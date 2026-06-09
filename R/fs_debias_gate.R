# =============================================================================
# fs_debias_gate.R
#
# Tier-2 fast de-biased gate for a selected forestsearch subgroup.
#
# After forestsearch() selects a harm subgroup H-hat, this computes a fast
# multiplier-bootstrap approximation of the bootstrap bias-corrected treatment
# effect and flags whether that de-biased estimate is still consistent with
# harm.  It is a leading-order approximation of the full bootstrap correction
# (forestsearch_bootstrap_dofuture(), "Tier 1") and does NOT replace it.
#
# The expensive quantity -- the per-subject treatment dfbeta for each candidate
# -- is exactly what the resample consistency engine already computes via
# .consistency_glm_pieces() / .consistency_cox_pieces().  Phase 1 re-derives
# those pieces from the candidate family (a few seconds at the ACTG175 scale);
# Phase 2 retains them in the consistency engine so the gate is near-zero cost.
#
# Closed form (see consistency_resample.R / the resampling-theory note):
#   beta*(b) = beta_hat + D(b),   D_g(b) = sum_i G_i(b) dfbeta_{g,i}
#   selection bias = mean_b D_{H*_b}(b)        (winner's curse of the re-selected
#                                               candidate on draw b)
#   de-biased beta = beta_hat(H-hat) - selection_bias - mean_b D_{H-hat}(b)
#
# Variance: infinitesimal jackknife (Leon et al. 2024, Eq. VInfJ_bc), computed
# from the same draws -- the centered multipliers play the role of the bootstrap
# multiplicities (K*_bi - Kbar*_i) and the per-draw residual is
#   r_b = (selection_bias + fixed_bias) - D_{H*_b}(b) - D_{H-hat}(b).
# This is the leading-order analogue of the Tier-1 IJ variance; with
# ci_method = "wald" the gate falls back to the subgroup robust SE.
# =============================================================================


#' Default near-null gate threshold for a measure
#'
#' @param effect_measure Character measure label (`"OR"`, `"RR"`, `"IRR"`,
#'   `"HR"`, `"RD"`, `"MD"`, ...) or `NULL`.
#' @return `1` for ratio measures (effect scale), `0` for differences.
#' @keywords internal
.fs_dg_gate_null <- function(effect_measure) {
  if (is.null(effect_measure)) return(1)
  if (effect_measure %in% c("OR", "RR", "IRR", "HR")) 1 else 0
}


#' Per-candidate influence pieces (Cox / GLM dispatch)
#'
#' Thin dispatch onto the existing forestsearch internals so the gate uses the
#' identical treatment `dfbeta`, `beta_hat`, and robust `sigma_D` the resample
#' consistency engine uses.
#'
#' @param df_sub Data frame for one candidate subgroup's members.
#' @param spec List with `outcome_type`, `effect_measure`, `treat.name`,
#'   `outcome.name`, `event.name`, `offset.name`, `adjust_covariates`,
#'   `adverse_outcome`.
#' @return Output of [.consistency_glm_pieces()] / [.consistency_cox_pieces()],
#'   or `NULL`.
#' @keywords internal
.fs_dg_pieces <- function(df_sub, spec) {
  if (identical(spec$outcome_type, "survival")) {
    .consistency_cox_pieces(
      df_sub, tte.name = spec$outcome.name, event.name = spec$event.name,
      treat.name = spec$treat.name, adjust_covariates = spec$adjust_covariates)
  } else {
    .consistency_glm_pieces(
      df_sub, outcome_type = spec$outcome_type, effect_measure = spec$effect_measure,
      treat.name = spec$treat.name, outcome.name = spec$outcome.name,
      offset.name = spec$offset.name, adjust_covariates = spec$adjust_covariates,
      adverse_outcome = if (is.null(spec$adverse_outcome)) TRUE else spec$adverse_outcome)
  }
}


#' Assemble the influence matrix from candidate membership lists
#'
#' For each candidate (a vector of row indices into `df`), fit the treatment
#' effect once and place its `dfbeta` into the candidate's rows (zero elsewhere).
#' Candidates whose fit fails or whose `dfbeta` length does not match are dropped.
#'
#' @param df Analysis data frame.
#' @param candidates Named list of integer row-index vectors.
#' @param spec See [.fs_dg_pieces()].
#' @return List with `B` (`N x S_kept`), `beta_hat`, `sigma_D`, `sizes`,
#'   `names`, `keep`, `log_scale`.
#' @keywords internal
.fs_dg_assemble <- function(df, candidates, spec) {
  N <- nrow(df); S <- length(candidates)
  B  <- matrix(0, N, S)
  bh <- rep(NA_real_, S); sdv <- rep(NA_real_, S); sz <- rep(NA_integer_, S)
  keep <- logical(S); log_scale <- TRUE
  for (g in seq_len(S)) {
    idx <- candidates[[g]]
    if (length(idx) < 6L) next
    pc <- tryCatch(.fs_dg_pieces(df[idx, , drop = FALSE], spec), error = function(e) NULL)
    if (is.null(pc) || length(pc$dfbeta) != length(idx)) next
    B[idx, g] <- pc$dfbeta
    bh[g] <- pc$beta_hat; sdv[g] <- pc$sigma_D; sz[g] <- length(idx)
    keep[g] <- TRUE; log_scale <- isTRUE(pc$log_scale)
  }
  list(B = B[, keep, drop = FALSE], beta_hat = bh[keep], sigma_D = sdv[keep],
       sizes = sz[keep], names = names(candidates)[keep], keep = keep,
       log_scale = log_scale)
}


#' Mean-zero unit-variance multiplier vectors (one column per draw)
#' @keywords internal
.fs_dg_multipliers <- function(n, draws, type) {
  switch(type,
    rademacher = matrix(sample(c(-1, 1), n * draws, replace = TRUE), n, draws),
    gaussian   = matrix(stats::rnorm(n * draws), n, draws),
    poisson    = matrix(stats::rpois(n * draws, 1) - 1, n, draws),
    stop("unknown multiplier: ", type))
}


#' Apply a selection rule among passing candidates on one draw
#' @keywords internal
.fs_dg_select <- function(beta, zcons, sizes, passers, rule, nbhd) {
  if (!length(passers)) return(NA_integer_)
  pick <- switch(rule,
    maxcons  = passers[which.max(zcons[passers])],
    maxeff   = passers[which.max(beta[passers])],
    maxSG    = passers[which.max(sizes[passers])],
    minSG    = passers[which.min(sizes[passers])],
    effMaxSG = { b <- passers[beta[passers] >= max(beta[passers]) - nbhd]
                 b[which.max(sizes[b])] },
    effMinSG = { b <- passers[beta[passers] >= max(beta[passers]) - nbhd]
                 b[which.min(sizes[b])] },
    stop("unknown reselection rule: ", rule))
  as.integer(pick)
}


#' Infinitesimal-jackknife variance of the de-biased (bagged) estimate
#'
#' Implements Leon et al. (2024) Eq. (VInfJ)/(VInfJ_bc) in multiplier form.
#' The multiplier matrix `Xi` supplies the centered bootstrap multiplicities
#' \eqn{K^*_{bi} - \bar K^*_i}, and `r` is the per-draw residual
#' \eqn{r_b = \hat\beta(\widehat H) - \eta^*_b(\widehat H^*_b) - \eta^*_b(\widehat H) - \hat\beta^*(\widehat H)},
#' i.e. `(selection_bias + fixed_bias) - D_{H*_b}(b) - D_{H}(b)`, evaluated only
#' on the draws in `ok` that produced a re-selected winner.
#'
#' @param Xi `N x draws` multiplier matrix (one column per draw).
#' @param r Length-`draws` residual vector (`NA` where no winner).
#' @param ok Integer indices of usable draws.
#' @return List with `tilde_V` (raw IJ), `hat_V` (Wager 2014 bias-corrected),
#'   and `B_ok` (number of usable draws).
#' @keywords internal
.fs_dg_ij_var <- function(Xi, r, ok) {
  if (length(ok) < 2L)
    return(list(tilde_V = NA_real_, hat_V = NA_real_, B_ok = length(ok)))
  Xk  <- Xi[, ok, drop = FALSE]            # N x B_ok
  rb  <- r[ok]
  Bok <- length(ok)
  Xc  <- Xk - rowMeans(Xk)                 # centered multiplicities over used draws
  cov_i   <- as.numeric(Xc %*% rb) / Bok   # cov_i = (1/B) sum_b (K* - Kbar) r_b
  tilde_V <- sum(cov_i^2)
  hat_V   <- tilde_V - (nrow(Xi) / Bok) * mean(rb^2)
  list(tilde_V = tilde_V, hat_V = hat_V, B_ok = Bok)
}


#' Resolve a de-biased SE from the IJ variance, with graceful fallback
#'
#' Prefers the bias-corrected IJ variance; if it is non-positive (too few
#' draws), falls back to the raw IJ variance, then to the subgroup robust SE.
#'
#' @param ij Output of [.fs_dg_ij_var()].
#' @param se_fallback Robust subgroup SE (`sigma_D`) used as last resort.
#' @return List with `se`, `var`, and `source`
#'   (`"ij"`, `"ij_raw"`, or `"wald_fallback"`).
#' @keywords internal
.fs_dg_se_from_ij <- function(ij, se_fallback) {
  if (is.finite(ij$hat_V) && ij$hat_V > 0)
    return(list(se = sqrt(ij$hat_V), var = ij$hat_V, source = "ij"))
  if (is.finite(ij$tilde_V) && ij$tilde_V > 0)
    return(list(se = sqrt(ij$tilde_V), var = ij$tilde_V, source = "ij_raw"))
  list(se = se_fallback, var = se_fallback^2, source = "wald_fallback")
}


#' Fast de-biased gate for a selected forestsearch subgroup
#'
#' Computes a multiplier-bootstrap approximation of the bootstrap
#' bias-corrected treatment effect for the selected subgroup and flags whether
#' it is still consistent with harm.  Reuses the resample consistency engine's
#' per-candidate treatment `dfbeta`; see the file header for the closed form.
#'
#' @param df Analysis data frame (the standardized `df.fs` inside
#'   `forestsearch()`); must contain the columns named in `spec`.
#' @param candidates Named list of integer row-index vectors, one per screened
#'   candidate subgroup (the family the selection rule chose among).
#' @param spec List with `outcome_type`, `effect_measure`, `treat.name`,
#'   `outcome.name`, `event.name`, `offset.name`, `adjust_covariates`,
#'   `adverse_outcome`.
#' @param selected_members Integer row indices of the observed selected subgroup
#'   (`which(grp.consistency$sg.harm.id == 1)`).
#' @param c_screen,c_consistency Screening and consistency thresholds on the
#'   **comparison scale** (log for ratio measures, identity for differences) --
#'   i.e. `forestsearch()`'s resolved `effect_threshold` / `consistency_threshold`.
#' @param p_star Consistency-rate cutoff (`pconsistency.threshold`).
#' @param t_gate Gate threshold on the **effect scale**; `NULL` uses the
#'   near-null default (`1` ratio, `0` difference). Set near the null, not at the
#'   screen, since the correction over-shrinks true effects.
#' @param gate `"point"` (de-biased point estimate clears `t_gate`) or `"ci"`
#'   (one-sided selection-adjusted lower bound clears `t_gate`).
#' @param reselection Bootstrap re-selection rule for the bias term; default
#'   `"maxcons"` (validated). `"effMaxSG"` etc. are available but approximate.
#' @param effect_neighborhood Band for the `eff*SG` re-selection rules.
#' @param draws,multiplier,seed Multiplier-bootstrap controls. `multiplier`
#'   defaults to `"poisson"` (mimics the nonparametric bootstrap).
#' @param include_complement Logical.  When `TRUE`, also de-bias the complement
#'   subgroup (everyone not in the selected subgroup).  The complement is not
#'   selected independently -- its bias is induced by selecting the subgroup --
#'   so on each multiplier draw the complement of the re-selected winner is
#'   tracked and de-biased with the complement's own influence.  Complements are
#'   fit only for candidates that win across draws (plus the selected one), so
#'   the extra cost is small.  Default `FALSE`.
#' @param ci_method `"ij"` (default) bases the **de-biased** CI on the
#'   infinitesimal-jackknife variance (Leon et al. 2024, Eq. VInfJ_bc), computed
#'   from the same multiplier draws -- the leading-order analogue of the Tier-1
#'   interval.  `"wald"` uses the subgroup robust SE (`sigma_D`).  The naive CI
#'   always uses the robust SE.
#' @return List with the selected index/label, `naive` and `debiased` estimates
#'   (effect scale, with approximate 95% CIs), `selection_bias`, `fixed_bias`,
#'   `selection_rate`, the `gate` settings, `harm_flag`, family/subgroup sizes,
#'   and `timing_seconds`. The `debiased` element carries `se_ij`, `se_wald`,
#'   `var_ij`, and `ij_source`; its CI uses the IJ SE under the default
#'   `ci_method = "ij"` (the Tier-1 analogue) and the robust SE under `"wald"`.
#'   When `include_complement = TRUE`, a `complement` element carries the
#'   complement subgroup's `naive`/`debiased` estimates and bias terms in the
#'   same form, including its own IJ variance.
#' @keywords internal
fs_debias_gate <- function(df, candidates, spec, selected_members,
                           c_screen, c_consistency = 0, p_star = 0.90,
                           t_gate = NULL, gate = c("point", "ci"),
                           reselection = c("maxcons", "maxeff", "maxSG",
                                           "minSG", "effMaxSG", "effMinSG"),
                           effect_neighborhood = 0.10,
                           draws = 2000L,
                           multiplier = c("poisson", "gaussian", "rademacher"),
                           include_complement = FALSE,
                           ci_method = c("ij", "wald"),
                           seed = NULL) {
  gate <- match.arg(gate); reselection <- match.arg(reselection)
  multiplier <- match.arg(multiplier); ci_method <- match.arg(ci_method)
  if (!is.null(seed)) set.seed(seed)
  df <- as.data.frame(df)
  if (!length(candidates)) candidates <- list()
  if (is.null(names(candidates)))
    names(candidates) <- paste0("cand", seq_along(candidates))

  # Ensure the observed selected subgroup is in the family and is the target.
  H_lab <- ".selected_H"
  hit <- if (length(candidates))
    which(vapply(candidates, function(ix) setequal(ix, selected_members), logical(1)))
  else integer(0)
  if (length(hit)) sel_lab <- names(candidates)[hit[1]]
  else { candidates[[H_lab]] <- selected_members; sel_lab <- H_lab }

  asm <- .fs_dg_assemble(df, candidates, spec)
  if (is.null(t_gate)) t_gate <- if (asm$log_scale) 1 else 0
  sel <- match(sel_lab, asm$names)

  gate_meta <- list(t_gate = t_gate, type = gate, reselection = reselection,
                    multiplier = multiplier, draws = as.integer(draws))
  if (is.na(sel)) {
    return(list(selected_index = NA_integer_, selected_label = sel_lab,
                harm_flag = NA, gate = gate_meta,
                note = "selected subgroup could not be fit in the reconstructed family",
                n_family = length(asm$names)))
  }

  B <- asm$B; bh <- asm$beta_hat; sdv <- asm$sigma_D; sz <- asm$sizes
  log_scale <- asm$log_scale
  to_eff    <- function(x) if (log_scale) exp(x) else x
  z      <- stats::qnorm((1 + p_star) / 2)
  z975   <- stats::qnorm(0.975)
  t_g    <- pmax(c_screen, c_consistency + z * sdv)

  t0 <- proc.time()
  Xi <- .fs_dg_multipliers(nrow(B), draws, multiplier)
  P  <- crossprod(B, Xi)                 # S x draws : D_g(b)
  beta_star <- bh + P
  sel_bias <- rep(NA_real_, draws)
  winner   <- rep(NA_integer_, draws)    # which candidate won on draw b
  for (b in seq_len(draws)) {
    bs <- beta_star[, b]
    pass <- which(bs >= t_g)
    if (!length(pass)) next
    s <- .fs_dg_select(bs, (bs - c_consistency) / sdv, sz, pass, reselection,
                       effect_neighborhood)
    if (!is.na(s)) { sel_bias[b] <- P[s, b]; winner[b] <- s }
  }
  timing <- as.numeric((proc.time() - t0)["elapsed"])

  selection_rate <- mean(!is.na(sel_bias))
  selection_bias <- mean(sel_bias, na.rm = TRUE)
  fixed_bias     <- mean(P[sel, ])
  fb             <- if (is.finite(fixed_bias)) fixed_bias else 0
  beta_naive <- bh[sel]
  beta_deb   <- beta_naive - selection_bias - fb
  se_wald    <- sdv[sel]

  # Infinitesimal-jackknife variance of the de-biased estimate (Eq. VInfJ_bc),
  # from the same draws: r_b = (selection_bias + fixed_bias) - D_{H*_b}(b) - D_H(b).
  r_H   <- (selection_bias + fb) - sel_bias - P[sel, ]
  ijH   <- .fs_dg_ij_var(Xi, r_H, which(is.finite(sel_bias)))
  se_ij <- .fs_dg_se_from_ij(ijH, se_wald)
  se    <- if (ci_method == "ij") se_ij$se else se_wald

  gate_cmp <- if (log_scale) log(t_gate) else t_gate
  ci_lo_1s <- beta_deb - stats::qnorm(0.95) * se
  flag <- if (gate == "point") (beta_deb >= gate_cmp) else (ci_lo_1s >= gate_cmp)

  # ---------------------------------------------------------------------------
  # Complement subgroup (optional).  The complement is induced by the selection,
  # not chosen independently, so its bias is the perturbation of the complement
  # of the re-selected winner on each draw.  Fit complements only for candidates
  # that win across draws (plus the selected one) to keep the cost small.
  # ---------------------------------------------------------------------------
  complement <- NULL
  if (isTRUE(include_complement)) {
    kept   <- candidates[asm$keep]              # aligns with asm columns
    Ncol   <- length(asm$names)
    Nall   <- nrow(df)
    winset <- sort(unique(c(winner[!is.na(winner)], sel)))
    Bc   <- matrix(0, Nall, Ncol)
    bh_c <- rep(NA_real_, Ncol); sdv_c <- rep(NA_real_, Ncol)
    for (w in winset) {
      comp_idx <- setdiff(seq_len(Nall), kept[[w]])
      if (length(comp_idx) < 6L) next
      pcc <- tryCatch(.fs_dg_pieces(df[comp_idx, , drop = FALSE], spec),
                      error = function(e) NULL)
      if (is.null(pcc) || length(pcc$dfbeta) != length(comp_idx)) next
      Bc[comp_idx, w] <- pcc$dfbeta
      bh_c[w] <- pcc$beta_hat; sdv_c[w] <- pcc$sigma_D
    }
    Pc <- crossprod(Bc, Xi)                      # Ncol x draws : D_{complement}(b)
    ok <- which(!is.na(winner))
    selb_c <- rep(NA_real_, draws)
    if (length(ok)) {
      vals <- Pc[cbind(winner[ok], ok)]
      vals[is.na(bh_c[winner[ok]])] <- NA_real_  # winner's complement not fit
      selb_c[ok] <- vals
    }
    selbias_c <- mean(selb_c, na.rm = TRUE)
    fixed_c   <- mean(Pc[sel, ])
    bnc       <- bh_c[sel]
    if (is.finite(bnc)) {
      sbc <- if (is.finite(selbias_c)) selbias_c else 0
      fcc <- if (is.finite(fixed_c)) fixed_c else 0
      bdc <- bnc - sbc - fcc
      sec <- sdv_c[sel]
      r_c   <- (sbc + fcc) - selb_c - Pc[sel, ]
      ijC   <- .fs_dg_ij_var(Xi, r_c, which(is.finite(selb_c)))
      se_ijc <- .fs_dg_se_from_ij(ijC, sec)
      sec_used <- if (ci_method == "ij") se_ijc$se else sec
      complement <- list(
        naive    = list(est = to_eff(bnc),
                        lower = to_eff(bnc - z975 * sec),
                        upper = to_eff(bnc + z975 * sec)),
        debiased = list(est = to_eff(bdc),
                        lower = to_eff(bdc - z975 * sec_used),
                        upper = to_eff(bdc + z975 * sec_used),
                        lower_1s = to_eff(bdc - stats::qnorm(0.95) * sec_used),
                        se = sec_used, se_ij = se_ijc$se, se_wald = sec,
                        var_ij = se_ijc$var, ij_source = se_ijc$source,
                        ij_draws = ijC$B_ok),
        selection_bias = selbias_c, fixed_bias = fixed_c,
        n = Nall - sz[sel])
    } else {
      complement <- list(note = "complement subgroup could not be fit")
    }
  }

  list(
    selected_index = sel, selected_label = asm$names[sel],
    measure = spec$effect_measure, log_scale = log_scale,
    ci_method = ci_method,
    naive    = list(est = to_eff(beta_naive),
                    lower = to_eff(beta_naive - z975 * se_wald),
                    upper = to_eff(beta_naive + z975 * se_wald)),
    debiased = list(est = to_eff(beta_deb),
                    lower = to_eff(beta_deb - z975 * se),
                    upper = to_eff(beta_deb + z975 * se),
                    lower_1s = to_eff(ci_lo_1s),
                    se = se, se_ij = se_ij$se, se_wald = se_wald,
                    var_ij = se_ij$var, ij_source = se_ij$source,
                    ij_draws = ijH$B_ok),
    selection_bias = selection_bias, fixed_bias = fixed_bias,
    selection_rate = selection_rate,
    complement = complement,
    gate = gate_meta, harm_flag = isTRUE(flag),
    n_family = length(asm$names), n_selected = sz[sel],
    timing_seconds = timing)
}
