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
#' @return List with the selected index/label, `naive` and `debiased` estimates
#'   (effect scale, with approximate 95% CIs), `selection_bias`, `fixed_bias`,
#'   `selection_rate`, the `gate` settings, `harm_flag`, family/subgroup sizes,
#'   and `timing_seconds`. The de-biased CI uses the subgroup's robust SE and is
#'   narrower than the Tier-1 infinitesimal-jackknife CI.  When
#'   `include_complement = TRUE`, a `complement` element carries the complement
#'   subgroup's `naive`/`debiased` estimates and bias terms in the same form.
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
                           seed = NULL) {
  gate <- match.arg(gate); reselection <- match.arg(reselection)
  multiplier <- match.arg(multiplier)
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
  beta_naive <- bh[sel]
  beta_deb   <- beta_naive - selection_bias - (if (is.finite(fixed_bias)) fixed_bias else 0)
  se <- sdv[sel]

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
      bdc <- bnc - (if (is.finite(selbias_c)) selbias_c else 0) -
                   (if (is.finite(fixed_c)) fixed_c else 0)
      sec <- sdv_c[sel]
      complement <- list(
        naive    = list(est = to_eff(bnc),
                        lower = to_eff(bnc - z975 * sec),
                        upper = to_eff(bnc + z975 * sec)),
        debiased = list(est = to_eff(bdc),
                        lower = to_eff(bdc - z975 * sec),
                        upper = to_eff(bdc + z975 * sec),
                        lower_1s = to_eff(bdc - stats::qnorm(0.95) * sec)),
        selection_bias = selbias_c, fixed_bias = fixed_c,
        n = Nall - sz[sel])
    } else {
      complement <- list(note = "complement subgroup could not be fit")
    }
  }

  list(
    selected_index = sel, selected_label = asm$names[sel],
    measure = spec$effect_measure, log_scale = log_scale,
    naive    = list(est = to_eff(beta_naive),
                    lower = to_eff(beta_naive - z975 * se),
                    upper = to_eff(beta_naive + z975 * se)),
    debiased = list(est = to_eff(beta_deb),
                    lower = to_eff(beta_deb - z975 * se),
                    upper = to_eff(beta_deb + z975 * se),
                    lower_1s = to_eff(ci_lo_1s)),
    selection_bias = selection_bias, fixed_bias = fixed_bias,
    selection_rate = selection_rate,
    complement = complement,
    gate = gate_meta, harm_flag = isTRUE(flag),
    n_family = length(asm$names), n_selected = sz[sel],
    timing_seconds = timing)
}
