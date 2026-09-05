# =============================================================================
# fs_mr_inference.R
#
# Multiplier resampling (MR) for a selected forestsearch subgroup.
#
# After forestsearch() selects a harm subgroup H-hat, this computes a fast
# multiplier-bootstrap approximation of the bootstrap bias-corrected treatment
# effect and flags whether that de-biased estimate is still consistent with
# harm.  It is a leading-order approximation of the full bootstrap (FB)
# correction (forestsearch_bootstrap_dofuture()) and does NOT replace it.
#
# The expensive quantity -- the per-subject treatment dfbeta for each candidate
# -- is exactly what the resample consistency engine already computes via
# .consistency_glm_pieces() / .consistency_cox_pieces().  Phase 1 re-derives
# those pieces from the candidate family (a few seconds at the ACTG175 scale);
# Phase 2 retains them in the consistency engine so MR is near-zero cost.
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
# This is the leading-order analogue of the FB IJ variance; with
# ci_method = "wald" MR falls back to the subgroup robust SE.
# =============================================================================


#' Default near-null harm-confirmation threshold for a measure
#'
#' The value `t_confirm` takes when left `NULL`: the null effect on the effect
#' scale, i.e. no effect at all rather than the screening threshold.
#'
#' @param effect_measure Character measure label (`"OR"`, `"RR"`, `"IRR"`,
#'   `"HR"`, `"RD"`, `"MD"`, ...) or `NULL`.
#' @return `1` for ratio measures (effect scale), `0` for differences.
#' @keywords internal
.fs_mr_confirm_null <- function(effect_measure) {
  if (is.null(effect_measure)) return(1)
  if (effect_measure %in% c("OR", "RR", "IRR", "HR")) 1 else 0
}


#' Per-candidate influence pieces (Cox / GLM dispatch)
#'
#' Thin dispatch onto the existing forestsearch internals so MR uses the
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
.fs_mr_pieces <- function(df_sub, spec) {
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
#' @param spec See [.fs_mr_pieces()].
#' @return List with `B` (`N x S_kept`), `beta_hat`, `sigma_D`, `sizes`,
#'   `names`, `keep`, `log_scale`.
#' @keywords internal
.fs_mr_assemble <- function(df, candidates, spec) {
  N <- nrow(df); S <- length(candidates)
  B  <- matrix(0, N, S)
  bh <- rep(NA_real_, S); sdv <- rep(NA_real_, S); sz <- rep(NA_integer_, S)
  keep <- logical(S); log_scale <- TRUE
  for (g in seq_len(S)) {
    idx <- candidates[[g]]
    if (length(idx) < 6L) next
    pc <- tryCatch(.fs_mr_pieces(df[idx, , drop = FALSE], spec), error = function(e) NULL)
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
.fs_mr_multipliers <- function(n, draws, type) {
  switch(type,
    rademacher = matrix(sample(c(-1, 1), n * draws, replace = TRUE), n, draws),
    gaussian   = matrix(stats::rnorm(n * draws), n, draws),
    poisson    = matrix(stats::rpois(n * draws, 1) - 1, n, draws),
    stop("unknown multiplier: ", type))
}


#' Apply a selection rule among passing candidates on one draw
#'
#' For the size-within-band rules (`effMaxSG`/`effMinSG`) the inclusion band is
#' built to match the search's `selection_rule`:
#'   * `"neighborhood"` - effect within `nbhd` of the max, on the NATURAL effect
#'     scale and multiplicative (`eff >= (1 - nbhd) * max(eff)`), matching
#'     `sort_subgroups()` exactly (G3);
#'   * `"pareto"` - the 2-D non-dominated set in (effect, size), via the same
#'     dominance core the search uses (`.pareto_dominated_xy()`) (G2);
#'   * `"both"` - the intersection of the two.
#' `beta` is on the working scale (log for ratio measures); `log_scale` controls
#' the conversion to the natural effect for the band.
#' @keywords internal
.fs_mr_select <- function(beta, zcons, sizes, passers, rule, nbhd,
                          selection_rule = "neighborhood", log_scale = TRUE) {
  if (!length(passers)) return(NA_integer_)
  .inband <- function() {
    eff <- if (log_scale) exp(beta[passers]) else beta[passers]  # natural effect
    sz  <- sizes[passers]
    # Band from the shared helper, the same one the consistency engine, DINA
    # and (since the frontier fix) GRF use.  This replaces a local copy that
    # reimplemented the same three rules.  The local copy differed in its NA
    # handling: it used a plain max() with no is.na() guard, so a single
    # non-finite effect made the threshold NA, made the whole band NA, and --
    # via the emptiness fallback below -- silently degraded the band to ALL
    # passers rather than to a sensible band.  The shared helper uses
    # max(na.rm = TRUE) and guards !is.na(), degrading one candidate instead
    # of the draw.  That is a bug fix, not a preference.
    ib <- .compute_inclusion_band(hr_vec = eff, n_vec = sz,
                                  selection_rule = selection_rule,
                                  effect_neighborhood = nbhd) == 1L
    # EMPTY-BAND FALLBACK -- deliberately here, at the call site, not inside
    # the shared helper.  The two callers want different things and the
    # difference is not reconcilable:
    #
    #   * In the consistency/DINA sort an all-zero band is harmless: -in_band
    #     is the leading sort key and the next key breaks the tie, so nothing
    #     is lost.
    #   * Here the band is a FILTER inside a draw.  Returning nothing loses
    #     that draw from sel_bias -- and it loses exactly the draws where the
    #     perturbed effects were extreme enough to empty the band, which
    #     biases the correction rather than protecting it.  Same failure mode
    #     as an over-tight admission floor.
    #
    # Keeping it visible at the call site records that it is a deliberate
    # choice for MR, not a property of the band.
    if (!any(ib)) ib <- rep(TRUE, length(passers))
    passers[ib]
  }
  pick <- switch(rule,
    maxcons  = passers[which.max(zcons[passers])],
    maxeff   = passers[which.max(beta[passers])],
    maxSG    = passers[which.max(sizes[passers])],
    minSG    = passers[which.min(sizes[passers])],
    effMaxSG = { b <- .inband(); b[which.max(sizes[b])] },
    effMinSG = { b <- .inband(); b[which.min(sizes[b])] },
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
.fs_mr_ij_var <- function(Xi, r, ok) {
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
#' @param ij Output of [.fs_mr_ij_var()].
#' @param se_fallback Robust subgroup SE (`sigma_D`) used as last resort.
#' @return List with `se`, `var`, and `source`
#'   (`"ij"`, `"ij_raw"`, or `"wald_fallback"`).
#' @keywords internal
.fs_mr_se_from_ij <- function(ij, se_fallback) {
  if (is.finite(ij$hat_V) && ij$hat_V > 0)
    return(list(se = sqrt(ij$hat_V), var = ij$hat_V, source = "ij"))
  if (is.finite(ij$tilde_V) && ij$tilde_V > 0)
    return(list(se = sqrt(ij$tilde_V), var = ij$tilde_V, source = "ij_raw"))
  list(se = se_fallback, var = se_fallback^2, source = "wald_fallback")
}


#' Multiplier resampling (MR) for a selected forestsearch subgroup
#'
#' Computes a multiplier-bootstrap approximation of the bootstrap
#' bias-corrected treatment effect for the selected subgroup and flags whether
#' it is still consistent with harm.  Reuses the resample consistency engine's
#' per-candidate treatment `dfbeta`; see the file header for the closed form.
#'
#' This is **post-selection inference on a completed analysis**: it runs after
#' [forestsearch()] has chosen the subgroup and cannot change that choice.  It
#' is the fast, refit-free counterpart of the full bootstrap (FB,
#' [forestsearch_bootstrap_dofuture()]), which re-runs the entire search in
#' every replicate; MR perturbs the influence contributions of a single fit
#' instead.  MR approximates FB to leading order and does not replace it.
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
#' @param admission The resolved admission set, as returned by
#'   `.fs_resolve_admission()`: a list with `effect_floor` (numeric on the
#'   **comparison scale** -- log for ratio measures -- or `NULL`) and
#'   `consistency` (a list with `c_cons` and `p_star`, or `NULL`).
#'
#'   This replaces the former `c_screen` / `c_consistency` / `p_star`
#'   arguments, from which MR used to rebuild
#'   `t_g <- pmax(c_screen, c_consistency + z * sigma_D)` at this site. That
#'   reconstruction was the defect: the identifier resolved its floors in one
#'   place and MR re-derived them in another, so the two could disagree about
#'   which candidates were ever in contention -- and under
#'   `sg_focus = "maxeff"` they did. MR is a linearization of the selection
#'   *map*, which is a ranking and a domain; passing the resolved domain makes
#'   the two equal by construction.
#'
#'   `NULL` means the floor does not apply, and is handled by a separate code
#'   path rather than by sentinel arithmetic. With neither floor present there
#'   is no admission filter at all.
#' @param t_confirm Harm-confirmation threshold on the **effect scale** (HR/OR/
#'   RR/IRR, or RD/MD for differences -- not the working log scale). `NULL` uses
#'   the near-null default: `1` for ratio measures, `0` for differences. Set it
#'   near the null rather than at the screening threshold, because the
#'   de-biasing correction over-shrinks true effects.
#' @param confirm_rule Which harm-confirmation rule to apply to the de-biased
#'   estimate: `"point"` requires the de-biased point estimate to clear
#'   `t_confirm`; `"ci"` requires the one-sided 95% selection-adjusted lower
#'   bound to clear it, and is therefore the stricter of the two.
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
#' @param return_reselection Logical. When `TRUE`, the return list gains a
#'   `reselection` element exposing the per-draw re-selection already computed
#'   for the bias term: `winner` (integer vector of length `draws`, the index
#'   of the re-selected candidate in the kept family on each draw, `NA` where
#'   no draw winner exists) and `p_hat` (named numeric over the kept family,
#'   `tabulate(winner) / draws` -- so `sum(p_hat) == selection_rate`, and the
#'   frequencies sum to 1 exactly when every draw produced a winner). The
#'   default `FALSE` reproduces the previous return object exactly; nothing in
#'   the arithmetic depends on this switch.
#' @param ci_method `"ij"` (default) bases the **de-biased** CI on the
#'   infinitesimal-jackknife variance (Leon et al. 2024, Eq. VInfJ_bc), computed
#'   from the same multiplier draws -- the leading-order analogue of the FB
#'   interval.  `"wald"` uses the subgroup robust SE (`sigma_D`).  The naive CI
#'   always uses the robust SE.  `"field"` computes everything the `"ij"` path
#'   computes -- the `debiased` element is identical -- and additionally runs
#'   the field-calibrated interval (method proposal,
#'   `dev/tasks/TASK_mr_field_vs_guohe_2026-09-05.md`), returned as a `field`
#'   element: Gaussian-multiplier perturbations of the shrunk candidate-effect
#'   field (`w = beta_hat` with the winner's entry replaced by the de-biased
#'   estimate) are pushed through the configured re-selection map, with an
#'   inner Monte Carlo estimating the selection drift, and the resulting
#'   `Lambda*` distribution is inverted around the de-biased estimate
#'   (basic-bootstrap form).  Drawn under `seed + 900000L` when `seed` is
#'   given, after -- never inside -- the main multiplier stream, so `"ij"` and
#'   `"wald"` output is unaffected byte-for-byte.
#' @param field_R_out,field_R_in Outer and inner Monte Carlo sizes for
#'   `ci_method = "field"` (defaults 1000 / 500); ignored otherwise.  The
#'   inner draws are shared across outer draws.
#' @param field_uniform Logical (default `FALSE`); only consulted under
#'   `ci_method = "field"`.  When `TRUE`, after the field block completes the
#'   gate additionally runs the uniform (kappa) calibration
#'   (`fs_mr_field_uniform()`, method proposal
#'   `dev/tasks/TASK_mr_field_uniform_2026-09-05.md`) and attaches its result
#'   as `field$uniform`: the smallest widening factor `kappa` in `[1, 2]`
#'   such that the widened two-sided interval attains 95% coverage uniformly
#'   over the **winner-profile protection family** -- hypothetical true
#'   fields equal to the shrunk field `w` with the winner's entry set to the
#'   runner-up level plus `delta` winner-SEs, `delta` on a 0-4 grid --
#'   computed from the trial's own influence structure (`Sigma_hat` via
#'   `B_eff`), restricted to the **mass-carrying candidate set** (smallest
#'   `M <= 12` candidates holding >= 99% of the re-selection mass, winner
#'   always included).  Guarantee as documented there: the one-sided field
#'   bound is uniformly valid; the two-sided interval widened by the
#'   reported `kappa` is uniformly valid over the winner-profile family; the
#'   plain quantile interval is approximate.  The sweep draws under
#'   `seed + 910000L` after -- never inside -- the field's own stream, so
#'   every other output, the field block included, is byte-identical whether
#'   or not it runs.
#' @return List with the selected index/label, `naive` and `debiased` estimates
#'   (effect scale, with approximate 95% CIs), `selection_bias`, `fixed_bias`,
#'   `selection_rate`, `mean_r`, `mean_r_c`, the `settings` actually used (`t_confirm`,
#'   `confirm_rule`, `reselection`, `selection_rule`, `multiplier`, `draws`),
#'   `harm_flag`, family/subgroup sizes, and `timing_seconds`. The `debiased`
#'   element carries `se_ij`, `se_wald`, `var_ij`, and `ij_source`; its CI uses
#'   the IJ SE under the default `ci_method = "ij"` (the FB analogue) and the
#'   robust SE under `"wald"`.
#'   When `include_complement = TRUE`, a `complement` element carries the
#'   complement subgroup's `naive`/`debiased` estimates and bias terms in the
#'   same form, including its own IJ variance.  The complement's de-biased CI
#'   follows the same SE convention as the winner's: the IJ SE under
#'   `ci_method = "ij"` and `"field"` (under `"field"` the complement, like the
#'   `debiased` element, is identical to the `"ij"` output), the robust SE
#'   under `"wald"`.
#'
#'   `mean_r` is the mean of the IJ residual \eqn{r_b} over `ok_H` -- exactly
#'   the draws entering the variance, not all `draws`, since on an excluded
#'   draw \eqn{r_b} is not a meaningful quantity. **The invariant is that it is
#'   zero by construction whenever both bias terms share a denominator**, which
#'   is the convention the package implements: `selection_bias` and
#'   `fixed_bias` both average over the draws that produced a winner, and the
#'   IJ runs on that same set. A non-zero value therefore means the two terms
#'   are being normalised differently somewhere -- the defect corrected in
#'   `dad0415`, whose signature was precisely a non-zero residual mean.
#'   `mean_r_c` is the same quantity for the complement over `use_c`, and is
#'   `NA_real_` when no complement was fit. Both are diagnostics of the
#'   correction's internal consistency rather than properties of the estimate,
#'   which is why they sit beside `selection_rate` and not inside `debiased`.
#'
#'   Note the mixed scales in the flat `debiased` list: `est`/`lower`/`upper`/
#'   `lower_1s` are on the **effect** scale, while `se`/`se_ij`/`se_wald`/
#'   `var_ij` and the two bias terms are on the **working** (log, for ratio
#'   measures) scale.
#'
#'   When the selected subgroup cannot be fit in the reconstructed family the
#'   return is a short variant carrying only `selected_index` (`NA`),
#'   `selected_label`, `harm_flag` (`NA`), `settings`, `note` and `n_family`;
#'   consumers must tolerate that shape.
#'
#'   Under `return_reselection = TRUE` the full return additionally carries a
#'   `reselection` element (see that argument); the short variant never does.
#'
#'   Under `ci_method = "field"` the full return additionally carries a
#'   `field` element: `lambda_mean`, `lambda_sd`/`se_field` and the
#'   `Lambda*` quantiles `q05/q25/q50/q75/q95/q025/q975` on the **working**
#'   scale; `n_out_used`, `n_in_used_mean`, `R_out`, `R_in`, `seed_offset`,
#'   `timing_seconds`; and, on the **effect** scale (the `debiased`
#'   convention), the second-order point estimate `est2` (de-biased estimate
#'   minus `lambda_mean`), the primary one-sided 95% bound `lower_1s`
#'   (de-biased minus `q95`), the two-sided quantile interval
#'   `lower_2s`/`upper_2s`, and the supplementary SE-type interval
#'   `lower_se`/`upper_se` around `est2`.
#' @section Alignment is assumed, not checked here:
#' This is an engine-level entry point. Its arguments are a candidate family,
#' a specification, and a selected membership vector -- the identifier
#' configuration that produced them is not visible, so the alignment
#' conditions MR requires (selection ranking on the inferential coefficient
#' \eqn{\hat\beta(g)}, and a fixed candidate family) cannot be verified at
#' this level and are assumed to have been established upstream.
#'
#' They are enforced by `.validate_mr_configuration()` at the three
#' configuration-visible entry points -- [forestsearch()] under
#' `mr_inference = TRUE`, and [forestsearch_bootstrap_dofuture()] and
#' [forestsearch_Kfold()] under `mr_in_replicates = TRUE`. Calling
#' `fs_mr_inference()` directly bypasses those guards: a family ranked on
#' DINA's native tau-hat, on GRF's doubly-robust score, or on a GRF
#' policy-tree objective will still produce numbers, but they do not de-bias
#' the reported effect.
#'
#' @seealso [forestsearch()] for the `mr_inference` switch and the
#'   vocabulary section; [mr_estimates_table()] to render the result;
#'   [forestsearch_bootstrap_dofuture()] for the full bootstrap (FB) this
#'   approximates; [fs_fdr_report()] which sweeps `c_confirm` thresholds.
#' @keywords internal
fs_mr_inference <- function(df, candidates, spec, selected_members,
                           admission,
                           t_confirm = NULL, confirm_rule = c("point", "ci"),
                           reselection = c("maxcons", "maxeff", "maxSG",
                                           "minSG", "effMaxSG", "effMinSG"),
                           effect_neighborhood = 0.10,
                           selection_rule = c("neighborhood", "pareto", "both"),
                           draws = 2000L,
                           multiplier = c("poisson", "gaussian", "rademacher"),
                           include_complement = FALSE,
                           ci_method = c("ij", "wald", "field"),
                           seed = NULL,
                           return_reselection = FALSE,
                           field_R_out = 1000L,
                           field_R_in = 500L,
                           field_uniform = FALSE) {
  confirm_rule <- match.arg(confirm_rule); reselection <- match.arg(reselection)
  selection_rule <- match.arg(selection_rule)
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

  asm <- .fs_mr_assemble(df, candidates, spec)
  if (is.null(t_confirm)) t_confirm <- if (asm$log_scale) 1 else 0
  sel <- match(sel_lab, asm$names)

  mr_settings <- list(t_confirm = t_confirm, confirm_rule = confirm_rule,
                      reselection = reselection,
                      selection_rule = selection_rule,
                      multiplier = multiplier, draws = as.integer(draws))
  if (is.na(sel)) {
    return(list(selected_index = NA_integer_, selected_label = sel_lab,
                harm_flag = NA, settings = mr_settings,
                note = "selected subgroup could not be fit in the reconstructed family",
                n_family = length(asm$names)))
  }

  B <- asm$B; bh <- asm$beta_hat; sdv <- asm$sigma_D; sz <- asm$sizes
  log_scale <- asm$log_scale
  to_eff    <- function(x) if (log_scale) exp(x) else x
  z975   <- stats::qnorm(0.975)

  # ---------------------------------------------------------------------------
  # ADMISSION SET
  # ---------------------------------------------------------------------------
  # MR linearizes the identifier's selection MAP, which is a ranking AND a
  # domain.  The domain is resolved once by .fs_resolve_admission() and carried
  # here; it is NOT reconstructed from raw thresholds at this site, because
  # reconstruction is what let MR's admission set drift from the identifier's.
  #
  # Absence is absence: a floor that does not apply is NULL, not -Inf.  The
  # three cases are three code paths, not one expression fed sentinels -- with
  # a sentinel, "no filter" is arithmetically indistinguishable from "a filter
  # that happens to admit everything", and only the second is testable.
  .has_effect <- !is.null(admission$effect_floor)
  .has_cons   <- !is.null(admission$consistency)
  c_cons <- if (.has_cons) admission$consistency$c_cons else NULL

  if (.has_effect && .has_cons) {
    z   <- stats::qnorm((1 + admission$consistency$p_star) / 2)
    t_g <- pmax(admission$effect_floor, c_cons + z * sdv)
    .admit <- function(bs) which(bs >= t_g)
  } else if (.has_effect) {
    t_g <- admission$effect_floor
    .admit <- function(bs) which(bs >= t_g)
  } else {
    # Unrestricted: every estimable candidate is admissible, matching an
    # identifier that applied no auxiliary selection condition.
    t_g <- NULL
    .admit <- function(bs) seq_along(bs)
  }

  # The consistency-standardized score is only defined when a consistency floor
  # exists, and only the "maxcons" rule consults it.  Ranking on it without one
  # would be the same class of mismatch this change removes, in the ranking
  # rather than the domain, so it is refused rather than silently centred at 0.
  if (identical(reselection, "maxcons") && !.has_cons) {
    stop("reselection = \"maxcons\" ranks on a consistency-standardized ",
         "statistic, but the resolved admission set has no consistency floor ",
         "(subgroup_method without a consistency screen, or ",
         "sg_focus = \"maxeff\").  MR would then rank on a quantity the ",
         "identifier never computed.", call. = FALSE)
  }
  .zcons <- function(bs) if (.has_cons) (bs - c_cons) / sdv else NULL

  t0 <- proc.time()
  Xi <- .fs_mr_multipliers(nrow(B), draws, multiplier)
  P  <- crossprod(B, Xi)                 # S x draws : D_g(b)
  beta_star <- bh + P
  sel_bias <- rep(NA_real_, draws)
  winner   <- rep(NA_integer_, draws)    # which candidate won on draw b
  for (b in seq_len(draws)) {
    bs <- beta_star[, b]
    pass <- .admit(bs)
    if (!length(pass)) next
    s <- .fs_mr_select(bs, .zcons(bs), sz, pass, reselection,
                       effect_neighborhood, selection_rule, log_scale)
    if (!is.na(s)) { sel_bias[b] <- P[s, b]; winner[b] <- s }
  }
  timing <- as.numeric((proc.time() - t0)["elapsed"])

  # BOTH BIAS TERMS ARE CONDITIONAL ON IDENTIFICATION, over the same draws.
  #
  # The defect was never the denominator on its own -- it was that the two terms
  # carried DIFFERENT ones: selection_bias over the draws that produced a
  # winner, fixed_bias over all B.  The residual then mixed differently
  # normalised quantities and mean(r) was not zero, which is the condition
  # Eq. 13's uncentered rbar2 needs in order to be Wager's centered v-hat.
  #
  # Averaging BOTH over `ok` repairs that.  The alternative repair -- both over
  # B, with D := 0 on a no-winner draw -- is equally centered, so centering does
  # not choose between them; they differ in ESTIMAND.  Conditioning is the
  # deliberate choice: the reported analysis exists only because a subgroup was
  # identified, beta(H-hat) is already a conditional estimand, and a bootstrap
  # replicate that identifies nothing contributes nothing to the full bootstrap
  # MR approximates.  It also avoids having to invent a convention for the
  # complement of an empty winner, whose complement is everyone.
  # Named ok_H, not ok: the complement block below rebinds `ok` to its own draw
  # set, and two different index sets under one name is the kind of quiet
  # aliasing this change exists to remove.
  ok_H <- which(is.finite(sel_bias))
  selection_rate <- mean(!is.na(sel_bias))
  selection_bias <- mean(sel_bias, na.rm = TRUE)   # over `ok_H`
  fixed_bias     <- mean(P[sel, ok_H])             # over `ok_H`, was all B
  fb             <- if (is.finite(fixed_bias)) fixed_bias else 0
  beta_naive <- bh[sel]
  beta_deb   <- beta_naive - selection_bias - fb
  se_wald    <- sdv[sel]

  # Infinitesimal-jackknife variance of the de-biased estimate (Eq. VInfJ_bc),
  # from the same draws: r_b = (selection_bias + fixed_bias) - D_{H*_b}(b) - D_H(b).
  # Evaluated on `ok`, the same draws both bias terms average over -- which is
  # what makes mean(r_b) identically zero there.
  r_H   <- (selection_bias + fb) - sel_bias - P[sel, ]
  # Exposure only -- reads the residual the IJ is about to consume, over the
  # same `ok_H`.  Zero by construction while both bias terms share a
  # denominator; a non-zero value means they do not, which is what F13 was.
  mean_r <- mean(r_H[ok_H])
  ijH   <- .fs_mr_ij_var(Xi, r_H, ok_H)
  se_ij <- .fs_mr_se_from_ij(ijH, se_wald)
  # "field" keeps the debiased element on the IJ interval (identical to "ij");
  # only "wald" switches it to the robust SE.
  se    <- if (ci_method == "wald") se_wald else se_ij$se

  t_cmp    <- if (log_scale) log(t_confirm) else t_confirm
  ci_lo_1s <- beta_deb - stats::qnorm(0.95) * se
  flag <- if (confirm_rule == "point") (beta_deb >= t_cmp) else (ci_lo_1s >= t_cmp)

  # ---------------------------------------------------------------------------
  # Complement subgroup (optional).  The complement is induced by the selection,
  # not chosen independently, so its bias is the perturbation of the complement
  # of the re-selected winner on each draw.  Fit complements only for candidates
  # that win across draws (plus the selected one) to keep the cost small.
  # ---------------------------------------------------------------------------
  complement <- NULL
  mean_r_c   <- NA_real_        # stays NA when no complement is fit
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
      pcc <- tryCatch(.fs_mr_pieces(df[comp_idx, , drop = FALSE], spec),
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
    # Same principle as the selected subgroup: both bias terms and the IJ share
    # one denominator.  Here that is the draws on which the complement's
    # perturbation EXISTS, which is_finite(selb_c) already identifies -- it
    # excludes two distinct events at once, and under this convention both are
    # excluded rather than one being substituted:
    #
    #   (a) the draw admitted no candidate      -- winner[b] is NA;
    #   (b) the winner's complement was not fit -- set NA at the line above.
    #
    # (b) is missing data, not a selection outcome, so it must never be
    # substituted with a value.  Conditioning makes (a) excluded too, so the
    # distinction needs no code: naming the shared set is enough.  It also
    # avoids the question of what the complement of an EMPTY winner is -- under
    # this convention that draw is excluded on both sides, and the complement is
    # always a genuine complement.
    use_c     <- which(is.finite(selb_c))
    selbias_c <- mean(selb_c[use_c])
    fixed_c   <- mean(Pc[sel, use_c])            # over `use_c`, was all B
    bnc       <- bh_c[sel]
    if (is.finite(bnc)) {
      sbc <- if (is.finite(selbias_c)) selbias_c else 0
      fcc <- if (is.finite(fixed_c)) fixed_c else 0
      bdc <- bnc - sbc - fcc
      sec <- sdv_c[sel]
      r_c   <- (sbc + fcc) - selb_c - Pc[sel, ]
      mean_r_c <- mean(r_c[use_c])          # exposure only; see mean_r above
      ijC   <- .fs_mr_ij_var(Xi, r_c, use_c)
      se_ijc <- .fs_mr_se_from_ij(ijC, sec)
      sec_used <- if (ci_method %in% c("ij", "field")) se_ijc$se else sec
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

  # ---------------------------------------------------------------------------
  # Field-calibrated interval (ci_method = "field") -- method proposal,
  # TASK_mr_field_vs_guohe_2026-09-05.  Add-only: this block runs only when
  # requested, after the main multiplier stream is fully consumed, under a
  # derived seed -- so the "ij"/"wald" paths above are untouched, RNG stream
  # included.
  #
  # Shrunk field w = beta_hat with the winner's entry replaced by the two-term
  # de-biased estimate (E1).  Gaussian-multiplier perturbations zeta = B' xi,
  # xi ~ N(0, I_n) -- no Cholesky, no explicit Sigma.  Per outer draw r:
  #   v_r = w + zeta*_r;  G_r = S(v_r);
  #   m-hat(v_r) = mean_j zeta'_{j, S(v_r + zeta'_j)}  (shared inner draws,
  #                draws with no winner skipped, the bias_sel convention);
  #   Lambda*_r  = zeta*_{r, G_r} - m-hat(v_r).
  # The interval inverts Lambda* around beta_deb (basic-bootstrap form).
  # S is the gate's own re-selection: under maxeff with no admission floor it
  # is a plain argmax (vectorized over inner draws; ties.method = "first"
  # matches which.max); any other configuration goes through .fs_mr_select
  # per draw, identically to the main loop above.
  # ---------------------------------------------------------------------------
  field <- NULL
  if (ci_method == "field") {
    t0f <- proc.time()
    if (!is.null(seed)) set.seed(as.integer(seed) + 900000L)
    Np <- nrow(B)
    Zo <- crossprod(B, matrix(stats::rnorm(Np * field_R_out), Np, field_R_out))
    Zi <- crossprod(B, matrix(stats::rnorm(Np * field_R_in), Np, field_R_in))
    w <- bh; w[sel] <- beta_deb
    fast <- identical(reselection, "maxeff") && is.null(t_g)
    sel_one <- function(v) {
      pass <- .admit(v)
      if (!length(pass)) return(NA_integer_)
      .fs_mr_select(v, .zcons(v), sz, pass, reselection,
                    effect_neighborhood, selection_rule, log_scale)
    }
    lam <- rep(NA_real_, field_R_out)
    n_in_used <- rep(NA_real_, field_R_out)
    ii <- seq_len(field_R_in)
    for (r in seq_len(field_R_out)) {
      v <- w + Zo[, r]
      G <- if (fast) which.max(v) else sel_one(v)
      if (is.na(G)) next
      if (fast) {
        win <- max.col(t(v + Zi), ties.method = "first")
        lam[r] <- Zo[G, r] - mean(Zi[cbind(win, ii)])
        n_in_used[r] <- field_R_in
      } else {
        wi <- vapply(ii, function(j) sel_one(v + Zi[, j]), integer(1))
        ok_in <- which(!is.na(wi))
        if (!length(ok_in)) next
        lam[r] <- Zo[G, r] - mean(Zi[cbind(wi[ok_in], ok_in)])
        n_in_used[r] <- length(ok_in)
      }
    }
    ok_f <- which(is.finite(lam))
    if (length(ok_f) >= 2L) {
      lf <- lam[ok_f]
      qs <- stats::quantile(lf, c(.05, .25, .50, .75, .95, .025, .975),
                            names = FALSE, type = 7)
      est2_w <- beta_deb - mean(lf)
      sd_f <- stats::sd(lf)
      field <- list(
        lambda_mean = mean(lf), lambda_sd = sd_f,
        q05 = qs[1], q25 = qs[2], q50 = qs[3], q75 = qs[4], q95 = qs[5],
        q025 = qs[6], q975 = qs[7],
        n_out_used = length(ok_f), n_in_used_mean = mean(n_in_used[ok_f]),
        est2 = to_eff(est2_w),
        lower_1s = to_eff(beta_deb - qs[5]),
        lower_2s = to_eff(beta_deb - qs[7]), upper_2s = to_eff(beta_deb - qs[6]),
        se_field = sd_f,
        lower_se = to_eff(est2_w - z975 * sd_f),
        upper_se = to_eff(est2_w + z975 * sd_f),
        R_out = as.integer(field_R_out), R_in = as.integer(field_R_in),
        seed_offset = 900000L,
        timing_seconds = as.numeric((proc.time() - t0f)["elapsed"]))

      # -- Uniform (kappa) calibration -- add-only, after the field's stream
      # is fully consumed; its own derived seed (+910000L) mirrors the field
      # block's (+900000L), so nothing above changes byte-for-byte whether or
      # not the sweep runs (TASK_mr_field_uniform_2026-09-05).
      if (isTRUE(field_uniform)) {
        p_hat_u <- tabulate(winner[!is.na(winner)],
                            nbins = length(asm$names)) / draws
        uni <- tryCatch(
          fs_mr_field_uniform(
            B = B, w = w, sel = sel, sigma_sel = sdv[sel],
            p_hat = p_hat_u, t_g = t_g,
            reselection = reselection, sz = sz,
            effect_neighborhood = effect_neighborhood,
            selection_rule = selection_rule, log_scale = log_scale,
            sdv = sdv, zcons_c = c_cons,
            seed = if (!is.null(seed)) as.integer(seed) + 910000L else NULL),
          error = function(e)
            list(note = paste("uniform sweep failed:", conditionMessage(e))))
        field$uniform <- if (is.null(uni$note)) list(
          kappa = uni$kappa, kappa_mcse = uni$kappa_mcse,
          M = uni$M, mass_covered = uni$mass_covered,
          minC1 = uni$minC1, C1 = uni$C1, C2_k1 = uni$C2_k1,
          C2_kstar = uni$C2_kstar, delta_grid = uni$delta_grid,
          n_kept = uni$n_kept,
          # The uniform two-sided interval for THIS analysis: the real
          # field's quantiles widened about its own mean (task step 6).
          lower_2u = to_eff(beta_deb - mean(lf) - uni$kappa * (qs[7] - mean(lf))),
          upper_2u = to_eff(beta_deb - mean(lf) - uni$kappa * (qs[6] - mean(lf))),
          R_rep = uni$R_rep, R_out = uni$R_out, R_in = uni$R_in,
          seed_offset = 910000L,
          timing_seconds = uni$timing_seconds) else uni
      }
    } else {
      field <- list(note = "fewer than 2 usable outer draws",
                    n_out_used = length(ok_f),
                    R_out = as.integer(field_R_out),
                    R_in = as.integer(field_R_in))
    }
  }

  out <- list(
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
    mean_r = mean_r, mean_r_c = mean_r_c,
    complement = complement,
    settings = mr_settings, harm_flag = isTRUE(flag),
    n_family = length(asm$names), n_selected = sz[sel],
    timing_seconds = timing)
  if (isTRUE(return_reselection)) {
    p_hat <- tabulate(winner[!is.na(winner)], nbins = length(asm$names)) / draws
    names(p_hat) <- asm$names
    out$reselection <- list(winner = winner, p_hat = p_hat)
  }
  if (!is.null(field)) out$field <- field
  out
}
