# test-mr-ij-residual-joint.R
#
# fs_mr_inference(ij_residual =) and field$joint
# (TASK_complement_refinements_2026-09-06).  Contracts:
#   1. add-only: with the defaults, every pre-existing element is identical to
#      a call made before the variants existed (the new elements excluded);
#      ij_residual = "winner"/"winner_floor" changes only the reported SE and
#      bounds, never the point estimates.
#   2. K = 1: se_ij_winner equals the naive SE within Monte Carlo error on both
#      blocks, two-term ~ 2x it (the 4 sigma^2 identity), winner_floor = naive.
#   3. joint: gamma in [alpha/2, alpha], achieved probability >= 1 - alpha,
#      pair bounds consistent with the marginal quantiles, Bonferroni pair
#      never tighter than the calibrated pair.

.mrj_data <- function(n = 400L, seed = 21L) {
  set.seed(seed)
  data.frame(time  = stats::rexp(n, 0.2), event = stats::rbinom(n, 1L, 0.8),
             treat = stats::rbinom(n, 1L, 0.5),
             x1 = stats::rnorm(n), x2 = stats::rnorm(n))
}
.mrj_spec <- list(outcome_type = "survival", effect_measure = "HR",
                  treat.name = "treat", outcome.name = "time",
                  event.name = "event", offset.name = NULL,
                  adjust_covariates = NULL, adverse_outcome = TRUE)
.mrj_call <- function(df, cands, sel, ..., draws = 600L) {
  fs_mr_inference(df = df, candidates = cands, spec = .mrj_spec,
                  selected_members = cands[[sel]],
                  admission = list(effect_floor = NULL, consistency = NULL),
                  reselection = "maxeff", draws = draws,
                  multiplier = "poisson", seed = 777L, ...)
}
.new_el <- c("se_ij_two_term", "se_ij_winner", "ij_source_winner",
             "se_ij_winner_floor", "ij_source_winner_floor",
             "lower_w", "upper_w", "lower_1s_w", "upper_1s_w",
             "lower_wf", "upper_wf", "lower_1s_wf", "upper_1s_wf")
.strip_new <- function(x) {
  x$timing_seconds <- NULL; x$ij_residual <- NULL
  x$debiased[.new_el] <- NULL
  if (is.list(x$complement$debiased)) x$complement$debiased[.new_el] <- NULL
  if (is.list(x$field)) { x$field$timing_seconds <- NULL; x$field$joint <- NULL
    if (is.list(x$field$complement)) x$field$complement$timing_seconds <- NULL }
  x
}

test_that("ij_residual variants are add-only and leave point estimates alone", {
  df <- .mrj_data()
  cands <- list(g1 = which(df$x1 > 0), g2 = which(df$x2 > 0),
                g3 = which(df$x1 + df$x2 > 0.3))
  bh <- vapply(cands, function(ix) unname(coef(survival::coxph(
    survival::Surv(time, event) ~ treat, data = df[ix, ]))), numeric(1))
  sel <- names(cands)[which.max(bh)]
  d0 <- .mrj_call(df, cands, sel, ci_method = "ij", include_complement = TRUE)
  d1 <- .mrj_call(df, cands, sel, ci_method = "ij", include_complement = TRUE,
                  ij_residual = "two_term")
  w  <- .mrj_call(df, cands, sel, ci_method = "ij", include_complement = TRUE,
                  ij_residual = "winner")
  wf <- .mrj_call(df, cands, sel, ci_method = "ij", include_complement = TRUE,
                  ij_residual = "winner_floor")
  expect_identical(.strip_new(d0), .strip_new(d1))
  expect_identical(d0$ij_residual, "two_term")
  # Default reports the two-term SE; the variants report theirs.
  expect_equal(d0$debiased$se_ij, d0$debiased$se_ij_two_term)
  expect_equal(w$debiased$se_ij, d0$debiased$se_ij_winner)
  expect_equal(wf$debiased$se_ij, d0$debiased$se_ij_winner_floor)
  expect_equal(w$complement$debiased$se_ij, d0$complement$debiased$se_ij_winner)
  expect_equal(wf$complement$debiased$se_ij, d0$complement$debiased$se_ij_winner_floor)
  expect_gte(wf$debiased$se_ij, d0$debiased$se_wald - 1e-12)
  expect_gte(wf$complement$debiased$se_ij, d0$complement$debiased$se_wald - 1e-12)
  # Point estimates and bias terms untouched by the variant.
  for (o in list(w, wf)) {
    expect_identical(o$debiased$est, d0$debiased$est)
    expect_identical(o$complement$debiased$est, d0$complement$debiased$est)
    expect_identical(o$selection_bias, d0$selection_bias)
    expect_identical(o$fixed_bias, d0$fixed_bias)
    expect_identical(o$naive, d0$naive)
  }
  # Bounds follow the reported SE.
  expect_equal(w$debiased$lower, d0$debiased$lower_w)
  expect_equal(w$debiased$upper, d0$debiased$upper_w)
  expect_equal(wf$complement$debiased$upper, d0$complement$debiased$upper_wf)
  # The harm block's exposed side is the lower bound: no upper_1s_w/_wf there.
  expect_true(all(setdiff(.new_el, c("upper_1s_w", "upper_1s_wf")) %in% names(d0$debiased)))
  expect_false(any(c("upper_1s_w", "upper_1s_wf") %in% names(d0$debiased)))
  expect_true(all(.new_el %in% names(d0$complement$debiased)))
})

test_that("K = 1: winner-only SE is the naive SE; two-term is ~2x; floor equals naive", {
  df <- .mrj_data(n = 500L, seed = 22L)
  cands <- list(S1 = which(df$x1 > 0))
  o <- .mrj_call(df, cands, "S1", ci_method = "ij", include_complement = TRUE,
                 draws = 4000L)
  for (blk in list(o$debiased, o$complement$debiased)) {
    expect_lt(abs(blk$se_ij_winner / blk$se_wald - 1), 0.10)
    expect_lt(abs(blk$se_ij_two_term / (2 * blk$se_wald) - 1), 0.12)
    expect_equal(blk$se_ij_winner_floor, max(blk$se_ij_winner, blk$se_wald))
  }
})

test_that("field$joint: gamma in range, achieved probability >= 1 - alpha, pairs consistent", {
  df <- .mrj_data(n = 400L, seed = 23L)
  cands <- list(g1 = which(df$x1 > 0), g2 = which(df$x2 > 0),
                g3 = which(df$x1 + df$x2 > 0.3))
  bh <- vapply(cands, function(ix) unname(coef(survival::coxph(
    survival::Surv(time, event) ~ treat, data = df[ix, ]))), numeric(1))
  sel <- names(cands)[which.max(bh)]
  off <- .mrj_call(df, cands, sel, ci_method = "field", include_complement = TRUE,
                   field_R_out = 400L, field_R_in = 100L, field_complement = TRUE)
  j <- off$field$joint
  expect_true(is.list(j)); expect_null(j$note)
  expect_true(j$gamma >= 0.025 - 1e-12 && j$gamma <= 0.05 + 1e-12)
  expect_gte(j$joint_prob, 0.95 - 1e-12)
  expect_equal(j$bonf_gamma, 0.025)
  expect_gte(j$bonf_joint_prob, 0.95 - 1e-12)
  # Calibrated pair is never tighter than... never wider than Bonferroni's.
  expect_gte(j$lower_H, j$bonf_lower_H - 1e-12)
  expect_lte(j$upper_Hc, j$bonf_upper_Hc + 1e-12)
  # Marginal one-sided bounds untouched and consistent with the pair at gamma = 0.05.
  expect_true(is.finite(off$field$lower_1s) && is.finite(off$field$complement$upper_1s))
  expect_lte(j$lower_H, off$field$lower_1s + 1e-12)
  expect_gte(j$upper_Hc, off$field$complement$upper_1s - 1e-12)
  expect_true(abs(j$corr) <= 1)
  expect_equal(j$n_joint_draws, off$field$complement$n_out_used)
  # The complement field list is unchanged by the joint addition.
  expect_false("joint" %in% names(off$field$complement))
  # No joint without the complement field.  field_complement is set explicitly
  # here: it has defaulted to TRUE since TASK_cert20_2026-09-08 Part D, so
  # omitting it no longer produces the "without" condition this checks.
  no <- .mrj_call(df, cands, sel, ci_method = "field", include_complement = TRUE,
                  field_R_out = 400L, field_R_in = 100L, field_complement = FALSE)
  expect_null(no$field$joint)
})
