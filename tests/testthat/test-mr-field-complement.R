# test-mr-field-complement.R
#
# fs_mr_inference(field_complement = TRUE): the complement's field block
# (TASK_mr_field_complement_2026-09-06).  Contracts checked:
#   1. add-only: with the flag on, every output other than field$complement
#      is identical to the flag-off call (the harm field's draws are shared,
#      not redrawn); the flag off reproduces the omitted-argument call.
#   2. K = 1: Lambda*c is centred with sd equal to the complement's naive SE,
#      so upper_1s is the naive one-sided upper bound (within Monte Carlo error).
#   3. bound identities at 1e-12 on the returned quantities.
#   4. fs_sim_bias_coverage(side = "upper") reads fld_Hc_up1s and scores
#      truth <= upper bound; the default side is unchanged.

.mfc_data <- function(n = 300L, seed = 11L) {
  set.seed(seed)
  data.frame(time  = stats::rexp(n, 0.2),
             event = stats::rbinom(n, 1L, 0.8),
             treat = stats::rbinom(n, 1L, 0.5),
             x1 = stats::rnorm(n), x2 = stats::rnorm(n))
}
.mfc_spec <- list(outcome_type = "survival", effect_measure = "HR",
                  treat.name = "treat", outcome.name = "time",
                  event.name = "event", offset.name = NULL,
                  adjust_covariates = NULL, adverse_outcome = TRUE)
.mfc_call <- function(df, cands, sel, ...) {
  fs_mr_inference(df = df, candidates = cands, spec = .mfc_spec,
                  selected_members = cands[[sel]],
                  admission = list(effect_floor = NULL, consistency = NULL),
                  reselection = "maxeff", draws = 400L,
                  multiplier = "poisson", seed = 4242L, ...)
}
.strip_t <- function(x) {
  x$timing_seconds <- NULL
  if (is.list(x$field)) { x$field$timing_seconds <- NULL
                          x$field$complement <- NULL }
  x
}

test_that("field_complement is add-only: everything else is byte-identical", {
  df <- .mfc_data()
  cands <- list(g1 = which(df$x1 > 0), g2 = which(df$x2 > 0),
                g3 = which(df$x1 + df$x2 > 0.3))
  # Selected: the candidate with the largest naive log-HR.
  bh <- vapply(cands, function(ix) unname(coef(survival::coxph(
    survival::Surv(time, event) ~ treat, data = df[ix, ]))), numeric(1))
  sel <- names(cands)[which.max(bh)]
  off  <- .mfc_call(df, cands, sel, ci_method = "field",
                    include_complement = TRUE, field_R_out = 300L,
                    field_R_in = 100L)
  off0 <- .mfc_call(df, cands, sel, ci_method = "field",
                    include_complement = TRUE, field_R_out = 300L,
                    field_R_in = 100L, field_complement = FALSE)
  on   <- .mfc_call(df, cands, sel, ci_method = "field",
                    include_complement = TRUE, field_R_out = 300L,
                    field_R_in = 100L, field_complement = TRUE)
  expect_null(off$field$complement)
  expect_identical(.strip_t(off), .strip_t(off0))
  expect_identical(.strip_t(on), .strip_t(off))
  expect_identical(names(on), names(off))
  # The complement field itself.
  fc <- on$field$complement
  expect_true(is.list(fc)); expect_null(fc$note)
  expect_true(all(c("est2", "upper_1s", "lower_1s", "lower_2s", "upper_2s",
                    "lower_se", "upper_se", "se_field", "lambda_mean",
                    "q05", "q25", "q50", "q75", "q95", "q025", "q975",
                    "n_out_used", "n_in_used_mean", "n_complement_fits",
                    "n_new_fits", "share_draws_new_fit",
                    "timing_seconds") %in% names(fc)))
  expect_true(is.finite(fc$upper_1s) && fc$upper_1s > 0)
  expect_gte(fc$upper_2s, fc$lower_2s)
  expect_gte(fc$upper_1s, fc$lower_1s)
  expect_gte(fc$n_complement_fits, 1L)
  expect_true(fc$share_draws_new_fit >= 0 && fc$share_draws_new_fit <= 1)
  # Bound identities (task 1b), at 1e-12 on the log scale.
  bt <- log(on$complement$debiased$est)
  expect_equal(log(fc$upper_1s), bt - fc$q05, tolerance = 1e-12)
  expect_equal(log(fc$lower_1s), bt - fc$q95, tolerance = 1e-12)
  expect_equal(log(fc$lower_2s), bt - fc$q975, tolerance = 1e-12)
  expect_equal(log(fc$upper_2s), bt - fc$q025, tolerance = 1e-12)
  expect_equal(log(fc$est2), bt - fc$lambda_mean, tolerance = 1e-12)
  expect_equal(log(fc$upper_se), log(fc$est2) + qnorm(0.975) * fc$se_field,
               tolerance = 1e-12)
  expect_equal(log(fc$lower_se), log(fc$est2) - qnorm(0.975) * fc$se_field,
               tolerance = 1e-12)
  # Off the field path, or without the complement, the flag is inert.
  ij <- .mfc_call(df, cands, sel, ci_method = "ij", include_complement = TRUE,
                  field_complement = TRUE)
  expect_null(ij$field)
  nc <- .mfc_call(df, cands, sel, ci_method = "field",
                  include_complement = FALSE, field_R_out = 300L,
                  field_R_in = 100L, field_complement = TRUE)
  expect_null(nc$complement); expect_null(nc$field$complement)
})

test_that("K = 1: the complement field is the naive one-sided upper bound", {
  df <- .mfc_data(n = 400L, seed = 12L)
  cands <- list(S1 = which(df$x1 > 0))
  on <- .mfc_call(df, cands, "S1", ci_method = "field",
                  include_complement = TRUE, field_R_out = 4000L,
                  field_R_in = 200L, field_complement = TRUE)
  fc <- on$field$complement
  sec <- on$complement$debiased$se_wald          # complement's naive SE
  expect_lt(abs(fc$lambda_mean), 4 * sec * sqrt(1 / 4000 + 1 / 200))
  expect_lt(abs(fc$lambda_sd / sec - 1), 0.05)
  expect_lt(abs(-fc$q05 / (qnorm(0.95) * sec) - 1), 0.06)
  # upper_1s vs the naive one-sided upper bound around beta-tilde^c.
  bt <- log(on$complement$debiased$est)
  expect_lt(abs(log(fc$upper_1s) - (bt + qnorm(0.95) * sec)), 0.06 * sec + 1e-12)
  expect_identical(fc$n_complement_fits, 1L)
  expect_identical(fc$n_new_fits, 0L)
  expect_equal(fc$share_draws_new_fit, 0)
})

test_that("fs_sim_bias_coverage(side = 'upper') scores the upper bound", {
  set.seed(3)
  n <- 200L
  tr <- exp(stats::rnorm(n, 0, 0.2))
  e  <- tr * exp(stats::rnorm(n, 0, 0.2))
  se <- rep(0.2, n)
  res <- data.frame(
    detected = 1L, betaHhat_Hc = tr, betaHhat_H = tr,
    nv_Hc_est = e, nv_Hc_lo = e * exp(-1.96 * se), nv_Hc_hi = e * exp(1.96 * se),
    nv_Hc_se = se,
    mr_Hc_est = e, mr_Hc_lo = e * exp(-1.96 * se), mr_Hc_hi = e * exp(1.96 * se),
    mr_Hc_se_ij = se,
    fld_Hc_est2 = e, fld_Hc_lo2s = e * exp(-1.96 * se),
    fld_Hc_hi2s = e * exp(1.96 * se), fld_Hc_se = se,
    fld_Hc_lo1s = e * exp(-1.645 * se), fld_Hc_up1s = e * exp(1.645 * se))
  up <- fs_sim_bias_coverage(res, block = "Hc", side = "upper")
  expect_setequal(up$estimator, c("naive", "mr", "fld"))
  f <- up[up$estimator == "fld", ]
  expect_equal(f$cov1, mean(tr <= res$fld_Hc_up1s))
  expect_equal(f$cov1_ref, pnorm(qnorm(0.95) * f$r + f$b))
  m <- up[up$estimator == "mr", ]
  expect_equal(m$cov1, mean(tr <= exp(log(e) + qnorm(0.95) * se)))
  lo <- fs_sim_bias_coverage(res, block = "Hc", side = "lower")
  fl <- lo[lo$estimator == "fld", ]
  expect_equal(fl$cov1, mean(tr >= res$fld_Hc_lo1s))
  expect_equal(fl$cov1_ref, pnorm(qnorm(0.95) * fl$r - fl$b))
  expect_identical(names(lo), names(up))
  # Without complement field columns the fld row is dropped with a message.
  res2 <- res[, !grepl("^fld_", names(res))]
  expect_message(d <- fs_sim_bias_coverage(res2, block = "Hc"), "dropping")
  expect_setequal(d$estimator, c("naive", "mr"))
})
