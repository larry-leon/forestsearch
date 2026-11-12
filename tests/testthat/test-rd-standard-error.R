# Acceptance tests for F13 (audit item A-055): the standard error of the
# binary risk difference in the unadjusted, unweighted logistic tier.
#
# When one arm is at 0% or 100% events, the identity-link binomial fit
# (tier 1 of .estimate_rd) cannot start, and the slice is estimated by
# logistic g-computation (tier 2, method_used "logistic_margins").  For the
# unadjusted two-arm model, which is saturated, the arm means of the fitted
# probabilities are the observed arm proportions p0 and p1, so the exact
# delta-method standard error of mean(pred1) - mean(pred0) is the two-sample
# binomial (Wald) standard error
#
#   se = sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0),
#
# with n0, n1 the arm sizes.  The estimate (the difference of the arm
# proportions), `converged` (TRUE) and `method_used` ("logistic_margins")
# do not change.  The slice must not be routed to the raw-proportion
# tier 3 (method_used "raw_means", converged = FALSE).
#
# Every expected value below is worked out by hand from the formula above.

# Two arms of 40 patients, with a given number of events in each arm.
.rdse_slice <- function(e0, e1, n = 40L) {
  data.frame(
    treat = rep(c(0L, 1L), each = n),
    y     = c(rep(1L, e0), rep(0L, n - e0), rep(1L, e1), rep(0L, n - e1))
  )
}

.rdse_fn <- function() {
  make_effect_estimator(outcome_type = "binary", treat.name = "treat",
                        outcome.name = "y", effect_measure = "RD")
}

# Tolerance for the ESTIMATE only, fixed before any run: under the
# quasi-separation of a 0% or 100% arm, glm() stops once the relative change
# in deviance is below its epsilon (1e-8), which leaves the degenerate arm's
# fitted probability of order 1e-8 away from 0 or 1.  The g-computation
# estimate therefore matches the raw difference of proportions to about
# 1e-8; 1e-6 bounds that with headroom.  The standard error is compared at
# the default tolerance, since the formula uses the observed proportions.
.rdse_est_tol <- 1e-6

test_that("RD, control 30% and treated 100%, 40 per arm: Wald SE", {
  # Control: 12 / 40 = 0.3 events; treated: 40 / 40 = 1.0 events.
  res <- .rdse_fn()(.rdse_slice(e0 = 12L, e1 = 40L))

  # By hand: p0 (1 - p0) / n0 = 0.3 * 0.7 / 40 = 0.21 / 40 = 0.00525;
  #          p1 (1 - p1) / n1 = 1.0 * 0.0 / 40 = 0;
  #          se = sqrt(0.00525) = 0.0724569 (0.0724569^2 = 0.0052500).
  expect_equal(res$se, sqrt(0.3 * 0.7 / 40 + 1.0 * 0.0 / 40))
  expect_equal(round(res$se, 5), 0.07246)

  # Estimate: difference of the arm proportions, 1.0 - 0.3 = 0.7.
  expect_equal(res$estimate, 1.0 - 0.3, tolerance = .rdse_est_tol)

  # converged stays TRUE and the method label stays the unadjusted tier-2
  # label (glm_effect_estimators.R, .estimate_rd tier 2), so the slice is
  # neither routed to tier 3 ("raw_means", converged = FALSE) nor relabelled.
  expect_true(res$converged)
  expect_identical(res$method_used, "logistic_margins")
  expect_identical(res$measure, "RD")
  expect_identical(c(res$n0, res$n1), c(40L, 40L))
})

test_that("RD, control 0% and treated 80%, 40 per arm: Wald SE", {
  # Control: 0 / 40 = 0.0 events; treated: 32 / 40 = 0.8 events.
  res <- .rdse_fn()(.rdse_slice(e0 = 0L, e1 = 32L))

  # By hand: p0 (1 - p0) / n0 = 0.0 * 1.0 / 40 = 0;
  #          p1 (1 - p1) / n1 = 0.8 * 0.2 / 40 = 0.16 / 40 = 0.004;
  #          se = sqrt(0.004) = 0.0632456 (0.0632456^2 = 0.0040000).
  expect_equal(res$se, sqrt(0.8 * 0.2 / 40 + 0.0 * 1.0 / 40))
  expect_equal(round(res$se, 5), 0.06325)

  # Estimate: difference of the arm proportions, 0.8 - 0.0 = 0.8.
  expect_equal(res$estimate, 0.8 - 0.0, tolerance = .rdse_est_tol)

  expect_true(res$converged)
  expect_identical(res$method_used, "logistic_margins")
  expect_identical(res$measure, "RD")
  expect_identical(c(res$n0, res$n1), c(40L, 40L))
})
