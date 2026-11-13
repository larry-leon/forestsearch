# Acceptance tests for B1 (CRAN preparation, part 5): the standard error of
# the binary risk difference in the covariate-adjusted logistic tier.
#
# When the identity-link binomial fit (tier 1 of .estimate_rd) fails, the
# risk difference is estimated by logistic g-computation (tier 2,
# method_used "logistic_gcomp_adj" when adjustment covariates are given):
#
#   rd(b) = mean(plogis(X1 %*% b)) - mean(plogis(X0 %*% b)),
#
# with X1, X0 the design matrices of every subject with treatment set to 1
# and to 0, and b the logistic coefficients.  Its delta-method standard
# error is sqrt(g' V g), with V the covariance matrix of the logistic fit and
# g the gradient of rd(b):
#
#   g = colMeans(p1 * (1 - p1) * X1) - colMeans(p0 * (1 - p0) * X0),
#
# p1, p0 the fitted probabilities under treatment and under control.  The
# old adjusted standard error was abs(mean(p1 * (1 - p1))) * se(treatment
# coefficient); it is kept only when sqrt(g' V g) is not a finite,
# non-negative number.
#
# Two deterministic routes reach tier 2 with covariates at this tree:
#   - a factor covariate with three levels: tier 1 supplies one start value
#     per covariate NAME (c(p0, p1 - p0, rep(0, n_adj))), the design has one
#     column per non-reference LEVEL, glm() stops on the start length, and
#     the slice goes to tier 2 with no degenerate arm;
#   - an arm at 100% events: the start value puts that arm's fitted
#     probability at exactly 1, which binomial() rejects, so tier 1 stops.
# Each test through make_effect_estimator() asserts method_used first, as the
# premise that the slice reached tier 2.
#
# Every expected value is worked out by hand or recomputed in the test by
# central finite differences on a separately fitted glm(); none is taken
# from the code under test.

# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

# 60 patients, 30 per arm; site a/b/c in turn; x on a 10-point grid in
# [-1.5, 1.5]; y a fixed pattern.  Every treat x site cell has both outcomes
# (control 10/30 events, treated 16/30), so the logistic fit is regular.
.rdadj_fac <- function() {
  i <- 1:60
  treat <- rep(c(0L, 1L), each = 30)
  site  <- factor(rep(c("a", "b", "c"), length.out = 60))
  x     <- rep(seq(-1.5, 1.5, length.out = 10), length.out = 60)
  y     <- as.integer((i * 7) %% 10 < 3 + 2 * treat + (site == "c"))
  data.frame(treat = treat, site = site, x = x, y = y)
}

.rdadj_fn <- function(adj) {
  make_effect_estimator(outcome_type = "binary", treat.name = "treat",
                        outcome.name = "y", effect_measure = "RD",
                        adjust_covariates = adj)
}

# Delta-method standard error recomputed independently: a separate glm(),
# design matrices built from an explicit right-hand side, the gradient of
# rd(b) by central finite differences, and V = vcov(fit, complete = FALSE)
# (estimable coefficients only; aliased columns are dropped from X1, X0).
#
# Step and tolerance, fixed before any run.  Central differences have
# truncation error h^2 / 6 * |f'''| and rounding error about u * |f| / h
# (u = 1.1e-16, |f| <= 1).  Along any coefficient, |f'''| <= 2 * (1/8) *
# max|x|^3, since |plogis'''| <= 1/8; the design columns here are 0/1 or in
# [-1.5, 3] (x2 = 2 * x reaches 3 but is aliased and dropped), so
# |f'''| <= 2 * 0.125 * 3.375 = 0.84.  With h = 1e-5: truncation <= 1.4e-11
# and rounding about 1e-11, against gradient entries of order 1e-2 to 1e-1,
# a relative error of at most about 1e-9 in g and in sqrt(g' V g).
# Tolerance 1e-6 (relative, testthat's expect_equal) leaves two to three
# orders of headroom.
.rdadj_h   <- 1e-5
.rdadj_tol <- 1e-6

.rdadj_fd_se <- function(data, rhs) {
  fit <- stats::glm(stats::update(rhs, y ~ .), data = data,
                    family = stats::binomial(link = "logit"))
  b   <- stats::coef(fit)
  est <- !is.na(b)
  d1 <- data; d1$treat <- 1L
  d0 <- data; d0$treat <- 0L
  X1 <- stats::model.matrix(rhs, d1)[, est, drop = FALSE]
  X0 <- stats::model.matrix(rhs, d0)[, est, drop = FALSE]
  rd <- function(beta) {
    mean(stats::plogis(X1 %*% beta)) - mean(stats::plogis(X0 %*% beta))
  }
  bb <- b[est]
  g <- vapply(seq_along(bb), function(j) {
    e <- numeric(length(bb)); e[j] <- .rdadj_h
    (rd(bb + e) - rd(bb - e)) / (2 * .rdadj_h)
  }, numeric(1))
  V <- stats::vcov(fit, complete = FALSE)
  sqrt(drop(t(g) %*% V %*% g))
}

# ---------------------------------------------------------------------------
# (a) Adjusted logistic fit: SE equals the finite-difference delta method
# ---------------------------------------------------------------------------

test_that("adjusted RD, logistic tier: SE is the delta-method sqrt(g' V g)", {
  d   <- .rdadj_fac()
  res <- .rdadj_fn(c("site", "x"))(d)

  # Premise: the slice reached the adjusted logistic tier (tier 1 stopped on
  # the start length: 4 start values for the 5 columns (Intercept), treat,
  # siteb, sitec, x).
  expect_identical(res$method_used, "logistic_gcomp_adj")
  expect_true(res$converged)

  expect_equal(res$se, .rdadj_fd_se(d, ~ treat + site + x),
               tolerance = .rdadj_tol)
})

test_that("adjusted RD, logistic tier, aliased covariate: estimable columns only", {
  # x2 = 2 * x is aliased with x: its coefficient is NA and vcov() has an NA
  # row and column for it.  The delta method uses the estimable coefficients
  # only, so the SE is the same finite-difference value on the reduced model.
  d <- .rdadj_fac()
  d$x2 <- 2 * d$x
  res <- .rdadj_fn(c("site", "x", "x2"))(d)

  expect_identical(res$method_used, "logistic_gcomp_adj")
  expect_equal(res$se, .rdadj_fd_se(d, ~ treat + site + x + x2),
               tolerance = .rdadj_tol)
  # Dropping x2 from the model gives the same fit, hence the same SE.
  expect_equal(res$se, .rdadj_fd_se(d, ~ treat + site + x),
               tolerance = .rdadj_tol)
})

# ---------------------------------------------------------------------------
# (b) No covariates: the same calculation is the two-sample binomial SE
# ---------------------------------------------------------------------------

test_that("g' V g with no covariates equals p1 (1 - p1) / n1 + p0 (1 - p0) / n0", {
  # Control 12 / 40 = 0.3 events; treated 15 / 30 = 0.5 events.
  d <- data.frame(
    treat = rep(c(0L, 1L), c(40L, 30L)),
    y     = c(rep(1L, 12L), rep(0L, 28L), rep(1L, 15L), rep(0L, 15L)))
  fit <- stats::glm(y ~ treat, data = d, family = stats::binomial())
  nd0 <- d; nd0$treat <- 0L
  nd1 <- d; nd1$treat <- 1L
  p0 <- stats::predict(fit, newdata = nd0, type = "response")
  p1 <- stats::predict(fit, newdata = nd1, type = "response")

  # By hand.  The two-arm model is saturated: p0, p1 are the arm proportions
  # for every subject, X0 rows (1, 0), X1 rows (1, 1), q = 1 - p, so
  #   g = (p1 q1 - p0 q0, p1 q1),
  #   V = solve([[n0 p0 q0 + n1 p1 q1, n1 p1 q1], [n1 p1 q1, n1 p1 q1]])
  #     : V_aa = 1 / (n0 p0 q0), V_at = -1 / (n0 p0 q0),
  #       V_tt = 1 / (n0 p0 q0) + 1 / (n1 p1 q1),
  #   g' V g = (g_a - g_t)^2 / (n0 p0 q0) + g_t^2 / (n1 p1 q1)
  #          = p0 q0 / n0 + p1 q1 / n1.
  # Here: 0.3 * 0.7 / 40 + 0.5 * 0.5 / 30 = 0.00525 + 0.0083333 = 0.0135833,
  #       sqrt(0.0135833) = 0.1165476.
  # Tolerance 1e-6, fixed before any run: the fitted probabilities equal the
  # arm proportions, and V = vcov(fit) the inverse information, only to
  # glm()'s convergence (relative deviance change below 1e-8); vcov() also
  # uses the IRLS weights from before glm()'s final update.  Measured after
  # the run: a relative difference of 1.9e-7.
  expect_equal(.fs_rd_gcomp_se(fit, nd0, nd1, p0, p1),
               sqrt(0.3 * 0.7 / 40 + 0.5 * 0.5 / 30), tolerance = 1e-6)
  expect_equal(round(.fs_rd_gcomp_se(fit, nd0, nd1, p0, p1), 5), 0.11655)
})

# ---------------------------------------------------------------------------
# (c) Through make_effect_estimator() on a slice with an arm at 100%
# ---------------------------------------------------------------------------

test_that("adjusted RD, treated arm at 100%: SE by hand, 40 per arm", {
  # Binary covariate z.  Control: z = 0, 20 patients, 4 events (0.2);
  # z = 1, 20 patients, 10 events (0.5).  Treated: z = 0, 10 patients;
  # z = 1, 30 patients; all 40 have the event.
  d <- data.frame(
    treat = rep(c(0L, 1L), each = 40L),
    z     = c(rep(0L, 20L), rep(1L, 20L), rep(0L, 10L), rep(1L, 30L)),
    y     = c(rep(1L, 4L), rep(0L, 16L), rep(1L, 10L), rep(0L, 10L),
              rep(1L, 40L)))
  res <- .rdadj_fn("z")(d)

  # Premise: tier 1 stopped (its start value puts the treated arm at
  # probability exactly 1) and the slice reached the adjusted logistic tier.
  expect_identical(res$method_used, "logistic_gcomp_adj")
  expect_true(res$converged)

  # By hand.  With the treated arm at 100%, the logistic MLE sends the
  # treatment coefficient towards +Inf (glm() stops at a large finite value),
  # every p1 tends to 1 and p1 (1 - p1) to 0, so the treated arm contributes
  # nothing to g' V g.  The score equations for the intercept and z then hold
  # on the control arm alone, where the model in z is saturated: p0(z) is the
  # control proportion at z, with variance p0z (1 - p0z) / n0z.  The estimate
  # averages p0(z) over all 80 patients with weights w_z, the share of
  # patients at z: w_0 = (20 + 10) / 80 = 0.375, w_1 = (20 + 30) / 80 = 0.625.
  #   SE^2 = w_0^2 p00 q00 / n00 + w_1^2 p01 q01 / n01
  #        = 0.140625 * 0.2 * 0.8 / 20 + 0.390625 * 0.5 * 0.5 / 20
  #        = 0.001125 + 0.0048828125 = 0.0060078125,
  #   SE   = sqrt(0.0060078125) = 0.0775101.
  # Estimate: 1 - (0.375 * 0.2 + 0.625 * 0.5) = 1 - 0.3875 = 0.6125.
  #
  # Tolerance 1e-6, fixed before any run: glm() stops once the relative
  # change in deviance is below 1e-8, which leaves the treated arm's fitted
  # probabilities of order 1e-8 away from 1 (as in test-rd-standard-error.R).
  # The terms this adds to g' V g are of that order times entries of g and V
  # of order 0.1 to 1 (the treatment coefficient's large variance enters
  # multiplied by its gradient squared, of order 1e-16), i.e. about 1e-9
  # against SE^2 = 6e-3: a relative error of at most about 1e-7 in the SE.
  expect_equal(res$se,
               sqrt(0.375^2 * 0.2 * 0.8 / 20 + 0.625^2 * 0.5 * 0.5 / 20),
               tolerance = 1e-6)
  expect_equal(round(res$se, 6), 0.07751)
  expect_equal(res$estimate, 0.6125, tolerance = 1e-6)
})

# ---------------------------------------------------------------------------
# Safety rule: a value that is not finite and non-negative keeps the old form
# ---------------------------------------------------------------------------

test_that("adjusted RD: a non-finite or negative delta-method value keeps the old SE", {
  d <- .rdadj_fac()
  # The old expression, by its definition, on a separately fitted glm():
  # abs(mean(p1 * (1 - p1))) * se(treatment coefficient).
  fit <- stats::glm(y ~ treat + site + x, data = d, family = stats::binomial())
  nd1 <- d; nd1$treat <- 1L
  p1  <- stats::predict(fit, newdata = nd1, type = "response")
  old <- abs(mean(p1 * (1 - p1))) * sqrt(diag(stats::vcov(fit)))[["treat"]]

  # The helper is replaced by one returning each kind of bad value; the
  # closure must then return the old expression.  Tolerance 1e-12, fixed
  # before any run: the same model on the same data with the same arithmetic,
  # so only rounding can differ.
  for (bad in list(NA_real_, NaN, Inf, -0.5)) {
    local_mocked_bindings(.fs_rd_gcomp_se = function(...) bad)
    res <- .rdadj_fn(c("site", "x"))(d)
    expect_identical(res$method_used, "logistic_gcomp_adj")
    expect_equal(res$se, old, tolerance = 1e-12)
  }
})

test_that("the delta-method helper never errors: failures return NA", {
  d   <- .rdadj_fac()
  fit <- stats::glm(y ~ treat + site + x, data = d, family = stats::binomial())
  # Not a fitted model.
  expect_identical(.fs_rd_gcomp_se(NULL, d, d, 0.5, 0.5), NA_real_)
  # Probabilities of the wrong length.
  expect_identical(.fs_rd_gcomp_se(fit, d, d, c(0.5, 0.5), c(0.5, 0.5)),
                   NA_real_)
  # A missing covariate value: the estimate is NA and the closure returns
  # without an error, with the old expression (also NA) as the SE.
  d$x[3] <- NA
  expect_no_error(res <- .rdadj_fn(c("site", "x"))(d))
  expect_identical(res$method_used, "logistic_gcomp_adj")
  expect_true(is.na(res$estimate))
  expect_true(is.na(res$se))
})

# ---------------------------------------------------------------------------
# Not estimable: the treatment coefficient is aliased (repair round)
# ---------------------------------------------------------------------------

test_that("adjusted RD, treatment effect not estimable: the SE stays NA, not 0", {
  # Control: 30 patients, z = 0, 1 in turn, an event at every third patient
  # (5 events at z = 0 and 5 at z = 1), so the control fit in z is regular.
  ctl <- data.frame(treat = 0L, z = rep(0:1, 15L),
                    y = as.integer(1:30 %% 3L == 0L))
  # Two slices on which no data inform the treatment coefficient:
  #   - no treated patient at all (one arm empty);
  #   - 10 treated patients whose outcomes are all missing (glm() drops
  #     them, so the fitted rows hold the control arm only).
  slices <- list(
    empty_arm   = ctl,
    trt_all_na  = rbind(ctl, data.frame(treat = 1L, z = rep(0:1, 5L),
                                        y = NA_integer_)))

  for (nm in names(slices)) {
    d <- slices[[nm]]
    # Premise 1: in a separately fitted glm() the treatment column is
    # constant on the fitted rows, so its coefficient is aliased (NA).
    fit <- stats::glm(y ~ treat + z, data = d, family = stats::binomial())
    expect_true(is.na(stats::coef(fit)[["treat"]]), info = nm)

    # Expected value, by definition: the risk difference is not estimable,
    # so it has no standard error.  The old expression gives NA there too,
    # since vcov() has an NA row and column for the aliased coefficient
    # (abs(avg_deriv) * NA = NA).  The delta method on the estimable columns
    # alone would give 0 (X1 and X0 then differ only in the dropped
    # treatment column, so g = 0), which is wrong.
    #
    # predict() on a rank-deficient fit may warn about non-estimable rows
    # (in the closure and here); that warning is not what is tested.
    nd0 <- d; nd0$treat <- 0L
    nd1 <- d; nd1$treat <- 1L
    p0 <- suppressWarnings(stats::predict(fit, newdata = nd0, type = "response"))
    p1 <- suppressWarnings(stats::predict(fit, newdata = nd1, type = "response"))
    expect_identical(.fs_rd_gcomp_se(fit, nd0, nd1, p0, p1), NA_real_,
                     info = nm)

    # Through make_effect_estimator(): premise 2, the slice reached the
    # adjusted logistic tier; then the SE is NA.
    res <- suppressWarnings(.rdadj_fn("z")(d))
    expect_identical(res$method_used, "logistic_gcomp_adj", info = nm)
    expect_true(is.na(res$se), info = nm)
  }
})
