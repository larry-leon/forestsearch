# =============================================================================
# The lm() influence path: .dfbeta_glm() and everything downstream of it
# =============================================================================
# REGRESSION TEST for a live defect.
#
# .dfbeta_glm() read `w <- fit$weights` unconditionally.  For a glm that is
# always a vector (the IRLS working weights); for an UNWEIGHTED lm() it is
# NULL, so `crossprod(X, w * X)` was non-conformable and the function errored.
#
# effect_measure = "MD" is the only measure that fits with lm(), so the whole
# continuous-outcome influence path went down with it, silently at every level:
#
#   .dfbeta_glm()            errors
#   .consistency_glm_pieces() catches it and returns NULL
#   consistency_resample()   returns all-NA rates
#   the consistency engine   falls back to literal splitting, ignoring
#                            consistency_method = "resample" without a message
#   fs_mr_inference()        can fit NO candidate, so it returns the short
#                            variant with n_family = 0 and no `debiased` element
#
# WHY THE ORIGINAL VALIDATION MISSED IT.  The validation grid in
# dev/identifier-alignment/dfbeta_glm_validation.qmd had a "gaussian identity"
# control, and it passed.  But that control was a glm(family = gaussian()),
# where fit$weights is a vector of ones.  It exercised the gaussian FAMILY and
# never the lm PATH -- the two coincide numerically and differ structurally,
# and only the structural difference mattered.  So the grid's most-exact row
# was the one least able to catch this.
#
# The lesson generalises past this function: a control that passes because the
# quantity is easy is not evidence that the code path is reachable.  Hence the
# end-to-end assertions at the bottom -- they run forestsearch() and inspect
# the MR object, rather than checking the leaf and inferring the rest.

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

.dfb_data <- local({
  set.seed(20260805L)
  n <- 300L
  d <- data.frame(Treat = rbinom(n, 1L, 0.5), x = rnorm(n))
  d$Y  <- 1 + 0.5 * d$Treat + 0.3 * d$x + rnorm(n)
  d$pw <- runif(n, 0.5, 3)
  d
})


# ---------------------------------------------------------------------------
# 1. The leaf: lm() no longer errors, and agrees with stats::dfbeta()
# ---------------------------------------------------------------------------

test_that(".dfbeta_glm() handles an unweighted lm(), which has NULL weights", {
  fit <- stats::lm(Y ~ Treat + x, data = .dfb_data)

  # The precondition that made this fail.  Pinned explicitly: if a future R
  # gives lm() a weights vector, this test's premise is gone and the reader
  # should know that rather than see an unexplained pass.
  expect_null(fit$weights)

  expect_no_error(db <- forestsearch:::.dfbeta_glm(fit))
  expect_true(is.matrix(db))
  expect_identical(nrow(db), nrow(.dfb_data))
  expect_true(all(is.finite(db)))
})

test_that(".dfbeta_glm() matches stats::dfbeta() wherever the two are comparable", {
  # Identity link with gaussian errors: the working residual IS the response
  # residual, so the R 4.6.0 change to weighted.residuals() does not bite and
  # the two computations must agree on EVERY R version.  These are the only
  # cases where agreement is a version-independent assertion -- for
  # glm(binomial) the two differ on R < 4.6.0 BY DESIGN, which is the whole
  # reason .dfbeta_glm() exists, so no agreement is asserted there.
  cases <- list(
    "lm(), no prior weights" =
      stats::lm(Y ~ Treat + x, data = .dfb_data),
    "lm(), prior weights" =
      stats::lm(Y ~ Treat + x, data = .dfb_data, weights = .dfb_data$pw),
    "glm(gaussian) -- the control that was blind to the lm path" =
      stats::glm(Y ~ Treat + x, data = .dfb_data, family = stats::gaussian())
  )

  for (nm in names(cases)) {
    fit <- cases[[nm]]
    expect_equal(forestsearch:::.dfbeta_glm(fit), stats::dfbeta(fit),
                 tolerance = 1e-8, ignore_attr = TRUE, label = nm)
  }
})


# ---------------------------------------------------------------------------
# 2. One level up: the MD pieces are produced, not swallowed
# ---------------------------------------------------------------------------

test_that(".consistency_glm_pieces() returns usable pieces for effect_measure = 'MD'", {
  pc <- forestsearch:::.consistency_glm_pieces(
    .dfb_data, outcome_type = "continuous", effect_measure = "MD",
    treat.name = "Treat", outcome.name = "Y")

  # Before the fix this was NULL -- the failure that propagated everywhere.
  expect_false(is.null(pc))

  expect_length(pc$dfbeta, nrow(.dfb_data))
  expect_true(is.finite(pc$beta_hat))
  expect_true(is.finite(pc$sigma_D) && pc$sigma_D > 0)
  expect_false(pc$log_scale)          # MD is an identity-scale measure
  expect_identical(pc$measure, "MD")

  # sigma_D is sqrt(sum(dfbeta^2)) and nothing else (manuscript Sec 2.2).
  expect_equal(pc$sigma_D, sqrt(sum(pc$dfbeta^2)), tolerance = 1e-12)

  # Sec 2.2 also requires the influences to sum to zero.
  expect_equal(sum(pc$dfbeta), 0, tolerance = 1e-8)

  # beta_hat is the treatment coefficient of the model actually fitted.
  expect_equal(pc$beta_hat,
               unname(stats::coef(stats::lm(Y ~ Treat, data = .dfb_data))[["Treat"]]),
               tolerance = 1e-10)
})

test_that("consistency_resample() produces a rate for MD instead of NA", {
  rr <- consistency_resample(
    .dfb_data, outcome_type = "continuous", effect_measure = "MD",
    comparison_threshold = 0, treat.name = "Treat", outcome.name = "Y",
    method = "closed")

  expect_false(is.na(rr$rate_closed))   # was NA_real_ before the fix
  expect_true(rr$rate_closed >= 0 && rr$rate_closed <= 1)
  expect_true(is.finite(rr$sigma_D) && rr$sigma_D > 0)
})


# ---------------------------------------------------------------------------
# 3. End to end: forestsearch() on a continuous outcome reaches a real MR object
# ---------------------------------------------------------------------------

test_that("MR runs end-to-end on a continuous outcome and de-biases a real family", {
  # This is the assertion the leaf checks could not make.  With the pre-fix
  # body every candidate is dropped in .fs_mr_assemble(), so fs_mr_inference()
  # returns the SHORT variant: n_family = 0, a `note`, and no `debiased`
  # element.  Verified against the pre-fix body -- it fails all three.
  dc <- .make_continuous_data(N = 400L, MD_harm = 2.0)

  args <- .fs_args_for(
    "continuous",
    confounders = c("age", "biomarker_hi", "biomarker", "sex"),
    extra = list(
      use_grf = FALSE, use_lasso = FALSE, quiet = TRUE,
      consistency_method = "resample",
      # MD is a DIFFERENCE, so the shared hr.* thresholds are on the mean
      # scale here, not the ratio scale: 1.25 would be a huge mean shift.
      hr.threshold   = 0.5,
      hr.consistency = 0.0,
      mr_inference      = TRUE,
      mr_inference_args = list(draws = 200L, seed = 11L)))

  fs <- .run_fs_capture(dc, args)$value
  skip_if(is.null(fs$sg.harm), "no subgroup identified on this fixture")

  mr <- fs$mr_inference
  expect_false(is.null(mr))

  # The three marks of the short variant, each absent.
  expect_null(mr$note)
  expect_gt(mr$n_family, 0L)
  expect_true("debiased" %in% names(mr))

  expect_false(is.na(mr$selected_index))
  expect_true(is.finite(mr$naive$est))
  expect_true(is.finite(mr$debiased$est))
  expect_true(is.finite(mr$debiased$se) && mr$debiased$se > 0)

  # The de-biasing actually engaged: a family was competed over on some draws.
  expect_gt(mr$selection_rate, 0)
  expect_true(is.finite(mr$selection_bias))

  # Identity scale is carried through -- MD must not be exponentiated.
  expect_false(mr$log_scale)
})
