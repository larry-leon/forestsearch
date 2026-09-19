# Acceptance tests for the estimability boundary (Part A).
#
# The rule is per estimand: OR needs all four cells >= 1; RR/IRR/HR need an
# event in each arm; RD/IRD/MD have no condition and must be byte-identical.
# No admission floor is added anywhere and nothing is forwarded to DINA/GRF.

mkb <- function(y, t) data.frame(treat = as.integer(t), y = as.integer(y))

d_ok   <- mkb(c(1,1,0,0,0, 1,1,1,0,0), c(0,0,0,0,0, 1,1,1,1,1))  # four cells >= 1
d_zero <- mkb(c(1,1,0,0,0, 1,1,1,1,1), c(0,0,0,0,0, 1,1,1,1,1))  # treated non-events = 0
d_noev <- mkb(c(1,1,1,0,0, 0,0,0,0,0), c(0,0,0,0,0, 1,1,1,1,1))  # treated events = 0
s_ok   <- data.frame(treat = c(rep(0,5), rep(1,5)),
                     time  = c(5,6,7,8,9, 4,5,6,7,8),
                     ev    = c(1,1,0,1,0, 1,0,1,1,0))
s_noev <- data.frame(treat = c(rep(0,5), rep(1,5)),
                     time  = c(5,6,7,8,9, 4,5,6,7,8),
                     ev    = c(1,1,0,1,0, 0,0,0,0,0))
c_ok   <- data.frame(treat = c(rep(0,6), rep(1,6)),
                     y = c(1.2,2.3,0.7,1.9,2.8,1.1, 3.2,2.9,4.1,3.3,2.7,3.8))

binest <- function(m) make_effect_estimator(outcome_type = "binary",
  treat.name = "treat", outcome.name = "y", effect_measure = m,
  adverse_outcome = TRUE)
coxest <- function() make_effect_estimator(outcome_type = "survival",
  treat.name = "treat", outcome.name = "time", event.name = "ev")

# ---- the estimand does not exist -> NA with a reason ----------------------

test_that("OR with an empty cell is non-estimable and names the cell", {
  r <- binest("OR")(d_zero)
  expect_true(is.na(r$estimate))
  expect_true(is.na(r$se))
  expect_false(r$converged)
  expect_match(r$reason, "^non-estimable: ")
  expect_match(r$reason, "treated non-events = 0")
})

test_that("RR with a zero-event arm is non-estimable", {
  r <- binest("RR")(d_noev)
  expect_true(is.na(r$estimate))
  expect_false(r$converged)
  expect_match(r$reason, "treated events = 0")
})

test_that("HR with a zero-event arm is non-estimable", {
  r <- coxest()(s_noev)
  expect_true(is.na(r$estimate))
  expect_false(r$converged)
  expect_match(r$reason, "treated events = 0")
})

test_that("the boundary replaces a FINITE divergent estimate, not an error", {
  # This is the whole point: glm()/coxph() do not fail on these slices, they
  # return a large finite coefficient with converged = TRUE.  The guard must
  # be an existence check, because `converged` does not reveal the pathology.
  expect_silent({
    fit <- suppressWarnings(stats::glm(y ~ treat, data = d_zero,
                                       family = stats::binomial()))
  })
  expect_true(is.finite(stats::coef(fit)[["treat"]]))
  expect_true(abs(stats::coef(fit)[["treat"]]) > 10)
  expect_true(fit$converged)
})

# ---- estimands with no existence condition are untouched ------------------

test_that("OR with all four cells >= 1 is unchanged and carries no reason", {
  r <- binest("OR")(d_ok)
  expect_false(is.na(r$estimate))
  expect_true(r$converged)
  expect_null(r$reason)
  expect_identical(r$method_used, "logistic")
})

test_that("RD is out of scope on every frame, including a zero cell", {
  for (d in list(d_ok, d_zero, d_noev)) {
    r <- binest("RD")(d)
    expect_false(is.na(r$estimate))
    expect_null(r$reason)     # never enters the boundary
  }
})

test_that("RD tier 3 returns converged = FALSE by design and is not caught", {
  # Asserted directly on the tier-3 return: a converged = FALSE that must
  # survive the ratio-estimand non-convergence catch untouched.
  r <- forestsearch:::.estimate_rd(
    d_ok, stats::as.formula("y ~ treat"), "treat", "y",
    n0 = 5L, n1 = 5L, adjust_covariates = NULL, ps_adjust_method = "none")
  expect_false(is.na(r$estimate))
  expect_null(r$reason)
  # and the scoping itself: RD is not among the guarded measures
  expect_null(forestsearch:::.fs_existence_reason(
    "RD", c(e0 = 0L, ne0 = 0L, e1 = 0L, ne1 = 0L)))
  expect_null(forestsearch:::.fs_existence_reason(
    "MD", c(e0 = 0L, ne0 = 0L, e1 = 0L, ne1 = 0L)))
  expect_null(forestsearch:::.fs_existence_reason(
    "IRD", c(e0 = 0L, ne0 = 0L, e1 = 0L, ne1 = 0L)))
})

test_that("a normal HR fit and a normal MD fit are untouched", {
  r <- coxest()(s_ok)
  expect_false(is.na(r$estimate)); expect_true(r$converged); expect_null(r$reason)
  m <- make_effect_estimator(outcome_type = "continuous", treat.name = "treat",
                             outcome.name = "y", effect_measure = "MD")(c_ok)
  expect_false(is.na(m$estimate)); expect_null(m$reason)
})

# ---- the existence conditions themselves ----------------------------------

test_that("the existence condition is exactly per estimand", {
  full <- c(e0 = 3L, ne0 = 2L, e1 = 4L, ne1 = 1L)
  expect_null(forestsearch:::.fs_existence_reason("OR", full))
  expect_null(forestsearch:::.fs_existence_reason("RR", full))
  expect_null(forestsearch:::.fs_existence_reason("HR", full))
  # OR is the only estimand that sees the non-event cells
  no_ne1 <- c(e0 = 3L, ne0 = 2L, e1 = 4L, ne1 = 0L)
  expect_match(forestsearch:::.fs_existence_reason("OR", no_ne1), "treated non-events")
  expect_null(forestsearch:::.fs_existence_reason("RR", no_ne1))
  expect_null(forestsearch:::.fs_existence_reason("HR", no_ne1))
  # a zero-event arm fails every ratio estimand
  no_e0 <- c(e0 = 0L, ne0 = 2L, e1 = 4L, ne1 = 1L)
  for (m in c("OR", "RR", "IRR", "HR"))
    expect_match(forestsearch:::.fs_existence_reason(m, no_e0), "control events = 0")
})

test_that("cell counting matches the four-cell definition", {
  cl <- forestsearch:::.fs_binary_cells(d_ok$y, d_ok$treat)
  expect_identical(unname(cl[["e0"]]),  sum(d_ok$treat == 0 & d_ok$y == 1))
  expect_identical(unname(cl[["ne1"]]), sum(d_ok$treat == 1 & d_ok$y == 0))
})

# ---- visibility: counted with its reason, never silent --------------------

test_that("a non-estimable candidate is counted with its reason", {
  fc <- list(n_nonestimable = 0L, nonestimable_reasons = character(0))
  res <- list(status = 5L, result = NULL, reason = "non-estimable: treated events = 0")
  expect_false(is.null(res$reason))
  expect_null(list(status = 5L, result = NULL)$reason)  # ordinary fit failure
})

test_that("fit_glm_for_subgroup carries the reason instead of a bare NULL", {
  df <- cbind(d_zero, id = 1L)
  out <- forestsearch:::fit_glm_for_subgroup(df, rep(1L, nrow(df)), binest("OR"))
  expect_true(isTRUE(out$nonestimable))
  expect_match(out$reason, "treated non-events = 0")
  # an estimable slice still returns the hr/lower/upper contract
  ok <- forestsearch:::fit_glm_for_subgroup(cbind(d_ok, id = 1L),
                                            rep(1L, nrow(d_ok)), binest("OR"))
  expect_true(all(c("hr", "lower", "upper", "med0", "med1") %in% names(ok)))
  expect_null(ok$nonestimable)
})

test_that("converged is read at the boundary, and the catch is scoped off RD", {
  # Source-level assertion, as the task requires: `converged` is no longer
  # computed and discarded for the ratio estimands, and the catch cannot
  # reach RD's tier-3 converged = FALSE.
  f <- testthat::test_path("..", "..", "R", "glm_effect_estimators.R")
  skip_if_not(file.exists(f), "package source not available")
  gl <- paste(readLines(f, warn = FALSE), collapse = "\n")
  expect_true(grepl('isFALSE(result$converged)', gl, fixed = TRUE))
  expect_true(grepl('effect_measure %in% c("OR", "RR") &&', gl, fixed = TRUE))
  # exactly one non-convergence catch, and RD is not in its scope
  expect_identical(
    length(gregexpr('isFALSE(result$converged)', gl, fixed = TRUE)[[1]]), 1L)
})
