# =============================================================================
# Continuous IPTW weighted fit (regression, 0.4.0.9000): the weights did not
# resolve inside .make_lm_estimator()'s closure, so every candidate failed.
# Closure-level identity against lm(weights = ...), and the path end to end.
# =============================================================================

test_that("continuous IPTW closure reproduces the weighted OLS fit", {
  set.seed(11)
  n <- 200L
  d <- data.frame(treat = rep(0:1, each = n / 2))
  d$y  <- 1 + 0.5 * d$treat + stats::rnorm(n)
  ps   <- stats::runif(n, 0.2, 0.8)
  d$sw <- ifelse(d$treat == 1L, mean(d$treat) / ps, (1 - mean(d$treat)) / (1 - ps))

  est <- make_effect_estimator(outcome_type = "continuous", treat.name = "treat",
                               outcome.name = "y", ps_adjust_method = "iptw")
  r   <- est(d)
  ref <- stats::lm(y ~ treat, data = d, weights = d$sw)   # d is visible here
  expect_true(isTRUE(r$converged))
  expect_identical(r$method_used, "ols")
  expect_identical(r$estimate, stats::coef(ref)[["treat"]])
  expect_identical(r$se, sqrt(diag(stats::vcov(ref)))[["treat"]])

  # Without an sw column the same closure is the unweighted fit (unchanged).
  d0 <- d; d0$sw <- NULL
  r0 <- est(d0)
  ref0 <- stats::lm(y ~ treat, data = d0)
  expect_identical(r0$estimate, stats::coef(ref0)[["treat"]])
  expect_identical(r0$se, sqrt(diag(stats::vcov(ref0)))[["treat"]])
})

test_that("continuous IPTW path is reachable end to end and yields fitted candidates", {
  skip_on_cran()   # the closure-level test above is the regression guard everywhere
  df <- .make_continuous_data(N = 300L, MD_harm = 1.5, seed = 42L)
  df$biomarker_hi <- as.numeric(as.character(df$biomarker_hi))  # numeric confounders
  df$sex <- as.numeric(as.character(df$sex))                    # for the PS model
  fit <- suppressWarnings(suppressMessages(forestsearch(
    df.analysis = df, confounders.name = c("age", "biomarker", "biomarker_hi", "sex"),
    outcome.name = "y", event.name = "y", treat.name = "treat", id.name = "id",
    outcome_type = "continuous", effect_measure = "MD", adverse_outcome = TRUE,
    effect.threshold = 0.5, consistency.threshold = 0.0, pconsistency.threshold = 0.5,
    fs.splits = 20, n.min = 30, maxk = 2, use_lasso = FALSE, use_grf = FALSE,
    sg_focus = "eff", is.RCT = FALSE, ps_method = "logistic", ps_adjust_method = "iptw",
    details = FALSE, quiet = TRUE, seedit = 1L, parallel_args = list(plan = "sequential"))))
  tab <- fit$find.grps$out.found$hr.subgroups
  expect_true(!is.null(tab) && nrow(tab) > 0L)
})
