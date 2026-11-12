# =============================================================================
# GLM bootstrap uses the parent fit's estimator settings (A-009).
#
# forestsearch_bootstrap_dofuture() rebuilt its GLM estimator without the
# parent's resolved adverse_outcome and ps_adjust_method, so H_obs, H_star,
# Hstar_obs, Hstar_star (and FSsg_tab) used adverse_outcome = TRUE and
# ps_adjust_method = "none" whatever the parent used.  The parent's own rebuild
# under active propensity adjustment also dropped adverse_outcome.
#
# Invariant: the bootstrap's estimate of the parent's subgroup on the original
# data (H_estimates$H0 = H_obs on the identity scale) equals the parent fit's
# own estimate for that subgroup (max_sg_est), sign included.
# =============================================================================

.a009_pa <- list(plan = "sequential", workers = 1L, show_message = FALSE)

.a009_data <- function() {
  df <- .make_continuous_data(N = 300L, MD_harm = 1.5, seed = 42L)
  df$biomarker_hi <- as.numeric(as.character(df$biomarker_hi))  # numeric
  df$sex <- as.numeric(as.character(df$sex))                    # PS confounders
  df
}

.a009_fit <- function(df, ...) {
  suppressWarnings(suppressMessages(forestsearch(
    df.analysis = df,
    confounders.name = c("age", "biomarker", "biomarker_hi", "sex"),
    outcome.name = "y", treat.name = "treat", id.name = "id",
    outcome_type = "continuous", effect_measure = "MD",
    effect.threshold = 0.5, consistency.threshold = 0,
    pconsistency.threshold = 0.5, fs.splits = 20, n.min = 30, maxk = 2,
    use_lasso = FALSE, use_grf = FALSE, sg_focus = "hr",
    details = FALSE, quiet = TRUE, seedit = 1L, parallel_args = .a009_pa, ...)))
}

.a009_boot <- function(fit) {
  suppressWarnings(suppressMessages(forestsearch_bootstrap_dofuture(
    fit, nb_boots = 3L, details = FALSE, show_three = FALSE,
    parallel_args = .a009_pa)))
}

.a009_bc_cols <- c("H_biasadj_1", "H_biasadj_2", "Hc_biasadj_1", "Hc_biasadj_2")


test_that("bootstrap H0 equals the parent's subgroup estimate (adverse_outcome FALSE and TRUE)", {
  skip_on_cran()
  df <- .a009_data()
  ben <- df; ben$y <- -ben$y          # higher y = better: the mirror twin

  fit_f <- .a009_fit(ben, is.RCT = TRUE)                       # unset -> FALSE
  fit_t <- .a009_fit(df,  is.RCT = TRUE, adverse_outcome = TRUE)
  expect_identical(fit_f$args_call_all$adverse_outcome, FALSE)
  expect_identical(fit_t$args_call_all$adverse_outcome, TRUE)
  expect_false(is.null(fit_f$sg.harm))
  expect_identical(fit_f$sg.harm, fit_t$sg.harm)

  bt_f <- .a009_boot(fit_f)
  bt_t <- .a009_boot(fit_t)
  expect_false(is.null(bt_f$H_estimates))
  expect_false(is.null(bt_t$H_estimates))

  # Parent's own estimate for its subgroup, sign included.
  expect_equal(bt_f$H_estimates$H0, fit_f$max_sg_est, tolerance = 1e-10)
  expect_equal(bt_t$H_estimates$H0, fit_t$max_sg_est, tolerance = 1e-10)
  expect_equal(bt_f$H_estimates$H0,
               fit_f$grp.consistency$out_sg$result$hr[1], tolerance = 1e-10)
  expect_gt(bt_f$H_estimates$H0, 0)   # planted harm reads as harm on both

  # The mirror twins give the same bootstrap output: the FALSE fit's
  # bootstrap is on the parent's (harm-oriented) scale, not the raw one.
  expect_equal(as.data.frame(bt_f$H_estimates), as.data.frame(bt_t$H_estimates),
               tolerance = 1e-10)
  expect_equal(as.data.frame(bt_f$Hc_estimates), as.data.frame(bt_t$Hc_estimates),
               tolerance = 1e-10)
  expect_equal(as.data.frame(bt_f$results)[, .a009_bc_cols],
               as.data.frame(bt_t$results)[, .a009_bc_cols], tolerance = 1e-10)
  expect_identical(bt_f$FSsg_tab[["MD (95% CI)"]], bt_t$FSsg_tab[["MD (95% CI)"]])
})


test_that("bootstrap estimator carries the parent's ps_adjust_method (iptw)", {
  skip_on_cran()
  df <- .a009_data()
  fit <- .a009_fit(df, is.RCT = FALSE, ps_method = "logistic",
                   ps_adjust_method = "iptw", adverse_outcome = TRUE)
  expect_identical(fit$args_call_all$ps_adjust_method, "iptw")
  expect_false(is.null(fit$sg.harm))
  expect_true("sw" %in% names(fit$df.est))

  H <- subset(fit$df.est, treat.recommend == 0)
  unweighted <- stats::coef(stats::lm(y ~ treat, data = H))[["treat"]]
  weighted   <- stats::coef(stats::lm(y ~ treat, data = H, weights = sw))[["treat"]]
  expect_equal(fit$max_sg_est, weighted, tolerance = 1e-10)
  expect_gt(abs(weighted - unweighted), 1e-6)   # the two are distinguishable

  bt <- .a009_boot(fit)
  expect_false(is.null(bt$H_estimates))
  expect_equal(bt$H_estimates$H0, fit$max_sg_est, tolerance = 1e-10)
})


test_that("parent's propensity-adjusted rebuild keeps adverse_outcome", {
  skip_on_cran()
  df <- .a009_data()
  ben <- df; ben$y <- -ben$y
  fit_t <- .a009_fit(df,  is.RCT = FALSE, ps_method = "logistic",
                     ps_adjust_method = "iptw", adverse_outcome = TRUE)
  fit_f <- .a009_fit(ben, is.RCT = FALSE, ps_method = "logistic",
                     ps_adjust_method = "iptw", adverse_outcome = FALSE)
  expect_false(is.null(fit_t$sg.harm))
  expect_identical(fit_f$sg.harm, fit_t$sg.harm)
  expect_equal(fit_f$max_sg_est, fit_t$max_sg_est, tolerance = 1e-10)
})
