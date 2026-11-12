# =============================================================================
# Displays use the fit's own orientation (F12, cran-prep part 4).
#
# plot_sg_glm_outcomes(), plot_subgroup_results_forestplot() (GLM rows) and
# mr_estimates_table(boot_results = NULL) rebuilt their effect estimator
# without the fit's adverse_outcome, so it fell back to make_effect_estimator()'s
# formal (TRUE) whatever the fit used.  For a continuous fit at the default
# orientation (adverse_outcome unset -> FALSE: higher y is better) the fit
# reports MD on -y (> 0 = treatment harm) while the displays showed MD on y:
# same magnitude, opposite sign.
#
# Invariant: each display's estimate for a set equals the fit's own estimate
# for that set, sign included.  The arm means (and their raw difference) stay
# raw; they are not checked here.
#
# Expected values are computed by hand from the data (no package code):
#   - unadjusted OLS of y on a 0/1 treatment gives the treatment coefficient
#     mean(y | T = 1) - mean(y | T = 0) and the pooled-variance standard
#     error sqrt(s2 * (1/n0 + 1/n1)), s2 = (SS0 + SS1) / (n0 + n1 - 2);
#   - orientation (?forestsearch, @param adverse_outcome): TRUE means higher
#     values indicate harm, FALSE (the continuous default) that higher values
#     indicate benefit; the package's convention is effect > 0 = treatment
#     harm, so under FALSE the fit's MD is the coefficient for -y, i.e.
#     -(mean(y | T = 1) - mean(y | T = 0)), with the same standard error.
#   The hand value for the selected subgroup is also checked against the fit's
#   own reported estimate (max_sg_est), which does not depend on the displays.
# Tolerances are the display rounding: half a unit of the last printed digit
# (labels at 3 decimals: 5e-4; the MR table at digits = 4: 5e-5 for the
# estimate, 1e-4 for its limits, whose z is the literal 1.96 rather than
# qnorm(0.975), a difference below 1e-5 at these standard errors).
# =============================================================================

.dor_pa <- list(plan = "sequential", workers = 1L, show_message = FALSE)

# Mirror twins: `harm` has a planted region where treatment RAISES y by 1.5
# (age <= 50 & biomarker_hi); `ben` is the same data with y negated, so under
# the continuous default (higher y = better) the same region is where
# treatment harms.
.dor_data <- function() {
  df <- .make_continuous_data(N = 300L, MD_harm = 1.5, seed = 42L)
  df$biomarker_hi <- as.numeric(as.character(df$biomarker_hi))
  df$sex <- as.numeric(as.character(df$sex))
  ben <- df; ben$y <- -ben$y
  list(harm = df, ben = ben)
}

.dor_fit <- function(df, ...) {
  suppressWarnings(suppressMessages(forestsearch(
    df.analysis = df,
    confounders.name = c("age", "biomarker", "biomarker_hi", "sex"),
    outcome.name = "y", treat.name = "treat", id.name = "id",
    outcome_type = "continuous", effect_measure = "MD",
    effect.threshold = 0.5, consistency.threshold = 0,
    pconsistency.threshold = 0.5, fs.splits = 20, n.min = 30, maxk = 2,
    use_lasso = FALSE, use_grf = FALSE, sg_focus = "hr", is.RCT = TRUE,
    mr_inference = TRUE,
    mr_inference_args = list(draws = 200L, seed = 11L, ci_method = "ij",
                             include_complement = FALSE),
    details = FALSE, quiet = TRUE, seedit = 1L, parallel_args = .dor_pa, ...)))
}

# One fit per orientation, built once for the file.
.dor_fixture <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      d <- .dor_data()
      cache <<- list(
        data  = d,
        fit_f = .dor_fit(d$ben),                           # unset -> FALSE
        fit_t = .dor_fit(d$harm, adverse_outcome = TRUE))
    }
    cache
  }
})

# Hand MD (see header): OLS coefficient and pooled SE, oriented.
.dor_md <- function(d, adverse) {
  y0 <- d$y[d$treat == 0]; y1 <- d$y[d$treat == 1]
  n0 <- length(y0); n1 <- length(y1)
  s2 <- (sum((y0 - mean(y0))^2) + sum((y1 - mean(y1))^2)) / (n0 + n1 - 2)
  raw <- mean(y1) - mean(y0)
  list(est = if (adverse) raw else -raw, se = sqrt(s2 * (1 / n0 + 1 / n1)))
}

# Sets of the fit: H (treat.recommend == 0), Hc (== 1), ITT (all of df.est).
.dor_sets <- function(fit) {
  d <- fit$df.est
  list(ITT = d, H = d[d$treat.recommend == 0, ], Hc = d[d$treat.recommend == 1, ])
}

# Numbers printed in a label such as "MD = -1.094 (-1.665, -0.523)".
.dor_nums <- function(s) {
  as.numeric(regmatches(s, gregexpr("-?[0-9]+\\.[0-9]+", s))[[1]])
}

.dor_check_glm_plot <- function(fit, adverse) {
  ref <- lapply(.dor_sets(fit), .dor_md, adverse = adverse)
  z <- stats::qnorm(0.975)                     # conf.level = 0.95 (default)
  p <- plot_sg_glm_outcomes(fit, outcome.name = "y", treat.name = "treat")
  eff <- p$effects
  lab <- c(ITT = "ITT", H = "Questionable", Hc = "Recommend")   # default names
  for (k in names(lab)) {
    got <- .dor_nums(eff$effect_label[as.character(eff$Subgroup) == lab[[k]]])
    want <- c(ref[[k]]$est, ref[[k]]$est - z * ref[[k]]$se,
              ref[[k]]$est + z * ref[[k]]$se)
    expect_length(got, 3L)
    expect_lt(max(abs(got - want)), 5e-4 + 1e-9, label = paste("label", k))
    expect_identical(sign(got[1]), sign(want[1]), label = paste("sign", k))
  }
}

.dor_check_forestplot <- function(fit, df, adverse) {
  ref <- lapply(.dor_sets(fit), .dor_md, adverse = adverse)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  fp <- suppressWarnings(suppressMessages(plot_subgroup_results_forestplot(
    fs_results = list(fs.est = fit), df_analysis = df,
    outcome.name = "y", event.name = NULL, treat.name = "treat")))
  dt <- fp$data; rt <- fp$row_types
  expect_identical(length(rt), nrow(dt))
  # Row types (est.scale = "hr"): "itt"; "posthoc" = the complement's
  # full-analysis row; "posthoc_complement" = the selected subgroup's.
  pick <- c(ITT = "itt", Hc = "posthoc", H = "posthoc_complement")
  for (k in names(pick)) {
    i <- which(rt == pick[[k]])
    expect_length(i, 1L)
    expect_equal(dt$est[i], ref[[k]]$est, tolerance = 1e-10,
                 label = paste("forest est", k))
    expect_equal(dt$se[i], ref[[k]]$se, tolerance = 1e-10,
                 label = paste("forest se", k))
  }
}

.dor_check_mr_table <- function(fit, adverse) {
  ref <- lapply(.dor_sets(fit), .dor_md, adverse = adverse)
  tbl <- mr_estimates_table(fit, digits = 4)
  d <- as.data.frame(tbl[["_data"]], stringsAsFactors = FALSE)
  expect_identical(d$Subgroup, c("Qstnbl", "Recmnd"))
  rows <- c(H = "Qstnbl", Hc = "Recmnd")
  for (k in names(rows)) {
    got <- .dor_nums(d[["MD (95% CI)"]][d$Subgroup == rows[[k]]])
    want <- c(ref[[k]]$est, ref[[k]]$est - 1.96 * ref[[k]]$se,
              ref[[k]]$est + 1.96 * ref[[k]]$se)
    expect_length(got, 3L)
    expect_lt(abs(got[1] - want[1]), 5e-5 + 1e-9, label = paste("MR est", k))
    expect_identical(sign(got[1]), sign(want[1]), label = paste("MR sign", k))
    # The selected row's limits are fs$mr_inference$naive's robust (sandwich)
    # limits (?mr_estimates_table), not the pooled-SE limits computed here;
    # the complement row (no include_complement) is the table's own OLS fit.
    if (k == "Hc") {
      expect_lt(max(abs(got[2:3] - want[2:3])), 1e-4, label = paste("MR CI", k))
    }
  }
}


test_that("continuous fit at the default orientation: the fixture is what the tests assume", {
  skip_on_cran()
  fx <- .dor_fixture()
  fit <- fx$fit_f
  expect_identical(fit$args_call_all$adverse_outcome, FALSE)  # forestsearch()'s rule
  expect_false(is.null(fit$sg.harm))
  expect_false(is.null(fit$mr_inference))
  expect_null(fit$mr_inference$complement)       # include_complement = FALSE
  # The fit's own estimate for its subgroup is the hand MD on -y, and the
  # planted harm reads as harm (> 0).
  h <- .dor_md(.dor_sets(fit)$H, adverse = FALSE)
  expect_equal(fit$max_sg_est, h$est, tolerance = 1e-10)
  expect_gt(h$est, 0)
  expect_equal(fit$mr_inference$naive$est, h$est, tolerance = 1e-10)
})

test_that("plot_sg_glm_outcomes(): effect annotations are the fit's own estimates, sign included", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  fx <- .dor_fixture()
  .dor_check_glm_plot(fx$fit_f, adverse = FALSE)
})

test_that("plot_subgroup_results_forestplot(): GLM rows are the fit's own estimates, sign included", {
  skip_on_cran()
  skip_if_not_installed("forestploter", "1.2.0")
  fx <- .dor_fixture()
  .dor_check_forestplot(fx$fit_f, fx$data$ben, adverse = FALSE)
})

test_that("mr_estimates_table(boot_results = NULL): both rows are the fit's own estimates, sign included", {
  skip_on_cran()
  skip_if_not_installed("gt")
  fx <- .dor_fixture()
  # The complement row is built by the table's own estimator (no
  # include_complement); the selected row is fs$mr_inference$naive.
  .dor_check_mr_table(fx$fit_f, adverse = FALSE)
})

test_that("no stored adverse_outcome: the displays use forestsearch()'s rule (continuous -> FALSE)", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("forestploter", "1.2.0")
  skip_if_not_installed("gt")
  fx <- .dor_fixture()
  fit <- fx$fit_f
  fit$args_call_all$adverse_outcome <- NULL      # an object not built by forestsearch()
  expect_false("adverse_outcome" %in% names(fit$args_call_all))
  .dor_check_glm_plot(fit, adverse = FALSE)
  .dor_check_forestplot(fit, fx$data$ben, adverse = FALSE)
  .dor_check_mr_table(fit, adverse = FALSE)
})

test_that("adverse_outcome = TRUE: the displays show the raw MD, which is the fit's own estimate", {
  skip_on_cran()
  skip_if_not_installed("ggplot2")
  skip_if_not_installed("forestploter", "1.2.0")
  skip_if_not_installed("gt")
  fx <- .dor_fixture()
  fit <- fx$fit_t
  expect_identical(fit$args_call_all$adverse_outcome, TRUE)
  expect_false(is.null(fit$sg.harm))
  h <- .dor_md(.dor_sets(fit)$H, adverse = TRUE)
  expect_equal(fit$max_sg_est, h$est, tolerance = 1e-10)
  .dor_check_glm_plot(fit, adverse = TRUE)
  .dor_check_forestplot(fit, fx$data$harm, adverse = TRUE)
  .dor_check_mr_table(fit, adverse = TRUE)
  # Mirror twins: the FALSE fit on -y selects the same subgroup, and its
  # displays carry the same oriented estimates as this fit's.
  expect_identical(fx$fit_f$sg.harm, fit$sg.harm)
})
