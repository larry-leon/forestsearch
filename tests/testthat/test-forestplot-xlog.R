# ============================================================================
# test-forestplot-xlog.R
#
# Regression guard for GitHub #7.  plot_subgroup_results_forestplot() passed
# `xlog = xlog` to forestploter::forest(), which has never had an `xlog`
# argument, so the HR/OR axis was linear whatever `xlog` said.  Since
# 0.4.0.9000 the call uses the forestploter >= 1.2.0 API: forest(style = ),
# set_xaxis(x_trans = ), set_labs(), and refline_gp in the theme.
#
# Fixture: the null-DGM DINA fit of test-no-subgroup-dina-grf-plot.R (no FS
# subgroup), with one reference subgroup in subgroup_list so the call
# reaches forestploter::forest().
#
# The plot object exposes forestploter's recipe as
# attr(p, "forest_recipe"); its $xaxis$x_trans records the axis scale.
# ============================================================================


.xlog_fixture <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      df <- .make_survival_data(N = 160L, HR_harm = 1.0, seed = 7L)
      args <- .fs_args_for("survival", confounders = c("age", "stage", "sex"),
                           extra = list(hr.threshold = 10.0,
                                        pconsistency.threshold = 0.99,
                                        subgroup_method = "dina"))
      fs <- suppressWarnings(do.call(forestsearch, c(list(df.analysis = df), args)))
      cache <<- fs
    }
    cache
  }
})

.xlog_plot <- function(fs, ...) {
  plot_subgroup_results_forestplot(
    fs_results    = list(fs.est = fs),
    df_analysis   = as.data.frame(fs$args_call_all$df.analysis),
    subgroup_list = list(sex1 = list(subset_expr = "sex == 1",
                                     name = "sex = 1", type = "reference")),
    outcome.name  = "time",
    event.name    = "event",
    treat.name    = "treat",
    ...
  )
}


# NOTE: forestploter's "Unknown arguments ignored" / "will be deprecated"
# messages fire once per R session, so this assertion is only meaningful as
# the FIRST forestploter plot drawn in the session.  Keep it the first test
# of this file and run the file in its own session.
test_that("forestplot call emits no forestploter deprecation / unknown-argument messages", {
  skip_on_cran()
  skip_if_not_installed("forestploter", "1.2.0")
  fs <- .xlog_fixture()

  msgs <- character()
  res <- withCallingHandlers(
    .xlog_plot(fs),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m))
      invokeRestart("muffleMessage")
    }
  )
  expect_s3_class(res$plot, "forestplot")
  bad <- grep("Unknown arguments ignored|will be deprecated", msgs, value = TRUE)
  expect_length(bad, 0L)
})


test_that("xlog = TRUE and xlog = FALSE both draw, and differ", {
  skip_on_cran()
  skip_if_not_installed("forestploter", "1.2.0")
  fs <- .xlog_fixture()

  p_log <- suppressMessages(.xlog_plot(fs, xlog = TRUE))$plot
  p_lin <- suppressMessages(.xlog_plot(fs, xlog = FALSE))$plot
  expect_s3_class(p_log, "forestplot")
  expect_s3_class(p_lin, "forestplot")
  expect_false(identical(p_log, p_lin))
})


test_that("xlog maps to set_xaxis(x_trans = ) in the forestploter recipe", {
  skip_on_cran()
  skip_if_not_installed("forestploter", "1.2.0")
  fs <- .xlog_fixture()

  p_log <- suppressMessages(.xlog_plot(fs, xlog = TRUE))$plot
  p_lin <- suppressMessages(.xlog_plot(fs, xlog = FALSE))$plot
  expect_identical(attr(p_log, "forest_recipe")$xaxis$x_trans, "log")
  expect_identical(attr(p_lin, "forest_recipe")$xaxis$x_trans, "none")
})
