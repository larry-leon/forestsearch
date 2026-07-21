# ============================================================================
# test-no-subgroup-dina-grf-plot.R
#
# Regression guard for the no-subgroup DINA/GRF fault on the PLOT entry
# points, surfaced by the completeness audit of the summary-table fix and
# parallel to the CV-path fix (981a6c1).
#
# subgroup_method = "dina"/"grf" return df.est = df (a populated frame
# WITHOUT a treat.recommend column) with sg.harm NULL when no subgroup is
# identified, unlike the consistency path (df.est NULL).  Two exported plot
# entry points gated their no-subgroup handling on is.null(df.est) alone,
# so a no-subgroup DINA/GRF fit slipped through:
#
#   * plot.forestsearch()                -> the intended graceful
#       "No subgroup identified -- nothing to plot" no-op was skipped, so
#       the call hard-errored downstream in plot_sg_results().
#   * plot_subgroup_results_forestplot() -> SECTION 3 entered and ran
#       subset(df_fs, treat.recommend == ...) on the absent column,
#       crashing with the cryptic NSE "object 'treat.recommend' not found".
#
# Both guards now key on sg.harm (the actual no-subgroup contract).  A
# subgroup-found fit has sg.harm non-NULL, so each condition reduces to the
# former is.null(df.est) test -- behaviour byte-identical for found fits.
#
# The strict plotters (plot_sg_results, plot_sg_weighted_km,
# plot_sg_glm_outcomes) already fail cleanly via an explicit
# required-column check and are intentionally left unchanged.
# ============================================================================


.make_null_dg_fit_plot <- function(method) {
  forcing <- if (method == "grf") {
    list(hr.threshold = 10.0, pconsistency.threshold = 0.99,
         subgroup_method = "grf", dmin.grf = 1e6)
  } else {
    list(hr.threshold = 10.0, pconsistency.threshold = 0.99,
         subgroup_method = "dina")
  }
  df <- .make_survival_data(N = 160L, HR_harm = 1.0, seed = 7L)
  args <- .fs_args_for("survival", confounders = c("age", "stage", "sex"),
                       extra = forcing)
  fs <- suppressWarnings(do.call(forestsearch, c(list(df.analysis = df), args)))
  if (!is.null(fs$sg.harm)) {
    skip(sprintf("null DGM unexpectedly identified a subgroup (%s)", method))
  }
  # Unified no-subgroup contract (Option i): df.est = NULL for all identifiers.
  expect_null(fs$df.est)
  fs
}


# ---------------------------------------------------------------------------
# plot.forestsearch(): graceful no-op, not a hard error
# ---------------------------------------------------------------------------

for (method in c("dina", "grf")) {

  test_that(sprintf("plot.forestsearch() is a graceful no-op for no-subgroup %s fit",
                    method), {
    skip_on_cran()
    fs <- .make_null_dg_fit_plot(method)

    for (ptype in c("combined", "km", "forest", "summary")) {
      res <- NULL
      # Emits the graceful message and returns the fit invisibly (never
      # errors).  Assign inside expect_message so the value can be checked.
      expect_message(
        res <- plot(fs, type = ptype),
        "No subgroup identified -- nothing to plot"
      )
      expect_identical(res, fs)
    }
  })
}


# ---------------------------------------------------------------------------
# plot_subgroup_results_forestplot(): skips FS rows, no cryptic NSE crash
# ---------------------------------------------------------------------------

for (method in c("dina", "grf")) {

  test_that(sprintf("plot_subgroup_results_forestplot() skips FS rows for no-subgroup %s fit",
                    method), {
    skip_on_cran()
    skip_if_not_installed("forestploter")
    fs <- .make_null_dg_fit_plot(method)
    df <- as.data.frame(fs$args_call_all$df.analysis)

    # Must not raise the pre-fix "object 'treat.recommend' not found" crash.
    # A clean return OR any error that is NOT the treat.recommend NSE
    # failure both prove the no-subgroup block was skipped.
    err <- tryCatch({
      suppressWarnings(suppressMessages(
        plot_subgroup_results_forestplot(
          fs_results  = list(fs.est = fs),
          df_analysis = df,
          outcome.name = "time",
          event.name   = "event",
          treat.name   = "treat")
      ))
      NULL
    }, error = function(e) conditionMessage(e))

    err_msg <- if (is.null(err)) "" else err
    expect_false(
      isTRUE(grepl("treat.recommend", err_msg, fixed = TRUE)),
      info = sprintf("forestplot still crashed on treat.recommend: %s", err_msg)
    )
  })
}
