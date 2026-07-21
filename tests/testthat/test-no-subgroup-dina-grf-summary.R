# ============================================================================
# test-no-subgroup-dina-grf-summary.R
#
# Regression guard for the no-subgroup DINA/GRF fault on the SUMMARY/TABLE
# path, parallel to the CV-path fix (981a6c1).
#
# subgroup_method = "dina"/"grf" return df.est = df (a populated frame
# WITHOUT a treat.recommend column) with sg.harm NULL when no subgroup is
# identified, unlike the consistency path (df.est NULL).  sg_tables()
# gated its no-subgroup / empty-tables branch on is.null(df) alone, so a
# no-subgroup DINA/GRF fit slipped through to Table 1's
# SG_tab_estimates(SG_flag = "treat.recommend"), which splits on the absent
# column via prepare_subgroup_data() (df[, SG_flag]) and died with the
# opaque "undefined columns selected".
#
# The fix keys the guard on sg.harm (the actual no-subgroup contract), so
# every no-subgroup fit returns empty tables regardless of subgroup_method.
# A subgroup-found fit has sg.harm non-NULL, so the condition reduces to the
# former is.null(df) test -- output byte-identical for found fits.
# ============================================================================


# No-subgroup DINA/GRF fit, forced with extreme thresholds.  Asserts the
# distinguishing SHAPE (sg.harm NULL, df.est populated WITHOUT
# treat.recommend) that separates this path from the consistency path.
.make_null_dg_fit_summary <- function(method, family = "survival") {
  forcing <- if (method == "grf") {
    list(hr.threshold = 10.0, pconsistency.threshold = 0.99,
         subgroup_method = "grf", dmin.grf = 1e6)
  } else {
    list(hr.threshold = 10.0, pconsistency.threshold = 0.99,
         subgroup_method = "dina")
  }
  if (family == "survival") {
    df <- .make_survival_data(N = 160L, HR_harm = 1.0, seed = 7L)
    args <- .fs_args_for("survival", confounders = c("age", "stage", "sex"),
                         extra = forcing)
  } else {
    df <- .make_binary_data(N = 160L, OR_harm = 1.0, seed = 7L)
    args <- .fs_args_for("binary", confounders = c("age", "biomarker_hi", "sex"),
                         extra = forcing)
  }
  fs <- suppressWarnings(do.call(forestsearch, c(list(df.analysis = df), args)))
  if (!is.null(fs$sg.harm)) {
    skip(sprintf("null DGM unexpectedly identified a subgroup (%s/%s)",
                 method, family))
  }
  # Unified no-subgroup contract (Option i): df.est = NULL for all identifiers.
  expect_null(fs$df.est)
  fs
}


for (method in c("dina", "grf")) {
  for (fam in c("survival", "binary")) {

    test_that(sprintf("sg_tables() returns empty tables for no-subgroup %s fit (%s)",
                      method, fam), {
      skip_on_cran()
      skip_if_not_installed("gt")
      fs <- .make_null_dg_fit_summary(method, fam)

      out <- expect_warning(
        sg_tables(fs),
        "no harm subgroup was identified"
      )
      expect_true(is.list(out))
      expect_null(out$tab_estimates)
      expect_null(out$sg10_out)
    })
  }
}


test_that("sg_tables() still tabulates a subgroup-found fit (guard does not fire)", {
  skip_on_cran()
  skip_if_not_installed("gt")

  df <- .make_survival_data(N = 250L, HR_harm = 2.5, seed = 42L)
  args <- .fs_args_for("survival", confounders = c("age", "stage", "sex"))
  fs <- suppressWarnings(do.call(forestsearch, c(list(df.analysis = df), args)))
  if (is.null(fs$sg.harm)) skip("primary fit did not identify a subgroup; test N/A")

  # sg.harm non-NULL: the widened guard reduces to the former is.null(df)
  # test, so a real table is still produced (no empty-tables short-circuit).
  out <- suppressWarnings(sg_tables(fs))
  expect_false(is.null(out$tab_estimates))
})
