# ============================================================================
# test-cv-no-subgroup-dina-grf.R
#
# No-subgroup CV coverage for the DINA and GRF selection modes
# (subgroup_method = "dina" / "grf").  Regression guard for the CV
# base-frame fault these scenarios surfaced:
#
#   Unlike the default consistency path (df.est NULL when no subgroup is
#   identified), .forestsearch_dina_select() and .forestsearch_grf_select()
#   return df.est = df -- a non-NULL frame WITHOUT a treat.recommend
#   column -- on their not-found paths.  .fs_cv_base_frame() used to gate
#   its ITT fallback on is.null(fs.est$df.est) alone, so a no-subgroup
#   DINA/GRF fit took the subgroup-found branch and every CV entry point
#   died with the cryptic 'undefined columns selected'.  The gate now also
#   keys on sg.harm (the actual "no subgroup identified" contract), which
#   leaves every subgroup-found path and every consistency-path fit
#   byte-identical.
#
# No-subgroup forcing uses EXTREME thresholds (hr.threshold = 10 for DINA's
# harm floor; dmin.grf = 1e6 for GRF's DR-score floor) so GRF's known
# near-threshold floating-point flip sensitivity cannot produce a find.
# ============================================================================


.cv_pa_dg <- list(plan = "sequential", workers = 1L, show_message = FALSE)

# Null-DGM no-subgroup fit for subgroup_method = "dina" / "grf".
# Asserts the no-subgroup premise (sg.harm NULL) before returning.
.make_null_fit_dg <- function(method, family) {
  forcing <- if (method == "grf") {
    list(subgroup_method = "grf", dmin.grf = 1e6)
  } else {
    list(subgroup_method = "dina", hr.threshold = 10.0)
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
  fs
}


for (method in c("dina", "grf")) {
  for (fam in c("survival", "binary")) {

    test_that(sprintf("tenfold and Kfold complete for no-subgroup %s fit (%s)",
                      method, fam), {
      skip_on_cran()
      fs <- .make_null_fit_dg(method, fam)

      # The shape this regression guards: DINA/GRF not-found fits carry a
      # NON-NULL df.est without treat.recommend (the consistency path
      # carries NULL).  Either shape must be handled; assert the current
      # one so a future contract change is caught deliberately.
      expect_false(is.null(fs$df.est))
      expect_false("treat.recommend" %in% names(fs$df.est))

      cv <- suppressWarnings(forestsearch_tenfold(
        fs.est = fs, sims = 2L, Kfolds = 3L,
        details = FALSE, parallel_args = .cv_pa_dg
      ))
      .expect_cv_return_shape(cv)
      expect_true(all(cv$fold_summary$any_found == 0L))
      expect_true(all(cv$find_summary == 0))
      expect_true(all(is.na(cv$sens_summary)))

      kf <- suppressWarnings(forestsearch_Kfold(
        fs.est = fs, Kfolds = 3L, parallel_args = .cv_pa_dg
      ))
      expect_s3_class(kf, "fs_kfold")
      # No original subgroup: every subject's original recommendation is
      # ITT, over the FULL analysis frame
      expect_true(all(kf$resCV$treat.recommend.original == 1.0))
      expect_identical(
        nrow(kf$resCV),
        nrow(as.data.frame(fs$args_call_all$df.analysis))
      )
    })
  }
}
