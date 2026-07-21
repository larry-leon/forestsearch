# ============================================================================
# test-no-subgroup-cv.R
#
# A no-subgroup primary fit (df.est NULL, sg.harm NULL) is a valid,
# completable input across the CV surface -- forestsearch_tenfold(),
# forestsearch_Kfold(), and forestsearch_KfoldOut() (outall FALSE and
# TRUE) -- for EVERY outcome family: survival, binary, continuous, count.
#
# Continuous and count are the critical families here: forestsearch()
# synthesizes a placeholder event column (".event_placeholder" = 1L) that
# exists in df.est but NOT in the raw analysis frame, so the no-subgroup
# fallback must re-create it (.fs_cv_base_frame()).  Count additionally
# exercises offset.name alongside the placeholder.
#
# Also covers the survival empty-arm summary-table guard: an all-ITT CV
# prediction frame yields a zero-member complementary arm, which
# analyze_subgroup() must summarize as an NA row, not a survfit() error.
#
# All tests use skip_on_cran(): each runs a full (small) forestsearch fit
# plus CV, taking several seconds.  Local devtools::test() runs them all.
# ============================================================================


# Null-DGM fit for a family: no real treatment-effect subgroup, and
# thresholds tight enough that nothing is identified.  Asserts the
# no-subgroup premise (df.est / sg.harm NULL) before returning.
.make_null_fit <- function(family) {
  null_extra <- list(hr.threshold = 10.0, pconsistency.threshold = 0.99)
  if (family == "survival") {
    df <- .make_survival_data(N = 160L, HR_harm = 1.0, seed = 7L)
    args <- .fs_args_for("survival", confounders = c("age", "stage", "sex"),
                         extra = null_extra)
  } else if (family == "binary") {
    df <- .make_binary_data(N = 160L, OR_harm = 1.0, seed = 7L)
    args <- .fs_args_for("binary", confounders = c("age", "biomarker_hi", "sex"),
                         extra = null_extra)
  } else if (family == "continuous") {
    df <- .make_continuous_data(N = 160L, MD_harm = 0.0, seed = 7L)
    args <- .fs_args_for("continuous", confounders = c("age", "biomarker_hi", "sex"),
                         extra = null_extra)
  } else {
    df <- .make_count_data(N = 180L, IRR_harm = 1.0, seed = 7L)
    args <- .fs_args_for("count", confounders = c("age", "biomarker_hi", "sex"),
                         extra = null_extra)
  }
  fs <- suppressWarnings(do.call(forestsearch, c(list(df.analysis = df), args)))
  if (!is.null(fs$sg.harm)) {
    skip(sprintf("null DGM unexpectedly identified a subgroup (%s)", family))
  }
  fs
}

.cv_pa <- list(plan = "sequential", workers = 1L, show_message = FALSE)


for (fam in c("survival", "binary", "continuous", "count")) {

  test_that(sprintf("tenfold completes as 0%%-identification for no-subgroup fit (%s)", fam), {
    skip_on_cran()
    fs <- .make_null_fit(fam)
    expect_null(fs$df.est)

    cv <- suppressWarnings(forestsearch_tenfold(
      fs.est = fs, sims = 2L, Kfolds = 3L,
      details = FALSE, parallel_args = .cv_pa
    ))
    .expect_cv_return_shape(cv)
    expect_true(all(cv$fold_summary$any_found == 0L))
    expect_true(all(cv$find_summary == 0))
    expect_true(all(is.na(cv$sens_summary)))
  })

  test_that(sprintf("Kfold completes for no-subgroup fit (%s)", fam), {
    skip_on_cran()
    fs <- .make_null_fit(fam)

    kf <- suppressWarnings(forestsearch_Kfold(
      fs.est = fs, Kfolds = 3L, parallel_args = .cv_pa
    ))
    expect_s3_class(kf, "fs_kfold")
    # No original subgroup: every subject's original recommendation is ITT
    expect_true(all(kf$resCV$treat.recommend.original == 1.0))
    # resCV covers the full analysis frame (placeholder column did not
    # silently drop subjects for continuous/count)
    expect_identical(
      nrow(kf$resCV),
      nrow(as.data.frame(fs$args_call_all$df.analysis))
    )
  })

  test_that(sprintf("KfoldOut handles no-subgroup CV result, outall FALSE and TRUE (%s)", fam), {
    skip_on_cran()
    fs <- .make_null_fit(fam)
    kf <- suppressWarnings(forestsearch_Kfold(
      fs.est = fs, Kfolds = 3L, parallel_args = .cv_pa
    ))

    out_f <- forestsearch_KfoldOut(kf, outall = FALSE)
    expect_named(out_f$find_metrics,
                 c("Any", "Exact", "At least 1", "Cov1", "Cov2",
                   "Cov 1 & 2", "Cov1 exact", "Cov2 exact"))
    expect_named(out_f$sens_metrics_original,
                 c("sens_H", "sens_Hc", "ppv_H", "ppv_Hc"))

    out_t <- forestsearch_KfoldOut(kf, outall = TRUE)
    # ITT-only tables: subgroup tables absent, ITT table doubles as tab_all
    expect_null(out_t$SG_tab_original)
    expect_null(out_t$SG_tab_Kfold)
    expect_identical(out_t$tab_all, out_t$itt_tab)
  })
}


test_that("SG_tab_estimates returns NA row (not error) for empty subgroup arm", {
  set.seed(1)
  df <- data.frame(
    tte   = rexp(80, 0.05),
    event = rbinom(80, 1L, 0.7),
    treat = rbinom(80, 1L, 0.5),
    treat.recommend = 1.0            # all-ITT: SG == 0 arm is empty
  )
  tab <- SG_tab_estimates(
    df, SG_flag = "treat.recommend",
    outcome.name = "tte", event.name = "event", treat.name = "treat",
    sg1_name = "Recommend", sg0_name = "Not recommend"
  )
  expect_identical(nrow(tab), 2L)
  # Empty arm: zero counts, NA estimates; both rows same column set
  expect_match(tab[1L, "n"], "^0 ")
  expect_true(is.na(tab[1L, "HR (95% CI)"]))
  # Non-empty arm still estimated
  expect_false(is.na(tab[2L, "HR (95% CI)"]))
})
