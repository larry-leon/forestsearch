# ============================================================================
# test-cv-no-subgroup-edges.R
#
# Tier-2 edge-case regression coverage for the no-subgroup CV path:
#
#   T2.4  est.scale = "1/hr": forestsearch_Kfold() SECTION-7 flips the
#         treatment indicator (treat2 = 1 - treat) and swaps the subgroup
#         labels; the NULL-df.est ITT default (treat.recommend.original =
#         1.0) must pass through that flip unharmed.
#
#   T2.5  adjust_covariates on the no-subgroup fallback: the CV base
#         frame is rebuilt from the raw analysis frame, and must retain
#         the covariate-adjustment columns -- including the busiest
#         compound case (continuous outcome + synthesized placeholder
#         event column + adjustment columns).
#
#   T2.6  Leave-one-out scale: the original no-subgroup bug surfaced at
#         LOO (Kfolds = N); both the clean completion AND the fatal
#         refit-error path (with correct LOO fold naming) must hold there,
#         not just at Kfolds = 3-5.
#
#   T2.7  Parallel/sequential equivalence: same seed, plan = "multisession"
#         (2 workers) vs "sequential" must produce identical() tenfold
#         no-subgroup results (fold_summary + metrics).  Requires the
#         package to be INSTALLED (multisession workers resolve
#         forestsearch functions from the installed namespace).
# ============================================================================


.cv_pa_t2 <- list(plan = "sequential", workers = 1L, show_message = FALSE)

.sens_names_t2 <- c("sens_H", "sens_Hc", "ppv_H", "ppv_Hc")
.find_names_t2 <- c("Any", "Exact", "At least 1", "Cov1", "Cov2",
                    "Cov 1 & 2", "Cov1 exact", "Cov2 exact")

# Null-DGM no-subgroup fit (consistency path); premise asserted.
.make_null_fit_t2 <- function(family = "survival", N = 160L, extra = list()) {
  null_extra <- utils::modifyList(
    list(hr.threshold = 10.0, pconsistency.threshold = 0.99), extra)
  if (family == "survival") {
    df <- .make_survival_data(N = N, HR_harm = 1.0, seed = 7L)
    args <- .fs_args_for("survival", confounders = c("age", "stage", "sex"),
                         extra = null_extra)
  } else {
    df <- .make_continuous_data(N = N, MD_harm = 0.0, seed = 7L)
    args <- .fs_args_for("continuous",
                         confounders = c("age", "biomarker_hi", "sex"),
                         extra = null_extra)
  }
  fs <- suppressWarnings(do.call(forestsearch, c(list(df.analysis = df), args)))
  if (!is.null(fs$sg.harm)) {
    skip(sprintf("null DGM unexpectedly identified a subgroup (%s)", family))
  }
  fs
}


# ----------------------------------------------------------------------------
# T2.4: est.scale = "1/hr" + no-subgroup
# ----------------------------------------------------------------------------

test_that("est.scale='1/hr' + no-subgroup: Kfold completes through the flip", {
  skip_on_cran()

  fs <- .make_null_fit_t2("survival", extra = list(est.scale = "1/hr"))
  expect_null(fs$df.est)
  expect_identical(fs$args_call_all$est.scale, "1/hr")

  kf <- suppressWarnings(forestsearch_Kfold(
    fs.est = fs, Kfolds = 3L, parallel_args = .cv_pa_t2
  ))
  expect_s3_class(kf, "fs_kfold")

  # SECTION-7 flip applied: treat2 = 1 - treat, labels swapped
  treat.name <- fs$args_call_all$treat.name
  expect_true("treat2" %in% names(kf$resCV))
  expect_identical(kf$resCV$treat2, 1 - kf$resCV[[treat.name]])
  expect_identical(kf$sg0.name, "Recommend")
  expect_identical(kf$sg1.name, "Not recommend")

  # The ITT default survives the flip: all-1 original recommendation over
  # the full analysis frame, and the no-subgroup metric shape is intact
  expect_true(all(kf$resCV$treat.recommend.original == 1.0))
  expect_identical(nrow(kf$resCV),
                   nrow(as.data.frame(fs$args_call_all$df.analysis)))
  expect_identical(names(kf$sens_summary), .sens_names_t2)
  expect_identical(names(kf$find_summary), .find_names_t2)
  expect_true(all(is.na(kf$sens_summary)))
  expect_true(all(kf$find_summary == 0))
})


# ----------------------------------------------------------------------------
# T2.5: adjust_covariates on the no-subgroup fallback
# ----------------------------------------------------------------------------

test_that("adjust_covariates + no-subgroup: Kfold completes, column carried (survival)", {
  skip_on_cran()

  fs <- .make_null_fit_t2("survival", extra = list(adjust_covariates = "noise"))
  expect_null(fs$df.est)

  kf <- suppressWarnings(forestsearch_Kfold(
    fs.est = fs, Kfolds = 3L, parallel_args = .cv_pa_t2
  ))
  expect_s3_class(kf, "fs_kfold")
  expect_true("noise" %in% names(kf$resCV))
  expect_true(all(kf$resCV$treat.recommend.original == 1.0))
  expect_identical(nrow(kf$resCV),
                   nrow(as.data.frame(fs$args_call_all$df.analysis)))
})


test_that("continuous + placeholder event + adjust_covariates carries through CV", {
  skip_on_cran()

  fs <- .make_null_fit_t2("continuous",
                          extra = list(adjust_covariates = "biomarker"))
  expect_null(fs$df.est)

  # The no-subgroup fallback frame must re-create the placeholder event
  # column AND retain the adjustment column
  event.name <- fs$args_call_all$event.name
  bf <- forestsearch:::.fs_cv_base_frame(fs)
  expect_true("biomarker" %in% names(bf))
  expect_true(event.name %in% names(bf))
  expect_true(all(bf[[event.name]] == 1L))
  expect_true(all(bf$treat.recommend.original == 1.0))

  kf <- suppressWarnings(forestsearch_Kfold(
    fs.est = fs, Kfolds = 3L, parallel_args = .cv_pa_t2
  ))
  expect_s3_class(kf, "fs_kfold")
  expect_true("biomarker" %in% names(kf$resCV))
  expect_identical(nrow(kf$resCV),
                   nrow(as.data.frame(fs$args_call_all$df.analysis)))

  # tenfold's per-sim CV frames must also carry the adjustment column
  tf <- suppressWarnings(forestsearch_tenfold(
    fs.est = fs, sims = 2L, Kfolds = 3L, details = FALSE,
    parallel_args = .cv_pa_t2, keep_resCV = TRUE
  ))
  .expect_cv_return_shape(tf)
  expect_true(all(tf$fold_summary$any_found == 0L))
  expect_true(all(vapply(tf$resCV_all,
                         function(d) "biomarker" %in% names(d),
                         logical(1))))
})


# ----------------------------------------------------------------------------
# T2.6: leave-one-out scale
# ----------------------------------------------------------------------------

test_that("LOO-scale (Kfolds = N) no-subgroup CV completes", {
  skip_on_cran()

  N <- 40L
  fs <- .make_null_fit_t2("survival", N = N, extra = list(n.min = 12L))
  expect_null(fs$df.est)

  kf <- suppressWarnings(forestsearch_Kfold(
    fs.est = fs, Kfolds = N, parallel_args = .cv_pa_t2
  ))
  expect_s3_class(kf, "fs_kfold")
  expect_identical(kf$Kfolds, N)
  expect_identical(nrow(kf$resCV), N)
  # The PRIMARY fit found no subgroup, so every subject's ORIGINAL
  # (full-data) recommendation is ITT -- this is the no-subgroup contract
  # the LOO-scale fallback must honour, and where the original bug bit.
  expect_true(all(kf$resCV$treat.recommend.original == 1.0))
  # (Individual n-1 training folds MAY still identify a subgroup at LOO
  # scale; that is expected live-machinery behaviour, not a contract
  # violation, so it is deliberately NOT asserted here.)
  expect_identical(names(kf$sens_summary), .sens_names_t2)
  expect_identical(names(kf$find_summary), .find_names_t2)
})


test_that("LOO-scale refit failure still names the correct fold", {
  skip_on_cran()

  N <- 40L
  fs <- .make_null_fit_t2("survival", N = N, extra = list(n.min = 12L))

  real_fs <- forestsearch
  calls <- 0L
  testthat::local_mocked_bindings(
    forestsearch = function(...) {
      calls <<- calls + 1L
      if (calls == 23L) stop("induced LOO refit failure XYZZY")
      real_fs(...)
    }
  )

  err <- tryCatch(
    suppressWarnings(forestsearch_Kfold(fs.est = fs, Kfolds = N,
                                        parallel_args = .cv_pa_t2)),
    error = function(e) e)
  expect_s3_class(err, "error")
  msg <- conditionMessage(err)
  fold_idx <- as.integer(sub(".*FAILED on fold (\\d+) of (\\d+).*", "\\1", msg))
  fold_tot <- as.integer(sub(".*FAILED on fold (\\d+) of (\\d+).*", "\\2", msg))
  expect_identical(fold_idx, 23L)
  expect_identical(fold_tot, 40L)
  expect_match(msg, "induced LOO refit failure XYZZY", fixed = TRUE)
})


# ----------------------------------------------------------------------------
# T2.7: parallel vs sequential equivalence on the no-subgroup path
# ----------------------------------------------------------------------------

test_that("multisession and sequential tenfold agree on the no-subgroup path", {
  skip_on_cran()

  # forestsearch_tenfold()'s foreach body references package-internal
  # helpers (e.g. .fs_adjust_vars) unqualified.  Against an installed
  # package a multisession worker resolves them from its loaded namespace;
  # under devtools::test()/load_all() the parent runs the DEV namespace but
  # each spawned worker loads the INSTALLED package, so the serialised
  # foreach closure (whose environment is the dev namespace) fails with
  # "could not find function".  That is a dev-load harness artifact, NOT a
  # code fault: verified separately that against the installed package the
  # sequential and multisession no-subgroup results are identical().  Skip
  # cleanly when dev-loaded so devtools::test() stays green; the lock still
  # runs under R CMD check on the built (installed) package.
  if (isTRUE(tryCatch(
        requireNamespace("pkgload", quietly = TRUE) &&
          pkgload::is_dev_package("forestsearch"),
        error = function(e) FALSE))) {
    skip("dev-loaded (devtools::test); multisession workers load the installed pkg -- lock verified on installed builds")
  }

  fs <- .make_null_fit_t2("survival")

  run_tf <- function(pa) {
    set.seed(20260720L)
    suppressWarnings(forestsearch_tenfold(
      fs.est = fs, sims = 3L, Kfolds = 3L,
      details = FALSE, parallel_args = pa
    ))
  }

  tf_seq <- run_tf(list(plan = "sequential", workers = 1L,
                        show_message = FALSE))
  tf_par <- run_tf(list(plan = "multisession", workers = 2L,
                        show_message = FALSE))

  # Everything except timing must agree exactly (locks parallel-RNG
  # reproducibility for the no-subgroup path)
  expect_identical(tf_par$fold_summary, tf_seq$fold_summary)
  expect_identical(tf_par$sens_out,     tf_seq$sens_out)
  expect_identical(tf_par$find_out,     tf_seq$find_out)
  expect_identical(tf_par$sens_summary, tf_seq$sens_summary)
  expect_identical(tf_par$find_summary, tf_seq$find_summary)
  expect_identical(tf_par$sims,         tf_seq$sims)
  expect_identical(tf_par$Kfolds,       tf_seq$Kfolds)
})
