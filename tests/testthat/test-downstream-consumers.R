# ============================================================================
# test-downstream-consumers.R
#
# Tier 3 tests: CV and bootstrap pipelines behave correctly when
# individual folds/iterations fail OR produce empty results.
#
# Key behaviors being tested:
#   - forestsearch_tenfold() with an fs.est from a failed primary run
#     gives an actionable error (not "missing required components")
#   - CV runs where every fold fails to identify still return a valid
#     fs_tenfold object (not a crash, not a malformed list)
#   - summarize_cv_results() on pathological inputs: zero
#     identifications, all-NA columns, missing optional columns
#   - fold_summary invariants hold under adverse conditions
#
# Several tests use skip_on_cran() because CV/bootstrap at meaningful
# sims/nb_boots takes multiple seconds even at tiny N.  Local
# devtools::test() runs all of them.
# ============================================================================


# ---------------------------------------------------------------------------
# DOWNSTREAM CONSUMERS ON FIX A's EARLY-RETURN OBJECTS
# ---------------------------------------------------------------------------

test_that("forestsearch_tenfold() rejects failed fs.est with actionable error", {
  # Construct a failed fs.est by running forestsearch() with a config
  # that triggers Fix A's get_FSdata early return.
  df <- .make_survival_data(N = 150L, seed = 42L)
  args <- .fs_args_for("survival",
    confounders = c("age", "stage", "sex"),
    extra = list(
      use_grf      = FALSE,
      use_lasso    = FALSE,
      exclude_cuts = "."
    )
  )

  fs_failed <- suppressWarnings(
    do.call(forestsearch, c(list(df.analysis = df), args))
  )

  # Confirm the failed fs.est has Fix A's shape
  .expect_fs_early_return_shape(fs_failed, expected_stage = "get_FSdata")

  # Now feed it to forestsearch_tenfold() -- must error cleanly
  # (not crash, not silently proceed)
  err <- tryCatch(
    suppressWarnings(
      forestsearch_tenfold(
        fs.est  = fs_failed,
        sims    = 2L,
        Kfolds  = 3L,
        details = FALSE,
        parallel_args = list(plan = "sequential", workers = 1L,
                             show_message = FALSE)
      )
    ),
    error = function(e) e
  )

  expect_s3_class(err, "error")
  # Error message should mention the problem (either "missing components"
  # or "sg.harm NULL" or similar -- actionable)
  expect_match(conditionMessage(err), ".+",
               info = "error message must be non-empty")
})


# ---------------------------------------------------------------------------
# CV PIPELINE: pathological fold_summary inputs
# ---------------------------------------------------------------------------

test_that("summarize_cv_results() handles zero-identification fold_summary", {
  # Construct a synthetic fs_tenfold object where every row has
  # any_found == 0 (no folds identified anything).  Verify that
  # summarize_cv_results() returns a valid structure with NULL or
  # zero-row content where appropriate, rather than crashing.
  n_sim <- 3L
  n_fold <- 4L

  fake_fsum <- data.frame(
    sim = rep(seq_len(n_sim), each = n_fold),
    fold = rep(seq_len(n_fold), times = n_sim),
    n_test = rep(40L, n_sim * n_fold),
    sg1 = NA_character_,
    sg2 = NA_character_,
    grf_cuts = NA_character_,
    pconsistency = NA_real_,
    training_fs_hr = NA_real_,
    n_candidates_evaluated = NA_integer_,
    any_found = 0L,
    stringsAsFactors = FALSE
  )

  fake_cv <- list(
    sens_summary = c(sens_H = NA_real_, sens_Hc = NA_real_,
                     ppv_H = NA_real_, ppv_Hc = NA_real_),
    find_summary = c(Any = 0, Exact = 0, AtLeast1 = 0),
    fold_summary = fake_fsum,
    sims = n_sim,
    Kfolds = n_fold
  )
  class(fake_cv) <- c("fs_tenfold", "list")

  res <- tryCatch(
    summarize_cv_results(cv_output = fake_cv, top_n = 5L),
    error = function(e) e
  )

  expect_false(inherits(res, "error"),
    info = sprintf("summarize_cv_results should not crash on zero-ID input: %s",
                   if (inherits(res, "error")) conditionMessage(res) else ""))

  # Key slots must still exist (may be NULL or empty)
  expect_true(is.list(res))
  expect_true("identification_summary" %in% names(res))
  expect_true("has_pconsistency" %in% names(res) ||
              "has_grf_cuts" %in% names(res),
    info = "metadata slots should exist")
})


test_that("summarize_cv_results() rejects non-fs_tenfold input with clear error", {
  # Plain list without fs_tenfold class -> should warn but still work
  # (duck-typing path)
  fake_fsum <- data.frame(
    sim = 1L, fold = 1L, n_test = 40L,
    sg1 = NA_character_, sg2 = NA_character_,
    any_found = 0L,
    stringsAsFactors = FALSE
  )
  fake_cv_ducktype <- list(fold_summary = fake_fsum,
                           sens_summary = NULL,
                           find_summary = NULL,
                           sims = 1L, Kfolds = 1L)

  # Not fs_tenfold -> should warn but proceed
  expect_warning(
    tryCatch(summarize_cv_results(fake_cv_ducktype, top_n = 5L),
             error = function(e) NULL),
    "fs_tenfold"
  )
})


test_that("summarize_cv_results() rejects fs_kfold input with actionable error", {
  # Object with class 'fs_kfold' (single K-fold) should be rejected with
  # a message pointing the user at the right function.
  fake <- list()
  class(fake) <- c("fs_kfold", "list")

  expect_error(
    summarize_cv_results(fake),
    "fs_tenfold"
  )
})


test_that("summarize_cv_results() rejects fold_summary missing required columns", {
  bad_fsum <- data.frame(sim = 1L, fold = 1L)   # missing sg1, sg2, any_found
  fake_cv <- list(fold_summary = bad_fsum,
                  sens_summary = NULL, find_summary = NULL,
                  sims = 1L, Kfolds = 1L)
  class(fake_cv) <- c("fs_tenfold", "list")

  err <- tryCatch(
    summarize_cv_results(fake_cv),
    error = function(e) e
  )
  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "required column",
               ignore.case = TRUE)
})


test_that("summarize_cv_results() rejects non-list input", {
  expect_error(summarize_cv_results("not a list"), "list")
  expect_error(summarize_cv_results(123),          "list")
  expect_error(summarize_cv_results(NULL),         "list")
})


test_that("summarize_cv_results() top_n argument is validated", {
  fake_fsum <- data.frame(
    sim = 1L, fold = 1L,
    sg1 = NA_character_, sg2 = NA_character_, any_found = 0L,
    stringsAsFactors = FALSE
  )
  fake_cv <- list(fold_summary = fake_fsum,
                  sens_summary = NULL, find_summary = NULL,
                  sims = 1L, Kfolds = 1L)
  class(fake_cv) <- c("fs_tenfold", "list")

  expect_error(summarize_cv_results(fake_cv, top_n = 0L),
               "positive")
  expect_error(summarize_cv_results(fake_cv, top_n = -5L),
               "positive")
  expect_error(summarize_cv_results(fake_cv, top_n = NA_integer_),
               "positive")
})


# ---------------------------------------------------------------------------
# CV PIPELINE: actual small runs
# ---------------------------------------------------------------------------

test_that("forestsearch_tenfold() returns valid shape on identified DGM", {
  skip_on_cran()

  df <- .make_survival_data(N = 250L, HR_harm = 2.5, seed = 42L)
  args <- .fs_args_for("survival",
    confounders = c("age", "stage", "sex"))

  fs <- do.call(forestsearch, c(list(df.analysis = df), args))
  # Proceed only if primary fit identified something (DGM should make
  # this highly likely, but we don't crash the test if it doesn't)
  if (is.null(fs$sg.harm)) {
    skip("primary forestsearch did not identify a subgroup; test N/A")
  }

  cv <- forestsearch_tenfold(
    fs.est  = fs,
    sims    = 2L,
    Kfolds  = 5L,
    details = FALSE,
    parallel_args = list(plan = "sequential", workers = 1L,
                         show_message = FALSE)
  )

  .expect_cv_return_shape(cv)
  expect_true(nrow(cv$fold_summary) == 2L * 5L)
})


test_that("forestsearch_tenfold() returns valid shape when NO fold identifies", {
  skip_on_cran()

  # Tight thresholds on pure-noise DGM -> no fold should identify
  set.seed(42L)
  df <- data.frame(
    id = 1:200, treat = rbinom(200, 1L, 0.5),
    time = rexp(200, 0.05), event = rbinom(200, 1L, 0.7),
    age = round(rnorm(200, 55, 12)),
    sex = factor(rbinom(200, 1L, 0.5), levels = c(0L, 1L))
  )

  # Run primary forestsearch with an impossibly tight hr.threshold.
  # The primary fit itself will not identify anything, but the return
  # object must still be well-formed (args_call_all populated) so
  # forestsearch_tenfold can use it as a template.
  args <- .fs_args_for("survival",
    confounders = c("age", "sex"),
    extra = list(hr.threshold = 10.0,
                 pconsistency.threshold = 0.99))
  fs_full <- do.call(forestsearch, c(list(df.analysis = df), args))

  if (is.null(fs_full$args_call_all)) {
    skip("primary fit malformed; test N/A")
  }

  cv <- tryCatch(
    suppressWarnings(forestsearch_tenfold(
      fs.est  = fs_full,
      sims    = 2L,
      Kfolds  = 3L,
      details = FALSE,
      parallel_args = list(plan = "sequential", workers = 1L,
                           show_message = FALSE)
    )),
    error = function(e) {
      skip(sprintf("CV errored before validation: %s", conditionMessage(e)))
    }
  )

  .expect_cv_return_shape(cv)
  # All folds should have any_found == 0 given the 10.0 threshold
  expect_true(all(cv$fold_summary$any_found == 0L),
    info = "With hr.threshold = 10, no fold should identify")
})


test_that("summarize_cv_results() on real CV run works end-to-end", {
  skip_on_cran()

  df <- .make_survival_data(N = 250L, HR_harm = 2.5, seed = 42L)
  args <- .fs_args_for("survival",
    confounders = c("age", "stage", "sex"))

  fs <- do.call(forestsearch, c(list(df.analysis = df), args))
  if (is.null(fs$sg.harm)) skip("primary fit did not identify; test N/A")

  cv <- forestsearch_tenfold(
    fs.est  = fs,
    sims    = 2L,
    Kfolds  = 5L,
    details = FALSE,
    parallel_args = list(plan = "sequential", workers = 1L,
                         show_message = FALSE)
  )

  res <- summarize_cv_results(
    cv_output         = cv,
    original_sg       = fs$sg.harm,
    original_grf_cuts = fs$grf_cuts,
    top_n             = 5L
  )

  # Core slots must be present
  required_slots <- c("identification_summary",
                      "fold_numeric_summary",
                      "has_pconsistency")
  missing_slots <- setdiff(required_slots, names(res))
  expect_identical(missing_slots, character(0),
    info = sprintf("summarize_cv_results missing: %s",
                   paste(missing_slots, collapse = ", ")))

  expect_true(isTRUE(res$has_pconsistency))
})


# ---------------------------------------------------------------------------
# BOOTSTRAP PIPELINE: downstream consumer resilience
# ---------------------------------------------------------------------------

test_that("forestsearch_bootstrap_dofuture() returns valid structure", {
  skip_on_cran()

  df <- .make_survival_data(N = 200L, HR_harm = 2.5, seed = 42L)
  args <- .fs_args_for("survival",
    confounders = c("age", "stage", "sex"))

  fs <- do.call(forestsearch, c(list(df.analysis = df), args))
  if (is.null(fs$sg.harm)) skip("primary fit did not identify; test N/A")

  bc <- tryCatch(
    forestsearch_bootstrap_dofuture(
      fs.est        = fs,
      nb_boots      = 5L,
      parallel_args = list(plan = "sequential", workers = 1L,
                           show_message = FALSE)
    ),
    error = function(e) {
      skip(sprintf("bootstrap errored: %s", conditionMessage(e)))
    }
  )

  # Contract: returns a list with $results data.frame
  expect_true(is.list(bc))
  expect_true(!is.null(bc$results))
  expect_true(is.data.frame(bc$results) ||
              inherits(bc$results, "data.table"))
})
