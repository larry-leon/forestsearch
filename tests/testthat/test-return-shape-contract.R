# ============================================================================
# test-return-shape-contract.R
#
# Tier 1 tests: forestsearch() must return objects satisfying the
# downstream CV/bootstrap consumer contract on EVERY code path, whether
# successful or failed.
#
# Key contract:
#   Success  -> object with $args_call_all, $df.est, $sg.harm (possibly NULL),
#               $outcome_type, $effect_measure, $threshold_config.
#   Failure  -> object with $args_call_all, $sg.harm = NULL, $error_log.
#
# Missing $args_call_all on ANY return path is a regression (this was
# Fix A in v0.2.0-dev; see forestsearch_main.R line ~1172).
# ============================================================================


# ---------------------------------------------------------------------------
# HAPPY PATHS: object should satisfy full return contract
# ---------------------------------------------------------------------------

test_that("survival happy path returns full consumer contract", {
  df <- .make_survival_data(N = 200L, HR_harm = 2.0, seed = 42L)
  args <- .fs_args_for("survival",
    confounders = c("age", "stage", "sex", "noise"))

  fs <- do.call(forestsearch, c(list(df.analysis = df), args))

  .expect_fs_full_return(fs)
  expect_identical(fs$outcome_type, "survival")
  expect_true(inherits(fs, "forestsearch"))
})


test_that("binary GLM happy path returns full consumer contract", {
  df <- .make_binary_data(N = 200L, OR_harm = 2.5, seed = 42L)
  args <- .fs_args_for("binary",
    confounders = c("age", "wtkg", "biomarker", "biomarker_hi", "sex"))

  fs <- do.call(forestsearch, c(list(df.analysis = df), args))

  .expect_fs_full_return(fs)
  expect_identical(fs$outcome_type,   "binary")
  expect_identical(fs$effect_measure, "OR")
})


test_that("continuous GLM happy path returns full consumer contract", {
  df <- .make_continuous_data(N = 200L, MD_harm = 1.5, seed = 42L)
  args <- .fs_args_for("continuous",
    confounders = c("age", "biomarker", "biomarker_hi", "sex"))

  fs <- do.call(forestsearch, c(list(df.analysis = df), args))

  .expect_fs_full_return(fs)
  expect_identical(fs$outcome_type,   "continuous")
  expect_identical(fs$effect_measure, "MD")
})


test_that("count GLM with offset happy path returns full consumer contract", {
  df <- .make_count_data(N = 250L, IRR_harm = 2.0, seed = 42L)
  args <- .fs_args_for("count",
    confounders = c("age", "biomarker", "biomarker_hi", "sex"))

  fs <- do.call(forestsearch, c(list(df.analysis = df), args))

  .expect_fs_full_return(fs)
  expect_identical(fs$outcome_type,   "count")
  expect_identical(fs$effect_measure, "IRR")
})


# ---------------------------------------------------------------------------
# SUCCESSFUL-BUT-NO-SUBGROUP PATH: sg.harm = NULL, everything else valid
# ---------------------------------------------------------------------------

test_that("successful run with no subgroup identified returns full contract", {
  # Pure-noise DGM: no treatment-covariate interaction exists
  set.seed(42L)
  df <- data.frame(
    id    = 1:200,
    treat = rbinom(200, 1L, 0.5),
    time  = rexp(200, 0.05),
    event = rbinom(200, 1L, 0.7),
    age   = round(rnorm(200, 55, 12)),
    sex   = factor(rbinom(200, 1L, 0.5), levels = c(0L, 1L))
  )
  # Use TIGHT thresholds that shouldn't be met by pure noise
  args <- .fs_args_for("survival",
    confounders = c("age", "sex"),
    extra = list(hr.threshold = 5.0, pconsistency.threshold = 0.99)
  )

  fs <- do.call(forestsearch, c(list(df.analysis = df), args))

  # Full contract still satisfied even though no subgroup found
  .expect_fs_full_return(fs)
  # Specifically: sg.harm is NULL, but the object is complete
  expect_null(fs$sg.harm)
})


# ---------------------------------------------------------------------------
# get_FSdata FAILURE PATH: Fix A regression test
# ---------------------------------------------------------------------------

test_that("get_FSdata failure early-return includes args_call_all [Fix A]", {
  # Force get_FSdata to fail deterministically by passing exclude_cuts
  # that matches every candidate cut.  This is cleaner than trying to
  # construct adversarial covariates, because it triggers the exact
  # n_confs == 0 path (Fix B's territory) -> get_FSdata throws ->
  # forestsearch() catches and returns its early-return object (Fix A's
  # territory).
  df <- .make_survival_data(N = 150L, seed = 42L)
  args <- .fs_args_for("survival",
    confounders = c("age", "stage", "sex"),
    extra = list(
      use_grf      = FALSE,
      use_lasso    = FALSE,
      exclude_cuts = "."     # matches everything
    )
  )

  # Run capturing warnings (the failure path warns before returning)
  res <- suppressWarnings(
    .run_fs_capture(df, args)
  )
  fs <- res$value

  # --- FIX A regression: args_call_all MUST be present on this path
  .expect_fs_early_return_shape(fs, expected_stage = "get_FSdata")

  # --- Message hints that the failure was at get_FSdata
  expect_match(fs$error_log$message, "FSdata|get_FSdata",
               ignore.case = TRUE)
})


# ---------------------------------------------------------------------------
# Return-object class
# ---------------------------------------------------------------------------

test_that("successful return has class 'forestsearch'", {
  df <- .make_binary_data(N = 150L, seed = 42L)
  args <- .fs_args_for("binary",
    confounders = c("age", "biomarker", "biomarker_hi", "sex"))

  fs <- do.call(forestsearch, c(list(df.analysis = df), args))

  expect_s3_class(fs, "forestsearch")
})


# ---------------------------------------------------------------------------
# outcome_type / effect_measure are preserved on return
# ---------------------------------------------------------------------------

test_that("outcome_type round-trips through forestsearch()", {
  for (otype in c("survival", "binary", "continuous", "count")) {
    df <- switch(otype,
      survival   = .make_survival_data(N = 150L),
      binary     = .make_binary_data(N = 150L),
      continuous = .make_continuous_data(N = 150L),
      count      = .make_count_data(N = 200L)
    )
    conf <- switch(otype,
      survival   = c("age", "stage", "sex", "noise"),
      c("age", "biomarker", "biomarker_hi", "sex")
    )
    args <- .fs_args_for(otype, confounders = conf)

    fs <- do.call(forestsearch, c(list(df.analysis = df), args))

    expect_identical(fs$outcome_type, otype,
      info = sprintf("outcome_type = %s should be preserved", otype))
  }
})


test_that("effect_measure round-trips for GLM outcomes", {
  for (otype in c("binary", "continuous", "count")) {
    em <- switch(otype, binary = "OR", continuous = "MD", count = "IRR")
    df <- switch(otype,
      binary     = .make_binary_data(N = 150L),
      continuous = .make_continuous_data(N = 150L),
      count      = .make_count_data(N = 200L)
    )
    args <- .fs_args_for(otype,
      confounders = c("age", "biomarker", "biomarker_hi", "sex"))

    fs <- do.call(forestsearch, c(list(df.analysis = df), args))

    expect_identical(fs$effect_measure, em,
      info = sprintf("effect_measure = %s should be preserved for %s",
                     em, otype))
  }
})


# ---------------------------------------------------------------------------
# Consistency: early-return shape mirrors subgroup.search early-return shape
# ---------------------------------------------------------------------------

test_that("all early-return paths use the same list shape", {
  # This is a structural test: the early-return from get_FSdata (line 1172
  # in forestsearch_main.R, patched by Fix A) and subgroup.search (line
  # 1356) should produce lists with the same minimum fields.
  #
  # We exercise get_FSdata failure via exclude_cuts = "." and verify the
  # shape fields are all present.  We don't have a clean way to force
  # only subgroup.search to fail at this test tier (requires get_FSdata
  # to succeed but subgroup.search to error), but the structural check
  # is encoded in .expect_fs_early_return_shape().
  df <- .make_survival_data(N = 100L, seed = 42L)
  args <- .fs_args_for("survival",
    confounders = c("age", "stage", "sex"),
    extra = list(
      use_grf      = FALSE,
      use_lasso    = FALSE,
      exclude_cuts = "."
    )
  )

  fs <- suppressWarnings(
    do.call(forestsearch, c(list(df.analysis = df), args))
  )

  # All three of these must be present, regardless of which stage failed
  required_fields <- c("sg.harm", "args_call_all", "error_log")
  missing_fields <- setdiff(required_fields, names(fs))
  expect_identical(missing_fields, character(0),
    info = sprintf("missing fields on early-return: %s",
                   paste(missing_fields, collapse = ", ")))
})
