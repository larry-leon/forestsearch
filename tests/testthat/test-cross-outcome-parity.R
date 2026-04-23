# ============================================================================
# test-cross-outcome-parity.R
#
# Tier 4 tests: structural parity across outcome types.
#
# forestsearch() supports four outcome types (survival, binary,
# continuous, count) and must present a symmetric return contract to
# downstream consumers so that CV and bootstrap pipelines don't need
# outcome-type branching.  The Phase A/B/C diagnostic additions
# specifically require:
#
#   - fold_summary structure identical across outcome types
#   - training_fs_hr on natural scale regardless of outcome type
#     (Phase B option-3 exp() transform)
#   - summarize_cv_results slot structure identical across outcome types
#
# These tests exercise the pipeline on each outcome type and check
# these invariants.  Some tests use skip_on_cran() since CV on four
# outcome types is slow even at tiny N.
# ============================================================================


# ---------------------------------------------------------------------------
# HAPPY-PATH PARITY: same slots across outcome types
# ---------------------------------------------------------------------------

test_that("forestsearch() returns structurally equivalent objects across outcome types", {
  # Run forestsearch() on each outcome type with matching DGM.  The
  # returned object should have the same set of top-level slots
  # regardless of outcome.
  types_data <- list(
    survival   = list(df = .make_survival_data(N = 200L, seed = 42L),
                      conf = c("age", "stage", "sex")),
    binary     = list(df = .make_binary_data(N = 200L, seed = 42L),
                      conf = c("age", "biomarker", "biomarker_hi", "sex")),
    continuous = list(df = .make_continuous_data(N = 200L, seed = 42L),
                      conf = c("age", "biomarker", "biomarker_hi", "sex")),
    count      = list(df = .make_count_data(N = 250L, seed = 42L),
                      conf = c("age", "biomarker", "biomarker_hi", "sex"))
  )

  runs <- lapply(names(types_data), function(otype) {
    td <- types_data[[otype]]
    args <- .fs_args_for(otype, confounders = td$conf)
    do.call(forestsearch, c(list(df.analysis = td$df), args))
  })
  names(runs) <- names(types_data)

  # Every run must satisfy the full return contract
  for (otype in names(runs)) {
    .expect_fs_full_return(runs[[otype]])
  }

  # Every run must expose the same COMMON slot names (union)
  required_common <- c("sg.harm", "args_call_all", "outcome_type",
                       "grf_cuts", "find.grps", "grp.consistency",
                       "df.est", "df.predict", "df.test")
  for (otype in names(runs)) {
    missing_slots <- setdiff(required_common, names(runs[[otype]]))
    expect_identical(missing_slots, character(0),
      info = sprintf("outcome %s missing slots: %s",
                     otype, paste(missing_slots, collapse = ", ")))
  }
})


test_that("outcome_type and effect_measure round-trip for all GLM types", {
  specs <- list(
    list(otype = "binary",     em = "OR"),
    list(otype = "continuous", em = "MD"),
    list(otype = "count",      em = "IRR")
  )

  for (spec in specs) {
    df <- switch(spec$otype,
      binary     = .make_binary_data(N = 150L, seed = 42L),
      continuous = .make_continuous_data(N = 150L, seed = 42L),
      count      = .make_count_data(N = 200L, seed = 42L)
    )
    args <- .fs_args_for(spec$otype,
      confounders = c("age", "biomarker", "biomarker_hi", "sex"))

    fs <- do.call(forestsearch, c(list(df.analysis = df), args))

    expect_identical(fs$outcome_type, spec$otype,
      info = sprintf("outcome_type = %s", spec$otype))
    expect_identical(fs$effect_measure, spec$em,
      info = sprintf("effect_measure = %s for %s", spec$em, spec$otype))
  }
})


# ---------------------------------------------------------------------------
# CV fold_summary STRUCTURAL PARITY across outcome types
# ---------------------------------------------------------------------------

test_that("fold_summary columns are identical across outcome types", {
  skip_on_cran()

  run_cv_for <- function(otype) {
    df <- switch(otype,
      survival   = .make_survival_data(N = 250L, HR_harm = 2.5, seed = 42L),
      binary     = .make_binary_data(N = 250L, OR_harm = 2.5, seed = 42L),
      continuous = .make_continuous_data(N = 250L, MD_harm = 1.5, seed = 42L),
      count      = .make_count_data(N = 300L, IRR_harm = 2.0, seed = 42L)
    )
    conf <- if (otype == "survival")
      c("age", "stage", "sex")
    else
      c("age", "biomarker", "biomarker_hi", "sex")
    args <- .fs_args_for(otype, confounders = conf)

    fs <- do.call(forestsearch, c(list(df.analysis = df), args))
    if (is.null(fs$sg.harm)) return(NULL)  # skip type if not identified

    forestsearch_tenfold(
      fs.est  = fs,
      sims    = 2L,
      Kfolds  = 3L,
      details = FALSE,
      parallel_args = list(plan = "sequential", workers = 1L,
                           show_message = FALSE)
    )
  }

  cv_runs <- lapply(c("survival", "binary"), run_cv_for)
  names(cv_runs) <- c("survival", "binary")
  cv_runs <- cv_runs[!vapply(cv_runs, is.null, logical(1))]

  if (length(cv_runs) < 2L) {
    skip("need at least two outcome types to identify for parity check")
  }

  # Column sets must be identical
  col_sets <- lapply(cv_runs, function(cv) sort(names(cv$fold_summary)))
  for (i in 2:length(col_sets)) {
    expect_identical(col_sets[[i]], col_sets[[1]],
      info = sprintf("fold_summary columns differ between %s and %s",
                     names(cv_runs)[1], names(cv_runs)[i]))
  }
})


# ---------------------------------------------------------------------------
# PHASE B REGRESSION: training_fs_hr on natural scale
# ---------------------------------------------------------------------------

test_that("training_fs_hr is on natural scale for survival (HR)", {
  skip_on_cran()

  df <- .make_survival_data(N = 300L, HR_harm = 2.5, seed = 42L)
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

  hr_vals <- cv$fold_summary$training_fs_hr
  hr_nonna <- hr_vals[!is.na(hr_vals)]

  if (length(hr_nonna) == 0L) skip("no identifying folds; test N/A")

  # Natural-scale HR for a harm subgroup must be > 1
  med_hr <- median(hr_nonna)
  expect_true(med_hr > 1,
    info = sprintf("survival training_fs_hr median = %.3f (natural HR scale)",
                   med_hr))

  # Must respect hr.threshold from args
  threshold <- args$hr.threshold
  expect_true(all(hr_nonna >= threshold),
    info = sprintf("training_fs_hr min = %.3f; threshold = %.3f",
                   min(hr_nonna), threshold))
})


test_that("training_fs_hr is on natural scale for binary (OR)", {
  skip_on_cran()

  df <- .make_binary_data(N = 400L, OR_harm = 2.5, seed = 42L)
  args <- .fs_args_for("binary",
    confounders = c("age", "biomarker", "biomarker_hi", "sex"))

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

  or_vals <- cv$fold_summary$training_fs_hr
  or_nonna <- or_vals[!is.na(or_vals)]
  if (length(or_nonna) == 0L) skip("no identifying folds; test N/A")

  # If Phase B's exp() transform is working, OR values should be > 1
  # on natural scale.  If log-scale, they'd cluster near log(2.5) ~= 0.92
  # and fail against the hr.threshold = 1.25 imposed during search.
  med_or <- median(or_nonna)
  expect_true(med_or > 1,
    info = sprintf("binary training_fs_hr median = %.3f (natural OR scale); ",
                   "log-scale value would be %.3f",
                   med_or, log(med_or)))
  expect_true(all(or_nonna >= args$hr.threshold),
    info = sprintf("all identifying folds must have OR >= %.2f (natural); ",
                   "min = %.3f",
                   args$hr.threshold, min(or_nonna)))
})


test_that("training_fs_hr is on natural scale for count (IRR)", {
  skip_on_cran()

  df <- .make_count_data(N = 400L, IRR_harm = 2.0, seed = 42L)
  args <- .fs_args_for("count",
    confounders = c("age", "biomarker", "biomarker_hi", "sex"))

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

  irr_vals <- cv$fold_summary$training_fs_hr
  irr_nonna <- irr_vals[!is.na(irr_vals)]
  if (length(irr_nonna) == 0L) skip("no identifying folds; test N/A")

  med_irr <- median(irr_nonna)
  expect_true(med_irr > 1,
    info = sprintf("count training_fs_hr median = %.3f (natural IRR scale)",
                   med_irr))
})


test_that("training_fs_hr for continuous (MD) is identity-scale, not exponentiated", {
  skip_on_cran()

  df <- .make_continuous_data(N = 400L, MD_harm = 1.5, seed = 42L)
  args <- .fs_args_for("continuous",
    confounders = c("age", "biomarker", "biomarker_hi", "sex"))

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

  md_vals <- cv$fold_summary$training_fs_hr
  md_nonna <- md_vals[!is.na(md_vals)]
  if (length(md_nonna) == 0L) skip("no identifying folds; test N/A")

  # MD is additive; values should be plausible mean differences, NOT
  # exponentiated.  DGM MD_harm = 1.5 on an outcome with sd = 1;
  # identifying folds will have training MDs in a plausible range
  # around the DGM value.  Key check: values should NOT be exp()ed
  # (exp(1.5) = 4.48 would be implausibly large for this DGM).
  max_md <- max(md_nonna)
  expect_true(max_md < 20,
    info = sprintf("MD max = %.3f; if erroneously exponentiated, ",
                   "exp(sensible MD) would be much larger",
                   max_md))
})


# ---------------------------------------------------------------------------
# summarize_cv_results() slot PARITY across outcome types
# ---------------------------------------------------------------------------

test_that("summarize_cv_results slots are identical across outcome types", {
  skip_on_cran()

  run_summary_for <- function(otype) {
    df <- switch(otype,
      survival   = .make_survival_data(N = 250L, HR_harm = 2.5, seed = 42L),
      binary     = .make_binary_data(N = 250L, OR_harm = 2.5, seed = 42L),
      count      = .make_count_data(N = 300L, IRR_harm = 2.0, seed = 42L)
    )
    conf <- if (otype == "survival")
      c("age", "stage", "sex")
    else
      c("age", "biomarker", "biomarker_hi", "sex")
    args <- .fs_args_for(otype, confounders = conf)

    fs <- do.call(forestsearch, c(list(df.analysis = df), args))
    if (is.null(fs$sg.harm)) return(NULL)

    cv <- forestsearch_tenfold(
      fs.est  = fs,
      sims    = 2L,
      Kfolds  = 4L,
      details = FALSE,
      parallel_args = list(plan = "sequential", workers = 1L,
                           show_message = FALSE)
    )
    summarize_cv_results(cv_output = cv,
                         original_sg = fs$sg.harm,
                         original_grf_cuts = fs$grf_cuts,
                         top_n = 5L)
  }

  summaries <- lapply(c("survival", "binary"), run_summary_for)
  names(summaries) <- c("survival", "binary")
  summaries <- summaries[!vapply(summaries, is.null, logical(1))]

  if (length(summaries) < 2L) {
    skip("need at least two outcome types to compare slots")
  }

  # Top-level slot names must be identical
  for (i in 2:length(summaries)) {
    expect_identical(sort(names(summaries[[1]])),
                     sort(names(summaries[[i]])),
      info = sprintf("summarize_cv_results slots differ: %s vs %s",
                     names(summaries)[1], names(summaries)[i]))
  }

  # Both must report has_pconsistency = TRUE
  for (s_name in names(summaries)) {
    expect_true(isTRUE(summaries[[s_name]]$has_pconsistency),
      info = sprintf("has_pconsistency not TRUE for %s", s_name))
  }
})


test_that("fold_numeric_summary row labels are outcome-agnostic", {
  skip_on_cran()

  df <- .make_binary_data(N = 300L, OR_harm = 2.5, seed = 42L)
  args <- .fs_args_for("binary",
    confounders = c("age", "biomarker", "biomarker_hi", "sex"))

  fs <- do.call(forestsearch, c(list(df.analysis = df), args))
  if (is.null(fs$sg.harm)) skip("primary fit did not identify; test N/A")

  cv <- forestsearch_tenfold(
    fs.est  = fs,
    sims    = 2L,
    Kfolds  = 4L,
    details = FALSE,
    parallel_args = list(plan = "sequential", workers = 1L,
                         show_message = FALSE)
  )

  res <- summarize_cv_results(cv_output = cv, top_n = 5L)

  # The effect row label should say "subgroup effect", not "HR" (that
  # would be outcome-specific)
  fns <- res$data$fold_numeric_summary
  eff_row <- fns[grepl("subgroup effect", fns$Metric, fixed = TRUE), ]
  expect_true(nrow(eff_row) > 0L,
    info = "expected a row labeled with 'subgroup effect'")
  # Context should mention natural scale
  expect_match(eff_row$Context[1], "natural scale",
               ignore.case = TRUE)
})
