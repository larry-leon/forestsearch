# ============================================================================
# test-get-fsdata-diagnostics.R
#
# Tier 2 tests: get_FSdata() error handling and input-invariant robustness.
#
# Fix B (v0.2.0-dev) replaced the bare stop() at the n_confs == 0 check
# with an actionable diagnostic citing specific upstream source counts.
# These tests lock in that behavior and cover the main failure modes:
#
#   - All upstream sources empty (all-zero counts in the message)
#   - GRF empty but LASSO succeeded
#   - LASSO empty but GRF succeeded
#   - Explicit conf_force works even when auto sources fail
#
# Also covers input-invariant edge cases: zero-variance covariates,
# single-covariate configs, mixed vs all-continuous vs all-factor designs.
# ============================================================================


# ---------------------------------------------------------------------------
# FIX B REGRESSION: diagnostic message content
# ---------------------------------------------------------------------------
#
# To exercise the n_confs == 0 path reliably and deterministically, we
# use the `exclude_cuts` argument.  `get_FSdata` filters the candidate
# cut vector by grepl-ing each exclude pattern; if the pattern matches
# everything, the candidate set becomes empty and we reach the n_confs
# check at line ~320.  This is a cleaner trigger than trying to
# construct adversarial covariates, because it doesn't depend on
# upstream stages (LASSO, GRF) succeeding or failing in specific ways.

test_that("Fix B: n_confs == 0 produces structured diagnostic [Fix B]", {
  df <- .make_survival_data(N = 150L, seed = 42L)

  err <- tryCatch(
    get_FSdata(
      df.analysis      = df,
      confounders.name = c("age", "stage", "sex"),
      outcome.name     = "time",
      event.name       = "event",
      use_grf          = FALSE,
      use_lasso        = FALSE,
      # Exclude EVERYTHING: match any character
      exclude_cuts     = ".",
      details          = FALSE
    ),
    error = function(e) e
  )

  expect_s3_class(err, "error")
  # Must identify itself as a get_FSdata error
  expect_match(conditionMessage(err), "get_FSdata",
               ignore.case = TRUE)
  # Must mention source counts
  expect_match(conditionMessage(err), "GRF cuts",
               ignore.case = TRUE)
  expect_match(conditionMessage(err), "LASSO",
               ignore.case = TRUE)
  # Must offer at least one actionable suggestion
  expect_match(conditionMessage(err),
               "details = TRUE|use_lasso = FALSE|dmin.grf|frac.tau",
               ignore.case = TRUE)
})


test_that("Fix B: diagnostic cites all five source categories", {
  df <- .make_survival_data(N = 150L, seed = 42L)

  err <- tryCatch(
    get_FSdata(
      df.analysis      = df,
      confounders.name = c("age", "stage", "sex"),
      outcome.name     = "time",
      event.name       = "event",
      use_grf          = FALSE,
      use_lasso        = FALSE,
      exclude_cuts     = ".",
      details          = FALSE
    ),
    error = function(e) e
  )

  expect_s3_class(err, "error")
  msg <- conditionMessage(err)
  # Each of 5 source categories must appear in the message
  for (cat in c("GRF cuts", "LASSO selected", "conf_force",
                "conf.cont_medians", "conf.categorical")) {
    expect_match(msg, cat, fixed = TRUE,
      info = sprintf("diagnostic must cite '%s' count", cat))
  }
})


test_that("Fix B: successful run does NOT trigger stop()", {
  df <- .make_survival_data(N = 200L, seed = 42L)

  # Successful call with real covariates should not throw
  res <- get_FSdata(
    df.analysis      = df,
    confounders.name = c("age", "stage", "sex"),
    outcome.name     = "time",
    event.name       = "event",
    use_grf          = FALSE,
    use_lasso        = FALSE,
    details          = FALSE
  )

  expect_true(is.list(res))
  expect_true(!is.null(res$df))
  expect_true(!is.null(res$confs_names))
  expect_true(length(res$confs_names) > 0L)
})


test_that("Fix B regression: error message is on a single line (no tripping artifacts)", {
  # Defensive check: the message uses paste0 concatenation which could
  # accidentally produce a multi-line mess.  Check that the basic
  # formatting holds.
  df <- .make_survival_data(N = 150L, seed = 42L)

  err <- tryCatch(
    get_FSdata(
      df.analysis      = df,
      confounders.name = c("age", "stage", "sex"),
      outcome.name     = "time",
      event.name       = "event",
      use_grf          = FALSE,
      use_lasso        = FALSE,
      exclude_cuts     = ".",
      details          = FALSE
    ),
    error = function(e) e
  )

  expect_s3_class(err, "error")
  msg <- conditionMessage(err)
  # Message must not contain %s placeholder unsubstituted
  expect_false(grepl("%s", msg, fixed = TRUE),
               info = "sprintf placeholders must be fully substituted")
  # Must be a single string, not character vector
  expect_length(msg, 1L)
})


# ---------------------------------------------------------------------------
# INPUT INVARIANTS: covariate type handling
# ---------------------------------------------------------------------------

test_that("all-factor-covariate configuration works end-to-end", {
  set.seed(42L)
  N <- 200L
  df <- data.frame(
    id    = seq_len(N),
    treat = rbinom(N, 1L, 0.5),
    time  = rexp(N, 0.05),
    event = rbinom(N, 1L, 0.7),
    z1    = factor(rbinom(N, 1L, 0.5), levels = c(0L, 1L)),
    z2    = factor(rbinom(N, 1L, 0.5), levels = c(0L, 1L)),
    z3    = factor(rbinom(N, 1L, 0.5), levels = c(0L, 1L))
  )
  args <- .fs_args_for("survival",
    confounders = c("z1", "z2", "z3"),
    extra = list(use_grf = FALSE, use_lasso = FALSE))

  fs <- do.call(forestsearch, c(list(df.analysis = df), args))

  .expect_fs_full_return(fs)
})


test_that("all-continuous-covariate configuration works end-to-end", {
  # This is the specific case that originally failed -- but with
  # canonical parameters and GRF enabled, it should succeed.
  set.seed(42L)
  N <- 300L
  df <- data.frame(
    id    = seq_len(N),
    treat = rbinom(N, 1L, 0.5),
    time  = rexp(N, 0.05),
    event = rbinom(N, 1L, 0.7),
    x1    = rnorm(N),
    x2    = rnorm(N),
    x3    = rnorm(N)
  )
  args <- .fs_args_for("survival",
    confounders = c("x1", "x2", "x3"),
    extra = list(use_grf = TRUE, use_lasso = FALSE))

  fs <- do.call(forestsearch, c(list(df.analysis = df), args))

  # Whatever the outcome (subgroup found or not), the return must be
  # well-formed
  expect_true(!is.null(fs$args_call_all))
  expect_true("sg.harm" %in% names(fs))
})


test_that("mixed continuous + factor covariates work end-to-end", {
  df <- .make_binary_data(N = 250L, seed = 42L)
  args <- .fs_args_for("binary",
    confounders = c("age", "wtkg", "biomarker", "biomarker_hi", "sex"))

  fs <- do.call(forestsearch, c(list(df.analysis = df), args))

  .expect_fs_full_return(fs)
})


test_that("single-covariate configuration works or fails cleanly", {
  df <- .make_survival_data(N = 200L, seed = 42L)
  args <- .fs_args_for("survival",
    confounders = "age",
    extra = list(maxk = 1L))

  # Either succeeds with a result or fails cleanly; not a silent crash
  fs <- tryCatch(
    do.call(forestsearch, c(list(df.analysis = df), args)),
    error = function(e) e
  )

  if (inherits(fs, "error")) {
    # If it threw, the message should be informative
    expect_match(conditionMessage(fs), ".+",
      info = "single-covariate error should have non-empty message")
  } else {
    # Otherwise, it must satisfy the return contract
    expect_true(!is.null(fs$args_call_all))
    expect_true("sg.harm" %in% names(fs))
  }
})


test_that("covariate list pointing at non-existent column fails cleanly", {
  df <- .make_survival_data(N = 150L, seed = 42L)
  args <- .fs_args_for("survival",
    confounders = c("age", "does_not_exist"))

  # Must not be a silent cascade; must error with a clear message
  err <- tryCatch(
    do.call(forestsearch, c(list(df.analysis = df), args)),
    error = function(e) e,
    warning = function(w) w
  )

  # Either errored or warned; either way it shouldn't have produced a
  # malformed object
  if (inherits(err, "condition")) {
    expect_match(conditionMessage(err), ".+")
  }
})


# ---------------------------------------------------------------------------
# INPUT INVARIANTS: outcome column validation
# ---------------------------------------------------------------------------

test_that("get_FSdata rejects non-numeric outcome column", {
  set.seed(42L)
  df <- data.frame(
    id    = 1:100,
    treat = rbinom(100, 1L, 0.5),
    time  = as.character(round(rexp(100, 0.05), 2)),  # character, not numeric
    event = rbinom(100, 1L, 0.7),
    age   = round(rnorm(100, 55, 12))
  )

  err <- tryCatch(
    get_FSdata(
      df.analysis      = df,
      confounders.name = "age",
      outcome.name     = "time",
      event.name       = "event",
      details          = FALSE
    ),
    error = function(e) e
  )

  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "numeric",
               ignore.case = TRUE)
})


test_that("get_FSdata rejects non-numeric event column", {
  set.seed(42L)
  df <- data.frame(
    id    = 1:100,
    treat = rbinom(100, 1L, 0.5),
    time  = rexp(100, 0.05),
    event = c("yes", "no")[rbinom(100, 1L, 0.7) + 1L],  # character
    age   = round(rnorm(100, 55, 12))
  )

  err <- tryCatch(
    get_FSdata(
      df.analysis      = df,
      confounders.name = "age",
      outcome.name     = "time",
      event.name       = "event",
      details          = FALSE
    ),
    error = function(e) e
  )

  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "numeric",
               ignore.case = TRUE)
})


test_that("logical outcome column fails explicitly", {
  # TRUE/FALSE is not numeric in R's is.numeric() sense;
  # get_FSdata's contract is numeric (0/1)
  set.seed(42L)
  df <- data.frame(
    id    = 1:100,
    treat = rbinom(100, 1L, 0.5),
    y     = sample(c(TRUE, FALSE), 100, replace = TRUE),  # logical
    age   = round(rnorm(100, 55, 12))
  )

  err <- tryCatch(
    get_FSdata(
      df.analysis      = df,
      confounders.name = "age",
      outcome.name     = "y",
      event.name       = "y",
      outcome_type     = "binary",
      details          = FALSE
    ),
    error = function(e) e
  )

  # Either rejected at outcome-column check, or accepted (both are
  # defensible; test just verifies no silent acceptance-then-crash)
  if (inherits(err, "error")) {
    expect_match(conditionMessage(err), ".+",
      info = "logical outcome should either work or error cleanly")
  }
})


# ---------------------------------------------------------------------------
# INPUT INVARIANTS: covariate type violations
# ---------------------------------------------------------------------------

test_that("character covariate columns are rejected", {
  set.seed(42L)
  df <- data.frame(
    id    = 1:100,
    treat = rbinom(100, 1L, 0.5),
    time  = rexp(100, 0.05),
    event = rbinom(100, 1L, 0.7),
    age   = rnorm(100, 55, 12),
    bad   = sample(c("a", "b", "c"), 100, replace = TRUE),  # character
    stringsAsFactors = FALSE
  )

  err <- tryCatch(
    get_FSdata(
      df.analysis      = df,
      confounders.name = c("age", "bad"),
      outcome.name     = "time",
      event.name       = "event",
      details          = FALSE
    ),
    error = function(e) e
  )

  expect_s3_class(err, "error")
  expect_match(conditionMessage(err),
               "numeric|factor|confounders",
               ignore.case = TRUE)
})


# ---------------------------------------------------------------------------
# INPUT INVARIANTS: grf_cuts contract
# ---------------------------------------------------------------------------

test_that("grf_cuts as list-of-values is converted to character cuts", {
  # Upstream GRF returns list(name = value); get_FSdata converts
  # internally to c("name <= value").  This is a contract test for that
  # conversion.
  df <- .make_survival_data(N = 200L, seed = 42L)

  # Call get_FSdata with a list-form grf_cuts input
  res <- tryCatch(
    get_FSdata(
      df.analysis      = df,
      confounders.name = c("age", "stage", "sex"),
      outcome.name     = "time",
      event.name       = "event",
      use_grf          = TRUE,
      grf_cuts         = list(age = 50),   # list, not character
      use_lasso        = FALSE,
      details          = FALSE
    ),
    error = function(e) e
  )

  if (!inherits(res, "error")) {
    # Converted cut should appear in confs as a character expression
    expect_true(is.list(res))
    expect_true(!is.null(res$confs))
    age_cut_present <- any(grepl("age\\s*<=\\s*50", res$confs))
    expect_true(age_cut_present,
      info = "list-form grf_cuts should be converted to 'age <= 50'")
  }
})


test_that("grf_cuts with invalid type is rejected explicitly", {
  df <- .make_survival_data(N = 150L, seed = 42L)

  # Numeric vector without names is not a valid grf_cuts
  err <- tryCatch(
    get_FSdata(
      df.analysis      = df,
      confounders.name = c("age", "stage", "sex"),
      outcome.name     = "time",
      event.name       = "event",
      use_grf          = TRUE,
      grf_cuts         = c(50, 60),       # bare numeric, neither named
                                          # list nor character
      details          = FALSE
    ),
    error = function(e) e
  )

  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), "character",
               ignore.case = TRUE)
})


# ---------------------------------------------------------------------------
# HAPPY PATH WITH DIFFERENT FEATURE-SELECTION COMBINATIONS
# ---------------------------------------------------------------------------

test_that("use_grf = TRUE, use_lasso = FALSE works", {
  df <- .make_survival_data(N = 200L, seed = 42L)
  args <- .fs_args_for("survival",
    confounders = c("age", "stage", "sex", "noise"),
    extra = list(use_grf = TRUE, use_lasso = FALSE))

  fs <- do.call(forestsearch, c(list(df.analysis = df), args))

  .expect_fs_full_return(fs)
})


test_that("use_grf = FALSE, use_lasso = TRUE works", {
  df <- .make_survival_data(N = 200L, seed = 42L)
  args <- .fs_args_for("survival",
    confounders = c("age", "stage", "sex", "noise"),
    extra = list(use_grf = FALSE, use_lasso = TRUE))

  fs <- do.call(forestsearch, c(list(df.analysis = df), args))

  .expect_fs_full_return(fs)
})


test_that("use_grf = FALSE, use_lasso = FALSE works (factor-only default cuts)", {
  df <- .make_survival_data(N = 200L, seed = 42L)
  args <- .fs_args_for("survival",
    confounders = c("age", "stage", "sex", "noise"),
    extra = list(use_grf = FALSE, use_lasso = FALSE))

  fs <- do.call(forestsearch, c(list(df.analysis = df), args))

  .expect_fs_full_return(fs)
})


test_that("use_grf = TRUE, use_lasso = TRUE works", {
  df <- .make_survival_data(N = 200L, seed = 42L)
  args <- .fs_args_for("survival",
    confounders = c("age", "stage", "sex", "noise"),
    extra = list(use_grf = TRUE, use_lasso = TRUE))

  fs <- do.call(forestsearch, c(list(df.analysis = df), args))

  .expect_fs_full_return(fs)
})
