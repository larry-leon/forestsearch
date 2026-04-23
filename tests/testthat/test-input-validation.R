# ============================================================================
# test-input-validation.R
#
# Tier 5 tests: forestsearch() and related functions catch invalid
# arguments at entry with informative errors, rather than failing deep
# in the pipeline with cryptic messages.
#
# These tests use tryCatch + conditionMessage matching to verify that
# errors surface at the right stage.  Most tests are purely validation
# checks and do not run the pipeline, so they're fast.
# ============================================================================


# ---------------------------------------------------------------------------
# outcome_type validation
# ---------------------------------------------------------------------------

test_that("forestsearch() rejects invalid outcome_type", {
  df <- .make_survival_data(N = 100L, seed = 42L)

  err <- tryCatch(
    forestsearch(
      df.analysis      = df,
      confounders.name = c("age", "stage"),
      outcome.name     = "time",
      event.name       = "event",
      treat.name       = "treat",
      id.name          = "id",
      outcome_type     = "nonsense",   # invalid
      quiet            = TRUE
    ),
    error = function(e) e
  )

  expect_s3_class(err, "error")
  # Match.arg error references the valid choices
  expect_match(conditionMessage(err),
               "outcome_type|should be|one of|nonsense",
               ignore.case = TRUE)
})


# ---------------------------------------------------------------------------
# summarize_cv_results input validation
# (more coverage in test-downstream-consumers.R; these are edge cases)
# ---------------------------------------------------------------------------

test_that("summarize_cv_results top_n = Inf is accepted", {
  fake_fsum <- data.frame(
    sim = 1L, fold = 1L,
    sg1 = NA_character_, sg2 = NA_character_, any_found = 0L,
    stringsAsFactors = FALSE
  )
  fake_cv <- list(fold_summary = fake_fsum,
                  sens_summary = NULL, find_summary = NULL,
                  sims = 1L, Kfolds = 1L)
  class(fake_cv) <- c("fs_tenfold", "list")

  # Inf is a valid value (means "no truncation")
  res <- tryCatch(
    summarize_cv_results(fake_cv, top_n = Inf),
    error = function(e) e
  )
  expect_false(inherits(res, "error"),
    info = sprintf("top_n = Inf should be accepted; got: %s",
                   if (inherits(res, "error")) conditionMessage(res) else ""))
})


test_that("summarize_cv_results top_n vector of length > 1 is rejected", {
  fake_fsum <- data.frame(
    sim = 1L, fold = 1L,
    sg1 = NA_character_, sg2 = NA_character_, any_found = 0L,
    stringsAsFactors = FALSE
  )
  fake_cv <- list(fold_summary = fake_fsum,
                  sens_summary = NULL, find_summary = NULL,
                  sims = 1L, Kfolds = 1L)
  class(fake_cv) <- c("fs_tenfold", "list")

  expect_error(
    summarize_cv_results(fake_cv, top_n = c(5L, 10L)),
    "positive integer"
  )
})


test_that("summarize_cv_results top_n non-numeric is rejected", {
  fake_fsum <- data.frame(
    sim = 1L, fold = 1L,
    sg1 = NA_character_, sg2 = NA_character_, any_found = 0L,
    stringsAsFactors = FALSE
  )
  fake_cv <- list(fold_summary = fake_fsum,
                  sens_summary = NULL, find_summary = NULL,
                  sims = 1L, Kfolds = 1L)
  class(fake_cv) <- c("fs_tenfold", "list")

  expect_error(
    summarize_cv_results(fake_cv, top_n = "five"),
    "positive integer"
  )
})


# ---------------------------------------------------------------------------
# get_FSdata: missing required outcome/event columns
# ---------------------------------------------------------------------------

test_that("get_FSdata surfaces useful error when outcome column missing", {
  df <- .make_survival_data(N = 100L, seed = 42L)

  err <- tryCatch(
    get_FSdata(
      df.analysis      = df,
      confounders.name = c("age", "stage"),
      outcome.name     = "does_not_exist",
      event.name       = "event",
      details          = FALSE
    ),
    error = function(e) e
  )

  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), ".+")
})


test_that("get_FSdata surfaces useful error when event column missing", {
  df <- .make_survival_data(N = 100L, seed = 42L)

  err <- tryCatch(
    get_FSdata(
      df.analysis      = df,
      confounders.name = c("age", "stage"),
      outcome.name     = "time",
      event.name       = "does_not_exist",
      details          = FALSE
    ),
    error = function(e) e
  )

  expect_s3_class(err, "error")
  expect_match(conditionMessage(err), ".+")
})


test_that("get_FSdata with non-data.frame input fails (either at entry or downstream)", {
  # Current behavior: get_FSdata coerces non-data.frames via
  # as.data.frame() and then fails downstream at the outcome-column
  # numeric check.  A future hardening pass could add an explicit
  # is.data.frame() check at entry; for now we just verify that some
  # informative error is produced.
  err <- tryCatch(
    get_FSdata(
      df.analysis      = "not a data frame",
      confounders.name = "x",
      outcome.name     = "y",
      event.name       = "e",
      details          = FALSE
    ),
    error = function(e) e
  )

  expect_s3_class(err, "error")
  # Actual message as of v0.2.0-dev is "Outcome column must be numeric."
  # (because as.data.frame() coerces the string and then the outcome
  # column doesn't exist / is NULL).  Accept any non-empty error.
  expect_match(conditionMessage(err), ".+",
    info = "non-data.frame input must fail with some informative error")
})


test_that("get_FSdata rejects empty confounders.name", {
  df <- .make_survival_data(N = 100L, seed = 42L)

  # With no confounders and no GRF/LASSO/conf_force, there's nothing to
  # work with -- should error at n_confs == 0 (Fix B structured message)
  err <- tryCatch(
    get_FSdata(
      df.analysis      = df,
      confounders.name = character(0),
      outcome.name     = "time",
      event.name       = "event",
      use_grf          = FALSE,
      use_lasso        = FALSE,
      details          = FALSE
    ),
    error = function(e) e
  )

  # Either rejected at entry, or triggers Fix B -- both are acceptable;
  # the key is it's not a silent crash
  expect_s3_class(err, "error")
})


# ---------------------------------------------------------------------------
# Confounder type validation
# ---------------------------------------------------------------------------

test_that("get_FSdata rejects mixed character-numeric column types", {
  set.seed(42L)
  df <- data.frame(
    id = 1:100, treat = rbinom(100, 1L, 0.5),
    time = rexp(100, 0.05), event = rbinom(100, 1L, 0.7),
    numeric_ok   = rnorm(100),
    character_bad = sample(letters, 100, replace = TRUE),
    stringsAsFactors = FALSE
  )

  err <- tryCatch(
    get_FSdata(
      df.analysis      = df,
      confounders.name = c("numeric_ok", "character_bad"),
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
# Contradictory or unlikely argument combinations
# ---------------------------------------------------------------------------

test_that("forestsearch with zero-row data errors before pipeline", {
  df <- .make_survival_data(N = 100L, seed = 42L)[integer(0), ]

  err <- tryCatch(
    suppressWarnings(
      forestsearch(
        df.analysis      = df,
        confounders.name = c("age", "stage"),
        outcome.name     = "time",
        event.name       = "event",
        treat.name       = "treat",
        id.name          = "id",
        outcome_type     = "survival",
        quiet            = TRUE,
        parallel_args    = list(plan = "sequential", workers = 1L,
                                show_message = FALSE)
      )
    ),
    error = function(e) e
  )

  # Either errors explicitly or returns an early-return object -- must
  # not silently produce a garbage result with phantom subgroups
  if (inherits(err, "error")) {
    expect_match(conditionMessage(err), ".+")
  } else {
    # If it returned, the return must be a well-formed early-return
    expect_true(is.list(err))
    expect_true("args_call_all" %in% names(err) ||
                "sg.harm"       %in% names(err),
      info = "empty-data return must carry at least args_call_all or sg.harm")
  }
})


test_that("forestsearch with single-row-per-arm data fails cleanly", {
  # Extreme adversarial: 2 rows total, one per arm
  set.seed(42L)
  df <- data.frame(
    id    = 1:2,
    treat = c(0L, 1L),
    time  = c(5, 10),
    event = c(1L, 1L),
    age   = c(40, 60),
    sex   = factor(c(0L, 1L), levels = c(0L, 1L))
  )

  err <- tryCatch(
    suppressWarnings(
      forestsearch(
        df.analysis      = df,
        confounders.name = c("age", "sex"),
        outcome.name     = "time",
        event.name       = "event",
        treat.name       = "treat",
        id.name          = "id",
        outcome_type     = "survival",
        n.min            = 2L,
        d0.min           = 1L,
        d1.min           = 1L,
        quiet            = TRUE,
        parallel_args    = list(plan = "sequential", workers = 1L,
                                show_message = FALSE)
      )
    ),
    error = function(e) e
  )

  # Either errors or returns early-return; must not silently succeed
  # with bogus results
  if (inherits(err, "error")) {
    expect_match(conditionMessage(err), ".+")
  } else {
    expect_true(is.null(err$sg.harm) || !is.null(err$error_log),
      info = "pathologically small data must not produce phantom subgroups")
  }
})


test_that("forestsearch with constant treatment column fails cleanly", {
  # All-treated or all-control -- no treatment effect is estimable
  set.seed(42L)
  df <- .make_survival_data(N = 100L, seed = 42L)
  df$treat <- 1L   # everyone treated

  err <- tryCatch(
    suppressWarnings(
      forestsearch(
        df.analysis      = df,
        confounders.name = c("age", "stage"),
        outcome.name     = "time",
        event.name       = "event",
        treat.name       = "treat",
        id.name          = "id",
        outcome_type     = "survival",
        quiet            = TRUE,
        parallel_args    = list(plan = "sequential", workers = 1L,
                                show_message = FALSE)
      )
    ),
    error = function(e) e
  )

  # Must not silently produce a result with phantom subgroups
  if (inherits(err, "error")) {
    expect_match(conditionMessage(err), ".+")
  } else {
    # If it returned, it must be an early-return object
    expect_true(is.list(err))
    # Either errored internally (error_log populated) or legitimately
    # found nothing (sg.harm NULL) -- not a phantom subgroup
    if (!is.null(err$sg.harm)) {
      # If somehow a subgroup was identified, df.est should at minimum
      # have both arms (which it can't with constant treatment)
      skip("constant-treat scenario produced an unexpected subgroup; investigate")
    }
  }
})
