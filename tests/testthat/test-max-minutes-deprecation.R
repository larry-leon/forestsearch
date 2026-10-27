# =============================================================================
# max.minutes deprecation (0.4.0.9000): the formal is kept in forestsearch()
# and subgroup.search(), never consulted, and warned about.  forestsearch()
# warns on a non-default value (args_call_all always carries max.minutes and
# is replayed by bootstrap/CV, so a bare !missing() trigger would fire there);
# subgroup.search() warns whenever it is passed.  forestsearch() does not
# forward it to subgroup.search().
# =============================================================================

.mm_df <- .make_continuous_data(N = 300L, MD_harm = 1.5, seed = 42L)

.mm_run <- function(...) {
  suppressMessages(forestsearch(
    df.analysis = .mm_df, confounders.name = c("age", "biomarker", "biomarker_hi", "sex"),
    outcome.name = "y", event.name = "y", treat.name = "treat", id.name = "id",
    outcome_type = "continuous", effect_measure = "MD", adverse_outcome = TRUE,
    effect.threshold = 0.5, consistency.threshold = 0.0, pconsistency.threshold = 0.5,
    fs.splits = 20, n.min = 30, maxk = 2, use_lasso = FALSE, use_grf = FALSE,
    sg_focus = "eff", is.RCT = TRUE, details = FALSE,
    quiet = TRUE, seedit = 1L, parallel_args = list(plan = "sequential"), ...))
}

# Collect every warning message raised while evaluating expr; return both.
.mm_capture <- function(expr) {
  w <- character(0)
  value <- withCallingHandlers(expr, warning = function(cond) {
    w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning")
  })
  list(value = value, warnings = w)
}

# Timing fields, as stripped by the existing identity tests
# (test-declaration-c0.R, test-declaration-calibration.R).
.mm_strip_t <- function(x) {
  if (is.list(x) && !is.data.frame(x)) {
    drop <- names(x) %in% c("timing_seconds", "minutes_all", "time_search")
    if (any(drop)) x <- x[!drop]
    for (i in seq_along(x)) if (!is.null(x[[i]])) x[[i]] <- .mm_strip_t(x[[i]])
  }
  x
}

# Saved-vs-saved: compare readRDS() copies (data.table internals).
.mm_rds <- function(x) {
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f), add = TRUE)
  saveRDS(x, f)
  readRDS(f)
}

.mm_default <- .mm_capture(.mm_run())
.mm_passed  <- .mm_capture(.mm_run(max.minutes = 5))

test_that("a default forestsearch() call emits no max.minutes warning", {
  expect_false(any(grepl("max.minutes", .mm_default$warnings, fixed = TRUE)))
  expect_false(any(grepl("subgroup.search():", .mm_default$warnings, fixed = TRUE)))
})

test_that("max.minutes = 5 emits exactly one deprecation warning", {
  hits <- grep("'max.minutes' is deprecated", .mm_passed$warnings, fixed = TRUE)
  expect_identical(length(hits), 1L)
  expect_true(startsWith(.mm_passed$warnings[hits], "forestsearch(): "))
  # not forwarded: subgroup.search() never sees it
  expect_false(any(grepl("subgroup.search():", .mm_passed$warnings, fixed = TRUE)))
})

test_that("max.minutes is inert: the fit is identical apart from the recorded argument", {
  a <- .mm_default$value
  b <- .mm_passed$value
  expect_identical(b$args_call_all$max.minutes, 5)
  a$args_call_all$max.minutes <- NULL
  b$args_call_all$max.minutes <- NULL
  expect_identical(.mm_rds(.mm_strip_t(b)), .mm_rds(.mm_strip_t(a)))
})

test_that("subgroup.search() warns once with its own prefix when passed max.minutes", {
  d <- .make_survival_data(N = 200L, seed = 42L)
  Z <- cbind(age_le50 = as.integer(d$age <= 50),
             stage3   = as.integer(d$stage == "3"),
             sex1     = as.integer(d$sex == "1"))
  call_ss <- function(...) {
    subgroup.search(Y = d$time, Event = d$event, Treat = d$treat, Z = Z,
                    n.min = 20, d0.min = 5, d1.min = 5, hr.threshold = 1.0,
                    maxk = 2, parallel_workers = 1, ...)
  }
  passed <- .mm_capture(suppressMessages(call_ss(max.minutes = 5)))
  hits <- grep("'max.minutes' is deprecated", passed$warnings, fixed = TRUE)
  expect_identical(length(hits), 1L)
  expect_true(startsWith(passed$warnings[hits], "subgroup.search(): "))
  dflt <- .mm_capture(suppressMessages(call_ss()))
  expect_false(any(grepl("max.minutes", dflt$warnings, fixed = TRUE)))
})

test_that("the max.minutes formals are kept through this release", {
  expect_true("max.minutes" %in% names(formals(forestsearch)))
  expect_true("max.minutes" %in% names(formals(subgroup.search)))
})
