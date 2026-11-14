# =============================================================================
# Core test of forestsearch(); runs on CRAN (no skip_on_cran()).
#
# One complete default survival forestsearch() call on survival::gbsg
# (hormonal therapy as the treatment), on a sequential plan.  It checks the
# structure documented in ?forestsearch (Value) and that a second call with
# the same seed returns the same result.  Structure checks and
# same-session comparisons only: no numeric value is typed into this file.
# The comparison uses a tolerance because some linear-algebra libraries (MKL) do not reproduce floating-point results bit for bit, even within a session.
# The default call fits no grf forest and starts no parallel workers.
# =============================================================================

.core_run <- function() {
  suppressMessages(forestsearch(
    df.analysis      = survival::gbsg,
    outcome.name     = "rfstime",
    event.name       = "status",
    treat.name       = "hormon",
    id.name          = "pid",
    confounders.name = c("age", "meno", "size", "grade", "nodes", "pgr", "er"),
    parallel_args    = list(plan = "sequential")
  ))
}

# The first fit, made once and shared by both blocks.
.core_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) fit <<- .core_run()
    fit
  }
})

# The run-time fields: the only parts of the object that differ between calls.
.core_strip_timing <- function(fs) {
  fs$minutes_all <- NULL
  fs$find.grps$time_search <- NULL
  fs
}

test_that("a default survival forestsearch() call returns the documented structure", {
  fs <- .core_fit()
  dat <- survival::gbsg

  expect_s3_class(fs, "forestsearch")
  expect_type(fs, "list")
  documented <- c("grp.consistency", "find.grps", "confounders.candidate",
                  "confounders.evaluated", "df.est", "df.predict", "df.test",
                  "minutes_all", "grf_res", "sg_focus", "sg.harm", "grf_cuts",
                  "prop_maxk", "max_sg_est", "grf_plot", "args_call_all",
                  "threshold_config", "family_status", "mr_inference",
                  "mr_harm_confirmed")
  expect_identical(setdiff(documented, names(fs)), character())

  # The call as recorded: every formal, the data, the default seed and focus
  expect_identical(setdiff(names(formals(forestsearch)), names(fs$args_call_all)),
                   character())
  expect_identical(fs$args_call_all$df.analysis, dat)
  expect_identical(fs$args_call_all$seedit, eval(formals(forestsearch)$seedit))
  expect_identical(fs$sg_focus, eval(formals(forestsearch)$sg_focus))

  # Thresholds: the survival path compares hazard ratios; the stored
  # screening and consistency values are their logs
  tc <- fs$threshold_config
  expect_identical(setdiff(c("outcome_type", "effect_measure", "screening",
                             "consistency", "screening_natural",
                             "consistency_natural", "scale", "pconsistency"),
                           names(tc)), character())
  expect_identical(tc$outcome_type, "survival")
  expect_identical(tc$screening_natural, fs$args_call_all$hr.threshold)
  expect_identical(tc$consistency_natural, fs$args_call_all$hr.consistency)
  expect_equal(tc$screening, log(tc$screening_natural))
  expect_equal(tc$consistency, log(tc$consistency_natural))
  expect_identical(tc$pconsistency, fs$args_call_all$pconsistency.threshold)

  # The search: candidate factors are columns of the analysis data, and each
  # label of the selected subgroup is one of the evaluated factors
  expect_true(is.list(fs$find.grps))
  expect_type(fs$confounders.candidate, "character")
  expect_type(fs$confounders.evaluated, "character")
  expect_true(all(fs$confounders.candidate %in% names(fs$df.est)))
  expect_type(fs$sg.harm, "character")
  expect_true(all(gsub("^!?\\{|\\}$", "", fs$sg.harm) %in% fs$confounders.evaluated))
  expect_true(is.numeric(fs$max_sg_est))
  expect_true(is.numeric(fs$prop_maxk))
  expect_true(is.numeric(fs$minutes_all))

  # The consistency evaluation
  gc <- fs$grp.consistency
  expect_identical(setdiff(c("out_sg", "sg_focus", "df_flag", "sg.harm", "sg.harm.id",
                             "algorithm", "n_candidates_evaluated", "n_passed"),
                           names(gc)), character())
  expect_identical(gc$sg_focus, fs$sg_focus)
  expect_identical(gc$sg.harm, fs$sg.harm)
  expect_true(gc$algorithm %in% c("twostage", "fixed"))
  expect_true(gc$n_passed <= gc$n_candidates_evaluated)

  # The returned data: one row per patient in input order; the subgroup and
  # its complement are both non-empty, and treat.recommend is the complement
  # of subgroup membership
  expect_s3_class(fs$df.est, "data.frame")
  expect_identical(fs$df.est$pid, dat$pid)
  member <- as.logical(gc$sg.harm.id)
  expect_true(any(member))
  expect_false(all(member))
  expect_true(is.numeric(fs$df.est$treat.recommend))
  expect_identical(as.logical(fs$df.est$treat.recommend), !member)
  expect_identical(gc$df_flag$id, fs$df.est$pid)
  expect_identical(gc$df_flag$treat.recommend, fs$df.est$treat.recommend)
  expect_null(fs$df.test)

  # Defaults off: no front end (no LASSO, grf or DINA), no multiplier resampling
  expect_identical(fs$family_status, "no-front-end")
  expect_null(fs$grf_res)
  expect_null(fs$mr_inference)
  expect_identical(fs$mr_harm_confirmed, NA)
})

test_that("a second call with the same seed returns the same result", {
  a <- .core_fit()
  b <- .core_run()
  # serialize() round trips drop data.table's session-only self-reference
  # pointers, so the comparison sees values only
  roundtrip <- function(x) unserialize(serialize(x, NULL))
  expect_equal(roundtrip(.core_strip_timing(b)), roundtrip(.core_strip_timing(a)), tolerance = 1e-8)
})
