# A-151: test_hte_from_forestsearch() reads fs_result$confounders.evaluated.
# On the consistency path those are subgroup-factor labels: a categorical
# confounder keeps its bare name ("bm"), a continuous one becomes cut labels
# such as "age <= 55" (get_fsdata(), default cut_type).  Each label is mapped
# to the df columns it refers to (a label that is exactly a column name is
# that column; otherwise all.vars() of the parsed label, intersected with
# names(df)), and the unique columns are tested.  A label that cannot be
# resolved is skipped with a warning naming it, and a call left with no
# covariate warns before returning NULL.

skip_on_cran()

.a151_data <- function() {
  set.seed(11); n <- 300
  df <- data.frame(id = seq_len(n), treat = rbinom(n, 1, 0.5),
                   age = round(rnorm(n, 55, 12)),
                   bmi = round(rnorm(n, 27, 4), 1),
                   bm  = rbinom(n, 1, 0.5))
  harm <- df$age <= 50 & df$bm == 1
  tt <- rexp(n, 0.05 * exp(log(2.5) * df$treat * harm))
  df$time <- pmin(tt, 60); df$event <- as.integer(tt <= 60)
  df
}

.a151_fit <- function(df, confs) {
  old_plan <- future::plan(future::sequential)
  on.exit(future::plan(old_plan), add = TRUE)
  suppressWarnings(suppressMessages(forestsearch(
    df.analysis = df, outcome.name = "time", event.name = "event",
    treat.name = "treat", id.name = "id", confounders.name = confs,
    use_lasso = FALSE, use_grf = FALSE,
    fs.splits = 50, n.min = 30, d0.min = 5, d1.min = 5, hr.threshold = 1.1,
    pconsistency.threshold = 0.5, details = FALSE, quiet = TRUE,
    plot.sg = FALSE,
    parallel_args = list(plan = "sequential", workers = 1L,
                         show_message = FALSE))))
}

# Evaluate expr, muffling and recording every warning message.
.a151_collect <- function(expr) {
  w <- character(0)
  value <- withCallingHandlers(expr, warning = function(cnd) {
    w <<- c(w, conditionMessage(cnd))
    invokeRestart("muffleWarning")
  })
  list(value = value, warnings = w)
}

test_that("continuous and categorical confounders are all used (no LASSO, no GRF)", {
  skip_on_cran()
  df <- .a151_data()
  fit <- .a151_fit(df, c("age", "bmi", "bm"))
  ev <- fit$confounders.evaluated

  # Preconditions (the case A-151 describes): age and bmi have 4 or more
  # distinct values, so under the default cut_type they appear only as cut
  # labels "var <= value" (process_conf_force_expr), never as bare names;
  # bm is binary, so it keeps its bare name.
  expect_true(any(startsWith(ev, "age <= ")))
  expect_true(any(startsWith(ev, "bmi <= ")))
  expect_false(any(c("age", "bmi") %in% ev))
  expect_true("bm" %in% ev)

  r <- .a151_collect(test_hte_from_forestsearch(fit, df, "time", "treat"))

  # Every label refers to a column of df (age, bmi or bm), so no warning.
  expect_length(r$warnings, 0)
  # Expected covariates: the three names passed to forestsearch(), by hand.
  # Order: get_fsdata() returns confs = c(conf.cont_cuts, conf.categorical),
  # with the cuts in confounders.name order, so the first label of each
  # variable appears in the order age, bmi, bm -- the order passed above
  # (continuous confounders listed first on purpose).
  expect_identical(r$value$covariates_used, c("age", "bmi", "bm"))
  # Three covariates -> ncol(X) = 3 (integer, from ncol()).
  expect_identical(r$value$n_covariates, 3L)
  # The source text is unchanged by the fix (spec): the start-SHA string.
  expect_identical(r$value$source,
                   "confounders.evaluated (GRF/LASSO selected)")
  # The returned test equals test_hte() called directly on those columns,
  # built as its documentation asks: a numeric matrix, one column each.
  X <- cbind(age = as.numeric(df$age), bmi = as.numeric(df$bmi),
             bm = as.numeric(df$bm))
  expect_equal(r$value$hte_test, test_hte(df$time, df$treat, X))
})

test_that("a fit with only continuous confounders returns a test", {
  skip_on_cran()
  df <- .a151_data()
  fit <- .a151_fit(df, c("age", "bmi"))
  ev <- fit$confounders.evaluated

  # Precondition: every evaluated label is a cut label, none is a column
  # name, so a column-name lookup would find nothing.
  expect_true(length(ev) > 0)
  expect_false(any(ev %in% names(df)))

  r <- .a151_collect(test_hte_from_forestsearch(fit, df, "time", "treat"))

  # By hand: every label in ev is a cut label "age <= v" or "bmi <= v"
  # (precondition above); all.vars() of the parsed label is "age" or "bmi",
  # both columns of df, so no label is skipped (no "cannot be resolved"
  # warning) and X gets two columns (no "No covariate found" warning, no
  # NULL).  hte_test is test_hte()'s return value, whose class test_hte()
  # sets to "hte_test".
  # Before the A-151 fix the labels were looked up as column names, none
  # matched, and the call returned NULL without a warning: the warning
  # count passed then too; the NULL and class checks did not.
  expect_length(r$warnings, 0)
  expect_false(is.null(r$value))
  expect_s3_class(r$value$hte_test, "hte_test")
  # By hand: the two confounders passed, in the order passed (cuts are
  # emitted in confounders.name order).
  expect_identical(r$value$covariates_used, c("age", "bmi"))
  expect_identical(r$value$n_covariates, 2L)
  X <- cbind(age = as.numeric(df$age), bmi = as.numeric(df$bmi))
  expect_equal(r$value$hte_test, test_hte(df$time, df$treat, X))
})

test_that("a label that cannot be resolved is skipped with a warning naming it", {
  df <- .a151_data()
  # A minimal result object: test_hte_from_forestsearch() reads only
  # $confounders.evaluated, $confounders.candidate and $find.grps.
  fs_like <- list(confounders.evaluated =
                    c("age <= 55", "zzz <= 1", "bmi <=", "bm"))
  r <- .a151_collect(test_hte_from_forestsearch(fs_like, df, "time", "treat"))

  # By hand: "age <= 55" -> age; "zzz <= 1" parses but zzz is not a column
  # of df; "bmi <=" does not parse; "bm" -> bm.  So two labels cannot be
  # resolved (one warning each, naming the label) and age, bm are used.
  expect_length(r$warnings, 2)
  expect_true(any(grepl("zzz <= 1", r$warnings, fixed = TRUE)))
  expect_true(any(grepl("bmi <=", r$warnings, fixed = TRUE)))
  expect_identical(r$value$covariates_used, c("age", "bm"))
  X <- cbind(age = as.numeric(df$age), bm = as.numeric(df$bm))
  expect_equal(r$value$hte_test, test_hte(df$time, df$treat, X))
})

test_that("no covariate left: a warning, then NULL", {
  df <- .a151_data()

  # One unresolvable label: by the spec, one warning naming the label and
  # one warning that no covariate remains, then NULL.  Code path: the
  # confounders.evaluated loop warns once for "zzz <= 1" (zzz is not a
  # column), leaving no covariate; the empty-X_cols branch then warns once
  # ("No covariate found ...", the message this fix defines).  2 in all.
  r <- .a151_collect(test_hte_from_forestsearch(
    list(confounders.evaluated = "zzz <= 1"), df, "time", "treat"))
  expect_null(r$value)
  expect_length(r$warnings, 2)
  expect_true(any(grepl("zzz <= 1", r$warnings, fixed = TRUE)))
  expect_true(any(grepl("No covariate found", r$warnings, fixed = TRUE)))

  # Candidate names that are not columns of df (the q1, q2 dummy names of
  # the consistency path).  Code path: the confounders.candidate branch
  # emits no per-label warning (it only sets covs), the column loop skips
  # q1 and q2 (not columns), and the empty-X_cols branch warns once.  So
  # exactly 1 warning, the no-covariate one, then NULL (not a silent NULL).
  r2 <- .a151_collect(test_hte_from_forestsearch(
    list(confounders.candidate = c("q1", "q2")), df, "time", "treat"))
  expect_null(r2$value)
  expect_length(r2$warnings, 1)
  expect_true(any(grepl("No covariate found", r2$warnings, fixed = TRUE)))
})

test_that("a label that is exactly a column name is that column (non-syntactic)", {
  df <- .a151_data()
  # forestsearch() accepts categorical confounders whose names are not
  # syntactic R names (evaluate_comparison() looks a bare name up as an
  # exact string), and stores them as bare labels in confounders.evaluated.
  # Add two such columns, plus columns ER and pos, which the parsed
  # expression ER - pos would refer to.  Values: deterministic functions of
  # the existing columns (no extra RNG draws).
  df[["bm status"]] <- 1L - df$bm
  df[["ER-pos"]] <- as.integer(df$bmi > 27)
  df$ER <- df$age / 10
  df$pos <- df$bmi / 10

  # By hand: each label is exactly a column name of df, so it refers to that
  # column -- as at the start SHA, which looked every label up as a column
  # name (the spec's unaffected case: labels that are all column names give
  # identical output).  It is not parsed, so "ER-pos" is not ER - pos and
  # "bm status" is not dropped as unparseable.  Nothing unresolved: no
  # warning.
  r <- .a151_collect(test_hte_from_forestsearch(
    list(confounders.evaluated = c("bm status", "ER-pos")), df,
    "time", "treat"))
  expect_length(r$warnings, 0)
  expect_identical(r$value$covariates_used, c("bm status", "ER-pos"))
  expect_identical(r$value$n_covariates, 2L)
  X <- cbind(`bm status` = as.numeric(df[["bm status"]]),
             `ER-pos` = as.numeric(df[["ER-pos"]]))
  expect_equal(r$value$hte_test, test_hte(df$time, df$treat, X))

  # Mixed with a cut label, by hand: "bm status" -> itself (column name);
  # "age <= 55" is not a column name, so it is parsed -> age; "ER-pos" ->
  # itself.  Unique columns in label order; no warning.
  r2 <- .a151_collect(test_hte_from_forestsearch(
    list(confounders.evaluated = c("bm status", "age <= 55", "ER-pos")),
    df, "time", "treat"))
  expect_length(r2$warnings, 0)
  expect_identical(r2$value$covariates_used, c("bm status", "age", "ER-pos"))
})
