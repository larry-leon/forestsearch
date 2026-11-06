# test-mrct-minor-helpers.R
#
# mrct_record(), fs_in_subgroup() and mrct_build_dgm() (TASK_wrapper_batch_2026-11-06,
# B5).

# ---- transplanted: the packages chunk of the MRCT consistency articles
.mmh_dd <- function() {
  dd <- c(file.path("..", "..", "inst", "extdata", "mrct-consistency"),
          file.path("inst", "extdata", "mrct-consistency"),
          system.file("extdata", "mrct-consistency", package = "forestsearch"))
  dd <- dd[nzchar(dd) & dir.exists(dd)][1]
  stopifnot(!is.na(dd), file.exists(file.path(dd, "PROVENANCE.md")))
  dd
}
.mmh_rd <- function(f, ...) read.csv(file.path(.mmh_dd(), f), stringsAsFactors = FALSE, check.names = FALSE, ...)

test_that("mrct_record() is the articles' rd(), with informative errors", {
  for (f in c("rep5_trial.csv", "mr_field_S2_rep5.csv", "predictions_full.csv"))
    expect_identical(mrct_record(f), .mmh_rd(f), label = f)
  expect_identical(mrct_record("s2_example_four_sets.csv", colClasses = c(hr_3dp = "character")),
                   .mmh_rd("s2_example_four_sets.csv", colClasses = c(hr_3dp = "character")))
  expect_error(mrct_record("nope.csv"), "no record 'nope.csv'")
  expect_error(mrct_record(c("a", "b")), "one file name")
  expect_error(forestsearch:::.mrct_record_dir(file.path(tempdir(), "no-such-dir")), "was not found")
  empty <- tempfile(); dir.create(empty); on.exit(unlink(empty, recursive = TRUE))
  expect_error(forestsearch:::.mrct_record_dir(empty), "no PROVENANCE.md")
})

.mmh_fit_memo <- new.env()
.mmh_fit <- function() {
  if (is.null(.mmh_fit_memo$fit)) {
    set.seed(1); n <- 500
    d <- data.frame(id = seq_len(n), region = stats::rbinom(n, 1, 0.25), treat = stats::rbinom(n, 1, 0.5),
                    x1 = stats::rnorm(n), x2 = stats::rnorm(n))
    d$time <- stats::rexp(n, 0.1 * ifelse(d$treat == 1, ifelse(d$x1 < -0.6, 3, 0.45), 1))
    d$event <- as.integer(d$time < 30); d$time <- pmin(d$time, 30)
    .mmh_fit_memo$d <- d
    .mmh_fit_memo$fit <- forestsearch(
      df.analysis = d[d$region == 0, ], df.test = d[d$region == 1, ], outcome.name = "time",
      event.name = "event", treat.name = "treat", id.name = "id", confounders.name = c("x1", "x2"),
      use_grf = FALSE, use_lasso = FALSE, maxk = 1, fs.splits = 200, sg_focus = "hrMaxSG",
      parallel_args = list(plan = "sequential"), quiet = TRUE, details = FALSE, plot.sg = FALSE,
      show_candidate_summary = FALSE)
  }
  .mmh_fit_memo
}

test_that("fs_in_subgroup() equals the raw idiom and the fit's own classification", {
  m <- .mmh_fit(); fs <- m$fit; trial <- m$d
  expect_false(is.null(fs$sg.harm))
  in_Hhat <- fs_in_subgroup(fs, trial)
  expect_identical(in_Hhat, get_dfpred(trial, sg.harm = fs$sg.harm, version = 2)$treat.recommend == 0)
  # transplanted: the case-sets chunk's stopifnot
  stopifnot(identical(in_Hhat, trial$id %in% c(fs$df.est$id[fs$df.est$treat.recommend == 0],
                                                fs$df.test$id[fs$df.test$treat.recommend == 0])))
  expect_gt(sum(in_Hhat), 0)
  expect_identical(fs_in_subgroup(fs, data.table::as.data.table(trial)), in_Hhat)
  none <- fs; none$sg.harm <- NULL
  expect_identical(fs_in_subgroup(none, trial), rep(FALSE, nrow(trial)))
  expect_error(fs_in_subgroup(list(sg.harm = "x"), trial), "forestsearch\\(\\) fit")
  expect_error(fs_in_subgroup(fs, as.list(trial)), "data frame")
})

.mmh_args <- function() {
  list(outcome_var = "y_sim", event_var = "event_sim", treatment_var = "treat_sim",
       continuous_vars = c("age", "wtkg", "karnof", "preanti", "cd40", "cd80"),
       factor_vars = c("hemo", "homo", "drugs", "race", "gender", "symptom"),
       x_pred = "cd40",
       spline_spec = list(knot = 300, zeta = 800, log_hrs = rep(log(0.55), 3)),
       region = list(observed = "region"),
       x3 = list(vars = "cd40", cuts = list(cd40 = 230), loghr = log(2 / 0.55)),
       n_super = 20000, expand = "copula", seed = 7L,
       continuous_vars_cens = c("age", "cd40"), factor_vars_cens = "symptom")
}
.mmh_strip <- function(x) { attr(x, "build_warnings") <- NULL; x }

test_that("mrct_build_dgm() captures the warnings of a build known to warn, equal to the direct call", {
  skip_if_not_installed("latentcor")
  # the application article's pre-build chunk: its censoring-model warnings
  trial_pre <- mrct_record("rep5_trial.csv")
  trial_pre <- trial_pre[, setdiff(names(trial_pre), c("id", "flag_harm"))]
  names(trial_pre) <- sub("^z_", "", names(trial_pre))
  w_pre <- character(0)
  direct <- withCallingHandlers(
    do.call(inject_mrct_structure, c(list(trial_pre), .mmh_args())),
    warning = function(cnd) { w_pre <<- c(w_pre, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  expect_no_warning(wrapped <- do.call(mrct_build_dgm, c(list(trial_pre), .mmh_args())))
  expect_gt(length(w_pre), 0L)
  expect_identical(attr(wrapped, "build_warnings"), w_pre)
  expect_true(isTRUE(all.equal(.mmh_strip(wrapped), direct, tolerance = 0)))
  expect_identical(wrapped$df_super, direct$df_super)
})

test_that("mrct_build_dgm(m = ) scales as the direct calls, and a clean build keeps no warnings", {
  skip_if_not_installed("speff2trial")
  # the simulations article's reader-data and build-fn chunks
  actg <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))
  actg$treat_mrct <- as.integer(actg$arms == 1L)
  actg$time_months <- actg$days / 30.4375
  rownames(actg) <- NULL
  latent_R <- as.matrix(mrct_record("latent_R.csv", row.names = 1))
  args <- list(outcome_var = "time_months", event_var = "cens", treatment_var = "treat_mrct",
               continuous_vars = c("age", "wtkg", "karnof", "preanti", "cd40", "cd80"),
               factor_vars = c("hemo", "homo", "drugs", "race", "gender", "symptom"),
               x_pred = "cd40", spline_spec = list(knot = 300, zeta = 800, log_hrs = rep(log(0.55), 3)),
               region = list(prevalence = 0.20, or_pred = 1, x1_vars = "symptom", or_x1 = 3,
                             loghr = -log(4), band = list(cut = 230, or = 1)),
               x3 = list(vars = "cd40", cuts = list(cd40 = 230), loghr = log(2 / 0.55)),
               region_treat_loghr = NULL, n_super = 20000, expand = "copula", seed = 7L,
               continuous_vars_cens = c("age", "cd40"), factor_vars_cens = "symptom",
               latent_R = latent_R)
  direct0 <- do.call(inject_mrct_structure, c(list(actg), args))
  wrapped0 <- do.call(mrct_build_dgm, c(list(actg), args))
  expect_identical(attr(wrapped0, "build_warnings"), character(0))
  expect_true(isTRUE(all.equal(.mmh_strip(wrapped0), direct0, tolerance = 0)))
  direct3 <- scale_dgm_hazard(direct0, m = 3)
  wrapped3 <- do.call(mrct_build_dgm, c(list(actg), args, list(m = 3)))
  expect_identical(attr(wrapped3, "build_warnings"), character(0))
  expect_true(isTRUE(all.equal(.mmh_strip(wrapped3), direct3, tolerance = 0)))
  expect_identical(wrapped3$model_params, direct3$model_params)
  expect_false(isTRUE(all.equal(wrapped3$model_params$mu, wrapped0$model_params$mu)))
  expect_error(do.call(mrct_build_dgm, c(list(actg), args, list(m = 0))), "positive")
})
