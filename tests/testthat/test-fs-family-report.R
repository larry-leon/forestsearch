# =============================================================================
# fs_family_report(): shape, resolution agreement with the engine, coverage of
# forestsearch()'s formals, status sensitivity, and the fitted-object path.
# =============================================================================

.FR_STATUSES <- c("deterministic", "disabled", "inert", "data-dependent",
                  "data-dependent (not disableable)")

# The drivers' own argument set (sim_fs_maxeffCons_mr_md40_*), minus data names.
.fr_driver_args <- function(...) {
  utils::modifyList(list(
    confounders.name = c("age", "preanti", "wtkg", "karnof", "cd40", "cd80",
                         "hemo", "homo", "drugs", "race", "gender", "symptom", "str2"),
    outcome_type = "continuous", effect_measure = "MD",
    effect.threshold = 30, consistency.threshold = 10,
    pconsistency.threshold = 0.90, fs.splits = 400L, n.min = 60L,
    d0.min = 12L, d1.min = 12L, maxk = 2L, vi.grf.min = -0.2,
    sg_focus = "maxeffCons", selection_rule = "neighborhood",
    effect_neighborhood = 0.10, stop_threshold = NULL,
    consistency_method = "resample",
    conf.cont_jcuts = list(age = 10, preanti = 10),
    use_lasso = FALSE, use_dina = FALSE, use_grf = FALSE, use_twostage = TRUE,
    is.RCT = TRUE, adverse_outcome = FALSE), list(...), keep.null = TRUE)
}

# forestsearch() formals the report deliberately does NOT classify.  Kept here,
# not in the function: a new formal that is neither classified nor listed
# below fails test 3, forcing a decision.
.FR_OUT_OF_SCOPE <- c(
  # data and column names
  "df.analysis", "outcome.name", "event.name", "treat.name", "id.name",
  "potentialOutcome.name", "flag_harm.name", "confounders.name", "df.predict",
  "df.test", "offset.name",
  # design / estimand, not family
  "is.RCT", "est.scale", "outcome_type", "effect_measure", "adverse_outcome",
  "overdispersion", "grf_count_transform", "adjust_covariates",
  "ps_method", "ps_adjust_method", "ps_hat",
  # supplied screening results / engine tuning (the engine rows carry the switch)
  "grf_res", "grf_cuts", "dina_res", "dina_cuts", "dina_args",
  "dina_select_statistic", "grf_depth", "grf_selection",
  "grf_select_statistic", "dmin.grf", "frac.tau", "tune_grf",
  "return_selected_cuts_only", "replace_med_grf",
  # execution, printing, plotting
  "parallel_args", "seedit", "details", "quiet", "show_candidate_summary",
  "max_print", "by.risk", "plot.sg", "plot.grf",
  # post-selection inference
  "mr_inference", "mr_inference_args")

test_that("shape and contract: statuses, real formals, prints, with and without data", {
  r <- fs_family_report(.fr_driver_args())
  expect_s3_class(r, "fs_family_report")
  expect_s3_class(r, "data.frame")
  expect_identical(names(r), c("stage", "arguments", "values", "status", "note"))
  expect_true(all(r$status %in% .FR_STATUSES))
  expect_false(anyDuplicated(r$stage) > 0L)
  args_named <- unique(unlist(strsplit(r$arguments[nzchar(r$arguments)], ", ")))
  expect_true(all(args_named %in% names(formals(forestsearch))),
              info = paste("not formals:", paste(setdiff(args_named, names(formals(forestsearch))), collapse = ", ")))
  expect_type(attr(r, "verdict"), "character")
  expect_identical(sum(attr(r, "status_counts")), nrow(r))
  expect_false(attr(r, "data_supplied"))
  expect_output(print(r), "Candidate family is data-dependent")
  expect_output(print(r), "Intrinsic to the method")
  # with data
  set.seed(1)
  df <- data.frame(age = round(rnorm(300, 40, 8)), preanti = round(rexp(300, 1 / 500)),
                   wtkg = round(rnorm(300, 75, 12)), karnof = sample(c(70, 80, 90, 100), 300, TRUE),
                   cd40 = round(rnorm(300, 350, 100)), cd80 = round(rnorm(300, 900, 300)),
                   hemo = rbinom(300, 1, .1), homo = rbinom(300, 1, .6), drugs = rbinom(300, 1, .15),
                   race = rbinom(300, 1, .3), gender = rbinom(300, 1, .8), symptom = rbinom(300, 1, .2),
                   str2 = rbinom(300, 1, .6))
  r2 <- fs_family_report(.fr_driver_args(), data = df)
  expect_true(attr(r2, "data_supplied"))
  expect_true(is.integer(attr(r2, "n_cut_columns")) && attr(r2, "n_cut_columns") > 0L)
  expect_true(attr(r2, "n_combinations") > attr(r2, "n_cut_columns"))
  expect_output(print(r2), "cut columns")
  # errors
  expect_error(fs_family_report(42), "named list")
  expect_error(fs_family_report(list(confounders.name = "age")), "outcome_type")
  expect_error(fs_family_report(list(not_a_formal = 1), outcome_type = "continuous"), "not forestsearch\\(\\) formals")
})

test_that("resolution agreement (drift guard): stop_threshold and the max_n_confounders cap match the engine", {
  df <- .make_continuous_data(N = 300L, MD_harm = 1.5, seed = 42L)
  # sg_focus = "eff" (alias of "hr"), stop_threshold left at its default:
  # forestsearch() resets it to NULL and re-syncs args_call_all
  # (forestsearch_main.R L1575-1609).  The report must resolve the same.
  fit <- suppressWarnings(suppressMessages(forestsearch(
    df.analysis = df, confounders.name = c("age", "biomarker", "biomarker_hi", "sex"),
    outcome.name = "y", event.name = "y", treat.name = "treat", id.name = "id",
    outcome_type = "continuous", effect_measure = "MD", adverse_outcome = TRUE,
    effect.threshold = 0.5, consistency.threshold = 0.0, pconsistency.threshold = 0.5,
    fs.splits = 20, n.min = 30, maxk = 2, use_lasso = FALSE, use_grf = FALSE,
    vi.grf.min = NULL, sg_focus = "eff", is.RCT = TRUE, details = FALSE,
    quiet = TRUE, seedit = 1L, parallel_args = list(plan = "sequential"))))
  expect_null(fit$args_call_all$stop_threshold)
  r <- fs_family_report(fit)
  es <- r[r$stage == "early stopping", ]
  expect_identical(es$status, "disabled")
  expect_match(es$values, "stop_threshold = NULL")
  # the report's resolved value equals the engine's recorded one
  rep_resolved <- if (grepl("stop_threshold = NULL", es$values)) NULL else NA
  expect_identical(rep_resolved, fit$args_call_all$stop_threshold)
  # max_n_confounders inert under vi.grf.min = NULL, on the fitted object and on its args
  expect_null(fit$args_call_all$vi.grf.min)
  expect_identical(r[r$stage == "confounder cap", "status"], "inert")
  expect_identical(r[r$stage == "GRF variable-importance ordering", "status"], "disabled")
})

test_that("coverage guard: every forestsearch() formal is classified or explicitly out of scope", {
  r <- fs_family_report(.fr_driver_args())
  classified <- unique(unlist(strsplit(r$arguments[nzchar(r$arguments)], ", ")))
  all_formals <- names(formals(forestsearch))
  expect_identical(setdiff(all_formals, c(classified, .FR_OUT_OF_SCOPE)), character(0))
  expect_identical(intersect(classified, .FR_OUT_OF_SCOPE), character(0))
  expect_identical(setdiff(.FR_OUT_OF_SCOPE, all_formals), character(0))
})

test_that("status sensitivity: switches flip their rows; the intrinsic rows never move", {
  row_status <- function(r, stage) r[r$stage == stage, "status"]
  r_null <- fs_family_report(.fr_driver_args(vi.grf.min = NULL))
  expect_identical(row_status(r_null, "GRF variable-importance ordering"), "disabled")
  expect_identical(row_status(r_null, "confounder cap"), "inert")
  r_num <- fs_family_report(.fr_driver_args(vi.grf.min = -0.2))
  expect_identical(row_status(r_num, "GRF variable-importance ordering"), "data-dependent")
  expect_identical(row_status(r_num, "confounder cap"), "data-dependent")
  expect_identical(row_status(fs_family_report(.fr_driver_args(use_lasso = TRUE)), "LASSO screen"), "data-dependent")
  expect_identical(row_status(fs_family_report(.fr_driver_args(use_lasso = FALSE)), "LASSO screen"), "disabled")
  expect_identical(row_status(fs_family_report(.fr_driver_args(use_grf = TRUE)), "DINA / GRF subgroup paths"), "data-dependent")
  # maxeff: floors and consistency relaxed, dedup becomes exact-membership but stays intrinsic
  r_me <- fs_family_report(.fr_driver_args(sg_focus = "maxeff"))
  expect_identical(row_status(r_me, "per-factor prevalence floor"), "disabled")
  expect_identical(row_status(r_me, "consistency screen"), "disabled")
  expect_identical(row_status(r_me, "effect screen"), "disabled")
  expect_identical(row_status(r_me, "candidate cap"), "inert")
  expect_match(r_me[r_me$stage == "near-duplicate removal", "values"], "exact-membership")
  # outcome type moves the per-arm row
  expect_identical(row_status(fs_family_report(.fr_driver_args(outcome_type = "count")), "per-arm floors"), "inert")
  expect_identical(row_status(fs_family_report(.fr_driver_args(outcome_type = "binary")), "per-arm floors"), "data-dependent")
  # early stopping: meaningful for maxeffCons only
  expect_identical(row_status(fs_family_report(.fr_driver_args(stop_threshold = 0.95)), "early stopping"), "data-dependent")
  expect_identical(row_status(fs_family_report(.fr_driver_args(sg_focus = "eff", stop_threshold = 0.95)), "early stopping"), "disabled")
  # the claim the function exists to make: under EVERY combination tried, the
  # cut-construction and near-duplicate rows are data-dependent and not disableable
  combos <- list(
    .fr_driver_args(),
    .fr_driver_args(vi.grf.min = NULL),
    .fr_driver_args(vi.grf.min = NULL, use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE),
    .fr_driver_args(sg_focus = "maxeff"),
    .fr_driver_args(sg_focus = "eff", stop_threshold = NULL),
    .fr_driver_args(sg_focus = "hrMaxSG"),
    .fr_driver_args(outcome_type = "survival", effect_measure = "HR", hr.threshold = 1.25, hr.consistency = 1.0),
    .fr_driver_args(outcome_type = "binary", effect_measure = "OR"),
    .fr_driver_args(minp = 0, max_subgroups_search = Inf, m1.threshold = Inf),
    .fr_driver_args(conf_force = c("age <= 34", "preanti <= 744.5"), conf.cont_jcuts = NULL))
  for (a in combos) {
    r <- fs_family_report(a)
    expect_identical(row_status(r, "cut construction"), "data-dependent (not disableable)",
                     info = paste("cut construction under", deparse(a$sg_focus), deparse(a$vi.grf.min)))
    expect_identical(row_status(r, "near-duplicate removal"), "data-dependent (not disableable)",
                     info = paste("dedup under", deparse(a$sg_focus)))
    expect_identical(row_status(r, "subgroup size floor"), "data-dependent (not disableable)")
    expect_identical(row_status(r, "redundancy (rmin)"), "data-dependent (not disableable)")
    expect_true(attr(r, "status_counts")[["data-dependent (not disableable)"]] >= 4L)
    expect_match(attr(r, "verdict"), "intrinsic to the method")
  }
})

test_that("fitted-object path gives the same table as its args_call_all", {
  df <- .make_continuous_data(N = 300L, MD_harm = 1.5, seed = 42L)
  fit <- suppressWarnings(suppressMessages(forestsearch(
    df.analysis = df, confounders.name = c("age", "biomarker", "biomarker_hi", "sex"),
    outcome.name = "y", event.name = "y", treat.name = "treat", id.name = "id",
    outcome_type = "continuous", effect_measure = "MD", adverse_outcome = TRUE,
    effect.threshold = 0.5, consistency.threshold = 0.0, pconsistency.threshold = 0.5,
    fs.splits = 20, n.min = 30, maxk = 2, use_lasso = FALSE, use_grf = FALSE,
    vi.grf.min = -0.2, sg_focus = "maxeffCons", is.RCT = TRUE, details = FALSE,
    quiet = TRUE, seedit = 1L, parallel_args = list(plan = "sequential"))))
  r_fit  <- fs_family_report(fit)
  r_args <- fs_family_report(fit$args_call_all, outcome_type = fit$args_call_all$outcome_type)
  expect_identical(as.data.frame(r_fit), as.data.frame(r_args))
  expect_identical(attr(r_fit, "verdict"), attr(r_args, "verdict"))
  expect_identical(attr(r_fit, "outcome_type"), "continuous")
})
