# test-mrct-default-confounders.R
#
# A-125 (TASK_cran_prep_code_fixes_2026-11-12, F9): when confounders.name is
# not given, mrct_region_sims() and mrct_case_study() take their default
# confounder set from dgm$df_super. For a mechanism built with a spline, the
# three internal spline columns <var>_treat, <var>_k and <var>_k_treat are
# not baseline covariates and are left out of that default; the region
# indicator stays; a covariate whose own name ends in "_treat"
# (z_prior_treat) stays, since the exclusion names the three columns exactly.

# A small synthetic seed trial with the columns create_dgm_for_mrct() requires
# (validate_mrct_data(): tte, event, treat, age, bm, male, histology,
# prior_treat, regA). All factors are 0/1.
.mdc_dgm <- function() {
  set.seed(11L); n <- 400L
  d <- data.frame(tte = stats::rexp(n, 0.05), event = stats::rbinom(n, 1L, 0.75),
                  treat = rep(0:1, n / 2L),
                  age = round(stats::rnorm(n, 62, 8)), bm = round(stats::runif(n, 0, 12), 1),
                  male = stats::rbinom(n, 1L, 0.6), histology = stats::rbinom(n, 1L, 0.4),
                  prior_treat = stats::rbinom(n, 1L, 0.3), regA = stats::rbinom(n, 1L, 0.3))
  suppressMessages(create_dgm_for_mrct(d, model_type = "alt", log_hrs = log(c(3, 1.25, 0.5)),
                                       verbose = FALSE))
}

# The same cheap forestsearch() settings for both functions, so that the two
# draw and search the same trial (seed 100, trial k = 1, n = 300).
.mdc_fs <- list(use_grf = FALSE, use_lasso = FALSE, fs.splits = 50)

test_that("A-125: the default confounders leave out the spline helper columns, in both functions", {
  skip_on_cran()                                     # about 5 s
  dgm <- .mdc_dgm()

  # Precondition: the mechanism has spline information on z_bm
  # (create_dgm_for_mrct() sets spline_spec var = "z_bm"), and
  # create_spline_variables() has put its three helper columns, named
  # z_bm + "_treat", "_k", "_k_treat", into df_super.
  expect_identical(dgm$model_params$spline_info$var, "z_bm")
  expect_true(all(c("z_bm_treat", "z_bm_k", "z_bm_k_treat") %in% names(dgm$df_super)))

  # Expected default, written out by hand: create_dgm_for_mrct() passes
  # continuous_vars = c("age", "bm") and factor_vars = c("male", "histology",
  # "prior_treat", "regA") to generate_aft_dgm_flex(), which carries each as a
  # z_ column of df_super, continuous first, then factors, in that order
  # (generate_aft_dgm_helpers.R: process_continuous_vars() then
  # process_factor_vars(), each looping over its vector; a two-level factor is
  # one 0/1 indicator column). Those six are the baseline covariates. The
  # region indicator z_regA stays (F9: "The region indicator stays").
  # z_prior_treat is a baseline covariate and stays, although its name ends
  # in "_treat". The three spline helpers above are not in the list.
  # The other columns of df_super are outside the default for another
  # reason: both functions keep only names matching "^z_|^ecog|^strat", so
  # the raw copies age, bm, male, histology, prior_treat and regA, and
  # zcens_age and zcens_prior_treat ("zcens_" does not match "^z_"), never
  # enter it, nor do id, y, treat, event, flag_harm, lin_pred_*, theta_* and
  # loghr_po. Among the matching names only the three spline helpers are
  # removed, by exact name.
  expected <- c("z_age", "z_bm", "z_male", "z_histology", "z_prior_treat", "z_regA")

  # mrct_region_sims(): confounders.name omitted (formal default NULL). With
  # details = TRUE the replicate prints the set it passes to forestsearch()
  # ("  Confounders: " followed by the names joined by ", ").
  out <- capture.output(suppressMessages(
    rs <- mrct_region_sims(dgm, n_sims = 1, n_sample = 300, region_var = "z_regA",
                           conf_force = character(0), fs_args = .mdc_fs,
                           details = TRUE, verbose_n_sims = 1,
                           parallel_args = list(plan = "sequential", workers = 1,
                                                show_message = FALSE),
                           seed = 100)))
  line <- grep("^\\s*Confounders:", out, value = TRUE)
  # Exactly one line: one replicate (n_sims = 1), and with details = TRUE the
  # line is printed for each replicate sim <= verbose_n_sims (= 1) that runs
  # forestsearch() (inside if (!is.na(hr_test)) in mrct_region_sims()).
  expect_length(line, 1L)
  rs_set <- strsplit(trimws(sub("^\\s*Confounders:", "", line)), ", ", fixed = TRUE)[[1]]
  expect_identical(rs_set, expected)
  expect_true("z_prior_treat" %in% rs_set)

  # mrct_case_study(): fs_args without confounders.name, so it takes the
  # default; call_args$fs_args is the forestsearch() argument list it used
  # (documented in its @return).
  cs <- suppressMessages(
    mrct_case_study(dgm, seed = 100, k = 1, n = 300, region_var = "z_regA",
                    fs_args = c(list(conf_force = character(0)), .mdc_fs),
                    mr_inference = FALSE))
  cs_set <- cs$call_args$fs_args$confounders.name
  expect_identical(cs_set, expected)
  expect_true("z_prior_treat" %in% cs_set)

  # The two defaults are the same, so the case study reproduces the
  # replicate's identifier (mrct_case_study() documentation: a run's direct
  # arguments and its fs_args together reproduce its identifier).
  # A regression guard, not a check of the fix: it passes before and after
  # the fix here, since this trial's identifier (!{z_male} on both trees)
  # uses no spline helper column.
  expect_identical(cs$Hhat, as.character(rs$sg_found))
})
