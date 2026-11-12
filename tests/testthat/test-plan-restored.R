# =============================================================================
# The caller's future plan is restored (code audit INV-2, plan half).
#
# Every exported function that sets a future plan puts the caller's plan back
# on exit, on a normal return and on an error.  The caller's plan here is a
# two-worker multisession plan, which none of the calls below sets for itself,
# so a function that resets to sequential, or leaves its own plan in place,
# fails the check.  setup_parallel_SGcons() and reset_workers() are left alone:
# setting the plan is their documented purpose.
#
# fpr_calibration() and fs_fdr_report() run their replicates in workers that
# load the installed package; their %dofuture% loops are stubbed out here
# (local_mocked_bindings), since what is tested is the plan handling around
# the loop.  mrct_region_sims() and the bootstrap run for real on a
# sequential plan of their own.
# =============================================================================

.pr_sig <- function() {
  list(class = class(future::plan()), workers = future::nbrOfWorkers())
}

# Evaluate `expr` under a caller plan of multisession(workers = 2); return the
# plan signature before and after the call and the value (or the error).  The
# plan in force when the helper is called is put back on exit.
.pr_under_caller_plan <- function(expr) {
  outer <- future::plan()
  on.exit(future::plan(outer), add = TRUE)
  future::plan(future::multisession, workers = 2L)
  before <- .pr_sig()
  value <- tryCatch(expr, error = function(e) e)
  list(before = before, after = .pr_sig(), value = value)
}

# A stand-in for the %dofuture% loop: one empty, non-declaring result per
# iteration, shaped for both fpr_calibration() and fs_fdr_report().
.pr_fake_dofuture <- function(obj, ex) {
  rep(list(list(found = 0L, flag = 0L, deb_hr = NA_real_, msg = NA_character_)),
      length(obj$args[[1L]]))
}

# Timing fields differ between any two runs; compare saved copies.
.pr_strip_t <- function(x) {
  if (is.list(x) && !is.data.frame(x)) {
    drop <- names(x) %in% c("timing_seconds", "minutes_all", "time_search")
    if (any(drop)) x <- x[!drop]
    for (i in seq_along(x)) if (!is.null(x[[i]])) x[[i]] <- .pr_strip_t(x[[i]])
  }
  x
}
.pr_rds <- function(x) {
  f <- tempfile(fileext = ".rds")
  on.exit(unlink(f), add = TRUE)
  saveRDS(x, f)
  readRDS(f)
}

.pr_seq <- list(plan = "sequential", workers = 1L, show_message = FALSE)

.pr_data <- .make_survival_data(N = 400L, HR_harm = 3.0)
.pr_fs <- function() {
  suppressWarnings(do.call(forestsearch, c(
    list(df.analysis = .pr_data),
    .fs_args_for("survival", confounders = c("age", "stage", "sex", "noise"),
                 extra = list(use_grf = FALSE, use_lasso = FALSE, quiet = TRUE)))))
}

.pr_ss <- function() {
  d <- .make_survival_data(N = 200L, seed = 42L)
  Z <- cbind(age_le50 = as.integer(d$age <= 50),
             stage3   = as.integer(d$stage == "3"),
             sex1     = as.integer(d$sex == "1"))
  suppressMessages(subgroup.search(
    Y = d$time, Event = d$event, Treat = d$treat, Z = Z,
    n.min = 20, d0.min = 5, d1.min = 5, hr.threshold = 1.0,
    maxk = 2, parallel_workers = 1))
}


test_that("subgroup.search() restores the caller's plan", {
  skip_on_cran()
  r <- .pr_under_caller_plan(.pr_ss())
  expect_false(inherits(r$value, "error"))
  expect_identical(r$after, r$before)
})

test_that("subgroup.search() restores the caller's plan when the search fails", {
  skip_on_cran()
  local_mocked_bindings(
    search_combinations_parallel = function(...) stop("search stopped"),
    .package = "forestsearch")
  r <- .pr_under_caller_plan(.pr_ss())
  expect_s3_class(r$value, "error")
  expect_match(conditionMessage(r$value), "search stopped", fixed = TRUE)
  expect_identical(r$after, r$before)
})

test_that("forestsearch() restores the caller's plan; its result does not depend on it", {
  skip_on_cran()
  r <- .pr_under_caller_plan(.pr_fs())
  expect_s3_class(r$value, "forestsearch")
  expect_identical(r$after, r$before)

  outer <- future::plan()
  future::plan(future::sequential)
  ref <- .pr_fs()
  future::plan(outer)
  expect_identical(.pr_rds(.pr_strip_t(r$value)), .pr_rds(.pr_strip_t(ref)))
})

test_that("forestsearch_bootstrap_dofuture() restores the caller's plan", {
  skip_on_cran()
  outer <- future::plan()
  future::plan(future::sequential)
  fs <- .pr_fs()
  future::plan(outer)
  skip_if(is.null(fs$sg.harm), "no subgroup identified on the primary fit")
  r <- .pr_under_caller_plan(suppressWarnings(suppressMessages(
    forestsearch_bootstrap_dofuture(fs.est = fs, nb_boots = 2L,
                                    parallel_args = .pr_seq))))
  expect_false(inherits(r$value, "error"))
  expect_identical(r$after, r$before)
})

test_that("mrct_region_sims() restores the caller's plan, also on an error", {
  skip_on_cran()
  dgm <- mrct_build_S2()
  call_rs <- function() {
    suppressWarnings(suppressMessages(mrct_region_sims(
      dgm = dgm, n_sims = 2L, n_sample = 600, region_var = "z_region",
      sg_focus = "hrMaxSG", maxk = 1, hr.threshold = 0.90, hr.consistency = 0.80,
      pconsistency.threshold = 0.90,
      confounders.name = analysis_covariates(dgm, "observed"),
      conf_force = character(0),
      fs_args = list(conf.cont_jcuts = list(z_cd40 = 10, z_wtkg = 10, z_preanti = 10),
                     use_grf = FALSE, use_lasso = FALSE, fs.splits = 50),
      sim_args = list(rand_ratio = 1, draw_treatment = TRUE, max_entry = 11,
                      cens_adjust = 50),
      analysis_time = 40, cens_adjust = 50, parallel_args = .pr_seq,
      details = FALSE, seed = 20261101L)))
  }
  r <- .pr_under_caller_plan(call_rs())
  expect_false(inherits(r$value, "error"))
  expect_identical(nrow(r$value), 2L)
  expect_identical(r$after, r$before)

  local_mocked_bindings(`%dofuture%` = function(obj, ex) stop("loop stopped"),
                        .package = "forestsearch")
  r_err <- .pr_under_caller_plan(call_rs())
  expect_s3_class(r_err$value, "error")
  expect_identical(r_err$after, r_err$before)
})

test_that("fpr_calibration(n_workers > 1) restores the caller's plan", {
  skip_on_cran()
  local_mocked_bindings(`%dofuture%` = .pr_fake_dofuture, .package = "forestsearch")
  d <- data.frame(id = 1:60, treat = rep(0:1, 30), y = rep(0:1, each = 30))
  r <- .pr_under_caller_plan(fpr_calibration(
    c1 = 1.25, c2 = 1.0, df.analysis = d, fs_params = list(treat.name = "treat"),
    nsims = 2L, n_workers = 3L, quiet = TRUE))
  expect_s3_class(r$value, "fpr_calibration")
  expect_identical(r$after, r$before)
})

test_that("fs_fdr_report() restores the caller's plan", {
  skip_on_cran()
  local_mocked_bindings(`%dofuture%` = .pr_fake_dofuture, .package = "forestsearch")
  fs_analysis <- list(args_call_all = list(hr.threshold = 1.25, hr.consistency = 1.0))
  r <- .pr_under_caller_plan(fs_fdr_report(
    fs_analysis, df_analysis = data.frame(x = seq_len(100)), null = "hom",
    nsims = 2L, n_super = 2000L, n_workers = 3L, quiet = TRUE))
  expect_s3_class(r$value, "fs_fdr_report")
  expect_identical(r$after, r$before)
})

test_that("setup_parallel_SGcons() still sets the plan (its documented purpose)", {
  skip_on_cran()
  outer <- future::plan()
  on.exit(future::plan(outer), add = TRUE)
  future::plan(future::sequential)
  setup_parallel_SGcons(list(plan = "multisession", workers = 2L, show_message = FALSE))
  expect_identical(future::nbrOfWorkers(), 2L)
})
