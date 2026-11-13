# mrct_select_case() and mrct_case_study() (TASK_fs-mrct-article-r2_2026-11-01,
# section 2). On the S2 build of helper-mrct-analytic.R with the campaign's
# call seed and arguments (mrct-consistency sims/121_ident_rlog4_full_call.R):
# the selection rule on replicates 1-6 returns replicate 5 (the rule of
# sims/150_s2_example_km.R), and the case study of replicate 5 reproduces the
# committed record of results/s2_example/ (fixtures/mrct_case_study/). The
# mrct_region_sims() workers run the installed build: install before running.
# Multiplier resampling is reduced (200 draws, field 50/25) where the check
# does not depend on it.

cs_fixture <- function(f) testthat::test_path("fixtures", "mrct_case_study", f)
cs_seed <- 613800000L   # the campaign's call seed of S2

cs_fs_args <- function(dgm, ...) c(list(
  sg_focus = "hrMaxSG", maxk = 2,
  hr.threshold = 0.90, hr.consistency = 0.80, pconsistency.threshold = 0.90,
  confounders.name = analysis_covariates(dgm, "observed"), conf_force = character(0),
  conf.cont_jcuts = list(z_cd40 = 10, z_wtkg = 10, z_preanti = 10),
  effect_neighborhood = 0.20, use_grf = FALSE, use_lasso = FALSE), list(...))
cs_sim_args <- list(rand_ratio = 1, draw_treatment = TRUE, max_entry = 11, cens_adjust = 50)
cs_mr_small <- list(ci_method = "field", include_complement = TRUE, field_complement = TRUE,
                    draws = 200, field_R_out = 50, field_R_in = 25)

cs_memo <- new.env()
cs_get <- function(key, expr) {
  if (is.null(cs_memo[[key]])) cs_memo[[key]] <- expr
  cs_memo[[key]]
}

# The campaign's mrct_region_sims() call on S2, replicates 1-6.
cs_runs <- function() cs_get("runs", {
  dgm <- mrct_build_S2()
  fa <- cs_fs_args(dgm)
  direct <- c("sg_focus", "maxk", "hr.threshold", "hr.consistency", "pconsistency.threshold",
              "confounders.name", "conf_force")
  res <- do.call(mrct_region_sims, c(list(dgm = dgm, n_sims = 6L, n_sample = 1000,
                                          region_var = "z_region"),
                                     fa[direct],
                                     list(fs_args = fa[setdiff(names(fa), direct)],
                                          sim_args = cs_sim_args, analysis_time = 40,
                                          cens_adjust = 50,
                                          parallel_args = list(plan = "multisession",
                                                               workers = 2L,
                                                               show_message = FALSE),
                                          details = FALSE, seed = cs_seed, record = "full")))
  future::plan(future::sequential)
  res
})

cs_rep5 <- function() cs_get("rep5", {
  dgm <- mrct_build_S2()
  mrct_case_study(dgm, seed = cs_seed, k = 5, n = 1000, region_var = "z_region",
                  fs_args = cs_fs_args(dgm), sim_args = cs_sim_args,
                  analysis_time = 40, cens_adjust = 50,
                  mr_inference = TRUE, mr_inference_args = cs_mr_small)
})

test_that("mrct_select_case() refuses a record without the estimators attribute", {
  res <- data.frame(sim = 1:2, any_found = c(1, 0))
  expect_error(mrct_select_case(res), "record = \"full\"")
  est <- data.frame(sim = 1L, estimator = "itt_all", hr = 0.8, hi = 0.9)
  attr(res, "estimators") <- est
  expect_error(mrct_select_case(res, conditions = "nope"), "conditions must be")
  expect_error(mrct_select_case(res, itt_ub = NA), "single numbers")
})

test_that("mrct_select_case() evaluates the three conditions per trial; NA never qualifies", {
  # five stand-in trials: 1 untriggered; 2 triggered, undeclared; 3 triggered,
  # declared, Ghat >= 1; 4 all three; 5 all three; 6 a missing region HR
  mk <- function(sim, hr_all, hi_all, hr_reg, hr_gh)
    data.frame(sim = sim, estimator = c("itt_all", "itt_region", "region_Ghat"),
               hr = c(hr_all, hr_reg, hr_gh), hi = c(hi_all, NA, NA))
  est <- rbind(mk(1, 0.8, 1.02, 1.2, 0.9), mk(2, 0.8, 0.95, 1.2, 1.2), mk(3, 0.8, 0.95, 1.1, 1.05),
               mk(4, 0.8, 0.95, 1.1, 0.9), mk(5, 0.7, 0.85, 1.3, 0.7), mk(6, 0.8, 0.95, NA, 0.9))
  res <- data.frame(sim = 6:1, any_found = c(1, 1, 1, 1, 0, 1), sg_found = "x")
  attr(res, "estimators") <- est
  first <- mrct_select_case(res)
  expect_identical(first$sim, 4L)
  expect_identical(attr(first, "n_meeting"), 2L)
  all3 <- mrct_select_case(res, which = "all")
  expect_identical(all3$sim, c(4L, 5L))
  expect_true(all(all3$trigger & all3$declared & all3$region_Ghat_lt1))
  expect_identical(all3$hr_itt_region, c(1.1, 1.3))
  expect_false("hr_region_Gstar" %in% names(all3))
  trig <- mrct_select_case(res, conditions = "trigger", which = "all")
  expect_identical(trig$sim, c(2L, 3L, 4L, 5L))       # 1: ub >= 1; 6: NA region HR
  # the bounds are arguments
  expect_identical(mrct_select_case(res, conditions = "trigger", which = "all", itt_ub = 0.9)$sim, 5L)
  expect_identical(mrct_select_case(res, conditions = "trigger", which = "all", region_itt = 1.25)$sim, 5L)
  expect_message(none <- mrct_select_case(res, conditions = "trigger", itt_ub = 0.5), "no trial")
  expect_identical(nrow(none), 0L)
})

test_that("on S2 with the campaign's seed, replicates 1-6 select replicate 5", {
  mrct_skip_if_no_build()
  res <- cs_runs()
  expect_identical(sort(res$sim), 1:6)
  sel <- mrct_select_case(res)
  expect_identical(nrow(sel), 1L)
  expect_identical(sel$sim, 5L)
  expect_true(sel$trigger && sel$declared && sel$region_Ghat_lt1)
  # the flags of replicates 1-5 are the committed ones
  tr <- read.csv(cs_fixture("s2_example_trials.csv"), stringsAsFactors = FALSE)
  allf <- mrct_select_case(res, conditions = "trigger", which = "all")
  flags <- mrct_select_case(res, conditions = "declared", which = "all")
  expect_identical(as.integer(1:5 %in% allf$sim), as.integer(tr$trigger))
  expect_identical(as.integer(1:5 %in% flags$sim), as.integer(tr$declared))
  every <- mrct_select_case(res, conditions = "trigger", which = "all", itt_ub = Inf, region_itt = -Inf)
  expect_identical(every$sim, 1:6)
  expect_equal(every$hr_itt_all[1:5], tr$hr_itt, tolerance = 1e-10)
  expect_equal(every$itt_ub[1:5], tr$itt_ub, tolerance = 1e-10)
  expect_equal(every$hr_itt_region[1:5], tr$hr_region, tolerance = 1e-10)
  expect_identical(every$sg_found[1:5], tr$sg_found)
  expect_true("hr_region_Gstar" %in% names(every))
})

test_that("mrct_case_study(k = 5) reproduces the committed four sets, H-hat and 2 x 2 tables", {
  mrct_skip_if_no_build()
  cs <- cs_rep5()
  expect_s3_class(cs, "mrct_case_study")
  rec <- read.csv(cs_fixture("s2_example_four_sets.csv"), stringsAsFactors = FALSE,
                  colClasses = c(hr_3dp = "character", ci_3dp = "character"))
  expect_identical(cs$sets$set, c("itt_all", "itt_region", "region_Gstar", "region_Ghat"))
  rec <- rec[match(cs$sets$set, rec$set), ]
  expect_identical(formatC(cs$sets$hr, format = "f", digits = 3), rec$hr_3dp)
  expect_identical(sprintf("%s, %s", formatC(cs$sets$lo, format = "f", digits = 3),
                           formatC(cs$sets$hi, format = "f", digits = 3)), rec$ci_3dp)
  expect_equal(cs$sets$hr, rec$hr, tolerance = 1e-10)
  for (k in c("n", "n_exp", "n_ctl", "events_exp", "events_ctl"))
    expect_identical(as.integer(cs$sets[[k]]), as.integer(rec[[k]]), info = k)
  tr <- read.csv(cs_fixture("s2_example_trials.csv"), stringsAsFactors = FALSE)
  expect_identical(cs$Hhat, tr$sg_found[tr$sim == 5])
  expect_true(cs$trigger)
  # the 2 x 2 tables of the committed definitions
  defs <- readLines(cs_fixture("s2_example_definitions.md"), encoding = "UTF-8")
  rec2x2 <- function(lab) {
    l <- defs[grepl(paste0("^- ", lab, "( \\(|:)"), defs)]
    vapply(c(TP = "TP", FP = "FP", FN = "FN", TN = "TN"),
           function(k) as.integer(sub(sprintf(".*%s ([0-9]+).*", k), "\\1", l)), 0L)
  }
  expect_identical(cs$confusion$region, rec2x2("Region"))
  expect_identical(cs$confusion$nonregion, rec2x2("Non-region"))
  # H-hat membership is the identifier's own classification
  ids_H <- c(cs$fs$df.est$id[cs$fs$df.est$treat.recommend == 0],
             cs$fs$df.test$id[cs$fs$df.test$treat.recommend == 0])
  expect_identical(cs$in_Hhat, cs$trial$id %in% ids_H)
  # equal to the simulation's own record of the trial
  est <- as.data.frame(attr(cs_runs(), "estimators"))
  e5 <- est[est$sim == 5, ]; e5 <- e5[match(cs$sets$set, e5$estimator), ]
  expect_equal(cs$sets$hr, e5$hr, tolerance = 1e-12)
  expect_equal(cs$sets$hi, e5$hi, tolerance = 1e-12)
  # MR on H-hat, seeded seed + k
  expect_false(is.null(cs$mr))
  expect_identical(cs$mr$ci_method, "field")
  expect_identical(as.integer(cs$call_args$fs_args$mr_inference_args$seed), cs_seed + 5L)
  tab <- mr_estimates_table(cs$fs, est.scale = "hr", digits = 3)
  expect_false(is.null(attr(tab, "mr_values")))
  expect_equal(attr(tab, "mr_values")$field_bound,
               c(cs$mr$field$lower_1s, cs$mr$field$complement$upper_1s), tolerance = 1e-12)
})

test_that("mrct_case_study() restores the caller's generator and stream", {
  mrct_skip_if_no_build()
  dgm <- mrct_build_S2()
  set.seed(11)
  before <- .Random.seed; kind <- RNGkind()
  cs <- mrct_case_study(dgm, seed = cs_seed, k = 1, n = 1000, region_var = "z_region",
                        fs_args = cs_fs_args(dgm, maxk = 1, fs.splits = 200),
                        sim_args = cs_sim_args, analysis_time = 40, cens_adjust = 50,
                        mr_inference = FALSE)
  expect_identical(RNGkind(), kind)
  expect_identical(.Random.seed, before)
  # the draw is the simulation's replicate 1
  e1 <- as.data.frame(attr(cs_runs(), "estimators"))
  e1 <- e1[e1$sim == 1 & e1$estimator == "itt_all", ]
  expect_equal(cs$sets$hr[cs$sets$set == "itt_all"], e1$hr, tolerance = 1e-12)
  expect_null(cs$mr)
})

test_that("without flag_harm the G* set and the 2 x 2 tables are omitted, with a message", {
  mrct_skip_if_no_build()
  dgm <- mrct_build_S2()
  dgm$df_super$flag_harm <- NULL
  expect_message(
    cs <- mrct_case_study(dgm, seed = cs_seed, k = 5, n = 1000, region_var = "z_region",
                          fs_args = cs_fs_args(dgm, maxk = 1, fs.splits = 200),
                          sim_args = cs_sim_args, analysis_time = 40, cens_adjust = 50,
                          mr_inference = FALSE),
    "no flag_harm")
  expect_false("flag_harm" %in% names(cs$trial))
  expect_identical(cs$sets$set, c("itt_all", "itt_region", "region_Ghat"))
  expect_null(cs$confusion)
  expect_output(print(cs), "MRCT case study")
  pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(cs))
})

test_that("print() and plot() run on the replicate-5 case study", {
  mrct_skip_if_no_build()
  cs <- cs_rep5()
  out <- capture.output(print(cs))
  expect_true(any(grepl("Region minus H\\* \\(G\\*\\)", out)))
  expect_true(any(grepl("Multiplier resampling \\(field", out)))
  expect_true(any(grepl("H-hat complement", out)))
  pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(cs))
  expect_no_error(plot(cs, which = "forest"))
  expect_no_error(plot(cs, which = "km"))
})
