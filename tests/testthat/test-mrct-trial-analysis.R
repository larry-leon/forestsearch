# test-mrct-trial-analysis.R
#
# mrct_trial_analysis() and mrct_four_sets() on the core they share with
# mrct_case_study() (TASK_mrct_trial_analysis_2026-11-06, D9). The replicate-5 gates
# against the committed record run once as gate artifacts, not here.

# ---- FX-fast: a small simulated cell and one case study, regenerated in the test
.mta_fast_dgm <- function() {
  set.seed(101L); n <- 500L
  d <- data.frame(age = stats::rnorm(n, 50, 10), bm = stats::rnorm(n, 25, 4), lab = stats::rnorm(n, 300, 100),
                  f1 = stats::rbinom(n, 1L, 0.4), f2 = stats::rbinom(n, 1L, 0.5), treat = stats::rbinom(n, 1L, 0.5))
  d$tte <- stats::rexp(n, 0.05 * exp(-0.5 * d$treat + 0.01 * (d$age - 50)))
  d$event <- as.integer(d$tte < 40); d$tte <- pmin(d$tte, 40)
  suppressWarnings(inject_mrct_structure(
    d, outcome_var = "tte", event_var = "event", treatment_var = "treat",
    continuous_vars = c("age", "bm", "lab"), factor_vars = c("f1", "f2"), x_pred = "lab",
    spline_spec = list(knot = 300, zeta = 800, log_hrs = rep(log(0.55), 3)),
    region = list(prevalence = 0.25, or_pred = 1, x1_vars = "f1", or_x1 = 3, loghr = -log(3),
                  band = list(cut = 230, or = 1)),
    x3 = list(vars = "lab", cuts = list(lab = 230), loghr = log(3 / 0.55)),
    n_super = 4000, expand = "bootstrap", seed = 7L,
    continuous_vars_cens = c("age", "lab"), factor_vars_cens = "f1"))
}
.mta_fast_fs_args <- function(dgm) {
  list(sg_focus = "hrMaxSG", maxk = 1, confounders.name = analysis_covariates(dgm, "observed"),
       conf_force = character(0), use_grf = FALSE, use_lasso = FALSE, fs.splits = 200,
       conf.cont_jcuts = list(z_lab = 5, z_age = 3, z_bm = 3))
}
.mta_memo <- new.env()
.mta_fast <- function() {
  if (is.null(.mta_memo$cs)) {
    .mta_memo$dgm <- .mta_fast_dgm()
    set.seed(99L)
    .mta_memo$before <- .Random.seed
    .mta_memo$cs <- mrct_case_study(.mta_memo$dgm, seed = 5000L, k = 3, n = 600, region_var = "z_region",
                                    fs_args = .mta_fast_fs_args(.mta_memo$dgm),
                                    sim_args = list(rand_ratio = 1, draw_treatment = TRUE, max_entry = 6, cens_adjust = 20),
                                    analysis_time = 30, cens_adjust = 20, mr_inference = FALSE)
    .mta_memo$after <- .Random.seed
  }
  .mta_memo
}

# ---- tiny engineered trials, one per stage: a hazard ratio by patient, exponential
# times, administrative censoring at 30
.mta_trial <- function(hr_fun, n = 800L, seed = 1L) {
  set.seed(seed)
  d <- data.frame(id = seq_len(n), region = stats::rbinom(n, 1L, 0.25), treat = stats::rbinom(n, 1L, 0.5),
                  x1 = stats::rnorm(n), x2 = stats::rnorm(n))
  d$harm <- as.integer(d$x1 < ifelse(d$region == 1, 0.8, -0.6))
  d$time <- stats::rexp(n, 0.1 * ifelse(d$treat == 1, hr_fun(d), 1))
  d$event <- as.integer(d$time < 30)
  d$time <- pmin(d$time, 30)
  d
}
.mta_scenarios <- list(
  complete = function(d) ifelse(d$harm == 1, 3, 0.45),
  no_declaration = function(d) ifelse(d$region == 1, 1.5, 0.5),
  no_trigger = function(d) rep(0.5, nrow(d)),
  not_positive = function(d) rep(1.1, nrow(d)))
.mta_fit <- function(d, ...) {
  mrct_trial_analysis(d, "time", "event", "treat", region.name = "region", id.name = "id",
                      fs_args = list(sg_focus = "hrMaxSG", confounders.name = c("x1", "x2"),
                                     use_grf = FALSE, use_lasso = FALSE, fs.splits = 200),
                      mr_inference = FALSE, ...)
}
.mta_cox <- function(d) {
  s <- summary(survival::coxph(survival::Surv(time, event) ~ treat, data = d))
  c(hr = unname(s$conf.int[1, "exp(coef)"]), lo = unname(s$conf.int[1, "lower .95"]),
    hi = unname(s$conf.int[1, "upper .95"]), se = unname(s$coefficients[1, "se(coef)"]))
}
.mta_sets_cols <- c("set", "label", "n", "n_exp", "n_ctl", "events_exp", "events_ctl", "hr", "lo", "hi", "se_loghr")

test_that("FX-fast: mrct_case_study() keeps its fields, class, sets and print layout", {
  m <- .mta_fast(); cs <- m$cs
  expect_identical(class(cs), "mrct_case_study")
  expect_identical(names(cs), c("trial", "fs", "sets", "in_Hhat", "confusion", "Hhat", "mr", "trigger",
                                "call_args", "warnings"))
  expect_identical(names(cs$call_args), c("seed", "k", "n", "region_var", "analysis_time", "sim_args", "fs_args"))
  expect_false(any(c("df.analysis", "df.test") %in% names(cs$call_args$fs_args)))
  expect_identical(names(cs$sets), .mta_sets_cols)
  expect_identical(cs$sets$set, c("itt_all", "itt_region", "region_Gstar", "region_Ghat"))
  expect_identical(cs$sets$label, c("All patients, ITT", "Region, ITT", "Region minus H* (G*)",
                                    "Region minus H-hat (G-hat)"))
  expect_identical(names(cs$confusion), c("nonregion", "region"))
  expect_identical(names(cs$confusion$region), c("TP", "FP", "FN", "TN"))
  expect_identical(m$after, m$before)                # the caller's stream is restored
  out <- capture.output(print(cs))
  expect_identical(out[1], "MRCT case study: trial 3 (seed 5000 + 3), n = 600, region 'z_region'")
  expect_identical(out[2], paste0("H-hat: ", cs$Hhat))
  expect_identical(out[3], paste0("Trigger (ITT ub < 1.0 and region ITT > 1.0): ", cs$trigger))
  expect_true(any(out == "H-hat against H* (TP, FP, FN, TN):"))
  expect_identical(out[length(out)], "Multiplier resampling: not run, or nothing declared.")
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(cs, which = "forest"))
  expect_no_error(plot(cs, which = "km"))
})

test_that("FX-fast: the three entries agree on the same patients", {
  m <- .mta_fast(); cs <- m$cs
  expect_true("id" %in% names(cs$trial))
  obs <- mrct_trial_analysis(cs$trial, "y_sim", "event_sim", "treat_sim", region.name = "z_region",
                             id.name = "id", truth.name = "flag_harm", fs_args = .mta_fast_fs_args(m$dgm),
                             mr_inference = FALSE, proceed_without_trigger = TRUE)
  expect_identical(class(obs), "mrct_trial_analysis")
  expect_identical(obs$Hhat, cs$Hhat)
  expect_identical(obs$in_Hhat, cs$in_Hhat)
  expect_identical(obs$confusion, cs$confusion)
  expect_equal(obs$sets, cs$sets, tolerance = 1e-12)
  fs4 <- mrct_four_sets(cs$trial, "y_sim", "event_sim", "treat_sim", region.name = "z_region",
                        in_Hhat = cs$in_Hhat, truth.name = "flag_harm")
  expect_identical(as.data.frame(fs4), cs$sets)
  expect_identical(fs4$confusion, cs$confusion)
})

test_that("each stage on an engineered trial, and never an error", {
  d <- .mta_trial(.mta_scenarios$not_positive)
  x <- .mta_fit(d)
  expect_gte(.mta_cox(d)[["hi"]], 1)
  expect_identical(x$stage, "not_positive")
  expect_false(x$trigger)
  expect_null(x$fs); expect_null(x$mr); expect_null(x$in_Hhat); expect_null(x$confusion)
  expect_identical(x$Hhat, NA_character_)
  expect_identical(x$sets$set, c("itt_all", "itt_region"))
  expect_equal(x$sets$hr[1], .mta_cox(d)[["hr"]], tolerance = 1e-12)
  expect_true(any(grepl("Stage: not_positive \\(the identifier was not run\\)", capture.output(print(x)))))

  d <- .mta_trial(.mta_scenarios$no_trigger)
  x <- .mta_fit(d)
  expect_lt(.mta_cox(d)[["hi"]], 1)
  expect_lte(.mta_cox(d[d$region == 1, ])[["hr"]], 1)
  expect_identical(x$stage, "no_trigger")
  expect_false(x$trigger)
  expect_null(x$fs)
  expect_identical(x$sets$set, c("itt_all", "itt_region"))

  d <- .mta_trial(.mta_scenarios$no_declaration)
  x <- .mta_fit(d)
  expect_lt(.mta_cox(d)[["hi"]], 1)
  expect_gt(.mta_cox(d[d$region == 1, ])[["hr"]], 1)
  expect_identical(x$stage, "no_declaration")
  expect_true(x$trigger)
  expect_null(x$fs$sg.harm)
  expect_identical(x$Hhat, "none")
  expect_identical(x$sets$set, c("itt_all", "itt_region", "region_Ghat"))
  expect_identical(x$sets$n[3], x$sets$n[2])               # region(G-hat) is the whole region
  expect_false(any(x$in_Hhat))

  d <- .mta_trial(.mta_scenarios$complete)
  x <- .mta_fit(d)
  expect_identical(x$stage, "complete")
  expect_true(x$trigger)
  expect_false(identical(x$Hhat, "none"))
  expect_identical(x$sets$set, c("itt_all", "itt_region", "region_Ghat"))
  expect_equal(x$sets$hr[3], .mta_cox(d[d$region == 1 & !x$in_Hhat, ])[["hr"]], tolerance = 1e-12)
  expect_identical(x$in_Hhat,
                   get_dfpred(d, sg.harm = x$fs$sg.harm, version = 2)$treat.recommend == 0)
})

test_that("proceed_without_trigger runs the identifier and still reports the trigger", {
  d <- .mta_trial(.mta_scenarios$no_trigger)
  x <- .mta_fit(d, proceed_without_trigger = TRUE)
  expect_false(is.null(x$fs))
  expect_false(x$trigger)
  expect_true(x$stage %in% c("no_declaration", "complete"))
  expect_identical(x$sets$set[3], "region_Ghat")
  d <- .mta_trial(.mta_scenarios$not_positive)
  x <- .mta_fit(d, proceed_without_trigger = TRUE)
  expect_false(is.null(x$fs))
  expect_false(x$trigger)
  expect_true(x$stage %in% c("no_declaration", "complete"))
})

test_that("truth absent: no G* row, no confusion, no error, no message; truth present: both", {
  d <- .mta_trial(.mta_scenarios$complete)
  expect_no_message(x <- .mta_fit(d))
  expect_false("region_Gstar" %in% x$sets$set)
  expect_null(x$confusion)
  expect_false(any(grepl("H-hat against H\\*", capture.output(print(x)))))
  y <- mrct_trial_analysis(d, "time", "event", "treat", region.name = "region", id.name = "id",
                           truth.name = "harm",
                           fs_args = list(sg_focus = "hrMaxSG", confounders.name = c("x1", "x2"),
                                          use_grf = FALSE, use_lasso = FALSE, fs.splits = 200),
                           mr_inference = FALSE)
  expect_identical(y$sets$set, c("itt_all", "itt_region", "region_Gstar", "region_Ghat"))
  expect_equal(y$sets$hr[3], .mta_cox(d[d$region == 1 & d$harm == 0, ])[["hr"]], tolerance = 1e-12)
  expect_identical(sum(y$confusion$region) + sum(y$confusion$nonregion), nrow(d))
  expect_identical(y$Hhat, x$Hhat)
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(y, which = "forest"))
  expect_no_error(plot(y, which = "km", time_unit = "months", analysis_time = 30))
})

test_that("with MR, print reports the evaluated construction set, by location", {
  skip_on_cran()                                     # about 9 s
  d <- .mta_trial(.mta_scenarios$complete)
  x <- mrct_trial_analysis(d, "time", "event", "treat", region.name = "region", id.name = "id",
                           fs_args = list(sg_focus = "hrMaxSG", confounders.name = c("x1", "x2"),
                                          use_grf = FALSE, use_lasso = FALSE, fs.splits = 200),
                           mr_inference_args = list(ci_method = "field", include_complement = TRUE,
                                                    field_complement = TRUE, draws = 100,
                                                    field_R_out = 30, field_R_in = 15, seed = 7L))
  expect_identical(x$stage, "complete")
  expect_identical(x$mr, x$fs$mr_inference)
  out <- capture.output(print(x, digits = 3))
  v <- mr_values(x$fs)
  f <- function(z) formatC(z, format = "f", digits = 3)
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "field-s 1-sided 95%", fixed = TRUE)
  expect_match(txt, "Bonferroni, field-s", fixed = TRUE)
  expect_match(txt, paste(">=", f(v$field_lower_1s[1])), fixed = TRUE)
  expect_match(txt, paste("<=", f(v$field_upper_1s[2])), fixed = TRUE)
  expect_match(txt, paste("<=", f(v$field_upper_1s_s[2])), fixed = TRUE)
  expect_match(txt, paste(">=", f(v$bonf_lower_H[1])), fixed = TRUE)
  expect_match(txt, paste("<=", f(v$bonf_upper_Hc[2])), fixed = TRUE)
  expect_false(grepl("significan|n\\.s\\.|winner", txt, ignore.case = TRUE))
})

test_that("the input is checked", {
  d <- .mta_trial(.mta_scenarios$complete, n = 50L)
  fa <- list(confounders.name = c("x1", "x2"))
  expect_error(mrct_trial_analysis(as.list(d), "time", "event", "treat", "id", "region", fs_args = fa), "must be a data frame")
  expect_error(mrct_trial_analysis(d, "nope", "event", "treat", "id", "region", fs_args = fa), "outcome.name 'nope' not found")
  expect_error(mrct_trial_analysis(d, "time", "event", "treat", "id", "region", truth.name = "nope", fs_args = fa), "truth.name 'nope'")
  expect_error(mrct_trial_analysis(d, "time", "event", "treat", "id", "region"), "confounders.name is required")
  bad <- d; bad$region[1] <- 2L
  expect_error(mrct_trial_analysis(bad, "time", "event", "treat", "id", "region", fs_args = fa), "must be 0/1")
  bad <- d; bad$region[1] <- NA
  expect_error(mrct_trial_analysis(bad, "time", "event", "treat", "id", "region", fs_args = fa), "missing values")
  expect_error(mrct_trial_analysis(d, "time", "event", "treat", "nope", "region", fs_args = fa), "id.name 'nope'")
  expect_error(mrct_four_sets(d, "time", "event", "treat", "region", in_Hhat = TRUE), "one non-missing value per row")
})

test_that("mrct_four_sets() on a hand-built data frame equals direct coxph() summaries", {
  d <- .mta_trial(.mta_scenarios$complete, seed = 4L)
  H <- d$x1 < -0.2
  fs4 <- mrct_four_sets(d, "time", "event", "treat", region.name = "region", in_Hhat = H, truth.name = "harm")
  s <- as.data.frame(fs4)
  expect_s3_class(fs4, "mrct_four_sets")
  expect_identical(class(s), "data.frame")
  expect_identical(names(s), .mta_sets_cols)
  expect_identical(s$set, c("itt_all", "itt_region", "region_Gstar", "region_Ghat"))
  idx <- list(rep(TRUE, nrow(d)), d$region == 1, d$region == 1 & d$harm == 0, d$region == 1 & !H)
  for (i in seq_along(idx)) {
    z <- d[idx[[i]], ]; cx <- .mta_cox(z)
    expect_equal(c(s$hr[i], s$lo[i], s$hi[i], s$se_loghr[i]), unname(cx), tolerance = 1e-12)
    expect_identical(as.integer(s$n[i]), nrow(z))
    expect_identical(as.integer(s$n_exp[i]), sum(z$treat == 1L))
    expect_identical(as.integer(s$n_ctl[i]), sum(z$treat == 0L))
    expect_identical(as.integer(s$events_exp[i]), sum(z$event[z$treat == 1L]))
    expect_identical(as.integer(s$events_ctl[i]), sum(z$event[z$treat == 0L]))
  }
  band <- d$harm == 1; reg <- d$region == 1
  expect_identical(fs4$confusion$region[["TP"]], sum(H[reg] & band[reg]))
  expect_identical(fs4$confusion$nonregion[["TN"]], sum(!H[!reg] & !band[!reg]))

  no_truth <- mrct_four_sets(d, "time", "event", "treat", region.name = "region", in_Hhat = H)
  expect_identical(no_truth$sets$set, c("itt_all", "itt_region", "region_Ghat"))
  expect_null(no_truth$confusion)
  nothing <- mrct_four_sets(d, "time", "event", "treat", region.name = "region")
  expect_identical(nothing$sets$n[3], nothing$sets$n[2])

  out <- capture.output(print(fs4))
  expect_identical(out[1], sprintf("MRCT four sets: n = %d, region 'region'", nrow(d)))
  expect_true(any(out == "H-hat against H* (TP, FP, FN, TN):"))
  grDevices::pdf(NULL); on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(fs4))
  expect_no_error(plot(fs4, which = "forest"))
  expect_no_error(plot(fs4, which = "km", time_unit = "months", analysis_time = 30))
  expect_invisible(plot(no_truth, which = "forest"))
})
