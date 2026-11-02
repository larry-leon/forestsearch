# mrct_region_sims(record = ) (TASK_fs-mrct-region-sims-record_2026-11-01,
# section 4): the point path is unchanged, and the full record agrees with the
# point estimates, with its own counts and with the targets recomputed on the
# redrawn trials. The call is mrct_rs_call() (helper-mrct-region-sims-record.R)
# on the S2 build (helper-mrct-analytic.R). The workers run the installed
# build: install before running this file.

rs_memo <- new.env()
rs_get <- function(key, expr) {
  if (is.null(rs_memo[[key]])) rs_memo[[key]] <- expr
  rs_memo[[key]]
}
rs_point_default <- function() rs_get("p0", mrct_rs_call(mrct_build_S2()))
rs_point <- function() rs_get("pp", mrct_rs_call(mrct_build_S2(), record = "point"))
rs_full <- function() rs_get("pf", mrct_rs_call(mrct_build_S2(), record = "full"))

# The trial of replicate r, redrawn as the workers draw it (L'Ecuyer stream,
# simulate_from_dgm(seed = seed + r)).
rs_redraw <- function(r, sim_args = list(rand_ratio = 1, draw_treatment = TRUE,
                                         max_entry = 11, cens_adjust = 50)) {
  kind <- RNGkind("L'Ecuyer-CMRG")
  on.exit(do.call(RNGkind, as.list(kind)))
  s <- do.call(simulate_from_dgm, utils::modifyList(sim_args, list(
    dgm = mrct_build_S2(), n = 1000, analysis_time = 40, cens_adjust = 50,
    seed = mrct_rs_seed + r)))
  data.table::as.data.table(s)
}

point_cols <- c("sim", "n_itt", "hr_itt", "hr_ittX", "n_train", "hr_train", "n_test",
                "hr_test", "any_found", "sg_found", "n_sg", "hr_sg", "POhr_sg", "prev_sg",
                "n_sg_train", "hr_sg_train", "POhr_sg_train", "hr_sg_null")
count_cols <- paste0(rep(c("TP", "FP", "FN", "TN"), 2), rep(c("_nonregion", "_region"), each = 4))
mr_cols <- c("mr_ci_method", "mr_debiased_H", "mr_ij_lower_H", "mr_ij_upper_H",
             "mr_field_lower_1s_H", "mr_debiased_Hc", "mr_ij_lower_Hc", "mr_ij_upper_Hc",
             "mr_field_upper_1s_Hc", "mr_field_upper_1s_s_Hc", "mr_bonf_lower_H",
             "mr_bonf_upper_Hc", "mr_p_hat_H", "mr_selection_rate")
estimator_names <- c("itt_all", "itt_nonregion", "itt_region", "region_Ghat",
                     "region_Gstar", "Hstar_nonregion")

# ---- 1. point path unchanged ------------------------------------------------

test_that("record = 'point' is the call without the argument and the pre-change call", {
  mrct_skip_if_no_build()
  p0 <- rs_point_default(); pp <- rs_point()
  expect_identical(names(pp), names(p0))
  expect_identical(as.list(pp), as.list(p0))
  expect_identical(names(p0), point_cols)
  expect_null(attr(pp, "estimators"))
  # fixture: the same call with the installed build of 6b5c85c0, whose
  # mrct_region_sims() body is that of 19b703b
  fx <- readRDS(test_path("fixtures", "mrct_region_sims_record", "point_6b5c85c0.rds"))
  expect_identical(names(pp), names(fx))
  expect_identical(as.list(pp), as.list(fx))
  # the full path leaves every point column as it is and appends the counts
  pf <- rs_full()
  expect_identical(names(pf), c(point_cols, count_cols, mr_cols))
  expect_identical(as.list(pf)[point_cols], as.list(pp))
})

# ---- 2. full path against the point estimates ---------------------------------

test_that("the estimators' hr equal hr_itt, hr_train, hr_test, hr_sg", {
  mrct_skip_if_no_build()
  pf <- rs_full(); est <- attr(pf, "estimators")
  expect_s3_class(est, "data.table")
  expect_identical(names(est), c("sim", "estimator", "n", "n_exp", "events_exp", "events_ctl",
                                 "hr", "lo", "hi", "se_loghr", "ahr", "cde", "declared"))
  expect_identical(est$estimator, rep(estimator_names, nrow(pf)))
  expect_identical(est$sim, rep(pf$sim, each = length(estimator_names)))
  pick <- function(k, v = "hr") est[[v]][est$estimator == k]
  expect_equal(pick("itt_all"), pf$hr_itt, tolerance = 1e-12)       # no strat column in S2
  expect_equal(pick("itt_nonregion"), pf$hr_train, tolerance = 1e-12)
  expect_equal(pick("itt_region"), pf$hr_test, tolerance = 1e-12)
  def <- !is.na(pf$hr_sg)
  expect_true(any(def))
  expect_equal(pick("region_Ghat")[def], pf$hr_sg[def], tolerance = 1e-12)
  expect_equal(pick("region_Ghat", "ahr")[def], pf$POhr_sg[def], tolerance = 1e-12)
  expect_identical(pick("itt_all", "declared"), as.integer(pf$any_found))
  # Wald limits and SE from the same fit
  ok <- !is.na(est$hr)
  expect_equal(log(est$hi[ok]) - log(est$hr[ok]), qnorm(0.975) * est$se_loghr[ok], tolerance = 1e-10)
  expect_equal(log(est$hr[ok]) - log(est$lo[ok]), qnorm(0.975) * est$se_loghr[ok], tolerance = 1e-10)
})

# ---- 3. counts ------------------------------------------------------------------

test_that("the counts add up to the set sizes", {
  mrct_skip_if_no_build()
  pf <- rs_full(); est <- attr(pf, "estimators")
  pick <- function(k, v) est[[v]][est$estimator == k]
  expect_true(all(est$n_exp + (est$n - est$n_exp) == est$n))
  expect_true(all(est$n_exp <= est$n & est$events_exp <= est$n_exp &
                    est$events_ctl <= est$n - est$n_exp))
  expect_identical(pick("itt_all", "n"), pf$n_itt)
  expect_identical(pick("itt_nonregion", "n"), pf$n_train)
  expect_identical(pick("itt_region", "n"), pf$n_test)
  expect_identical(pick("itt_nonregion", "n") + pick("itt_region", "n"), pick("itt_all", "n"))
  with(pf, {
    expect_identical(TP_nonregion + FP_nonregion + FN_nonregion + TN_nonregion, n_train)
    expect_identical(TP_region + FP_region + FN_region + TN_region, n_test)
    expect_identical(pick("region_Ghat", "n"), n_test - (TP_region + FP_region))
    expect_identical(pick("region_Gstar", "n"), n_test - (TP_region + FN_region))
    expect_identical(pick("Hstar_nonregion", "n"), TP_nonregion + FN_nonregion)
    dec <- any_found == 1
    expect_identical(pick("region_Ghat", "n")[dec], n_sg[dec])
    expect_identical(pick("region_Ghat", "n")[!dec], n_test[!dec])
    expect_true(all(TP_nonregion[!dec] + FP_nonregion[!dec] + TP_region[!dec] + FP_region[!dec] == 0))
  })
})

# ---- 4. targets and sets on the redrawn trials -------------------------------

test_that("ahr, cde and the sets agree with the trial frame", {
  mrct_skip_if_no_build()
  pf <- rs_full(); est <- attr(pf, "estimators")
  for (r in pf$sim) {
    s <- rs_redraw(r)
    reg <- s$z_region == 1; band <- s$flag_harm == 1
    k <- which(pf$sim == r)
    H <- if (pf$any_found[k] == 1) {
      get_dfpred(as.data.frame(s), strsplit(pf$sg_found[k], " & ", fixed = TRUE)[[1]], version = 2)$treat.recommend == 0
    } else rep(FALSE, nrow(s))
    sets <- list(itt_all = rep(TRUE, nrow(s)), itt_nonregion = !reg, itt_region = reg,
                 region_Ghat = reg & !H, region_Gstar = reg & !band, Hstar_nonregion = !reg & band)
    e <- est[est$sim == r, ]
    for (j in seq_along(sets)) {
      d <- s[sets[[j]]]
      expect_identical(e$estimator[j], names(sets)[j])
      expect_identical(e$n[j], nrow(d))
      expect_identical(e$events_exp[j], sum(d$event_sim[d$treat_sim == 1]))
      expect_equal(e$ahr[j], exp(mean(d$loghr_po)), tolerance = 1e-12)
      expect_equal(e$cde[j], mean(exp(d$theta_1)) / mean(exp(d$theta_0)), tolerance = 1e-12)
    }
    expect_identical(pf$TP_region[k], sum(H[reg] & band[reg]))
    expect_identical(pf$FN_nonregion[k], sum(!H[!reg] & band[!reg]))
  }
})

test_that("without flag_harm the G* and H* rows are omitted and the counts are NA", {
  mrct_skip_if_no_build()
  s <- rs_redraw(1L); s[, flag_harm := NULL]
  f <- forestsearch:::.mrct_full_record(s, "z_region", NULL, 1L, 0L)
  expect_identical(f$estimators$estimator, estimator_names[1:4])
  expect_identical(names(f$counts), count_cols)
  expect_true(all(is.na(unlist(f$counts))))
})

# ---- 5. entry_delay passes through sim_args --------------------------------------

test_that("sim_args entry_delay reaches the replicate: region entry in [6, 17], non-region in [0, 11]", {
  mrct_skip_if_no_build()
  sa <- list(rand_ratio = 1, draw_treatment = TRUE, max_entry = 11, cens_adjust = 50,
             entry_delay = list(var = "z_region", delay = 6))
  rd <- mrct_rs_call(mrct_build_S2(), n_sims = 1L, sim_args = sa, record = "full")
  est <- attr(rd, "estimators")
  s <- rs_redraw(1L, sa); s0 <- rs_redraw(1L)
  reg <- s$z_region == 1
  # with cens_adjust = 50 the only censoring is administrative: c_time = 40 - entry
  entry <- 40 - s$c_time
  expect_true(all(entry[reg] >= 6 - 1e-9 & entry[reg] <= 17 + 1e-9))
  expect_true(all(entry[!reg] >= -1e-9 & entry[!reg] <= 11 + 1e-9))
  expect_true(min(entry[reg]) > 1e-9)
  # the replicate's record is that of the delayed trial, not of the undelayed one
  e_reg <- est[est$estimator == "itt_region", ]
  expect_identical(e_reg$events_exp + e_reg$events_ctl, sum(s$event_sim[reg]))
  expect_false(identical(sum(s$event_sim[reg]), sum(s0$event_sim[s0$z_region == 1])))
  hr_reg <- unname(exp(coef(survival::coxph(survival::Surv(y_sim, event_sim) ~ treat_sim, data = s[reg]))))
  expect_equal(rd$hr_test, hr_reg, tolerance = 1e-12)
})
