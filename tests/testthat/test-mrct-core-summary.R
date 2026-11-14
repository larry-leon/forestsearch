# test-mrct-core-summary.R
#
# mrct_core_summary() (TASK_wrapper_batch_2026-11-06, B1): on a hand-built record the
# rows equal the demo-summary chunk of the simulations article, transplanted below,
# and the counted values worked out by hand. The gate on a real campaign ran once.

skip_on_cran()

# Twelve trials: 3 not positive; 2 positive, region <= 1; 7 triggered, of which 5
# declared (Ghat < 1, < 1, >= 1, >= 1, NA), one of them with an undefined Sens.
.mcs_record <- function() {
  sim <- 12:1                                   # not in order: the rows are matched by sim
  mk <- function(sim, hi, hr_reg, hr_gh, hr_gs, n_gh) {
    data.frame(sim = sim, estimator = c("itt_all", "itt_region", "region_Ghat", "region_Gstar"),
               n = c(100, 25, n_gh, 20), hr = c(0.8, hr_reg, hr_gh, hr_gs), hi = c(hi, NA, NA, NA))
  }
  est <- rbind(
    mk(1, 1.05, 1.2, 1.2, 0.9, 25), mk(2, 1.01, 0.9, 0.9, 0.8, 25), mk(3, NA, 1.1, 1.1, 0.8, 25),
    mk(4, 0.95, 0.9, 0.9, 0.8, 25), mk(5, 0.90, 1.0, 1.0, 0.8, 25),
    mk(6, 0.95, 1.2, 0.9, 0.8, 20), mk(7, 0.92, 1.3, 0.7, 0.9, 15), mk(8, 0.90, 1.1, 1.0, 1.2, 22),
    mk(9, 0.91, 1.4, 1.2, 0.9, 24), mk(10, 0.93, 1.2, NA, 1.1, 10), mk(11, 0.94, 1.1, 1.1, 0.7, 25),
    mk(12, 0.96, 1.5, 1.5, 1.3, 25))
  res <- data.frame(sim = sim, any_found = as.numeric(sim %in% c(3, 4, 6, 7, 8, 9, 10)))
  tp <- ifelse(res$sim == 9, 0L, 8L); fn <- ifelse(res$sim == 9, 0L, 2L)
  res$TP_nonregion <- tp; res$FN_nonregion <- fn; res$FP_nonregion <- 1L; res$TN_nonregion <- 9L
  res$TP_region <- 4L; res$FN_region <- 1L; res$FP_region <- 0L; res$TN_region <- 5L
  attr(res, "estimators") <- est
  res
}

test_that("every row equals the article's demo-summary definitions on the same object", {
  demo <- .mcs_record(); est <- attr(demo, "estimators"); nsims_demo <- nrow(demo)
  pc <- function(x, d = 1) formatC(100 * x, format = "f", digits = d)
  # ---- transplanted: vignettes/articles/mrct-consistency-simulations.qmd, demo-summary
  g_est <- function(k) { y <- est[est$estimator == k, ]; y[match(demo$sim, y$sim), ] }
  itt <- g_est("itt_all"); rgn <- g_est("itt_region"); gh <- g_est("region_Ghat"); gs <- g_est("region_Gstar")
  lt1 <- function(v) !is.na(v) & v < 1; gt1 <- function(v) !is.na(v) & v > 1
  A <- lt1(itt$hi); K <- A & gt1(rgn$hr); dec <- demo$any_found == 1; DK <- dec & K
  ss <- function(tp, fn) tp / (tp + fn)
  sens_non <- with(demo, ss(TP_nonregion, FN_nonregion)); spec_non <- with(demo, ss(TN_nonregion, FP_nonregion))
  sens_reg <- with(demo, ss(TP_region, FN_region)); spec_reg <- with(demo, ss(TN_region, FP_region))
  pk <- function(w, d = K) if (sum(d)) pc(sum(w & d) / sum(d)) else "\u2013"
  av <- function(v) if (sum(DK)) pc(mean(v[DK], na.rm = TRUE)) else "\u2013"
  demonstration <- c(pc(mean(A)), sprintf("%d of %d", sum(K), nsims_demo), pk(dec),
                     paste(av(sens_non), "/", av(spec_non)), paste(av(sens_reg), "/", av(spec_reg)),
                     pk(dec & lt1(gh$hr)), pk(lt1(gs$hr)), pk(dec & !lt1(gh$hr)))
  # ---- the function
  s <- mrct_core_summary(demo)
  v <- stats::setNames(s$display, s$quantity)
  ours <- c(v[["power"]], v[["trigger"]], v[["declared"]],
            paste(v[["sens_nonregion"]], "/", v[["spec_nonregion"]]),
            paste(v[["sens_region"]], "/", v[["spec_region"]]),
            v[["declared_Ghat_lt1"]], v[["Gstar_lt1"]], v[["declared_Ghat_ge1"]])
  expect_identical(ours, demonstration)
})

test_that("the counts, denominators, estimates and SEs by hand", {
  s <- mrct_core_summary(.mcs_record())
  expect_identical(class(s), "data.frame")
  expect_identical(names(s), c("quantity", "label", "estimate", "se", "count", "denominator", "display"))
  expect_identical(s$quantity, c("power", "trigger", "declared", "sens_nonregion", "spec_nonregion",
                                 "sens_region", "spec_region", "declared_Ghat_lt1", "Gstar_lt1",
                                 "declared_Ghat_ge1", "Ghat_share"))
  r <- function(q) s[s$quantity == q, ]
  # A: trials 4-12 (hi < 1); K: 6-12 (region > 1)
  expect_identical(c(r("power")$count, r("power")$denominator), c(9, 12))
  expect_equal(r("power")$estimate, 75, tolerance = 1e-12)
  expect_equal(r("power")$se, 100 * sqrt(0.75 * 0.25 / 12), tolerance = 1e-12)
  expect_identical(c(r("trigger")$count, r("trigger")$denominator), c(7, 12))
  expect_identical(r("trigger")$display, "7 of 12")
  expect_true(is.na(r("trigger")$se))
  # declared among K: 6, 7, 8, 9, 10
  expect_identical(c(r("declared")$count, r("declared")$denominator), c(5, 7))
  expect_equal(r("declared")$se, 100 * sqrt((5 / 7) * (2 / 7) / 7), tolerance = 1e-12)
  # declared and Ghat < 1: 6, 7; a missing HR (10) counts as >= 1, as in the article
  expect_identical(c(r("declared_Ghat_lt1")$count, r("declared_Ghat_ge1")$count), c(2, 3))
  # G* < 1 among K: 6, 7, 9, 11
  expect_identical(c(r("Gstar_lt1")$count, r("Gstar_lt1")$denominator), c(4, 7))
  # Sens non-region over 6, 7, 8, 10: trial 9 has TP + FN = 0, undefined
  expect_identical(r("sens_nonregion")$denominator, 4)
  expect_equal(r("sens_nonregion")$estimate, 80, tolerance = 1e-12)
  expect_equal(r("sens_nonregion")$se, 0, tolerance = 1e-12)
  expect_identical(r("spec_nonregion")$denominator, 5)
  # the region's share in Ghat over 6, 7, 8, 9, 10
  expect_equal(r("Ghat_share")$estimate, 100 * mean(c(20, 15, 22, 24, 10) / 25), tolerance = 1e-12)
  expect_true(is.na(r("Ghat_share")$se))
})

test_that("an empty denominator, no truth, and a result without the record", {
  res <- .mcs_record(); est <- attr(res, "estimators")
  est$hi[est$estimator == "itt_all"] <- 1.5                     # no positive trial
  attr(res, "estimators") <- est
  s <- mrct_core_summary(res)
  r <- s[s$quantity == "declared", ]
  expect_true(is.na(r$estimate)); expect_identical(r$display, "–"); expect_identical(r$denominator, 0)
  expect_equal(s$estimate[s$quantity == "power"], 0, tolerance = 1e-12)

  res <- .mcs_record(); est <- attr(res, "estimators")
  attr(res, "estimators") <- est[est$estimator != "region_Gstar", ]
  res[grep("^(TP|FP|FN|TN)_", names(res))] <- NA_integer_
  s <- mrct_core_summary(res)
  expect_true(all(is.na(s$estimate[s$quantity %in% c("sens_nonregion", "spec_region", "Gstar_lt1")])))
  expect_identical(s$count[s$quantity == "declared"], 5)

  expect_error(mrct_core_summary(data.frame(sim = 1)), "record = \"full\"")
  bad <- data.frame(sim = 1); attr(bad, "estimators") <- data.frame(sim = 1)
  expect_error(mrct_core_summary(bad), "lacks the columns")
})

test_that("a real record = \"full\" campaign: equal to the article's rows", {
  skip_on_cran()                                     # about 5 s
  set.seed(1); n <- 400
  d <- data.frame(age = stats::rnorm(n, 50, 10), lab = stats::rnorm(n, 300, 100),
                  f1 = stats::rbinom(n, 1, 0.4), treat = stats::rbinom(n, 1, 0.5))
  d$tte <- stats::rexp(n, 0.05 * exp(-0.5 * d$treat))
  d$event <- as.integer(d$tte < 40); d$tte <- pmin(d$tte, 40)
  dgm <- suppressWarnings(inject_mrct_structure(
    d, outcome_var = "tte", event_var = "event", treatment_var = "treat",
    continuous_vars = c("age", "lab"), factor_vars = "f1", x_pred = "lab",
    spline_spec = list(knot = 300, zeta = 800, log_hrs = rep(log(0.55), 3)),
    region = list(prevalence = 0.25, or_pred = 1, x1_vars = "f1", or_x1 = 3, loghr = -log(3),
                  band = list(cut = 230, or = 1)),
    x3 = list(vars = "lab", cuts = list(lab = 230), loghr = log(3 / 0.55)),
    n_super = 2000, expand = "bootstrap", seed = 7L,
    continuous_vars_cens = c("age", "lab"), factor_vars_cens = "f1"))
  dgm <- scale_dgm_hazard(dgm, m = 3)
  demo <- mrct_region_sims(
    dgm, n_sims = 8, n_sample = 800, region_var = "z_region", sg_focus = "hrMaxSG",
    confounders.name = analysis_covariates(dgm, "observed"), conf_force = character(0),
    fs_args = list(use_grf = FALSE, use_lasso = FALSE),
    sim_args = list(rand_ratio = 1, draw_treatment = TRUE, max_entry = 6, cens_adjust = 20),
    analysis_time = 30, cens_adjust = 20, parallel_args = list(plan = "sequential"),
    details = FALSE, verbose_n_sims = 0, seed = 5000L, record = "full")
  est <- as.data.frame(attr(demo, "estimators"))
  g_est <- function(k) { y <- est[est$estimator == k, ]; y[match(demo$sim, y$sim), ] }
  itt <- g_est("itt_all"); rgn <- g_est("itt_region"); gh <- g_est("region_Ghat")
  A <- !is.na(itt$hi) & itt$hi < 1; K <- A & !is.na(rgn$hr) & rgn$hr > 1
  s <- mrct_core_summary(demo)
  expect_identical(s$count[s$quantity == "power"], as.numeric(sum(A)))
  expect_identical(s$count[s$quantity == "trigger"], as.numeric(sum(K)))
  expect_gt(sum(K), 0)
  expect_identical(s$count[s$quantity == "declared"], as.numeric(sum(demo$any_found == 1 & K)))
})
