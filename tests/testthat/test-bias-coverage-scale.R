# fs_sim_bias_coverage(scale = ) -- add-only (TASK_continuous_field_mac_2026-09-07).
#   1. the default reproduces the log-scale output exactly (byte-identical
#      to a call written before the argument existed);
#   2. scale = "identity" scores difference measures on the natural scale:
#      bias = mean(est - target), SD/SE untransformed, one-sided bound
#      est -/+ z * SE, and non-positive estimates are NOT dropped;
#   3. the Gaussian reference is evaluated at the same (b, r) formula on
#      either scale, and the output columns are identical.

.bcs_frame <- function(n = 300L, seed = 7L, shift = 0) {
  set.seed(seed)
  tr <- shift + stats::rnorm(n, 0, 20)          # MD-scale targets (any sign)
  e  <- tr + 5 + stats::rnorm(n, 0, 15)          # estimates with bias +5
  se <- rep(15, n)
  data.frame(
    detected = 1L, betaHhat_H = tr, betaHhat_Hc = tr,
    nv_H_est = e, nv_H_lo = e - 1.96 * se, nv_H_hi = e + 1.96 * se, nv_H_se = se,
    mr_H_est = e, mr_H_lo = e - 1.96 * se, mr_H_hi = e + 1.96 * se, mr_H_se_ij = se,
    fld_H_est2 = e, fld_H_lo2s = e - 1.96 * se, fld_H_hi2s = e + 1.96 * se,
    fld_H_se = se, fld_H_lo1s = e - 1.645 * se, fld_H_up1s = e + 1.645 * se)
}

test_that("scale = 'log' (the default) is unchanged and is the default", {
  set.seed(3)
  n <- 200L
  tr <- exp(stats::rnorm(n, 0, 0.2)); e <- tr * exp(stats::rnorm(n, 0, 0.2))
  se <- rep(0.2, n)
  res <- data.frame(
    detected = 1L, betaHhat_H = tr,
    nv_H_est = e, nv_H_lo = e * exp(-1.96 * se), nv_H_hi = e * exp(1.96 * se), nv_H_se = se,
    mr_H_est = e, mr_H_lo = e * exp(-1.96 * se), mr_H_hi = e * exp(1.96 * se), mr_H_se_ij = se,
    fld_H_est2 = e, fld_H_lo2s = e * exp(-1.96 * se), fld_H_hi2s = e * exp(1.96 * se),
    fld_H_se = se, fld_H_lo1s = e * exp(-1.645 * se))
  d <- fs_sim_bias_coverage(res, block = "H")
  l <- fs_sim_bias_coverage(res, block = "H", scale = "log")
  expect_identical(d, l)
  # Hand computation on the log scale for the MR row.
  m <- l[l$estimator == "mr", ]
  expect_equal(m$bias_log, mean(log(e) - log(tr)))
  expect_equal(m$sd_emp, stats::sd(log(e)))
  expect_equal(m$cov1, mean(tr >= exp(log(e) - qnorm(0.95) * se)))
  expect_equal(m$cov1_ref, pnorm(qnorm(0.95) * m$r - m$b))
  expect_error(fs_sim_bias_coverage(res, block = "H", scale = "sqrt"))
})

test_that("scale = 'identity' scores on the natural scale", {
  res <- .bcs_frame()
  tr <- res$betaHhat_H; e <- res$mr_H_est; se <- res$mr_H_se_ij
  id <- fs_sim_bias_coverage(res, block = "H", scale = "identity")
  expect_setequal(id$estimator, c("naive", "mr", "fld"))
  m <- id[id$estimator == "mr", ]
  expect_equal(m$n, nrow(res))
  expect_equal(m$bias_log, mean(e - tr))
  expect_equal(m$sd_emp, stats::sd(e))
  expect_equal(m$se_mean, mean(se))
  expect_equal(m$b, mean(e - tr) / stats::sd(e))
  expect_equal(m$r, mean(se) / stats::sd(e))
  expect_equal(m$cov2, mean(tr >= res$mr_H_lo & tr <= res$mr_H_hi))
  expect_equal(m$cov1, mean(tr >= e - qnorm(0.95) * se))
  expect_equal(m$cov1_ref, pnorm(qnorm(0.95) * m$r - m$b))
  expect_equal(m$cov2_ref, pnorm(qnorm(0.975) * m$r - m$b) - pnorm(-qnorm(0.975) * m$r - m$b))
  f <- id[id$estimator == "fld", ]
  expect_equal(f$cov1, mean(tr >= res$fld_H_lo1s))
  # Upper side: bound est + z * SE, coverage target <= bound, reference at +b.
  up <- fs_sim_bias_coverage(res, block = "H", side = "upper", scale = "identity")
  mu <- up[up$estimator == "mr", ]
  expect_equal(mu$cov1, mean(tr <= e + qnorm(0.95) * se))
  expect_equal(mu$cov1_ref, pnorm(qnorm(0.95) * mu$r + mu$b))
  fu <- up[up$estimator == "fld", ]
  expect_equal(fu$cov1, mean(tr <= res$fld_H_up1s))
  expect_identical(names(id),
                   names(suppressWarnings(fs_sim_bias_coverage(res, block = "H"))))
})

test_that("scale = 'identity' keeps non-positive estimates; 'log' drops them", {
  res <- .bcs_frame(shift = -30)                 # most estimates negative
  expect_true(mean(res$mr_H_est <= 0) > 0.5)
  id <- fs_sim_bias_coverage(res, block = "H", scale = "identity")
  m  <- id[id$estimator == "mr", ]
  expect_equal(m$n, nrow(res))
  expect_true(is.finite(m$bias_log) && is.finite(m$cov1))
  # The log path only sees the positive pairs, so its n is strictly smaller.
  lg <- suppressWarnings(fs_sim_bias_coverage(res, block = "H", scale = "log"))
  ml <- lg[lg$estimator == "mr", ]
  expect_lt(ml$n, m$n)
})
