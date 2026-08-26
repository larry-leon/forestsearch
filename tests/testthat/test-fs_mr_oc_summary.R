# Tests for fs_mr_oc_summary().
#
# Deterministic and fast: synthetic payloads with hand-computable answers, plus
# every error path.  Parity against the batch documents' inline machinery on the
# committed payloads is checked separately by
# dev/glm-continuous-sims/parity_fs_mr_oc_summary.R, which needs those payloads
# and so cannot live in the suite.

make_payload <- function(n_sims = 10, detected = rep(TRUE, n_sims),
                         adverse = FALSE, nb_boots = NULL) {
  s <- seq_len(n_sims)
  r <- data.frame(
    sim_id    = s,
    status    = ifelse(detected, "DETECTED", "NO-DETECTION"),
    detected  = as.integer(detected),
    n_harm    = rep(70L, n_sims),
    betaHhat_H  = rep(-30, n_sims),      # RAW scale
    betaHhat_Hc = rep(-20, n_sims),
    or_H_est = rep(40, n_sims), or_H_lo = rep(30, n_sims),
    or_H_hi  = rep(50, n_sims), or_H_se = rep(5, n_sims),
    nv_H_est = rep(60, n_sims), nv_H_lo = rep(50, n_sims),
    nv_H_hi  = rep(70, n_sims), nv_H_se = rep(5, n_sims),
    mr_H_est = rep(35, n_sims), mr_H_lo = rep(25, n_sims),
    mr_H_hi  = rep(45, n_sims), mr_H_se_ij = rep(5, n_sims),
    or_Hc_est = rep(26, n_sims), or_Hc_lo = rep(20, n_sims),
    or_Hc_hi  = rep(32, n_sims), or_Hc_se = rep(3, n_sims),
    nv_Hc_est = rep(22, n_sims), nv_Hc_lo = rep(16, n_sims),
    nv_Hc_hi  = rep(28, n_sims), nv_Hc_se = rep(3, n_sims),
    mr_Hc_est = rep(24, n_sims), mr_Hc_lo = rep(18, n_sims),
    mr_Hc_hi  = rep(30, n_sims), mr_Hc_se_ij = rep(3, n_sims),
    sens = rep(0.4, n_sims), spec = rep(0.8, n_sims),
    ppv  = rep(0.5, n_sims), npv  = rep(0.7, n_sims),
    stringsAsFactors = FALSE
  )
  list(results = r,
       truth = list(effect_Q = -40, effect_Qc = -26, beta_inter = -14,
                    prevalence_Q = 0.3446, marg_H = -36, marg_Hc = -22),
       meta = list(n_sample = 700L, n_sims = n_sims, adverse_outcome = adverse,
                   nb_boots = nb_boots, mr_draws = 5000L,
                   sg_focus = "maxeffCons", pkg_version = "test"))
}

test_that("orientation bridges raw betaHhat to the oriented estimator scale", {
  oc <- fs_mr_oc_summary(make_payload())
  expect_equal(oc$targets$orient, -1)
  expect_equal(oc$targets$struct_H, 40)      # -1 * -40
  expect_equal(oc$targets$marg_H, 36)
  H <- oc$estimation[oc$estimation$block == "H", ]
  # oriented beta(Hhat) = 30; naive est 60 => bias 30
  expect_equal(H$bias_beta[H$estimator == "naive"], 30)
  expect_equal(H$bias_struct[H$estimator == "naive"], 20)   # 60 - 40
  expect_equal(H$bias_marg[H$estimator == "naive"], 24)     # 60 - 36
  expect_equal(H$bias_oracle[H$estimator == "naive"], 20)   # 60 - 40
})

test_that("adverse_outcome flips the orientation", {
  oc <- fs_mr_oc_summary(make_payload(adverse = TRUE))
  expect_equal(oc$targets$orient, 1)
  expect_equal(oc$targets$struct_H, -40)
})

test_that("coverage is computed against each of the four targets", {
  H <- fs_mr_oc_summary(make_payload())$estimation
  H <- H[H$block == "H", ]
  mr <- H[H$estimator == "MR", ]
  expect_equal(mr$cov_beta, 1)     # 30 in [25, 45]
  expect_equal(mr$cov_struct, 1)   # 40 in [25, 45]
  expect_equal(mr$cov_marg, 1)     # 36 in [25, 45]
  nv <- H[H$estimator == "naive", ]
  expect_equal(nv$cov_beta, 0)     # 30 not in [50, 70]
})

test_that("width and se_hat are means of the stored quantities", {
  H <- fs_mr_oc_summary(make_payload())$estimation
  H <- H[H$block == "H", ]
  expect_equal(unique(H$width), 20)
  expect_equal(unique(H$se_hat), 5)
})

test_that("MR takes its standard error from the _se_ij column", {
  pl <- make_payload()
  pl$results$mr_H_se_ij <- 9
  H <- fs_mr_oc_summary(pl)$estimation
  expect_equal(H$se_hat[H$block == "H" & H$estimator == "MR"], 9)
})

test_that("FB is dropped unless the run carried a full bootstrap", {
  expect_false("FB" %in% fs_mr_oc_summary(make_payload())$estimation$estimator)
  pl <- make_payload(nb_boots = 200L)
  pl$results$fb_H_est <- 33; pl$results$fb_H_lo <- 20
  pl$results$fb_H_hi  <- 46; pl$results$fb_H_se <- 6
  expect_true("FB" %in% fs_mr_oc_summary(pl)$estimation$estimator)
})

test_that("the two identification conventions agree at full detection", {
  id <- fs_mr_oc_summary(make_payload())$identification
  expect_equal(id$sens[1], id$sens[2])
  expect_equal(id$detection[1], 1)
})

test_that("the two conventions diverge when detection is incomplete", {
  det <- c(rep(TRUE, 5), rep(FALSE, 5))
  id  <- fs_mr_oc_summary(make_payload(detected = det))$identification
  expect_equal(id$detection[1], 0.5)
  expect_equal(id$sens[id$convention == "conditional"], 0.4)
  expect_equal(id$sens[id$convention == "unconditional"], 0.2)   # 0.4 * 0.5
  expect_equal(id$spec[id$convention == "unconditional"], 0.9)   # (.8+1)/2
  expect_equal(id$npv[id$convention == "unconditional"], 0.85)
  expect_equal(id$n_used, c(5, 10))
})

test_that("only detected replicates enter the estimation table", {
  det <- c(rep(TRUE, 5), rep(FALSE, 5))
  e <- fs_mr_oc_summary(make_payload(detected = det))$estimation
  expect_true(all(e$n == 5))
})

test_that("blocks and estimators can be selected", {
  e <- fs_mr_oc_summary(make_payload(), blocks = "H")$estimation
  expect_identical(unique(e$block), "H")
  e2 <- fs_mr_oc_summary(make_payload(), estimators = c("oracle", "MR"))$estimation
  expect_identical(sort(unique(e2$estimator)), c("MR", "oracle"))
})

test_that("digits rounds output without affecting the computation", {
  a <- fs_mr_oc_summary(make_payload())$estimation
  b <- fs_mr_oc_summary(make_payload(), digits = 1)$estimation
  expect_equal(round(a$avg, 1), b$avg)
})

test_that("malformed input is rejected", {
  expect_error(fs_mr_oc_summary(list(results = 1, truth = 1)), "meta")
  expect_error(fs_mr_oc_summary("no-such-file.rds"), "no such payload")
  expect_error(fs_mr_oc_summary(make_payload(), estimators = "nope"),
               "unknown or unavailable")
})

test_that("print returns its argument invisibly", {
  oc <- fs_mr_oc_summary(make_payload())
  expect_output(out <- print(oc), "operating-characteristics")
  expect_identical(out, oc)
})
