# Tests for fs_mr_oc_summary().
#
# Deterministic and fast: synthetic payloads with hand-computable answers, plus
# every error path.  Parity against the batch documents' inline machinery on the
# committed payloads is checked separately by
# dev/glm-continuous-sims/parity_fs_mr_oc_summary.R, which needs those payloads
# and so cannot live in the suite.

skip_on_cran()

make_payload <- function(n_sims = 10, detected = rep(TRUE, n_sims),
                         adverse = FALSE, nb_boots = NULL) {
  s <- seq_len(n_sims)
  r <- data.frame(
    sim_id    = s,
    status    = ifelse(detected, "DETECTED", "NO-DETECTION"),
    detected  = as.integer(detected),
    n_harm    = rep(70L, n_sims),
    n_true    = rep(241L, n_sims),       # true |H|; |Hc|/N = 459/700
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
  # non-detection contributes npv = |Hc|/N = 1 - 241/700 = 459/700
  expect_equal(id$npv[id$convention == "unconditional"], (0.7 + 459/700) / 2)
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


# ---- NPV for a replicate with no identified subgroup ------------------------
# Leon et al. (2024), Section 3: with no subgroup identified the estimated
# complement is everyone, so spec = sens(Hc) = 1 and npv = ppv(Hc) = |Hc|/N,
# the share of that replicate's subjects truly in the complement (1 under
# the null).  Hand-built payloads; no search, no simulation.

.nofind_payload <- function(null = FALSE) {
  pl <- make_payload(n_sims = 4, detected = c(TRUE, TRUE, FALSE, FALSE))
  pl$meta$n_sample <- 200L
  r <- pl$results
  if (null) {
    r$n_true <- rep(0L, 4)
    r$sens   <- NA_real_                    # 0/0 under the null
    r$spec   <- c(0.85, 0.95, NA, NA)
    r$ppv    <- c(0, 0, NA, NA)
    r$npv    <- c(1, 1, NA, NA)             # found-case npv when Q is empty
  } else {
    r$n_true <- c(40L, 45L, 50L, 60L)       # |Hc|/N = 0.75, 0.70 for 3, 4
    r$sens   <- c(0.6, 0.5, NA, NA)
    r$spec   <- c(0.85, 0.95, NA, NA)
    r$ppv    <- c(0.5, 0.7, NA, NA)
    r$npv    <- c(0.9, 0.8, NA, NA)         # the MR harness stores NA
  }
  pl$results <- r
  pl
}

test_that("no-find, alternative: a non-detected replicate contributes npv = |Hc|/N and spec = 1", {
  pl <- .nofind_payload()
  id <- fs_mr_oc_summary(pl)$identification
  u  <- id[id$convention == "unconditional", ]
  expect_equal(u$npv,  (0.9 + 0.8 + 150/200 + 140/200) / 4)
  expect_equal(u$spec, (0.85 + 0.95 + 1 + 1) / 4)
  expect_equal(u$sens, (0.6 + 0.5 + 0 + 0) / 4)
  expect_equal(u$ppv,  (0.5 + 0.7 + 0 + 0) / 4)

  # one non-detected replicate on its own: exactly |Hc|/N and spec 1
  one <- pl; one$results <- pl$results[3, , drop = FALSE]
  u1 <- fs_mr_oc_summary(one)$identification
  u1 <- u1[u1$convention == "unconditional", ]
  expect_equal(u1$npv,  150 / 200)
  expect_equal(u1$spec, 1)

  # the detected-only summary is unchanged
  c0 <- id[id$convention == "conditional", ]
  expect_equal(c0$npv,  (0.9 + 0.8) / 2)
  expect_equal(c0$spec, (0.85 + 0.95) / 2)
})

test_that("no-find, null: the same input gives npv = 1 and spec = 1 for a non-detected replicate", {
  pl <- .nofind_payload(null = TRUE)
  id <- fs_mr_oc_summary(pl)$identification
  u  <- id[id$convention == "unconditional", ]
  expect_equal(u$npv,  1)
  expect_equal(u$spec, (0.85 + 0.95 + 1 + 1) / 4)

  one <- pl; one$results <- pl$results[3:4, , drop = FALSE]
  u1 <- fs_mr_oc_summary(one)$identification
  u1 <- u1[u1$convention == "unconditional", ]
  expect_equal(u1$npv,  1)
  expect_equal(u1$spec, 1)
})

test_that("no-find: unconditional npv is NA when |Hc|/N is not available", {
  pl <- .nofind_payload(); pl$results$n_true <- NULL
  id <- fs_mr_oc_summary(pl)$identification
  expect_true(is.na(id$npv[id$convention == "unconditional"]))
  expect_equal(id$npv[id$convention == "conditional"], (0.9 + 0.8) / 2)
  expect_equal(id$spec[id$convention == "unconditional"], (0.85 + 0.95 + 1 + 1) / 4)

  pl2 <- .nofind_payload(); pl2$meta$n_sample <- NULL
  id2 <- fs_mr_oc_summary(pl2)$identification
  expect_true(is.na(id2$npv[id2$convention == "unconditional"]))

  pl3 <- .nofind_payload(); pl3$results$n_true[4] <- NA_integer_
  id3 <- fs_mr_oc_summary(pl3)$identification
  expect_true(is.na(id3$npv[id3$convention == "unconditional"]))
})

test_that("no-find: unconditional npv agrees with the run_simulation_analysis() extractors", {
  # Four replicates of N = 200, scored by .extract_fs_estimates_gen(): two
  # with an identified subgroup (treat.recommend == 0), two with none.
  N <- 200L
  mk <- function(nH) data.frame(
    id = seq_len(N), y = seq_len(N) / N, ev = rep(c(1L, 1L, 0L), length.out = N),
    tr = rep(c(0L, 1L), length.out = N),
    flag_harm = c(rep(1L, nH), rep(0L, N - nH)))
  score <- function(df, rec = NULL) {
    fs_full <- if (is.null(rec)) NULL else list(df.est = data.frame(treat.recommend = rec))
    .extract_fs_estimates_gen(
      df = df, fs_res = NULL, fs_full = fs_full, dgm = NULL,
      cox_formula = NULL, cox_formula_adj = NULL,
      outcome_name = "y", event_name = "ev", treat_name = "tr",
      harm_col = "flag_harm", analysis = "FS", verbose = FALSE)
  }
  nH  <- c(40L, 45L, 50L, 60L)
  rec <- list(c(rep(0L, 30), rep(1L, N - 30)),            # 30 of 40 caught
              c(rep(1L, 20), rep(0L, 60), rep(1L, N - 80)), # 25 of 45, 35 false
              NULL, NULL)
  sc <- lapply(1:4, function(k) suppressWarnings(score(mk(nH[k]), rec[[k]])))
  ext_npv  <- vapply(sc, function(o) o$npv,  numeric(1))
  ext_spec <- vapply(sc, function(o) o$spec, numeric(1))
  expect_identical(vapply(sc, function(o) o$any.H, integer(1)), c(1L, 1L, 0L, 0L))

  # The MR harness's record of the same replicates: the found-case values for
  # detected ones, NA for the non-detected (its .classify() returns NA).
  pl <- make_payload(n_sims = 4, detected = c(TRUE, TRUE, FALSE, FALSE))
  pl$meta$n_sample <- N
  pl$results$n_true <- nH
  pl$results$npv  <- c(ext_npv[1:2], NA, NA)
  pl$results$spec <- c(ext_spec[1:2], NA, NA)
  id <- fs_mr_oc_summary(pl)$identification
  u  <- id[id$convention == "unconditional", ]
  expect_equal(u$npv,  mean(ext_npv))
  expect_equal(u$spec, mean(ext_spec))
})

test_that("no-find: a payload without a status column keeps its unconditional npv", {
  # No status column: no replicate is marked, nothing is filled, and the
  # stored values are averaged, as before.
  pl <- .nofind_payload(); pl$results$status <- NULL
  id <- fs_mr_oc_summary(pl)$identification
  expect_equal(id$npv[id$convention == "unconditional"], (0.9 + 0.8) / 2)

  # a non-numeric n_sample gives NA, not an error
  pl2 <- .nofind_payload(); pl2$meta$n_sample <- "200"
  id2 <- fs_mr_oc_summary(pl2)$identification
  expect_true(is.na(id2$npv[id2$convention == "unconditional"]))
})
