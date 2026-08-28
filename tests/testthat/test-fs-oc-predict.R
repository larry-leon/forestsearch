# Tests for fs_oc_family_enumerate() and fs_oc_predict().
#
# Small and fast: a 2000-row synthetic glm_dgm with two continuous covariates
# and one binary factor, and a hand-built M = 3 family.  Draw counts are kept
# to a few thousand; nothing here depends on Monte-Carlo precision beyond
# what a seed pins.
#
# The fidelity of fs_oc_predict() to the prediction document's chunk is NOT
# tested here (it needs the tracked payload); see
# dev/glm-continuous-sims/fidelity_fs_oc_predict_2026-08-28.R.

# ---- fixtures ---------------------------------------------------------------

make_oc_dgm <- function(N = 2000L, seed = 1L) {
  set.seed(seed)
  age <- round(stats::rnorm(N, 35, 9))
  pre <- round(stats::rexp(N, 1 / 500))
  V   <- factor(stats::rbinom(N, 1L, 0.42), levels = 0:1)
  inQ <- as.integer(age > 34 & pre <= 745)
  mu0 <- 40 + 0.2 * age
  structure(list(
    df_super = data.frame(age = age, preanti = pre, V = V,
                          mu0 = mu0, mu1 = mu0 - 26 - 14 * inQ,
                          flag_harm = inQ),
    outcome_type = "continuous", effect_measure = "MD",
    model_params = list(sigma = 127.5)),
    class = c("glm_dgm", "list"))
}

oc_args <- list(confounders.name = c("age", "preanti", "V"),
                conf.cont_jcuts = list(age = 4, preanti = 4),
                n.min = 60, effect.threshold = 30, consistency.threshold = 10)

make_hand_family <- function(piQ = 0.34) {
  fam <- structure(list(
    lab = c("Q", "P", "D"), Pg = c(piQ, 0.45, 0.31),
    PQg = c(1, 0.28 / 0.45, 1),
    sens_g = c(1, 0.28 / piQ, 0.31 / piQ),
    spec_g = c(1, 1 - 0.17 / (1 - piQ), 1),
    ovl = matrix(c(piQ, 0.28, 0.31,
                   0.28, 0.45, 0.28 * 0.31 / piQ,
                   0.31, 0.28 * 0.31 / piQ, 0.31), 3, 3),
    M = 3L, PQ = piQ), class = c("fs_oc_family", "list"))
  fam$beta_g <- 26 + 14 * fam$PQg
  fam$se_g   <- 13.7 * sqrt(2) * sqrt(piQ / fam$Pg)
  fam
}

oc_quantities <- c("det_rate", "EnH", "Esens", "Espec", "Eppv", "Enpv",
                   "EbetaH", "Enaive_bias", "mass_below")


# ---- 1. interface contract --------------------------------------------------

test_that("fs_oc_predict runs a hand-built M = 3 family and returns finite quantities", {
  fam <- make_hand_family()
  p <- fs_oc_predict(family = fam, n = 500, c1 = 30, c2 = 10,
                     consistency_method = "split", draws = 5000, seed = 7)
  expect_s3_class(p, "fs_oc_predict")
  for (nm in oc_quantities) {
    expect_true(is.numeric(p[[nm]]) && length(p[[nm]]) == 1L && is.finite(p[[nm]]),
                info = nm)
  }
  expect_true(p$det_rate >= 0 && p$det_rate <= 1)
  expect_equal(sum(p$sel_c), 1)
  expect_length(p$P1, 3L)
  expect_length(p$p_sel, 3L)
  expect_true(all(p$P1 >= 0 & p$P1 <= 1))
  expect_equal(p$det_rate_se, sqrt(p$det_rate * (1 - p$det_rate) / 5000))
  expect_equal(p$M, 3L)
  expect_identical(p$settings$consistency_method, "split")
  expect_output(print(p), "Predicted operating characteristics")
})


# ---- 2. family construction -------------------------------------------------

test_that("fs_oc_family_enumerate reports population proportions coherently", {
  dgm <- make_oc_dgm()
  fam <- fs_oc_family_enumerate(dgm, oc_args, n = 500)
  expect_s3_class(fam, "fs_oc_family")
  expect_equal(fam$M, length(fam$lab))
  expect_equal(fam$M, ncol(fam$memb))
  expect_equal(fam$Pg, unname(colMeans(fam$memb)))
  expect_true(isSymmetric(fam$ovl))
  expect_equal(unname(diag(fam$ovl)), fam$Pg)
  expect_true(all(fam$PQg >= 0 & fam$PQg <= 1))
  expect_true(all(fam$sens_g >= 0 & fam$sens_g <= 1))
  expect_true(all(fam$spec_g >= 0 & fam$spec_g <= 1))
  expect_true(all(fam$Pg >= 60 / 500))
  expect_equal(fam$PQ, mean(dgm$df_super$flag_harm))
  # direct recomputation of purity and sensitivity from the membership
  inQ <- dgm$df_super$flag_harm == 1
  expect_equal(fam$PQg, unname(colMeans(fam$memb & inQ)) / fam$Pg)
  expect_equal(fam$sens_g, unname(colMeans(fam$memb & inQ)) / mean(inQ))
  # no two candidates share a membership vector
  keys <- apply(fam$memb, 2L, function(v) paste(which(v), collapse = ","))
  expect_false(any(duplicated(keys)))
  expect_equal(fam$counts[["M"]], fam$M)
  expect_output(print(fam), "Population-enumerated candidate family")
})

test_that("the scale table drives beta_g and se_g", {
  dgm <- make_oc_dgm()
  fam <- fs_oc_family_enumerate(dgm, oc_args, n = 500)
  sc  <- fs_dgm_scale(dgm)
  reg <- sc$regions
  iQ <- match("Q", reg$region); iQc <- match("Qc", reg$region)
  tauQc <- abs(reg$m_tau[iQc]); bint <- abs(reg$m_tau[iQ]) - tauQc
  seQ1000 <- sqrt(reg$V_eff[iQ] / (1000 * reg$P_g[iQ]))
  expect_equal(fam$beta_g, tauQc + bint * fam$PQg)
  expect_equal(fam$se_g, seQ1000 * sqrt(2) * sqrt(reg$P_g[iQ] / fam$Pg))
})


# ---- 3. conf.cont_jcuts is honoured ------------------------------------------

test_that("conf.cont_jcuts sets the number of cuts and places them at population quantiles", {
  dgm <- make_oc_dgm()
  args10 <- utils::modifyList(oc_args, list(conf.cont_jcuts = list(age = 10, preanti = 4)))
  args4  <- oc_args
  fam10 <- fs_oc_family_enumerate(dgm, args10, n = 500)
  fam4  <- fs_oc_family_enumerate(dgm, args4,  n = 500)
  age10 <- fam10$cuts[grepl("^age", fam10$cuts)]
  age4  <- fam4$cuts[grepl("^age", fam4$cuts)]
  expect_length(age10, 10L)
  expect_length(age4, 4L)
  # cut values equal the k/(J+1) population quantiles of df_super$age
  vals10 <- as.numeric(sub("^age <= ", "", age10))
  vals4  <- as.numeric(sub("^age <= ", "", age4))
  expect_equal(vals10, unname(round(stats::quantile(dgm$df_super$age, (1:10) / 11), 1)))
  expect_equal(vals4,  unname(round(stats::quantile(dgm$df_super$age, (1:4) / 5), 1)))
})


# ---- 4. determinism ---------------------------------------------------------

test_that("the family is deterministic and the prediction is seed-reproducible", {
  dgm <- make_oc_dgm()
  f1 <- fs_oc_family_enumerate(dgm, oc_args, n = 500)
  f2 <- fs_oc_family_enumerate(dgm, oc_args, n = 500)
  expect_identical(f1, f2)
  p1 <- fs_oc_predict(family = f1, n = 500, c1 = 30, c2 = 10,
                      consistency_method = "split", draws = 4000, seed = 11)
  p2 <- fs_oc_predict(family = f1, n = 500, c1 = 30, c2 = 10,
                      consistency_method = "split", draws = 4000, seed = 11)
  p3 <- fs_oc_predict(family = f1, n = 500, c1 = 30, c2 = 10,
                      consistency_method = "split", draws = 4000, seed = 12)
  expect_identical(p1[oc_quantities], p2[oc_quantities])
  expect_identical(p1$P1, p2$P1)
  expect_false(identical(p1[oc_quantities], p3[oc_quantities]))
})


# ---- 5. override semantics ---------------------------------------------------

test_that("explicit c1/c2 override forestsearch_args and det_rate is monotone in c1", {
  fam <- make_hand_family()
  args <- list(effect.threshold = 30, consistency.threshold = 10)
  p_def <- fs_oc_predict(forestsearch_args = args, family = fam, n = 500,
                         consistency_method = "split", draws = 4000, seed = 3)
  expect_equal(p_def$settings$c1, 30)
  expect_equal(p_def$settings$c2, 10)
  p_ovr <- fs_oc_predict(forestsearch_args = args, family = fam, n = 500,
                         c1 = 45, c2 = 20, consistency_method = "split",
                         draws = 4000, seed = 3)
  expect_equal(p_ovr$settings$c1, 45)
  expect_equal(p_ovr$settings$c2, 20)
  expect_false(identical(p_def$det_rate, p_ovr$det_rate))
  # same draws (same seed) across a c1 ladder: declaration can only shrink
  ladder <- c(20, 30, 40, 50, 60)
  rates <- vapply(ladder, function(cc)
    fs_oc_predict(family = fam, n = 500, c1 = cc, c2 = 10,
                  consistency_method = "split", draws = 4000, seed = 3)$det_rate,
    numeric(1))
  expect_true(all(diff(rates) <= 0))
  # n override: EnH scales with n at fixed draws
  p700 <- fs_oc_predict(family = fam, n = 700, c1 = 30, c2 = 10,
                        consistency_method = "split", draws = 4000, seed = 3)
  expect_equal(p700$EnH / p_def$EnH, 700 / 500)
  # missing thresholds are an error
  expect_error(fs_oc_predict(family = fam, n = 500, consistency_method = "split"),
               "c1")
})


# ---- 6. size guard ----------------------------------------------------------

test_that("max_M below the realized M raises the size guard", {
  dgm <- make_oc_dgm()
  fam <- fs_oc_family_enumerate(dgm, oc_args, n = 500)
  expect_gt(fam$M, 5L)
  expect_error(fs_oc_family_enumerate(dgm, oc_args, n = 500, max_M = 5L),
               "above max_M = 5")
  expect_error(fs_oc_predict(dgm, oc_args, n = 500, consistency_method = "split",
                             draws = 100, max_M = 5L),
               "above max_M = 5")
})


# ---- 7. the resample gate is refused, not approximated ----------------------

test_that("consistency_method = 'resample' stops with the sigma_D explanation", {
  fam <- make_hand_family()
  expect_error(fs_oc_predict(family = fam, n = 500, c1 = 30, c2 = 10),
               "sigma_D")
  expect_error(fs_oc_predict(family = fam, n = 500, c1 = 30, c2 = 10,
                             consistency_method = "resample"),
               "not implemented")
})


# ---- floors -----------------------------------------------------------------

test_that("n.min = NULL resolves to max(60, ceiling(n.min.frac * n)) and the size floor moves", {
  dgm <- make_oc_dgm()
  args_null <- utils::modifyList(oc_args, list(n.min = NULL, n.min.frac = 0.2))
  args_null["n.min"] <- list(NULL)
  fam <- fs_oc_family_enumerate(dgm, args_null, n = 500)
  expect_equal(fam$args_used$n.min, 100L)
  expect_true(all(fam$Pg >= 100 / 500))
  fam60 <- fs_oc_family_enumerate(dgm, oc_args, n = 500)
  expect_gt(fam60$M, fam$M)
})

test_that("a malformed family is rejected", {
  fam <- make_hand_family()
  fam$ovl <- fam$ovl[1:2, 1:2]
  expect_error(fs_oc_predict(family = fam, n = 500, c1 = 30, c2 = 10,
                             consistency_method = "split", draws = 10),
               "M x M")
  expect_error(fs_oc_predict(family = list(), n = 500, c1 = 30, c2 = 10,
                             consistency_method = "split", draws = 10),
               "fs_oc_family")
})
