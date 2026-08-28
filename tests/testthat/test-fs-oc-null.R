# Tests for the null branch of fs_oc_family_enumerate() and its propagation
# through fs_oc_predict() / fs_oc_grid() / fs_oc_invert().

.null_dgm <- function(N = 2000L, seed = 1L, model = "null") {
  set.seed(seed)
  age <- round(stats::rnorm(N, 35, 9)); pre <- round(stats::rexp(N, 1 / 500))
  V   <- factor(stats::rbinom(N, 1L, 0.42), levels = 0:1)
  mu0 <- 40 + 0.2 * age
  d <- list(df_super = data.frame(age = age, preanti = pre, V = V, mu0 = mu0,
                                  mu1 = mu0 - 26, flag_harm = 0L),
            outcome_type = "continuous", effect_measure = "MD",
            model_params = list(sigma = 127.5))
  if (!is.null(model)) d$model <- model
  structure(d, class = c("glm_dgm", "list"))
}
.alt_dgm <- function(N = 2000L, seed = 1L) {
  set.seed(seed)
  age <- round(stats::rnorm(N, 35, 9)); pre <- round(stats::rexp(N, 1 / 500))
  V   <- factor(stats::rbinom(N, 1L, 0.42), levels = 0:1)
  inQ <- as.integer(age > 34 & pre <= 745); mu0 <- 40 + 0.2 * age
  structure(list(df_super = data.frame(age = age, preanti = pre, V = V, mu0 = mu0,
                                       mu1 = mu0 - 26 - 14 * inQ, flag_harm = inQ),
                 outcome_type = "continuous", effect_measure = "MD",
                 model_params = list(sigma = 127.5), model = "alt"),
            class = c("glm_dgm", "list"))
}
.null_args <- list(confounders.name = c("age", "preanti", "V"),
                   conf.cont_jcuts = list(age = 4, preanti = 4), n.min = 60,
                   effect.threshold = 30, consistency.threshold = 10)


test_that("the null branch sets the family fields as specified", {
  dgm <- .null_dgm()
  fam <- fs_oc_family_enumerate(dgm, .null_args, n = 500)
  expect_true(isTRUE(fam$null))
  expect_true(all(fam$PQg == 0))
  expect_true(all(is.na(fam$sens_g)))
  expect_equal(fam$spec_g, 1 - fam$Pg)
  expect_equal(fam$PQ, 0)
  expect_equal(length(unique(fam$beta_g)), 1L)
  expect_equal(fam$beta_g[1], 26)                     # |effect_Qc|, oriented
  # se_g from the whole-population effective variance at (n, Pg)
  sc <- fs_dgm_scale(dgm, regions = list(S = rep(TRUE, nrow(dgm$df_super))))
  expect_equal(fam$se_g, sqrt(sc$regions$V_eff / 1000) * sqrt(1000 / 500) * sqrt(1 / fam$Pg))
  expect_equal(fam$scale$regions$region, "S")
  # enumeration and floors are untouched by the branch: same M as the
  # alternative built on the same covariates
  fam_alt <- fs_oc_family_enumerate(.alt_dgm(), .null_args, n = 500)
  expect_equal(fam$M, fam_alt$M)
  expect_identical(fam$lab, fam_alt$lab)
  expect_identical(fam$ovl, fam_alt$ovl)
  expect_output(print(fam), "NULL DGM")
  # a null DGM with no model field is still detected from the flag
  fam2 <- fs_oc_family_enumerate(.null_dgm(model = NULL), .null_args, n = 500)
  expect_true(isTRUE(fam2$null))
  expect_equal(fam2$se_g, fam$se_g)
})

test_that("an inconsistent model/flag pair is refused", {
  bad <- .null_dgm(); bad$df_super$flag_harm[1:50] <- 1L
  expect_error(fs_oc_family_enumerate(bad, .null_args, n = 500), "inconsistent")
  bad2 <- .alt_dgm(); bad2$df_super$flag_harm[] <- 0L
  expect_error(fs_oc_family_enumerate(bad2, .null_args, n = 500), "inconsistent")
})

test_that("NA sensitivity propagates only to E[sens]; declaration is monotone in c1", {
  fam <- fs_oc_family_enumerate(.null_dgm(), .null_args, n = 500)
  for (g in c("resample", "split")) {
    p <- fs_oc_predict(family = fam, n = 500, c1 = 30, c2 = 10,
                       consistency_method = g, draws = 5000, seed = 3)
    expect_true(is.na(p$Esens))
    for (nm in c("det_rate", "EnH", "Espec", "Eppv", "Enpv", "EbetaH",
                 "Enaive_bias", "mass_below")) {
      expect_true(is.finite(p[[nm]]), info = paste(g, nm))
    }
    expect_equal(p$Eppv, 0)
    expect_equal(p$Enpv, 1)
    expect_equal(p$EbetaH, 26)
    expect_equal(p$mass_below, 1)                  # every true mean (26) < c1 = 30
    expect_equal(sum(p$sel_c), 1)
  }
  G <- fs_oc_grid(family = fam, n = 500, c1 = c(10, 20, 30, 45, 60, 80),
                  c2 = 10, consistency_method = c("resample", "split"),
                  draws = 5000, seed = 3)
  for (g in c("resample", "split")) {
    r <- G$table[G$table$consistency_method == g, ]
    expect_true(all(diff(r$det_rate[order(r$c1)]) <= 0), info = g)
    expect_true(all(is.na(r$Esens)))
  }
  iv <- fs_oc_invert(family = fam, n = 500, target = c(0.05, 0.10), solve_for = "c1",
                     c2 = 10, consistency_method = "resample", draws = 5000, seed = 3)
  expect_s3_class(iv, "fs_oc_invert_list")
  tb <- attr(iv, "table")
  expect_true(all(tb$attainable))
  expect_true(all(tb$achieved >= tb$target))
  expect_true(tb$value[1] > tb$value[2])          # smaller type-I error needs a higher c1
})

test_that("the alternative path is untouched by the null branch", {
  ref_path <- testthat::test_path("..", "..", "dev", "glm-continuous-sims",
                                  "prerefactor_reference_2026-08-29.rds")
  skip_if_not(file.exists(ref_path), "pre-refactor reference not present (dev/ not shipped)")
  ref <- readRDS(ref_path)
  fam <- fs_oc_family_enumerate(.alt_dgm(), .null_args, n = 500)
  expect_false(isTRUE(fam$null))
  strip <- function(p) { p$family <- NULL; p }
  expect_identical(strip(fs_oc_predict(family = fam, n = 500, c1 = 30, c2 = 10,
                                       consistency_method = "resample", draws = 2e4, seed = 20260829)),
                   ref$fam_resample)
  expect_identical(strip(fs_oc_predict(family = fam, n = 500, c1 = 30, c2 = 10,
                                       consistency_method = "split", draws = 2e4, seed = 20260829)),
                   ref$fam_split)
})
