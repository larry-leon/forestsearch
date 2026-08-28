# Tests for fs_oc_grid() and fs_oc_invert(), and the refactor guard on
# fs_oc_predict().  Small families, a few thousand draws.

# ---- fixtures ---------------------------------------------------------------

.grid_hand_family <- function(piQ = 0.34, n = 500) {
  fam <- structure(list(
    lab = c("Q", "P", "D"), Pg = c(piQ, 0.45, 0.31),
    PQg = c(1, 0.28 / 0.45, 1),
    sens_g = c(1, 0.28 / piQ, 0.31 / piQ),
    spec_g = c(1, 1 - 0.17 / (1 - piQ), 1),
    ovl = matrix(c(piQ, 0.28, 0.31,
                   0.28, 0.45, 0.28 * 0.31 / piQ,
                   0.31, 0.28 * 0.31 / piQ, 0.31), 3, 3),
    M = 3L, PQ = piQ, n = n), class = c("fs_oc_family", "list"))
  fam$beta_g <- 26 + 14 * fam$PQg
  fam$se_g   <- 13.7 * sqrt(2) * sqrt(piQ / fam$Pg)
  fam
}

.grid_dgm <- function(N = 2000L, seed = 1L) {
  set.seed(seed)
  age <- round(stats::rnorm(N, 35, 9)); pre <- round(stats::rexp(N, 1 / 500))
  V   <- factor(stats::rbinom(N, 1L, 0.42), levels = 0:1)
  inQ <- as.integer(age > 34 & pre <= 745); mu0 <- 40 + 0.2 * age
  structure(list(df_super = data.frame(age = age, preanti = pre, V = V, mu0 = mu0,
                                       mu1 = mu0 - 26 - 14 * inQ, flag_harm = inQ),
                 outcome_type = "continuous", effect_measure = "MD",
                 model_params = list(sigma = 127.5)), class = c("glm_dgm", "list"))
}
.grid_args <- list(confounders.name = c("age", "preanti", "V"),
                   conf.cont_jcuts = list(age = 4, preanti = 4), n.min = 60,
                   effect.threshold = 30, consistency.threshold = 10)


# ---- 1. refactor guard (GATE) -------------------------------------------------

test_that("fs_oc_predict() is identical to the 0.2.4 pre-refactor reference on both gates", {
  ref_path <- testthat::test_path("..", "..", "dev", "glm-continuous-sims",
                                  "prerefactor_reference_2026-08-29.rds")
  skip_if_not(file.exists(ref_path), "pre-refactor reference not present (dev/ not shipped)")
  ref <- readRDS(ref_path)
  hand <- .grid_hand_family(); hand$n <- NULL          # the reference's hand family had no n
  fam  <- fs_oc_family_enumerate(.grid_dgm(), .grid_args, n = 500)
  expect_equal(fam$M, ref$fam_M)
  strip <- function(p) { p$family <- NULL; p }
  expect_identical(strip(fs_oc_predict(family = hand, n = 500, c1 = 30, c2 = 10,
                                       consistency_method = "resample", draws = 2e4, seed = 20260829)),
                   ref$hand_resample)
  expect_identical(strip(fs_oc_predict(family = hand, n = 500, c1 = 30, c2 = 10,
                                       consistency_method = "split", draws = 2e4, seed = 20260829)),
                   ref$hand_split)
  expect_identical(strip(fs_oc_predict(family = fam, n = 500, c1 = 30, c2 = 10,
                                       consistency_method = "resample", draws = 2e4, seed = 20260829)),
                   ref$fam_resample)
  expect_identical(strip(fs_oc_predict(family = fam, n = 500, c1 = 30, c2 = 10,
                                       consistency_method = "split", draws = 2e4, seed = 20260829)),
                   ref$fam_split)
})


# ---- 2. grid / predict identity (GATE) ------------------------------------------

test_that("fs_oc_grid() at one point with block = Inf is identical() to fs_oc_predict()", {
  fam <- .grid_hand_family()
  for (g in c("resample", "split")) {
    p <- fs_oc_predict(family = fam, n = 500, c1 = 30, c2 = 10,
                       consistency_method = g, draws = 8000, seed = 21)
    G <- fs_oc_grid(family = fam, n = 500, c1 = 30, c2 = 10,
                    consistency_method = g, draws = 8000, block = Inf, seed = 21)
    expect_identical(G$results[[1]], p, info = g)
    expect_equal(nrow(G$table), 1L)
    expect_identical(G$table$det_rate, p$det_rate)
  }
})


# ---- 3. blocking invariance ---------------------------------------------------

test_that("blocked draws agree with the one-block path to MC precision (not bit-for-bit)", {
  # The RNG stream is laid out differently across blocks, so the two are
  # different Monte-Carlo estimates of the same quantity: agreement within a
  # few MC SEs holds; identity does not, and is not claimed.
  fam <- .grid_hand_family()
  for (g in c("resample", "split")) {
    G1 <- fs_oc_grid(family = fam, n = 500, c1 = 30, c2 = 10,
                     consistency_method = g, draws = 20000, block = Inf, seed = 4)
    Gb <- fs_oc_grid(family = fam, n = 500, c1 = 30, c2 = 10,
                     consistency_method = g, draws = 20000, block = 4000, seed = 4)
    se <- G1$table$det_rate_se
    expect_lt(abs(G1$table$det_rate - Gb$table$det_rate), 4 * se)
    expect_lt(abs(G1$table$EnH - Gb$table$EnH) / G1$table$EnH, 0.02)
    expect_lt(abs(G1$table$Eppv - Gb$table$Eppv), 0.02)
    expect_equal(sum(Gb$results[[1]]$sel_c), 1)
    expect_equal(Gb$table$block, 4000)
  }
})


# ---- 4. monotonicity ------------------------------------------------------------

test_that("det_rate is non-increasing along c1 and along c2 on one draw set", {
  fam <- .grid_hand_family()
  G <- fs_oc_grid(family = fam, n = 500, c1 = c(10, 20, 30, 40, 50, 60),
                  c2 = c(0, 10, 20, 30), consistency_method = c("resample", "split"),
                  draws = 6000, seed = 8)
  tb <- G$table
  for (g in c("resample", "split")) {
    for (cc in unique(tb$c2)) {
      r <- tb[tb$consistency_method == g & tb$c2 == cc, ]
      expect_true(all(diff(r$det_rate[order(r$c1)]) <= 0), info = paste(g, "c2 =", cc))
    }
    for (cc in unique(tb$c1)) {
      r <- tb[tb$consistency_method == g & tb$c1 == cc, ]
      expect_true(all(diff(r$det_rate[order(r$c2)]) <= 0), info = paste(g, "c1 =", cc))
    }
  }
  expect_equal(nrow(tb), 6 * 4 * 2)
  expect_s3_class(summary(G), "summary.fs_oc_grid")
  expect_output(print(G), "Threshold sweep")
})


# ---- 5. inversion round-trip ----------------------------------------------------

test_that("fs_oc_invert() returns a c1 whose grid rate matches the target", {
  fam <- .grid_hand_family()
  for (g in c("resample", "split")) {
    inv <- fs_oc_invert(family = fam, n = 500, target = 0.4, solve_for = "c1",
                        c2 = 10, consistency_method = g, draws = 10000, seed = 15)
    expect_true(inv$attainable)
    expect_true(is.finite(inv$value))
    G <- fs_oc_grid(family = fam, n = 500, c1 = inv$value, c2 = 10,
                    consistency_method = g, draws = 10000, block = Inf, seed = 15)
    expect_identical(G$table$det_rate, inv$achieved)          # same draws, same gate
    expect_gte(inv$achieved, 0.4)                              # the step the target falls on
    expect_lt(inv$achieved - 0.4, 0.01)                        # and not far above it
    expect_lt(inv$next_step_rate, 0.4)
    expect_equal(inv$achieved_se, sqrt(inv$achieved * (1 - inv$achieved) / 10000))
    expect_output(print(inv), "Declaration-rate inversion")
  }
  # solve for c2 at fixed c1
  inv2 <- fs_oc_invert(family = fam, n = 500, target = 0.5, solve_for = "c2",
                       c1 = 30, consistency_method = "split", draws = 10000, seed = 15)
  expect_true(inv2$attainable)
  G2 <- fs_oc_grid(family = fam, n = 500, c1 = 30, c2 = inv2$value,
                   consistency_method = "split", draws = 10000, block = Inf, seed = 15)
  expect_identical(G2$table$det_rate, inv2$achieved)
})


# ---- 6. ceiling ----------------------------------------------------------------

test_that("a target above the ceiling returns NA with the ceiling reported", {
  fam <- .grid_hand_family()
  inv <- fs_oc_invert(family = fam, n = 500, target = 0.999, solve_for = "c1",
                      c2 = 10, consistency_method = "resample", draws = 10000, seed = 2)
  expect_false(inv$attainable)
  expect_true(is.na(inv$value))
  expect_true(is.finite(inv$ceiling) && inv$ceiling < 0.999)
  expect_match(inv$binding, "c2")
  # the ceiling is the grid rate at c1 = -Inf
  G <- fs_oc_grid(family = fam, n = 500, c1 = -Inf, c2 = 10,
                  consistency_method = "resample", draws = 10000, block = Inf, seed = 2)
  expect_identical(G$table$det_rate, inv$ceiling)
  expect_output(print(inv), "infeasible")
})


# ---- 7. family reuse guard ---------------------------------------------------

test_that("a family built at one n is refused at another n", {
  fam <- .grid_hand_family(n = 500)
  expect_error(fs_oc_grid(family = fam, n = 700, c1 = 30, c2 = 10, draws = 100),
               "built at n = 500")
  expect_error(fs_oc_grid(family = fam, n = c(500, 700), c1 = 30, c2 = 10, draws = 100),
               "built at n = 500")
  expect_error(fs_oc_invert(family = fam, n = 700, target = 0.5, c2 = 10, draws = 100),
               "built at n = 500")
  fam_no_n <- fam; fam_no_n$n <- NULL
  expect_error(fs_oc_grid(family = fam_no_n, n = 500, c1 = 30, c2 = 10, draws = 100),
               "unrecorded")
  # enumerated families carry their n and are accepted at it
  dgm <- .grid_dgm()
  f5 <- fs_oc_family_enumerate(dgm, .grid_args, n = 500)
  G <- fs_oc_grid(family = f5, n = 500, c1 = 30, c2 = 10, draws = 500, seed = 1)
  expect_equal(unique(G$table$M), f5$M)
  expect_equal(nrow(G$table), 2L)          # both gates by default
  # and per-n enumeration inside the grid gives a different M at another n
  G2 <- fs_oc_grid(dgm, .grid_args, n = c(500, 700), c1 = 30, c2 = 10,
                   consistency_method = "resample", draws = 500, seed = 1)
  expect_equal(nrow(G2$table), 2L)
  expect_false(G2$families[["500"]]$M == G2$families[["700"]]$M)
  expect_equal(G2$families[["700"]]$floor, 60 / 700)
  expect_true(all(c("enumerate_secs", "draw_secs", "sweep_secs") %in% names(G2$timing)))
})
