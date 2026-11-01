# The analytic layer (R/ahr_analytic.R) against the committed mrct-consistency
# records (fixtures/mrct_analytic/, commits in fixtures/README.md), and the
# S2 build of the records (TASK_fs-mrct-analytic-layer_2026-11-01, section 5,
# tests 1-3). The builders are in helper-mrct-analytic.R.

tol_ref <- 1e-9   # probabilities absolute; other quantities relative to max(1, |ref|)
near <- function(ours, ref, tol = tol_ref) abs(ours - ref) <= tol * pmax(1, abs(ref))

# ---- 1. build identity ----------------------------------------------------

test_that("the S2 build matches the committed build.csv row for S2_rmlog4", {
  mrct_skip_if_no_build()
  dgm <- mrct_build_S2()
  bc <- read.csv(mrct_fixture("build_S2_rmlog4.csv"))
  expect_equal(nrow(bc), 1L)
  pop <- mrct_pop_S2()
  e_reg <- ahr_limit(pop, "region", b = 0.55, h = 2, variance = FALSE)
  e_non <- ahr_limit(pop, "nonregion", b = 0.55, h = 2, variance = FALSE)
  ds <- dgm$df_super; reg <- ds$z_region == 1; H <- ds$flag_harm == 1; rg <- dgm$mrct$region
  ours <- c(
    h = 2, r = -log(4), mu = unname(dgm$model_params$mu), mu_fitted = unname(dgm$mrct$mu_fitted),
    tau = dgm$model_params$tau, m = dgm$mrct$hazard_multiplier,
    b0_z_region = unname(dgm$model_params$b0["z_region"]), b0_treat = unname(dgm$model_params$b0["treat"]),
    b0_treat_harm = unname(dgm$model_params$b0["treat_harm"]),
    prognostic_loghr_recorded = rg$prognostic_loghr, region_alpha0 = rg$alpha0,
    region_coef_cd40 = unname(rg$coefs["cd40"]), region_coef_symptom = unname(rg$coefs["symptom"]),
    region_coef_cd4_le_230 = log(rg$band$or), prevalence = mean(reg),
    Hstar_share_region = mean(H[reg]), Hstar_share_nonregion = mean(H[!reg]),
    event_prob_40_region = e_reg$d, event_prob_40_nonregion = e_non$d,
    events_n1000_region = 1000 * mean(reg) * e_reg$d, events_n1000_nonregion = 1000 * mean(!reg) * e_non$d,
    mean_theta0_region = mean(ds$theta_0[reg]), mean_theta0_nonregion = mean(ds$theta_0[!reg]),
    n_build_warnings = length(attr(dgm, "build_warnings")))
  # every numeric column of the record that describes the build; the campaign's
  # bookkeeping (seed, truth_seed, replicates) and the simulated trials'
  # prevalence_trials are not properties of the build
  num <- names(bc)[vapply(bc, is.numeric, NA)]
  expect_setequal(setdiff(num, names(ours)), c("prevalence_trials", "seed", "truth_seed", "replicates"))
  for (k in names(ours)) expect_equal(ours[[k]], bc[[k]], tolerance = 1e-8, label = k)
  expect_equal(attr(dgm, "m"), 3)
})

# ---- 2. the layer at no delay ---------------------------------------------

test_that("ahr_design() and ahr_limit() reproduce the S2 predictions at r = -log 4", {
  mrct_skip_if_no_build()
  pop <- mrct_pop_S2()
  pr <- read.csv(mrct_fixture("predictions_S2_rmlog4.csv"), stringsAsFactors = FALSE)
  expect_true(all(pr$b == 0.55 & pr$h == 2 & pr$n == 1000))
  dz <- ahr_design(pop, 0.55, 2, 1000, sets = c("all", "nonregion", "region", "region_Gstar"))
  pk <- c(p_A = "p_itt_ub_lt1", p_B = "p_region_itt_gt1", p_A_B = "p_trigger",
          p_Cstar = "p_region_Gstar_lt1", p_A_B_Cstar = "p_trigger_and_Gstar_lt1")
  for (k in names(pk)) {
    ref <- pr$value[pr$quantity == k]
    expect_length(ref, 1L)
    expect_lte(abs(dz$prob[[pk[[k]]]] - ref), tol_ref)
  }
  st <- dz$sets
  for (s in st$set) {
    row <- st[st$set == s, ]
    for (m in list(c("se_model_loghr", "se_model"), c("se_robust_loghr", "se_robust"), c("events", "events"))) {
      ref <- pr$value[pr$quantity == m[1] & pr$set == s]
      expect_length(ref, 1L)
      expect_true(near(row[[m[2]]], ref), label = paste(m[1], s))
    }
  }
  # ahr_limit() on the six sets of the record
  setmap <- c(all = "all", nonregion = "nonregion", region = "region", region_Gstar = "region_Gstar",
              band_nonregion = "Hstar_nonregion", band_region = "Hstar_region")
  for (s in names(setmap)) {
    e <- ahr_limit(pop, setmap[[s]], b = 0.55, h = 2, variance = FALSE)
    got <- c(marghr_40 = e$hr, marghr_inf = e$marghr_inf, ahr = e$ahr, cde = e$cde,
             share_of_population = e$share, Hstar_share = e$band_share, event_prob_40 = e$d)
    for (k in names(got)) {
      ref <- pr$value[pr$quantity == k & pr$set == s]
      expect_length(ref, 1L)
      expect_true(near(got[[k]], ref), label = paste(k, s))
    }
  }
})

test_that("ahr_design() reproduces the grid.csv row (b 0.55, h 2, n 1,000, r = -log 2)", {
  mrct_skip_if_no_build()
  pop <- mrct_pop_S2(-log(2))
  gt <- read.csv(mrct_fixture("grid_b055_h2_n1000_rmlog2.csv"), stringsAsFactors = FALSE)
  expect_equal(nrow(gt), 1L)
  expect_true(abs(gt$r + log(2)) < 1e-12 && gt$b == 0.55 && gt$h == 2 && gt$n == 1000)
  dz <- ahr_design(pop, 0.55, 2, 1000, sets = c("all", "nonregion", "region", "region_Gstar"))
  s <- function(set, col) dz$sets[[col]][dz$sets$set == set]
  probs <- c(p_itt_ub_lt1 = "p_itt_ub_lt1", p_region_itt_gt1 = "p_region_itt_gt1", p_trigger = "p_trigger",
             p_region_Gstar_lt1 = "p_region_Gstar_lt1", p_oracle_joint = "p_trigger_and_Gstar_lt1")
  for (k in names(probs)) expect_lte(abs(dz$prob[[probs[[k]]]] - gt[[k]]), tol_ref, label = k)
  hs <- ahr_limit(pop, "Hstar_region", 0.55, 2, variance = FALSE)
  ea <- ahr_limit(pop, "all", 0.55, 2, variance = FALSE)
  er <- ahr_limit(pop, "region", 0.55, 2, variance = FALSE)
  got <- c(marghr40_all = s("all", "marghr"), marghr40_region = s("region", "marghr"),
           marghr40_nonregion = s("nonregion", "marghr"), marghr40_Gstar = s("region_Gstar", "marghr"),
           marghr40_Hstar_region = hs$hr, marghrinf_all = ea$marghr_inf, marghrinf_region = er$marghr_inf,
           events_all = s("all", "events"), events_region = s("region", "events"),
           events_Gstar = s("region_Gstar", "events"),
           se_model_all = s("all", "se_model"), se_model_region = s("region", "se_model"),
           se_robust_all = s("all", "se_robust"), se_robust_region = s("region", "se_robust"))
  for (k in names(got)) expect_true(near(got[[k]], gt[[k]]), label = k)
})

# ---- 3. the layer with a delay --------------------------------------------

test_that("ahr_design() at the two committed roots reproduces roots.csv", {
  mrct_skip_if_no_build()
  rt <- read.csv(mrct_fixture("roots.csv"), stringsAsFactors = FALSE)
  roots <- rt[startsWith(rt$row, "D1") | startsWith(rt$row, "D2"), ]
  expect_equal(nrow(roots), 2L)
  for (i in seq_len(nrow(roots))) {
    x <- roots[i, ]
    pop <- mrct_pop_S2(x$r, entry_delay = list(var = "z_region", delay = x$a))
    dz <- ahr_design(pop, 0.55, 2, 1000)
    s <- function(set, col) dz$sets[[col]][dz$sets$set == set]
    for (k in c("p_itt_ub_lt1", "p_region_itt_gt1", "p_trigger"))
      expect_lte(abs(dz$prob[[k]] - x[[k]]), tol_ref, label = paste(x$row, k))
    expect_lte(abs(dz$prob[["p_trigger_and_Gstar_lt1"]] - x$p_trigger_Gstar_lt1), tol_ref)
    got <- c(marghr40_region = s("region", "marghr"), marghr40_all = s("all", "marghr"),
             events_region = s("region", "events"), events_all = s("all", "events"))
    for (k in names(got)) expect_true(near(got[[k]], x[[k]]), label = paste(x$row, k))
  }
})

test_that("ahr_limit(region) with a region delay a equals no delay at horizon 40 - a", {
  mrct_skip_if_no_build()
  a <- 6
  e6 <- ahr_limit(mrct_pop_S2(-log(2), entry_delay = list(var = "z_region", delay = a)),
                  "region", b = 0.55, h = 2)
  e34 <- ahr_limit(mrct_pop_S2(-log(2)), "region", b = 0.55, h = 2, horizon = 40 - a)
  expect_lte(abs(e6$beta - e34$beta), 1e-6)
  expect_lte(abs(e6$i1 / e34$i1 - 1), 1e-6)
  expect_lte(abs(e6$d / e34$d - 1), 1e-6)
  expect_lte(abs(e6$var_robust / e34$var_robust - 1), 1e-6)
  # the non-region is untouched by the region's delay
  n6 <- ahr_limit(mrct_pop_S2(-log(2), entry_delay = list(var = "z_region", delay = a)),
                  "nonregion", b = 0.55, h = 2, variance = FALSE)
  n0 <- ahr_limit(mrct_pop_S2(-log(2)), "nonregion", b = 0.55, h = 2, variance = FALSE)
  expect_lte(abs(n6$beta - n0$beta), 1e-6)
})

test_that("ahr_solve() finds the committed D1 root on its bracket", {
  mrct_skip_if_no_build()
  rt <- read.csv(mrct_fixture("roots.csv"), stringsAsFactors = FALSE)
  d1 <- rt[startsWith(rt$row, "D1"), ]
  dgm <- mrct_build_S2(-log(2))
  f <- function(a) ahr_design(ahr_population(dgm, 40, 11, entry_delay = list(var = "z_region", delay = a)),
                              0.55, 2, 1000)$prob[["p_trigger"]]
  sol <- ahr_solve(f, 0.25, c(11, 12), scan = c(11, 12), tol = 1e-8)
  expect_equal(sol$n_roots, 1L)
  expect_lte(abs(sol$roots - d1$a), 1e-6)
  expect_lte(abs(sol$values - 0.25), 1e-6)
})

# ---- interfaces -----------------------------------------------------------

test_that("ahr_solve() returns every bracketed root and exact scan-point roots", {
  f <- function(x, k) (x - 1) * (x - 3) * (x - k)
  s <- ahr_solve(f, 0, c(0, 6), k = 5, n_scan = 13, tol = 1e-12)
  expect_equal(s$n_roots, 3L)
  expect_equal(s$roots, c(1, 3, 5), tolerance = 1e-8)
  s2 <- ahr_solve(function(x) x - 2, 0, c(0, 4), scan = 0:4)
  expect_identical(s2$roots, 2)
  s3 <- ahr_solve(function(x) x^2 + 1, 0, c(-1, 1))
  expect_equal(s3$n_roots, 0L)
  expect_error(ahr_solve(f, 0, c(2, 1), k = 5), "lower < upper")
})

test_that("ahr_population() and ahr_limit() validate their inputs", {
  mrct_skip_if_no_build()
  dgm <- mrct_build_S2()
  pop <- mrct_pop_S2()
  expect_s3_class(pop, "ahr_population")
  expect_identical(pop$N, nrow(dgm$df_super))
  expect_output(print(pop), "ahr_population")
  expect_error(ahr_population(dgm, 40, 11, entry_delay = list(var = "nope", delay = 1)), "entry_delay")
  expect_error(ahr_population(dgm, 40, 11, entry_delay = list(var = "z_region", delay = 30)), "positive")
  expect_error(ahr_population(dgm, 40, 11, entry_delay = 1:3), "length")
  bad <- dgm; bad$df_super$theta_0 <- bad$df_super$theta_0 + 1
  expect_error(ahr_population(bad, 40, 11), "not consistent")
  expect_error(ahr_limit(pop, "nope", 0.55, 2), "set must be one of")
  expect_error(ahr_limit(pop, "all", b = 0.55), "both b and h")
  # a logical set equal to a standard set agrees with it (own grid vs cached grid)
  e1 <- ahr_limit(pop, pop$region, 0.55, 2, variance = FALSE)
  e2 <- ahr_limit(pop, "region", 0.55, 2, variance = FALSE)
  expect_lte(abs(e1$beta - e2$beta), 1e-6)
  # the DGM's own effects (psi = loghr_po) equal (b, h) = (0.55, 2) on this build
  e3 <- ahr_limit(pop, "region", variance = FALSE)
  expect_lte(abs(e3$beta - e2$beta), 1e-6)
})
