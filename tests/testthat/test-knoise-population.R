# Tests for the population-noise scheme in the DGM builders:
# setup_gbsg_dgm() / generate_glm_dgm() gain k_random_noise / noise_seed and
# draw noise ONCE onto the super-population at construction (scheme of commit
# 2a4787bc, previously inline in the gbsg_redux sweep drivers).

# Survival builder shared by several tests (small n_super for speed).
.make_gbsg_noise_dgm <- function(...) {
  setup_gbsg_dgm(model = "alt", n_super = 1000L, verbose = FALSE, ...)
}

# Minimal GLM fixture: no test previously exercised generate_glm_dgm, so this
# builds on .make_binary_data() from helper-synthetic-dgm.R.
.make_glm_noise_dgm <- function(...) {
  generate_glm_dgm(
    data          = .make_binary_data(N = 200L),
    factor_vars   = c("biomarker_hi", "sex"),
    outcome_var   = "y",
    treatment_var = "treat",
    outcome_type  = "binary",
    subgroup_vars = c("biomarker_hi", "sex"),
    subgroup_cuts = list(biomarker_hi = 1L, sex = 1L),
    model         = "alt",
    k_inter       = 0.5,
    n_super       = 1000L,
    verbose       = FALSE,
    ...
  )
}

# -- 1. Survival: default k_random_noise = 0 is noise-free -------------------

test_that("setup_gbsg_dgm default has no noise columns and 'none' scheme", {
  dgm <- .make_gbsg_noise_dgm()
  expect_false(any(grepl("^noise[0-9]+$", names(dgm$df_super))))
  expect_false(any(grepl("^noise[0-9]+$", names(dgm$df_super_rand))))
  expect_identical(dgm$noise_names, character(0))
  expect_identical(dgm$noise_scheme, "none")
  expect_identical(dgm$noise_seed, NA_integer_)
})

# -- 2. Survival: k = 3 lands on BOTH frames, identical values ---------------

test_that("setup_gbsg_dgm k_random_noise = 3 puts identical noise on both frames", {
  dgm <- .make_gbsg_noise_dgm(k_random_noise = 3L)
  expect_identical(dgm$noise_names, c("noise1", "noise2", "noise3"))
  expect_true(all(dgm$noise_names %in% names(dgm$df_super)))
  expect_true(all(dgm$noise_names %in% names(dgm$df_super_rand)))
  for (nm in dgm$noise_names) {
    expect_identical(dgm$df_super[[nm]], dgm$df_super_rand[[nm]])
    expect_length(dgm$df_super[[nm]], 1000L)
  }
  expect_identical(dgm$noise_scheme, "population")
  expect_identical(dgm$noise_seed, 20260807L)
})

# -- 3. Survival: reproducible by noise_seed, sensitive to it ----------------

test_that("setup_gbsg_dgm noise is reproducible and seed-sensitive", {
  d_a <- .make_gbsg_noise_dgm(k_random_noise = 3L)
  d_b <- .make_gbsg_noise_dgm(k_random_noise = 3L)
  d_c <- .make_gbsg_noise_dgm(k_random_noise = 3L, noise_seed = 1L)
  for (nm in d_a$noise_names) {
    expect_identical(d_a$df_super[[nm]], d_b$df_super[[nm]])
    expect_false(identical(d_a$df_super[[nm]], d_c$df_super[[nm]]))
  }
  expect_identical(d_c$noise_seed, 1L)
})

# -- 4. Survival: non-noise columns golden vs the k = 0 build ----------------

test_that("setup_gbsg_dgm k_random_noise perturbs nothing but the noise columns", {
  d0 <- .make_gbsg_noise_dgm()
  d3 <- .make_gbsg_noise_dgm(k_random_noise = 3L)
  expect_identical(d3$df_super[names(d0$df_super)], d0$df_super)
  expect_identical(d3$df_super_rand[names(d0$df_super_rand)], d0$df_super_rand)
  expect_setequal(names(d3$df_super), c(names(d0$df_super), d3$noise_names))
})

# -- 5. Survival: inheritance into trials and the eval frame -----------------

test_that("noise columns inherit into simulate_from_dgm and fs_build_eval_frame", {
  dgm <- .make_gbsg_noise_dgm(k_random_noise = 3L)

  sim <- simulate_from_dgm(dgm, n = 200, seed = 1)
  expect_true(all(dgm$noise_names %in% names(sim)))
  for (nm in dgm$noise_names) {
    expect_true(all(sim[[nm]] %in% dgm$df_super[[nm]]))
  }

  dgm <- compute_dgm_cde(dgm)
  ev <- fs_build_eval_frame(dgm)
  expect_true(all(dgm$noise_names %in% names(ev)))
  for (nm in dgm$noise_names) {
    expect_true(all(ev[[nm]] %in% dgm$df_super[[nm]]))
  }
})

# -- 6. GLM mirror of 1-4, inheritance via simulate_from_glm_dgm -------------

test_that("generate_glm_dgm default has no noise columns and 'none' scheme", {
  g <- .make_glm_noise_dgm()
  expect_false(any(grepl("^noise[0-9]+$", names(g$df_super))))
  expect_identical(g$noise_names, character(0))
  expect_identical(g$noise_scheme, "none")
  expect_identical(g$noise_seed, NA_integer_)
})

test_that("generate_glm_dgm k_random_noise = 3 draws population noise", {
  g <- .make_glm_noise_dgm(k_random_noise = 3L)
  expect_identical(g$noise_names, c("noise1", "noise2", "noise3"))
  expect_true(all(g$noise_names %in% names(g$df_super)))
  for (nm in g$noise_names) expect_length(g$df_super[[nm]], 1000L)
  expect_identical(g$noise_scheme, "population")
  expect_identical(g$noise_seed, 20260807L)
})

test_that("generate_glm_dgm noise is reproducible and seed-sensitive", {
  g_a <- .make_glm_noise_dgm(k_random_noise = 3L)
  g_b <- .make_glm_noise_dgm(k_random_noise = 3L)
  g_c <- .make_glm_noise_dgm(k_random_noise = 3L, noise_seed = 1L)
  for (nm in g_a$noise_names) {
    expect_identical(g_a$df_super[[nm]], g_b$df_super[[nm]])
    expect_false(identical(g_a$df_super[[nm]], g_c$df_super[[nm]]))
  }
})

test_that("generate_glm_dgm k_random_noise perturbs nothing but the noise columns", {
  g0 <- .make_glm_noise_dgm()
  g3 <- .make_glm_noise_dgm(k_random_noise = 3L)
  expect_identical(g3$df_super[names(g0$df_super)], g0$df_super)
  expect_setequal(names(g3$df_super), c(names(g0$df_super), g3$noise_names))
})

test_that("noise columns inherit into simulate_from_glm_dgm", {
  g <- .make_glm_noise_dgm(k_random_noise = 3L)
  sim <- simulate_from_glm_dgm(g, n = 200, seed = 1)
  expect_true(all(g$noise_names %in% names(sim)))
  for (nm in g$noise_names) {
    expect_true(all(sim[[nm]] %in% g$df_super[[nm]]))
  }
})
