# inject_mrct_structure() observed-region mode on ACTG 175.
# TASK_fs-inject-observed-region_2026-11-02, section 3.

obs_seed_data <- function() {
  d <- mrct_seed_data()                     # helper-mrct-analytic.R
  set.seed(20261102L)
  d$reg_obs <- as.integer(d$symptom == 1 | stats::rbinom(nrow(d), 1L, 0.15) == 1L)
  d
}

obs_build <- function(region = list(observed = "reg_obs"), expand = "bootstrap",
                      seed_data = obs_seed_data(), n_super = 3000L, seed = 11L, ...) {
  inject_mrct_structure(
    seed_data,
    outcome_var = "time_months", event_var = "cens", treatment_var = "treat_mrct",
    continuous_vars = c("age", "wtkg", "karnof", "cd40"),
    factor_vars = c("hemo", "homo", "drugs", "symptom"),
    x_pred = "cd40",
    spline_spec = list(knot = 300, zeta = 800, log_hrs = rep(log(0.55), 3)),
    region = region,
    x3 = list(vars = "cd40", cuts = list(cd40 = 230), loghr = log(2 / 0.55)),
    n_super = n_super, expand = expand, seed = seed,
    continuous_vars_cens = c("age", "cd40"), factor_vars_cens = "symptom", ...)
}

obs_memo <- new.env()
obs_dgm <- function() {
  if (is.null(obs_memo$boot)) obs_memo$boot <- obs_build()
  obs_memo$boot
}

test_that("observed-region build succeeds and records the mode", {
  skip_if_not_installed("speff2trial")
  dgm <- obs_dgm()
  expect_s3_class(dgm, "aft_dgm_flex")
  r <- dgm$mrct$region
  expect_true(isTRUE(r$observed))
  expect_identical(r$name, "reg_obs")
  expect_identical(r$z_name, "z_reg_obs")
  expect_true(r$z_name %in% names(dgm$df_super))
  zr <- dgm$df_super[[r$z_name]]
  expect_true(all(zr %in% c(0, 1)))
  expect_setequal(unique(zr), c(0, 1))
  expect_equal(r$prevalence_pop, mean(zr))
  for (f in c("prevalence_target", "alpha0", "coefs", "scale", "band"))
    expect_null(r[[f]], label = f)
  expect_true("region_model" %in% names(dgm$mrct))
  expect_null(dgm$mrct$region_model)
  expect_identical(dgm$mrct$classes$X1, character(0))
})

test_that("with expand = 'bootstrap' each df_super row keeps its source row's region", {
  skip_if_not_installed("speff2trial")
  dgm <- obs_dgm()
  d <- obs_seed_data()
  use <- c("time_months", "cens", "treat_mrct", "age", "wtkg", "karnof", "cd40",
           "hemo", "homo", "drugs", "symptom", "reg_obs")
  d <- d[stats::complete.cases(d[, use]), ]
  set.seed(11L)                               # expand_covariates() seed = seed
  idx <- sample.int(nrow(d), 3000L, replace = TRUE)
  ## alignment: the covariates are the source rows'
  expect_equal(dgm$df_super$z_cd40, d$cd40[idx])
  expect_equal(dgm$df_super$z_age, d$age[idx])
  ## the region is the source row's
  expect_equal(dgm$df_super$z_reg_obs, as.numeric(d$reg_obs[idx]))
})

test_that("loghr = NULL estimates the region log-HR; numeric loghr fixes it", {
  skip_if_not_installed("speff2trial")
  dgm <- obs_dgm()
  pl <- dgm$mrct$region$prognostic_loghr
  expect_true(is.numeric(pl) && length(pl) == 1L && is.finite(pl))
  expect_equal(pl, unname(dgm$model_params$b0["z_reg_obs"]))

  fx <- obs_build(region = list(observed = "reg_obs", loghr = -log(3)))
  expect_equal(fx$mrct$region$prognostic_loghr, -log(3))
  expect_equal(unname(fx$model_params$b0["z_reg_obs"]), -log(3))
})

test_that("planted-model elements and a one-level region column are errors", {
  skip_if_not_installed("speff2trial")
  bad <- list(prevalence = 0.2, or_pred = 2, x1_vars = "symptom", or_x1 = 3,
              scale = "rank", band = list(cut = 230, or = 2))
  for (nm in names(bad)) {
    reg <- c(list(observed = "reg_obs"), bad[nm])
    expect_error(obs_build(region = reg, n_super = 200L), nm, fixed = TRUE,
                 label = nm)
  }
  d1 <- obs_seed_data(); d1$reg_obs <- 1L
  expect_error(obs_build(seed_data = d1, n_super = 200L), "one level")
})

test_that("ahr_population() accepts the DGM and ahr_limit() on the region is finite", {
  skip_if_not_installed("speff2trial")
  dgm <- obs_dgm()
  pop <- ahr_population(dgm, analysis_time = 40, max_entry = 11)
  expect_s3_class(pop, "ahr_population")
  expect_equal(mean(pop$region), mean(dgm$df_super$z_reg_obs))
  e <- ahr_limit(pop, "region")
  for (f in c("beta", "hr", "ahr", "cde", "var_robust", "d"))
    expect_true(all(is.finite(e[[f]])), label = f)
})

test_that("observed region with expand = 'copula' enters the latent correlation", {
  skip_if_not_installed("speff2trial")
  skip_if_not_installed("latentcor")
  dgm <- obs_build(expand = "copula")
  expect_true(isTRUE(dgm$mrct$region$observed))
  R <- dgm$mrct$latent_R
  expect_true("reg_obs" %in% rownames(R))
  expect_true(all(dgm$df_super$z_reg_obs %in% c(0, 1)))
  ## a latent_R without the region column fails the existing check
  keep <- setdiff(rownames(R), "reg_obs")
  expect_error(obs_build(expand = "copula", n_super = 200L,
                         latent_R = R[keep, keep]),
               "latent_R names must match")
})
