# Smoke test for the MRCT structure-injection module (R/mrct_inject_dgm.R)
# on survival::cgd, time to first serious infection.
# TASK_fs-add-mrct-inject-module (2026-10-02), section 3.

test_that("inject_mrct_structure builds a usable aft_dgm_flex on cgd", {
  cgd <- survival::cgd
  cgd1 <- cgd[cgd$enum == 1, ]
  cgd1 <- data.frame(
    tte      = cgd1$tstop - cgd1$tstart,
    event    = cgd1$status,
    treat    = as.numeric(cgd1$treat == "rIFN-g"),
    age      = cgd1$age,
    height   = cgd1$height,
    weight   = cgd1$weight,
    female   = as.numeric(cgd1$sex == "female"),
    autosom  = as.numeric(cgd1$inherit == "autosomal"),
    steroids = cgd1$steroids,
    propylac = cgd1$propylac
  )

  dgm <- inject_mrct_structure(
    seed_data = cgd1, outcome_var = "tte", event_var = "event",
    treatment_var = "treat", continuous_vars = c("age", "height", "weight"),
    factor_vars = c("female", "autosom", "steroids", "propylac"), x_pred = "age",
    spline_spec = list(knot = 12, zeta = 25, log_hrs = log(c(0.35, 0.80, 1.30))),
    region = list(prevalence = 0.20, or_pred = 10, x1_vars = "female", or_x1 = 3,
                  loghr = 0),
    n_super = 5000, expand = "copula", seed = 1
  )

  expect_s3_class(dgm, "aft_dgm_flex")
  expect_equal(nrow(dgm$df_super), 5000)
  expect_true("z_region" %in% names(dgm$df_super))
  prev <- mean(dgm$df_super$z_region)
  expect_gte(prev, 0.17)
  expect_lte(prev, 0.23)

  truth <- mrct_truth_by_region(dgm)
  expect_equal(nrow(truth), 3)
  expect_true(all(is.finite(truth$AHR)))
  expect_gt(attr(truth, "smd_x_pred"), 0.3)

  sim <- simulate_from_dgm(dgm, n = 500, analysis_time = 500, max_entry = 120,
                           seed = 1)
  expect_equal(nrow(sim), 500)
  expect_true(all(sim$event_sim %in% c(0, 1)))

  met <- sim_region_metrics(sim, "z_region")
  expect_true(is.finite(met$overall["hr"]))

  covs <- analysis_covariates(dgm, "unobserved")
  expect_false("z_age" %in% covs)
  expect_false("z_region" %in% covs)
})

# TASK_fs-region-band-term (2026-10-03): band term in the region logit.
test_that("region band term concentrates the X3 band in the region on cgd", {
  cgd <- survival::cgd
  cgd1 <- cgd[cgd$enum == 1, ]
  cgd1 <- data.frame(
    tte      = cgd1$tstop - cgd1$tstart,
    event    = cgd1$status,
    treat    = as.numeric(cgd1$treat == "rIFN-g"),
    age      = cgd1$age,
    height   = cgd1$height,
    weight   = cgd1$weight,
    female   = as.numeric(cgd1$sex == "female"),
    autosom  = as.numeric(cgd1$inherit == "autosomal"),
    steroids = cgd1$steroids,
    propylac = cgd1$propylac
  )

  dgm <- inject_mrct_structure(
    seed_data = cgd1, outcome_var = "tte", event_var = "event",
    treatment_var = "treat", continuous_vars = c("age", "height", "weight"),
    factor_vars = c("female", "autosom", "steroids", "propylac"), x_pred = "age",
    spline_spec = list(knot = 12, zeta = 25, log_hrs = rep(log(0.70), 3)),
    region = list(prevalence = 0.20, or_pred = 1,
                  band = list(cut = 10, or = 20)),
    x3 = list(vars = "age", cuts = list(age = 10), loghr = log(5)),
    n_super = 5000, expand = "copula", seed = 1
  )

  ds <- dgm$df_super
  prev <- mean(ds$z_region)
  expect_gte(prev, 0.17)
  expect_lte(prev, 0.23)

  expect_true(all(ds$flag_harm == (ds$z_age <= 10)))

  in_band <- ds$z_age <= 10
  expect_gt(mean(in_band[ds$z_region == 1]), mean(in_band[ds$z_region == 0]))

  expect_equal(dgm$mrct$region$band, list(var = "age", cut = 10, or = 20))
  expect_equal(dgm$mrct$region_model$band, dgm$mrct$region$band)
})
