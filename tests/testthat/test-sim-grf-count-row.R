# run_simulation_analysis(): the GRF row of a count simulation
# (TASK_cran_prep_finalize_2026-11-13, Step 1).  For outcome_type = "count",
# grf.subg.harm.glm() accepts only effect_measure "log_IRR" or "log_RR", and
# the forest fit does not use the value.  .run_grf_analysis_gen() passes
# "log_IRR" for the simulation's "IRR" and "IRD"; the GRF row's effects are
# computed by the simulation's extractor on the simulation's own scale.
# Before that change every count GRF row failed with "GRF failed: For
# outcome_type = 'count', effect_measure must be 'log_IRR' or 'log_RR'".

skip_on_cran()

test_that("count GRF rows complete for IRR and IRD with an exposure offset", {
  skip_if_not_installed("grf")
  skip_if_not_installed("policytree")
  old_cores <- Sys.getenv("_R_CHECK_LIMIT_CORES_", unset = NA)
  Sys.setenv("_R_CHECK_LIMIT_CORES_" = "TRUE")   # grf num.threads = 2
  on.exit(if (is.na(old_cores)) Sys.unsetenv("_R_CHECK_LIMIT_CORES_")
          else Sys.setenv("_R_CHECK_LIMIT_CORES_" = old_cores), add = TRUE)

  for (em in c("IRR", "IRD")) {
    dgm <- generate_glm_dgm(
      data = .make_count_data(N = 300L), factor_vars = c("biomarker_hi", "sex"),
      continuous_vars = "age", outcome_var = "y", treatment_var = "treat",
      outcome_type = "count", effect_measure = em, offset_var = "ftime",
      subgroup_vars = c("biomarker_hi", "sex"),
      subgroup_cuts = list(biomarker_hi = 1L, sex = 1L),
      model = "alt", k_inter = 2.0, n_super = 2000L, seed = 8316951L,
      verbose = FALSE)

    warns <- character(0)
    expect_no_error(
      r <- withCallingHandlers(
        suppressMessages(run_simulation_analysis(
          sim_id = 1L, dgm = dgm, n_sample = 600,
          confounders_base = c("biomarker_hi", "sex", "age"),
          run_fs = FALSE, run_fs_grf = FALSE, run_grf = TRUE,
          fs_params = list(outcome_type = "count", effect_measure = em,
                           offset.name = "ftime"),
          grf_params = list(dmin.grf = 0), verbose = FALSE)),
        warning = function(cnd) {
          warns <<- c(warns, conditionMessage(cnd))
          invokeRestart("muffleWarning")
        }))

    expect_false(any(grepl("effect_measure must be", warns, fixed = TRUE)),
                 info = em)
    expect_false(any(grepl("GRF failed", warns, fixed = TRUE)), info = em)
    expect_identical(r$analysis, "GRF", info = em)
  }
})
