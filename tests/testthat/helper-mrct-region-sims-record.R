# The small mrct_region_sims() call of test-mrct-region-sims-record.R, on the
# S2 build of helper-mrct-analytic.R: the campaign's specification
# (mrct-consistency sims/121_ident_rlog4_full_call.R) at maxk = 1 and
# fs.splits = 200, 4 replicates, 2 workers, fixed seed. Also sourced by the
# script that wrote fixtures/mrct_region_sims_record/point_6b5c85c0.rds, so
# the fixture and the test make the same call.

mrct_rs_n_sims <- 4L
mrct_rs_seed <- 20261101L

mrct_rs_call <- function(dgm, ..., n_sims = mrct_rs_n_sims, seed = mrct_rs_seed,
                         sim_args = list(rand_ratio = 1, draw_treatment = TRUE,
                                         max_entry = 11, cens_adjust = 50)) {
  args <- list(
    dgm = dgm, n_sims = n_sims, n_sample = 1000, region_var = "z_region",
    sg_focus = "hrMaxSG", maxk = 1,
    hr.threshold = 0.90, hr.consistency = 0.80, pconsistency.threshold = 0.90,
    confounders.name = analysis_covariates(dgm, "observed"),
    conf_force = character(0),
    fs_args = list(conf.cont_jcuts = list(z_cd40 = 10, z_wtkg = 10, z_preanti = 10),
                   effect_neighborhood = 0.20, use_grf = FALSE, use_lasso = FALSE,
                   fs.splits = 200),
    sim_args = sim_args, analysis_time = 40, cens_adjust = 50,
    parallel_args = list(plan = "multisession", workers = 2, show_message = FALSE),
    details = FALSE, seed = seed)
  do.call(mrct_region_sims, c(args, list(...)))
}
