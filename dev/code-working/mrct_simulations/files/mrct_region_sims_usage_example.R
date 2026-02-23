# =============================================================================
# Example: Replicating a specific forestsearch() configuration inside
#          mrct_region_sims() using the new fs_args passthrough
# =============================================================================
#
# The original manual call was:
#
#   fs <- forestsearch(dfa,
#     confounders.name = confounders.name,
#     outcome.name = "y_sim", treat.name = "treat_sim",
#     event.name = "event_sim", id.name = "id",
#     potentialOutcome.name = "loghr_po",
#     df.test = as.data.frame(df_AP),
#     flag_harm.name = NULL,
#     hr.threshold = 0.9, hr.consistency = 0.80,
#     pconsistency.threshold = 0.80,
#     stop_threshold = NULL,
#     sg_focus = "minSG", max_subgroups_search = 30,
#     showten_subgroups = TRUE, details = TRUE,
#     conf_force = c("z_age <= 65", "z_bm <= 0", "z_bm <= 1",
#                    "z_bm <= 2", "z_bm <= 5"),
#     cut_type = "default", use_grf = TRUE, plot.grf = TRUE,
#     use_lasso = TRUE,
#     maxk = 1, n.min = 60, d0.min = 12, d1.min = 12,
#     plot.sg = TRUE, by.risk = 6,
#     use_twostage = FALSE,
#     twostage_args = list(n.splits.screen = 50, batch.size = 25,
#                          conf.level = 0.99),
#     parallel_args = list(plan = "callr", workers = 100,
#                          show_message = TRUE)
#   )
#
# Below maps each parameter into the mrct_region_sims() call.
# =============================================================================

library(forestsearch)
library(doFuture)
library(doRNG)

# ---- Step 1: Create the DGM ----

dgm_spline <- generate_aft_dgm_flex(
  data = df.case,
  continuous_vars = c("age", "bm"),
  factor_vars = c("male", "histology", "prior_treat", "AP"),
  set_beta_spec = list(
    set_var = c("z_AP"),
    beta_var = -log(5)
  ),
  continuous_vars_cens = c("age"),
  factor_vars_cens = c("prior_treat"),
  cens_type = "weibull",
  outcome_var = "tte",
  event_var = "event",
  treatment_var = "treat",
  subgroup_vars = NULL,
  subgroup_cuts = NULL,
  k_inter = 0.0,
  model = "alt",
  spline_spec = list(
    var = "z_bm",
    knot = 5,
    zeta = 10,
    log_hrs = log(c(2, 1.25, 0.5))
  ),
  verbose = TRUE,
  standardize = FALSE
)


# ---- Step 2: Run simulations ----
#
# Parameter mapping:
#
#   EXPOSED directly by mrct_region_sims():
#     hr.threshold, hr.consistency, pconsistency.threshold,
#     sg_focus, maxk, confounders.name, conf_force,
#     analysis_time, cens_adjust, region_var
#
#   PASSED VIA fs_args (the new passthrough):
#     use_grf, use_lasso, cut_type, n.min, d0.min, d1.min,
#     max_subgroups_search, use_twostage, twostage_args
#
#   PASSED VIA sim_args:
#     rand_ratio, draw_treatment
#
#   FORCED INTERNALLY (appropriate for batch sims):
#     showten_subgroups = FALSE, details = FALSE, plot.sg = FALSE,
#     plot.grf = FALSE, parallel_args = list(plan = "sequential")
#     (inner forestsearch runs sequentially; outer loop is parallelized)

results <- mrct_region_sims(
  dgm       = dgm_spline,
  n_sims    = 1000,
  n_sample  = 500,
  region_var = "z_AP",

  # --- Parameters exposed directly ---
  sg_focus               = "minSG",
  maxk                   = 1,
  hr.threshold           = 0.90,
  hr.consistency         = 0.80,
  pconsistency.threshold = 0.80,   # NOTE: default is 0.90; your call uses 0.80

  confounders.name = c("z_age", "z_bm", "z_male", "ecog",
                        "z_histology", "z_prior_treat", "strat"),
  conf_force       = c("z_age <= 65", "z_bm <= 0", "z_bm <= 1",
                        "z_bm <= 2", "z_bm <= 5"),

  # --- ForestSearch passthrough (fs_args) ---
  # These flow directly to forestsearch() inside each replicate.
  # Explicit mrct_region_sims params (above) take precedence.
  fs_args = list(
    use_grf              = TRUE,
    use_lasso            = TRUE,
    cut_type             = "default",
    n.min                = 60,
    d0.min               = 12,
    d1.min               = 12,
    max_subgroups_search = 30,
    use_twostage         = FALSE,
    twostage_args        = list(
      n.splits.screen = 50,
      batch.size      = 25,
      conf.level      = 0.99
    )
  ),

  # --- simulate_from_dgm passthrough (sim_args) ---
  sim_args = list(
    rand_ratio     = 1,
    draw_treatment = FALSE   # Retain original treatment assignment from df.case
  ),

  # --- Simulation infrastructure ---
  analysis_time = 60,
  cens_adjust   = 0,
  parallel_args = list(plan = "callr", workers = 100, show_message = TRUE),
  details       = TRUE,
  seed          = 20240101
)


# ---- Step 3: Summarize ----
cat("Subgroup identification rate:",
    round(100 * mean(results$any_found), 2), "%\n")

temp <- summaryout_mrct(
  pop_summary = NULL,
  mrct_sims   = results,
  showtable   = TRUE
)

plot_out <- SGplot_estimates(
  df             = results,
  label_training = "Non-AP, ITT",
  label_itt      = "Overall, ITT",
  label_testing  = "AP, ITT",
  label_sg       = "AP, identified subgroup"
)
print(plot_out$plot_estimates)
