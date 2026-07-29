# =============================================================================
# Phase 0 -- cell definitions
# =============================================================================
# The six `consistency` (fs_*) cells in quarto/simulations/gbsg_redux that have
# a pooled 500-replicate bundle.  Every knob below was read out of the cell's
# own batch .qmd; the six .qmd files are IDENTICAL apart from the three fields
# that vary here (target_hr_harm, n_sample, k_random_noise) -- verified by
# diffing the setup chunks.

PHASE0_CELLS <- list(
  list(cell = "h10_knoise0_n500",  target_hr_harm = 1.0, n_sample = 500L,  k_random_noise = 0L),
  list(cell = "h10_knoise0_n1000", target_hr_harm = 1.0, n_sample = 1000L, k_random_noise = 0L),
  list(cell = "h15_knoise0_n500",  target_hr_harm = 1.5, n_sample = 500L,  k_random_noise = 0L),
  list(cell = "h15_knoise0_n1000", target_hr_harm = 1.5, n_sample = 1000L, k_random_noise = 0L),
  list(cell = "h15_knoise3_n1500", target_hr_harm = 1.5, n_sample = 1500L, k_random_noise = 3L),
  list(cell = "h20_knoise3_n1500", target_hr_harm = 2.0, n_sample = 1500L, k_random_noise = 3L)
)
names(PHASE0_CELLS) <- vapply(PHASE0_CELLS, `[[`, character(1), "cell")

# Knobs shared by every cell (verbatim from the .qmd setup chunks).
PHASE0_COMMON <- list(
  seed_base          = 8316951L,
  dgm_model          = "alt",
  analysis_time      = 84,
  cens_adjust        = log(1.5),
  n_super            = 100000L,
  subgroup_method    = "consistency",
  consistency_method = "resample",
  sg_focus           = "eff",
  selection_rule     = "neighborhood",
  hr_threshold       = 0.90,
  hr_consistency     = 0.80,
  pconsistency       = 0.90,
  fs_splits          = 400L,
  maxk               = 2L,
  n_min              = NULL,
  d0_min             = 10L,
  d1_min             = 10L,
  use_lasso          = FALSE,
  use_dina           = FALSE,
  use_grf            = FALSE,
  use_twostage       = TRUE,
  fs_conf_force      = c("meno == 0", "er <= 0", "pgr <= 0"),
  fs_conf.cont_jcuts = list(er = 10),
  confounders_base   = c("er", "age", "meno", "pgr", "nodes", "size", "grade"),
  outcome_name       = "y_sim",
  event_name         = "event_sim",
  treat_name         = "treat_sim",
  id_name            = "id",
  harm_col           = "flag_harm"
)

SIM_DIR <- "quarto/simulations/gbsg_redux"

phase0_bundle_path <- function(cell) {
  file.path(SIM_DIR, sprintf("fs_t1_t2_m1_%s_combined_1_500.rds",
                             sub("^h", "h", cell)))
}
