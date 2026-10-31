# ---- transplant (verbatim chunk bodies; fences removed) ----
# source: forestsearch/quarto/applications/actg175/analysis_actg175_binary_multimethod_fixed_family.qmd
# forestsearch commit 43b051b6 (2026-08-16 12:07; the render that wrote the payload at 2026-08-16 16:50:40 UTC ran this text — its HTML is committed in the same commit; identical to HEAD; differs from fd1ce74f, the last commit before the build, in fs_selection_rule both->neighborhood)
#   chunk 'fs-grf-select': fence lines 1440-1491
# ----------------------------------------------------------

## ==== chunk 'fs-grf-select' ====
t0 <- proc.time()

fs_grf <- forestsearch(
  # --- Required data + roles ----------------------------------------------
  df.analysis         = actg_df,
  confounders.name    = confounders.name,
  outcome.name        = outcome.name,
  event.name          = outcome.name,
  treat.name          = treat.name,
  id.name             = id.name,
  is.RCT              = TRUE,
  seedit              = analysis_seed,   # also seeds the forest
  # --- GLM parameters -----------------------------------------------------
  outcome_type        = fs_outcome_type,
  effect_measure      = fs_effect_measure,
  adverse_outcome     = TRUE,
  # --- GRF-selection mode -------------------------------------------------
  subgroup_method     = "grf",
  grf_selection       = fs_grf_demo_selection,   # "tree" | "frontier"
  grf_select_statistic = fs_grf_select_statistic, # frontier only
  grf_depth           = fs_grf_demo_depth,
  dmin.grf            = fs_grf_demo_dmin,          # DR-score harm floor (RD scale)
  frac.tau            = 0.8,
  # --- Selection criteria -------------------------------------------------
  sg_focus            = fs_grf_sg_focus,          # -> GRF frontier/tree rule
  selection_rule      = "neighborhood",
  effect_neighborhood = 0.10,
  n.min               = 60L,
  # --- Tier-2 de-biased gate ----------------------------------------------
  mr_inference         = run_mr,
  mr_inference_args    = list(draws = mr_draws,
                             family_native_neighborhood = fs_grf_family_native_nbhd),
  # --- Output -------------------------------------------------------------
  details             = TRUE
  # NOTE: subgroup_method = "grf" bypasses LASSO, DINA and the consistency
  # search; use_lasso / use_dina / use_twostage / fs.splits / maxk and related
  # parameters do not apply in this mode.
)

timings$fs_grf_select <- (proc.time() - t0)["elapsed"]

if (!is.null(fs_grf$sg.harm) && length(fs_grf$sg.harm) > 0) {
  cat("GRF-selected subgroup (H):",
      fmt_sg(fs_grf$sg.harm), "\n")
  cat("Subgroup size:", sum(fs_grf$df.est$treat.recommend == 0),
      sprintf("(%.1f%% of ITT)\n",
              100 * mean(fs_grf$df.est$treat.recommend == 0)))
} else {
  cat("GRF selection: no subgroup met the criteria.\n")
}

