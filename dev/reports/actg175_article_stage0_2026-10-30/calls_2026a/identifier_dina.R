# ---- transplant (verbatim chunk bodies; fences removed) ----
# source: forestsearch/quarto/applications/actg175/analysis_actg175_binary_multimethod_fixed_family.qmd
# forestsearch commit 43b051b6 (2026-08-16 12:07; the render that wrote the payload at 2026-08-16 16:50:40 UTC ran this text — its HTML is committed in the same commit; identical to HEAD; differs from fd1ce74f, the last commit before the build, in fs_selection_rule both->neighborhood)
#   chunk 'dina-config': fence lines 990-1001
#   chunk 'fs-dina-select': fence lines 1218-1267
# ----------------------------------------------------------

## ==== chunk 'dina-config' ====
#| code-fold: false

# Shared DINA *fit* settings (binary outcome => family mapped to "binomial").
# Selection parameters (sg_focus, or_threshold) live in the setup chunk under
# the fs_dina_* knobs so they can be tuned independently; (B) and (C) read them.
dina_seed   <- 8316951L
dina_nfolds <- 2L
# Harm/benefit floor on the natural-parameter (log-OR) link scale for the
# binomial family.
m_diff_dina <- log(fs_dina_or_threshold)

## ==== chunk 'fs-dina-select' ====
t0 <- proc.time()

fs_dina <- forestsearch(
  # --- Required data + roles ----------------------------------------------
  df.analysis         = actg_df,
  confounders.name    = confounders.name,
  outcome.name        = outcome.name,
  event.name          = outcome.name,
  treat.name          = treat.name,
  id.name             = id.name,
  is.RCT              = TRUE,
  seedit              = analysis_seed,
  # --- GLM parameters -----------------------------------------------------
  outcome_type        = fs_outcome_type,
  effect_measure      = fs_effect_measure,
  adverse_outcome     = TRUE,
  # --- DINA-selection mode ------------------------------------------------
  subgroup_method     = "dina",
  dina_args           = list(seed = dina_seed, n_folds = dina_nfolds),  # fit tuning only
  dina_select_statistic = fs_dina_select_statistic,
  # --- Selection criteria (passed through to dina_subgroup) ---------------
  sg_focus            = fs_dina_sg_focus,
  selection_rule      = "neighborhood",
  effect_neighborhood = 0.10,
  n.min               = 60L,
  effect.threshold    = fs_dina_or_threshold,  # -> log-OR harm floor
  # --- Tier-2 de-biased gate ----------------------------------------------
  mr_inference         = run_mr,
  mr_inference_args    = list(draws = mr_draws,
                             family_native_neighborhood = fs_dina_family_native_nbhd),
  # --- Output -------------------------------------------------------------
  details             = TRUE
  # NOTE: subgroup_method = "dina" bypasses GRF, LASSO and the consistency
  # search; use_grf / use_lasso / use_twostage / fs.splits / maxk and related
  # parameters do not apply in this mode.
)

timings$fs_dina_select <- (proc.time() - t0)["elapsed"]

if (!is.null(fs_dina$sg.harm) && length(fs_dina$sg.harm) > 0) {
  cat("DINA-selected subgroup (H):",
      fmt_sg(fs_dina$sg.harm), "\n")
  cat("Subgroup size:", sum(fs_dina$df.est$treat.recommend == 0),
      sprintf("(%.1f%% of ITT)\n",
              100 * mean(fs_dina$df.est$treat.recommend == 0)))
} else {
  cat("DINA selection: no subgroup met the criteria.\n")
}

