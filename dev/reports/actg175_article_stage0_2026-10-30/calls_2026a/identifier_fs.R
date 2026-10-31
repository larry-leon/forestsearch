# ---- transplant (verbatim chunk bodies; fences removed) ----
# source: forestsearch/quarto/applications/actg175/analysis_actg175_binary_multimethod_fixed_family.qmd
# forestsearch commit 43b051b6 (2026-08-16 12:07; the render that wrote the payload at 2026-08-16 16:50:40 UTC ran this text — its HTML is committed in the same commit; identical to HEAD; differs from fd1ce74f, the last commit before the build, in fs_selection_rule both->neighborhood)
#   chunk 'forestsearch-main': fence lines 551-629
# ----------------------------------------------------------

## ==== chunk 'forestsearch-main' ====
#| fig-width: 10
#| fig-height: 7

t0 <- proc.time()

fs <- forestsearch(
  # --- Required data + roles ----------------------------------------------
  df.analysis               = actg_df,
  confounders.name          = confounders.name,
  outcome.name              = outcome.name,
  event.name                = outcome.name,    # binary: outcome IS the event
  treat.name                = treat.name,
  id.name                   = id.name,
  is.RCT                    = TRUE,
  seedit                    = analysis_seed,
  # --- GLM parameters -----------------------------------------------------
  outcome_type              = fs_outcome_type,
  effect_measure            = fs_effect_measure,
  adverse_outcome           = TRUE,
  # --- Screening engines (candidate generation) ---------------------------
  use_lasso                 = FALSE,
  use_grf                   = FALSE,
  grf_depth                 = fs_grf_depth,
  dmin.grf                  = fs_dmin_grf,
  frac.tau                  = 0.80,
  vi.grf.min                = -0.2,
  return_selected_cuts_only = TRUE,
  tune_grf                  = FALSE,
  use_dina                  = FALSE,
  dina_args                 = fs_dina_args,
  # --- Candidate-cut construction -----------------------------------------
  cut_type                  = "default",
  cont.cutoff               = 4L,
  conf.cont_jcuts           = list(cd40 = 10, wtkg = 10),
  # --- Selection engine ---------------------------------------------------
  subgroup_method           = fs_subgroup_method,        # consistency | dina | grf
  grf_selection             = fs_grf_selection,          # used when subgroup_method = "grf"
  grf_select_statistic      = fs_grf_select_statistic,   # frontier only
  dina_select_statistic     = fs_dina_select_statistic,  # used when subgroup_method = "dina"
  sg_focus                  = fs_sg_focus,
  selection_rule            = fs_selection_rule,
  effect_neighborhood       = 0.10,
  # --- Tier-2 de-biased gate ----------------------------------------------
  mr_inference               = run_mr,
  mr_inference_args          = list(draws = mr_draws,
                                   family_native_neighborhood = fs_family_native_nbhd),
  # --- Subgroup constraints -----------------------------------------------
  n.min                     = 60L,
  d0.min                    = 10L,
  d1.min                    = 10L,
  maxk                      = 2L,
  # --- Thresholds (OR scale) ----------------------------------------------
  effect.threshold          = fs_or_threshold,
  consistency.threshold     = fs_or_consistency,
  pconsistency.threshold    = fs_pconsistency,
  consistency_method        = fs_consistency_method,   # "resample" (gate requires it)
  # --- Consistency search -------------------------------------------------
  fs.splits                 = 500L,
  use_twostage              = TRUE,
  # --- Output / diagnostics -----------------------------------------------
  show_candidate_summary    = TRUE,
  details                   = TRUE,
  quiet                     = FALSE,
  plot.sg                   = FALSE,   # KM not applicable; use plot_sg_glm_outcomes()
  # --- Parallel -----------------------------------------------------------
  parallel_args = list(
    plan = "multisession",
    workers = n_cores,
    show_message = TRUE
  )
)

plan("sequential"); gc(verbose = FALSE)
timings$forestsearch <- (proc.time() - t0)["elapsed"]

cat("\nForestSearch completed in",
    round(timings$forestsearch, 1), "seconds\n")

