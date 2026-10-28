# ---------------------------------------------------------------------------
# Table S2 of León & Anderson (2026a), FS row (MR column): main forestsearch() fit -- verbatim transplant.
# Source: forestsearch quarto/applications/gbsg/analysis_gbsg_survival_multimethod.qmd
# Commit: cf4d64328b6288ce2ec009f2c2fb15689de6a4cb (2026-08-15 21:12 -0700), the last
#   commit touching the .qmd before the Table S2 payload (built 2026-08-15 22:50:30 PDT,
#   Mac-Studio-3.local, forestsearch 0.2.0); render html committed in 43b051b6.
#   At HEAD the .qmd differs from cf4d6432 only at lines 144-145 (run_cv / run_loo
#   FALSE -> TRUE, 6ad5a724); no transplanted chunk below is affected.
# Chunks (name:fence lines, inclusive): forestsearch-main:687-809
# Chunk bodies are copied verbatim; only these header comments and the
# '## ---- chunk: <name> ----' separators were added. Transplanted 2026-10-28.
# ---------------------------------------------------------------------------

## ---- chunk: forestsearch-main (lines 688-808) ----
#| fig-width: 10
#| fig-height: 7

t0 <- proc.time()

fs <- forestsearch(
  # --- Required data + roles ----------------------------------------------
  df.analysis,
  outcome.name              = outcome.name,
  event.name                = event.name,
  treat.name                = treat.name,
  id.name                   = id.name,
  confounders.name          = confounders.name,
  potentialOutcome.name     = NULL,           # default NULL
  flag_harm.name            = NULL,           # default NULL
  df.predict                = NULL,           # default NULL
  df.test                   = NULL,           # default NULL
  is.RCT                    = TRUE,           # default TRUE
  seedit                    = 8316951,        # default 8316951
  est.scale                 = "hr",           # default "hr"
  # --- Screening engines (candidate generation) ---------------------------
  use_lasso                 = FALSE,           # default TRUE
  use_grf                   = FALSE,           # default TRUE
  grf_res                   = NULL,           # default NULL
  grf_cuts                  = NULL,           # default NULL
  grf_depth                 = fs_grf_depth,   # set in setup chunk (default 2)
  dmin.grf                  = fs_dmin_grf,     # set in setup chunk (default 0.0)
  frac.tau                  = 0.80,           # default 0.8
  vi.grf.min                = -0.2,           # default -0.2
  return_selected_cuts_only = TRUE,           # default TRUE
  tune_grf                  = FALSE,          # default FALSE
  use_dina                  = FALSE,          # default FALSE
  dina_res                  = NULL,           # default NULL
  dina_cuts                 = NULL,           # default NULL
  dina_args                 = fs_dina_args,   # extra DINA tuning (setup chunk)
  # --- Candidate-cut construction -----------------------------------------
  max_n_confounders         = 1000,           # default 1000
  conf_force                = fs_conf_force,           # default NULL
  defaultcut_names          = NULL,           # default NULL
  cut_type                  = "default",      # default "default"
  exclude_cuts              = NULL,           # default NULL
  replace_med_grf           = FALSE,          # default FALSE
  cont.cutoff               = 4,              # default 4
  conf.cont_medians         = NULL,           # default NULL
  conf.cont_medians_force   = NULL,           # default NULL
  conf.cont_jcuts           = list(er = 10, pgr = 10),           # default NULL
  collapse_cuts = TRUE,
  # --- Selection engine ---------------------------------------------------
  subgroup_method           = fs_subgroup_method,  # set in setup chunk: consistency | dina | grf
  grf_selection             = fs_grf_selection,    # set in setup chunk (used when subgroup_method = "grf")
  grf_select_statistic      = fs_grf_select_statistic,  # set in setup chunk (frontier only)
  dina_select_statistic     = fs_dina_select_statistic,  # set in setup chunk (used when subgroup_method = "dina")
  sg_focus                  = fs_sg_focus,         # set in setup chunk
  selection_rule            = "neighborhood",  # default "neighborhood"
  effect_neighborhood       = 0.10,            # default 0.10
  # --- MR ----------------------------------------------
  mr_inference               = run_mr,           # set in setup chunk
  mr_inference_args          = list(draws = mr_draws,
                                   family_native_neighborhood = fs_family_native_nbhd),
  # --- Subgroup constraints -----------------------------------------------
  n.min                     = 60,                       # default 60
  d0.min                    = 10,                       # default 10
  d1.min                    = 10,                       # default 10
  maxk                      = 2,                        # default 2
  # --- Thresholds (survival / HR scale) -----------------------------------
  hr.threshold              = fs_hr_threshold,    # set in setup chunk
  hr.consistency            = fs_hr_consistency,  # set in setup chunk
  effect.threshold          = NULL,           # default NULL (alias of hr.threshold)
  consistency.threshold     = NULL,           # default NULL (alias of hr.consistency)
  pconsistency.threshold    = 0.90,           # default 0.90
  consistency_method        = fs_consistency_method,
  m1.threshold              = Inf,            # default Inf
  minp                      = 0.025,          # default 0.025
  # --- Consistency search -------------------------------------------------
  fs.splits                 = 1000,           # default 1000
  stop_threshold            = NULL,           # NULL is what this configuration RESOLVES to, not a
                                              #   preference.  This call's sg_focus ("effMaxSG")
                                              #   normalizes to "hrMaxSG", and forestsearch_main.R:1403
                                              #   resets stop_threshold to NULL for that focus: a
                                              #   neighborhood band needs every candidate evaluated,
                                              #   so a truncated scan is unsound.  Any non-NULL value
                                              #   reaches the same NULL but warns on every render;
                                              #   passing NULL skips the reset branch outright (its
                                              #   guard is !is.null(stop_threshold)).  Named here so
                                              #   the resolved setting stays inspectable.
  use_twostage              = TRUE,           # default TRUE
  twostage_args             = list(),         # default list()
  max.minutes               = 3,              # default 3
  # --- Output / diagnostics -----------------------------------------------
  show_candidate_summary    = TRUE,           # default FALSE
  details                   = TRUE,           # default FALSE
  quiet                     = FALSE,          # default FALSE
  by.risk                   = 12,             # default 12
  plot.sg                   = TRUE,           # default FALSE
  plot.grf                  = FALSE,          # default FALSE
  # --- Parallel -----------------------------------------------------------
  parallel_args = list(
    plan = "multisession",
    workers = n_cores,
    show_message = TRUE
  )
)

plan("sequential")
timings$forestsearch <- (proc.time() - t0)["elapsed"]

# Fit-integrity guard.  NOT stopifnot(!is.null(fs$sg.harm)): a NULL subgroup is
# a legitimate outcome this document handles explicitly (see the *-no-subgroup
# chunks).  Assert only that the fit returned and captured its arguments, which
# every downstream provenance table depends on.
stopifnot(!is.null(fs), !is.null(fs$args_call_all))

# Guard for the GRF-dependent diagnostics below, read from the FIT rather than
# from the call, so it tracks whatever the fit actually resolved to.  With
# use_grf = FALSE no GRF cut is ever produced, and the bootstrap GRF-cut
# stability table would render a single "(no GRF cut) -- 100%" row presented as
# a result.
fs_used_grf <- isTRUE(fs$args_call_all$use_grf)

cat("\nForestSearch completed in",
    round(timings$forestsearch, 1), "seconds\n")

