# ---------------------------------------------------------------------------
# Table S2 of León & Anderson (2026a), DINA row (MR column): forestsearch(subgroup_method = 'dina') -- verbatim transplant.
# Source: forestsearch quarto/applications/gbsg/analysis_gbsg_survival_multimethod.qmd
# Commit: cf4d64328b6288ce2ec009f2c2fb15689de6a4cb (2026-08-15 21:12 -0700), the last
#   commit touching the .qmd before the Table S2 payload (built 2026-08-15 22:50:30 PDT,
#   Mac-Studio-3.local, forestsearch 0.2.0); render html committed in 43b051b6.
#   At HEAD the .qmd differs from cf4d6432 only at lines 144-145 (run_cv / run_loo
#   FALSE -> TRUE, 6ad5a724); no transplanted chunk below is affected.
# Chunks (name:fence lines, inclusive): fs-dina-select:1499-1544
# Chunk bodies are copied verbatim; only these header comments and the
# '## ---- chunk: <name> ----' separators were added. Transplanted 2026-10-28.
# ---------------------------------------------------------------------------

## ---- chunk: fs-dina-select (lines 1500-1543) ----
t0 <- proc.time()

fs_dina <- forestsearch(
  # --- Required data + roles ----------------------------------------------
  df.analysis,
  outcome.name        = outcome.name,
  event.name          = event.name,
  treat.name          = treat.name,
  id.name             = id.name,
  confounders.name    = confounders.name,
  df.predict          = NULL,            # default NULL
  df.test             = NULL,            # default NULL
  is.RCT              = TRUE,            # default TRUE
  seedit              = 8316951,         # default 8316951
  # --- DINA-selection mode ------------------------------------------------
  subgroup_method     = "dina",          # default "consistency"
  dina_args           = list(seed = dina_seed, n_folds = dina_nfolds),  # fit tuning only
  dina_select_statistic = fs_dina_select_statistic,  # setup-chunk selector (top-level)
  dina_res            = NULL,            # default NULL (reuse a fitted dina object)
  # --- Selection criteria (passed through to dina_subgroup) ---------------
  sg_focus            = fs_dina_sg_focus,  # set in setup chunk
  selection_rule      = "neighborhood",  # default "neighborhood"
  effect_neighborhood = 0.10,            # default 0.10
  n.min               = 60,              # default 60
  hr.threshold        = fs_dina_hr_threshold,  # set in setup chunk (-> log-HR floor)
  # --- MR ----------------------------------------------
  mr_inference         = run_mr,           # set in setup chunk
  mr_inference_args    = list(draws = mr_draws,
                             family_native_neighborhood = fs_dina_family_native_nbhd),
  # --- Output -------------------------------------------------------------
  details             = TRUE             # default FALSE
  # NOTE: subgroup_method = "dina" bypasses GRF, LASSO and the consistency
  # search; use_grf / use_lasso / use_twostage / fs.splits / maxk /
  # stop_threshold and related parameters do not apply in this mode.
)

timings$fs_dina_select <- (proc.time() - t0)["elapsed"]
stopifnot(!is.null(fs_dina), !is.null(fs_dina$args_call_all))

cat("DINA-selected subgroup (H):",
    paste(fs_dina$sg.harm, collapse = " & "), "\n")
cat("Subgroup size:", sum(fs_dina$df.est$treat.recommend == 0),
    sprintf("(%.1f%% of ITT)\n",
            100 * mean(fs_dina$df.est$treat.recommend == 0)))

