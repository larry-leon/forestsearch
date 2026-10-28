# ---------------------------------------------------------------------------
# Table S2 of León & Anderson (2026a), GRF row (MR column): forestsearch(subgroup_method = 'grf') -- verbatim transplant.
# Source: forestsearch quarto/applications/gbsg/analysis_gbsg_survival_multimethod.qmd
# Commit: cf4d64328b6288ce2ec009f2c2fb15689de6a4cb (2026-08-15 21:12 -0700), the last
#   commit touching the .qmd before the Table S2 payload (built 2026-08-15 22:50:30 PDT,
#   Mac-Studio-3.local, forestsearch 0.2.0); render html committed in 43b051b6.
#   At HEAD the .qmd differs from cf4d6432 only at lines 144-145 (run_cv / run_loo
#   FALSE -> TRUE, 6ad5a724); no transplanted chunk below is affected.
# Chunks (name:fence lines, inclusive): fs-grf-select:1740-1787
# Chunk bodies are copied verbatim; only these header comments and the
# '## ---- chunk: <name> ----' separators were added. Transplanted 2026-10-28.
# ---------------------------------------------------------------------------

## ---- chunk: fs-grf-select (lines 1741-1786) ----
t0 <- proc.time()

fs_grf <- forestsearch(
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
  seedit              = 8316951,         # default 8316951 (also seeds the forest)
  # --- GRF-selection mode -------------------------------------------------
  subgroup_method     = "grf",           # default "consistency"
  grf_selection       = fs_grf_demo_selection,  # set in setup chunk ("tree" | "frontier")
  grf_select_statistic = fs_grf_select_statistic,  # set in setup chunk (frontier only)
  grf_depth           = fs_grf_demo_depth,      # set in setup chunk (tree / DR-pair depth)
  dmin.grf            = fs_grf_demo_dmin,        # set in setup chunk (DR-score harm floor)
  frac.tau            = 0.8,             # default 0.8 (RMST truncation fraction)
  grf_res             = NULL,            # default NULL (reuse a fitted forest)
  # --- Selection criteria -------------------------------------------------
  sg_focus            = fs_grf_sg_focus,  # set in setup chunk (-> GRF frontier/tree rule)
  selection_rule      = "neighborhood",  # default "neighborhood"
  effect_neighborhood = 0.10,            # default 0.10
  n.min               = 60,              # default 60
  # --- MR ----------------------------------------------
  mr_inference         = run_mr,           # set in setup chunk
  mr_inference_args    = list(draws = mr_draws,
                             family_native_neighborhood = fs_grf_family_native_nbhd),
  # --- Output -------------------------------------------------------------
  details             = TRUE             # default FALSE
  # NOTE: subgroup_method = "grf" bypasses LASSO, DINA and the consistency
  # search; use_lasso / use_dina / use_twostage / fs.splits / maxk /
  # hr.threshold / pconsistency.threshold do not apply in this mode.
)

timings$fs_grf_select <- (proc.time() - t0)["elapsed"]
stopifnot(!is.null(fs_grf), !is.null(fs_grf$args_call_all))

cat("GRF-selected subgroup (H):",
    paste(fs_grf$sg.harm, collapse = " & "), "\n")
cat("Subgroup size:", sum(fs_grf$df.est$treat.recommend == 0),
    sprintf("(%.1f%% of ITT)\n",
            100 * mean(fs_grf$df.est$treat.recommend == 0)))

