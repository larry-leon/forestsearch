# fixture_binary.R -- the binary (OR) fixture for Stage 0b.
#
# TRANSPLANTED (CC, 2026-10-26) from the committed ACTG175 binary application:
# quarto/applications/actg175/analysis_actg175_binary_multimethod_fixed_family.qmd,
# chunks `parallel-setup` (the fs_* settings the call references),
# `analysis-roles`, `data-treat`, `data-outcomes`, `data-covariates` (data
# preparation, verbatim) and `forestsearch-main` (the forestsearch() parameter
# block, verbatim), with exactly these overrides for the per-replicate shape:
#   mr_inference = FALSE (mr_inference_args dropped), details = FALSE,
#   quiet = TRUE, plot.sg = FALSE, plot.grf = FALSE,
#   parallel_args = list(plan = "sequential").
# use_lasso / use_grf / vi.grf.min / use_twostage / consistency_method and
# every threshold are kept as the document has them.

fixture_binary <- function() {
  if (!requireNamespace("speff2trial", quietly = TRUE))
    stop("speff2trial is required for the ACTG175 fixtures")

  # ---- parallel-setup: settings referenced by the call (verbatim) ----
  fs_sg_focus              <- "effMaxSG"   # GLM-natural alias of "hrMaxSG"
  fs_selection_rule        <- "neighborhood"       # neighborhood ∩ Pareto (effMaxSG/effMinSG)
  fs_or_threshold          <- 1.0
  fs_or_consistency        <- 1.0
  fs_pconsistency          <- 0.90
  fs_consistency_method    <- "resample"
  fs_subgroup_method <- "consistency"  # "consistency" | "dina" | "grf"
  fs_grf_selection   <- "frontier"     # only if subgroup_method == "grf": "tree" | "frontier"
  fs_grf_select_statistic <- "effect"  # only if subgroup_method=="grf" & grf_selection=="frontier": "dr" | "effect".
  fs_grf_depth       <- 2L             # GRF depth (screening + GRF selection)
  fs_dmin_grf        <- 0              # GRF DR-score screening floor (RD scale for binary; 0 = permissive)
  fs_dina_args       <- list()         # extra args forwarded to DINA selection
  fs_dina_select_statistic <- "effect" # only if subgroup_method=="dina": "dina" | "effect".

  # ---- analysis-roles (verbatim) ----
  analysis_seed <- 8316951L
  set.seed(analysis_seed)
  cont_vars <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80")
  bin_vars  <- c("hemo", "homo", "drugs", "race", "gender", "symptom")
  confounders.name <- c(cont_vars, bin_vars)
  id.name   <- "id"
  treat.name   <- "treat"
  outcome.name <- "y_binary_neg"   # FS adverse-scale outcome (1 = no improvement)
  fs_outcome_type   <- "binary"
  fs_effect_measure <- "OR"

  # ---- data-treat / data-outcomes / data-covariates (verbatim) ----
  actg_df <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))
  actg_df$id <- seq_len(nrow(actg_df))
  actg_df$treat <- ifelse(actg_df$arms == 1L, 1L, 0L)
  actg_df <- actg_df[!is.na(actg_df$cd420), ]
  actg_df$y_binary_pos <- as.integer(actg_df$cd420 > actg_df$cd40)
  actg_df$y_binary_neg <- 1L - actg_df$y_binary_pos
  for (v in bin_vars) actg_df[[v]] <- as.numeric(actg_df[[v]])

  # ---- forestsearch-main (verbatim parameter block, overrides marked) ----
  args <- list(
    df.analysis               = actg_df,
    confounders.name          = confounders.name,
    outcome.name              = outcome.name,
    event.name                = outcome.name,    # binary: outcome IS the event
    treat.name                = treat.name,
    id.name                   = id.name,
    is.RCT                    = TRUE,
    seedit                    = analysis_seed,
    outcome_type              = fs_outcome_type,
    effect_measure            = fs_effect_measure,
    adverse_outcome           = TRUE,
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
    cut_type                  = "default",
    cont.cutoff               = 4L,
    conf.cont_jcuts           = list(cd40 = 10, wtkg = 10),
    subgroup_method           = fs_subgroup_method,        # consistency | dina | grf
    grf_selection             = fs_grf_selection,          # used when subgroup_method = "grf"
    grf_select_statistic      = fs_grf_select_statistic,   # frontier only
    dina_select_statistic     = fs_dina_select_statistic,  # used when subgroup_method = "dina"
    sg_focus                  = fs_sg_focus,
    selection_rule            = fs_selection_rule,
    effect_neighborhood       = 0.10,
    mr_inference              = FALSE,                     # override (document: run_mr = TRUE; mr_inference_args dropped)
    n.min                     = 60L,
    d0.min                    = 10L,
    d1.min                    = 10L,
    maxk                      = 2L,
    effect.threshold          = fs_or_threshold,
    consistency.threshold     = fs_or_consistency,
    pconsistency.threshold    = fs_pconsistency,
    consistency_method        = fs_consistency_method,   # "resample" (gate requires it)
    fs.splits                 = 500L,
    use_twostage              = TRUE,
    show_candidate_summary    = TRUE,
    details                   = FALSE,                     # override (document: TRUE)
    quiet                     = TRUE,                      # override (document: FALSE)
    plot.sg                   = FALSE,   # KM not applicable; use plot_sg_glm_outcomes()
    plot.grf                  = FALSE,                     # override (document: absent)
    parallel_args             = list(plan = "sequential")  # override (document: multisession, n_cores)
  )
  list(label = "binary: ACTG175 OR, committed application (fixed family, effMaxSG, eps = 0.10)",
       df = actg_df, args = args,
       source = "quarto/applications/actg175/analysis_actg175_binary_multimethod_fixed_family.qmd, chunks parallel-setup + analysis-roles + data-* + forestsearch-main")
}
