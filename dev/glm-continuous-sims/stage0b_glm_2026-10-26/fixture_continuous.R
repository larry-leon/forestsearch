# fixture_continuous.R -- the continuous (MD) fixture for Stage 0b.
#
# TRANSPLANTED (CC, 2026-10-26) from the committed ACTG175 continuous
# application: quarto/applications/actg175/analysis_actg175_continuous_intervals.qmd,
# chunks `data-prep` (data preparation, verbatim) and `anchor` (the
# forestsearch() parameter block, verbatim), with exactly these overrides for
# the per-replicate shape:
#   mr_inference = FALSE (mr_inference_args dropped), details = FALSE,
#   quiet = TRUE, plot.sg = FALSE, plot.grf = FALSE,
#   parallel_args = list(plan = "sequential").
# seedit = params$seed is resolved to the document's YAML params value 8316951.
# use_lasso / use_grf / use_twostage / consistency_method and every threshold
# are kept as the document has them.

fixture_continuous <- function() {
  if (!requireNamespace("speff2trial", quietly = TRUE))
    stop("speff2trial is required for the ACTG175 fixtures")
  params <- list(seed = 8316951)

  # ---- data-prep (verbatim) ----
  cont_vars <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80")
  bin_vars  <- c("hemo", "homo", "drugs", "race", "gender", "symptom")
  confounders.name <- c(cont_vars, bin_vars)

  actg_df <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))
  actg_df$id <- seq_len(nrow(actg_df))
  actg_df$treat <- ifelse(actg_df$arms == 1L, 1L, 0L)
  actg_df <- actg_df[!is.na(actg_df$cd420), ]
  actg_df$y_decline <- actg_df$cd40 - actg_df$cd420
  for (v in bin_vars) actg_df[[v]] <- as.factor(actg_df[[v]])

  N <- nrow(actg_df)
  stopifnot(N == 1083L)

  # ---- anchor (verbatim parameter block, overrides marked) ----
  args <- list(
    df.analysis      = actg_df,
    confounders.name = confounders.name,
    outcome.name     = "y_decline",
    treat.name       = "treat",
    id.name          = "id",
    outcome_type     = "continuous",
    effect_measure   = "MD",
    adverse_outcome  = TRUE,
    seedit           = params$seed,
    sg_focus         = "effMaxSG",
    effect_neighborhood = 0.20,
    selection_rule   = "neighborhood",
    consistency_method = "resample",
    effect.threshold       = 10,
    consistency.threshold  = 5,
    pconsistency.threshold = 0.90,
    use_twostage     = TRUE,
    conf.cont_jcuts  = list(age = 10, preanti = 10, wtkg = 10,
                            karnof = 10, cd40 = 10, cd80 = 10),
    cut_type         = "default",
    maxk             = 2L,
    n.min            = 60L,
    d0.min           = 10L,
    d1.min           = 10L,
    fs.splits        = 500L,
    use_lasso        = FALSE,
    use_grf          = FALSE,
    use_dina         = FALSE,
    is.RCT           = TRUE,
    parallel_args    = list(plan = "sequential"),   # override (document: same)
    details          = FALSE,                       # override (document: same)
    quiet            = TRUE,                        # override (document: same)
    mr_inference     = FALSE,                       # override (document: TRUE; mr_inference_args dropped)
    plot.sg          = FALSE,                       # override (document: absent)
    plot.grf         = FALSE                        # override (document: absent)
  )
  list(label = "continuous: ACTG175 MD, committed application (effMaxSG, eps = 0.20)",
       df = actg_df, args = args,
       source = "quarto/applications/actg175/analysis_actg175_continuous_intervals.qmd, chunks data-prep + anchor")
}
