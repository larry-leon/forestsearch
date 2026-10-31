# ---- transplant (verbatim chunk bodies; fences removed) ----
# source: fs-glms-interpretable/quarto/actg175/analysis_actg175_binary_mr.qmd
# fs-glms-interpretable HEAD 19c6ac41c1198dbc513f5c6978d0d6c5df1fc157 (file last changed d97ca1f / 935ae87, 2026-09-23)
#   chunk 'fit': fence lines 264-297
# ----------------------------------------------------------

## ==== chunk 'fit' ====
mr_draws <- 5000L
t0 <- proc.time()[["elapsed"]]
fit <- forestsearch(
  df.analysis = df,
  outcome.name = "y_binary_neg", event.name = "y_binary_neg", treat.name = "treat", id.name = "id",
  confounders.name = confounders,
  is.RCT = TRUE, seedit = 8316951,
  outcome_type = "binary", effect_measure = "OR", adverse_outcome = TRUE,
  use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE,
  max_n_confounders = 1000,
  cut_type = "default", cont.cutoff = 4, conf.cont_jcuts = list(cd40 = 10, wtkg = 10),
  collapse_cuts = TRUE,
  subgroup_method = "consistency", sg_focus = "effMaxSG",
  selection_rule = "neighborhood", effect_neighborhood = 0.20,
  n.min = 60, d0.min = 10, d1.min = 10, maxk = 2, max_subgroups_search = Inf,
  effect.threshold = 1.00, consistency.threshold = 1.00, pconsistency.threshold = 0.90,
  consistency_method = "resample", m1.threshold = Inf, minp = 0.025, fs.splits = 1000,
  stop_threshold = NULL,
  mr_inference = TRUE,
  mr_inference_args = list(ci_method = "field", draws = mr_draws, include_complement = TRUE,
                           confirm_rule = "point", field_uniform = FALSE,
                           field_complement = TRUE, ij_residual = "two_term",
                           field_decompose = TRUE, field_scale_complement = "selected",
                           field_recovery = FALSE, return_reselection = TRUE,
                           keep_declaration_field = TRUE, keep_field_matrix = FALSE,
                           declaration_c0 = c(0.70, 0.75, 0.80, 0.85)))
fit_sec <- proc.time()[["elapsed"]] - t0
if (fit_sec > 900) stop(sprintf("the fit took %.0f s, over the 15-minute cap", fit_sec))
g <- fit$mr_inference
cat(sprintf("forestsearch %s; fit including MR: %.1f s (multiplier stage %.1f s; field blocks %.1f + %.1f s)\n",
            packageVersion("forestsearch"), fit_sec, g$timing_seconds,
            g$field$timing_seconds, g$field$complement$timing_seconds))

