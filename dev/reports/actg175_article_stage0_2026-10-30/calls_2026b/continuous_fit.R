# ---- transplant (verbatim chunk bodies; fences removed) ----
# source: fs-glms-interpretable/quarto/actg175/analysis_actg175_continuous_mr.qmd
# fs-glms-interpretable HEAD 19c6ac41c1198dbc513f5c6978d0d6c5df1fc157 (file last changed d97ca1f / 935ae87, 2026-09-23)
#   chunk 'fit': fence lines 282-315
# ----------------------------------------------------------

## ==== chunk 'fit' ====
mr_draws <- 5000L
t0 <- proc.time()[["elapsed"]]
fit <- forestsearch(
  df.analysis = df, confounders.name = confounders,
  outcome.name = "y_decline", treat.name = "treat", id.name = "id",      # lines 95-97
  outcome_type = "continuous", effect_measure = "MD",                    # lines 98-99
  adverse_outcome = TRUE,                                                # line 100
  seedit = 8316951,                                                      # line 101 (params$seed)
  sg_focus = "effMaxSG", effect_neighborhood = 0.20,                     # THIS WORKSTREAM'S DEPARTURE (§1.1)
  selection_rule = "neighborhood",                                       # line 103
  consistency_method = "resample",                                       # line 104
  effect.threshold = 10, consistency.threshold = 5,                      # lines 105-106
  pconsistency.threshold = 0.90,                                         # line 107
  use_twostage = TRUE,                                                   # line 108
  conf.cont_jcuts = list(age = 10, preanti = 10, wtkg = 10,              # lines 109-110
                         karnof = 10, cd40 = 10, cd80 = 10),
  cut_type = "default", maxk = 2L, n.min = 60L,                          # lines 111-113
  d0.min = 10L, d1.min = 10L, fs.splits = 500L,                          # lines 114-116
  use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE, is.RCT = TRUE,   # lines 117-120
  mr_inference = TRUE,                                                   # line 128
  mr_inference_args = list(ci_method = "field", draws = mr_draws,        # lines 129-134
                           include_complement = TRUE, field_complement = TRUE,
                           field_scale_complement = "selected",
                           return_reselection = TRUE,
                           keep_declaration_field = TRUE, keep_field_matrix = FALSE,
                           declaration_c0 = c(-5, -10, -15)))
fit_sec <- proc.time()[["elapsed"]] - t0
if (fit_sec > 900) stop(sprintf("the fit took %.0f s, over the 15-minute cap", fit_sec))
g <- fit$mr_inference
cat(sprintf("forestsearch %s; fit including MR: %.1f s (multiplier stage %.1f s; field blocks %.1f + %.1f s)\n",
            packageVersion("forestsearch"), fit_sec, g$timing_seconds,
            g$field$timing_seconds, g$field$complement$timing_seconds))

