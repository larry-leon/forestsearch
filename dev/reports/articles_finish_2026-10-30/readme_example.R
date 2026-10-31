library(forestsearch)
library(survival)

df <- within(gbsg, {
  id          <- seq_len(nrow(gbsg))
  time_months <- rfstime / 30.4375
  grade3      <- ifelse(grade == "3", 1, 0)
})

fit <- forestsearch(
  df.analysis = df,
  outcome.name = "time_months", event.name = "status",
  treat.name = "hormon", id.name = "id",
  confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
  is.RCT = TRUE, est.scale = "hr", seedit = 8316951,
  # the candidate family
  conf.cont_jcuts = list(er = 10, pgr = 10), conf_force = c("er <= 0", "pgr <= 0"),
  n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
  # the thresholds and the selection rule
  hr.threshold = 1.00, hr.consistency = 1.00, pconsistency.threshold = 0.90,
  stop_threshold = NULL,
  sg_focus = "effMaxSG", selection_rule = "neighborhood", effect_neighborhood = 0.20,
  # selection-adjusted inference
  mr_inference = TRUE,
  mr_inference_args = list(ci_method = "field", draws = 5000L, include_complement = TRUE,
                           field_scale_complement = "selected", field_decompose = TRUE,
                           keep_declaration_field = TRUE,
                           declaration_c0 = c(0.70, 0.75, 0.80, 0.85)),
  parallel_args = list(plan = "sequential")
)

print(fit)
summary(fit)

# Would a screen protecting a uniform benefit (HR 0.70 to 0.85) have declared it?
fs_declaration_calibration(fit, alpha = 0.10, c0 = c(0.70, 0.75, 0.80, 0.85))
cat("REGION", paste(fit$sg.harm, collapse = " & "), "\n")
cat(sprintf("CHECK corrected %.3f field_lower %.3f complement %.3f field_s %.3f\n",
            fit$mr_inference$debiased$est, fit$mr_inference$field$lower_1s,
            fit$mr_inference$complement$debiased$est, fit$mr_inference$field$complement$upper_1s_s))
