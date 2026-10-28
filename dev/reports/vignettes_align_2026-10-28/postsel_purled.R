## ----setup--------------------------------------------------------------------
#| message: false
#| warning: false
knitr::opts_chunk$set(comment = "#>", fig.align = "center")
library(forestsearch)
library(survival)


## ----data---------------------------------------------------------------------
df.analysis <- gbsg
df.analysis <- within(df.analysis, {
  id          <- seq_len(nrow(df.analysis))
  time_months <- rfstime / 30.4375
  grade3      <- ifelse(grade == "3", 1, 0)
})
confounders.name <- c("age", "meno", "size", "grade3", "nodes", "pgr", "er")
cat("Sample size:", nrow(df.analysis),
    "| events:", sum(df.analysis$status),
    sprintf("(%.1f%%)\n", 100 * mean(df.analysis$status)))


## ----fit----------------------------------------------------------------------
t0 <- proc.time()
fs <- forestsearch(
  df.analysis,
  outcome.name = "time_months", event.name = "status",
  treat.name = "hormon", id.name = "id",
  confounders.name = confounders.name,
  is.RCT = TRUE, seedit = 8316951, est.scale = "hr",
  use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE,
  conf.cont_jcuts = list(er = 10, pgr = 10), collapse_cuts = TRUE,
  conf_force = c("er <= 0", "pgr <= 0"),
  subgroup_method = "consistency", sg_focus = "effMaxSG",
  selection_rule = "neighborhood", effect_neighborhood = 0.20,
  mr_inference = TRUE,
  mr_inference_args = list(draws = 5000L, field_recovery = TRUE),
  n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
  hr.threshold = 1.0, hr.consistency = 1.0, pconsistency.threshold = 0.90,
  stop_threshold = NULL,
  details = FALSE, quiet = TRUE, plot.sg = FALSE,
  parallel_args = list(plan = "sequential"))
cat(sprintf("fit time: %.1f seconds\n", (proc.time() - t0)[["elapsed"]]))


## ----subgroup-----------------------------------------------------------------
cat("Selected subgroup:", paste(fs$sg.harm, collapse = " & "), "\n")


## ----print--------------------------------------------------------------------
print(fs)


## ----summary------------------------------------------------------------------
summary(fs)


## ----locations----------------------------------------------------------------
g <- fs$mr_inference
cat(sprintf("unadjusted HR on Hhat       %.3f\n", g$naive$est))
cat(sprintf("corrected (two-term) HR     %.3f\n", g$debiased$est))
cat(sprintf("field one-sided lower       %.3f\n", g$field$lower_1s))
cat(sprintf("field-s one-sided upper Hc  %.3f  (vs 0.80 / 0.85)\n",
            g$field$complement$upper_1s_s))


## ----phat---------------------------------------------------------------------
ph  <- g$reselection$p_hat
lab <- g$selected_label
cat(sprintf("p-hat(Hhat) = %.3f over a family of %d candidates\n",
            unname(ph[[lab]]), g$n_family))
print(round(sort(ph, decreasing = TRUE)[1:5], 3))


## ----recovery-----------------------------------------------------------------
rc <- g$field$recovery
cat(sprintf("sens_H  = %.3f   mean share of Hhat retained by a re-selection
", rc$sens_H))
cat(sprintf("ppv_H   = %.3f   mean share of a re-selection that lies in Hhat
", rc$ppv_H))
cat(sprintf("containment q10 / q50 / q90 = %.3f / %.3f / %.3f
", rc$q10, rc$q50, rc$q90))
cat(sprintf("share of draws containing all of Hhat = %.3f  (over %d draws)
",
            rc$share_equal_1, rc$n_used))


## ----session------------------------------------------------------------------
sessionInfo()

