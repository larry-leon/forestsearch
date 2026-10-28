## ----setup--------------------------------------------------------------------
#| message: false
#| warning: false
knitr::opts_chunk$set(comment = "#>", fig.align = "center")
library(forestsearch)
library(survival)


## ----data---------------------------------------------------------------------
df <- survival::gbsg
df <- within(df, {
  id          <- seq_len(nrow(df))
  time_months <- rfstime / 30.4375
  grade3      <- ifelse(grade == "3", 1, 0)
})
confounders <- c("age", "meno", "size", "grade3", "nodes", "pgr", "er")
stopifnot(nrow(df) == 686, sum(df$status) == 299)
whole <- coxph(Surv(time_months, status) ~ hormon, data = df)
sprintf("Whole trial: %d patients, %d events, hazard ratio %.3f (%.3f, %.3f)",
        nrow(df), sum(df$status), exp(coef(whole)), exp(confint(whole))[1],
        exp(confint(whole))[2])


## ----fit-fs-------------------------------------------------------------------
#| message: false
#| results: hide
t0 <- proc.time()[["elapsed"]]
fit_fs <- forestsearch(
  df.analysis = df,
  outcome.name = "time_months", event.name = "status", treat.name = "hormon", id.name = "id",
  confounders.name = confounders,
  is.RCT = TRUE, est.scale = "hr", seedit = 8316951,
  conf.cont_jcuts = list(er = 10, pgr = 10), conf_force = c("er <= 0", "pgr <= 0"),
  n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
  subgroup_method = "consistency", consistency_method = "resample",
  hr.threshold = 1.00, hr.consistency = 1.00, pconsistency.threshold = 0.90,
  stop_threshold = NULL,
  sg_focus = "effMaxSG", selection_rule = "neighborhood", effect_neighborhood = 0.10,
  mr_inference = TRUE,
  mr_inference_args = list(ci_method = "field", draws = 5000L, field_recovery = TRUE),
  quiet = TRUE, parallel_args = list(plan = "sequential"))
sec_fs <- proc.time()[["elapsed"]] - t0


## ----fit-dina-----------------------------------------------------------------
#| message: false
#| results: hide
t0 <- proc.time()[["elapsed"]]
fit_dina <- forestsearch(
  df.analysis = df,
  outcome.name = "time_months", event.name = "status", treat.name = "hormon", id.name = "id",
  confounders.name = confounders,
  is.RCT = TRUE, est.scale = "hr", seedit = 8316951,
  subgroup_method = "dina", dina_select_statistic = "effect",
  dina_args = list(seed = 8316951, n_folds = 2),
  hr.threshold = 0.90, n.min = 60,
  sg_focus = "effMaxSG", selection_rule = "neighborhood", effect_neighborhood = 0.10,
  mr_inference = TRUE, mr_inference_args = list(ci_method = "ij", draws = 5000L),
  quiet = TRUE, parallel_args = list(plan = "sequential"))
sec_dina <- proc.time()[["elapsed"]] - t0


## ----fit-grf------------------------------------------------------------------
#| message: false
#| results: hide
t0 <- proc.time()[["elapsed"]]
fit_grf <- forestsearch(
  df.analysis = df,
  outcome.name = "time_months", event.name = "status", treat.name = "hormon", id.name = "id",
  confounders.name = confounders,
  is.RCT = TRUE, est.scale = "hr", seedit = 8316951,
  subgroup_method = "grf", grf_selection = "frontier", grf_select_statistic = "effect",
  grf_depth = 2L, dmin.grf = 0.0, frac.tau = 0.8, n.min = 60,
  sg_focus = "effMaxSG", selection_rule = "neighborhood", effect_neighborhood = 0.10,
  mr_inference = TRUE, mr_inference_args = list(ci_method = "ij", draws = 5000L),
  quiet = TRUE, parallel_args = list(plan = "sequential"))
sec_grf <- proc.time()[["elapsed"]] - t0


## ----selected-----------------------------------------------------------------
fits <- list(FS = fit_fs, DINA = fit_dina, GRF = fit_grf)
secs <- c(FS = sec_fs, DINA = sec_dina, GRF = sec_grf)
for (m in names(fits)) {
  gm <- fits[[m]]$mr_inference
  cat(sprintf("%-4s %-28s N %3d   fit %5.1f s (multiplier stage %.1f s)   family %d\n",
              m, paste(fits[[m]]$sg.harm, collapse = " & "), gm$n_selected, secs[[m]],
              gm$timing_seconds, gm$n_family))
}


## ----mr-fs--------------------------------------------------------------------
g <- fit_fs$mr_inference
cat(sprintf("Selection bias (log scale): %.3f over a family of %d candidates\n",
            g$selection_bias, g$n_family))
cat(sprintf("Region:     naive %.2f (%.2f, %.2f)  ->  corrected %.2f (%.2f, %.2f)\n",
            g$naive$est, g$naive$lower, g$naive$upper,
            g$debiased$est, g$debiased$lower, g$debiased$upper))
cc <- g$complement
cat(sprintf("Complement: naive %.2f (%.2f, %.2f)  ->  corrected %.2f (%.2f, %.2f)\n",
            cc$naive$est, cc$naive$lower, cc$naive$upper,
            cc$debiased$est, cc$debiased$lower, cc$debiased$upper))


## ----bootstrap----------------------------------------------------------------
#| message: false
#| warning: false
B <- c(FS = 25L, DINA = 25L, GRF = 10L)
boots <- list(); boot_secs <- c()
for (m in names(fits)) {
  t0 <- proc.time()[["elapsed"]]
  boots[[m]] <- forestsearch_bootstrap_dofuture(
    fs.est = fits[[m]], nb_boots = B[[m]], seed = 8316951L,
    mr_in_replicates = FALSE, parallel_args = list(plan = "sequential"))
  boot_secs[m] <- proc.time()[["elapsed"]] - t0
}
round(boot_secs, 1)


## ----boot-summary-------------------------------------------------------------
summ_fs <- summarize_bootstrap_results(sgharm = fit_fs$sg.harm, boot_results = boots$FS,
                                       create_plots = FALSE, est.scale = "hr")
summ_fs$table


## ----mr-table-----------------------------------------------------------------
mr_estimates_table(fit_fs, boot_results = boots$FS, est.scale = "hr", digits = 2)


## ----three-columns------------------------------------------------------------
three_columns <- function(method, fit, boot) {
  g <- fit$mr_inference; cc <- g$complement
  n_H <- g$n_selected; n_Hc <- cc$n
  f2 <- function(e, l, u) sprintf("%.2f (%.2f, %.2f)", e, l, u)
  data.frame(
    method = method, region = c("H", "Hc"),
    n = c(sprintf("%d (%d%%)", n_H, round(100 * n_H / (n_H + n_Hc))),
          sprintf("%d (%d%%)", n_Hc, round(100 * n_Hc / (n_H + n_Hc)))),
    naive = c(f2(g$naive$est, g$naive$lower, g$naive$upper),
              f2(cc$naive$est, cc$naive$lower, cc$naive$upper)),
    bootstrap = c(f2(boot$H_estimates$H2, boot$H_estimates$H2_lower,
                     boot$H_estimates$H2_upper),
                  f2(boot$Hc_estimates$H2, boot$Hc_estimates$H2_lower,
                     boot$Hc_estimates$H2_upper)),
    multiplier = c(f2(g$debiased$est, g$debiased$lower, g$debiased$upper),
                   f2(cc$debiased$est, cc$debiased$lower, cc$debiased$upper)))
}
side_by_side <- rbind(three_columns("FS", fit_fs, boots$FS),
                      three_columns("DINA", fit_dina, boots$DINA),
                      three_columns("GRF", fit_grf, boots$GRF))
boot_label <- sprintf("Full bootstrap (B = %d / %d / %d)", B[["FS"]], B[["DINA"]], B[["GRF"]])
knitr::kable(side_by_side, row.names = FALSE, align = "llrrrr",
             col.names = c("Method", "SG", "n (%)", "Naive", boot_label,
                           "Multiplier (5,000 draws)"),
             caption = paste("This vignette's run: hazard ratios with 95% intervals.",
                             "The bootstrap column is an illustration at a small B."))


## ----table-s2-----------------------------------------------------------------
table_s2 <- data.frame(
  method = c("FS", "FS", "DINA", "DINA", "GRF", "GRF"),
  region = rep(c("H", "Hc"), 3),
  n = c("61 (9%)", "625 (91%)", "89 (13%)", "597 (87%)", "82 (12%)", "604 (88%)"),
  naive = c("2.54 (1.25, 5.17)", "0.61 (0.47, 0.79)", "1.60 (0.88, 2.91)",
            "0.63 (0.48, 0.82)", "1.95 (1.04, 3.67)", "0.61 (0.47, 0.80)"),
  bootstrap = c("1.94 (0.89, 4.21)", "0.63 (0.43, 0.93)", "0.89 (0.38, 2.07)",
                "0.66 (0.40, 1.11)", "1.04 (0.49, 2.24)", "0.66 (0.39, 1.13)"),
  multiplier = c("1.44 (0.62, 3.36)", "0.64 (0.40, 1.03)", "1.11 (0.45, 2.75)",
                 "0.66 (0.40, 1.09)", "1.05 (0.47, 2.34)", "0.66 (0.41, 1.06)"))
knitr::kable(table_s2, row.names = FALSE, align = "llrrrr",
             col.names = c("Method", "SG", "n (%)", "Naive", "Full bootstrap (B = 1,000)",
                           "Multiplier"),
             caption = "Table S2 of León and Anderson (2026a), as published.")


## ----reading-values-----------------------------------------------------------
#| echo: false
r2 <- function(x) sprintf("%.2f", x)
gd <- fit_dina$mr_inference; gg <- fit_grf$mr_inference


## ----cost---------------------------------------------------------------------
cost <- data.frame(
  method = names(fits),
  mr_seconds = sapply(fits, function(f) f$mr_inference$timing_seconds),
  bootstrap_B = B[names(fits)],
  bootstrap_seconds = boot_secs[names(fits)])
cost$seconds_per_replicate <- cost$bootstrap_seconds / cost$bootstrap_B
cost$projected_B_1000 <- round(cost$seconds_per_replicate * 1000)
knitr::kable(cost, row.names = FALSE, digits = 2,
             col.names = c("Method", "MR (s)", "B", "bootstrap (s)", "s per replicate",
                           "projected at B = 1,000 (s)"),
             caption = "Single-core timings from this run. The projection is linear in B.")


## ----phat---------------------------------------------------------------------
ph  <- g$reselection$p_hat
cat(sprintf("p-hat(H) = %.3f over a family of %d candidates\n",
            unname(ph[[g$selected_label]]), g$n_family))
print(round(sort(ph, decreasing = TRUE)[1:5], 3))


## ----recovery-----------------------------------------------------------------
rc <- g$field$recovery
cat(sprintf("sens_H = %.3f   mean share of H retained by a re-selection\n", rc$sens_H))
cat(sprintf("ppv_H  = %.3f   mean share of a re-selection that lies in H\n", rc$ppv_H))
cat(sprintf("containment q10 / q50 / q90 = %.3f / %.3f / %.3f\n", rc$q10, rc$q50, rc$q90))
cat(sprintf("share of draws containing all of H = %.3f  (over %d draws)\n",
            rc$share_equal_1, rc$n_used))


## ----session------------------------------------------------------------------
sessionInfo()

