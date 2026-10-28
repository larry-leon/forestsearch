## ----setup--------------------------------------------------------------------
#| message: false
#| warning: false
knitr::opts_chunk$set(comment = "#>", fig.align = "center")
library(forestsearch)
library(survival)


## ----data---------------------------------------------------------------------
df <- survival::gbsg
df <- within(df, {
  id          <- seq_len(nrow(df))              # the fit's subject identifier
  time_months <- rfstime / 30.4375              # recurrence-free survival, days to months
  grade3      <- ifelse(grade == "3", 1, 0)     # grade 3 against grade 1 or 2
})
confounders <- c("age", "meno", "size", "grade3", "nodes", "pgr", "er")
stopifnot(nrow(df) == 686, sum(df$status) == 299)

whole <- coxph(Surv(time_months, status) ~ hormon, data = df)
whole_hr <- c(est = unname(exp(coef(whole))), exp(confint(whole))[1, ])
c(patients = nrow(df), hormonal_therapy = sum(df$hormon == 1),
  none = sum(df$hormon == 0), events = sum(df$status))
sprintf("Whole-trial hazard ratio %.3f (%.3f, %.3f), hormonal therapy vs none",
        whole_hr[["est"]], whole_hr[[2]], whole_hr[[3]])


## ----fit-fs-------------------------------------------------------------------
#| message: false
#| results: hide
t0 <- proc.time()[["elapsed"]]
fit_fs <- forestsearch(
  df.analysis = df,
  outcome.name = "time_months", event.name = "status", treat.name = "hormon", id.name = "id",
  confounders.name = confounders,
  is.RCT = TRUE, est.scale = "hr", seedit = 8316951,
  # the candidate family
  conf.cont_jcuts = list(er = 10, pgr = 10), conf_force = c("er <= 0", "pgr <= 0"),
  n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
  # the identifier, its thresholds and the selection rule
  subgroup_method = "consistency", consistency_method = "resample",
  hr.threshold = 1.00, hr.consistency = 1.00, pconsistency.threshold = 0.90,
  stop_threshold = NULL,                      # evaluate every candidate; no early stop
  sg_focus = "effMaxSG", selection_rule = "neighborhood", effect_neighborhood = 0.20,
  # selection-adjusted inference
  mr_inference = TRUE,
  mr_inference_args = list(ci_method = "field", draws = 5000L, include_complement = TRUE,
                           field_scale_complement = "selected", field_decompose = TRUE,
                           keep_declaration_field = TRUE,
                           declaration_c0 = c(0.70, 0.75, 0.80, 0.85)),
  quiet = TRUE, parallel_args = list(plan = "sequential"))
fit_sec_fs <- proc.time()[["elapsed"]] - t0


## ----fs-selected--------------------------------------------------------------
g_fs <- fit_fs$mr_inference
cat(sprintf("Identified subgroup: %s\n", paste(fit_fs$sg.harm, collapse = " & ")))
cat(sprintf("Fit including inference: %.1f s\n", fit_sec_fs))
cat(sprintf("  multiplier stage %.1f s; field blocks %.1f + %.1f s\n", g_fs$timing_seconds,
            g_fs$field$timing_seconds, g_fs$field$complement$timing_seconds))


## ----fs-funnel----------------------------------------------------------------
counts <- fit_fs$find.grps$filter_counts
funnel <- data.frame(
  stage = c("enumerated (at most two cuts, size-eligible)",
            "fittable (ten events in each arm, Cox fit obtained)",
            "clearing the screening threshold, HR > 1.0",
            "passing the consistency screen, p_cons >= 0.90",
            "re-selection family for inference"),
  candidates = c(counts$n_evaluated, counts$n_passed_cox, counts$n_passed_hr,
                 fit_fs$grp.consistency$n_passed, g_fs$n_family))
knitr::kable(funnel, row.names = FALSE, format.args = list(big.mark = ","))


## ----fs-candidates------------------------------------------------------------
cand <- as.data.frame(fit_fs$grp.consistency$out_sg$result)
cand$candidate <- ifelse(cand$M.2 == "", cand$M.1, paste(cand$M.1, cand$M.2, sep = " & "))
cand <- cand[order(-cand$N), ]
knitr::kable(cand[, c("candidate", "N", "hr", "Pcons")], digits = 3,
             row.names = FALSE, col.names = c("candidate", "N", "HR", "p_cons"),
             caption = "The candidates passing the consistency screen, largest first.")
band <- (1 - 0.20) * max(cand$hr)
within_band <- cand$hr >= band
sel <- cand[within_band, ][which.max(cand$N[within_band]), ]
cat(sprintf("Largest HR over the family: %.2f (N %d). Band: HR >= %.2f.\n",
            max(cand$hr), cand$N[which.max(cand$hr)], band))
cat(sprintf("max N(eps) takes the largest candidate within the band: N %d, HR %.2f.\n",
            sel$N, sel$hr))
cat(sprintf("max N alone would take N %d (HR %.2f); max beta the N %d.\n",
            max(cand$N), cand$hr[which.max(cand$N)], cand$N[which.max(cand$hr)]))
stopifnot(sel$N == g_fs$n_selected)   # the rule, re-applied, picks the fit's own selection


## ----table1-helper------------------------------------------------------------
table1_column <- function(fit, data) {
  g  <- fit$mr_inference; cc <- g$complement
  fb <- g$field; fc <- fb$complement; js <- fb$joint_s
  in_H <- fit$df.est$treat.recommend[match(data$id, fit$df.est$id)] == 0
  stopifnot(!anyNA(in_H), sum(in_H) == g$n_selected, sum(!in_H) == cc$n)
  f3 <- function(x) sprintf("%.3f", x)
  ci <- function(e, l, u) sprintf("%s (%s, %s)", f3(e), f3(l), f3(u))
  c(paste(fit$sg.harm, collapse = " & "),
    sprintf("%d (%d)", sum(in_H), sum(data$status[in_H])),
    ci(g$naive$est, g$naive$lower, g$naive$upper),
    f3(g$debiased$est),
    f3(fb$lower_1s),
    sprintf("%d (%d)", sum(!in_H), sum(data$status[!in_H])),
    ci(cc$naive$est, cc$naive$lower, cc$naive$upper),
    f3(cc$debiased$est),
    f3(fc$upper_1s_s),
    sprintf("%s / %s", f3(js$bonf_lower_H), f3(js$bonf_upper_Hc)))
}
table1_rows <- c("Region", "N (events)", "Unadjusted, 95% CI", "Corrected",
                 "Field, lower bound", "Complement N (events)", "Unadjusted, 95% CI",
                 "Corrected", "Field-s, upper bound", "Bonferroni pair")
knitr::kable(data.frame(row = table1_rows, FS = table1_column(fit_fs, df)),
             col.names = c("", "FS"), row.names = FALSE,
             caption = "Forest search: the identified region and its complement.")


## ----fs-supporting------------------------------------------------------------
fb <- g_fs$field; fc <- fb$complement; js <- fb$joint_s
supporting <- data.frame(
  quantity = c("mean of the simulated selection root, region (log scale)",
               "second-order estimate, region",
               "second-order estimate, complement (field-s)",
               "scale ratio rho_c (field-s vs the unstudentized complement bound)",
               "Bonferroni pair: simulated joint coverage",
               "calibrated split gamma (Bonferroni's own alpha/2 is 0.025)",
               "correlation of the aligned roots"),
  value = c(sprintf("%.3f", fb$lambda_mean), sprintf("%.3f", fb$est2),
            sprintf("%.3f", fc$est2_s), sprintf("%.3f", fc$decomp_fields$scale_ratio_c),
            sprintf("%.4f", js$bonf_joint_prob), sprintf("%.3f", js$gamma),
            sprintf("%.3f", js$corr)))
knitr::kable(supporting, row.names = FALSE)


## ----fit-dina-----------------------------------------------------------------
#| message: false
#| results: hide
t0 <- proc.time()[["elapsed"]]
fit_dina <- forestsearch(
  df.analysis = df,
  outcome.name = "time_months", event.name = "status", treat.name = "hormon", id.name = "id",
  confounders.name = confounders,
  is.RCT = TRUE, est.scale = "hr", seedit = 8316951,
  conf.cont_jcuts = list(er = 10, pgr = 10), conf_force = c("er <= 0", "pgr <= 0"),
  n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
  subgroup_method = "dina", dina_select_statistic = "effect",     # DINA's own
  hr.threshold = 1.00, stop_threshold = NULL,
  sg_focus = "effMaxSG", selection_rule = "neighborhood", effect_neighborhood = 0.20,
  mr_inference = TRUE,
  mr_inference_args = list(ci_method = "field", draws = 5000L, include_complement = TRUE,
                           field_scale_complement = "selected", field_decompose = TRUE),
  quiet = TRUE, parallel_args = list(plan = "sequential"))
fit_sec_dina <- proc.time()[["elapsed"]] - t0


## ----fit-grf------------------------------------------------------------------
#| message: false
#| results: hide
t0 <- proc.time()[["elapsed"]]
fit_grf <- forestsearch(
  df.analysis = df,
  outcome.name = "time_months", event.name = "status", treat.name = "hormon", id.name = "id",
  confounders.name = confounders,
  is.RCT = TRUE, est.scale = "hr", seedit = 8316951,
  conf.cont_jcuts = list(er = 10, pgr = 10), conf_force = c("er <= 0", "pgr <= 0"),
  n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
  subgroup_method = "grf", grf_selection = "frontier",             # GRF's own
  grf_select_statistic = "effect", grf_depth = 2L, dmin.grf = 0.0, # GRF's own
  hr.threshold = 1.00, stop_threshold = NULL,
  sg_focus = "effMaxSG", selection_rule = "neighborhood", effect_neighborhood = 0.20,
  mr_inference = TRUE,
  mr_inference_args = list(ci_method = "field", draws = 5000L, include_complement = TRUE,
                           field_scale_complement = "selected", field_decompose = TRUE),
  quiet = TRUE, parallel_args = list(plan = "sequential"))
fit_sec_grf <- proc.time()[["elapsed"]] - t0


## ----dina-grf-selected--------------------------------------------------------
g_dina <- fit_dina$mr_inference
g_grf  <- fit_grf$mr_inference
cat(sprintf("DINA: %s  (fit %.1f s; re-selection family %d)\n",
            paste(fit_dina$sg.harm, collapse = " & "), fit_sec_dina,
            g_dina$n_family))
cat(sprintf("GRF:  %s  (fit %.1f s; re-selection family %d)\n",
            paste(fit_grf$sg.harm, collapse = " & "), fit_sec_grf, g_grf$n_family))


## ----agreement----------------------------------------------------------------
member <- function(fit) {
  fit$df.est$treat.recommend[match(df$id, fit$df.est$id)] == 0
}
m <- data.frame(FS = member(fit_fs), DINA = member(fit_dina), GRF = member(fit_grf))
partition <- aggregate(cbind(N = 1, events = df$status), by = m, FUN = sum)
knitr::kable(partition[order(-partition$N), ], row.names = FALSE,
             caption = "The three-way partition of the 686 patients by the three selections.")


## ----table1-------------------------------------------------------------------
table1 <- data.frame(row = table1_rows,
                     FS = table1_column(fit_fs, df),
                     DINA = table1_column(fit_dina, df),
                     GRF = table1_column(fit_grf, df))
caption1 <- sprintf(paste(
  "Selection-adjusted analysis of the three identified regions, hazard-ratio",
  "scale. Unadjusted: two-sided 95%% Wald interval; field and field-s: one-sided",
  "at nominal 0.95; Bonferroni pair: region / complement, joint 0.95.",
  "Whole trial: %d patients, %d events, hazard ratio %.3f (%.3f, %.3f)."),
  nrow(df), sum(df$status), whole_hr[["est"]], whole_hr[[2]], whole_hr[[3]])
knitr::kable(table1, col.names = c("", "FS", "DINA", "GRF"), row.names = FALSE,
             align = "lrrr", caption = caption1)


## ----reading-values-----------------------------------------------------------
#| echo: false
r3 <- function(x) sprintf("%.3f", x)
ev_region <- vapply(m, function(z) sum(df$status[z]), numeric(1))


## ----declaration--------------------------------------------------------------
c0_grid <- c(0.70, 0.75, 0.80, 0.85)
cal    <- fs_declaration_calibration(fit_fs, alpha = 0.10, c0 = c0_grid)
cal_c2 <- fs_declaration_calibration(fit_fs, alpha = 0.10)   # c0 = NULL: the claim threshold
tb <- cal$c0$table
fld <- g_fs$declaration_field
T_H   <- cal$T_hat[[fld$family_id[fld$meta$selected_index]]]
T_max <- max(cal$T_hat)
kappa <- c(tb$kappa_hat, cal_c2$kappa_hat)
declaration <- data.frame(
  c0 = c(sprintf("HR %.2f%s", tb$c0, ifelse(tb$c0 == 0.75, " (primary)", "")),
         "HR 1.00 (claim threshold)"),
  kappa   = sprintf("%.3f", kappa),
  implied = sprintf("%.3f", 2 * pnorm(kappa) - 1),
  fw      = sprintf("%.3f", c(tb$fw_size, cal_c2$fw_size)),
  admits  = c(lengths(cal$c0$admitted_calibrated),
              length(cal_c2$admitted_calibrated)))
knitr::kable(declaration, row.names = FALSE, align = "lrrrr",
             col.names = c("protected level c0", "kappa-hat", "implied level",
                           "FW size", "admits"),
             caption = paste(
               "Declaration calibration at alpha = 0.10. kappa-hat: the calibrated",
               "cutoff; implied level: 2*Phi(kappa-hat) - 1; FW size: the family-wise",
               "size of the p* = 0.90 screen at c0; admits: candidates the calibrated",
               "rule admits."))
cat(sprintf("T of the declared region: %.3f; largest T in the family: %.3f\n", T_H, T_max))
cat(sprintf("The p* = 0.90 screen's cutoff, as the package applies it: %.3f\n", cal$z_pstar))


## ----session------------------------------------------------------------------
sessionInfo()

