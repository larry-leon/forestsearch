# ---- transplant (verbatim chunk bodies; fences removed) ----
# source: fs-glms-interpretable/quarto/actg175/analysis_actg175_continuous_mr.qmd
# fs-glms-interpretable HEAD 19c6ac41c1198dbc513f5c6978d0d6c5df1fc157 (file last changed d97ca1f / 935ae87, 2026-09-23)
#   chunk 'selection': fence lines 357-400
#   chunk 'extract': fence lines 421-480
#   chunk 'payload': fence lines 773-845
# ----------------------------------------------------------

## ==== chunk 'selection' ====
label_fs <- paste(fit$sg.harm, collapse = " & ")
est <- fit$df.est
in_H <- est$treat.recommend[match(df$id, est$id)] == 0
stopifnot(!anyNA(in_H), identical(in_H, subgroup_member(label_fs, df)))
n_H <- sum(in_H); n_Hc <- sum(!in_H)
md_H <- md_est(df[in_H, ])$md
cand <- as.data.frame(fit$grp.consistency$out_sg$result)
pcons_H <- cand$Pcons[1]

anchor <- data.frame(
  check = c("subgroup", "n", "oriented MD", "Pcons"),
  recorded = c(describe_subgroup("!{cd40 <= 507} & {gender}"), "79", "71.10", "0.91"),
  this_fit = c(describe_subgroup(label_fs), n_H, sprintf("%.2f", md_H), sprintf("%.2f", pcons_H)))
knitr::kable(anchor, row.names = FALSE,
             caption = "Gate 3.1: the effMaxSG anchor of REPORT_actg175_applied_effmaxsg_stage0_2026-09-16.md §6.")
if (!identical(label_fs, "!{cd40 <= 507} & {gender}") || n_H != 79L ||
    abs(md_H - 71.10) > 0.005 || abs(pcons_H - 0.91) > 0.005)
  stop("the fit does not reproduce the recorded effMaxSG anchor")

fcn <- fit$find.grps$filter_counts
funnel <- data.frame(
  stage = c("candidates enumerated", "passing the variance and prevalence checks", "passing redundancy",
            "passing the per-arm minimum (10 / 10)", "clearing the screen (MD ≥ 10, oriented)",
            "evaluated for consistency", "passing consistency (c2 = 5, p* = 0.90)",
            "on the Pareto frontier", "MR candidate family"),
  n = c(fcn$n_evaluated, fcn$n_passed_prevalence, fcn$n_passed_redundancy, fcn$n_passed_sample_size,
        fcn$n_passed_hr, fit$grp.consistency$n_candidates_evaluated, fit$grp.consistency$n_passed,
        nrow(fit$grp.consistency$out_sg$pareto_frontier), g$n_family))
knitr::kable(funnel, col.names = c("", "count"), row.names = FALSE, caption = "The candidate funnel.")

# the consistency-passing candidates. On this path the package's `hr` column is the effect on the
# working (identity) scale, with no transform to undo, and `E` repeats `N`: a continuous outcome
# has no events.
mcols <- grep("^M\\.", names(cand), value = TRUE)
cand_def <- apply(cand[, mcols, drop = FALSE], 1, function(z) paste(z[nzchar(z) & !is.na(z)], collapse = " & "))
cand_tab <- data.frame(
  subgroup = vapply(seq_along(cand_def), function(i) check_description(cand_def[i], cand$N[i], df), ""),
  N = cand$N, `oriented MD` = sprintf("%.1f", cand$hr), Pcons = sprintf("%.2f", cand$Pcons),
  check.names = FALSE)
knitr::kable(cand_tab, row.names = FALSE,
             caption = sprintf("All %d candidates passing consistency, ordered as the package returns them; the first row is the max N(ε) selection. Mean differences in CD4 cells per mm³, oriented.", nrow(cand)))
stopifnot(abs(cand$hr[1] - md_H) < 1e-8, cand$N[1] == n_H, identical(unname(cand_def[1]), label_fs))

## ==== chunk 'extract' ====
cc <- g$complement
fb <- g$field
fc <- fb$complement
js <- fb$joint_s

bounds <- list(
  field_lower_H     = fb$lower_1s,
  field_s_upper_Hc  = fc$upper_1s_s,
  bonf_lower_H      = js$bonf_lower_H,
  bonf_upper_Hc     = js$bonf_upper_Hc,
  bonf_joint_prob   = js$bonf_joint_prob,
  calib_gamma       = js$gamma,
  calib_joint_prob  = js$joint_prob,
  calib_lower_H     = js$lower_H,
  calib_upper_Hc    = js$upper_Hc,
  joint_corr        = js$corr,
  n_joint_draws     = js$n_joint_draws,
  # the applied document does not pass field_decompose (package default FALSE), so the
  # complement's decomposition fields, and with them the scale ratio rho-c, are not computed
  rho_c             = if (is.null(fc$decomp_fields)) NA_real_ else fc$decomp_fields$scale_ratio_c,
  n_out_H           = fb$n_out_used,
  n_in_H            = fb$n_in_used_mean,
  n_out_Hc          = fc$n_out_used,
  field_sec_H       = fb$timing_seconds,
  field_sec_Hc      = fc$timing_seconds,
  alpha             = js$alpha,
  # est2 = two-term - lambda_mean on the WORKING scale, here the identity scale, so lambda_mean
  # is itself a number of CD4 cells and the correction is a subtraction, not a ratio
  field_est2_H            = fb$est2,
  field_lambda_mean_H     = fb$lambda_mean,
  field_s_est2_Hc         = fc$est2_s,
  field_s_lambda_mean_Hc  = fc$lambda_mean_s,
  field_est2_Hc           = fc$est2,
  field_lambda_mean_Hc    = fc$lambda_mean)
stopifnot(is.finite(bounds$field_est2_H), is.finite(bounds$field_lambda_mean_H),
          is.finite(bounds$field_s_est2_Hc), is.finite(bounds$field_s_lambda_mean_Hc),
          abs(bounds$field_est2_H - (g$debiased$est - bounds$field_lambda_mean_H)) < 1e-9,
          abs(bounds$field_s_est2_Hc - (cc$debiased$est - bounds$field_s_lambda_mean_Hc)) < 1e-9)

whole <- reference_md[1, ]
ci <- function(l, u) sprintf("%.1f–%.1f", l, u)
summary_table <- data.frame(
  row = c("**Reference effects**", "Whole trial",
          "**Identified subgroup**", "Unadjusted", "Two-term", "Field, lower bound",
          "**Complement**", "Unadjusted", "Two-term", "Field-s, upper bound"),
  N = c("", whole$n, "", n_H, n_H, n_H, "", n_Hc, n_Hc, n_Hc),
  MD = c("", f1(whole$md), "",
         f1(g$naive$est), f1(g$debiased$est), f1(bounds$field_est2_H), "",
         f1(cc$naive$est), f1(cc$debiased$est), f1(bounds$field_s_est2_Hc)),
  CI = c("", ci(whole$lo, whole$hi), "",
         ci(g$naive$lower, g$naive$upper), paste(ci(g$debiased$lower, g$debiased$upper), "†"),
         paste(f1(bounds$field_lower_H), "‡"), "",
         ci(cc$naive$lower, cc$naive$upper), paste(ci(cc$debiased$lower, cc$debiased$upper), "†"),
         paste(f1(bounds$field_s_upper_Hc), "‡")),
  stringsAsFactors = FALSE)
knitr::kable(summary_table, col.names = c("", "N", "MD", "95% CI"), align = "lrrl", row.names = FALSE,
             caption = sprintf("Mean difference in CD4 decline to week 20, **CD4 cells per mm³**: positive means a worse CD4 outcome under ZDV + ddI, negative means benefit from it. Whole trial: %d on ZDV + ddI, %d on ddI. Identified subgroup: %s.",
                               whole$n_trt, whole$n_ctrl, check_description(label_fs, n_H, df)))

## ==== chunk 'payload' ====
arm_means <- local({
  h <- df[in_H, ]
  list(treated_n = sum(h$treat == 1), control_n = sum(h$treat == 0),
       treated_mean_decline = mean(h$y_decline[h$treat == 1]),
       control_mean_decline = mean(h$y_decline[h$treat == 0]))
})
table1 <- data.frame(
  method = "FS", region = c("H", "Hc"),
  n = c(n_H, n_Hc), pct = round(100 * c(n_H, n_Hc) / nrow(df)),
  naive_est = c(g$naive$est, cc$naive$est), naive_lo = c(g$naive$lower, cc$naive$lower),
  naive_hi = c(g$naive$upper, cc$naive$upper),
  fb_est = NA_real_, fb_lo = NA_real_, fb_hi = NA_real_,
  mr_est = c(g$debiased$est, cc$debiased$est), mr_lo = c(g$debiased$lower, cc$debiased$lower),
  mr_hi = c(g$debiased$upper, cc$debiased$upper),
  stringsAsFactors = FALSE)

si <- sessionInfo()
payload <- list(
  table  = table1,
  labels = list(fs = label_fs, dina = NA_character_, grf = NA_character_,
                dina_standalone = list(found = FALSE, label = NA_character_, n = NA_integer_)),
  bounds = bounds,
  meta   = list(
    B = NA_integer_, mr_draws = mr_draws, n_total = nrow(df),
    cores = 1L,
    machine = list(sysname = unname(Sys.info()[["sysname"]]), arch = unname(Sys.info()[["machine"]]),
                   cores_used = 1L, blas = si$BLAS, r_version = R.version.string,
                   timestamp = format(Sys.time(), tz = "UTC", usetz = TRUE)),
    outcome = "y_decline", event = NA_character_, covariates = confounders,
    endpoint = "continuous: CD4 decline to week 20, oriented (positive = worse under ZDV + ddI)",
    specification_source = "forestsearch quarto/applications/actg175/analysis_actg175_continuous_oc.qmd (lines cited in the document); sg_focus and effect_neighborhood are this workstream's",
    settings = list(sg_focus = "effMaxSG", effect_neighborhood = 0.20, selection_rule = "neighborhood",
                    outcome_type = "continuous", effect_measure = "MD", adverse_outcome = TRUE,
                    effect.threshold = 10, consistency.threshold = 5, pconsistency.threshold = 0.90,
                    consistency_method = "resample", n.min = 60, d0.min = 10, d1.min = 10,
                    maxk = 2, fs.splits = 500, use_twostage = TRUE,
                    conf.cont_jcuts = list(age = 10, preanti = 10, wtkg = 10,
                                           karnof = 10, cd40 = 10, cd80 = 10),
                    seedit = 8316951, mr_settings = g$settings,
                    field_R_out = fb$R_out, field_R_in = fb$R_in),
    thresholds_received = received,
    mr_harm_confirmed = fit$mr_harm_confirmed,
    filter_counts = funnel,
    anchor_check = anchor,
    fit_sec = fit_sec,
    timings = list(fs = list(label = label_fs, n_family = g$n_family, selection_bias = g$selection_bias,
                             boot_sec = NA_real_, mr_sec = g$timing_seconds, speedup = NA_real_))),
  extras = list(concordance = NULL, arm_means = arm_means,
                reference_md = reference_md, candidates = cand_tab,
                counts_cascade = cascade, summary_table = summary_table),
  est_scale = "md_decline_oriented",
  built_at = Sys.time(),
  forestsearch_version = as.character(packageVersion("forestsearch")),
  # declaration calibration (§ Declaration calibration): added after every existing element
  declaration_calibration = list(alpha10 = cal, alpha05 = cal_05, boundary = cal_c2, c0_grid = c0_grid,
                                 c0_primary = c0_primary, T_by_candidate = T_by_candidate, max_T = max_T,
                                 T_H = T_H, id_H = id_H, applied_admitted = applied, table = declcal_table,
                                 family_size_prereduction = cal$n_family_prereduction,
                                 # the build that wrote this payload (TASK_declcal_applications_fsglms_2026-09-22_v2)
                                 build = list(
                                   r_version = R.version.string, platform = R.version$platform,
                                   os = utils::sessionInfo()$running,
                                   blas = extSoftVersion()[["BLAS"]], lapack = La_library(),
                                   forestsearch_version = as.character(packageVersion("forestsearch")),
                                   forestsearch_built = packageDescription("forestsearch")$Built,
                                   forestsearch_packaged = packageDescription("forestsearch")$Packaged,
                                   note = "Built under the forestsearch recorded in `build$forestsearch_built`.")))

payload_file <- "analysis_actg175_continuous_mr_payload.rds"
saveRDS(payload, payload_file)
cat("Wrote", payload_file, "\n")

