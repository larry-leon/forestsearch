# ---- transplant (verbatim chunk bodies; fences removed) ----
# source: fs-glms-interpretable/quarto/actg175/analysis_actg175_binary_mr.qmd
# fs-glms-interpretable HEAD 19c6ac41c1198dbc513f5c6978d0d6c5df1fc157 (file last changed d97ca1f / 935ae87, 2026-09-23)
#   chunk 'selection': fence lines 335-361
#   chunk 'extract': fence lines 435-498
#   chunk 'payload': fence lines 796-866
# ----------------------------------------------------------

## ==== chunk 'selection' ====
label_fs <- paste(fit$sg.harm, collapse = " & ")
est <- fit$df.est
in_H <- est$treat.recommend[match(df$id, est$id)] == 0
stopifnot(!anyNA(in_H), identical(in_H, subgroup_member(label_fs, df)))
n_H <- sum(in_H); n_Hc <- sum(!in_H)
or_H <- or_est(df[in_H, ])$or
fcn <- fit$find.grps$filter_counts
funnel <- data.frame(
  stage = c("candidates enumerated", "fittable", "clearing the OR screen (≥ 1.00)",
            "evaluated for consistency", "passing consistency (p* = 0.90)", "MR candidate family"),
  n = c(fcn$n_evaluated, fcn$n_passed_cox, fcn$n_passed_hr,
        fit$grp.consistency$n_candidates_evaluated, fit$grp.consistency$n_passed, g$n_family))
knitr::kable(funnel, col.names = c("", "count"), row.names = FALSE, caption = "The candidate funnel.")

# the consistency-passing candidates. On the GLM path this table's `hr` column is the log odds ratio.
cand <- as.data.frame(fit$grp.consistency$out_sg$result)
mcols <- grep("^M\\.", names(cand), value = TRUE)
cand_def <- apply(cand[, mcols, drop = FALSE], 1, function(z) paste(z[nzchar(z) & !is.na(z)], collapse = " & "))
cand_tab <- data.frame(
  subgroup = vapply(seq_along(cand_def), function(i) check_description(cand_def[i], cand$N[i], df), ""),
  N = cand$N, `no improvement` = cand$E, OR = sprintf("%.3f", exp(cand$hr)), Pcons = sprintf("%.2f", cand$Pcons),
  check.names = FALSE)
knitr::kable(cand_tab, row.names = FALSE,
             caption = "Candidates passing the consistency screen, ordered as the package returns them; the first row is the max N(ε) selection. The package's own table carries the effect as a log odds ratio on this path; it is exponentiated here.")
stopifnot(abs(exp(cand$hr[1]) - or_H) < 1e-8, cand$N[1] == n_H, identical(unname(cand_def[1]), label_fs))

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
  rho_c             = fc$decomp_fields$scale_ratio_c,
  n_out_H           = fb$n_out_used,
  n_in_H            = fb$n_in_used_mean,
  n_out_Hc          = fc$n_out_used,
  field_sec_H       = fb$timing_seconds,
  field_sec_Hc      = fc$timing_seconds,
  alpha             = js$alpha,
  # the field's second-order point estimates: est2 = exp(two-term log-OR - lambda_mean), lambda_mean the mean of
  # the simulated root on the log scale (fs_mr_inference(): field$est2, field$lambda_mean; the complement block
  # carries the same pair for the unstudentized root and, with the _s suffix, for field-s, the reported product)
  field_est2_H            = fb$est2,
  field_lambda_mean_H     = fb$lambda_mean,
  field_s_est2_Hc         = fc$est2_s,
  field_s_lambda_mean_Hc  = fc$lambda_mean_s,
  field_est2_Hc           = fc$est2,
  field_lambda_mean_Hc    = fc$lambda_mean)
stopifnot(is.finite(bounds$field_est2_H), is.finite(bounds$field_lambda_mean_H),
          is.finite(bounds$field_s_est2_Hc), is.finite(bounds$field_s_lambda_mean_Hc),
          abs(bounds$field_est2_H - exp(log(g$debiased$est) - bounds$field_lambda_mean_H)) < 1e-10,
          abs(bounds$field_s_est2_Hc - exp(log(cc$debiased$est) - bounds$field_s_lambda_mean_Hc)) < 1e-10)

ev_H  <- sum(df$y_binary_neg[in_H]); ev_Hc <- sum(df$y_binary_neg[!in_H])
stopifnot(ev_H + ev_Hc == sum(df$y_binary_neg))   # the two regions' adverse outcomes make up the analysis set's
h3 <- function(x) sprintf("%.3f", x)
ci <- function(l, u) sprintf("%s–%s", h3(l), h3(u))
whole <- reference_or[1, ]
summary_table <- data.frame(
  row = c("**Reference effects**", "Whole trial",
          "**Identified subgroup**", "Unadjusted", "Two-term", "Field, lower bound",
          "**Complement**", "Unadjusted", "Two-term", "Field-s, upper bound"),
  N = c("", whole$n, "", n_H, n_H, n_H, "", n_Hc, n_Hc, n_Hc),
  Events = c("", whole$ev_ctrl + whole$ev_trt, "", ev_H, "—", "—", "", ev_Hc, "—", "—"),
  OR = c("", h3(whole$or), "",
         h3(g$naive$est), h3(g$debiased$est), h3(bounds$field_est2_H), "",
         h3(cc$naive$est), h3(cc$debiased$est), h3(bounds$field_s_est2_Hc)),
  CI = c("", ci(whole$lo, whole$hi), "",
         ci(g$naive$lower, g$naive$upper), paste(ci(g$debiased$lower, g$debiased$upper), "†"),
         paste(h3(bounds$field_lower_H), "‡"), "",
         ci(cc$naive$lower, cc$naive$upper), paste(ci(cc$debiased$lower, cc$debiased$upper), "†"),
         paste(h3(bounds$field_s_upper_Hc), "‡")),
  stringsAsFactors = FALSE)
knitr::kable(summary_table,
             col.names = c("", "N", "No CD4 improvement", "OR", "95% CI"), align = "lrrrl", row.names = FALSE,
             caption = sprintf("Odds ratios of NO CD4 improvement, ZDV + ddI against ddI; OR above 1 is harm. Whole trial: %d on ZDV + ddI, %d on ddI. Identified subgroup: %s.",
                               whole$n_trt, whole$n_ctrl, check_description(label_fs, n_H, df)))

## ==== chunk 'payload' ====
harm_rates <- local({
  h <- df[in_H, ]
  et <- h$y_binary_neg[h$treat == 1]; ec <- h$y_binary_neg[h$treat == 0]
  list(treated_n = length(et), control_n = length(ec),
       treated_events = sum(et == 1), control_events = sum(ec == 1),
       treated_pct = 100 * mean(et == 1), control_pct = 100 * mean(ec == 1))
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
    outcome = "y_binary_neg", event = "y_binary_neg", covariates = confounders,
    endpoint = "binary: no CD4 improvement at week 20 (adverse scale)",
    settings = list(sg_focus = "effMaxSG", effect_neighborhood = 0.20, selection_rule = "neighborhood",
                    outcome_type = "binary", effect_measure = "OR", adverse_outcome = TRUE,
                    effect.threshold = 1.00, consistency.threshold = 1.00, pconsistency.threshold = 0.90,
                    consistency_method = "resample", max_subgroups_search = Inf,
                    conf.cont_jcuts = list(cd40 = 10, wtkg = 10),
                    seedit = 8316951, mr_settings = g$settings,
                    field_R_out = fb$R_out, field_R_in = fb$R_in),
    thresholds_received = received,
    mr_harm_confirmed = fit$mr_harm_confirmed,
    filter_counts = funnel,
    fit_sec = fit_sec,
    timings = list(fs = list(label = label_fs, n_family = g$n_family, selection_bias = g$selection_bias,
                             boot_sec = NA_real_, mr_sec = g$timing_seconds, speedup = NA_real_))),
  extras = list(concordance = NULL, harm_rates = harm_rates,
                reference_or = reference_or, agreement = agreement,
                counts_cascade = cascade, summary_table = summary_table,
                candidates = cand_tab),
  est_scale = "or",
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

payload_file <- "analysis_actg175_binary_mr_payload.rds"
saveRDS(payload, payload_file)
cat("Wrote", payload_file, "\n")

