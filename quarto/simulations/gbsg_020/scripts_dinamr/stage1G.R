# STAGE 1 battery for the grfmr completion (TASK_grfmr_completion_2026-09-12):
# the 5-replicate smoke at (12.4%, HR 1.00, n 500) on the GRF path.
# stage1_checks.R's battery, repointed at the GRF smoke, plus what the kickoff
# names for this smoke: admitted_n present and populated on the null path; the
# selection rate and the admitted_n distribution; realized prevalence; and how
# many non-detections the smoke produced and whether admitted_n is NA or 0 on
# them.  A non-detection is NOT a failure and is not chased here.
#
# Uses gate2.R's CORRECTED identity (field-s inverted around the same bdc), not
# the bonf-vs-raw pair.  Exit status 1 on any failure (stop-on-failure).
QMD_DIR <- Sys.getenv("DINAMR_QMD_DIR", unset = "..")
f <- file.path(QMD_DIR, "results/grf_effMaxSG_fb_mr_field_m1_h100_knoise0_n500_nb20_grfmrsmk_res_1_5.rds")
b <- readRDS(f); r <- b$results; m <- b$meta
fail <- 0L; pass <- 0L
P <- function(lab, ok, extra = "") { if (isTRUE(ok)) pass <<- pass + 1L else fail <<- fail + 1L
  cat(sprintf("  %-52s %-8s %s\n", lab, if (isTRUE(ok)) "PASS" else "**FAIL**", extra)) }
cat("=== STAGE 1 SMOKE: 12.4% (Z1Q unset), HR 1.00, n 500, 5 replicates, GRF ===\n")
cat(sprintf("stem: %s\n", basename(f)))
cat(sprintf("meta: method=%s focus=%s nbhd=%.2f z1q=%.2f super-prevalence=%.5f hr=%.2f n=%d workers=%s ver=%s host=%s campaign=%s\n",
    m$subgroup_method, m$sg_focus, m$effect_neighborhood, m$harm_z1_quantile, m$harm_prevalence_super,
    m$target_hr_harm, m$n_sample, format(m$n_workers), m$forestsearch_version, m$hostname, m$campaign_tag))
P("meta: subgroup_method grf", identical(m$subgroup_method, "grf"))
P("meta: target_hr_harm 1.00", isTRUE(all.equal(m$target_hr_harm, 1.00)))
P("meta: sg_focus effMaxSG / eps 0.20",
  identical(m$sg_focus, "effMaxSG") && isTRUE(all.equal(m$effect_neighborhood, 0.20)))
P("meta: field knobs + ij two_term",
  isTRUE(m$field_complement) && isTRUE(m$field_decompose) && isTRUE(m$field_recovery) &&
  identical(m$field_scale_complement, "selected") && identical(m$ij_residual, "two_term"))
P("5 rows", nrow(r) == 5L, sprintf("(%d)", nrow(r)))
ce <- sum(grepl("CONFIG", as.character(r$status), ignore.case = TRUE))
P("no CONFIG-ERROR", ce == 0L, sprintf("(%d)", ce))

cat("\n--- per replicate ---\n")
print(data.frame(sim_id = r$sim_id, status = r$status, detected = r$detected,
                 n_family = r$n_family, admitted_n = r$admitted_n, n_sel = r$n_sel,
                 n_true = r$n_true, secs = round(r$fit_mr_secs, 2)), row.names = FALSE)

D <- r[r$detected %in% 1L, , drop = FALSE]; ND <- r[!(r$detected %in% 1L), , drop = FALSE]
sr <- mean(r$detected %in% 1L)
cat(sprintf("\n>> SELECTION RATE : %.4f (%d / %d)  -- a selection rate, nothing more\n", sr, nrow(D), nrow(r)))
A <- r$admitted_n[is.finite(r$admitted_n)]
if (length(A)) cat(sprintf(">> ADMITTED_N     : min %g q25 %g MED %g q75 %g p90 %g max %g  (finite on %d of %d rows)\n",
    min(A), stats::quantile(A,.25,names=FALSE), stats::median(A), stats::quantile(A,.75,names=FALSE),
    stats::quantile(A,.90,names=FALSE), max(A), length(A), nrow(r)))
K <- r$n_family[is.finite(r$n_family)]
if (length(K)) cat(sprintf(">> N_FAMILY (ENUMERATED POOL, not the qualified set): min %g med %g max %g\n",
    min(K), stats::median(K), max(K)))
cat(sprintf(">> NON-DETECTIONS : %d%s\n", nrow(ND),
    if (nrow(ND)) sprintf("  -- admitted_n on them: NA %d | == 0 %d | > 0 %d ; status %s ; err_msg %s",
        sum(is.na(ND$admitted_n)), sum(ND$admitted_n %in% 0L),
        sum(!is.na(ND$admitted_n) & ND$admitted_n > 0L),
        paste(unique(as.character(ND$status)), collapse = "/"),
        paste(unique(as.character(ND$err_msg)), collapse = " | ")) else ""))
cat(sprintf(">> REALIZED PREVALENCE: super-population %.5f | trial mean %.5f (per replicate %s)\n",
    m$harm_prevalence_super, mean(r$n_true)/m$n_sample,
    paste(sprintf("%.4f", r$n_true/m$n_sample), collapse = " ")))
P("realized prevalence ~ 12.4% (super-population in [0.12, 0.13])",
  m$harm_prevalence_super >= 0.12 && m$harm_prevalence_super <= 0.13)

cat("\n--- presence and population (detected replicates) ---\n")
P("admitted_n column present", "admitted_n" %in% names(r))
P("admitted_n populated on the null path (finite on >= 1 detected row)",
  nrow(D) > 0L && any(is.finite(D$admitted_n)))
recov <- c("fld_recov_sens_H","fld_recov_ppv_H","fld_recov_sens_Hc","fld_recov_npv_Hc",
           "fld_recov_q10","fld_recov_q50","fld_recov_q90","fld_recov_share1","fld_recov_n_used")
phat <- c("p_hat_H","p_hat_sum","p_hat_top1")
rhoc <- c("fld_Hc_scale_ratio","fld_Hc_scale_sel","fld_Hc_scale_win","fld_Hc_scale_cv")
pres <- function(cols) all(cols %in% names(r)) &&
  !any(vapply(cols, function(k) all(is.na(D[[k]])), logical(1)))
P("nine recovery columns present & populated", pres(recov))
P("p-hat block present & populated", pres(phat))
P("rho-c / scale block present & populated", pres(rhoc))

fin <- c("nv_H_est","nv_Hc_est","mr_H_est","mr_Hc_est","mr_H_lo","mr_H_hi","mr_Hc_lo","mr_Hc_hi",
         "fld_H_est2","fld_H_lo1s","fld_H_se","fld_Hc_est2","fld_Hc_up1s","fld_Hc_se",
         "fld_Hc_est2_s","fld_Hc_up1s_s","fld_Hc_se_s",
         "fld_joint_gamma","fld_joint_bonf_loH","fld_joint_bonf_upHc",
         "fld_joint_s_gamma","fld_joint_s_bonf_loH","fld_joint_s_bonf_upHc",
         "betaHhat_H","betaHhat_Hc","sens","spec","ppv","npv", phat, rhoc, recov)
fin <- fin[fin %in% names(r)]
nf <- fin[vapply(fin, function(k) any(!is.finite(D[[k]])), logical(1))]
P("every construction finite on detected replicates", length(nf) == 0L,
  if (length(nf)) paste("non-finite:", paste(nf, collapse = ", ")) else sprintf("(%d quantities)", length(fin)))
P("classification set defined (sens/spec/PPV/NPV finite)",
  all(is.finite(c(D$sens, D$spec, D$ppv, D$npv))))

cat("\n--- invariants ---\n")
P("fld_H_lo1s <= fld_H_est2",           all(D$fld_H_lo1s <= D$fld_H_est2))
P("fld_Hc_est2 <= fld_Hc_up1s",         all(D$fld_Hc_est2 <= D$fld_Hc_up1s))
P("fld_Hc_est2_s <= fld_Hc_up1s_s",     all(D$fld_Hc_est2_s <= D$fld_Hc_up1s_s))
P("joint : bonf_loH <= fld_H_est2",     all(D$fld_joint_bonf_loH <= D$fld_H_est2))
P("joint : fld_Hc_est2 <= bonf_upHc",   all(D$fld_Hc_est2 <= D$fld_joint_bonf_upHc))
P("jointS: bonf_loH <= fld_H_est2",     all(D$fld_joint_s_bonf_loH <= D$fld_H_est2))
P("jointS: fld_Hc_est2_s <= bonf_upHc", all(D$fld_Hc_est2_s <= D$fld_joint_s_bonf_upHc))
P("IJ    : mr_H_lo <= est <= mr_H_hi",  all(D$mr_H_lo <= D$mr_H_est & D$mr_H_est <= D$mr_H_hi))
P("IJ    : mr_Hc_lo <= est <= mr_Hc_hi",all(D$mr_Hc_lo <= D$mr_Hc_est & D$mr_Hc_est <= D$mr_Hc_hi))
P("rho-c : scale_ratio > 0",            all(D$fld_Hc_scale_ratio > 0))
P("recov : sens_H, ppv_H in [0,1]",     all(D$fld_recov_sens_H >= 0 & D$fld_recov_sens_H <= 1 &
                                            D$fld_recov_ppv_H >= 0 & D$fld_recov_ppv_H <= 1))
P("p-hat : 0 <= p_hat_H <= 1, <= sum, <= top1",
  all(D$p_hat_H >= 0 & D$p_hat_H <= 1 & D$p_hat_H <= D$p_hat_sum + 1e-12 & D$p_hat_H <= D$p_hat_top1 + 1e-12))
P("admitted_n >= 1 on every detected replicate", all(D$admitted_n >= 1L, na.rm = TRUE))
g1 <- D$fld_joint_gamma; g2 <- D$fld_joint_s_gamma
P("gamma (joint)   in [0.025, 0.05]", all(g1 >= 0.025 & g1 <= 0.05), sprintf("[%.5f, %.5f]", min(g1), max(g1)))
P("gamma (joint-s) in [0.025, 0.05]", all(g2 >= 0.025 & g2 <= 0.05), sprintf("[%.5f, %.5f]", min(g2), max(g2)))
i1 <- max(abs((log(D$fld_Hc_est2_s) + D$fld_Hc_lam_mean_s) - (log(D$fld_Hc_est2) + D$fld_Hc_lam_mean)))
P("identity: field-s inverted around the same bdc", i1 <= 1e-12, sprintf("max |diff| = %.3g", i1))

cat("\n--- structural (reported as such, never a failure) ---\n")
for (k in c("n_cons_qual","band_n","p_star"))
  cat(sprintf("  STRUCTURAL %-12s %s\n", k,
    if (!k %in% names(r)) "not a recorder column (admission-set term; NULL on GRF)"
    else if (all(is.na(r[[k]]))) "present, all-NA -- structural on GRF (no consistency screen)"
    else "POPULATED (unexpected on GRF)"))
cat(sprintf("\n  oracle target log(1.00) = 0; planted marg_H %.4f, marg_Hc %.4f\n",
    b$truth$marg_H, b$truth$marg_Hc))
cat(sprintf("  on detected: sens %s | spec %s | ppv %s | npv %s\n",
    paste(sprintf("%.3f", D$sens), collapse=" "), paste(sprintf("%.3f", D$spec), collapse=" "),
    paste(sprintf("%.3f", D$ppv), collapse=" "), paste(sprintf("%.3f", D$npv), collapse=" ")))
cat(sprintf("  field lower bound (HR): %s ; theta(Hhat): %s\n",
    paste(sprintf("%.3f", D$fld_H_lo1s), collapse=" "), paste(sprintf("%.3f", D$betaHhat_H), collapse=" ")))
cat(sprintf("\nSTAGE 1: %d passes, %d failures -- %s\n", pass, fail, if (fail) "FAIL: STOP" else "PASS"))
if (fail) quit(status = 1L)
