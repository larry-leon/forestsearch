chk <- function(path, label) {
  b <- readRDS(path); r <- b$results; m <- b$meta
  cat("\n================ ", label, " ================\n", sep="")
  cat(sprintf("stem: %s\n", basename(path)))
  cat(sprintf("meta: method=%s focus=%s nbhd=%.2f z1q=%.2f realized super-prevalence=%.5f hr=%.2f n=%d workers=%d ver=%s host=%s\n",
      m$subgroup_method, m$sg_focus, m$effect_neighborhood, m$harm_z1_quantile,
      m$harm_prevalence_super, m$target_hr_harm, m$n_sample, m$n_workers, m$forestsearch_version, m$hostname))
  cat(sprintf("meta: field_complement=%s field_decompose=%s field_scale_complement=%s field_recovery=%s ij_residual=%s fb_mode=%s er_jcuts=%s\n",
      m$field_complement, m$field_decompose, m$field_scale_complement, m$field_recovery, m$ij_residual, m$fb_mode, m$er_jcuts))
  cat(sprintf("rows=%d cols=%d  detected=%d/%d\n", nrow(r), ncol(r), sum(r$detected %in% 1L), nrow(r)))
  cat(sprintf("realized trial prevalence n_true/n: %s  (mean %.4f)\n",
      paste(sprintf("%.4f", r$n_true/m$n_sample), collapse=" "), mean(r$n_true)/m$n_sample))
  cat(sprintf("PROPOSED-FAMILY SIZE per replicate (n_family): %s\n", paste(r$n_family, collapse=" ")))
  cat(sprintf("n_sel: %s   fit_mr_secs: %s\n", paste(r$n_sel, collapse=" "), paste(sprintf("%.2f", r$fit_mr_secs), collapse=" ")))

  D <- r[r$detected %in% 1L, , drop=FALSE]

  # --- PRESENCE: the nine recovery columns, p-hat, rho-c ---
  recov <- c("fld_recov_sens_H","fld_recov_ppv_H","fld_recov_sens_Hc","fld_recov_npv_Hc",
             "fld_recov_q10","fld_recov_q50","fld_recov_q90","fld_recov_share1","fld_recov_n_used")
  phat  <- c("p_hat_H","p_hat_sum","p_hat_top1")
  rhoc  <- c("fld_Hc_scale_ratio","fld_Hc_scale_sel","fld_Hc_scale_win","fld_Hc_scale_cv")
  pres <- function(cols, nm) {
    miss <- cols[!cols %in% names(r)]
    allna <- cols[cols %in% names(r)][vapply(cols[cols %in% names(r)], function(k) all(is.na(D[[k]])), logical(1))]
    cat(sprintf("%-26s column PRESENT: %d/%d%s | all-NA on detected: %s\n", nm,
        sum(cols %in% names(r)), length(cols),
        if (length(miss)) paste0("  MISSING: ", paste(miss, collapse=",")) else "",
        if (length(allna)) paste(allna, collapse=",") else "none"))
    length(miss) == 0L && length(allna) == 0L
  }
  ok_recov <- pres(recov, "nine recovery columns")
  ok_phat  <- pres(phat,  "p-hat block")
  ok_rhoc  <- pres(rhoc,  "rho-c / scale block")

  # --- FINITENESS on detected replicates ---
  fin_cols <- c("nv_H_est","nv_Hc_est","mr_H_est","mr_Hc_est","mr_H_lo","mr_H_hi","mr_Hc_lo","mr_Hc_hi",
                "fld_H_est2","fld_H_lo1s","fld_H_se","fld_Hc_est2","fld_Hc_up1s","fld_Hc_se",
                "fld_Hc_est2_s","fld_Hc_up1s_s","fld_Hc_se_s",
                "fld_joint_gamma","fld_joint_bonf_loH","fld_joint_bonf_upHc",
                "fld_joint_s_gamma","fld_joint_s_bonf_loH","fld_joint_s_bonf_upHc",
                "betaHhat_H","betaHhat_Hc", phat, rhoc, recov)
  fin_cols <- fin_cols[fin_cols %in% names(r)]
  nonfin <- fin_cols[vapply(fin_cols, function(k) any(!is.finite(D[[k]])), logical(1))]
  cat(sprintf("FINITENESS over %d constructions on %d detected replicates: %s\n", length(fin_cols), nrow(D),
      if (length(nonfin)) paste("NON-FINITE in:", paste(nonfin, collapse=", ")) else "ALL FINITE"))

  # --- INTERVAL INVARIANTS ---
  inv <- list(
    "harm  : fld_H_lo1s <= fld_H_est2"          = all(D$fld_H_lo1s   <= D$fld_H_est2),
    "compl : fld_Hc_est2 <= fld_Hc_up1s"        = all(D$fld_Hc_est2  <= D$fld_Hc_up1s),
    "_s    : fld_Hc_est2_s <= fld_Hc_up1s_s"    = all(D$fld_Hc_est2_s<= D$fld_Hc_up1s_s),
    "joint : bonf_loH <= fld_H_est2"            = all(D$fld_joint_bonf_loH   <= D$fld_H_est2),
    "joint : fld_Hc_est2 <= bonf_upHc"          = all(D$fld_Hc_est2          <= D$fld_joint_bonf_upHc),
    "jointS: bonf_loH <= fld_H_est2"            = all(D$fld_joint_s_bonf_loH <= D$fld_H_est2),
    "jointS: fld_Hc_est2_s <= bonf_upHc_s"      = all(D$fld_Hc_est2_s        <= D$fld_joint_s_bonf_upHc),
    "IJ    : mr_H_lo <= mr_H_est <= mr_H_hi"    = all(D$mr_H_lo <= D$mr_H_est & D$mr_H_est <= D$mr_H_hi),
    "IJ    : mr_Hc_lo <= mr_Hc_est <= mr_Hc_hi" = all(D$mr_Hc_lo <= D$mr_Hc_est & D$mr_Hc_est <= D$mr_Hc_hi),
    "rho-c : fld_Hc_scale_ratio > 0"            = all(D$fld_Hc_scale_ratio > 0),
    "recov : 0 <= sens_H <= 1"                  = all(D$fld_recov_sens_H >= 0 & D$fld_recov_sens_H <= 1),
    "recov : 0 <= ppv_H  <= 1"                  = all(D$fld_recov_ppv_H  >= 0 & D$fld_recov_ppv_H  <= 1),
    "p-hat : 0 <= p_hat_H <= 1"                 = all(D$p_hat_H >= 0 & D$p_hat_H <= 1))
  for (nm in names(inv)) cat(sprintf("  INV %-42s %s\n", nm, if (isTRUE(inv[[nm]])) "OK" else "**VIOLATED**"))

  # --- GAMMA in range ---
  g1 <- D$fld_joint_gamma; g2 <- D$fld_joint_s_gamma
  cat(sprintf("GAMMA joint   : %s   in [0.025,0.05]: %s\n", paste(sprintf("%.5f",g1),collapse=" "),
      all(g1 >= 0.025 & g1 <= 0.05)))
  cat(sprintf("GAMMA joint-s : %s   in [0.025,0.05]: %s\n", paste(sprintf("%.5f",g2),collapse=" "),
      all(g2 >= 0.025 & g2 <= 0.05)))

  # --- BOUND <-> QUANTILE identity (<= 1e-12) ---
  d1 <- max(abs(D$fld_joint_bonf_loH - D$fld_joint_loH), na.rm=TRUE)
  d2 <- max(abs(D$fld_joint_bonf_upHc - D$fld_joint_upHc), na.rm=TRUE)
  cat(sprintf("BOUND<->QUANTILE max |bonf - raw| : loH %.3g  upHc %.3g   (<=1e-12: %s)\n",
      d1, d2, (d1 <= 1e-12) && (d2 <= 1e-12)))

  # --- structurally-NA columns ---
  sn <- c("n_cons_qual","band_n","p_star")
  for (k in sn) cat(sprintf("STRUCTURAL %-12s %s\n", k,
      if (!k %in% names(r)) "not a recorder column (admission-set term)"
      else if (all(is.na(r[[k]]))) "present, all-NA (structural on DINA: no consistency screen)"
      else paste("POPULATED:", paste(r[[k]], collapse=" "))))

  # --- recovery values ---
  cat("RECOVERY means on detected: ")
  cat(paste(sprintf("%s=%.4f", sub("fld_recov_","",recov), colMeans(D[,recov,drop=FALSE], na.rm=TRUE)), collapse="  "), "\n")
  cat(sprintf("p_hat_H: %s\n", paste(sprintf("%.4f", D$p_hat_H), collapse=" ")))
  cat(sprintf("rho-c (fld_Hc_scale_ratio): %s\n", paste(sprintf("%.4f", D$fld_Hc_scale_ratio), collapse=" ")))
  invisible(list(ok_recov=ok_recov, ok_phat=ok_phat, ok_rhoc=ok_rhoc, r=r, m=m))
}
a <- chk("results/dina_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_nb20_dinamrsmk_res_1_5.rds",
         "SMOKE 1 -- 12.4% (Z1Q unset), HR 1.50, n 500")
b <- chk("results/dina_effMaxSG_fb_mr_field_m1_h150_knoise0_n1500_z1q60_nb20_dinamrsmk_res_1_5.rds",
         "SMOKE 2 -- 31% (Z1Q=0.60), HR 1.50, n 1500")
