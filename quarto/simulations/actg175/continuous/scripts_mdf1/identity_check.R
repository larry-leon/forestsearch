# Stage 1b identities: 5-replicate template renders (ij and field) vs the committed bundles.
cells <- list(
  list(tag = "md40 n500",  new = "fs_maxeffCons_mr_field_md40_knoise0_n500_%s_d5000/fs_maxeffCons_mr_field_md40_knoise0_n500_%s_res_1_5.rds",
       old = "fs_maxeffCons_mr_md40_knoise0_n500_s1000_d5000/fs_maxeffCons_mr_md40_knoise0_n500_res_1_1000.rds"),
  list(tag = "md40 n700",  new = "fs_maxeffCons_mr_field_md40_knoise0_n700_%s_d5000/fs_maxeffCons_mr_field_md40_knoise0_n700_%s_res_1_5.rds",
       old = "fs_maxeffCons_mr_md40_knoise0_n700_s1000_d5000/fs_maxeffCons_mr_md40_knoise0_n700_res_1_1000.rds"),
  list(tag = "md120 n500", new = "fs_maxeffCons_mr_field_md120_knoise0_n500_%s_d5000/fs_maxeffCons_mr_field_md120_knoise0_n500_%s_res_1_5.rds",
       old = "fs_maxeffCons_mr_md120_knoise0_n500_s1000_d5000/fs_maxeffCons_mr_md120_knoise0_n500_res_1_1000.rds"),
  list(tag = "null n500",  new = "fs_maxeffCons_mr_field_mdnull_knoise0_n500_%s_d5000/fs_maxeffCons_mr_field_mdnull_knoise0_n500_%s_res_1_5.rds",
       old = "fs_maxeffCons_mr_mdnull_knoise0_n500_s1000_d5000/fs_maxeffCons_mr_mdnull_knoise0_n500_res_1_1000.rds"))
relmax <- function(a, z) {                      # max relative diff over numeric columns; NA==NA ok
  out <- c()
  for (k in names(a)) {
    x <- a[[k]]; y <- z[[k]]
    if (!is.numeric(x)) next
    both_na <- is.na(x) & is.na(y)
    if (any(is.na(x) != is.na(y))) { out[k] <- Inf; next }
    d <- abs(x - y) / pmax(abs(y), 1e-300); d[both_na] <- 0
    out[k] <- max(d)
  }
  out
}
tol <- 1e-8; all_ok <- TRUE
for (cl in cells) {
  old <- readRDS(file.path("mr_md_harm", cl$old)); ro <- old$results[old$results$sim_id %in% 1:5, ]
  ro <- ro[order(ro$sim_id), ]
  pre <- setdiff(names(ro), c("fit_mr_secs", "fb_secs", "mr_msg"))   # twin's 59 columns minus timings/messages
  cat(sprintf("\n=== %s === committed pkg %s | truth effect_Q %s\n", cl$tag, old$meta$pkg_version, format(old$truth$effect_Q, digits = 12)))
  for (ci in c("idij", "idfield")) {
    nb <- readRDS(file.path("mr_md_harm", sprintf(cl$new, ci, ci))); rn <- nb$results[order(nb$results$sim_id), ]
    stopifnot(identical(rn$sim_id, ro$sim_id))
    rm <- relmax(rn[, pre], ro[, pre]); rm <- rm[is.finite(rm) | is.infinite(rm)]
    ok_num <- all(rm <= tol)
    # Rules compared as SETS of terms (package vintages order conjunction terms differently).
    .terms <- function(v) lapply(strsplit(v, " & ", fixed = TRUE), function(t) sort(t))
    rules_same <- identical(.terms(rn$sg_def), .terms(ro$sg_def))
    if (!identical(rn$sg_def, ro$sg_def)) cat(sprintf("  NOTE: sg_def strings differ on sim %s (term order only: %s)\n",
        paste(rn$sim_id[rn$sg_def != ro$sg_def], collapse = ","), rules_same))
    ok_chr <- rules_same && identical(rn$status, ro$status) && identical(rn$n_harm, ro$n_harm) &&
              identical(rn$ij_source, ro$ij_source) && identical(rn$betaHhat_status, ro$betaHhat_status)
    ok_truth <- isTRUE(all.equal(nb$truth, old$truth, tolerance = 1e-10))
    cat(sprintf("  [%s] pre-existing numeric cols: max rel diff %.2e over %d cols -> %s | rule-sets/status/n_harm/ij_source identical: %s | truth equal: %s | timing mean fit_mr %.2f s\n",
                ci, max(rm), length(rm), if (ok_num) "PASS" else "FAIL", ok_chr, ok_truth, mean(rn$fit_mr_secs)))
    if (!ok_num) print(sort(rm[rm > tol], decreasing = TRUE))
    all_ok <- all_ok && ok_num && ok_chr && ok_truth
    if (ci == "idfield") {
      d <- rn[rn$detected %in% 1L, ]
      fin <- function(v) all(is.finite(v))
      ff <- fin(d$fld_H_est2) && fin(d$fld_H_lo1s) && fin(d$fld_H_lo2s) && fin(d$fld_H_hi2s) && fin(d$fld_H_se)
      fc <- fin(d$fld_Hc_est2) && fin(d$fld_Hc_up1s) && fin(d$fld_Hc_lo2s) && fin(d$fld_Hc_hi2s) && fin(d$fld_Hc_se)
      # Bound identities at 1e-12 (absolute, MD scale): inversion about beta-tilde.
      idn <- c(lo1s = max(abs(d$fld_H_lo1s - (d$mr_H_est - d$fld_H_q95))),
               lo2s = max(abs(d$fld_H_lo2s - (d$mr_H_est - d$fld_H_q975))),
               hi2s = max(abs(d$fld_H_hi2s - (d$mr_H_est - d$fld_H_q025))),
               est2 = max(abs(d$fld_H_est2 - (d$mr_H_est - d$fld_H_lam_mean))),
               lo_se = max(abs(d$fld_H_lo_se - (d$fld_H_est2 - qnorm(.975) * d$fld_H_se))),
               c_up1s = max(abs(d$fld_Hc_up1s - (d$mr_Hc_est - d$fld_Hc_q05))),
               c_lo1s = max(abs(d$fld_Hc_lo1s - (d$mr_Hc_est - d$fld_Hc_q95))),
               c_lo2s = max(abs(d$fld_Hc_lo2s - (d$mr_Hc_est - d$fld_Hc_q975))),
               c_hi2s = max(abs(d$fld_Hc_hi2s - (d$mr_Hc_est - d$fld_Hc_q025))),
               c_est2 = max(abs(d$fld_Hc_est2 - (d$mr_Hc_est - d$fld_Hc_lam_mean))))
      gam_ok <- all(d$fld_joint_gamma >= 0.025 - 1e-12 & d$fld_joint_gamma <= 0.05 + 1e-12)
      ph_ok  <- all(is.finite(d$p_hat_H) & d$p_hat_H >= 0 & d$p_hat_H <= 1)
      # IJ columns identical between the ij and field renders (same seed stream).
      ri <- readRDS(file.path("mr_md_harm", sprintf(cl$new, "idij", "idij"))); ri <- ri$results[order(ri$results$sim_id), ]
      ij_same <- max(relmax(rn[, pre], ri[, pre]))
      cat(sprintf("  [field] detected %d/5 | field finite H: %s, Hc: %s | bound identities max abs %.2e (%s) | gamma in [0.025,0.05]: %s (values %s) | p_hat finite in [0,1]: %s (values %s) | IJ cols vs ij render: max rel %.1e | nout min %d | Hc nfit mean %.0f | secs field/complement mean %.2f/%.2f\n",
                  nrow(d), ff, fc, max(idn), if (max(idn) <= 1e-12) "PASS" else "FAIL", gam_ok,
                  paste(format(d$fld_joint_gamma), collapse = ","), ph_ok, paste(format(round(d$p_hat_H, 3)), collapse = ","),
                  ij_same, min(d$fld_H_nout), mean(d$fld_Hc_nfit), mean(d$fld_H_secs), mean(d$fld_Hc_secs)))
      cat(sprintf("          field lo1s: %s | Hc up1s: %s | oriented betaHhat_H: %s\n",
                  paste(sprintf("%.1f", d$fld_H_lo1s), collapse = ","), paste(sprintf("%.1f", d$fld_Hc_up1s), collapse = ","),
                  paste(sprintf("%.1f", -d$betaHhat_H), collapse = ",")))
      all_ok <- all_ok && ff && fc && max(idn) <= 1e-12 && gam_ok && ph_ok && ij_same <= tol
    }
  }
}
cat(sprintf("\nOVERALL 1b IDENTITIES: %s\n", if (all_ok) "PASS" else "FAIL"))
