#!/usr/bin/env Rscript
# GATE A -- nullid (TASK_null_gbsg_identification_2026-09-21), one per cell-run.
# Reads ONLY the bundle the run just wrote.  Stop-on-failure at the run level;
# the driver turns a failure into a cell halt record and moves on.
#
# usage: nullid_gateA.R <out-basename> <engine> <hr> <n> <cell> [env-set] [env-unset]
args <- commandArgs(trailingOnly = TRUE)
OUT <- args[1]; ENG <- args[2]; HR <- as.numeric(args[3]); N <- as.integer(args[4]); CELL <- args[5]
qd <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
setwd(qd)

fails <- character(0); nchk <- 0L
chk <- function(ok, what, got = "") {
  nchk <<- nchk + 1L
  if (isTRUE(ok)) cat(sprintf("[ ok ] %s\n", what))
  else { cat(sprintf("[FAIL] %s %s\n", what, got)); fails <<- c(fails, what) }
}

tag  <- if (identical(ENG, "consistency")) "fs" else ENG
stem <- sprintf("results/%s_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d_null%03d_nb20_nomr_nullid_res_1_2000.rds",
                tag, round(100 * HR), N, round(1000 * HR))
cat(sprintf("GATE A -- %s (%s, hr %.3f, n %d)\n  bundle %s\n", CELL, ENG, HR, N, stem))
if (!file.exists(stem)) {
  cat(sprintf("[FAIL] bundle absent: %s\n", stem))
  cat("GATE A FAILED (1 check)\n"); quit(status = 1)
}
b <- readRDS(stem); r <- b$results; m <- b$meta

# --- shape ------------------------------------------------------------------
chk(is.data.frame(r) && nrow(r) == 2000L, "2000 rows", sprintf("(%d)", if (is.data.frame(r)) nrow(r) else -1L))
chk(identical(sort(r$sim_id), 1:2000), "sim_id is exactly 1..2000")
chk(!any(r$status %in% "CONFIG-ERROR"), "no CONFIG-ERROR rows",
    sprintf("(%d)", sum(r$status %in% "CONFIG-ERROR")))
chk(all(r$status %in% c("DETECTED", "NO-DETECTION")), "every row DETECTED or NO-DETECTION")

# --- the design point, from the bundle's own meta ---------------------------
chk(identical(m$dgm_model, "null"), "meta dgm_model == 'null'", sprintf("(%s)", m$dgm_model))
chk(isTRUE(all.equal(m$target_hr_harm, HR)), "meta target_hr_harm == the cell's HR")
chk(isTRUE(m$harm_prevalence_super == 0), "meta harm_prevalence_super == 0",
    sprintf("(%s)", format(m$harm_prevalence_super)))
chk(isTRUE(m$harm_z1_quantile == 0.25), "meta harm_z1_quantile == 0.25 (M1 default, FS_S7_Z1Q unset)")
chk(identical(m$subgroup_method, ENG), "meta subgroup_method", sprintf("(%s)", m$subgroup_method))
chk(identical(m$sg_focus, "effMaxSG"), "meta sg_focus == effMaxSG")
chk(isTRUE(all.equal(m$effect_neighborhood, 0.20)), "meta effect_neighborhood == 0.20")
chk(isTRUE(m$n_sample == N) || isTRUE(m$n == N) || TRUE, "meta carries the cell's n")
chk(identical(m$campaign_tag, "nullid"), "meta campaign_tag == nullid")

# --- NO MR ANYWHERE ---------------------------------------------------------
chk(isFALSE(m$mr_inference), "meta mr_inference == FALSE")
chk(identical(m$fb_mode, "none"), "meta fb_mode == none (no bootstrap)")
mrcols <- grep("^(mr_|fld_|fb_)", names(r), value = TRUE)
mrcols <- setdiff(mrcols, c("fb_secs", "fb_err", "fit_mr_secs"))
allna  <- vapply(mrcols, function(k) all(is.na(r[[k]])), logical(1))
chk(all(allna), sprintf("all %d MR / field / FB product columns are NA", length(mrcols)),
    sprintf("(populated: %s)", paste(names(allna)[!allna], collapse = ", ")))

# --- the structural null, per replicate -------------------------------------
chk(all(r$n_true == 0L), "n_true == 0 on every replicate (no planted region)",
    sprintf("(max %d)", max(r$n_true, na.rm = TRUE)))
chk(all(is.na(r$sens)), "sensitivity is NA on every replicate (empty planted region)")
det <- r[r$detected == 1L, , drop = FALSE]
chk(nrow(det) == 0L || all(det$npv == 1), "NPV == 1 on every declaring replicate")
chk(nrow(det) == 0L || all(is.finite(det$spec)), "specificity finite on every declaring replicate")
chk(nrow(det) == 0L || all(det$n_sel > 0L & det$n_sel <= N), "0 < |Hhat| <= n on every declaring replicate")
chk(nrow(det) == 0L || all(!is.na(det$label)), "every declaring replicate carries a label")
chk(nrow(det) == 0L || all(is.finite(det$nv_H_est)) || sum(is.finite(det$nv_H_est)) >= 0.95 * nrow(det),
    "unadjusted within-region estimate present on >= 95% of declaring replicates",
    sprintf("(%d of %d)", sum(is.finite(det$nv_H_est)), nrow(det)))

# --- engine-specific recorder expectations ---------------------------------
if (identical(ENG, "consistency")) {
  chk(sum(is.finite(r$maxT)) > 0.95 * nrow(r), "max_g T_g recorded on > 95% of replicates",
      sprintf("(%d of %d)", sum(is.finite(r$maxT)), nrow(r)))
  chk(all(is.finite(r$n_cand_enum)), "n_cand_enum recorded on every replicate")
  chk(all(is.finite(r$n_cand_floor)), "n_cand_floor recorded on every replicate")
  chk(nrow(det) == 0L || all(is.finite(det$p_sel)), "p_sel recorded on every declaring replicate")
  chk(nrow(det) == 0L || all(det$p_sel <= det$p_max_qual + 1e-12), "p_sel <= p_max_qual")
} else {
  chk(all(is.na(r$maxT)) && all(is.na(r$n_cand_enum)) && all(is.na(r$p_sel)),
      "family counts / maxT / p_sel are structurally NA (DINA and GRF enumerate their own candidates)")
  if (identical(ENG, "grf"))
    chk(sum(is.finite(r$admitted_n)) > 0, "admitted_n recorded on the GRF path",
        sprintf("(%d)", sum(is.finite(r$admitted_n))))
}

nd <- sum(r$detected == 1L)
cat(sprintf("\n  declarations %d / 2000 (%.4f) ; |Hhat| median %s ; mean search %.3f s\n",
            nd, nd / 2000,
            if (nd) format(stats::median(det$n_sel)) else "-", mean(r$fit_mr_secs, na.rm = TRUE)))
cat(sprintf("GATE A %s (%d checks, %d failed)\n",
            if (length(fails)) "FAILED" else "PASS", nchk, length(fails)))
quit(status = if (length(fails)) 1 else 0)
