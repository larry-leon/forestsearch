#!/usr/bin/env Rscript
# GATE A -- nullc125 (TASK_null_gbsg_thresholds_2026-09-21_v2), one per run.
# nullid_gateA.R's invariants, plus the threshold and ITT checks of Step 2.2.
# Reads ONLY the bundle the run just wrote.  Invariants only -- never a guessed
# rate, count or share: a strict screen that rarely declares is an outcome.
#
# One departure from nullid_gateA.R, recorded in the report: its check
# "unadjusted within-region estimate present on >= 95% of declaring
# replicates" was a guessed share.  It is replaced by the invariant behind it
# (est, SE and one-sided bound are finite together on declaring rows, and all
# NA on non-declaring rows); the coverage is printed, not gated.
#
# usage: nullthr_gateA.R <engine> <hr> <n> <cell> <c1> <c2> <campaign> <reps> [quickrun TRUE|FALSE]
args <- commandArgs(trailingOnly = TRUE)
ENG <- args[1]; HR <- as.numeric(args[2]); N <- as.integer(args[3]); CELL <- args[4]
C1 <- as.numeric(args[5]); C2 <- as.numeric(args[6]); CAMP <- args[7]; REPS <- as.integer(args[8])
QR <- identical(if (length(args) >= 9) args[9] else "FALSE", "TRUE")
qd <- normalizePath(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)[1])), ".."))
setwd(qd)

fails <- character(0); nchk <- 0L
chk <- function(ok, what, got = "") {
  nchk <<- nchk + 1L
  if (isTRUE(ok)) cat(sprintf("[ ok ] %s\n", what))
  else { cat(sprintf("[FAIL] %s %s\n", what, got)); fails <<- c(fails, what) }
}

tag  <- if (identical(ENG, "consistency")) "fs" else ENG
thr_tag <- if (abs(C1 - 0.90) < 1e-12 && abs(C2 - 0.80) < 1e-12) "" else
  sprintf("_c%03dc%03d", round(100 * C1), round(100 * C2))
stem <- sprintf("results/%s_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d_null%03d_nb20_nomr%s_%s%s_res_1_%d.rds",
                tag, round(100 * HR), N, round(1000 * HR), thr_tag, CAMP,
                if (QR) "_quickrun" else "", REPS)
cat(sprintf("GATE A -- %s (%s, hr %.3f, n %d, c1 %.2f, c2 %.2f, %d reps)\n  bundle %s\n",
            CELL, ENG, HR, N, C1, C2, REPS, stem))
if (!file.exists(stem)) {
  cat(sprintf("[FAIL] bundle absent: %s\n", stem))
  cat("GATE A FAILED (1 check)\n"); quit(status = 1)
}
b <- readRDS(stem); r <- b$results; m <- b$meta

# --- shape ------------------------------------------------------------------
chk(is.data.frame(r) && nrow(r) == REPS, sprintf("%d rows", REPS),
    sprintf("(%d)", if (is.data.frame(r)) nrow(r) else -1L))
chk(identical(sort(r$sim_id), seq_len(REPS)), sprintf("sim_id is exactly 1..%d", REPS))
chk(!any(r$status %in% "CONFIG-ERROR"), "no CONFIG-ERROR rows",
    sprintf("(%d)", sum(r$status %in% "CONFIG-ERROR")))
chk(all(r$status %in% c("DETECTED", "NO-DETECTION")), "every row DETECTED or NO-DETECTION")
chk(identical(r$status == "DETECTED", r$detected == 1L), "status DETECTED exactly where detected == 1")

# --- the design point, from the bundle's own meta ---------------------------
chk(identical(m$dgm_model, "null"), "meta dgm_model == 'null'", sprintf("(%s)", m$dgm_model))
chk(isTRUE(all.equal(m$target_hr_harm, HR)), "meta target_hr_harm == the cell's HR")
chk(isTRUE(m$harm_prevalence_super == 0), "meta harm_prevalence_super == 0",
    sprintf("(%s)", format(m$harm_prevalence_super)))
chk(isTRUE(m$harm_z1_quantile == 0.25), "meta harm_z1_quantile == 0.25 (M1 default, FS_S7_Z1Q unset)")
chk(identical(m$subgroup_method, ENG), "meta subgroup_method", sprintf("(%s)", m$subgroup_method))
chk(identical(m$sg_focus, "effMaxSG"), "meta sg_focus == effMaxSG")
chk(isTRUE(all.equal(m$effect_neighborhood, 0.20)), "meta effect_neighborhood == 0.20")
chk(isTRUE(m$er_jcuts == 10L), "meta er_jcuts == 10 (FS_S7_ER_JCUTS unset)")
chk(identical(m$campaign_tag, CAMP), sprintf("meta campaign_tag == %s", CAMP),
    sprintf("(%s)", m$campaign_tag))

# --- thresholds (Step 2.2) --------------------------------------------------
chk(isTRUE(all.equal(m$c1, C1)), sprintf("meta c1 == %.2f", C1), sprintf("(%s)", format(m$c1)))
chk(isTRUE(all.equal(m$c2, C2)), sprintf("meta c2 == %.2f", C2), sprintf("(%s)", format(m$c2)))
chk(isTRUE(all.equal(m$pstar, 0.90)), "meta pstar == 0.90 (the literal)", sprintf("(%s)", format(m$pstar)))
chk(isTRUE(all.equal(m$dmin_grf, 0.0)), "meta dmin_grf == 0.0", sprintf("(%s)", format(m$dmin_grf)))
chk(isTRUE(all.equal(m$dina_m_diff, log(C1))), "meta dina_m_diff == log(c1)",
    sprintf("(%s)", format(m$dina_m_diff)))
rc1 <- r$c1_resolved; rc2 <- r$c2_resolved
chk(all(is.na(rc1) | abs(rc1 - C1) < 1e-12), "resolved c1 == the knob on every replicate where recorded",
    sprintf("(values: %s)", paste(unique(rc1[!is.na(rc1)]), collapse = ",")))
chk(all(is.na(rc2) | abs(rc2 - C2) < 1e-12), "resolved c2 == the knob on every replicate where recorded",
    sprintf("(values: %s)", paste(unique(rc2[!is.na(rc2)]), collapse = ",")))
cat(sprintf("       (resolved thresholds recorded on %d / %d replicates for c1, %d / %d for c2)\n",
            sum(!is.na(rc1)), REPS, sum(!is.na(rc2)), REPS))
chk(all(is.finite(r$itt_est)), "itt_est finite on every replicate",
    sprintf("(%d of %d)", sum(is.finite(r$itt_est)), REPS))
chk(all(is.finite(r$itt_se) & r$itt_se > 0), "itt_se finite and positive on every replicate")

# --- NO MR ANYWHERE ---------------------------------------------------------
chk(isFALSE(m$mr_inference), "meta mr_inference == FALSE")
chk(all(r$mr_ok == 0L), "mr_ok == 0 on every replicate (MR never ran)")
chk(identical(m$fb_mode, "none"), "meta fb_mode == none (no bootstrap)")
mrcols <- grep("^(mr_|fld_|fb_)", names(r), value = TRUE)
# mr_ok, fb_secs, fb_err and fit_mr_secs match the prefix but are not MR
# PRODUCTS (nullid_gateA.R, fixed in 1f5c0076).
mrcols <- setdiff(mrcols, c("mr_ok", "fb_secs", "fb_err", "fit_mr_secs"))
allna  <- vapply(mrcols, function(k) all(is.na(r[[k]])), logical(1))
chk(all(allna), sprintf("all %d MR / field / FB product columns are NA", length(mrcols)),
    sprintf("(populated: %s)", paste(names(allna)[!allna], collapse = ", ")))

# --- the structural null, per replicate -------------------------------------
chk(all(r$n_true == 0L), "n_true == 0 on every replicate (no planted region)",
    sprintf("(max %d)", max(r$n_true, na.rm = TRUE)))
chk(all(is.na(r$sens)), "sensitivity is NA on every replicate (empty planted region)")
det <- r[r$detected == 1L, , drop = FALSE]
nod <- r[r$detected == 0L, , drop = FALSE]
chk(nrow(det) == 0L || all(det$npv == 1), "NPV == 1 on every declaring replicate")
chk(nrow(det) == 0L || all(is.finite(det$spec)), "specificity finite on every declaring replicate")
chk(nrow(det) == 0L || all(det$n_sel > 0L & det$n_sel <= N), "0 < |Hhat| <= n on every declaring replicate")
chk(nrow(det) == 0L || all(!is.na(det$label)), "every declaring replicate carries a label")
chk(nrow(det) == 0L || (identical(is.finite(det$nv_H_est), is.finite(det$nv_H_se)) &&
                        identical(is.finite(det$nv_H_est), is.finite(det$nv_H_lo1s))),
    "unadjusted estimate, SE and one-sided bound are finite together on declaring replicates")
chk(nrow(nod) == 0L || all(is.na(nod$nv_H_est)), "no unadjusted estimate on non-declaring replicates")
cat(sprintf("       (unadjusted within-region estimate present on %d / %d declaring replicates)\n",
            sum(is.finite(det$nv_H_est)), nrow(det)))

# --- engine-specific recorder expectations ---------------------------------
if (identical(ENG, "consistency")) {
  chk(identical(is.finite(r$maxT), r$n_cand_floor > 0L),
      sprintf("max_g T_g is finite exactly when the screened family is non-empty (%d of %d)",
              sum(is.finite(r$maxT)), nrow(r)))
  chk(sum(r$n_cand_floor == 0L & r$detected == 1L) == 0L,
      "no replicate declares a region when nothing cleared the effect floor",
      sprintf("(%d)", sum(r$n_cand_floor == 0L & r$detected == 1L)))
  chk(all(is.finite(r$n_cand_enum)), "n_cand_enum recorded on every replicate")
  chk(all(is.finite(r$n_cand_floor)), "n_cand_floor recorded on every replicate")
  chk(nrow(det) == 0L || all(is.finite(det$p_sel)), "p_sel recorded on every declaring replicate")
  chk(nrow(det) == 0L || all(det$p_sel <= det$p_max_qual + 1e-12), "p_sel <= p_max_qual")
  chk(nrow(det) == 0L || all(det$p_sel >= 0.90 - 1e-12), "p_sel >= p* on every declaring replicate")
} else {
  chk(all(is.na(r$maxT)) && all(is.na(r$n_cand_enum)) && all(is.na(r$p_sel)),
      "family counts / maxT / p_sel are structurally NA (DINA and GRF enumerate their own candidates)")
  # nullid asserted admitted_n was recorded on SOME replicate (> 0), a count.
  # The invariant: a declared region was admitted, so admitted_n >= 1 there.
  if (identical(ENG, "grf"))
    chk(nrow(det) == 0L || all(is.finite(det$admitted_n) & det$admitted_n >= 1L),
        "admitted_n >= 1 on every declaring GRF replicate",
        sprintf("(recorded on %d of %d replicates)", sum(is.finite(r$admitted_n)), nrow(r)))
}

nd <- sum(r$detected == 1L)
cat(sprintf("\n  declarations %d / %d (%.4f) ; |Hhat| median %s ; mean search %.3f s\n",
            nd, REPS, nd / REPS,
            if (nd) format(stats::median(det$n_sel)) else "-", mean(r$fit_mr_secs, na.rm = TRUE)))
cat(sprintf("GATE A %s (%d checks, %d failed)\n",
            if (length(fails)) "FAILED" else "PASS", nchk, length(fails)))
quit(status = if (length(fails)) 1 else 0)
