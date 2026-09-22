#!/usr/bin/env Rscript
# GATE A -- nullmr (TASK_null_gbsg_mr_2026-09-21, Step 2.2), one per run.
# nullthr_gateA.R's invariants, except that its MR-absence checks invert:
#   - every MR product is NA wherever nothing was declared (and mr_ok == 0);
#   - wherever mr_ok == 1 the five evaluated products are finite (naive, IJ
#     two-term, field lower on beta(Hhat), field-s upper on beta(Hhat^c), the
#     field-s Bonferroni pair);
#   - every full-bootstrap product is NA;
#   - meta c1 0.90 / c2 0.80 / p* 0.90, and the resolved thresholds equal them
#     on EVERY replicate (not only where recorded);
#   - betaHhat_H and betaHhat_Hc finite on every declaring replicate.
# The mr_ok rate among declaring replicates is printed, not gated.
# Reads ONLY the bundle the run just wrote.  Invariants only.
#
# Three nullthr checks read columns the template fills ONLY on its MR-off
# branch (template :1376-1414): nv_H_lo1s (:1408-1409), p_sel and p_max_qual
# (:1391-1395).  With MR on they are structurally NA, so those checks become
# "NA on every replicate"; the unadjusted-estimate pairing keeps est and SE.
#
# usage: nullmr_gateA.R <engine> <hr> <n> <cell> <c1> <c2> <campaign> <reps> [quickrun TRUE|FALSE]
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
stem <- sprintf("results/%s_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d_null%03d_nb20%s_%s%s_res_1_%d.rds",
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
chk(all(!is.na(rc1) & abs(rc1 - C1) < 1e-12), "resolved c1 == c1 on every replicate",
    sprintf("(values: %s)", paste(unique(rc1[!is.na(rc1)]), collapse = ",")))
chk(all(!is.na(rc2) & abs(rc2 - C2) < 1e-12), "resolved c2 == c2 on every replicate",
    sprintf("(values: %s)", paste(unique(rc2[!is.na(rc2)]), collapse = ",")))
cat(sprintf("       (resolved thresholds recorded on %d / %d replicates for c1, %d / %d for c2)\n",
            sum(!is.na(rc1)), REPS, sum(!is.na(rc2)), REPS))
chk(all(is.finite(r$itt_est)), "itt_est finite on every replicate",
    sprintf("(%d of %d)", sum(is.finite(r$itt_est)), REPS))
chk(all(is.finite(r$itt_se) & r$itt_se > 0), "itt_se finite and positive on every replicate")

# --- MR ON (Step 2.2: the MR-absence checks, inverted) ----------------------
chk(isTRUE(m$mr_inference), "meta mr_inference == TRUE")
chk(identical(m$ci_method, "field"), "meta ci_method == field")
chk(isTRUE(m$mr_draws == 5000L), "meta mr_draws == 5000 (template :545)")
chk(isTRUE(m$field_complement) && isTRUE(m$field_decompose) &&
    identical(m$field_scale_complement, "selected") && identical(m$ij_residual, "two_term") &&
    isFALSE(m$field_uniform),
    "meta MR knobs: field_complement TRUE, field_decompose TRUE, scale_complement selected, ij two_term, uniform FALSE")
chk(identical(m$field_recovery, !identical(ENG, "consistency")),
    sprintf("meta field_recovery == %s (the committed %s setting)",
            !identical(ENG, "consistency"), if (identical(ENG, "consistency")) "FS" else "DINA / GRF"),
    sprintf("(%s)", format(m$field_recovery)))
chk(identical(m$fb_mode, "none"), "meta fb_mode == none (no bootstrap)")
chk(all(r$mr_ok %in% c(0L, 1L)), "mr_ok is 0 or 1 on every replicate")
chk(all(r$mr_ok[r$detected == 0L] == 0L), "mr_ok == 0 on every non-declaring replicate")
# MR products: every mr_* / fld_* column (mr_ok and the timing column
# fit_mr_secs are not products; nullid_gateA.R, 1f5c0076), plus the
# MR-only recorder fields n_family and p_hat_*.
mrcols <- grep("^(mr_|fld_)", names(r), value = TRUE)
mrcols <- setdiff(mrcols, c("mr_ok", "fit_mr_secs"))
mrcols <- c(mrcols, "n_family", grep("^p_hat_", names(r), value = TRUE))
nod0 <- r[r$detected == 0L, , drop = FALSE]
popn <- vapply(mrcols, function(k) any(!is.na(nod0[[k]])), logical(1))
chk(!any(popn), sprintf("all %d MR / field product columns are NA on every non-declaring replicate (%d rows)",
                        length(mrcols), nrow(nod0)),
    sprintf("(populated: %s)", paste(names(popn)[popn], collapse = ", ")))
ok1 <- r[r$mr_ok == 1L, , drop = FALSE]
FIVE <- list(
  naive          = c("nv_H_est", "nv_H_lo", "nv_H_hi", "nv_H_se", "nv_Hc_est", "nv_Hc_lo", "nv_Hc_hi", "nv_Hc_se"),
  "IJ two-term"  = c("mr_H_est", "mr_H_lo", "mr_H_hi", "mr_H_se_ij", "mr_Hc_est", "mr_Hc_lo", "mr_Hc_hi", "mr_Hc_se_ij"),
  "field lower"  = "fld_H_lo1s",
  "field-s upper" = "fld_Hc_up1s_s",
  "Bonferroni pair" = c("fld_joint_s_bonf_loH", "fld_joint_s_bonf_upHc"))
for (p in names(FIVE)) {
  bad <- vapply(FIVE[[p]], function(k) if (is.null(ok1[[k]])) -1L else sum(!is.finite(ok1[[k]])), integer(1))
  chk(all(bad == 0L), sprintf("%s finite wherever mr_ok == 1 (%s; %d rows)", p,
                              paste(FIVE[[p]], collapse = ", "), nrow(ok1)),
      sprintf("(non-finite: %s)", paste(names(bad)[bad != 0L], bad[bad != 0L], sep = "=", collapse = " ")))
}
fbcols <- setdiff(grep("^fb_", names(r), value = TRUE), c("fb_secs", "fb_err"))
fbpop  <- vapply(fbcols, function(k) any(!is.na(r[[k]])), logical(1))
chk(length(fbcols) > 0L && !any(fbpop), sprintf("all %d full-bootstrap product columns are NA", length(fbcols)),
    sprintf("(populated: %s)", paste(names(fbpop)[fbpop], collapse = ", ")))

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
chk(nrow(det) == 0L || identical(is.finite(det$nv_H_est), is.finite(det$nv_H_se)),
    "unadjusted estimate and SE are finite together on declaring replicates")
chk(all(is.na(r$nv_H_lo1s)), "nv_H_lo1s NA on every replicate (filled on the MR-off branch only, template :1408-1409)")
chk(nrow(det) == 0L || all(is.finite(det$betaHhat_H) & is.finite(det$betaHhat_Hc)),
    "betaHhat_H and betaHhat_Hc finite on every declaring replicate",
    sprintf("(%d / %d finite on both)", sum(is.finite(det$betaHhat_H) & is.finite(det$betaHhat_Hc)), nrow(det)))
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
  chk(all(is.na(r$p_sel)) && all(is.na(r$p_max_qual)),
      "p_sel / p_max_qual NA on every replicate (filled on the MR-off branch only, template :1391-1395)")
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
nmr <- sum(r$detected == 1L & r$mr_ok == 1L)
cat(sprintf("\n  mr_ok among declaring replicates: %d / %d (%s) -- printed, not gated\n",
            nmr, nd, if (nd) sprintf("%.4f", nmr / nd) else "-"))
cat(sprintf("\n  declarations %d / %d (%.4f) ; |Hhat| median %s ; mean search %.3f s\n",
            nd, REPS, nd / REPS,
            if (nd) format(stats::median(det$n_sel)) else "-", mean(r$fit_mr_secs, na.rm = TRUE)))
cat(sprintf("GATE A %s (%d checks, %d failed)\n",
            if (length(fails)) "FAILED" else "PASS", nchk, length(fails)))
quit(status = if (length(fails)) 1 else 0)
