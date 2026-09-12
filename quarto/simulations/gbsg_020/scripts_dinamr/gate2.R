# Portability header added when this script was committed: the campaign ran with
# SCRATCH = the session scratchpad and QMD_DIR = quarto/simulations/gbsg_020.
SCRATCH <- Sys.getenv("DINAMR_SCRATCH", unset = ".")
QMD_DIR <- Sys.getenv("DINAMR_QMD_DIR", unset = "..")

# ===== GATE 2, per cell, campaign dinamr =====================================
# Completeness; detection rate; proposed-family distribution; realized prevalence;
# finiteness of every product on detected replicates; interval invariants;
# gamma in [0.025, 0.05]; bound<->quantile identities <= 1e-12; structurally-NA
# columns reported AS SUCH; and AMENDMENT 3 -- the same-draws assertions against
# the committed FS comparator that shares the DGM draws (n_true identical() on
# all 2,000 rows; truth by all.equal at tolerance 1e-8, NOT identical(), because
# the FS bundles are Linux-produced and cross-machine BLAS moves truth at ~1e-16).
TOL_TRUTH <- 1e-8
R <- file.path(QMD_DIR, "results/")
dstem <- function(hr, n, blk) sprintf("%sdina_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d%s_nb20_dinamr",
                                      R, round(100*hr), n, if (blk == "B") "_z1q60" else "")
# The committed FS comparator that shares the DGM draws, per cell.
# Block A (12.4%): p12ext at HR 1.50/1.00, tier2 at HR 1.75 -- maxeffCons, eps 0.10.
# Block B (31%)  : cert20 at n 1000/1500, e1stud at n 500  -- effMaxSG, eps 0.20.
fscomp <- function(hr, n, blk) {
  if (blk == "A") {
    camp <- if (abs(hr - 1.75) < 1e-9) "tier2" else if (abs(hr - 1.00) < 1e-9 && n == 500L) "tier2" else "p12ext"
    sprintf("%sfs_maxeffCons_fb_mr_field_m1_h%03d_knoise0_n%d_%s_combined_1_2000.rds", R, round(100*hr), n, camp)
  } else {
    # Block C addendum (TASK_dinamr_blockC_grfprobe_2026-09-11): e1stud was run
    # only at HR 1.50 and 1.75, so at HR 1.00 the n = 500 / 31% comparator is
    # cert20, which IS on disk and is criterion-matched (effMaxSG, eps 0.20).
    camp <- if (n == 500L && abs(hr - 1.00) > 1e-9) "e1stud" else "cert20"
    sprintf("%sfs_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d_z1q60_nb20_%s_combined_1_2000.rds", R, round(100*hr), n, camp)
  }
}
gate2 <- function(hr, n, blk) {
  st <- dstem(hr, n, blk); cf <- paste0(st, "_combined_1_2000.rds")
  nm <- sprintf("Block %s -- HR %.2f, n %d", blk, hr, n)
  cat("\n########## GATE 2: ", nm, " ##########\n", sep = "")
  if (!file.exists(cf)) { cat("  NOT ON DISK (deferred or failed):", basename(cf), "\n"); return(invisible(NULL)) }
  b <- readRDS(cf); r <- b$results; m <- b$meta
  P <- function(lab, ok, extra = "") cat(sprintf("  %-46s %-6s %s\n", lab, if (isTRUE(ok)) "PASS" else "**FAIL**", extra))

  # --- completeness ---
  P("2,000 rows", nrow(r) == 2000L, sprintf("(%d)", nrow(r)))
  P("sim_id 1..2000 exactly, no duplicates",
    identical(sort(r$sim_id), 1:2000) && !any(duplicated(r$sim_id)))
  ce <- if ("status" %in% names(r)) sum(grepl("CONFIG", as.character(r$status), ignore.case = TRUE)) else 0L
  P("no CONFIG-ERROR replicate", ce == 0L, sprintf("(%d)", ce))
  # batch metas carry n_workers / forestsearch_version
  bts <- Sys.glob(paste0(st, "_res_*.rds")); bm <- lapply(bts, function(f) readRDS(f)$meta)
  P("meta: subgroup_method == 'dina'", identical(m$subgroup_method, "dina"))
  P("meta: sg_focus == 'effMaxSG'", identical(m$sg_focus, "effMaxSG"))
  P("meta: effect_neighborhood == 0.20", isTRUE(all.equal(m$effect_neighborhood, 0.20)))
  P("meta: field_complement TRUE", isTRUE(m$field_complement))
  P("meta: field_decompose TRUE", isTRUE(m$field_decompose))
  P("meta: field_scale_complement 'selected'", identical(m$field_scale_complement, "selected"))
  P("meta: field_recovery TRUE", isTRUE(m$field_recovery))
  P("meta: ij_residual 'two_term'", identical(m$ij_residual, "two_term"))
  P("meta: campaign_tag 'dinamr'", identical(m$campaign_tag, "dinamr"))
  P("meta: 2 seed-disjoint batches", length(bts) == 2L, sprintf("(%d)", length(bts)))
  P("meta: n_workers recorded in every batch",
    all(vapply(bm, function(x) !is.null(x$n_workers), TRUE)),
    sprintf("(%s)", paste(unique(vapply(bm, function(x) as.character(x$n_workers %||% NA), "")), collapse = "/")))
  P("meta: forestsearch_version recorded",
    all(vapply(bm, function(x) !is.null(x$forestsearch_version), TRUE)),
    sprintf("(%s)", paste(unique(vapply(bm, function(x) x$forestsearch_version %||% NA_character_, "")), collapse = "/")))
  cat(sprintf("  meta seed_base %s | host %s | R %s\n", m$seed_base,
              paste(unique(vapply(bm, function(x) x$hostname %||% NA_character_, "")), collapse="/"),
              paste(unique(vapply(bm, function(x) x$r_version %||% NA_character_, "")), collapse="/")))

  # --- detection, family, prevalence (recorded prominently) ---
  det <- mean(r$detected %in% 1L); K <- r$n_family[is.finite(r$n_family)]
  cat(sprintf("  >> DETECTION RATE           : %.4f (%d / %d)\n", det, sum(r$detected %in% 1L), nrow(r)))
  cat(sprintf("  >> PROPOSED-FAMILY SIZE     : min %g  q10 %g  med %g  q90 %g  max %g  (mean %.1f, CV %.3f)\n",
      min(K), quantile(K,.1,names=FALSE), median(K), quantile(K,.9,names=FALSE), max(K), mean(K), sd(K)/mean(K)))
  cat(sprintf("  >> REALIZED PREVALENCE      : super-population %.5f | trial mean %.5f\n",
      m$harm_prevalence_super, mean(r$n_true)/m$n_sample))

  D <- r[r$detected %in% 1L, , drop = FALSE]
  recov <- c("fld_recov_sens_H","fld_recov_ppv_H","fld_recov_sens_Hc","fld_recov_npv_Hc",
             "fld_recov_q10","fld_recov_q50","fld_recov_q90","fld_recov_share1","fld_recov_n_used")
  fin <- c("nv_H_est","nv_Hc_est","mr_H_est","mr_Hc_est","mr_H_lo","mr_H_hi","mr_Hc_lo","mr_Hc_hi",
           "fld_H_est2","fld_H_lo1s","fld_H_se","fld_Hc_est2","fld_Hc_up1s","fld_Hc_se",
           "fld_Hc_est2_s","fld_Hc_up1s_s","fld_Hc_se_s","fld_Hc_scale_ratio",
           "fld_joint_gamma","fld_joint_bonf_loH","fld_joint_bonf_upHc",
           "fld_joint_s_gamma","fld_joint_s_bonf_loH","fld_joint_s_bonf_upHc",
           "betaHhat_H","betaHhat_Hc","p_hat_H","p_hat_sum","p_hat_top1", recov)
  fin <- fin[fin %in% names(r)]
  nf <- fin[vapply(fin, function(k) any(!is.finite(D[[k]])), logical(1))]
  P("every product finite on detected replicates", length(nf) == 0L,
    if (length(nf)) paste("non-finite:", paste(nf, collapse=", ")) else sprintf("(%d quantities)", length(fin)))
  P("nine recovery columns present & populated",
    all(recov %in% names(r)) && !any(vapply(recov, function(k) all(is.na(D[[k]])), logical(1))))
  P("p-hat block present & populated",
    all(c("p_hat_H","p_hat_sum","p_hat_top1") %in% names(r)) && !all(is.na(D$p_hat_H)))
  P("rho-c (fld_Hc_scale_ratio) present & populated",
    "fld_Hc_scale_ratio" %in% names(r) && !all(is.na(D$fld_Hc_scale_ratio)))

  # --- interval invariants ---
  P("invariant harm  : fld_H_lo1s <= fld_H_est2",       all(D$fld_H_lo1s <= D$fld_H_est2))
  P("invariant compl : fld_Hc_est2 <= fld_Hc_up1s",     all(D$fld_Hc_est2 <= D$fld_Hc_up1s))
  P("invariant _s    : fld_Hc_est2_s <= fld_Hc_up1s_s", all(D$fld_Hc_est2_s <= D$fld_Hc_up1s_s))
  P("invariant joint : bonf_loH <= fld_H_est2",         all(D$fld_joint_bonf_loH <= D$fld_H_est2))
  P("invariant joint : fld_Hc_est2 <= bonf_upHc",       all(D$fld_Hc_est2 <= D$fld_joint_bonf_upHc))
  P("invariant jointS: bonf_loH_s <= fld_H_est2",       all(D$fld_joint_s_bonf_loH <= D$fld_H_est2))
  P("invariant jointS: fld_Hc_est2_s <= bonf_upHc_s",   all(D$fld_Hc_est2_s <= D$fld_joint_s_bonf_upHc))
  P("invariant IJ    : mr_H_lo <= est <= mr_H_hi",      all(D$mr_H_lo <= D$mr_H_est & D$mr_H_est <= D$mr_H_hi))
  P("invariant IJ    : mr_Hc_lo <= est <= mr_Hc_hi",    all(D$mr_Hc_lo <= D$mr_Hc_est & D$mr_Hc_est <= D$mr_Hc_hi))
  P("invariant rho-c : scale_ratio > 0",                all(D$fld_Hc_scale_ratio > 0))

  # --- gamma ---
  g1 <- D$fld_joint_gamma; g2 <- D$fld_joint_s_gamma
  P("gamma (joint)   in [0.025, 0.05]", all(g1 >= 0.025 & g1 <= 0.05),
    sprintf("[%.5f, %.5f]", min(g1), max(g1)))
  P("gamma (joint-s) in [0.025, 0.05]", all(g2 >= 0.025 & g2 <= 0.05),
    sprintf("[%.5f, %.5f]", min(g2), max(g2)))

  # --- bound <-> quantile identities (the cert20 form) -----------------------
  # THE identity that matters for the studentized block, quoted from
  # REPORT_cert20_2026-09-08: field-s is inverted around the SAME bdc as the
  # unscaled complement field, log(est2_s) + lam_mean_s == log(est2) + lam_mean.
  i1 <- max(abs((log(D$fld_Hc_est2_s) + D$fld_Hc_lam_mean_s) -
                (log(D$fld_Hc_est2)   + D$fld_Hc_lam_mean)), na.rm = TRUE)
  P("identity: field-s inverted around the same bdc", i1 <= 1e-12,
    sprintf("max |diff| = %.3g", i1))
  # The Bonferroni joint bound equals the raw joint bound EXACTLY WHEN gamma
  # sits at its 0.025 floor, and is a different quantity otherwise -- so the
  # identity is gated on the floor rows, and the share at the floor is reported
  # beside the gamma range, as cert20 did.
  atfl <- D$fld_joint_gamma == 0.025; atfs <- D$fld_joint_s_gamma == 0.025
  j1 <- if (any(atfl)) max(abs(D$fld_joint_bonf_loH[atfl]  - D$fld_joint_loH[atfl]),
                           abs(D$fld_joint_bonf_upHc[atfl] - D$fld_joint_upHc[atfl])) else 0
  j2 <- if (any(atfs)) max(abs(D$fld_joint_s_bonf_loH[atfs]  - D$fld_joint_s_loH[atfs]),
                           abs(D$fld_joint_s_bonf_upHc[atfs] - D$fld_joint_s_upHc[atfs])) else 0
  P("identity: joint bonf == raw where gamma at floor", max(j1, j2) <= 1e-12,
    sprintf("max |diff| = %.3g (share at floor %.3f / %.3f)", max(j1, j2), mean(atfl), mean(atfs)))
  # p-hat validity
  pv <- all(D$p_hat_H >= 0 & D$p_hat_H <= 1) && all(D$p_hat_sum >= 0) &&
        all(D$p_hat_H <= D$p_hat_sum + 1e-12) && all(D$p_hat_H <= D$p_hat_top1 + 1e-12)
  P("p-hat validity (0<=p<=1, p_H<=sum, p_H<=top1)", pv)
  # classification, recorded as cert20 did
  cat(sprintf("  >> CLASSIFICATION           : sens %.4f spec %.4f ppv %.4f npv %.4f | mean |Hhat| %.1f\n",
      mean(D$sens), mean(D$spec), mean(D$ppv), mean(D$npv), mean(D$n_sel)))

  # --- NON-DETECTIONS, split as far as the bundle allows --------------------
  ND <- r[!(r$detected %in% 1L), , drop = FALSE]
  if (nrow(ND)) {
    cat(sprintf("  >> NON-DETECTIONS           : %d (%.4f). n_family NA %d | == 0 %d | > 0 %d | err_msg %d | status %s\n",
        nrow(ND), nrow(ND)/nrow(r), sum(is.na(ND$n_family)),
        sum(!is.na(ND$n_family) & ND$n_family == 0), sum(!is.na(ND$n_family) & ND$n_family > 0),
        sum(!is.na(ND$err_msg)), paste(unique(as.character(ND$status)), collapse = "/")))
    cat("     SPLIT (empty proposed family vs candidates proposed but none admitted): NOT AVAILABLE\n")
    cat("     The recorder returns its all-NA record at `if (!found) { rec$status <- \"NO-DETECTION\"; return(rec) }`\n")
    cat("     BEFORE n_family is written; n_family is read from the MR gate object, which does not\n")
    cat("     exist when nothing is selected.  So n_family is NA on every non-detection -- never 0,\n")
    cat("     never positive -- and the two causes are not separable from the committed columns.\n")
  } else cat("  >> NON-DETECTIONS           : none\n")

  # --- structurally-NA, reported AS SUCH (never as a failure) ---
  for (k in c("n_cons_qual","band_n","p_star"))
    cat(sprintf("  STRUCTURAL %-12s %s\n", k,
      if (!k %in% names(r)) "not a recorder column (admission-set term; NULL on DINA)"
      else if (all(is.na(r[[k]]))) "present, all-NA -- STRUCTURAL on DINA (no consistency screen), NOT a failure"
      else "POPULATED (unexpected on DINA)"))

  # --- AMENDMENT 3: same-draws assertions ---
  fp <- fscomp(hr, n, blk)
  cat("  --- Amendment 3 (same-draws vs the committed FS comparator) ---\n")
  if (!file.exists(fp)) {
    cat("    FS comparator not on disk:", basename(fp), " -- assertion not evaluable.\n")
  } else {
    fb <- readRDS(fp); fr <- fb$results
    same_rows <- nrow(fr) == nrow(r) && identical(sort(fr$sim_id), sort(r$sim_id))
    o <- order(r$sim_id); of <- order(fr$sim_id)
    nt_ok <- same_rows && identical(r$n_true[o], fr$n_true[of])
    ae <- all.equal(b$truth, fb$truth, tolerance = TOL_TRUTH)
    tr_ok <- isTRUE(ae)
    cat(sprintf("    comparator: %s  (focus %s, eps %s, campaign %s)\n", basename(fp),
        fb$meta$sg_focus %||% NA, format(fb$meta$effect_neighborhood %||% NA), fb$meta$campaign_tag %||% NA))
    cat(sprintf("    n_true identical() on all %d rows : %s%s\n", nrow(r),
        if (nt_ok) "YES" else "NO",
        if (!nt_ok && same_rows) sprintf("   (%d of %d rows differ)", sum(r$n_true[o] != fr$n_true[of]), nrow(r)) else ""))
    cat(sprintf("    truth all.equal(tol = %g)          : %s%s\n", TOL_TRUTH,
        if (tr_ok) "YES" else "NO", if (!tr_ok) paste0("  -> ", paste(ae, collapse="; ")) else ""))
    cat(sprintf("    truth identical() [reported, not asserted]: %s\n", identical(b$truth, fb$truth)))
    ta <- unlist(b$truth); tb <- unlist(fb$truth)
    cat(sprintf("    truth max |abs diff| %.4g ; max |rel diff| %.4g\n",
        max(abs(ta - tb)), max(abs(ta - tb) / pmax(abs(tb), .Machine$double.xmin))))
    FK <- fr$n_family[is.finite(fr$n_family)]
    cat(sprintf("    FS enumerated family (same cell): min %g med %g q90 %g max %g CV %.4f | detection %.4f\n",
        min(FK), stats::median(FK), stats::quantile(FK,.9,names=FALSE), max(FK),
        stats::sd(FK)/mean(FK), mean(fr$detected %in% 1L)))
    cat(sprintf("    DINA family (this cell)        : min %g med %g q90 %g max %g CV %.4f | detection %.4f\n",
        min(K), stats::median(K), stats::quantile(K,.9,names=FALSE), max(K),
        stats::sd(K)/mean(K), det))
    if (!nt_ok || !tr_ok)
      cat("    >> DRAWS DO NOT MATCH. Recorded as a FINDING ABOUT THE DGM PATH, not a cell failure.\n")
  }
  invisible(NULL)
}
`%||%` <- function(a,b) if (is.null(a) || length(a)==0 || all(is.na(a))) b else a
args <- commandArgs(trailingOnly = TRUE)
blocks <- if (length(args)) args else c("A","B","C")
if ("A" %in% blocks) for (hr in c(1.50,1.75)) for (n in c(500L,1000L,1500L)) gate2(hr, n, "A")
if ("B" %in% blocks) for (hr in c(1.50,1.75)) for (n in c(500L,1000L,1500L)) gate2(hr, n, "B")
if ("C" %in% blocks) for (blk in c("A","B")) for (n in c(500L,1000L,1500L)) gate2(1.00, n, blk)
