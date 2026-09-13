# Portability header, as gate2.R: P12X20_QMD_DIR overrides the directory holding
# results/.  P12X20_TAG / P12X20_Z1Q_TAG exist only so the checker can be
# dry-run against a committed 31% bundle; the campaign leaves both unset.
QMD_DIR <- Sys.getenv("P12X20_QMD_DIR", unset = "..")
TAG     <- Sys.getenv("P12X20_TAG", unset = "p12x20")
Z1Q_TAG <- Sys.getenv("P12X20_Z1Q_TAG", unset = "")
WORKERS <- Sys.getenv("P12X20_WORKERS", unset = NA)

# ===== GATE 2, per cell, campaign p12x20 (FS, effMaxSG, eps 0.20, 12.4%) =====
# gate2.R repointed at the FS cells (TASK_p12x20_partA_2026-09-12_v2 §5-§6), as
# gate2G.R repointed it at GRF.  Same content -- completeness; detection;
# realized prevalence; finiteness of every product on detected replicates; the
# interval invariants; gamma in [0.025, 0.05]; the CORRECTED bound<->quantile
# identity log(est2_s) + lam_mean_s == log(est2) + lam_mean; the Bonferroni
# identity gated on the gamma-at-floor rows (NOT the fld_joint_bonf_* vs
# fld_joint_* pair on all rows); p-hat validity -- with these changes:
#
#   ADDED.  §5 combine assertions first: 2,000 rows, no duplicate sim_id, the
#   two batches' sim_id sets disjoint and their union 1:2000.
#   ADDED.  §6 same-draws in TWO directions: against the committed 12.4% FS
#   bundle (p12ext / tier2, gate2.R's block-A map) AND against dinamr and
#   grfmr at 12.4%.  Each: n_true identical() on all 2,000 rows; truth
#   all.equal() at 1e-8 (asserted); identical() on truth and the maximum
#   absolute and relative discrepancy beside it (reported).
#   ADDED.  Payload size: > 100 MB is a hard-stop failure; > 50 MB is flagged.
#   META.  subgroup_method 'consistency', campaign 'p12x20', and the Stage 1
#   §1b knobs read off BOTH batch metas.  field_recovery is asserted FALSE
#   (FS_S7_FIELD_RECOV unset, as cert20), so the nine recovery columns are
#   REPORTED, not required, and are excluded from the finiteness list.
#   CONSISTENCY SCREEN.  FS computes Pcons: n_cons_qual / band_n / p_star are
#   reported as present / populated, never labelled structural.
#
# Exit status 1 if any check fails (or the checker errors); the runner halts.
TOL_TRUTH <- 1e-8
`%||%` <- function(a,b) if (is.null(a) || length(a)==0 || all(is.na(a))) b else a
R <- file.path(QMD_DIR, "results/")
fstem <- function(hr, n) sprintf("%sfs_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d%s_nb20_%s",
                                 R, round(100*hr), n, Z1Q_TAG, TAG)
# gate2.R's block-A (12.4%) FS map, unchanged: tier2 at HR 1.75 and at HR 1.00
# n 500, p12ext otherwise -- maxeffCons, eps 0.10.
fscomp <- function(hr, n) {
  camp <- if (abs(hr - 1.75) < 1e-9) "tier2" else if (abs(hr - 1.00) < 1e-9 && n == 500L) "tier2" else "p12ext"
  sprintf("%sfs_maxeffCons_fb_mr_field_m1_h%03d_knoise0_n%d_%s_combined_1_2000.rds", R, round(100*hr), n, camp)
}
engcomp <- function(hr, n, eng) {
  camp <- c(dina = "dinamr", grf = "grfmr")[[eng]]
  sprintf("%s%s_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d%s_nb20_%s_combined_1_2000.rds",
          R, eng, round(100*hr), n, Z1Q_TAG, camp)
}

N_RUN <- 0L; N_PASS <- 0L; FAILED <- character(0)
P <- function(lab, ok, extra = "") {
  N_RUN <<- N_RUN + 1L
  if (isTRUE(ok)) N_PASS <<- N_PASS + 1L else FAILED <<- c(FAILED, lab)
  cat(sprintf("  %-58s %-8s %s\n", lab, if (isTRUE(ok)) "PASS" else "**FAIL**", extra))
}

same_draws <- function(b, r, fp, label) {
  cat(sprintf("  --- same-draws: %s ---\n", label))
  P(sprintf("[%s] comparator on disk", label), file.exists(fp), sprintf("(%s)", basename(fp)))
  if (!file.exists(fp)) return(invisible(NULL))
  fb <- readRDS(fp); fr <- fb$results
  cat(sprintf("    comparator meta: method %s, focus %s, eps %s, campaign %s, host %s\n",
      fb$meta$subgroup_method %||% NA, fb$meta$sg_focus %||% NA, format(fb$meta$effect_neighborhood %||% NA),
      fb$meta$campaign_tag %||% NA, fb$meta$hostname %||% NA))
  same_rows <- nrow(fr) == nrow(r) && identical(sort(fr$sim_id), sort(r$sim_id))
  P(sprintf("[%s] same sim_id set (%d rows)", label, nrow(r)), same_rows, sprintf("(comparator %d rows)", nrow(fr)))
  o <- order(r$sim_id); of <- order(fr$sim_id)
  nt_ok <- same_rows && identical(r$n_true[o], fr$n_true[of])
  P(sprintf("[%s] n_true identical() on all rows", label), nt_ok,
    if (!nt_ok && same_rows) sprintf("(%d of %d rows differ)", sum(r$n_true[o] != fr$n_true[of]), nrow(r)) else "")
  ae <- all.equal(b$truth, fb$truth, tolerance = TOL_TRUTH)
  P(sprintf("[%s] truth all.equal(tol = %g)", label, TOL_TRUTH), isTRUE(ae),
    if (!isTRUE(ae)) paste(ae, collapse = "; ") else "")
  ta <- unlist(b$truth); tb <- unlist(fb$truth)
  cat(sprintf("    truth identical() [reported, not asserted]: %s\n", identical(b$truth, fb$truth)))
  cat(sprintf("    truth max |abs diff| %.4g ; max |rel diff| %.4g\n",
      max(abs(ta - tb)), max(abs(ta - tb) / pmax(abs(tb), .Machine$double.xmin))))
  if (!nt_ok || !isTRUE(ae))
    cat("    >> DRAWS DO NOT MATCH. A FINDING ABOUT THE DGM PATH: record and STOP (task §6).\n")
  invisible(NULL)
}

gate2 <- function(hr, n, cell) {
  st <- fstem(hr, n); cf <- paste0(st, "_combined_1_2000.rds")
  cat(sprintf("\n########## GATE 2 (p12x20): %s -- HR %.2f, n %d ##########\n", cell, hr, n))
  cat(sprintf("  bundle: %s\n", basename(cf)))
  P("combined payload on disk", file.exists(cf))
  if (!file.exists(cf)) return(invisible(NULL))
  b <- readRDS(cf); r <- b$results; m <- b$meta

  # --- §5 combine assertions ---
  cat("  --- combine assertions (task §5) ---\n")
  f1 <- paste0(st, "_res_1_1000.rds"); f2 <- paste0(st, "_res_1001_2000.rds")
  bts <- Sys.glob(paste0(st, "_res_*.rds"))
  P("exactly 2 batch files, res_1_1000 and res_1001_2000",
    length(bts) == 2L && file.exists(f1) && file.exists(f2), sprintf("(%d)", length(bts)))
  s1 <- if (file.exists(f1)) readRDS(f1)$results$sim_id else integer(0)
  s2 <- if (file.exists(f2)) readRDS(f2)$results$sim_id else integer(0)
  P("2,000 rows", nrow(r) == 2000L, sprintf("(%d)", nrow(r)))
  P("no duplicate sim_id", !any(duplicated(r$sim_id)), sprintf("(%d duplicated)", sum(duplicated(r$sim_id))))
  P("batch sim_id sets disjoint", length(intersect(s1, s2)) == 0L,
    sprintf("(batch1 %d, batch2 %d, overlap %d)", length(s1), length(s2), length(intersect(s1, s2))))
  P("batch sim_id union == 1:2000", identical(sort(union(s1, s2)), 1:2000))
  P("combined sim_id == 1:2000", identical(sort(r$sim_id), 1:2000))
  ce <- if ("status" %in% names(r)) sum(grepl("CONFIG", as.character(r$status), ignore.case = TRUE)) else 0L
  P("no CONFIG-ERROR replicate", ce == 0L, sprintf("(%d)", ce))

  # --- meta: both batch metas carry the Stage 1 §1b knob set ---
  cat("  --- meta (both batches) ---\n")
  bm <- lapply(bts, function(f) readRDS(f)$meta)
  mk <- function(key, want) {
    got <- vapply(bm, function(x) paste(format(x[[key]] %||% "<absent>"), collapse = "|"), "")
    ok <- length(bm) == 2L && all(vapply(bm, function(x) isTRUE(all.equal(x[[key]], want)), TRUE))
    P(sprintf("meta: %s == %s", key, paste(format(want), collapse = "|")), ok, sprintf("(%s)", paste(got, collapse = " / ")))
  }
  mk("subgroup_method", "consistency"); mk("sg_focus", "effMaxSG"); mk("effect_neighborhood", 0.20)
  mk("selection_rule", "neighborhood"); mk("stop_threshold", "NULL")
  mk("harm_z1_quantile", if (nzchar(Z1Q_TAG)) 0.60 else 0.25); mk("er_jcuts", 10L)
  mk("mr_inference", TRUE); mk("ci_method", "field"); mk("mr_draws", 5000L)
  mk("field_uniform", FALSE); mk("field_complement", TRUE); mk("field_decompose", TRUE)
  mk("field_scale_complement", "selected"); mk("field_recovery", FALSE); mk("ij_residual", "two_term")
  mk("fb_mode", "none"); mk("nb_boots", 0L); mk("seed_base", 8316951L); mk("campaign_tag", TAG)
  mk("target_hr_harm", hr); mk("n_sample", as.integer(n))
  if (!is.na(WORKERS)) mk("n_workers", as.integer(WORKERS))
  P("meta: forestsearch_version recorded", all(vapply(bm, function(x) !is.null(x$forestsearch_version), TRUE)),
    sprintf("(%s)", paste(unique(vapply(bm, function(x) x$forestsearch_version %||% NA_character_, "")), collapse = "/")))
  cat(sprintf("  meta seed_base %s | host %s | R %s | built_at %s\n", m$seed_base,
              paste(unique(vapply(bm, function(x) x$hostname %||% NA_character_, "")), collapse="/"),
              paste(unique(vapply(bm, function(x) x$r_version %||% NA_character_, "")), collapse="/"),
              paste(vapply(bm, function(x) format(x$built_at %||% NA), ""), collapse=" / ")))

  # --- detection, family, prevalence (recorded) ---
  det <- mean(r$detected %in% 1L); K <- r$n_family[is.finite(r$n_family)]
  cat(sprintf("  >> %-26s: %.4f (%d / %d)\n", if (abs(hr - 1.00) < 1e-9) "SELECTION RATE" else "DETECTION RATE",
      det, sum(r$detected %in% 1L), nrow(r)))
  if (length(K))
    cat(sprintf("  >> N_FAMILY                 : min %g  q10 %g  med %g  q90 %g  max %g  (mean %.1f, CV %.3f, n %d)\n",
        min(K), quantile(K,.1,names=FALSE), median(K), quantile(K,.9,names=FALSE), max(K), mean(K), sd(K)/mean(K), length(K)))
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
           "betaHhat_H","betaHhat_Hc","p_hat_H","p_hat_sum","p_hat_top1")
  miss <- setdiff(fin, names(r))
  P("every finiteness column present", length(miss) == 0L, if (length(miss)) paste("absent:", paste(miss, collapse=", ")) else "")
  fin <- fin[fin %in% names(r)]
  nf <- fin[vapply(fin, function(k) any(!is.finite(D[[k]])), logical(1))]
  P("every product finite on detected replicates", length(nf) == 0L,
    if (length(nf)) paste("non-finite:", paste(nf, collapse=", ")) else sprintf("(%d quantities)", length(fin)))
  cat(sprintf("  recovery columns [reported; FS_S7_FIELD_RECOV unset per §1b]: present %d/9 ; populated on detected %d/9\n",
      sum(recov %in% names(r)), sum(vapply(recov, function(k) k %in% names(r) && !all(is.na(D[[k]])), logical(1)))))
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
  P("gamma (joint)   in [0.025, 0.05]", all(g1 >= 0.025 & g1 <= 0.05), sprintf("[%.5f, %.5f]", min(g1), max(g1)))
  P("gamma (joint-s) in [0.025, 0.05]", all(g2 >= 0.025 & g2 <= 0.05), sprintf("[%.5f, %.5f]", min(g2), max(g2)))

  # --- bound <-> quantile identities (the corrected, cert20 form) ---
  i1 <- max(abs((log(D$fld_Hc_est2_s) + D$fld_Hc_lam_mean_s) -
                (log(D$fld_Hc_est2)   + D$fld_Hc_lam_mean)), na.rm = TRUE)
  P("identity: field-s inverted around the same bdc", i1 <= 1e-12, sprintf("max |diff| = %.3g", i1))
  atfl <- D$fld_joint_gamma == 0.025; atfs <- D$fld_joint_s_gamma == 0.025
  j1 <- if (any(atfl)) max(abs(D$fld_joint_bonf_loH[atfl]  - D$fld_joint_loH[atfl]),
                           abs(D$fld_joint_bonf_upHc[atfl] - D$fld_joint_upHc[atfl])) else 0
  j2 <- if (any(atfs)) max(abs(D$fld_joint_s_bonf_loH[atfs]  - D$fld_joint_s_loH[atfs]),
                           abs(D$fld_joint_s_bonf_upHc[atfs] - D$fld_joint_s_upHc[atfs])) else 0
  P("identity: joint bonf == raw where gamma at floor", max(j1, j2) <= 1e-12,
    sprintf("max |diff| = %.3g (share at floor %.3f / %.3f)", max(j1, j2), mean(atfl), mean(atfs)))
  pv <- all(D$p_hat_H >= 0 & D$p_hat_H <= 1) && all(D$p_hat_sum >= 0) &&
        all(D$p_hat_H <= D$p_hat_sum + 1e-12) && all(D$p_hat_H <= D$p_hat_top1 + 1e-12)
  P("p-hat validity (0<=p<=1, p_H<=sum, p_H<=top1)", pv)
  cat(sprintf("  >> CLASSIFICATION           : sens %.4f spec %.4f ppv %.4f npv %.4f | mean |Hhat| %.1f\n",
      mean(D$sens), mean(D$spec), mean(D$ppv), mean(D$npv), mean(D$n_sel)))

  # --- non-detections ---
  ND <- r[!(r$detected %in% 1L), , drop = FALSE]
  if (nrow(ND)) {
    cat(sprintf("  >> NON-DETECTIONS           : %d (%.4f). status %s | err_msg %d | n_family NA %d\n",
        nrow(ND), nrow(ND)/nrow(r), paste(unique(as.character(ND$status)), collapse = "/"),
        sum(!is.na(ND$err_msg)), sum(is.na(ND$n_family))))
  } else cat("  >> NON-DETECTIONS           : none\n")

  # --- consistency-screen columns: FS computes Pcons; reported as present ---
  for (k in c("n_cons_qual","band_n","p_star"))
    cat(sprintf("  CONSISTENCY-SCREEN %-12s %s\n", k,
      if (!k %in% names(r)) "not a recorder column"
      else sprintf("present, non-NA %d / %d rows (detected: %d / %d)", sum(!is.na(r[[k]])), nrow(r),
                   sum(!is.na(D[[k]])), nrow(D))))

  # --- same-draws, two directions (task §6) ---
  cat("  --- same-draws direction 1: the committed 12.4% FS bundle ---\n")
  same_draws(b, r, fscomp(hr, n), "FS p12ext/tier2")
  cat("  --- same-draws direction 2: dinamr and grfmr at 12.4% ---\n")
  same_draws(b, r, engcomp(hr, n, "dina"), "dinamr")
  same_draws(b, r, engcomp(hr, n, "grf"), "grfmr")

  # --- payload size ---
  cat("  --- payload size (100 MB hard stop; 50 MB flag) ---\n")
  for (f in c(f1, f2, cf)) if (file.exists(f)) {
    sz <- file.size(f)
    P(sprintf("size <= 100 MB: %s", basename(f)), sz <= 100 * 1024^2,
      sprintf("(%d B)%s", sz, if (sz > 50 * 1024^2) "  ** FLAG: over 50 MB **" else ""))
  }
  invisible(NULL)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("usage: Rscript gate2F.R <hr> <n> <cell>")
res <- tryCatch({ gate2(as.numeric(args[1]), as.integer(args[2]), args[3]); TRUE },
                error = function(e) { P("checker ran without error", FALSE, conditionMessage(e)); FALSE })
cat(sprintf("\nGATE_COUNTS run=%d passed=%d failed=%d\n", N_RUN, N_PASS, N_RUN - N_PASS))
if (length(FAILED)) cat("FAILED:", paste(FAILED, collapse = " || "), "\n")
quit(status = if (N_RUN > N_PASS) 1L else 0L)
