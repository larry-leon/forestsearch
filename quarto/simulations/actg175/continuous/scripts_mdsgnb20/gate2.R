# ===== GATE 2, per cell, campaign mdsgnb20 (FS, effMaxSG, eps 0.20, ACTG175 continuous / MD) =====
# Transplant of quarto/simulations/gbsg_020/scripts_p12x20/gate2F.R (the committed
# Gate 2 of p12x20: combine assertions, meta, finiteness, invariants, gamma,
# bound<->quantile identities, same-draws in two directions, payload size),
# pointed at the MD bundles with the checks of TASK_md_field_rerun_2026-09-15 §2.3:
#   - the combined bundle has 2,000 rows with sim_id exactly 1-2000, and the batch
#     files match it on every column;
#   - meta carries the campaign rule, field_scale_complement = "selected",
#     pkg_version 0.3.5 and host pop-os (both batch metas);
#   - same draws as mdf1, both directions: n_true identical and the oracle columns
#     within 1e-8 relative on all 2,000 sim_ids (null cell: complement oracle only);
#   - §1.6(b) field-s checks on every replicate whose complement field block was filled;
#   - MR failures on declared replicates <= max(20, 2 x mdf1's count in the cell).
# Identities are on the MD (identity) scale: est2_s + lam_mean_s == est2 + lam_mean
# (gate2F's log() form is the ratio-scale analogue).
# usage: Rscript gate2.R <md: 40|120|null> <n> <cell-label>; exit 1 on any failure.
QMD_DIR <- Sys.getenv("MDSG_DIR", unset = "..")
TAG     <- Sys.getenv("MDSG_TAG", unset = "mdsgnb20")
WORKERS <- Sys.getenv("MDSG_WORKERS", unset = NA)
HOST    <- Sys.getenv("MDSG_HOST", unset = "pop-os")
TOL <- 1e-8
`%||%` <- function(a,b) if (is.null(a) || length(a)==0 || all(is.na(a))) b else a
mdtok <- function(md) if (md == "null") "mdnull" else sprintf("md%s", md)
fstem <- function(md, n) sprintf("fs_effMaxSG_mr_field_%s_knoise0_n%d_nb20_%s", mdtok(md), n, TAG)
fdir  <- function(st) file.path(QMD_DIR, "mr_md_harm", paste0(st, "_d5000"))
mdf1p <- function(md, n) { st <- sprintf("fs_maxeffCons_mr_field_%s_knoise0_n%d_mdf1", mdtok(md), n)
  file.path(QMD_DIR, "mr_md_harm", paste0(st, "_d5000"), paste0(st, "_combined_1_2000.rds")) }

N_RUN <- 0L; N_PASS <- 0L; FAILED <- character(0)
P <- function(lab, ok, extra = "") {
  N_RUN <<- N_RUN + 1L
  if (isTRUE(ok)) N_PASS <<- N_PASS + 1L else FAILED <<- c(FAILED, lab)
  cat(sprintf("  %-62s %-8s %s\n", lab, if (isTRUE(ok)) "PASS" else "**FAIL**", extra))
}
relmax <- function(x, y) { d <- abs(x - y) / pmax(abs(y), 1e-300); d[is.na(x) & is.na(y)] <- 0; d[xor(is.na(x), is.na(y))] <- Inf; max(d) }

same_draws <- function(r, fr, label, null_cell) {
  cat(sprintf("  --- same-draws: %s ---\n", label))
  same_rows <- nrow(fr) == nrow(r) && identical(sort(fr$sim_id), sort(r$sim_id))
  P(sprintf("[%s] same sim_id set (%d rows)", label, nrow(r)), same_rows, sprintf("(comparator %d rows)", nrow(fr)))
  o <- order(r$sim_id); of <- order(fr$sim_id)
  nt_ok <- same_rows && identical(r$n_true[o], fr$n_true[of])
  P(sprintf("[%s] n_true identical() on all rows", label), nt_ok,
    if (!nt_ok && same_rows) sprintf("(%d of %d rows differ)", sum(r$n_true[o] != fr$n_true[of]), nrow(r)) else "")
  orc <- if (null_cell) c("or_Hc_est","or_Hc_lo","or_Hc_hi","or_Hc_se") else c("or_H_est","or_H_lo","or_H_hi","or_H_se","or_Hc_est","or_Hc_lo","or_Hc_hi","or_Hc_se")
  mx <- if (same_rows) max(vapply(orc, function(k) relmax(r[[k]][o], fr[[k]][of]), numeric(1))) else Inf
  P(sprintf("[%s] oracle columns (%s) <= %g relative on all rows", label, if (null_cell) "complement only" else "H and Hc"), mx <= TOL, sprintf("(max %.3g)", mx))
  invisible(NULL)
}

gate2 <- function(md, n, cell) {
  null_cell <- md == "null"
  st <- fstem(md, n); dr <- fdir(st); cf <- file.path(dr, paste0(st, "_combined_1_2000.rds"))
  cat(sprintf("\n########## GATE 2 (%s): %s -- md %s, n %d ##########\n", TAG, cell, md, n))
  cat(sprintf("  bundle: %s\n", cf))
  P("combined payload on disk", file.exists(cf))
  if (!file.exists(cf)) return(invisible(NULL))
  b <- readRDS(cf); r <- b$results; m <- b$meta

  # --- combine assertions (§2.3) ---
  cat("  --- combine assertions ---\n")
  f1 <- file.path(dr, paste0(st, "_res_1_1000.rds")); f2 <- file.path(dr, paste0(st, "_res_1001_2000.rds"))
  bts <- Sys.glob(file.path(dr, paste0(st, "_res_*.rds")))
  P("exactly 2 batch files, res_1_1000 and res_1001_2000", length(bts) == 2L && file.exists(f1) && file.exists(f2), sprintf("(%d)", length(bts)))
  b1 <- if (file.exists(f1)) readRDS(f1) else NULL; b2 <- if (file.exists(f2)) readRDS(f2) else NULL
  s1 <- if (!is.null(b1)) b1$results$sim_id else integer(0); s2 <- if (!is.null(b2)) b2$results$sim_id else integer(0)
  P("2,000 rows", nrow(r) == 2000L, sprintf("(%d)", nrow(r)))
  P("combined sim_id == 1:2000", identical(sort(r$sim_id), 1:2000))
  P("batch sim_id sets 1:1000 and 1001:2000", identical(sort(s1), 1:1000) && identical(sort(s2), 1001:2000), sprintf("(batch1 %d, batch2 %d)", length(s1), length(s2)))
  if (!is.null(b1) && !is.null(b2)) {
    rb <- rbind(b1$results, b2$results); rb <- rb[order(rb$sim_id), ]; rc <- r[order(r$sim_id), ]
    rownames(rb) <- NULL; rownames(rc) <- NULL
    same_cols <- identical(names(rb), names(rc))
    colsame <- if (same_cols) vapply(names(rc), function(k) identical(rb[[k]], rc[[k]]), logical(1)) else FALSE
    P("batch files match the combined bundle on every column", same_cols && all(colsame),
      if (same_cols && !all(colsame)) paste("differ:", paste(names(rc)[!colsame], collapse = ",")) else sprintf("(%d columns)", ncol(rc)))
  } else P("batch files match the combined bundle on every column", FALSE, "(batch file missing)")
  ce <- sum(grepl("CONFIG", as.character(r$status), ignore.case = TRUE))
  P("no CONFIG-ERROR replicate", ce == 0L, sprintf("(%d)", ce))

  # --- meta: both batch metas carry the campaign rule and constructions ---
  cat("  --- meta (both batches) ---\n")
  bm <- lapply(bts, function(f) readRDS(f)$meta)
  mk <- function(key, want) {
    got <- vapply(bm, function(x) paste(format(x[[key]] %||% "<absent>"), collapse = "|"), "")
    ok <- length(bm) == 2L && all(vapply(bm, function(x) isTRUE(all.equal(x[[key]], want)), TRUE))
    P(sprintf("meta: %s == %s", key, paste(format(want), collapse = "|")), ok, sprintf("(%s)", paste(got, collapse = " / ")))
  }
  mk("subgroup_method", "consistency"); mk("sg_focus", "effMaxSG"); mk("effect_neighborhood", 0.20)
  mk("selection_rule", "neighborhood"); mk("consistency_method", "resample")
  mk("ci_method", "field"); mk("mr_draws", 5000L); mk("field_uniform", FALSE); mk("field_complement", TRUE)
  mk("field_scale_complement", "selected"); mk("ij_residual", "two_term"); mk("return_reselection", TRUE)
  mk("fb_mode", "none"); mk("seed_base", 8316951L); mk("campaign_tag", TAG); mk("n_sample", as.integer(n))
  mk("null_cell", null_cell); mk("effect_threshold", 30); mk("consistency_threshold", 10)
  mk("pkg_version", "0.3.5"); mk("hostname", HOST)
  if (!is.na(WORKERS)) mk("n_workers", as.integer(WORKERS))
  cat(sprintf("  meta seed_base %s | host %s | R %s | built_at %s\n", m$seed_base,
              paste(unique(vapply(bm, function(x) x$hostname %||% NA_character_, "")), collapse="/"),
              paste(unique(vapply(bm, function(x) x$r_version %||% NA_character_, "")), collapse="/"),
              paste(vapply(bm, function(x) format(x$built_at %||% NA), ""), collapse=" / ")))

  # --- detection, MR failures on declared replicates ---
  det <- r$detected %in% 1L
  cat(sprintf("  >> DECLARATION RATE         : %.4f (%d / %d)\n", mean(det), sum(det), nrow(r)))
  D <- r[det, , drop = FALSE]
  mrfail <- sum(!(D$mr_ok %in% 1L) | !is.finite(D$mr_H_est))
  om <- readRDS(mdf1p(md, n)); Dm <- om$results[om$results$detected %in% 1L, ]
  mrfail_m <- sum(!(Dm$mr_ok %in% 1L) | !is.finite(Dm$mr_H_est))
  P(sprintf("MR failures on declared replicates <= max(20, 2 x mdf1's %d)", mrfail_m), mrfail <= max(20L, 2L * mrfail_m), sprintf("(%d)", mrfail))

  # --- finiteness of every product on filled replicates; field-s wiring (§1.6(b)) ---
  fin <- c("nv_H_est","nv_Hc_est","mr_H_est","mr_Hc_est","mr_H_lo","mr_H_hi","mr_Hc_lo","mr_Hc_hi",
           "fld_H_est2","fld_H_lo1s","fld_H_se","betaHhat_H","betaHhat_Hc","p_hat_H")
  miss <- setdiff(fin, names(r))
  P("every finiteness column present", length(miss) == 0L, if (length(miss)) paste("absent:", paste(miss, collapse=", ")) else "")
  Dm1 <- D[D$mr_ok %in% 1L & is.finite(D$fld_H_est2), ]
  nf <- fin[vapply(fin, function(k) any(!is.finite(Dm1[[k]])), logical(1))]
  P("harm products finite on declared replicates with a field block", length(nf) == 0L,
    if (length(nf)) paste("non-finite:", paste(nf, collapse=", ")) else sprintf("(%d rows)", nrow(Dm1)))
  f <- D[is.finite(D$fld_Hc_est2), ]
  sC <- c("fld_Hc_est2_s","fld_Hc_up1s_s","fld_Hc_lo1s_s","fld_Hc_lo2s_s","fld_Hc_hi2s_s","fld_Hc_lo_se_s","fld_Hc_hi_se_s","fld_Hc_se_s","fld_Hc_lam_mean_s")
  jC <- c("fld_joint_s_gamma","fld_joint_s_prob","fld_joint_s_loH","fld_joint_s_upHc","fld_joint_s_bonf_loH","fld_joint_s_bonf_upHc","fld_joint_s_bonf_prob","fld_joint_s_corr","fld_joint_s_n")
  P("nine fld_Hc_*_s and nine fld_joint_s_* columns present", all(c(sC, jC) %in% names(r)))
  P(sprintf("fld_Hc_*_s finite on all %d filled replicates", nrow(f)), nrow(f) > 0 && all(vapply(sC, function(k) all(is.finite(f[[k]])), logical(1))))
  P(sprintf("fld_joint_s_* finite on all %d filled replicates", nrow(f)), nrow(f) > 0 && all(vapply(jC, function(k) all(is.finite(f[[k]])), logical(1))))
  cat(sprintf("  complement field filled on %d of %d declared replicates (%d notes)\n", nrow(f), nrow(D), sum(!is.na(D$fld_Hc_note))))

  # --- interval invariants ---
  P("invariant harm  : fld_H_lo1s <= fld_H_est2",       all(Dm1$fld_H_lo1s <= Dm1$fld_H_est2))
  P("invariant compl : fld_Hc_est2 <= fld_Hc_up1s",     all(f$fld_Hc_est2 <= f$fld_Hc_up1s))
  P("invariant _s    : fld_Hc_lo1s_s <= fld_Hc_up1s_s", all(f$fld_Hc_lo1s_s <= f$fld_Hc_up1s_s))
  P("invariant _s    : fld_Hc_lo2s_s <= fld_Hc_hi2s_s", all(f$fld_Hc_lo2s_s <= f$fld_Hc_hi2s_s))
  P("invariant _s    : fld_Hc_est2_s <= fld_Hc_up1s_s", all(f$fld_Hc_est2_s <= f$fld_Hc_up1s_s))
  P("invariant joint : bonf_loH <= fld_H_est2",         all(f$fld_joint_bonf_loH <= f$fld_H_est2))
  P("invariant jointS: fld_Hc_est2_s <= bonf_upHc_s",   all(f$fld_Hc_est2_s <= f$fld_joint_s_bonf_upHc))
  P("invariant IJ    : mr_H_lo <= est <= mr_H_hi",      all(Dm1$mr_H_lo <= Dm1$mr_H_est & Dm1$mr_H_est <= Dm1$mr_H_hi))

  # --- gamma; identities (identity scale) ---
  g1 <- f$fld_joint_gamma; g2 <- f$fld_joint_s_gamma
  P("gamma (joint)   in [0.025, 0.05]", all(g1 >= 0.025 - 1e-12 & g1 <= 0.05 + 1e-12), sprintf("[%.5f, %.5f]", min(g1), max(g1)))
  P("gamma (joint-s) in [0.025, 0.05]", all(g2 >= 0.025 - 1e-12 & g2 <= 0.05 + 1e-12), sprintf("[%.5f, %.5f]", min(g2), max(g2)))
  i1 <- max(abs((f$fld_Hc_est2_s + f$fld_Hc_lam_mean_s) - (f$fld_Hc_est2 + f$fld_Hc_lam_mean)))
  P("identity: field-s inverted around the same beta-tilde^c", i1 <= 1e-9, sprintf("max |diff| = %.3g", i1))
  agree <- f$fld_joint_n == f$fld_joint_s_n
  dj <- if (any(agree)) max(abs(f$fld_joint_bonf_loH[agree] - f$fld_joint_s_bonf_loH[agree])) else 0
  P("identity: Bonferroni harm bound joint == joint_s where draw counts agree", dj <= 1e-12, sprintf("(%d of %d agree; max |diff| %.3g)", sum(agree), nrow(f), dj))
  i2 <- max(abs(Dm1$fld_H_lo1s - (Dm1$mr_H_est - Dm1$fld_H_q95)), abs(f$fld_Hc_up1s - (f$mr_Hc_est - f$fld_Hc_q05)))
  P("identity: lo1s = beta-tilde - q95; up1s = beta-tilde^c - q05", i2 <= 1e-9, sprintf("max |diff| = %.3g", i2))
  pv <- all(Dm1$p_hat_H >= 0 & Dm1$p_hat_H <= 1)
  P("p-hat in [0, 1]", pv, sprintf("(mean %.3f, share < 0.5: %.3f)", mean(Dm1$p_hat_H), mean(Dm1$p_hat_H < 0.5)))
  cat(sprintf("  >> CLASSIFICATION           : sens %.4f ppv %.4f | mean |Hhat| (n_harm) %.1f | mdf1 %.1f\n",
      mean(D$sens, na.rm = TRUE), mean(D$ppv, na.rm = TRUE), mean(D$n_harm, na.rm = TRUE), mean(Dm$n_harm, na.rm = TRUE)))

  # --- same-draws, both directions (§2.3) ---
  fr <- om$results
  same_draws(r, fr, "mdsgnb20 -> mdf1", null_cell)
  same_draws(fr, r, "mdf1 -> mdsgnb20", null_cell)

  # --- payload size ---
  for (fp in c(f1, f2, cf)) if (file.exists(fp)) {
    sz <- file.size(fp)
    P(sprintf("size <= 100 MB: %s", basename(fp)), sz <= 100 * 1024^2, sprintf("(%d B)%s", sz, if (sz > 50 * 1024^2) "  ** FLAG: over 50 MB **" else ""))
  }
  cat(sprintf("  timing: fit_mr_secs mean %.1f median %.1f p90 %.1f max %.1f | fld_H_secs mean %.1f | fld_Hc_secs mean %.2f\n",
              mean(r$fit_mr_secs), median(r$fit_mr_secs), quantile(r$fit_mr_secs, .9, names = FALSE), max(r$fit_mr_secs),
              mean(r$fld_H_secs, na.rm = TRUE), mean(r$fld_Hc_secs, na.rm = TRUE)))
  invisible(NULL)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("usage: Rscript gate2.R <md> <n> <cell>")
res <- tryCatch({ gate2(args[1], as.integer(args[2]), args[3]); TRUE },
                error = function(e) { P("checker ran without error", FALSE, conditionMessage(e)); FALSE })
cat(sprintf("\nGATE_COUNTS run=%d passed=%d failed=%d\n", N_RUN, N_PASS, N_RUN - N_PASS))
if (length(FAILED)) cat("FAILED:", paste(FAILED, collapse = " || "), "\n")
quit(status = if (N_RUN > N_PASS) 1L else 0L)
