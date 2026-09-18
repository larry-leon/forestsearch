# ===== GATE 2, per cell, ACTG175 binary / OR campaigns orfs | orgrf | ordina =====
# Transplant of ../../continuous/scripts_mddina/gate2.R (combine assertions, meta, finiteness,
# invariants, gamma, identities, same-draws in two directions, payload size), pointed at the OR
# bundles (stem <fs|grf|dina>_effMaxSG_mr_field_or0XX_n<N>_nb20_<tag>) with
# TASK_actg175_binary_campaign_2026-09-17 Stage 2's per-cell Gate 2:
#   - 2,000 rows with sim_id exactly 1-2000, and the batch files matching the combined bundle on
#     every column;
#   - meta carries the rule (effMaxSG, eps 0.20, neighborhood), the thresholds (0.90 / 0.80 /
#     0.90 with adverse_outcome = TRUE), target_or_h, the truths, field_scale_complement,
#     pkg_version and host;
#   - SAME DRAWS ACROSS IDENTIFIERS: for orgrf and ordina, the data-level columns -- the truths,
#     the oracle columns and the seed -- are identical to orfs's in the same cell, in BOTH
#     directions, within 1e-8 relative;
#   - the §1.5(d) construction checks on every declared replicate;
#   - MR failures on declared replicates at most 40.
# EVERY coverage number the orgrf and ordina campaigns produce is coverage of beta(H-hat)
# CONDITIONAL ON THE PROPOSED FAMILY: GRF's and DINA's families are generated from fitted
# surfaces, so the fixed-family condition does not hold.  FS's family is the prespecified cut
# grid.  Comparisons across the three are descriptive, not a contest.
# n_cons_qual / band_n are FS-only (no consistency screen on GRF or DINA) and are reported as
# structurally NA there, not as failures; admitted_n is GRF's forest-qualified count / DINA's
# admitted count and is NA on FS.
# usage: Rscript gate2.R <target: 0.75|1.0|1.5> <n> <cell-label> <tag: orfs|orgrf|ordina>
QMD_DIR <- Sys.getenv("ORSG_DIR", unset = "..")
WORKERS <- Sys.getenv("ORSG_WORKERS", unset = NA)
HOST    <- Sys.getenv("ORSG_HOST", unset = "pop-os")
PKG     <- Sys.getenv("ORSG_PKG",  unset = "0.3.5")
MRFAIL  <- as.integer(Sys.getenv("ORSG_MRFAIL", unset = "40"))
TOL <- 1e-8
`%||%` <- function(a,b) if (is.null(a) || length(a)==0 || all(is.na(a))) b else a
dtok  <- function(target) sprintf("or%03d", as.integer(round(100 * as.numeric(target))))
mtag  <- function(tag) switch(tag, orfs = "fs", orgrf = "grf", ordina = "dina",
                              stop("unknown campaign tag: ", tag))
mfull <- function(tag) switch(tag, orfs = "consistency", orgrf = "grf", ordina = "dina")
fstem <- function(target, n, tag)
  sprintf("%s_effMaxSG_mr_field_%s_n%d_nb20_%s", mtag(tag), dtok(target), n, tag)
fdir  <- function(st) file.path(QMD_DIR, "mr_or_harm", paste0(st, "_d5000"))
fcomb <- function(target, n, tag) { st <- fstem(target, n, tag)
  file.path(fdir(st), paste0(st, "_combined_1_2000.rds")) }

N_RUN <- 0L; N_PASS <- 0L; FAILED <- character(0)
P <- function(lab, ok, extra = "") {
  N_RUN <<- N_RUN + 1L
  if (isTRUE(ok)) N_PASS <<- N_PASS + 1L else FAILED <<- c(FAILED, lab)
  cat(sprintf("  %-70s %-8s %s\n", lab, if (isTRUE(ok)) "PASS" else "**FAIL**", extra))
}
relmax <- function(x, y) { d <- abs(x - y) / pmax(abs(y), 1e-300)
  d[is.na(x) & is.na(y)] <- 0; d[xor(is.na(x), is.na(y))] <- Inf; max(d) }

# The data-level columns: rule-independent by construction (the seed is the table lookup, n_true
# the true-region size, and the oracle refits on the TRUE region), so two identifiers running the
# same cell MUST agree on them exactly.
DATA_COLS <- c("seed", "n_true",
               "or_H_est","or_H_lo","or_H_hi","or_H_se",
               "or_Hc_est","or_Hc_lo","or_Hc_hi","or_Hc_se")
TRUTH_KEYS <- c("or_causal", "marg_H", "marg_Hc", "cde_H", "cde_Hc",
                "prevalence_Q", "beta_inter")

same_draws <- function(r, fr, tr, ftr, label) {
  cat(sprintf("  --- same draws as orfs: %s ---\n", label))
  same_rows <- nrow(fr) == nrow(r) && identical(sort(fr$sim_id), sort(r$sim_id))
  P(sprintf("[%s] same sim_id set (%d rows)", label, nrow(r)), same_rows,
    sprintf("(comparator %d rows)", nrow(fr)))
  o <- order(r$sim_id); of <- order(fr$sim_id)
  if (same_rows) {
    mx <- vapply(DATA_COLS, function(k) relmax(as.numeric(r[[k]][o]), as.numeric(fr[[k]][of])), numeric(1))
    P(sprintf("[%s] the data-level columns (%d: seed, n_true, the 8 oracle columns) <= %g relative",
              label, length(DATA_COLS), TOL), max(mx) <= TOL,
      sprintf("(max %.3g%s)", max(mx),
              if (max(mx) > TOL) paste0(" on ", paste(names(mx)[mx > TOL], collapse = ",")) else ""))
    P(sprintf("[%s] seed identical() on all rows", label),
      identical(as.integer(r$seed[o]), as.integer(fr$seed[of])))
  } else {
    P(sprintf("[%s] the data-level columns <= %g relative", label, TOL), FALSE, "(row sets differ)")
    P(sprintf("[%s] seed identical() on all rows", label), FALSE, "(row sets differ)")
  }
  tmx <- vapply(TRUTH_KEYS, function(k) relmax(tr[[k]] %||% NA_real_, ftr[[k]] %||% NA_real_), numeric(1))
  P(sprintf("[%s] the truths (%s) <= %g relative", label, paste(TRUTH_KEYS, collapse = ","), TOL),
    max(tmx) <= TOL, sprintf("(max %.3g)", max(tmx)))
  invisible(NULL)
}

gate2 <- function(target, n, cell, tag) {
  st <- fstem(target, n, tag); dr <- fdir(st); cf <- fcomb(target, n, tag)
  meth <- mfull(tag)
  cat(sprintf("\n########## GATE 2 (%s): %s -- target_or_h %s, n %d, identifier %s ##########\n",
              tag, cell, target, n, meth))
  cat(if (identical(meth, "consistency"))
        "  FS's candidate family is the prespecified cut grid: the fixed-family condition holds.\n" else
        "  Every coverage figure of this campaign is coverage of beta(H-hat) CONDITIONAL ON THE PROPOSED FAMILY.\n")
  cat(sprintf("  bundle: %s\n", cf))
  P("combined payload on disk", file.exists(cf))
  if (!file.exists(cf)) return(invisible(NULL))
  b <- readRDS(cf); r <- b$results; m <- b$meta; tr <- b$truth

  cat("  --- combine assertions ---\n")
  f1 <- file.path(dr, paste0(st, "_res_1_1000.rds")); f2 <- file.path(dr, paste0(st, "_res_1001_2000.rds"))
  bts <- Sys.glob(file.path(dr, paste0(st, "_res_*.rds")))
  P("exactly 2 batch files, res_1_1000 and res_1001_2000",
    length(bts) == 2L && file.exists(f1) && file.exists(f2), sprintf("(%d)", length(bts)))
  b1 <- if (file.exists(f1)) readRDS(f1) else NULL; b2 <- if (file.exists(f2)) readRDS(f2) else NULL
  s1 <- if (!is.null(b1)) b1$results$sim_id else integer(0)
  s2 <- if (!is.null(b2)) b2$results$sim_id else integer(0)
  P("2,000 rows", nrow(r) == 2000L, sprintf("(%d)", nrow(r)))
  P("combined sim_id == 1:2000 exactly", identical(sort(as.integer(r$sim_id)), 1:2000))
  P("batch sim_id sets 1:1000 and 1001:2000",
    identical(sort(as.integer(s1)), 1:1000) && identical(sort(as.integer(s2)), 1001:2000),
    sprintf("(batch1 %d, batch2 %d)", length(s1), length(s2)))
  if (!is.null(b1) && !is.null(b2)) {
    rb <- rbind(b1$results, b2$results); rb <- rb[order(rb$sim_id), ]; rc <- r[order(r$sim_id), ]
    rownames(rb) <- NULL; rownames(rc) <- NULL
    same_cols <- identical(names(rb), names(rc))
    colsame <- if (same_cols) vapply(names(rc), function(k) identical(rb[[k]], rc[[k]]), logical(1)) else FALSE
    P("batch files match the combined bundle on every column", same_cols && all(colsame),
      if (same_cols && !all(colsame)) paste("differ:", paste(names(rc)[!colsame], collapse = ","))
      else sprintf("(%d columns)", ncol(rc)))
  } else P("batch files match the combined bundle on every column", FALSE, "(batch file missing)")
  ce <- sum(grepl("CONFIG", as.character(r$status), ignore.case = TRUE))
  P("no CONFIG-ERROR replicate", ce == 0L, sprintf("(%d)", ce))

  cat("  --- meta (both batches) ---\n")
  bm <- lapply(bts, function(f) readRDS(f)$meta)
  mk <- function(key, want) {
    got <- vapply(bm, function(x) paste(format(x[[key]] %||% "<absent>"), collapse = "|"), "")
    ok <- length(bm) == 2L && all(vapply(bm, function(x) isTRUE(all.equal(x[[key]], want)), TRUE))
    P(sprintf("meta: %s == %s", key, paste(format(want), collapse = "|")), ok,
      sprintf("(%s)", paste(got, collapse = " / ")))
  }
  mk("subgroup_method", meth)
  mk("sg_focus", "effMaxSG"); mk("effect_neighborhood", 0.20); mk("selection_rule", "neighborhood")
  mk("effect_threshold", 0.90); mk("consistency_threshold", 0.80); mk("pconsistency", 0.90)
  mk("adverse_outcome", TRUE); mk("outcome_type", "binary"); mk("effect_measure", "OR")
  mk("target_or_h", as.numeric(target)); mk("design_tag", dtok(target)); mk("dgm_model", "alt")
  mk("sg_quantile", 0.70); mk("n_super", 100000L); mk("eval_seed", 20260628L)
  mk("ci_method", "field"); mk("mr_draws", 5000L); mk("field_uniform", FALSE)
  mk("field_complement", TRUE); mk("field_scale_complement", "selected")
  mk("ij_residual", "two_term"); mk("return_reselection", TRUE); mk("fb_mode", "none")
  mk("seed_base", 8316951L); mk("campaign_tag", tag); mk("n_sample", as.integer(n))
  mk("k_random_noise", 0L); mk("consistency_method", "resample")
  if (identical(meth, "grf")) {
    mk("dmin_grf", 0.0); mk("grf_selection", "frontier")
    mk("grf_depth", 2L); mk("grf_select_statistic", "effect")
  }
  if (identical(meth, "dina")) { mk("dina_select_statistic", "effect"); mk("dina_args", "list()") }
  mk("pkg_version", PKG); mk("hostname", HOST)
  if (!is.na(WORKERS)) mk("n_workers", as.integer(WORKERS))
  P("meta carries the truths", all(vapply(paste0("truth_", c("or_causal","marg_H","marg_Hc","cde_H","cde_Hc","prevalence_Q","beta_inter")),
                                          function(k) is.finite(m[[k]] %||% NA_real_), logical(1))),
    sprintf("(marg_H %.10f | marg_Hc %.10f | cde_H %.10f | cde_Hc %.10f | prev %.6f)",
            m$truth_marg_H, m$truth_marg_Hc, m$truth_cde_H, m$truth_cde_Hc, m$truth_prevalence_Q))
  cat(sprintf("  meta seed_base %s | seed_scheme %s | host %s | R %s | pkg_commit %s | built_at %s\n",
              m$seed_base, m$seed_scheme,
              paste(unique(vapply(bm, function(x) x$hostname %||% NA_character_, "")), collapse="/"),
              paste(unique(vapply(bm, function(x) x$r_version %||% NA_character_, "")), collapse="/"),
              m$pkg_commit %||% NA_character_,
              paste(vapply(bm, function(x) format(x$built_at %||% NA), ""), collapse=" / ")))

  det <- r$detected %in% 1L
  cat(sprintf("  >> DECLARATION RATE         : %.4f (%d / %d)\n", mean(det), sum(det), nrow(r)))
  D <- r[det, , drop = FALSE]
  mrfail <- sum(!(D$mr_ok %in% 1L) | !is.finite(D$mr_H_est))
  P(sprintf("MR failures on declared replicates <= %d", MRFAIL), mrfail <= MRFAIL, sprintf("(%d)", mrfail))

  # --- identifier accounting ---
  A <- r$admitted_n[is.finite(r$admitted_n)]; K <- D$n_family[is.finite(D$n_family)]
  if (length(K)) cat(sprintf("  >> N_FAMILY (MR's kept family K): min %g  med %g  p90 %g  max %g\n",
      min(K), stats::median(K), stats::quantile(K, .9, names = FALSE), max(K)))
  if (length(A)) { q <- stats::quantile(A, c(.25,.50,.75,.90), names = FALSE)
    cat(sprintf("  >> ADMITTED_N               : min %g  q25 %g  MED %g  q75 %g  p90 %g  max %g  (mean %.1f, n %d); == 0 on %d row(s)\n",
        min(A), q[1], q[2], q[3], q[4], max(A), mean(A), length(A), sum(A == 0L))) }
  if (identical(meth, "dina")) {
    PR <- D$dina_proposed_n[is.finite(D$dina_proposed_n)]
    if (length(PR)) { q <- stats::quantile(PR, c(.25,.50,.75,.90), names = FALSE)
      cat(sprintf("  >> DINA_PROPOSED_N          : min %g  q25 %g  MED %g  q75 %g  p90 %g  max %g  (mean %.1f, n %d); searched %s\n",
          min(PR), q[1], q[2], q[3], q[4], max(PR), mean(PR), length(PR),
          paste(unique(range(D$dina_searched_n, na.rm = TRUE)), collapse = "-"))) }
    E2 <- c("dina_searched_n", "dina_proposed_n", "dina_tau_min", "admitted_n")
    P("DINA proposal fields filled on every declared replicate",
      all(E2 %in% names(r)) && all(vapply(E2, function(k) all(is.finite(D[[k]])), logical(1))),
      sprintf("(%s)", paste(sprintf("%s %d/%d", E2, vapply(E2, function(k) sum(is.finite(D[[k]])), 1L), nrow(D)), collapse = ", ")))
    # DINA's tau-hat is on the LINK scale: forestsearch() derives the floor as
    # m_diff = log(hr.threshold) for every non-Gaussian family
    # (R/forestsearch_helpers.R:1434-1437), and .dina_collect_candidates() drops
    # mean_tau < m_diff (R/dina_subgroup.R:748-749).  The floor AS APPLIED is the
    # OR-scale threshold on the harm side, asserted on the log scale.
    P("every proposed candidate at oriented tau-hat >= log(0.90), the OR effect threshold on the harm side",
      all(D$dina_tau_min >= log(0.90) - 1e-12),
      sprintf("(min %.8f = OR %.6f; floor log(0.90) = %.8f)",
              min(D$dina_tau_min), exp(min(D$dina_tau_min)), log(0.90)))
    P("1 <= admitted_n <= dina_proposed_n on every declared replicate",
      all(D$admitted_n >= 1L & D$admitted_n <= D$dina_proposed_n))
  }
  for (k in c("n_cons_qual","band_n"))
    cat(sprintf("  STRUCTURAL %-12s %s\n", k,
        if (all(is.na(r[[k]]))) sprintf("present, all-NA -- %s on %s, NOT a failure",
                                        if (identical(meth, "consistency")) "UNEXPECTED" else "STRUCTURAL (no consistency screen)", meth)
        else sprintf("POPULATED (expected on FS%s)", if (identical(meth, "consistency")) "" else " only -- unexpected here")))
  P("n_family finite on every declared replicate with a gate", all(is.finite(D$n_family[D$mr_ok %in% 1L])))
  P("p_hat_H and p_hat_sum recorded on declared replicates with a gate",
    all(is.finite(D$p_hat_H[D$mr_ok %in% 1L])) && all(is.finite(D$p_hat_sum[D$mr_ok %in% 1L])))
  P("p-hat validity (0<=p<=1, p_H<=sum)",
    all(D$p_hat_H >= 0 & D$p_hat_H <= 1, na.rm = TRUE) &&
    all(D$p_hat_H <= D$p_hat_sum + 1e-12, na.rm = TRUE))
  wm <- r$warn_msg[!is.na(r$warn_msg)]
  cat(sprintf("  >> WARNINGS                 : %d of %d rows carry warn_msg; distinct: %s\n", length(wm), nrow(r),
      if (length(wm)) paste(unique(gsub(" \\[x[0-9]+\\]", "", unlist(strsplit(wm, " | ", fixed = TRUE)))), collapse = " || ") else "none"))
  P("zero factor-comparison warnings", sum(grepl("not meaningful for factors", wm, fixed = TRUE)) == 0L)

  # --- §1.5(d) the constructions, on every declared replicate ---
  fin <- c("nv_H_est","nv_Hc_est","mr_H_est","mr_Hc_est","mr_H_lo","mr_H_hi","mr_Hc_lo","mr_Hc_hi",
           "fld_H_est2","fld_H_lo1s","fld_H_se","betaHhat_H","betaHhat_Hc","p_hat_H",
           "C_dagger_H","C_dagger_Hc","C_ddagger_H","C_ddagger_Hc")
  miss <- setdiff(fin, names(r))
  P("every finiteness column present (incl. C_dagger_* / C_ddagger_*)", length(miss) == 0L,
    if (length(miss)) paste("absent:", paste(miss, collapse=", ")) else "")
  Dm1 <- D[D$mr_ok %in% 1L & is.finite(D$fld_H_est2), , drop = FALSE]
  nf <- fin[vapply(fin, function(k) any(!is.finite(Dm1[[k]])), logical(1))]
  P("harm products finite on declared replicates with a field block", length(nf) == 0L,
    if (length(nf)) paste("non-finite:", paste(nf, collapse=", ")) else sprintf("(%d rows)", nrow(Dm1)))
  f <- D[is.finite(D$fld_Hc_est2), , drop = FALSE]
  sC <- c("fld_Hc_est2_s","fld_Hc_up1s_s","fld_Hc_lo1s_s","fld_Hc_lo2s_s","fld_Hc_hi2s_s",
          "fld_Hc_lo_se_s","fld_Hc_hi_se_s","fld_Hc_se_s","fld_Hc_lam_mean_s")
  jC <- c("fld_joint_s_gamma","fld_joint_s_prob","fld_joint_s_loH","fld_joint_s_upHc",
          "fld_joint_s_bonf_loH","fld_joint_s_bonf_upHc","fld_joint_s_bonf_prob",
          "fld_joint_s_corr","fld_joint_s_n")
  P("nine fld_Hc_*_s and nine fld_joint_s_* columns present", all(c(sC, jC) %in% names(r)))
  P(sprintf("fld_Hc_*_s finite on all %d filled replicates", nrow(f)),
    nrow(f) > 0 && all(vapply(sC, function(k) all(is.finite(f[[k]])), logical(1))))
  P(sprintf("fld_joint_s_* finite on all %d filled replicates", nrow(f)),
    nrow(f) > 0 && all(vapply(jC, function(k) all(is.finite(f[[k]])), logical(1))))
  cat(sprintf("  complement field filled on %d of %d declared replicates (%d notes)\n",
              nrow(f), nrow(D), sum(!is.na(D$fld_Hc_note))))

  P("invariant harm  : fld_H_lo1s <= fld_H_est2",       all(Dm1$fld_H_lo1s <= Dm1$fld_H_est2))
  P("invariant compl : fld_Hc_est2 <= fld_Hc_up1s",     all(f$fld_Hc_est2 <= f$fld_Hc_up1s))
  P("invariant _s    : fld_Hc_lo1s_s <= fld_Hc_up1s_s", all(f$fld_Hc_lo1s_s <= f$fld_Hc_up1s_s))
  P("invariant _s    : fld_Hc_lo2s_s <= fld_Hc_hi2s_s", all(f$fld_Hc_lo2s_s <= f$fld_Hc_hi2s_s))
  P("invariant _s    : fld_Hc_est2_s <= fld_Hc_up1s_s", all(f$fld_Hc_est2_s <= f$fld_Hc_up1s_s))
  P("invariant joint : bonf_loH <= fld_H_est2",         all(f$fld_joint_bonf_loH <= f$fld_H_est2))
  P("invariant jointS: fld_Hc_est2_s <= bonf_upHc_s",   all(f$fld_Hc_est2_s <= f$fld_joint_s_bonf_upHc))
  P("invariant IJ    : mr_H_lo <= est <= mr_H_hi",      all(Dm1$mr_H_lo <= Dm1$mr_H_est & Dm1$mr_H_est <= Dm1$mr_H_hi))
  # The OR scale, with degenerate logistic fits accounted for (§1.5(d)).  Under complete
  # separation in an arm the logistic MLE diverges (the study's .logit_or_ci() guards only the
  # OVERALL >= 5 events / >= 5 non-events), beta-hat and SE blow up together, and exp()
  # underflows the lower bound to 0 / overflows the upper to Inf.  The bound is still positive
  # mathematically; the double cannot hold it.  So: ESTIMATES must be strictly positive, a
  # bound is a failure only when NEGATIVE, and zero / infinite bounds are counted and reported
  # as a degenerate-fit diagnostic.  Such rows are already dropped from every coverage figure
  # by the finiteness masks, here and in the committed study.
  PEST <- intersect(c("or_H_est","or_Hc_est","nv_H_est","nv_Hc_est","mr_H_est","mr_Hc_est",
                      "fld_H_est2","fld_Hc_est2","fld_Hc_est2_s","betaHhat_H","betaHhat_Hc",
                      "C_dagger_H","C_dagger_Hc","C_ddagger_H","C_ddagger_Hc"), names(r))
  PBND <- intersect(c("or_H_lo","or_H_hi","or_Hc_lo","or_Hc_hi","nv_H_lo","nv_H_hi","nv_Hc_lo","nv_Hc_hi",
                      "mr_H_lo","mr_H_hi","mr_Hc_lo","mr_Hc_hi",
                      "fld_H_lo1s","fld_H_lo2s","fld_H_hi2s",
                      "fld_Hc_up1s","fld_Hc_lo1s","fld_Hc_lo2s","fld_Hc_hi2s",
                      # field-s BOUNDS only: fld_Hc_se_s and fld_Hc_lam_mean_s are log-OR
                      # quantities (R/fs_mr_inference.R:480-488), routinely negative, not bounds.
                      "fld_Hc_up1s_s","fld_Hc_lo1s_s","fld_Hc_lo2s_s","fld_Hc_hi2s_s","fld_Hc_lo_se_s","fld_Hc_hi_se_s",
                      "fld_joint_loH","fld_joint_upHc","fld_joint_bonf_loH","fld_joint_bonf_upHc",
                      "fld_joint_s_loH","fld_joint_s_upHc","fld_joint_s_bonf_loH","fld_joint_s_bonf_upHc"), names(r))
  ne <- vapply(PEST, function(k) { x <- r[[k]]; sum(is.finite(x) & x <= 0) }, integer(1))
  P(sprintf("every ESTIMATE is a positive OR (%d columns)", length(PEST)), all(ne == 0L),
    if (any(ne > 0L)) paste(sprintf("%s %d", names(ne)[ne > 0], ne[ne > 0]), collapse = ", ") else "")
  nb <- vapply(PBND, function(k) { x <- r[[k]]; sum(is.finite(x) & x < 0) }, integer(1))
  P(sprintf("no NEGATIVE bound (%d columns)", length(PBND)), all(nb == 0L),
    if (any(nb > 0L)) paste(sprintf("%s %d", names(nb)[nb > 0], nb[nb > 0]), collapse = ", ") else "")
  dg <- vapply(PBND, function(k) { x <- r[[k]]; sum((is.finite(x) & x == 0) | is.infinite(x)) }, integer(1))
  dor <- sum(is.finite(r$or_H_est) & (!is.finite(r$or_H_hi) | r$or_H_lo %in% 0), na.rm = TRUE)
  dorc <- sum(is.finite(r$or_Hc_est) & (!is.finite(r$or_Hc_hi) | r$or_Hc_lo %in% 0), na.rm = TRUE)
  cat(sprintf("  >> DEGENERATE BOUNDS (separation) : %s | oracle rows affected: H %d, Hc %d of %d\n",
      if (any(dg > 0L)) paste(sprintf("%s %d", names(dg)[dg > 0], dg[dg > 0]), collapse = ", ") else "none",
      dor, dorc, nrow(r)))

  g1 <- f$fld_joint_gamma; g2 <- f$fld_joint_s_gamma
  P("gamma (joint)   in [0.025, 0.05]", all(g1 >= 0.025 - 1e-12 & g1 <= 0.05 + 1e-12),
    sprintf("[%.5f, %.5f]", min(g1), max(g1)))
  P("gamma (joint-s) in [0.025, 0.05]", all(g2 >= 0.025 - 1e-12 & g2 <= 0.05 + 1e-12),
    sprintf("[%.5f, %.5f]", min(g2), max(g2)))
  # est2 = to_eff(beta_deb - lambda_mean) with to_eff = exp on a ratio measure
  # (R/fs_mr_inference.R:931, :480-488): lambda_mean and the Lambda* quantiles are
  # on the WORKING (log-OR) scale, est2 and the bounds on the EFFECT (OR) scale.
  # The MD checker's additive identities are identity-scale specializations; here
  # they are written on the log scale.
  i1 <- max(abs((log(f$fld_Hc_est2_s) + f$fld_Hc_lam_mean_s) -
                (log(f$fld_Hc_est2)   + f$fld_Hc_lam_mean)))
  P("identity: field-s inverted around the same beta-tilde^c (log scale)", i1 <= 1e-9, sprintf("max |diff| = %.3g", i1))
  i1b <- max(abs(log(f$fld_Hc_est2) + f$fld_Hc_lam_mean - log(f$mr_Hc_est)),
             abs(log(f$fld_Hc_est2_s) + f$fld_Hc_lam_mean_s - log(f$mr_Hc_est)))
  P("identity: log(est2) + lambda_mean = log(beta-tilde^c), field and field-s", i1b <= 1e-9, sprintf("max |diff| = %.3g", i1b))
  agree <- f$fld_joint_n == f$fld_joint_s_n
  dj <- if (any(agree)) max(abs(f$fld_joint_bonf_loH[agree] - f$fld_joint_s_bonf_loH[agree])) else 0
  P("identity: Bonferroni harm bound joint == joint_s where draw counts agree", dj <= 1e-12,
    sprintf("(%d of %d agree; max |diff| %.3g)", sum(agree), nrow(f), dj))
  i2 <- max(abs(log(Dm1$fld_H_lo1s) - (log(Dm1$mr_H_est) - Dm1$fld_H_q95)),
            abs(log(f$fld_Hc_up1s)   - (log(f$mr_Hc_est)  - f$fld_Hc_q05)))
  P("identity: log(lo1s) = log(beta-tilde) - q95; log(up1s) = log(beta-tilde^c) - q05", i2 <= 1e-9, sprintf("max |diff| = %.3g", i2))
  # The per-replicate target columns are the cell's population constants.
  P("C_dagger_* / C_ddagger_* equal the truth table on every row",
    max(abs(r$C_dagger_H - tr$marg_H), abs(r$C_dagger_Hc - tr$marg_Hc),
        abs(r$C_ddagger_H - tr$cde_H), abs(r$C_ddagger_Hc - tr$cde_Hc)) <= 1e-12)
  cat(sprintf("  >> p-hat(H): mean %.3f, share < 0.5: %.3f | CLASSIFICATION: sens %.4f ppv %.4f | mean |Hhat| %.1f\n",
      mean(Dm1$p_hat_H), mean(Dm1$p_hat_H < 0.5), mean(D$sens, na.rm = TRUE),
      mean(D$ppv, na.rm = TRUE), mean(D$n_harm, na.rm = TRUE)))
  cat(sprintf("  >> BOUNDS: mean fld_H_lo1s %.4f | share >= 1.0 %.3f | mean fld_Hc_up1s_s %.4f | share <= 1.0 %.3f\n",
      mean(Dm1$fld_H_lo1s), mean(Dm1$fld_H_lo1s >= 1.0),
      mean(f$fld_Hc_up1s_s), mean(f$fld_Hc_up1s_s <= 1.0)))

  # --- same draws across identifiers, both directions (orgrf / ordina vs orfs) ---
  if (!identical(tag, "orfs")) {
    fp <- fcomb(target, n, "orfs")
    P("orfs comparator on disk", file.exists(fp), basename(fp))
    if (file.exists(fp)) {
      om <- readRDS(fp); fr <- om$results
      same_draws(r, fr, tr, om$truth, sprintf("%s -> orfs", tag))
      same_draws(fr, r, om$truth, tr, sprintf("orfs -> %s", tag))
      cat(sprintf("  orfs comparator (same cell): declaration %.4f | mean |Hhat| %.1f | mean fld_H_lo1s %.4f\n",
                  mean(fr$detected %in% 1L), mean(fr$n_harm, na.rm = TRUE),
                  mean(fr$fld_H_lo1s, na.rm = TRUE)))
    }
  } else cat("  same-draws: this IS the orfs campaign (the reference for the other two).\n")

  for (fp2 in c(f1, f2, cf)) if (file.exists(fp2)) {
    sz <- file.size(fp2)
    P(sprintf("size <= 100 MB: %s", basename(fp2)), sz <= 100 * 1024^2,
      sprintf("(%d B)%s", sz, if (sz > 50 * 1024^2) "  ** FLAG: over 50 MB **" else ""))
  }
  cat(sprintf("  timing: fit_mr_secs mean %.1f median %.1f p90 %.1f max %.1f | id_secs mean %.2f median %.2f max %.2f | fld_H_secs mean %.1f | fld_Hc_secs mean %.2f\n",
              mean(r$fit_mr_secs), stats::median(r$fit_mr_secs),
              stats::quantile(r$fit_mr_secs, .9, names = FALSE), max(r$fit_mr_secs),
              mean(r$id_secs, na.rm = TRUE), stats::median(r$id_secs, na.rm = TRUE), max(r$id_secs, na.rm = TRUE),
              mean(r$fld_H_secs, na.rm = TRUE), mean(r$fld_Hc_secs, na.rm = TRUE)))
  invisible(NULL)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) stop("usage: Rscript gate2.R <target> <n> <cell> <tag>")
res <- tryCatch({ gate2(args[1], as.integer(args[2]), args[3], args[4]); TRUE },
                error = function(e) { P("checker ran without error", FALSE, conditionMessage(e)); FALSE })
cat(sprintf("\nGATE_COUNTS run=%d passed=%d failed=%d\n", N_RUN, N_PASS, N_RUN - N_PASS))
if (length(FAILED)) cat("FAILED:", paste(FAILED, collapse = " || "), "\n")
quit(status = if (N_RUN > N_PASS) 1L else 0L)
