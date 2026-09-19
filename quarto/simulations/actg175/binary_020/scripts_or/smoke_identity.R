# Smoke checks for the ACTG175 binary/OR campaigns (TASK_actg175_binary_campaign_2026-09-17 §1.5).
# Transplant of ../../continuous/scripts_mddina/smoke_identity.R with the OR fields and meta, in
# four modes matching the task's four renders:
#   recipe : §1.5(a) THE DATA RECIPE -- the knobs at the STUDY's rule (maxeffCons, eps 0.10),
#            compared against the committed study bundle mr_sweep/<run_tag>/fs_mr_n<n>_res.rds on
#            sim_id 1..nsims: the truths, the oracle estimates and every other RULE-INDEPENDENT
#            data-level column (seed, n_true) within 1e-8 relative.  A difference means the port
#            is unfaithful and the checker exits 1.  The selection agreement under the same rule
#            is REPORTED beside, not gated (the MR construction differs: field vs ij).
#   fs     : §1.5(b) FS under the CAMPAIGN rule (effMaxSG, eps 0.20) -- it runs, declares, and its
#            meta carries the rule, the thresholds, target_or_h, the truths,
#            field_scale_complement = "selected", pkg_version and host.
#   grf    : §1.5(c) GRF under the campaign rule -- zero factor-comparison warnings, zero
#            NA-membership candidates.
#   dina   : §1.5(c) DINA under the campaign rule -- the same, plus: the proposal and admission
#            floors, AS APPLIED, are the OR-scale threshold on the harm side (dina_tau_min >= the
#            oriented effect threshold), and 1 <= admitted_n <= dina_proposed_n.
# §1.5(d) THE CONSTRUCTIONS runs in EVERY mode on every declared replicate: the nine fld_Hc_*_s
# and nine fld_joint_s_* columns finite, fld_Hc_lo1s_s <= fld_Hc_up1s_s, the Bonferroni harm bound
# identical between joint and joint_s where the draw counts agree, every bound a positive OR, and
# p-hat(Hhat) recorded.
# Facts reported in every mode: the declared count, each identifier's selection on sim_id 1 beside
# S0 §7's fit, fit_mr_secs and the field seconds.
# EVERY coverage figure on the GRF and DINA paths is coverage of beta(H-hat) CONDITIONAL ON THE
# PROPOSED FAMILY (their families are generated from fitted surfaces).  FS's family is the
# prespecified cut grid.
# usage: Rscript smoke_identity.R <target: 0.75|1.0|1.5> <n> <tag> <nsims> <mode: recipe|fs|grf|dina>
suppressPackageStartupMessages(library(forestsearch))
args <- commandArgs(TRUE)
if (length(args) != 5L) stop("usage: Rscript smoke_identity.R <target> <n> <tag> <nsims> <mode>")
target <- args[1]; n <- as.integer(args[2]); tag <- args[3]
nsims <- as.integer(args[4]); mode <- args[5]
stopifnot(mode %in% c("recipe", "fs", "grf", "dina"))
setwd(Sys.getenv("ORSG_DIR", unset = "~/Documents/GitHub/forestsearch/quarto/simulations/actg175/binary_020"))
HOST  <- Sys.getenv("ORSG_HOST", unset = "pop-os")
PKG   <- Sys.getenv("ORSG_PKG",  unset = "0.3.5")
TOL   <- 1e-8
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a

design_tag <- sprintf("or%03d", as.integer(round(100 * as.numeric(target))))
# Under the STUDY's rule eps is 0.10, which the template leaves untagged; the campaign rule tags _nb20.
nbhd_tag <- if (identical(mode, "recipe")) "" else "_nb20"
focus_tg <- if (identical(mode, "recipe")) "maxeffCons" else "effMaxSG"  # fs_focus_tag(): "maxeffCons" is unmapped on the consistency path, so it passes through
meth     <- switch(mode, recipe = "fs", fs = "fs", grf = "grf", dina = "dina")
method_f <- switch(mode, recipe = "consistency", fs = "consistency", grf = "grf", dina = "dina")
stem <- sprintf("%s_%s_mr_field_%s_n%d%s_%s", meth, focus_tg, design_tag, n, nbhd_tag, tag)
bp   <- file.path("mr_or_harm", paste0(stem, "_d5000"), sprintf("%s_res_1_%d.rds", stem, nsims))
# The committed study bundle for this (identifier, n): the (a) comparator, and the reference
# selection quoted beside every mode's sim_id 1 fact.
RUNTAG <- "maxeffCons_actg175_or075_seedtab_s1000"
op <- file.path("mr_sweep", RUNTAG, sprintf("%s_mr_n%d_res.rds", meth, n))

cat(sprintf("== SMOKE %s: target_or_h=%s n=%d tag=%s | %s vs committed %s, sim_id 1-%d ==\n",
            mode, target, n, tag, basename(bp), basename(op), nsims))
ok_all <- TRUE
chk <- function(cond, msg) { cat(sprintf("  [%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg))
                             ok_all <<- ok_all && isTRUE(cond) }
relmax <- function(x, y) { d <- abs(x - y) / pmax(abs(y), 1e-300)
  d[is.na(x) & is.na(y)] <- 0; d[xor(is.na(x), is.na(y))] <- Inf; max(d) }

chk(file.exists(bp), paste("smoke bundle exists:", bp))
chk(file.exists(op), paste("committed study bundle exists:", op))
if (!file.exists(bp)) { cat(sprintf("SMOKE %s: FAIL (no bundle)\n", mode)); quit(status = 1L) }
b <- readRDS(bp); r1 <- b$results; m <- b$meta; tr <- b$truth
o <- if (file.exists(op)) readRDS(op) else NULL
ro <- if (!is.null(o)) o$results[o$results$sim_id %in% seq_len(nsims), ] else NULL
r1 <- r1[order(r1$sim_id), ]; if (!is.null(ro)) ro <- ro[order(ro$sim_id), ]

chk(identical(as.integer(r1$sim_id), seq_len(nsims)), sprintf("carries sim_id 1-%d (%d rows)", nsims, nrow(r1)))
chk(all(r1$status %in% c("DETECTED", "NO-DETECTION")),
    sprintf("no CONFIG-ERROR rows (status: %s)",
            paste(names(table(r1$status)), table(r1$status), collapse = ", ")))
if (any(!is.na(r1$err_msg)))
  cat("  err_msg present on", sum(!is.na(r1$err_msg)), "row(s):",
      paste(unique(stats::na.omit(r1$err_msg)), collapse = " || "), "\n")

# ---- truths (every mode): the template's DGM must reproduce the committed truth table ----
if (!is.null(o)) {
  tk <- c("or_causal", "marg_H", "marg_Hc", "cde_H", "cde_Hc")
  # The committed study planted the region at prevalence 9.632% (sg_quantile 0.70).
  # TASK_binary_study_redesign_2026-09-18 raised it to 14.917% (sg_quantile 0.62850),
  # so under the design of record the committed truths describe a SUPERSEDED
  # super-population and are not a comparator for anything.  The gate therefore
  # applies only when the two carry the SAME planted prevalence; otherwise the
  # comparison is REPORTED with both prevalences named, never gated.
  # The committed study's truth list carries no prevalence, so the planted
  # region is compared through meta$sg_quantile, which both bundles carry.
  .sgq_here  <- m$sg_quantile %||% NA_real_
  .sgq_there <- o$meta$sg_quantile %||% NA_real_
  same_prev  <- isTRUE(abs(.sgq_here - .sgq_there) <= 1e-10)
  if (identical(design_tag, "or075") && same_prev) {
    tmx <- max(vapply(tk, function(k) relmax(tr[[k]], o$truth[[k]]), numeric(1)))
    chk(tmx <= TOL, sprintf("truth targets match the committed study within %g relative (max %.3g); %s",
                            TOL, tmx, paste(sprintf("%s %.10f", tk, unlist(tr[tk])), collapse = " | ")))
  } else if (identical(design_tag, "or075")) {
    cat(sprintf("  SUPERSEDED BY DESIGN CHANGE: the committed study plants at sg_quantile %s; this DGM of record plants at %s, prevalence(H) %.6f. The committed truths are not a comparator for it, and are not gated. Truths here: %s\n",
                format(.sgq_there), format(.sgq_here), tr$prevalence_Q %||% NA_real_,
                paste(sprintf("%s %.10f", tk, unlist(tr[tk])), collapse = " | ")))
  } else {
    cat(sprintf("  truths at %s (no committed comparator at this design point): %s\n",
                design_tag, paste(sprintf("%s %.10f", tk, unlist(tr[tk])), collapse = " | ")))
  }
  cat(sprintf("  prevalence(H) %.6f | k_inter %.7f\n", tr$prevalence_Q, tr$beta_inter))
}

# ---- (a) the data recipe: RULE-INDEPENDENT data-level columns ----------------------------
# Strictly rule-independent: the per-replicate seed, the true-region size, and the oracle
# quadruples (refit on the TRUE region, so the realized rule cannot touch them).  The template
# records the oracle under the MD naming or_*; the study uses ora_*.
if (identical(mode, "recipe")) {
  if (is.null(ro)) { chk(FALSE, "committed comparator readable"); }
  else {
    chk(nrow(ro) == nsims && identical(as.integer(ro$sim_id), seq_len(nsims)),
        sprintf("committed comparator carries sim_id 1-%d", nsims))
    chk(identical(as.integer(r1$seed), as.integer(ro$seed)),
        sprintf("seed identical on every row (sim_id 1 seed %d, committed %d)", r1$seed[1], ro$seed[1]))
    chk(identical(as.integer(r1$n_true), as.integer(ro$n_true)),
        sprintf("n_true identical on every row (sim_id 1: %d vs %d)", r1$n_true[1], ro$n_true[1]))
    pairs <- rbind(
      data.frame(here = paste0("or_H_",  c("est","lo","hi","se")),
                 there = paste0("ora_H_", c("est","lo","hi","se")), stringsAsFactors = FALSE),
      data.frame(here = paste0("or_Hc_",  c("est","lo","hi","se")),
                 there = paste0("ora_Hc_", c("est","lo","hi","se")), stringsAsFactors = FALSE))
    # The committed study computes the oracle only AFTER its detection return
    # (maxeffCons_mr_coverage_sweep_or075.qmd :628-631), so its ora_* columns are
    # NA on NO-DETECTION rows.  This template follows the MD template and computes
    # the oracle BEFORE the search, so it is filled on every replicate.  The
    # comparison therefore runs on the rows the committed bundle HAS an oracle on;
    # the extra rows are reported, and are strictly more information on a
    # rule-independent quantity, not a disagreement about a value.
    kk <- is.finite(as.numeric(ro$ora_H_est))
    omx <- max(vapply(seq_len(nrow(pairs)), function(i)
      relmax(as.numeric(r1[[pairs$here[i]]])[kk], as.numeric(ro[[pairs$there[i]]])[kk]), numeric(1)))
    chk(omx <= TOL, sprintf("oracle columns (or_* here vs ora_* committed, 8 columns) within %g relative on all %d rows the committed bundle scores (max %.3g)", TOL, sum(kk), omx))
    extra <- which(!kk & is.finite(r1$or_H_est))
    cat(sprintf("  REPORTED: the committed bundle leaves ora_* NA on %d row(s) (%s -- all NO-DETECTION, its oracle runs after the detection return); this template fills the oracle on all %d (rule-independent, computed before the search)\n",
                length(extra), paste(r1$sim_id[extra], collapse = ","), nrow(r1)))
    cat(sprintf("  FACT sim_id 1 oracle: or_H_est %.10f (committed ora_H_est %.10f) | or_Hc_est %.10f (committed %.10f)\n",
                r1$or_H_est[1], ro$ora_H_est[1], r1$or_Hc_est[1], ro$ora_Hc_est[1]))
    # REPORTED, not gated: the selection under the same rule, and the naive estimate.
    same_sel <- mapply(function(a, z) identical(sort(strsplit(a %||% "", " & ")[[1]]),
                                                sort(strsplit(z %||% "", " & ")[[1]])),
                       r1$sg_def, ro$sg_def)
    nvd <- relmax(as.numeric(r1$nv_H_est), as.numeric(ro$nv_H_est))
    cat(sprintf("  REPORTED (rule-DEPENDENT, not a gate): same realized rule on %d of %d rows; status agrees on %d; naive OR max rel diff %.3g\n",
                sum(same_sel %in% TRUE), nsims,
                sum((r1$status == ro$status) %in% TRUE), nvd))
    dif <- which(!(same_sel %in% TRUE))
    for (i in dif)
      cat(sprintf("    sim %d: committed [%s] | here [%s]\n", r1$sim_id[i],
                  ro$sg_def[i] %||% "<none>", r1$sg_def[i] %||% "<none>"))
  }
  chk(identical(m$sg_focus, "maxeffCons") && isTRUE(all.equal(m$effect_neighborhood, 0.10)),
      sprintf("meta carries the STUDY's rule (sg_focus %s, eps %s)", m$sg_focus, format(m$effect_neighborhood)))
}

# ---- (b)/(c) meta under the campaign rule -------------------------------------------------
if (!identical(mode, "recipe")) {
  cat(sprintf("  meta: method %s | focus %s/%s | eps %s | rule %s | ci %s | scalec %s | draws %s | thresholds %s/%s/%s | adverse %s | target_or_h %s | pkg %s | host %s | workers %s | R %s\n",
              m$subgroup_method, m$sg_focus, m$focus_tag, format(m$effect_neighborhood),
              m$selection_rule, m$ci_method, m$field_scale_complement, format(m$mr_draws),
              format(m$effect_threshold), format(m$consistency_threshold), format(m$pconsistency),
              m$adverse_outcome, format(m$target_or_h), m$pkg_version, m$hostname,
              format(m$n_workers), m$r_version))
  chk(identical(m$subgroup_method, method_f) &&
      identical(m$sg_focus, "effMaxSG") && isTRUE(all.equal(m$effect_neighborhood, 0.20)) &&
      identical(m$selection_rule, "neighborhood"),
      sprintf("meta carries the CAMPAIGN rule (identifier %s, effMaxSG, eps 0.20, neighborhood)", method_f))
  chk(isTRUE(all.equal(m$effect_threshold, 0.90)) &&
      isTRUE(all.equal(m$consistency_threshold, 0.80)) &&
      isTRUE(all.equal(m$pconsistency, 0.90)) && isTRUE(m$adverse_outcome),
      "meta carries the study's thresholds 0.90 / 0.80 / 0.90 with adverse_outcome = TRUE")
  chk(isTRUE(all.equal(m$target_or_h, as.numeric(target))) &&
      identical(m$design_tag, design_tag) && identical(m$dgm_model, "alt"),
      sprintf("meta carries target_or_h %s and design_tag %s on the calibrated alt branch", target, design_tag))
  chk(all(vapply(c("truth_marg_H","truth_marg_Hc","truth_cde_H","truth_cde_Hc",
                   "truth_or_causal","truth_prevalence_Q"),
                 function(k) is.finite(m[[k]] %||% NA_real_), logical(1))),
      "meta carries the truths (marg_H/Hc, cde_H/Hc, or_causal, prevalence)")
  chk(identical(m$field_scale_complement, "selected") && identical(m$ci_method, "field") &&
      isTRUE(m$field_complement) && identical(m$ij_residual, "two_term") &&
      isTRUE(m$return_reselection) && identical(as.integer(m$mr_draws), 5000L),
      "meta carries the MR construction set (field, complement, field-s selected, two_term, reselection, 5000 draws)")
  chk(identical(m$pkg_version, PKG) && identical(m$hostname, HOST),
      sprintf("meta carries pkg_version %s and host %s", PKG, HOST))
  # Identifier arguments, as applied.
  if (identical(mode, "grf"))
    chk(isTRUE(all.equal(m$dmin_grf, 0.0)) && identical(m$grf_selection, "frontier") &&
        identical(as.integer(m$grf_depth), 2L) && identical(m$grf_select_statistic, "effect"),
        sprintf("meta carries the study's GRF arguments (dmin %s, %s, depth %s, %s)",
                format(m$dmin_grf), m$grf_selection, format(m$grf_depth), m$grf_select_statistic))
  if (identical(mode, "dina"))
    chk(identical(m$dina_select_statistic, "effect") && identical(m$dina_args, "list()"),
        sprintf("meta carries the study's DINA arguments (%s, %s)", m$dina_select_statistic, m$dina_args))
}

D <- r1[r1$detected %in% 1L, , drop = FALSE]
cat(sprintf("  FACT declared: %d of %d | NO-DETECTION %d | CONFIG-ERROR %d\n",
            nrow(D), nsims, sum(r1$detected %in% 0L), sum(is.na(r1$detected))))
if (!identical(mode, "recipe")) chk(nrow(D) > 0L, sprintf("%s declares on at least one replicate", method_f))

# ---- warnings: zero factor-comparison warnings; zero NA-membership candidates (§1.5(c)) ----
wm <- r1$warn_msg[!is.na(r1$warn_msg)]
cat(sprintf("  WARNINGS (distinct, with row counts, on %d of %d rows):\n", length(wm), nsims))
for (w in unique(wm)) cat("    ", w, "  [rows:", sum(wm == w), "]\n")
if (!length(wm)) cat("     none\n")
nfac <- sum(grepl("not meaningful for factors", wm, fixed = TRUE))
nna  <- sum(grepl("NA|missing", wm) & grepl("member|subgroup|candidate", wm))
chk(nfac == 0L, sprintf("zero factor-comparison warnings (%d rows carry one)", nfac))
chk(nna == 0L, sprintf("zero NA-membership candidate warnings (%d rows carry one)", nna))

# ---- sim_id 1 selection, beside S0 §7's fit ------------------------------------------------
S0FIT <- list(
  `consistency-500`  = "{race} & !{karnof <= 95} (n 73), oriented log-OR 0.8303, MR family 2238, field lower_1s 0.2937, field-s upper_1s_s 1.1332",
  `grf-500`          = "{age <= 28} & {preanti <= 777.2} (n 89), oriented log-OR 0.6484, MR family 1051, field lower_1s 0.2920, field-s upper_1s_s 1.1375",
  `dina-500`         = "{age <= 28} & {preanti <= 777.2} (n 89), oriented log-OR 0.6484, MR family 1157, field lower_1s 0.3028, field-s upper_1s_s 1.1313",
  `consistency-2000` = "!{preanti <= 741} & !{wtkg <= 78} (n 145), oriented log-OR 0.3844, MR family 2978, field lower_1s 0.2824, field-s upper_1s_s 0.7565",
  `grf-2000`         = "{preanti > 0} & {wtkg <= 64.5} (n 276), oriented log-OR 0.1354, MR family 1345, field lower_1s 0.2127, field-s upper_1s_s 0.7459",
  `dina-2000`        = "{wtkg >= 68.04} & {cd40 <= 214} (n 129), oriented log-OR 0.3589, MR family 70, field lower_1s 0.3015, field-s upper_1s_s 0.7501")
s1 <- r1[r1$sim_id == 1L, ]
cat(sprintf("  FACT sim_id 1 [%s]: status %s | sg_def [%s] | n_sel %s | n_harm %s | naive OR %s (oriented log-OR %s) | MR family %s | fld_H_lo1s %s | fld_Hc_up1s_s %s | p_hat_H %s | sens %s ppv %s\n",
            method_f, s1$status, s1$sg_def %||% "<none>", format(s1$n_sel), format(s1$n_harm),
            format(round(s1$nv_H_est, 6)),
            if (is.finite(s1$nv_H_est) && s1$nv_H_est > 0) format(round(log(s1$nv_H_est), 4)) else "NA",
            format(s1$n_family), format(round(s1$fld_H_lo1s, 4)),
            format(round(s1$fld_Hc_up1s_s, 4)), format(round(s1$p_hat_H, 3)),
            format(round(s1$sens, 3)), format(round(s1$ppv, 3))))
k0 <- sprintf("%s-%d", method_f, n)
if (!is.null(S0FIT[[k0]]))
  cat(sprintf("  FACT S0 section 7 fit (same design, effMaxSG eps 0.20, field, one worker), quoted: %s\n", S0FIT[[k0]]))
if (!is.null(ro)) {
  os1 <- ro[ro$sim_id == 1L, ]
  cat(sprintf("  FACT committed study payload (maxeffCons, eps 0.10, ci_method ij), sim_id 1: status %s | sg_def [%s] | n_sel %s | naive OR %s | t2_secs %s\n",
              os1$status, os1$sg_def %||% "<none>", format(os1$n_sel),
              format(round(os1$nv_H_est, 6)), format(round(os1$t2_secs, 1))))
}

# ---- DINA floors, as applied (§1.5(c)) -----------------------------------------------------
if (identical(mode, "dina") && nrow(D)) {
  E2 <- c("dina_searched_n", "dina_proposed_n", "dina_tau_min", "admitted_n")
  chk(all(vapply(E2, function(k) all(is.finite(D[[k]])), logical(1))),
      sprintf("DINA proposal fields filled on every declared replicate (%s)",
              paste(sprintf("%s %d/%d", E2, vapply(E2, function(k) sum(is.finite(D[[k]])), 1L), nrow(D)), collapse = ", ")))
  # DINA's tau-hat is on the LINK scale: forestsearch() derives the proposal floor
  # as m_diff = log(hr.threshold) for every non-Gaussian family
  # (R/forestsearch_helpers.R:1434-1437) and .dina_collect_candidates() drops any
  # candidate with mean_tau < m_diff (R/dina_subgroup.R:748-749).  So the floor AS
  # APPLIED is the OR-scale effect threshold on the harm side, asserted here on the
  # log scale and reported back as an OR.
  FLOOR <- as.numeric(m$effect_threshold %||% 0.90)
  chk(all(D$dina_tau_min >= log(FLOOR) - 1e-12),
      sprintf("every proposed candidate at oriented tau-hat >= log(effect threshold %.2f) = %.8f, i.e. the OR-scale threshold on the harm side (min over replicates %.8f = OR %.6f)",
              FLOOR, log(FLOOR), min(D$dina_tau_min), exp(min(D$dina_tau_min))))
  chk(all(D$admitted_n >= 1L & D$admitted_n <= D$dina_proposed_n),
      "1 <= admitted_n <= dina_proposed_n on every declared replicate")
  cat(sprintf("  FACT DINA family sizes (declared): searched %s | proposed min %d median %.0f max %d | admitted min %d median %.0f max %d | tau_min min %.6f\n",
              paste(unique(range(D$dina_searched_n)), collapse = "-"),
              min(D$dina_proposed_n), stats::median(D$dina_proposed_n), max(D$dina_proposed_n),
              min(D$admitted_n), stats::median(D$admitted_n), max(D$admitted_n), min(D$dina_tau_min)))
  cat(sprintf("  FACT DINA floor as applied: m_diff = log(%.2f) = %.8f (link scale); smallest proposed tau-hat %.8f = OR %.6f\n",
              FLOOR, log(FLOOR), min(D$dina_tau_min), exp(min(D$dina_tau_min))))
}
if (identical(mode, "grf") && nrow(D)) {
  A <- r1$admitted_n[is.finite(r1$admitted_n)]
  cat(sprintf("  FACT GRF admitted_n: %s (== 0 on %d of %d rows)\n",
              if (length(A)) sprintf("min %g median %g max %g (n %d)", min(A), stats::median(A), max(A), length(A)) else "none finite",
              sum(A == 0L), length(A)))
}

# ---- §1.5(d) THE CONSTRUCTIONS, every mode, every declared replicate -----------------------
sC <- c("fld_Hc_est2_s","fld_Hc_up1s_s","fld_Hc_lo1s_s","fld_Hc_lo2s_s","fld_Hc_hi2s_s",
        "fld_Hc_lo_se_s","fld_Hc_hi_se_s","fld_Hc_se_s","fld_Hc_lam_mean_s")
jC <- c("fld_joint_s_gamma","fld_joint_s_prob","fld_joint_s_loH","fld_joint_s_upHc",
        "fld_joint_s_bonf_loH","fld_joint_s_bonf_upHc","fld_joint_s_bonf_prob",
        "fld_joint_s_corr","fld_joint_s_n")
chk(all(c(sC, jC) %in% names(r1)), "nine fld_Hc_*_s and nine fld_joint_s_* columns present")
f <- D[is.finite(D$fld_Hc_est2), , drop = FALSE]
chk(nrow(f) > 0 && all(vapply(sC, function(k) all(is.finite(f[[k]])), logical(1))),
    sprintf("nine fld_Hc_*_s finite on all %d declared replicates with a complement field block", nrow(f)))
chk(nrow(f) > 0 && all(vapply(jC, function(k) all(is.finite(f[[k]])), logical(1))),
    sprintf("nine fld_joint_s_* finite on all %d such replicates", nrow(f)))
cat(sprintf("  complement field filled on %d of %d declared replicates (%d carry a degenerate note)\n",
            nrow(f), nrow(D), sum(!is.na(D$fld_Hc_note))))
chk(nrow(f) > 0 && all(f$fld_Hc_lo1s_s <= f$fld_Hc_up1s_s), "fld_Hc_lo1s_s <= fld_Hc_up1s_s on every such replicate")
chk(nrow(f) > 0 && all(f$fld_Hc_lo2s_s <= f$fld_Hc_hi2s_s), "fld_Hc_lo2s_s <= fld_Hc_hi2s_s on every such replicate")
agree <- f$fld_joint_n == f$fld_joint_s_n
dj <- if (any(agree)) max(abs(f$fld_joint_bonf_loH[agree] - f$fld_joint_s_bonf_loH[agree])) else 0
chk(dj <= 1e-12, sprintf("Bonferroni harm bound identical between joint and joint_s where the draw counts agree (%d of %d rows agree; max |diff| %.2e)",
                         sum(agree), nrow(f), dj))
# est2 = to_eff(beta_deb - lambda_mean) with to_eff = exp on a ratio measure
# (R/fs_mr_inference.R:931, :480-488), so the "same beta-tilde^c" identity is
# log(est2) + lambda_mean = log(beta-tilde^c), checked here on the log scale for
# BOTH the unstudentized field and field-s.  (The MD template's additive form is
# an identity-scale specialization and does not hold for an OR.)
idn <- if (nrow(f)) max(abs((log(f$fld_Hc_est2_s) + f$fld_Hc_lam_mean_s) -
                            (log(f$fld_Hc_est2)   + f$fld_Hc_lam_mean))) else NA_real_
idb <- if (nrow(f)) max(abs(log(f$fld_Hc_est2) + f$fld_Hc_lam_mean - log(f$mr_Hc_est))) else NA_real_
ids <- if (nrow(f)) max(abs(log(f$fld_Hc_est2_s) + f$fld_Hc_lam_mean_s - log(f$mr_Hc_est))) else NA_real_
chk(is.finite(idn) && idn <= 1e-9, sprintf("field-s inverted around the same beta-tilde^c, on the log scale (max |diff| %.2e)", idn))
chk(is.finite(idb) && idb <= 1e-9, sprintf("identity: log(est2) + lambda_mean = log(beta-tilde^c) (max |diff| %.2e)", idb))
chk(is.finite(ids) && ids <= 1e-9, sprintf("identity: log(est2_s) + lambda_mean_s = log(beta-tilde^c) (max |diff| %.2e)", ids))
# The OR scale, with degenerate logistic fits accounted for.  ESTIMATES must be strictly
# positive; a bound fails only when NEGATIVE.  Under complete separation in an arm the logistic
# MLE diverges and exp() underflows the lower bound to 0 / overflows the upper to Inf -- the
# bound is still positive mathematically, so those are counted and reported, not failed.
PEST <- intersect(c("or_H_est","or_Hc_est","nv_H_est","nv_Hc_est","mr_H_est","mr_Hc_est",
                    "fld_H_est2","fld_Hc_est2","fld_Hc_est2_s","betaHhat_H","betaHhat_Hc"), names(r1))
PBND <- intersect(c("or_H_lo","or_H_hi","or_Hc_lo","or_Hc_hi","nv_H_lo","nv_H_hi","nv_Hc_lo","nv_Hc_hi",
                    "mr_H_lo","mr_H_hi","mr_Hc_lo","mr_Hc_hi",
                    "fld_H_lo1s","fld_H_lo2s","fld_H_hi2s",
                    "fld_Hc_up1s","fld_Hc_lo1s","fld_Hc_lo2s","fld_Hc_hi2s",
                    # field-s BOUNDS only: fld_Hc_se_s and fld_Hc_lam_mean_s are log-OR
                    # quantities (R/fs_mr_inference.R:480-488), routinely negative, not bounds.
                    "fld_Hc_up1s_s","fld_Hc_lo1s_s","fld_Hc_lo2s_s","fld_Hc_hi2s_s","fld_Hc_lo_se_s","fld_Hc_hi_se_s",
                    "fld_joint_loH","fld_joint_upHc","fld_joint_bonf_loH","fld_joint_bonf_upHc",
                    "fld_joint_s_loH","fld_joint_s_upHc","fld_joint_s_bonf_loH","fld_joint_s_bonf_upHc"), names(r1))
ne <- vapply(PEST, function(k) { x <- r1[[k]]; sum(is.finite(x) & x <= 0) }, integer(1))
chk(all(ne == 0L), sprintf("every ESTIMATE is a positive OR (%d columns%s)", length(PEST),
    if (any(ne > 0L)) paste0("; ", paste(sprintf("%s %d", names(ne)[ne > 0], ne[ne > 0]), collapse = ", ")) else ""))
nb <- vapply(PBND, function(k) { x <- r1[[k]]; sum(is.finite(x) & x < 0) }, integer(1))
chk(all(nb == 0L), sprintf("no NEGATIVE bound (%d columns%s)", length(PBND),
    if (any(nb > 0L)) paste0("; ", paste(sprintf("%s %d", names(nb)[nb > 0], nb[nb > 0]), collapse = ", ")) else ""))
dg <- vapply(PBND, function(k) { x <- r1[[k]]; sum((is.finite(x) & x == 0) | is.infinite(x)) }, integer(1))
cat(sprintf("  DEGENERATE BOUNDS (separation): %s\n",
    if (any(dg > 0L)) paste(sprintf("%s %d", names(dg)[dg > 0], dg[dg > 0]), collapse = ", ") else "none"))
Dm1 <- D[D$mr_ok %in% 1L & is.finite(D$fld_H_est2), , drop = FALSE]
chk(nrow(Dm1) > 0 && all(is.finite(Dm1$p_hat_H)),
    sprintf("p-hat(Hhat) recorded on all %d declared replicates with a field block (mean %.3f, share < 0.5 %.3f)",
            nrow(Dm1), mean(Dm1$p_hat_H), mean(Dm1$p_hat_H < 0.5)))
chk(nrow(Dm1) > 0 && all(Dm1$fld_H_lo1s <= Dm1$fld_H_est2), "fld_H_lo1s <= fld_H_est2 on every such replicate")
i2h <- if (nrow(Dm1)) max(abs(log(Dm1$fld_H_lo1s) - (log(Dm1$mr_H_est) - Dm1$fld_H_q95))) else NA_real_
i2c <- if (nrow(f))   max(abs(log(f$fld_Hc_up1s)  - (log(f$mr_Hc_est)  - f$fld_Hc_q05)))  else NA_real_
chk(is.finite(i2h) && i2h <= 1e-9, sprintf("identity: log(lo1s) = log(beta-tilde) - q95 (max |diff| %.2e)", i2h))
chk(is.finite(i2c) && i2c <= 1e-9, sprintf("identity: log(up1s) = log(beta-tilde^c) - q05 (max |diff| %.2e)", i2c))
chk(nrow(f) > 0 && all(f$fld_Hc_est2 <= f$fld_Hc_up1s), "fld_Hc_est2 <= fld_Hc_up1s on every such replicate")
chk(nrow(Dm1) > 0 && all(Dm1$mr_H_lo <= Dm1$mr_H_est & Dm1$mr_H_est <= Dm1$mr_H_hi),
    "mr_H_lo <= mr_H_est <= mr_H_hi on every such replicate")
g1 <- f$fld_joint_gamma; g2 <- f$fld_joint_s_gamma
if (length(g1)) chk(all(g1 >= 0.025 - 1e-12 & g1 <= 0.05 + 1e-12),
                    sprintf("gamma (joint) in [0.025, 0.05] ([%.5f, %.5f])", min(g1), max(g1)))
if (length(g2)) chk(all(g2 >= 0.025 - 1e-12 & g2 <= 0.05 + 1e-12),
                    sprintf("gamma (joint-s) in [0.025, 0.05] ([%.5f, %.5f])", min(g2), max(g2)))
mrfail <- sum(!(D$mr_ok %in% 1L) | !is.finite(D$mr_H_est))
cat(sprintf("  MR failures on declared replicates: %d of %d\n", mrfail, nrow(D)))

cat(sprintf("  FACT timing: fit_mr_secs mean %.1f median %.1f max %.1f | id_secs mean %.2f median %.2f max %.2f | fld_H_secs mean %.2f median %.2f | fld_Hc_secs mean %.2f median %.2f\n",
            mean(r1$fit_mr_secs, na.rm = TRUE), stats::median(r1$fit_mr_secs, na.rm = TRUE), max(r1$fit_mr_secs, na.rm = TRUE),
            mean(r1$id_secs, na.rm = TRUE), stats::median(r1$id_secs, na.rm = TRUE), max(r1$id_secs, na.rm = TRUE),
            mean(r1$fld_H_secs, na.rm = TRUE), stats::median(r1$fld_H_secs, na.rm = TRUE),
            mean(r1$fld_Hc_secs, na.rm = TRUE), stats::median(r1$fld_Hc_secs, na.rm = TRUE)))
if (!identical(method_f, "consistency"))
  cat("  READING: every coverage figure on this path is coverage of beta(H-hat) CONDITIONAL ON THE PROPOSED FAMILY.\n")
cat(sprintf("SMOKE %s target_or_h=%s n=%d: %s\n", mode, target, n, if (ok_all) "PASS" else "FAIL"))
quit(status = if (ok_all) 0 else 1)
