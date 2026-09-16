# Smoke checks for campaign mdgrf (TASK_md_grf_2026-09-16 §1.6).  Transplant of
# scripts_mdsgnb20/smoke_identity.R (its pairing proof and classification, :39-79; its field-s
# wiring checks, :94-107), with the comparator moved from mdf1 to mdsgnb20's committed combined
# bundle and two modes:
#   identity : §1.6(a) FS regression -- the knob at its default, every mdsgnb20 recorder column
#              except *_secs compared on sim_id 1..n (numeric <= 1e-8 relative, NA == NA; sg_def /
#              covs as term sets; other character columns identical); messages / FB reported apart;
#              the template's NEW columns reported apart; zero selection flips is the gate.
#   grf      : §1.6(b) -- n_true identical and the oracle columns <= 1e-8 vs mdsgnb20; meta;
#              the E3 fields filled on every declared replicate; the field-s wiring; p-hat recorded;
#              facts: declared count, sim_id 1's selection, fit_mr_secs and the GRF fit time.
# usage: Rscript smoke_identity.R <md: 40|120|null> <n> <tag> <nsims> <mode: identity|grf>
suppressPackageStartupMessages(library(forestsearch))
args <- commandArgs(TRUE)
md <- args[1]; n <- as.integer(args[2]); tag <- args[3]; nsims <- as.integer(args[4]); mode <- args[5]
stopifnot(mode %in% c("identity", "grf"))
setwd(Sys.getenv("MDSG_DIR", unset = "~/Documents/GitHub/forestsearch/quarto/simulations/actg175/continuous"))
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a
mdtok <- if (md == "null") "mdnull" else sprintf("md%s", md)
meth  <- if (mode == "grf") "grf" else "fs"
stem  <- sprintf("%s_effMaxSG_mr_field_%s_knoise0_n%d_nb20_%s", meth, mdtok, n, tag)
bp    <- file.path("mr_md_harm", paste0(stem, "_d5000"), sprintf("%s_res_1_%d.rds", stem, nsims))
ostem <- sprintf("fs_effMaxSG_mr_field_%s_knoise0_n%d_nb20_mdsgnb20", mdtok, n)
op    <- file.path("mr_md_harm", paste0(ostem, "_d5000"), paste0(ostem, "_combined_1_2000.rds"))
cat(sprintf("== SMOKE %s: md=%s n=%d tag=%s | %s vs mdsgnb20 sim_id 1-%d ==\n", mode, md, n, tag, basename(bp), nsims))
ok_all <- TRUE
chk <- function(cond, msg) { cat(sprintf("  [%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg)); ok_all <<- ok_all && isTRUE(cond) }
chk(file.exists(bp), paste("smoke bundle exists:", bp)); chk(file.exists(op), paste("mdsgnb20 bundle exists:", op))
b <- readRDS(bp); o <- readRDS(op); r1 <- b$results; ro <- o$results
ro <- ro[ro$sim_id %in% seq_len(nsims), ]; r1 <- r1[order(r1$sim_id), ]; ro <- ro[order(ro$sim_id), ]
chk(identical(r1$sim_id, seq_len(nsims)) && identical(ro$sim_id, seq_len(nsims)), sprintf("both carry sim_id 1-%d", nsims))
chk(all(r1$status %in% c("DETECTED", "NO-DETECTION")), sprintf("no CONFIG-ERROR rows (status: %s)", paste(names(table(r1$status)), table(r1$status), collapse = ", ")))
msgc <- c("mr_msg", "err_msg", "fb_err"); fbc <- grep("^fb_", names(ro), value = TRUE)
pre  <- setdiff(names(ro), c(grep("_secs$", names(ro), value = TRUE), msgc, fbc))
newc <- setdiff(names(r1), names(ro))
cat(sprintf("  new columns in this template, reported apart (%d): %s\n", length(newc), paste(newc, collapse = ", ")))
num  <- pre[vapply(pre, function(k) is.numeric(ro[[k]]), logical(1))]
chr  <- setdiff(pre, num)
.terms <- function(v) lapply(strsplit(as.character(v), " & ", fixed = TRUE), sort)
rel <- sapply(num, function(k) { x <- r1[[k]]; y <- ro[[k]]; d <- abs(x - y) / pmax(abs(y), 1e-300); d[is.na(x) & is.na(y)] <- 0; d[xor(is.na(x), is.na(y))] <- Inf; d })
rel <- matrix(rel, nrow = nsims, dimnames = list(NULL, num))
if (mode == "identity") {
  chk(all(pre %in% names(r1)), sprintf("all %d paired mdsgnb20 columns present in the smoke bundle", length(pre)))
  rowmax <- apply(rel, 1, max)
  chr_ok <- sapply(chr, function(k) { x <- r1[[k]]; y <- ro[[k]]
    if (k %in% c("sg_def", "covs")) (mapply(identical, .terms(x), .terms(y)) | (is.na(x) & is.na(y)))
    else (x == y) | (is.na(x) & is.na(y)) })
  chr_ok <- matrix(chr_ok %in% TRUE, nrow = nsims, dimnames = list(NULL, chr))
  same_row <- (rowmax <= 1e-8) & apply(chr_ok, 1, all)
  enum <- r1$sim_id[!same_row]
  cat(sprintf("  pairing: %d of %d rows identical on %d numeric + %d character columns (max rel diff among identical rows %.2e)\n",
              sum(same_row), nsims, length(num), length(chr), if (any(same_row)) max(rowmax[same_row]) else 0))
  ix1 <- match(enum, r1$sim_id); ixo <- match(enum, ro$sim_id)
  .same <- function(k) { a <- r1[[k]][ix1]; z <- ro[[k]][ixo]; d <- abs(a - z) / pmax(abs(z), 1e-300) <= 1e-8; d[is.na(a) & is.na(z)] <- TRUE; d %in% TRUE }
  if (length(enum)) {
    samp_same <- (r1$n_harm[ix1] == ro$n_harm[ixo] | (is.na(r1$n_harm[ix1]) & is.na(ro$n_harm[ixo]))) %in% TRUE & r1$status[ix1] == ro$status[ixo] & .same("nv_H_est") & .same("nv_Hc_est") & .same("or_H_est") & .same("sens")
    mr_same   <- .same("mr_H_est") & .same("mr_H_se_ij") & .same("mr_Hc_est")
    tgt_same  <- .same("betaHhat_H") & .same("betaHhat_Hc")
    cls <- ifelse(!samp_same, ifelse(r1$status[ix1] != ro$status[ixo], "DETECTION FLIP (status differs)", "SELECTION FLIP (sample membership differs)"),
           ifelse(!mr_same, "same selection, MR numerics differ",
           ifelse(!tgt_same, "label tie, super-population target moves", "label tie, no numeric consequence")))
    for (i in seq_along(enum)) {
      dcols <- c(num[rel[ix1[i], ] > 1e-8], chr[!chr_ok[ix1[i], ]])
      cat(sprintf("    sim %d -> %s | committed [%s] | here [%s] | max rel %.2e | differing: %s\n", enum[i], cls[i], ro$sg_def[ixo[i]], r1$sg_def[ix1[i]], rowmax[ix1[i]], paste(dcols, collapse = ",")))
    }
    n_flip <- sum(!samp_same)
    cat(sprintf("  classification: %d selection/detection flips; %d MR-numerics; %d label ties with target move; %d pure label ties\n",
                n_flip, sum(samp_same & !mr_same), sum(samp_same & mr_same & !tgt_same), sum(samp_same & mr_same & tgt_same)))
  } else { n_flip <- 0L; cat("  classification: no enumerated rows\n") }
  chk(n_flip == 0L, sprintf("zero selection flips (%d enumerated rows, all in a mdsgnb20 Gate 2 class)", length(enum)))
  for (k in msgc) cat(sprintf("  reported: %s differs on %d of %d rows\n", k, sum(!((r1[[k]] == ro[[k]]) | (is.na(r1[[k]]) & is.na(ro[[k]])))), nsims))
  cat(sprintf("  reported: FB columns (%d) -- finite fb_H_est: mdsgnb20 %d rows, smoke %d\n", length(fbc), sum(is.finite(ro$fb_H_est)), sum(is.finite(r1$fb_H_est))))
  cat(sprintf("  reported (new columns on the FS path): n_cons_qual finite %d | band_n finite %d | n_family finite %d | admitted_n finite %d | p_hat_sum finite %d | warn_msg non-NA %d | id_secs finite %d\n",
              sum(is.finite(r1$n_cons_qual)), sum(is.finite(r1$band_n)), sum(is.finite(r1$n_family)), sum(is.finite(r1$admitted_n)), sum(is.finite(r1$p_hat_sum)), sum(!is.na(r1$warn_msg)), sum(is.finite(r1$id_secs))))
  chk(identical(b$meta$subgroup_method, "consistency") && grepl("^fs_effMaxSG_", stem), "stem and meta carry the consistency engine (knob at its default)")
  if (any(!is.na(r1$warn_msg))) cat("  FS-path warnings (previously suppressed):", paste(unique(r1$warn_msg[!is.na(r1$warn_msg)]), collapse = " || "), "\n")
} else {
  # ---- (b) same draws: n_true identical; oracle columns <= 1e-8 vs mdsgnb20 ----
  chk(identical(r1$n_true, ro$n_true), "n_true identical on every row")
  orc <- grep("^or_(H|Hc)_", num, value = TRUE)
  orm <- apply(rel[, orc, drop = FALSE], 1, max)
  chk(all(orm <= 1e-8), sprintf("oracle columns (%s) within 1e-8 relative on every row (max %.2e)", paste(orc, collapse = ","), max(orm)))
  m <- b$meta
  cat(sprintf("  meta: subgroup_method %s | dmin_grf %s | grf_selection %s | grf_depth %s | grf_select_statistic %s | focus %s | nbhd %s | rule %s | ci %s | scalec %s | pkg %s | host %s | workers %d | R %s\n",
              m$subgroup_method, format(m$dmin_grf), m$grf_selection, format(m$grf_depth), m$grf_select_statistic, m$sg_focus, format(m$effect_neighborhood), m$selection_rule, m$ci_method, m$field_scale_complement, m$pkg_version, m$hostname, m$n_workers, m$r_version))
  chk(identical(m$subgroup_method, "grf") && isTRUE(all.equal(m$dmin_grf, 30)) && identical(m$sg_focus, "effMaxSG") && isTRUE(all.equal(m$effect_neighborhood, 0.20)) &&
      identical(m$selection_rule, "neighborhood") && identical(m$ci_method, "field") && identical(m$field_scale_complement, "selected") &&
      identical(m$pkg_version, "0.3.5") && identical(m$hostname, "pop-os") && grepl("^grf_effMaxSG_", stem),
      "meta carries identifier grf, dmin.grf 30, focus effMaxSG, band 0.20, the rule, ci_method field, field_scale_complement selected, pkg 0.3.5, host pop-os; stem grf_effMaxSG_")
  D <- r1[r1$detected %in% 1L, ]
  cat(sprintf("  FACT declared: %d of %d | NO-DETECTION %d | admitted_n on non-detections: %s\n", nrow(D), nsims, sum(r1$detected %in% 0L),
              paste(r1$admitted_n[r1$detected %in% 0L], collapse = ",")))
  s1 <- r1[r1$sim_id == 1L, ]
  cat(sprintf("  FACT sim_id 1: sg_def [%s] | n_sel %s | n_harm %s | nv_H_est %.4f | admitted_n %s | n_family %s | sens %.3f ppv %.3f\n",
              s1$sg_def, s1$n_sel, s1$n_harm, s1$nv_H_est, s1$admitted_n, s1$n_family, s1$sens, s1$ppv))
  cat("  FACT S0.3 (quoted): {preanti <= 792.8} & {cd40 > 364}, n 156 / 344, oriented MD 72.04, admitted_n 234, 1,257 enumerated, sens 0.346 ppv 0.404\n")
  chk(all(is.finite(D$admitted_n)) && all(D$admitted_n >= 1L), sprintf("admitted_n filled and >= 1 on every declared replicate (min %s max %s)", min(D$admitted_n), max(D$admitted_n)))
  chk(all(is.finite(r1$admitted_n)), "admitted_n recorded on every replicate (0L where the admitted set was empty)")
  chk(all(is.finite(D$n_family[D$mr_ok %in% 1L])), sprintf("n_family filled on every declared replicate with a gate (values %s)", paste(range(D$n_family, na.rm = TRUE), collapse = "-")))
  chk(all(is.finite(D$p_hat_H[D$mr_ok %in% 1L])) && all(is.finite(D$p_hat_sum[D$mr_ok %in% 1L])), sprintf("p_hat_H and p_hat_sum recorded (p_hat_H mean %.3f)", mean(D$p_hat_H, na.rm = TRUE)))
  cat(sprintf("  STRUCTURAL n_cons_qual all-NA: %s | band_n all-NA: %s (GRF has no consistency screen)\n", all(is.na(r1$n_cons_qual)), all(is.na(r1$band_n))))
  cat(sprintf("  FACT the floor as applied: meta dmin_grf = %s on the DR-score harm effect; effect_threshold = %s on the harm-oriented MD (admission)\n", format(m$dmin_grf), format(m$effect_threshold)))
  wm <- r1$warn_msg[!is.na(r1$warn_msg)]
  cat(sprintf("  WARNINGS (verbatim, distinct, with per-row counts on %d of %d rows):\n", length(wm), nsims))
  for (w in unique(wm)) cat("    ", w, "  [rows:", sum(wm == w), "]\n")
  cat(sprintf("  FACT timing: fit_mr_secs mean %.1f median %.1f max %.1f | id_secs (GRF fit, incl. re-selection) mean %.2f median %.2f max %.2f | fld_H_secs mean %.1f | fld_Hc_secs mean %.2f\n",
              mean(r1$fit_mr_secs), median(r1$fit_mr_secs), max(r1$fit_mr_secs), mean(r1$id_secs, na.rm = TRUE), median(r1$id_secs, na.rm = TRUE), max(r1$id_secs, na.rm = TRUE),
              mean(r1$fld_H_secs, na.rm = TRUE), mean(r1$fld_Hc_secs, na.rm = TRUE)))
  cat(sprintf("  FACT vs mdsgnb20 (sim 1-%d): mean |Hhat| (n_harm) %.1f vs %.1f | sens %.3f vs %.3f | ppv %.3f vs %.3f | declared %d vs %d\n", nsims,
              mean(r1$n_harm, na.rm = TRUE), mean(ro$n_harm, na.rm = TRUE), mean(r1$sens, na.rm = TRUE), mean(ro$sens, na.rm = TRUE), mean(r1$ppv, na.rm = TRUE), mean(ro$ppv, na.rm = TRUE), sum(r1$detected %in% 1L), sum(ro$detected %in% 1L)))
}
# ---- field-s wiring on every replicate whose complement field block was filled (scripts_mdsgnb20/smoke_identity.R:94-107) ----
f <- r1[is.finite(r1$fld_Hc_est2), ]
sC <- c("fld_Hc_est2_s","fld_Hc_up1s_s","fld_Hc_lo1s_s","fld_Hc_lo2s_s","fld_Hc_hi2s_s","fld_Hc_lo_se_s","fld_Hc_hi_se_s","fld_Hc_se_s","fld_Hc_lam_mean_s")
jC <- c("fld_joint_s_gamma","fld_joint_s_prob","fld_joint_s_loH","fld_joint_s_upHc","fld_joint_s_bonf_loH","fld_joint_s_bonf_upHc","fld_joint_s_bonf_prob","fld_joint_s_corr","fld_joint_s_n")
chk(all(c(sC, jC) %in% names(r1)), "nine fld_Hc_*_s and nine fld_joint_s_* columns present")
chk(nrow(f) > 0 && all(vapply(sC, function(k) all(is.finite(f[[k]])), logical(1))), sprintf("nine fld_Hc_*_s finite on all %d filled replicates", nrow(f)))
chk(nrow(f) > 0 && all(vapply(jC, function(k) all(is.finite(f[[k]])), logical(1))), sprintf("nine fld_joint_s_* finite on all %d filled replicates", nrow(f)))
chk(all(f$fld_Hc_lo1s_s <= f$fld_Hc_up1s_s) && all(f$fld_Hc_lo2s_s <= f$fld_Hc_hi2s_s), "fld_Hc_lo1s_s <= fld_Hc_up1s_s and fld_Hc_lo2s_s <= fld_Hc_hi2s_s")
agree <- f$fld_joint_n == f$fld_joint_s_n
dj <- if (any(agree)) max(abs(f$fld_joint_bonf_loH[agree] - f$fld_joint_s_bonf_loH[agree])) else 0
chk(dj <= 1e-12, sprintf("Bonferroni harm bound identical between joint and joint_s where draw counts agree (%d of %d rows agree; max |diff| %.2e)", sum(agree), nrow(f), dj))
idn <- max(abs((f$fld_Hc_est2_s + f$fld_Hc_lam_mean_s) - (f$fld_Hc_est2 + f$fld_Hc_lam_mean)))
chk(idn <= 1e-9, sprintf("field-s inverted around the same beta-tilde^c (max |diff| %.2e)", idn))
chk(identical(b$meta$field_scale_complement, "selected"), sprintf("meta field_scale_complement = %s", b$meta$field_scale_complement %||% "<absent>"))
cat(sprintf("SMOKE %s md=%s n=%d: %s\n", mode, md, n, if (ok_all) "PASS" else "FAIL"))
quit(status = if (ok_all) 0 else 1)
