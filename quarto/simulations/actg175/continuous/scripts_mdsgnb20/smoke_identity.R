# Smoke identity and field-s wiring checks (TASK_md_field_rerun_2026-09-15 §1.6).
# Transplant of the pairing proof and classification in scripts_mdf1/gate2_check.R
# (:22-53; the script that wrote mdf1's gate2_flips.txt), repointed at a smoke
# batch bundle vs sim_id 1..n of mdf1's committed combined bundle, with:
#   - every recorder column except *_secs compared (numeric: <= 1e-8 relative,
#     NA == NA; sg_def / covs as term sets; other character columns identical);
#     message (mr_msg, err_msg, fb_err) and joined FB columns (fb_*) are
#     REPORTED separately, not paired -- mdf1's Gate 2 convention excluded
#     "timings, messages and the FB columns" (REPORT_continuous_field_gate2_2026-09-07.md:5);
#   - classification exactly as gate2_check.R:41-51 (selection flip / detection
#     flip / MR-numerics / label tie with target move / pure label tie);
#   - §1.6(b) field-s wiring checks on every replicate whose complement field
#     block was filled;
#   - mode "rule": §1.6(c) -- n_true identical, oracle columns <= 1e-8, and the
#     |Hhat| grew / stayed / shrank count against mdf1 (a fact, not a gate).
# usage: Rscript smoke_identity.R <md: 40|120|null> <n> <focus> <nbhd> <tag> <nsims> <mode: identity|rule>
suppressPackageStartupMessages(library(forestsearch))
args <- commandArgs(TRUE)
md <- args[1]; n <- as.integer(args[2]); focus <- args[3]; nbhd <- as.numeric(args[4]); tag <- args[5]
nsims <- as.integer(args[6]); mode <- args[7]
stopifnot(mode %in% c("identity", "rule"))
setwd(Sys.getenv("MDSG_DIR", unset = "~/Documents/GitHub/forestsearch/quarto/simulations/actg175/continuous"))
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a
focus_tag <- forestsearch::fs_focus_tag("consistency", focus)
nbhd_tag  <- if (abs(nbhd - 0.10) < 1e-12) "" else sprintf("_nb%02d", round(100 * nbhd))
mdtok <- if (md == "null") "mdnull" else sprintf("md%s", md)
stem <- sprintf("fs_%s_mr_field_%s_knoise0_n%d%s_%s", focus_tag, mdtok, n, nbhd_tag, tag)
bp   <- file.path("mr_md_harm", paste0(stem, "_d5000"), sprintf("%s_res_1_%d.rds", stem, nsims))
ostem <- sprintf("fs_maxeffCons_mr_field_%s_knoise0_n%d_mdf1", mdtok, n)
op   <- file.path("mr_md_harm", paste0(ostem, "_d5000"), paste0(ostem, "_combined_1_2000.rds"))
cat(sprintf("== SMOKE %s: md=%s n=%d focus=%s nbhd=%.2f tag=%s | %s vs mdf1 sim_id 1-%d ==\n", mode, md, n, focus, nbhd, tag, basename(bp), nsims))
ok_all <- TRUE
chk <- function(cond, msg) { cat(sprintf("  [%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg)); ok_all <<- ok_all && isTRUE(cond) }
chk(file.exists(bp), paste("smoke bundle exists:", bp)); chk(file.exists(op), paste("mdf1 bundle exists:", op))
b <- readRDS(bp); o <- readRDS(op); r1 <- b$results; ro <- o$results
ro <- ro[ro$sim_id %in% seq_len(nsims), ]; r1 <- r1[order(r1$sim_id), ]; ro <- ro[order(ro$sim_id), ]
chk(identical(r1$sim_id, seq_len(nsims)) && identical(ro$sim_id, seq_len(nsims)), sprintf("both carry sim_id 1-%d", nsims))
chk(all(r1$status %in% c("DETECTED", "NO-DETECTION")), sprintf("no CONFIG-ERROR rows (status: %s)", paste(names(table(r1$status)), table(r1$status), collapse = ", ")))
# ---- pairing set: every mdf1 recorder column except *_secs; messages / FB reported apart ----
msgc <- c("mr_msg", "err_msg", "fb_err"); fbc <- grep("^fb_", names(ro), value = TRUE)
pre  <- setdiff(names(ro), c(grep("_secs$", names(ro), value = TRUE), msgc, fbc))
chk(all(pre %in% names(r1)), sprintf("all %d paired mdf1 columns present in the smoke bundle", length(pre)))
num  <- pre[vapply(pre, function(k) is.numeric(ro[[k]]), logical(1))]
chr  <- setdiff(pre, num)
.terms <- function(v) lapply(strsplit(as.character(v), " & ", fixed = TRUE), sort)
rel <- sapply(num, function(k) { x <- r1[[k]]; y <- ro[[k]]; d <- abs(x - y) / pmax(abs(y), 1e-300); d[is.na(x) & is.na(y)] <- 0; d[xor(is.na(x), is.na(y))] <- Inf; d })
rel <- matrix(rel, nrow = nsims, dimnames = list(NULL, num))
rowmax <- apply(rel, 1, max)
chr_ok <- sapply(chr, function(k) { x <- r1[[k]]; y <- ro[[k]]
  if (k %in% c("sg_def", "covs")) (mapply(identical, .terms(x), .terms(y)) | (is.na(x) & is.na(y)))
  else (x == y) | (is.na(x) & is.na(y)) })
chr_ok <- matrix(chr_ok %in% TRUE, nrow = nsims, dimnames = list(NULL, chr))
same_row <- (rowmax <= 1e-8) & apply(chr_ok, 1, all)
if (mode == "identity") {
  enum <- r1$sim_id[!same_row]
  cat(sprintf("  pairing: %d of %d rows identical on %d numeric + %d character columns (max rel diff among identical rows %.2e); rule-string-order-only differences on %d rows\n",
              sum(same_row), nsims, length(num), length(chr), if (any(same_row)) max(rowmax[same_row]) else 0,
              sum((apply(chr_ok[, c("sg_def", "covs"), drop = FALSE], 1, all) & (r1$sg_def != ro$sg_def)) %in% TRUE)))
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
  chk(n_flip == 0L, sprintf("zero selection flips (%d enumerated rows, all in a mdf1 Gate 2 class)", length(enum)))
  # messages / FB, reported apart (not paired; mdf1 convention)
  for (k in msgc) cat(sprintf("  reported: %s differs on %d of %d rows\n", k, sum(!((r1[[k]] == ro[[k]]) | (is.na(r1[[k]]) & is.na(ro[[k]])))), nsims))
  cat(sprintf("  reported: FB columns (%d) -- mdf1 has finite fb_H_est on %d rows (joined), smoke on %d\n", length(fbc), sum(is.finite(ro$fb_H_est)), sum(is.finite(r1$fb_H_est))))
} else {
  # ---- (c) same draws under the campaign rule: n_true identical; oracle columns <= 1e-8 ----
  chk(identical(r1$n_true, ro$n_true), "n_true identical on every row")
  orc <- grep("^or_(H|Hc)_", num, value = TRUE)
  orm <- apply(rel[, orc, drop = FALSE], 1, max)
  chk(all(orm <= 1e-8), sprintf("oracle columns (%s) within 1e-8 relative on every row (max %.2e)", paste(orc, collapse = ","), max(orm)))
  d <- r1$n_harm - ro$n_harm
  cat(sprintf("  FACT |Hhat| vs mdf1 (n_harm): grew %d, stayed %d, shrank %d, undefined %d | mean |Hhat| here %.1f vs mdf1 %.1f | same rule string on %d rows\n",
              sum(d > 0, na.rm = TRUE), sum(d == 0, na.rm = TRUE), sum(d < 0, na.rm = TRUE), sum(is.na(d)),
              mean(r1$n_harm, na.rm = TRUE), mean(ro$n_harm, na.rm = TRUE), sum(mapply(identical, .terms(r1$sg_def), .terms(ro$sg_def)))))
  cat(sprintf("  FACT n_sel: grew %d, stayed %d, shrank %d\n", sum((r1$n_sel - ro$n_sel) > 0, na.rm = TRUE), sum((r1$n_sel - ro$n_sel) == 0, na.rm = TRUE), sum((r1$n_sel - ro$n_sel) < 0, na.rm = TRUE)))
  cat(sprintf("  FACT stem: %s | meta sg_focus %s | effect_neighborhood %s | selection_rule %s\n", stem, b$meta$sg_focus, format(b$meta$effect_neighborhood), b$meta$selection_rule))
  chk(identical(b$meta$sg_focus, focus) && isTRUE(all.equal(b$meta$effect_neighborhood, nbhd)) && grepl(sprintf("^fs_%s_", focus_tag), stem) && (nbhd_tag == "" || grepl(nbhd_tag, stem, fixed = TRUE)), "stem and meta carry the rule")
}
# ---- (b) field-s wiring on every replicate whose complement field block was filled ----
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
cat(sprintf("  meta: pkg %s | host %s | R %s | workers %d | ci %s | focus %s | nbhd %s | rule %s | scalec %s | built %s\n", b$meta$pkg_version, b$meta$hostname, b$meta$r_version, b$meta$n_workers, b$meta$ci_method, b$meta$sg_focus, format(b$meta$effect_neighborhood), b$meta$selection_rule, b$meta$field_scale_complement, format(b$meta$built_at)))
cat(sprintf("  timing: fit_mr_secs mean %.1f median %.1f max %.1f | fld_H_secs mean %.1f | fld_Hc_secs mean %.2f | detected %d/%d\n", mean(r1$fit_mr_secs), median(r1$fit_mr_secs), max(r1$fit_mr_secs), mean(r1$fld_H_secs, na.rm = TRUE), mean(r1$fld_Hc_secs, na.rm = TRUE), sum(r1$detected %in% 1L), nsims))
cat(sprintf("SMOKE %s md=%s n=%d: %s\n", mode, md, n, if (ok_all) "PASS" else "FAIL"))
quit(status = if (ok_all) 0 else 1)
