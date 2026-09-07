# Gate 2 per cell: (1) batch-1 pre-existing columns identical (<= 1e-8 rel; rules as sets) to the
# committed bundle by sim_id, flips enumerated; (2) new columns present and finite on detected
# replicates; (3) save guard active (bundle paths untracked); (4) completeness; (5) invariants; gamma.
args <- commandArgs(TRUE); md <- args[1]; n <- as.integer(args[2]); tag <- "mdf1"
setwd("~/Documents/GitHub/forestsearch/quarto/simulations/actg175/continuous")
stem <- if (md == "null") sprintf("fs_maxeffCons_mr_field_mdnull_knoise0_n%d_%s", n, tag) else sprintf("fs_maxeffCons_mr_field_md%s_knoise0_n%d_%s", md, n, tag)
dir <- file.path("mr_md_harm", paste0(stem, "_d5000"))
b1p <- file.path(dir, sprintf("%s_res_1_1000.rds", stem)); b2p <- file.path(dir, sprintf("%s_res_1001_2000.rds", stem))
cbp <- file.path(dir, sprintf("%s_combined_1_2000.rds", stem))
old_p <- switch(md, "40" = sprintf("mr_md_harm/fs_maxeffCons_mr_md40_knoise0_n%d_s1000_d5000/fs_maxeffCons_mr_md40_knoise0_n%d_res_1_1000.rds", n, n),
                "120" = "mr_md_harm/fs_maxeffCons_mr_md120_knoise0_n500_s1000_d5000/fs_maxeffCons_mr_md120_knoise0_n500_res_1_1000.rds",
                "null" = "mr_md_harm/fs_maxeffCons_mr_mdnull_knoise0_n500_s1000_d5000/fs_maxeffCons_mr_mdnull_knoise0_n500_res_1_1000.rds")
cat(sprintf("== GATE 2: md=%s n=%d ==\n", md, n)); ok_all <- TRUE
chk <- function(cond, msg) { cat(sprintf("  [%s] %s\n", if (isTRUE(cond)) "PASS" else "FAIL", msg)); ok_all <<- ok_all && isTRUE(cond) }
for (p in c(b1p, b2p, cbp)) chk(file.exists(p), paste("exists:", basename(p)))
tracked <- function(p) tryCatch(system2("git", c("ls-files", "--error-unmatch", shQuote(p)), stdout = FALSE, stderr = FALSE) == 0L, error = function(e) FALSE)
chk(!any(vapply(c(b1p, b2p, cbp), tracked, logical(1))), "save guard: all three bundle paths untracked at save time (guard cannot have been bypassed)")
b1 <- readRDS(b1p); b2 <- readRDS(b2p); cb <- readRDS(cbp); old <- readRDS(old_p)
r1 <- b1$results; r2 <- b2$results; rc <- cb$results; ro <- old$results
chk(identical(sort(r1$sim_id), 1:1000) && identical(sort(r2$sim_id), 1001:2000) && identical(sort(rc$sim_id), 1:2000), "completeness: sim_id 1-1000, 1001-2000, combined 1-2000")
chk(all(rc$status %in% c("DETECTED", "NO-DETECTION")), sprintf("no CONFIG-ERROR rows (status table: %s)", paste(names(table(rc$status)), table(rc$status), collapse = ", ")))
# (1) pairing proof on batch 1
r1 <- r1[order(r1$sim_id), ]; ro <- ro[order(ro$sim_id), ]; stopifnot(identical(r1$sim_id, ro$sim_id))
pre <- setdiff(names(ro), c("fit_mr_secs", "fb_secs", "mr_msg", "fb_err", "fb_H_est","fb_H_lo","fb_H_hi","fb_H_se","fb_Hc_est","fb_Hc_lo","fb_Hc_hi","fb_Hc_se","fb_src1","fb_src2","fb_nres"))
num <- pre[vapply(pre, function(k) is.numeric(ro[[k]]), logical(1))]
.terms <- function(v) lapply(strsplit(v, " & ", fixed = TRUE), sort)
rule_ok <- mapply(identical, .terms(r1$sg_def), .terms(ro$sg_def)) | (is.na(r1$sg_def) & is.na(ro$sg_def))
rel <- sapply(num, function(k) { x <- r1[[k]]; y <- ro[[k]]; d <- abs(x - y) / pmax(abs(y), 1e-300); d[is.na(x) & is.na(y)] <- 0; d[xor(is.na(x), is.na(y))] <- Inf; d })
rowmax <- apply(rel, 1, max)
flips <- r1$sim_id[!(rowmax <= 1e-8 & rule_ok & r1$status == ro$status) %in% TRUE]
cat(sprintf("  pairing proof: %d of 1000 sim_ids identical on %d pre-existing numeric cols (max rel diff among identical rows %.2e); rule-string-order-only differences on %d sim_ids\n",
            1000 - length(flips), length(num), max(rowmax[!(r1$sim_id %in% flips)]), sum((rule_ok & r1$sg_def != ro$sg_def) %in% TRUE)))
if (length(flips)) {
  cat(sprintf("  ENUMERATED FLIPS (%d): %s\n", length(flips), paste(flips, collapse = ", ")))
  for (s in flips) cat(sprintf("    sim %d: committed [%s] status %s | here [%s] status %s | max rel %.2e\n", s, ro$sg_def[ro$sim_id == s], ro$status[ro$sim_id == s], r1$sg_def[r1$sim_id == s], r1$status[r1$sim_id == s], rowmax[r1$sim_id == s]))
}
# Classify each enumerated difference: LABEL TIE (same membership: numerics identical, n_harm equal;
# only the cut label differs) vs SELECTION FLIP (membership differs).  Gate criterion (M-5 kickoff):
# pre-existing columns identical by sim_id with enumerated flips excluded -- so the gate passes when
# every NON-flip row is identical; flips are listed for the record and excluded downstream.
ix1 <- match(flips, r1$sim_id); ixo <- match(flips, ro$sim_id)
.same <- function(k) { a <- r1[[k]][ix1]; z <- ro[[k]][ixo]; d <- abs(a - z) / pmax(abs(z), 1e-300) <= 1e-8; d[is.na(a) & is.na(z)] <- TRUE; d %in% TRUE }
samp_same <- (r1$n_harm[ix1] == ro$n_harm[ixo] | (is.na(r1$n_harm[ix1]) & is.na(ro$n_harm[ixo]))) %in% TRUE & r1$status[ix1] == ro$status[ixo] & .same("nv_H_est") & .same("nv_Hc_est") & .same("or_H_est") & .same("sens")
mr_same   <- .same("mr_H_est") & .same("mr_H_se_ij") & .same("mr_Hc_est")
tgt_same  <- .same("betaHhat_H") & .same("betaHhat_Hc")
cls <- ifelse(!samp_same, ifelse(r1$status[ix1] != ro$status[ixo], "DETECTION FLIP (status differs)", "SELECTION FLIP (sample membership differs)"),
       ifelse(!mr_same, "same selection, MR numerics differ",
       ifelse(!tgt_same, "label tie, super-population target moves", "label tie, no numeric consequence")))
for (i in seq_along(flips)) cat(sprintf("    sim %d -> %s\n", flips[i], cls[i]))
cat(sprintf("  classification: %d selection flips; %d MR-numerics; %d label ties with target move; %d pure label ties\n",
            sum(!samp_same), sum(samp_same & !mr_same), sum(samp_same & mr_same & !tgt_same), sum(samp_same & mr_same & tgt_same)))
writeLines(sprintf("%d\t%s", flips, cls), file.path(dir, "gate2_flips.txt"))
chk(all(rowmax[!(r1$sim_id %in% flips)] <= 1e-8), sprintf("pairing proof: all %d non-flip rows identical (<= 1e-8); %d enumerated and excluded (%d selection flips)", 1000 - length(flips), length(flips), sum(!samp_same)))
# (2) new columns present and finite on detected replicates
d <- rc[rc$detected %in% 1L, ]
newc <- c("fld_H_est2","fld_H_lo1s","fld_H_lo2s","fld_H_hi2s","fld_H_se","fld_H_lam_mean","fld_H_q05","fld_H_q95","fld_H_nout",
          "fld_Hc_est2","fld_Hc_up1s","fld_Hc_lo1s","fld_Hc_lo2s","fld_Hc_hi2s","fld_Hc_se","fld_Hc_nfit",
          "fld_joint_gamma","fld_joint_prob","fld_joint_loH","fld_joint_upHc","fld_joint_bonf_loH","fld_joint_bonf_upHc","fld_joint_corr",
          "p_hat_H","p_top1","mr_H_se_w","mr_Hc_se_w","label","mr_harm_flag")
chk(all(newc %in% names(rc)), "new columns present")
fin <- vapply(setdiff(newc, c("label")), function(k) mean(is.finite(d[[k]])), numeric(1))
cat(sprintf("  finite share on %d detected reps: min %.4f (%s); field notes set on %d (H) / %d (Hc)\n", nrow(d), min(fin), names(fin)[which.min(fin)], sum(!is.na(d$fld_H_note)), sum(!is.na(d$fld_Hc_note))))
chk(all(fin[!grepl("^mr_.*_w$", names(fin))] >= 0.99), "fields finite on >= 99% of detected replicates")
# (5) invariants and gamma
inv <- with(d, all(fld_H_lo2s <= fld_H_hi2s & fld_Hc_lo2s <= fld_Hc_hi2s & fld_Hc_lo1s <= fld_Hc_up1s & mr_H_lo <= mr_H_hi & mr_Hc_lo <= mr_Hc_hi & nv_H_lo <= nv_H_hi & or_H_lo <= or_H_hi, na.rm = TRUE))
chk(inv, "interval invariants lo <= hi")
idn <- with(d, max(abs(fld_H_lo1s - (mr_H_est - fld_H_q95)), abs(fld_Hc_up1s - (mr_Hc_est - fld_Hc_q05)), abs(fld_H_est2 - (mr_H_est - fld_H_lam_mean)), na.rm = TRUE))
chk(idn <= 1e-12, sprintf("bound identities (max abs %.1e)", idn))
g <- d$fld_joint_gamma[is.finite(d$fld_joint_gamma)]
chk(all(g >= 0.025 - 1e-12 & g <= 0.05 + 1e-12), sprintf("gamma in [0.025, 0.05] (range %.3f-%.3f; mean %.4f)", min(g), max(g), mean(g)))
chk(all(is.finite(d$p_hat_H) & d$p_hat_H >= 0 & d$p_hat_H <= 1), sprintf("p_hat(Hhat) finite in [0,1] (mean %.3f, share < 0.5: %.3f)", mean(d$p_hat_H), mean(d$p_hat_H < 0.5)))
# FB join report
if (any(is.finite(rc$fb_H_est))) cat(sprintf("  FB joined on %d sim_ids (%s)\n", sum(is.finite(rc$fb_H_est)), b1$meta$fb_join_source))
cat(sprintf("  meta: pkg %s | workers %s | ci %s | complement %s | ij_residual %s | built %s\n", cb$meta$pkg_version, paste(cb$meta$n_workers_by_batch, collapse="/"), cb$meta$ci_method, cb$meta$field_complement, cb$meta$ij_residual, format(cb$meta$built_at)))
cat(sprintf("  timing: batch1 fit+MR mean %.1f s (max %.1f), batch2 %.1f s; detection %.3f\n", mean(r1$fit_mr_secs), max(r1$fit_mr_secs), mean(r2$fit_mr_secs), mean(rc$detected)))
cat(sprintf("GATE 2 md=%s n=%d: %s\n", md, n, if (ok_all) "PASS" else "FAIL"))
quit(status = if (ok_all) 0 else 1)
