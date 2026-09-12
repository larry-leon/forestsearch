# PART C -- FS extraction for the manuscript (TASK_grfmr_campaign_2026-09-11).
#
# READING COMMITTED BUNDLES ONLY.  No re-run, no new simulation, no recorder
# change.  Every number below comes from columns that are already in the
# committed FS comparator bundles.
#
# Coverage: the committed FS grid at both prevalences, harm and null --
#   12.4% : p12ext / tier2, sg_focus = maxeffCons, eps 0.10
#   31%   : cert20 / e1stud, sg_focus = effMaxSG,  eps 0.20
# The comparator's sg_focus and eps are carried beside every row, because the
# two prevalences are NOT run under the same criterion.
#
# WHAT IS AND IS NOT DERIVABLE.
#   * sens / spec / ppv / npv are recorded PER REPLICATE (template
#     .classify(), sim_fs_...template.qmd:811-822) as tp/(tp+fn), tn/(tn+fp),
#     tp/(tp+fp), tn/(tn+fn) over the replicate's own subjects.  The 2x2 itself
#     is not a recorder column, but it is EXACTLY recoverable: n_sel = tp + fp,
#     n_true = tp + fn, and tp = sens * n_true = ppv * n_sel, so
#     tp, fp, fn, tn = n_sample - tp - fp - fn all follow.  The reconstruction
#     is CHECKED here (sens*n_true vs ppv*n_sel) and the max discrepancy is
#     reported; a Wilson interval is then put on the POOLED subject-level
#     counts, which is the quantity a Wilson interval is defined for.  The
#     replicate-mean rate is reported beside it, unbracketed, because that is
#     the form the DINA tables use and a Wilson interval on a mean of
#     proportions would not be one.
#   * mean |Hhat| is n_sel; |H| is n_true.  Both are recorder columns.
#   * Bound location uses the SAME columns as the DINA location tables
#     (blockA_rest.R:63-76): fld_H_lo1s, betaHhat_H, fld_H_est2, nv_H_est,
#     truth$marg_H.
#   * Every rate and every location figure is computed on the DETECTED
#     replicates -- the set on which the quantities exist.  n_eval is stated.
SCRATCH <- Sys.getenv("DINAMR_SCRATCH", unset = ".")
QMD_DIR <- Sys.getenv("DINAMR_QMD_DIR", unset = "..")
R <- file.path(QMD_DIR, "results/")

wilson <- function(x, n, conf = 0.95) {
  if (!is.finite(n) || n <= 0) return(c(NA_real_, NA_real_))
  z <- stats::qnorm(1 - (1 - conf)/2); p <- x/n
  c((p + z^2/(2*n) - z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n),
    (p + z^2/(2*n) + z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n))
}
fmtw <- function(x, n) sprintf("%.4f [%.4f, %.4f]", x/n, wilson(x,n)[1], wilson(x,n)[2])

# --- the committed FS grid --------------------------------------------------
grid <- rbind(
  expand.grid(hr = c(1.50,1.75,1.00), n = c(500L,1000L,1500L), prev = "12.4%",
              stringsAsFactors = FALSE),
  expand.grid(hr = c(1.50,1.75,1.00), n = c(500L,1000L,1500L), prev = "31%",
              stringsAsFactors = FALSE))
fsfile <- function(hr, n, prev) {
  if (prev == "12.4%") {
    camp <- if (abs(hr-1.75) < 1e-9) "tier2" else if (abs(hr-1.00) < 1e-9 && n == 500L) "tier2" else "p12ext"
    sprintf("%sfs_maxeffCons_fb_mr_field_m1_h%03d_knoise0_n%d_%s_combined_1_2000.rds", R, round(100*hr), n, camp)
  } else {
    camp <- if (n == 500L && abs(hr-1.00) > 1e-9) "e1stud" else "cert20"
    sprintf("%sfs_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d_z1q60_nb20_%s_combined_1_2000.rds", R, round(100*hr), n, camp)
  }
}
lab <- function(g) sprintf("%s HR %.2f n %d", g$prev, g$hr, g$n)

B <- list()
for (i in seq_len(nrow(grid))) {
  f <- fsfile(grid$hr[i], grid$n[i], grid$prev[i])
  if (!file.exists(f)) { cat("MISSING:", basename(f), "\n"); next }
  b <- readRDS(f)
  B[[lab(grid[i,])]] <- list(b = b, f = basename(f), g = grid[i,])
}
cat(sprintf("Bundles read: %d of %d grid cells.\n\n", length(B), nrow(grid)))

## ===========================================================================
cat("## 1. CLASSIFICATION AGAINST THE PLANTED REGION\n\n")
cat("Detected replicates only.  Pooled rates put a Wilson interval on the\n")
cat("subject-level 2x2 summed over replicates; rep-mean is the mean of the\n")
cat("per-replicate rates, the form the DINA tables use.\n\n")
C1 <- list(); recon_max <- 0
for (k in names(B)) {
  b <- B[[k]]$b; r <- b$results; m <- b$meta
  d <- r[r$detected %in% 1L & is.finite(r$sens), , drop = FALSE]
  if (!nrow(d)) next
  tp1 <- d$sens * d$n_true; tp2 <- d$ppv * d$n_sel
  recon_max <- max(recon_max, max(abs(tp1 - tp2), na.rm = TRUE))
  tp <- round(tp1); fp <- d$n_sel - tp; fn <- d$n_true - tp
  tn <- m$n_sample - tp - fp - fn
  TP <- sum(tp); FP <- sum(fp); FN <- sum(fn); TN <- sum(tn)
  C1[[k]] <- data.frame(
    cell = k, campaign = sub("^.*_nb20_|^.*knoise0_n[0-9]+_", "", sub("_combined.*$","",B[[k]]$f)),
    focus = m$sg_focus, eps = m$effect_neighborhood,
    n_eval = nrow(d), detection = mean(r$detected %in% 1L),
    sens_pool = fmtw(TP, TP+FN), spec_pool = fmtw(TN, TN+FP),
    ppv_pool  = fmtw(TP, TP+FP), npv_pool  = fmtw(TN, TN+FN),
    sens_rep = mean(d$sens), spec_rep = mean(d$spec),
    ppv_rep = mean(d$ppv), npv_rep = mean(d$npv),
    mean_Hhat = mean(d$n_sel), mean_H = mean(d$n_true),
    ratio_Hhat_H = mean(d$n_sel)/mean(d$n_true),
    med_ratio_paired = stats::median(d$n_sel/d$n_true),
    stringsAsFactors = FALSE)
}
C1 <- do.call(rbind, C1)
cat("### 1a. pooled rates with Wilson intervals\n\n")
print(C1[, c("cell","campaign","focus","eps","n_eval","detection",
             "sens_pool","spec_pool","ppv_pool","npv_pool")], row.names = FALSE)
cat("\n### 1b. replicate-mean rates, and |Hhat| against |H|\n\n")
print(C1[, c("cell","focus","eps","sens_rep","spec_rep","ppv_rep","npv_rep",
             "mean_Hhat","mean_H","ratio_Hhat_H","med_ratio_paired")],
      row.names = FALSE, digits = 4)
cat(sprintf("\n2x2 reconstruction check: max |sens*n_true - ppv*n_sel| over every row of every cell = %.3g\n",
            recon_max))
cat("(zero to floating point means the reconstruction is exact, not an approximation)\n")

## ===========================================================================
cat("\n\n## 2. BOUND LOCATION ON THE HR SCALE\n\n")
cat("Same columns as the DINA location tables (scripts_dinamr/blockA_rest.R:63-76):\n")
cat("fld_H_lo1s (field lower bound), betaHhat_H (realized theta(Hhat) on the\n")
cat("super-population), fld_H_est2, nv_H_est, truth$marg_H.  Detected replicates only.\n\n")
C2 <- list()
for (k in names(B)) {
  b <- B[[k]]$b; r <- b$results; m <- b$meta
  d <- r[r$detected %in% 1L & is.finite(r$fld_H_lo1s) & is.finite(r$betaHhat_H), , drop = FALSE]
  if (!nrow(d)) next
  s100 <- sum(d$fld_H_lo1s >= 1.00); s125 <- sum(d$fld_H_lo1s >= 1.25)
  C2[[k]] <- data.frame(
    cell = k, focus = m$sg_focus, eps = m$effect_neighborhood, n_eval = nrow(d),
    med_bound = stats::median(d$fld_H_lo1s),
    med_theta = stats::median(d$betaHhat_H),
    med_est2  = stats::median(d$fld_H_est2),
    med_naive = stats::median(d$nv_H_est),
    bound_minus_theta = stats::median(d$fld_H_lo1s) - stats::median(d$betaHhat_H),
    bound_over_theta  = stats::median(d$fld_H_lo1s) / stats::median(d$betaHhat_H),
    med_paired_ratio  = stats::median(d$fld_H_lo1s / d$betaHhat_H),
    planted_marg_H = unname(b$truth$marg_H),
    share_ge_100 = fmtw(s100, nrow(d)),
    share_ge_125 = fmtw(s125, nrow(d)),
    stringsAsFactors = FALSE)
}
C2 <- do.call(rbind, C2)
print(C2[, c("cell","focus","eps","n_eval","med_bound","med_theta","med_est2","med_naive",
             "planted_marg_H")], row.names = FALSE, digits = 4)
cat("\n### 2b. the gap, three ways, and the location shares with Wilson intervals\n\n")
print(C2[, c("cell","focus","eps","bound_minus_theta","bound_over_theta",
             "med_paired_ratio","share_ge_100","share_ge_125")],
      row.names = FALSE, digits = 4)

saveRDS(list(classification = C1, location = C2, grid = grid),
        file.path(SCRATCH, "fs_extraction.rds"))
cat("\n\nWritten: fs_extraction.rds\n")
