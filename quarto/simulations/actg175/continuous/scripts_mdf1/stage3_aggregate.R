# Stage 3 cross-cell aggregation from the combined bundles: Table-2 layout rows per cell/block,
# bound-location shares, joint pair, regime diagnostics, p-hat, and the display table
# (fs_sim_bias_coverage scale = "identity").  Writes markdown fragments + the display PNG.
suppressPackageStartupMessages({ library(forestsearch); library(ggplot2) })
args <- commandArgs(TRUE); tag <- if (length(args)) args[1] else "mdf1"
setwd("~/Documents/GitHub/forestsearch/quarto/simulations/actg175/continuous")
out_dir <- Sys.getenv("OUT_DIR", "."); dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
cells <- list(list(lab = "md40 n500",  md = "40",  n = 500L), list(lab = "md120 n500", md = "120", n = 500L),
              list(lab = "null n500",  md = "null", n = 500L), list(lab = "md40 n700",  md = "40",  n = 700L))
stem_of <- function(cl) if (cl$md == "null") sprintf("fs_maxeffCons_mr_field_mdnull_knoise0_n%d_%s", cl$n, tag) else sprintf("fs_maxeffCons_mr_field_md%s_knoise0_n%d_%s", cl$md, cl$n, tag)
path_of <- function(cl) { s <- stem_of(cl); f <- Sys.glob(file.path("mr_md_harm", paste0(s, "_d5000"), paste0(s, "_combined_*.rds"))); if (length(f)) f[1] else NA_character_ }
z95 <- qnorm(0.95)
wilson <- function(x, n, z = qnorm(0.975)) { if (!is.finite(x) || n <= 0) return(c(NA, NA)); ctr <- (x + z^2/(2*n))/(1 + z^2/n); hw <- z*sqrt(x*(1-x)/n + z^2/(4*n^2))/(1 + z^2/n); c(ctr - hw, ctr + hw) }
fmtw <- function(x, n) { w <- wilson(x, n); sprintf("%.3f (%.3f, %.3f)", x, w[1], w[2]) }
mfin <- function(x) { x <- x[is.finite(x)]; if (length(x)) mean(x) else NA_real_ }
sdfin <- function(x) { x <- x[is.finite(x)]; if (length(x) > 1) sd(x) else NA_real_ }
rows <- list(); bl <- list(H = list(), Hc = list()); jt <- list(); rg <- list(); disp <- list(); status <- list()
for (cl in cells) {
  p <- path_of(cl); if (is.na(p)) { status[[cl$lab]] <- "NOT RUN (deferred)"; next }
  b <- readRDS(p); r <- b$results; d <- r[r$detected %in% 1L, ]
  or <- -1; d$bH <- or * d$betaHhat_H; d$bHc <- or * d$betaHhat_Hc
  status[[cl$lab]] <- sprintf("%d reps, %d detected (%.3f), pkg %s", nrow(r), nrow(d), mean(r$detected), b$meta$pkg_version)
  struct <- list(H = or * b$truth$effect_Q, Hc = or * b$truth$effect_Qc)
  for (sfx in c("H", "Hc")) {
    side <- if (sfx == "H") "lower" else "upper"; sgn <- if (side == "lower") -1 else 1
    beta <- d[[if (sfx == "H") "bH" else "bHc"]]
    ests <- list(naive = list(e = d[[paste0("nv_", sfx, "_est")]], lo = d[[paste0("nv_", sfx, "_lo")]], hi = d[[paste0("nv_", sfx, "_hi")]], se = d[[paste0("nv_", sfx, "_se")]]),
                 oracle = list(e = d[[paste0("or_", sfx, "_est")]], lo = d[[paste0("or_", sfx, "_lo")]], hi = d[[paste0("or_", sfx, "_hi")]], se = d[[paste0("or_", sfx, "_se")]]),
                 `MR (IJ)` = list(e = d[[paste0("mr_", sfx, "_est")]], lo = d[[paste0("mr_", sfx, "_lo")]], hi = d[[paste0("mr_", sfx, "_hi")]], se = d[[paste0("mr_", sfx, "_se_ij")]]),
                 `MR (field)` = list(e = d[[paste0("fld_", sfx, "_est2")]], lo = d[[paste0("fld_", sfx, "_lo2s")]], hi = d[[paste0("fld_", sfx, "_hi2s")]], se = d[[paste0("fld_", sfx, "_se")]]))
    for (k in names(ests)) {
      ec <- ests[[k]]; e <- ec$e; se <- ec$se
      tgt <- if (k == "oracle") rep(struct[[sfx]], length(e)) else beta
      ctr <- if (k == "MR (field)") d[[paste0("mr_", sfx, "_est")]] else e
      b1 <- if (k == "MR (field)") d[[paste0("fld_", sfx, if (side == "upper") "_up1s" else "_lo1s")]] else e + sgn * z95 * se
      ok <- is.finite(e) & is.finite(tgt); ok2 <- ok & is.finite(ec$lo) & is.finite(ec$hi); ok1 <- ok & is.finite(b1)
      bias <- mean(e[ok] - tgt[ok]); sde <- sdfin(e); sem <- mfin(se)
      c2 <- mean(tgt[ok2] >= ec$lo[ok2] & tgt[ok2] <= ec$hi[ok2]); c1 <- if (side == "lower") mean(tgt[ok1] >= b1[ok1]) else mean(tgt[ok1] <= b1[ok1])
      rows[[length(rows) + 1]] <- data.frame(cell = cl$lab, block = sfx, estimator = k, n = sum(ok), bias = bias, bias_sd = bias / sde, SD = sde, SE = sem, SE_SD = sem / sde,
                                             cov2 = c2, cov2_w = fmtw(c2, sum(ok2)), cov1 = c1, cov1_w = fmtw(c1, sum(ok1)), side = side,
                                             halfwidth = mfin((ec$hi - ec$lo) / 2), margin1 = mfin(sgn * (b1 - ctr)), stringsAsFactors = FALSE)
      if (k != "oracle") {
        thr <- if (sfx == "H") c(0, 10, 20, 30, 40) else c(0, 10, 20, 30)
        bb <- b1[is.finite(b1)]
        sh <- vapply(thr, function(t) if (side == "lower") mean(bb >= t) else mean(bb <= t), numeric(1))
        bl[[sfx]][[length(bl[[sfx]]) + 1]] <- cbind(data.frame(cell = cl$lab, block = sfx, estimator = k, n = length(bb), mean = mean(bb), median = median(bb), q10 = unname(quantile(bb, .1)), q90 = unname(quantile(bb, .9)), stringsAsFactors = FALSE),
                                      as.data.frame(as.list(setNames(sh, sprintf(if (side == "lower") "P(L>=%g)" else "P(U<=%g)", thr))), check.names = FALSE))
      }
    }
  }
  # joint pair
  dj <- d[is.finite(d$fld_joint_gamma), ]
  pairs <- list(`separate 95% field bounds` = list(lo = dj$fld_H_lo1s, up = dj$fld_Hc_up1s), `Bonferroni (gamma = 0.025)` = list(lo = dj$fld_joint_bonf_loH, up = dj$fld_joint_bonf_upHc), `calibrated gamma` = list(lo = dj$fld_joint_loH, up = dj$fld_joint_upHc))
  for (nm in names(pairs)) { pp <- pairs[[nm]]; cH <- dj$bH >= pp$lo; cC <- dj$bHc <= pp$up
    jt[[length(jt) + 1]] <- data.frame(cell = cl$lab, pair = nm, n = nrow(dj), joint = fmtw(mean(cH & cC), nrow(dj)), cov_H = mean(cH), cov_Hc = mean(cC), margin_H = mean(dj$mr_H_est - pp$lo), margin_Hc = mean(pp$up - dj$mr_Hc_est), stringsAsFactors = FALSE) }
  # regime diagnostics
  rg[[length(rg) + 1]] <- data.frame(cell = cl$lab, n_det = nrow(d), mean_pH = mean(d$n_harm), p_hat_mean = mean(d$p_hat_H), p_hat_med = median(d$p_hat_H), p_hat_lt05 = mean(d$p_hat_H < 0.5), p_hat_ge09 = mean(d$p_hat_H >= 0.9),
                                     gamma_mean = mean(dj$fld_joint_gamma), corr_lam = mean(dj$fld_joint_corr), nfit_mean = mean(d$fld_Hc_nfit), share_newfit = mean(d$fld_Hc_share_newfit),
                                     sd_btc_naive = sdfin(d$mr_Hc_est) / mfin(d$nv_Hc_se), lamc_naive = mfin(d$fld_Hc_se) / mfin(d$nv_Hc_se), ij_sd_H = mfin(d$mr_H_se_ij) / sdfin(d$mr_H_est), ij_sd_Hc = mfin(d$mr_Hc_se_ij) / sdfin(d$mr_Hc_est),
                                     fit_secs = mean(r$fit_mr_secs), field_secs = mfin(d$fld_H_secs), comp_secs = mfin(d$fld_Hc_secs), stringsAsFactors = FALSE)
  # display rows (identity scale)
  ro <- r; ro$betaHhat_H <- or * r$betaHhat_H; ro$betaHhat_Hc <- or * r$betaHhat_Hc
  dh <- fs_sim_bias_coverage(ro, block = "H", estimators = c("mr", "fld"), scale = "identity"); dh$block <- "H"; dh$cell <- cl$lab
  dc <- fs_sim_bias_coverage(ro, block = "Hc", estimators = c("mr", "fld"), side = "upper", scale = "identity"); dc$block <- "Hc"; dc$cell <- cl$lab
  disp[[length(disp) + 1]] <- rbind(dh, dc)
}
T2 <- do.call(rbind, rows); BL_H <- do.call(rbind, bl$H); BL_Hc <- do.call(rbind, bl$Hc); JT <- do.call(rbind, jt); RG <- do.call(rbind, rg); DS <- do.call(rbind, disp)
saveRDS(list(T2 = T2, BL_H = BL_H, BL_Hc = BL_Hc, JT = JT, RG = RG, DS = DS, status = status), file.path(out_dir, "stage3_tables.rds"))
md <- function(df, digits = 3) { df2 <- df; for (k in names(df2)) if (is.numeric(df2[[k]])) df2[[k]] <- formatC(df2[[k]], digits = digits, format = "f"); paste(c(paste0("| ", paste(names(df2), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(df2)), collapse = "|"), "|"), apply(df2, 1, function(x) paste0("| ", paste(x, collapse = " | "), " |"))), collapse = "\n") }
sink(file.path(out_dir, "stage3_tables.md"))
cat("## status\n\n"); for (n in names(status)) cat(sprintf("- %s: %s\n", n, status[[n]]))
cat("\n## Table-2 layout (both blocks)\n\n", md(T2[, c("cell","block","estimator","n","bias","bias_sd","SD","SE","SE_SD","cov2_w","cov1_w","side","halfwidth","margin1")], 3), "\n")
cat("\n## Compact\n\n"); cmp <- data.frame(cell = T2$cell, block = T2$block, estimator = T2$estimator, `bias (MD | SD)` = sprintf("%+.2f | %+.3f", T2$bias, T2$bias_sd), SD = sprintf("%.2f", T2$SD), SE = sprintf("%.2f", T2$SE), `SE/SD` = sprintf("%.3f", T2$SE_SD), `two-sided` = T2$cov2_w, `one-sided` = paste(T2$cov1_w, T2$side), check.names = FALSE); cat(md(cmp), "\n")
cat("\n## Bound location, harm block (one-sided 95% LOWER bound)\n\n", md(BL_H, 3), "\n")
cat("\n## Bound location, complement block (one-sided 95% UPPER bound)\n\n", md(BL_Hc, 3), "\n")
cat("\n## Joint pair\n\n", md(JT, 3), "\n")
cat("\n## Regime diagnostics\n\n", md(RG, 3), "\n")
cat("\n## Display (identity scale)\n\n", md(DS[, c("cell","block","estimator","n","bias_log","sd_emp","se_mean","b","r","cov1","cov1_ref","cov2","cov2_ref")], 3), "\n")
sink()
# display figure: harm block (lower) and complement (upper), all cells
DS$cell <- factor(DS$cell)
pH <- fs_plot_bias_coverage(DS[DS$block == "H", ], side = "lower")
pC <- fs_plot_bias_coverage(DS[DS$block == "Hc", ], side = "upper")
ggsave(file.path(out_dir, "bias_coverage_display_H.png"), pH, width = 13, height = 4.6, dpi = 130)
ggsave(file.path(out_dir, "bias_coverage_display_Hc.png"), pC, width = 13, height = 4.6, dpi = 130)
cat("wrote", file.path(out_dir, "stage3_tables.md"), "\n")
