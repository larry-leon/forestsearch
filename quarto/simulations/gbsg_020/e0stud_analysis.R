# E0 analysis (TASK_field_studentize_stage1_e0_2026-09-08): tables 1-4 on the
# e0stud bundle (nb20-A HR 1.75 n500 configuration + field_decompose, sim_id
# 1-200).  Every number printed verbatim into the report.  Tertiles within the
# detected replicates; Wilson intervals for shares; SEs for means.
# Usage: Rscript e0_analysis.R <e0stud.rds> <out.md>
args <- commandArgs(trailingOnly = TRUE)
b <- readRDS(args[1]); out <- args[2]
r <- b$results
d <- r[r$detected %in% 1L & is.finite(r$fld_Hc_scale_ratio) & is.finite(r$fld_Hc_se) &
         is.finite(r$nv_Hc_se) & is.finite(r$p_hat_H), , drop = FALSE]
n_det <- sum(r$detected %in% 1L)
d$rho <- d$fld_Hc_scale_ratio
d$ratio_H <- d$n_harm / d$n_true
d$fld_over_nv <- d$fld_Hc_se / d$nv_Hc_se
d$corr_over_nv <- d$rho * d$fld_over_nv
wilson <- function(x, n, z = qnorm(0.975)) {
  p <- x / n; den <- 1 + z^2 / n; c0 <- (p + z^2 / (2 * n)) / den
  h <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den; c(lo = c0 - h, hi = c0 + h)
}
tert <- function(v) {
  br <- quantile(v, c(0, 1/3, 2/3, 1), names = FALSE)
  list(g = cut(v, breaks = unique(br), include.lowest = TRUE, labels = FALSE), br = br)
}
f3 <- function(x) formatC(x, digits = 3, format = "f")
f2 <- function(x) formatC(x, digits = 2, format = "f")
sink(out)
cat(sprintf("Detected replicates: %d of %d rows; analysed (every input finite): %d\n\n", n_det, nrow(r), nrow(d)))
cat(sprintf("Range of sim_id: %d-%d. Sanity: n_true range %d-%d; n_harm range %d-%d; p_hat_H range %.3f-%.3f.\n\n",
            min(d$sim_id), max(d$sim_id), min(d$n_true), max(d$n_true), min(d$n_harm), max(d$n_harm), min(d$p_hat_H), max(d$p_hat_H)))

# ---- Table 1: correlation and mean rho by tertile --------------------------
cat("## Table 1 -- corr(rho^c, fld_Hc_se / nv_Hc_se) and mean rho^c by tertile\n\n")
pe <- cor.test(d$rho, d$fld_over_nv, method = "pearson")
sp <- suppressWarnings(cor.test(d$rho, d$fld_over_nv, method = "spearman"))
cat(sprintf("| statistic | value | 95%% CI / p |\n|---|---|---|\n"))
cat(sprintf("| Pearson corr(rho^c, fld_Hc_se/nv_Hc_se) | %s | [%s, %s] |\n", f3(pe$estimate), f3(pe$conf.int[1]), f3(pe$conf.int[2])))
cat(sprintf("| Spearman corr(rho^c, fld_Hc_se/nv_Hc_se) | %s | p = %.2g |\n", f3(sp$estimate), sp$p.value))
cat(sprintf("| Pearson corr(rho^c, p_hat_H) | %s | |\n", f3(cor(d$rho, d$p_hat_H))))
cat(sprintf("| Spearman corr(rho^c, p_hat_H) | %s | |\n", f3(cor(d$rho, d$p_hat_H, method = "spearman"))))
cat(sprintf("| Pearson corr(rho^c, \\|Hhat\\|/\\|H\\|) | %s | |\n", f3(cor(d$rho, d$ratio_H))))
cat(sprintf("| Spearman corr(rho^c, \\|Hhat\\|/\\|H\\|) | %s | |\n\n", f3(cor(d$rho, d$ratio_H, method = "spearman"))))

strat_tab <- function(v, what) {
  tt <- tert(v); g <- tt$g
  cat(sprintf("| %s tertile | n | range | mean rho^c (SE) | median rho^c | share rho^c > 1 [Wilson] | mean fld/nv (SE) | mean rho^c x fld/nv (SE) | sqrt(mean lam^2)/sqrt(mean nSE^2) | corrected RMS form |\n", what))
  cat("|---|---|---|---|---|---|---|---|---|---|\n")
  for (k in sort(unique(g))) {
    dk <- d[g == k, ]; n <- nrow(dk); x <- sum(dk$rho > 1); w <- wilson(x, n)
    cat(sprintf("| T%d | %d | [%s, %s] | %s (%s) | %s | %s [%s, %s] | %s (%s) | %s (%s) | %s | %s |\n",
                k, n, f3(min(v[g == k])), f3(max(v[g == k])),
                f3(mean(dk$rho)), f3(sd(dk$rho) / sqrt(n)), f3(median(dk$rho)),
                f2(x / n), f2(w[1]), f2(w[2]),
                f3(mean(dk$fld_over_nv)), f3(sd(dk$fld_over_nv) / sqrt(n)),
                f3(mean(dk$corr_over_nv)), f3(sd(dk$corr_over_nv) / sqrt(n)),
                f3(sqrt(mean(dk$fld_Hc_se^2) / mean(dk$nv_Hc_se^2))),
                f3(sqrt(mean((dk$rho * dk$fld_Hc_se)^2) / mean(dk$nv_Hc_se^2)))))
  }
  cat("\n")
}
cat("**By p-hat(Hhat) tertile:**\n\n"); strat_tab(d$p_hat_H, "p-hat")
cat("**By |Hhat|/|H| tertile:**\n\n"); strat_tab(d$ratio_H, "|Hhat|/|H|")

# ---- Table 2: corrected vs uncorrected ratio by tertile (compact) ----------
cat("## Table 2 -- (rho^c x fld_Hc_se)/nv_Hc_se beside fld_Hc_se/nv_Hc_se, by tertile\n\n")
comp <- function(v, what) {
  tt <- tert(v); g <- tt$g
  cat(sprintf("| %s tertile | n | uncorrected mean fld/nv (SE) | corrected mean rho x fld/nv (SE) | uncorrected RMS | corrected RMS | mean nv_Hc_se | mean fld_Hc_se | mean rho x fld_Hc_se |\n", what))
  cat("|---|---|---|---|---|---|---|---|---|\n")
  for (k in sort(unique(g))) {
    dk <- d[g == k, ]; n <- nrow(dk)
    cat(sprintf("| T%d | %d | %s (%s) | %s (%s) | %s | %s | %s | %s | %s |\n", k, n,
                f3(mean(dk$fld_over_nv)), f3(sd(dk$fld_over_nv) / sqrt(n)),
                f3(mean(dk$corr_over_nv)), f3(sd(dk$corr_over_nv) / sqrt(n)),
                f3(sqrt(mean(dk$fld_Hc_se^2) / mean(dk$nv_Hc_se^2))),
                f3(sqrt(mean((dk$rho * dk$fld_Hc_se)^2) / mean(dk$nv_Hc_se^2))),
                f3(mean(dk$nv_Hc_se)), f3(mean(dk$fld_Hc_se)), f3(mean(dk$rho * dk$fld_Hc_se))))
  }
  n <- nrow(d)
  cat(sprintf("| all | %d | %s (%s) | %s (%s) | %s | %s | %s | %s | %s |\n\n", n,
              f3(mean(d$fld_over_nv)), f3(sd(d$fld_over_nv) / sqrt(n)),
              f3(mean(d$corr_over_nv)), f3(sd(d$corr_over_nv) / sqrt(n)),
              f3(sqrt(mean(d$fld_Hc_se^2) / mean(d$nv_Hc_se^2))),
              f3(sqrt(mean((d$rho * d$fld_Hc_se)^2) / mean(d$nv_Hc_se^2))),
              f3(mean(d$nv_Hc_se)), f3(mean(d$fld_Hc_se)), f3(mean(d$rho * d$fld_Hc_se))))
}
cat("**By |Hhat|/|H| tertile (A2's stratification):**\n\n"); comp(d$ratio_H, "|Hhat|/|H|")
cat("**By p-hat(Hhat) tertile:**\n\n"); comp(d$p_hat_H, "p-hat")
# per-replicate regression of the ratios on each other
fit_u <- lm(fld_over_nv ~ ratio_H, data = d); fit_c <- lm(corr_over_nv ~ ratio_H, data = d)
cat(sprintf("Slope of fld/nv on |Hhat|/|H|: %s (SE %s); slope of rho x fld/nv on |Hhat|/|H|: %s (SE %s).\n",
            f3(coef(fit_u)[2]), f3(sqrt(vcov(fit_u)[2, 2])), f3(coef(fit_c)[2]), f3(sqrt(vcov(fit_c)[2, 2]))))
cat(sprintf("SD across replicates: fld_Hc_se %s, rho x fld_Hc_se %s, nv_Hc_se %s; corr(fld_Hc_se, nv_Hc_se) %s, corr(rho x fld_Hc_se, nv_Hc_se) %s.\n\n",
            f3(sd(d$fld_Hc_se)), f3(sd(d$rho * d$fld_Hc_se)), f3(sd(d$nv_Hc_se)),
            f3(cor(d$fld_Hc_se, d$nv_Hc_se)), f3(cor(d$rho * d$fld_Hc_se, d$nv_Hc_se))))

# ---- Table 3: per-draw stability of s_G ------------------------------------
cat("## Table 3 -- fld_Hc_scale_cv (CV of the outer winners' complement scales, per replicate)\n\n")
q <- quantile(d$fld_Hc_scale_cv, c(.10, .25, .50, .75, .90, .99), names = FALSE)
cat("| mean | q10 | q25 | q50 | q75 | q90 | q99 | max |\n|---|---|---|---|---|---|---|---|\n")
cat(sprintf("| %s | %s | %s | %s | %s | %s | %s | %s |\n\n", f3(mean(d$fld_Hc_scale_cv)), f3(q[1]), f3(q[2]), f3(q[3]), f3(q[4]), f3(q[5]), f3(q[6]), f3(max(d$fld_Hc_scale_cv))))
tt <- tert(d$p_hat_H)
cat("| p-hat tertile | mean CV | q50 CV | q90 CV | mean scale_sel | mean scale_win |\n|---|---|---|---|---|---|\n")
for (k in sort(unique(tt$g))) { dk <- d[tt$g == k, ]
  cat(sprintf("| T%d | %s | %s | %s | %s | %s |\n", k, f3(mean(dk$fld_Hc_scale_cv)), f3(median(dk$fld_Hc_scale_cv)),
              f3(quantile(dk$fld_Hc_scale_cv, .9)), f3(mean(dk$fld_Hc_scale_sel)), f3(mean(dk$fld_Hc_scale_win)))) }
cat("\n")
cat(sprintf("Cross-check: scale_sel vs nv_Hc_se -- corr %s; mean scale_sel / mean nv_Hc_se = %s (scale_sel is the influence-norm sqrt(sum dfbeta^2), the sandwich-variance building block).\n\n",
            f3(cor(d$fld_Hc_scale_sel, d$nv_Hc_se)), f3(mean(d$fld_Hc_scale_sel) / mean(d$nv_Hc_se))))

# ---- Table 4: context row ---------------------------------------------------
cat("## Table 4 -- context: rho^c summary against the cell's committed average deficit\n\n")
n <- nrow(d); x <- sum(d$rho > 1); w <- wilson(x, n)
lam2 <- mean(d$fld_Hc_se^2); nse2 <- mean(d$nv_Hc_se^2)
cat("| quantity | value |\n|---|---|\n")
cat(sprintf("| mean rho^c (SE) | %s (%s) |\n", f3(mean(d$rho)), f3(sd(d$rho) / sqrt(n))))
cat(sprintf("| median rho^c | %s |\n", f3(median(d$rho))))
cat(sprintf("| q10 / q90 rho^c | %s / %s |\n", f3(quantile(d$rho, .1)), f3(quantile(d$rho, .9))))
cat(sprintf("| share rho^c > 1 [Wilson] | %s [%s, %s] (%d of %d) |\n", f2(x / n), f2(w[1]), f2(w[2]), x, n))
cat(sprintf("| these 200: lambda^2 / naive SE^2 (mean of squares) | %s |\n", f3(lam2 / nse2)))
cat(sprintf("| these 200: lambda / naive SE = sqrt of that | %s |\n", f3(sqrt(lam2 / nse2))))
cat(sprintf("| these 200: implied mean rho^c = naive SE / lambda | %s |\n", f3(sqrt(nse2 / lam2))))
cat(sprintf("| these 200: corrected (rho x lambda)^2 / naive SE^2 | %s |\n", f3(mean((d$rho * d$fld_Hc_se)^2) / nse2)))
cat(sprintf("| these 200: mean per-replicate fld/nv | %s |\n", f3(mean(d$fld_over_nv))))
cat(sprintf("| these 200: mean per-replicate rho x fld/nv | %s |\n", f3(mean(d$corr_over_nv))))
cat("| committed cell (A1, nb20-A HR 1.75 n500, n = 1999): lambda^2 / naive SE^2 | 0.894 |\n")
cat("| committed cell: lambda / naive SE = sqrt(0.894) | 0.945 |\n")
cat("| committed cell: implied mean rho^c = 1 / 0.945 | 1.058 |\n\n")
sink()
cat("written", out, "\n")
