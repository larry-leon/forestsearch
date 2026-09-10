# p12ext findings (TASK_p12ext_2026-09-09, Stage 3): the discriminator table,
# the mechanism, the nulls, the one-sided products and the recovery
# diagnostics, as markdown, verbatim into REPORT_p12ext_2026-09-09.md.
# Definitions follow summary_p12ext.qmd (this directory) and, for the
# fixed p-hat bands, REPORT_fixedphat_ij2s_2026-09-09.md Part G.
# Usage: Rscript p12ext_findings.R <out.md>
suppressPackageStartupMessages(library(forestsearch))
out <- commandArgs(trailingOnly = TRUE)[1]
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a

f3 <- function(x) formatC(x, digits = 3, format = "f")
f4 <- function(x) formatC(x, digits = 4, format = "f")
wil <- function(p, n, z = qnorm(.975)) { if (!is.finite(p) || n == 0) return(c(NA, NA))
  c <- (p + z^2/(2*n))/(1+z^2/n); h <- z*sqrt(p*(1-p)/n + z^2/(4*n^2))/(1+z^2/n); c(c - h, c + h) }
wfmt <- function(p, n) { w <- wil(p, n); sprintf("%s [%s, %s]", f3(p), f3(w[1]), f3(w[2])) }
tert <- function(v) { br <- quantile(v, c(0, 1/3, 2/3, 1), names = FALSE)
  cut(v, breaks = unique(br), include.lowest = TRUE, labels = FALSE) }
# Fixed p-hat bands, exactly as REPORT_fixedphat_ij2s_2026-09-09 Part G
PB <- c(0, 0.05, 0.10, 0.20, 0.35, 0.55, 1.0)
PBLAB <- c("[0, 0.05)", "[0.05, 0.10)", "[0.10, 0.20)", "[0.20, 0.35)", "[0.35, 0.55)", "[0.55, 1.0]")
pband <- function(v) cut(v, breaks = PB, include.lowest = TRUE, right = FALSE, labels = PBLAB)

# ---- the nine cells of the extended grid ------------------------------------
G <- list(
  list(key="h150_n500",  hr=1.50, n=500,  camp="p12ext", kind="harm"),
  list(key="h150_n1000", hr=1.50, n=1000, camp="p12ext", kind="harm"),
  list(key="h150_n1500", hr=1.50, n=1500, camp="p12ext", kind="harm"),
  list(key="h175_n500",  hr=1.75, n=500,  camp="tier2",  kind="harm"),
  list(key="h175_n1000", hr=1.75, n=1000, camp="tier2",  kind="harm"),
  list(key="h175_n1500", hr=1.75, n=1500, camp="tier2",  kind="harm"),
  list(key="h100_n500",  hr=1.00, n=500,  camp="tier2",  kind="null"),
  list(key="h100_n1000", hr=1.00, n=1000, camp="p12ext", kind="null"),
  list(key="h100_n1500", hr=1.00, n=1500, camp="p12ext", kind="null"))
gpath <- function(g) sprintf("results/fs_maxeffCons_fb_mr_field_m1_h%03d_knoise0_n%d_%s_combined_1_2000.rds",
                             round(100*g$hr), g$n, g$camp)
lab <- function(g) sprintf("HR %.2f, n = %d%s [%s]", g$hr, g$n,
                           if (identical(g$kind,"null")) " (null)" else "", g$camp)
for (i in seq_along(G)) { G[[i]]$path <- gpath(G[[i]]); G[[i]]$label <- lab(G[[i]])
  G[[i]]$have <- file.exists(G[[i]]$path) }
B <- lapply(G, function(g) if (g$have) readRDS(g$path) else NULL)
names(B) <- vapply(G, `[[`, "", "key")

# ---- the 31% comparator rows (cert20; a DIFFERENT identifier and band) ------
G31 <- list(
  list(hr=1.50, n=500), list(hr=1.50, n=1000), list(hr=1.50, n=1500),
  list(hr=1.75, n=500), list(hr=1.75, n=1000), list(hr=1.75, n=1500),
  list(hr=1.00, n=500), list(hr=1.00, n=1000), list(hr=1.00, n=1500))
# The 31% profile at eps 0.20 was produced by two campaigns: `cert20` covers
# every cell except HR 1.50 / 1.75 at n = 500, which come from `p30sgnb20`
# (the same effMaxSG eps 0.20 J 10 configuration, arm A).  Both are carried and
# the supplying campaign is named in the table.
for (i in seq_along(G31)) {
  cand <- sprintf("results/fs_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d_z1q60_nb20_%s_combined_1_2000.rds",
                  round(100*G31[[i]]$hr), G31[[i]]$n, c("cert20", "p30sgnb20"))
  hit <- cand[file.exists(cand)]
  G31[[i]]$path <- if (length(hit)) hit[1] else cand[1]
  G31[[i]]$camp <- if (length(hit)) sub("^.*z1q60_nb20_(.*)_combined.*$", "\\1", hit[1]) else NA_character_
  G31[[i]]$have <- length(hit) > 0 }

det <- function(b) { r <- b$results; r[r$detected %in% 1L, ] }
# IJ two-term two-sided coverage of beta(Hhat) on the harm block
ij2s <- function(d) mean(d$betaHhat_H >= d$mr_H_lo & d$betaHhat_H <= d$mr_H_hi)
# miss split: interval entirely ABOVE the target (lo > target) vs entirely BELOW (hi < target)
miss_above <- function(d) mean(d$mr_H_lo  > d$betaHhat_H)
miss_below <- function(d) mean(d$mr_H_hi  < d$betaHhat_H)

# ---- definitions transplanted verbatim from summary_p12ext.qmd --------------
subst_s <- function(r) { for (s in c("est2", "up1s", "lo1s", "lo2s", "hi2s", "lo_se", "hi_se", "se", "lam_mean")) r[[paste0("fld_Hc_", s)]] <- r[[paste0("fld_Hc_", s, "_s")]]; r }
err_sd <- function(r, block, est) { d <- r[r$detected %in% 1L, ]
  e <- switch(est, naive = d[[paste0("nv_", block, "_est")]], mr = d[[paste0("mr_", block, "_est")]], fld = d[[paste0("fld_", block, "_est2")]])
  ok <- is.finite(e) & is.finite(d[[paste0("betaHhat_", block)]]); sd(log(e[ok]) - log(d[[paste0("betaHhat_", block)]][ok])) }
cov_block <- function(r, block, side) {
  t1 <- fs_sim_bias_coverage(r, block = block, estimators = c("naive", "mr", "fld"), side = side)
  t1$sd_err <- vapply(as.character(t1$estimator), function(e) err_sd(r, block, e), 0)
  t1$construction <- c(naive = "naive", mr = "IJ two-term", fld = "field")[as.character(t1$estimator)]
  if (block == "Hc") { rs <- subst_s(r); t2 <- fs_sim_bias_coverage(rs, block = "Hc", estimators = "fld", side = side)
    t2$sd_err <- err_sd(rs, "Hc", "fld"); t2$construction <- "field-s"; t1 <- rbind(t1, t2) }
  t1$se_over_sd_err <- t1$se_mean / t1$sd_err
  t1 }

sink(out)

## ===========================================================================
cat("## 1. Standard tables\n\n")
cat("Nine cells: the five new `p12ext` cells and the four committed `tier2` cells read as they\n")
cat("stand. Rows naive / field / field-s / IJ two-term on both blocks; the harm block has no\n")
cat("field-s row (it is untouched by that construction). Winner-only and winner-floor excluded.\n")
cat("`bias` is on the log scale; `b` = bias / marginal SD and `bias/err SD` = bias / error SD are\n")
cat("the two SD-unit readings; `r` = mean SE / marginal SD; `SE/err SD` = mean SE / error SD.\n")
cat("One-sided coverage is on the **exposed** side (lower on the harm block, upper on the\n")
cat("complement); Wilson 95% limits on both the one-sided and the two-sided rates.\n\n")
COVALL <- list()
cat("| cell | block | construction | n | bias (log) | marginal SD | error SD | b | bias/err SD | mean SE | r | SE/err SD | one-sided [Wilson] | two-sided [Wilson] | Gaussian ref |\n")
cat("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n")
for (g in G) { if (!g$have) next
  r <- B[[g$key]]$results
  tb <- rbind(cbind(block = "Hhat (lower)",   cov_block(r, "H",  "lower")),
              cbind(block = "Hhat^c (upper)", cov_block(r, "Hc", "upper")))
  tb$construction <- factor(tb$construction, levels = c("naive", "field", "field-s", "IJ two-term"))
  tb <- tb[order(tb$block, tb$construction), ]
  tb$cell <- g$label
  COVALL[[g$key]] <- tb
  for (i in seq_len(nrow(tb))) { w2 <- wil(tb$cov2[i], tb$n[i])
    cat(sprintf("| %s | %s | %s | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s [%s, %s] | %s [%s, %s] | %s |\n",
      g$label, tb$block[i], as.character(tb$construction[i]), tb$n[i],
      f4(tb$bias_log[i]), f4(tb$sd_emp[i]), f4(tb$sd_err[i]),
      f3(tb$b[i]), f3(tb$bias_log[i]/tb$sd_err[i]), f4(tb$se_mean[i]), f3(tb$r[i]),
      f3(tb$se_over_sd_err[i]),
      f3(tb$cov1[i]), f3(tb$cov1_wilson_lo[i]), f3(tb$cov1_wilson_hi[i]),
      f3(tb$cov2[i]), f3(w2[1]), f3(w2[2]), f3(tb$cov1_ref[i]))) } }
cat("\n### 1a. Across cells\n\n")
CA <- do.call(rbind, COVALL)
cat("| block | construction | cells | one-sided min | mean | max | two-sided mean | mean b | mean r | mean SE/err SD |\n|---|---|---|---|---|---|---|---|---|---|\n")
for (bl in unique(CA$block)) for (co in levels(CA$construction)) {
  d <- CA[CA$block == bl & CA$construction == co, ]; if (!nrow(d)) next
  cat(sprintf("| %s | %s | %d | %s | %s | %s | %s | %s | %s | %s |\n", bl, co, nrow(d),
      f3(min(d$cov1)), f3(mean(d$cov1)), f3(max(d$cov1)), f3(mean(d$cov2)),
      f3(mean(d$b)), f3(mean(d$r)), f3(mean(d$se_over_sd_err)))) }
cat("\n")

## ===========================================================================
cat("## 2. The discriminator table\n\n")
cat("IJ two-term **two-sided** coverage of beta(Hhat) on the **harm block**, at the M1 default\n")
cat("prevalence (12.4%), focus `maxeffCons`, J = 10. Wilson 95% limits. The 31% rows beside are the\n")
cat("31% cells -- prevalence 0.3065 (`FS_S7_Z1Q=0.60`), focus `effMaxSG` at eps = 0.20, campaigns\n")
cat("`cert20` and (at n = 500) `p30sgnb20`, named in each entry -- a\n")
cat("different identifier and band, carried for contrast only.\n\n")
cat("| n | HR 1.50 (12.4%) [Wilson] | HR 1.75 (12.4%) [Wilson] | HR 1.50 minus HR 1.75 | HR 1.50 (31%) | HR 1.75 (31%) |\n")
cat("|---|---|---|---|---|---|\n")
DISC <- list()
for (nn in c(500, 1000, 1500)) {
  k150 <- sprintf("h150_n%d", nn); k175 <- sprintf("h175_n%d", nn)
  c150 <- if (!is.null(B[[k150]])) { d <- det(B[[k150]]); c(ij2s(d), nrow(d)) } else c(NA, 0)
  c175 <- if (!is.null(B[[k175]])) { d <- det(B[[k175]]); c(ij2s(d), nrow(d)) } else c(NA, 0)
  g31a <- Filter(function(z) z$hr == 1.50 && z$n == nn && z$have, G31)
  g31b <- Filter(function(z) z$hr == 1.75 && z$n == nn && z$have, G31)
  s31a <- if (length(g31a)) { d <- det(readRDS(g31a[[1]]$path)); sprintf("%s (%s)", wfmt(ij2s(d), nrow(d)), g31a[[1]]$camp) } else "-"
  s31b <- if (length(g31b)) { d <- det(readRDS(g31b[[1]]$path)); sprintf("%s (%s)", wfmt(ij2s(d), nrow(d)), g31b[[1]]$camp) } else "-"
  cat(sprintf("| %d | %s | %s | %s | %s | %s |\n", nn,
      if (is.na(c150[1])) "-" else wfmt(c150[1], c150[2]),
      if (is.na(c175[1])) "-" else wfmt(c175[1], c175[2]),
      if (is.na(c150[1]) || is.na(c175[1])) "-" else f3(c150[1] - c175[1]), s31a, s31b))
  DISC[[as.character(nn)]] <- list(c150 = c150, c175 = c175)
}
cat("\n")
cat("**The decay at each HR, as a trajectory in n** (12.4% prevalence, harm block, IJ two-term two-sided):\n\n")
for (h in c("h150", "h175")) {
  v <- vapply(c(500,1000,1500), function(nn) { k <- sprintf("%s_n%d", h, nn)
    if (is.null(B[[k]])) NA_real_ else ij2s(det(B[[k]])) }, 0)
  cat(sprintf("- HR %s: %s -> %s -> %s at n = 500 -> 1000 -> 1500; total change %s\n",
      sub("h", "", sub("h1", "1.", h)), f4(v[1]), f4(v[2]), f4(v[3]),
      if (all(is.finite(v))) f4(v[3] - v[1]) else "-"))
}
cat("\n")

## ===========================================================================
cat("## 3. The mechanism on the harm block\n\n")
cat("Bias is on the log scale, `mean(log(mr_H_est) - log(betaHhat_H))`; SE is `mean(mr_H_se_ij)`;\n")
cat("the two-sided miss is split by **side**: *above* = the interval lies entirely above the\n")
cat("target (`mr_H_lo > betaHhat_H`), *below* = entirely below it (`mr_H_hi < betaHhat_H`).\n\n")
cat("### 3a. Trajectory in n, per HR\n\n")
cat("| cell | n detected | bias (log) | mean SE | bias / SE | marginal SD | error SD | two-sided cov | miss above | miss below |\n")
cat("|---|---|---|---|---|---|---|---|---|---|\n")
MECH <- list()
for (g in G) { if (!g$have || g$kind == "null") next
  d <- det(B[[g$key]]); ok <- is.finite(d$mr_H_est) & is.finite(d$betaHhat_H)
  e <- log(d$mr_H_est[ok]) - log(d$betaHhat_H[ok])
  bi <- mean(e); se <- mean(d$mr_H_se_ij[ok]); sd_m <- sd(log(d$mr_H_est[ok])); sd_e <- sd(e)
  cv <- ij2s(d); ma <- miss_above(d); mb <- miss_below(d)
  cat(sprintf("| %s | %d | %s | %s | %s | %s | %s | %s | %s | %s |\n", g$label, nrow(d),
      f4(bi), f4(se), f3(bi/se), f4(sd_m), f4(sd_e), wfmt(cv, nrow(d)), f4(ma), f4(mb)))
  MECH[[g$key]] <- c(bias = bi, se = se, ratio = bi/se, cov = cv, above = ma, below = mb)
}
cat("\n### 3b. The same quantities by fixed p-hat band\n\n")
cat("Bands are fixed (not tertiles), identical across cells, so the same band is the same\n")
cat("re-selection regime everywhere. Bands with fewer than 30 replicates are flagged.\n\n")
cat("| cell | p-hat band | n | flag | bias (log) | mean SE | bias / SE | two-sided cov [Wilson] | miss above | miss below |\n")
cat("|---|---|---|---|---|---|---|---|---|---|\n")
for (g in G) { if (!g$have || g$kind == "null") next
  d <- det(B[[g$key]]); d <- d[is.finite(d$p_hat_H), ]; bnd <- pband(d$p_hat_H)
  for (bl in PBLAB) { dk <- d[!is.na(bnd) & bnd == bl, , drop = FALSE]
    if (nrow(dk) == 0) { cat(sprintf("| %s | %s | 0 | EMPTY | - | - | - | - | - | - |\n", g$label, bl)); next }
    ok <- is.finite(dk$mr_H_est) & is.finite(dk$betaHhat_H)
    e <- log(dk$mr_H_est[ok]) - log(dk$betaHhat_H[ok]); bi <- mean(e); se <- mean(dk$mr_H_se_ij[ok])
    cat(sprintf("| %s | %s | %d | %s | %s | %s | %s | %s | %s | %s |\n", g$label, bl, nrow(dk),
        if (nrow(dk) < 30) "**n < 30**" else "", f4(bi), f4(se), f3(bi/se),
        wfmt(ij2s(dk), nrow(dk)), f4(miss_above(dk)), f4(miss_below(dk)))) } }
cat("\n")

## ===========================================================================
cat("## 4. The nulls (HR 1.00)\n\n")
cat("Detection rate is over all 2,000 replicates; every other quantity is among detections.\n")
cat("`reaches 1.00` is the share of detected replicates whose harm **lower** bound is >= 1.00.\n\n")
cat("| cell | detection [Wilson] | n detected | harm field lower one-sided [Wilson] | field-s complement upper one-sided [Wilson] | harm lower mean | share lower >= 1.00 | share lower >= 0.95 | share lower >= 0.85 | Hc upper mean | share upper < 0.85 | share upper < 0.80 |\n")
cat("|---|---|---|---|---|---|---|---|---|---|---|---|\n")
for (g in G) { if (!g$have || g$kind != "null") next
  r <- B[[g$key]]$results; dr <- mean(r$detected %in% 1L); d <- det(B[[g$key]])
  covH  <- mean(d$betaHhat_H  >= d$fld_H_lo1s)
  covHc <- mean(d$betaHhat_Hc <= d$fld_Hc_up1s_s)
  lo <- d$fld_H_lo1s; up <- d$fld_Hc_up1s_s
  cat(sprintf("| %s | %s | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n", g$label,
      wfmt(dr, nrow(r)), nrow(d), wfmt(covH, nrow(d)), wfmt(covHc, nrow(d)),
      f3(mean(lo)), f3(mean(lo >= 1.00)), f3(mean(lo >= 0.95)), f3(mean(lo >= 0.85)),
      f3(mean(up)), f3(mean(up < 0.85)), f3(mean(up < 0.80)))) }
cat("\n")

## ===========================================================================
cat("## 5. One-sided products across the extended grid\n\n")
cat("Harm block: the **field** one-sided lower bound on beta(Hhat) (`fld_H_lo1s`).\n")
cat("Complement block: the **field-s** one-sided upper bound on beta(Hhat^c) (`fld_Hc_up1s_s`).\n")
cat("Certified ranges carried as reference lines only: harm field lower **0.944-0.974**,\n")
cat("field-s complement upper **0.912-0.960**. No criterion is pre-registered for this campaign.\n\n")
cat("| cell | n | harm field lower [Wilson] | in 0.944-0.974 | field-s Hc upper [Wilson] | in 0.912-0.960 |\n")
cat("|---|---|---|---|---|---|\n")
for (g in G) { if (!g$have) next
  d <- det(B[[g$key]]); cH <- mean(d$betaHhat_H >= d$fld_H_lo1s); cC <- mean(d$betaHhat_Hc <= d$fld_Hc_up1s_s)
  cat(sprintf("| %s | %d | %s | %s | %s | %s |\n", g$label, nrow(d), wfmt(cH, nrow(d)),
      if (cH >= 0.944 && cH <= 0.974) "yes" else "no", wfmt(cC, nrow(d)),
      if (cC >= 0.912 && cC <= 0.960) "yes" else "no")) }
cat("\n### 5a. By p-hat tertile (within cell), both blocks\n\n")
cat("| cell | tertile | range | n | harm field lower [Wilson] | field-s Hc upper [Wilson] |\n|---|---|---|---|---|---|\n")
for (g in G) { if (!g$have) next
  d <- det(B[[g$key]]); d <- d[is.finite(d$p_hat_H), ]; tt <- tert(d$p_hat_H)
  for (t in sort(unique(tt))) { dk <- d[tt == t, ]
    cat(sprintf("| %s | T%d | [%s, %s] | %d | %s | %s |\n", g$label, t,
        f3(min(dk$p_hat_H)), f3(max(dk$p_hat_H)), nrow(dk),
        wfmt(mean(dk$betaHhat_H >= dk$fld_H_lo1s), nrow(dk)),
        wfmt(mean(dk$betaHhat_Hc <= dk$fld_Hc_up1s_s), nrow(dk)))) } }
cat("\n### 5b. By fixed p-hat band, both blocks\n\n")
cat("| cell | p-hat band | n | flag | harm field lower [Wilson] | field-s Hc upper [Wilson] |\n|---|---|---|---|---|---|\n")
for (g in G) { if (!g$have) next
  d <- det(B[[g$key]]); d <- d[is.finite(d$p_hat_H), ]; bnd <- pband(d$p_hat_H)
  for (bl in PBLAB) { dk <- d[!is.na(bnd) & bnd == bl, , drop = FALSE]
    if (nrow(dk) == 0) { cat(sprintf("| %s | %s | 0 | EMPTY | - | - |\n", g$label, bl)); next }
    cat(sprintf("| %s | %s | %d | %s | %s | %s |\n", g$label, bl, nrow(dk),
        if (nrow(dk) < 30) "**n < 30**" else "",
        wfmt(mean(dk$betaHhat_H >= dk$fld_H_lo1s), nrow(dk)),
        wfmt(mean(dk$betaHhat_Hc <= dk$fld_Hc_up1s_s), nrow(dk)))) } }
cat("\n")

## ===========================================================================
cat("## 6. Recovery diagnostics at campaign scale\n\n")
cat("Recorded only on the five `p12ext` cells (`field_recovery = TRUE`); the four `tier2` cells\n")
cat("predate the knob and carry no `fld_recov_*` columns. Descriptive: no construction reads them.\n")
cat("`sens_H` is the mean share of the identified patients that the outer draws' re-selections\n")
cat("retain; `q10/q50/q90` and `share1` describe the per-draw containment behind that mean.\n")
cat("`sens` (truth-referenced) is the recorder's sensitivity of Hhat against the planted region.\n\n")
cat("| cell | n | sens_H mean [sd] | q10 | q50 | q90 | ppv_H mean [sd] | sens_Hc mean | npv_Hc mean | share1 mean | n_used mean | all finite | in [0,1] |\n")
cat("|---|---|---|---|---|---|---|---|---|---|---|---|---|\n")
for (g in G) { if (!g$have) next
  d <- det(B[[g$key]]); if (!("fld_recov_sens_H" %in% names(d)) || all(is.na(d$fld_recov_sens_H))) next
  q <- quantile(d$fld_recov_sens_H, c(.10,.50,.90), na.rm = TRUE)
  fin <- all(is.finite(d$fld_recov_sens_H)) && all(is.finite(d$fld_recov_ppv_H)) &&
         all(is.finite(d$fld_recov_sens_Hc)) && all(is.finite(d$fld_recov_npv_Hc))
  inr <- all(d$fld_recov_sens_H >= 0 & d$fld_recov_sens_H <= 1 &
             d$fld_recov_ppv_H  >= 0 & d$fld_recov_ppv_H  <= 1, na.rm = TRUE)
  cat(sprintf("| %s | %d | %s [%s] | %s | %s | %s | %s [%s] | %s | %s | %s | %s | %s | %s |\n",
      g$label, nrow(d), f3(mean(d$fld_recov_sens_H)), f3(sd(d$fld_recov_sens_H)),
      f3(q[1]), f3(q[2]), f3(q[3]),
      f3(mean(d$fld_recov_ppv_H)), f3(sd(d$fld_recov_ppv_H)),
      f3(mean(d$fld_recov_sens_Hc)), f3(mean(d$fld_recov_npv_Hc)),
      f3(mean(d$fld_recov_share1)), f3(mean(d$fld_recov_n_used)), fin, inr)) }
cat("\n**Containment quantiles across replicates** (the per-replicate `q10`/`q50`/`q90` columns, summarised):\n\n")
cat("| cell | q10: mean [q25, q75] | q50: mean [q25, q75] | q90: mean [q25, q75] | share1 mean | share1 q90 |\n|---|---|---|---|---|---|\n")
for (g in G) { if (!g$have) next
  d <- det(B[[g$key]]); if (!("fld_recov_q10" %in% names(d)) || all(is.na(d$fld_recov_q10))) next
  s <- function(v) sprintf("%s [%s, %s]", f3(mean(v)), f3(quantile(v,.25)), f3(quantile(v,.75)))
  cat(sprintf("| %s | %s | %s | %s | %s | %s |\n", g$label,
      s(d$fld_recov_q10), s(d$fld_recov_q50), s(d$fld_recov_q90),
      f3(mean(d$fld_recov_share1)), f3(quantile(d$fld_recov_share1, .90)))) }
cat("\n**Correlations** (Pearson, on detected replicates):\n\n")
cat("| cell | corr(sens_H, p-hat) | corr(ppv_H, p-hat) | corr(sens_H, sens truth-ref) | corr(ppv_H, ppv truth-ref) | corr(sens_H, ppv_H) | corr(share1, p-hat) |\n|---|---|---|---|---|---|---|\n")
for (g in G) { if (!g$have) next
  d <- det(B[[g$key]]); if (!("fld_recov_sens_H" %in% names(d)) || all(is.na(d$fld_recov_sens_H))) next
  cc <- function(a, b) { ok <- is.finite(a) & is.finite(b); if (sum(ok) < 3) NA else cor(a[ok], b[ok]) }
  cat(sprintf("| %s | %s | %s | %s | %s | %s | %s |\n", g$label,
      f3(cc(d$fld_recov_sens_H, d$p_hat_H)), f3(cc(d$fld_recov_ppv_H, d$p_hat_H)),
      f3(cc(d$fld_recov_sens_H, d$sens)),    f3(cc(d$fld_recov_ppv_H, d$ppv)),
      f3(cc(d$fld_recov_sens_H, d$fld_recov_ppv_H)), f3(cc(d$fld_recov_share1, d$p_hat_H)))) }
cat("\n")

## ===========================================================================
cat("## Gate 2 record, per cell\n\n")
cat("| cell | rows | sim_id 1-2000 | dups | CONFIG-ERROR | detected | realized prevalence (mean n_true/n) | non-finite in any block | invariants | gamma in [0.025,0.05] | joint >= 0.95 - 2/n_joint | max bound<->quantile dev | p-hat valid | recovery finite / in [0,1] | n_workers | version |\n")
cat("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n")
for (g in G) { if (!g$have) next
  b <- B[[g$key]]; r <- b$results; m <- b$meta; d <- det(b)
  ids_ok <- identical(sort(as.integer(r$sim_id)), 1:2000)
  dups <- sum(duplicated(r$sim_id))
  cfg <- sum(grepl("CONFIG-ERROR", paste(r$status, r$err_msg), fixed = TRUE))
  blocks <- c("nv_H_est","nv_H_lo","nv_H_hi","nv_H_se","nv_Hc_est","nv_Hc_lo","nv_Hc_hi","nv_Hc_se",
    "fld_H_est2","fld_H_lo2s","fld_H_hi2s","fld_H_lo1s","fld_H_se","fld_H_lam_mean",
    "fld_Hc_est2","fld_Hc_up1s","fld_Hc_lo2s","fld_Hc_hi2s","fld_Hc_se","fld_Hc_lam_mean",
    "fld_Hc_est2_s","fld_Hc_up1s_s","fld_Hc_lo2s_s","fld_Hc_hi2s_s","fld_Hc_se_s",
    "fld_Hc_scale_ratio","fld_joint_gamma","fld_joint_prob","fld_joint_loH","fld_joint_upHc",
    "mr_H_est","mr_H_lo","mr_H_hi","mr_H_se_ij","mr_Hc_est","mr_Hc_lo","mr_Hc_hi","mr_Hc_se_ij",
    "betaHhat_H","betaHhat_Hc","p_hat_H")
  blocks <- intersect(blocks, names(d))
  nf <- sum(vapply(blocks, function(cn) sum(!is.finite(d[[cn]])), 0L))
  invs <- all(d$fld_H_lo2s <= d$fld_H_hi2s, na.rm=TRUE) && all(d$fld_Hc_lo2s <= d$fld_Hc_hi2s, na.rm=TRUE) &&
          all(d$fld_Hc_lo2s_s <= d$fld_Hc_hi2s_s, na.rm=TRUE) && all(d$mr_H_lo <= d$mr_H_hi, na.rm=TRUE) &&
          all(d$mr_Hc_lo <= d$mr_Hc_hi, na.rm=TRUE) && all(d$fld_H_lo1s >= d$fld_H_lo2s - 1e-12, na.rm=TRUE) &&
          all(d$fld_Hc_up1s <= d$fld_Hc_hi2s + 1e-12, na.rm=TRUE)
  gok <- all(d$fld_joint_gamma >= 0.025 & d$fld_joint_gamma <= 0.05, na.rm=TRUE)
  jok <- all(d$fld_joint_prob >= 0.95 - 2/d$fld_joint_n, na.rm=TRUE)
  bdH <- log(d$fld_H_est2) + d$fld_H_lam_mean; bdC <- log(d$fld_Hc_est2) + d$fld_Hc_lam_mean
  dev <- max(c(abs(log(d$fld_H_lo2s)  - (bdH - d$fld_H_q975)),
               abs(log(d$fld_H_hi2s)  - (bdH - d$fld_H_q025)),
               abs(log(d$fld_H_lo1s)  - (bdH - d$fld_H_q95)),
               abs(log(d$fld_Hc_lo2s) - (bdC - d$fld_Hc_q975)),
               abs(log(d$fld_Hc_hi2s) - (bdC - d$fld_Hc_q025)),
               abs(log(d$fld_Hc_up1s) - (bdC - d$fld_Hc_q05))), na.rm = TRUE)
  pok <- all(d$p_hat_H >= 0 & d$p_hat_H <= 1, na.rm=TRUE) && all(d$p_hat_sum <= 1 + 1e-8, na.rm=TRUE)
  hasrec <- "fld_recov_sens_H" %in% names(d) && !all(is.na(d$fld_recov_sens_H))
  rok <- if (!hasrec) "n/a (knob predates cell)" else
    sprintf("%s / %s", all(is.finite(d$fld_recov_sens_H)) && all(is.finite(d$fld_recov_ppv_H)),
            all(d$fld_recov_sens_H >= 0 & d$fld_recov_sens_H <= 1 &
                d$fld_recov_ppv_H  >= 0 & d$fld_recov_ppv_H  <= 1))
  nw <- m$n_workers
  if (is.null(nw)) { bf <- sub("_combined_1_2000", "_res_1_1000", g$path)
    nw <- if (file.exists(bf)) readRDS(bf)$meta$n_workers else NA }
  cat(sprintf("| %s | %d | %s | %d | %d | %d | %s | %d | %s | %s | %s | %s | %s | %s | %s | %s |\n",
      g$label, nrow(r), ids_ok, dups, cfg, nrow(d), f4(mean(r$n_true/m$n_sample)),
      nf, invs, gok, jok, format(dev, digits = 3, scientific = TRUE), pok, rok,
      as.character(nw %||% NA), as.character(m$forestsearch_version))) }
cat("\n**Meta knobs as set (per cell):**\n\n")
cat("| cell | focus | z1q | prev_super | J | nbhd | field_complement | field_decompose | field_scale_complement | field_recovery | ij_residual | fb | campaign | seed_base | n_batches |\n")
cat("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n")
for (g in G) { if (!g$have) next; m <- B[[g$key]]$meta
  cat(sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n", g$label,
      m$sg_focus, f3(m$harm_z1_quantile), f4(m$harm_prevalence_super), as.character(m$er_jcuts),
      f3(m$effect_neighborhood), as.character(m$field_complement), as.character(m$field_decompose),
      as.character(m$field_scale_complement), as.character(m$field_recovery %||% NA),
      as.character(m$ij_residual), paste(unique(m$fb_mode_by_batch), collapse="/"),
      as.character(m$campaign_tag), as.character(m$seed_base), as.character(m$n_batches))) }
cat("\n")
sink()
cat("wrote", out, "\n")
