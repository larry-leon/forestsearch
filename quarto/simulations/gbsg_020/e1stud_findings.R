# E1 findings (TASK_field_studentize_e1_2026-09-08, Stage 3): the standard
# tables and the five pre-registered findings as markdown, verbatim into
# REPORT_field_studentize_e1_2026-09-08.md.  Same definitions as
# summary_e1stud.qmd (this directory).  Usage: Rscript e1stud_findings.R <out.md>
suppressPackageStartupMessages(library(forestsearch))
out <- commandArgs(trailingOnly = TRUE)[1]
cells <- list(
  list(key = "eps20_h150", label = "effMaxSG eps 0.20, HR 1.50", kind = "band", comp = "nb20_p30sgnb20", e1 = "nb20_e1stud", focus = "effMaxSG", h = "h150"),
  list(key = "eps20_h175", label = "effMaxSG eps 0.20, HR 1.75", kind = "band", comp = "nb20_p30sgnb20", e1 = "nb20_e1stud", focus = "effMaxSG", h = "h175"),
  list(key = "eps30_h150", label = "effMaxSG eps 0.30 (stress), HR 1.50", kind = "band", comp = "nb30_banddial", e1 = "nb30_e1stud", focus = "effMaxSG", h = "h150"),
  list(key = "eps30_h175", label = "effMaxSG eps 0.30 (stress), HR 1.75", kind = "band", comp = "nb30_banddial", e1 = "nb30_e1stud", focus = "effMaxSG", h = "h175"),
  list(key = "maxSG_h175", label = "maxSG, HR 1.75", kind = "end", comp = "banddial", e1 = "e1stud", focus = "maxSG", h = "h175"),
  list(key = "minSG_h175", label = "minSG, HR 1.75", kind = "end", comp = "banddial", e1 = "e1stud", focus = "minSG", h = "h175"))
path <- function(cl, which) sprintf("results/fs_%s_fb_mr_field_m1_%s_knoise0_n500_z1q60_%s_combined_1_2000.rds", cl$focus, cl$h, cl[[which]])
B <- lapply(cells, function(cl) if (file.exists(path(cl, "e1"))) list(e1 = readRDS(path(cl, "e1")), comp = readRDS(path(cl, "comp")), cl = cl) else NULL)
B <- B[!vapply(B, is.null, TRUE)]
f3 <- function(x) formatC(x, digits = 3, format = "f"); f4 <- function(x) formatC(x, digits = 4, format = "f")
wil <- function(p, n, z = qnorm(.975)) { c <- (p + z^2/(2*n))/(1+z^2/n); h <- z*sqrt(p*(1-p)/n + z^2/(4*n^2))/(1+z^2/n); c(c - h, c + h) }
wfmt <- function(p, n) { w <- wil(p, n); sprintf("%s [%s, %s]", f3(p), f3(w[1]), f3(w[2])) }
subst_s <- function(r) { for (s in c("est2", "up1s", "lo1s", "lo2s", "hi2s", "lo_se", "hi_se", "se", "lam_mean")) r[[paste0("fld_Hc_", s)]] <- r[[paste0("fld_Hc_", s, "_s")]]; r }
err_sd <- function(r, block, est) { d <- r[r$detected %in% 1L, ]
  e <- switch(est, naive = d[[paste0("nv_", block, "_est")]], mr = d[[paste0("mr_", block, "_est")]], fld = d[[paste0("fld_", block, "_est2")]])
  ok <- is.finite(e) & is.finite(d[[paste0("betaHhat_", block)]]); sd(log(e[ok]) - log(d[[paste0("betaHhat_", block)]][ok])) }
tert <- function(v) { br <- quantile(v, c(0, 1/3, 2/3, 1), names = FALSE); cut(v, breaks = unique(br), include.lowest = TRUE, labels = FALSE) }
sink(out)
# ---- Gate 2 summary per cell -------------------------------------------------
cat("## Gate 2 per cell (pairing identity to the committed comparator)\n\n")
cat("| cell | rows | detected | pre-existing non-timing columns identical | truth identical | _s finite on detected | mean rho^c | joint_s at the alpha/2 fallback (joint) |\n|---|---|---|---|---|---|---|---|\n")
timing_cols <- c("fb_secs", "fit_mr_secs", "fld_H_secs", "fld_Hc_secs", "fld_H_uniform_secs")
for (b in B) { rc <- b$comp$results; rf <- b$e1$results; rc <- rc[order(rc$sim_id), ]; rf <- rf[order(rf$sim_id), ]
  pre <- setdiff(names(rc), timing_cols); id <- vapply(pre, function(cn) identical(rc[[cn]], rf[[cn]]), TRUE); det <- rf$detected %in% 1L
  fin <- all(is.finite(rf$fld_Hc_up1s_s[det])) && all(is.finite(rf$fld_joint_s_gamma[det]))
  cat(sprintf("| %s | %d | %d | %d / %d | %s | %s | %s | %s (%s) |\n", b$cl$label, nrow(rf), sum(det), sum(id), length(pre), identical(b$comp$truth, b$e1$truth), fin,
              f3(mean(rf$fld_Hc_scale_ratio[det])), f3(mean(abs(rf$fld_joint_s_gamma[det] - 0.025) < 1e-9)), f3(mean(abs(rf$fld_joint_gamma[det] - 0.025) < 1e-9)))) }
cat("\n")
# ---- Constructions table ------------------------------------------------------
cat("## Constructions per cell (identical replicates): naive / field / field-s / IJ two-term\n\n")
cat("| cell | block | construction | n | bias (log) | marginal SD | error SD | mean SE | b | r = SE/marg SD | SE/error SD | one-sided cov [Wilson] | two-sided cov | Gaussian ref |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n")
COV <- list()
for (b in B) { r <- b$e1$results
  for (bl in c("H", "Hc")) { side <- if (bl == "H") "lower" else "upper"
    t1 <- fs_sim_bias_coverage(r, block = bl, estimators = c("naive", "mr", "fld"), side = side); t1$sd_err <- vapply(as.character(t1$estimator), function(e) err_sd(r, bl, e), 0)
    t1$construction <- c(naive = "naive", mr = "IJ two-term", fld = "field")[as.character(t1$estimator)]
    if (bl == "Hc") { rs <- subst_s(r); t2 <- fs_sim_bias_coverage(rs, block = "Hc", estimators = "fld", side = side); t2$sd_err <- err_sd(rs, "Hc", "fld"); t2$construction <- "field-s"; t1 <- rbind(t1, t2) }
    t1 <- t1[order(match(t1$construction, c("naive", "field", "field-s", "IJ two-term"))), ]
    for (i in seq_len(nrow(t1))) cat(sprintf("| %s | %s | %s | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n", b$cl$label, if (bl == "H") "Hhat (lower)" else "Hhat^c (upper)", t1$construction[i], t1$n[i],
      f3(t1$bias_log[i]), f3(t1$sd_emp[i]), f3(t1$sd_err[i]), f3(t1$se_mean[i]), f3(t1$b[i]), f3(t1$r[i]), f3(t1$se_mean[i] / t1$sd_err[i]), wfmt(t1$cov1[i], t1$n[i]), f3(t1$cov2[i]), f3(t1$cov1_ref[i])))
    COV[[length(COV) + 1]] <- cbind(cell = b$cl$label, block = bl, t1) } }
COV <- do.call(rbind, COV); cat("\n")
cat("**Across cells** (one-sided coverage: min / mean / max; mean two-sided; mean b; mean r; mean SE/error SD):\n\n| block | construction | cells | cov1 min | cov1 mean | cov1 max | cov2 mean | b mean | r mean | SE/error SD mean |\n|---|---|---|---|---|---|---|---|---|---|\n")
for (bl in c("H", "Hc")) for (cn in c("naive", "field", "field-s", "IJ two-term")) { d <- COV[COV$block == bl & COV$construction == cn, ]; if (!nrow(d)) next
  cat(sprintf("| %s | %s | %d | %s | %s | %s | %s | %s | %s | %s |\n", if (bl == "H") "Hhat (lower)" else "Hhat^c (upper)", cn, nrow(d), f3(min(d$cov1)), f3(mean(d$cov1)), f3(max(d$cov1)), f3(mean(d$cov2)), f3(mean(d$b)), f3(mean(d$r)), f3(mean(d$se_mean / d$sd_err)))) }
cat("\n")
# ---- Finding 1 ----------------------------------------------------------------
cat("## Finding 1 -- complement one-sided upper coverage: field / field-s / IJ two-term (Wilson)\n\n")
cat("| cell | n | field | field-s | IJ two-term | field-s - field | flips to cover / to miss | mean rho^c | lambda-SD^c/nSE unscaled -> studentized | field-s Wilson lower >= 0.92 | point >= 0.92 |\n|---|---|---|---|---|---|---|---|---|---|---|\n")
for (b in B) { r <- b$e1$results; d <- r[r$detected %in% 1L & is.finite(r$betaHhat_Hc) & is.finite(r$fld_Hc_up1s_s), ]; n <- nrow(d)
  cf <- mean(d$betaHhat_Hc <= d$fld_Hc_up1s); cs <- mean(d$betaHhat_Hc <= d$fld_Hc_up1s_s); ci <- mean(d$betaHhat_Hc <= exp(log(d$mr_Hc_est) + qnorm(.95) * d$mr_Hc_se_ij))
  cat(sprintf("| %s | %d | %s | %s | %s | %+.3f | %s / %s | %s | %s -> %s | %s | %s |\n", b$cl$label, n, wfmt(cf, n), wfmt(cs, n), wfmt(ci, n), cs - cf,
              f3(mean(d$betaHhat_Hc > d$fld_Hc_up1s & d$betaHhat_Hc <= d$fld_Hc_up1s_s)), f3(mean(d$betaHhat_Hc <= d$fld_Hc_up1s & d$betaHhat_Hc > d$fld_Hc_up1s_s)),
              f3(mean(d$fld_Hc_scale_ratio)), f3(sqrt(mean(d$fld_Hc_se^2) / mean(d$nv_Hc_se^2))), f3(sqrt(mean(d$fld_Hc_se_s^2) / mean(d$nv_Hc_se^2))), wil(cs, n)[1] >= 0.92, cs >= 0.92)) }
cat("\n")
# ---- Finding 2 ----------------------------------------------------------------
cat("## Finding 2 -- shape: complement upper coverage by |Hhat|/|H| tertile and by p-hat tertile (field -> field-s)\n\n")
cat("| cell | stratification | T1 range / T2 / T3 | n | mean rho^c T1 / T2 / T3 | fld/nv T1 / T2 / T3 | flds/nv T1 / T2 / T3 | field T1 / T2 / T3 | field-s T1 / T2 / T3 | field-s Wilson T1 / T2 / T3 |\n|---|---|---|---|---|---|---|---|---|---|\n")
for (b in B) { r <- b$e1$results; d <- r[r$detected %in% 1L & is.finite(r$betaHhat_Hc) & is.finite(r$fld_Hc_up1s_s), ]
  for (what in c("|Hhat|/|H|", "p-hat")) { v <- if (what == "p-hat") d$p_hat_H else d$n_harm / d$n_true; g <- tert(v); ks <- sort(unique(g))
    st <- function(fun) paste(vapply(ks, function(t) fun(d[g == t, ], v[g == t]), ""), collapse = " / ")
    cat(sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n", b$cl$label, what,
      st(function(dk, vk) sprintf("[%.2f, %.2f]", min(vk), max(vk))), st(function(dk, vk) sprintf("%d", nrow(dk))),
      st(function(dk, vk) f3(mean(dk$fld_Hc_scale_ratio))), st(function(dk, vk) f3(mean(dk$fld_Hc_se / dk$nv_Hc_se))), st(function(dk, vk) f3(mean(dk$fld_Hc_se_s / dk$nv_Hc_se))),
      st(function(dk, vk) f3(mean(dk$betaHhat_Hc <= dk$fld_Hc_up1s))), st(function(dk, vk) f3(mean(dk$betaHhat_Hc <= dk$fld_Hc_up1s_s))),
      st(function(dk, vk) { w <- wil(mean(dk$betaHhat_Hc <= dk$fld_Hc_up1s_s), nrow(dk)); sprintf("[%.3f, %.3f]", w[1], w[2]) }))) } }
cat("\n")
# ---- Finding 4: bound locations and the joint pair --------------------------------
cat("## Finding 4 -- bound locations (HR scale) and the joint pair\n\n")
cat("| cell | construction | mean beta(Hhat) | H lower mean | share >= 0.85 | share >= 0.95 | mean beta(Hhat^c) | Hc upper mean | share < 0.85 | share < 0.80 | margin_Hc (log) |\n|---|---|---|---|---|---|---|---|---|---|---|\n")
for (b in B) { d <- b$e1$results; d <- d[d$detected %in% 1L, ]
  rows <- list(naive = list(lo = exp(log(d$nv_H_est) - qnorm(.95) * d$nv_H_se), up = exp(log(d$nv_Hc_est) + qnorm(.95) * d$nv_Hc_se)),
               field = list(lo = d$fld_H_lo1s, up = d$fld_Hc_up1s), `field-s` = list(lo = d$fld_H_lo1s, up = d$fld_Hc_up1s_s),
               `IJ two-term` = list(lo = exp(log(d$mr_H_est) - qnorm(.95) * d$mr_H_se_ij), up = exp(log(d$mr_Hc_est) + qnorm(.95) * d$mr_Hc_se_ij)))
  for (nm in names(rows)) { p <- rows[[nm]]; ok <- is.finite(p$lo) & is.finite(p$up)
    cat(sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n", b$cl$label, nm, f3(mean(d$betaHhat_H)), f3(mean(p$lo[ok])), f3(mean(p$lo[ok] >= 0.85)), f3(mean(p$lo[ok] >= 0.95)),
                f3(mean(d$betaHhat_Hc)), f3(mean(p$up[ok])), f3(mean(p$up[ok] < 0.85)), f3(mean(p$up[ok] < 0.80)), f3(mean(log(p$up[ok]) - log(d$mr_Hc_est[ok]))))) } }
cat("\n**Joint pair** (coverage of (beta(Hhat) >= lower_H, beta(Hhat^c) <= upper_Hc), Wilson; margins in log units):\n\n")
cat("| cell | pair | n | joint cov [Wilson] | cov H | cov Hc | margin H | margin Hc | mean gamma | mean corr |\n|---|---|---|---|---|---|---|---|---|---|\n")
for (b in B) { r <- b$e1$results; d <- r[r$detected %in% 1L & is.finite(r$betaHhat_H) & is.finite(r$betaHhat_Hc) & is.finite(r$fld_joint_gamma) & is.finite(r$fld_joint_s_gamma), ]
  pairs <- list(`separate 95% (field)` = list(lo = d$fld_H_lo1s, up = d$fld_Hc_up1s, g = NA, cr = NA), `separate 95% (field-s)` = list(lo = d$fld_H_lo1s, up = d$fld_Hc_up1s_s, g = NA, cr = NA),
                `Bonferroni (joint)` = list(lo = d$fld_joint_bonf_loH, up = d$fld_joint_bonf_upHc, g = 0.025, cr = d$fld_joint_corr), `Bonferroni (joint_s)` = list(lo = d$fld_joint_s_bonf_loH, up = d$fld_joint_s_bonf_upHc, g = 0.025, cr = d$fld_joint_s_corr),
                `calibrated (joint)` = list(lo = d$fld_joint_loH, up = d$fld_joint_upHc, g = d$fld_joint_gamma, cr = d$fld_joint_corr), `calibrated (joint_s)` = list(lo = d$fld_joint_s_loH, up = d$fld_joint_s_upHc, g = d$fld_joint_s_gamma, cr = d$fld_joint_s_corr))
  for (nm in names(pairs)) { p <- pairs[[nm]]; cH <- d$betaHhat_H >= p$lo; cC <- d$betaHhat_Hc <= p$up; jc <- mean(cH & cC)
    cat(sprintf("| %s | %s | %d | %s | %s | %s | %s | %s | %s | %s |\n", b$cl$label, nm, nrow(d), wfmt(jc, nrow(d)), f3(mean(cH)), f3(mean(cC)), f3(mean(log(d$mr_H_est) - log(p$lo))), f3(mean(log(p$up) - log(d$mr_Hc_est))),
                if (all(is.na(p$g))) "-" else f3(mean(p$g)), if (all(is.na(p$cr))) "-" else f3(mean(p$cr)))) } }
cat("\n")
# ---- Finding 5 ----------------------------------------------------------------
cat("## Finding 5 -- se_field_s / (rho^c x se_field) per cell (R1 vs the global rescale; informational)\n\n")
cat("| cell | n | mean | sd | q01 | q10 | q50 | q90 | q99 | min | max | scale CV mean / q90 | mean rho^c | share rho^c > 1 | corr(se_s, nSE) / corr(se, nSE) | SD across reps: se_s / se / nSE |\n|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|\n")
for (b in B) { r <- b$e1$results; d <- r[r$detected %in% 1L & is.finite(r$fld_Hc_se_s) & is.finite(r$fld_Hc_scale_ratio), ]
  v <- d$fld_Hc_se_s / (d$fld_Hc_scale_ratio * d$fld_Hc_se); q <- quantile(v, c(.01, .10, .50, .90, .99))
  cat(sprintf("| %s | %d | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s / %s | %s | %s | %s / %s | %s / %s / %s |\n", b$cl$label, nrow(d), f4(mean(v)), f4(sd(v)), f4(q[1]), f4(q[2]), f4(q[3]), f4(q[4]), f4(q[5]), f4(min(v)), f4(max(v)),
              f3(mean(d$fld_Hc_scale_cv)), f3(quantile(d$fld_Hc_scale_cv, .9)), f3(mean(d$fld_Hc_scale_ratio)), f3(mean(d$fld_Hc_scale_ratio > 1)),
              f3(cor(d$fld_Hc_se_s, d$nv_Hc_se)), f3(cor(d$fld_Hc_se, d$nv_Hc_se)), f4(sd(d$fld_Hc_se_s)), f4(sd(d$fld_Hc_se)), f4(sd(d$nv_Hc_se)))) }
cat("\n")
# ---- identification context --------------------------------------------------
cat("## Identification context per cell (unchanged from the comparators by pairing)\n\n| cell | detected | mean |Hhat| | |Hhat|/|H| median | sens | spec | p-hat mean | complement fits/rep | fit+MR s/rep (e1stud) |\n|---|---|---|---|---|---|---|---|---|\n")
for (b in B) { r <- b$e1$results; d <- r[r$detected %in% 1L, ]
  cat(sprintf("| %s | %d | %.1f | %s | %s | %s | %s | %.0f | %.1f |\n", b$cl$label, nrow(d), mean(d$n_harm), f3(median(d$n_harm / d$n_true)), f3(mean(d$sens, na.rm = TRUE)), f3(mean(d$spec, na.rm = TRUE)), f3(mean(d$p_hat_H)), mean(d$fld_Hc_nfit, na.rm = TRUE), mean(d$fit_mr_secs, na.rm = TRUE))) }
sink(); cat("written", out, "\n")
