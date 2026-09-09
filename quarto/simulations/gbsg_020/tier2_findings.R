## tier2_findings.R -- the Stage 3 tables of TASK_tier2_mac_2026-09-08, verbatim
## source: the four campaign `tier2` combined bundles (12.5% prevalence, focus maxeffCons).
## Winner-only and winner-floor excluded everywhere.  Definitions match summary_tier2.qmd.
suppressMessages(library(forestsearch))
options(width = 200)
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || all(is.na(a))) b else a
z95 <- qnorm(0.95); z975 <- qnorm(0.975)
wil <- function(x, n, z = z975) { c <- (x + z^2/(2*n))/(1+z^2/n)
  h <- z*sqrt(x*(1-x)/n + z^2/(4*n^2))/(1+z^2/n); c(c - h, c + h) }
fmtw <- function(p, n) sprintf("%.3f [%.3f, %.3f]", p, wil(p, n)[1], wil(p, n)[2])
md <- function(df, digits = 3) {
  df <- as.data.frame(df)
  for (j in seq_along(df)) if (is.numeric(df[[j]])) df[[j]] <- formatC(df[[j]], format = "f", digits = digits)
  hdr <- paste0("| ", paste(names(df), collapse = " | "), " |")
  sep <- paste0("|", paste(rep("---", ncol(df)), collapse = "|"), "|")
  rows <- apply(df, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |"))
  cat(hdr, sep, rows, sep = "\n"); cat("\n\n")
}

cells <- list(
  list(key="h175_n500",  label="HR 1.75, n = 500",       kind="harm", hr=1.75, n=500,
       comp="results/fs_maxeffCons_fb_mr_field_m1_h175_knoise0_n500_s7c_combined_1_2000.rds",  comp_tag="s7c",
       f="results/fs_maxeffCons_fb_mr_field_m1_h175_knoise0_n500_tier2_combined_1_2000.rds"),
  list(key="h175_n1000", label="HR 1.75, n = 1000",      kind="harm", hr=1.75, n=1000,
       comp=NA_character_, comp_tag="none at this cell",
       f="results/fs_maxeffCons_fb_mr_field_m1_h175_knoise0_n1000_tier2_combined_1_2000.rds"),
  list(key="h175_n1500", label="HR 1.75, n = 1500",      kind="harm", hr=1.75, n=1500,
       comp=NA_character_, comp_tag="none exact (map1c h150 n1500 is a different HR)",
       f="results/fs_maxeffCons_fb_mr_field_m1_h175_knoise0_n1500_tier2_combined_1_2000.rds"),
  list(key="h100_n500",  label="HR 1.00, n = 500 (null)", kind="null", hr=1.00, n=500,
       comp="results/fs_maxeffCons_fb_mr_field_m1_h100_knoise0_n500_s7c_combined_1_2000.rds",  comp_tag="s7c null",
       f="results/fs_maxeffCons_fb_mr_field_m1_h100_knoise0_n500_tier2_combined_1_2000.rds"))
cells <- Filter(function(cl) file.exists(cl$f), cells)
B <- lapply(cells, function(cl) readRDS(cl$f)); names(B) <- vapply(cells, `[[`, "", "key")
cl_f <- setNames(lapply(cells, `[[`, "f"), names(B))
CMP <- lapply(cells, function(cl) if (!is.na(cl$comp) && file.exists(cl$comp)) readRDS(cl$comp) else NULL)
names(CMP) <- names(B)
labs <- setNames(vapply(cells, `[[`, "", "label"), names(B))
kind <- setNames(vapply(cells, `[[`, "", "kind"), names(B))
ctag <- setNames(vapply(cells, `[[`, "", "comp_tag"), names(B))
det  <- lapply(B, function(b) b$results[b$results$detected %in% 1L, , drop = FALSE])

## ---------------------------------------------------------------- GATE 2 ----
cat("## Gate 2 per cell\n\n")
G2 <- do.call(rbind, lapply(names(B), function(k) {
  r <- B[[k]]$results; m <- B[[k]]$meta; d <- det[[k]]
  fin <- function(cols) all(vapply(intersect(cols, names(d)), function(cc) all(is.finite(as.numeric(d[[cc]]))), TRUE))
  inv_ok <- all(d$fld_H_lo2s <= d$fld_H_lo1s + 1e-12, d$fld_H_lo1s <= d$fld_H_est2 + 1e-12,
                d$fld_H_est2 <= d$fld_H_hi2s + 1e-12,
                d$fld_Hc_lo2s <= d$fld_Hc_lo1s + 1e-12, d$fld_Hc_lo1s <= d$fld_Hc_est2 + 1e-12,
                d$fld_Hc_est2 <= d$fld_Hc_up1s + 1e-12, d$fld_Hc_up1s <= d$fld_Hc_hi2s + 1e-12,
                d$fld_Hc_lo2s_s <= d$fld_Hc_lo1s_s + 1e-12, d$fld_Hc_lo1s_s <= d$fld_Hc_est2_s + 1e-12,
                d$fld_Hc_est2_s <= d$fld_Hc_up1s_s + 1e-12, d$fld_Hc_up1s_s <= d$fld_Hc_hi2s_s + 1e-12)
  ## bound <-> quantile identities (within bundle, machine-local)
  id1 <- max(abs(log(d$fld_Hc_up1s) - (log(d$fld_Hc_est2) + (d$fld_Hc_lam_mean - d$fld_Hc_q05))), na.rm = TRUE)
  id2 <- max(abs(log(d$fld_H_lo1s)  - (log(d$fld_H_est2)  + (d$fld_H_lam_mean  - d$fld_H_q95))),  na.rm = TRUE)
  ## the Gate 2 bar is per replicate: 0.95 - 2/n_joint with THAT replicate's own n_joint
  bar   <- 0.95 - 2/d$fld_joint_n;   bad   <- sum(d$fld_joint_prob   < bar - 1e-12)
  bar_s <- 0.95 - 2/d$fld_joint_s_n; bad_s <- sum(d$fld_joint_s_prob < bar_s - 1e-12)
  ## machine: the combine meta drops hostname, so read it from this cell's batch bundles
  hosts <- unique(vapply(Sys.glob(sub("_combined_1_2000\\.rds$", "_res_*.rds", cl_f[[k]])),
                         function(z) readRDS(z)$meta$hostname %||% NA_character_, ""))
  cc <- CMP[[k]]
  nt <- if (is.null(cc)) NA else {
    rc <- cc$results[order(cc$results$sim_id), ]; rf <- r[order(r$sim_id), ]
    identical(as.integer(rc$n_true), as.integer(rf$n_true)) }
  tr <- if (is.null(cc)) NA_real_ else max(vapply(names(cc$truth), function(nm)
    abs(as.numeric(B[[k]]$truth[[nm]]) - as.numeric(cc$truth[[nm]]))/abs(as.numeric(cc$truth[[nm]])), 0))
  data.frame(cell = labs[[k]], rows = nrow(r),
    sim_id = sprintf("%d-%d", min(r$sim_id), max(r$sim_id)),
    dups = anyDuplicated(r$sim_id),
    config_err = sum(grepl("CONFIG", as.character(r$status)) | grepl("CONFIG", as.character(r$err_msg %||% ""))),
    detected = sum(r$detected %in% 1L),
    prevalence = mean(r$n_true)/m$n_sample,
    all_finite = fin(c("fld_H_est2","fld_H_lo1s","fld_H_se","fld_Hc_est2","fld_Hc_up1s","fld_Hc_se",
                       "fld_Hc_est2_s","fld_Hc_up1s_s","fld_Hc_se_s","fld_Hc_scale_ratio",
                       "fld_joint_gamma","fld_joint_s_gamma","fld_joint_bonf_loH","fld_joint_s_bonf_upHc",
                       "mr_H_se_ij","mr_Hc_se_ij","betaHhat_H","betaHhat_Hc","p_hat_H","fld_Hc_nfit")),
    invariants = inv_ok,
    gamma_range = sprintf("[%.4f, %.4f]", min(d$fld_joint_gamma), max(d$fld_joint_gamma)),
    gamma_s_range = sprintf("[%.4f, %.4f]", min(d$fld_joint_s_gamma), max(d$fld_joint_s_gamma)),
    joint_n_range = sprintf("%d-%d", min(d$fld_joint_n), max(d$fld_joint_n)),
    joint_below_bar = bad, joint_s_below_bar = bad_s,
    bound_id_max = max(id1, id2),
    n_true_identical = nt, truth_max_reldiff = if (is.na(tr)) NA_character_ else sprintf("%.2e", tr),
    machine = paste(hosts, collapse = ","), version = m$forestsearch_version,
    knobs = sprintf("%s/%s/%s/%s/J=%d/z1q=%.2f", m$field_complement, m$field_decompose,
                    m$field_scale_complement, m$ij_residual, m$er_jcuts, m$harm_z1_quantile),
    stringsAsFactors = FALSE) }))
md(G2, 6)

## ------------------------------------------------- 1. STANDARD TABLES -------
cat("## 1. Constructions per cell (identical replicates): naive / field / field-s / IJ two-term\n\n")
subst_s <- function(r) { for (s in c("est2","up1s","lo1s","lo2s","hi2s","lo_se","hi_se","se","lam_mean"))
  r[[paste0("fld_Hc_", s)]] <- r[[paste0("fld_Hc_", s, "_s")]]; r }
err_sd <- function(r, block, est) { d <- r[r$detected %in% 1L, ]
  e <- switch(est, naive = d[[paste0("nv_",block,"_est")]], mr = d[[paste0("mr_",block,"_est")]],
              fld = d[[paste0("fld_",block,"_est2")]])
  ok <- is.finite(e) & is.finite(d[[paste0("betaHhat_",block)]])
  sd(log(e[ok]) - log(d[[paste0("betaHhat_",block)]][ok])) }
cov_block <- function(r, block, side) {
  t1 <- fs_sim_bias_coverage(r, block = block, estimators = c("naive","mr","fld"), side = side)
  t1$sd_err <- vapply(as.character(t1$estimator), function(e) err_sd(r, block, e), 0)
  t1$construction <- c(naive="naive", mr="IJ two-term", fld="field")[as.character(t1$estimator)]
  if (block == "Hc") { rs <- subst_s(r)
    t2 <- fs_sim_bias_coverage(rs, block = "Hc", estimators = "fld", side = side)
    t2$sd_err <- err_sd(rs, "Hc", "fld"); t2$construction <- "field-s"; t1 <- rbind(t1, t2) }
  t1$se_over_sd_err <- t1$se_mean / t1$sd_err; t1 }
COV <- do.call(rbind, lapply(names(B), function(k) { r <- B[[k]]$results
  rbind(cbind(cell = labs[[k]], block = "Hhat (lower)",   cov_block(r, "H",  "lower")),
        cbind(cell = labs[[k]], block = "Hhat^c (upper)", cov_block(r, "Hc", "upper"))) }))
COV$construction <- factor(COV$construction, levels = c("naive","field","field-s","IJ two-term"))
COV <- COV[order(match(COV$cell, labs), COV$block, COV$construction), ]
COV$cov1_wilson <- sprintf("%.3f [%.3f, %.3f]", COV$cov1, COV$cov1_wilson_lo, COV$cov1_wilson_hi)
md(COV[, c("cell","block","construction","n","bias_log","sd_emp","sd_err","se_mean","b","r",
           "se_over_sd_err","cov1_wilson","cov2","cov1_ref")], 3)

cat("**Across cells** (one-sided coverage min / mean / max; mean two-sided; mean b; mean r; mean SE/error SD):\n\n")
ACR <- do.call(rbind, lapply(split(COV, list(COV$block, COV$construction), drop = TRUE), function(d)
  data.frame(block = d$block[1], construction = as.character(d$construction[1]), cells = nrow(d),
             cov1_min = min(d$cov1), cov1_mean = mean(d$cov1), cov1_max = max(d$cov1),
             cov2_mean = mean(d$cov2), b_mean = mean(d$b), r_mean = mean(d$r),
             se_over_sd_err_mean = mean(d$se_over_sd_err), stringsAsFactors = FALSE)))
ACR <- ACR[order(ACR$block, ACR$construction), ]; md(ACR, 3)

## --------------------------------------- 2. THE DOMINATED-REGIME CHECK ------
cat("## 2. The dominated-regime check: rho^c, the SD ratios, field vs field-s\n\n")
D2 <- do.call(rbind, lapply(names(B), function(k) { d <- det[[k]]
  d <- d[is.finite(d$betaHhat_Hc) & is.finite(d$fld_Hc_up1s_s), ]
  rc <- d$fld_Hc_scale_ratio
  cf <- mean(d$betaHhat_Hc <= d$fld_Hc_up1s); cs <- mean(d$betaHhat_Hc <= d$fld_Hc_up1s_s)
  wf <- wil(cf, nrow(d)); ws <- wil(cs, nrow(d))
  data.frame(cell = labs[[k]], n = nrow(d),
             rho_c_mean = mean(rc), rho_c_q10 = quantile(rc, .10, names = FALSE),
             rho_c_q90 = quantile(rc, .90, names = FALSE), rho_c_share_gt1 = mean(rc > 1),
             lamSD_over_nSE = sqrt(mean(d$fld_Hc_se^2))/sqrt(mean(d$nv_Hc_se^2)),
             lamSDs_over_nSE = sqrt(mean(d$fld_Hc_se_s^2))/sqrt(mean(d$nv_Hc_se^2)),
             field = cf, field_lo = wf[1], field_hi = wf[2],
             field_s = cs, field_s_lo = ws[1], field_s_hi = ws[2],
             diff = cs - cf, wilson_overlap = (ws[1] <= wf[2]) && (wf[1] <= ws[2]),
             flip_to_cover = mean(d$betaHhat_Hc > d$fld_Hc_up1s & d$betaHhat_Hc <= d$fld_Hc_up1s_s),
             flip_to_miss  = mean(d$betaHhat_Hc <= d$fld_Hc_up1s & d$betaHhat_Hc > d$fld_Hc_up1s_s),
             stringsAsFactors = FALSE) }))
md(D2, 4)

## ------------------------------------- 3. AGAINST THE COMMITTED RECORDS -----
cat("## 3. Against the committed unscaled-field comparators (pairing by DGM draws, not fitted values)\n\n")
D3 <- do.call(rbind, lapply(names(B), function(k) { d <- det[[k]]
  cf <- mean(d$betaHhat_Hc <= d$fld_Hc_up1s, na.rm = TRUE); wf <- wil(cf, nrow(d))
  cc <- CMP[[k]]
  if (is.null(cc)) return(data.frame(cell = labs[[k]], comparator = ctag[[k]], n = nrow(d),
    tier2_field = cf, tier2_lo = wf[1], tier2_hi = wf[2],
    comp_field = NA_real_, comp_lo = NA_real_, comp_hi = NA_real_, diff = NA_real_,
    wilson_overlap = NA, max_absdiff_up1s = NA_character_, max_absdiff_betaHc = NA_character_,
    stringsAsFactors = FALSE))
  dc <- cc$results[cc$results$detected %in% 1L, ]
  cg <- mean(dc$betaHhat_Hc <= dc$fld_Hc_up1s, na.rm = TRUE); wg <- wil(cg, nrow(dc))
  ## informational only (never asserted across machines): row-wise agreement of the
  ## unscaled field bound and of beta(Hhat^c) on the paired sim_ids
  rf <- B[[k]]$results[order(B[[k]]$results$sim_id), ]; rc <- cc$results[order(cc$results$sim_id), ]
  ok <- rf$detected %in% 1L & rc$detected %in% 1L
  mx_up <- max(abs(rf$fld_Hc_up1s[ok] - rc$fld_Hc_up1s[ok]), na.rm = TRUE)
  mx_bt <- max(abs(rf$betaHhat_Hc[ok] - rc$betaHhat_Hc[ok]), na.rm = TRUE)
  data.frame(cell = labs[[k]], comparator = ctag[[k]], n = nrow(d),
             tier2_field = cf, tier2_lo = wf[1], tier2_hi = wf[2],
             comp_field = cg, comp_lo = wg[1], comp_hi = wg[2], diff = cf - cg,
             wilson_overlap = (wf[1] <= wg[2]) && (wg[1] <= wf[2]),
             max_absdiff_up1s = sprintf("%.2e", mx_up), max_absdiff_betaHc = sprintf("%.2e", mx_bt),
             stringsAsFactors = FALSE) }))
md(D3, 4)

## ------------------------------------------- 4. BY P-HAT TERTILE ------------
cat("## 4. By p-hat tertile, both blocks (observed with Wilson; Gaussian-implied beside)\n\n")
tert <- function(v) { br <- quantile(v, c(0, 1/3, 2/3, 1), names = FALSE)
  cut(v, breaks = unique(br), include.lowest = TRUE, labels = FALSE) }
gauss <- function(err, se, side) { b <- mean(err)/sd(err); rr <- mean(se)/sd(err)
  if (side == "lower") pnorm(z95*rr - b) else pnorm(z95*rr + b) }
D4 <- do.call(rbind, lapply(names(B), function(k) { d <- det[[k]]
  d <- d[is.finite(d$p_hat_H) & is.finite(d$betaHhat_H) & is.finite(d$betaHhat_Hc), ]
  g <- tert(d$p_hat_H)
  do.call(rbind, lapply(sort(unique(g)), function(t) { dk <- d[g == t, ]
    cH  <- mean(dk$betaHhat_H  >= dk$fld_H_lo1s)
    cC  <- mean(dk$betaHhat_Hc <= dk$fld_Hc_up1s)
    cCs <- mean(dk$betaHhat_Hc <= dk$fld_Hc_up1s_s)
    wH <- wil(cH, nrow(dk)); wC <- wil(cC, nrow(dk)); wCs <- wil(cCs, nrow(dk))
    eH  <- log(dk$fld_H_est2)   - log(dk$betaHhat_H)
    eC  <- log(dk$fld_Hc_est2)  - log(dk$betaHhat_Hc)
    eCs <- log(dk$fld_Hc_est2_s)- log(dk$betaHhat_Hc)
    data.frame(cell = labs[[k]], tertile = sprintf("T%d [%.3f, %.3f]", t, min(dk$p_hat_H), max(dk$p_hat_H)),
      n = nrow(dk), p_hat_mean = mean(dk$p_hat_H), rho_c = mean(dk$fld_Hc_scale_ratio),
      harm_field = cH, harm_lo = wH[1], harm_hi = wH[2], harm_gauss = gauss(eH, dk$fld_H_se, "lower"),
      comp_field = cC, comp_lo = wC[1], comp_hi = wC[2], comp_gauss = gauss(eC, dk$fld_Hc_se, "upper"),
      comp_field_s = cCs, comp_s_lo = wCs[1], comp_s_hi = wCs[2],
      comp_s_gauss = gauss(eCs, dk$fld_Hc_se_s, "upper"), stringsAsFactors = FALSE) })) }))
md(D4, 3)

## --------------------------------------- 5. IDENTIFICATION AND n ------------
cat("## 5. Identification and n\n\n")
D5 <- do.call(rbind, lapply(names(B), function(k) { r <- B[[k]]$results; d <- det[[k]]; tr <- B[[k]]$truth
  rel <- d$n_harm / d$n_true
  data.frame(cell = labs[[k]], n = B[[k]]$meta$n_sample, rows = nrow(r),
    detection = mean(r$detected %in% 1L),
    n_true_mean = mean(r$n_true), n_sel_mean = mean(d$n_sel),
    relsize_median = median(rel), relsize_q90 = quantile(rel, .90, names = FALSE),
    sens = mean(d$sens, na.rm = TRUE), spec = mean(d$spec, na.rm = TRUE),
    ppv = mean(d$ppv, na.rm = TRUE), npv = mean(d$npv, na.rm = TRUE),
    betaHhat_H_mean = mean(d$betaHhat_H), planted_marg_H = as.numeric(tr$marg_H),
    planted_cde_H = as.numeric(tr$cde_H),
    betaHhat_Hc_mean = mean(d$betaHhat_Hc), planted_marg_Hc = as.numeric(tr$marg_Hc),
    planted_cde_Hc = as.numeric(tr$cde_Hc),
    p_hat_mean = mean(d$p_hat_H), n_family_mean = mean(d$n_family),
    comp_fits_mean = mean(d$fld_Hc_nfit), share_new_fit = mean(d$fld_Hc_share_newfit),
    stringsAsFactors = FALSE) }))
md(D5, 3)

## ----------------------------------- 6. PRE-REGISTERED ACCEPTANCE CRITERIA --
cat("## 6. The pre-registered acceptance criteria, evaluated as findings\n\n")
D6 <- do.call(rbind, lapply(names(B), function(k) { d <- det[[k]]; nn <- nrow(d)
  hf <- mean(d$betaHhat_H >= d$fld_H_lo1s); whf <- wil(hf, nn)
  fs <- mean(d$betaHhat_Hc <= d$fld_Hc_up1s_s); wfs <- wil(fs, nn)
  ij2H <- mean(d$betaHhat_H >= d$mr_H_lo & d$betaHhat_H <= d$mr_H_hi)
  ij2C <- mean(d$betaHhat_Hc >= d$mr_Hc_lo & d$betaHhat_Hc <= d$mr_Hc_hi)
  jb  <- mean(d$betaHhat_H >= d$fld_joint_bonf_loH & d$betaHhat_Hc <= d$fld_joint_bonf_upHc)
  jbs <- mean(d$betaHhat_H >= d$fld_joint_s_bonf_loH & d$betaHhat_Hc <= d$fld_joint_s_bonf_upHc)
  isnull <- kind[[k]] == "null"
  data.frame(cell = labs[[k]], n = nn,
    harm_field = hf, harm_wilson_lo = whf[1],
    harm_met = if (isnull) "recorded (null)" else if (whf[1] >= 0.94) "MET" else "NOT MET",
    field_s = fs, field_s_wilson_lo = wfs[1],
    field_s_met = if (isnull) "recorded (null)" else if (fs >= 0.93) "MET" else "NOT MET",
    ij2_H = ij2H, ij2_Hc = ij2C,
    ij2_met = if (isnull) "recorded (null)" else if (min(ij2H, ij2C) >= 0.93) "MET" else "NOT MET",
    joint_bonf = jb, joint_s_bonf = jbs,
    joint_met = if (isnull) "recorded (null)" else if (jb >= 0.93) "MET" else "NOT MET",
    stringsAsFactors = FALSE) }))
md(D6, 4)

cat("## Reading lines\n\n")
for (i in seq_len(nrow(D2))) cat(sprintf(
 "- %s (n %d): rho^c mean %.4f [q10 %.4f, q90 %.4f], share>1 %.3f; lambda-SD^c/naive %.4f, se_field_s/naive %.4f; field %s vs field-s %s (diff %+.4f, Wilson overlap %s)\n",
 D2$cell[i], D2$n[i], D2$rho_c_mean[i], D2$rho_c_q10[i], D2$rho_c_q90[i], D2$rho_c_share_gt1[i],
 D2$lamSD_over_nSE[i], D2$lamSDs_over_nSE[i], fmtw(D2$field[i], D2$n[i]), fmtw(D2$field_s[i], D2$n[i]),
 D2$diff[i], D2$wilson_overlap[i]))
cat("\n")
for (i in seq_len(nrow(D6))) cat(sprintf(
 "- %s: harm field %s [%s]; field-s %s [%s]; IJ two-sided H %.3f / Hc %.3f [%s]; joint Bonferroni %.3f (joint_s %.3f) [%s]\n",
 D6$cell[i], fmtw(D6$harm_field[i], D6$n[i]), D6$harm_met[i], fmtw(D6$field_s[i], D6$n[i]), D6$field_s_met[i],
 D6$ij2_H[i], D6$ij2_Hc[i], D6$ij2_met[i], D6$joint_bonf[i], D6$joint_s_bonf[i], D6$joint_met[i]))
