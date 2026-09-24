# declcalc0_findings.R -- tables for REPORT_declcal_c0_campaign_2026-09-22
# (TASK_declcal_c0_campaign_2026-09-22, section 6).  Transplant of
# declcal_findings.R: same helpers (wilson, fmt_w, med_iqr), same cell-status
# rule (a cell whose payload is missing, incomplete or tripped contributes no
# rate).  Reads results/declcalc0_{inull,power}_<cell>_res_1_2000.rds, the
# committed results/declcal_{inull,power}_<cell>_res_1_2000.rds and
# logs/declcal_fixedk_practical.txt; writes logs/declcalc0_findings.txt
# (markdown tables) and declcalc0_findings.rds.  Run from scripts_dinamr/.
wilson <- function(x, n, z = qnorm(0.975)) {
  p <- x / n; d <- 1 + z^2 / n
  c((p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / d,
    (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / d)
}
fmt_w <- function(x) { x <- x[!is.na(x)]; n <- length(x); k <- sum(x); w <- wilson(k, n)
  sprintf("%.4f [%.4f, %.4f]", k / n, w[1], w[2]) }
med_iqr <- function(x, d = 0) { q <- stats::quantile(x, c(.5, .25, .75), na.rm = TRUE)
  sprintf(paste0("%.", d, "f (%.", d, "f-%.", d, "f)"), q[1], q[2], q[3]) }
cells <- list(inull = c("B1", "B2", "B3", "B4", "B5", "B6"),
              power = c("C1", "C2", "C3", "C4"))
dgm_txt <- c(null0657 = "uniform benefit HR 0.657", null0721 = "uniform benefit HR 0.721",
             alt_h150 = "planted harm HR 1.5", alt_h200 = "planted harm HR 2.0")
c0_grid <- c(0.70, 0.75, 0.80, 0.85)
sfx <- sprintf("c%03d", as.integer(round(100 * c0_grid)))
id_cols <- c("max_T_pre", "max_T_post", "G_pre", "G_post", "declared_conv",
             "declared_conv_exact", "Mstar_q90", "Mstar_q95", "Mstar_q99",
             "kappa_hat_05", "kappa_hat_10", "declared_cal05", "declared_cal10")

# Every column this script reads must be present: a payload from the fixed
# drivers (TASK_declcal_consumers_2026-09-24_v3) no longer carries the
# exact-cutoff columns, and a missing column would otherwise reach sprintf()
# as a zero-length or NA value without any error.
need_cols <- function(df, cols, what) {
  miss <- setdiff(cols, names(df))
  if (length(miss)) stop(sprintf("%s lacks column(s) read by this script: %s", what,
                                 paste(miss, collapse = ", ")), call. = FALSE)
  invisible(TRUE)
}
c0_read <- c("kappa_hat_05", "kappa_hat_10", "declared_cal05", "declared_cal10",
             "n_admitted_cal05", "fw_1621", "fw_1645")
cols_res <- unique(c("status", "dgm", "n", "rep", "declared_conv", "declared_cal05", "declared_cal10",
                     "kappa_hat_05", "kappa_hat_10", "max_T_pre", "max_T_post", "n_admitted_cal05",
                     "alpha_FW_hat_1621", "alpha_FW_hat_1645", "pconsistency_digits", "wall_sec", id_cols,
                     as.vector(outer(c0_read, sfx, paste, sep = "_"))))
cols_aux <- c("rep", "search_declared", "c0_mono_draws", "c0_kappa_le")
cols_ref <- unique(c("rep", "status", id_cols, "declared_cal05", "declared_cal10", "max_T_post"))
AX <- readRDS("../results/declcal_c0approx_res.rds")$approx   # committed at e9de87c7
# Settable p* (TASK_declcal_consumers_2026-09-24_v3), derived from the stored
# kappa_hat by the package's own inverse (.fs_decl_settable_table(), built on
# .fs_decl_settable()): pstar_settable, the smallest p* on the fit's
# pconsistency.digits grid whose rounded-screen cutoff reaches kappa_hat (NA:
# none <= 1 at those digits), and the finer pair (pstar_fine, digits_fine), the
# smallest digits at which a settable p* sits within 0.01 z of kappa_hat.
stab <- function(k, d = dig) { d <- unique(as.integer(d)); stopifnot(length(d) == 1L, !is.na(d))
  forestsearch:::.fs_decl_settable_table(k, d) }
set_txt <- function(t) sprintf("%s ; none <= 1: %.4f ; fine %.5f [digits %s]",
  if (any(t$pstar_achievable)) med_iqr(t$pstar_settable[t$pstar_achievable], 2) else "none",
  mean(!t$pstar_achievable), median(t$pstar_fine),
  paste(unique(range(t$digits_fine)), collapse = "-"))
ptxt <- function(k) { t <- stab(k)
  if (t$pstar_achievable) sprintf("%.2f", t$pstar_settable) else
    sprintf("none <= 1 (fine %.5f [digits %d])", t$pstar_fine, t$digits_fine) }
FOOT <- paste0("Settable p* = the smallest p* on the pconsistency.digits grid whose rounded-screen cutoff ",
               "reaches kappa_hat (forestsearch:::.fs_decl_settable()); 'none <= 1' when no p* <= 1 reaches ",
               "it at those digits; pstar_fine [digits_fine] is the settable p* at the smallest digits within 0.01 z ",
               "of kappa_hat. Setting that p* at those digits reproduces a cutoff at or above kappa_hat.")
P <- list(); R0 <- list(); status <- list()
for (blk in names(cells)) for (cl in cells[[blk]]) {
  f <- sprintf("../results/declcalc0_%s_%s_res_1_2000.rds", blk, cl)
  R0[[cl]] <- readRDS(sprintf("../results/declcal_%s_%s_res_1_2000.rds", blk, cl))
  need_cols(R0[[cl]]$results, cols_ref, sprintf("declcal_%s_%s $results", blk, cl))
  if (!file.exists(f)) { status[[cl]] <- "not run"; next }
  p <- readRDS(f)
  need_cols(p$results, cols_res, paste0(f, " $results")); need_cols(p$aux, cols_aux, paste0(f, " $aux"))
  st <- p$meta$cell_status; fin <- isTRUE(p$meta$final)
  status[[cl]] <- sprintf("%s%s ; %d rows ; %s ; elapsed %.0f s ; median wall %.2f s ; max wall %.1f s",
                          st, if (fin) "" else " (not final)", nrow(p$results),
                          paste(names(table(p$results$status)), table(p$results$status), collapse = ", "),
                          p$meta$elapsed_s, median(p$results$wall_sec), max(p$results$wall_sec))
  if (identical(st, "complete") && fin) P[[cl]] <- p
}
ok_rows <- function(cl) { r <- P[[cl]]$results; r[r$status == "ok", ] }
dig <- unique(unlist(lapply(P, function(p) unique(p$results$pconsistency_digits))))
stopifnot(length(dig) == 1L, !is.na(dig))
# the executed (rounded) cutoff as the payload recorded it
z_exec <- function(cl) { z <- P[[cl]]$meta$z_round; if (is.null(z)) z <- unique(P[[cl]]$results$z_pstar)
  stopifnot(length(z) == 1L, is.finite(z)); z }
out <- character(0); say <- function(...) { s <- sprintf(...); cat(s, "\n"); out <<- c(out, s) }
say("## Cell status"); for (cl in names(status)) say("- %s: %s", cl, status[[cl]])

# ---- gates --------------------------------------------------------------------
say("\n## Gates")
say("| cell | ok reps | identity: reps disagreeing on any of %d columns | Mstar_c0 monotone in c0 (all draws) | kappa_hat_05_c0 <= kappa_hat_05 | fidelity (declared_conv = search indicator) | abort / error |",
    length(id_cols))
say("|---|---|---|---|---|---|---|")
gates <- list()
for (cl in names(P)) {
  r <- P[[cl]]$results; a <- P[[cl]]$aux[match(r$rep, P[[cl]]$aux$rep), ]
  ref <- R0[[cl]]$results; m <- ref[match(r$rep, ref$rep), ]
  bad <- rep(FALSE, nrow(r))
  for (v in id_cols) { x <- r[[v]]; y <- m[[v]]
    bad <- bad | !((is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & x == y)) }
  ok <- r$status == "ok"
  gates[[cl]] <- c(n_ok = sum(ok), id_bad = sum(bad), mono = sum(a$c0_mono_draws %in% TRUE),
                   kap = sum(a$c0_kappa_le %in% TRUE), fid = sum(r$declared_conv[ok] == a$search_declared[ok]),
                   fail = sum(r$status %in% c("abort_time", "error")))
  say("| %s | %d | %d | %d / %d | %d / %d | %d / %d | %d |", cl, sum(ok), sum(bad),
      sum(a$c0_mono_draws %in% TRUE), nrow(r), sum(a$c0_kappa_le %in% TRUE), nrow(r),
      gates[[cl]][["fid"]], sum(ok), gates[[cl]][["fail"]])
}

# ---- 6.1 primary tables -------------------------------------------------------
for (al in c("05", "10")) {
  say("\n## 6.1 Primary table, calibrated alpha = 0.%s", al)
  say("| cell | DGM | n | p* = 0.90 as executed | %s | c0 = c2 (committed) |",
      paste(sprintf("cal alpha 0.%s, c0 %.2f", al, c0_grid), collapse = " | "))
  say("|---|---|---|---|%s---|", strrep("---|", length(c0_grid)))
  for (cl in names(P)) {
    r <- ok_rows(cl); ref <- R0[[cl]]$results; ref <- ref[ref$status == "ok", ]
    say("| %s | %s | %d | %s | %s | %s |", cl, dgm_txt[[r$dgm[1]]], r$n[1], fmt_w(r$declared_conv),
        paste(vapply(sfx, function(s) fmt_w(r[[sprintf("declared_cal%s_%s", al, s)]]), ""), collapse = " | "),
        fmt_w(ref[[sprintf("declared_cal%s", al)]]))
    fmt_r <- function(p) fmt_w(c(rep(1L, round(p * 2000)), rep(0L, 2000 - round(p * 2000))))
    say("| %s approx | plug-in median kappa_hat(c0) at n, max_T_pre | %d | - | %s | - |", cl, r$n[1],
        paste(vapply(AX$rates[sprintf("%s|%.2f|%s", cl, c0_grid, al)], fmt_r, ""), collapse = " | "))
  }
  say("\n\"approx\" rows: REPORT_declcal_c0_approx_2026-09-22 (e9de87c7) -- each rate is mean(max_T_pre >= median kappa_hat(c0) over 40 captures at that n, B 2000), a fixed cutoff, not the per-replicate calibrated rule of the exact row above it.")
}

# ---- 6.2 calibration quantities per c0 ---------------------------------------
say("\n## 6.2 The calibration's own quantities per c0")
say("| cell | c0 | kappa_hat_05 (median, IQR) | settable p*: median (IQR) ; share none <= 1 ; median pstar_fine [digits_fine] | mean fw_1645 | mean fw_1621 | n_admitted_cal05 median (all reps) | declaring reps (cal05) | share of declaring reps with n_admitted_cal05 = 1 |")
say("|---|---|---|---|---|---|---|---|---|")
q62 <- list()
for (cl in names(P)) {
  r <- ok_rows(cl)
  for (j in seq_along(sfx)) {
    s <- sfx[j]; dec <- r[[paste0("declared_cal05_", s)]] == 1L
    na <- r[[paste0("n_admitted_cal05_", s)]]
    sh1 <- if (any(dec)) sprintf("%.4f (%d / %d)", mean(na[dec] == 1L), sum(na[dec] == 1L), sum(dec)) else "no declarations"
    say("| %s | %.2f | %s | %s | %.4f | %.4f | %d | %d | %s |", cl, c0_grid[j],
        med_iqr(r[[paste0("kappa_hat_05_", s)]], 3), set_txt(stab(r[[paste0("kappa_hat_05_", s)]])),
        mean(r[[paste0("fw_1645_", s)]]), mean(r[[paste0("fw_1621_", s)]]),
        as.integer(median(na)), sum(dec), sh1)
    q62[[paste(cl, s)]] <- data.frame(cell = cl, n = r$n[1], c0 = c0_grid[j],
      kappa_med = median(r[[paste0("kappa_hat_05_", s)]]), n_dec = sum(dec),
      n_dec_adm1 = sum(na[dec] == 1L), n_dec_adm_gt1 = sum(na[dec] > 1L))
  }
  say("| %s | c2 = 1.00 (committed) | %s | %s | %.4f | %.4f | %d | %d | %s |", cl,
      med_iqr(r$kappa_hat_05, 3), set_txt(stab(r$kappa_hat_05)), mean(r$alpha_FW_hat_1645),
      mean(r$alpha_FW_hat_1621), as.integer(median(r$n_admitted_cal05)), sum(r$declared_cal05),
      if (any(r$declared_cal05 == 1)) sprintf("%.4f (%d / %d)", mean(r$n_admitted_cal05[r$declared_cal05 == 1] == 1),
                                             sum(r$n_admitted_cal05[r$declared_cal05 == 1] == 1), sum(r$declared_cal05)) else "no declarations")
}
q62 <- do.call(rbind, q62)
say("\n%s", FOOT)

# ---- 6.3 the trade, c0 as rows ---------------------------------------------
fk <- readLines("logs/declcal_fixedk_practical.txt")
fk <- grep("^\\| [BC][0-9] \\| [0-9]+ \\| 2\\.0000 \\|", fk, value = TRUE)
fk <- do.call(rbind, lapply(strsplit(fk, "\\|"), function(z) { z <- trimws(z)
  data.frame(cell = z[2], pre = as.numeric(sub(" .*", "", z[6])), post = as.numeric(sub(" .*", "", z[7]))) }))
# cross-check the parsed k = 2.0 rates against the committed payloads
fk_chk <- vapply(fk$cell, function(cl) { r <- R0[[cl]]$results
  isTRUE(all.equal(mean(ifelse(is.na(r$max_T_post), FALSE, r$max_T_post >= 2.0)),
                   fk$post[fk$cell == cl], tolerance = 1e-9)) }, logical(1))
Bc <- intersect(cells$inull, names(P)); Cc <- intersect(cells$power, names(P))
trade_row <- function(lab, rate_fun) {
  rb <- vapply(Bc, rate_fun, 0); rc <- vapply(Cc, rate_fun, 0)
  w <- which.max(rb)
  sprintf("| %s | %.4f (%s) | %s |", lab, rb[w], names(rb)[w],
          paste(sprintf("%.4f", rc[c("C1", "C2", "C3", "C4")]), collapse = " | "))
}
for (al in c("05", "10")) {
  say("\n## 6.3 The trade, calibrated alpha = 0.%s: worst uniform-benefit false-declaration rate (B1-B6) and planted-harm power", al)
  say("| screen | worst B rate (cell) | HR 1.5 n 1000 (C1) | HR 1.5 n 1500 (C2) | HR 2.0 n 1000 (C3) | HR 2.0 n 1500 (C4) |")
  say("|---|---|---|---|---|---|")
  for (j in seq_along(sfx)) {
    say(trade_row(sprintf("calibrated, c0 = %.2f", c0_grid[j]),
                  function(cl) mean(ok_rows(cl)[[sprintf("declared_cal%s_%s", al, sfx[j])]])))
    say(trade_row(sprintf("  approx (plug-in median kappa_hat), c0 = %.2f", c0_grid[j]),
                  function(cl) AX$rates[[sprintf("%s|%.2f|%s", cl, c0_grid[j], al)]]))
  }
  say(trade_row("calibrated, c0 = c2 = 1.00 (committed)",
                function(cl) { r <- R0[[cl]]$results; mean(r[[sprintf("declared_cal%s", al)]][r$status == "ok"]) }))
  say(trade_row("fixed k = 2.0, p* 0.9545, post-reduction (committed practical)",
                function(cl) fk$post[fk$cell == cl]))
  say(trade_row("fixed k = 2.0, p* 0.9545, pre-reduction", function(cl) fk$pre[fk$cell == cl]))
  say(trade_row("p* = 0.90 as executed", function(cl) mean(ok_rows(cl)$declared_conv)))
}
say("\nfixed-k rates parsed from logs/declcal_fixedk_practical.txt (k = 2.0000 rows); post-reduction column reproduced from the committed payloads in %d / %d cells.",
    sum(fk_chk), length(fk_chk))

# ---- kappa_hat_c0 configuration invariance -----------------------------------
for (al in c("05", "10")) {
  say("\n## kappa_hat_%s_c0 across configurations, grouped by n; exact (per-replicate, B 500) beside approx (40 captures, B 2000)", al)
  say("| n | c0 | cells | exact: median of per-cell medians (min-max) | exact: settable p* | approx median kappa_hat | approx settable p* |")
  say("|---|---|---|---|---|---|---|")
  for (nn in c(500L, 1000L, 1500L)) {
    cc <- names(P)[vapply(names(P), function(cl) P[[cl]]$meta$n == nn, logical(1))]
    for (j in seq_len(length(c0_grid) + 1L)) {
      col <- if (j <= length(sfx)) sprintf("kappa_hat_%s_%s", al, sfx[j]) else sprintf("kappa_hat_%s", al)
      lab <- if (j <= length(sfx)) sprintf("%.2f", c0_grid[j]) else "c2 = 1.00 (committed)"
      akey <- if (j <= length(sfx)) sprintf("%d|%.2f|%s", nn, c0_grid[j], al) else sprintf("%d|c2 (1.00, unshifted)|%s", nn, al)
      v <- vapply(cc, function(cl) median(ok_rows(cl)[[col]]), 0)
      ak <- AX$median_kappa[[akey]]
      say("| %d | %s | %s | %.3f (%.3f-%.3f) | %s | %.3f | %s |", nn, lab, paste(cc, collapse = ","),
          median(v), min(v), max(v), ptxt(median(v)), ak, ptxt(ak))
    }
  }
}
say("\n%s", FOOT)

# ---- n_admitted == 1 share, pooled ------------------------------------------
say("\n## Share of declaring replicates (cal05) with n_admitted_cal05_c0 = 1, pooled")
say("| c0 | B cells: declaring reps | B: share = 1 | C cells: declaring reps | C: share = 1 | C: median n_admitted among declaring |")
say("|---|---|---|---|---|---|")
for (j in seq_along(sfx)) {
  pool <- function(cc) { d <- do.call(c, lapply(cc, function(cl) { r <- ok_rows(cl)
    r[[paste0("n_admitted_cal05_", sfx[j])]][r[[paste0("declared_cal05_", sfx[j])]] == 1L] })); d }
  b <- pool(Bc); cx <- pool(Cc)
  say("| %.2f | %d | %s | %d | %s | %s |", c0_grid[j], length(b),
      if (length(b)) sprintf("%.4f", mean(b == 1L)) else "-", length(cx),
      if (length(cx)) sprintf("%.4f", mean(cx == 1L)) else "-",
      if (length(cx)) sprintf("%.0f", median(cx)) else "-")
}

# ---- fw_c0 against realized rates in the B cells ---------------------------
say("\n## fw_c0 (Eq. 8 at the shifted field) beside realized p* = 0.90 rates, B cells")
say("| cell | uniform HR | c0 | mean fw_1621_c0 | realized pre-family rate at the executed cutoff (max_T_pre >= z_pstar) | diff | paired MC SE | mean fw_1621_c0 | executed rate (post, rounded) | diff | paired MC SE |")
say("|---|---|---|---|---|---|---|---|---|---|---|")
for (cl in Bc) {
  r <- ok_rows(cl); hr <- if (r$dgm[1] == "null0657") 0.657 else 0.721
  pre_ex <- as.integer(r$max_T_pre >= z_exec(cl))
  for (j in seq_along(sfx)) {
    f21 <- r[[paste0("fw_1621_", sfx[j])]]
    d1 <- f21 - pre_ex; d2 <- f21 - r$declared_conv
    say("| %s | %.3f | %.2f | %.4f | %.4f | %+.4f | %.4f | %.4f | %.4f | %+.4f | %.4f |", cl, hr, c0_grid[j],
        mean(f21), mean(pre_ex), mean(d1), sd(d1) / sqrt(nrow(r)),
        mean(f21), mean(r$declared_conv), mean(d2), sd(d2) / sqrt(nrow(r)))
  }
}
# ---- the three free checks (declcal task section 7.3, per c0) ----------------
say("\n## Free check (a): Eq. 8 at the shifted field as an estimator, B cells (all c0)")
say("(the fw_c0 table above: mean fw_1621_c0 vs the as-executed rate, mean fw_1621_c0 vs the pre-family rate at the executed cutoff, paired MC SE = sd/sqrt(2000); Block A is not re-run)")
say("\n## Free check (b): how strict the calibration is -- pstar_fine at alpha = 0.05 (the settable p* at digits_fine), per c0")
say("| cell | c0 | min | 5%% | 25%% | 50%% | 75%% | 95%% | max | digits_fine range | fraction settable at the fit's digits with p* > 0.90 | fraction with no settable p* <= 1 at the fit's digits | fraction < p*(executed cutoff) |")
say("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
pr <- c(0, .05, .25, .5, .75, .95, 1)
for (cl in names(P)) { r <- ok_rows(cl)
  for (j in seq_along(sfx)) { t <- stab(r[[paste0("kappa_hat_05_", sfx[j])]]); v <- t$pstar_settable
    say("| %s | %.2f | %s | %s | %.4f | %.4f | %.4f |", cl, c0_grid[j], paste(sprintf("%.5f", stats::quantile(t$pstar_fine, pr, type = 7)), collapse = " | "),
        paste(unique(range(t$digits_fine)), collapse = "-"), mean(!is.na(v) & v > 0.90), mean(is.na(v)), mean(r[[paste0("kappa_hat_05_", sfx[j])]] < z_exec(cl))) } }
say("(last column: share of replicates whose kappa_hat_05_c0 is below the as-executed cutoff z = 1.621, i.e. the calibrated rule is looser than p* 0.90 as executed)")
say("\n## Free check (c): can the calibrated rule declare where the conventional one did not?")
say("| cell | c0 | cal05 = 1 & conv = 0 | cal10 = 1 & conv = 0 | cal05 = 1 & exact = 0 | cal10 = 1 & exact = 0 | min kappa_hat_05_c0 | min kappa_hat_10_c0 | cal05 & !conv: via family / via cutoff | cal10 & !conv: via family / via cutoff |")
say("|---|---|---|---|---|---|---|---|---|---|")
fc <- list()
for (cl in names(P)) { r <- ok_rows(cl)
  for (j in seq_along(sfx)) { s5 <- r[[paste0("declared_cal05_", sfx[j])]]; s10 <- r[[paste0("declared_cal10_", sfx[j])]]
    x <- c(sum(s5 == 1 & r$declared_conv == 0), sum(s10 == 1 & r$declared_conv == 0),
           sum(s5 == 1 & r$declared_conv_exact == 0), sum(s10 == 1 & r$declared_conv_exact == 0))
    k5 <- r[[paste0("kappa_hat_05_", sfx[j])]]; k10 <- r[[paste0("kappa_hat_10_", sfx[j])]]
    post_ge <- function(k) !is.na(r$max_T_post) & r$max_T_post >= k
    v5 <- s5 == 1 & r$declared_conv == 0; v10 <- s10 == 1 & r$declared_conv == 0
    sp <- c(sum(v5 & !post_ge(k5)), sum(v5 & post_ge(k5)), sum(v10 & !post_ge(k10)), sum(v10 & post_ge(k10)))
    fc[[paste(cl, sfx[j])]] <- c(x, sp)
    say("| %s | %.2f | %d | %d | %d | %d | %.3f | %.3f | %d / %d | %d / %d |", cl, c0_grid[j], x[1], x[2], x[3], x[4],
        min(k5), min(k10), sp[1], sp[2], sp[3], sp[4]) } }
tot <- Reduce(`+`, fc)
say("\nPooled over the 40 (cell, c0) rows: cal05 & !conv %d (via family %d, via cutoff %d); cal10 & !conv %d (via family %d, via cutoff %d).",
    tot[1], tot[5], tot[6], tot[2], tot[7], tot[8])
say("via family: max_T_post < kappa_hat (the admitting subgroup is in the pre-reduction family but not the post-reduction family the executed screen evaluated); via cutoff: max_T_post >= kappa_hat yet the rounded p* 0.90 rule did not admit, i.e. kappa_hat below the executed cutoff 1.621.")
say("(conv = declared_conv, the rounded post-reduction rule as executed; exact = declared_conv_exact, max_T_post >= 1.6449)")
writeLines(out, "logs/declcalc0_findings.txt")
saveRDS(list(status = status, gates = gates, q62 = q62, fk = fk, free_c = fc), "declcalc0_findings.rds")
