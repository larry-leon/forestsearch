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

AX <- readRDS("../results/declcal_c0approx_res.rds")$approx   # committed at e9de87c7
ptxt <- function(k) sprintf("%.5f", 2 * pnorm(k) - 1)
FOOT <- paste0("Implied p* = 2*pnorm(kappa) - 1. A literal re-run of the fixed-p* screen at that p* needs ",
               "pconsistency.digits >= 5 (the implied p* are quoted to 5 decimals). At the default ",
               "pconsistency.digits = 2 the consistency proportion is rounded to 2 decimals before it is ",
               "compared with p*, so e.g. p* = 0.99586 acts as 'rounded proportion = 1.00' and the executed ",
               "cutoff is not kappa_hat.")
P <- list(); R0 <- list(); status <- list()
for (blk in names(cells)) for (cl in cells[[blk]]) {
  f <- sprintf("../results/declcalc0_%s_%s_res_1_2000.rds", blk, cl)
  R0[[cl]] <- readRDS(sprintf("../results/declcal_%s_%s_res_1_2000.rds", blk, cl))
  if (!file.exists(f)) { status[[cl]] <- "not run"; next }
  p <- readRDS(f)
  st <- p$meta$cell_status; fin <- isTRUE(p$meta$final)
  status[[cl]] <- sprintf("%s%s ; %d rows ; %s ; elapsed %.0f s ; median wall %.2f s ; max wall %.1f s",
                          st, if (fin) "" else " (not final)", nrow(p$results),
                          paste(names(table(p$results$status)), table(p$results$status), collapse = ", "),
                          p$meta$elapsed_s, median(p$results$wall_sec), max(p$results$wall_sec))
  if (identical(st, "complete") && fin) P[[cl]] <- p
}
ok_rows <- function(cl) { r <- P[[cl]]$results; r[r$status == "ok", ] }
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
    say("| %s approx | plug-in median kappa_hat(c0) at n, max_T_pre | %d | %.4f | %s | %.4f |", cl, r$n[1],
        AX$rates[[sprintf("%s|conv|%s", cl, al)]],
        paste(sprintf("%.4f", AX$rates[sprintf("%s|%.2f|%s", cl, c0_grid, al)]), collapse = " | "),
        AX$rates[[sprintf("%s|c2|%s", cl, al)]])
  }
  say("\n\"approx\" rows: REPORT_declcal_c0_approx_2026-09-22 (e9de87c7) -- each rate is mean(max_T_pre >= median kappa_hat(c0) over 40 captures at that n, B 2000), a fixed cutoff, not the per-replicate calibrated rule of the exact row above it.")
}

# ---- 6.2 calibration quantities per c0 ---------------------------------------
say("\n## 6.2 The calibration's own quantities per c0")
say("| cell | c0 | kappa_hat_05 (median, IQR) | implied p* (median, IQR) | mean fw_1645 | mean fw_1621 | n_admitted_cal05 median (all reps) | declaring reps (cal05) | share of declaring reps with n_admitted_cal05 = 1 |")
say("|---|---|---|---|---|---|---|---|---|")
q62 <- list()
for (cl in names(P)) {
  r <- ok_rows(cl)
  for (j in seq_along(sfx)) {
    s <- sfx[j]; dec <- r[[paste0("declared_cal05_", s)]] == 1L
    na <- r[[paste0("n_admitted_cal05_", s)]]
    sh1 <- if (any(dec)) sprintf("%.4f (%d / %d)", mean(na[dec] == 1L), sum(na[dec] == 1L), sum(dec)) else "no declarations"
    say("| %s | %.2f | %s | %s | %.4f | %.4f | %d | %d | %s |", cl, c0_grid[j],
        med_iqr(r[[paste0("kappa_hat_05_", s)]], 3), med_iqr(r[[paste0("pstar_implied_05_", s)]], 5),
        mean(r[[paste0("fw_1645_", s)]]), mean(r[[paste0("fw_1621_", s)]]),
        as.integer(median(na)), sum(dec), sh1)
    q62[[paste(cl, s)]] <- data.frame(cell = cl, n = r$n[1], c0 = c0_grid[j],
      kappa_med = median(r[[paste0("kappa_hat_05_", s)]]), n_dec = sum(dec),
      n_dec_adm1 = sum(na[dec] == 1L), n_dec_adm_gt1 = sum(na[dec] > 1L))
  }
  say("| %s | c2 = 1.00 (committed) | %s | %s | %.4f | %.4f | %d | %d | %s |", cl,
      med_iqr(r$kappa_hat_05, 3), med_iqr(r$pstar_implied_05, 5), mean(r$alpha_FW_hat_1645),
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
  say("| n | c0 | cells | exact: median of per-cell medians (min-max) | exact implied p* | approx median kappa_hat | approx implied p* |")
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

# ---- fw_1645_c0 against realized rates in the B cells ------------------------
say("\n## fw_c0 (Eq. 8 at the shifted field) beside realized p* = 0.90 rates, B cells")
say("| cell | uniform HR | c0 | mean fw_1645_c0 | realized pre-family exact rate (max_T_pre >= 1.6449) | diff | paired MC SE | mean fw_1621_c0 | executed rate (post, rounded) | diff | paired MC SE |")
say("|---|---|---|---|---|---|---|---|---|---|---|")
for (cl in Bc) {
  r <- ok_rows(cl); hr <- if (r$dgm[1] == "null0657") 0.657 else 0.721
  pre_ex <- as.integer(r$max_T_pre >= qnorm(0.95))
  for (j in seq_along(sfx)) {
    f45 <- r[[paste0("fw_1645_", sfx[j])]]; f21 <- r[[paste0("fw_1621_", sfx[j])]]
    d1 <- f45 - pre_ex; d2 <- f21 - r$declared_conv
    say("| %s | %.3f | %.2f | %.4f | %.4f | %+.4f | %.4f | %.4f | %.4f | %+.4f | %.4f |", cl, hr, c0_grid[j],
        mean(f45), mean(pre_ex), mean(d1), sd(d1) / sqrt(nrow(r)),
        mean(f21), mean(r$declared_conv), mean(d2), sd(d2) / sqrt(nrow(r)))
  }
}
writeLines(out, "logs/declcalc0_findings.txt")
saveRDS(list(status = status, gates = gates, q62 = q62, fk = fk), "declcalc0_findings.rds")
