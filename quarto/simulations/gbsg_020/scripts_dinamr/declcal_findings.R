# declcal_findings.R -- Stage 2 tables for REPORT_declaration_calibration_evaluation
# (TASK_declcal_CAMPAIGN_2026-09-22_v2, section 7).  Reads the committed payloads
# results/declcal_{bnull,inull,power}_<cell>_res_1_2000.rds; writes
# logs/declcal_findings.txt (markdown tables) and declcal_findings.rds.
# A cell whose payload is missing, incomplete or tripped the per-cell gate
# contributes no rate (task section 9).
wilson <- function(x, n, z = qnorm(0.975)) {
  p <- x / n; d <- 1 + z^2 / n
  c((p + z^2 / (2 * n) - z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / d,
    (p + z^2 / (2 * n) + z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))) / d)
}
fmt_w <- function(x) { x <- x[!is.na(x)]; n <- length(x); k <- sum(x); w <- wilson(k, n)
  sprintf("%.4f [%.4f, %.4f]", k / n, w[1], w[2]) }
med_iqr <- function(x, d = 0) { q <- stats::quantile(x, c(.5, .25, .75), na.rm = TRUE)
  sprintf(paste0("%.", d, "f (%.", d, "f-%.", d, "f)"), q[1], q[2], q[3]) }
cells <- list(
  declcal_bnull = c("A1", "A2", "A3"),
  declcal_inull = c("B1", "B2", "B3", "B4", "B5", "B6"),
  declcal_power = c("C1", "C2", "C3", "C4"))
dgm_txt <- c(null1000 = "complete null", null0657 = "uniform benefit HR 0.657",
             null0721 = "uniform benefit HR 0.721", alt_h150 = "planted harm HR 1.5",
             alt_h200 = "planted harm HR 2.0")
P <- list(); status <- list()
for (camp in names(cells)) for (cl in cells[[camp]]) {
  f <- sprintf("../results/%s_%s_res_1_2000.rds", camp, cl)
  if (!file.exists(f)) { status[[cl]] <- "not run"; next }
  p <- readRDS(f)
  st <- p$meta$cell_status; fin <- isTRUE(p$meta$final)
  status[[cl]] <- sprintf("%s%s ; %d rows ; %s ; elapsed %.0f s", st, if (fin) "" else " (not final)",
                          nrow(p$results), paste(names(table(p$results$status)), table(p$results$status), collapse = ", "),
                          p$meta$elapsed_s)
  if (identical(st, "complete") && fin) P[[cl]] <- p
}
out <- character(0); say <- function(...) { s <- sprintf(...); cat(s, "\n"); out <<- c(out, s) }
say("## Cell status"); for (cl in names(status)) say("- %s: %s", cl, status[[cl]])

say("\n## 7.1 Primary table")
say("| cell | DGM | n | conventional (as executed) | conventional (exact z) | calibrated alpha = 0.05 | calibrated alpha = 0.10 |")
say("|---|---|---|---|---|---|---|")
tab <- list()
for (cl in names(P)) {
  r <- P[[cl]]$results; r <- r[r$status == "ok", ]
  say("| %s | %s | %d | %s | %s | %s | %s |", cl, dgm_txt[[r$dgm[1]]], r$n[1],
      fmt_w(r$declared_conv), fmt_w(r$declared_conv_exact), fmt_w(r$declared_cal05), fmt_w(r$declared_cal10))
  tab[[cl]] <- data.frame(cell = cl, n_ok = nrow(r), conv = mean(r$declared_conv),
                          exact = mean(r$declared_conv_exact), cal05 = mean(r$declared_cal05),
                          cal10 = mean(r$declared_cal10))
}

say("\n## 7.2 The calibration's own quantities")
say("| cell | G_pre (median, IQR) | G_post (median, IQR) | kappa_hat_05 (median, IQR) | implied p* (median, IQR) | mean alpha_FW_hat at 1.6449 | mean alpha_FW_hat at 1.621 | median n_band |")
say("|---|---|---|---|---|---|---|---|")
for (cl in names(P)) {
  r <- P[[cl]]$results; r <- r[r$status == "ok", ]
  say("| %s | %s | %s | %s | %s | %.4f | %.4f | %d |", cl, med_iqr(r$G_pre), med_iqr(r$G_post),
      med_iqr(r$kappa_hat_05, 3), med_iqr(r$pstar_implied_05, 5),
      mean(r$alpha_FW_hat_1645), mean(r$alpha_FW_hat_1621), as.integer(median(r$n_band)))
}
say("\nn_band distribution (candidates of the post-reduction family with closed-form rate in [0.895, 0.900)):")
say("| cell | min | q25 | median | q75 | q95 | max | share of reps with n_band > 0 | reps where rounding alone declared (conv 1, exact 0) |")
say("|---|---|---|---|---|---|---|---|---|")
for (cl in names(P)) {
  r <- P[[cl]]$results; r <- r[r$status == "ok", ]
  q <- stats::quantile(r$n_band, c(0, .25, .5, .75, .95, 1), type = 1)
  say("| %s | %s | %.4f | %d |", cl, paste(q, collapse = " | "), mean(r$n_band > 0),
      sum(r$declared_conv == 1 & r$declared_conv_exact == 0))
}

say("\n## 7.3 (a) Eq. 8 as an estimator, Block A")
say("| cell | mean alpha_FW_hat_1621 | as-executed rate | diff | paired MC SE | mean alpha_FW_hat_1645 | exact-z rate | diff | paired MC SE |")
say("|---|---|---|---|---|---|---|---|---|")
for (cl in intersect(c("A1", "A2", "A3"), names(P))) {
  r <- P[[cl]]$results; r <- r[r$status == "ok", ]
  d1 <- r$alpha_FW_hat_1621 - r$declared_conv; d2 <- r$alpha_FW_hat_1645 - r$declared_conv_exact
  say("| %s | %.4f | %.4f | %+.4f | %.4f | %.4f | %.4f | %+.4f | %.4f |", cl,
      mean(r$alpha_FW_hat_1621), mean(r$declared_conv), mean(d1), sd(d1) / sqrt(nrow(r)),
      mean(r$alpha_FW_hat_1645), mean(r$declared_conv_exact), mean(d2), sd(d2) / sqrt(nrow(r)))
}
say("\n## 7.3 (b) How strict the calibration is: pstar_implied_05")
say("| cell | min | 5%% | 25%% | 50%% | 75%% | 95%% | max | fraction > 0.90 |")
say("|---|---|---|---|---|---|---|---|---|")
for (cl in names(P)) {
  r <- P[[cl]]$results; r <- r[r$status == "ok", ]
  q <- stats::quantile(r$pstar_implied_05, c(0, .05, .25, .5, .75, .95, 1), type = 1)
  say("| %s | %s | %.4f |", cl, paste(sprintf("%.5f", q), collapse = " | "), mean(r$pstar_implied_05 > 0.90))
}
say("\n## 7.3 (c) Calibrated declares where conventional did not")
say("| cell | cal05 = 1 & conv = 0 | cal10 = 1 & conv = 0 | reps with kappa_hat_10 < 1.6449 |")
say("|---|---|---|---|")
for (cl in names(P)) {
  r <- P[[cl]]$results; r <- r[r$status == "ok", ]
  say("| %s | %d | %d | %d |", cl, sum(r$declared_cal05 == 1 & r$declared_conv == 0),
      sum(r$declared_cal10 == 1 & r$declared_conv == 0), sum(r$kappa_hat_10 < qnorm(0.95)))
}
say("\n## Diagnostics: candidate size behind declarations (sg_size_argmax = size of the max-T candidate)")
say("| cell | n.min | sg_size_argmax median (IQR), all reps | ... on cal05 declarations | sg_size_declared median (IQR), conv declarations | mean wall_sec |")
say("|---|---|---|---|---|---|")
for (cl in names(P)) {
  r <- P[[cl]]$results; r <- r[r$status == "ok", ]
  say("| %s | %s | %s | %s | %s | %.2f |", cl, sub("_.*", "", sub("nmin", "", r$floors_id[1])),
      med_iqr(r$sg_size_argmax), if (any(r$declared_cal05 == 1)) med_iqr(r$sg_size_argmax[r$declared_cal05 == 1]) else "none",
      med_iqr(r$sg_size_declared), mean(r$wall_sec))
}
say("\n## Fidelity gate (declared_conv vs the search's own indicator, every ok replicate)")
for (cl in names(P)) {
  r <- P[[cl]]$results; a <- P[[cl]]$aux[match(r$rep, P[[cl]]$aux$rep), ]
  ok <- r$status == "ok"
  say("- %s: %d of %d agree ; replay_check TRUE on %d ; n_unmatched total %d ; floors %s ; digits %s",
      cl, sum(r$declared_conv[ok] == a$search_declared[ok]), sum(ok), sum(a$replay_check[ok] %in% TRUE),
      sum(a$n_unmatched[ok], na.rm = TRUE), paste(unique(r$floors_id), collapse = ","),
      paste(unique(r$pconsistency_digits), collapse = ","))
}
writeLines(out, "logs/declcal_findings.txt")
saveRDS(list(tab = do.call(rbind, tab), status = status), "declcal_findings.rds")
