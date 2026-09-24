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
cols_res <- c("status", "dgm", "n", "rep", "declared_conv", "declared_conv_exact",
              "declared_cal05", "declared_cal10", "G_pre", "G_post", "kappa_hat_05", "kappa_hat_10",
              "pconsistency_digits", "alpha_FW_hat_1645", "alpha_FW_hat_1621", "n_band",
              "floors_id", "sg_size_argmax", "sg_size_declared", "wall_sec")
cols_aux <- c("rep", "search_declared", "replay_check", "n_unmatched")
# the executed (rounded) cutoff as the payload recorded it
z_exec <- function(cl) { z <- P[[cl]]$meta$z_round; if (is.null(z)) z <- unique(P[[cl]]$results$z_pstar)
  if (length(z) != 1L || !is.finite(z))
    stop(sprintf("%s records no single executed cutoff (meta$z_round or results$z_pstar)", cl), call. = FALSE)
  z }
cells <- list(
  declcal_bnull = c("A1", "A2", "A3"),
  declcal_inull = c("B1", "B2", "B3", "B4", "B5", "B6"),
  declcal_power = c("C1", "C2", "C3", "C4"))
# Settable p* (TASK_declcal_consumers_2026-09-24_v3), derived from the stored
# kappa_hat by the package's own inverse (.fs_decl_settable_table(), built on
# .fs_decl_settable()): pstar_settable, the smallest p* on the fit's
# pconsistency.digits grid whose rounded-screen cutoff reaches kappa_hat (NA:
# none <= 1 at those digits), and the finer pair (pstar_fine, digits_fine), the
# smallest digits at which a settable p* sits within 0.01 z of kappa_hat.
stab <- function(k, d) { d <- unique(as.integer(d)); stopifnot(length(d) == 1L, !is.na(d))
  forestsearch:::.fs_decl_settable_table(k, d) }
set_txt <- function(t) sprintf("%s ; none <= 1: %.4f ; fine %.5f [digits %s]",
  if (any(t$pstar_achievable)) med_iqr(t$pstar_settable[t$pstar_achievable], 2) else "none",
  mean(!t$pstar_achievable), median(t$pstar_fine),
  paste(unique(range(t$digits_fine)), collapse = "-"))
dgm_txt <- c(null1000 = "complete null", null0657 = "uniform benefit HR 0.657",
             null0721 = "uniform benefit HR 0.721", alt_h150 = "planted harm HR 1.5",
             alt_h200 = "planted harm HR 2.0")
P <- list(); status <- list()
for (camp in names(cells)) for (cl in cells[[camp]]) {
  f <- sprintf("../results/%s_%s_res_1_2000.rds", camp, cl)
  if (!file.exists(f)) { status[[cl]] <- "not run"; next }
  p <- readRDS(f)
  need_cols(p$results, cols_res, paste0(f, " $results")); need_cols(p$aux, cols_aux, paste0(f, " $aux"))
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
say("| cell | G_pre (median, IQR) | G_post (median, IQR) | kappa_hat_05 (median, IQR) | settable p*: median (IQR) at the fit's digits ; share none <= 1 ; median pstar_fine [digits_fine] | mean alpha_FW_hat at 1.6449 | mean alpha_FW_hat at 1.621 | median n_band |")
say("|---|---|---|---|---|---|---|---|")
for (cl in names(P)) {
  r <- P[[cl]]$results; r <- r[r$status == "ok", ]
  say("| %s | %s | %s | %s | %s | %.4f | %.4f | %d |", cl, med_iqr(r$G_pre), med_iqr(r$G_post),
      med_iqr(r$kappa_hat_05, 3), set_txt(stab(r$kappa_hat_05, r$pconsistency_digits)),
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
say("\n## 7.3 (b) How strict the calibration is: pstar_fine at alpha = 0.05 (the settable p* at digits_fine)")
say("| cell | min | 5%% | 25%% | 50%% | 75%% | 95%% | max | digits_fine range | fraction settable at the fit's digits with p* > 0.90 | fraction with no settable p* <= 1 at the fit's digits |")
say("|---|---|---|---|---|---|---|---|---|---|---|")
for (cl in names(P)) {
  r <- P[[cl]]$results; r <- r[r$status == "ok", ]
  t <- stab(r$kappa_hat_05, r$pconsistency_digits); ps <- t$pstar_settable
  q <- stats::quantile(t$pstar_fine, c(0, .05, .25, .5, .75, .95, 1), type = 1)
  say("| %s | %s | %s | %.4f | %.4f |", cl, paste(sprintf("%.5f", q), collapse = " | "),
      paste(unique(range(t$digits_fine)), collapse = "-"), mean(!is.na(ps) & ps > 0.90), mean(is.na(ps)))
}
say("\n## 7.3 (c) Calibrated declares where conventional did not")
say("| cell | cal05 = 1 & conv = 0 | cal10 = 1 & conv = 0 | reps with kappa_hat_10 < executed cutoff z_pstar |")
say("|---|---|---|---|")
for (cl in names(P)) {
  r <- P[[cl]]$results; r <- r[r$status == "ok", ]
  say("| %s | %d | %d | %d |", cl, sum(r$declared_cal05 == 1 & r$declared_conv == 0),
      sum(r$declared_cal10 == 1 & r$declared_conv == 0), sum(r$kappa_hat_10 < z_exec(cl)))
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
