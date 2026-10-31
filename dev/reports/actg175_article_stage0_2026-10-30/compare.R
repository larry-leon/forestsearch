# compare.R -- gates A/B/C and the recorded comparisons, from rungs/<rung>_result.rds.
# Usage: Rscript dev/reports/actg175_article_stage0_2026-10-30/compare.R <rung>
REC <- "dev/reports/actg175_article_stage0_2026-10-30"
rung <- commandArgs(TRUE)[1]
res <- readRDS(file.path(REC, "rungs", paste0(rung, "_result.rds")))
x <- res$cells
ex <- read.csv(file.path(REC, "manuscript_source_2026b", "applications_metrics.csv"))
out <- list()
row <- function(cell, ours, printed, dp, payload = NA_real_) {
  ok <- if (is.na(printed)) NA else identical(sprintf("%.*f", dp, ours), sprintf("%.*f", dp, printed))
  out[[length(out) + 1L]] <<- data.frame(cell = cell, ours = sprintf("%.10g", ours),
    ours_printed = sprintf("%.*f", dp, ours), table = if (is.na(printed)) "" else sprintf("%.*f", dp, printed),
    match = ok, payload = sprintf("%.10g", payload),
    abs_diff_payload = if (is.na(payload)) NA_real_ else abs(ours - payload),
    all_equal_1e8 = if (is.na(payload)) NA else isTRUE(all.equal(ours, payload, tolerance = 1e-8)))
}
if (rung %in% c("M-bin", "M-cont")) {
  bin <- rung == "M-bin"; dp <- if (bin) 3 else 1
  an <- if (bin) "actg175_binary" else "actg175_continuous"
  H <- ex[ex$analysis == an & ex$region == "H", ]; Hc <- ex[ex$analysis == an & ex$region == "Hc", ]
  T16 <- if (bin) list(whole = c(0.665, 0.520, 0.850), N = 72, ev = 37, naive = c(3.584, 1.315, 9.774),
                       mr = c(1.773, 0.604, 5.206), field = 0.686, cN = 1011, cev = 393,
                       cnaive = c(0.587, 0.454, 0.759), cmr = c(0.620, 0.372, 1.033), fs = 0.786,
                       bonf = c(0.597, 0.817), sb = 0.704, fam = 2794, wN = 1083, wev = 430)
         else list(whole = c(-27.6, -43.1, -12.1), N = 79, ev = NA, naive = c(71.1, -5.1, 147.3),
                   mr = c(22.1, -77.5, 121.6), field = -38.3, cN = 1004, cev = NA,
                   cnaive = c(-36.2, -51.4, -20.9), cmr = c(-32.7, -62.2, -3.2), fs = -17.6,
                   bonf = c(-51.2, -15.8), sb = 47.9, fam = 4935, wN = 1083, wev = NA)
  row("whole N", x$whole_n, T16$wN, 0, H$n_total)
  if (bin) row("whole events", x$whole_events, T16$wev, 0, H$events_total)
  row("whole est", x$whole_est, T16$whole[1], dp, H$whole_trial_est)
  row("whole lo", x$whole_lo, T16$whole[2], dp, H$whole_trial_lo)
  row("whole hi", x$whole_hi, T16$whole[3], dp, H$whole_trial_hi)
  row("H N", x$N_H, T16$N, 0, H$n)
  if (bin) row("H events", x$events_H, T16$ev, 0, H$events)
  row("H unadj est", x$naive_est, T16$naive[1], dp, H$naive_est)
  row("H unadj lo", x$naive_lo, T16$naive[2], dp, H$naive_lo)
  row("H unadj hi", x$naive_hi, T16$naive[3], dp, H$naive_hi)
  row("H corrected est", x$mr_est, T16$mr[1], dp, H$mr_est)
  row("H corrected lo", x$mr_lo, T16$mr[2], dp, H$mr_lo)
  row("H corrected hi", x$mr_hi, T16$mr[3], dp, H$mr_hi)
  row("H field lower", x$field_lower_H, T16$field, dp, H$field_bound)
  row("Hc N", x$N_Hc, T16$cN, 0, Hc$n)
  if (bin) row("Hc events", x$events_Hc, T16$cev, 0, Hc$events)
  row("Hc unadj est", x$c_naive_est, T16$cnaive[1], dp, Hc$naive_est)
  row("Hc unadj lo", x$c_naive_lo, T16$cnaive[2], dp, Hc$naive_lo)
  row("Hc unadj hi", x$c_naive_hi, T16$cnaive[3], dp, Hc$naive_hi)
  row("Hc corrected est", x$c_mr_est, T16$cmr[1], dp, Hc$mr_est)
  row("Hc corrected lo", x$c_mr_lo, T16$cmr[2], dp, Hc$mr_lo)
  row("Hc corrected hi", x$c_mr_hi, T16$cmr[3], dp, Hc$mr_hi)
  row("Hc field-s upper", x$field_s_upper_Hc, T16$fs, dp, Hc$field_bound)
  row("Bonferroni H lower", x$bonf_lower_H, T16$bonf[1], dp, H$bonf_bound)
  row("Bonferroni Hc upper", x$bonf_upper_Hc, T16$bonf[2], dp, Hc$bonf_bound)
  row("selection bias", x$selection_bias, T16$sb, dp, H$selection_bias)
  row("re-selection family", x$n_family, T16$fam, 0, H$n_family)
  # region by membership: the paper's region, evaluated on the data the rung built
  dfp <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L)); dfp <- dfp[!is.na(dfp$cd420), ]
  ref <- if (bin) dfp$wtkg > 86 & dfp$cd40 > 380 else dfp$cd40 > 507 & dfp$gender == 1
  cat("REGION label (this fit):", x$region, "\n")
  cat("REGION label (extract): ", H$region_def, "\n")
  cat("REGION membership identical to Table S16's region:", identical(unname(res$membership), ref), "\n")
}
if (grepl("^M-S3-", rung)) {
  m <- sub("^M-S3-", "", rung)
  P <- readRDS(file.path(REC, "manuscript_source_2026a", "analysis_actg175_binary_multimethod_fixed_family_payload.rds"))
  H <- P$table[P$table$method == m & P$table$region == "H", ]; Hc <- P$table[P$table$method == m & P$table$region == "Hc", ]
  S3 <- list(FS = list(n = 72, nv = c(3.58, 1.31, 9.77), mr = c(1.79, 0.63, 5.08), cnv = c(0.59, 0.45, 0.76), cmr = c(0.62, 0.37, 1.03)),
             DINA = list(n = 98, nv = c(1.94, 0.84, 4.47), mr = c(0.95, 0.29, 3.10), cnv = c(0.60, 0.46, 0.78), cmr = c(0.63, 0.38, 1.06)),
             GRF = list(n = 95, nv = c(2.25, 0.97, 5.20), mr = c(1.10, 0.43, 2.80), cnv = c(0.58, 0.45, 0.76), cmr = c(0.62, 0.37, 1.04)))[[m]]
  row("H n", x$n_H, S3$n, 0, H$n)
  row("H naive est", x$naive_est, S3$nv[1], 2, H$naive_est); row("H naive lo", x$naive_lo, S3$nv[2], 2, H$naive_lo)
  row("H naive hi", x$naive_hi, S3$nv[3], 2, H$naive_hi)
  row("H MR est", x$mr_est, S3$mr[1], 2, H$mr_est); row("H MR lo", x$mr_lo, S3$mr[2], 2, H$mr_lo)
  row("H MR hi", x$mr_hi, S3$mr[3], 2, H$mr_hi)
  row("Hc n", x$n_Hc, 1083 - S3$n, 0, Hc$n)
  row("Hc naive est", x$c_naive_est, S3$cnv[1], 2, Hc$naive_est); row("Hc naive lo", x$c_naive_lo, S3$cnv[2], 2, Hc$naive_lo)
  row("Hc naive hi", x$c_naive_hi, S3$cnv[3], 2, Hc$naive_hi)
  row("Hc MR est", x$c_mr_est, S3$cmr[1], 2, Hc$mr_est); row("Hc MR lo", x$c_mr_lo, S3$cmr[2], 2, Hc$mr_lo)
  row("Hc MR hi", x$c_mr_hi, S3$cmr[3], 2, Hc$mr_hi)
  row("re-selection family", x$n_family, NA, 0, P$meta$timings[[tolower(m)]]$n_family)
  row("selection bias", x$selection_bias, NA, 3, P$meta$timings[[tolower(m)]]$selection_bias)
  cat("REGION label (this fit):", x$region, "\n")
  cat("REGION label (payload): ", P$labels[[tolower(m)]], "\n")
  dfp <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L)); dfp <- dfp[!is.na(dfp$cd420), ]; dfp$id <- seq_len(nrow(dfp))
  ref <- switch(m, FS = dfp$wtkg > 86 & dfp$cd40 > 380, DINA = dfp$preanti >= 849.4 & dfp$cd40 >= 338,
                GRF = dfp$wtkg > 84.3696 & dfp$cd40 > 368.2)
  cat("REGION membership identical to Table S3's region:", identical(sort(dfp$id[ref]), sort(res$membership_ids)), "\n")
}
tab <- do.call(rbind, out)
options(width = 220); print(tab, row.names = FALSE, right = FALSE)
cat("\nPRINTED-PRECISION MATCH:", sum(tab$match, na.rm = TRUE), "/", sum(!is.na(tab$match)),
    "| payload all.equal(1e-8):", sum(tab$all_equal_1e8, na.rm = TRUE), "/", sum(!is.na(tab$all_equal_1e8)),
    "| max abs diff vs payload:", signif(max(tab$abs_diff_payload, na.rm = TRUE), 4), "\n")
write.csv(tab, file.path(REC, "rungs", paste0(rung, "_compare.csv")), row.names = FALSE)
