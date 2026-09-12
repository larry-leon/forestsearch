# Gate 1 projection for TASK_grfmr_completion_2026-09-12: the two deferred
# 31% HR 1.75 cells and the six HR 1.00 cells.
#
# projectG.R's machinery, unchanged in method:
#   1. From the per-replicate cost DISTRIBUTION, not a mean: each cell's
#      2,000-replicate total is bootstrap-resampled (B = 4,000) from measured
#      per-replicate `fit_mr_secs`.
#   2. Overhead charged PER BATCH RENDER (3 per cell), at projectC.R's measured
#      per-block values -- not 3 x uniform.
#   3. CALIBRATED on grfmr's realized-over-projected ratio, 0.8962 (Part A, ten
#      cells, every cell 0.864-0.932: REPORT_grfmr_2026-09-11.md).  Unlike
#      projectG.R, which reported its calibration beside the projection, this
#      kickoff calibrates ON it, so the 0.896 factor is applied and the
#      uncalibrated figure is printed beside it.
#
# Cost basis per cell (stated in the `basis` column, never hidden):
#   - 12.4% HR 1.00 n 500 : MEASURED -- the null-corner probe (13.50 s median,
#     the cheapest of the five).
#   - every other HR 1.00 cell: no probe exists.  Costed from the HR 1.50 pool at
#     its own prevalence/n (n 1000 = n500/n1500 draws pooled), exactly as the
#     reference projection of ~4.9 h did.  This is conservative: at the one
#     coordinate where both exist the null probe runs 13.50 s against 14.18 s.
#     A second, null-scaled basis (HR 1.50 pool x the measured null/harm ratio
#     at 12.4% n 500) is printed beside it for the record, not used.
#   - 31% HR 1.75 n 1000 / n 1500: the HR 1.50 pool at 31%, as in projectG.R;
#     Part A showed the HR 1.75 cells' realized ratios indistinguishable from
#     HR 1.50's at matched coordinates.
#
# Ceiling 9 h wall; hard timeout 10 h (kickoff, GATE 1).
SCRATCH <- Sys.getenv("DINAMR_SCRATCH", unset = ".")
QMD_DIR <- Sys.getenv("DINAMR_QMD_DIR", unset = "..")
RES     <- file.path(QMD_DIR, "results")
CEILING_H <- 9; TIMEOUT_H <- 10
W <- 12L
CALIB <- 0.8962                                           # grfmr Part A realized / projected
OVERHEAD_BY_BLOCK <- c("12.4%" = 122.8, "31%" = 196.6)   # projectC.R, measured
RENDERS_PER_CELL  <- 3L
set.seed(8316951)
B <- 4000L

gprobe <- function(hr, n, z1q)
  file.path(RES, sprintf("grf_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d%s_nb20_grfprobe_res_1_36.rds",
                         round(100*hr), n, if (z1q) "_z1q60" else ""))
spec <- data.frame(hr = c(1.50,1.50,1.50,1.50,1.00), n = c(500L,1500L,500L,1500L,500L),
                   z1q = c(FALSE,FALSE,TRUE,TRUE,FALSE))
P <- lapply(seq_len(nrow(spec)), function(i) {
  b <- readRDS(gprobe(spec$hr[i], spec$n[i], spec$z1q[i])); r <- b$results
  list(hr = spec$hr[i], prev = if (spec$z1q[i]) "31%" else "12.4%", n = spec$n[i],
       s = r$fit_mr_secs[is.finite(r$fit_mr_secs)], det = mean(r$detected %in% 1L),
       W = as.integer(b$meta$n_workers))
})
stopifnot(length(P) == 5L, all(vapply(P, function(p) p$W, integer(1)) == W))
cat("=== GRF COST PROBES (36 replicates, 12 workers) ===\n")
print(do.call(rbind, lapply(P, function(p) data.frame(prev = p$prev, hr = p$hr, n = p$n,
  selection = round(p$det, 4), s_q10 = unname(quantile(p$s,.1)), s_med = median(p$s),
  s_q90 = unname(quantile(p$s,.9)), s_max = max(p$s)))), row.names = FALSE, digits = 4)

get_s <- function(prev, hr, n) {
  i <- which(vapply(P, function(p) p$prev == prev && p$hr == hr && p$n == n, logical(1)))
  if (length(i)) P[[i]]$s else NULL }
pool150 <- function(prev, n) if (n == 1000L) c(get_s(prev, 1.50, 500L), get_s(prev, 1.50, 1500L)) else get_s(prev, 1.50, n)
null_ratio <- median(get_s("12.4%", 1.00, 500L)) / median(get_s("12.4%", 1.50, 500L))

cells <- data.frame(order = 1:8,
  prev = c("31%","31%","12.4%","12.4%","12.4%","31%","31%","31%"),
  hr   = c(1.75, 1.75, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00),
  n    = c(1500L, 1000L, 500L, 1000L, 1500L, 500L, 1000L, 1500L),
  tag  = c("A31_h175_n1500","A31_h175_n1000","C124_h100_n500","C124_h100_n1000","C124_h100_n1500",
           "C31_h100_n500","C31_h100_n1000","C31_h100_n1500"), stringsAsFactors = FALSE)
boot_h <- function(s, ovh) { tot <- replicate(B, sum(sample(s, 2000L, replace = TRUE)))
  (tot / W + RENDERS_PER_CELL * ovh) / 3600 }
G <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
  pv <- cells$prev[i]; hr <- cells$hr[i]; n <- cells$n[i]; ovh <- unname(OVERHEAD_BY_BLOCK[pv])
  meas <- get_s(pv, hr, n)
  if (!is.null(meas)) { s <- meas; basis <- "measured (null-corner probe)" }
  else { s <- pool150(pv, n)
         basis <- paste0("HR 1.50 pool", if (n == 1000L) ", n500/n1500 pooled" else "",
                         if (hr == 1.75) " (HR 1.75 unprobed)" else " (HR 1.00 unprobed here)") }
  h <- boot_h(s, ovh)
  alt <- if (hr == 1.00 && is.null(meas)) median(boot_h(s * null_ratio, ovh)) * CALIB else NA_real_
  data.frame(order = i, cell = cells$tag[i], prev = pv, hr = hr, n = n, basis = basis, ovh_s = ovh,
             h_uncal = median(h), h_q05 = unname(quantile(h,.05)) * CALIB, h_med = median(h) * CALIB,
             h_q95 = unname(quantile(h,.95)) * CALIB, h_null_scaled = alt, stringsAsFactors = FALSE)
}))
G$cum_h <- cumsum(G$h_med)
cat(sprintf("\nnull/harm cost ratio at 12.4%% n 500 (medians): %.4f\n", null_ratio))
cat(sprintf("\n=== THE EIGHT CELLS, 2,000 replicates (bootstrap B = %d, %d workers), CALIBRATED x %.4f ===\n", B, W, CALIB))
cat("h_uncal = uncalibrated median (compute + 3 x ovh_s); h_q05/h_med/h_q95 = calibrated.\n")
cat("h_null_scaled = the alternative, null-scaled basis for the unprobed null cells -- reported, NOT used.\n")
print(G, row.names = FALSE, digits = 4)

tot <- sum(G$h_med)
cat("\n=== GATE 1 TOTAL ===\n")
cat(sprintf("two harm cells : %.3f h   (reference ~1.9 h)\n", sum(G$h_med[G$hr == 1.75])))
cat(sprintf("six null cells : %.3f h   (reference ~4.9 h)\n", sum(G$h_med[G$hr == 1.00])))
cat(sprintf("eight cells    : %.3f h calibrated (90%% band %.3f-%.3f h); uncalibrated %.3f h   (reference ~6.8 h)\n",
            tot, sum(G$h_q05), sum(G$h_q95), sum(G$h_uncal)))
cat(sprintf("ceiling %.0f h wall; hard timeout %.0f h (kickoff, GATE 1)\n", CEILING_H, TIMEOUT_H))
cat(sprintf("even uncalibrated the eight read %.3f h: headroom %.3f h under the ceiling\n",
            sum(G$h_uncal), CEILING_H - sum(G$h_uncal)))
fit <- G$cum_h <= CEILING_H; k <- sum(fit)
cat(sprintf("cells that fit under the ceiling in run order: %d of 8 (cumulative %.3f h)\n", k, if (k) G$cum_h[k] else 0))
if (k < 8L) { cat("DEFERRED, from the tail:\n"); print(G[!fit, c("order","cell","h_med")], row.names = FALSE) } else cat("none deferred.\n")
cat(sprintf("\nGATE 1: %s -- run %d, defer %d.  Replicate count NOT reduced.\n", if (k) "GO" else "NO-GO", k, 8L - k))
saveRDS(list(cells = G, total_h = tot, total_uncal_h = sum(G$h_uncal), calib = CALIB,
             ceiling_h = CEILING_H, timeout_h = TIMEOUT_H, n_fit = k, null_ratio = null_ratio,
             overhead_by_block = OVERHEAD_BY_BLOCK), file.path(SCRATCH, "projectionGC.rds"))
