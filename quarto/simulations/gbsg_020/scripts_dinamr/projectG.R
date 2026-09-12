# Gate 1 projection for TASK_grfmr_campaign_2026-09-11, PART A (twelve harm cells).
#
# Same machinery as projectC.R, repointed at the GRF probes:
#   1. Projects from the per-replicate cost DISTRIBUTION, not a mean: each
#      cell's 2,000-replicate total is bootstrap-resampled from the measured
#      per-replicate `fit_mr_secs` of the relevant probe, so every cell figure
#      carries a 90% band.
#   2. Charges the MEASURED per-render overhead per batch render (3 renders per
#      cell: batch 1, batch 1001, combine), not "3 x uniform".  The per-block
#      overhead constants are projectC.R's, measured on the eleven completed
#      dinamr cells by walls.R.
#   3. Calibrates on the dinamr realized-over-projected record: Block A 1.27
#      before the overhead correction and 0.949 after; Block B 0.92.  The
#      correction that produced 0.949 is the one applied here, so NO blanket
#      multiplier is used; the calibration is reported as a check beside the
#      projection, not folded into it.
#
# THE GRF COST PROFILE IS THE OPPOSITE OF DINA'S.  The probes measured
# 13.5-19.8 s median per replicate, rising with n and with prevalence and FLAT
# in family size (|rho| <= 0.171).  So cost is interpolated ACROSS n within a
# prevalence -- never across family size, which carries no cost signal here.
#
# Ceiling 9 h wall; hard timeout 10 h (kickoff, "GATE 1 - compute go/no-go";
# Larry's available window is 9-10 h).
#
# HR 1.75 is not probed.  The probes are at HR 1.50 and HR 1.00 only, so the
# HR 1.75 cells are costed at their own prevalence/n from the HR 1.50 pool and
# that assumption is stated in the output rather than hidden in a constant.
# This is exactly what the original dinamr Gate 1 did for its HR 1.75 cells.
SCRATCH <- Sys.getenv("DINAMR_SCRATCH", unset = ".")
QMD_DIR <- Sys.getenv("DINAMR_QMD_DIR", unset = "..")
RES     <- file.path(QMD_DIR, "results")
CEILING_H <- 9; TIMEOUT_H <- 10
W <- 12L
OVERHEAD_BY_BLOCK <- c("12.4%" = 122.8, "31%" = 196.6)   # projectC.R, measured
RENDERS_PER_CELL  <- 3L                                   # batch 1, batch 1001, combine
set.seed(8316951)
B <- 4000L

gprobe <- function(hr, n, z1q)
  file.path(RES, sprintf("grf_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d%s_nb20_grfprobe_res_1_36.rds",
                         round(100*hr), n, if (z1q) "_z1q60" else ""))

## ---- 1. the GRF probes, as MEASURED ----------------------------------------
spec <- data.frame(hr = c(1.50,1.50,1.50,1.50,1.00),
                   n  = c(500L,1500L,500L,1500L,500L),
                   z1q= c(FALSE,FALSE,TRUE,TRUE,FALSE))
P <- lapply(seq_len(nrow(spec)), function(i) {
  f <- gprobe(spec$hr[i], spec$n[i], spec$z1q[i])
  if (!file.exists(f)) return(NULL)
  b <- readRDS(f); r <- b$results
  list(hr = spec$hr[i], prev = if (spec$z1q[i]) "31%" else "12.4%", n = spec$n[i],
       s = r$fit_mr_secs[is.finite(r$fit_mr_secs)],
       K = r$n_family[is.finite(r$n_family)],
       # rho(s, K) is computed on the rows where BOTH are finite: n_family is
       # NA on an undetected replicate (the HR 1.00 probe has one), so the two
       # vectors above are not always the same length.
       sK = { ii <- is.finite(r$fit_mr_secs) & is.finite(r$n_family)
              data.frame(s = r$fit_mr_secs[ii], K = r$n_family[ii]) },
       det = mean(r$detected %in% 1L), reps = nrow(r),
       W = as.integer(b$meta$n_workers))
})
P <- Filter(Negate(is.null), P)
stopifnot(length(P) == 5L, all(vapply(P, function(p) p$W, integer(1)) == W))
cat("=== GRF COST PROBES (36 replicates each, 12 workers; TASK_dinamr_blockC_grfprobe Part B) ===\n")
tab <- do.call(rbind, lapply(P, function(p) data.frame(
  prev = p$prev, hr = p$hr, n = p$n, reps = p$reps, detection = p$det,
  K_med = median(p$K), K_min = min(p$K), K_max = max(p$K),
  s_q10 = unname(quantile(p$s,.10)), s_med = median(p$s),
  s_q90 = unname(quantile(p$s,.90)), s_max = max(p$s), s_mean = mean(p$s))))
print(tab, row.names = FALSE, digits = 4)
cat("\nwall against family size, per probe (the DINA stratifier carries no cost signal here):\n")
for (p in P) cat(sprintf("  %-6s HR %.2f n %-5d : Pearson rho(s, K) = %+0.3f   Spearman = %+0.3f   (%d rows both finite)\n",
  p$prev, p$hr, p$n, suppressWarnings(cor(p$sK$s, p$sK$K)),
  suppressWarnings(cor(p$sK$s, p$sK$K, method = "spearman")), nrow(p$sK)))

## ---- 2. the per-replicate cost pool for each of the twelve cells ------------
# Within a prevalence: n 500 and n 1500 are measured at HR 1.50; n 1000 mixes
# the two pools 50/50, interpolating the DISTRIBUTION and not only its mean
# (projectC.R's device).  HR 1.75 reuses its prevalence/n pool from HR 1.50.
pool_of <- function(prev, n) {
  ps  <- Filter(function(p) p$prev == prev && p$hr == 1.50, P)
  s5  <- ps[[which(vapply(ps, function(p) p$n, integer(1)) == 500L)]]$s
  s15 <- ps[[which(vapply(ps, function(p) p$n, integer(1)) == 1500L)]]$s
  if (n == 500L) s5 else if (n == 1500L) s15 else c(s5, s15)
}
basis_of <- function(hr, n) paste0(
  if (n == 1000L) "n500/n1500 pooled draws" else "measured",
  if (hr != 1.50) "; HR 1.50 pool (HR 1.75 unprobed)" else "")

boot_cell <- function(s, ovh) {
  tot <- replicate(B, sum(sample(s, 2000L, replace = TRUE)))
  (tot / W + RENDERS_PER_CELL * ovh) / 3600
}

# Run order of the kickoff; defer from the TAIL.
cells <- data.frame(
  order = 1:12,
  prev  = rep(c("12.4%","31%","12.4%","31%"), each = 3),
  hr    = rep(c(1.50,1.50,1.75,1.75), each = 3),
  n     = rep(c(500L,1000L,1500L), times = 4),
  stringsAsFactors = FALSE)
G <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
  ovh <- unname(OVERHEAD_BY_BLOCK[cells$prev[i]])
  h <- boot_cell(pool_of(cells$prev[i], cells$n[i]), ovh)
  data.frame(order = cells$order[i], prev = cells$prev[i], hr = cells$hr[i],
             n = cells$n[i], basis = basis_of(cells$hr[i], cells$n[i]),
             ovh_s = ovh, h_q05 = unname(quantile(h,.05)), h_med = median(h),
             h_q95 = unname(quantile(h,.95)), stringsAsFactors = FALSE)
}))
G$cum_h <- cumsum(G$h_med)
cat(sprintf("\n=== PART A, twelve cells at 2,000 replicates (bootstrap B = %d, %d workers) ===\n", B, W))
cat(sprintf("h_med / h_q05 / h_q95 already carry the measured per-render overhead (%d renders x ovh_s).\n",
            RENDERS_PER_CELL))
print(G, row.names = FALSE, digits = 4)

## ---- 3. calibration against the dinamr realized record ---------------------
cat("\n=== CALIBRATION (reported, NOT applied) ===\n")
cat("dinamr realized-over-projected: Block A 1.27 BEFORE the overhead correction,\n")
cat("0.949 AFTER; Block B 0.92.  The corrected projection is the one used above, so\n")
cat("the calibration factor that applies to it is 0.949/0.92, i.e. within +-8% either way.\n")
cat(sprintf("At the worst of those, the twelve-cell total would read %.3f h (x1.27 would be %.3f h).\n",
            sum(G$h_med)*0.949, sum(G$h_med)*1.27))

## ---- 4. the go/no-go and the defer decision --------------------------------
tot <- sum(G$h_med)
cat("\n=== GATE 1 TOTAL, PART A ===\n")
cat(sprintf("twelve cells : %.3f h  (90%% band %.3f-%.3f h)\n", tot, sum(G$h_q05), sum(G$h_q95)))
cat(sprintf("ceiling %.0f h wall; hard timeout %.0f h (kickoff: Larry's window is 9-10 h)\n",
            CEILING_H, TIMEOUT_H))
fit <- G$cum_h <= CEILING_H
k <- sum(fit)
cat(sprintf("cells that fit under the ceiling in run order : %d of 12 (cumulative %.3f h)\n",
            k, if (k) G$cum_h[k] else 0))
if (k < 12L) {
  cat("DEFERRED, from the tail of the run order:\n")
  for (i in which(!fit))
    cat(sprintf("  order %2d : %-6s HR %.2f n %-5d  (%.3f h)\n",
                G$order[i], G$prev[i], G$hr[i], G$n[i], G$h_med[i]))
} else cat("none deferred.\n")
cat(sprintf("\nGATE 1: %s -- run %d cell(s), defer %d.  Replicate count is NOT reduced.\n",
            if (k > 0L) "GO" else "NO-GO", k, 12L - k))
saveRDS(list(cells = G, total_h = tot, ceiling_h = CEILING_H, timeout_h = TIMEOUT_H,
             n_fit = k, overhead_by_block = OVERHEAD_BY_BLOCK, probes = tab),
        file.path(SCRATCH, "projectionG.rds"))
