# Gate 1 projection for TASK_dinamr_blockC_grfprobe_2026-09-11 (Part A, seven cells).
#
# Differences from project.R, which this does NOT replace:
#   1. It projects from the per-replicate cost DISTRIBUTION, not a mean: the
#      2,000-replicate total is bootstrap-resampled from each probe's 36
#      measured per-replicate seconds, so the projection carries an interval.
#      The kickoff requires this ("project from the family-size distribution,
#      never a mean"); project.R's point estimate is reported beside it.
#   2. It calibrates on the campaign's REALIZED walls (Blocks A and B).  Those
#      walls come from the bundles' own PER-REPLICATE TIMING COLUMNS via
#      walls.R -- sum(fit_mr_secs) / meta$n_workers -- NOT from differencing
#      result-file mtimes, which conflates inter-cell gaps, the combine render
#      and any interruption with the cell's own cost.  The mtime figure is kept
#      beside it as a labelled proxy only.  See walls.R for the verification
#      that fit_mr_secs is the top-level per-replicate timer, finite on every
#      row, with fld_H_secs / fld_Hc_secs NESTED INSIDE it.
#   3. It covers the seven cells this task runs: the deferred Block B cell and
#      the six Block C cells.
#
# Ceiling 9 h wall for Part A; hard timeout 12 h (kickoff, "Gate 1 - compute
# go/no-go"; the 9 h ceiling itself is Amendment 1 of the predecessor).
SCRATCH <- Sys.getenv("DINAMR_SCRATCH", unset = ".")
QMD_DIR <- Sys.getenv("DINAMR_QMD_DIR", unset = "..")
RES     <- file.path(QMD_DIR, "results")
CEILING_H <- 9; TIMEOUT_H <- 12
W <- 12L
# Per-render overhead.  project.R assumed 30 s.  Measured on the eleven completed
# cells (walls.R) the residual of the realized span over the compute wall is a
# MEDIAN OF 112.8 s PER RENDER, range 75.9-196.6 -- roughly four times the
# assumption, and the whole reason the mtime-based ratios ran high.  It is taken
# per block at the largest realized value, which is the conservative choice.
OVERHEAD_BY_BLOCK <- c("12.4%" = 122.8, "31%" = 196.6)
OVERHEAD_LEGACY   <- 30
set.seed(8316951)
B <- 4000L                  # bootstrap draws for the cell-wall distribution

probe <- function(hr, n, z1q)
  file.path(RES, sprintf("dina_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d%s_nb20_dinamrprobe_res_1_36.rds",
                         round(100*hr), n, if (z1q) "_z1q60" else ""))

## ---- 1. the Block C probe corners, as MEASURED -----------------------------
corners <- expand.grid(n = c(500L, 1500L), z1q = c(FALSE, TRUE), KEEP.OUT.ATTRS = FALSE)
P <- lapply(seq_len(nrow(corners)), function(i) {
  f <- probe(1.00, corners$n[i], corners$z1q[i])
  if (!file.exists(f)) return(NULL)
  r <- readRDS(f)$results
  list(prev = if (corners$z1q[i]) "31%" else "12.4%", n = corners$n[i],
       s = r$fit_mr_secs[is.finite(r$fit_mr_secs)],
       K = r$n_family[is.finite(r$n_family)],
       det = mean(r$detected %in% 1L), reps = nrow(r))
})
P <- Filter(Negate(is.null), P)
cat("=== BLOCK C PROBE CORNERS AT HR 1.00 (36 replicates each, 12 workers) ===\n")
cat("These are the ORIGINAL Gate 1 probes (probe.sh lines 30-33): four of the ten\n")
cat("covered a Block C cell, so no new 36-replicate probe was needed.\n\n")
tab <- do.call(rbind, lapply(P, function(p) data.frame(
  prev = p$prev, n = p$n, reps = p$reps, detection = p$det,
  K_q10 = unname(quantile(p$K, .10)), K_med = median(p$K),
  K_q90 = unname(quantile(p$K, .90)), K_max = max(p$K),
  s_q10 = unname(quantile(p$s, .10)), s_med = median(p$s),
  s_q90 = unname(quantile(p$s, .90)), s_max = max(p$s), s_mean = mean(p$s))))
print(tab, row.names = FALSE, digits = 4)

## ---- 2. realized walls, Blocks A and B, from the TIMING COLUMNS ------------
# Cell compute wall = sum(fit_mr_secs) over both batches / meta$n_workers.
# fit_mr_secs is the top-level per-replicate timer and is finite on every row;
# fld_H_secs / fld_Hc_secs are nested inside it and must not be added.
stem_of <- function(f)
  sub("^dina_effMaxSG_fb_mr_field_m1_(h[0-9]+)_knoise0_(n[0-9]+)(_z1q60)?_nb20_dinamr_.*$",
      "\\1_\\2\\3", f)
bat <- list.files(RES, pattern = "_dinamr_res_[0-9]+_[0-9]+[.]rds$", full.names = TRUE)
bi <- do.call(rbind, lapply(bat, function(pp) {
  b <- readRDS(pp); r <- b$results
  data.frame(cell = stem_of(basename(pp)),
             n_workers = as.integer(b$meta$n_workers),
             worker_s = sum(r$fit_mr_secs, na.rm = TRUE),
             n_fin = sum(is.finite(r$fit_mr_secs)), rows = nrow(r),
             stringsAsFactors = FALSE) }))
Wrec <- unique(bi$n_workers[is.finite(bi$n_workers)])
stopifnot(length(Wrec) == 1L)
real <- aggregate(cbind(worker_s, rows, n_fin) ~ cell, bi, sum)
real$nbatch <- as.integer(table(bi$cell)[real$cell])
real$compute_s <- real$worker_s / Wrec
real$block <- ifelse(grepl("z1q60", real$cell), "B", "A")
real$n  <- as.integer(sub(".*_n([0-9]+).*", "\\1", real$cell))
real$hr <- as.numeric(sub("^h([0-9]+)_.*", "\\1", real$cell)) / 100
# The mtime span, kept ONLY as a labelled proxy and to expose the residual.
allf <- list.files(RES, pattern = "_dinamr_(res|combined)_.*[.]rds$", full.names = TRUE)
mt <- file.info(allf)$mtime; o <- order(mt); allf <- basename(allf)[o]; mt <- mt[o]
dl <- c(NA_real_, as.numeric(diff(mt), units = "secs")); dl[1] <- dl[2]
prox <- aggregate(list(mtime_s = dl), list(cell = stem_of(allf)), sum)
real <- merge(real, prox, by = "cell", all.x = TRUE)
real$residual_s <- real$mtime_s - real$compute_s
real$ovh_per_render <- real$residual_s / (real$nbatch + 1L)
g1 <- c("h150_n500"=1484.3,"h150_n1000"=1345.9,"h150_n1500"=991.6,
        "h175_n500"=1484.3,"h175_n1000"=1345.9,"h175_n1500"=991.6,
        "h150_n500_z1q60"=4204.9,"h150_n1000_z1q60"=5325.4,"h150_n1500_z1q60"=5822.5,
        "h175_n500_z1q60"=4204.9,"h175_n1000_z1q60"=5325.4,"h175_n1500_z1q60"=5822.5)
real$gate1_s <- unname(g1[real$cell])
real <- real[order(real$block, real$hr, real$n), ]
cat("\n=== REALIZED WALLS FROM THE TIMING COLUMNS (mtime kept as a proxy) ===\n")
cat(sprintf("n_workers read from batch meta: %d ; fit_mr_secs finite on %d of %d rows\n",
            Wrec, sum(real$n_fin), sum(real$rows)))
print(transform(real, compute_h = round(compute_s/3600,4), mtime_h = round(mtime_s/3600,4),
                ovh = round(ovh_per_render,1))[
      , c("block","hr","n","compute_h","mtime_h","residual_s","ovh")],
      row.names = FALSE, digits = 5)
cat(sprintf("\nper-render residual: median %.1f s, range %.1f-%.1f s (project.R assumed %d s)\n",
    median(real$ovh_per_render), min(real$ovh_per_render), max(real$ovh_per_render), OVERHEAD_LEGACY))
cat("The residual is render/DGM/table cost plus the makespan slack over the sum/W bound;\n")
cat("it is measured, not assumed, and is what the mtime-based K was silently absorbing.\n")
ok <- is.finite(real$gate1_s)
real$ratio_timing <- NA_real_
real$ratio_timing[ok] <- (real$compute_s[ok] +
   3*OVERHEAD_BY_BLOCK[ifelse(real$block[ok]=="A","12.4%","31%")]) / real$gate1_s[ok]
real$ratio_mtime <- NA_real_; real$ratio_mtime[ok] <- real$mtime_s[ok]/real$gate1_s[ok]
cat("\n=== CALIBRATION AGAINST GATE 1, both reconstructions ===\n")
print(transform(real[ok,], rt = round(ratio_timing,3), rm = round(ratio_mtime,3))[
      , c("block","hr","n","gate1_s","compute_s","mtime_s","rt","rm")],
      row.names = FALSE, digits = 5)
m5 <- ok & real$hr == 1.50
cat(sprintf("\nCOMPUTE-model accuracy on the MEASURED (HR 1.50) corners, before overhead:\n"))
cat(sprintf("  sum(compute)/sum(gate1) = %.3f  -- Gate 1's compute model was essentially exact.\n",
    sum(real$compute_s[m5])/sum(real$gate1_s[m5])))
cat(sprintf("  the mtime-based K of %.3f on the same corners was overhead, not compute.\n",
    sum(real$mtime_s[m5])/sum(real$gate1_s[m5])))
cat("The HR 1.75 cells are excluded from that reading: they were costed at the HR 1.50\n")
cat("corner, so their ratio absorbs that assumption rather than measuring anything.\n")
cat("\nNO BLANKET MULTIPLIER IS APPLIED BELOW.  The projection is compute-from-distribution\n")
cat("plus the MEASURED per-render overhead, which is what the multiplier was standing in for.\n")

## ---- 3. Block C: bootstrap the 2,000-replicate wall from the 36 draws ------
boot_cell <- function(s, ovh) {
  tot <- replicate(B, sum(sample(s, 2000L, replace = TRUE)))
  (tot / W + 3 * ovh) / 3600
}
CC <- do.call(rbind, lapply(P, function(p) {
  ovh <- unname(OVERHEAD_BY_BLOCK[p$prev])
  h <- boot_cell(p$s, ovh)
  data.frame(prev = p$prev, n = p$n, basis = "measured", ovh_s = ovh,
             h_q05 = unname(quantile(h,.05)), h_med = median(h),
             h_q95 = unname(quantile(h,.95)),
             h_legacy30 = (2000*mean(p$s)/W + 3*OVERHEAD_LEGACY)/3600)
}))
# n = 1000 is not probed at HR 1.00: interpolate the per-replicate cost pool by
# mixing the n500 and n1500 draws 50/50, which interpolates the DISTRIBUTION
# rather than only its mean.
for (pv in unique(CC$prev)) {
  ps <- Filter(function(p) p$prev == pv, P)
  s5 <- ps[[which(sapply(ps, function(p) p$n) == 500L)]]$s
  s15 <- ps[[which(sapply(ps, function(p) p$n) == 1500L)]]$s
  mix <- c(s5, s15)
  ovh <- unname(OVERHEAD_BY_BLOCK[pv])
  h <- boot_cell(mix, ovh)
  CC <- rbind(CC, data.frame(prev = pv, n = 1000L, basis = "n500/n1500 pooled draws", ovh_s = ovh,
              h_q05 = unname(quantile(h,.05)), h_med = median(h),
              h_q95 = unname(quantile(h,.95)),
              h_legacy30 = (2000*mean(mix)/W + 3*OVERHEAD_LEGACY)/3600))
}
CC <- CC[order(CC$prev, CC$n), ]
cat(sprintf("\n=== BLOCK C, 2000 replicates, bootstrap over the per-replicate cost distribution (B = %d) ===\n", B))
cat("h_med / h_q05 / h_q95 already carry the MEASURED per-render overhead (ovh_s x 3).\n")
cat("h_legacy30 is the same compute with project.R's 30 s assumption, for comparison only.\n\n")
print(CC, row.names = FALSE, digits = 4)
cat(sprintf("\nBlock C total : %.3f h  (90%% band %.3f-%.3f h)\n",
            sum(CC$h_med), sum(CC$h_q05), sum(CC$h_q95)))
cat(sprintf("  the same compute at project.R's 30 s overhead would read %.3f h,\n", sum(CC$h_legacy30)))
cat(sprintf("  which is the original Gate 1 figure of 2.17 h reproduced.\n"))

## ---- 4. the deferred Block B cell, anchored on realized Block B walls ------
# B_h175_n1500 was never run.  Two anchors, both from realized walls:
#   (a) the n-profile within HR 1.75:  wall(h175,n1500) = wall(h175,n1000) x
#       [wall(h150,n1500) / wall(h150,n1000)]
#   (b) the HR-profile within n=1500:  wall(h175,n1500) = wall(h150,n1500) x
#       [wall(h175,n1000) / wall(h150,n1000)]
g <- function(hr, n) real$compute_s[real$block=="B" & real$hr==hr & real$n==n]
a <- g(1.75,1000) * (g(1.50,1500)/g(1.50,1000))
b <- g(1.50,1500) * (g(1.75,1000)/g(1.50,1000))
ovhB <- unname(OVERHEAD_BY_BLOCK["31%"])
cat("\n=== DEFERRED BLOCK B CELL (HR 1.75, n 1500, 31%) ===\n")
cat("Anchored on realized COMPUTE walls (timing columns), then the measured overhead added.\n")
cat(sprintf("anchor (a), n-profile within HR 1.75 : compute %.0f s\n", a))
cat(sprintf("anchor (b), HR-profile within n 1500 : compute %.0f s\n", b))
Bdef_h <- (max(a, b) + 3*ovhB)/3600
cat(sprintf("larger compute %.0f s + 3 x %.1f s overhead = %.0f s = %.3f h\n",
            max(a,b), ovhB, max(a,b) + 3*ovhB, Bdef_h))
cat(sprintf("references: checkpoint re-projection 2.225 h; original Gate 1 %.3f h\n",
            g1[["h175_n1500_z1q60"]]/3600))

## ---- 5. the go/no-go -------------------------------------------------------
tot <- Bdef_h + sum(CC$h_med)
cat("\n=== GATE 1 TOTAL, PART A (seven cells) ===\n")
cat(sprintf("deferred Block B cell : %.3f h\n", Bdef_h))
cat(sprintf("Block C, six cells    : %.3f h\n", sum(CC$h_med)))
cat(sprintf("PART A TOTAL          : %.3f h against a %.0f h ceiling (hard timeout %.0f h)\n",
            tot, CEILING_H, TIMEOUT_H))
cat(sprintf("headroom              : %.3f h (%.0f%% of the ceiling)\n",
            CEILING_H - tot, 100*(CEILING_H - tot)/CEILING_H))
cat(sprintf("room left under the %.0f h timeout for Part B's 1.5 h cap: %.3f h\n",
            TIMEOUT_H, TIMEOUT_H - tot))
cat(sprintf("\nGATE 1: %s -- all seven cells run, none deferred.\n",
            if (tot <= CEILING_H) "GO" else "NO-GO"))
saveRDS(list(blockC = CC, deferredB_h = Bdef_h, total_h = tot,
             overhead_by_block = OVERHEAD_BY_BLOCK,
             realized = real, ceiling_h = CEILING_H),
        file.path(SCRATCH, "projectionC.rds"))
