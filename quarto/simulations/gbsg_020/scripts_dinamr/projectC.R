# Gate 1 projection for TASK_dinamr_blockC_grfprobe_2026-09-11 (Part A, seven cells).
#
# Differences from project.R, which this does NOT replace:
#   1. It projects from the per-replicate cost DISTRIBUTION, not a mean: the
#      2,000-replicate total is bootstrap-resampled from each probe's 36
#      measured per-replicate seconds, so the projection carries an interval.
#      The kickoff requires this ("project from the family-size distribution,
#      never a mean"); project.R's point estimate is reported beside it.
#   2. It calibrates on the campaign's REALIZED walls (Blocks A and B, read from
#      the committed bundles' mtimes), not on the probe surface alone.
#   3. It covers the seven cells this task runs: the deferred Block B cell and
#      the six Block C cells.
#
# Ceiling 9 h wall for Part A; hard timeout 12 h (kickoff, "Gate 1 - compute
# go/no-go"; the 9 h ceiling itself is Amendment 1 of the predecessor).
SCRATCH <- Sys.getenv("DINAMR_SCRATCH", unset = ".")
QMD_DIR <- Sys.getenv("DINAMR_QMD_DIR", unset = "..")
RES     <- file.path(QMD_DIR, "results")
CEILING_H <- 9; TIMEOUT_H <- 12
W <- 12L; OVERHEAD <- 30    # workers; per-render overhead (s), as project.R
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

## ---- 2. realized walls, Blocks A and B, from the committed bundles ---------
# Cell wall = batch(1) + batch(1001) + combine.  Reconstructed from mtimes; the
# very first cell's batch-1 duration has no predecessor timestamp and is taken
# equal to its own batch-2 duration (the two batches are the same size).
fs <- list.files(RES, pattern = "_dinamr_(res|combined)_.*[.]rds$", full.names = TRUE)
mt <- file.info(fs)$mtime
o  <- order(mt); fs <- basename(fs)[o]; mt <- mt[o]
dl <- c(NA_real_, as.numeric(diff(mt), units = "secs"))
dl[1] <- dl[2]                                   # first batch-1, imputed
cellof <- sub("^dina_effMaxSG_fb_mr_field_m1_(h[0-9]+)_knoise0_(n[0-9]+)(_z1q60)?_nb20_dinamr_.*$",
              "\\1_\\2\\3", fs)
real <- aggregate(list(wall_s = dl), list(cell = cellof), sum)
real$block <- ifelse(grepl("z1q60", real$cell), "B", "A")
real$n  <- as.integer(sub(".*_n([0-9]+).*", "\\1", real$cell))
real$hr <- as.numeric(sub("^h([0-9]+)_.*", "\\1", real$cell)) / 100
# project.R's Gate 1 point projections, in seconds (its printed cell_wall_s).
g1 <- c("h150_n500"=1484.3,"h150_n1000"=1345.9,"h150_n1500"=991.6,
        "h175_n500"=1484.3,"h175_n1000"=1345.9,"h175_n1500"=991.6,
        "h150_n500_z1q60"=4204.9,"h150_n1000_z1q60"=5325.4,"h150_n1500_z1q60"=5822.5,
        "h175_n500_z1q60"=4204.9,"h175_n1000_z1q60"=5325.4,"h175_n1500_z1q60"=5822.5)
real$gate1_s <- unname(g1[real$cell]); real$ratio <- real$wall_s / real$gate1_s
real <- real[order(real$block, real$hr, real$n), ]
cat("\n=== REALIZED WALLS vs GATE 1 (Blocks A and B, from the committed bundles) ===\n")
print(transform(real, wall_h = round(wall_s/3600,3), gate1_h = round(gate1_s/3600,3),
                ratio = round(ratio,3))[, c("block","hr","n","wall_h","gate1_h","ratio")],
      row.names = FALSE)
# The HR 1.75 cells were costed at the HR 1.50 corner, so their ratio absorbs
# that assumption.  Block C's corners are MEASURED at HR 1.00, so the honest
# calibration band is the measured-corner (HR 1.50) cells only.
meas <- real[real$hr == 1.50, ]
K_meas <- sum(meas$wall_s) / sum(meas$gate1_s)
K_all  <- sum(real$wall_s) / sum(real$gate1_s)
K_A    <- with(real[real$block=="A",], sum(wall_s)/sum(gate1_s))
K_B    <- with(real[real$block=="B",], sum(wall_s)/sum(gate1_s))
cat(sprintf("\nK on MEASURED corners (HR 1.50 cells only) : %.3f  (per-cell %.3f-%.3f)\n",
            K_meas, min(meas$ratio), max(meas$ratio)))
cat(sprintf("K, Block A (6 cells) : %.3f     K, Block B (5 cells) : %.3f     K, all 11 : %.3f\n",
            K_A, K_B, K_all))
cat("The kickoff's prior ratios (A 1.27, B 0.92) are realized-over-CHECKPOINT for B;\n")
cat("realized-over-GATE-1 is what is tabulated here.\n")
K_USE <- K_A   # the conservative choice: the largest block-level calibration

## ---- 3. Block C: bootstrap the 2,000-replicate wall from the 36 draws ------
boot_cell <- function(s) {
  tot <- replicate(B, sum(sample(s, 2000L, replace = TRUE)))
  (tot / W + 3 * OVERHEAD) / 3600
}
CC <- do.call(rbind, lapply(P, function(p) {
  h <- boot_cell(p$s)
  data.frame(prev = p$prev, n = p$n, basis = "measured",
             h_q05 = unname(quantile(h,.05)), h_med = median(h),
             h_q95 = unname(quantile(h,.95)),
             h_point = (2000*mean(p$s)/W + 3*OVERHEAD)/3600)
}))
# n = 1000 is not probed at HR 1.00: interpolate the per-replicate cost pool by
# mixing the n500 and n1500 draws 50/50, which interpolates the DISTRIBUTION
# rather than only its mean.
for (pv in unique(CC$prev)) {
  ps <- Filter(function(p) p$prev == pv, P)
  s5 <- ps[[which(sapply(ps, function(p) p$n) == 500L)]]$s
  s15 <- ps[[which(sapply(ps, function(p) p$n) == 1500L)]]$s
  mix <- c(s5, s15)
  h <- boot_cell(mix)
  CC <- rbind(CC, data.frame(prev = pv, n = 1000L, basis = "n500/n1500 pooled draws",
              h_q05 = unname(quantile(h,.05)), h_med = median(h),
              h_q95 = unname(quantile(h,.95)),
              h_point = (2000*mean(mix)/W + 3*OVERHEAD)/3600))
}
CC <- CC[order(CC$prev, CC$n), ]
CC$h_cal <- CC$h_med * K_USE
cat(sprintf("\n=== BLOCK C, 2000 replicates, bootstrap over the per-replicate cost distribution (B = %d) ===\n", B))
print(CC, row.names = FALSE, digits = 4)
cat(sprintf("\nBlock C uncalibrated median total : %.3f h  (90%% band %.3f-%.3f h)\n",
            sum(CC$h_med), sum(CC$h_q05), sum(CC$h_q95)))
cat(sprintf("Block C calibrated  (x %.3f)      : %.3f h\n", K_USE, sum(CC$h_cal)))
cat(sprintf("Reference to beat/correct (original Gate 1, project.R): 2.17 h\n"))

## ---- 4. the deferred Block B cell, anchored on realized Block B walls ------
# B_h175_n1500 was never run.  Two anchors, both from realized walls:
#   (a) the n-profile within HR 1.75:  wall(h175,n1500) = wall(h175,n1000) x
#       [wall(h150,n1500) / wall(h150,n1000)]
#   (b) the HR-profile within n=1500:  wall(h175,n1500) = wall(h150,n1500) x
#       [wall(h175,n1000) / wall(h150,n1000)]
g <- function(hr, n) real$wall_s[real$block=="B" & real$hr==hr & real$n==n]
a <- g(1.75,1000) * (g(1.50,1500)/g(1.50,1000))
b <- g(1.50,1500) * (g(1.75,1000)/g(1.50,1000))
cat("\n=== DEFERRED BLOCK B CELL (HR 1.75, n 1500, 31%) ===\n")
cat(sprintf("anchor (a), n-profile within HR 1.75 : %.0f s = %.3f h\n", a, a/3600))
cat(sprintf("anchor (b), HR-profile within n 1500 : %.0f s = %.3f h\n", b, b/3600))
Bdef_h <- max(a, b)/3600
cat(sprintf("taken as the LARGER of the two       : %.3f h\n", Bdef_h))
cat(sprintf("references: checkpoint re-projection 2.225 h; original Gate 1 %.3f h\n",
            g1[["h175_n1500_z1q60"]]/3600))

## ---- 5. the go/no-go -------------------------------------------------------
tot <- Bdef_h + sum(CC$h_cal)
cat("\n=== GATE 1 TOTAL, PART A (seven cells) ===\n")
cat(sprintf("deferred Block B cell : %.3f h\n", Bdef_h))
cat(sprintf("Block C, six cells    : %.3f h  (calibrated)\n", sum(CC$h_cal)))
cat(sprintf("PART A TOTAL          : %.3f h against a %.0f h ceiling (hard timeout %.0f h)\n",
            tot, CEILING_H, TIMEOUT_H))
cat(sprintf("headroom              : %.3f h (%.0f%% of the ceiling)\n",
            CEILING_H - tot, 100*(CEILING_H - tot)/CEILING_H))
cat(sprintf("room left under the %.0f h timeout for Part B's 1.5 h cap: %.3f h\n",
            TIMEOUT_H, TIMEOUT_H - tot))
cat(sprintf("\nGATE 1: %s -- all seven cells run, none deferred.\n",
            if (tot <= CEILING_H) "GO" else "NO-GO"))
saveRDS(list(blockC = CC, deferredB_h = Bdef_h, total_h = tot, K = K_USE,
             realized = real, ceiling_h = CEILING_H),
        file.path(SCRATCH, "projectionC.rds"))
