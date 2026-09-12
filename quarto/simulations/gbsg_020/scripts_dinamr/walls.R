# Realized cell walls from the bundles' own PER-REPLICATE TIMING COLUMNS.
# TASK_dinamr_blockC_grfprobe_2026-09-11, Gate 1 item (3).
#
# WHY THIS REPLACES mtime DIFFERENCING.  The first Gate 1 pass reconstructed
# realized walls by differencing the results files' mtimes.  That conflates
# whatever sat between two writes -- inter-cell gaps, the combine render,
# parallel writes, any interruption -- with the cell's own cost.  The bundles do
# carry per-replicate timing, so the walls are reconstructed from that instead
# and the mtime figure is kept beside it only as corroboration, labelled a proxy.
#
# WHAT THE COLUMNS ARE.  Verified on the committed bundles:
#   * fit_mr_secs is the per-replicate TOP-LEVEL worker time and is finite on
#     EVERY replicate (1000/1000 on the batch checked), detected or not.
#   * fld_H_secs and fld_Hc_secs are NESTED INSIDE it -- fld_H_secs <=
#     fit_mr_secs on 878 of 878 rows where both are finite (median ratio 0.829,
#     max 0.940), and fld_H_secs + fld_Hc_secs <= fit_mr_secs on all of them.
#     So they must NOT be added to fit_mr_secs; that would double-count.
#   * fb_secs and fld_H_uniform_secs are identically zero on this campaign
#     (FS_S7_FB=none, field_uniform off).
# There is no per-replicate WALL column, so the compute wall is the summed
# worker-seconds divided by the recorded n_workers.  meta$n_workers is recorded
# in every BATCH meta (12 here); the combined meta does not carry it.
#
#   compute_wall = sum(fit_mr_secs) / n_workers
#   cell_wall    = compute_wall + per-render overhead x 3 (two batches + combine)
#
# The per-render overhead is not itself recorded; it is ESTIMATED here as the
# residual of the mtime span over the compute wall, and reported, not assumed.
SCRATCH <- Sys.getenv("DINAMR_SCRATCH", unset = ".")
QMD_DIR <- Sys.getenv("DINAMR_QMD_DIR", unset = "..")
RES     <- file.path(QMD_DIR, "results")

stem_of <- function(f)
  sub("^dina_effMaxSG_fb_mr_field_m1_(h[0-9]+)_knoise0_(n[0-9]+)(_z1q60)?_nb20_dinamr_.*$",
      "\\1_\\2\\3", f)

bat <- list.files(RES, pattern = "_dinamr_res_[0-9]+_[0-9]+[.]rds$", full.names = TRUE)
if (!length(bat)) stop("no dinamr batch bundles on disk")
info <- do.call(rbind, lapply(bat, function(p) {
  b <- readRDS(p); m <- b$meta; r <- b$results
  data.frame(file = basename(p), cell = stem_of(basename(p)),
             start = m$sim_id_start, nsims = m$n_sims,
             n_workers = as.integer(m$n_workers %||% NA),
             built_at = as.POSIXct(m$built_at),
             mtime = file.info(p)$mtime,
             worker_s = sum(r$fit_mr_secs, na.rm = TRUE),
             n_fin = sum(is.finite(r$fit_mr_secs)), rows = nrow(r),
             stringsAsFactors = FALSE) }))
`%||%` <- function(a,b) if (is.null(a) || length(a)==0 || all(is.na(a))) b else a

cat("=== BATCH TIMING, read from the bundles ===\n")
print(info[order(info$built_at), c("cell","start","rows","n_fin","n_workers","worker_s","built_at")],
      row.names = FALSE, digits = 6)
if (any(info$n_fin != info$rows))
  cat("\nNOTE: fit_mr_secs is not finite on every row of every batch; the sums below\n",
      "     therefore understate those cells.  Affected: ",
      paste(unique(info$cell[info$n_fin != info$rows]), collapse = ", "), "\n", sep = "")

W <- unique(info$n_workers[is.finite(info$n_workers)])
cat(sprintf("\nn_workers recorded in batch meta: %s\n", paste(W, collapse = "/")))
if (length(W) != 1L) stop("batches disagree about n_workers; not reconstructing walls")

## ---- per-cell: the timing-column wall -------------------------------------
CELL <- aggregate(cbind(worker_s, rows) ~ cell, info, sum)
CELL$nbatch <- as.integer(table(info$cell)[CELL$cell])
CELL$compute_wall_s <- CELL$worker_s / W
CELL$block <- ifelse(grepl("z1q60", CELL$cell), "B", "A")
CELL$n  <- as.integer(sub(".*_n([0-9]+).*", "\\1", CELL$cell))
CELL$hr <- as.numeric(sub("^h([0-9]+)_.*", "\\1", CELL$cell)) / 100

## ---- the mtime proxy, kept for corroboration only -------------------------
allf <- list.files(RES, pattern = "_dinamr_(res|combined)_.*[.]rds$", full.names = TRUE)
mt <- file.info(allf)$mtime; o <- order(mt); allf <- basename(allf)[o]; mt <- mt[o]
dl <- c(NA_real_, as.numeric(diff(mt), units = "secs")); dl[1] <- dl[2]  # first batch-1, imputed
PROXY <- aggregate(list(mtime_wall_s = dl), list(cell = stem_of(allf)), sum)
CELL <- merge(CELL, PROXY, by = "cell", all.x = TRUE)
CELL$overhead_s <- CELL$mtime_wall_s - CELL$compute_wall_s
CELL <- CELL[order(CELL$block, CELL$hr, CELL$n), ]

cat("\n=== REALIZED CELL WALLS ===\n")
cat("compute_wall_s = sum(fit_mr_secs) / n_workers -- the bundles' own timing.\n")
cat("mtime_wall_s   = PROXY, the span between consecutive result-file writes.\n")
cat("overhead_s     = the residual: render/DGM/table cost plus any gap the proxy swept up.\n\n")
print(transform(CELL,
        compute_h = round(compute_wall_s/3600, 4),
        mtime_h   = round(mtime_wall_s/3600, 4),
        ovh_per_render_s = round(overhead_s/pmax(CELL$nbatch + 1L, 1L), 1))[
      , c("block","hr","n","nbatch","rows","worker_s","compute_h","mtime_h","overhead_s","ovh_per_render_s")],
      row.names = FALSE, digits = 5)
cat(sprintf("\ntotal compute wall over %d cells: %.3f h ; mtime proxy total: %.3f h ; residual %.3f h\n",
    nrow(CELL), sum(CELL$compute_wall_s)/3600, sum(CELL$mtime_wall_s)/3600,
    sum(CELL$overhead_s)/3600))
cat(sprintf("per-render overhead implied: median %.1f s, range %.1f-%.1f s (project.R assumed 30 s)\n",
    median(CELL$overhead_s/(CELL$nbatch+1)), min(CELL$overhead_s/(CELL$nbatch+1)),
    max(CELL$overhead_s/(CELL$nbatch+1))))

## ---- recalibrate Gate 1 on the timing-column walls -------------------------
g1 <- c("h150_n500"=1484.3,"h150_n1000"=1345.9,"h150_n1500"=991.6,
        "h175_n500"=1484.3,"h175_n1000"=1345.9,"h175_n1500"=991.6,
        "h150_n500_z1q60"=4204.9,"h150_n1000_z1q60"=5325.4,"h150_n1500_z1q60"=5822.5,
        "h175_n500_z1q60"=4204.9,"h175_n1000_z1q60"=5325.4,"h175_n1500_z1q60"=5822.5)
CELL$gate1_s <- unname(g1[CELL$cell])
CELL$ratio_timing <- CELL$mtime_wall_s * NA
ok <- is.finite(CELL$gate1_s)
CELL$ratio_timing[ok] <- (CELL$compute_wall_s[ok] + 3*30) / CELL$gate1_s[ok]
CELL$ratio_mtime[ok]  <- CELL$mtime_wall_s[ok] / CELL$gate1_s[ok]
cat("\n=== CALIBRATION AGAINST GATE 1, both reconstructions ===\n")
print(transform(CELL[ok, ], ratio_timing = round(ratio_timing,3), ratio_mtime = round(ratio_mtime,3))[
      , c("block","hr","n","gate1_s","compute_wall_s","mtime_wall_s","ratio_timing","ratio_mtime")],
      row.names = FALSE, digits = 5)
Ktim <- sum(CELL$compute_wall_s[ok] + 3*30) / sum(CELL$gate1_s[ok])
Kmt  <- sum(CELL$mtime_wall_s[ok]) / sum(CELL$gate1_s[ok])
mA <- ok & CELL$block=="A"; mB <- ok & CELL$block=="B"; m5 <- ok & CELL$hr==1.50
cat(sprintf("\nK (timing columns): all %.3f | Block A %.3f | Block B %.3f | measured HR 1.50 corners %.3f\n",
    Ktim, sum(CELL$compute_wall_s[mA]+3*30)/sum(CELL$gate1_s[mA]),
    sum(CELL$compute_wall_s[mB]+3*30)/sum(CELL$gate1_s[mB]),
    sum(CELL$compute_wall_s[m5]+3*30)/sum(CELL$gate1_s[m5])))
cat(sprintf("K (mtime proxy)   : all %.3f | Block A %.3f | Block B %.3f | measured HR 1.50 corners %.3f\n",
    Kmt, sum(CELL$mtime_wall_s[mA])/sum(CELL$gate1_s[mA]),
    sum(CELL$mtime_wall_s[mB])/sum(CELL$gate1_s[mB]),
    sum(CELL$mtime_wall_s[m5])/sum(CELL$gate1_s[m5])))
saveRDS(CELL, file.path(SCRATCH, "walls.rds"))
