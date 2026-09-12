# Portability header added when this script was committed.  The campaign ran with
# SCRATCH = the session scratchpad and QMD_DIR = quarto/simulations/gbsg_020.
# Both are overridable; the defaults work from a clone.
SCRATCH <- Sys.getenv("DINAMR_SCRATCH", unset = ".")
QMD_DIR <- Sys.getenv("DINAMR_QMD_DIR", unset = "..")

# Gate 1 projection from the probe corners.
# Cost model: per-replicate WORKER-seconds (fit_mr_secs) is what parallelises;
# the document overhead (DGM build, tables, render) is a fixed per-render cost.
R <- file.path(QMD_DIR, "results", "")   # trailing "" keeps the %s concatenation below valid
stem <- function(hr, n, z1q) sprintf("%sdina_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d%s_nb20_dinamrprobe_res_1_36.rds",
                                     R, round(100*hr), n, if (z1q) "_z1q60" else "")
grid <- expand.grid(hr = c(1.50, 1.00), n = c(500L, 1000L, 1500L), z1q = c(FALSE, TRUE),
                    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
out <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  f <- stem(grid$hr[i], grid$n[i], grid$z1q[i])
  if (!file.exists(f)) return(NULL)
  b <- readRDS(f); r <- b$results
  s <- r$fit_mr_secs[is.finite(r$fit_mr_secs)]
  K <- r$n_family[is.finite(r$n_family)]
  data.frame(block = if (grid$z1q[i]) "B (31%)" else "A (12.4%)", hr = grid$hr[i], n = grid$n[i],
             reps = nrow(r), detection = mean(r$detected %in% 1L),
             K_med = stats::median(K), K_q90 = unname(stats::quantile(K, .9)), K_max = max(K),
             s_med = stats::median(s), s_mean = mean(s), s_q90 = unname(stats::quantile(s, .9)),
             s_max = max(s), s_total = sum(s)) }))
out <- out[order(out$block, out$hr, out$n), ]
cat("=== PROBE CORNERS (36 replicates each, 12 workers) ===\n")
print(out, row.names = FALSE, digits = 4)
cat("\n")
# Per-cell projection: 2000 replicates.  The parallel wall is driven by the SUM of
# per-replicate seconds divided by the workers (the skew averages out over 2000/12 = 167
# rounds per worker), plus a measured per-render overhead, times two batches + a combine.
W <- 12
OVERHEAD <- 30    # s per render (DGM build + tables + quarto), measured below
out$worker_s_per_rep <- out$s_mean
out$cell_wall_s <- 2000 * out$s_mean / W + 3 * OVERHEAD   # 2 batches + 1 combine render
# HR 1.75 is not probed: cost it at the HR 1.50 corner (both detect near-always; the
# family is the admission-set candidate table, which the effect target barely moves).
h175 <- out[out$hr == 1.50, ]; h175$hr <- 1.75; h175$assumed <- "costed at the HR 1.50 corner"
out$assumed <- "measured"
h100_mid <- do.call(rbind, lapply(unique(out$block), function(bk) {
  d <- out[out$block == bk & out$hr == 1.00, ]
  if (!nrow(d) || 1000L %in% d$n) return(NULL)
  x <- d[order(d$n), ]
  y <- approx(x$n, x$s_mean, xout = 1000)$y
  z <- x[1, ]; z$n <- 1000L; z$s_mean <- y; z$s_med <- NA; z$s_q90 <- NA; z$s_max <- NA
  z$K_med <- NA; z$K_q90 <- NA; z$K_max <- NA; z$detection <- NA; z$reps <- NA
  z$cell_wall_s <- 2000 * y / W + 3 * OVERHEAD; z$assumed <- "interpolated n500<->n1500"
  z }))
ALL <- rbind(out, h175, h100_mid)
ALL$wall_h <- ALL$cell_wall_s / 3600
ALL$campaign_block <- ifelse(ALL$hr == 1.00, "C (null)", ifelse(grepl("^A", ALL$block), "A", "B"))
ALL <- ALL[order(ALL$campaign_block, ALL$block, ALL$hr, ALL$n), ]
cat("=== PER-CELL PROJECTION, 2000 replicates, 12 workers ===\n")
print(ALL[, c("campaign_block","block","hr","n","s_mean","cell_wall_s","wall_h","assumed")],
      row.names = FALSE, digits = 4)
cat("\n=== BLOCK TOTALS ===\n")
tot <- aggregate(wall_h ~ campaign_block, ALL, sum)
print(tot, row.names = FALSE, digits = 4)
cat(sprintf("\nGRAND TOTAL (all 18 cells): %.2f h\n", sum(ALL$wall_h)))
cat(sprintf("A + B only (12 cells)     : %.2f h\n", sum(ALL$wall_h[ALL$campaign_block %in% c("A","B")])))
cat(sprintf("A only (6 cells)          : %.2f h\n", sum(ALL$wall_h[ALL$campaign_block == "A"])))
cat(sprintf("B only (6 cells)          : %.2f h\n", sum(ALL$wall_h[ALL$campaign_block == "B"])))
cat(sprintf("C only (6 cells)          : %.2f h\n", sum(ALL$wall_h[ALL$campaign_block == "C (null)"])))
saveRDS(ALL, file.path(SCRATCH, "projection.rds"))
