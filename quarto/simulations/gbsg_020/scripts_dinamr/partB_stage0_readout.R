# Part B initial measurement (TASK_partB_measurement_2026-09-12) -- the Stage 0
# readout.  READS COMMITTED BUNDLES ONLY: no render, no simulation, no MR-off run.
#
# Stage 0c found that the template hard-codes mr_inference = TRUE
# (sim_fs_maxeffCons_fb_mr_field_m1_template.qmd:1017) with no FS_S7_* knob, so
# Stage M did not run and there is NO measured MR-off cost.  Everything below
# is MR-ON data from the committed bundles, and the two cost figures it derives
# are BOUNDS, labelled as such:
#   * CEILING  -- fit_mr_secs itself.  MR off cannot cost more than MR on.
#   * FIELD-EXCLUDED BOUND -- fit_mr_secs - fld_H_secs - fld_Hc_secs.  Both
#     field timers are nested inside fit_mr_secs and disjoint from each other
#     (walls.R:14-16: fld_H_secs + fld_Hc_secs <= fit_mr_secs on every row), so
#     this removes the field block only.  The rest of MR (the de-biasing draws,
#     the IJ, the reselection record) stays in, so MR-off cost is BELOW this
#     bound, and 1 - bound/fit is a LOWER bound on the MR share.  NA field
#     timers (non-detected rows, where no field ran) count as 0.
# Neither is a substitute for the Stage M measurement.
#
# usage (from scripts_dinamr/): Rscript partB_stage0_readout.R
SCRATCH <- Sys.getenv("DINAMR_SCRATCH", unset = ".")
QMD_DIR <- Sys.getenv("DINAMR_QMD_DIR", unset = "..")
RES <- file.path(QMD_DIR, "results")
W <- 12L

rd <- function(stem) {
  f <- file.path(RES, paste0(stem, c("_res_1_1000.rds", "_res_1001_2000.rds")))
  stopifnot(all(file.exists(f)))
  b <- lapply(f, readRDS)
  r <- do.call(rbind, lapply(b, `[[`, "results"))
  z <- function(x) { x[!is.finite(x)] <- 0; x }
  r$fit_minus_field <- r$fit_mr_secs - z(r$fld_H_secs) - z(r$fld_Hc_secs)
  list(r = r, meta = b[[1]]$meta, bytes = sum(file.size(f)))
}
q4 <- function(x) round(unname(quantile(x, c(.25, .5, .75, .9), na.rm = TRUE)), 2)

## ---- 1. the three comparators at 31% HR 1.50 n 500, effMaxSG eps 0.20 -----
# FS: the DESIGNATED comparator (fs_extraction.R / current_status.md section
# 2.1) is e1stud.  p30sgnb20 is the same cell and criterion, shown beside it.
comp <- c(consistency   = "fs_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_z1q60_nb20_e1stud",
          "consistency (p30sgnb20)" = "fs_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_z1q60_nb20_p30sgnb20",
          dina          = "dina_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_z1q60_nb20_dinamr",
          grf           = "grf_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_z1q60_nb20_grfmr")
cat("=== 1. COMPARATORS at 31% HR 1.50 n 500 (MR ON, 2,000 replicates each) ===\n")
C1 <- do.call(rbind, lapply(names(comp), function(e) {
  x <- rd(comp[[e]]); r <- x$r; m <- x$meta
  s <- q4(r$fit_mr_secs); d <- q4(r$fit_minus_field)
  data.frame(engine = e, host = m$hostname, workers = m$n_workers, R = m$r_version,
             focus = m$sg_focus, eps = m$effect_neighborhood, reps = nrow(r),
             sim_ids = paste(range(r$sim_id), collapse = "-"),
             detected = mean(r$detected %in% 1L),
             on_q25 = s[1], on_med = s[2], on_q75 = s[3], on_p90 = s[4],
             bnd_q25 = d[1], bnd_med = d[2], bnd_q75 = d[3], bnd_p90 = d[4],
             share_lower_bound = round(1 - median(r$fit_minus_field) / median(r$fit_mr_secs), 3),
             kb_per_rep = round(x$bytes / nrow(r) / 1024, 3),
             stringsAsFactors = FALSE)
}))
print(C1, row.names = FALSE)

## ---- 2. is FS cost flat across sg_focus?  same cell, same host -------------
cat("\n=== 2. FS ACROSS sg_focus at 31% HR 1.50 n 500 (all pop-os, 100 workers, MR ON) ===\n")
foc <- c("maxeffCons (p30)"     = "fs_maxeffCons_fb_mr_field_m1_h150_knoise0_n500_z1q60_p30",
         "effMaxSG 0.10 (p30sg)" = "fs_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_z1q60_p30sg",
         "effMaxSG 0.20 (p30sgnb20)" = "fs_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_z1q60_nb20_p30sgnb20",
         "effMaxSG 0.20 (e1stud)" = "fs_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_z1q60_nb20_e1stud",
         "effMaxSG 0.30 (banddial)" = "fs_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_z1q60_nb30_banddial",
         "maxSG (banddial)"     = "fs_maxSG_fb_mr_field_m1_h150_knoise0_n500_z1q60_banddial",
         "minSG (banddial)"     = "fs_minSG_fb_mr_field_m1_h150_knoise0_n500_z1q60_banddial")
C2 <- do.call(rbind, lapply(names(foc), function(k) {
  x <- rd(foc[[k]]); r <- x$r
  data.frame(focus = k, host = x$meta$hostname, workers = x$meta$n_workers,
             on_med = round(median(r$fit_mr_secs), 2),
             bnd_med = round(median(r$fit_minus_field), 2),
             field_med = round(median(r$fit_mr_secs - r$fit_minus_field), 2),
             share_lower_bound = round(1 - median(r$fit_minus_field) / median(r$fit_mr_secs), 3),
             stringsAsFactors = FALSE)
}))
print(C2, row.names = FALSE)
cat("Not covered by any committed bundle: maxeff (consistency) and effMinSG.\n")

## ---- 3. host factor: the same FS cell on pop-os/100 and on the Mac/12 -------
cat("\n=== 3. HOST FACTOR, FS maxeffCons 12.4% n 500 (sim_id 1-1000, same seeds) ===\n")
hf <- list(c("h100", "s7c", "tier2"), c("h175", "s7c", "tier2"))
C3 <- do.call(rbind, lapply(hf, function(h) {
  one <- function(camp) {
    f <- file.path(RES, sprintf("fs_maxeffCons_fb_mr_field_m1_%s_knoise0_n500_%s_res_1_1000.rds", h[1], camp))
    b <- readRDS(f); r <- b$results
    z <- function(x) { x[!is.finite(x)] <- 0; x }
    c(host = b$meta$hostname, W = b$meta$n_workers,
      med = median(r$fit_mr_secs), bnd = median(r$fit_mr_secs - z(r$fld_H_secs) - z(r$fld_Hc_secs)))
  }
  a <- one(h[2]); m <- one(h[3])
  data.frame(cell = h[1], popos = sprintf("%s/%s", a["host"], a["W"]), mac = sprintf("%s/%s", m["host"], m["W"]),
             on_med_popos = round(as.numeric(a["med"]), 2), on_med_mac = round(as.numeric(m["med"]), 2),
             ratio_on = round(as.numeric(a["med"]) / as.numeric(m["med"]), 2),
             bnd_popos = round(as.numeric(a["bnd"]), 2), bnd_mac = round(as.numeric(m["bnd"]), 2),
             ratio_bnd = round(as.numeric(a["bnd"]) / as.numeric(m["bnd"]), 2))
}))
print(C3, row.names = FALSE)

## ---- 4. per-engine 18-cell cost, MR ON, from every committed cell -----------
# FS uses the designated grid (fs_extraction.R's fsfile()): at 12.4% maxeffCons
# eps 0.10 (tier2 / p12ext), at 31% effMaxSG eps 0.20 (e1stud / cert20).
fs_stem <- function(hr, n, prev) {
  if (prev == "12.4%") {
    camp <- if (abs(hr - 1.75) < 1e-9 || (abs(hr - 1) < 1e-9 && n == 500L)) "tier2" else "p12ext"
    sprintf("fs_maxeffCons_fb_mr_field_m1_h%03d_knoise0_n%d_%s", round(100 * hr), n, camp)
  } else {
    camp <- if (n == 500L && abs(hr - 1) > 1e-9) "e1stud" else "cert20"
    sprintf("fs_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d_z1q60_nb20_%s", round(100 * hr), n, camp)
  }
}
eng_stem <- function(e, hr, n, prev)
  sprintf("%s_effMaxSG_fb_mr_field_m1_h%03d_knoise0_n%d%s_nb20_%smr", e, round(100 * hr), n,
          if (prev == "31%") "_z1q60" else "", e)
grid <- expand.grid(prev = c("12.4%", "31%"), hr = c(1.00, 1.50, 1.75), n = c(500L, 1000L, 1500L),
                    stringsAsFactors = FALSE)
C4 <- do.call(rbind, lapply(c("consistency", "dina", "grf"), function(e) {
  do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
    g <- grid[i, ]
    st <- if (e == "consistency") fs_stem(g$hr, g$n, g$prev) else eng_stem(e, g$hr, g$n, g$prev)
    x <- rd(st); r <- x$r
    data.frame(engine = e, prev = g$prev, hr = g$hr, n = g$n, host = x$meta$hostname,
               workers = x$meta$n_workers, reps = nrow(r),
               on_med = round(median(r$fit_mr_secs), 2),
               on_worker_h = sum(r$fit_mr_secs) / 3600,
               bnd_worker_h = sum(r$fit_minus_field) / 3600,
               stringsAsFactors = FALSE)
  }))
}))
cat("\n=== 4. PER-CELL, MR ON, 2,000 replicates (worker-hours = sum of per-replicate seconds / 3600) ===\n")
print(transform(C4, on_worker_h = round(on_worker_h, 2), bnd_worker_h = round(bnd_worker_h, 2)), row.names = FALSE)

## ---- 5. what the bounds imply for 288 cell-runs -- NOT the Gate P projection
runs <- c(consistency = 6L, dina = 5L, grf = 5L)
agg <- aggregate(cbind(on_worker_h, bnd_worker_h) ~ engine, C4, sum)
agg$hosts <- tapply(paste0(C4$host, "/", C4$workers), C4$engine,
                    function(h) paste(names(table(h)), table(h), sep = " x", collapse = "; "))[agg$engine]
agg$criteria <- runs[agg$engine]
cat("\n=== 5. BOUNDS FOR THE PART B SWEEP (criteria x 18 cells), NOT A PROJECTION ===\n")
cat("One criterion's 18-cell cost is taken from the committed effMaxSG-family grid and\n")
cat("assumed flat across sg_focus (section 2 supports this for the field-excluded bound\n")
cat("on FS; it does NOT for the MR-on ceiling).  Wall = worker-hours / 12.\n")
cat("FS hours are mostly POP-OS / 100-WORKER seconds (section 3: 4-6x the Mac's per replicate).\n\n")
C5 <- do.call(rbind, lapply(c(2000L, 1000L, 500L), function(R) {
  s <- R / 2000
  data.frame(reps = R, engine = agg$engine, cell_runs = 18L * agg$criteria,
             ceiling_worker_h = round(agg$on_worker_h * agg$criteria * s, 1),
             ceiling_wall_h_12w = round(agg$on_worker_h * agg$criteria * s / W, 1),
             bound_worker_h = round(agg$bnd_worker_h * agg$criteria * s, 1),
             bound_wall_h_12w = round(agg$bnd_worker_h * agg$criteria * s / W, 1),
             hosts = agg$hosts, stringsAsFactors = FALSE)
}))
print(C5, row.names = FALSE)
tot <- aggregate(cbind(ceiling_wall_h_12w, bound_wall_h_12w) ~ reps, C5, sum)
cat("\nAll three engines, wall hours at 12 workers:\n"); print(tot[order(-tot$reps), ], row.names = FALSE)

## ---- 6. Monte Carlo resolution --------------------------------------------
wilson_hw <- function(p, n) { z <- qnorm(.975); z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / (1 + z^2 / n) }
C6 <- data.frame(reps = c(2000L, 1000L, 500L))
C6$hw_p50 <- round(wilson_hw(.50, C6$reps), 4)
C6$hw_p90 <- round(wilson_hw(.90, C6$reps), 4)
C6$thinnest_detected <- round(0.344 * C6$reps)
C6$hw_p50_on_detected <- round(wilson_hw(.50, C6$thinnest_detected), 4)
C6$hw_p90_on_detected <- round(wilson_hw(.90, C6$thinnest_detected), 4)
cat("\n=== 6. MONTE CARLO RESOLUTION (Wilson 95% half-width) ===\n")
cat("thinnest cell: DINA, 12.4%, HR 1.00, n 1500, detection 0.344\n")
print(C6, row.names = FALSE)

saveRDS(list(comparators = C1, fs_focus = C2, host_factor = C3, per_cell = C4,
             bounds = C5, bounds_total = tot, mc = C6),
        file.path(SCRATCH, "partB_stage0_readout.rds"))
