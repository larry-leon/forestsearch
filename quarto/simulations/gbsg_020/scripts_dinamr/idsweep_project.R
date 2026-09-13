# idsweep Gate 1 projection, the per-cell re-projection, and realized walls
# (TASK_idsweep_2026-09-12).
#
# usage (from scripts_dinamr/):
#   Rscript idsweep_project.R gate1                    the upfront go/no-go
#   Rscript idsweep_project.R next <i> <elapsed_s>     before cell i: exit 0 GO, 10 DEFER
#   Rscript idsweep_project.R walls                    realized against projected
#
# BASIS.  Per-replicate cost per (engine, criterion) from the pBoc smoke
# (partBoc_table.rds, 12.4% HR 1.50 n 500, MR off, 12 workers).  A cell's cost
# is that basis scaled by the engine's committed 18-cell cost profile
# (partB_stage0_readout.rds $per_cell: the MR-on field-excluded bound; the
# only committed per-cell profile) -- applied to COMPUTE ONLY.  Render overhead
# is per render and does not scale with per-replicate cost.  One render per
# cell-run: 500 replicates run as one batch, with no combine render.
#
# The consistency profile mixes hosts (14 pop-os/100-worker cells, 4
# Mac/12-worker cells); the multiplier in REPORT_partB_enabling divides the
# mixed mean by a pop-os value.  Here pop-os cells are put on the Mac scale by
# the committed host factor (mean ratio_bnd, $host_factor).  Both variants are
# printed.
#
# RE-PROJECTION.  Before each cell the cost of every remaining cell is
# re-estimated from the REALIZED sweep: per (engine, criterion), the realized
# per-replicate cost at the nearest completed anchor cell (same prevalence and
# n > same n > same prevalence > any), scaled by the profile ratio between the
# target and the anchor; overhead is the realized median per render.  Cell i is
# started only if elapsed + its projected wall <= the 13 h ceiling; otherwise
# it and every later cell are deferred (from the tail).
SCRATCH <- Sys.getenv("DINAMR_SCRATCH", unset = ".")
QMD_DIR <- Sys.getenv("DINAMR_QMD_DIR", unset = "..")
RES <- file.path(QMD_DIR, "results")
W <- 12L; REPS <- 500L; TAG <- "idsweep"
CEILING_S <- 13 * 3600; TIMEOUT_S <- 16 * 3600
args <- commandArgs(trailingOnly = TRUE); MODE <- if (length(args)) args[1] else "gate1"

cells <- read.table(file.path(SCRATCH, "idsweep.cells"), col.names = c("z", "n", "hr", "tag"),
                    colClasses = c("character", "integer", "numeric", "character"))
cells$prev <- ifelse(cells$z == "-", "12.4%", "31%")
runs <- rbind(data.frame(engine = "consistency", focus = c("effMaxSG", "effMinSG", "maxeffCons", "maxeff", "maxSG", "minSG")),
              data.frame(engine = "dina", focus = c("effMaxSG", "effMinSG", "maxSG", "minSG", "maxeffCons")),
              data.frame(engine = "grf",  focus = c("effMaxSG", "effMinSG", "maxSG", "minSG", "maxeffCons")))
runs$fo <- ifelse(runs$engine != "consistency" & runs$focus == "maxeffCons", "eff", runs$focus)
ENG <- c("consistency", "dina", "grf")

oc <- readRDS(file.path(SCRATCH, "partBoc_table.rds"))
basis <- merge(runs, oc$per_run[, c("engine", "focus", "mean_s")], by = c("engine", "focus"))
OVH0 <- oc$overhead_s
s0 <- readRDS(file.path(SCRATCH, "partB_stage0_readout.rds"))
HF <- mean(s0$host_factor$ratio_bnd)
pc <- s0$per_cell
pc$bnd_norm <- ifelse(pc$engine == "consistency" & pc$host == "pop-os", pc$bnd_worker_h / HF, pc$bnd_worker_h)
P <- function(e, prev, hr, n, norm = TRUE) {
  s <- pc[pc$engine == e & pc$prev == prev & abs(pc$hr - hr) < 1e-9 & pc$n == n, ]
  stopifnot(nrow(s) == 1L); if (norm) s$bnd_norm else s$bnd_worker_h }
Hs <- function(e, prev, hr, n)
  pc$host[pc$engine == e & pc$prev == prev & abs(pc$hr - hr) < 1e-9 & pc$n == n]
hms <- function(s) sprintf("%.2f h", s / 3600)

wall_of <- function(out) {
  lg <- file.path(SCRATCH, "logs", paste0(out, ".log"))
  if (!file.exists(lg)) return(NA_real_)
  w <- grep("^WALL_SECONDS=.* RC=0 ", readLines(lg, warn = FALSE), value = TRUE)
  if (!length(w)) NA_real_ else as.numeric(sub("^WALL_SECONDS=([0-9]+).*", "\\1", tail(w, 1)))
}
stem_of <- function(e, focus, ci) {
  z <- cells$z[ci]; band <- focus %in% c("effMaxSG", "effMinSG")
  sprintf("%s/%s_%s_fb_mr_field_m1_h%03d_knoise0_n%d%s%s_nomr_%s_res_1_%d.rds", RES,
          if (e == "consistency") "fs" else e, forestsearch::fs_focus_tag(e, focus),
          round(100 * cells$hr[ci]), cells$n[ci], if (z == "-") "" else "_z1q60",
          if (band) "_nb20" else "", TAG, REPS)
}
# Realized per run for every COMPLETED cell (all 16 renders RC 0 with bundles).
realized <- function() {
  out <- list()
  for (ci in seq_len(nrow(cells))) {
    rr <- lapply(seq_len(nrow(runs)), function(j) {
      e <- runs$engine[j]; f <- runs$focus[j]
      o <- sprintf("idsweep_%s_%s_%s", cells$tag[ci], e, runs$fo[j]); p <- stem_of(e, f, ci); w <- wall_of(o)
      if (!file.exists(p) || !is.finite(w)) return(NULL)
      x <- readRDS(p)$results$fit_mr_secs
      data.frame(ci = ci, engine = e, focus = f, mean_s = mean(x, na.rm = TRUE),
                 compute_s = sum(x, na.rm = TRUE) / W, wall_s = w)
    })
    if (all(!vapply(rr, is.null, logical(1)))) out[[length(out) + 1L]] <- do.call(rbind, rr)
  }
  if (length(out)) do.call(rbind, out) else NULL
}
# Projected wall for cell ci from anchors (realized R, or the smoke basis).
project_cell <- function(ci, R, ovh) {
  tg <- cells[ci, ]
  per <- lapply(ENG, function(e) {
    b <- basis[basis$engine == e, ]
    done <- if (is.null(R)) integer(0) else unique(R$ci[R$engine == e])
    if (length(done)) {
      sc <- vapply(done, function(a) { ca <- cells[a, ]
        if (ca$prev == tg$prev && ca$n == tg$n) 0 else if (ca$n == tg$n) 1 else if (ca$prev == tg$prev) 2 else 3 }, 0)
      pref <- done[sc == min(sc)]
      hr15 <- pref[abs(cells$hr[pref] - 1.5) < 1e-9]
      a <- if (length(hr15)) max(hr15) else max(pref)
      ra <- R[R$ci == a & R$engine == e, ]
      ms <- ra$mean_s[match(b$focus, ra$focus)]
      # Same-host profile ratio only.  The consistency profile's Mac and pop-os
      # cells scale with n very differently (Mac x1.34, pop-os x4.5 from n 500
      # to 1500), so a cross-host ratio is not a cost ratio; across hosts the
      # realized anchor is carried flat (the same-host pop-os HR 1.75 / 1.50
      # pairs sit at 1.02-1.05).
      ht <- Hs(e, tg$prev, tg$hr, tg$n); ha <- Hs(e, cells$prev[a], cells$hr[a], cells$n[a])
      ratio <- if (identical(ht, ha)) P(e, tg$prev, tg$hr, tg$n, FALSE) / P(e, cells$prev[a], cells$hr[a], cells$n[a], FALSE) else 1
      anchor <- paste0(cells$tag[a], if (identical(ht, ha)) "" else "(flat:cross-host)")
    } else {
      ms <- b$mean_s
      ratio <- P(e, tg$prev, tg$hr, tg$n) / P(e, "12.4%", 1.50, 500L)
      anchor <- "pBoc smoke"
    }
    data.frame(engine = e, anchor = anchor, ratio = ratio,
               compute_s = sum(ms) * ratio * REPS / W, overhead_s = nrow(b) * ovh)
  })
  x <- do.call(rbind, per); x$cell <- tg$tag; x$wall_s <- x$compute_s + x$overhead_s; x
}

if (MODE == "gate1") {
  cat("=== GATE 1: idsweep, 288 cell-runs, 500 replicates, MR off, 12 workers ===\n")
  cat(sprintf("Ceiling %s wall (planning bound).  Hard timeout %s (kill bound).\n\n", hms(CEILING_S), hms(TIMEOUT_S)))
  cat("--- 1. Reconciling REPORT_partB_enabling's 500-replicate projection ---\n")
  pj <- oc$projection[oc$projection$reps == 500L, ]
  ck <- readRDS(file.path(SCRATCH, "partBoc_checks.rds"))
  pj$mult <- unname(ck$multiplier[pj$engine])
  pj$as_reported <- pj$compute_h * pj$mult + pj$overhead_h
  pj$mult_on_wall <- pj$wall_h * pj$mult
  print(transform(pj[, c("engine", "compute_h", "overhead_h", "wall_h", "mult", "as_reported", "mult_on_wall")],
                  mult = round(mult, 3), as_reported = round(as_reported, 2), mult_on_wall = round(mult_on_wall, 2)),
        row.names = FALSE)
  cat(sprintf("  multiplier on compute only (the report's method): %.2f h\n", sum(pj$as_reported)))
  cat(sprintf("  multiplier on the whole wall (unrounded): %.2f h; on the rounded 3.2/1.4/2.9: %.2f h\n",
              sum(pj$mult_on_wall), sum(round(pj$wall_h, 1) * pj$mult)))
  cat("  -> 11.1 h is the arithmetic the report's method defines; ~13.6 h multiplies the per-render\n")
  cat("     overhead as well, which does not scale with per-replicate cost.  Overhead in both assumes\n")
  cat(sprintf("     2 renders per cell-run (batch + combine): %d renders x %.1f s = %.2f h.  idsweep runs one\n",
              288L * 2L, OVH0, 288 * 2 * OVH0 / 3600))
  cat(sprintf("     render per cell-run (500 replicates, one batch, no combine): 288 x %.1f s = %.2f h.\n\n",
              OVH0, 288 * OVH0 / 3600))
  cat("--- 2. The consistency multiplier mixes hosts ---\n")
  cm <- pc[pc$engine == "consistency", c("prev", "hr", "n", "host", "bnd_worker_h", "bnd_norm")]
  print(transform(cm, bnd_worker_h = round(bnd_worker_h, 2), bnd_norm = round(bnd_norm, 2)), row.names = FALSE)
  mult <- sapply(ENG, function(e) { s <- pc[pc$engine == e, ]
    c(as_reported = mean(s$bnd_worker_h) / P(e, "12.4%", 1.50, 500L, FALSE),
      host_normalized = mean(s$bnd_norm) / P(e, "12.4%", 1.50, 500L, TRUE)) })
  cat(sprintf("  host factor (mean ratio_bnd, pop-os/100 over Mac/12): %.3f\n", HF))
  print(round(mult, 3))
  cat("  DINA and GRF profiles are all-Mac, so their multipliers are unchanged.  All three profiles are\n")
  cat("  MR-ON cost (MR is ~0.85-0.95 of DINA/GRF per-replicate cost at the smoke cell), so they are an\n")
  cat("  approximation for identification-only cost; the re-projection replaces them with realized cost.\n\n")
  cat("--- 3. Per-criterion basis (smoke, s per replicate) ---\n")
  print(transform(basis[, c("engine", "focus", "fo", "mean_s")], mean_s = round(mean_s, 2)), row.names = FALSE)
  cat("\n--- 4. Corrected projection, per cell in run order (host-normalized profile; compute x profile\n")
  cat(sprintf("        ratio; 16 renders per cell x %.1f s, the 30-replicate overhead -- unmeasured at 500) ---\n", OVH0))
  PR <- do.call(rbind, lapply(seq_len(nrow(cells)), function(ci) project_cell(ci, NULL, OVH0)))
  agg <- aggregate(cbind(compute_s, overhead_s, wall_s) ~ cell, PR, sum)
  agg <- agg[match(cells$tag, agg$cell), ]; agg$cum_h <- cumsum(agg$wall_s) / 3600
  print(transform(agg, compute_s = round(compute_s), overhead_s = round(overhead_s), wall_s = round(wall_s),
                  cum_h = round(cum_h, 2)), row.names = FALSE)
  byE <- aggregate(cbind(compute_s, overhead_s, wall_s) ~ engine, PR, sum)
  byE[, -1] <- round(byE[, -1] / 3600, 2); cat("\nby engine (h):\n"); print(byE, row.names = FALSE)
  tot <- sum(PR$wall_s)
  cat(sprintf("\nTOTAL (corrected basis): %.2f h  (compute %.2f h + overhead %.2f h)\n", tot / 3600,
              sum(PR$compute_s) / 3600, sum(PR$overhead_s) / 3600))
  # sensitivities
  unn <- sum(vapply(seq_len(nrow(cells)), function(ci) {
    tg <- cells[ci, ]; sum(vapply(ENG, function(e) sum(basis$mean_s[basis$engine == e]) * REPS / W *
      P(e, tg$prev, tg$hr, tg$n, FALSE) / P(e, "12.4%", 1.50, 500L, FALSE), 0)) }, 0)) + 288 * OVH0
  cat(sprintf("  with the consistency profile as reported (not host-normalized): %.2f h\n", unn / 3600))
  SENS <- data.frame(overhead_s_per_render = c(OVH0, 30, 60, 90, 120))
  SENS$total_h <- round((sum(PR$compute_s) + 288 * SENS$overhead_s_per_render) / 3600, 2)
  SENS$cells_within_13h <- vapply(SENS$overhead_s_per_render, function(o)
    sum(cumsum(agg$compute_s + 16 * o) <= CEILING_S), integer(1))
  cat("  render-overhead sensitivity (the 500-replicate MR-off render is unmeasured):\n")
  print(SENS, row.names = FALSE)
  fit <- sum(agg$cum_h * 3600 <= CEILING_S)
  cat(sprintf("\nGATE 1: %d of 18 cells fit under the %s ceiling on the corrected basis -> %s\n", fit,
              hms(CEILING_S), if (fit == 18L) "GO, all 288" else sprintf("GO, with cells %d-18 projected to defer", fit + 1L)))
  cat("Re-projected before every cell from realized cost; the first cell measures the 500-replicate render overhead.\n")
  saveRDS(list(reconcile = pj, multipliers = mult, host_factor = HF, basis = basis, per_engine_cell = PR,
               per_cell = agg, total_h = tot / 3600, total_unnormalized_h = unn / 3600, sensitivity = SENS,
               ceiling_s = CEILING_S, timeout_s = TIMEOUT_S), file.path(SCRATCH, "idsweep_gate1.rds"))

} else if (MODE == "next") {
  ci <- as.integer(args[2]); el <- as.numeric(args[3])
  R <- realized()
  if (!is.null(R)) {
    ovr <- R$wall_s - R$compute_s; ovh <- median(ovr)
    src <- sprintf("realized median over %d renders (%d cells)", nrow(R), length(unique(R$ci)))
  } else { ovh <- OVH0; src <- "smoke's 30-replicate median (no cell realized yet)" }
  cat(sprintf("=== RE-PROJECTION before cell %d (%s) at elapsed %s ===\n", ci, cells$tag[ci], hms(el)))
  cat(sprintf("render overhead: %.1f s per render -- %s\n", ovh, src))
  if (!is.null(R)) {
    re <- aggregate(cbind(compute_s, wall_s) ~ ci, R, sum)
    cat("realized cells (s):\n"); print(data.frame(cell = cells$tag[re$ci], compute_s = round(re$compute_s),
                                                   wall_s = re$wall_s), row.names = FALSE)
  }
  rem <- ci:nrow(cells)
  PR <- do.call(rbind, lapply(rem, function(k) project_cell(k, R, ovh)))
  agg <- aggregate(cbind(compute_s, overhead_s, wall_s) ~ cell, PR, sum)
  agg <- agg[match(cells$tag[rem], agg$cell), ]
  agg$finish_h <- (el + cumsum(agg$wall_s)) / 3600
  anc <- aggregate(anchor ~ cell, PR, function(a) paste(unique(a), collapse = "/"))
  agg$anchors <- anc$anchor[match(agg$cell, anc$cell)]
  print(transform(agg, compute_s = round(compute_s), overhead_s = round(overhead_s), wall_s = round(wall_s),
                  finish_h = round(finish_h, 2)), row.names = FALSE)
  nxt <- agg$wall_s[1]; go <- el + nxt <= CEILING_S
  fits <- sum(agg$finish_h * 3600 <= CEILING_S)
  cat(sprintf("projected: all remaining finish at %.2f h; %d of %d remaining fit under the ceiling in order\n",
              tail(agg$finish_h, 1), fits, length(rem)))
  f <- file.path(SCRATCH, "idsweep_reproject.rds")
  hist <- if (file.exists(f)) readRDS(f) else list()
  hist[[cells$tag[ci]]] <- list(elapsed_s = el, overhead_s = ovh, projected = agg, per_engine = PR, go = go,
                                at = Sys.time())
  saveRDS(hist, f)
  cat(sprintf("DECISION cell %d %s: elapsed %s + projected %.0f s = %s vs ceiling %s -> %s\n", ci, cells$tag[ci],
              hms(el), nxt, hms(el + nxt), hms(CEILING_S), if (go) "GO" else "DEFER (this cell and every later cell)"))
  quit(status = if (go) 0L else 10L)

} else if (MODE == "walls") {
  R <- realized(); stopifnot(!is.null(R))
  g1 <- readRDS(file.path(SCRATCH, "idsweep_gate1.rds"))
  hist <- if (file.exists(file.path(SCRATCH, "idsweep_reproject.rds"))) readRDS(file.path(SCRATCH, "idsweep_reproject.rds")) else list()
  dl <- file.path(SCRATCH, "logs", "idsweep.driver.log")
  cd <- if (file.exists(dl)) grep("^CELL DONE: ", readLines(dl, warn = FALSE), value = TRUE) else character(0)
  cw <- setNames(as.numeric(sub(".*wall=([0-9]+)s.*", "\\1", cd)), sub("^CELL DONE: (\\S+).*", "\\1", cd))
  re <- aggregate(cbind(compute_s, wall_s) ~ ci, R, sum)
  re$cell <- cells$tag[re$ci]
  re$overhead_s <- re$wall_s - re$compute_s
  re$cell_done_wall_s <- unname(cw[re$cell])
  re$gate1_proj_s <- g1$per_cell$wall_s[match(re$cell, g1$per_cell$cell)]
  re$reproj_s <- vapply(re$cell, function(k) if (is.null(hist[[k]])) NA_real_ else hist[[k]]$projected$wall_s[1], 0)
  re$real_over_gate1 <- re$cell_done_wall_s / re$gate1_proj_s
  re$real_over_reproj <- re$cell_done_wall_s / re$reproj_s
  print(transform(re[, c("cell", "compute_s", "overhead_s", "wall_s", "cell_done_wall_s", "gate1_proj_s", "reproj_s",
                          "real_over_gate1", "real_over_reproj")],
                  compute_s = round(compute_s), overhead_s = round(overhead_s), gate1_proj_s = round(gate1_proj_s),
                  reproj_s = round(reproj_s), real_over_gate1 = round(real_over_gate1, 3),
                  real_over_reproj = round(real_over_reproj, 3)), row.names = FALSE)
  byE <- aggregate(cbind(compute_s, wall_s) ~ engine, R, sum); byE$overhead_s <- byE$wall_s - byE$compute_s
  g1E <- aggregate(cbind(compute_s, overhead_s, wall_s) ~ engine, g1$per_engine_cell[g1$per_engine_cell$cell %in% re$cell, ], sum)
  cat("\nby engine over the completed cells (h): realized vs Gate 1\n")
  print(data.frame(engine = byE$engine, real_compute_h = round(byE$compute_s / 3600, 2),
                   real_overhead_h = round(byE$overhead_s / 3600, 2), real_wall_h = round(byE$wall_s / 3600, 2),
                   gate1_compute_h = round(g1E$compute_s / 3600, 2), gate1_wall_h = round(g1E$wall_s / 3600, 2)),
        row.names = FALSE)
  ovr <- R$wall_s - R$compute_s
  cat(sprintf("\nrender overhead per render, realized: median %.1f s [%.0f, %.0f] over %d renders\n",
              median(ovr), min(ovr), max(ovr), length(ovr)))
  pcrit <- aggregate(mean_s ~ engine + focus, R, mean)
  cat("per-criterion mean s per replicate, averaged over completed cells:\n"); print(transform(pcrit, mean_s = round(mean_s, 2)), row.names = FALSE)
  cat(sprintf("\nTOTAL: realized %.2f h (sum of CELL DONE walls, %d cells) vs Gate 1 %.2f h for the same cells; Gate 1 all-18 %.2f h\n",
              sum(re$cell_done_wall_s, na.rm = TRUE) / 3600, nrow(re), sum(re$gate1_proj_s) / 3600, g1$total_h))
  saveRDS(list(per_cell = re, per_run = R, by_engine = byE), file.path(SCRATCH, "idsweep_walls.rds"))
}
