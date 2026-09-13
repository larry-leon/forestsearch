# Part B OC table (TASK_partB_enabling_2026-09-12, Part 2) from the sixteen
# `pBoc` smoke bundles.  Reads bundles and render logs only.
#
# Definitions (all rates over the run's 30 replicates unless stated):
#   rate        detected / 30, Wilson 95%.  At HR 1.50 this is a detection rate.
#   sens..npv   replicate means over DETECTED replicates of the template's
#               per-replicate .classify() rates (the form the DINA tables use;
#               no interval on a mean of proportions).
#   |Hhat|,|H|  mean n_sel and mean n_true over detected replicates; ratio of
#               the two means.
#   n_family    NA with MR off by construction (MR's own fitted family; see the
#               template's record_replicate() comment).  Printed as NA.
#   admitted_n  GRF only: median [min, max] over replicates where recorded.
#   prevalence  mean(n_true / n) over all 30 replicates.
#   wall        median per-replicate forestsearch() seconds (fit_mr_secs, which
#               with MR off times the identification alone); the run's total
#               wall is the render's WALL_SECONDS from its log.
#
# usage (from scripts_dinamr/): Rscript partBoc_table.R
SCRATCH <- Sys.getenv("DINAMR_SCRATCH", unset = ".")
QMD_DIR <- Sys.getenv("DINAMR_QMD_DIR", unset = "..")
RES <- file.path(QMD_DIR, "results")
W <- 12L; N <- 500L; REPS <- 30L

runs <- rbind(
  data.frame(engine = "consistency",
             focus = c("effMaxSG", "effMinSG", "maxeffCons", "maxeff", "maxSG", "minSG")),
  data.frame(engine = "dina", focus = c("effMaxSG", "effMinSG", "maxSG", "minSG", "maxeffCons")),
  data.frame(engine = "grf",  focus = c("effMaxSG", "effMinSG", "maxSG", "minSG", "maxeffCons")))
runs$tag  <- mapply(function(e, f) forestsearch::fs_focus_tag(e, f), runs$engine, runs$focus)
runs$band <- runs$focus %in% c("effMaxSG", "effMinSG")
stem <- function(e, tag, band) sprintf("%s_%s_fb_mr_field_m1_h150_knoise0_n500%s_nomr_pBoc_res_1_30.rds",
                                       if (e == "consistency") "fs" else e, tag, if (band) "_nb20" else "")
wilson <- function(x, n) { z <- qnorm(.975); p <- x / n
  c((p + z^2/(2*n) - z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n),
    (p + z^2/(2*n) + z*sqrt(p*(1-p)/n + z^2/(4*n^2))) / (1 + z^2/n)) }
wall_of <- function(e, f) {
  lg <- file.path(SCRATCH, "logs", sprintf("partBoc_%s_%s.log", e, f))
  if (!file.exists(lg)) return(NA_real_)
  w <- grep("^WALL_SECONDS=", readLines(lg), value = TRUE)
  if (!length(w)) NA_real_ else as.numeric(sub("^WALL_SECONDS=([0-9]+).*", "\\1", tail(w, 1)))
}
mrng <- function(x) { x <- x[is.finite(x)]
  if (!length(x)) "NA" else sprintf("%g [%g, %g]", median(x), min(x), max(x)) }

B <- list(); rows <- list()
for (i in seq_len(nrow(runs))) {
  e <- runs$engine[i]; fo <- runs$focus[i]
  p <- file.path(RES, stem(e, runs$tag[i], runs$band[i]))
  if (!file.exists(p)) { cat("MISSING:", basename(p), "\n"); next }
  b <- readRDS(p); r <- b$results; B[[paste(e, fo)]] <- r
  stopifnot(nrow(r) == REPS, identical(b$meta$mr_inference, FALSE))
  d <- r$detected %in% 1L; nd <- sum(d); wi <- wilson(nd, REPS)
  m <- function(cc) if (nd) mean(r[[cc]][d], na.rm = TRUE) else NA_real_
  rows[[length(rows) + 1L]] <- data.frame(
    engine = e,
    sg_focus = if (e != "consistency" && fo == "maxeffCons") "eff (= maxeff = maxeffCons)" else fo,
    eps = if (runs$band[i]) sprintf("%.2f", b$meta$effect_neighborhood) else "inert",
    rate = sprintf("%.3f [%.3f, %.3f]", nd / REPS, wi[1], wi[2]),
    sens = round(m("sens"), 3), spec = round(m("spec"), 3),
    ppv = round(m("ppv"), 3), npv = round(m("npv"), 3),
    mean_Hhat = round(m("n_sel"), 1), mean_H = round(m("n_true"), 1),
    ratio = round(m("n_sel") / m("n_true"), 3),
    n_family = mrng(r$n_family),
    admitted_n = if (e == "grf") mrng(r$admitted_n) else "--",
    prevalence = round(mean(r$n_true / N), 4),
    med_rep_s = round(median(r$fit_mr_secs, na.rm = TRUE), 2),
    total_wall_s = wall_of(e, fo),
    status = paste(names(table(r$status)), table(r$status), sep = ":", collapse = " "),
    stringsAsFactors = FALSE)
}
OC <- do.call(rbind, rows)
cat("=== THE OC TABLE: 12.4%, HR 1.50, n 500, 30 replicates, MR off ===\n")
op <- options(width = 250); print(OC, row.names = FALSE); options(op)

cat("\n=== NA COLUMNS ===\n")
for (k in names(B)) {
  r <- B[[k]]; d <- r$detected %in% 1L
  idc <- c("label", "sg_def", "n_sel", "n_cons_qual", "band_n", "admitted_n", "n_family",
           "sens", "spec", "ppv", "npv", "betaHhat_H")
  na <- idc[vapply(idc, function(cc) any(d) && all(is.na(r[[cc]][d])), logical(1))]
  cat(sprintf("%-24s all-NA on detected rows: %s\n", k, if (length(na)) paste(na, collapse = ", ") else "<none>"))
}

cat("\n=== SAME SUBGROUP ACROSS CRITERIA, SAME SEEDS (within engine; sg_def string match, both detected) ===\n")
AG <- list()
for (e in unique(runs$engine)) {
  ks <- grep(paste0("^", e, " "), names(B), value = TRUE)
  M <- matrix(NA_integer_, length(ks), length(ks), dimnames = list(sub("^\\S+ ", "", ks), sub("^\\S+ ", "", ks)))
  for (a in seq_along(ks)) for (bb in seq_along(ks)) {
    ra <- B[[ks[a]]]; rb <- B[[ks[bb]]]
    stopifnot(identical(ra$sim_id, rb$sim_id))
    both <- ra$detected %in% 1L & rb$detected %in% 1L
    M[a, bb] <- sum(both & ra$sg_def == rb$sg_def, na.rm = TRUE)
  }
  cat(sprintf("\n%s -- cell [i,j] = replicates (of 30) where both selected the identical rule:\n", e))
  print(M); AG[[e]] <- M
  for (a in seq_along(ks)) for (bb in seq_along(ks)) if (a < bb) {
    ra <- B[[ks[a]]]; rb <- B[[ks[bb]]]
    both <- ra$detected %in% 1L & rb$detected %in% 1L
    if (sum(both) && M[a, bb] == sum(both))
      cat(sprintf("  NEVER DIFFER: %s vs %s -- identical on all %d jointly detected replicates\n",
                  rownames(M)[a], rownames(M)[bb], sum(both)))
  }
}
# The truth is the same data on every run (same seeds, same cell).
nt <- vapply(B, function(r) paste(r$n_true, collapse = ","), character(1))
cat(sprintf("\nn_true identical across all %d runs (same DGM draws): %s\n", length(B), length(unique(nt)) == 1L))

cat("\n=== PROJECTION: 288 cell-runs from the measured per-replicate seconds ===\n")
cat("Basis: this ONE cell (12.4%, HR 1.50, n 500).  Compute = sum over the engine's criteria of\n")
cat("mean per-replicate seconds x reps x 18 cells / 12 workers.  Render overhead per render =\n")
cat("median over the 16 runs of (WALL_SECONDS - sum(fit_mr_secs)/12); renders per cell-run =\n")
cat("ceiling(reps/1000) batches + 1 combine.\n\n")
per <- do.call(rbind, lapply(names(B), function(k) {
  r <- B[[k]]; e <- sub(" .*", "", k); fo <- sub("^\\S+ ", "", k)
  data.frame(engine = e, focus = fo, mean_s = mean(r$fit_mr_secs, na.rm = TRUE),
             med_s = median(r$fit_mr_secs, na.rm = TRUE),
             compute_s = sum(r$fit_mr_secs, na.rm = TRUE) / W, wall_s = wall_of(e, fo))
}))
per$overhead_s <- per$wall_s - per$compute_s
ovh <- median(per$overhead_s, na.rm = TRUE)
print(transform(per, mean_s = round(mean_s, 2), med_s = round(med_s, 2),
                compute_s = round(compute_s, 1), overhead_s = round(overhead_s, 1)), row.names = FALSE)
cat(sprintf("\nmedian per-render overhead: %.1f s\n", ovh))
PJ <- do.call(rbind, lapply(c(2000L, 1000L, 500L), function(R) {
  do.call(rbind, lapply(c("consistency", "dina", "grf"), function(e) {
    s <- per[per$engine == e, ]
    renders <- (ceiling(R / 1000) + 1L) * 18L * nrow(s)
    comp_h <- sum(s$mean_s) * R * 18 / W / 3600
    data.frame(reps = R, engine = e, cell_runs = 18L * nrow(s),
               compute_h = round(comp_h, 2), overhead_h = round(renders * ovh / 3600, 2),
               wall_h = round(comp_h + renders * ovh / 3600, 2))
  }))
}))
print(PJ, row.names = FALSE)
TOT <- aggregate(cbind(compute_h, overhead_h, wall_h) ~ reps, PJ, sum)
cat("\nAll three engines:\n"); print(TOT[order(-TOT$reps), ], row.names = FALSE)

saveRDS(list(oc = OC, agreement = AG, per_run = per, overhead_s = ovh,
             projection = PJ, projection_total = TOT),
        file.path(SCRATCH, "partBoc_table.rds"))
