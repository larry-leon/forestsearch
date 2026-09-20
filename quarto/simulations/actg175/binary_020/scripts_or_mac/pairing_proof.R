# Cross-machine pairing proof for TASK_binary_stage2_v3_2026-09-19's report.
# For a Mac orgrf/ordina cell, compare the DATA-ONLY columns per replicate against
# the pop-os orfs bundle at the SAME design point and n: seed and n_true (the true-
# region size) and the 8 oracle columns (the oracle refits on the TRUE region), all
# rule-independent by construction.  Tolerance-based, not identical(): the machines
# differ in BLAS/LAPACK and R version, so the same draws agree to floating point,
# not bitwise.  Reports the matched fraction, as scripts_mdf1's pairing proof did.
# usage: Rscript pairing_proof.R <qmd-dir> <campaign: orgrf|ordina> <dtok: or075|or100|or150> <n>
a <- commandArgs(TRUE); QD <- a[1]; TAG <- a[2]; DT <- a[3]; N <- as.integer(a[4])
TOL <- 1e-8
mtag <- function(t) switch(t, orgrf = "grf", ordina = "dina", orfs = "fs")
bun <- function(tag) { st <- sprintf("%s_effMaxSG_mr_field_%s_n%d_nb20_%s", mtag(tag), DT, N, tag)
  file.path(QD, "mr_or_harm", paste0(st, "_d5000"), sprintf("%s_combined_1_1000.rds", st)) }
fm <- bun(TAG); fp <- bun("orfs")
if (!file.exists(fm) || !file.exists(fp)) { cat(sprintf("SKIP %s_%s_n%d: bundle missing\n", TAG, DT, N)); quit(status = 0) }
m <- readRDS(fm); p <- readRDS(fp); r <- m$results; f <- p$results
o <- order(r$sim_id); of <- order(f$sim_id)
COLS <- c("seed","n_true","or_H_est","or_H_lo","or_H_hi","or_H_se","or_Hc_est","or_Hc_lo","or_Hc_hi","or_Hc_se")
cat(sprintf("\n=== %s_%s_n%d  vs  orfs_%s_n%d ===\n", TAG, DT, N, DT, N))
cat(sprintf("  Mac    : host %s, R %s, pkg %s\n", m$meta$hostname, m$meta$r_version, m$meta$pkg_version))
cat(sprintf("  pop-os : host %s, R %s, pkg %s\n", p$meta$hostname, p$meta$r_version, p$meta$pkg_version))
if (!(nrow(r) == nrow(f) && identical(sort(r$sim_id), sort(f$sim_id)))) {
  cat("  ROW SETS DIFFER -- pairing is by construction (seed table) and NOT verified here\n"); quit(status = 0) }
tot <- 0L; mx_all <- 0
for (k in COLS) {
  x <- as.numeric(r[[k]][o]); y <- as.numeric(f[[k]][of])
  d <- abs(x - y) / pmax(abs(y), 1e-300); d[is.na(x) & is.na(y)] <- 0; d[xor(is.na(x), is.na(y))] <- Inf
  cat(sprintf("  %-10s matched %4d/%4d (%.4f)  max rel %.3g\n", k, sum(d <= TOL), length(d), mean(d <= TOL), max(d)))
  tot <- tot + sum(d <= TOL); mx_all <- max(mx_all, max(d))
}
cat(sprintf("  ALL %d data-only columns: matched %d/%d (%.4f) at tol %g; max rel %.3g\n",
            length(COLS), tot, length(COLS) * nrow(r), tot / (length(COLS) * nrow(r)), TOL, mx_all))
cat(sprintf("  seed identical() on all rows: %s | or_H_est identical(): %s (floating point, not bitwise)\n",
            identical(as.integer(r$seed[o]), as.integer(f$seed[of])),
            identical(as.numeric(r$or_H_est[o]), as.numeric(f$or_H_est[of]))))
cat(sprintf("  discriminating contrast -- declaration rate: %s %.4f vs orfs %.4f\n",
            TAG, mean(r$detected %in% 1L), mean(f$detected %in% 1L)))
