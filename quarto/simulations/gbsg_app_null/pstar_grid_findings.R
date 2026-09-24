# =============================================================================
# pstar_grid_findings.R -- read-out and post-condition checks for
# TASK_gbsg_pstar_grid_2026-09-23 (driver run_gbsg_pstar_grid.R)
#
# Usage (from the package root):
#   Rscript quarto/simulations/gbsg_app_null/pstar_grid_findings.R \
#     > quarto/simulations/gbsg_app_null/logs/pstar_grid_findings.txt
#
# Reads only the committed payloads in results/.  Pcons is reported to two
# decimals (pconsistency.digits = 2), so the declaration-rate CDF is evaluated
# on the 0.01 grid, compared on the integer percent scale to avoid
# floating-point edge effects: rate(p) = mean(round(100 * max Pcons) >= 100 p),
# with a replicate whose out_sg is empty counted as not declaring.
# =============================================================================

rd <- "quarto/simulations/gbsg_app_null/results"
gA <- readRDS(file.path(rd, "gbsg_pstar_gateA.rds"))
gB <- readRDS(file.path(rd, "gbsg_pstar_gateB.rds"))
gC <- readRDS(file.path(rd, "gbsg_pstar_gateC.rds"))
cl <- lapply(1:2, function(k) readRDS(file.path(rd, sprintf("gbsg_pstar_cell%d.rds", k))))

wilson <- function(x, n, alpha = 0.10) {
  z <- stats::qnorm(1 - alpha / 2); p <- x / n
  den <- 1 + z^2 / n
  ctr <- (p + z^2 / (2 * n)) / den
  hw <- z / den * sqrt(p * (1 - p) / n + z^2 / (4 * n^2))
  c(lo = max(0, ctr - hw), hi = min(1, ctr + hw))
}
pc_int <- function(r) ifelse(is.na(r$max_pcons), -1L, as.integer(round(100 * r$max_pcons)))
rate_at <- function(r, p) {
  x <- sum(pc_int(r) >= round(100 * p)); n <- nrow(r)
  ci <- wilson(x, n)
  c(p = p, count = x, n = n, rate = x / n, lo = ci[["lo"]], hi = ci[["hi"]])
}
# FW_0.10(0.75) of the GBSG application fit at the aligned (rounded) screen,
# 0.6672: REPORT_declcal_rounding_alignment_2026-09-23 Gate 2 (ac860c6d),
# re-checked as mean(Mstar_c0 > z_pstar) on that report's capture
# (TASK_declcal_consumers_2026-09-24_v3). The application payload in
# fs-glms-interpretable still carries the exact-cutoff value until it re-runs.
fw_app_075 <- 0.6672
pc <- list()
fail <- function(id, ok, note = "") {
  pc[[id]] <<- ok
  cat(sprintf("  [%s] %-4s %s\n", if (ok) "PASS" else "FAIL", id, note))
}

cat("=== Pins ===\n")
for (k in 1:2) cat(sprintf("cell %d: git_head %s | forestsearch %s built %s | %s | workers %d\n",
                           k, substr(cl[[k]]$git_head, 1, 8), cl[[k]]$forestsearch_version,
                           cl[[k]]$forestsearch_built, cl[[k]]$R, cl[[k]]$n_workers))
cat(sprintf("gates: A %s | B %s | C %s\n", substr(gA$git_head, 1, 8),
            substr(gB$git_head, 1, 8), substr(gC$git_head, 1, 8)))

cat("\n=== Truth (uncensored population marginal Cox HR, on df_source) ===\n")
cat(sprintf("k_treat %.6f (recalibrated: %s) | overall %.4f | patient-level HR %.4f\n",
            gA$k_treat_used, !is.null(gA$recalibration), gA$hr_source, gA$hr_patient))
h <- gA$family$hr_true_uncensored
cat(sprintf("family M = %d (task route) | MR route %d vs application 1,744\n",
            length(h), gA$family$counts[["mr_route"]]))
print(gA$family$counts)
cat(sprintf("per-candidate HR: min %.4f median %.4f max %.4f | > 0.75: %d of %d (%.1f%%)\n",
            min(h), stats::median(h), max(h), sum(h > 0.75), length(h), 100 * mean(h > 0.75)))

cat("\n=== Cost ===\n")
cores <- parallel::detectCores()
walls <- vapply(cl, `[[`, numeric(1), "wall_seconds")
cost <- data.frame(
  run = c("Gate B (p* 0.90 + floor)", "Gate C (floor)", "Cell 1 (floor)", "Cell 2 (floor)"),
  replicates = c(2L * nrow(gB$results_floor), nrow(gC$results),
                 nrow(cl[[1]]$results), nrow(cl[[2]]$results)),
  workers = c(gB$n_workers, gC$n_workers, cl[[1]]$n_workers, cl[[2]]$n_workers),
  wall_s = round(c(gB$wall_090 + gB$wall_floor, gC$wall_seconds, walls), 1))
cost$wall_min <- round(cost$wall_s / 60, 2)
cost$s_per_rep_per_worker <- round(cost$wall_s * cost$workers / cost$replicates, 2)
print(cost, row.names = FALSE)
cat(sprintf("machine %d logical cores; cells total %.1f s = %.1f min for %d replicates\n",
            cores, sum(walls), sum(walls) / 60, sum(cost$replicates[3:4])))
cat(sprintf("mean per-replicate fit time: cell 1 %.2f s, cell 2 %.2f s\n",
            mean(cl[[1]]$results$seconds), mean(cl[[2]]$results$seconds)))

grid <- c(0.90, 0.95, 0.96, 0.97, 0.98, 0.99)
steps <- seq(0.50, 1.00, by = 0.01)
for (k in 1:2) {
  r <- cl[[k]]$results
  cat(sprintf("\n=== Cell %d (%s) ===\n", k,
              if (k == 1) "draw_treatment = TRUE, rand_ratio = 246/440"
              else "draw_treatment = FALSE"))
  cat(sprintf("replicates %d | errored %d | denominator %d\n",
              nrow(r), sum(r$errored), nrow(r)))
  if (any(r$errored)) print(r[r$errored, c("sim", "error_msg")])
  cat(sprintf("treated: mean fraction %.4f | range %d-%d | event rate mean %.4f (SD %.4f)\n",
              mean(r$n_treated / r$n), min(r$n_treated), max(r$n_treated),
              mean(r$event_rate), stats::sd(r$event_rate)))
  cdf <- t(vapply(steps, function(p) rate_at(r, p), numeric(6)))
  cat("CDF of declaration rate vs p* (0.01 steps):\n")
  print(data.frame(p = sprintf("%.2f", cdf[, "p"]), count = cdf[, "count"],
                   rate = sprintf("%.4f", cdf[, "rate"])), row.names = FALSE)
  first <- which(cdf[, "rate"] <= 0.10)[1]
  if (is.na(first)) {
    cat("rate never falls to 0.10 on [0.50, 1.00]\n")
  } else {
    br <- cdf[max(1, first - 1):first, , drop = FALSE]
    cat(sprintf("first p* with rate <= 0.10: %.2f (rate %.4f, 90%% Wilson [%.4f, %.4f], %d / %d)\n",
                cdf[first, "p"], cdf[first, "rate"], cdf[first, "lo"], cdf[first, "hi"],
                cdf[first, "count"], nrow(r)))
    cat("bracketing 0.01 grid points:\n")
    print(round(br, 4))
  }
  cat(sprintf("rate at 0.90: %.4f vs FW_0.10(0.75) %.4f and super-population 0.329\n",
              rate_at(r, 0.90)[["rate"]], fw_app_075))
}

cat("\n=== Grid table (cells x p*; rate [90% Wilson] count/B) ===\n")
tab <- do.call(rbind, lapply(1:2, function(k) {
  r <- cl[[k]]$results
  vapply(grid, function(p) {
    x <- rate_at(r, p)
    sprintf("%.4f [%.4f, %.4f] %d/%d", x[["rate"]], x[["lo"]], x[["hi"]],
            x[["count"]], x[["n"]])
  }, character(1))
}))
dimnames(tab) <- list(c("Cell 1", "Cell 2"), sprintf("%.2f", grid))
print(t(tab), quote = FALSE)

cat("\n=== Consistency with the gates ===\n")
f <- gB$results_floor
cols <- c("sim", "declared", "max_pcons", "n_candidates_total", "n_passed",
          "sg.def", "n_treated", "events")
same200 <- isTRUE(all.equal(f[, cols], cl[[1]]$results[seq_len(nrow(f)), cols],
                            check.attributes = FALSE))
cat(sprintf("cell 1 rows 1-%d identical to the Gate B floor run: %s\n", nrow(f), same200))
same10 <- isTRUE(all.equal(gC$results[, cols], cl[[1]]$results[seq_len(nrow(gC$results)), cols],
                           check.attributes = FALSE))
cat(sprintf("cell 1 rows 1-%d identical to the Gate C run: %s\n", nrow(gC$results), same10))
seq_vs_full <- isTRUE(all.equal(gC$seq_rerun[, cols], cl[[1]]$results[1, cols],
                                check.attributes = FALSE))
cat(sprintf("sequential L'Ecuyer re-run of replicate 1 equals cell 1 replicate 1: %s\n",
            seq_vs_full))

cat("\n=== Post-conditions ===\n")
fail("PC1", isTRUE(gA$gateA_pass) && abs(gA$hr_source / 0.75 - 1) < 0.01 &&
       isTRUE(gA$gate[["flag_harm_zero_source"]]),
     sprintf("df_source marginal HR %.4f (%+.2f%%), flag_harm identically 0",
             gA$hr_source, 100 * (gA$hr_source / 0.75 - 1)))
set_all <- c(gB$results_090$settings_ok, gB$results_floor$settings_ok,
             gC$results$settings_ok, cl[[1]]$results$settings_ok,
             cl[[2]]$results$settings_ok)
fail("PC2", all(set_all %in% TRUE),
     sprintf("args_call_all settings incl. stop_threshold = NULL on %d / %d fits",
             sum(set_all %in% TRUE), length(set_all)))
fail("PC3", all(gB$checks[c("exact_agreement", "n_candidates_match")]),
     sprintf("Gate B: %d disagreements, %d n_candidates_total mismatches over %d",
             length(gB$disagree), length(gB$ncand_mismatch), nrow(gB$results_floor)))
fail("PC4", !is.null(gC$wall_seconds) && all(walls < 7200),
     sprintf("Gate C wall %.1f s recorded; cells %.0f s and %.0f s < 7,200 s",
             gC$wall_seconds, walls[1], walls[2]))
fail("PC5", all(vapply(cl, function(x) sum(x$results$errored) == 0L &&
                         nrow(x$results) == 5000L, logical(1))),
     sprintf("errors %d and %d; denominators %d and %d",
             sum(cl[[1]]$results$errored), sum(cl[[2]]$results$errored),
             nrow(cl[[1]]$results), nrow(cl[[2]]$results)))
fail("PC6", all(cl[[2]]$results$n_treated == 246L) &&
       abs(mean(cl[[1]]$results$n_treated / 686) - 0.359) <= 0.01,
     sprintf("cell 2 treated == 246 on %d / %d; cell 1 mean fraction %.4f",
             sum(cl[[2]]$results$n_treated == 246L), nrow(cl[[2]]$results),
             mean(cl[[1]]$results$n_treated / 686)))
fail("PC7", isTRUE(gC$repro_seq_lecuyer) && seq_vs_full,
     "replicate 1 re-run sequentially under L'Ecuyer-CMRG equals the parallel result")
ev <- vapply(cl, function(x) mean(x$results$event_rate), numeric(1))
fail("PC8", all(abs(ev - 0.436) <= 0.02),
     sprintf("event rate %.4f / %.4f vs 0.436", ev[1], ev[2]))
cat("  PC9 (catalogue pin = HEAD) and PC10 (write locations, no R/ change) are checked at commit time.\n")
cat(sprintf("\nALL COMPUTABLE POST-CONDITIONS PASS: %s\n", all(unlist(pc))))
