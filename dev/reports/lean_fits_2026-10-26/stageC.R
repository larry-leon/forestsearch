# stageC.R -- timing on the installed A+F build, and the bootstrap payloads.
# Task: dev/tasks/cc_task_lean_fits_2026-10-26.md, section 6.
# Run from this directory: Rscript stageC.R
source("timing_common.R")

pd <- packageDescription("forestsearch")[c("Version", "Packaged", "Repository", "Built")]
print(pd)
saveRDS(list(sessionInfo = sessionInfo(), packageDescription = pd),
        file.path(OUT, "stageC_sessionInfo.rds"))

tab <- run_timing("C")
print(tab)

SUMMARY <- file.path(OUT, "stageC_summary.txt")
cat(sprintf("Stage C summary -- lean fits A+F, installed forestsearch %s (Built %s)\n", pd$Version, pd$Built),
    file = SUMMARY)
sline <- function(...) { s <- paste0(...); cat(s, "\n"); cat(s, "\n", file = SUMMARY, append = TRUE, sep = "") }
failed <- character(0)
gate <- function(ok, label) {
  sline(if (isTRUE(ok)) "PASS" else "FAIL", ": ", label)
  if (!isTRUE(ok)) failed <<- c(failed, label)
}

# ---- sections 3 and 4: gate lines copied in from the committed summaries -----
for (s in c(sec3 = "identA_summary.txt", sec4 = "identF_summary.txt")) {
  for (l in grep("^(PASS|FAIL): ", readLines(s), value = TRUE))
    sline(sub("^(PASS|FAIL): ", "\\1: [", l), if (s == "identA_summary.txt") " (section 3)]" else " (section 4)]")
}
# ---- section 5: test outcomes (counts parsed from the test logs) -------------
totals <- function(f) {
  t <- grep("\\[ FAIL [0-9]+ \\| WARN [0-9]+ \\| SKIP [0-9]+ \\| PASS [0-9]+ \\]", readLines(f), value = TRUE)
  lapply(regmatches(t, gregexpr("[0-9]+", t)), as.integer)
}
g <- totals("tests_gating_log.txt")[[1]]
gate(g[1] == 0L && g[4] > 0L,
     sprintf("[section 5] test-lean-fits.R (gating): failures %d, warnings %d, skips %d, passes %d", g[1], g[2], g[3], g[4]))
ng <- totals("tests_nongating_log.txt")
for (i in seq_along(ng))
  sline(sprintf("INFO: [section 5] %s (non-gating): failures %d, warnings %d, skips %d, passes %d",
                c("test-fs-family-report.R", "test-vi-deprecation.R")[i],
                ng[[i]][1], ng[[i]][2], ng[[i]][3], ng[[i]][4]))

# ---- section 6 gates: bootstrap payloads, saved against saved ----------------
for (fx_name in names(FIXTURES)) {
  bA <- readRDS(file.path(OUT, sprintf("stageA_%s_boot.rds", fx_name)))
  bC <- readRDS(file.path(OUT, sprintf("stageC_%s_boot.rds", fx_name)))
  gate(identical(bA$names_found, bC$names_found),
       sprintf("[section 6] %s bootstrap: same components found (%s)", fx_name, paste(bC$names_found, collapse = ", ")))
  for (k in bA$names_found)
    gate(identical(bC[[k]], bA[[k]]), sprintf("[section 6] %s bootstrap: %s identical to Stage A", fx_name, k))
}

# ---- timing table ------------------------------------------------------------
tA <- read.csv(file.path(OUT, "stageA_timing.csv")); tC <- tab
pred <- c(survival = "1.51 s", continuous = "0", binary = "3.13 s")
pct <- function(a, c) 100 * (c - a) / a
for (fx_name in names(FIXTURES)) {
  a <- tA[tA$fixture == fx_name, ]; c <- tC[tC$fixture == fx_name, ]
  save_rep <- 100 * (a$rep_mean_s - c$rep_mean_s) / a$rep_mean_s
  bar <- if (fx_name %in% c("survival", "binary") && save_rep < 30) "; FINDING: below the 30% bar" else ""
  sline(sprintf(paste0("FINDING: %s | call median s: %.2f -> %.2f (%+.1f%%) | replicate mean s: %.2f -> %.2f (%+.1f%%; saving %.2f s) ",
                       "| boot outer s: %.1f -> %.1f (%+.1f%%) | Stage-0b predicted saving per replicate: %s%s"),
                fx_name, a$call_median_s, c$call_median_s, pct(a$call_median_s, c$call_median_s),
                a$rep_mean_s, c$rep_mean_s, pct(a$rep_mean_s, c$rep_mean_s), a$rep_mean_s - c$rep_mean_s,
                a$boot_outer_s, c$boot_outer_s, pct(a$boot_outer_s, c$boot_outer_s), pred[[fx_name]], bar))
}

sline(if (length(failed)) sprintf("STAGE C GATES: %d FAILED", length(failed)) else "STAGE C GATES: ALL PASS")
if (length(failed)) quit(status = 1)
