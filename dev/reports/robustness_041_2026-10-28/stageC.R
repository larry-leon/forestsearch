# stageC.R -- identity against Stage A, records and timing on the installed
# robustness build.  Task: dev/tasks/cc_task_robustness_041_2026-10-28.md, section 8.
# Transplant of ../lean_fits_2026-10-26/stageC.R (sections-3/4 copy-in and the
# prediction table dropped; gate/summary mechanics kept).
# Run from this directory: Rscript stageC.R
source("timing_common.R")

pd <- packageDescription("forestsearch")[c("Version", "Packaged", "Repository", "Built")]
print(pd)
saveRDS(list(sessionInfo = sessionInfo(), packageDescription = pd),
        file.path(OUT, "stageC_sessionInfo.rds"))

tab <- run_timing("C")
print(tab)

SUMMARY <- file.path(OUT, "stageC_summary.txt")
cat(sprintf("Stage C summary -- robustness 0.4.1, installed forestsearch %s (Built %s)\n", pd$Version, pd$Built),
    file = SUMMARY)
sline <- function(...) { s <- paste0(...); cat(s, "\n"); cat(s, "\n", file = SUMMARY, append = TRUE, sep = "") }
failed <- character(0)
gate <- function(ok, label) {
  sline(if (isTRUE(ok)) "PASS" else "FAIL", ": ", label)
  if (!isTRUE(ok)) failed <<- c(failed, label)
}

# ---- 1. end to end, saved against saved ---------------------------------------
E2E <- c("sg.harm", "membership", "out_sg_result", "hr.subgroups", "filter_counts",
         "n_candidates", "n_fits", "n_warnings", "warnings", "args_call_all")
for (fx_name in names(FIXTURES)) {
  a <- readRDS(file.path(OUT, sprintf("stageA_%s_ii.rds", fx_name)))
  c <- readRDS(file.path(OUT, sprintf("stageC_%s_ii.rds", fx_name)))
  for (k in E2E)
    gate(identical(c[[k]], a[[k]]), sprintf("[1] %s end to end: %s identical to Stage A", fx_name, k))
}

# ---- 2. bootstrap payloads, saved against saved -------------------------------
for (fx_name in names(FIXTURES)) {
  bA <- readRDS(file.path(OUT, sprintf("stageA_%s_boot.rds", fx_name)))
  bC <- readRDS(file.path(OUT, sprintf("stageC_%s_boot.rds", fx_name)))
  gate(identical(bA$names_found, bC$names_found),
       sprintf("[2] %s bootstrap: same components found (%s)", fx_name, paste(bC$names_found, collapse = ", ")))
  for (k in bA$names_found)
    gate(identical(bC[[k]], bA[[k]]), sprintf("[2] %s bootstrap: %s identical to Stage A", fx_name, k))
}

# ---- 3. test outcomes (counts parsed from the test logs) ----------------------
totals <- function(f) {
  t <- grep("\\[ FAIL [0-9]+ \\| WARN [0-9]+ \\| SKIP [0-9]+ \\| PASS [0-9]+ \\]", readLines(f), value = TRUE)
  as.integer(regmatches(t[length(t)], gregexpr("[0-9]+", t[length(t)]))[[1]])
}
for (f in c("tests_gating_lean_fits.txt", "tests_gating_iptw.txt")) {
  g <- totals(f)
  gate(g[1] == 0L && g[4] > 0L,
       sprintf("[3] %s (gating): failures %d, warnings %d, skips %d, passes %d", f, g[1], g[2], g[3], g[4]))
}
for (f in c("tests_nongating_family_report.txt", "tests_nongating_vi_deprecation.txt")) {
  g <- totals(f)
  sline(sprintf("INFO: [3] %s (non-gating): failures %d, warnings %d, skips %d, passes %d", f, g[1], g[2], g[3], g[4]))
}

# ---- 4. reproduce: before (copied in) and after (re-executed) -----------------
pre <- grep("^(PASS|FAIL): closure-level", readLines("reproduce_iptw_log.txt"), value = TRUE)
gate(length(pre) == 1L && startsWith(pre, "PASS"), sprintf("[4] before the fix (reproduce_iptw_log.txt): %s", pre))
set.seed(11)
n <- 200L
d <- data.frame(treat = rep(0:1, each = n / 2))
d$y  <- 1 + 0.5 * d$treat + stats::rnorm(n)
ps   <- stats::runif(n, 0.2, 0.8)
d$sw <- ifelse(d$treat == 1L, mean(d$treat) / ps, (1 - mean(d$treat)) / (1 - ps))
est <- make_effect_estimator(outcome_type = "continuous", treat.name = "treat",
                             outcome.name = "y", ps_adjust_method = "iptw")
r <- est(d)
sline(sprintf("INFO: [4] closure after the fix: estimate = %s, se = %s, converged = %s, method_used = %s",
              format(r$estimate), format(r$se), r$converged, r$method_used))
gate(!is.na(r$estimate), "[4] closure-level weighted fit now returns an estimate")

# ---- 5. timing (informational) ------------------------------------------------
tA <- read.csv(file.path(OUT, "stageA_timing.csv")); tC <- tab
pct <- function(a, c) 100 * (c - a) / a
for (fx_name in names(FIXTURES)) {
  a <- tA[tA$fixture == fx_name, ]; c <- tC[tC$fixture == fx_name, ]
  sline(sprintf(paste0("INFO: [5] %s | call median s: %.2f -> %.2f (%+.1f%%) | replicate mean s: %.2f -> %.2f (%+.1f%%) ",
                       "| boot outer s: %.1f -> %.1f (%+.1f%%)"),
                fx_name, a$call_median_s, c$call_median_s, pct(a$call_median_s, c$call_median_s),
                a$rep_mean_s, c$rep_mean_s, pct(a$rep_mean_s, c$rep_mean_s),
                a$boot_outer_s, c$boot_outer_s, pct(a$boot_outer_s, c$boot_outer_s)))
}

# ---- 6. probe (informational) -------------------------------------------------
cv <- forestsearch:::.fs_cox_lean_conventions()
sline(sprintf("INFO: [6] probe: .fs_cox_lean_state$ok = %s; conventions: ties = %s, nocenter = c(%s); survival %s",
              format(forestsearch:::.fs_cox_lean_state$ok), cv$ties, paste(cv$nocenter, collapse = ", "),
              as.character(packageVersion("survival"))))

sline(if (length(failed)) sprintf("STAGE C GATES: %d GATE(S) DID NOT PASS", length(failed)) else "STAGE C GATES: ALL PASS")
if (length(failed)) quit(status = 1)
