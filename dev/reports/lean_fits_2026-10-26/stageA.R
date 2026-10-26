# stageA.R -- timing baseline and bootstrap payloads on the deprecation build
# (no edits yet).  Task: dev/tasks/cc_task_lean_fits_2026-10-26.md, section 2.
# Run from this directory: Rscript stageA.R
source("timing_common.R")

pd <- packageDescription("forestsearch")[c("Version", "Packaged", "Repository", "Built")]
print(pd)
saveRDS(list(sessionInfo = sessionInfo(), packageDescription = pd),
        file.path(OUT, "stageA_sessionInfo.rds"))

tab <- run_timing("A")
print(tab)

# GATE: the survival bootstrap payload reproduces the committed deprecation record
# (saved against saved).
bA <- readRDS(file.path(OUT, "stageA_survival_boot.rds"))
bR <- readRDS("../vi_deprecation_2026-10-26/stageC_survival_boot.rds")
ok <- identical(bA$names_found, bR$names_found)
cat(sprintf("%s: survival bootstrap: same components found (%s)\n",
            if (ok) "PASS" else "FAIL", paste(bA$names_found, collapse = ", ")))
for (k in bR$names_found) {
  r <- identical(bA[[k]], bR[[k]])
  cat(sprintf("%s: survival bootstrap: %s identical to committed stageC_survival_boot.rds\n",
              if (r) "PASS" else "FAIL", k))
  ok <- ok && r
}
if (!ok) { cat("GATE FAIL: survival bootstrap payload\n"); quit(status = 1) }
cat("GATE PASS: survival bootstrap payload identical to the committed record\n")
