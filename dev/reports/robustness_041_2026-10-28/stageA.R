# stageA.R -- timing baseline, end-to-end and bootstrap payloads on master
# (no edits yet).  Task: dev/tasks/cc_task_robustness_041_2026-10-28.md, section 2.
# Transplant of ../lean_fits_2026-10-26/stageA.R minus its deprecation-record gate.
# Run from this directory: Rscript stageA.R
source("timing_common.R")

pd <- packageDescription("forestsearch")[c("Version", "Packaged", "Repository", "Built")]
print(pd)
saveRDS(list(sessionInfo = sessionInfo(), packageDescription = pd),
        file.path(OUT, "stageA_sessionInfo.rds"))
print(sessionInfo())

tab <- run_timing("A")
print(tab)

# INFORMATIONAL (not a gate): Stage A against the committed lean-fits records,
# saved against saved, component by component.
for (fx_name in names(FIXTURES)) {
  a <- readRDS(file.path(OUT, sprintf("stageA_%s_ii.rds", fx_name)))
  r <- readRDS(sprintf("../lean_fits_2026-10-26/identF_%s_ii.rds", fx_name))
  for (k in union(names(r), names(a)))
    cat(sprintf("INFO: %s end to end: %s %s vs identF_%s_ii.rds\n", fx_name, k,
                if (identical(a[[k]], r[[k]])) "identical" else "DIFFERS", fx_name))
  bA <- readRDS(file.path(OUT, sprintf("stageA_%s_boot.rds", fx_name)))
  bR <- readRDS(sprintf("../lean_fits_2026-10-26/stageC_%s_boot.rds", fx_name))
  for (k in union(bR$names_found, bA$names_found))
    cat(sprintf("INFO: %s bootstrap: %s %s vs stageC_%s_boot.rds\n", fx_name, k,
                if (identical(bA[[k]], bR[[k]])) "identical" else "DIFFERS", fx_name))
}
