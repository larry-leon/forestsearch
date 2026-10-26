# stageA.R -- baselines before any edit (installed package = CRAN 0.4.0).
# Task: dev/tasks/cc_task_vi_deprecation_2026-10-26.md, section 2.
# Run from this directory: Rscript stageA.R
source("vi_common.R")

pd <- packageDescription("forestsearch")[c("Version", "Packaged", "Repository", "Built")]
print(pd)
saveRDS(list(sessionInfo = sessionInfo(), packageDescription = pd),
        file.path(OUT, "stageA_sessionInfo.rds"))

tab <- list()
for (fx_name in names(FIXTURES)) {
  fx <- FIXTURES[[fx_name]]()
  for (run in c("i", "ii")) {
    t0 <- proc.time()[["elapsed"]]
    rc <- run_collect(variant_args(fx$args, run))
    x <- extract_run(rc, fx$args$id.name)
    x$elapsed <- proc.time()[["elapsed"]] - t0
    if (fx_name == "survival" && run == "ii") {
      bt <- run_boot(rc$fs)
      saveRDS(extract_boot(bt), file.path(OUT, "stageA_survival_boot.rds"))
      cat("bootstrap components found:", paste(names(bt), collapse = ", "), "\n")
    }
    saveRDS(x, out_file("A", fx_name, run))
    tab[[length(tab) + 1L]] <- summary_row(fx_name, run, x)
    cat(sprintf("[%s %s] n_fits = %d, warnings = %d, %.1fs\n", fx_name, run,
                x$n_fits, x$n_warnings, x$elapsed))
    if (x$n_warnings) print(x$warnings)
  }
}
tab <- do.call(rbind, tab)
print(tab, right = FALSE)
writeLines(capture.output(print(tab, right = FALSE)), file.path(OUT, "stageA_table.txt"))

# GATE: binary (i) exactly one forest fit, every other run zero.
expected <- ifelse(tab$fixture == "binary" & tab$run == "i", 1L, 0L)
if (!identical(as.integer(tab$n_fits), as.integer(expected))) {
  cat("GATE FAIL: forest-fit counts do not match the premise\n")
  quit(status = 1)
}
cat("GATE PASS: binary (i) n_fits = 1, all other runs 0\n")
