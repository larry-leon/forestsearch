# stageC.R -- verification against the installed branch build.
# Task: dev/tasks/cc_task_vi_deprecation_2026-10-26.md, section 4.
# Run from this directory: Rscript stageC.R   (identity runs + gates only;
# the test files are run separately under their own timeouts).
source("vi_common.R")

pd <- packageDescription("forestsearch")[c("Version", "Packaged", "Repository", "Built")]
print(pd)
saveRDS(list(sessionInfo = sessionInfo(), packageDescription = pd),
        file.path(OUT, "stageC_sessionInfo.rds"))

SUMMARY <- file.path(OUT, "stageC_summary.txt")
cat(sprintf("Stage C summary -- installed forestsearch %s (Built %s)\n", pd$Version, pd$Built),
    file = SUMMARY)
sline <- function(...) { s <- paste0(...); cat(s, "\n"); cat(s, "\n", file = SUMMARY, append = TRUE, sep = "") }
failed <- character(0)
gate <- function(ok, label) {
  sline(if (isTRUE(ok)) "PASS" else "FAIL", ": ", label)
  if (!isTRUE(ok)) failed <<- c(failed, label)
}
n_depr <- function(x) sum(grepl("deprecated", x$warnings))
cmp_components <- function(a, b, tag) {
  for (k in IDENTITY_COMPONENTS) gate(identical(a[[k]], b[[k]]), sprintf("%s: %s identical", tag, k))
}

A <- list(); C <- list(); tab <- list()
for (fx_name in names(FIXTURES)) {
  fx <- FIXTURES[[fx_name]]()
  for (run in c("i", "ii")) {
    t0 <- proc.time()[["elapsed"]]
    rc <- run_collect(variant_args(fx$args, run))
    x <- extract_run(rc, fx$args$id.name)
    x$elapsed <- proc.time()[["elapsed"]] - t0
    if (fx_name == "survival" && run == "ii") {
      bt <- run_boot(rc$fs)
      saveRDS(extract_boot(bt), file.path(OUT, "stageC_survival_boot.rds"))
    }
    saveRDS(x, out_file("C", fx_name, run))
    key <- paste(fx_name, run)
    # Compare saved against saved: a data.table read back by readRDS() has a
    # nil .internal.selfref pointer, so identical() on an in-memory table vs a
    # reloaded one is FALSE even when the contents are identical (run 1 of
    # this script hit exactly that; see stageC_summary_run1_selfref_artifact.txt).
    C[[key]] <- readRDS(out_file("C", fx_name, run))
    A[[key]] <- readRDS(out_file("A", fx_name, run))
    tab[[length(tab) + 1L]] <- summary_row(fx_name, run, x)
    cat(sprintf("[%s %s] n_fits = %d, warnings = %d, %.1fs\n", fx_name, run,
                x$n_fits, x$n_warnings, x$elapsed))
    if (x$n_warnings) print(x$warnings)
  }
}
tab <- do.call(rbind, tab)
print(tab, right = FALSE)
writeLines(capture.output(print(tab, right = FALSE)), file.path(OUT, "stageC_table.txt"))

# ---- survival and continuous, both runs; binary (ii) ---------------------------
for (key in c("survival i", "survival ii", "continuous i", "continuous ii", "binary ii")) {
  cmp_components(C[[key]], A[[key]], sprintf("%s vs Stage A", key))
  gate(n_depr(C[[key]]) == 0L, sprintf("%s: zero deprecation warnings", key))
  if (key == "binary ii") gate(C[[key]]$n_warnings == 0L, "binary ii: zero warnings")
  gate(C[[key]]$n_fits == 0L, sprintf("%s: zero forest fits", key))
}

# ---- binary (i): vi.grf.min = -0.2 ---------------------------------------------
bi <- C[["binary i"]]
gate(n_depr(bi) == 1L, "binary i: exactly one warning matching 'deprecated'")
gate(bi$n_fits == 0L, sprintf("binary i: zero forest fits (Stage A had %d)", A[["binary i"]]$n_fits))
cmp_components(bi, C[["binary ii"]], "binary i vs Stage C binary ii")
cmp_components(bi, A[["binary ii"]], "binary i vs Stage A binary ii")
sline("FINDING: binary i warnings (", bi$n_warnings, "): ", paste(bi$warnings, collapse = " | "))

# Report, not gate: against the old screened run, Stage A (i).
ai <- A[["binary i"]]
sline("FINDING: binary i vs Stage A (i) [report only]: setequal(sg.harm clauses) = ",
      setequal(bi$sg.harm, ai$sg.harm))
sline("FINDING: binary i vs Stage A (i) [report only]: identical(sorted membership) = ",
      identical(bi$membership, ai$membership))
cand_norm <- function(h) {
  h <- as.data.frame(h)
  qcols <- grep("^q[0-9]+\\.[01]$", names(h), value = TRUE)
  lab <- apply(as.matrix(h[, qcols, drop = FALSE]), 1L, function(r) paste(sort(qcols[r == 1]), collapse = " & "))
  keep <- setdiff(names(h), c("grp", qcols))
  out <- data.frame(h[, keep, drop = FALSE], label = lab, check.names = FALSE, stringsAsFactors = FALSE)
  out <- out[order(out$HR, out$K, out$label), , drop = FALSE]
  rownames(out) <- NULL
  out
}
cn_b <- cand_norm(bi$hr.subgroups); cn_a <- cand_norm(ai$hr.subgroups)
sline("FINDING: binary i vs Stage A (i) [report only]: candidate table equal after sorting by (HR, K, label) = ",
      isTRUE(all.equal(cn_b, cn_a, check.attributes = FALSE)),
      sprintf(" (identical = %s; rows %d vs %d; label = sorted active q-columns, grp dropped)",
              identical(cn_b, cn_a), nrow(cn_b), nrow(cn_a)))
sline("FINDING: binary i vs Stage A (i) [report only]: sg.harm clause order differs = ",
      !identical(bi$sg.harm, ai$sg.harm), sprintf(" (Stage C: %s; Stage A (i): %s)",
      paste(bi$sg.harm, collapse = " & "), paste(ai$sg.harm, collapse = " & ")))
for (k in setdiff(IDENTITY_COMPONENTS, c("sg.harm", "membership", "hr.subgroups")))
  sline("FINDING: binary i vs Stage A (i) [report only]: ", k, " identical = ", identical(bi[[k]], ai[[k]]))

# ---- survival nb = 20 bootstrap payload -----------------------------------------
bA <- readRDS("stageA_survival_boot.rds"); bC <- readRDS(file.path(OUT, "stageC_survival_boot.rds"))
gate(identical(bA$names_found, bC$names_found),
     sprintf("survival bootstrap: same components found (%s)", paste(bC$names_found, collapse = ", ")))
for (k in bA$names_found) gate(identical(bC[[k]], bA[[k]]), sprintf("survival bootstrap: %s identical", k))

sline(if (length(failed)) sprintf("STAGE C IDENTITY GATES: %d FAILED", length(failed))
      else "STAGE C IDENTITY GATES: ALL PASS")
if (length(failed)) quit(status = 1)
