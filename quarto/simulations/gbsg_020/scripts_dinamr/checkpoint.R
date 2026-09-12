# Portability header added when this script was committed: the campaign ran with
# SCRATCH = the session scratchpad and QMD_DIR = quarto/simulations/gbsg_020.
SCRATCH <- Sys.getenv("DINAMR_SCRATCH", unset = dirname(normalizePath(sys.frame(1)$ofile %||% ".", mustWork = FALSE)))
if (!nzchar(SCRATCH) || is.na(SCRATCH)) SCRATCH <- "."
QMD_DIR <- Sys.getenv("DINAMR_QMD_DIR", unset = "..")
`%||%` <- function(a, b) if (is.null(a)) b else a

# ===== AMENDMENT 2 CHECKPOINT, run ONCE after Block A ========================
# Re-project Block B from Block A's REALIZED walls against the remaining budget
# to the 13 h ceiling (Amendment 4).  Block C is PINNED DEFERRED regardless of
# what this shows (Larry, checkpoint pin) -- this decides Block B's cells only.
CEILING_H <- 13; TIMEOUT_H <- 16
log <- "/private/tmp/claude-501/-Users-larryleon-Documents-GitHub-forestsearch/c9c94bd3-80ce-4b91-9acd-4d47ca95cfaa/tasks/bft741io9.output"
L <- readLines(log)
done <- grep("^CELL DONE", L, value = TRUE)
realized <- data.frame(
  cell = sub("^CELL DONE: ([A-Za-z0-9_]+).*", "\\1", done),
  wall_s = as.numeric(sub(".*wall=([0-9]+)s.*", "\\1", done)), stringsAsFactors = FALSE)
# Gate 1's per-cell projections (s), from the measured probe surface.
proj <- c(A_h150_n500 = 1484.3, A_h150_n1000 = 1345.9, A_h150_n1500 = 991.6,
          A_h175_n500 = 1484.3, A_h175_n1000 = 1345.9, A_h175_n1500 = 991.6,
          B_h150_n500 = 4204.9, B_h150_n1000 = 5325.4, B_h150_n1500 = 5822.5,
          B_h175_n500 = 4204.9, B_h175_n1000 = 5325.4, B_h175_n1500 = 5822.5)
realized$projected_s <- unname(proj[realized$cell])
realized$ratio <- realized$wall_s / realized$projected_s
cat("=== BLOCK A: REALIZED vs GATE 1 PROJECTION ===\n")
print(transform(realized, wall_min = round(wall_s/60,1), proj_min = round(projected_s/60,1),
                ratio = round(ratio,3))[, c("cell","wall_s","projected_s","wall_min","proj_min","ratio")],
      row.names = FALSE)
A_real <- sum(realized$wall_s); A_proj <- sum(realized$projected_s)
K <- A_real / A_proj
cat(sprintf("\nBlock A realized total : %.0f s = %.3f h  (%d of 6 cells)\n", A_real, A_real/3600, nrow(realized)))
cat(sprintf("Block A projected total: %.0f s = %.3f h\n", A_proj, A_proj/3600))
cat(sprintf("CALIBRATION FACTOR K   : %.4f  (realized / projected; per-cell range %.3f-%.3f)\n",
            K, min(realized$ratio), max(realized$ratio)))
# Per-n calibration, since the n-profile is what drives Block B.
realized$n <- as.integer(sub(".*_n", "", realized$cell))
byn <- aggregate(cbind(wall_s, projected_s) ~ n, realized, sum)
byn$K_n <- byn$wall_s / byn$projected_s
cat("\nCalibration by n (Block B's cost is n-driven):\n"); print(byn, row.names = FALSE, digits = 4)
Bcells <- names(proj)[grepl("^B_", names(proj))]
BB <- data.frame(cell = Bcells, projected_s = unname(proj[Bcells]), stringsAsFactors = FALSE)
BB$n <- as.integer(sub(".*_n", "", BB$cell))
BB$K_used <- byn$K_n[match(BB$n, byn$n)]
BB$K_used[is.na(BB$K_used)] <- K
BB$reproj_s <- BB$projected_s * BB$K_used
BB$reproj_h <- BB$reproj_s / 3600
cat("\n=== BLOCK B RE-PROJECTION (Block A's realized walls, calibrated by n) ===\n")
print(transform(BB, reproj_min = round(reproj_s/60,1))[, c("cell","n","projected_s","K_used","reproj_s","reproj_min","reproj_h")],
      row.names = FALSE, digits = 4)
budget_h <- CEILING_H - A_real/3600
cat(sprintf("\nCeiling %.0f h  -  Block A realized %.3f h  =  REMAINING BUDGET %.3f h\n", CEILING_H, A_real/3600, budget_h))
cat(sprintf("Block B re-projected total: %.3f h\n", sum(BB$reproj_h)))
cat(sprintf("A + B re-projected total  : %.3f h against the %.0f h ceiling\n", A_real/3600 + sum(BB$reproj_h), CEILING_H))
fits <- sum(BB$reproj_h) <= budget_h
cat(sprintf("\nDOES BLOCK B FIT IN FULL?  %s  (headroom %.3f h = %.1f%%)\n",
            if (fits) "YES" else "NO", budget_h - sum(BB$reproj_h),
            100*(budget_h - sum(BB$reproj_h))/budget_h))
if (!fits) {
  # Stated defer order within B: n = 1500 cells first, then n = 1000.
  ord <- c("B_h150_n1500","B_h175_n1500","B_h150_n1000","B_h175_n1000")
  keep <- BB; dropped <- character(0)
  for (cl in ord) { if (sum(keep$reproj_h) <= budget_h) break
    dropped <- c(dropped, cl); keep <- keep[keep$cell != cl, ] }
  cat("\nDEFER PER THE STATED ORDER (B's n = 1500 first, then n = 1000):\n")
  cat("  DEFER: ", if (length(dropped)) paste(dropped, collapse=", ") else "none", "\n")
  cat("  RUN  : ", paste(keep$cell, collapse=", "), sprintf("  (%.3f h)\n", sum(keep$reproj_h)))
  writeLines(keep$cell, file.path(SCRATCH, "blockB.run")
} else {
  cat("\nRUN ALL SIX BLOCK B CELLS.\n")
  writeLines(BB$cell, file.path(SCRATCH, "blockB.run")
}
cat(sprintf("\nBLOCK C: PINNED DEFERRED to a follow-up session, whatever this shows (Larry, checkpoint pin).\n"))
cat(sprintf("Hard timeout stands at %.0f h.\n", TIMEOUT_H))
