# =============================================================================
# Verify the batch_size selection defect is gone
# =============================================================================
# Reproduces the configuration that originally surfaced it (the efficiency
# evaluation's P2 sweep), not the Phase 3/4 acceptance fixture -- those used
# different conf_force and candidate counts, and GBSG collapses most foci to the
# same subgroup under them.
#
# BEFORE the fix, with sg_focus = "hr" (the default) and stop_threshold
# inherited at 0.95:
#     batch_size = 1   ->  {er <= 0} & {size <= 35}   (1 candidate evaluated)
#     batch_size >= 12 ->  {er <= 0} & {pgr <= 33}    (8 candidates evaluated)
# Same data, same statistical settings, different answer.
#
# AFTER the fix, "hr" no longer early-stops, so every candidate is evaluated
# regardless of batch size and the selection must be identical across all three.
#
# Run from the repository root:
#   Rscript dev/sg-focus-work/R/verify_batchsize_fix.R

suppressPackageStartupMessages({
  library(forestsearch)
  library(survival)
})

df <- as.data.frame(survival::gbsg)
df$grade3 <- as.integer(df$grade == 3)
df$id     <- seq_len(nrow(df))

CONFS <- c("age", "meno", "size", "grade3", "nodes", "pgr", "er")

run_one <- function(bs, focus = "hr") {
  fs <- forestsearch(
    df.analysis          = df,
    confounders.name     = CONFS,
    outcome.name         = "rfstime",
    event.name           = "status",
    treat.name           = "hormon",
    id.name              = "id",
    sg_focus             = focus,
    fs.splits            = 400,
    max_subgroups_search = 200,
    maxk                 = 2,
    details              = FALSE,
    quiet                = TRUE,
    parallel_args = list(plan = "sequential", workers = 1L,
                         batch_size = as.integer(bs))
  )
  gc_ <- fs$grp.consistency
  list(
    sg    = if (is.null(fs$sg.harm)) NA_character_
            else paste(sort(as.character(fs$sg.harm)), collapse = " & "),
    n_sel = if (is.null(gc_$sg.harm.id)) NA_integer_ else sum(gc_$sg.harm.id == 1L),
    n_ev  = gc_$n_candidates_evaluated,
    n_tot = gc_$n_candidates_total,
    early = isTRUE(gc_$early_stop_triggered),
    stopt = gc_$stop_threshold
  )
}

cat("=== sg_focus = \"hr\" (default): selection must NOT vary with batch_size ===\n\n")
cat(sprintf("%11s %10s %8s %9s  %s\n",
            "batch_size", "evaluated", "early", "n_sel", "selected subgroup"))

sizes <- c(1L, 8L, 64L)
res <- lapply(sizes, run_one)
for (i in seq_along(sizes)) {
  r <- res[[i]]
  cat(sprintf("%11d %6d/%-3d %8s %9s  %s\n", sizes[i], r$n_ev, r$n_tot,
              r$early, format(r$n_sel), r$sg))
}

sgs <- vapply(res, function(r) r$sg, character(1))
ok  <- length(unique(sgs)) == 1L
cat(sprintf("\nstop_threshold in force: %s  (expected NULL -- reset for \"hr\")\n",
            if (is.null(res[[1]]$stopt)) "NULL" else format(res[[1]]$stopt)))
cat(sprintf("early stopping fired    : %s  (expected FALSE)\n",
            paste(vapply(res, function(r) r$early, logical(1)), collapse = " ")))
cat(sprintf("\n%s\n", if (ok)
  "PASS - selection identical across all batch sizes." else
  "FAIL - selection VARIES with batch_size. The defect is not fixed."))

# ---- maxeffCons: early stopping is active, and must still not change the answer
cat("\n=== sg_focus = \"maxeffCons\": early-stops, yet must still be invariant ===\n\n")
cat(sprintf("%11s %10s %8s %9s  %s\n",
            "batch_size", "evaluated", "early", "n_sel", "selected subgroup"))
res2 <- lapply(sizes, run_one, focus = "maxeffCons")
for (i in seq_along(sizes)) {
  r <- res2[[i]]
  cat(sprintf("%11d %6d/%-3d %8s %9s  %s\n", sizes[i], r$n_ev, r$n_tot,
              r$early, format(r$n_sel), r$sg))
}
sgs2 <- vapply(res2, function(r) r$sg, character(1))
ok2  <- length(unique(sgs2)) == 1L
cat(sprintf("\n%s\n", if (ok2)
  paste("PASS - selection identical although n_candidates_evaluated differs,",
        "\n       which is expected: a larger batch scores more candidates",
        "\n       before the stop condition is noticed.") else
  "FAIL - selection VARIES with batch_size under maxeffCons."))

cat(sprintf("\n%s\n", if (ok && ok2) "BOTH CHECKS PASS" else "*** CHECK FAILED ***"))
