# =============================================================================
# Configuration for the forestsearch efficiency evaluation
# =============================================================================
# Edit these; nothing else in the harness should need changing.

EVAL <- list(

  ## ---- what to run -------------------------------------------------------
  run_glm        = FALSE,     # TRUE adds the ACTG175 binary/OR gates + micro
  quick          = FALSE,     # TRUE shrinks reps/sweeps for a smoke run

  ## ---- worker counts to sweep (P stages) ---------------------------------
  ## RESOLVED ON-MACHINE: this box has 128 physical cores, so the auto-derived
  ## sweep topped out at W = 128 = detectCores(). A multisession (PSOCK) pool of
  ## 128 CANNOT be created: R's default connection limit is 128 and ~4 are
  ## already in use, so 128 nodes fail ("only 124 connections left"). The
  ## package's own default, .default_parallel_workers() = floor(0.80 * 128) =
  ## 102, stays under that ceiling and is the realistic W for a real analysis,
  ## so the sweep is pinned to feasible values with 102 as the top.
  workers_sweep  = c(1L, 8L, 32L, 64L, 102L),
  workers_main   = 102L,      # = .default_parallel_workers() on this machine

  ## ---- batch_size sweep (P2, the headline) -------------------------------
  ## NULL -> c(1, W/8, W/2, W, 2W) deduplicated
  batch_sweep    = NULL,

  ## ---- forestsearch() call used for E / P2 -------------------------------
  ## Kept deliberately close to a real analysis. Adjust to match the
  ## configuration whose runtime you actually care about.
  ## NOTE (resolved on-machine): the first forestsearch() formal is
  ## `df.analysis`, not `df`; the consistency-split count is `fs.splits`
  ## (not `n.splits`) and the subgroup cap is `max_subgroups_search` (not
  ## `stop_Kgroups`). n.splits/stop_Kgroups are subgroup.consistency()'s
  ## internal names; forestsearch() has no `...`, so passing them errors.
  fs_args = list(
    confounders.name     = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
    outcome.name         = "rfstime",
    treat.name           = "hormon",
    event.name           = "status",
    fs.splits            = 400,
    max_subgroups_search = 200,
    maxk                 = 2,
    details              = FALSE
  ),

  ## ---- microbenchmark settings -------------------------------------------
  micro_times    = 50L,
  gate_reps      = 40L,       # random subgroups per equivalence gate
  n_splits_gate  = 200L,

  ## ---- output -------------------------------------------------------------
  results_dir    = "results",
  seed           = 8316951
)

if (isTRUE(EVAL$quick)) {
  EVAL$micro_times   <- 10L
  EVAL$gate_reps     <- 10L
  EVAL$n_splits_gate <- 60L
  EVAL$fs_args$n.splits <- 100
}
