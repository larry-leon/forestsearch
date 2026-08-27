# Search + consistency reproducibility across worker counts.
#
# The consistency stage pre-generates one L'Ecuyer-CMRG stream per candidate,
# indexed by GLOBAL candidate position (.make_candidate_rng_seeds), and hands
# future_lapply the streams for each batch (candidate_seeds[batch_indices]).
# That makes every candidate's resampling RNG a function of (master seed,
# global index) only -- invariant to worker count and batch size.  This test
# pins the invariant at workers 1 vs 2.
#
# Sensitivity: batch size is min(2 * workers, n_candidates), so workers 1 and
# 2 partition the candidates differently (batches of 2 vs 4).  If the stream
# indexing regresses to WITHIN-BATCH position (the historical defect: a scalar
# future.seed per batch call), candidate i's draws change with the batch
# partition and the per-candidate Pcons values diverge between worker counts.
# The fixture is sized so many candidates carry mid-range Pcons -- a
# first-candidate-only comparison would miss the defect, because the first
# element of the first batch receives stream 1 under either indexing.
#
# consistency_method = "split" is used deliberately: the production default
# "resample" evaluates the closed-form multiplier rate and consumes NO RNG,
# so it is worker-invariant structurally and cannot detect a stream-indexing
# regression.  "split" performs literal random sample splits -- the path the
# per-candidate streams exist to govern, and the production fallback whenever
# the closed form is inestimable on a candidate.

.repro_df <- function(n = 360, seed = 42) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n); x3 <- rnorm(n)
  b1 <- rbinom(n, 1, 0.5); b2 <- rbinom(n, 1, 0.4)
  treat <- rbinom(n, 1, 0.5)
  harm <- as.integer(x1 > 0 & b1 == 1)
  y <- 50 + 5 * x2 + treat * (2 - 8 * harm) + rnorm(n, sd = 4)
  data.frame(id = seq_len(n), y = y, treat = treat,
             x1 = x1, x2 = x2, x3 = x3, b1 = b1, b2 = b2)
}

.repro_run <- function(workers) {
  fs <- suppressWarnings(forestsearch(
    df.analysis = .repro_df(),
    confounders.name = c("x1", "x2", "x3", "b1", "b2"),
    outcome.name = "y", treat.name = "treat", id.name = "id",
    outcome_type = "continuous", effect_measure = "MD",
    effect.threshold = 2, consistency.threshold = 2,
    pconsistency.threshold = 0.05, fs.splits = 40L,
    n.min = 30, d0.min = 5, d1.min = 5, maxk = 2L,
    sg_focus = "maxeffCons", consistency_method = "split",
    stop_threshold = NULL,   # evaluate every candidate; no early stop
    is.RCT = TRUE, adverse_outcome = FALSE,
    details = FALSE, quiet = TRUE, seedit = 20260827,
    parallel_args = list(plan = "multisession", workers = workers,
                         show_message = FALSE),
    mr_inference = FALSE))
  future::plan(future::sequential)
  fs
}

# Strip data.table attributes (.internal.selfref is an externalptr, which
# identical() compares by address) so identical() tests values bitwise.
.repro_flat <- function(x) {
  x <- as.data.frame(x)
  rownames(x) <- NULL
  x
}

test_that("search + consistency results are invariant to worker count", {
  skip_on_cran()   # two full forestsearch runs (~10 s) with multisession pools

  fs1 <- .repro_run(1L)
  fs2 <- .repro_run(2L)

  # The fixture must actually exercise the multi-batch regime, or the test
  # cannot see a stream-indexing regression: candidates must outnumber the
  # workers-2 batch size (4), and several mid-range Pcons rows must exist.
  n_cand <- fs1$grp.consistency$n_candidates_total
  expect_gt(n_cand, 4L)
  pc1 <- fs1$grp.consistency$out_sg$result
  expect_gt(nrow(pc1), 3L)
  expect_true(any(as.numeric(pc1$Pcons) > 0.1 & as.numeric(pc1$Pcons) < 0.9))

  # Search candidate table: deterministic fits, no RNG -- must match exactly.
  expect_identical(.repro_flat(fs1$find.grps$out.found$hr.subgroups),
                   .repro_flat(fs2$find.grps$out.found$hr.subgroups))

  # Per-candidate consistency table, every Pcons included: this is the
  # RNG-dependent surface the per-candidate streams make worker-invariant.
  expect_identical(.repro_flat(pc1),
                   .repro_flat(fs2$grp.consistency$out_sg$result))

  # The selected subgroup string.
  expect_identical(fs1$sg.harm, fs2$sg.harm)
})
