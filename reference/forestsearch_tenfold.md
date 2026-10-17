# ForestSearch Repeated K-Fold Cross-Validation

This function performs multiple independent K-fold cross-validations to
assess the variability in subgroup identification. Each simulation:

1.  Randomly shuffles the data

2.  Performs K-fold CV

3.  Records sensitivity and agreement metrics

Results are summarized across all simulations.

## Usage

``` r
forestsearch_tenfold(
  fs.est,
  sims,
  Kfolds = 10,
  details = TRUE,
  seed = 8316951L,
  parallel_args = list(plan = "multisession", workers = 6, show_message = TRUE),
  keep_resCV = FALSE,
  mr_in_replicates = FALSE
)
```

## Arguments

- fs.est:

  List. ForestSearch results object from
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md).

- sims:

  Integer. Number of simulation repetitions.

- Kfolds:

  Integer. Number of folds per simulation (default: 10).

- details:

  Logical. Print progress details (default: TRUE).

- seed:

  Integer. Base random seed for fold shuffling. Default 8316951L. Each
  simulation uses seed + 1000 \* ksim for reproducibility.

- parallel_args:

  List. Parallelization configuration (elements `plan` / `workers` /
  `show_message`; see **Parallelization Strategy**). Here
  `parallel_args` controls **simulation-level** parallelism: the `sims`
  repetitions are distributed across workers, while each simulation runs
  its K-fold CV – and every inner
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  refit – sequentially (so workers are not oversubscribed by nested
  parallelism). As in
  [`forestsearch_Kfold`](https://larry-leon.github.io/forestsearch/reference/forestsearch_Kfold.md),
  the plan is taken from `parallel_args`, not from any ambient
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  set by the caller (which is saved and restored on exit). For example,
  to run the simulations on four workers:
  `parallel_args = list(plan = "multisession", workers = 4)`.

- keep_resCV:

  Logical. If `TRUE`, preserve the full per-subject cross-validation
  prediction data frame for every simulation in the returned `resCV_all`
  slot. Useful for post-hoc metric computation, reproducibility
  debugging, or auditing which subjects landed in which fold. Memory
  cost scales as `sims * nrow(df.est) * ~9 columns` (for a 686-subject
  trial with 200 sims, ~10 MB; for a 10,000-subject trial with 1,000
  sims, ~1 GB). Default: `FALSE`. Regardless of this flag, a compact
  `fold_summary` (one row per sim x fold) is always returned; see
  Return.

- mr_in_replicates:

  Logical. Whether multiplier resampling (MR) runs *inside* each CV
  fold. **Default `FALSE`**.

  `mr_inference` reaches each fold through the reconstructed call in
  `fs.est$args_call_all`, so an analysis fitted with
  `mr_inference = TRUE` would otherwise run MR in every fold of every
  simulation – `sims * Kfolds` runs – and nothing consumes the result.
  Under the default the flag is stripped from the per-fold call, so that
  work is not done at all. The top-level analysis on `fs.est` is
  unaffected.

  `TRUE` lets `mr_inference` propagate unchanged into each fold *and*
  retains the results, returned as `$mr_replicates` – a list with one
  element per successful simulation (named `sim<k>` after its `sim_id`,
  so a dropped simulation shows as a gap rather than a silent re-index),
  each itself a list of length `Kfolds` holding that fold's
  `mr_inference` object or `NULL`. MR therefore runs in the folds only
  if the original analysis set `mr_inference = TRUE`.

  **The retained objects are not aggregable into an estimate for the
  original analysis.** Each fold re-runs the search on its training
  split, selects its own subgroup, and de-biases *that* subgroup – so
  the collection describes MR's behaviour across fold-specific
  selections, not a better estimate of the original quantity. Averaging
  them, or pooling their intervals, is an error.

## Value

List with components:

- sens_summary:

  Named vector of median sensitivity metrics across simulations

- find_summary:

  Named vector of median subgroup-finding metrics

- sens_out:

  Matrix of sensitivity metrics (sims x metrics)

- find_out:

  Matrix of finding metrics (sims x metrics)

- fold_summary:

  Data frame with one row per (sim, fold) combination. Columns: `sim`,
  `fold`, `n_test`, `sg1`, `sg2`, `grf_cuts`, `pconsistency`,
  `training_fs_hr`, `n_candidates_evaluated`, `any_found`. The
  `grf_cuts` column records the GRF policy-tree cut expressions returned
  for each training fold (collapsed with " \| " when multiple cuts are
  returned; `NA` when GRF was not used, failed, or returned no cuts).
  The `pconsistency` column records the consistency probability
  (`Pcons`) achieved by the identified subgroup on each fold (`NA_real_`
  when no subgroup was identified or the training ForestSearch call
  errored). The `training_fs_hr` column records the in-sample effect
  estimate for the identified subgroup on the training fold, on its
  NATURAL SCALE regardless of outcome type (hazard ratio for survival;
  odds ratio, rate ratio, or risk ratio for GLM ratio measures after
  internal exponentiation; mean/risk difference for GLM additive
  measures); optimistically biased relative to any independent estimate
  and surfaced for diagnostic comparison only, with `NA_real_` when no
  subgroup was identified. The `n_candidates_evaluated` column records
  the integer count of candidate subgroups actually evaluated for
  consistency on this training fold (populated whenever the consistency
  stage ran, even if zero candidates met the threshold and thus no
  subgroup was identified; `NA_integer_` when the training call errored
  or consistency evaluation did not run). Always returned (compact;
  cheap). Lets you tabulate, for example, which subgroup was identified
  in each fold of each simulation, the empirical distribution of GRF cut
  choices across the full sim x fold grid, the relationship between
  GRF's cut and the final identified subgroup, or the near-miss
  consistency values among folds that did not surface a subgroup.

- resCV_all:

  List of length `sims`; element *i* is the per-subject `resCV` data
  frame from simulation *i*. Only populated when `keep_resCV = TRUE`;
  otherwise `NULL`.

- timing_minutes:

  Total execution time

- sims:

  Number of simulations run

- Kfolds:

  Number of folds per simulation

- mr_replicates:

  `NULL` under the default `mr_in_replicates = FALSE`. Otherwise a
  nested list: one element per successful simulation (named `sim<k>`),
  each a list of length `Kfolds` holding that fold's `mr_inference`
  object (or `NULL`), computed against *that fold's* candidate family
  and selected subgroup. Not aggregable – see the `mr_in_replicates`
  parameter.

## Details

Runs repeated K-fold cross-validation simulations for ForestSearch and
summarizes subgroup identification stability across repetitions.

## Parallelization Strategy

Unlike the single K-fold function which parallelizes across folds, this
function parallelizes across simulations for better efficiency when
running many repetitions. Each simulation runs its K-fold CV
sequentially.

## See also

[`forestsearch_Kfold`](https://larry-leon.github.io/forestsearch/reference/forestsearch_Kfold.md)
for single K-fold CV
[`forestsearch_KfoldOut`](https://larry-leon.github.io/forestsearch/reference/forestsearch_KfoldOut.md)
for summarizing CV results

## Examples

``` r
if (FALSE) { # \dontrun{
# Run 100 repetitions of 10-fold CV
tenfold_results <- forestsearch_tenfold(
  fs.est = fs_result,
  sims = 100,
  Kfolds = 10,
  parallel_args = list(plan = "multisession", workers = 6),
  details = TRUE
)

# View summary
print(tenfold_results$sens_summary)
print(tenfold_results$find_summary)
} # }
```
