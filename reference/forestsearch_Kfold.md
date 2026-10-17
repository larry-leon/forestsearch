# ForestSearch K-Fold Cross-Validation

This function assesses the stability and reproducibility of ForestSearch
subgroup identification through cross-validation. For each fold:

1.  Train ForestSearch on (K-1) folds

2.  Apply the identified subgroup to the held-out fold

3.  Compare predictions to the original full-data analysis

## Usage

``` r
forestsearch_Kfold(
  fs.est,
  Kfolds = .fs_cv_nobs(fs.est),
  seedit = 8316951L,
  parallel_args = list(plan = "multisession", workers = 6, show_message = TRUE),
  sg0.name = "Not recommend",
  sg1.name = "Recommend",
  details = FALSE,
  mr_in_replicates = FALSE
)
```

## Arguments

- fs.est:

  List. ForestSearch results object from
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md).
  Must contain `df.est` (data frame) and `args_call_all` (list of
  arguments).

- Kfolds:

  Integer. Number of folds (default: `nrow(fs.est$df.est)` for LOO).

- seedit:

  Integer. Random seed for fold assignment (default: 8316951).

- parallel_args:

  List. Parallelization configuration with elements:

  - `plan`: Character. One of "multisession", "multicore", "sequential"

  - `workers`: Integer. Number of parallel workers

  - `show_message`: Logical. Show parallel setup messages

  `parallel_args` controls **fold-level** parallelism only: the fold
  loop is distributed across workers from this argument, while each
  fold's inner
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  refit is pinned to `plan = "sequential"` (deliberate
  nested-parallelism avoidance, so fold-level workers are not
  oversubscribed by inner-search workers). The plan is taken from
  `parallel_args`, not from any ambient
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  set by the caller (which is saved and restored on exit). For example,
  to run the folds on four workers:
  `parallel_args = list(plan = "multisession", workers = 4)`.

- sg0.name:

  Character. Label for subgroup 0 (default: "Not recommend").

- sg1.name:

  Character. Label for subgroup 1 (default: "Recommend").

- details:

  Logical. Print progress details (default: FALSE).

- mr_in_replicates:

  Logical. Whether multiplier resampling (MR) runs *inside* each CV
  fold. **Default `FALSE`**.

  `mr_inference` reaches each fold through the reconstructed call in
  `fs.est$args_call_all`, so an analysis fitted with
  `mr_inference = TRUE` would otherwise run MR in every fold – and
  nothing consumes the result. Under the default the flag is stripped
  from the per-fold call, so that work is not done at all. The top-level
  analysis on `fs.est` is unaffected.

  `TRUE` lets `mr_inference` propagate unchanged into each fold *and*
  retains the result, returned as `$mr_replicates` – a list of length
  `Kfolds` whose `k`-th element is fold `k`'s `mr_inference` object, or
  `NULL` where that fold identified no subgroup. MR therefore runs in
  the folds only if the original analysis set `mr_inference = TRUE`.

  **The retained objects are not aggregable into an estimate for the
  original analysis.** Each fold re-runs the search on its training
  split, selects its own subgroup, and de-biases *that* subgroup – so
  the collection describes MR's behaviour across fold-specific
  selections, not a better estimate of the original quantity. Averaging
  them, or pooling their intervals, is an error.

## Value

List with components:

- resCV:

  Data frame with CV predictions for each observation

- cv_args:

  Arguments used for CV ForestSearch calls

- timing_minutes:

  Execution time in minutes

- prop_SG_found:

  Percentage of folds where a subgroup was found

- sg_analysis:

  Original subgroup definition from full-data analysis

- sg0.name, sg1.name:

  Subgroup labels

- Kfolds:

  Number of folds used

- sens_summary:

  Named vector of sensitivity metrics (sens_H, sens_Hc, ppv_H, ppv_Hc)

- find_summary:

  Named vector of subgroup-finding metrics (Any, Exact, etc.)

- mr_replicates:

  `NULL` under the default `mr_in_replicates = FALSE`. Otherwise a list
  of length `Kfolds`, element `k` holding fold `k`'s `mr_inference`
  object (or `NULL`), each computed against *that fold's* candidate
  family and selected subgroup. Not aggregable – see the
  `mr_in_replicates` parameter.

## Details

Performs K-fold cross-validation for ForestSearch, evaluating subgroup
identification and agreement between training and test sets.

## Cross-Validation Types

- **Leave-One-Out (LOO)**: When `Kfolds = nrow(df)`, each observation is
  held out once. Most thorough but computationally intensive.

- **K-Fold**: When `Kfolds < nrow(df)`, data is split into K roughly
  equal folds. Good balance of bias-variance tradeoff.

## Output Metrics

The returned `resCV` data frame contains:

- `treat.recommend`: Prediction from CV model

- `treat.recommend.original`: Prediction from full-data model

- `cvindex`: Fold assignment

- `sg1`, `sg2`: Subgroup definitions found in each fold

## See also

[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
for initial subgroup identification
[`forestsearch_KfoldOut`](https://larry-leon.github.io/forestsearch/reference/forestsearch_KfoldOut.md)
for summarizing CV results
[`forestsearch_tenfold`](https://larry-leon.github.io/forestsearch/reference/forestsearch_tenfold.md)
for repeated K-fold simulations

## Examples

``` r
if (FALSE) { # \dontrun{
# Run initial ForestSearch
fs_result <- forestsearch(
  df.analysis = trial_data,
  outcome.name = "time",
  event.name = "status",
  treat.name = "treatment",
  confounders.name = c("age", "biomarker")
)

# Run 10-fold cross-validation
cv_results <- forestsearch_Kfold(
  fs.est = fs_result,
  Kfolds = 10,
  parallel_args = list(plan = "multisession", workers = 4),
  details = TRUE
)

# Summarize results
cv_summary <- forestsearch_KfoldOut(cv_results, outall = TRUE)
} # }
```
