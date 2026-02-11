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
  Kfolds = nrow(fs.est$df.est),
  seedit = 8316951L,
  parallel_args = list(plan = "multisession", workers = 6, show_message = TRUE),
  sg0.name = "Not recommend",
  sg1.name = "Recommend",
  details = FALSE
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

- sg0.name:

  Character. Label for subgroup 0 (default: "Not recommend").

- sg1.name:

  Character. Label for subgroup 1 (default: "Recommend").

- details:

  Logical. Print progress details (default: FALSE).

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
