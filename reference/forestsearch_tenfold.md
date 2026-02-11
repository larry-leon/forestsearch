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
  parallel_args = list(plan = "multisession", workers = 6, show_message = TRUE)
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

  List. Parallelization configuration.

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

- timing_minutes:

  Total execution time

- sims:

  Number of simulations run

- Kfolds:

  Number of folds per simulation

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
