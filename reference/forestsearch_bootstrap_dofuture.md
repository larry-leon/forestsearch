# ForestSearch Bootstrap with doFuture Parallelization

Orchestrates bootstrap analysis for ForestSearch using doFuture
parallelization. Implements bias correction methods to adjust for
optimism in subgroup selection.

## Usage

``` r
forestsearch_bootstrap_dofuture(
  fs.est,
  nb_boots,
  details = FALSE,
  show_three = FALSE,
  parallel_args = list()
)
```

## Arguments

- fs.est:

  List. ForestSearch results object from
  [`forestsearch`](https://github.com/larry-leon/forestsearch/reference/forestsearch.md).
  Must contain `df.est` (data frame) and `args_call_all` (list of
  arguments).

- nb_boots:

  Integer. Number of bootstrap samples (recommend 500-1000).

- details:

  Logical. If `TRUE`, prints detailed progress information. Default:
  `FALSE`.

- show_three:

  Logical. If `TRUE`, shows verbose output for first 3 bootstrap
  iterations for debugging. Default: `FALSE`.

- parallel_args:

  List. Parallelization configuration with elements:

  - `plan`: Character. One of "multisession", "multicore", "callr", or
    "sequential"

  - `workers`: Integer. Number of parallel workers

  - `show_message`: Logical. Show parallel setup messages

  If empty list, inherits settings from original forestsearch call.

## Value

List with the following components:

- results:

  Data.table with bias-corrected estimates for each bootstrap iteration

- SG_CIs:

  List of confidence intervals for H and Hc (raw and bias-corrected)

- FSsg_tab:

  Formatted table of subgroup estimates

- Ystar_mat:

  Matrix (nb_boots x n) of bootstrap sample indicators

- H_estimates:

  Detailed estimates for subgroup H

- Hc_estimates:

  Detailed estimates for subgroup Hc

- summary:

  (If create_summary=TRUE) Enhanced summary with tables and diagnostics

## Bias Correction Methods

Two bias correction approaches are implemented:

1.  **Method 1 (Simple Optimism)**: \$\$H\_{adj1} = H\_{obs} -
    (H^\*\_{\*} - H^\*\_{obs})\$\$ where \\H^\*\_{\*}\\ is the new
    subgroup HR on bootstrap data and \\H^\*\_{obs}\\ is the new
    subgroup HR on original data.

2.  **Method 2 (Double Bootstrap)**: \$\$H\_{adj2} = 2 \times H\_{obs} -
    (H\_{\*} + H^\*\_{\*} - H^\*\_{obs})\$\$ where \\H\_{\*}\\ is the
    original subgroup HR on bootstrap data.

## Variable Naming Convention

- `H`: Original subgroup (harm/questionable, treat.recommend == 0)

- `Hc`: Complement subgroup (recommend, treat.recommend == 1)

- `_obs`: Estimate from original data

- `_star`: Estimate from bootstrap data

- `_biasadj_1`: Bias correction method 1

- `_biasadj_2`: Bias correction method 2

## Performance

Typical runtime: 1-5 seconds per bootstrap iteration. For 1000
bootstraps with 6 workers, expect 3-10 minutes total. Memory usage
scales with dataset size and number of workers.

## Requirements

- Original `fs.est` must have identified a valid subgroup

- Requires packages: `data.table`, `foreach`, `doFuture`, `survival`

- For plots: requires `ggplot2`

## See also

[`forestsearch`](https://github.com/larry-leon/forestsearch/reference/forestsearch.md)
for initial subgroup identification
[`bootstrap_results`](https://github.com/larry-leon/forestsearch/reference/bootstrap_results.md)
for the core bootstrap worker function
[`build_cox_formula`](https://github.com/larry-leon/forestsearch/reference/build_cox_formula.md)
for Cox formula construction
[`fit_cox_models`](https://github.com/larry-leon/forestsearch/reference/fit_cox_models.md)
for Cox model fitting

## Examples

``` r
if (FALSE) { # \dontrun{
# Run ForestSearch
fs_result <- forestsearch(
  df.analysis = mydata,
  outcome.name = "time",
  event.name = "status",
  treat.name = "treatment",
  confounders.name = c("age", "sex", "stage")
)

# Run bootstrap with bias correction
boot_results <- forestsearch_bootstrap_dofuture(
  fs.est = fs_result,
  nb_boots = 1000,
  parallel_args = list(
    plan = "multisession",
    workers = 6,
    show_message = TRUE
  ),
  create_summary = TRUE,
  create_plots = TRUE
)

# View results
print(boot_results$FSsg_tab)
print(boot_results$summary$table)

# Check success rate
mean(!is.na(boot_results$results$H_biasadj_2))
} # }
```
