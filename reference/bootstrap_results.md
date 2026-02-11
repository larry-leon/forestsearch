# Bootstrap Results for ForestSearch with Bias Correction

Runs bootstrap analysis for ForestSearch, fitting Cox models and
computing bias-corrected estimates and valid CIs (see vignette for
references)

## Usage

``` r
bootstrap_results(
  fs.est,
  df_boot_analysis,
  cox.formula.boot,
  nb_boots,
  show_three,
  H_obs,
  Hc_obs,
  seed = 8316951L
)
```

## Arguments

- fs.est:

  List. ForestSearch results object from
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md).
  Must contain:

  - `df.est`: Data frame with analysis data including `treat.recommend`

  - `confounders.candidate`: Character vector of confounder names

  - `args_call_all`: List of original forestsearch call arguments

- df_boot_analysis:

  Data frame. Bootstrap analysis data with same structure as
  `fs.est$df.est`. Must contain columns for outcome, event, treatment,
  and the `treat.recommend` flag.

- cox.formula.boot:

  Formula. Cox model formula for bootstrap, typically created by
  [`build_cox_formula`](https://larry-leon.github.io/forestsearch/reference/build_cox_formula.md).
  Should be of form `Surv(outcome, event) ~ treatment`.

- nb_boots:

  Integer. Number of bootstrap samples to generate (e.g., 500-1000).
  More iterations provide better bias correction but increase
  computation time.

- show_three:

  Logical. If `TRUE`, prints detailed progress for the first three
  bootstrap iterations for debugging purposes. Default: `FALSE`.

- H_obs:

  Numeric. Observed log hazard ratio for subgroup H (harm/questionable
  group, `treat.recommend == 0`) from original sample. Used as reference
  for bias correction.

- Hc_obs:

  Numeric. Observed log hazard ratio for subgroup H^c
  (complement/recommend, `treat.recommend == 1`) from original sample.
  Used as reference for bias correction.

- seed:

  Integer. Random seed for reproducibility. Default 8316951L. Must match
  the seed used in
  [`bootstrap_ystar`](https://larry-leon.github.io/forestsearch/reference/bootstrap_ystar.md)
  to ensure bootstrap index alignment.

## Value

Data.table with one row per bootstrap iteration and columns:

- boot_id:

  Integer. Bootstrap iteration number (1 to `nb_boots`)

- H_biasadj_1:

  Bias-corrected estimate for H using method 1:
  `H_obs - (Hstar_star - Hstar_obs)`

- H_biasadj_2:

  Bias-corrected estimate for H using method 2:
  `2*H_obs - (H_star + Hstar_star - Hstar_obs)`

- Hc_biasadj_1:

  Bias-corrected estimate for H^c using method 1

- Hc_biasadj_2:

  Bias-corrected estimate for H^c using method 2

- max_sg_est:

  Numeric. Maximum subgroup hazard ratio found

- L:

  Integer. Number of candidate factors evaluated

- max_count:

  Integer. Maximum number of factor combinations

- events_H_0:

  Integer. Number of events in control arm of original subgroup H on
  bootstrap sample

- events_H_1:

  Integer. Number of events in treatment arm of original subgroup H on
  bootstrap sample

- events_Hc_0:

  Integer. Number of events in control arm of original subgroup H^c on
  bootstrap sample

- events_Hc_1:

  Integer. Number of events in treatment arm of original subgroup H^c on
  bootstrap sample

- events_Hstar_0:

  Integer. Number of events in control arm of new subgroup H\* on
  original data

- events_Hstar_1:

  Integer. Number of events in treatment arm of new subgroup H\* on
  original data

- events_Hcstar_0:

  Integer. Number of events in control arm of new subgroup H^c\* on
  original data

- events_Hcstar_1:

  Integer. Number of events in treatment arm of new subgroup H^c\* on
  original data

- tmins_search:

  Numeric. Minutes spent on subgroup search in this iteration

- tmins_iteration:

  Numeric. Total minutes for this bootstrap iteration

- Pcons:

  Numeric. Consistency p-value for top subgroup

- hr_sg:

  Numeric. Hazard ratio for top subgroup

- N_sg:

  Integer. Sample size of top subgroup

- E_sg:

  Integer. Number of events in top subgroup

- K_sg:

  Integer. Number of factors defining top subgroup

- g_sg:

  Numeric. Subgroup group ID

- m_sg:

  Numeric. Subgroup index

- M.1:

  Character. First factor label

- M.2:

  Character. Second factor label

- M.3:

  Character. Third factor label

- M.4:

  Character. Fourth factor label

- M.5:

  Character. Fifth factor label

- M.6:

  Character. Sixth factor label

- M.7:

  Character. Seventh factor label

Rows where no valid subgroup was found will have `NA` for bias
corrections. The returned object has a "timing" attribute with summary
statistics.

## Note

This function is designed to be called within a `foreach` loop with
`%dofuture%` operator. It requires:

- All functions in
  [`get_bootstrap_exports`](https://larry-leon.github.io/forestsearch/reference/get_bootstrap_exports.md)
  to be available in the parallel workers

- Packages listed in `BOOTSTRAP_REQUIRED_PACKAGES` to be installed

- Proper parallel backend setup via
  [`setup_parallel_SGcons`](https://larry-leon.github.io/forestsearch/reference/setup_parallel_SGcons.md)

## Bias Correction Methods

Two bias correction approaches are implemented:

1.  **Method 1 (Simple Optimism)**: \$\$H\_{adj1} = H\_{obs} -
    (H^\*\_{\*} - H^\*\_{obs})\$\$ where \\H^\*\_{\*}\\ is the new
    subgroup HR on bootstrap data and \\H^\*\_{obs}\\ is the new
    subgroup HR on original data.

2.  **Method 2 (Double Bootstrap)**: \$\$H\_{adj2} = 2 \times H\_{obs} -
    (H\_{\*} + H^\*\_{\*} - H^\*\_{obs})\$\$ where \\H\_{\*}\\ is the
    original subgroup HR on bootstrap data.

where:

- `H_obs`: Original subgroup HR on original data

- `H_star`: Original subgroup HR on bootstrap data

- `Hstar_obs`: New subgroup (found in bootstrap) HR on original data

- `Hstar_star`: New subgroup (found in bootstrap) HR on bootstrap data

## Computational Details

- Uses `doFuture` backend for parallel execution (configured externally)

- Sets reproducible seeds: `8316951 + boot * 100` for each iteration

- Each bootstrap iteration runs full ForestSearch pipeline including
  variable selection, subgroup search, and consistency evaluation

- Sequential execution within each bootstrap prevents nested
  parallelization

- Failed bootstrap iterations generate warnings but don't stop execution

- Confounders are removed from bootstrap data to force fresh variable
  selection

## Bootstrap Configuration

Each bootstrap iteration modifies ForestSearch arguments to:

- **Suppress output**: `details`, `showten_subgroups`, `plot.sg`,
  `plot.grf` all set to `FALSE`

- **Force re-selection**: `grf_res` and `grf_cuts` set to `NULL`

- **Prevent nested parallel**: `parallel_args$plan = "sequential"`,
  `workers = 1`

## Performance Considerations

- Typical runtime: 1-5 seconds per bootstrap iteration

- For 1000 bootstraps with 6 workers: ~3-10 minutes total

- Memory usage scales with dataset size and number of workers

- Consider reducing `nb_boots` for initial testing (e.g., 100)

## Error Handling

The function gracefully handles three failure modes:

1.  Bootstrap sample creation fails: Returns row with all `NA`

2.  ForestSearch fails to run: Warns and returns row with all `NA`

3.  ForestSearch runs but finds no subgroup: Returns row with all `NA`

All three cases ensure the foreach loop can still combine results via
`rbind`.

## See also

[`forestsearch_bootstrap_dofuture`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md)
for the wrapper function that sets up parallelization and calls this
function
[`build_cox_formula`](https://larry-leon.github.io/forestsearch/reference/build_cox_formula.md)
for creating the Cox formula
[`fit_cox_models`](https://larry-leon.github.io/forestsearch/reference/fit_cox_models.md)
for initial Cox model fitting
[`get_Cox_sg`](https://larry-leon.github.io/forestsearch/reference/get_Cox_sg.md)
for Cox model fitting on subgroups
[`get_dfRes`](https://larry-leon.github.io/forestsearch/reference/get_dfRes.md)
for processing bootstrap results into confidence intervals
[`bootstrap_ystar`](https://larry-leon.github.io/forestsearch/reference/bootstrap_ystar.md)
for generating the Ystar matrix

## Examples

``` r
if (FALSE) { # \dontrun{
# Typically called via forestsearch_bootstrap_dofuture()
# Manual usage for debugging:

# 1. Fit initial ForestSearch model
fs_result <- forestsearch(
  df.analysis = mydata,
  outcome.name = "time",
  event.name = "status",
  treat.name = "treatment",
  confounders.name = c("age", "sex", "stage")
)

# 2. Build Cox formula
cox_formula <- build_cox_formula("time", "status", "treatment")

# 3. Get observed estimates
cox_fits <- fit_cox_models(fs_result$df.est, cox_formula)

# 4. Set up parallel backend
library(doFuture)
registerDoFuture()
plan(multisession, workers = 6)

# 5. Run bootstrap (note: this is already parallelized internally)
boot_results <- bootstrap_results(
  fs.est = fs_result,
  df_boot_analysis = fs_result$df.est,
  cox.formula.boot = cox_formula,
  nb_boots = 100,
  show_three = TRUE,
  H_obs = cox_fits$H_obs,
  Hc_obs = cox_fits$Hc_obs
)

# 6. Check results
summary(boot_results)

# Proportion of bootstraps that found a subgroup
mean(!is.na(boot_results$H_biasadj_2))
} # }
```
