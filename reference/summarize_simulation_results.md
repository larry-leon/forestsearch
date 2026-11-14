# Summarize Simulation Results

Creates a summary table of operating characteristics across all
simulations. Includes both HR and AHR metrics.

## Usage

``` r
summarize_simulation_results(
  results,
  analyses = NULL,
  digits = 2,
  digits_hr = 3
)
```

## Arguments

- results:

  data.table with simulation results from run_simulation_analysis

- analyses:

  Character vector. Analysis methods to include. Default: NULL, meaning
  all analyses in results

- digits:

  Integer. Decimal places for proportions. Default: 2

- digits_hr:

  Integer. Decimal places for hazard ratios. Default: 3

## Value

Data frame with summary statistics

## Examples

``` r
if (FALSE) { # \dontrun{
dgm <- setup_gbsg_dgm()
# run_simulation_analysis() runs one replicate per call; ten are stacked.
# confounders_base names the DGM's numeric candidate covariates, as in the
# simulations article (the default v1-v7 columns are factors).
sim_res <- data.table::rbindlist(lapply(1:10, function(i)
  run_simulation_analysis(sim_id = i, dgm = dgm, n_sample = 500,
    confounders_base = c("z1", "z2", "z3", "z4", "z5", "size", "grade3"))),
  fill = TRUE)
summarize_simulation_results(sim_res)
} # }
```
