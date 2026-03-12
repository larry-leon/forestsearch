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

  Character vector. Analysis methods to include. Default: all

- digits:

  Integer. Decimal places for proportions. Default: 2

- digits_hr:

  Integer. Decimal places for hazard ratios. Default: 3

## Value

Data frame with summary statistics

## Examples

``` r
if (FALSE) { # \dontrun{
dgm <- create_gbsg_dgm()
sim_res <- run_simulation_analysis(dgm, nsim = 10)
summarize_simulation_results(sim_res)
} # }
```
