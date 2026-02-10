# Format Operating Characteristics Results as GT Table

Creates a formatted gt table from simulation operating characteristics
results.

## Usage

``` r
format_oc_results(
  results,
  analyses = NULL,
  metrics = "all",
  digits = 3,
  digits_hr = 3,
  title = "Operating Characteristics Summary",
  subtitle = NULL,
  use_gt = TRUE
)
```

## Arguments

- results:

  data.table or data.frame. Simulation results from
  [`run_simulation_analysis`](https://github.com/larry-leon/forestsearch/reference/run_simulation_analysis.md)
  or combined results from multiple simulations.

- analyses:

  Character vector. Analysis methods to include. Default: NULL (all
  analyses in results)

- metrics:

  Character vector. Metrics to display. Options include: "detection",
  "classification", "hr_estimates", "subgroup_size", "all". Default:
  "all"

- digits:

  Integer. Decimal places for proportions. Default: 3

- digits_hr:

  Integer. Decimal places for hazard ratios. Default: 3

- title:

  Character. Table title. Default: "Operating Characteristics Summary"

- subtitle:

  Character. Table subtitle. Default: NULL

- use_gt:

  Logical. Return gt table if TRUE, data.frame if FALSE. Default: TRUE

## Value

A gt table object (if use_gt = TRUE and gt package available) or
data.frame

## Details

The function summarizes simulation results across multiple metrics:

- **Detection**: Proportion of simulations finding a subgroup (any.H)

- **Classification**: Sen, spec, PPV, NPV

- **HR Estimates**: Mean hazard ratios in H and Hc subgroups

- **Subgroup Size**: Average, min, max sizes
