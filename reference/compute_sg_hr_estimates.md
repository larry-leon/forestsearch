# Compute Hazard Ratio Estimates for Subgroups

Internal function to compute Cox model hazard ratio estimates with
confidence intervals for ITT, H, and Hc subgroups.

## Usage

``` r
compute_sg_hr_estimates(
  df,
  df_H,
  df_Hc,
  outcome.name,
  event.name,
  treat.name,
  conf.level = 0.95,
  verbose = FALSE
)
```

## Arguments

- df:

  Full analysis data frame

- df_H:

  Data frame for H subgroup

- df_Hc:

  Data frame for Hc subgroup

- outcome.name:

  Character. Outcome variable name

- event.name:

  Character. Event indicator name

- treat.name:

  Character. Treatment variable name

- conf.level:

  Numeric. Confidence level

- verbose:

  Logical. Print messages

## Value

Data frame with HR estimates
