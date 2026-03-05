# Compute Hazard Ratio for a Single Subgroup

Internal helper function to compute HR and CI for a subgroup. Uses
robust (sandwich) standard errors for consistency with cox_summary().

## Usage

``` r
compute_sg_hr(
  df,
  sg_name,
  outcome.name,
  event.name,
  treat.name,
  E.name,
  C.name,
  z_alpha = qnorm(0.975),
  conf.level = 0.95
)
```

## Arguments

- df:

  Data frame for the subgroup.

- sg_name:

  Character. Name of the subgroup.

- outcome.name:

  Character. Name of survival time variable.

- event.name:

  Character. Name of event indicator variable.

- treat.name:

  Character. Name of treatment variable.

- E.name:

  Character. Label for experimental arm.

- C.name:

  Character. Label for control arm.

- z_alpha:

  Numeric. Z-multiplier for CI (default: qnorm(0.975) for 95% CI).

- conf.level:

  Numeric. Confidence level for intervals (default: 0.95).

## Value

Data frame with single row of HR estimates, or NULL if model fails.
