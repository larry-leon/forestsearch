# Compute Summary Statistics for Subgroups

Internal function to compute summary statistics for each subgroup.

## Usage

``` r
compute_sg_summary(
  df,
  df_H,
  df_Hc,
  outcome.name,
  event.name,
  treat.name,
  sg0_name,
  sg1_name
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

- sg0_name:

  Character. Label for H subgroup

- sg1_name:

  Character. Label for Hc subgroup

## Value

Data frame with summary statistics
