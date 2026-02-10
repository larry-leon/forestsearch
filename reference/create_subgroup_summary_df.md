# Create Subgroup Summary Data Frame for Forest Plot

Creates a data frame suitable for forestploter from multiple subgroup
analyses. This is a more flexible alternative for complex subgroup
configurations.

## Usage

``` r
create_subgroup_summary_df(
  df_analysis,
  subgroups,
  outcome.name,
  event.name,
  treat.name,
  E.name = "E",
  C.name = "C",
  fs_bc_list = NULL,
  fs_kfold_list = NULL,
  conf.level = 0.95
)
```

## Arguments

- df_analysis:

  Data frame. The analysis dataset.

- subgroups:

  Named list of subgroup definitions.

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

- fs_bc_list:

  List. Named list of bootstrap results for each subgroup.

- fs_kfold_list:

  List. Named list of k-fold results for each subgroup.

- conf.level:

  Numeric. Confidence level for intervals (default: 0.95).

## Value

Data frame with HR estimates for all subgroups.
