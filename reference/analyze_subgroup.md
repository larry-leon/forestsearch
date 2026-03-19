# Analyze subgroup for summary table (OPTIMIZED)

Analyzes a subgroup and returns formatted results for summary table.
Uses optimized cox_summary() and reduces redundant calculations.

## Usage

``` r
analyze_subgroup(
  df_sub,
  outcome.name,
  event.name,
  treat.name,
  strata.name,
  subgroup_name,
  hr_a,
  potentialOutcome.name,
  return_medians,
  N
)
```

## Arguments

- df_sub:

  Data frame for subgroup.

- outcome.name:

  Character. Name of outcome variable.

- event.name:

  Character. Name of event indicator variable.

- treat.name:

  Character. Name of treatment variable.

- strata.name:

  Character. Name of strata variable (optional).

- subgroup_name:

  Character. Subgroup name.

- hr_a:

  Character. Adjusted hazard ratio (optional).

- potentialOutcome.name:

  Character. Name of potential outcome variable (optional).

- return_medians:

  Logical. Use medians or RMST.

- N:

  Integer. Total sample size.

## Value

Character vector of results.
