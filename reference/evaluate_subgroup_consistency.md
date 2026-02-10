# Evaluate Single Subgroup for Consistency (Fixed-Sample)

Evaluates a single subgroup for consistency across random splits using a
fixed number of splits.

## Usage

``` r
evaluate_subgroup_consistency(
  m,
  index.Z,
  names.Z,
  df,
  found.hrs,
  n.splits,
  hr.consistency,
  pconsistency.threshold,
  pconsistency.digits = 2,
  maxk,
  confs_labels,
  details = FALSE
)
```

## Arguments

- m:

  Integer. Index of the subgroup to evaluate.

- index.Z:

  Data.table or matrix. Factor indicators for all subgroups.

- names.Z:

  Character vector. Names of factor columns.

- df:

  Data.frame. Original data with Y, Event, Treat, id columns.

- found.hrs:

  Data.table. Subgroup hazard ratio results.

- n.splits:

  Integer. Number of random splits for consistency evaluation.

- hr.consistency:

  Numeric. Minimum HR threshold for consistency.

- pconsistency.threshold:

  Numeric. Minimum proportion of splits meeting consistency.

- pconsistency.digits:

  Integer. Rounding digits for consistency proportion.

- maxk:

  Integer. Maximum number of factors in a subgroup.

- confs_labels:

  Character vector. Labels for confounders.

- details:

  Logical. Print details during execution.

## Value

Named numeric vector with consistency results, or NULL if criteria not
met.
