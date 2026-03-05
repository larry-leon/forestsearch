# Evaluate Consistency (Two-Stage Algorithm)

Evaluates a single subgroup for consistency using a two-stage approach:
Stage 1 screens with fewer splits, Stage 2 uses sequential batched
evaluation with early stopping for efficient evaluation.

## Usage

``` r
evaluate_consistency_twostage(
  m,
  index.Z,
  names.Z,
  df,
  found.hrs,
  hr.consistency,
  pconsistency.threshold,
  pconsistency.digits = 2,
  maxk,
  confs_labels,
  details = FALSE,
  n.splits.screen = 30,
  screen.threshold = NULL,
  n.splits.max = 400,
  batch.size = 20,
  conf.level = 0.95,
  min.valid.screen = 10
)
```

## Arguments

- m:

  Integer. Index of subgroup to evaluate.

- index.Z:

  data.table or matrix. Factor indicators for all subgroups.

- names.Z:

  Character vector. Names of factor columns.

- df:

  data.frame. Original data with Y, Event, Treat, id columns.

- found.hrs:

  data.table. Subgroup hazard ratio results.

- hr.consistency:

  Numeric. Minimum HR threshold for consistency.

- pconsistency.threshold:

  Numeric. Final consistency threshold.

- pconsistency.digits:

  Integer. Rounding digits for output.

- maxk:

  Integer. Maximum number of factors in a subgroup.

- confs_labels:

  Character vector. Labels for confounders.

- details:

  Logical. Print progress details.

- n.splits.screen:

  Integer. Number of splits for Stage 1 (default 30).

- screen.threshold:

  Numeric. Screening threshold for Stage 1 (default auto-calculated).

- n.splits.max:

  Integer. Maximum total splits (default 400).

- batch.size:

  Integer. Splits per batch in Stage 2 (default 20).

- conf.level:

  Numeric. Confidence level for early stopping (default 0.95).

- min.valid.screen:

  Integer. Minimum valid splits in Stage 1 (default 10).

## Value

Named numeric vector with consistency results, or NULL if not met.
