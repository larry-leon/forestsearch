# Interpret and Diagnose Search Configuration

Generates a human-readable diagnostic of how ForestSearch will interpret
the user's parameter combination. Prints via
[`message()`](https://rdrr.io/r/base/message.html) and returns an
invisible list of diagnostics.

## Usage

``` r
interpret_search_config(
  outcome_type,
  effect_measure,
  adverse_outcome,
  effect_threshold,
  consistency_threshold,
  use_lasso,
  use_grf,
  outcome.name,
  event.name,
  treat.name,
  offset.name = NULL,
  subgroup_method = "consistency",
  quiet = FALSE
)
```

## Arguments

- outcome_type:

  Character. One of `"survival"`, `"binary"`, `"continuous"`, `"count"`.

- effect_measure:

  Character. `"HR"`, `"OR"`, `"RR"`, `"RD"`, `"MD"`, `"IRR"`, `"IRD"`.

- adverse_outcome:

  Logical. As resolved by forestsearch().

- effect_threshold:

  Numeric. Resolved screening threshold (on log scale for ratio
  measures, identity for others).

- consistency_threshold:

  Numeric. Resolved consistency threshold.

- use_lasso:

  Logical.

- use_grf:

  Logical. Whether GRF candidate-cut generation is on; see
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md).

- outcome.name:

  Character. Name of outcome column.

- event.name:

  Character. Name of event column.

- treat.name:

  Character. Name of treatment column.

- offset.name:

  Character or NULL.

- subgroup_method:

  Character, the resolved `subgroup_method`. `"dina"` and `"grf"` return
  before the consistency stage and never read the consistency threshold,
  so the echoed value is annotated as unused on those paths. Display
  only; nothing is computed from it. Default `"consistency"`.

- quiet:

  Logical. If TRUE, suppress output.

## Value

Invisible list with diagnostic fields.
