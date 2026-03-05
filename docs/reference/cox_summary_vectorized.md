# Cox model summary for subgroup - vectorized version

Efficiently processes multiple subgroups at once. Useful when analyzing
many subgroups (e.g., in cross-validation).

## Usage

``` r
cox_summary_vectorized(
  data,
  outcome_col,
  event_col,
  treat_col,
  strata_col = NULL,
  subgroup_col = "subgroup",
  return_format = c("formatted", "numeric")
)
```

## Arguments

- data:

  Data frame with columns for Y, E, Treat, and optionally Strata.

- outcome_col:

  Character. Name of outcome column.

- event_col:

  Character. Name of event column.

- treat_col:

  Character. Name of treatment column.

- strata_col:

  Character. Name of strata column (optional).

- subgroup_col:

  Character. Name of subgroup indicator column.

- return_format:

  Character. "formatted" or "numeric".

## Value

Data frame with one row per subgroup and HR results.
