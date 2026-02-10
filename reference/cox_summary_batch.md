# Batch Cox summaries with caching

For repeated calls with the same data structure but different subsets,
this version pre-processes the data structure once.

## Usage

``` r
cox_summary_batch(
  Y,
  E,
  Treat,
  Strata = NULL,
  subset_indices,
  return_format = c("formatted", "numeric")
)
```

## Arguments

- Y:

  Numeric vector of outcome (full dataset).

- E:

  Numeric vector of event indicators (full dataset).

- Treat:

  Numeric vector of treatment indicators (full dataset).

- Strata:

  Vector of strata (optional, full dataset).

- subset_indices:

  List of integer vectors, each defining a subset to analyze.

- return_format:

  Character. "formatted" or "numeric".

## Value

List of results, one per subset.
