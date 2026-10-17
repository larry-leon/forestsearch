# HTE Test Using ForestSearch-Selected Covariates

HTE Test Using ForestSearch-Selected Covariates

## Usage

``` r
test_hte_from_forestsearch(
  fs_result,
  df,
  outcome.name,
  treat.name,
  regression = c("ols", "logistic", "poisson"),
  offset = NULL,
  poly_order = 1L
)
```

## Arguments

- fs_result:

  A `forestsearch` result object.

- df:

  Data frame used in the
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  call.

- outcome.name:

  Character. Outcome variable name.

- treat.name:

  Character. Treatment variable name.

- regression:

  Character: "ols", "logistic", or "poisson".

- offset:

  Numeric offset vector (Poisson only).

- poly_order:

  Integer. Default: 1.

## Value

List with hte_test, covariates_used, source.
