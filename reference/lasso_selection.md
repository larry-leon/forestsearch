# LASSO selection for Cox model

Performs LASSO variable selection using Cox regression.

## Usage

``` r
lasso_selection(
  df,
  confounders.name,
  outcome.name,
  event.name,
  seedit = 8316951
)
```

## Arguments

- df:

  Data frame.

- confounders.name:

  Character vector of confounder names.

- outcome.name:

  Character. Name of outcome variable.

- event.name:

  Character. Name of event indicator variable.

- seedit:

  Integer. Random seed.

## Value

List with selected, omitted variables, coefficients, lambda, and fits.
