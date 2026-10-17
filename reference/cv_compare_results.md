# Compare Multiple CV Results

Creates a comparison table from multiple cross-validation runs with
different configurations.

## Usage

``` r
cv_compare_results(
  cv_list,
  metrics = c("all", "finding", "agreement"),
  show_percentages = TRUE,
  digits = 1,
  use_gt = TRUE
)
```

## Arguments

- cv_list:

  Named list of cv_result objects from
  [`forestsearch_tenfold()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_tenfold.md)
  or
  [`forestsearch_Kfold()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_Kfold.md).

- metrics:

  Character vector. Which metrics to include. Options: "finding",
  "agreement", "all". Default: "all".

- show_percentages:

  Logical. Display as percentages. Default: TRUE.

- digits:

  Integer. Decimal places. Default: 1.

- use_gt:

  Logical. Return gt table if TRUE. Default: TRUE.

## Value

A gt table or data.frame comparing CV results across configurations.

## Examples

``` r
if (FALSE) { # \dontrun{
# Compare CV results from different configurations
cv_comparison <- cv_compare_results(
  cv_list = list(
    "maxk=1" = cv_maxk1,
    "maxk=2" = cv_maxk2,
    "maxk=3" = cv_maxk3
  )
)
} # }
```
