# Create Compact CV Summary Text

Generates a compact text string summarizing CV results, suitable for
annotations in plots or reports.

## Usage

``` r
cv_summary_text(
  cv_result,
  est.scale = "hr",
  include_finding = TRUE,
  include_agreement = TRUE
)
```

## Arguments

- cv_result:

  List. Result from
  [`forestsearch_tenfold()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_tenfold.md)
  or
  [`forestsearch_Kfold()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_Kfold.md).

- est.scale:

  Character. "hr" or "1/hr" to determine label orientation. Default:
  "hr".

- include_finding:

  Logical. Include subgroup finding rate. Default: TRUE.

- include_agreement:

  Logical

## Value

Character string with formatted CV metrics.
