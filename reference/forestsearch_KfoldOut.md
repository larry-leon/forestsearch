# ForestSearch K-Fold Cross-Validation Output Summary

Summarizes cross-validation results for ForestSearch, including subgroup
agreement and performance metrics.

## Usage

``` r
forestsearch_KfoldOut(res, details = FALSE, outall = FALSE)
```

## Arguments

- res:

  List. Result object from ForestSearch cross-validation, must contain
  elements: `cv_args`, `sg_analysis`, `sg0.name`, `sg1.name`, `Kfolds`,
  `resCV`.

- details:

  Logical. Print details during execution (default: FALSE).

- outall:

  Logical. If TRUE, returns all summary tables; if FALSE, returns only
  metrics (default: FALSE).

## Value

If `outall=FALSE`, a list with `sens_metrics_original` and
`find_metrics`. If `outall=TRUE`, a list with summary tables and
metrics.
