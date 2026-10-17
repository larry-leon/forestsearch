# ForestSearch K-Fold Cross-Validation Output Summary

Summarizes cross-validation results for ForestSearch, including subgroup
agreement and performance metrics.

## Usage

``` r
forestsearch_KfoldOut(res, details = FALSE, outall = FALSE, digits = 4)
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

- digits:

  Integer. Decimal places for construction-time numeric formatting in
  the per-fold summary tables built when `outall = TRUE` (arm summaries,
  Diff, effect-estimate CI). Default 4 to provide round-down headroom
  for downstream display reformatting (e.g.,
  [`cv_summary_tables`](https://larry-leon.github.io/forestsearch/reference/cv_summary_tables.md)`(digits = ...)`).
  Threaded through
  [`SG_tab_estimates_glm`](https://larry-leon.github.io/forestsearch/reference/SG_tab_estimates_glm.md)
  for GLM outcomes. Ignored when `outall = FALSE`.

## Value

If `outall=FALSE`, a list with `sens_metrics_original` and
`find_metrics`. If `outall=TRUE`, a list with summary tables and
metrics.

## Examples

``` r
if (FALSE) { # \dontrun{
library(survival)
df <- survival::gbsg
df$grade3 <- as.integer(df$grade == "3")
fs <- forestsearch(df,
  confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
  outcome.name = "rfstime", treat.name = "hormon", event.name = "status")
oob <- forestsearch_Kfold(fs.est = fs, Kfolds = nrow(df))
result <- forestsearch_KfoldOut(oob, outall = TRUE)
} # }
```
