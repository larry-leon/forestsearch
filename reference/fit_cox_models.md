# Fit Cox Models for Subgroups

Fits Cox models for two subgroups defined by treatment recommendation.

## Usage

``` r
fit_cox_models(df, formula, treat.name = NULL)
```

## Arguments

- df:

  Data frame.

- formula:

  Cox model formula.

- treat.name:

  Character or `NULL`. Passed to
  [`get_Cox_sg`](https://larry-leon.github.io/forestsearch/reference/get_Cox_sg.md)
  to extract the treatment coefficient by name when `formula` is
  adjusted for covariates. `NULL` (default) preserves the treatment-only
  behaviour.

## Value

List with HR and SE for each subgroup.

## Examples

``` r
if (FALSE) { # \dontrun{
library(survival)
df <- data.frame(
  tte            = gbsg$rfstime / 30.4375,
  event          = gbsg$status,
  treat          = gbsg$hormon,
  treat.recommend = as.integer(gbsg$er > 0)
)
formula <- build_cox_formula("tte", "event", "treat")
fit_cox_models(df, formula)
} # }
```
