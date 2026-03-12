# Fit Cox Models for Subgroups

Fits Cox models for two subgroups defined by treatment recommendation.

## Usage

``` r
fit_cox_models(df, formula)
```

## Arguments

- df:

  Data frame.

- formula:

  Cox model formula.

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
