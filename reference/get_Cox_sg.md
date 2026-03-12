# Fit Cox Model for Subgroup

Fits a Cox model for a subgroup and returns estimate and standard error.

## Usage

``` r
get_Cox_sg(df_sg, cox.formula, est.loghr = TRUE, cox_initial = log(1))
```

## Arguments

- df_sg:

  Data frame for subgroup.

- cox.formula:

  Cox model formula.

- est.loghr:

  Logical. Is estimate on log(HR) scale?

- cox_initial:

  Optional pre-fitted Cox model object to use instead of fitting a new
  model. Default NULL

## Value

List with estimate and standard error.

## Details

Function is utilized throughout codebase

## Examples

``` r
if (FALSE) { # \dontrun{
library(survival)
df <- data.frame(
  tte   = gbsg$rfstime / 30.4375,
  event = gbsg$status,
  treat = gbsg$hormon
)
formula <- build_cox_formula("tte", "event", "treat")
result  <- get_Cox_sg(df, cox.formula = formula)
exp(result$est_obs)  # hazard ratio
} # }
```
