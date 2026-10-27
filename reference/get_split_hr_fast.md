# Fast Cox Model HR Estimation

Fits a minimal Cox model to estimate hazard ratio with reduced overhead.
Disables robust variance, model matrix storage, and other extras for
speed.

## Usage

``` r
get_split_hr_fast(df, cox_init = 0, adjust_covariates = NULL)
```

## Arguments

- df:

  data.frame or data.table with Y, Event, Treat columns.

- cox_init:

  Numeric. Initial value for coefficient (default 0). Ignored when
  `adjust_covariates` is supplied (the default coxph starting values are
  used so the init length matches the number of model coefficients).

- adjust_covariates:

  This option is implemented but has not been tested against the
  manuscript analyses or the simulation studies of León and Anderson
  (2026a, 2026b), which use the unadjusted treatment-only fits.
  Character vector or `NULL`. Additional terms appended to the model
  right-hand side after `Treat`. Terms are pasted verbatim, so
  `"strata(x1)"` produces a stratified Cox model. The referenced columns
  must be present in `df`. `NULL` (default) fits the unadjusted
  `Surv(Y, Event) ~ Treat` model.

## Value

Numeric. Estimated hazard ratio for `Treat`, or NA if model fails.

## Examples

``` r
# \donttest{
set.seed(42)
df <- data.frame(
  Y     = rexp(80),
  Event = rbinom(80, 1, 0.6),
  Treat = rep(0:1, 40)
)
get_split_hr_fast(df)
#> [1] 0.7923594
# }
```
