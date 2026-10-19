# Build Cox Model Formula

Constructs a Cox model formula from variable names, optionally adjusted
for additional covariates.

## Usage

``` r
build_cox_formula(
  outcome.name,
  event.name,
  treat.name,
  adjust_covariates = NULL
)
```

## Arguments

- outcome.name:

  Character. Name of outcome variable.

- event.name:

  Character. Name of event indicator variable.

- treat.name:

  Character. Name of treatment variable.

- adjust_covariates:

  Character vector or `NULL`. Additional model terms appended to the
  right-hand side after the treatment term. Terms are pasted verbatim,
  so survival modelling functions such as
  [`strata()`](https://rdrr.io/pkg/survival/man/strata.html) and
  [`pspline()`](https://rdrr.io/pkg/survival/man/pspline.html) may be
  used directly, e.g. `adjust_covariates = "strata(x1)"` produces a
  stratified Cox model rather than including `x1` as a linear covariate.
  `NULL` (default) reproduces the unadjusted, treatment-only formula.

## Value

An R formula object for Cox regression. When `adjust_covariates` is
supplied, the treatment term is always placed first on the right-hand
side.

## Examples

``` r
build_cox_formula("time_months", "event", "treat")
#> Surv(time_months, event) ~ treat
#> <environment: 0x55c82bab8f28>
build_cox_formula("time_months", "event", "treat",
                  adjust_covariates = c("strata(site)", "age"))
#> Surv(time_months, event) ~ treat + strata(site) + age
#> <environment: 0x55c82baf8cc0>
```
