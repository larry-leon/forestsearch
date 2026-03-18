# Build Cox Model Formula

Constructs a Cox model formula from variable names.

## Usage

``` r
build_cox_formula(outcome.name, event.name, treat.name)
```

## Arguments

- outcome.name:

  Character. Name of outcome variable.

- event.name:

  Character. Name of event indicator variable.

- treat.name:

  Character. Name of treatment variable.

## Value

An R formula object for Cox regression.

## Examples

``` r
build_cox_formula("time_months", "event", "treat")
#> Surv(time_months, event) ~ treat
#> <environment: 0x55d10ba76698>
```
