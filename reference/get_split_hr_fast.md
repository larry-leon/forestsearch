# Fast Cox Model HR Estimation

Fits a minimal Cox model to estimate hazard ratio with reduced overhead.
Disables robust variance, model matrix storage, and other extras for
speed.

## Usage

``` r
get_split_hr_fast(df, cox_init = 0)
```

## Arguments

- df:

  data.frame or data.table with Y, Event, Treat columns.

- cox_init:

  Numeric. Initial value for coefficient (default 0).

## Value

Numeric. Estimated hazard ratio, or NA if model fails.

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
#>     Treat 
#> 0.7923594 
# }
```
