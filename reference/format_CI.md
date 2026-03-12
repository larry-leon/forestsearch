# Format Confidence Interval for Estimates

Formats confidence interval for estimates.

## Usage

``` r
format_CI(estimates, col_names)
```

## Arguments

- estimates:

  Data frame or data.table of estimates.

- col_names:

  Character vector of column names for estimate, lower, upper.

## Value

Character string formatted as \\estimate (lower, upper)\\.

## Examples

``` r
# \donttest{
library(data.table)
est <- data.table(hr = 1.58, lower = 0.86, upper = 2.90)
format_CI(est, col_names = c("hr", "lower", "upper"))
#> [1] "1.58 (0.86,2.9)"
# }
```
