# Format Hazard Ratio and Confidence Interval

Formats a hazard ratio and confidence interval for display.

## Usage

``` r
hrCI_format(hrest)
```

## Arguments

- hrest:

  Numeric vector with HR, lower, and upper confidence limits.

## Value

Character string formatted as \\HR (lower, upper)\\.

## Examples

``` r
hrCI_format(c(1.45, 0.98, 2.14))
#> [1] "1.45 (0.98, 2.14)"
```
