# Round half away from zero

Base R [`round()`](https://rdrr.io/r/base/Round.html) uses
round-half-to-even ("banker's rounding"), so `round(0.5) == 0` and
`round(74.5) == 74`. Cut coarsening rounds representative thresholds
half *up* (away from zero) so that, e.g., a cluster centroid of 74.5
maps to 75, matching the intuitive "nearest integer" rule.

## Usage

``` r
collapse_cuts_round(x, digits = 0L)
```

## Arguments

- x:

  Numeric vector.

- digits:

  Integer. Number of decimal places.

## Value

Numeric vector rounded half away from zero.
