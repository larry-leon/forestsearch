# Format a numeric threshold for a cut expression

Renders a (already-rounded) numeric threshold as a clean string with no
trailing zeros or floating-point noise, for use in "var \<= value" cut
expressions. Integer-valued thresholds print without a decimal point.

## Usage

``` r
collapse_cuts_fmt(v, digits = 0L)
```

## Arguments

- v:

  Numeric scalar.

- digits:

  Integer. Decimal places used during rounding.

## Value

Character scalar.
