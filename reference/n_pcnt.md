# Calculate n and percent

Returns count and percent for a vector relative to a denominator.

## Usage

``` r
n_pcnt(x, denom)
```

## Arguments

- x:

  Vector of values.

- denom:

  Denominator for percent calculation.

## Value

Character string formatted as \\n (percent%)\\.

## Examples

``` r
n_pcnt(1:30, 100)
#> [1] "30 (30.0%)"
```
