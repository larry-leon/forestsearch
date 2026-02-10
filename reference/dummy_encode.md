# Dummy-code a data frame (numeric pass-through, factors expanded)

Dummy-code a data frame (numeric pass-through, factors expanded)

## Usage

``` r
dummy_encode(df)
```

## Arguments

- df:

  Data frame with numeric and/or factor columns.

## Value

Data frame with numeric columns unchanged and factor columns expanded
via
[`acm.disjctif`](https://github.com/larry-leon/forestsearch/reference/acm.disjctif.md).
