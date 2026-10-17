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
[`acm.disjctif`](https://larry-leon.github.io/forestsearch/reference/acm.disjctif.md).

## Examples

``` r
df <- data.frame(age = c(40, 55, 70), grade = factor(c("1", "2", "3")))
dummy_encode(df)
#>   age grade.1 grade.2 grade.3
#> 1  40       1       0       0
#> 2  55       0       1       0
#> 3  70       0       0       1
```
