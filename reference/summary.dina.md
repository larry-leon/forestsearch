# Summary method for a DINA fit.

Builds a coefficient table with estimate, standard error, Wald
z-statistic, and two-sided p-value for each coefficient.

## Usage

``` r
# S3 method for class 'dina'
summary(object, ...)
```

## Arguments

- object:

  an object of class `"dina"`.

- ...:

  unused.

## Value

an object of class `"summary.dina"`, printable via
[`print()`](https://rdrr.io/r/base/print.html).
