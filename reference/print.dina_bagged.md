# Print a bagged DINA fit.

Extends `print.dina()` with bag count and per-bag cross-fitting info.

## Usage

``` r
# S3 method for class 'dina_bagged'
print(x, digits = max(3L, getOption("digits") - 3L), ...)
```

## Arguments

- x:

  an object of class `"dina_bagged"`.

- digits:

  number of digits to show for coefficients.

- ...:

  unused.

## Value

invisibly returns `x`.
