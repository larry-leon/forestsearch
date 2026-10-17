# Predicted DINA values for new covariate data.

Returns the estimated treatment-effect contrast \\\hat\tau(x) =
\hat\beta_0 + x^\top \hat\beta\_{1:d}\\, on the natural parameter scale
of the fitted family (mean difference for Gaussian, log odds ratio for
binomial, log mean ratio for Poisson, log hazard ratio for Cox).

## Usage

``` r
# S3 method for class 'dina'
predict(object, newdata, ...)
```

## Arguments

- object:

  an object of class `"dina"`.

- newdata:

  a numeric matrix or data frame with the same number of columns as the
  training `X` and columns in the same order.

- ...:

  unused.

## Value

numeric vector of predicted treatment-effect contrasts, one per row of
`newdata`.
