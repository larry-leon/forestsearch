# Target Estimate and Standard Error for Bootstrap

Calculates target estimate and standard error for bootstrap samples.

## Usage

``` r
get_targetEst(x, ystar, cov_method = "standard", cov_trim = 0)
```

## Arguments

- x:

  Numeric vector of estimates.

- ystar:

  Matrix of bootstrap samples.

- cov_method:

  Character. Covariance method ("standard" or "nocorrect").

- cov_trim:

  Numeric. Trimming proportion for covariance (default: 0.0).

## Value

List with target estimate, standard errors, and correction term.
