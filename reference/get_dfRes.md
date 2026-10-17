# Bootstrap Confidence Interval and Bias Correction Results

Calculates confidence intervals and bias-corrected estimates for
bootstrap results.

## Usage

``` r
get_dfRes(
  Hobs,
  seHobs,
  H1_adj,
  H2_adj = NULL,
  ystar,
  cov_method = "standard",
  cov_trim = 0,
  est.scale = "hr",
  est.loghr = TRUE
)
```

## Arguments

- Hobs:

  Numeric. Observed estimate.

- seHobs:

  Numeric. Standard error of observed estimate.

- H1_adj:

  Numeric. Bias-corrected estimate 1.

- H2_adj:

  Numeric. Bias-corrected estimate 2 (optional).

- ystar:

  Matrix of bootstrap samples.

- cov_method:

  Character. Covariance method ("standard" or "nocorrect").

- cov_trim:

  Numeric. Trimming proportion for covariance (default: 0.0).

- est.scale:

  Character. "hr" or "1/hr".

- est.loghr:

  Logical. Is estimate on log(HR) scale?

## Value

Data.table with confidence intervals and estimates.

## Examples

``` r
if (FALSE) { # \dontrun{
  # get_dfRes is called internally by forestsearch_bootstrap_dofuture().
  # To use directly, supply bootstrap sample matrix from bootstrap_ystar():
  #   ystar <- bootstrap_ystar(df, nb_boots = 100, seed = 42L)
  #   res <- get_dfRes(Hobs = log(1.5), seHobs = 0.3,
  #                    H1_adj = 0.05, H2_adj = 0.03,
  #                    ystar = ystar)
} # }
```
