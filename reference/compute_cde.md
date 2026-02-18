# Compute CDE from theta_0 and theta_1

Computes Controlled Direct Effect as the ratio of average hazard
contributions on the natural scale:
`CDE(S) = mean(exp(theta_1[S])) / mean(exp(theta_0[S]))`.

## Usage

``` r
compute_cde(df, subset_indicator = NULL)
```

## Arguments

- df:

  Data frame with `theta_0` and `theta_1` columns.

- subset_indicator:

  Optional logical/integer vector for subsetting. If provided, only rows
  where `subset_indicator == 1` are used.

## Value

Numeric CDE value, or `NA_real_` if columns are missing.
