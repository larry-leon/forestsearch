# Expand a seed covariate frame to n_super rows via a Gaussian copula.

Expand a seed covariate frame to n_super rows via a Gaussian copula.

## Usage

``` r
expand_covariates(
  seed_df,
  vars,
  n_super,
  seed = 20260903L,
  integer_vars = NULL,
  method = c("copula", "bootstrap")
)
```

## Arguments

- seed_df:

  data.frame of the seed trial (only `vars` are used)

- vars:

  covariates to expand (numeric or 0/1 binary)

- n_super:

  rows to generate

- seed:

  RNG seed

- integer_vars:

  variables to round back to integers (auto-detected if NULL)

- method:

  "copula" (default) or "bootstrap" (row resampling, as the package does
  by default – no novel covariate patterns)

## Value

data.frame with n_super rows and columns `vars`
