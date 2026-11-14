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
  method = c("copula", "bootstrap"),
  latent_R = NULL
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

- latent_R:

  NULL (default), or a latent correlation matrix to use in place of a
  fitted one. It must be numeric, square, with row and column names
  equal (as sets) to `vars`; it is reordered to `vars`, and checked
  symmetric with unit diagonal and positive definite. A supplied matrix
  is used as given and no fit is attempted (latentcor is not called), so
  a DGM can be rebuilt on any platform from a stored
  `dgm$mrct$latent_R`. Ignored when `method = "bootstrap"`.

## Value

data.frame with n_super rows and columns `vars`

## Details

Latent correlation (copula). With `latent_R = NULL` the latent
correlation is fitted by
[`latentcor::latentcor()`](https://rdrr.io/pkg/latentcor/man/latentcor.html)
(types "bin"/"con", method "original", nearPD) when latentcor is
installed. When latentcor is not installed, or its call fails, the
normal-score (van der Waerden) correlation of the seed rows, projected
by [`Matrix::nearPD()`](https://rdrr.io/pkg/Matrix/man/nearPD.html), is
used instead and a warning names the fallback and the reason. The two
give different super-populations from the same seed, so a DGM built
without latentcor is not the one built with it; pass `latent_R` to make
the build platform-independent. The matrix used is returned as
`attr(, "latent_R")`.
