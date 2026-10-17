# Validate the consistency approximation against literal splitting (Cox)

Calls the package's
[`run_single_consistency_split()`](https://larry-leon.github.io/forestsearch/reference/run_single_consistency_split.md)
`R_true` times and compares the Monte Carlo rate to the closed-form and
multiplier approximations from
[`consistency_resample()`](https://larry-leon.github.io/forestsearch/reference/consistency_resample.md).
Survival path only; for the GLM paths use the standalone GLM validator.

## Usage

``` r
consistency_resample_compare(
  df,
  hr.consistency = 1,
  R_true = 400L,
  draws = 1000L,
  multiplier = "rademacher",
  tte.name = "Y",
  event.name = "Event",
  treat.name = "Treat",
  cox_init = 0,
  adjust_covariates = NULL,
  seed = NULL
)
```

## Arguments

- df:

  data.frame/data.table with the relevant outcome columns.

- hr.consistency:

  Numeric; survival HR consistency threshold (also the default GLM
  threshold when `consistency_threshold` is `NULL`). Default 1.0. This
  is `c2`, an **effect** threshold on the half-sample effect –
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)'s
  `consistency.threshold` / `hr.consistency`. It is not a proportion:
  the fraction of splits that must clear it is `pconsistency.threshold`
  (`p*`), which this function does not take – it returns the rate rather
  than thresholding it.

- R_true:

  Integer; number of literal random splits used for the truth.

- draws:

  Integer; multiplier draws for the Monte Carlo path.

- multiplier:

  Character; `"rademacher"` (default), `"normal"`, `"poisson"`.

- tte.name, event.name, treat.name:

  Character; survival/treatment column names.

- cox_init:

  Numeric; Cox warm-start.

- adjust_covariates:

  Character vector or `NULL`; adjustment terms.

- seed:

  Integer or `NULL`; RNG seed for the Monte Carlo draws.

## Value

One-row data.frame with `n`, `d`, `beta_hat`, `sigma_D`, `rate_true`,
`valid_true`, `rate_closed`, `rate_mc`, `err_closed`, `err_mc`.

## See also

[`consistency_resample()`](https://larry-leon.github.io/forestsearch/reference/consistency_resample.md),
[`run_single_consistency_split()`](https://larry-leon.github.io/forestsearch/reference/run_single_consistency_split.md)
