# Fit Effect Models for Both Subgroups (GLM path)

GLM counterpart to
[`fit_cox_models`](https://larry-leon.github.io/forestsearch/reference/fit_cox_models.md).
Calls the estimator closure on H (treat.recommend == 0) and Hc
(treat.recommend == 1) subgroups.

## Usage

``` r
fit_effect_models(df, estimator_fn)
```

## Arguments

- df:

  Data frame with a `treat.recommend` column.

- estimator_fn:

  Closure from
  [`make_effect_estimator()`](https://larry-leon.github.io/forestsearch/reference/make_effect_estimator.md).

## Value

List with `H_obs`, `seH_obs`, `Hc_obs`, `seHc_obs`.
