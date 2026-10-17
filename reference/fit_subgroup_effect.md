# Estimate Subgroup Effect via Estimator Closure

GLM counterpart to
[`get_Cox_sg`](https://larry-leon.github.io/forestsearch/reference/get_Cox_sg.md).
Calls the pre-built estimator closure on a data slice and returns
`est_obs` and `se_obs` in the same format that
[`get_Cox_sg()`](https://larry-leon.github.io/forestsearch/reference/get_Cox_sg.md)
returns.

## Usage

``` r
fit_subgroup_effect(df_sg, estimator_fn)
```

## Arguments

- df_sg:

  Data frame. The subgroup data slice.

- estimator_fn:

  Closure from
  [`make_effect_estimator()`](https://larry-leon.github.io/forestsearch/reference/make_effect_estimator.md).

## Value

List with `est_obs` (effect estimate) and `se_obs` (standard error).
