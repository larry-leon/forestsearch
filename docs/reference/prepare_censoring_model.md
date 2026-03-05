# Prepare Censoring Model Parameters

Constructs the censoring model object and appends per-subject
counterfactual censoring linear predictors (`lin_pred_cens_0`,
`lin_pred_cens_1`) to the super-population data frame.

## Usage

``` r
prepare_censoring_model(
  df_work,
  cens_type,
  cens_params,
  df_super,
  select_censoring = TRUE,
  verbose = TRUE
)
```

## Arguments

- df_work:

  Working data frame (output of `prepare_working_dataset`).

- cens_type:

  Character. `"weibull"` or `"uniform"`.

- cens_params:

  Named list of user-supplied censoring parameters.

- df_super:

  Super-population data frame; receives `lin_pred_cens_0` and
  `lin_pred_cens_1` columns.

- select_censoring:

  Logical. If `TRUE` (default), fits the censoring distribution from
  observed data using AIC-based `survreg` model comparison. If `FALSE`,
  uses `cens_params` directly with no model fitting. See
  [`generate_aft_dgm_flex`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md)
  for the required `cens_params` structure under each combination of
  `select_censoring` and `cens_type`.

- verbose:

  Logical. If `TRUE` (default), prints the censoring model comparison
  table and recommendation. Set to `FALSE` to suppress all censoring
  model selection output.

## Value

A named list:

- cens_model:

  List of censoring distribution parameters stored in
  `dgm$model_params$censoring`.

- df_super:

  Updated super-population data frame with `lin_pred_cens_0` and
  `lin_pred_cens_1` appended. These hold covariate contributions only
  (\\\gamma_c' X\\); the intercept is excluded.

## Details

### Linear predictor convention

`lin_pred_cens_0` and `lin_pred_cens_1` store the **covariate
contribution only** — i.e. \\\gamma_c' X\\, with the intercept \\\mu_c\\
excluded. This matches the convention used for the outcome model
(`lin_pred_0`, `lin_pred_1` = \\\gamma' X\\, no intercept) computed in
[`calculate_linear_predictors()`](https://larry-leon.github.io/forestsearch/reference/calculate_linear_predictors.md).

[`simulate_from_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md)
reconstructs the full log-censoring time as: \$\$\log C = \mu_c +
\delta + \tau_c \epsilon + \gamma_c' X\$\$ where \\\mu_c\\ =
`params$censoring$mu`, \\\delta\\ = `cens_adjust`, \\\tau_c\\ =
`params$censoring$tau`, and \\\gamma_c' X\\ = `lin_pred_cens_{0|1}`.

When `select_censoring = TRUE`, `predict(survreg, type = "linear")`
returns the full linear predictor \\\mu_c + \gamma_c' X\\. The stored
intercept \\\mu_c\\ is therefore subtracted before writing
`lin_pred_cens_*`, so that
[`simulate_from_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md)
can add `params$censoring$mu` exactly once. Omitting this subtraction
causes \\\mu_c\\ to be counted twice, producing astronomically large
censoring times and universal censoring.

When `select_censoring = FALSE` with a Weibull/lognormal `cens_type`,
the intercept-only model has zero covariate contribution, so
`lin_pred_cens_0 = lin_pred_cens_1 = 0`. Storing `mu` instead of 0
causes the same double-counting.
