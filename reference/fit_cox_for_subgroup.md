# Fit Cox Model for Subgroup

Fit Cox Model for Subgroup

## Usage

``` r
fit_cox_for_subgroup(
  yy,
  dd,
  tt,
  id.x,
  df_clean = NULL,
  adjust_covariates = NULL
)
```

## Arguments

- yy, dd, tt, id.x:

  Numeric vectors: outcome, event, treatment, and the subgroup
  membership indicator (1 = in subgroup).

- df_clean:

  Data frame or `NULL`. Cleaned analysis frame row-aligned with the
  vectors, used only to source raw `adjust_covariates` columns.

- adjust_covariates:

  Character vector or `NULL`. Additional Cox terms (e.g. `"strata(x1)"`)
  appended after `Treat`. `NULL` (default) fits the treatment-only
  model.
