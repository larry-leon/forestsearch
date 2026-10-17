# Evaluate a Single Factor Combination with Status Tracking

Tests whether a specific combination meets all criteria and returns a
status code indicating how far the evaluation progressed.

## Usage

``` r
evaluate_combination_with_status(
  covs.in,
  yy,
  dd,
  tt,
  zz,
  n.min,
  d0.min,
  d1.min,
  hr.threshold,
  minp,
  rmin,
  kk,
  estimator_fn = NULL,
  df_clean = NULL,
  outcome_type = NULL,
  adjust_covariates = NULL,
  disable_effect_floor = FALSE
)
```

## Arguments

- covs.in:

  Numeric vector. Factor selection indicators.

- yy:

  Numeric vector. Outcome values.

- dd:

  Numeric vector. Event indicators.

- tt:

  Numeric vector. Treatment indicators.

- zz:

  Matrix. Factor indicators.

- n.min:

  Integer. Minimum sample size.

- d0.min:

  Integer. Minimum control events (binary/survival) or ignored
  (continuous).

- d1.min:

  Integer. Minimum treatment events (binary/survival) or ignored
  (continuous).

- hr.threshold:

  Numeric. Effect threshold, on the scale this function compares on: the
  natural hazard-ratio scale on the survival path, and the resolved
  comparison scale (link for OR/RR/IRR, identity for RD/IRD/MD) on the
  GLM paths. See
  [`subgroup.search`](https://larry-leon.github.io/forestsearch/reference/subgroup.search.md).

- minp:

  Numeric. Minimum prevalence.

- rmin:

  Integer. Minimum size reduction.

- kk:

  Integer. Combination index.

- adjust_covariates:

  Character vector or `NULL`. Cox adjustment terms passed to
  `fit_cox_for_subgroup` on the survival path.

## Value

List with:

- status:

  Integer status code: 0 = failed variance check, 1 = passed variance,
  failed prevalence, 2 = passed prevalence, failed redundancy, 3 =
  passed redundancy, failed per-arm filter, 4 = passed per-arm filter,
  failed sample size, 5 = passed sample size, failed model fit, 6 =
  passed model fit, failed effect threshold, 7 = passed all criteria
  (success)

- result:

  Result row if successful, NULL otherwise
