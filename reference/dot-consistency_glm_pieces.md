# Per-subject GLM pieces for the consistency approximation

Fits an unweighted, unadjusted treatment-effect model and returns the
treatment coefficient, its per-subject `dfbeta` influence contributions,
and the split-perturbation SD.

## Usage

``` r
.consistency_glm_pieces(
  df,
  outcome_type,
  effect_measure = NULL,
  treat.name = "Treat",
  outcome.name = "Y",
  offset.name = NULL,
  adjust_covariates = NULL,
  adverse_outcome = TRUE
)
```

## Arguments

- df:

  data.frame/data.table with `outcome.name`, `treat.name`, any
  `adjust_covariates`, and `offset.name` (rate measures).

- outcome_type:

  Character; `"binary"`, `"continuous"`, or `"count"`.

- effect_measure:

  Character or `NULL`; resolved to the outcome-type default when `NULL`
  (binary `"OR"`, continuous `"MD"`, count `"IRR"`).

- treat.name, outcome.name, offset.name:

  Character; column names.

- adjust_covariates:

  Character vector or `NULL`.

- adverse_outcome:

  Logical; if `FALSE`, binary `Y` is flipped to `1 - Y` and continuous
  `Y` negated, matching the estimator so that a larger coefficient
  consistently means treatment harm.

## Value

List with `beta_hat`, `dfbeta`, `sigma_D`, `n`, `d`, `log_scale`,
`measure`; or `NULL` if unsupported/degenerate/non-convergent.

## Details

Supported: `"OR"` (logistic), `"RR"` (log-binomial with Poisson
fallback), `"RD"` (identity-link binomial), `"MD"` (OLS), `"IRR"`
(Poisson with `offset(log(time))`). `"IRD"` is not a single coefficient
and returns `NULL`, and the caller falls back to splitting.

## Propensity adjustment is not represented here

This function has **no** `ps_method` or `ps_adjust_method` argument and
reads no `sw` / `ps_hat` / `ips_covar` column, so under
`ps_method != "none"` it fits the effect **as if unadjusted** rather
than declining to fit it. It cannot do otherwise: nothing in its
signature tells it which case it is in.

The consequence is that the identifier ranks candidates on the IPTW or
g-computation adjusted effect built by
[`make_effect_estimator()`](https://larry-leon.github.io/forestsearch/reference/make_effect_estimator.md),
while every quantity derived here – \\\hat\beta(g)\\, \\db\_{g,i}\\,
\\\sigma\_{D,g}\\, and therefore the whole multiplier-resampling
correction – describes the unadjusted functional.
[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
also forbids `adjust_covariates` alongside `ps_method != "none"`, so the
model fitted here is fully unadjusted, not differently adjusted.

Earlier versions of this documentation claimed that propensity-adjusted
effects "return `NULL` (the caller falls back to splitting)". That was
never true and was not implementable as written. Whether to extend the
influence path, add the argument purely so the claim can be honoured, or
refuse the combination outright is unresolved.
