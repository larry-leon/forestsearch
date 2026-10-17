# Standard within-subgroup treatment-effect model for a DINA subgroup

Given a subgroup discovered by
[`dina_subgroup`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md),
fit the *standard* treatment-effect model on the patients in that
subgroup: a Cox proportional-hazards model when `sg$family == "cox"`, or
a generalized linear model otherwise. The reported effect is the
treatment coefficient on the natural-parameter scale (log-HR for Cox,
log-OR for binomial, log-IRR for Poisson, mean difference for Gaussian),
directly comparable to the DINA effect returned by
[`dina_subgroup()`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md)
but obtained from a conventional within-subgroup regression rather than
from the cross-fitted CATE surface.

## Usage

``` r
dina_subgroup_refit(
  sg,
  df,
  treatment,
  outcome,
  covariates,
  status = NULL,
  strata = NULL,
  confounders = "none",
  alpha = 0.05
)
```

## Arguments

- sg:

  a `"dina_subgroup"` object from
  [`dina_subgroup`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md)
  with `sg$found == TRUE`.

- df:

  data frame the subgroup was identified on; must contain the outcome,
  treatment, covariate, and any `strata` columns.

- treatment:

  character(1); name of the binary treatment column.

- outcome:

  character(1); name of the response column. For `family = "cox"` this
  is the event/censoring time column.

- covariates:

  character vector of the DINA covariate names. Used only to construct
  the automatic adjustment set when `confounders = NULL`; ignored
  otherwise.

- status:

  character(1) or `NULL`; for Cox only, the event indicator column name.
  Required when `sg$family == "cox"`.

- strata:

  `NULL` (default) or a character vector of column names entered as Cox
  [`strata()`](https://rdrr.io/pkg/survival/man/strata.html) terms. Only
  valid for `family = "cox"`.

- confounders:

  `"none"` (default, unadjusted), `NULL` (automatic; see Details), or a
  character vector of column names.

- alpha:

  confidence level for the Wald interval. Default `0.05` for a 95% CI.

## Value

An object of class `"dina_subgroup_refit"`, a list with components:

- effect:

  treatment coefficient on the natural-parameter scale (the point
  estimate).

- se:

  Wald standard error of `effect`.

- ci:

  named length-2 numeric vector (`lower`, `upper`); Wald `1 - alpha` CI
  on the natural-parameter scale.

- effect_scale:

  character label for the scale: one of `"log-HR"`, `"log-OR"`,
  `"log-IRR"`, `"MD"`.

- ratio_scale:

  logical; `TRUE` when `effect` is a log-ratio (Cox/binomial/Poisson) so
  that [`exp()`](https://rdrr.io/r/base/Log.html) maps it to HR/OR/IRR.

- signature:

  character; the subgroup rule, e.g. `"nodes >= 12"`.

- formula:

  the fitted model formula.

- confounders_used, strata_used:

  the resolved adjustment set and stratification columns actually used.

- n_subgroup, n_treated, n_control:

  subgroup size and per-arm counts.

- n_events, n_events_treated, n_events_control:

  event counts within the subgroup (Cox only; `NA` otherwise).

- family, alpha:

  echoed inputs.

- model:

  the fitted `coxph` / `glm` object.

- call:

  the matched call.

## Details

The model is unadjusted by default (`~ treatment`); optional covariate
adjustment and stratification are available via `confounders` and
`strata`.

## Adjustment set (the confounders argument)

Three modes:

- `"none"` (default) – unadjusted: `~ treatment`.

- `NULL` – automatic: the DINA covariates supplied in `covariates`, with
  the subgroup-defining covariate (`sg$covariate`) omitted. Omitting it
  avoids adjusting on a variable whose range is restricted by the
  subgroup definition.

- a character vector – exactly those columns. A warning is issued if it
  includes `sg$covariate`.

## Inference caveat

The Wald confidence interval treats the subgroup definition as
*pre-specified*. It does **not** adjust for the data-driven selection of
the signature; it is a within-signature (conditional) interval, not a
selection-adjusted one. A selection-adjusted interval requires a
bias-corrected bootstrap of the full search procedure (Leon et al. 2024,
Section 3.2) and is not provided here.

## See also

[`dina_subgroup`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup.md)
for the subgroup search and its DINA (BLP-analog) effect;
[`dina_subgroup_bootstrap`](https://larry-leon.github.io/forestsearch/reference/dina_subgroup_bootstrap.md)
for bootstrap inference that reports this within-subgroup model
alongside the DINA effect.

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- dina(df, outcome = "time", treatment = "trt",
            covariates = c("age", "nodes", "er"),
            family = "cox", status = "status", seed = 1L)
sg  <- dina_subgroup(fit, df, covariates = c("age", "nodes", "er"),
                     m_diff = 0, n_min = 60L)

# Unadjusted within-subgroup Cox model (default):
dina_subgroup_refit(sg, df, treatment = "trt", outcome = "time",
                    covariates = c("age", "nodes", "er"),
                    status = "status")

# Adjusted for DINA covariates except the subgroup-defining one:
dina_subgroup_refit(sg, df, treatment = "trt", outcome = "time",
                    covariates = c("age", "nodes", "er"),
                    status = "status", confounders = NULL)
} # }
```
