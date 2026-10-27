# Resampling approximation to the splitting consistency rate

Approximates the forestsearch consistency rate for a candidate subgroup
without repeated sample splits and refits. The subgroup is fit once; the
random-split pair of half-sample treatment effects is represented as
`beta_hat +/- D`, a multiplier sum of the per-subject treatment `dfbeta`
contributions (see the file header). Handles the Cox (survival) path and
the GLM paths (binary OR/RR/RD, continuous MD, count IRR).

## Usage

``` r
consistency_resample(
  df,
  hr.consistency = 1,
  consistency_threshold = NULL,
  comparison_threshold = NULL,
  outcome_type = c("survival", "binary", "continuous", "count"),
  effect_measure = NULL,
  method = c("both", "closed", "mc"),
  multiplier = c("rademacher", "normal", "poisson"),
  draws = 1000L,
  tte.name = "Y",
  event.name = "Event",
  treat.name = "Treat",
  outcome.name = "Y",
  offset.name = NULL,
  cox_init = 0,
  adjust_covariates = NULL,
  adverse_outcome = TRUE,
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

- consistency_threshold:

  Numeric or `NULL`; GLM consistency threshold on the natural scale
  (OR/RR/IRR, or RD/MD/IRD). `NULL` uses `hr.consistency`.

- comparison_threshold:

  Numeric or `NULL`; the consistency threshold already resolved to the
  comparison scale `c` (`log(threshold)` for ratio measures, identity
  for identity measures). When supplied it is used directly and the
  internal transform is bypassed. This is the value
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  threads through for GLM outcomes (it resolves the scale upstream), so
  callers avoid a double log-transform.

- outcome_type:

  Character; `"survival"` (default), `"binary"`, `"continuous"`, or
  `"count"`.

- effect_measure:

  Character or `NULL`; GLM effect measure (resolved to the outcome-type
  default when `NULL`: binary `"OR"`, continuous `"MD"`, count `"IRR"`).
  Ignored for survival.

- method:

  Character; `"closed"`, `"mc"`, or `"both"` (default).

- multiplier:

  Character; `"rademacher"` (default), `"normal"`, `"poisson"`.

- draws:

  Integer; multiplier draws for the Monte Carlo path.

- tte.name, event.name, treat.name:

  Character; survival/treatment column names.

- outcome.name, offset.name:

  Character; GLM outcome and (rate) offset column names.

- cox_init:

  Numeric; Cox warm-start.

- adjust_covariates:

  This option is implemented but has not been tested against the
  manuscript analyses or the simulation studies of León and Anderson
  (2026a, 2026b), which use the unadjusted treatment-only fits.
  Character vector or `NULL`; adjustment terms.

- adverse_outcome:

  Logical; GLM adverse-outcome convention (see
  [`.consistency_glm_pieces()`](https://larry-leon.github.io/forestsearch/reference/dot-consistency_glm_pieces.md)).

- seed:

  Integer or `NULL`; RNG seed for the Monte Carlo draws.

## Value

List with `rate_closed`, `rate_mc`, `beta_hat`, `sigma_D`, `n`, `d`,
`delta`, `measure`; all-`NA` rates (with available summaries) if the fit
fails or the configuration is unsupported.

## Details

`multiplier = "rademacher"` (default) reproduces the 50/50 Bernoulli
split of
[`run_single_consistency_split()`](https://larry-leon.github.io/forestsearch/reference/run_single_consistency_split.md);
`"normal"` is the Lin (1997) Gaussian multiplier and `"poisson"` a
centred unit-Poisson (DBP) multiplier.

## See also

[`run_single_consistency_split()`](https://larry-leon.github.io/forestsearch/reference/run_single_consistency_split.md),
[`consistency_resample_compare()`](https://larry-leon.github.io/forestsearch/reference/consistency_resample_compare.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# survival (Cox)
consistency_resample(sub, hr.consistency = 1.0)
# binary odds ratio
consistency_resample(sub, outcome_type = "binary", effect_measure = "OR",
                     consistency_threshold = 1.0,
                     outcome.name = "Y", treat.name = "Treat")
# count incidence-rate ratio
consistency_resample(sub, outcome_type = "count", effect_measure = "IRR",
                     consistency_threshold = 1.0,
                     outcome.name = "events", offset.name = "fu_time")
} # }
```
