# Truth by region on the super-population: prevalence, x_pred shift, causal AHR / CDE, and a large-sample Cox HR on UNCENSORED potential outcomes. Compare HR_cox_PO with the censored-trial HR from sim_region_metrics(): a gap between them isolates the prognostic-shift / follow-up-truncation mechanism (the Feb-2026 deck) from covariate shift (AHR itself differing by region).

Truth by region on the super-population: prevalence, x_pred shift,
causal AHR / CDE, and a large-sample Cox HR on UNCENSORED potential
outcomes. Compare HR_cox_PO with the censored-trial HR from
sim_region_metrics(): a gap between them isolates the prognostic-shift /
follow-up-truncation mechanism (the Feb-2026 deck) from covariate shift
(AHR itself differing by region).

## Usage

``` r
mrct_truth_by_region(dgm, seed = 1L)
```

## Arguments

- dgm:

  An MRCT data-generating mechanism returned by
  [`inject_mrct_structure()`](https://larry-leon.github.io/forestsearch/reference/inject_mrct_structure.md).
  Uses its super-population `df_super`, the region indicator and
  effect-modifier names in `mrct`, and the AFT parameters `mu` and `tau`
  in `model_params`.

- seed:

  Integer seed passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) before drawing the
  common error term used to generate the uncensored potential outcomes.
  Default `1L`.

## Value

A data frame with one row per stratum (`"ALL"`, `"NON-REGION"`,
`"REGION"`) and columns `stratum`, `n`, `prev`, `x_pred_mean`, `AHR`,
`CDE` and `HR_cox_PO`. Attribute `"smd_x_pred"` holds the standardised
mean difference of the effect modifier between region and non-region.
