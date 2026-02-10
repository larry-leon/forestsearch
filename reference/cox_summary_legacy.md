# Cox model summary for subgroup

Called in analyze_subgroup() \<– SG_tab_estimates

## Usage

``` r
cox_summary_legacy(Y, E, Treat, Strata)
```

## Arguments

- Y:

  Numeric vector of outcome.

- E:

  Numeric vector of event indicators.

- Treat:

  Numeric vector of treatment indicators.

- Strata:

  Vector of strata (optional).

## Value

Character string with formatted HR and CI.

## Details

Calculates hazard ratio and confidence interval for a subgroup using Cox
regression.
