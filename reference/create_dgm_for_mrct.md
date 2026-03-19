# Create Data Generating Mechanism for MRCT Simulations

Wrapper function to create a data generating mechanism (DGM) for MRCT
simulation scenarios using
[`generate_aft_dgm_flex`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md).

## Usage

``` r
create_dgm_for_mrct(
  df_case,
  model_type = c("alt", "null"),
  log_hrs = NULL,
  confounder_var = NULL,
  confounder_effect = NULL,
  include_regA = TRUE,
  verbose = FALSE
)
```

## Arguments

- df_case:

  Data frame containing case study data

- model_type:

  Character. Either "alt" (alternative hypothesis with heterogeneous
  treatment effects) or "null" (uniform treatment effect)

- log_hrs:

  Numeric vector. Log hazard ratios for spline specification. If NULL,
  defaults are used based on model_type

- confounder_var:

  Character. Name of a confounder variable to include with a forced
  prognostic effect. Default: NULL (no forced effect)

- confounder_effect:

  Numeric. Log hazard ratio for confounder_var effect. Only used if
  confounder_var is specified

- include_regA:

  Logical. Include regA as a factor in the model. Default: TRUE

- verbose:

  Logical. Print detailed output. Default: FALSE

## Value

An object of class "aft_dgm_flex" for use with
[`simulate_from_dgm`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md)
and
[`mrct_region_sims`](https://larry-leon.github.io/forestsearch/reference/mrct_region_sims.md)

## Details

### Model Types

- alt:

  Alternative hypothesis: Treatment effect varies by biomarker level
  (heterogeneous treatment effect). Default log_hrs create HR ranging
  from 2.0 (harm) to 0.5 (benefit) across biomarker range

- null:

  Null hypothesis: Uniform treatment effect regardless of biomarker
  level. Default log_hrs = log(0.7) uniformly

### Confounder Effects

By default, NO prognostic confounder effect is forced. The
confounder_var and confounder_effect parameters allow optionally
specifying ANY baseline covariate to have a fixed prognostic effect in
the outcome model.

The regA variable (region indicator) is included as a factor by default
but without a forced effect - its coefficient is estimated from data.

## See also

[`generate_aft_dgm_flex`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md)
for underlying DGM creation
[`mrct_region_sims`](https://larry-leon.github.io/forestsearch/reference/mrct_region_sims.md)
for running simulations with the DGM
