# Calibrate k_treat for Target Overall Hazard Ratio

Finds the treatment effect multiplier (`k_treat`) that achieves a target
overall hazard ratio in a DGM built by
[`generate_aft_dgm_flex`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md).
The root is found by `uniroot` over `k_treat` with every other argument
to
[`generate_aft_dgm_flex()`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md)
held fixed.

## Usage

``` r
calibrate_k_treat(
  target_hr_overall,
  base_args,
  k_treat_range = c(-5, 5),
  tol = 1e-06,
  tol_rel = 2.5,
  use_ahr = FALSE,
  verbose = FALSE
)
```

## Arguments

- target_hr_overall:

  Numeric. Target overall hazard ratio (must be positive).

- base_args:

  Named list of arguments to
  [`generate_aft_dgm_flex()`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md)
  excluding `k_treat`. Every other parameter (data, covariates, model,
  subgroup_vars, seed, n_super, etc.) is taken from this list. See
  Details below for reproducibility considerations.

- k_treat_range:

  Numeric vector of length 2. Search range for `k_treat` passed to
  `uniroot`. Default `c(-5, 5)` is typically wide enough; widen if the
  boundary HRs do not bracket the target (uniroot will error otherwise).

- tol:

  Numeric. Tolerance for root finding. Default `1e-6`.

- tol_rel:

  Numeric. Maximum tolerated relative error, in percent; the call stops
  if the achieved HR is not within `tol_rel` percent of the target.
  Default `2.5`.

- use_ahr:

  Logical. If `TRUE`, calibrate to the average hazard ratio
  (`dgm$hazard_ratios$AHR`) instead of the Cox-based overall HR
  (`dgm$hazard_ratios$overall`). Default `FALSE`.

- verbose:

  Logical. Print diagnostic information. Default `FALSE`.

## Value

Numeric scalar. The calibrated `k_treat` value. Stops with an error if
the target cannot be reached – the search interval fails to bracket it
even after automatic extension, or the achieved HR lies outside
`tol_rel` percent of the target.

## Details

This is the natural sibling of
[`calibrate_k_inter`](https://larry-leon.github.io/forestsearch/reference/calibrate_k_inter.md)
(which targets the harm-subgroup HR under `model = "alt"`) and
[`calibrate_cens_adjust`](https://larry-leon.github.io/forestsearch/reference/calibrate_cens_adjust.md)
(which targets the censoring rate). It is typically used with
`model = "null"` to fix the overall HR for a uniform-effect DGM, but
also works with `model = "alt"` when the marginal HR is the calibration
target rather than the subgroup-specific HR.

Reproducibility requires the `uniroot` objective to be deterministic in
`k_treat`. This is achieved by including an integer `seed` field in
`base_args` so that every call to
[`generate_aft_dgm_flex()`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md)
during the search uses the same random-number stream. If
`base_args$seed` is `NULL` the objective becomes stochastic and
`uniroot` may not converge; the function emits a warning in this case.

## See also

[`calibrate_k_inter`](https://larry-leon.github.io/forestsearch/reference/calibrate_k_inter.md),
[`calibrate_cens_adjust`](https://larry-leon.github.io/forestsearch/reference/calibrate_cens_adjust.md),
[`generate_aft_dgm_flex`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(survival)
data(gbsg)
gbsg$time_months <- gbsg$rfstime / 30.4375

base_args <- list(
  data            = gbsg,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars     = c("meno", "grade"),
  outcome_var     = "time_months",
  event_var       = "status",
  treatment_var   = "hormon",
  model           = "null",
  n_super         = 5000,
  seed            = 99,
  verbose         = FALSE
)

# Calibrate k_treat so that the overall HR is exactly 0.70
k <- calibrate_k_treat(target_hr_overall = 0.70,
                       base_args         = base_args,
                       verbose           = TRUE)

# Verify
dgm <- do.call(generate_aft_dgm_flex,
               c(base_args, list(k_treat = k)))
print(dgm$hazard_ratios$overall)   # ~ 0.70

# Calibrate to AHR instead
k_ahr <- calibrate_k_treat(target_hr_overall = 0.70,
                           base_args         = base_args,
                           use_ahr           = TRUE)
dgm_ahr <- do.call(generate_aft_dgm_flex,
                   c(base_args, list(k_treat = k_ahr)))
print(dgm_ahr$hazard_ratios$AHR)   # ~ 0.70
} # }
```
