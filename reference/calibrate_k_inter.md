# Calibrate k_inter for Target Subgroup Hazard Ratio

Finds the interaction effect multiplier (k_inter) that achieves a target
hazard ratio in the harm subgroup.

## Usage

``` r
calibrate_k_inter(
  target_hr_harm,
  model = "alt",
  k_treat = 1,
  cens_type = "weibull",
  k_inter_range = c(-100, 100),
  tol = 1e-06,
  tol_rel = 2.5,
  use_ahr = FALSE,
  verbose = FALSE,
  ...
)
```

## Arguments

- target_hr_harm:

  Numeric. Target hazard ratio for the harm subgroup

- model:

  Character. Model type ("alt" only). Default: "alt"

- k_treat:

  Numeric. Treatment effect multiplier. Default: 1

- cens_type:

  Character. Censoring type. Default: "weibull"

- k_inter_range:

  Numeric vector of length 2. Search range for k_inter. Default: c(-100,
  100)

- tol:

  Numeric. Tolerance for root finding. Default: 1e-6

- tol_rel:

  Numeric. Maximum tolerated relative error, in percent; the call stops
  if the achieved HR is not within `tol_rel` percent of the target.
  Default `2.5`.

- use_ahr:

  Logical. If TRUE, calibrate to AHR instead of Cox-based HR. Default:
  FALSE

- verbose:

  Logical. Print diagnostic information. Default: FALSE

- ...:

  Additional arguments passed to the internal GBSG builder that
  [`setup_gbsg_dgm`](https://larry-leon.github.io/forestsearch/reference/setup_gbsg_dgm.md)
  also calls: `k_z3`, `z1_quantile`, `n_super`, `use_rand_params` and
  `seed`, with the meanings given there.

## Value

Numeric value of k_inter that achieves the target HR

## Details

This function uses `uniroot` to find the k_inter value such that the
empirical HR (or AHR) in the harm subgroup equals target_hr_harm.

## Superseded

This function is kept so that published simulations can be reproduced;
it calibrates `k_inter` for
[`setup_gbsg_dgm()`](https://larry-leon.github.io/forestsearch/reference/setup_gbsg_dgm.md).
New work should use
[`find_k_inter_for_target_hr()`](https://larry-leon.github.io/forestsearch/reference/find_k_inter_for_target_hr.md),
which calibrates `k_inter` for
[`generate_aft_dgm_flex()`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md)
on any survival dataset with a 0/1 treatment column. It targets the
Cox-based harm-subgroup hazard ratio, as `use_ahr = FALSE` does here;
there is no packaged equivalent of `use_ahr = TRUE`. The calibrated
quantity does not depend on the censoring model, so the censoring step
of the GBSG wrapper (see
[`setup_gbsg_dgm`](https://larry-leon.github.io/forestsearch/reference/setup_gbsg_dgm.md))
does not affect the calibrated value. This function will be deprecated
in a future release.

## Examples

``` r
if (FALSE) { # \dontrun{
# Find k_inter for HR = 1.5 in harm subgroup
k <- calibrate_k_inter(target_hr_harm = 1.5, verbose = TRUE)

# Verify
dgm <- setup_gbsg_dgm(model = "alt", k_inter = k, verbose = TRUE)
print(dgm$hr_H_true)  # Should be close to 1.5

# Calibrate to AHR instead
k_ahr <- calibrate_k_inter(target_hr_harm = 1.5, use_ahr = TRUE, verbose = TRUE)
dgm_ahr <- setup_gbsg_dgm(model = "alt", k_inter = k_ahr, verbose = TRUE)
print(dgm_ahr$AHR_H_true)  # Should be close to 1.5
} # }
```
