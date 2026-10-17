# Calibrate GLM Interaction for a Target Subgroup Effect Size

Searches over a grid of interaction multipliers (`k_inter`) to find the
value that produces a target effect size in the subgroup \\Q\\. Returns
a modified `"glm_dgm"` object with the calibrated interaction and
updated true effects.

## Usage

``` r
calibrate_glm_interaction(
  data,
  factor_vars,
  continuous_vars = NULL,
  outcome_var,
  treatment_var,
  target_effect,
  outcome_type = c("binary", "continuous", "count"),
  effect_measure = NULL,
  offset_var = NULL,
  subgroup_vars = NULL,
  subgroup_cuts = NULL,
  k_treat = 1,
  adverse_outcome = FALSE,
  k_inter_range = c(-10, 10),
  grid_step = 0.05,
  tol_rel = 2.5,
  n_super = 5000L,
  seed = 8316951L,
  verbose = FALSE
)
```

## Arguments

- data:

  The source data frame (same as passed to
  [`generate_glm_dgm`](https://larry-leon.github.io/forestsearch/reference/generate_glm_dgm.md)).

- factor_vars:

  Character vector of factor variable names.

- continuous_vars:

  Character vector of continuous prognostic variable names. Passed
  through to
  [`generate_glm_dgm`](https://larry-leon.github.io/forestsearch/reference/generate_glm_dgm.md).
  Default `NULL`.

- outcome_var:

  Character string naming the outcome variable.

- treatment_var:

  Character string naming the treatment variable.

- target_effect:

  Numeric. Target effect size in Q on the scale determined by
  `effect_measure` (e.g., OR = 2.0 for a binary outcome).

- outcome_type:

  Character. One of `"binary"`, `"continuous"`, or `"count"`.

- effect_measure:

  Character. Effect measure. Default `NULL` (selects default per
  `outcome_type`).

- offset_var:

  Character or `NULL`. Name of the exposure / follow-up time column for
  count outcomes with an offset (IRR, IRD). Default `NULL`.

- subgroup_vars:

  Character vector of subgroup-defining variables.

- subgroup_cuts:

  Named list of cutpoint specifications.

- k_treat:

  Numeric. Treatment effect scaling factor passed to
  [`generate_glm_dgm`](https://larry-leon.github.io/forestsearch/reference/generate_glm_dgm.md).
  Default `1` (preserve fitted treatment effect). Set to `0` to zero out
  the ITT effect.

- adverse_outcome:

  Logical. When `TRUE` for binary outcomes with a ratio effect measure
  (OR, RR), the calibration searches for `1 / target_effect` on the
  actual Y scale because
  [`generate_glm_dgm`](https://larry-leon.github.io/forestsearch/reference/generate_glm_dgm.md)
  negates `k_inter` internally. The returned DGM stores
  `adverse_outcome = TRUE` so that downstream functions (e.g.,
  [`run_simulation_analysis`](https://larry-leon.github.io/forestsearch/reference/run_simulation_analysis.md))
  can auto-coordinate the analysis direction. Default `FALSE`.

- k_inter_range:

  Numeric vector of length 2. Search range for `k_inter`. Default
  `c(0, 10)`.

- grid_step:

  Numeric. Grid resolution. Default `0.05`.

- tol_rel:

  Numeric. Maximum tolerated relative error, in percent; the call stops
  (fails loudly) if the achieved OR is not within `tol_rel` percent of
  the target. Default `2.5`.

- n_super:

  Integer. Super-population size. Default `5000L`.

- seed:

  Integer. Random seed. Default `8316951L`.

- verbose:

  Logical. Print calibration progress. Default `FALSE`.

## Value

An object of class `"glm_dgm"` with the calibrated `k_inter` and updated
`hazard_ratios`.

## Details

The function constructs a grid of `k_inter` values, calls
[`generate_glm_dgm`](https://larry-leon.github.io/forestsearch/reference/generate_glm_dgm.md)
for each, and selects the value whose subgroup effect is closest to
`target_effect`.

For binary outcomes with `effect_measure = "OR"`, the target is on the
switched-treatment OR scale (OR \> 1 = treatment increases the outcome).

## See also

[`generate_glm_dgm`](https://larry-leon.github.io/forestsearch/reference/generate_glm_dgm.md),
[`simulate_from_glm_dgm`](https://larry-leon.github.io/forestsearch/reference/simulate_from_glm_dgm.md)

## Examples

``` r
if (FALSE) { # \dontrun{
dgm <- calibrate_glm_interaction(
  data          = actg_df,
  factor_vars   = paste0("z", 1:12),
  outcome_var   = "y_binary",
  treatment_var = "treat",
  target_effect = 2.0,
  outcome_type  = "binary",
  subgroup_vars = c("z1", "z2"),
  subgroup_cuts = list(z1 = 1L, z2 = 1L),
  verbose       = TRUE
)
print(dgm)
} # }
```
