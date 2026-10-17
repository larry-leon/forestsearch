# Generate GLM-Based Data Generating Mechanism

Creates a data generating mechanism (DGM) for non-survival outcomes
using generalised linear models.

## Usage

``` r
generate_glm_dgm(
  data,
  factor_vars,
  continuous_vars = NULL,
  outcome_var,
  treatment_var,
  outcome_type = c("binary", "continuous", "count"),
  effect_measure = NULL,
  offset_var = NULL,
  subgroup_vars = NULL,
  subgroup_cuts = NULL,
  model = c("alt", "null"),
  k_treat = 1,
  k_inter = 0,
  adverse_outcome = FALSE,
  n_super = 5000L,
  seed = 8316951L,
  verbose = FALSE,
  k_random_noise = 0L,
  noise_seed = 20260807L
)
```

## Arguments

- data:

  A data frame containing the source dataset.

- factor_vars:

  Character vector of factor/categorical variable names to include as
  covariates. These are the candidate confounders for ForestSearch
  (e.g., `paste0("z", 1:12)`).

- continuous_vars:

  Character vector of continuous variable names to include as prognostic
  covariates in the baseline GLM. These improve the model fit but do not
  participate in subgroup definition. Consistent with
  [`generate_aft_dgm_flex`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md),
  which accepts both `continuous_vars` and `factor_vars`. Default `NULL`
  (factor-only model, backward compatible).

- outcome_var:

  Character string naming the outcome variable.

- treatment_var:

  Character string naming the treatment variable (must be 0/1 coded).

- outcome_type:

  Character. One of `"binary"`, `"continuous"`, or `"count"`.

- effect_measure:

  Character. Effect measure for the GLM. Default is `NULL`, which
  selects `"OR"` for binary, `"MD"` for continuous, `"IRR"` for count
  outcomes.

- offset_var:

  Character or `NULL`. Name of the exposure / follow-up time column for
  count outcomes with an offset. Required when `effect_measure` is
  `"IRR"` or `"IRD"`. The offset enters the Poisson GLM as
  `offset(log(offset_var))`. Default `NULL`.

- subgroup_vars:

  Character vector of variable names defining the subgroup. Default
  `NULL` (no subgroup).

- subgroup_cuts:

  Named list of cutpoint specifications for subgroup variables. Uses the
  same flexible format as
  [`generate_aft_dgm_flex`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md):
  fixed values, quantiles, or functions.

- model:

  Character. `"alt"` for alternative hypothesis with treatment–subgroup
  interaction, `"null"` for uniform treatment effect. Default `"alt"`.

- k_treat:

  Numeric. Scaling factor for the fitted treatment coefficient on the
  linear predictor scale. `k_treat = 1` (default) preserves the fitted
  treatment effect as-is. `k_treat > 1` amplifies the treatment main
  effect, useful when the source dataset has a weak treatment effect and
  a stronger ITT is desired for simulation. `k_treat < 1` attenuates the
  effect. Analogous to the `k_treat` parameter in
  [`generate_aft_dgm_flex`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md).
  Note: the null DGM (`model = "null"`) represents a *homogeneous*
  treatment effect (no HTE), not a zero treatment effect. `k_treat`
  scales the magnitude of that homogeneous effect.

- k_inter:

  Numeric. Direct additive shift on the linear predictor for Q members
  under treatment. The interpretation depends on `outcome_type`:

  binary

  :   Shift on the log-odds scale. Positive increases P(Y=1) for Q under
      treatment.

  continuous

  :   Shift on the mean scale (identity link). Positive increases
      \\E(Y)\\ for Q under treatment.

  count

  :   Shift on the log-rate scale. Positive increases \\E(Y)\\ for Q
      under treatment (multiplicative: rate multiplied by
      `exp(k_inter)`).

  Default `0` (no interaction, equivalent to `model = "null"`).

- adverse_outcome:

  Logical. If `TRUE` and `outcome_type = "binary"`, the interaction
  shift `k_inter` is negated internally so that a positive `k_inter`
  consistently amplifies the treatment contrast on the beneficial
  (1 - Y) scale. The flag is stored in the returned object for
  downstream use by
  [`calibrate_glm_interaction`](https://larry-leon.github.io/forestsearch/reference/calibrate_glm_interaction.md)
  and
  [`run_simulation_analysis`](https://larry-leon.github.io/forestsearch/reference/run_simulation_analysis.md).
  Default `FALSE`.

- n_super:

  Integer. Size of the super-population. Default `5000L`.

- seed:

  Integer. Random seed for super-population sampling. Default
  `8316951L`.

- verbose:

  Logical. Print diagnostic information. Default `FALSE`.

- k_random_noise:

  Integer. Number of standard-normal noise columns (`noise1`, `noise2`,
  ...) drawn once onto the super-population at construction. Default
  `0L` (no noise columns; output identical to previous versions).

- noise_seed:

  Integer. Seed for the population-level noise draw, used only when
  `k_random_noise > 0`. Default `20260807L`, the committed
  population-noise seed (commit `2a4787bc`).

## Value

An object of class `c("glm_dgm", "list")` with:

- `df_super`:

  Super-population data frame with potential outcome columns: `p0`, `p1`
  (binary probabilities), or `mu0`, `mu1` (conditional means for
  continuous / count outcomes), and `flag_harm`.

- `df_source`:

  The observed analysis frame in observed order (every source patient
  exactly once), carrying the same potential-outcome columns as
  `df_super` computed by the identical prediction calls, plus
  `flag_harm` and the discretised covariates. This is the frozen panel
  consumed by `simulate_from_glm_dgm(baseline = "fixed")`.

- `hazard_ratios`:

  Named list with `harm_subgroup`, `no_harm_subgroup`, and `overall` –
  effect estimates on the scale determined by `effect_measure`. Field
  name retained for compatibility with
  [`get_dgm_hr()`](https://larry-leon.github.io/forestsearch/reference/get_dgm_hr.md)
  and reporting functions.

- `outcome_type`:

  Character: `"binary"`, `"continuous"`, or `"count"`.

- `effect_measure`:

  Character: the effect measure used.

- `model_params`:

  List with fitted coefficients, family, `k_inter`, residual SD
  (continuous), and offset variable name (count).

- `subgroup_info`:

  List with subgroup definition, size, and proportion.

- `model_type`:

  Character: `"alt"` or `"null"`.

- `noise_names`:

  Character vector of noise column names (`character(0)` when
  `k_random_noise = 0`).

- `noise_scheme`:

  `"population"` when noise was drawn, otherwise `"none"`.

- `noise_seed`:

  The seed used for the noise draw, or `NA_integer_` when none was
  drawn.

## Details

Fits a GLM to the source data, creates a super-population with
individual-level potential outcomes, and stores true subgroup effects.
The returned object can be passed to
[`simulate_from_glm_dgm`](https://larry-leon.github.io/forestsearch/reference/simulate_from_glm_dgm.md)
for simulation studies, and is compatible with
[`run_simulation_analysis`](https://larry-leon.github.io/forestsearch/reference/run_simulation_analysis.md)
(which dispatches on class).

Noise columns requested via `k_random_noise` are inert N(0,1) population
attributes intended as candidate confounders; the builder never puts
them in the outcome model. They are drawn once onto `df_super` after
sampling, so simulated trials inherit them by row resampling.

## See also

[`simulate_from_glm_dgm`](https://larry-leon.github.io/forestsearch/reference/simulate_from_glm_dgm.md),
[`calibrate_glm_interaction`](https://larry-leon.github.io/forestsearch/reference/calibrate_glm_interaction.md),
[`generate_aft_dgm_flex`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(speff2trial)
actg <- subset(ACTG175, arms %in% c(1, 3))
actg$treat  <- 1L - as.integer(actg$arms == 1)
actg$y      <- as.integer(actg$cd420 > actg$cd40)
actg$z1     <- as.factor(ifelse(actg$age > 34, 1L, 0L))
actg$z2     <- as.factor(ifelse(actg$preanti <= 745, 1L, 0L))

dgm <- generate_glm_dgm(
  data          = actg,
  factor_vars   = c("z1", "z2"),
  outcome_var   = "y",
  treatment_var = "treat",
  outcome_type  = "binary",
  subgroup_vars = c("z1", "z2"),
  subgroup_cuts = list(z1 = 1L, z2 = 1L),
  model         = "alt",
  k_inter       = 1.0,   # additive log-odds shift for Q under treatment
  verbose       = TRUE
)
} # }
```
