# Set Up a GBSG-Based AFT Data Generating Mechanism

Creates a GBSG-based data generating mechanism that is fully compatible
with
[`simulate_from_dgm`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md).
It builds the mechanism with the package's internal GBSG builder and
returns it as an object of class
`c("aft_dgm_flex", "gbsg_dgm", "list")`.

## Usage

``` r
setup_gbsg_dgm(
  model = c("alt", "null"),
  k_treat = 1,
  k_inter = 1,
  k_z3 = 1,
  z1_quantile = 0.25,
  n_super = 5000L,
  cens_type = c("weibull", "uniform"),
  use_rand_params = FALSE,
  seed = 8316951L,
  verbose = FALSE,
  k_random_noise = 0L,
  noise_seed = 20260807L
)
```

## Arguments

- model:

  Character. Either "alt" for alternative hypothesis with heterogeneous
  treatment effects, or "null" for uniform treatment effect. Default:
  "alt"

- k_treat:

  Numeric. Treatment effect multiplier applied to the treatment
  coefficient from the fitted AFT model. Values \> 1 strengthen the
  treatment effect. Default: 1

- k_inter:

  Numeric. Interaction effect multiplier for the treatment-subgroup
  interaction (z1 \* z3). Only used when model = "alt". Higher values
  create more heterogeneity between HR(H) and HR(Hc). Default: 1

- k_z3:

  Numeric. Effect multiplier for the z3 (menopausal status) coefficient.
  Default: 1

- z1_quantile:

  Numeric. Quantile threshold for z1 (estrogen receptor). Observations
  with ER \<= quantile are coded as z1 = 1. Default: 0.25

- n_super:

  Integer. Size of super-population for empirical HR estimation.
  Default: 5000

- cens_type:

  Character. Censoring distribution type: "weibull" or "uniform".
  Default: "weibull"

- use_rand_params:

  Logical. If TRUE, modifies confounder coefficients using estimates
  from randomized subset (meno == 0). Default: FALSE

- seed:

  Integer. Random seed for super-population generation. Default: 8316951

- verbose:

  Logical. Print diagnostic information. Default: FALSE

- k_random_noise:

  Integer. Number of standard-normal noise columns (`noise1`, `noise2`,
  ...) drawn once onto the super-population at construction. Default
  `0L` (no noise columns; output identical to previous versions).

- noise_seed:

  Integer. Seed for the population-level noise draw, used only when
  `k_random_noise > 0`. Default `20260807L`, the committed
  population-noise seed (commit `2a4787bc`).

## Value

An object of class `c("aft_dgm_flex", "gbsg_dgm", "list")` with the
builder's fields (`df_super_rand`; the true hazard ratios `hr_causal`,
`hr_H_true`, `hr_Hc_true`, `AHR`, `AHR_H_true`, `AHR_Hc_true` and the
list `hazard_ratios`; `model_params`, `cens_params`, `subgroup_info`,
`effect_modifiers`, `analysis_vars`, `model_type`, `n_super` and `seed`)
plus:

- `df_super`:

  Super-population data frame with
  [`simulate_from_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md)-compatible
  column names.

- `model_params$tau`:

  Copy of `model_params$sigma`.

- `model_params$censoring`:

  Sub-list with `type`, `mu`, `tau` for the censoring model.

- `noise_names`:

  Character vector of noise column names (`character(0)` when
  `k_random_noise = 0`).

- `noise_scheme`:

  `"population"` when noise was drawn, otherwise `"none"`.

- `noise_seed`:

  The seed used for the noise draw, or `NA_integer_` when none was
  drawn.

## Details

Internally the function calls that builder and then:

1.  Adds a `df_super` field with column names aligned to
    [`simulate_from_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md)
    conventions (`lin_pred_1`, `lin_pred_0`, `lin_pred_cens_1`,
    `lin_pred_cens_0`, `flag_harm`).

2.  Adds a `model_params$tau` field (= `model_params$sigma`) and a
    `model_params$censoring` sub-list.

3.  Sets class to `c("aft_dgm_flex", "gbsg_dgm", "list")`.

The original `df_super_rand` field is kept so that
[`compute_dgm_cde()`](https://larry-leon.github.io/forestsearch/reference/compute_dgm_cde.md)
and `print.gbsg_dgm` continue to work.

**Subgroup.** The harm subgroup H is `z1 = 1` and `z3 = 1`: `z1` is 1
when the estrogen receptor value `er` is at or below its `z1_quantile`
sample quantile in
[`survival::gbsg`](https://rdrr.io/pkg/survival/man/gbsg.html) (the 25th
percentile by default), and `z3` is 1 for premenopausal patients
(`meno == 0`). Its indicator is `flag.harm` in `df_super_rand` and
`flag_harm` in `df_super`. Under `model = "null"` the indicator is 0 for
every patient.

**Model.** A Weibull accelerated failure time model is fitted with
[`survival::survreg()`](https://rdrr.io/pkg/survival/man/survreg.html)
to [`survival::gbsg`](https://rdrr.io/pkg/survival/man/gbsg.html), with
time `rfstime / 30.4375` (months), event `status` and treatment
`treat = hormon`. Its covariates are `treat`, `z1`, `z2`, `z3`, `z4` and
`z5`, where `z2`, `z4` and `z5` are the quartile groups of `age`, `pgr`
and `nodes`, entered as scores 1 to 4; under `model = "alt"` the
interaction `zh = treat * z1 * z3` is added. `k_treat` and `k_z3`
multiply the fitted `treat` and `z3` coefficients. With AFT coefficients
`gamma` and scale `sigma`, the hazard-scale coefficients are
`-gamma / sigma` (`model_params$b_hr`).

**Interaction (`k_inter`).** Under `model = "alt"`, `k_inter` multiplies
the fitted `zh` coefficient, `gamma[zh] <- k_inter * gamma[zh]` (the
default `k_inter = 1` keeps the fitted value). Each patient's log hazard
ratio (`loghr_po = theta_1 - theta_0` in `df_super_rand`) is then
`-(gamma[treat] + gamma[zh]) / sigma` in H and `-gamma[treat] / sigma`
in its complement, so `k_inter = 0` gives every patient the same hazard
ratio. The exponentials of these two values are the stored `AHR_H_true`
and `AHR_Hc_true`. The Cox-based `hr_H_true` and `hr_Hc_true` are fitted
to stacked potential outcomes of the `n_super` patients drawn with
replacement from the GBSG rows, and differ from them; at `k_inter = 0`
they still differ slightly from each other. Under `model = "null"`,
`k_inter` is not used.

Noise columns requested via `k_random_noise` are inert N(0,1) population
attributes intended as candidate confounders; the builder never puts
them in the outcome model. They are drawn once onto the super-population
(both `df_super_rand` and `df_super`), so simulated trials and
evaluation frames inherit them by row resampling.

## Superseded

This function is kept so that published simulations can be reproduced.
New work should use
[`generate_aft_dgm_flex()`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md)
with
[`simulate_from_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md),
which take any survival dataset with a 0/1 treatment column. On the GBSG
data they give the same outcome model, true hazard ratios and `k_inter`
calibration once the caller has built the GBSG covariates, at `k_z3 = 1`
and `use_rand_params = FALSE`
([`generate_aft_dgm_flex()`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md)
has neither argument); there is no drop-in replacement call. The
censoring step differs: here (with the default `cens_type = "weibull"`)
a Weibull censoring model is fitted to the resampled super-population,
in which treatment is the arm assigned for simulation, not the observed
treatment, and its covariate term is the fitted term multiplied by the
ratio of the event-model scale to the censoring-model scale. The true
hazard ratios stored in the object and the `k_inter` calibration are not
affected; simulated censoring times are, and with them the observed
times and event indicators and anything computed from them, such as the
beta(Hhat) targets scored on a frame from
[`fs_build_eval_frame()`](https://larry-leon.github.io/forestsearch/reference/fs_build_eval_frame.md).
The "Extreme Subgroup Effects" article on the package website shows the
general machinery on the GBSG data, with a uniform treatment effect.
This function will be deprecated in a future release.

## See also

[`generate_aft_dgm_flex`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md),
[`simulate_from_dgm`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md),
[`compute_dgm_cde`](https://larry-leon.github.io/forestsearch/reference/compute_dgm_cde.md)

## Examples

``` r
if (FALSE) { # \dontrun{
dgm <- setup_gbsg_dgm(model = "alt", k_inter = 2, verbose = FALSE)
dgm <- compute_dgm_cde(dgm)
print(dgm)
sim <- simulate_from_dgm(dgm, n = 400, seed = 1)
} # }
```
