# Create GBSG-Based AFT Data Generating Mechanism

Creates a data generating mechanism (DGM) for survival simulations based
on the German Breast Cancer Study Group (GBSG) dataset. Supports
heterogeneous treatment effects via treatment-subgroup interactions.

## Usage

``` r
create_gbsg_dgm(
  model = c("alt", "null"),
  k_treat = 1,
  k_inter = 1,
  k_z3 = 1,
  z1_quantile = 0.25,
  n_super = DEFAULT_N_SUPER,
  cens_type = c("weibull", "uniform"),
  use_rand_params = FALSE,
  seed = SEED_BASE,
  verbose = FALSE
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

## Value

A list of class "gbsg_dgm" containing:

- df_super_rand:

  Data frame with randomized super-population including potential
  outcomes (theta_0, theta_1, loghr_po)

- hr_H_true:

  Empirical hazard ratio in harm subgroup (Cox-based)

- hr_Hc_true:

  Empirical hazard ratio in complement subgroup (Cox-based)

- hr_causal:

  Overall causal (ITT) hazard ratio (Cox-based)

- AHR:

  Overall average hazard ratio (from loghr_po)

- AHR_H_true:

  Average hazard ratio in harm subgroup

- AHR_Hc_true:

  Average hazard ratio in complement subgroup

- hazard_ratios:

  List matching generate_aft_dgm_flex output format

- model_params:

  List with AFT model parameters (mu, sigma, gamma, etc.)

- cens_params:

  List with censoring model parameters

- subgroup_info:

  List with subgroup definitions and true factor names

- analysis_vars:

  Character vector of analysis variable names

- model_type:

  Character indicating "alt" or "null"

## Details

This version is aligned with
[`generate_aft_dgm_flex()`](https://github.com/larry-leon/forestsearch/reference/generate_aft_dgm_flex.md)
and
[`calculate_hazard_ratios()`](https://github.com/larry-leon/forestsearch/reference/calculate_hazard_ratios.md)
methodology, computing individual-level potential outcomes and average
hazard ratios (AHR).

### Subgroup Definition

The harm subgroup H is defined as: z1 = 1 AND z3 = 1, where:

- z1: Low estrogen receptor (ER \<= 25th percentile by default)

- z3: Premenopausal status (meno == 0)

### Model Specification

The AFT model uses covariates: treat, z1, z2, z3, z4, z5, and (for
"alt") the interaction zh = treat \* z1 \* z3.

### Interaction Effect (k_inter)

The k_inter parameter modifies the zh coefficient in the AFT model:

    gamma[zh] <- k_inter * gamma[zh]

This affects the hazard ratio for the harm subgroup:

- HR(H) = exp(-gamma\[treat\]/sigma - gamma\[zh\]/sigma)

- HR(Hc) = exp(-gamma\[treat\]/sigma)

When k_inter = 0, HR(H) = HR(Hc) (no heterogeneity).

### Alignment with generate_aft_dgm_flex

This function now computes:

- theta_0: Log-hazard contribution under control

- theta_1: Log-hazard contribution under treatment

- loghr_po: Individual causal log hazard ratio (theta_1 - theta_0)

- AHR metrics: exp(mean(loghr_po)) for overall and subgroups

## See also

[`simulate_from_gbsg_dgm`](https://github.com/larry-leon/forestsearch/reference/simulate_from_gbsg_dgm.md)
for generating data from the DGM
[`calibrate_k_inter`](https://github.com/larry-leon/forestsearch/reference/calibrate_k_inter.md)
for finding k_inter to achieve target HR

## Examples

``` r
if (FALSE) { # \dontrun{
# Alternative hypothesis with default parameters
dgm_alt <- create_gbsg_dgm(model = "alt", verbose = TRUE)

# Null hypothesis
dgm_null <- create_gbsg_dgm(model = "null", verbose = TRUE)

# Custom subgroup HR via k_inter
dgm_custom <- create_gbsg_dgm(
  model = "alt",
  k_treat = 1.2,
  k_inter = 2.0,
  verbose = TRUE
)

# Access AHR metrics (aligned with generate_aft_dgm_flex)
dgm_alt$hazard_ratios$AHR_harm
dgm_alt$hazard_ratios$AHR_no_harm
} # }
```
