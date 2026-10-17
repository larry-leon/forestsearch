# Build an aft_dgm_flex with synthetic regions and paper-style effect structure.

Build an aft_dgm_flex with synthetic regions and paper-style effect
structure.

## Usage

``` r
inject_mrct_structure(
  seed_data,
  outcome_var,
  event_var,
  treatment_var,
  continuous_vars,
  factor_vars,
  x_pred,
  spline_spec,
  region = list(),
  x3 = NULL,
  region_treat_loghr = NULL,
  n_super = 5000L,
  expand = c("copula", "bootstrap"),
  continuous_vars_cens = NULL,
  factor_vars_cens = NULL,
  cens_type = "weibull",
  seed = 8316951L,
  verbose = FALSE
)
```

## Arguments

- seed_data:

  seed trial data.frame (one row per patient)

- outcome_var, event_var, treatment_var:

  column names in seed_data

- continuous_vars, factor_vars:

  covariates entering the AFT model (factor_vars must be 0/1 numeric or
  2-level factors for clean expansion)

- x_pred:

  region-associated effect modifier (X2), a continuous_var

- spline_spec:

  list(knot, zeta, log_hrs) on the NATURAL scale of x_pred: true log-HR
  is piecewise linear through (0, `log_hrs[1]`), (knot, `log_hrs[2]`),
  (zeta, `log_hrs[3]`). log_hrs constant =\> no heterogeneity (paper's
  beta1 = 0).

- region:

  list(name = "region", prevalence = 0.2, or_pred = 10, x1_vars =
  character(0), or_x1 = 1, loghr = 0, scale = "rank") – loghr is the
  prognostic log-HR of Region (no TE modification; 0 for the paper's
  DGM; -log(5) reproduces the Feb-2026 deck's strongly prognostic
  region); scale = "rank" (ECDF) or "minmax" (the paper's `[0,1]`
  scaling); optional band = list(cut, or) (and optional var, default
  `x_pred`) adds `log(or) * 1(var <= cut)` to the region logit, see
  make_region_model()

- x3:

  NULL, or list(vars = c(...), cuts = list(...), loghr = ...) for a
  balanced effect modifier: cuts as in generate_aft_dgm_flex (use fixed
  values), loghr = log-HR change (treatment x subgroup).

- region_treat_loghr:

  NULL, or a log-HR for a direct treat x Region interaction (the
  residual Region-U-TE pathway). Mutually exclusive with x3 (the package
  supports one interaction term).

- n_super:

  size of the expanded super-population

- expand:

  "copula" or "bootstrap"

- continuous_vars_cens, factor_vars_cens, cens_type:

  passed through

- seed:

  RNG seed

- verbose:

  logical

## Value

aft_dgm_flex with df_super = expanded population (+ `z_<region>`), and a
`$mrct` element documenting the injected structure.
