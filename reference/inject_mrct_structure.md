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
  region_loghr_shift = 0,
  region_harm_loghr_shift = 0,
  n_super = 5000L,
  expand = c("copula", "bootstrap"),
  continuous_vars_cens = NULL,
  factor_vars_cens = NULL,
  cens_type = "weibull",
  seed = 8316951L,
  verbose = FALSE,
  latent_R = NULL
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
  make_region_model(). Observed-region mode:
  `region = list(observed = "<col>", loghr = NULL)` takes the region
  from `seed_data[[col]]`, a 0/1 column (a 2-level factor is recoded as
  the other covariates are) that must have both levels in the complete
  seed rows. The region is then named `col` (`z_<col>` in `df_super`;
  `col` must not also be in `continuous_vars` or `factor_vars`), it is
  expanded with the covariates, and make_region_model() is not called.
  `prevalence`, `or_pred`, `x1_vars`, `or_x1`, `scale` and `band`
  describe the planted model, and giving any of them with `observed` is
  an error. In this mode `loghr = NULL` (the default) means the region's
  prognostic log-HR is estimated by the AFT fit; a number fixes it, as
  in the planted mode.

- x3:

  NULL, or list(vars = c(...), cuts = list(...), loghr = ...) for a
  balanced effect modifier: cuts as in generate_aft_dgm_flex (use fixed
  values), loghr = log-HR change (treatment x subgroup).

- region_treat_loghr:

  NULL, or a log-HR for a direct treat x Region interaction (the
  residual Region-U-TE pathway). Mutually exclusive with x3 (the package
  supports one interaction term).

- region_loghr_shift:

  Log hazard-ratio shift \\\gamma_R\\ added to the treatment effect of
  every region patient (default `0`, no shift). With
  `region_harm_loghr_shift` it gives the region-dependent effect
  \$\$\psi(X) = \psi\_{base}(X) + \gamma_R 1\\X \in R\\ +
  \gamma\_{H^\*R} 1\\X \in H^\*\\ 1\\X \in R\\,\$\$ where
  \\\psi\_{base}\\ is the effect built from `spline_spec` and `x3` and
  \\H^\*\\ is the `x3` subgroup (`flag_harm`). When the base effect is
  \\\log b\\ outside \\H^\*\\ and \\\log h\\ inside it, hazard ratios
  \\b_R\\ and \\h_R\\ inside the region are obtained with \\\gamma_R =
  \log(b_R/b)\\ and \\\gamma\_{H^\*R} = \log(h_R/h) - \log(b_R/b)\\. The
  shift is applied to the treated potential outcome (`theta_1`,
  `lin_pred_1`, `loghr_po`, and `lin_pred_obs` of treated rows) after
  the AFT fit, so the control arm and the prognostic structure are
  unchanged. Not available with `region_treat_loghr`.

- region_harm_loghr_shift:

  Log hazard-ratio shift \\\gamma\_{H^\*R}\\ added for patients in both
  \\H^\*\\ and the region (default `0`). A non-zero value needs `x3`
  (which defines \\H^\*\\).

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

- latent_R:

  NULL (default), or the latent correlation matrix of the copula
  expansion, passed to
  [`expand_covariates()`](https://larry-leon.github.io/forestsearch/reference/expand_covariates.md):
  named rows and columns matching the expanded covariates
  (`continuous_vars`, `factor_vars`, the censoring covariates, the `x3`
  variables and, with `region$observed`, the region column), symmetric
  and positive definite (a matrix without the observed region column
  fails the check). When supplied no latentcor fit is made, so the
  super-population is the same on every platform; e.g. pass
  `dgm$mrct$latent_R` of a stored build to rebuild it. With NULL the
  matrix is fitted by latentcor when installed; otherwise (or if the fit
  fails) the normal-score correlation is used with a warning (see
  [`expand_covariates()`](https://larry-leon.github.io/forestsearch/reference/expand_covariates.md)).
  The matrix used is stored in `dgm$mrct$latent_R` either way.

## Value

aft_dgm_flex with df_super = expanded population (+ `z_<region>`), and a
`$mrct` element documenting the injected structure.
`$mrct$region$observed` is `TRUE` in the observed-region mode, where
`prevalence_target`, `alpha0`, `coefs`, `scale`, `band` and
`$mrct$region_model` are `NULL`, `$mrct$classes$X1` is `character(0)`,
and `prognostic_loghr` is the log-HR used (the fixed `loghr`, or the AFT
estimate `model_params$b0[z_name]`). When `region_loghr_shift` or
`region_harm_loghr_shift` is non-zero, `$mrct$region_effect` records
both shifts; with both at `0` (the default) the object is exactly the
one built without them.

## Details

An observed region (`region$observed`) is carried through the expansion
with the covariates; nothing about its membership is modelled.
