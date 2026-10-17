# GLM Treatment-Effect Profile for Continuous Biomarker Interaction

Estimates a treatment effect profile across a continuous biomarker using
a generalized linear model with natural cubic spline interactions.
Supports binary, continuous, and count or rate outcomes. Standard errors
are computed via the delta method applied to the full model covariance,
with optional sandwich variance estimation for modified-Poisson models.

## Usage

``` r
glm_effect_profile(
  df,
  outcome.name,
  treat.name = "treat",
  z.name = "bm",
  effect_measure = c("log_OR", "log_RR", "log_IRR", "RD", "MD"),
  offset.name = NULL,
  overdispersion = c("none", "quasi", "negbin"),
  spline_df = 3L,
  z_by = 0.05,
  z_quantile = 0.95,
  alpha = 0.05,
  conf.level = NULL,
  strata.name = NULL,
  show_plot = FALSE,
  verbose = FALSE
)
```

## Arguments

- df:

  Data frame containing the analysis data.

- outcome.name:

  Character. Name of the outcome variable. Binary outcomes must be coded
  0/1; count outcomes must be non-negative integers; continuous outcomes
  may be any numeric.

- treat.name:

  Character. Name of the binary treatment indicator (0 = control, 1 =
  treated). Default: `"treat"`.

- z.name:

  Character. Name of the continuous biomarker (effect modifier).
  Default: `"bm"`.

- effect_measure:

  Character. Effect measure to estimate on the link scale. Allowed
  values: `"log_OR"` (log odds ratio), `"log_RR"` (log relative risk,
  modified Poisson), `"log_IRR"` (log incidence rate ratio, requires
  offset), `"RD"` (risk difference), `"MD"` (mean difference). Default:
  `"log_OR"`.

- offset.name:

  Character or `NULL`. Column name for person-time or exposure. When
  provided, `log(offset)` enters the linear predictor; required for
  `effect_measure = "log_IRR"`. Default: `NULL`.

- overdispersion:

  Character. Overdispersion correction for Poisson family models.
  `"none"` uses standard Poisson (model-based SEs or sandwich if
  `effect_measure = "log_RR"`); `"quasi"` uses quasi-Poisson; `"negbin"`
  uses negative binomial via
  [`MASS::glm.nb`](https://rdrr.io/pkg/MASS/man/glm.nb.html) (requires
  MASS). Default: `"none"`.

- spline_df:

  Integer. Degrees of freedom for natural spline basis. Default: `3`.

- z_by:

  Numeric. Increment for the biomarker prediction grid. Default: `0.05`.

- z_quantile:

  Numeric. Upper quantile of the biomarker used as the grid endpoint
  (avoids extrapolation into sparse tails). Default: `0.95`.

- alpha:

  Numeric. Two-sided significance level for confidence intervals.
  Default: `0.05`.

- conf.level:

  Numeric or `NULL`. Confidence level; overrides `alpha` when supplied.
  Default: `NULL`.

- strata.name:

  Character or `NULL`. Optional stratification variable entered as a
  factor covariate. Default: `NULL`.

- show_plot:

  Logical. Display a base-R diagnostic plot of the estimated treatment
  effect profile. Default: `FALSE`.

- verbose:

  Logical. Print diagnostic messages. Default: `FALSE`.

## Value

An object of class `"glm_effect_profile"`: a named list containing
`z_profile` (biomarker grid), `est` (point estimates on the link scale),
`lower` and `upper` (confidence bounds), `se` (standard errors),
`effect_measure`, `family_used`, `overdispersion`, `dispersion` (quasi
or negbin dispersion estimate, else NA), `model_fit`, `spline_basis`,
and `alpha`.

## See also

[`glm_cs_fit`](https://larry-leon.github.io/forestsearch/reference/glm_cs_fit.md)
for a simpler family-based interface,
[`cox_cs_fit`](https://larry-leon.github.io/forestsearch/reference/cox_cs_fit.md)
for survival outcomes,
[`grf.subg.harm.glm`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.glm.md)
for GRF-based subgroup identification.
