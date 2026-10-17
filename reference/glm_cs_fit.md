# Fit GLM with Spline Treatment-Biomarker Interaction

Estimates treatment effects as a function of a continuous covariate
using a generalized linear model with natural cubic splines. This is the
GLM analog of
[`cox_cs_fit`](https://larry-leon.github.io/forestsearch/reference/cox_cs_fit.md)
for binary, continuous, and count outcomes.

## Usage

``` r
glm_cs_fit(
  df,
  outcome_name = "event",
  treat_name = "treat",
  z_name = "bm",
  family = "binomial",
  offset_name = NULL,
  alpha = 0.05,
  spline_df = 3,
  z_by = 0.05,
  z_quantile = 0.95,
  z_max = Inf,
  z_window = 0,
  show_plot = FALSE,
  plot_params = NULL,
  truebeta_name = NULL,
  verbose = TRUE
)

# S3 method for class 'glm_cs_fit'
print(x, ...)
```

## Arguments

- df:

  Data frame containing outcome data.

- outcome_name:

  Character. Name of the outcome variable. Default: `"event"`.

- treat_name:

  Character. Name of the treatment variable (0/1). Default: `"treat"`.

- z_name:

  Character. Name of the continuous biomarker covariate. Default:
  `"bm"`.

- family:

  Character or family object. One of `"binomial"` (default),
  `"poisson"`, or `"gaussian"`. Determines the link function and the
  scale of the treatment-effect profile: `"binomial"` uses logit link
  (profile is log-OR); `"poisson"` uses log link (profile is log-IRR,
  requires `offset_name`); `"gaussian"` uses identity link (profile is
  mean difference).

- offset_name:

  Character or `NULL`. Name of the offset variable (e.g., log follow-up
  time for Poisson models). The offset is applied as `log(offset_var)`
  internally. Required when `family = "poisson"`.

- alpha:

  Numeric. Significance level for confidence intervals (two-sided).
  Default: 0.05 (95 percent CI).

- spline_df:

  Integer. Degrees of freedom for the natural spline. Default: 3.

- z_by:

  Numeric. Increment for the biomarker prediction grid. Default: 0.05.

- z_quantile:

  Numeric (0–1). Upper quantile for the prediction grid. Default: 0.95
  (5th to 95th percentile).

- z_max:

  Numeric. Maximum z value for predictions. Default: `Inf` (no
  truncation).

- z_window:

  Numeric. Half-width for counting observations near each grid point.
  Default: 0.0.

- show_plot:

  Logical. Display a base-R diagnostic plot. Default: `FALSE`.

- plot_params:

  List. Optional plot parameter overrides (see
  [`cox_cs_fit`](https://larry-leon.github.io/forestsearch/reference/cox_cs_fit.md)
  for format).

- truebeta_name:

  Character or `NULL`. Name of a column containing true log-effect
  values (for simulation validation).

- verbose:

  Logical. Print diagnostic information. Default: `TRUE`.

- x:

  A `glm_cs_fit` object.

- ...:

  Additional arguments (unused).

## Value

A list of class `"glm_cs_fit"` containing:

- z_profile:

  Numeric vector. Biomarker grid values.

- loghr_est:

  Numeric vector. Point estimates of the log treatment effect (log-OR,
  log-IRR, or MD) at each grid point.

- loghr_lower:

  Numeric vector. Lower confidence bounds.

- loghr_upper:

  Numeric vector. Upper confidence bounds.

- se_loghr:

  Numeric vector. Standard errors (delta method).

- counts_profile:

  Integer vector. Observation counts near each grid point.

- glm_primary:

  Numeric. Overall treatment effect from the no-interaction GLM.

- model_fit:

  The fitted `glm` object.

- spline_basis:

  The natural spline basis object.

- family:

  Character. The GLM family used.

- effect_label:

  Character. Human-readable label for the effect scale (e.g., "log(OR)",
  "log(IRR)", "Mean Difference").

- lrt_pvalue:

  Numeric. P-value from the likelihood ratio test comparing the
  interaction model to the main-effects-only model.

- alpha:

  Numeric. Significance level used.

- ci_level:

  Numeric. Confidence level (1 - alpha).

Invisibly returns `x`.

## Details

### Model Structure

For binary outcomes (`family = "binomial"`):
\$\$\text{logit}(P(Y=1\|Z,A)) = \beta_0 A + f(Z) + g(Z) \cdot A\$\$

For count outcomes with offset (`family = "poisson"`):
\$\$\log(E\[Y\|Z,A\]) = \beta_0 A + f(Z) + g(Z) \cdot A + \log(t)\$\$

For continuous outcomes (`family = "gaussian"`): \$\$E\[Y\|Z,A\] =
\beta_0 A + f(Z) + g(Z) \cdot A\$\$

The treatment-effect profile \\\beta(Z) = \beta_0 + g(Z)\\ is estimated
on the link scale (log-OR, log-IRR, or identity) with pointwise
delta-method confidence intervals.

## See also

[`cox_cs_fit`](https://larry-leon.github.io/forestsearch/reference/cox_cs_fit.md)
for the survival analog,
[`glm_effect_profile`](https://larry-leon.github.io/forestsearch/reference/glm_effect_profile.md)
for the extended interface with sandwich SEs and overdispersion
correction.

## Examples

``` r
if (FALSE) { # \dontrun{
# Binary outcome: log-OR profile over a biomarker
set.seed(42)
df <- data.frame(
  event = rbinom(500, 1, 0.3),
  treat = rbinom(500, 1, 0.5),
  bm    = rnorm(500, 0, 1)
)
result <- glm_cs_fit(df, outcome_name = "event", z_name = "bm",
                      family = "binomial")

# Poisson outcome with offset: log-IRR profile
df$ftime <- rexp(500, 0.01) + 1
result_irr <- glm_cs_fit(df, outcome_name = "event", z_name = "bm",
                           family = "poisson", offset_name = "ftime")
} # }
```
