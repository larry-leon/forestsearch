# Generate Synthetic Survival Data using AFT Model with Flexible Subgroups

Creates a data generating mechanism (DGM) for survival data using an
Accelerated Failure Time (AFT) model with Weibull distribution. Supports
flexible subgroup definitions and treatment-subgroup interactions.

## Usage

``` r
generate_aft_dgm_flex(
  data,
  continuous_vars,
  factor_vars,
  continuous_vars_cens = NULL,
  factor_vars_cens = NULL,
  set_beta_spec = list(set_var = NULL, beta_var = NULL),
  outcome_var,
  event_var,
  treatment_var = NULL,
  subgroup_vars = NULL,
  subgroup_cuts = NULL,
  draw_treatment = FALSE,
  model = "alt",
  k_treat = 1,
  k_inter = 1,
  n_super = 5000,
  select_censoring = TRUE,
  cens_type = "weibull",
  cens_params = list(),
  seed = 8316951,
  verbose = TRUE,
  standardize = FALSE,
  spline_spec = NULL
)
```

## Arguments

- data:

  A data.frame containing the input dataset to base the simulation on

- continuous_vars:

  Character vector of continuous variable names to be standardized and
  included as covariates

- factor_vars:

  Character vector of factor/categorical variable names to be converted
  to dummy variables (largest value as reference)

- continuous_vars_cens:

  Character vector of continuous variable names to be used for censoring
  model. If NULL, uses same as continuous_vars. Default NULL

- factor_vars_cens:

  Character vector of factor variable names to be used for censoring
  model. If NULL, uses same as factor_vars. Default NULL

- set_beta_spec:

  List with elements 'set_var' and 'beta_var' for manually setting
  specific beta coefficients. Default list(set_var = NULL, beta_var =
  NULL)

- outcome_var:

  Character string specifying the name of the outcome/time variable

- event_var:

  Character string specifying the name of the event/status variable (1 =
  event, 0 = censored)

- treatment_var:

  Character string specifying the name of the treatment variable. If
  NULL, treatment will be randomly simulated with 50/50 allocation

- subgroup_vars:

  Character vector of variable names defining the subgroup. Default is
  NULL (no subgroups)

- subgroup_cuts:

  Named list of cutpoint specifications for subgroup variables. See
  Details section for flexible specification options

- draw_treatment:

  Logical indicating whether to redraw treatment assignment in
  simulation. Default is FALSE (use original assignments)

- model:

  Character string: "alt" for alternative model with subgroup effects,
  "null" for null model without subgroup effects. Default is "alt"

- k_treat:

  Numeric treatment effect modifier. Values \>1 increase treatment
  effect, \<1 decrease it. Default is 1 (no modification)

- k_inter:

  Numeric interaction effect modifier for treatment-subgroup
  interaction. Default is 1 (no modification)

- n_super:

  Integer specifying size of super population to generate. Default is
  5000

- select_censoring:

  Logical. If `TRUE` (default), fits the censoring distribution to the
  observed censoring times in `data` using `survreg` with AIC-based
  selection among Weibull and log-normal models (with and without
  covariates). If `FALSE`, no model is fitted; the censoring
  distribution is specified entirely by `cens_params`. Default `TRUE`.

- cens_type:

  Character string specifying censoring distribution type: `"weibull"`
  or `"uniform"`. Controls which parametric family is considered when
  `select_censoring = TRUE`, and determines the required structure of
  `cens_params` when `select_censoring = FALSE`. Default `"weibull"`.

- cens_params:

  Named list of censoring distribution parameters. Interpretation
  depends on `select_censoring` and `cens_type`:

  `select_censoring = TRUE`

  :   Ignored; all parameters are estimated from data.

  `select_censoring = FALSE, cens_type = "uniform"`

  :   Must supply `min` and `max`. If either is absent, defaults to
      `0.5 * min(y)` and `1.5 * max(y)` with a message.

  `select_censoring = FALSE, cens_type = "weibull"`

  :   Must supply `mu` (log-scale location) and `tau` (scale).
      Optionally supply `type` (`"weibull"` or `"lognormal"`); defaults
      to `"weibull"`. Censoring is treated as intercept-only (no
      covariate or treatment dependence):
      `lin_pred_cens_0 = lin_pred_cens_1 = mu`.

  Default [`list()`](https://rdrr.io/r/base/list.html).

- seed:

  Integer random seed for reproducibility. Default is 8316951

- verbose:

  Logical indicating whether to print diagnostic information during
  execution. Default is TRUE

- standardize:

  Logical indicating whether to standardize continuous variables.
  Default is FALSE

- spline_spec:

  List specifying spline configuration for treatment effect. Must
  include 'var' (variable name), 'knot', 'zeta', and 'log_hrs' (vector
  of length 3). Default NULL (no spline)

## Details

### Subgroup Cutpoint Specifications

The `subgroup_cuts` parameter accepts multiple flexible specifications:

#### Fixed Value

    subgroup_cuts = list(er = 20)  # er <= 20

#### Quantile-based

    subgroup_cuts = list(
      er = list(type = "quantile", value = 0.25)  # er <= 25th percentile
    )

#### Function-based

    subgroup_cuts = list(
      er = list(type = "function", fun = median)  # er <= median
    )

#### Range

    subgroup_cuts = list(
      age = list(type = "range", min = 40, max = 60)  # 40 <= age <= 60
    )

#### Greater than

    subgroup_cuts = list(
      nodes = list(type = "greater", quantile = 0.75)  # nodes > 75th percentile
    )

#### Multiple values (for categorical)

    subgroup_cuts = list(
      grade = list(type = "multiple", values = c(2, 3))  # grade in (2, 3)
    )

#### Custom function

    subgroup_cuts = list(
      er = list(
        type = "custom",
        fun = function(x) x <= quantile(x, 0.3) | x >= quantile(x, 0.9)
      )
    )

### Model Structure

The AFT model with Weibull distribution is specified as: \$\$\log(T) =
\mu + \gamma' X + \sigma \epsilon\$\$

Where:

- \\T\\ is the survival time

- \\\mu\\ is the intercept

- \\\gamma\\ contains the covariate effects

- \\X\\ includes treatment, covariates, and treatment x subgroup
  interaction

- \\\sigma\\ is the scale parameter

- \\\epsilon\\ follows an extreme value distribution

### Interaction Term

The model creates a SINGLE interaction term representing the treatment
effect modification when ALL subgroup conditions are simultaneously
satisfied. This is not multiple separate interactions but one combined
indicator.

## References

Leon, L.F., et al. (2024). Statistics in Medicine.

Kalbfleisch, J.D. and Prentice, R.L. (2002). The Statistical Analysis of
Failure Time Data (2nd ed.). Wiley.

## Author

Your Name

## Examples

``` r
if (FALSE) { # \dontrun{
df <- survival::gbsg
dgm <- generate_aft_dgm_flex(
  data            = df,
  outcome_var     = "rfstime",
  event_var       = "status",
  treatment_var   = "hormon",
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars     = "meno",
  model           = "null",
  verbose         = FALSE
)
str(dgm)
} # }
```
