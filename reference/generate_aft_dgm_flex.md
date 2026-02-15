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
  select_censoring = FALSE,
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

  Logical indicating whether to fit censoring model from data. If FALSE,
  censoring is not modeled. Default FALSE

- cens_type:

  Character string specifying censoring distribution: "weibull" or
  "uniform". Default is "weibull"

- cens_params:

  List of parameters for censoring distribution. For uniform: list(min =
  value, max = value). For Weibull: fitted from data

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

## Value

An object of class `c("aft_dgm_flex", "list")` containing:

- df_super:

  Data frame with the super population including all covariates, linear
  predictors, and potential outcomes

- model_params:

  List containing model parameters:

  mu

  :   Intercept from AFT model

  sigma

  :   Scale parameter

  gamma

  :   Vector of regression coefficients

  b_weibull

  :   Weibull parameterization coefficients

  b0_weibull

  :   Weibull baseline hazard coefficients

  censoring

  :   Censoring distribution parameters

- subgroup_info:

  List with subgroup information:

  vars

  :   Variables used to define subgroup

  cuts

  :   Cutpoint specifications used

  definitions

  :   Human-readable subgroup definitions

  size

  :   Number of observations in subgroup

  proportion

  :   Proportion of observations in subgroup

- hazard_ratios:

  List of true hazard ratios:

  overall

  :   Overall treatment HR

  harm_subgroup

  :   HR within subgroup (if model="alt")

  no_harm_subgroup

  :   HR outside subgroup (if model="alt")

- analysis_vars:

  List of variable classifications for analysis

- model_type:

  Character: "alt" or "null"

- n_super:

  Size of super population

- seed:

  Random seed used

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

- \\X\\ includes treatment, covariates, and treatment×subgroup
  interaction

- \\\sigma\\ is the scale parameter

- \\\epsilon\\ follows an extreme value distribution

### Interaction Term

The model creates a SINGLE interaction term representing the treatment
effect modification when ALL subgroup conditions are simultaneously
satisfied. This is not multiple separate interactions but one combined
indicator.

## References

Kalbfleisch, J.D. and Prentice, R.L. (2002). The Statistical Analysis of
Failure Time Data (2nd ed.). Wiley.

## See also

[`simulate_from_dgm`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md)
for generating simulated data from the DGM
[`find_quantile_for_proportion`](https://larry-leon.github.io/forestsearch/reference/find_quantile_for_proportion.md)
for finding quantiles that achieve target subgroup proportions

## Author

Your Name

## Examples

``` r
if (FALSE) { # \dontrun{
library(survival)
data(cancer)

# Example 1: Simple fixed cutpoints
dgm1 <- generate_aft_dgm_flex(
  data = gbsg,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = c("er", "meno"),
  subgroup_cuts = list(
    er = 20,         # Fixed value
    meno = 0         # Factor level
  ),
  model = "alt",
  verbose = TRUE
)

# Example 2: Quantile-based cutpoints
dgm2 <- generate_aft_dgm_flex(
  data = gbsg,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = c("er", "pgr", "age"),
  subgroup_cuts = list(
    er = list(type = "quantile", value = 0.25),
    pgr = list(type = "function", fun = median),
    age = list(type = "range", min = 40, max = 60)
  ),
  model = "alt",
  k_inter = 2,  # Double the interaction effect
  verbose = TRUE
)

# Print summary
print(dgm2)
} # }
```
