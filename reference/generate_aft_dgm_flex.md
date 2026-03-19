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

## Value

A named list of class `aft_dgm` containing:

- data:

  Simulated trial data frame with outcome, event, and treatment columns.

- model_params:

  Model parameters used for data generation (coefficients, dispersion,
  spline info if applicable).

- subgroup_info:

  Subgroup definition and membership indicators, if a heterogeneous
  treatment effect was specified.

- censoring_info:

  Censoring model parameters and observed censoring rate.

- call_args:

  Arguments used in the call, for reproducibility.

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
# \donttest{
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
#> List of 8
#>  $ df_super     :'data.frame':   5000 obs. of  33 variables:
#>   ..$ id             : int [1:5000] 1 2 3 4 5 6 7 8 9 10 ...
#>   ..$ y              : int [1:5000] 191 2172 195 286 600 1730 1264 970 624 1193 ...
#>   ..$ treat          : int [1:5000] 0 1 0 1 0 0 0 0 0 0 ...
#>   ..$ event          : num [1:5000] 1 0 0 1 1 0 0 0 1 1 ...
#>   ..$ z_age          : int [1:5000] 45 59 51 34 53 49 50 66 59 51 ...
#>   ..$ z_size         : int [1:5000] 10 8 30 30 75 21 40 28 27 35 ...
#>   ..$ z_nodes        : int [1:5000] 1 2 1 12 19 5 1 2 20 1 ...
#>   ..$ z_pgr          : int [1:5000] 14 181 119 0 375 80 80 488 9 6 ...
#>   ..$ z_er           : int [1:5000] 3 0 44 5 107 152 21 298 2 1 ...
#>   ..$ z_meno         : num [1:5000] 0 1 0 0 0 1 1 1 1 1 ...
#>   ..$ pid            : int [1:5000] 1102 821 761 359 884 1469 1393 1131 987 94 ...
#>   ..$ age            : int [1:5000] 45 59 51 34 53 49 50 66 59 51 ...
#>   ..$ meno           : int [1:5000] 0 1 0 0 0 1 1 1 1 1 ...
#>   ..$ size           : int [1:5000] 10 8 30 30 75 21 40 28 27 35 ...
#>   ..$ grade          : int [1:5000] 2 2 2 3 3 2 2 2 3 3 ...
#>   ..$ nodes          : int [1:5000] 1 2 1 12 19 5 1 2 20 1 ...
#>   ..$ pgr            : int [1:5000] 14 181 119 0 375 80 80 488 9 6 ...
#>   ..$ er             : int [1:5000] 3 0 44 5 107 152 21 298 2 1 ...
#>   ..$ zcens_age      : int [1:5000] 45 59 51 34 53 49 50 66 59 51 ...
#>   ..$ zcens_size     : int [1:5000] 10 8 30 30 75 21 40 28 27 35 ...
#>   ..$ zcens_nodes    : int [1:5000] 1 2 1 12 19 5 1 2 20 1 ...
#>   ..$ zcens_pgr      : int [1:5000] 14 181 119 0 375 80 80 488 9 6 ...
#>   ..$ zcens_er       : int [1:5000] 3 0 44 5 107 152 21 298 2 1 ...
#>   ..$ zcens_meno     : num [1:5000] 0 1 0 0 0 1 1 1 1 1 ...
#>   ..$ flag_harm      : num [1:5000] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..$ lin_pred_1     : num [1:5000] 0.547 0.743 0.664 -0.118 0.182 ...
#>   ..$ lin_pred_0     : num [1:5000] 0.263 0.459 0.38 -0.402 -0.102 ...
#>   ..$ lin_pred_obs   : num [1:5000] 0.263 0.743 0.38 -0.118 -0.102 ...
#>   ..$ theta_0        : num [1:5000] -0.362 -0.633 -0.524 0.554 0.141 ...
#>   ..$ theta_1        : num [1:5000] -0.755 -1.025 -0.916 0.162 -0.251 ...
#>   ..$ loghr_po       : num [1:5000] -0.392 -0.392 -0.392 -0.392 -0.392 ...
#>   ..$ lin_pred_cens_0: num [1:5000] 0 0 0 0 0 0 0 0 0 0 ...
#>   ..$ lin_pred_cens_1: num [1:5000] 0 0 0 0 0 0 0 0 0 0 ...
#>  $ model_params :List of 6
#>   ..$ mu         : Named num 7.52
#>   .. ..- attr(*, "names")= chr "(Intercept)"
#>   ..$ tau        : num 0.725
#>   ..$ gamma      : Named num [1:7] 0.28429 0.0075 -0.00626 -0.03902 0.00194 ...
#>   .. ..- attr(*, "names")= chr [1:7] "treat" "z_age" "z_size" "z_nodes" ...
#>   ..$ b0         : Named num [1:7] -0.39223 -0.01035 0.00864 0.05384 -0.00268 ...
#>   .. ..- attr(*, "names")= chr [1:7] "treat" "z_age" "z_size" "z_nodes" ...
#>   ..$ censoring  :List of 4
#>   .. ..$ mu   : Named num 7.45
#>   .. .. ..- attr(*, "names")= chr "(Intercept)"
#>   .. ..$ tau  : num 0.425
#>   .. ..$ gamma: Named num(0) 
#>   .. .. ..- attr(*, "names")= chr(0) 
#>   .. ..$ type : chr "weibull"
#>   ..$ spline_info: NULL
#>  $ subgroup_info:List of 5
#>   ..$ vars       : NULL
#>   ..$ cuts       : NULL
#>   ..$ definitions: list()
#>   ..$ size       : num 0
#>   ..$ proportion : num 0
#>  $ hazard_ratios:List of 7
#>   ..$ overall    : Named num 0.75
#>   .. ..- attr(*, "names")= chr "treat"
#>   ..$ AHR        : num 0.676
#>   ..$ AHR_harm   : num NaN
#>   ..$ AHR_no_harm: num 0.676
#>   ..$ CDE        : num 0.676
#>   ..$ CDE_harm   : num NA
#>   ..$ CDE_no_harm: num 0.676
#>  $ analysis_vars:List of 6
#>   ..$ continuous: chr [1:5] "age" "size" "nodes" "pgr" ...
#>   ..$ factor    : chr "meno"
#>   ..$ covariates: chr [1:6] "z_age" "z_size" "z_nodes" "z_pgr" ...
#>   ..$ treatment : chr "treat"
#>   ..$ outcome   : chr "y_sim"
#>   ..$ event     : chr "event_sim"
#>  $ model_type   : chr "null"
#>  $ n_super      : num 5000
#>  $ seed         : num 8316951
#>  - attr(*, "class")= chr [1:2] "aft_dgm_flex" "list"
# }
```
