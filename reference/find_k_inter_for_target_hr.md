# Find k_inter Value to Achieve Target Harm Subgroup Hazard Ratio

Uses numerical root-finding to determine the interaction parameter
(k_inter) that achieves a specified target hazard ratio in the harm
subgroup. This is the most efficient method for single target
calibration.

## Usage

``` r
find_k_inter_for_target_hr(
  target_hr_harm,
  data,
  continuous_vars,
  factor_vars,
  outcome_var,
  event_var,
  treatment_var,
  subgroup_vars,
  subgroup_cuts,
  k_treat = 1,
  k_inter_range = c(-10, 10),
  tol = 0.001,
  n_super = 5000,
  verbose = TRUE
)
```

## Arguments

- target_hr_harm:

  Numeric value specifying the target hazard ratio for the harm
  subgroup. Must be positive.

- data:

  A data.frame containing the dataset to use for model fitting.

- continuous_vars:

  Character vector of continuous variable names to be standardized and
  included as covariates.

- factor_vars:

  Character vector of factor/categorical variable names to be converted
  to dummy variables.

- outcome_var:

  Character string specifying the name of the outcome/time variable.

- event_var:

  Character string specifying the name of the event/status variable (1 =
  event, 0 = censored).

- treatment_var:

  Character string specifying the name of the treatment variable.

- subgroup_vars:

  Character vector of variable names defining the subgroup.

- subgroup_cuts:

  Named list of cutpoint specifications for subgroup variables. See
  [`generate_aft_dgm_flex`](https://github.com/larry-leon/forestsearch/reference/generate_aft_dgm_flex.md)
  for details on flexible specifications.

- k_treat:

  Numeric value for treatment effect modifier. Default is 1 (no
  modification).

- k_inter_range:

  Numeric vector of length 2 specifying the search range for k_inter.
  Default is c(-10, 10).

- tol:

  Numeric value specifying tolerance for root finding convergence.
  Default is 0.001.

- n_super:

  Integer specifying size of super population for hazard ratio
  calculation. Default is 5000.

- verbose:

  Logical indicating whether to print progress information. Default is
  TRUE.

## Value

A list of class "k_inter_result" containing:

- k_inter:

  Numeric value of optimal k_inter parameter

- achieved_hr_harm:

  Numeric value of achieved hazard ratio in harm subgroup

- target_hr_harm:

  Numeric value of target hazard ratio (for reference)

- error:

  Numeric value of absolute error between achieved and target HR

- dgm:

  Object of class "aft_dgm_flex" containing the final DGM

- convergence:

  Integer number of iterations to convergence

- method:

  Character string "root-finding" indicating method used

## Details

This function uses the `uniroot` algorithm to solve the equation:
\$\$HR\_{harm}(k\_{inter}) - HR\_{target} = 0\$\$

The algorithm typically converges within 5-10 iterations and achieves
high precision (within the specified tolerance). If the root-finding
fails, the function evaluates the boundaries and provides diagnostic
information.

## See also

[`sensitivity_analysis_k_inter`](https://github.com/larry-leon/forestsearch/reference/sensitivity_analysis_k_inter.md)
for sensitivity analysis
[`generate_aft_dgm_flex`](https://github.com/larry-leon/forestsearch/reference/generate_aft_dgm_flex.md)
for DGM generation

## Examples

``` r
if (FALSE) { # \dontrun{
library(survival)
data(cancer)

# Find k_inter for target HR = 2.0 in harm subgroup
result <- find_k_inter_for_target_hr(
  target_hr_harm = 2.0,
  data = gbsg,
  continuous_vars = c("age", "er", "pgr"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = c("er", "meno"),
  subgroup_cuts = list(
    er = list(type = "quantile", value = 0.25),
    meno = 0
  ),
  k_treat = 1.0,
  verbose = TRUE
)

cat("Optimal k_inter:", result$k_inter, "\n")
cat("Achieved HR:", result$achieved_hr_harm, "\n")
} # }
```
