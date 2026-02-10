# Sensitivity Analysis of Hazard Ratios to k_inter

Analyzes how the interaction parameter k_inter affects hazard ratios in
different populations (overall, harm subgroup, no-harm subgroup).

## Usage

``` r
sensitivity_analysis_k_inter(
  k_inter_range = c(-5, 5),
  n_points = 21,
  plot = TRUE,
  ...
)
```

## Arguments

- k_inter_range:

  Numeric vector of length 2 specifying the range of k_inter values to
  analyze. Default is c(-5, 5).

- n_points:

  Integer number of points to evaluate within the range. Default is 21.

- plot:

  Logical indicating whether to create visualization plots. Default is
  TRUE.

- ...:

  Additional arguments passed to
  [`generate_aft_dgm_flex`](https://github.com/larry-leon/forestsearch/reference/generate_aft_dgm_flex.md).

## Value

A data.frame of class "k_inter_sensitivity" with columns:

- k_inter:

  Numeric k_inter value

- hr_harm:

  Numeric hazard ratio in harm subgroup

- hr_no_harm:

  Numeric hazard ratio in no-harm subgroup

- hr_overall:

  Numeric overall hazard ratio

- subgroup_size:

  Integer size of harm subgroup

## Details

This function evaluates the hazard ratios at evenly spaced points across
the k_inter range. If plot = TRUE, it creates a 4-panel visualization
showing:

1.  Harm subgroup HR vs k_inter

2.  All HRs (overall, harm, no-harm) vs k_inter

3.  Ratio of HRs (harm/no-harm) showing effect modification

4.  Table of key values

## Examples

``` r
if (FALSE) { # \dontrun{
# Analyze sensitivity to k_inter
sensitivity_results <- sensitivity_analysis_k_inter(
  k_inter_range = c(-2, 2),
  n_points = 11,
  data = gbsg,
  continuous_vars = c("age", "er", "pgr"),
  factor_vars = c("meno", "grade"),
  outcome_var = "rfstime",
  event_var = "status",
  treatment_var = "hormon",
  subgroup_vars = c("er", "meno"),
  subgroup_cuts = list(er = 20, meno = 0),
  model = "alt",
  plot = TRUE
)

# Results show relationship between k_inter and HRs
print(sensitivity_results)
} # }
```
