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
  [`generate_aft_dgm_flex`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md).

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
