# Validate k_inter Effect on HR Heterogeneity

Test function to verify that k_inter properly modulates the difference
between HR(H) and HR(Hc), and that AHR metrics align with Cox-based HRs.

## Usage

``` r
validate_k_inter_effect(
  k_inter_values = c(-2, -1, 0, 1, 2, 3),
  verbose = TRUE,
  ...
)
```

## Arguments

- k_inter_values:

  Numeric vector of k_inter values to test. Default: c(-2, -1, 0, 1, 2,
  3)

- verbose:

  Logical. Print results. Default: TRUE

- ...:

  Additional arguments passed to the internal GBSG builder that
  [`setup_gbsg_dgm`](https://larry-leon.github.io/forestsearch/reference/setup_gbsg_dgm.md)
  also calls: `k_treat`, `k_z3`, `z1_quantile`, `n_super`, `cens_type`,
  `use_rand_params` and `seed`, with the meanings given there.

## Value

Data frame with k_inter, hr_H, hr_Hc, AHR_H, AHR_Hc, and ratio columns

## Note

This function builds its mechanisms through the superseded GBSG wrapper
(see
[`setup_gbsg_dgm`](https://larry-leon.github.io/forestsearch/reference/setup_gbsg_dgm.md)).
That wrapper's censoring step, in which the covariate term of the
censoring model is multiplied by the ratio of the event-model scale to
the censoring-model scale, does not affect the hazard ratios this
function reports.

## Examples

``` r
if (FALSE) { # \dontrun{
# Test k_inter effect
results <- validate_k_inter_effect()

# k_inter = 0 should give hr_H approximately equals hr_Hc (ratio approximately 1)
} # }
```
