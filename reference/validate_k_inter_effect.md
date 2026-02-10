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

  Additional arguments passed to create_gbsg_dgm

## Value

Data frame with k_inter, hr_H, hr_Hc, AHR_H, AHR_Hc, and ratio columns

## Examples

``` r
if (FALSE) { # \dontrun{
# Test k_inter effect
results <- validate_k_inter_effect()

# k_inter = 0 should give hr_H approximately equals hr_Hc (ratio approximately 1)
} # }
```
