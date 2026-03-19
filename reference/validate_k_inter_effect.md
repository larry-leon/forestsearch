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
# \donttest{
# Test k_inter effect
results <- validate_k_inter_effect()
#> Testing k_inter effect on HR heterogeneity...
#> 
#> k_inter  HR(H)    HR(Hc)   AHR(H)   AHR(Hc)  Ratio(Cox) Ratio(AHR)
#> ---------------------------------------------------------------------- 
#> -2.0     0.1336   0.6612   0.0884   0.5848   0.2021     0.1512    
#> -1.0     0.3033   0.6612   0.2274   0.5848   0.4587     0.3888    
#> 0.0      0.6552   0.6612   0.5848   0.5848   0.9909     1.0000    
#> 1.0      1.3873   0.6612   1.5041   0.5848   2.0982     2.5721    
#> 2.0      2.9651   0.6612   3.8687   0.5848   4.4846     6.6157    
#> 3.0      6.6375   0.6612   9.9507   0.5848   10.0387    17.0162   
#> 
#> PASS: k_inter = 0 gives Cox ratio ~= 1 (no heterogeneity)
#> PASS: k_inter = 0 gives AHR ratio ~= 1 (no heterogeneity)
#> 
#> AHR vs Cox HR alignment:
#>   k_inter = -2.0: HR(H) vs AHR(H) diff = 0.0452
#>   k_inter = -1.0: HR(H) vs AHR(H) diff = 0.0759
#>   k_inter = 0.0: HR(H) vs AHR(H) diff = 0.0704
#>   k_inter = 1.0: HR(H) vs AHR(H) diff = 0.1168
#>   k_inter = 2.0: HR(H) vs AHR(H) diff = 0.9036
#>   k_inter = 3.0: HR(H) vs AHR(H) diff = 3.3132

# k_inter = 0 should give hr_H approximately equals hr_Hc (ratio approximately 1)
# }
```
