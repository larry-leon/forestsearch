# Calibrate k_inter for Target Subgroup Hazard Ratio

Finds the interaction effect multiplier (k_inter) that achieves a target
hazard ratio in the harm subgroup.

## Usage

``` r
calibrate_k_inter(
  target_hr_harm,
  model = "alt",
  k_treat = 1,
  cens_type = "weibull",
  k_inter_range = c(-100, 100),
  tol = 1e-06,
  use_ahr = FALSE,
  verbose = FALSE,
  ...
)
```

## Arguments

- target_hr_harm:

  Numeric. Target hazard ratio for the harm subgroup

- model:

  Character. Model type ("alt" only). Default: "alt"

- k_treat:

  Numeric. Treatment effect multiplier. Default: 1

- cens_type:

  Character. Censoring type. Default: "weibull"

- k_inter_range:

  Numeric vector of length 2. Search range for k_inter. Default: c(-100,
  100)

- tol:

  Numeric. Tolerance for root finding. Default: 1e-6

- use_ahr:

  Logical. If TRUE, calibrate to AHR instead of Cox-based HR. Default:
  FALSE

- verbose:

  Logical. Print diagnostic information. Default: FALSE

- ...:

  Additional arguments passed to
  [`create_gbsg_dgm`](https://larry-leon.github.io/forestsearch/reference/create_gbsg_dgm.md)

## Value

Numeric value of k_inter that achieves the target HR

## Details

This function uses `uniroot` to find the k_inter value such that the
empirical HR (or AHR) in the harm subgroup equals target_hr_harm.

## Examples

``` r
# \donttest{
# Find k_inter for HR = 1.5 in harm subgroup
k <- calibrate_k_inter(target_hr_harm = 1.5, verbose = TRUE)
#> Calibrating k_inter to achieve HR(H) = 1.5000
#> Search range: [-100.0, 100.0]
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> HR(H) at k_inter = -100.0: 0.0000
#> HR(H) at k_inter = 100.0: 2013671808.8162
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> 
#> === Calibration Result ===
#> Found k_inter = 1.111548
#> Achieved HR(H) = 1.5000 (target: 1.5000)
#> Error: 0.000017
#> Iterations: 24

# Verify
dgm <- setup_gbsg_dgm(model = "alt", k_inter = k, verbose = FALSE)
print(dgm)
#> GBSG-Based AFT Data Generating Mechanism (Aligned)
#> ===================================================
#> 
#> Model type: alt 
#> Super-population size: 5000 
#> 
#> Effect Modifiers:
#>   k_treat: 1 
#>   k_inter: 1.111548 
#>   k_z3: 1 
#> 
#> Hazard Ratios (Cox-based, stacked PO):
#>   Overall (causal): 0.7103 
#>   Harm subgroup (H): 1.5 
#>   Complement (Hc): 0.6612 
#>   Ratio HR(H)/HR(Hc): 2.2686 
#> 
#> Average Hazard Ratios (from loghr_po):
#>   AHR (overall): 0.6681 
#>   AHR_harm (H): 1.6713 
#>   AHR_no_harm (Hc): 0.5848 
#>   Ratio AHR(H)/AHR(Hc): 2.8579 
#> 
#> Subgroup definition: z1 == 1 & z3 == 1 (low ER & premenopausal) 
#> ER threshold: 8 (quantile = 0.25)
#> Subgroup size: 634 (12.7%)
#> Analysis variables: v1, v2, v3, v4, v5, v6, v7 

# Calibrate to AHR instead
k_ahr <- calibrate_k_inter(target_hr_harm = 1.5, use_ahr = TRUE, verbose = TRUE)
#> Calibrating k_inter to achieve AHR(H) = 1.5000
#> Search range: [-100.0, 100.0]
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> AHR(H) at k_inter = -100.0: 0.0000
#> AHR(H) at k_inter = 100.0: 62477087940906073029797013052975443083264.0000
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> 
#> === Calibration Result ===
#> Found k_inter = 0.997107
#> Achieved AHR(H) = 1.5000 (target: 1.5000)
#> Error: 0.000000
#> Iterations: 20
dgm_ahr <- setup_gbsg_dgm(model = "alt", k_inter = k_ahr, verbose = FALSE)
print(dgm_ahr)
#> GBSG-Based AFT Data Generating Mechanism (Aligned)
#> ===================================================
#> 
#> Model type: alt 
#> Super-population size: 5000 
#> 
#> Effect Modifiers:
#>   k_treat: 1 
#>   k_inter: 0.997107 
#>   k_z3: 1 
#> 
#> Hazard Ratios (Cox-based, stacked PO):
#>   Overall (causal): 0.7064 
#>   Harm subgroup (H): 1.3843 
#>   Complement (Hc): 0.6612 
#>   Ratio HR(H)/HR(Hc): 2.0937 
#> 
#> Average Hazard Ratios (from loghr_po):
#>   AHR (overall): 0.659 
#>   AHR_harm (H): 1.5 
#>   AHR_no_harm (Hc): 0.5848 
#>   Ratio AHR(H)/AHR(Hc): 2.5651 
#> 
#> Subgroup definition: z1 == 1 & z3 == 1 (low ER & premenopausal) 
#> ER threshold: 8 (quantile = 0.25)
#> Subgroup size: 634 (12.7%)
#> Analysis variables: v1, v2, v3, v4, v5, v6, v7 
# }
```
