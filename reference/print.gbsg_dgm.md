# Print Method for gbsg_dgm Objects

Print Method for gbsg_dgm Objects

## Usage

``` r
# S3 method for class 'gbsg_dgm'
print(x, ...)
```

## Arguments

- x:

  A gbsg_dgm object

- ...:

  Additional arguments (unused)

## Value

Invisibly returns `x`.

## Examples

``` r
# \donttest{
dgm <- setup_gbsg_dgm(model = "alt", verbose = FALSE)
print(dgm)
#> GBSG-Based AFT Data Generating Mechanism (Aligned)
#> ===================================================
#> 
#> Model type: alt 
#> Super-population size: 5000 
#> 
#> Effect Modifiers:
#>   k_treat: 1 
#>   k_inter: 1 
#>   k_z3: 1 
#> 
#> Hazard Ratios (Cox-based, stacked PO):
#>   Overall (causal): 0.7065 
#>   Harm subgroup (H): 1.3873 
#>   Complement (Hc): 0.6612 
#>   Ratio HR(H)/HR(Hc): 2.0982 
#> 
#> Average Hazard Ratios (from loghr_po):
#>   AHR (overall): 0.6592 
#>   AHR_harm (H): 1.5041 
#>   AHR_no_harm (Hc): 0.5848 
#>   Ratio AHR(H)/AHR(Hc): 2.5721 
#> 
#> Subgroup definition: z1 == 1 & z3 == 1 (low ER & premenopausal) 
#> ER threshold: 8 (quantile = 0.25)
#> Subgroup size: 634 (12.7%)
#> Analysis variables: v1, v2, v3, v4, v5, v6, v7 
# }
```
