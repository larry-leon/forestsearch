# LASSO selection for Cox model

Performs LASSO variable selection using Cox regression.

## Usage

``` r
lasso_selection(
  df,
  confounders.name,
  outcome.name,
  event.name,
  seedit = 8316951
)
```

## Arguments

- df:

  Data frame.

- confounders.name:

  Character vector of confounder names.

- outcome.name:

  Character. Name of outcome variable.

- event.name:

  Character. Name of event indicator variable.

- seedit:

  Integer. Random seed.

## Value

List with selected, omitted variables, coefficients, lambda, and fits.

## Examples

``` r
# \donttest{
library(survival)
df <- survival::gbsg
df$grade3 <- as.integer(df$grade == "3")
lasso_selection(df,
  confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
  outcome.name = "rfstime", event.name = "status")
#> $selected
#> [1] "size"   "grade3" "nodes"  "pgr"   
#> 
#> $omitted
#> [1] "age"  "meno" "er"  
#> 
#> $coefficients
#>          age         meno         size       grade3        nodes          pgr 
#>  0.000000000  0.000000000  0.005433058  0.178235777  0.049670354 -0.001813976 
#>           er 
#>  0.000000000 
#> 
#> $lambda_min
#> [1] 0.01843202
#> 
#> $cvfit
#> 
#> Call:  glmnet::cv.glmnet(x = x, y = y, family = "cox", alpha = 1) 
#> 
#> Measure: Partial Likelihood Deviance 
#> 
#>      Lambda Index Measure     SE Nonzero
#> min 0.01843    26   5.709 0.2035       4
#> 1se 0.18866     1   5.811 0.2006       0
#> 
#> $fit
#> 
#> Call:  glmnet::glmnet(x = x, y = y, family = "cox", alpha = 1, lambda = lambda_min) 
#> 
#>   Df %Dev  Lambda
#> 1  4 2.44 0.01843
#> 
# }
```
