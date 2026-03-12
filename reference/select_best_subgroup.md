# Select best subgroup based on criterion

Identifies the optimal subgroup according to the specified criterion

## Usage

``` r
select_best_subgroup(values, sg.criterion, dmin.grf, n.max)
```

## Arguments

- values:

  Data frame. Node metrics from policy trees

- sg.criterion:

  Character. "mDiff" for maximum difference, "Nsg" for largest size

- dmin.grf:

  Numeric. Minimum difference threshold

- n.max:

  Integer. Maximum allowed subgroup size (total sample size)

## Value

Data frame row with best subgroup or NULL if none found

## Examples

``` r
vals <- data.frame(diff = c(8.5, 6.2, 3.1), Nsg = c(120, 95, 80))
select_best_subgroup(values = vals, sg.criterion = "mDiff",
                     dmin.grf = 6, n.max = 500)
#>   diff Nsg
#> 1  8.5 120
```
