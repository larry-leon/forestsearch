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
