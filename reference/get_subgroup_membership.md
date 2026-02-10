# Get subgroup membership vector

Returns a vector indicating subgroup membership (1 if all selected
factors are present, 0 otherwise).

## Usage

``` r
get_subgroup_membership(zz, covs.in)
```

## Arguments

- zz:

  Matrix or data frame of subgroup factor indicators.

- covs.in:

  Numeric vector indicating which factors are selected (1 = included).

## Value

Numeric vector of subgroup membership (1/0).
