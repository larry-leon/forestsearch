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

## Examples

``` r
zz <- matrix(c(1, 0, 1, 1, 0, 1), nrow = 3,
             dimnames = list(NULL, c("er_le0", "grade3")))
get_subgroup_membership(zz, covs.in = c(1, 1))  # both factors required
#> [1]  TRUE FALSE  TRUE
get_subgroup_membership(zz, covs.in = c(1, 0))  # first factor only
#> [1]  TRUE FALSE  TRUE
```
