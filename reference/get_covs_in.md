# Get indicator vector for selected subgroup factors

Returns a vector indicating which factors are included in a subgroup
combination.

## Usage

``` r
get_covs_in(
  kk,
  maxk,
  L,
  counts_1factor,
  index_1factor,
  counts_2factor = NULL,
  index_2factor = NULL,
  counts_3factor = NULL,
  index_3factor = NULL
)
```

## Arguments

- kk:

  Integer. Index of the combination.

- maxk:

  Integer. Maximum number of factors in a combination.

- L:

  Integer. Number of subgroup factors.

- counts_1factor:

  Integer. Number of single-factor combinations.

- index_1factor:

  Matrix of indices for single-factor combinations.

- counts_2factor:

  Integer. Number of two-factor combinations.

- index_2factor:

  Matrix of indices for two-factor combinations.

- counts_3factor:

  Integer. Number of three-factor combinations.

- index_3factor:

  Matrix of indices for three-factor combinations.

## Value

Numeric vector indicating selected factors (1 = included, 0 = not
included).
