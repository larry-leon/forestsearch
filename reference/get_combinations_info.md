# Get all combinations of subgroup factors up to maxk

Generates all possible combinations of subgroup factors up to a
specified maximum size.

## Usage

``` r
get_combinations_info(L, maxk)
```

## Arguments

- L:

  Integer. Number of subgroup factors.

- maxk:

  Integer. Maximum number of factors in a combination.

## Value

List with `max_count` (total combinations) and `indices_list` (indices
for each k).
