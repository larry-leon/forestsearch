# Build candidate-membership lists from a (v1,d1,c1,v2,d2,c2) table

Build candidate-membership lists from a (v1,d1,c1,v2,d2,c2) table

## Usage

``` r
.fs_mr_family_from_table(df, tab, op_right = ">=", n_min = 1L, digits = 17L)
```

## Arguments

- df:

  Analysis data frame (for membership evaluation).

- tab:

  Data frame with columns `v1`, `d1`, `c1`, `v2`, `d2`, `c2`.

- op_right:

  Operator used when direction is not `"left"`: `">="` (DINA) or `">"`
  (GRF).

- n_min:

  Minimum subgroup size; smaller candidates are dropped.

- digits:

  Threshold label formatting (matches the package convention).

## Value

Named list of integer row-index vectors (one per candidate).
