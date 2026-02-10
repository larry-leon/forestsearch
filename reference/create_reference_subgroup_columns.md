# Create Reference Subgroup Indicator Columns

Creates indicator columns (0/1) in the data frame for each reference
subgroup based on the provided subset expressions.

## Usage

``` r
create_reference_subgroup_columns(df, ref_subgroups, verbose = FALSE)
```

## Arguments

- df:

  Data frame to modify.

- ref_subgroups:

  Named list of reference subgroup definitions. Each element should have
  `subset_expr` and optionally `label` and `color`.

- verbose:

  Logical. Print diagnostic messages.

## Value

List with modified df, cols, labels, and colors vectors.
