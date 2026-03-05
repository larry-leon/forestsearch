# Create Subgroup Indicator from Factor Definitions

Parses factor definitions (e.g., "v1.1", "grade3.1") and creates a
binary indicator for subgroup membership.

## Usage

``` r
create_subgroup_indicator(df, sg_factors)
```

## Arguments

- df:

  Data frame containing the variables

- sg_factors:

  Character vector of factor definitions

## Value

Integer vector (1 = in subgroup, 0 = not in subgroup)
