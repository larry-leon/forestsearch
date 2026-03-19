# Prepare subgroup data for analysis

Splits a data frame into two subgroups based on a flag and treatment
scale.

## Usage

``` r
prepare_subgroup_data(df, SG_flag, est.scale, treat.name)
```

## Arguments

- df:

  Data frame.

- SG_flag:

  Character. Name of subgroup flag variable.

- est.scale:

  Character. Effect scale ("hr" or "1/hr").

- treat.name:

  Character. Name of treatment variable.

## Value

List with subgroup data frames and treatment variable name.
