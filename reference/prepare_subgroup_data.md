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

## Examples

``` r
df <- data.frame(
  treat = c(0, 1, 0, 1, 0),
  sg_flag = c(1, 1, 0, 0, 1)
)
result <- prepare_subgroup_data(df, SG_flag = "sg_flag",
                                 est.scale = "hr", treat.name = "treat")
nrow(result$df_1)
#> [1] 3
```
