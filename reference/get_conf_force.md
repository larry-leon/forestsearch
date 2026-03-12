# Get forced cut expressions for variables

For each variable in `conf.force.names`, returns cut expressions if
continuous.

## Usage

``` r
get_conf_force(df, conf.force.names, cont.cutoff = 4)
```

## Arguments

- df:

  Data frame.

- conf.force.names:

  Character vector of variable names.

- cont.cutoff:

  Integer. Cutoff for continuous.

## Value

Character vector of cut expressions.

## Examples

``` r
df <- data.frame(age = c(45, 60, 35), size = c(20, 35, 15), grade = c(1, 2, 3))
get_conf_force(df, conf.force.names = c("age", "size"))
#> NULL
```
