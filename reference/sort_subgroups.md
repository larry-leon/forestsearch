# Sort Subgroups by Focus

Sorts a data.table of subgroup results according to the specified focus.

## Usage

``` r
sort_subgroups(result_new, sg_focus)
```

## Arguments

- result_new:

  A data.table of subgroup results.

- sg_focus:

  Sorting focus: "hr", "hrMaxSG", "maxSG", "hrMinSG", "minSG".

## Value

A sorted data.table.

## Examples

``` r
# \donttest{
library(data.table)
dt <- data.table(Pcons = c(0.92, 0.95, 0.88),
                 hr    = c(1.8, 1.5, 2.1),
                 N     = c(80, 120, 60),
                 K     = c(1, 2, 1))
sort_subgroups(dt, sg_focus = "hr")
#>    Pcons    hr     N     K
#>    <num> <num> <num> <num>
#> 1:  0.95   1.5   120     2
#> 2:  0.92   1.8    80     1
#> 3:  0.88   2.1    60     1
# }
```
