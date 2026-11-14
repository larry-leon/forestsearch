# Plot Method for ForestSearch Forest Plot

Plot Method for ForestSearch Forest Plot

## Usage

``` r
# S3 method for class 'fs_forestplot'
plot(x, ...)
```

## Arguments

- x:

  An fs_forestplot object

- ...:

  Additional arguments (ignored)

## Value

The input `x`, invisibly.

## Examples

``` r
if (FALSE) { # \dontrun{
gbsg <- survival::gbsg
gbsg$grade3 <- as.integer(gbsg$grade == 3)
fs <- forestsearch(gbsg, id.name = "pid",
  confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
  outcome.name = "rfstime", treat.name = "hormon", event.name = "status")
fp <- plot_subgroup_results_forestplot(list(fs.est = fs), gbsg,
  outcome.name = "rfstime", event.name = "status", treat.name = "hormon")
plot(fp)
} # }
```
