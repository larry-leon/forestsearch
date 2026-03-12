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

## Examples

``` r
if (FALSE) { # \dontrun{
fs <- forestsearch(gbsg,
  confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
  outcome.name = "rfstime", treat.name = "hormon", event.name = "status")
fp <- plot_subgroup_results_forestplot(list(fs.est = fs), gbsg,
  outcome.name = "rfstime", event.name = "status", treat.name = "hormon")
plot(fp)
} # }
```
