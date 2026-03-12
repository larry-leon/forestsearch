# Plot Method for fs_sg_plot Objects

Plot Method for fs_sg_plot Objects

## Usage

``` r
# S3 method for class 'fs_sg_plot'
plot(x, which = 1, ...)
```

## Arguments

- x:

  An fs_sg_plot object

- which:

  Character or integer. Which plot to display. Default: 1 (first
  available)

- ...:

  Additional arguments passed to plot functions

## Examples

``` r
if (FALSE) { # \dontrun{
fs <- forestsearch(gbsg,
  confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
  outcome.name = "rfstime", treat.name = "hormon", event.name = "status")
p <- plot_subgroup(fs)
plot(p)
} # }
```
