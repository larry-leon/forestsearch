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

## Value

`x`, invisibly. No plot is drawn: a message names the
[`plot_sg_results()`](https://larry-leon.github.io/forestsearch/reference/plot_sg_results.md)
call that produces the requested plot. An error is raised when a
character `which` does not name an available plot.

## Examples

``` r
if (FALSE) { # \dontrun{
gbsg <- survival::gbsg
gbsg$grade3 <- as.integer(gbsg$grade == 3)
fs <- forestsearch(gbsg, id.name = "pid",
  confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
  outcome.name = "rfstime", treat.name = "hormon", event.name = "status")
p <- plot_sg_results(fs, outcome.name = "rfstime", event.name = "status",
  treat.name = "hormon")
plot(p)
} # }
```
