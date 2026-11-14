# Print Method for fs_sg_plot Objects

Print Method for fs_sg_plot Objects

## Usage

``` r
# S3 method for class 'fs_sg_plot'
print(x, ...)
```

## Arguments

- x:

  An fs_sg_plot object

- ...:

  Additional arguments (unused)

## Value

`x`, invisibly. Called for its side effect of printing the subgroup
definition, summary statistics and hazard ratio estimates.

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
print(p)
} # }
```
