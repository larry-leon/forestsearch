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

## Examples

``` r
if (FALSE) { # \dontrun{
fs <- forestsearch(gbsg,
  confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
  outcome.name = "rfstime", treat.name = "hormon", event.name = "status")
p <- plot_subgroup(fs)
print(p)
} # }
```
