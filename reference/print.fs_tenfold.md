# Print Method for Repeated K-Fold CV Results

Print Method for Repeated K-Fold CV Results

## Usage

``` r
# S3 method for class 'fs_tenfold'
print(x, ...)
```

## Arguments

- x:

  An fs_tenfold object

- ...:

  Additional arguments (ignored)

## Examples

``` r
if (FALSE) { # \dontrun{
fs <- forestsearch(gbsg,
  confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
  outcome.name = "rfstime", treat.name = "hormon", event.name = "status")
kfold <- forestsearch_tenfold(fs.est = fs, sims = 5)
print(kfold)
} # }
```
