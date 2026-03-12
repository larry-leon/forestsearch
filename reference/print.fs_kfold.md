# Print Method for K-Fold CV Results

Print Method for K-Fold CV Results

## Usage

``` r
# S3 method for class 'fs_kfold'
print(x, ...)
```

## Arguments

- x:

  An fs_kfold object

- ...:

  Additional arguments (ignored)

## Examples

``` r
if (FALSE) { # \dontrun{
fs <- forestsearch(gbsg,
  confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
  outcome.name = "rfstime", treat.name = "hormon", event.name = "status")
oob <- forestsearch_Kfold(fs.est = fs)
print(oob)
} # }
```
