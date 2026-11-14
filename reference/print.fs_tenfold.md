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

## Value

`x`, invisibly. Called for its side effect of printing the repeated
K-fold summary.

## Examples

``` r
if (FALSE) { # \dontrun{
gbsg <- survival::gbsg
gbsg$grade3 <- as.integer(gbsg$grade == 3)
fs <- forestsearch(gbsg, id.name = "pid",
  confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
  outcome.name = "rfstime", treat.name = "hormon", event.name = "status")
kfold <- forestsearch_tenfold(fs.est = fs, sims = 5)
print(kfold)
} # }
```
