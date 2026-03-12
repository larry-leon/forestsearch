# Print method for survreg_comparison objects

Print method for survreg_comparison objects

## Usage

``` r
# S3 method for class 'multi_survreg_comparison'
print(x, ...)
```

## Arguments

- x:

  A survreg_comparison object

- ...:

  Additional arguments (not used)

## Value

Invisibly returns the input object

## Examples

``` r
if (FALSE) { # \dontrun{
library(survival)
df <- survival::gbsg
res <- compare_multiple_survreg(df, outcome.name = "rfstime",
                                 event.name = "status")
print(res)
} # }
```
