# Summary method for cox_ahr_cde objects

Summary method for cox_ahr_cde objects

## Usage

``` r
# S3 method for class 'cox_ahr_cde'
summary(object, ...)
```

## Arguments

- object:

  A `cox_ahr_cde` object from
  [`cox_ahr_cde_analysis`](https://larry-leon.github.io/forestsearch/reference/cox_ahr_cde_analysis.md).

- ...:

  Additional arguments (not used).

## Value

Invisibly returns the input object.

## Examples

``` r
if (FALSE) { # \dontrun{
library(survival)
df <- survival::gbsg
df$grade3 <- as.integer(df$grade == "3")
res <- cox_ahr_cde_analysis(df,
  outcome.name = "rfstime", event.name = "status", treat.name = "hormon",
  confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"))
summary(res)
} # }
```
