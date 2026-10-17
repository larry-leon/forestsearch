# Process forced cut expression for a variable

Evaluates a cut expression (e.g., "age \<= mean(age)") and returns the
expression with the value.

## Usage

``` r
process_conf_force_expr(expr, df)
```

## Arguments

- expr:

  Character string of the cut expression.

- df:

  Data frame.

## Value

Character string with evaluated value.

## Examples

``` r
df <- data.frame(age = c(40, 55, 70, 35), size = c(20, 30, 25, 15))
process_conf_force_expr("age <= mean(age)", df)
#> [1] "age <= 50"
process_conf_force_expr("size <= median(size)", df)
#> [1] "size <= 22.5"
```
