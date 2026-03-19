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
