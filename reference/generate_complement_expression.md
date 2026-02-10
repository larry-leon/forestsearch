# Generate Complement Expression

Creates the logical complement of a subgroup expression. Handles common
patterns like "var \<= x" -\> "var \> x".

## Usage

``` r
generate_complement_expression(expr)
```

## Arguments

- expr:

  Character vector of expressions to negate.

## Value

Character string with negated expression.
