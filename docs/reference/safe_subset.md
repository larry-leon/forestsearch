# Subset a data frame using an expression string

Thin wrapper around
[`safe_eval_expr`](https://larry-leon.github.io/forestsearch/reference/safe_eval_expr.md)
that uses the logical result to subset rows.

## Usage

``` r
safe_subset(df, expr)
```

## Arguments

- df:

  Data frame.

- expr:

  Character. Subset expression (e.g., `"BM > 1 & tmrsize > 19"`).

## Value

Subset of `df`, or `NULL` on failure.
