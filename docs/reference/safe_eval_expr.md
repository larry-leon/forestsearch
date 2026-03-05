# Evaluate an expression string in a data-frame scope

Parses and evaluates `expr` in a restricted environment containing only
the columns of `df` (parent:
[`baseenv()`](https://rdrr.io/r/base/environment.html)). This isolates
evaluation from the global environment, reducing scope for unintended
side effects.

## Usage

``` r
safe_eval_expr(df, expr)
```

## Arguments

- df:

  Data frame providing column names as variables.

- expr:

  Character. Expression to evaluate (e.g., `"BM > 1 & tmrsize > 19"`).

## Value

Result of evaluating `expr`, or `NULL` on failure.
