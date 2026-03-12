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

## Note

`eval(parse())` is used intentionally here.
[`evaluate_comparison`](https://larry-leon.github.io/forestsearch/reference/evaluate_comparison.md)
handles only single comparisons (e.g., `"er <= 0"`); this function is
needed for the compound logical expressions produced by the ForestSearch
subgroup enumeration algorithm (e.g., `"er <= 0 & nodes > 3"`).
Evaluation is sandboxed: the environment contains only the columns of
`df` with [`baseenv()`](https://rdrr.io/r/base/environment.html) as
parent, so neither the global environment nor any package namespace is
in scope. No user-supplied strings are evaluated; only
internally-constructed subgroup definition strings reach this function.

## See also

[`evaluate_comparison`](https://larry-leon.github.io/forestsearch/reference/evaluate_comparison.md)
for the single-comparison operator-dispatch alternative that avoids
`eval(parse())`.
