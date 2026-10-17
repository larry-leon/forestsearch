# Generate J-quantile cut expressions for a continuous variable

For a continuous variable, returns `J` cut expressions of the form
`"X <= qj(X, k, J + 1)"` for `k = 1, ..., J`. These `J` cut points are
placed at the (k/(J+1))-th empirical quantiles of X and partition its
range into `J + 1` non-overlapping intervals \$\$\[\min(X), c_1),\\
\[c_1, c_2),\\ \ldots,\\ \[c_J, \max(X)\]\$\$ where \\c_k\\ is the
(k/(J+1))-th quantile of X. Note the cut expressions themselves are
nested half-spaces (each is a subset of the next); ForestSearch consumes
them as binary candidate factors and combines them via intersection
during the search.

## Usage

``` r
cut_var_jq(x, J)
```

## Arguments

- x:

  Character. Variable name.

- J:

  Integer \>= 1. Number of binary cut expressions to emit. The resulting
  partition has `J + 1` non-overlapping intervals. Cuts are placed at
  the (k/(J+1))-th empirical quantiles for `k = 1, ..., J` (no cut at
  the boundaries).

## Value

Character vector of `J` cut expressions.

## Details

Expressions are emitted in deferred form (with literal `qj(...)` calls
inside the string) so that they are correctly recomputed when the same
expression is processed against a different data subset (e.g., a
bootstrap replicate). They are subsequently resolved to literal numerics
by
[`process_conf_force_expr()`](https://larry-leon.github.io/forestsearch/reference/process_conf_force_expr.md).
