# Bias-vs-coverage display (three panels)

Draws the standard bias-vs-coverage display from a stacked
[`fs_sim_bias_coverage()`](https://larry-leon.github.io/forestsearch/reference/fs_sim_bias_coverage.md)
table carrying a `cell` column: (1) one-sided coverage against `b` with
Gaussian-reference curves `Phi(z1 * r - b)` for the given `curves`
values of `r`; (2) two-sided coverage against `b` with curves
`Phi(z2 * r - b) - Phi(-z2 * r - b)`; (3) observed coverage against the
Gaussian reference at each point's own `(b, r)`, with the diagonal.
Estimators are distinguished by marker shape, the nominal level by a
dashed line, and cells by text labels.

## Usage

``` r
fs_plot_bias_coverage(
  tbl,
  labels = TRUE,
  curves = c(1, 1.25, 1.5, 2),
  level = 0.95,
  side = c("lower", "upper")
)
```

## Arguments

- tbl:

  Stacked output of
  [`fs_sim_bias_coverage()`](https://larry-leon.github.io/forestsearch/reference/fs_sim_bias_coverage.md)
  with an added `cell` column (one factor level per simulation cell).

- labels:

  Draw cell labels next to the points (default `TRUE`).

- curves:

  Gaussian-reference `r` values for panels 1-2 (default
  `c(1, 1.25, 1.5, 2)`).

- level:

  Nominal level for the reference lines (default 0.95).

- side:

  `"lower"` (default) or `"upper"`: which one-sided bound the table's
  `cov1` was computed for (pass the same value given to
  [`fs_sim_bias_coverage()`](https://larry-leon.github.io/forestsearch/reference/fs_sim_bias_coverage.md)).
  Only panel 1's reference curves and title change: `Phi(z1 * r + b)`
  under `"upper"`.

## Value

A patchwork object of the three ggplot panels.
