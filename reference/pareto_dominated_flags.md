# Pareto Dominance Flags

Internal helper. For each row of `result_dt`, returns `TRUE` if the row
is dominated in (`hr`, `N`) space (both maximized) by another row, and
`FALSE` otherwise. Rows with `NA` `hr` or `N` are flagged as dominated.

## Usage

``` r
pareto_dominated_flags(result_dt, effect_log_scale = FALSE)
```

## Arguments

- result_dt:

  Data.table of candidate subgroups with columns `hr` and `N`.

- effect_log_scale:

  Logical. If `TRUE`, `hr` is exponentiated before comparison.

## Value

Logical vector of length `nrow(result_dt)`.

## Details

Used by both
[`compute_pareto_frontier()`](https://larry-leon.github.io/forestsearch/reference/compute_pareto_frontier.md)
(which filters by the negation of this vector) and the
`selection_rule = "pareto"` branches of
[`sort_subgroups()`](https://larry-leon.github.io/forestsearch/reference/sort_subgroups.md)
/
[`sort_subgroups_preview()`](https://larry-leon.github.io/forestsearch/reference/sort_subgroups_preview.md)
(which use it as an inclusion indicator).
