# Combined Pareto-Frontier Plot Across Configurations Sharing a Passing Set

For two or more
[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
fits whose consistency-passing candidate sets are **identical**,
produces a single Pareto-frontier plot annotated with one
`S<i>: <combo_label>` marker per configuration at each selected
subgroup. If two configurations pick the same subgroup, the markers
stack into a single multi-line label (e.g.\\
`"S1: hrMaxSG/pareto\nS2: hrMaxSG/both"`).

## Usage

``` r
plot_pareto_combined(
  fs_list,
  combo_labels = NULL,
  ci_table_list = NULL,
  show_band = TRUE,
  xlim = NULL,
  tolerance = 1e-06,
  verbose = TRUE
)
```

## Arguments

- fs_list:

  A list of `forestsearch` result objects. Length must be \>= 2.

- combo_labels:

  Character vector of length `length(fs_list)`. Used in the
  `S<i>: <label>` annotations. Default `NULL` = auto, e.g.\\
  `"config 1"`, `"config 2"`, ...

- ci_table_list:

  Optional list of CI tables from
  [`compute_frontier_cis`](https://larry-leon.github.io/forestsearch/reference/compute_frontier_cis.md),
  parallel to `fs_list`. When provided and the tables match (within
  tolerance), the CI bars are drawn using the first table.

- show_band:

  Logical. Draw the effect-band shading if applicable. Default `TRUE`.

- xlim:

  Numeric vector of length 2 or `NULL`. Default `NULL` = auto-range.

- tolerance:

  Numeric. Per-cell tolerance for the value equality check (`hr`, `N`,
  `E`, `K`). Default `1e-6`.

- verbose:

  Logical. Emit a warning naming the specific equality-check failure
  mode when the sets don't match. Default `TRUE`.

## Value

A `ggplot` object, or `NULL` (with a warning) if the passing sets don't
satisfy the equality criterion.

## Details

The combined plot is meaningful only when the passing sets match. When
they don't, the function returns `NULL` with a warning; the side-by-side
view (one panel per configuration) is the appropriate alternative.

## Equality criterion

Two passing sets are considered equal when, after sorting by subgroup
definition string:

- they have the same number of rows;

- the same set of subgroup definitions (concatenated `M.*` columns),
  compared as sets (order-independent);

- within each matched definition, the `hr`, `N`, `E`, and `K` values
  agree within `tolerance`.

`Pcons` is deliberately excluded from the value check. It can
legitimately drift across rules because the preview sort (which depends
on `selection_rule`) changes the queue order, which changes the
random-split state each candidate consumes. Drift of up to ~0.10 between
runs on the same subgroup is expected and does not indicate a real
disagreement. Note also that the internal candidate id `m` can differ
across configurations for the same reason; `m` is NOT used for equality.

## See also

[`compare_selection_rules`](https://larry-leon.github.io/forestsearch/reference/compare_selection_rules.md),
[`plot_pareto_frontier`](https://larry-leon.github.io/forestsearch/reference/plot_pareto_frontier.md),
[`pareto_frontier_table`](https://larry-leon.github.io/forestsearch/reference/pareto_frontier_table.md).

## Examples

``` r
if (FALSE) { # \dontrun{
out <- compare_selection_rules(
  df.analysis    = actg_df,
  sg_focus       = c("hrMaxSG", "hrMaxSG"),
  selection_rule = c("pareto",  "both"),
  ...
)

# Auto-attached by compare_selection_rules():
print(out$plot_combined)

# Or call directly:
p <- plot_pareto_combined(
  fs_list       = out$fs,
  combo_labels  = out$combos$label,
  ci_table_list = out$ci_tab
)
print(p)
} # }
```
