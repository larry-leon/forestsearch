# Plot Distribution of Identified Subgroups

Bar chart of subgroups identified across simulations, filtered to those
appearing in at least `min_pct` of the found simulations. Supports two
column schemas (MRCT and `run_simulation_analysis`), top-\\K\\ capping
with a pooled "Other" bar, automatic threshold halving when no label
clears the initial cut, and an `n_bars` attribute for adaptive
figure-height computation.

## Usage

``` r
plot_sg_distribution(
  results,
  min_pct = 5,
  title = "Distribution of Identified Subgroups",
  wrap_width = 25,
  any_col = "any_found",
  label_col = "sg_found",
  top_k = NULL,
  show_other = TRUE,
  min_floor_pct = 0.5,
  max_halvings = 0L,
  bar_label_inside = FALSE,
  placeholder_on_empty = FALSE
)
```

## Arguments

- results:

  `data.frame` or `data.table` with at least two columns: a
  per-replicate detection flag (0/1) and an identified-subgroup label
  string. Column names are configurable via `any_col` and `label_col`.
  The defaults (`"any_found"`, `"sg_found"`) match
  [`mrct_region_sims()`](https://larry-leon.github.io/forestsearch/reference/mrct_region_sims.md);
  pass `any_col = "any.H"`, `label_col = "sg.def"` to plot
  [`run_simulation_analysis()`](https://larry-leon.github.io/forestsearch/reference/run_simulation_analysis.md)
  output.

- min_pct:

  Numeric. Minimum percentage threshold for display (0-100). Default:
  `5`.

- title:

  Character. Plot title.

- wrap_width:

  Integer. Character width for wrapping long subgroup labels. Default:
  `25`.

- any_col:

  Character. Column name of the per-replicate detection flag. Default:
  `"any_found"`.

- label_col:

  Character. Column name of the identified-subgroup label string.
  Default: `"sg_found"`.

- top_k:

  Integer or `NULL`. Maximum number of bars to display; additional
  passing labels are pooled into a single "Other (\\k\\ labels)" bar
  (see `show_other`). Prevents unreadable charts under diffuse null
  hypotheses where finds are fragmented across many singleton labels.
  Default: `NULL` (no cap, i.e. the original behaviour).

- show_other:

  Logical. When `top_k` is set and extra labels are pooled, attach the
  "Other" bar. Default: `TRUE`. Set to `FALSE` to suppress the pooled
  bar entirely; its count and percent are still reported in the subtitle
  for honesty.

- min_floor_pct:

  Numeric. Lower bound on automatic threshold halving. Default: `0.5`.

- max_halvings:

  Integer. Maximum number of halvings applied to `min_pct` when no label
  clears the initial threshold. Default: `0` (no halving, i.e. the
  original behaviour). Set to `2` to match the extended-diagnostic
  behaviour used in the benefit-search vignettes.

- bar_label_inside:

  Logical. Render bar labels inside bars (`hjust = 1.05`) instead of
  outside (`hjust = -0.1`). Prevents small-count labels from clipping
  off the right edge. Default: `FALSE` (original outside-of-bar
  placement).

- placeholder_on_empty:

  Logical. When no subgroups clear the threshold (even after halving),
  return a minimal "no data" placeholder plot with
  `attr(p, "n_bars") = 0L` instead of `invisible(NULL)`. Default:
  `FALSE` (the original behaviour). Set to `TRUE` to simplify adaptive-
  figure-height rendering in vignettes.

## Value

A `ggplot2` object, or `invisible(NULL)` when no subgroups are found and
`placeholder_on_empty = FALSE`. When a plot is returned it carries two
attributes:

- `n_bars`:

  Integer. Number of bars in the plot (0 in the placeholder case). Use
  [`sgdist_fig_height()`](https://larry-leon.github.io/forestsearch/reference/sgdist_fig_height.md)
  to compute an adaptive figure height for knitr.

- `effective_min_pct`:

  Numeric. The threshold in effect after any automatic halving.

## See also

[`sgdist_fig_height`](https://larry-leon.github.io/forestsearch/reference/sgdist_fig_height.md)
for computing an adaptive figure height from `n_bars`, and
[`run_simulation_analysis`](https://larry-leon.github.io/forestsearch/reference/run_simulation_analysis.md)
/
[`mrct_region_sims`](https://larry-leon.github.io/forestsearch/reference/mrct_region_sims.md)
for the upstream simulation engines.

## Examples

``` r
if (FALSE) { # \dontrun{
# MRCT / mrct_region_sims() output (default column schema)
fs <- forestsearch(gbsg,
  confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
  outcome.name = "rfstime", treat.name = "hormon", event.name = "status")
plot_sg_distribution(fs$grp.consistency$result_new)

# run_simulation_analysis() output with top-K cap and auto-halving
fs_rows <- results_alt[results_alt$analysis == "FS", ]
p <- plot_sg_distribution(
  fs_rows,
  any_col              = "any.H",
  label_col            = "sg.def",
  top_k                = 15L,
  show_other           = FALSE,
  max_halvings         = 2L,
  bar_label_inside     = TRUE,
  placeholder_on_empty = TRUE,
  title                = "HTE: Subgroups Identified (H-hat)"
)
fig_h <- sgdist_fig_height(attr(p, "n_bars"))
} # }
```
