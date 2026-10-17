# Plot the Pareto Frontier of Candidate Subgroups

Produces a 2D scatter of candidate subgroups in (effect, N) space, with
the Pareto frontier drawn as a step polyline and the selected subgroup
highlighted. When `ci_table` is supplied (the return value of
[`compute_frontier_cis`](https://larry-leon.github.io/forestsearch/reference/compute_frontier_cis.md)),
horizontal error bars for the split-derived 95\\

## Usage

``` r
plot_pareto_frontier(
  fs,
  ci_table = NULL,
  show_band = FALSE,
  effect_neighborhood = NULL,
  label_members = TRUE,
  point_size = 3,
  xlim = NULL,
  xlim_trim = NULL
)
```

## Arguments

- fs:

  A `forestsearch` result object.

- ci_table:

  Optional `data.table` from
  [`compute_frontier_cis`](https://larry-leon.github.io/forestsearch/reference/compute_frontier_cis.md).
  When supplied, horizontal error bars depict the split-derived 95\\
  Default `NULL` (no error bars).

- show_band:

  Logical. If `TRUE` and `sg_focus` is `"hrMaxSG"` or `"hrMinSG"`, shade
  the effect-size neighborhood band on the right of the plot. Default
  `FALSE`.

- effect_neighborhood:

  Numeric. Override of the band width; defaults to the value used in the
  original `forestsearch` call (or `0.10` if not recoverable).

- label_members:

  Logical. If `TRUE`, label each frontier member with its short cut
  description. Default `TRUE`.

- point_size:

  Numeric. Size of frontier-member points. Default `3`.

- xlim:

  Optional numeric vector of length 2 controlling the x-axis range,
  e.g.\\ `c(0.5, 3)`. By default the axis auto-expands to include
  split-CI bars when supplied, which can make individual points look
  crowded. Pass an explicit `xlim` (or use `xlim_trim = TRUE`) to zoom
  in.

- xlim_trim:

  Logical or `NULL`. If `TRUE` (and `xlim` is `NULL`), the x-axis is
  trimmed to the range of frontier point estimates with a small padding,
  ignoring CI bars. If `NULL` (default), auto-decides: `TRUE` when the
  effect measure is on a ratio (log) scale (`"OR"`, `"RR"`, `"IRR"`),
  where wide split-CI tails on the natural scale routinely dwarf the
  point-estimate range and force the axis to span orders of magnitude;
  `FALSE` otherwise. Pass `TRUE` or `FALSE` explicitly to override the
  auto behavior.

## Value

A `ggplot` object.

## Details

The frontier polyline is drawn as a downward step from large-effect /
small-N candidates to small-effect / large-N candidates. Points off the
frontier (if any are present in the `fs$grp.consistency$out_sg$result`
table beyond those on the frontier) are not currently plotted; this
function only displays the frontier and its selected member.

For `sg_focus = "hrMinSG"`, the selected subgroup may not lie on the
Pareto frontier (the focus deliberately prefers small subgroups, which
are typically N-dominated). In that case the selected marker appears off
the polyline.

## See also

[`pareto_frontier_table`](https://larry-leon.github.io/forestsearch/reference/pareto_frontier_table.md),
[`compute_frontier_cis`](https://larry-leon.github.io/forestsearch/reference/compute_frontier_cis.md),
[`frontier_member_flags`](https://larry-leon.github.io/forestsearch/reference/frontier_member_flags.md).

## Examples

``` r
if (FALSE) { # \dontrun{
p <- plot_pareto_frontier(fs)
print(p)

# With split-derived 95% CIs:
ci_tab <- compute_frontier_cis(fs, n_splits = 1000)
p2 <- plot_pareto_frontier(fs, ci_table = ci_tab)
print(p2)
} # }
```
