# Apply a selection rule among passing candidates on one draw

For the size-within-band rules (`effMaxSG`/`effMinSG`) the inclusion
band is built to match the search's `selection_rule`:

- `"neighborhood"` - effect within `nbhd` of the max, on the NATURAL
  effect scale and multiplicative (`eff >= (1 - nbhd) * max(eff)`),
  matching
  [`sort_subgroups()`](https://larry-leon.github.io/forestsearch/reference/sort_subgroups.md)
  exactly (G3);

- `"pareto"` - the 2-D non-dominated set in (effect, size), via the same
  dominance core the search uses (`.pareto_dominated_xy()`) (G2);

- `"both"` - the intersection of the two. `beta` is on the working scale
  (log for ratio measures); `log_scale` controls the conversion to the
  natural effect for the band.

## Usage

``` r
.fs_mr_select(
  beta,
  zcons,
  sizes,
  passers,
  rule,
  nbhd,
  selection_rule = "neighborhood",
  log_scale = TRUE
)
```
