# Output Subgroup Consistency Results

Returns the top subgroup(s) and recommended treatment flags.

## Usage

``` r
sg_consistency_out(
  df,
  result_new,
  sg_focus,
  index.Z,
  names.Z,
  details = FALSE,
  plot.sg = FALSE,
  by.risk = 12,
  confs_labels,
  is_glm = FALSE,
  selection_rule = "neighborhood",
  effect_neighborhood = 0.1,
  effect_log_scale = FALSE
)
```

## Arguments

- df:

  Data.frame. Original analysis data.

- result_new:

  Data.table. Sorted subgroup results.

- sg_focus:

  Character. Sorting focus criterion.

- index.Z:

  Matrix. Subgroup factor indicators.

- names.Z:

  Character vector. Factor column names.

- details:

  Logical. Print details.

- plot.sg:

  Logical. Plot subgroup curves.

- by.risk:

  Numeric. Risk interval for plotting.

- confs_labels:

  Character vector. Human-readable labels.

- is_glm:

  Logical. If `TRUE`, suppresses Kaplan-Meier survival plots (which are
  not meaningful for GLM outcomes). Default `FALSE`.

- selection_rule:

  Character. Rule defining the candidate inclusion set for `"hrMaxSG"` /
  `"hrMinSG"`. See
  [`sort_subgroups`](https://larry-leon.github.io/forestsearch/reference/sort_subgroups.md).
  Default `"neighborhood"`.

- effect_neighborhood:

  Numeric in `[0, 1)`. Relative tolerance for `"hrMaxSG"`/`"hrMinSG"`
  selection. See
  [`sort_subgroups`](https://larry-leon.github.io/forestsearch/reference/sort_subgroups.md).
  Default `0.10`.

- effect_log_scale:

  Logical. If `TRUE`, the `hr` column stores log-scale values and is
  exponentiated for the neighborhood test and for Pareto-frontier
  dominance. Default `FALSE`.

## Value

List with elements:

- result:

  Sorted candidate table (top row = selected subgroup).

- pareto_frontier:

  Data.table of non-dominated candidates on (effect, N), both maximized.
  See
  [`compute_pareto_frontier`](https://larry-leon.github.io/forestsearch/reference/compute_pareto_frontier.md).
  May be `NULL` if computation failed.

- sg.harm:

  Factor-level cut names defining the selected subgroup.

- sg.harm_label:

  Human-readable subgroup labels.

- df_flag:

  Per-subject treatment-recommendation flags.

- sg.harm.id:

  Per-subject subgroup-membership indicator.

## Examples

``` r
if (FALSE) { # \dontrun{
# sg_consistency_out is called internally by forestsearch().
# See forestsearch() for the standard entry point.
fs <- forestsearch(gbsg,
  confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
  outcome.name = "rfstime", treat.name = "hormon", event.name = "status")
fs$grp.consistency$out_sg
} # }
```
