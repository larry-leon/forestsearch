# Tidy summary of selected subgroups across a comparison run

Reduces a `forestsearch_comparison` object to a tidy data.frame with one
row per combo, summarizing each configuration's selected subgroup.
Intended for post-hoc aggregation across analyses (e.g.\\ side-by-side
summary of a continuous-outcome comparison and a binary-outcome
comparison in a single document).

## Usage

``` r
extract_selected_subgroups(x, analysis_label = NULL)
```

## Arguments

- x:

  A `forestsearch_comparison` object as returned by
  [`compare_selection_rules`](https://larry-leon.github.io/forestsearch/reference/compare_selection_rules.md).

- analysis_label:

  Optional character scalar identifying the parent analysis (e.g.\\
  `"continuous"` or `"binary"`). When supplied, it appears as the first
  column. Useful when row-binding tidy summaries from multiple analyses.

## Value

A `data.frame` with columns:

- `analysis`:

  If `analysis_label` is supplied; otherwise omitted.

- `combo`:

  Integer index, `1..n_combos`.

- `combo_label`:

  Pretty label (e.g.\\ `"effMaxSG / pareto"`).

- `sg_focus`, `selection_rule`:

  The two comparison axes.

- `outcome_type`, `effect_measure`:

  Recorded per row from `fs$outcome_type` / `fs$effect_measure`, so
  mixed-effect-measure runs remain comparable.

- `subgroup`:

  Selected subgroup definition – `NA_character_` if no subgroup was
  identified; otherwise `paste(fs$sg.harm, collapse = " & ")`.

- `N_H`, `N_Hc`:

  Sample sizes in the harm subgroup H (`treat.recommend == 0`) and its
  complement Hc.

- `effect`:

  The winner's effect estimate on the natural scale (\\\exp\\ applied to
  the stored `hr` column for ratio-scale measures `OR`/`RR`/`IRR`;
  passed through otherwise).

- `Pcons`:

  The winner's consistency probability.

- `K`:

  Number of factors in the winner's subgroup definition.

- `n_passed`:

  Number of candidates that passed consistency (=
  `nrow(out_sg$result)`).

- `on_frontier`:

  Logical – did the selected subgroup lie on the Pareto frontier? `NA`
  when the frontier is unavailable.

- `error`:

  The error message for failed combos (`NA_character_` on success).

## Examples

``` r
if (FALSE) { # \dontrun{
cmp_cont <- compare_selection_rules(...)   # continuous analysis
cmp_bin  <- compare_selection_rules(...)   # binary analysis

tidy_cont <- extract_selected_subgroups(cmp_cont, analysis_label = "continuous")
tidy_bin  <- extract_selected_subgroups(cmp_bin,  analysis_label = "binary")

all_selected <- rbind(tidy_cont, tidy_bin)
saveRDS(all_selected, "selected_subgroups_all.rds")
} # }
```
