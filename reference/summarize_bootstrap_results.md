# Enhanced Bootstrap Results Summary

Creates comprehensive output including formatted table with subgroup
footnote, diagnostic plots, bootstrap quality metrics, and detailed
timing analysis.

## Usage

``` r
summarize_bootstrap_results(
  sgharm,
  boot_results,
  create_plots = FALSE,
  est.scale = "hr"
)
```

## Arguments

- sgharm:

  The selected subgroup object from forestsearch results. Can be:

  - Character vector of factor definitions (e.g., c("{age\>=50}",
    "{nodes\>=3}"))

  - List with `sgharm` element containing factor definitions

  - List with `sg.harm_label` element (human-readable labels)

- boot_results:

  List. Output from forestsearch_bootstrap_dofuture()

- create_plots:

  Logical. Generate diagnostic plots (default: FALSE)

- est.scale:

  Character. "hr" or "1/hr" for effect scale

## Value

List with components:

- table:

  gt table with treatment effects and subgroup footnote

- diagnostics:

  List of bootstrap quality metrics

- diagnostics_table_gt:

  gt table of diagnostics

- plots:

  List of ggplot2 diagnostic plots (if create_plots=TRUE)

- timing:

  List of timing analysis (if timing data available)

- subgroup_summary:

  List from summarize_bootstrap_subgroups()

## Details

The `table` output includes a footnote displaying the identified
subgroup definition, analogous to the `tab_estimates` table from
[`sg_tables`](https://larry-leon.github.io/forestsearch/reference/sg_tables.md).
This is achieved by extracting the subgroup definition from `sgharm` and
passing it to
[`format_bootstrap_table`](https://larry-leon.github.io/forestsearch/reference/format_bootstrap_table.md).

## See also

[`format_bootstrap_table`](https://larry-leon.github.io/forestsearch/reference/format_bootstrap_table.md)
for table creation
[`sg_tables`](https://larry-leon.github.io/forestsearch/reference/sg_tables.md)
for analogous main analysis tables
[`summarize_bootstrap_subgroups`](https://larry-leon.github.io/forestsearch/reference/summarize_bootstrap_subgroups.md)
for subgroup stability analysis
