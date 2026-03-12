# Format Bootstrap Results Table with gt

Creates a publication-ready table from ForestSearch bootstrap results,
with bias-corrected confidence intervals, informative formatting, and
optional subgroup definition footnote.

## Usage

``` r
format_bootstrap_table(
  FSsg_tab,
  nb_boots,
  est.scale = "hr",
  boot_success_rate = NULL,
  sg_definition = NULL,
  title = NULL,
  subtitle = NULL
)
```

## Arguments

- FSsg_tab:

  Data frame or matrix from forestsearch_bootstrap_dofuture()\$FSsg_tab

- nb_boots:

  Integer. Number of bootstrap iterations performed

- est.scale:

  Character. "hr" or "1/hr" for effect scale

- boot_success_rate:

  Numeric. Proportion of bootstraps that found subgroups

- sg_definition:

  Character. Subgroup definition string to display as footnote (e.g.,
  "{age\>=50} & {nodes\>=3}"). If NULL, no subgroup footnote is added.

- title:

  Character. Custom title (optional)

- subtitle:

  Character. Custom subtitle (optional)

## Value

A gt table object

## Examples

``` r
if (FALSE) { # \dontrun{
fs <- forestsearch(gbsg,
  confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
  outcome.name = "rfstime", treat.name = "hormon", event.name = "status")
fs_bc <- forestsearch_bootstrap_dofuture(fs, nb_boots = 100)
summaries <- summarize_bootstrap_results(fs$sg.harm, fs_bc)
format_bootstrap_table(summaries$table_data)
} # }
```
