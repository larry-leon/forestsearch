# Identification structure across a sweep

Applies
[`fs_identification_structure()`](https://larry-leon.github.io/forestsearch/reference/fs_identification_structure.md)
across a named collection of bundles – typically one per trial size –
and returns a single object the plot and table helpers consume. This is
the package-level replacement for the inline
`fs_identification_figures*.qmd` fragments.

## Usage

``` r
fs_identification_summary(
  bundles,
  anchor,
  partner,
  proxy,
  confounders = NULL,
  sweep_label = "Trial size (n)",
  sg_def_col = "sg_def",
  covs_col = "covs"
)
```

## Arguments

- bundles:

  A named list. Each element is a results data frame, a bundle carrying
  `results`, or a path to an `.rds`. Names label the sweep points (for
  example `"500"` and `"1000"`) and are used on the x axis in the order
  supplied.

- anchor, partner, proxy:

  Character. The planted-structure triple.

- confounders:

  Character vector. The analysis covariate pool.

- sweep_label:

  Character. Axis label for the sweep dimension.

- sg_def_col, covs_col:

  Character. Passed through.

## Value

An object of class `fs_identification`: a list with `per_point` (the raw
aggregates), tidy `involvement` and `structure` data frames, an
`accuracy` data frame, and the triple it was built with.

## Examples

``` r
if (FALSE) { # \dontrun{
id <- fs_identification_summary(
  list(`500` = "results/fs_..._n500_combined_1_500.rds",
       `1000` = "results/fs_..._n1000_combined_1_500.rds"),
  anchor = "er", partner = "meno", proxy = "age",
  confounders = c("er", "age", "meno", "pgr", "nodes", "size", "grade"))
plot_fs_identification_involvement(id)
plot_fs_identification_structure(id)
fs_identification_table(id)
} # }
```
