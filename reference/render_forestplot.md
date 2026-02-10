# Render ForestSearch Forest Plot

Renders a forest plot from
[`plot_subgroup_results_forestplot()`](https://github.com/larry-leon/forestsearch/reference/plot_subgroup_results_forestplot.md).

## Usage

``` r
render_forestplot(x, newpage = TRUE)
```

## Arguments

- x:

  An fs_forestplot object from
  [`plot_subgroup_results_forestplot()`](https://github.com/larry-leon/forestsearch/reference/plot_subgroup_results_forestplot.md).

- newpage:

  Logical. Call grid.newpage() before drawing. Default: TRUE.

## Value

Invisibly returns the grob object.

## Details

To control plot sizing, create a custom theme using
[`create_forest_theme()`](https://github.com/larry-leon/forestsearch/reference/create_forest_theme.md)
and pass it to
[`plot_subgroup_results_forestplot()`](https://github.com/larry-leon/forestsearch/reference/plot_subgroup_results_forestplot.md):

`my_theme <- create_forest_theme(base_size = 14, row_padding = c(6, 4))`

`result <- plot_subgroup_results_forestplot(..., theme = my_theme)`

`render_forestplot(result)`

## Examples

``` r
if (FALSE) { # \dontrun{
# For larger plot, use theme parameter in plot_subgroup_results_forestplot:
large_theme <- create_forest_theme(
  base_size = 14,
  row_padding = c(6, 4),
  cv_fontsize = 11
)

result <- plot_subgroup_results_forestplot(
  fs_results = list(fs.est = fs, fs_bc = fs_bc),
  df_analysis = df.analysis,
  outcome.name = "time",
  event.name = "status",
  treat.name = "treatment",
  theme = large_theme
)

render_forestplot(result)
} # }
```
