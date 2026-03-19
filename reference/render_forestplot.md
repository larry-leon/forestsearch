# Render ForestSearch Forest Plot

Renders a forest plot from
[`plot_subgroup_results_forestplot()`](https://larry-leon.github.io/forestsearch/reference/plot_subgroup_results_forestplot.md).

## Usage

``` r
render_forestplot(x, newpage = TRUE)
```

## Arguments

- x:

  An fs_forestplot object from
  [`plot_subgroup_results_forestplot()`](https://larry-leon.github.io/forestsearch/reference/plot_subgroup_results_forestplot.md).

- newpage:

  Logical. Call grid.newpage() before drawing. Default: TRUE.

## Value

Invisibly returns the grob object.

## Details

To control plot sizing, create a custom theme using
[`create_forest_theme()`](https://larry-leon.github.io/forestsearch/reference/create_forest_theme.md)
and pass it to
[`plot_subgroup_results_forestplot()`](https://larry-leon.github.io/forestsearch/reference/plot_subgroup_results_forestplot.md):

`my_theme <- create_forest_theme(base_size = 14, row_padding = c(6, 4))`

`result <- plot_subgroup_results_forestplot(..., theme = my_theme)`

`render_forestplot(result)`
