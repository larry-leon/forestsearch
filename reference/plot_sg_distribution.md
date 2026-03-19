# Plot Distribution of Identified Subgroups

Bar chart of subgroups identified across simulations, filtered to those
appearing in at least `min_pct` of the found simulations.

## Usage

``` r
plot_sg_distribution(
  results,
  min_pct = 5,
  title = "Distribution of Identified Subgroups",
  wrap_width = 25
)
```

## Arguments

- results:

  data.table from
  [`mrct_region_sims()`](https://larry-leon.github.io/forestsearch/reference/mrct_region_sims.md)

- min_pct:

  Numeric. Minimum percentage threshold for display (default: 5)

- title:

  Character. Plot title. Default: "Distribution of Identified Subgroups"

- wrap_width:

  Integer. Character width for wrapping long subgroup labels. Default:
  25

## Value

A ggplot2 object
