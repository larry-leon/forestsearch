# Recommended figure height for a forest panel

Converts a panel's row count to the `fig.height` (inches) used by the
vignette chunks: `rows * row_height + overhead`, rounded to 0.1.

## Usage

``` r
forest_height(
  x,
  panel = c("single", "combo", "highrisk"),
  row_height = 0.45,
  overhead = 1.5
)
```

## Arguments

- x:

  A `"subgroup_sims_summary"` object.

- panel:

  `"single"`, `"combo"`, or `"highrisk"`.

- row_height:

  Inches per row (default 0.45).

- overhead:

  Fixed inches for title, x-axis, and footnote (default 1.5).

## Value

A single numeric, rounded to one decimal.

## See also

[`plot.subgroup_sims_summary()`](https://larry-leon.github.io/forestsearch/reference/plot.subgroup_sims_summary.md)

## Examples

``` r
if (FALSE) { # \dontrun{
fh_single <- forest_height(S, "single")
} # }
```
