# Adaptive Figure Height for a Subgroup-Distribution Plot

Computes a reasonable knitr `fig.height` (in inches) for a
[`plot_sg_distribution()`](https://larry-leon.github.io/forestsearch/reference/plot_sg_distribution.md)
panel based on the number of bars, keeping the visual quality consistent
across scenarios with very different numbers of distinct subgroup
labels.

## Usage

``` r
sgdist_fig_height(
  n_bars,
  floor_in = 4,
  ceil_in = 10,
  per_bar = 0.35,
  base_in = 3
)
```

## Arguments

- n_bars:

  Integer. Number of bars in the plot, typically retrieved via
  `attr(p, "n_bars")` on a plot returned by
  [`plot_sg_distribution()`](https://larry-leon.github.io/forestsearch/reference/plot_sg_distribution.md).

- floor_in:

  Numeric. Minimum height in inches. Default: `4.0`.

- ceil_in:

  Numeric. Maximum height in inches. Default: `10.0`.

- per_bar:

  Numeric. Height contribution per bar in inches. Default: `0.35`.

- base_in:

  Numeric. Base height added to `per_bar * n_bars` before applying the
  floor / ceiling. Default: `3.0`.

## Value

Numeric scalar, a height in inches in `[floor_in, ceil_in]`. Suitable
for the `fig.height` chunk option.

## See also

[`plot_sg_distribution`](https://larry-leon.github.io/forestsearch/reference/plot_sg_distribution.md)

## Examples

``` r
if (FALSE) { # \dontrun{
p <- plot_sg_distribution(results, placeholder_on_empty = TRUE)
h <- sgdist_fig_height(attr(p, "n_bars"))
# Use `h` as the fig.height chunk option.
} # }
```
