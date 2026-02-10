# Plot Detection Probability Curve

Creates a visualization of the detection probability curve.

## Usage

``` r
plot_detection_curve(
  curve_data,
  add_reference_lines = TRUE,
  add_threshold_line = TRUE,
  title = NULL,
  ...
)
```

## Arguments

- curve_data:

  A data.frame from
  [`generate_detection_curve`](https://github.com/larry-leon/forestsearch/reference/generate_detection_curve.md)
  or with columns: theta, probability.

- add_reference_lines:

  Logical. Add horizontal reference lines at 0.05, 0.10, 0.80. Default:
  TRUE

- add_threshold_line:

  Logical. Add vertical line at hr_threshold. Default: TRUE

- title:

  Character. Plot title. Default: auto-generated

- ...:

  Additional arguments passed to plot()

## Value

Invisibly returns the input data.

## Examples

``` r
if (FALSE) { # \dontrun{
curve_data <- generate_detection_curve(n_sg = 60, prop_cens = 0.2)
plot_detection_curve(curve_data)
} # }
```
