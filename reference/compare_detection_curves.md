# Compare Detection Curves Across Sample Sizes

Generates and compares detection probability curves for multiple
subgroup sample sizes.

## Usage

``` r
compare_detection_curves(
  n_sg_values,
  prop_cens = 0.3,
  hr_threshold = 1.25,
  hr_consistency = 1,
  theta_range = c(0.5, 3),
  n_points = 40L,
  verbose = TRUE
)
```

## Arguments

- n_sg_values:

  Integer vector. Subgroup sample sizes to compare.

- prop_cens:

  Numeric. Proportion censored. Default: 0.3

- hr_threshold:

  Numeric. HR threshold. Default: 1.25

- hr_consistency:

  Numeric. HR consistency threshold. Default: 1.0

- theta_range:

  Numeric vector of length 2. Range of HR values. Default: c(0.5, 3.0)

- n_points:

  Integer. Number of points per curve. Default: 40

- verbose:

  Logical. Print progress. Default: TRUE

## Value

A data.frame with all curves combined, including n_sg as a factor.

## Examples

``` r
if (FALSE) { # \dontrun{
comparison <- compare_detection_curves(
  n_sg_values = c(40, 60, 80, 100),
  prop_cens = 0.2
)

# Plot with ggplot2
library(ggplot2)
ggplot(comparison, aes(x = theta, y = probability, color = factor(n_sg))) +
  geom_line(linewidth = 1) +
  labs(x = "True HR", y = "P(Detect)", color = "n_sg") +
  theme_minimal()
} # }
```
