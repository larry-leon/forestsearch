# Generate Detection Probability Curve

Computes detection probability across a range of hazard ratios to create
a power-like curve for subgroup detection.

## Usage

``` r
generate_detection_curve(
  theta_range = c(0.5, 3),
  n_points = 50L,
  n_sg,
  prop_cens = 0.3,
  hr_threshold = 1.25,
  hr_consistency = 1,
  include_reference = TRUE,
  method = "cubature",
  verbose = TRUE
)
```

## Arguments

- theta_range:

  Numeric vector of length 2. Range of HR values to evaluate. Default:
  c(0.5, 3.0)

- n_points:

  Integer. Number of points to evaluate. Default: 50

- n_sg:

  Integer. Subgroup sample size.

- prop_cens:

  Numeric. Proportion censored (0-1). Default: 0.3

- hr_threshold:

  Numeric. HR threshold for detection. Default: 1.25

- hr_consistency:

  Numeric. HR consistency threshold. Default: 1.0

- include_reference:

  Logical. Include reference HR values (0.5, 0.75, 1.0). Default: TRUE

- method:

  Character. Integration method. Default: "cubature"

- verbose:

  Logical. Print progress. Default: TRUE

## Value

A data.frame with columns:

- theta:

  Hazard ratio values

- probability:

  Detection probability

- n_sg:

  Subgroup size (repeated)

- prop_cens:

  Censoring proportion (repeated)

- hr_threshold:

  Detection threshold (repeated)

## Examples

``` r
if (FALSE) { # \dontrun{
# Generate detection curve
curve_data <- generate_detection_curve(
  n_sg = 60,
  prop_cens = 0.2,
  hr_threshold = 1.25
)

# Plot
plot_detection_curve(curve_data)
} # }
```
