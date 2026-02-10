# Create Sample Size Table for Multiple Scenarios

Generates a table of required sample sizes for different combinations of
true hazard ratios and censoring proportions.

## Usage

``` r
create_sample_size_table(
  theta_values,
  prop_cens_values,
  target_power = 0.8,
  hr_threshold = 1.25,
  verbose = TRUE
)
```

## Arguments

- theta_values:

  Numeric vector. True hazard ratios to evaluate.

- prop_cens_values:

  Numeric vector. Censoring proportions to evaluate.

- target_power:

  Numeric. Target detection probability. Default: 0.80

- hr_threshold:

  Numeric. HR threshold. Default: 1.25

- verbose:

  Logical. Print progress. Default: TRUE

## Value

A data.frame with columns: theta, prop_cens, n_required, achieved_power

## Examples

``` r
if (FALSE) { # \dontrun{
ss_table <- create_sample_size_table(
  theta_values = c(1.5, 1.75, 2.0, 2.5),
  prop_cens_values = c(0.2, 0.3, 0.4),
  target_power = 0.80
)
print(ss_table)
} # }
```
