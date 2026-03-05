# Create Bootstrap Diagnostic Plots

Generates diagnostic visualization plots for bootstrap analysis.

## Usage

``` r
create_bootstrap_diagnostic_plots(
  results,
  H_estimates,
  Hc_estimates,
  overall_timing = NULL
)
```

## Arguments

- results:

  Data frame with bootstrap results

- H_estimates:

  List with H subgroup estimates

- Hc_estimates:

  List with Hc subgroup estimates

- overall_timing:

  List with overall timing information (optional)

## Value

List of ggplot2 objects
