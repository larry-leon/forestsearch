# Plot Forest Plot of Hazard Ratios

Creates a forest plot showing hazard ratios with confidence intervals.

## Usage

``` r
plot_sg_forest(hr_estimates, sg0_name, sg1_name, colors, title = NULL, ...)
```

## Arguments

- hr_estimates:

  Data frame with HR estimates

- sg0_name:

  Character. Label for H subgroup

- sg1_name:

  Character. Label for Hc subgroup

- colors:

  List. Color specifications

- title:

  Character. Plot title

- ...:

  Additional arguments

## Value

Invisible NULL (creates plot as side effect)
