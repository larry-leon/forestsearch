# Plot Kaplan-Meier Survival Curves for Subgroups

Creates side-by-side Kaplan-Meier survival curves for the H and Hc
subgroups.

## Usage

``` r
plot_sg_km(
  df_H,
  df_Hc,
  outcome.name,
  event.name,
  treat.name,
  by.risk,
  sg0_name,
  sg1_name,
  treat_labels,
  colors,
  show_ci = TRUE,
  show_logrank = TRUE,
  show_hr = TRUE,
  hr_estimates = NULL,
  conf.level = 0.95,
  title = NULL,
  ...
)
```

## Arguments

- df_H:

  Data frame for H subgroup

- df_Hc:

  Data frame for Hc subgroup

- outcome.name:

  Character. Outcome variable name

- event.name:

  Character. Event indicator name

- treat.name:

  Character. Treatment variable name

- by.risk:

  Numeric. Risk table interval

- sg0_name:

  Character. Label for H subgroup

- sg1_name:

  Character. Label for Hc subgroup

- treat_labels:

  Named character vector. Treatment labels

- colors:

  List. Color specifications

- show_ci:

  Logical. Show confidence intervals

- show_logrank:

  Logical. Show log-rank p-value

- show_hr:

  Logical. Show HR annotation

- hr_estimates:

  Data frame with HR estimates

- conf.level:

  Numeric. Confidence level

- title:

  Character. Plot title

- ...:

  Additional arguments

## Value

Invisible NULL (creates plot as side effect)
