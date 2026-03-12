# Plot Spline Treatment Effect Function

Plot Spline Treatment Effect Function

## Usage

``` r
plot_spline_treatment_effect(dgm_result, add_points = TRUE)
```

## Arguments

- dgm_result:

  Result object from generate_aft_dgm_flex with spline

- add_points:

  Logical; add observed data points. Default TRUE

## Examples

``` r
if (FALSE) { # \dontrun{
library(survival)
df <- survival::gbsg
dgm <- generate_aft_dgm_flex(df, outcome.name = "rfstime",
                              event.name = "status", treat.name = "hormon",
                              confounders.name = c("age", "meno", "nodes"))
plot_spline_treatment_effect(dgm)
} # }
```
