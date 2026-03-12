# RMST calculation for subgroup

Calculates restricted mean survival time (RMST) for a subgroup.

## Usage

``` r
rmst_calculation(
  df,
  tte.name = "tte",
  event.name = "event",
  treat.name = "treat"
)
```

## Arguments

- df:

  Data frame.

- tte.name:

  Character. Name of time-to-event variable.

- event.name:

  Character. Name of event indicator variable.

- treat.name:

  Character. Name of treatment variable.

## Value

List with tau, RMST, RMST for treatment, RMST for control.

## Examples

``` r
if (FALSE) { # \dontrun{
library(survival)
df <- data.frame(
  tte   = gbsg$rfstime / 30.4375,
  event = gbsg$status,
  treat = gbsg$hormon
)
rmst_calculation(df, tte.name = "tte", event.name = "event",
                 treat.name = "treat")
} # }
```
