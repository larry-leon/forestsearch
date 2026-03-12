# Subgroup summary table estimates

Returns a summary table of subgroup estimates (HR, RMST, medians, etc.).

## Usage

``` r
SG_tab_estimates(
  df,
  SG_flag,
  outcome.name = "tte",
  event.name = "event",
  treat.name = "treat",
  strata.name = NULL,
  hr_1a = NA,
  hr_0a = NA,
  potentialOutcome.name = NULL,
  sg1_name = NULL,
  sg0_name = NULL,
  draws = 0,
  details = FALSE,
  return_medians = TRUE,
  est.scale = "hr"
)
```

## Arguments

- df:

  Data frame.

- SG_flag:

  Character. Subgroup flag variable.

- outcome.name:

  Character. Name of outcome variable.

- event.name:

  Character. Name of event indicator variable.

- treat.name:

  Character. Name of treatment variable.

- strata.name:

  Character. Name of strata variable (optional).

- hr_1a:

  Character. Adjusted HR for subgroup 1 (optional).

- hr_0a:

  Character. Adjusted HR for subgroup 0 (optional).

- potentialOutcome.name:

  Character. Name of potential outcome variable (optional).

- sg1_name:

  Character. Name for subgroup 1.

- sg0_name:

  Character. Name for subgroup 0.

- draws:

  Integer. Number of draws for resampling (optional).

- details:

  Logical. Print details.

- return_medians:

  Logical. Use medians or RMST.

- est.scale:

  Character. Effect scale ("hr" or "1/hr").

## Value

Data frame of subgroup summary estimates.

## Examples

``` r
if (FALSE) { # \dontrun{
library(survival)
df <- data.frame(
  tte          = gbsg$rfstime / 30.4375,
  event        = gbsg$status,
  treat        = gbsg$hormon,
  treat.recommend = as.integer(gbsg$er > 0)
)
SG_tab_estimates(df, SG_flag = "ITT",
                 outcome.name = "tte",
                 event.name   = "event",
                 treat.name   = "treat")
} # }
```
