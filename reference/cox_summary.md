# Cox model summary for subgroup (OPTIMIZED)

Called in analyze_subgroup() \<– SG_tab_estimates

## Usage

``` r
cox_summary(
  Y,
  E,
  Treat,
  Strata = NULL,
  use_strata = !is.null(Strata),
  return_format = c("formatted", "numeric")
)
```

## Arguments

- Y:

  Numeric vector of outcome.

- E:

  Numeric vector of event indicators.

- Treat:

  Numeric vector of treatment indicators.

- Strata:

  Vector of strata (optional).

- use_strata:

  Logical. Whether to use strata in the model (default: TRUE if Strata
  provided).

- return_format:

  Character. "formatted" (default) or "numeric" for downstream use.

## Value

Character string with formatted HR and CI (or numeric vector if
return_format="numeric").

## Details

Calculates hazard ratio and confidence interval for a subgroup using Cox
regression. Optimized version with reduced overhead and better error
handling.

## Examples

``` r
# \donttest{
library(survival)
cox_summary(
  Y     = gbsg$rfstime / 30.4375,
  E     = gbsg$status,
  Treat = gbsg$hormon
)
#> [1] "0.69 (0.54, 0.89)"
# }
```
