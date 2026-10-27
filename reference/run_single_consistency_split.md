# Run Single Consistency Split

Performs one random 50/50 split and evaluates whether both halves meet
the HR consistency threshold.

## Usage

``` r
run_single_consistency_split(
  df.x,
  N.x,
  hr.consistency,
  cox_init = 0,
  estimator_fn = NULL,
  consistency_threshold = NULL,
  adjust_covariates = NULL
)
```

## Arguments

- df.x:

  data.table. Subgroup data with columns Y, Event, Treat.

- N.x:

  Integer. Number of observations in subgroup.

- hr.consistency:

  Numeric. `c2`, **an effect threshold** – the minimum effect required
  in *each* half of a split for that split to count as consistent.
  **Not** a proportion and **not** the consistency rate.
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  calls this `consistency.threshold` / `hr.consistency`.

- cox_init:

  Numeric. Initial value for Cox model (log HR).

- estimator_fn:

  Closure or `NULL`. Effect-estimator closure from
  [`make_effect_estimator`](https://larry-leon.github.io/forestsearch/reference/make_effect_estimator.md)
  for GLM outcomes. Default `NULL` (uses Cox model).

- consistency_threshold:

  Numeric or `NULL`. Threshold for GLM consistency evaluation. Default
  `NULL` (uses `hr.consistency`).

- adjust_covariates:

  This option is implemented but has not been tested against the
  manuscript analyses or the simulation studies of León and Anderson
  (2026a, 2026b), which use the unadjusted treatment-only fits.
  Character vector or `NULL`. Additional Cox model terms (e.g.
  `"strata(x1)"`) passed to
  [`get_split_hr_fast`](https://larry-leon.github.io/forestsearch/reference/get_split_hr_fast.md)
  on the survival path. Ignored on the GLM path. Default `NULL`
  (unadjusted).

## Value

Numeric. 1 if both splits meet threshold, 0 if not, NA if error.

## Examples

``` r
# \donttest{
library(data.table)
set.seed(1)
df <- data.table(
  Y     = rexp(100),
  Event = rbinom(100, 1, 0.55),
  Treat = rep(0:1, 50)
)
run_single_consistency_split(df, N.x = 100, hr.consistency = 1.0)
#> [1] 1
# }
```
