# Subgroup Search for Treatment Effect Heterogeneity (Improved, Parallelized)

Searches for subgroups with treatment effect heterogeneity using
combinations of candidate factors. Evaluates subgroups for minimum
prevalence, event counts, and hazard ratio threshold. Parallelizes the
main search loop.

## Usage

``` r
subgroup.search(
  Y,
  Event,
  Treat,
  ID = NULL,
  Z,
  n.min = 30,
  d0.min = 15,
  d1.min = 15,
  hr.threshold = 1,
  max.minutes = 30,
  minp = 0.05,
  rmin = 5,
  details = FALSE,
  maxk = 2,
  parallel_workers = parallel::detectCores(),
  estimator_fn = NULL,
  df_analysis = NULL,
  effect_threshold = NULL,
  effect_measure = NULL,
  outcome_type = NULL,
  adjust_covariates = NULL,
  disable_effect_floor = FALSE
)
```

## Arguments

- Y:

  Numeric vector of outcome (e.g., time-to-event).

- Event:

  Numeric vector of event indicators (0/1).

- Treat:

  Numeric vector of treatment group indicators (0/1).

- ID:

  Optional vector of subject IDs.

- Z:

  Matrix or data frame of candidate subgroup factors (binary
  indicators).

- n.min:

  Integer. Minimum subgroup size.

- d0.min:

  Integer. Minimum per-arm filter. For survival and binary outcomes:
  minimum number of **events** in the control arm within the candidate
  subgroup. Ignored for continuous outcomes (only `n.min` applies).

- d1.min:

  Integer. Same as `d0.min` for the treatment arm.

- hr.threshold:

  Numeric. `c1`, the **screening threshold on a candidate subgroup's own
  effect**: a candidate is retained only if its fitted effect exceeds
  this value. The scale is the one the caller compares on, and it
  differs by path. On the survival path `effect_threshold` is `NULL`, so
  this value reaches the comparison unchanged and is compared against
  the fitted hazard ratio on the **natural** scale
  ([`fit_cox_for_subgroup()`](https://larry-leon.github.io/forestsearch/reference/fit_cox_for_subgroup.md)
  returns `exp(coef)`): 1.25 means HR \>= 1.25, not log(1.25). On the
  GLM paths the value arrives already on the comparison scale that
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  resolved from the effect measure – link scale for the ratio measures
  (OR, RR, IRR), identity for RD, IRD and MD – and is compared against
  the estimator's estimate on that same scale. This is
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)'s
  `effect.threshold` / `hr.threshold`. It is not the per-split threshold
  (`c2`, `hr.consistency`) and not the consistency rate (`p*`,
  `pconsistency.threshold`); neither is applied here – this function
  performs the screening stage only.

- max.minutes:

  Numeric. **Currently inert; scheduled for deprecation in v0.3.0.**
  Previously intended as a wall-clock budget for the combination search,
  it is threaded to `search_combinations_parallel()` and no path
  consults it, so it has no effect on behaviour. Search scope is
  governed by `maxk` and the candidate-factor controls. Retained because
  it is part of a public signature; this description now matches
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)'s,
  which has documented it as inert for some time.

- minp:

  Numeric. Minimum prevalence rate for each factor.

- rmin:

  Integer. Minimum required reduction in sample size when adding a
  factor.

- details:

  Logical. Print details during execution.

- maxk:

  Integer. Maximum number of factors in a subgroup.

- parallel_workers:

  Integer. Number of parallel workers (default: all available cores).

- estimator_fn:

  Closure or `NULL`. Effect-estimator closure from
  [`make_effect_estimator`](https://larry-leon.github.io/forestsearch/reference/make_effect_estimator.md)
  for GLM outcomes. Default `NULL` (uses Cox model).

- df_analysis:

  Data frame or `NULL`. Full analysis dataset, passed to the GLM
  estimator for within-subgroup evaluation.

- effect_threshold:

  Numeric or `NULL`. GLM effect threshold for subgroup selection.
  Default `NULL` (uses `hr.threshold`).

- effect_measure:

  Character or `NULL`. Effect measure name (e.g., `"OR"`, `"HR"`,
  `"RD"`). Used for display labels in the filtering summary. Default
  `NULL`.

- outcome_type:

  Character or `NULL`. One of `"survival"`, `"binary"`, `"continuous"`,
  `"count"`. Controls d0.min/d1.min interpretation: events per arm for
  binary/survival, ignored for continuous. Default `NULL` (inferred from
  `estimator_fn` presence).

- adjust_covariates:

  Character vector or `NULL`. Additional Cox model terms (e.g.
  `"strata(x1)"` or `"x2"`) used to adjust the candidate-search Cox
  scorer on the survival path. Referenced columns are sourced from
  `df_analysis`. Ignored on the GLM path. Default `NULL`
  (treatment-only).

- disable_effect_floor:

  Logical. When `TRUE`, the Status-6 effect threshold (`hr.threshold`)
  is not applied, so the full estimable candidate family is retained
  regardless of effect size. Set by
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  when `sg_focus = "maxeff"` so the effect maximiser ranges over an
  unconditional family. Default `FALSE`.

## Value

List with found subgroups, maximum HR, search time, configuration info,
and filtering statistics.

## Examples

``` r
if (FALSE) { # \dontrun{
library(survival)
df <- survival::gbsg
df$grade3 <- as.integer(df$grade == "3")
Z <- df[, c("age", "meno", "size", "grade3", "nodes", "pgr", "er")]
res <- subgroup.search(Y = df$rfstime, Event = df$status,
                       Treat = df$hormon, Z = Z,
                       hr.threshold = 1.25)
} # }
```
