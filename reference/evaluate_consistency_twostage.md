# Evaluate Consistency (Two-Stage Algorithm)

Evaluates a single subgroup for consistency using a two-stage approach:
Stage 1 screens with fewer splits, Stage 2 uses sequential batched
evaluation with early stopping for efficient evaluation.

## Usage

``` r
evaluate_consistency_twostage(
  m,
  index.Z,
  names.Z,
  df,
  found.hrs,
  hr.consistency,
  pconsistency.threshold,
  pconsistency.digits = 2,
  maxk,
  confs_labels,
  details = FALSE,
  n.splits.screen = 30,
  screen.threshold = NULL,
  n.splits.max = 400,
  batch.size = 20,
  conf.level = 0.95,
  min.valid.screen = 10,
  estimator_fn = NULL,
  consistency_threshold = NULL,
  adjust_covariates = NULL,
  consistency_method = "split",
  glm_resample_spec = NULL
)
```

## Arguments

- m:

  Integer. Index of subgroup to evaluate.

- index.Z:

  data.table or matrix. Factor indicators for all subgroups.

- names.Z:

  Character vector. Names of factor columns.

- df:

  data.frame. Original data with Y, Event, Treat, id columns.

- found.hrs:

  data.table. Subgroup hazard ratio results.

- hr.consistency:

  Numeric. `c2`, **an effect threshold** – the minimum effect required
  in *each* half of a split for that split to count as consistent.
  **Not** a proportion and **not** the consistency rate.
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  calls this `consistency.threshold` / `hr.consistency`.

- pconsistency.threshold:

  Numeric. `p*`, **a proportion** in `[0, 1]` – the fraction of splits
  that must be consistent, where one split's consistency is decided by
  `hr.consistency` (`c2`). `c2` sets the bar; `p*` counts how often it
  is met.

- pconsistency.digits:

  Integer. Number of decimal places to which the consistency proportion
  is rounded **before** it is compared with `pconsistency.threshold`.
  The rounded value decides admission and is the value the selection
  sort orders on, so this setting affects which subgroups are admitted,
  not only what is displayed. At the default of 2, a threshold of `p*`
  admits any candidate whose unrounded proportion is at or above the
  largest 0.01 grid point strictly below `p*` plus half a step (0.005) –
  so `p* = 0.99` admits from about 0.985. Raising the setting makes
  admission finer-grained and produces fewer exact ties for the
  selection rule to break. Under `consistency_method = "resample"` the
  proportion is a continuous closed-form quantity, so the rounding is a
  deliberate coarsening rather than a formatting step. Default: 2

- maxk:

  Integer. Maximum number of factors in a subgroup.

- confs_labels:

  Character vector. Labels for confounders.

- details:

  Logical. Print progress details.

- n.splits.screen:

  Integer. Number of splits for Stage 1 (default 30).

- screen.threshold:

  Numeric. Screening threshold for Stage 1 (default auto-calculated).

- n.splits.max:

  Integer. Maximum total splits (default 400).

- batch.size:

  Integer. Splits per batch in Stage 2 (default 20).

- conf.level:

  Numeric. Confidence level for early stopping (default 0.95).

- min.valid.screen:

  Integer. Minimum valid splits in Stage 1 (default 10).

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
  `"strata(x1)"`) used when scoring survival subgroups. Referenced
  columns must be present in `df`. Ignored on the GLM path. Default
  `NULL` (unadjusted).

- consistency_method:

  Character. `"split"` (default) or `"resample"`; see
  [`subgroup.consistency`](https://larry-leon.github.io/forestsearch/reference/subgroup.consistency.md).
  When `"resample"` engages it bypasses the two-stage split screening
  (the rate is computed from a single fit); an
  unsupported/non-convergent GLM rate falls through to the Stage 1/2
  splitting.

- glm_resample_spec:

  List or `NULL`. GLM resampling specification threaded from
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md);
  see
  [`subgroup.consistency`](https://larry-leon.github.io/forestsearch/reference/subgroup.consistency.md).
  `NULL` on the survival path.

## Value

Named numeric vector with consistency results, or NULL if not met.

## Examples

``` r
if (FALSE) { # \dontrun{
# evaluate_consistency_twostage() is called internally by forestsearch()
# when use_twostage = TRUE. Use forestsearch() as the entry point:
fs <- forestsearch(gbsg,
  confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
  outcome.name = "rfstime", treat.name = "hormon", event.name = "status",
  use_twostage = TRUE)
} # }
```
