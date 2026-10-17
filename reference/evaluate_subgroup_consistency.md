# Evaluate Subgroup Consistency (Single-Stage)

Evaluates a single candidate subgroup for treatment-effect consistency
using repeated random data splits (the standard, non-two-stage
algorithm). For each split the subgroup effect is scored and compared
against the harm threshold; the consistency rate is the proportion of
splits meeting it. Called internally by
[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
(via
[`subgroup.consistency`](https://larry-leon.github.io/forestsearch/reference/subgroup.consistency.md));
[`evaluate_consistency_twostage`](https://larry-leon.github.io/forestsearch/reference/evaluate_consistency_twostage.md)
is the screened, early-stopping counterpart used when
`use_twostage = TRUE`.

## Usage

``` r
evaluate_subgroup_consistency(
  m,
  index.Z,
  names.Z,
  df,
  found.hrs,
  n.splits,
  hr.consistency,
  pconsistency.threshold,
  pconsistency.digits = 2,
  maxk,
  confs_labels,
  details = FALSE,
  estimator_fn = NULL,
  consistency_threshold = NULL,
  adjust_covariates = NULL,
  consistency_method = "split",
  glm_resample_spec = NULL,
  skip_consistency = FALSE
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

- n.splits:

  Integer. Number of random splits used to evaluate consistency.

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

- estimator_fn:

  Closure or `NULL`. Effect-estimator closure from
  [`make_effect_estimator`](https://larry-leon.github.io/forestsearch/reference/make_effect_estimator.md)
  for GLM outcomes. Default `NULL` (uses Cox model).

- consistency_threshold:

  Numeric or `NULL`. Threshold for GLM consistency evaluation. Default
  `NULL` (uses `hr.consistency`).

- adjust_covariates:

  Character vector or `NULL`. Additional Cox model terms (e.g.
  `"strata(x1)"`) used when scoring survival subgroups. Referenced
  columns must be present in `df`. Ignored on the GLM path. Default
  `NULL` (unadjusted).

- consistency_method:

  Character. `"split"` (default) or `"resample"`; see
  [`subgroup.consistency`](https://larry-leon.github.io/forestsearch/reference/subgroup.consistency.md).
  When `"resample"` engages it bypasses the split evaluation (the rate
  is computed from a single fit); an unsupported/non-convergent GLM rate
  falls through to the split evaluation.

- glm_resample_spec:

  List or `NULL`. GLM resampling specification threaded from
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md);
  see
  [`subgroup.consistency`](https://larry-leon.github.io/forestsearch/reference/subgroup.consistency.md).
  `NULL` on the survival path.

- skip_consistency:

  Logical. When `TRUE`, skip the consistency estimation entirely and
  return the effect-carry result row (effect from `found.hrs`,
  `Pcons = NA`) without forming the subgroup subset or refitting. Used
  by the `sg_focus = "maxeff"` fast path, which needs effects for the
  argmax but consistency for the selected winner only. Default `FALSE`.

## Value

Named numeric vector with consistency results, or `NULL` if the subgroup
does not meet the criteria.

## Examples

``` r
if (FALSE) { # \dontrun{
# evaluate_subgroup_consistency() is called internally by forestsearch().
# Use forestsearch() as the entry point:
fs <- forestsearch(gbsg,
  confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
  outcome.name = "rfstime", treat.name = "hormon", event.name = "status")
} # }
```
