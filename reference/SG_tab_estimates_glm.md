# Subgroup Summary Table for GLM Outcomes

GLM counterpart to
[`SG_tab_estimates`](https://larry-leon.github.io/forestsearch/reference/SG_tab_estimates.md).
Produces a summary data frame with per-arm outcome rates and the GLM
effect estimate for each subgroup.

## Usage

``` r
SG_tab_estimates_glm(
  df,
  SG_flag,
  outcome.name,
  treat.name,
  estimator_fn = NULL,
  effect_measure = "RD",
  outcome_type = "binary",
  effect_a_1 = NA,
  effect_a_0 = NA,
  sg1_name = "Recommend",
  sg0_name = "Questionable",
  est.scale = "hr",
  digits = 4
)
```

## Arguments

- df:

  Data frame with a `treat.recommend` column (or ITT).

- SG_flag:

  Character. "ITT" for the full population, or the name of the subgroup
  flag column (e.g., "treat.recommend").

- outcome.name:

  Character. Name of outcome variable.

- treat.name:

  Character. Name of treatment variable.

- estimator_fn:

  Closure from
  [`make_effect_estimator()`](https://larry-leon.github.io/forestsearch/reference/make_effect_estimator.md),
  or NULL.

- effect_measure:

  Character. Effect measure label (e.g., "RD", "OR").

- outcome_type:

  Character. One of `"binary"`, `"continuous"`, or `"count"`. Controls
  column labels and formatting: binary uses `"Rate(C)"`/`"Rate(T)"` with
  percentages; continuous and count use `"Mean(C)"`/`"Mean(T)"` with raw
  values. Default `"binary"`.

- effect_a_1:

  Character. Adjusted effect for subgroup 1 (optional).

- effect_a_0:

  Character. Adjusted effect for subgroup 0 (optional).

- sg1_name:

  Character. Label for subgroup 1 (treat.recommend == 1).

- sg0_name:

  Character. Label for subgroup 0 (treat.recommend == 0).

- est.scale:

  Character. Effect scale ("hr" or "1/hr").

- digits:

  Integer. Decimal places for numeric formatting; passed through to
  [`analyze_subgroup_glm`](https://larry-leon.github.io/forestsearch/reference/analyze_subgroup_glm.md).
  Default 4 to provide round-down headroom for downstream display
  reformatting (e.g., `summarize_bootstrap_results(digits = ...)`).

## Value

Data frame. One row when `SG_flag = "ITT"`; two rows (sg0 then sg1) when
`SG_flag = "treat.recommend"`. Columns: `Subgroup`, `n`, `n1`, `Rate(C)`
or `Mean(C)` (binary vs continuous/count), `Rate(T)` or `Mean(T)`,
`Diff`, and the effect-measure header `"\{effect_measure\} (95% CI)"`.
When `effect_a_*` is supplied, an additional `"\{effect_measure\}*"`
column carries the bias-corrected estimate.
