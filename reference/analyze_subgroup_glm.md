# Analyze subgroup for GLM outcomes

GLM counterpart to
[`analyze_subgroup`](https://larry-leon.github.io/forestsearch/reference/analyze_subgroup.md).
Computes per-arm outcome rates, sample sizes, and the treatment effect
estimate (RD, OR, etc.) using the estimator closure.

## Usage

``` r
analyze_subgroup_glm(
  df_sub,
  outcome.name,
  treat.name,
  subgroup_name,
  effect_a = NA,
  estimator_fn = NULL,
  effect_measure = "RD",
  outcome_type = "binary",
  N,
  digits = 4
)
```

## Arguments

- df_sub:

  Data frame for the subgroup.

- outcome.name:

  Character. Name of outcome variable.

- treat.name:

  Character. Name of treatment variable.

- subgroup_name:

  Character. Label for this subgroup.

- effect_a:

  Character. Adjusted effect estimate string (optional).

- estimator_fn:

  Closure from
  [`make_effect_estimator()`](https://larry-leon.github.io/forestsearch/reference/make_effect_estimator.md),
  or NULL.

- effect_measure:

  Character. Effect measure label (e.g., "RD", "OR").

- outcome_type:

  Character. One of `"binary"`, `"continuous"`, or `"count"`. Controls
  formatting: binary displays per-arm means as percentages (`"Rate"`);
  continuous and count display raw means (`"Mean"`). Default `"binary"`.

- N:

  Integer. Total sample size (for percentage calculation).

- digits:

  Integer. Number of decimal places for arm-mean / Diff /
  effect-estimate / CI numeric formatting. Default 4 to provide headroom
  for downstream display reformatting (e.g., via
  `summarize_bootstrap_results(digits = ...)`, which can round down but
  cannot recover precision lost at construction). Note this does not
  affect the count-with-percent-of-total formatting in the `n` and `n1`
  columns, which remain at 1 decimal.

## Value

One-row data frame of formatted subgroup results. Column names are
placeholders (`Subgroup`, `n`, `n1`, `rate0`, `rate1`, `Diff`, `effect`,
optionally `effect_a`) and are overwritten by the calling function
([`SG_tab_estimates_glm`](https://larry-leon.github.io/forestsearch/reference/SG_tab_estimates_glm.md))
with the user-facing labels (e.g., `"Rate(C)"` vs `"Mean(C)"`, the
effect-measure header). Returning a data frame – rather than a named
character vector via `c(...)` – ensures type stability across the ITT
and subgroup call paths and prevents downstream rendering failures in
consumers that go through
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) (e.g.,
gt).
