# Compute GLM Effect Estimate for a Single Subgroup

Fits a GLM within the specified subgroup and returns a one-row data
frame containing the effect estimate, confidence interval, sample size,
and event count. Supports binary, continuous, and count (with offset)
outcomes.

## Usage

``` r
create_glm_row(
  df,
  outcome.name,
  treat.name = "treat",
  effect_measure = "log_OR",
  offset.name = NULL,
  overdispersion = "none",
  conf.level = 0.95,
  min_arm_n = 5L,
  min_arm_events = 3L,
  verbose = FALSE
)
```

## Arguments

- df:

  Data frame for the subgroup (already subset).

- outcome.name:

  Character. Name of the outcome column.

- treat.name:

  Character. Name of the treatment indicator.

- effect_measure:

  Character. Effect measure; see
  [`grf.subg.harm.glm`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.glm.md)
  for options.

- offset.name:

  Character or `NULL`. Name of the person-time column. When supplied,
  `log(offset)` is added as an offset term. Default: `NULL`.

- overdispersion:

  Character. One of `"none"`, `"quasi"`, `"negbin"`. Default: `"none"`.

- conf.level:

  Numeric. Confidence level. Default: `0.95`.

- min_arm_n:

  Integer. Minimum observations per arm; returns `NA` row if either arm
  is below this threshold. Default: `5L`.

- min_arm_events:

  Integer. For binary/count outcomes, minimum events per arm; returns
  `NA` row if check fails. Default: `3L`.

- verbose:

  Logical. Default: `FALSE`.

## Value

A one-row data frame with columns: `N`, `n_treat`, `n_control`,
`events_treat`, `events_control` (for binary/count), `est`, `lower`,
`upper`, `se`, `effect_measure`.

## See also

[`grf.subg.harm.glm`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.glm.md)
