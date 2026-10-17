# Summarize an extreme-subgroups simulation study

Computes, from the raw matrices of a
[`run_subgroup_sims()`](https://larry-leon.github.io/forestsearch/reference/run_subgroup_sims.md)
result (or a loaded RDS payload coerced with `class<-`), every statistic
used by the vignettes' results table and forest panels: 1st/50th/99th
percentile matrices for the estimate and its upper bound, convergence
counts, mean subgroup N, tail probabilities, unconditional median
display strings, the formatted `results_tbl`, category labels and
colours, validity masks, the single/combination panel index vectors, and
the high-risk panel (threshold filter, ITT anchor, ordering).

## Usage

``` r
# S3 method for class 'subgroup_sims'
summary(
  object,
  hr_true = object$hr_true,
  probs = c(0.01, 0.5, 0.99),
  est_thresholds = NULL,
  ub_thresholds = NULL,
  est_label = NULL,
  ub_label = NULL,
  single_cats = c("ITT", "Clinical", "Continuous"),
  combo_cats = c("Interaction", "3-way", "Random"),
  pr_thresh_highrisk = 0.1,
  itt_name = "All Patients",
  cat_cols = c(ITT = "black", Clinical = "steelblue", Continuous = "darkgreen",
    Interaction = "purple4", `3-way` = "violetred3", Random = "darkorange"),
  ...
)
```

## Arguments

- object:

  A `"subgroup_sims"` object.

- hr_true:

  True effect on the fitter's estimate scale, for reference lines (a
  hazard ratio for Cox fits, a mean difference for
  [`subgroup_glm()`](https://larry-leon.github.io/forestsearch/reference/subgroup_glm.md)
  fits); defaults to the value stored on `object`.

- probs:

  Quantile probabilities for the ECI matrices; the default
  `c(0.01, 0.50, 0.99)` matches the vignettes and downstream code
  indexes column 2 as the median.

- est_thresholds:

  Length-2 numeric `c(low, high)`: the estimate tails are
  `Pr(est < low)` and `Pr(est > high)`. `NULL` resolves via the `effect`
  metadata, then the HR legacy `c(0.5, 1.0)`. `NA` entries disable that
  tail (see above).

- ub_thresholds:

  Length-2 numeric: the upper-bound tails are
  `Pr(UB >= ub_thresholds[1])` and `Pr(UB >= ub_thresholds[2])`; the
  first also drives the high-risk panel. Resolution as for
  `est_thresholds` (HR legacy `c(2, 3)`).

- est_label, ub_label:

  Display labels for the estimate and its upper bound in the table
  headers (legacy `"HR"` / `"UB"`).

- single_cats, combo_cats:

  Category labels routed to the single-variable and combination forest
  panels.

- pr_thresh_highrisk:

  High-risk panel inclusion threshold on the first upper-bound tail
  probability, default `0.10`.

- itt_name:

  Name of the ITT row, always included in and anchored at the top of the
  high-risk panel.

- cat_cols:

  Named colour map for categories.

- ...:

  Unused.

## Value

An object of class `"subgroup_sims_summary"`: a list with fields named
as in the vignettes (`hr_q`, `ub_q`, `n_valid`, `pct_valid`, `mean_n`,
`pr_hr_lt050`, `pr_hr_gt1`, `pr_ub_ge2`, `pr_ub_ge3`, `mhr_uncond`,
`mub_uncond`, `results_tbl`, `sg_names`, `cat`, `cat_cols`, `ok`,
`ok_ub`, `idx_hr_single`, `idx_hr_combo`, `idx_ub_single`,
`idx_ub_combo`, `n_single`, `n_combo`, `highrisk`, `n_sims`, `design`,
`hr_true`). Legacy calls return exactly this set; effect-aware calls
append `effect` (the result's metadata, possibly `NULL` when only
explicit overrides were given), `thresholds` (`list(est =, ub =)` as
resolved), and `labels` (`list(est =, ub =)`).

## Details

Field names are structural and retained across outcome types (`hr_q`,
`pr_hr_lt050`, `pr_ub_ge2`, ...), exactly as `sim_hrs` and
`hr.threshold` serve generic duty elsewhere in the package: on a GLM
result they hold the fitter's estimate-scale statistics. Estimates
follow the estimator layer's harm-positive normalization, so the upper
tail (`Pr(est > null)`, `Pr(UB >= t)`) always points toward harm
whatever the raw outcome's direction.

Threshold and label resolution: explicit arguments win; otherwise the
result's `effect` metadata supplies them
([`subgroup_glm()`](https://larry-leon.github.io/forestsearch/reference/subgroup_glm.md)
stamps scale-aware defaults from Phase 4.4 – ratio measures carry the
legacy tails `c(0.5, 1)` / `c(2, 3)`, identity measures
`c(NA, null_value)` / `c(NA, NA)` – with labels from the effect
measure); otherwise the HR legacy values `c(0.5, 1.0)` / `c(2, 3)` with
the vignettes' literal headers. An `NA` threshold degrades gracefully:
its probability field is `NA` and its table column renders `"-"` under a
placeholder header.

Structurally empty subgroups (all-`NA` columns) yield `NA` medians and
`NaN` means/proportions, exactly as in the vignettes; the formatted
table renders them as `"-"`.

## See also

[`run_subgroup_sims()`](https://larry-leon.github.io/forestsearch/reference/run_subgroup_sims.md),
[`subgroup_glm()`](https://larry-leon.github.io/forestsearch/reference/subgroup_glm.md)
