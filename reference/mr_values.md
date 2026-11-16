# Multiplier-resampling (MR) inference quantities as a tidy data frame

Extracts the MR quantities of a
[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
fit (its `mr_inference` element) or of a bare
[`fs_mr_inference()`](https://larry-leon.github.io/forestsearch/reference/fs_mr_inference.md)
result into a plain data frame, one row per region: `"H"` (the selected
subgroup) and, when the complement block exists, `"Hc"` (its
complement). This is pure extraction: every value already exists in
`mr_inference` and nothing is recomputed, apart from the optional
[`log()`](https://rdrr.io/r/base/Log.html) of `scale = "log"`.

## Usage

``` r
mr_values(x, scale = c("effect", "log"))
```

## Arguments

- x:

  A
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  fit with a non-`NULL` `mr_inference` element, or the list returned by
  [`fs_mr_inference()`](https://larry-leon.github.io/forestsearch/reference/fs_mr_inference.md)
  (including its short variant, which yields an `"H"` row that is `NA`
  apart from `n_family`).

- scale:

  `"effect"` (the default) returns every column exactly as stored.
  `"log"` applies [`log()`](https://rdrr.io/r/base/Log.html) to the
  columns listed as *effect scale* below, for ratio measures (`HR`,
  `OR`, `RR`, `IRR`); it is an error for identity measures (`MD`, `RD`).
  The other columns are returned as stored under both settings.

## Value

A plain `data.frame` with a fixed set of columns, whatever the fit
computed: a quantity the fit did not compute is `NA`, and no column is
ever dropped. The `"Hc"` row is present whenever `complement` is
non-`NULL`; if the complement could not be fit, that row is all `NA`.
Below, *working scale* means log for ratio measures and the effect
itself for identity measures. Paths are relative to `mr_inference`.

- `row`: `"H"` or `"Hc"`.

- `n`: Patients in the region: `n_selected` on `"H"`, `complement$n` on
  `"Hc"`. Count. `NA` in the short return, or when the complement could
  not be fit.

- `naive_est`, `naive_lower`, `naive_upper`: The unadjusted estimate
  with its two-sided 95% robust interval: `naive$est`, `naive$lower`,
  `naive$upper`, and `complement$naive` on `"Hc"`. Effect scale. `NA` in
  the short return, or when the complement could not be fit.

- `debiased_est`: The two-term de-biased estimate: `debiased$est`, and
  `complement$debiased$est` on `"Hc"`. Effect scale. `NA` as
  `naive_est`.

- `ij_lower`, `ij_upper`: The two-sided 95% two-term IJ interval:
  `debiased$lower`, `debiased$upper`, and `complement$debiased` on
  `"Hc"`. Effect scale. Also `NA` under `ci_method = "wald"` or
  `ij_residual != "two_term"`, when those two stored limits are not this
  interval.

- `ij_se`: The two-term IJ standard error: `debiased$se_ij_two_term`,
  and `complement$debiased$se_ij_two_term` on `"Hc"`. Working scale,
  under both settings of `scale`. It is the subgroup robust SE when the
  IJ variance is unusable (`debiased$ij_source` is `"wald_fallback"`).
  `NA` as `naive_est`.

- `field_est2`: The field second-order point estimate: `field$est2`, and
  `field$complement$est2` on `"Hc"`. Effect scale. `NA` unless
  `ci_method = "field"`, and on `"Hc"` unless the complement field ran
  (`field_complement = TRUE` with `include_complement = TRUE`), or when
  that block returned a note instead of numbers.

- `field_est2_s`: The studentized field-s counterpart on `"Hc"`:
  `field$complement$est2_s`. Effect scale. `NA` on `"H"`, and on `"Hc"`
  as `field_est2` or under `field_scale_complement = "none"`.

- `field_lower_1s`: The field one-sided 95% lower bound on `"H"`:
  `field$lower_1s`. Effect scale. `NA` on `"Hc"`, and unless
  `ci_method = "field"`.

- `field_upper_1s`: The field one-sided 95% upper bound on `"Hc"`:
  `field$complement$upper_1s`. Effect scale. `NA` on `"H"`, and on
  `"Hc"` as `field_est2`.

- `field_upper_1s_s`: The field-s one-sided 95% upper bound on `"Hc"`:
  `field$complement$upper_1s_s`. Effect scale. `NA` as `field_est2_s`.

- `bonf_lower_H`, `bonf_upper_Hc`: The Bonferroni pair (each bound at
  97.5%, jointly 95%), from the field's joint block:
  `field$joint$bonf_lower_H` on `"H"`, `field$joint$bonf_upper_Hc` on
  `"Hc"`. Effect scale. `NA` on the other row, and when `field$joint` is
  absent, which it is unless the complement field ran, or has no bounds.

- `bonf_lower_H_s`, `bonf_upper_Hc_s`: The same pair from the field-s
  joint block, `field$joint_s`. Effect scale. `NA` on the other row, and
  when `field$joint_s` is absent, as for `bonf_lower_H`, or under
  `field_scale_complement = "none"`.

- `selection_bias`: The selection-bias term: `selection_bias`. Working
  scale, under both settings of `scale`. Fit-level: `"H"` row only. `NA`
  in the short return.

- `p_hat`: The selected candidate's re-selection frequency:
  `reselection$p_hat[[selected_index]]`. A probability. `"H"` row only.
  `NA` under `return_reselection = FALSE`, or in the short return.

- `n_family`: The candidate family's size: `n_family`. Count. `"H"` row
  only.

- `n_selected`: Patients in the selected subgroup: `n_selected`. Count.
  `"H"` row only. `NA` in the short return.

- `timing_seconds`: Seconds in the main multiplier stage:
  `timing_seconds`. `"H"` row only. `NA` in the short return.

- `field_timing_seconds`: Seconds in the field block:
  `field$timing_seconds`. `"H"` row only. `NA` unless
  `ci_method = "field"`.

- `field_complement_timing_seconds`: Seconds in the complement field
  block: `field$complement$timing_seconds`. `"H"` row only. `NA` as
  `field_est2` on `"Hc"`.

## Details

Only the evaluated construction set is exposed: the naive estimate, the
two-term infinitesimal-jackknife (IJ) de-biased estimate, the field
one-sided bounds (lower on the selected subgroup, upper on its
complement, plain and studentized) and the Bonferroni pair, with the
fit-level diagnostics.

## See also

[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
for the `mr_inference` switch;
[`fs_mr_inference()`](https://larry-leon.github.io/forestsearch/reference/fs_mr_inference.md)
for the computation;
[`mr_estimates_table()`](https://larry-leon.github.io/forestsearch/reference/mr_estimates_table.md)
to render the estimates as a table.

## Examples

``` r
set.seed(11)
n <- 300
d <- data.frame(time = rexp(n, 0.2), event = rbinom(n, 1, 0.8),
                treat = rbinom(n, 1, 0.5), x1 = rnorm(n), x2 = rnorm(n))
family <- list(g1 = which(d$x1 > 0), g2 = which(d$x2 > 0),
               g3 = which(d$x1 + d$x2 > 0.3))
spec <- list(outcome_type = "survival", effect_measure = "HR",
             treat.name = "treat", outcome.name = "time",
             event.name = "event", offset.name = NULL,
             adjust_covariates = NULL, adverse_outcome = TRUE)
fit <- fs_mr_inference(
  df = d, candidates = family, spec = spec,
  selected_members = family$g1,
  admission = list(effect_floor = NULL, consistency = NULL),
  reselection = "maxeff", draws = 300L, seed = 4242L,
  ci_method = "field", include_complement = TRUE,
  field_R_out = 200L, field_R_in = 50L)
mr_values(fit)
#>   row   n naive_est naive_lower naive_upper debiased_est  ij_lower ij_upper
#> 1   H 151 0.9279521   0.6438566    1.337402    0.8627701 0.4005693 1.858286
#> 2  Hc 149 0.7760111   0.5378543    1.119621    0.8694547 0.4515283 1.674206
#>       ij_se field_est2 field_est2_s field_lower_1s field_upper_1s
#> 1 0.3914672  0.8563254           NA      0.6222182             NA
#> 2 0.3343062  0.9039651    0.9050143             NA       1.230306
#>   field_upper_1s_s bonf_lower_H bonf_upper_Hc bonf_lower_H_s bonf_upper_Hc_s
#> 1               NA    0.5935985            NA      0.5935985              NA
#> 2         1.230768           NA      1.416401             NA        1.408442
#>   selection_bias p_hat n_family n_selected timing_seconds field_timing_seconds
#> 1     0.07958234  0.46        3        151          0.004                0.014
#> 2             NA    NA       NA         NA             NA                   NA
#>   field_complement_timing_seconds
#> 1                           0.013
#> 2                              NA
mr_values(fit, scale = "log")[, c("row", "naive_est", "ij_se", "selection_bias")]
#>   row  naive_est     ij_se selection_bias
#> 1   H -0.0747752 0.3914672     0.07958234
#> 2  Hc -0.2535885 0.3343062             NA
```
