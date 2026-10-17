# Assert the exported selected-subgroup membership matches forestsearch's own

Integrity check for
[`fs_export_maxeff_family()`](https://larry-leon.github.io/forestsearch/reference/fs_export_maxeff_family.md):
the materialized membership column of the selected candidate must equal,
row for row, the per-subject membership forestsearch itself computed
(`fit$grp.consistency$sg.harm.id`). A mismatch means the label-to-dummy
reconstruction is wrong, so this stops with an error rather than letting
a silently different subgroup be de-biased.
[`fs_to_guohe()`](https://larry-leon.github.io/forestsearch/reference/fs_to_guohe.md)
runs this check by default (`verify = TRUE`).

## Usage

``` r
fs_assert_membership(exp, fit)
```

## Arguments

- exp:

  The list returned by
  [`fs_export_maxeff_family()`](https://larry-leon.github.io/forestsearch/reference/fs_export_maxeff_family.md).

- fit:

  The `forestsearch` object `exp` was exported from.

## Value

`TRUE`, invisibly, if the memberships agree; otherwise an error is
thrown.

## See also

[`fs_export_maxeff_family()`](https://larry-leon.github.io/forestsearch/reference/fs_export_maxeff_family.md),
[`fs_to_guohe()`](https://larry-leon.github.io/forestsearch/reference/fs_to_guohe.md)

## Examples

``` r
set.seed(11)
n <- 200
df <- data.frame(id = seq_len(n), treat = rbinom(n, 1, 0.5),
                 age = round(rnorm(n, 55, 12)),
                 bm = factor(rbinom(n, 1, 0.5)))
harm <- df$age <= 50 & df$bm == "1"
tt <- rexp(n, 0.05 * exp(log(2.5) * df$treat * harm))
df$time <- pmin(tt, 60)
df$event <- as.integer(tt <= 60)
fit <- forestsearch(
  df.analysis = df,
  outcome.name = "time", event.name = "event",
  treat.name = "treat", id.name = "id",
  confounders.name = c("age", "bm"),
  sg_focus = "maxeff", use_grf = FALSE, use_lasso = FALSE,
  fs.splits = 50, n.min = 30, d0.min = 5, d1.min = 5,
  hr.threshold = 1.1, hr.consistency = 1.0,
  pconsistency.threshold = 0.5, maxk = 2,
  details = FALSE, plot.sg = FALSE,
  parallel_args = list(plan = "sequential", workers = 1L,
                       show_message = FALSE)
)
#> Warning: sg_focus = 'maxeff' selects the effect maximiser with no auxiliary conditions; overriding: pconsistency.threshold -> 0, stop_threshold -> NULL, use_twostage -> FALSE, minp -> 0. The subgroup-search effect floor (hr.threshold) is disabled, redundancy pruning is relaxed to exact-duplicate removal (rmin -> 0), and consistency-stage duplicate removal collapses identical-membership candidates to their fewest-cut representative -- so the argmax ranges over the full estimable family with no arbitrary tie-break.
#> 
#> [forestsearch] Subgroup Identification Configuration
#>   Outcome type:   survival (Cox PH)
#>   Effect measure: HR (hazard ratio)
#>   Scale:          log (ratio)
#> 
#>   Screening:      HR >= 1.10  (log: 0.0953)
#>   Consistency:    HR >= 1.00  (log: 0.0000)
#>   Consistency rate threshold: 0%
#> 
#>   Interpretation:
#>     Candidates: subgroup hazard ratio >= 1.10
#>     Splits:     each split HR >= 1.00
exp_fam <- fs_export_maxeff_family(fit)
fs_assert_membership(exp_fam, fit)
```
