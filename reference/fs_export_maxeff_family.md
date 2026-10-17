# Materialize the maxeff candidate family as 0/1 membership columns

Reconstructs, for every candidate in the deduplicated maxeff family kept
on a `forestsearch(sg_focus = "maxeff")` fit, its per-subject 0/1
membership on the estimation frame `fit$df.est`. Each candidate's cut
labels (e.g. `"{er <= 0}"`, or a negated cut `"!{size <= 20}"`) are
mapped back to the aligned dummy columns of `fit$df.est` via
`fit$confounders.evaluated` / `fit$confounders.candidate`, and a
candidate's membership is the intersection (AND) over its cuts. The
result is the materialized family that
[`guohe_algorithm3()`](https://larry-leon.github.io/forestsearch/reference/guohe_algorithm3.md)
consumes;
[`fs_to_guohe()`](https://larry-leon.github.io/forestsearch/reference/fs_to_guohe.md)
calls this internally.

## Usage

``` r
fs_export_maxeff_family(fit, prefix = "sg_")
```

## Arguments

- fit:

  A `forestsearch` object fitted with `sg_focus = "maxeff"`.

- prefix:

  Character prefix for the generated membership column names (`sg_0001`,
  `sg_0002`, ... by default).

## Value

A list with elements `data` (`fit$df.est` with one appended 0/1
membership column per candidate), `candidates` (the appended column
names), `selected` (the column name of the forestsearch-selected
candidate; the family table is sorted selected-first), and `family` (the
candidate family table `fit$grp.consistency$out_sg$result`).

## See also

[`fs_to_guohe()`](https://larry-leon.github.io/forestsearch/reference/fs_to_guohe.md),
[`fs_assert_membership()`](https://larry-leon.github.io/forestsearch/reference/fs_assert_membership.md)

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
head(exp_fam$candidates)
#> [1] "sg_0001" "sg_0002" "sg_0003" "sg_0004" "sg_0005" "sg_0006"
exp_fam$selected
#> [1] "sg_0001"
```
