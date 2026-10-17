# Run Guo and He Algorithm 3 on a forestsearch maxeff fit

Bridge from a `forestsearch(sg_focus = "maxeff")` fit to the Guo and He
(2021) Algorithm 3 de-biasing estimator
[`guohe_algorithm3()`](https://larry-leon.github.io/forestsearch/reference/guohe_algorithm3.md).
This is the "post-hoc identified subgroups" application of Guo and He:
forestsearch enumerates a large candidate family and selects the argmax;
Guo and He supplies a de-biased effect and one-sided bound for the
selected subgroup, correcting the winner's-curse optimism of that
selection. The candidate family is materialized from the fit with
[`fs_export_maxeff_family()`](https://larry-leon.github.io/forestsearch/reference/fs_export_maxeff_family.md)
and, by default, verified against forestsearch's own selected-subgroup
membership with
[`fs_assert_membership()`](https://larry-leon.github.io/forestsearch/reference/fs_assert_membership.md)
before any de-biasing runs.

## Usage

``` r
fs_to_guohe(
  fit,
  time = "time_months",
  event = "status",
  treatment = "hormon",
  B = 2000L,
  r = 0.03,
  level = 0.05,
  seed = NULL,
  min_events = 5L,
  parallel = FALSE,
  verify = TRUE
)
```

## Arguments

- fit:

  A `forestsearch` object fitted with `sg_focus = "maxeff"`.

- time, event, treatment:

  Names of the time, event, and 0/1 treatment columns on `fit$df.est`.
  The defaults match the GBSG analyses (`"time_months"`, `"status"`,
  `"hormon"`); set them to the names used in the forestsearch call.

- B:

  Number of bootstrap resamples passed to
  [`guohe_algorithm3()`](https://larry-leon.github.io/forestsearch/reference/guohe_algorithm3.md).

- r:

  Shrinkage tuning parameter, strictly between 0 and 0.5; see
  [`guohe_algorithm3()`](https://larry-leon.github.io/forestsearch/reference/guohe_algorithm3.md)
  and
  [`guohe_adaptive_r()`](https://larry-leon.github.io/forestsearch/reference/guohe_adaptive_r.md).

- level:

  One-sided level; the two-sided interval has coverage `1 - level`.

- seed:

  Optional integer for reproducible resampling.

- min_events:

  Minimum events for a candidate to be estimable.

- parallel:

  Logical; passed to
  [`guohe_algorithm3()`](https://larry-leon.github.io/forestsearch/reference/guohe_algorithm3.md).
  Requires a
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  to have been set by the caller.

- verify:

  Logical; if `TRUE` (default), run
  [`fs_assert_membership()`](https://larry-leon.github.io/forestsearch/reference/fs_assert_membership.md)
  and stop on any mismatch between the reconstructed and the original
  selected-subgroup membership.

## Value

A list with elements `gh` (the full `"guohe_a3"` object), `export` (the
[`fs_export_maxeff_family()`](https://larry-leon.github.io/forestsearch/reference/fs_export_maxeff_family.md)
result), `selected` (the Guo and He argmax candidate name), and the
convenience scalars `naive_hr`, `debiased_hr`, and `bound_hr`
(hazard-ratio scale).

## Details

The de-biasing engine is an independent implementation of the published
method, validated by reproducing the simulation study of Section 5 of
Guo and He (2021).

The Guo and He selection is `orient = +1` (the maximum-hazard-ratio,
i.e. harm, subgroup), so its argmax coincides with the forestsearch
maxeff selection. If the two nevertheless disagree (possible when
candidates are dropped as non-estimable), a warning is issued and the
returned Guo and He objects refer to the Guo and He argmax.

## References

Guo, X. and He, X. (2021). Inference on Selected Subgroups in Clinical
Trials. *Journal of the American Statistical Association*, **116**(535),
1498–1506.
[doi:10.1080/01621459.2020.1740096](https://doi.org/10.1080/01621459.2020.1740096)

## See also

[`fs_export_maxeff_family()`](https://larry-leon.github.io/forestsearch/reference/fs_export_maxeff_family.md),
[`fs_assert_membership()`](https://larry-leon.github.io/forestsearch/reference/fs_assert_membership.md),
[`guohe_algorithm3()`](https://larry-leon.github.io/forestsearch/reference/guohe_algorithm3.md),
[`guohe_adaptive_r()`](https://larry-leon.github.io/forestsearch/reference/guohe_adaptive_r.md)

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
gh <- fs_to_guohe(fit, time = "time", event = "event",
                  treatment = "treat", B = 50, seed = 3)
gh$naive_hr
#> [1] 1.285522
gh$debiased_hr
#> [1] 0.8041826
gh$bound_hr
#> [1] 0.5205394
```
