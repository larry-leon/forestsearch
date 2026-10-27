# Cross-validated choice of the shrinkage tuning parameter r

Implements Algorithm 2 of Guo and He (2021), Section 2.5: a v-fold
cross-validated choice of the shrinkage tuning parameter `r` used by
[`guohe_algorithm3()`](https://larry-leon.github.io/forestsearch/reference/guohe_algorithm3.md).
For each candidate `r` and each fold, the bias-reduced estimator is
computed on the training folds and compared, on the held-out fold,
against every candidate subgroup's effect estimate with its sampling
variance subtracted; the `r` minimising the resulting CV objective is
selected.

## Usage

``` r
guohe_adaptive_r(
  data,
  outcome = c("survival", "binary", "continuous"),
  treatment,
  candidates,
  time = NULL,
  event = NULL,
  y = NULL,
  orient = -1,
  r_grid = c(0.03, 0.1, 0.2, 0.3, 0.4, 0.45),
  v = 5L,
  B = 200L,
  level = 0.05,
  seed = NULL,
  min_events = 5L,
  refit = TRUE,
  adjust_covariates = NULL,
  fast = NULL,
  parallel = FALSE
)
```

## Arguments

- data:

  A data frame.

- outcome:

  One of `"survival"`, `"binary"`, `"continuous"`.

- treatment:

  Name of the 0/1 treatment column.

- candidates:

  Character vector of names of 0/1 membership columns forming the
  enumerated family.

- time, event:

  Survival only: names of the time and event columns.

- y:

  Binary or continuous only: name of the outcome column.

- orient:

  Either `-1` or `+1`; sign orienting the effect of interest.

- r_grid:

  Numeric vector of candidate values in (0, 0.5).

- v:

  Number of cross-validation folds.

- B:

  Number of bootstrap resamples.

- level:

  One-sided level; the two-sided interval has coverage `1 - level`.

- seed:

  Optional integer for reproducible resampling.

- min_events:

  Minimum events (survival) or minimum per-arm count for a candidate to
  be estimable.

- refit:

  Logical; if `TRUE`, also return a full-data fit at the selected `r`.

- adjust_covariates:

  This option is implemented but has not been tested against the
  manuscript analyses or the simulation studies of León and Anderson
  (2026a, 2026b), which use the unadjusted treatment-only fits. Optional
  character vector of covariate names for the within-subgroup model,
  matching the argument of the same name in
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)'s
  `.estimate_or()` / `.make_lm_estimator()`. Default `NULL` is the
  unadjusted model. Supplying it DISABLES the closed-form fast path,
  because an adjusted coefficient is a partial regression coefficient,
  not a 2x2 log-odds-ratio or a difference of arm means.

  NOTE:
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  has a SECOND adjustment axis, `ps_adjust_method` (`"iptw"` weighted
  fitting, `"dr_gcomp"` Bang-Robins G-computation). Those modes are not
  implemented in this comparator, and the closed forms would be invalid
  under them too – IPTW fits a weighted logistic, and DR reports a
  marginal log-OR from G-computed probabilities; neither equals the 2x2
  statistic. If the comparator is extended to those modes, the `fast_ok`
  gate must be extended with them.

- fast:

  Logical or `NULL`. Closed-form evaluation of the within-subgroup
  effect, available only for UNADJUSTED binary and continuous outcomes.
  These are EXACT, not approximations: with a single binary covariate
  the logistic model is saturated, so its MLE is the 2x2 log-odds-ratio,
  and the Gaussian coefficient is the difference of arm means. `NULL`
  (default) auto-enables it whenever it is valid. `TRUE` errors if the
  configuration makes it invalid (survival, or any `adjust_covariates`),
  rather than silently falling back.

- parallel:

  Logical; if `TRUE`, distribute the bootstrap model fitting with
  `future.apply`. Requires a
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  to have been set by the caller. Results are identical to
  `parallel = FALSE` for the same `seed`, because resample indices are
  drawn serially before any distribution. Do NOT enable this inside an
  already-parallel outer loop.

## Value

An object of class `"guohe_ar"`: the selected `r` (`r_hat`), the CV
objective over the grid (`objective`), the per-candidate objective at
the selected `r` (`per_candidate`), bookkeeping fields (`r_grid`, `v`,
`B`, `n`, `n_candidates`, `fast`, `outcome`, `orient`, `level`,
`adjust_covariates`), and, when `refit = TRUE`, the refitted full-data
[`guohe_algorithm3()`](https://larry-leon.github.io/forestsearch/reference/guohe_algorithm3.md)
result (`fit`).

## Details

This is an independent implementation of the published method, validated
by reproducing the simulation study of Section 5 of Guo and He (2021).

## References

Guo, X. and He, X. (2021). Inference on Selected Subgroups in Clinical
Trials. *Journal of the American Statistical Association*, **116**(535),
1498–1506.
[doi:10.1080/01621459.2020.1740096](https://doi.org/10.1080/01621459.2020.1740096)

## Examples

``` r
set.seed(1)
n <- 160
df <- data.frame(
  treat = rbinom(n, 1, 0.5),
  y     = rbinom(n, 1, 0.3),
  g1    = rbinom(n, 1, 0.5),
  g2    = rbinom(n, 1, 0.5),
  g3    = rbinom(n, 1, 0.5)
)
cv <- guohe_adaptive_r(
  df, outcome = "binary", treatment = "treat",
  candidates = c("g1", "g2", "g3"), y = "y",
  orient = +1, r_grid = c(0.05, 0.25), v = 2, B = 25, seed = 2
)
cv$r_hat
#> [1] 0.05
```
