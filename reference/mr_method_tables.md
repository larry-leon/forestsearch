# Naive, bootstrap and multiplier estimates, side by side, and their cost

Plain data frames for the comparison of a
[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
analysis's naive estimates with its full-bootstrap and
multiplier-resampling corrections, method by method; the caller styles
them.

## Usage

``` r
mr_method_comparison(fits, boots, est.scale = "hr", digits = 2)

mr_cost_table(fits, boots, boot_seconds = NULL, project_B = 1000L)
```

## Arguments

- fits:

  A named list of
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  fits with multiplier resampling and the complement de-biased. The
  names label the methods.

- boots:

  A list of
  [`forestsearch_bootstrap_dofuture()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md)
  results with the same names.

- est.scale:

  `"hr"` (the default), or `"1/hr"`, which inverts every ratio (limits
  swap) and is accepted for ratio measures only.

- digits:

  Decimals of the estimates. Default 2.

- boot_seconds:

  The wall-clock seconds of each bootstrap, a vector named by method.
  `NULL` (the default) reads what each result stores, the time of its
  replicate loop (`attr(boot$results, "timing")$total_minutes`), which
  leaves out the call's setup and summaries.

- project_B:

  The number of replicates to project the cost to. Default 1000.

## Value

`mr_method_comparison()`: a plain data frame with two rows per method,
`H` and `Hc`: `method`, `region`, `n` (patients, with their percent),
`naive` (estimate with its robust 95% interval), `bootstrap` (the point
estimate), and `multiplier` (the two-term de-biased estimate with its IJ
95% interval). The interval reads `NA` on a fit whose stored limits are
not the two-term IJ interval; see
[`mr_values()`](https://larry-leon.github.io/forestsearch/reference/mr_values.md).

`mr_cost_table()`: a plain data frame with one row per method: `method`,
`mr_seconds`, `bootstrap_B`, `bootstrap_seconds`,
`seconds_per_replicate`, and `projected_B_<project_B>`.

## Details

**What the comparison shows.** The table compares displays; it does not
adjudicate one quantity. The bootstrap column carries the point estimate
only: at a small `B` that point is an illustration, and its interval is
not a result. For families generated from a model (the DINA and GRF
paths), the full bootstrap, which re-runs the search in every replicate,
and multiplier resampling, which perturbs one fixed family, target
different estimands.

`mr_cost_table()`'s `mr_seconds` is the multiplier stage as the fit
records it (`mr_inference$timing_seconds`), and the projection is linear
in `B` at the worker count the bootstrap ran on.

## See also

[`mr_values()`](https://larry-leon.github.io/forestsearch/reference/mr_values.md),
[`mr_estimates_table()`](https://larry-leon.github.io/forestsearch/reference/mr_estimates_table.md),
[`forestsearch_bootstrap_dofuture()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md).

## Examples

``` r
# \donttest{
set.seed(1)
n <- 500
d <- data.frame(id = seq_len(n), treat = rbinom(n, 1, 0.5), x1 = rnorm(n), x2 = rnorm(n))
d$time <- rexp(n, 0.1 * ifelse(d$treat == 1, ifelse(d$x1 < -0.6, 3, 0.45), 1))
d$event <- as.integer(d$time < 30)
d$time <- pmin(d$time, 30)
fit <- forestsearch(
  df.analysis = d, outcome.name = "time", event.name = "event", treat.name = "treat",
  id.name = "id", confounders.name = c("x1", "x2"), use_grf = FALSE, use_lasso = FALSE,
  maxk = 1, fs.splits = 200, sg_focus = "hrMaxSG", mr_inference = TRUE,
  mr_inference_args = list(draws = 200, field_R_out = 50, field_R_in = 25),
  parallel_args = list(plan = "sequential"), quiet = TRUE, details = FALSE,
  plot.sg = FALSE, show_candidate_summary = FALSE)
boot <- forestsearch_bootstrap_dofuture(fit, nb_boots = 2,
                                        parallel_args = list(plan = "sequential"))
#> Warning: nb_boots < 100 may produce unreliable confidence intervals. Recommend at least 500 for production use.
#> Bootstrap parallel config: Using user-provided settings
#> System max cores available: 4
#> Bootstrap will use:  workers with 'sequential' plan
#> Bootstrap iteration 1/2
#> Bootstrap iteration 2/2
mr_method_comparison(list(FS = fit), list(FS = boot))
#>   method region         n             naive bootstrap        multiplier
#> 1     FS      H 134 (27%) 3.56 (2.45, 5.16)      5.07 3.62 (1.50, 8.74)
#> 2     FS     Hc 366 (73%) 0.54 (0.43, 0.67)      0.58 0.52 (0.29, 0.93)
mr_cost_table(list(FS = fit), list(FS = boot))
#>   method mr_seconds bootstrap_B bootstrap_seconds seconds_per_replicate
#> 1     FS      0.007           2             3.594                 1.797
#>   projected_B_1000
#> 1             1797
# }
```
