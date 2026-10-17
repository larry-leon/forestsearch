# Design-time feasibility of a DGM's planted region

Draws replicates through the DGM's own generator and reports, per sample
size, how often the planted harm region \\Q\\ would be **undeclarable**
(too small for the search's own size test), how often it falls under the
per-arm events floor, and how often the analysis estimand does not exist
on it.

## Usage

``` r
fs_dgm_feasibility(
  dgm,
  n = c(500, 750, 1000, 2000),
  n.min = 60,
  d0.min = 10,
  d1.min = 10,
  n_rep = 200L,
  tolerance = 0.05,
  effect_measure = NULL,
  seed = NULL,
  rand_ratio = 1
)

# S3 method for class 'fs_dgm_feasibility'
print(x, ...)
```

## Arguments

- dgm:

  An object of class `"glm_dgm"` – the interface the
  operating-characteristics family accepts
  ([`fs_oc_predict`](https://larry-leon.github.io/forestsearch/reference/fs_oc_predict.md),
  [`fs_oc_grid`](https://larry-leon.github.io/forestsearch/reference/fs_oc_grid.md),
  which validate `dgm$df_super` and its `flag_harm` column). Survival
  DGMs are not accepted; see Note.

- n:

  Numeric vector of sample sizes to evaluate.

- n.min, d0.min, d1.min:

  The search's own floors, passed so the report describes the campaign
  that will actually run. Defaults are the package defaults (60, 10,
  10). Nothing is imposed: these are read, not applied.

- n_rep:

  Integer. Replicates per sample size (default 200).

- tolerance:

  Numeric. The largest undeclarable share that still counts as feasible
  (default 0.05).

- effect_measure:

  Character or `NULL`. The estimand whose existence condition is
  checked. `NULL` (default) reads `dgm$effect_measure`.

- seed:

  Integer or `NULL`. Base seed; replicate `i` of sample size `n` is
  drawn at `seed + i`.

- rand_ratio:

  Numeric. Passed to
  [`simulate_from_glm_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_glm_dgm.md).

- x:

  An `fs_dgm_feasibility` object.

- ...:

  Ignored.

## Value

An object of class `"fs_dgm_feasibility"`: a list with `$table` (one row
per `n`), `$feasible` (a single logical – `TRUE` when every undeclarable
share is at or below `tolerance`), and `$args`.

## Details

The motivating case: at `n = 500` a region averaging 48.5 subjects is at
or below `n.min = 60` in 96\\ structurally unable to recover a region of
the planted size – what it declares is a larger overlapping region. That
is a property of the design, knowable before any replicate is run, and
this function is how to know it.

## Note

**Survival DGMs are not supported, and what they would need is this:**
the OC family's contract is `inherits(dgm, "glm_dgm")` with a `df_super`
carrying `flag_harm`.
[`setup_gbsg_dgm()`](https://larry-leon.github.io/forestsearch/reference/setup_gbsg_dgm.md)
produces a different object whose own generator is
`simulate_from_dgm(dgm, n, analysis_time, cens_adjust, seed)` – a
different signature, and it takes two censoring arguments that have no
GLM counterpart and that change the per-arm event counts this function
reports. Supporting it means a generator-dispatch layer plus passing
`analysis_time` and `cens_adjust` through; that is a separate change and
is deliberately not guessed at here.

## What is counted

Per `n`, over `n_rep` draws:

- `share_undeclarable`:

  share with \\\|Q\| \le\\ `n.min`. The comparison is `<=`, exactly as
  [`subgroup.search()`](https://larry-leon.github.io/forestsearch/reference/subgroup.search.md)
  tests it (`nx <= n.min` rejects), so a region must *exceed* `n.min` to
  be declarable.

- `share_under_events`:

  share with fewer than `d0.min` events in the control arm or `d1.min`
  in the treated arm. Binary and survival only – the floor is skipped
  entirely for continuous and count outcomes.

- `share_nonestimable`:

  share on which the estimand does not exist, under the same
  per-estimand condition the estimator boundary applies: OR needs all
  four cells \\\ge 1\\; RR, IRR and HR need an event in each arm; RD,
  IRD and MD have no condition and this share is always 0.

Together with per-arm cell summaries (mean, 5th and 95th percentile, and
minimum of each cell across replicates).

## Drawing through the DGM's own generator

Replicates are drawn with
[`simulate_from_glm_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_glm_dgm.md)
– the same function the campaign templates call – never a
re-implementation, so the population this reports on is the population a
campaign would run on.

This function **does not change the RNG kind**. It seeds each replicate
with [`set.seed()`](https://rdrr.io/r/base/Random.html) under whatever
generator is current and restores the caller's RNG state on exit. The
caller's own obligation is the mirror of that one: **build and calibrate
the DGM before any replicate switches the RNG kind.** Calibrating a DGM
after a switch to `"L'Ecuyer-CMRG"` yields a *different*
super-population, and every count here would then describe a population
no campaign ever used.

## See also

[`simulate_from_glm_dgm`](https://larry-leon.github.io/forestsearch/reference/simulate_from_glm_dgm.md),
[`fs_oc_predict`](https://larry-leon.github.io/forestsearch/reference/fs_oc_predict.md)

## Examples

``` r
# \donttest{
set.seed(1)
N <- 400
d <- data.frame(treat = rbinom(N, 1, 0.5),
                age   = round(rnorm(N, 55, 12)),
                sex   = factor(rbinom(N, 1, 0.5), levels = 0:1))
d$y <- rbinom(N, 1, plogis(-0.5 + 0.3 * d$treat))
dgm <- generate_glm_dgm(data = d, factor_vars = "sex",
                        continuous_vars = "age", outcome_var = "y",
                        treatment_var = "treat", outcome_type = "binary",
                        effect_measure = "OR", subgroup_vars = "age",
                        subgroup_cuts = list(age = list(type = "greater",
                                                        quantile = 0.55)),
                        k_inter = 1, n_super = 2000, seed = 8316951)
feas <- fs_dgm_feasibility(dgm, n = c(500, 1000), n_rep = 50)
feas
#> DGM feasibility -- can the planted region be declared at this n?
#>   outcome_type = binary; effect_measure = OR; 50 replicates per n
#>   floors read (not imposed): n.min = 60, d0.min = 10, d1.min = 10
#> 
#>        n size_mean size_q05 size_min   undeclarable   under_events  non-estimable
#>      500     223.0    206.3      197           0.0%          0.0%           0.0%
#>     1000     444.9    419.0      390           0.0%          0.0%           0.0%
#> 
#>   FEASIBLE: every undeclarable share is at or below the 5% tolerance.
feas$feasible
#> [1] TRUE
# }
```
