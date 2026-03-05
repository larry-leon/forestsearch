# Calibrate Censoring Adjustment to Match DGM Reference Distribution

Uses root-finding to select a value of `cens_adjust` for
[`simulate_from_dgm`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md)
such that a chosen censoring summary statistic in the simulated data
matches the corresponding statistic from the DGM reference data
(`dgm$df_super`).

## Usage

``` r
calibrate_cens_adjust(
  dgm,
  target = c("rate", "km_median"),
  n = 1000,
  rand_ratio = 1,
  analysis_time = 48,
  max_entry = 24,
  seed = 42,
  interval = c(-3, 3),
  tol = 1e-04,
  n_eval = 2000,
  verbose = TRUE,
  ...
)
```

## Arguments

- dgm:

  An `"aft_dgm_flex"` object from
  [`generate_aft_dgm_flex`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md).

- target:

  Character. Calibration target: `"rate"` (default) or `"km_median"`.

- n:

  Integer. Sample size passed to
  [`simulate_from_dgm`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md).
  Default `1000`.

- rand_ratio:

  Numeric. Randomisation ratio passed to
  [`simulate_from_dgm`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md).
  Default `1`.

- analysis_time:

  Numeric. Calendar analysis time passed to
  [`simulate_from_dgm`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md).
  Must be on the DGM time scale. Default `48`.

- max_entry:

  Numeric. Maximum staggered entry time passed to
  [`simulate_from_dgm`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md).
  Default `24`.

- seed:

  Integer. Base random seed. Each evaluation of the objective function
  uses this seed for reproducibility. Default `42`.

- interval:

  Numeric vector of length 2. Search interval for `cens_adjust` on the
  log scale. Default `c(-3, 3)` (corresponding roughly to a 20-fold
  decrease/increase in censoring times).

- tol:

  Numeric. Root-finding tolerance. Default `1e-4`.

- n_eval:

  Integer. Sample size used inside the objective function during
  root-finding. Smaller values are faster but noisier; increase for
  precision. Default `2000`.

- verbose:

  Logical. Print search progress and final result. Default `TRUE`.

- ...:

  Additional arguments passed to
  [`simulate_from_dgm`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md)
  (e.g. `strata_rand`, `time_eos`).

## Value

A named list with elements:

- `cens_adjust`:

  Calibrated `cens_adjust` value.

- `target`:

  Calibration target used.

- `ref_value`:

  Reference metric value from `dgm$df_super`.

- `sim_value`:

  Achieved metric value in simulated data at the calibrated
  `cens_adjust`.

- `residual`:

  Absolute difference between `sim_value` and `ref_value`.

- `iterations`:

  Number of `uniroot` iterations.

- `diagnostic`:

  Output of
  [`check_censoring_dgm`](https://larry-leon.github.io/forestsearch/reference/check_censoring_dgm.md)
  at the calibrated value (invisibly).

## Details

Two calibration targets are supported:

- `"rate"`:

  Overall censoring rate (proportion censored). Finds `cens_adjust` such
  that `mean(event_sim == 0)` in simulated data equals
  `mean(event == 0)` in `dgm$df_super`.

- `"km_median"`:

  KM-based median censoring time, estimated by reversing the event
  indicator so censored observations become the "event" of interest.
  Finds `cens_adjust` such that the simulated KM median matches the
  reference KM median.

### How the objective function works

At each candidate `cens_adjust` value, the objective function:

1.  Calls
    [`simulate_from_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md)
    with `n = n_eval` and the candidate `cens_adjust`.

2.  Calls
    [`check_censoring_dgm()`](https://larry-leon.github.io/forestsearch/reference/check_censoring_dgm.md)
    with `verbose = FALSE` to extract the target metric.

3.  Returns `sim_metric - ref_metric`.

`uniroot` finds the zero crossing, i.e. the `cens_adjust` at which
simulated and reference metrics are equal.

### Monotonicity

The objective is monotone in `cens_adjust` for both targets:

- Larger `cens_adjust` → longer censoring times → lower censoring rate
  and higher KM median.

- Smaller `cens_adjust` → shorter censoring times → higher censoring
  rate and lower KM median.

If `uniroot` fails (the target lies outside the search interval), the
boundary values are printed and a wider `interval` should be tried.

### Stochastic noise

Because the objective function involves simulation, there is Monte Carlo
noise. Setting a fixed `seed` and a sufficiently large `n_eval` (\>=
2000) reduces noise enough for reliable root-finding. The `tol` argument
controls the root-finding tolerance on the `cens_adjust` scale (not the
metric scale).

## See also

[`simulate_from_dgm`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md),
[`check_censoring_dgm`](https://larry-leon.github.io/forestsearch/reference/check_censoring_dgm.md),
[`generate_aft_dgm_flex`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md)

## Examples

``` r
if (FALSE) { # \dontrun{
library(survival)

# Build DGM on months scale
gbsg$time_months <- gbsg$rfstime / 30.4375

dgm <- generate_aft_dgm_flex(
  data            = gbsg,
  continuous_vars = c("age", "size", "nodes", "pgr", "er"),
  factor_vars     = c("meno", "grade"),
  outcome_var     = "time_months",
  event_var       = "status",
  treatment_var   = "hormon",
  subgroup_vars   = c("er", "meno"),
  subgroup_cuts   = list(er = 20, meno = 0)
)

# Calibrate so simulated censoring rate matches reference
cal_rate <- calibrate_cens_adjust(
  dgm           = dgm,
  target        = "rate",
  n             = 1000,
  analysis_time = 84,
  max_entry     = 24
)
cat("Calibrated cens_adjust (rate):", cal_rate$cens_adjust, "\n")

# Calibrate to KM median censoring time instead
cal_km <- calibrate_cens_adjust(
  dgm           = dgm,
  target        = "km_median",
  n             = 1000,
  analysis_time = 84,
  max_entry     = 24
)
cat("Calibrated cens_adjust (km_median):", cal_km$cens_adjust, "\n")

# Use calibrated value in simulation
sim <- simulate_from_dgm(
  dgm           = dgm,
  n             = 1000,
  analysis_time = 84,
  max_entry     = 24,
  cens_adjust   = cal_rate$cens_adjust,
  seed          = 123
)
mean(sim$event_sim)   # event rate
mean(sim$event_sim == 0)  # censoring rate — should match ref
} # }
```
