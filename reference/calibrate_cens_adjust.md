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
# \donttest{
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
#> 
#> === Subgroup Definitions ===
#>    er <= 20 
#>     Proportion in subgroup: 0.397 
#>    meno <= 0 
#>     Proportion in subgroup: 0.423 
#> 
#> Overall subgroup (all conditions met):
#>   Size: 134 out of 686 
#>   Proportion: 0.195 
#> 
#> === Model Parameters (AFT, log(T)) ===
#> Intercept (mu): 4.057 
#> Scale (tau): 0.721 
#> Treatment effect: 0.321 
#> Interaction effect: -0.313 
#> 
#> Overall subgroup in super-population (all conditions met):
#>   Proportion: 0.197 
#> 
#> === Hazard Ratios (super popln) ===
#> Overall HR: 0.768 
#> Causal AHR: 0.698 
#> Harm subgroup HR: 0.987 
#> Harm subgroup AHR: 0.99 
#> No-harm subgroup HR: 0.733 
#> No-harm subgroup AHR: 0.641 
#> 
#>  ====================================================================== 
#>  CENSORING MODEL SELECTION:
#>  MULTIPLE SURVIVAL REGRESSION MODEL COMPARISONS
#> ====================================================================== 
#> 
#> CONVERGENCE SUMMARY ( 4 / 4  converged):
#> -------------------------------------------------- 
#> 
#> Weibull (weibull):
#>   Status:  ✓ Converged 
#>   Iterations:  9 
#>   Log-likelihood:  -1838.404 
#>   Scale parameter:  0.418 
#> 
#> LogNormal (lognormal):
#>   Status:  ✓ Converged 
#>   Iterations:  4 
#>   Log-likelihood:  -2000.309 
#>   Scale parameter:  0.8707 
#> 
#> Weibull0 (weibull):
#>   Status:  ✓ Converged 
#>   Iterations:  7 
#>   Log-likelihood:  -1843.717 
#>   Scale parameter:  0.4246 
#> 
#> LogNormal0 (lognormal):
#>   Status:  ✓ Converged 
#>   Iterations:  6 
#>   Log-likelihood:  -2004.326 
#>   Scale parameter:  0.8768 
#> 
#> MODEL COMPARISON TABLE:
#> -------------------------------------------------- 
#>       Model     AIC delta_AIC AIC_weight     BIC delta_BIC
#>    Weibull0 3691.43      0.00      0.976 3700.50      0.00
#>     Weibull 3698.81      7.37      0.024 3748.65     48.15
#>  LogNormal0 4012.65    321.22      0.000 4021.71    321.22
#>   LogNormal 4022.62    331.18      0.000 4072.46    371.96
#> 
#> MODEL RANKINGS:
#> -------------------------------------------------- 
#> By AIC:  Weibull0 > Weibull > LogNormal0 > LogNormal 
#> By BIC:  Weibull0 > Weibull > LogNormal0 > LogNormal 
#> 
#> EVIDENCE ASSESSMENT:
#> -------------------------------------------------- 
#>   Strength:  Strong evidence 
#>   Details:  Top model has 97.6% of AIC weight 
#> 
#> FINAL RECOMMENDATION:
#> -------------------------------------------------- 
#>   Clear winner: Weibull0 (best by both AIC and BIC)
#>   Evidence: Strong evidence - Top model has 97.6% of AIC weight 
#> ====================================================================== 
#> 

# Calibrate so simulated censoring rate matches reference
cal_rate <- calibrate_cens_adjust(
  dgm           = dgm,
  target        = "rate",
  n             = 1000,
  analysis_time = 84,
  max_entry     = 24
)
#> 
#> -------------------------------------------------------
#>   calibrate_cens_adjust()
#> -------------------------------------------------------
#>   Target metric    : rate
#>   Reference value  : 0.5532
#>   Search interval  : [-3.00, 3.00]
#>   n_eval per call  : 2000
#>   Tolerance        : 1.0e-04
#> -------------------------------------------------------
#>   f(interval[1] = -3.00) = +0.4363
#>   f(interval[2] = 3.00) = -0.1597
#>   Searching...
#> 
#> =========================================================
#>   Censoring Diagnostic: DGM Reference vs Simulated
#> =========================================================
#>   DGM censoring type  : weibull
#>   DGM mu_cens         : 4.0312  [exp = 56.33 time units]
#>   DGM tau_cens        : 0.4246
#>   Reference n (super) : 5000
#>   Simulated n         : 1000
#> 
#> --- 1. Censoring Rates ---
#>    Group Ref_rate Sim_rate    Diff
#>  Overall    55.3%    55.1% -0.2 pp
#>    Arm 0    51.7%    49.0% -2.7 pp
#>    Arm 1    61.8%    61.2% -0.6 pp
#> 
#> --- 2. Censoring-Time Quantiles (censored subjects only) ---
#>     Ref: y[event == 0]  |  Sim: c_time[event_sim == 0]
#>  Quantile   Ref   Sim Ratio
#>       25% 28.19 30.78 1.092
#>       50% 47.05 46.19 0.982
#>       75% 61.04 61.66 1.010
#>       90% 71.00 71.27 1.004
#> 
#> --- 3. KM Median Censoring Time (reversed event indicator) ---
#>    Group Ref_median Sim_median Ratio
#>  Overall   54.04517   53.37382 0.988
#>    Arm 0   51.25257   53.37382 1.041
#>    Arm 1   57.69199   53.39228 0.925
#> 
#> --- 4. Implied Time-Scale Check ---
#>   exp(params$mu)          = 57.83  [median outcome time, DGM scale]
#>   exp(params$censoring$mu)= 56.33  [median censoring time, DGM scale]
#>   If these values are implausibly large relative to your
#>   analysis_time, the DGM was likely built on a different time
#>   scale (e.g. days vs months).
#> 
#> [OK] No substantial discrepancies detected.
#> =========================================================
#> 
#> 
#>   === Calibration Result ===
#>   cens_adjust      : 0.097831
#>   Reference rate      : 0.5532
#>   Simulated rate      : 0.5510
#>   Residual         : 0.0022
#>   uniroot iters    : 18
#> -------------------------------------------------------
#> 
cat("Calibrated cens_adjust (rate):", cal_rate$cens_adjust, "\n")
#> Calibrated cens_adjust (rate): 0.09783093 

# Calibrate to KM median censoring time instead
cal_km <- calibrate_cens_adjust(
  dgm           = dgm,
  target        = "km_median",
  n             = 1000,
  analysis_time = 84,
  max_entry     = 24
)
#> 
#> -------------------------------------------------------
#>   calibrate_cens_adjust()
#> -------------------------------------------------------
#>   Target metric    : km_median
#>   Reference value  : 54.0452
#>   Search interval  : [-3.00, 3.00]
#>   n_eval per call  : 2000
#>   Tolerance        : 1.0e-04
#> -------------------------------------------------------
#>   f(interval[1] = -3.00) = -51.6931
#>   f(interval[2] = 3.00) = +17.6392
#>   Searching...
#> 
#> =========================================================
#>   Censoring Diagnostic: DGM Reference vs Simulated
#> =========================================================
#>   DGM censoring type  : weibull
#>   DGM mu_cens         : 4.0312  [exp = 56.33 time units]
#>   DGM tau_cens        : 0.4246
#>   Reference n (super) : 5000
#>   Simulated n         : 1000
#> 
#> --- 1. Censoring Rates ---
#>    Group Ref_rate Sim_rate    Diff
#>  Overall    55.3%    54.8% -0.5 pp
#>    Arm 0    51.7%    48.6% -3.1 pp
#>    Arm 1    61.8%    61.0% -0.8 pp
#> 
#> --- 2. Censoring-Time Quantiles (censored subjects only) ---
#>     Ref: y[event == 0]  |  Sim: c_time[event_sim == 0]
#>  Quantile   Ref   Sim Ratio
#>       25% 28.19 31.42 1.114
#>       50% 47.05 47.20 1.003
#>       75% 61.04 62.37 1.022
#>       90% 71.00 71.53 1.008
#> 
#> --- 3. KM Median Censoring Time (reversed event indicator) ---
#>    Group Ref_median Sim_median Ratio
#>  Overall   54.04517   54.42848 1.007
#>    Arm 0   51.25257   54.46035 1.063
#>    Arm 1   57.69199   54.42848 0.943
#> 
#> --- 4. Implied Time-Scale Check ---
#>   exp(params$mu)          = 57.83  [median outcome time, DGM scale]
#>   exp(params$censoring$mu)= 56.33  [median censoring time, DGM scale]
#>   If these values are implausibly large relative to your
#>   analysis_time, the DGM was likely built on a different time
#>   scale (e.g. days vs months).
#> 
#> [OK] No substantial discrepancies detected.
#> =========================================================
#> 
#> 
#>   === Calibration Result ===
#>   cens_adjust      : 0.117984
#>   Reference km_median : 54.0452
#>   Simulated km_median : 54.4285
#>   Residual         : 0.3833
#>   uniroot iters    : 8
#> -------------------------------------------------------
#> 
cat("Calibrated cens_adjust (km_median):", cal_km$cens_adjust, "\n")
#> Calibrated cens_adjust (km_median): 0.1179835 

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
#> [1] 0.458
mean(sim$event_sim == 0)  # censoring rate — should match ref
#> [1] 0.542
# }
```
