# Run Null Calibration Simulation

Runs ForestSearch under the null hypothesis at multiple sample sizes to
calibrate \\L\_{\text{eff}}\\, then returns a calibration object for
predicting corrected FPR.

## Usage

``` r
run_null_calibration(
  sim_null_fn,
  fs_args,
  d_eff_fn,
  N_values = c(200, 500, 1000),
  c1 = 1.25,
  c2 = 1,
  n_sims = 100L,
  prop_harm = 0.3,
  n_min = 60L,
  seed_base = 42L,
  verbose = TRUE
)
```

## Arguments

- sim_null_fn:

  Function. Generates one H0 dataset. Signature: `sim_null_fn(n, seed)`
  returning a data.frame suitable for
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md).

- fs_args:

  Named list. Arguments passed to
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  (excluding `df.analysis`, `hr.threshold`, `hr.consistency`, and
  `seedit`, which are set by the wrapper).

- d_eff_fn:

  Function. Computes \\d\_{\text{eff}}\\ for a subgroup. Signature:
  `d_eff_fn(n_sg)` returning numeric. Example:
  `function(n) d_eff_binary(n, 0.30)`.

- N_values:

  Integer vector. Sample sizes to simulate. At least 2 required.
  Default: `c(200, 500, 1000)`.

- c1:

  Numeric. Screening threshold for calibration. Default: 1.25.

- c2:

  Numeric. Consistency threshold for calibration. Default: 1.0.

- n_sims:

  Integer. Simulations per sample size. Default: 100.

- prop_harm:

  Numeric. Expected proportion in harm subgroup (for computing d_eff).
  Default: 0.30.

- n_min:

  Integer. Minimum subgroup size. Default: 60.

- seed_base:

  Integer. Base seed. Default: 42.

- verbose:

  Logical. Print progress. Default: TRUE.

## Value

An `"leff_calibration"` object.

## Details

For each N in `N_values`, the function:

1.  Generates `n_sims` datasets via `sim_null_fn()`

2.  Runs
    [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
    on each

3.  Computes the simulated FPR = proportion with any H found

4.  Computes \\P_1\\ via the GLM approximation

5.  Calls
    [`calibrate_L_eff()`](https://larry-leon.github.io/forestsearch/reference/calibrate_L_eff.md)
    to fit the power-law model

ForestSearch is run with `parallel_args = list(plan = "sequential")` to
avoid nested parallelism. The outer simulation loop is sequential. For
faster execution, wrap in `foreach %dofuture%` externally.

## Examples

``` r
# \donttest{
# Binary DGM under H0
sim_null <- function(n, seed) {
  set.seed(seed)
  data.frame(
    id = seq_len(n),
    treat = rbinom(n, 1, 0.5),
    bm1 = as.factor(rbinom(n, 1, 0.70)),
    bm2 = as.factor(rbinom(n, 1, 0.50)),
    age = round(rnorm(n, 55, 10)),
    ecog = as.factor(sample(0:1, n, TRUE, c(0.6, 0.4))),
    progressed = rbinom(n, 1, 0.30)
  )
}

fs_args <- list(
  confounders.name = c("bm1", "bm2", "age", "ecog"),
  outcome.name = "progressed",
  event.name = "progressed",
  treat.name = "treat",
  id.name = "id",
  outcome_type = "binary",
  effect_measure = "OR",
  adverse_outcome = TRUE,
  pconsistency.threshold = 0.90,
  fs.splits = 200,
  n.min = 60,
  maxk = 2,
  use_lasso = TRUE,
  use_grf = TRUE,
  use_twostage = TRUE,
  is.RCT = TRUE,
  details = FALSE,
  plot.sg = FALSE,
  parallel_args = list(plan = "sequential",
                       workers = 1,
                       show_message = FALSE)
)

cal <- run_null_calibration(
  sim_null_fn = sim_null,
  fs_args = fs_args,
  d_eff_fn = function(n) d_eff_binary(n, 0.30),
  N_values = c(200, 500, 1000),
  n_sims = 50,
  verbose = TRUE
)
#> Null calibration: 3 N values x 50 sims = 150 tasks
#> Warning: Error in get_FSdata: get_FSdata: no candidate cut expressions could be constructed from the supplied inputs.  Upstream source counts: GRF cuts (0), DINA cuts (0), LASSO selected (1), conf_force (0), conf.cont_medians (0), conf.categorical (0).  Check that: (1) GRF is succeeding for this outcome type (rerun forestsearch() with details = TRUE to inspect); (2) LASSO is not shrinking all covariates to zero (try use_lasso = FALSE); (3) candidate covariates are numeric or factor and have sufficient variability; (4) GRF-related thresholds (dmin.grf, frac.tau, vi.grf.min) are appropriate for the outcome type and sample size (frac.tau is survival-specific and may be ignored for GLM).
#> Warning: FSdata failure - returning NULL result
#> Warning: Error in get_FSdata: get_FSdata: no candidate cut expressions could be constructed from the supplied inputs.  Upstream source counts: GRF cuts (0), DINA cuts (0), LASSO selected (1), conf_force (0), conf.cont_medians (0), conf.categorical (0).  Check that: (1) GRF is succeeding for this outcome type (rerun forestsearch() with details = TRUE to inspect); (2) LASSO is not shrinking all covariates to zero (try use_lasso = FALSE); (3) candidate covariates are numeric or factor and have sufficient variability; (4) GRF-related thresholds (dmin.grf, frac.tau, vi.grf.min) are appropriate for the outcome type and sample size (frac.tau is survival-specific and may be ignored for GLM).
#> Warning: FSdata failure - returning NULL result
#>   N =  200: FPR = 0.240 (12/50) | P1 = 0.227 | 17.5 sec
#> Warning: Error in get_FSdata: get_FSdata: no candidate cut expressions could be constructed from the supplied inputs.  Upstream source counts: GRF cuts (0), DINA cuts (0), LASSO selected (1), conf_force (0), conf.cont_medians (0), conf.categorical (0).  Check that: (1) GRF is succeeding for this outcome type (rerun forestsearch() with details = TRUE to inspect); (2) LASSO is not shrinking all covariates to zero (try use_lasso = FALSE); (3) candidate covariates are numeric or factor and have sufficient variability; (4) GRF-related thresholds (dmin.grf, frac.tau, vi.grf.min) are appropriate for the outcome type and sample size (frac.tau is survival-specific and may be ignored for GLM).
#> Warning: FSdata failure - returning NULL result
#>   N =  500: FPR = 0.600 (30/50) | P1 = 0.195 | 30.0 sec
#>   N = 1000: FPR = 0.640 (32/50) | P1 = 0.153 | 47.9 sec
#> 
#> L_eff Calibration
#>   Model: L_eff = 0.313 * (N / 60)^1.109
#>   R-squared: 0.945
#>   Calibrated at N = 200, 500, 1000
#> 
#>   Implied L_eff:
#>     N =  200: L_eff = 1.07 (fitted: 1.19)
#>     N =  500: L_eff = 4.22 (fitted: 3.29)
#>     N = 1000: L_eff = 6.15 (fitted: 7.09)
print(cal)
#> L_eff Calibration
#>   Model: L_eff = 0.313 * (N / 60)^1.109
#>   R-squared: 0.945
#>   Calibrated at N = 200, 500, 1000
#> 
#>   Implied L_eff:
#>     N =  200: L_eff = 1.07 (fitted: 1.19)
#>     N =  500: L_eff = 4.22 (fitted: 3.29)
#>     N = 1000: L_eff = 6.15 (fitted: 7.09)
# }
```
