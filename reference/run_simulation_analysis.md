# Run One Simulation Replicate

General replacement for the legacy `run_simulation_analysis()` that was
coupled to
[`simulate_from_gbsg_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_gbsg_dgm.md)
and GBSG-specific column names. This version calls
[`simulate_from_dgm`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md)
and accepts explicit column-name parameters, making it applicable to any
DGM built with
[`generate_aft_dgm_flex`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md).

## Usage

``` r
run_simulation_analysis(
  sim_id,
  dgm,
  n_sample,
  analysis_time = Inf,
  cens_adjust = 0,
  max_follow = NULL,
  muC_adj = NULL,
  confounders_base = c("v1", "v2", "v3", "v4", "v5", "v6", "v7"),
  n_add_noise = 0L,
  outcome_name = "y_sim",
  event_name = "event_sim",
  treat_name = "treat_sim",
  harm_col = "flag_harm",
  run_fs = TRUE,
  run_fs_grf = TRUE,
  run_grf = TRUE,
  fs_params = list(),
  grf_params = list(),
  cox_formula = NULL,
  cox_formula_adj = NULL,
  n_sims_total = NULL,
  seed_base = 8316951L,
  verbose = FALSE,
  verbose_n = NULL,
  debug = FALSE
)
```

## Arguments

- sim_id:

  Integer. Simulation replicate index (used as seed offset).

- dgm:

  An `"aft_dgm_flex"` object from
  [`generate_aft_dgm_flex`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md)
  or
  [`setup_gbsg_dgm`](https://larry-leon.github.io/forestsearch/reference/setup_gbsg_dgm.md).

- n_sample:

  Integer. Per-replicate sample size.

- analysis_time:

  Numeric. Calendar time of analysis on the DGM time scale. Use `Inf`
  (default) for no administrative censoring — equivalent to the legacy
  `max_follow = Inf`.

- cens_adjust:

  Numeric. Log-scale shift to censoring times passed to
  `simulate_from_dgm(cens_adjust = ...)`. Replaces legacy `muC_adj`.
  Default `0`.

- max_follow:

  **Deprecated.** Use `analysis_time` instead. If supplied, its value is
  forwarded to `analysis_time` with a warning. Retained for backward
  compatibility with legacy scripts.

- muC_adj:

  **Deprecated.** Use `cens_adjust` instead. If supplied, its value is
  forwarded to `cens_adjust` with a warning. Retained for backward
  compatibility with legacy scripts.

- confounders_base:

  Character vector of base confounder names.

- n_add_noise:

  Integer. Number of independent N(0,1) noise variables to append.
  Default `0L`.

- outcome_name:

  Name of the observed time column in simulated data. Default `"y_sim"`.

- event_name:

  Name of the event indicator column. Default `"event_sim"`.

- treat_name:

  Name of the treatment column. Default `"treat_sim"`.

- harm_col:

  Name of the true-subgroup indicator column. Default `"flag_harm"`.

- run_fs:

  Logical. Run ForestSearch (LASSO). Default `TRUE`.

- run_fs_grf:

  Logical. Run ForestSearch (LASSO + GRF). Default `TRUE`.

- run_grf:

  Logical. Run standalone GRF. Default `TRUE`.

- fs_params:

  Named list of ForestSearch parameter overrides.

- grf_params:

  Named list of GRF parameter overrides.

- cox_formula:

  Optional Cox formula for unadjusted ITT.

- cox_formula_adj:

  Optional adjusted Cox formula.

- n_sims_total:

  Integer. Total simulations (for progress messages).

- seed_base:

  Integer. Base seed; replicate seed = `seed_base + sim_id`. Default
  `8316951L`.

- verbose:

  Logical. Print progress messages. Default `FALSE`.

- verbose_n:

  Integer. If set, only print for `sim_id <= verbose_n`. Default `NULL`.

- debug:

  Logical. Print detailed debug output. Default `FALSE`.

## Value

A `data.table` with one row per analysis method, containing subgroup
size, HR, AHR, CDE, and classification metrics.

## See also

[`simulate_from_dgm`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md),
[`generate_aft_dgm_flex`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md),
[`setup_gbsg_dgm`](https://larry-leon.github.io/forestsearch/reference/setup_gbsg_dgm.md)
