# Run Single Simulation Analysis

Executes ForestSearch and/or GRF analysis on a single simulated dataset.
This is the core function called within a simulation loop.

## Usage

``` r
run_simulation_analysis(
  sim_id,
  dgm,
  n_sample,
  max_follow = Inf,
  muC_adj = 0,
  confounders_base = c("v1", "v2", "v3", "v4", "v5", "v6", "v7"),
  n_add_noise = 0L,
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

  Integer. Simulation index for seed offset and tracking

- dgm:

  A DGM object from
  [`create_gbsg_dgm`](https://github.com/larry-leon/forestsearch/reference/create_gbsg_dgm.md)
  or similar

- n_sample:

  Integer. Sample size for simulation

- max_follow:

  Numeric. Maximum follow-up time. Default: Inf

- muC_adj:

  Numeric. Censoring adjustment. Default: 0

- confounders_base:

  Character vector. Base confounder names

- n_add_noise:

  Integer. Number of noise variables to add. Default: 0

- run_fs:

  Logical. Run ForestSearch with LASSO variable selection. Default:
  TRUE. Analysis label: "FS"

- run_fs_grf:

  Logical. Run ForestSearch with LASSO + GRF variable selection.
  Default: TRUE. Analysis label: "FSlg"

- run_grf:

  Logical. Run standalone GRF analysis (grf.subg.harm.survival).
  Default: TRUE. Analysis label: "GRF"

- fs_params:

  List. ForestSearch parameters (overrides all defaults including
  use_lasso/use_grf). User-provided values take precedence over
  analysis-type defaults. For example,
  `fs_params = list(hr.threshold = 1.5, use_twostage = TRUE)` will apply
  to both FS and FSlg analyses.

- grf_params:

  List. GRF parameters for standalone GRF analysis (overrides defaults).
  Accepts all parameters for
  [`grf.subg.harm.survival()`](https://github.com/larry-leon/forestsearch/reference/grf.subg.harm.survival.md):
  n.min, dmin.grf, frac.tau, maxdepth, RCT, sg.criterion, seedit,
  outcome.name, event.name, treat.name, id.name. User-provided values
  take precedence over defaults.

- cox_formula:

  Formula. Cox model formula for estimation

- cox_formula_adj:

  Formula. Adjusted Cox model formula

- n_sims_total:

  Integer. Total simulations (for progress display)

- seed_base:

  Integer. Base random seed. Default: 8316951

- verbose:

  Logical. Print progress. Default: FALSE

- verbose_n:

  Integer. Only print verbose output for first N simulations. Default:
  NULL (print for all simulations when verbose = TRUE)

- debug:

  Logical. Print detailed debugging information. Default: FALSE

## Value

A data.table with analysis results for all requested methods, including
both HR and AHR metrics. Contains columns:

- sim:

  Simulation ID

- sizeH_true, propH_true:

  True harm subgroup size/proportion in sample

- analysis:

  Analysis method: "FS", "FSlg", or "GRF"

- any.H:

  1 if subgroup identified, 0 otherwise

- size.H, size.Hc:

  Size of identified H and complement

- hr.H.true, hr.H.hat:

  True and estimated HR in identified H

- hr.Hc.true, hr.Hc.hat:

  True and estimated HR in identified Hc

- ahr.H.true, ahr.H.hat:

  True and estimated AHR in identified H

- sens, spec, ppv, npv:

  Classification metrics

## Details

Aligned with create_gbsg_dgm() and generate_aft_dgm_flex() output
structures.

### Analysis Methods

The function can run up to three analysis types:

- **FS**: ForestSearch with LASSO variable selection only (default:
  use_lasso = TRUE, use_grf = FALSE)

- **FSlg**: ForestSearch with LASSO + GRF variable selection (default:
  use_lasso = TRUE, use_grf = TRUE)

- **GRF**: Standalone GRF-based subgroup identification using
  grf.subg.harm.survival()

### Parameter Merging Order

Parameters are merged in the following order (later values override
earlier):

1.  [`default_fs_params()`](https://github.com/larry-leon/forestsearch/reference/default_fs_params.md) -
    package defaults

2.  Analysis-type-specific defaults (use_lasso/use_grf for FS vs FSlg)

3.  User's `fs_params` - **final authority**

This means if you pass `fs_params = list(use_grf = TRUE)`, it will
override the FS analysis default of `use_grf = FALSE`.

### Output Column Naming Convention

The output distinguishes between:

- **True subgroup from DGM**: sizeH_true, propH_true (known from data
  generation)

- **Identified subgroup**: size.H, hr.H.hat (estimated by analysis)

## Examples

``` r
if (FALSE) { # \dontrun{
# Create DGM (aligned version)
dgm <- create_gbsg_dgm(model = "alt", k_inter = 2, verbose = TRUE)

# Run single simulation with LASSO only
result <- run_simulation_analysis(
  sim_id = 1,
  dgm = dgm,
  n_sample = 500,
  confounders_base = c("v1", "v2", "v3", "v4", "v5", "v6", "v7"),
  run_fs = TRUE,      # LASSO only
  run_fs_grf = FALSE, # Skip LASSO+GRF
  run_grf = FALSE,    # Skip standalone GRF
  verbose = TRUE
)

# Run all three analysis types
result_all <- run_simulation_analysis(
  sim_id = 1,
  dgm = dgm,
  n_sample = 500,
  run_fs = TRUE,
  run_fs_grf = TRUE,
  run_grf = TRUE,
  verbose = TRUE
)
# result_all has 3 rows: one for FS, one for FSlg, one for GRF

# With use_twostage = TRUE for faster analysis
result_fast <- run_simulation_analysis(
  sim_id = 1,
  dgm = dgm,
  n_sample = 500,
  fs_params = list(use_twostage = TRUE),
  verbose = TRUE
)
} # }
```
