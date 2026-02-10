# MRCT Regional Subgroup Simulation

Simulates multi-regional clinical trials and evaluates ForestSearch
subgroup identification. Splits data by region into training and testing
populations, identifies subgroups using ForestSearch on training data,
and evaluates performance on the testing region.

## Usage

``` r
mrct_region_sims(
  dgm,
  n_sims,
  n_sample = NULL,
  region_var = "z_regA",
  sg_focus = "minSG",
  maxk = 1,
  hr.threshold = 0.9,
  hr.consistency = 0.8,
  pconsistency.threshold = 0.9,
  confounders.name = NULL,
  conf_force = NULL,
  analysis_time = 60,
  cens_adjust = 0,
  parallel_args = list(plan = "multisession", workers = NULL, show_message = TRUE),
  details = FALSE,
  seed = NULL
)
```

## Arguments

- dgm:

  Data generating mechanism object from
  [`generate_aft_dgm_flex`](https://github.com/larry-leon/forestsearch/reference/generate_aft_dgm_flex.md)

- n_sims:

  Integer. Number of simulations to run

- n_sample:

  Integer. Sample size per simulation. If NULL (default), uses the
  entire super-population from dgm

- region_var:

  Character. Name of the region indicator variable used to split data
  into training (region_var == 0) and testing (region_var == 1)
  populations. Default: "z_regA"

- sg_focus:

  Character. Subgroup selection criterion passed to
  [`forestsearch`](https://github.com/larry-leon/forestsearch/reference/forestsearch.md):
  "minSG", "hr", or "maxSG". Default: "minSG"

- maxk:

  Integer. Maximum number of factors in subgroup combinations (1 or 2).
  Default: 1

- hr.threshold:

  Numeric. Hazard ratio threshold for subgroup identification. Default:
  0.90

- hr.consistency:

  Numeric. Consistency threshold for hazard ratio. Default: 0.80

- pconsistency.threshold:

  Numeric. Probability threshold for consistency. Default: 0.90

- confounders.name:

  Character vector. Confounder variable names for ForestSearch. If NULL,
  automatically extracted from dgm

- conf_force:

  Character vector. Forced cuts to consider in ForestSearch. Default:
  c("z_age \<= 65", "z_bm \<= 0", "z_bm \<= 1", "z_bm \<= 2", "z_bm \<=
  5")

- analysis_time:

  Numeric. Time of analysis for administrative censoring. Default: 60

- cens_adjust:

  Numeric. Adjustment factor for censoring rate on log scale. Default: 0

- parallel_args:

  List. Parallel processing configuration with components:

  - plan: "multisession", "multicore", "callr", or "sequential"

  - workers: Number of workers (NULL for auto-detect)

  - show_message: Logical for progress messages

- details:

  Logical. Print detailed progress information. Default: FALSE

- seed:

  Integer. Base random seed for reproducibility. Default: NULL

## Value

A data.table with simulation results containing:

- sim:

  Simulation index

- n_itt:

  ITT sample size

- hr_itt:

  ITT hazard ratio (stratified if strat variable present)

- hr_ittX:

  ITT hazard ratio stratified by region

- n_train:

  Training (non-region A) sample size

- hr_train:

  Training population hazard ratio

- n_test:

  Testing (region A) sample size

- hr_test:

  Testing population hazard ratio

- any_found:

  Indicator: 1 if subgroup identified, 0 otherwise

- sg_found:

  Character description of identified subgroup

- n_sg:

  Subgroup sample size

- hr_sg:

  Subgroup hazard ratio

- POhr_sg:

  Potential outcome hazard ratio in subgroup

- prev_sg:

  Subgroup prevalence (proportion of testing population)

- hr_sg_null:

  Subgroup HR when found, NA otherwise

## Details

### Simulation Process

For each simulation:

1.  Sample from super-population using
    [`simulate_from_dgm`](https://github.com/larry-leon/forestsearch/reference/simulate_from_dgm.md)

2.  Split by region_var into training and testing populations

3.  Estimate HRs in ITT, training, and testing populations

4.  Run
    [`forestsearch`](https://github.com/larry-leon/forestsearch/reference/forestsearch.md)
    on training population

5.  Apply identified subgroup to testing population

6.  Calculate subgroup-specific estimates

### Region Variable

The region_var parameter is used ONLY for splitting data into
training/testing populations. It does not imply any prognostic effect.
To include prognostic confounder effects, specify them when creating the
DGM using
[`create_dgm_for_mrct`](https://github.com/larry-leon/forestsearch/reference/create_dgm_for_mrct.md)
or
[`generate_aft_dgm_flex`](https://github.com/larry-leon/forestsearch/reference/generate_aft_dgm_flex.md).

## See also

[`forestsearch`](https://github.com/larry-leon/forestsearch/reference/forestsearch.md)
for subgroup identification algorithm
[`generate_aft_dgm_flex`](https://github.com/larry-leon/forestsearch/reference/generate_aft_dgm_flex.md)
for DGM creation
[`simulate_from_dgm`](https://github.com/larry-leon/forestsearch/reference/simulate_from_dgm.md)
for data simulation
[`create_dgm_for_mrct`](https://github.com/larry-leon/forestsearch/reference/create_dgm_for_mrct.md)
for MRCT-specific DGM wrapper
[`summaryout_mrct`](https://github.com/larry-leon/forestsearch/reference/summaryout_mrct.md)
for summarizing simulation results

## Examples

``` r
if (FALSE) { # \dontrun{
# Create DGM for alternative hypothesis
dgm_alt <- create_dgm_for_mrct(
  df_case = df_case,
  model_type = "alt",
  log_hrs = log(c(3, 1.25, 0.50)),
  verbose = TRUE
)

# Run simulations
results <- mrct_region_sims(
  dgm = dgm_alt,
  n_sims = 100,
  region_var = "z_regA",
  sg_focus = "minSG",
  parallel_args = list(plan = "multisession", workers = 4),
  details = TRUE
)

# Summarize results
cat("Subgroup identification rate:", mean(results$any_found) * 100, "%\n")
} # }
```
