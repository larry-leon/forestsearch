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
  fs_args = list(),
  sim_args = list(rand_ratio = 1, draw_treatment = TRUE),
  analysis_time = 60,
  cens_adjust = 0,
  parallel_args = list(plan = "multisession", workers = NULL, show_message = TRUE),
  details = FALSE,
  verbose_n_sims = 2L,
  seed = NULL,
  record = c("point", "full")
)
```

## Arguments

- dgm:

  Data generating mechanism object from
  [`generate_aft_dgm_flex`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md)

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
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md):
  "minSG", "hr", or "maxSG". Default: "minSG"

- maxk:

  Integer. Maximum number of factors in subgroup combinations (1 or 2).
  Default: 1

- hr.threshold:

  Numeric. `c1`, the screening threshold on a candidate subgroup's own
  hazard ratio, passed to
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md).
  Default: 0.90

- hr.consistency:

  Numeric. `c2`, **an effect threshold** – the minimum hazard ratio
  required in each half of a split. Not a proportion and not the
  consistency rate. Default: 0.80

- pconsistency.threshold:

  Numeric. `p*`, **a proportion** in `[0, 1]` – the fraction of splits
  that must clear `hr.consistency` (`c2`). Not an effect threshold.
  Default: 0.90

- confounders.name:

  Character vector. Confounder variable names for ForestSearch. If NULL,
  automatically extracted from dgm: the columns of `dgm$df_super` whose
  names start with `z_`, `ecog` or `strat`, except that, for a mechanism
  built with a spline (`dgm$model_params$spline_info`), the three
  internal columns of the spline variable `<var>` are left out by their
  exact names `<var>_treat`, `<var>_k` and `<var>_k_treat` (for
  [`create_dgm_for_mrct`](https://larry-leon.github.io/forestsearch/reference/create_dgm_for_mrct.md):
  `z_bm_treat`, `z_bm_k`, `z_bm_k_treat`). A covariate such as
  `z_prior_treat` and the region indicator stay in the set.

- conf_force:

  Character vector. Forced cuts to consider in ForestSearch. Default:
  NULL, which resolves to c("z_age \<= 65", "z_bm \<= 0", "z_bm \<= 1",
  "z_bm \<= 2", "z_bm \<= 5")

- fs_args:

  Named list. Additional arguments passed directly to
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  inside each simulation replicate. Use this to control parameters not
  exposed by `mrct_region_sims` (e.g., `use_grf`, `use_lasso`,
  `cut_type`, `d0.min`, `d1.min`, `n.min`, `max_subgroups_search`,
  `use_twostage`, `twostage_args`). Parameters already in the
  `mrct_region_sims` signature (`hr.threshold`, `hr.consistency`,
  `pconsistency.threshold`, `sg_focus`, `maxk`, `confounders.name`,
  `conf_force`) take precedence over values in `fs_args`. Default:
  list() (uses forestsearch defaults)

- sim_args:

  Named list. Additional arguments passed to
  [`simulate_from_dgm`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md)
  inside each replicate (e.g., `rand_ratio`, `draw_treatment`,
  `max_entry`, or `entry_delay = list(var = "z_region", delay = a)` for
  a region enrolling `a` time units later, redrawn every replicate). If
  `dgm$df_super` has a column `entrytime`, it is passed as
  `entry_var = "entrytime"` (one fixed entry time per super-population
  row) and takes precedence over `entry_delay`. Parameters already in
  the `mrct_region_sims` signature (`analysis_time`, `cens_adjust`) take
  precedence. Default: list(rand_ratio = 1, draw_treatment = TRUE)

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

- verbose_n_sims:

  Integer. When `details = TRUE`, print full ForestSearch diagnostics
  (including internal output) for only the first `verbose_n_sims`
  simulation replicates. Set to 0 to suppress per-sim output, or `Inf`
  to print all. Default: 2

- seed:

  Integer. Base random seed for reproducibility. Default: NULL

- record:

  Character. `"point"` (default) returns the point estimates only,
  exactly as before this argument existed. `"full"` adds, from the same
  simulated trials (nothing is redrawn), the per-replicate record
  described under *Full record* below.

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

  Subgroup hazard ratio in testing population

- POhr_sg:

  Potential outcome hazard ratio in subgroup (testing)

- prev_sg:

  Subgroup prevalence (proportion of testing population)

- n_sg_train:

  Subgroup sample size in training population

- hr_sg_train:

  Subgroup hazard ratio in training population

- POhr_sg_train:

  Potential outcome hazard ratio in subgroup (training)

- hr_sg_null:

  Subgroup HR when found, NA otherwise

With `record = "full"`, the same table followed by the columns
`TP_nonregion`, `FP_nonregion`, `FN_nonregion`, `TN_nonregion`,
`TP_region`, `FP_region`, `FN_region`, `TN_region`, then the
multiplier-resampling columns `mr_ci_method`, `mr_debiased_H`,
`mr_ij_lower_H`, `mr_ij_upper_H`, `mr_field_lower_1s_H`,
`mr_debiased_Hc`, `mr_ij_lower_Hc`, `mr_ij_upper_Hc`,
`mr_field_upper_1s_Hc`, `mr_field_upper_1s_s_Hc`, `mr_bonf_lower_H`,
`mr_bonf_upper_Hc`, `mr_p_hat_H`, `mr_selection_rate`, and the attribute
`"estimators"` (see *Full record*).

## Details

### Simulation Process

For each simulation:

1.  Sample from super-population using
    [`simulate_from_dgm`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md)

2.  Split by region_var into training and testing populations

3.  Estimate HRs in ITT, training, and testing populations

4.  Run
    [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
    on training population

5.  Apply identified subgroup to testing population

6.  Calculate subgroup-specific estimates

### Region Variable

The region_var parameter is used ONLY for splitting data into
training/testing populations. It does not imply any prognostic effect.
To include prognostic confounder effects, specify them when creating the
DGM using
[`create_dgm_for_mrct`](https://larry-leon.github.io/forestsearch/reference/create_dgm_for_mrct.md)
or
[`generate_aft_dgm_flex`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md).

### Full record

With `record = "full"` the returned table is the `"point"` table with
twenty-two columns added (eight counts, then fourteen
multiplier-resampling columns) and a second table attached as an
attribute, `est <- attr(res, "estimators")`; `res$hr_itt`, `nrow(res)`
and every `"point"` column are unchanged. (A data.table subset such as
`res[order(sim)]` drops the attribute: take it first.)

`estimators` has one row per trial and patient set: `itt_all` (all
patients), `itt_nonregion` (`region_var == 0`), `itt_region`
(`region_var == 1`), `region_Ghat` (the region outside the declared harm
subgroup \\\hat H\\; the whole region in an undeclared trial),
`region_Gstar` (the region with `flag_harm == 0`) and `Hstar_nonregion`
(the non-region patients with `flag_harm == 1`). Columns: `sim`,
`estimator`, `n`, `n_exp`, `events_exp`, `events_ctl`; the unadjusted
two-sample Cox fit `Surv(y_sim, event_sim) ~ treat_sim` on the set,
`hr`, `lo`, `hi` (Wald 95\\ `se_loghr` (a fit with an arm without events
is kept; a fit that errors, e.g. on an empty set, gives `NA`); the
targets on the same patients, `ahr` = `exp(mean(loghr_po))` and `cde` =
`mean(exp(theta_1)) / mean(exp(theta_0))` (`NA` when the simulated trial
lacks those columns); and `declared` (0/1, the trial's `any_found`).
`itt_all` is unadjusted, so it equals `hr_itt` only when the trial has
no `strat` column.

The added columns of the per-trial table are the 2 x 2 counts of \\\hat
H\\ membership against `flag_harm == 1` over the non-region patients and
over the region patients. \\\hat H\\ is evaluated from the selected
definition with
[`get_dfpred`](https://larry-leon.github.io/forestsearch/reference/get_dfpred.md)
(the route of `df.test`); it is empty in an undeclared trial, whose
counts are kept.

When the simulated trial has no `flag_harm` column, the `region_Gstar`
and `Hstar_nonregion` rows are omitted and the counts are `NA`, with a
message.

After the counts come fourteen columns read from the trial's
[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
result `fs$mr_inference` (multiplier resampling,
`fs_args = list(mr_inference = TRUE, ...)`; see
[`fs_mr_inference`](https://larry-leon.github.io/forestsearch/reference/fs_mr_inference.md)),
on the effect (HR) scale:

- `mr_ci_method`:

  `ci_method`.

- `mr_debiased_H`, `mr_ij_lower_H`, `mr_ij_upper_H`:

  `debiased$est`, `$lower`, `$upper`: the de-biased estimate on \\\hat
  H\\ and its interval (IJ SE under `"ij"` and `"field"`, robust SE
  under `"wald"`).

- `mr_field_lower_1s_H`:

  `field$lower_1s`.

- `mr_debiased_Hc`, `mr_ij_lower_Hc`, `mr_ij_upper_Hc`:

  `complement$debiased$est`, `$lower`, `$upper`.

- `mr_field_upper_1s_Hc`, `mr_field_upper_1s_s_Hc`:

  `field$complement$upper_1s` and its studentized companion
  `field$complement$upper_1s_s`.

- `mr_bonf_lower_H`, `mr_bonf_upper_Hc`:

  `field$joint$bonf_lower_H`, `field$joint$bonf_upper_Hc`.

- `mr_p_hat_H`:

  `reselection$p_hat[selected_index]`.

- `mr_selection_rate`:

  `selection_rate`.

Each is `NA` when its element is absent: every one in a trial without
`mr_inference` or without a declared subgroup (nothing to de-bias), or
with the short MR return; the complement columns without
`include_complement = TRUE`; the field columns unless
`ci_method = "field"`; the complement field and Bonferroni columns
without `field_complement = TRUE`; `mr_p_hat_H` without
`return_reselection = TRUE`. Under `record = "point"` none of these
columns is added.

## See also

[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
for subgroup identification algorithm
[`generate_aft_dgm_flex`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md)
for DGM creation
[`simulate_from_dgm`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md)
for data simulation
[`create_dgm_for_mrct`](https://larry-leon.github.io/forestsearch/reference/create_dgm_for_mrct.md)
for MRCT-specific DGM wrapper
[`summaryout_mrct`](https://larry-leon.github.io/forestsearch/reference/summaryout_mrct.md)
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
