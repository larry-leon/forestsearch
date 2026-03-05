# ForestSearch: Exploratory Subgroup Identification

Identifies subgroups with differential treatment effects in clinical
trials using a combination of Generalized Random Forests (GRF), LASSO
variable selection, and exhaustive combinatorial search with
split-sample validation.

## Usage

``` r
forestsearch(
  df.analysis,
  outcome.name = "tte",
  event.name = "event",
  treat.name = "treat",
  id.name = "id",
  potentialOutcome.name = NULL,
  flag_harm.name = NULL,
  confounders.name = NULL,
  parallel_args = list(plan = "callr", workers = 6, show_message = TRUE),
  df.predict = NULL,
  df.test = NULL,
  is.RCT = TRUE,
  seedit = 8316951,
  est.scale = "hr",
  use_lasso = TRUE,
  use_grf = TRUE,
  grf_res = NULL,
  grf_cuts = NULL,
  max_n_confounders = 1000,
  grf_depth = 2,
  dmin.grf = 12,
  frac.tau = 0.6,
  return_selected_cuts_only = TRUE,
  conf_force = NULL,
  defaultcut_names = NULL,
  cut_type = "default",
  exclude_cuts = NULL,
  replace_med_grf = FALSE,
  cont.cutoff = 4,
  conf.cont_medians = NULL,
  conf.cont_medians_force = NULL,
  n.min = 60,
  hr.threshold = 1.25,
  hr.consistency = 1,
  sg_focus = "hr",
  fs.splits = 1000,
  m1.threshold = Inf,
  pconsistency.threshold = 0.9,
  stop_threshold = 0.95,
  showten_subgroups = FALSE,
  d0.min = 12,
  d1.min = 12,
  max.minutes = 3,
  minp = 0.025,
  details = FALSE,
  maxk = 2,
  by.risk = 12,
  plot.sg = FALSE,
  plot.grf = FALSE,
  max_subgroups_search = 10,
  vi.grf.min = -0.2,
  use_twostage = TRUE,
  twostage_args = list()
)
```

## Arguments

- df.analysis:

  Data frame. Analysis dataset with required columns.

- outcome.name:

  Character. Name of time-to-event outcome variable. Default "tte".

- event.name:

  Character. Name of event indicator (1=event, 0=censored). Default
  "event".

- treat.name:

  Character. Name of treatment variable (1=treatment, 0=control).
  Default "treat".

- id.name:

  Character. Name of subject ID variable. Default "id".

- potentialOutcome.name:

  Character. Name of potential outcome variable (optional).

- flag_harm.name:

  Character. Name of true harm flag for simulations (optional).

- confounders.name:

  Character vector. Names of candidate subgroup-defining variables.

- parallel_args:

  List. Parallel processing configuration:

  plan

  :   Character. One of "multisession", "multicore", "callr",
      "sequential"

  workers

  :   Integer. Number of parallel workers

  show_message

  :   Logical. Show parallel setup messages

- df.predict:

  Data frame. Prediction dataset (optional).

- df.test:

  Data frame. Test dataset (optional).

- is.RCT:

  Logical. Is this a randomized controlled trial? Default TRUE.

- seedit:

  Integer. Random seed. Default 8316951.

- est.scale:

  Character. Estimation scale ("hr" or "rmst"). Default "hr".

- use_lasso:

  Logical. Use LASSO for variable selection. Default TRUE.

- use_grf:

  Logical. Use GRF for variable importance. Default TRUE.

- grf_res:

  GRF results object (optional, for reuse).

- grf_cuts:

  List. Custom GRF cut points (optional).

- max_n_confounders:

  Integer. Maximum confounders to consider. Default 1000.

- grf_depth:

  Integer. GRF tree depth. Default 2.

- dmin.grf:

  Integer. Minimum events for GRF. Default 12.

- frac.tau:

  Numeric. Fraction of tau for RMST. Default 0.6.

- return_selected_cuts_only:

  Logical. If TRUE (default), GRF returns only cuts from the tree depth
  that identified the selected subgroup meeting `dmin.grf`. If FALSE
  returns all cuts from all fitted trees (depths 1 through `grf_depth`).
  See
  [`grf.subg.harm.survival`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.survival.md)
  for details.

- conf_force:

  Character vector. Variables to force include (optional).

- defaultcut_names:

  Character vector. Default cut variable names (optional).

- cut_type:

  Character. Cut type ("default" or "custom"). Default "default".

- exclude_cuts:

  Character vector. Variables to exclude from cutting (optional).

- replace_med_grf:

  Logical. Replace median with GRF cuts. Default FALSE.

- cont.cutoff:

  Integer. Cutoff for continuous vs categorical. Default 4.

- conf.cont_medians:

  Named numeric vector. Median values for continuous variables
  (optional).

- conf.cont_medians_force:

  Named numeric vector. Forced median values (optional).

- n.min:

  Integer. Minimum subgroup size. Default 60.

- hr.threshold:

  Numeric. Minimum HR for candidate subgroups. Default 1.25.

- hr.consistency:

  Numeric. Minimum HR for consistency validation. Default 1.0.

- sg_focus:

  Character. Subgroup selection focus. One of "hr", "hrMaxSG", "maxSG",
  "hrMinSG", "minSG". Default "hr".

- fs.splits:

  Integer. Number of splits for consistency evaluation (or maximum
  splits when `use_twostage = TRUE`). Default 1000.

- m1.threshold:

  Numeric. Maximum median survival threshold. Default Inf.

- pconsistency.threshold:

  Numeric. Minimum consistency proportion. Default 0.90.

- stop_threshold:

  Numeric in `[0, 1]` or `NULL`. Early stopping threshold for
  consistency evaluation. When a candidate subgroup's consistency
  probability (`Pcons`) meets or exceeds this threshold, evaluation
  stops early — remaining candidates are skipped. Set to `NULL` to
  disable early stopping and force evaluation of all candidates up to
  `max_subgroups_search`. Default `0.95`.

  **Note:** Values \> 1.0 are not permitted. To disable early stopping,
  use `stop_threshold = NULL`, not a value above 1.

  Automatically reset to `NULL` (with a warning) when `sg_focus` is
  `"hrMaxSG"` or `"hrMinSG"`, as these compound criteria require
  comparing HR *and* size across all candidates.

- showten_subgroups:

  Logical. Show top 10 subgroups. Default FALSE.

- d0.min:

  Integer. Minimum control arm events. Default 12.

- d1.min:

  Integer. Minimum treatment arm events. Default 12.

- max.minutes:

  Numeric. Maximum search time in minutes. Default 3.

- minp:

  Numeric. Minimum prevalence threshold. Default 0.025.

- details:

  Logical. Print progress details. Default FALSE.

- maxk:

  Integer. Maximum number of factors per subgroup. Default 2.

- by.risk:

  Integer. Risk table interval. Default 12.

- plot.sg:

  Logical. Plot subgroup survival curves. Default FALSE.

- plot.grf:

  Logical. Plot GRF results. Default FALSE.

- max_subgroups_search:

  Integer. Maximum subgroups to evaluate. Default 10.

- vi.grf.min:

  Numeric. Minimum GRF variable importance. Default -0.2.

- use_twostage:

  Logical. Use two-stage sequential consistency algorithm for improved
  performance. Default FALSE for backward compatibility. When TRUE,
  `fs.splits` becomes the maximum number of splits, and early stopping
  is enabled. See Details.

- twostage_args:

  List. Parameters for two-stage algorithm (only used when
  `use_twostage = TRUE`):

  n.splits.screen

  :   Integer. Splits for Stage 1 screening. Default 30.

  screen.threshold

  :   Numeric. Consistency threshold for Stage 1. Default is
      automatically calculated to provide ~2.5 SE margin.

  batch.size

  :   Integer. Splits per batch in Stage 2. Default 20.

  conf.level

  :   Numeric. Confidence level for early stopping. Default 0.95.

  min.valid.screen

  :   Integer. Minimum valid Stage 1 splits. Default 10.

## Value

A list of class "forestsearch" containing:

- grp.consistency:

  Consistency evaluation results including:

  - out_sg: Selected subgroup based on sg_focus

  - sg_focus: Focus criterion used

  - df_flag: Treatment recommendations

  - algorithm: "twostage" or "fixed"

  - n_candidates_evaluated: Number evaluated

  - n_passed: Number passing threshold

- find.grps:

  Subgroup search results

- confounders.candidate:

  Candidate confounders considered

- confounders.evaluated:

  Confounders after variable selection

- df.est:

  Analysis data with treatment recommendations

- df.predict:

  Prediction data with recommendations (if provided)

- df.test:

  Test data with recommendations (if provided)

- minutes_all:

  Total computation time

- grf_res:

  GRF results object

- sg_focus:

  Subgroup focus criterion used

- sg.harm:

  Selected subgroup definition

- grf_cuts:

  GRF cut points used

- prop_maxk:

  Proportion of max combinations searched

- max_sg_est:

  Maximum subgroup HR estimate

- grf_plot:

  GRF plot object (if plot.grf = TRUE)

- args_call_all:

  All arguments for reproducibility

## Details

**Algorithm Overview:**

1.  **Variable Selection**: GRF identifies variables with treatment
    effect heterogeneity; LASSO selects most predictive

2.  **Subgroup Discovery**: Exhaustive search over factor combinations
    up to `maxk`

3.  **Consistency Validation**: Split-sample validation ensures
    reproducibility

4.  **Selection**: Choose subgroup based on `sg_focus` criterion

**Two-Stage Consistency Algorithm:** When `use_twostage = TRUE`, the
consistency evaluation uses an optimized algorithm that can provide
3-10x speedup:

- **Stage 1**: Quick screening with `n.splits.screen` splits eliminates
  clearly non-viable candidates

- **Stage 2**: Sequential batched evaluation with early stopping for
  candidates passing Stage 1

The two-stage algorithm is recommended for:

- Exploratory analyses with many candidate subgroups

- Large `fs.splits` values (\>200)

- Iterative model development

For final regulatory submissions, `use_twostage = FALSE` may be
preferred for exact reproducibility.

## References

- FDA Guidance for Industry: Enrichment Strategies for Clinical Trials

- Athey & Imbens (2016). Recursive partitioning for heterogeneous causal
  effects. PNAS.

- Wager & Athey (2018). Estimation and inference of heterogeneous
  treatment effects using random forests. JASA.

## See also

[`subgroup.consistency`](https://larry-leon.github.io/forestsearch/reference/subgroup.consistency.md)
for consistency evaluation details
[`forestsearch_bootstrap_dofuture`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md)
for bootstrap inference
[`forestsearch_Kfold`](https://larry-leon.github.io/forestsearch/reference/forestsearch_Kfold.md)
for cross-validation

Package website: <https://larry-leon.github.io/forestsearch/>

Source code: <https://github.com/larry-leon/forestsearch>

## Examples

``` r
if (FALSE) { # \dontrun{
# Example 1: Standard analysis (backward compatible)
result <- forestsearch(
  df.analysis = trial_data,
  sg_focus = "hr",
  hr.threshold = 1.25,
  pconsistency.threshold = 0.90,
  fs.splits = 400,
  details = TRUE
)

# Example 2: Fast exploratory analysis with two-stage
result_fast <- forestsearch(
  df.analysis = trial_data,
  sg_focus = "maxSG",
  hr.threshold = 1.15,
  pconsistency.threshold = 0.85,
  fs.splits = 500,
  use_twostage = TRUE,
  details = TRUE
)

# Example 3: Two-stage with custom parameters
result_custom <- forestsearch(
  df.analysis = trial_data,
  sg_focus = "hr",
  hr.threshold = 1.3,
  pconsistency.threshold = 0.95,
  fs.splits = 600,
  use_twostage = TRUE,
  twostage_args = list(
    n.splits.screen = 50,
    batch.size = 25,
    conf.level = 0.99
  ),
  parallel_args = list(plan = "multisession", workers = 4),
  details = TRUE
)
} # }
```
