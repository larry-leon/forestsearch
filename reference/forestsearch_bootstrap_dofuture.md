# ForestSearch Bootstrap with doFuture Parallelization

Orchestrates bootstrap analysis for ForestSearch using doFuture
parallelization. Implements bias correction methods to adjust for
optimism in subgroup selection.

## Usage

``` r
forestsearch_bootstrap_dofuture(
  fs.est,
  nb_boots,
  seed = 8316951L,
  details = FALSE,
  show_three = FALSE,
  parallel_args = list(),
  digits = 4,
  mr_in_replicates = FALSE
)
```

## Arguments

- fs.est:

  List. ForestSearch results object from
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md).
  Must contain `df.est` (data frame) and `args_call_all` (list of
  arguments).

- nb_boots:

  Integer. Number of bootstrap samples (recommend 500-1000).

- seed:

  Integer. Random seed for reproducibility of bootstrap sample
  generation. Default `8316951L`. The value is passed to both
  [`bootstrap_ystar`](https://larry-leon.github.io/forestsearch/reference/bootstrap_ystar.md)
  (which constructs the \\n \times B\\ bootstrap index matrix) and
  [`bootstrap_results`](https://larry-leon.github.io/forestsearch/reference/bootstrap_results.md)
  (which re-runs the ForestSearch algorithm on each replicate); both
  calls must use the same seed to ensure bootstrap index alignment. Set
  to `NULL` for a non-reproducible run.

- details:

  Logical. If `TRUE`, prints detailed progress information. Default:
  `FALSE`.

- show_three:

  Logical. If `TRUE`, shows verbose output for first 3 bootstrap
  iterations for debugging. Default: `FALSE`.

- parallel_args:

  List. Parallelization configuration with elements:

  - `plan`: Character. One of "multisession", "multicore", "callr", or
    "sequential"

  - `workers`: Integer. Number of parallel workers

  - `show_message`: Logical. Show parallel setup messages

  If empty list, inherits settings from original forestsearch call.

- digits:

  Integer. Decimal places for construction-time numeric formatting in
  `FSsg_tab` (arm summaries, Diff, effect-estimate CI, bias-corrected
  CI). Default 4 to provide round-down headroom for downstream display
  reformatting (e.g.,
  [`summarize_bootstrap_results`](https://larry-leon.github.io/forestsearch/reference/summarize_bootstrap_results.md)`(digits = ...)`,
  which can round down but cannot recover precision lost at
  construction). Threaded through
  [`format_CI`](https://larry-leon.github.io/forestsearch/reference/format_CI.md)
  and
  [`SG_tab_estimates_glm`](https://larry-leon.github.io/forestsearch/reference/SG_tab_estimates_glm.md).

- mr_in_replicates:

  Logical. Whether multiplier resampling (MR) runs *inside* each
  bootstrap replicate. **Default `FALSE`**.

  `mr_inference` reaches each replicate through the reconstructed call
  in `fs.est$args_call_all`, so an analysis fitted with
  `mr_inference = TRUE` would otherwise run MR in every replicate – and
  nothing consumes the result. Under the default the flag is stripped
  from the replicate call, so that work is not done at all. This does
  *not* affect the top-level analysis: `fs.est` keeps whatever de-biased
  estimate it was fitted with.

  `TRUE` lets `mr_inference` propagate unchanged into each replicate
  *and* retains the result, returned as `$mr_replicates` – a list of
  length `nb_boots` whose `b`-th element is that replicate's
  `mr_inference` object, or `NULL` where the replicate failed or
  identified no subgroup. MR therefore runs in the replicates only if
  the original analysis set `mr_inference = TRUE`; running it and
  discarding the result is never a reachable combination.

  **The retained objects are not aggregable into an estimate for the
  original analysis.** Each replicate re-runs the search, selects its
  own subgroup, and de-biases *that* subgroup – so the collection
  describes MR's sampling behaviour across resampled selections, not a
  better estimate of the original quantity. Averaging them, or pooling
  their intervals, is an error. The bias-corrected estimate for the
  original analysis is the one this function already returns in
  `SG_CIs`, and the original subgroup's own MR estimate is on `fs.est`.

## Value

An `fs_bootstrap` object (a list with class `c("fs_bootstrap", "list")`)
containing:

- results:

  Data.table with bias-corrected estimates for each bootstrap iteration

- SG_CIs:

  List of confidence intervals for H and Hc (raw and bias-corrected)

- FSsg_tab:

  Formatted table of subgroup estimates

- Ystar_mat:

  Matrix (nb_boots x n) of bootstrap sample indicators

- H_estimates:

  Detailed estimates for subgroup H

- Hc_estimates:

  Detailed estimates for subgroup Hc

- mr_replicates:

  `NULL` under the default `mr_in_replicates = FALSE`. Otherwise a list
  of length `nb_boots`, element `b` holding replicate `b`'s
  `mr_inference` object (or `NULL` for a failed or no-subgroup
  replicate), each computed against *that replicate's* candidate family
  and selected subgroup. See the `mr_in_replicates` parameter for why
  these must not be averaged.

- summary:

  (If create_summary=TRUE) Enhanced summary with tables and diagnostics

- nb_boots:

  Integer. Number of bootstrap iterations requested. Used by
  [`print.fs_bootstrap`](https://larry-leon.github.io/forestsearch/reference/print.fs_bootstrap.md)
  to compute identification percentages without re-inspecting `results`.
  Added in v0.2.0.

- original_sg:

  Character vector. The subgroup identified by the primary
  [`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  call (`fs.est$sg.harm`), carried forward so that
  [`print.fs_bootstrap`](https://larry-leon.github.io/forestsearch/reference/print.fs_bootstrap.md)
  can report exact- and partial-match rates without retaining a
  reference to `fs.est`. Added in v0.2.0.

- outcome_type:

  Character. Outcome type of the underlying analysis: one of
  `"survival"`, `"binary"`, `"continuous"`, or `"count"`. Added in
  v0.2.0.

- effect_measure:

  Character. Effect measure label (`"HR"` for survival; `"OR"`, `"RR"`,
  `"RD"`, `"IRR"`, `"IRD"`, or `"MD"` for GLM outcomes). Added in
  v0.2.0.

- est.scale:

  Character. Estimation scale for confidence intervals (`"hr"` or
  `"1/hr"`). Added in v0.2.0.

## Bias Correction Methods

Two bias correction approaches are implemented:

1.  **Method 1 (Simple Optimism)**: \$\$H\_{adj1} = H\_{obs} -
    (H^\*\_{\*} - H^\*\_{obs})\$\$ where \\H^\*\_{\*}\\ is the new
    subgroup HR on bootstrap data and \\H^\*\_{obs}\\ is the new
    subgroup HR on original data.

2.  **Method 2 (Double Bootstrap)**: \$\$H\_{adj2} = 2 \times H\_{obs} -
    (H\_{\*} + H^\*\_{\*} - H^\*\_{obs})\$\$ where \\H\_{\*}\\ is the
    original subgroup HR on bootstrap data.

## Variable Naming Convention

- `H`: Original subgroup (harm/questionable, treat.recommend == 0)

- `Hc`: Complement subgroup (recommend, treat.recommend == 1)

- `_obs`: Estimate from original data

- `_star`: Estimate from bootstrap data

- `_biasadj_1`: Bias correction method 1

- `_biasadj_2`: Bias correction method 2

## Valid for every identifier configuration

Unlike multiplier resampling, the full bootstrap places no alignment
requirement on the identifier. It replays the entire pipeline on each
resample – refitting the front ends, re-enumerating the candidate family
and re-running selection – rather than linearizing the selection event,
so neither of MR's two conditions (selection ranking on the inferential
coefficient, and a fixed candidate family) bears on its validity. In the
equivalence result those conditions state, the bootstrap is the
reference standard, not a party to it.

It is therefore valid for all six identifier configurations, including
the three that rank on a native statistic: DINA under
`dina_select_statistic = "dina"`, GRF under
`grf_select_statistic = "dr"`, and GRF under `grf_selection = "tree"`.
For those, it correctly estimates the selection bias that ranking on the
native statistic induces on the reported effect \\\hat\beta\\. That is a
coherent quantity, and the estimate of it is valid – it simply sits
outside the framework the manuscript develops, which is about
identifiers whose selection map and reported effect are the same
functional.

The one place an alignment check does apply here is
`mr_in_replicates = TRUE`, which asks for MR *inside* each replicate;
see that parameter.

## Performance

Typical runtime: 1-5 seconds per bootstrap iteration. For 1000
bootstraps with 6 workers, expect 3-10 minutes total. Memory usage
scales with dataset size and number of workers.

## Requirements

- Original `fs.est` must have identified a valid subgroup

- Requires packages: `data.table`, `foreach`, `doFuture`, `survival`

- For plots: requires `ggplot2`

## See also

[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
for initial subgroup identification, and its vocabulary section for how
the full bootstrap (FB) implemented here relates to multiplier
resampling (MR)
[`fs_mr_inference`](https://larry-leon.github.io/forestsearch/reference/fs_mr_inference.md)
for MR, the fast refit-free approximation to this function
[`mr_estimates_table`](https://larry-leon.github.io/forestsearch/reference/mr_estimates_table.md)
to render MR estimates in this function's table layout
[`fs_fdr_report`](https://larry-leon.github.io/forestsearch/reference/fs_fdr_report.md)
for harm-confirmation threshold sweeps
[`bootstrap_results`](https://larry-leon.github.io/forestsearch/reference/bootstrap_results.md)
for the core bootstrap worker function
[`build_cox_formula`](https://larry-leon.github.io/forestsearch/reference/build_cox_formula.md)
for Cox formula construction
[`fit_cox_models`](https://larry-leon.github.io/forestsearch/reference/fit_cox_models.md)
for Cox model fitting

## Examples

``` r
if (FALSE) { # \dontrun{
# Run ForestSearch
fs_result <- forestsearch(
  df.analysis = mydata,
  outcome.name = "time",
  event.name = "status",
  treat.name = "treatment",
  confounders.name = c("age", "sex", "stage")
)

# Run bootstrap with bias correction
boot_results <- forestsearch_bootstrap_dofuture(
  fs.est = fs_result,
  nb_boots = 1000,
  parallel_args = list(
    plan = "multisession",
    workers = 6,
    show_message = TRUE
  ),
  create_summary = TRUE,
  create_plots = TRUE
)

# View results
print(boot_results$FSsg_tab)
print(boot_results$summary$table)

# Check success rate
mean(!is.na(boot_results$results$H_biasadj_2))
} # }
```
