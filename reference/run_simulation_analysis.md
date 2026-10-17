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
  methods = NULL,
  fs_params = list(),
  grf_params = list(),
  cox_formula = NULL,
  cox_formula_adj = NULL,
  n_sims_total = NULL,
  seed_base = 8316951L,
  verbose = FALSE,
  verbose_n = NULL,
  debug = FALSE,
  keep = character(0),
  keep_first_n = Inf
)
```

## Arguments

- sim_id:

  Integer. Simulation replicate index (used as seed offset).

- dgm:

  An `"aft_dgm_flex"` object from
  [`generate_aft_dgm_flex`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md)
  or
  [`setup_gbsg_dgm`](https://larry-leon.github.io/forestsearch/reference/setup_gbsg_dgm.md),
  or a `"glm_dgm"` object from
  [`generate_glm_dgm`](https://larry-leon.github.io/forestsearch/reference/generate_glm_dgm.md).

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

- methods:

  Optional specification of analysis arms that generalizes the `run_fs`
  / `run_fs_grf` / `run_grf` booleans. When non-`NULL` the booleans are
  ignored and **every** arm is run via
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  (through `.run_fs_analysis_gen()`), so all arms — including `dina-*`
  and `grf-*` — inherit the bias-corrected post-selection inference. Two
  forms are accepted:

  - A character vector of canonical `"<engine>-<focus>"` labels, e.g.
    `c("fs-effMaxSG", "dina-effMaxSG", "grf-tree")`, expanded internally
    to per-arm
    [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
    overrides. Engines are `"fs"` (`subgroup_method = "consistency"`),
    `"dina"`, and `"grf"` (`grf_selection = "frontier"`, or `"tree"` for
    the `grf-tree` label). Defaults applied during expansion:
    `selection_rule = "neighborhood"`; for `fs-*` arms `use_grf = TRUE`,
    `use_dina = TRUE`, `use_lasso = FALSE`, and `n.min = NULL` (the
    sample-size-adaptive rule).

  - A named list of explicit
    [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
    argument-override lists, label = name, e.g.
    `list("fs-effMaxSG" = list(subgroup_method = "consistency", sg_focus = "effMaxSG", use_grf = TRUE, use_dina = TRUE, use_lasso = FALSE), "grf-tree" = list(subgroup_method = "grf", grf_selection = "tree"))`.

  Each arm's overrides are merged on top of
  [`default_sim_params()`](https://larry-leon.github.io/forestsearch/reference/default_sim_params.md)
  overlaid with `fs_params`; per-arm overrides take precedence. The
  merge preserves an explicit `n.min = NULL`. Default `NULL` preserves
  the legacy boolean behavior.

- fs_params:

  Named list of ForestSearch parameter overrides. Any element of
  [`default_sim_params()`](https://larry-leon.github.io/forestsearch/reference/default_sim_params.md)
  can be overridden, including `parallel_args` (see *Parallel
  Processing* section below). Passing `NULL` or an empty list uses
  defaults throughout.

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

- keep:

  Character vector. Optional names of heavy diagnostic objects to attach
  as a list in `attr(result, "diagnostics")`. Allowed values (any
  combination):

  `"fs_full"`

  :   The full
      [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
      result object (minus the `df.est` analysis data frame, which is
      stripped to control memory). Includes `sg.harm`,
      `grp.consistency`, `confounders.candidate`,
      `confounders.evaluated`, GRF cuts, and timing.

  `"df.est"`

  :   The analysis data frame with `treat.recommend` column attached.
      Heavy; use sparingly.

  `"candidate_table"`

  :   The data table of candidate subgroups evaluated
      (`fs_full$grp.consistency$out_sg$result`). Medium-weight; useful
      for diagnosing which candidates were near-threshold.

  `"grf_full"`

  :   The full `grf.subg.harm.*()` result object.

  Default `character(0)` (attach nothing; current behavior).

- keep_first_n:

  Integer. When `keep` is non-empty, attach the requested diagnostics
  only if `sim_id <= keep_first_n`. Default `Inf` (attach for all
  replicates). Use this to capture full details for, e.g., the first 20
  simulations without bloating the result object for long production
  runs.

## Value

A `data.table` with one row per analysis method. Always-present scalar
columns cover: true-subgroup size/effect from the DGM, per-method
detection flag (`any.H`), estimated subgroup size (`size.H`, `size.Hc`),
effect estimates (`hr.H.hat`, `hr.Hc.hat`, `hr.itt`, and their true
counterparts), classification metrics against the DGM-stored subgroup
(`sens`, `spec`, `ppv`, `npv`), and pairwise concordance between methods
(`agree.*.*`, `kappa.*.*`). As of v0.2.0 the following ForestSearch
diagnostic columns are also always populated (`NA` for GRF rows):

- `sg.def`:

  Character. Cut-expression string of the identified subgroup, formed as
  `paste(fs_full$sg.harm, collapse = " & ")`. Empty string when no
  subgroup was identified. Read across all replicates to diagnose which
  cuts the method is selecting (e.g., `table(results$sg.def)`).

- `sg.n_factors`:

  Integer. Number of factors in the identified subgroup definition.

- `n_candidates_evaluated`:

  Integer. Number of candidate subgroups actually evaluated for
  consistency (before or at the early-stop candidate).

- `n_candidates_total`:

  Integer. Total number of candidate subgroups available.

- `n_passed`:

  Integer. Number of candidates meeting the consistency threshold.

- `consistency_algorithm`:

  Character. `"fixed"` or `"twostage"`.

- `early_stop_triggered`:

  Logical. Whether evaluation stopped early due to `stop_threshold`.

- `fs_minutes`:

  Numeric. Wall time for the FS run in minutes.

The analogous GRF row carries `sg.def` (from `grf_result$sg.harm.id`)
and `grf_selected_depth`. When `keep` is non-empty and
`sim_id <= keep_first_n`, the heavy objects are attached as
`attr(result, "diagnostics")` — a named list with one element per keeper
(`fs_full`, `grf_full`, etc.).

## GLM Parameters

GLM-specific parameters (`outcome_type`, `effect_measure`,
`offset.name`) must be passed inside `fs_params`, not as top-level
arguments. They route only to the estimation step
(`.extract_fs_estimates_gen`, `.extract_grf_estimates_gen`), not to
[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
itself, which uses Cox PH for subgroup identification in v0.1.x. Passing
these as top-level arguments will result in them being silently ignored.

## Parallel Processing

`run_simulation_analysis()` is designed to be called inside a
`foreach() %dofuture%` loop (one replicate per worker). In that idiom,
*outer* parallelism is provided by `%dofuture%` (one replicate per
worker) and the *inner*
[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
call should run sequentially within each worker. Running both layers
multisession produces nested parallelism: each outer worker tries to
spawn its own pool of inner workers, which `parallelly` rejects with a
300\\

To prevent this,
[`default_sim_params()`](https://larry-leon.github.io/forestsearch/reference/default_sim_params.md)
sets `parallel_args = list(plan = "sequential")` as the default for the
inner
[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
call. Users who call `run_simulation_analysis()` once at the top level
(e.g., for interactive debugging) and want the inner pipeline to run in
parallel can opt back in by passing
`fs_params = list(parallel_args = list(plan = "multisession", workers = N))`.

## See also

[`simulate_from_dgm`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md),
[`generate_aft_dgm_flex`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md),
[`setup_gbsg_dgm`](https://larry-leon.github.io/forestsearch/reference/setup_gbsg_dgm.md)

## Examples

``` r
if (FALSE) { # \dontrun{
dgm <- setup_gbsg_dgm(model = "null", verbose = FALSE)

# Inner forestsearch() runs sequentially by default; safe for both
# interactive use and embedding inside foreach() %dofuture% loops.
result <- run_simulation_analysis(sim_id = 1, dgm = dgm, n_sample = 500)

# Opt back into inner multisession for a one-off top-level call:
result <- run_simulation_analysis(
  sim_id    = 1,
  dgm       = dgm,
  n_sample  = 500,
  fs_params = list(parallel_args = list(plan = "multisession",
                                        workers = 4))
)
} # }
```
