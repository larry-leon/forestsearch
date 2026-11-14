# Redraw and examine one trial of an MRCT simulation

Redraws trial `k` of a
[`mrct_region_sims`](https://larry-leon.github.io/forestsearch/reference/mrct_region_sims.md)
run exactly as the simulation drew it, runs
[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
on its non-region patients (the region's patients as the test set) with
the run's arguments and multiplier resampling, and analyses four patient
sets: all patients, the region, region(\\G^\*\\) (the region without the
true harmed subgroup, when the trial has `flag_harm`) and region(\\\hat
G\\) (the region without the declared subgroup \\\hat H\\).

## Usage

``` r
mrct_case_study(
  dgm,
  seed,
  k,
  n = NULL,
  region_var = "z_regA",
  fs_args = list(),
  sim_args = list(rand_ratio = 1, draw_treatment = TRUE),
  analysis_time = 60,
  cens_adjust = 0,
  mr_inference = TRUE,
  mr_inference_args = list(ci_method = "field", include_complement = TRUE,
    field_complement = TRUE),
  entry_delay = NULL
)
```

## Arguments

- dgm:

  The DGM of the run
  ([`generate_aft_dgm_flex`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md),
  [`inject_mrct_structure`](https://larry-leon.github.io/forestsearch/reference/inject_mrct_structure.md)).

- seed:

  The run's base seed (`mrct_region_sims(seed = )`).

- k:

  The trial's index (`sim`).

- n:

  The trial size (`mrct_region_sims(n_sample = )`); `NULL` (default)
  uses `nrow(dgm$df_super)`, as
  [`mrct_region_sims()`](https://larry-leon.github.io/forestsearch/reference/mrct_region_sims.md)
  does.

- region_var:

  The region indicator column. Default `"z_regA"`.

- fs_args:

  Named list of
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  arguments. It may also hold the arguments that
  [`mrct_region_sims()`](https://larry-leon.github.io/forestsearch/reference/mrct_region_sims.md)
  takes directly (`sg_focus`, `maxk`, `hr.threshold`, `hr.consistency`,
  `pconsistency.threshold`, `confounders.name`, `conf_force`); those not
  given take
  [`mrct_region_sims()`](https://larry-leon.github.io/forestsearch/reference/mrct_region_sims.md)'s
  defaults, so a run's direct arguments and its `fs_args` together
  reproduce its identifier. The default `confounders.name` is that of
  [`mrct_region_sims()`](https://larry-leon.github.io/forestsearch/reference/mrct_region_sims.md):
  for a mechanism built with a spline, it leaves out the three internal
  columns of the spline variable `<var>` by their exact names
  `<var>_treat`, `<var>_k` and `<var>_k_treat`; a covariate such as
  `z_prior_treat` and the region indicator stay. The data, the outcome,
  treatment, event, id and potential-outcome names, the display flags
  and the (sequential) parallel plan are set here, as in
  [`mrct_region_sims()`](https://larry-leon.github.io/forestsearch/reference/mrct_region_sims.md).

- sim_args:

  Named list of
  [`simulate_from_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md)
  arguments, as in
  [`mrct_region_sims()`](https://larry-leon.github.io/forestsearch/reference/mrct_region_sims.md).
  Default `list(rand_ratio = 1, draw_treatment = TRUE)`.

- analysis_time:

  Calendar time of the analysis. Default 60.

- cens_adjust:

  Censoring adjustment on the log scale. Default 0.

- mr_inference:

  Logical. Run multiplier resampling on \\\hat H\\
  ([`fs_mr_inference`](https://larry-leon.github.io/forestsearch/reference/fs_mr_inference.md)).
  Default `TRUE`.

- mr_inference_args:

  List passed to `forestsearch(mr_inference_args = )`; `seed` defaults
  to `seed + k`. Default: the field method with the complement de-biased
  and given its own field.

- entry_delay:

  `NULL` (default) or `list(var = <column>, delay = a)`, passed to
  [`simulate_from_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md)
  (overrides `sim_args$entry_delay`).

## Value

An object of class `mrct_case_study`, a list with

- `trial`:

  The redrawn trial (data.frame).

- `fs`:

  The
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  result (its `mr_inference` element holds the MR result);
  [`mr_estimates_table`](https://larry-leon.github.io/forestsearch/reference/mr_estimates_table.md)`(x$fs)`
  gives the inference table.

- `sets`:

  The four-set table: `set`, `label`, `n`, `n_exp`, `n_ctl`,
  `events_exp`, `events_ctl`, and the unadjusted two-sample Cox fit
  `Surv(y_sim, event_sim) ~ treat_sim` on the set: `hr`, `lo`, `hi`
  (Wald 95\\ `itt_all`, `itt_region`, `region_Gstar` (when `flag_harm`
  exists) and `region_Ghat` (the whole region when nothing is declared).

- `in_Hhat`:

  Logical, per trial row: membership of \\\hat H\\ by
  [`get_dfpred`](https://larry-leon.github.io/forestsearch/reference/get_dfpred.md)
  on the selected definition.

- `confusion`:

  When `flag_harm` exists, a list of two integer vectors
  `c(TP, FP, FN, TN)` of \\\hat H\\ against \\H^\*\\ (`flag_harm == 1`),
  `nonregion` and `region`; else `NULL`.

- `Hhat`:

  The definition as
  [`mrct_region_sims()`](https://larry-leon.github.io/forestsearch/reference/mrct_region_sims.md)
  records it (`"none"` when nothing is declared).

- `mr`:

  `fs$mr_inference`.

- `trigger`:

  Logical: all-patient upper limit \< 1.0 and region hazard ratio \>
  1.0.

- `call_args`:

  `seed`, `k`, `n`, `region_var`, `analysis_time`, the
  [`simulate_from_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md)
  and
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  argument lists (without data).

- `warnings`:

  Warnings raised by
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md),
  muffled and kept here.

## Reproducing a simulated trial

[`mrct_region_sims()`](https://larry-leon.github.io/forestsearch/reference/mrct_region_sims.md)
draws trial `k` in a worker with `simulate_from_dgm(seed = seed + k)`,
under the L'Ecuyer-CMRG generator that the workers run;
`mrct_case_study()` switches to that generator for the draw (the same
[`simulate_from_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md)
arguments: the `sim_args` defaults overridden by `n`, `analysis_time`,
`cens_adjust` and `seed + k`, and `entry_var = "entrytime"` when
`dgm$df_super` has that column) and restores the caller's generator, and
the caller's random-number state, on exit. The identifier seeds its own
random steps (data preparation and the consistency splits), so it
returns the simulation's \\\hat H\\ under any generator. Multiplier
resampling runs after the draw under the caller's generator with
`seed = seed + k` unless `mr_inference_args` gives a seed. (Inside
[`mrct_region_sims()`](https://larry-leon.github.io/forestsearch/reference/mrct_region_sims.md)
the resampling, if requested through `fs_args`, runs under the workers'
generator, so its MR columns need not equal this object's.)

## See also

[`mrct_select_case`](https://larry-leon.github.io/forestsearch/reference/mrct_select_case.md),
[`mrct_region_sims`](https://larry-leon.github.io/forestsearch/reference/mrct_region_sims.md),
[`mr_estimates_table`](https://larry-leon.github.io/forestsearch/reference/mr_estimates_table.md).

## Examples

``` r
if (FALSE) { # \dontrun{
k <- mrct_select_case(res)$sim
cs <- mrct_case_study(dgm, seed = 613800000, k = k, n = 1000,
                      region_var = "z_region",
                      fs_args = list(sg_focus = "hrMaxSG", maxk = 2,
                                     conf_force = character(0)),
                      sim_args = list(rand_ratio = 1, draw_treatment = TRUE,
                                      max_entry = 11),
                      analysis_time = 40, cens_adjust = 50)
cs
plot(cs)
mr_estimates_table(cs$fs)
} # }
```
