# Run the extreme-subgroups simulation study

Draws `n_sims` synthetic trials from a calibrated DGM via
[`simulate_from_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md),
fits `fit` in every pre-defined subgroup of every trial, and returns the
raw result matrices plus full provenance. With the default `fit`,
`benchmarks`, seeds, and skip rule, the matrices are bit-identical to
those produced by the `extreme_subgroups` vignettes for the matching
`baseline`, and the returned object is byte-compatible with their RDS
payload: `saveRDS(sims, path)` yields a file the design-comparison memo
loads unchanged.

## Usage

``` r
run_subgroup_sims(
  dgm,
  subgroups,
  n_sims,
  fit = subgroup_cox(survival::Surv(y_sim, event_sim) ~ treat_sim),
  baseline = c("resample", "fixed"),
  n = NULL,
  analysis_time,
  max_entry,
  cens_adjust,
  cutpoints = list(),
  benchmarks = benchmark_spec(),
  min_n = 5L,
  workers = NULL,
  seed_base = 0L,
  rand_seed_offset = 1000000L,
  hr_true = NULL,
  k_treat = NULL,
  future_packages = c("forestsearch", "survival"),
  validate = TRUE,
  verbose = FALSE
)
```

## Arguments

- dgm:

  A DGM from
  [`generate_aft_dgm_flex()`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md)
  (survival) or
  [`generate_glm_dgm()`](https://larry-leon.github.io/forestsearch/reference/generate_glm_dgm.md)
  (GLM). For `baseline = "fixed"` (survival only) it must carry
  `df_source`.

- subgroups:

  List of subgroup definitions, each a list with character scalars `id`
  (an expression over trial columns, exposed cut-points, `flag_itt`, and
  benchmark columns), `name`, and `grp`.

- n_sims:

  Number of simulated trials.

- fit:

  Per-subgroup analysis function `data -> c(estimate, upper)`; see
  [`subgroup_cox()`](https://larry-leon.github.io/forestsearch/reference/subgroup_cox.md).
  Default is an unstratified Cox fit on the simulated columns (the
  vignettes pass a grade-stratified formula); for a `"glm_dgm"` the
  default is
  [`subgroup_glm()`](https://larry-leon.github.io/forestsearch/reference/subgroup_glm.md)
  constructed from the DGM's `outcome_type` / `effect_measure` (and, for
  binary, its `adverse_outcome`).

- baseline:

  `"resample"` (random-X: draw `n` patients from the super-population
  per trial) or `"fixed"` (conditional-on-X: every trial is the
  `df_source` panel). Survival only; a `"glm_dgm"` always resamples.

- n:

  Patients per trial. Required for `"resample"` and for a `"glm_dgm"`;
  must be `NULL` for `"fixed"`, where N is `nrow(dgm$df_source)` by
  construction.

- analysis_time, max_entry, cens_adjust:

  Passed to
  [`simulate_from_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md)
  (administrative censoring time, recruitment window, calibrated
  censoring shift). Survival-only: supplying any of them with a
  `"glm_dgm"` is an error.

- cutpoints:

  Named list of scalars exposed as columns of each simulated trial so
  `id` expressions like `"age <= age_med"` resolve.

- benchmarks:

  A
  [`benchmark_spec()`](https://larry-leon.github.io/forestsearch/reference/benchmark_spec.md),
  or `NULL` to skip benchmark generation entirely (no RNG is consumed
  for it).

- min_n:

  Subgroups with fewer rows are skipped for that trial (row stays `NA`),
  default `5L` as in the vignettes.

- workers:

  Parallel workers. `NULL` (default) reuses a caller-established
  non-sequential
  [`future::plan()`](https://future.futureverse.org/reference/plan.html)
  when one exists – no worker-pool spawn or teardown, the mode the
  vignettes use after their `parallel-setup` chunk – and otherwise
  creates a multisession plan on `ceiling(0.90 * availableCores())`,
  restored on exit. A numeric value forces a plan for this call
  (restored on exit): a fraction in (0, 1) of
  [`availableCores()`](https://parallelly.futureverse.org/reference/availableCores.html),
  an integer count, or `1` for sequential. Results are identical under
  every setting (seeding is per-iteration); only wall time differs.
  `sim_config` records `workers`, `plan_reused`, and `t_plan_secs`.

- seed_base, rand_seed_offset:

  Seed scheme; defaults (`0L`, `1e6L`) reproduce the vignettes.

- hr_true, k_treat:

  Optional provenance values stored verbatim in the result (the
  vignettes store `0.70` and the calibrated `k_treat`). `hr_true` is the
  true effect on the fitter's estimate scale: a hazard ratio for Cox
  fits, a mean difference for
  [`subgroup_glm()`](https://larry-leon.github.io/forestsearch/reference/subgroup_glm.md)
  fits.

- future_packages:

  Packages loaded on workers, default `c("forestsearch", "survival")` as
  in the vignettes.

- validate:

  Run
  [`validate_subgroups()`](https://larry-leon.github.io/forestsearch/reference/validate_subgroups.md)
  first (default `TRUE`).

- verbose:

  Print progress lines.

## Value

An object of class `"subgroup_sims"`: a named list with `design`,
`n_sims`, matrices `sim_hrs` / `sim_ubs` / `sim_ns` (`n_sims` x
subgroups, columns named by subgroup `name`), `subgroups`, `sim_config`,
`cens_adjust`, `k_treat`, `hr_true`, `created`, `r_version`, and
`forestsearch_version`. For a `"glm_dgm"`: `design = "glm"`,
`cens_adjust = NULL`, and an additional `effect` field holding the
fitter's scale metadata (see
[`subgroup_glm()`](https://larry-leon.github.io/forestsearch/reference/subgroup_glm.md));
survival results are constructed exactly as before, with no new fields.

## Details

Per-iteration seeding is explicit (trial `ss` uses `seed_base + ss` for
the DGM draw and `seed_base + ss + rand_seed_offset` for benchmark
membership), so results are independent of the parallel schedule: any
`workers` setting, including sequential, gives identical matrices.

GLM dispatch: when `dgm` inherits `"glm_dgm"`
([`generate_glm_dgm()`](https://larry-leon.github.io/forestsearch/reference/generate_glm_dgm.md)),
trials are drawn by
[`simulate_from_glm_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_glm_dgm.md)
(`n` is required; the survival-only `analysis_time` / `max_entry` /
`cens_adjust` must not be supplied; `baseline = "fixed"` uses the stored
df_source panel), the default `fit` becomes
[`subgroup_glm()`](https://larry-leon.github.io/forestsearch/reference/subgroup_glm.md)
constructed from the DGM's own `outcome_type` and `effect_measure`
(binary DGMs also pass their `adverse_outcome`; a continuous DGM
resolves to the same fitter as before), and the result additionally
carries the fitter's `effect` scale metadata with `design = "glm"`. The
`sim_hrs` / `sim_ubs` matrices then hold that fitter's
`(estimate, upper bound)` pairs – the field names are retained for
structural compatibility, exactly as `hr.threshold` serves
generic-threshold duty for GLM outcomes in
[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md).

## See also

[`summary.subgroup_sims()`](https://larry-leon.github.io/forestsearch/reference/summary.subgroup_sims.md),
[`benchmark_spec()`](https://larry-leon.github.io/forestsearch/reference/benchmark_spec.md),
[`subgroup_cox()`](https://larry-leon.github.io/forestsearch/reference/subgroup_cox.md),
[`subgroup_glm()`](https://larry-leon.github.io/forestsearch/reference/subgroup_glm.md),
[`validate_subgroups()`](https://larry-leon.github.io/forestsearch/reference/validate_subgroups.md)

## Examples

``` r
if (FALSE) { # \dontrun{
sims <- run_subgroup_sims(
  dgm = dgm_uniform, subgroups = subgroups, n_sims = 1000,
  fit = subgroup_cox(survival::Surv(y_sim, event_sim) ~ treat_sim +
                       survival::strata(grade)),
  baseline = "fixed",
  analysis_time = 84, max_entry = 24,
  cens_adjust = cal_uniform$cens_adjust,
  cutpoints = list(age_med = median(gbsg$age)),
  hr_true = 0.70, k_treat = k_treat_uniform
)
saveRDS(sims, "results/extreme_sims_fixed_1000_payload.rds")
} # }
```
