# Operating False-Discovery Rate of a Fitted ForestSearch Analysis

Characterizes the operating false-discovery rate (FDR) of a **fitted**
[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
analysis and tests whether requiring harm confirmation on the de-biased
effect reduces it. This is a **measurement** tool, not a calibration
tool: it does not solve for thresholds that control an FDR target. For
the \\(c_1, c_2)\\ pair the analysis actually used, it reports two
quantities under a null:

## Usage

``` r
fs_fdr_report(
  fs_analysis,
  df_analysis,
  null = c("hom", "harm"),
  c_confirm = c(0.9, 1, 1.25),
  nsims = 1000L,
  n = NULL,
  confounders = NULL,
  n_super = 25000L,
  analysis_time = 84,
  cens_adjust = log(1.5),
  event_min = 5L,
  n_workers = NULL,
  seed_base = 8316951L,
  quiet = FALSE
)
```

## Arguments

- fs_analysis:

  A fitted
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  object. Must carry `args_call_all` (the stored call arguments),
  including `hr.threshold` and `hr.consistency`.

- df_analysis:

  The analyst's data frame. Used only to default the simulated sample
  size (`n = nrow(df_analysis)`); the null data are generated from the
  GBSG-AFT DGM, not resampled from this frame.

- null:

  Character vector; which null(s) to simulate. `"hom"` is the
  homogeneous null (`k_inter = 0`: a real, uniform treatment effect, no
  modifier); `"harm"` is the null-harm-region null
  (`target_hr_harm = 1.0`: a real modifier tuned so the flagged region
  has HR = 1.0 while treatment stays active elsewhere). Default both.

- c_confirm:

  Numeric vector of de-biased-HR harm-confirmation thresholds at which
  to report B, i.e. the values `t_confirm` is swept over. Default
  `c(0.9, 1.0, 1.25)`. On real data harm confirmation is inert; here it
  annotates the null's operating characteristics only.

- nsims:

  Integer; null replicates per null. Default 1000.

- n:

  Integer or `NULL`; simulated trial size. Default `nrow(df_analysis)`.
  The published reference baseline used `n = 700`; pass that to
  reproduce it exactly.

- confounders:

  Character vector or `NULL`; candidate covariates. Default `NULL`
  inherits `confounders.name` from the fitted object (intersected per
  replicate with the DGM columns), falling back to the GBSG covariate
  set when the analysis carried none.

- n_super:

  Integer; super-population size for the DGM. Default 25000.

- analysis_time:

  Numeric; administrative censoring time (months). Default 84.

- cens_adjust:

  Numeric; censoring shift passed to
  [`simulate_from_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md).
  Default `log(1.5)`.

- event_min:

  Integer; minimum events for a replicate to be scored. Replicates below
  this (or with a single treatment arm) are counted as visible failures,
  never silently dropped. Default 5.

- n_workers:

  Integer or `NULL`; parallel workers for the outer replicate loop.
  Default `NULL` uses ~90% of physical cores. The inner
  [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  call always runs sequentially
  (`parallel_args = list(plan = "sequential")`) to avoid nested
  parallelism, which breaks forestsearch.

- seed_base:

  Integer; base seed. Each null uses a disjoint seed block and each
  replicate an explicit seed, so results are plan-invariant. Default
  8316951.

- quiet:

  Logical; suppress progress and per-null messages. Default `FALSE`.

## Value

An object of class `"fs_fdr_report"`: a list with

- `fdr`:

  Data frame, one row per null x `c_confirm`, with `A`, `A_lo`, `A_hi`
  (raw FDR and 95% Wilson CI), `B`, `B_lo`, `B_hi` (confirmed FDR and
  CI), `reduction` (`A - B`), and the accounting columns `n_valid`,
  `n_declared`, `n_deb_na`, `fs_errors`, `unhandled`.

- `hr_structure`:

  Data frame of the realized hazard ratios each DGM encodes – the
  structural check that the nulls are different worlds.

- `verdict`:

  `"CLEAN"` if no replicate failed across any null, else `"WARNING"`.

- `n_errors`:

  Total failed replicates (forestsearch errors + unhandled).

- `c1`, `c2`, `n`, `nsims`, `nulls`:

  The inherited thresholds and run settings.

- `call_args`:

  Reproducibility metadata.

## Details

- **A** (raw FDR): \\A = P(\text{declare any harm subgroup} \mid c_1,
  c_2)\\.

- **B(g)** (confirmed FDR): \\B(g) = P(\text{declare AND de-biased } HR
  \ge g)\\. Because harm confirmation can only remove declarations,
  \\B(g) \le A\\ always; the gap \\A - B(g)\\ is what the de-biasing
  step buys at confirmation threshold \\g\\.

All thresholds and identifier knobs are inherited from
`fs_analysis$args_call_all`, so the FDR reported is that of the analysis
as run. Multiplier resampling is forced on internally
(`mr_inference = TRUE`) regardless of the fitted setting, because B
requires the per-replicate de-biased HR.

### What "false discovery" means here

Under `"hom"` the treatment effect is real but uniform, so any declared
subgroup is a manufactured heterogeneity – the winner's-curse-of-search
rate. Under `"harm"` a real modifier exists and the flagged region sits
exactly at HR = 1.0, so a small upward fluctuation reads as harm; this
is the more adversarial null. The declared cut need not overlap
`flag_harm`: under either null any declaration is a false discovery.

### The confirmation accounting (why n_deb_na matters)

B(g) counts a declaration only when its de-biased HR is finite and
clears `g`; a declaration whose de-biased HR could not be computed
counts as unconfirmed. Such cases are counted in `n_deb_na` and
reported, so that `A - B` is never mistaken for pure de-biasing benefit
when part of it is de-biased-HR unavailability.

### Inheritance and rebinding

Identifier knobs (thresholds, `subgroup_method`, `consistency_method`,
`conf_force`, `conf.cont_jcuts`, `maxk`, MR arguments, ...) transfer
verbatim from `args_call_all`. The data-binding arguments – the five
column names and the confounder set – are rebound to the DGM, because
[`simulate_from_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md)
emits fixed names (`y_sim / event_sim / treat_sim / id / flag_harm`).
`mr_inference` is forced `TRUE` so B is computable.

## See also

[`forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
(and its vocabulary section),
[`fs_mr_inference`](https://larry-leon.github.io/forestsearch/reference/fs_mr_inference.md)
for the MR step this forces on,
[`mr_estimates_table`](https://larry-leon.github.io/forestsearch/reference/mr_estimates_table.md),
[`forestsearch_bootstrap_dofuture`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md)
for the full bootstrap,
[`fpr_calibration`](https://larry-leon.github.io/forestsearch/reference/fpr_calibration.md),
[`setup_gbsg_dgm`](https://larry-leon.github.io/forestsearch/reference/setup_gbsg_dgm.md),
[`calibrate_k_inter`](https://larry-leon.github.io/forestsearch/reference/calibrate_k_inter.md),
[`simulate_from_dgm`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md)

## Examples

``` r
if (FALSE) { # \dontrun{
# A fitted survival/Cox analysis (thresholds and knobs live on the object)
fit <- forestsearch(
  df.analysis = gbsg_df,
  outcome.name = "time", event.name = "status",
  treat.name = "treat", id.name = "id",
  confounders.name = c("er", "age", "meno", "pgr", "nodes", "size", "grade"),
  is.RCT = TRUE, est.scale = "hr",
  subgroup_method = "consistency", consistency_method = "resample",
  hr.threshold = 1.25, hr.consistency = 1.0,
  pconsistency.threshold = 0.90, maxk = 2,
  mr_inference = TRUE)

# Characterize its operating FDR and whether harm confirmation helps
rep <- fs_fdr_report(fit, gbsg_df, nsims = 1000, n_workers = 100)
rep                       # prints the realized-HR table and the A/B table
rep$fdr                   # one row per null x c_confirm

# Reproduce the published reference (n = 700): A ~ 0.138 (hom), 0.278 (harm)
fs_fdr_report(fit, gbsg_df, n = 700, nsims = 1000, n_workers = 100)
} # }
```
