# forestsearch 0.2.0

## New: `mean_r` and `mean_r_c` on the MR object

`fs_mr_inference()` now returns two diagnostics beside `selection_rate`:

- `mean_r` -- the mean of the infinitesimal-jackknife residual over the draws
  the IJ actually uses (`ok_H`), not over all `draws`.
- `mean_r_c` -- the same for the complement, `NA_real_` when no complement was
  fit.

**No estimate changes.** This is exposure of a quantity the correction already
formed; no computation was touched.

They exist to make one invariant checkable: the residual mean is **zero by
construction whenever both bias terms share a denominator**, which is the
convention the package implements -- `selection_bias` and `fixed_bias` both
average over the draws that produced a winner, and the IJ runs on that same
set. A non-zero value means the two terms are being normalised differently
somewhere. That was exactly the defect corrected earlier in this cycle, whose
signature was a non-zero residual mean; until now nothing asserted its absence
directly.

## Breaking: `max_subgroups_search` now defaults to `Inf`

`forestsearch()` capped the candidate pool at the top **200** subgroups before
consistency evaluation. The default is now `Inf` -- no truncation.

The cap has no counterpart in the method. The candidate family is every
conjunction of at most two conditions, less empty intersections and
size/event failures, and selection ranges over that family. Truncation applied
a *preview* ordering and discarded everything below the cut, so the subgroup
that is optimal under the actual selection criterion could be dropped before
it was ever evaluated -- most acutely for the size foci `"maxSG"` / `"minSG"`,
where small high-effect subgroups sort last, and for the band foci
`"effMaxSG"` / `"effMinSG"`.

**This changes FS results wherever the cap previously bound**, i.e. wherever
the post-screen pool exceeded 200 and the truncation warning fired. Where the
pool was smaller than the cap, nothing changes. Runs that relied on the cap for
runtime should set a finite value explicitly, knowing it can change which
subgroup is selected:

```r
forestsearch(..., max_subgroups_search = 200)   # previous behaviour
```

`sg_focus = "maxeff"` already forced `Inf` and is unaffected.

## Identifier alignment: new defaults (BREAKING) and an MR configuration guard

Multiplier resampling (MR) linearizes the selection event, which is valid only
when two conditions hold: selection ranks on the same within-subgroup effect
the analysis reports, and the candidate family is fixed rather than estimated
from the outcome data.  Several package defaults did not satisfy them.  Five
defaults changed and a hard configuration check was added.

### Breaking: `use_lasso` and `use_grf` now default to `FALSE`

`forestsearch()` previously ran a prognostic Cox-LASSO prefilter and a GRF
importance screen before enumerating candidates.  Both are now off by default.

**This changes the candidate family on every `forestsearch()` call with
`subgroup_method = "consistency"`, not only under MR.**  A different family
can yield a different selected subgroup, so results from 0.1.0 and earlier
will not necessarily reproduce at the new defaults.

Reproducing Leon et al. (2024, *Statistics in Medicine*, doi:10.1002/sim.10163),
which used the prognostic lasso prefilter, now requires setting the flags
explicitly:

```r
forestsearch(..., use_lasso = TRUE, use_grf = TRUE)
```

That configuration remains fully supported.  It is simply no longer the
default, because a model-based front end makes the family data-dependent.

### Changed defaults: three `match.arg()` re-orderings

The permitted values of all three are unchanged; only which one is the default
moved.  Each previously defaulted to its engine's native statistic, which is
not the effect the analysis reports.

| argument | was | now |
|---|---|---|
| `dina_select_statistic` | `"dina"` (subgroup-mean tau-hat) | `"effect"` |
| `grf_selection` | `"tree"` (policy-tree objective) | `"frontier"` |
| `grf_select_statistic` | `"dr"` (mean doubly-robust score) | `"effect"` |

`use_dina` is unchanged (`FALSE`).

### New: MR configurations are validated, not silently accepted

A misaligned configuration combined with MR previously ran to completion and
returned a correction that did not de-bias the reported effect.  It now raises
an error naming the offending argument and the setting that resolves it, at
three entry points:

* `forestsearch()` under `mr_inference = TRUE`;
* `forestsearch_bootstrap_dofuture()` under `mr_in_replicates = TRUE`;
* `forestsearch_Kfold()` under `mr_in_replicates = TRUE`.

Rejected: `dina_select_statistic = "dina"`, `grf_select_statistic = "dr"`,
`grf_selection = "tree"`, and -- for `subgroup_method = "consistency"` only --
any of `use_lasso` / `use_grf` / `use_dina` left `TRUE`.  That last check is
**not** applied to `subgroup_method` of `"dina"` or `"grf"`: their families
regenerate under resampling by construction, which is a property of the
methods, not a configuration error.

Errors only -- no warnings, no silent coercion.  All six configurations remain
reachable; the three misaligned ones now require explicit opt-in and are
unavailable under MR.

### New: `$family_status` on the returned object

`forestsearch()` now records the candidate-family status as a character
scalar, so a sweep can tabulate it without parsing text: `"no-front-end"`,
`"conditional-removable"` (a front end is on, and turning it off would reach
`"no-front-end"`), or `"conditional-inherent"` (DINA/GRF).  The first level is
deliberately not called `"fixed"`: it records that no fitted model shapes the
family on the observed data, which is weaker than the manuscript's fixed
family, that additionally requires resample-invariant cut locations.  Shown by
`print()` and `summary()`.  Descriptive only; it never raises a condition.

### Unaffected: the full bootstrap

`forestsearch_bootstrap_dofuture()` remains valid for all six configurations.
It replays the whole pipeline per resample rather than linearizing it, so the
alignment conditions do not bear on its validity -- it is the reference
standard those conditions are stated against.  See its documentation.

## Terminology: "gate" and "tier" are gone

The post-selection-inference surface was named after a "gate" while mostly
doing something else, and the two de-biasing routes were called "Tier 1" and
"Tier 2", which said nothing about what either one does.  Both vocabularies
are retired.  **None of the renamed names appeared in CRAN 0.1.0** -- the
whole family was added after that release -- so these are renames of
unreleased names, not breaking changes for any CRAN user, and no deprecation
forwards are provided.  `?forestsearch` now carries a *Vocabulary* section
defining every term.

* **"Tier 1" is now the full bootstrap (FB)** --
  `forestsearch_bootstrap_dofuture()`.  It resamples subjects and re-runs the
  whole search per replicate.
* **"Tier 2" is now multiplier resampling (MR)** -- `fs_mr_inference()`.  It
  resamples nothing; it perturbs the influence contributions of a single set
  of fits.  MR is post-selection inference on a completed analysis and cannot
  change which subgroup is identified.

| removed | replacement |
|---|---|
| `debias_gate` (`forestsearch()` argument) | `mr_inference` |
| `debias_gate_args` | `mr_inference_args` |
| `out$debias_gate` | `out$mr_inference` |
| `out$harm_flag_debiased` | `out$mr_harm_confirmed` |
| `gate_estimates_table()` | `mr_estimates_table()` |
| `fs_debias_gate()` (internal) | `fs_mr_inference()` |
| `mr_inference_args$gate` | `mr_inference_args$confirm_rule` |
| `mr_inference_args$t_gate` | `mr_inference_args$t_confirm` |
| `c_gate` (`fs_fdr_report()` argument and `rep$fdr` column) | `c_confirm` |
| `out$mr_inference$gate` | `out$mr_inference$settings` |
| `out$mr_inference$gate$type` | `out$mr_inference$settings$confirm_rule` |
| `.fs_dg_*`, `.fs_apply_debias_gate()` (internal) | `.fs_mr_*`, `.fs_apply_mr()` |

The permitted values of `confirm_rule` are unchanged (`"point"`, `"ci"`);
only the argument name moved.  Three files were renamed to match their main
function: `R/fs_mr_inference.R`, `R/fs_mr_inference_methods.R`,
`R/mr_estimates_table.R`.

Unrelated uses of the word were deliberately left alone: the calibrators'
hard tolerance gates (a genuine pass/fail that `stop()`s), and the
`tier1`/`tier2` risk-difference fallback ladder in
`glm_effect_estimators.R`.

**`mr_harm_confirmed` is three-valued and `NA` is not a negative result.**
`TRUE` = MR ran and harm is confirmed; `FALSE` = MR ran and it is not;
`NA` = MR did not run or could not be computed -- including the default
`mr_inference = FALSE`, and the case where no subgroup was identified.
Because `isTRUE(NA)` is `FALSE`, code that reads this field with `isTRUE()`
alone reports "harm not confirmed" for an analysis that never ran.  Branch on
`is.na()` first.

## Known limitations of the MR step

Documented here because they are measurable and were not previously written
down.  None is changed by the rename; all are behaviour, and all remain open.

* **MR is silently skipped for GLM outcomes under the default
  `consistency_method = "split"`.**  The consistency branch runs MR only when
  `consistency_method == "resample"` with a GLM effect, or on the survival/Cox
  path.  For a binary/count/continuous outcome at the default
  `consistency_method`, setting `mr_inference = TRUE` does nothing, emits no
  warning, and leaves `mr_harm_confirmed = NA`.  Pass
  `consistency_method = "resample"` to get MR on a GLM outcome.
* **MR's re-selection family does not fully replay the identifier's
  admissibility.**  The per-arm event minima (`d0.min` / `d1.min`) and
  `max_subgroups_search` are not re-applied, so MR re-selects over a superset
  of the candidates the search actually chose among.  This inflates the
  selection-bias term, i.e. over-corrects the de-biased estimate.
* **MR is recomputed inside every bootstrap replicate and CV fold, and the
  result is discarded.**  `mr_inference` propagates through
  `args_call_all`, and neither `forestsearch_bootstrap_dofuture()` nor
  `forestsearch_Kfold()` overrides it.  This is wasted work only -- the
  reported result comes from the original analysis -- but at the default
  `draws = 2000` it is not cheap.

`tests/testthat/test-mr-inference.R` is new and covers this surface, which
previously had no regression tests at all; the first limitation above is
pinned there as present behaviour.

## Behaviour Changes

**These change which subgroup a search returns.**  0.1.0 is on CRAN; an
analysis re-run under 0.2.0 can select a different subgroup than it did
before.  Anyone reproducing published results should pin the version or
re-record the selection.

* **`sg_focus = "hr"` (and its aliases `"eff"`, `"maxcons"`) no longer stops
  early.**  Early stopping halts the candidate scan at the first candidate
  reaching `stop_threshold` and then applies the focus's selection key to that
  prefix.  That is sound only when the key is fully determined by the
  enumeration order, and `Pcons` never is -- it is unknown until a candidate is
  evaluated.  Any key containing a `Pcons` term can therefore be beaten by a
  candidate the scan never reached.  This affects `"hr"` (`Pcons` primary) and
  also `"maxSG"` / `"minSG"`, where `Pcons` breaks ties in subgroup size and a
  truncated scan can halt before a tied candidate with higher `Pcons` is seen.
  `stop_threshold` is now reset to `NULL` for all of these, with a warning when
  it was supplied explicitly.

  A visible consequence: the selected subgroup no longer depends on
  `parallel_args$batch_size`.  Before this change, that performance knob could
  alter a scientific result.

* **`stop_threshold` now defaults to `pconsistency.threshold`** instead of a
  fixed `0.95`.  Previously the two defaults were `0.95` and `0.90`, so a
  candidate with `Pcons` in `[0.90, 0.95)` qualified but did not halt the scan.
  The default is a promise, so it tracks a user-supplied floor rather than
  reverting to the old gap.  This is an efficiency change, not a correctness
  fix: where early stopping is sound at all, the prefix sort recovers the same
  winner for any `stop_threshold >= pconsistency.threshold`.

## New Features

* **New `sg_focus = "maxeffCons"`**: the effect maximiser among candidates
  clearing `pconsistency.threshold`.  It sits between the two existing
  effect-oriented rules -- `"maxeff"` maximises effect with no consistency
  condition at all, while `"hr"`/`"maxcons"` maximises consistency and uses
  effect only as a tiebreak.  `"maxeffCons"` is the only focus for which early
  stopping is sound, since its key `(-hr, K)` is exactly the enumeration order.
  Multiplier resampling replays this rule correctly (it maps to MR's
  `"maxeff"`, which already ranges only over consistency-qualifying
  candidates).

* **New `sg_focus = "maxcons"`**, an alias for `"hr"`.  `"hr"` selects the
  consistency argmax with effect only as a tiebreak, which is the opposite of
  what the name suggests to most readers; `"maxcons"` names the rule it
  actually implements.  `"hr"` and `"eff"` are unchanged and keep their
  meaning.

* `sort_subgroups()` and `sort_subgroups_preview()` now accept the alias
  vocabulary (`"eff"`, `"effMaxSG"`, `"effMinSG"`, `"maxcons"`) when called
  directly, as the other entry points already did.

* `parallel_args$batch_size` is now documented in `forestsearch()`.  It was
  accepted and validated but appeared in no roxygen block.


### GLM Extension

* Subgroup identification for binary, continuous, and count outcomes
  using odds ratios (OR), risk differences (RD), relative risks (RR),
  incidence rate ratios (IRR), incidence rate differences (IRD), and
  mean differences (MD).

* New `outcome_type` and `effect_measure` parameters in `forestsearch()`
  control the GLM pipeline.  Identification uses outcome-type-aware
  LASSO screening and GLM-based subgroup evaluation; estimation uses
  pluggable estimator closures via `make_effect_estimator()`.

* `generate_glm_dgm()`, `calibrate_glm_interaction()`, and
  `simulate_from_glm_dgm()` provide a simulation suite for binary,
  continuous, and count DGMs with calibrated treatment-by-subgroup
  interactions.

* `run_simulation_analysis()` dispatches automatically on DGM class:
  `"glm_dgm"` objects route to `simulate_from_glm_dgm()` and
  `grf.subg.harm.glm()` without user intervention.

### Consistency Resampling Approximation

* New `consistency_method` parameter in `forestsearch()` and
  `subgroup.consistency()`.  `"split"` (default) runs the original
  repeated 50/50 split-and-refit consistency calculation and is
  byte-identical to prior behavior; `"resample"` replaces it with a
  multiplier (influence-function) approximation that returns the rate
  from a single subgroup fit.  Each candidate's half-split treatment
  effects are represented as `beta_hat +/- D`, where `D` is a
  multiplier sum of the per-subject treatment `dfbeta` contributions,
  giving the closed form `2 * Phi((beta_hat - c) / sigma_D) - 1` with
  `sigma_D` the robust (sandwich) SE of the treatment coefficient and
  `c` the consistency threshold on the comparison scale.

* The `"resample"` path applies to the Cox (survival) outcome and to
  GLM outcomes whose effect is a single model coefficient (OR, RR, RD,
  MD, IRR).  Configurations it cannot represent that way --- IRD,
  propensity-adjusted effects, or a non-convergent fit on a given
  candidate --- fall back to literal splitting automatically, so a
  consistency rate is always produced.  When `use_twostage = TRUE`,
  `"resample"` bypasses the two-stage split screening entirely.

* The method propagates through the bootstrap automatically (each
  resample iteration inherits `consistency_method`), where the saving
  compounds because the inner `forestsearch()` runs sequentially.

* New exported `consistency_resample()` (Cox + GLM rate from a single
  fit) and `consistency_resample_compare()` (validates the
  approximation against literal splitting on the survival path).

### GRF Standalone Subgroup Identification

* `grf.subg.harm.glm()` provides standalone GRF-based subgroup
  identification for binary, continuous, and count outcomes using
  `grf::causal_forest()` and `policytree::policy_tree()`.

* `adverse_outcome` parameter controls the Y-flip before
  `causal_forest()`.  When `TRUE` (the default for binary and count
  outcomes), the outcome is flipped (`1 - Y` for binary, `-Y` for
  continuous) so that the most-negative-CATE subgroup corresponds to
  the harm subgroup.  This aligns the GRF sign convention with
  `forestsearch()`'s internal GRF pre-screening and with
  `grf::causal_survival_forest()`'s RMST convention.

* `tune_grf` parameter added to `forestsearch()`,
  `grf.subg.harm.glm()`, and `grf.subg.harm.survival()`.
  When `TRUE`, enables cross-validated hyperparameter tuning
  (`tune.parameters = "all"`) in `causal_forest()`, tuning
  `min.node.size`, `mtry`, `sample.fraction`, `alpha`, and
  `imbalance.penalty`.  Default `FALSE` preserves existing behavior.
  Most impactful for observational data with propensity score
  adjustment; defaults are near-optimal for RCTs (Dandl et al., 2024).

* `dmin.grf` is on the **risk difference** scale for binary outcomes
  (e.g., `dmin.grf = 0.10` requires at least a 10 percentage point
  absolute difference in event rates).  Calibration documents are
  provided for both synthetic 4-confounder settings and the ACTG175
  12-confounder setting.

### Selection Criteria

* `sg_focus = "hrMaxSG"` and `"hrMinSG"` redefined to use
  effect-size neighborhood selection.  Among candidates whose
  effect is within `effect_neighborhood` (default 10%) of the
  maximum effect, `"hrMaxSG"` picks the largest sample size and
  `"hrMinSG"` picks the smallest.  Setting `effect_neighborhood = 0`
  reduces these to a strict max-effect filter.  Applies to all
  outcome types; for GLM ratio measures (OR, IRR), the neighborhood
  test is on the natural scale.

* New `effect_neighborhood` parameter (default `0.10`) in
  `forestsearch()` and `subgroup.consistency()`.

* New `selection_rule` parameter (default `"neighborhood"`) controls
  the candidate-inclusion logic for `"hrMaxSG"` / `"hrMinSG"`.  One of:
    - `"neighborhood"`: within `effect_neighborhood` of the maximum
      effect (the legacy v0.2.0 behavior).
    - `"pareto"`: Pareto-non-dominated set on (effect, N), both
      maximized.
    - `"both"`: intersection of `"neighborhood"` and `"pareto"`.
  Must be `"neighborhood"` for single-criterion focus values
  (`"hr"`, `"maxSG"`, `"minSG"`).

* GLM-natural `sg_focus` vocabulary added as aliases for the existing
  Cox-flavored names:
    - `"eff"`      is an alias for `"hr"`
    - `"effMaxSG"` is an alias for `"hrMaxSG"`
    - `"effMinSG"` is an alias for `"hrMinSG"`
  Both vocabularies produce identical results.  The `"eff*"` forms
  read more naturally in GLM contexts (continuous MD, binary OR/RR/RD,
  count IRR) where there is no hazard ratio.  Old code using the
  `"hr*"` forms continues to work without changes.

* New canonical threshold names `effect.threshold` and
  `consistency.threshold` in `forestsearch()`.  The legacy
  `hr.threshold` and `hr.consistency` continue to work; the new
  names take precedence when both are provided.  Like the
  `sg_focus` aliases, the new names read more naturally in GLM
  contexts.

* Frontier-preserving preview sort.  The candidate triage that
  determines which subgroups enter consistency evaluation now
  guarantees that the full Pareto frontier appears in the top-K
  candidates, regardless of `selection_rule`.  Previously,
  restrictive rules such as `selection_rule = "both"` could crowd
  low-N frontier members out of the top-K by filling slots with
  higher-N dominated candidates, producing different post-consistency
  diagnostics across rules even on the same data.  Affects candidate
  *filtering* only; the post-consistency winner selection in
  `sort_subgroups()` is unchanged, so the selected subgroup is
  unaffected.

### Dual-View Candidate Diagnostic

* `forestsearch()` and `subgroup.consistency()` accept
  `show_candidate_summary = TRUE` to print two diagnostic tables
  during a run:
    - **Pre-consistency preview**: all candidates entering
      consistency evaluation, with Frontier and InBand flags
      computed from the candidate HR/N values.
    - **Post-consistency summary**: passing candidates with Pcons,
      Frontier, InBand, and Selected flags.
  Together the two views make the rule's filter visible end-to-end.
  Column headers use abbreviated forms (`P` for Pcons, `OF` for
  on-frontier, `IB` for in-band, `S` for selected) with a legend
  printed below each table.

### Pareto-Frontier Diagnostic Suite

* `pareto_frontier_table()` renders the frontier as a formatted
  `gt` table or returns it as a `data.table`.  Works uniformly for
  survival (HR) and GLM (OR, RR, IRR, RD, IRD, MD) outcomes;
  effect-column label and scale handling derived from
  `fs$effect_measure`.  Selected subgroup is marked with a ★ and
  optionally highlighted.  The `include_dominated = TRUE` option
  extends the view to all passing candidates (not just the
  frontier), with `on_frontier` and `in_band` flag columns.

* `plot_pareto_frontier()` renders the frontier as a `ggplot`
  scatter with optional ε-band shading and split-derived CI bars.
  Shared style conventions across all frontier plots (theme_bw,
  steelblue step polyline, `#D55E00` highlight for the selected
  point).

* `plot_pareto_combined()` composes multiple `forestsearch` fits
  onto a single Pareto plot when their consistency-passing sets
  are identical.  Each selected subgroup is annotated with one
  or more `S1: <combo_label>` markers naming the configuration(s)
  that picked it; multiple combos picking the same subgroup yield
  stacked multi-line labels.  Returns `NULL` with a clear warning
  when the equality precondition fails (same row count, same
  subgroup-definition set, hr/N/E/K agreement within tolerance);
  the side-by-side panel composition is the natural alternative.

* `compute_frontier_cis()` computes three CIs per frontier member:
  a naive CI (model-based or robust SE), a half-jackknife split CI
  (Shao 1996), and an FSBC-mimic CI.  Pluggable into the table and
  plot functions via `ci_table = ...`.

* `explain_pareto_selection()` produces a verbal account of why the
  selected subgroup wins on the configured lexicographic criterion
  relative to other non-dominated candidates.  Format and verbosity
  configurable; supports markdown rendering via `results = "asis"`
  in Quarto chunks.

* `frontier_member_flags()` returns a per-subject membership matrix
  indicating which frontier subgroups each subject belongs to;
  useful for downstream within-frontier comparisons.

* `compare_selection_rules()` wrapper runs `forestsearch()` across
  multiple `(sg_focus, selection_rule)` tuple combinations with all
  other parameters held fixed.  Captures each run's stdout, builds
  per-combo Pareto plots, composes them via `patchwork`
  side-by-side, and auto-builds a `plot_pareto_combined()` view
  when the passing sets match.  Returns a structured
  `forestsearch_comparison` object with per-combo fs, ci_tab, plot,
  console, and diagnostics slots.

### Other

* Propensity score adjustment for observational studies via stabilized
  IPTW with GRF, LASSO, or logistic regression PS estimation.

* Forest plot (`plot_subgroup_results_forestplot()`) extended to
  support GLM effect measures with auto-detected display defaults,
  clinically meaningful axis limits, and extreme CI capping.

* `plot_sg_glm_outcomes()`: bar charts (binary) and
  point-and-errorbar (continuous) outcome visualization.

* `glm_effect_profile()`: delta-method treatment effect profiles
  across continuous biomarkers with natural cubic spline interactions.

## Breaking Changes (vs. v0.1.0)

* **`sg_focus = "hrMaxSG"` / `"hrMinSG"` selection rule changed.**
  Previously these used lexicographic sorts (size primary, effect
  secondary); now they use effect-size neighborhood selection (see
  *Selection Criteria* above).  Results from these focus values will
  differ from v0.1.0 runs on the same data.  Single-criterion focus
  values (`"hr"`, `"maxSG"`, `"minSG"`) are unchanged.

* **`showten_subgroups` argument removed.**  Renamed to
  `show_candidate_summary`.  The old name no longer works; callers
  passing `showten_subgroups = TRUE` will see an
  `unused argument` error.  The replacement triggers two
  diagnostic views (preview + summary) rather than the legacy fixed
  "top 10" display; see *Dual-View Candidate Diagnostic* above.

* **`pareto_frontier_table()` `digits` default lowered.**  Master
  `digits` argument now defaults to `2L` (was `3L`).  Per-column
  overrides (`digits_effect`, `digits_pcons`, `digits_ci`) inherit
  from the master when not explicitly set.  Rendered tables are more
  compact; pass `digits = 3L` to recover the prior precision.

## Bug Fixes

* **Fixed an incoherence between the two multiplier-resampling bias terms'
  denominators.** `selection_bias` averaged over the draws that produced a
  re-selected winner while `fixed_bias` averaged over all `B` draws. The
  infinitesimal-jackknife residual combines the two, so it mixed differently
  normalised quantities and its mean over the draws the IJ used was not zero --
  precisely the condition the finite-`B` correction needs in order to be the
  centered quantity it is meant to be. Both terms, and the IJ, are now taken
  over the same draws: those on which a selection occurred.

  This conditions the correction on identification, which is deliberate. The
  reported analysis exists only because a subgroup was identified, the selected
  subgroup's effect is already a conditional estimand, and a bootstrap replicate
  that identifies nothing contributes nothing to the full bootstrap that MR
  approximates.

  **Wherever `selection_rate < 1` the de-biased point estimate shifts by exactly
  the old residual mean** -- the incoherence was surfacing directly as a bias of
  that size. Intervals are essentially unchanged, and **`selection_rate = 1`
  analyses are unaffected to every digit**. `$mr_inference$selection_rate` is
  the field that tells you which case you are in. Measured on ACTG175 at
  `selection_rate = 0.865`:

  | quantity | before | after |
  |---|---:|---:|
  | `selection_bias` | 0.6623937 | 0.6623937 (identical) |
  | `fixed_bias` | 0.0053309 | 0.0908747 |
  | de-biased OR | 1.52988 | 1.40445 |
  | `tilde_V` | 0.2965604 | 0.2965604 (identical) |
  | `se_ij` | 0.4190578 | 0.4217818 |

  Note that `fixed_bias` is now a **conditional** expectation and is no longer
  mean-zero: the draws that clear the admission floor are exactly those where
  the perturbation pushed effects up, so the conditioning event correlates with
  the quantity being averaged. It shrinks monotonically to zero as
  `selection_rate` approaches 1 -- the behaviour of a selection effect rather
  than an artefact. Both the shipped DINA configuration and GRF on ACTG175 ran
  at `selection_rate = 1`.

* Fixed cross-validation reusing a cached `grf_res` / `dina_res` inside every
  training fold.  `forestsearch_Kfold()` and `forestsearch_tenfold()` rebuild
  each fold's call from `args_call_all`, which carries all of `forestsearch()`'s
  formals -- including a user-supplied cached fit.  The bootstrap has always
  nulled `grf_res`, `grf_cuts`, `dina_res` and `dina_cuts` for exactly this
  reason, and CV nulled `ps_hat` on the same reasoning while leaving the other
  four.  A fold therefore ranked over a candidate family estimated on **all**
  the data, including the fold being held out, so `sensCV`, `ppvCV` and
  exact-match were optimistic by construction.  All four are now nulled at both
  entry points.

  **This changes CV metrics for anyone passing a cached fit**, and the effect is
  not small.  Measured on ACTG175 (binary, OR, `Kfolds = 3`), cached versus
  uncached:

  | method | metric | cached, before | uncached | cached, after |
  |---|---|---:|---:|---:|
  | GRF | `ppv_H` | 0.6486 | 0.4896 | 0.4896 |
  | GRF | exact match | 0.3333 | 0.0000 | 0.0000 |
  | DINA | `sens_H` | 0.4688 | 0.2188 | 0.2188 |
  | DINA | `ppv_H` | 0.3000 | 0.1538 | 0.1538 |

  After the fix the cached and uncached runs agree exactly on every metric.
  Runs that passed `grf_res = NULL` / `dina_res = NULL` (the default) are
  unaffected -- there was nothing to reuse.

* Fixed the influence path erroring on every `lm()` fit, which disabled
  multiplier resampling (MR) for continuous outcomes entirely.  The internal
  one-step `dfbeta` read `fit$weights` unconditionally; an unweighted `lm()`
  has `NULL` there rather than a vector of ones, so the computation was
  non-conformable and failed.  `effect_measure = "MD"` is the only measure
  fitted with `lm()`, and the failure was silent at every level:
  `consistency_resample()` returned `NA` rates, the consistency engine fell
  back to literal splitting while `consistency_method = "resample"` was in
  force, and `fs_mr_inference()` could fit no candidate at all -- so
  `mr_inference = TRUE` on a continuous outcome produced no correction and
  reported nothing.  **Any MD analysis run with `mr_inference = TRUE` between
  the introduction of the one-step `dfbeta` and this release has no MR result
  and should be re-run.**  Binary, count and survival outcomes are unaffected:
  all fit with `glm()` or `coxph()`, where the weights are always present.

* Fixed `grf.subg.harm.glm()` ignoring `adverse_outcome`: the Y-flip
  before `causal_forest()` was accepted but never applied, causing the
  policy tree to identify the complement of the true subgroup in
  benefit-search scenarios.

* Fixed `run_simulation_analysis()` hardcoding `simulate_from_dgm()`:
  added class-based dispatch so `"glm_dgm"` objects route to
  `simulate_from_glm_dgm()`.

* Fixed `run_simulation_analysis()` GRF dispatch: binary/continuous
  outcomes now route to `grf.subg.harm.glm()` instead of
  `grf.subg.harm.survival()`.  Also removes `event.name` from
  `grf_args` (not accepted by `grf.subg.harm.glm()`) and propagates
  `adverse_outcome` with the same default as `forestsearch()`.

* Fixed `valid_pnames` whitelist in `run_simulation_analysis()`:
  added `outcome_type`, `effect_measure`, `offset.name`, and
  `tune_grf` so GLM parameters pass through to `forestsearch()`.

* Removed stale duplicate `grf.subg.harm.glm()` from `grf_main.R`
  that conflicted with the canonical implementation in
  `grf_subg_harm_glm.R`.

* Fixed `merge(by = "id")` ignoring `id.name` parameter in
  `forestsearch()` output.

* Fixed `geom_errorbarh(height = 0, ...)` deprecation warning under
  `ggplot2 >= 3.5` in `plot_pareto_frontier()` and
  `plot_pareto_combined()`.  Replaced with `geom_linerange()` which
  produces an identical visual without the warning.

* Fixed non-ASCII characters in `plot_pareto_combined.R` for R CMD
  check compliance: em-dashes in user-visible string literals
  replaced with ASCII hyphens; ε and ≥ glyphs replaced with R
  unicode escapes (`\u03b5`, `\u2265`).  Plot output is unchanged.

* Fixed bare `<<-` assignment leakage inside nested
  `capture.output()` / `tryCatch()` in `compare_selection_rules()`.
  Replaced with an explicit holder environment; the prior pattern
  could silently write to the global environment under certain
  R session configurations.

## Internal

* `fit_causal_forest()` and `fit_causal_forest_glm()` accept
  `tune_grf` and pass `tune.parameters` to `grf::causal_forest()`
  and `grf::causal_survival_forest()`.

* Bootstrap bias correction extended for GLM outcomes with
  log-scale / identity-scale handling.

* Cross-validation extended to handle GLM parameters and propensity
  score re-estimation.

* New internal helper `.normalize_sg_focus()` in
  `forestsearch_helpers.R` translates the `"eff*"` vocabulary to
  the canonical `"hr*"` form at entry to `forestsearch()` and
  `subgroup.consistency()`.  Downstream code is keyed on the
  canonical form; the aliases are purely user-facing.

* New internal helper `extract_candidate_diagnostics()` in
  `compare_selection_rules.R` slices the `CANDIDATE EVALUATION
  PREVIEW` and `CANDIDATE EVALUATION SUMMARY` blocks out of
  captured forestsearch stdout for focused presentation.

* Pareto frontier on (effect, N) -- both maximized -- attached to
  consistency results as `out_sg$pareto_frontier`.  Post-hoc
  diagnostic listing non-dominated alternatives to the selected
  subgroup; not used for selection.

* `Pcons` excluded from the value-tolerance check in
  `plot_pareto_combined()`'s equality precondition.  Pcons can
  legitimately drift across rules because the preview-sort change
  alters queue order and therefore the random-split state each
  candidate consumes; equality is keyed on subgroup definitions
  and hr/N/E/K only.

* GLM consistency resampling is wired through an explicit
  `glm_resample_spec` threaded from `forestsearch()` to the
  consistency evaluators (rather than introspecting the estimator
  closure).  `consistency_resample()` gained a `comparison_threshold`
  argument so the engine passes the already-resolved comparison-scale
  threshold without a second log-transform.  New internal helpers
  `.glm_resample_supported()` and `.consistency_via_splits()` gate the
  supported measures and provide the shared split fallback.

## References

* Dandl, S., Haslinger, C., Hothorn, T., Seibold, H., Sverdrup, E.,
  Wager, S., & Zeileis, A. (2024). What makes forest-based
  heterogeneous treatment effect estimators work? *Annals of Applied
  Statistics*, 18(1), 506-528.

* Rehill, P. (2025). How do applied researchers use the causal forest?
  A methodological review. *International Statistical Review*.

* Shao, J. (1996). Bootstrap model selection. *Journal of the
  American Statistical Association*, 91(434), 655-665.

# forestsearch 0.1.0

* Initial CRAN release.
