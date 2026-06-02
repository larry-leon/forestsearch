# forestsearch 0.2.0

## New Features

### GRF frontier selection (experimental)

* New `grf_selection` argument (one of `"tree"` (default) or `"frontier"`) for
  the GRF subgroup routines `grf.subg.harm.survival()` and
  `grf.subg.harm.glm()`, and for `forestsearch(subgroup_method = "grf")`.
  `"tree"` is the standard policy-tree selection and is unchanged.

* `grf_selection = "frontier"` is an experimental alternative that consumes the
  same doubly-robust scores as the policy tree but selects differently: it
  enumerates single-covariate thresholds (and depth-2 covariate-pair
  conjunctions when `grf_depth >= 2`), scores each candidate subgroup by its
  doubly-robust harm-effect and size, takes the Pareto frontier, and selects one
  subgroup under `frontier_rule` (`"effMaxSG"` (default), `"eff"`, or
  `"maxSG"`), with `effect_neighborhood` giving the relative band for
  `"effMaxSG"`. The selected subgroup is always a single conjunction. Within
  `forestsearch()`, the frontier rule is taken from `sg_focus` (which the tree
  path ignores). Implemented for survival and GLM outcomes; membership flows
  through the existing estimation and bootstrap machinery unchanged.

* This mode is provided for comparison and exploration. In internal benchmarks
  the GRF policy tree matched or outperformed the frontier on harm recovery
  (the tree's global value objective aligns well with isolating harm), so the
  tree remains the recommended default.

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

### DINA Estimator (Difference in Natural Parameters)

* `dina_fit()` implements the DINA estimator of heterogeneous treatment
  effects on the natural-parameter scale (Gao & Hastie 2021) for
  Gaussian, binomial, Poisson, and Cox outcomes.  Uses an orthogonalized
  (Robinson / R-learner) estimating equation with cross-fitting and a
  sandwich variance.  Propensity and baseline nuisance functions are
  pluggable: built-in logistic / GLM / Cox fits, a numeric constant, or
  user-supplied closures (`propensity_method`, `baseline_method =
  list(eta0, eta1)` or `list(nu)`).  S3 methods `coef`, `vcov`,
  `predict`, `confint`, `summary`, and `print` follow the conventional
  contract with `(Intercept), X1..Xd` naming.

* `dina_subgroup()` searches the cross-fit DINA surface for the largest
  subgroup whose mean effect meets a threshold, returning the discovered
  signature and a BLP-analog effect (`a*'beta`).

* `dina_subgroup_bootstrap()` provides bootstrap inference for the
  discovered subgroup: a percentile CI on the BLP-analog effect (with
  `a*` fixed from the original subgroup) and, when `refit = TRUE`, a
  bootstrap CI for the within-subgroup standard model refit within the
  fixed signature on each resample.  Both CIs are conditional on the
  discovered signature; neither adjusts for signature selection.

* `dina_subgroup_refit()` fits the standard within-subgroup
  treatment-effect model (a plain Cox model for survival, a GLM
  otherwise) comparing the arms inside a discovered subgroup --- the
  conventional clinical estimate, distinct from the BLP-analog.
  Unadjusted by default, with optional confounder adjustment
  (`"none"` / `NULL`-automatic / explicit) and Cox stratification.

* **Reference validation.**  `dina_fit()` was validated against the
  authors' reference implementation (`github.com/ZijunGao/DINA`) on
  synthetic data across all four families: coefficients, sandwich
  standard errors, and predicted HTE surfaces agree to numerical
  precision under matched nuisance learners and cross-fitting folds.
  See the QC document `quarto/qc/dina_reference_validation.qmd`.



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

## References

* Gao, Z., & Hastie, T. (2021). Estimating heterogeneous treatment
  effects for general responses. *arXiv preprint* arXiv:2103.04277.

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
