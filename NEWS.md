# forestsearch 0.2.0

## New Features

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

## Internal

* `fit_causal_forest()` and `fit_causal_forest_glm()` accept
  `tune_grf` and pass `tune.parameters` to `grf::causal_forest()`
  and `grf::causal_survival_forest()`.

* Bootstrap bias correction extended for GLM outcomes with
  log-scale / identity-scale handling.

* Cross-validation extended to handle GLM parameters and propensity
  score re-estimation.

## References

* Dandl, S., Haslinger, C., Hothorn, T., Seibold, H., Sverdrup, E.,
  Wager, S., & Zeileis, A. (2024). What makes forest-based
  heterogeneous treatment effect estimators work? *Annals of Applied
  Statistics*, 18(1), 506-528.

* Rehill, P. (2025). How do applied researchers use the causal forest?
  A methodological review. *International Statistical Review*.

# forestsearch 0.1.0

* Initial CRAN release.
