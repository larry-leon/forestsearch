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
# forestsearch NEWS

## forestsearch 0.2.0

### New functions

* `generate_aft_dgm_flex()` — general data-generating model (DGM) builder.
  Accepts any survival dataset and fits an accelerated failure time (AFT)
  super-population model with user-specified treatment effect heterogeneity
  parameters. This is the recommended starting point when building a DGM
  based on a dataset other than GBSG.

* `simulate_from_dgm()` — general simulator for drawing trial replicates from
  an `aft_dgm_flex` DGM. Supersedes `simulate_from_gbsg_dgm()` for new code.
  Column names in the returned data frame use underscore notation (`y_sim`,
  `event_sim`, `treat_sim`, `flag_harm`).

* `run_simulation_analysis()` (general version) — simulation wrapper that
  calls `simulate_from_dgm()` and accepts explicit column-name parameters,
  making it applicable to any DGM built with `generate_aft_dgm_flex()` or
  `setup_gbsg_dgm()`. The GBSG dataset is one application of this general
  pipeline rather than a separate code path.

* `setup_gbsg_dgm()` — the recommended entry point for all GBSG-based
  simulation work. Encodes the data preparation and subgroup definition from
  León et al. (2024) and returns an `aft_dgm_flex`-compatible object accepted
  by `simulate_from_dgm()` and `run_simulation_analysis()`. Existing scripts
  using `create_gbsg_dgm()` can migrate with a one-line change:
  `dgm <- setup_gbsg_dgm(model = "alt", k_inter = k, seed = seed)`.

### Superseded functions

`create_gbsg_dgm()` is superseded by `setup_gbsg_dgm()`. It remains fully
functional and continues to produce correct results; no existing GBSG
simulation scripts need to change. The distinction is that `setup_gbsg_dgm()`
returns an object of class `c("aft_dgm_flex", "gbsg_dgm")` compatible with
the general pipeline, whereas `create_gbsg_dgm()` returns only `"gbsg_dgm"`.
A `.Deprecated()` signal is emitted to encourage migration in new code.

`simulate_from_gbsg_dgm()` is superseded by `simulate_from_dgm()` for new
code. Column names in the output change from dot-notation to underscore
notation — see the mapping table below. Pass `analysis_time = Inf` to match
the legacy `max_follow = Inf` behaviour.

  | Legacy column | General column |
  |---------------|----------------|
  | `y.sim`       | `y_sim`        |
  | `event.sim`   | `event_sim`    |
  | `treat`       | `treat_sim`    |
  | `flag.harm`   | `flag_harm`    |

### Superseded parameters

* `run_simulation_analysis(max_follow)` → use `analysis_time`. If supplied,
  `max_follow` is forwarded to `analysis_time` with a warning.

* `run_simulation_analysis(muC_adj)` → use `cens_adjust`. If supplied,
  `muC_adj` is forwarded to `cens_adjust` with a warning.

### Bug fixes

The following bugs were discovered and fixed during the general pipeline
migration. All affected code paths were exercised by GBSG factor variables
(`v1`–`v7`) stored as `factor()` rather than `numeric()`.

* `lasso_selection()` (`get_FSdata_helpers.R`): `as.matrix()` on a data frame
  containing factor columns produced a character matrix that `cv.glmnet()`
  rejected. Factor columns with all-numeric levels are now coerced via
  `as.integer(as.character(.))` before matrix conversion.

* `process_conf_force_expr()` (`get_FSdata_helpers.R`): `mean()` applied to a
  factor column returned `NA`. Factor columns are now coerced to numeric before
  `mean()`, `median()`, and `quantile()` calls.

* `evaluate_comparison()` (`forestsearch_helpers.R`): the `<=` / `>=`
  operator applied to a factor column triggered an `Ops.factor` warning and
  returned `NA`. Factor columns are now coerced to numeric before comparison.

* `forestsearch()` (`forestsearch_main.R`): `df[, conf.screen]` dropped to a
  vector when `conf.screen` had length 1, causing `dummy()` to error on a
  non-data-frame input. Fixed by adding `drop = FALSE`.

* `default_grf_params_gen()` (`run_simulation_analysis.R`): `maxdepth` was
  initialised to `4`, exceeding the maximum of `3` accepted by
  `grf.subg.harm.survival()`. Corrected to `2` (matching the legacy default).

* `default_grf_params_gen()` (`run_simulation_analysis.R`): `sg.criterion`
  was set to `"hr"`, which is not a valid value. Corrected to `"mDiff"`
  (matching the legacy default).

### Internal changes

* `create_gbsg_dgm()` and `simulate_from_gbsg_dgm()` are now thin public
  wrappers that call `.create_gbsg_dgm_()` and `.simulate_from_gbsg_dgm_()`
  internally. This prevents warning spam in functions that call these
  in loops or binary searches (`calibrate_k_inter()`, `get_dgm_with_output()`,
  `validate_k_inter_effect()`).

* `compute_dgm_cde()` now resolves the super-population data frame from
  `dgm$df_super_rand` (GBSG DGMs) or `dgm$df_super` (general `aft_dgm_flex`
  DGMs), making it compatible with both class hierarchies.

* `globals.R`: added `"sim_id"` to `utils::globalVariables()` to suppress a
  spurious `R CMD check` NOTE from `run_simulation_analysis.R`.

---

## forestsearch 0.1.0

* Initial release.
