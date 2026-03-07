# forestsearch NEWS

## forestsearch 0.2.0

### New functions

* `generate_aft_dgm_flex()` — general data-generating model (DGM) builder.
  Replaces the GBSG-specific `create_gbsg_dgm()` as the primary DGM
  constructor. Accepts any survival dataset and fits an accelerated failure
  time (AFT) super-population model with user-specified treatment effect
  heterogeneity parameters.

* `simulate_from_dgm()` — general simulator for drawing trial replicates from
  an `aft_dgm_flex` DGM. Replaces `simulate_from_gbsg_dgm()`. Column names
  in the returned data frame use underscore notation (`y_sim`, `event_sim`,
  `treat_sim`, `flag_harm`).

* `run_simulation_analysis()` (general version) — simulation wrapper that
  calls `simulate_from_dgm()` and accepts explicit column-name parameters,
  making it applicable to any DGM built with `generate_aft_dgm_flex()`. The
  GBSG dataset is now one application of this general pipeline rather than a
  separate code path.

* `setup_gbsg_dgm()` — convenience bridge function. Wraps
  `create_gbsg_dgm()` (see Deprecated below) and reshapes its output to the
  `aft_dgm_flex` class expected by `simulate_from_dgm()` and
  `run_simulation_analysis()`. Existing GBSG-based simulation scripts can
  adopt the general pipeline with a one-line change:
  `dgm <- setup_gbsg_dgm(model = "alt", k_inter = k, seed = seed)`.

### Deprecated functions

The following functions are retained and fully functional but will be removed
in a future version. Each emits a deprecation warning on first call with a
concrete migration example.

* `create_gbsg_dgm()` → use `generate_aft_dgm_flex()` or `setup_gbsg_dgm()`.

* `simulate_from_gbsg_dgm()` → use `simulate_from_dgm(analysis_time = Inf)`.
  Note: the new function defaults to `analysis_time = 48` (staggered-entry
  administrative censoring); pass `analysis_time = Inf` to match the legacy
  `max_follow = Inf` behaviour. Column names in the result also change — see
  the mapping table below.

  | Legacy column | General column |
  |---------------|----------------|
  | `y.sim`       | `y_sim`        |
  | `event.sim`   | `event_sim`    |
  | `treat`       | `treat_sim`    |
  | `flag.harm`   | `flag_harm`    |

### Deprecated parameters

* `run_simulation_analysis(max_follow)` → use `analysis_time`. If supplied,
  `max_follow` is silently forwarded to `analysis_time` with a warning.

* `run_simulation_analysis(muC_adj)` → use `cens_adjust`. If supplied,
  `muC_adj` is silently forwarded to `cens_adjust` with a warning.

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
  internally. This prevents deprecation-warning spam in functions that call
  these internally in loops or binary searches (`calibrate_k_inter()`,
  `get_dgm_with_output()`, `validate_k_inter_effect()`).

* `compute_dgm_cde()` now resolves the super-population data frame from
  `dgm$df_super_rand` (GBSG DGMs) or `dgm$df_super` (general `aft_dgm_flex`
  DGMs), making it compatible with both class hierarchies.

* `globals.R`: added `"sim_id"` to `utils::globalVariables()` to suppress a
  spurious `R CMD check` NOTE from `run_simulation_analysis.R`.

---

## forestsearch 0.1.0

* Initial release.
