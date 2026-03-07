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
