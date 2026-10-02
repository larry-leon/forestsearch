# TASK (forestsearch): add the MRCT structure-injection module to R/

Date: 2026-10-02
Owner: Larry (approves) · Spec: chat · Execution: Claude Code (CC)
Repo: forestsearch, branch `feature/glm-extension` (the current local checkout). Do not switch branches.

## 0. First action

1. Copy `~/Downloads/TASK_fs-add-mrct-inject-module.md` to `dev/tasks/` and commit it.
2. Copy `~/Downloads/mrct_inject_dgm.R` to `R/mrct_inject_dgm.R`. Do not commit until G1 passes.
3. CC never pushes. Larry pushes; the project repo pins the resulting SHA.

Gates are stop-on-failure: on a failed gate, stop, record the failure in the REPORT (§4), end the task.

## 1. R/ change — called out separately

| Kind | Item | Classification |
|---|---|---|
| New file | `R/mrct_inject_dgm.R` — `expand_covariates`, `make_region_model`, `inject_mrct_structure`, `analysis_covariates`, `mrct_truth_by_region`, `plot_region_modifier`, `sim_region_metrics` | **New code.** Moves nothing, changes no existing function's behaviour, changes no method. |
| New file | `tests/testthat/test-mrct_inject_dgm.R` (§3) | New test only. |
| Edit | `DESCRIPTION`: bump the dev version (fourth component); `latentcor` → Suggests; `mvtnorm`, `Matrix` → Imports if not already there | Metadata only. |
| Regenerated | `NAMESPACE`, `man/` for the seven functions via `devtools::document()` | Mechanical. |

No existing R/ file is edited.

## 2. Edits to the module on import (behaviour unchanged)

- Delete `.fs_get()` and call the package-internal functions directly: `generate_aft_dgm_flex`, `prepare_working_dataset`, `define_subgroups`, `create_spline_variables`, `calculate_linear_predictors`, `prepare_censoring_model`, `calculate_hazard_ratios`. Verify each exists in the current R/ before editing (read the source; names, not descriptions).
- Add `@export` to the seven functions. Add `@importFrom mvtnorm rmvnorm`, `@importFrom Matrix nearPD`, and `@importFrom` lines for the `stats`, `graphics`, `grDevices`, `utils` functions used (`quantile`, `pnorm`, `plogis`, `uniroot`, `rbinom`, `ecdf`, `coef`, `vcov`, `rexp`, `complete.cases`, `setNames`, `modifyList`, `png`, `dev.off`, `adjustcolor`, `par`, `plot`, `hist`, `abline`, `legend`). Keep `latentcor` behind `requireNamespace()`.
- Leave all function signatures, defaults and numerics as supplied.

## 3. Smoke test (`tests/testthat/test-mrct_inject_dgm.R`) — runs on `survival::cgd`, which ships with R

Build on cgd, time to first serious infection (`enum == 1` rows of `survival::cgd`; `tte = tstop - tstart`, `event = status`, `treat = 1{treat == "rIFN-g"}`; covariates `age, height, weight` continuous, `female, autosom, steroids, propylac` binary), with:

```
inject_mrct_structure(seed_data = cgd1, outcome_var = "tte", event_var = "event",
  treatment_var = "treat", continuous_vars = c("age","height","weight"),
  factor_vars = c("female","autosom","steroids","propylac"), x_pred = "age",
  spline_spec = list(knot = 12, zeta = 25, log_hrs = log(c(0.35, 0.80, 1.30))),
  region = list(prevalence = 0.20, or_pred = 10, x1_vars = "female", or_x1 = 3, loghr = 0),
  n_super = 5000, expand = "copula", seed = 1)
```

Assertions: object inherits `aft_dgm_flex`; `nrow(df_super) == 5000`; `z_region` present with mean in [0.17, 0.23]; `mrct_truth_by_region()` returns three rows with finite AHR and `attr(,"smd_x_pred") > 0.3`; `simulate_from_dgm(dgm, n = 500, analysis_time = 500, max_entry = 120, seed = 1)` returns 500 rows with `event_sim` in {0,1}; `sim_region_metrics(sim, "z_region")$overall["hr"]` is finite; `analysis_covariates(dgm, "unobserved")` excludes `z_age` and `z_region`. Keep the test under ~10 s (n_super = 5000, no plotting).

## 4. Gate and report

**G1.** `devtools::document()` clean; `R CMD INSTALL .` succeeds; `library(forestsearch)` exports the seven functions and still exports `mrct_region_sims`, `simulate_from_dgm`, `generate_aft_dgm_flex`, `check_censoring_dgm`, `calibrate_cens_adjust`, and `fs_mr_inference` (confirm the last from the current NAMESPACE — the project's analysis task depends on it); `devtools::test()` passes the new test; `R CMD check --no-manual --no-vignettes` has 0 errors (warnings/notes recorded, not gating; pre-existing notes identified as such).

Commit on G1 pass. `REPORT_fs-add-mrct-inject-module.md` beside the other `REPORT_*` files: gate outcome, R and package versions, check output summary, test timing, and **the commit SHA on `feature/glm-extension`** — this SHA is what `mrct-consistency/deps.R` pins next.
