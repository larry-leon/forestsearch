# REPORT: add the MRCT structure-injection module to R/

Date: 2026-10-02
Task: `dev/tasks/TASK_fs-add-mrct-inject-module.md` (committed as received, `41aac1a7`)
Branch: `feature/glm-extension`. Not pushed.

## Gate outcome

**G1 (as amended): PASS.**

Amendment from Larry in chat: `R CMD check` is not run. G1 is `devtools::document()`, install,
`devtools::test(filter = "mrct_inject_dgm")`, and the export list.

| Check | Result |
|---|---|
| `devtools::document()` | Clean on the final run (see "Edits beyond the spec" for the four link warnings on the first run) |
| `R CMD INSTALL .` | `* DONE (forestsearch)` |
| Export list via `library(forestsearch)` | All 13 present: the seven new functions, plus `mrct_region_sims`, `simulate_from_dgm`, `generate_aft_dgm_flex`, `check_censoring_dgm`, `calibrate_cens_adjust`, `fs_mr_inference` |
| `devtools::test(filter = "mrct_inject_dgm")` | `[ FAIL 0 | WARN 3 | SKIP 0 | PASS 13 ]` |
| `R CMD check` | Not run (amendment) |

## Commit SHA to pin

**`2b4e77624155a220790831e77d713d9cbfe2a5cc`** (`2b4e7762`) on `feature/glm-extension`.

This is the SHA for `mrct-consistency/deps.R`. This report is committed separately, after that SHA.

## Versions

- R 4.6.1 (2026-06-24)
- forestsearch 0.3.5.9001 (bumped from 0.3.5.9000)
- roxygen2 8.1.0, devtools 2.5.2
- mvtnorm 1.4.2, Matrix 1.7.6
- latentcor: not installed on this machine, so the smoke test exercised the normal-score
  correlation fallback in `expand_covariates()`, not the `latentcor` path.

## Test timing

`tests/testthat/test-mrct_inject_dgm.R`: 0.63 s test time, 0.83 s wall (limit was ~10 s).

## Test warnings (not failures)

The test passes with three warnings, all the same message from `survival::survreg.fit`:
"Ran out of iterations and did not converge". All three come from the Weibull censoring-model
fit in `prepare_censoring_model()` (`R/generate_aft_dgm_helpers.R:312`) on the cgd seed data:

- twice inside `generate_aft_dgm_flex()` (`R/generate_aft_dgm_main.R:358` and `:393`);
- once from the module's own call (`R/mrct_inject_dgm.R:357`).

The test does not suppress them. The censoring model for the cgd-based DGM is therefore taken
from a non-converged fit. This was not investigated further; it is in existing code and the
supplied call, both outside this task's scope.

## What changed

| Kind | Item |
|---|---|
| New file | `R/mrct_inject_dgm.R` |
| New file | `tests/testthat/test-mrct_inject_dgm.R` |
| Edit | `DESCRIPTION`: version 0.3.5.9001; `Matrix`, `mvtnorm` added to Imports; `latentcor` added to Suggests |
| Regenerated | `NAMESPACE`; seven new `man/*.Rd` |
| Extra edit | `R/fs_mr_inference.R` (see next heading) |

Module edits per section 2:

- `.fs_get()` deleted; the seven internals are called directly. Each was confirmed by name in
  the current `R/` before editing (`generate_aft_dgm_flex`, `prepare_working_dataset`,
  `define_subgroups`, `create_spline_variables`, `calculate_linear_predictors`,
  `prepare_censoring_model`, `calculate_hazard_ratios`).
- `@export` added to the seven functions; `@importFrom` added for `mvtnorm::rmvnorm`,
  `Matrix::nearPD` and the listed `stats`/`utils`/`graphics`/`grDevices` functions.
- `latentcor` stays behind `requireNamespace()`.
- Signatures, defaults and numerics are as supplied.

## Extra edit: `fs_mr_inference` exported

Authorised by Larry in chat; not part of the original task, which said no existing R/ file is
edited.

- File: `R/fs_mr_inference.R`. One line added: `#' @export` above the function.
  `@keywords internal` is kept.
- No behaviour change. The API surface is widened by one function.
- Before this edit `fs_mr_inference` was defined but absent from `NAMESPACE`, so the G1 export
  check would have failed.

## Edits beyond the spec (documentation text only)

- The first `document()` run gave four "Could not resolve link to topic" warnings, because
  roxygen markdown read `log_hrs[1]`, `[2]`, `[3]` and `[0,1]` in `@param` text as links. They
  are now in backticks. No code changed.
- The header comment line "Requires the forestsearch R sources (or the installed package) to be
  loaded." was removed along with `.fs_get()`, since it no longer applies.

## Side issues, not fixed

- `@return` of `inject_mrct_structure` contains `z_<region>`; roxygen renders `<region>` as
  HTML-only, so it drops out of the text and PDF help.
- `analysis_covariates`, `mrct_truth_by_region`, `plot_region_modifier` and
  `sim_region_metrics` have no `@param`/`@return` tags. `R CMD check` will flag undocumented
  arguments when it is next run.
- `stats::qnorm`, `stats::cor`, `stats::sd`, `survival::coxph` and `survival::Surv` are used
  with `::` in the module but were not in the spec's `@importFrom` list, so none was added.
- `graphics` and `grDevices` are imported in `NAMESPACE` (already before this task) but are not
  listed in `DESCRIPTION` Imports.
