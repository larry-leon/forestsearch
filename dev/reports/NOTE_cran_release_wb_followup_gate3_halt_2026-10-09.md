# NOTE — CRAN release Stage 3: halted at Gate 3 (2026-10-09)

Task: `dev/tasks/TASK_cran_release_wb_followup_2026-10-09.md`. The task's gates stop on failure. **Gate 3 failed** (`FAIL 8` in the full suite), so work stopped before Step 4. Nothing beyond the task doc is committed, and nothing was uploaded.

## State

- Stage 2 HEAD: `f43c4299b14144ccb79b3b218570be8c735b973d`.
- Task-doc commit: `50b1b9a5` (task doc only).
- Preconditions: P1 clean; P2 `Version: 0.4.0` and Stage 2 report present; P3 pandoc, quarto and pdflatex on `PATH`. All passed.
- Uncommitted working tree (Steps 1 and 3):
  - `R/mrct_inject_dgm.R`
  - `man/inject_mrct_structure.Rd`
  - 7 × `tests/testthat/test-*.R`

## Step 1 — HTML NOTE (Gate 1 PASSED)

- `man/inject_mrct_structure.Rd:69` before:
  `aft_dgm_flex with df_super = expanded population (+ z_\if{html}{\out{<region>}}), and`
- `R/mrct_inject_dgm.R:255` before:
  `#' @return aft_dgm_flex with df_super = expanded population (+ z_<region>), and`
- Cause: roxygen's markdown mode passed the raw `<region>` through as an HTML tag (`\out{}`).
- Fix: wrapped it in backticks, giving `` `z_<region>` ``. The Rd now has `\verb{z_<region>}`.
- Scan of all `man/*.Rd`: rendered every page with `tools::Rd2HTML` and checked for unknown tags. No other hits; the only `\out{<...>}` uses are standard `<div>`, and `<h5>` is valid HTML.
- Gate 1 results:
  - `Rd2HTML` output has 1 × `&lt;region&gt;` and 0 × raw `<region>`.
  - `NAMESPACE` is byte-identical.
  - The only `R/` diff line is roxygen (`#'`).

## Step 2 — Per-file test timing (`NOT_CRAN` unset, fresh R process each)

Total across 57 files: **151.41 s**.

| File | Elapsed (s) | Step 3 |
|---|---:|---|
| test-mr-inference.R | 85.867 | skipped |
| test-declaration-calibration.R | 13.618 | skipped |
| test-declaration-c0.R | 12.021 | skipped |
| test-return-shape-contract.R | 3.888 | skipped |
| test-binary-default-or.R | 3.101 | skipped |
| test-directive-c.R | 2.799 | skipped |
| test-grf-inclusion-band.R | 2.691 | skipped |
| test-get-fsdata-diagnostics.R | 2.54 | runs |
| test-cross-outcome-parity.R | 2.519 | runs |
| test-dfbeta-lm-influence-path.R | 2.345 | runs |
| test-threshold-pair-directive-a.R | 2.175 | runs |
| test-betaHhat-contract.R | 1.841 | runs |
| test-mr-closedform-fixture.R | 1.351 | runs |
| test-glm-pipeline-integration.R | 1.336 | runs |
| test-subgroup-sims.R | 1.289 | runs |
| test-fs-oc-null.R | 1.289 | runs |
| test-fs-oc-grid.R | 1.288 | runs |
| test-fs-oc-predict.R | 1.173 | runs |
| test-fs-family-report.R | 1.065 | runs |
| test-threshold-sync.R | 0.719 | runs |
| test-mr-field-complement.R | 0.69 | runs |
| test-input-validation.R | 0.618 | runs |
| test-mr-ij-residual-joint.R | 0.549 | runs |
| test-mr-inclusion-band.R | 0.539 | runs |
| test-mrct_inject_dgm.R | 0.427 | runs |
| test-knoise-population.R | 0.348 | runs |
| test-fs-oc-signed.R | 0.332 | runs |
| test-fs_dgm_scale.R | 0.267 | runs |
| test-mr-admission-rounded.R | 0.249 | runs |
| test-fs-dgm-feasibility.R | 0.233 | runs |
| test-declcal-rounding-alignment.R | 0.212 | runs |
| test-downstream-consumers.R | 0.176 | runs |
| test-mr-admission-set.R | 0.173 | runs |
| test-sg-focus-dispatch.R | 0.156 | runs |
| test-fs_mr_oc_summary.R | 0.151 | runs |
| test-sg-focus-transparency.R | 0.149 | runs |
| test-estimability-boundary.R | 0.138 | runs |
| test-dina-proposal-orientation.R | 0.128 | runs |
| test-cv-mixed-and-fatal.R | 0.121 | runs |
| test-bias-coverage-scale.R | 0.097 | runs |
| test-sg-focus-selection.R | 0.089 | runs |
| test-fs-focus-tag.R | 0.074 | runs |
| test-no-subgroup-cv.R | 0.064 | runs |
| test-binary-default-or-entry-points.R | 0.064 | runs |
| test-mr-residual-centering.R | 0.056 | runs |
| test-grf-membership-coding.R | 0.045 | runs |
| test-dina-bootstrap-sg-focus-whitelist.R | 0.042 | runs |
| test-fs_sym_root.R | 0.041 | runs |
| test-save-leff-calibration.R | 0.036 | runs |
| test-no-subgroup-dina-grf-summary.R | 0.032 | runs |
| test-cv-no-subgroup-edges.R | 0.032 | runs |
| test-check-thread-cap.R | 0.032 | runs |
| test-no-subgroup-unified-contract.R | 0.03 | runs |
| test-no-subgroup-dina-grf-plot.R | 0.03 | runs |
| test-cv-no-subgroup-dina-grf.R | 0.029 | runs |
| test-no-subgroup-bootstrap.R | 0.028 | runs |
| test-search-reproducibility.R | 0.018 | runs |

## Step 3 — CRAN skips

- Added a file-level `skip_on_cran()` as the first executable line of 7 files:
  - `test-mr-inference.R`
  - `test-declaration-calibration.R`
  - `test-declaration-c0.R`
  - `test-return-shape-contract.R`
  - `test-binary-default-or.R`
  - `test-directive-c.R`
  - `test-grf-inclusion-band.R`
- Each diff is exactly `+skip_on_cran()`; no test bodies were touched.
- `test-save-leff-calibration.R` and `test-check-thread-cap.R` are not skipped.
- Remaining total projected from Step 2: 151.41 − 123.97 = 27.44 s.
- Re-timed after the edits: **28.11 s**. This meets the ≤ 30 s target.

## Gate 3 — full suite (`NOT_CRAN=true`): FAILED

- `devtools::test()` gave `FAIL 8 | ERR 0 | SKIP 3 | PASS 5866`.
- The same failure counts appear at `50b1b9a5` (task-doc commit, no Stage 3 edits), checked in a temporary git worktree. **The failures were already there and are not caused by Stage 3.**

| File | Fails | Test | Nature |
|---|---:|---|---|
| `test-fs-oc-grid.R` | 4 | `fs_oc_predict() is identical to the 0.2.4 pre-refactor reference on both gates` | `identical()` against a saved `dev/` reference; values differ in the last few digits (e.g. `Esens` 0.91580692796093999 vs 0.91580692796094010) |
| `test-fs-oc-null.R` | 2 | `the alternative path is untouched by the null branch` | same kind of last-digit difference (e.g. `EnH` 86.308365156362925 vs …896) |
| `test-search-reproducibility.R` | 2 | `search + consistency results are invariant to worker count` | one of the two runs returned `sg.harm = NULL` and the other `c('!{x1 <= 0.6}', '{b1}')`; `M.2` differs |

Likely CRAN impact (not verified, because Step 4 was not run):

- The two `fs-oc` tests are guarded by `skip_if_not(file.exists(ref_path))`. `dev/` is not shipped in the build, so they should skip under `R CMD check`.
- `test-search-reproducibility.R` already calls `skip_on_cran()`.

## Not done

- Step 4: local `R CMD check --as-cran`.
- Step 5: commits and `check_win_devel` upload.
- Step 6: `cran-comments.md`.
- The Stage 3 report.

## Decision needed

1. **Proceed (recommended):** continue with Steps 4–6 and record Gate 3 as a deviation, citing the identical failure counts at `50b1b9a5`. This is defensible because none of these tests run in the CRAN check.
2. **Fix first:** fix the 8 failures before proceeding. That needs a separate authorization, since it means changing test assertions (e.g. switching `identical` to a tolerance comparison) or investigating the worker-count difference in `R/`.
