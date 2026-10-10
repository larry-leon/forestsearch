# REPORT — CRAN release, Stage 3: win-builder follow-up (forestsearch 0.4.0)

- Date: 2026-10-09
- Task: `dev/tasks/TASK_cran_release_wb_followup_2026-10-09.md`
- Branch: `feature/glm-extension`
- Outcome: all steps done. Gate 3 was amended (see §1). No push, no tag, no CRAN submission.

## 0. Preconditions and commits

- P1: `git status --porcelain` was empty.
- P2: `Version: 0.4.0`; Stage 2 report present. Stage 2 HEAD: `f43c4299b14144ccb79b3b218570be8c735b973d`.
- P3: all three tools found on `PATH`:
  - pandoc: `/Applications/RStudio.app/.../quarto/bin/tools/aarch64/pandoc`
  - quarto: `/Applications/RStudio.app/.../quarto/bin/quarto`
  - pdflatex: `~/Library/TinyTeX/bin/universal-darwin/pdflatex`
- Commits:
  - `50b1b9a5`: task doc
  - `d82cb0ce`: Step 1 docs fix
  - `315bd2ba`: Step 3 test skips
  - `8ac745af`: check log
  - `64283b0d`: cran-comments
  - this report, committed together with the Gate 3 halt note

## 1. Gate 3 amendment (Larry's disposition)

- **Amendment:** Gate 3 now reads "FAIL 0 except the 8 pre-existing failures enumerated in the halt note (identical at `50b1b9a5`)". Disposition option 1: continue with Steps 4–6.
- Halt note: `dev/reports/NOTE_cran_release_wb_followup_gate3_halt_2026-10-09.md`.
- The full suite was run once (`NOT_CRAN=true`): `FAIL 8 | ERR 0 | SKIP 3 | PASS 5866`.
- The same failure counts appear at `50b1b9a5` (a clean worktree without the Stage 3 edits).
- No task was opened for these failures; they are recorded below as findings.

### Findings — 8 pre-existing test failures

- **F1.** `test-fs-oc-grid.R`, 4 failures. Test: `fs_oc_predict() is identical to the 0.2.4 pre-refactor reference on both gates`.
  - It uses `identical()` against a saved `dev/` reference, and the values differ in the last few digits.
  - Example: `Esens` 0.91580692796093999 vs 0.91580692796094010.
- **F2.** `test-fs-oc-null.R`, 2 failures. Test: `the alternative path is untouched by the null branch`.
  - Same kind of last-digit difference.
  - Example: `EnH` 86.308365156362925 vs 86.308365156362896.
- **F3.** `test-search-reproducibility.R`, 2 failures. Test: `search + consistency results are invariant to worker count`.
  - One run returned `sg.harm = NULL`; the other returned `c('!{x1 <= 0.6}', '{b1}')`. `M.2` also differs.
- **CRAN impact:** none of the 8 runs in the CRAN check.
  - F1 and F2 skip when the `dev/` reference is absent, and `dev/` is not shipped in the build.
  - F3 already calls `skip_on_cran()`.
  - The Step 4 check confirms this: tests passed.

## 2. Step 1 — HTML NOTE (Gate 1 passed)

- `man/inject_mrct_structure.Rd:69` before:
  `aft_dgm_flex with df_super = expanded population (+ z_\if{html}{\out{<region>}}), and`
- `R/mrct_inject_dgm.R:255` before:
  `#' @return aft_dgm_flex with df_super = expanded population (+ z_<region>), and`
- Cause: roxygen's markdown mode passed `<region>` through as raw HTML (`\out{}`).
- Fix: wrapped it in backticks (`` `z_<region>` ``). The Rd now has `\verb{z_<region>}`.
- Scan of all `man/*.Rd` (every page rendered with `tools::Rd2HTML`): no other non-HTML tags. The only `\out{<...>}` uses are `<div>`.
- Gate 1 results:
  - `Rd2HTML` output has 1 × `&lt;region&gt;` and 0 × raw `<region>`.
  - `NAMESPACE` is byte-identical.
  - The only `R/` diff line is roxygen.

## 3. Step 2 — test timing (`NOT_CRAN` unset, fresh R process per file)

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

## 4. Step 3 — CRAN skips

- Added a file-level `skip_on_cran()` as the first executable line of 7 files:
  - `test-mr-inference.R`
  - `test-declaration-calibration.R`
  - `test-declaration-c0.R`
  - `test-return-shape-contract.R`
  - `test-binary-default-or.R`
  - `test-directive-c.R`
  - `test-grf-inclusion-band.R`
- `test-save-leff-calibration.R` and `test-check-thread-cap.R` are not skipped.
- Remaining total projected from Step 2: 27.44 s. Re-timed after the edits: **28.11 s**, which meets the ≤ 30 s target.
- Full suite: see §1 (amended gate passed).

## 5. Step 4 — local release check (Gate 4 passed)

- Temp library under the scratchpad:
  - `weightedsurv` 0.1.0 (CRAN) and `latentcor` 2.0.2 (CRAN), both built from source.
  - 26 missing dependencies installed as CRAN macOS binaries, same as the Stage 1 deviation (no gfortran on this machine).
  - Nothing went into the user library.
- Build: `R CMD build` with vignettes. Exit 0, 54 s; tarball 3,349,100 bytes.
- Check environment:
  ```
  R_LIBS=<scratch>/rlib
  _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE
  _R_CHECK_LIMIT_CORES_=TRUE
  OMP_THREAD_LIMIT=2
  _R_CHECK_EXAMPLE_TIMING_CPU_TO_ELAPSED_THRESHOLD_=2.5
  _R_CHECK_TEST_TIMING_CPU_TO_ELAPSED_THRESHOLD_=2.5
  ```
- Check command: `R CMD check --as-cran`. Wall time 199 s (Stage 1: 331 s).
- Platform: R 4.5.2, aarch64-apple-darwin20.
- **Status line: `Status: 1 NOTE`** (0 errors, 0 warnings).
- **Tests:** `checking tests ... [19s/17s] OK`; `testthat.Rout` reports elapsed 17.127 s. In Stage 1 the tests took 133 s.
- The NOTE, verbatim:
  ```
  * checking HTML version of manual ... NOTE
  Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.
  Please obtain a recent version of HTML Tidy by downloading a binary
  release or compiling the source code from <https://www.html-tidy.org/>.
  ```
  This comes from the old local HTML Tidy, the same as in Stage 1. Because of it, the local check cannot confirm the `<region>` fix; win-builder round 2 will.
- INFO (not counted as a NOTE): installed size 8.6 Mb (R 2.0, doc 4.8, help 1.6).
- Vignette rebuild: `[23s/27s] OK`.
- `00check.log` saved to `dev/reports/cran_wbfollowup_00check_2026-10-09.log`. The temp library and scratch build/check directories were deleted.

## 6. Step 5 — win-builder R-devel upload (one upload)

The check log was committed (`8ac745af`) before the upload. Console output of `devtools::check_win_devel(pkg = ".")`, verbatim:

```
Building windows version of forestsearch (0.4.0)
ℹ Using R-devel with win-builder.r-project.org.
── R CMD build ─────────────────────────────────────────────────────────────────
* checking for file ‘/Users/larryleon/Documents/GitHub/forestsearch/DESCRIPTION’ ... OK
* preparing ‘forestsearch’:
* checking DESCRIPTION meta-information ... OK
* installing the package to build vignettes
* creating vignettes ... OK
* checking for LF line-endings in source and make files and shell scripts
* checking for empty or unneeded directories
Removed empty directory ‘forestsearch/tests/testthat/_snaps’
* building ‘forestsearch_0.4.0.tar.gz’

────────────────────────────────────────────────────────────────────────────────
ℹ Check <'larry.leon.05@post.harvard.edu'> for the results in 15-30 mins
  (~05:37 PM).
```

I did not wait for the result email.

## 7. Step 6 — cran-comments.md (`64283b0d`)

- **Test environments:**
  - local macOS aarch64, R 4.5.2
  - win-builder round 1: R-release 4.6.1 and R-devel r90650
  - win-builder round 2: R-devel, pending
  - mac-builder: unavailable (HTTP 502 ×2, 2026-10-09)
- **R CMD check results:**
  - `0 errors | 0 warnings | 0 notes` expected after the fixes
  - the round-1 HTML-validation NOTE (`inject_mrct_structure.Rd`) is stated as fixed
  - local status quoted verbatim as `Status: 1 NOTE`, the old-Tidy NOTE
- Added one line saying long-running tests are skipped on CRAN and the full suite runs locally.
- `## Submission` and `## Reverse dependencies` are unchanged.

## 8. Post-conditions

- **PC1:** passed, with one addition. The files changed since `50b1b9a5` are:
  - `R/mrct_inject_dgm.R`, `man/inject_mrct_structure.Rd`
  - 7 test files
  - `cran-comments.md`
  - the check log
  - this report
  - **Addition:** the Gate 3 halt note, committed with this report because the Gate 3 amendment cites it.
- **PC2:** passed. `NAMESPACE` is byte-identical to `f43c4299`; `Version: 0.4.0`.
- **PC3:** passed. Each of the 7 test-file diffs is exactly `+skip_on_cran()` (1 line added, 0 removed).
- **PC4:** passed. The tree is clean after the report commit. No tarball, no `.Rcheck`, no temp library in the repo; the scratch `rlib` and `build` directories are deleted.
