READY FOR WIN-BUILDER

# REPORT — Windows check failures and check time, then one release check (2026-11-13)

Task: `dev/tasks/TASK_cran_windows_fixes_2026-11-13.md`. Session `fs-cran-windows-fixes`, host
`pop-os`, branch `master`, version 0.5.0. R 4.6.2 (2026-10-29), testthat 3.3.2, survival 3.8.12,
Pop!_OS 24.04 LTS. No agents were used. Every test run used the scratch build
`R_LIBS=/tmp/cran9/lib` (`git archive` of `0f1a2730`'s package paths, `R CMD INSTALL`).

READY because: the CRAN-mode suite has 0 failures and took 24.1 s (driver wall clock, as in the
finalize task; 25.4 s from the shell, which adds R start-up), and the release check has 0 errors,
0 warnings and 0 NOTEs (`Status: OK`).

| Step | Result | Commit | Key numbers |
|---|---|---|---|
| 0 Start | gate PASS | `0f1a2730` (task file) | branch `master`, `Version: 0.5.0`, `git status --porcelain` empty; scratch install 15.9 s |
| 1 Censoring fit on Linux | facts below | none | fit does not converge (30 of 30 iterations), does not stop; 12 warnings per build, all "Ran out of iterations and did not converge" |
| 2 `test-ahr-tables.R` | file-level `skip_on_cran()` + 1 comment line | `0d0d3f77` | +3 lines; 3 blocks / 46 expectations pass with `NOT_CRAN=true` |
| 3 Stored/typed references | 26 comparisons to `tolerance = 1e-12` | `1803b163` | the win-builder line (`test-mrct-region-sims-record.R:53`) + 25 lines in 7 files still run on CRAN |
| 4 Check time | 17 files skipped on CRAN, 3 rounds | `8ec06c42` | CRAN-mode total 102.6 s (projected after Step 2) → 24.1 s measured; 0 failures; 25 changed files pass with `NOT_CRAN=true` |
| 5 Release check | **Status: OK**: 0 errors, 0 warnings, 0 NOTEs | `47d5f728` (`cran-comments.md`) | tests 26 s / 26 s; vignette re-build 76 s / 75 s; check 5 min 52.4 s; build 1 min 44.8 s; tarball 1,661,046 bytes |
| 6 Report | this file + check log | (this commit) | |

## Step 1 — the censoring fit on Linux (facts)

Probe: `.at_dgm(3)` as in the test file, with `suppressWarnings()` replaced by a recorder and a
`trace()` exit hook on `prepare_censoring_model()` (in memory only; `R/` unchanged).

1. The fit does not converge and does not stop: `survreg(cens_formula, data = df_work, dist = "weibull")` ran 30 of 30 iterations (`maxiter` 30); so did the other three AIC candidates (lognormal, and the two `~ 1` fits).
2. Warnings: 12 per build, every one `Ran out of iterations and did not converge` from `survreg.fit(...)`: 3 `prepare_censoring_model()` calls (Steps 7 and 7b of `generate_aft_dgm_flex()`, and `inject_mrct_structure()`'s refit) × 4 fits. No other warning. `.at_dgm(2)` gives the same 12 (the censoring fit does not use `h`).
3. `df_work`: 500 rows; censoring events `sum(1 - event)` = 128. `cens_formula`: `Surv(y, 1 - event) ~ treat + zcens_age + zcens_lab + zcens_f1` (not standardised: age mean 49.4, sd 9.7; lab mean 299.3, sd 95.2; treat and f1 binary).
4. The Weibull fit with covariates: (Intercept) 2.277, treat 0.375, zcens_age 0.00557, zcens_lab −0.000475, zcens_f1 0.144; scale 0.0216; loglik −167.09 (intercept only, survreg's own) and −60,725.5 (with covariates).
5. The AIC choice is Weibull0 (`~ 1`): intercept 3.6859 (exp = 39.88), scale 0.00426, loglik 56.28. `compare_multiple_survreg()` prints "Models converged: 4 / 4".
6. Ill-conditioning: all 128 censoring events are at one time, y = 40 (the fixture's `pmin(tte, 40)`; 1 distinct value). The censoring distribution is a point mass, the likelihood keeps rising as the scale goes to 0, and every fit stops at the limit (scale 0.004 to 0.035). The covariate fit's final loglik is far below the intercept-only value. No covariate is near-constant; 128 events is not few.
7. Same message as the simulations article. Its own-trial build `build_own(h = 2)` (chunk `own-build`, `rep5_trial.csv`) shows it 3 times. There it comes from the same Weibull-with-covariates fit, at 30 iterations, once in each of the 3 calls; the other three fits converge in 9–18 iterations. Censoring times there are spread (595 distinct, 29.0 to 40.0), and the AIC choice is LogNormal0.
8. The application article builds the same `rep5_trial.csv` DGM (`mrct-consistency.qmd:287`) and prints the count (`build_warnings` 3) rather than the messages. No article builds from the fixture's configuration.
9. Same configuration (seed data cut at 40 by `tte < 40` and `pmin(tte, 40)`, then `inject_mrct_structure()`): the examples of `ahr_truths_table()`/`ahr_design_table()` (`man/ahr_tables.Rd`), `mrct_build_dgm()` (`man/mrct_build_dgm.Rd`) and `mrct_core_summary()` (`man/mrct_core_summary.Rd`, in `\donttest{}`), all `set.seed(1)`, n = 400; the tests `test-mrct-trial-analysis.R` (`.mta_fast_dgm()`, n = 500, seed 101, five covariates) and `test-mrct-core-summary.R:114` (n = 400, seed 1, already `skip_on_cran()`).

## Files skipped on CRAN by this task

All are a top-level `skip_on_cran()` after the file's header comment (+2 lines; +3 with the comment
line in Step 2). Times are CRAN-mode wall clock per file: the finalize task's run
(`/tmp/cran8/out/suite_cran.rds`), or this task's round 1/2 where marked.

| # | File | Seconds | Expectations | Step |
|---|---|---|---|---|
| 1 | `test-ahr-tables.R` | 0.5 | 46 | 2 (Windows error) |
| 2 | `test-mrct-region-sims-record.R` | 21.4 | 172 | 4 |
| 3 | `test-mr-estimates-table-nullboot.R` | 9.8 | 33 | 4 |
| 4 | `test-betaHhat-exact.R` | 8.7 | 38 | 4 |
| 5 | `test-ahr-analytic.R` | 8.3 | 208 | 4 |
| 6 | `test-dfbeta-lm-influence-path.R` | 4.3 | 31 | 4 |
| 7 | `test-cross-outcome-parity.R` | 3.3 | 38 | 4 |
| 8 | `test-threshold-remap-messages.R` | 3.3 | 27 | 4 |
| 9 | `test-missing-membership-messages.R` | 3.0 | 60 | 4 |
| 10 | `test-mrct-region-interaction.R` | 2.9 | 201 | 4 |
| 11 | `test-get-fsdata-diagnostics.R` | 2.6 | 63 | 4 |
| 12 | `test-betaHhat-contract.R` | 2.3 | 154 | 4 |
| 13 | `test-mr-closedform-fixture.R` | 2.2 | 41 | 4 |
| 14 | `test-subgroup-sims.R` | 2.1 | 78 | 4 |
| 15 | `test-mrct-trial-analysis.R` | 1.8 | 122 | 4 |
| 16 | `test-mrct-minor-helpers.R` | 1.6 | 27 | 4 |
| 17 | `test-mrct-dgm-args.R` | 1.8 (round 1; 0.8 in finalize) | 47 | 4, round 2 |
| 18 | `test-fs-oc-predict.R` | 1.6 (round 2) | 89 | 4, round 3 |

- `test-mrct-dgm-args.R` was 0.8 s in the finalize run. Once the earlier files were skipped, it was
  the only CRAN-run file that built the memoised S2 DGM (`helper-mrct-analytic.R`), so the build
  cost moved to it (1.8 s).
- Projection (finalize per-file times): 103.1 s (finalize, after its two skips) − 0.5 (Step 2) =
  102.6 s; − files 2–16 = 24.8 s.

## Step 4 — the CRAN-mode runs

Driver: the finalize task's `suite_run.R` (what `tests/testthat.R` → `test_check()` does, plus a
per-file timer), run from an export of the working tree's `tests/` with `NOT_CRAN` unset,
`_R_CHECK_LIMIT_CORES_=TRUE`, `R_LIBS=/tmp/cran9/lib`. Reporter evidence at start and end:
`on_cran()` TRUE, `NOT_CRAN` unset, `interactive()` FALSE.

| Round | Skipped (cumulative) | Files | Blocks | Expectations | Failed | Errors | Skipped blocks | Warnings | Driver wall | Shell real |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 15 + Step 2 | 112 | 465 | 3,477 | 0 | 0 | 107 | 8 | 27.6 s | 28.8 s |
| 2 | + `test-mrct-dgm-args.R` | 112 | 459 | 3,430 | 0 | 0 | 107 | 8 | 26.0 s | 27.2 s |
| 3 (final) | + `test-fs-oc-predict.R` | 112 | 447 | 3,341 | 0 | 0 | 107 | 8 | **24.1 s** | 25.4 s |

- CRAN-mode total: **102.6 s before** (projected, after Step 2), **24.1 s after** (round 3, driver
  wall clock, the measure the finalize task used; the shell figure adds R start-up). The finalize
  report's last measured figure was 103.1 s.
- Still run on CRAN after round 3: 63 of 112 files run at least one block (340 blocks, 3,234
  expectations in them); 32 files are skipped at file level ("On CRAN"; 14 before this task); 17
  files run no block for other reasons (in-block skips). Skip records: 126 "On CRAN", 9 threshold
  probe unavailable, 2 pre-refactor reference not present, 1 baseline CSV not in this tree, 1
  package source not available.
- Slowest file left: `test-mr_values.R`, 1.6 s.
- The tests run in round 3 are byte-identical to HEAD `8ec06c42`'s `tests/` (`diff -r`).
- Step 4.4, `NOT_CRAN=true`, the 25 files changed in Steps 2–4, each once (`load_all()` of the HEAD
  export, sequential plan, finalize's `run_tests.R`): 25 files, 195 blocks, 1,743 expectations,
  0 failures, 0 errors, 5 skipped blocks (3 "simulation modules not present" in
  `test-betaHhat-contract.R`, 2 "primary fit did not identify; test N/A" in
  `test-cross-outcome-parity.R`; none "On CRAN"; `on_cran()` FALSE); 2 min 49 s.

## Step 3 — comparisons changed to `tolerance = 1e-12`

Each `expect_identical(a, b[, label])` became `expect_equal(a, b[, label], tolerance = 1e-12)`; no
expected value changed.

- 3.1: `test-mrct-region-sims-record.R:53` (`pp` against the stored fixture `point_6b5c85c0.rds`).
  Line 46 (`pp` against the call made in the same session) stays `expect_identical`.
- 3.2 survey: an extra CRAN-mode run of the suite (Steps 2 + 4 + 3.1 applied, 33 s) with a tracer on
  `expect_identical()` recorded 1,290 calls at 454 source lines; 163 lines compared a double. Each
  was read. Static greps of the running files for `tolerance = 0`, `all.equal(..., tolerance = 0)`,
  `expect_true(identical(...))`, `== <decimal>` and hashed values added one line,
  `test-mr-method-tables.R:91` (`tolerance = 0`, same-session reference, left).
- Changed (25): a double computed by the package or the test, against a value typed into the test.
  No comparison run on CRAN takes its reference from a file: the three that would
  (`test-fs-oc-grid.R`, `test-fs-oc-null.R`, `test-threshold-pair-directive-a.R`, all reading `dev/`)
  skip because `dev/` is not in the tarball.

| File | Lines |
|---|---|
| `test-cv-mixed-and-fatal.R` | 125, 142 (`find_metrics` rows, `rep(0, 8L)`) |
| `test-fs-dgm-feasibility.R` | 59, 62, 93, 119 (shares 1, 0, 0, 0) |
| `test-fs-oc-signed.R` | 90 (`m_tau` 0), 103 (`PQg` 0), 104 (`beta_g` 26) |
| `test-mr-method-tables.R` | 97 (30), 98 (1.5), 99 (750) |
| `test-mr-residual-centering.R` | 98 (1) |
| `test-mrct-core-summary.R` | 94 (power estimate 0) |
| `test-npv-empty-complement.R` | 63, 65, 67, 70 (`50 / 200`), 71, 108, 128, 129, 130, 132, 134 |

- Left: same-session pairs; counts and denominators; signs, codes and indicators
  (`.dina_tau_sign()`, `.grf_code_column()`, `treat.recommend.original == 1.0`); inputs passed
  through (`tab$c0`, `max.minutes`, `hazard_multiplier`, `bundle$C` after `saveRDS()`, the
  threshold pair's supplied `c2`). `test-cv-mixed-and-fatal.R:198`
  (`c(NaN, NA, 0, NA)`) is left because `expect_equal(..., tolerance = )` does not tell `NaN` from
  `NA` (checked), and that difference is what the line tests.

## Step 5 — release check

- `R CMD build ~/Documents/GitHub/forestsearch` from `/tmp/cran9/build` (HEAD `8ec06c42`), then
  `R CMD check --as-cran forestsearch_0.5.0.tar.gz`. `NOT_CRAN` and `OMP_THREAD_LIMIT` unset;
  `R_LIBS=/tmp/cran9/lib`. Nothing was written into the repository (`git status --porcelain` showed
  only the Step 5 edit afterwards).
- Tarball 1,661,046 bytes (finalize: 1,660,903); no `legacy_vignette` or `Rplots.pdf` entry; its
  `tests/` equal HEAD's (`diff -rq`). Installed size 4,760 KB.
- NOTEs: **none**. `Status: OK`. The companion file
  `REPORT_cran_windows_fixes_2026-11-13_check_log.txt` is this check's `00check.log`.

Timing lines (verbatim from `00check.log`):

```
* checking CRAN incoming feasibility ... [6s/15s] OK
* checking whether package ‘forestsearch’ can be installed ... [19s/19s] OK
* checking R code for possible problems ... [44s/44s] OK
* checking examples ... [34s/34s] OK
* checking examples with --run-donttest ... [60s/61s] OK
* checking tests ... [26s/26s] OK
  Running ‘testthat.R’ [26s/25s]
* checking re-building of vignette outputs ... [76s/75s] OK
* checking PDF version of manual ... [13s/13s] OK
* checking HTML version of manual ... [18s/17s] OK
```

| Item | This check | Finalize run 2 |
|---|---|---|
| Tests (`testthat.R`) | 26 s / 25 s | 93 s / 103 s |
| Vignette re-build | 76 s / 75 s | 74 s / 73 s |
| `R CMD check` total (shell `time`) | 5 min 52.4 s | 6 min 56.4 s |
| `R CMD build` (shell `time`) | 1 min 44.8 s | 1 min 41.0 s |

- `cran-comments.md` (`47d5f728`): the counts line stays `0 errors | 0 warnings | 0 notes` (this
  check's counts); one line added under "R CMD check results": "Long-running tests are skipped on
  CRAN with `skip_on_cran()`; the full suite runs locally." The win-builder and R-hub lines are
  unchanged.

## Order and departures

- Step 4 was done before Step 3: Step 3.2 surveys the files that still run on CRAN, which needs
  Step 4's skips. The Step 3 commit comes first in history. `test-mrct-region-sims-record.R` has a
  hunk in each: the Step 3 hunk was staged alone (`git update-index`), the skip went with Step 4.
- One extra CRAN-mode run (33 s, with the `expect_identical()` tracer) was made for the Step 3.2
  survey, besides the three timing rounds. No `NOT_CRAN=true` full-suite run was made.
- Step 4.3 needed two of its three allowed extra rounds.
- `/tmp/cran8/out` still held the finalize per-file list, so the suite was not re-timed before the
  rule.

## Commits (this session)

Start SHA: `0f1a2730c7e863ed5473131bb109f422ae97b42b`, the task-file commit; its parent `418b172d` was
HEAD at the start.

| SHA | Commit |
|---|---|
| `0f1a2730` | dev: task CRAN Windows fixes (as received) |
| `0d0d3f77` | test: skip test-ahr-tables.R on CRAN (Windows fixes, Step 2) |
| `1803b163` | test: tolerance 1e-12 on comparisons with stored or typed references (Windows fixes, Step 3) |
| `8ec06c42` | test: skip the slowest CRAN-mode files on CRAN (Windows fixes, Step 4) |
| `47d5f728` | docs: cran-comments.md, check results and the skip_on_cran() line (Windows fixes, Step 5) |
| (this) | dev: CRAN Windows fixes report and check log (Step 6) |

## Final checks

Run with the report and the check log staged, before the report commit; re-run after it (result
in the closing message).

```
1. paths touched by this session's commits (5 commits since 418b172d, plus the staged report commit): 29 files; outside the allowed set: none
   tests/testthat/ 25, cran-comments.md 1, task file 1, report + check log 2
2. DESCRIPTION: Version: 0.5.0
3. git status --porcelain: only the two staged report files
4. merges since 418b172d: 0; tags containing 0f1a2730: 0; remote-tracking branches containing 0f1a2730: 0; branch master
RESULT: PASS
```

Scratch: `/tmp/cran9/lib` deleted after the final checks; `/tmp/cran9/build` kept (tarball and
`.Rcheck`), with `/tmp/cran9/out` (drivers, logs, `.rds` results), `/tmp/cran9/src`,
`/tmp/cran9/run_final` and `/tmp/cran9/run_survey`.
