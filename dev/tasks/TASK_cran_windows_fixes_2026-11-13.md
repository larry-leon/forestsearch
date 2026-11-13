# TASK — Windows check failures and check time, then one release check (2026-11-13)

Session name: `fs-cran-windows-fixes`

Host `pop-os`, repository `~/Documents/GitHub/forestsearch`, branch **`master`**, version 0.5.0.
Read this whole file first.

## Background

win-builder checked the 0.5.0 tarball on 2026-11-13 (R-release 4.6.2 and R-devel r90881). Both
ended `Status: 1 ERROR`, with the same four test failures and nothing else:

1. `test-ahr-tables.R:31`, `:47`, `:68` (three errors). `.at_dgm(3)` →
   `inject_mrct_structure()` → `generate_aft_dgm_flex()` → `prepare_censoring_model()` →
   `survival::survreg(cens_formula, data = df_work, dist = "weibull")` stops with
   `initial iteration failed (use starting estimates?)`. On Linux the same call does not stop.
2. `test-mrct-region-sims-record.R:53` (one failure).
   `expect_identical(as.list(pp), as.list(fx))` differs in the last bit on Windows
   (`hr_train[3]` 0.82157082749266153 against 0.82157082749266142; `POhr_sg[4]`
   0.55000000000000004 against 0.55000000000000016).

Check time on Windows: 678 s (R-devel) and 741 s (R-release), of which the tests took 225 s and
230 s. CRAN's incoming check flags an overall check time above 10 minutes on Windows. The last
local figure for the tests as CRAN runs them is 103 s (finalize report), so Windows runs them
about 2.2 times slower.

## Rules

- **Tests only.** No file under `R/`, `man/`, `vignettes/` or `inst/` changes, and neither does
  `DESCRIPTION`, `NAMESPACE` or `NEWS.md`. No expected value in any test is edited.
- **Lean run.** The main session does everything. No agents. One short report, plus the raw
  check log as its only companion file.
- **No full test suite with `NOT_CRAN=true`.** Step 5 runs the one release check Larry
  authorized for this task.
- Work on `master`. Commits are path-limited, one per step. **No push, pull, merge or tag.**
- No install into the user or site library. Every test run puts a scratch build first:
  `R_LIBS=/tmp/cran9/lib`. Scratch: `/tmp/cran9`.
- Only Step 0 has a stop gate. After it, record a failure and go on.

## Step 0 — start

1. Copy this file into `dev/tasks/` and commit it. Record the start SHA.
2. **Gate (the only stop):** branch `master`; `DESCRIPTION` says `Version: 0.5.0`;
   `git status --porcelain` is empty.
3. `git archive` HEAD's package paths and `R CMD INSTALL` them into `/tmp/cran9/lib`.

## Step 1 — what the censoring fit does on Linux (read-only facts)

Run the fixture of `tests/testthat/test-ahr-tables.R` (`.at_dgm(3)`) on this machine with
warnings shown, not suppressed, and look at the `survreg` fit inside
`prepare_censoring_model()`. Report, in at most 15 lines:

- whether the fit converges, the number of iterations, and every warning verbatim;
- the number of rows, the number of censoring events, the terms of `cens_formula`, and the
  fitted coefficients and scale;
- anything that makes the fit ill-conditioned (a near-constant covariate, very few events,
  coefficients running off);
- whether this is the same warning that the `mrct-consistency-simulations` article shows
  ("Ran out of iterations and did not converge"), and which articles or examples use this
  configuration.

Facts only. Change nothing in `R/`, and propose nothing.

## Step 2 — `test-ahr-tables.R`: skip on CRAN

The fixture's censoring fit is platform-sensitive and stops on Windows. Add `skip_on_cran()`
so that the tests that build `.at_dgm()` do not run on CRAN (a file-level skip if every test
in the file depends on it; otherwise in each such test). Add one comment line giving the
reason and the win-builder date. Nothing else in the file changes. One commit.

## Step 3 — comparisons with stored references: allow last-bit differences

1. `test-mrct-region-sims-record.R:53`: the comparison of `pp` with the stored reference `fx`
   becomes `expect_equal(..., tolerance = 1e-12)`. The comparison with the call made in the
   same session stays `expect_identical`.
2. After Step 4, look through the test files that **still run on CRAN** for the same kind of
   comparison: a floating-point result required to be `identical` (or equal with zero
   tolerance) to a value stored in a file or typed into the test. Change each to
   `tolerance = 1e-12` and list it. Comparisons of two values computed in the same session
   are left alone.
3. One commit. No expected value changes.

## Step 4 — check time: fewer tests on CRAN

1. Take the per-file times of the suite as CRAN runs it from the finalize task
   (`/tmp/cran8/out`, or `dev/reports/REPORT_cran_prep_finalize_2026-11-13.md`). If the
   per-file list is no longer available, time the suite once in CRAN mode, as that task did.
2. **Larry's rule for this task:** add a file-level `skip_on_cran()` to the slowest files, one
   by one, until the projected CRAN-mode total is **25 seconds or less**. Nothing else in
   those files changes.
3. Run the suite once in CRAN mode (`testthat::on_cran()` TRUE, `_R_CHECK_LIMIT_CORES_=TRUE`,
   `R_LIBS=/tmp/cran9/lib`). Record files, blocks, expectations, failures, skips and the wall
   clock. If the total is still above 25 s, skip the next slowest file and run again, at most
   three rounds.
4. Run the files changed in Steps 2–4 once with `NOT_CRAN=true` and confirm that each still
   runs and passes locally (changed files only, not the whole suite).
5. One commit. In the report: the files skipped by this step with their times, and the
   CRAN-mode total before and after.

## Step 5 — the release check (once)

1. From `/tmp/cran9/build`, run `R CMD build` on the repository and
   `R CMD check --as-cran` on the tarball. Do not set `NOT_CRAN` or `OMP_THREAD_LIMIT`.
   Nothing is written into the repository.
2. Record: errors, warnings, every NOTE verbatim; the tests time, the vignette re-build time
   and the total; the tarball size.
3. Copy `00check.log` to `dev/reports/REPORT_cran_windows_fixes_2026-11-13_check_log.txt`.
4. `cran-comments.md`, only these edits, one commit:
   - "R CMD check results": the counts from this check;
   - add one line: long-running tests are skipped on CRAN with `skip_on_cran()`; the full
     suite runs locally;
   - the win-builder and R-hub lines stay as they are.

## Step 6 — report and close

1. Write `dev/reports/REPORT_cran_windows_fixes_2026-11-13.md`. Short:
   - first line: READY FOR WIN-BUILDER or NOT READY (ready means: the CRAN-mode suite has 0
     failures and takes 25 s or less, and the release check has 0 errors and 0 warnings);
   - top table: one row per step — result, commit, key numbers;
   - Step 1's facts;
   - the files skipped on CRAN by this task, and how many files and expectations still run on
     CRAN;
   - Step 3's list;
   - every NOTE verbatim and the timing lines;
   - this session's commits. Facts only.
2. Commit the report and the check log (path-limited). Copy the report to `~/Downloads`.
3. **Final checks** (report the result): this session's commits touch only files under
   `tests/testthat/`, `cran-comments.md`, this task file, the report and the check log;
   `DESCRIPTION` still says `Version: 0.5.0`; `git status --porcelain` is empty; nothing was
   pushed, merged or tagged.
4. Delete `/tmp/cran9/lib`. Keep `/tmp/cran9/build`.
5. Print the first line and the top table as the closing message.
