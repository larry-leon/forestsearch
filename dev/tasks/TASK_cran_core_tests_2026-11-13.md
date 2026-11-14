# TASK — run only a small core set of `forestsearch()` tests on CRAN (2026-11-13)

Session name: `fs-cran-core-tests`

Host `pop-os`, repository `~/Documents/GitHub/forestsearch`, branch **`master`**, version 0.5.0.
Read this whole file first.

## Why

Larry's decision (2026-11-13): CRAN runs only a few vital tests of the main function
`forestsearch()`. Every other test runs locally only. The full suite stays as it is and keeps
running locally before each release.

Background: the Windows check ran 678 s (R-devel) on a quiet server before the last task, with
225 s in the tests, and CRAN flags a Windows check above 10 minutes. After that task, 63 test
files still run at least one block on CRAN (24.1 s here); all of them passed on win-builder
R-devel and R-release on 2026-11-13.

## Rules

- **Tests only.** No file under `R/`, `man/`, `vignettes/` or `inst/` changes, and neither does
  `DESCRIPTION`, `NAMESPACE` or `NEWS.md`. No expected value in any existing test is edited,
  and no existing test is removed.
- **Lean run.** The main session does everything. No agents. One short report, plus the raw
  check log as its only companion file.
- **No full test suite with `NOT_CRAN=true`.** Step 4 runs the one release check Larry
  authorized for this task.
- Work on `master`. Commits are path-limited. **No push, pull, merge or tag.**
- No install into the user or site library. Every test run puts a scratch build first:
  `R_LIBS=/tmp/cran10/lib`. Scratch: `/tmp/cran10`.
- Only Step 0 has a stop gate. After it, record a failure and go on.

## Step 0 — start

1. Copy this file into `dev/tasks/` and commit it. Record the start SHA.
2. **Gate (the only stop):** branch `master`; `DESCRIPTION` says `Version: 0.5.0`;
   `git status --porcelain` is empty.
3. `git archive` HEAD's package paths and `R CMD INSTALL` them into `/tmp/cran10/lib`.

## Step 1 — choose the core set

Start from the test files that currently run at least one block in CRAN mode (the last task's
round 3 results in `/tmp/cran9/out`; if they are gone, run the suite once in CRAN mode).

**A file stays on CRAN only if all of these hold:**

1. the tests it runs on CRAN exercise `forestsearch()` itself: the search, the consistency
   evaluation, the returned object, reproducibility, or the main function's argument checks;
2. in CRAN mode it starts no parallel workers and uses at most 2 threads;
3. it compares nothing with a value stored in a file;
4. it takes 2 seconds or less in CRAN mode.

**Size:** the kept files total **5 seconds or less** in CRAN mode. If they total more, drop
the slowest until they fit.

**One end-to-end test is required.** The core set must contain a test that makes a complete
default survival `forestsearch()` call on a small data set, on a sequential plan, and checks
that the result has the documented structure and that a second call with the same seed gives
an identical result. If no kept file has such a test running on CRAN, add it as a new file
`tests/testthat/test-core-forestsearch.R`:

- structure checks and same-session comparisons only; no numeric value typed into the test;
- sequential plan, at most 2 threads, 3 seconds or less;
- it runs on CRAN (no `skip_on_cran()`), and it passes with `NOT_CRAN=true` as well.

## Step 2 — everything else is local only

Add a file-level `skip_on_cran()` (after the header comment, the form already used) to every
other file that currently runs a block on CRAN. Nothing else in those files changes. One
commit for Steps 1 and 2 together.

## Step 3 — confirm

1. Run the suite once in CRAN mode (`testthat::on_cran()` TRUE, `_R_CHECK_LIMIT_CORES_=TRUE`,
   `R_LIBS=/tmp/cran10/lib`). Record: the files that run a block, blocks, expectations,
   failures, and the wall clock. Expected: 0 failures, only the core files run, 6 seconds or
   less in total.
2. Run the changed files once with `NOT_CRAN=true` (changed files only) and confirm that each
   still runs and passes locally.

## Step 4 — the release check (once)

1. From `/tmp/cran10/build`, run `R CMD build` on the repository and
   `R CMD check --as-cran` on the tarball. Do not set `NOT_CRAN` or `OMP_THREAD_LIMIT`.
   Nothing is written into the repository.
2. Record: errors, warnings, every NOTE verbatim; the tests time, the vignette re-build time
   and the total; the tarball size.
3. Copy `00check.log` to `dev/reports/REPORT_cran_core_tests_2026-11-13_check_log.txt`.

## Step 5 — `cran-comments.md`

Under "R CMD check results", replace the sentence about long-running tests with:

> On CRAN a small core set of tests of the main function `forestsearch()` is run. The full
> test suite is run locally before each release.

Set the counts line from this check. Change nothing else in the file. One commit.

## Step 6 — report and close

1. Write `dev/reports/REPORT_cran_core_tests_2026-11-13.md`. Short:
   - first line: READY FOR WIN-BUILDER or NOT READY (ready means: the CRAN-mode suite has 0
     failures and takes 6 seconds or less, and the release check has 0 errors and 0 warnings);
   - **core set table:** each file kept on CRAN — what its CRAN-run tests check (one line),
     CRAN-mode seconds, expectations; and the total;
   - the new file, if one was added, quoted in full;
   - the number of files newly skipped on CRAN;
   - the check result, every NOTE verbatim, the timing lines;
   - this session's commits. Facts only.
2. Commit the report and the check log (path-limited). Copy the report to `~/Downloads`.
3. **Final checks** (report the result): this session's commits touch only files under
   `tests/testthat/`, `cran-comments.md`, this task file, the report and the check log;
   `DESCRIPTION` still says `Version: 0.5.0`; `git status --porcelain` is empty; nothing was
   pushed, merged or tagged.
4. Delete `/tmp/cran10/lib`. Keep `/tmp/cran10/build`.
5. Print the first line and the core set table as the closing message.
