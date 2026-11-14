READY FOR WIN-BUILDER

# REPORT — a small core set of `forestsearch()` tests on CRAN (2026-11-13)

Task: `dev/tasks/TASK_cran_core_tests_2026-11-13.md`. Session `fs-cran-core-tests`, host `pop-os`,
branch `master`, version 0.5.0. R 4.6.2 (2026-10-29), testthat 3.3.2. No agents were used. Every test run
used the scratch build `R_LIBS=/tmp/cran10/lib` (`git archive` of `df0d2e44`'s package paths,
`R CMD INSTALL`, 16.3 s). Start SHA `df0d2e44` (the task commit; HEAD before it `37b865fa`). Step 0 gate:
PASS (`master`, `Version: 0.5.0`, `git status --porcelain` empty).

READY because the CRAN-mode suite has 0 failures and took 4.6 s, and the release check has 0 errors and
0 warnings (`Status: OK`, 0 NOTEs).

## Core set

| File | What its CRAN-run tests check | CRAN-mode seconds | Expectations |
|---|---|---|---|
| `test-core-forestsearch.R` (new) | one complete default survival `forestsearch()` call on `survival::gbsg`, sequential plan: the return has the structure documented in `?forestsearch`; a second call with the same seed returns an identical object (timing fields aside) | 1.33 | 42 |
| **Total** | 1 file, 2 blocks | **1.33** | **42** |

Seconds and expectations are from the Step 3 CRAN-mode run. The file runs with no parallel worker and
one thread (CPU 1.27 s for 1.28 s wall; the default call fits no grf forest).

**Why no existing file is kept (Step 1).** Start: the 63 files that ran a block in CRAN mode (round 3 of
the Windows-fixes task, `/tmp/cran9/out`; tests unchanged since). An instrumented CRAN-mode pass traced
`forestsearch()` calls, `plan()` changes, worker launches and grf `num.threads` per test block.
Criterion 1 fails for all 63:

- 53 files make no `forestsearch()` call in a CRAN-run block.
- 10 files call it, but each also runs CRAN blocks that test another function:
  - `test-vi-deprecation.R`: `fs_family_report()`
  - `test-input-validation.R`: `get_FSdata()`, `summarize_cv_results()`
  - `test-glm-pipeline-integration.R`: `grf.subg.harm.glm()`, `get_FSdata()`
  - `test-threshold-pair-directive-a.R`: `.fs_resolve_threshold_pair()`, `make_effect_estimator()`
  - `test-fs-family-report.R`: `fs_family_report()`
  - `test-max-minutes-deprecation.R`: `subgroup.search()`
  - `test-lasso-intercept.R`: `lasso_selection()`, `get_FSdata()`
  - `test-frontier-cis-orientation.R`: `compute_frontier_cis()`
  - `test-fs-declaration-table.R`: `fs_declaration_table()`
  - `test-downstream-consumers.R`: `forestsearch_tenfold()`, `summarize_cv_results()`

Criteria 2 and 4 held for all 63. Every plan set had 1 worker, so future ran it in-process. There was no
PSOCK, fork or callr launch, and grf `num.threads` was 2. CPU/elapsed was at most 1.62 per file, and the
slowest file took 1.59 s. Criterion 3 was not needed.

## The new file

`tests/testthat/test-core-forestsearch.R`, in full. It has no numeric literal. The only digits in it are
the `3` in `expect_s3_class`. Results:

- CRAN mode: 1.28 s.
- `NOT_CRAN=true`, installed build: 1.30 s.
- `NOT_CRAN=true` under `load_all()`, which is not byte-compiled: 3.08 s.
- Negative control: a second call with another `seedit` fails the identity block.

```r
# =============================================================================
# Core test of forestsearch(); runs on CRAN (no skip_on_cran()).
#
# One complete default survival forestsearch() call on survival::gbsg
# (hormonal therapy as the treatment), on a sequential plan.  It checks the
# structure documented in ?forestsearch (Value) and that a second call with
# the same seed returns an identical object.  Structure checks and
# same-session comparisons only: no numeric value is typed into this file.
# The default call fits no grf forest and starts no parallel workers.
# =============================================================================

.core_run <- function() {
  suppressMessages(forestsearch(
    df.analysis      = survival::gbsg,
    outcome.name     = "rfstime",
    event.name       = "status",
    treat.name       = "hormon",
    id.name          = "pid",
    confounders.name = c("age", "meno", "size", "grade", "nodes", "pgr", "er"),
    parallel_args    = list(plan = "sequential")
  ))
}

# The first fit, made once and shared by both blocks.
.core_fit <- local({
  fit <- NULL
  function() {
    if (is.null(fit)) fit <<- .core_run()
    fit
  }
})

# The run-time fields: the only parts of the object that differ between calls.
.core_strip_timing <- function(fs) {
  fs$minutes_all <- NULL
  fs$find.grps$time_search <- NULL
  fs
}

test_that("a default survival forestsearch() call returns the documented structure", {
  fs <- .core_fit()
  dat <- survival::gbsg

  expect_s3_class(fs, "forestsearch")
  expect_type(fs, "list")
  documented <- c("grp.consistency", "find.grps", "confounders.candidate",
                  "confounders.evaluated", "df.est", "df.predict", "df.test",
                  "minutes_all", "grf_res", "sg_focus", "sg.harm", "grf_cuts",
                  "prop_maxk", "max_sg_est", "grf_plot", "args_call_all",
                  "threshold_config", "family_status", "mr_inference",
                  "mr_harm_confirmed")
  expect_identical(setdiff(documented, names(fs)), character())

  # The call as recorded: every formal, the data, the default seed and focus
  expect_identical(setdiff(names(formals(forestsearch)), names(fs$args_call_all)),
                   character())
  expect_identical(fs$args_call_all$df.analysis, dat)
  expect_identical(fs$args_call_all$seedit, eval(formals(forestsearch)$seedit))
  expect_identical(fs$sg_focus, eval(formals(forestsearch)$sg_focus))

  # Thresholds: the survival path compares hazard ratios; the stored
  # screening and consistency values are their logs
  tc <- fs$threshold_config
  expect_identical(setdiff(c("outcome_type", "effect_measure", "screening",
                             "consistency", "screening_natural",
                             "consistency_natural", "scale", "pconsistency"),
                           names(tc)), character())
  expect_identical(tc$outcome_type, "survival")
  expect_identical(tc$screening_natural, fs$args_call_all$hr.threshold)
  expect_identical(tc$consistency_natural, fs$args_call_all$hr.consistency)
  expect_equal(tc$screening, log(tc$screening_natural))
  expect_equal(tc$consistency, log(tc$consistency_natural))
  expect_identical(tc$pconsistency, fs$args_call_all$pconsistency.threshold)

  # The search: candidate factors are columns of the analysis data, and each
  # label of the selected subgroup is one of the evaluated factors
  expect_true(is.list(fs$find.grps))
  expect_type(fs$confounders.candidate, "character")
  expect_type(fs$confounders.evaluated, "character")
  expect_true(all(fs$confounders.candidate %in% names(fs$df.est)))
  expect_type(fs$sg.harm, "character")
  expect_true(all(gsub("^!?\\{|\\}$", "", fs$sg.harm) %in% fs$confounders.evaluated))
  expect_true(is.numeric(fs$max_sg_est))
  expect_true(is.numeric(fs$prop_maxk))
  expect_true(is.numeric(fs$minutes_all))

  # The consistency evaluation
  gc <- fs$grp.consistency
  expect_identical(setdiff(c("out_sg", "sg_focus", "df_flag", "sg.harm", "sg.harm.id",
                             "algorithm", "n_candidates_evaluated", "n_passed"),
                           names(gc)), character())
  expect_identical(gc$sg_focus, fs$sg_focus)
  expect_identical(gc$sg.harm, fs$sg.harm)
  expect_true(gc$algorithm %in% c("twostage", "fixed"))
  expect_true(gc$n_passed <= gc$n_candidates_evaluated)

  # The returned data: one row per patient in input order; the subgroup and
  # its complement are both non-empty, and treat.recommend is the complement
  # of subgroup membership
  expect_s3_class(fs$df.est, "data.frame")
  expect_identical(fs$df.est$pid, dat$pid)
  member <- as.logical(gc$sg.harm.id)
  expect_true(any(member))
  expect_false(all(member))
  expect_true(is.numeric(fs$df.est$treat.recommend))
  expect_identical(as.logical(fs$df.est$treat.recommend), !member)
  expect_identical(gc$df_flag$id, fs$df.est$pid)
  expect_identical(gc$df_flag$treat.recommend, fs$df.est$treat.recommend)
  expect_null(fs$df.test)

  # Defaults off: no front end (no LASSO, grf or DINA), no multiplier resampling
  expect_identical(fs$family_status, "no-front-end")
  expect_null(fs$grf_res)
  expect_null(fs$mr_inference)
  expect_identical(fs$mr_harm_confirmed, NA)
})

test_that("a second call with the same seed returns an identical result", {
  a <- .core_fit()
  b <- .core_run()
  # serialize() round trips drop data.table's session-only self-reference
  # pointers, so the comparison sees values only
  roundtrip <- function(x) unserialize(serialize(x, NULL))
  expect_identical(roundtrip(.core_strip_timing(b)), roundtrip(.core_strip_timing(a)))
})
```

## Files newly skipped on CRAN: 63

Each of the 63 files gains `skip_on_cran()` and a blank line after its header comment, the form already
used. Nothing else in them changes: +126/−0 lines, all in commit `9cf2023b`.

## Step 3 — confirm

1. CRAN mode (`on_cran()` TRUE, `NOT_CRAN` unset, `_R_CHECK_LIMIT_CORES_=TRUE`, test_check-equivalent
   driver on an export of HEAD's `tests/`):
   - 113 files; 1 runs a block (`test-core-forestsearch.R`).
   - 70 blocks: 2 run, 68 skipped.
   - 42 expectations pass (110 counted with the skips).
   - 0 failures, 0 errors, 0 warnings.
   - Wall clock 4.6 s for the driver, 4.9 s from the shell.
2. `NOT_CRAN=true`, the 64 changed files once each (`load_all()` of the HEAD export, sequential plan):
   - 64 files, 381 blocks, 3,489 expectations.
   - 0 failures, 0 errors.
   - 3 skips, all of them references that are not shipped: `test-fs-oc-grid.R`, `test-fs-oc-null.R` and
     `test-threshold-pair-directive-a.R`.
   - 1 min 54 s.

## Step 4 — release check (once)

From `/tmp/cran10/build`: `R CMD build` on the repository, then `R CMD check --as-cran
forestsearch_0.5.0.tar.gz`. `NOT_CRAN` and `OMP_THREAD_LIMIT` were not set.

- **Status: OK.** 0 errors, 0 warnings, 0 NOTEs (no NOTE to quote).
- Tests: `checking tests ... OK` prints no time, because the step is under the 10 s print threshold.
  `testthat.Rout` reads `proc.time()` elapsed 4.741 s (user 4.489 s). The log reads
  `[ FAIL 0 | WARN 0 | SKIP 163 | PASS 42 ]`.
- Vignettes: `checking re-building of vignette outputs ... [78s/77s] OK`.
- Total: check 5 min 39.0 s, build 1 min 44.9 s (shell `time`).
- Tarball `forestsearch_0.5.0.tar.gz`: 1,662,444 bytes. Its `tests/testthat/` is identical to HEAD's.
- Check log: `dev/reports/REPORT_cran_core_tests_2026-11-13_check_log.txt`, a copy of `00check.log`.

`R CMD build` logged "Removed empty directory ‘forestsearch/tests/testthat/_snaps’". That empty,
untracked directory was made in the repository by this session's own `test_file()` runs of the new
file (testthat's snapshot reporter). It was removed with `rmdir` after the check.

## Step 5 — `cran-comments.md`

The sentence about long-running tests is replaced by the task's two sentences. The counts line is from
this check, `0 errors | 0 warnings | 0 notes`, which is unchanged. Nothing else changes (commit
`cc2bf933`).

## Commits (this session, `master`, not pushed)

| Commit | Content |
|---|---|
| `df0d2e44` | dev: task CRAN core tests (as received) |
| `9cf2023b` | test: core CRAN test of forestsearch(); every other CRAN-run file skips on CRAN (Steps 1 and 2) — 64 files under `tests/testthat/` |
| `cc2bf933` | docs: cran-comments.md, the core CRAN test set (Step 5) |
| (this commit) | dev: this report and the check log (Step 6) |

## Final checks (run with this report staged)

1. Paths touched by this session's commits plus the staged report commit: 68; outside the allowed set: none.
   By area: tests/testthat/ 64, cran-comments.md 1, task file 1, report + check log 2.
   Under R/, man/, vignettes/, inst/, DESCRIPTION, NAMESPACE, NEWS.md: 0.
2. DESCRIPTION: Version: 0.5.0.
3. git status --porcelain: 2 entries, all staged additions of the report and check log: yes.
4. Not pushed: origin/master = 37b865fa (last updated 2026-11-13 17:36:04 -0800, before this session); master is ahead 3 commits.
   Merges in this session's commits: 0. Tags on them: 0 (latest tag v0.4.0, 2026-10-09). No pull or merge was run.

Step 6.4 follows this commit: delete `/tmp/cran10/lib`, keep `/tmp/cran10/build`.
