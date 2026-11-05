# TASK_fs-estimability-fix_2026-11-04 — reconcile the estimability-boundary invariant test

**Repo**: forestsearch (Pop!_OS box). **Authorization**: Larry, 2026-11-04 ("proceed with A and B"); this document is the task of record.
**First action**: copy this file into `dev/tasks/` and commit that copy alone, before any other change.
**Git**: `export GIT_TERMINAL_PROMPT=0` at session start. Never push; no GitHub credentials. Short-lived branch from `master`, fast-forwarded, deleted (check `git branch -r` first for GitHub-Desktop auto-published copies).
**R/ callout (binding)**: this task is expected to touch `tests/` (plus `DESCRIPTION`/`NEWS.md` and `dev/tasks/`) ONLY. It must not change the behaviour of anything under `R/`. If the investigation concludes the correct fix is in `R/`, STOP at the decision gate below — report; do not edit.
**Amendments**: any chat amendment is appended verbatim to the committed copy as "Amendment N" in its own commit.

## Context

During TASK_fs-dgm-region-interaction_2026-11-04, the full suite surfaced one pre-existing failure, unrelated to that task's diff:

- `tests/testthat/test-estimability-boundary.R:163`
- Suspected origin: commit `f6735639`, which added a second `isFALSE(result$converged)` condition to `R/glm_effect_estimators.R`.

Larry waived the failure to allow that task's fast-forward and has now authorized this fix. The failure blocks a green suite and therefore the 0.4.2 release preconditions.

The coordinates above come from the prior task's report. Re-confirm every one of them from source. If the failure does not reproduce, or the line/commit differ, record the actuals in NOTES and proceed against the actuals.

## Objective

A fully green `devtools::test()` via a test-side reconciliation: the invariant test updated to describe the *intended* behaviour of `R/glm_effect_estimators.R`, with that intent established from source before any edit is made.

## Steps

1. **Baseline.** Run `devtools::test()`; record the failure verbatim (file:line, expectation, observed). Run `devtools::check()`; record baseline errors/warnings/notes for later comparison.
2. **Read the sources.** The failing expectation and its setup in `tests/testthat/test-estimability-boundary.R`; the function(s) it exercises in `R/glm_effect_estimators.R`.
3. **Read the originating commit.** `git log -1 --format='%H %ad %s' f6735639` and `git show f6735639 -- R/glm_effect_estimators.R`. Classify the second `isFALSE(result$converged)`:
   (i) literally redundant — same object, same branch, a no-op duplicate (then the failure has another cause: find it before touching anything);
   (ii) semantically distinct — different object or branch, a deliberate second check;
   (iii) a behaviour change — results now classified differently than the test's invariant assumes.
4. **Root cause first.** Before any edit, write to NOTES one paragraph naming the exact mechanism by which the line-163 expectation fails, with the commit evidence.
5. **Decision gate.**
   - **Branch A — the `R/` behaviour is intentional and correct** (the invariant predates f6735639 and no longer describes intended behaviour): update the TEST ONLY — the minimal edit that makes the invariant describe intended behaviour — with a comment in the test citing f6735639. Continue.
   - **Branch B — the `R/` code is wrong** (accidental duplicate with semantic effect, or any other defect): **STOP.** No edit to `R/` and no edit to the test. Write findings and a proposed minimal `R/` fix to NOTES, leave the branch with only the task-doc commit, fast-forward nothing beyond it, and report. An `R/` behaviour change is proposal-only per protocol; it would be a new task.
6. *(Branch A only)* Bump `Version:` to `0.4.1.9002`; `NEWS.md` entry under the dev version: test-suite reconciliation of the estimability-boundary invariant with `glm_effect_estimators` behaviour; no functional change.
7. Full `devtools::test()` and `devtools::check()`.

## Gates (stop on first failure; machine-checkable)

| Gate | Check |
|---|---|
| G0 | Task doc committed alone into `dev/tasks/` as this task's first commit |
| G1 | Baseline failure and check baseline recorded in NOTES before any edit |
| G2 | Root-cause statement in NOTES, with f6735639 evidence, before any edit |
| G3 | Decision gate recorded (Branch A or B). Branch B ⇒ stop and report |
| G4 | `git diff --stat master...HEAD` touches only `tests/`, `DESCRIPTION`, `NEWS.md`, `dev/tasks/` |
| G5 | `devtools::test()`: 0 failures (pre-existing skips allowed) |
| G6 | `devtools::check()`: 0 errors, 0 warnings; notes no worse than the G1 baseline |

On green: fast-forward `master`, delete the branch, report with the standard table (gate | status | evidence). A green report needs no chat round trip. NOTES to `~/Downloads/NOTES_fs-estimability-fix_2026-11-04.md`.
