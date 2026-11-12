# TASK — empty identified set: detection flag, PPV, and the CV agreement metrics, 2026-11-11

Complete on its own. Run from the forestsearch repo root. Authorized by Larry 2026-11-11.

Supersedes `TASK_npv_nofind_instances_2026-11-11_v2.md` (do not run that one). The NPV
part of v2 is already done on branch `fix/npv-nofind-instances`
(`dev/reports/REPORT_npv_nofind_instances_2026-11-11.md`). This task does the rest, on
that same branch.

## Purpose

- **Kind of change: changes behaviour. It does not change the method.**
- The definition (León et al. 2024, Section 3): when no subgroup is identified, the
  estimated complement is everyone. Under an alternative that replicate has `sens` = 0,
  `ppv` = 0 (convention for 0/0), `spec` = 1, `npv` = |Hc|/N. Under the null: `spec` = 1,
  `npv` = 1, and `sens` and `ppv` stay NA, as the no-find path stores them now.
- Two places still depart from it:
  1. **Part 1.** In `R/run_simulation_analysis.R`, a replicate can have an empty
     identified set without taking the no-find path. It is then flagged as detected, and
     its `ppv` is NA where the no-find path stores 0
     (`dev/reports/REPORT_a003_fix_2026-11-11.md`, open items 3 and 4).
  2. **Part 2.** In `forestsearch_KfoldOut()`, when the folds identify nothing and the
     full-data analysis has a subgroup, the agreement sensitivity for H is reported as 1
     for zero agreement (`REPORT_npv_nofind_instances_2026-11-11.md`, open item 4).

## Rules

- No `git fetch`, `pull`, `push` or `merge`. Every `git add` names explicit paths.
- No package install into any library. No render, no simulation or application re-run.
- `quarto/`, `vignettes/` and the two sibling repositories: not touched.
- Evidence is the current source. Quote file:line.
- Gates are stop-on-failure, never stop-to-ask. On a failure, write what failed into the
  report, commit the report, and stop.
- Nothing else changes: no refactor, no renaming, no other fix.

## Step 0 — branch and record

1. `git status --porcelain` is empty. **Gate.**
2. Working branch: check out `fix/npv-nofind-instances`. If it no longer exists (merged
   and deleted), create a branch of that name from `master`.
3. **Gates on the working branch:**
   - `grep -cE 'npv +<- +sum\(!true_H\) / nrow\(df\)' R/run_simulation_analysis.R`
     returns 4 or more;
   - `grep -c 'hc_share' R/fs_mr_oc_summary.R` returns 1 or more.
4. Copy `~/Downloads/TASK_empty_set_nofind_2026-11-11.md` to `dev/tasks/` and commit it on
   the working branch (`dev: add empty-set no-find task`).
5. Record the working branch, its HEAD SHA before step 4 (the start SHA), and the
   `DESCRIPTION` version (expected 0.4.1.9005).

## Part 1 — an empty identified set is a no-find on every path (`R/run_simulation_analysis.R`)

In `.extract_fs_estimates_gen` and `.extract_grf_estimates_gen`, two known cases:

- FS "site 2": the result has rows but no `M.*` factor is recovered. `any.H` was set to 1
  earlier (near `:1041`) and is not reset, while the block stores size 0 and no-find
  classification.
- The found path with an identified set of size 0 (e.g. rows present, `treat.recommend`
  all 1). The found-case formula gives `ppv` = 0/0 = NA.

**Rule E.** Whenever no subject is in the identified set, by whatever path, the replicate
is scored exactly as the no-find path scores it for the same data: `any.H` = 0, the same
sizes, the same classification values, and the same values for every other per-replicate
output.

- First trace, from the source, how each case can arise in the real pipeline, and put the
  trace in the report.
- **Guard.** If the trace shows that a case can be reached while a non-empty identified
  set exists for that replicate (membership available from another object), do not change
  that case: report it as `BLOCKED` with the trace.
- Make the smallest edit that satisfies Rule E. A replicate with a non-empty identified
  set must give exactly the outputs it gives now.
- List every output whose value changes for an empty-set replicate: name, old value, new
  value.

## Part 2 — CV agreement metrics when the folds identify nothing (`R/forestsearch_cross_validation.R`)

Site: the single-level branch of `forestsearch_KfoldOut()` (near `:1409-1420`), where the
table of fold assignment against full-data assignment has only one level of the fold
assignment.

**Rule C.** When the full-data assignment has both levels, compute the four agreement
values from their two-level definitions with a missing fold level counted as zero
subjects, and with 0/0 = 0:

| Fold assignment | `sens_H` | `ppv_H` | `sens_Hc` | `ppv_Hc` |
|---|---|---|---|---|
| Nobody in H (folds identify nothing) | 0 | 0 | 1 | share of subjects in the full-data complement |
| Everybody in H | 1 | share of subjects in the full-data H | 0 | 0 |

- When the full-data assignment itself has a single level, leave the code as it is.
- Use the source's own names and output order for the four values.
- Update the roxygen or comments at the site if they state the old values. Do not touch
  any other CV code.

## Common to both parts

- If any roxygen changed, run `devtools::document()` once.
- `DESCRIPTION`: 0.4.1.9005 → 0.4.1.9006.
- `NEWS.md`, development section, one entry per part:
  - Part 1: `run_simulation_analysis()`: a replicate whose identified set is empty is now
    scored as a no-find on every path (`any.H` = 0; `ppv` = 0 under an alternative). The
    detection rate and the all-replicate PPV average can change where such replicates
    occur.
  - Part 2: `forestsearch_KfoldOut()`: when the folds identify no subgroup and the
    full-data analysis has one, the agreement metrics are now 0, 0, 1 and the full-data
    complement share (the sensitivity for H was reported as 1).

**Gates after the edits** (diffs are against the start SHA):

1. `git diff <start> --stat -- R/` lists only `R/run_simulation_analysis.R` and
   `R/forestsearch_cross_validation.R`.
2. Every hunk in `git diff <start> -- R/` maps to Rule E or Rule C; the report gives the
   mapping.
3. If `document()` ran: `git diff <start> --stat -- man/` lists only the Rd files of
   functions whose roxygen changed.
4. The two Step 0 greps still pass, and
   `grep -c 'spec <- sum(!true_H)' R/run_simulation_analysis.R` returns 0.

## Part 3 — tests

- Part 1, in `tests/testthat/test-nofind-classification.R`, with hand-built inputs:
  1. FS site 2 input: `any.H` is 0 and every output equals the site 1 output for the same
     data;
  2. found path with an empty identified set, FS and GRF: every output equals the no-find
     output for the same data (`ppv` 0 under an alternative; NA under the null);
  3. found path with a non-empty identified set, FS and GRF: `any.H` is 1 and the four
     classification values equal the found-case formulas computed inline.
- Part 2:
  1. update `tests/testthat/test-cv-mixed-and-fatal.R:142` (it pins the old values) to the
     Rule C values; list old and new;
  2. add both rows of the Rule C table as expectations, and one two-level case whose four
     values are unchanged from the start SHA.
- If a case cannot be built by hand without running a search, do not run one: say so.
- No search, no simulation, no parallel plan in any new test. Hard cap: each test file
  runs in under 2 minutes.

## Part 4 — tests and check

1. Run the changed test files. **Gate:** 0 failures.
2. Then, once: the full suite and `devtools::check()`, on the same check route as
   `dev/reports/REPORT_npv_nofind_instances_2026-11-11.md` §6.
   - Baseline, from that report (do not re-run it): suite 608 blocks, 81 files, 8,218
     expectations, 0 failures, 3 skips, 38 warnings; check 0 errors / 0 warnings /
     0 notes.
   - Hard caps: suite 30 minutes, `check()` 60 minutes. On a cap, stop the run, record the
     wall clock, and stop the task.
   - Pre-set decision: a failing expectation that pins an old value this task changes on
     purpose (`any.H` = 1 or `ppv` = NA for an empty identified set; the old CV values in
     the single-level branch) is updated and listed. **Gate:** any other failure, or any
     new check ERROR, WARNING or NOTE → stop.
3. Commit on the working branch by explicit path
   (`fix: empty identified set is a no-find on every path; CV agreement when folds find nothing`).

## Part 5 — report and close

1. Write `dev/reports/REPORT_empty_set_nofind_2026-11-11.md`, in bullets and tables:
   - provenance;
   - **a summary table for Larry**: one row per case, with columns
     `Case | Fixed / BLOCKED | Can the real pipeline reach it? (how) | Outputs changed (old → new)`;
   - the Part 1 trace;
   - `git diff <start> -- R/` verbatim with the hunk mapping, and the gate outputs;
   - updated and new expectations;
   - test, suite and check results beside the baseline, with measured wall clocks;
   - OPEN ITEMS.
2. Commit the report on the working branch (`dev: empty-set no-find report`).
3. **Gates:** `git diff --name-only <start>..HEAD` lists only
   `R/run_simulation_analysis.R`, `R/forestsearch_cross_validation.R`, `man/` files only
   if roxygen changed, test files the report names, `DESCRIPTION`, `NEWS.md`, this task
   document and the report. `git status --porcelain` is empty.
4. Copy the report to `~/Downloads/REPORT_empty_set_nofind_2026-11-11.md`.
5. Print the summary table, the gate results and the branch HEAD SHA. Do not push or merge.
