# TASK — NPV for a replicate with no subgroup: all remaining instances, 2026-11-11

Complete on its own. Run from the forestsearch repo root. Authorized by Larry 2026-11-11:
"all instances need to be fixed".

## Purpose

- **Kind of change: changes behaviour. It does not change the method.**
- The definition (León et al. 2024, Section 3): when no subgroup is identified, the
  estimated complement is everyone. So for that replicate:
  - `sens` = 0 and `ppv` = 0 (convention for 0/0);
  - `spec` = sens(Hc) = 1;
  - `npv` = ppv(Hc) = |Hc|/N, the share of that replicate's subjects truly in the
    complement. Under the null this is 1.
- The A-003 fix (`dev/reports/REPORT_a003_fix_2026-11-11.md`) made the extractors in
  `R/run_simulation_analysis.R` follow this. Other places still score `npv` = 1 for such a
  replicate under an alternative. That report's open items 1 and 2 name two of them.
- This task finds every remaining instance (Part A) and fixes the ones in the package
  (Part B). Instances outside the package are listed, not changed.

## Rules

- No `git fetch`, `pull`, `push` or `merge`. Every `git add` names explicit paths.
- No package install into any library. No render, no simulation or application re-run.
- fs-post-selection and fs-glms-interpretable: read-only, no branch change.
- `quarto/` and `vignettes/` in forestsearch: read-only in this task.
- Evidence is the current source. Quote file:line.
- Gates are stop-on-failure, never stop-to-ask. On a failure, write what failed into the
  report, commit the report, and stop.

## Step 0 — record and branch

1. On `master`, `git status --porcelain` is empty. **Gate.**
2. **Gate:** the A-003 fix is in `master`:
   `grep -cE 'npv +<- +sum\(!true_H\) / nrow\(df\)' R/run_simulation_analysis.R` returns 4.
   If not, stop: branch `fix/a003-nofind-spec-npv` has to be merged first.
3. Copy `~/Downloads/TASK_npv_nofind_instances_2026-11-11.md` to `dev/tasks/` and commit it
   on `master` (`dev: add NPV no-find instances task`).
4. Record `master` SHA and `DESCRIPTION` version (expected 0.4.1.9004).
5. **Gate:** branch `fix/npv-nofind-instances` does not exist. Create it from `master`.

## Part A — inventory (read-only)

Search tracked files in forestsearch `R/`, `tests/`, `man/`, `vignettes/` and `quarto/`
for every site that sets, fills, documents or tests the value of NPV for a replicate with
no identified subgroup (also written as not detected, not declared, `any.H == 0`,
`detected == 0`). Names to search include `npv`, `NPV`, `ppv_sgc`, `ppv(Hc)`, `ppv_Hc`.
At each such site also read the value given to specificity; it should be 1.

Known starting points (verify each; they are not the full list):

- `R/fs_mr_oc_summary.R`: roxygen near `:56`; `.fs_oc_identification` fill near
  `:280-309`; `man/fs_mr_oc_summary.Rd`.
- `R/simulation_tables.R`: `build_classification_table()` fill and comment near
  `:360-372`.
- `tests/testthat/test-fs_mr_oc_summary.R:108`.
- `quarto/simulations/actg175/binary_methods/actg175_binary_methods_template.qmd:357-360`
  (`oc_agg_from_results`) and its copies.
- The `.classify()` helpers in the multiplier-resampling simulation templates (e.g.
  `quarto/simulations/actg175/binary_020/sim_fs_mr_field_or_template.qmd:533-543`).

Then, read-only, the same search in:

- fs-glms-interpretable (`~/Documents/GitHub/fs-glms-interpretable`): start at
  `dev/payloads/R/identifiers_metrics.R:359-362`; include `manuscript/` and `quarto/`.
- fs-post-selection (`~/Documents/GitHub/fs-post-selection`): `manuscript_rev1/` and its
  generators.

One table, one row per site (copies of one template may share a row, with a count):

| Repo | File:line | Kind (code / roxygen / comment / test / quarto source) | Value given to NPV and to spec for a no-subgroup replicate | Can the NPV value be reached under an alternative? (yes / no, and why) | Is \|Hc\|/N for that replicate available at the site? (from what) | Reported outputs it feeds (table, figure, payload) |
|---|---|---|---|---|---|---|

Also state: how many tracked payloads store unconditional NPV computed with the old
value, and whether any tracked document displays it.

## Part B — fix (forestsearch `R/`, `man/`, `tests/`, `DESCRIPTION`, `NEWS.md` only)

Apply to every Part A row whose file is under `R/`:

- **R1 — the value.** For a replicate with no identified subgroup, NPV is |Hc|/N for that
  replicate. Specificity is 1. `sens` and `ppv` are not changed.
- **R2 — where |Hc|/N comes from.** Use quantities the site already has: a per-replicate
  true-subgroup size or flag, sizes in the same row, or a stored value already computed by
  the found-case formula. Do not add arguments, columns or inputs.
- **R3 — not available.** If the site cannot obtain |Hc|/N from what it already has, leave
  its code unchanged and list it in the report as `BLOCKED`, with exactly what it lacks.
  Do not substitute another value and do not leave the fill out.
- Nothing else changes: no refactor, no renaming, no change to detection flags, `sens`,
  `ppv` or any conditional (detected-only) summary.
- **Roxygen and comments.** Where text states the old convention (NPV = 1 for a
  non-detected replicate, or "matches `build_classification_table()` / Table 1" where that
  is no longer true), correct it to the definition above. If a site is `BLOCKED`, its text
  must describe what the code does, not what it should do. If any roxygen changed, run
  `devtools::document()` once.
- `DESCRIPTION`: 0.4.1.9004 → 0.4.1.9005.
- `NEWS.md`, development section, one entry naming each function changed and stating:
  a replicate with no identified subgroup now contributes NPV = |Hc|/N (it was 1);
  specificity, sensitivity, PPV and detected-only summaries are unchanged; stored payloads
  are not regenerated.

**Gates after the edits:**

1. `git diff master --stat -- R/` lists only files that have a Part A row.
2. Every hunk in `git diff master -- R/` maps to a Part A row; the report gives the
   mapping.
3. If `document()` ran: `git diff master --stat -- man/` lists only the Rd files of
   functions whose roxygen changed.

## Part C — tests

- Update `tests/testthat/test-fs_mr_oc_summary.R:108` (and any other expectation that pins
  NPV = 1 for a non-detected replicate under an alternative) to the new value. List each:
  file:line, old value, new value.
- Add, in the test file of each function changed, for a hand-built input with at least one
  detected and one non-detected replicate under an alternative:
  1. the non-detected replicate contributes NPV = |Hc|/N and spec = 1;
  2. the same input under the null gives NPV = 1 and spec = 1.
- Add one agreement expectation if it can be built by hand without a search: on the same
  replicates, the unconditional NPV from `fs_mr_oc_summary()` equals the all-replicate mean
  of `npv` as the `R/run_simulation_analysis.R` extractors score those replicates. If it
  cannot be built by hand, say so in the report.
- No search, no simulation, no parallel plan in any test. Hard cap: each test file runs in
  under 2 minutes.

## Part D — tests and check

1. Run the changed and new test files. **Gate:** 0 failures.
2. Then, once: the full suite and `devtools::check()`, on the same check route as
   `dev/reports/REPORT_a003_fix_2026-11-11.md` §4.
   - Baseline, from that report (do not re-run it): suite 600 blocks, 8,183 expectations,
     0 failures, 3 skips, 38 warnings; check 0 errors / 0 warnings / 0 notes.
   - Hard caps: suite 30 minutes, `check()` 60 minutes. On a cap, stop the run, record the
     wall clock, and stop the task.
   - Pre-set decision: a failing expectation that pins NPV = 1 for a non-detected replicate
     under an alternative is updated and listed. **Gate:** any other failure, or any new
     check ERROR, WARNING or NOTE → stop.
3. Commit on `fix/npv-nofind-instances` by explicit path
   (`fix: NPV for a replicate with no subgroup is |Hc|/N in all package summaries`).

## Part E — report and close

1. Write `dev/reports/REPORT_npv_nofind_instances_2026-11-11.md`, in bullets and tables:
   - provenance;
   - the Part A table and the stored-payload statement;
   - **a summary table for Larry**: one row per instance, with columns
     `Instance | Fixed here / BLOCKED / outside the package | Reported output affected | What addressing it needs (source edit, re-render, re-run)`;
   - `git diff master -- R/` verbatim with the hunk-to-row mapping, and the gate outputs;
   - updated and new expectations;
   - test, suite and check results beside the baseline, with measured wall clocks;
   - OPEN ITEMS.
2. Commit the report on the branch (`dev: NPV no-find instances report`).
3. **Gates:** `git diff --name-only master...HEAD` lists only files under `R/`, `man/` and
   `tests/testthat/` that the report names, plus `DESCRIPTION`, `NEWS.md` and the report.
   `git status --porcelain` is empty here and unchanged in the two sibling repositories.
4. Copy the report to `~/Downloads/REPORT_npv_nofind_instances_2026-11-11.md`.
5. Print the summary table, the gate results and the branch HEAD SHA. Do not push or merge.
