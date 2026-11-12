# TASK — CRAN preparation, part 2: the fixed seed, and a usage inventory of the legacy GBSG wrapper, 2026-11-11

Complete on its own. Run from the forestsearch repo root. Authorized by Larry 2026-11-11.
This task runs unattended and never stops to ask.

It continues on branch `cran-prep`, after `TASK_cran_prep_2026-11-11.md`.

## Purpose

- **Part 1 (changes behaviour, not the method).** `subgroup.consistency()` calls
  `set.seed(8316951)` and leaves the caller's random-number state changed. After this
  task, a call to an exported forestsearch function leaves the caller's state as it found
  it, and every value the package computes is unchanged.
- **Part 2 (read-only).** Which of the package's own vignettes and articles rely on the
  state that seed leaves behind.
- **Part 3 (read-only).** Where the legacy GBSG simulation wrapper is still used, so Larry
  can decide whether to remove it.

## Rules

- No `git fetch`, `pull`, `push` or `merge`. Every `git add` names explicit paths.
- No install into the user or site library.
- fs-post-selection and fs-glms-interpretable: read-only, no branch change.
- `vignettes/` and `quarto/`: read-only in this task.
- Evidence is the current source. Quote file:line.
- Gates are stop-on-failure for their own part, never stop-to-ask. If Part 1 fails a gate,
  revert it completely, report it as FAILED, and go on to Parts 2 and 3. The report is
  always written.
- No timing estimates. Report measured wall clocks.

## Step 0 — branch and record

1. `git status --porcelain` is empty. **Gate.**
2. If branch `cran-prep` exists, check it out. If
   `dev/reports/REPORT_cran_prep_2026-11-11.md` is not on it (or the branch does not
   exist), the first CRAN task has not been run: run it now, completely, with all its own
   gates, from `dev/tasks/TASK_cran_prep_2026-11-11.md` if present, else
   `~/Downloads/TASK_cran_prep_2026-11-11.md`. **Gate:** it ends with its report
   committed on `cran-prep`. If neither file exists, stop this task.
3. Copy `~/Downloads/TASK_cran_prep_seed_legacy_2026-11-11.md` to `dev/tasks/` and commit
   it on `cran-prep` (`dev: add CRAN preparation task, part 2`).
4. Record the HEAD SHA before step 3 (the start SHA) and the `DESCRIPTION` version.

## Part 1 — the caller's random-number state is restored

1. **Inventory.** List every `set.seed(` call in package code under `R/` (not roxygen
   examples): function, file:line, whether the seed is fixed or comes from an argument,
   whether the previous state is restored, and which exported functions reach it.
2. **Scope.** Change only the fixed-seed call in `subgroup.consistency()`
   (`R/subgroup_consistency_main.R` near `:498`) and what is needed around it. Seeds that
   come from a user-supplied argument are listed, not changed.
3. **Rule S.**
   - After any exported forestsearch function that reaches that call returns, the
     caller's `.Random.seed` is what it was on entry (absent again if it was absent).
   - **Every value the package computes is identical to the start SHA**, given the same
     inputs and the same caller state on entry. In particular, the draws made inside one
     exported call after the fixed seed is set must come from the same stream as now.
   - So the state is saved and restored at the outermost exported call only. A nested call
     must not restore part-way through an outer call. Use the smallest mechanism that
     achieves this.
4. **Equivalence check, in scratch** (`/tmp/seed_check_2026-11-11/`; start-SHA tree from
   `git archive`, branch tree; each loaded with `pkgload::load_all()` in its own
   `Rscript --vanilla` process; sequential plan; small synthetic data; hard cap 15 minutes
   in total). In both trees call `set.seed(123)` immediately before each call below, and
   compare the results with `identical()` after dropping timing fields:
   1. `forestsearch()`, survival, defaults;
   2. the same with `mr_inference = TRUE` and no seed in `mr_inference_args`;
   3. the same with a `df.test`;
   4. `forestsearch_bootstrap_dofuture()` on the fit of 1, small B;
   5. one K-fold cross-validation call, small K;
   6. `forestsearch()` on one GLM outcome.
   **Gate:** all six identical, and in the branch tree `.Random.seed` after each call
   equals `.Random.seed` before it.
5. **The consequence, recorded not gated.** In both trees run: `set.seed(123)`, then
   `forestsearch()`, then `forestsearch_bootstrap_dofuture()` with no seed between the
   two. Report whether the bootstrap output differs between the trees, and whether in each
   tree it repeats exactly when the three lines are run twice.
6. **Tests.** Add `tests/testthat/test-rng-state-restored.R`: the state after
   `forestsearch()` equals the state before it, with and without `mr_inference`; the same
   for a direct `subgroup.consistency()` call if it can be called directly; two
   consecutive identical calls give identical results. `skip_on_cran()` where a search
   runs; under 2 minutes.
7. `DESCRIPTION`: raise the development version by one. `NEWS.md`: one entry saying that
   forestsearch functions no longer leave the caller's random-number state changed; that
   their own results are unchanged; and that a later call that draws random numbers
   without its own seed (for example a bootstrap run straight after `forestsearch()`) now
   depends on the caller's seed, so set one.
8. **Suite and check**, once, on the check route of
   `dev/reports/REPORT_cran_prep_2026-11-11.md`; baseline = that report's final figures
   (do not re-run it). Hard caps: suite 30 minutes, `check()` 60 minutes.
   - Pre-set decision: a failing expectation that pins values of a call made after
     `forestsearch()` with no explicit seed is corrected by adding an explicit
     `set.seed()` before that later call and updating the pinned values. List each:
     file:line, old, new.
   - Any other failure, or a new check ERROR, WARNING or NOTE: revert Part 1, report it
     as FAILED with the failing output.
9. Commit by explicit path
   (`fix: forestsearch functions restore the caller's random-number state`).

## Part 2 — vignettes and articles that rely on the old state (read-only)

For every tracked file under `vignettes/` (the vignettes CRAN builds, and
`vignettes/articles/`): each place where a call that draws random numbers follows a call
reaching the fixed seed, with no explicit seed (a `set.seed()` or a seed argument) between
them.

| File:line | Later call | Seeded? | Displayed outputs that would vary from one render to the next | Built by CRAN, or pkgdown only |
|---|---|---|---|---|

No edit. State which rows are in a vignette that CRAN builds.

## Part 3 — the legacy GBSG simulation wrapper (read-only)

The wrapper is everything defined in `R/sim_aft_gbsg.R` and `R/setup_gbsg_dgm.R`,
including `create_gbsg_dgm()`, `setup_gbsg_dgm()`, `simulate_from_gbsg_dgm()`,
`calibrate_k_inter()`, `validate_k_inter_effect()`, `get_dgm_with_output()` and the
internal functions behind them.

1. **Per function:** file:line; exported or internal; present in the CRAN 0.4.0 release
   (from the release commit or tag; `CANNOT DETERMINE` if neither can be found);
   deprecation status.
2. **Callers, per function, with counts and one example each:**
   - other package functions under `R/` (name them; include `fs_fdr_report()` and
     `R/betaHhat_truth.R`);
   - `tests/`;
   - `vignettes/` and `vignettes/articles/` (name every file);
   - `quarto/`, split into current directories and archive / legacy / broken directories;
   - `inst/`, `_pkgdown.yml`, `README`;
   - fs-post-selection: `manuscript_rev1/` and its generators;
   - fs-glms-interpretable: `quarto/`, `manuscript/`, `dev/payloads/`.
3. **One summary table for Larry:**

   | Function | Exported | In CRAN 0.4.0 | Called by package code | Called by a vignette or article | Called by current manuscript work | Called only by archive or legacy documents |
   |---|---|---|---|---|---|---|

4. **What removal would involve:** the files, exports, help pages and tests that would
   go; each dependent package function that would have to change or go with it; and, for
   each current document that calls the wrapper, whether
   `generate_aft_dgm_flex()` with `simulate_from_dgm()` offers the same step (yes / no /
   unknown, from reading only). Describe; change nothing.

## Part 4 — report and close

1. Write `dev/reports/REPORT_cran_prep_seed_legacy_2026-11-11.md`, in bullets and tables:
   provenance (and whether Step 0 ran the first CRAN task); a status line per part;
   Part 1's inventory, diff, equivalence results, the recorded consequence, tests, suite
   and check beside the baseline with wall clocks; the Part 2 table; the Part 3 tables;
   OPEN ITEMS.
2. Commit it on `cran-prep` by explicit path (`dev: CRAN preparation report, part 2`).
3. **Gates:** `git status --porcelain` is empty; the two sibling repositories are
   unchanged; `git diff <start> --stat` lists only files the report names, and nothing
   under `vignettes/` or `quarto/`.
4. Copy the report to `~/Downloads/REPORT_cran_prep_seed_legacy_2026-11-11.md`.
5. Print the status lines, the Part 3 summary table, the gate results and the branch HEAD
   SHA. Do not push or merge.
