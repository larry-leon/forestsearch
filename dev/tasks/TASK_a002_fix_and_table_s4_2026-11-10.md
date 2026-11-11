# TASK — A-002 fix (arm B) on a forestsearch branch, then Table S4 under arm B, 2026-11-10

Two parts. Part B starts only if Part A ends green. Gates are stop-on-failure, never
stop-to-ask. Do not push or merge in either repository; Larry does that.

- **Part A** (forestsearch): implement arm B on a new branch. **Changes behaviour**
  (bootstrap bias-corrected estimates on the consistency path); does not change the method —
  it makes the code compute León et al. (2024) Eq. 7 as written.
- **Part B** (fs-post-selection): re-run the simulation behind supplement Table S4 with the
  Part A branch, and compare against the committed Table S4. Nothing committed is
  overwritten.

Background: `dev/reports/REPORT_a002_magnitude_v6_2026-11-10.md` on forestsearch branch
`fs-a002-check` (arm A = merge on id, arm B = row-wise; v1 report §3 holds the two-line diff).

---

## Part A — forestsearch branch

### A0. Pre-conditions

1. forestsearch repo on `master`, `git status --porcelain` empty (**gate**).
2. Copy this file to `dev/tasks/` and commit it on `master`
   (`dev: add A-002 fix and Table S4 task`).
3. **Gate:** branch `fix/a002-bootstrap-rowwise` does not exist. Create it from `master`.

### A1. The change (R/forestsearch_main.R only)

Apply exactly the two-line change recorded in
`dev/reports/REPORT_a002_magnitude_2026-11-10.md` §3 (on branch `fs-a002-check`; read it with
`git show fs-a002-check:dev/reports/REPORT_a002_magnitude_2026-11-10.md`), at the
consistency-path merge block (`df.est_out` and `df.predict_out`):

```r
-      df.est_out <- merge(df, temp, by.x = id.name, by.y = "id", all.x = TRUE)
+      df.est_out <- df
+      df.est_out$treat.recommend <- temp$treat.recommend
...
-        df.predict_out <- merge(df.predict, temp, by.x = id.name, by.y = "id", all.x = TRUE)
+        df.predict_out <- get_dfpred(df.predict, sg.harm, version = 2)
```

- **Gate:** after the edit, `diff` of `R/forestsearch_main.R` against
  `/tmp/a002_2026-11-10/fs_patched/R/forestsearch_main.R` is empty (if that file no longer
  exists, the diff against `master` must be exactly these two hunks).
- Add a three-line code comment above the `df.est_out` line: flags are attached by row
  position (the flag table is built row-wise from `df`) and `df.predict` is evaluated by the
  rule, because a bootstrap resample repeats ids and a merge on id would compute Eq. 7's
  observed-data term on the resample.
- `DESCRIPTION`: version 0.4.1.9002 → 0.4.1.9003.
- `NEWS.md`, development section, one entry: bootstrap bias correction on the consistency
  path now evaluates the bootstrap-selected subgroup on the original data and on the
  bootstrap sample as Eq. 7 defines (previously both were computed on the resample through a
  merge on repeated ids); bias-corrected H/Hc estimates from
  `forestsearch_bootstrap_dofuture()` after a consistency-path fit change; naive estimates,
  selection, DINA/GRF and multiplier resampling are unchanged.
- Nothing else under `R/` changes. No roxygen change, so no `document()`.

### A2. Regression test

New file `tests/testthat/test-bootstrap-rowwise-flags.R`, transplanted from the A-002
reproducer R01 (`dev/reports/REPORT_code_audit_2026-11-09.md` appendix, branch
`fs-code-audit`): `set.seed(1)`, n = 600, the same data, `maxk = 1`, sequential plan.
Expectations:

1. `nrow(fs$df.est) == nrow(boot)` (resample size, not sum(k^2)).
2. `nrow(fs$df.predict) == nrow(d)` and no `NA` in `fs$df.predict$treat.recommend`.
3. `fs$df.predict$treat.recommend` equals `get_dfpred(d, fs$sg.harm, version = 2)$treat.recommend`.
4. On the original data (unique ids, `df.predict` NULL), `df.est$treat.recommend` sorted by
   id equals the value the merge gave before the change (compute the merge inline in the
   test for comparison).

Add `skip_on_cran()`.

### A3. Tests and check

1. Run the new test file and every test file whose name contains `bootstrap` (task's own
   tests). **Gate:** 0 failures.
2. Then, once, the full suite and `devtools::check()` (branch touches `R/`). Baseline on
   record: check 0/0/0; suite 0 failures, 3 skips (REPORT_wrapper_batch_2026-11-06).
   - Pre-set decision: a failure whose expectation pins a bootstrap bias-corrected value
     (H/Hc `H1`/`H2`, `H_biasadj_*`, `Hc_biasadj_*`) after a consistency-path fit is updated
     to the new output; list each such file:line, old value and new value in the report.
   - **Gate:** any other failure, or any new check ERROR/WARNING/NOTE → stop.
3. Commit on `fix/a002-bootstrap-rowwise`
   (`fix: bootstrap Eq. 7 terms use row-wise subgroup flags (A-002)`).
4. Record the branch HEAD SHA for Part B.

---

## Part B — Table S4 under arm B (fs-post-selection)

### B0. Locate

1. fs-post-selection repo: `git status --porcelain` empty (**gate**). Record branch and SHA.
2. Find the committed payload behind Table S4:
   `_payloads/fs_maxeffCons_fb_mr_m1_h10_knoise0_n500_combined_1_500.rds`, read by
   `_sim_estimation_fs_h10.qmd`. Record its path and its `meta` (n_sims, n_sample, nb_boots,
   mr_draws, seeds, any recorded package version / SHA / wall clock).
3. Find the document or script that generates that payload (grep the repo for
   `fs_maxeffCons_fb_mr` and for `combined_`). **Gate:** exactly one generator; otherwise stop
   and list candidates.
4. Create branch `a002-table-s4` from the current branch. Copy this task file into the
   repo's task folder (whichever folder holds this repo's task docs; if none, `dev/tasks/`)
   and commit it.

### B1. Run

1. Copy the generator to scratch `/tmp/a002_tableS4/` as a script. Transplant only; marked
   changes:
   - `library(forestsearch)` → `pkgload::load_all("<forestsearch repo path>")` with the
     forestsearch checkout on `fix/a002-bootstrap-rowwise` at the Part A SHA. Leave that
     checkout untouched for the duration.
   - Outer parallel plan → `future::plan(future::multicore, workers = 100)` (forked workers
     run the loaded branch code; multisession workers would load the installed package).
     Any inner plan the generator sets per simulation stays as it is.
   - Output payload name: the committed stem with `_armB` before `_combined`, written to
     scratch.
   - Nothing else: same DGM, n_sims, n_sample, nb_boots, mr_draws, seeds, chunking.
2. **Fork gate** before the run: a multicore future returns TRUE for
   `grepl("get_dfpred(df.predict, sg.harm", paste(deparse(forestsearch::forestsearch), collapse = ""), fixed = TRUE)`.
3. Run with `Rscript --vanilla` from a terminal. **Hard cap: 10 h wall clock.** On the cap,
   stop and report the wall clock and replicates completed.

### B2. Compare and render

1. Render a copy of `_sim_estimation_fs_h10.qmd` as `_sim_estimation_fs_h10_armB.qmd`
   (only change: `rds_candidates` → the `_armB` payload) to HTML. Do not edit the original
   fragment.
2. Comparison table (report), every cell of Table S4: committed vs arm B vs difference, for
   both blocks (H, Hc) and all four estimators, plus the caption quantities (detected
   replicates, detection rate, mean |Ĥ|, sens/spec/PPV/NPV, κ).
3. Record, without stopping: whether the non-FB quantities (detection, classification,
   oracle, naive, MR rows) are identical to the committed payload. Only the FB rows should
   move. If non-FB quantities differ, say so plainly: the committed payload was produced
   under a different package state, and the FB difference is then not attributable to A-002
   alone.

### B3. Record and close

1. Write `REPORT_a002_table_s4_2026-11-10.md` in the repo's reports folder (whichever
   folder holds this repo's reports; if none, `dev/reports/`): Part A summary (SHA,
   tests, check, any updated pinned values), Part B header (generator, settings, wall
   clock), gate checklist, the comparison table, and one line on whether only FB moved. No
   recommendation.
2. Commit on `a002-table-s4`: the report, `_sim_estimation_fs_h10_armB.qmd`, and the
   `_armB` payload in `_payloads/` if it is under 50 MB (otherwise keep it in scratch and
   record its path and size). Never modify the committed payload or fragment.
   **Gate:** `git diff --name-status` against the starting branch lists only added files.
3. Copy to `~/Downloads`: the report and the rendered `_sim_estimation_fs_h10_armB.html`.
4. Print the comparison table and the gate checklist.
