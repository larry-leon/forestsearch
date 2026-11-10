# TASK — A-002 magnitude check (bootstrap id merge), 2026-11-10

Complete on its own. No change to `R/`, `tests/`, `man/`, `vignettes/`, `inst/`,
`quarto/`, `DESCRIPTION` or `NAMESPACE`. No package install into any library.

## Purpose

- Audit finding A-002 (`dev/reports/REPORT_code_audit_2026-11-09.md` on branch
  `fs-code-audit`): on the consistency path, `forestsearch()` attaches `treat.recommend` to
  `df.est` and `df.predict` by `merge()` on the id column. Inside a bootstrap replicate the
  resample repeats ids, so `df.est` gets sum(k^2) rows and `df.predict` holds in-bag subjects
  only. The replicate's "new subgroup on original data" estimate is then a bootstrap-sample
  estimate.
- This task measures, on one committed application, how much that changes the bias-corrected
  estimates. It compares the same method under two versions of the code: as is (arm A) and
  with flags attached row-wise (arm B). Same data, settings, seed and B.
- Decision rule, set by Larry: **material** if the H2 point estimate of H or of Hc moves by
  0.05 or more on the hazard-ratio scale; otherwise **immaterial**.

## Rules

- Gates are stop-on-failure, never stop-to-ask. On a failure, write what failed into the
  report and stop.
- Verify from the current source. Do not use earlier reports or descriptions as evidence of
  what the code does.
- Do not run the test suite, `check()`, `document()`, `build()` or any render.
- Run everything with `Rscript --vanilla` from a scratch directory under `tempdir()` (or
  `/tmp/a002_<date>`), with `future::plan(future::sequential)` so that both arms run the code
  loaded by `pkgload::load_all()` in the R session (multisession workers would load the
  installed package instead).
- Hard caps: arm A ≤ 30 min wall clock, arm B ≤ 30 min. On hitting a cap, stop the run,
  record the wall clock, and stop the task.
- Do not push.

## Step 0 — record and branch

1. On `master`, `git status --porcelain` must be empty. **Gate.**
2. Copy `~/Downloads/TASK_a002_magnitude_2026-11-10.md` to `dev/tasks/` and commit it on
   `master` (`dev: add A-002 magnitude task`).
3. **Gate:** branch `fs-a002-check` must not exist. Create it from `master` and check it out.
4. Record `master` SHA and `DESCRIPTION` version for the report.

## Step 1 — the application

1. Find the committed GBSG survival application that runs `forestsearch()` on the consistency
   path (default `subgroup_method`) followed by `forestsearch_bootstrap_dofuture()`. Look
   first for `quarto/applications/gbsg/analysis_gbsg_survival_effMaxSG.qmd`. If that file
   is absent or has no bootstrap call, use the one GBSG survival document under
   `quarto/applications/` that has both calls. **Gate:** zero, or more than one, candidate
   → stop and list what was found.
2. Copy into a scratch script, verbatim, everything the two calls need: data preparation, the
   `forestsearch()` call, the `forestsearch_bootstrap_dofuture()` call with its `nb_boots`
   and `seed`. Change only: `parallel_args` plans to `"sequential"`; nothing else. Do not
   edit the `.qmd`.
3. Record the document path, the `forestsearch()` arguments, `nb_boots` and `seed`.

## Step 2 — the patched copy (arm B code)

1. Export the tree to scratch: `git archive HEAD | tar -x -C <scratch>/fs_patched`. Nothing
   is written into the repo.
2. In `<scratch>/fs_patched/R/forestsearch_main.R`, change only the consistency-path block
   that builds the returned `df.est` and `df.predict` by `merge(..., by.x = id.name,
   by.y = "id", ...)` (around lines 3700–3735 at the audited SHA; locate it by reading):
   - `df.est`: attach `treat.recommend` from the flag table by row position, not by id. First
     confirm by reading `extract_subgroup()` (`subgroup_consistency_helpers.R`) that the flag
     table has one row per row of `df`, in the same order. **Gate:** if it does not, stop.
   - `df.predict`: evaluate the selected subgroup definition on `df.predict`'s own rows,
     transplanted from the way the DINA/GRF branches of the same file build their
     `df.predict` (the audit cites `get_dfpred` at `forestsearch_helpers.R:1664-1666`).
     Copy that code; do not write new logic.
   - Leave `df.test` and every other line unchanged.
3. Save `diff -u` of the original and patched file for the report.

## Step 3 — equivalence gate (no bootstrap yet)

1. Fit the application's `forestsearch()` call once with the repo tree loaded (arm A code) and
   once with `<scratch>/fs_patched` loaded (arm B code), each in its own `Rscript` process.
2. **Gate:** same `sg.harm`; same `treat.recommend` for every subject in `df.est` and in
   `df.predict` after sorting by id; same row counts. On original data the ids are unique,
   so the two codes must agree here. Any difference → stop.

## Step 4 — the two bootstrap runs

1. Arm A: repo tree loaded, the application's `forestsearch()` then
   `forestsearch_bootstrap_dofuture()` exactly as recorded in Step 1. Save the bootstrap
   object to scratch (`.rds`).
2. Arm B: the same script with `<scratch>/fs_patched` loaded.
3. Record each arm's wall clock.
4. **Gate:** per replicate, the selected subgroup is identical in the two arms
   (`results$M.1`, `results$M.2` and the found / not-found indicator, compared row by row).
   Selection happens before the changed code, so any difference means the comparison is not
   clean → stop and report how many replicates differ.

## Step 5 — report

Write `dev/reports/REPORT_a002_magnitude_2026-11-10.md`:

1. Header: `master` SHA, version, application document, `forestsearch()` arguments,
   `nb_boots`, `seed`, wall clock per arm.
2. Gate checklist, pass/fail.
3. The patch: the `diff -u` from Step 2.
4. Results, one table, HR scale:

   | Subgroup | Quantity | Arm A (as is) | Arm B (row-wise) | Difference (B − A) |
   |---|---|---|---|---|
   | H | naive (H0) |  |  |  |
   | H | H2 estimate |  |  |  |
   | H | H2 95% CI |  |  |  |
   | Hc | naive (H0) |  |  |  |
   | Hc | H2 estimate |  |  |  |
   | Hc | H2 95% CI |  |  |  |

   Also report the number of replicates with a subgroup found (same in both arms by the
   Step 4 gate). H1 is not reported.
5. Verdict, one line: `MATERIAL` if |Difference| ≥ 0.05 for the H2 estimate of H or of Hc,
   otherwise `IMMATERIAL`, stating both differences.
6. No recommendation, no proposed fix.

## Step 6 — close

1. Commit the report on `fs-a002-check` (`dev: A-002 magnitude report`).
2. **Gate:** `git diff --name-status master..fs-a002-check` lists exactly one line, `A`, the
   report. `git status --porcelain` empty.
3. Do not merge, push or delete the branch. Leave scratch files where they are.
4. Copy `dev/reports/REPORT_a002_magnitude_2026-11-10.md` to `~/Downloads`.
5. Print the results table, the verdict and the gate checklist.
