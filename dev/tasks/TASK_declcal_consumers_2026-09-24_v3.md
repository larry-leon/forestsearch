# TASK v3 — Update the consumers of the changed declaration-calibration output

**Supersedes** v2 (`TASK_declcal_consumers_2026-09-24_v2.md`, never run, never committed) and v1
(`dev/tasks/TASK_declcal_consumers_2026-09-23.md`, `1b11967c`, which stopped at its inventory by rule).
v3 is v2 with four corrections: the standing check policy, Gate 2's residual classes, a test fixture that
can actually discriminate, and no catalogue writes in `gbsg_020`.

**Repo:** `larry-leon/forestsearch`, branch `feature/glm-extension`, the Pop!_OS clone. Nothing is written to
`fs-glms-interpretable`.

**Purpose.** The 2026-09-23 alignment (`96f84ad8`) removed `pstar_implied` from the c0 table and moved `fw_size`
to the screen's rounded threshold (`7713942e` aligned MR admission the same way). Scripts in this repo still read
the removed field, hold old values, or compute their own version of the threshold. The applications re-run in
`fs-glms-interpretable` is blocked behind this.

**`R/` CALLOUT — scripts and tests only. No `R/` change, not even temporarily.** The OC predictor's exact-cutoff
gate (`R/fs_oc_predict.R:279`, `R/fs_oc_grid.R:582`) is out of scope. If any other fix turns out to need `R/`,
stop and report.

**Kind:** maintenance. **Compute:** negligible. No campaign and no simulation. **Hard abort at 2 h** for the whole
task.

**Git:** never `fetch`, `pull` or `push`. Every `git add` names explicit paths. Pre-existing untracked files are
never staged.

---

## 0. First action

1. Copy this file to `dev/tasks/TASK_declcal_consumers_2026-09-24_v3.md`, `git add` that path, and commit
   `docs(tasks): declcal consumers v3 (2026-09-24)`.
2. Record HEAD, R version, platform and the installed forestsearch build time. Assert
   `exists(".fs_decl_settable", envir = asNamespace("forestsearch"), inherits = FALSE)` and stop if it is false.
3. Run `git status --short`, record the pre-existing untracked files, and never stage them.

---

## 1. Take the inventory from the record — do not redo the search

Read `dev/reports/REPORT_declcal_consumers_2026-09-23.md` (`4ef26719`) and work from it. Verify each line number
against current source as you edit, since the report is a day old.

---

## 2. The three dispositions (decided)

**2a. Committed payloads storing the exact-cutoff value.** For each row, derive the settable pair from the stored
`kappa_hat` with `.fs_decl_settable()` (`R/fs_declaration_calibration.R:205`), reached the way the script already
loads forestsearch. Keep no second inverse anywhere. Where a script stores the column for future runs, store the
derived settable values and drop the exact-cutoff column. Do not re-run any campaign.

**2b. Run scripts are fixed; archived findings are annotated.** Classify every inventory site, per file:

- **Run scripts** are anything invoked to produce new output: campaign drivers, verification scripts and
  read-out scripts. **Fix them.**
- **Archived findings** are committed outputs and theory notes. **Do not rewrite them.** Add one header line
  saying the values were computed at the exact cutoff, before the 2026-09-23 alignment, and naming this task.
  Change nothing else.
- A file that is both counts as a run script and is fixed. Its committed output gets the header line.

**2c. The live bug, fixed regardless.** `dev/verification/report_values_c0.R:12` calls
`f(t05$pstar_implied[i])`, which prints nothing. Treat `report_values_c0.R:13, 22` and
`report_values.R:12-13, 15, 29-38, 49-50` the same way. `report_values.R:50` rebuilds the admission rule inline
from `z_pstar`; replace it with the package's threshold, not a local reconstruction.

---

## 3. The stale test

`tests/testthat/test-declaration-calibration.R` test 1 (`:111-124`) rebuilds the screen without rounding:
`rate <- 2 * pnorm(T) - 1; rate >= dc$p_star`.

- **Oracle.** State the screen's rounded rule literally, `round(rate, digits) >= dc$p_star`, at the fixture's
  `pconsistency.digits` (default 2). It must be independent of `.fs_pcons_eff()`, so that a regression inside
  the helper is caught.
- **Fixture.** The fixture must contain at least one candidate on which the rounded and unrounded rules
  disagree. Without one, the test cannot tell the two rules apart. If the fixture has no such candidate,
  adjust it and record the change. Add an assertion that such a candidate exists, with a message saying the
  test loses its power without it.
- **Old-rule failure.** In a scratch copy **outside the repo**, swap the oracle back to the unrounded rule and
  run it. It must fail. Put the failure output in the report.
- Record test 1's old and new expectations.

---

## 4. Verification — every gate stops on failure

- For each fixed run script, run its read-out path on a committed payload. Use a small single fit only where
  no payload path exists. Run no campaign.
- **Gate 1 — no silent emptiness and no wrong number.** For each fixed consumer, no formatted field is
  zero-length. Every implied-p\* or threshold value it reports must equal what the calibration returns for the
  same input.
- **Gate 2 — no second derivation left.** Re-run the inventory's own search patterns over the tree outside
  `R/`. List every remaining hit in the report, assigning each to exactly one class:
  - (a) an annotated archived file, with the header line present;
  - (b) a test in `tests/testthat/` that states the screen's rule as a deliberate oracle, including §3's
    band assertion;
  - (c) `tests/testthat/test-fs-oc-predict.R`, which pins the out-of-scope OC gate and is untouched pending
    Larry's decision.

  The gate fails on any hit in a run script, and on any hit that fits none of (a)–(c).
- **Gate 3 — this task's own test file only.** Run `devtools::test(filter = "declaration-calibration")` with a
  hard cap of 10 min. Do not run the full suite, `R CMD check` or a vignette build; those run only when Larry
  asks.
- **Gate 4 — no `R/` change.** `git status --short -- R/` must be empty at the end.

---

## 5. Report, catalogue and commits

Write the report to `dev/reports/REPORT_declcal_consumers_2026-09-24_v3.md`:

1. **One table with one row per inventory site:** file, line, defect kind, class (run script / archived / out of
   scope), what changed (or "header added"), and the Gate 1 result.
2. **Gate 2:** the residual hits by class, then one sentence of the form "Outside `R/`, no run script derives the
   threshold; remaining hits: n(a) archived, n(b) test oracles, n(c) OC pin."
3. **Test 1:** old and new expectations, the fixture change (if any) and the old-rule failure output.

**Catalogue:** do **not** edit `quarto/simulations/gbsg_020/status_curated.md` and do not regenerate
`quarto/simulations/gbsg_020/current_status.md`. The Section 5 re-run on the Mac Studio commits both files at its
closeout, and a second edit here would conflict when the two clones are reconciled. Instead, list in the report
the catalogue entries this task owes (the fixed drivers and the annotated archived files), so the next
`gbsg_020` closeout picks them up. Any other simulation directory touched follows its own closeout rule as usual.

Commit with explicit paths, in this order: task doc; run-script fixes; archived-file header lines; test fix;
report.

---

## POST-CONDITIONS (machine-checkable)

1. Every inventory site is dispositioned: fixed, annotated, or out of scope with the reason.
2. No consumer references `pstar_implied`.
3. Gate 1 holds for every fixed consumer.
4. Gate 2: every hit outside `R/` falls in class (a), (b) or (c), and none is in a run script.
5. Gate 3 passes. Test 1's fixture contains a candidate on which the two rules disagree, and the old-rule
   failure is shown.
6. No `R/` file is modified.
7. `gbsg_020/status_curated.md` and `gbsg_020/current_status.md` are unmodified.
8. Nothing is written to `fs-glms-interpretable`.

---

## OUT OF SCOPE

- The OC predictor's exact-cutoff gate (`R/fs_oc_predict.R:279`, `R/fs_oc_grid.R:582`) and the test that pins
  it (`test-fs-oc-predict.R:214`). This is a separate decision.
- Any `R/` change.
- Any re-run of a declcal campaign.
- The GBSG and ACTG 175 application documents, which live in the other repo.
- Manuscript edits.
- The full test suite and CRAN checks.
