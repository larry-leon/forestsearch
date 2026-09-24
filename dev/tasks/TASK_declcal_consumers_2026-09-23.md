# TASK — Update the consumers of the changed declaration-calibration output

**Repo:** `larry-leon/forestsearch`, branch `feature/glm-extension` — run CC from the forestsearch clone.
Nothing is written to `fs-glms-interpretable`.

**Purpose.** The 2026-09-23 alignment (`96f84ad8`) removed `pstar_implied` from the c0 table, replacing it with
the settable-pair columns, and changed `fw_size` to the screen's rounded threshold. Scripts in this repo still
read the removed field and hold the old value.

**Why this blocks.** The failure is quiet, not loud. `tb$pstar_implied` returns `NULL` with no warning, and
`sprintf()` on a zero-length argument returns `character(0)` — so a line **disappears from the output** rather
than erroring. A report then prints one fewer row than it used to and nothing says why. Until this is fixed,
any declaration-calibration re-run produces output that cannot be trusted to be complete.

**`R/` CALLOUT — expected to be scripts only.** If a fix turns out to require an `R/` change, stop and report
rather than making it.

**Kind:** maintenance. **Compute: negligible.** Hard abort at 1 h.

---

## 0. First action

1. Copy this file to `dev/tasks/TASK_declcal_consumers_2026-09-23.md`; `git add` that path; commit
   `docs(tasks): update declaration-calibration consumers (2026-09-23)`.
2. Record HEAD, R version, platform.
3. `git status --short`: record pre-existing untracked files; never stage them.

---

## 1. Find every consumer — from source

Search the repo for references to the removed and changed fields: `pstar_implied`, `fw_size`, `z_pstar`,
`pconsistency_digits`, and any hard-coded `0.651`, `0.6508` or `0.9936`.

Known starting points, to be verified rather than assumed:

- `quarto/simulations/gbsg_020/scripts_dinamr/declcalc0_run.R` (reported at `:230`, `:315-321`, `:414`)
- `quarto/simulations/gbsg_020/scripts_dinamr/declcal_c0approx_run.R`
- `dev/verification/report_values*.R`

Report every hit with its file, line, and what it does with the value: reads it, formats it, stores it in a
payload schema, or compares against it.

---

## 2. Fix

For each consumer, replace the removed field with the settable-pair columns the c0 table now carries —
`pstar_settable`, `pstar_achievable`, `pcons_eff_settable`, `z_eff_settable`, `z_gap`, `digits_fine`,
`pstar_fine`, `z_gap_fine` — choosing whichever the consumer's purpose actually needs, and say why in the
report. Do not carry all eight into a display that needed one number.

Where a consumer reproduces the screen's rounded admission inline (`declcalc0_run.R:414` does), **call
`.fs_pcons_eff()` instead**. The 2026-09-23 task made "the effective threshold is computed in exactly one
place" a post-condition; a second copy in a script is the drift that produced this whole sequence.

Where a hard-coded old value appears, replace it with a value read from the payload, or update it and record
that it is now pinned to a specific commit.

---

## 3. The stale test, folded in

`test-declaration-calibration.R` test 1 still compares against the **unrounded** rule. It passes only because
its fixture happens to contain no candidate between 0.895 and 0.90 — so it would fail spuriously if the fixture
changed, and, worse, it would not catch a regression reverting the alignment.

Update it to the rounded rule, as the other three constants in that file and in `test-declaration-c0.R` were
updated. Report its old and new expectation. Confirm the updated test **can** fail on the old rule, as
`test-mr-admission-rounded.R` does.

---

## 4. Verification

- For each fixed consumer, demonstrate its read-out path produces **non-empty, non-`NA` output** on a small
  fit or a committed payload. Do not run a full campaign.
- **Gate 1 — no silent emptiness.** Assert that no field the consumer formats is zero-length. This is the
  failure mode the task exists to remove, so assert it rather than eyeball it.
- **Gate 2 — test suite.** `devtools::test()` passes; every updated test listed with old and new expectation.
- **Gate 3 — `R CMD check --as-cran`.** Baseline is **1 NOTE** ("Version contains large components
  (0.3.5.9000)"), deliberately open as the dev-version marker. Passes if the modified tree produces exactly
  that NOTE and nothing else. **Do not bump the version.** Skip this gate if no file under `R/` or `tests/`
  changed, and say so.
- **Gate 4 — no `R/` change.** If one proved necessary, stop and report instead.

---

## 5. Record and commits

Report to `dev/reports/REPORT_declcal_consumers_2026-09-23.md`.

Commits, explicit paths, in order: task doc; the script fixes; the test fix; the report.

---

## POST-CONDITIONS (machine-checkable)

1. Every hit from §1 is listed with file, line and disposition — fixed, or left with a stated reason.
2. No consumer references `pstar_implied`.
3. The effective-threshold expression still appears in exactly one place; scripts call `.fs_pcons_eff()`.
4. Gate 1: no formatted field is zero-length in any fixed consumer's demonstrated output.
5. `devtools::test()` passes; test 1's old and new expectations reported; the updated test shown to fail on
   the old rule.
6. No `R/` file modified.
7. Nothing written to `fs-glms-interpretable`.

---

## OUT OF SCOPE

No change to `R/`. No re-run of any declcal campaign — this task makes the consumers correct, it does not
regenerate their outputs. No manuscript edits. No change to the GBSG or ACTG 175 application documents, which
live in the other repo.
