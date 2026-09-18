# CC TASK — Directive B: the binary default estimand becomes OR

**Opened:** 2026-09-18
**Repository:** forestsearch (the package).
**Authorized by:** Larry, 2026-09-18.
**Evidence base:** `quarto/simulations/actg175/binary_020/REPORT_criterion_and_defaults_audit_2026-09-18.md`
(`ce5f1e84`), §3.1, §3.8, Part C. Every cited line is a pointer — **re-verify from current source.**
**Compute:** about six single fits, seeded, bootstrap and CV off. Target under two minutes; **hard abort at
ten** — if exceeded, stop and report. Nothing else runs.

## Classification

| Change | Class |
|---|---|
| Binary default `effect_measure`: `"RD"` → `"OR"` at the live resolution site | **changes behaviour** (a default) |
| Deleting the second, unreachable resolution site | removes existing code (after proof, step 3) |
| Fixing the wrong scale claim in `subgroup_search.R`'s `hr.threshold` roxygen | documentation only |

**No method change.** No selection logic is touched. Explicit `effect_measure` of any value keeps exactly its
current behaviour. Continuous `MD` and count `IRR` defaults are unchanged. `adverse_outcome` is not touched
at either of its sites. No threshold validation, no derivation, no c2/c1 rule — that is Directive A, not this
task.

## Rules

- Copy this document into `dev/tasks/` and commit it before anything else.
- **No install, reinstall, or change to any R library.** Work via `devtools::load_all()`;
  `rcmdcheck::rcmdcheck(args = "--as-cran")` is permitted (temporary library).
- **Do not push.**
- Gates are stop-on-failure, not stop-to-ask.
- Find everything by search, never by the audit's line numbers.

## Steps

### Step 0 — task document into the record

Commit this file to `dev/tasks/` alone.

### Step 1 — source verification on a clean tracked tree

1. Locate every site where `effect_measure` is resolved from `outcome_type`. The audit found exactly two in
   `forestsearch_main.R` (`:1334-1341` live; `:1731-1738` inside the estimator-closure block) and no third
   anywhere in `R/`. Confirm the count by search over all of `R/`.
2. **Prove or refute the audit's unreachability claim for the second site:** the first site runs
   unconditionally for `outcome_type != "survival"` and leaves `effect_measure` non-`NULL`, so the second
   site's `is.null(effect_measure)` condition can never hold on any path through `forestsearch()`. Trace it
   from source. Also check whether any other function calls into the code containing the second site with
   `effect_measure` still `NULL`.
3. Record where `args_call_all` is captured relative to the first resolution site, and therefore whether a
   bootstrap/CV replay carries the resolved estimand or re-resolves it.
4. Locate the campaign template's `forestsearch()` call (the audit's Part C: `sim_fs_mr_field_or_template.qmd`,
   `effect_measure = "OR"` at `:307`/`:657`) — the real file for check 6.

**Gate 1:** if the site count is not two, or the reachability analysis contradicts the audit, stop and report.

### Step 2 — baseline capture, before any edit

On the clean tracked tree (untracked campaign artefacts are acceptable; record the SHA), run a resolution
probe — no fits — recording resolved `effect_measure`, resolved screening and consistency thresholds, and
their comparison scale, for every cell:

| Cell |
|---|
| binary, `effect_measure` unset |
| binary, `= "OR"` explicit |
| binary, `= "RD"` explicit |
| binary, `= "RR"` explicit |
| each of the four, under `subgroup_method = "consistency"`, `"dina"`, `"grf"` |
| survival, unset |
| continuous, unset and `= "MD"` explicit |
| count, unset and `= "IRR"` explicit |

Then the fit-level baseline: seeded, bootstrap off, on a fixed-seed simulated binary dataset defined in the
test file (and one survival, one MD dataset). Fits: binary unset; binary `"OR"` explicit; binary `"RD"`
explicit; survival default; MD default. Record selection digests.

Commit probe results + digests with the SHA. **No source edit yet.**

### Step 3 — the change

1. At the live site: `binary = "RD"` becomes `binary = "OR"`. Nothing else in the `switch` changes.
2. If step 1 **proved** the second site unreachable: delete that duplicate block, in its **own commit**,
   classified removes-code, citing the proof in the commit message. If not proven: change it identically to
   the live site **in the same commit as** the default change, and record that consolidation is deferred.
3. The roxygen for `effect_measure` (wherever it documents the per-outcome defaults) is updated to state the
   binary default is `OR`, and the docs-task threshold-vocabulary table's RD line is checked: `RD` rows
   describe **explicit** `RD` now — adjust wording only where it asserts RD is the binary default.

### Step 4 — the documentation rider (own commit, doc-only)

Fix `R/subgroup_search.R` `hr.threshold` roxygen (the audit's finding: `:25` claims log scale "for ratio
measures (OR, HR)"). Correct statement, verified from source first: on the survival path the comparison is
against the fitted HR on the **natural** scale (`effect_threshold` is `NULL` there, so the natural
`hr.threshold` reaches the comparison); on GLM paths the passed value arrives already on the comparison
(link or identity) scale. The rider obeys the docs task's rule: only roxygen lines change in this commit.

### Step 5 — re-run the probe and diff

**Gate 5:** the only cells that differ from baseline are binary-with-`effect_measure`-unset (all three
identifiers): resolved measure RD → OR, thresholds 0.05/0.0 → log(1.25)/log(1.0), scale identity → log.
Every other cell — every explicit cell, survival, MD, IRR — is **identical**. Any other difference fails the
task.

### Step 6 — fits

Re-run the step 2 fits. **Gate 6:** binary `"OR"` explicit, binary `"RD"` explicit, survival, MD — digests
**identical** to baseline. Binary unset — completes on the OR path; record its (expectedly different)
selection beside the baseline's, both kept in the report.

### Step 7 — tests, NEWS, check

testthat tests for the acceptance checks below (failing pre-change where applicable — state that this was
confirmed). `devtools::document()`. `NEWS.md`: one behaviour-change entry naming who is affected — binary
callers who inherit `effect_measure` — and that explicit values are untouched. Then
`rcmdcheck --as-cran`: **no new** finding relative to a pre-change run at this task's starting SHA (the
docs-task tree; record both sets).

### Step 8 — report

`REPORT_binary_default_or_2026-09-18.md` beside the existing reports in
`quarto/simulations/actg175/binary_020/`: step 1's proof, the SHA-pinned probe diff, fit digests with wall
clock, each acceptance check's result, commit list with classifications, findings with no task attached.
Commit. **Do not push.**

## Acceptance checks

1. `outcome_type = "binary"` with `effect_measure` unset resolves to `"OR"`, with ratio-scale thresholds
   `log(1.25)` / `log(1.0)` — not `0.05` / `0.0` — under all three `subgroup_method` values.
2. Binary unset and binary `"OR"` explicit resolve identically.
3. Binary `"RD"` explicit resolves exactly as at baseline (`0.05` / `0.0`, identity scale).
4. Survival, MD, IRR cells: byte-identical resolution to baseline.
5. Exactly **one** `effect_measure` resolution site remains (or two identical ones, if step 1 refuted
   unreachability), proven by the same search as step 1.
6. The campaign template's call re-resolves to the same estimand and thresholds as at baseline (it passes
   `"OR"` explicitly) — tested against the real file.
7. A replayed `args_call_all` list (built the way the bootstrap builds it) resolves the same estimand as its
   parent fit, for binary unset — before and after.
8. Identifier agreement (audit §3.8): after the change, a default binary call screens on the same estimand
   under `consistency`, `dina` and `grf`. State in the report that DINA already screened on log-OR, so B
   aligns FS and GRF to it — DINA's own behaviour is unchanged.
9. The `hr.threshold` roxygen no longer claims log scale for HR; rider commit is roxygen-only (assert
   mechanically as in the docs task).
10. `NEWS.md` entry present; check findings did not grow.

## Out of scope

Directive A entirely (validation, derivation, c2 ≤ c1 documentation); the `args_call_all` sync for the legacy
threshold spellings (audit §3.9); `adverse_outcome` at both sites; floors and admission criteria; Directive C;
anything in fs-glms-interpretable.
