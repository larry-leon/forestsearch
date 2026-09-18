# CC TASK — threshold sync: replicates resolve what the parent fit resolved

**Opened:** 2026-09-18 · **Repository:** forestsearch · **Authorized by:** Larry, 2026-09-18.
**Evidence base:** `REPORT_criterion_and_defaults_audit_2026-09-18.md` §3.9 (pointers, not facts — re-verify
from source; line numbers there predate the docs and B commits and HAVE drifted).
**Compute:** none. No fit, no bootstrap, no CV is executed. Probes construct and resolve argument lists only.
**Testing (standing policy):** ONLY this task's own acceptance test files, hard abort at 3 minutes. No R CMD
check, no full suite.
**Every step has a hard cap. Nothing runs open-ended.**

## The defect being fixed

Resolved thresholds live in locals and are never written back. `args_call_all` captures the formals, and
bootstrap (`bootstrap_analysis_dofuture.R`, `args_FS_boot`) and CV replay it with every formal supplied — so
`missing(hr.threshold)` / `missing(hr.consistency)` are FALSE inside a replicate and `user_set_*` is wrongly
TRUE. Consequence today (audit §3.9): on identity-scale estimands a replicate resolves different thresholds
than its parent fit (e.g. consistency 1.0 instead of 0.0). Ratio measures and survival are currently
unaffected only by accident of the ratio branch — and Directive A's derivation cannot land safely until
replicates provably resolve what the parent resolved.

## Classification

| Change | Class |
|---|---|
| Write resolved thresholds back so replays resolve identically | **changes behaviour** (bugfix on replicate/CV resolution; parent-fit resolution must NOT change) |
| Rider (separable): binary default RD → OR at the three exported estimation entry points | **changes behaviour** (a default), zero committed exposure per B's inventory |

## Rules

Copy this document into `dev/tasks/` and commit it first. No install; `devtools::load_all()` only. Do not
push. Gates are stop-on-failure. Find everything by search. Report ≤ ~120 lines, cite file + function (line
numbers optional), filed at `dev/reports/REPORT_threshold_sync_2026-09-18.md` (create `dev/reports/`; this
starts the convention for package-level reports).

## Steps

### Step 1 — source verification + exposure inventory (cap: read-only, no limit needed)

1. Locate at current HEAD: the alias merge; the `user_set_*` detection; the estimand resolution block; the
   locals holding resolved values; `.sync_args_call_all()` (signature and what it writes); where
   `args_call_all` is captured; how `args_FS_boot` is built and replayed; the CV equivalent.
2. **Exposure inventory (Gate 1):** every committed in-repo caller that (a) uses an identity-scale estimand
   (RD, IRD, MD), AND (b) leaves thresholds at defaults or supplies them via the legacy spellings only, AND
   (c) runs bootstrap or CV. **If any exists, STOP and report** — Larry decides pin-vs-scope before any edit.
   (Expected empty: the OR campaign uses the NULL-defaulted spellings; the MD campaigns pass thresholds
   explicitly — verify the spelling from their files.)

### Step 2 — baseline: the replicate-equality probe (committed before any edit)

For every cell of: estimand {survival, binary-unset(→OR), binary "RD", binary "OR", MD, IRR, IRD} ×
threshold supply {none, legacy spellings, new spellings} × {c1 only, both, neither} — construct the parent
argument list, resolve it; construct the replicate list exactly as `args_FS_boot` is built from
`args_call_all`, resolve it; record both. **The violation table (cells where replicate ≠ parent) is the bug,
committed as the baseline with the HEAD SHA.** Expect violations only on identity-scale cells; if ratio or
survival cells violate at baseline, STOP and report.

### Step 3 — the sync

Implement write-back so that after resolution, a replay resolves identically. Suggested mechanism (CC may
improve; the post-conditions govern, not the mechanism): sync the resolved values, on the natural scale, into
the NULL-defaulted spellings (`effect.threshold`, `consistency.threshold`) in `args_call_all` via
`.sync_args_call_all()` — those spellings are `is.null()`-detected and therefore wrapper-safe, and the alias
merge then makes the replicate honour them. Do not alter the resolution logic itself.

### Step 4 — the rider (SEPARABLE; skip on Larry's word without touching steps 1–3)

At the three exported entry points from B's report (F2 / caller inventory): `make_effect_estimator()`,
`consistency_resample()`, `consistency_resample_compare()` — the internal binary default becomes `"OR"`,
matching `forestsearch()`. Own commit. Exposure gate: re-run B's inventory grep and confirm still no
committed caller reaches those defaults; cite B's report. Update the three `@param` docs accordingly.

### Step 5 — re-probe and gates

Re-run the step 2 probe on the edited tree.
- **Gate 5a:** zero violations — every cell's replicate resolution equals its parent's.
- **Gate 5b (parent invariance):** every cell's PARENT resolution is byte-identical to baseline. The sync may
  change only what replays see, never what the original fit resolves.
- **Gate 5c:** ratio and survival replicate cells are byte-identical to baseline (they were already equal;
  they must stay exactly as committed work produced them).

### Step 6 — acceptance tests (hard abort: 3 minutes)

As testthat files under `tests/`, run via `testthat::test_file` + `load_all()`:
1. Replicate-equality across the step 2 matrix (the probe, as a test).
2. Parent-fit resolution unchanged across the matrix.
3. A replayed list built the bootstrap's way carries non-NULL `effect.threshold` / `consistency.threshold`
   equal to the parent's resolved natural-scale values.
4. Rider (if kept): the three entry points resolve binary→OR when unset; explicit values untouched.
Record wall clock.

### Step 7 — NEWS and report

`NEWS.md`: one bugfix entry (replicates/CV now resolve the parent's thresholds; who was affected:
identity-scale analyses using defaults or legacy spellings with bootstrap/CV) and one entry for the rider if
kept. Then the report: gate outcomes, the baseline violation table, the inventory result, findings (no tasks
attached). Commit. **Do not push.**

## Out of scope

Directive A (validation, derivation, c2 ≤ c1 doc) — next task, on top of this one. Directive C (DINA guard,
frontier warnings, display defaults). Floors. Anything in fs-glms-interpretable.
