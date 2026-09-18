# CC TASK — Directive A: c2 > c1 fails loudly; silent c2 derives 0.80 · c1

**Opened:** 2026-09-18 · **Repository:** forestsearch · **Authorized by:** Larry, 2026-09-18.
**Prerequisite, already landed:** the threshold sync (`REPORT_threshold_sync_2026-09-18.md`) — replicates
resolve what the parent resolved, so a derived c2 written back through the same mechanism is what every
replay sees. Verify SECTION 2B-ii exists at HEAD before starting.
**Evidence base:** the audit report and the sync report in `dev/reports/` — pointers, re-verify from source.
**Compute:** at most two seeded micro-fits inside the acceptance tests (smallest data that completes a
binary-OR and a survival fit), bootstrap/CV off. **Hard abort 5 minutes.** Everything else is
argument-list resolution.
**Testing:** ONLY this task's own acceptance test files, hard abort 3 minutes (fits under their own cap).
No R CMD check, no full suite.

## Scope — dispositions already taken (these govern)

1. **HR and binary OR only.** RD, IRD, MD and IRR: no error, no derivation, resolution unchanged, byte-identical
   to baseline in every probe cell.
2. **The error and the derivation act only where the consistency stage runs** — `subgroup_method = "consistency"`.
   Under `"dina"` and `"grf"` c2 is inert: no error, no derivation, exactly as today.
3. **c2 > c1 is a `stop()`**, message naming values and the spellings the caller used, e.g.
   `c2 > c1 not allowed for FS: consistency.threshold = 1.00 exceeds effect.threshold = 0.90`. c2 = c1 passes.
4. **c1 supplied with c2 not supplied derives c2 = 0.80 × c1 on the ratio scale** — an additive shift on the
   log scale, never `0.8 × log(c1)`. Worked values: 1.25 → 1.00, 1.00 → 0.80, 0.90 → 0.72. Announced with a
   `message()` naming the derived value and the c1 it came from; written back through the sync so
   `args_call_all` carries it and every replay resolves it.
5. **Explicit beats derived.** A supplied c2 (either spelling) is never overridden.
6. **Both spellings of one quantity supplied and disagreeing is an error** naming both.
7. **No method change.** Selection logic untouched.

## Classification

| Change | Class |
|---|---|
| `stop()` on c2 > c1 (consistency path, HR + binary OR) | changes behaviour (validation) |
| Derivation of silent c2 (same scope), announced and synced | changes behaviour |
| Deleting the branch made unreachable by the validation | removes code (after proof) |
| roxygen: the c2 ≤ c1 requirement and the derivation, documented in the same commit that makes them true | documentation |
| Rider (separable): reorder `make_effect_estimator()`'s binary `match.arg` choices so `"OR"` is first (sync report S3 — latent, unreachable today) | alignment, no behaviour change |

## Rules

Copy this document into `dev/tasks/` and commit it first. No install; `devtools::load_all()` only. Do not
push. Gates stop-on-failure. Find everything by search. Report ≤ ~120 lines, file + function citations, to
`dev/reports/REPORT_directive_A_2026-09-18.md`.

## Steps

### Step 1 — source verification

Locate at HEAD: the `user_set_*` detection and alias merge; the ratio-branch resolution; SECTION 2B-ii (the
sync); the consistency-stage entry condition (`has_subgroups`, the `any(hr_values > check_threshold)` guard)
and exactly what runs when it is false; where `subgroup_method` is known at resolution time. Confirm
"c2 not supplied" is detectable there for both spellings (the sync's own analysis says it is — `is.null` on
the new spellings, `missing()` on the legacy ones, evaluated in `forestsearch()`'s own frame before any
replay). **Gate 1: stop if any of this differs from the sync report's description.**

### Step 2 — baseline probe, before any edit

Extend the sync's probe matrix (reuse `helper-threshold-sync.R`; keep one copy) with derivation and
validation cells, and record parent + replicate resolution for every cell at the pre-edit SHA, committed:

| Cell family | Expectation after the change |
|---|---|
| HR / binary-OR: c1 supplied (each spelling, several values incl. 0.90), c2 not | c2 = 0.80·c1, announced, parent = replicate |
| HR / binary-OR: both supplied, c2 < c1 and c2 = c1 | unchanged, no message |
| HR / binary-OR: both supplied, c2 > c1, `subgroup_method = "consistency"` | `stop()`, message as specified |
| same, `subgroup_method = "dina"` / `"grf"` | **no error, no derivation** — byte-identical to baseline |
| c2 supplied, c1 not | honored against default c1; error rule applies |
| both spellings of one quantity, disagreeing | error naming both |
| every RD / IRD / MD / IRR cell from the sync matrix | byte-identical to baseline |
| survival + binary-OR defaults (nothing supplied) | byte-identical — the default pair already satisfies the rule |

### Step 3 — implement

Validation first, then derivation, one commit each. The derivation writes the derived c2 (natural scale)
through the existing sync so replays inherit it. The `message()` fires once per fit, in the parent frame
only — verify it cannot fire per replicate (the replay carries the derived value explicitly, so the
derivation branch is not re-entered; assert this in a test).

### Step 4 — delete the unreachable branch

With c2 ≤ c1 enforced on the consistency path, the interval (c1, c2) is empty and the skip branch at the
consistency-stage entry cannot be reached by any legal call. Prove it from source (the audit's finding 11:
the entry condition binds only when c2 > c1), then delete, own commit. **Gate: if any legal call still
reaches it — including dina/grf paths, where no error fires — do not delete; report.**

### Step 5 — documentation, in the same commits as the behaviour

The roxygen sentences deferred from the docs task now become true and land here: c2 must be ≤ c1 on the
consistency path (error otherwise); a silent c2 is derived as 0.80 × c1 for ratio estimands, announced.
Update the threshold-vocabulary section's wording where it must now mention the derivation. `document()`.

### Step 6 — re-probe, gates

- **Gate 6a:** every cell matches the Step 2 expectation table; any other movement fails.
- **Gate 6b:** parent-resolution invariance for every non-derivation, non-error cell.
- **Gate 6c:** replicate = parent everywhere, including derived-c2 cells.

### Step 7 — acceptance tests

testthat files (3-minute abort): the full Step 2 matrix as assertions; the error message text (both
spellings); derivation announced exactly once per fit; the two micro-fits (5-minute cap): a binary-OR fit
with c1 = 0.90, c2 unset completes with the consistency stage receiving 0.72, and a c2 > c1 call errors
**before** any model is fit; rider (if kept): choice-order test.

### Step 8 — NEWS and report

`NEWS.md`: one entry for the validation + derivation (who is affected: callers who set c1 without c2 —
previously silent 1.0, now derived 0.80·c1; and callers with c2 > c1 — previously degenerate, now an error).
Report to `dev/reports/`, gates + findings, no tasks attached. Commit. **Do not push.**

## Out of scope

Directive C (DINA guard, frontier warnings, display defaults, dina/grf display annotation) — next task.
RD / IRD / MD / IRR behaviour. Floors. fs-glms-interpretable. `fpr_calibration()`'s own c2 ≤ c1 check
(already present and already documented).
