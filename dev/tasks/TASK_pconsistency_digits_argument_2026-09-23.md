# TASK — Expose `pconsistency.digits` as a `forestsearch()` argument

**Repo:** `larry-leon/forestsearch`, branch `feature/glm-extension` — run CC from the forestsearch clone.

**Purpose.** `subgroup.consistency()` and the consistency evaluators already take `pconsistency.digits`
(default 2), but `forestsearch()` has no such argument, no `...`, and does not pass one through. Every
`forestsearch()` run is therefore pinned to 2 digits with no way to change it. This task makes the setting
reachable.

**`R/` CALLOUT — this ADDS AN ARGUMENT. It does not change behaviour and does not change the method.**
At the default of 2 every result must be bit-identical to the current package. Gate 1 enforces that. Nothing
about the rounding design is altered: the rounded value continues to decide admission, and the selection sort
continues to order on it.

**Kind:** package change. No simulation. **Compute: negligible** — test suite and `R CMD check` only. Hard
abort at 30 min.

---

## 0. First action

1. Copy this file to `dev/tasks/TASK_pconsistency_digits_argument_2026-09-23.md`; `git add` that path; commit
   `docs(tasks): add pconsistency.digits argument task (2026-09-23)`.
2. Record HEAD, installed forestsearch version and build date, R version, platform.
3. `git status --short`: record pre-existing untracked files; never stage them.

No `git fetch` / `pull` / `push`. Explicit paths on every `git add`.

---

## 1. Baseline capture (before any edit)

Run the GBSG application as `quarto/simulations/gbsg_app_null/run_gbsg_app_null.R` (`e6477a22`) sets it up, at
its own `p* = 0.90`, on the **unmodified** package. Save the full `grp.consistency$out_sg` result table to a
scratch `.rds` under `~/Downloads` (not the repo). This is the comparator for Gate 1 and is not committed.

---

## 2. The change

Read the current source before editing; do not work from this description.

- **`R/forestsearch_main.R`** — add a formal argument `pconsistency.digits = 2` to `forestsearch()`, in a
  position consistent with the surrounding threshold arguments, and pass it to the `subgroup.consistency()`
  call (currently around `:3490-3512`). Include it in `args_call_all` so a fit records what it received.
- **No change to `subgroup.consistency()`, the evaluators, or `.consistency_via_splits()`** — they already
  accept and honour the argument.
- **No change to any rounding logic, comparison, or selection sort.**

If passing the argument requires touching the `.make_eval_*` factories, report that rather than assuming; they
are documented as passing it through already.

---

## 3. Roxygen — correct the misleading wording

The existing roxygen describes this as digits "for output" (`R/subgroup_consistency_helpers.R:1176, 1304, 1584`
and `R/subgroup_consistency_main.R:188`). That is wrong: the rounded value decides admission against
`pconsistency.threshold` and is what the selection sort orders on.

Rewrite those descriptions, and write the new `forestsearch()` entry, to state plainly:

- the consistency proportion is rounded to this many digits **before** it is compared with
  `pconsistency.threshold`, so the setting affects which subgroups are admitted, not only what is displayed;
- at the default of 2, a threshold of `p*` admits any candidate whose unrounded proportion is at or above
  the nearest 0.01 grid point below `p*` plus half a step — so `p* = 0.99` admits from about 0.985;
- raising the setting makes admission finer-grained and produces fewer exact ties for the selection rule to
  break;
- under `consistency_method = "resample"` the proportion is a continuous closed-form quantity, so the rounding
  is a deliberate coarsening rather than a formatting step.

Do not editorialise beyond this; state the behaviour.

---

## 4. NEWS

Add a NEWS.md entry under the current development version: `forestsearch()` gains `pconsistency.digits`
(default 2), previously reachable only on the internal consistency functions; default behaviour unchanged.

---

## 5. Gates (stop on failure)

- **Gate 1 — behaviour unchanged at the default.** Re-run Step 1's GBSG fit on the modified package, with
  `pconsistency.digits` **omitted**. The `out_sg` result table must be identical to the Step 1 baseline —
  every row, every column, including `Pcons` — and the declared subgroup must be `{er <= 0} & {pgr <= 26}`,
  N 75, HR 2.22. Repeat with `pconsistency.digits = 2` passed explicitly; same requirement.
- **Gate 2 — the argument arrives.** `args_call_all` records `pconsistency.digits` for both runs above.
- **Gate 3 — the pass-through is real.** Run the same fit with `pconsistency.digits = 6`. Report the `Pcons`
  of the declared subgroup and of the maximum-effect candidate, and whether either is on the 0.01 grid. Report
  the declared subgroup, N and HR. **No expected value is asserted here** — record what comes back.
- **Gate 4 — test suite.** `devtools::test()` passes. Report `test-search-reproducibility.R` explicitly: it
  compares whole `Pcons` tables across worker counts, so if it passes only because rounding masks differences,
  that would show as a failure once digits are raised — run that file additionally with digits = 6 and report
  the outcome as a finding either way. **A failure at digits = 6 is not a gate failure**; a failure at the
  default is.
- **Gate 5 — `R CMD check`.** Clean, no new NOTEs, WARNINGs or ERRORs against the pre-change baseline.

---

## 6. Commits

Explicit paths, in order: task doc; `R/` change plus roxygen; regenerated `man/` pages; NEWS.md; the
verification record beside the existing `REPORT_*` files.

Do not commit the Step 1 baseline `.rds`.

---

## POST-CONDITIONS (machine-checkable)

1. Gate 1: `out_sg` identical to baseline with the argument omitted, and again with it set to 2.
2. Gate 2: `args_call_all` carries `pconsistency.digits` on every fit.
3. Gate 3: digits = 6 runs without error; the reported values are recorded.
4. Gate 4: `devtools::test()` passes at the default. The digits = 6 outcome for
   `test-search-reproducibility.R` is reported.
5. Gate 5: `R CMD check` clean against baseline.
6. `man/` regenerated and consistent with the roxygen.
7. Files modified are confined to `R/forestsearch_main.R`, the roxygen in the consistency files named in §3,
   `man/`, `NEWS.md`, `dev/tasks/` and the verification record.
8. No change to any rounding, comparison or selection logic.

---

## OUT OF SCOPE

No change to the rounding design — the rounded value continues to decide admission, by decision on 2026-09-23.
No change to `fs_declaration_calibration.R`, to κ̂ or FŴ, or to MR admission. No re-run of the GBSG p\* grid
campaign. No simulation of any kind. No change to the default, which stays 2.
