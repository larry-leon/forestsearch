# TASK — Align MR admission with the search's rounded rule, and measure it on GBSG

**Repo:** `larry-leon/forestsearch`, branch `feature/glm-extension` — run CC from the forestsearch clone.
Nothing is written to `fs-glms-interpretable`; its application `.qmd` may be read **read-only** for the MR
settings if they are not recorded in this repo.

**Purpose.** `R/fs_mr_inference.R:660-661` builds MR's admission floor from `qnorm((1 + p_star)/2)` — the exact
cutoff — while the search admits on `round(Pcons, pconsistency.digits) >= p_star`. MR therefore re-selects
under a different admission rule from the one that produced Ĥ. Theorem 1's correction evaluates the optimism
under the map the search actually ran, so this is a premise violation, not a conservatism preference, and it
runs in the under-adjusting direction. This task aligns MR and measures the effect on the GBSG application.

**`R/` CALLOUT — this CHANGES BEHAVIOUR, and it changes published results.** MR's admitted set changes for
candidates whose standardized statistic falls in `[z_eff, z_exact)`. Corrected estimates, field bounds and the
Bonferroni pair can all move. **Forest search only:** GRF and DINA carry no consistency screen, so `p_star`
must not enter their admission at all — Gate 2 makes that a hard requirement rather than an expectation.

**Kind:** package change plus one applied re-run. No simulation. **Compute:** six MR fits on n = 686 (three
identifiers, before and after), plus the test suite and `R CMD check`. Report the first MR fit's wall clock
before continuing. **Hard abort at 3 h.**

---

## 0. First action

1. Copy this file to `dev/tasks/TASK_mr_admission_alignment_2026-09-23.md`; `git add` that path; commit
   `docs(tasks): align MR admission with the rounded screen (2026-09-23)`.
2. Record HEAD, installed forestsearch version and build date, R version, platform.
3. `git status --short`: record pre-existing untracked files; never stage them.

No `git fetch` / `pull` / `push`. Explicit paths on every `git add`.

---

## 1. Baseline (before any edit) — the paired comparison depends on this

On the **unmodified** package, run the GBSG application with MR on, for **all three identifiers**, at the
application's own settings and **a fixed, recorded seed**.

- Take the MR settings from wherever they are recorded in this repo. If they are not, read
  `~/Documents/GitHub/fs-glms-interpretable/quarto/gbsg/analysis_gbsg_mr.qmd` **read-only** for the
  `forestsearch()` call, as the 2026-09-23 null task did. Quote the settings in the report.
- **The seed must be identical before and after.** With the same seed the multiplier draws are the same
  vectors, so every difference is attributable to the admission threshold and not to Monte Carlo noise. A
  comparison on different draws is void.
- Record per identifier: the region (rule, N, events), unadjusted HR, corrected estimate, field lower bound,
  field-s upper bound, Bonferroni pair, selection bias on the log scale, and the re-selection family size.
- Save to a scratch `.rds` under `~/Downloads`, not the repo, not committed.
- Report the first fit's wall clock before running the rest.

---

## 2. The change

Read `R/fs_mr_inference.R` before editing; do not work from this description.

- MR's admission floor must be built from the **effective** threshold of the rounded rule, not from
  `qnorm((1 + p_star)/2)`.
- **Use the existing helper `.fs_pcons_eff(p_star, digits)`** in `R/fs_declaration_calibration.R`. Package
  internals are visible across files, so call it — do not move it and do not reimplement the expression. The
  2026-09-23 alignment task made "computed in exactly one place" a post-condition; this task inherits it.
- **`digits`** comes from the same source the calibration uses: the fit's `pconsistency.digits`, falling back to
  the `subgroup.consistency()` default. Determine from source whether MR has the argument directly in scope
  (it runs inside `forestsearch()`) or must read `args_call_all`, and record which.
- **Apply on both consistency paths.** Under `consistency_method = "split"` the closed form is an
  approximation to `k / n_valid`, but MR already used that approximation; substituting the effective threshold
  corrects the rounding without changing the approximation. If reading the source shows MR's floor is built
  differently on the two paths, report that instead of forcing it.
- **Do not change `R/fs_declaration_calibration.R`.** Do not change any selection rule, the screen, or the
  rounding design.

---

## 3. The GBSG re-run

Repeat Step 1 exactly — same three identifiers, same settings, **same seed** — on the modified package.

Report one table, identifiers × quantities, with before, after, and the difference:

- region (rule, N, events), unadjusted HR, corrected estimate, field lower bound, field-s upper bound,
  Bonferroni pair (region / complement), selection bias, re-selection family size.

---

## 4. Draw-level diagnostic — how much the change actually bites

This is what tells Larry whether the Section 5 simulation re-run is a refresh or a finding, so it matters more
than the headline numbers.

For forest search, on the shared draws:

- The number and share of outer draws on which the **admitted set** differs between `z_exact` and `z_eff`.
- The number and share on which the **re-selected winner** differs.
- How many candidates of the re-selection family have an observed `T` in `[z_eff, z_exact)`, with their rules,
  N and HR. On the observed fit we expect `{pgr <= 0} & {er <= 9}` (HR 1.683928, `Pcons` 0.898655) to be one —
  **report what is actually found**, and say whether it matches.

If the winner differs on very few draws, say so plainly: that is the result, not a disappointment.

---

## 5. Gates (stop on failure)

- **Gate 1 — the change is real.** At least one of the reported FS quantities differs between Step 1 and
  Step 3, or, if none does, the Step 4 diagnostic shows why (no candidate in the band on any draw). A silent
  no-op with a populated band is a failure.
- **Gate 2 — GRF and DINA are untouched.** Their Step 3 outputs must be **identical** to Step 1, every
  reported quantity. They carry no consistency screen, so `p_star` must not reach their admission. A
  difference means the change leaked into a path it must not touch: stop and report.
- **Gate 3 — the helper is not duplicated.** `grep` confirms the effective-threshold expression appears in
  exactly one place in `R/`; MR calls the helper.
- **Gate 4 — test suite.** `devtools::test()` passes. Any test encoding MR's exact-cutoff admission must be
  updated to the aligned rule and **listed individually** with its old and new expectation.
- **Gate 5 — `R CMD check --as-cran`.** The baseline is **1 NOTE** ("Version contains large components
  (0.3.5.9000)"), deliberately open as the dev-version marker. The gate passes if the modified tree produces
  exactly that NOTE and nothing else. **Do not bump the version.**

---

## 6. NEWS, record and commits

- NEWS.md under the current development version: MR admission now uses the screen's rounded threshold;
  corrected estimates and field bounds can move for forest search; GRF and DINA unaffected.
- Report to `dev/reports/REPORT_mr_admission_alignment_2026-09-23.md`.
- Commits, explicit paths, in order: task doc; `R/` change and tests; regenerated `man/` if roxygen changed;
  NEWS.md; the report.

---

## POST-CONDITIONS (machine-checkable)

1. The same seed was used in Steps 1 and 3; it is recorded, and the multiplier draws are shared.
2. Gate 2: GRF and DINA outputs identical before and after, every reported quantity.
3. Gate 3: the effective-threshold expression appears in exactly one place in `R/`.
4. Step 4's three diagnostics are reported with counts and shares, not adjectives.
5. `devtools::test()` passes; every updated test listed with its old and new expectation.
6. `R CMD check --as-cran` produces exactly the one dev-version NOTE.
7. `R/fs_declaration_calibration.R` is unmodified.
8. Files modified are confined to `R/fs_mr_inference.R`, `tests/testthat/`, `man/`, `NEWS.md`, `dev/tasks/`
   and `dev/reports/`.
9. Nothing written to `fs-glms-interpretable`.

---

## OUT OF SCOPE

No re-run of the Section 5 simulation cells — that is a separate decision once the GBSG magnitude is known. No
ACTG 175 re-run. No change to the rounding design, the screen, κ̂, FŴ, or any selection rule. No fix to the
consumers that read the removed `pstar_implied` (`declcalc0_run.R`, `dev/verification/report_values*.R`) —
separate and not blocking here. No manuscript edits: Section 4.6's Eq. (8), Section 2.3.3's Step 3 and
Section 4.4 all state the exact cutoff and now need revising, but that drafting belongs to the manuscript chat.
