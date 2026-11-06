# TASK: `mrct_trial_analysis()` + `mrct_four_sets()` on a shared core (review proposals P1+P2, Option B)

**Date:** 2026-11-06
**Origin:** REVIEW_wrapper_functions_2026-11-05.md (P1 flagship + P2), Option **B** approved by
Larry 2026-11-06: shared core — `mrct_case_study()` is rebuilt internally on the same analysis
core the new real-trial entry point uses, so the simulated and observed routes provably run
identical code. Standing constraint: **current functionality preserved; overall improvement** —
enforced by the identity gates below.
**Compute:** no simulation campaign. The heaviest gate call is one S2 build + one replicate-5
case study with the MR field call (the application article's whole run is ~5–6 min); total gate
compute well under 15 min. No compute go/no-go applies.
**Git:** commit the task doc on `master`, then create work branch `fs-mrct-trial` from it and
make every other change there (the P5 pattern). Commit as you go. **Do not push.**

---

## Objective

1. **`mrct_trial_analysis()`** (new export): the regional-consistency pipeline of the
   application article run on an **observed** trial — ITT fit, region fit and the trigger; the
   identifier on the non-region patients with the region as `df.test`; Ĥ applied; the
   four-set Cox summary; MR inference carried on the fit. The real-trial counterpart of
   `mrct_case_study()`.
2. **`mrct_four_sets()`** (new export): the four-set constructor with `print()` and
   `plot(which = c("forest", "km"))` decoupled from the `mrct_case_study` class — the plotting
   the application article currently re-implements by hand because the methods are class-bound.
3. **Shared core (Option B):** the analysis/report portion of `mrct_case_study()` is
   **extracted** (transplanted, not re-authored) into an internal core; `mrct_case_study()`
   keeps its public signature, return structure, class vector, `print()` and `plot()` output
   **exactly** as they are, and calls the core; `mrct_trial_analysis()` is a second public
   entry to the same core for observed data.

## R/ impact callout (per protocol)

- **New code:** `mrct_trial_analysis()`, `mrct_four_sets()` + methods, tests, roxygen, NEWS.
- **Moves existing code: YES** — the analysis-generic portion of `mrct_case_study()` (and its
  `print`/`plot` workhorse code) moves into the shared core. This is the task's only touch to
  existing functions, and it is gated to be behaviour-preserving (Gate 1).
- **Changes behaviour:** none. The case-study route must be output-identical (objects, printed
  text, plotted quantities) before and after.
- **Changes method:** none.

## Out of scope

- Any vignette or article edit (the articles move to the wrappers in a later editing pass).
- `mrct_region_sims()` and its summarization (review P3) — separate task.
- Any change to `forestsearch()`, `fs_mr_inference()`, the `ahr_*` layer, or any `R/` file
  other than the one(s) housing `mrct_case_study()` and its methods.
- Any new statistic or any addition of fields to `mrct_case_study()`'s return object.

---

## Stage 0 — source map and baseline capture

Verify from source; the review, the articles and prior summaries are not evidence of current
structure.

1. Record branch + HEAD; run the full test suite and `devtools::check()` once; record the
   baseline (errors / warnings / notes, failures, the skip list, expectation count).
2. Read the file(s) housing `mrct_case_study()` and its `print`/`plot` methods. Produce a
   **line-partition map**: which lines are simulation-specific (DGM redraw, seed = call seed +
   k bookkeeping, `call_args`, generator save/restore) and which are analysis-generic (the
   identifier call on non-region with `df.test`, Ĥ membership, the four-set Cox summaries,
   confusion tables, object assembly, print/plot internals). File:line references required.
   Record also: the exact class vector of the return object; how the fs call is seeded on the
   case-study route; the sets' row set, ordering and labels; which components are conditional
   (MR on/off, which = forest/km).
3. Confirm the gate data dependencies are in the repo: `inst/extdata/mrct-consistency/`
   (`rep5_trial.csv`, `mr_field_S2_rep5.csv`, `s2_example_four_sets.csv`,
   `s2_example_trials.csv`, `s2_example_definitions.md`, `latent_R.csv`, `PROVENANCE.md`) and
   that the S2 build can run with the shipped `latent_R` (no live latentcor fit), per the
   simulations article's `build-fn` chunk — transplant that build code for the fixtures.
4. **Baseline capture:** at baseline HEAD, run and `saveRDS()` (outside the package tree, in
   the task's scratch area):
   - **FX-fast:** a small synthetic DGM (small `inject_mrct_structure()` build on a reduced
     dataset or the smallest committed-precedent equivalent), one `mrct_case_study()` call with
     `mr_inference = FALSE` — seconds. Save the full return object and `capture.output(print(x))`.
   - **FX-rep5:** the S2 build (shipped `latent_R`) + `mrct_case_study(dgm_S2, seed = <S2
     campaign seed>, k = 5, …, mr_inference = TRUE)` with the simulations article's arguments —
     the article's own reproduction fixture. Save the object and print text.

**Gate 0 (stop-on-failure):** the partition map is complete, the class vector and seeding are
recorded, the gate data are present, and both baseline captures saved. If the analysis code is
interleaved with simulation bookkeeping in a way that cannot be separated without changing
behaviour: STOP and report — do not improvise a riskier extraction.

---

## Decisions (pre-set; no mid-task questions)

- **D1 — names and files.** `mrct_trial_analysis()` and `mrct_four_sets()`, in new files
  following the repo's mrct file conventions; the internal core named per repo convention
  (e.g. an unexported `.mrct_analyze_trial()`). The core is the **moved** case-study code.
- **D2 — `mrct_case_study()` frozen.** Public signature, return fields, class vector, print
  text and plotted quantities unchanged; no new fields, even harmless ones. Its `print`/`plot`
  become thin shims over the shared workhorses, dispatching on the existing class.
- **D3 — `mrct_trial_analysis()` conventions mirror `mrct_case_study()`.** Variable-name
  arguments and the `fs_args` shape exactly as the case study takes them, plus:
  `region.name` (the observed region column), optional `truth.name` (a flag_harm-style column:
  adds the G* row and the confusion tables; absent → neither, no error), `mr_inference` /
  `mr_inference_args` as on the case study, and `proceed_without_trigger = FALSE`.
- **D4 — stage semantics, never `stop()`.** The ITT and region fits are always computed.
  `$stage` ∈ {`not_positive`, `no_trigger`, `no_declaration`, `complete`}. Below `complete`
  the object carries the fits and the sets computable so far; the identifier runs on a
  non-triggered trial only when `proceed_without_trigger = TRUE`.
- **D5 — sets.** Row set, ordering and labels exactly those the current `mrct_case_study()`
  produces; the G* row and confusion only with truth. `mrct_four_sets()` returns the same
  table shape with `as.data.frame()`, `print()`, `plot()`.
- **D6 — print with MR present** reports the evaluated construction set only — naive,
  two-term IJ, field lower on Ĥ, field / field-s upper on Ĥᶜ, Bonferroni pair — via
  `mr_values()` where convenient. Bounds are labeled by location (no significance framing at
  HR = 1.0, no "significant"/"n.s." wording anywhere in methods output).
- **D7 — new class.** `mrct_trial_analysis()` returns its own class (not
  `"mrct_case_study"`); both classes' methods call the shared workhorses.
- **D8 — docs.** Roxygen for both exports with small runnable examples (`\donttest` if > ~5 s);
  `mrct_case_study()` roxygen edited only where it describes internals that moved (text only);
  NEWS entries; `_pkgdown.yml` index entries for both new exports beside the mrct family (in
  scope this time — the P5 lesson), `pkgdown::check_pkgdown()` clean.
- **D9 — tests committed** (`testthat`): FX-fast identity (regenerated in-test, not the RDS);
  `mrct_trial_analysis()` stage-semantics on tiny engineered trials (one per stage); a
  truth-absent run (no G* row, no confusion, no error); `mrct_four_sets()` on a hand-built
  data frame against directly computed `coxph()` summaries at 1e-12. Anything over ~5 s wrapped
  in `skip_on_cran()`. The rep5 gates run once as gate artifacts (report), not as suite tests.

---

## Stage 1 — the refactor (Gate 1: case-study identity)

Extract the core per the Stage-0 map; rebuild `mrct_case_study()` on it.

**Gate 1 (stop-on-failure):**
1. **FX-fast:** refactored return object vs baseline RDS — `all.equal(tolerance = 1e-12)` with
   attributes and class checked, and `identical()` on every non-floating component; print text
   byte-identical; both `plot()` types execute cleanly.
2. **FX-rep5:** the refactored route reproduces the committed record at the simulations
   article's own tolerances — **transplant the article's check code** (`find-check`,
   `reproduction`, `mr-repro` chunks), do not re-author it: redrawn trial equal to
   `rep5_trial.csv` (per-column `all.equal`, tol 1e-12); Ĥ string identical to
   `s2_example_trials.csv`; four sets' n and per-arm events identical and HRs equal to
   `s2_example_four_sets.csv` at the article's precision; both confusion 2×2s identical to
   `s2_example_definitions.md`; MR values equal to `mr_field_S2_rep5.csv` at 1e-8. Also
   object-vs-baseline-RDS equality as in (1).

## Stage 2 — the new entry points (Gate 2: route equivalence)

Implement `mrct_trial_analysis()` and `mrct_four_sets()` on the core, with D3–D7 and the tests
of D9.

**Gate 2 (stop-on-failure):**
1. **Route equivalence on replicate 5:** `mrct_trial_analysis()` run on the shipped
   `rep5_trial.csv` **as an observed trial** (region = `z_region`, `truth.name = "flag_harm"`),
   with the application article's identifier arguments and the record's call seed (mirror the
   article's `case-fs` call, including its MR seeding — Stage 0 records how), must equal the
   case-study route: Ĥ identical, four sets at 1e-12, confusion identical, MR values at 1e-8
   against `mr_field_S2_rep5.csv`. This is Option B's claim — real route ≡ simulated route on
   the same patients — demonstrated on the committed record.
2. D9 tests all green (stages, truth-absent, four-sets vs direct `coxph`).

## Stage 3 — preservation (Gate 3)

1. `git diff` confinement: changed existing files are only the one(s) housing
   `mrct_case_study()` and its methods (plus `NAMESPACE`, `NEWS.md`, `_pkgdown.yml`,
   regenerated `man/`); every other existing `R/` file untouched.
2. Full existing test suite green; skip list unchanged from baseline.
3. `devtools::check()`: 0 errors, 0 warnings, notes identical to or fewer than baseline.
4. `pkgdown::check_pkgdown()` clean.

## Stage 4 — report

`REPORT_mrct_trial_analysis_2026-11-06.md`, beside the P5 report (same location convention):
branch/HEAD and commits; the Stage-0 partition map (file:line); baseline-capture inventory;
**verbatim** Gate 1–3 outputs (the rep5 comparisons with their largest differences, as the
simulations article prints them); a before/after sketch of `mrct_case_study()`'s structure
(lines moved vs lines remaining); OPEN ITEMS for anything non-blocking. No push.

---

## Amendment 1 (Larry, 2026-11-06, in chat; arrived after Stages 0-4 had run)

Amendment to TASK_mrct_trial_analysis_2026-11-06.md, effective immediately:
1. Baseline (Stage 0): do NOT run the full suite or check(). Verify by
   git diff that nothing R CMD check sees (R/, tests/, man/, vignettes/,
   DESCRIPTION, NAMESPACE) changed since the P5 gated commit 90eed9b9;
   then reuse REPORT_mr_values_2026-11-06.md's final figures as the
   baseline and record that reuse in the report.
2. Inner loop: while developing, run ONLY the targeted test files and the
   FX-fast fixture. Never run the full suite or check() between edits.
3. The heavy FX-rep5 fixture runs exactly three times total: baseline
   capture, Gate 1, Gate 2. Do not re-run it per edit.
4. Full suite + check(): ONCE, at Stage 3 only.
All gates, tolerances and deliverables otherwise unchanged.
