# TASK: `mr_values()` — exported tidy accessor for `fit$mr_inference` (review proposal P5)

**Date:** 2026-11-06
**Origin:** REVIEW_wrapper_functions_2026-11-05.md (proposal P5), approved by Larry 2026-11-06
with the standing constraint: **current package functionality is preserved; the change is an
overall improvement.** That constraint is enforced by the gates below, not by judgment calls.
**Compute:** none — test fixtures are seconds-scale; no simulation go/no-go applies.
**Git:** work on the current checkout's branch; record branch and HEAD in the report. Commit as
you go. **Do not push** (Larry pushes via GitHub Desktop).

---

## Objective

One new exported function that extracts the multiplier-resampling inference quantities from a
`forestsearch()` fit (or a standalone `fs_mr_inference()` result) into a tidy data frame — one
row per region (Ĥ, and Ĥᶜ when the complement was computed). Pure extraction: every value
already exists in `mr_inference`; nothing new is computed (an optional log transform of
ratio-measure columns excepted).

The accessor exposes **only the evaluated construction set**: naive, two-term IJ (de-biased),
field (one-sided lower on Ĥ), field-s / field (one-sided upper on Ĥᶜ, whichever the fit
computed), and the Bonferroni pair — plus the fit-level diagnostics (selection bias, p̂(Ĥ),
family size, timing). **No winner-variant quantities** in any form.

## R/ impact callout (per protocol)

- **New code:** `R/mr_values.R` (one exported function), tests, roxygen, NEWS entry.
- **Moves existing code:** none.
- **Changes behaviour:** none. `mr_estimates_table()`, `fs_mr_inference()`,
  `forestsearch()` and every other existing `R/` file are untouched. (The option of having
  `mr_estimates_table()` delegate its `"mr_values"` attribute to the new accessor is
  **explicitly out of scope** for this task.)
- **Changes method:** none.

## Out of scope

- Any edit to an existing `R/` file's logic (regenerated `NAMESPACE`/`man/` and `NEWS.md` are
  the only permitted touches outside the new files).
- Any vignette or article edit.
- Winner-only / winner-floor quantities.
- Any new statistic.

---

## Stage 0 — verify structure from source, record baseline

Prior summaries, the review document and the project KB are **not** evidence of current
structure. Read the executing code:

1. Record the current branch and `HEAD`, then run the existing test suite and
   `devtools::check()` once, and record the baseline counts (errors / warnings / notes, and
   which notes). This baseline is what Gate 2 compares against.
2. Read `R/fs_mr_inference.R` (and the `forestsearch()` pass-through that assembles
   `fit$mr_inference`, wherever it lives) and the source of `mr_estimates_table()` (including
   how its `"mr_values"` attribute is built).
3. Produce a **field map** with file:line references: for each quantity in D4 below, the exact
   source path inside `mr_inference`, the condition under which it exists (`ci_method`,
   `include_complement`, `field_complement`, `field_scale_complement`, `return_reselection`,
   `outcome_type` / `effect_measure`, …), and the scale it is stored on (effect scale vs
   log/link — note that point estimates appear to be stored on the effect scale while
   `selection_bias` is on the log/link scale for ratio measures; **verify, do not assume**).
   Both `joint` and `joint_s` Bonferroni blocks and both `upper_1s` and `upper_1s_s`
   complement bounds have appeared in committed article code — map whichever of these the
   current source actually produces, under which options.

**Gate 0 (stop-on-failure):** every D4 quantity has either a verified source path or a declared
absence condition (→ `NA` in the output). If a specified quantity has no source under any
configuration, or the structure contradicts this spec in a way that changes the API surface:
STOP and report — do not improvise an API.

---

## Decisions (pre-set; no mid-task questions)

- **D1 — name and input.** `mr_values(x, scale = "effect")` in new file `R/mr_values.R`.
  `x` is a `forestsearch` fit (uses `x$mr_inference`) or a bare `fs_mr_inference()` result;
  detect which, error informatively when neither / when `mr_inference` is absent.
- **D2 — shape.** Plain `data.frame`. Row `"H"` always; row `"Hc"` when the complement block
  exists. Fixed column superset; `NA` where a component was not computed (never drop columns).
  Fit-level scalars (selection bias, family size, p̂, timing) sit on the `"H"` row, `NA` on
  `"Hc"`.
- **D3 — scale.** `scale = c("effect", "log")`. `"effect"` (default) returns every column
  exactly as stored — no transformation. `"log"` applies `log()` to the point/interval/bound
  columns of ratio measures (HR/OR) and errors informatively for identity measures (MD).
  Columns whose stored scale differs from the rest (e.g. `selection_bias`) are returned as
  stored under both settings, with the scale stated per column in the roxygen. Final
  column-by-column scale statements come from the Stage-0 field map.
- **D4 — quantities** (exact column names finalized at Stage 0 from the verified source names;
  record the mapping in the report):
  n of the region; naive est / lower / upper; de-biased (two-term) est / lower / upper and its
  IJ SE; the field point estimate(s) where present; field one-sided lower bound (H row); field
  one-sided upper bound(s) for the complement — `upper_1s` and/or `upper_1s_s`, whichever the
  fit computed, as separate columns; Bonferroni pair (from `joint` and/or `joint_s`, as
  separate columns if both exist); selection bias; p̂(Ĥ) (when re-selection was returned);
  family size; n_selected; timing.
- **D5 — fixtures and tests** (`testthat`):
  - **F1 (full-coverage):** a small synthetic survival dataset with a small hand-built
    candidate family run through `fs_mr_inference()` directly (the standalone-call pattern of
    the committed Guo–He article is the transplant precedent), with `ci_method = "field"`,
    complement and complement-field on, re-selection returned, a few hundred draws, fixed seed.
  - **F2 (NA-fill):** the same data with `ci_method = "ij"` and no complement.
  - **F3 (fallback, only if the standalone call cannot produce some mapped component):** a
    minimal `forestsearch()` fit on a tiny dataset; CC's construction, deterministic seed.
  - Fixtures must each run in well under a minute locally; wrap anything over ~5 s in
    `skip_on_cran()`.
- **D6 — identity assertions (the substance of Gate 1):** for every fixture, every non-`NA`
  cell of `mr_values()` equals the value read directly from the raw `mr_inference` paths
  (`all.equal`, tolerance 1e-12); and where `mr_estimates_table()` runs on the fixture, the
  overlapping quantities equal its `"mr_values"` attribute. `NA` cells must correspond to
  genuinely absent components, asserted explicitly.
- **D7 — docs.** Roxygen with a runnable example mirroring F1 (`\donttest` if over ~5 s);
  the roxygen states, per column, what it is, its source component, when it is `NA`, and its
  scale. `NEWS.md` entry under the development version.

---

## Stage 1 — implement

Write `R/mr_values.R`, roxygen, tests per D1–D7; `devtools::document()`; run the new tests.

**Gate 1 (stop-on-failure):**
1. All D6 identity tests green, with the verbatim test output captured for the report.
2. `git diff --stat` confirms the change set is confined to: `R/mr_values.R`,
   `tests/testthat/test-mr_values.R` (and fixture helpers), regenerated `NAMESPACE` and
   `man/mr_values.Rd`, `NEWS.md`. Any line changed in an existing `R/` file fails the gate.

## Stage 2 — preservation check

Run the **full** existing test suite and `devtools::check()`.

**Gate 2 (stop-on-failure):**
1. Existing test suite: green, no skips introduced beyond baseline.
2. `check()`: 0 errors, 0 warnings, notes identical to or fewer than the Stage-0 baseline
   (no new notes).

## Stage 3 — report

`REPORT_mr_values_2026-11-06.md`, placed beside the repo's existing function-task `REPORT_*`
files (follow the location convention found at Stage 0; transplant-first). Contents:

- branch + HEAD at start and end; commit list;
- the Gate 0 field map (quantity → source path, condition, scale; file:line);
- fixture definitions and the **verbatim** identity-test output;
- baseline vs final `check()` summary, side by side;
- the improvement demonstration (non-gating): the ~25-line `s16_column`-style raw-path
  extraction and its `mr_values()`-based equivalent (~5 lines) run on F1, with their outputs
  shown equal — included in the report only, no article is edited;
- OPEN ITEMS block for documentation gaps (per protocol, these do not stop the run; work
  failures do).

No push. Larry reviews the report and pushes.
