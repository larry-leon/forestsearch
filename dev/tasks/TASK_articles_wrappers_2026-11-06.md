# TASK: article editing pass — the four articles moved onto the wrappers (review §6)

**Date:** 2026-11-06
**Origin:** REVIEW_wrapper_functions_2026-11-05.md §6, approved by Larry 2026-11-06. Closes the
wrapper workstream: P5, P1+P2 and the B1–B5 batch are merged and installed; this task points
the articles at them. Standing constraint: **every published number the articles print is
reproduced; code shrinks, output does not move.**
**Prerequisite (Gate 0 stop condition):** the wrapper batch is merged into `master` and the
installed library copy exports the batch functions. If either fails, STOP and report.
**Compute:** no simulation campaign and (expected) no R code change. The cost is renders:
two baseline control renders plus one final render of the four articles, ≈ 60–90 min wall in
total, plus single-article re-renders while editing. No compute go/no-go applies.
**Git:** task doc committed on `master`; work branch `fs-articles-wrappers`; **one commit per
article.** Do not push.

## R/ impact callout (per protocol)

**None.** No `R/` file, test, or `man/` page changes. The diff is confined to the four `.qmd`
files (plus any rendered-output or freeze files the repo's own conventions track, identified
at Stage 0). New code: none. Moves existing code: none. Changes behaviour: none. Changes
method: none.

## Scope — the four articles and exactly these replacements

### A1 `vignettes/articles/mrct-consistency.qmd` (application article)
- **Steps 1–5 stay unpacked** — the article is the record of every call; that purpose is
  unchanged.
- Replace the hand-coded plotting (`fig-case-forest`, `case-km`) with the decoupled `plot()`
  methods on the four-sets/analysis objects, and update the one sentence that currently
  explains the hand-coding ("drawn from these objects with the same code, kept short").
- Add a short closing subsection **"The one-call form"** after step 5: the same trial through
  `mrct_trial_analysis()` with `mr_inference = TRUE` and the record's seed, followed by a
  `stopifnot` that its Ĥ, four sets and MR values equal the step-by-step objects (render cost
  ≈ +1 min).
- Retarget the "*With your trial*" notes at `mrct_trial_analysis()` where they currently
  describe hand-assembly; the step-level notes about argument meanings stay.
- `rd()`/`dd` block → `mrct_record()`; the `in_Hhat` idiom → `fs_in_subgroup()` **keeping**
  the existing `stopifnot` identity checks around it.
- Extend the install-check `needed` vector with every new function the article now uses.

### A2 `vignettes/articles/mrct-consistency-simulations.qmd`
- `bridge-case` → `mrct_trial_analysis()` on the shipped `rep5_trial.csv` (the one-call
  form), with the downstream objects read from its result; **every existing committed-record
  `stopifnot` in the article stays and must still pass.**
- `own-summary` and `demo-summary` → `mrct_core_summary()`; the demonstration-vs-campaign
  side-by-side table keeps its exact row labels and the committed campaign column.
- `own-analytic`, `analytic-cells`, the inline `design_tab()`, `design-n700` →
  `ahr_truths_table()` / `ahr_design_table()`, keeping the checks against
  `predictions_full.csv` / `predictions_S1.csv` / `predictions_n700.csv` at 1e-9.
- `build_own` / `build_cell` → `mrct_build_dgm()` (warning capture via its attribute; the
  build-warning sentences keep rendering the same counts and messages).
- `rd()` block → `mrct_record()`; install-check `needed` vector extended.

### A3 `vignettes/articles/actg175-glm-outcomes.qmd`
- `s16_column()` → built on `mr_values()`; the parts that read the data rather than
  `mr_inference` (N, event counts, region membership) stay, with `in_region()` →
  `fs_in_subgroup()` and its membership `stopifnot` kept.
- Part 2: `three_columns()`/`side_by_side` → `mr_method_comparison()`; the cost chunk →
  `mr_cost_table()`. The captions' framing (bootstrap point only at small B; the two routes
  answer different questions for model-generated families) stays verbatim.
- `s17_column()` → `fs_declaration_table()`, with the published-S17 comparison table
  untouched.

### A4 `vignettes/articles/survival-post-selection.qmd`
- `three_columns()` → `mr_method_comparison()`; the cost chunk → `mr_cost_table()`.
- Nothing else: the short narrative `fit$mr_inference` reads in its prose sections stay —
  they are few lines and instructive, and chasing every `$` read is not this task.

## Out of scope

Any other article or vignette (`vignette("forestsearch")`, the 2024 vignettes, Guo–He,
paper_simulations, the theory documents); any `R/`, test or `man/` change; any change to the
articles' readings, published-table comparisons, or References; file moves between
`vignettes/` and `vignettes/articles/`.

---

## Stage 0 — locate, inventory, baseline control renders

1. Confirm the prerequisite: `master` contains the batch merge; the installed copy exports
   the nine batch functions plus `mr_values`, `mrct_trial_analysis`, `mrct_four_sets`.
2. Locate the four files; record for each whether it is a CRAN vignette or an Rbuildignored
   website article (this decides the final preservation step), and whether the repo tracks
   any rendered output or freeze files for them (follow the convention found; add nothing
   new).
3. **Protected-check inventory:** for each article, list every assertion against a shipped
   record, committed value or published table (file:line). These may not be deleted or
   weakened. Internal scaffolding assertions may go only together with the scaffolding they
   tested, each listed in the report with the batch-gate that now covers the same identity.
4. **Baseline control renders:** render each article **twice** at the unedited tree. Extract
   the numeric tokens from each rendered output; the tokens that differ between the two
   control renders (timings, dates, sessionInfo, wall-clock sentences) form that article's
   **volatile set**. Everything outside it is frozen.

**Gate 0 (stop-on-failure):** prerequisite holds; both control renders of every article
succeed (all internal checks pass on the unedited tree); volatile sets recorded; protected
inventory complete.

## Decisions (pre-set)

- **D1:** only the replacements listed above; prose edits only where a sentence describes
  code that changed or the install-check vectors — every text edit quoted in the report.
- **D2:** protected checks retained verbatim; the A1 one-call subsection's new `stopifnot`
  is the only new check added.
- **D3:** the render-equality gate compares numeric tokens, volatile sets excluded; zero
  differences required per article.
- **D4:** inner loop re-renders only the article being edited; the full four-article render
  happens once, at Gate R.
- **D5:** grep count of `stopifnot(` per article must be ≥ baseline minus the listed
  scaffolding removals; the report shows the per-article before/after counts.

## Stage 1 — edits (one commit per article, A1–A4)

## Stage 2 — Gate R: render equality (stop-on-failure, per article)

Render all four once at the final tree. For each article: (a) the render succeeds — every
retained committed-record check passed; (b) numeric tokens equal the baseline exactly outside
the volatile set; (c) A1's one-call `stopifnot` passed. Any difference outside a volatile
set: stop, fix, re-render that article only, and re-run Gate R for it.

## Stage 3 — preservation

1. **Diff confinement:** changed files are exactly the four `.qmd` (plus any convention
   files from Stage 0 item 2). Any `R/`, `tests/`, `man/`, `NAMESPACE` or `DESCRIPTION`
   entry in the diff fails the gate.
2. If all four files are outside the built package (Rbuildignored articles): the test suite
   and `check()` are **not run** — record the `.Rbuildignore` evidence instead. If any file
   is a CRAN vignette: one full `check()` at the end, 0/0/0 vs the reused baseline (the
   batch report's figures, after the usual empty-git-diff verification).

## Stage 4 — report

`REPORT_articles_wrappers_2026-11-06.md` beside the prior three. Per article: chunks
replaced with before/after line counts; every prose edit quoted; removed scaffolding
assertions with their covering batch gate; the volatile set; the Gate R token comparison
(counts compared, zero non-volatile differences); the protected-check inventory with
before/after `stopifnot` counts. Plus the Stage 3 record and the commit list. No push.

## Amendment 2 (2026-11-06, received in chat, verbatim)

Amendment 2 to TASK_articles_wrappers_2026-11-06.md: track ALL four HTML
renders, survival-post-selection.html included.
1. Check whether it is untracked by an explicit .gitignore rule or simply
   never added. If a rule excludes it, lift that rule for this file only
   (one-line .gitignore edit, permitted by this amendment).
2. After Gate R passes, commit its Gate-R render in its own commit, same
   pattern as the other three.
3. Stage 3 confinement is amended accordingly: the permitted change set
   now includes the four tracked HTML files (and .gitignore if item 1
   required it).
4. Record the finding (rule vs omission) and this amendment in the report.
Everything else unchanged.
