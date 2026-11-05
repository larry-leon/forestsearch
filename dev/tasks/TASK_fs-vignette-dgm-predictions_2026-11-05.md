# TASK_fs-vignette-dgm-predictions_2026-11-05

**Dispositions on record (Larry, 2026-11-05 12:59)**: the mrct-consistency article gains, in its Overview section, a "The data-generating mechanism" subsection and, before the analysis, the steps producing the report's Table 2 analytic predictions (S2 column only, rows (1)–(7)) **from `rep5_trial.csv`**; exact reproduction of Table 2 is not expected — the rep5-based build may differ, and the article says so.

**Repo**: `forestsearch` only. **No `R/` changes of any kind** — this task edits one vignette file. If any function or argument the insert uses is missing from the installed dev package, **STOP** and record it (an `R/` change is a separate proposal, never done here).
**Branch**: short-lived from `master`, fast-forward merge, delete branch (check `git branch -r` first). Record `master` HEAD before work (last known c6a18311, version 0.4.1.9002; verify from source).
**Compute**: none beyond the article's own render (the insert adds one copula build + quadrature, seconds; no new simulation).
**Session start**: `export GIT_TERMINAL_PROMPT=0`. Credentials-failed fetch recorded, not a stop.
**Files from ~/Downloads**: this task doc and `vignette-insert_dgm-predictions_2026-11-05.qmd` (the insert, verbatim content). **First action**: copy both into `dev/tasks/` and commit them together as the task commit.
**Gates** are stop-on-failure. `NOTES_fs-vignette-dgm-predictions_2026-11-05.md` to `~/Downloads` at close.

## Objective

Insert two subsections into `vignettes/articles/mrct-consistency.qmd` (verify path from source): `## The data-generating mechanism {#sec-dgm-overview}` and `## Analytic predictions before the analysis {#sec-predictions}`, verbatim from the insert file, plus one row in the Overview's labels table. Render the article green.

## 1. The insert

- Location: immediately **after** the closing `:::` of the callout-important block "The region in this article" (end of the "The trigger and the question" subsection) and **before** `## The pipeline, for one trial`.
- Content: everything between the `INSERT-BEGIN` and `INSERT-END` markers of `vignette-insert_dgm-predictions_2026-11-05.qmd`, verbatim (markers and the file's leading comment excluded).
- The insert is self-contained: it uses only helpers from the article's `packages` chunk (`rd`, `f3`, `near`), new object names (`trial_pre`, `dgm_pre`, `pop_pre`, `p_pre`, `rows_pre`, `rows_rep`, `w_pre`, `t_pre`, `r_pre`), and shipped records (`rep5_trial.csv`, `predictions_full.csv`). It must not disturb any later chunk; in particular the own-sims section's `own`/`dgm_own` objects are untouched and its build still runs (the article now performs that build twice — once here, once there; record both timings in NOTES, no de-duplication in this task).
- Section ids `#sec-dgm-overview` and `#sec-predictions` must not collide with existing ids (`#sec-dgm` exists and is distinct; verify by grep).

## 2. The labels table

In the Overview's "Labels used throughout." table, add one row after the `joint` row, verbatim:

```
| consistency / non-consistency | the manuscript's names: a positive trial with region ITT $\le$ 1.0; the trigger |
```

No other prose changes anywhere in the article.

## 3. Gates (stop-on-failure)

1. `quarto render` of the article completes; every `stopifnot` in the article passes — including the insert's: the three identities of the predictions rows at 1e-12, and the record check that the report column prints exactly `0.844, 0.326, 0.602, 0.241, 0.187, 0.054, 0.790` (this ties the shipped `predictions_full.csv` to the manuscript's Table 2 S2 column; if it fails, STOP — do not adjust the expected strings).
2. The rep5-based column is **not** gated against any values (differences from the report column are expected); record it in NOTES.
3. `git diff` is confined to `vignettes/articles/mrct-consistency.qmd` (plus the `dev/tasks/` copies). `DESCRIPTION`, `R/`, `NAMESPACE`, records under `inst/extdata/` all untouched.
4. Rendered HTML: the two new subsections appear in the TOC before "The pipeline, for one trial"; all cross-references (`@sec-dgm-overview`, `@sec-predictions`, `@sec-own-sims`, `@sec-case-campaign`, `@sec-sims`, `@sec-dgm`) resolve — no "?@sec-" in the rendered HTML.
5. Article render wall-clock recorded (prior total ≈ 5 minutes per the article's own note; the insert adds one build).

## 4. Close-out

Commit sequence: task commit (task doc + insert in `dev/tasks/`); the edit commit; fast-forward `master`; delete branch. NOTES to `~/Downloads` with: HEADs before/after; the rendered §"Analytic predictions before the analysis" table (both columns, all seven rows) and `r_pre` (the estimated regional shift) and `exp(r_pre)`; the insert build's warning count and seconds, and the own-sims build's seconds (the duplicate-build cost); total render minutes; any flags.

## Out of scope

Any `R/`, records, or DESCRIPTION change; any rewording of existing article text beyond the one labels-table row; harmonizing the article's vocabulary ("trigger", "joint", "oracle joint") to the manuscript's named events beyond the bridge sentences the insert carries — that is a separate decision for Larry; the mrct-consistency repo.
