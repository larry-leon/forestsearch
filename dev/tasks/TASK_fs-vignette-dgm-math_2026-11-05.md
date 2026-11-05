# TASK_fs-vignette-dgm-math_2026-11-05

**Dispositions on record (Larry, 2026-11-05 13:25)**: in the mrct-consistency article, (i) the "report (Table 2, S2)" comparison column is removed — the article is a standalone document; (ii) the "The data-generating mechanism" subsection carries the mathematical details of the manuscript's data-generating-mechanism block (manuscript pp. 4–5: the AFT law, the $\theta_0/\psi$ split, (M1)–(M4), and the region-dependent extension), adapted to the article's own notation ($\beta$ for the coefficient vector and $\mathbf 1\{X \in R\}$, as its appendix already uses; no manuscript macros; the proposition reference replaced by a pointer to the article's appendix).

**Repo**: `forestsearch` only. **No `R/` changes of any kind.** One vignette file is edited. **STOP** if anything the replacement uses is missing.
**Branch**: short-lived from `master` (now carrying TASK_fs-vignette-dgm-predictions's merge; record HEAD, verify from source), fast-forward, delete branch (`git branch -r` first).
**Compute**: none beyond the article's own render (unchanged computations; the record-comparison read is removed).
**Session start**: `export GIT_TERMINAL_PROMPT=0`. Credentials-failed fetch recorded, not a stop.
**Files from ~/Downloads**: this task doc and `vignette-replace_dgm-math_2026-11-05.qmd`. **First action**: copy both into `dev/tasks/` and commit them together as the task commit.
**Gates** stop-on-failure. `NOTES_fs-vignette-dgm-math_2026-11-05.md` to `~/Downloads` at close.

## 1. The replacement

In `vignettes/articles/mrct-consistency.qmd`, replace everything from the line
`## The data-generating mechanism {#sec-dgm-overview}` up to (not including) the line
`## The pipeline, for one trial` with the content between the `REPLACE-BEGIN`/`REPLACE-END` markers of the replacement file, verbatim. The `pre-build` chunk inside the replacement is byte-identical to the current one (verify by diff); the changes are:

- §"The data-generating mechanism": the mathematical block (three display equations, the (M1)–(M4) paragraph, the region-dependent extension display with $\gamma_R$, $\gamma_{H^*\!R}$, the base-model null, and the S2 parameter paragraph).
- §"Analytic predictions before the analysis": the `pre-table` chunk loses `rows_rep`, the `predictions_full.csv` read and the record-value gate; the table has one value column ("from rep5_trial"); the caption and prose no longer reference the report's Table 2 or its values; the three identity `stopifnot`s remain.

No other changes anywhere in the article — the labels-table row added by the previous task stays.

## 2. Gates (stop-on-failure)

1. `quarto render` completes; every `stopifnot` in the article passes (the removed record gate is **gone by design** — its absence is correct, do not re-add it).
2. Grep gates on the rendered article source/HTML: zero occurrences of "Table 2" and of "the report (Table 2"; `predictions_full.csv` is no longer read in the Overview section (it is still read later in @sec-own-sims and @sec-analytic — those are untouched).
3. The three display equations and the extension display render (no raw `$$` or LaTeX errors in the HTML); cross-references `@sec-dgm`, `@sec-appendix`, `@sec-own-sims`, `@sec-case-campaign`, `@sec-sims` resolve; no "?@sec-".
4. `git diff` confined to `vignettes/articles/mrct-consistency.qmd` plus the `dev/tasks/` copies; `R/`, `DESCRIPTION`, `NAMESPACE`, `inst/extdata/` untouched.
5. Render wall-clock recorded.

## 3. Close-out

Commit sequence: task commit (both files in `dev/tasks/`); the edit commit; fast-forward `master`; delete branch. NOTES with: HEADs; the rendered one-column predictions table and $\hat r$; confirmation the record gate and column are gone; render minutes; flags.

## Out of scope

Everything else in the article; `R/`; the mrct-consistency repo; vignette-wide vocabulary harmonization.
