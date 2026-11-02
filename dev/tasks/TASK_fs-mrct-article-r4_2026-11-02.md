# TASK: MRCT article r4 — own-trial section with the observed-region mode (forestsearch, 2026-11-02)

**Repository:** forestsearch. **Branch:** `mrct-article-r3` (continue on it). **Touches:** `vignettes/articles/mrct-consistency.qmd` only. Nothing under `R/`.

**Prerequisite:** `master` contains `eaacb9e6` (the observed-region mode, merged from `fs-inject-observed-region`). If not, stop and say so in the notes.

## 0. First actions

1. Check out `mrct-article-r3`; merge `master` into it (`git merge master`); stop on a conflict.
2. Install the package from this tree.
3. Copy this file to `dev/tasks/TASK_fs-mrct-article-r4_2026-11-02.md` and commit it alone ("task: MRCT article r4").

## 1. What changes

Rewrite `#sec-own-sims` ("Setting up simulations from your own trial") so that replicate 5, read from `rep5_trial.csv` and treated as the analyst's trial, is built with `inject_mrct_structure(region = list(observed = ...))`, and the analytic layer runs on that build. Everything else in the article stays as it is.

### Step 1 (evaluated): the DGM from the trial, region observed

- Replace the `generate_aft_dgm_flex()` build with one visible `inject_mrct_structure()` call on the data read from `rep5_trial.csv`: `region = list(observed = <the region column>)`, `loghr` left at NULL (estimated), the same covariates as the report's build, `x_pred = "cd40"`, `spline_spec` with $b = 0.55$ everywhere, `x3` for $H^*$ = CD4 ≤ 230 with $h = 2$, `expand = "copula"`, `n_super = 20000`, a stated `seed`.
- **Column names.** `inject_mrct_structure()` adds the `z_` prefix itself. Read the CSV's column names; if they already carry `z_` (e.g. `z_cd40`, `z_region`), strip it in one visible line in the chunk before the call, so the build's columns are `z_cd40`, not `z_z_cd40`.
- **No `scale_dgm_hazard()` here.** The trial's own data carry its hazard scale (replicate 5 was drawn with all hazards × 3); the AFT fit to it estimates that scale. Say so in one sentence.
- Report the estimated regional log-HR (`dgm$mrct$region$prognostic_loghr`) beside −log 4, the value in the mechanism that generated the trial, and the region's share of the super-population beside its share of the trial.
- Keep the note on any censoring-model warnings if they still occur (all censoring in the trial is administrative); show them (`#| warning: true`) rather than suppressing.

### Step 3 (evaluated): the analytic layer on that build

`ahr_population()` with `analysis_time = 40`, `max_entry = 11`; `ahr_limit()` for the region, the non-region patients and `"region_Gstar"`; `ahr_design()` at $n = 1{,}000$ giving P(trigger) and the oracle joint. One short table. One sentence comparing them with S2's (@sec-analytic), as a reading: the build from one trial of 1,000 is not S2's build from ACTG 175 and is not expected to reproduce it.

### Step 4 (eval: false, as now): the campaign

Use `analysis_covariates(dgm_own, "observed")` for the confounders in place of the explicit twelve-covariate list (the observed-mode record supplies `z_name`, so the region column is excluded).

### Text

Remove the sentences stating the two gaps (no copula expansion for an observed region; `ahr_population()` accepting only planted builds). Add one sentence to Step 1: with `region$observed`, the region is carried through the expansion with the covariates and nothing about membership is modelled. Callouts 4.3a–4.3d unchanged.

## 2. Gate and notes

- Render by the same route as r3 (render, then `build_article()`, restoring the local HTML). Stop on a render error or a failed check.
- Commit on `mrct-article-r3`. Do not push.
- Write `~/Downloads/NOTES_fs-mrct-article-r4_2026-11-02.md`: status, commits, render time, the estimated regional log-HR, the step-3 table, any error text.
