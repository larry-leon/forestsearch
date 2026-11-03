# TASK: MRCT article — interpretive edits from the review packet (forestsearch, 2026-11-02)

**Repository:** forestsearch. **Branch:** create `article-review-edits` from `master` (at `fa3b1397` or later). **Touches:** `vignettes/articles/mrct-consistency.qmd` and its tracked `.html`. Nothing under `R/`. No chunk computation changes.

## 0. First action

Copy this file to `dev/tasks/TASK_fs-mrct-article-review_2026-11-02.md` and commit it alone ("task: MRCT article review edits").

## 1. Edits (verbatim)

Each edit names text in the current file. Find it by the quoted text (line wrapping may differ); if a quoted text is not found exactly once, stop and say so in the notes.

**E1 — "explain" → "account for" (Overview, "The trigger and the question").**
Replace `can explain the regional result` with `can account for the regional result`.
Replace `and is not evidence of a region-specific effect.` with `and does not by itself point to a region-specific effect.`

**E2 — oracle as a benchmark, not a ceiling.**
In the core-OC reading (Simulation results), replace `against an oracle ceiling of` with `against an oracle benchmark of`.
In callout 4.3d, replace `the oracle joint as the ceiling that removing $H^*$ exactly would give` with `the oracle joint as the benchmark that removing $H^*$ exactly would give under the stated mechanism (a benchmark, not a bound: in a given trial $\hat G$ can do better)`.

**E3 — Sens / Spec labelled as simulation recovery metrics.**
In the core-OC table caption, replace `(Sens / Spec averaged over the declared triggered trials)` with `(Sens / Spec: recovery of the simulation's $H^*$, averaged over the declared triggered trials; not available in a real trial)`.

**E4 — what a move below 1.0 shows, and what it does not.** Insert this callout in the case-study section (`#sec-case`), directly after the paragraph or table that reports region($\hat G$) for replicate 5 (the four-sets result), before the bounds:

```markdown
::: {.callout-warning}
## What the move below 1.0 shows, and what it does not

Removing $\hat H$ from the region and refitting the Cox model compares two different patient sets; it
does not decompose the regional hazard ratio. The region's estimate is not a weighted average of the
estimates in $\hat H$ and in $\hat G$: removing patients changes the risk sets, the events and the
censoring mixture, and each Cox estimate converges to the limit of its own set at the analysis time
(@sec-analytic), which can differ from that set's AHR or CDE. A regional estimate that moves from above
1.0 to below 1.0 is therefore compatible with a trial-wide subgroup accounting for the regional result.
It does not show that $\hat H$ caused it, does not validate $\hat H$ as a clinically actionable
subgroup, and, read with the null cells (@sec-sims), can arise when there is no harmed subgroup at all.
:::
```

(If `@sec-analytic` or `@sec-sims` is not the current label, use the label of the analytic-layer and simulation-results sections and say so in the notes.)

## 2. Gate and notes

- Render by the r3/r4 route (render, then `build_article()`, restore the local HTML). Stop on a render error or a failed check.
- Commit the `.qmd` and re-rendered `.html` together ("article: interpretive edits from review"), then `git checkout master && git merge --ff-only article-review-edits`. Do not push.
- Write `~/Downloads/NOTES_fs-mrct-article-review_2026-11-02.md`: status, commits, new master SHA, where E4 sits, any text not found.
