# Split specification — TASK_fs-vignette-split_2026-11-05

Companion to the task doc. Part A is new content, verbatim. Parts B–C are the per-document edit
tables. "App" = `mrct-consistency.qmd` (Sections 1–2); "Sims" = the new
`mrct-consistency-simulations.qmd` (Sections 3–6), both under `vignettes/articles/`.

---

## A. New content, verbatim

### A.1 Sims front matter (replaces the YAML of the moved content)

```yaml
---
title: "MRCT regional consistency: simulations"
subtitle: "The data-generating mechanism, the analytic layer and the campaign behind the case study"
author: "ForestSearch Package"
date: today
format:
  html:
    embed-resources: true
    toc: true
    toc-depth: 3
    number-sections: true
    code-fold: true
    code-overflow: wrap
    fig-width: 9
    fig-height: 6
vignette: >
  %\VignetteIndexEntry{MRCT regional consistency: simulations}
  %\VignetteEngine{quarto::html}
  %\VignetteEncoding{UTF-8}
---
```

Then, copied **byte-identical from the app article**: the `setup` chunk (including `nsims_demo`),
the CSS block, and the `packages` chunk (helpers `rd`, `f2`–`f4`, `pc`, `big`, `near`,
`t_article`). The `install-check` chunk is copied with its `needed` vector replaced by the Sims
list of §B.3 below. The "Installing forestsearch" callout is copied unchanged; the "Rendering this
article yourself" callout is copied with the download URL's filename changed to
`mrct-consistency-simulations.qmd`.

### A.2 Sims opening paragraph (after the install/render callouts, before the first section)

> This article is the simulation side of the regional-consistency strategy: how the report's
> data-generating mechanism and cells were built, the analytic layer that predicts the design's
> behaviour before any simulation, the committed campaign's results, and how the case study of the
> application article ([mrct-consistency](mrct-consistency.html)) was found. It executes
> independently of that article: the case study's objects, which the checks below compare against,
> are recomputed first (@sec-case-objects).

### A.3 The bridge section (first section of Sims, before "Setting up simulations from your own trial")

~~~
# The case study's objects, recomputed {#sec-case-objects}

The application article analyses one trial, replicate 5 of S2, from its shipped data
(`rep5_trial.csv`) to its multiplier-resampling bounds, with every call explained. The simulation
sections below check themselves against that analysis — the redraw of @sec-find must reproduce the
trial, its $\hat H$ and its four Cox fits, and @sec-case-sim reads its MR values — so this chunk
recomputes the same objects with the same calls, condensed; the application article is the
walkthrough.

```{r bridge-case}
<CONTENT: the following chunks of the app article, concatenated in order, with their display
 lines removed and their stopifnot checks kept —
 1. `case-read` in full (trial, covs, seed_rep5, and its stopifnot);
 2. from `case-itt`: the fit_itt / itt_ub lines and `stopifnot(itt_ub < 1)` (drop the printed c());
 3. from `case-region`: reg, fit_reg, hr_reg, trigger and `stopifnot(trigger)` (drop the c());
 4. `case-fs` in full (the forestsearch() call with mr_inference and seed_rep5; keep t_mr, mr,
    Hhat_def and the stopifnot; drop the final c() and Hhat_def display lines);
 5. from `case-sets`: in_Hhat and its stopifnot, cox_set(), Hstar, set_idx, four, its near()
    stopifnot, and set_label (drop the kable);
 6. from `mr-table`: mr_tab, mv, cc, fl and the near() stopifnot (drop the printed mr_tab);
 then one closing display line:
 c(patients = nrow(trial), seconds_mr = round(t_mr, 1)); Hhat_def>
```

$\hat H$ is `` `r Hhat_def` ``; the four sets and the bounds are read and discussed in the
application article.
~~~

(The chunk is assembled from the app article's committed code so the two articles cannot drift on
the call; no line of new computation is introduced.)

### A.4 App closing line (end of the app article, after the campaign-in-brief reading)

> The article took `r "`r round((proc.time()[[\"elapsed\"]] - t_article) / 60, 1)`"` minutes to run.

### A.5 App intro addition (end of the opening paragraph "This article reproduces…")

> The simulation design that generated the case study, and the campaign's results, are in the
> companion article [mrct-consistency-simulations](mrct-consistency-simulations.html).

---

## B. The app article (`mrct-consistency.qmd`)

1. **Keep**: front matter; `setup` (drop the `nsims_demo` lines and their comment — nothing in
   §§1–2 uses it); CSS; `packages`; install callout; `install-check` with the App `needed` list of
   §B.3; render callout (URL unchanged); §1 Overview in full; §2 The case study in full, including
   "The S2 campaign in brief". **Remove**: everything from `# Setting up simulations from your own
   trial {#sec-own-sims}` to the end of the file. Append A.4; insert A.5.
2. **Cross-reference rewrites** (each old fragment occurs once unless noted; rewrite exactly):

| Where | Old | New |
|---|---|---|
| §1.1 callout | `planted in a simulated population by a logistic model (@sec-dgm)` | `planted in a simulated population by a logistic model (the simulations article)` |
| §1.1 callout | `is carried into the data-generating mechanism instead (@sec-own-sims)` | `is carried into the data-generating mechanism instead (the simulations article)` |
| §1.2 end | `are in @sec-dgm onward.` | `are in the companion article [mrct-consistency-simulations](mrct-consistency-simulations.html).` |
| §1.3 caption | `are in @sec-case-campaign and @sec-sims.` | `are in @sec-case-campaign and in the simulations article.` |
| §1.3 prose | `against the generating $-\log 4$ = … of @sec-dgm-overview, and its covariates` | keep `@sec-dgm-overview` (internal) — no change |
| §1.3 prose | `a build from the full construction of @sec-dgm gives slightly different ones` | `a build from the full construction (the simulations article) gives slightly different ones` |
| §1.3 prose | `explained in @sec-own-sims` | `explained in the simulations article` |
| §1.3 prose | `(@sec-own-sims walks through them)` | `(the simulations article walks through them)` |
| how-to-use tip | `@sec-own-sims shows how to set up simulations around your trial, treating replicate 5 as if it were yours, with its region observed. @sec-dgm onward is the report's own simulation design, which generated replicate 5, and its results.` | `The companion article [mrct-consistency-simulations](mrct-consistency-simulations.html) shows how to set up simulations around your trial, and contains the report's own simulation design, which generated replicate 5, and its results.` |
| §2 step 3 | `The call is the one `mrct_case_study()` makes (@sec-find), written out` | `The call is the one `mrct_case_study()` makes (see the simulations article), written out` |
| §2 forest text | `from an `mrct_case_study` object only (@sec-find)` | `from an `mrct_case_study` object only (simulations article)` |
| §2 warning callout | `at the analysis time (@sec-analytic)` | `at the analysis time (simulations article)` |
| §2 warning callout | `read with the null cells (@sec-sims)` | `read with the null cells of the simulations article` |
| §2 end of step 5 | `that needs the simulation design, and is in @sec-case-sim.` | `that needs the simulation design, and is in the simulations article.` |
| §2.6 | `The campaign ran 10,000 trials of the same design (@sec-sims)` | `The campaign ran 10,000 trials of the same design (the simulations article)` |

3. **`needed` lists** (`install-check`): App:
   `c("inject_mrct_structure", "ahr_population", "ahr_design", "get_dfpred", "compute_frontier_cis", "mr_estimates_table")`.
   Sims:
   `c("inject_mrct_structure", "scale_dgm_hazard", "ahr_population", "ahr_limit", "ahr_design", "ahr_solve", "mrct_region_sims", "mrct_select_case", "mrct_case_study", "analysis_covariates", "get_dfpred", "mr_estimates_table")`.
   Gate: every name is in `getNamespaceExports("forestsearch")` at render (the chunk stops itself).

---

## C. The sims article (`mrct-consistency-simulations.qmd`)

1. **Assembly order**: A.1 front matter; setup/CSS/packages/install callout/install-check (Sims
   list)/render callout per A.1; A.2 opening paragraph; A.3 bridge section; then, moved from the
   app article **unchanged except the rewrites below**: `# Setting up simulations from your own
   trial {#sec-own-sims}` …; `# How replicate 5 was generated…` (…@sec-dgm, cells, @sec-analytic,
   @sec-roots, @sec-find, @sec-case-sim); `# Simulation results {#sec-sims}` (…@sec-campaign-code,
   @sec-demo); `# Appendix: the analytic tools {#sec-appendix}`.
2. **Cross-reference rewrites**:

| Where | Old | New |
|---|---|---|
| §own-sims step 4 | `the identifier's arguments of @sec-case (shown, not run` | `the identifier's arguments of @sec-case-objects (shown, not run` |
| callout 4.3c | `the identifier's arguments you used in @sec-case` | `the identifier's arguments you used in @sec-case-objects` |
| #sec-find step 3 heading | `**3. Redraw it, and check it against @sec-case.**` | `**3. Redraw it, and check it against @sec-case-objects.**` |
| #sec-find step 3 | `its $\hat H$ and four sets those of @sec-case; the article stops otherwise.` | `its $\hat H$ and four sets those of @sec-case-objects; the article stops otherwise.` |
| #sec-find step 4 | `as in @sec-case` | `as in the application article` |
| #sec-find end | `the report's case example and the trial of @sec-case.` | `the report's case example and the trial of @sec-case-objects.` |
| #sec-case-sim intro | `rather than in @sec-case.` | `rather than in the application article.` |
| #sec-case-sim | `in @sec-rep5-mr` (the MR-draws sentence of #sec-find step 3: `the MR draws are made once, in @sec-rep5-mr`) | `the MR draws are made once, in @sec-case-objects` |
| §sims mr-field intro | `The MR field call of @sec-rep5-mr was run` | `The MR field call of @sec-case-objects was run` |

3. **Rule for any `@sec-` occurrence not in these tables**: after the split, every `@sec-` in each
   article must reference an id defined in the same article. A remaining cross-document reference
   is rewritten by the rule — in App, to "the simulations article" (with the relative link where a
   sentence names the article); in Sims, to `@sec-case-objects` when it denotes the case-study
   analysis/arguments/objects, otherwise to "[the application article](mrct-consistency.html)" —
   and every such extra rewrite is recorded in NOTES.

4. **pkgdown**: if `_pkgdown.yml` lists articles explicitly, add the new article beside the
   existing one; otherwise nothing. Record which.
