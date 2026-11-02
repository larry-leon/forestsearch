# TASK: MRCT article r3 — replicate 5 as the walkthrough; own-trial simulation setup (forestsearch, 2026-11-02)

**Supersedes** `TASK_fs-mrct-article-tone_2026-11-02.md` (delivered, never launched; do not use it).

**Repository:** forestsearch. **Branch:** create `mrct-article-r3` from `master` (at `0b1f5ba7` or later).
**Touches:** `vignettes/articles/mrct-consistency.qmd`; one new data file and its entry in `vignettes/articles/mrct-consistency/PROVENANCE.md`; one generation script under `dev/`. **Nothing under `R/`** (see §5 if the own-trial section turns out to need it).

## 0. First action

Copy this file to `dev/tasks/TASK_fs-mrct-article-r3_2026-11-02.md` and commit it alone ("task: MRCT article r3").

## 1. Goal

The article's main purpose becomes: walk the analyst through the analysis of one trial, replicate 5 of S2, from its data to its bounds, so they can follow the same steps with their own trial. The reader starts from the shipped replicate 5 dataset and needs none of the simulation setup to do so. The simulation material stays, in two places: (i) how an analyst would set up the analogous simulations from a trial with an observed region, using replicate 5 as if it were their trial; (ii) the report's own simulation design (the current DGM build, cells, analytic layer, campaign, case selection), now presented as how replicate 5 was generated.

## 2. Text edits in the opening (verbatim)

**(A)** In the opening paragraph, delete from `It is condensed from the report's chapters` through `the report's derivations are cited, not repeated.` The paragraph then ends at `is applied inside it.`

**(B)** Replace `estimate; no regulator codifies it as a threshold for harm.` with `estimate.`

**(C)** Replace `Two explanations fit a trigger:` with `We consider two explanations for the trigger:`

**(D)** Insert after the paragraph that ends `**signal generation, not confirmation**.` (in the file this wraps across two lines):

```markdown
::: {.callout-important}
## The region in this article

ACTG 175 was not a multiregional trial and has no region variable. The region used here is
hypothetical: about 20% of patients, planted in a simulated population by a logistic model
(@sec-dgm). The way a region enters an analysis differs between real data and simulation:

- **In a real analysis** the region is an observed column of the trial data and nothing about it is
  modelled. The trial is split into the region's patients and the rest; the identifier runs on the
  rest; the declared rule is applied to the region.
- **In simulation** the region has to be generated together with the patients. Here, with no observed
  region, it is planted by a model. With a trial that has an observed region, the observed region is
  carried into the data-generating mechanism instead (@sec-own-sims).

Replicate 5 of S2, the case study, is one simulated trial. From the point at which its data are read
(@sec-case), it is analysed exactly as a real trial would be.
:::
```

## 3. The replicate 5 dataset

Write `vignettes/articles/mrct-consistency/rep5_trial.csv`: replicate 5 of S2 as drawn (the trial `mrct_case_study()` returns as `cs$trial` for `k = 5`, or the equivalent `simulate_from_dgm()` call — read the source for which is exact), one row per patient, with: the patient id, time in months, event indicator, treatment (1 = experimental), the region indicator, every covariate the identifier is given (the `analysis_covariates(dgm_S2, "observed")` set), and `flag_harm` (membership in H*). Keep the column names the analysis uses. Generate it with a script `dev/article/make_rep5_trial.R` (committed; the article does not read `dev/`), and add the file to `PROVENANCE.md` (commit, generating script, row count, the replicate's call seed).

## 4. Article structure after r3

Keep every existing chunk's computation and every check against shipped records; move and re-head sections as below. Section labels in braces are new where marked.

1. **Overview** (`#sec-overview`): opening edited per §2; the trigger and the question; the region note (D); the pipeline for one trial; the labels table. Then the roadmap callout of §4.1.
2. **The case study: replicate 5 of S2** `{#sec-case}` (new; the main walkthrough). See §4.2.
3. **Setting up simulations from your own trial** `{#sec-own-sims}` (new). See §4.3.
4. **How replicate 5 was generated: the report's simulation design** — the current "The data-generating mechanism, built from ACTG 175" (`#sec-dgm`, with its "To use your own trial" list kept), the cells, the analytic layer (`#sec-analytic`, `#sec-roots`) and "Finding a case study" (`#sec-find`). In `#sec-find`, step 3 (`mrct_case_study()`) becomes the check that the simulation route returns the same trial and the same $\hat H$ as §4.2; run it with `mr_inference = FALSE` so the MR draws are made once in the article (in §4.2), and compare $\hat H$ and the four hazard ratios to §4.2's.
5. **Simulation results** (`#sec-sims`) unchanged.
6. **Appendix** unchanged.

Fix every cross-reference the move affects (`@sec-rep5`, `@sec-rep5-mr`, `@sec-find`, `@sec-demo`): the current `# Replicate 5 of S2 {#sec-rep5}` material is absorbed into §4.2; keep `#sec-rep5-mr` as the label of the inference subsection there.

### 4.1 Roadmap callout (verbatim; at the end of the Overview)

```markdown
::: {.callout-tip}
## How to use this article

@sec-case is the analysis of one trial, replicate 5 of S2, read from a shipped dataset and taken from
the first Cox fit to the bounds on the declared subgroup. With your own trial, follow the same steps
with your data frame and your column names; nothing in that section depends on how the trial arose.
@sec-own-sims shows how to set up simulations around your trial, treating replicate 5 as if it were
yours, with its region observed. @sec-dgm onward is the report's own simulation design, which
generated replicate 5, and its results.
:::
```

### 4.2 The case study section (`#sec-case`)

Start by reading `rep5_trial.csv` (and nothing from the simulation sections). Then the pipeline's five steps, each as a visible call of the fundamental function — not a wrapper:

1. ITT Cox fit (`survival::coxph`), ITT ub.
2. Region Cox fit: the trigger.
3. `forestsearch()` on the non-region patients with the region's patients as the set the declared rule is applied to, with the campaign's arguments **and** the MR arguments, in one call. Read the exact argument names and values from `mrct_case_study()`'s source (how it calls `forestsearch()`), and write the call out in full. Same seed (`seed_S2 + 5`) under R's default generator; the replicate's call seed is written as a literal or read from the shipped record, not from the simulation sections.
4. Apply $\hat H$ to the region; the Cox fit on $\hat G$; the four sets table and the forest plot and Kaplan–Meier panels (use the package's plotting functions on these objects if `plot(cs)` needs a `mrct_case_study` object; if no exported function draws them from these objects, draw them with the same code `plot.mrct_case_study` uses, kept short).
5. The Pareto frontier, `mr_estimates_table()`, and the reading of the bounds.

Each step gets one line of plain direction, of the form "With your trial: …" (which column or argument to change). Place the bounds-reading callout below after step 5 (verbatim):

```markdown
::: {.callout-note}
## Reading the bounds

`mr_estimates_table()` needs only the `forestsearch()` result run with `mr_inference = TRUE`. Read the
field lower bound for $\hat H$ by its location against the hazard ratios that would matter clinically
in your setting: a lower bound well below 1.0 leaves benefit in $\hat H$ open; one near or above 0.90
says that, after accounting for the selection, little benefit remains there. Read the field upper
bound for $\hat H^c$ the same way: below 1.0 says the non-region patients outside $\hat H$ benefit.
Use the Bonferroni pair when both statements are made together. The naive estimate in $\hat H$ is
selected and is not a bound. The re-selection share $\hat p(\hat H)$ is a stability diagnostic: a
small share says another rule would often have been declared from similar data.

These bounds are for sets among the non-region patients. The region's fits on the whole region and
on $\hat G$ are ordinary Wald fits on sets defined by a fixed rule, with the region's few events; read
their intervals by location as well.
:::
```

Then a subsection **"Because this trial is simulated"**: $\hat H$ against $H^*$ (the 2 × 2 tables, from `flag_harm`), region($G^*$), and the four sets beside their analytic limits (this subsection may use `dgm_S2` and the analytic layer, so it moves after, or refers forward to, `#sec-dgm`; say so in one sentence). Open it with: `With a real trial none of this is available: $H^*$ is unknown, so there is no $G^*$, no oracle and no sensitivity or specificity.`

Keep the existing reproduction checks (four hazard ratios, counts, $\hat H$, 2 × 2 tables, MR values against `mr_field_S2_rep5.csv`) here, against §4.2's own results; the article still stops on a mismatch.

### 4.3 Setting up simulations from your own trial (`#sec-own-sims`)

Treat `rep5_trial.csv` as the analyst's trial, with its region **observed**. Outline, in this order, each step as a call where the package supports it:

1. Build a DGM from the trial in which the observed region is carried as a patient characteristic (with its prognostic effect estimated from the data), rather than planted by a membership model.
2. State $H^*$ and the cells: $h$ inside $H^*$, $b$ outside, null cells included (callout 4.3a).
3. The analytic layer on that DGM: P(trigger) at the trial size, the four truths of the sets (callout 4.3b).
4. The campaign: `mrct_region_sims(record = "full")` with the trial's region variable (callout 4.3c).
5. Case selection and the summaries conditional on the trigger (callout 4.3d).

**First read the source** of `inject_mrct_structure()`, `create_dgm_for_mrct()`, `generate_aft_dgm_flex()`, `ahr_population()` and `mrct_region_sims()` and establish whether step 1 can be done with existing exported functions (e.g. the region indicator passed as a covariate of the AFT fit and of the copula expansion, with no membership model). Then:

- **If it can:** write the steps as calls on `rep5_trial.csv`, evaluated, with one cell (S2's $h = 2$, $b = 0.55$, $H^*$ = CD4 ≤ 230) and the analytic layer's P(trigger) at $n = 1{,}000$; the campaign call is shown with `#| eval: false` (no campaign is run here).
- **If it cannot:** do not change anything under `R/`. Write step 1 as prose naming exactly what is missing (which function, which argument), steps 2–5 as `#| eval: false` calls with the gap marked, and record the finding in the notes file (§6). A proposal for the missing piece comes back to chat; it is not part of this task.

Callouts (verbatim):

```markdown
::: {.callout-note}
## 4.3a — Choosing the cells

A cell is one data-generating mechanism: your trial's build with a stated $H^*$, a hazard ratio $h$
inside it and $b$ outside it, in every region. Choose $h$ as the sizes of harm you want the strategy
to be able to explain, and keep cells in which there is nothing harmful to find ($h = b$; $h$ between
$b$ and 1; $h = 1$): those are where the strategy explains a trigger that has no trial-wide
explanation, and the rate at which it does so is a result to report.
:::
```

```markdown
::: {.callout-note}
## 4.3b — The analytic layer

`ahr_population()` on each cell's mechanism, with your analysis time and entry period, then
`ahr_limit()` on the sets you care about gives each set's four truths: read margHR at the analysis
time as what the set's Cox estimate converges to, and the AHR and CDE as the effect in the set.
`ahr_design()` at your trial size gives P(ITT ub < 1.0), P(region ITT > 1.0), P(trigger) and the
oracle joint before any simulation; if P(trigger) is very small, the campaign will have few triggered
trials to summarise. `ahr_solve()` places one free parameter so that P(trigger) reaches a target
(@sec-roots). These are computed by quadrature and take seconds per cell.
:::
```

```markdown
::: {.callout-note}
## 4.3c — The campaign

Give `mrct_region_sims()` your cell's mechanism, your trial size, your region variable and the
identifier's arguments you used in @sec-case; `sim_args` and `analysis_time` carry the entry, window
and dropout of the design. Keep `record = "full"`: case selection reads the per-trial Cox fits it
stores. Give each cell its own seed and record it: a replicate is redrawn from the mechanism, the call
seed and its index. Before a full campaign, time a short run at the worker count you will use; the
time per trial on a lightly loaded machine does not carry over to all workers busy.
:::
```

```markdown
::: {.callout-note}
## 4.3d — Reading the results

Summarise the campaign conditional on the trigger: a team runs the joint analysis only when the trial
is positive and the region's estimate is above 1.0, so the triggered trials are the denominator that
matters. Report for each cell the share of triggered trials with a declaration, the joint and, beside
it, the oracle joint as the ceiling that removing $H^*$ exactly would give; then the same rows in the
null cells. `mrct_select_case()` picks a case by a stated rule before it is looked at, and
`mrct_case_study()` redraws and examines it (@sec-find).
:::
```

## 5. What stays fixed

Every number, table and check in the current simulation sections and results; the collapsed "Reproducing this report's build exactly" chunk; the appendix. Terminology as in the article now (trigger, joint, Ĥᶜ, "built from ACTG 175").

## 6. Gate and notes

- Render the article by the same route as `dev/tasks/TASK_fs-mrct-article-r2_2026-11-01.md`. Stop on a render error or a failed reproduction check.
- Commits on `mrct-article-r3` (the data file and script may be their own commit). Do not push.
- Write `~/Downloads/NOTES_fs-mrct-article-r3_2026-11-02.md`: status, commits, render time, the `forestsearch()` call as written in §4.2, the finding of §4.3 step 1 (which function carries an observed region, or what is missing), any error text.
