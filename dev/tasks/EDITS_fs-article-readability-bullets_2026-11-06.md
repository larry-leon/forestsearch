# EDITS_fs-article-readability-bullets_2026-11-06 — companion to the task doc

Seven old → new replacements in `vignettes/articles/mrct-consistency.qmd`. Match each `old` with
whitespace-normalized exact matching (one occurrence); write `new` with the list structure as
given, prose re-wrapped at ~100 characters. Sentences not inside a block are untouched.

## R1 — §1.1 "The trigger and the question", the two explanations

**old**

We consider two explanations for the trigger: the region differs in its treatment effect, or the trial contains a subgroup that is harmed, present in every region, and the region's estimate is pulled up by it (with the region's prognosis and the finite follow-up window shaping how far).

**new**

We consider two explanations for the trigger:

- **a region-specific effect** — the region differs in its treatment effect; or
- **a trial-wide harmed subgroup** — present in every region, with the region's estimate pulled up by it (the region's prognosis and the finite follow-up window shaping how far).

(The following "The question is whether …" sentences are untouched.)

## R2 — §1.2, the S2 cell's specification

**old**

In S2, the cell that produced replicate 5: $b = 0.55$, $h = 2$, $r = -\log 4$, the region about 20% of patients (planted by a logistic model with odds ratio 3 on symptomatic status and no CD4 term), all hazards multiplied by 3, and a trial is $n = 1{,}000$ patients, 1:1, entry uniform over 11 months, analysed at 40 months, with administrative censoring only.

**new**

In S2, the cell that produced replicate 5:

- **effects** — $b = 0.55$ outside $H^*$, $h = 2$ inside it;
- **region** — about 20% of patients, planted by a logistic model with odds ratio 3 on symptomatic status and no CD4 term, with prognostic shift $r = -\log 4$;
- **hazard scale** — all hazards multiplied by 3;
- **the trial** — $n = 1{,}000$ patients, 1:1, entry uniform over 11 months, analysed at 40 months, with administrative censoring only.

(The following companion-link sentence is untouched.)

## R3 — §1.3, the three calls

**old**

Before a single Cox model is fitted, the design's event probabilities can be computed with no simulation. Three calls do it: rebuild the mechanism **from `rep5_trial.csv`**, with its region taken as observed and S2's stated effects ($h = 2$ in $H^*$ = {CD4 $\le$ 230}, $b = 0.55$ outside; the region's prognostic shift $r$ is estimated from the trial); `ahr_population()` at the design's window; `ahr_design()` at the trial size.

**new**

Before a single Cox model is fitted, the design's event probabilities can be computed with no simulation. Three calls do it:

1. rebuild the mechanism **from `rep5_trial.csv`**, with its region taken as observed and S2's stated effects ($h = 2$ in $H^*$ = {CD4 $\le$ 230}, $b = 0.55$ outside; the region's prognostic shift $r$ is estimated from the trial);
2. `ahr_population()` at the design's window;
3. `ahr_design()` at the trial size.

## R4 — §1.3, the seven design probabilities

**old**

Seven design probabilities follow from four of them: with $A$ = {ITT ub < 1.0}, $B$ = {region ITT > 1.0} and $C^*$ = {region($G^*$) < 1.0}, consistency is $A \setminus B$, non-consistency (the trigger) is $A \cap B$, oracle consistency (the oracle joint) is $A \cap B \cap C^*$, and row (7) adds the two routes to consistency.

**new**

Seven design probabilities follow from four of them, built from three events:

- $A$ = {ITT ub < 1.0} — the trial is positive;
- $B$ = {region ITT > 1.0} — the region's point estimate is unfavourable;
- $C^*$ = {region($G^*$) < 1.0} — removing $H^*$ exactly leaves the region's estimate below 1.0.

The named rows:

- **consistency** = $A \setminus B$;
- **non-consistency (the trigger)** = $A \cap B$;
- **oracle consistency (the oracle joint)** = $A \cap B \cap C^*$;
- **row (7)** adds the two routes to consistency.

## R5 — §1.3, the paragraph after the predictions table

**old**

This build estimates the region's prognostic shift from the trial: $\hat r$ = `r f3(r_pre)` ($e^{\hat r}$ = `r f3(exp(r_pre))`), against the generating $-\log 4$ = `r f3(-log(4))` of @sec-dgm-overview, and its covariates are one trial of 1,000 re-expanded by the copula — so these values carry the sampling of that one trial, and a build from the full construction (the simulations article) gives slightly different ones. Consistency, non-consistency and oracle consistency are the names used in the accompanying manuscript; this article's trigger is the non-consistency event and its oracle joint the oracle consistency, as the rows state. Any build warnings are the censoring model's fit and are explained in the simulations article. *With your trial:* the same three calls, with your columns and your stated $H^*$, $b$ and $h$ (the simulations article walks through them); if P(non-consistency) is very small at your design, a campaign will have few triggered trials to summarise.

**new**

This build estimates the region's prognostic shift from the trial: $\hat r$ = `r f3(r_pre)` ($e^{\hat r}$ = `r f3(exp(r_pre))`), against the generating $-\log 4$ = `r f3(-log(4))` of @sec-dgm-overview, and its covariates are one trial of 1,000 re-expanded by the copula — so these values carry the sampling of that one trial, and a build from the full construction (the simulations article) gives slightly different ones. Any build warnings are the censoring model's fit and are explained in the simulations article.

Consistency, non-consistency and oracle consistency are the names used in the accompanying manuscript; this article's trigger is the non-consistency event and its oracle joint the oracle consistency, as the rows state.

*With your trial:* the same three calls, with your columns and your stated $H^*$, $b$ and $h$ (the simulations article walks through them); if P(non-consistency) is very small at your design, a campaign will have few triggered trials to summarise.

(All sentences unchanged; the warnings sentence moves beside the build reading; the inline-R
expressions are identical.)

## R6 — §2, the case-study columns

**old**

The columns are the ones the analysis uses: the patient `id`, the time `y_sim` in months, the event indicator `event_sim`, the treatment `treat_sim` (1 = experimental), the region indicator `z_region`, the twelve baseline covariates the identifier is given (`z_age` to `z_symptom`, on their natural scales), and `flag_harm`, membership of $H^*$. A real trial has no `flag_harm`; here it is used only for the region($G^*$) row of step 4 and in the simulations article.

**new**

The columns are the ones the analysis uses:

- `id` — the patient;
- `y_sim`, `event_sim` — the time in months, and the event indicator;
- `treat_sim` — the treatment (1 = experimental);
- `z_region` — the region indicator;
- `z_age` to `z_symptom` — the twelve baseline covariates the identifier is given, on their natural scales;
- `flag_harm` — membership of $H^*$. A real trial has no `flag_harm`; here it is used only for the region($G^*$) row of step 4 and in the simulations article.

## R7 — §2, step 3: the identifier's settings

**old**

`forestsearch()` runs on the non-region patients, with the region's patients passed as the test set (`df.test`) to which the declared rule is applied, with the campaign's arguments: candidate cuts on the observed covariates with deciles for CD4, weight and prior ART, and at most two factors per rule. Two thresholds do different work: the screen admits a candidate whose full-data hazard ratio is at least 0.90 (`hr.threshold`), and the consistency requirement then asks the candidate to show a hazard ratio of at least 0.80 in both halves of a random split of the data, in at least 90% of the splits (`hr.consistency`, `pconsistency.threshold`). A candidate must also meet the identifier's subgroup-size and per-arm event minimums, at their defaults here. Among the candidates that pass, the declared subgroup is the largest within an effect neighbourhood of 0.20 of the largest effect.

**new**

`forestsearch()` runs on the non-region patients, with the region's patients passed as the test set (`df.test`) to which the declared rule is applied, with the campaign's arguments:

- **candidates** — cuts on the observed covariates, with deciles for CD4, weight and prior ART, and at most two factors per rule;
- **the screen** — a candidate is admitted when its full-data hazard ratio is at least 0.90 (`hr.threshold`);
- **the consistency requirement** — the candidate must show a hazard ratio of at least 0.80 in both halves of a random split of the data, in at least 90% of the splits (`hr.consistency`, `pconsistency.threshold`);
- **admissibility** — the identifier's subgroup-size and per-arm event minimums, at their defaults here;
- **selection** — among the candidates that pass, the declared subgroup is the largest within an effect neighbourhood of 0.20 of the largest effect.

(The following "Multiplier resampling (MR, …)" sentence and the *With your trial* aside are
untouched.)
