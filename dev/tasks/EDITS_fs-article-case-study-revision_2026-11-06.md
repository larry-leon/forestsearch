# EDITS_fs-article-case-study-revision_2026-11-06 — companion to the task doc

Five old → new replacements in `vignettes/articles/mrct-consistency.qmd`. Match `old` with
whitespace-normalized exact matching (one occurrence); write `new` re-wrapped at ~100 characters.
Sentences shown around the blocks (the MR sentence after E3; the *With your trial* asides) are
context only where not inside a block. The anchors were taken from the pre-wrapper source of
2026-11-06; any mismatch is a G1 STOP.

## E1 — §1 Overview, "The pipeline, for one trial", closing paragraph

**old**

$\hat H$ is learned on the non-region patients only and applied to the region, whose outcomes the identifier never sees: the region is a held-out test set, and its Cox fits on the whole region and on $\hat G$ are ordinary Wald fits on sets defined by a fixed rule. The non-region estimate inside $\hat H$ is selected; @sec-rep5-mr bounds it by multiplier resampling.

**new**

$\hat H$ is learned on the non-region patients only and applied to the region, whose outcomes the identifier never sees: the region is a held-out test set, and its Cox fits on the whole region and on $\hat G$ are ordinary Wald fits on sets defined by a fixed rule. The hold-out removes one selection, not two: the region's fits carry no bias from a search their outcomes never entered, but the pipeline itself runs only when the region's estimate is already above 1.0, so everything after step 2 is conditioned on an unfavourable regional draw. That conditioning is why the strategy stays signal generation: the null cells of the simulations article show how often it "explains" a trigger with nothing to find. The non-region estimate inside $\hat H$ is selected; @sec-rep5-mr bounds it by multiplier resampling.

## E2 — Step 2 prose (the region's Cox fit, and the trigger)

**old**

The region's estimate is `r f3(hr_reg)`, above 1.0, from `r sum(trial$event_sim[reg])` events: the trigger. *With your trial:* your region column in `reg`; if region ITT $\le$ 1.0 the pipeline stops here.

**new**

The region's estimate is `r f3(hr_reg)` (95% CI `r f3(fit_reg$conf.int[1, "lower .95"])`, `r f3(fit_reg$conf.int[1, "upper .95"])`), above 1.0, from `r sum(trial$event_sim[reg])` events: the trigger. The interval is wide, and a point estimate above 1.0 is not evidence that treatment harms the region's patients; it is the second half of the trigger, and what it triggers is an investigation, not a conclusion. *With your trial:* your region column in `reg`; if region ITT $\le$ 1.0 the pipeline stops here.

(If the wrapper migration reworded the *With your trial* aside, the full block will miss. Fall back
to matching the `old` text up through "events: the trigger." only, apply the `new` text up through
"not a conclusion.", and leave the aside as it stands; STOP only if that reduced match also fails
to occur exactly once.)

## E3 — Step 3 opening prose (the identifier's settings)

**old**

`forestsearch()` runs on the non-region patients, with the region's patients passed as the test set (`df.test`) to which the declared rule is applied, with the campaign's arguments: candidate cuts on the observed covariates with deciles for CD4, weight and prior ART, at most two factors, hazard-ratio threshold 0.90 for screening and 0.80 for consistency at probability 0.90, the largest consistent subgroup within an effect neighbourhood of 0.20 of the largest effect.

**new**

`forestsearch()` runs on the non-region patients, with the region's patients passed as the test set (`df.test`) to which the declared rule is applied, with the campaign's arguments: candidate cuts on the observed covariates with deciles for CD4, weight and prior ART, and at most two factors per rule. Two thresholds do different work: the screen admits a candidate whose full-data hazard ratio is at least 0.90 (`hr.threshold`), and the consistency requirement then asks the candidate to show a hazard ratio of at least 0.80 in both halves of a random split of the data, in at least 90% of the splits (`hr.consistency`, `pconsistency.threshold`). A candidate must also meet the identifier's subgroup-size and per-arm event minimums, at their defaults here. Among the candidates that pass, the declared subgroup is the largest within an effect neighbourhood of 0.20 of the largest effect.

(The following "Multiplier resampling (MR, …)" sentence is untouched. The minimums sentence is
conditional on G4's documentation check.)

## E4 — Step 4 callout "What the move below 1.0 shows, and what it does not", closing sentences

**old**

A regional estimate that moves from above 1.0 to below 1.0 is therefore compatible with a trial-wide subgroup accounting for the regional result. It does not show that $\hat H$ caused it, does not validate $\hat H$ as a clinically actionable subgroup, and, read with the null cells of the simulations article, can arise when there is no harmed subgroup at all.

**new**

A regional estimate that moves from above 1.0 to below 1.0 is therefore compatible with a trial-wide subgroup accounting for the regional result. It does not show that $\hat H$ caused it, does not validate $\hat H$ as a clinically actionable subgroup, and, read with the null cells of the simulations article, can arise when there is no harmed subgroup at all. Nor is the distance between region($\hat G$) and region($G^*$) a measure of how well $\hat H$ recovered $H^*$: they are different patient sets, and each estimate converges to its own set's limit (the simulations article sets $\hat H$ against $H^*$ directly).

## E5a — Step 5 prose (reading the bounds)

**old**

Read by location: the naive hazard ratio of $\hat H$ is `r f3(mr$naive$est)` and its de-biased value `r f3(mr$debiased$est)`, but the field lower bound is `r f3(fl$lower_1s)`, `r if (fl$lower_1s < 0.90) "below 0.90" else "above 0.90 but below 1.0"`: once the selection is accounted for, the data do not rule out a modest benefit in $\hat H$. The field upper bound for $\hat H^c$ is `r f3(fl$complement$upper_1s)`, below 1.0: the non-region patients outside $\hat H$ benefit. $\hat H$ is re-selected on a share `r f3(mr$reselection$p_hat[[mr$selected_index]])` of the draws.

**new**

Read by location: the naive hazard ratio of $\hat H$ is `r f3(mr$naive$est)` and its de-biased value `r f3(mr$debiased$est)`, but the field lower bound is `r f3(fl$lower_1s)`, `r if (fl$lower_1s < 0.90) "below 0.90" else "above 0.90 but below 1.0"`: once the selection is accounted for, the data do not rule out a modest benefit in $\hat H$ — and a lower bound does not establish harm there either. The field upper bound for $\hat H^c$ is `r f3(fl$complement$upper_1s)`, below 1.0: the non-region patients outside $\hat H$ benefit. The two bounds are one-sided statements about two different sets, not an interval for the contrast between them, so together they do not establish that the two effects are separated; the Bonferroni pair in the table is the joint statement. $\hat H$ is re-selected on a share `r f3(mr$reselection$p_hat[[mr$selected_index]])` of the draws.

## E5b — Step 5 callout "Reading the bounds", two sentences

**old**

The naive estimate in $\hat H$ is selected and is not a bound. The re-selection share $\hat p(\hat H)$ is a stability diagnostic: a small share says another rule would often have been declared from similar data.

**new**

The naive estimate in $\hat H$ is selected and is not a bound. The field bounds are not ordinary intervals computed after treating $\hat H$ as fixed: they are selection-adjusted by construction, bounds over the multiplier draws of the whole candidate family. The re-selection share $\hat p(\hat H)$ is a stability diagnostic for the exact declared rule: a small share says another rule would often have been declared from similar data. It does not say whether neighbouring rules or the broader patient pattern are stable, and it is not evidence for or against a harmed subgroup existing.
