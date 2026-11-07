# EDITS_fs-article-ms-alignment_2026-11-07 — companion to the task doc

Targets `vignettes/articles/mrct-consistency.qmd`. Matching and application per the task doc's
convention. Everything not inside a block is untouched.

## (A) — §1.1, the two explanations (old → new)

**old**

We consider two explanations for the trigger:

- **a region-specific effect** — the region differs in its treatment effect; or
- **a trial-wide harmed subgroup** — present in every region, with the region's estimate pulled up by it (the region's prognosis and the finite follow-up window shaping how far).

**new**

We consider two explanations for the trigger:

- **a trial-wide harmed subgroup** — the explanation this article's analysis pursues: a subgroup present in every region, with the region's estimate pulled up by it (the region's prognosis and the finite follow-up window shaping how far);
- **a region-specific effect** — the region differs in its treatment effect. The analysis approach does not pursue this explanation; it enters through the extension of @sec-dgm-overview, whose configurations with no regional treatment benefit ($b_R = 1.0$) quantify how often the strategy explains a trigger that a region-specific effect produced — the false-positive rate of the subgroup pathway (the accompanying manuscript's region-interaction cells).

## (B) — §1.1 callout, S2's first occurrence (old → new)

**old**

Replicate 5 of S2, the case study, is one simulated trial. From the point at which its data are read (@sec-case), it is analysed exactly as a real trial would be.

**new**

Replicate 5 of S2 (the simulation cell defined in @sec-dgm-overview), the case study, is one simulated trial. From the point at which its data are read (@sec-case), it is analysed exactly as a real trial would be.

## (C) — §1.4, new callout (insert before the paragraph beginning "$\hat H$ is learned on the non-region patients only")

::: {.callout-tip}
## The pipeline in one call

@sec-case runs the five steps one at a time, each call visible. `mrct_trial_analysis()` runs the same pipeline in one call — @sec-one-call shows it on the case study and checks it against the step-by-step objects.
:::

## (D) — §2.6, delimited section replacement

Keep the line `## The S2 campaign in brief {#sec-case-campaign}`; replace everything after it, up
to (not including) the `## The one-call form` heading line, with:

---BEGIN SECTION CONTENT---

Replicate 5 is one trial of S2. The campaign ran 10,000 trials of the same design (the simulations article). The first table gives the campaign-observed counterparts of the analytic rows of @sec-predictions — unconditional, over all trials; the second follows the pipeline of this section, conditional on the trigger where its denominators say so. The numbers are read from the campaign's records, shipped with the package.

```{r case-campaign}
cc_core <- mrct_record("core_full_numeric.csv")
cc_core <- cc_core[cc_core$setting == "S2_rmlog4", ]
cc_disp <- mrct_record("core_n1000.csv")
stopifnot(nrow(cc_core) == 1, cc_core$replicates == 10000)
# the trial-positive count, from the display record's first row for S2 ("...; N of M)")
r1 <- cc_disp$S2[cc_disp$row == "r1"]
n_A <- as.numeric(gsub(",", "", sub("^.*; ([0-9,]+) of [0-9,]+\\)$", "\\1", r1)))
reps <- cc_core$replicates
n_trig <- cc_core$n_trigger
n_orac <- cc_core$n_Gstar_lt1
n_dcl  <- cc_core$n_declared_Ghat_lt1
stopifnot(length(n_A) == 1, !is.na(n_A), n_trig <= n_A, n_orac <= n_trig, n_dcl <= n_trig)
counts <- c(n_A, n_A - n_trig, n_trig, n_orac, n_trig - n_orac,
            (n_A - n_trig) + n_orac, n_dcl, (n_A - n_trig) + n_dcl)
obs <- counts / reps
# identities by construction: (3)+(4)=(1), (5)+(6)=(4), (3)+(5)=(7), (3)+(8)=(9)
stopifnot(near(obs[2] + obs[3], obs[1], 1e-12), near(obs[4] + obs[5], obs[3], 1e-12),
          near(obs[2] + obs[4], obs[6], 1e-12), near(obs[2] + obs[7], obs[8], 1e-12))
lab_obs <- c("(1) P(ITT ub < 1.0): the trial is positive",
             "(3) P(consistency): positive, region ITT ≤ 1.0",
             "(4) P(non-consistency): the trigger",
             "(5) P(oracle consistency): the oracle joint",
             "(6) P(oracle non-consistency)",
             "(7) P(consistency or oracle consistency)",
             "(8) P̂(declared consistency)",
             "(9) P̂(consistency or declared consistency)")
knitr::kable(data.frame(row = append(lab_obs, "(2) P(region ITT > 1.0)", after = 1),
                        `campaign-observed` = append(sprintf("%s (%s)", f3(obs), big(counts)), "—", after = 1),
                        check.names = FALSE),
             row.names = FALSE, align = "lr",
             caption = "S2 at n = 1,000, campaign-observed over 10,000 trials; rows as in @sec-predictions, identities asserted. Row (2) is analytic-only: the committed core records carry the pipeline from the positive trial onward.")
```

```{r case-campaign-pipeline}
cc_mr <- mrct_record("mr_field_core.csv")
cc_ent <- mrct_record("entry_cells_core_numeric.csv")
cc_ev <- function(q) cc_ent$S2[cc_ent$quantity == q]
stopifnot(sum(cc_mr$row == "lb_100") == 1,
          sum(cc_ent$quantity %in% c("trigger_n", "joint", "Ghat_region")) == 3)
camp <- list(
  reps = reps, n_trig = n_trig, n_decl = cc_core$n_declared, n_joint = n_dcl,
  lb100 = as.numeric(cc_mr$S2[cc_mr$row == "lb_100"]),
  ghat_reg = cc_ev("Ghat_region"))
camp$p_trig <- 100 * camp$n_trig / camp$reps
camp$p_decl <- 100 * camp$n_decl / camp$n_trig
camp$p_joint <- 100 * camp$n_joint / camp$n_trig
# the two records agree on S2
stopifnot(camp$n_trig == cc_ev("trigger_n"), abs(camp$p_joint - cc_ev("joint")) <= 1e-9)
# the reading's verbal claims
stopifnot(round(camp$p_trig / 25) == 1, abs(camp$p_joint - 50) < 5)
pc1 <- function(x) formatC(x, format = "f", digits = 1)
knitr::kable(
  data.frame(
    step = c("trigger: ITT ub < 1.0 and region ITT > 1.0", "subgroup declared",
             "Ĥ field lower bound (one-sided 95%) > 1.0", "declared, region(Ĝ) < 1.0",
             "Ĝ, % of region patients"),
    S2 = c(sprintf("%s (%s)", pc1(camp$p_trig), big(camp$n_trig)),
           sprintf("%s (%s)", formatC(camp$p_decl, format = "f", digits = 2), big(camp$n_decl)),
           pc1(camp$lb100),
           sprintf("%s (%s)", pc1(camp$p_joint), big(camp$n_joint)),
           pc1(camp$ghat_reg)),
    of = c(sprintf("all trials (%s)", big(camp$reps)), "triggered trials", "declared triggered trials",
           "triggered trials", "declared triggered trials")),
  align = "lrl",
  caption = "Conditional on the trigger: the pipeline of this section, as often as it went each way (each row's denominator stated).")
```

Read unconditionally: the trial is positive in `r f3(obs[1])` of trials; consistency holds outright in `r f3(obs[2])`, and the trigger occurs in `r f3(obs[3])` — about one S2 trial in four. The oracle resolution holds in `r f3(obs[4])` of all trials (row 5) and the declared subgroup's in `r f3(obs[7])` (row 8), so the two routes together give `r f3(obs[6])` by the oracle pathway and `r f3(obs[8])` by the declared pathway (rows 7 and 9). Conditional on the trigger: a subgroup was declared in `r big(camp$n_decl)` of `r big(camp$n_trig)` triggered trials, region(Ĝ) was below 1.0 in about half (`r pc1(camp$p_joint)`%), and Ĝ kept on average `r pc1(camp$ghat_reg)`% of the region's patients. The field lower bound for Ĥ was above 1.0 in `r pc1(camp$lb100)`% of the declared triggered trials; in the others it lies below 1.0, and it is read by its location, as for replicate 5 in step 5. Replicate 5 is one of the triggered trials with region(Ĝ) < 1.0.

---END SECTION CONTENT---

(The article's timing sentence "The article took … minutes to run." currently sits in §2.7 and is
untouched. `near`, `f3`, `big` are the setup chunk's helpers.)

## (E1) — §1.3, the events intro and bullets (old → new)

**old**

Seven design probabilities follow from four of them, built from three events:

- $A$ = {ITT ub < 1.0} — the trial is positive;
- $B$ = {region ITT > 1.0} — the region's point estimate is unfavourable;
- $C^*$ = {region($G^*$) < 1.0} — removing $H^*$ exactly leaves the region's estimate below 1.0.

**new**

Seven design probabilities follow from four of them, built from three events (the manuscript's full event notation is reproduced in the [appendix](#sec-event-notation)):

- $A$ = {ITT ub < 1.0} — the trial is positive (the manuscript's $\mathrm{ITT}^+$);
- $B$ = {region ITT > 1.0} — the region's point estimate is unfavourable (the manuscript's $\mathrm{R}^-$);
- $C^*$ = {region($G^*$) < 1.0} — removing $H^*$ exactly leaves the region's estimate below 1.0 (the manuscript's $\mathrm{R}(G^*)^+$).

## (E2) — append at the end of the file

# Appendix: the manuscript's event notation {#sec-event-notation .unnumbered}

The accompanying manuscript collects the design's events in its Table 1, reproduced here; a superscript $+$ ($-$) marks a point estimate or bound in (opposite to) the direction of consistency. This article's $A$, $B$ and $C^*$ are $\mathrm{ITT}^+$, $\mathrm{R}^-$ and $\mathrm{R}(G^*)^+$.

| Event | Definition | Reading |
|---|---|---|
| *Primitive events* | | |
| $\mathrm{ITT}^+$ | {ITT 95% upper bound < 1.0} | The trial is positive |
| $\mathrm{R}^-$, $\mathrm{R}^+$ | {region ITT > 1.0}, {region ITT $\le$ 1.0} | Region point estimate opposite to, in, the direction of consistency |
| $\mathrm{R}(G^*)^-$, $\mathrm{R}(G^*)^+$ | {region($G^*$) > 1.0}, {region($G^*$) $\le$ 1.0} | Likewise, for the region outside the true subgroup $H^*$ |
| $\mathrm{R}(\hat G)^+$ | {subgroup declared, region($\hat G$) < 1.0} | The declared resolution is in the consistent direction |
| *Named events* | | |
| consistency | $\mathrm{ITT}^+ \cap \mathrm{R}^+$ | Positive trial, region in the consistent direction — the standard pathway |
| non-consistency | $\mathrm{ITT}^+ \cap \mathrm{R}^-$ | Positive trial, region opposite — the trigger for subgroup identification |
| oracle consistency | $\mathrm{ITT}^+ \cap \mathrm{R}^- \cap \mathrm{R}(G^*)^+$ | Non-consistency resolved by the true subgroup |
| oracle non-consistency | $\mathrm{ITT}^+ \cap \mathrm{R}^- \cap \mathrm{R}(G^*)^-$ | Non-consistency persisting under the oracle resolution |
| declared consistency | $\mathrm{ITT}^+ \cap \mathrm{R}^- \cap \mathrm{R}(\hat G)^+$ | Non-consistency resolved by the declared subgroup |

: The manuscript's Table 1: event notation.
