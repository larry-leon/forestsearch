# Interpretable subgroup identification with forestsearch

``` r

knitr::opts_chunk$set(comment = "#>", fig.align = "center")
library(forestsearch)
library(survival)
```

A confirmatory trial is powered for an effect in its whole population,
and the indication it supports describes a population in which efficacy
and a favourable benefit–risk balance have been concluded. Regulatory
guidance asks that this conclusion be examined across relevant subgroups
rather than assumed, and that a sub-population in which it does not hold
be left out of the indication; it credits a subgroup finding only when
credible — tied to prior plausibility and replication — and warns that
an estimate singled out for its extremity is biased by the selection
(European Medicines Agency 2019). Approvals have been restricted to a
subgroup after a positive result in the whole population (Amatya et
al. 2021).

A subgroup of concern is reached by searching the trial’s covariates
after the data are seen, or by examining many prespecified subgroups and
singling one out. Either way it is chosen for the extremity of its
estimate, and its effect is reported through the analysis the trial
would have run had that subgroup been the only one: the treatment
coefficient of the trial’s working model fitted within it, a Cox log
hazard ratio for a time-to-event outcome. So chosen, the subgroup
carries an estimate that overstates its effect, and the same selection
distorts the estimate for the patients outside it. **forestsearch** is
built for this situation.

**Identification — what the package searches for, and how**

- **Interpretable candidates.** Every candidate subgroup is a
  conjunction of at most two covariate cuts, so whether a patient
  belongs is read from the patient’s own measurements.
- **Three identifiers, one form.** Forest search (León et al. 2024), the
  causal forest of generalized random forests (GRF; Athey et al. 2019,
  Cui et al. 2023) and the difference-in-natural-parameters estimator
  (DINA; Gao and Hastie 2025) each score a family of candidates by an
  estimated treatment effect, admit those whose standard-analysis
  estimate $`\hat\beta(g)`$ clears a threshold common to all three, and
  return one admitted candidate as a selection rule directs — or none,
  if none is admitted. Forest search enumerates the candidates
  exhaustively; DINA and GRF propose them from a surface fitted to the
  outcome, and every proposal is re-ranked on $`\hat\beta(g)`$ against
  the same threshold.
- **The orientation.** The search runs in the direction of a sub-optimal
  effect — the package’s “harm” direction, a hazard ratio above one for
  a survival outcome — the direction in which the trial’s conclusion of
  benefit would not hold for the patients concerned.
- **The consistency screen, in closed form.** Forest search additionally
  requires a candidate’s effect to replicate across random half-splits.
  That screen is complementary-pairs stability selection applied to a
  smooth estimating-function coefficient, so the replication rate
  follows in closed form from one fit of the working model within the
  candidate: nothing is drawn, and the search is deterministic given the
  data.
- **The selection rule is a modelling choice.** The rule moves the
  returned region along a frontier of size against effect — the largest
  admitted region, the largest within a band of the strongest effect,
  the strongest effect itself — and so determines what identification
  and inference are about.
- **The targets.** The region returned is $`\hat H`$; the targets are
  $`\beta(\hat H)`$ and, for the patients outside it,
  $`\beta(\hat H^{c})`$: the working model’s coefficient in each
  population, its least-false value whether or not the model is correct,
  on the scale on which the trial was powered and on which regulators
  and clinicians read its results.

**Inference — what is reported about the identified region, and about
the rest of the trial**

- **The correction.** Refit-free multiplier resampling over the
  candidate family removes the optimism of selection from the estimate
  in $`\hat H`$ and in $`\hat H^{c}`$ (León and Anderson 2026a).
- **One-sided bounds, and a claim about both groups.** The *field*
  inverts the law of the residual the correction leaves into a one-sided
  lower bound on $`\beta(\hat H)`$; *field-s* gives the one-sided upper
  bound on $`\beta(\hat H^{c})`$, relevant to a benefit claim for the
  patients outside the subgroup; the Bonferroni pair joins the two in a
  two-subgroup claim at a joint level. Coverage is asymptotically exact
  under separated selection, characterized under local competition, and
  conditional on the family proposed when the family comes from the
  outcomes (León and Anderson 2026b).
- **Declaration control.** Whether a subgroup should be declared at all
  is governed separately: the law of the maximum over the enumerated
  family gives the family-wise rate of falsely declaring a subgroup when
  every candidate’s effect is at least a clinically specified benefit —
  a diagnostic of the rate the chosen thresholds imply and, where a
  capped rate is wanted, the calibrated cutoff that holds it (León and
  Anderson 2026b).
- **Reading.** Each estimate and bound is read by its location against
  the trial’s own effect and against clinically meaningful effect sizes,
  not as a test at the null. Whether and how to act on a finding is a
  broader judgement of benefit and risk, which this inference is meant
  to inform.

**This vignette** runs the analysis reported as Table 1 of León and
Anderson (2026b), on the German Breast Cancer Study Group (GBSG) trial
that ships with **survival**, at the manuscript’s settings: the three
identifiers, each by one call; the selection-adjusted analysis of each
identified region and its complement; and the declaration calibration at
a protected level. The computation takes about a minute on one core; the
fit times are printed as they are measured.

## The data

The analysis set is the 686 patients with complete prognostic data, out
of 720 randomized: 246 received hormonal therapy and 440 did not, with
299 recurrences or deaths. Seven covariates enter the search: age,
menopausal status, tumour size, grade, the number of positive lymph
nodes, and the progesterone and estrogen receptor concentrations in
fmol/l. Randomization status is masked in the distributed data; hormonal
therapy is analysed as randomized.

``` r

df <- survival::gbsg
df <- within(df, {
  id          <- seq_len(nrow(df))              # the fit's subject identifier
  time_months <- rfstime / 30.4375              # recurrence-free survival, days to months
  grade3      <- ifelse(grade == "3", 1, 0)     # grade 3 against grade 1 or 2
})
confounders <- c("age", "meno", "size", "grade3", "nodes", "pgr", "er")
stopifnot(nrow(df) == 686, sum(df$status) == 299)

whole <- coxph(Surv(time_months, status) ~ hormon, data = df)
whole_hr <- c(est = unname(exp(coef(whole))), exp(confint(whole))[1, ])
c(patients = nrow(df), hormonal_therapy = sum(df$hormon == 1),
  none = sum(df$hormon == 0), events = sum(df$status))
```

    #>         patients hormonal_therapy             none           events
    #>              686              246              440              299

``` r

sprintf("Whole-trial hazard ratio %.3f (%.3f, %.3f), hormonal therapy vs none",
        whole_hr[["est"]], whole_hr[[2]], whole_hr[[3]])
```

    #> [1] "Whole-trial hazard ratio 0.695 (0.544, 0.888), hormonal therapy vs none"

Every bound below is read against the trial’s own effect, a hazard ratio
of 0.695, and against clinically meaningful effect sizes — not as a test
at a hazard ratio of one.

## The specification

Candidates are conjunctions of at most two covariate cuts (`maxk = 2`)
with more than 60 patients (`n.min`); forest search additionally
requires ten events in each arm (`d0.min`, `d1.min`). Continuous
covariates take the package’s default cuts, except that the two
receptors carry ten candidate cuts each (`conf.cont_jcuts`) with the
value zero — no measured receptor — forced in as a cut of its own
(`conf_force`).

Three thresholds govern the search, on the hazard-ratio scale with harm
upward:

- the **screening threshold** $`c_1`$ = `hr.threshold`, which a
  candidate’s estimated hazard ratio must clear to be admitted;
- the **per-split threshold** $`c_2`$ = `hr.consistency`, the claim
  threshold at which an effect stops being beneficial, which both halves
  of a random half-split must clear for the split to count as
  consistent;
- the **consistency level** $`p^\star`$ = `pconsistency.threshold`, the
  proportion of consistent splits a candidate must reach.

The applications of León and Anderson (2026b) set $`c_1 = c_2 = 1.0`$
and $`p^\star = 0.90`$. The splitting-consistency rate is *not* computed
by drawing splits. With `consistency_method = "resample"` (the default)
the rate is evaluated in closed form, as Section 2.2 of that paper gives
it: to first order the two half-sample estimates are a mirror pair about
the full-sample fit, so the probability that both clear $`c_2`$ is
``` math
p_{\mathrm{cons}}(g) \;=\; 2\,\Phi\!\left(\frac{\hat\beta(g) - \log c_2}{\sigma_D(g)}\right) - 1,
```
with $`\sigma_D(g)`$ the scale of the half-sample difference, obtained
from the fit’s influence function. Despite its name, `"resample"` draws
nothing. Literal random half-splitting (`consistency_method = "split"`,
`fs.splits` refits per candidate) is retained only as the fallback for a
candidate whose closed form cannot be evaluated — a non-convergent fit,
for instance — and `fs.splits` has no effect otherwise. Here every
candidate has at least ten events in each arm and the Cox fits converge,
so the closed form is evaluable throughout and the search is
deterministic given the data. The rate is rounded to
`pconsistency.digits = 2` decimals before it is compared with
$`p^\star`$.

The **selection rule** is $`\max_{\mathcal{A}} N(\varepsilon)`$: among
the admitted candidates $`\mathcal{A}`$, the largest whose hazard ratio
is within a fraction $`\varepsilon = 0.20`$ of the largest hazard ratio
over the family (`sg_focus = "effMaxSG"`,
`selection_rule = "neighborhood"`, `effect_neighborhood = 0.20`).
Section 2.1 of León and Anderson (2026b) sets this rule beside the
alternatives ($`\max_{\mathcal{A}} N`$,
$`\max_{\mathcal{A}} \hat\beta`$, and their smallest-subgroup
counterparts), which Figure 1 of that paper draws on this trial’s
candidate frontier.

**Inference** (`mr_inference = TRUE`): from 5,000 multiplier draws over
every size-eligible candidate the identifier proposes, the corrected
(two-term) estimate of the effect in the identified region $`\hat H`$
and in its complement $`\hat H^{c}`$ (León and Anderson 2026a); the
nominal 95% one-sided *field* lower bound on $`\beta(\hat H)`$, the
*field-s* upper bound on $`\beta(\hat H^{c})`$, the studentized
complement bound, and the Bonferroni pair, a joint 95% claim on both
(León and Anderson 2026b). The field’s own Monte Carlo runs at 1,000
outer and 500 inner draws, the package’s fixed setting on this path. The
only randomness in the analysis is in these draws; they are seeded by
`seedit`, so the vignette’s numbers are reproducible.

## Forest search

One call. Identification, the unadjusted estimates and every
selection-adjusted quantity come out of it together.
`parallel_args = list(plan = "sequential")` keeps the run on one core;
the package’s default plan is parallel.

``` r

t0 <- proc.time()[["elapsed"]]
fit_fs <- forestsearch(
  df.analysis = df,
  outcome.name = "time_months", event.name = "status", treat.name = "hormon", id.name = "id",
  confounders.name = confounders,
  is.RCT = TRUE, est.scale = "hr", seedit = 8316951,
  # the candidate family
  conf.cont_jcuts = list(er = 10, pgr = 10), conf_force = c("er <= 0", "pgr <= 0"),
  n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
  # the identifier, its thresholds and the selection rule
  subgroup_method = "consistency", consistency_method = "resample",
  hr.threshold = 1.00, hr.consistency = 1.00, pconsistency.threshold = 0.90,
  stop_threshold = NULL,                      # evaluate every candidate; no early stop
  sg_focus = "effMaxSG", selection_rule = "neighborhood", effect_neighborhood = 0.20,
  # selection-adjusted inference
  mr_inference = TRUE,
  mr_inference_args = list(ci_method = "field", draws = 5000L, include_complement = TRUE,
                           field_scale_complement = "selected", field_decompose = TRUE,
                           keep_declaration_field = TRUE,
                           declaration_c0 = c(0.70, 0.75, 0.80, 0.85)),
  quiet = TRUE, parallel_args = list(plan = "sequential"))
fit_sec_fs <- proc.time()[["elapsed"]] - t0
```

``` r

g_fs <- fit_fs$mr_inference
cat(sprintf("Identified subgroup: %s\n", paste(fit_fs$sg.harm, collapse = " & ")))
```

    #> Identified subgroup: {er <= 0} & {pgr <= 26}

``` r

cat(sprintf("Fit including inference: %.1f s\n", fit_sec_fs))
```

    #> Fit including inference: 33.7 s

``` r

cat(sprintf("  multiplier stage %.1f s; field blocks %.1f + %.1f s\n", g_fs$timing_seconds,
            g_fs$field$timing_seconds, g_fs$field$complement$timing_seconds))
```

    #>   multiplier stage 1.0 s; field blocks 16.5 + 2.4 s

The identified region is `{er <= 0} & {pgr <= 26}`: no measured estrogen
receptor and progesterone receptor at most 26 fmol/l. Two arguments in
the call are not the defaults and are there for a reason:
`stop_threshold = NULL` has the consistency stage evaluate every
screened candidate rather than stop at the first that passes, and
`field_decompose = TRUE` returns the supporting quantities shown below.
Everything else that is not the default is the specification itself.

### The search in numbers

The fit records how the candidates were whittled down. Supplementary
Table S13 of León and Anderson (2026b) reports the same funnel.

``` r

counts <- fit_fs$find.grps$filter_counts
funnel <- data.frame(
  stage = c("enumerated (at most two cuts, size-eligible)",
            "fittable (ten events in each arm, Cox fit obtained)",
            "clearing the screening threshold, HR > 1.0",
            "passing the consistency screen, p_cons >= 0.90",
            "re-selection family for inference"),
  candidates = c(counts$n_evaluated, counts$n_passed_cox, counts$n_passed_hr,
                 fit_fs$grp.consistency$n_passed, g_fs$n_family))
knitr::kable(funnel, row.names = FALSE, format.args = list(big.mark = ","))
```

| stage                                               | candidates |
|:----------------------------------------------------|-----------:|
| enumerated (at most two cuts, size-eligible)        |      2,211 |
| fittable (ten events in each arm, Cox fit obtained) |      1,410 |
| clearing the screening threshold, HR \> 1.0         |        121 |
| passing the consistency screen, p_cons \>= 0.90     |          9 |
| re-selection family for inference                   |      1,744 |

The re-selection family is larger than the screened set: multiplier
resampling re-selects over every size-eligible candidate under the
screening threshold and the consistency screen, without the per-arm
event minimum, so the family contains the one the search ranked
(Supplementary Section S5).

The candidates that pass the screen, with their size, hazard ratio and
consistency rate, are the points of Figure 1 of the paper. The selection
rule acts on them.

``` r

cand <- as.data.frame(fit_fs$grp.consistency$out_sg$result)
cand$candidate <- ifelse(cand$M.2 == "", cand$M.1, paste(cand$M.1, cand$M.2, sep = " & "))
cand <- cand[order(-cand$N), ]
knitr::kable(cand[, c("candidate", "N", "hr", "Pcons")], digits = 3,
             row.names = FALSE, col.names = c("candidate", "N", "HR", "p_cons"),
             caption = "The candidates passing the consistency screen, largest first.")
```

| candidate                   |   N |    HR | p_cons |
|:----------------------------|----:|------:|-------:|
| {er \<= 0}                  |  82 | 1.951 |   0.96 |
| {er \<= 0} & {pgr \<= 26}   |  75 | 2.222 |   0.99 |
| {pgr \<= 0} & {er \<= 9}    |  71 | 1.684 |   0.90 |
| {er \<= 0} & {pgr \<= 16}   |  70 | 1.918 |   0.95 |
| {er \<= 0} & {pgr \<= 8}    |  66 | 2.129 |   0.97 |
| {er \<= 0} & {pgr \<= 2}    |  62 | 1.886 |   0.92 |
| {er \<= 0} & {size \<= 35}  |  61 | 2.537 |   0.99 |
| {er \<= 0} & {nodes \<= 7}  |  61 | 2.335 |   0.97 |
| {er \<= 0} & !{size \<= 20} |  61 | 2.054 |   0.95 |

The candidates passing the consistency screen, largest first. {.table
.caption-top}

``` r

band <- (1 - 0.20) * max(cand$hr)
within_band <- cand$hr >= band
sel <- cand[within_band, ][which.max(cand$N[within_band]), ]
cat(sprintf("Largest HR over the family: %.2f (N %d). Band: HR >= %.2f.\n",
            max(cand$hr), cand$N[which.max(cand$hr)], band))
```

    #> Largest HR over the family: 2.54 (N 61). Band: HR >= 2.03.

``` r

cat(sprintf("max N(eps) takes the largest candidate within the band: N %d, HR %.2f.\n",
            sel$N, sel$hr))
```

    #> max N(eps) takes the largest candidate within the band: N 75, HR 2.22.

``` r

cat(sprintf("max N alone would take N %d (HR %.2f); max beta the N %d.\n",
            max(cand$N), cand$hr[which.max(cand$N)], cand$N[which.max(cand$hr)]))
```

    #> max N alone would take N 82 (HR 1.95); max beta the N 61.

``` r

stopifnot(sel$N == g_fs$n_selected)   # the rule, re-applied, picks the fit's own selection
```

All of the passing candidates restrict the estrogen receptor, and the
rules differ only in how far that region is trimmed. The
$`\max_{\mathcal{A}}
N(\varepsilon)`$ rule takes the 75 patients with no measured estrogen
receptor and progesterone receptor at most 26 fmol/l;
$`\max_{\mathcal{A}} N`$ would take the 82 with no measured estrogen
receptor; the effect-maximising rules the 61 with, in addition, tumour
size at most 35 mm.

## The selection-adjusted analysis of a region

The nine quantities of Table 1 are read from `fit$mr_inference`. There
is no single accessor for them, so the vignette defines one;
`print(fit)` and `summary(fit)` report the bounds and the diagnostics in
their own layout.

``` r

table1_column <- function(fit, data) {
  g  <- fit$mr_inference; cc <- g$complement
  fb <- g$field; fc <- fb$complement; js <- fb$joint_s
  in_H <- fit$df.est$treat.recommend[match(data$id, fit$df.est$id)] == 0
  stopifnot(!anyNA(in_H), sum(in_H) == g$n_selected, sum(!in_H) == cc$n)
  f3 <- function(x) sprintf("%.3f", x)
  ci <- function(e, l, u) sprintf("%s (%s, %s)", f3(e), f3(l), f3(u))
  c(paste(fit$sg.harm, collapse = " & "),
    sprintf("%d (%d)", sum(in_H), sum(data$status[in_H])),
    ci(g$naive$est, g$naive$lower, g$naive$upper),
    f3(g$debiased$est),
    f3(fb$lower_1s),
    sprintf("%d (%d)", sum(!in_H), sum(data$status[!in_H])),
    ci(cc$naive$est, cc$naive$lower, cc$naive$upper),
    f3(cc$debiased$est),
    f3(fc$upper_1s_s),
    sprintf("%s / %s", f3(js$bonf_lower_H), f3(js$bonf_upper_Hc)))
}
table1_rows <- c("Region", "N (events)", "Unadjusted, 95% CI", "Corrected",
                 "Field, lower bound", "Complement N (events)", "Unadjusted, 95% CI",
                 "Corrected", "Field-s, upper bound", "Bonferroni pair")
knitr::kable(data.frame(row = table1_rows, FS = table1_column(fit_fs, df)),
             col.names = c("", "FS"), row.names = FALSE,
             caption = "Forest search: the identified region and its complement.")
```

|                       | FS                        |
|:----------------------|:--------------------------|
| Region                | {er \<= 0} & {pgr \<= 26} |
| N (events)            | 75 (41)                   |
| Unadjusted, 95% CI    | 2.222 (1.169, 4.221)      |
| Corrected             | 1.276                     |
| Field, lower bound    | 0.634                     |
| Complement N (events) | 611 (258)                 |
| Unadjusted, 95% CI    | 0.605 (0.465, 0.788)      |
| Corrected             | 0.645                     |
| Field-s, upper bound  | 0.804                     |
| Bonferroni pair       | 0.509 / 0.840             |

Forest search: the identified region and its complement. {.table
.caption-top}

The unadjusted hazard ratio in the region is opposite in direction to
the trial’s; the correction removes 0.543 on the log scale and leaves
1.276. The field lower bound is the one-sided 95% bound on the region’s
effect; the field-s upper bound the one-sided 95% bound on the
complement’s; the Bonferroni pair is both at joint 95%, each side at
0.975.

The quantities that support the bounds, as in Supplementary Table S18 of
the paper:

``` r

fb <- g_fs$field; fc <- fb$complement; js <- fb$joint_s
supporting <- data.frame(
  quantity = c("mean of the simulated selection root, region (log scale)",
               "second-order estimate, region",
               "second-order estimate, complement (field-s)",
               "scale ratio rho_c (field-s vs the unstudentized complement bound)",
               "Bonferroni pair: simulated joint coverage",
               "calibrated split gamma (Bonferroni's own alpha/2 is 0.025)",
               "correlation of the aligned roots"),
  value = c(sprintf("%.3f", fb$lambda_mean), sprintf("%.3f", fb$est2),
            sprintf("%.3f", fc$est2_s), sprintf("%.3f", fc$decomp_fields$scale_ratio_c),
            sprintf("%.4f", js$bonf_joint_prob), sprintf("%.3f", js$gamma),
            sprintf("%.3f", js$corr)))
knitr::kable(supporting, row.names = FALSE)
```

| quantity                                                          | value  |
|:------------------------------------------------------------------|:-------|
| mean of the simulated selection root, region (log scale)          | -0.022 |
| second-order estimate, region                                     | 1.305  |
| second-order estimate, complement (field-s)                       | 0.647  |
| scale ratio rho_c (field-s vs the unstudentized complement bound) | 1.002  |
| Bonferroni pair: simulated joint coverage                         | 0.9495 |
| calibrated split gamma (Bonferroni’s own alpha/2 is 0.025)        | 0.025  |
| correlation of the aligned roots                                  | -0.026 |

The second-order estimate subtracts from the corrected estimate the mean
of the simulated selection root; the bound subtracts a quantile of the
same law and is not centred on it. A scale ratio near one means field-s
and the unstudentized complement bound agree on these data.

## The two model-based identifiers

DINA and GRF propose candidates from a surface fitted to the outcome and
have them re-ranked on the same Cox coefficient $`\hat\beta(g)`$ against
the same screening threshold; neither applies the consistency screen.
Supplementary Section S5 of León and Anderson (2026b) sets the three
identifiers side by side, with every setting used here. The selection
rule and the inference are the same. The calls differ from forest
search’s only in the identifier’s own arguments.

``` r

t0 <- proc.time()[["elapsed"]]
fit_dina <- forestsearch(
  df.analysis = df,
  outcome.name = "time_months", event.name = "status", treat.name = "hormon", id.name = "id",
  confounders.name = confounders,
  is.RCT = TRUE, est.scale = "hr", seedit = 8316951,
  conf.cont_jcuts = list(er = 10, pgr = 10), conf_force = c("er <= 0", "pgr <= 0"),
  n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
  subgroup_method = "dina", dina_select_statistic = "effect",     # DINA's own
  hr.threshold = 1.00, stop_threshold = NULL,
  sg_focus = "effMaxSG", selection_rule = "neighborhood", effect_neighborhood = 0.20,
  mr_inference = TRUE,
  mr_inference_args = list(ci_method = "field", draws = 5000L, include_complement = TRUE,
                           field_scale_complement = "selected", field_decompose = TRUE),
  quiet = TRUE, parallel_args = list(plan = "sequential"))
fit_sec_dina <- proc.time()[["elapsed"]] - t0
```

``` r

t0 <- proc.time()[["elapsed"]]
fit_grf <- forestsearch(
  df.analysis = df,
  outcome.name = "time_months", event.name = "status", treat.name = "hormon", id.name = "id",
  confounders.name = confounders,
  is.RCT = TRUE, est.scale = "hr", seedit = 8316951,
  conf.cont_jcuts = list(er = 10, pgr = 10), conf_force = c("er <= 0", "pgr <= 0"),
  n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
  subgroup_method = "grf", grf_selection = "frontier",             # GRF's own
  grf_select_statistic = "effect", grf_depth = 2L, dmin.grf = 0.0, # GRF's own
  hr.threshold = 1.00, stop_threshold = NULL,
  sg_focus = "effMaxSG", selection_rule = "neighborhood", effect_neighborhood = 0.20,
  mr_inference = TRUE,
  mr_inference_args = list(ci_method = "field", draws = 5000L, include_complement = TRUE,
                           field_scale_complement = "selected", field_decompose = TRUE),
  quiet = TRUE, parallel_args = list(plan = "sequential"))
fit_sec_grf <- proc.time()[["elapsed"]] - t0
```

``` r

g_dina <- fit_dina$mr_inference
g_grf  <- fit_grf$mr_inference
cat(sprintf("DINA: %s  (fit %.1f s; re-selection family %d)\n",
            paste(fit_dina$sg.harm, collapse = " & "), fit_sec_dina,
            g_dina$n_family))
```

    #> DINA: {grade3 >= 1} & {pgr <= 20}  (fit 10.3 s; re-selection family 30)

``` r

cat(sprintf("GRF:  %s  (fit %.1f s; re-selection family %d)\n",
            paste(fit_grf$sg.harm, collapse = " & "), fit_sec_grf, g_grf$n_family))
```

    #> GRF:  {er <= 0}  (fit 26.5 s; re-selection family 858)

GRF returns the 82 patients with no measured estrogen receptor, a region
that contains forest search’s; DINA returns the 113 who are grade 3 with
progesterone receptor at most 20 fmol/l, a different classification. For
the two model-based identifiers the candidate family is generated from a
surface fitted to the outcome, so their bounds are conditional on the
family proposed (Supplementary Section S1.9 of León and Anderson 2026b);
forest search’s family is enumerated from the covariates alone.

``` r

member <- function(fit) {
  fit$df.est$treat.recommend[match(df$id, fit$df.est$id)] == 0
}
m <- data.frame(FS = member(fit_fs), DINA = member(fit_dina), GRF = member(fit_grf))
partition <- aggregate(cbind(N = 1, events = df$status), by = m, FUN = sum)
knitr::kable(partition[order(-partition$N), ], row.names = FALSE,
             caption = "The three-way partition of the 686 patients by the three selections.")
```

| FS    | DINA  | GRF   |   N | events |
|:------|:------|:------|----:|-------:|
| FALSE | FALSE | FALSE | 526 |    212 |
| FALSE | TRUE  | FALSE |  78 |     42 |
| TRUE  | FALSE | TRUE  |  40 |     22 |
| TRUE  | TRUE  | TRUE  |  35 |     19 |
| FALSE | FALSE | TRUE  |   7 |      4 |

The three-way partition of the 686 patients by the three selections.
{.table .caption-top}

## Table 1

``` r

table1 <- data.frame(row = table1_rows,
                     FS = table1_column(fit_fs, df),
                     DINA = table1_column(fit_dina, df),
                     GRF = table1_column(fit_grf, df))
caption1 <- sprintf(paste(
  "Selection-adjusted analysis of the three identified regions, hazard-ratio",
  "scale. Unadjusted: two-sided 95%% Wald interval; field and field-s: one-sided",
  "at nominal 0.95; Bonferroni pair: region / complement, joint 0.95.",
  "Whole trial: %d patients, %d events, hazard ratio %.3f (%.3f, %.3f)."),
  nrow(df), sum(df$status), whole_hr[["est"]], whole_hr[[2]], whole_hr[[3]])
knitr::kable(table1, col.names = c("", "FS", "DINA", "GRF"), row.names = FALSE,
             align = "lrrr", caption = caption1)
```

|  | FS | DINA | GRF |
|:---|---:|---:|---:|
| Region | {er \<= 0} & {pgr \<= 26} | {grade3 \>= 1} & {pgr \<= 20} | {er \<= 0} |
| N (events) | 75 (41) | 113 (61) | 82 (45) |
| Unadjusted, 95% CI | 2.222 (1.169, 4.221) | 1.315 (0.746, 2.317) | 1.951 (1.037, 3.673) |
| Corrected | 1.276 | 1.076 | 1.111 |
| Field, lower bound | 0.634 | 0.638 | 0.510 |
| Complement N (events) | 611 (258) | 573 (238) | 604 (254) |
| Unadjusted, 95% CI | 0.605 (0.465, 0.788) | 0.641 (0.489, 0.839) | 0.615 (0.472, 0.801) |
| Corrected | 0.645 | 0.658 | 0.663 |
| Field-s, upper bound | 0.804 | 0.842 | 0.846 |
| Bonferroni pair | 0.509 / 0.840 | 0.572 / 0.887 | 0.404 / 0.874 |

Selection-adjusted analysis of the three identified regions,
hazard-ratio scale. Unadjusted: two-sided 95% Wald interval; field and
field-s: one-sided at nominal 0.95; Bonferroni pair: region /
complement, joint 0.95. Whole trial: 686 patients, 299 events, hazard
ratio 0.695 (0.544, 0.888). {.table .caption-top}

This is Table 1 of León and Anderson (2026b). The calls, the seed and
the draw counts are those of the manuscript’s analysis, so the published
values are reproduced to the digit; a different draw count would move
the corrected estimates and the bounds in the third decimal and leave
the regions and the unadjusted estimates unchanged.

### Reading the table

Each estimate is read by its location against the trial’s own effect,
0.695; the bounds are reported beside it and are not read as tests at a
hazard ratio of one. A corrected estimate above 0.90 is read as
suggestive of limited benefit, and where its lower bound lies below 1.0
a benefit cannot be ruled out.

- **Forest search.** The unadjusted hazard ratio in its region is 2.222,
  opposite in direction to the trial’s. The correction leaves 1.276,
  above 0.90, so suggestive of limited benefit; with a field lower bound
  of 0.634, below 1.0, a benefit cannot be ruled out. The complement,
  89% of the trial, has a corrected estimate of 0.645, a larger benefit
  than the trial’s, with a field-s upper bound of 0.804.
- **DINA and GRF** read the same way: corrected estimates of 1.076 and
  1.111 in their regions, with lower bounds below 1.0, and of 0.658 and
  0.663 in their complements.
- **The two-subgroup claim.** Forest search’s Bonferroni pair is 0.509
  for the region and 0.840 for the complement. The complement’s benefit,
  a hazard reduction of at least 16%, survives the joint adjustment; the
  region excludes only hazard ratios below 0.509, which is no claim of
  harm. A claim that the two subgroups differ is not supported; the
  claim of benefit outside the region is.
- **The substantive reading** is the same under all three partitions: a
  small region whose unfavourable unadjusted estimate the correction
  largely removes, and a complement, most of the trial, with a corrected
  benefit larger than the trial’s. With 41 to 61 events in the regions,
  a claim of harm at a conventional level is not realistic at this
  trial’s size; the bounds are reported, and how much statistical
  assurance to require of them is the analyst’s decision in their own
  setting.

## Declaration at a protected level

Every region above was declared at the identifiers’ default settings,
forest search’s by the conventional screen at $`p^\star = 0.90`$.
Section 4 of León and Anderson (2026b) asks a different question: would
a screen calibrated to protect a clinically meaningful uniform benefit
have declared it? The answer is read for forest search, whose family is
enumerated. At a protected level $`c_0`$ and family-wise level
$`\alpha`$, the calibrated cutoff $`\hat\kappa_\alpha(c_0)`$ is the
$`(1 - \alpha)`$ quantile of the maximum over the family of each
candidate’s standardized statistic
``` math
T(g) \;=\; \frac{\hat\beta(g) - \log c_2}{\sigma_D(g)},
```
with every candidate shifted to sit at $`c_0`$, under the law of the
fit’s own multiplier draws. The calibrated rule admits a candidate when
$`T(g) \ge \hat\kappa_\alpha(c_0)`$; it holds the family-wise
probability of declaring any harm subgroup at $`\alpha`$ whenever every
candidate’s true effect is $`c_0`$ or better.
[`fs_declaration_calibration()`](https://larry-leon.github.io/forestsearch/reference/fs_declaration_calibration.md)
reads it from the draws the fit above already made
(`keep_declaration_field = TRUE`, with the grid passed as
`declaration_c0`); nothing is redrawn and nothing else in the fit
changes.

``` r

c0_grid <- c(0.70, 0.75, 0.80, 0.85)
cal    <- fs_declaration_calibration(fit_fs, alpha = 0.10, c0 = c0_grid)
cal_c2 <- fs_declaration_calibration(fit_fs, alpha = 0.10)   # c0 = NULL: the claim threshold
tb <- cal$c0$table
fld <- g_fs$declaration_field
T_H   <- cal$T_hat[[fld$family_id[fld$meta$selected_index]]]
T_max <- max(cal$T_hat)
kappa <- c(tb$kappa_hat, cal_c2$kappa_hat)
declaration <- data.frame(
  c0 = c(sprintf("HR %.2f%s", tb$c0, ifelse(tb$c0 == 0.75, " (primary)", "")),
         "HR 1.00 (claim threshold)"),
  kappa   = sprintf("%.3f", kappa),
  implied = sprintf("%.3f", 2 * pnorm(kappa) - 1),
  fw      = sprintf("%.3f", c(tb$fw_size, cal_c2$fw_size)),
  admits  = c(lengths(cal$c0$admitted_calibrated),
              length(cal_c2$admitted_calibrated)))
knitr::kable(declaration, row.names = FALSE, align = "lrrrr",
             col.names = c("protected level c0", "kappa-hat", "implied level",
                           "FW size", "admits"),
             caption = paste(
               "Declaration calibration at alpha = 0.10. kappa-hat: the calibrated",
               "cutoff; implied level: 2*Phi(kappa-hat) - 1; FW size: the family-wise",
               "size of the p* = 0.90 screen at c0; admits: candidates the calibrated",
               "rule admits."))
```

| protected level c0        | kappa-hat | implied level | FW size | admits |
|:--------------------------|----------:|--------------:|--------:|-------:|
| HR 0.70                   |     2.578 |         0.990 |   0.558 |      0 |
| HR 0.75 (primary)         |     2.726 |         0.994 |   0.667 |      0 |
| HR 0.80                   |     2.891 |         0.996 |   0.769 |      0 |
| HR 0.85                   |     3.039 |         0.998 |   0.844 |      0 |
| HR 1.00 (claim threshold) |     3.503 |         1.000 |   0.974 |      0 |

Declaration calibration at alpha = 0.10. kappa-hat: the calibrated
cutoff; implied level: 2*Phi(kappa-hat) - 1; FW size: the family-wise
size of the p* = 0.90 screen at c0; admits: candidates the calibrated
rule admits. {.table .caption-top}

``` r

cat(sprintf("T of the declared region: %.3f; largest T in the family: %.3f\n", T_H, T_max))
```

    #> T of the declared region: 2.438; largest T in the family: 2.564

``` r

cat(sprintf("The p* = 0.90 screen's cutoff, as the package applies it: %.3f\n", cal$z_pstar))
```

    #> The p* = 0.90 screen's cutoff, as the package applies it: 1.621

The screen’s cutoff is the $`p^\star = 0.90`$ rule as the package
applies it, after the rate is rounded to `pconsistency.digits = 2`,
which is why it sits a little below the 1.645 of an unrounded 0.90.

At the primary protected level, a hazard ratio of 0.75, the calibrated
cutoff is 2.726, an implied screening level of 0.994; the region’s
statistic is 2.438 and the largest in the family 2.564, so the
calibrated screen admits no candidate at that level, and likewise at
every protected level on the grid and at the claim threshold. The
FW-size column is the diagnostic of the paper’s equation (8): the
probability that the $`p^\star = 0.90`$ screen as executed declares some
harm subgroup if every candidate carried a uniform benefit of $`c_0`$ —
a property of this family and covariate design, not the trial’s realized
error rate; at the claim threshold it is a boundary size. To run the
screen at the implied level, set `pconsistency.threshold` to it and
`pconsistency.digits` to at least 5, so the level is not rounded away;
the package’s `pstar_settable` and `pstar_fine` columns in
`cal$c0$table` give what the digit grid can express. The supplement’s
Section S1.8.4 adds a realized-rate simulation on this design, which is
not run here. The selection and every bound above are those of the fit
as run; the calibration reads the fit and alters none of them.

## Where to go next

- The [post-selection
  article](https://larry-leon.github.io/forestsearch/articles/survival-post-selection.html)
  reads the bounds in more depth on the same trial: the re-selection
  frequency p̂ as a diagnostic, membership agreement, and what each
  construction does and does not say.
- The package website carries the original workflow of León et
  al. (2024) — the full argument listing, the bootstrap bias correction
  and the cross-validation — as the article *Exploratory subgroup
  identification with forestsearch (the 2024 workflow)*, and the
  methodology and simulation articles.
- [`?forestsearch`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  documents every argument, with examples for this analysis and for the
  binary and continuous ACTG 175 analyses of León and Anderson (2026b),
  and a section on the options not tested against the manuscript
  analyses.

## References

Amatya, A. K., Fiero, M. H., Bloomquist, E. W., Sinha, A. K., Lemery, S.
J., Singh, H., Ibrahim, A., Donoghue, M., Fashoyin-Aje, L. A., de Claro,
R. A., Gormley, N. J., Amiri-Kordestani, L., Sridhara, R., Theoret, M.
R., Kluetz, P. G., Pazdur, R., Beaver, J. A. and Tang, S. (2021).
Subgroup analyses in oncology trials: Regulatory considerations and case
examples. *Clinical Cancer Research* 27(21), 5753–5756.

Athey, S., Tibshirani, J. and Wager, S. (2019). Generalized random
forests. *The Annals of Statistics* 47(2), 1148–1178.

Cui, Y., Kosorok, M. R., Sverdrup, E., Wager, S. and Zhu, R. (2023).
Estimating heterogeneous treatment effects with right-censored data via
causal survival forests. *Journal of the Royal Statistical Society
Series B* 85(2), 179–211.

European Medicines Agency (2019). Guideline on the investigation of
subgroups in confirmatory clinical trials. EMA/CHMP/539146/2013,
Committee for Medicinal Products for Human Use.

Gao, Z. and Hastie, T. (2025). Estimating heterogeneous treatment
effects for general responses. *Biometrics* 81(4), ujaf162.

León, L. F. and Anderson, K. M. (2026a). Inference for standard trial
estimands after data-driven subgroup discovery. arXiv:2609.38361.
<https://doi.org/10.48550/arXiv.2609.38361>

León, L. F. and Anderson, K. M. (2026b). Interpretable subgroup
identification across survival and generalized-linear outcomes.
arXiv:2610.28379. <https://doi.org/10.48550/arXiv.2610.28379>

León, L. F., Jemielita, T., Guo, Z., Marceau West, R. and Anderson, K.
M. (2024). Exploratory subgroup identification in the heterogeneous Cox
model: A relatively simple procedure. *Statistics in Medicine*.
<https://doi.org/10.1002/sim.10163>

``` r

sessionInfo()
```

    #> R version 4.6.2 (2026-10-29)
    #> Platform: x86_64-pc-linux-gnu
    #> Running under: Ubuntu 24.04.5 LTS
    #>
    #> Matrix products: default
    #> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
    #> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    #>
    #> locale:
    #>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8
    #>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8
    #>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C
    #> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C
    #>
    #> time zone: UTC
    #> tzcode source: system (glibc)
    #>
    #> attached base packages:
    #> [1] stats     graphics  grDevices utils     datasets  methods   base
    #>
    #> other attached packages:
    #> [1] future_1.80.0      survival_3.8-12    forestsearch_0.4.0
    #>
    #> loaded via a namespace (and not attached):
    #>  [1] gt_1.4.0             generics_0.1.4       xml2_1.6.0
    #>  [4] shape_1.4.6.1        lattice_0.23-1       listenv_1.1.0
    #>  [7] digest_0.6.39        magrittr_2.0.5       grf_2.6.1
    #> [10] evaluate_1.0.5       grid_4.6.2           RColorBrewer_1.1-3
    #> [13] iterators_1.0.14     policytree_1.2.5     mvtnorm_1.4-2
    #> [16] fastmap_1.2.0        foreach_1.5.2        jsonlite_2.0.0
    #> [19] Matrix_1.7-6         glmnet_5.2           scales_1.4.0
    #> [22] weightedsurv_0.1.0   codetools_0.2-20     cli_3.6.6
    #> [25] rlang_1.3.0          parallelly_1.48.0    future.apply_1.21.0
    #> [28] splines_4.6.2        yaml_2.3.12          otel_0.2.0
    #> [31] tools_4.6.2          parallel_4.6.2       doFuture_1.4.0
    #> [34] dplyr_1.2.1          ggplot2_4.0.3        globals_0.20.0
    #> [37] vctrs_0.7.3          R6_2.6.1             lifecycle_1.0.5
    #> [40] randomForest_4.7-1.2 fs_2.1.0             pkgconfig_2.0.3
    #> [43] progressr_1.1.0      pillar_1.11.1        gtable_0.3.6
    #> [46] data.table_1.18.6.1  glue_1.8.1           Rcpp_1.1.2
    #> [49] xfun_0.62            tibble_3.3.1         tidyselect_1.2.1
    #> [52] knitr_1.52           farver_2.1.2         patchwork_1.3.2
    #> [55] htmltools_0.5.9      rmarkdown_2.32       compiler_4.6.2
    #> [58] S7_0.2.2
