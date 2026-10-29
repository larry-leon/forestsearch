# Post-selection inference for an identified survival subgroup

> **A website article.** In forestsearch 0.4.0 this was a CRAN vignette;
> from 0.4.1 it is an article on the package website, because the full
> analysis takes a few minutes to run, too long for CRAN’s package
> checks. The text is otherwise unchanged; where it says “vignette”,
> read “article”.

``` r

knitr::opts_chunk$set(comment = "#>", fig.align = "center")
library(forestsearch)
library(survival)
```

A subgroup discovered in a trial is reported through the standard
analysis the trial would have run had that subgroup been the only one:
the treatment coefficient of the working model fitted within it, a Cox
log hazard ratio here. Selection by an effect-driven rule — rewarding
large estimated effects subject to screening and size criteria —
inflates that coefficient. This is the winner’s curse, and it
invalidates a naive interval whether or not the identifier’s own
internal estimand was de-biased by cross-fitting (León and Anderson
2026a). The same selection distorts the estimate for the patients
outside the subgroup.

**What is corrected, and how.** The target is the standard-analysis
effect of the selected region, with the candidate family held fixed:
$`\beta(\hat H)`$ for the region and $`\beta(\hat H^{c})`$ for its
complement. Two routes remove the optimism of selection.

- **The full bootstrap** re-runs the whole pipeline — enumerate, screen,
  select — on each resample, and corrects the observed estimate by the
  bootstrap estimate of its selection bias. It is the procedure of León
  et al. (2024). It is costly, because every replicate repeats the
  search.
- **Multiplier resampling** (MR) never refits. Once the candidate family
  and every candidate’s effect and influence function are in hand,
  resampling the search reduces to perturbing all candidate effects
  jointly with one set of multiplier weights and re-selecting; the
  correction is read from the re-selected effects, and an
  infinitesimal-jackknife (IJ) interval comes from the same draws. León
  and Anderson (2026a) show the IJ interval’s coverage under standard
  regularity and an explicit condition on the limiting competition
  between candidates, and that MR matches the full bootstrap on a fixed
  family.
- **Model-generated families.** For DINA and GRF the candidate family is
  proposed from a surface fitted to the outcome, so a bootstrap
  regenerates it. There MR remains valid for the same conditional target
  — the effect of the selected region given the family proposed — which
  the bootstrap reproduces only under further conditions.

**This vignette** runs the GBSG analysis of León and Anderson (2026a)
for the three identifiers — forest search (FS), DINA and GRF — at the
manuscript’s settings, with the multiplier correction in full (5,000
draws) and the full bootstrap implemented and run at a small number of
replicates to show the machinery and its cost. The manuscript’s Table S2
(Supplementary Section S6), computed with B = 1,000 bootstrap
replicates, is reproduced for comparison.
[`vignette("forestsearch")`](https://larry-leon.github.io/forestsearch/articles/forestsearch.md)
runs the companion analysis of León and Anderson (2026b) on the same
trial — the search, the one-sided bounds and the declaration
calibration; this vignette is about the two-sided correction and the
bootstrap it replaces.

## The data

``` r

df <- survival::gbsg
df <- within(df, {
  id          <- seq_len(nrow(df))
  time_months <- rfstime / 30.4375
  grade3      <- ifelse(grade == "3", 1, 0)
})
confounders <- c("age", "meno", "size", "grade3", "nodes", "pgr", "er")
stopifnot(nrow(df) == 686, sum(df$status) == 299)
whole <- coxph(Surv(time_months, status) ~ hormon, data = df)
sprintf("Whole trial: %d patients, %d events, hazard ratio %.3f (%.3f, %.3f)",
        nrow(df), sum(df$status), exp(coef(whole)), exp(confint(whole))[1],
        exp(confint(whole))[2])
```

    #> [1] "Whole trial: 686 patients, 299 events, hazard ratio 0.695 (0.544, 0.888)"

## The settings

The analysis of León and Anderson (2026a), Section 6, and the calls
below follow the document that produced its Table S2:

- **Candidates** are conjunctions of at most two covariate cuts with at
  least 60 patients and 10 events in each arm. Forest search carries ten
  candidate cuts on each receptor with zero forced in as a cut of its
  own; DINA and GRF propose their candidates from a fitted surface on
  the covariates as given.
- **Forest search** screens at a hazard ratio of 1.0 and requires
  splitting consistency at $`p^\star = 0.90`$, evaluated in closed form
  (`consistency_method = "resample"`); `stop_threshold = NULL` has every
  screened candidate evaluated.
- **The selection rule** is $`\max_{\mathcal{A}} N(\varepsilon)`$ at
  $`\varepsilon = 0.10`$: the largest admitted region whose effect is
  within 10% of the largest effect over the family. The companion
  vignette uses $`\varepsilon = 0.20`$ and returns a larger
  forest-search region (N 75); at $`\varepsilon = 0.10`$ the rule
  returns the 61 patients with no measured estrogen receptor and tumour
  size at most 35 mm — the effect-maximising region of that vignette’s
  candidate table.
- **DINA** is re-ranked on the Cox coefficient against a floor of 0.90,
  with a two-fold cross-fit surface; **GRF** uses the frontier of its
  covariate grid, re-ranked on the Cox coefficient, at depth 2 with a
  doubly-robust-score floor of 0.
- **Multiplier resampling** at 5,000 draws, reported as the two-term
  corrected estimate with its IJ two-sided 95% interval
  (`ci_method = "ij"` for DINA and GRF). The forest-search fit asks for
  `ci_method = "field"` with `field_recovery = TRUE`: the two-term
  estimate and the IJ interval are the same under both settings, and the
  field adds the one-sided bounds of the companion vignette and the
  membership-recovery diagnostics used below.

## The three fits

The fits run on one core here
(`parallel_args = list(plan = "sequential")`); the bootstrap below uses
up to two, the limit a CRAN check imposes. The package’s default plan is
parallel.

``` r

t0 <- proc.time()[["elapsed"]]
fit_fs <- forestsearch(
  df.analysis = df,
  outcome.name = "time_months", event.name = "status", treat.name = "hormon", id.name = "id",
  confounders.name = confounders,
  is.RCT = TRUE, est.scale = "hr", seedit = 8316951,
  conf.cont_jcuts = list(er = 10, pgr = 10), conf_force = c("er <= 0", "pgr <= 0"),
  n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
  subgroup_method = "consistency", consistency_method = "resample",
  hr.threshold = 1.00, hr.consistency = 1.00, pconsistency.threshold = 0.90,
  stop_threshold = NULL,
  sg_focus = "effMaxSG", selection_rule = "neighborhood", effect_neighborhood = 0.10,
  mr_inference = TRUE,
  mr_inference_args = list(ci_method = "field", draws = 5000L, field_recovery = TRUE),
  quiet = TRUE, parallel_args = list(plan = "sequential"))
sec_fs <- proc.time()[["elapsed"]] - t0
```

``` r

t0 <- proc.time()[["elapsed"]]
fit_dina <- forestsearch(
  df.analysis = df,
  outcome.name = "time_months", event.name = "status", treat.name = "hormon", id.name = "id",
  confounders.name = confounders,
  is.RCT = TRUE, est.scale = "hr", seedit = 8316951,
  subgroup_method = "dina", dina_select_statistic = "effect",
  dina_args = list(seed = 8316951, n_folds = 2),
  hr.threshold = 0.90, n.min = 60,
  sg_focus = "effMaxSG", selection_rule = "neighborhood", effect_neighborhood = 0.10,
  mr_inference = TRUE, mr_inference_args = list(ci_method = "ij", draws = 5000L),
  quiet = TRUE, parallel_args = list(plan = "sequential"))
sec_dina <- proc.time()[["elapsed"]] - t0
```

``` r

t0 <- proc.time()[["elapsed"]]
fit_grf <- forestsearch(
  df.analysis = df,
  outcome.name = "time_months", event.name = "status", treat.name = "hormon", id.name = "id",
  confounders.name = confounders,
  is.RCT = TRUE, est.scale = "hr", seedit = 8316951,
  subgroup_method = "grf", grf_selection = "frontier", grf_select_statistic = "effect",
  grf_depth = 2L, dmin.grf = 0.0, frac.tau = 0.8, n.min = 60,
  sg_focus = "effMaxSG", selection_rule = "neighborhood", effect_neighborhood = 0.10,
  mr_inference = TRUE, mr_inference_args = list(ci_method = "ij", draws = 5000L),
  quiet = TRUE, parallel_args = list(plan = "sequential"))
sec_grf <- proc.time()[["elapsed"]] - t0
```

``` r

fits <- list(FS = fit_fs, DINA = fit_dina, GRF = fit_grf)
secs <- c(FS = sec_fs, DINA = sec_dina, GRF = sec_grf)
for (m in names(fits)) {
  gm <- fits[[m]]$mr_inference
  cat(sprintf("%-4s %-28s N %3d   fit %5.1f s (multiplier stage %.1f s)   family %d\n",
              m, paste(fits[[m]]$sg.harm, collapse = " & "), gm$n_selected, secs[[m]],
              gm$timing_seconds, gm$n_family))
}
```

    #> FS   {er <= 0} & {size <= 35}     N  61   fit  31.7 s (multiplier stage 1.0 s)   family 1744
    #> DINA {grade3 >= 1} & {pgr <= 10}  N  89   fit   1.4 s (multiplier stage 0.3 s)   family 84
    #> GRF  {er <= 0}                    N  82   fit   9.7 s (multiplier stage 0.4 s)   family 858

## The multiplier correction

For each identifier, `fit$mr_inference` holds the unadjusted (naive)
estimate of the selected region, the selection bias the draws estimate,
and the two-term corrected estimate with its IJ interval; `$complement`
holds the same for the patients outside the region. The family size is
the number of candidates the draws re-select over.

``` r

g <- fit_fs$mr_inference
cat(sprintf("Selection bias (log scale): %.3f over a family of %d candidates\n",
            g$selection_bias, g$n_family))
```

    #> Selection bias (log scale): 0.555 over a family of 1744 candidates

``` r

cat(sprintf("Region:     naive %.2f (%.2f, %.2f)  ->  corrected %.2f (%.2f, %.2f)\n",
            g$naive$est, g$naive$lower, g$naive$upper,
            g$debiased$est, g$debiased$lower, g$debiased$upper))
```

    #> Region:     naive 2.54 (1.25, 5.17)  ->  corrected 1.44 (0.62, 3.36)

``` r

cc <- g$complement
cat(sprintf("Complement: naive %.2f (%.2f, %.2f)  ->  corrected %.2f (%.2f, %.2f)\n",
            cc$naive$est, cc$naive$lower, cc$naive$upper,
            cc$debiased$est, cc$debiased$lower, cc$debiased$upper))
```

    #> Complement: naive 0.61 (0.47, 0.79)  ->  corrected 0.64 (0.40, 1.03)

The correction moves the forest-search region’s estimate from an
apparent harm to a value near one, and widens the interval: the IJ
interval accounts for the selection’s own variability, which a naive
interval ignores. The complement moves the other way and by much less —
it is most of the trial, and the selection barely touches it.

## The full bootstrap

[`forestsearch_bootstrap_dofuture()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md)
draws a nonparametric bootstrap sample, re-runs the identifier on it —
the whole search, the screen and the selection, with the front-end
surfaces refitted for DINA and GRF — and records the selected region’s
effect on the resample and on the original data. The bias-corrected
estimate is the observed estimate less the bootstrap estimate of its
selection bias, averaged over replicates (the correction of León et al.
2024); its interval uses an infinitesimal-jackknife variance of the
bootstrap mean with a finite-B correction. Multiplier resampling is not
run inside the replicates (`mr_in_replicates = FALSE`, the default), so
each replicate costs one search.

The manuscript used B = 1,000. Here B is small, to show the
implementation and its cost: the point estimates are usable as an
illustration, the intervals are not — the finite-B variance needs B in
the hundreds, the package warns below 100, and at small B the variance
can fall back to its uncorrected form. Read the interval columns of the
bootstrap tables below as what the machinery produces at this B, not as
results.

Replicates are independent, so the bootstrap parallelises across them:
[`forestsearch_bootstrap_dofuture()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md)
dispatches them through the **future** framework, and
`parallel_args = list(plan = "multisession", workers = k)` runs them on
k processes. The document behind Table S2 ran its B = 1,000 on 11
workers. This vignette uses at most two, the limit a CRAN check imposes
([`future::availableCores()`](https://parallelly.futureverse.org/reference/availableCores.html)
returns what the host allows); for a real run, set `n_workers` below to
the cores you have.

``` r

n_workers <- min(2L, future::availableCores())   # set to your cores for a real run
boot_plan <- if (n_workers > 1L) list(plan = "multisession", workers = n_workers) else
  list(plan = "sequential")
B <- c(FS = 25L, DINA = 25L, GRF = 10L)
boots <- list(); boot_secs <- c()
for (m in names(fits)) {
  t0 <- proc.time()[["elapsed"]]
  boots[[m]] <- forestsearch_bootstrap_dofuture(
    fs.est = fits[[m]], nb_boots = B[[m]], seed = 8316951L,
    mr_in_replicates = FALSE, parallel_args = boot_plan)
  boot_secs[m] <- proc.time()[["elapsed"]] - t0
}
round(boot_secs, 1)
```

    #>   FS DINA  GRF
    #> 71.1 31.8 34.0

The package’s own summaries, shown for forest search. The first is the
bootstrap table: the naive estimate beside the bootstrap-corrected one,
for the region and its complement.

``` r

summ_fs <- summarize_bootstrap_results(sgharm = fit_fs$sg.harm, boot_results = boots$FS,
                                       create_plots = FALSE, est.scale = "hr")
```

    #>
    #> ===============================================================
    #>            BOOTSTRAP ANALYSIS SUMMARY
    #> ===============================================================
    #>
    #> IDENTIFIED SUBGROUP:
    #> -------------------------------------------------------------
    #>   H: {er <= 0} & {size <= 35}
    #>
    #> BOOTSTRAP SUCCESS METRICS:
    #> -------------------------------------------------------------
    #>   Total iterations:              25
    #>   Successful subgroup ID:        23 (92.0%)
    #>   Failed to find subgroup:       2 (8.0%)
    #>
    #> TIMING ANALYSIS:
    #> -------------------------------------------------------------
    #> Overall:
    #>   Total bootstrap time:          1.17 minutes (0.02 hours)
    #>   Average per iteration:         0.05 min (2.8 sec)
    #>   Projected for 1000 boots:      46.65 min (0.78 hrs)

``` r

summ_fs$table
```

[TABLE]

[`mr_estimates_table()`](https://larry-leon.github.io/forestsearch/reference/mr_estimates_table.md)
puts the multiplier correction into the same layout: it takes the
bootstrap table as its template — the descriptive columns and the naive
estimate — and replaces the corrected column by the two-term estimate
with its IJ interval, so the two corrections can be read side by side in
one format.

``` r

mr_estimates_table(fit_fs, boot_results = boots$FS, est.scale = "hr", digits = 2)
```

[TABLE]

## Naive, bootstrap and multiplier side by side

The three columns for the three identifiers, assembled from the fit
objects and the bootstrap objects. The bootstrap column carries the
point estimate only: at this B its interval is not a result (the
package’s own table above shows what the machinery returns for it).

``` r

three_columns <- function(method, fit, boot) {
  g <- fit$mr_inference; cc <- g$complement
  n_H <- g$n_selected; n_Hc <- cc$n
  f2 <- function(e, l, u) sprintf("%.2f (%.2f, %.2f)", e, l, u)
  data.frame(
    method = method, region = c("H", "Hc"),
    n = c(sprintf("%d (%d%%)", n_H, round(100 * n_H / (n_H + n_Hc))),
          sprintf("%d (%d%%)", n_Hc, round(100 * n_Hc / (n_H + n_Hc)))),
    naive = c(f2(g$naive$est, g$naive$lower, g$naive$upper),
              f2(cc$naive$est, cc$naive$lower, cc$naive$upper)),
    bootstrap = sprintf("%.2f", c(boot$H_estimates$H2, boot$Hc_estimates$H2)),
    multiplier = c(f2(g$debiased$est, g$debiased$lower, g$debiased$upper),
                   f2(cc$debiased$est, cc$debiased$lower, cc$debiased$upper)))
}
side_by_side <- rbind(three_columns("FS", fit_fs, boots$FS),
                      three_columns("DINA", fit_dina, boots$DINA),
                      three_columns("GRF", fit_grf, boots$GRF))
boot_label <- sprintf("Full bootstrap estimate (B = %d / %d / %d)",
                      B[["FS"]], B[["DINA"]], B[["GRF"]])
knitr::kable(side_by_side, row.names = FALSE, align = "llrrrr",
             col.names = c("Method", "SG", "n (%)", "Naive", boot_label,
                           "Multiplier (5,000 draws)"),
             caption = paste("This vignette's run: hazard ratios with 95% intervals;",
                             "the bootstrap estimate at a small B, without an interval."))
```

| Method | SG | n (%) | Naive | Full bootstrap estimate (B = 25 / 25 / 10) | Multiplier (5,000 draws) |
|:---|:---|---:|---:|---:|---:|
| FS | H | 61 (9%) | 2.54 (1.25, 5.17) | 1.88 | 1.44 (0.62, 3.36) |
| FS | Hc | 625 (91%) | 0.61 (0.47, 0.79) | 0.63 | 0.64 (0.40, 1.03) |
| DINA | H | 89 (13%) | 1.60 (0.88, 2.91) | 0.83 | 1.11 (0.45, 2.75) |
| DINA | Hc | 597 (87%) | 0.63 (0.48, 0.82) | 0.66 | 0.66 (0.40, 1.09) |
| GRF | H | 82 (12%) | 1.95 (1.04, 3.67) | 1.29 | 1.05 (0.47, 2.34) |
| GRF | Hc | 604 (88%) | 0.61 (0.47, 0.80) | 0.73 | 0.66 (0.41, 1.06) |

This vignette’s run: hazard ratios with 95% intervals; the bootstrap
estimate at a small B, without an interval. {.table .caption-top}

Table S2 of León and Anderson (2026a), as published, with the full
bootstrap at B = 1,000:

``` r

table_s2 <- data.frame(
  method = c("FS", "FS", "DINA", "DINA", "GRF", "GRF"),
  region = rep(c("H", "Hc"), 3),
  n = c("61 (9%)", "625 (91%)", "89 (13%)", "597 (87%)", "82 (12%)", "604 (88%)"),
  naive = c("2.54 (1.25, 5.17)", "0.61 (0.47, 0.79)", "1.60 (0.88, 2.91)",
            "0.63 (0.48, 0.82)", "1.95 (1.04, 3.67)", "0.61 (0.47, 0.80)"),
  bootstrap = c("1.94 (0.89, 4.21)", "0.63 (0.43, 0.93)", "0.89 (0.38, 2.07)",
                "0.66 (0.40, 1.11)", "1.04 (0.49, 2.24)", "0.66 (0.39, 1.13)"),
  multiplier = c("1.44 (0.62, 3.36)", "0.64 (0.40, 1.03)", "1.11 (0.45, 2.75)",
                 "0.66 (0.40, 1.09)", "1.05 (0.47, 2.34)", "0.66 (0.41, 1.06)"))
knitr::kable(table_s2, row.names = FALSE, align = "llrrrr",
             col.names = c("Method", "SG", "n (%)", "Naive", "Full bootstrap (B = 1,000)",
                           "Multiplier"),
             caption = "Table S2 of León and Anderson (2026a), as published.")
```

| Method | SG | n (%) | Naive | Full bootstrap (B = 1,000) | Multiplier |
|:---|:---|---:|---:|---:|---:|
| FS | H | 61 (9%) | 2.54 (1.25, 5.17) | 1.94 (0.89, 4.21) | 1.44 (0.62, 3.36) |
| FS | Hc | 625 (91%) | 0.61 (0.47, 0.79) | 0.63 (0.43, 0.93) | 0.64 (0.40, 1.03) |
| DINA | H | 89 (13%) | 1.60 (0.88, 2.91) | 0.89 (0.38, 2.07) | 1.11 (0.45, 2.75) |
| DINA | Hc | 597 (87%) | 0.63 (0.48, 0.82) | 0.66 (0.40, 1.11) | 0.66 (0.40, 1.09) |
| GRF | H | 82 (12%) | 1.95 (1.04, 3.67) | 1.04 (0.49, 2.24) | 1.05 (0.47, 2.34) |
| GRF | Hc | 604 (88%) | 0.61 (0.47, 0.80) | 0.66 (0.39, 1.13) | 0.66 (0.41, 1.06) |

Table S2 of León and Anderson (2026a), as published. {.table
.caption-top}

- **Naive and multiplier columns.** The regions, their sizes, the naive
  estimates and the multiplier estimates and intervals are those of the
  published table at its precision: the calls, the seed and the draw
  count are the manuscript’s.
- **The bootstrap column** at this B is an illustration of the
  implementation. At B = 1,000 the manuscript finds the two corrections
  agree closely for GRF (1.04 and 1.05), remain apart for forest search
  (1.94 and 1.44, both well below the naive 2.54), and fall on opposite
  sides of one for DINA (0.89 and 1.11) — the two routes estimate the
  same selection bias from different perturbations of the search, and
  DINA’s family is itself regenerated by the bootstrap.
- **Reading, by location.** The forest-search region’s corrected
  estimates, 1.44 by multiplier resampling and 1.94 by the manuscript’s
  bootstrap, sit far below its naive 2.54: most of the apparent harm was
  selection. The complement’s corrected estimate, 0.64, is a benefit at
  least as large as the trial’s own 0.695 for the patients outside the
  region. Neither is read as a test at a hazard ratio of one.

### The cost of the two corrections

``` r

cost <- data.frame(
  method = names(fits),
  mr_seconds = sapply(fits, function(f) f$mr_inference$timing_seconds),
  bootstrap_B = B[names(fits)],
  bootstrap_seconds = boot_secs[names(fits)])
cost$seconds_per_replicate <- cost$bootstrap_seconds / cost$bootstrap_B
cost$projected_B_1000 <- round(cost$seconds_per_replicate * 1000)
plan_label <- if (n_workers > 1L) sprintf("%d workers", n_workers) else "one process"
knitr::kable(cost, row.names = FALSE, digits = 2,
             col.names = c("Method", "MR (s)", "B", "bootstrap (s)", "s per replicate",
                           "projected at B = 1,000 (s)"),
             caption = sprintf(paste(
               "Wall-clock timings from this run: the multiplier stage under the fits'",
               "sequential plan, the bootstrap on %s. The projection is linear in B",
               "at that worker count."), plan_label))
```

| Method | MR (s) |   B | bootstrap (s) | s per replicate | projected at B = 1,000 (s) |
|:-------|-------:|----:|--------------:|----------------:|---------------------------:|
| FS     |   0.99 |  25 |         71.12 |            2.84 |                       2845 |
| DINA   |   0.25 |  25 |         31.82 |            1.27 |                       1273 |
| GRF    |   0.35 |  10 |         33.98 |            3.40 |                       3398 |

Wall-clock timings from this run: the multiplier stage under the fits’
sequential plan, the bootstrap on 2 workers. The projection is linear in
B at that worker count. {.table .caption-top}

The multiplier correction is seconds; the bootstrap is minutes at the B
its intervals need, and the ratio widens with B. The seconds are those
of this machine and this package version: the MR column is the whole
multiplier stage as the fit records it (`timing_seconds`); the bootstrap
column is the whole
[`forestsearch_bootstrap_dofuture()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md)
call on 2 workers, in wall-clock seconds, so the per-replicate figure
and the projection fall as workers are added. The manuscript’s figures —
72 s to about 10 minutes at B = 1,000 on 11 workers, 0.1 to 0.4 s for
the multiplier stage — were recorded under forestsearch 0.2.0 on another
machine, and are not comparable second for second; the order of
magnitude between the two routes is the point.

## Diagnostics

### The re-selection frequency

$`\hat p(\hat H)`$ is the frequency with which the re-selection picks
exactly the identified region across the multiplier draws. It is a
recorded diagnostic — computed after the correction, read by no
construction, and estimated from the same draws — so it is a flag, not
an instrument.

``` r

ph  <- g$reselection$p_hat
cat(sprintf("p-hat(H) = %.3f over a family of %d candidates\n",
            unname(ph[[g$selected_label]]), g$n_family))
```

    #> p-hat(H) = 0.237 over a family of 1744 candidates

``` r

print(round(sort(ph, decreasing = TRUE)[1:5], 3))
```

    #>  q1.1 & q27.1  q1.1 & q18.1  q1.1 & q31.1  q1.1 & q26.0 q10.0 & q30.0
    #>         0.237         0.140         0.104         0.058         0.052

A small $`\hat p`$ says the exact boundary of the region is not stable
under perturbation; by itself it is neither evidence for nor against the
correction.

### What $`\hat p`$ alone cannot say

$`\hat p`$ counts exact re-selections. Over a family in which most
candidates differ from $`\hat H`$ by a single cut, a small exact-match
rate is compatible with the draws re-selecting near-twins of $`\hat H`$
or unrelated regions, and $`\hat p`$ does not distinguish them.
Membership agreement does: each draw’s re-selection is compared with
$`\hat H`$ as a set of patients, in the vocabulary of the
cross-validation agreement metrics of
[`forestsearch_Kfold()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_Kfold.md)
(`sens_H`, `ppv_H`). These diagnostics are computed from the field
draws, which the forest-search fit asked for.

``` r

rc <- g$field$recovery
cat(sprintf("sens_H = %.3f   mean share of H retained by a re-selection\n", rc$sens_H))
```

    #> sens_H = 0.543   mean share of H retained by a re-selection

``` r

cat(sprintf("ppv_H  = %.3f   mean share of a re-selection that lies in H\n", rc$ppv_H))
```

    #> ppv_H  = 0.465   mean share of a re-selection that lies in H

``` r

cat(sprintf("containment q10 / q50 / q90 = %.3f / %.3f / %.3f\n", rc$q10, rc$q50, rc$q90))
```

    #> containment q10 / q50 / q90 = 0.000 / 0.656 / 0.967

``` r

cat(sprintf("share of draws containing all of H = %.3f  (over %d draws)\n",
            rc$share_equal_1, rc$n_used))
```

    #> share of draws containing all of H = 0.020  (over 986 draws)

Read together: the median draw retains 66% of the region’s members, at
least a tenth of the draws retain none of them, and 2% of draws contain
all of it. Exact re-selection is rare; what moves from draw to draw is
the boundary. These diagnostics measure re-selection within the fixed
candidate family under multiplier perturbation; the bootstrap’s and the
cross-validation’s recovery rates
([`forestsearch_bootstrap_dofuture()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md),
[`forestsearch_Kfold()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_Kfold.md))
measure re-discovery from resampled data with the family rebuilt, which
is stronger and costs a full re-run each time.

## Scope

The coverage result of León and Anderson (2026a) is for a fixed,
enumerated family — forest search’s. DINA and GRF propose their families
from the outcomes, and a bootstrap regenerates them; for those
identifiers the multiplier correction targets the effect of the selected
region conditional on the family proposed, and the results above are
empirical evaluations of how the correction behaves on such families,
not instances of the theorem. The one-sided bounds, the two-subgroup
claim and the declaration calibration that build on this correction are
the subject of León and Anderson (2026b) and of
[`vignette("forestsearch")`](https://larry-leon.github.io/forestsearch/articles/forestsearch.md).

## Where to go next

- [`vignette("forestsearch")`](https://larry-leon.github.io/forestsearch/articles/forestsearch.md):
  the same trial under the identification paper — the search, the three
  identifiers at the $`\max_{\mathcal{A}} N(\varepsilon = 0.20)`$ rule,
  the one-sided field bounds, the two-subgroup claim and the declaration
  calibration.
- The package website’s article *Exploratory subgroup identification
  with forestsearch (the 2024 workflow)*: the full bootstrap at a
  working B with its diagnostics and plots, and K-fold cross-validation.
- [`?forestsearch_bootstrap_dofuture`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md),
  [`?summarize_bootstrap_results`](https://larry-leon.github.io/forestsearch/reference/summarize_bootstrap_results.md),
  [`?mr_estimates_table`](https://larry-leon.github.io/forestsearch/reference/mr_estimates_table.md).

## References

Efron, B. (2014). Estimation and accuracy after model selection.
*Journal of the American Statistical Association* 109(507), 991–1007.

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

    #> R version 4.6.1 (2026-06-24)
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
    #> [1] future_1.80.0      survival_3.8-6     forestsearch_0.4.0
    #>
    #> loaded via a namespace (and not attached):
    #>  [1] gt_1.4.0             sass_0.4.10          generics_0.1.4
    #>  [4] xml2_1.6.0           shape_1.4.6.1        lattice_0.22-9
    #>  [7] listenv_1.1.0        digest_0.6.39        magrittr_2.0.5
    #> [10] grf_2.6.1            evaluate_1.0.5       grid_4.6.1
    #> [13] RColorBrewer_1.1-3   iterators_1.0.14     policytree_1.2.5
    #> [16] mvtnorm_1.4-2        fastmap_1.2.0        foreach_1.5.2
    #> [19] jsonlite_2.0.0       Matrix_1.7-5         glmnet_5.2
    #> [22] scales_1.4.0         weightedsurv_0.1.0   codetools_0.2-20
    #> [25] cli_3.6.6            rlang_1.3.0          litedown_0.11
    #> [28] parallelly_1.48.0    future.apply_1.21.0  commonmark_2.0.0
    #> [31] splines_4.6.1        base64enc_0.1-6      withr_3.0.3
    #> [34] yaml_2.3.12          otel_0.2.0           tools_4.6.1
    #> [37] parallel_4.6.1       doFuture_1.4.0       dplyr_1.2.1
    #> [40] ggplot2_4.0.3        globals_0.20.0       vctrs_0.7.3
    #> [43] R6_2.6.1             lifecycle_1.0.5      randomForest_4.7-1.2
    #> [46] fs_2.1.0             pkgconfig_2.0.3      progressr_1.1.0
    #> [49] pillar_1.11.1        gtable_0.3.6         data.table_1.18.6.1
    #> [52] glue_1.8.1           Rcpp_1.1.2           xfun_0.62
    #> [55] tibble_3.3.1         tidyselect_1.2.1     knitr_1.52
    #> [58] farver_2.1.2         patchwork_1.3.2      htmltools_0.5.9
    #> [61] rmarkdown_2.32       compiler_4.6.1       S7_0.2.2
    #> [64] markdown_2.1
