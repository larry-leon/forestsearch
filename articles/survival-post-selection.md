# Post-selection inference for an identified survival subgroup

``` r

knitr::opts_chunk$set(comment = "#>", fig.align = "center")
library(forestsearch)
library(survival)
```

[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
identifies a subgroup by searching over candidate subgroups, so a
treatment effect estimated within the identified subgroup is affected by
the selection: a confidence interval computed as if the subgroup had
been prespecified tends to under-cover. With `mr_inference = TRUE`,
[`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
reports selection-adjusted estimates and bounds for the identified
subgroup and its complement, based on multiplier resampling of the
candidate effects. The methodology, its assumptions and its simulation
operating characteristics are described in León and Anderson (2026);
this vignette shows how to obtain and read the reported quantities on a
public dataset.

## The data and the fit

The GBSG node-positive breast-cancer trial shipped with **survival**:
686 patients, hormonal therapy against none, recurrence-free survival.
The covariate set is the usual seven; `er` and `pgr` are forced in as
candidate receptor cutpoints, and continuous covariates get 10 quantile
cuts each.

``` r

df.analysis <- gbsg
df.analysis <- within(df.analysis, {
  id          <- seq_len(nrow(df.analysis))
  time_months <- rfstime / 30.4375
  grade3      <- ifelse(grade == "3", 1, 0)
})
confounders.name <- c("age", "meno", "size", "grade3", "nodes", "pgr", "er")
cat("Sample size:", nrow(df.analysis),
    "| events:", sum(df.analysis$status),
    sprintf("(%.1f%%)\n", 100 * mean(df.analysis$status)))
```

    #> Sample size: 686 | events: 299 (43.6%)

The fit below uses the consistency search without a LASSO/GRF/DINA
screen, `sg_focus = "effMaxSG"` (the largest subgroup within 20% of the
strongest observed effect), 200 consistency splits and 500 multiplier
draws. `mr_inference = TRUE` attaches the post-selection results to the
fit as `fs$mr_inference`; `field_recovery = TRUE` adds the
membership-agreement diagnostics used below. Bootstrap and
cross-validation are not run here; see
[`vignette("forestsearch")`](https://larry-leon.github.io/forestsearch/articles/forestsearch.md).

``` r

t0 <- proc.time()
fs <- forestsearch(
  df.analysis,
  outcome.name = "time_months", event.name = "status",
  treat.name = "hormon", id.name = "id",
  confounders.name = confounders.name,
  is.RCT = TRUE, seedit = 8316951, est.scale = "hr",
  use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE,
  conf.cont_jcuts = list(er = 10, pgr = 10), collapse_cuts = TRUE,
  conf_force = c("er <= 0", "pgr <= 0"),
  subgroup_method = "consistency", sg_focus = "effMaxSG",
  selection_rule = "neighborhood", effect_neighborhood = 0.20,
  mr_inference = TRUE,
  mr_inference_args = list(draws = 500L, field_recovery = TRUE),
  n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
  hr.threshold = 1.0, hr.consistency = 1.0, pconsistency.threshold = 0.80,
  fs.splits = 200, max.minutes = 3,
  details = FALSE, quiet = TRUE, plot.sg = FALSE,
  parallel_args = list(plan = "sequential"))
cat(sprintf("fit time: %.1f seconds\n", (proc.time() - t0)[["elapsed"]]))
```

    #> fit time: 32.8 seconds

## The identified subgroup

``` r

cat("Selected subgroup:", paste(fs$sg.harm, collapse = " & "), "\n")
```

    #> Selected subgroup: {er <= 0} & {pgr <= 26}

The subgroup is a conjunction of at most `maxk = 2` binary cutpoints
drawn from the candidate family. `family_status` records how the
candidate family was built; `"no-front-end"` means the family was
enumerated directly rather than proposed by a model fitted to the data.

## The reported quantities

[`print()`](https://rdrr.io/r/base/print.html) reports them; nothing
extra needs to be called.

``` r

print(fs)
```

    #> ForestSearch Results
    #> ====================
    #>
    #> Selected Subgroup:
    #>   Definition: {er <= 0} & {pgr <= 26}
    #>   sg_focus: hrMaxSG
    #>   N: 75
    #>   HR: 2.222
    #>   Pcons: 0.99
    #>   Algorithm: twostage
    #>   Candidate family: no-front-end
    #>   Admission set: effect floor 0; consistency floor 0 at p* = 0.8
    #>   Candidates evaluated: 120
    #>   Candidates passed: 16
    #>
    #> Post-selection inference (selection-adjusted; Leon & Anderson 2026):
    #>   Harm subgroup H:        one-sided 95% lower bound on HR   0.609
    #>   Complement Hc:          one-sided 95% upper bound on HR   0.812   [field-s]
    #>   Joint (Bonferroni):     H lower 0.500, Hc upper 0.848  (gamma 0.025 each side) [field-s]
    #>   Two-sided (IJ, secondary): H (0.549, 2.900)
    #>   Re-selection frequency  p-hat(H) = 0.006
    #>
    #>   One-sided bounds are the primary products; the two-sided IJ interval is
    #>   secondary. p-hat(H) is a descriptive diagnostic. Read each bound by its
    #>   location relative to a clinically meaningful effect size. Methods and
    #>   simulation operating characteristics: Leon & Anderson (2026),
    #>   arXiv:2609.38361.
    #>
    #> Computation time: 0.08 minutes

- **Identified subgroup:** the one-sided 95% lower bound on the hazard
  ratio (the *field* construction).
- **Complement:** the one-sided 95% upper bound on the hazard ratio
  (*field-s*, the studentized complement field), relevant to a benefit
  claim.
- **Both together:** the Bonferroni pair, each side at level 0.025.
- **Two-sided interval:** the infinitesimal-jackknife (IJ) two-term
  interval, reported as a secondary summary.
- **p-hat(H):** the re-selection frequency, a descriptive diagnostic
  (below).

The field bounds require `ci_method = "field"` (the default); with
`ci_method = "ij"` only the two-sided interval is reported. León and
Anderson (2026) describe each construction and report its performance in
simulation.

[`summary()`](https://rdrr.io/r/base/summary.html) adds the standard
errors and the leading re-selection frequencies:

``` r

summary(fs)
```

    #> ForestSearch Summary
    #> ====================
    #>
    #> Analysis Parameters:
    #>   sg_focus: hrMaxSG
    #>   hr.threshold: 1
    #>   hr.consistency: 1
    #>   pconsistency.threshold: 0.8
    #>   n.min: 60
    #>   fs.splits: 200
    #>   maxk: 2
    #>   use_twostage: TRUE
    #>   use_lasso: FALSE
    #>   use_grf: FALSE
    #>   use_dina: FALSE
    #>   candidate family: no-front-end
    #>   admission set: effect floor 0; consistency floor 0 at p* = 0.8
    #>
    #> Variable Selection:
    #>   Candidate confounders: 33
    #>   Confounders evaluated: 33
    #>
    #> Search Space:
    #>   Proportion of max combinations searched: 1
    #>   Maximum subgroup HR estimate: 2.537
    #>
    #> Consistency Evaluation:
    #>   Algorithm: twostage
    #>   Candidates evaluated: 120
    #>   Candidates passed: 16
    #>   Stage 1 screening splits:
    #>   Stage 2 batch size:
    #>
    #> Selected Subgroup:
    #>   Definition: {er <= 0} & {pgr <= 26}
    #>   Sample size: 75
    #>   Hazard ratio: 2.222
    #>   Consistency: 99 %
    #>   Estimation data: n = 686 (harm = 75 , complement = 611 )
    #>
    #> Post-selection inference (selection-adjusted; Leon & Anderson 2026):
    #>   Harm subgroup H:        one-sided 95% lower bound on HR   0.609
    #>   Complement Hc:          one-sided 95% upper bound on HR   0.812   [field-s]
    #>   Joint (Bonferroni):     H lower 0.500, Hc upper 0.848  (gamma 0.025 each side) [field-s]
    #>   Two-sided (IJ, secondary): H (0.549, 2.900)
    #>
    #>   Standard errors (log scale):
    #>     field (H)                 se_field    0.411
    #>     field-s (Hc)              se_field_s  0.135   (naive complement SE 0.134)
    #>     IJ two-term (H)           se_ij       0.425
    #>   Re-selection frequency  p-hat(H) = 0.006
    #>     top-3 re-selection mass:  q1.1 & q18.1 0.182 | q1.1 & q27.1 0.156 | q10.0 & q30.0 0.056
    #>     p_hat_sum = 1.000 over a family of 1744 candidates
    #>     membership recovery: sens_H = 0.559 (mean share of Hhat retained over 998 draws)
    #>
    #>   One-sided bounds are the primary products; the two-sided IJ interval is
    #>   secondary. p-hat(H) is a descriptive diagnostic. Read each bound by its
    #>   location relative to a clinically meaningful effect size. Methods and
    #>   simulation operating characteristics: Leon & Anderson (2026),
    #>   arXiv:2609.38361.
    #>
    #> Computation time: 0.08 minutes

## How to read a bound: by location

These bounds are read **by location against a clinically meaningful
effect size**, never as significance at the null.

The harm-side question is *how much harm the data can rule out*. A
one-sided lower bound in the high 0.8s or above, against a naive
subgroup estimate well above 1, says the adjusted evidence does not
establish clinically meaningful harm – the naive estimate was largely
selection.

The complement-side question is a benefit claim: *is the one-sided upper
bound below a threshold that would matter?* Name the threshold you are
reading against – for example 0.80 or 0.85 – and say whether the bound
clears it. A threshold quoted this way is a reading aid, not a decision
rule.

``` r

g <- fs$mr_inference
cat(sprintf("naive HR on Hhat            %.3f\n", g$naive$est))
```

    #> naive HR on Hhat            2.222

``` r

cat(sprintf("de-biased HR on Hhat        %.3f\n", g$debiased$est))
```

    #> de-biased HR on Hhat        1.262

``` r

cat(sprintf("field one-sided lower       %.3f\n", g$field$lower_1s))
```

    #> field one-sided lower       0.609

``` r

cat(sprintf("field-s one-sided upper Hc  %.3f  (vs 0.80 / 0.85)\n",
            g$field$complement$upper_1s_s))
```

    #> field-s one-sided upper Hc  0.812  (vs 0.80 / 0.85)

The distance between the naive estimate and the de-biased one is what
selection cost; the spread across the constructions is what the
adjustments pay for it.

## p-hat as a diagnostic

`p_hat(Hhat)` is the frequency with which the gate’s own re-selection
map picks the same subgroup across the multiplier draws. It is a
**recorded diagnostic**: it is computed after the bounds are formed, and
**no construction reads it**. It is also estimated from the same draws
that build the bound, so it is a flag, not an independent instrument –
which is a reason not to let any interval depend on it.

``` r

ph  <- g$reselection$p_hat
lab <- g$selected_label
cat(sprintf("p-hat(Hhat) = %.3f over a family of %d candidates\n",
            unname(ph[[lab]]), g$n_family))
```

    #> p-hat(Hhat) = 0.006 over a family of 1744 candidates

``` r

print(round(sort(ph, decreasing = TRUE)[1:5], 3))
```

    #>  q1.1 & q18.1  q1.1 & q27.1 q10.0 & q30.0          q1.1 q30.1 & q32.0
    #>         0.182         0.156         0.056         0.054         0.052

A small p-hat indicates that the exact boundary of the identified
subgroup is not stable under perturbation; by itself it is neither
evidence for nor against the bounds. León and Anderson (2026) discuss
how the bias correction behaves across values of p-hat.

### What p-hat alone cannot say

p-hat counts **exact** re-selections of `Hhat`. Over a large family in
which most candidates differ from `Hhat` by a single cutpoint, a small
exact-match rate is compatible with two very different situations: the
draws re-selected near-twins of `Hhat`, or they re-selected unrelated
regions. p-hat does not distinguish them. **Membership agreement** does,
and it is computable from the same draws. Each draw’s re-selection is
compared with `Hhat` as a set of patients, in the vocabulary of the
cross-validation agreement metrics
([`forestsearch_Kfold()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_Kfold.md)’s
`sens_H` / `ppv_H`).

``` r

rc <- g$field$recovery
cat(sprintf("sens_H  = %.3f   mean share of Hhat retained by a re-selection
", rc$sens_H))
```

    #> sens_H  = 0.559   mean share of Hhat retained by a re-selection

``` r

cat(sprintf("ppv_H   = %.3f   mean share of a re-selection that lies in Hhat
", rc$ppv_H))
```

    #> ppv_H   = 0.563   mean share of a re-selection that lies in Hhat

``` r

cat(sprintf("containment q10 / q50 / q90 = %.3f / %.3f / %.3f
", rc$q10, rc$q50, rc$q90))
```

    #> containment q10 / q50 / q90 = 0.000 / 0.547 / 1.000

``` r

cat(sprintf("share of draws containing all of Hhat = %.3f  (over %d draws)
",
            rc$share_equal_1, rc$n_used))
```

    #> share of draws containing all of Hhat = 0.244  (over 998 draws)

Read together, the two kinds of number say different things about the
same draws: exact re-selection can be rare while the typical
re-selection still retains much of the identified subgroup. That pattern
points to an unstable boundary rather than an unstable region.

One limit. These diagnostics measure re-selection **within the fixed
candidate family**, under multiplier perturbation of the candidate
effects. The bootstrap
([`forestsearch_bootstrap_dofuture()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md))
and cross-validation
([`forestsearch_Kfold()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_Kfold.md))
recovery rates measure something stronger: re-**discovery** of a
subgroup from resampled data, with the search re-run and the family
rebuilt each time. The two are related but neither substitutes for the
other, and this one is free in every analysis while those cost a full
re-run.

## References

León, L. F. and Anderson, K. M. (2026). Inference for standard trial
estimands after data-driven subgroup discovery. arXiv:2609.38361.
<https://arxiv.org/abs/2609.38361>

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
    #>  [1] gt_1.4.0             generics_0.1.4       xml2_1.6.0
    #>  [4] shape_1.4.6.1        lattice_0.22-9       listenv_1.1.0
    #>  [7] digest_0.6.39        magrittr_2.0.5       grf_2.6.1
    #> [10] evaluate_1.0.5       grid_4.6.1           RColorBrewer_1.1-3
    #> [13] iterators_1.0.14     policytree_1.2.5     mvtnorm_1.4-2
    #> [16] fastmap_1.2.0        foreach_1.5.2        jsonlite_2.0.0
    #> [19] Matrix_1.7-5         glmnet_5.1           scales_1.4.0
    #> [22] weightedsurv_0.1.0   codetools_0.2-20     cli_3.6.6
    #> [25] rlang_1.3.0          parallelly_1.48.0    future.apply_1.21.0
    #> [28] splines_4.6.1        yaml_2.3.12          otel_0.2.0
    #> [31] tools_4.6.1          parallel_4.6.1       doFuture_1.4.0
    #> [34] dplyr_1.2.1          ggplot2_4.0.3        globals_0.20.0
    #> [37] vctrs_0.7.3          R6_2.6.1             lifecycle_1.0.5
    #> [40] randomForest_4.7-1.2 fs_2.1.0             pkgconfig_2.0.3
    #> [43] progressr_1.1.0      pillar_1.11.1        gtable_0.3.6
    #> [46] data.table_1.18.6.1  glue_1.8.1           Rcpp_1.1.2
    #> [49] xfun_0.62            tibble_3.3.1         tidyselect_1.2.1
    #> [52] knitr_1.52           farver_2.1.2         patchwork_1.3.2
    #> [55] htmltools_0.5.9      rmarkdown_2.32       compiler_4.6.1
    #> [58] S7_0.2.2
