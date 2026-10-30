
# forestsearch <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

Interpretable subgroup identification in randomized clinical trials, with
post-selection inference for the effect the trial reports.

**forestsearch** identifies a subgroup defined by explicit covariate cuts,
such as "no measured estrogen receptor and progesterone receptor at most
26 fmol/l". It reports the trial's standard analysis for that subgroup and
for the patients outside it, corrected for the selection that produced the
subgroup. Survival outcomes use the Cox hazard ratio. Binary, continuous and
count outcomes use the coefficient of a generalized linear model: an odds
ratio, relative risk, risk difference, mean difference or incidence-rate
ratio.

The methods are described in León and Anderson (2026a, 2026b), and in León
et al. (2024) for the original procedure (see Citation below).

## How it works

1. **Candidates.** The search enumerates every conjunction of at most two
   covariate cuts (`maxk = 2`) that meets minimum size and event counts. The
   candidates come from the covariates alone; no outcome model decides them.
2. **Screen.** Each candidate receives the trial's own analysis within the
   subgroup. A candidate is admitted when its effect clears the screening
   threshold and its splitting-consistency rate reaches `p*`. The rate is the
   probability that both halves of a random half-split clear the
   consistency threshold, and it is evaluated in closed form from one fit
   (`consistency_method = "resample"`, the default).
3. **Select.** A stated rule chooses among the admitted candidates. For
   example, `sg_focus = "effMaxSG"` with `effect_neighborhood = 0.20` takes
   the largest subgroup whose effect is within 20% of the strongest. The
   search can also return no subgroup.
4. **Correct.** With `mr_inference = TRUE`, multiplier resampling re-applies
   the selection rule to jointly perturbed candidate effects, without
   refitting. It reports a de-biased estimate for the subgroup and for its
   complement, a one-sided selection-adjusted bound for each, and a
   Bonferroni pair for a joint claim on both.
5. **Calibrate.** `fs_declaration_calibration()` estimates how often the
   screen, as run, would declare a subgroup if every patient benefited by a
   clinically specified amount. It also gives the cutoff that would hold that
   rate at a nominal level.

Model-based identifiers can propose the candidates instead: DINA
(`subgroup_method = "dina"`) and causal forests (`subgroup_method = "grf"`).
Their regions are re-ranked on the same coefficient and receive the same
inference, conditional on the candidate family they propose. The full
bootstrap correction and the cross-validation of León et al. (2024) remain
available.

Bounds are read by where they lie against the trial's own effect and against
clinically meaningful effect sizes, not as tests at a hazard ratio of one.

## Installation

From CRAN:

``` r
install.packages("forestsearch")
```

The development version from GitHub:

``` r
# install.packages("pak")
pak::pak("larry-leon/forestsearch")
```

## Quick start

This is the forest-search analysis of the German Breast Cancer Study Group
trial in Table 1 of León and Anderson (2026b). `vignette("forestsearch")` runs
it in full, with the two model-based identifiers, and reads every number.

``` r
library(forestsearch)
library(survival)

df <- within(gbsg, {
  id          <- seq_len(nrow(gbsg))
  time_months <- rfstime / 30.4375
  grade3      <- ifelse(grade == "3", 1, 0)
})

fit <- forestsearch(
  df.analysis = df,
  outcome.name = "time_months", event.name = "status",
  treat.name = "hormon", id.name = "id",
  confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
  is.RCT = TRUE, est.scale = "hr", seedit = 8316951,
  # the candidate family
  conf.cont_jcuts = list(er = 10, pgr = 10), conf_force = c("er <= 0", "pgr <= 0"),
  n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
  # the thresholds and the selection rule
  hr.threshold = 1.00, hr.consistency = 1.00, pconsistency.threshold = 0.90,
  stop_threshold = NULL,
  sg_focus = "effMaxSG", selection_rule = "neighborhood", effect_neighborhood = 0.20,
  # selection-adjusted inference
  mr_inference = TRUE,
  mr_inference_args = list(ci_method = "field", draws = 5000L, include_complement = TRUE,
                           field_scale_complement = "selected", field_decompose = TRUE,
                           keep_declaration_field = TRUE,
                           declaration_c0 = c(0.70, 0.75, 0.80, 0.85)),
  parallel_args = list(plan = "sequential")
)

print(fit)
summary(fit)

# Would a screen protecting a uniform benefit (HR 0.70 to 0.85) have declared it?
fs_declaration_calibration(fit, alpha = 0.10, c0 = c(0.70, 0.75, 0.80, 0.85))
```

## Main functions

| Function | Purpose |
|---|---|
| `forestsearch()` | Candidate enumeration, screening, selection and (with `mr_inference = TRUE`) selection-adjusted inference, in one call |
| `fs_mr_inference()` | Multiplier-resampling inference on a supplied candidate family |
| `fs_declaration_calibration()` | Family-wise declaration rate of the screen at a protected level, and the calibrated cutoff |
| `mr_estimates_table()` | The multiplier-resampling estimates in the layout of the bootstrap tables |
| `guohe_algorithm3()` | The de-biased estimate and one-sided bootstrap bound of Guo and He (2021), for comparison |
| `forestsearch_bootstrap_dofuture()`, `summarize_bootstrap_results()` | The full bootstrap bias correction of León et al. (2024) |
| `forestsearch_Kfold()`, `forestsearch_tenfold()`, `cv_metrics_tables()` | Cross-validation of the identification |
| `plot_subgroup_results_forestplot()`, `plot_sg_weighted_km()` | Forest plots and weighted Kaplan-Meier curves |
| `generate_aft_dgm_flex()`, `create_gbsg_dgm()`, `run_simulation_analysis()` | Data-generating mechanisms and the simulation framework |
| `mrct_region_sims()`, `create_dgm_for_mrct()` | Regional consistency in multi-regional trials |

## Documentation

- `vignette("forestsearch")`: the analysis above, end to end.
- The package website, <https://larry-leon.github.io/forestsearch/>, has
  articles on post-selection inference for the same trial (the full
  bootstrap beside multiplier resampling), on the comparison with Guo and He
  (2021), on treatment-effect definitions and their causal reading, on the
  simulation framework, and on the original workflow of León et al. (2024).
- `?forestsearch` documents every argument, with examples for this analysis
  and for the binary and continuous ACTG 175 analyses of León and Anderson
  (2026b).

## Citation

If you use **forestsearch** in your work, please cite:

``` bibtex
@article{Leon2024,
  author  = {Le{\'o}n, Larry F. and Jemielita, Thomas and Guo, Zifang
             and Marceau West, Rachel and Anderson, Keaven M.},
  title   = {Exploratory subgroup identification in the heterogeneous
             {C}ox model: A relatively simple procedure},
  journal = {Statistics in Medicine},
  year    = {2024},
  doi     = {10.1002/sim.10163}
}

@misc{LeonAnderson2026a,
  author        = {Le{\'o}n, Larry F. and Anderson, Keaven M.},
  title         = {Inference for standard trial estimands after data-driven
                   subgroup discovery},
  year          = {2026},
  eprint        = {2609.38361},
  archivePrefix = {arXiv},
  primaryClass  = {stat.ME},
  doi           = {10.48550/arXiv.2609.38361}
}

@misc{LeonAnderson2026b,
  author        = {Le{\'o}n, Larry F. and Anderson, Keaven M.},
  title         = {Interpretable subgroup identification across survival and
                   generalized-linear outcomes},
  year          = {2026},
  eprint        = {2610.28379},
  archivePrefix = {arXiv},
  primaryClass  = {stat.ME},
  doi           = {10.48550/arXiv.2610.28379}
}
```

## Related work

- Guo, X. & He, X. (2021). Inference on selected subgroups in clinical
  trials. *JASA*.
- Athey, S. & Imbens, G. (2016). Recursive partitioning for
  heterogeneous causal effects. *PNAS*.
- Wager, S. & Athey, S. (2018). Estimation and inference of
  heterogeneous treatment effects using random forests. *JASA*.
- Efron, B. (2014). Estimation and accuracy after model selection.
  *JASA*.

## License

MIT © Larry F. Leon
