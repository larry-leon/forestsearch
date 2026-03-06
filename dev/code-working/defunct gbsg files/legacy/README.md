
# forestsearch <img src="man/figures/logo.png" align="right" height="139" alt="" />

<!-- badges: start -->
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/larry-leon/forestsearch/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/larry-leon/forestsearch/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

Exploratory subgroup identification in clinical trials with survival
endpoints using bootstrap bias correction.

**forestsearch** implements the Forest Search methodology described in
Leon et al. (2024),
*"Exploratory subgroup identification in the heterogeneous Cox model,"*
[*Statistics in Medicine*](https://doi.org/10.1002/sim.10163). The
package screens all possible subgroups (formed by combinations of
baseline factors) against Cox hazard ratio thresholds indicative of
harm, selects the subgroup "maximally consistent with harm" via
split-sample evaluation, and provides bias-corrected treatment effect
estimates with valid confidence intervals through an infinitesimal
jackknife bootstrap procedure.

## Overview

In randomized clinical trials, important subgroups based on patient
characteristics may not be anticipated or well understood prior to the
trial. Forest Search addresses the need for exploratory, post-hoc
subgroup identification by:

1. **Constructing candidate factors** via LASSO, Generalized Random
   Forests (GRF), and quartile-based binary splits.
2. **Enumerating** all one- and two-factor subgroup combinations meeting
   minimum sample size and event criteria.
3. **Screening** candidates against a hazard ratio threshold (default
   HR >= 1.25) and evaluating **split-sample consistency** across
   repeated random 50/50 splits (default *R* = 400).
4. **Selecting** the subgroup with the highest consistency rate (>=
   90%) as the estimated detrimental subgroup.
5. **Correcting** for selection bias via bootstrap bias correction with
   infinitesimal jackknife variance estimation.
6. **Validating** stability through N-fold (leave-one-out) and K-fold
   cross-validation.

By reversing the treatment indicator, the same framework identifies
subgroups with substantial **benefit**.

## Installation

You can install the development version of **forestsearch** from
GitHub:

``` r
# install.packages("pak")
pak::pak("larry-leon/forestsearch")
```

Or using `remotes`:

``` r
# install.packages("remotes")
remotes::install_github("larry-leon/forestsearch")
```

The package depends on
[weightedsurv](https://github.com/larry-leon/weightedsurv), which is
installed automatically from GitHub via the `Remotes` field in the
DESCRIPTION.

## Quick Start

``` r
library(forestsearch)

# --- Step 1: Run ForestSearch ---
fs <- forestsearch(
  df.analysis     = trial_data,
  outcome.name    = "tte",
  event.name      = "event",
  treat.name      = "treat",
  confounders.name = c("age", "biomarker", "stage", "region"),
  sg_focus        = "maxSG",
  hr.threshold    = 1.25,
  hr.consistency  = 1.0,
  pconsistency.threshold = 0.90,
  fs.splits       = 400,
  use_lasso       = TRUE,
  use_grf         = TRUE,
  parallel_args   = list(plan = "multisession", workers = 4)
)

print(fs)
summary(fs)

# --- Step 2: Bootstrap Bias Correction ---
fs_bc <- forestsearch_bootstrap_dofuture(
  fs.est    = fs,
  nb_boots  = 2000,
  parallel_args = list(plan = "multisession", workers = 6)
)

# Publication-ready tables
summaries <- summarize_bootstrap_results(fs_bc)
summaries$main_table_gt

# --- Step 3: Cross-Validation ---
# K-fold CV (repeated 10-fold)
fs_cv <- forestsearch_tenfold(
  fs.est  = fs,
  sims    = 200,
  Kfolds  = 10,
  parallel_args = list(plan = "multisession", workers = 6)
)

cv_metrics_tables(fs_cv)

# N-fold (leave-one-out) CV
fs_oob <- forestsearch_Kfold(
  fs.est  = fs,
  Kfolds  = nrow(trial_data),
  parallel_args = list(plan = "multisession", workers = 6)
)

# --- Step 4: Visualization ---
# Forest plot with bias-corrected estimates and CV metrics
plot_subgroup_results_forestplot(
  fs_results  = list(fs.est = fs, fs_bc = fs_bc, fs_kfold = fs_cv),
  df_analysis = trial_data,
  outcome.name = "tte",
  event.name   = "event",
  treat.name   = "treat"
)

# Weighted Kaplan-Meier curves by subgroup
plot_sg_weighted_km(fs, fs_bc = fs_bc)
```

## Key Features

### Core Algorithm

| Feature | Description |
|---|---|
| `forestsearch()` | Main search engine: LASSO/GRF candidate selection, exhaustive enumeration, consistency evaluation |
| `subgroup.consistency()` | Split-sample consistency with optional two-stage sequential early stopping |
| `select_best_subgroup()` | Multiple selection criteria: `"hr"`, `"maxSG"`, `"minSG"`, `"hrMaxSG"`, `"hrMinSG"` |

### Bootstrap Inference

| Feature | Description |
|---|---|
| `forestsearch_bootstrap_dofuture()` | Parallelised bootstrap with `doFuture`; full algorithm mimicked in each replicate |
| Bias correction | Two-source correction on the log-HR scale (Leon et al., Eq. 6--7) |
| IJ variance | Infinitesimal jackknife variance estimation for valid confidence intervals |
| `summarize_bootstrap_results()` | Publication-ready `gt` tables with diagnostics |

### Cross-Validation

| Feature | Description |
|---|---|
| `forestsearch_Kfold()` | N-fold (leave-one-out) CV for algorithmic stability assessment |
| `forestsearch_tenfold()` | Repeated K-fold CV (e.g., 200 x 10-fold) |
| `cv_metrics_tables()` | sensCV, ppvCV, and exact-match correspondence metrics |

### Visualization

| Feature | Description |
|---|---|
| `plot_subgroup_results_forestplot()` | Forest plot with bias-corrected HRs, reference subgroups, and CV annotations |
| `plot_sg_weighted_km()` | Weighted Kaplan-Meier curves for identified subgroups |
| `plot_km_band_forestsearch()` | Multi-subgroup KM comparison bands |
| `plot_spline_treatment_effect()` | Spline-based treatment effect curves |
| `plot_detection_curve()` | Power/detection probability across effect sizes |

### Simulation Framework

| Feature | Description |
|---|---|
| `generate_aft_dgm_flex()` | Flexible AFT data-generating mechanism with configurable heterogeneity |
| `create_gbsg_dgm()` | GBSG-based DGM for replicating paper simulations |
| `run_simulation_analysis()` | Full simulation wrapper: FS, FS+GRF, and standalone GRF |
| `summarize_simulation_results()` | Operating characteristics tables (type-1 error, power, sensitivity, PPV) |

### Multi-Regional Clinical Trials

| Feature | Description |
|---|---|
| `mrct_region_sims()` | Regional consistency evaluation in MRCTs |
| `create_dgm_for_mrct()` | DGM with region-level heterogeneity |

## Methodology at a Glance

Forest Search identifies the subgroup *maximally consistent with harm*
through the following procedure:

```
Candidate Factors ──► Enumerate Subgroups ──► HR Screening (≥ 1.25)
                                                     │
                                              Consistency Evaluation
                                              (R = 400 random 50/50 splits)
                                                     │
                                              Select Ĥ (consistency ≥ 90%)
                                                     │
                                              Bootstrap Bias Correction
                                              (B = 300–2000 replicates)
```

The **bias-corrected estimator** accounts for two sources of optimism:

```
β̂*(Ĥ) = β̂(Ĥ) − (1/B) Σ [η*_b(Ĥ*_b) + η*_b(Ĥ)]
```

where η\*\_b(Ĥ\*\_b) and η\*\_b(Ĥ) capture the discrepancies between
bootstrap and observed data estimates for the bootstrap-selected and
observed subgroups, respectively. Confidence intervals use the
**infinitesimal jackknife** variance estimator.

Simulation studies (20,000 replicates) demonstrate:

- **Type-1 error** ≤ 3% for FS with LASSO (FS_l), even with noise
  factors
- **Power** 71--96% for detecting HR = 2.0 subgroups (N = 300--700)
- **Conservative estimation**: bias-corrected estimates underestimate
  harm and overestimate benefit---a favorable property for exploratory
  analyses
- **Oracle CI coverage** ≥ 95% across all scenarios

## Requirements

- **R** ≥ 4.1.0
- Key dependencies: `survival`, `data.table`, `glmnet`, `grf`,
  `doFuture`, `future`, `ggplot2`, `gt`
- Optional: `forestploter` (forest plots), `DiagrammeR` (diagrams),
  `cubature` (numerical integration)

## Documentation

- **Getting Started vignette**: `vignette("forestsearch")` --- end-to-end
  walkthrough with the GBSG breast cancer trial
- **Methodology vignette**: `vignette("methodology")` --- detailed
  statistical framework, algorithm, simulation results, and applications
- **pkgdown site**: <https://larry-leon.github.io/forestsearch/>
- **Function reference**: organised across 16 topic areas (core
  algorithm, bootstrap, CV, visualization, simulation, MRCT, and more)

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
```

## Related Work

- Athey, S. & Imbens, G. (2016). Recursive partitioning for
  heterogeneous causal effects. *PNAS*.
- Wager, S. & Athey, S. (2018). Estimation and inference of
  heterogeneous treatment effects using random forests. *JASA*.
- Efron, B. (2014). Estimation and accuracy after model selection.
  *JASA*.

## License

MIT © Larry F. Leon
