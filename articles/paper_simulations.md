# Simulation Studies for Evaluating ForestSearch Performance

## 1 Introduction

This vignette demonstrates how to conduct simulation studies to evaluate
the performance of ForestSearch for identifying subgroups with
differential treatment effects. The simulation framework allows you to:

- Generate synthetic clinical trial data with known treatment effect
  heterogeneity
- Evaluate subgroup identification rates (power)
- Assess classification accuracy (sens, spec, PPV, NPV)
- Compare different analysis methods (ForestSearch, GRF)
- Estimate Type I error under null hypothesis
- Track and summarize computational timings across simulations

### 1.1 Simulation Framework Overview

The simulation workflow consists of four main steps:

``` mermaid
flowchart LR
    A[Create DGM] --> B[Simulate Trials]
    B --> C[Run Analyses]
    C --> D[Summarize Results]
```

1.  **Create DGM**: Define a data generating mechanism with specified
    treatment effects
2.  **Simulate Trials**: Generate multiple simulated datasets
3.  **Run Analyses**: Apply ForestSearch (and optionally GRF) to each
    dataset
4.  **Summarize Results**: Aggregate operating characteristics across
    simulations

### 1.2 Key Updates in This Version

The simulation framework has been aligned with
[`generate_aft_dgm_flex()`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md)
methodology:

| Feature | Description |
|----|----|
| **Individual Potential Outcomes** | `theta_0`, `theta_1`, `loghr_po` columns |
| **Average Hazard Ratios (AHR)** | Alternative to Cox-based HR estimation |
| **Stacked PO for Cox HR** | Same epsilon for causal inference |
| **`use_twostage` Parameter** | Faster exploratory analysis option |
| **Backward Compatible** | Works with old and new DGM formats |

## 2 Setup

``` r
# Core packages
library(forestsearch)
library(weightedsurv)

library(data.table)
library(survival)
library(ggplot2)
library(gt)

# Parallel processing
library(foreach)
library(doFuture)
library(future)

# Source simulation functions (if not yet in package)
# source("sim_aft_gbsg_refactored.R")
# source("oc_analyses_gbsg.R")
```

### 2.1 Initializing the Timing Framework

A structured timing framework tracks computational cost at each stage of
the simulation study. We initialize a named list that accumulates
elapsed times for DGM creation, calibration, validation, simulation
loops under H1 and H0, summarization, and table formatting.

``` r
# ── Master timing list ──────────────────────────────────────────────────────
# Each entry stores elapsed wall-clock seconds for one stage.
# Entries are added cumulatively as we proceed through the vignette.
timings <- list()

t_vignette_start <- proc.time()
```

## 3 Creating a Data Generating Mechanism

The simulation framework uses the German Breast Cancer Study Group
(GBSG) dataset as a template for realistic covariate distributions and
censoring patterns.

### 3.1 Understanding the DGM

The
[`create_gbsg_dgm()`](https://larry-leon.github.io/forestsearch/reference/create_gbsg_dgm.md)
function creates a data generating mechanism (DGM) based on an
Accelerated Failure Time (AFT) model with Weibull distribution. Key
features:

- **Covariates**: Age, estrogen receptor, menopausal status,
  progesterone receptor, nodes
- **Treatment effect heterogeneity**: Specified via interaction terms
- **Subgroup definition**: H = {low estrogen receptor AND premenopausal}
- **Censoring**: Weibull or uniform censoring model

#### 3.1.1 New Output Structure (Aligned with generate_aft_dgm_flex)

The DGM now includes:

    dgm$hazard_ratios <- list(
      overall        = hr_causal,      # Cox-based overall HR
      AHR            = AHR,            # Average HR from loghr_po

      AHR_harm       = AHR_H_true,     # AHR in harm subgroup
      AHR_no_harm    = AHR_Hc_true,    # AHR in complement
      harm_subgroup  = hr_H_true,      # Cox-based HR in H
      no_harm_subgroup = hr_Hc_true    # Cox-based HR in Hc
    )

The super-population data (`dgm$df_super_rand`) now contains:

| Column     | Description                                            |
|------------|--------------------------------------------------------|
| `theta_0`  | Log-hazard contribution under control                  |
| `theta_1`  | Log-hazard contribution under treatment                |
| `loghr_po` | Individual causal log hazard ratio (theta_1 - theta_0) |

### 3.2 Alternative Hypothesis (Heterogeneous Treatment Effect)

Under the alternative hypothesis, we create a DGM where the treatment
effect varies across patient subgroups:

``` r
t0 <- proc.time()[3]

# Create DGM with heterogeneous treatment effect
# HR in harm subgroup (H) will be > 1 (treatment harmful)
# HR in complement (H^c) will be < 1 (treatment beneficial)

dgm_alt <- create_gbsg_dgm(
  model = "alt",
  k_treat = 1.0,
  k_inter = 2.0,      # Interaction effect multiplier
  k_z3 = 1.0,
  z1_quantile = 0.25, # ER threshold at 25th percentile
  n_super = 5000,
  cens_type = "weibull",
  seed = 8316951,
  verbose = TRUE
)

# Examine the DGM (print method now shows both HR and AHR metrics)
print(dgm_alt)
```

    GBSG-Based AFT Data Generating Mechanism (Aligned)
    ===================================================

    Model type: alt
    Super-population size: 5000

    Effect Modifiers:
      k_treat: 1
      k_inter: 2
      k_z3: 1

    Hazard Ratios (Cox-based, stacked PO):
      Overall (causal): 0.7331
      Harm subgroup (H): 2.9651
      Complement (Hc): 0.6612
      Ratio HR(H)/HR(Hc): 4.4846

    Average Hazard Ratios (from loghr_po):
      AHR (overall): 0.7431
      AHR_harm (H): 3.8687
      AHR_no_harm (Hc): 0.5848
      Ratio AHR(H)/AHR(Hc): 6.6157

    Subgroup definition: z1 == 1 & z3 == 1 (low ER & premenopausal)
    ER threshold: 8 (quantile = 0.25)
    Subgroup size: 634 (12.7%)
    Analysis variables: v1, v2, v3, v4, v5, v6, v7 

``` r
timings$dgm_creation <- proc.time()[3] - t0
```

#### 3.2.1 Accessing Hazard Ratios (New Aligned Format)

``` r
# Traditional access (backward compatible)
cat("Cox-based HRs:\n")
```

    Cox-based HRs:

``` r
cat("  HR(H):", round(dgm_alt$hr_H_true, 4), "\n")
```

      HR(H): 2.9651 

``` r
cat("  HR(Hc):", round(dgm_alt$hr_Hc_true, 4), "\n")
```

      HR(Hc): 0.6612 

``` r
cat("  HR(overall):", round(dgm_alt$hr_causal, 4), "\n")
```

      HR(overall): 0.7331 

``` r
# New AHR metrics (aligned with generate_aft_dgm_flex)
cat("\nAverage Hazard Ratios (from loghr_po):\n")
```

    Average Hazard Ratios (from loghr_po):

``` r
cat("  AHR(H):", round(dgm_alt$AHR_H_true, 4), "\n")
```

      AHR(H): 3.8687 

``` r
cat("  AHR(Hc):", round(dgm_alt$AHR_Hc_true, 4), "\n")
```

      AHR(Hc): 0.5848 

``` r
cat("  AHR(overall):", round(dgm_alt$AHR, 4), "\n")
```

      AHR(overall): 0.7431 

``` r
# Using hazard_ratios list (unified access)
cat("\nVia hazard_ratios list:\n")
```

    Via hazard_ratios list:

``` r
cat("  harm_subgroup:", round(dgm_alt$hazard_ratios$harm_subgroup, 4), "\n")
```

      harm_subgroup: 2.9651 

``` r
cat("  AHR_harm:", round(dgm_alt$hazard_ratios$AHR_harm, 4), "\n")
```

      AHR_harm: 3.8687 

#### 3.2.2 Examining Individual-Level Treatment Effects

``` r
# The super-population now includes individual log hazard ratios
df_super <- dgm_alt$df_super_rand

cat("Individual-level potential outcomes:\n")
```

    Individual-level potential outcomes:

``` r
cat("  theta_0 (control log-hazard) range:", 
    round(range(df_super$theta_0), 3), "\n")
```

      theta_0 (control log-hazard) range: -0.891 1.76 

``` r
cat("  theta_1 (treatment log-hazard) range:", 
    round(range(df_super$theta_1), 3), "\n")
```

      theta_1 (treatment log-hazard) range: -1.427 2.909 

``` r
cat("  loghr_po (individual log-HR) range:", 
    round(range(df_super$loghr_po), 3), "\n")
```

      loghr_po (individual log-HR) range: -0.537 1.353 

``` r
# Verify AHR calculation
ahr_manual <- exp(mean(df_super$loghr_po))
cat("\nAHR verification:\n")
```

    AHR verification:

``` r
cat("  exp(mean(loghr_po)) =", round(ahr_manual, 4), "\n")
```

      exp(mean(loghr_po)) = 0.7431 

``` r
cat("  dgm$AHR =", round(dgm_alt$AHR, 4), "\n")
```

      dgm$AHR = 0.7431 

``` r
# Distribution of individual treatment effects
cat("\nIndividual HR distribution:\n")
```

    Individual HR distribution:

``` r
individual_hr <- exp(df_super$loghr_po)
cat("  Mean:", round(mean(individual_hr), 4), "\n")
```

      Mean: 1.0012 

``` r
cat("  Median:", round(median(individual_hr), 4), "\n")
```

      Median: 0.5848 

``` r
cat("  SD:", round(sd(individual_hr), 4), "\n")
```

      SD: 1.0928 

### 3.3 Calibrating for a Target Hazard Ratio

Often, you want to specify the exact hazard ratio in the harm subgroup.
Use
[`calibrate_k_inter()`](https://larry-leon.github.io/forestsearch/reference/calibrate_k_inter.md)
to find the interaction parameter that achieves this.

#### 3.3.1 Calibrate to Cox-based HR (Default)

``` r
t0 <- proc.time()[3]

# Find k_inter for Cox-based HR = 2.0 in harm subgroup
k_inter_cox <- calibrate_k_inter(
  target_hr_harm = 2.0,
  model = "alt",
  k_treat = 1.0,
  cens_type = "weibull",
  use_ahr = FALSE,  # Default: calibrate to Cox-based HR
  verbose = TRUE
)

# Create DGM with calibrated k_inter
dgm_calibrated_cox <- create_gbsg_dgm(
  model = "alt",
  k_treat = 1.0,
  k_inter = k_inter_cox,
  verbose = TRUE
)

cat("\nVerification (Cox-based):\n")
```

    Verification (Cox-based):

``` r
cat("Achieved HR(H):", round(dgm_calibrated_cox$hr_H_true, 3), "\n")
```

    Achieved HR(H): 2 

``` r
cat("HR(H^c):", round(dgm_calibrated_cox$hr_Hc_true, 3), "\n")
```

    HR(H^c): 0.661 

``` r
cat("Overall HR:", round(dgm_calibrated_cox$hr_causal, 3), "\n")
```

    Overall HR: 0.722 

``` r
timings$calibrate_cox <- proc.time()[3] - t0
```

#### 3.3.2 Calibrate to AHR (New Option)

``` r
t0 <- proc.time()[3]

# Alternatively, calibrate to Average Hazard Ratio
k_inter_ahr <- calibrate_k_inter(
  target_hr_harm = 2.0,
  model = "alt",
  k_treat = 1.0,
  cens_type = "weibull",
  use_ahr = TRUE,  # NEW: calibrate to AHR instead
  verbose = TRUE
)

# Create DGM with AHR-calibrated k_inter
dgm_calibrated_ahr <- create_gbsg_dgm(
  model = "alt",
  k_treat = 1.0,
  k_inter = k_inter_ahr,
  verbose = TRUE
)

cat("\nVerification (AHR-based):\n")
```

    Verification (AHR-based):

``` r
cat("Achieved AHR(H):", round(dgm_calibrated_ahr$AHR_H_true, 3), "\n")
```

    Achieved AHR(H): 2 

``` r
cat("AHR(H^c):", round(dgm_calibrated_ahr$AHR_Hc_true, 3), "\n")
```

    AHR(H^c): 0.585 

``` r
cat("Overall AHR:", round(dgm_calibrated_ahr$AHR, 3), "\n")
```

    Overall AHR: 0.683 

``` r
timings$calibrate_ahr <- proc.time()[3] - t0
```

#### 3.3.3 Compare Cox HR vs AHR Calibration

``` r
# Compare the two calibration approaches
cat("Comparison of calibration methods:\n")
```

    Comparison of calibration methods:

``` r
cat(sprintf("%-20s %-12s %-12s\n", "Metric", "Cox-calib", "AHR-calib"))
```

    Metric               Cox-calib    AHR-calib   

``` r
cat(sprintf("%-20s %-12.4f %-12.4f\n", "k_inter", k_inter_cox, k_inter_ahr))
```

    k_inter              1.4947       1.3016      

``` r
cat(sprintf("%-20s %-12.4f %-12.4f\n", "HR(H)", 
            dgm_calibrated_cox$hr_H_true, dgm_calibrated_ahr$hr_H_true))
```

    HR(H)                2.0000       1.7233      

``` r
cat(sprintf("%-20s %-12.4f %-12.4f\n", "AHR(H)", 
            dgm_calibrated_cox$AHR_H_true, dgm_calibrated_ahr$AHR_H_true))
```

    AHR(H)               2.4001       2.0000      

### 3.4 Validating k_inter Effect on Heterogeneity

Use
[`validate_k_inter_effect()`](https://larry-leon.github.io/forestsearch/reference/validate_k_inter_effect.md)
to verify the interaction parameter properly modulates treatment effect
heterogeneity:

``` r
t0 <- proc.time()[3]

# Test k_inter effect on HR heterogeneity
# k_inter = 0 should give ratio ~ 1 (no heterogeneity)
validation_results <- validate_k_inter_effect(
  k_inter_values = c(-2, -1, 0, 1, 2, 3),
  verbose = TRUE
)
```

    Testing k_inter effect on HR heterogeneity...

    k_inter  HR(H)    HR(Hc)   AHR(H)   AHR(Hc)  Ratio(Cox) Ratio(AHR)
    ----------------------------------------------------------------------
    -2.0     0.1336   0.6612   0.0884   0.5848   0.2021     0.1512
    -1.0     0.3033   0.6612   0.2274   0.5848   0.4587     0.3888
    0.0      0.6552   0.6612   0.5848   0.5848   0.9909     1.0000
    1.0      1.3873   0.6612   1.5041   0.5848   2.0982     2.5721
    2.0      2.9651   0.6612   3.8687   0.5848   4.4846     6.6157
    3.0      6.6375   0.6612   9.9507   0.5848   10.0387    17.0162

    PASS: k_inter = 0 gives Cox ratio ~= 1 (no heterogeneity)
    PASS: k_inter = 0 gives AHR ratio ~= 1 (no heterogeneity)

    AHR vs Cox HR alignment:
      k_inter = -2.0: HR(H) vs AHR(H) diff = 0.0452
      k_inter = -1.0: HR(H) vs AHR(H) diff = 0.0759
      k_inter = 0.0: HR(H) vs AHR(H) diff = 0.0704
      k_inter = 1.0: HR(H) vs AHR(H) diff = 0.1168
      k_inter = 2.0: HR(H) vs AHR(H) diff = 0.9036
      k_inter = 3.0: HR(H) vs AHR(H) diff = 3.3132

``` r
timings$validation <- proc.time()[3] - t0
```

### 3.5 Null Hypothesis (Uniform Treatment Effect)

For Type I error evaluation, create a DGM with uniform treatment effect:

``` r
t0 <- proc.time()[3]

# Create null DGM (no treatment effect heterogeneity)
dgm_null <- create_gbsg_dgm(
  model = "null",
  k_treat = 1.0,
  verbose = TRUE
)

cat("\nNull hypothesis HRs:\n")
```

    Null hypothesis HRs:

``` r
cat("Overall HR:", round(dgm_null$hr_causal, 3), "\n")
```

    Overall HR: 0.722 

``` r
cat("HR(H^c):", round(dgm_null$hr_Hc_true, 3), "\n")
```

    HR(H^c): 0.722 

``` r
cat("AHR(H^c):", round(dgm_null$AHR_Hc_true, 3), "\n")
```

    AHR(H^c): 0.654 

``` r
cat("AHR:", round(dgm_null$AHR, 3), "\n")
```

    AHR: 0.654 

``` r
timings$dgm_null <- proc.time()[3] - t0
```

## 4 Simulating Trial Data

### 4.1 Single Trial Simulation

Use
[`simulate_from_gbsg_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_gbsg_dgm.md)
to generate a single simulated trial:

``` r
# Use the Cox-calibrated DGM for simulations
dgm_calibrated <- dgm_calibrated_cox

# Simulate a single trial
sim_data <- simulate_from_gbsg_dgm(
  dgm = dgm_calibrated,
  n = 700,
  rand_ratio = 1,        # 1:1 randomization
  sim_id = 1,
  max_follow = 84,       # 84 months administrative censoring
  muC_adj = log(1.5)     # Censoring adjustment
)

# Examine the data
cat("Simulated trial:\n")
```

    Simulated trial:

``` r
cat("  N =", nrow(sim_data), "\n")
```

      N = 700 

``` r
cat("  Events =", sum(sim_data$event.sim), 
    "(", round(100 * mean(sim_data$event.sim), 1), "%)\n")
```

      Events = 376 ( 53.7 %)

``` r
cat("  Harm subgroup size =", sum(sim_data$flag.harm),
    "(", round(100 * mean(sim_data$flag.harm), 1), "%)\n")
```

      Harm subgroup size = 86 ( 12.3 %)

``` r
# Quick survival analysis
fit_itt <- coxph(Surv(y.sim, event.sim) ~ treat, data = sim_data)
cat("  Estimated ITT HR =", round(exp(coef(fit_itt)), 3), "\n")
```

      Estimated ITT HR = 0.741 

#### 4.1.1 Examining Individual-Level Effects in Simulated Data

``` r
# The simulated data now includes loghr_po
if ("loghr_po" %in% names(sim_data)) {
  cat("\nIndividual treatment effects in simulated trial:\n")
  
  # Compute AHR in simulated data by subgroup
  ahr_H_sim <- exp(mean(sim_data$loghr_po[sim_data$flag.harm == 1]))
  ahr_Hc_sim <- exp(mean(sim_data$loghr_po[sim_data$flag.harm == 0]))
  ahr_overall_sim <- exp(mean(sim_data$loghr_po))
  
  cat("  AHR(H) in sim:", round(ahr_H_sim, 3), "\n")
  cat("  AHR(Hc) in sim:", round(ahr_Hc_sim, 3), "\n")
  cat("  AHR(overall) in sim:", round(ahr_overall_sim, 3), "\n")
} else {
  cat("\nNote: loghr_po not available in simulated data\n")
}
```

    Individual treatment effects in simulated trial:
      AHR(H) in sim: 2.4
      AHR(Hc) in sim: 0.585
      AHR(overall) in sim: 0.696 

### 4.2 Examining Covariate Structure

``` r
dfcount <- df_counting(
  df = sim_data,
  by.risk = 6,
  tte.name = "y.sim", 
  event.name = "event.sim", 
  treat.name = "treat"
)
plot_weighted_km(dfcount, conf.int = TRUE, show.logrank = TRUE, 
                 ymax = 1.05, xmed.fraction = 0.775, ymed.offset = 0.125)
```

![](data:image/png;base64,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AggggAACCCCAAAIIxFaAAGhsDy07hgACCCCAAAIIIIAAAggggAACCCCAAAIEQDkHEEAAAQQQQAABBBBAAAEEEEAAAQQQQCC2AgRAY3to2TEEEEAAAQQQQAABBBBAAAEEEEAAAQQQIADKOYAAAggggAACCCCAAAIIIIAAAggggAACsRUgABrbQ8uOIYAAAggggAACCCCAAAIIIIAAAggggAABUM4BBBBAAAEEEEAAAQQQQAABBBBAAAEEEIitAAHQ2B5adgwBBBBAAAEEEEAAAQQQQAABBBBAAAEECIByDiCAAAIIIIAAAggggAACCCCAAAIIIIBAbAUIgMb20LJjCCCAAAIIIIAAAggggAACCCCAAAIIIEAAlHMAAQQQQAABBBBAAAEEEEAAAQQQQAABBGIrQAA0toeWHUMAAQQQQAABBBBAAAEEEEAAAQQQQAABAqCcAwgggAACCCCAAAIIIIAAAggggAACCCAQWwECoLE9tOwYAggggAACCCCAAAIIIIAAAggggAACCBAA5RxAAAEEEEAAAQQQQAABBBBAAAEEEEAAgdgKEACN7aFlxxBAAAEEEEAAAQQQQAABBBBAAAEEEECAACjnAAIIIIAAAggggAACCCCAAAIIIIAAAgjEVoAAaGwPLTuGAAIIIIAAAggggAACCCCAAAIIIIAAAgRAOQcQQAABBBBAAAEEEEAAAQQQQAABBBBAILYCBEBje2jZMQQQQAABBBBAAAEEEEAAAQQQQAABBBAgAMo5gAACCCCAAAIIIIAAAggggAACCCCAAAKxFSAAGttDy44hgAACCCCAAAIIIIAAAggggAACCCCAAAFQzgEEEEAAAQQQQAABBBBAAAEEEEAAAQQQiK0AAdDYHlp2DAEEEEAAAQQQQAABBBBAAAEEEEAAAQQIgHIOIIAAAggggAACCCCAAAIIIIAAAggggEBsBQiAxvbQsmMIIIAAAggggAACCCCAAAIIIIAAAgggQACUcwABBBBAAAEEEEAAAQQQQAABBBBAAAEEYitAADS2h5YdQwABBBBAAAEEEEAAAQQQQAABBBBAAAECoJwDCCCAAAIIIIAAAggggAACCCCAAAIIIBBbAQKgsT207BgCCCCAAAIIIIAAAggggAACCCCAAAIIEADlHEAAAQQQQAABBBBAAAEEEEAAAQQQQACB2AoQAI3toWXHEEAAAQQQQAABBBBAAAEEEEAAAQQQQIAAKOcAAggggAACCCCAAAIIIIAAAggggAACCMRWgABobA8tO4YAAggggAACCCCAAAIIIIAAAggggAACBEA5BxBAAAEEEEAAAQQQQAABBBBAAAEEEEAgtgIEQGN7aNkxBBBAAAEEEEAAAQQQQAABBBBAAAEEECAAyjmAAAIIIIAAAggggAACCCCAAAIIIIAAArEVIAAa20PLjiGAAAIIIIAAAggggAACCCCAAAIIIIAAAVDOAQQQQAABBBBAAAEEEEAAAQQQQAABBBCIrQAB0NgeWnYMAQQQQAABBBBAAAEEEEAAAQQQQAABBAiAcg4ggAACCCCAAAIIIIAAAggggAACCCCAQGwFCIDG9tCyYwgggAACCCCAAAIIIIAAAggggAACCCBAAJRzAAEEEEAAAQQQQAABBBBAAAEEEEAAAQRiK0AANLaHlh1DAAEEEEAAAQQQQAABBBBAAAEEEEAAAQKgnAMIIIAAAggggAACCCCAAAIIIIAAAgggEFsBAqCxPbTsGAIIIIAAAggggAACCCCAAAIIIIAAAggQAOUcQAABBBBAAAEEEEAAAQQQQAABBBBAAIHYChAAje2hZccQQAABBBBAAAEEEEAAAQQQQAABBBBAgAAo5wACCCCAAAIIIIAAAggggAACCCCAAAIIxFaAAGhsDy07hgACCCCAAAIIIIAAAggggAACCCCAAAIEQDkHEEAAAQQQQAABBBBAAAEEEEAAAQQQQCC2AgRAY3to2TEEEEAAAQQQQAABBBBAAAEEEEAAAQQQIADKOYAAAggggAACCCCAAAIIIIAAAggggAACsRUgABrbQ8uOIYAAAggggAACCCCAAAIIIIAAAggggAABUM4BBBBAAAEEEEAAAQQQQAABBBBAAAEEEIitAAHQ2B5adgwBBBBAAAEEEEAAAQQQQAABBBBAAAEECIByDiCAAAIIIIAAAggggAACCCCAAAIIIIBAbAUIgMb20LJjCCCAAAIIIIAAAggggAACCCCAAAIIIEAAlHMAAQQQQAABBBBAAAEEEEAAAQQQQAABBGIrQAA0toeWHUMAAQQQQAABBBBAAAEEEEAAAQQQQAABAqCcAwgggAACCCCAAAIIIIAAAggggAACCCAQWwECoLE9tOwYAggggAACCCCAAAIIIIAAAggggAACCBAA5RxAAAEEEEAAAQQQQAABBBBAAAEEEEAAgdgKEACN7aFlxxBAAAEEEEAAAQQQQAABBBBAAAEEEECAACjnAAIIIIAAAggggAACCCCAAAIIIIAAAgjEVoAAaGwPLTuGAAIIIIAAAggggAACCCCAAAIIIIAAAgRAOQcQQAABBBBAAAEEEEAAAQQQQAABBBBAILYCBEBje2jZMQQQQAABBBBAAAEEEEAAAQQQQAABBBAgAMo5gAACCCCAAAIIIIAAAggggAACCCCAAAKxFSAAGttDy44hgAACCCCAAAIIIIAAAggggAACCCCAAAFQzgEEEEAAAQQQQAABBBBAAAEEEEAAAQQQiK0AAdDYHlp2DAEEEEAAAQQQQAABBBBAAAEEEEAAAQQIgHIOIIAAAggggAACCCCAAAIIIIAAAggggEBsBQiAxvbQsmMIIIAAAggggAACCCCAAAIIIIAAAgggQACUcwABBBBAAAEEEEAAAQQQQAABBBBAAAEEYitAADS2h5YdQwABBBBAAAEEEEAAAQQQQAABBBBAAIEqCBBoToEff/zRbrjhBquvr2/O1bIuBBBAAAEEEEAAAQQQQAABBBBAYLEU6Nixo5155pnWtm3bxXL7W8JGEwBtCUdhCdoGBT+vvPLKJWiP2VUEEEAAAQQQQAABBBBAAAEEEECgcQJrrLGGHXDAAY1byBI8NwHQJfjgL4pdD5mf/fr1s2233XZRbALrRAABBBBAAAEEEEAAAQQQQAABBBYLgXvuucc+/PBDatI28mgRAG0kILOXJqDg54ABA0qbmbkQQAABBBBAAAEEEEAAAQQQQACBJUBAwU+9KI0ToBOkxvkxNwIIIIAAAggggAACCCCAAAIIIIAAAgi0YAECoC344LBpCCCAAAIIIIAAAggggAACCCCAAAIIINA4AQKgjfNjbgQQQAABBBBAAAEEEEAAAQQQQAABBBBowQIEQFvwwWHTEEAAAQQQQAABBBBAAAEEEEAAAQQQQKBxAgRAG+fH3AgggAACCCCAAAIIIIAAAggggAACCCDQggUIgLbgg8OmIYAAAggggAACCCCAAAIIIIAAAggggEDjBAiANs6PuRFAAAEEEEAAAQQQQAABBBBAAAEEEECgBQsQAG3BB4dNQwABBBBAAAEEEEAAAQQQQAABBBBAAIHGCRAAbZwfcyOAAAIIIIAAAggggAACCCCAAAIIIIBACxYgANqCDw6bhgACCCCAAAIIIIAAAggggAACCCCAAAKNEyAA2jg/5kYAAQQQQAABBBBAAAEEEEAAAQQQQACBFixAALQFHxw2DQEEEEAAAQQQQAABBBBAAAEEEEAAAQQaJ0AAtHF+zI0AAggggAACCCCAAAIIIIAAAggggAACLViAAGgLPjhsGgIIIIAAAggggAACCCCAAAIIIIAAAgg0ToAAaOP8mBsBBBBAAAEEEEAAAQQQQAABBBBAAAEEWrAAAdAWfHDYNAQQQAABBBBAAAEEEEAAAQQQQAABBBBonAAB0Mb5MTcCCCCAAAIIIIAAAggggAACCCCAAAIItGCBqha8bWwaAotc4LHHHrOnn37adtttN9t7770X+fYsig044YQTbM6cOfaPf/zD2rVrV5ZN+H//7//ZaaedZi+//LKttNJKGZc5e/ZsO/HEE22ttdayU045JeM05RxY7H5+9tlndv3119tvfvMbGzBgQDk3JWVZmdYzf/58++1vf2trrLGG/fOf/0yZPi5fnnrqKXvttdds/PjxttFGG9m2227rrUvZv1KWVco8pWwb8yCAAAIIIIAAAi1d4NZbb7V33nnHdtppJ9t///1zbu7f/vY3++mnn+yII46wLbfc0k87cuRIGzhwYNb5qqurrUOHDtarVy9/jbvyyitnnbY5R/zyyy9255132tdff21LL720358ddtjB2rZtW9BmfPXVV/bmm28WNO32229vffr0yTntX//6Vxs1apT93//9n8ms2DJ37lzvu+6669p1112XMntj9zVlYe6L9v3VV1+1jz76yFq3bm3rrbeeHXzwwVZbW5s+aaO+l+uaXfeejzzyiH388cc2duxYW3XVVe2QQw6xFVZYoaDty3fv+sQTT9jJJ59sL7zwgl92QQtlovILJCgINKOACxQl3FmccMG0Zlxr6as677zz/Paee+65pS9kMZ+zrq7OG7gfgrLsiQtoJbp165Y4+uijcy7vsMMO8+t1P5Q5pyvXyGL30/2I+e3beeedy7UJGZeTbT2PPvqoX7/7Ec043+I6cMaMGQl3ce33TX8rwquioiLhgr1F7VYpyyplnqI2iokRQAABBBBAAIHFTEDX47omO/XUU/NuuQt0+Wld0DQ5rXugn7ymC9d22d5ramoSd999d3LeRfVB19ou+aPBdrvgYWL06NEFbdbNN9/cYP5s+33vvffmXOaNN96YXJauV0spl1xyScIFIBPDhg1Lmb0c+xpdoO71q6qqktsb9tkFeBMumBydtOTP5bxm//TTTxMu6abB9mofbrnlloK2Md+967x58xLrr79+Yquttkq4ZJaClhmdyAVj/fa1hP83otu1uH2mCrz7v5GCAALNJ6DMzylTptiFF16YcaXKbtTTzbvuuivj+KYauN1229mOO+5Y0tPUptqmXMvt16+fbb755nbMMcfY9OnTc026WI3T/rz44ou24YYb+kwDPYG97777rH379qYs3WLOi1KWVco8ixUwG4sAAggggAACCCxCAdXw+d///pfy+u9//2u33367KbtSNc8OPfRQcwHBRbaVgwcPtoMOOshmzpxpV111lSmDVbWy9ttvP3PBMl8zqb6+Pu/2qSabsv6yvbbZZhu/DGW/brzxxlmXpwzU008/Pev4QkYoI1PZuX/605+sd+/eyVnKta9hgfJygXLr1KmTKXNYmaVffPGF9e3b14YOHWr77ruvP8Zh+lLfy3XNPnXqVF/bUw66t3rrrbdswoQJPlNTTscdd5y98cYbWTez0HtXl8xhl112mc8IVgYvZREJLG4RW7Z38RYgA3TxO37FZkbm2sPPP//cP7k6/vjjM07mfhwTm266acrTt+bKAM24QTkGZsvMzDFLSaNyrcdV+fBWl19+eUnLbmkz6fi3atUq4aoVJSZOnJiyea65BL+vrtp/yvBsX0pZVinzZFs/wxFAAAEEEEAAgbgIlDMD1AU4c7L84Q9/8Nd8rgmynNM15UhXzd9vg2oDRouy+FzShB93//33R0cV/dkFBhM9e/b0176PP/541vnl5ZqDSrRp08av14WNEqVkgP7+979PVFZWJn7++eeUdZVzX13AOLHMMsv4fdJ9SrRMmzYt4YKifh9cskN0VNGfy3nNfuaZZ/ptUgZo+rnpAt8+C9g1f5Ax61fbUey9q45l165dEy6Bpaj9JgO0KK6sE5MBuogCz6w2HgJqR0VP5F555RX/ZLDQvXJ/LP1TRPcj6mfRctRuiJ4ghaJx0WFffvmlff/992F0yrueXH3wwQf2+uuvm/tRSxkX/aLl6alqKO6HyD/R+vbbb839lQiDC37X8vTS9hdS9ERQxV3YNJhc7a26agH27rvv2iabbGJnn312g2lKGaD9HTJkiOnJ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URX5nQ/5o62w2jauen2NJeanjte0CozeebHb3m4LOm+Kzp/r86obzrKt98pf3T3XMtLHqX3RkVP1WlhVXx0mrdqlnf3GtRXazvVaT0EAAQQQQAABBBBAAAEEEEAAAQSKFSCiUKxYmaYfP368bbzxxsmljRgxwpZddtnkd31QNqeyQJ977jl7//33bfvtt7err77atthiC9P8jzzyiJ1xxhm++rXaqLzsssvoTCZFcNF+UbA6BKznux6LlBGqTM7wKiUouvPOO9tOO+1k33zzjR122GE2eeIEe/XJR7PuaP3sWTbw/lVcsLThJLec39XeebZdzp7mG861YMi7z9XZjCkLqrwruKry5Xu1rtOlBeP1b2VVwnquVG+uT6iMpW2H+S44m9puavqEc+cnbMi4qf5VV11perWtrnKvXz+7oOiCYfpeZZX0MJ9OyHcEEEAAAQQQQAABBBBAAAEElngBAqAt/BT4z3/+YwcddJDP/lTG6I477ugzRdXmZLQoMNq3b9/oID63IAEFqNWsgV6hRIOi06ZN850qpR/XMG30XVmmq622ml177bU+EB4dFz6//+GH9pzLDn77mcdTfB6oAABAAElEQVRs3W12tE1+t0cYlXw/8IyJ9hvXc7vr3L3g8taTbV1bpG1cr/UV9sA1XVLm++jltqZXMeWcW0fbqhvMLmiWGfXzTK9xFomyps2pwGi3trXWvV1r/+roqtBTEEAAAQQQQAABBBBAAAEEEEBgyRYgANrCj7+qUz/77LN2/vnn20033eSrVEeDZOr1XZmfVH1v4Qcyw+ZFg6LdunXzUygrdPLkyf41depUU5A0W9l8882zjbINNtjAvvj8c/vpp5/sP5dfYG07dHLBzq1Tpu/QZb5tuUdx1dg7dp1nY0ZW26xpC6vTz5nVygdEq2tdJ081C6Oparpz6eXmWuu6hcPCBowcWu2r319/2tLWeZl5NmDQGOvSPUOaapihwPfpLkA6bNIM/9IsrasqrKtrU5R2RQsEZDIEEEAAAQQQQAABBBBAAAEEYijQKuFKDPcrtrv0/fff+w6R2rRpY6uuuqr17t3bFEhbXIoCcx9//LH169fPHn00e9XtxWV/mnI7FfxUZqjaEh07dqxv97OY9T300EN21VVX+VlWWms9W2PTLf25snW//W2pnssXs6ic095zaWd7+cEOdvA5E2zHA6bm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3SrwyvMnQ+QoQv2hH1QWJt+mTRtJSEgQq7OVp2vp0aOHzovZrVs3T82qvA/FnFDrxES9VtYBzo/2GI87wzWhBoqzaIjtiJDt1KmTu0N13126dJGmTZu6bWN2QJBu3769dOzY0WxyOzefi2mANAzG8xjXPmTIEIeiVikpKTpk3rSvyrxPnz5yyimnVOUQtiUBnxAIUhXLynzREyrBw4vTJGru37+/IHku/vCTk5MFbwV27dol+I94zpw5OpEwKovB7rnnHnnwwQd9MQz24ecEUCQLb5HOO+88+fTTT/18tBweCZAACZAACZAACZAACZAACQQeATgamWLFgTd6jtgfCECnmTdvntZqTj/9dIGX6MKFC2s8NOSKXbBgQY378WUH0KwgzPurXX755fLOO+/I22+/LWPGjPHXYfr9uHziAYqrxJsReH+i6lpBQYH8/vvv2sPTEICXp6tCR8gRgj8sGgmQAAmQAAmQAAmQAAmQAAmQAAmQQM0JwAEJHpSoBE4jAV8SGDVqlCQlJVW5SzjLvfDCCy51oSp3xgNIoBoEfCaA4twTJ07UrtFQp+HlZzVn8RMVwyB8QjANCwuzNuUyCZAACZAACZAACZAACZAACZAACZBADQggTJoCaA0ANvJDTQi88Yw086FDh1ariDW8kiGAlpaWyo4dO3R4fiNHzMuvYwI+FUAx9l69emnvT4Q3r1mzRtauXaun9PR0Qa6Krl27CvJiXHbZZYK8HjQSIAESIAESIAESIAESIAESIAESIAHfEoiLi/Nth+ytURGYNGmSjBw5Utq2bauv++abbxZ4f5r16sKAJ+gll1wic+fOlRYtWlS3Gx5HAlUm4HMBFCOAR+fFF19c5cHwABIgARIgARIgARIgARIgARIgARIggZoTMMV7cnJyxExIV0cjAW8III2hVex0XvemD2sbpGWAgApnOXwPMzMzBYWannrqKRk4cKCceOKJ1uZcJgGfE6gVAdTno2SHJEACJEACJEACJEACJEACJEACJEACVSIAgQmTMaSmy83NtQuiEEYpiho6nNc2gSlTpuiI4S1btuhToZg2BFGsUwCtbfrsv9YEUOR1WLZsmaxbt0727t0rJSUl0rp1a+nYsaMMHjxYwsPDSZ8ESIAESIAESIAESIAESIAESIAESKCOCMCLz5Momp2dLWlpaXU0Gp6mMRO49tprpaysTCOAfkQjgdom4HMBFLk+H3roIXnnnXdk3759LseP/3DPPfdc3Q6CKI0ESIAESIAESIAESIAESIAESIAESKDuCTiLovAITUlJkf379+uCNXU/Ip6xIRNAbRh4fEJsp5FAXRII9uXJFi9erIscIYeDO/ET50Ouh7ffflt69uwpr7/+ui+HwL5IgARIgARIgARIgARIgARIgARIgASqSQC5Q+Go1KdPH50DEgIpjQR8RWD69Okyb948+frrr33VJfshAa8I+Ox/stWrV8uIESMkIyNDnxj/aQ4fPlw6d+4sHTp0EPynuWPHDtm+fbt89913kpqaqnOPjB8/XhISEuT888/3asBsRAIkQAIkQAIkQAIkQAIkQAIkQAIkULsE8AyPIjhJSUly4MAB7eSEnKE0EqgOgfz8fJk1a5ZkZWXpw61h78gFisrwTZs2leTkZAkKCtKpE5E+kUYCviLgMwH09ttvt4ufyOVw3333OVQMsw4Y/2m+9NJLcs899+iEy2PHjpWhQ4dKXFyctRmXSYAESIAESIAESIAESIAESIAESIAE6pFAcHCwJCYm6gnP8oj2hCBqFbDqcXg8dYAQgNMcIoFdGdIubN26Ve9asWKFnq9du1aLoK7acxsJVIeATwTQDRs2aK9ODABi5quvvupxLDExMQLBtGXLlnLVVVdp4RRvAlARjEYCJEACJEACJEACJEACJEACJEACJOB/BPAs36lTJx3liahOFDyGZx+NBCojcMwxx8iMGTPsHqBov3z5cq0lwfMT0cQoinTbbbfpCOK+fftW1qXe/95778n8+fPl/9k7E7g7pvv/f20RW0KIBCEiCCGEJChKa99rp7WkaqcqtOpHil9F1Vqp8iN2EqWWEmvttdUSe+yxJREhKRIklizzn8/M/9xn7r1zl5nnPs8z9973eb2e586cOefM+b7nzJkz3znn+x05cqR17dq1qjwkak4CNVGAvvbaawG9hRZayP76179WTXLo0KHBFOinnnrKHnvsMRSgVZMjIQQgAAEIQAACEIAABCAAAQhAoGMI6N2/R48ettxyy9l7772XWw3aMbXhrPVAQDOJtfI3GmQyUSYSZTbx7bffth9++MHWWmstW2SRRQLFuuIqhUcffdTefPNN++CDD6xapWmlMjnemARqogCVhziFAQMGJF7Gvtlmm5kUoPICRoAABCAAAQhA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W2utteyjjz6yr7/+2jzPawCpECFrBGqiAFXjHDRokK288so2dOhQ22uvvWyBBRbImqzUBwIQgAAEIAABCEAAAhCAAAQgAAEIQCBjBDp37hwoQeVA+5tvvgnsgs6cOdNmzZqVsZpSnXolUBMFqDy7v/rqq8Gflrzvvffe9cqDekMAAhCAAAQgAAEIQAACEIAABCAAAQh0AAHZBe3SpUvw16tXL9NqYzlJcn/ffvttB9SKUzYCgZooQGXr04Wdd97ZbfILAQhAAAIQgAAEIAABCEAAAhCAAAQgAIFUBGRqcZlllgn+VMCcOXNyylApRb///vtU5ZKp+QjURAHav3//HDlNUW6WoC8R1113nd100002YcKE4CbcaKONbLPNNrOddtrJBg8eXDMU9957r/39738PzqNzyUZGv379bO2117ajjz7aBg4cWLNzURAEIAABCEAAAhCAAAQgAAEIQAACEMgaAfmcWXbZZYM/1W369Ok2ceJE09J5AgTKEaiJAnTzzTe3Pn362Icffmhjx461SZMm2SqrrFLuvHV/bMqUKbbDDjvY66+/nifLgw8+aPobMWKEXXvttXbggQfmHU+688EHH9iRRx5pDz/8cFHWadOmmcwPXH311Xbsscfa2WefbUsssURROiIgAAEIQAACEIAABCAAAQhAAAIQgECjEejevbstueSS9v7779vs2bMbTTzkqSGBmihANSX50UcfDTy/jxs3zgYMGGB//OMfbZNNNrE11lgjp5mvYb07tChNs9YMT6f81OzLXXfd1VZaaSV74okn7I477jDZpTj44INNM2KlnEwTvvvuO9t9991t/PjxQfbll1/eDjjgANOMW93YL774YjArdO7cuXbxxRfbjBkz7Prrr09zKvJAAAIQgAAEIAABCEAAAhCAAAQgAIG6I7DYYosFepLJkyfbZ599Vnf1p8LtQ6AmClApBM866yxbZ5117J133gmWgp9wwgk5Cbp27Rpo5HMRMRsnnnii6a8ewplnnmmvvfZaUNV9993XRo8eHSxJV4Rmaz711FO2yy67BMrPYcOG2R577GErrrhiYtF+97vf5ZSf22+/faDs7NatW145YqZjuslvuOGGQBGLE6o8ROxAAAIQgAAEIAABCEAAAhCAAAQg0MAE5Dypd+/eJv2TFKE4S2rgi51StAVT5svLpoalZdiyhyllaGHQ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AAQgAAEIQAACDUwABWgDX1xEgwAEIAABCEAAAhCAAAQgAAEIQAACEIBAsxNAAdrsLQD5IQABCEAAAhCAAAQgAAEIQAACEIAABCDQwARQgDbwxUU0CEAAAhCAAAQgAAEIQAACEIAABCAAAQg0OwEUoM3eApAfAhCAAAQgAAEIQAACEIAABCAAAQhAAAINTAAFaANfXESDAAQgAAEIQAACEIAABCAAAQhAAAIQgECzE0AB2uwtAPkhAAEIQAACEIAABCAAAQhAAAIQgAAEINDABFCAtuLizps3z66++mrbZpttrHfv3rbMMsvY9ttvb2eeeaa98MILrSi5OOv7779vhx12mG2wwQbWpUsX69+/f7B/7bXX2qxZs4ozEAMBCEAAAhCAAAQgAAEIQAACEIAABCAAAQjYwjBIR2DKlCm2ww472Ouvv55XwIMPPmj6GzFihEk5eeCBB+YdT7NzwQUX2Kmnnmpz5szJZX/rrbdMf1LAXnnllXbvvfcGCthcAjYgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABFKBp2sBXX31lO+20U075OXDgQNt1111tpZVWsieeeMLuuOMO+/bbb+3ggw+2mTNn2rHHHpvmNEEeKVFPOumkYLtz586233772Y9+9CObPHmy3XPPPfbqq6/aM888Y1tuuaU99NBD1qNHj9TnIiMEIAABCEAAAhCAAAQgAAEIQAACEIAABBqNADNAU1xRLXF/7bXXgpz77ruvjR492jp16hTsH3nkkfbUU0/ZLrvsEig/hw0bZnvssYetuOKKic80ffr0nPK0a9euNnbs2EDR6Qo644wzAiXrzTffbOPHjw+W3l966aXuML8QgAAEIAABCEAAAhCAAAQgAAEIQAACEGh6AtgATdgEvvjiCxs1alSQa5VVVslTfrqiNt98cxszZkywO3fuXLviiivcoUS/I0eODGaSKtO5556bp/xU3CKLLGI33nijDRkyRLt2ww032Ndffx1s8w8CEIAABCAAAQhAAAIQgAAEIAABCEAAAhAwQwGasBXceuut9s033wS5jjrqqNzMz8JiNAO0X79+QbQUoFH7nYVpS+1fc801waGlllrKDjnkkNhkCy64oJ144onBMdVLSlACBCAAAQhAAAIQgAAEIAABCEAAAhCAAAQgEBJAAZqwJcjepgvbbbed24z9lXd4halTp9pjjz0Wm6ZU5EcffWSffvppcFj2Pd0S+7j0W2+9tS2wwALBIc0IJUAAAhCAAAQgAAEIQAACEIAABCAAAQhAAAIhARSgCVvCs88+G+TQzMv11luvbO71118/d7zQW3zuQIkNdx4dlpOlcqF79+62wgorBEmSnqdcuRyDAAQgAAEIQAACEIAABCAAAQhAAAIQgEC9E0ABmvAKvvfee0EOeXyXDc5yQTZCXXj77bfdZlW/EyZMyKXr06dPbrvUhjuXbIB+8sknpZIRDwEIQAACEIAABCAAAQhAAAIQgAAEIACBpiKAF/gEl3vWrFk2b968IEePHj0q5tTMTBfkPClJ+Oqrr3LJ05wrjdf53AljNr788kubP39+zJFkUXIKRYAABCAAAQhAAAIQgAAEIAABCEAAAhCAQHsRQAGagPTMmTNzqRdbbLHcdqmNaJrZs2eXShYb357niq1AJPKYY46xyy67LBLT+s1p06a1vhBKgAAEIAABCEAAAhCAAAQgAAEIQAACEIBABQIoQCsAih7W8nIXOnfu7DZL/i666KK5Y0kVoO15rlwlS2ysvPLKttxyy5nneSVSVB+tma1z5swxt2S/+pykhAAEIAABCEAAAhCAAAQgAAEIQAACEIBAcgIoQBMwi9r8rGYpdzRNNQrTaFXa81zR88Ztn3LKKaa/WoQTTjjBRo4caRtvvHEtiqMMCEAAAhCAAAQgAAEIQAACEIAABCAAAQiUJYATpLJ48g8uueSSuYjvvvsut11qI5qma9eupZLFxrfnuWIrQCQEIAABCEAAAhCAAAQgAAEIQAACEIAABBqAAArQBBdxqaWWyqWOOinKRRZsRNN06dKl4Gj53Wj6aDmlckXTRPOWSk88BCAAAQhAAAIQgAAEIAABCEAAAhCAAASagQAK0ARXWU6NnHf1yZMnV8wZTdOzZ8+K6aMJ+vbtm9uNlpOLLNhwaRZeeGFbdtllC46yCwEIQAACEIAABCAAAQhAAAIQgAAEIACB5iSAAjThde/fv3+QQzMup0+fXjb3hAkTcseHDBmS265mw51Had97772yWeRUaOLEiUGaAQMGWFJ7o2UL5yAEIAABCEAAAhCAAAQgAAEIQAACEIAABOqYAArQhBcv6rzniSeeKJv7ySefzB2P5stFltlYf/31zXmRr3Se559/3r7//vugtKTnKVMFDkEAAhCAAAQgAAEIQAACEIAABCAAAQhAoO4JoABNeAn33nvvXI7Ro0fntgs3Jk2aZI8//ngQPXjwYOvRo0dhkrL7coK0ww47BGlef/11e/nll0umv+GGG3LHdtlll9w2GxCAAAQgAAEIQAACEIAABCAAAQhAAAIQaHYCKEATtoCBAwfaoEGDglx33XWXjRkzpqiEb7/91g477DDT0nSFk08+uSjN3Llz7YUXXsj9ubTRhCrDhSOPPNJmzJjhdnO/9913n1177bXB/rrrrms77bRT7hgbEIAABCAAAQhAAAIQgAAEIAABCEAAAhBodgIoQFO0gEsuucQWWGAB8zzPhg4dan/605/sww8/DBSeTz31VDBz86GHHgpK3mSTTWzPPfcsOsvnn39usgvq/qZNm1aURrM53SzQcePG2U9/+lN79NFH7bvvvrMpU6bYxRdfbLvvvntw3gUXXNDOOeecoF5FBREBAQhAAAIQgAAEIAABCEAAAhCAAAQgAIEmJYACNMWFl1Lz5ptvtiWWWMLmz59vf/jDH2y11VYL9n/84x+bs9m5+uqr2913321STqYNN954o2233XZB9ldeecW23npr69Kli/Xq1cuOP/743CzTCy+80Hbeeee0pyEfBCAAAQhAAAIQgAAEIAABCEAAAhCAAAQakkB6zVxD4qheqH333deeffZZk33PhRZaKMjolrF36tTJhg0bFhxfbrnlqi80JmW3bt3s/vvvt1NPPdW0reDOo215fb/33nuD82mfAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCLQQWLhlk62kBGRzU0vTZ8+ebZqdKcdHmgnar18/69q1a9ni5BRJS+irCZpBqmX2bqm9HCIttthituaaa1qfPn1aNcO0mvOTBgIQgAAEIAABCEAAAhCAAAQgAAEIQAAC9UoABWgNrtziiy9um266afBXg+LKFiGFp/4IEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAQGUCLIGvzIgUEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAQJ0SYAZonV64eq+2lvH/4x//qHcxqD8EIAABCEAAAhCAAAQgAAEIQAACEGgzAh999FGbld1MBaMAbaarnSFZb7jhBtMfAQIQgAAEIAABCEAAAhCAAAQgAAEIQAACbUkABWhb0qXskgQGDhwYOIsqmYADEIAABCAAAQhAAAIQgAAEIAABCECgyQk8++yzNnHixCan0HrxUYC2niElpCAwdOhQGzZsWIqcZIEABCAAAQhAAAIQgAAEIAABCEAAAs1B4KCDDkIBWoNLjROkGkCkCAhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQCCbBFCAZvO6UCsIQAACEIAABCAAAQhAAAIQgAAEIAABCECgBgRQgNYAIkVAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAAC2SSAAjSb14VaQQACEIAABCAAAQhAAAIQgAAEIAABCEAAAjUggAK0BhApAgIQgAAEIAABCEAAAhCAAAQgAAEIQAACEMgmARSg2bwu1AoCEIAABCAAAQhAAAIQgAAEIAABCEAAAhCoAQEUoDWASBEQgAAEIAABCEAAAhCAAAQgAAEIQAACEIBANgmgAM3mdaFWEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAQA0ILFyDMigCApkj8P33ZmPGmL32mtkPP5htuqnZlluarbJKfFX/8x+zd96JP6ZY5dt66+Ljn31mdv31Yd7u3c0228xsq63MlliiOG3WYpIyiqv/uHFmV1xhtsMOZnvtFZfCbPx4s/vvN3vvPbMF/U8ua61ltuuuZn37xqfPUmwaRmpvt9xi9vrrZlOmmPXubTZwoNmee4byx8k3f77Z7bebvfCC2aefhmzUjjbfPC51tuLSyCuukvfll82mTzdbc02zAw8sfX9K4jR5skWK2lQikLSvSNtvJz1PpXq35/HW1v3uu8Nn1iWXmPXsWbrmH3xg9s9/mr39ttlii5mts47ZL35h1qVL6TxZPKJ+49e/Nlt3XbPjjy9dw1rIe/rpZlOnmv3f/5ktskjpc2XtSDWMvv3W7Oaby9f8Zz8z69YtP80rr5ipzb3/vtkaa4TPtB/9yKxTp/x0Wd+rhpFkULokY0/lSfMMVb6shWoZqd733GP2+ONmn39uNnhwOD5XHxMX0ra9uLKyFpdm7Kf3moceCseYSy8d3lN77FF6fJk1mautz9//bjZ7duXUXbua7bNPfro092F+CfW3p2f1v/9t9tJLZp07h+8dBxxgtuii1clSzftcdSVlL1XaPqQZ21H2rl4b1cgjQKAdCQwbNszzm7J30UUXtdlZX3rJ81Zd1fPPk/+3yCKe9/e/x592m23y0xbm3W234nx33OF5Sy5ZnG/99T3v00+L02cpJg2jwvp//bXn9e0byv/73xce9bx58zzvqKM8b8EFixktuqjnjRhRnCdLMWkYPf205622WrG8ak+DBnne228XSzh5sucNGRKf58ADPe+HH4rzZCUmjbyvvup5665bLO/CC3veVVfFS5YmT3xJHRf72mued+65nnf44Z535JGe3wd63nvvFdfnrbc878orq/tz+W+8sbr0t9xSfL4sxKTtK5L222nPU8+MonVXG3TPrAkTokfyt//6V8/r3Ln4Hl1xRc975JH8tFnfO/jgUI4DDihd01rIe9llLbxmzy59riweqYbRU0+1yFc4PnL7emZGw29+E//833xzz/vyy2jK7G9Xw6i9xgxZpVUNI90b225b3JY0Trz00njJ0rS9+JKyFZt07Dd3rueddprnLbRQMT+9d8ycmS35WlubFVYoltP1NdHfNdfMP1Oa+zC/hPrb03hSY+goF23rHe2ddyrLU+l9rnIJ2U6Rpg/Jajs60H8xlB5l9OjR2Yae8dr5twcBAu1HoK0VoJ9/7nmrrBI+BHbf3fMefNDzPvrI8846y/OWWcbzFljA8+KUAMstF+Y59ljP06C98E8KiWgYPz58QdRA5IILPG/KFM/Ty+U++4Tl9OuXXcVVWkZR+bV96KEtD9s4BagUnHoA66X55ps977PPPO/jjz3vkktaXsKz2n+nYTRtmue5AdvQoZ73+uvhgPT55z3PKWrWWcfzvv22heT8+Z73k5+EnAYPDhX0aksjR3per15h/CmntKTP0lYaeb/6yvNWWimU60c/8rz//MfzvvgivE/XWMPz9JHiiSfypUyTJ7+Ejt1LqnQbNSrkUziQjdt3H3Rcu4tLE40rfFHoWDItZ0/bVyTtt9Oep6WmHbfV2rrrJahPn5a2VUoBKgWnnpO6F086yfOkkNc9ud9+Yd5llw2fdx1Horoz676TssC1/1IK0FrIqw9biy/ecq4mULLEAABAAElEQVR6UYBWy0jE9dwWSz2nCsdHbl/PdxfOPz9M36NH+OFn0iTPu+cez9t++zB+vfXCj6QufVZ/q2XUXmOGLHKqlpHqro+6akf6IPzss543fbrn3XST53XtGsZff32xhEnbXnEJ2YtJM/b7y19CRppAoG3dbx984Hkaf6vP3nnn+rinqr0a6r9d31L4e9xxnrfYYiGPo49uKTHNfdiSuz63XF+r8dDVV4fvWm+8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NtE667n2wQdmvXvHO9eKptV2mrZXWEYW95OO/eRk1Y0V5ZhN/AqdsWZRzvauUzPeh199ZfbGG+E779JLtzfx7J8vTR+StXZ00EEH2ZgxY2z06NHmK0OzDz2jNUQBmtEL06jVak8FaKMyRC4IQAACEIAABCAAAQhAAAIQgAAEmoMACtDaXGff3y4BAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQgEBjEkAB2pjXFakgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABnwAKUJoBBCAAAQhAAAIQgAAEIAABCEAAAhCAAAQg0LAEUIA27KVFMAhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAAFKG0AAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQaFgCCzesZAjW1AS+//57GzNmjL322mv2ww8/2KabbmpbbrmlrbLKKrFc/vOf/9g777wTe0yRyrf11lsXHf/ss8/s+uuvD/J2797dNttsM9tqq61siSWWKEqbtYikjFT/NHmict99990Br0suucR69uwZPZTJ7dbKq/y//vWvbd1117Xjjz++rIz33HOPPf744/b555/b4MGDg/a6zjrrlM2ThYO6v2655RZ7/fXXbcqUKda7d28bOHCg7bnnnrbggvHf2ObPn2+33367vfDCC/bpp59a3759g/tm8803LylSa69FyYLb6cD48ePt/vvvt/feey/gstZaa9muu+4ayB5XhTR9khiJ68svv2zTp0+3Nddc0w488MCS/V7ceTsyLikj1TVNHnF64IEH7MUXX7SPPvoouAY/+9nPbP311+9I8ROfu9r+VM/Bhx56KLhHl156adN9tscee5S8P9NyTSxAO2Sotg9O2ibS3J/tIG6qUyRh1B7jqlRCtHGmahmpGkmf5Wmeh20sbk2KT9PvRE9cbf8WzVMv23//+99t9uzZFavbtWtX22efffLSvfLKKyY277//vq2xxhpBf/6jH/3IOnXqlJeuEXa+/fZbu/nmm8uKomd3t27ditK8/fbb9u9//9teeukl69y5czAuPeCAA2zRRRctSttoEePGjbMrrrjCdthhB9trr71ixVOfVs/jxVihykRK3iTPr9a0vTLV4FBWCHgECLQjgWHDhnl+2/cuuuiiNjur/7DzVl111eA8Opf7W2SRRTx/0BF73m222SaXzqWP/u62225F+e644w5vySWXLMrnv0R7vlKnKH2WItIwSpMnKrM/GM7xmjBhQvRQJrdbK6+EOvjgg4P24Q+6SsroD4K9bbfdtqgd+cpD79JLLy2ZLwsHnn76aW+11VYrqrvunUGDBnn+ALSompMnT/aGDBkSm8dX1nm+QrUoTy2uRVGh7RQxb94876ijjvJ0PaN9irb9gbg3YsSI2Jok7ZNeffVVz1e0F51j4YUX9q666qrYc2QlMg2jNHkkr6/wDNpm4bXQ9TnllFOygqRiParpT+fOneuddtpp3kILLVTULvScmjlzZtF50nItKigjEdX0wWnaRNL7MyM4YqtRDaM0fXCzMUrzLE/zPIy9iBmKTNPvFFa/mv6tME897a+wwgpFfXLhM0n7/kfMPLF+85vfxI4l/I9a3pdffpmXthF2nnrqqYqc1DcVBr1jauxTyNT/2O75k10KkzfU/tdff+1JTsn++9//Pla2eh4vxgpUITLN8ytt26tQlVYf1nuSru3o0aNbXVYzF2DNLDyytz+BtlaA+rPnPH+2ZtA57L777t6DDz4YvPCeddZZ3jLLLOMtsMACnj9brUjw5ZZbLshz7LHHehpgFP5deeWVeXn8WUee/0UxeLG84IILPH/mm6cBm/+lNiinX79+sYqcvEI6aCcNozR5ouJpwNGnT5+AjTrurCtAWyuvlAhSPLjBVzkFqHuYSWH47LPPev7MPe+mm27y/C//QX5/hnEUZWa2p02b5rlB/NChQz1/BmigUHn++ec99+Lrz2D1/K+ouTr7M128n/zkJ4Fc/izX4IOE7p2RI0d6vXr1CuILlVCtvRa5k3fQhhScagcrrrii589k8PxZ497HH3/s+bOgcx8E4gYySfqkr776yltppZWC8/gzQTx/dpr3xRdfBP2fP0PE08efJ554ooMIVD5tGkZp8ujFoEePHgGn/fffP8fpuuuu87p06RLEl/pIVlmK9ktRbX/6l7/8JZBJinZtq9198MEHwUuRnoU777yzp74qGtJwjebPyna1fXDaNpHk/swKk8J6VMsobR/cTIzENumzPM3zsPAaZnE/Tb8TlaPa/i2ap962NT4sfM9w+8cdd5y32GKLBX330UcfnRPt/PPPD+L0DDv33HO9SZMmef5sY2/77bcP4tdbb72i/jyXuU43NE7S+EnjRcen8FfPtWhwnNT/XH311cGY64033gie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AAEEEEAAAQQQQAABBBBAAAFfBAiA+rInyQcCCCCAAAIIIIAAAggggAACCCCAAAII5AkQAM0jYQICCCCAAAIIIIAAAggggAACCCCAAAII+CLQ0peMkA8EmiLw888/m0mTJpm//OUvZvbs2eY///mP6d+/v9liiy3MkCFDTN++fZuyeG+++9JLL5l7773X/Otf/zI6rkOvXr3MlltuaU466STTrl07b/JaTkYGDx5spkyZYsaOHWuOPPLIcr5iPv/8c3P55ZeHfnqsffzxx2bttdc2PXv2NFtttZU54ogjTOvWrctaVhYSPfvss2G+OnToYBYsWFDWJi9cuNBcccUV5umnnzZvv/12+L3u3bubdddd1xx44IFm2223LWs5tZrowQcfNLfeemtY1ugxoPtb97/m77DDDjMbb7xxozZdy62tt97afPLJJ2a33XYz1157baOWUwtfuueee8zkyZPD/f/++++blVZayay//vpmv/32M/vuu69Zcskly9rMd955x1xwwQXh+abjnTt3NptvvnlYtu+9995mmWWWKWs5tZgoKaNFixaF178XX3wxLNvnzJlj1lxzzbBsP/zww812221Xi9lv9DY15jxJ65xtdCZS/mJjjHwtt++//36j+7+coW3btuH1vZy0Nk1jrO13a+1vEmWJ7/fjSRjpfvfpflzPr7vvvrtRh/OIESNMv379Cn7Xp3I7LaMXXnjBjBs3zrz11ltm5syZ4e+49dZbz/ziF78wJ554olluueUK2mZhYprXpKuvvtqcc845IcO0adPMOuuskwUSttEKBAwI1LnA3Llzgw022CCQc6Lgv5YtWwY333xznSsFwVVXXRW0atWqoJHadezYMXjggQfqxmn8+PHOQgKgZeV7zJgxQfv27d33Ch1zEggLnnjiibKWV+uJvvjii6BHjx5hfldcccWyNvfJJ58MunbtWtJIgmDBV199VdbyaimRBOCCQYMGlcybBPaCo48+Ovj2228r3vTf//73btn77LNPxd+vhS9I8Db45S9/6fJR6Bz51a9+FXz//fexm3vJJZeULLM222yzQI/RrA1JGn344YeB/Hgs6b3HHnsEP/zwQ9aYim5vJedJ2uds0Y2s8oxKjHRTfS63tSwtVA4VmiYPgSvec5VaV7yCZvpCEmWJ7/fjSRjp7vTtfnz06NFln2O5591tt92Wd4T7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7duHwc7olzRwOXfuXDN8+HBzwgknmO233958/PHH0STheMuWLc1jjz1mtt566zAwqq/Vf/3113npNMh65513mmKv3b/11ltmv/32M9OnT8/7rp2wzz77mGuuucZ06NDBTuIvAggggAACCCCAQMoC1ABNGZjFI4AAAggggAACCFQm0Lt3b3PooYeamTNnhrUw9dtbbrml2WCDDcIFrbbaahUtcMaMGWHAVGuKrrrqqmannXYyK6ywgnn00UfNa6+9Zn766SejgckrrrjCHHzwwea7774Lg6V9+vQJA6ia7rPPPgs7WzriiCOMBjqXWmqpBtvw1FNPhYHRH374IZzet29fs+OOO5q11lrLvPHGG2bKlClhfm6//Xbz/PPPm1dffdW0a9euwTL4gAACCCCAAAIIIJCOADVA03FlqQgggAACCCCAAAJNFJg8eXJYe1MXM27cOKPBx9yhnBqg9jv6uvzo0aPtR6Ovqw8dOtRMnTrVTdOapPr6+7bbbuumaa3Rfv36ha/Q68R7773X7LLLLm6+BlD1FXcNtOpwySWXGH1Nf4kl/tfc/sKFC81hhx1mJk2aFKY59thj817ZD2fwHwIIIIAAAggggEDiAv+7K0t80SwQAQQQQAABBBBAAIHaENh5553Nueee22BjWrVqZU488cQG00aNGtUg+Kkztdaotv9ph3feeceOhn+vvfZaF/w84IADwtfto8FPTbT00kubiRMnhu2R6mcN6L799ts6yoAAAggggAACCCCQsgAB0JSBWTwCCCCAAAIIIIBA9QVOPvlk06JFi7wNWXfddd00nT9ixAj3OTpiO13SadrmaHQYP368+xitYeomRka0ZqgOWvv01ltvjcxhFAEEEEAAAQQQQCAtAQKgacmyXAQQQAABBBBAAIGaEdh0000LbsuKK67opmvbom3atHGfoyPLLLOM+/jjjz+68Z9//tm8++674Wd9HV9ri5Yaor3Nz549u1RS5iGAAAIIIIAAAggkJEAnSAlBshgEEEAAAQQQQACB2hTQDou0B/dCQ7RWaKnOlXJfabfL+uCDD4wNiH700Uemc+fOdlbsXwKgsUQkQAABBBBAAAEEEhEgAJoIIwtBAAEEEEAAAQQQqFWBaO3NUtu45JJLlppdcN6cOXPc9EWLFrmOktzEEiO5bYmWSMosBBBAAAEEEEAAgSYIEABtAh5f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)

``` r
create_summary_table(
  data = sim_data, 
  treat_var = "treat", 
  table_title = "Characteristics by Treatment Arm",
  vars_continuous = c("z1", "z2", "size", "z3", "z4", "z5"),
  vars_categorical = c("flag.harm", "grade3"),
  font_size = 12
)
```

| Characteristics by Treatment Arm |  |  |  |  |  |
|----|----|----|----|----|----|
| Characteristic |  | Control (n=350) | Treatment (n=350) | P-value¹ | SMD² |
| z1 | Mean (SD) | 0.3 (0.4) | 0.2 (0.4) | 0.729 | 0.03 |
| z2 | Mean (SD) | 2.4 (1.1) | 2.5 (1.1) | 0.157 | 0.11 |
| size | Mean (SD) | 29.2 (14.3) | 30.1 (17.0) | 0.461 | 0.06 |
| z3 | Mean (SD) | 0.4 (0.5) | 0.4 (0.5) | 0.249 | 0.09 |
| z4 | Mean (SD) | 2.6 (1.1) | 2.5 (1.1) | 0.428 | 0.06 |
| z5 | Mean (SD) | 2.4 (1.1) | 2.4 (1.1) | 0.563 | 0.04 |
| flag.harm |  |  |  | 0.908 | 0.00 |
|  | 0 | 306 (87.4%) | 308 (88.0%) |  |  |
|  | 1 | 44 (12.6%) | 42 (12.0%) |  |  |
| grade3 |  |  |  | 0.530 | 0.02 |
|  | 0 | 273 (78.0%) | 265 (75.7%) |  |  |
|  | 1 | 77 (22.0%) | 85 (24.3%) |  |  |
| ¹ P-values: t-test for continuous, chi-square/Fisher's exact for categorical/binary variables |  |  |  |  |  |
| ² SMD = Standardized mean difference (Cohen's d for continuous, Cramer's V for categorical) |  |  |  |  |  |

## 5 Running Simulation Studies

### 5.1 Setting Up Parallel Processing

For efficient simulation studies, use parallel processing:

``` r
# Configure parallel backend
n_workers <- min(parallel::detectCores() - 1, 120)

plan(multisession, workers = n_workers)
registerDoFuture()

cat("Using", n_workers, "parallel workers\n")
```

    Using 13 parallel workers

### 5.2 Define Simulation Parameters

``` r
# Simulation settings
sim_config_alt <- list(
  n_sims = 100,          # Number of simulations (use 500-1000 for final)
  n_sample = 700,         # Sample size per trial
  max_follow = 84,        # Maximum follow-up (months)
  seed_base = 8316951,
  muC_adj = log(1.5)
)

sim_config_null <- list(
  n_sims = 100,          # More simulations for Type I error estimation
  n_sample = 700,         # Sample size per trial
  max_follow = 84,        # Maximum follow-up (months)
  seed_base = 8316951,
  muC_adj = log(1.5)
)

# ForestSearch parameters (now includes use_twostage)
fs_params <- list(
  outcome.name = "y.sim",
  event.name = "event.sim",
  treat.name = "treat",
  id.name = "id",
  use_lasso = TRUE,
  use_grf = TRUE,
  hr.threshold = 1.25,
  hr.consistency = 1.0,
  pconsistency.threshold = 0.90,
  fs.splits = 400,
  n.min = 60,
  d0.min = 12,
  d1.min = 12,
  maxk = 2,
  by.risk = 12,
  vi.grf.min = -0.2,
  # NEW: Two-stage algorithm option
  use_twostage = TRUE,      # Set TRUE for faster exploratory analysis
  twostage_args = list()     # Optional tuning parameters
)


# Confounders for analysis
confounders_base <- c("z1", "z2", "z3", "z4", "z5", "size", "grade3")
```

#### 5.2.1 Two-Stage Algorithm Option

The `use_twostage` parameter enables a faster two-stage search
algorithm:

``` r
# Fast configuration with two-stage algorithm
fs_params_fast <- modifyList(fs_params, list(
  use_twostage = TRUE,
  twostage_args = list(
    n.splits.screen = 30,    # Stage 1 screening splits
    batch.size = 20,         # Stage 2 batch size
    conf.level = 0.95        # Early stopping confidence
  )
))

cat("Standard search: use_twostage =", fs_params$use_twostage, "\n")
```

    Standard search: use_twostage = TRUE 

``` r
cat("Fast search: use_twostage =", fs_params_fast$use_twostage, "\n")
```

    Fast search: use_twostage = TRUE 

### 5.3 Running Alternative Hypothesis Simulations

``` r
cat("Running", sim_config_alt$n_sims, "simulations under H1...\n")
```

    Running 100 simulations under H1...

``` r
start_time <- Sys.time()
t0 <- proc.time()[3]

results_alt <- foreach(
  sim = 1:sim_config_alt$n_sims,
  .combine = rbind,
  .errorhandling = "remove",
  .options.future = list(
    packages = c("forestsearch", "survival", "data.table"),
    seed = TRUE
  )
) %dofuture% {
  run_simulation_analysis(
    sim_id = sim,
    dgm = dgm_calibrated,
    n_sample = sim_config_alt$n_sample,
    max_follow = sim_config_alt$max_follow,
    muC_adj = sim_config_alt$muC_adj,
    confounders_base = confounders_base,
    cox_formula_adj = survival::Surv(y.sim, event.sim) ~ treat + z1 + z2 + z3,
    n_add_noise = 0L,
    run_fs = TRUE,
    run_fs_grf = FALSE,
    run_grf = TRUE,
    fs_params = fs_params,
    verbose = TRUE,
    debug = FALSE,
    verbose_n = 3  # Only print first 3 simulations
  )
}
```

    === Two-Stage Consistency Evaluation Enabled ===
    Stage 1 screening splits: 30
    Maximum total splits: 400
    Batch size: 20
    ================================================

    GRF stage for cut selection with dmin, tau = 12 0.6
    tau, maxdepth = 48.53742 2
       leaf.node control.mean control.size control.se depth
    1          2        -5.33       520.00       1.11     1
    2          3         1.02       180.00       2.32     1
    11         4        -6.71       426.00       1.21     2
    21         5         3.44       142.00       2.62     2
    4          7        -7.11        78.00       3.18     2
    GRF subgroup NOT found
    GRF cuts identified: 0
    # of continuous/categorical characteristics 4 3
    Continuous characteristics: z2 z4 z5 size
    Categorical characteristics: z1 z3 grade3
    ## Prior to lasso: z2 z4 z5 size
    #### Lasso selection results
    7 x 1 sparse Matrix of class "dgCMatrix"
                    s0
    z1      0.28638272
    z2      .
    z3      0.04573174
    z4     -0.34834939
    z5      0.45235142
    size    .
    grade3  0.01493646
    Cox-LASSO selected: z1 z3 z4 z5 grade3
    Cox-LASSO not selected: z2 size
    ### End Lasso selection
    ## After lasso: z4 z5
    Default cuts included from Lasso: z4 <= mean(z4) z4 <= median(z4) z4 <= qlow(z4) z4 <= qhigh(z4) z5 <= mean(z5) z5 <= median(z5) z5 <= qlow(z5) z5 <= qhigh(z5)
    Categorical after Lasso: z1 z3 grade3
    Factors per GRF:
    Initial GRF cuts included
    Factors included per GRF (not in lasso)

    ===== CONSOLIDATED CUT EVALUATION (IMPROVED) =====
    Evaluating 10 cut expressions once and caching...
    Cut evaluation summary:
      Total cuts:  10
      Valid cuts:  10
      Errors:  0
    ✓ All 10 factors validated as 0/1
    ===== END CONSOLIDATED CUT EVALUATION =====

    # of candidate subgroup factors= 10
     [1] "z4 <= 2.5" "z4 <= 3"   "z4 <= 2"   "z5 <= 2.4" "z5 <= 2"   "z5 <= 1"
     [7] "z5 <= 3"   "z1"        "z3"        "grade3"
    Number of possible configurations (<= maxk): maxk = 2 , # combinations = 210
    Events criteria: control >= 12 , treatment >= 12
    Sample size criteria: n >= 60
    Subgroup search completed in 0.02 minutes

    --- Filtering Summary ---
      Combinations evaluated: 210
      Passed variance check: 189
      Passed prevalence (>= 0.025 ): 189
      Passed redundancy check: 178
      Passed event counts (d0>= 12 , d1>= 12 ): 160
      Passed sample size (n>= 60 ): 158
      Cox model fit successfully: 158
      Passed HR threshold (>= 1.25 ): 3
    -------------------------

    Found 3 subgroup candidate(s)
    # of candidate subgroups (meeting all criteria) = 3
    Random seed set to: 8316951
    Two-stage parameters:
      n.splits.screen: 30
      screen.threshold: 0.763
      batch.size: 20
      conf.level: 0.95
    Removed 1 near-duplicate subgroups
    # of unique initial candidates: 2
    # Restricting to top stop_Kgroups = 10
    # of candidates to evaluate: 2
    # Early stop threshold: 0.95
    Parallel config: workers = 6 , batch_size = 1
    Batch 1 / 2 : candidates 1 - 1
    Batch 2 / 2 : candidates 2 - 2
    Evaluated 2 of 2 candidates (complete)
    1 subgroups passed consistency threshold
    SG focus = hr
    Seconds and minutes forestsearch overall = 10.04 0.1673
    Consistency algorithm used: twostage
    tau, maxdepth = 48.53742 2
       leaf.node control.mean control.size control.se depth
    1          2        -5.33       520.00       1.11     1
    2          3         1.02       180.00       2.32     1
    11         4        -6.71       426.00       1.21     2
    21         5         3.44       142.00       2.62     2
    4          7        -7.11        78.00       3.18     2
    GRF subgroup NOT found

    === Two-Stage Consistency Evaluation Enabled ===
    Stage 1 screening splits: 30
    Maximum total splits: 400
    Batch size: 20
    ================================================

    GRF stage for cut selection with dmin, tau = 12 0.6
    tau, maxdepth = 49.03326 2
      leaf.node control.mean control.size control.se depth
    1         2        -5.56       534.00       1.18     1
    2         3         2.39       166.00       2.32     1
    3         4        -5.10       262.00       1.63     2
    4         5        10.27        96.00       2.89     2
    5         6        -7.05       336.00       1.52     2

    Selected subgroup:
      leaf.node control.mean control.size control.se depth
    4         5        10.27        96.00       2.89     2

    GRF subgroup found
    Terminating node at max.diff (sg.harm.id):
    [1] "z1 <= 0"

    All splits (from all trees):
    [1] "z1 <= 0"    "z2 <= 2"    "size <= 58"
    GRF cuts identified: 3
      Cuts: z1 <= 0, z2 <= 2, size <= 58
      Selected tree depth: 2
    # of continuous/categorical characteristics 4 3
    Continuous characteristics: z2 z4 z5 size
    Categorical characteristics: z1 z3 grade3
    ## Prior to lasso: z2 z4 z5 size
    #### Lasso selection results
    7 x 1 sparse Matrix of class "dgCMatrix"
                   s0
    z1      .
    z2      .
    z3      .
    z4     -0.2924157
    z5      0.4106361
    size    .
    grade3  .
    Cox-LASSO selected: z4 z5
    Cox-LASSO not selected: z1 z2 z3 size grade3
    ### End Lasso selection
    ## After lasso: z4 z5
    Default cuts included from Lasso: z4 <= mean(z4) z4 <= median(z4) z4 <= qlow(z4) z4 <= qhigh(z4) z5 <= mean(z5) z5 <= median(z5) z5 <= qlow(z5) z5 <= qhigh(z5)
    Categorical after Lasso:
    Factors per GRF: z1 <= 0 z2 <= 2 size <= 58
    Initial GRF cuts included z1 <= 0 z2 <= 2 size <= 58
    Factors included per GRF (not in lasso) z1 <= 0 z2 <= 2 size <= 58

    ===== CONSOLIDATED CUT EVALUATION (IMPROVED) =====
    Evaluating 10 cut expressions once and caching...
    Cut evaluation summary:
      Total cuts:  10
      Valid cuts:  10
      Errors:  0
    ✓ All 10 factors validated as 0/1
    ===== END CONSOLIDATED CUT EVALUATION =====

    # of candidate subgroup factors= 10
     [1] "z2 <= 2"    "size <= 58" "z4 <= 2.5"  "z4 <= 3"    "z4 <= 2"
     [6] "z5 <= 2.5"  "z5 <= 2"    "z5 <= 1.8"  "z5 <= 3"    "z1 <= 0"
    Number of possible configurations (<= maxk): maxk = 2 , # combinations = 210
    Events criteria: control >= 12 , treatment >= 12
    Sample size criteria: n >= 60
    Subgroup search completed in 0.01 minutes

    --- Filtering Summary ---
      Combinations evaluated: 210
      Passed variance check: 189
      Passed prevalence (>= 0.025 ): 189
      Passed redundancy check: 174
      Passed event counts (d0>= 12 , d1>= 12 ): 143
      Passed sample size (n>= 60 ): 142
      Cox model fit successfully: 142
      Passed HR threshold (>= 1.25 ): 9
    -------------------------

    Found 9 subgroup candidate(s)
    # of candidate subgroups (meeting all criteria) = 9
    Random seed set to: 8316951
    Two-stage parameters:
      n.splits.screen: 30
      screen.threshold: 0.763
      batch.size: 20
      conf.level: 0.95
    Removed 3 near-duplicate subgroups
    # of unique initial candidates: 6
    # Restricting to top stop_Kgroups = 10
    # of candidates to evaluate: 6
    # Early stop threshold: 0.95
    Parallel config: workers = 6 , batch_size = 1
    Batch 1 / 6 : candidates 1 - 1

    ==================================================
    EARLY STOP TRIGGERED (batch 1 )
      Candidate: 1 of 6
      Pcons: 1 >= 0.95
    ==================================================

    Evaluated 1 of 6 candidates (early stop)
    1 subgroups passed consistency threshold
    SG focus = hr
    Seconds and minutes forestsearch overall = 4.088 0.0681
    Consistency algorithm used: twostage
    tau, maxdepth = 49.03326 2
      leaf.node control.mean control.size control.se depth
    1         2        -5.56       534.00       1.18     1
    2         3         2.39       166.00       2.32     1
    3         4        -5.10       262.00       1.63     2
    4         5        10.27        96.00       2.89     2
    5         6        -7.05       336.00       1.52     2

    Selected subgroup:
      leaf.node control.mean control.size control.se depth
    4         5        10.27        96.00       2.89     2

    GRF subgroup found
    Terminating node at max.diff (sg.harm.id):
    [1] "z1 <= 0"

    All splits (from all trees):
    [1] "z1 <= 0"    "z2 <= 2"    "size <= 58"

    === Two-Stage Consistency Evaluation Enabled ===
    Stage 1 screening splits: 30
    Maximum total splits: 400
    Batch size: 20
    ================================================

    GRF stage for cut selection with dmin, tau = 12 0.6
    tau, maxdepth = 49.57713 2
       leaf.node control.mean control.size control.se depth
    1          2        -6.22       494.00       1.21     1
    2          3         5.44       206.00       2.11     1
    21         5        -6.09       123.00       2.07     2
    3          6        -6.92       385.00       1.46     2
    4          7         7.19       143.00       2.55     2

    Selected subgroup:
      leaf.node control.mean control.size control.se depth
    4         7         7.19       143.00       2.55     2

    GRF subgroup found
    Terminating node at max.diff (sg.harm.id):
    [1] "z1 <= 0"

    All splits (from all trees):
    [1] "z1 <= 0"    "z5 <= 1"    "size <= 18"
    GRF cuts identified: 3
      Cuts: z1 <= 0, z5 <= 1, size <= 18
      Selected tree depth: 2
    # of continuous/categorical characteristics 4 3
    Continuous characteristics: z2 z4 z5 size
    Categorical characteristics: z1 z3 grade3
    ## Prior to lasso: z2 z4 z5 size
    #### Lasso selection results
    7 x 1 sparse Matrix of class "dgCMatrix"
                    s0
    z1      .
    z2      .
    z3      0.07273592
    z4     -0.40237672
    z5      0.35008925
    size    .
    grade3  .
    Cox-LASSO selected: z3 z4 z5
    Cox-LASSO not selected: z1 z2 size grade3
    ### End Lasso selection
    ## After lasso: z4 z5
    Default cuts included from Lasso: z4 <= mean(z4) z4 <= median(z4) z4 <= qlow(z4) z4 <= qhigh(z4) z5 <= mean(z5) z5 <= median(z5) z5 <= qlow(z5) z5 <= qhigh(z5)
    Categorical after Lasso: z3
    Factors per GRF: z1 <= 0 z5 <= 1 size <= 18
    Initial GRF cuts included z1 <= 0 z5 <= 1 size <= 18
    Factors included per GRF (not in lasso) z1 <= 0 z5 <= 1 size <= 18

    ===== CONSOLIDATED CUT EVALUATION (IMPROVED) =====
    Evaluating 11 cut expressions once and caching...
    Cut evaluation summary:
      Total cuts:  11
      Valid cuts:  10
      Errors:  0
    Dropping variables (cut only has 1 level): z4 <= 4
    Total cuts after dropping invalid: 10
    ✓ All 10 factors validated as 0/1
    ===== END CONSOLIDATED CUT EVALUATION =====

    # of candidate subgroup factors= 10
     [1] "z5 <= 1"    "size <= 18" "z4 <= 2.5"  "z4 <= 2"    "z4 <= 1"
     [6] "z5 <= 2.5"  "z5 <= 2"    "z5 <= 3"    "z3"         "z1 <= 0"
    Number of possible configurations (<= maxk): maxk = 2 , # combinations = 210
    Events criteria: control >= 12 , treatment >= 12
    Sample size criteria: n >= 60
    Subgroup search completed in 0.01 minutes

    --- Filtering Summary ---
      Combinations evaluated: 210
      Passed variance check: 189
      Passed prevalence (>= 0.025 ): 189
      Passed redundancy check: 178
      Passed event counts (d0>= 12 , d1>= 12 ): 162
      Passed sample size (n>= 60 ): 154
      Cox model fit successfully: 154
      Passed HR threshold (>= 1.25 ): 14
    -------------------------

    Found 14 subgroup candidate(s)
    # of candidate subgroups (meeting all criteria) = 14
    Random seed set to: 8316951
    Two-stage parameters:
      n.splits.screen: 30
      screen.threshold: 0.763
      batch.size: 20
      conf.level: 0.95
    Removed 4 near-duplicate subgroups
    # of unique initial candidates: 10
    # Restricting to top stop_Kgroups = 10
    # of candidates to evaluate: 10
    # Early stop threshold: 0.95
    Parallel config: workers = 6 , batch_size = 1
    Batch 1 / 10 : candidates 1 - 1

    ==================================================
    EARLY STOP TRIGGERED (batch 1 )
      Candidate: 1 of 10
      Pcons: 1 >= 0.95
    ==================================================

    Evaluated 1 of 10 candidates (early stop)
    1 subgroups passed consistency threshold
    SG focus = hr
    Seconds and minutes forestsearch overall = 4.251 0.0709
    Consistency algorithm used: twostage
    tau, maxdepth = 49.57713 2
       leaf.node control.mean control.size control.se depth
    1          2        -6.22       494.00       1.21     1
    2          3         5.44       206.00       2.11     1
    21         5        -6.09       123.00       2.07     2
    3          6        -6.92       385.00       1.46     2
    4          7         7.19       143.00       2.55     2

    Selected subgroup:
      leaf.node control.mean control.size control.se depth
    4         7         7.19       143.00       2.55     2

    GRF subgroup found
    Terminating node at max.diff (sg.harm.id):
    [1] "z1 <= 0"

    All splits (from all trees):
    [1] "z1 <= 0"    "z5 <= 1"    "size <= 18"

``` r
timings$sims_alt_elapsed <- proc.time()[3] - t0
runtime_alt <- difftime(Sys.time(), start_time, units = "mins")
timings$sims_alt_wall <- as.numeric(runtime_alt) * 60  # store in seconds
cat("Completed in", round(runtime_alt, 1), "minutes\n")
```

    Completed in 1.1 minutes

``` r
cat("Results:", nrow(results_alt), "rows\n")
```

    Results: 200 rows

### 5.4 Running Null Hypothesis Simulations

``` r
cat("Running", sim_config_null$n_sims, "simulations under H0...\n")
```

    Running 100 simulations under H0...

``` r
start_time <- Sys.time()
t0 <- proc.time()[3]

results_null <- foreach(
  sim = 1:sim_config_null$n_sims,
  .combine = rbind,
  .errorhandling = "remove",
  .options.future = list(
    packages = c("forestsearch", "survival", "data.table"),
    seed = TRUE
  )
) %dofuture% {
  
  run_simulation_analysis(
    sim_id = sim,
    dgm = dgm_null,
    n_sample = sim_config_null$n_sample,
    max_follow = sim_config_null$max_follow,
    muC_adj = sim_config_null$muC_adj,
    confounders_base = confounders_base,
    cox_formula_adj = survival::Surv(y.sim, event.sim) ~ treat + z1 + z2 + z3,
    n_add_noise = 0L,
    run_fs = TRUE,
    run_fs_grf = FALSE,
    run_grf = TRUE,
    fs_params = fs_params,
    verbose = TRUE,
    verbose_n = 3  # Only print first 3 simulations

  )
}
```

    === Two-Stage Consistency Evaluation Enabled ===
    Stage 1 screening splits: 30
    Maximum total splits: 400
    Batch size: 20
    ================================================

    GRF stage for cut selection with dmin, tau = 12 0.6
    tau, maxdepth = 47.91247 2
      leaf.node control.mean control.size control.se depth
    2         3        -4.78       695.00       1.00     1
    1         4        -5.09       568.00       1.10     2
    GRF subgroup NOT found
    GRF cuts identified: 0
    # of continuous/categorical characteristics 4 3
    Continuous characteristics: z2 z4 z5 size
    Categorical characteristics: z1 z3 grade3
    ## Prior to lasso: z2 z4 z5 size
    #### Lasso selection results
    7 x 1 sparse Matrix of class "dgCMatrix"
                      s0
    z1      0.0078483334
    z2      0.0324623135
    z3      .
    z4     -0.3338944143
    z5      0.4663620842
    size    0.0002302206
    grade3  .
    Cox-LASSO selected: z1 z2 z4 z5 size
    Cox-LASSO not selected: z3 grade3
    ### End Lasso selection
    ## After lasso: z2 z4 z5 size
    Default cuts included from Lasso: z2 <= mean(z2) z2 <= median(z2) z2 <= qlow(z2) z2 <= qhigh(z2) z4 <= mean(z4) z4 <= median(z4) z4 <= qlow(z4) z4 <= qhigh(z4) z5 <= mean(z5) z5 <= median(z5) z5 <= qlow(z5) z5 <= qhigh(z5) size <= mean(size) size <= median(size) size <= qlow(size) size <= qhigh(size)
    Categorical after Lasso: z1
    Factors per GRF:
    Initial GRF cuts included
    Factors included per GRF (not in lasso)

    ===== CONSOLIDATED CUT EVALUATION (IMPROVED) =====
    Evaluating 15 cut expressions once and caching...
    Cut evaluation summary:
      Total cuts:  15
      Valid cuts:  14
      Errors:  0
    Dropping variables (cut only has 1 level): z2 <= 4
    Total cuts after dropping invalid: 14
    ✓ All 14 factors validated as 0/1
    ===== END CONSOLIDATED CUT EVALUATION =====

    # of candidate subgroup factors= 14
     [1] "z2 <= 2.5"    "z2 <= 2"      "z4 <= 2.5"    "z4 <= 3"      "z4 <= 2"
     [6] "z5 <= 2.4"    "z5 <= 2"      "z5 <= 1"      "z5 <= 3"      "size <= 29.1"
    [11] "size <= 25"   "size <= 20"   "size <= 35"   "z1"
    Number of possible configurations (<= maxk): maxk = 2 , # combinations = 406
    Events criteria: control >= 12 , treatment >= 12
    Sample size criteria: n >= 60
    Subgroup search completed in 0.03 minutes

    --- Filtering Summary ---
      Combinations evaluated: 406
      Passed variance check: 373
      Passed prevalence (>= 0.025 ): 373
      Passed redundancy check: 354
      Passed event counts (d0>= 12 , d1>= 12 ): 327
      Passed sample size (n>= 60 ): 316
      Cox model fit successfully: 316
      Passed HR threshold (>= 1.25 ): 4
    -------------------------

    Found 4 subgroup candidate(s)
    # of candidate subgroups (meeting all criteria) = 4
    Random seed set to: 8316951
    Two-stage parameters:
      n.splits.screen: 30
      screen.threshold: 0.763
      batch.size: 20
      conf.level: 0.95
    Removed 2 near-duplicate subgroups
    # of unique initial candidates: 2
    # Restricting to top stop_Kgroups = 10
    # of candidates to evaluate: 2
    # Early stop threshold: 0.95
    Parallel config: workers = 6 , batch_size = 1
    Batch 1 / 2 : candidates 1 - 1
    Batch 2 / 2 : candidates 2 - 2
    Evaluated 2 of 2 candidates (complete)
    No subgroups found meeting consistency threshold
    Seconds and minutes forestsearch overall = 11.296 0.1883
    Consistency algorithm used: twostage
    tau, maxdepth = 47.91247 2
      leaf.node control.mean control.size control.se depth
    2         3        -4.78       695.00       1.00     1
    1         4        -5.09       568.00       1.10     2
    GRF subgroup NOT found

    === Two-Stage Consistency Evaluation Enabled ===
    Stage 1 screening splits: 30
    Maximum total splits: 400
    Batch size: 20
    ================================================

    GRF stage for cut selection with dmin, tau = 12 0.6
    tau, maxdepth = 49.50766 2
      leaf.node control.mean control.size control.se depth
    1         2        -4.62       689.00       1.07     1
    2         5        -6.72       124.00       2.32     2
    3         6        -5.20       513.00       1.27     2
    GRF subgroup NOT found
    GRF cuts identified: 0
    # of continuous/categorical characteristics 4 3
    Continuous characteristics: z2 z4 z5 size
    Categorical characteristics: z1 z3 grade3
    ## Prior to lasso: z2 z4 z5 size
    #### Lasso selection results
    7 x 1 sparse Matrix of class "dgCMatrix"
                   s0
    z1      .
    z2      .
    z3     -0.1081342
    z4     -0.2654241
    z5      0.4210064
    size    .
    grade3  .
    Cox-LASSO selected: z3 z4 z5
    Cox-LASSO not selected: z1 z2 size grade3
    ### End Lasso selection
    ## After lasso: z4 z5
    Default cuts included from Lasso: z4 <= mean(z4) z4 <= median(z4) z4 <= qlow(z4) z4 <= qhigh(z4) z5 <= mean(z5) z5 <= median(z5) z5 <= qlow(z5) z5 <= qhigh(z5)
    Categorical after Lasso: z3
    Factors per GRF:
    Initial GRF cuts included
    Factors included per GRF (not in lasso)

    ===== CONSOLIDATED CUT EVALUATION (IMPROVED) =====
    Evaluating 8 cut expressions once and caching...
    Cut evaluation summary:
      Total cuts:  8
      Valid cuts:  8
      Errors:  0
    ✓ All 8 factors validated as 0/1
    ===== END CONSOLIDATED CUT EVALUATION =====

    # of candidate subgroup factors= 8
    [1] "z4 <= 2.5" "z4 <= 3"   "z4 <= 2"   "z5 <= 2.5" "z5 <= 2"   "z5 <= 1.8"
    [7] "z5 <= 3"   "z3"
    Number of possible configurations (<= maxk): maxk = 2 , # combinations = 136
    Events criteria: control >= 12 , treatment >= 12
    Sample size criteria: n >= 60
    Subgroup search completed in 0.01 minutes

    --- Filtering Summary ---
      Combinations evaluated: 136
      Passed variance check: 117
      Passed prevalence (>= 0.025 ): 117
      Passed redundancy check: 106
      Passed event counts (d0>= 12 , d1>= 12 ): 98
      Passed sample size (n>= 60 ): 98
      Cox model fit successfully: 98
      Passed HR threshold (>= 1.25 ): 0
    -------------------------

    NO subgroup candidates found
    tau, maxdepth = 49.50766 2
      leaf.node control.mean control.size control.se depth
    1         2        -4.62       689.00       1.07     1
    2         5        -6.72       124.00       2.32     2
    3         6        -5.20       513.00       1.27     2
    GRF subgroup NOT found

    === Two-Stage Consistency Evaluation Enabled ===
    Stage 1 screening splits: 30
    Maximum total splits: 400
    Batch size: 20
    ================================================

    GRF stage for cut selection with dmin, tau = 12 0.6
    tau, maxdepth = 46.7094 2
      leaf.node control.mean control.size control.se depth
    1         2         2.55       101.00       2.51     1
    2         3        -5.13       599.00       1.05     1
    3         4         3.54        91.00       2.44     2
    4         5        -6.03       449.00       1.13     2
    5         6        -7.26       105.00       2.90     2
    GRF subgroup NOT found
    GRF cuts identified: 0
    # of continuous/categorical characteristics 4 3
    Continuous characteristics: z2 z4 z5 size
    Categorical characteristics: z1 z3 grade3
    ## Prior to lasso: z2 z4 z5 size
    #### Lasso selection results
    7 x 1 sparse Matrix of class "dgCMatrix"
                      s0
    z1     -0.2046970204
    z2      .
    z3     -0.0162003130
    z4     -0.3904744742
    z5      0.3456587351
    size    .
    grade3  0.0002404273
    Cox-LASSO selected: z1 z3 z4 z5 grade3
    Cox-LASSO not selected: z2 size
    ### End Lasso selection
    ## After lasso: z4 z5
    Default cuts included from Lasso: z4 <= mean(z4) z4 <= median(z4) z4 <= qlow(z4) z4 <= qhigh(z4) z5 <= mean(z5) z5 <= median(z5) z5 <= qlow(z5) z5 <= qhigh(z5)
    Categorical after Lasso: z1 z3 grade3
    Factors per GRF:
    Initial GRF cuts included
    Factors included per GRF (not in lasso)

    ===== CONSOLIDATED CUT EVALUATION (IMPROVED) =====
    Evaluating 10 cut expressions once and caching...
    Cut evaluation summary:
      Total cuts:  10
      Valid cuts:  9
      Errors:  0
    Dropping variables (cut only has 1 level): z4 <= 4
    Total cuts after dropping invalid: 9
    ✓ All 9 factors validated as 0/1
    ===== END CONSOLIDATED CUT EVALUATION =====

    # of candidate subgroup factors= 9
    [1] "z4 <= 2.5" "z4 <= 2"   "z4 <= 1"   "z5 <= 2.5" "z5 <= 2"   "z5 <= 3"
    [7] "z1"        "z3"        "grade3"
    Number of possible configurations (<= maxk): maxk = 2 , # combinations = 171
    Events criteria: control >= 12 , treatment >= 12
    Sample size criteria: n >= 60
    Subgroup search completed in 0.01 minutes

    --- Filtering Summary ---
      Combinations evaluated: 171
      Passed variance check: 154
      Passed prevalence (>= 0.025 ): 154
      Passed redundancy check: 146
      Passed event counts (d0>= 12 , d1>= 12 ): 141
      Passed sample size (n>= 60 ): 135
      Cox model fit successfully: 135
      Passed HR threshold (>= 1.25 ): 0
    -------------------------

    NO subgroup candidates found
    tau, maxdepth = 46.7094 2
      leaf.node control.mean control.size control.se depth
    1         2         2.55       101.00       2.51     1
    2         3        -5.13       599.00       1.05     1
    3         4         3.54        91.00       2.44     2
    4         5        -6.03       449.00       1.13     2
    5         6        -7.26       105.00       2.90     2
    GRF subgroup NOT found

``` r
timings$sims_null_elapsed <- proc.time()[3] - t0
runtime_null <- difftime(Sys.time(), start_time, units = "mins")
timings$sims_null_wall <- as.numeric(runtime_null) * 60
cat("Completed in", round(runtime_null, 1), "minutes\n")
```

    Completed in 0.8 minutes

## 6 Summarizing Results

### 6.1 Operating Characteristics Summary

``` r
t0 <- proc.time()[3]

# Summarize alternative hypothesis results
summary_alt <- summarize_simulation_results(results_alt)
print(summary_alt)
```

                      FS     GRF
    any.H          0.890   0.770
    sens           0.850   0.920
    spec           0.980   0.970
    ppv            0.880   0.850
    npv            0.980   0.990
    Avg(#H)       86.000  98.000
    minH          61.000  62.000
    maxH         148.000 224.000
    Avg(#Hc)     624.000 625.000
    minHc        552.000 476.000
    maxHc        700.000 700.000
    hat(H*)        2.378     NaN
    hat(hat[H])    2.383   2.266
    hat(Hc*)       0.631     NaN
    hat(hat[Hc])   0.638   0.628
    hat(H*)all     2.378     NaN
    hat(Hc*)all    0.631     NaN
    hat(ITT)       0.746   0.746
    hat(ITTadj)    0.761   0.761

``` r
timings$summarize_alt <- proc.time()[3] - t0
```

``` r
t0 <- proc.time()[3]

# Summarize null hypothesis results
summary_null <- summarize_simulation_results(results_null)
print(summary_null)
```

                      FS     GRF
    any.H          0.070   0.100
    sens             NaN     NaN
    spec           0.850   0.890
    ppv            0.000   0.000
    npv            1.000   1.000
    Avg(#H)      106.000  75.000
    minH          81.000  61.000
    maxH         154.000  99.000
    Avg(#Hc)     693.000 692.000
    minHc        546.000 601.000
    maxHc        700.000 700.000
    hat(H*)          NaN     NaN
    hat(hat[H])    1.739   1.385
    hat(Hc*)       0.782     NaN
    hat(hat[Hc])   0.687   0.713
    hat(H*)all       NaN     NaN
    hat(Hc*)all    0.782     NaN
    hat(ITT)       0.697   0.697
    hat(ITTadj)    0.688   0.688

``` r
timings$summarize_null <- proc.time()[3] - t0
```

## 7 Theoretical Subgroup Detection Rate Approximation

The function
[`compute_detection_probability()`](https://larry-leon.github.io/forestsearch/reference/compute_detection_probability.md)
provides an analytical approximation based on asymptotic normal theory:

``` r
#| label: theoretical-power
#| fig-width: 8
#| fig-height: 6

# =============================================================================
# Theoretical Detection Probability Analysis
# =============================================================================

# Calculate expected subgroup characteristics
n_sg_expected <- sim_config_alt$n_sample * mean(dgm_calibrated$df_super_rand$flag.harm)
prop_cens <- mean(results_alt$p.cens)  # Censoring proportion

cat("=== Subgroup Characteristics ===\n")
```

    === Subgroup Characteristics ===

``` r
cat("Expected subgroup size (n_sg):", round(n_sg_expected), "\n")
```

    Expected subgroup size (n_sg): 89 

``` r
cat("Censoring proportion:", round(prop_cens, 3), "\n")
```

    Censoring proportion: 0.454 

``` r
cat("True HR in H:", round(dgm_calibrated$hr_H_true, 3), "\n")
```

    True HR in H: 2 

``` r
cat("HR threshold:", fs_params$hr.threshold, "\n")
```

    HR threshold: 1.25 

``` r
# -----------------------------------------------------------------------------
# Single-Point Detection Probability
# -----------------------------------------------------------------------------

# True H is dgm_calibrated$hr_H_true
# However we want at plim of observed estimate
#plim_hr_hatH <- c(summary_alt[c("hat(hat[H])"),1])

dgm_calibrated$hr_H_true
```

       treat
    1.999999 

``` r
# Compute detection probability at the true HR
prob_detect <- compute_detection_probability(
 theta = dgm_calibrated$hr_H_true,
  n_sg = round(n_sg_expected),
  prop_cens = prop_cens,
  hr_threshold = fs_params$hr.threshold,
  hr_consistency = fs_params$hr.consistency,
  method = "cubature"
)

# Compare theoretical to empirical (alternative)
cat("\n=== Detection Probability Comparison ===\n")
```

    === Detection Probability Comparison ===

``` r
cat("Theoretical FS (asymptotic):", round(prob_detect, 3), "\n")
```

    Theoretical FS (asymptotic): 0.899 

``` r
cat("Empirical FS:", round(mean(results_alt[analysis == "FS"]$any.H), 3), "\n")
```

    Empirical FS: 0.89 

``` r
cat("Empirical FSlg:", round(mean(results_alt[analysis == "FSlg"]$any.H), 3), "\n")
```

    Empirical FSlg: NaN 

``` r
if ("GRF" %in% results_alt$analysis) {
  cat("Empirical GRF:", round(mean(results_alt[analysis == "GRF"]$any.H), 3), "\n")
}
```

    Empirical GRF: 0.77 

``` r
# Null 

#plim_hr_itt <- c(summary_alt[c("hat(ITT)all"),1])

# Calculate at min SG size
# Compute detection probability at the true HR
prob_detect_null <- compute_detection_probability(
 theta = dgm_null$hr_causal,
  n_sg = fs_params$n.min,
  prop_cens = prop_cens,
  hr_threshold = fs_params$hr.threshold,
  hr_consistency = fs_params$hr.consistency,
  method = "cubature"
)


# Compare theoretical to empirical (alternative)
cat("\n=== Detection Probability Comparison ===\n")
```

    === Detection Probability Comparison ===

``` r
cat("Under the null calculate at min SG size:", fs_params$n.min,"\n")
```

    Under the null calculate at min SG size: 60 

``` r
cat("Theoretical FS at min(SG) (asymptotic):", round(prob_detect_null, 6), "\n")
```

    Theoretical FS at min(SG) (asymptotic): 0.039765 

``` r
cat("Empirical FS:", round(mean(results_null[analysis == "FS"]$any.H), 6), "\n")
```

    Empirical FS: 0.07 

``` r
cat("Empirical FSlg:", round(mean(results_null[analysis == "FSlg"]$any.H), 6), "\n")
```

    Empirical FSlg: NaN 

``` r
if ("GRF" %in% results_null$analysis) {
  cat("Empirical GRF:", round(mean(results_null[analysis == "GRF"]$any.H), 6), "\n")
}
```

    Empirical GRF: 0.1 

``` r
prop_cens <- mean(results_null$p.cens)  # Censoring proportion
cat("Censoring proportion:", round(prop_cens, 3), "\n")
```

    Censoring proportion: 0.467 

``` r
# -----------------------------------------------------------------------------
# Generate Full Detection Curve
# -----------------------------------------------------------------------------

# Generate detection probability curve across HR values
detection_curve <- generate_detection_curve(
  theta_range = c(0.5, 3.0),
  n_points = 50,
  n_sg = round(n_sg_expected),
  prop_cens = prop_cens,
  hr_threshold = fs_params$hr.threshold,
  hr_consistency = fs_params$hr.consistency,
  include_reference = TRUE,
  verbose = FALSE
)

# -----------------------------------------------------------------------------
# Visualization
# -----------------------------------------------------------------------------

# Plot detection curve with empirical overlay
plot_detection_curve(
  detection_curve,
  add_reference_lines = TRUE,
  add_threshold_line = TRUE,
  title = sprintf(
    "Detection Probability Curve (n=%d, cens=%.0f%%, threshold=%.2f)",
    round(n_sg_expected), 100 * prop_cens, fs_params$hr.threshold
  )
)

# Add empirical results as points
empirical_rates <- c(
  FS = mean(results_alt[analysis == "FS"]$any.H),
  FSlg = mean(results_alt[analysis == "FSlg"]$any.H)
)
if ("GRF" %in% results_alt$analysis) {
  empirical_rates["GRF"] <- mean(results_alt[analysis == "GRF"]$any.H)
}

# Mark the true HR and empirical detection rates
points(
  x = rep(dgm_calibrated$hr_H_true, length(empirical_rates)),
  y = empirical_rates,
  pch = c(16, 17, 18)[1:length(empirical_rates)],
  col = c("blue", "darkgreen", "purple")[1:length(empirical_rates)],
  cex = 1.5
)

# Add vertical line at true HR
abline(v = dgm_calibrated$hr_H_true, lty = 2, col = "blue", lwd = 1)

# Legend for empirical points
legend(
  "topleft",
  legend = c(
    sprintf("H true = %.2f", dgm_calibrated$hr_H_true),
    paste(names(empirical_rates), "=", round(empirical_rates, 3))
  ),
  pch = c(NA, 16, 17, 18)[1:(length(empirical_rates) + 1)],
  lty = c(2, rep(NA, length(empirical_rates))),
  col = c("blue", "blue", "darkgreen", "purple")[1:(length(empirical_rates) + 1)],
  cex = 0.8,
  bty = "n"
)
```

![](data:image/png;base64,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BCAAAQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAACEEAARh+AAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAhCoFAEEYJVqTioDAQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAACEIAABCCAAIw+AAEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAAQgUCkCCMAq1ZxUBgIQgAAEIAABCEAAAhCAAAQgAAEIQAACEIAABCAAAQhAAAEYfQACEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQKBSBCavVG2oDAQgAIEBIvDCCy+YCy64oFDuI0eONFNPPbX7m2uuucyiiy5qZp555kJxCdQ9As8995y56KKLMhkecMABZvjw4Rn/oeDx8ssvm3PPPbfPqk422WSOke/rM844o5lzzjnNYostZqaccso+4/dqgDvuuMP85S9/aSje7LPPbnbeeecGv165KLu8RZ+RouHErZWwvcK5U+V4//33zZ///Gdz/fXXm4cffti89tpr5vXXXzdf+tKXzMILL2wWWmghd/z6179uZptttk5lSzolErjyyivNDjvsUM/hvPPOM1tuuWX9uionn376qXn00UfNE088YZ588kmjMdEcc8xh5p9/fjPffPO5fjts2LCWqlur1cy///1v849//MM89NBD5t133zVf/OIXzXLLLWdWXnllo+8MDgJlELjssssy3/pf/OIXmfHLBx98YN544422izDVVFMZjZNSTuOun/3sZ+aBBx5wz9jo0aPNIossYjbddFP3l4qT8tN3Q2nIfeELX3BpTTvttKmgPen3zDPPmIsvvjhTtoMPPjjjh0eawGBmWObY8NlnnzXjxo3LQDvwwAN75jtz//33m7POOquhjMcdd5yZYYYZGvw6efHUU0+Zf/3rX0bHp59+2swyyyxmgQUWMF/5ylfccYoppmgpuxdffNHou96q01jCu4MOOsicfPLJ7nLEiBHmnnvucWMMf58jBCAAgaYE7EsIBwEIQAACbRKwE80a0fX7b8kll6xdeOGFtQ8//LDNknQu+vjx42tHHnlk7eabb+5coj2SUpG63XLLLcn27KU26jZO+0MjyaRo37eTLrUf//jHtZdeeqnbRe9IfvbHZqb+Sy+9dEfSLiORsstb9BkpGk4MWgkbMyvyXMdxeuHaCr7cu3a66abL9K/Us2UnLWsnnXRS7eOPP+6F4lOGHAJWYFOzEzf1NrULXWoTJkzICT14vW+99dbaUkstVa9nqs8uu+yyNY2TirqbbrqpNmbMmNw0rVC4dv7559es4K1okoSDQCECjz32WM0u1Mn0vTfffDMTX+P2VH9v1e/b3/52Jm15/O53v6vNNNNMuXko3jvvvJOMG3ped911DWkcdthh4e0BPy/y7f7Tn/7UUAfPmHdA8eYbzAzbGRv2RUi/c31/Co8fffRRX1G7dt8KsDNlfP7550vJ/+67766ts846mfxCNnbxX80KDQvnbxdzNU0vTDs8twspG/LQ+9kupqmnpXLiIAABCBQlwNI5+4bFQQACEBhoAlqVqVXiVhDmVk8PdHmuuuoqp21wyCGHGK1wrZKrct16vZ3eeustc/zxx5sllljCabj0enkp3+AhMFifa737tZpW71o7kVkIuDRhfvSjH5llllnG/P3vfy8Uh0DdJ3D44Yc7TSif80477WRaXTHt4/bqcffddzfSLNHq9GZOq7SltfWDH/ygWTA33thjjz3MWmut5bRB8wJrJfmOO+5oNt54Y2OFinnB8IdASwQ++eQTs/322xsrkGkpXhmB1ce33nprpwWcl/6vf/1ro/dMX84KvOpBpP2111571a8H+mSwfrsHmhv5Dx0CV199tbELTbpS4V133dVpWf/xj39smp80U7fbbjvzjW98w9hFjU3D6uaDDz7YZ5giAaSBtsYaa9SDqpzStMdBAAIQKEIAAVgRSoSBAAQg0CUC//znP93A8/bbb+9Sjo3ZPPLII25gudlmmxmZyqiSq3LdBls7ybTbqquuatTfcRBoh8Bgfq5l6lCTBzIX1x+nCQUJFezq5f5EJ06JBP7f//t/RibTvJO5vv/93//1l5U4yuzzGWecUbgudnWmOfXUU5OmpnwiYvTLX/6ysJmka6+91myyySZGggscBNolIJNivbKo4NBDD20QxH3rW99yC+ROP/10M/nkn+9ioWuZh8tzekbuvffe+m0Jv6affvr69UCdDOZv90AxI9+hR+C2225zQvlu1PyUU07JmFnsK18J5rRoQN/3Zs6bX20Wpug9LbwJnRaEyYQ4DgIQgEBfBBCA9UWI+xCAAAS6TMCaWTFbbLGFefXVV7ucs3GaOfEeR10vREkZSuuoqnUrCVmpyWrfDGsO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paa2CBAgQIAAAQIECBAgQIAAAQIECBAgQIAAAQIEahcQAKv9FhgAAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQIBAlQICYFVqaosAAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQKB2AQGw2m+BARAgQIAAAQIECBAgQIAAAQIECBAgQIAAAQIECFQpIABWpaa2CBAgQIAAAQIECBAgQIAAAQIECBAgQIAAAQIEahcQAKv9FhgAAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQIBAlQICYFVqaosAAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQKB2AQGw2m+BARAgQIAAAQIECBAgQIAAAQIECBAgQIAAAQIECFQpIABWpaa2CBAgQIAAAQIECBAgQIAAAQIECBAgQIAAAQIEahcQAKv9FhgAAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQIBAlQICYFVqaosAAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQKB2AQGw2m+BARAgQIAAAQIECBAgQIAAAQIECBAgQIAAAQIECFQpIABWpaa2CBAgQIAAAQIECBAgQIAAAQIECBAgQIAAAQIEahcQAKv9FhgAAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQIBAlQICYFVqaosAAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQKB2AQGw2m+BARAgQIAAAQIECBAgQIAAAQIECBAgQIAAAQIECFQpIABWpaa2CBAgQIAAAQIECBAgQIAAAQIECBAgQIAAAQIEahcQAKv9FhgAAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQIBAlQICYFVqaosAAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQKB2AQGw2m+BARAgQIAAAQIECBAgQIAAAQIECBAgQIAAAQIECFQpIABWpaa2CBAgQIAAAQIECBAgQIAAAQIECBAgQIAAAQIEahcQAKv9FhgAAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQIBAlQICYFVqaosAAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQKB2AQGw2m+BARAgQIAAAQIECBAgQIAAAQIECBAgQIAAAQIECFQpIABWpaa2CBAgQIAAAQIECBAgQIAAAQIECBAgQIAAAQIEahcQAKv9FhgAAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQIBAlQICYFVqaosAAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQKB2AQGw2m+BARAgQIAAAQIECBAgQIAAAQIECBAgQIAAAQIECFQpIABWpaa2CBAgQIAAAQIECBAgQIAAAQIECBAgQIAAAQIEahcQAKv9FhgAAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQIBAlQICYFVqaosAAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQKB2AQGw2m+BARAgQIAAAQIECBAgQIAAAQIECBAgQIAAAQIECFQpIABWpaa2CBAgQIAAAQIECBAgQIAAAQIECBAgQIAAAQIEahcQAKv9FhgAAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQIBAlQICYFVqaosAAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQKB2AQGw2m+BARAgQIAAAQIECBAgQIAAAQIECBAgQIAAAQIECFQpIABWpaa2CBAgQIAAAQIECBAgQIAAAQIECBAgQIAAAQIEahcQAKv9FhgAAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQIBAlQICYFVqaosAAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQKB2AQGw2m+BARAgQIAAAQIECBAgQIAAAQIECBAgQIAAAQIECFQpIABWpaa2CBAgQIAAAQIECBAgQIAAAQIECBAgQIAAAQIEahcQAKv9FhgAAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQIBAlQICYFVqaosAAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQKB2AQGw2m+BARAgQIAAAQIECBAgQIAAAQIECBAgQIAAAQIECFQpIABWpaa2CBAgQIAAAQIECBAgQIAAAQIECBAgQIAAAQIEahcQAKv9FhgAAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQIBAlQJTVtmYtggQIECAAAECBAgQIECAAIGBCTz00EPp+eefT8OHD0/Dhg0b2MmTae0bbrghffrpp2nttddOU07pq43J9GXgsgkQIECAAAECvQrIAOuVx0ECBAgQIECAAAECPQW+853vpNlnn32Cfy644IKejU+iexZeeOHS64033hjQVTZbjxs3bkDnqpzSjTfeWNoff/zx40Uy99xzl20034/u2/PMM09aYokl0iqrrJK23377NGbMmDTY92zs2LEtr/HSSy8tr+HMM89sWaeOna+88kpabbXV0qhRo9KQIUNaDqHdNV100UXlNZ177rktz61j5yeffJJefPHFjnZ9+eWXp/XWWy8de+yxHe1H4wQIECBAgAABAhOvgP8mNfHeOyMnQIAAAQIECBCoSeCtt95Kr7322gT3/uGHH05wGxNLA6+//np6++238+FG1sZAShXWA+lvUqsbAajC8P333x+vy4vzP/744wGd+5e//CVFoOmggw5KEahZYYUVBnT+QCvHtR166KEpAssvv/xyj9P/85//lA4ffPBBj+N17dh5553z90aMfY455ugyjL6uKT5Dinv7WbmmO+64I+2yyy5pu+22S7vuumuX66nyyWGHHZYiqHnwwQenDTfcMEWQXSFAgAABAgQIECDQLCAA1qxhmwABAgQIECBAgEA/BOKL1vfee69lzfjC+vHHH8+PTTPNNHk2TMuK2c5ZZ5213SH7CXxmBb73ve+lqaaaquX4Isj0wgsvpKeeeioVmX7PPPNMnqlz9913p3nnnbfleVXsHDFiRLr//vvTzDPPXEVzg9LGZZddlq6++uo011xzpT333LNHnxPbNV177bVp5MiRPa6jEzsiWPjjH/84D7Buu+226fbbb09TTDFFJ7rSJgECBAgQIECAwEQqIAA2kd44wyZAgAABAgQIEKhP4PTTT2/b+X333ZeWW265/PhXvvKVdM8997St6wCBiVHgjDPOSNNPP32vQ//vf/+bYkq+COpEsDiCYQcccEC68MILez1vQg6++eab+entgiBf//rXU0ybF2WppZbKH+v855133ikzpPbZZ5807bTT9hhOX9fU44Sad0R2bFHa3YfieBWPo0ePTieccEKKrLN4Xe64445VNKsNAgQIECBAgACBSUTAGmCTyI10GQQIECBAgAABAgQIEPisCESGWKwBduqpp5ZD+s1vflNOg1nuHMSNYcOGpY022ij/ieB03SUCNjE16EwzzZRPF1j3eCbG/v/v//6vDHodfvjhA56mc2K8ZmMmQIAAAQIECBDov4AAWP+t1CRAgAABAgQIECBAgACBAQjEdIlFtthHH32UnnzyyQGcPelWjQy5E088Mb/ADTbYIE099dST7sV2+Mo233zzvIeYejOCrAoBAgQIECBAgACBQsAUiIWERwIECBAgQIAAAQI1Ctx6663pkksuyUdw5JFHppgeLR5/+9vf5gGEb37zm/l0cvPMM0/afffdU6y1NNtss6VDDjmk7aj/9re/pdNOOy0/vt566+XrMLWqHF/Gn3POOfl0jQ8++GD6xz/+kWKds6WXXjqtueaaKb6gnxjKY489lq666qr0xBNP5D/PPvtsbrTgggum+IkvypdffvkelxLrVf3iF7/osb+vHTvttFPLqfTGdxzR31577ZVPGbjqqqum73//+/k0ghdccEGKNhdffPE822WTTTbpstbRuHHjUtSJ18/f//739PnPfz7F+THl35ZbbtnXZXT0+Oc+97l8rbtYGy9K8xR53TseX7c99tgjffDBB+lf//pX3mRMuVhMhbfQQgulOB7lr3/9azr//PPz7U033TStscYa+Xb3f6KtyFy766678vX8XnvttbToooumJZdcMq2++up5Bln3cwb6/OKLL04vvfRSflqMpXvp7zV1Py+e33LLLenKK6/M38/PPfdcWnbZZdPKK6+c1l577bTCCiu0OqXLvvH5PPjjH/+YB5/ivVSUeD0+/PDD+dO4nrgXzeWTTz5JY8aMSfGZU7xnwz6y86JufP7EfewrOBj3JepH38cee2z67ne/29yNbQIECBAgQIAAgclZoKEQIECAAAECBAgQIFCZwL333tvI/r7If5ZYYol+t3vyySeX52WBq0Y2XVv5vGjv7LPPztubYYYZ8mPZF8W9tn/11VeXbRx00EEt62ZBh8YyyyxT1iv6an7cbLPNGtlaRC3P7+/OGWecsezj1Vdf7e9peb3msWSBvx7nfvzxx40sWNjIvigv+2g+p9jOAkON7Iv4RvYFf5c2suBjr+cV53d/DN/mMqHjiLZmmWWWfCzbbrtt47jjjusxrnhdfPrpp2W3WUZVY7755utRrxhrFvhsZEGG8vhhhx1WnjuQjSmnnLJsIwsw9fvUp59+ujwvxpQFZHqcO6Fuza+t4rqLxywIWPZ33nnnlWM5/vjjy/3NG5deemnjy1/+clmvaKf5cf3112+88sorzacNeLt4z2XTHzayzLge5/f3ms4666xyrNlaWI399tuvka29Ve5rHne8/i+77LIefTXvGN/Pg/Bs7qv7drzHmkv0kwWjez0n2sgC8Y0sENl8asvtn/70p2Vbt912W8s6dhIgQIAAAQIECEx+AjLAst+qFQIECBAgQIAAAQKfJYHIAoosrCiRzROZEtNNN12KzJ8qS/ZFcVp33XXThx9+mDe73HLL5Vki8847b3rkkUfSjTfemE9ZlwUF0p133pkio2zo0KFVDqGStrIAXbriiivytmI9pZh2L7JIpplmmtwxsuhi7OEYmV6RERN1ipIFnXrNcosMlGgrTO677778tFh7KAtwFk3kjxM6jubGYqrACy+8sNwVmVRZ4CttscUWZfZXvEZWWWWV9M9//jOvF9f87W9/O80///zpgQceyLPhrr/++nT77beX7QzmRngfddRRZZdf+tKXUhbAK58XGxPqttVWW+Wv4V//+tfp3Xffze9VkfnWPeuo6LPV43XXXZdiLEVZccUV02qrrZYi6zJeP/E6iqytqBf3/v7778+PFfX7+5gFz/L7E/VHjhyZYr207mV8runAAw/M11jLAmB5ZmJkb2bBxfSHP/whPfroo/nrP7KjsgB6Wmeddbp3mSbk82D48OF5tla8biMDLUpkIRbvkTnnnLPsL/tPAvmxyFyMkgUp80zTeH1EhmBkAkaGXByP9uI1//jjj6d4D7QrG2+8cYo1wKLEfRoxYkS7qvYTIECAAAECBAhMTgKTX8zPFRMgQIAAAQIECBDonEAVGWDZ3yONBRZYoJEFLxqR8RQZEJGZUpQqMsAiCyqbUq/MmjjmmGMaWcCi6CJ/zAJjjSyQUNYZPXp0l+MDedKc0ZIF+BrZl/X9/gmP4qd7Blg2dVqZ8bLIIou0zFSLjKn999+/bOOrX/3qQIae180CCI3iGiIb6qabburSRlXjKDLAiuuNjLXImsqCXI0sINYYO3Zs2W8WEC2vacMNN2x0z8yK87Lp4co60eahhx5anj+QjYFmgEWGTzZdYJe+swBVjy6rcouGs8Bt3l8Ytiq9ZYCFbxaAyc/PAlKNLFDao4nIgozsr+LeZFOD9qjTnx3NGXknnXRSr6f0dU3NGWAxrniNZlMKdmkz3utZQK0cd9yX7qWqz4Nm48hqbVVGjRpVjiWbarJVlUYW3G3MMcccZb0s6NiyXrEzPruK12gWyC92eyRAgAABAgQIEJjMBdJkfv0unwABAgQIECBAgEClAlUEwLJMh8Ydd9zRdlxVBMDii/fii/ytt966bV9xYKWVVsrrRmAgy8TotW67g0XwqOhzfB+7B8B23nnn8jqyjLV23efTHhZu8UV58zSCbU/634GY7i6mHizGnK2r1uOUqsbRHADLsnR69FPsiPtQjCdbC66RrZ1UHOrymK0J1oip74q6VQTAYvq+CCJ2/8ky6xpZxlUZKCz6jMcddtihy7iKJ1W5RXt9BYuagzPdp0DMMgJLoz333LMYXo/HmK4wgtPFtUWQeqAl3m/F+d2n0ezeVl/X1BwAGzJkSKPd9H/Z+miNLIs07zfqdX/9V/V50GzcKgCWZdA14nMkrj/LsOt+uV2ex7SmhVN/XrfFezQ+P7NMsi5teUKAAAECBAgQIDB5CrSfQyD7TVMhQIAAAQIECBAgQGDwBbLAQj61XSd7PuWUU8rms3Whyu1WG1kWUr47yxJJF110Uasqte3bZZdd0uWXX56y9bLSWmut1XYcWdArnxowKsS0cFkgo23d5gPvv/9+ytbRKqekDIsdd9yxuUq+3Ylx7Lbbbj36KXbENG9F2X333dO0005bPO3yGNMiNk/r1+XgeD6J6RVjKsjuPzEl4FNPPZVPw1c0na1Plq655pp0+umnF7u6PHbCrUsH/Xxyww035DW/8IUvpH333bftWTFdYRaMKY/HtIsDLTGdYlGy9caKzQl+zDK72k79F1ODFtNBxms/popsLoP1eRDTL2bZk+nMM89MRx99dPMQemwX440D8T7sqxSWMVVoTFeqECBAgAABAgQIELAGmNcAAQIECBAgQIAAgc+YQAQtOllibaZnn3027yICFNlUY712t8IKK5THs4yicnt8N2JtsVh3q78l1iZrVxZbbLEUP72VWHMp1iUq1sqKuhEEi7W9eivhtPnmm5frfsV6TdlUkS1P6cQ4ensd3HrrreU4Yp2q3koE8KoMXC666KL52nTRZwQb/v3vf6fXXnstRYA0SqxBdcABB+TrrDUHMfKD3f7phFu3Lvp8+sYbb6RsesO8Xqz7NfPMM/d6TngWJZvCsdjs92P0V5RYX6yqEoHz3srcc8+dsukR8yrvvPNOinXsogzm50H0GYG6+GlX4r15zz33pCIoGfViX1+l2bLZuK/zHCdAgAABAgQIEJh0BQTAJt1768oIECBAgAABAgQmUoEFF1ywoyPP1tcpM6BefvnlFF+M97dUEQBbeumlU7beUn+77He9bJq39Oc//zll0wOmGGdkIz355JPp9ddf79FGNgFIj33dd0QGVjbFXb47m94vRbZPNr1a92o9nlcxjmzawhTByXYlgnpF6ev+9XW8aKe/jxGcmH766btUj+Bitm5WOuKII2Ka/XTiiSembEq6MuuoS+U2T6pwa9N0r7vj9VKUbMrBYrPtY2SJRZAsxhuvsbjeCPr1txTBmcjaG0gguK/2I8urtxJZkEWJwGVR6vw8iCzCyCgs3q/h+fTTT5efT8UY+/N+LTLA4pzCuDjfIwECBAgQIECAwOQp8P9/A548r99VEyBAgAABAgQIEPjMCVQdsOh+gfEFc1HGjRuXsnV5iqd9Pj7zzDN91hnsCs8//3z6+c9/ni644IKUrRHWsvvIDomsl/jpT4lMr2ytr7xq3I/rrruuR9CneztVjiMChDHdXrvS/AV/Xxl8zYGBdu1N6P4I5Bx++OEp+oppDd9+++20zTbb5FMzjho1qtfmq3TrtaM2ByPoUpT+WsVrIgJgMTVfPPY3kPXhhx+m+Iky55xzFt1W8jjNNNOMVzt1fB6MGTMmHXXUUW2nKozX/vzzz58HsPt7Uc2fm3FPFAIECBAgQIAAAQICYF4DBAgQIECAAAECBD5jAv3JMooh95UV0W7asObsnZjecIsttui3QG9BmX43UmHFyGAbMWJEuUZXNP3FL34xDR8+PC211FJpySWXTHGN8RhTKUbGSZTeMnYuu+yytPfee+f1Ys2iCH71Fayoehx9vQYiC6koEdTrLQDT7nVQnF/l40477ZRefPHFPBgW7W699dYppnJcfvnlW3ZTtVvLTvrY2Zw5FYG7/pTmAEsxlWB/zosgVbyHYrrICD5/Fspgfx5EsHr//fcvLz2yHeM1Eu/ReM/GT7ynb7/99rT++uvn9Xp7vxYNNXs2vz+K4x4JECBAgAABAgQmPwEBsMnvnrtiAgQIECBAgACBSUSgWHOp3eW0y3bqvi7Trrvu2q6Jz/z+zTbbrAx+xRpIv/rVr1JMV9iqxFpVRWmeAq7YF48xheIPfvCDPLgYQaiLL744xZSNfZWqx9FXf7PPPnt69NFH82qRQdVbACymuBvMcvDBB6ebb7453XXXXXlG3gYbbJDuvffelmvNDbZbK4fmKUfHjh3bqkqXffG+i8BdlKFDh6YhQ4Z0Od7bkwjkRHZfBAmb16Tr7ZxOHxvMz4N4XRTBr5iSMabM3G677XLH7tfZn/dr8znNnn1lRTafZ5sAAQIECBAgQGDSFeh7AvtJ99pdGQECBAgQIECAAIGJUqDIwiqmUmt3Ec8991x5qDlbbNZZZ00zzjhjfuzhhx/Os1HKii02Pvjgg3TTTTfl05H11WeL0zu2680338wDVtFBZOFExki74FdkQTWvm9UqKyrWIRo5cmQ5jWKsabXeeuv1Of6qx9Fnh1mFyJIpyiOPPFJstnxsnuKuZYWKd0Zg49xzz01TTz113nIEiyIzrHupw637GOJ5ZB8VGUZh2fxeaVU/6hQB1Dh3oKUIzsR0nTGFYt1lMD8PIpuyKAcddFDaa6+9Wga/ok4EdovS6v1aHCseBcAKCY8ECBAgQIAAAQKFgABYIeGRAAECBAgQIECAwEQiUEy5Fl/4Nk/F1jz8+ML4+uuvb97VZXvFFVfMn0dA66STTupyrPuTM844I62zzjppkUUWSVtuuWX3w7U9v+eee8pgxUorrZSmm266tmO59tpr03vvvVce/+STT8rt2AjLddddt/SMdax+9KMfdanT7kmV42jXR/f9G220Ubnr2GOPLR3Knf/biGDOySef3H13x58vuuiiKQIcRbnmmmvSFVdcUTzNH6t2i6n0onS/t/nOXv6Zdtppy4BirAfWfZzdT421zoqy4YYbFpv9fiwCYHFCc1C2VQPje02t2uptX1WfB8V4o69W9yGyAovyjW98o9js8fjRRx+lyy+/vNzfqq3y4P82mi37mrK0+7meEyBAgAABAgQITJoCAmCT5n11VQQIECBAgAABApOwQDElXwQ39txzz5ZXuu2226YHH3yw5bHYedRRR6Xiy+oDDzywnE6v+wlPPvlkijV7ivLDH/6w2Kz9ca655irH8MQTT6T40rxViWkNu69zFtk3RYkgYKw19Mwzz+S7IqhxwgknFIf7fKxqHH121FQhAn7LLLNMviey+NoFMU877bTUV4ZYU7OVbkZ2T6zrVJTddtstNa+xVbVbsZbVu+++W2bxFX339XjccceVVfbdd9/00ksvlc+bNyKDqQiQRYbbqFGjmg/3a3uVVVYp6zUHhMqdTRsTck1NzfS5WdXnQTHe6PD111/v0W/zPX/ooYd6HI8dEewK12K9vtjX/H6N561KYRnBr2HDhrWqYh8BAgQIECBAgMBkJiAANpndcJdLgAABAgQIECAw8QvssMMO5UWcf/75+fRyMUVhTHkYU89FMCf295YFEYGJnXfeOW8nMqOWW265FJlEMQ1gBNYiUBFtrb766mVW1MYbb5yav7wvB1HTRmQZxVpYUWKdq0033TSfqnHcuHH5FHX3339/imkMv/Wtb6WYxrG5NE+Xts0226Q777wzPzzbbLOlTTbZJF1yySXpzDPPTKecckqeQRVZVN1/Lr300vycqsbRPL6+tmPKvhhfEcSMbLVYyy3WBYvp+SIgOHr06BSZbMX0fn21WfXxmArxrLPOSrGWWpTI0PnJT35SdlO1W6ytFSWyH2Pqyl/+8pfpnHPOKfvrbWONNdZI8fqOElNGRpD5wgsvzNfqioBMBBn33nvvfIrMYorEo48+Os0///y9NdvyWLw/i3LHHXcUmy0fJ+SaWjbYZmdVnwfFeKObCCLHel8RXIzPlSjxeVKU8DzxxBPzz63YF+/Jq666KsWacVdeeWVRLX9sfr92OfC/JxGoL+qEb12v+VZjs48AAQIECBAgQKBGgeyXd4UAAQIECBAgQIAAgYoE7r333kb2633+s8QSS/S71Sy4Up6XBQ16PS8LcDSyL3nL+kV/zY8LLbRQI8toKutk09H1aDMLfDW23377sk5x/pAhQ3rsW3PNNRtZFkaPNvq7I1tzrGzz1Vdf7e9peb1iXPHYfQy/+93vGtmX3WXbUWeGGWZoDB06tNyXBWIaWeClkQX4yn3ZtI7lGOadd95yf3NffW1nQZKyjSrGEY3NMsss+Vjmnnvusu3eNrJp4hrTTDNNl/HH9RZjz6b3a5x66qnl88MOO6y35toea24zXjcDKVlwruw/7tWtt95anl6VWzTY/B4qrj9b36rs67zzzivHcfzxx5f7i40sgNLIppYs6xRtNF977MsCeo2DDz64OG28HhdYYIG8n74+I/q6pvisKMYZdXsrWWCprDt27NgeVav4PMgCzY155pmn7KcYWxZIzvvLgomNESNG9DieBeu7vI8XX3zxxp/+9KdGFlDL62aZY404t11pdsimfm1XzX4CBAgQIECAAIHJTEAGWPYbuUKAAAECBAgQIEBgYhKI7IZYU+mYY45JWcCky9BjfbCRI0emyCzpaxqwmK4s1vfKghBp+PDhZTZR81SC2Rf1eRbNDTfckGLKt89aWWuttdItt9ySll122XJokdH2zjvvpCyQl2ecRCZYTPEWmSVFGTNmTLFZyWNd44i1wOJeR5ZbUSIDKsqCCy6Y4r41HyvqDOZjFnRLWVAk7zL7eztlQddySrsq3SLbbb/99kvNWUhvvPFGmcHY1zXPPPPM+bpTF110UYrXfZFFVHhGtl1kF912223pZz/7WV/N9Xq8yAKLjL0XX3yxbd0Jvaa2Dbc4UMXnQayn9vvf/z7PKJ1qqqnKXh5//PF8O7IBY23CffbZJ2WB2/L4yy+/nGeeRkbnEUcckR544IG06qqrpnh9RIkpKbPAab7d6p/IgI0S19Db2mKtzrWPAAECBAgQIEBg0hWYIgJ+k+7luTICBAgQIECAAAECk75ATP8X631FwCumMiumnBvolcc6OzF1XkxXlmXO5NO7ZZlI493eQPufkPrxZ83zzz+fr+P11ltvpUUWWSQtvPDCKcvemZBmB3xuneOIQErcv1jT7Gtf+1qK6RwnllK1W1jEWmAReGtel2ogHnF+rJ/2wgsvpFi7KsuqzN8XA2mjXd2Ysm+xxRbLp6s85JBD0gEHHNCuarm/imsqG+vHxoR+HsRUpE899VSKwGIEtrp/LkWgOqabjLX3ZppppjwI3z2g349h5lMfxv2JwH2sMxfTKioECBAgQIAAAQIEQkAAzOuAAAECBAgQIECAAAECBAgMssDmm2+eYh25bArO9Oyzz5YZZ4M8jIm+u1jrbY899sgzPiOYFkF7hQABAgQIECBAgEAImALR64AAAQIECBAgQIAAAQIECAyyQEzXGNMsZutxpZtvvnmQe590ujv77LPzi9lyyy0Fvyad2+pKCBAgQIAAAQKVCAiAVcKoEQIECBAgQIAAAQIECBAg0H+BWHdv4403zk848sgj+3+imqVArCcW66jFemKxrphCgAABAgQIECBAoFlAAKxZwzYBAgQIECBAgAABAgQIEBgkgVivKta/igywG2+8cZB6nTS6+eSTT9Lee++dX8xhhx2W5ptvvknjwlwFAQIECBAgQIBAZQICYJVRaogAAQIECBAgQIAAAQIECPRfYPbZZ08nn3xyfkIEcz799NP+nzyZ1zzvvPPSY489llZeeeU0evToyVzD5RMgQIAAAQIECLQSEABrpWIfAQIECBAgQIAAAQIECBAYBIFRo0al3XffPc0222zp7rvvHoQeJ40ubr/99rTmmmumCIR97nO+2pg07qqrIECAAAECBAhUKzBFIyvVNqk1AgQIECBAgAABAgQIECBAgAABAgQIECBAgAABAvUJ+G9S9dnrmQABAgQIECBAgAABAgQIECBAgAABAgQIECBAoAMCAmAdQNUkAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQIBAfQICYPXZ65kAAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQKADAgJgHUDVJAECBAgQIECAAAECBAgQIECAAAECBAgQIECAQH0CAmD12euZAAECBAgQIECAAAECBAgQIECAAAECBAgQIECgAwICYB1A1SQBAgQIECBAgAABAgQIECBAgAABAgQIECBAgEB9AgJg9dnrmQABAgQIECBAgAABAgQIECBAgAABAgQIECBAoAMCAmAdQNUkAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQIBAfQICYPXZ65kAAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQKADAgJgHUDVJAECBAgQIECAAAECBAgQIECAAAECBAgQIECAQH0CAmD12euZAAECBAgQIECAAAECBAgQIECAAAECBAgQIECgAwICYB1A1SQBAgQIECBAgAABAgQIECBAgACDDHMUAAAECElEQVQBAgQIECBAgEB9AgJg9dnrmQABAgQIECBAgAABAgQIECBAgAABAgQIECBAoAMCAmAdQNUkAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQIBAfQICYPXZ65kAAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQKADAgJgHUDVJAECBAgQIECAAAECBAgQIECAAAECBAgQIECAQH0CAmD12euZAAECBAgQIECAAAECBAgQIECAAAECBAgQIECgAwICYB1A1SQBAgQIECBAgAABAgQIECBAgAABAgQIECBAgEB9AgJg9dnrmQABAgQIECBAgAABAgQIECBAgAABAgQIECBAoAMCAmAdQNUkAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQIBAfQICYPXZ65kAAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQKADAgJgHUDVJAECBAgQIECAAAECBAgQIECAAAECBAgQIECAQH0CAmD12euZAAECBAgQIECAAAECBAgQIECAAAECBAgQIECgAwICYB1A1SQBAgQIECBAgAABAgQIECBAgAABAgQIECBAgEB9AgJg9dnrmQABAgQIECBAgAABAgQIECBAgAABAgQIECBAoAMCAmAdQNUkAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQIBAfQICYPXZ65kAAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQKADAgJgHUDVJAECBAgQIECAAAECBAgQIECAAAECBAgQIECAQH0CAmD12euZAAECBAgQIECAAAECBAgQIECAAAECBAgQIECgAwICYB1A1SQBAgQIECBAgAABAgQIECBAgAABAgQIECBAgEB9AgJg9dnrmQABAgQIECBAgAABAgQIECBAgAABAgQIECBAoAMCAmAdQNUkAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQIBAfQICYPXZ65kAAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQKADAgJgHUDVJAECBAgQIECAAAECBAgQIECAAAECBAgQIECAQH0CAmD12euZAAECBAgQIECAAAECBAgQIECAAAECBAgQIECgAwICYB1A1SQBAgQIECBAgAABAgQIECBAgAABAgQIECBAgEB9AgJg9dnrmQABAgQIECBAgAABAgQIECBAgAABAgQIECBAoAMCAmAdQNUkAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQIBAfQICYPXZ65kAAQIECBAgQIAAAQIECBAgQIAAAQIECBAgQKADAgJgHUDVJAECBAgQIECAAAECBAgQIECAAAECBAgQIECAQH0CAmD12euZAAECBAgQIECAAAECBAgQIECAAAECBAgQIECgAwL/D2xpubZBaQ4tAAAAAElFTkSuQmCC)

### 7.1 AHR Metrics in Results (New)

The aligned analysis functions now compute AHR estimates in addition to
Cox-based HRs:

``` r
# Check for AHR columns in results
ahr_cols <- grep("ahr", names(results_alt), value = TRUE)
cat("AHR columns in results:", paste(ahr_cols, collapse = ", "), "\n\n")
```

    AHR columns in results:  

``` r
if (length(ahr_cols) > 0) {
  # Summarize AHR estimates
  results_found <- results_alt[results_alt$any.H == 1, ]
  
  if (nrow(results_found) > 0 && "ahr.H.hat" %in% names(results_found)) {
    cat("AHR estimates (when subgroup found):\n")
    cat("  Mean AHR(H) estimated:", round(mean(results_found$ahr.H.hat, na.rm = TRUE), 3), "\n")
    cat("  Mean AHR(Hc) estimated:", round(mean(results_found$ahr.Hc.hat, na.rm = TRUE), 3), "\n")
    cat("  True AHR(H):", round(dgm_calibrated$AHR_H_true, 3), "\n")
    cat("  True AHR(Hc):", round(dgm_calibrated$AHR_Hc_true, 3), "\n")
  }
}
```

### 7.2 Formatted Tables

``` r
# Format operating characteristics for H1
format_oc_results(
  results = results_alt,
  title = "Operating Characteristics (Alternative Hypothesis)",
  subtitle = sprintf("n = %d, %d simulations, HR(H) = %.2f",
                     sim_config_alt$n_sample,
                     sim_config_alt$n_sims,
                     dgm_calibrated$hr_H_true),
  use_gt = TRUE
)
```

[TABLE]

``` r
# Format operating characteristics for H0
format_oc_results(
  results = results_null,
  title = "Type I Error (Null Hypothesis)",
  subtitle = sprintf("n = %d, %d simulations, HR(overall) = %.2f",
                     sim_config_null$n_sample,
                     sim_config_null$n_sims,
                     dgm_null$hr_causal),
  use_gt = TRUE
)
```

[TABLE]

### 7.3 Key Metrics

``` r
# Extract key metrics
cat("=== KEY OPERATING CHARACTERISTICS ===\n\n")
```

    === KEY OPERATING CHARACTERISTICS ===

``` r
cat("Alternative Hypothesis (H1):\n")
```

    Alternative Hypothesis (H1):

``` r
for (analysis in unique(results_alt$analysis)) {
  res <- results_alt[results_alt$analysis == analysis, ]
  cat(sprintf("  %s: Power = %.3f, Sens = %.3f, Spec = %.3f, PPV = %.3f\n",
              analysis,
              mean(res$any.H),
              mean(res$sens, na.rm = TRUE),
              mean(res$spec, na.rm = TRUE),
              mean(res$ppv, na.rm = TRUE)))
}
```

      FS: Power = 0.830, Sens = 0.880, Spec = 0.977, PPV = 0.869
      GRF: Power = 0.830, Sens = 0.880, Spec = 0.977, PPV = 0.869

``` r
cat("\nNull Hypothesis (H0):\n")
```

    Null Hypothesis (H0):

``` r
for (analysis in unique(results_null$analysis)) {
  res <- results_null[results_null$analysis == analysis, ]
  cat(sprintf("  %s: Type I Error = %.4f\n",
              analysis,
              mean(res$any.H)))
}
```

      FS: Type I Error = 0.0850
      GRF: Type I Error = 0.0850

## 8 Using `format_oc_results()`

The
[`format_oc_results()`](https://larry-leon.github.io/forestsearch/reference/format_oc_results.md)
function is a flexible tool for creating publication-quality tables from
simulation results. It accepts the raw `data.table` produced by
[`run_simulation_analysis()`](https://larry-leon.github.io/forestsearch/reference/run_simulation_analysis.md)
(or combined via `rbind` / `rbindlist` across simulations) and returns
either a `gt` table object or a plain `data.frame`.

### 8.1 Function Signature

``` r
format_oc_results(
  results,
  analyses = NULL,
  metrics = "all",
  digits = 3,
  digits_hr = 3,
  title = "Operating Characteristics Summary",
  subtitle = NULL,
  use_gt = TRUE
)
```

### 8.2 Key Parameters

| Parameter | Default | Description |
|----|----|----|
| `results` | (required) | A `data.table` or `data.frame` with columns `analysis`, `any.H`, `sens`, `spec`, `ppv`, `npv`, `hr.H.hat`, etc., as produced by [`run_simulation_analysis()`](https://larry-leon.github.io/forestsearch/reference/run_simulation_analysis.md). |
| `analyses` | `NULL` | Character vector of analysis labels to include, e.g., `c("FS", "FSlg", "GRF")`. When `NULL`, all unique values of `results$analysis` are included. |
| `metrics` | `"all"` | Which metric groups to show: `"detection"`, `"classification"`, `"hr_estimates"`, `"subgroup_size"`, or `"all"`. |
| `digits` | `3` | Decimal places for proportions (detection rate, sens, spec, PPV, NPV). |
| `digits_hr` | `3` | Decimal places for hazard ratio estimates. |
| `title` | `"Operating Characteristics Summary"` | Table title shown in the `gt` header. |
| `subtitle` | `NULL` | Optional subtitle (e.g., sample size and true HR). |
| `use_gt` | `TRUE` | If `TRUE` and the `gt` package is available, returns a styled `gt` table; otherwise returns a `data.frame`. |

### 8.3 Output Structure

When `metrics = "all"`, the returned table contains one row per analysis
method with the following columns:

| Column | Description |
|----|----|
| `Detection` | Proportion of simulations that identified any subgroup (`any.H`) |
| `Sen` | Mean sensitivity across all simulations |
| `Spec` | Mean specificity |
| `PPV` | Mean positive predictive value |
| `NPV` | Mean negative predictive value |
| `HR_H_hat` | Mean estimated HR in identified subgroup (when found) |
| `HR_Hc_hat` | Mean estimated HR in complement (when found) |
| `HR_ITT` | Mean overall (ITT) HR across all simulations |
| `Size_H_mean` | Mean size of identified subgroup (when found) |

### 8.4 Usage Patterns

**Pattern 1: Quick summary of a single simulation run**

``` r
# After running simulations
format_oc_results(results_alt, title = "H1 Results")
```

**Pattern 2: Compare specific analysis methods**

``` r
# Compare only FS and GRF
format_oc_results(
  results_alt,
  analyses = c("FS", "GRF"),
  title = "FS vs GRF Comparison"
)
```

**Pattern 3: Focus on classification metrics only**

``` r
format_oc_results(
  results_alt,
  metrics = "classification",
  title = "Classification Performance"
)
```

**Pattern 4: Extract as data.frame for custom processing**

``` r
# Get raw summary data.frame for further manipulation
oc_df <- format_oc_results(results_alt, use_gt = FALSE)

# Use for custom plotting or LaTeX export
```

**Pattern 5: Multiple scenarios side-by-side**

``` r
# Create tables for different sample sizes or HR targets
for (scenario in names(results_list)) {
  cat("\n", scenario, "\n")
  print(format_oc_results(
    results_list[[scenario]],
    title = scenario,
    subtitle = sprintf("n = %d", sample_sizes[[scenario]])
  ))
}
```

**Pattern 6: Combine with
[`summarize_simulation_results()`](https://larry-leon.github.io/forestsearch/reference/summarize_simulation_results.md)
for detailed output**

[`summarize_simulation_results()`](https://larry-leon.github.io/forestsearch/reference/summarize_simulation_results.md)
returns a `data.frame` with row names like `any.H`, `sens`, `spec`, etc.
(from
[`summarize_single_analysis()`](https://larry-leon.github.io/forestsearch/reference/summarize_single_analysis.md)),
and one column per analysis method. This is useful for programmatic
access to the raw summary statistics, while
[`format_oc_results()`](https://larry-leon.github.io/forestsearch/reference/format_oc_results.md)
is oriented toward presentation-quality tables:

``` r
# Detailed numeric summary (for computation)
summary_df <- summarize_simulation_results(results_alt)

# Presentation table (for reports)
format_oc_results(results_alt, title = "Publication Table")
```

## 9 Advanced Topics

### 9.1 Comparing Standard vs Two-Stage Algorithm

``` r
# Run simulations with two-stage algorithm for comparison
results_twostage <- foreach(
  sim = 1:100,
  .combine = rbind,
  .options.future = list(
    packages = c("forestsearch", "survival", "data.table"),
    seed = TRUE
  )
) %dofuture% {
  
  run_simulation_analysis(
    sim_id = sim,
    dgm = dgm_calibrated,
    n_sample = sim_config_alt$n_sample,
    confounders_base = confounders_base,
    run_fs = TRUE,
    run_fs_grf = FALSE,
    run_grf = FALSE,
    fs_params = fs_params_fast,  # use_twostage = TRUE
    verbose = FALSE
  )
}

# Compare detection rates
cat("Standard algorithm power:", round(mean(results_alt$any.H[results_alt$analysis == "FS"]), 3), "\n")
cat("Two-stage algorithm power:", round(mean(results_twostage$any.H), 3), "\n")
```

### 9.2 Adding Noise Variables

Test ForestSearch robustness by including irrelevant noise variables:

``` r
# Run simulations with noise variables
results_noise <- foreach(
  sim = 1:sim_config_alt$n_sims,
  .combine = rbind,
  .errorhandling = "remove",
  .options.future = list(
    packages = c("forestsearch", "survival", "data.table"),
    seed = TRUE
  )
) %dofuture% {
  
  run_simulation_analysis(
    sim_id = sim,
    dgm = dgm_calibrated,
    n_sample = sim_config_alt$n_sample,
    confounders_base = confounders_base,
    n_add_noise = 10,  # Add 10 noise variables
    run_fs = TRUE,
    fs_params = fs_params,
    verbose = FALSE
  )
}

# Compare detection rates
cat("Without noise:", round(mean(results_alt$any.H), 3), "\n")
cat("With 10 noise vars:", round(mean(results_noise$any.H), 3), "\n")
```

### 9.3 Sensitivity Analysis: Varying Parameters

``` r
# Test different HR thresholds
thresholds <- c(1.10, 1.25, 1.50)
results_by_thresh <- list()

for (thresh in thresholds) {
  results_by_thresh[[as.character(thresh)]] <- foreach(
    sim = 1:100,
    .combine = rbind,
    .options.future = list(
      packages = c("forestsearch", "survival", "data.table"),
      seed = TRUE
    )
  ) %dofuture% {
    
    run_simulation_analysis(
      sim_id = sim,
      dgm = dgm_calibrated,
      n_sample = sim_config_alt$n_sample,
      confounders_base = confounders_base,
      run_fs = TRUE,
      fs_params = modifyList(fs_params, list(hr.threshold = thresh)),
      verbose = FALSE
    )
  }
  results_by_thresh[[as.character(thresh)]]$threshold <- thresh
}

# Combine and summarize
combined <- rbindlist(results_by_thresh)
combined[, .(power = mean(any.H), ppv = mean(ppv, na.rm = TRUE)), 
         by = .(threshold, analysis)]
```

## 10 Saving Results

``` r
# Save simulation results for later use
save_simulation_results(
  results = results_alt,
  dgm = dgm_calibrated,
  summary_table = summary_alt,
  runtime_hours = as.numeric(runtime_alt) / 60,
  output_file = "forestsearch_simulation_alt.Rdata",
  # Include AHR metrics in saved output
  ahr_metrics = list(
    AHR_H_true = dgm_calibrated$AHR_H_true,
    AHR_Hc_true = dgm_calibrated$AHR_Hc_true,
    AHR = dgm_calibrated$AHR
  )
)

save_simulation_results(
  results = results_null,
  dgm = dgm_null,
  summary_table = summary_null,
  runtime_hours = as.numeric(runtime_null) / 60,
  output_file = "forestsearch_simulation_null.Rdata"
)
```

## 11 Computational Timing Summary

Tracking wall-clock time for every stage of a simulation study is
essential for planning larger experiments and for reporting
reproducibility information.

``` r
# Finalize total vignette time
timings$total <- (proc.time() - t_vignette_start)["elapsed"]

# ── Build timing data.frame ─────────────────────────────────────────────────
timing_df <- data.frame(
  Stage = c(
    "DGM creation (H1)",
    "Calibrate k_inter (Cox)",
    "Calibrate k_inter (AHR)",
    "Validate k_inter",
    "DGM creation (H0)",
    "Simulations H1",
    "Simulations H0",
    "Summarize H1",
    "Summarize H0",
    "Total vignette"
  ),
  Seconds = c(
    timings$dgm_creation,
    timings$calibrate_cox,
    timings$calibrate_ahr,
    timings$validation,
    timings$dgm_null,
    timings$sims_alt_elapsed,
    timings$sims_null_elapsed,
    timings$summarize_alt,
    timings$summarize_null,
    timings$total
  ),
  stringsAsFactors = FALSE
)

timing_df$Minutes <- timing_df$Seconds / 60
timing_df$Pct <- 100 * timing_df$Seconds / timings$total

# ── Present as gt table ─────────────────────────────────────────────────────
gt(timing_df) |>
  tab_header(
    title = "Computational Timing Summary",
    subtitle = sprintf(
      "%d H1 + %d H0 simulations, %d workers",
      sim_config_alt$n_sims,
      sim_config_null$n_sims,
      n_workers
    )
  ) |>
  fmt_number(columns = Seconds, decimals = 1) |>
  fmt_number(columns = Minutes, decimals = 2) |>
  fmt_number(columns = Pct, decimals = 1) |>
  cols_label(
    Stage   = "Stage",
    Seconds = "Time (sec)",
    Minutes = "Time (min)",
    Pct     = "% of Total"
  ) |>
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(rows = Stage == "Total vignette")
  ) |>
  tab_footnote(
    footnote = sprintf("Parallel backend: %d workers via future::multisession.",
                       n_workers),
    locations = cells_column_labels(columns = Seconds)
  )
```

| Computational Timing Summary |  |  |  |
|----|----|----|----|
| 100 H1 + 100 H0 simulations, 13 workers |  |  |  |
| Stage | Time (sec)¹ | Time (min) | % of Total |
| DGM creation (H1) | 0.0 | 0.00 | 0.0 |
| Calibrate k_inter (Cox) | 2.4 | 0.04 | 1.7 |
| Calibrate k_inter (AHR) | 1.0 | 0.02 | 0.7 |
| Validate k_inter | 0.2 | 0.00 | 0.2 |
| DGM creation (H0) | 0.0 | 0.00 | 0.0 |
| Simulations H1 | 67.4 | 1.12 | 49.6 |
| Simulations H0 | 50.6 | 0.84 | 37.2 |
| Summarize H1 | 0.0 | 0.00 | 0.0 |
| Summarize H0 | 0.0 | 0.00 | 0.0 |
| Total vignette | 136.1 | 2.27 | 100.0 |
| ¹ Parallel backend: 13 workers via future::multisession. |  |  |  |

#### 11.0.1 Per-Simulation Timing

If per-simulation timing is available in the results (via a
`time_elapsed` column or similar), it can be summarized to characterize
the distribution of runtimes across individual simulations:

``` r
# If run_simulation_analysis() stored per-sim elapsed time:
if ("time_elapsed" %in% names(results_alt)) {
  per_sim_stats <- results_alt[
    , .(
      mean_sec  = mean(time_elapsed, na.rm = TRUE),
      median_sec = median(time_elapsed, na.rm = TRUE),
      sd_sec    = sd(time_elapsed, na.rm = TRUE),
      min_sec   = min(time_elapsed, na.rm = TRUE),
      max_sec   = max(time_elapsed, na.rm = TRUE)
    ),
    by = analysis
  ]
  
  gt(per_sim_stats) |>
    tab_header(title = "Per-Simulation Timing (seconds)") |>
    fmt_number(columns = 2:6, decimals = 1)
}
```

Alternatively, the overall wall-clock time divided by the number of
simulations gives the effective per-simulation cost accounting for
parallel overhead:

``` r
cat("Effective per-simulation cost:\n")
```

    Effective per-simulation cost:

``` r
cat(sprintf("  H1: %.1f sec/sim (wall) across %d sims on %d workers\n",
            timings$sims_alt_elapsed / sim_config_alt$n_sims,
            sim_config_alt$n_sims, n_workers))
```

      H1: 0.7 sec/sim (wall) across 100 sims on 13 workers

``` r
cat(sprintf("  H0: %.1f sec/sim (wall) across %d sims on %d workers\n",
            timings$sims_null_elapsed / sim_config_null$n_sims,
            sim_config_null$n_sims, n_workers))
```

      H0: 0.5 sec/sim (wall) across 100 sims on 13 workers

## 12 Complete Example Script

Here’s a minimal self-contained script for running a simulation study:

``` r
# ===========================================================================
# Complete ForestSearch Simulation Study - Minimal Example (Aligned)
# ===========================================================================

library(forestsearch)
library(data.table)
library(survival)
library(foreach)
library(doFuture)

# --- Configuration ---
N_SIMS <- 500
N_SAMPLE <- 500
TARGET_HR_HARM <- 1.5

# --- Setup parallel processing ---
plan(multisession, workers = 4)
registerDoFuture()

# --- Create DGM ---
# Option 1: Calibrate to Cox-based HR
k_inter <- calibrate_k_inter(target_hr_harm = TARGET_HR_HARM, 
                             use_ahr = FALSE, verbose = TRUE)

# Option 2: Calibrate to AHR instead
# k_inter <- calibrate_k_inter(target_hr_harm = TARGET_HR_HARM, 
#                              use_ahr = TRUE, verbose = TRUE)

dgm <- create_gbsg_dgm(model = "alt", k_inter = k_inter, verbose = TRUE)

# Verify hazard ratios (new aligned output)
cat("\nDGM Summary:\n")
cat("  Cox HR(H):", round(dgm$hr_H_true, 3), "\n")
cat("  AHR(H):", round(dgm$AHR_H_true, 3), "\n")
cat("  Cox HR(Hc):", round(dgm$hr_Hc_true, 3), "\n")
cat("  AHR(Hc):", round(dgm$AHR_Hc_true, 3), "\n")

# --- Run simulations ---
confounders <- c("v1", "v2", "v3", "v4", "v5", "v6", "v7")

results <- foreach(
  sim = 1:N_SIMS, 
  .combine = rbind,
  .options.future = list(
    packages = c("forestsearch", "survival", "data.table"),
    seed = TRUE
  )
) %dofuture% {
  run_simulation_analysis(
    sim_id = sim,
    dgm = dgm,
    n_sample = N_SAMPLE,
    max_follow = 60,
    confounders_base = confounders,
    run_fs = TRUE,
    run_fs_grf = TRUE,
    run_grf = FALSE,
    fs_params = list(
      hr.threshold = 1.25, 
      fs.splits = 300, 
      maxk = 2,
      use_twostage = FALSE  # Set TRUE for faster analysis
    )
  )
}

# --- Summarize ---
summary_table <- summarize_simulation_results(results)
print(summary_table)

# --- Display formatted table ---
format_oc_results(results = results, title = sprintf("Operating Characteristics (n=%d)", N_SAMPLE))

# --- Report AHR metrics (new) ---
results_found <- results[results$any.H == 1, ]
if (nrow(results_found) > 0 && "ahr.H.hat" %in% names(results_found)) {
  cat("\nAHR Estimates:\n")
  cat("  True AHR(H):", round(dgm$AHR_H_true, 3), "\n")
  cat("  Mean estimated AHR(H):", round(mean(results_found$ahr.H.hat, na.rm = TRUE), 3), "\n")
}
```

## 13 Summary

This vignette demonstrated the complete workflow for evaluating
ForestSearch performance through simulation:

| Step | Function | Purpose |
|----|----|----|
| 1\. Create DGM | [`create_gbsg_dgm()`](https://larry-leon.github.io/forestsearch/reference/create_gbsg_dgm.md) | Define data generating mechanism |
| 2\. Calibrate | [`calibrate_k_inter()`](https://larry-leon.github.io/forestsearch/reference/calibrate_k_inter.md) | Achieve target subgroup HR (Cox or AHR) |
| 3\. Validate | [`validate_k_inter_effect()`](https://larry-leon.github.io/forestsearch/reference/validate_k_inter_effect.md) | Verify heterogeneity control |
| 4\. Simulate | [`simulate_from_gbsg_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_gbsg_dgm.md) | Generate trial data |
| 5\. Analyze | [`run_simulation_analysis()`](https://larry-leon.github.io/forestsearch/reference/run_simulation_analysis.md) | Run ForestSearch/GRF |
| 6\. Summarize | [`summarize_simulation_results()`](https://larry-leon.github.io/forestsearch/reference/summarize_simulation_results.md) | Aggregate metrics |
| 7\. Display | [`format_oc_results()`](https://larry-leon.github.io/forestsearch/reference/format_oc_results.md) | Create gt tables |

**Key metrics to report:**

- **Power** (H1) / **Type I Error** (H0): Subgroup detection rate
- **Sensitivity**: P(identified \| true harm subgroup)
- **Specificity**: P(not identified \| true complement)
- **PPV**: P(true harm \| identified)
- **NPV**: P(true complement \| not identified)

**New aligned features:**

- **AHR metrics**: Alternative to Cox-based HR (from `loghr_po`)
- **`use_ahr` calibration**: Calibrate to AHR instead of Cox HR
- **`use_twostage`**: Faster two-stage search algorithm option
- **Individual effects**: Access `theta_0`, `theta_1`, `loghr_po` per
  subject

## 14 Appendix A: Publication-Quality Tables

This appendix demonstrates how to build publication-quality tables from
[`summarize_simulation_results()`](https://larry-leon.github.io/forestsearch/reference/summarize_simulation_results.md)
output that match the structure and content of the tables in Leon et
al. (2024). The two target table formats are:

1.  **Table 1 (Classification)**: Subgroup identification and
    classification rates across multiple data generation scenarios and
    analysis methods (cf. Table 4 of Leon et al.).
2.  **Table 2 (Estimation)**: Estimation properties including bias,
    coverage, and confidence interval metrics for the identified
    subgroup hazard ratios.

### 14.1 Package Functions for Table Construction

The `forestsearch` package provides three exported functions for
building publication-quality simulation tables. These are defined in
`R/simulation_tables.R` and available after loading the package:

- [`build_classification_table()`](https://larry-leon.github.io/forestsearch/reference/build_classification_table.md):
  Constructs a grouped `gt` table of subgroup identification and
  classification rates across scenarios, matching the layout of Table 4
  in Leon et al. (2024).
- [`build_estimation_table()`](https://larry-leon.github.io/forestsearch/reference/build_estimation_table.md):
  Summarizes estimation properties (Avg, SD, Min, Max, relative bias)
  for HR estimates in identified subgroups, matching Table 5 of Leon et
  al. (2024).
- [`render_reference_table()`](https://larry-leon.github.io/forestsearch/reference/render_reference_table.md):
  Renders a pre-assembled data frame of reference results as a styled
  `gt` table with consistent formatting.

See
[`?build_classification_table`](https://larry-leon.github.io/forestsearch/reference/build_classification_table.md),
[`?build_estimation_table`](https://larry-leon.github.io/forestsearch/reference/build_estimation_table.md),
and
[`?render_reference_table`](https://larry-leon.github.io/forestsearch/reference/render_reference_table.md)
for full documentation.

### 14.2 Generating Tables from Current Simulation Results

Using the package table functions, we construct tables from the
simulation results obtained in this vignette. These demonstrate the
table format; for a full replication of the Leon et al. (2024) tables,
one would run multiple scenarios (e.g., M1: N=700, M2: N=500, M3: N=300)
with 20,000 simulations each.

#### 14.2.1 Table 1: Classification Rates

``` r
# ── Assemble scenario list from the current vignette results ─────────────
scenario_list <- list(
  null = list(
    results    = results_null,
    label      = "M",
    n_sample   = sim_config_null$n_sample,
    dgm        = dgm_null,
    hypothesis = "null"
  ),
  alt = list(
    results    = results_alt,
    label      = "M",
    n_sample   = sim_config_alt$n_sample,
    dgm        = dgm_calibrated,
    hypothesis = "alt"
  )
)

# ── Build and display the classification table ───────────────────────────
build_classification_table(
  scenario_results = scenario_list,
  analyses = sort(unique(c(
    unique(results_null$analysis),
    unique(results_alt$analysis)
  ))),
  digits = 2,
  title = "Subgroup Identification and Classification Rates",
  n_sims = sim_config_alt$n_sims
)
```

| Subgroup Identification and Classification Rates |  |  |
|----|----|----|
| Across 100 simulations per scenario |  |  |
|  | FS | GRF |
| M Null: N=700, theta(ITT) = 0.72 |  |  |
| any(H) | 0.07 | 0.10 |
| sens(Hc) | 0.85 | 0.89 |
| ppv(Hc) | 1.00 | 1.00 |
| avg\|H\| | 106.00 | 75.00 |
| M Alt: N=700, p_H=13%, theta(H)=2, theta(Hc)=0.66, theta(ITT)=0.72 |  |  |
| any(H) | 0.89 | 0.77 |
| sens(H) | 0.85 | 0.92 |
| sens(Hc) | 0.98 | 0.97 |
| ppv(H) | 0.88 | 0.85 |
| ppv(Hc) | 0.98 | 0.99 |
| avg\|H\| | 86.00 | 98.00 |

#### 14.2.2 Table 2: Estimation Properties

``` r
# ── Build estimation table for the preferred analysis method ──────────────
# Uses the alternative-hypothesis results where subgroups are identified.
# If "FSlg" is not present, fall back to "FS".

est_analysis <- if ("FSlg" %in% unique(results_alt$analysis)) "FSlg" else "FS"

build_estimation_table(
  results = results_alt,
  dgm     = dgm_calibrated,
  analysis_method = est_analysis,
  digits  = 2,
  title   = "Estimation Properties for Identified Subgroup"
)
```

| Estimation Properties for Identified Subgroup |  |  |  |  |  |
|----|----|----|----|----|----|
| FS: 89 estimable realizations (B = 300 bootstraps) |  |  |  |  |  |
|  | Avg | SD | Min | Max | Rel Bias (%) |
| H: 89 estimable, avg \|H\| = 86, theta(H) = 2 |  |  |  |  |  |
| theta-hat(H-hat) | 2.38 | 0.92 | 1.37 | 8.16 | 19.17 |
| Hc: avg \|Hc\| = 614, theta(Hc) = 0.66 |  |  |  |  |  |
| theta-hat(Hc-hat) | 0.64 | 0.08 | 0.45 | 0.84 | -3.54 |

#### 14.2.3 Producing Multi-Scenario Tables (Full Replication)

To produce a table matching Table 4 of Leon et al. (2024) with three
model scenarios (M1, M2, M3), noise / no-noise comparisons, and six
analysis methods, the following pattern extends the above approach:

``` r
# ===========================================================================
# Full replication: M1 (N=700), M2 (N=500), M3 (N=300)
# Each with no-noise and 10-noise variants
# Six analysis methods: FS_l, FS_lg, GRF, GRF_60, VT(24), VT(36)
# ===========================================================================

# Assume results have been run and stored in a named list:
# all_results <- list(
#   "M1_null_nonoise" = list(results = ..., label = "M1", ...),
#   "M1_alt_nonoise"  = list(results = ..., label = "M1", ...),
#   "M1_null_noise"   = list(results = ..., label = "M1", ...),
#   ...
# )

# ── Build two-panel table (no noise | with noise) ──────────────────────
# Panel 1: No additional noise factors
table_nonoise <- build_classification_table(
  scenario_results = all_results[grep("nonoise", names(all_results))],
  analyses = c("FS", "FSlg", "GRF", "GRF60", "VT24", "VT36"),
  title = "Analysis with No Additional Noise Factors",
  n_sims = 20000
)

# Panel 2: With additional noise factors
table_noise <- build_classification_table(
  scenario_results = all_results[grep("noise$", names(all_results))],
  analyses = c("FS", "FSlg", "GRF", "GRF60", "VT24", "VT36"),
  title = "Analysis with Additional Noise Factors",
  n_sims = 20000
)

# Display
table_nonoise
table_noise

# ── For combined LaTeX output ────────────────────────────────────────────
# Use gt::as_latex() or gt::gtsave() to export:
# gt::gtsave(table_nonoise, filename = "table4_nonoise.tex")
# gt::gtsave(table_noise,   filename = "table4_noise.tex")
```

#### 14.2.4 Notes on Matching the Leon et al. (2024) Table Format

The LaTeX tables in the reference use specific formatting conventions:

| Feature | LaTeX Table | [`build_classification_table()`](https://larry-leon.github.io/forestsearch/reference/build_classification_table.md) Equivalent |
|----|----|----|
| Bold for inflated Type I error | `\bf{0.25}` | Use [`tab_style()`](https://gt.rstudio.com/reference/tab_style.html) with conditional logic on `any(H) > 0.05` |
| Superscript footnote markers | `$^{a}$` | Use [`tab_footnote()`](https://gt.rstudio.com/reference/tab_footnote.html) with [`cells_body()`](https://gt.rstudio.com/reference/cells_body.html) locations |
| Multi-column headers | `\multicolumn{6}{c}{}` | Use [`tab_spanner()`](https://gt.rstudio.com/reference/tab_spanner.html) in gt |
| Model scenario group rows | `\multicolumn{13}{c}{}` | Use `groupname_col` in gt |
| Horizontal rules | `\midrule` | Handled by [`tab_style()`](https://gt.rstudio.com/reference/tab_style.html) borders |

For the estimation properties table (Table 2), the three-row structure
per subgroup (oracle estimate, plugin estimate, bias-corrected estimate)
maps directly to the `make_est_row()` helper with different input
columns from the bootstrap results.

## 15 Appendix B: Reference Tables from Leon et al. (2024)

This appendix reproduces the published simulation results from Leon et
al. (2024) Tables 4 and 5 as `gt` tables, enabling direct visual
comparison with the vignette’s simulation output. The values below were
digitized from the LaTeX source of the original manuscript. All results
are based on 20,000 simulations (classification) or 1,000 simulations
(estimation).

The
[`render_reference_table()`](https://larry-leon.github.io/forestsearch/reference/render_reference_table.md)
function from the package accepts a pre-assembled data frame and returns
a styled `gt` object with the same visual conventions used by
[`build_classification_table()`](https://larry-leon.github.io/forestsearch/reference/build_classification_table.md).

### 15.1 Table B1: Classification Rates (Table 4 of Leon et al.)

``` r
# ── Digitized data from Tables.tex (no-noise panel) ─────────────────────
ref_class_nonoise <- data.frame(
  Scenario = c(
    # M1 Null
    rep("M1 Null: N=700, theta(ITT)=0.7", 4),
    # M1 Alt
    rep("M1 Alt: N=700, p_H=13%, theta(H)=2, theta(Hc)=0.65, theta(ITT)=0.71", 6),
    # M2 Null
    rep("M2 Null: N=500, theta(ITT)=0.69", 4),
    # M2 Alt
    rep("M2 Alt: N=500, p_H=20%, theta(H)=2, theta(Hc)=0.69, theta(ITT)=0.79", 6),
    # M3 Null
    rep("M3 Null: N=300, theta(ITT)=0.55", 4),
    # M3 Alt
    rep("M3 Alt: N=300, p_H=30%, theta(H)=2, theta(Hc)=0.56, theta(ITT)=0.74", 6)
  ),
  Metric = c(
    "any(H)", "sens(Hc)", "ppv(Hc)", "avg|H|",
    "any(H)", "sens(H)", "sens(Hc)", "ppv(H)", "ppv(Hc)", "avg|H|",
    "any(H)", "sens(Hc)", "ppv(Hc)", "avg|H|",
    "any(H)", "sens(H)", "sens(Hc)", "ppv(H)", "ppv(Hc)", "avg|H|",
    "any(H)", "sens(Hc)", "ppv(Hc)", "avg|H|",
    "any(H)", "sens(H)", "sens(Hc)", "ppv(H)", "ppv(Hc)", "avg|H|"
  ),
  FSl = c(
    0.02, 1, 1, 114,
    0.77, 0.72, 0.99, 0.69, 0.96, 94,
    0.02, 1, 1, 114,
    0.92, 0.84, 0.98, 0.84, 0.96, 102,
    0.00, 1, 1, 76,
    0.89, 0.73, 0.97, 0.80, 0.90, 82
  ),
  FSlg = c(
    0.03, 1, 1, 99,
    0.86, 0.82, 0.99, 0.80, 0.98, 92,
    0.03, 0.99, 1, 100,
    0.96, 0.88, 0.98, 0.88, 0.97, 101,
    0.00, 1, 1, 75,
    0.92, 0.78, 0.97, 0.84, 0.92, 83
  ),
  GRF = c(
    0.25, 0.97, 1, 88,
    0.94, 0.84, 0.97, 0.78, 0.98, 102,
    0.23, 0.96, 1, 87,
    0.98, 0.87, 0.94, 0.79, 0.97, 116,
    0.05, 0.99, 1, 74,
    0.97, 0.87, 0.93, 0.83, 0.95, 96
  ),
  GRF60 = c(
    0.05, 0.99, 1, 78,
    0.72, 0.66, 0.98, 0.61, 0.96, 99,
    0.05, 0.99, 1, 76,
    0.83, 0.73, 0.94, 0.66, 0.94, 115,
    0.01, 1, 1, 70,
    0.82, 0.72, 0.93, 0.68, 0.90, 97
  ),
  `VT(24)` = c(
    0.03, 1, 1, 78,
    0.49, 0.46, 0.99, 0.44, 0.93, 92,
    0.03, 1, 1, 76,
    0.66, 0.59, 0.98, 0.59, 0.91, 102,
    0.01, 1, 1, 70,
    0.61, 0.49, 0.97, 0.53, 0.83, 82
  ),
  `VT(36)` = c(
    0.04, 1, 1, 79,
    0.47, 0.42, 0.99, 0.41, 0.93, 93,
    0.04, 0.99, 1, 80,
    0.64, 0.56, 0.98, 0.56, 0.91, 103,
    0.02, 1, 1, 71,
    0.63, 0.52, 0.97, 0.55, 0.85, 86
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

render_reference_table(
  ref_class_nonoise,
  title = "Reference: Classification Rates (No Noise Factors)",
  subtitle = "Leon et al. (2024) Table 4, left panel --- 20,000 simulations"
)
```

| Reference: Classification Rates (No Noise Factors) |  |  |  |  |  |  |
|----|----|----|----|----|----|----|
| Leon et al. (2024) Table 4, left panel --- 20,000 simulations |  |  |  |  |  |  |
|  | FSl | FSlg | GRF | GRF60 | VT(24) | VT(36) |
| M1 Null: N=700, theta(ITT)=0.7 |  |  |  |  |  |  |
| any(H) | 0.02 | 0.03 | 0.25 | 0.05 | 0.03 | 0.04 |
| sens(Hc) | 1.00 | 1.00 | 0.97 | 0.99 | 1.00 | 1.00 |
| ppv(Hc) | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| avg\|H\| | 114.00 | 99.00 | 88.00 | 78.00 | 78.00 | 79.00 |
| M1 Alt: N=700, p_H=13%, theta(H)=2, theta(Hc)=0.65, theta(ITT)=0.71 |  |  |  |  |  |  |
| any(H) | 0.77 | 0.86 | 0.94 | 0.72 | 0.49 | 0.47 |
| sens(H) | 0.72 | 0.82 | 0.84 | 0.66 | 0.46 | 0.42 |
| sens(Hc) | 0.99 | 0.99 | 0.97 | 0.98 | 0.99 | 0.99 |
| ppv(H) | 0.69 | 0.80 | 0.78 | 0.61 | 0.44 | 0.41 |
| ppv(Hc) | 0.96 | 0.98 | 0.98 | 0.96 | 0.93 | 0.93 |
| avg\|H\| | 94.00 | 92.00 | 102.00 | 99.00 | 92.00 | 93.00 |
| M2 Null: N=500, theta(ITT)=0.69 |  |  |  |  |  |  |
| any(H) | 0.02 | 0.03 | 0.23 | 0.05 | 0.03 | 0.04 |
| sens(Hc) | 1.00 | 0.99 | 0.96 | 0.99 | 1.00 | 0.99 |
| ppv(Hc) | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| avg\|H\| | 114.00 | 100.00 | 87.00 | 76.00 | 76.00 | 80.00 |
| M2 Alt: N=500, p_H=20%, theta(H)=2, theta(Hc)=0.69, theta(ITT)=0.79 |  |  |  |  |  |  |
| any(H) | 0.92 | 0.96 | 0.98 | 0.83 | 0.66 | 0.64 |
| sens(H) | 0.84 | 0.88 | 0.87 | 0.73 | 0.59 | 0.56 |
| sens(Hc) | 0.98 | 0.98 | 0.94 | 0.94 | 0.98 | 0.98 |
| ppv(H) | 0.84 | 0.88 | 0.79 | 0.66 | 0.59 | 0.56 |
| ppv(Hc) | 0.96 | 0.97 | 0.97 | 0.94 | 0.91 | 0.91 |
| avg\|H\| | 102.00 | 101.00 | 116.00 | 115.00 | 102.00 | 103.00 |
| M3 Null: N=300, theta(ITT)=0.55 |  |  |  |  |  |  |
| any(H) | 0.00 | 0.00 | 0.05 | 0.01 | 0.01 | 0.02 |
| sens(Hc) | 1.00 | 1.00 | 0.99 | 1.00 | 1.00 | 1.00 |
| ppv(Hc) | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| avg\|H\| | 76.00 | 75.00 | 74.00 | 70.00 | 70.00 | 71.00 |
| M3 Alt: N=300, p_H=30%, theta(H)=2, theta(Hc)=0.56, theta(ITT)=0.74 |  |  |  |  |  |  |
| any(H) | 0.89 | 0.92 | 0.97 | 0.82 | 0.61 | 0.63 |
| sens(H) | 0.73 | 0.78 | 0.87 | 0.72 | 0.49 | 0.52 |
| sens(Hc) | 0.97 | 0.97 | 0.93 | 0.93 | 0.97 | 0.97 |
| ppv(H) | 0.80 | 0.84 | 0.83 | 0.68 | 0.53 | 0.55 |
| ppv(Hc) | 0.90 | 0.92 | 0.95 | 0.90 | 0.83 | 0.85 |
| avg\|H\| | 82.00 | 83.00 | 96.00 | 97.00 | 82.00 | 86.00 |

Reference classification rates from Leon et al. (2024) Table 4 —
Analysis without additional noise factors.

``` r
# ── Digitized data from Tables.tex (noise panel) ────────────────────────
ref_class_noise <- data.frame(
  Scenario = c(
    rep("M1 Null: N=700, theta(ITT)=0.7", 4),
    rep("M1 Alt: N=700, p_H=13%, theta(H)=2, theta(Hc)=0.65, theta(ITT)=0.71", 6),
    rep("M2 Null: N=500, theta(ITT)=0.69", 4),
    rep("M2 Alt: N=500, p_H=20%, theta(H)=2, theta(Hc)=0.69, theta(ITT)=0.79", 6),
    rep("M3 Null: N=300, theta(ITT)=0.55", 4),
    rep("M3 Alt: N=300, p_H=30%, theta(H)=2, theta(Hc)=0.56, theta(ITT)=0.74", 6)
  ),
  Metric = c(
    "any(H)", "sens(Hc)", "ppv(Hc)", "avg|H|",
    "any(H)", "sens(H)", "sens(Hc)", "ppv(H)", "ppv(Hc)", "avg|H|",
    "any(H)", "sens(Hc)", "ppv(Hc)", "avg|H|",
    "any(H)", "sens(H)", "sens(Hc)", "ppv(H)", "ppv(Hc)", "avg|H|",
    "any(H)", "sens(Hc)", "ppv(Hc)", "avg|H|",
    "any(H)", "sens(H)", "sens(Hc)", "ppv(H)", "ppv(Hc)", "avg|H|"
  ),
  FSl = c(
    0.02, 1, 1, 126,
    0.71, 0.64, 0.98, 0.60, 0.95, 96,
    0.03, 0.99, 1, 117,
    0.89, 0.77, 0.97, 0.77, 0.95, 103,
    0.00, 1, 1, 76,
    0.88, 0.68, 0.96, 0.76, 0.88, 80
  ),
  FSlg = c(
    0.11, 0.99, 1, 91,
    0.83, 0.74, 0.98, 0.71, 0.97, 93,
    0.14, 0.98, 1, 88,
    0.96, 0.81, 0.96, 0.81, 0.95, 101,
    0.02, 1, 1, 74,
    0.93, 0.71, 0.95, 0.78, 0.89, 81
  ),
  GRF = c(
    0.61, 0.92, 1, 94,
    0.94, 0.66, 0.93, 0.60, 0.95, 106,
    0.60, 0.89, 1, 89,
    0.99, 0.70, 0.88, 0.62, 0.92, 118,
    0.13, 0.97, 1, 74,
    0.96, 0.73, 0.88, 0.70, 0.89, 95
  ),
  GRF60 = c(
    0.27, 0.97, 1, 81,
    0.71, 0.52, 0.96, 0.47, 0.94, 101,
    0.32, 0.95, 1, 80,
    0.86, 0.58, 0.89, 0.51, 0.90, 119,
    0.07, 0.98, 1, 71,
    0.87, 0.62, 0.87, 0.59, 0.85, 96
  ),
  `VT(24)` = c(
    0.04, 1, 1, 79,
    0.44, 0.37, 0.99, 0.36, 0.92, 92,
    0.04, 0.99, 1, 77,
    0.56, 0.44, 0.97, 0.43, 0.88, 101,
    0.01, 1, 1, 70,
    0.51, 0.36, 0.95, 0.39, 0.79, 83
  ),
  `VT(36)` = c(
    0.06, 0.99, 1, 81,
    0.42, 0.34, 0.99, 0.33, 0.92, 93,
    0.06, 0.99, 1, 79,
    0.53, 0.40, 0.97, 0.40, 0.87, 102,
    0.02, 1, 1, 72,
    0.53, 0.37, 0.95, 0.40, 0.80, 85
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

render_reference_table(
  ref_class_noise,
  title = "Reference: Classification Rates (With Noise Factors)",
  subtitle = "Leon et al. (2024) Table 4, right panel --- 20,000 simulations"
)
```

| Reference: Classification Rates (With Noise Factors) |  |  |  |  |  |  |
|----|----|----|----|----|----|----|
| Leon et al. (2024) Table 4, right panel --- 20,000 simulations |  |  |  |  |  |  |
|  | FSl | FSlg | GRF | GRF60 | VT(24) | VT(36) |
| M1 Null: N=700, theta(ITT)=0.7 |  |  |  |  |  |  |
| any(H) | 0.02 | 0.11 | 0.61 | 0.27 | 0.04 | 0.06 |
| sens(Hc) | 1.00 | 0.99 | 0.92 | 0.97 | 1.00 | 0.99 |
| ppv(Hc) | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| avg\|H\| | 126.00 | 91.00 | 94.00 | 81.00 | 79.00 | 81.00 |
| M1 Alt: N=700, p_H=13%, theta(H)=2, theta(Hc)=0.65, theta(ITT)=0.71 |  |  |  |  |  |  |
| any(H) | 0.71 | 0.83 | 0.94 | 0.71 | 0.44 | 0.42 |
| sens(H) | 0.64 | 0.74 | 0.66 | 0.52 | 0.37 | 0.34 |
| sens(Hc) | 0.98 | 0.98 | 0.93 | 0.96 | 0.99 | 0.99 |
| ppv(H) | 0.60 | 0.71 | 0.60 | 0.47 | 0.36 | 0.33 |
| ppv(Hc) | 0.95 | 0.97 | 0.95 | 0.94 | 0.92 | 0.92 |
| avg\|H\| | 96.00 | 93.00 | 106.00 | 101.00 | 92.00 | 93.00 |
| M2 Null: N=500, theta(ITT)=0.69 |  |  |  |  |  |  |
| any(H) | 0.03 | 0.14 | 0.60 | 0.32 | 0.04 | 0.06 |
| sens(Hc) | 0.99 | 0.98 | 0.89 | 0.95 | 0.99 | 0.99 |
| ppv(Hc) | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| avg\|H\| | 117.00 | 88.00 | 89.00 | 80.00 | 77.00 | 79.00 |
| M2 Alt: N=500, p_H=20%, theta(H)=2, theta(Hc)=0.69, theta(ITT)=0.79 |  |  |  |  |  |  |
| any(H) | 0.89 | 0.96 | 0.99 | 0.86 | 0.56 | 0.53 |
| sens(H) | 0.77 | 0.81 | 0.70 | 0.58 | 0.44 | 0.40 |
| sens(Hc) | 0.97 | 0.96 | 0.88 | 0.89 | 0.97 | 0.97 |
| ppv(H) | 0.77 | 0.81 | 0.62 | 0.51 | 0.43 | 0.40 |
| ppv(Hc) | 0.95 | 0.95 | 0.92 | 0.90 | 0.88 | 0.87 |
| avg\|H\| | 103.00 | 101.00 | 118.00 | 119.00 | 101.00 | 102.00 |
| M3 Null: N=300, theta(ITT)=0.55 |  |  |  |  |  |  |
| any(H) | 0.00 | 0.02 | 0.13 | 0.07 | 0.01 | 0.02 |
| sens(Hc) | 1.00 | 1.00 | 0.97 | 0.98 | 1.00 | 1.00 |
| ppv(Hc) | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 | 1.00 |
| avg\|H\| | 76.00 | 74.00 | 74.00 | 71.00 | 70.00 | 72.00 |
| M3 Alt: N=300, p_H=30%, theta(H)=2, theta(Hc)=0.56, theta(ITT)=0.74 |  |  |  |  |  |  |
| any(H) | 0.88 | 0.93 | 0.96 | 0.87 | 0.51 | 0.53 |
| sens(H) | 0.68 | 0.71 | 0.73 | 0.62 | 0.36 | 0.37 |
| sens(Hc) | 0.96 | 0.95 | 0.88 | 0.87 | 0.95 | 0.95 |
| ppv(H) | 0.76 | 0.78 | 0.70 | 0.59 | 0.39 | 0.40 |
| ppv(Hc) | 0.88 | 0.89 | 0.89 | 0.85 | 0.79 | 0.80 |
| avg\|H\| | 80.00 | 81.00 | 95.00 | 96.00 | 83.00 | 85.00 |

Reference classification rates from Leon et al. (2024) Table 4 —
Analysis with additional noise factors.

### 15.2 Table B2: Estimation Properties (Table 5 of Leon et al.)

The estimation properties table reports bias, variability, and coverage
for the FSlg method across all three DGM scenarios. Three estimators are
shown per subgroup:

- **theta-hat(H)**: Oracle estimate using the true subgroup membership.
- **theta-hat(H-hat)**: Plugin estimate using the identified subgroup.
- **theta-hat\*(H-hat)**: Bootstrap bias-corrected estimate.

Bias metrics use three reference points: oracle, double-dagger, and
dagger (see Leon et al. Sections 4–5 for definitions).

``` r
# ── Digitized data from Tables.tex (estimation table) ───────────────────
ref_est <- data.frame(
  Scenario = c(
    # M1 H-hat
    rep("M1 H-hat: 839 estimable, avg |H|=89, \u03b8\u2020(H)=2, \u03b8\u2021(H)=2.25", 3),
    # M1 Hc-hat
    rep("M1 Hc-hat: avg |Hc|=611, \u03b8\u2020(Hc)=0.65, \u03b8\u2021(Hc)=0.60", 3),
    # M2 H-hat
    rep("M2 H-hat: 949 estimable, avg |H|=101, \u03b8\u2020(H)=2, \u03b8\u2021(H)=2.61", 3),
    # M2 Hc-hat
    rep("M2 Hc-hat: avg |Hc|=399, \u03b8\u2020(Hc)=0.69, \u03b8\u2021(Hc)=0.64", 3),
    # M3 H-hat
    rep("M3 H-hat: 924 estimable, avg |H|=90, \u03b8\u2020(H)=2, \u03b8\u2021(H)=2.56", 3),
    # M3 Hc-hat
    rep("M3 Hc-hat: avg |Hc|=210, \u03b8\u2020(Hc)=0.56, \u03b8\u2021(Hc)=0.49", 3)
  ),
  Estimator = rep(c(
    "$\\hat{\\theta}(H)$",
    "$\\hat{\\theta}(\\widehat{H})$",
    "$\\hat{\\theta}^{*}(\\widehat{H})$"
  ), 6),
  Avg = c(
    2.22, 2.18, 1.80,
    0.65, 0.65, 0.66,
    2.34, 2.39, 1.96,
    0.69, 0.69, 0.71,
    2.29, 2.48, 1.95,
    0.55, 0.59, 0.62
  ),
  SD = c(
    0.58, 0.53, 0.48,
    0.08, 0.08, 0.08,
    0.60, 0.58, 0.52,
    0.10, 0.10, 0.11,
    0.61, 0.62, 0.52,
    0.11, 0.13, 0.14
  ),
  SD_hat = c(
    0.57, 0.57, 0.53,
    0.07, 0.07, 0.11,
    0.57, 0.61, 0.59,
    0.10, 0.10, 0.14,
    0.60, 0.73, 0.69,
    0.11, 0.11, 0.17
  ),
  Min = c(
    1.06, 1.40, 1.07,
    0.44, 0.44, 0.45,
    1.10, 1.41, 1.11,
    0.43, 0.43, 0.44,
    1.00, 1.45, 1.11,
    0.25, 0.28, 0.28
  ),
  Max = c(
    6.20, 6.08, 4.82,
    0.99, 0.90, 0.92,
    5.75, 5.75, 4.95,
    1.01, 1.05, 1.12,
    6.97, 6.97, 5.83,
    1.10, 1.35, 1.41
  ),
  `Bias_oracle (%)` = c(
    0.00, -0.54, -18.55,
    0.00, -0.26, 1.41,
    0.00, 3.09, -15.95,
    0.00, 0.47, 3.49,
    0.00, 10.21, -13.66,
    0.00, 6.79, 12.62
  ),
  `Bias_ddagger (%)` = c(
    -1.12, 14.13, -6.28,
    8.05, 2.84, 4.55,
    -10.27, 8.99, -11.09,
    7.52, -1.82, 1.12,
    -10.64, 12.62, -11.58,
    11.31, -9.69, -4.76
  ),
  `Bias_dagger (%)` = c(
    11.21, 9.17, -10.04,
    0.93, 0.64, 2.33,
    17.17, 19.40, -2.05,
    0.04, 0.50, 3.56,
    14.34, 23.97, -2.39,
    -1.32, 5.14, 10.93
  ),
  Length = c(
    2.35, 2.32, 2.21,
    0.29, 0.29, 0.43,
    2.33, 2.48, 2.45,
    0.38, 0.38, 0.56,
    2.47, 3.04, 2.96,
    0.45, 0.45, 0.68
  ),
  `Cov_oracle` = c(
    1.00, 0.98, 0.95,
    1.00, 1.00, 1.00,
    1.00, 0.99, 0.99,
    1.00, 1.00, 1.00,
    1.00, 0.99, 1.00,
    1.00, 0.99, 1.00
  ),
  `Cov_ddagger` = c(
    0.97, 0.93, 0.87,
    0.89, 0.87, 0.96,
    0.93, 0.92, 0.85,
    0.92, 0.83, 0.94,
    0.94, 0.95, 0.89,
    0.92, 0.76, 0.93
  ),
  `Cov_dagger` = c(
    0.96, 0.97, 0.91,
    0.94, 0.94, 0.99,
    0.92, 0.93, 0.97,
    0.95, 0.94, 0.98,
    0.95, 0.95, 0.97,
    0.94, 0.92, 0.97
  ),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

# ── Render with gt ───────────────────────────────────────────────────────
gt_ref_est <- gt::gt(ref_est, groupname_col = "Scenario") |>
  gt::tab_header(
    title = "Reference: Estimation Properties (FSlg)",
    subtitle = "Leon et al. (2024) Table 5 --- 1,000 simulations, B=300 bootstraps"
  ) |>
  gt::cols_label(
    Estimator          = "",
    SD_hat             = gt::md("$\\widehat{\\text{SD}}$"),
    `Bias_oracle (%)`  = gt::md("$\\hat{b}^{\\text{oracle}}$ (%)"),
    `Bias_ddagger (%)` = gt::md("$\\hat{b}^{\\ddagger}$ (%)"),
    `Bias_dagger (%)`  = gt::md("$b^{\\dagger}$ (%)"),
    Cov_oracle         = gt::md("$\\hat{C}^{\\text{oracle}}$"),
    Cov_ddagger        = gt::md("$\\hat{C}^{\\ddagger}$"),
    Cov_dagger         = gt::md("$C^{\\dagger}$")
  ) |>
  gt::fmt_markdown(columns = Estimator) |>
  gt::tab_spanner(
    label = gt::md("Relative Bias (%)"),
    columns = c(`Bias_oracle (%)`, `Bias_ddagger (%)`, `Bias_dagger (%)`)
  ) |>
  gt::tab_spanner(
    label = "CI Coverage",
    columns = c(Cov_oracle, Cov_ddagger, Cov_dagger)
  ) |>
  gt::tab_style(
    style = gt::cell_text(size = "small"),
    locations = gt::cells_body()
  ) |>
  gt::tab_style(
    style = gt::cell_text(weight = "bold", size = "small"),
    locations = gt::cells_row_groups()
  ) |>
  gt::tab_style(
    style = gt::cell_text(size = "small"),
    locations = gt::cells_column_labels()
  ) |>
  gt::tab_style(
    style = gt::cell_text(size = "small"),
    locations = gt::cells_column_spanners()
  ) |>
  gt::tab_options(
    table.font.size = gt::px(11),
    row_group.padding = gt::px(4)
  ) |>
  gt::tab_footnote(
    footnote = gt::md(paste(
      "$\\hat{\\theta}(H)$ = oracle estimate using true subgroup;",
      "$\\hat{\\theta}(\\widehat{H})$ = plugin estimate using identified subgroup;",
      "$\\hat{\\theta}^{*}(\\widehat{H})$ = bootstrap bias-corrected estimate."
    ))
  ) |>
  gt::tab_footnote(
    footnote = gt::md(paste(
      "$\\hat{b}^{\\text{oracle}}$: bias relative to oracle HR;",
      "$\\hat{b}^{\\ddagger}$: bias relative to $\\theta^{\\ddagger}$ (expected HR in $\\widehat{H}$);",
      "$b^{\\dagger}$: bias relative to $\\theta^{\\dagger}$ (true causal HR)."
    ))
  )

gt_ref_est
```

[TABLE]

Reference estimation properties from Leon et al. (2024) Table 5 — FSlg
method, 1,000 simulations, B=300 bootstraps.

### 15.3 Comparing Vignette Results with Published Benchmarks

When running the full 20,000-simulation study, you can place the
vignette’s
[`build_classification_table()`](https://larry-leon.github.io/forestsearch/reference/build_classification_table.md)
output alongside the reference tables above to assess reproducibility.
Key comparisons:

- **Type I error control**: The `any(H)` row under null scenarios should
  be close to the reference values (e.g., FSl: 0.02, FSlg: 0.03 for M1).
- **Power**: The `any(H)` row under alternative scenarios (e.g., FSlg:
  0.86 for M1 Alt) serves as the detection benchmark.
- **Classification accuracy**: `sens(H)` and `ppv(H)` values indicate
  how well the method recovers the true subgroup boundary.
- **Estimation bias**: The bias-corrected estimator `theta-hat*(H-hat)`
  should reduce the upward bias seen in the plugin `theta-hat(H-hat)`.

Note that exact replication requires the same DGM calibration, censoring
rate, and random seeds. The vignette’s reduced simulation count
(`n_sims = 100`) will show higher Monte Carlo variability than the
published 20,000-simulation results.

## 16 Session Info

``` r
sessionInfo()
```

    R version 4.5.1 (2025-06-13)
    Platform: aarch64-apple-darwin20
    Running under: macOS Tahoe 26.2

    Matrix products: default
    BLAS:   /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRblas.0.dylib
    LAPACK: /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.1

    locale:
    [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

    time zone: America/Los_Angeles
    tzcode source: internal

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base

    other attached packages:
    [1] doFuture_1.2.0          future_1.69.0           foreach_1.5.2
    [4] gt_1.3.0                ggplot2_4.0.2           survival_3.8-6
    [7] data.table_1.18.2.1     weightedsurv_0.1.0      forestsearch_0.0.0.9000

    loaded via a namespace (and not attached):
     [1] sass_0.4.10          generics_0.1.4       xml2_1.5.2
     [4] shape_1.4.6.1        stringi_1.8.7        lattice_0.22-9
     [7] cubature_2.1.4-1     listenv_0.10.0       digest_0.6.39
    [10] magrittr_2.0.4       grf_2.5.0            evaluate_1.0.5
    [13] grid_4.5.1           RColorBrewer_1.1-3   iterators_1.0.14
    [16] policytree_1.2.3     fastmap_1.2.0        jsonlite_2.0.0
    [19] Matrix_1.7-4         glmnet_4.1-10        scales_1.4.0
    [22] codetools_0.2-20     cli_3.6.5            rlang_1.1.7
    [25] litedown_0.9         parallelly_1.46.1    future.apply_1.20.1
    [28] commonmark_2.0.0     splines_4.5.1        base64enc_0.1-6
    [31] withr_3.0.2          yaml_2.3.12          otel_0.2.0
    [34] tools_4.5.1          parallel_4.5.1       dplyr_1.2.0
    [37] globals_0.19.0       vctrs_0.7.1          R6_2.6.1
    [40] lifecycle_1.0.5      stringr_1.6.0        randomForest_4.7-1.2
    [43] fs_1.6.6             pkgconfig_2.0.3      progressr_0.18.0
    [46] pillar_1.11.1        gtable_0.3.6         glue_1.8.0
    [49] Rcpp_1.1.1           xfun_0.56            tibble_3.3.1
    [52] tidyselect_1.2.1     rstudioapi_0.18.0    knitr_1.51
    [55] farver_2.1.2         htmltools_0.5.9      rmarkdown_2.30
    [58] compiler_4.5.1       S7_0.2.1             markdown_2.0        

## 17 References

1.  Leon LF, Marceau-West CT, He W, et al. (2024). “Identifying Patient
    Subgroups with Differential Treatment Effects: A Forest Search
    Approach.” *Statistics in Medicine*.

2.  Athey S, Imbens GW. (2016). “Recursive partitioning for
    heterogeneous causal effects.” *PNAS*, 113(27):7353-7360.

3.  Wager S, Athey S. (2018). “Estimation and inference of heterogeneous
    treatment effects using random forests.” *JASA*, 113(523):1228-1242.
