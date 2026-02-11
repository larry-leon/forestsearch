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
  n_sims = 1000,          # Number of simulations (use 500-1000 for final)
  n_sample = 700,         # Sample size per trial
  max_follow = 84,        # Maximum follow-up (months)
  seed_base = 8316951,
  muC_adj = log(1.5)
)

sim_config_null <- list(
  n_sims = 1000,          # More simulations for Type I error estimation
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

    Running 1000 simulations under H1...

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
      return_selected_cuts_only = TRUE: using cuts from selected tree only
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
    Seconds and minutes forestsearch overall = 10.743 0.1791
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
      return_selected_cuts_only = TRUE: using cuts from selected tree only
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

    Cuts from selected tree (depth = 2 ):
    [1] "z2 <= 2"    "z1 <= 0"    "size <= 58"
    GRF cuts identified: 3
      Cuts: z2 <= 2, z1 <= 0, size <= 58
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
    Factors per GRF: z2 <= 2 z1 <= 0 size <= 58
    Initial GRF cuts included z2 <= 2 z1 <= 0 size <= 58
    Factors included per GRF (not in lasso) z2 <= 2 z1 <= 0 size <= 58

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
    Seconds and minutes forestsearch overall = 4.749 0.0792
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
      return_selected_cuts_only = TRUE: using cuts from selected tree only
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

    Cuts from selected tree (depth = 2 ):
    [1] "z5 <= 1"    "size <= 18" "z1 <= 0"
    GRF cuts identified: 3
      Cuts: z5 <= 1, size <= 18, z1 <= 0
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
    Factors per GRF: z5 <= 1 size <= 18 z1 <= 0
    Initial GRF cuts included z5 <= 1 size <= 18 z1 <= 0
    Factors included per GRF (not in lasso) z5 <= 1 size <= 18 z1 <= 0

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
    Subgroup search completed in 0.02 minutes

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
    Seconds and minutes forestsearch overall = 5.074 0.0846
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

    Completed in 10.3 minutes

``` r
cat("Results:", nrow(results_alt), "rows\n")
```

    Results: 2000 rows

### 5.4 Running Null Hypothesis Simulations

``` r
cat("Running", sim_config_null$n_sims, "simulations under H0...\n")
```

    Running 1000 simulations under H0...

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
      return_selected_cuts_only = TRUE: using cuts from selected tree only
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
    Seconds and minutes forestsearch overall = 10.937 0.1823
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
      return_selected_cuts_only = TRUE: using cuts from selected tree only
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
      return_selected_cuts_only = TRUE: using cuts from selected tree only
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

    Completed in 8.4 minutes

## 6 Summarizing Results

### 6.1 Operating Characteristics Summary

``` r
t0 <- proc.time()[3]

# Summarize alternative hypothesis results
summary_alt <- summarize_simulation_results(results_alt)
print(summary_alt)
```

                      FS     GRF
    any.H          0.880   0.750
    sens           0.860   0.890
    spec           0.980   0.970
    ppv            0.890   0.840
    npv            0.980   0.980
    Avg(#H)       86.000  98.000
    minH          61.000  60.000
    maxH         226.000 224.000
    Avg(#Hc)     624.000 626.000
    minHc        474.000 476.000
    maxHc        700.000 700.000
    hat(H*)        2.275     NaN
    hat(hat[H])    2.306   2.118
    hat(Hc*)       0.636     NaN
    hat(hat[Hc])   0.641   0.632
    hat(H*)all     2.275     NaN
    hat(Hc*)all    0.636     NaN
    hat(ITT)       0.744   0.744
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
    any.H          0.060   0.070
    sens             NaN     NaN
    spec           0.860   0.890
    ppv            0.000   0.000
    npv            1.000   1.000
    Avg(#H)       99.000  77.000
    minH          61.000  60.000
    maxH         186.000 121.000
    Avg(#Hc)     694.000 695.000
    minHc        514.000 579.000
    maxHc        700.000 700.000
    hat(H*)          NaN     NaN
    hat(hat[H])    1.763   1.507
    hat(Hc*)       0.776     NaN
    hat(hat[Hc])   0.687   0.677
    hat(H*)all       NaN     NaN
    hat(Hc*)all    0.776     NaN
    hat(ITT)       0.700   0.700
    hat(ITTadj)    0.693   0.693

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

    Censoring proportion: 0.453 

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

    Empirical FS: 0.88 

``` r
cat("Empirical FSlg:", round(mean(results_alt[analysis == "FSlg"]$any.H), 3), "\n")
```

    Empirical FSlg: NaN 

``` r
if ("GRF" %in% results_alt$analysis) {
  cat("Empirical GRF:", round(mean(results_alt[analysis == "GRF"]$any.H), 3), "\n")
}
```

    Empirical GRF: 0.754 

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

    Theoretical FS at min(SG) (asymptotic): 0.039618 

``` r
cat("Empirical FS:", round(mean(results_null[analysis == "FS"]$any.H), 6), "\n")
```

    Empirical FS: 0.061 

``` r
cat("Empirical FSlg:", round(mean(results_null[analysis == "FSlg"]$any.H), 6), "\n")
```

    Empirical FSlg: NaN 

``` r
if ("GRF" %in% results_null$analysis) {
  cat("Empirical GRF:", round(mean(results_null[analysis == "GRF"]$any.H), 6), "\n")
}
```

    Empirical GRF: 0.067 

``` r
prop_cens <- mean(results_null$p.cens)  # Censoring proportion
cat("Censoring proportion:", round(prop_cens, 3), "\n")
```

    Censoring proportion: 0.463 

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

![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAABsAAAASACAYAAABBfog7AAAEDmlDQ1BrQ0dDb2xvclNwYWNlR2VuZXJpY1JHQgAAOI2NVV1oHFUUPpu5syskzoPUpqaSDv41lLRsUtGE2uj+ZbNt3CyTbLRBkMns3Z1pJjPj/KRpKT4UQRDBqOCT4P9bwSchaqvtiy2itFCiBIMo+ND6R6HSFwnruTOzu5O4a73L3PnmnO9+595z7t4LkLgsW5beJQIsGq4t5dPis8fmxMQ6dMF90A190C0rjpUqlSYBG+PCv9rt7yDG3tf2t/f/Z+uuUEcBiN2F2Kw4yiLiZQD+FcWyXYAEQfvICddi+AnEO2ycIOISw7UAVxieD/Cyz5mRMohfRSwoqoz+xNuIB+cj9loEB3Pw2448NaitKSLLRck2q5pOI9O9g/t/tkXda8Tbg0+PszB9FN8DuPaXKnKW4YcQn1Xk3HSIry5ps8UQ/2W5aQnxIwBdu7yFcgrxPsRjVXu8HOh0qao30cArp9SZZxDfg3h1wTzKxu5E/LUxX5wKdX5SnAzmDx4A4OIqLbB69yMesE1pKojLjVdoNsfyiPi45hZmAn3uLWdpOtfQOaVmikEs7ovj8hFWpz7EV6mel0L9Xy23FMYlPYZenAx0yDB1/PX6dledmQjikjkXCxqMJS9WtfFCyH9XtSekEF+2dH+P4tzITduTygGfv58a5VCTH5PtXD7EFZiNyUDBhHnsFTBgE0SQIA9pfFtgo6cKGuhooeilaKH41eDs38Ip+f4At1Rq/sjr6NEwQqb/I/DQqsLvaFUjvAx+eWirddAJZnAj1DFJL0mSg/gcIpPkMBkhoyCSJ8lTZIxk0TpKDjXHliJzZPO50dR5ASNSnzeLvIvod0HG/mdkmOC0z8VKnzcQ2M/Yz2vKldduXjp9bleLu0ZWn7vWc+l0JGcaai10yNrUnXLP/8Jf59ewX+c3Wgz+B34Df+vbVrc16zTMVgp9um9bxEfzPU5kPqUtVWxhs6OiWTVW+gIfywB9uXi7CGcGW/zk98k/kmvJ95IfJn/j3uQ+4c5zn3Kfcd+AyF3gLnJfcl9xH3OfR2rUee80a+6vo7EK5mmXUdyfQlrYLTwoZIU9wsPCZEtP6BWGhAlhL3p2N6sTjRdduwbHsG9kq32sgBepc+xurLPW4T9URpYGJ3ym4+8zA05u44QjST8ZIoVtu3qE7fWmdn5LPdqvgcZz8Ww8BWJ8X3w0PhQ/wnCDGd+LvlHs8dRy6bLLDuKMaZ20tZrqisPJ5ONiCq8yKhYM5cCgKOu66Lsc0aYOtZdo5QCwezI4wm9J/v0X23mlZXOfBjj8Jzv3WrY5D+CsA9D7aMs2gGfjve8ArD6mePZSeCfEYt8CONWDw8FXTxrPqx/r9Vt4biXeANh8vV7/+/16ffMD1N8AuKD/A/8leAvFY9bLAAAAOGVYSWZNTQAqAAAACAABh2kABAAAAAEAAAAaAAAAAAACoAIABAAAAAEAAAbAoAMABAAAAAEAAASAAAAAAKPtvSoAAEAASURBVHgB7J0HvCRFubdrWdhAvJIWFBZQQECiSBKRJDkKKChBCVdJKpIUBMkgYUEJEiQv4AVJCn4gGOCSQYLAFck56pJhl12gv/r3UrM11dUzPaHPmenz1O93TndXV3yqurum3nrfGpZYZ3AQgAAEIAABCEAAAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQqAiB6SpSD6oBAQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAACEIAABCAAgZQAAjA6AgQgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAQKUIIACrVHNSGQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAACEIAABCAAAQRg9AEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAAQgAIFKEUAAVqnmpDIQgAAEIAABCEAAAhCAAAQgAAEIQAACEIAABCAAAQhAAAIIwOgDEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAAClSKAAKxSzUllIAABCEAAAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEEIDRByAAAQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAACEIAABCpFAAFYpZqTykAAAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCCAAow9AAAIQgAAEIAABCEAAAhCAAAQgAAEIQAACEIAABCAAAQhUigACsEo1J5WBAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAhBAAEYfgAAEIAABCEAAAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQqBQBBGCVak4qAwEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAAQggACMPgABCEAAAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIFApAgjAKtWcVAYCEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAABGH0AAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCECgUgQQgFWqOakMBCAAAQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAACEIAAAjD6AAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAQKUIIACrVHNSGQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAACEIAABCAAAQRg9AEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAAQgAIFKEUAAVqnmpDIQgAAEIAABCEAAAhCAAAQgAAEIQAACEIAABCAAAQhAAAIIwOgDEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAAClSKAAKxSzUllIAABCEAAAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEEIDRByAAAQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAACEIAABCpFAAFYpZqTykAAAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCCAAow9AAAIQgAAEIAABCEAAAhCAAAQgAAEIQAACEIAABCAAAQhUigACsEo1J5WBAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAhBAAEYfgAAEIAABCEAAAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQqBQBBGCVak4qAwEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAAQggACMPgABCEAAAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIFApAgjAKtWcVAYCEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAABGH0AAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCECgUgQQgFWqOakMBCAAAQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAACEIAAAjD6AAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAQKUIIACrVHNSGQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAACEIAABCAAAQRg9AEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAAQgAIFKEUAAVqnmpDIQgAAEIAABCEAAAhCAAAQgAAEIQAACEIAABCAAAQhAAAIIwOgDEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAAClSKAAKxSzUllIAABCEAAAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEEIDRByAAAQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAACEIAABCpFAAFYpZqTykAAAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCCAAow9AAAIQgAAEIAABCEAAAhCAAAQgAAEIQAACEIAABCAAAQhUigACsEo1J5WBAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAhBAAEYfgAAEIAABCEAAAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQqBQBBGCVak4qAwEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAAQggACMPgABCEAAAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIFApAgjAKtWcVAYCEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAABGH0AAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCECgUgQQgFWqOakMBCAAAQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAACEIAAAjD6AAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAQKUIIACrVHNSGQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAACEIAABCAAAQRg9AEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAAQgAIFKEUAAVqnmpDIQgAAEIAABCEAAAhCAAAQgAAEIQAACEIAABCAAAQhAAAIIwOgDEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAAClSKAAKxSzUllIAABCEAAAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEEIDRByAAAQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAACEIAABCpFAAFYpZqTykAAAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCCAAow9AAAIQgAAEIAABCEAAAhCAAAQgAAEIQAACEIAABCAAAQhUigACsEo1J5WBAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAhBAAEYfgAAEIAABCEAAAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQqBQBBGCVak4qAwEIQAACEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAAQggACMPgABCEAAAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIFApAgjAKtWcVAYCEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAABGH0AAhCAAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIAABCECgUgSmr1RtqAwEIACBQSLwwgsvmPPPP79Q7iNGjDAzzTRT+rfAAguYpZZaysw555yF4hJo4Ag899xz5sILL8xkeMABB5jhw4dn/IeCx8svv2zOOeecplWdbrrpUkaur88+++xm/vnnN0svvbQZPXp00/i9GuD22283f/3rX+uKN++885qdd965zq9XLsoub9FnpGg4cWslbK9w7lY53nvvPfOXv/zFXHfddeahhx4yr732mpkwYYL5zGc+Y5ZYYgmz+OKLp8fVVlvNzDPPPN3KlnRKJHDFFVeY73znO7Uczj33XPPNb36zdl2Vk48//tg88sgj5vHHHzdPPPGE0ZhovvnmM4sssohZeOGF0347bNiwlqqbJIl56qmnzD/+8Q/z4IMPmnfeecd8+tOfNiuuuKJZddVVjb4zOAiUQeDSSy/NfOt/+ctftjR+eeONN8z111+f9mG9y/XO/vznP2+WXXZZ89nPfrZwsTXuOuGEE8z999+fPmNjxowxSy65pNliiy3Sv6IJ6buhNOT+67/+K01rlllmKRp90MM988wz5qKLLsqU46CDDsr44REn0M8MyxwbPvvss2b8+PEZaAceeGDPfGfuu+8+c+aZZ9aV8Re/+IX51Kc+VefXzYsnn3zS/Otf/zI6Pv3002auueYyiy66aPoe03HkyJEtZffiiy8afddbdRpLOPezn/3M/OpXv0ovZ5hhBnPPPfekYwx3nyMEIACBhgTsSwgHAQhAAAIdErATzRrRtf233HLLJRdccEHywQcfdFiS7kWfOHFicsQRRyQ33XRT9xLtkZSK1O3mm2+OtmcvtdFA47Q/NKJMivZ9O+mS/OQnP0leeumlgS56V/KzPzYz9V9++eW7knYZiZRd3qLPSNFwYtBK2JBZkec6jNML11bwlb5rZ5111kz/ij1bdtIyOemkk5IPP/ywF4pPGXIIWIFNYiduam1qF7okkyZNygndv9633HJL8sUvfrFWz1ifXWGFFRKNk4q6G2+8MRk7dmxumlYonJx33nmJFbwVTZJwEChE4NFHH03sQp1M37MCrULx7UR9YoXciZ2czaShZ8MKbpMdd9yx0Djo97//fTLHHHNE01Fa3/rWt5K33367abmuvfbaujQOPfTQpnEGMkCRb/ef//znujq49wzvgOIt1c8MOxkbNiOk37muP/nHKVOmNIs6YPetADtTxueff76U/O++++5k/fXXz+Tns7GL/xIrNCycv13M1TA9P23/3C6krMtD72e7mKaWlsqJgwAEIFCUAEvn7BsWBwEIQGCwCWhVplaJW0FYunp6sMtz5ZVXptoGBx98sHn//fcHuzhdzb/KdesqqBISe/PNN82xxx6broCWhgsOAt0i0K/Ptd790grQu9ZOZBbCIU2YH//4x+ZLX/qSueuuuwrFIdDAEzjssMNSTSiX80477dTyimkXt1ePu+++u5FmiVanN3JapS2trR/+8IeNgqXjjT322MOsu+66qTZoXmCtJLdCBLPZZpsZK1TMC4Y/BFoi8NFHH5kddtjBWIFMS/FcYPcduuyyy4ydPHfedUdpS1rhbfrebzQOUh/fZpttUi3gugS8i9/+9rdG75lmzgq8akGk/bXXXnvVrgf7xDGr4u+NwWZL/tUgcNVVVxm70GRAKrPrrrumWtbSXm3kpJm6/fbbm69+9avGLmpsFDS998ADDzQNUySANM/WXnvtWlCVU5r2OAhAAAJFCCAAK0KJMBCAAAQGiMA///nPdOB52223DVCO9dk8/PDD6cByyy23NDKVUSVX5br1WzvJHNAaa6xh1N9xEOiEQD8/1zJ1qMkDmYtrx2lCQUIFu3q5nejEKZHA//3f/xmZTHNO5vq+//3vu8tKHGX2+fTTTy9cF7s605xyyilRU1MuETH69a9/XdhM0jXXXGM233xzI8EFDgKdEpBJsXYXFfzud78zW2+9tZEp2yJOCxlkxvCtt96KBj/kkEPqBHFf//rX0wVyp512mpl++mm7WOha5uHynJ6Rv//977XbEn7NNttstevBOunnb/dgMSPfoUfg1ltvTYXyA1Hzk08+OWNmsVm+Esxp0YC+742cM7/a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ussvv3yaw6WXXmpuvvnmjnLTirt99tknTWP22Wc3Cy20UEfpERkCEIAABCAAAQhAAAIQgAAEIFCEwLPPTmfOPXcW88EHRUITBgIQgAAEIAABCEAAAoNPAAFYyW1wwAEHpCYLZSZis802M2eeeaZpR8X7gQceMOuuu67RUW7XXXctueQk3yqBYcOGpVHcsdX4hIcABCAQEuB9EhLh2icwbPjIqZfDR/nenEOgjoB7j7hj3U0uIAABCLRA4MQTZzJHHDGHuf766VuIRVAIQAACWQJuXOKO2RD4QAACEChGwL1H3LFYLEINJQKMXEtubZlAPProo81+++1n3nrrrVRwpfPVV1/dLLvssqkW15gxY8zo0aPNqFGj0k2FJSzTfg3PP/+8eeKJJ8z//u//mocffrhWUgnCjjjiiNo1J71BYMYZZzTSzNMRBwEIQKATArPMMovR4G2mmWbqJBniVp3Akv9ths1g+8hi21W9ptSvAwKzzjqrGT58uJl55pk7SIWoEIAABIw1LzTiEwyfLMAACgQgAIE2CXzqU59K505GjuR90iZCokEAAp8QmGOOOdLtgmaYYQaYQCBKAAFYFEt3Pffdd1+jh3GPPfYwEydONO+884659tpr079Wc1p//fXNxRdfbKabDuW9VtmVHV5tokmmstzvf//7VEC65JJLms9//vNNs7niiivSMBK0ur3omkYiQEsEZGf48ccfT//0bKtdFl100Y6FoBKWS9vzxRdfNEsttZT5whe+UOiZbzdeS5Um8IAQKPt9MiCVIJPSCQwbNbsxy+xZej5k0N8EJPwqc3zS33QoPQQg0AoB9xuUFdatUCMsBCAQIzBixAijPxwEIACBTglIkI4wvVOK1Y6PFGWA2nfHHXc0zz77rDnwwAPNPPPM01KueohlPvGaa64x1113Xapl1FICBK4Ege22285stdVW5sorryxUH4XVn/pMK+7jjz82Z5xxhrnxxhtbiTbkwl5++eWpwGuJJZZIn89tttnGLLfccmbeeec1v/rVr9refPO4444zY8eONWussYbZdtttzdJLL50+8+PGjWvIuN14DRPlJgQgAAEIQAACEIAABCAAAQhAAAIQgAAEIACBPiWABtgANtxcc81ljjrqqPTvmWeeMXfeeWeqOSJzh9LckGaY1DVlokYrdaW1o8n1ZZZZBrM1A9hOQz2rddZZx/z1r381l1566VBHkVv/008/3ey+++7pfT2fa665ZmrC9I477jA333yz2Wuvvcx9991nLrjggtw0YjcOPvhgc+SRR6am7zbaaKNU+HX//febG264wUiTdMKECalJ1TBuu/HCdLiGAAQgAAEIQAACEIAABCAAAQhAAAIQgAAEIFAVAgjABqklF1xwQaM/HAR6jcAjjzzSa0XqqfJoX7699947LdP+++9vjjnmmDrzhOecc47ZZZddzIUXXmi+/vWvm80337xQ+WXyUMIvOcWVxp9zl1xyifnOd76T5qU0V1hhBXcrNZXYTrxaApxAAAIQgAAEIAABCEAAAhCAAAQgAAEIQAACEKggAUwgVrBRqRIEIFAeAe3FNmnSpHSvL2l0ur0QXI4777yz2WCDDdLLq666ynk3PcqEodyXv/zlOuGX/L797W+bLbbYQqfmtNNOS4/uX7vxXHyOEIAABCAAAQhAAAIQgAAEIAABCEAAAhCAAASqSAABWBVblToNCgHtnSVTlh9++OGg5N9ppo899lhqau+DDz5Ik3rwwQfT63/+85/p9auvvppeS1NJ7pZbbjGHH364+dOf/mQmT56c+t16661pmBdffDG9Dv+98sor6X2ZCcxzDz/8sBk/fny6X55MDcqsYC8xvffee9Oir7rqqmb66eNKtCuuuGIa5qGHHsqrZp1/kiTm+uuvT/223377unvuQkIwOZmmVD+TazdeGpl/PU3go48+SttZRxwE8ggkE/9jkvt/ZZJJb+QFwR8C6TdU3w2NU3AQgAAEOiHg3iMag+IgAAEIdEJAcwjaDoT3SScUiQsBCIiA5jG1rRAOAnkE4rO3eaHxH1QCU6ZMMRJCODfffPO5U449QOD99983b7zxhlE7zTnnnD1QotaKcNZZZ5lx48bVIrn96r73ve+ZM888M93bauutt05N+q288srmpz/9aS2s9qu69tprzXe/+13z5JNPmrPPPttIEyp0f/vb31Jtpnnnnde89NJLdbfFbtdddzWXXXZZnb8uZPLvoosuSrWuMjdzPA499FBz/vnn59zNeu+xxx5mv/32y94IfGSO8Nxzz60J/YLb6eWzzz6bHueff/7Y7YyfhI+qv9zXvva1zH15rL322qm/tM/+8Y9/mK9+9aum3XjRDPDsKQLvvvtu2ic00fSpT32qp8pGYXqIwMNnm+Tuw82w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0+UaAhCAAAQgAAEIQAACEIAABCAAAQhAAAIiIMHTpEkffvI3xTuX3xQzceLU47Qw064nTpxi5/figqqpgqy8e76/BF4fpcKvoi2y5ZZfsFu2YCawKK+hFA4B2FBqbeoKgQEiIG2umNO+RqGbf/75Uy/t8aX9w7rlpE224oorFk7OlaNoBJV1vfXWM0888UQqsLrsssvM+uuvXzR6XTgnAJOn9n3Snk+he+utt1JTkfKXtppcu/HSyPyDAAQgAAEIQAACEIAABCAAAQhAAAIQ6GkCmmOTwEmCJff3/vvufKr/tOup/v61O/eFVUpn6rU71gu7JHzqN3fFFf9nnn/+TTP//P/Vb0WnvCUTQABWMmCSh8BQIiANLTmZBYy5l156KeO96KKLpn5PP/10qt0lDbOYu//++83ss89u92j7TKF9TLQnl/7KcM8995xZc801jY7aZPOPf/yjWXrppdvOSvbPZd5UG3f+6U9/Mrvttlsmreuvvz71GzlypFl22WXT83bjZRLHAwIQgAAEIAABCEAAAhCAAAQgAAEIQKAlAtJoeu+9yUZCJh39c+fnju6+C+OEWRJuvf/+5DohlxNaKUw/CqNagtilwLPPPtrMOedMXUqNZKpEAAFYlVqTukCgCwTcPiGNzBHmZeOEVw8++GAmiNK76qqrMv7a00waTzL9N27cOHPkkUdmwtx8882pwEkq2GeffXZbZgYzibbpoTJsvfXWqfBLe6GpbN3YE22PPfYwP//5z83FF18cFYBdcsklaYmlZTZixIha6duNV0uAEwhAAAIQgAAEIAABCEAAAhCAAAQgUFECElK9++4Hdt5pcnp8910dp56/884H6bkTSsWPEmxNFXA5YZYEVvL7+OOkotR6o1p2txQ7BzbcjBw5/Sd/8fMxY2a2i+BXs1unzNAbBacUPUUAAVhPNQeF6WcC2odJe0LNNttsZsYZZ+zbqsw888xp2e+++26z3XbbtbSB5CqrrGJuueUWc95555nvfOc75itf+Uqa1vvvv2923HFH88ADD2S4SKtLgp/99tvPHH/88al5v+23374W7pVXXjE77bRTan947rnnNttuu23t3mCcnHPOOebOO+9MuRx11FF2lc775pFHHokWZYYZZjALL7xw7d7LL79sjj322PT6W9/6lllppZVq93bffXdzzDHHmNtuu80cccQR5uCD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70Hr38XmIRE7eTSKTp7jq/gb5Oqrr+7svPPO1eeO+Axz4YUX1ouo1n0+qSisjCEQQyFIBKa1wJ/+9KfOi1/84urBGgGsMuhVBsZimXtt9XUa5AvrOPjQQw/t5LH2R9QVXybV64n1vffeu3jo963QDgIEFiuBJr+wLk/8jDPO6OSJ60c8P7qfXfFl9t13310eYkmAwBQQaPp5Ev+0k+cBG/EsiX/YmTFjxohtT3nKUzp5iOYpIOgUCEwfgfgnmO6/MwZ9feutt46AGuQLpjjA55MRbF4QmDICk/088flkytw6ToRAT4H4nmKjjTYa8TklAl3df4M87nGP65xwwgk9y4iNPp/0pbGjJmAIxPwXgDS9BZZbbrmiS+373ve+FOuRymHDYn3mzJkp995K++yzT7xcqHTAAQekK6+8Mm211VZVOeWk80suuWTKD/+U/7MqHX300SleSwQIEOgnsOOOO6ZLL720mN+rnMS+fHblD43FMyv2515n/YqwnQABAmmVVVZJMadHzPG17LLLFiL5b4WUe5QW6yuuuGIx3EgMfRbzCUoECCw+Atdcc82kN9bnk0knVyGBSRGY7OeJzyeTcllVQqA1gfie4oorriimiCn/Bsk9w6q/QfI/9xZDJcaz5y1vectCt9Pnk4UmXKwLeEQEwxbrM9B4Ag0L3HTTTSkP+ZNiYsZ11lknrbnmmin/F0LDtaRijNvrrrsu3XDDDWn11VdPMW5+1CkRIEBgogL3339/ivkFY/zrtdZaK6277rrVF9kTLUt+AgSmr0D80Rlz+8Tnk5hLcJNNNik+o0xfEWdOgMDCCPh8sjB6jiVAoBTw+aSUsCQwNQXiPZ6niEnXXntt8TfIBhtsUHwfG0GwYSSfT4ahumiXKQC2aF8frSNAgAABAgQIECBAgAABAgQIECBAgAABAgQIEJigQPPdWibYANkJECBAgAABAgQIECBAgAABAgQIECBAgAABAgQINCkgANakprIIECBAgAABAgQIECBAgAABAgQIECBAgAABAgRaFxAAa/0SaAABAgQIECBAgAABAgQIECBAgAABAgQIECBAgECTAgJgTWoqiwABAgQIECBAgAABAgQIECBAgAABAgQIECBAoHUBAbDWL4EGECBAgAABAgQIECBAgAABAgQIECBAgAABAgQINCkgANakprIIECBAgAABAgQIECBAgAABAgQIECBAgAABAgRaFxAAa/0SaAABAgQIECBAgAABAgQIECBAgAABAgQIECBAgECTAgJgTWoqiwABAgQIECBAgAABAgQIECBAgAABAgQIECBAoHUBAbDWL4EGECBAgAABAgQIECBAgAABAgQIECBAgAABAgQINCkgANakprIIECBAgAABAgQIECBAgAABAgQIECBAgAABAgRaFxAAa/0SaAABAgQIECBAgAABAgQIECBAgAABAgQIECBAgECTAgJgTWoqiwABAgQIECBAgAABAgQIECBAgAABAgQIECBAoHUBAbDWL4EGECBAgAABAgQIECBAgAABAgQIECBAgAABAgQINCkgANakprIIECBAgAABAgQIECBAgAABAgQIECBAgAABAgRaFxAAa/0SaAABAgQIECBAgAABAgQIECBAgAABAgQIECBAgECTAgJgTWoqiwABAgQIECBAgAABAgQIECBAgAABAgQIECBAoHUBAbDWL4EGECBAgAABAgQIECBAgAABAgQIECBAgAABAgQINCkgANakprIIECBAgAABAgQIECBAgAABAgQIECBAgAABAgRaFxAAa/0SaAABAgQIECBAgAABAgQIECBAgAABAgQIECBAgECTAgJgTWoqiwABAgQIECBAgAABAgQIECBAgAABAgQIECBAoHUBAbDWL4EGECBAgAABAgQIECBAgAABAgQIECBAgAABAgQINCkgANakprIIECBAgAABAgQIECBAgAABAgQIECBAgAABAgRaFxAAa/0SaAABAgQIECBAgAABAgQIECBAgAABAgQIECBAgECTAgJgTWoqiwABAgQIECBAgAABAgQIECBAgAABAgQIECBAoHUBAbDWL4EGECBAgAABAgQIECBAgAABAgQIECBAgAABAgQINCkgANakprIIECBAgAABAgQIECBAgAABAgQIECBAgAABAgRaFxAAa/0SaAABAgQIECBAgAABAgQIECBAgAABAgQIECBAgECTAgJgTWoqiwABAgQIECBAgAABAgQIECBAgAABAgQIECBAoHUBAbDWL4EGECBAgAABAgQIECBAgAABAgQIECBAgAABAgQINCkgANakprIIECBAgAABAgQIECBAgAABAgQIECBAgAABAgRaFxAAa/0SaAABAgQIECBAgAABAgQIECBAgAABAgQIECBAgECTAgJgTWoqiwABAgQIECBAgAABAgQIECBAgAABAgQIECBAoHUBAbDWL4EGECBAgAABAgQIECBAgAABAgQIECBAgAABAgQINCkgANakprIIECBAgAABAgQIECBAgAABAgQIECBAgAABAgRaFxAAa/0SaAABAgQIECBAgAABAgQIECBAgAABAgQIECBAgECTAgJgTWoqiwABAgQIECBAgAABAgQIECBAgAABAgQIECBAoHUBAbDWL4EGECBAgAABAgQIECBAgAABAgQIECBAgAABAgQINCkgANakprIIECBAgAABAgQIECBAgAABAgQIECBAgAABAgRaFxAAa/0SaAABAgQIECBAgAABAgQIECBAgAABAgQIECBAgECTAgJgTWoqiwABAgQIECBAgAABAgQIECBAgAABAgQIECBAoHUBAbDWL4EGECBAgAABAgQIECBAgAABAgQIECBAgAABAgQINCkgANakprIIECBAgAABAgQIECBAgAABAgQIECBAgAABAgRaFxAAa/0SaAABAgQIECBAgAABAgQIECBAgAABAgQIECBAgECTAgJgTWoqiwABAgQIECBAgAABAgQIECBAgAABAgQIECBAoHUBAbDWL4EGECBAgAABAgQIECBAgAABAgQIECBAgAABAgQINCkgANakprIIECBAgAABAgQIECBAgAABAgQIECBAgAABAgRaFxAAa/0SaAABAgQIECBAgAABAgQIECBAgAABAgQIECBAgECTAgJgTWoqiwABAgQIECBAgAABAgQIECBAgAABAgQIECBAoHUBAbDWL4EGECBAgAABAgQIECBAgAABAgQIECBAgAABAgQINCkgANakprIIECBAgAABAgQIECBAgAABAgQIECBAgAABAgRaFxAAa/0SaAABAgQIECBAgAABAgQIECBAgAABAgQIECBAgECTAks2WZiyCBAgQIAAAQIECBAgQIAAgYkJXH311emWW25JM2fOTKuvvvrEDp6muc8555z08MMPp5e85CVpySV9tTFNbwOnTYAAAQIECBAYU0APsDF57CRAgAABAgQIECAwWuA1r3lNWnnllRf655RTThld+BTdsu6661Zed99994TOsm794IMPTuhYmVM699xzK/sjjzxygUhWXXXVqoz69eheX2211dKGG26Ytthii/T2t789zZs3L032Nbv55pt7nuMZZ5xRncPxxx/fM08bG//whz+k5z3veWnOnDlpxowZPZvQ75xOO+206pxOPPHEnse2sXH+/PnptttuG2rVZ511Vtpmm23SEUccMdR6FE6AAAECBAgQILD4Cvg3qcX32mk5AQIECBAgQIBASwL33HNPuvPOOxe69gceeGChy1hcCrjrrrvSvffeWzQ3em1MJDVhPZH6plreCECVhn//+98X6PTi+H/9618TOvbHP/5xikDTwQcfnCJQs9lmm03o+IlmjnP70Ic+lCKwfPvtt486/B//+EflcP/994/a39aGPfbYo3hvRNuf9KQnjWjGeOcUz5Dy2i4q53TJJZekPffcM+26665pr732GnE+Tb6YO3duiqDmIYcckrbddtsUQXaJAAECBAgQIECAQF1AAKyuYZ0AAQIECBAgQIDAAALxRevf/va3njnjC+trr7222Lf00ksXvWF6ZswbV1xxxX67bCewyAq8/vWvT0sttVTP9kWQ6dZbb03XX399Knv63XDDDUVPnZ/97GdpjTXW6HlcExtnzZqVrrzyyrT88ss3UdyklHHmmWemr3/962mVVVZJ++2336g6F7dz+uY3v5lmz5496jyGsSGChf/93/9dBFh32WWXdPHFF6dHPOIRw6hKmQQIECBAgAABAoupgADYYnrhNJsAAQIECBAgQKA9geOOO65v5VdccUXadNNNi/1Pe9rT0mWXXdY3rx0EFkeBz372s+kxj3nMmE3/5z//mWJIvgjqRLA4gmEHHnhgOvXUU8c8bmF2/vnPfy4O7xcEef7zn59i2LxIz3jGM4plm7/uu+++qofUAQcckJZZZplRzRnvnEYd0PKG6B1bpn7XodzfxHLfffdNRx11VIpeZ3Ff7r777k0UqwwCBAgQIECAAIEpImAOsClyIZ0GAQIECBAgQIAAAQIEFhWB6CEWc4B9+tOfrpr05S9/uRoGs9o4iSurr7562m677YqfCE63nSJgE0ODLrfccsVwgW23Z3Gs/3GPe1wV9PrIRz4y4WE6F8dz1mYCBAgQIECAAIHBBQTABreSkwABAgQIECBAgAABAgQmIBDDJZa9xR566KF03XXXTeDoqZs1esgdffTRxQm++tWvTo961KOm7skO+cx22mmnooYYejOCrBIBAgQIECBAgACBUsAQiKWEJQECBAgQIECAAIEWBS688MJ0+umnFy346Ec/mmJ4tFh++9vfLgIIL37xi4vh5FZbbbW0zz77pJhraaWVVkof/OAH+7b6F7/4RfrMZz5T7N9mm22KeZh6ZY4v40844YRiuMarrroq/e53v0sxz9kzn/nMtPXWW6f4gn5xSNdcc0362te+ln7zm98UPzfeeGNhtPbaa6f4iS/Kn/WsZ406lZiv6hOf+MSo7eNteMc73tFzKL0FbUfU9573vKcYMnDLLbdMb3jDG4phBE855ZQUZW6wwQZFb5cddthhxFxHDz74YIo8cf/83//9X3rkIx+Z4vgY8m/nnXce7zSGun+JJZYo5rqLufEi1YfI6654Qd3e/e53p/vvvz/96U9/KoqMIRfLofDWWWedFPsj/eQnP0knn3xysb7jjjumF7zgBcV6968oK3qu/fSnPy3m87vzzjvT+uuvnzbaaKO01VZbFT3Iuo+Z6OsvfelL6fe//31xWLSlOw16Tt3HxesLLrggffWrXy3ezzfddFPaZJNN0nOe85z0kpe8JG222Wa9DhmxbUGeBz/4wQ+K4FO8l8oU9+Mvf/nL4mWcT1yLepo/f36aN29eimdO+Z4N++idF3nj+RPXcbzgYFyXyB91H3HEEel1r3tdvRrrBAgQIECAAAEC01mgIxEgQIAAAQIECBAg0JjA5Zdf3sl/XxQ/G2644cDlHnvssdVxOXDVycO1Va/L8r7whS8U5T32sY8t9uUviscs/+tf/3pVxsEHH9wzbw46dDbeeOMqX1lXffna1762k+ci6nn8oBsf//jHV3Xccccdgx5W5Ku3JQf+Rh37r3/9q5ODhZ38RXlVR/2Ycj0Hhjr5i/hO/oJ/RBk5+DjmceXx3cvwraeFbUeUtcIKKxRt2WWXXTof//jHR7Ur7ouHH364qjb3qOqsueaao/KVbc2Bz04OMlT7586dWx07kZUll1yyKiMHmAY+9Le//W11XLQpB2RGHbuwbvV7qzzvcpmDgFV9J510UtWWI488stpeXznjjDM6T3nKU6p8ZTn15Stf+crOH/7wh/phE14v33N5+MNO7hk36vhBz+nzn/981dY8F1bn/e9/fyfPvVVtq7c77v8zzzxzVF31DQv6PAjPel3d6/Eeq6eoJwejxzwmysiB+E4ORNYP7bn+vve9ryrroosu6pnHRgIECBAgQIAAgeknoAdY/lQtESBAgAABAgQIEFiUBKIXUPTCihS9eaKnxKMf/egUPX+aTPmL4vSyl70sPfDAA0Wxm266adFLZI011ki/+tWv0rnnnlsMWZeDAunSSy9N0aNs2WWXbbIJjZSVA3TpK1/5SlFWzKcUw+5FL5Kll166cIxedNH2cIyeXtEjJvKUKQedxuzlFj1QoqwwueKKK4rDYu6hHOAsiyiWC9uOemExVOCpp55abYqeVDnwld74xjdWvb/iHtliiy3SH//4xyJfnPOrXvWqtNZaa6Wf//znRW+4s88+O1188cVVOZO5Et6HHXZYVeUTn/jElAN41etyZWHd3vKWtxT38Be/+MX017/+tbhWZc+37l5HZZ29lt/61rdStKVMm2++eXre856Xotdl3D9xH0WvrcgX1/7KK68s9pX5B13m4FlxfSL/7NmzU8yX1p0W5JwOOuigYo61HAAreiZG780cXEzf//73069//evi/o/eUTmAnl760pd2V5kW5nkwc+bMordW3LfRAy1S9EIs3yNPfvKTq/ryPwkU+6LnYqQcpCx6msb9ET0Eoydg9JCL/VFe3PPXXnttivdAv7T99tunmAMsUlynWbNm9ctqOwECBAgQIECAwHQSmH4xP2dMgAABAgQIECBAYHgCTfQAy3+PdJ761Kd2cvCiEz2eogdE9EwpUxM9wKIXVB5Sr+o18bGPfayTAxZlFcUyB8Y6OZBQ5dl3331H7J/Ii3qPlhzg6+Qv6wf+CY/yp7sHWB46rerxst566/XsqRY9pj7wgQ9UZfznf/7nRJpe5M0BhE55DtEb6rzzzhtRRlPtKHuAlecbPdai11QOcnVyQKxz8803V/XmgGh1Tttuu22nu2dWHJeHh6vyRJkf+tCHquMnsjLRHmDRwycPFzii7hygGlVlU25RcA7cFvWFYa80Vg+w8M0BmOL4HJDq5EDpqCKiF2T0/iqvTR4adFSeQTbUe+Qdc8wxYx4y3jnVe4BFu+IezUMKjigz3us5oFa1O65Ld2rqeVA3jl6tvdKcOXOqtuShJntl6eTgbudJT3pSlS8HHXvmKzfGs6u8R3Mgv9xsSYAAAQIECBAgMM0F0jQ/f6dPgAABAgQIECBAoFGBJgJguadD55JLLunbriYCYPHFe/lF/lvf+ta+dcWOZz/72UXeCAzknhhj5u23swwelXUu6LI7ALbHHntU55F7rPWrvhj2sHSLL8rrwwj2PejfO2K4uxh6sGxznldt1CFNtaMeAMu9dEbVU26I61C2J88F18lzJ5W7RizznGCdGPquzNtEACyG74sgYvdP7lnXyT2uqkBhWWcsd9tttxHtKl805RbljRcsqgdnuodAzD0CK6P99tuvbN6oZQxXGMHp8twiSD3RFO+38vjuYTS7yxrvnOoBsBkzZnT6Df+X50fr5F6kRb2Rr/v+b+p5UDfuFQDLPeg68RyJ88897LpPd8TrGNa0dBrkvi3fo/H8zD3JRpTlBQECBAgQIECAwPQU6D+GQP6kKREgQIAAAQIECBAgMPkCObBQDG03zJo/9alPVcXneaGq9V4ruRdSsTn3EkmnnXZaryytbdtzzz3TWWedlfJ8WelFL3pR33bkoFcxNGBkiGHhciCjb976jr///e8pz6NVDUkZFrvvvns9S7E+jHbsvffeo+opN8Qwb2XaZ5990jLLLFO+HLGMYRHrw/qN2LmAL2J4xRgKsvsnhgS8/vrri2H4yqLz/GTpG9/4RjruuOPKTSOWw3AbUcGAL84555wi53/8x3+k9773vX2PiuEKczCm2h/DLk40xXCKZcrzjZWrC73MPbv6Dv0XQ4OWw0HGvR9DRdbTZD0PYvjF3HsyHX/88enwww+vN2HUetne2BHvw/FSaRlDhcZwpRIBAgQIECBAgAABc4C5BwgQIECAAAECBAgsYgIRtBhmirmZbrzxxqKKCFDkocbGrG6zzTar9uceRdX6gq7E3GIx79agKeYm65ee/vSnp/gZK8WcSzEvUTlXVuSNIFjM7TVWCqeddtqpmvcr5mvKQ0X2PGQY7RjrPrjwwgurdsQ8VWOlCOA1Gbhcf/31i7npos4INvzlL39Jd955Z4oAaaSYg+rAAw8s5lmrBzGKnV2/huHWVcW4L+++++6Uhzcs8sW8X8svv/yYx4RnmfIQjuXqwMuor0wxv1hTKQLnY6VVV1015eERiyz33XdfinnsIk3m8yDqjEBd/PRL8d687LLLUhmUjHyxbbxUt6wbj3ec/QQIECBAgAABAlNXQABs6l5bZ0aAAAECBAgQILCYCqy99tpDbXmeX6fqAXX77ben+GJ80NREAOyZz3xmyvMtDVrlwPnyMG/pRz/6UcrDA6ZoZ/RGuu6669Jdd901qow8AMiobd0bogdWHuKu2JyH90vR2ycPr9adbdTrJtqRhy1MEZzslyKoV6bxrt94+8tyBl1GcOIxj3nMiOwRXMzzZqVDDz00htlPRx99dMpD0lW9jkZk7vOiCbc+RY+5Oe6XMuUhB8vVvsvoJRZBsmhv3GNxvhH0GzSVwZnotTeRQPB45Ucvr7FS9IIsUwQuy9Tm8yB6EUaPwvL9Gp6//e1vq+dT2cZB3q9lD7A4pjQuj7ckQIAAAQIECBCYngL//xPw9Dx/Z02AAAECBAgQIEBgkRNoOmDRfYLxBXOZHnzwwZTn5Slfjru84YYbxs0z2RluueWW9OEPfzidcsopKc8R1rP66B0SvV7iZ5AUPb3yXF9F1rge3/rWt0YFfbrLabIdESCM4fb6pfoX/OP14KsHBvqVt7DbI5DzkY98JEVdMazhvffem972trcVQzPOmTNnzOKbdBuzoj47I+hSpkGt4p6IAFgMzRfLQQNZDzzwQIqfSE9+8pPLahtZLr300gtUThvPg3nz5qXDDjus71CFce+vtdZaRQB70JOqPzfjmkgECBAgQIAAAQIEBMDcAwQIECBAgAABAgQWMYFBehlFk8frFdFv2LB6750Y3vCNb3zjwAJjBWUGLqTBjNGDbdasWdUcXVH0E57whDRz5sz0jGc8I2200UYpzjGWMZRi9DiJNFaPnTPPPDPtv//+Rb6YsyiCX+MFK5pux3j3QPRCKlME9cYKwPS7D8rjm1y+4x3vSLfddlsRDIty3/rWt6YYyvFZz3pWz2qadutZyTgb6z2nInA3SKoHWMqhBAc5LoJU8R6K4SIj+LwopMl+HkSw+gMf+EB16tHbMe6ReI/GezZ+4j198cUXp1e+8pVFvrHer2VBdc/6+6Pcb0mAAAECBAgQIDD9BATApt81d8YECBAgQIAAAQJTRKCcc6nf6fTr7dQ9L9Nee+3Vr4hFfvtrX/vaKvgVcyB97nOfSzFcYa8Uc1WVqT4EXLktljGE4pve9KYiuBhBqC996UsphmwcLzXdjvHqW3nlldOvf/3rIlv0oBorABZD3E1mOuSQQ9L555+ffvrTnxY98l796lenyy+/vOdcc5Pt1suhPuTozTff3CvLiG3xvovAXaRll102zZgxY8T+sV5EICd690WQsD4n3VjHDHvfZD4P4r4og18xJGMMmbnrrrsWjt3nOcj7tX5M3XO8XpH146wTIECAAAECBAhMXYHxB7CfuufuzAgQIECAAAECBAgslgJlL6xyKLV+J3HTTTdVu+q9xVZcccX0+Mc/vtj3y1/+suiNUmXssXL//fen8847rxiObLw6exw+tE1//vOfi4BVVBC9cKLHSL/gV/SCqs+b1atXVMxDNHv27GoYxZjTapttthm3/U23Y9wKc4boJVOmX/3qV+Vqz2V9iLueGRreGIGNE088MT3qUY8qSo5gUfQM605tuHW3IV5H76Oyh1FY1t8rvfJHnjKAGsdONJXBmRiuM4ZQbDtN5vMgelOW6eCDD07vec97ega/Ik8EdsvU6/1a7iuXAmClhCUBAgQIECBAgEApIABWSlgSIECAAAECBAgQWEwEyiHX4gvf+lBs9ebHF8Znn312fdOI9c0337x4HQGtY445ZsS+7hef/exn00tf+tK03nrrpZ133rl7d2uvL7vssipY8exnPzs9+tGP7tuWb37zm+lvf/tbtX/+/PnVeqyE5cte9rLKM+axete73jUiT78XTbajXx3d27fbbrtq0xFHHFE5VBv/vRLBnGOPPbZ789Bfr7/++ikCHGX6xje+kb7yla+UL4tl024xlF6k7mtbbBzj1zLLLFMFFGM+sO52dh8ac52Vadttty1XB16WAbA4oB6U7VXAgp5Tr7LG2tbU86Bsb9TV6zpEr8AyvfCFLyxXRy0feuihdNZZZ1Xbe5VV7fz3St1yvCFLu4/1mgABAgQIECBAYGoKCIBNzevqrAgQIECAAAECBKawQDkkXwQ39ttvv55nussuu6Srrrqq577YeNhhh6Xyy+qDDjqoGk6v+4DrrrsuxZw9ZXrnO99Zrra+XGWVVao2/OY3v0nxpXmvFMMads9zFr1vyhRBwJhr6IYbbig2RVDjqKOOKnePu2yqHeNWVMsQAb+NN9642BK9+PoFMT/zmc+k8XqI1YptdDV698S8TmXae++9U32Orabdyrms/vrXv1a9+Mq6x1t+/OMfr7K8973vTb///e+r1/WV6MFUBsiih9ucOXPquwda32KLLap89YBQtbG2sjDnVCtm3NWmngdle6PCu+66a1S99Wt+9dVXj9ofGyLYFa7lfH2xrf5+jde9UmkZwa/VV1+9VxbbCBAgQIAAAQIEppmAANg0u+BOlwABAgQIECBAYPEX2G233aqTOPnkk4vh5WKIwhjyMIaei2BObB+rF0QEJvbYY4+inOgZtemmm6boSRTDAEZgLQIVUdZWW21V9YrafvvtU/3L+6oRLa1EL6OYCytSzHO14447FkM1Pvjgg8UQdVdeeWWKYQxf/vKXpxjGsZ7qw6W97W1vS5deemmxe6WVVko77LBDOv3009Pxxx+fPvWpTxU9qKIXVffPGWecURzTVDvq7RtvPYbsi/aVQczorRZzucW8YDE8XwQE99133xQ92crh/cYrs+n9MRTi5z//+RRzqUWKHjr/8z//U1XTtFvMrRUpej/G0JWf/OQn0wknnFDVN9bKC17wghT3d6QYMjKCzKeeemoxV1cEZCLIuP/++xdDZJZDJB5++OFprbXWGqvYnvvi/VmmSy65pFztuVyYc+pZYJ+NTT0PyvZGNRFEjvm+IrgYz5VI8TwpU3geffTRxXMrtsV78mtf+1qKOeO++tWvltmKZf39OmLHv19EoL7ME75t3fO92mYbAQIECBAgQIBAiwL5w7tEgAABAgQIECBAgEBDApdffnknf7wvfjbccMOBS83Bleq4HDQY87gc4OjkL3mr/GV99eU666zTyT2aqjx5OLpRZebAV+ftb397lac8fsaMGaO2bb311p3cC2NUGYNuyHOOVWXecccdgx5W5CvbFcvuNnz3u9/t5C+7q7Ijz2Mf+9jOsssuW23LgZhODrx0coCv2paHdazasMYaa1Tb63WNt56DJFUZTbQjClthhRWKtqy66qpV2WOt5GHiOksvvfSI9sf5lm3Pw/t1Pv3pT1ev586dO1ZxfffVy4z7ZiIpB+eq+uNaXXjhhdXhTblFgfX3UHn+eX6rqq6TTjqpaseRRx5ZbS9XcgClk4eWrPKUZdTPPbblgF7nkEMOKQ9boOVTn/rUop7xnhHjnVM8K8p2Rt6xUg4sVXlvvvnmUVmbeB7kQHNntdVWq+op25YDyUV9OZjYmTVr1qj9OVg/4n28wQYbdH74wx92ckCtyJt7jnXi2H6p7pCHfu2XzXYCBAgQIECAAIFpJqAHWP5ELhEgQIAAAQIECBBYnASid0PMqfSxj30s5YDJiKbH/GCzZ89O0bNkvGHAYriymN8rByHSzJkzq95E9aEE8xf1RS+ac845J8WQb4taetGLXpQuuOCCtMkmm1RNix5t9913X8qBvKLHSfQEiyHeomdJmebNm1euNrJsqx0xF1hc6+jlVqboARVp7bXXTnHd6vvKPJO5zEG3lIMiRZX57+2Ug67VkHZNukVvt/e///2p3gvp7rvvrnowjnfOyy+/fDHv1GmnnZbivi97EZWe0dsuehdddNFF6X//93/HK27M/WUvsOixd9ttt/XNu7Dn1LfgHjuaeB7EfGrf+973ih6lSy21VFXLtddeW6xHb8CYm/CAAw5IOXBb7b/99tuLnqfRo/PQQw9NP//5z9OWW26Z4v6IFENS5sBpsd7rV/SAjRTnMNbcYr2OtY0AAQIECBAgQGDqCjwiAn5T9/ScGQECBAgQIECAAIGpLxDD/8V8XxHwiqHMyiHnJnrmMc9ODJ0Xw5XlnjPF8G65J9IClzfR+hcmf/xZc8sttxTzeN1zzz1pvfXWS+uuu27KvXcWptgJH9tmOyKQEtcv5jT7r//6rxTDOS4uqWm3sIi5wCLwVp+XaiIecXzMn3brrbemmLsq96os3hcTKaNf3hiy7+lPf3oxXOUHP/jBdOCBB/bLWm1v4pyqwgZYWdjnQQxFev3116cILEZgq/u5FIHqGG4y5t5bbrnliiB8d0B/gGYWQx/p2mdDAAAHiklEQVTG9YnAfcwzF8MqSgQIECBAgAABAgRCQADMfUCAAAECBAgQIECAAAECBCZZYKeddkoxj1wegjPdeOONVY+zSW7GYl9dzPX27ne/u+jxGcG0CNpLBAgQIECAAAECBELAEIjuAwIECBAgQIAAAQIECBAgMMkCMVxjDLOY5+NK559//iTXPnWq+8IXvlCczM477yz4NXUuqzMhQIAAAQIECDQiIADWCKNCCBAgQIAAAQIECBAgQIDA4AIx7972229fHPDRj3508APlrARiPrGYRy3mE4t5xSQCBAgQIECAAAECdQEBsLqGdQIECBAgQIAAAQIECBAgMEkCMV9VzH8VPcDOPffcSap1alQzf/78tP/++xcnM3fu3LTmmmtOjRNzFgQIECBAgAABAo0JCIA1RqkgAgQIECBAgAABAgQIECAwuMDKK6+cjj322OKACOY8/PDDgx88zXOedNJJ6ZprrknPec5z0r777jvNNZw+AQIECBAgQIBALwEBsF4qthEgQIAAAQIECBAgQIAAgUkQmDNnTtpnn33SSiutlH72s59NQo1To4qLL744bb311ikCYUss4auNqXFVnQUBAgQIECBAoFmBR3RyarZIpREgQIAAAQIECBAgQIAAAQIECBAgQIAAAQIECBBoT8C/SbVnr2YCBAgQIECAAAECBAgQIECAAAECBAgQIECAAIEhCAiADQFVkQQIECBAgAABAgQIECBAgAABAgQIECBAgAABAu0JCIC1Z69mAgQIECBAgAABAgQIECBAgAABAgQIECBAgACBIQgIgA0BVZEECBAgQIAAAQIECBAgQIAAAQIECBAgQIAAAQLtCQiAtWevZgIECBAgQIAAAQIECBAgQIAAAQIECBAgQIAAgSEICIANAVWRBAgQIECAAAECBAgQIECAAAECBAgQIECAAAEC7QkIgLVnr2YCBAgQIECAAAECBAgQIECAAAECBAgQIECAAIEhCAiADQFVkQQIECBAgAABAgQIECBAgAABAgQIECBAgAABAu0JCIC1Z69mAgQIECBAgAABAgQIECBAgAABAgQIECBAgACBIQgIgA0BVZEECBAgQIAAAQIECBAgQIAAAQIECBAgQIAAAQLtCQiAtWevZgIECBAgQIAAAQIECBAgQIAAAQIECBAgQIAAgSEICIANAVWRBAgQIECAAAECBAgQIECAAAECBAgQIECAAAEC7QkIgLVnr2YCBAgQIECAAAECBAgQIECAAAECBAgQIECAAIEhCAiADQFVkQQIECBAgAABAgQIECBAgAABAgQIECBAgAABAu0JCIC1Z69mAgQIECBAgAABAgQIECBAgAABAgQIECBAgACBIQgIgA0BVZEECBAgQIAAAQIECBAgQIAAAQIECBAgQIAAAQLtCQiAtWevZgIECBAgQIAAAQIECBAgQIAAAQIECBAgQIAAgSEICIANAVWRBAgQIECAAAECBAgQIECAAAECBAgQIECAAAEC7QkIgLVnr2YCBAgQIECAAAECBAgQIECAAAECBAgQIECAAIEhCAiADQFVkQQIECBAgAABAgQIECBAgAABAgQIECBAgAABAu0JCIC1Z69mAgQIECBAgAABAgQIECBAgAABAgQIECBAgACBIQgIgA0BVZEECBAgQIAAAQIECBAgQIAAAQIECBAgQIAAAQLtCQiAtWevZgIECBAgQIAAAQIECBAgQIAAAQIECBAgQIAAgSEICIANAVWRBAgQIECAAAECBAgQIECAAAECBAgQIECAAAEC7QkIgLVnr2YCBAgQIECAAAECBAgQIECAAAECBAgQIECAAIEhCAiADQFVkQQIECBAgAABAgQIECBAgAABAgQIECBAgAABAu0JCIC1Z69mAgQIECBAgAABAgQIECBAgAABAgQIECBAgACBIQgIgA0BVZEECBAgQIAAAQIECBAgQIAAAQIECBAgQIAAAQLtCQiAtWevZgIECBAgQIAAAQIECBAgQIAAAQIECBAgQIAAgSEICIANAVWRBAgQIECAAAECBAgQIECAAAECBAgQIECAAAEC7QkIgLVnr2YCBAgQIECAAAECBAgQIECAAAECBAgQIECAAIEhCAiADQFVkQQIECBAgAABAgQIECBAgAABAgQIECBAgAABAu0JCIC1Z69mAgQIECBAgAABAgQIECBAgAABAgQIECBAgACBIQgIgA0BVZEECBAgQIAAAQIECBAgQIAAAQIECBAgQIAAAQLtCQiAtWevZgIECBAgQIAAAQIECBAgQIAAAQIECBAgQIAAgSEICIANAVWRBAgQIECAAAECBAgQIECAAAECBAgQIECAAAEC7QkIgLVnr2YCBAgQIECAAAECBAgQIECAAAECBAgQIECAAIEhCAiADQFVkQQIECBAgAABAgQIECBAgAABAgQIECBAgAABAu0JCIC1Z69mAgQIECBAgAABAgQIECBAgAABAgQIECBAgACBIQgIgA0BVZEECBAgQIAAAQIECBAgQIAAAQIECBAgQIAAAQLtCQiAtWevZgIECBAgQIAAAQIECBAgQIAAAQIECBAgQIAAgSEI/D9mTSTrJ7NjegAAAABJRU5ErkJggg==)

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

      FS: Power = 0.817, Sens = 0.875, Spec = 0.978, PPV = 0.870
      GRF: Power = 0.817, Sens = 0.875, Spec = 0.978, PPV = 0.870

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

      FS: Type I Error = 0.0640
      GRF: Type I Error = 0.0640

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
| 1000 H1 + 1000 H0 simulations, 13 workers |  |  |  |
| Stage | Time (sec)¹ | Time (min) | % of Total |
| DGM creation (H1) | 0.0 | 0.00 | 0.0 |
| Calibrate k_inter (Cox) | 2.2 | 0.04 | 0.2 |
| Calibrate k_inter (AHR) | 0.9 | 0.01 | 0.1 |
| Validate k_inter | 0.3 | 0.00 | 0.0 |
| DGM creation (H0) | 0.0 | 0.00 | 0.0 |
| Simulations H1 | 618.4 | 10.31 | 54.3 |
| Simulations H0 | 502.5 | 8.38 | 44.1 |
| Summarize H1 | 0.0 | 0.00 | 0.0 |
| Summarize H0 | 0.0 | 0.00 | 0.0 |
| Total vignette | 1,139.2 | 18.99 | 100.0 |
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

      H1: 0.6 sec/sim (wall) across 1000 sims on 13 workers

``` r
cat(sprintf("  H0: %.1f sec/sim (wall) across %d sims on %d workers\n",
            timings$sims_null_elapsed / sim_config_null$n_sims,
            sim_config_null$n_sims, n_workers))
```

      H0: 0.5 sec/sim (wall) across 1000 sims on 13 workers

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
| Across 1,000 simulations per scenario |  |  |
|  | FS | GRF |
| M Null: N=700, theta(ITT) = 0.72 |  |  |
| any(H) | 0.06 | 0.07 |
| sens(Hc) | 0.86 | 0.89 |
| ppv(Hc) | 1.00 | 1.00 |
| avg\|H\| | 99.00 | 77.00 |
| M Alt: N=700, p_H=13%, theta(H)=2, theta(Hc)=0.66, theta(ITT)=0.72 |  |  |
| any(H) | 0.88 | 0.75 |
| sens(H) | 0.86 | 0.89 |
| sens(Hc) | 0.98 | 0.97 |
| ppv(H) | 0.89 | 0.84 |
| ppv(Hc) | 0.98 | 0.98 |
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
| FS: 880 estimable realizations (B = 300 bootstraps) |  |  |  |  |  |
|  | Avg | SD | Min | Max | Rel Bias (%) |
| H: 880 estimable, avg \|H\| = 86, theta(H) = 2 |  |  |  |  |  |
| theta-hat(H-hat) | 2.31 | 0.63 | 1.32 | 8.16 | 15.31 |
| Hc: avg \|Hc\| = 614, theta(Hc) = 0.66 |  |  |  |  |  |
| theta-hat(Hc-hat) | 0.64 | 0.07 | 0.45 | 0.91 | -3.04 |

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
(`n_sims = 1000`) will show higher Monte Carlo variability than the
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
