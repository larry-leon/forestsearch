# ForestSearch GLM Extension: Document Map and Naming System

## Naming Convention

```
{role}_{outcome}_{topic}.qmd
```

| Component | Values |
|---|---|
| **role** | `theory`, `calibration`, `selection`, `validation` |
| **outcome** | `glm` (all), `binary`, `survival`, `count`, `continuous` |
| **topic** | descriptive slug |


---

## Complete Inventory: Current → Proposed Names

### Theory (T)

| ID | Current Name | Proposed Name | Lines |
|---|---|---|---|
| T1 | `glm_detection_probability.qmd` | `theory_glm_detection_probability.qmd` | 784 |
| T2 | `null_approximation_correction.qmd` | `theory_glm_leff_correction.qmd` | 418 |


### Calibration (C) — expensive, render once per DGM

| ID | Current Name | Proposed Name | Lines | Output |
|---|---|---|---|---|
| C1 | `binary_threshold_calibration.qmd` | `calibration_binary_leff.qmd` | 974 | C=0.220, alpha=1.298 |
| C2 | `cox_vs_glm_approximation.qmd` | `calibration_survival_fourway.qmd` | 1077 | Validation only |
| C3 | `cox-vs-glm_approximation_null-grid.qmd` | `calibration_survival_leff_grid.qmd` | 793 | C=0.029, alpha=0.882 |
| -- | *MISSING* | `calibration_count_leff.qmd` | -- | TBD |


### Selection (S) — fast, uses prior L_eff

| ID | Current Name | Proposed Name | Lines | Outcome |
|---|---|---|---|---|
| S1 | `binary_threshold_calibration_null-search.qmd` | `selection_binary_frontier.qmd` | 1140 | Binary |
| S2 | `cox_vs_glm_approximation_null-search_3.qmd` | `selection_survival_frontier.qmd` | 907 | Survival |
| S3 | `suggest_thresholds_guide.qmd` | `selection_glm_guide.qmd` | 859 | All four |
| -- | *MISSING* | `selection_count_frontier.qmd` | -- | Count |


### Validation (V)

| ID | Current Name | Proposed Name | Lines |
|---|---|---|---|
| V1 | `glm_simulation_study2.qmd` | `validation_glm_simulation_study.qmd` | 919 |
| V2 | `hte_tests_vignette.qmd` | `validation_hte_tests_crump.qmd` | 1510 |


### R Source Files

| Current Name | Proposed Name | Lines | Purpose |
|---|---|---|---|
| `R/suggest_thresholds.R` | (keep) | 173 | Threshold frontier computation |
| `R/test_hte_crump.R` | (keep) | 419 | Crump et al. HTE tests |


---

## Per-Outcome Gap Analysis

### Survival (Cox PH)

| Pipeline Stage | Document | Status |
|---|---|---|
| d_eff formula | T1: `d_eff_survival(n_H, prop_cens)` | Done |
| L_eff calibration | C3: `calibration_survival_leff_grid.qmd` | Done (C=0.029, alpha=0.882) |
| Four-way validation | C2: `calibration_survival_fourway.qmd` | Done (Cox vs Poisson vs approx vs L_eff) |
| Threshold frontier | S2: `selection_survival_frontier.qmd` | Done (approximation-driven) |
| Guide entry | S3: `selection_glm_guide.qmd` | Done |
| Simulation OC | Not in V1 (survival uses separate GBSG vignette) | Done (main package vignette) |
| HTE pre-screen | V2: Poisson+offset on GBSG | Done |

**Missing:** Nothing. Survival pipeline is complete.


### Binary (Logistic)

| Pipeline Stage | Document | Status |
|---|---|---|
| d_eff formula | T1: `d_eff_binary(n_sg, p_event)` | Done |
| L_eff calibration | C1: `calibration_binary_leff.qmd` | Done (C=0.220, alpha=1.298) |
| Threshold frontier | S1: `selection_binary_frontier.qmd` | Done (approximation-driven) |
| Guide entry | S3: `selection_glm_guide.qmd` | Done |
| Simulation OC | V1: B-H0, B-H1-S, B-H0-L, B-H1-L | Done |
| HTE pre-screen | V2: OLS + logistic on simulation DGM | Done |

**Missing:** Nothing. Binary pipeline is complete.


### Count (Poisson + offset)

| Pipeline Stage | Document | Status |
|---|---|---|
| d_eff formula | T1: `d_eff_count(total_events)` | Done |
| L_eff calibration | *MISSING* | Uses binary L_eff as approximation |
| Threshold frontier | *MISSING* | Uses `suggest_thresholds()` with binary L_eff |
| Guide entry | S3: `selection_glm_guide.qmd` | Done (but with approximate L_eff) |
| Simulation OC | V1: K-H0, K-H1 | Done (K-H0 FPR ~7% with binary L_eff) |
| HTE pre-screen | V2: OLS + Poisson on simulation DGM | Done |

**Missing (HIGH priority):**

1. `calibration_count_leff.qmd` — Run ForestSearch under H0 with
   the count DGM (4 confounders: eos_low, smoking, male, age) at
   multiple N values (500, 700, 1000, 1500).  Derive count-specific
   C and alpha.  Structure follows C1 exactly.

2. `selection_count_frontier.qmd` — Once count L_eff is calibrated,
   produce an approximation-driven frontier.  Structure follows S1.

**Current workaround:** Binary L_eff (C=0.220, alpha=1.298) applied
to count outcomes.  This brought K-H0 FPR from ~27% (L_eff=1) to
~7% (with binary L_eff).  The approximation is reasonable because
both DGMs have 4 confounders, but is unvalidated.


### Continuous (Gaussian)

| Pipeline Stage | Document | Status |
|---|---|---|
| d_eff formula | T1: `d_eff_continuous(n_sg, sigma_y)` | Done |
| L_eff calibration | Not needed | FPR ~1% without correction |
| Threshold frontier | Not needed | `suggest_thresholds()` with L_eff=1 |
| Guide entry | S3: `selection_glm_guide.qmd` | Done |
| Simulation OC | V1: C-H0, C-H1 | Done (C-H0 FPR = 1%) |
| HTE pre-screen | V2: OLS on simulation DGM | Done |

**Missing:** Nothing.  The high d_eff for continuous outcomes
(large n_H / sigma^2) makes P1 so small that L_eff correction
is unnecessary.  The 1% FPR under H0 confirms this.

**Note:** Continuous thresholds are on the **mean-difference
scale** (c1=0.12, c2=0.12), qualitatively different from the
ratio-scale thresholds for binary/count/survival.  This is
documented in S3.


---

## Dependency Graph

```
T1 (d_eff, all families) ─────────────────────┐
T2 (L_eff correction) ────────────────────────┤
                                               │
C1 (binary L_eff) ─────────┐                  │
   C=0.220, alpha=1.298    │                  │
     S1 (binary frontier)  │                  │
                            ├── S3 (guide) ───┤
C3 (survival L_eff+grid) ─┐│                  │
   C=0.029, alpha=0.882   ││                  │
     S2 (survival frontier)││                  │
                            ││                  │
C2 (survival four-way) ───┘│                  │
                            │                  │
[calibration_count_leff] ──┘│  (MISSING)       │
     [selection_count_frontier]  (MISSING)     │
                            │                  │
V1 (simulation study) ─────┘                  │
V2 (HTE tests + GBSG) ───────────────────────┘
```


---

## Proposed File Layout

```
quarto/
|
|  # Theory
+-- theory_glm_detection_probability.qmd          (T1)
+-- theory_glm_leff_correction.qmd                (T2)
|
|  # Calibration (expensive, render once)
+-- calibration_binary_leff.qmd                   (C1)
+-- calibration_survival_fourway.qmd              (C2)
+-- calibration_survival_leff_grid.qmd            (C3)
+-- [calibration_count_leff.qmd]                  MISSING
|
|  # Selection (fast, uses prior L_eff)
+-- selection_binary_frontier.qmd                 (S1)
+-- selection_survival_frontier.qmd               (S2)
+-- selection_glm_guide.qmd                       (S3)
+-- [selection_count_frontier.qmd]                MISSING
|
|  # Validation
+-- validation_glm_simulation_study.qmd           (V1)
+-- validation_hte_tests_crump.qmd                (V2)
|
|  # Persisted results
+-- _output/
    +-- binary_leff_calibration.rds
    +-- survival_leff_calibration.rds
    +-- binary_frontier.rds
    +-- survival_frontier.rds
```


---

## Renaming Checklist

If you choose to rename, do these in order to avoid breaking
cross-references:

| Step | Old | New |
|---|---|---|
| 1 | `glm_detection_probability.qmd` | `theory_glm_detection_probability.qmd` |
| 2 | `null_approximation_correction.qmd` | `theory_glm_leff_correction.qmd` |
| 3 | `binary_threshold_calibration.qmd` | `calibration_binary_leff.qmd` |
| 4 | `cox_vs_glm_approximation.qmd` | `calibration_survival_fourway.qmd` |
| 5 | `cox-vs-glm_approximation_null-grid.qmd` | `calibration_survival_leff_grid.qmd` |
| 6 | `binary_threshold_calibration_null-search.qmd` | `selection_binary_frontier.qmd` |
| 7 | `cox_vs_glm_approximation_null-search_3.qmd` | `selection_survival_frontier.qmd` |
| 8 | `suggest_thresholds_guide.qmd` | `selection_glm_guide.qmd` |
| 9 | `glm_simulation_study2.qmd` | `validation_glm_simulation_study.qmd` |
| 10 | `hte_tests_vignette.qmd` | `validation_hte_tests_crump.qmd` |

After renaming, update cross-references in:
- C_prior / alpha_prior comments citing source documents
- output_dir paths in the selection documents
- Any `source()` calls pointing to .qmd filenames


---

## Reading Order for a New Reader

1. T1: Theory -- d_eff framework
2. T2: Theory -- L_eff correction
3. S3: Guide -- practical threshold selection (all outcomes)
4. S1: Selection -- binary deep-dive
5. S2: Selection -- survival deep-dive
6. V1: Validation -- simulation study
7. V2: Validation -- HTE pre-screening + GBSG analysis

Calibration documents (C1, C2, C3) can be skipped by readers
who accept the prior L_eff values.
