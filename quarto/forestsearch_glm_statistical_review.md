# ForestSearch GLM Extension: Comprehensive Statistical Review

**Author:** Review of 14 Quarto documents by Larry F. León  
**Date:** April 2026  
**Package:** `forestsearch` (v0.1.0 on CRAN; GLM extension targeting v0.2.0)

---

## 1. Executive Summary

These 14 documents form a rigorous, interconnected body of work extending the ForestSearch methodology (León et al., 2024, *Statistics in Medicine*) from survival endpoints to general GLM outcomes — binary (logistic), continuous (Gaussian), and count (Poisson) responses. The extension rests on three foundational pillars:

1. **A unified detection probability theory** via the effective information $d_{\text{eff}}$, which generalizes the survival score-statistic variance $8/d$ to arbitrary GLM families.

2. **A multiplicity correction** via $L_{\text{eff}}$, the effective number of approximately independent candidate subgroups, which accounts for the gap between per-subgroup and procedure-level false positive rates.

3. **Simulation-calibrated threshold selection**, where an approximation-guided frontier identifies candidate $(c_1, c_2)$ pairs, validated by Monte Carlo simulation.

The documents are organized into four functional roles: **Theory** (T1–T2), **Calibration** (C1–C3), **Selection** (S1–S3), and **Validation** (V1–V2), with two additional literature reviews and a biomarker comparison study.

---

## 2. Theoretical Foundations

### 2.1 The Section 2.1 Detection Probability (T1: `theory_glm_detection_probability.qmd`)

The core approximation from León et al. (2024, eq. 3) gives the probability that ForestSearch identifies a specific subgroup $H$ with true effect $\beta = \log(\theta^{\dagger}(H))$:

$$
p(c_1, c_2;\, d,\, \theta^{\dagger}(H)) = P\bigl(W_1 + W_2 \geq 2k_1,\;\min(W_1, W_2) \geq k_2\bigr)
$$

where $W_1, W_2 \stackrel{\text{iid}}{\sim} N(\beta,\, 8/d)$ are the split-half log-effect estimates, $k_1 = \log(c_1)$, and $k_2 = \log(c_2)$.

The two conditions encode ForestSearch's dual criteria:

- **Screening**: the average effect across halves exceeds $c_1$ (i.e., $\hat{\beta}_s = (W_1 + W_2)/2 \geq k_1$)
- **Consistency**: each half independently exceeds $c_2$ (i.e., $\min(W_1, W_2) \geq k_2$)

#### Generalization to GLM via $d_{\text{eff}}$

The key insight is that for any canonical GLM with balanced arms ($n_1 = n_0 = n_H/2$), the Fisher information for the treatment coefficient takes the form $I(\beta) = (n_H/4)\cdot\kappa$, yielding:

$$
\operatorname{Var}(\hat{\beta}_{\text{full}}) = \frac{4}{d_{\text{eff}}}, \qquad
\operatorname{Var}(\hat{\beta}_{\text{half}}) = \frac{8}{d_{\text{eff}}}
$$

The family-specific $d_{\text{eff}}$ derivations are:

| Family | Link | Fisher Info Factor $\kappa$ | $d_{\text{eff}}$ | Effect Scale |
|--------|------|------|------|------|
| Cox PH | $\log$ | $(1-p_c)$ | $n_H(1-p_c) = d$ (events) | $\log$ |
| Poisson + offset | $\log$ | $\bar{\mu}$ (mean rate) | $D$ (total events) | $\log$ |
| Binomial | logit | $\bar{p}(1-\bar{p})$ | $n_H \bar{p}(1-\bar{p})$ | $\log$ |
| Gaussian | identity | $1/\sigma^2$ | $n_H / \sigma^2$ | identity |
| Quasi-Poisson | $\log$ | $\bar{\mu}/\hat{\phi}$ | $D/\hat{\phi}$ | $\log$ |

**Critical distinction for continuous outcomes**: thresholds are on the identity scale ($k_1 = c_1$, $k_2 = c_2$), not the log scale, since the effect $\beta$ is a mean difference rather than a log-ratio.

The document provides a unified implementation `compute_detection_probability_glm()` with companion helpers `d_eff_survival()`, `d_eff_binary()`, `d_eff_count()`, and `d_eff_continuous()`, validated against the existing Cox-specific `compute_detection_probability()` with numerical agreement to four decimal places.

#### Cross-Family Comparison

For $n_H = 100$ with $c_1 = 1.25$, $c_2 = 1.0$:

| Outcome | $d_{\text{eff}}$ | $P(\theta=1.0)$ | $P(\theta=2.0)$ | 80% Power at |
|---------|------------------|------------------|------------------|--------------|
| Survival (45% cens) | 55 | ~15% | ~82% | $\theta \approx 1.9$ |
| Poisson ($\mu=0.8$) | 80 | ~9% | ~92% | $\theta \approx 1.6$ |
| Binary ($p=0.30$) | 21 | ~23% | ~50% | $\theta \approx 2.8$ |
| Binary ($p=0.50$) | 25 | ~21% | ~57% | $\theta \approx 2.5$ |

This reveals that binary outcomes have substantially lower $d_{\text{eff}}$ than survival for the same $n_H$, explaining why the paper's survival-calibrated defaults ($c_1 = 1.25$, $c_2 = 1.0$) yield inadequate FPR control for binary data.


### 2.2 The $L_{\text{eff}}$ Multiplicity Correction (T2: `theory_glm_leff_correction.qmd`)

The per-subgroup approximation $P_1$ models a single pre-specified subgroup. Under $H_0$ the procedure evaluates many candidates, so the procedure-level FPR is:

$$
\boxed{\,
P_{\text{null}}^{\text{corrected}}(N) = 1 - (1 - P_1)^{L_{\text{eff}}(N)}
\,}
$$

where $L_{\text{eff}}$ is the effective number of approximately independent candidates. This is **not** the raw candidate count but the number of effectively independent "rolls of the die."

#### Why $L_{\text{eff}}$ Grows with $N$

Three mechanisms drive this growth:

1. **More subgroup combinations survive `n.min`**: at $N = 200$ with `n.min = 60`, a combination needing $\geq 30\%$ prevalence is excluded; at $N = 1000$ the same combination enters the candidate set.
2. **GRF generates more split points** on continuous variables with larger samples.
3. **LASSO retains more variables** due to increased power for spurious prognostic signals.

#### Power-Law Model

The implied $L_{\text{eff}}$ is recovered from simulation via:

$$
L_{\text{eff}} = \frac{\log(1 - \text{FPR}_{\text{sim}})}{\log(1 - P_1)}
$$

and modeled as a power law:

$$
L_{\text{eff}}(N) = C \cdot \left(\frac{N}{n_{\min}}\right)^{\alpha}
$$

Key property: $L_{\text{eff}}$ depends on the candidate generation process (confounders, GRF, LASSO, `maxk`, `n.min`) but **not** on the thresholds $(c_1, c_2)$. The thresholds control $P_1$; $L_{\text{eff}}$ controls how many times $P_1$ is "rolled." This threshold-independence is validated empirically across multiple $(c_1, c_2)$ pairs in the calibration documents.

#### Why the Correction Was Unnecessary in the Original Paper

For the GBSG survival DGM with $N = 300$–700 and `n.min = 60`, the events-based filters (`d0.min`, `d1.min`) are more restrictive, giving $L_{\text{eff}} \approx 0.3$–1.5. The per-subgroup approximation already slightly overestimates the FPR, making the defaults conservative. For binary outcomes at larger $N$, $L_{\text{eff}}$ can reach 5–15, making the correction essential.

---

## 3. Calibration Studies

### 3.1 Binary Outcome Calibration (C1: `calibration_binary_leff.qmd`)

This document serves as the primary calibration source for binary outcomes:

- **DGM**: Logistic model with 4 confounders (bm1, bm2, age, ECOG), $p_{\text{event}} = 0.30$, 30% harm prevalence, $\text{OR}_H = 2.0$.
- **Simulation**: ForestSearch under $H_0$ and $H_1$ across $N \in \{300, 500, 700\}$ and a grid of $(c_1, c_2)$ pairs.
- **Calibration output**: $L_{\text{eff}}$ model fitted via `calibrate_L_eff()` and saved as `.rds` for downstream use.

The document demonstrates the full workflow: (1) survival benchmark reproduction, (2) binary $d_{\text{eff}}$ calculation showing the information deficit, (3) simulation engine with parallel execution, (4) approximation vs. simulation comparison, (5) $L_{\text{eff}}$ calibration and corrected FPR.

Key findings: power agreement (Sim/Approx $\approx 1.0$ for $N \geq 500$) vs. FPR disagreement (Sim/Approx reaching 5–15× without correction), resolved to $\approx 1.0$ with the $L_{\text{eff}}$ correction.

### 3.2 Survival Four-Way Validation (C2: `calibration_survival_fourway.qmd`)

This is the most computationally intensive document, running four parallel calculations to validate the GLM extension through the Poisson–Cox equivalence:

| Label | Method | Source |
|-------|--------|--------|
| **(A)** | ForestSearch Poisson+offset simulation | GLM pipeline |
| **(B)** | GLM approximation: $d_{\text{eff}} = D$ | GLM theory |
| **(C)** | ForestSearch Cox PH simulation | Established pipeline |
| **(D)** | Cox approximation: `compute_detection_probability()` | Package |

The theoretical basis: under proportional hazards with piecewise-constant baseline, $d_{\text{eff}}^{\text{Poisson}} = D = d^{\text{Cox}}$ (Laird & Olivier, 1981; Whitehead, 1980). Therefore (B) = (D) validates the variance formula, and (A) = (C) validates the GLM pipeline.

Validation confirms: near-identical detection rates, classification metrics (sensitivity, specificity, PPV), and per-simulation concordance between Cox and Poisson pipelines ($>95\%$ agreement). The $L_{\text{eff}}$ calibration for the GBSG DGM (7 confounders) is also computed and saved.

### 3.3 Survival $L_{\text{eff}}$ Grid (C3: `calibration_survival_leff_grid.qmd`)

A faster grid-sweep variant of C2, designed to verify $L_{\text{eff}}$ threshold-independence across a broader $(c_1, c_2)$ grid with fewer replicates. Acts as a complementary check; C2 is the primary calibration source.

---

## 4. Threshold Selection Framework

### 4.1 The `suggest_thresholds()` Function (S3: `suggest_thresholds_guide.qmd`, `selection_glm_guide.qmd`)

The function implements approximation-driven threshold selection:

$$
\text{FPR}_{\text{corr}} = 1 - (1 - P_1)^{L_{\text{eff}}(N)} \leq \text{target}
$$

It returns a **frontier table**: for each $c_1$ (ascending), the minimum $c_2$ that achieves the FPR target. The function handles both ratio-scale effects (HR, OR, IRR with $k_j = \log(c_j)$) and difference-scale effects (MD with $k_j = c_j$ directly).

#### Recommended Starting Thresholds

| Outcome | $c_1$ | $c_2$ | Scale | Notes |
|---------|-------|-------|-------|-------|
| Survival | 1.25 | 1.00 | HR | Paper defaults; safe for GBSG-like settings |
| Binary ($p \approx 0.30$) | 1.55–1.60 | 1.50–1.55 | OR | Requires $L_{\text{eff}}$ correction |
| Count (high rate, $\mu \geq 3$) | 1.10–1.25 | 0.90–1.00 | IRR | High $d_{\text{eff}}$; defaults may suffice |
| Count (low rate, $\mu < 1$) | 1.25–1.40 | 1.00–1.15 | IRR | Similar to survival |
| Continuous ($\sigma = 1$) | 0.03–0.10 | 0.01–0.05 | Mean diff | Cohen's $d$ scale when standardized |

### 4.2 Binary Frontier Selection (S1: `selection_binary_frontier.qmd`)

Implements the full threshold selection pipeline for binary outcomes:

1. Loads $L_{\text{eff}}$ calibration from C1's `.rds` bundle
2. Computes the corrected FPR frontier via `suggest_thresholds()`
3. Validates selected pairs with ForestSearch $H_0$ simulation
4. Generates power–FPR tradeoff plots
5. Produces optimal pair recommendations

The document demonstrates the parameter-inheritance pattern: S1 loads C1's DGM configuration (confounders, event rate, ForestSearch parameters) to ensure consistency.

### 4.3 Survival Frontier Selection (S2: `selection_survival_frontier.qmd`)

The survival analog of S1, loading calibration from C2's `.rds` and applying `suggest_thresholds()` to the GBSG DGM. Confirms that the paper defaults lie safely within the frontier for survival settings.

---

## 5. Pre-Screening for Treatment Effect Heterogeneity

### 5.1 Crump et al. (2008) HTE Tests (V1: `validation_hte_tests_crump.qmd`)

This document implements and validates a pre-screening step: before running ForestSearch, test whether HTE exists at all. Three exported functions are provided:

**Test 1 — Zero CATE ($H_0: \tau(x) = 0\;\forall x$):**

$$
Q = (\hat{\xi}_{1,K} - \hat{\xi}_{0,K})'\,\hat{V}^{-1}\,(\hat{\xi}_{1,K} - \hat{\xi}_{0,K}) \;\xrightarrow{d}\; \chi^2(K) \quad\text{under}\; H_0
$$

where $\hat{\xi}_{w,K}$ are the sieve regression coefficients fitted separately in each treatment arm using a polynomial basis $P_K(x) = (1, x_1, \ldots, x_d)'$, and $\hat{V}$ is the HC sandwich variance of their difference.

**Test 2 — Constant CATE ($H_0': \tau(x) = \tau\;\forall x$):**

$$
Q' = (\hat{\xi}_{1,1,K} - \hat{\xi}_{0,1,K})'\,\hat{V}_{11}^{-1}\,(\hat{\xi}_{1,1,K} - \hat{\xi}_{0,1,K}) \;\xrightarrow{d}\; \chi^2(K-1) \quad\text{under}\; H_0'
$$

The slope coefficients (indices $2:K$) are tested while the intercept (which captures the overall ATE) is excluded.

**Key distinction**: the standard ATE test ($H_0'': E[\tau(X)] = 0$) averages over subpopulations and can miss opposing effects. The Crump tests detect heterogeneity even when harm and benefit cancel on average.

The implementation extends the original OLS framework to logistic and Poisson+offset regression using the GLM sandwich:

$$
\hat{\Omega}_w = N_w \cdot A^{-1} B\, A^{-1}
$$

where $A$ is the expected Fisher information with working weights and $B$ is the empirical score outer product. For Poisson: $w_i = \mu_i$ and $s_i = P_i(Y_i - \mu_i)$. For logistic: $w_i = \mu_i(1 - \mu_i)$ and $s_i = P_i(Y_i - \mu_i)$.

The document includes a 12-scenario Monte Carlo study covering size calibration (6 null scenarios) and power (6 alternative scenarios) across OLS, logistic, and Poisson+offset specifications, with variance calibration diagnostics and GBSG real-data analysis.

---

## 6. Causal Foundations for Survival Estimation

### 6.1 Cox/AFT Causal Review (V2: `cox_causal_review.qmd`)

This document provides a comprehensive treatment of causal estimands for survival data, synthesizing Hernán (2010), Aalen, Cook & Røysland (2015), Martinussen (2022), Fay & Li (2024), and Knudsen et al. (2025).

#### The Weibull AFT/Cox Dual Representation

The ForestSearch DGM uses the Weibull AFT model:

$$
\log(T_i) = \mu + \mathbf{X}_i^{\top}\boldsymbol{\gamma} + \sigma\,\varepsilon_i, \qquad \varepsilon_i \sim \text{Extreme Value}
$$

which is equivalently a Cox PH model with the AFT-to-hazard transformation:

$$
\boldsymbol{\beta}_0 = -\boldsymbol{\gamma}\,/\,\sigma
$$

This duality means potential outcomes are **simulated** on the AFT scale (guaranteeing exact exchangeability) while effects are **interpreted** on the HR scale (connecting to clinical convention).

#### Estimand Taxonomy (Fay & Li, 2024)

| Level | Estimand | Definition | Causal Status |
|-------|----------|------------|---------------|
| Individual | $\text{HR}_i(t)$ | $\lambda_i^{(1)}(t)/\lambda_i^{(0)}(t)$ | Valid under AFT |
| Population | AHR | $\exp\bigl(E[\log(\text{HR}_i)]\bigr)$ | Valid (geometric mean) |
| Population | Marginal Causal HR | $\bar{\lambda}^{(1)}(t)/\bar{\lambda}^{(0)}(t)$ | Varies with $t$ |
| Population | CDE | $E[T_i(0)]/E[T_i(1)]$ | Natural-scale complement |

The **Average Hazard Ratio (AHR)** is the primary simulation estimand because it inherits dual individual-/population-level causal validity from the AFT scale-change parameter and is deterministic given the super-population. The CDE provides a complementary natural-scale perspective.

#### The Built-In Selection Bias Problem

The instantaneous Cox HR at time $t$ conditions on survival to $t$, which is affected by treatment assignment. Even in a perfectly randomized trial, the risk sets at $t > 0$ are no longer exchangeable — sicker treated patients have been "selected out," creating a spurious attenuation of the observed HR over time. This is resolved by the AFT framework, where the causal parameter is the time-scale change, not the instantaneous hazard ratio.

---

## 7. Biomarker Evaluation and Clinical Motivation

### 7.1 Continuous Biomarker Comparison (V2: `biomarker_comparison.qmd`)

Compares three methods for evaluating treatment effect modification by continuous biomarkers:

1. **`cox_cs_fit()`** — natural cubic spline interaction in a Cox model, providing a smooth log-HR profile $\hat{\beta}(z)$ over the biomarker range with pointwise confidence bands
2. **STEPP** — Subpopulation Treatment Effect Pattern Plot using overlapping sliding windows of patients, with permutation-based tests for HTE
3. **MFPI** — Multivariable Fractional Polynomial Interaction using data-driven power transforms of the biomarker

Simulation study (5000 replicates, N = 700) evaluates bias, coverage, and CI width across 8 biomarker values. STEPP permutation power study (500 replicates × 999 permutations) compares test power to the spline Wald test.

### 7.2 GLM HTE Literature Summary (`glm_hte_literature_summary.qmd`)

Comprehensive survey organized by outcome type:

- **Binary**: Tumor response (RECIST), ACR20/50/70, virologic suppression. Key methods: Virtual Twins, mCART. DINA estimand = log-OR.
- **Continuous**: DAS28 change, viral load change, FEV1 change. Analyzed via ANCOVA/MMRM. DINA estimand = mean difference.
- **Count**: COPD exacerbations (MATINEE, GALATHEA/TERRANOVA, BOREAS), MS relapses, HF hospitalizations. Standard model: negative binomial. DINA estimand = log-IRR.

### 7.3 Count Data Framework (`count_data_hte_summary.qmd`)

Detailed treatment of count/rate outcomes for the Poisson extension:

- **Clinical DGM**: Respiratory exacerbation trial ($N = 1000$, control rate $\sim$1.2/year, harm subgroup IRR $\approx 2.7$, complement IRR $\approx 0.61$)
- **GRF screening**: $\log(Y + 0.5)$ variance-stabilizing transform before causal forest splits
- **Overdispersion**: quasi-Poisson ($\hat{\phi} > 1$) and negative binomial both produce identical point estimates for log-IRR; only SEs differ
- **Poisson–Cox connection**: under proportional hazards with constant baseline, $\text{IRR}_{\text{Poisson}} \approx \text{HR}_{\text{Cox}}$

---

## 8. Document Architecture and Workflow

The documents follow a four-role taxonomy with a strict dependency chain:

```
Theory (T1, T2) → Calibration (C1, C2, C3) → Selection (S1, S2, S3) → Validation (V1, V2)
```

Data flows via `.rds` bundles:

- C1 saves `calibration_binary_leff.rds` → S1 loads
- C2 saves `calibration_survival_leff.rds` → S2 loads
- Bundles contain: $L_{\text{eff}}$ parameters, DGM config, ForestSearch settings

Each document uses a centralized parameter block in the setup chunk, with short-name aliases in an `eval-params` chunk. The `tryCatch` fallback pattern ensures documents render even without cached calibration data.

---

## 9. Suggested Improvements and Future Work

### 9.1 Immediate Priorities (v0.2.0 Blocking)

1. **Count outcome routing**: The open architectural question of whether `family = Poisson` routes through `outcome_type = "binary"` with `effect_measure = "IRR"` or has a dedicated `"count"` value needs resolution. The documents consistently use `outcome_type = "count"` and `effect_measure = "IRR"`, which suggests a dedicated value is preferred. GRF screening for true count outcomes (integers > 1) should be validated separately.

2. **NAMESPACE exports**: Functions referenced across documents — `grf.subg.harm.glm()`, `glm_cs_fit()`, `plot_sg_glm_outcomes()`, `suggest_thresholds()`, `calibrate_L_eff()`, `save_leff_calibration()`, `load_leff_calibration()` — need formal `@export` tags and NAMESPACE registration.

3. **Continuous outcome validation gap**: The continuous outcome ($d_{\text{eff}} = n_H/\sigma^2$) has no dedicated calibration study or simulation validation. The `suggest_thresholds_guide.qmd` uses it in examples, but $L_{\text{eff}}$ is uncalibrated.  A minimal C-level document with a continuous DGM (e.g., DAS28 change-from-baseline) would close this gap.

### 9.2 Methodological Extensions

4. **$L_{\text{eff}}$ as a function of `maxk` and confounder count**: The power-law model $L_{\text{eff}} = C(N/n_{\min})^{\alpha}$ treats $C$ and $\alpha$ as constants for a given confounder configuration. A meta-model predicting $C$ and $\alpha$ from `maxk`, $|\text{confounders}|$, and the number of binary vs. continuous confounders would reduce the need for per-DGM calibration simulations. The two existing calibrations (survival: 7 confounders, $C \approx 0.029$, $\alpha \approx 0.88$; binary: 4 confounders, $C \approx 0.22$, $\alpha \approx 1.30$) suggest that fewer confounders paradoxically increase $L_{\text{eff}}$ at a given $N$, possibly because events-based filters (`d0.min`, `d1.min`) in survival settings prune more candidates.

5. **Negative binomial support**: The count data framework discusses quasi-Poisson and negative binomial but the ForestSearch pipeline currently supports only Poisson+offset. Adding `MASS::glm.nb()` as an alternative for overdispersed count data would be clinically relevant for COPD/MS trials. The $d_{\text{eff}}$ adjustment is $D/\hat{\phi}$, which is already documented for quasi-Poisson.

6. **Non-proportional effects**: The detection probability theory assumes constant treatment effect within a subgroup. When effects are time-varying (non-proportional hazards) or dose-dependent, the split-half variance $8/d_{\text{eff}}$ may not capture the estimation uncertainty accurately. The CHR (Cumulative Hazard Ratio) discussed in the causal review could serve as an alternative estimand for non-PH settings.

7. **Integration of Crump HTE pre-tests**: The pre-screening tests (`test_zero_cate`, `test_constant_cate`, `test_hte`) are fully implemented but not yet integrated into the `forestsearch()` pipeline. A natural extension would be a `test_hte_from_forestsearch()` wrapper that applies the test using ForestSearch's own screening-selected covariates, with the option to abort the subgroup search if $p > 0.10$.

### 9.3 Simulation and Validation Gaps

8. **GLM-native DGM**: The current binary and count simulations use `forestsearch()` with `outcome_type = "binary"` or `"count"`, but the DGMs are ad-hoc R functions rather than package-level generators like `setup_gbsg_dgm()`. A formal `setup_glm_dgm()` function with calibration infrastructure (analogous to `calibrate_k_inter()` for survival) would standardize the simulation workflow.

9. **Bootstrap validation for GLM**: The bootstrap bias correction step (`forestsearch_bootstrap_dofuture()`) has been validated extensively for survival but lacks a dedicated GLM simulation study showing that bootstrap-corrected estimates are closer to truth than raw estimates for binary/count/continuous outcomes.

10. **Cross-validation performance**: Similarly, `forestsearch_Kfold()` for GLM outcomes needs a simulation study characterizing its expected coverage and bias properties across outcome types.

### 9.4 Documentation and Packaging

11. **Unified parameter naming**: The documents use `hr.threshold` and `hr.consistency` for all outcome types, which is confusing for non-HR effects. An alias system where `effect.threshold` maps to `hr.threshold` internally would improve the API for GLM users without breaking backward compatibility.

12. **Label generalization**: Several utility functions (`build_estimation_table()`, `compare_metrics()`, `format_oc_results()`) contain survival-centric labels (e.g., "HR", "Events"). An `effect_label` parameter with outcome-specific defaults would make outputs correct for GLM analyses.

13. **Vignette packaging**: The 14 `.qmd` documents collectively form a comprehensive methodological supplement. A condensed CRAN vignette covering (a) the $d_{\text{eff}}$ framework, (b) a binary outcome worked example, and (c) the `suggest_thresholds()` workflow would serve as the user-facing entry point, with the full document set available via `pkgdown`.

### 9.5 Additional Theoretical Directions

14. **Adaptive $L_{\text{eff}}$**: Currently $L_{\text{eff}}$ is calibrated from null simulations and treated as fixed. An online estimate using the observed number of candidate subgroups (from `find.grps$out.found`) could provide a data-dependent correction without requiring a separate calibration simulation.

15. **Bayesian detection probability**: The frequentist $P_1$ framework could be extended to incorporate prior information on the effect size distribution. A Bayesian detection probability integrating over a prior $\pi(\theta)$ would be useful for Phase II trials where external information about treatment effect heterogeneity is available.

16. **Multi-endpoint extension**: Many clinical trials collect multiple endpoints (e.g., binary response + continuous biomarker change + time-to-event). A joint detection probability accounting for the correlation structure across endpoints would allow simultaneous subgroup identification, avoiding the current practice of running ForestSearch separately for each endpoint.

---

## 10. Summary Assessment

This body of work represents a rigorous and well-structured extension of the ForestSearch methodology. The mathematical foundation is sound: the $d_{\text{eff}}$ generalization follows naturally from Fisher information theory, and the $L_{\text{eff}}$ correction addresses a real and well-characterized gap between per-subgroup and procedure-level error rates. The four-way Poisson–Cox validation provides a compelling internal consistency check.

The main areas for strengthening are: (a) closing the continuous-outcome calibration gap, (b) formalizing the GLM-native DGM and bootstrap validation, and (c) packaging the extensive document set into a user-accessible vignette hierarchy. The theoretical framework is ready for v0.2.0; the remaining work is primarily engineering and validation.

---

## References

- León, L.F., Jemielita, T., Guo, Z., Marceau West, R., Anderson, K.M. (2024). Exploratory subgroup identification in the heterogeneous Cox model. *Statistics in Medicine*, 43, 4253–4275. doi:10.1002/sim.10163
- Crump, R.K., Hotz, V.J., Imbens, G.W., Mitnik, O.A. (2008). Nonparametric tests for treatment effect heterogeneity. *Review of Economics and Statistics*, 90(3), 389–405.
- Fay, M.P. and Li, F. (2024). Causal estimands for time-to-event outcomes. *Lifetime Data Analysis*, 30, 588–612.
- Hernán, M.A. (2010). The hazards of hazard ratios. *Epidemiology*, 21(1), 13–15.
- Aalen, O.O., Cook, R.J., Røysland, K. (2015). Does Cox analysis of a randomized survival study yield a causal treatment effect? *Lifetime Data Analysis*, 21, 579–593.
- Jennison, C. and Turnbull, B.W. (1984). Repeated confidence intervals for group sequential clinical trials. *Controlled Clinical Trials*, 5(1), 33–45.
- McCullagh, P. and Nelder, J.A. (1989). *Generalized Linear Models*, 2nd ed. Chapman and Hall.
- Laird, N. and Olivier, D. (1981). Covariance analysis of censored survival data using log-linear analysis techniques. *JASA*, 76(374), 231–240.
- DINA (2025). Estimating heterogeneous treatment effects for general responses. *Biometrics*, 81(4).
