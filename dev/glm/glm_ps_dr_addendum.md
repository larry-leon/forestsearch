# GLM Extension — Propensity-Score & Doubly-Robust Addendum

**Context:** forestsearch GLM extension, `feature/glm-extension` branch
**References:** Bang & Robins (2005, Biometrics); Funk et al. (2011, Am J Epidemiol)
**Code basis:** `BangRobins_TreatmentEffectEstimationFunctions.R` (Larry Leon)
**Date:** 2026-03-24

---

## 1. Motivation

ForestSearch operates in three settings where propensity-score methods add
value beyond the current unadjusted estimators:

1. **Covariate adjustment within subgroups.** Even in a randomized trial,
   conditioning on the subgroup membership indicator (which is a function of
   covariates) can induce residual covariate imbalance within the identified
   subgroup.  Covariate-adjusted estimators recover efficiency.

2. **Multi-regional clinical trials (MRCT).** Region-level treatment allocation
   ratios may differ from the nominal 1:1, and region is itself a stratification
   factor.  A PS model that includes region provides automatic re-weighting.

3. **Non-randomized comparisons.** ForestSearch's consistency-based subgroup
   identification is applicable to observational data (e.g., the cross-trial
   Pola vs. RCHOP comparison in Larry's presentation), where confounding
   adjustment is essential.

In all three cases the relevant estimand is the **risk difference (RD)** —
computed via g-computation (standardization) over the covariate distribution,
not from the model coefficient.

---

## 2. Three Estimator Tiers

The `adjustment` parameter threads through `make_effect_estimator()` and
controls the estimation strategy.  Each tier inherits the same consistency
criterion and bootstrap infrastructure.

### Tier 0: Unadjusted (`adjustment = "none"`)

Current default.  The closure fits `Y ~ treat` within each data slice and
returns the treatment coefficient (on log or identity scale depending on
`effect_measure`).  This is sufficient for randomized trials with the
unadjusted RD, OR, etc.

### Tier 1: PS-Covariate / Bang & Robins DR (`adjustment = "ps_covariate"`)

From Bang & Robins (2005), as implemented in `BR.TE.fit()`:

```
1.  Fit PS model:  π̂(L) = P(treat = 1 | L)   via logistic regression
2.  Compute inverse PS:
      inv_ps_i = 1/π̂_i        if treat_i = 1
      inv_ps_i = 1/(1 - π̂_i)  if treat_i = 0
3.  Fit augmented outcome model:
      Y ~ treat + confounders + inv_ps    [family per effect_measure]
4.  Predict for each arm:
      μ̂₁(L_i) = predict(fit, treat = 1, L_i, inv_ps = 1/π̂_i)
      μ̂₀(L_i) = predict(fit, treat = 0, L_i, inv_ps = 1/(1-π̂_i))
5.  Risk difference:  RD = mean(μ̂₁) - mean(μ̂₀)
```

This is the simplest DR approach.  The estimate is consistent if **either**
the PS model or the outcome model is correctly specified.  It maps directly to
the existing `BR.TE.fit(..., get.main = TRUE)` code path.

**Key implementation detail:** The inverse PS is added as a *linear* covariate
in the outcome model.  Larry's code also implements non-linear versions via
GAM splines (`s(psE)`, `s(psC)` — v5/v6 in `BR.TE.fit()`), but the linear
version is the starting point for forestsearch.

### Tier 2: AIPW / Funk & Davidian (`adjustment = "aipw"`)

From Funk et al. (2011), as implemented in v3 of `BR.TE.fit()`:

```
1.  Fit PS model:  π̂(L) = P(treat = 1 | L)
2.  Fit outcome regression:  μ̂(Z, L) via GLM
3.  Predict potential outcomes:
      μ̂₁(L_i) = predict(fit, treat = 1, L_i)
      μ̂₀(L_i) = predict(fit, treat = 0, L_i)
4.  Compute AIPW influence-function estimates:

    φ̂₁_i = Z_i Y_i / π̂_i  −  (Z_i − π̂_i)/π̂_i · μ̂₁(L_i)

    φ̂₀_i = (1−Z_i) Y_i / (1−π̂_i)  +  (Z_i − π̂_i)/(1−π̂_i) · μ̂₀(L_i)

5.  RD_AIPW = mean(φ̂₁) − mean(φ̂₀)
```

This is the influence-function-based DR estimator.  It has the same double-
robustness property as the PS-covariate approach but also achieves
the semiparametric efficiency bound when both models are correct.

---

## 3. Architecture: Where the PS Model Is Fitted

### Option A: Fit once, store as column (recommended for forestsearch)

Fit `π̂(L)` once on the full `df.analysis` before entering the consistency
loop.  Store `inv_ps` as a column.  The estimator closure picks it up from
each data slice.

**Pros:**
- Consistent with the principle that the PS model reflects the full-study
  covariate structure, not a random half of a subgroup
- Much faster (one `glm()` fit vs. 2 × R × n_subgroups fits)
- Exactly replicates Larry's existing code pattern (`data$inv.ps <- g.DR`)

**Cons:**
- The PS model is not re-estimated in the bootstrap, which means the
  bootstrap does not capture uncertainty in the PS model.  This is a known
  limitation addressed in the AIPW variance formula.

**Implementation:** In `forestsearch()`, if `adjustment != "none"`:
- Fit `ps_model <- glm(treat ~ confounders, family = binomial, data = df.analysis)`
- Compute `pihat <- fitted(ps_model)`
- Truncate: `pihat <- pmax(pmin(pihat, 1 - ps_truncate), ps_truncate)`
  (default `ps_truncate = 0.01`)
- Store: `df.analysis$inv_ps <- ifelse(treat == 1, 1/pihat, 1/(1-pihat))`
- Pass `ps.name = "inv_ps"` to `make_effect_estimator()`

### Option B: Fit within each closure call

Re-fit the PS model on each data slice.  More statistically rigorous for the
bootstrap but prohibitively expensive for the consistency loop (R = 400 splits
× n_subgroups × 2 halves = tens of thousands of PS model fits).

**Recommendation:** Option A for the consistency loop; consider Option B only
for the final bootstrap bias-correction step where computational budget is
larger.

---

## 4. Estimator Closure Signatures

### PS-Covariate Estimator

```r
.make_glm_ps_covariate_estimator <- function(
    treat.name,
    outcome.name,
    confounders.name,
    ps.name        = "inv_ps",
    effect_measure = "RD",
    family_fn      = stats::binomial(link = "identity"),
    ...
)
```

The closure receives a `data_slice` that already contains the `inv_ps` column.
It fits the augmented model, predicts under each arm, and returns the
g-computed RD (or OR, etc.).

### AIPW Estimator

```r
.make_aipw_estimator <- function(
    treat.name,
    outcome.name,
    confounders.name,
    pihat.name     = ".pihat",
    effect_measure = "RD",
    family_fn      = stats::binomial(link = "logit"),
    ...
)
```

The closure fits the outcome model on the data slice, predicts potential
outcomes, then applies the AIPW formula.  Needs the raw `π̂` values
(not just `inv_ps`) to compute the influence-function terms.

---

## 5. Consistency Criterion with Adjusted Estimators

The consistency criterion operates identically regardless of the adjustment
method.  The closure returns `estimate` (a scalar RD, log-OR, etc.) and the
split loop checks:

```
both_halves_consistent = (est_half1 >= threshold) & (est_half2 >= threshold)
```

For the PS-covariate and AIPW estimators, `estimate` is the g-computed RD
rather than a model coefficient, which is the correct target for clinical
interpretation.

---

## 6. Weight Truncation

Extreme inverse-PS weights destabilise GLM fits in small subgroup halves.
The `ps_truncate` parameter (default 0.01) clips `π̂` to `[0.01, 0.99]`,
mirroring the `wt.truncate` parameter in `BR.TE.fit()`.  Larry's code
uses `x.truncate()` for quantile-based truncation; we default to
fixed-bound truncation as it is more predictable in small samples but
offer the quantile option via `ps_truncate_method = c("fixed", "quantile")`.

---

## 7. Updated `make_effect_estimator()` Signature

```r
make_effect_estimator <- function(
    outcome_type,
    treat.name,
    outcome.name,
    event.name       = NULL,
    offset.name      = NULL,
    effect_measure   = NULL,
    adjustment       = c("none", "regression", "ps_covariate", "aipw"),
    confounders.name = NULL,   # required if adjustment != "none"
    ps.name          = "inv_ps",
    pihat.name       = ".pihat",
    robust_se        = TRUE,
    ...
)
```

The `adjustment` parameter defaults to `"none"`.  When set to
`"ps_covariate"` or `"aipw"`, `confounders.name` is required.

For `"regression"`, the closure fits `Y ~ treat + confounders` and
g-computes the RD — no PS model is involved, but the confounders
improve efficiency.  This corresponds to the OR estimator in
`BR.TE.fit()`.

---

## 8. Relationship to Existing `BR.TE.fit()` Estimators

| `BR.TE.fit()` estimator | forestsearch `adjustment` | Description |
|---|---|---|
| `delta.naive` | `"none"` | Raw group means / unadjusted GLM coefficient |
| `delta.OR` | `"regression"` | Outcome regression with g-computation |
| `delta.DR` (v1) | `"ps_covariate"` | BR: outcome + inv.ps as covariate |
| `delta.DR3` (v3) | `"aipw"` | Funk & Davidian AIPW influence function |
| `delta.HT` | not implemented | Pure IPW (Horvitz-Thompson) — less stable |
| `delta.match` | not implemented | PS matching — different paradigm |

Versions 2, 4, 5, 6 from `BR.TE.fit()` (separate arm-specific PS terms,
GAM smoothing of PS) are available as future extensions but are not
included in the initial implementation.

---

## 9. ForestSearch-Specific Considerations

### 9.1 RCT vs. Observational Data

In a properly randomised trial, the true propensity score is known
(`π = 0.5` for 1:1 allocation, or the allocation ratio for stratified
randomisation).  Using the **known** propensity score is valid and avoids
estimating a nuisance model.  The `ps_known` argument allows the user to
supply the true π directly:

```r
forestsearch(
  ...,
  adjustment  = "ps_covariate",
  ps_known    = 0.5              # or a vector of known π per subject
)
```

When `ps_known` is supplied, no PS model is fitted; `inv_ps` is computed
directly from the known value.  For observational data, `ps_known = NULL`
(default) triggers PS estimation from `confounders.name`.

### 9.2 Bootstrap Bias Correction

The Forest Search bootstrap already re-runs the entire pipeline on each
bootstrap replicate.  If `adjustment = "ps_covariate"`, the PS model is
re-fitted within each bootstrap replicate (Option B), automatically
capturing the uncertainty in the PS model.  This uses Option A (pre-fitted)
for the consistency loop (speed) and Option B (re-fitted) for the bootstrap
(correctness).

### 9.3 Cross-Validation

The K-fold and N-fold CV paths also re-run the pipeline on each fold.
The PS model is re-fitted per fold, providing proper out-of-sample
calibration.

---

## 10. Implementation Phases (Updated)

| Phase | Additions for PS/DR methods |
|---|---|
| Phase 1 | `"regression"` adjustment only (no PS model needed) |
| Phase 2 | `"ps_covariate"` with pre-computed `inv_ps`; `ps_truncate`; threading `confounders.name` |
| Phase 3 | `"aipw"` influence-function estimator; `ps_known` for RCTs |
| Phase 4 | Bootstrap re-estimation of PS model within each replicate |
| Phase 5 | Validation: compare `"ps_covariate"` output to `BR.TE.fit()` on identical data |

Phase 1 is the natural first step because it exercises the g-computation
prediction machinery (predict under each arm, average) without any PS
modelling, and it directly corresponds to `delta.OR` in `BR.TE.fit()`.

---

## 11. References

- Bang H, Robins JM (2005). Doubly robust estimation in missing data and
  causal inference models. *Biometrics* 61:962–972.
- Funk MJ, Westreich D, Wiesen C, Stürmer T, Brookhart MA, Davidian M
  (2011). Doubly robust estimation of causal effects. *Am J Epidemiol*
  173:761–767.
- Robins JM, Rotnitzky A, Zhao LP (1994). Estimation of regression
  coefficients when some regressors are not always observed. *JASA*
  89:846–866.
- Tsiatis AA, Davidian M, Zhang M, Lu X (2008). Covariate adjustment for
  two-sample treatment comparisons in randomized clinical trials:
  a principled yet flexible approach. *Stat Med* 27:4658–4677.
- León LF et al. (2024). Exploratory subgroup identification in the
  heterogeneous Cox model. *Stat Med*. doi:10.1002/sim.10163.
