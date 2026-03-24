# GLM Extension Design for `forestsearch` — Version 2

**Document type:** Architecture design & implementation roadmap
**Scope:** Extending `forestsearch` to binary, continuous, and count/rate outcomes
**Status:** Pre-implementation proposal (revision of v1)
**Date:** 2026-03-23

---

## 1. What Changed from Version 1

Version 1 proposed OR (logistic regression) as the default binary effect measure.
This revision makes three substantive changes:

1. **Risk difference (RD) is now the primary estimand for binary outcomes.**
   The ICH E9(R1) addendum emphasises the estimand framework, and the risk
   difference is the most directly interpretable treatment effect measure for
   binary endpoints — it answers "how many additional patients per 100 are
   harmed by treatment in this subgroup?" This aligns with the clinical
   question that motivates Forest Search.

2. **Count/rate outcomes via Poisson regression with log-time offset** are added
   as a new capability.  Specifically, `glm(event ~ treat, family = poisson,
   offset = log(time))` approximates the Cox model when the baseline hazard is
   locally constant within each subject's observation window.  This provides a
   **built-in sanity check**: running the Forest Search pipeline with
   `effect_measure = "IRR"` (incidence rate ratio) on survival data should
   reproduce the subgroup identification and effect estimates from the existing
   `outcome_type = "survival"` path.

3. **The effect-measure taxonomy is expanded.** Binary outcomes now support five
   measures: RD, OR, RR, IRR, and IRD — all under `outcome_type = "binary"`.
   The rate-based measures (IRR, IRD) require an `offset.name` argument
   pointing to the follow-up time column.

---

## 2. Complete Effect-Measure Taxonomy

### 2.1 Binary Outcomes (`outcome_type = "binary"`)

| Effect Measure | Link | Model | Scale | `offset.name` required? |
|---|---|---|---|---|
| **RD** (risk difference) | identity | `binomial(link = "identity")` | identity (e.g., 0.08) | No |
| OR (odds ratio) | logit | `binomial(link = "logit")` | log-OR | No |
| RR (risk ratio) | log | `binomial(link = "log")` / modified Poisson | log-RR | No |
| IRR (incidence rate ratio) | log | `poisson(link = "log")` | log-IRR | **Yes** |
| IRD (incidence rate difference) | — | Computed from Poisson-predicted rates | identity | **Yes** |

**Default:** `effect_measure = "RD"`.

### 2.2 Continuous Outcomes (`outcome_type = "continuous"`)

| Effect Measure | Model | Scale |
|---|---|---|
| **MD** (mean difference) | `lm()` | identity |

### 2.3 Survival Outcomes (`outcome_type = "survival"`)

Unchanged from v0.1.0; the Cox pipeline is not touched.

---

## 3. The Poisson-Offset Bridge to Cox

### 3.1 Theoretical Basis

For subject $i$ with binary event indicator $D_i$ and follow-up time $T_i$,
the Poisson regression model with log-time offset is:

$$\log E(D_i \mid A_i) = \log(T_i) + \beta_0 + \beta_1 A_i$$

which implies $E(D_i) = T_i \exp(\beta_0 + \beta_1 A_i)$, so $\exp(\beta_1)$
is the **incidence rate ratio** (IRR).  Under the exponential hazard assumption
(constant baseline hazard), this equals the Cox hazard ratio.  More generally,
with Weibull or other parametric hazards the approximation improves as sample
size increases and the model includes covariates that capture the shape of the
hazard function.

### 3.2 Role as Sanity Check

Given a survival dataset with columns `tte`, `event`, `treat`, and covariates:

```r
# Existing survival path
fs_surv <- forestsearch(
  df.analysis  = trial_data,
  outcome.name = "tte",
  event.name   = "event",
  treat.name   = "treat",
  outcome_type = "survival",
  ...
)

# Poisson-offset path (should approximate the above)
fs_pois <- forestsearch(
  df.analysis   = trial_data,
  outcome.name  = "event",
  treat.name    = "treat",
  outcome_type  = "binary",
  effect_measure = "IRR",
  offset.name   = "tte",
  ...
)

# Compare: fs_surv and fs_pois should identify similar subgroups
# and produce similar log-HR / log-IRR estimates
```

This is not expected to give *identical* results — the Cox model is
semiparametric and the Poisson model is parametric — but the subgroup
identification and effect estimates should agree closely, especially in
larger samples and when the true hazard is approximately constant.

### 3.3 When the Approximation Breaks Down

The Poisson-offset approximation is known to diverge from the Cox partial
likelihood when: the hazard function is strongly time-varying (non-PH is
extreme), individual follow-up times are very heterogeneous (heavy
administrative censoring at different times), or event rates are very high
(> 50%) such that the rare-event Poisson assumption fails.  These are
precisely the settings where the user should prefer the Cox path.

---

## 4. Estimator Closure Design

### 4.1 Updated `make_effect_estimator()` Signature

```r
make_effect_estimator <- function(
    outcome_type,
    treat.name,
    outcome.name,
    event.name     = NULL,    # survival only
    offset.name    = NULL,    # IRR / IRD only (log of this col used as offset)
    effect_measure = NULL,    # NULL -> default per outcome_type
    ...
)
```

### 4.2 New Estimator: `.make_poisson_rate_estimator()`

Returns a closure `function(data_slice)` that fits:

```r
glm(
  outcome ~ treat,
  family = poisson(link = "log"),
  offset = log(time),
  data   = data_slice
)
```

The treatment coefficient is the log-IRR.  SE is model-based by default;
if `sandwich` is available and `robust_se = TRUE`, the robust sandwich SE
is used (recommended for the modified-Poisson interpretation where the
Poisson variance assumption is known to be wrong for binary data).

### 4.3 Incidence Rate Difference (IRD)

IRD is computed from Poisson-predicted rates rather than via an identity-link
Poisson model (which has severe convergence problems).  The procedure:

1. Fit `glm(event ~ treat, family = poisson(link = "log"), offset = log(time))`
2. Predict rate for each arm: $\hat\lambda_a = \exp(\hat\beta_0 + \hat\beta_1 a) / \bar{T}_a$
3. $\widehat{IRD} = \hat\lambda_1 - \hat\lambda_0$
4. SE via delta method from the Poisson `vcov`

This avoids the convergence problems of the identity-link Poisson while still
providing a difference-scale measure suitable for the consistency criterion.

---

## 5. Consistency Criterion Mapping (Updated)

| Outcome | Measure | Screening threshold arg | Consistency threshold arg | Default screen | Default consistency |
|---|---|---|---|---|---|
| Survival | HR | `hr.threshold` | `hr.consistency` | 1.25 | 1.0 |
| Binary | **RD** | `rd.threshold` | `rd.consistency` | 0.05 | 0.0 |
| Binary | OR | `or.threshold` | `or.consistency` | 1.5 | 1.0 |
| Binary | RR | `rr.threshold` | `rr.consistency` | 1.5 | 1.0 |
| Binary | IRR | `irr.threshold` | `irr.consistency` | 1.25 | 1.0 |
| Binary | IRD | `ird.threshold` | `ird.consistency` | 0.01 | 0.0 |
| Continuous | MD | `md.threshold` | `md.consistency` | — | 0.0 |

For log-scale measures (HR, OR, RR, IRR) the consistency criterion compares
`log(estimate) >= log(threshold)`.  For identity-scale measures (RD, IRD, MD)
it compares `estimate >= threshold` directly.

The **RD consistency default of 0.0** means "both random halves show any
positive risk difference (treatment worse than control)" — the direct analogue
of HR > 1.0 in the survival setting.

---

## 6. Risk Difference Implementation Notes

### 6.1 Why RD as Default

In the context of Forest Search's subgroup identification for harm, the risk
difference answers the most actionable clinical question: "What is the absolute
excess event probability attributable to treatment in this subgroup?"  A finding
of RD = 0.12 in a subgroup means 12 additional events per 100 patients treated —
directly communicable to clinicians and regulators.

By contrast, an OR of 2.0 in a subgroup with 5% baseline risk (RD ≈ 0.05) and
an OR of 2.0 in a subgroup with 40% baseline risk (RD ≈ 0.16) have very different
clinical implications despite the identical odds ratio.

### 6.2 Convergence Strategy for RD

The identity-link binomial is notoriously fragile.  The estimator closure uses a
three-tier fallback:

1. **Identity-link binomial** with starting values `c(p0, p1 - p0)`.
2. **Margins from logistic regression:** Fit logistic, predict probabilities for
   each arm, take the difference.  SE via delta method from logistic `vcov`.
3. **Raw means fallback:** `p1 - p0` with delta-method SE
   `sqrt(p1(1-p1)/n1 + p0(1-p0)/n0)`.

The fallback is flagged in the return list (`converged = FALSE`,
`method_used = "logistic_margins"` or `"raw_means"`) so that summaries can
report which subgroups required the approximation.

---

## 7. GRF Candidate Selection

`grf::causal_forest()` estimates CATEs on the **risk-difference scale** for
binary outcomes and the **mean-difference scale** for continuous outcomes.
This directly aligns with the RD-first philosophy.  No change to
`grf.subg.harm.glm()` is needed.

For rate-based analyses (IRR/IRD), GRF is still run on the binary event
indicator without the time offset, since `causal_forest()` does not accept
offsets.  The candidate factors from GRF are then passed to the Poisson
consistency loop, which does incorporate the offset.

---

## 8. LASSO Candidate Selection

The LASSO path in `get_FSdata_helpers.R` currently uses
`glmnet(family = "cox")`.  The branching by `outcome_type`:

| `outcome_type` | `glmnet` family | Outcome vector |
|---|---|---|
| `"survival"` | `"cox"` | `Surv(time, event)` |
| `"binary"` (RD, OR, RR) | `"binomial"` | event indicator (0/1) |
| `"binary"` (IRR, IRD) | `"poisson"` | event indicator; `offset = log(time)` |
| `"continuous"` | `"gaussian"` | continuous outcome |

`glmnet::cv.glmnet()` supports `offset` and all these families natively.

---

## 9. Simulation DGMs for Testing

### 9.1 Binary Trial with Known Subgroup (`simulate_binary_trial()`)

Already in v1.  Unchanged except `in_subgroup` is now `TRUE`/`FALSE` logical
with associated `true_rd` (risk difference) recorded in the output.

### 9.2 Continuous Trial (`simulate_continuous_trial()`)

Unchanged from v1.

### 9.3 **NEW:** Rate-Based Trial (`simulate_rate_trial()`)

Generates survival-like data with binary event indicator and follow-up time,
suitable for both Cox and Poisson-offset analysis.  Parameters control the true
incidence rate ratio in the subgroup.  The DGM is:

```
h_i(t) = lambda_0 * exp(beta_treat * A_i + beta_sg * A_i * 1(i in H) + Z_i' gamma)
T_i ~ Exp(h_i)          (exponential survival)
C_i ~ Unif(c_lo, c_hi)  (administrative censoring)
event_i = 1(T_i <= C_i)
tte_i   = min(T_i, C_i)
```

Under this DGM, `exp(beta_treat)` is the ITT hazard ratio (≈ IRR) and
`exp(beta_treat + beta_sg)` is the subgroup-specific HR (≈ IRR).  The
exponential baseline ensures the Poisson-offset approximation is exact,
providing the cleanest sanity check.

---

## 10. Updated Implementation Phases

| Phase | Deliverables | Effort |
|---|---|---|
| **Phase 1** | Refactor `make_effect_estimator()` with all 7 estimators (HR, RD, OR, RR, IRR, IRD, MD); unit smoke tests on simulated data | ~2 days |
| **Phase 2** | Thread `outcome_type` + `effect_measure` + `offset.name` through `forestsearch()`, `subgroup.consistency()`, `subgroup.search()`, bootstrap, CV | ~2 days |
| **Phase 3** | `grf.subg.harm.glm()`; LASSO `family` branching; `simulate_rate_trial()`; Cox vs. Poisson sanity check script | ~1.5 days |
| **Phase 4** | Visualization: `plot_subgroup_binary()` (RD bar charts, forest plots); `plot_subgroup_rates()` (rate-based comparisons) | ~1.5 days |
| **Phase 5** | Vignette; `devtools::check(args = "--as-cran")`; win-builder | ~1 day |

---

## 11. Target User-Facing API

```r
# ── Binary outcome, risk difference (primary use case) ───────────────────────
fs_bin <- forestsearch(
  df.analysis    = trial_data,
  outcome.name   = "response",
  treat.name     = "treat",
  confounders.name = c("age", "biomarker", "stage"),
  outcome_type   = "binary",
  effect_measure = "RD",             # default for binary

  rd.threshold   = 0.05,             # screen: RD >= 5 percentage points
  rd.consistency = 0.0,              # split consistency: RD >= 0
  pconsistency.threshold = 0.90,
  use_grf        = TRUE
)

# ── Poisson-offset sanity check against Cox ──────────────────────────────────
fs_rate <- forestsearch(
  df.analysis    = trial_data,
  outcome.name   = "event",
  treat.name     = "treat",
  confounders.name = c("age", "biomarker", "stage"),
  outcome_type   = "binary",
  effect_measure = "IRR",
  offset.name    = "tte",            # log(tte) used as Poisson offset
  irr.threshold  = 1.25,
  irr.consistency = 1.0,
  pconsistency.threshold = 0.90,
  use_grf        = TRUE
)

# ── Compare with existing survival pipeline ──────────────────────────────────
fs_surv <- forestsearch(
  df.analysis    = trial_data,
  outcome.name   = "tte",
  event.name     = "event",
  treat.name     = "treat",
  outcome_type   = "survival",       # existing pipeline, unchanged
  hr.threshold   = 1.25,
  hr.consistency = 1.0,
  pconsistency.threshold = 0.90,
  use_grf        = TRUE
)
```

---

## 12. DESCRIPTION Changes

```
Suggests:
    sandwich        # robust SE for modified Poisson (RR) and rate models (IRR)
```

Updated Title and Description to mention binary and rate-based outcomes:

```
Title: Exploratory Subgroup Identification in Clinical Trials
Description: Implements the Forest Search methodology for exploratory subgroup
    identification in randomised controlled trials. Supports survival endpoints
    (Cox proportional hazards, AFT), binary outcomes (risk difference, odds
    ratio, risk ratio), and rate-based analyses (Poisson regression with
    exposure-time offset for incidence rate ratios) with bootstrap bias
    correction. The Poisson-offset path approximates the Cox model, providing
    a built-in cross-validation between the survival and GLM pipelines.
    ...
```

No new hard `Imports` required — `glm()`, `lm()`, `poisson()`, `binomial()`
are base R; `glmnet` and `grf` are already imported.

---

## 13. Key Design Principles (Unchanged from v1, Reinforced)

- **Backward compatibility:** `outcome_type` defaults to `"survival"`;
  every existing call works unchanged.
- **Single codebase:** No package duplication.
- **Convergence robustness:** Every GLM estimator has a graceful fallback.
- **Tidyverse style:** `snake_case`, roxygen2, CRAN-compliant examples.
- **Bootstrap agnosticism:** The IJ bootstrap operates on a scalar log-scale
  (or identity-scale) estimate; it does not need to know the model.

---

## 14. References

- León et al. (2024). *Statistics in Medicine.* doi:10.1002/sim.10163
- Zou G (2004). *A modified Poisson regression approach.* Am J Epidemiol 159(7):702–6.
- Laird N, Olivier D (1981). *Covariance analysis of censored survival data
  using log-linear analysis techniques.* JASA 76:231–240.
- Whitehead J (1980). *Fitting Cox's regression model to survival data using
  GLIM.* Applied Statistics 29:268–275.
- ICH E9(R1) (2019). *Addendum on estimands and sensitivity analysis.*
