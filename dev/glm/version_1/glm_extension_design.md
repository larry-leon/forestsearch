# GLM Extension Design for `forestsearch`

**Document type:** Architecture design & implementation roadmap  
**Scope:** Extending `forestsearch` to binary and continuous (GLM) outcomes  
**Status:** Pre-implementation proposal  

---

## 1. Executive Summary

The `forestsearch` pipeline is composed of five logically separable stages. Of these,
**stages 2–4 (enumeration, evaluation, and bootstrap) are largely outcome-agnostic** at
their structural level — they operate on subgroup membership indicators and a scalar
effect estimate returned from a pluggable model function. Only stages 1 (candidate
factor selection) and 5 (visualization) are deeply survival-specific.

**Recommendation: extend within the current package** using a strategy-pattern
("pluggable estimator") architecture rather than creating a parallel codebase. This
avoids duplicating the combinatorial search engine, the bootstrap/CV infrastructure,
and the pkgdown site. The GLM layer adds two new functional areas:
`glm_effect_estimator` (internal) and `grf.subg.harm.glm` (exported), plus
`outcome_type` threading through the existing top-level functions.

The full implementation is moderately large (~500–800 lines of new R code) but
follows well-established patterns already present in the package.

---

## 2. Current Architecture: Stage-by-Stage Assessment

### Stage 1 — Candidate Factor Selection (`get_FSdata`, `grf.subg.harm.survival`)

| Component | Survival-specific? | Notes |
|---|---|---|
| LASSO via `glmnet` | **Yes** — uses `family = "cox"` | Needs `family = "binomial"` / `"gaussian"` branch |
| GRF via `causal_survival_forest()` | **Yes** | Needs `causal_forest()` branch for binary/continuous |
| Quartile/median binary cuts | No | Fully reusable |
| Factor enumeration | No | Fully reusable |

### Stage 2 — Subgroup Enumeration (`subgroup.search`, combinatorial helpers)

**Fully reusable.** The combinatorial enumeration operates on binary subgroup
membership indicators `J_1, ..., J_L` and is entirely independent of outcome type.
Minimum sample-size (`n.min`) and minimum-event (`d0.min`, `d1.min`) filters are the
only outcome-adjacent logic. For GLM outcomes, event-count filters are replaced by a
simple minimum-arm-size criterion.

### Stage 3 — Subgroup Evaluation / Consistency Criterion (`subgroup.consistency`)

**Structurally reusable, effect estimator must be plugged in.** The consistency loop:

```
for r in 1:R:
    split subgroup 50/50
    estimate effect in each half
    flag consistent if min(effect_1, effect_2) >= threshold
```

This loop is outcome-agnostic. What changes is the function that maps a data slice to
a scalar log-effect estimate. Currently this is an anonymous Cox call. The refactor
wraps this in an internal `effect_estimator_fn` argument.

### Stage 4 — Bootstrap Bias Correction (`forestsearch_bootstrap_dofuture`)

**Structurally reusable.** The IJ bootstrap procedure applies identically to any
scalar log-scale effect — log-HR, log-OR, log-RR, or mean difference. The correction:

    beta_corrected = beta_obs - mean(eta_b_selected + eta_b_obs)

does not assume a specific model. `outcome_type` threading is sufficient.

### Stage 5 — Visualization

**Not reusable.** Weighted KM curves, RMST displays, and survival-specific forest
plots all require `Surv` objects. A GLM extension needs its own visualization layer:
bar charts of response rates by arm, forest plots of odds ratios or mean differences,
and optionally risk difference tables.

---

## 3. Effect Measure Design Decisions

### 3.1 Binary Outcomes

Three effect measures are clinically relevant. Each maps to a GLM via a different link
function:

| Effect Measure | Link | `glm()` call | Scale for consistency |
|---|---|---|---|
| Odds Ratio (OR) | logit | `family = binomial(link = "logit")` | log-OR ≥ log(threshold) |
| Risk Ratio (RR) | log | `family = binomial(link = "log")` | log-RR ≥ log(threshold) |
| Risk Difference (RD) | identity | `family = binomial(link = "identity")` | RD ≥ threshold |

**Recommended default: log-OR** (logistic regression). Reasons:

1. Numerically most stable for binary outcomes; log-binomial (RR) frequently fails
   to converge when event rates approach 0 or 1 — a common situation in subgroups.
2. Odds ratios are the standard summary in most clinical trial regulatory submissions
   (Europe, FDA) for binary endpoints.
3. Modified Poisson (RR via Poisson with robust SE) is a well-known workaround for
   RR convergence issues and can be offered as an option.

**Convergence safeguard:** wrap each `glm()` call in `tryCatch`; if log-binomial
fails to converge, silently fall back to logistic and flag the subgroup.

**Consistency threshold mapping:** For OR, the consistency criterion becomes:

    both splits have log-OR >= log(or.threshold)

analogous to the Cox HR criterion. Default `or.threshold = 1.0` for consistency
(any excess odds), `or.screening.threshold = 1.5` for initial screening.

### 3.2 Continuous Outcomes

| Effect Measure | Model | Scale |
|---|---|---|
| Mean Difference (MD) | Linear regression (`lm`) | MD ≥ threshold |
| Standardised MD | Linear regression + SD scaling | SMD ≥ threshold |

**Recommended: raw mean difference** (unstandardised), with the user supplying a
clinically meaningful `md.threshold` for screening and `md.consistency = 0` for the
split criterion (analogous to `hr.consistency = 1.0` on the log scale).

For continuous outcomes the standard `lm()` coefficient for treatment is the within-
subgroup mean difference (experimental − control under 1:1 randomisation), which is
interpretable without further transformation.

---

## 4. GRF Candidate Selection for GLM Outcomes

The existing `grf.subg.harm.survival()` function uses `grf::causal_survival_forest()`
with an RMST-based tau criterion. The GLM analogue uses `grf::causal_forest()`:

```r
# Binary / continuous
cf <- grf::causal_forest(
  X = X_matrix,
  Y = outcome_vector,   # binary 0/1 or continuous numeric
  W = treat_vector,
  num.trees = num.trees,
  ...
)
tau_hat <- predict(cf)$predictions   # CATE estimates (risk diff scale for binary)
```

For binary outcomes `causal_forest()` estimates the individual treatment effect on the
**risk-difference scale** (ATE estimand). Policy trees (`policytree::policy_tree()`)
operate the same way regardless of outcome type, so the tree-based factor selection
logic from `grf.subg.harm.survival()` is directly portable.

A new exported function `grf.subg.harm.glm()` with arguments mirroring
`grf.subg.harm.survival()` but without `event.name` / RMST parameters covers this.

**LASSO** is straightforward: pass `family = "binomial"` or `family = "gaussian"` to
`glmnet::cv.glmnet()`. The existing LASSO selection code only needs a `family`
argument threaded through.

---

## 5. Implementation Architecture

### 5.1 New `outcome_type` Parameter

Thread `outcome_type = c("survival", "binary", "continuous")` through:

- `forestsearch()` — top-level dispatcher
- `subgroup.search()` — passes to effect estimator
- `subgroup.consistency()` — chooses effect estimator function
- `forestsearch_bootstrap_dofuture()` — passes to each bootstrap replicate
- `forestsearch_tenfold()` / `forestsearch_Kfold()` — same

For `outcome_type = "survival"` the behaviour is 100% unchanged (no regression).

### 5.2 Internal `effect_estimator_fn` Pattern

Introduce an internal function factory:

```r
#' @noRd
make_effect_estimator <- function(
    outcome_type,
    treat.name,
    outcome.name,
    event.name    = NULL,    # survival only
    effect_measure = NULL,   # "OR", "RR", "RD", "MD"
    ...
) {
  switch(outcome_type,
    survival   = .make_cox_estimator(treat.name, outcome.name, event.name, ...),
    binary     = .make_glm_binary_estimator(treat.name, outcome.name,
                                             effect_measure, ...),
    continuous = .make_lm_estimator(treat.name, outcome.name, ...)
  )
}
```

Each `.make_*_estimator()` returns a **closure** of the form:

```r
function(data_slice) {
  # Returns named list:
  #   estimate  : scalar log-effect (log-OR, log-HR, mean diff)
  #   se        : standard error
  #   converged : logical
  #   n0, n1    : arm sample sizes
}
```

This single interface means `subgroup.consistency()` and the bootstrap loop need
**no further modification** — they call `estimator_fn(data_slice)` and test
`estimate >= log(threshold)` (or `estimate >= threshold` for RD/MD).

### 5.3 New and Modified Files

#### New files

| File | Contents |
|---|---|
| `R/glm_effect_estimators.R` | `make_effect_estimator()`, `.make_glm_binary_estimator()`, `.make_lm_estimator()`, convergence fallback logic |
| `R/grf_subg_harm_glm.R` | `grf.subg.harm.glm()` — causal forest candidate selection for binary/continuous |
| `R/glm_simulate.R` | `simulate_binary_trial()`, `simulate_continuous_trial()` — GLM data-generating mechanisms for testing and vignettes |
| `R/plot_glm_results.R` | `plot_subgroup_binary()`, `plot_subgroup_continuous()` — bar/forest plots for GLM summaries |
| `R/summarize_glm_bootstrap.R` | `summarize_glm_bootstrap_results()` — `gt` tables for OR/RR/RD/MD estimates |

#### Modified files (minimal, additive changes only)

| File | Change |
|---|---|
| `R/forestsearch_main.R` | Add `outcome_type`, `effect_measure` args; dispatch to `make_effect_estimator()` |
| `R/subgroup_search.R` | Thread `estimator_fn` argument; replace inline Cox call |
| `R/subgroup_consistency.R` | Accept `estimator_fn`; replace inline Cox call |
| `R/get_FSdata_helpers.R` | Add `family` argument to LASSO section |
| `R/bootstrap_dofuture_main.R` | Thread `outcome_type` + `estimator_fn` into worker |
| `R/forestsearch_cross_validation.R` | Same as bootstrap |
| `DESCRIPTION` | No new hard dependencies needed (all GLM/lm functions are base R or already-imported `glmnet`/`grf`); add `Suggests: sandwich` for robust SE in RR via modified Poisson |

---

## 6. Consistency Criterion Mapping by Outcome Type

| Outcome | Effect scale | Screening threshold arg | Consistency threshold arg | Default screen | Default consistency |
|---|---|---|---|---|---|
| Survival | log-HR | `hr.threshold` | `hr.consistency` | 1.25 | 1.0 |
| Binary (OR) | log-OR | `or.threshold` | `or.consistency` | 1.5 | 1.0 |
| Binary (RR) | log-RR | `rr.threshold` | `rr.consistency` | 1.5 | 1.0 |
| Binary (RD) | RD | `rd.threshold` | `rd.consistency` | 0.05 | 0.0 |
| Continuous | MD | `md.threshold` | `md.consistency` | — | 0.0 |

The split-consistency criterion remains structurally identical:

    p_consistency = mean(min(effect_1r, effect_2r) >= consistency_threshold, r = 1..R)

For continuous outcomes, `md.threshold = 0` for consistency is equivalent to "both
halves show a positive treatment effect" — the analogue of HR > 1.0 for harm.

---

## 7. Minimum Sample Size Criteria

For GLM outcomes the event-count filters (`d0.min`, `d1.min`) are replaced by
arm-size filters only:

- `n0.min`: minimum control arm subjects per subgroup (default 30)
- `n1.min`: minimum experimental arm subjects per subgroup (default 30)

A warning should be issued if response rates are very low (< 5%) for binary
outcomes and the log-binomial link is selected, advising the user to consider
`effect_measure = "OR"`.

---

## 8. DESCRIPTION Changes

```
# Add to Suggests (not Imports — GLM uses base R):
sandwich,       # robust SE for modified Poisson (RR estimation)

# Update Title / Description to mention GLM outcomes
Title: Exploratory Subgroup Identification in Clinical Trials
Description: Implements the Forest Search methodology for exploratory subgroup
    identification in randomised controlled trials. Supports survival endpoints
    (Cox proportional hazards, AFT) and generalised linear model outcomes
    (binary: odds ratio, risk ratio, risk difference; continuous: mean difference).
    ...
```

No new hard `Imports` are needed because:
- `glm()` and `lm()` are base R
- `glmnet` is already imported (LASSO)
- `grf::causal_forest` is already imported

---

## 9. Key Design Principles

### 9.1 Backward Compatibility

`outcome_type` defaults to `"survival"` and `effect_measure` defaults to `NULL`
(resolved to `"HR"` internally). Every existing call to `forestsearch()` works
unchanged with no deprecation.

### 9.2 Convergence Robustness

Log-binomial GLMs frequently fail to converge, especially in small subgroups. The
estimator closure must:

1. Attempt `glm(..., family = binomial("log"))` first.
2. On `warning("did not converge")` or `tryCatch` error, fall back to
   `glm(..., family = binomial("logit"))` and compute the approximate RR from the OR
   via the formula `RR ≈ OR / ((1 - p0) + p0 * OR)` where `p0` is the control-arm
   response rate.
3. Record `converged = FALSE` in the return list and propagate a flag so the final
   summary can note which subgroups used the fallback.

### 9.3 Bootstrap with GLM

The IJ bootstrap loop calls `estimator_fn(data_slice)` in each replicate — identical
to the survival case. The `seed` argument already exposed in
`forestsearch_bootstrap_dofuture()` continues to control RNG. No structural change to
the bootstrap worker is required.

### 9.4 Tidyverse Style Compliance

All new functions follow the existing package conventions:
- `snake_case` function and argument names
- `@param`, `@return`, `@examples` roxygen blocks for all exported functions
- `\dontrun{}` for examples requiring substantial computation
- No `library()` calls in function bodies
- `data.table` for any in-function data manipulation consistent with existing code

---

## 10. Suggested Implementation Order

| Phase | Deliverables | Effort |
|---|---|---|
| **Phase 1** | Refactor `subgroup.consistency()` to accept `estimator_fn`; implement `make_effect_estimator()` and `.make_glm_binary_estimator()` for OR; unit tests with simulated 2-arm binary data | ~2 days |
| **Phase 2** | Thread `outcome_type` through `forestsearch()`, bootstrap, and CV; implement `grf.subg.harm.glm()`; LASSO `family` branching | ~2 days |
| **Phase 3** | Continuous outcome estimator (`lm`); RR / RD estimators with convergence fallback; `summarize_glm_bootstrap_results()` | ~1 day |
| **Phase 4** | Visualization layer (`plot_subgroup_binary()`, `plot_subgroup_continuous()`); simulation DGMs for testing; vignette | ~2 days |
| **Phase 5** | Full `devtools::check(args = "--as-cran")` pass; win-builder | ~0.5 day |

---

## 11. Illustrative User-Facing API (Target State)

```r
# --- Binary outcome (logistic / OR) ---
fs_bin <- forestsearch(
  df.analysis       = trial_data,
  outcome.name      = "response",       # binary 0/1
  treat.name        = "treat",
  confounders.name  = c("age", "biomarker", "stage"),
  outcome_type      = "binary",
  effect_measure    = "OR",             # "OR" | "RR" | "RD"
  or.threshold      = 1.5,              # screening: OR >= 1.5
  or.consistency    = 1.0,             # split consistency: OR >= 1.0
  pconsistency.threshold = 0.90,
  fs.splits         = 400,
  use_lasso         = TRUE,
  use_grf           = TRUE
)

# --- Continuous outcome (mean difference) ---
fs_cont <- forestsearch(
  df.analysis       = trial_data,
  outcome.name      = "biomarker_change",
  treat.name        = "treat",
  confounders.name  = c("age", "baseline_score", "region"),
  outcome_type      = "continuous",
  md.threshold      = 2.0,             # screening: MD >= 2.0 units
  md.consistency    = 0.0,             # split consistency: MD >= 0
  pconsistency.threshold = 0.90,
  fs.splits         = 400
)

# Bootstrap and CV work identically
fs_bc <- forestsearch_bootstrap_dofuture(fs_bin, nb_boots = 1000)
summaries <- summarize_glm_bootstrap_results(fs_bc)
```

---

## 12. Alternative: Separate Package (`forestsearch.glm`)

For completeness, the case for a separate package:

**Pros:**
- Clean separation of survival vs. GLM dependencies in DESCRIPTION
- Independent version lifecycle
- No risk of breaking existing survival API

**Cons:**
- Duplicates the entire combinatorial search engine (~5,000 lines)
- Duplicates the bootstrap/CV infrastructure
- Two pkgdown sites, two CRAN submissions, two maintenance burdens
- Users must install two packages for mixed-endpoint trials (e.g., OS + ORR)

**Conclusion:** The separation cost is high and the benefit is low. The shared
infrastructure is the valuable part of `forestsearch`; it should not be duplicated.
The strategy-pattern extension keeps the codebase unified and the user API coherent.

---

## 13. References

- León et al. (2024). *Exploratory subgroup identification in the heterogeneous Cox model.* Statistics in Medicine. https://doi.org/10.1002/sim.10163  
- Tibshirani et al. (2023). grf: Generalized Random Forests. R package. https://grf-labs.github.io/grf/  
- Shirvaikar & Dominici (2024). *Targeting relative risk heterogeneity with causal forests.* arXiv:2309.15793  
- Zou G (2004). *A modified Poisson regression approach to prospective studies with binary data.* Am J Epidemiol 159(7):702–6.  
- Ballarini et al. (2021). *Subgroup identification in clinical trials: an overview of available methods.* Biometrical Journal.  
- Dandl et al. (2024). *Model-based forests for subgroup identification.* Statistics in Medicine.
