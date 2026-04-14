# Deep-Dive Transparency Review: `actg175_simulations_2.qmd`

**Reviewer:** Claude (Anthropic)
**Date:** 2026-04-04
**Scope:** Transparency, adoptability, and clarity for users adapting the
ForestSearch simulation framework to their own trial data.

---

## Executive Summary

The document does a commendable job of centralising simulation parameters
in the setup chunk, providing a side-by-side GBSG-vs-ACTG175 comparison,
and walking through the `generate_aft_dgm_flex()` pipeline step by step.
However, several aspects reduce its effectiveness as a *template* that
other users can confidently adapt to their own datasets.  The issues
cluster into five themes:

1. **Hidden "why" behind parameter choices** — many values are stated but
   not justified
2. **Unexplained bridging between data-prep and DGM** — column renaming
   and factor coding create silent coupling
3. **The adjusted Cox formula is a silent contract** — `cox_formula_adj`
   must match the DGM but nothing enforces or documents this
4. **Calibration workflow is ad hoc** — the manual grid search works but
   disguises what `calibrate_k_inter()` would do if it supported non-GBSG
   data
5. **Key ForestSearch parameter choices lack guidance** — thresholds are
   stated but the document doesn't explain when a user should change them

Below are the findings organised by document section, each with a
severity tag (**Critical** / **High** / **Medium** / **Low**) and a
suggested fix.

---

## 1. Setup Chunk (Lines 22–94)

### 1.1 Cutpoint rationale is missing [HIGH]

The cutpoints `actg_age_cut = 34`, `actg_wtkg_cut = 75`, and
`actg_cd80_cut = 1292` are stated as constants but their clinical or
data-driven origin is only hinted at in a single sentence (line 148):

> "The subgroup definition for ACTG175 is derived from the GRF screening
> analysis in `actg175_analysis_v1.Rnw`"

A user adapting this to their own trial has no way to replicate this step
because `actg175_analysis_v1.Rnw` is not included or described.

**Suggested fix:** Add a subsection (e.g. §1.3 "Choosing Discretisation
Cutpoints") that explicitly describes the process:

```markdown
### Choosing Discretisation Cutpoints

ForestSearch operates on binary (factor) covariates.  Continuous
variables must be discretised before building the DGM.  There are
three common strategies:

1. **Clinical cutpoints** — use thresholds defined by domain knowledge
   (e.g., CD4 < 200 for AIDS-defining immunosuppression).
2. **Data-driven medians** — split at the sample median; appropriate
   when no clinical threshold exists.  Used here for z3 (Karnofsky)
   and z4 (CD4 baseline).
3. **GRF-guided splits** — run a preliminary GRF analysis on the
   observed data to identify the variable-importance-weighted split
   points that maximise treatment effect heterogeneity.  This was
   the approach for z1 (age ≤ 34) and z5 (CD8 > 1292).

For your own dataset, we recommend starting with (1) where available,
falling back to (2), and optionally using (3) as a sensitivity check.
```

### 1.2 `dgm_subgroup_cuts` semantics are unclear [HIGH]

The line:

```r
dgm_subgroup_cuts <- list(z1 = 1L, z5 = 1L)
```

…means "subjects with `z1 == 1` AND `z5 == 1` are in the harm
subgroup."  But the user must know that `z1 = 1` means "age ≤ 34" (not
"age > 34"), which depends on the factor coding done 50 lines later.
This coupling is invisible.

**Suggested fix:** Add an inline annotation:

```r
# Within-level values that define membership in H.
# z1 = 1 means age <= actg_age_cut  (young)
# z5 = 1 means CD8 > actg_cd80_cut  (high CD8)
dgm_subgroup_cuts <- list(z1 = 1L, z5 = 1L)
```

And in the prose, add: "The values in `dgm_subgroup_cuts` must match the
factor coding in the data preparation chunk — `1L` means the subject
falls into the 'at-risk' level of each variable."

### 1.3 `cox_formula_adj` uses only 3 of 12 confounders [MEDIUM]

```r
cox_formula_adj <- survival::Surv(y_sim, event_sim) ~
  treat_sim + z1 + z2 + z3
```

The DGM has 12 binary confounders (`z1`–`z12`), but the adjusted Cox
formula uses only `z1 + z2 + z3`.  There is no explanation of:

- *Why* only three?  Is this by design (parsimony, avoiding
  overadjustment) or an oversight?
- Which confounders were selected and how?
- The consequence of this mismatch for the adjusted ITT HR.

This matters because `run_simulation_analysis()` uses this formula to
compute `hr.adj.itt`.  If a user copies this pattern and includes all
12 confounders, they may get different results without understanding why.

**Suggested fix:** Add a comment block:

```r
# Adjusted Cox formula for the simulation.
# We include only the subgroup-defining variables (z1) and strong
# prognostic confounders (z2, z3) to mimic a parsimonious clinical
# analysis.  The full set z1–z12 could be used for a saturated model
# but may over-adjust in small simulated samples.
#
# IMPORTANT: This formula must reference the simulated column names
# (y_sim, event_sim, treat_sim) and confounder names that exist in
# the super-population (z1, z2, ...).
```

### 1.4 `sim_analysis_time = 1200` needs time-scale anchoring [MEDIUM]

The value 1200 days (~3.3 years) is stated but not motivated.  Users
adapting to their own dataset need to know: *how should I choose
`analysis_time`?*

**Suggested fix:** Add guidance:

```r
# Administrative censoring time on the DGM's time scale.
# Must use the SAME units as the outcome variable in the source data.
# ACTG175 uses days, so analysis_time is in days.
# Rule of thumb: set to the observed median follow-up time,
# or the planned analysis calendar time of the trial.
sim_analysis_time <- 1200   # ~3.3 years, close to ACTG175 median follow-up
```

### 1.5 `fs_vi_grf_min = -0.2` is unexplained [LOW]

A negative GRF variable-importance threshold is unusual — it means all
variables pass the screen.  If this is intentional (disable GRF
screening), say so.  If it's a permissive default, explain when a user
should tighten it.


---

## 2. Data Preparation (Lines 178–226)

### 2.1 Factor coding direction is not self-documenting [HIGH]

The factor coding block (lines 198–219) creates `z1`–`z12` as factors,
but a reader must mentally track which level is `1` and which is `0`.
For example:

```r
actg_df$z5 <- as.factor(ifelse(actg_df$cd80 > actg_cd80_cut, 1L, 0L))
```

This means `z5 = 1` → "high CD8", which in the harm subgroup means
treatment is harmful for patients with high CD8.  But `z1 = 1` → "young
age", so the directionality is inconsistent (sometimes "1" = high,
sometimes "1" = low).

**Suggested fix:** Add a summary comment block after the factor creation:

```r
# ── Factor coding summary ──────────────────────────────────────
# z1  = 1 ⇒ age ≤ 34       (young)
# z2  = 1 ⇒ weight ≤ 75 kg (low weight)
# z3  = 1 ⇒ Karnofsky ≤ median (lower performance)
# z4  = 1 ⇒ CD4 ≤ median   (lower baseline CD4)
# z5  = 1 ⇒ CD8 > 1292     (high CD8)
# z6  = 1 ⇒ hemophilia = yes
# z7  = 1 ⇒ homosexual activity = yes
# z8  = 1 ⇒ IV drug use = yes
# z9  = 1 ⇒ non-white
# z10 = 1 ⇒ male
# z11 = 1 ⇒ prior antiretroviral = yes
# z12 = 1 ⇒ symptomatic = yes
```

### 2.2 `flag_harm` construction uses programmatic lookup but is hard to verify [MEDIUM]

Lines 222–225 build `flag_harm` from `dgm_subgroup_vars` and
`dgm_subgroup_cuts` using `[[` indexing:

```r
actg_df$flag_harm <- as.integer(
  actg_df[[dgm_subgroup_vars[1]]] == dgm_subgroup_cuts[[dgm_subgroup_vars[1]]] &
  actg_df[[dgm_subgroup_vars[2]]] == dgm_subgroup_cuts[[dgm_subgroup_vars[2]]]
)
```

This is general-purpose but opaque.  The reader cannot see at a glance
what the subgroup actually is.  It also silently compares a *factor*
column against an *integer* value (`1L`), which works because R coerces,
but hides a type mismatch.

**Suggested fix:** Add an explicit assertion and human-readable echo:

```r
# Verify subgroup definition is as expected
stopifnot(
  identical(dgm_subgroup_vars, c("z1", "z5")),
  identical(dgm_subgroup_cuts, list(z1 = 1L, z5 = 1L))
)
cat("Harm subgroup H: z1 == 1 (age ≤", actg_age_cut,
    ") AND z5 == 1 (CD8 >", actg_cd80_cut, ")\n")
```


---

## 3. DGM Construction (Lines 458–663)

### 3.1 No time-scale consistency check [CRITICAL]

The `extreme_subgroups` vignette has an excellent warning about
time-scale mismatch:

> "The single most common simulation failure is building the DGM on days
> and then passing `analysis_time` in months."

This document uses `dfa$y <- dfa$time_days` (days) and
`sim_analysis_time = 1200` (days), which is correct.  But there is no
diagnostic check to *verify* consistency.

**Suggested fix:** Add a diagnostic block after creating `dgm_alt`:

```r
# ── Time-scale sanity check ────────────────────────────────────
# exp(mu) is the implied median event time on the DGM time scale.
# It should be comparable to the observed median event time.
observed_median <- median(dfa$y[dfa$event == 1])
dgm_median      <- exp(dgm_alt$model_params$mu)

cat(sprintf("Observed median event time: %.0f days\n", observed_median))
cat(sprintf("DGM implied median:         %.0f days\n", dgm_median))
cat(sprintf("analysis_time:              %.0f days\n", sim_analysis_time))

if (sim_analysis_time < dgm_median * 0.5)
  warning("analysis_time is much shorter than DGM median — ",
          "expect heavy administrative censoring!")
```

### 3.2 `compute_dgm_cde()` is called without explanation [HIGH]

Line 500:
```r
dgm_alt <- compute_dgm_cde(dgm_alt)
```

This enriches the DGM with controlled direct effect (CDE) metrics, but
the document never explains what CDEs are, why they matter, or what
fields they add.  A user adapting this template might omit this call and
then get errors downstream when `build_estimation_table()` tries to
access CDE columns.

**Suggested fix:** Add a brief paragraph:

```markdown
`compute_dgm_cde()` adds individual-level potential-outcome quantities
(`theta_0`, `theta_1`, `loghr_po`) to the super-population.  These are
used by `build_estimation_table()` and `interpret_estimation_table()` to
compute CDE bias metrics.  Always call this after creating any
alternative-hypothesis DGM.
```

### 3.3 Calibration explains *what* but not *why* [HIGH]

The calibration section (lines 515–569) does a grid search over
`k_inter` to hit `cal_target_hr = 2.0`.  But the document never
explains:

- **Why HR(H) = 2.0?**  Is this a conventional choice?  A power
  calculation target?
- **What does `k_inter` control?**  It multiplies the treatment-by-
  subgroup interaction in the AFT linear predictor.  The connection
  between `k_inter` and the resulting HR(H) is non-linear and
  dataset-dependent.
- **Why is the grid search necessary?**  The existing
  `calibrate_k_inter()` wraps `setup_gbsg_dgm()`, so it doesn't work
  for non-GBSG data.  This gap should be flagged as a known limitation
  and future enhancement.

**Suggested fix:** Add a prose paragraph before the chunk:

```markdown
### Why calibrate?

The DGM's `k_inter` parameter controls the strength of the
treatment-by-subgroup interaction on the AFT (log-time) scale.
Because the mapping from `k_inter` to the Cox-model hazard ratio
HR(H) is non-linear and depends on the specific covariate
distribution, we cannot set `k_inter` analytically to achieve a
target HR.  Instead, we search over a grid of candidate values.

A target HR(H) of 2.0 is used here because it represents a clinically
meaningful level of harm (treatment doubles the hazard in the
subgroup) and matches the simulation scenarios in León et al. (2024).
Users should choose their target based on the effect size they want
to power the simulation to detect.

Note: `calibrate_k_inter()` automates this for GBSG-based DGMs.
For non-GBSG datasets, the manual grid search shown here is required.
This is a known limitation and will be addressed in a future release.
```

### 3.4 Null DGM doesn't document its `k_inter` [MEDIUM]

Lines 643–657 create the null DGM with `model = "null"`.  The `k_inter`
argument is absent, which means it defaults to `k_inter = 1` inside
`generate_aft_dgm_flex()`.  But this is implicit — a user might wonder
whether `k_inter` matters under the null.

**Suggested fix:** Add a comment:

```r
# Under model = "null", the treatment effect is uniform (no HTE),
# so k_inter is irrelevant and can be omitted (defaults to 1).
```


---

## 4. Simulation Machinery (Lines 690–811)

### 4.1 `confounders_actg` silently controls ForestSearch's search space [CRITICAL]

Line 701:
```r
confounders_actg <- dgm_factor_vars   # paste0("z", 1:12)
```

This vector is passed to `run_simulation_analysis()` as
`confounders_base`, which in turn becomes the set of candidate
variables for ForestSearch and GRF.  **This is the single most
important parameter for adoptability** — if a user's covariates are
named `age_bin`, `cd4_low`, etc., they must set `confounders_base`
accordingly.  But the document treats it as a one-liner.

**Suggested fix:** Add a dedicated subsection:

```markdown
### Defining the Candidate Variable Set

`confounders_base` tells `run_simulation_analysis()` which columns
in the simulated data to pass to ForestSearch and GRF as candidate
subgroup-defining variables.  These must:

1. Be present in the super-population (`dgm$df_super`)
2. Be binary factors (0/1) for ForestSearch
3. Match the names used in `factor_vars` when building the DGM

If your dataset has different variable names, update both
`factor_vars` (in the DGM call) and `confounders_base` (in the
simulation call) to match.
```

### 4.2 `fs_params` list silently overrides defaults [HIGH]

The `fs_params` list (lines 704–726) sets many ForestSearch parameters.
A user adapting this document has no guidance on which of these they
should change and which are safe defaults.

**Suggested fix:** Annotate each parameter with its role and when to
change it:

```r
fs_params <- list(
  # ── Column name mappings (must match simulate_from_dgm output) ──
  outcome.name           = "y_sim",
  event.name             = "event_sim",
  treat.name             = "treat_sim",
  id.name                = "id",

  # ── Variable screening (usually keep defaults) ──────────────────
  use_lasso              = TRUE,    # LASSO pre-screening
  use_grf                = TRUE,    # GRF pre-screening

  # ── Subgroup detection thresholds ───────────────────────────────
  # c1: HR must exceed this to flag a subgroup.
  # Increase for more conservative detection; decrease for more power.
  hr.threshold           = 1.25,
  # c2: Bootstrap-replicate HR consistency threshold.
  hr.consistency         = 1.0,
  # Proportion of bootstrap replicates exceeding c2.
  pconsistency.threshold = 0.90,

  # ── Computational settings ─────────────────────────────────────
  fs.splits              = 400L,   # Bootstrap replicates
  n.min                  = 60L,    # Minimum subgroup size (both arms)
  d0.min                 = 12L,    # Minimum events, control arm
  d1.min                 = 12L,    # Minimum events, treatment arm
  maxk                   = 2L,     # Maximum subgroup-defining variables
  by.risk                = 12L,    # Risk-table granularity (time units)

  # ...
)
```

### 4.3 `use_twostage = TRUE` is set without explanation [MEDIUM]

Line 721 enables two-stage estimation, but the document never mentions
what this does or why it matters.

### 4.4 Inner `future::plan("sequential")` is unexplained [MEDIUM]

Inside the `%dofuture%` body (line 752):
```r
future::plan("sequential")
```

This resets the plan *within each worker* to prevent nested parallelism.
Users unfamiliar with the `future` framework will not understand why this
is there, and may remove it, causing hard-to-debug failures.

**Suggested fix:** Add a comment:

```r
# Reset plan inside workers to avoid nested parallelism.
# Each worker runs ForestSearch sequentially; the outer loop
# parallelises across simulation replicates.
future::plan("sequential")
```


---

## 5. Results and Interpretation (Lines 814–911)

### 5.1 `format_oc_results()` metrics are not defined [HIGH]

Lines 822–832 call:
```r
format_oc_results(
  results_alt,
  metrics = c("detection", "classification", "cde_estimates"),
  ...
)
```

The metrics `"detection"`, `"classification"`, and `"cde_estimates"` are
magic strings.  The document never defines what each metric group
contains (e.g., detection rate, sensitivity, specificity, PPV for
"classification").  A user reading the output table has no glossary.

**Suggested fix:** Add a brief metrics glossary before the results
section:

```markdown
### Metrics Glossary

- **detection**: proportion of simulations where ForestSearch identifies
  *any* subgroup (`any.H = 1`).
- **classification**: among simulations with a detected subgroup,
  sensitivity (overlap with true H), specificity (non-overlap with Hc),
  and positive predictive value (PPV).
- **cde_estimates**: bias and RMSE of the controlled direct effect (CDE)
  estimator relative to the true DGM values.
```

### 5.2 `build_estimation_table()` survival-centric labels [MEDIUM]

Per your known deferred issues, `build_estimation_table()` uses labels
like "HR(H)" and "AHR(Hc)" which are survival-specific.  This is fine
for this document but should be noted for users who may later run
GLM-extension simulations.

### 5.3 `build_classification_table()` input structure is undocumented [MEDIUM]

Lines 883–910 build a `scenario_list` with specific field names
(`results`, `label`, `n_sample`, `dgm`, `hypothesis`).  This structure
is not documented in `?build_classification_table`, so a user must
reverse-engineer it from the vignette.


---

## 6. Detection Curve (Lines 914–964)

### 6.1 `compute_detection_probability()` parameters need context [MEDIUM]

The call uses `method = "cubature"` without explaining what this is or
what alternatives exist.  The `hr_threshold` and `hr_consistency`
parameters are passed from `fs_params` but the connection to the
detection probability formula is not made explicit.

### 6.2 Base R `abline()` mixed with ggplot [LOW]

Line 963 uses `abline()` on what is presumably a base-R plot from
`plot_detection_curve()`.  If a user expects ggplot (given the rest of
the document uses ggplot), this will be confusing.


---

## 7. KM Curve Construction (Lines 372–419)

### 7.1 GBSG KM plot uses `data.frame(survfit(...))` which may fail [LOW]

Line 377 wraps `survfit()` output in `data.frame()`.  This relies on an
implicit coercion method that may not produce the expected columns in all
versions of the `survival` package.  The ACTG175 side does manual
extraction (lines 393–402), which is more robust but inconsistent with
the GBSG side.


---

## 8. Cross-Cutting Transparency Issues

### 8.1 No "Your Dataset Here" checklist [CRITICAL]

The document's stated purpose is to serve as a template, but there is no
consolidated checklist a user can follow.  The information is scattered
across prose, code comments, and the comparison table.

**Suggested fix:** Add a §1.3 checklist:

```markdown
### Adaptation Checklist

To substitute your own trial data for ACTG175, you must address
each of the following:

- [ ] **Load your data** and create `id`, `y` (time), `event`, and
      `treat` (0/1) columns.
- [ ] **Choose confounders**: identify the candidate variables for
      subgroup search.
- [ ] **Discretise continuous variables**: create binary `z1`, `z2`,
      … columns.  Document your cutpoint rationale.
- [ ] **Define the harm subgroup**: set `dgm_subgroup_vars` and
      `dgm_subgroup_cuts` consistently with your factor coding.
- [ ] **Set time-scale parameters**: ensure `sim_analysis_time` is
      in the same units as your outcome variable.
- [ ] **Update `cox_formula_adj`**: include confounders appropriate
      for your study.
- [ ] **Update `confounders_base`**: must list all `z*` names passed
      to ForestSearch.
- [ ] **Run time-scale diagnostic**: verify `exp(dgm$model_params$mu)`
      is close to your observed median event time.
- [ ] **Calibrate `k_inter`**: run the grid search targeting a
      clinically meaningful HR(H).
- [ ] **Set simulation count**: use ≥ 500 replicates for stable OCs.
```

### 8.2 No censoring diagnostic [HIGH]

The document simulates data and reports event rates but never compares
the simulated censoring distribution against the observed data.
`check_censoring_dgm()` exists in the package for exactly this purpose
but is not called.

**Suggested fix:** Add a censoring diagnostic after the single-trial
simulation:

```r
check_censoring_dgm(
  dgm           = dgm_calibrated,
  n             = single_n,
  analysis_time = sim_analysis_time,
  seed          = single_seed
)
```

### 8.3 `nsims_alt = 44` and `nsims_null = 44` are too low without context [MEDIUM]

The Limitations section (line 1060) acknowledges this, but 44 is such
an unusual number that it invites confusion.  A user might wonder if
this is meaningful.

**Suggested fix:** Add an inline comment in the setup chunk:

```r
# Demonstration run; increase to >= 500 for stable operating
# characteristics (Monte Carlo SE for 50% detection ≈ 2.2% at 500 sims).
nsims_alt   <- 44L
nsims_null  <- 44L
```


---

## Summary Priority Matrix

| # | Issue | Severity | Section |
|---|-------|----------|---------|
| 8.1 | No adaptation checklist | Critical | Cross-cutting |
| 3.1 | No time-scale consistency check | Critical | DGM |
| 4.1 | `confounders_base` under-documented | Critical | Simulation |
| 1.1 | Cutpoint rationale missing | High | Setup |
| 1.2 | `dgm_subgroup_cuts` semantics unclear | High | Setup |
| 3.2 | `compute_dgm_cde()` unexplained | High | DGM |
| 3.3 | Calibration why/what missing | High | DGM |
| 4.2 | `fs_params` no guidance on what to change | High | Simulation |
| 5.1 | Metrics not defined | High | Results |
| 8.2 | No censoring diagnostic | High | Cross-cutting |
| 2.1 | Factor coding direction not summarised | High | Data prep |
| 1.3 | `cox_formula_adj` only 3 of 12 confounders | Medium | Setup |
| 1.4 | `analysis_time` needs time-scale anchoring | Medium | Setup |
| 2.2 | `flag_harm` construction opaque | Medium | Data prep |
| 3.4 | Null DGM `k_inter` implicit | Medium | DGM |
| 4.3 | `use_twostage` unexplained | Medium | Simulation |
| 4.4 | Inner `future::plan("sequential")` unexplained | Medium | Simulation |
| 5.2 | Survival-centric labels | Medium | Results |
| 5.3 | `build_classification_table()` input undocumented | Medium | Results |
| 6.1 | Detection probability method not explained | Medium | Detection |
| 8.3 | `nsims = 44` unusual without context | Medium | Cross-cutting |
| 1.5 | `fs_vi_grf_min = -0.2` unexplained | Low | Setup |
| 6.2 | Mixed base-R / ggplot | Low | Detection |
| 7.1 | KM plot construction inconsistent | Low | Comparison |
