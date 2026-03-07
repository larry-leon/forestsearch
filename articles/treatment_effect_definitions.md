# Treatment Effect Definitions in ForestSearch

## Introduction

ForestSearch implements multiple approaches for quantifying treatment
effects in survival analysis, each grounded in the potential outcomes
framework (Rubin, 1974). When identifying subgroups that experience
differential treatment effects, the choice of estimand matters:
different definitions capture different aspects of treatment effect
heterogeneity and carry distinct assumptions about the underlying
data-generating process. Understanding these distinctions is critical
for proper interpretation of ForestSearch results and for designing
simulation studies that faithfully evaluate operating characteristics.

This vignette defines the three treatment effect measures used
throughout the ForestSearch package, details their computation from an
Accelerated Failure Time (AFT) model with Weibull errors, clarifies the
distinction between DGM-level (super-population) and simulation-level
(per-replicate) hazard ratios, and discusses when each measure is most
informative. A summary comparison table is provided at the end for quick
reference.

## Notation and Setup

### Potential Outcomes Framework

For subject $i$ in a two-arm randomized trial with treatment
$A_{i} \in \{ 0,1\}$, the potential outcomes framework posits a pair of
survival times:

- $T_{i}(1)$: survival time had subject $i$ received treatment
  ($A = 1$),
- $T_{i}(0)$: survival time had subject $i$ received control ($A = 0$).

Only one potential outcome is observed for each subject. Let
$\mathbf{X}_{i}$ denote baseline covariates.

### AFT Model Specification

ForestSearch’s data-generating mechanisms are built on a Weibull AFT
model:

$$\log T_{i}(a) = \mu + \mathbf{X}_{i}(a)^{\top}{\mathbf{γ}} + \sigma\,\varepsilon_{i},\qquad\varepsilon_{i} \sim \text{Extreme Value}$$

where $\mu$ is the intercept, $\sigma > 0$ is the scale parameter, and
$\mathbf{γ}$ are regression coefficients on the log-time scale. For the
alternative-hypothesis model with heterogeneous treatment effects, the
covariate vector $\mathbf{X}_{i}(a)$ includes both the treatment
indicator and a treatment-by-subgroup interaction:

$$\mathbf{X}_{i}(a)^{\top}{\mathbf{γ}} = a \cdot \gamma_{\text{treat}} + \mathbf{Z}_{i}^{\top}{\mathbf{γ}}_{Z} + a \cdot H_{i} \cdot \gamma_{\text{int}}$$

where $H_{i} \in \{ 0,1\}$ flags membership in the harm subgroup and
${\mathbf{γ}}_{Z}$ are prognostic covariate effects.

### Hazard-Scale Coefficients

The AFT-to-hazard-scale transformation yields:

$${\mathbf{β}}_{0} = - {\mathbf{γ}}/\sigma$$

Under this parameterisation the conditional hazard function is

$$h\left( t \mid \mathbf{X}_{i} \right) = h_{0}(t)\exp\left( \mathbf{X}_{i}^{\top}{\mathbf{β}}_{0} \right)$$

and the individual log-hazard contributions under each treatment arm
are:

$$\theta_{i}(a) = \mathbf{X}_{i}(a)^{\top}{\mathbf{β}}_{0}$$

These are stored in the ForestSearch super-population data frame as
`theta_1` (treatment) and `theta_0` (control).

## Treatment Effect Definitions

### Marginal (Causal) Hazard Ratio

#### Definition

The **marginal hazard ratio** for a subgroup $\mathcal{S}$ is the
population-level treatment effect obtained by fitting a Cox proportional
hazards model to stacked potential outcomes within $\mathcal{S}$:

$$\text{HR}_{\text{marg}}(\mathcal{S}) = \exp\left( {\widehat{\beta}}_{\text{Cox}} \mid \text{data} = \left\{ \left( T_{i}(1),1,1 \right),\left( T_{i}(0),1,0 \right):i \in \mathcal{S} \right\} \right)$$

This is the coefficient from a Cox model fit to a dataset where each
subject contributes two rows—one under treatment, one under control—both
with event indicator 1 (since potential outcomes are fully observed in
the super-population).

#### Computation in ForestSearch

In
[`calculate_hazard_ratios()`](https://larry-leon.github.io/forestsearch/reference/calculate_hazard_ratios.md)
and
[`setup_gbsg_dgm()`](https://larry-leon.github.io/forestsearch/reference/setup_gbsg_dgm.md),
the marginal HR is computed as follows:

1.  Generate a common set of error terms:
    $\varepsilon_{i} = \log\left( \text{Exp}(1) \right)$ for
    $i = 1,\ldots,n$.

2.  Compute potential survival times:

    - $\log T_{i}(1) = \mu + \sigma\varepsilon_{i} + \mathbf{X}_{i}(1)^{\top}{\mathbf{γ}}$
    - $\log T_{i}(0) = \mu + \sigma\varepsilon_{i} + \mathbf{X}_{i}(0)^{\top}{\mathbf{γ}}$

3.  Stack into a single dataset of $2n$ rows with columns `time`,
    `event = 1`, `treat`, `flag_harm`.

4.  Fit Cox model:
    `coxph(Surv(time, event) ~ treat_sim, data = df_stacked)`.

5.  Extract:
    $\text{HR}_{\text{marg}} = \exp\left( {\widehat{\beta}}_{\text{treat}} \right)$.

#### Key Properties

The marginal HR uses the same $\varepsilon_{i}$ for both potential
outcomes, encoding the causal constraint that each subject’s latent
frailty is identical regardless of treatment assignment. This yields a
proper causal estimand under the Rubin framework. However, the Cox model
fit to stacked data targets a population-averaged effect and may deviate
from the conditional (individual-level) hazard ratio when the
proportional hazards assumption does not hold exactly.

#### Package Variables

| Variable                   | Description                                 |
|:---------------------------|:--------------------------------------------|
| `hr_causal` / `hr_overall` | Overall marginal HR                         |
| `hr_H_true`                | Marginal HR in harm subgroup $\mathcal{H}$  |
| `hr_Hc_true`               | Marginal HR in complement $\mathcal{H}^{c}$ |

### Average Hazard Ratio (AHR)

#### Definition

The **average hazard ratio** (AHR) for a subgroup $\mathcal{S}$ is
defined directly from the individual-level potential-outcome log hazard
ratios:

$$\text{AHR}(\mathcal{S}) = \exp\left( \frac{1}{|\mathcal{S}|}\sum\limits_{i \in \mathcal{S}}\left\lbrack \theta_{i}(1) - \theta_{i}(0) \right\rbrack \right) = \exp\left( \frac{1}{|\mathcal{S}|}\sum\limits_{i \in \mathcal{S}}\log\text{HR}_{i}^{(\text{po})} \right)$$

where $\log\text{HR}_{i}^{(\text{po})} = \theta_{i}(1) - \theta_{i}(0)$
is the individual causal log hazard ratio.

#### Computation in ForestSearch

The AHR is computed from the `loghr_po` column in the super-population
data frame:

1.  For each subject, compute potential-outcome log-hazards:
    - $\theta_{i}(0) = \mathbf{X}_{i}(0)^{\top}{\mathbf{β}}_{0}$
    - $\theta_{i}(1) = \mathbf{X}_{i}(1)^{\top}{\mathbf{β}}_{0}$
    - $\text{loghr\_po}_{i} = \theta_{i}(1) - \theta_{i}(0)$
2.  Aggregate over the subgroup:
    $\text{AHR}(\mathcal{S}) = \exp\left( {\overline{\text{loghr\_po}}}_{\mathcal{S}} \right)$

#### Key Properties

The AHR is a deterministic function of the model coefficients and
covariate distribution—it does not depend on the random error
$\varepsilon_{i}$. This makes it fully reproducible across simulation
replications (given the same super-population). Because it operates on
the log-hazard scale and averages individual effects, it captures
heterogeneity in treatment effect across the covariate distribution. For
a subgroup defined by a biomarker threshold $z \geq c$, the AHR curve
$\text{AHR}(c)$ provides a smooth summary of how the treatment effect
changes with the threshold.

#### Package Variables

| Variable                 | Description                                         |
|:-------------------------|:----------------------------------------------------|
| `loghr_po`               | Individual log HR: $\theta_{i}(1) - \theta_{i}(0)$  |
| `AHR` / `ahr`            | Overall average hazard ratio (DGM-level)            |
| `ahr_H` / `AHR_H_true`   | AHR in harm subgroup (DGM truth)                    |
| `ahr_Hc` / `AHR_Hc_true` | AHR in complement (DGM truth)                       |
| `ahr.H.hat`              | Per-replicate AHR in identified $\widehat{H}$       |
| `ahr.Hc.hat`             | Per-replicate AHR in identified ${\widehat{H}}^{c}$ |
| `ahr.H.true`             | Per-replicate AHR in true $H$                       |
| `ahr.Hc.true`            | Per-replicate AHR in true $H^{c}$                   |

### Controlled Direct Effect (CDE)

#### Definition

The **controlled direct effect** (CDE) hazard ratio for a subgroup
$\mathcal{S}$ is defined as the ratio of average exponentiated
log-hazards under treatment versus control:

$$\text{CDE}(\mathcal{S}) = \frac{\frac{1}{|\mathcal{S}|}\sum\limits_{i \in \mathcal{S}}\exp\left( \theta_{i}(1) \right)}{\frac{1}{|\mathcal{S}|}\sum\limits_{i \in \mathcal{S}}\exp\left( \theta_{i}(0) \right)} = \frac{{\overline{\exp\left( \theta_{i}(1) \right)}}_{\mathcal{S}}}{{\overline{\exp\left( \theta_{i}(0) \right)}}_{\mathcal{S}}}$$

This is the ratio of average hazard *contributions* (on the natural
scale) under treatment versus control. The term “controlled direct
effect” reflects that each $\theta_{i}(a)$ controls for baseline
covariates through the linear predictor.

#### Computation in ForestSearch

In
[`cox_ahr_cde_analysis()`](https://larry-leon.github.io/forestsearch/reference/cox_ahr_cde_analysis.md),
[`plot_subgroup_effects()`](https://larry-leon.github.io/forestsearch/reference/plot_subgroup_effects.md),
and
[`compute_dgm_cde()`](https://larry-leon.github.io/forestsearch/reference/compute_dgm_cde.md),
the CDE is computed from the `theta_1` and `theta_0` columns:

1.  For each subject, retrieve:
    - $\theta_{i}(0) = \mathbf{X}_{i}(0)^{\top}{\mathbf{β}}_{0}$
    - $\theta_{i}(1) = \mathbf{X}_{i}(1)^{\top}{\mathbf{β}}_{0}$
2.  Compute subgroup CDE:
    $\text{CDE}(\mathcal{S}) = {\overline{\exp\left( \theta(1) \right)}}_{\mathcal{S}}/{\overline{\exp\left( \theta(0) \right)}}_{\mathcal{S}}$

#### Key Properties

Like the AHR, the CDE is deterministic given the super-population (no
dependence on $\varepsilon_{i}$). However, by averaging on the hazard
(exponential) scale rather than the log-hazard scale, the CDE gives more
weight to subjects with larger absolute hazard contributions. Jensen’s
inequality implies that in general
$\text{CDE}(\mathcal{S}) \neq \text{AHR}(\mathcal{S})$, with the
discrepancy increasing as within-subgroup heterogeneity grows.

#### Package Variables

| Variable                                   | Description                                         |
|:-------------------------------------------|:----------------------------------------------------|
| `theta_0`                                  | Log-hazard contribution under control               |
| `theta_1`                                  | Log-hazard contribution under treatment             |
| `cde_H` / `cde_results[["recommend"]]`     | DGM-level CDE in harm subgroup (truth)              |
| `cde_Hc` / `cde_results[["questionable"]]` | DGM-level CDE in complement (truth)                 |
| `cde_results[["overall"]]`                 | Overall CDE                                         |
| `cde.H.hat`                                | Per-replicate CDE in identified $\widehat{H}$       |
| `cde.Hc.hat`                               | Per-replicate CDE in identified ${\widehat{H}}^{c}$ |
| `cde.H.true`                               | Per-replicate CDE in true $H$                       |
| `cde.Hc.true`                              | Per-replicate CDE in true $H^{c}$                   |

## Theoretical Relationships

### When All Three Agree

Under a constant treatment effect model (no heterogeneity), the
individual log hazard ratio is identical for all subjects:

$$\log\text{HR}_{i}^{(\text{po})} = \beta_{\text{treat}}\quad\forall\, i$$

In this case, all three measures coincide:
$\text{HR}_{\text{marg}} = \text{AHR} = \text{CDE} = \exp\left( \beta_{\text{treat}} \right)$.

### Sources of Divergence

When treatment effects are heterogeneous, three sources of divergence
arise.

First, **Jensen’s inequality** separates AHR from CDE:

$$\exp\left( \overline{\log\text{HR}_{i}} \right) \neq \overline{\exp\left( \theta_{i}(1) \right)}/\overline{\exp\left( \theta_{i}(0) \right)}$$

unless all $\log\text{HR}_{i}$ are identical.

Second, the marginal HR from the stacked Cox model can differ from the
AHR because the Cox partial likelihood implicitly weights subjects by
their risk-set contributions, creating a form of **informative
censoring** in the stacked framework.

Third, the marginal HR depends on the realized $\varepsilon_{i}$ draw,
introducing **Monte Carlo variability** absent from the AHR and CDE.

### Practical Implications for ForestSearch

In the ForestSearch simulation framework, the AHR is the primary target
estimand for characterising treatment effect heterogeneity because it is
deterministic and directly interpretable as the geometric mean hazard
ratio across subjects. The marginal HR serves as a validation target,
confirming that the AFT model’s causal structure produces sensible
population-level effects. The CDE provides a complementary perspective
that is particularly useful when comparing treatment effects across
subgroups defined by continuous biomarker thresholds.

## Computation Under the ForestSearch AFT Model

### Step-by-Step Derivation

Given a fitted AFT model with parameters $(\mu,\sigma,{\mathbf{γ}})$:

**Step 1: Hazard-scale transformation**

$${\mathbf{β}}_{0} = - {\mathbf{γ}}/\sigma$$

**Step 2: Build potential-outcome design matrices**

For each subject $i$, construct $\mathbf{X}_{i}(1)$ (treatment arm) and
$\mathbf{X}_{i}(0)$ (control arm). The treatment indicator is set to 1
or 0 respectively; all other covariates remain at their observed values.
Interaction terms (e.g., $A \times H$) are computed accordingly.

**Step 3: Individual log-hazard contributions**

$$\theta_{i}(a) = \mathbf{X}_{i}(a)^{\top}{\mathbf{β}}_{0},\quad a \in \{ 0,1\}$$

$$\log\text{HR}_{i}^{(\text{po})} = \theta_{i}(1) - \theta_{i}(0)$$

**Step 4: Treatment effect summaries** for subgroup $\mathcal{S}$

$$\text{AHR}(\mathcal{S}) = \exp\left( \frac{1}{|\mathcal{S}|}\sum\limits_{i \in \mathcal{S}}\log\text{HR}_{i}^{(\text{po})} \right)$$

$$\text{CDE}(\mathcal{S}) = \frac{\sum\limits_{i \in \mathcal{S}}\exp\left( \theta_{i}(1) \right)}{\sum\limits_{i \in \mathcal{S}}\exp\left( \theta_{i}(0) \right)}$$

$$\text{HR}_{\text{marg}}(\mathcal{S}):\text{Cox model on stacked potential outcomes}$$

### GBSG Example: Harm Subgroup

For the GBSG-based DGM with alternative hypothesis, the harm subgroup is
$\mathcal{H} = \{ z_{1} = 1\;\text{AND}\; z_{3} = 1\}$ (low estrogen
receptor AND premenopausal). The theoretical hazard ratios from the AFT
coefficients are:

$$\text{HR}(\mathcal{H}) = \exp\left( \beta_{0}\left\lbrack \text{treat} \right\rbrack + \beta_{0}\left\lbrack \text{zh} \right\rbrack \right)$$

$$\text{HR}\left( \mathcal{H}^{c} \right) = \exp\left( \beta_{0}\left\lbrack \text{treat} \right\rbrack \right)$$

The `k_inter` parameter modulates the interaction coefficient:
$\left. \gamma\left\lbrack \text{zh} \right\rbrack\leftarrow k_{\text{inter}} \times \gamma\left\lbrack \text{zh} \right\rbrack \right.$,
providing direct control over the magnitude of treatment effect
heterogeneity.

## DGM-Level vs. Simulation-Level Hazard Ratios

The variable names `hr.H.true` and `hr.Hc.true` appear in two distinct
contexts within the ForestSearch codebase. Both use the same
computational core—a Cox model fit within the true subgroup defined by
`flag_harm`—but they operate on different data and serve different
purposes. Understanding this distinction is essential for correctly
interpreting simulation operating characteristics.

### DGM-Level Truth (Super-Population)

When
[`setup_gbsg_dgm()`](https://larry-leon.github.io/forestsearch/reference/setup_gbsg_dgm.md)
or
[`calculate_hazard_ratios()`](https://larry-leon.github.io/forestsearch/reference/calculate_hazard_ratios.md)
constructs the data-generating mechanism, the marginal HR is computed
from the **stacked potential outcomes** in a large super-population
(typically $n_{\text{super}} = 5,000$). Every subject contributes
uncensored survival times under *both* treatment and control, using the
same random error $\varepsilon_{i}$ for both arms:

``` r
# In setup_gbsg_dgm() — stacked potential outcomes
df_po <- data.frame(
  time      = c(T_1, T_0),
  event     = 1L,
  treat     = c(rep(1, n_super), rep(0, n_super)),
  flag_harm = rep(flag_harm, 2)
)

hr_H_true <- exp(coxph(
  Surv(time, event) ~ treat_sim,
  data = subset(df_po, flag_harm == 1)
)$coefficients)

hr_Hc_true <- exp(coxph(
  Surv(time, event) ~ treat_sim,
  data = subset(df_po, flag_harm == 0)
)$coefficients)
```

These DGM-level values are essentially **fixed constants** for a given
DGM specification. They define the ground truth against which estimator
performance is evaluated and converge to the theoretical values:

$$\text{HR}_{\text{DGM}}(\mathcal{H}) = \exp\!(\beta_{0}\left\lbrack \text{treat} \right\rbrack + \beta_{0}\left\lbrack \text{zh} \right\rbrack),\qquad\text{HR}_{\text{DGM}}\left( \mathcal{H}^{c} \right) = \exp\!(\beta_{0}\left\lbrack \text{treat} \right\rbrack)$$

as $\left. n_{\text{super}}\rightarrow\infty \right.$.

**Stored as:** `dgm$hr_H_true`, `dgm$hr_Hc_true`, `dgm$hr_causal`;
equivalently in `dgm$hazard_ratios$harm_subgroup`,
`dgm$hazard_ratios$no_harm_subgroup`, `dgm$hazard_ratios$overall`.

### Simulation-Level Truth (Per-Replicate)

During each simulation replicate, `run_fs_analysis()` and
[`run_grf_analysis()`](https://larry-leon.github.io/forestsearch/reference/run_grf_analysis.md)
in `oc_analyses_gbsg_refactored.R` recompute the true subgroup HRs on
the **simulated sample**—with censoring, finite sample size, and
randomized (not stacked) treatment:

``` r
# In run_fs_analysis() — observed data, single treatment arm per subject
hr.H.true <- exp(coxph(
  Surv(y_sim, event_sim) ~ treat_sim,
  data = subset(df, flag_harm == 1)
)$coefficients)

hr.Hc.true <- exp(coxph(
  Surv(y_sim, event_sim) ~ treat_sim,
  data = subset(df, flag_harm == 0)
)$coefficients)
```

These per-replicate values are **random variables** that fluctuate
across simulation iterations due to sampling variability, censoring
patterns, and the particular treatment randomization.

**Stored as:** columns `hr.H.true` and `hr.Hc.true` in the simulation
results `data.table`, one row per replicate.

### Side-by-Side Comparison

| **DGM-Level vs. Simulation-Level True Hazard Ratios**                                                             |                                                  |                                                  |
|-------------------------------------------------------------------------------------------------------------------|--------------------------------------------------|--------------------------------------------------|
| Both use `coxph(Surv(time, event) ~ treat_sim)` within `flag_harm` subgroups, but on fundamentally different data |                                                  |                                                  |
|                                                                                                                   | **DGM-Level (Super-Population)**                 | **Simulation-Level (Per-Replicate)**             |
| Data source                                                                                                       | Stacked potential outcomes from super-population | Observed simulated trial data                    |
| Potential outcomes                                                                                                | Both T(1) and T(0) observed per subject          | Only T(A_i) observed per subject                 |
| Censoring                                                                                                         | None (event = 1 for all rows)                    | Yes (administrative + random censoring)          |
| Sample size                                                                                                       | 2 x n_super (e.g. 10,000 rows)                   | n_sample (e.g. 300 -- 700)                       |
| Treatment assignment                                                                                              | Both arms per subject (causal framework)         | Randomized (one arm per subject)                 |
| Subgroup membership                                                                                               | Known flag_harm (truth)                          | Known flag_harm (truth)                          |
| Nature of quantity                                                                                                | Essentially fixed (large N constant)             | Random (varies across replicates)                |
| Primary purpose                                                                                                   | Define ground truth for bias calculation         | Per-replicate reference for estimator comparison |
| Location in code                                                                                                  | setup_gbsg_dgm(), calculate_hazard_ratios()      | run_fs_analysis(), run_grf_analysis()            |
| Storage                                                                                                           | dgm\$hr_H_true, dgm\$hr_Hc_true                  | results\$hr.H.true, results\$hr.Hc.true          |

[ Code](#collapse-dgmvssimtable)

``` r
dgm_sim_df <- data.frame(
  Aspect = c(
    "Data source",
    "Potential outcomes",
    "Censoring",
    "Sample size",
    "Treatment assignment",
    "Subgroup membership",
    "Nature of quantity",
    "Primary purpose",
    "Location in code",
    "Storage"
  ),
  DGM_Level = c(
    "Stacked potential outcomes from super-population",
    "Both T(1) and T(0) observed per subject",
    "None (event = 1 for all rows)",
    "2 x n_super (e.g. 10,000 rows)",
    "Both arms per subject (causal framework)",
    "Known flag_harm (truth)",
    "Essentially fixed (large N constant)",
    "Define ground truth for bias calculation",
    "setup_gbsg_dgm(), calculate_hazard_ratios()",
    "dgm$hr_H_true, dgm$hr_Hc_true"
  ),
  Sim_Level = c(
    "Observed simulated trial data",
    "Only T(A_i) observed per subject",
    "Yes (administrative + random censoring)",
    "n_sample (e.g. 300 -- 700)",
    "Randomized (one arm per subject)",
    "Known flag_harm (truth)",
    "Random (varies across replicates)",
    "Per-replicate reference for estimator comparison",
    "run_fs_analysis(), run_grf_analysis()",
    "results$hr.H.true, results$hr.Hc.true"
  ),
  stringsAsFactors = FALSE
)

gt(dgm_sim_df) |>
  cols_label(
    Aspect = "",
    DGM_Level = md("**DGM-Level (Super-Population)**"),
    Sim_Level = md("**Simulation-Level (Per-Replicate)**")
  ) |>
  cols_width(
    Aspect ~ px(200),
    DGM_Level ~ px(340),
    Sim_Level ~ px(340)
  ) |>
  tab_header(
    title = md("**DGM-Level vs. Simulation-Level True Hazard Ratios**"),
    subtitle = md(
      paste0(
        "Both use `coxph(Surv(time, event) ~ treat_sim)` within `flag_harm` ",
        "subgroups, but on fundamentally different data"
      )
    )
  ) |>
  tab_style(
    style = cell_text(weight = "bold", size = "small"),
    locations = cells_body(columns = Aspect)
  ) |>
  tab_style(
    style = cell_text(size = "small"),
    locations = cells_body(columns = c(DGM_Level, Sim_Level))
  ) |>
  tab_style(
    style = cell_text(font = "monospace", size = "x-small"),
    locations = cells_body(
      columns = c(DGM_Level, Sim_Level),
      rows = c(9, 10)
    )
  ) |>
  tab_style(
    style = list(
      cell_fill(color = "#E8F4E8")
    ),
    locations = cells_body(
      columns = c(DGM_Level, Sim_Level),
      rows = c(7, 8)
    )
  ) |>
  tab_style(
    style = cell_borders(
      sides = "bottom",
      color = "#B0BEC5",
      weight = px(1)
    ),
    locations = cells_body(rows = c(6, 8))
  ) |>
  tab_options(
    table.font.size = px(13),
    heading.title.font.size = px(16),
    heading.subtitle.font.size = px(12),
    column_labels.font.weight = "bold",
    table.border.top.color = "#2C3E50",
    table.border.bottom.color = "#2C3E50"
  ) |>
  tab_style(
    style = list(
      cell_fill(color = "#2C3E50"),
      cell_text(color = "white", weight = "bold")
    ),
    locations = cells_column_labels()
  )
```

### Role in Operating Characteristics

In the simulation results produced by `run_oc_gbsg_sims()` and
summarised by
[`build_estimation_table()`](https://larry-leon.github.io/forestsearch/reference/build_estimation_table.md)
and
[`format_oc_results()`](https://larry-leon.github.io/forestsearch/reference/format_oc_results.md),
multiple HR columns work together to quantify estimation performance:

| Column       | Meaning                                                                                                                       |
|:-------------|:------------------------------------------------------------------------------------------------------------------------------|
| `hr.H.true`  | Per-replicate Cox HR in the **true** H, using observed data                                                                   |
| `hr.H.hat`   | Cox HR in the **identified** $\widehat{H}$, using observed data                                                               |
| `hr.H.bc`    | Bootstrap bias-corrected HR in $\widehat{H}$                                                                                  |
| `ahr.H.hat`  | AHR = $\exp\left( \overline{\text{loghr\_po}} \right)$ in the **identified** $\widehat{H}$                                    |
| `cde.H.hat`  | CDE = $\overline{\exp\left( \theta_{1} \right)}/\overline{\exp\left( \theta_{0} \right)}$ in the **identified** $\widehat{H}$ |
| `hr.Hc.true` | Per-replicate Cox HR in the **true** $H^{c}$, using observed data                                                             |
| `hr.Hc.hat`  | Cox HR in the **identified** ${\widehat{H}}^{c}$, using observed data                                                         |
| `hr.Hc.bc`   | Bootstrap bias-corrected HR in ${\widehat{H}}^{c}$                                                                            |
| `ahr.Hc.hat` | AHR in the **identified** ${\widehat{H}}^{c}$                                                                                 |
| `cde.Hc.hat` | CDE in the **identified** ${\widehat{H}}^{c}$                                                                                 |

Relative bias is computed against the DGM-level constant
$\theta_{H}^{\dagger} = {\mathtt{d}\mathtt{g}\mathtt{m}\mathtt{\$}\mathtt{h}\mathtt{r}\mathtt{\_}\mathtt{H}\mathtt{\_}\mathtt{t}\mathtt{r}\mathtt{u}\mathtt{e}}$
(the marginal/causal HR truth, using the $\dagger$ notation from León et
al., 2024):

$$b^{\dagger}(\%) = 100 \times \frac{\overline{\mathtt{h}\mathtt{r}\mathtt{.}\mathtt{H}\mathtt{.}\mathtt{h}\mathtt{a}\mathtt{t}} - \theta_{H}^{\dagger}}{\theta_{H}^{\dagger}}$$

When CDE truth values ($\theta^{\ddagger}$) are available,
[`build_estimation_table()`](https://larry-leon.github.io/forestsearch/reference/build_estimation_table.md)
reports a second bias column measuring distance from the controlled
direct effect:

$$b^{\ddagger}(\%) = 100 \times \frac{\overline{\mathtt{h}\mathtt{r}\mathtt{.}\mathtt{H}\mathtt{.}\mathtt{h}\mathtt{a}\mathtt{t}} - \theta_{H}^{\ddagger}}{\theta_{H}^{\ddagger}}$$

The notation aligns across both summary functions: $H$ denotes the true
(oracle) subgroup, $\widehat{H}$ denotes the identified (data-driven)
subgroup, and $*$ is reserved exclusively for bootstrap bias-corrected
estimates (${\widehat{\theta}}^{*}$)—it is never used to denote “true
subgroup.”

This decomposition allows the simulation to separate two sources of
error: (1) **subgroup misidentification** (the identified $\widehat{H}$
may not match the true $H$), and (2) **estimation variability** (even if
$\widehat{H} = H$, the Cox coefficient estimate fluctuates in finite
samples).

The per-replicate `hr.H.true` helps diagnose the second source. When
`mean(hr.H.true)` across replicates closely matches `dgm$hr_H_true`, the
simulation is confirming that the DGM produces the intended treatment
effect in finite samples (on average). Any gap between `mean(hr.H.hat)`
and `mean(hr.H.true)` is attributable to subgroup misidentification,
while the gap between `mean(hr.H.true)` and `dgm$hr_H_true` reflects
finite-sample estimation noise.

### Companion AHR and CDE Metrics

For each context, parallel **Average Hazard Ratio** and **Controlled
Direct Effect** values are available:

**DGM-level AHR:** Deterministic, computed from individual `loghr_po`
values in the super-population.

``` r
AHR_H_true  <- exp(mean(df_super$loghr_po[df_super$flag_harm == 1]))
AHR_Hc_true <- exp(mean(df_super$loghr_po[df_super$flag_harm == 0]))
```

**DGM-level CDE:** Deterministic, computed from individual `theta_0` and
`theta_1` values in the super-population.

``` r
CDE_H_true  <- mean(exp(df_super$theta_1[df_super$flag_harm == 1])) /
               mean(exp(df_super$theta_0[df_super$flag_harm == 1]))
CDE_Hc_true <- mean(exp(df_super$theta_1[df_super$flag_harm == 0])) /
               mean(exp(df_super$theta_0[df_super$flag_harm == 0]))
```

These are computed by
[`compute_dgm_cde()`](https://larry-leon.github.io/forestsearch/reference/compute_dgm_cde.md)
and stored in the DGM’s `hazard_ratios` list.

**Simulation-level AHR:** Can be computed per replicate when `loghr_po`
is carried through to the simulated sample (available via
[`simulate_from_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md)).
Stored as `ahr.H.hat` and `ahr.Hc.hat` in the results `data.table`.

**Simulation-level CDE:** Similarly computed per replicate from
`theta_0` and `theta_1` in the simulated sample. Stored as `cde.H.hat`
and `cde.Hc.hat` in the results `data.table`.

The DGM-level AHR and Cox-based marginal HR generally agree closely but
can diverge when within-subgroup treatment effects are highly
heterogeneous, as discussed in Section @ref(theoretical-relationships).
The CDE generally tracks the AHR but diverges due to Jensen’s
inequality, with the gap increasing as within-subgroup heterogeneity
grows.

In
[`build_estimation_table()`](https://larry-leon.github.io/forestsearch/reference/build_estimation_table.md),
all four estimator types—Cox HR ($\widehat{\theta}$), bootstrap
bias-corrected Cox HR (${\widehat{\theta}}^{*}$), AHR
($\widehat{a}\text{hr}$), and CDE ($\theta^{\ddagger}$)—can appear as
rows, each with appropriate bias columns. AHR and CDE rows appear
conditionally when their source columns are present in the simulation
results.

## Summary Comparison Table

The following `gt` table provides a side-by-side comparison of all three
treatment effect measures.

| **Comparison of Treatment Effect Definitions**                                                                                                     |                                                                            |                                                                               |                                                                          |
|----------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------|-------------------------------------------------------------------------------|--------------------------------------------------------------------------|
| Marginal HR, Average Hazard Ratio, and Controlled Direct Effect                                                                                    |                                                                            |                                                                               |                                                                          |
|                                                                                                                                                    | **Marginal HR**                                                            | **AHR**                                                                       | **CDE**                                                                  |
| **Definition**                                                                                                                                     |                                                                            |                                                                               |                                                                          |
| Formal name                                                                                                                                        | Marginal (causal) hazard ratio                                             | Average hazard ratio                                                          | Controlled direct effect (hazard ratio)                                  |
| Notation                                                                                                                                           | HR_marg(S), theta-dagger                                                   | AHR(S)                                                                        | CDE(S), theta-ddagger                                                    |
| Formula (subgroup S)                                                                                                                               | exp(beta_Cox) from stacked PO Cox model                                    | exp( mean(loghr_po_i) ) for i in S                                            | mean(exp(theta_1_i)) / mean(exp(theta_0_i))                              |
| **Statistical Properties**                                                                                                                         |                                                                            |                                                                               |                                                                          |
| Averaging scale                                                                                                                                    | Hazard ratio (Cox partial likelihood)                                      | Log-hazard (then exponentiated)                                               | Hazard (natural scale ratio)                                             |
| Depends on epsilon_i                                                                                                                               | Yes                                                                        | No                                                                            | No                                                                       |
| Reproducible across seeds                                                                                                                          | No (Monte Carlo variability)                                               | Yes (deterministic)                                                           | Yes (deterministic)                                                      |
| Accounts for covariate heterogeneity                                                                                                               | Averaged out (population-level)                                            | Geometric mean preserves individual effects                                   | Ratio of mean hazards preserves scale                                    |
| Sensitive to outlier hazard contributions                                                                                                          | No (rank-based)                                                            | Moderate (log scale dampens extremes)                                         | Yes (exponential scale amplifies large values)                           |
| **Computation**                                                                                                                                    |                                                                            |                                                                               |                                                                          |
| Standard estimation method                                                                                                                         | Cox PH on stacked potential outcomes                                       | Arithmetic mean of individual log HRs                                         | Ratio of average exponentiated log-hazards                               |
| Primary ForestSearch use                                                                                                                           | Estimation rows in build_estimation_table(); DGM validation                | Estimation row in build_estimation_table(); primary subgroup characterization | Estimation row in build_estimation_table(); biomarker-threshold analysis |
| Key package variable(s)                                                                                                                            | hr_causal, hr_H_true, hr_Hc_true, hr.H.hat, hr.Hc.hat                      | loghr_po, ahr.H.hat, ahr.Hc.hat, ahr_H, ahr_Hc                                | theta_0, theta_1, cde.H.hat, cde.Hc.hat, cde_H, cde_Hc                   |
| **Interpretation & Use**                                                                                                                           |                                                                            |                                                                               |                                                                          |
| Key package function(s)                                                                                                                            | calculate_hazard_ratios(), setup_gbsg_dgm(), build_estimation_table()      | cox_ahr_cde_analysis(), build_estimation_table()                              | cox_ahr_cde_analysis(), compute_dgm_cde(), build_estimation_table()      |
| Interpretation                                                                                                                                     | Population-average HR as estimated by Cox regression on potential outcomes | Geometric mean of individual causal hazard ratios across subjects             | Ratio of average hazard contributions under treatment vs. control        |
| Based on the potential outcomes framework within the ForestSearch AFT data-generating mechanism. See León et al. (2024), *Statistics in Medicine*. |                                                                            |                                                                               |                                                                          |

[ Code](#collapse-comparisontable)

``` r
# Build comparison data frame
comparison_df <- data.frame(
  Property = c(
    "Formal name",
    "Notation",
    "Formula (subgroup S)",
    "Averaging scale",
    "Depends on epsilon_i",
    "Reproducible across seeds",
    "Accounts for covariate heterogeneity",
    "Sensitive to outlier hazard contributions",
    "Standard estimation method",
    "Primary ForestSearch use",
    "Key package variable(s)",
    "Key package function(s)",
    "Interpretation"
  ),
  Marginal_HR = c(
    "Marginal (causal) hazard ratio",
    "HR_marg(S), theta-dagger",
    "exp(beta_Cox) from stacked PO Cox model",
    "Hazard ratio (Cox partial likelihood)",
    "Yes",
    "No (Monte Carlo variability)",
    "Averaged out (population-level)",
    "No (rank-based)",
    "Cox PH on stacked potential outcomes",
    "Estimation rows in build_estimation_table(); DGM validation",
    "hr_causal, hr_H_true, hr_Hc_true, hr.H.hat, hr.Hc.hat",
    "calculate_hazard_ratios(), setup_gbsg_dgm(), build_estimation_table()",
    "Population-average HR as estimated by Cox regression on potential outcomes"
  ),
  AHR = c(
    "Average hazard ratio",
    "AHR(S)",
    "exp( mean(loghr_po_i) ) for i in S",
    "Log-hazard (then exponentiated)",
    "No",
    "Yes (deterministic)",
    "Geometric mean preserves individual effects",
    "Moderate (log scale dampens extremes)",
    "Arithmetic mean of individual log HRs",
    "Estimation row in build_estimation_table(); primary subgroup characterization",
    "loghr_po, ahr.H.hat, ahr.Hc.hat, ahr_H, ahr_Hc",
    "cox_ahr_cde_analysis(), build_estimation_table()",
    "Geometric mean of individual causal hazard ratios across subjects"
  ),
  CDE = c(
    "Controlled direct effect (hazard ratio)",
    "CDE(S), theta-ddagger",
    "mean(exp(theta_1_i)) / mean(exp(theta_0_i))",
    "Hazard (natural scale ratio)",
    "No",
    "Yes (deterministic)",
    "Ratio of mean hazards preserves scale",
    "Yes (exponential scale amplifies large values)",
    "Ratio of average exponentiated log-hazards",
    "Estimation row in build_estimation_table(); biomarker-threshold analysis",
    "theta_0, theta_1, cde.H.hat, cde.Hc.hat, cde_H, cde_Hc",
    "cox_ahr_cde_analysis(), compute_dgm_cde(), build_estimation_table()",
    "Ratio of average hazard contributions under treatment vs. control"
  ),
  stringsAsFactors = FALSE
)

# Create publication-quality gt table
gt(comparison_df) |>
  cols_label(
    Property = "",
    Marginal_HR = md("**Marginal HR**"),
    AHR = md("**AHR**"),
    CDE = md("**CDE**")
  ) |>
  cols_width(
    Property ~ px(260),
    Marginal_HR ~ px(280),
    AHR ~ px(280),
    CDE ~ px(280)
  ) |>
  tab_header(
    title = md("**Comparison of Treatment Effect Definitions**"),
    subtitle = md(
      "Marginal HR, Average Hazard Ratio, and Controlled Direct Effect"
    )
  ) |>
  tab_source_note(
    source_note = md(
      paste0(
        "Based on the potential outcomes framework within the ",
        "ForestSearch AFT data-generating mechanism. ",
        "See Le\u00f3n et al. (2024), *Statistics in Medicine*."
      )
    )
  ) |>
  tab_row_group(
    label = md("**Definition**"),
    rows = 1:3
  ) |>
  tab_row_group(
    label = md("**Statistical Properties**"),
    rows = 4:8
  ) |>
  tab_row_group(
    label = md("**Computation**"),
    rows = 9:11
  ) |>
  tab_row_group(
    label = md("**Interpretation & Use**"),
    rows = 12:13
  ) |>
  row_group_order(
    groups = c("**Definition**", "**Statistical Properties**",
               "**Computation**", "**Interpretation & Use**")
  ) |>
  tab_style(
    style = list(
      cell_fill(color = "#F0F4F8"),
      cell_text(weight = "bold", size = "small")
    ),
    locations = cells_row_groups()
  ) |>
  tab_style(
    style = cell_text(weight = "bold", size = "small"),
    locations = cells_body(columns = Property)
  ) |>
  tab_style(
    style = cell_text(size = "small"),
    locations = cells_body(columns = c(Marginal_HR, AHR, CDE))
  ) |>
  tab_style(
    style = list(
      cell_fill(color = "#E8F4E8")
    ),
    locations = cells_body(
      columns = c(Marginal_HR, AHR, CDE),
      rows = c(5, 6)
    )
  ) |>
  tab_style(
    style = cell_borders(
      sides = "bottom",
      color = "#B0BEC5",
      weight = px(1)
    ),
    locations = cells_body(rows = c(3, 8, 11))
  ) |>
  tab_options(
    table.font.size = px(13),
    heading.title.font.size = px(16),
    heading.subtitle.font.size = px(13),
    column_labels.font.weight = "bold",
    table.border.top.color = "#2C3E50",
    table.border.bottom.color = "#2C3E50",
    row_group.border.top.color = "#B0BEC5",
    row_group.border.bottom.color = "#B0BEC5",
    source_notes.font.size = px(11)
  ) |>
  tab_style(
    style = list(
      cell_fill(color = "#2C3E50"),
      cell_text(color = "white", weight = "bold")
    ),
    locations = cells_column_labels()
  )
```

## Detailed Formulas: Quick-Reference Card

| **Treatment Effect Computation Pipeline**      |                                                                             |                                                                             |
|------------------------------------------------|-----------------------------------------------------------------------------|-----------------------------------------------------------------------------|
| From AFT model fit to subgroup-level estimands |                                                                             |                                                                             |
| **Step**                                       | **Mathematical Formula**                                                    | **R Implementation**                                                        |
| 1\. AFT model                                  | log(T_i) = mu + X_i' gamma + sigma \* epsilon_i                             | survreg(Surv(y, event) ~ ..., dist = 'weibull')                             |
| 2\. Hazard-scale transform                     | beta_0 = -gamma / sigma                                                     | b0 \<- -gamma / tau                                                         |
| 3\. Individual log-hazards                     | theta_i(a) = X_i(a)' beta_0                                                 | theta_0 \<- X_control %\*% b0; theta_1 \<- X_treat %\*% b0                  |
| 4\. Individual causal log HR                   | loghr_po_i = theta_i(1) - theta_i(0)                                        | loghr_po \<- theta_1 - theta_0                                              |
| 5a. AHR (subgroup S)                           | AHR(S) = exp( mean(loghr_po_i, i in S) )                                    | exp(mean(df\$loghr_po\[idx\]))                                              |
| 5b. CDE (subgroup S)                           | CDE(S) = mean(exp(theta_1_i)) / mean(exp(theta_0_i))                        | mean(exp(df\$theta_1\[idx\])) / mean(exp(df\$theta_0\[idx\]))               |
| 5c. Marginal HR (subgroup S)                   | HR_marg(S): exp(beta) from coxph() on stacked PO data                       | exp(coxph(Surv(time,event)~treat_sim, data=df_stacked)\$coef)               |
| 6\. Theoretical HR (harm subgroup)             | HR(H) = exp( beta_0\[treat\] + beta_0\[zh\] )                               | exp(b0\['treat'\] + b0\['zh'\])                                             |
| 7\. Per-replicate truth (sim-level)            | hr.H.true: coxph() on observed data, flag_harm == 1                         | exp(coxph(Surv(y_sim,event_sim)~treat_sim, subset(df, flag_harm==1))\$coef) |
| 8a. Per-replicate estimate (identified, Cox)   | hr.H.hat: coxph() on observed data, sg_hat == 1                             | exp(coxph(Surv(y_sim,event_sim)~treat_sim, subset(df, sg_hat==1))\$coef)    |
| 8b. Per-replicate AHR (identified)             | ahr.H.hat: exp(mean(loghr_po_i)) for i in sg_hat == 1                       | exp(mean(df\$loghr_po\[df\$sg_hat == 1\]))                                  |
| 8c. Per-replicate CDE (identified)             | cde.H.hat: mean(exp(theta_1_i)) / mean(exp(theta_0_i)) for i in sg_hat == 1 | mean(exp(df\$theta_1\[sg_hat\])) / mean(exp(df\$theta_0\[sg_hat\]))         |

[ Code](#collapse-pipelinetable)

``` r
formula_df <- data.frame(
  Step = c(
    "1. AFT model",
    "2. Hazard-scale transform",
    "3. Individual log-hazards",
    "4. Individual causal log HR",
    "5a. AHR (subgroup S)",
    "5b. CDE (subgroup S)",
    "5c. Marginal HR (subgroup S)",
    "6. Theoretical HR (harm subgroup)",
    "7. Per-replicate truth (sim-level)",
    "8a. Per-replicate estimate (identified, Cox)",
    "8b. Per-replicate AHR (identified)",
    "8c. Per-replicate CDE (identified)"
  ),
  Formula = c(
    "log(T_i) = mu + X_i' gamma + sigma * epsilon_i",
    "beta_0 = -gamma / sigma",
    "theta_i(a) = X_i(a)' beta_0",
    "loghr_po_i = theta_i(1) - theta_i(0)",
    "AHR(S) = exp( mean(loghr_po_i, i in S) )",
    "CDE(S) = mean(exp(theta_1_i)) / mean(exp(theta_0_i))",
    "HR_marg(S): exp(beta) from coxph() on stacked PO data",
    "HR(H) = exp( beta_0[treat] + beta_0[zh] )",
    "hr.H.true: coxph() on observed data, flag_harm == 1",
    "hr.H.hat: coxph() on observed data, sg_hat == 1",
    "ahr.H.hat: exp(mean(loghr_po_i)) for i in sg_hat == 1",
    "cde.H.hat: mean(exp(theta_1_i)) / mean(exp(theta_0_i)) for i in sg_hat == 1"
  ),
  ForestSearch_Code = c(
    "survreg(Surv(y, event) ~ ..., dist = 'weibull')",
    "b0 <- -gamma / tau",
    "theta_0 <- X_control %*% b0; theta_1 <- X_treat %*% b0",
    "loghr_po <- theta_1 - theta_0",
    "exp(mean(df$loghr_po[idx]))",
    "mean(exp(df$theta_1[idx])) / mean(exp(df$theta_0[idx]))",
    "exp(coxph(Surv(time,event)~treat_sim, data=df_stacked)$coef)",
    "exp(b0['treat'] + b0['zh'])",
    "exp(coxph(Surv(y_sim,event_sim)~treat_sim, subset(df, flag_harm==1))$coef)",
    "exp(coxph(Surv(y_sim,event_sim)~treat_sim, subset(df, sg_hat==1))$coef)",
    "exp(mean(df$loghr_po[df$sg_hat == 1]))",
    "mean(exp(df$theta_1[sg_hat])) / mean(exp(df$theta_0[sg_hat]))"
  ),
  stringsAsFactors = FALSE
)

gt(formula_df) |>
  cols_label(
    Step = md("**Step**"),
    Formula = md("**Mathematical Formula**"),
    ForestSearch_Code = md("**R Implementation**")
  ) |>
  cols_width(
    Step ~ px(220),
    Formula ~ px(380),
    ForestSearch_Code ~ px(380)
  ) |>
  tab_header(
    title = md("**Treatment Effect Computation Pipeline**"),
    subtitle = "From AFT model fit to subgroup-level estimands"
  ) |>
  tab_style(
    style = cell_text(font = "monospace", size = "small"),
    locations = cells_body(columns = ForestSearch_Code)
  ) |>
  tab_style(
    style = cell_text(size = "small"),
    locations = cells_body(columns = c(Step, Formula))
  ) |>
  tab_style(
    style = cell_text(weight = "bold", size = "small"),
    locations = cells_body(columns = Step)
  ) |>
  tab_style(
    style = list(
      cell_fill(color = "#FFF8E1")
    ),
    locations = cells_body(rows = 5:7)
  ) |>
  tab_style(
    style = list(
      cell_fill(color = "#E8F4E8")
    ),
    locations = cells_body(rows = 9:12)
  ) |>
  tab_options(
    table.font.size = px(13),
    column_labels.font.weight = "bold"
  ) |>
  tab_style(
    style = list(
      cell_fill(color = "#2C3E50"),
      cell_text(color = "white", weight = "bold")
    ),
    locations = cells_column_labels()
  )
```

## Worked Example

The following code demonstrates how the three measures are computed from
a single DGM object, highlighting where they agree and diverge.

``` r
library(forestsearch)

# Create a GBSG-based DGM with heterogeneous treatment effects
dgm <- setup_gbsg_dgm(
  model   = "alt",
  k_treat = 1.0,
  k_inter = 2.0,
  verbose = TRUE
)

# --- AHR (deterministic, from loghr_po) ---
df <- dgm$df_super_rand
ahr_overall <- exp(mean(df$loghr_po))
ahr_harm    <- exp(mean(df$loghr_po[df$flag_harm == 1]))
ahr_noharm  <- exp(mean(df$loghr_po[df$flag_harm == 0]))

cat("AHR (overall):", round(ahr_overall, 4), "\n")
cat("AHR (harm):   ", round(ahr_harm, 4), "\n")
cat("AHR (Hc):     ", round(ahr_noharm, 4), "\n")

# --- CDE (deterministic, from theta_0 and theta_1) ---
cde_overall <- mean(exp(df$theta_1)) / mean(exp(df$theta_0))
cde_harm    <- mean(exp(df$theta_1[df$flag_harm == 1])) /
               mean(exp(df$theta_0[df$flag_harm == 1]))
cde_noharm  <- mean(exp(df$theta_1[df$flag_harm == 0])) /
               mean(exp(df$theta_0[df$flag_harm == 0]))

cat("CDE (overall):", round(cde_overall, 4), "\n")
cat("CDE (harm):   ", round(cde_harm, 4), "\n")
cat("CDE (Hc):     ", round(cde_noharm, 4), "\n")

# Or equivalently via compute_dgm_cde():
dgm <- compute_dgm_cde(dgm)
cat("\nCDE from DGM (harm):", round(dgm$hazard_ratios$CDE_harm, 4), "\n")
cat("CDE from DGM (Hc):  ", round(dgm$hazard_ratios$CDE_no_harm, 4), "\n")

# --- Marginal HR: DGM-level (stochastic, from stacked Cox model) ---
cat("\n--- DGM-Level Marginal HR (super-population truth) ---\n")
cat("HR_marg (overall):", round(dgm$hr_causal, 4), "\n")
cat("HR_marg (harm):   ", round(dgm$hr_H_true, 4), "\n")
cat("HR_marg (Hc):     ", round(dgm$hr_Hc_true, 4), "\n")

# --- Simulation-level truth (per-replicate, observed data) ---
cat("\n--- Simulation-Level HR (single replicate) ---\n")
sim_data <- simulate_from_dgm(dgm, n = 500, seed = 1)

hr_H_sim <- exp(survival::coxph(
  survival::Surv(y_sim, event_sim) ~ treat_sim,
  data = subset(sim_data, flag_harm == 1)
)$coefficients)

hr_Hc_sim <- exp(survival::coxph(
  survival::Surv(y_sim, event_sim) ~ treat_sim,
  data = subset(sim_data, flag_harm == 0)
)$coefficients)

cat("hr.H.true  (this replicate):", round(hr_H_sim, 4), "\n")
cat("hr.Hc.true (this replicate):", round(hr_Hc_sim, 4), "\n")
cat("\nDGM vs sim-level gap (H):  ",
    round(hr_H_sim - dgm$hr_H_true, 4), "(finite-sample noise)\n")

# --- Biomarker-threshold analysis ---
results <- cox_ahr_cde_analysis(
  df          = dgm$df_super_rand,
  z_name      = "z_er",
  hr_threshold = 0.7,
  plot_style  = "grid",
  verbose     = TRUE
)
```

## Guidance for Choosing an Estimand

The choice among these three measures depends on the analytic objective:

**Use AHR** as the default estimand for ForestSearch simulation studies.
It is deterministic, directly interpretable, and aligns with the
geometric-mean interpretation of treatment effect heterogeneity. The AHR
is a primary target when evaluating operating characteristics
(sensitivity, PPV, bias) across simulation replications.

**Use marginal HR** for validation purposes. The stacked-Cox approach
provides a familiar, model-based estimand that confirms the DGM’s causal
structure produces the intended population-level effects. It is also the
natural target for analyses of observed trial data where potential
outcomes are not directly available. In
[`build_estimation_table()`](https://larry-leon.github.io/forestsearch/reference/build_estimation_table.md),
the naive Cox HR ($\widehat{\theta}$) and bias-corrected
(${\widehat{\theta}}^{*}$) rows report bias against the marginal truth
($\theta^{\dagger}$) and, when available, against the CDE truth
($\theta^{\ddagger}$).

**Use CDE** as a complementary estimand that provides a natural-scale
perspective on treatment effect heterogeneity. The CDE appears as an
estimation row ($\theta^{\ddagger}\left( \widehat{H} \right)$) in
[`build_estimation_table()`](https://larry-leon.github.io/forestsearch/reference/build_estimation_table.md)
alongside Cox HR and AHR rows, enabling direct comparison of estimation
properties across all three estimand types. When exploring treatment
effects across continuous biomarker thresholds via
[`cox_ahr_cde_analysis()`](https://larry-leon.github.io/forestsearch/reference/cox_ahr_cde_analysis.md),
the CDE’s natural-scale averaging can reveal patterns that are
attenuated on the log scale, particularly when hazard contributions are
highly variable.

## References

- León, L.F., Jemielita, T., Guo, Z., Marceau West, R., & Anderson, K.M.
  (2024). Exploratory subgroup identification in the heterogeneous Cox
  model: A relatively simple procedure. *Statistics in Medicine*.
  <doi:10.1002/sim.10163>.

- Rubin, D.B. (1974). Estimating causal effects of treatments in
  randomized and nonrandomized studies. *Journal of Educational
  Psychology*, 66(5), 688–701.

- Hernán, M.A. & Robins, J.M. (2020). *Causal Inference: What If*. Boca
  Raton: Chapman & Hall/CRC.

- Kalbfleisch, J.D. & Prentice, R.L. (2002). *The Statistical Analysis
  of Failure Time Data* (2nd ed.). Wiley.
