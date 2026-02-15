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
Accelerated Failure Time (AFT) model with Weibull errors, and discusses
when each measure is most informative. A summary comparison table is
provided at the end for quick reference.

## Notation and Setup

### Potential Outcomes Framework

For subject $`i`$ in a two-arm randomized trial with treatment
$`A_i \in \{0, 1\}`$, the potential outcomes framework posits a pair of
survival times:

- $`T_i(1)`$: survival time had subject $`i`$ received treatment
  ($`A = 1`$),
- $`T_i(0)`$: survival time had subject $`i`$ received control
  ($`A = 0`$).

Only one potential outcome is observed for each subject. Let
$`\mathbf{X}_i`$ denote baseline covariates.

### AFT Model Specification

ForestSearch’s data-generating mechanisms are built on a Weibull AFT
model:

``` math

\log T_i(a) = \mu + \mathbf{X}_i(a)^{\top} \boldsymbol{\gamma}
+ \sigma \, \varepsilon_i,
\qquad \varepsilon_i \sim \text{Extreme Value}
```

where $`\mu`$ is the intercept, $`\sigma > 0`$ is the scale parameter,
and $`\boldsymbol{\gamma}`$ are regression coefficients on the log-time
scale. For the alternative-hypothesis model with heterogeneous treatment
effects, the covariate vector $`\mathbf{X}_i(a)`$ includes both the
treatment indicator and a treatment-by-subgroup interaction:

``` math

\mathbf{X}_i(a)^{\top} \boldsymbol{\gamma}
= a \cdot \gamma_{\text{treat}} + \mathbf{Z}_i^{\top} \boldsymbol{\gamma}_Z
+ a \cdot H_i \cdot \gamma_{\text{int}}
```

where $`H_i \in \{0, 1\}`$ flags membership in the harm subgroup and
$`\boldsymbol{\gamma}_Z`$ are prognostic covariate effects.

### Hazard-Scale Coefficients

The AFT-to-hazard-scale transformation yields:

``` math

\boldsymbol{\beta}_0 = -\boldsymbol{\gamma} / \sigma
```

Under this parameterisation the conditional hazard function is

``` math

h(t \mid \mathbf{X}_i) = h_0(t) \exp\left(\mathbf{X}_i^{\top}
\boldsymbol{\beta}_0\right)
```

and the individual log-hazard contributions under each treatment arm
are:

``` math

\theta_i(a) = \mathbf{X}_i(a)^{\top} \boldsymbol{\beta}_0
```

These are stored in the ForestSearch super-population data frame as
`theta_1` (treatment) and `theta_0` (control).

## Treatment Effect Definitions

### Marginal (Causal) Hazard Ratio

#### Definition

The **marginal hazard ratio** for a subgroup $`\mathcal{S}`$ is the
population-level treatment effect obtained by fitting a Cox proportional
hazards model to stacked potential outcomes within $`\mathcal{S}`$:

``` math

\text{HR}_{\text{marg}}(\mathcal{S})
= \exp\left(\hat{\beta}_{\text{Cox}}  \mid 
\text{data} = \left\{\left(T_i(1), 1, 1\right), \left(T_i(0), 1, 0\right)
: i \in \mathcal{S}\right\}\right)
```

This is the coefficient from a Cox model fit to a dataset where each
subject contributes two rows—one under treatment, one under control—both
with event indicator 1 (since potential outcomes are fully observed in
the super-population).

#### Computation in ForestSearch

In
[`calculate_hazard_ratios()`](https://larry-leon.github.io/forestsearch/reference/calculate_hazard_ratios.md)
and
[`create_gbsg_dgm()`](https://larry-leon.github.io/forestsearch/reference/create_gbsg_dgm.md),
the marginal HR is computed as follows:

1.  Generate a common set of error terms:
    $`\varepsilon_i = \log(\text{Exp}(1))`$ for $`i = 1, \ldots, n`$.

2.  Compute potential survival times:

    - $`\log T_i(1) = \mu + \sigma \varepsilon_i + \mathbf{X}_i(1)^{\top} \boldsymbol{\gamma}`$
    - $`\log T_i(0) = \mu + \sigma \varepsilon_i + \mathbf{X}_i(0)^{\top} \boldsymbol{\gamma}`$

3.  Stack into a single dataset of $`2n`$ rows with columns `time`,
    `event = 1`, `treat`, `flag_harm`.

4.  Fit Cox model:
    `coxph(Surv(time, event) ~ treat, data = df_stacked)`.

5.  Extract:
    $`\text{HR}_{\text{marg}} = \exp(\hat{\beta}_{\text{treat}})`$.

#### Key Properties

The marginal HR uses the same $`\varepsilon_i`$ for both potential
outcomes, encoding the causal constraint that each subject’s latent
frailty is identical regardless of treatment assignment. This yields a
proper causal estimand under the Rubin framework. However, the Cox model
fit to stacked data targets a population-averaged effect and may deviate
from the conditional (individual-level) hazard ratio when the
proportional hazards assumption does not hold exactly.

#### Package Variables

| Variable                   | Description                                  |
|:---------------------------|:---------------------------------------------|
| `hr_causal` / `hr_overall` | Overall marginal HR                          |
| `hr_H_true`                | Marginal HR in harm subgroup $`\mathcal{H}`$ |
| `hr_Hc_true`               | Marginal HR in complement $`\mathcal{H}^c`$  |

### Average Hazard Ratio (AHR)

#### Definition

The **average hazard ratio** (AHR) for a subgroup $`\mathcal{S}`$ is
defined directly from the individual-level potential-outcome log hazard
ratios:

``` math

\text{AHR}(\mathcal{S})
= \exp\left(\frac{1}{|\mathcal{S}|} \sum_{i \in \mathcal{S}}
\left[\theta_i(1) - \theta_i(0)\right]\right)
= \exp\left(\frac{1}{|\mathcal{S}|} \sum_{i \in \mathcal{S}}
\log \text{HR}_{i}^{(\text{po})}\right)
```

where $`\log \text{HR}_{i}^{(\text{po})} = \theta_i(1) - \theta_i(0)`$
is the individual causal log hazard ratio.

#### Computation in ForestSearch

The AHR is computed from the `loghr_po` column in the super-population
data frame:

1.  For each subject, compute potential-outcome log-hazards:
    - $`\theta_i(0) = \mathbf{X}_i(0)^{\top} \boldsymbol{\beta}_0`$
    - $`\theta_i(1) = \mathbf{X}_i(1)^{\top} \boldsymbol{\beta}_0`$
    - $`\text{loghr\_po}_i = \theta_i(1) - \theta_i(0)`$
2.  Aggregate over the subgroup:
    $`\text{AHR}(\mathcal{S}) = \exp\left(\overline{\text{loghr\_po}}_{\mathcal{S}}\right)`$

#### Key Properties

The AHR is a deterministic function of the model coefficients and
covariate distribution—it does not depend on the random error
$`\varepsilon_i`$. This makes it fully reproducible across simulation
replications (given the same super-population). Because it operates on
the log-hazard scale and averages individual effects, it captures
heterogeneity in treatment effect across the covariate distribution. For
a subgroup defined by a biomarker threshold $`z \geq c`$, the AHR curve
$`\text{AHR}(c)`$ provides a smooth summary of how the treatment effect
changes with the threshold.

#### Package Variables

| Variable | Description |
|:---|:---|
| `loghr_po` | Individual log HR: $`\theta_i(1) - \theta_i(0)`$ |
| `AHR` | Overall average hazard ratio |
| `AHR_H_true` / `AHR_harm` | AHR in harm subgroup |
| `AHR_Hc_true` / `AHR_no_harm` | AHR in complement |

### Controlled Direct Effect (CDE)

#### Definition

The **controlled direct effect** (CDE) hazard ratio for a subgroup
$`\mathcal{S}`$ is defined as the ratio of average exponentiated
log-hazards under treatment versus control:

``` math

\text{CDE}(\mathcal{S})
= \frac{\frac{1}{|\mathcal{S}|}
\sum_{i \in \mathcal{S}} \exp\left(\theta_i(1)\right)}
{\frac{1}{|\mathcal{S}|}
\sum_{i \in \mathcal{S}} \exp\left(\theta_i(0)\right)}
= \frac{\overline{\exp(\theta_i(1))}_{\mathcal{S}}}
{\overline{\exp(\theta_i(0))}_{\mathcal{S}}}
```

This is the ratio of average hazard *contributions* (on the natural
scale) under treatment versus control. The term “controlled direct
effect” reflects that each $`\theta_i(a)`$ controls for baseline
covariates through the linear predictor.

#### Computation in ForestSearch

In
[`cox_ahr_cde_analysis()`](https://larry-leon.github.io/forestsearch/reference/cox_ahr_cde_analysis.md)
and
[`plot_subgroup_effects()`](https://larry-leon.github.io/forestsearch/reference/plot_subgroup_effects.md),
the CDE is computed from the `theta_1` and `theta_0` columns:

1.  For each subject, retrieve:
    - $`\theta_i(0) = \mathbf{X}_i(0)^{\top} \boldsymbol{\beta}_0`$
    - $`\theta_i(1) = \mathbf{X}_i(1)^{\top} \boldsymbol{\beta}_0`$
2.  Compute subgroup CDE:
    $`\text{CDE}(\mathcal{S}) = \overline{\exp(\theta(1))}_{\mathcal{S}} / \overline{\exp(\theta(0))}_{\mathcal{S}}`$

#### Key Properties

Like the AHR, the CDE is deterministic given the super-population (no
dependence on $`\varepsilon_i`$). However, by averaging on the hazard
(exponential) scale rather than the log-hazard scale, the CDE gives more
weight to subjects with larger absolute hazard contributions. Jensen’s
inequality implies that in general
$`\text{CDE}(\mathcal{S}) \neq \text{AHR}(\mathcal{S})`$, with the
discrepancy increasing as within-subgroup heterogeneity grows.

#### Package Variables

| Variable                        | Description                             |
|:--------------------------------|:----------------------------------------|
| `theta_0`                       | Log-hazard contribution under control   |
| `theta_1`                       | Log-hazard contribution under treatment |
| `cde_results[["overall"]]`      | Overall CDE                             |
| `cde_results[["recommend"]]`    | CDE in recommended subgroup             |
| `cde_results[["questionable"]]` | CDE in questionable subgroup            |

## Theoretical Relationships

### When All Three Agree

Under a constant treatment effect model (no heterogeneity), the
individual log hazard ratio is identical for all subjects:

``` math

\log \text{HR}_i^{(\text{po})} = \beta_{\text{treat}} \quad
\forall\, i
```

In this case, all three measures coincide:
$`\text{HR}_{\text{marg}} = \text{AHR} = \text{CDE} = \exp(\beta_{\text{treat}})`$.

### Sources of Divergence

When treatment effects are heterogeneous, three sources of divergence
arise.

First, **Jensen’s inequality** separates AHR from CDE:

``` math

\exp\left(\overline{\log \text{HR}_i}\right) \neq
\overline{\exp(\theta_i(1))} / \overline{\exp(\theta_i(0))}
```

unless all $`\log \text{HR}_i`$ are identical.

Second, the marginal HR from the stacked Cox model can differ from the
AHR because the Cox partial likelihood implicitly weights subjects by
their risk-set contributions, creating a form of **informative
censoring** in the stacked framework.

Third, the marginal HR depends on the realized $`\varepsilon_i`$ draw,
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

Given a fitted AFT model with parameters
$`(\mu, \sigma, \boldsymbol{\gamma})`$:

**Step 1: Hazard-scale transformation**

``` math

\boldsymbol{\beta}_0 = -\boldsymbol{\gamma} / \sigma
```

**Step 2: Build potential-outcome design matrices**

For each subject $`i`$, construct $`\mathbf{X}_i(1)`$ (treatment arm)
and $`\mathbf{X}_i(0)`$ (control arm). The treatment indicator is set to
1 or 0 respectively; all other covariates remain at their observed
values. Interaction terms (e.g., $`A \times H`$) are computed
accordingly.

**Step 3: Individual log-hazard contributions**

``` math

\theta_i(a) = \mathbf{X}_i(a)^{\top} \boldsymbol{\beta}_0,
\quad a \in \{0, 1\}
```

``` math

\log \text{HR}_i^{(\text{po})} = \theta_i(1) - \theta_i(0)
```

**Step 4: Treatment effect summaries** for subgroup $`\mathcal{S}`$

``` math

\text{AHR}(\mathcal{S})
= \exp\left(\frac{1}{|\mathcal{S}|} \sum_{i \in \mathcal{S}}
\log \text{HR}_i^{(\text{po})}\right)
```

``` math

\text{CDE}(\mathcal{S})
= \frac{\sum_{i \in \mathcal{S}} \exp(\theta_i(1))}
{\sum_{i \in \mathcal{S}} \exp(\theta_i(0))}
```

``` math

\text{HR}_{\text{marg}}(\mathcal{S}):
\text{Cox model on stacked potential outcomes}
```

### GBSG Example: Harm Subgroup

For the GBSG-based DGM with alternative hypothesis, the harm subgroup is
$`\mathcal{H} = \{z_1 = 1 \;\text{AND}\; z_3 = 1\}`$ (low estrogen
receptor AND premenopausal). The theoretical hazard ratios from the AFT
coefficients are:

``` math

\text{HR}(\mathcal{H}) = \exp\left(\beta_0[\text{treat}] +
\beta_0[\text{zh}]\right)
```

``` math

\text{HR}(\mathcal{H}^c) = \exp\left(\beta_0[\text{treat}]\right)
```

The `k_inter` parameter modulates the interaction coefficient:
$`\gamma[\text{zh}] \leftarrow k_{\text{inter}} \times \gamma[\text{zh}]`$,
providing direct control over the magnitude of treatment effect
heterogeneity.

## Summary Comparison Table

The following `gt` table provides a side-by-side comparison of all three
treatment effect measures.

| **Comparison of Treatment Effect Definitions** |  |  |  |
|----|----|----|----|
| Marginal HR, Average Hazard Ratio, and Controlled Direct Effect |  |  |  |
|  | **Marginal HR** | **AHR** | **CDE** |
| **Definition** |  |  |  |
| Formal name | Marginal (causal) hazard ratio | Average hazard ratio | Controlled direct effect (hazard ratio) |
| Notation | HR_marg(S) | AHR(S) | CDE(S) |
| Formula (subgroup S) | exp(beta_Cox) from stacked PO Cox model | exp( mean(loghr_po_i) ) for i in S | mean(exp(theta_1_i)) / mean(exp(theta_0_i)) |
| **Statistical Properties** |  |  |  |
| Averaging scale | Hazard ratio (Cox partial likelihood) | Log-hazard (then exponentiated) | Hazard (natural scale ratio) |
| Depends on epsilon_i | Yes | No | No |
| Reproducible across seeds | No (Monte Carlo variability) | Yes (deterministic) | Yes (deterministic) |
| Accounts for covariate heterogeneity | Averaged out (population-level) | Geometric mean preserves individual effects | Ratio of mean hazards preserves scale |
| Sensitive to outlier hazard contributions | No (rank-based) | Moderate (log scale dampens extremes) | Yes (exponential scale amplifies large values) |
| **Computation** |  |  |  |
| Standard estimation method | Cox PH on stacked potential outcomes | Arithmetic mean of individual log HRs | Ratio of average exponentiated log-hazards |
| Primary ForestSearch use | Validation of DGM causal structure | Primary estimand for subgroup characterization | Complementary biomarker-threshold analysis |
| Key package variable(s) | hr_causal, hr_H_true, hr_Hc_true | loghr_po, AHR, AHR_H_true, AHR_Hc_true | theta_0, theta_1, cde_results |
| **Interpretation & Use** |  |  |  |
| Key package function(s) | calculate_hazard_ratios(), create_gbsg_dgm() | cox_ahr_cde_analysis(), plot_subgroup_effects() | cox_ahr_cde_analysis(), plot_subgroup_effects() |
| Interpretation | Population-average HR as estimated by Cox regression on potential outcomes | Geometric mean of individual causal hazard ratios across subjects | Ratio of average hazard contributions under treatment vs. control |
| Based on the potential outcomes framework within the ForestSearch AFT data-generating mechanism. See León et al. (2024), *Statistics in Medicine*. |  |  |  |

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
    "HR_marg(S)",
    "exp(beta_Cox) from stacked PO Cox model",
    "Hazard ratio (Cox partial likelihood)",
    "Yes",
    "No (Monte Carlo variability)",
    "Averaged out (population-level)",
    "No (rank-based)",
    "Cox PH on stacked potential outcomes",
    "Validation of DGM causal structure",
    "hr_causal, hr_H_true, hr_Hc_true",
    "calculate_hazard_ratios(), create_gbsg_dgm()",
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
    "Primary estimand for subgroup characterization",
    "loghr_po, AHR, AHR_H_true, AHR_Hc_true",
    "cox_ahr_cde_analysis(), plot_subgroup_effects()",
    "Geometric mean of individual causal hazard ratios across subjects"
  ),
  CDE = c(
    "Controlled direct effect (hazard ratio)",
    "CDE(S)",
    "mean(exp(theta_1_i)) / mean(exp(theta_0_i))",
    "Hazard (natural scale ratio)",
    "No",
    "Yes (deterministic)",
    "Ratio of mean hazards preserves scale",
    "Yes (exponential scale amplifies large values)",
    "Ratio of average exponentiated log-hazards",
    "Complementary biomarker-threshold analysis",
    "theta_0, theta_1, cde_results",
    "cox_ahr_cde_analysis(), plot_subgroup_effects()",
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

| **Treatment Effect Computation Pipeline** |  |  |
|----|----|----|
| From AFT model fit to subgroup-level estimands |  |  |
| **Step** | **Mathematical Formula** | **R Implementation** |
| 1\. AFT model | log(T_i) = mu + X_i' gamma + sigma \* epsilon_i | survreg(Surv(y, event) ~ ..., dist = 'weibull') |
| 2\. Hazard-scale transform | beta_0 = -gamma / sigma | b0 \<- -gamma / tau |
| 3\. Individual log-hazards | theta_i(a) = X_i(a)' beta_0 | theta_0 \<- X_control %\*% b0; theta_1 \<- X_treat %\*% b0 |
| 4\. Individual causal log HR | loghr_po_i = theta_i(1) - theta_i(0) | loghr_po \<- theta_1 - theta_0 |
| 5a. AHR (subgroup S) | AHR(S) = exp( mean(loghr_po_i, i in S) ) | exp(mean(df\$loghr_po\[idx\])) |
| 5b. CDE (subgroup S) | CDE(S) = mean(exp(theta_1_i)) / mean(exp(theta_0_i)) | mean(exp(df\$theta_1\[idx\])) / mean(exp(df\$theta_0\[idx\])) |
| 5c. Marginal HR (subgroup S) | HR_marg(S): exp(beta) from coxph() on stacked PO data | exp(coxph(Surv(time,event)~treat, data=df_stacked)\$coef) |
| 6\. Theoretical HR (harm subgroup) | HR(H) = exp( beta_0\[treat\] + beta_0\[zh\] ) | exp(b0\['treat'\] + b0\['zh'\]) |

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
    "6. Theoretical HR (harm subgroup)"
  ),
  Formula = c(
    "log(T_i) = mu + X_i' gamma + sigma * epsilon_i",
    "beta_0 = -gamma / sigma",
    "theta_i(a) = X_i(a)' beta_0",
    "loghr_po_i = theta_i(1) - theta_i(0)",
    "AHR(S) = exp( mean(loghr_po_i, i in S) )",
    "CDE(S) = mean(exp(theta_1_i)) / mean(exp(theta_0_i))",
    "HR_marg(S): exp(beta) from coxph() on stacked PO data",
    "HR(H) = exp( beta_0[treat] + beta_0[zh] )"
  ),
  ForestSearch_Code = c(
    "survreg(Surv(y, event) ~ ..., dist = 'weibull')",
    "b0 <- -gamma / tau",
    "theta_0 <- X_control %*% b0; theta_1 <- X_treat %*% b0",
    "loghr_po <- theta_1 - theta_0",
    "exp(mean(df$loghr_po[idx]))",
    "mean(exp(df$theta_1[idx])) / mean(exp(df$theta_0[idx]))",
    "exp(coxph(Surv(time,event)~treat, data=df_stacked)$coef)",
    "exp(b0['treat'] + b0['zh'])"
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
dgm <- create_gbsg_dgm(
  model   = "alt",
  k_treat = 1.0,
  k_inter = 2.0,
  verbose = TRUE
)

# --- AHR (deterministic, from loghr_po) ---
df <- dgm$df_super_rand
ahr_overall <- exp(mean(df$loghr_po))
ahr_harm    <- exp(mean(df$loghr_po[df$flag.harm == 1]))
ahr_noharm  <- exp(mean(df$loghr_po[df$flag.harm == 0]))

cat("AHR (overall):", round(ahr_overall, 4), "\n")
cat("AHR (harm):   ", round(ahr_harm, 4), "\n")
cat("AHR (Hc):     ", round(ahr_noharm, 4), "\n")

# --- CDE (deterministic, from theta_0 and theta_1) ---
cde_overall <- mean(exp(df$theta_1)) / mean(exp(df$theta_0))
cde_harm    <- mean(exp(df$theta_1[df$flag.harm == 1])) /
               mean(exp(df$theta_0[df$flag.harm == 1]))
cde_noharm  <- mean(exp(df$theta_1[df$flag.harm == 0])) /
               mean(exp(df$theta_0[df$flag.harm == 0]))

cat("CDE (overall):", round(cde_overall, 4), "\n")
cat("CDE (harm):   ", round(cde_harm, 4), "\n")
cat("CDE (Hc):     ", round(cde_noharm, 4), "\n")

# --- Marginal HR (stochastic, from stacked Cox model) ---
cat("HR_marg (overall):", round(dgm$hr_causal, 4), "\n")
cat("HR_marg (harm):   ", round(dgm$hr_H_true, 4), "\n")
cat("HR_marg (Hc):     ", round(dgm$hr_Hc_true, 4), "\n")

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
is the primary target when evaluating operating characteristics
(sensitivity, PPV, bias) across simulation replications.

**Use marginal HR** for validation purposes. The stacked-Cox approach
provides a familiar, model-based estimand that confirms the DGM’s causal
structure produces the intended population-level effects. It is also the
natural target for analyses of observed trial data where potential
outcomes are not directly available.

**Use CDE** for complementary analysis when exploring treatment effect
profiles across continuous biomarker thresholds. The CDE’s natural-scale
averaging can reveal patterns that are attenuated on the log scale,
particularly when hazard contributions are highly variable.

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
