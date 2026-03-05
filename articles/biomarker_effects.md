# Biomarker Treatment Effects in the Weibull AFT/Cox Framework

## Introduction

This vignette develops the mathematical framework for
**biomarker-dependent treatment effects** used in the ForestSearch
simulation machinery. The central idea is a piecewise-linear spline
model embedded in a Weibull regression that is simultaneously an
accelerated failure time (AFT) model and a Cox proportional hazards
model. This dual representation allows treatment effects to be
*specified* as causal log-hazard-ratios on the Cox scale while survival
times are *generated* on the AFT scale.

Three treatment effect estimands are defined for summarizing
biomarker-modulated heterogeneity:

- the **Average Hazard Ratio (AHR)** — a geometric-mean summary on the
  log-hazard scale,
- the **Controlled Direct Effect (CDE)** — a ratio of average hazards on
  the natural scale (Aalen, Cook & Roysland, 2015),
- the **Marginal (causal) HR** — a population-level Cox estimate from
  stacked potential outcomes.

The first two are deterministic functions of the data-generating
mechanism (DGM); the third serves as a stochastic validation target.
Their theoretical relationships, conditions for agreement and
divergence, and causal standing are developed below.

## The Weibull AFT/Cox Model

### Weibull Distribution Essentials

For the Weibull distribution with shape parameter $\nu > 0$ and scale
parameter $\theta > 0$, the density, survival, hazard, and cumulative
hazard functions are: $$\begin{aligned}
{f(t)} & {= \nu\,\theta^{- \nu}\, t^{\nu - 1}\,\exp\!( - (t/\theta)^{\nu}),} \\
{S(t)} & {= \exp\!( - (t/\theta)^{\nu}),} \\
{\lambda(t)} & {= \nu\,\theta^{- \nu}\, t^{\nu - 1},} \\
{\Lambda(t)} & {= (t/\theta)^{\nu}.}
\end{aligned}$$

For $\nu \in (0,1)$ the hazard is strictly decreasing; for $\nu > 1$ it
is strictly increasing. A useful distributional identity is that the
cumulative hazard evaluated at the random variable $T$ is
unit-exponential:

**Probability-Integral Property.**
$$\Lambda(T) = (T/\theta)^{\nu}\; \sim \;{Exp}(1).$$ This follows
because
$\Pr\left( \Lambda(T) \leq t \right) = \Pr\left( T \leq \Lambda^{- 1}(t) \right) = F\left( \Lambda^{- 1}(t) \right) = 1 - e^{- t}$.

### The Dual AFT / Cox Representation

Writing $Q = \Lambda(T) \sim {Exp}(1)$ gives
$\log Q = \nu\left( \log T - \log\theta \right)$, which rearranges to
the **AFT form**:

$$\log(T) = \log(\theta) + \tau\,\varepsilon,$$

where $\tau = 1/\nu$ is the AFT scale parameter and
$\varepsilon = \log Q$ has the standard extreme-value distribution with
density $f_{\varepsilon}(x) = \exp\left( x - e^{x} \right)$.

Incorporating a covariate vector $L$, the **Cox (hazard)
parameterization** is:
$$\lambda(t;\, L) = \lambda_{0}(t)\,\exp(L\prime\beta),$$ and the
corresponding **AFT parameterization** is:
$$\log(T) = \log(\theta) + L\prime\gamma + \tau\,\varepsilon,$$ where
the two coefficient vectors are linked by the **AFT-to-hazard
transformation**:

**AFT / Hazard-Scale Transformation.** $$\beta = - \,\gamma\,/\,\tau.$$
R’s `survreg` estimates
$\left( \log\widehat{\theta},\;\widehat{\gamma},\;\widehat{\tau} \right)$
on the AFT scale; hazard-ratio coefficients $\widehat{\beta}$ are
obtained via (4).

This transformation is the linchpin of the framework: potential outcomes
are *simulated* on the AFT scale (3), while treatment effects are
*interpreted* on the hazard-ratio scale via $\beta$.

## Biomarker-Dependent Treatment Effects

### The Two-Phase Spline Model

To induce treatment effects that vary with a continuous biomarker $Z$,
the Cox linear predictor is specified as a piecewise-linear spline with
a single interior knot at $Z = k$:

$$L\prime\beta\;:=\;\beta_{1}\, X + \beta_{2}\, Z + \beta_{3}\, Z\! X + \beta_{4}\,(Z - k)\,\mathbf{1}(Z > k) + \beta_{5}\,(Z - k)\,\mathbf{1}(Z > k)\, X,$$

where $X \in \{ 0,1\}$ denotes the treatment indicator. The five terms
have the following roles:

|  Parameter  | Term                           | Role                                                                      |
|:-----------:|:-------------------------------|:--------------------------------------------------------------------------|
| $\beta_{1}$ | $X$                            | Main treatment effect (intercept of log-HR)                               |
| $\beta_{2}$ | $Z$                            | Prognostic biomarker effect (both arms)                                   |
| $\beta_{3}$ | $ZX$                           | Biomarker $\times$ treatment interaction (slope of log-HR for $Z \leq k$) |
| $\beta_{4}$ | $(Z - k)\mathbf{1}(Z > k)$     | Spline term: change in prognostic slope above $k$                         |
| $\beta_{5}$ | $(Z - k)\mathbf{1}(Z > k)\, X$ | Spline term: change in treatment-effect slope above $k$                   |

### The Causal Log-Hazard-Ratio Function

Under the potential outcomes framework, let $\lambda_{x,z}(t)$ denote
the hazard function for a subject with biomarker level $Z = z$ had they
followed treatment $X = x$. The **causal log-hazard-ratio** at biomarker
level $z$ is:

$$\psi^{0}(z)\;:=\;\log\!(\lambda_{1,z}(t)\,/\,\lambda_{0,z}(t))\; = \;\beta_{1}^{0} + \beta_{3}^{0}\, z + \beta_{5}^{0}\,(z - k)\,\mathbf{1}(z > k).$$

Because the baseline hazard $\lambda_{0}(t)$ cancels in the ratio,
$\psi^{0}(z)$ is free of $t$ — a direct consequence of the proportional
hazards structure conditional on $Z$. Expanding (6) piecewise:

$$\psi^{0}(z) = \begin{cases}
{\beta_{1}^{0} + \beta_{3}^{0}\, z,} & {z \leq k,} \\
{\beta_{1}^{0} + \beta_{3}^{0}\, z + \beta_{5}^{0}\,(z - k),} & {z > k.}
\end{cases}$$

### Anchor-Point Parameterization

Rather than specifying
$\left( \beta_{1}^{0},\beta_{3}^{0},\beta_{5}^{0} \right)$ directly, the
framework specifies the log-HR at three interpretable **anchor points**
— $\psi^{0}(0)$, $\psi^{0}(k)$, and $\psi^{0}(\zeta)$ for some
$\zeta > k$ — and solves:

$$\beta_{1}^{0} = \psi^{0}(0),\qquad\beta_{3}^{0} = \frac{\psi^{0}(k) - \psi^{0}(0)}{k},\qquad\beta_{5}^{0} = \frac{\psi^{0}(\zeta) - \psi^{0}(0) - \beta_{3}^{0}\,\zeta}{\zeta - k}.$$

**Verification.** These formulas are verified by substitution. For
instance, evaluating $\psi^{0}(\zeta)$ from (6) gives
$\psi^{0}(\zeta) = \beta_{1}^{0} + \beta_{3}^{0}\,\zeta + \beta_{5}^{0}(\zeta - k)$;
solving for $\beta_{5}^{0}$ recovers the third expression in (7)
exactly. The same check holds at $z = 0$ and $z = k$, confirming
internal consistency.

### Connection to Potential-Outcome Means

There is an important link between the log-HR $\psi^{0}(z)$ and the
conditional means of the potential log-survival times. Define
$m(x,z):=E\left\lbrack \log T \mid \text{Treat} = x,\, Z = z \right\rbrack$
under the AFT model (3). Then:

$$\psi^{0}(z) = \frac{m(0,z) - m(1,z)}{\tau}.$$

**Derivation.** On the AFT scale the treatment-related terms in
$\log(T)$ contribute $L\prime\gamma = - \tau(L\prime\beta)$ to $m(x,z)$.
Hence the difference in conditional means between control ($x = 0$) and
treatment ($x = 1$) at biomarker level $z$ is:

$$m(1,z) - m(0,z) = \gamma_{1} + \gamma_{3}\, z + \gamma_{5}\,(z - k)\,\mathbf{1}(z > k) = - \tau\,\psi^{0}(z),$$

which gives $\psi^{0}(z) = \{ m(0,z) - m(1,z)\}/\tau$ as required. This
confirms that $\psi^{0}(z)$ measures the *causal difference* in expected
log-survival between the two treatment arms, scaled by $\tau$.

**Sign Convention.** Throughout this vignette, $\psi^{0}(z) < 0$ means
treatment is **beneficial** (HR $< 1$), while $\psi^{0}(z) > 0$ means
treatment is **detrimental** (HR $> 1$). This follows the standard Cox
model convention where negative log-HR favors the experimental arm.

### Extension to Prognostic Factors

When an additional baseline prognostic factor $W$ is included, the
linear predictor becomes:

$$l_{x,z,w} = \beta_{1}^{0}\, x + \beta_{2}^{0}\, z + \beta_{3}^{0}\, zx + \beta_{4}^{0}\,(z - k)\,\mathbf{1}(z > k) + \beta_{5}^{0}\,(z - k)\,\mathbf{1}(z > k)\, x + \beta_{w}^{0}\, w.$$

The prognostic term $\beta_{w}^{0}\, w$ does not interact with
treatment, so the causal log-HR remains:
$$\psi^{0}(z,w) = \frac{m(0,z,w) - m(1,z,w)}{\tau} = \psi^{0}(z),$$
i.e., it depends on $z$ but not on $w$. When prognostic factors *do*
interact with treatment, $\psi^{0}$ would depend on $(z,w)$ jointly, and
the AHR/CDE definitions below would average over both.

### Working Example

A representative configuration for an oncology biomarker setting:

| Anchor | $z$  |       $\psi^{0}(z)$        |     HR | Interpretation       |
|:------:|:----:|:--------------------------:|-------:|:---------------------|
|  Low   | $0$  |   $\log(3) \approx 1.10$   | $3.00$ | Strongly detrimental |
|  Knot  | $5$  | $\log(1.25) \approx 0.22$  | $1.25$ | Modestly detrimental |
|  High  | $10$ | $\log(0.5) \approx - 0.69$ | $0.50$ | Strongly beneficial  |

This profile produces a biomarker-modulated treatment effect that
transitions from harm at low biomarker levels to substantial benefit at
high levels, crossing the null (HR $= 1$) between $z = 5$ and $z = 10$.
The overall average hazard ratio across the full biomarker distribution
is approximately $0.74$, indicating net benefit in the population.

## Treatment Effect Estimands

Three estimands are used to summarize treatment effects across
biomarker-defined sub-populations. Each captures a different aspect of
treatment effect heterogeneity.

### Average Hazard Ratio (AHR)

The **biomarker AHR** is the exponentiated mean of $\psi^{0}(Z)$ across
a biomarker-threshold sub-population:

$$\text{AHR}\left( z^{+} \right):=\exp\!\{ E_{Z \geq z}\;\psi^{0}(Z)\},\qquad\text{AHR}\left( z^{-} \right):=\exp\!\{ E_{Z \leq z}\;\psi^{0}(Z)\}.$$

The AHR is the **geometric mean** of the individual causal hazard ratios
$\exp\left( \psi^{0}\left( Z_{i} \right) \right)$ across the
sub-population. In ForestSearch’s super-population data, this is
computed from the `loghr_po` column as
$\text{AHR}(\mathcal{S}) = \exp\!({\overline{\text{loghr\_po}}}_{\mathcal{S}})$.

**Key properties.**

- *Deterministic*: depends only on the model coefficients and covariate
  distribution, not on the random error $\varepsilon_{i}$.
- *Reproducible*: identical across simulation replications for a fixed
  super-population.
- *Dual causal status*: as a geometric-mean ratio of potential-outcome
  hazards, the AHR achieves both individual-level and population-level
  causal interpretability under the Fay & Li (2024) taxonomy. This is
  because $\log$ linearizes the ratio, so
  $E\left\lbrack \log\left( \text{HR}_{i} \right) \right\rbrack = E\left\lbrack \log\lambda_{1,Z} \right\rbrack - E\left\lbrack \log\lambda_{0,Z} \right\rbrack$,
  and the order of comparison and summarization is irrelevant.

### Controlled Direct Effect (CDE)

The **CDE** is the ratio of average exponentiated log-hazards under
treatment versus control (Aalen, Cook & Roysland, 2015):

$$\text{CDE}\left( z^{+} \right):=\frac{{\bar{\theta}}^{1}(z + )}{{\bar{\theta}}^{0}(z + )},\qquad\text{CDE}\left( z^{-} \right):=\frac{{\bar{\theta}}^{1}(z - )}{{\bar{\theta}}^{0}(z - )},$$

where
${\bar{\theta}}^{x}(z + ) = E_{Z \geq z}\,\exp\left( l_{x,Z} \right)$
and $l_{x,z}$ is the Cox linear predictor at treatment $x$ and biomarker
$z$. In ForestSearch, this is computed from the `theta_1` and `theta_0`
columns:
$\text{CDE}(\mathcal{S}) = {\overline{\exp\left( \theta_{i}(1) \right)}}_{\mathcal{S}}\;/\;{\overline{\exp\left( \theta_{i}(0) \right)}}_{\mathcal{S}}$.

**Key properties.**

- *Deterministic* (like the AHR): no dependence on $\varepsilon_{i}$.
- *Natural-scale averaging*: by averaging on the hazard (exponential)
  scale rather than the log-hazard scale, the CDE gives more weight to
  subjects with larger absolute hazard contributions.
- *Differs from AHR under heterogeneity*: Jensen’s inequality implies
  $\text{CDE}(\mathcal{S}) \neq \text{AHR}(\mathcal{S})$ in general,
  with the discrepancy growing as within-subgroup treatment-effect
  variability increases.

### Marginal (Causal) Hazard Ratio

The **marginal HR** for a subgroup $\mathcal{S}$ is obtained by fitting
a Cox model to *stacked potential outcomes* — a dataset where each
subject contributes two rows (one under treatment, one under control),
both with event indicator $1$:

$$\text{HR}_{\text{marg}}(\mathcal{S}) = \exp\left( {\widehat{\beta}}_{\text{Cox}} \right),$$

where ${\widehat{\beta}}_{\text{Cox}}$ is the coefficient from
`coxph(Surv(time, event) ~ treat)` fit to the stacked data within
$\mathcal{S}$.

**Key properties.**

- *Stochastic*: depends on the realized error terms $\varepsilon_{i}$,
  introducing Monte Carlo variability across simulation replications.
- *Population-level*: the Cox partial likelihood implicitly weights
  subjects by their risk-set contributions, targeting a
  population-averaged effect.
- *Validation role*: serves as the DGM-level ground truth against which
  per-replicate estimator bias is assessed.

### When the Three Estimands Agree — and When They Diverge

Under a **constant treatment effect** (no heterogeneity), the individual
log-HR is identical for all subjects, and all three measures coincide:
$$\text{HR}_{\text{marg}} = \text{AHR} = \text{CDE} = \exp\left( \beta_{\text{treat}} \right).$$

Under **heterogeneous treatment effects**, three sources of divergence
arise:

1.  **Jensen’s inequality** separates AHR from CDE:
    $\exp\!(\overline{\log\text{HR}_{i}}) \neq \overline{\exp\left( \theta_{i}(1) \right)}\;/\;\overline{\exp\left( \theta_{i}(0) \right)}$
    unless all $\log\text{HR}_{i}$ are identical.

2.  **Cox partial-likelihood weighting** separates the marginal HR from
    the AHR, because risk-set membership creates implicit non-uniform
    weights.

3.  **Monte Carlo variability** affects the marginal HR (which depends
    on $\varepsilon_{i}$) but not the AHR or CDE.

The AHR is the recommended **primary estimand** for simulation studies
because it is deterministic, directly interpretable as a geometric-mean
HR, and achieves dual causal status. The marginal HR serves as a
validation target, and the CDE provides a complementary natural-scale
perspective.

## Causal Foundations

### Why the AFT Model Serves as the Causal DGM

Both Aalen, Cook & Roysland (2015) and Fay & Li (2024) identify the AFT
model as a natural causal foundation for survival analysis. The Weibull
AFT scale-change parameter $\exp\left( - \beta_{\text{treat}} \right)$
satisfies:
$$\exp\left( - \beta_{\text{treat}} \right) = \gamma_{1}/\gamma_{0},$$
where $\gamma_{x}$ is the geometric mean of the potential survival time
$T_{i}(x)$. Crucially, this geometric-mean ratio is simultaneously:

- **Individual-level causal**: it equals
  $\exp\!(E\left\lbrack \log T_{i}(1) - \log T_{i}(0) \right\rbrack)$,
  with comparison formed *within* individuals before averaging.
- **Population-level causal**: it is identifiable from marginal
  distributions in a randomized trial.

This dual status — unique among common survival ratio estimands — is why
the ForestSearch DGM is built on the Weibull AFT. The biomarker log-HR
$\psi^{0}(z) = \{ m(0,z) - m(1,z)\}/\tau$ inherits this property at each
biomarker level.

### The AHR as a Causally Valid Functional

The ForestSearch AHR is defined from individual potential-outcome log-HR
differences:
$$\text{AHR}(\mathcal{S}) = \exp\!\left( \frac{1}{|\mathcal{S}|}\sum\limits_{i \in \mathcal{S}}\left\lbrack \theta_{i}(1) - \theta_{i}(0) \right\rbrack \right),$$
where $\theta_{i}(x) = \mathbf{X}_{i}(x)\prime{\mathbf{β}}_{0}$. Because
this is a geometric mean of individual causal log-hazard differences, it
falls within the class of estimands that achieves both individual-level
and population-level causal status under the Fay & Li (2024) taxonomy.

This property is not shared by the instantaneous Cox HR
$e^{\widehat{\beta}}$, which suffers from conditioning-set distortion at
each event time (the “collider” argument of Hernan 2010; Aalen et
al. 2015; Martinussen 2022): among survivors at time $t > 0$, the
treatment and frailty are no longer independent even in a randomized
trial.

### Proportional Hazards and Marginal Non-PH

The Weibull DGM generates data satisfying proportional hazards
*conditional on the biomarker*. Under conditional PH, the Cox HR is a
valid population-level causal estimand (Fay & Li 2024; Beyersmann,
Schmoor & Schumacher 2025). However, *marginally* — averaging over the
biomarker distribution — PH need not hold when treatment effects are
heterogeneous. This is precisely why the framework reports AHR and CDE
(marginal summaries) rather than a single Cox HR as the primary
characterization of biomarker-modulated treatment effects.

The Abrahamowicz et al. (2025) simulations further support this design:
frailty-induced bias in Cox HR estimates is substantial only under very
strong unmeasured susceptibility. Since the ForestSearch DGM has no
unmodeled frailty, the Cox HR within simulations targets the correct
conditional parameter, while the AHR provides the appropriate marginal
summary for real-data applications where unmeasured heterogeneity may
exist.

### Cox HR as an Operational Tool for Subgroup Selection

In the ForestSearch workflow, per-replicate Cox HR estimates and
thresholds (e.g., $\widehat{\beta} \geq \log(0.90)$) serve as
**operational selection criteria** for identifying candidate subgroups.
This is distinct from their role as causal estimands. The recent
literature (Edelmann 2025; Martinussen 2022; Fay & Li 2024) notes that
the instantaneous Cox HR at time $t$ is a causal effect for the
*baseline* population but not necessarily for the population at risk at
$t$.

Within the ForestSearch framework, this distinction is addressed as
follows:

- The *true* treatment effects are defined through $\psi^{0}(z)$, a
  baseline causal quantity.
- Subgroup *selection* uses the Cox HR as an efficient screening tool.
- Subgroup *evaluation* uses the AHR and CDE, which avoid the
  conditioning-set problem.
- Bootstrap bias correction (`hr.H.bc`) and the decomposition of error
  into subgroup misidentification versus finite-sample noise provide
  further safeguards.

## Computation Pipeline

The complete pipeline from model specification to estimand computation
is summarized below.

| Step | Operation                            | Formula                                                                                                                                           |
|:----:|:-------------------------------------|:--------------------------------------------------------------------------------------------------------------------------------------------------|
|  1   | Fit Weibull AFT                      | $\log T_{i} = \mu + \mathbf{X}_{i}\prime{\mathbf{γ}} + \sigma\varepsilon_{i}$                                                                     |
|  2   | Hazard-scale transform               | ${\mathbf{β}}_{0} = - {\mathbf{γ}}/\sigma$                                                                                                        |
|  3   | Individual log-hazards               | $\theta_{i}(x) = \mathbf{X}_{i}(x)\prime{\mathbf{β}}_{0}$                                                                                         |
|  4   | Individual causal log-HR             | $\text{loghr\_po}_{i} = \theta_{i}(1) - \theta_{i}(0)$                                                                                            |
|  5a  | AHR (subgroup $\mathcal{S}$)         | $\text{AHR}(\mathcal{S}) = \exp({\overline{\text{loghr\_po}}}_{\mathcal{S}})$                                                                     |
|  5b  | CDE (subgroup $\mathcal{S}$)         | $\text{CDE}(\mathcal{S}) = {\overline{\exp\left( \theta_{1} \right)}}_{\mathcal{S}}\;/\;{\overline{\exp\left( \theta_{0} \right)}}_{\mathcal{S}}$ |
|  5c  | Marginal HR (subgroup $\mathcal{S}$) | $\exp\left( {\widehat{\beta}}_{\text{Cox}} \right)$ from stacked potential outcomes                                                               |

In ForestSearch, Steps 1–4 are performed by
[`generate_aft_dgm_flex()`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md),
which stores ${\mathbf{β}}_{0}$ in `dgm$model_params$b0` and the
individual-level quantities (`loghr_po`, `lin_pred_0`, `lin_pred_1`) in
`dgm$df_super`. Steps 5a–5b are computed by
[`cox_ahr_cde_analysis()`](https://larry-leon.github.io/forestsearch/reference/cox_ahr_cde_analysis.md);
Step 5c is computed within
[`calculate_hazard_ratios()`](https://larry-leon.github.io/forestsearch/reference/calculate_hazard_ratios.md).

**Consistency Across the ForestSearch Vignettes.** The mathematical
pipeline above — AFT model $\rightarrow$ hazard-scale transformation (4)
$\rightarrow$ individual log-HR $\rightarrow$ AHR/CDE aggregation — is
consistent across the MRCT analysis document, the
`treatment_effect_definitions` vignette, and the
`causal_effects_brief_review` vignette. All three use the same sign
conventions ($\left. \psi^{0}(z) < 0\Leftrightarrow \right.$
beneficial), the same transformation $\beta = - \gamma/\tau$, and the
same averaging operations for AHR and CDE.

The correspondences are:

| This vignette             | `treatment_effect_definitions`                                                                                      | `causal_effects_brief_review`             |
|:--------------------------|:--------------------------------------------------------------------------------------------------------------------|:------------------------------------------|
| $\psi^{0}(z)$             | `loghr_po` $= \theta_{i}(1) - \theta_{i}(0)$                                                                        | Individual-level AFT causal quantity      |
| $\text{AHR}(\mathcal{S})$ | $\exp\left( {\overline{\text{loghr\_po}}}_{\mathcal{S}} \right)$                                                    | Geometric-mean HR with dual causal status |
| $\text{CDE}(\mathcal{S})$ | ${\overline{\exp\left( \theta_{1} \right)}}_{\mathcal{S}}/{\overline{\exp\left( \theta_{0} \right)}}_{\mathcal{S}}$ | Natural-scale complement to AHR           |
| $\beta = - \gamma/\tau$   | ${\mathbf{β}}_{0} = - {\mathbf{γ}}/\sigma$                                                                          | AFT-to-hazard transformation              |

## References

- Aalen, O.O., Cook, R.J. & Roysland, K. (2015). Does Cox analysis of a
  randomized survival study yield a causal treatment effect? *Lifetime
  Data Analysis*, 21(4), 579–593.

- Abrahamowicz, M., Beauchamp, M.-E., Roberts, E.K. & Taylor, J.M.G.
  (2025). Revisiting the hazards of hazard ratios through simulations
  and case studies. *European Journal of Epidemiology*, 40, 611–629.

- Beyersmann, J., Schmoor, C. & Schumacher, M. (2025). Hazards
  constitute key quantities for analyzing, interpreting and
  understanding time-to-event data. *Biometrical Journal*, 67, e70057.

- Edelmann, D. (2025). Revisiting hazard ratios: Can we define causal
  estimands for time-dependent treatment effects? *Biometrical Journal*,
  67, e70100.

- Fay, M.P. & Li, F. (2024). Causal interpretation of the hazard ratio
  in randomized clinical trials. *Clinical Trials*, 21(5), 623–635.

- Hernan, M.A. (2010). The hazards of hazard ratios. *Epidemiology*,
  21(1), 13–15.

- Leon, L.F., Jemielita, T., Guo, Z., Marceau West, R. & Anderson, K.M.
  (2024). Exploratory subgroup identification in the heterogeneous Cox
  model. *Statistics in Medicine*. <doi:10.1002/sim.10163>.

- Martinussen, T. (2022). Causality and the Cox regression model.
  *Annual Review of Statistics and Its Application*, 9, 249–259.

- Prentice, R.L. & Aragaki, A.K. (2022). Intention-to-treat comparisons
  in randomized trials. *Statistical Science*, 37(3), 380–393.
