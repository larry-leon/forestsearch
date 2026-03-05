# Causal Interpretation of Hazard Ratios in Survival Analysis

------------------------------------------------------------------------

> **Purpose.** This document provides a self-contained conceptual review
> of the causal interpretation of hazard ratios in survival analysis,
> drawing on five foundational papers:
>
> | Paper | Primary Contribution |
> |----|----|
> | Aalen, Cook & Røysland (2015) | Foundational critique: Cox HR is not causally interpretable under unmodelled frailty |
> | Martinussen (2022) | Mathematical elaboration; hazard differences equally affected; causal HR defined |
> | Knudsen et al. (2025) | Markov property impossibility under time-varying treatment |
> | Prentice & Aragaki (2022) | Defense of hazard functionals (AHR, RMST) as causally valid in RCTs |
> | **Fay & Li (2024)** | **Individual- vs. population-level estimand taxonomy; CHR as cleanest non-PH causal estimand** |
>
> The first three papers develop the theoretical limitations of the Cox
> HR as a causal estimand; the fourth clarifies what hazard modeling
> *can* validly deliver; the fifth provides a unifying taxonomy that
> reconciles all four. Together they frame the estimand choices
> implemented in ForestSearch.

------------------------------------------------------------------------

## A Foundational Taxonomy: Individual-Level vs. Population-Level Estimands

Before examining what is and is not causally interpretable about the
hazard ratio, Fay & Li (2024) establish a taxonomy that clarifies *why*
the HR causes difficulty — a framework not made explicit in the earlier
literature.

### The Core Distinction

A causal estimand requires comparing two sets of potential outcomes
$`\{Y_i(0) : i \in \Omega\}`$ and $`\{Y_i(1) : i \in \Omega\}`$ over a
common comparison set $`\Omega`$ (a superpopulation satisfying trial
entry criteria). The crucial question is: in what *order* does one apply
the operations of comparison and summarization?

- **Individual-level estimand**: Compare first within individuals, then
  summarize across individuals. The average treatment effect
  $`\Delta_{ATE} = E[Y_i(1) - Y_i(0)]`$ is the canonical example — each
  individual’s causal difference is formed, then averaged.
- **Population-level estimand**: Summarize each arm first, then compare
  the summaries. The ratio of means $`E[Y_i(1)] / E[Y_i(0)]`$ is a
  population-level estimand.

**When these coincide.** Some estimands are both simultaneously. The ATE
for uncensored outcomes is equivalent whether one “compares then
summarizes” or “summarizes then compares,” because expectation is
linear. Similarly, the survival difference $`S_1(t^*) - S_0(t^*)`$ and
the RMST difference are both individual-level and population-level
estimands. This dual status is precisely what makes them the most
*straightforwardly* causal quantities.

**When they diverge.** For ratio estimands the order matters: under
general potential outcome distributions
$`E[Y_i(1)/Y_i(0)] \neq E[Y_i(1)]/E[Y_i(0)]`$. The left-hand
individual-level ratio depends on the unobservable joint distribution of
$`(Y_i(0), Y_i(1))`$ and is not identifiable from a randomized trial.
Only the population-level ratio $`\rho = E[Y_i(1)]/E[Y_i(0)]`$ is
identifiable. Crucially, $`\rho = 2`$ does *not* imply that the
treatment doubles each individual’s outcome on average — it implies only
that the average outcome under treatment is double the average under
control.

**The geometric mean exception.** One ratio estimand that *is* both
individual-level and population-level is the geometric mean ratio for
positive outcomes:
``` math

\exp\!\left(E_F\!\left[\log\frac{Y_i(1)}{Y_i(0)}\right]\right)
= \exp\!\left(E_{F_1}[\log Y_i(1)] - E_{F_0}[\log Y_i(0)]\right)
= \gamma_1/\gamma_0,
```
where $`\gamma_z`$ is the geometric mean of $`Y_i(z)`$. The order of
comparison and summarization is irrelevant here because $`\log`$
linearises the ratio, making the geometric mean ratio identifiable *and*
individually interpretable. This is precisely why the AFT scale-change
parameter has both individual-level and population-level causal status
(see Section 6).

### Classification of Key Survival Estimands

| **Causal Estimand Taxonomy** |  |  |  |  |
|----|----|----|----|----|
| Classification of common survival estimands by Fay & Li (2024) framework |  |  |  |  |
| **Estimand** | **Individual-level?** | **Population-level?** | **Identifiable from RCT?** | **Key Source** |
| Survival difference: S1(t\*) - S0(t\*) | Yes | Yes | Yes | Fay & Li (2024) |
| RMST difference: integral S1(t)dt - integral S0(t)dt | Yes | Yes | Yes | Fay & Li (2024); Prentice & Aragaki (2022) |
| AFT scale-change: exp(-beta) | Yes | Yes | Yes | Fay & Li (2024); Aalen et al. (2015) |
| Cumulative hazard ratio (CHR): log S1(t) / log S0(t) | No | Yes | Yes | Fay & Li (2024); Wei & Schaubel (2008) |
| Cox AHR \[theta_CoxAHR\]: limiting partial likelihood value | No | Yes (under PH) | Yes | Fay & Li (2024) |
| Instantaneous hazard ratio theta(t): lambda1(t)/lambda0(t) | No | Problematic | Yes (but conditioning sets differ) | Aalen (2015); Hernan (2010); Fay & Li (2024) |
| Constant Cox HR theta under PH: lambda1(t)/lambda0(t) = theta | No | Yes | Yes (under PH) | Fay & Li (2024); Prentice & Aragaki (2022) |
| Marginal causal HR \[theta_mcHR\]: conditioning on joint survival | Yes | Yes | No -- requires joint (Y0,Y1) | Martinussen et al. (2020); Fay & Li (2024) |
| Green = both individual- and population-level. Yellow = population-level only. Red = identifiability issues. |  |  |  |  |

------------------------------------------------------------------------

## The Core Problem: Why the Cox HR Is Not Causally Interpretable as an HR

### Setup and Notation

Consider a two-arm randomized trial with survival outcome $`T`$,
treatment indicator $`A \in \{0, 1\}`$, and unmeasured frailty $`Z`$
that influences $`T`$ but is independent of $`A`$ at baseline:
$`Z \perp\!\!\!\perp A`$ at $`t = 0`$.

The **potential outcomes** $`T^0`$ and $`T^1`$ represent what would be
observed under control and treatment respectively. A population causal
effect exists if $`P(T^1 > t)`$ and $`P(T^0 > t)`$ differ. The Cox model
for the marginal hazard is

``` math

\lambda(t \mid A = a) = \lambda_0(t)\, e^{\beta a}.
```

Under randomization $`\lambda_{T^a}(t) = \lambda(t \mid A = a)`$, so Cox
regression yields a **consistent estimator** of $`\beta`$ when the model
is correctly specified. But consistency of estimation does not imply
causal interpretability of the estimand.

### Three Sources of Interpretive Difficulty (Fay & Li 2024)

Fay & Li identify three distinct properties of $`\lambda_z(t)`$ that
complicate its causal interpretation, beyond the collider argument:

1.  **It is a derivative** — a mathematical limit, not directly a
    probability. One can work around this with discrete hazards, but it
    means $`\theta(t)`$ has no simple expectation-based definition.

2.  **It is a local effect** — conveying information only at time $`t`$,
    unlike $`S_z(t)`$ which summarises cumulative experience from 0 to
    $`t`$.

3.  **Its conditioning sets are not comparable.** The hazard at $`t`$ in
    arm $`z`$ conditions on the set
    $`\Omega_z(t) = \{i : Y_i(z) \geq t,\, i \in \Omega\}`$. When the
    treatment has an effect, $`\Omega_0(t) \neq \Omega_1(t)`$ for any
    $`t > 0`$. This violates the identical comparison set requirement
    (equation

    1.  of Fay & Li) that defines a causal estimand.

### The Mechanism: Collider Activation

Even though $`Z \perp\!\!\!\perp A`$ at baseline, among survivors at
time $`t > 0`$

``` math

Z \not{\!\perp\!\!\!\perp} A \mid T \geq t,
```

unless the joint hazard $`h(t, A, Z)`$ is additive in $`A`$ and $`Z`$ —
a condition the Cox multiplicative model violates. Survival $`S_t`$ is a
**collider** on the path $`A \to S_t \leftarrow Z \to S_{t+\Delta}`$.
Conditioning opens a non-causal association between $`A`$ and $`Z`$ that
grows throughout follow-up. Risk sets at every event time beyond the
first are therefore no longer from a randomized comparison.

The joint distribution among survivors is
``` math

P(Z = z, A = a \mid T \geq t) \;\propto\;
\exp\!\left(-\int_0^t h(s, a, z)\, ds\right) P(Z = z)\, P(A = a),
```
which factors only when $`h`$ is additive. Under the Cox multiplicative
structure it does not, so $`e^{\hat\beta}`$ converges to a well-defined
parameter that conflates the true treatment effect with the growing
frailty imbalance between arms. Hernan (2010) termed this “built-in
selection bias.”

### Mathematical Elaboration (Martinussen 2022)

Martinussen locates the root cause in the act of conditioning on
survival at any $`t > 0`$, regardless of model structure. For a
piecewise Cox model with change point $`\nu`$,

``` math

\frac{\lambda(t;\, A=1)}{\lambda(t;\, A=0)} =
\begin{cases} e^{\beta}, & t \leq \nu \\ 1, & t > \nu, \end{cases}
```

the HR declining to 1 after $`\nu`$ is fully consistent with a
beneficial treatment effect at *all* times — including $`t > \nu`$ —
when examined through the frailty-conditional HR $`\text{HR}_Z(t)`$.
These two conclusions are contradictory despite both models being
correctly specified and both estimators consistent. The contradiction
arises purely from selection. Key results:

- A genuinely causal HR conditions on joint survival under both
  potential outcomes:
  ``` math

  \text{HR}(t) =
  \frac{P(t \leq T^1 < t+h \mid T^0 \geq t,\, T^1 \geq t) / h}
       {P(t \leq T^0 < t+h \mid T^0 \geq t,\, T^1 \geq t) / h}.
  ```
  When $`e^{\beta} < 1`$ and the joint distribution follows a frailty
  copula, $`\text{HR}(t) < e^{\beta} < 1`$ — the true causal effect is
  *stronger* than the Cox HR implies. This quantity is not identifiable
  without untestable assumptions. In the Fay & Li taxonomy, this is the
  marginal causal HR ($`\theta_{mcHR}`$): an individual-level and
  population-level estimand, but not identifiable from a randomized
  trial because it requires the joint distribution $`(T^0, T^1)`$.

- **Hazard differences are equally affected.** The Aalen additive
  hazards model produces $`\psi(t) = 0`$ for $`t > \nu`$ for the same
  selection reason — the problem lies in the interpretation of the
  hazard function itself, not in its multiplicative structure.

- The preferred causal quantity is the **survival function contrast**
  $`\delta(t) = P(T^1 > t) - P(T^0 > t)`$, which does not condition on
  survival past $`t = 0`$ and admits efficient doubly-robust estimation
  via the one-step estimator.

### The Vaccine Waning Efficacy Problem (Fay & Li 2024)

Fay & Li’s vaccine trial example (motivated by BNT162b2 COVID-19 data)
makes the practical stakes vivid. They construct two data-generating
models that produce **identical** population cumulative incidence curves
for both trial arms:

- **Model 1** (heterogeneous population, constant individual efficacy):
  95% healthy individuals (hazard ratio: vaccine/placebo = 0.04) and 5%
  frail individuals (same 0.04 ratio). The individual vaccine efficacy
  is constant at 96% for every person, but frail individuals are
  depleted from the placebo arm faster, making the *observed* population
  hazard ratio decline over time.

- **Model 2** (homogeneous population, waning individual efficacy):
  Everyone has genuinely declining protection over time.

Both models generate the same pair of Kaplan-Meier curves. Without
information about individual frailty, the two models cannot be
distinguished from time-to-event data alone. This has profound policy
implications: under Model 1, the research imperative is to identify and
protect frail individuals; under Model 2, it is to develop booster
schedules. A declining observed HR over time does *not* tell us which
world we are in.

**ForestSearch relevance.** This aliasing between population-level
time-varying effects and individual-level heterogeneity is exactly the
problem ForestSearch’s subgroup identification machinery is designed to
address: by conditioning on baseline covariates $`X_i`$, the algorithm
recovers individual-level treatment effect variation that would
otherwise be absorbed into a time-varying population HR. The Fay & Li
example formalises why such conditioning is both theoretically necessary
and practically important.

------------------------------------------------------------------------

## The Cox HR Under and Without Proportional Hazards

### Under Proportional Hazards: A Valid Population-Level Causal Estimand

Fay & Li (2024) identify one setting where the constant Cox HR does
achieve clean population-level causal status. When
$`\lambda_1(t)/\lambda_0(t) = \theta`$ for all $`t`$, the PH assumption
implies $`S_1(t) = S_0(t)^\theta`$, so

``` math

\theta = \frac{\log S_1(t)}{\log S_0(t)} \quad \text{for all } t > 0.
```

This is the **cumulative hazard ratio (CHR)** at any time $`t`$, and it
is a comparison of two survival functions — both population-level causal
estimands from the randomized trial. Under PH the three interpretive
difficulties of $`\theta(t)`$ all vanish: there are no limits, no local
effects, and the conditioning issue is resolved because $`\theta`$
equals $`\log S_1(t)/\log S_0(t)`$ at every $`t`$.

**Important caveat.** Even under PH, $`\theta`$ is *not* an
individual-level causal estimand. A hazard ratio of $`\theta = 0.5`$
does not imply that the treatment halves each individual’s hazard. It
implies that the population hazard rate in the treatment arm is half
that in the control arm at every time $`t`$. The individual-level causal
HR requires untestable assumptions about the joint distribution
$`(T^0, T^1)`$.

### Without Proportional Hazards: The Cumulative Hazard Ratio Is Preferred

When $`\theta(t)`$ is not constant, the various weighted average HRs
(simple AHR, geometric AHR, Cox AHR) all represent different parameters
with different interpretational properties. Fay & Li’s key
recommendation — also endorsed by Wei & Schaubel (2008) — is:

> Under non-proportional hazards, the **cumulative hazard ratio**
> $`\theta_{CHR}(t) = \log S_1(t) / \log S_0(t)`$ is the most
> straightforward causal estimand among those that reduce to $`\theta`$
> under PH.

The CHR achieves population-level causal status at every $`t`$ without
the proportional hazards assumption, because it is a direct contrast of
identifiable survival functions. It removes all three interpretive
difficulties: no limits, no local effects, and both $`S_1(t)`$ and
$`S_0(t)`$ summarise cumulative experience from 0 to $`t`$ on the same
superpopulation $`\Omega`$.

By contrast, the **Cox AHR** $`\theta_{CoxAHR}`$ (the limiting value of
the partial likelihood estimator) is interpretionally complex, depends
on the censoring distribution under non-PH, and lacks the clean causal
derivation the CHR enjoys.

The **loss ratio** $`\theta_{LR}(t)`$, while attractive as a summary,
suffers the same conditioning-set problem as $`\theta(t)`$: its
conditioning sets select different proportions from each arm when there
is a treatment effect.

------------------------------------------------------------------------

## A Practical Defense: What the Cox Model *Can* Legitimately Deliver

### Prentice & Aragaki (2022): The ITT Perspective

Prentice & Aragaki offer a critical counterpoint from decades of
experience with the Women’s Health Initiative (WHI) hormone therapy
trials. They argue that hazard rate contrasts retain valid causal
interpretations in randomized trials — not as instantaneous HR
quantities at a point in time, but as **functionals of hazard
functions** that inherit causal standing directly from randomization.

Hazard rates are the fundamental identifiable quantities in censored
failure time data:
``` math

\lambda(t; z) = \lim_{\Delta\to 0}
\frac{P(T \in [t,t+\Delta) \mid T \geq t,\, A=z)}{\Delta}.
```
Under independent censoring this is consistently estimable via the
Nelson-Aalen estimator, and so is any smooth functional thereof. The
survival function $`F(t; z)`$, the RMST, and the nonparametric **average
hazard ratio** (AHR) all derive a **causal interpretation directly from
randomization** — the independence of $`A`$ from all pre-randomization
risk factors ensures that differences in outcome patterns can be
attributed to the treatment. Prentice & Aragaki explicitly cite
Martinussen et al. (2020) for the view that “various functional
contrasts derived from hazard rates have a causal interpretation using
counterfactuals.”

Fay & Li (2024) provide the formal framework that explains *why*: these
functionals are either individual-level estimands (RMST difference,
survival difference) or population-level estimands whose conditioning
sets are not distorted by differential depletion (CHR, cumulative
survival comparisons).

### Reconciling All Five Papers

The apparent tension between the Aalen/Martinussen/Knudsen critique, the
Prentice & Aragaki defense, and the Fay & Li taxonomy resolves when the
object of inference is precisely specified:

| **Reconciling the Five Papers** |  |  |
|----|----|----|
| What the Cox model can and cannot deliver causally |  |  |
| **Claim** | **Source** | **Resolution** |
| Cox HR e^beta is NOT a causal HR | Aalen (2015); Martinussen (2022); Hernan (2010); Fay & Li (2024) | e^beta consistently estimates a parameter that mixes treatment effect with selection. The conditioning sets Omega0(t) != Omega1(t) violate the common comparison set requirement for causal estimands (Fay & Li). Martinussen shows the true causal HR conditions on joint survival and is stronger than the Cox HR implies. |
| Under PH, constant theta IS a valid pop-level causal estimand | Fay & Li (2024); Prentice & Aragaki (2022) | When PH holds, theta = log S1(t)/log S0(t) -- a comparison of identifiable survival functions -- resolving all three interpretive difficulties. Still not individual-level. |
| CHR = log S1(t)/log S0(t) is causally valid under non-PH | Fay & Li (2024); Wei & Schaubel (2008) | The CHR directly contrasts survival functions without conditioning on survival at intermediate times. Recommended over Cox AHR under non-PH. |
| RMST difference and survival difference are both individual- and population-level | Fay & Li (2024) | Because E\[S1(t) - S0(t)\] = E\[I(Y(1) \> t) - I(Y(0) \> t)\], these estimands compare-then-summarize at the individual level and are also identifiable as population summaries. |
| AFT exp(-beta) is both individual-level and population-level | Fay & Li (2024); Aalen et al. (2015) | The geometric mean log ratio exp(E\[log Y(1)\] - E\[log Y(0)\]) is both individually and population interpretable -- the log transformation makes it collapsible. |
| Time-varying HRs -\> 1 do not imply effect waning | Martinussen (2022); Fay & Li (2024) | A correctly specified piecewise HR declining to 1 reflects growing frailty imbalance, not biological effect waning. The Fay & Li vaccine example shows two indistinguishable worlds. Baseline covariate conditioning (GRF/LASSO) is the correct remedy. |
| Hazard function functionals are causally valid in RCTs | Prentice & Aragaki (2022); Martinussen et al. (2020) | Randomization confers causal validity on hazard function estimates; functionals (AHR, RMST, survival contrasts) inherit this. Cox is a valid estimation engine. |
| Markov property cannot hold under frailty | Knudsen et al. (2025) | Any non-degenerate frailty violates the Markov property required for causal interpretation of Cox MSMs with time-varying treatment. |

The key distinction is between:

- **The Cox partial likelihood parameter** $`e^{\hat\beta}`$ —
  consistently estimated but not interpretable as a causal HR. This is
  what Aalen et al., Martinussen, and Fay & Li critique.
- **Functionals of the estimated hazard function** (CHR, AHR, RMST,
  survival contrasts) — causally valid in randomized trials because
  randomization confers validity on the underlying hazard estimates,
  even when those estimates come from Cox fits. This is what Prentice &
  Aragaki defend, and Fay & Li formalise by classifying each
  functional’s individual- and population-level status.

### What the Cox Model Can Legitimately Do

Combining all five papers, the following uses of Cox regression are
defensible:

1.  **Hypothesis testing** via the log-rank score test: valid under any
    model specification.

2.  **Under PH**: the Cox partial likelihood estimator consistently
    estimates $`\theta`$, a valid population-level causal estimand
    equalling $`\log S_1(t)/\log S_0(t)`$ at every $`t`$.

3.  **CHR estimation**: plugging Kaplan-Meier estimates into
    $`\log \hat{S}_1(t)/\log \hat{S}_0(t)`$ provides the cleanest
    population-level causal summary under non-PH (Fay & Li; Wei &
    Schaubel).

4.  **Flexible AHR and RMST estimation**: Cox-based survival function
    estimates feed into the nonparametric AHR and RMST contrasts that
    are causally valid in randomized trials (Prentice & Aragaki).

5.  **Covariate adjustment**: adding baseline covariates to the Cox
    model moves the estimand progressively closer to an individual-level
    causal quantity. If enough covariates capture frailty variation, the
    estimated conditional HR approaches the individual-level causal HR
    (Fay & Li).

6.  **Descriptive clinical communication**: “persons assigned to CEE+MPA
    had about an 80% increase in CHD incidence over year 1” is a valid
    and useful ITT statement, provided it is not interpreted as a causal
    HR at that specific time point in isolation (Prentice & Aragaki).

------------------------------------------------------------------------

## The AFT Model as a Causal Ground Truth

### Collapsibility and Dual Causal Status

Both Aalen et al. (2015) and Fay & Li (2024) identify the AFT model as
the natural causal alternative to the Cox model, with Fay & Li providing
the clearest formal explanation of *why*. The Weibull AFT can be
written:

``` math

\log T_i = \mu + \mathbf{X}_i^\top \boldsymbol{\gamma} + \sigma\, \varepsilon_i,
\qquad \varepsilon_i \sim \text{Extreme Value}.
```

The scale-change parameter $`\exp(-\beta_{\text{treat}})`$ satisfies
$`S_1(t) = S_0(t\, e^{-\beta})`$. Because failure times are positive,
the treatment effect can be expressed as a ratio as in equation (4) of
Fay & Li:

``` math

\exp(-\beta) = \gamma_1/\gamma_0,
```

where $`\gamma_z`$ is the geometric mean of $`Y_i(z)`$. The geometric
mean ratio is simultaneously individual-level (it equals
$`\exp(E[\log Y_i(1) - \log Y_i(0)])`$) and population-level (it is
identifiable from marginal distributions), making $`\exp(-\beta)`$ the
**only common survival estimand to achieve both simultaneously**.

This collapsibility — invariance to marginalising over the frailty
distribution $`Z`$ — ensures the marginal AFT treatment effect is
consistent regardless of how heterogeneous the population is.
Empirically, Aalen et al. (2015) show AFT estimators are unbiased across
a wide range of censoring rates and frailty variances, while Cox
estimators carry substantial variable bias. Prentice & Aragaki do not
dispute AFT collapsibility, but argue that Cox modeling with flexible HR
specification and careful functional derivation can recover causally
valid summaries in practice. The two positions are not contradictory:
the AFT provides the superior causal DGM, while the Cox model provides a
semiparametrically efficient estimation engine for the AHR and RMST
functionals that flow from it.

This is why ForestSearch builds its DGM on a Weibull AFT.

### AFT-to-Hazard Scale Transformation

Under the Weibull AFT parameterisation, the hazard-scale coefficients
are

``` math

\boldsymbol{\beta}_0 = -\boldsymbol{\gamma} / \sigma,
```

so the individual conditional hazard is
$`h(t \mid \mathbf{X}_i) = h_0(t)\, \exp(\mathbf{X}_i^\top \boldsymbol{\beta}_0)`$.
This links the AFT causal ground truth to the hazard-scale estimands
(AHR, CDE, CHR, marginal HR), enabling ForestSearch to compute all from
AFT-fitted parameters while retaining the individual- and
population-level causal validity of the underlying DGM.

------------------------------------------------------------------------

## The Average Hazard Ratio: A Causally Valid Functional

### Definition and Estimation

The nonparametric AHR of Kalbfleisch & Prentice (1981), prominently
featured in Prentice & Aragaki (2022) and discussed in Fay & Li (2024),
is:

``` math

\widehat{\text{AHR}}(t) =
\left[1 - \hat{F}(t;\, z=0)\right]^{-1}
\int_0^t \hat{F}(s;\, z=0)\, d\hat{\Lambda}(s;\, z=1).
```

This weights each instantaneous hazard ratio by the control-group
failure probability at that moment — a time-averaged summary causally
attributable to randomization group assignment. Fay & Li note that while
this AHR is not a straightforward causal estimand from a
conditioning-set standpoint, it is “at least in one example similar in
magnitude to the ratio of cumulative hazards” — i.e., it tracks the
clean CHR estimand in practice.

### AHR vs. CHR vs. RMST: A Practical Ordering

Fay & Li (2024) and Prentice & Aragaki (2022) together suggest the
following ordering of causal cleanness for non-PH settings:

| Estimand | Causal Status | PH required? | Sensitive to early effects? |
|----|----|----|----|
| RMST difference | Individual & population | No | Lower (rare outcomes) |
| Survival difference at $`t^*`$ | Individual & population | No | Depends on $`t^*`$ |
| CHR $`\log S_1(t)/\log S_0(t)`$ | Population | No | Yes |
| Nonparametric AHR | Population (approximately) | No | Yes |
| Cox HR under PH | Population | Yes (must verify) | Yes |
| Instantaneous HR $`\theta(t)`$ | Problematic | N/A | Yes (but misleading) |

### AHR vs. RMST: Empirical Evidence from WHI

Prentice & Aragaki’s analyses of the CEE+MPA trial provide instructive
real-data evidence for the relative sensitivity of these two estimands.
For **breast cancer**:

- The **AHR** detected a statistically significant elevation by **6
  years** from randomization (AHR $`\approx`$ 1.20, 95% CI 1.00–1.46).
- The **RMST** contrast did not show a clear reduction until
  approximately **12 years**, because the absolute risk difference
  accumulates slowly for a relatively rare outcome.

Crucially, a breast cancer risk elevation was the *principal trigger*
for the early stoppage of the CEE+MPA trial following 5.6 years of
follow-up. Had RMST been the primary reporting tool, the trial-stopping
signal would not have been detectable at that point. This real-data
finding validates the theoretical argument that AHR functionals are more
sensitive to early or moderate effects in prevention-trial settings.

### Relationship to ForestSearch’s AHR Estimand

ForestSearch defines the AHR at the **individual level** from the DGM:

``` math

\text{AHR}(S) = \exp\!\left(\frac{1}{|S|} \sum_{i \in S}
(\theta_i(1) - \theta_i(0))\right),
```

where $`\theta_i(a) = \mathbf{X}_i(a)^\top \boldsymbol{\beta}_0`$ are
individual log-hazards from the AFT DGM. This geometric mean of
individual causal log-hazard differences is — by the Fay & Li framework
— precisely the kind of estimand that achieves both individual-level and
population-level status, because it averages log-ratios of potential
outcomes and is therefore collapsible. The ForestSearch AHR is not
subject to the conditioning-set problem because it is computed directly
from potential outcomes without conditioning on observed survival.

The distinction from the Prentice & Aragaki AHR matters:

- The **Prentice & Aragaki AHR** is a population-level quantity weighted
  by control-group failure time, estimable from observed data.
- The **ForestSearch AHR** is a DGM-level quantity averaging over
  individual potential-outcome log-hazards, providing the simulation
  ground truth for subgroup causal effect.

Both avoid the collider bias problem: the Prentice AHR through
integration over follow-up time, the ForestSearch AHR through direct use
of potential outcomes from the AFT DGM without conditioning on any
observed survival.

------------------------------------------------------------------------

## Extension to Time-Varying Treatment (Knudsen et al. 2025)

### The Markov Property and Its Impossibility

When treatment $`A(t)`$ evolves over time, Cox MSMs require the **Markov
property** — that the marginal hazard depends only on current treatment
$`A(t)`$, not on treatment history $`\overline{A}(t)`$. Knudsen et
al. (2025) prove this requires a degenerate (point-mass) frailty
distribution. With complete-data hazard
$`\lambda(t \mid A(t), Z) = Z\, h(t, A(t))`$,

``` math

\lambda(t \mid \overline{A}(t)) =
h(t, A(t))\, E\!\left(Z \mid \overline{A}(t), T \geq t\right),
```

and the conditional expectation of $`Z`$ depends on the entire
cumulative hazard history
$`H(t, \overline{A}(t)) = \int_0^t h(s, A(s))\, ds`$, not just $`A(t)`$.
For Markovianity, the Laplace transform of $`Z`$ must be exponential —
corresponding to a point-mass distribution. Since unmeasured
heterogeneity (genetic, environmental, behavioural) is universal, the
Markov property **cannot realistically hold in any real dataset**.

### Implications for Cox MSMs and Extended Kaplan-Meier

Marginal structural Cox models estimate $`\beta`$ in
$`\lambda^a(t) = \lambda_0(t)\, e^{\beta a(t)}`$ using IPTW weighting.
Even with perfect confounding control, causal interpretation of
$`e^{\hat\beta}`$ as a survival probability ratio requires the Markov
property. Knudsen et al. (2025) demonstrate with a concrete
illness-death model that the rate-based survival estimate can
substantially underestimate the true causal survival benefit: at
$`t = 3`$ the rate-based estimate is 0.14 vs. the true causal contrast
of 0.22 — a 57% underestimate.

The extended Kaplan-Meier curve (Simon & Makuch 1984) estimates
$`\exp(-\text{cumulative rate})`$, not an intervention-specific survival
probability. Under the Cox rate model, the ratio of log-extended-KM
curves estimates the rate ratio — but the rate ratio has no clear causal
interpretation without the Markov property.

**Connection to Fay & Li.** The Knudsen et al. finding is the
time-varying treatment analogue of the Fay & Li conditioning-set
problem: in both cases, the hazard at $`t`$ conditions on a set of
survivors whose composition has been distorted by differential treatment
effects up to $`t`$. Fay & Li address this for fixed baseline treatment;
Knudsen et al. extend it to the time-varying case where the distortion
is actively controlled by the treatment regime being studied.

### Relevance to ForestSearch

ForestSearch’s DGM uses **fixed baseline treatment assignment**, placing
it in the Fay & Li baseline-treatment setting where these complications
do not arise. No Markov property is invoked, and the potential outcome
approach used to define the ForestSearch AHR and CDE directly avoids the
conditioning-set distortion identified by Fay & Li.

Any extension of ForestSearch to **time-varying treatment subgroup
identification** would need to contend with Knudsen et al.’s results
directly: subgroup-specific Cox MSM HRs would not be causally
interpretable even with perfect IPTW weighting, and the estimand choice
(rate ratio vs. intervention-specific survival contrast) would become
critical.

------------------------------------------------------------------------

## Connections to ForestSearch Estimands

The five papers collectively motivate the specific estimand choices in
ForestSearch.

| **Theory–Practice Mapping** |  |
|----|----|
| From causal inference literature to ForestSearch estimand design |  |
| **Theoretical Issue** | **ForestSearch Design Response** |
| Cox HR e^beta mixes treatment effect and selection bias; conditioning sets are not comparable (Aalen 2015; Martinussen 2022; Fay & Li 2024) | AHR and CDE are computed directly from AFT potential outcomes; no risk-set conditioning involved. Stacked PO marginal HR bypasses partial likelihood risk sets entirely; observed Cox HR treated as biased comparator |
| Under PH, theta = log S1(t)/log S0(t) achieves valid population-level causal status (Fay & Li 2024) | When Cox models are used for AHR/RMST estimation in applied output, proportionality should be assessed and CHR reported as primary summary |
| CHR is cleanest estimand under non-PH; recommended over Cox AHR (Fay & Li 2024; Wei & Schaubel 2008) | ForestSearch applied output should report CHR alongside AHR; the nonparametric AHR tracks CHR closely in practice |
| AFT exp(-beta) is both individual-level and population-level; the only survival ratio estimand with dual status (Fay & Li 2024; Aalen 2015) | Weibull AFT DGM provides ground truth with both individual- and population-level causal validity; beta0 = -gamma/sigma links to all hazard-scale estimands |
| RMST and survival differences are both individual- and population-level causal estimands (Fay & Li 2024) | RMST difference is a natural complement to AHR in ForestSearch output; less sensitive to early effects but cleanly dual-causal |
| Declining observed HR cannot distinguish genuine effect waning from frailty aliasing (Martinussen 2022; Fay & Li 2024) | GRF/LASSO subgroup identification conditions on baseline X_i, which Fay & Li show moves estimation toward individual-level effects and resolves the aliasing |
| Hazard differences share the same causal limitation as HRs (Martinussen 2022) | CDE = mean(exp(theta_i(1))) / mean(exp(theta_i(0))) operates on exponential scale, avoiding log compression; parallels the Prentice & Aragaki AHR weighting structure |
| Adding baseline covariates moves Cox HR toward individual-level causal estimand (Fay & Li 2024) | GRF baseline covariate conditioning is justified by Fay & Li: adding covariates moves the estimated HR toward individual-level causal truth |
| Markov property fails under any non-degenerate frailty in time-varying treatment (Knudsen 2025) | Baseline-only treatment keeps ForestSearch in the Aalen et al. / Fay & Li safe zone; no Markov assumption required |

### Estimand Selection Guide

**Primary simulation estimand: ForestSearch AHR.** Computed from AFT
potential outcomes as the geometric mean individual log-hazard
difference, it is the only hazard-ratio-scale estimand that achieves
both individual-level and population-level causal status (Fay & Li
framework). It is unaffected by conditioning-set distortion,
deterministic given the superpopulation, and more sensitive to early
effects than RMST. From the WHI empirical evidence, this sensitivity
advantage is clinically material: AHR detected the breast cancer risk
elevation at 6 years while RMST required 12.

**Applied reporting: CHR alongside AHR.** For applied output — subgroup
hazard ratio summaries reported to clinical trialists — the cumulative
hazard ratio $`\log \hat{S}_1(t)/\log \hat{S}_0(t)`$ is the most
defensible population-level causal estimand under non-PH (Fay & Li; Wei
& Schaubel). The nonparametric AHR of Prentice & Aragaki tracks this
closely in practice and adds useful time-averaged narrative.

**Complement: RMST difference.** Both individual- and population-level
causal (Fay & Li). Less sensitive than AHR/CHR for rare outcomes over
short follow-up (WHI breast cancer example), but provides clinically
intuitive absolute time-to-event framing that is entirely unambiguous
causally.

**Use the CDE** as a complementary natural-scale perspective. By
averaging $`\exp(\theta_i(1))`$ and $`\exp(\theta_i(0))`$ separately
before taking a ratio, the CDE parallels the Prentice & Aragaki AHR
weighting structure on the hazard scale, and is particularly informative
when hazard contributions are highly variable across subjects.

**Use the marginal HR (stacked PO)** for validation. It provides the
closest DGM-level analogue to what a correctly specified
population-level Cox model would target. Bias in the naive Cox HR
(`hr.H.hat`) relative to `dgm$hr_H_true` quantifies the combined effect
of finite-sample noise, censoring, and the selection bias documented by
Aalen et al.

**GRF baseline covariate conditioning is theoretically justified.** Fay
& Li demonstrate that conditioning on more baseline covariates moves the
Cox HR estimand continuously from a population-level toward an
individual-level causal quantity. ForestSearch’s GRF/LASSO covariate
conditioning is not merely a statistical procedure — it is the mechanism
by which individual-level heterogeneous treatment effects become
estimable and by which the frailty aliasing illustrated in the vaccine
efficacy example is resolved.

**Avoid interpreting per-replicate Cox HR causally.** The `hr.H.hat`
column in simulation results estimates the marginal Cox AHR, which under
non-PH depends on the censoring distribution, is not cleanly
individual-level, and conflates treatment effect with selection. Its
divergence from `dgm$hr_H_true` quantifies the combined selection bias,
finite-sample noise, and censoring distortion.

------------------------------------------------------------------------

## Key Takeaways

1.  **Causal estimands fall into individual-level, population-level, or
    both (Fay & Li 2024).** The RMST difference, survival difference,
    and AFT $`\exp(-\beta)`$ are both simultaneously. The Cox HR is at
    best population-level only (and only under PH).

2.  **The instantaneous Cox HR $`\theta(t)`$ fails the common comparison
    set requirement** of causal inference: when there is a treatment
    effect, $`\Omega_0(t) \neq \Omega_1(t)`$ for $`t > 0`$. This is the
    Fay & Li formalisation of the Aalen et al. collider argument.

3.  **Under proportional hazards, the constant Cox HR achieves valid
    population-level causal status** as $`\log S_1(t)/\log S_0(t)`$.
    This is Fay & Li’s most constructive result and reconciles the
    Aalen/Martinussen critique with the Prentice & Aragaki defense.

4.  **Under non-proportional hazards, the cumulative hazard ratio is the
    preferred causal estimand** (Fay & Li; Wei & Schaubel). It directly
    contrasts identifiable survival functions and retains
    population-level causal validity without a PH assumption.

5.  **The AFT $`\exp(-\beta)`$ is the only common survival ratio
    estimand that is both individual-level and population-level** (Fay &
    Li). This geometric mean property underpins ForestSearch’s choice of
    AFT as the causal DGM.

6.  **A declining observed HR over time cannot distinguish genuine
    effect waning from constant individual-level effects in a
    heterogeneous population.** The Fay & Li vaccine example and
    Martinussen’s piecewise Cox model both demonstrate this aliasing.
    Baseline covariate conditioning (GRF/LASSO in ForestSearch) is the
    correct remedy.

7.  **AHR is more sensitive than RMST to early effects for rare
    outcomes** (Prentice & Aragaki). The WHI breast cancer example —
    where AHR detected risk elevation at 6 years while RMST required 12
    — provides real-data validation that AHR/CHR-scale estimands are the
    right primary tools for subgroup detection in prevention-trial
    settings.

8.  **Hazard differences share the same causal limitation as hazard
    ratios** (Martinussen 2022). The problem lies in interpreting the
    hazard function itself, not in the multiplicative structure of the
    Cox model.

9.  **The Markov property is always violated under frailty** (Knudsen et
    al. 2025), disqualifying Cox MSMs with time-varying treatment from a
    clean causal interpretation. ForestSearch’s baseline-only treatment
    design avoids this problem entirely.

10. **Cox regression is a valid estimation engine but a problematic
    interpretive frame.** All five papers support the same division of
    labour: use Cox models for efficient, flexible estimation of hazard
    functions; communicate causal effects through functionals (CHR, AHR,
    RMST, survival contrasts) rather than as point-in-time HRs.

------------------------------------------------------------------------

## References

- Aalen, O.O., Cook, R.J., & Røysland, K. (2015). Does Cox analysis of a
  randomized survival study yield a causal treatment effect? *Lifetime
  Data Analysis*, 21(4), 579–593.
  <https://doi.org/10.1007/s10985-015-9335-y>

- Fay, M.P., & Li, F. (2024). Causal interpretation of the hazard ratio
  in randomized clinical trials. *Clinical Trials*, 21(5), 623–635.
  <https://doi.org/10.1177/17407745241243308>

- Hernan, M.A. (2010). The hazards of hazard ratios. *Epidemiology*,
  21(1), 13–15.

- Hernan, M.A., & Robins, J.M. (2020). *Causal Inference: What If*. Boca
  Raton: Chapman & Hall/CRC.

- Kalbfleisch, J.D., & Prentice, R.L. (1981). Estimation of the average
  hazard ratio. *Biometrika*, 68(1), 105–112.

- Knudsen, M.B., Gabriel, E.E., Martinussen, T., Rytgaard, H.C.W., &
  Sjolander, A. (2025). On the limitations for causal inference in Cox
  models with time-varying treatment. *arXiv:2504.01524*.

- Martinussen, T. (2022). Causality and the Cox regression model.
  *Annual Review of Statistics and Its Application*, 9, 249–259.
  <https://doi.org/10.1146/annurev-statistics-040320-114441>

- Martinussen, T., Vansteelandt, S., & Andersen, P.K. (2020). Subtleties
  in the interpretation of hazard contrasts. *Lifetime Data Analysis*,
  26(4), 833–855.

- Pearl, J. (2009). *Causality: Models, Reasoning, and Inference* (2nd
  ed.). Cambridge University Press.

- Prentice, R.L., & Aragaki, A.K. (2022). Intention-to-treat comparisons
  in randomized trials. *Statistical Science*, 37(3), 380–393.
  <https://doi.org/10.1214/21-STS830>

- Rubin, D.B. (1974). Estimating causal effects of treatments in
  randomized and nonrandomized studies. *Journal of Educational
  Psychology*, 66(5), 688–701.

- Simon, R., & Makuch, R.W. (1984). A non-parametric graphical
  representation of the relationship between survival and the occurrence
  of an event. *Statistics in Medicine*, 3(1), 35–44.

- Sjolander, A. (2020). A cautionary note on extended Kaplan-Meier
  curves for time-varying covariates. *Epidemiology*, 31(4), 517.

- Thomas, S.J., et al. (2021). Safety and efficacy of the BNT162b2 mRNA
  Covid-19 vaccine through 6 months. *New England Journal of Medicine*,
  385, 1761–1773.

- Uno, H., Wittes, J., Fu, H., et al. (2015). Alternatives to hazard
  ratios for comparing the efficacy or safety of therapies in
  noninferiority studies. *Annals of Internal Medicine*, 163(2),
  127–134.

- Wei, G., & Schaubel, D.E. (2008). Estimating cumulative treatment
  effects in the presence of nonproportional hazards. *Biometrics*,
  64(3), 724–732.
