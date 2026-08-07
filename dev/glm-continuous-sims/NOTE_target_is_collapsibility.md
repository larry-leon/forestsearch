# The conditional estimand for the continuous pathway

Corrects three claims in `HANDOFF_glm_continuous_simulations.md`. Each was
established by measurement; the scripts are in `design-checks/`.

## 1. The oracle is not 24–32x sharper

The handoff reports a ratio of 24–32 across $M$ and $\rho$ and reads it as a
property of the oracle. It is a property of the draw budgets. Q1 compares
$\mathrm{sd}(o)/\sqrt{10}$ at `nb = 2e6` against $\mathrm{sd}(p)$ at `nb = 2e4`;
the ratio is $\sqrt{10 \cdot 2\times10^{6} / 2\times10^{4}} = \sqrt{1000} = 31.6$
by construction, and the 24–32 spread is noise in two standard deviations each
estimated from ten observations.

At matched budget the ratio is 1.0 for every configuration (0.90, 1.12, 1.02,
0.93, 0.99, 1.08), and $\mathrm{sd}$ scales as $nb^{-0.517}$ against a predicted
$-0.5$.

The oracle's real advantage is cost per draw — one $M$-dimensional normal against
a pipeline replicate that re-scores the family — which is plausible and remains
**unmeasured**. It needs a timing comparison against `forestsearch()`.

This converts a reported property into a design requirement, which is more
useful: **fix the oracle budget as a multiple of $B$.** For the oracle's own
noise to sit an order of magnitude below the finite-$B$ component,
$nb_{\mathrm{oracle}} \ge 100B$; at $B = 20{,}000$ that is $2\times10^{6}$, the
value the script reached by accident.

What the ratio table obscured is the substantive dependence:

| $M$ | $\rho = 0$ | $\rho = 0.7$ |
|---|---|---|
| 2 | 0.5623 | 0.3069 |
| 5 | 1.1625 | 0.6355 |
| 20 | 1.8668 | 1.0214 |

Monotone in $M$, roughly halved by $\rho = 0.7$. This is the trend bearing on the
family-growth aim.

## 2. Collapsibility, not independence

The handoff derives $\beta(\widehat H) = \tau \cdot \mathbb P(H_{\mathrm{true}}
\mid \widehat H)$ as a ratio of rectangle probabilities and states three limits:
correlated covariates need `mvtnorm::pmvnorm`, GRF disjunctions need
inclusion–exclusion, and the form requires conjunction rules.

**All three dissolve.** The operative property is collapsibility of the mean
difference. Under randomisation ($A \perp X$), for *any* region $R$,

$$\beta(R) = \mathbb E[Y \mid A = 1, R] - \mathbb E[Y \mid A = 0, R]
           = \mathbb E[\tau(X) \mid R],$$

because the prognostic part $\mathbb E[m(X) \mid R]$ is common to both arms and
cancels. No independence between covariates, no rectangle shape, no step-shaped
CATE.

Verified at $n = 3\times10^{6}$ with $\mathrm{cor}(X_1, X_2) = 0.68$, a graded
nonlinear CATE and a strong nonlinear prognostic term:

| region | $\mathbb E[\tau \mid R]$ | fitted MD | diff |
|---|---|---|---|
| conjunction | 27.7149 | 27.7281 | 0.0131 |
| disjunction | 22.4899 | 22.4618 | 0.0281 |
| disjunction $\wedge$ conjunction | 21.9552 | 21.9220 | 0.0332 |
| non-rectangle $\{X_1 + X_2 > 1.2\}$ | 23.9313 | 23.9555 | 0.0241 |
| negation | 24.0787 | 24.0845 | 0.0057 |

Two controls. Under correlation the product-of-marginals form gives 0.42857
against a truth of 0.28609 — off by 0.14, not Monte Carlo error. And on the OR
scale the analogous move fails (marginal 2.4226 against mean individual 2.4596),
as non-collapsibility requires. The second is the required negative control:
without it the MD result could be an artifact any effect measure would pass.

**Consequence for the architecture.** `df_super` is a fixed finite population, so
$\beta(\widehat H)$ is an exact finite mean over it. No evaluation frame, no
`eval_seed`, no model fit, and no residual Monte Carlo noise — not "negligible",
zero. The $\theta^{\dagger}$ sanity gate upgrades from an approximate agreement
test to an exact identity.

## 3. The handoff's formula omits the complement effect

Read from `R/generate_glm_dgm.R:414–427` and `R/simulate_from_glm_dgm.R:97–120`,
the continuous DGM is

$$\mu_{1,i} = \mu_{0,i} + k_{\mathrm{treat}}\beta_{\mathrm{treat}}
              + \beta_{\mathrm{inter}}\mathbb 1\{i \in Q\},$$

with $\beta_{\mathrm{treat}}$ constant because the base fit carries no
treatment-by-covariate interactions — `R/generate_glm_dgm.R:310–323`, where
`rhs <- c(treatment_var, factor_vars, continuous_vars)` builds the additive
formula. The subject-level effect is therefore
two-valued: $\delta := k_{\mathrm{treat}}\beta_{\mathrm{treat}}$ outside $Q$, and
$\delta + \beta_{\mathrm{inter}}$ inside. Hence

$$\beta(\widehat H) = \delta + \beta_{\mathrm{inter}}\cdot
                      \mathbb P(Q \mid \widehat H).$$

The handoff's $\tau \cdot \mathbb P(H_{\mathrm{true}} \mid \widehat H)$ sets the
complement effect to zero. It is not zero: `build_actg175_glm_dgm.R` states that
the complement inherits the fitted ACTG175 treatment effect. Measured against a
fitted MD at $n = 4\times10^{5}$:

| rule | fitted | exact | diff | handoff | diff |
|---|---|---|---|---|---|
| exact rule | 46.3264 | 46.5000 | 0.1736 | 40.0000 | 6.3264 |
| too narrow | 46.4176 | 46.5000 | 0.0824 | 40.0000 | 6.4176 |
| too wide | 16.1600 | 16.2555 | 0.0955 | 9.7555 | 6.4045 |
| one cut only | 16.5708 | 16.6188 | 0.0480 | 10.1188 | 6.4520 |
| disjunction | 10.2819 | 10.2956 | 0.0137 | 3.7956 | 6.4863 |
| disjoint from truth | 6.5193 | 6.5000 | 0.0193 | 0.0000 | 6.5193 |

Every handoff difference is $\delta$; every exact difference is Monte Carlo
error. The row the handoff called most important is the one most clearly wrong: a
rule capturing none of the true subgroup has $\beta = \delta$, not $0$.
`check_betaHhat_md.R` passed only because its synthetic DGM had a zero complement
effect — a control that passed because the quantity was easy.

$\delta$ needs no new computation: it is the DGM's own complement effect,
reproduced to `0.00e+00`.

### The numbers above are fixture-specific; the identity is not

**Correction.** The table was measured on a fixture with $\delta = 6.5$,
$\beta_{\mathrm{inter}} = 40$ and $\text{effect}_Q = 46.5$. That is *not* the
DGM the continuous vignette calibrates. Running
`calibrate_glm_interaction()` with the vignette's own settings
(`quarto/simulations/actg175/continuous/actg175_continuous_simulations.qmd`,
target MD$(Q) = 40$, `n_super = 5000`) gives

$$\delta = -26.255236, \qquad
  \beta_{\mathrm{inter}} = +66.255236, \qquad
  \text{effect}_Q = 40.$$

Two consequences, one cosmetic and one not.

- **$\delta$ is negative in the real fixture.** So the handoff's form does not
  run *low*; it runs **high**, by $26.26$. The phrase "low by $\delta$" is an
  artifact of the illustrative fixture's positive $\delta$ and should be read
  as "**differs from the truth by exactly $\delta$, in whichever direction
  $\delta$ points**". The same correction applies to the `design-checks` table
  in `README.md`.
- **The identity itself is untouched and now independently confirmed.** On the
  calibrated fixture, $\text{exact} - \text{handoff} = \delta$ on all six test
  regions to $0$ or $7.1\times10^{-15}$ — floating point, not method error.
  See `verification/acceptance_betaHhat_md.qmd`, criterion 3.

The structural claim of this section — that the handoff's formula omits the
complement effect and is wrong by exactly $\delta$ — stands. Only the
illustrative magnitudes and the direction word were fixture-specific.

## What remains unverified

Everything above is measured against the *law* or against package source read
directly. **Nothing has been run against `forestsearch()`.** The first
implementation step should be small and should measure.
