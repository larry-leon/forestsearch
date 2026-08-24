# The Crump machinery as the bridge from the MD detection theory to binary and survival

Companion script: `check_binary_rd_process.R` (this directory; base R +
mvtnorm, sources the package's `R/test_hte_crump.R` verbatim, seed
20260824, ~5 s, 8 PASS / 0 FAIL on both the drafting sandbox and the
primary machine).

## 1. Question and verdict

Whether the package's Crump (2008) heterogeneity-test machinery
(`R/test_hte_crump.R`; validated in
`validation_hte_tests_crump.qmd`) can be utilized to extend the OLS
analytic tools — the candidate-effect-process detection/selection theory
verified for continuous/MD — to binary and survival outcomes.

**Verdict: yes, in a strong and precise sense.** The Crump machinery is not
an alternative to the candidate-effect process; it is built from the same
mathematical object. Its two computational primitives — per-arm linear or
GLM projections of the outcome, and sandwich (score outer-product)
variances — are exactly the Level 3 plug-in layer of the extended process;
its linear-probability validation is the license for running the binary
extension on the risk-difference scale; and its Poisson-plus-offset branch
is the package's own chosen bridge to survival ("the same model underlying
the IRR pipeline"; validation doc §Application to Survival). The pieces the
Crump machinery does not supply — per-candidate mixture means, region-wise
thresholds, and the selection/argmax layer — are exactly the pieces the MD
theory already provides. The first rung of the bridge (binary/RD, plus the
joint law with the actual package gate) is verified numerically in the
companion script; the remaining rungs are specified with their governing
quantities named.

## 2. What the machinery is, from source

`test_zero_cate` / `test_constant_cate` (`R/test_hte_crump.R`) fit the
outcome on a polynomial basis $P_K(X)$ **separately by arm** —
`regression = "ols"`, `"logistic"`, or `"poisson"` (with `offset`) — form
the arm-coefficient difference $\Delta = \hat\beta_1 - \hat\beta_0$, and
test the quadratic form $Q = \Delta^\top \widehat V^{-1} \Delta$ against
$\chi^2$ (all components for zero CATE; the slope block for constant
CATE), with the Crump eq. 3.12 normalization $(Q - K)/\sqrt{2K}$ also
reported. $\widehat V = \Omega_0/N_0 + \Omega_1/N_1$ with HC0 sandwiches:
OLS residual-based (`sandwich_var_ols`) or GLM score-based
(`sandwich_var_glm`, scores $P_i\,(d\mu/d\eta)(Y_i-\mu_i)/\mathrm{Var}(\mu_i)$).
`test_hte_from_forestsearch` runs the tests on FS-selected covariates as a
pipeline pre-screen (validation doc §Connection: the constant-CATE test as
the formal justification gate before `forestsearch()`).

Two facts from the validation document anchor the extension. For binary
outcomes, the OLS branch is a linear probability model whose target is the
**risk-difference** CATE, with the HC sandwich absorbing the
$p(x)(1-p(x))$ heteroskedasticity automatically (§Application to Binary).
For survival, the Poisson-plus-offset branch models
$\log E[Y_i] = P_K(X_i)^\top\xi_w + \log t_i$ with event indicator $Y$ and
follow-up $t$ — the classical low-rate approximation to Cox, and the model
already underlying the package's IRR pathway (§Application to Survival).

## 3. Three roles in the extension

**Role 1 — the license and vehicle for binary (load-bearing).** Under OLS,
every arm-projection coefficient is linear in $Y$. The within-region
difference in proportions $\hat\beta_{\mathrm{RD}}(g) = \hat p_1(g) -
\hat p_0(g)$ is the same class of object with region-indicator
coefficients $c_i(g)$ — the candidate-effect process **is** the
indicator-basis Crump statistic, one region at a time. Consequently:

- *Means transfer as Tier A.* The risk difference is collapsible:
  $\beta_{\mathrm{RD}}(g) = \mathbb E[\tau_{\mathrm{RD}}(X)\mid g]$ for
  every region, enumerable exactly on the DGM's probability surfaces
  (and, for a two-valued risk model, again a two-point mixture). The
  known Tier A $\to$ Tier B degradation applies only to the OR *reporting*
  scale, not to the RD *detection* vehicle.
- *Covariance transfers as Tier B with one change.* Homoskedastic
  $\sigma^2$ becomes per-subject $\mathrm{Var}(Y_i) = p_{W_i}(X_i)\,
  (1 - p_{W_i}(X_i))$:
  $$
  \operatorname{Cov}\big(\hat\beta(g), \hat\beta(g')\mid D\big)
  = \sum_{a\in\{0,1\}} \sum_{i \in g\cap g',\, W_i = a}
    \frac{p_a(X_i)\,(1 - p_a(X_i))}{n_{g,a}\, n_{g',a}},
  $$
  exact and enumerable; `sandwich_var_ols` is precisely its
  single-dataset plug-in at the basis level (the Level 3 layer,
  pre-built).
- *The law becomes CLT, well-governed.* Level 2's exact Gaussianity
  becomes a bounded-summand CLT whose quality is set by the smallest
  arm-by-region binomial counts — the April document's
  $d_{\mathrm{eff}} = n_H\,\bar p(1-\bar p)$ in new clothing.
- *Every functional transfers unchanged.* Declaration, split-half $P_1$,
  and selection (including maxeffCons and zcons variants) are the same
  `pmvnorm`/Gaussian-level computations on the new mean vector and
  covariance.

**Role 2 — the influence-function template for the nonlinear scales.**
The package's binary screening events live on the log-OR estimate (ratio
thresholds are logged at `R/forestsearch_main.R:1734`; the comparison at
`R/subgroup_search.R:626` is against that estimate), so the faithful
binary process is the per-region logistic treatment coefficient — a
delta-method-normal vector whose covariance is the overlap sum of
influence contributions. `sandwich_var_glm` *is* the arm-level version of
that algebra; restricting the basis to region indicators gives the
candidate version, with Tier B exact log-OR means on `df_super` and one
new ladder rung (per-region MLE linearization, governed by region event
counts). For survival, the Poisson-plus-offset lane makes the per-region
arm fit closed-form — $\hat\xi_a(g) = \log(D_{g,a}/T_{g,a})$, events over
exposure — so the candidate log-IRR process has means computable exactly
and information carried by events. The natural covariance conjecture,
consistent with the April $d_{\mathrm{eff}} = $ events for survival, is
**event-overlap**:
$$
\operatorname{corr}\big(\hat\beta(g), \hat\beta(g')\big) \;\approx\;
\frac{D_{g\cap g'}}{\sqrt{D_g\, D_{g'}}},
$$
reducing to person-overlap under homogeneous rates. The Cox-faithful lane
replaces Poisson scores by partial-likelihood influence terms (dfbetas) —
objects the MR machinery already computes — giving the realized-family
Tier B covariance for Cox from existing internals.

**Role 3 — the gate itself becomes a predictable operating
characteristic.** Because the OLS Crump vector $\Delta$ is linear in $Y$,
the pair (candidate process, $\Delta$) is jointly Gaussian with
cross-covariance from the same $\sum_i c\,c'\,\mathrm{Var}(Y_i)$ sum. The
validation document's recommended workflow — gate on the constant-CATE
test, then search — therefore has computable joint operating
characteristics: $P(\text{gate rejects})$,
$P(\text{gate} \wedge \text{family declares})$, and the conditional
selection law given the gate are all functionals of one known law. This is
verified below with the *actual* package test in the loop.

## 4. What the Crump machinery does not supply

Per-candidate means and floors (the mixture/Tier B layer carries those);
the selection/argmax distribution (the Rung C layer); Tier A means on
non-collapsible scales (OR, HR — Tier B as established); family
regeneration and nomination (Level 5, unchanged); and, on the survival
Poisson lane, censoring beyond exposure-offset low-rate fidelity — the
Cox-IF lane and the approximation-error measurement own that boundary.

## 5. Verified first rung (B1) and the remaining ladder

`check_binary_rd_process.R` — on-record verification: primary machine
(pop-os), 2026-08-24, R 4.6.1 (2026-06-24), mvtnorm 1.4.2, seed 20260824, **8 PASS / 0 FAIL**,
2.7 s, output digit-identical (runtime aside) to the drafting sandbox
(R 4.3.3, mvtnorm 1.2.4, 4.4 s); family = the committed MD check's 8-profile geometry
(M = 14 after the 0.10 floor); logistic DGM with benefit outside $Q$ and
harm inside (RD means $[-0.086, 0.211]$, $\beta_{\mathrm{RD}}$ at the
$Q$-rule $= 0.211$), floor $c_1 = 0.10$, $n = 600$, 4000 data-level
replicates.

| Check | Result |
|---|---|
| BR1 RD means: conditional identity | within 2.03 SE; imbalance vs population 0.020 (informational) |
| BR2 heteroskedastic overlap covariance $\equiv C\,\mathrm{diag}(v)\,C^\top$ | $1.9\times 10^{-15}$; empirical corr within MC; homoskedastic approx off up to 6.8% |
| BR3 Rung A on RD, alt | computed 0.9545 vs MC 0.9565 |
| BR3 Rung A on RD, null (no effect, floor 0.10) | computed 0.4039 vs MC 0.4215 — bounded-summand CLT holds at family level |
| BR4 Crump $\Delta$: computed mean/cross-corr vs empirical | within 1.44 SE; max cross-corr dev 0.028 |
| BR4 $P(\text{gate})$, Gaussian (true $V$) vs data (package test, plug-in $\widehat V$) | 0.6801 vs 0.6913 |
| BR4 $P(\text{gate} \wedge \text{declare})$ | 0.6739 vs 0.6857 |
| BR4 coupling | $P(\text{declare}\mid\text{gate}) = 0.991$ vs marginal 0.954 |

Remaining rungs, specified: **B2** (binary/log-OR): per-region logistic
process, IF-overlap covariance vs MC, Tier B means, linearization error
vs smallest region event counts. **S-P** (survival, Poisson+offset):
closed-form per-region log-rate contrasts under an exponential/PWE DGM
with censoring as exposure; verify the event-overlap correlation
conjecture; measure the Cox-approximation error against censoring
fraction. **S-C** (survival, Cox-faithful): dfbeta-overlap covariance
from one fitted family (MR internals) vs MC — the realized-family route.
Sequencing follows the committed workstream order: MD (done) $\to$ binary
$\to$ survival.

## 6. Relation to the committed detection-theory artifacts

The consolidated document
(`quarto/simulations/actg175/continuous/analytic_verification_and_prediction_md_harm.qmd`)
carries the ladder this NOTE extends: Levels 0–4 transfer with the
substitutions above (heteroskedastic weights; delta-method rung for
nonlinear scales; event-carried information for survival); Levels 5–6
(regeneration, gate composition, MR fidelity) are outcome-agnostic and
carry over verbatim. The binary/survival boundary rows of that document's
Part V are the ones this program moves.
