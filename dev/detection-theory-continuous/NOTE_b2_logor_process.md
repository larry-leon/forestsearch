# B2: the log-OR candidate-effect process as a transform of the exact proportion core

Companion script: `check_binary_logor_process.R` (this directory; base R +
mvtnorm, standalone, seed 20260825, ~6 s, 14 PASS / 0 FAIL in the drafting
sandbox). Fulfills the B2 rung specified in
`NOTE_crump_bridge_binary_survival.md` §5.

## 1. Claim and verdict

The binary pathway's screening scale is the log odds ratio: the `"OR"`
effect measure fits logistic regression and **returns the log-OR**
(`R/glm_effect_estimators.R:272`), ratio-scale floors are mapped to the log
scale (`R/forestsearch_main.R:1803-1804`), the screen compares that
estimate to the floor (`R/subgroup_search.R:626`; disabled only under
`sg_focus = "maxeff"`), and the binary default `adverse_outcome = TRUE`
makes positive log-OR the harm direction with no negation
(`R/glm_effect_estimators.R:98`). For the unadjusted within-region fit the
treatment coefficient equals $\mathrm{logit}(\hat p_1) -
\mathrm{logit}(\hat p_0)$ exactly (saturated two-group model; verified at
$2.4\times 10^{-12}$ against `glm`).

**Verdict: the log-OR candidate process is the logit-difference transform
of an exactly tractable linear core — B2 = B1 $\circ$ logit — and the only
new approximation rung is the transform's curvature, measured at
$\le 0.005$ in family declaration rates down to region-arm effective
counts of 3.4.** Every functional of the MD/RD theory (declaration,
split-half, selection) transfers by composing with the transform.

## 2. The construction

**The core (exact, B1-class).** For a family $\mathcal G$ of $M$ regions,
stack the $2M$ region-arm event proportions
$\big(\hat p_0(g), \hat p_1(g)\big)_{g}$. Each is linear in $Y$ with
coefficients $\mathbb 1\{i \in g, W_i = a\}/n_{g,a}$, so conditional on
the design the core has exact means $\bar p_{g,a} =
n_{g,a}^{-1}\sum_{i\in g, W=a} p_a(X_i)$ and exact heteroskedastic
overlap covariance $\Sigma_{\mathrm{prop}} = C\,\mathrm{diag}\big(p(1-p)\big)\,C^\top$,
with the cross-arm block identically zero (disjoint subjects; verified at
machine zero). This is the same object verified for the RD pathway (B1);
the two binary scales share one core.

**The transform.** $\hat\theta(g) = \mathrm{logit}\,\hat p_1(g) -
\mathrm{logit}\,\hat p_0(g)$, i.e. $\hat\theta = h(\hat p)$ with $h$
acting coordinatewise. The delta law is
$$
\hat\theta \;\approx\; N_M\!\big(h(\bar p),\; J\,\Sigma_{\mathrm{prop}}\,J^\top\big),
\qquad
J_{g,(g,a)} = \frac{(-1)^{1-a}}{\bar p_{g,a}\,(1 - \bar p_{g,a})},
$$
a diagonal-per-(region, arm) Jacobian. The **governing quantity** for the
new rung is the region-arm effective count
$n_{g,a}\,\bar p_{g,a}(1-\bar p_{g,a})$ — the halves of the April
theory's binary $d_{\mathrm{eff}}$ — which controls both the curvature
error and the boundary event $\hat p \in \{0, 1\}$ (log-OR non-finite;
frequency $2\times 10^{-5}$ at $n = 600$ and $3\times 10^{-4}$ at
$n = 220$ in the verification DGM).

**Tier B means and the non-collapsibility gap, quantified.** The exact
reference is the marginal region log-OR
$\theta(g) = \mathrm{logit}\,p_1(g) - \mathrm{logit}\,p_0(g)$ with
$p_a(g)$ the enumerated region-arm risks — Tier B, since the mixture
identity fails on this scale. The would-be Tier A answer
$\mathbb E[\theta_{\mathrm{CATE}}(X)\mid g]$ differs by up to $0.021$
across the verification family (and by $0.0009$ at the homogeneous
$Q$-rule): the Tier A $\to$ Tier B degradation is now a measured number,
not only a structural statement.

**Three-way error decomposition (the check architecture).** Any computed
rate can be evaluated three ways: (i) the delta-Gaussian `pmvnorm`
functional; (ii) *transform-exact* Gaussian-level Monte Carlo — sample the
exact core law, apply $h$, then evaluate the event; (iii) data-level Monte
Carlo. Then (iii)$-$(ii) isolates the core's binomial-vs-Gaussian CLT and
(ii)$-$(i) isolates linearization. Both components are measured
separately below; conflating them is how a good approximation gets
misjudged.

## 3. Verified results (drafting sandbox, seed 20260825, 14 PASS / 0 FAIL)

DGM and family identical to B1 (8 enumerable profiles, $M = 14$, logistic
arm risks with benefit outside $Q$ and harm inside; region log-OR range
$[-0.350, 0.900]$, $\theta(Q\text{-rule}) = 0.899$); floors
$c_1 = \log 1.5$, $c_2 = \log 1.15$; $n = 600$, 4000 data replicates,
$2\times 10^5$ Gaussian-level draws.

| Check | Result |
|---|---|
| BL1 non-collapsibility gap | max 0.021 over family; 0.0009 at the $Q$-rule |
| BL2 core exactness | cross-arm block $= 0$ exactly; means within 1.28 SE; corr within MC |
| BL3 delta law of $\hat\theta$ | means within tolerance (delta bias $O(1/n)$); max corr dev 0.040; min effective count 7.2 |
| BL4 Rung A, alt | delta 0.9530 \| transform-exact 0.9578 \| data 0.9570 — linearization 0.005, core CLT 0.0008 |
| BL4 Rung A, null | 0.4022 \| 0.4071 \| 0.4170 — linearization 0.005, core CLT 0.010 |
| BL5 selection | max dev $6.3\times 10^{-3}$; sum identity $2.6\times 10^{-5}$; $Q$-rule 0.495 vs near-dup 0.430 |
| BL6 $n = 220$ stress | min arm 15, min effective count 3.4; linearization 0.005, core CLT 0.001; non-finite freq $3\times 10^{-4}$ |
| BL7 estimator tie | `glm` treatment coefficient $=$ logit difference at $2.4\times 10^{-12}$ |

The near-duplicate competition structure of the MD theory reappears
unchanged on the log-OR scale (BL5), as it must: selection is a property
of the joint law, and the transform preserves the law's geometry to the
measured linearization tolerance.

## 4. What this settles and what remains

Settled: the binary pathway's *executing* scale has a verified
semi-analytic prediction layer — declaration, consistency, and selection
functionals computable from Tier B region-arm risks and the exact core
covariance, with the single new rung (curvature) measured and governed by
a named quantity. The RD lane (B1) and the log-OR lane (B2) are two
transforms of one core, so a design can be predicted on both scales from
one set of enumerated inputs.

Remaining, unchanged in scope: the population-level layer for log-OR
composes the additive covariance decomposition of the continuous document
(§Level 4) on the *core* and then transforms — derivable, deferred; the
multi-split $\hat p_{\mathrm{cons}}$ composition; family regeneration;
and the survival lanes **S-P** (Poisson-plus-offset with the
event-overlap conjecture) and **S-C** (Cox dfbeta overlap), next per the
Crump-bridge NOTE's ladder.
