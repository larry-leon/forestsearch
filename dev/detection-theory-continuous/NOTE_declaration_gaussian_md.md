# Declaration and selection rates as Gaussian functionals: the continuous/MD case

Companion script: `check_declaration_gaussian_md.R` (this directory). Base R +
`mvtnorm`, seeded, standalone; runtime ~3 s.

## 1. Question and scope

Whether the detection-rate side of the operating characteristics — screening
declaration, split-half consistency, and subgroup selection — admits closed-form
(or deterministic-numerics) computation in the continuous/MD family, over a
**fixed enumerated candidate family**. The continuous case is the feasibility
filter: the mean structure is exact (collapsibility), the variance is known (no
link curvature), and the estimator is linear in $Y$ (joint Gaussianity is exact
conditional on design, not a delta-method limit). If closed forms fail here they
fail in the nonlinear families; if they work, binary and survival become the
same architecture with approximate ingredients.

Out of scope, per the session-starter closures: the MD/MR harness decisions,
regenerated families (DINA/GRF), benefit-side $\hat G$ notation, and any change
to committed simulation machinery. This NOTE and its script are additive.

## 2. Conventions verified from source

Branch `feature/glm-extension`, HEAD `3c40467`.

1. **Oriented scale.** For `outcome_type = "continuous"` with
   `adverse_outcome = FALSE` (the default), the estimator closure negates $Y$ so
   that MD $> 0$ means harm (`R/glm_effect_estimators.R:814-821`). All
   quantities below live on this oriented scale; the reference calibration
   `cal_target_md = -40` (raw) is $+40$ oriented.
2. **Estimator.** The within-subgroup OLS treatment coefficient of
   `lm(y ~ treat)` — identically the difference in arm means (same file,
   `.make_lm_estimator`).
3. **Events are defined on point estimates.** The screening floor compares the
   estimate to the threshold (`R/subgroup_search.R:626`), and the consistency
   re-filter does the same (`R/subgroup_consistency_main.R:528-531`); no test
   statistic is formed (`dev/glm-continuous-sims/NOTE_threshold_sign_md.md`).
   Estimated variance therefore enters the theory only through plug-in
   $\widehat\Sigma$, never through the event definitions.
4. **Mean structure.** The DGM's subject-level effect is two-valued:
   $\tau(X) = \delta + \beta_{\mathrm{inter}}\,\mathbb 1\{X \in Q\}$ with
   $\delta \ne 0$ (`dev/glm-continuous-sims/NOTE_target_is_collapsibility.md`
   §3; vignette calibration $\delta = -26.255236$, target MD$(Q) = 40$
   oriented). Under randomization, for **any** region $g$,
   $\beta(g) = \delta + \beta_{\mathrm{inter}}\,P(Q \mid g)$ — collapsibility,
   requiring no covariate independence and no rectangle shape (same NOTE, §2).
5. **Thresholds.** The committed discriminating cell $c_1 = 30$ (screening),
   $c_2 = 10$ (consistency), oriented scale
   (`dev/glm-continuous-sims/DGM_actg175_grounded_cells.md`).

## 3. The candidate-effect process

Fix a candidate family $\mathcal G = \{g_1, \dots, g_M\}$ (conjunction rules
above a prevalence floor). For a design $D = \{(X_i, A_i)\}_{i=1}^n$ write
$n_{g,a}$ for arm counts within $g$ and let
$$\widehat\beta(g) \;=\; \bar Y_1(g) - \bar Y_0(g)
  \;=\; \textstyle\sum_i c_i(g)\, Y_i, \qquad
  c_i(g) = \mathbb 1\{i \in g\}\Big(\tfrac{A_i}{n_{g,1}} -
  \tfrac{1-A_i}{n_{g,0}}\Big).$$

**Exact conditional law.** Conditional on $D$, with Gaussian noise
$\varepsilon_i \sim N(0, \sigma^2)$, the vector
$\widehat\beta = (\widehat\beta(g_1), \dots, \widehat\beta(g_M))$ is *exactly*
multivariate normal:
$$\widehat\beta \mid D \;\sim\; N_M\big(m(D),\, \sigma^2 C C^\top\big),
  \qquad m(D) = C\,E[Y \mid D],$$
with $C$ the $M \times n$ coefficient matrix. The conditional mean decomposes
as $m(D) = a_0 + \delta\, a_1 + \beta_{\mathrm{inter}}\, a_2$ (prognostic
imbalance, treatment, and realized-$P(Q\mid g, \text{arm } 1)$ components),
linear in $(\delta, \beta_{\mathrm{inter}})$ — so scenario families (null,
borderline, alternative) reuse one geometry.

**Covariance = overlap.** Expanding $CC^\top$,
$$\operatorname{Cov}\big(\widehat\beta(g), \widehat\beta(g')\mid D\big)
  = \sigma^2\Big(\tfrac{n_{\cap,1}}{n_{g,1}n_{g',1}} +
                 \tfrac{n_{\cap,0}}{n_{g,0}n_{g',0}}\Big),$$
where $n_{\cap,a}$ counts the intersection by arm. Under balanced arms this is
$4\sigma^2 n_{g\cap g'}/(n_g n_{g'})$ and the correlation is the overlap
fraction $P(g \cap g')/\sqrt{P(g)P(g')}$, free of $\sigma$.

**Population limit.** Averaging over designs, $m(D) \to \beta(\cdot)$ (the
collapsibility means) and the variances acquire enumerable inflation terms; see
§6, Check 6d, for exactly how far the naive population plug-in carries.

## 4. Declaration and selection as Gaussian functionals

With $Z = \widehat\beta$ on the oriented scale:

- **Rung A (screening exceedance).**
  $P(\max_m Z_m \ge c_1) = 1 - \Phi_M(c_1 \mathbf 1;\, m, \Sigma)$ — a
  rectangle probability, evaluated by `mvtnorm::pmvnorm` (Genz–Bretz
  quasi-Monte Carlo: deterministic-in-principle numerics on the limiting law,
  not data simulation). Tail rates require tightened tolerances
  (`abseps` $\le 10^{-5}$); the default $10^{-3}$ is 40%+ relative error at
  $P \approx 0.002$.
- **Rung B (idealized split-half consistency).** Two disjoint halves give a
  $2M$-block Gaussian with zero cross-half covariance. Per candidate,
  $\text{pass}_m = \{W_1 + W_2 \ge 2c_1,\ \min(W_1, W_2) \ge c_2\}$; the
  family event $\{\exists m: \text{pass}_m\}$ is a union of half-space
  intersections — not a rectangle — evaluated by Gaussian-level Monte Carlo on
  the exact $2M$ law (no data, no refits). The per-candidate probability is a
  one-dimensional integral,
  $P_1 = \int_{c_2}^\infty \phi_{m_1,s_1}(w)\,
  \bar\Phi_{m_2,s_2}\big(\max(c_2,\, 2c_1 - w)\big)\, dw$ — the
  Leon et al. 2024 (SIM, §2.1) object, per candidate with its own
  $(m_j, s_j)$.
- **Rung C (selection identity).** For the argmax rule,
  $P(\widehat H = g_m,\ \text{declared}) =
  P\big(Z_m \ge c_1,\ Z_m - Z_j \ge 0\ \forall j \ne m\big)$ — a rectangle in
  the difference transform, one `pmvnorm` call per candidate. The identity
  $\sum_m P(\widehat H = g_m, \text{declared}) = P(\text{Rung A})$ is an
  internal consistency constraint (verified to $3\times 10^{-7}$).

**Idealizations relative to forestsearch, stated.** (i) Rung B is a single
two-way split with the $2k_1/k_2$ screen; the real procedure's
$\hat p_{\mathrm{cons}}$ aggregates many random splits — its composition uses
the stability closed form already in `_02_identification` and is a deferred
rung, not a gap in the present objects. (ii) Rung C's argmax-over-all is the
`sg_focus`-restricted maximum over *admitted* candidates in the real procedure;
restricting the argmax to the admitted set is the same transform on a subset.
(iii) Candidate nomination (GRF/LASSO screening into the family) is unmodeled
here exactly as in the April L_eff document; the fixed-family computation is
conditional on the family.

## 5. Reconciliation with prior art

- **`design-checks/oracle_design.R`.** Its oracle draws the same
  $M$-dimensional Gaussian competition and evaluates the *selection-bias*
  functional (mean noise of the argmax — the winner's-curse magnitude). This
  NOTE adds the *declaration and selection-probability* functionals of the
  identical limiting object. One law, two families of functionals; nothing is
  duplicated.
- **`quarto/guides/theory_glm_detection_probability.qmd`.** Its $P_1$ with
  $d_{\mathrm{eff}} = n_H/\sigma^2$ and identity-scale thresholds is the
  $M = 1$ (per-candidate, competition-free) case of Rung B; the split-half
  variance $8/d_{\mathrm{eff}} = 8\sigma^2/n_H$ agrees with the two-arm
  half-sample variance above.
- **`quarto/guides/theory_glm_leff_correction.qmd`.** There
  $L_{\mathrm{eff}}$ is *fitted* from a reference simulation to repair the
  independence composite $1 - (1 - P_1)^{L}$. Here the family rate is computed
  from the joint law, and an implied count becomes an *output*. Two structural
  caveats emerge (Check 7): with prevalence-heterogeneous candidate variances a
  single-$P_1$ normalization is ill-defined in the deep tail (the April
  setting's near-homogeneous $d_{\mathrm{eff}}$ masked this), and the April
  growth of $L_{\mathrm{eff}}$ in $N$ is a family-*regeneration* phenomenon
  that the fixed-family computation deliberately isolates away.

## 6. Checks and initial verification

On-record verification: primary machine (pop-os), 2026-08-23, R 4.6.1 (2026-06-24), mvtnorm 1.4.2,
seed 20260823, total runtime 2.4 s; **16 PASS / 0 FAIL** -- numerically
identical line for line (runtime aside) to the drafting-stage sandbox run
(R 4.3.3, mvtnorm 1.2.4, 3.2 s). Synthetic family:
three binary covariates ($x_1, x_2$ dependent; $P(x_3{=}1) = 0.95$ to
manufacture a near-duplicate), $M = 14$ after a 0.10 prevalence floor;
$Q = \{x_1{=}1 \wedge x_2{=}1\}$, $\delta = -26.26$,
$\beta_{\mathrm{inter}} = 66.26$ ($\tau_Q = +40$, the vignette-calibration
shape), $\sigma = 70$ as an ACTG175-magnitude stand-in (the harness residual
sd is data-fitted, not a quotable constant), $c_1/c_2 = 30/10$.

| # | Derivation checked | Result (sandbox) |
|---|---|---|
| 1 | Family-wide collapsibility, exact enumeration | max dev $1.1\times10^{-14}$ |
| 2 | $CC^\top \equiv$ two-arm overlap formula; empirical mean/corr | $1.9\times10^{-15}$; within MC bands (balanced approx off $\le 2\%$, informational) |
| 3 | Rung A `pmvnorm` vs data MC, 3 scenarios | alt $0.8297$ vs $0.8243$; borderline $0.5281$ vs $0.5400$; null $0.0262$ vs $0.0288$ — all PASS |
| 4 | Rung B family (Gaussian-level vs data) and $P_1$ integral | $0.8206$ vs $0.8187$; per-candidate max dev $7.5\times10^{-4}$ |
| 5 | Rung C selection `pmvnorm` vs MC; sum identity | max dev $7.8\times10^{-3}$; identity to $3\times10^{-7}$ |
| 6a | Near-collinearity (realized corr $0.9803$) | `pmvnorm` error $3\times10^{-8}$; exact duplicates require de-duplication |
| 6b | $n = 150$, min arm 10, lognormal noise (skew $\approx 6.2$) | Gaussian computation within $0.013$ (Gaussian noise) and $0.016$ (skewed) of MC |
| 6c | Plug-in $\widehat\Sigma$ from one dataset | rate $0.8415$ vs $0.8297$ truth; lm SEs inflate $0.94$–$1.09\times$ (absorb prognostic + within-$g$ CATE spread) |
| 6d | Population level, random designs | alt: corrected within $0.007$ of MC (PASS). Borderline diagnostic: MC $0.5545$ vs computed $0.5360$ |
| 7 | Implied $L_{\mathrm{eff}}$, null, $N \in \{300, 600, 1200\}$ | $L_{\mathrm{eff}}(\max)$: $2.33, 1.80, 1.08$; family/union: $0.671, 0.913, 0.989$ |

**Findings of substance.**

1. **Competition is first-order.** At the alternative, the per-candidate
   $P_1$-style rate for $Q$ is $0.8078$, but
   $P(\text{select } Q) = 0.3743$ with $P(\text{select dup}) = 0.4553$: the
   near-duplicate cannibalizes selection while the union
   ($0.8296$) matches the family rate. Per-candidate detection probabilities
   are not selection rates; only the joint law separates the two. This is
   consistent with the tie-dynamics observations in the committed MD grid
   reports.
2. **The design-conditional law is exact; the population plug-in has a known
   second-order limit.** Naive population variance
   ($4\sigma^2/(nP(g))$) needs the enumerable correction
   $\sigma^2 \to \sigma^2 + \operatorname{Var}(\mu_0\mid g) +
   \operatorname{Var}(\tau\mid g)/2 + \operatorname{Cov}(\mu_0, \tau\mid g)$,
   and at the adversarial borderline point a residual ($-0.0185$) remains that
   is a *correlation* effect: imbalance components of $Q$ and its
   near-duplicate are far less correlated than their $\sigma$ components, so
   the $\sigma$-overlap correlation overstates joint dependence precisely where
   both sit at the threshold. The full additive covariance decomposition is
   enumerable and is the identified deferred item.
3. **CLT robustness is adequate at realistic sizes.** With minimum arm counts
   of 10 and noise skewness $\approx 6.2$, the Gaussian computation tracks MC
   within $0.016$ — the linear estimator's averaging does the work.

## 7. Deferred items

Additive population covariance decomposition (item 2 above); composition of
the multi-split $\hat p_{\mathrm{cons}}$ stability screen; admitted-set
(`sg_focus`) restriction of the argmax; nomination/regenerated families;
adjusted and IPTW estimator variants (covariance changes, stays computable);
transfer conditions to binary (delta-method means, $d_{\mathrm{eff}}$
variances) and survival.

## 8. What this buys the manuscripts

For `fs-glms-interpretable` `_04_operating_characteristics`, computed rates
from the joint law can *corroborate* simulated declaration and selection rates
on fixed families — the language remains "consistent with / corroborated by";
the computation does not replace simulation calibration. The JRSS-B §7 posture
(rates are properties of the identification procedure, not the correction) is
unchanged and, if anything, sharpened: the procedure-level rates are now
exhibited as explicit functionals of one limiting Gaussian object, the same
object whose selection-bias functional `oracle_design.R` already measures.
