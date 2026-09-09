# MEMO — Xu & Guo (2026), "In-Sample Evaluation of Subgroups Identified by Generic Machine Learning": reading, mapping onto forestsearch/MR, and go/no-go for the two comparisons

**Date:** 2026-09-08. Author: the Linux MR-field chat. **Source:** arXiv:2605.03141v1 (stat.ME, submitted 4 May 2026; Shuoxun Xu, Xinzhou Guo, HKUST), read from the PDF text and the arXiv TeX source tarball. Every quoted specification below is verbatim from the main text; derived quantities are marked.
**Status:** report and go/no-go, per Larry's instruction — no task document drafted.

## 0. Code and supplement availability (the "3-day job" gate)

- The arXiv source tarball contains only `arxiv-in-sample-subgroup.tex`, `refs.bib`, and four figure PDFs (`PCI_linear_10`, `PCI_linear_1`, `PCI_rpart_10`, `PCI_loess_1`). **No supplement and no R code are in the arXiv upload.**
- The paper states: "Supplementary Material A–E and the R code to reproduce the numerical results referenced in Sections 2.1, 2.2, 4 and 5 are available online" and "Simulation code reproducing the numerical results in Section 4 is provided as Supplementary Material" — i.e. a journal supplement not yet public. No GitHub link anywhere in the text.
- **Consequence:** reproducing *their method* is off the table (as Larry intended). But the two things we actually need — their simulation DGM and their ACTG175 specification — are **fully stated in the main text**. Neither comparison needs their code.

## 1. Their method, from source (what "conditional adaptive perturbation" is)

- **Target.** PISA(c) := E(Y(1) − Y(0) | D̂(Z) ≥ c): the average treatment effect over the post-hoc subgroup of subjects whose *estimated* CATE clears a constant threshold c (main text: "we focus on a constant c"; estimated-c variants are in their Supplement C). The subgroup is a superlevel set of a CATE estimate D̂ fit once on the whole dataset by any method, "parametric or nonparametric, correctly specified or misspecified, and even black-box."
- **Algorithm 1.** Compute the pivotal statistic on a subset I ⊆ [n] of size m: PISA_I(c) = Σ_{i∈I} ψ(O_i, π̂^{(−q(i))}, ĥ^{(−q(i))}) 1{D̂(Z_i) ≥ c} / Σ_{i∈I} 1{D̂(Z_i) ≥ c}, with ψ the doubly-robust (AIPW-type) influence function and cross-fitted nuisances (π by logistic regression, h by B-spline in their runs). For j = 1..M: draw V_i ~ N(1, 1) i.i.d. over I, recompute the V-weighted ratio, and take empirical quantiles of m^{1/2}(PISA*_I − PISA_I); the interval scales those quantiles back by m^{−1/2}. **D̂ is held fixed throughout** — "keep D̂ fixed to address selection bias and stabilize distribution"; there is no re-identification inside the perturbation.
- **Adaptive m.** m ≤ n chosen by a data-adaptive rule (their Algorithm 2, Supplement E, constant C = 1/2): m = n gives full efficiency under regularity; m = o(n) restores validity under nonregularity (point masses / non-smooth boundary of {D ≥ c}). Their M = 2000 perturbations.
- **Theory.** Validity for PISA(c) regardless of regularity; "triple robustness": full √n efficiency "as long as any two of the three components, subgroup identification, outcome regression and propensity score, converge in sufficiently fast nonparametric rates." Their own argument for why sample-splitting fails: it targets PISA_I(c) on the identification half, which "generally differs from PISA(c)."
- **What "conditional" means here:** conditional on the fitted D̂ (a perturbation scheme), *not* conditional on a detection event. There is no "no subgroup found" state in the framework — a threshold set always exists (possibly empty).

## 2. Mapping onto forestsearch / MR / field

| Dimension | Xu & Guo (2026) | forestsearch + MR/field (ours) |
|---|---|---|
| Target | PISA(c): effect in {D̂(Z) ≥ c}, a curve in c; depends on the chosen D̂ | β(Ĥ), β(Ĥᶜ): conditional estimands of one selected covariate-cut subgroup and its complement |
| Identifier | Any CATE learner, superlevel set; not interpretable as cuts | Enumerated conjunctions of covariate cuts + consistency screen; interpretable; **can return nothing** |
| Selection-bias handling | Subsample-scale (m-out-of-n) perturbation with D̂ fixed; stabilizes the nonregular law | Shared multiplier broadcast to the whole family with **re-selection** per draw; two-term de-bias; field inversion |
| Nonregularity | Formally covered (m = o(n)), at a width cost | No formal theorem; certified by campaign record; the field's harm side robust across dial/grid/pick |
| Non-detection | Not a state | Central: guarantees conditional on a declaration (our §5 level-dimension item) |
| Nuisances | DR with cross-fitted π̂, ĥ (observational-ready) | Standard unadjusted trial analysis (covariate adjustment is a closed line) |
| Output | Pointwise CI curve; no named subgroup | Named Ĥ, Ĥᶜ, corrected effects, one-sided bounds, joint pair |
| Compute | One D̂ fit + M weighted ratios (cheap) | Family enumeration + screen + multiplier field (heavier, but no refits per draw) |

The two lines answer different questions on the same data: "how large is the effect among those the model rates above c?" versus "which describable region does the search return, and what is its effect?" A comparison is legitimate only as *what each flags on the same specification*, and as *coverage/width relative to each method's own oracle on the same DGM* — never as one number against another.

## 3. Their simulations, verbatim (Section 4) — the material for Part 2

**DGM.** n = 1000. Treatment G ~ Bernoulli((1 + e^{−0.5(Z1+Z2)})^{−1}) (confounded by design); Y = h(0, Z) + G·D(Z) + e, e ~ N(0, 0.4²); Z = (Z1, Z2). For (C), (D): Z2 = T·U + (1−T)·X with T ~ Bernoulli(0.5), U ~ U[−1, 1], P(X = −1) = P(X = −0.8) = P(X = 0.8) = P(X = 1) = 1/16, P(X = 0) = 3/4 (a point mass; nonregular). π estimated by logistic regression, h by B-spline. M = 2000 perturbations; **500 repetitions**.

**Table 1 (theirs).**

| Setting | h(0, Z) | D(Z) | F_Z | D̂ |
|---|---|---|---|---|
| (A) | 1{Z2 ≥ 0} − 0.06 | 1{Z2 ≥ 0} − 0.06 | (U[−1,1], U[−1,1]) | parametric (linear) |
| (B) | 1{Z2 ≥ 0} − 0.06 | 1{Z2 ≥ 0} − 0.06 | (U[−1,1], U[−1,1]) | nonparametric (B-spline) |
| (C) | 1{Z2 ≥ 0}·\|Z2\|^{1/3} | [Z2 − 0.95·sign(Z2)]·1{\|Z2\| ≤ 0.95} | (U[−1,1], Mixture) | parametric |
| (D) | 1{Z2 ≥ 0}·\|Z2\|^{1/3} | [Z2 − 0.95·sign(Z2)]·1{\|Z2\| ≤ 0.95} | (U[−1,1], Mixture) | nonparametric |

**Table 2 (theirs): ECP / CIL, two-sided 95% for PISA(c).**

| Setting | Naive | SS | Oracle | m = n | m = n/2 | m = n/4 | m = n/8 | Alg. 2 (adaptive) |
|---|---|---|---|---|---|---|---|---|
| (A) ECP / CIL | 91.2 / 0.15 | 74.6 / 0.22 | 94.6 / 0.17 | 93.4 / 0.15 | 94.8 / 0.22 | 95.4 / 0.31 | 94.4 / 0.43 | 93.2 / 0.16 |
| (B) | 91.0 / 0.17 | 69.6 / 0.23 | 93.4 / 0.14 | 93.4 / 0.17 | 95.2 / 0.23 | 94.0 / 0.33 | 94.2 / 0.47 | 93.6 / 0.17 |
| (C) | 88.5 / 0.31 | 78.7 / 0.44 | 94.8 / 0.29 | 86.9 / 0.17 | 89.4 / 0.24 | 92.2 / 0.34 | 94.2 / 0.48 | 94.2 / 0.44 |
| (D) | 78.2 / 0.28 | 57.1 / 0.35 | 96.6 / 0.26 | 80.8 / 0.16 | 89.8 / 0.22 | 93.4 / 0.32 | 93.2 / 0.47 | 93.8 / 0.42 |

(The threshold c used in the simulations is not printed in the main text.)

**Where they perform very well:** the regular settings (A), (B) — adaptive m ≈ n, coverage 93.2–93.6, CIL 0.16–0.17 at or below the oracle's 0.14–0.17 (ratio ≈ 0.94–1.21, derived). **Where they are more challenged:** the nonregular settings (C), (D) — m = n collapses to 86.9 / 80.8; the adaptive rule restores 94.2 / 93.8 **at CIL 0.44 / 0.42 against oracle 0.29 / 0.26, i.e. ≈ 1.5–1.6× the oracle width** (derived). That width-for-validity trade under nonregularity is the comparison point for us.

**Structure of their truth, derived (for our identification metrics).** (A)/(B): benefit region B = {Z2 ≥ 0}, prevalence 0.50, effect +0.94; complement effect −0.06; boundary at the median of U[−1, 1] — exactly representable by a quantile cut. (C)/(D): D is a sawtooth — for Z2 ∈ [−0.95, 0): D = Z2 + 0.95 ∈ (0, 0.95] (benefit, largest just below 0); for Z2 ∈ (0, 0.95]: D = Z2 − 0.95 ∈ [−0.95, 0) (harm); D(0) = 0 (sign(0) = 0); ±0.8 → ∓0.15; |Z2| ≥ 0.95 → 0. Under the mixture: benefit region prevalence ≈ 0.269 with mean effect ≈ +0.44, harm region ≈ 0.269 with ≈ −0.44, null ≈ 0.46 of which **0.375 sits at the single point Z2 = 0 on the benefit/harm boundary**. That point mass is the nonregularity — and it is also what a cut-based family must sort to one side of a cut (e.g. "Z2 ≤ q" for the largest negative grid value captures most of B; "Z2 ≤ 0" would swallow the null mass). Signal strength under randomized assignment (derived): (A)/(B) effect/SE ≈ 26; (C)/(D) ≈ 9 — identification is not power-limited in either; the challenge is structural.

## 4. Their ACTG175 application, verbatim (Section 5) — the material for Part 1

- **Sample and contrast:** "all the 1046 patients who receive the treatment, ZDV+ddI, or the control, ZDV+zal"; outcome "CD4 count at 20 ± 5 weeks after receiving the therapy." Randomized trial treated as a special case of their observational framework (DR "for efficiency improvement").
- **Covariates (10):** age (years); weight (kg); homosexual activity (yes/no); history of IV drug use (yes/no); Karnofsky score (0–100); history of zidovudine use (ZDV in the 30 days prior to treatment initiation, yes/no); race (white vs non-white); gender (female vs male); antiretroviral history (classical vs experienced); antiretroviral history stratification (1 = classical, 2 = > 1 but ≤ 52 weeks prior ART, 3 = > 52 weeks). In `speff2trial` terms: `age, wtkg, homo, drugs, karnof, z30, race, gender, str2, strat`; outcome `cd420`; arms 1 (ZDV+ddI) vs 2 (ZDV+zal).
- **Working models (from the figure files and text):** 10-covariate linear regression; age-only linear regression; 10-covariate recursive partitioning (rpart); age-only local polynomial (loess). D̂ fit once on all 1046.
- **Findings (qualitative; figures only, no table, no named subgroup):** under the linear D̂ the proposed interval is wider than the naive in both the 10-covariate and age-only cases; the gap and the instability of the sample-split estimate "suggest nonregularity in this study"; near c ≈ 50 the proposed interval covers 0 for the age-only subgroup where the naive and sample-split do not; the linear-model PISA curve is non-monotone, read as misspecification; under the nonparametric D̂ the curve is more monotone, and "the proposed interval around c = 50 in the left panel [10-covariate rpart] suggests that ZDV+ddI is particularly preferred in the 10-covariate post-hoc subgroup at c = 50." Their "subgroup A" is therefore {rpart-D̂(Z) ≥ 50 CD4 cells}, whose membership rule the paper does not print.

## 5. Go / no-go

**Part 1 — application at their specification: GO.** Everything needed is pinned: arms, outcome, the ten covariates with codings, and their reference threshold (50 cells). Deltas from our S7 ACTG175 analysis, all knob-level: control arm ZDV+zal (arms 2) instead of ddI; outcome `cd420`, MD estimand (`effect_measure = "MD"`, confirm-null 0) instead of the no-improvement indicator/OR; covariate set = their ten (drops cd40, cd80, preanti, hemo, symptom relative to ours; adds strat) — note this removes baseline CD4, the anchor of our OR-path subgroup, which is the point of "same specification"; **benefit orientation** ("a larger treatment effect is preferred": Ĝ/Ĝᶜ), with the harm-direction run reported alongside since our identifier can return nothing in either direction; bounds read by location against 50 cells, mirroring their c = 50, never as significance at 0. Output: what we identify (a named cut-defined Ĝ with |Ĝ|, naive MD, MR-corrected MD, field/field-s one-sided lower bound on Ĝ and upper on Ĝᶜ, joint pair, CV stability, NPV with sens/spec/PPV where a reference exists), placed beside their statement of what their curve shows at c = 50. Cost: one rendered application document; ≈ half a CC session. Public data only (cough excluded by standing rule).

**Part 2 — cross-simulation, indirect: GO, with two labeled adaptations.** The DGM is regenerable from Table 1 and the paragraph above without their code. Adaptations, stated in the record rather than hidden: (i) **treatment randomized**, P(G = 1) = 1/2, replacing their logistic propensity — our method is built for standard randomized-trial analyses and covariate adjustment is a closed line; their DR estimator handles confounding via π̂, so removing it puts both methods on the selection problem alone (an *indirect* comparison, as intended); (ii) our target is β(Ĝ)/β(Ĝᶜ) for the identified cut subgroup, not PISA(c); the comparable statistics are **coverage of our own conditional target and interval width relative to our own oracle (refit on the true benefit region)**, set beside their Table 2 ECP and CIL/oracle ratios (≈ 1.0 regular, ≈ 1.5–1.6 nonregular). Identification metrics (detection, |Ĝ|/|B|, sens/spec/PPV/NPV, best-representable sub-region on the grid, as nb20 did for er ≤ 59) against the derived truth regions of §3. Four settings × 500 replicates at n = 1000, matching their replicate count; a two-covariate family is tiny, so per-replicate cost should be far below the survival campaigns — to be projected at Stage 1, expected well under an hour at 100 workers. The predictions to pre-register: (A)/(B) near-deterministic identification and near-oracle width (both methods "very well"); (C)/(D) the point mass at Z2 = 0 and the sawtooth boundary are where their width pays ≈ 55% over oracle — the question is whether a cut-based identifier with the field's shared-multiplier correction lands closer to oracle width at ≥ 0.94 coverage, or whether the boundary mass costs us in identification (Ĝ absorbing the null mass) instead.

## 6. Decision items (for Larry; defaults in brackets)

- **X-1** Part 1 scope: [the benefit-oriented Ĝ analysis at their spec, with the harm-direction run reported alongside; one document in `fs-glms-interpretable`'s applications].
- **X-2** Part 2 adaptations: [randomized assignment; β(Ĝ) target; their four settings at n = 1000, 500 replicates; both benefit- and harm-oriented identification in (C)/(D) since the DGM plants both].
- **X-3** Grid for Z2 in (C)/(D): [the standard J = 10 quantile grid with de-duplication of the point-mass quantiles (the Karnofsky precedent); report the best-representable sub-region of B as nb20 did]. Alternative: an evenly spaced grid — a deliberate stress-design choice, labeled if taken.
- **X-4** Constructions reported: [naive, MR (IJ two-term), field, field-s once E1 clears it; FB none; winner-only/winner-floor excluded].
- **X-5** Stream and machine: the continuous/MD field path lives on `feature/glm-extension-mac`; [run both parts on the Mac stream, or merge the Mac branch into `feature/glm-extension` first (handoff §5 says it adds files only)] — Larry's call; either way after the E1 review closes, one CC task per session.
- **X-6** Sequencing: [Part 1 first (application, half a session), Part 2 second (its own session with Gate 1 compute go)].

## 7. What this memo does not do

No task document and no kickoff (per instruction); no reproduction of their D̂, their perturbation, or PISA(c); no bearing on the submitted paper — the positioning sentence for S3.5 belongs to the revision stream, which may want: "Xu and Guo (2026) extend the de-biased-inference line to subgroups identified by generic machine learning as superlevel sets of a fitted CATE, holding the learner fixed and perturbing at an adaptive subsample scale; the target and the absence of a no-subgroup state differ from those here."
