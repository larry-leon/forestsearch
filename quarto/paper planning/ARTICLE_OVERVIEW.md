# Article overview — procedure-agnostic post-selection inference for subgroup treatment effects

**Working title.** *Post-selection inference for identified subgroups: a wild-bootstrap
recasting of the León (2024) de-biased estimator and infinitesimal-jackknife interval
for difference-in-natural-parameter estimands.*

**One-paragraph thesis.** León et al. (2024) introduced, for the forest-search procedure
in the Cox model, a bootstrap bias-corrected subgroup effect with an infinitesimal-jackknife
(IJ) confidence interval. We show this post-selection inference is neither tied to forest
search nor to survival outcomes. It is a general method for any subgroup treatment effect
that is a **difference in natural parameters (DINA)** — Cox log-hazard-ratio, logistic
log-odds-ratio, Poisson log-incidence-rate-ratio, Gaussian mean difference — and for any
selection procedure. The estimand layer (DINA, after Gao & Hastie 2025) supplies a single
additive scale on which each subject contributes a `dfbeta` influence; the inference layer
is a shared-multiplier **wild / martingale bootstrap** of those influences, which to first
order reproduces León's nonparametric bootstrap de-biasing and IJ variance exactly. The two
layers compose so that screening, selection, bias correction, and interval estimation all run
from one influence matrix, turning post-selection inference into a fast wrapper around any
subgroup identifier.

---

## Contributions

1. **Decoupling inference from the selection procedure.** The de-biased estimator and IJ
   variance depend on the identification procedure only through a selection map applied to
   candidate estimates. They therefore apply to *any* procedure — forest search, GRF /
   policy trees, virtual twins, prespecified-then-adaptive — via the literal bootstrap
   re-run (Tier-1).

2. **A DINA characterization of the admissible estimand class.** The method covers exactly
   the treatment effects expressible as a difference in natural parameters of a Cox or GLM
   model, unifying survival and GLM outcomes on one canonical scale with a common
   `dfbeta` influence representation.

3. **A wild-bootstrap recasting (the fast tier).** Replacing each bootstrap re-search with a
   shared-multiplier perturbation of the candidate influences yields a refit-free computation
   (one influence-matrix product + multiplier draws). With centered-Poisson multipliers this
   *is* the infinitesimal-jackknife linearization of the nonparametric bootstrap.

4. **An exact-correspondence result.** To first order the fast estimator equals León's Eq. 7
   bias correction — point estimate *and* IJ interval — for both the identified subgroup and
   its complement, with the only approximation the `dfbeta` linearization. The selection map
   enters both tiers identically, so its non-smoothness is not a source of disagreement.

5. **A single-influence-matrix pipeline.** The same N×S influence matrix drives the
   closed-form consistency screen, the family-level (false-discovery) routine, and the
   post-selection de-biasing — so the whole identify → screen → select → correct → interval
   workflow needs one fit per candidate and no refitting.

6. **A clean home for modern HTE identifiers (DINA, GRF).** Aligning their selection to the
   inferential effect coerces them into the effect-rank-over-a-family template the correction
   requires, so the same de-biasing and IJ interval wrap them soundly — with the residual
   data-driven-family caveat made explicit.

---

## The two layers, stated precisely (for the framing section)

**Estimand layer — DINA (Gao & Hastie 2025).** The per-arm natural (canonical) parameter and
the treatment effect as their difference. Gives (i) a common additive scale across subgroups,
(ii) the `dfbeta` as the influence on that scale, (iii) the closed-form consistency rate.
*Boundary:* the inference machinery needs only a regular asymptotically linear (RAL) estimator
with a computable influence function; DINA is the natural, essentially complete instantiation
for Cox/GLM and the one that also yields the closed-form screen. Non-DINA-but-RAL effects
(RMST difference, win ratio) could use the same multiplier-bootstrap inference but forgo the
closed-form screen and the additive scale — note as scope/future.

**Inference layer — wild / multiplier (martingale) bootstrap of the influence functions.**
The Cox score is a martingale integral (Lin–Wei–Ying 1993; weird bootstrap, Dobler & Pauly
2017); the GLM score a sum of mean-zero terms; in both the multiplier bootstrap of the
`dfbeta` is the same object. A single shared multiplier vector per draw, broadcast to all
candidates, reproduces the candidates' overlap covariance and hence the bootstrap's *selection*
distribution. Centered-Poisson weights = the nonparametric bootstrap multiplicity law;
Rademacher/Gaussian share the first two moments (wild-bootstrap family) and give the same
leading-order correction.

**Composition.** A DINA estimand is RAL with `dfbeta` = the score-based influence on the
natural-parameter scale, so the estimand layer hands the inference layer exactly the objects
it needs. One influence matrix; three uses; two tiers.

**Two tiers of generality (state both so the claim is exact).**
- *Tier-1 (literal bootstrap re-run):* defined for any procedure that can be re-run on a
  resample; the universal definition of the de-biased estimator + IJ variance; cost = R
  re-searches.
- *Tier-2 (multiplier / wild bootstrap):* additionally requires an enumerable candidate family
  with influences and a selection map that ranks on the **inferential effect** β̂(g); then it is a
  first-order-exact, cheap surrogate for Tier-1. Procedures that natively rank on another statistic
  (DINA's subgroup-mean τ̂, GRF's mean DR score) are accommodated by re-ranking their qualifying
  candidates on β̂(g) (`select_statistic = "effect"`), which makes them candidate *generators* with
  a forest-search-identical selection layer — exact conditional on the proposed family. A
  data-driven family carries a conditional-vs-unconditional estimand distinction (family
  regeneration), to be quantified against a literal Tier-1.

---

## Proposed section spine

1. **Introduction.** Post-selection (winner's-curse) bias in subgroup identification; the gap
   that inference is usually tied to one procedure and one outcome; the procedure-agnostic,
   DINA-general proposal.
2. **Subgroup treatment effects as DINA estimands.** Cox/GLM as differences in natural
   parameters; the M-estimator / `dfbeta` influence representation; the common scale.
3. **Selection and its bias.** The selection map; the winner's curse; León's de-biased
   estimator and IJ variance restated *generally* (not procedure-specific).
4. **The wild / martingale multiplier bootstrap of influences.** Shared-multiplier
   linearization; why one draw mimics a bootstrap re-search (the resample-as-reweighting →
   dfbeta-as-refit → shared-vector-as-correlated-competition → Poisson-as-bootstrap-weight
   chain); the exact-correspondence proposition (point + IJ, Ĥ and Ĥᶜ).
5. **Two tiers and computation.** Literal re-run vs multiplier surrogate; the one-influence-matrix
   pipeline; cost and where the fast tier is available.
6. **The closed-form consistency screen and family routine** (the other two uses of the
   influence matrix) — supporting, can be condensed or referenced.
7. **Simulation study.** León (2024) Table 2 design (M1–M3, targets θ̂(H)/θ‡/θ†, bias +
   coverage + SD calibration) extended with the Tier-2 de-biased + IJ block; paired Tier-1↔Tier-2
   agreement; the small-nₛ degradation frontier; GLM (OR) closure. *(See the simulation hand-off.)*
8. **Applications.** GBSG/Cox (de-biased HR with IJ CI; harm subgroup + complement) and
   ACTG175/OR (de-biased OR), as concrete procedure-agnostic demonstrations.
9. **Discussion.** RAL beyond DINA; small-nₛ linearization limits; use as an in-fit screen vs
   final-report bootstrap; applicability across identification procedures.

---

## What is already written vs to-draft (maps to `consistency_resampling_theory.qmd`)

| Article piece | Status / source |
|---|---|
| DINA scale, `dfbeta` representation, common-scale schema | Written — `@sec-glm`, `@sec-schema`, `@eq-dina` |
| Closed-form consistency screen | Written — `@sec-criterion`, `@eq-closed`, `@sec-leon-eq3` |
| Family / false-discovery routine | Written — `@sec-if-design`, `@eq-fdr` |
| Multiplier-bootstrap equivalences (Rademacher/Gaussian/Poisson) | Written — `@sec-multiplier`, consequence 2 |
| Selection map + winner's curse | Written — `@sec-debias`, `@eq-selop` |
| Shared-multiplier linearization + "why one draw mimics a re-search" | Written — `@sec-debias` intuition subsection, `@eq-linlemma`, `@eq-perturb`, `@eq-cov`, `@eq-reselect` |
| De-biased estimator (Ĥ and Ĥᶜ) | Written — `@eq-biasterms`, `@eq-debiased`, `@eq-debiased-c` |
| Exact correspondence to León Eq. 7 (point) | Written — `@sec-debias-corr`, `@eq-leon7` |
| IJ variance in multiplier form (interval) | Written — `@sec-debias-ij`, `@eq-ij`, `@eq-ijmult` |
| Procedure-agnostic framing + selection-criterion alignment (Tier-1 universal / Tier-2 fast; DINA & GRF) | Written — `@sec-debias-general`, `@eq-selfunc` |
| Estimand-class boundary (RAL ⊇ DINA); conditional-vs-unconditional family | Written — `@sec-debias-general` (Scope) |
| Simulation study | **To run + write** — hand-off prepared |
| Applications writeups | Partial — GBSG/ACTG175 analyses exist as Quarto docs |

---

## Headline results to feature

- **Correspondence (point):** β̃(Ĥ) = β̂(Ĥ) − selection_bias − fixed_bias equals León Eq. 7
  to first order, via η\*_b(g) = Σ_{i∈g}(K\*_bi − 1)·dfbeta_{g,i} = D_g(b) (Ĥ and Ĥᶜ).
- **Correspondence (interval):** the IJ variance (León Eq. VInfJ_bc) is computable from the
  same multiplier draws; the de-biased IJ interval is wider than the naive robust-SE interval
  because it carries the selection variability.
- **Validation already in hand:** ACTG175/OR (T2 2.91 vs T1 2.60); GBSG/Cox (T2 HR 1.81
  (0.71, 4.67) vs T1 1.93 (0.88, 4.24); IJ SE 0.48 vs robust 0.36; ≈637× faster). Point
  tracks; T2 IJ SE ~20% conservative, largest at small nₛ — to be characterized in simulation.

---

## Positioning relative to prior work

- **León et al. (2024):** origin of the de-biased Cox estimator + IJ; this paper generalizes
  (any procedure, Cox *and* GLM) and provides the fast wild-bootstrap computation with the
  correspondence guarantee.
- **Harrell et al. (1996):** the optimism-correction lineage; note that León retains the
  second bias term η\*_b(Ĥ) that Harrell omits, and that this term carries the IJ residual.
- **Efron (2014); Wager, Hastie & Efron (2014):** the IJ-for-bagging variance the interval
  uses.
- **Lin–Wei–Ying (1993); Parzen–Wei–Ying (1994); Jin–Ying–Wei (2001); Dobler & Pauly (2017):**
  the martingale / wild bootstrap of the estimating function that the inference layer is an
  instance of.
- **Gao & Hastie (2025):** the DINA estimand framework that characterizes the admissible class.
- **GRF / virtual twins (Athey–Imbens; Wager–Athey; Foster et al.):** alternative *selection*
  procedures the same inference can wrap.

---

## Open decisions to settle before drafting

- Lead framing: "post-selection inference, procedure-agnostic" vs "wild-bootstrap recasting of
  León" — both true; pick the title emphasis.
- How much of the consistency-screen / FDR material to carry (it is the same influence matrix,
  but is supporting rather than central to the *inference* thesis) — full section vs appendix
  vs reference to the methods note.
- Whether the GLM extension is a co-equal contribution or a generalization remark.
- Target venue and length (drives §6 depth and whether the schema/`@sec-schema` is included).
- Author list and division of labor (Larry runs all simulations; co-authors do not).
