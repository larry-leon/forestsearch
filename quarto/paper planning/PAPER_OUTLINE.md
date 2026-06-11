# Paper outline (rigor-first) — procedure-agnostic post-selection inference

**Working title.** *Procedure-agnostic post-selection inference for identified subgroups:
a multiplier-bootstrap recasting of the León (2024) de-biased estimator and
infinitesimal-jackknife interval.*

**Organizing principle.** Lead with the claim that is *rigorously* established and state its
scope honestly. The feasible, rigorous target is the **conditional-on-family** correspondence
for procedures that select on the inferential Cox/GLM effect. Everything broader (data-driven
family regeneration; faithful linearization of native non-effect selection) is positioned as
characterized scope and simulation, not as a proven claim.

---

## The questions, ranked by how rigorously the paper can address them

| # | Question | Rigor level | Where |
|---|----------|-------------|-------|
| **Q1** | For a procedure that ranks on the inferential effect β̂(g) over a **fixed** candidate family, do the fast (Tier-2) de-biased estimator and IJ interval equal León's Eq. 7 (Tier-1)? | **Proven** — first-order identity under the dfbeta linearization (point + interval, Ĥ and Ĥᶜ). | Methods (correspondence) |
| **Q2** | Can DINA / GRF / Lasso-based identifiers be brought under Q1? | **Proven, conditional** — yes *iff* their selection is re-ranked on β̂(g) (`select_statistic = "effect"`), **conditional on the proposed family**. This is the paper's central feasible contribution. | Methods (alignment) |
| **Q3** | Does the fast gate remain valid when the candidate family is itself data-driven (DINA surface / GRF forest regenerated per resample)? | **Characterized, not closed** — the fast gate targets the *conditional* estimand and omits the family-generation component; sign and size quantified empirically against a literal Tier-1 that regenerates the family. | Scope + simulation |
| **Q4** | What are the finite-sample operating characteristics (bias, coverage, SD calibration) across models and separation/size regimes? | **Empirical** — León (2024) Table 2 design extended with the Tier-2 + IJ block and a paired Tier-1↔Tier-2 panel. | Simulation |
| **Q5** | Faithful fast de-biasing of a *native* (non-effect) selection rule. | **Out of scope (named as future work)** — requires the cross-influence between the native statistic and the effect under shared multipliers; a strictly heavier two-influence-matrix construction. | Discussion |

The paper's spine is Q1 → Q2 (the contribution), with Q3–Q5 as honestly-bounded scope.

---

## Introduction (the requested framing)

1. **Post-selection bias in subgroup identification.** A subgroup chosen *because* its effect
   is extreme has an optimistically biased effect (winner's curse). Reporting that effect with
   a naive model CI under-covers.

2. **León et al. (2024): a two-term bootstrap correction with an IJ interval.** On the
   coefficient scale, León's bias-corrected estimator subtracts the average of *two* bootstrap
   discrepancies (their Eq. 7),
   β̃(Ĥ) = β̂(Ĥ) − (1/B) Σ_b [ η\*_b(Ĥ\*_b) + η\*_b(Ĥ) ],  η\*_b(g) = β̂\*_b(g) − β̂(g),
   and estimates its variance by the infinitesimal jackknife (Efron 2014; Wager–Hastie–Efron
   2014), treating the correction as a bagged estimator.

3. **The term Harrell omits.** The optimism correction of Harrell, Lee & Mark (1996) uses only
   the first discrepancy, η\*_b(Ĥ\*_b) — the bootstrap-vs-observed effect at the *re-selected*
   subgroup, i.e. the selection optimism. León retains a *second* term, η\*_b(Ĥ), the
   discrepancy at the *observed selected* subgroup, which corrects the estimator's own
   finite-sample bias there and makes the corrected estimator a proper bagged estimator (it is
   also the term that supplies the IJ residual). Naming this gap explicitly is the entry point.

4. **This paper.** We recast León's procedure as a **multiplier (wild / martingale) bootstrap
   of subgroup influence functions**. The recasting (i) makes the inference *procedure-agnostic*,
   (ii) extends it from Cox to any **difference-in-natural-parameter (DINA)** estimand
   (Gao & Hastie 2025) — log-HR, log-OR, log-IRR, mean difference, (iii) yields a fast,
   refit-free computation that equals León's bootstrap to first order, and (iv) makes precise
   what "procedure-agnostic" *requires* — and therefore how it *constrains* model-based
   selectors.

---

## What "procedure-agnostic" means — and how it limits model-based selectors

State the boundary plainly, because it is the paper's most important honest point.

- **The wrapper condition.** The de-biasing and IJ variance attach to any procedure whose
  **selection** is a functional of the per-candidate *inferential effects* over a family,
  Ĥ = 𝒮({β̂(g) : g ∈ ℱ}). Tier-1 (literal re-run) realizes this for *any* re-runnable
  procedure; Tier-2 (multiplier) realizes it *exactly to first order* whenever selection ranks
  on β̂(g).

- **The limitation on model-based inputs.** Model-based identifiers — **GRF** (ranks regions
  by a mean doubly-robust score), **DINA** (ranks by subgroup-mean per-patient effect τ̂), and
  **Lasso** (selects a rule along a penalized path) — may use their internal model to *propose
  candidates*, but their **native selection statistic is not the Cox/GLM effect**. Expressed
  through the resampling de-biasing, their selection must be re-ranked on β̂(g)
  (`select_statistic = "effect"`); this demotes the model to a **candidate generator** and
  makes the selection layer identical to forest search. What is given up is the native ranking's
  own optimality (GRF honesty, DINA cross-fit efficiency, Lasso sparsity) *as the selection
  criterion* — in exchange for an exactly correctable, validly covered effect.

- **Why this is forced, not chosen** (full argument in the derivation note): the bias term is an
  *effect* discrepancy, so the bootstrap re-selection that the bias averages over must also be
  governed by the effect; a native-statistic re-selection picks different winners on the same
  resample and averages the optimism over the wrong distribution.

---

## Section spine

1. **Introduction** (above): post-selection bias; León's two-term correction and the Harrell
   gap; the procedure-agnostic, DINA-general, fast-computation proposal; the limitation.
2. **Subgroup effects as DINA estimands.** Cox/GLM as differences in natural parameters; the
   M-estimator / `dfbeta` influence; one additive scale.
3. **León's correction and its IJ variance, restated generally** (not procedure-specific);
   the two-term decomposition = selection optimism + estimation bias; the Harrell contrast.
4. **The multiplier-bootstrap recasting.** Resample = reweighting; refit = `dfbeta` dot
   product; shared multiplier = correlated competition; centered-Poisson = bootstrap weight.
   The **correspondence proposition** (Q1): Tier-2 = Eq. 7 (point + IJ), Ĥ and Ĥᶜ.
5. **Selection-criterion alignment (Q2 — the contribution).** Why native DINA/GRF/Lasso
   selection breaks the correspondence; the precise condition (rank on β̂(g)); the worked
   example; procedures as candidate generators; **conditional-on-family** exactness.
6. **Scope (Q3, Q5).** Fixed vs data-driven family; the conditional estimand the paper
   certifies; the unconditional gap as an empirical quantity; native-selection faithful
   linearization as future work.
7. **Simulation (Q4).** León Table 2 (M1–M3; targets θ̂(H)/θ‡/θ†; bias, coverage, SD
   calibration) + Tier-2 block + paired Tier-1↔Tier-2 panel + the family-regeneration
   experiment + GLM/OR closure.
8. **Applications.** GBSG/Cox and ACTG175/OR, procedure-agnostic demonstrations.
9. **Discussion.** Procedure-agnostic use in practice; small-nₛ linearization limits;
   conditional vs unconditional targets; RAL beyond DINA.

*(Sections 2 and the closed-form consistency screen can be condensed or referenced to the
companion methods note, since they support but are not central to the inference thesis.)*

---

## Headline results to feature

- **Correspondence (point):** β̃(Ĥ) = β̂(Ĥ) − bias_sel − bias_fix equals León Eq. 7 to first
  order via η\*_b(g) = Σ_{i∈g}(K\*_bi − 1)·dfbeta_{g,i} = D_g(b) (Ĥ and Ĥᶜ).
- **Correspondence (interval):** León's IJ variance from the same multiplier draws; the
  de-biased IJ interval is wider than the naive robust-SE interval (carries selection variance).
- **Alignment proposition (Q2):** the correspondence holds for DINA/GRF/Lasso *iff* selection
  ranks on β̂(g); the worked example shows the per-draw winner flip that breaks it otherwise.
- **Validation in hand:** ACTG175/OR (T2 2.91 vs T1 2.60); GBSG/Cox (T2 HR 1.81 (0.71, 4.67)
  vs T1 1.93 (0.88, 4.24); IJ SE 0.48 vs robust 0.36; ≈637× faster). Point tracks; T2 IJ SE
  ~20% conservative at small nₛ.

---

## Scope & limitations (own them up front)

- The certified claim is **conditional on the proposed candidate family**; the unconditional
  (family-regenerating) target is empirical, not proven.
- The linearization is loosest at **small subgroup size nₛ** (inherited from the closed-form
  screen's regime).
- Faithful fast de-biasing of a **native (non-effect) selection rule** is out of scope; it
  needs the native↔effect cross-influence under shared multipliers.

---

## Prior work

León et al. (2024) — origin; Harrell, Lee & Mark (1996) — optimism correction (the omitted
term); Efron (2014), Wager–Hastie–Efron (2014) — IJ-for-bagging variance; Lin–Wei–Ying (1993),
Parzen–Wei–Ying (1994), Jin–Ying–Wei (2001), Dobler & Pauly (2017) — martingale / wild
bootstrap of the estimating function; Gao & Hastie (2025) — DINA estimand; Wager–Athey (2018),
Athey–Tibshirani–Wager (2019), Foster et al. (2011) — GRF / virtual-twins identifiers the
inference can wrap.

---

## Open decisions

- Title emphasis: "procedure-agnostic post-selection inference" vs "multiplier-bootstrap
  recasting of León."
- Depth of the DINA-scale / consistency-screen material (full sections vs reference to the
  methods note).
- GLM extension as co-equal contribution vs generalization remark.
- Target venue and length.
