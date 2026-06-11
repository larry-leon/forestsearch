# Hand-off — paper drafting: procedure-agnostic post-selection inference

**Phase.** Methodology, theory, the rigor-first outline, the full mathematical derivation, and an
overview deck are all complete and validated. The next phase is **drafting the actual paper
sections**. This hand-off captures the state so a fresh chat can begin writing immediately.

**Author / setting.** Larry León, author/maintainer of `forestsearch` (CRAN v0.1.0; v0.2.0 dev on
`feature/glm-extension`) and `weightedsurv`. The paper formalizes and generalizes the
post-selection de-biasing in León et al. 2024 (*Stat Med*, DOI 10.1002/sim.10163). Solo
simulation work; Pop!_OS ~127 cores; RStudio + Quarto + GitHub Desktop.

---

## 1. What to bring into the new chat

Upload these (they are the canonical current versions; the new chat will not otherwise have
them). Confirm each is the latest before editing — the most recent in-chat artifact supersedes the
project KB for any touched file.

**Drafting sources (primary):**
- `PAPER_OUTLINE.md` — the rigor-first outline; the spine for drafting (questions ranked Q1–Q5,
  intro framing, section spine, headline results, scope, prior work).
- `selection_alignment_derivation.qmd` — the full mathematical derivation (setup → León Eq. 7 →
  Tier-2 linearization → effect-vs-native selection → alignment proposition → DINA statistic
  explicit → worked logistic/OR example → IJ → conditional-on-family scope; plus
  `@sec-multipliers`, the full treatment of the multiplier K\*(b,i)). The long form to condense
  into the methods section.
- `consistency_resampling_theory.qmd` — the methods/theory note. Contains `@sec-debias`
  (intuition for why one multiplier draw stands in for a bootstrap re-search), `@sec-debias-corr`
  (the León Eq. 7 correspondence, `@eq-leon7`), `@sec-debias-ij` (the IJ interval), and
  `@sec-debias-general` (procedure-agnostic use + selection-criterion alignment for DINA/GRF).
- `references.bib` — bibliography (León, Harrell, Efron, Wager–Hastie–Efron, LWY/PWY/JYW
  martingale-bootstrap, Gao–Hastie DINA, Wager–Athey & Athey–Tibshirani–Wager GRF, Foster VT).
- `Leon_SIM_2024.pdf` — the foundational paper (already in project knowledge).

**Supporting:**
- `ARTICLE_OVERVIEW.md` — two-layer thesis (DINA estimand layer + multiplier-bootstrap inference
  layer), contributions (1–6), section spine, prior-work positioning.
- `SIMULATION_HANDOFF.md` — the simulation plan (León 2024 Table 2 design extended with the
  Tier-2 + IJ block, a paired Tier-1↔Tier-2 panel, the family-regeneration experiment, GLM/OR
  closure).
- `paper_overview.pptx` / `paper_overview.pdf` — 13-slide overview deck (talk-ready; mirrors the
  paper's argument).

**Implementation (reference only, if a methods detail needs checking against code):**
- `fs_debias_gate.R` (the Tier-2 gate + IJ variance), `forestsearch_main.R` (SECTION 9B hook),
  `gate_estimates_table.R`.

---

## 2. The paper at a glance

**Working title.** *Procedure-agnostic post-selection inference for identified subgroups: a
multiplier-bootstrap recasting of the León (2024) de-biased estimator and infinitesimal-jackknife
interval.*

**Two-layer thesis.** (i) An *estimand* layer: subgroup effects as differences in natural
parameters (DINA; Gao–Hastie 2025), giving one additive scale across Cox/log-HR, logistic/log-OR,
Poisson/log-IRR, Gaussian/mean-difference. (ii) An *inference* layer: León's bootstrap de-biasing
+ IJ variance, recast as a multiplier (wild / martingale) bootstrap of subgroup influence
functions — procedure-agnostic, fast (refit-free), and exact to first order.

**The rigor-ranked questions (the organizing principle):**

| # | Question | Status |
|---|----------|--------|
| Q1 | Fixed family, effect-ranked selection: does Tier-2 equal León Eq. 7 (point + IJ)? | **Proven** (first-order identity, Ĥ and Ĥᶜ) |
| Q2 | Can DINA/GRF/Lasso be brought under Q1? | **Proven, conditional on the proposed family** — iff selection is re-ranked on β̂(g). *Central feasible contribution.* |
| Q3 | Validity when the family is itself data-driven (regenerated per resample)? | **Characterized, open** — fast gate targets the conditional estimand; family-generation component to be quantified vs literal Tier-1. |
| Q4 | Finite-sample operating characteristics (bias, coverage, SD calibration)? | **Empirical** — simulation (León Table 2 extended). |
| Q5 | Faithful fast de-biasing of a *native* (non-effect) selection rule? | **Out of scope** — needs the native↔effect cross-influence under shared multipliers. |

The spine is **Q1 → Q2 (the contribution), with Q3–Q5 as honestly-bounded scope.**

**Contributions.** (1) procedure-agnostic post-selection inference; (2) Cox→GLM via the DINA
scale; (3) fast, refit-free, first-order-exact computation; (4) a single-influence-matrix pipeline
(screen + family routine + de-biasing from one N×S matrix); (5) the explicit selection-criterion
condition; (6) a clean theoretical home for DINA/GRF.

---

## 3. Core technical results (so the new chat has the substance)

**Tier-1 vs Tier-2.** Tier-1 = the literal nonparametric bootstrap that re-runs the whole search
per resample (León Eq. 7 bias correction + IJ variance). Tier-2 = a fast multiplier (wild)
bootstrap that reuses per-candidate dfbeta with no refit; one shared multiplier vector per draw
broadcasts to all candidates and re-selects the winner per draw.

**The multiplier K\*(b,i).** The count of subject i in resample b; the count vector *is* the
resample. Refit = reweight (Σ_i K\*(b,i)·ψ_i = 0; exact for GLM, first-order via the dfbeta for
Cox). Only the centered weight G_i(b) = K\*(b,i) − 1 enters (baseline 1 reproduces the original
fit). Then β̂\*_b(g) − β̂(g) ≈ Σ_{i∈g} G_i(b)·dfbeta_{g,i} = D_g(b). G has mean 0, variance 1 —
the IJ. Weight laws (Multinomial / Poisson(1) / Rademacher ±1 / Gaussian / Dirichlet) agree to
first order; **Rademacher ±1 is the 50/50 split** underlying the consistency criterion;
centered-Poisson is the IJ of the bootstrap. One vector drives all S candidates → one matrix
P = Bᵀ Ξ (data in B fixed, randomness in Ξ).

**The correspondence (Q1).** Tier-2 reproduces León's Eq. 7 to first order, for the selected
subgroup Ĥ *and* its complement Ĥᶜ, with the IJ interval computed from the same draws. The only
approximation is the dfbeta linearization (exact for mean-difference, loosest at small subgroup
size nₛ). The selection map enters Tier-1 and Tier-2 identically, so selection non-smoothness is
**not** a source of Tier-1↔Tier-2 disagreement.

**León's two terms vs Harrell.** León Eq. 7 subtracts the average of two discrepancies:
η\*_b(Ĥ\*_b) (selection optimism — the Harrell–Lee–Mark 1996 term) **plus** η\*_b(Ĥ) (the
discrepancy at the observed selected subgroup — the estimator's own bias there, which Harrell
omits). The second term makes β̃ a proper bagged estimator and supplies the IJ residual. This is
the introduction's hook.

**Selection-criterion alignment (Q2).** The correspondence holds whenever Ĥ = 𝒮({β̂(g): g∈ℱ}) is
a functional of the per-candidate *inferential effects*. DINA ranks on subgroup-mean τ̂
(T_g = ā_g^⊤β̂_D, with influence dt_{g,i} = ā_g^⊤IF_D ≠ dfbeta_{g,i}); GRF ranks on a mean
DR-score; Lasso on a penalized path — none equal β̂(g). So native selection re-selects different
winners on the same draw, and an effect-based gate averages the optimism over the wrong
distribution. Re-ranking the procedure's qualifying candidates on β̂(g) (`select_statistic =
"effect"`) makes T = β̂, E = D, one influence matrix, gate exact — DINA/GRF/Lasso become
**candidate generators**. **Validity-over-bias argument:** effect-ranking selects a more extreme
subgroup (larger raw bias) but an *exactly correctable* one; native ranking's smaller bias is one
the effect correction cannot see. Exact **conditional on the proposed family**.

**Worked example (the per-draw winner flip).** Two candidates: β̂(g₁)=0.90, β̂(g₂)=1.00 (effect
ranks g₂); τ̂: T₁=0.55, T₂=0.50 (native ranks g₁). On one draw D₁=+0.05, D₂=+0.30 (effect),
E₁=+0.12, E₂=+0.02 (native): native picks g₁ and credits +0.05; the effect gate picks g₂ and
credits +0.30. The gate over-states the native procedure's optimism (over-shrinks). Under
`select_statistic="effect"` the procedure selects g₂ and the gate is the faithful linearization.

**Conditional vs unconditional family (Q3).** For an enumerated family (forest search) ℱ is fixed
and the fast gate is faithful. For DINA/GRF ℱ is data-driven (the surface/forest regenerates per
resample as ℱ\*_b); the fast gate omits the family-generation component of bias and variance.
Direction: forest-search Tier-2 IJ was mildly *conservative*; an unaccounted family-regeneration
component pushes the DINA/GRF Tier-2 variance toward *under*-coverage — net unsigned a priori,
to be measured against a literal Tier-1 that regenerates ℱ per resample.

**Validation in hand.** ACTG175/OR (T2 2.91 vs T1 2.60). GBSG/Cox: selected {er≤0 & size≤35}
n=61; de-biased HR T2 1.81 (0.71, 4.67) vs T1 1.93 (0.88, 4.24); complement 0.635 vs 0.633; IJ SE
0.48 vs robust 0.36; full family 324 vs 28 survivors; 0.34s vs 217.9s (≈637×). Points track;
Tier-2 IJ SE ~20% larger (mildly conservative at small nₛ) — to characterize in simulation.

---

## 4. The introduction framing (Larry's explicit asks)

1. Lead with León's de-biasing + IJ variance and **name the bias term Harrell does not address**
   (η\*_b(Ĥ), the estimation bias at the selected subgroup; Harrell 1996 corrects only selection
   optimism).
2. State precisely what **"procedure-agnostic"** means (the inference wraps any selection that is
   a functional of the per-candidate effects) — and **how it limits model-based inputs**: GRF,
   DINA screening, and Lasso may generate candidates, but expressed via the resampling de-biasing
   their selection must rank on β̂(g); the native ranking's own optimality is given up as the
   selection criterion in exchange for an exactly correctable, validly covered effect.
3. The feasible rigorous target is **conditional on the proposed family**; say so plainly and
   position the unconditional/native cases as scope.

---

## 5. Drafting plan (section spine, with status)

1. **Introduction** — per §4 above. *To write.*
2. **Subgroup effects as DINA estimands** — Cox/GLM as natural-parameter differences; the
   M-estimator/dfbeta influence; one additive scale. *Condense from the theory note / Gao–Hastie.*
3. **León's correction and its IJ variance, restated generally** — the two-term decomposition
   (selection optimism + estimation bias); the Harrell contrast; the bagged-estimator/IJ
   structure. *Draft from `@sec-debias` + `selection_alignment_derivation.qmd` §León.*
4. **The multiplier-bootstrap recasting** — resample = reweighting; refit = dfbeta dot product;
   shared multiplier = correlated competition; centered-Poisson = bootstrap weight. **The
   correspondence proposition (Q1).** *Draft from `@sec-multipliers` + `@sec-debias-corr` +
   `@sec-debias-ij`.*
5. **Selection-criterion alignment (Q2 — the contribution)** — why native DINA/GRF/Lasso breaks
   the correspondence; the precise condition; the worked example; procedures as candidate
   generators; conditional-on-family exactness. *Condense from `selection_alignment_derivation.qmd`
   + `@sec-debias-general`.*
6. **Scope (Q3, Q5)** — fixed vs data-driven family; the conditional estimand certified; the
   unconditional gap as an empirical quantity; native-selection faithful linearization as future
   work. *Draft from the scope subsections.*
7. **Simulation (Q4)** — León Table 2 + Tier-2 block + paired Tier-1↔Tier-2 panel +
   family-regeneration experiment + GLM/OR closure. *From `SIMULATION_HANDOFF.md`; needs the runs.*
8. **Applications** — GBSG/Cox, ACTG175/OR; procedure-agnostic demonstrations. *Numbers in hand.*
9. **Discussion** — procedure-agnostic use in practice; small-nₛ linearization limits; conditional
   vs unconditional targets; RAL beyond DINA. *To write.*

Sections 2 and the closed-form consistency screen can be condensed or referenced to a companion
methods note (they support but are not central to the inference thesis).

---

## 6. Open decisions and remaining work

- **Title emphasis:** "procedure-agnostic post-selection inference" vs "multiplier-bootstrap
  recasting of León."
- **Depth of the DINA-scale / consistency-screen material:** full sections vs reference to the
  methods note.
- **GLM extension:** co-equal contribution vs generalization remark.
- **Target venue and length.**
- **Simulation (Q3 / Q4):** the family-regeneration experiment and the operating-characteristics
  study still need to be run (production scale ~1000 sims, B=300; ~12 h wall-clock projected on
  ~121 workers). This is the main empirical gap before submission.

---

## 7. Conventions to carry forward

- **Style/build:** tidyverse style; roxygen2 with markdown ON (literal `% < > &`, never
  Rd-escape; `@section` titles plain); CRAN compliance; Quarto-first two-file pattern; Conventional
  Commits.
- **Deliverables:** always provide the files containing functions/material under discussion as
  downloadable artifacts; **every response that creates or modifies a file ends with
  `present_files`.**
- **Scope discipline:** modify only what is explicitly requested; ask for / read actual file
  contents before editing — never infer; no opportunistic fixes; do not editorialize on
  methodology or suggest changing covariate sets/goals unless asked.
- **Source of truth:** the most recent in-chat artifact supersedes the project KB for any touched
  file; uploads at `/mnt/user-data/uploads` are read-only (copy to outputs before editing).
- **Method facts (do not relitigate):** H2 is the canonical Tier-1 estimate; `consistency_method =
  "resample"` throughout (do not mix with `"split"` in one study); `sg_focus = "eff"` normalizes
  to `"hr"`, use `reselection = "maxcons"` for faithful Tier-1/Tier-2 correspondence;
  `devtools::install()` (not `load_all()`) so parallel workers see changes; classification metrics
  = sensitivity/specificity/PPV/NPV (not Jaccard); benefit notation Ĝ/Ĝᶜ deferred to the benefit
  vignette.
- **Slide/Quarto rendering gotchas (learned this arc):** β̂/β̃/τ̂ combining diacritics render in
  Calibri but **drop in Cambria** under the LibreOffice QA renderer — set math-bearing slide text
  in Calibri (Ĥ = U+0124 precomposed is safe). No Unicode subscript b or g exists — use small-font
  runs or functional notation D(g), K\*(b,i), dfbeta(g,i). pptxgenjs: no `#` in hex; fresh shadow
  object per shape; no accent stripes/edge bars; equations as Unicode text (no LaTeX).
- **Notation:** β̂(g) within-subgroup effect; Ĥ/Ĥᶜ selected subgroup and complement; D_g(b) effect
  perturbation; η\*_b(g) effect discrepancy; T_g native statistic; E_g(b) native perturbation;
  𝒮 selection map; ℱ candidate family; Ξ multiplier matrix; B_eff influence matrix.

---

## 8. Suggested first message for the new chat

> Starting to draft the paper. Files attached: PAPER_OUTLINE.md (spine),
> selection_alignment_derivation.qmd (full math), consistency_resampling_theory.qmd (methods
> note), ARTICLE_OVERVIEW.md, SIMULATION_HANDOFF.md, references.bib, Leon_SIM_2024.pdf, and the
> overview deck. Let's begin with the **Introduction** — lead with León's two-term de-biasing +
> IJ variance and the bias term Harrell (1996) omits, then state what "procedure-agnostic" means
> and how it constrains GRF/DINA/Lasso, and frame the conditional-on-family correspondence as the
> feasible rigorous target. Draft it as a Quarto section I can drop into the manuscript.

---

## 9. Key references

León et al. 2024 (*Stat Med*, 10.1002/sim.10163) — origin. Harrell, Lee & Mark 1996 — optimism
correction (the omitted term). Efron 2014; Wager, Hastie & Efron 2014 — IJ-for-bagging variance.
Lin–Wei–Ying 1993; Parzen–Wei–Ying 1994; Jin–Ying–Wei 2001; Dobler & Pauly 2017 — martingale /
wild bootstrap of the estimating function. Gao & Hastie 2025 (*Biometrics*) — DINA estimand.
Wager & Athey 2018; Athey, Tibshirani & Wager 2019 — GRF. Foster, Taylor & Ruberg 2011 — virtual
twins. CAPITAL (Cai, Lu et al. 2022) — DINA 2-D extension (planned).
