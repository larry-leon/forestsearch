# Proof of concept — interval constructions for the MR de-biased estimate

Date: 2026-09-05. Status: exploratory; anything adopted from here is a method proposal for the gate, not a queue item.

## Question

The paper's interval uses the two-term residual r_b = (bias_sel + bias_fix) − D_{Ĥ*_b}(b) − D_Ĥ(b) in the infinitesimal jackknife (Eq. 11–12). Its limit variance is (p + e_Ĥ)ᵀΣ(p + e_Ĥ), which is 4σ² when one candidate dominates — a doubled standard error exactly where there is least to correct. Which constructions fix the inflation where it is largest without breaking the tie regime, and at what cost?

## Setup

Normal-means model (Y = τ·A + ε, ε ~ N(0,1)), within-candidate mean differences, dfbeta influences, centred-Poisson multipliers, B = 4,000 draws, pure argmax, fixed family. The influence linearization is exact here, so differences between constructions are due to the constructions. The point estimate is always the paper's β̃ = β̂(Ĥ) − bias_sel − bias_fix. 1,000 replicates per scenario (S5: 400); coverage Monte Carlo SE ≈ 0.007 (S5: 0.011). SE of a candidate effect ≈ 0.141 (n_g = 200).

Scenarios: S1 K = 1 (no competition; the complement-subgroup case); S2 K = 2 disjoint separated by 4 SE; S3 K = 2 separated by 1.5 SE; S4 K = 2 tie; S5 K = 10 disjoint tie; S6 K = 10 nested (n = 400, unequal candidate sizes) tie.

Constructions (all around β̃ unless stated):

| label | interval |
|---|---|
| naive | β̂(Ĥ) ± 1.96 σ̂_{D,Ĥ} — reference, not de-biased |
| paper | ± 1.96 √V̂, two-term residual (current gate) |
| bag | ± 1.96 √V̂_bag, residual r_b = bias_sel − D_{Ĥ*_b}(b) (winner term only; the total-derivative IJ) |
| bagfloor | as bag with V = max(V̂_bag, σ̂²_{D,Ĥ}) |
| fieldSE | ± 1.96 sd(Λ*), Λ* = ζ*_{G(β̂+ζ*)} − m̂(β̂+ζ*), ζ* ~ N(0, Σ̂), m̂ by inner Monte Carlo on the plug-in field |
| fieldQ | [β̃ − q_{.975}(Λ*), β̃ − q_{.025}(Λ*)] — basic-bootstrap form, carries the bias adjustment and the non-Gaussian shape |
| tieQ | as fieldQ with Λ* simulated on an exact-tie field (least-favourable bias) |

Everything the field constructions need — the effect vector β̂ and Σ̂ = B_effᵀB_eff — the gate already computes; the extra cost is K × 500 × 300 Gaussian draws per analysis.

## Results

Sanity: retained optimism 0.296 (S4) and 0.293 (S5) against the tie constant 0.293; 0.368 on the nested family — the same excess seen on the real §5.2 replicates (0.35–0.37), so that excess is a property of unequal candidate variances under nesting, not of the survival model. Sign reversal at 1.5 SE → 4 SE separation reproduced.

Two-sided 95% coverage of θ_Ĥ (half-width in parentheses; naive half-width 0.277):

| scenario | paper | bag | bagfloor | fieldSE | fieldQ | tieQ |
|---|---|---|---|---|---|---|
| S1 K=1 | 1.000 (0.554) | 0.946 (0.277) | 0.948 (0.279) | 0.945 (0.276) | 0.942 (0.274) | 0.945 |
| S2 sep 4 SE | 0.999 (0.532) | 0.902 (0.258) | 0.945 (0.277) | 0.948 (0.288) | 0.943 (0.286) | 0.884 |
| S3 sep 1.5 SE | 1.000 (0.473) | 0.873 (0.215) | 0.949 (0.276) | 0.948 (0.280) | 0.931 (0.276) | 0.911 |
| S4 K=2 tie | 1.000 (0.457) | 0.890 (0.204) | 0.969 (0.277) | 0.975 (0.274) | 0.952 (0.270) | 0.949 |
| S5 K=10 tie | 0.998 (0.331) | 0.642 (0.107) | 0.955 (0.276) | 0.975 (0.264) | 0.962 (0.260) | 0.952 |
| S6 nested tie | 1.000 (0.551) | 0.934 (0.266) | 0.961 (0.311) | 0.960 (0.302) | 0.949 (0.302) | 0.948 |

SE / empirical SD of β̃: paper 1.57–2.00; bag 0.51–1.00; bagfloor 0.97–1.31; fieldSE 0.99–1.26.

One-sided 95% lower bound (the Guo–He convention), coverage: paper 0.98–1.00; bag 0.63–0.95; bagfloor 0.925–0.951; fieldQ 0.940–0.948. Half-width in SE units: paper 1.9–3.3; bagfloor 1.65; fieldQ 1.63–1.87; Guo–He on the real replicates ≈ 1.4.

Check on the real Guo–He replicates (11 cells, survival, B = 5,000): the floor construction can be read directly from the stored columns as β̃ ± z·σ̂_{D,Ĥ}. One-sided coverage 0.922–0.940 (paper 0.97–1.00; G&H 0.946–0.958), two-sided 0.922–0.962, half-width 0.30 against the paper's 0.36–0.59 and G&H's 0.20–0.31. The under-coverage sits in the moderately separated t35 cells (2.2–2.7 SE), where β̃'s SD runs 10% above the naive SE — the m(μ+ζ) term's variance, which the field calibration captures and the floor does not.

## Reading

1. The winner-only IJ (bag) is the correct first-order variance and is exact when one candidate dominates — S1 and S2 return the naive interval. It under-covers in the local regime (down to 0.64 at K = 10) because the delta method does not hold there, which is what Theorem 2(ii)'s non-Gaussian limit already says.
2. Flooring it at the naive variance (bagfloor) repairs the local regime at zero cost: two-sided coverage 0.945–0.969 everywhere in the PoC, half-width at the naive level, i.e. 40–50% narrower than the current interval with no loss of coverage. Its weakness is the moderately separated regime on the real replicates (0.92–0.94), where β̃ is slightly noisier than β̂ — a one-to-three-point shortfall in exchange for removing a factor of 1.7–2 in width.
3. The plug-in field calibration (fieldSE / fieldQ) is the principled version: SE/SD 0.99–1.26, two-sided 0.93–0.975, one-sided 0.94–0.95, and it carries the retained-bias adjustment into one-sided bounds. Its residual weakness is S3 (0.931 two-sided), consistent with the plug-in field being optimistic at the winner so that the simulated competition is more separated than the truth. A shrunk field — de-biasing the winner's entry by bias_sel before simulating, or a Guo–He-style shrinkage of the field — is the natural next step and would bring this into the Guo–He family of methods without leaving the influence framework.
4. The tie-field variant is wrong under separation (0.88) and is dropped.

## Recommendation as a method proposal

- Minimal change to the gate: drop D_Ĥ(b) from the residual and floor V̂ at σ̂²_{D,Ĥ}. One-line change in arithmetic; the current V̂ can be retained as an optional `ci_method` for reproducibility of the paper's tables.
- Preferred: add a `ci_method = "field"` that simulates Λ* on (β̂, Σ̂) with a shrunk winner entry and reports both the SE-type and quantile-type intervals; evaluate on the Guo–He replicates (the bundles hold everything except Ξ, so the gate needs to run once more per replicate — MR cost only, ~2–10 s).
- Either way, report the paper's interval, the alternative, and the naive interval side by side in the Stage 3 record so the width/coverage trade is visible.

Caveats: normal-means model with exact linearization; pure argmax, no screen or consistency; K ≤ 10; the complement subgroup is represented only by the K = 1 mechanism; 1,000 replicates (S5: 400).

Files: `poc_ci.py` (simulation), `poc_ci_results_2026-09-05.csv` (all metrics), `mr_vs_guohe_waldfloor_check_2026-09-05.csv` (floor check on the real replicates).
