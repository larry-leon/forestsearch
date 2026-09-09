# REVIEW — E1 (campaign `e1stud`): the R1 "field-s" construction — formal close

**Date:** 2026-09-08. Reviewer: the Linux MR-field chat. Chat-side review record; companion to `REVIEW_partB_banddial_2026-09-08.md` and `REVIEW_partsAC_close_2026-09-08.md`.
**Sources (numbers verbatim):** `REPORT_field_studentize_e1_stage1_2026-09-08.md` (3880022e), `REPORT_field_studentize_e1_2026-09-08.md` (…57e00f69), `PROPOSAL_complement_field_scale_2026-09-08_v2.md` (§2, §5, §6), `REPORT_field_studentize_e0_2026-09-08.md`, `REPORT_complement_variance_2026-09-07.md`. Commit range 9a6e21cc..57e00f69, unpushed. Conventions: winner-only/winner-floor excluded; bounds by location; Wilson intervals; marginal and error-scale SDs both shown (Part A, A0).

## 1. Verdict

- **Execution: PASS.** Stage 0 quotes match HEAD; G2a (both knobs off) 131/131 identical, 158 = 136 + 22 columns, new ones NA; G2b (scale on) 131/131 identical, `_s` invariants and `joint_s` γ in range, R1-vs-global 0.9925–1.0012; G2c (e0stud regression) 135/135 identical including the four scale columns; Gate 1 projection 2 h 52 m vs 2 h 53 m realized; six of six cells Gate 2 PASS at 131/131 to the committed comparators; harm block, unscaled complement, joint pair and identification unchanged by pairing. Two `R/` files as named; classification "adds code; byte-identical defaults; new outputs only on the enabled path" holds.
- **Invariance argument, upgraded by the record:** `sel ∈ winset`, so the multiplier stage already attempted the selected complement's fit; the relocated ensure-fit can only populate `Bc[, sel]` when that fit failed, and then `bdc` is NA and the function has already returned. Provably inert, not merely gate-inert.
- **Pre-registered success criterion (proposal §6 E1: field-s Ĥᶜ one-sided upper coverage ≥ 0.92 with Wilson support in all four band cells): NOT MET.** Point estimates 0.912 / 0.919 / 0.921 / 0.919; Wilson lower limits 0.899–0.909; one point estimate ≥ 0.92.
- **Construction validated on its own terms:** SE/error-SD 0.96–0.97 (field 0.89–0.92; naive 0.99–1.00); λ-SDᶜ/naive SE 0.997–1.002 (field 0.917–0.947); per-replicate corr(se_s, naive SE) 0.94–0.96 vs 0.34–0.42; |Ĥ|/|H| shape flattened in all four band cells; never worse than field in total in any cell.

## 2. Standard table 1 — constructions per cell (identical replicates)

Columns: bias (log scale, vs the conditional target), marginal SD (across replicates, carries target spread), error SD (of the de-biased error), mean SE, r = SE / marginal SD, SE / error SD, one-sided coverage on the exposed side [Wilson], two-sided coverage.

### Ĥ block (lower bound) — unchanged between field and field-s by construction

| Cell | Construction | bias | marg SD | err SD | SE | r | SE/errSD | 1s cov | 2s cov |
|---|---|---|---|---|---|---|---|---|---|
| ε 0.20 HR 1.50 | naive | +0.430 | 0.238 | 0.281 | 0.250 | 1.049 | 0.887 | 0.484 [0.462, 0.506] | 0.619 |
| ε 0.20 HR 1.50 | field | −0.056 | 0.299 | 0.314 | 0.324 | 1.081 | 1.029 | 0.974 [0.967, 0.981] | 0.925 |
| ε 0.20 HR 1.50 | IJ two-term | −0.003 | 0.279 | 0.303 | 0.360 | 1.289 | 1.186 | 0.985 [0.979, 0.990] | 0.981 |
| ε 0.20 HR 1.75 | naive | +0.365 | 0.246 | 0.294 | 0.243 | 0.985 | 0.825 | 0.573 [0.551, 0.595] | 0.694 |
| ε 0.20 HR 1.75 | field | −0.078 | 0.312 | 0.329 | 0.323 | 1.035 | 0.983 | 0.970 [0.962, 0.977] | 0.911 |
| ε 0.20 HR 1.75 | IJ two-term | −0.037 | 0.291 | 0.317 | 0.356 | 1.227 | 1.125 | 0.985 [0.979, 0.989] | 0.972 |
| ε 0.30 HR 1.50 | naive | +0.299 | 0.224 | 0.244 | 0.213 | 0.949 | 0.873 | 0.633 [0.612, 0.654] | 0.735 |
| ε 0.30 HR 1.50 | field | −0.093 | 0.270 | 0.267 | 0.298 | 1.103 | 1.116 | 0.980 [0.973, 0.986] | 0.937 |
| ε 0.30 HR 1.50 | IJ two-term | −0.068 | 0.254 | 0.258 | 0.334 | 1.311 | 1.294 | 0.995 [0.991, 0.997] | 0.987 |
| ε 0.30 HR 1.75 | naive | +0.254 | 0.239 | 0.254 | 0.210 | 0.879 | 0.825 | 0.680 [0.660, 0.700] | 0.777 |
| ε 0.30 HR 1.75 | field | −0.097 | 0.286 | 0.277 | 0.297 | 1.038 | 1.074 | 0.976 [0.968, 0.982] | 0.925 |
| ε 0.30 HR 1.75 | IJ two-term | −0.082 | 0.271 | 0.268 | 0.333 | 1.231 | 1.243 | 0.994 [0.990, 0.997] | 0.978 |
| maxSG HR 1.75 | naive | +0.046 | 0.063 | 0.095 | 0.133 | 2.118 | 1.400 | 0.947 [0.936, 0.956] | 0.974 |
| maxSG HR 1.75 | field | −0.005 | 0.087 | 0.131 | 0.126 | 1.455 | 0.961 | 0.947 [0.936, 0.956] | 0.956 |
| maxSG HR 1.75 | IJ two-term | +0.002 | 0.080 | 0.123 | 0.221 | 2.773 | 1.799 | 0.999 [0.997, 1.000] | 0.996 |
| minSG HR 1.75 | naive | +0.326 | 0.217 | 0.272 | 0.320 | 1.476 | 1.177 | 0.804 [0.786, 0.821] | 0.881 |
| minSG HR 1.75 | field | +0.022 | 0.284 | 0.326 | 0.315 | 1.110 | 0.967 | 0.936 [0.924, 0.946] | 0.922 |
| minSG HR 1.75 | IJ two-term | +0.074 | 0.258 | 0.305 | 0.391 | 1.514 | 1.280 | 0.971 [0.963, 0.977] | 0.985 |

### Ĥᶜ block (upper bound)

| Cell | Construction | bias | marg SD | err SD | SE | r | SE/errSD | 1s cov | 2s cov |
|---|---|---|---|---|---|---|---|---|---|
| ε 0.20 HR 1.50 | naive | −0.123 | 0.156 | 0.143 | 0.143 | 0.912 | 0.997 | 0.781 [0.762, 0.798] | 0.849 |
| ε 0.20 HR 1.50 | field | −0.031 | 0.159 | 0.147 | 0.135 | 0.852 | 0.923 | 0.897 [0.883, 0.910] | 0.921 |
| ε 0.20 HR 1.50 | **field-s** | −0.031 | 0.159 | 0.147 | 0.142 | 0.896 | 0.970 | **0.912 [0.899, 0.924]** | 0.938 |
| ε 0.20 HR 1.50 | IJ two-term | −0.044 | 0.157 | 0.145 | 0.257 | 1.636 | 1.776 | 0.995 [0.991, 0.997] | 1.000 |
| ε 0.20 HR 1.75 | naive | −0.102 | 0.162 | 0.144 | 0.143 | 0.884 | 0.991 | 0.826 [0.809, 0.842] | 0.880 |
| ε 0.20 HR 1.75 | field | −0.019 | 0.165 | 0.149 | 0.136 | 0.822 | 0.911 | 0.912 [0.899, 0.924] | 0.918 |
| ε 0.20 HR 1.75 | **field-s** | −0.019 | 0.165 | 0.149 | 0.143 | 0.868 | 0.961 | **0.919 [0.907, 0.931]** | 0.933 |
| ε 0.20 HR 1.75 | IJ two-term | −0.030 | 0.163 | 0.147 | 0.258 | 1.586 | 1.761 | 0.995 [0.991, 0.998] | 0.999 |
| ε 0.30 HR 1.50 | naive | −0.137 | 0.162 | 0.153 | 0.153 | 0.944 | 1.001 | 0.774 [0.755, 0.792] | 0.854 |
| ε 0.30 HR 1.50 | field | −0.031 | 0.167 | 0.156 | 0.141 | 0.845 | 0.900 | 0.904 [0.890, 0.916] | 0.916 |
| ε 0.30 HR 1.50 | **field-s** | −0.032 | 0.167 | 0.157 | 0.152 | 0.915 | 0.974 | **0.921 [0.909, 0.932]** | 0.936 |
| ε 0.30 HR 1.50 | IJ two-term | −0.047 | 0.164 | 0.154 | 0.266 | 1.625 | 1.726 | 0.996 [0.992, 0.998] | 1.000 |
| ε 0.30 HR 1.75 | naive | −0.113 | 0.164 | 0.153 | 0.152 | 0.928 | 0.996 | 0.807 [0.790, 0.824] | 0.875 |
| ε 0.30 HR 1.75 | field | −0.021 | 0.171 | 0.158 | 0.141 | 0.828 | 0.894 | 0.903 [0.890, 0.916] | 0.914 |
| ε 0.30 HR 1.75 | **field-s** | −0.022 | 0.171 | 0.158 | 0.153 | 0.894 | 0.965 | **0.919 [0.906, 0.930]** | 0.941 |
| ε 0.30 HR 1.75 | IJ two-term | −0.033 | 0.167 | 0.155 | 0.268 | 1.600 | 1.723 | 0.996 [0.993, 0.998] | 0.999 |
| maxSG HR 1.75 | naive | −0.191 | 0.393 | 0.406 | 0.355 | 0.903 | **0.876** | 0.776 [0.757, 0.794] | 0.844 |
| maxSG HR 1.75 | field | −0.041 | 0.392 | 0.395 | 0.361 | 0.921 | 0.915 | 0.925 [0.913, 0.936] | 0.940 |
| maxSG HR 1.75 | **field-s** | −0.039 | 0.393 | 0.395 | 0.358 | 0.912 | 0.906 | 0.926 [0.914, 0.937] | 0.940 |
| maxSG HR 1.75 | IJ two-term | −0.064 | 0.389 | 0.394 | 0.621 | 1.598 | 1.576 | 0.995 [0.991, 0.997] | 0.999 |
| minSG HR 1.75 | naive | −0.049 | 0.139 | 0.136 | 0.130 | 0.932 | 0.958 | 0.879 [0.864, 0.893] | 0.924 |
| minSG HR 1.75 | field | −0.013 | 0.137 | 0.134 | 0.129 | 0.940 | 0.967 | 0.929 [0.917, 0.940] | 0.940 |
| minSG HR 1.75 | **field-s** | −0.013 | 0.137 | 0.134 | 0.130 | 0.945 | 0.972 | **0.929 [0.917, 0.939]** | 0.943 |
| minSG HR 1.75 | IJ two-term | −0.019 | 0.137 | 0.134 | 0.251 | 1.830 | 1.881 | 0.997 [0.994, 0.999] | 1.000 |

## 3. Standard table 2 — across cells (six cells)

| Block | Construction | 1s cov min / mean / max | 2s cov mean | b (bias/marg SD) mean | r mean | SE/errSD mean |
|---|---|---|---|---|---|---|
| Ĥ (lower) | naive | 0.484 / 0.687 / 0.947 | 0.780 | +1.322 | 1.243 | 0.998 |
| Ĥ (lower) | field | 0.936 / 0.964 / 0.980 | 0.929 | −0.183 | 1.137 | 1.022 |
| Ĥ (lower) | IJ two-term | 0.971 / 0.988 / 0.999 | 0.983 | −0.066 | 1.557 | 1.321 |
| Ĥᶜ (upper) | naive | 0.774 / 0.807 / 0.879 | 0.871 | −0.632 | 0.917 | 0.970 |
| Ĥᶜ (upper) | field | 0.897 / 0.912 / 0.929 | 0.925 | −0.136 | 0.868 | 0.918 |
| Ĥᶜ (upper) | **field-s** | **0.912 / 0.921 / 0.929** | 0.939 | −0.137 | 0.905 | 0.958 |
| Ĥᶜ (upper) | IJ two-term | 0.995 / 0.996 / 0.997 | 1.000 | −0.209 | 1.646 | 1.741 |

## 4. The five pre-registered findings, verified

1. **Coverage.** Gains +1.6 / +0.7 / +1.8 / +1.6 points to 0.912–0.921; flips to cover 1.0–1.8% vs to miss 0.1–0.3%; Wilson lower 0.899–0.909 — **criterion not met**; field and field-s Wilson intervals overlap in every cell. IJ 0.995–0.996. Two-sided rises 0.914–0.921 → 0.933–0.941.
2. **Shape.** |Ĥ|/|H| axis flattened: top-tertile gains 2.5–3.6 points, bottom unchanged (ρᶜ ≈ 1.00–1.02 there), studentized per-replicate ratio 0.993–1.003 in every tertile. **p̂ axis still falls:** field-s T1 / T2 / T3 = 0.936 / 0.917 / 0.884, 0.954 / 0.928 / 0.877, 0.945 / 0.932 / 0.887, 0.936 / 0.937 / 0.884; the high-p̂ tertile sits at 0.877–0.887 under both constructions where ρᶜ is 1.01–1.04 and the unscaled ratio 0.97–0.99 — not a scale effect.
3. **Ends.** maxSG 0.925 → 0.926, minSG 0.929 → 0.929 (Wilson overlap). Inside maxSG a redistribution: tertiles 0.916 / 0.925 / 0.934 → 0.888 / 0.937 / 0.956; the bottom tertile (ρᶜ 0.83, studentized SD *smaller*) loses 2.8 points. minSG inert.
4. **Locations.** Band cells: field-s Ĥᶜ upper mean 0.997 / 1.010 / 0.977 / 0.981 (field 0.987 / 0.999 / 0.961 / 0.965; IJ 1.159–1.207); share < 0.85 0.156–0.214 (field 0.178–0.252; IJ ≤ 0.028); margin 0.246–0.266 log (field 0.234–0.248; IJ 0.42–0.44). Joint Bonferroni pair: `joint_s` 0.942 / 0.939 / 0.949 / 0.951–0.952 vs `joint` 0.932 / 0.933 / 0.935 / 0.935 — the ε 0.30 pairs reach the 0.95 target; ends unchanged (0.939 vs 0.941–0.942; 0.930 vs 0.929–0.930). Harm-side locations identical by construction.
5. **R1 vs global.** 0.996–0.999 (SD 0.0035–0.0057) in band cells, 0.9999 at minSG, **0.969 (SD 0.032, min 0.81) at maxSG** where the winner-scale CV is 0.26 — the per-draw and global forms coincide wherever winners' scales are homogeneous and separate only under the pure-size pick.

## 5. Where the residual sits (the finding of the campaign)

- **From variance to location.** The Gaussian reference for field-s (0.900–0.910 vs observed 0.912–0.921) decomposes the shortfall from 0.95 into ≈ 0.6 points from the remaining SE/SD deficit (SE/err-SD 0.965 → Φ(1.588) ≈ 0.944) and ≈ 3 points from the residual **bias** of the corrected complement estimate (−0.019 to −0.032 log-HR, −0.12 to −0.19 marginal-SD units; still optimistic).
- **Localized to the stable-pick stratum.** Field-s covers 0.936–0.954 in the low-p̂ tertile and 0.877–0.887 in the high-p̂ tertile, under both constructions. **The analysis-time flag reverses:** for the unscaled field the unstable pick (low p̂) was the danger; for field-s the danger is the stable pick (high p̂).
- **Hypothesis (not established):** when the observed winner re-selects itself, the field's correction — the mean of Λ*ᶜ — collapses toward zero, while the complement's real optimism does not; the two-term correction, drawn from the same field, shrinks with it. If so, the residual is under-correction concentrated where selection variability is smallest — the "conditioning on the observed field" issue the proposal's §2(iii) set out of scope, and adjacent to the handoff §5 level-dimension item.
- **The decisive check is zero compute:** on the e1stud (or nb20) bundles, the *means* of a (naive error), c (two-term correction), e (de-biased error) and `fld_Hc_lam_mean` (the field's own correction) by p̂ tertile. Part A printed variances only; the document reads any bundle set via `FS_SUMCV_GLOBS`, so this is one added chunk. Prediction under the hypothesis: mean(a) stays ≈ −0.10 in T3 while both corrections shrink, leaving mean(e) ≈ −0.05 to −0.09 (Φ(1.645 − 0.06/0.15) ≈ 0.89 reproduces the observed 0.88); alternative: mean(a) → 0 in T3 and the story is elsewhere.
- **Caveat from maxSG (outside the adoptable range).** The premise "error ≈ naive SE" (Part A: 0.98–1.01 in the band regimes) fails under the pure-size pick: naive SE / error SD **0.876**. Field-s pins to the naive SE and therefore inherits its calibration — exactly as good as the naive SE, no better; the bottom-tertile loss at maxSG is that inheritance made visible. A documentation caveat for the construction, not a blocker.

## 6. Decision items (recommendations marked)

- **D-1 Adoption of field-s.** Options: (a) promote field-s to the documented one-sided complement product, field kept as the unscaled comparator, IJ two-term as the conservative option; (b) report field and field-s side by side until the location residual is understood; (c) keep the documented rule as is, field-s as an add-beside diagnostic. **Recommendation: (a), with the shortfall stated** — ≈ 3 points below 0.95 at n = 500 in the band regimes, concentrated in the high-p̂ stratum; field-s dominates field on every metric that matters (coverage, two-sided, joint pair, locations, error-scale calibration), is the principled scaling, costs nothing, and adopting it does not claim the complement problem solved. The documented rule becomes: field-s bound with its shortfall stated against p̂ (high p̂ = caution) and the per-replicate `se_field_s`/naive SE identity; two-term IJ when a conservative statement is required. Larry's call.
- **D-2 The zero-compute localization** (§5, the stratified means; document-level; one CC session, no compute). All five reviews of this line are closed, so the standing "no new tasks" condition is lifted. **Recommendation: yes, next, after the Mac merge.** The task document is ready to draft on Larry's word.
- **D-3 E2.** Not proposed. Any location repair waits on D-2's answer; the proposal's own scope statement stands (the residual at concentrated picks is out of the scale repair's scope).
- **D-4 Housekeeping.** Push 9a6e21cc..57e00f69 before the Mac merge; the maxSG caveat and the p̂-flag reversal go into the next handoff.

## 7. Status of the line

Reviews closed: Part A, Part B, Part C, E0, E1. Constructions in the package: `field` (documented), `field-s` (add-beside, defaults off, validated), `field_decompose` diagnostics (add-beside, defaults off). Interim documented rule unchanged until D-1 is decided.
