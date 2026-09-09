# REVIEW — Part B (campaign `banddial`): the identifier dial on the J = 10 harm cells

**Date:** 2026-09-08. Reviewer: the Linux MR-field chat (per `HANDOFF_mr_field_linux_2026-09-07.md`). Chat-side review record; not a repo record unless Larry routes it there.
**Sources (all numbers verbatim from these):** `REPORT_banddial_2026-09-07.md`, `REPORT_banddial_gate2_2026-09-07.md`, task `dev/tasks/TASK_complement_variance_banddial_2026-09-07.md` (4b245516), `REPORT_nb20_2026-09-07.md` for the paired reference settings. Commit range for all three parts: 4b245516..cef19e9d (unpushed by CC beyond c133902e). Derived quantities (deltas, exchange rates) are marked; nothing else is computed outside the records.
**Conventions:** bounds read by location, never as significance at 1.0; Wilson 95% intervals; "SD units" and SD ratios are the reports' marginal-SD convention, with Part A's error-scale values beside them where the record prints both; winner-only and winner-floor excluded throughout.

## Verdict

**PASS — Part B review closed.** Six cells, all Gate 2 PASS, 2 h 36 m against the 5 h ceiling, none deferred or dropped; same DGM draws and same family K as nb20 arm A asserted per cell; nesting minSG ⊆ ε 0.10 ⊆ ε 0.20 ⊆ ε 0.30 ⊆ maxSG holds on the full run; dial table, constructions tables, and across-settings table internally consistent on every traced value; γ ∈ [0.025, 0.030] with the joint invariants inside gates; T-2's ε = 0.40 condition fired (maxSG and ε 0.30 differ materially) and was correctly not run uninstructed. Scope as kicked off (three settings); the maxSG/minSG gate maps re-quoted from installed source at Stage 0 and the template's focus-guard widening made explicitly, covered by the exact knob-inert identity against nb20 arm A.

## 1. The frontier (identification; paired replicates, n = 1,999 per cell; true |H| = 153)

Values HR 1.50 / HR 1.75 where they differ.

| Setting | mean \|Ĥ\| | \|Ĥ\|/\|H\| median | share ≥ \|H\| | sens | spec | PPV | NPV | β(Ĥ) vs 1.499 / 1.746 | β(Ĥᶜ) vs 0.721 |
|---|---|---|---|---|---|---|---|---|---|
| minSG | 62.5 / 61.8 | 0.41 / 0.40 | 0.001 | 0.248 / 0.249 | 0.929 / 0.931 | 0.603 / 0.613 | 0.737 / 0.738 | 1.163 / 1.304 | 0.874 / 0.915 |
| effMaxSG ε 0.10 | 97.3 / 99.3 | 0.60 / 0.62 | 0.064 / 0.069 | 0.450 / 0.506 | 0.917 / 0.936 | 0.697 / 0.767 | 0.794 / 0.815 | 1.252 / 1.512 | 0.835 / 0.848 |
| effMaxSG ε 0.20 | 127.1 / 128.9 | 0.79 / 0.82 | 0.264 / 0.302 | 0.590 / 0.660 | 0.893 / 0.919 | 0.706 / 0.779 | 0.837 / 0.866 | 1.256 / 1.515 | 0.808 / 0.808 |
| effMaxSG ε 0.30 | 168.9 / 166.8 | 1.09 / 1.09 | 0.634 / 0.647 | 0.737 / 0.792 | 0.838 / 0.868 | 0.682 / 0.740 | 0.883 / 0.909 | 1.228 / 1.454 | 0.779 / 0.774 |
| maxSG | 367.3 / 401.6 | 2.52 / 2.75 | 0.989 / 0.999 | 0.915 / 0.948 | 0.344 / 0.260 | 0.402 / 0.372 | 0.915 / 0.936 | 0.977 / 1.014 | 0.762 / 0.758 |

**Exchange rate, derived** (Δsensitivity / Δspecificity per step, HR 1.50 / HR 1.75): ε 0.10 → 0.20: **5.8 / 9.1**; ε 0.20 → 0.30: **2.7 / 2.6**; ε 0.30 → maxSG: **0.36 / 0.26**. minSG → ε 0.10: +0.202 sens for −0.012 spec at HR 1.50; at HR 1.75, ε 0.10 improves both (+0.257 sens, +0.005 spec) — **minSG is frontier-dominated**.

Frontier readings (the record's, confirmed): ε 0.30 is the first setting whose median region reaches the planted size, overshooting on the top third (q90 1.43–1.48); PPV peaks at ε 0.20; β(Ĥ) is flat ε 0.10 → 0.20, dips first at ε 0.30, collapses to ≈ 1.0 under maxSG (73–80% of the sample — the overall population's effect); β(Ĥᶜ) moves monotonically toward θ†(Hᶜ) = 0.721 from above, reaching 0.76 under maxSG; minSG sits on the n.min floor (61–64 on 99% of replicates). The two ends are boundary markers, not subgroup rules at this design.

## 2. Standard summary table — Ĥ block (per cell × setting × estimator)

Bias on the log-HR scale against β(Ĥ); "mSD" = marginal SD units, "eSD" = error-SD units (Part A convention), as the record prints them; r = mean SE / marginal SD; 1s = one-sided lower (the exposed side); 2s = two-sided. Absolute SDs/SEs are not printed in the record — ratios are reported as given. Field bias in mSD units printed in the record only at ε 0.30 (−0.3); field point estimate is est₂.

### HR 1.50 n500

| Setting | naive: bias log (mSD; eSD) | naive 1s cov | field: bias log | field 1s cov (r) | field 2s cov | IJ two-term 1s cov (r) | harm SD(β̃)/naive SE marg (err) |
|---|---|---|---|---|---|---|---|
| minSG | +0.384 (+2.03; +1.53) | 0.752 (0.733, 0.771) | +0.03 | 0.942 (1.18) | 0.94 | 0.974 (1.62) | 0.73 (0.90) |
| ε 0.10 | +0.548 (+2.48; +1.98) | 0.381 (0.360, 0.402) | −0.00 | 0.970 (1.12) | 0.94 | 0.978 (1.41) | 0.94 (1.06) |
| ε 0.20 | +0.430 (+1.81; +1.53) | 0.484 (0.462, 0.506) | −0.06 | 0.974 (1.08) | 0.93 | 0.985 (1.29) | 1.10 (1.19) |
| ε 0.30 | +0.299 (+1.33; +1.23) | 0.633 (0.612, 0.654) | −0.09 | 0.980 (1.10) | 0.94 | 0.995 (1.31) | 1.18 (1.19) |
| maxSG | +0.082 (+1.47; +0.90) | 0.943 (0.932, 0.952) | +0.00 | 0.949 (1.69) | 0.95 | 0.999 (3.33) | 0.46 (0.87) |

### HR 1.75 n500

| Setting | naive: bias log (mSD; eSD) | naive 1s cov | field: bias log | field 1s cov (r) | field 2s cov | IJ two-term 1s cov (r) | harm SD(β̃)/naive SE marg (err) |
|---|---|---|---|---|---|---|---|
| minSG | +0.326 (+1.50; +1.20) | 0.804 (0.786, 0.821) | +0.02 | 0.936 (1.11) | 0.92 | 0.971 (1.51) | 0.80 (0.95) |
| ε 0.10 | +0.483 (+2.12; +1.67) | 0.474 (0.452, 0.496) | −0.02 | 0.964 (1.08) | 0.93 | 0.982 (1.34) | 1.00 (1.13) |
| ε 0.20 | +0.365 (+1.48; +1.24) | 0.573 (0.551, 0.595) | −0.08 | 0.970 (1.04) | 0.91 | 0.985 (1.23) | 1.18 (1.28) |
| ε 0.30 | +0.254 (+1.07; +1.00) | 0.680 (0.660, 0.700) | −0.10 | 0.976 (1.04) | 0.93 | 0.994 (1.23) | 1.27 (1.26) |
| maxSG | +0.046 (+0.74; +0.49) | 0.947 (0.936, 0.956) | −0.01 | 0.947 (1.46) | 0.96 | 0.999 (2.77) | 0.59 (0.92) |

Readings: the field's one-sided lower bound is inside **0.936–0.980 across the whole dial** (|Ĥ| 62 → 402, two pure-size picks) — the harm-side robustness result extends from band and grid (nb20) to the pick itself. Its retained bias grows with the band (−0.00 → −0.06/−0.08 → −0.09/−0.10 log-HR): the correction over-shoots as the region overshoots. The naive lower bound's one-sided coverage of 0.38–0.68 in the band settings is the sharpest price-of-selection number on record; its 0.94 under maxSG reflects β(Ĥ) ≈ 1.0 with a near-deterministic pick, not adequacy. The IJ two-term covers 0.971–0.999 at r 1.2–1.6, inflating to 2.8–3.3 under maxSG where the marginal SD collapses (pick nearly constant).

## 3. Standard summary table — Ĥᶜ block (per cell × setting × estimator)

Field bias on Ĥᶜ: the record states "−0.01 to −0.04 everywhere" (range; per-setting values not printed). 1s = one-sided upper (the exposed side). Marginal SD ratio carries target spread (Part A); error ratio beside it.

### HR 1.50 n500

| Setting | naive: bias log (mSD) | naive 1s cov | field 1s cov (r) | field 2s cov | IJ two-term 1s cov (r) | λ-SDᶜ/naive SE | SD(β̃ᶜ)/naive SE marg (err) |
|---|---|---|---|---|---|---|---|
| minSG | −0.056 (−0.40) | 0.867 (0.852, 0.882) | 0.923 (0.94) | 0.94 | 0.997 (1.82) | 0.99 | 1.06 (1.04) |
| ε 0.10 | −0.112 (−0.76) | 0.793 (0.775, 0.810) | 0.913 (0.89) | 0.94 | 0.996 (1.70) | 0.97 | 1.08 (1.01) |
| ε 0.20 | −0.123 (−0.79) | 0.781 (0.762, 0.798) | 0.897 (0.85) | 0.92 | 0.995 (1.64) | 0.95 | 1.10 (1.01) |
| ε 0.30 | −0.137 (−0.84) | 0.774 (0.755, 0.792) | 0.904 (0.85) | 0.92 | 0.996 (1.63) | 0.92 | 1.07 (1.00) |
| maxSG | −0.225 (−0.67) | 0.726 (0.706, 0.745) | 0.928 (0.97) | 0.94 | 0.995 (1.59) | 1.02 | 1.00 (1.02) |

### HR 1.75 n500

| Setting | naive: bias log (mSD) | naive 1s cov | field 1s cov (r) | field 2s cov | IJ two-term 1s cov (r) | λ-SDᶜ/naive SE | SD(β̃ᶜ)/naive SE marg (err) |
|---|---|---|---|---|---|---|---|
| minSG | −0.049 (−0.35) | 0.879 (0.864, 0.893) | 0.929 (0.94) | 0.94 | 0.997 (1.83) | 0.99 | 1.06 (1.03) |
| ε 0.10 | −0.098 (−0.65) | 0.819 (0.801, 0.835) | 0.919 (0.86) | 0.94 | 0.996 (1.65) | 0.97 | 1.12 (1.01) |
| ε 0.20 | −0.102 (−0.63) | 0.826 (0.809, 0.842) | 0.912 (0.82) | 0.92 | 0.995 (1.59) | 0.95 | 1.14 (1.02) |
| ε 0.30 | −0.113 (−0.69) | 0.807 (0.790, 0.824) | 0.903 (0.83) | 0.91 | 0.996 (1.60) | 0.93 | 1.09 (1.02) |
| maxSG | −0.191 (−0.49) | 0.776 (0.757, 0.794) | 0.925 (0.92) | 0.94 | 0.995 (1.60) | 0.99 | 1.04 (1.05) |

Readings: the field's complement upper coverage is **0.897–0.919 in the band settings and 0.923–0.929 at the two ends, never reaching 0.95**; it tracks λ-SDᶜ/naive SE (0.92–0.97 in the bands, 0.99–1.02 at the ends), which falls as the spread of |Ĥ|/|H| widens and returns to the naive SE where the pick stops varying (minSG top-3 mass 0.74–0.78; maxSG argmax 0.39–0.54). **Part A's mechanism corroborated from both ends of the dial**: the shortfall follows the pick's variability, not the complement's size or level. The error ratio is 1.00–1.05 at every setting — the marginal diagnostic's movement remains target spread. Residual 2–3 points below 0.95 at the ends, where the scale is right, is a separate smaller phenomenon (connects to the handoff §5 level-dimension open item, not to the scale mechanism). The IJ two-term covers 0.995–0.997 with upper bounds 1.16–1.37 (1.97–2.32 under maxSG) that rule out nothing.

## 4. Bound locations and the joint pair

| Cell | Setting | Ĥ field lower: mean (≥ 0.85 / ≥ 0.95) | Ĥᶜ field upper: mean (< 0.85 / < 0.80) | β(Ĥ) / β(Ĥᶜ) | Joint Bonferroni (calibrated) |
|---|---|---|---|---|---|
| HR 1.50 | minSG | 0.71 (0.21 / 0.12) | 1.08 (0.05 / 0.02) | 1.16 / 0.87 | 0.939 (0.938) |
| HR 1.50 | ε 0.10 | 0.71 (0.22 / 0.13) | 1.02 (0.13 / 0.06) | 1.25 / 0.84 | 0.940 (0.939) |
| HR 1.50 | ε 0.20 | 0.69 (0.20 / 0.12) | 0.99 (0.19 / 0.11) | 1.26 / 0.81 | 0.932 (0.932) |
| HR 1.50 | ε 0.30 | 0.68 (0.17 / 0.09) | 0.96 (0.25 / 0.16) | 1.23 / 0.78 | 0.935 (0.935) |
| HR 1.50 | maxSG | 0.80 (0.16 / 0.03) | 1.30 (0.14 / 0.10) | 0.98 / 0.76 | 0.937 (0.936) |
| HR 1.75 | minSG | 0.79 (0.32 / 0.21) | 1.13 (0.03 / 0.01) | 1.30 / 0.92 | 0.930 (0.929) |
| HR 1.75 | ε 0.10 | 0.84 (0.40 / 0.28) | 1.04 (0.11 / 0.06) | 1.51 / 0.85 | 0.940 (0.940) |
| HR 1.75 | ε 0.20 | 0.82 (0.37 / 0.26) | 1.00 (0.18 / 0.10) | 1.52 / 0.81 | 0.933 (0.933) |
| HR 1.75 | ε 0.30 | 0.81 (0.34 / 0.22) | 0.97 (0.25 / 0.15) | 1.45 / 0.77 | 0.935 (0.935) |
| HR 1.75 | maxSG | 0.83 (0.27 / 0.08) | 1.39 (0.14 / 0.10) | 1.01 / 0.76 | 0.942 (0.941) |

The joint Bonferroni pair holds at 0.930–0.942 at every setting (Wilson upper limits 0.940–0.952), calibrated = Bonferroni throughout (corr(Λ*, Λ*ᶜ) +0.01 to +0.05); γ off the Bonferroni floor on 44–50% of maxSG replicates (max 0.030) — informational (record's side issue 1). Along the band the two blocks' claim locations move in opposite directions: the harm-side location is strongest at ε 0.10 (share ≥ 0.85: 0.22 / 0.40) and weakens as β(Ĥ) dilutes; the complement-side benefit claim strengthens monotonically (upper mean 1.02/1.04 → 0.96/0.97, below 0.85 on a quarter of replicates at ε 0.30). Purity of the harm claim vs reach of the pair.

## 5. Checks performed

Gate 2 record ↔ Stage 3 report: identification block identical per cell (mean |Ĥ|, sens, spec, γ, p̂, band_n, complement fits, walls). Same draws / same K vs nb20 arm A asserted at Gate 2 per cell; 1,999 of 2,000 detected under every setting, the same undetected replicate. Nesting asserted on the full run (|Ĥ| ranges quoted). Internal spot checks across the report's four tables: consistent on every traced value. Pooled metas carry `effect_neighborhood` / `er_jcuts` (Part C live in production). Walls: 31–32 min (ε 0.30, as arm A — the band costs nothing), 23 min (size rules); projection 2.7 h, realized 2 h 36 m. Not in the record (so not in the tables above): absolute SDs/SEs, per-setting field bias on Ĥᶜ, IJ bias, naive r.

## 6. Decision items now with Larry (listed, not proposed)

1. **The band choice**, on the frontier above. The record's trade statement between ε 0.20 and 0.30: 14 sensitivity points and a doubling of the share reaching the planted size, against 5 specificity points, 2–4 PPV points, the first β(Ĥ) dilution, and 1–2 points of complement under-coverage on the HR 1.75 side; nothing on the constructions breaks at either.
2. **A fill-in ε (0.35 / 0.40)** to map the accelerating 0.30 → maxSG stretch — T-2's condition fired; not run, correctly, without instruction.
3. **Table conventions** (record's side issue 3): whether the standard tables adopt error-scale columns (SD(e), λ-SDᶜ/SD(e)) beside the marginal, given Part A's finding that marginal SD carries target spread on both blocks.
4. **The repair proposal go** is confirmed on the review side (A4 instability + banddial's two-end corroboration); drafting gated on A's formal close (committed record text) and the `R/` source of `.fs_mr_field_complement()`. Realistic target: close the band settings' 0.90–0.91 toward 0.92–0.93; the residual to 0.95 at concentrated-pick settings is the §5 level-dimension item.
5. Housekeeping: stale memory-monitor loop pid 2086198 (record's side issue 2; harmless); confirm the mid-run push to c133902e was Larry's, leaving cef19e9d to push.

## Status of the sibling reviews

Part A: decision made (repair warranted) and corroborated here; **formal close pending** the committed `REPORT_complement_variance_2026-09-07.md` (truncated in transcript; the A1 Var(e) = Var(a) + Var(c) − 2Cov(a, c) identity-check line not yet sighted). Part C: verified backward (re-combine identity exact, seven cells; save-disabled render wrote nothing) and forward (banddial pooled metas carry the knobs); **formal sign-off pending** the committed `REPORT_template_hygiene_2026-09-07.md`.
