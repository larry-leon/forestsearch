# GATE 2 RECORD — `effMaxSG` on the prevalence-30% cells, Stage 2 (campaign `p30sg`)

**Task:** `dev/tasks/TASK_p30_effMaxSG_2026-09-06.md` (70bb63a3). Records: Stage 0 631a9105 (Gate 0 PASS; alignment quoted) · Stage 1 / Gate 1 PASS 059eda40.
**Run:** 2026-09-06/07, unattended per Larry's pre-authorization, 100 workers, committed template at 059eda40 driven by `FS_S7_*` env only (`FS_S7_FOCUS=effMaxSG`, `FS_S7_Z1Q=0.60`, `FS_S7_FIELD_COMPLEMENT=TRUE`, `FS_S7_IJ_RESIDUAL=two_term`, `FS_S7_FB=none`), forestsearch 0.3.5 (8fbb3bdc installed), seeds `8316951 + sim_id`, two seed-disjoint batches of 1,000 then combine per cell, fail-fast per cell, `.refuse_if_tracked()` live on every save. Driver wall 22:46:26 → 00:53:16 (**2 h 07 m** for the four cells vs the 3 h ceiling; 4 h hard timeout untouched). No failed batches; no cell deferred or dropped; no mid-run changes.

## GATE 2: ALL FOUR CELLS PASS

Check script `gate2_checks_SG.R` (session scratchpad; per-cell records `gate2SG_c*.log`). Gated per cell: completeness (2,000 rows, sim_id 1–2,000, no duplicates, no CONFIG-ERROR, `meta$harm_z1_quantile = 0.60`, `meta$sg_focus = effMaxSG`); **same draws as `p30`** — `n_true` identical on all 2,000 rows and the `truth` object `identical()` to the p30 bundle's (the seed assertion); realized prevalence (distributional criterion, as p30); every MR / winner-only / field / complement-field / joint / β(Ĥ)/β(Ĥᶜ) quantity finite on every detected replicate; interval invariants; `gamma ∈ [0.025, 0.05]`; achieved joint probability ≥ 0.95 − 2/n_joint; bound↔quantile identities ≤ 1e-12; sens/spec/|Ĥ| recorded.

| # | Cell | Wall | Detections (= p30) | Same draws | Selection differs from p30 (both detected) | mean |Ĥ| p30sg / p30 (true) | sens p30sg / p30 | spec | Non-finite | γ mean (share at 0.025) | corr mean | complement fits (share new) | fit+MR s/rep |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | HR 1.00 n500 | 26 min | 1841 (92.0%) | yes | 959 / 1841 | 94.0 / 79.7 (153) | 0.301 / 0.251 | 0.862 | none | 0.0251 (0.947) | +0.073 | 530 (0.052) | 61.0 |
| 2 | HR 1.50 n500 | 30 min | 1999 (99.95%) | yes | 1166 / 1999 | 97.3 / 79.6 (153) | 0.450 / 0.357 | 0.917 | none | 0.0251 (0.903) | +0.046 | 484 (0.070) | 71.0 |
| 3 | HR 1.75 n500 | 31 min | 1999 (99.95%) | yes | 1126 / 1999 | 99.3 / 80.1 (153) | 0.506 / 0.395 | 0.936 | none | 0.0251 (0.893) | +0.028 | 457 (0.082) | 74.1 |
| 4 | HR 1.00 n1000 | 39 min | 1909 (95.5%) | yes | 1213 / 1909 | 170.3 / 132.0 (306) | 0.320 / 0.233 | 0.895 | none | 0.0251 (0.946) | +0.074 | 578 (0.065) | 100.1 |

Notes.

- Detection counts equal p30's exactly in every cell: `effMaxSG` and `maxeffCons` range over the same screened family (Stage 0), so a replicate detects under one iff it detects under the other; only the selected candidate differs (on 52–64% of detected replicates).
- `effMaxSG` never returns a *smaller* Ĥ than `maxeffCons` on the same replicate (Stage 3 table: share with |Ĥ_eff| < |Ĥ_cons| = 0.000 in every cell) — the neighbourhood rule can only move the pick to a larger in-band candidate or leave it.
- Complement fits per replicate rise to 457–578 (from 355–419 under `maxeffCons`) as anticipated; the complement block stays at 3.4–3.5 s/rep (7.3 s at n = 1000); cell walls 26–39 min vs p30's 21–36.
- Field draw usage healthy (min n_out 854–968); every draw winner's complement fit (`n_out_dropped_unfit = 0`); γ at the Bonferroni floor on 89–95% of replicates, never above 0.027; corr(Λ*, Λ*ᶜ) +0.03 to +0.07.
- Meta: every combined bundle records `campaign_tag = p30sg`, `sg_focus = effMaxSG`, `focus_tag = effMaxSG`, `harm_z1_quantile = 0.60`, `harm_prevalence_super = 0.3065`, `ci_method = field`, `field_complement = TRUE`, `ij_residual = two_term`, forestsearch 0.3.5; stems `fs_effMaxSG_…_z1q60_p30sg`.

Driver: `stage2SG_driver.sh` (session scratchpad); per-render logs `s2SG_*.log`; driver log `stage2SG_driver.log`.
