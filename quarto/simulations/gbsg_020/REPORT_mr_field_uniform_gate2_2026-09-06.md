# GATE 2 RECORD — Field uniform calibration Stage 2 (campaigns s7u + uniform)

**Task:** `dev/tasks/TASK_mr_field_uniform_2026-09-05.md` (1313f597), scope per the Gate 1 adjudication (5577b079): H3 M_cap 40 / mass target 0.99, κ grid [1, 3], H4 unchanged (300/300/150), H5 = h100/h175 sims 1–500 (campaign `s7u`) + the four Guo–He check cells (campaign `uniform`). Implementation b1316d22 + 3ad6ee13; 2b driver 8975af1c.
**H6 go decision:** measured at the raised cap (5 reps/cell, campaign `h6smoke`, committed evidence bundles): κ sweep 149–166 s/rep lightly loaded → projection ≈ 4.1 h ≤ 6 h → **compute pre-authorized, H4 unchanged**. **Actual wall 7 h 15 m** (17:27:59 → 00:43:19, 2026-09-05/06), inside the 8 h hard timeout but above the 6 h projection: the s7u batches ran ~×1.6 over projection (κ at load 583/651 s vs the ~300–490 s estimate) and t7's nested family at cap 40 took 192.7 min alone. Recorded as a projection miss, not a gate breach (the go condition was on the projection; no mid-run changes were made).

## GATE 2: ALL GREEN

**2a (s7u, sims 1–500, `field_uniform = TRUE`, FB none):**

| Cell | Completeness | Detections | Non-uniform identity vs s7 | κ* finite | κ* (mean/SD/range) | M / mass | minC₁ (mean/min) | κ secs |
|---|---|---|---|---|---|---|---|---|
| h100 | 500/500, 0 CONFIG-ERROR | 329 | **64/64 columns identical** (FB-join + wall excluded) | all 329 | 1.62 / 0.14 / [1.15, 2.12], 0 ceiling hits | mean 39.8, cap-40 on 323/329; mass mean 0.927, min 0.779, ≥ 0.99 on 8 | 0.925 / 0.880 | mean 583 s |
| h175 | 500/500 | 477 | **64/64 identical** | all 477 | 1.56 / 0.13 / [1.07, 2.06], 0 ceiling hits | mean 39.0, cap-40 on 431/477; mass mean 0.949, min 0.809, ≥ 0.99 on 51 | 0.925 / 0.892 | mean 651 s |

mcse(κ*) mean 0.096 / 0.093 at the H4 counts (recorded per replicate). **H3 note:** even at cap 40 the 0.99 mass target binds against the cap on most forestsearch replicates (mass means 0.93–0.95); recorded per replicate as instructed.

**2b (Guo–He check set, campaign `uniform`, driver `mr_field_uniform_vs_guohe_run.R`):** every cell 2000/2000 replicates, 0 errored, **0 pairing mismatches** (18 pre-existing columns — MR (IJ) est/SE/bias terms and the full field block — `identical()` per replicate against the committed 2026-09-05 field bundles), 0 uniform NAs.

| Cell | Elapsed | κ* (mean/SD/range) | M range | mass ≥ 0.99 | minC₁ min | κ secs |
|---|---|---|---|---|---|---|
| t35_beta2_03 | 26.1 min | 1.26 / 0.08 / [1.03, 1.54] | 1–2 | 2000/2000 | 0.887 | 76 |
| t35_beta2_05 | 25.7 min | 1.23 / 0.10 / [1.03, 1.51] | 1–2 | 2000/2000 | 0.887 | 74 |
| t6_k10 | 81.6 min | 1.56 / 0.09 / [1.21, 1.92] | 3–10 | 2000/2000 | 0.883 | 237 |
| t7_beta2_00 | 192.7 min | 1.15 / 0.06 / [1.01, 1.39] | 12–40 (cap on 1872/2000) | 141/2000 | 0.880 | 560 |

Check script: session scratchpad `uniform_gate2_checks.R`; driver logs `s7u_*.log`.
