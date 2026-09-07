# REPORT — `effMaxSG` band 0.20 at 31% prevalence: Stage 1 (Knobs, identities, smoke, projection)

**Task:** `dev/tasks/TASK_p30sg_nb20_2026-09-07.md` (095d8c4d); Stage 0 record 24911dc8 (Gate 0 PASS). Decisions S-1–S-4 as stated; winner-only and winner-floor excluded from every table, figure and report line.
**Date:** 2026-09-07. Executor: Claude Code, unattended per Larry's pre-authorization. No `R/` changes (suite run for the record).

---

## GATE 1: PASS — knob-inert identity exact against `p30sg` (114/114 shared non-p̂ columns, worst relative difference 0.0, `truth` identical, on three cells); arm A smoke: the ε = 0.20 band admits at least as many candidates as ε = 0.10 on every replicate and |Ĥ| is never smaller than p30sg's on any replicate; arm B smoke: family size K reported (1,920–2,290), prevalence unchanged (`n_true` identical to p30sg on every shared replicate), fields finite, γ in range, bound identities ≤ 1.7e-16 in every cell; the observed Ĥ reproduced by the gate's re-selection map on the unperturbed effects on 3/3 alignment fits (ε = 0.20; J = 10 and J = 20); suite 0 fail / 5025 pass / 3 skip / 32 warn. Projection ≈ 1.1 h (arm A) + 4.1 h (arm B) ≈ 5.2 h at 100 workers, inside the 7 h ceiling (9 h hard timeout). **Compute go per S-5: arm A both cells, arm B all five cells; none deferred at the start** (the driver re-projects arm B from realized walls before each cell and defers any cell whose projected finish would cross the ceiling).

## 1a — Template change (document-level, add-only; defaults reproduce `p30sg` exactly)

- `effect_neighborhood <- .env_num("FS_S7_NBHD", 0.10)` (validated in [0, 1)) and `er_jcuts <- .env_int("FS_S7_ER_JCUTS", 10L)` read next to the focus knob, before the stem; the later literal pin becomes a pointer comment; `fs_conf.cont_jcuts <- list(er = er_jcuts)`. `forestsearch()` already receives `effect_neighborhood` and forwards it to the gate as `nbhd` (Stage 0 0a).
- Stem tags `_nb%02d` and `_j%d` only when non-default (empty at defaults, every committed stem unchanged). Campaign tags written alphanumerically per the template's guard: `p30sgnb20`, `p30sgnb20j20` (Stage 0, Naming).
- `return_reselection = TRUE` in `mr_inference_args`; new recorder columns `n_family` (gate's kept family K), `n_cons_qual` (consistency-qualifying candidates, from `grp.consistency$out_sg$result`), `band_n` (candidates with hr ≥ (1 − ε)·max hr on the observed effects), `p_hat_H`, `p_hat_sum`, `p_hat_top1..3`, `p_hat_top_labels`; a render-side "Re-selection regime" callout after the subgroup-size callout (renders nothing on older bundles).
- Knobs echo line gains `nbhd=`, `er_jcuts=`; both join the bundle meta and the combine poolability gate.

## 1b — Identities and smoke (all PASS; `smoke_nb20_check.R`, `align_nb20.R`, session scratchpad; 5 workers, sequential renders, sims 1–5 at the p30sg seeds)

**Knob-inert identity** (all knobs at default, `FS_S7_FOCUS=effMaxSG`, `FS_S7_Z1Q=0.60`): campaign `nbinert` (HR 1.00 n500) and, for the band counts below, campaign `nb10smoke` (HR 1.50 n500, HR 1.75 n500) — on all three cells **114/114** shared non-p̂ columns identical to the committed `p30sg` bundle (worst relative difference 0.0; timing and `fb_*` excluded; the nine new columns have no p30sg counterpart), `truth` `identical()`, `meta$sg_focus = effMaxSG`, `effect_neighborhood = 0.10`, `er_jcuts = 10`. p̂ finite on every detected replicate; top-3 mass 0.16–0.65 (≤ 1); Σp̂ 0.998–1.000.

| Cell (defaults) | p̂(Ĥ) per replicate | top-3 mass | K | consistency-qualifying | band (ε 0.10) |
|---|---|---|---|---|---|
| HR 1.00 n500 (`nbinert`) | 0.154, 0.401, —, 0.031, 0.159 | 0.359, 0.568, —, 0.162, 0.407 | 1233, 1297, —, 1280, 1303 | 7, 28, —, 13, 23 | 4, 1, —, 3, 6 |
| HR 1.50 n500 (`nb10smoke`) | 0.146, 0.307, 0.093, 0.059, 0.194 | 0.365, 0.471, 0.251, 0.299, 0.595 | 1233, 1297, 1212, 1280, 1303 | 118, 227, 48, 337, 89 | 4, 1, 2, 6, 2 |
| HR 1.75 n500 (`nb10smoke`) | 0.142, 0.227, 0.074, 0.254, 0.215 | 0.341, 0.466, 0.257, 0.466, 0.654 | 1233, 1297, 1212, 1280, 1303 | 244, 388, 128, 508, 169 | 2, 1, 3, 1, 2 |

**Arm A smoke** (campaign `nb20smoke`: ε = 0.20, J = 10), paired to p30sg by sim_id (same data: `n_true` identical):

| Cell | band ε 0.20 vs ε 0.10 (per rep) | ≥ every rep | |Ĥ| ε 0.20 vs p30sg | ≥ every rep | fields finite | γ | bound ids | p̂(Ĥ) | top-3 mass | fit+MR s/rep (light) | complement fits |
|---|---|---|---|---|---|---|---|---|---|---|---|
| HR 1.50 n500 | 8, 1, 14, 28, 3 vs 4, 1, 2, 6, 2 | yes | 161, 69, 113, 180, 87 vs 84, 69, 72, 153, 68 | yes | yes | in range | 5.6e-17 | 0.013, 0.186, 0.064, 0.037, 0.159 | 0.33, 0.34, 0.23, 0.28, 0.47 | 36.3 | 637 |
| HR 1.75 n500 | 15, 3, 15, 3, 3 vs 2, 1, 3, 1, 2 | yes | 173, 127, 167, 153, 87 vs 84, 69, 69, 68, 68 | yes | yes | in range | 1.2e-16 | 0.040, 0.086, 0.041, 0.201, 0.180 | 0.26, 0.35, 0.25, 0.41, 0.50 | 37.7 | 600 |

Early read, not a check: the wider band moves the pick to the largest of 3–28 in-band candidates, and |Ĥ| now reaches 153–180 (the planted size at n = 500 is ~153) on 5 of 10 replicates; p̂(Ĥ) falls to 0.01–0.20 — the selected candidate is rarely the re-selection argmax.

**Arm B smoke** (campaign `nb20j20smoke`: ε = 0.20, J = 20; five cells):

| Cell | detected | K (family) | consistency-qualifying | band | |Ĥ| (p30sg J = 10 beside) | same draws | prevalence mean | fields finite | γ | bound ids | p̂(Ĥ) | fit+MR s/rep (light) | complement fits |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HR 1.00 n500 | 4/5 | 1943, 2018, —, 2125, 2155 | 8, 40, —, 13, 30 | 5, 1, —, 9, 7 | 84, 73, —, 193, 102 (77, 69, —, 134, 102) | yes | 0.326 | yes | in range | 1.7e-16 | 0.085, 0.153, 0.005, 0.221 | 35.7 | 969 |
| HR 1.50 n500 | 5/5 | 1943, 2018, 1920, 2125, 2155 | 153, 323, 64, 515, 153 | 8, 1, 18, 38, 4 | 159, 73, 166, 174, 87 (84, 69, 72, 153, 68) | yes | 0.326 | yes | in range | 1.1e-16 | 0.019, 0.107, 0.016, 0.012, 0.161 | 49.9 | 875 |
| HR 1.75 n500 | 5/5 | 1943, 2018, 1920, 2125, 2155 | 354, 560, 176, 804, 284 | 15, 10, 28, 6, 4 | 182, 160, 176, 157, 112 (84, 69, 69, 68, 68) | yes | 0.326 | yes | in range | 1.2e-16 | 0.012, 0.038, 0.005, 0.054, 0.088 | 53.5 | 821 |
| HR 1.50 n1000 (new) | 5/5 | 2155, 2290, 2167, 2159, 2051 | 463, 412, 740, 517, 793 | 22, 21, 13, 29, 24 | 228, 246, 265, 328, 262 (no reference) | n/a | 0.312 | yes | in range | 1.4e-16 | 0.063, 0.036, 0.070, 0.026, 0.030 | 70.1 | 833 |
| HR 1.00 n1000 | 5/5 | 2155, 2290, 2167, 2159, 2051 | 57, 75, 95, 60, 264 | 16, 13, 3, 30, 25 | 160, 229, 115, 263, 262 (144, 154, 115, 199, 175) | yes | 0.312 | yes | in range | 1.1e-16 | 0.112, 0.028, 0.199, 0.030, 0.013 | 61.9 | 1024 |

The J = 20 family is ×1.6–1.7 the J = 10 family (K 1,920–2,290 vs 1,212–1,303; Stage 0's probe), the consistency-qualifying set ×1.3–1.7, and the same replicates detect (sim 3 of the null n500 cell is undetected under every arm, as in p30sg). Prevalence is the design's on the shared cells (`n_true` identical) and 0.312 on the five new-cell replicates against 0.3065.

**Alignment in numbers** (standalone fits of HR 1.50 n500 under `effMaxSG`, ε = 0.20, `return_reselection = TRUE`; the candidate table from `fs$grp.consistency$out_sg$result`):

| sim | J | consistency-qualifying | in band (hr ≥ 0.8·max) | largest in-band N (ties) | identifier top row = gate map | observed |Ĥ| = top-row N | Σp̂ = selection rate | p̂(Ĥ) / argmax p̂ |
|---|---|---|---|---|---|---|---|---|
| 1 | 10 | 247 | 20 | 159 (1) | yes | yes | 1.000 = 1.000 | 0.018 / 0.132 |
| 1 | 20 | 374 | 2 | 98 (1) | yes | yes | 1.000 = 1.000 | 0.122 / 0.203 |
| 4 | 20 | 66 | 8 | 93 (1) | yes | yes | 1.000 = 1.000 | 0.021 / 0.101 |

On every fit the identifier's `order(−in_band, −N, −Pcons, −hr, K)` top row is the gate map's `which.max(N[in_band])` pick with no exact-N ties in the band (the Stage 0 residual never engaged), and the observed Ĥ has that N. Caveat as at p30sg Stage 1: standalone fits reproduce the harness up to the ambient RNG, so their picks need not coincide with the smoke bundles' on the same sim ids; the alignment is a self-consistency of each fit.

**Suite:** 0 fail / 5025 pass / 3 skip / 32 warn (no `R/` change; run for the record).

## 1c — Projection (100 workers)

Light-load fit+MR per replicate against the p30sg smoke's (26.1 / 35.8 / 37.8 / 41.3 s for HR 1.00 n500 / 1.50 n500 / 1.75 n500 / 1.00 n1000): arm A 36.3 / 37.7 s (×1.0 — the band costs nothing; complement fits +30%); arm B 35.7 / 49.9 / 53.5 s at n = 500 (×1.37–1.42) and 61.9 s for the null n1000 cell (×1.50), 70.1 s for the harm 1.5 n1000 cell. Scaling the p30sg walls (26 / 30 / 31 / 39 min): **arm A ≈ 30 + 31 ≈ 1.0–1.1 h**; **arm B ≈ 36 + 42 + 44 + 66 + 59 ≈ 4.1 h**; total ≈ 5.2 h, inside the 7 h ceiling (9 h hard timeout). Order: arm A (HR 1.50 n500 → 1.75 n500), then arm B (HR 1.00 n500 → 1.50 n500 → 1.75 n500 → 1.50 n1000 → 1.00 n1000); two seed-disjoint batches then combine per cell; fail-fast per cell; Gate 2 per cell after each combine; the driver (`stage2_nb20_driver.sh`, session scratchpad) re-projects each arm B cell from the realized walls (arm A mean × 1.4 until the first arm B n500 cell lands, then the realized arm B n500 mean; ×1.6 for n = 1000) and defers, listing it, any cell whose projected finish would cross the ceiling; every render caps its own timeout at the 9 h hard stop; every render is the committed template driven by `FS_S7_*` env only; the save guard is live.
