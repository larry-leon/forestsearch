# REPORT — `effMaxSG` on the prevalence-30% cells: Stage 1 (Knob, identities, smoke, projection)

**Task:** `dev/tasks/TASK_p30_effMaxSG_2026-09-06.md` (70bb63a3); Stage 0 record 631a9105 (Gate 0 PASS; alignment quoted). Q1–Q4 at defaults (winner-only evaluated; winner-floor excluded).
**Date:** 2026-09-06. Executor: Claude Code, unattended per Larry's pre-authorization. No `R/` changes.

---

## GATE 1: PASS — knob-inert identity exact against `p30` (114/114 columns, truth identical), the `effMaxSG` smokes pass in every cell (selection differs from `maxeffCons` on 2–3 of 5 replicates; fields finite; γ in range; bound identities 2.2e-16), the observed Ĥ is reproduced by the gate's re-selection map on the unperturbed effects on 5/5 fitted replicates with no in-band size ties, suite green (0 fail / 5025 pass / 3 skip / 32 warn), projection ≈ 2 h at 100 workers (< 3 h ceiling, 4 h hard timeout). Compute go: all four cells run, none deferred.

## 1a — Template change (document-level, add-only)

`sg_focus <- .env_chr("FS_S7_FOCUS", "maxeffCons")` with `stopifnot(sg_focus %in% c("maxeffCons", "effMaxSG"))`, replacing the literal at the same place (before the stem, which already follows `focus_tag`); the knobs echo line gains `focus=`; `sg_focus` joins the combine poolability gate. Nothing else changes — `selection_rule = "neighborhood"`, `effect_neighborhood = 0.10`, the thresholds and every MR knob were already explicit and forwarded, and `forestsearch()` derives the gate's re-selection rule from `sg_focus`. Stems: `fs_effMaxSG_fb_mr_field_m1_…_z1q60_p30sg` (the focus tag names the identifier, as it always has).

## 1b — Identities and smoke (all PASS; `p30sg_smoke_check.R`, `p30sg_align3/4.R`, session scratchpad)

**Knob-inert identity** (campaign `sginert`: HR 1.00 n500, `FS_S7_Z1Q=0.60`, sims 1–5 at the p30 seeds, default focus): all **114** shared columns identical to the committed `p30` bundle (worst relative difference 0.0; `fb_*`/wall-clock excluded), `truth` `identical()`, `meta$sg_focus = maxeffCons`.

**`effMaxSG` smoke** (campaign `p30sgsmoke`, `FS_S7_FOCUS=effMaxSG`, `FS_S7_Z1Q=0.60`, 5 replicates per cell at the p30 seeds, 5 workers, sequential renders), paired to `p30` by sim_id (same data: `n_true` identical on every replicate):

| Cell | selection differs from maxeffCons | detected | |Ĥ| effMaxSG vs maxeffCons (per replicate) | sensitivity effMaxSG vs maxeffCons | fields finite | γ | bound ids | complement fits | fit+MR s/rep (light) |
|---|---|---|---|---|---|---|---|---|---|
| HR 1.00 n500 | 3/5 | 4/5 | 77,69,—,134,102 vs 63,69,—,64,68 | 0.14,0.00,—,0.37,0.27 vs 0.12,0.00,—,0.23,0.36 | yes | in range | 2.2e-16 | 556 | 26.1 |
| HR 1.50 n500 | 3/5 | 5/5 | 84,69,72,**153**,68 vs 84,69,64,68,66 | 0.39,0.00,0.16,**0.89**,0.36 vs 0.39,0.00,0.23,0.27,0.43 | yes | in range | 2.2e-16 | 500 | 35.8 |
| HR 1.75 n500 | 2/5 | 5/5 | 84,69,69,68,68 vs 84,69,64,68,66 | 0.39,0.00,0.43,0.27,0.36 vs 0.39,0.00,0.23,0.27,0.43 | yes | in range | 2.2e-16 | 479 | 37.8 |
| HR 1.00 n1000 | 3/5 | 5/5 | 144,154,115,**199**,175 vs 107,154,115,109,141 | 0.36,0.32,0.12,0.27,0.30 vs 0.18,0.32,0.12,0.32,0.19 | yes | in range | 2.2e-16 | 567 | 41.3 |

The knob reaches the identifier (selections differ), `meta$sg_focus = effMaxSG` in every bundle, and the winner-only columns (`mr_*_se_w`) are finite on every detected replicate. Early read, not a check: `effMaxSG` returns a larger Ĥ on some replicates (up to the full 153 at n = 500 and 199 at n = 1000) and the same or a same-sized one on others; sensitivity moves in both directions at n = 5 — Stage 3 will say. Complement fits per replicate rise to 480–570 (from 355–420 under `maxeffCons`), as anticipated at Gate 0.

**Alignment in numbers** (task 1b; standalone fits of the HR 1.00 n500 cell under `effMaxSG` with `return_reselection = TRUE`, sims 1, 2, 3, 5, 6; the candidate table read from `fs$grp.consistency$out_sg$result`):

| sim | family (consistency-qualifying) | in band (hr ≥ 0.9·max) | largest in-band N (ties) | identifier top row = gate map | observed |Ĥ| = top-row N | p̂ sum = selection rate | p̂(Ĥ) / argmax p̂ |
|---|---|---|---|---|---|---|---|
| 1 | 71 | 6 | 123 (1) | yes | 123 | 0.999 = 0.999 | 0.038 / 0.131 |
| 2 | 21 | 2 | 91 (1) | yes | 91 | 0.996 = 0.996 | 0.149 / 0.149 |
| 3 | 57 | 2 | 76 (1) | yes | 76 | 1.000 = 1.000 | 0.079 / 0.222 |
| 5 | 12 | 1 | 64 (1) | yes | 64 | 0.994 = 0.994 | 0.077 / 0.077 |
| 6 | 23 | 3 | 91 (1) | yes | 91 | 1.000 = 1.000 | 0.221 / 0.221 |

On every fitted replicate the identifier's `order(−in_band, −N, −Pcons, −hr, K)` top row is the gate map's `which.max(N[in_band])` pick, with **no exact-N ties among in-band candidates** (so the residual tie-break difference of Stage 0 never engaged), and the observed Ĥ has exactly that N; the re-selection frequencies sum to the selection rate. **Caveat recorded:** these standalone fits reproduce the template's harness only up to the ambient RNG, so their selections do not coincide with the smoke bundles' on the same sim ids (e.g. sim 1: 123 here vs 77 in the smoke); the alignment property is a self-consistency of each fit and does not depend on which fit is examined. Regime read from p̂: under `effMaxSG` at 30% the selection is **competitive** — p̂(Ĥ) 0.04–0.22, the argmax candidate often another (0.13–0.22) — the tie regime the task's cell 1 was designed to expose.

**Suite:** 0 fail / 5025 pass / 3 skip / 32 warn (no `R/` change; run for the record).

## 1c — Projection (100 workers)

Light-load fit+MR 26–41 s/rep (vs 20–34 in the p30 smoke: the larger complement-fit count adds a few seconds). Anchors: p30 walls 21 / 23 / 24 / 36 min → **≈ 1.8–2.0 h for the four cells**, inside the 3 h ceiling (4 h hard timeout). Order HR 1.00 n500 → 1.50 n500 → 1.75 n500 → 1.00 n1000; campaign `p30sg`; two seed-disjoint batches then combine; per-cell fail-fast; the save guard live (`stage2SG_driver.sh`, session scratchpad; every render is the committed template driven by `FS_S7_*` env only).
