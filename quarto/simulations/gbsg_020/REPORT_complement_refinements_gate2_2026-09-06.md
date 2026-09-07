# GATE 2 RECORD — Complement refinements, Stage 2 (campaigns `s7w` / `map1w`)

**Task:** `dev/tasks/TASK_complement_refinements_2026-09-06.md` (95e6c8c5). Records: Stage 0 291df669 · Stage 1 / Gate 1 PASS 8fbb3bdc.
**Run:** 2026-09-06, unattended under K-4, 100 workers, committed template at 8fbb3bdc driven by `FS_S7_*` env only (`FS_S7_FIELD_COMPLEMENT=TRUE`, `FS_S7_IJ_RESIDUAL=two_term`, `FS_S7_FB=none`), forestsearch 0.3.5 (8fbb3bdc installed), seeds `8316951 + sim_id`, two seed-disjoint batches of 1,000 then combine per cell, fail-fast per cell. Driver wall 15:32:46 → 18:42:44 for cells 1–7, plus the cell-7 re-run 18:43:29 → 19:18:33 (**3 h 46 m** total vs the 5 h ceiling; 7 h hard timeout untouched). No failed batches; no cell deferred or dropped; no mid-run changes to settings.

## GATE 2: ALL SEVEN CELLS PASS

Check script `gate2_checks_K.R` (session scratchpad; per-cell console records `gate2K_*.log`). Gated per cell: completeness; identity of **every existing column** (`fb_*`, wall-clock excluded) to the committed `s7c`/`map1c` bundle at ≤ 1e-12 relative; all 21 new columns finite on every detected replicate; interval invariants (`lo_w ≤ hi_w`, `lo_wf ≤ hi_wf`, both blocks; the calibrated pair never wider than Bonferroni's); `gamma ∈ [0.025, 0.05]`; achieved joint probability ≥ 0.95 − 2/n_joint (one draw of quantile discreteness per tail at the Bonferroni floor; the achieved value reported); `winner_floor ≥ naive` SE.

| # | Cell (campaign) | Wall | Detections | Existing columns identical | New columns finite | γ mean (share at 0.025) | joint prob mean (min; max shortfall, draws) | corr(Λ*, Λ*ᶜ) | Harm SE two-term / winner / floor (= naive) | Complement SE two-term / winner / floor (= naive) |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | h100 n500 (`s7w`) | 16 min | 1361 (68.0%) | **93/93, worst 0.0** | 21/21 | 0.0250 (0.968) | 0.9499 (0.9479; 1.90) | +0.047 | 0.399 / 0.175 / 0.324 | 0.255 / 0.123 / 0.139 |
| 2 | h175 n500 (`s7w`) | 20 min | 1900 (95.0%) | **93/93, worst 0.0** | 21/21 | 0.0251 (0.924) | 0.9502 (0.9480; 1.90) | +0.014 | 0.393 / 0.179 / 0.310 | 0.257 / 0.124 / 0.139 |
| 3 | h150 n500 (`map1w`) | 19 min | 1822 (91.1%) | **93/93, worst 0.0** | 21/21 | 0.0251 (0.943) | 0.9501 (0.9477; 1.90) | +0.030 | 0.395 / 0.177 / 0.315 | 0.256 / 0.123 / 0.139 |
| 4 | h150 n1500 (`map1w`) | 57 min | 1976 (98.8%) | **93/93, worst 0.0** | 21/21 | 0.0251 (0.901) | 0.9503 (0.9480; 1.90) | −0.049 | 0.246 / 0.110 / 0.183 | 0.151 / 0.073 / 0.080 |
| 5 | h075 n500 (`map1w`) | 14 min | 1042 (52.1%) | **93/93, worst 0.0** | 21/21 | 0.0250 (0.971) | 0.9498 (0.9477; 1.90) | +0.044 | 0.402 / 0.175 / 0.329 | 0.255 / 0.123 / 0.140 |
| 6 | h100 n1000 (`map1w`) | 26 min | 1319 (66.0%) | **93/93, worst 0.0** | 21/21 | 0.0250 (0.968) | 0.9499 (0.9479; 1.90) | +0.038 | 0.298 / 0.125 / 0.252 | 0.179 / 0.087 / 0.097 |
| 7 | h175 knoise3 (`map1w`) | 35 min | 1945 (97.2%) | **93/93, worst 0.0** | 21/21 | 0.0251 (0.898) | 0.9503 (0.9480; 1.85) | +0.038 | 0.395 / 0.176 / 0.319 | 0.256 / 0.124 / 0.138 |

Notes.

- Every existing column (the 64 s7-era columns plus the 8 uniform and the 30 complement-field columns of the `s7c`/`map1c` bundles — 93 gated after exclusions) is exactly equal on all 2,000 replicates in every cell; detection counts equal the committed bundles'. The additions changed nothing they were not meant to change.
- `gamma` sits at the Bonferroni floor on 90–97% of replicates and never exceeds 0.026: with the harm and complement draws nearly uncorrelated (mean corr −0.05 to +0.05; a candidate's harm and complement influences have disjoint supports), the independence solution 1 − √0.95 = 0.0253 rounds to the grid's 0.025, so the calibrated pair *is* the Bonferroni pair here. The achieved joint probability's sub-0.95 values (min 0.9477, at most 1.9 draws of ~1,000 under 0.95) are the quantile discreteness of two empirical 0.025 tails, not a defect.
- The winner-floor binds on essentially every replicate in both blocks: on the harm side the winner-only SE (0.11–0.18) sits well under the naive SE (0.18–0.33); on the complement side it sits just under (0.073–0.124 vs 0.080–0.140), so `winner_floor` = naive there too.
- Cost: `fit_mr_secs` equals the s7c/map1c values within noise (the additions are arithmetic on existing draws); cell walls match the complement run's to the minute.

## Incident: cell 7 first ran under the wrong campaign tag (corrected; no committed data lost)

The driver was derived from the complement task's by text substitution; the substitution `map1c → map1w` missed the last `CELLS` line (it ends with a closing quote), so cell 7 first rendered under `FS_S7_CAMPAIGN=map1c` and **overwrote the six committed `map1c` knoise3 files** (two batch bundles, the combined bundle, three HTMLs) in the working tree — a protocol breach (committed bundles are read-only). Detected at the driver's end from the campaign tag in the event line; the six files were **restored from git** (`git checkout --`, working tree verified clean of tracked modifications) and cell 7 was re-run under `map1w` with settings otherwise identical (35 min, 18:43 → 19:18). The mis-tagged bundles are kept in the session scratchpad (`mistagged_knoise3/`) as evidence only; nothing from them is used. The committed `map1c` files are byte-identical to their committed state (git-tracked; 0 modified).

Driver: `stage2K_driver.sh` (corrected in place after the incident) and `rerun_knoise3_map1w.sh`; per-render logs `s2K_*.log`; driver log `stage2K_driver.log`.
