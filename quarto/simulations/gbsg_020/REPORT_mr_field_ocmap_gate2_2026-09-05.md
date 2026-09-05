# GATE 2 RECORD — OC map Stage 2 (campaign map1)

**Task:** `dev/tasks/TASK_mr_field_ocmap_2026-09-05.md` (0ef90878); Stage 0 26718386; Stage 1 1533b99d.
**Run:** 2026-09-05, unattended, 100 workers, committed template at b1316d22 via `FS_S7_*` env renders only (`FS_S7_UNIFORM` off everywhere), FB none, campaign `map1`. Driver wall 14:00:37 → 16:25:12 (**2 h 25 m** cumulative vs the 5 h ceiling; hard timeout 6 h untouched). Cell walls: c1 19.3 m · c2 14.0 m · c3 25.0 m · c4 33.3 m · c5 53.0 m. No failed batches; no cell deferred or dropped.

## GATE 2: ALL GREEN (per cell)

| Cell | Completeness | Detections | Field NAs (detected) | min n_out | Invariants | Anchored identity |
|---|---|---|---|---|---|---|
| c1 h150 n500 | 2000/2000, 0 CONFIG-ERROR | 1822 (91.1%) | 0 | 817 | PASS | **PASS**: 500 shared sims vs `fb_mr_m1_h15_knoise0_n500_combined_1_500` (0.2.2), 0 flips, worst 5.5e-15 |
| c2 h075 n500 | 2000/2000 | 1042 (52.1%) | 0 | 653 | PASS | no anchor (H-Q5) |
| c3 h100 n1000 | 2000/2000 | 1319 (66.0%) | 0 | 626 | PASS | **PASS**: 499/500 shared sims vs the four `fb_mr_m1_h10_knoise0_n1000_res_*` batch bundles, worst 3.6e-15; **1 cross-vintage flip enumerated** |
| c4 h175 knoise3 | 2000/2000 | 1945 (97.2%) | 0 | 874 | PASS | no anchor |
| c5 h150 n1500 | 2000/2000 | 1976 (98.8%) | 0 | 747 | PASS | no anchor |

**Enumerated flip (the 2026-09-05 convention):** c3 sim 309 — anchor (0.2.2) selected `{grade}` (n = 108), fresh (0.3.5) selects `!{er <= 0} & {grade}` (n = 103); oracle exact (rel diff 0), detection equal. Same knife-edge class as h10 sims 393/780 in the s7 arc; reported individually, excluded from the identity denominator.

Check script: session scratchpad `ocmap_gate2_checks.R`; driver `ocmap_stage2_driver.sh` with per-render logs `m1_*.log`.
