# REPORT — OC map: Stage 1 (Smoke and projection)

**Task:** `dev/tasks/TASK_mr_field_ocmap_2026-09-05.md` (0ef90878); Stage 0 record 26718386.
**Date:** 2026-09-05, unattended. Campaign `map1smoke`, 5 replicates per cell at the committed seeds, 5 workers, `FS_S7_FB=none`, `FS_S7_UNIFORM` unset (FALSE).

## GATE 1 (pre-authorized): ALL FIVE CELLS PASS 1a AND FIT THE CEILING — all run in Stage 2

### 1a — Smoke (all PASS)

| Cell | Render | Detected (of 5) | Anchored identity (31 cols, sims 1–5) | Field | fit+MR+field s/rep (light load) |
|---|---|---|---|---|---|
| 1 h150 n500 | OK, no CONFIG-ERROR | 5 | **PASS, worst 6.6e-16** (vs `fb_mr_m1_h15_knoise0_n500_combined_1_500`) | finite, invariants hold | 20–23 (field 8–9) |
| 2 h075 n500 | OK | 3 | no anchor (H-Q5) | finite, invariants hold | 19–21 (field 8) |
| 3 h100 n1000 | OK | 5 | **PASS, worst 8.9e-16** (vs `fb_mr_m1_h10_knoise0_n1000_res_1_20`) | finite, invariants hold | 28–30 (field 9) |
| 4 h175 knoise3 | OK | 5 | no anchor | finite, invariants hold | 32–38 (field 11–12) |
| 5 h150 n1500 | OK | 5 | no anchor | finite, invariants hold | 35–37 (field 9–10) |

Uniform columns are all-NA on every render (the uniform task's knob stays off, its state untouched). No cross-vintage flips appeared in the anchored windows (0.2.2 anchors, ulp-level agreement only).

### 1b — Projection (100 workers, 2,000 replicates per cell)

Loaded per-replicate cost estimated as smoke light-load cost × the s7-measured load-inflation factor (×1.4 at h100-like detection, ×1.9 at h175-like; the 0.2.x-vintage sweep timings overstate the current package by the 0.3.4/0.3.5 search speedup and are not used):

| Cell | Loaded s/rep (est.) | Wall est. (2 batches + combine) | Cumulative |
|---|---|---|---|
| 1 h150 n500 | ~37 | ~25 min | 0.4 h |
| 2 h075 n500 | ~19 (detection ~60%, undetected reps cheap) | ~20 min | 0.8 h |
| 3 h100 n1000 | ~49 | ~35 min | 1.4 h |
| 4 h175 knoise3 | ~66 | ~45 min | 2.1 h |
| 5 h150 n1500 | ~68 | ~47 min | **2.9 h** |

Cumulative ≈ **2.9–3.5 h** with render overheads — inside the 5 h ceiling (H-Q6), hard timeout 6 h. **No cell deferred, none dropped.** Stage 2 order: 1 → 2 → 3 → 4 → 5, two seed-disjoint batches each then combine, fail-fast per cell (a failed batch stops that cell, the next proceeds).
