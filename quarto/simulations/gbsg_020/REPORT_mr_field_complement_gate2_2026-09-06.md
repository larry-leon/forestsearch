# GATE 2 RECORD — Complement field, Stage 2 (campaigns `s7c` / `map1c`)

**Task:** `dev/tasks/TASK_mr_field_complement_2026-09-06.md` (2c519337). Records: Stage 0 92d9ed20 · Stage 1 / Gate 1 PASS ef3e609a (5e540de2).
**Run:** 2026-09-06, unattended under H-C5, 100 workers, committed template at ef3e609a driven by `FS_S7_*` env only (`FS_S7_FIELD_COMPLEMENT=TRUE`, `FS_S7_FB=none`, `FS_S7_UNIFORM` unset), forestsearch 0.3.5 (ef3e609a installed), seeds `8316951 + sim_id`, two seed-disjoint batches of 1,000 then combine per cell, fail-fast per cell. Driver wall 12:04:44 → 15:13:56 (**3 h 09 m** for all seven cells vs the 6 h ceiling; 8 h hard timeout untouched). No failed batches; no cell deferred or dropped; no mid-run changes.

## GATE 2: ALL SEVEN CELLS PASS

Check script `gate2_checks.R` (session scratchpad; per-cell console records `gate2_*.log`). Gated per cell: completeness (2,000 rows, sim_id 1–2,000, no duplicates, no CONFIG-ERROR); identity of **every shared non-complement-field column** (`fb_*` and wall-clock columns excluded) to the committed `s7`/`map1` bundle at ≤ 1e-12 relative — the pairing proof; complement-field NA count (zero or documented by `fld_Hc_note`); interval invariants (`lo2s ≤ hi2s`, `lo1s ≤ up1s`, `lo_se ≤ hi_se`) and the bound↔quantile identities at ≤ 1e-12.

| # | Cell (campaign) | Wall | Detections | Shared columns identical | `sg_def` identical | Complement field finite / notes | n_out mean (min) | nfit mean (max) | share_newfit mean (max) | Complement s/rep (q90) | Harm field s/rep | fit+MR s/rep |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | h100 n500 (`s7c`) | 16 min | 1361 (68.0%) | **64/64, worst 0.0** | 2000/2000 | 1361/1361, 0 | 976 (697) | 414.7 (572) | 0.048 (0.579) | 2.90 (3.89) | 14.6 | 41.3 |
| 2 | h175 n500 (`s7c`) | 19 min | 1900 (95.0%) | **64/64, worst 0.0** | 2000/2000 | 1900/1900, 0 | 993 (766) | 382.0 (539) | 0.104 (0.663) | 3.01 (3.85) | 16.0 | 45.9 |
| 3 | h150 n500 (`map1c`) | 19 min | 1822 (91.1%) | **72/72, worst 0.0** | 2000/2000 | 1822/1822, 0 | 989 (817) | 394.6 (565) | 0.080 (0.803) | 3.04 (3.89) | 15.7 | 44.8 |
| 4 | h150 n1500 (`map1c`) | 57 min | 1976 (98.8%) | **72/72, worst 0.0** | 2000/2000 | 1976/1976, 0 | 994 (747) | 347.1 (572) | 0.202 (0.932) | 14.14 (19.56) | 35.9 | 161.6 |
| 5 | h075 n500 (`map1c`) | 14 min | 1042 (52.1%) | **72/72, worst 0.0** | 2000/2000 | 1042/1042, 0 | 966 (653) | 419.9 (562) | 0.042 (0.382) | 2.81 (3.79) | 13.9 | 39.3 |
| 6 | h100 n1000 (`map1c`) | 26 min | 1319 (66.0%) | **72/72, worst 0.0** | 2000/2000 | 1319/1319, 0 | 972 (626) | 416.8 (555) | 0.069 (0.692) | 6.73 (10.04) | 20.7 | 79.9 |
| 7 | h175 knoise3 (`map1c`) | 35 min | 1945 (97.2%) | **72/72, worst 0.0** | 2000/2000 | 1945/1945, 0 | 998 (874) | 603.2 (839) | 0.092 (0.915) | 6.45 (8.69) | 24.7 | 90.9 |

Notes.

- The `s7` bundles predate the uniform task and carry 64 shared columns; the `map1` bundles carry 72 (the eight `fld_H_kappa/…` columns, all-NA on both sides). Every shared column — detection, memberships, naive/oracle/MR (IJ) for both blocks, the harm field block, β(Ĥ)/β(Ĥᶜ) attachments — is exactly equal on all 2,000 replicates in every cell; `sg_def` strings are identical on 2,000/2,000 in every cell (same package vintage on both sides, so none of the s7-era cross-vintage flips appear).
- Detection counts equal the committed bundles' in every cell (the pairing is replicate-for-replicate).
- Interval invariants hold everywhere; bound identities ≤ 3.3e-16 in every cell. `n_out` minima equal the committed harm field's (the complement drops no outer draw the harm field kept — `n_out_dropped_unfit` is zero throughout, every draw winner's complement having been fit).
- Cost: the complement block is 2.8–3.0 s per detected replicate at n = 500 (≈ 7% of fit+MR), 6.5–6.7 s at n = 1,000 / knoise3, 14.1 s at n = 1,500 (≈ 9%). Cell walls are within 0–4 min of the committed runs'. Complement fits per replicate average 350–600 (max 839), of which the multiplier stage had already fit most: the share of outer + inner draw-winner readings needing a new fit averages 4–20% (n = 1,500 highest).
- Meta: every combined bundle records `campaign_tag` `s7c`/`map1c`, `ci_method = field`, `field_complement = TRUE`, forestsearch 0.3.5, `fb_mode` none/none.

Driver: `stage2_driver.sh` (session scratchpad); per-render logs `s2_*.log`; driver log `stage2_driver.log` (timestamps above).
