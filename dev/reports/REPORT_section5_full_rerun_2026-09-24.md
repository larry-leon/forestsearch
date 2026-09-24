# REPORT — Section 5 full re-run under the aligned package (2026-09-24)

Task: `dev/tasks/TASK_section5_full_rerun_2026-09-24.md` (committed `24c18738`). Payloads and logs: `5fec9198`.

**Run status: COMPLETE.** 18 of 18 cells, 2,000 replicates each, 0 errors, no overrun, no memory stop.
**Answer to §4:** declaration and selected region are identical on every replicate of every cell, but the
re-run **does** move three published range endpoints by one replicate each (0.0005) at four decimals:
field-s upper 0.9125 → 0.9130, `joint` 0.9320 → 0.9315, `joint_s` 0.9395 → 0.9400 (§4.3).

## 0. Provenance

- HEAD at start `fdd2f500`; task doc `24c18738`. R 4.5.2, aarch64-apple-darwin20, Mac-Studio-3, 14 physical cores, 36 GB.
- **0a.** `7713942e`, `96f84ad8`, `06ac5391`: all ancestors of HEAD.
- **0b.** `R_LIBS` unset; `R/`, `DESCRIPTION`, `NAMESPACE` clean; `devtools::install(dependencies = FALSE, upgrade = FALSE)`.
  Verified from inside a `multisession` worker (pid ≠ main): `pconsistency.digits` in `formals(fs_mr_inference)` TRUE,
  `.fs_pcons_eff` in its body TRUE. Library `/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/forestsearch`,
  version 0.3.5.9000, built `2026-09-24 04:00:02 UTC`. The same check ran from all 13 workers before every render
  (`s5rerun/logs/*.buildcheck`, all `OK=TRUE`).
- Pre-existing and left out of every commit: 9 staged files under `quarto/simulations/actg175/binary_020/Claude outputs/`,
  and the modified `quarto/simulations/gbsg_020/REPORT_null_gbsg_identification_2026-09-21.md`.
- No `R/` file modified.

## 1. Configuration

- Cells, parameters and source stems from `mrs5sweep/cells.txt`; seeds `8316951 + sim_id`, `sim_id` 1–2000 (`seed_base` 8316951 in every meta).
- Knob set: exactly `scripts_p12x20/campaign_p12x20.sh`'s (= `mrs5sweep/runcell.sh`'s): `FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20
  FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none
  FS_S7_RETURN_RESEL=TRUE`, `FS_S7_Z1Q=0.60` on the 31% cells, every other inherited `FS_S7_*` unset. Unchanged template and `scripts_dinamr/render.sh`.
- Structure: per cell `batch 1 (1000)` → `batch 1001 (1000)` → `combine`, as `campaign_p12x20.sh`. Bundles are
  `results/fs_effMaxSG_fb_mr_field_m1_h{HHH}_knoise0_n{N}[_z1q60]_nb20_s5rerun_{res_1_1000,res_1001_2000,combined_1_2000}.rds`,
  the committed name pattern with the tag changed, so a figure re-run is a tag change (the `p12x20` bundles sit in
  `p12x20_2026-09-12/`, these in `results/`, as `cert20`/`e1stud`). Renders `s5rerun_<cell>_<mode>_<start>.html` in the `gbsg_020/` root.
- Tags: `s5rerun` (the run), `s5rerunchk` (Step 2 check), `s5smk1500` (memory smoke). No committed payload modified or overwritten.
- **Deviations** (runner `s5rerun/runcell.sh`, driver `s5rerun/driver.sh`):
  1. Mac paths; the installed library instead of a scratch `R_LIBS`; `render.sh` invoked with `zsh` (its shebang).
  2. Workers 13 on every cell (the template caps at physical cores − 1). Committed: 64 (`p12x20`), 100 (`cert20`, `e1stud`).
  3. Empty-array expansions guarded for macOS bash 3.2 (first launch of the Step 2 check failed on this before rendering).
  4. `s5rerun/record_cell.R` failed on A7 (`meta$n_workers` is NULL on combined metas); made NA-tolerant and A7's record
     backfilled by hand (noted in `driver.log`). Worker count per cell is from the driver.
  5. `caffeinate -is -w <driver pid>` attached at 21:35 on 2026-09-23 (the driver was launched under `caffeinate -i`); noted in `driver.log`.

## 2. Step 2 check — A7, 200 replicates, against `mrs5probe` A7

| quantity | result |
|---|---|
| declaration | identical on all 200 (126 / 200, 0.6300 both) |
| selected region (label, `n_sel`) | 0 changed |
| MR top label | 0 changed |
| errors / `mr_ok` false | 0 / 0 |
| `n_true` | identical on all rows; truth agrees within 1e-8 |
| max \|Δ\| corrected estimate `mr_H_est` | 4.4e-16 |
| max \|Δ\| field lower `fld_H_lo1s` | 6.1e-16 |
| max \|Δ\| field-s upper `fld_Hc_up1s_s` | 4.4e-16 |
| max \|Δ\| Bonferroni `joint_s` lower / upper | 1.2e-15 / 4.4e-16 |
| max \|Δ\| any other non-timing numeric column | 4.3e-15 (`itt_est`) |

The differences are machine precision. Wall 192 s at 13 workers (loop 166 s). Log `s5rerun/logs/s5rerunchk_A7_*`, pair `s5rerun/pair/A7_chk_vs_mrs5probe.rds`.

## 3. The run

- Memory smoke `s5smk1500` (31% HR 1.75 n 1500, 13 replicates, 13 workers): peak 13.48 GB total R RSS over 14 processes,
  ≥ 75% system memory free, i.e. ~1 GB per process; 13 workers projected inside the 24 GB rule, so n > 500 ran at 13.
- Whole run: peak total R RSS 17.18 GB, minimum system free 65% (`s5rerun/logs/s5rerun_mem.log`, 30 s samples). The 20% stop never fired.
- Guard 3× projection (projection = committed wall × 1.32); never reached.

| cell | source | HR | n | wall (s) | projection (s) | workers | replicates | errors | declaration rate |
|---|---|---|---|---|---|---|---|---|---|
| A7 | p12x20 | 1.00 | 500 | 1936 | 1728 | 13 | 2000 | 0 | 0.6805 |
| P31_h100_n500 | cert20 | 1.00 | 500 | 2679 | 2133 | 13 | 2000 | 0 | 0.9205 |
| A1 | p12x20 | 1.50 | 500 | 2439 | 2079 | 13 | 2000 | 0 | 0.9110 |
| A4 | p12x20 | 1.75 | 500 | 2529 | 2136 | 13 | 2000 | 0 | 0.9500 |
| P31_h150_n500 | e1stud | 1.50 | 500 | 3161 | 2453 | 13 | 2000 | 0 | 0.9995 |
| P31_h175_n500 | e1stud | 1.75 | 500 | 3291 | 2509 | 13 | 2000 | 0 | 0.9995 |
| A8 | p12x20 | 1.00 | 1000 | 1967 | 2327 | 13 | 2000 | 0 | 0.6595 |
| P31_h100_n1000 | cert20 | 1.00 | 1000 | 2961 | 3204 | 13 | 2000 | 0 | 0.9545 |
| A2 | p12x20 | 1.50 | 1000 | 2699 | 2892 | 13 | 2000 | 0 | 0.9740 |
| A5 | p12x20 | 1.75 | 1000 | 2800 | 2921 | 13 | 2000 | 0 | 0.9950 |
| P31_h150_n1000 | cert20 | 1.50 | 1000 | 3614 | 3522 | 13 | 2000 | 0 | 1.0000 |
| P31_h175_n1000 | cert20 | 1.75 | 1000 | 3774 | 3556 | 13 | 2000 | 0 | 1.0000 |
| A9 | p12x20 | 1.00 | 1500 | 1976 | 3163 | 13 | 2000 | 0 | 0.6240 |
| P31_h100_n1500 | cert20 | 1.00 | 1500 | 3131 | 4941 | 13 | 2000 | 0 | 0.9590 |
| A3 | p12x20 | 1.50 | 1500 | 2831 | 4425 | 13 | 2000 | 0 | 0.9880 |
| A6 | p12x20 | 1.75 | 1500 | 2900 | 4450 | 13 | 2000 | 0 | 0.9990 |
| P31_h150_n1500 | cert20 | 1.50 | 1500 | 3895 | 5436 | 13 | 2000 | 0 | 1.0000 |
| P31_h175_n1500 | cert20 | 1.75 | 1500 | 4075 | 5474 | 13 | 2000 | 0 | 1.0000 |
| **total** | | | | **52658** | 59349 | | 36000 | 0 | |

Cells not run: none.

## 4. Read-out — against the committed published bundles

### 4.1 Per cell (paired by `sim_id`, 2,000 replicates; shifts on declaring replicates)

"bound crossings" counts declaring-replicate readings of the field lower, field-s upper and the two `joint_s` Bonferroni
bounds that land on the other side of 0.75 or 1.25; "coverage-indicator changes" counts replicates whose covers / misses
flips on any of those four bounds.

| cell | source | HR | n | decl. rate before / after | est shift mean / med / p05 / p95 | field lower shift mean / med / p05 / p95 | field-s upper shift mean / med / p05 / p95 | region chg | bound crossings 0.75/1.25 | coverage-indicator changes | errors |
|---|---|---|---|---|---|---|---|---|---|---|---|
| A7 | p12x20 | 1.00 | 500 | 0.6805 / 0.6805 | +0.00089 / +0.00066 / -0.00021 / +0.00279 | +0.00055 / +0.00041 / -0.00015 / +0.00192 | +0.00048 / +0.00018 / -0.00080 / +0.00276 | 0 | 11 | 6 | 0 |
| P31_h100_n500 | cert20 | 1.00 | 500 | 0.9205 / 0.9205 | +0.00035 / +0.00014 / -0.00015 / +0.00156 | +0.00021 / +0.00008 / -0.00010 / +0.00103 | +0.00019 / +0.00001 / -0.00070 / +0.00188 | 0 | 1 | 2 | 0 |
| A1 | p12x20 | 1.50 | 500 | 0.9110 / 0.9110 | +0.00035 / +0.00015 / -0.00044 / +0.00176 | +0.00019 / +0.00009 / -0.00038 / +0.00121 | +0.00027 / +0.00005 / -0.00068 / +0.00203 | 0 | 3 | 1 | 0 |
| A4 | p12x20 | 1.75 | 500 | 0.9500 / 0.9500 | +0.00019 / +0.00002 / -0.00050 / +0.00137 | +0.00011 / +0.00002 / -0.00039 / +0.00102 | +0.00020 / +0.00001 / -0.00055 / +0.00176 | 0 | 3 | 5 | 0 |
| P31_h150_n500 | e1stud | 1.50 | 500 | 0.9995 / 0.9995 | +0.00003 / +0.00000 / -0.00012 / +0.00030 | +0.00001 / +0.00000 / -0.00007 / +0.00017 | +0.00003 / +0.00000 / -0.00012 / +0.00019 | 0 | 1 | 2 | 0 |
| P31_h175_n500 | e1stud | 1.75 | 500 | 0.9995 / 0.9995 | +0.00000 / +0.00000 / -0.00006 / +0.00007 | +0.00000 / +0.00000 / -0.00003 / +0.00004 | +0.00001 / +0.00000 / -0.00002 / +0.00002 | 0 | 0 | 2 | 0 |
| A8 | p12x20 | 1.00 | 1000 | 0.6595 / 0.6595 | +0.00080 / +0.00055 / -0.00012 / +0.00240 | +0.00064 / +0.00041 / -0.00017 / +0.00230 | +0.00036 / +0.00016 / -0.00043 / +0.00188 | 0 | 7 | 4 | 0 |
| P31_h100_n1000 | cert20 | 1.00 | 1000 | 0.9545 / 0.9545 | +0.00023 / +0.00007 / -0.00011 / +0.00113 | +0.00017 / +0.00004 / -0.00009 / +0.00088 | +0.00013 / +0.00001 / -0.00045 / +0.00121 | 0 | 5 | 2 | 0 |
| A2 | p12x20 | 1.50 | 1000 | 0.9740 / 0.9740 | +0.00015 / +0.00000 / -0.00040 / +0.00111 | +0.00011 / +0.00001 / -0.00039 / +0.00097 | +0.00013 / +0.00002 / -0.00035 / +0.00116 | 0 | 6 | 1 | 0 |
| A5 | p12x20 | 1.75 | 1000 | 0.9950 / 0.9950 | +0.00002 / +0.00000 / -0.00038 / +0.00054 | -0.00000 / +0.00000 / -0.00032 / +0.00044 | +0.00008 / +0.00000 / -0.00018 / +0.00065 | 0 | 3 | 2 | 0 |
| P31_h150_n1000 | cert20 | 1.50 | 1000 | 1.0000 / 1.0000 | +0.00000 / +0.00000 / -0.00001 / +0.00001 | +0.00000 / +0.00000 / -0.00000 / +0.00001 | -0.00000 / +0.00000 / -0.00001 / +0.00001 | 0 | 0 | 0 | 0 |
| P31_h175_n1000 | cert20 | 1.75 | 1000 | 1.0000 / 1.0000 | -0.00000 / +0.00000 / -0.00000 / +0.00000 | -0.00000 / +0.00000 / -0.00000 / +0.00000 | +0.00000 / +0.00000 / -0.00000 / +0.00000 | 0 | 0 | 0 | 0 |
| A9 | p12x20 | 1.00 | 1500 | 0.6240 / 0.6240 | +0.00089 / +0.00074 / -0.00013 / +0.00241 | +0.00080 / +0.00061 / -0.00025 / +0.00263 | +0.00033 / +0.00018 / -0.00057 / +0.00154 | 0 | 12 | 5 | 0 |
| P31_h100_n1500 | cert20 | 1.00 | 1500 | 0.9590 / 0.9590 | +0.00022 / +0.00008 / -0.00007 / +0.00095 | +0.00019 / +0.00006 / -0.00006 / +0.00086 | +0.00012 / +0.00001 / -0.00048 / +0.00118 | 0 | 4 | 6 | 0 |
| A3 | p12x20 | 1.50 | 1500 | 0.9880 / 0.9880 | +0.00006 / +0.00000 / -0.00042 / +0.00077 | +0.00006 / +0.00000 / -0.00039 / +0.00078 | +0.00010 / +0.00001 / -0.00031 / +0.00087 | 0 | 4 | 0 | 0 |
| A6 | p12x20 | 1.75 | 1500 | 0.9990 / 0.9990 | -0.00003 / +0.00000 / -0.00037 / +0.00025 | -0.00003 / +0.00000 / -0.00032 / +0.00024 | +0.00005 / +0.00000 / -0.00013 / +0.00040 | 0 | 6 | 0 | 0 |
| P31_h150_n1500 | cert20 | 1.50 | 1500 | 1.0000 / 1.0000 | +0.00000 / +0.00000 / -0.00000 / +0.00000 | +0.00000 / +0.00000 / -0.00000 / +0.00000 | +0.00000 / +0.00000 / -0.00000 / +0.00000 | 0 | 0 | 0 | 0 |
| P31_h175_n1500 | cert20 | 1.75 | 1500 | 1.0000 / 1.0000 | -0.00000 / +0.00000 / +0.00000 / +0.00000 | -0.00000 / +0.00000 / -0.00000 / +0.00000 | -0.00000 / +0.00000 / -0.00000 / +0.00000 | 0 | 0 | 0 | 0 |

### 4.2 Per-cell coverage (field lower on β(Ĥ), field-s upper on β(Ĥᶜ); detected replicates; `fs_sim_bias_coverage`)

| cell | field lower cov. before | after | field-s upper cov. before | after |
|---|---|---|---|---|
| A7 | 0.9618 | 0.9618 | 0.9331 | 0.9346 |
| P31_h100_n500 | 0.9734 | 0.9734 | 0.9115 | 0.9125 |
| A1 | 0.9704 | 0.9709 | 0.9336 | 0.9336 |
| A4 | 0.9716 | 0.9716 | 0.9321 | 0.9326 |
| P31_h150_n500 | 0.9745 | 0.9745 | 0.9125 | 0.9130 |
| P31_h175_n500 | 0.9700 | 0.9700 | 0.9195 | 0.9200 |
| A8 | 0.9318 | 0.9318 | 0.9553 | 0.9553 |
| P31_h100_n1000 | 0.9560 | 0.9560 | 0.9277 | 0.9277 |
| A2 | 0.9476 | 0.9476 | 0.9528 | 0.9533 |
| A5 | 0.9437 | 0.9437 | 0.9442 | 0.9447 |
| P31_h150_n1000 | 0.9585 | 0.9585 | 0.9420 | 0.9420 |
| P31_h175_n1000 | 0.9525 | 0.9525 | 0.9420 | 0.9420 |
| A9 | 0.9391 | 0.9375 | 0.9519 | 0.9519 |
| P31_h100_n1500 | 0.9666 | 0.9666 | 0.9270 | 0.9270 |
| A3 | 0.9580 | 0.9580 | 0.9474 | 0.9474 |
| A6 | 0.9585 | 0.9585 | 0.9550 | 0.9550 |
| P31_h150_n1500 | 0.9615 | 0.9615 | 0.9465 | 0.9465 |
| P31_h175_n1500 | 0.9610 | 0.9610 | 0.9460 | 0.9460 |

### 4.3 Published ranges and figures (`s5rerun/ranges.R` → `s5rerun/ranges.md`)

The "before" rows reproduce `REPORT_fs_products_reconciliation_2026-09-12.md` (designated-comparator set) and the FS
rows of `REPORT_grfmr_completion_2026-09-12.md` exactly, so the definitions match the published ones.

## Per harm cell (4 dp), before -> after

| cell | source | HR | n | field lower | field-s upper | joint | joint_s |
|---|---|---|---|---|---|---|---|
| A1 | p12x20 | 1.50 | 500 | **0.9704 -> 0.9709** | 0.9336 | **0.9517 -> 0.9523** | 0.9555 |
| A4 | p12x20 | 1.75 | 500 | 0.9716 | **0.9321 -> 0.9326** | 0.9474 | 0.9495 |
| P31_h150_n500 | e1stud | 1.50 | 500 | 0.9745 | **0.9125 -> 0.9130** | **0.9320 -> 0.9315** | **0.9420 -> 0.9425** |
| P31_h175_n500 | e1stud | 1.75 | 500 | 0.9700 | **0.9195 -> 0.9200** | 0.9335 | **0.9395 -> 0.9400** |
| A2 | p12x20 | 1.50 | 1000 | 0.9476 | **0.9528 -> 0.9533** | 0.9487 | 0.9502 |
| A5 | p12x20 | 1.75 | 1000 | 0.9437 | **0.9442 -> 0.9447** | 0.9387 | **0.9402 -> 0.9407** |
| P31_h150_n1000 | cert20 | 1.50 | 1000 | 0.9585 | 0.9420 | 0.9460 | 0.9500 |
| P31_h175_n1000 | cert20 | 1.75 | 1000 | 0.9525 | 0.9420 | 0.9445 | 0.9485 |
| A3 | p12x20 | 1.50 | 1500 | 0.9580 | 0.9474 | 0.9494 | 0.9514 |
| A6 | p12x20 | 1.75 | 1500 | 0.9585 | 0.9550 | 0.9580 | 0.9575 |
| P31_h150_n1500 | cert20 | 1.50 | 1500 | 0.9615 | 0.9465 | 0.9440 | 0.9505 |
| P31_h175_n1500 | cert20 | 1.75 | 1500 | 0.9610 | 0.9460 | 0.9425 | 0.9490 |

## Ranges over the 12 harm cells, before / after

| cell set | field lower | field-s upper | joint | joint_s |
|---|---|---|---|---|
| (a) Section 5 set, before | 0.9437-0.9745 | 0.9125-0.9550 | 0.9320-0.9580 | 0.9395-0.9575 |
| (a) Section 5 set, after | 0.9437-0.9745 | 0.9130-0.9550 | 0.9315-0.9580 | 0.9400-0.9575 |
| (b) designated-comparator set, before | 0.9410-0.9745 | 0.9125-0.9605 | 0.9320-0.9635 | 0.9395-0.9640 |
| (b) designated-comparator set, after | 0.9410-0.9745 | 0.9130-0.9605 | 0.9315-0.9635 | 0.9400-0.9640 |

## HR 1.00 cells: field lower bound location (detected replicates), before / after

| cell | source | n | arm | n_eval | median lower | median theta(Hhat) | bound/theta | paired ratio | share >= 1.00 | share >= 1.25 |
|---|---|---|---|---|---|---|---|---|---|---|
| A7 | p12x20 | 500 | before | 1361 | 0.4204 | 0.6982 | 0.6021 | 0.6067 | 0.0044 | 0.0000 |
| A7 | p12x20 | 500 | after | 1361 | 0.4210 | 0.6982 | 0.6030 | 0.6074 | 0.0044 | 0.0000 |
| P31_h100_n500 | cert20 | 500 | before | 1841 | 0.4591 | 0.8241 | 0.5571 | 0.5637 | 0.0060 | 0.0005 |
| P31_h100_n500 | cert20 | 500 | after | 1841 | 0.4591 | 0.8241 | 0.5572 | 0.5639 | 0.0060 | 0.0005 |
| A8 | p12x20 | 1000 | before | 1319 | 0.5167 | 0.7357 | 0.7023 | 0.7088 | 0.0045 | 0.0008 |
| A8 | p12x20 | 1000 | after | 1319 | 0.5173 | 0.7357 | 0.7031 | 0.7098 | 0.0045 | 0.0008 |
| P31_h100_n1000 | cert20 | 1000 | before | 1909 | 0.5512 | 0.8494 | 0.6489 | 0.6558 | 0.0105 | 0.0016 |
| P31_h100_n1000 | cert20 | 1000 | after | 1909 | 0.5515 | 0.8494 | 0.6493 | 0.6562 | 0.0105 | 0.0016 |
| A9 | p12x20 | 1500 | before | 1248 | 0.5801 | 0.7706 | 0.7529 | 0.7569 | 0.0040 | 0.0000 |
| A9 | p12x20 | 1500 | after | 1248 | 0.5809 | 0.7706 | 0.7539 | 0.7580 | 0.0040 | 0.0000 |
| P31_h100_n1500 | cert20 | 1500 | before | 1918 | 0.6032 | 0.8709 | 0.6926 | 0.7058 | 0.0094 | 0.0010 |
| P31_h100_n1500 | cert20 | 1500 | after | 1918 | 0.6033 | 0.8709 | 0.6927 | 0.7065 | 0.0094 | 0.0010 |

### 4.4 Reading

- **Unchanged:** every declaration, every selected region, every declaration / selection rate, every bound-location
  share (≥ 1.00, ≥ 1.25) at the HR 1.00 cells, and the field-lower harm range (0.9410–0.9745 designated; 0.9437–0.9745 Section 5 set).
- **Changed at four decimals — three range endpoints, one replicate each:**
  - field-s upper low end **0.9125 → 0.9130** (31% HR 1.50 n 500, `e1stud`): the published 0.9125–0.9605 becomes **0.9130–0.9605**.
  - `joint` low end **0.9320 → 0.9315** (same cell): 0.9320–0.9635 becomes **0.9315–0.9635**.
  - `joint_s` low end **0.9395 → 0.9400** (31% HR 1.75 n 500, `e1stud`): 0.9395–0.9640 becomes **0.9400–0.9640**.
- **At three decimals** it depends on the convention the manuscript uses at the low end. The certification record's
  "0.912–0.960" and "0.939–0.963" read the low endpoints downward; under that convention they become **0.913** and **0.940**,
  and `joint` 0.932 becomes **0.931**. Under round-half-up all three read the same before and after (0.913, 0.932, 0.940).
- **Per-cell four-decimal figures also move** where the manuscript quotes them per cell: 14 per-cell coverage values
  (9 field lower / field-s upper across the 18 cells in §4.2, 5 `joint` / `joint_s` on the harm cells in §4.3), by ≤ 0.0016;
  the largest is 12.4% HR 1.00 n 1500 field lower 0.9391 → 0.9375, a `p12x20` cell. Also
  the HR 1.00 median-bound figures (e.g. 31% n 1000 bound/θ 0.6489 → 0.6493, paired ratio 0.6558 → 0.6562).
- **The sweep predicted no change in any stated result; at four decimals the full run does not bear that out.**
  The shift distributions match the sweep's (5th–95th percentiles within ±0.003 on the estimate and both field bounds; the largest single
  shift is 0.020, 12.4% HR 1.75 n 500 sim 710 on the `joint_s` Bonferroni lower bound, against the sweep's 0.015).
  A coverage endpoint moves in steps of one replicate, though, and at 2,000 replicates one replicate crossing its target
  moves it by 0.0005. Every range moves by 0.0005 at most, and no coverage value by more than 0.0016.

## 5. Post-conditions

| # | condition | status |
|---|---|---|
| 1 | 0a ancestors; 0b install, worker-verified, path / version / build recorded | PASS (§0) |
| 2 | cells / parameters / seeds from `mrs5sweep/cells.txt`; seeds `8316951 + sim_id` | PASS |
| 3 | new tag; no committed payload modified | PASS (`s5rerun`; commit `5fec9198` adds 263 files, modifies none) |
| 4 | Step 2 comparison, declaration and region match | PASS (200 / 200, 0 region changes) |
| 5 | bundle structure matches the committed campaign's | PASS (res_1_1000 / res_1001_2000 / combined_1_2000, same stem pattern) |
| 6 | errors 0 per cell | PASS (0 of 36,000) |
| 7 | cells not run listed; complete or partial stated | PASS (none; complete) |
| 8 | no `R/` file modified | PASS |
| 9 | catalogue pin == HEAD at commit time | see the closeout commit (`scripts_dinamr/check_current_status.sh`) |
