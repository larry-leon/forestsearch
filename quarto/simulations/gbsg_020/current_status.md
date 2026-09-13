# current_status — `quarto/simulations/gbsg_020`

- **Pin:** `d62e1391` on `feature/glm-extension` — HEAD at the time this file was committed; the closeout commit that adds this file is its child.
- **Updated:** 2026-09-13, after `TASK_idsweep_2026-09-12`: the Part B identification sweep, 288 cell-runs, complete (§2.6), with Gate I Amendment 1 for review. Before that came `TASK_partB_enabling_2026-09-12` (the `FS_S7_MR` knob, Gate T3 PASS, the 16-run OC smoke) and `TASK_partB_measurement_2026-09-12`, which stopped at Stage 0c.
- **Purpose:** a catalog of what has been run in this directory and where the payloads are, so a chat or workstream on another machine can be brought up to speed by attaching this one file. It points at authoritative files; it does not restate their numbers.
- **Maintenance:** regenerated as the closeout step of every task that touches this directory. The pin above must equal HEAD at commit time. §3 is produced by `scripts_dinamr/status_inventory.R` from the directory.

---

## 1. The DGM, in one block

- GBSG-based survival simulation. One template: `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd`.
- Cells are (prevalence × target HR × n): prevalence 12.4% (`FS_S7_Z1Q` unset) or 31% (`FS_S7_Z1Q=0.60`); HR 1.50, 1.75 or 1.00; n = 500, 1000, 1500. 2,000 replicates per cell unless stated.
- **HR 1.00 is not a global null.** `FS_S7_HR` calibrates `k_inter` to a target Cox HR inside the planted region. The region rule depends on `FS_S7_Z1Q` alone, and `dgm_model <- "alt"` is a literal, so the harness's global-null path is unreachable. At HR 1.00 the planted region carries HR 1.00 against a benefiting complement (HR 0.657 at 12.4%, 0.721 at 31%). Describe it as **differentially null against a benefiting complement**, never with harm vocabulary.
- Seeds are `8316951 + sim_id` throughout, so cells at matched coordinates share DGM draws across campaigns. Verified per cell on both GRF and DINA grids: `n_true` identical on all 2,000 rows, truth agreeing to ≤ 8.9e-15 (cross-machine BLAS; the largest observed is 8.882e-15, on the 31% HR 1.50 / 1.75 cells).

## 2. Campaigns

### 2.1 FS comparators (the reference grid)

| Campaign | Prevalence | Criterion | Cells |
|---|---|---|---|
| `p12ext` | 12.4% | `maxeffCons` ε 0.10 | HR 1.50 (all n); HR 1.00 at n 1000, 1500 |
| `tier2` | 12.4% | `maxeffCons` ε 0.10 | HR 1.75 (all n); HR 1.00 at n 500 |
| `e1stud` | 31% | `effMaxSG` ε 0.20 | HR 1.50 and 1.75 at n 500 |
| `cert20` | 31% | `effMaxSG` ε 0.20 | all others |

- All 18 grid cells are covered, and all 18 designated comparator bundles are verified present. `gate2G.R` now asserts per cell that the designated comparator resolves. The table names the **designated** comparator, the one `gate2G.R` and `fs_extraction.R` resolve to. Earlier FS campaigns (`map1`, `s7`, and others) also hold bundles at some of these cells; those are not the comparator and must not be substituted.
- **Criterion is matched to DINA and GRF at 31% and not at 12.4%.** Any 12.4% cross-identifier gap carries a criterion confound on top of identifier, family construction and detection set.
  - **At 12.4% the confound is the rule, not ε.** `maxeffCons` never reads `effect_neighborhood` on any engine. The "ε 0.10" in the table is the recorded default, and it was inert on those nine cells (`REPORT_partB_measurement_2026-09-12.md`, Stage 0b). Only the informational `band_n` column reads it.
- FS one-sided products, recomputed from the committed bundles across **all 12 harm cells**, both prevalences:
  - **field lower on β(Ĥ): 0.941–0.975.** The certification record's 0.944–0.974 is the same quantity on a **smaller cell set**. It predates `p12ext`, which supplies the three 12.4% HR 1.50 cells, and 12.4% HR 1.50 n 1000 (0.9410) is the only harm cell below 0.944. Both are right for their set; quote the cell set with the range.
  - **field-s upper on β(Ĥᶜ): 0.9125–0.9605.** Identical to the record's "0.912–0.960": the endpoints are the same numbers under a different rounding convention, not a disagreement.
  - **Bonferroni joint: unscaled `joint` 0.932–0.964; studentized `joint_s` 0.940–0.964.** The record's "0.939–0.963" is **`joint_s`**, whose 9-cell range is 0.9395–0.9640. `joint` and `joint_s` are different constructions and must be named when quoted.
  - The **six HR 1.00 cells** are a separate cell set; their field-lower values are listed per cell, with the row set named, in `REPORT_grfmr_completion_2026-09-12.md`. Do not fold them into the harm ranges above.
- Classification and bound-location extraction: `REPORT_fs_extraction_2026-09-11.md` (all 18 cells; reading of committed bundles, no re-run).

### 2.2 `dinamr` — DINA, complete

- Engine `dina`, `effMaxSG`, ε 0.20, effect floor log(0.90), no consistency term. 18 of 18 cells, 2,000 replicates.
- Reports: `REPORT_dinamr_*`, `REPORT_dinamr_blockC_2026-09-11.md`. Summary: `summary_dinamr.qmd` / `.html` (18 cells). Review: `REVIEW_dinamr_blockC_2026-09-11.md`.
- DINA's `n_family` is the **surface-proposed candidate set**: small, volatile and shrinking with n. CV runs 0.296–1.458 across all 18 cells; the 1.008–1.458 band quoted in `REPORT_dinamr_blockC` is the **six HR 1.00 cells only**, where it is most extreme.
- **Detection falls with n**, sharply at 12.4%: 0.875 → 0.839 → 0.799 on HR 1.50 harm, 0.908 → 0.907 → 0.898 on HR 1.75 harm, and 0.714 → 0.525 → 0.344 at the null. The overall range across the 18 cells is 0.344–1.000. Every DINA n-trend conditions on a shrinking detection set; say so when reading one.
- On DINA, `admitted_n` is not recorded (the column postdates this campaign), and `n_family` and p̂ are strongly confounded (ρ = −0.47), so the joint count table is needed to read either stratification.

### 2.3 `grfmr` — GRF, complete

- Engine `grf`, `effMaxSG`, ε 0.20, `grf_selection = "frontier"`, `grf_select_statistic = "effect"`, `dmin.grf = 0.0`. **18 of 18 cells**, 2,000 replicates, matching the DINA grid cell for cell.
  - Ten harm cells under `TASK_grfmr_campaign_2026-09-11` (7.939 h).
  - The two deferred 31% HR 1.75 cells and the six HR 1.00 cells under `TASK_grfmr_completion_2026-09-12` (6.425 h against a 6.744 h projection). Nothing deferred, nothing dropped.
- Reports: `REPORT_grfmr_2026-09-11.md`, `REPORT_grfmr_completion_2026-09-12.md`. Per-cell tables: `TABLES_grfmr_completion_2026-09-12.md` (18 cells; supersedes for coverage the 10-cell `TABLES_grfmr_percell_2026-09-12.md`, whose ten cells it reproduces). Summary: `summary_grfmr.qmd` / `.html` (18 of 18, no guarded skips). Reviews: `REVIEW_grfmr_2026-09-12.md` (the ten Part A cells) and `REVIEW_grfmr_completion_2026-09-12.md` (the eight completion cells; accepts the report, no re-run).
- **GRF's `n_family` is the outcome-independent enumerated pool**, not the qualified set. At each n it is identical quantile-for-quantile across both prevalences and all three hazard ratios. The outcome-dependent quantity is **`admitted_n`**, with ρ(`admitted_n`, `n_family`) between +0.016 and +0.054 across all 18 cells.
  - **Stratify GRF on `admitted_n`, never on `n_family`.**
  - `admitted_n` and p̂ carry independent information on GRF (joint counts near-uniform at harm and null cells alike), the opposite of DINA.
- **Detection on the harm cells is 1.0000 at all six 31% cells and 0.9975–1.0000 at 12.4%**, flat in n, so no GRF harm n-trend carries a detection-conditioning caveat.
- **At the HR 1.00 cells the rate is a `selection_rate`.** It is 0.9970 / 0.9980 / 0.9935 at 31%, at 1 within its interval and flat in n. At 12.4% it is 0.9840 / 0.9575 / 0.8675, falling with n through non-detections. **The 12.4% null cells are the only GRF cells whose rate moves with n by more than 0.01**, so a 12.4% null n-trend carries a detection-conditioning qualification. Bound location, not the selection rate, carries the null-cell question: the ≥ 1.00 and ≥ 1.25 shares are in `REPORT_grfmr_completion_2026-09-12.md` with FS beside GRF.
- **Non-detections are separable on GRF.** At the null cells, 13 of 405 non-detections carry `admitted_n = 0` (the empty-admitted-set path) and 392 carry `admitted_n` NA (the re-selection returned before the count). All are `NO-DETECTION`; none is a CONFIG-ERROR.
- Source trace (Gate 0a, 13 sites): GRF's candidates are enumerated from quantiles of X subject to `n_min`, with DR scores entering only the effect column. MR re-evaluates the full enumerated pool, not the forest-qualified subset, and re-applies qualification per draw. This is recorded as a fact about the construction. It sits against handoff §3's description of GRF; that discrepancy is noted, unresolved, and blocks nothing.
- Two floors on two scales: `dmin.grf` is a DR-score pre-filter in RMST units; the binding effect-scale floor on the re-selection path is `hr.threshold = 0.90`, the same floor DINA carries. GRF is not "unfloored" relative to DINA.

### 2.4 `grfprobe` — GRF cost and mechanism

- Five 36-replicate probes. Cost 13.5–19.8 s per replicate, rising with n and prevalence, **flat in family size** (|ρ| ≤ 0.171).
- **Two mechanisms, kept distinct.**
  - **The frontier band cannot empty at `dmin.grf = 0.0`:** the eligible maximum is never negative (`R/subgroup_consistency_helpers.R:784–785`, `R/grf_subgroup_labels.R:358`, `:377–383`). Measured empty on 0 of 180 probe replicates.
  - **The admission floor can empty, and did** (`R/forestsearch_helpers.R:1642–1651`, `admitted_n <- 0L`): on **13 of the 12,000 HR 1.00 replicates**, all at null cells and all non-detections, and on none of the 24,000 harm replicates.
  - Calling the empty-selection branch "unreachable" without that qualifier is wrong: the band cannot empty, the floor can. Stated the same way in `summary_grfmr.qmd` ("The engine").

### 2.5 Part B — identification only: enabled and smoke run (the sweep is §2.6)

**The template (committed `8fd89e1d`).**
- **`FS_S7_MR`:** default `TRUE`, the former literal; `FALSE` gives identification and classification only. It tags the stem `_nomr` and is recorded in meta as `mr_inference`.
- **The focus guard** admits the six Part B criteria: `effMaxSG`, `maxeffCons`, `effMinSG`, `maxSG`, `minSG`, `maxeff`.
- **`FS_S7_NBHD`**, and a non-default `selection_rule`, are rejected unless the focus is `effMaxSG` or `effMinSG`.
- **With MR off**, the recorder reads `label`, `n_sel`, `n_cons_qual` and `band_n` from the `forestsearch()` result. **`n_family` is NA**: it is MR's own fitted family and cannot be recovered without an `R/` change.

**Gate T3: PASS.**
- **Unset:** all 163 non-timing columns and `truth` are `identical()` to the pre-change run.
- **MR off:** the identification and classification columns are `identical()` to MR on.
- **Supplementary:** the DINA and GRF on/off pairs pass too.

**OC smoke `pBoc`.**
- **Setup:** 12.4%, HR 1.50, n 500, 30 replicates, 16 runs, 318 s total.
- **Detection:** set by the engine and flat across criteria.
- **Size and classification:** move with the rule, `maxSG` > `effMaxSG` > `eff` > `effMinSG` ≈ `minSG`.
- **Criterion pairs:** none is identical. On consistency, `maxeffCons`/`maxeff` agree on 27 of 29 replicates and `effMinSG`/`minSG` on 26 of 29.
- **Against committed MR-on bundles:** the selections are identical to `dinamr` / `grfmr` / `p12ext` on sim_id 1–30.
- **Projected 288 cell-runs at 12 workers:** 23.3 / 12.3 / 7.5 h flat, or 37.8 / 19.6 / 11.1 h adjusted for how cost varies across cells, at 2,000 / 1,000 / 500 replicates.
- **Details:** `REPORT_partB_enabling_2026-09-12.md`.

### 2.6 `idsweep` — the Part B identification sweep, complete

- **What it is.** 288 cell-runs: three engines, six criteria, 18 cells.
  - 500 replicates per run: sim_id 1–500, one batch per run, no combine render.
  - MR off (`_nomr`), tag `idsweep`, `campaign.sh`'s knob set.
  - **Identification and classification only.** No MR, coverage, bias or bounds.
  - It is a self-contained campaign: no earlier bundle is used as data.
- **Criteria.**
  - Consistency: `effMaxSG` ε 0.20, `effMinSG` ε 0.20, `maxeffCons`, `maxeff`, `maxSG`, `minSG`.
  - DINA and GRF: `effMaxSG` ε 0.20, `effMinSG` ε 0.20, `maxSG`, `minSG`, `eff`. Here `eff` stands for `eff` = `maxeff` = `maxeffCons`, run once.
- **Outcome.**
  - 18 of 18 cells; nothing deferred, nothing dropped.
  - 12.01 h wall from the start, against a 10.21 h Gate 1 projection. The miss is render overhead: 44.5 s per render against the 16.4 s assumed.
  - Gate A: 288 of 288 runs. Gate I: 18 of 18 cells.
- **Gate I Amendment 1** (applied unattended; for review).
  - Part 1 stopped at cell 1. DINA `maxSG` NPV was NA because the selection was the whole trial, so NPV is 0/0.
  - Under the amendment, a classification rate may be NA on a detected replicate only where its denominator is zero, and each such case is counted.
  - This occurs on DINA `maxSG` only: 12 cells, 2–299 of 500 replicates, and 159–299 at the six 31% harm cells.
- **Detection is set by the engine.**
  - It is flat across criteria in every cell on DINA and GRF.
  - On consistency it is flat except `maxeff`, which detects or selects on 500 of 500 in every cell.
- **Size ordering** holds in all 18 cells on every engine: `maxSG` > `effMaxSG` > `eff` / `maxeffCons` (> `maxeff` on consistency) > `effMinSG` > `minSG`.
- **Criterion agreement** (identical rule over jointly detected replicates).
  - Consistency `maxeffCons` / `maxeff`: 0.869–1.000.
  - Consistency `effMinSG` / `minSG`: 0.208–0.942. It falls with n and with prevalence at the harm cells.
  - Every DINA and GRF pair is lower (medians ≤ 0.457).
- **Same draws.** `n_true` is `identical()` on sim_id 1–500 to every covering committed bundle at every coordinate.
- **Files.**
  - `REPORT_idsweep_2026-09-12.md`
  - `summary_idsweep.qmd` / `.html`
  - `scripts_dinamr/idsweep*`, `scripts_dinamr/logs/idsweep*` (includes `idsweep_findings.txt` and `idsweep_walls.txt`)
  - `results/*_nomr_idsweep_res_1_500.rds` (288)
  - `idsweep_*.html` (288, 1.09 GB)

## 3. Payload inventory

Regenerated from the directory by `scripts_dinamr/status_inventory.R`.

| path pattern (first match wins) | tracked/disk | total | largest single file | what it is |
|---|---|---|---|---|
| `results/*dinamr*.rds` | 66/66 | 50.54 MB | `dina_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_z1q60_nb20_dinamr_combined_1_2000.rds` 1.54 MB | `dinamr` per-replicate bundles + metas (batch and combined) |
| `results/*grfmr*.rds` | 55/55 | 55.17 MB | `grf_effMaxSG_fb_mr_field_m1_h100_knoise0_n1000_z1q60_nb20_grfmr_combined_1_2000.rds` 1.55 MB | `grfmr` per-replicate bundles + metas (batch and combined), and the `grfmrsmk` Stage 1 smoke |
| `results/*grfprobe*.rds` | 5/5 | 155 KB | `grf_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_nb20_grfprobe_res_1_36.rds` 31 KB | `grfprobe` cost probes, 36 replicates each |
| `results/*_nomr_idsweep_*.rds` | 288/288 | 16.92 MB | `grf_effMaxSG_fb_mr_field_m1_h100_knoise0_n500_z1q60_nb20_nomr_idsweep_res_1_500.rds` 69 KB | `idsweep` per-replicate bundles (MR off, 500 replicates, one batch per run; all three engines) |
| `results/fs_*.rds` | 324/324 | 149.18 MB | `fs_maxeffCons_fb_mr_field_m1_h150_knoise0_n1000_p12ext_combined_1_2000.rds` 1.51 MB | FS bundles — the comparator grid plus every earlier FS campaign |
| `results/*.rds` | 28/28 | 1.02 MB | `grf_eff_fb_mr_m1_h10_knoise0_n500_combined_1_500.rds` 150 KB | other bundles in `results/` |
| `mr_sweep/**` | 210/210 | 15.19 MB | `grf_mr_n500_res.rds` 130 KB | `mr_sweep/` — an earlier seed-table sweep, superseded, kept for provenance |
| `scripts_dinamr/logs/*` | 725/725 | 2.98 MB | `idsweep.driver.part2.log` 203 KB | per-render, driver and gate logs — `WALL_SECONDS` / `CELL DONE wall=` |
| `scripts_dinamr/*.R` | 33/33 | 254 KB | `gate2G.R` 17 KB | drivers, checkers, projections, extractions (R) |
| `scripts_dinamr/*.sh` | 13/13 | 27 KB | `idsweep.sh` 6 KB | render/campaign drivers and the closeout checker (shell) |
| `scripts_dinamr/*.py` | 2/2 | 13 KB | `transplant_grfmr.py` 12 KB | transplant / chunk-diff helpers (Python) |
| `scripts_dinamr/*.cells` | 8/8 | 2 KB | `idsweep.cells` 0 KB | cell lists, one line per cell |
| `scripts_dinamr/*.rds` | 15/15 | 129 KB | `grfmr_tables.rds` 62 KB | saved derived objects (projections, walls, extracted tables) |
| `scripts_dinamr/*.md` | 1/1 | 10 KB | `README.md` 10 KB | `README.md` — the standing rules, and what each script is |
| `scripts_dinamr/*.txt` | 1/1 | 5 KB | `grf_mechanism_output.txt` 5 KB | captured script output (GRF mechanism probe) |
| `dinamr_*.html` | 54/54 | 235.69 MB | `dinamr_C31_h100_n500_combine_1.html` 4.39 MB | `dinamr` batch and combine renders |
| `grfmr_*.html` | 54/54 | 236.03 MB | `grfmr_C124_h100_n1500_batch_1001.html` 4.40 MB | `grfmr` batch and combine renders |
| `grfprobe_*.html` | 5/5 | 21.56 MB | `grfprobe_g_p124_h150_n1500.html` 4.33 MB | `grfprobe` renders |
| `idsweep_*.html` | 288/288 | 1089.06 MB | `idsweep_p124_h100_n1500_dina_minSG.html` 3.81 MB | `idsweep` batch renders, one per cell-run |
| `probe_*.html` | 10/10 | 43.14 MB | `probe_p31_h100_n500.html` 4.33 MB | Gate 1 cost-probe renders (`dinamr` era) |
| `summary_*.html` | 15/15 | 55.85 MB | `summary_grfmr.html` 7.87 MB | summary rendered outputs, all campaigns |
| `summary_*.qmd` | 15/15 | 398 KB | `summary_grfmr.qmd` 93 KB | summary sources, all campaigns |
| `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` | 1/1 | 165 KB | `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` 165 KB | **the** template `dinamr` / `grfmr` / the FS grid all render |
| `sim_*.qmd` | 53/53 | 4.19 MB | `sim_fs_maxeffCons_fb_mr_m1_h10_knoise0_n500_batch_1_100.qmd` 98 KB | other simulation templates (earlier campaigns and variants) |
| `sim_*.html` | 52/52 | 150.72 MB | `sim_fs_maxeffCons_fb_mr_m1_h10_knoise0_n500_batch_1_1000.html` 3.33 MB | renders of those other templates |
| `fs_*.html` | 245/245 | 973.62 MB | `fs_maxeffCons_fb_mr_field_m1_h175_knoise0_n1000_tier2_combine_1_2000.html` 4.34 MB | FS campaign batch/combine renders (`p12ext`, `tier2`, `e1stud`, `cert20`, earlier) |
| `*.html` | 53/53 | 180.74 MB | `t3gate_t3grf.html` 4.31 MB | remaining renders (smoke, gate, dflt, compare; includes the `t3gate_*` and `partBoc_*` renders) |
| `REPORT_*.md` | 68/68 | 1.15 MB | `REPORT_fixedphat_ij2s_2026-09-09.md` 79 KB | REPORT documents |
| `TABLES_*.md` | 2/2 | 83 KB | `TABLES_grfmr_completion_2026-09-12.md` 53 KB | TABLES documents |
| `REVIEW_*.md` | 3/3 | 30 KB | `REVIEW_dinamr_blockC_2026-09-11.md` 13 KB | REVIEW documents |
| `current_status.md` | 1/1 | 28 KB | `current_status.md` 28 KB | this file — the directory's catalog at a pin |
| `*.md` | 1/1 | 2 KB | `payload_runbook_mr_only_20260819.md` 2 KB | other notes in the directory |
| `*.qmd` | 21/21 | 1.44 MB | `gate_d2_cim_unset.qmd` 156 KB | remaining `.qmd` |
| `*.R` | 4/4 | 63 KB | `p12ext_findings.R` 24 KB | top-level ad-hoc R scripts |
| `**` | 1/1 | 5 KB | `compare_1_20_vs_1_500.csv` 5 KB | everything else |
| **total** | **2720/2720** | **3286 MB** | `summary_grfmr.html` 7.87 MB | every file, each counted once |

Sizes are **apparent size** (`st_size`), in MiB/KiB, not disk usage; `du` reports block-allocated size and reads larger for many small files. **Every file is counted exactly once**: the rules are applied first-match-wins and the rows sum to the total. **Files over 50 MB: 0; over 100 MB: 0.** The gitignored `_gateT_pre_template_files/` and `.DS_Store` are excluded throughout. The `current_status.md` row shows this file's size at the pin, before this regeneration.

**Where to start, by question**

| If you need | Read |
|---|---|
| Per-cell coverage, bias, SDs, misses — DINA | `summary_dinamr.html`, or the Block reports |
| Per-cell numbers — GRF, all 18 cells | `TABLES_grfmr_completion_2026-09-12.md`, `summary_grfmr.html` |
| GRF at the HR 1.00 cells: selection rate and bound location, FS beside GRF | `REPORT_grfmr_completion_2026-09-12.md` |
| FS classification and bound location | `REPORT_fs_extraction_2026-09-11.md` |
| Why a quoted FS range differs from the certification record | `REPORT_fs_products_reconciliation_2026-09-12.md` |
| How a number was computed | `scripts_dinamr/` — the summary's own chunks are authoritative over any re-implementation |
| Cost and wall provenance | `scripts_dinamr/logs/`, `projectC.R`, `walls.R`, `wallsGC.R` |
| Part B: how `sg_focus` resolves per engine, MR-on cost bounds, host factor, Monte Carlo resolution | `REPORT_partB_measurement_2026-09-12.md`, `scripts_dinamr/partB_stage0_readout.R` |
| Part B: the `FS_S7_MR` change, Gate T3, the OC table, criterion agreement, the 288-cell-run projection | `REPORT_partB_enabling_2026-09-12.md`, `scripts_dinamr/partBoc_table.R`, `partBoc_checks.R` |
| Part B sweep: detection or selection, subgroup size, classification and criterion agreement for three engines × six criteria × 18 cells (MR off) | `summary_idsweep.html`, `REPORT_idsweep_2026-09-12.md`, `scripts_dinamr/logs/idsweep_findings.txt` |
| Raw per-replicate rows | `results/*<campaign>*.rds` |

## 4. Reading conventions that must travel with these numbers

- Read every bound **by its location** against clinically meaningful effect sizes. Never frame a result as significance at HR = 1.00.
- **Selection rate is not an error rate.** At the differentially-null cells FS, DINA and GRF all select frequently and admissibly; GRF at 31% selects at 1 within its interval. What speaks to an unsupported claim is the share of lower bounds reaching HR 1.00 and 1.25.
- Coverage and bias as one comparative table (cells × estimators), then a short plain-language reading, not narrative prose.
- Marginal SD and error SD side by side. A Gaussian reference on the marginal SD understates the prediction by up to 7.5 points where the target moves with the estimate; both forms are in the summaries with their formulas stated.
- Wilson intervals on every rate. A Wilson interval belongs on pooled subject-level counts or on a rate over replicates; a replicate-mean rate goes beside it without one.
- Cross-identifier comparisons are descriptive, never a ranking, and carry the confound: identifier, family construction, detection set, and at 12.4% the selection criterion.
- Evaluated estimator set: naive, oracle, IJ two-term, field, and Bonferroni for two-subgroup claims. Winner-only and winner-floor IJ variants are closed.
- `n_cons_qual`, `band_n` and `p_star` are structurally NA on DINA and GRF, not failures, and not on account of HR 1.00.
- On DINA and GRF, `eff` stands for `eff` = `maxeff` = `maxeffCons`: all three resolve to the plain effect argmax there, and the stem tags `eff`. Say so in every table.
- **A classification rate can be undefined.** NPV is NA when the selection is the whole trial, and PPV, sensitivity and specificity likewise when their denominator is zero. In `idsweep` this happens on DINA `maxSG` only. Report the count beside the mean; never read it as a recording failure.

## 5. Superseded or known-wrong — do not quote

The first three entries name files that live **outside this repository** (the `fs_glms_interpretable` workstream and the manuscript build). They are listed because their numbers circulate alongside this directory's, not because they are here. Only the entries after them are about files in `gbsg_020`.

- `claude/forestsearch_audit_spec.md` conditional-coverage figures: built from an earlier submission (different title, 35 + 85 pp.). The current build is 30 + 75 pp., dated 2026-08-20.
- Supplement §8.3 prose describing DINA/GRF conditional coverage as lower than FS's or recovering with n: does not agree with the adjacent tables and Figures S5/S7. Quote the tables and figures, not the prose.
- `BRIEF_dinamr_for_fs_glms_interpretable_2026-09-11.md` (v1): §8 listed FS classification metrics and bound-location shares as unavailable; both exist. Superseded by `BRIEF_fs_identifier_for_fs_glms_interpretable_2026-09-12.md`.
- `REVIEW_grfmr_2026-09-12.md` v1 framed its §2 as a blocking decision; withdrawn in the committed version.
- Earlier versions of this file: "`admitted_n` never 0 across 20,000 campaign replicates" is true of the harm cells only (see §2.4); "ρ(`admitted_n`, `n_family`) ≤ +0.042" was the ten-cell bound (now +0.054 over 18); "`grfmr` partial, 10 of 12, six null cells never run" is superseded by §2.3.
- **A reconciliation of the FS one-sided product ranges (2026-09-12; per-cell values in `REPORT_fs_products_reconciliation_2026-09-12.md`).** An earlier pass of this file "corrected" the certification figures; **two of those corrections were wrong and are withdrawn**. The certification records are correct as written and were not edited. What actually differs:
  - **field lower on β(Ĥ)** — a genuine **cell-set** difference. `NOTE_survival_products_2026-09-09.md` reports 0.944–0.974 over its harm cells. Its evidence list is `cert20` / `tier2` / `fixedphat_ij2s` / `field_studentize_e1` / `cimethod_flip` and does **not** include `REPORT_p12ext_2026-09-09`, the campaign that supplies the three 12.4% HR 1.50 cells. Excluding those, the committed bundles give exactly 0.944–0.974; including them gives 0.941–0.975, the 0.9410 coming from 12.4% HR 1.50 n 1000. *(The note says "ten harm cells"; nine is what reproduces the range, and no tenth bundle on disk carries field columns. Unreconciled, and it does not move the range.)*
  - **field-s upper on β(Ĥᶜ)** — **no difference at all.** The value is 0.9125–0.9605; "0.912–0.960" and "0.913–0.961" are the same endpoints rounded differently. The note's per-cell quotes reproduce exactly: 0.912 / 0.942 / 0.947 at 31% HR 1.50, 0.919 / 0.942 / 0.946 at HR 1.75, 0.941 / 0.956 / 0.961 at 12.4% HR 1.75.
  - **Bonferroni joint** — **different construction, not a different number.** The record's 0.939–0.963 is the **studentized** pair `fld_joint_s_bonf_*` (9-cell range 0.9395–0.9640). The 0.932–0.964 quoted against it was the **unscaled** pair `fld_joint_bonf_*`. Always name which.
  - Neither `SUMMARY_survival_properties_2026-09-10.md` nor `REVIEW_certification_2026-09-09.md` is in this repository, so the "0.941–0.980" attributed to them could not be checked. **0.980 does not reproduce from any field-lower computation on the committed bundles** (the nearest 0.98 in the record is the IJ two-sided at 31%, 0.971–0.981).
- **ε 0.10 as part of the 12.4% FS criterion (2026-09-12; `REPORT_partB_measurement_2026-09-12.md`, Stage 0b).** Reported; the records are not edited.
  - **What is wrong.** Several records treat ε 0.10 as operative on the `tier2` / `p12ext` cells:
    - `REPORT_dinamr_blockC_2026-09-11.md:596`: "at half the band width".
    - `summary_grfmr.qmd` captions: `fs_matched` / `criterion_matched` gated on "sg_focus and eps".
    - `REPORT_grfmr_completion_2026-09-12.md:257` and `REPORT_grfmr_2026-09-11.md:418`: "`maxeffCons` ε 0.10 against `effMaxSG` ε 0.20".
    - `fs_extraction.R`'s eps column.
  - **Why.** ε is inert under `maxeffCons`.
  - **Quote it as:** the 12.4% FS comparator differs from DINA and GRF in the **rule** (`maxeffCons`, no band).

## 6. Not derivable from the committed columns

- **DINA non-detection causes.** The recorder returns before `n_family` is written, so an empty proposal and a proposal with nothing admitted are indistinguishable. (On GRF they are separable through `admitted_n`: 0 versus NA.)
- **Per-subject membership**, so no classification metric finer than the four rates, and no agreement between two identifiers' Ĥ on the same draw.
- **2×2 counts** are not stored, but are exactly recoverable from the recorded rates with `n_sel` and `n_true` (checked to 5.68e-14).
- **`dinamr` render logs** were written to a session scratchpad that no longer exists. The compute wall is rebuildable from bundle timing columns (`walls.R`); the per-render overhead for that campaign is not.

## 7. Open work in this directory

- **Part B sweep `idsweep`: complete (§2.6); two items for review.**
  - **Gate I Amendment 1**, applied unattended after the first Gate I stop. It allows a classification rate to be NA only where its denominator is zero, and each case is counted. It needs Larry's acceptance or reversal; rationale and negative test are in `REPORT_idsweep_2026-09-12.md`.
  - **Push size.** The sweep adds 288 renders totalling 1.09 GB (largest 3.81 MB; nothing over 50 MB) to the pushable range.
  - Still not done: identification with MR on, and any criterion choice. The sweep carries no acceptance criterion and no recommendation.
- A criterion-matched FS comparator at 12.4% **with MR products**. None is committed.
  - `idsweep` now supplies FS identification and classification at 12.4% under `effMaxSG` ε 0.20 as well as `maxeffCons`, at every cell, MR off.
  - The MR-on comparison still needs compute.
- Whether DINA and GRF should default to the field constructions, now informed by two complete 18-cell grids.
- Whether to pin a commit in `forestsearch_version`.
- The frontier-filter asymmetry decision (GRF's band as a filter with no empty-band fallback, where DINA uses a sort key and MR's `.inband()` has a "never empty" fallback) remains open; §2.4 states what can and cannot empty.
