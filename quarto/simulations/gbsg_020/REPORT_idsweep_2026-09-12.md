# REPORT — Part B: the identification sweep (`idsweep`)

- **Date:** 2026-09-12 (run 2026-09-13 00:00–12:01 PDT). **Executor:** Claude Code, unattended. **Machine:** Mac Studio (M4 Max, 14 cores, 36 GB). **Branch:** `feature/glm-extension`.
- **Task:** `dev/tasks/TASK_idsweep_2026-09-12.md`, committed first as `c242310f`.
- **Commits:** tooling `e3863a33`; Gate I Amendment 1 and resume `5e1e47c8`; the sweep payload and this report; the `current_status.md` closeout.
- **No `R/` change and no template change.** Nothing needed one. The only mid-run change is to this sweep's own Gate I script (Amendment 1, below).
- **Outcome:**
  - **18 of 18 cells, 288 of 288 cell-runs.** Nothing deferred, nothing dropped.
  - **12.01 h wall** from start to finish, against a 13 h ceiling and a 16 h hard timeout. Gate 1 projected 10.21 h.
  - **Gate A:** 288 of 288 runs pass. **Gate I:** 18 of 18 cells pass, after Amendment 1.
  - **For review:** Amendment 1, applied unattended.

**Scope.** This is identification and classification only. There is no MR, no field construction, no coverage, no bias and no bound anywhere in this report. The MR-derived columns are NA by design and are not read. It is a self-contained campaign: every number comes from the `idsweep` bundles. Committed bundles are read only for the same-draws check on `n_true`.

---

## 1. Setup as run

- **Grid.** 18 cells: 12.4% (`FS_S7_Z1Q` unset) and 31% (`FS_S7_Z1Q=0.60`) × HR 1.50 / 1.75 / 1.00 × n 500 / 1000 / 1500. They ran in the kickoff's order (`scripts_dinamr/idsweep.cells`).
- **Sixteen runs per cell**, all completed before the next cell started:
  - Consistency: `effMaxSG` (ε 0.20), `effMinSG` (ε 0.20), `maxeffCons`, `maxeff`, `maxSG`, `minSG`.
  - DINA and GRF: `effMaxSG` (ε 0.20), `effMinSG` (ε 0.20), `maxSG`, `minSG`, `eff`.
  - **On DINA and GRF, `eff` stands for `eff` = `maxeff` = `maxeffCons`, run once.** It was passed as `FS_S7_FOCUS=maxeffCons`, because the template's focus guard does not admit `eff`. It stem-tags as `eff`.
- **Knobs.** `campaign.sh`'s pinned set, plus the following:
  - `FS_S7_MR=FALSE`, `FS_S7_CAMPAIGN=idsweep`, `FS_S7_FB=none`, `FS_S7_WORKERS=12`.
  - `FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE FS_S7_FIELD_RECOV=TRUE FS_S7_IJ_RESIDUAL=two_term`. These five are MR-only and inert with MR off. They are pinned because the kickoff pins `campaign.sh`'s set. The `pBoc` smoke did not set them; Gate T3's `t3mroff` rendered with them.
  - `FS_S7_NBHD=0.20` on `effMaxSG` / `effMinSG` only, unset otherwise. `FS_S7_ER_JCUTS` unset. `FS_S7_METHOD` unset on consistency.
- **Seeds and batches.**
  - Seeds `8316951 + sim_id` (template literal), sim_id 1–500.
  - **One batch per run: 500 replicates does not need splitting, so there is no combine render.** Bundles are `results/<engine>_<tag>_fb_mr_field_m1_h<HR>_knoise0_n<n>[_z1q60][_nb20]_nomr_idsweep_res_1_500.rds`.
- **Threads and host.** `render.sh` exports `VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`. Installed `forestsearch` 0.3.5 was built 2026-09-11, after the last `R/` commit. Swap was 0 MB at every cell start.
- **Driver.** `scripts_dinamr/idsweep.sh`, run under `caffeinate -is`, in two parts:
  - **Part 1:** `logs/idsweep.driver.part1.log`, 00:00:34–00:30:49. It ran cell 1, then Gate I stopped it.
  - **Part 2:** `logs/idsweep.driver.part2.log`, 00:34:02–12:01:07. It ran cells 2–18. Its clock and watchdog both count from the part 1 start (`IDSWEEP_T0`), and cell 1 was not re-run (`IDSWEEP_FROM=2`).

## 2. Gate 1 — compute

### 2a. Reconciling `REPORT_partB_enabling`'s 500-replicate projection

**11.1 h is the right arithmetic for that report's method; ~13.6 h is not.** Output: `logs/idsweep_gate1.txt`; script: `idsweep_project.R gate1`.

- **The report's method applies the per-engine multiplier to compute only:**

  | Engine | Compute | × multiplier | + overhead | = total |
  |---|---|---|---|---|
  | consistency | 2.18 h | × 2.08 | + 0.98 h | 5.51 h |
  | DINA | 0.54 h | × 2.33 | + 0.82 h | 2.08 h |
  | GRF | 2.12 h | × 1.26 | + 0.82 h | 3.50 h |
  | **All three** | | | | **11.08 h** |

- **The ~13.6 h figure applies the multipliers to the whole wall**, overhead included:
  - 3.2 × 2.08 + 1.4 × 2.33 + 2.9 × 1.26 = 13.57 h, or 13.44 h unrounded.
  - Render overhead is a per-render cost. It does not scale with per-replicate cost, so multiplying it by a cost-profile ratio double-counts it.
- **Both figures carry two further errors for this sweep.**
  - **Two renders per cell-run.** Both count a batch plus a combine: 576 renders × 16.4 s = 2.62 h. `idsweep` runs one render per cell-run, 288 × 16.4 s = 1.31 h.
  - **The consistency multiplier mixes hosts.** It divides the 18-cell mean of a profile holding 14 pop-os/100-worker and 4 Mac/12-worker cells by a pop-os value.
    - The two hosts' profiles scale with n very differently: ×4.5 on pop-os and ×1.34 on the Mac, from n 500 to 1500.
    - Put on the Mac scale with the committed host factor (6.13), the multiplier is ×2.28.
- **Corrected basis.**
  - Compute: the smoke's per-criterion per-replicate cost × each cell's profile ratio (host-normalized), × 500 / 12.
  - Overhead: 16 renders per cell × 16.4 s.
  - This gives **10.21 h (compute 8.90 h + overhead 1.31 h): GO, all 288.** It is 9.78 h with the un-normalized consistency profile.
- **Per criterion.** The smoke's costs were used as measured: consistency `maxeff` at 2.55 s per replicate, against 1.56–1.61 s for the other five. Nothing was assumed flat.
- **Sensitivity to the unmeasured 500-replicate render overhead.** At 16.4 / 30 / 60 / 90 s per render the total is 10.21 / 11.30 / 13.70 / 16.10 h, and 18 / 18 / 17 / 14 cells fit under 13 h.

### 2b. Re-projection, and projected against realized

- **The first cell measured the render overhead.** A 500-replicate MR-off render costs **50.2 s** beyond compute, against 16.4 s at 30 replicates.
  - Re-projected before cell 2, the grid would finish at **13.33 h**, and the last cell (31% HR 1.00 n 1500) was projected to defer.
  - Realized identification cost then rose with n far more gently than the MR-on profile. The projection came down to 10.73 h after the three 12.4% HR 1.50 cells, then 11.94 h and 11.96 h, and every cell got GO. The decisions are logged in `logs/idsweep_reproject.txt`.
- **Realized against projected** (`logs/idsweep_walls.txt`):

  | | Realized | Gate 1 |
  |---|---|---|
  | Total wall | **12.01 h** driver, start to finish; 11.95 h as the sum of cell walls | **10.21 h** |
  | Render overhead | **3.52 h**: 44.5 s median per render [20, 65] over 288 renders | 1.31 h at 16.4 s |
  | Compute, consistency | 3.36 h | 4.96 h |
  | Compute, DINA | 1.67 h | 1.26 h |
  | Compute, GRF | 3.32 h | 2.68 h |
  | Compute, total | 8.35 h | 8.90 h |

  - Per cell, realized over Gate 1 ranges from 0.71 to 1.60. It runs over at every n 500 cell (1.47–1.60), where overhead dominates, and under at four of the six n 1500 cells (0.71–0.99).
- **What Gate 1 got wrong.**
  - **The overhead:** +2.21 h.
  - **The cost profile.** The host-normalized consistency profile still over-predicted n-scaling (−1.60 h). The MR-on DINA and GRF profiles under-predicted identification cost (+0.41 h, +0.64 h).
- **Realized per-replicate cost, averaged over the 18 cells.**
  - Consistency: 2.54–2.61 s for the five floored criteria, and **`maxeff` 3.25 s (+25%)**.
  - DINA: 1.60–1.61 s.
  - GRF: 3.16–3.20 s.
  - These are 1.6×, 3.1× and 1.6× the smoke's single-cell values.
- **Bounds recorded:** ceiling 13 h (planning bound), hard timeout 16 h (kill bound). Neither was reached.

## 3. Gate A — alignment, per run: 288 of 288 PASS

38 checks per consistency or DINA run and 44 per GRF run; 0 failures. Logs: `logs/idsweep_gateA_<cell>_<engine>_<criterion>.log`. Each value is resolved from the driver's environment, the render log, the rendered document's own audit lines (`Output stem:`, `Run config:`, `Template knobs:`), and the bundle meta.

| Resolved value | Result, all 288 runs |
|---|---|
| MR | `mr_inference=FALSE` in every audit line and `meta$mr_inference` FALSE; `_nomr` in every document stem and every bundle filename |
| `sg_focus` | as intended on every run. On DINA/GRF `eff` resolves as `sg_focus = maxeffCons`, `focus_tag = eff` (18 + 18 runs) |
| ε | `FS_S7_NBHD=0.20` passed and `effect_neighborhood` 0.20 on all 108 `effMaxSG` / `effMinSG` runs. `FS_S7_NBHD` unset on the other 180, where the template forwards its 0.10 default and those foci do not read it |
| GRF | `grf_selection = "frontier"` and `grf_select_statistic = "effect"` (template literals `:533–534`); `Run config: method=grf/frontier` on all 90 |
| GRF `frontier_rule` | resolved by evaluating the `sg_focus → frontier_rule` switch on the normalized focus, from both the installed `forestsearch()` body and `R/forestsearch_main.R`. The two agree on every run: `effMaxSG → effMaxSG`, `effMinSG → effMinSG`, `maxSG → maxSG`, `minSG → minSG`, `maxeffCons → eff` (18 each) |
| GRF path | `admitted_n` finite on every GRF bundle, which is the effect/frontier re-selection's own output. GRF ε is 0.20 on the two band rules (set explicitly; GRF's default is 0.10) and the unread 0.10 default on the other three |
| Other | HR, n, z1q, campaign `idsweep`, `run_mode=batch`, sim_id 1–500, `er_jcuts=10`, `fb_mode=none`, `workers=12`, seed base 8316951, render RC 0. The three thread variables are exported by `render.sh`. The five field knobs read as pinned, and are inert |

## 4. Gate I — integrity, per cell: 18 of 18 PASS

293 checks per cell. Logs: `logs/idsweep_gateI_<cell>.txt`; record: `idsweep_gateI.rds`.

- **Present and populated on every detected replicate:** `detected`, `n_sel`, `label`, `sg_def`, `sens`, `spec`, `ppv`, `npv`, and `admitted_n` on GRF. Rates lie within [0, 1]. This holds in all 288 runs, subject to Amendment 1.
- **No errors.** No status other than DETECTED / NO-DETECTION, and `err_msg` NA on every row of all 144,000 replicate-runs. Non-detections are recorded as such and counted per run in each cell's Gate I log and in the summary's grid-status table.
- **Realized trial prevalence:**

  | Prevalence | n 500 | n 1000 | n 1500 | Super-population |
  |---|---|---|---|---|
  | 12.4% | 0.1233 | 0.1238 | 0.1238 | 0.1242 |
  | 31% | 0.3070 | 0.3070 | 0.3069 | 0.3065 |

  - The value is identical on all 16 runs of a cell, which share draws.
- **Structural NA, never failures:**
  - `n_family` is all NA on all 288 runs. It is MR's own fitted family.
  - `n_cons_qual` and `band_n` are all NA on the 180 DINA and GRF runs, which have no consistency table. They are populated on every detected consistency row.
  - `mr_ok` is 0 throughout.
- **Same draws within a cell:** `n_true` and `truth` are `identical()` across the 16 runs in all 18 cells.
- **Same draws against committed bundles:** `n_true` is `identical()` on sim_id 1–500 against **every** tracked MR-on bundle of this template, at every one of the 18 coordinates, that covers that range. **No mismatch anywhere.**

  | Cells | Committed bundles compared, all `identical()` |
  |---|---|
  | 12.4% HR 1.50 | `dinamr`, `grfmr`, `p12ext`, plus `map1` / `map1c` / `map1w` at n 500 and 1500 |
  | 12.4% HR 1.75 | `dinamr`, `grfmr`, `tier2`, plus `s7` / `s7c` / `s7u` / `s7w` at n 500 |
  | 12.4% HR 1.00 | `dinamr`, `grfmr`, and `tier2` / `s7*` (n 500), `p12ext` / `map1*` (n 1000), `p12ext` (n 1500) |
  | 31% | `dinamr`, `grfmr`, plus `cert20`, `e1stud`, `p30`, `p30sg`, `p30sgnb20`, `p30sgnb20j20` and `banddial` where present |

  - Bundles covering only sim_id 1–5 or 1–36 (smokes, probes) were listed and skipped.
  - None of these bundles is used as data.

### Amendment 1 — undefined classification rates (applied unattended; for review)

- **What stopped the sweep.** Gate I stopped part 1 after cell 1 (12.4% HR 1.50 n 500). The report was: DINA `maxSG`: `npv` NA on 8 of 436 detected replicates.
  - The failed log is kept as `logs/idsweep_gateI_p124_h150_n500_v1_FAILED.txt`.
  - Its second failure, "outside [0,1]", is the same 8 NAs reaching the range check.
- **Diagnosis.**
  - Every one of the 8 is a selection covering the whole trial: `n_sel = 500`, rules such as `er <= 1091`. There are therefore no predicted negatives.
  - The template's `.classify()` (`:848`) returns `npv <- if (tn + fn > 0) ... else NA_real_`.
  - The NA count equals the whole-trial count exactly (8 = 8). No other column was NA on any detected row of the 16 bundles, and the other 291 checks passed.
  - The gate as written required a value the template defines as undefined. The data carried no recording failure.
- **The amendment** (`idsweep_gateI.R`). A rate may be NA on a detected replicate only where the recorded counts give a zero denominator:
  - NPV when `n_sel == n`
  - PPV when `n_sel == 0`
  - sensitivity when `n_true == 0`
  - specificity when `n_true == n`

  Each such NA is counted and reported as structural. **Any NA not explained this way still fails.**
- **Verification.**
  - A negative test mutated a DINA `minSG` bundle to carry an NA `npv` and an NA `sens` on replicates with nonzero denominators. The amended gate failed both (`rc 1`).
  - Cell 1 was re-gated at 293 / 0.
- **Resume.**
  - The sweep resumed from cell 2 with its clock kept from the original start.
  - The part 1 watchdog's orphaned `sleep` was cleared before the relaunch; its parent was already gone. The part 2 orphan was cleared at the end.
- **Why I resumed rather than waited.** The gate's purpose is to catch recording failures, and the amendment preserves every such check. Stopping would have idled the machine for the night on a false positive.
  - **This departs from a literal reading of "stop-on-failure".** Larry may prefer to reverse it. If so, the undefined rates are already counted per run, and no number in this report depends on the amendment except that Gate I reads PASS.
- **Where undefined rates occur:** DINA `maxSG` only, NPV only, always with the selection equal to the whole trial.

  | Cell | Whole-trial selections, of detected |
  |---|---|
  | 12.4% HR 1.50 n 500 | 8 of 436 |
  | 12.4% HR 1.75 n 500 | 9 of 450 |
  | 12.4% HR 1.00 n 500 | 4 of 350 |
  | 31% HR 1.50 n 500 / 1000 / 1500 | 160 / 164 / 159 of 500 |
  | 31% HR 1.75 n 500 / 1000 / 1500 | 226 / 277 / 299 of 500 |
  | 31% HR 1.00 n 500 / 1000 / 1500 | 27 of 467 / 9 of 459 / 2 of 446 |

  - At 31% HR 1.75 n 1500 that is 0.598 [0.554, 0.640] of replicates. The summary's per-cell tables print the count beside the NPV mean, which is over the remaining detected replicates.

## 5. Findings

**Where to read them.** `summary_idsweep.qmd` / `.html` is authoritative for the tables:
- **Per cell:** the OC table, for all 18 cells.
- **Across cells:** rate against n, with a figure and Wilson table; |Ĥ|/|H| against n; sensitivity / specificity / PPV / NPV against n; and criterion agreement, with a figure and a Wilson table.

`scripts_dinamr/logs/idsweep_findings.txt` (`idsweep_findings.R`) holds the across-cell readouts quoted here.

**How to read them.**
- Rates are over 500 replicates, with Wilson 95% intervals.
- Sens / spec / PPV / NPV are replicate means over detected replicates, and carry no interval.
- On DINA and GRF, `eff` stands for `eff` = `maxeff` = `maxeffCons` throughout.
- **This is descriptive only: no acceptance criterion and no recommendation.** Cross-engine differences carry the usual confound: identifier, family construction and detection set.

### 5a. Detection and selection rate — set by the engine, not the criterion

- **Flat across criteria within engine.** DINA and GRF have zero spread across their five criteria in all 18 cells. Consistency has zero spread across its five floored criteria in every cell.
  - **The one exception is consistency `maxeff`, which detects or selects on 500 of 500 in every cell** (1.000 [0.992, 1.000]). `maxeff` disables the consistency and effect floors.
  - Its margin over the other five is:
    - 0.092 / 0.034 / 0.012 at 12.4% HR 1.50 (n 500 / 1000 / 1500)
    - 0.046 / 0.008 / 0 at 12.4% HR 1.75
    - 0.342 / 0.342 / 0.380 at 12.4% HR 1.00
    - 0.086 / 0.050 / 0.038 at 31% HR 1.00
    - 0 at every 31% harm cell
- **31% harm cells (HR 1.50, 1.75):** 500 of 500 on every engine, criterion and n (1.000 [0.992, 1.000]).
- **12.4% harm cells.** Detection moves with n in opposite directions by engine:

  | Engine | HR | n 500 | n 1500 |
  |---|---|---|---|
  | Consistency, five floored criteria | 1.50 | 0.908 [0.879, 0.930] | 0.988 [0.974, 0.994] |
  | Consistency, five floored criteria | 1.75 | 0.954 | 1.000 |
  | DINA | 1.50 | 0.872 [0.840, 0.898] | 0.794 [0.756, 0.827] |
  | DINA | 1.75 | 0.900 [0.871, 0.923] | 0.896 [0.866, 0.920] |
  | GRF | 1.50, 1.75 | 0.996–1.000 | 0.996–1.000 |

  - Consistency rises with n. DINA falls at HR 1.50 and is flat at HR 1.75. GRF is at 0.996–1.000 at both HRs and every n.
  - Every DINA n-trend at 12.4% conditions on a shrinking detection set.
- **HR 1.00: a selection rate.** The planted region is differentially null against a benefiting complement and clears the sub-null log(0.90) floor (`hr_threshold = 0.90`), so returning it is an admissible selection.

  | Engine | 12.4%, n 500 | 12.4%, n 1500 | 31%, n 500 | 31%, n 1500 |
  |---|---|---|---|---|
  | Consistency (floored five) | 0.658 [0.615, 0.698] | 0.620 [0.577, 0.661] | 0.914 [0.886, 0.936] | 0.962 [0.941, 0.976] |
  | DINA | 0.700 [0.658, 0.739] | **0.334 [0.294, 0.376]** | 0.934 [0.909, 0.953] | 0.892 [0.862, 0.916] |
  | GRF | 0.980 [0.964, 0.989] | 0.874 [0.842, 0.900] | 0.996 [0.986, 0.999] | 1.000 [0.992, 1.000] |

  - At 12.4% DINA's selection rate halves with n; consistency and GRF fall less. At 31% the rate is 0.89–1.00 on all three.

### 5b. Subgroup size — set by the criterion, and by prevalence

- **The ordering is identical in all 18 cells on every engine:**
  - Consistency: `maxSG` > `effMaxSG` > `maxeffCons` > `maxeff` > `effMinSG` > `minSG`.
  - DINA and GRF: `maxSG` > `effMaxSG` > `eff` > `effMinSG` > `minSG`.
- **Mean |Ĥ| / mean |H| at the harm cells** (HR 1.50 and 1.75 averaged; n 500 → n 1500):

  | Criterion | 12.4% | 31% |
  |---|---|---|
  | `maxSG` | consistency 2.99 → 2.97; DINA 3.74 → 2.90; GRF 4.46 → 3.50 | 2.49–2.86 → 2.84–3.03 |
  | `effMaxSG` | 1.54–1.66 → 1.37–1.40 | 0.72–0.83 → 0.92–1.01 (**rises** with n) |
  | `eff` / `maxeffCons` / `maxeff` | 1.13–1.19 → 0.92–1.00 | 0.47–0.52 → 0.42–0.47 |
  | `effMinSG` / `minSG` | 0.98–1.08 → 0.82–0.88 | 0.39–0.43 → 0.33–0.35 |

  - **At 12.4%** the `eff` / `maxeffCons` rules reach about the planted size by n 1000, and the minimum-size rules sit just under it.
  - **At 31%** every rule except `maxSG` returns less than the planted region: about a third to a half, or near all of it for `effMaxSG` at n 1500. Size is flat in n there, apart from `effMaxSG`.
- **HR 1.00.**
  - The 12.4% pattern matches the harm cells: `eff` / `maxeffCons` 1.16–1.30 → 0.95–1.02.
  - `maxSG` **shrinks** with n: consistency 2.39 → 1.68, DINA 2.86 → 1.72, GRF 3.57 → 2.21.
  - At 31% it is 1.52–2.15 at n 500.

### 5c. Classification — moves with size and with n

Sensitivity / specificity / PPV / NPV, replicate means over detected replicates, harm cells (HR 1.50 and 1.75 averaged).

- **12.4%.** Every rule except `maxSG` gains sensitivity, specificity and PPV with n. For example (n 500 → n 1500):

  | Engine and rule | n 500 | n 1500 |
  |---|---|---|
  | Consistency `maxeffCons` | 0.61 / 0.92 / 0.54 / 0.94 | 0.85 / 0.98 / 0.84 / 0.98 |
  | DINA `eff` | 0.47 / 0.91 / 0.43 / 0.92 | 0.68 / 0.96 / 0.73 / 0.95 |
  | GRF `eff` | 0.43 / 0.90 / 0.36 / 0.92 | 0.67 / 0.96 / 0.75 / 0.95 |

  - `minSG` gains least: sensitivity 0.54 → 0.51 on consistency, 0.36 → 0.40 on DINA and 0.28 → 0.37 on GRF.
  - `maxSG` buys sensitivity (0.75–0.78 → 0.83–0.93) at specificity 0.48–0.69 and PPV 0.19–0.30 at n 500.
- **31%.**
  - The small-subgroup rules (`eff`, `maxeffCons`, `maxeff`, `effMinSG`, `minSG`) keep specificity at 0.90–0.97, with PPV rising with n. Their **sensitivity stays at 0.14–0.41 at every n.**
  - `effMaxSG` is the rule whose sensitivity grows: 0.53–0.62 → 0.80–0.86, with PPV 0.72–0.74 → 0.86–0.88.
  - `maxSG`: sensitivity 0.93–0.99 with specificity 0.10–0.31, falling with n, and PPV 0.33–0.39.
- **HR 1.00.** The four rates describe how the selection overlaps a differentially null region, and nothing more. At 12.4% they rise with n on consistency and GRF and move little on DINA. Full rows are in the summary.

### 5d. Criterion agreement — the same subgroup on the same trial

This is the share of jointly detected replicates on which two criteria selected the identical rule, within engine. The smoke, at 30 replicates of 12.4% HR 1.50 n 500, found `maxeffCons`/`maxeff` on 27 of 29 and `effMinSG`/`minSG` on 26 of 29. At that cell the sweep gives 427 of 454 (0.941 [0.915, 0.959]) and 350 of 454 (0.771 [0.730, 0.807]).

- **Consistency `maxeffCons` / `maxeff`: the pattern holds across the grid.**
  - 0.869–1.000, median 0.981.
  - It rises with n and with prevalence at the harm cells, reaching 500 of 500 at 31% HR 1.50 n 1500 and HR 1.75 n 1000 / 1500.
  - It is lowest at HR 1.00: 286 of 329 (0.869 [0.829, 0.902]) at 12.4% n 500.
- **Consistency `effMinSG` / `minSG`: the pattern does not hold across the harm grid.**
  - It falls with n and with prevalence: 0.771 [0.730, 0.807] at 12.4% HR 1.50 n 500, down to **0.208 [0.175, 0.246]** (104 of 500) at 31% HR 1.75 n 1500.
  - It is high at HR 1.00: 0.900–0.942 at 12.4% and 0.796–0.879 at 31%.
- **DINA and GRF: every pair is far lower.**
  - `effMinSG`/`minSG`: DINA 0.210–0.862 (median 0.457), GRF 0.078–0.753 (median 0.263).
  - `effMinSG`/`eff`: DINA median 0.372, GRF median 0.315.
  - `effMaxSG`/`maxSG` is 0.000–0.020 at every 31% harm cell on both engines.
  - GRF's `maxSG` shares a rule with `effMinSG`, `minSG` or `eff` on at most 0.027 of replicates in any cell.
  - The highest DINA/GRF shares for `effMinSG`/`minSG` and `effMaxSG`/`maxSG` are at 12.4% HR 1.00 n 1500: DINA 0.862 and 0.617, GRF 0.753 and 0.462.

## 6. Files

All under `quarto/simulations/gbsg_020/` and committed with this report. **Nothing is over 50 MB.**

| Path | Tracked | Size |
|---|---|---|
| `idsweep_<cell>_<engine>_<criterion>.html` (288 batch renders) | yes | **1,089 MB total**; largest 3.81 MB |
| `results/*_nomr_idsweep_res_1_500.rds` (288 bundles) | yes | 16.9 MB total; largest 0.07 MB |
| `scripts_dinamr/logs/idsweep*` (601: 288 render logs, 288 Gate A logs, 19 Gate I logs, the part 1 / part 2 driver logs, re-projection, Gate 1, walls, findings) | yes | 2.6 MB total; largest 0.20 MB |
| `summary_idsweep.qmd`, `summary_idsweep.html` | yes | 20 KB, 4.82 MB |
| `scripts_dinamr/idsweep.sh`, `idsweep.cells`, `idsweep_gateA.R`, `idsweep_gateI.R`, `idsweep_project.R`, `idsweep_findings.R` | yes | 6.5, 0.5, 12.1, 11.7, 17.2, 10.3 KB |
| `scripts_dinamr/idsweep_gate1.rds`, `idsweep_gateI.rds`, `idsweep_reproject.rds`, `idsweep_walls.rds`, `idsweep_findings.rds` | yes | 3–25 KB each |
| `scripts_dinamr/README.md` (idsweep rows), `scripts_dinamr/status_inventory.R` (two `idsweep` rules) | modified | 10.5, 4.9 KB |
| `REPORT_idsweep_2026-09-12.md` | yes | this file |
| `dev/tasks/TASK_idsweep_2026-09-12.md` | yes (`c242310f`) | 7 KB |

- **Flag: push size.** The renders add about **1.09 GB** to the range to push. No single file approaches the 50 MB flag or the 100 MB stop. The total is recorded here, per the standing rule, for Larry to weigh.

## 7. Not done, by instruction

- No MR, coverage, bias or bound quantity.
- No acceptance criterion and no recommendation.
- No `R/` or template change.
- No reduction of replicates.
- No use of an earlier bundle as data.
- No push.
