# REPORT — Part B enabling change, and the OC smoke

- **Date:** 2026-09-12. **Executor:** Claude Code, unattended. **Machine:** Mac Studio. **Branch:** `feature/glm-extension`.
- **Task:** `dev/tasks/TASK_partB_enabling_2026-09-12.md`, committed first as `46167967`.
- **Outcome:**
  - The template change is made and committed as `8fd89e1d`.
  - **Gate T3 passed both halves.**
  - The 16-run OC smoke ran in **318 s**, against a 90-minute cap.
  - No sweep cell was run. **No `R/` change**: nothing needed one.

---

## Part 1 — the change

The only file changed is `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` (+89 / −8). Line numbers below are post-change.

### 1a. `FS_S7_MR`

- **Knob.** `mr_inference_on <- .env_chr("FS_S7_MR", "TRUE")`, validated to `TRUE` / `FALSE` (`:378–380`).
  - It is read beside the prevalence knob, because `FALSE` tags the stem.
  - `forestsearch(mr_inference = mr_inference_on)` (`:1054`) replaces the literal `TRUE`.
- **Stem.** `mr_tag` is `""` when MR is on and `"_nomr"` when it is off (`:449–450`). Every existing stem is unchanged, and MR-off bundles can never share a stem or `combine_glob` with MR-on ones.
- **Echo.**
  - The knob audit line gains `mr_inference=%s` (`:728`).
  - The settings readout gains the row "mr_inference (FS_S7_MR)" (`:1771`). When MR is off it reads "FALSE (identification only; every MR setting below is inert)".
- **Meta.**
  - Batch meta records `mr_inference = mr_inference_on` (`:1610`).
  - Combined meta records `mr_inference = bundles[[1]]$meta$mr_inference %||% NA` (`:1728`).
  - Both are record-only and follow the template's provenance convention. The field is not added to the poolability key, because the stem token already separates the two states.

### 1b. The focus guard

- **Six criteria.** `stopifnot(sg_focus %in% c("effMaxSG", "maxeffCons", "effMinSG", "maxSG", "minSG", "maxeff"))` (`:333`). It is still a guard: anything else stops the render.
- **ε.** Supplying `FS_S7_NBHD` with a focus other than `effMaxSG` / `effMinSG` is an error (`:350`).
  - The check is on whether the variable is set, not on its value. Unset, the 0.10 default is still forwarded, and it equals `forestsearch()`'s own default. So every committed non-band render is unchanged.
  - **Note:** this is template policy, as the task specifies. The package itself treats ε as inert for those foci rather than rejecting it.
- **`selection_rule`.** A value other than `"neighborhood"` with a non-band focus is an error (`:555`). This mirrors what `.validate_selection_rule()` already requires (`R/subgroup_consistency_helpers.R:743`), but fails before any replicate runs.

### 1c. The recorder with MR off

The MR-off block sits at `:1105–1141`. It runs only when `FS_S7_MR=FALSE`; the MR-on path is textually unchanged.

| Column | MR off: read from | Why it is the same quantity |
|---|---|---|
| `n_sel` | `sum(fs.est$grp.consistency$sg.harm.id == 1L)` | MR's `n_selected` is the size of the family member whose membership equals `which(sg.harm.id == 1)` (`R/fs_mr_inference.R:582–590`, `:1046`). On all 123,892 detected committed rows (FS, DINA, GRF), `n_sel == n_harm`. |
| `label`, consistency | `fs.est$grp.consistency$out_sg$sg.harm`: the factor codes `names.Z[indexm == 1]` (`R/subgroup_consistency_helpers.R:913`, returned at `:1078`), joined with `" & "` | MR names its enumerated family `paste(colnames(Z)[covs.in == 1], collapse = " & ")` (`R/forestsearch_main.R:3367–3368`). |
| `label`, DINA / GRF | `fs.est$sg.harm` with the braces stripped | `.fs_mr_family_from_table()` spells the conjunction that way (`R/fs_mr_inference_methods.R:62–66`). It equals MR's label on all 66,676 committed DINA/GRF rows. |
| `n_cons_qual`, `band_n` | `fs.est$grp.consistency$out_sg$result` | It is the identifier's own table, the same read the MR-on block makes. It remains structurally NA on DINA and GRF, as before. |
| `n_family` | **left NA** | See below. |

- **`n_family` cannot be recovered without an `R/` change.**
  - It is the size of **MR's own fitted family**. That family is enumerated only inside the MR branch: `R/forestsearch_main.R:3352–3371` for consistency, `:2330` for DINA and `:2520` for GRF.
  - It is then filtered by per-candidate fits in `.fs_mr_assemble()` (`R/fs_mr_inference.R:80–97`), which drops candidates with fewer than 6 members or a failed fit.
  - No field of the result carries it. Rebuilding it in the template would re-implement MR.
- **MR-derived columns** (`mr_*`, `nv_*`, `fld_*`, `p_hat_*`, `mr_harm_flag`, `ij_source`) are NA when MR is off, by design. `mr_ok` is 0.

## Gate T3 — PASS, both halves

Standing identity cell: consistency, `effMaxSG` ε 0.20, HR 1.50, n 500, `FS_S7_Z1Q=0.60`, sim_id 1–5. Driver: `scripts_dinamr/t3gate.sh`; checker: `t3gate.R`; output: `logs/t3gate_T3.txt`.

- **Unset (`t3pre`, before the edit, vs `t3post`, after): PASS.**
  - 5 rows each, and the same 168 columns in the same order.
  - **All 163 non-timing columns `identical()`, and `truth` `identical()`.**
  - `meta$mr_inference` is absent before and `TRUE` after; it is record-only.
- **MR off (`t3post` vs `t3mroff`, `FS_S7_MR=FALSE`): PASS.**
  - Detected: 5/5 in both runs.
  - `n_sel`, `label`, `detected`, `sens`, `spec`, `ppv`, `npv`, `n_cons_qual` and `band_n` are present and populated on every detected row.
  - **The 30 identification and classification columns are `identical()` on/off**: `sim_id`, `detected`, `status`, `n_harm`, `n_true`, `label`, `sg_def`, `covs`, `err_msg`, `n_sel`, `n_cons_qual`, `band_n`, `admitted_n`, `sens`–`npv`, `betaHhat_*`, `nH_eval`, `nHc_eval`, `or_*`.
  - `truth` is `identical()`.
  - `n_family` is NA on every MR-off row. This is the documented exemption; with MR on it read 1233 / 1297 / 1212 / 1280 / 1303.
  - All 123 MR-derived columns are NA with MR off, and `mr_ok` is 0.
- **Supplementary, beyond the gate as specified: DINA and GRF on/off pairs at the same cell, 0 failures.** Their `label` route differs from consistency's, so it needed its own test. All 30 identification and classification columns are `identical()` on/off.

**A second, independent confirmation from the smoke.** At the smoke cell (12.4%, HR 1.50, n 500), the MR-off selections on sim_id 1–30 were compared with the **committed MR-on bundles** on the same seeds (`scripts_dinamr/partBoc_checks.R`):
- **DINA `effMaxSG` vs `dinamr` and GRF `effMaxSG` vs `grfmr`:** `detected`, `status`, `sg_def`, `n_sel`, `n_harm`, `n_true`, `label`, `sens`–`npv`, `betaHhat_*` and GRF's `admitted_n` are all `identical()`. `truth` is `identical()`.
- **FS `maxeffCons` vs `p12ext`** (built on pop-os, R 4.6.1): everything is `identical()` except `betaHhat_H` on one replicate and `betaHhat_Hc` on another, each differing by **1.1e-16**. `truth` agrees to 1e-12. That is cross-machine floating point; the selection is identical.

---

## Part 2 — the OC smoke

**Setup.**
- **Cell:** 12.4% (`FS_S7_Z1Q` unset), HR 1.50, n 500, sim_id 1–30.
- **Run settings:** MR off, FB none, 12 workers, thread variables at 1. Tag `pBoc`, and every stem carries `_nomr`.
- **Driver:** `scripts_dinamr/partBoc.sh`; table: `partBoc_table.R`, output `logs/partBoc_table.txt`.
- **Wall:** 16 of 16 runs completed with RC 0 and none were skipped, **318 s in total**.

### The OC table

**Column definitions:**
- **Rate:** detected / 30, with a Wilson 95% interval.
- **sens / spec / PPV / NPV:** replicate means over detected replicates.
- **|Ĥ| and |H|:** mean `n_sel` and mean `n_true` over detected replicates.
- **Prevalence:** mean `n_true / n` over all 30 replicates.
- **Per-rep s:** median `fit_mr_secs`. With MR off this times the identification alone.
- **Total s:** the render's `WALL_SECONDS`.

| Engine | `sg_focus` | ε | Rate [Wilson 95%] | Sens | Spec | PPV | NPV | Mean \|Ĥ\| | Mean \|H\| | Ratio | `n_family` median [range] | `admitted_n` median [range] | Prevalence | Per-rep s (median) | Total s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| consistency | `effMaxSG` | 0.20 | 0.967 [0.833, 0.994] | 0.540 | 0.842 | 0.342 | 0.928 | 102.7 | 62.3 | 1.649 | NA | — | 0.1237 | 1.59 | 21 |
| consistency | `effMinSG` | 0.20 | 0.967 [0.833, 0.994] | 0.496 | 0.921 | 0.479 | 0.929 | 65.9 | 62.3 | 1.059 | NA | — | 0.1237 | 1.50 | 19 |
| consistency | `maxeffCons` | inert | 0.967 [0.833, 0.994] | 0.521 | 0.907 | 0.463 | 0.929 | 72.9 | 62.3 | 1.170 | NA | — | 0.1237 | 1.51 | 20 |
| consistency | `maxeff` | inert | 1.000 [0.886, 1.000] | 0.533 | 0.907 | 0.463 | 0.931 | 73.5 | 61.8 | 1.189 | NA | — | 0.1237 | 2.64 | 23 |
| consistency | `maxSG` | inert | 0.967 [0.833, 0.994] | 0.632 | 0.725 | 0.283 | 0.937 | 159.7 | 62.3 | 2.565 | NA | — | 0.1237 | 1.55 | 19 |
| consistency | `minSG` | inert | 0.967 [0.833, 0.994] | 0.478 | 0.920 | 0.465 | 0.926 | 65.5 | 62.3 | 1.052 | NA | — | 0.1237 | 1.54 | 20 |
| dina | `effMaxSG` | 0.20 | 0.933 [0.787, 0.982] | 0.425 | 0.831 | 0.260 | 0.914 | 100.9 | 61.8 | 1.634 | NA | — | 0.1237 | 0.33 | 18 |
| dina | `effMinSG` | 0.20 | 0.933 [0.787, 0.982] | 0.310 | 0.901 | 0.311 | 0.903 | 62.8 | 61.8 | 1.017 | NA | — | 0.1237 | 0.31 | 18 |
| dina | `maxSG` | inert | 0.933 [0.787, 0.982] | 0.721 | 0.639 | 0.240 | 0.949 | 203.0 | 61.8 | 3.288 | NA | — | 0.1237 | 0.31 | 18 |
| dina | `minSG` | inert | 0.933 [0.787, 0.982] | 0.298 | 0.903 | 0.304 | 0.901 | 60.9 | 61.8 | 0.986 | NA | — | 0.1237 | 0.33 | 18 |
| dina | `eff` (stands for `eff` = `maxeff` = `maxeffCons`) | inert | 0.933 [0.787, 0.982] | 0.367 | 0.897 | 0.330 | 0.910 | 68.2 | 61.8 | 1.104 | NA | — | 0.1237 | 0.32 | 18 |
| grf | `effMaxSG` | 0.20 | 1.000 [0.886, 1.000] | 0.511 | 0.829 | 0.293 | 0.925 | 106.8 | 61.8 | 1.728 | NA | 124.5 [4, 387] | 0.1237 | 2.15 | 21 |
| grf | `effMinSG` | 0.20 | 1.000 [0.886, 1.000] | 0.329 | 0.901 | 0.315 | 0.905 | 64.1 | 61.8 | 1.037 | NA | 124.5 [4, 387] | 0.1237 | 2.12 | 21 |
| grf | `maxSG` | inert | 1.000 [0.886, 1.000] | 0.733 | 0.510 | 0.191 | 0.929 | 260.5 | 61.8 | 4.213 | NA | 124.5 [4, 387] | 0.1237 | 2.12 | 22 |
| grf | `minSG` | inert | 1.000 [0.886, 1.000] | 0.320 | 0.906 | 0.321 | 0.904 | 61.0 | 61.8 | 0.987 | NA | 124.5 [4, 387] | 0.1237 | 2.11 | 21 |
| grf | `eff` (stands for `eff` = `maxeff` = `maxeffCons`) | inert | 1.000 [0.886, 1.000] | 0.415 | 0.891 | 0.337 | 0.915 | 73.3 | 61.8 | 1.185 | NA | 124.5 [4, 387] | 0.1237 | 2.11 | 21 |

**Statuses.**
- Consistency: DETECTED 29, NO-DETECTION 1 for every criterion except `maxeff`, which detected 30.
- DINA: DETECTED 28, NO-DETECTION 2 for all five.
- GRF: DETECTED 30 for all five.
- No CONFIG-ERROR anywhere.
- `n_true` is identical across all 16 runs, which confirms they used the same DGM draws.

**How to read it** (30 replicates, one cell, descriptive only):
- **Whether a subgroup is found is set by the engine, not the criterion.**
  - Detection is flat across the criteria within each engine.
  - The one exception is consistency `maxeff`, on sim_id 21. `maxeff` disables the consistency and effect floors, and returned `{er <= 20} & !{meno}` (86 patients, sens 1.00, PPV 0.57) where every floored rule found nothing.
- **What is found moves with the criterion, as designed.**
  - Size runs `maxSG` > `effMaxSG` > `eff` / `maxeffCons` > `effMinSG` ≈ `minSG` on every engine. The mean |Ĥ|/|H| ratio is 2.6–4.2, 1.6–1.7, 1.1–1.2 and 1.0 respectively.
  - Sensitivity rises with size (0.30–0.50 → 0.63–0.73), while specificity and PPV fall (PPV 0.30–0.48 → 0.19–0.28).
  - The difference between the engines shows up at the ends: GRF's `maxSG` region is the largest at 4.2×, and DINA's `minSG` / `effMinSG` sensitivity is the lowest at about 0.30.
- **Realized prevalence** is 0.1237 on every run, against the 12.4% design.

### NA columns, and why

- **`n_family`: NA on all 16 runs.**
  - With MR off it cannot be recovered (Part 1c).
  - For reference, the committed MR-on bundles at this cell on sim_id 1–30 give `n_family` for their own criterion only:

    | Bundle (criterion) | `n_family` median [range] |
    |---|---|
    | `p12ext` (FS, `maxeffCons`) | 1230 [1112, 1333] |
    | `dinamr` (DINA, `effMaxSG`) | 212.5 [16, 1695] |
    | `grfmr` (GRF, `effMaxSG`) | 778 [729, 853] |

  - With MR off, `admitted_n` is GRF's family-size stratifier. Its value does not depend on the criterion: identical across GRF's five runs, median 124.5 [4, 387].
- **`n_cons_qual`, `band_n`: NA on every DINA and GRF run.** This is structural: neither engine has a consistency table, as before. On consistency both are populated.
- **`admitted_n`: NA on consistency and DINA**, by design. It is a GRF-path column.
- **MR-derived columns:** NA on all runs, by design.
- **Nothing else is NA on a detected row.**

### Same subgroup across criteria, same seeds

Cell [i, j] is the number of replicates, out of 30, where both criteria selected the identical rule; the diagonal is the detected count. Comparison is within engine only, because the rule strings are engine-specific.

| consistency | `effMaxSG` | `effMinSG` | `maxeffCons` | `maxeff` | `maxSG` | `minSG` |
|---|---|---|---|---|---|---|
| `effMaxSG` | 29 | 5 | 7 | 7 | 9 | 4 |
| `effMinSG` | 5 | 29 | 14 | 13 | 4 | 26 |
| `maxeffCons` | 7 | 14 | 29 | 27 | 5 | 12 |
| `maxeff` | 7 | 13 | 27 | 30 | 5 | 11 |
| `maxSG` | 9 | 4 | 5 | 5 | 29 | 4 |
| `minSG` | 4 | 26 | 12 | 11 | 4 | 29 |

| dina | `effMaxSG` | `effMinSG` | `maxSG` | `minSG` | `eff` |
|---|---|---|---|---|---|
| `effMaxSG` | 28 | 4 | 6 | 2 | 6 |
| `effMinSG` | 4 | 28 | 1 | 14 | 13 |
| `maxSG` | 6 | 1 | 28 | 1 | 1 |
| `minSG` | 2 | 14 | 1 | 28 | 6 |
| `eff` | 6 | 13 | 1 | 6 | 28 |

| grf | `effMaxSG` | `effMinSG` | `maxSG` | `minSG` | `eff` |
|---|---|---|---|---|---|
| `effMaxSG` | 30 | 0 | 1 | 0 | 2 |
| `effMinSG` | 0 | 30 | 0 | 6 | 6 |
| `maxSG` | 1 | 0 | 30 | 0 | 1 |
| `minSG` | 0 | 6 | 0 | 30 | 2 |
| `eff` | 2 | 6 | 1 | 2 | 30 |

- **No pair of criteria is identical at this cell.** Every pair differs on at least two jointly detected replicates.
- **Two consistency pairs are near-duplicates:**
  - **`maxeffCons` vs `maxeff`: identical on 27 of the 29 jointly detected replicates.** They differ on sim_id 9 and 15, where `maxeff`'s unfloored argmax picked a different rule. They also differ in detection on sim_id 21.
  - **`effMinSG` vs `minSG`: identical on 26 of 29.**
  - At this cell, the sweep will therefore show these two pairs differing mainly through a few replicates. That bears on whether the grid needs both members of each pair.
- **The same pairs are far apart on DINA and GRF.** `effMinSG` / `minSG` match on 14 of 28 (DINA) and 6 of 30 (GRF). The band rules and `eff` share a rule on at most 6 replicates on either engine.
- **GRF's five criteria almost never coincide.** No off-diagonal cell exceeds 6.

### Projected cost of the full 288 cell-runs

**Method** (`partBoc_table.R`):
- Compute = the sum, over the engine's criteria, of mean per-replicate seconds × replicates × 18 cells ÷ 12 workers.
- Render overhead: **16.4 s per render**, the median over the 16 runs of (`WALL_SECONDS` − compute).
- Renders per cell-run = ⌈replicates / 1,000⌉ batches plus one combine.

**Measured per-replicate means at this cell (s):**

| Engine | Criterion | Mean s per replicate |
|---|---|---|
| consistency | `effMaxSG` / `effMinSG` / `maxeffCons` / `maxSG` / `minSG` | 1.61 / 1.57 / 1.56 / 1.58 / 1.57 |
| consistency | **`maxeff`** | **2.55 (+63%)**: it evaluates every candidate with no truncation |
| DINA | all five | 0.52 |
| GRF | all five | 2.03–2.04 |

- **The flat-across-`sg_focus` assumption holds on DINA and GRF, and on FS for every criterion except `maxeff`.** The projection uses each criterion's own measured cost, so this is already accounted for.

**Flat (this one cell's cost applied to all 18 cells), wall hours at 12 workers:**

| Replicates per cell | consistency (108) | DINA (90) | GRF (90) | **All 288** |
|---|---|---|---|---|
| 2,000 | 10.18 (compute 8.70 + overhead 1.47) | 3.40 (2.17 + 1.23) | 9.72 (8.50 + 1.23) | **23.30** |
| 1,000 | 5.33 (4.35 + 0.98) | 1.90 (1.09 + 0.82) | 5.07 (4.25 + 0.82) | **12.30** |
| 500 | 3.16 (2.18 + 0.98) | 1.36 (0.54 + 0.82) | 2.94 (2.12 + 0.82) | **7.46** |

**Adjusted for how cost varies across the 18 cells (an approximation; `partBoc_checks.R`):**
- **Method.** Each engine's compute is multiplied by the ratio of its committed 18-cell mean to its value at this cell. The quantity used is the field-excluded bound from `REPORT_partB_measurement`.
- **Multipliers:** **consistency ×2.08, DINA ×2.33, GRF ×1.26.**
  - DINA's cost is about 3–4× higher at 31% than at 12.4%.
  - GRF's rises gently with n.
  - FS's profile comes mostly from pop-os bundles, so its multiplier is the least reliable.

| Replicates per cell | consistency | DINA | GRF | **All 288** |
|---|---|---|---|---|
| 2,000 | 19.56 | 6.28 | 11.96 | **37.80** |
| 1,000 | 10.02 | 3.36 | 6.18 | **19.56** |
| 500 | 5.51 | 2.08 | 3.50 | **11.09** |

**Caveats on both tables:**
- **The per-render overhead was measured on 30-row renders.** A 1,000-row render builds larger tables. The MR-on campaigns measured 76–197 s per render, but those include the MR tables, so the MR-off overhead at 1,000 rows is unmeasured. At 3–5 h of total overhead it is the smaller term in either table.
- **The flat table is the direct answer to the task's method; read the adjusted table as the more realistic figure.** The adjustment assumes MR-off cost varies across cells the way the committed field-excluded bound does.

### What MR cost at this cell, as a by-product

The smoke gives measured MR-off per-replicate seconds, and the committed MR-on bundles at the same cell and seeds give the MR-on seconds. Medians over sim_id 1–30, one cell, a small sample:

| Engine | MR off (s) | MR on (s) | MR share |
|---|---|---|---|
| DINA (`effMaxSG`, Mac against Mac) | 0.32 | 7.12 | **≈ 0.95** |
| GRF (`effMaxSG`, Mac against Mac) | 2.11 | 14.16 | **≈ 0.85** |
| FS (`maxeffCons`) | 1.51 | 45.17 (`p12ext`, pop-os) | **not comparable** across hosts |

- **Bundle size:** 16 MR-off bundles total 105 KB, about **0.22 KB per replicate** against 0.63–0.80 KB with MR on.

---

## Files

**Committed in `8fd89e1d`, the enabling change:**

| Path under `quarto/simulations/gbsg_020/` | Size |
|---|---|
| `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` (modified) | 162 KB |
| `scripts_dinamr/t3gate.sh`, `scripts_dinamr/t3gate.R` | 2 KB, 6 KB |
| `scripts_dinamr/logs/t3gate_*` (7 render logs, 7 driver logs, `t3gate_T3.txt`) | < 60 KB |
| `results/*_t3{pre,post,mroff,dina,dinaoff,grf,grfoff}_res_1_5.rds` (7) | 2.8–6.6 KB each |
| `t3gate_*.html` (7) | 3.9–4.5 MB each |

**Committed with this report:**

| Path under `quarto/simulations/gbsg_020/` | Size |
|---|---|
| `scripts_dinamr/partBoc.sh`, `partBoc_table.R`, `partBoc_checks.R` | 3 KB, 7 KB, 4 KB |
| `scripts_dinamr/partBoc_table.rds`, `partBoc_checks.rds` | < 5 KB each |
| `scripts_dinamr/logs/partBoc_*.log` (16 renders), `partBoc.driver.log`, `partBoc_table.txt`, `partBoc_checks.txt` | < 60 KB total |
| `results/*_nomr_pBoc_res_1_30.rds` (16) | 105 KB total |
| `partBoc_*.html` (16) | 60.3 MB total; largest 3.8 MB |
| `scripts_dinamr/README.md` (rows added), this report | — |

- **Nothing is over 50 MB.**

**Not done, by instruction:** no sweep cell, no replicate-count decision, and no `R/` change.
