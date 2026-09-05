# GATE 2 RECORD — MR (field) Section 7 Stage 2 (campaign s7), as amended

**Task:** `dev/tasks/TASK_mr_field_section7_2026-09-05.md` (bb516569). Amendment: Larry's Stage 2 STOP resolution (option 1; template knob `fb_join_skip`, commit 98f608c5). Prior records: Stage 2 STOP 9d39673d (diagnosis of the vintage flips).
**Date:** 2026-09-05. Run: template `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` at 100 workers (Stage 1 calibration), 14 h ceiling, fail-fast, campaign `s7`, forestsearch 0.3.5 (fix a6702fd8 installed). Wall: batches 10:59:23 → 11:34:23 (~8.0 / 8.0 / 9.5 / 9.5 min), combines +33 s — **35.5 min total**, no failures, no mid-run changes.

## GATE 2: ALL GREEN

**Amended wording applied:** naive / oracle / detection / MR (IJ) identical (≤ 1e-12) to the committed bundles on sim_id 1–1,000 *except on the enumerated flip replicates (h100: 393, 780), reported individually*.

### Completeness

| Batch | rows / sim_id | detections | CONFIG-ERROR |
|---|---|---|---|
| h100 res_1_1000 | 1000 / 1–1000 | 675 | 0 |
| h100 res_1001_2000 | 1000 / 1001–2000 | 686 | 0 |
| h175 res_1_1000 | 1000 / 1–1000 | 950 | 0 |
| h175 res_1001_2000 | 1000 / 1001–2000 | 950 | 0 |

Combined: h100 2000 rows, 1361 detections (68.0%); h175 2000 rows, 1900 detections (95.0%); sim_id 1–2000, no duplicates, both cells.

### Identity vs committed bundles (sims 1–1000)

- **h100:** all 25 gated columns exact-zero on the 998 non-flip replicates; status vectors identical. The strict scan beyond the gate (memberships n_sel/n_harm/n_true/sens/spec/ppv/npv and the β(Ĥ)/β(Ĥᶜ) truth attachments) is ALSO exact on all 998.
- **h175:** all 25 gated columns exact-zero on all 1,000; status identical; memberships exact on all 1,000.

**Enumerated flip replicates (h100), reported individually per the amendment:**

| sim | fresh (0.3.5) | committed (0.2.0) | oracle | detection | FB |
|---|---|---|---|---|---|
| 393 | `{pgr <= 0}`, n=64, naive 1.3421 | `{er <= 121} & {pgr <= 0}`, n=62, naive 1.3974 | exact (rel diff 0) | equal | absent (skipped) |
| 780 | `{pgr <= 0} & {er <= 17}`, n=62, naive 1.3826 | `{pgr <= 38} & {er <= 0}`, n=62, naive 1.4257 | exact (rel diff 0) | equal | absent (skipped) |

**Two additional findings surfaced by the strict scan, outside the gated column families (h175, β(Ĥ) truth attachment only):** sims **38** and **721** carry a different *cut value* in the realized rule string across vintages — sim 38: fresh `{age <= 53} & {pgr <= 7}` vs committed `{age <= 52} & {pgr <= 7}` (betaHhat_H rel diff 0.0407); sim 721: fresh `{er <= 3} & {age <= 52}` vs committed `{er <= 3} & {age <= 51}` (0.0123). Every trial-level column — naive, oracle, MR (IJ), memberships — is exact on both sims: the two cut spellings select the identical trial subgroup and differ only when the rule is evaluated on the 100k super-population by `fs_attach_betaHhat`. Same vintage-flip phenomenon as 393/780, at the cut-representation rather than the candidate-selection level. Not gate failures (the gate names naive/oracle/detection/MR (IJ)); carried to the Stage 3 flips section.

### FB join accounting (key cell)

- Joined rows with finite `fb_H_est`: **328** = committed FB detections (329) − 1 skipped-and-detected in the joined range (sim 393; sim 780 lies outside the committed FB range). ✔
- FB present only on sim_id ≤ 500 and never on a skip sim; all 8 `fb_*` columns exact vs the committed bundle on every joined row; sims 501–2,000 recorded FB-absent. ✔ The render prints both selections for each skipped sim (self-documenting).

### MR (field) completeness

- h100: 0 NA `est2` on 1,361 detections; 0 degenerate notes; min n_out 697/1000.
- h175: 0 NA `est2` on 1,900 detections; 0 degenerate notes; min n_out 766/1000.
- Interval invariant (lo2s ≤ hi2s) holds everywhere.

### Meta

Both combined bundles: n_sims 2000, 2 batches, `ci_method = field`, campaign `s7`, forestsearch 0.3.5; h100 `fb_mode` by batch join/none with `fb_join_skip` "393,780" recorded on the joining batch; h175 none/none.

Check script: session scratchpad `gate2_checks.R` (gated columns per the amended wording; strict scan reported, not gated). Driver: `stage2_driver.sh`; per-render logs `s2_*.log`.
