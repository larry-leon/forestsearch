# REPORT — GRF on the ACTG175 continuous (MD) design: campaign `mdgrf`, Stage 1 resumed after the membership fix (smoke, calibration → Gate 1)

Date: 2026-09-16 (UTC 2026-09-17). Machine: `pop-os` (64 physical cores, 251 GB; R 4.6.1, reference BLAS). Branch `feature/glm-extension`. Task: `dev/tasks/TASK_md_grf_resume_2026-09-16.md` (committed as received, `b62ae1aa`), executing `dev/tasks/TASK_md_grf_2026-09-16.md` ("the GRF task") from §1.6 with substitutions S1–S7. The first Stage 1 stopped at §1.6(c): `REPORT_md_grf_stage1_2026-09-16.md` (`a4c063bf`). The fix: `REPORT_grf_dina_fixes_2026-09-16.md` (P1 `0cd33f7b`). Every render ran with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`. The governing constraint stands: GRF's candidate family is generated from a fitted surface, so every coverage figure this campaign produces is coverage of the estimand conditional on the proposed family, and comparisons with FS are descriptive.

## S1 — preconditions — GATE PASS

```
pop-os
feature/glm-extension
019be60f
fix closeout 019be60f in HEAD
record: P1 landed          (REPORT_grf_dina_fixes_2026-09-16.md: "**P1 landed** (`0cd33f7b`)")
064fce91                   (the fix record's last R/ commit)
no R/ DESCRIPTION NAMESPACE change since 064fce91
Built: R 4.6.1; ; 2026-09-17 04:47:31 UTC; unix      (the fix record's final Built)
1788826 R 4.6.1; ; 2026-09-17 04:47:31 UTC; unix     (doFuture worker)
1788825 R 4.6.1; ; 2026-09-17 04:47:31 UTC; unix     (doFuture worker)
tracked mods: none
procs: none
```
First commit: `b62ae1aa Add TASK_md_grf_resume_2026-09-16 as received`.

## S2 — not repeated — GATE PASS

- Template: `894da993 MD template E1-E4 (TASK_md_grf_2026-09-16 §1.5, transplanted from the survival m1 template): …` — `git diff --quiet 894da993..HEAD -- sim_fs_maxeffCons_mr_field_md_template.qmd` succeeds.
- Scripts: `f0b9c844 scripts_mdgrf (TASK_md_grf_2026-09-16 §1.6, transplants of scripts_mdsgnb20): …` — `git diff --quiet f0b9c844..HEAD -- scripts_mdgrf/` succeeds.
- The GRF task's §1.3 (Gate 0), §1.4 (survival reference) and §1.5 (edits) stand as recorded in `REPORT_md_grf_stage1_2026-09-16.md` §1.3–§1.5.

## S3 — the stopped smoke, cleaned

Deleted (all untracked; `git ls-files` on them: 0), and nothing else: `mr_md_harm/fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdgrfsmokefs_d5000/`, `mr_md_harm/grf_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdgrfsmoke_d5000/`, and `logs_mdgrf/` `smoke.sh`, `smoke_driver.log`, `fs.log`, `grf.log`, `fs.peak_mb`, `grf.peak_mb`, `smoke_fs.html`, `smoke_grf.html`, `smoke_identity_fs.txt`, `smoke_identity_grf.txt` (the first record's list). The smoke driver was then re-created with the same content (finding R4).

## §1.6 Smoke — (a) PASS, (b) PASS, S5 gate PASS

Renders (`logs_mdgrf/smoke.sh`; md40 n500, sim_id 1–20, 20 workers, `FS_MD_FOCUS=effMaxSG FS_MD_NBHD=0.20 FS_MD_FIELD_SCALEC=selected FS_MD_CI=field FS_MD_FB=none`; template `894da993`; installed `Built 2026-09-17 04:47:31 UTC`):
```
SMOKE fs: rc=0 wall_s=91 peak_mb=18089 2026-09-17T04:53:04Z
SMOKE grf: rc=0 wall_s=81 peak_mb=15066 2026-09-17T04:54:25Z
SMOKE DONE
```

**(a) FS regression — S4** (`scripts_mdgrf/smoke_identity.R 40 500 mdgrfsmokefs 20 identity`):
```
== SMOKE identity: md=40 n=500 tag=mdgrfsmokefs | fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdgrfsmokefs_res_1_20.rds vs mdsgnb20 sim_id 1-20 ==
  [PASS] smoke bundle exists: mr_md_harm/fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdgrfsmokefs_d5000/fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdgrfsmokefs_res_1_20.rds
  [PASS] mdsgnb20 bundle exists: mr_md_harm/fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdsgnb20_d5000/fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdsgnb20_combined_1_2000.rds
  [PASS] both carry sim_id 1-20
  [PASS] no CONFIG-ERROR rows (status: DETECTED 20)
  new columns in this template, reported apart (16): n_family, n_cons_qual, band_n, admitted_n, p_hat_sum, fld_recov_sens_H, fld_recov_ppv_H, fld_recov_sens_Hc, fld_recov_npv_Hc, fld_recov_q10, fld_recov_q50, fld_recov_q90, fld_recov_share1, fld_recov_n_used, warn_msg, id_secs
  [PASS] all 137 paired mdsgnb20 columns present in the smoke bundle
  pairing: 20 of 20 rows identical on 126 numeric + 11 character columns (max rel diff among identical rows 0.00e+00)
  classification: no enumerated rows
  [PASS] zero selection flips (0 enumerated rows, all in a mdsgnb20 Gate 2 class)
  reported: mr_msg differs on 0 of 20 rows
  reported: err_msg differs on 0 of 20 rows
  reported: fb_err differs on 0 of 20 rows
  reported: FB columns (13) -- finite fb_H_est: mdsgnb20 0 rows, smoke 0
  reported (new columns on the FS path): n_cons_qual finite 20 | band_n finite 20 | n_family finite 20 | admitted_n finite 0 | p_hat_sum finite 20 | warn_msg non-NA 0 | id_secs finite 20
  [PASS] stem and meta carry the consistency engine (knob at its default)
  [PASS] nine fld_Hc_*_s and nine fld_joint_s_* columns present
  [PASS] nine fld_Hc_*_s finite on all 20 filled replicates
  [PASS] nine fld_joint_s_* finite on all 20 filled replicates
  [PASS] fld_Hc_lo1s_s <= fld_Hc_up1s_s and fld_Hc_lo2s_s <= fld_Hc_hi2s_s
  [PASS] Bonferroni harm bound identical between joint and joint_s where draw counts agree (20 of 20 rows agree; max |diff| 0.00e+00)
  [PASS] field-s inverted around the same beta-tilde^c (max |diff| 1.78e-15)
  [PASS] meta field_scale_complement = selected
SMOKE identity md=40 n=500: PASS
```
FS after the reinstall is bit-identical to `mdsgnb20` on 20 of 20 replicates (max relative difference 0). **PASS.**

**(b) GRF** (`… mdgrfsmoke 20 grf`):
```
== SMOKE grf: md=40 n=500 tag=mdgrfsmoke | grf_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdgrfsmoke_res_1_20.rds vs mdsgnb20 sim_id 1-20 ==
  [PASS] smoke bundle exists: mr_md_harm/grf_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdgrfsmoke_d5000/grf_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdgrfsmoke_res_1_20.rds
  [PASS] mdsgnb20 bundle exists: mr_md_harm/fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdsgnb20_d5000/fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdsgnb20_combined_1_2000.rds
  [PASS] both carry sim_id 1-20
  [PASS] no CONFIG-ERROR rows (status: DETECTED 20)
  new columns in this template, reported apart (16): n_family, n_cons_qual, band_n, admitted_n, p_hat_sum, fld_recov_sens_H, fld_recov_ppv_H, fld_recov_sens_Hc, fld_recov_npv_Hc, fld_recov_q10, fld_recov_q50, fld_recov_q90, fld_recov_share1, fld_recov_n_used, warn_msg, id_secs
  [PASS] n_true identical on every row
  [PASS] oracle columns (or_H_est,or_H_lo,or_H_hi,or_H_se,or_Hc_est,or_Hc_lo,or_Hc_hi,or_Hc_se) within 1e-8 relative on every row (max 0.00e+00)
  meta: subgroup_method grf | dmin_grf 30 | grf_selection frontier | grf_depth 2 | grf_select_statistic effect | focus effMaxSG | nbhd 0.2 | rule neighborhood | ci field | scalec selected | pkg 0.3.5 | host pop-os | workers 20 | R 4.6.1
  [PASS] meta carries identifier grf, dmin.grf 30, focus effMaxSG, band 0.20, the rule, ci_method field, field_scale_complement selected, pkg 0.3.5, host pop-os; stem grf_effMaxSG_
  FACT declared: 20 of 20 | NO-DETECTION 0 | admitted_n on non-detections: 
  FACT sim_id 1: sg_def [{preanti <= 792.80000000000018} & {cd40 > 364}] | n_sel 156 | n_harm 156 | nv_H_est 72.0396 | admitted_n 474 | n_family 1257 | sens 0.346 ppv 0.404
  FACT S0.3 (quoted): {preanti <= 792.8} & {cd40 > 364}, n 156 / 344, oriented MD 72.04, admitted_n 234, 1,257 enumerated, sens 0.346 ppv 0.404
  [PASS] admitted_n filled and >= 1 on every declared replicate (min 256 max 1172)
  [PASS] admitted_n recorded on every replicate (0L where the admitted set was empty)
  [PASS] n_family filled on every declared replicate with a gate (values 1185-1372)
  [PASS] p_hat_H and p_hat_sum recorded (p_hat_H mean 0.084)
  STRUCTURAL n_cons_qual all-NA: TRUE | band_n all-NA: TRUE (GRF has no consistency screen)
  FACT the floor as applied: meta dmin_grf = 30 on the DR-score harm effect; effect_threshold = 30 on the harm-oriented MD (admission)
  WARNINGS (verbatim, distinct, with per-row counts on 0 of 20 rows):
  FACT timing: fit_mr_secs mean 32.5 median 32.5 max 35.6 | id_secs (GRF fit, incl. re-selection) mean 4.04 median 4.04 max 4.32 | fld_H_secs mean 18.9 | fld_Hc_secs mean 1.26
  FACT vs mdsgnb20 (sim 1-20): mean |Hhat| (n_harm) 119.7 vs 120.3 | sens 0.301 vs 0.309 | ppv 0.432 vs 0.423 | declared 20 vs 20
  [PASS] nine fld_Hc_*_s and nine fld_joint_s_* columns present
  [PASS] nine fld_Hc_*_s finite on all 20 filled replicates
  [PASS] nine fld_joint_s_* finite on all 20 filled replicates
  [PASS] fld_Hc_lo1s_s <= fld_Hc_up1s_s and fld_Hc_lo2s_s <= fld_Hc_hi2s_s
  [PASS] Bonferroni harm bound identical between joint and joint_s where draw counts agree (20 of 20 rows agree; max |diff| 0.00e+00)
  [PASS] field-s inverted around the same beta-tilde^c (max |diff| 7.11e-15)
  [PASS] meta field_scale_complement = selected
SMOKE grf md=40 n=500: PASS
```
**PASS.** Facts (S5): sim_id 1 selects `{preanti <= 792.8} & {cd40 > 364}`, n 156 — the same as Stage 0 S0.3 (`{preanti <= 792.8} & {cd40 > 364}`, n 156 / 344, admitted 234, 1,257 enumerated) and as the fix record's F4 after P1 (the same selection, n 156, admitted 474, 1,257 enumerated, 0 NA). Here `admitted_n` is 474 and MR's kept family `n_family` is 1,257 — the full enumerated pool (the first Stage 1 recorded 612, finding F3 there). Declared 20 of 20; `fit_mr_secs` mean 32.5 / median 32.5 / max 35.6 s; GRF fit (`id_secs`) mean 4.04 / median 4.04 / max 4.32 s.

**S5 gate — factor-comparison warnings and NA memberships.**
- Captured-warnings column: `warn_msg` is NA on 20 of 20 GRF rows ("WARNINGS … on 0 of 20 rows") — zero factor-comparison warnings.
- NA-membership candidates: the recorder has no such column and the template is frozen (S2), so each replicate was regenerated with the template's own chunks and argument block and re-identified with MR off (`logs_mdgrf/na_membership_check.R`, untracked), counting `.grf_evaluate_subgroup()` NA memberships over every enumerated candidate and checking the selection and `admitted_n` against the smoke bundle (finding R3):
```
forestsearch 0.3.5 Built R 4.6.1; ; 2026-09-17 04:47:31 UTC; unix 
sim  1: enumerated 1257 | NA membership 0 | factor warnings 0 | admitted 474 (bundle 474) | selection matches bundle: TRUE
sim  2: enumerated 1372 | NA membership 0 | factor warnings 0 | admitted 481 (bundle 481) | selection matches bundle: TRUE
sim  3: enumerated 1274 | NA membership 0 | factor warnings 0 | admitted 900 (bundle 900) | selection matches bundle: TRUE
sim  4: enumerated 1185 | NA membership 0 | factor warnings 0 | admitted 277 (bundle 277) | selection matches bundle: TRUE
sim  5: enumerated 1187 | NA membership 0 | factor warnings 0 | admitted 348 (bundle 348) | selection matches bundle: TRUE
sim  6: enumerated 1193 | NA membership 0 | factor warnings 0 | admitted 560 (bundle 560) | selection matches bundle: TRUE
sim  7: enumerated 1260 | NA membership 0 | factor warnings 0 | admitted 472 (bundle 472) | selection matches bundle: TRUE
sim  8: enumerated 1188 | NA membership 0 | factor warnings 0 | admitted 737 (bundle 737) | selection matches bundle: TRUE
sim  9: enumerated 1189 | NA membership 0 | factor warnings 0 | admitted 999 (bundle 999) | selection matches bundle: TRUE
sim 10: enumerated 1198 | NA membership 0 | factor warnings 0 | admitted 975 (bundle 975) | selection matches bundle: TRUE
sim 11: enumerated 1240 | NA membership 0 | factor warnings 0 | admitted 1172 (bundle 1172) | selection matches bundle: TRUE
sim 12: enumerated 1249 | NA membership 0 | factor warnings 0 | admitted 586 (bundle 586) | selection matches bundle: TRUE
sim 13: enumerated 1188 | NA membership 0 | factor warnings 0 | admitted 256 (bundle 256) | selection matches bundle: TRUE
sim 14: enumerated 1276 | NA membership 0 | factor warnings 0 | admitted 412 (bundle 412) | selection matches bundle: TRUE
sim 15: enumerated 1202 | NA membership 0 | factor warnings 0 | admitted 487 (bundle 487) | selection matches bundle: TRUE
sim 16: enumerated 1218 | NA membership 0 | factor warnings 0 | admitted 275 (bundle 275) | selection matches bundle: TRUE
sim 17: enumerated 1202 | NA membership 0 | factor warnings 0 | admitted 305 (bundle 305) | selection matches bundle: TRUE
sim 18: enumerated 1211 | NA membership 0 | factor warnings 0 | admitted 1071 (bundle 1071) | selection matches bundle: TRUE
sim 19: enumerated 1191 | NA membership 0 | factor warnings 0 | admitted 846 (bundle 846) | selection matches bundle: TRUE
sim 20: enumerated 1196 | NA membership 0 | factor warnings 0 | admitted 522 (bundle 522) | selection matches bundle: TRUE
TOTAL: NA-membership candidates 0 | factor-comparison warnings 0 | selections matching the bundle 20 of 20
```
**GATE PASS**: zero warnings and zero NA-membership candidates on every replicate; 20 of 20 selections and admitted counts match the bundle.

*GATE §1.6:* (a), (b) and the S5 replacement of (c) hold.

## §1.7 Calibration

md40 n700, knob `grf`, `FS_MD_DMIN_GRF=30`, the campaign knobs, `FS_MD_WORKERS` = 16, 32, 63 with 3 × W replicates (tags `mdgrfcal16/32/63`), memory sampled every 5 s (`logs_mdgrf/calib.sh`; summary `logs_mdgrf/calib_summary.R`, the `mdsgnb20` method of `REPORT_md_field_rerun_stage1_2026-09-15.md` §1.7):
```
CALIB W=16 reps=48 rc=0 wall_s=157 peak_mb=15169 2026-09-17T04:59:03Z
CALIB W=32 reps=96 rc=0 wall_s=182 peak_mb=29036 2026-09-17T05:02:05Z
CALIB W=63 reps=189 rc=0 wall_s=264 peak_mb=56580 2026-09-17T05:06:29Z
CALIB DONE
total wall 603 s
```
| W | replicates | render wall (s) | fit_mr_secs mean / median / p90 / max (s) | ratio to 16-worker mean | GRF fit (id_secs) mean / median / p90 (s) | fld_H_secs / fld_Hc_secs mean (s) | peak summed RSS (MB) | replicates per minute | declared | rows with warnings |
|---|---|---|---|---|---|---|---|---|---|---|
| 16 | 48 | 157 | 37.47 / 37.52 / 39.40 / 40.15 | 1.000 | 4.32 / 4.30 / 4.58 | 20.46 / 1.73 | 15169 | 18.34 | 48 | 0 |
| 32 | 96 | 182 | 38.66 / 39.06 / 40.84 / 44.31 | 1.032 | 4.51 / 4.50 / 4.94 | 21.14 / 1.77 | 29036 | 31.65 | 96 | 0 |
| 63 | 189 | 264 | 46.73 / 47.00 / 54.87 / 60.93 | 1.247 | 5.93 / 6.00 / 7.50 | 25.38 / 2.10 | 56580 | 42.95 | 189 | 0 |

Fixed render overhead (16-worker render, wall minus 3 rounds x mean): 44.6 s. Wall-based loop cost per replicate: W 16: 2.342 s; W 32: 1.431 s; W 63: 1.161 s.

Projection at W = 63 (loop cost 1.161 s per replicate at md40 n700, scaled per cell by mdsgnb20 fit_mr_secs ratios md40_n500 0.6626, md120_n500 0.8017, null_n500 0.6352, md40_n700 1.0000; overhead 45 s per batch render; combine render 90 s):

| cell | mdsgnb20 fit_mr_secs mean (s) | ratio | batch of 1,000 (s) | batch (min) | cell: 2 batches + combine (min) |
|---|---|---|---|---|---|
| md40_n500 | 56.24 | 0.6626 | 814 | 13.6 | 28.6 |
| md120_n500 | 68.04 | 0.8017 | 975 | 16.3 | 34.0 |
| null_n500 | 53.91 | 0.6352 | 782 | 13.0 | 27.6 |
| md40_n700 | 84.88 | 1.0000 | 1206 | 20.1 | 41.7 |

**Projection: 7913 s = 131.9 min = 2.20 h.  Ceiling (1.5x): 11870 s = 197.8 min = 3.30 h.  Per-render timeout (2 x the longest projected batch, at least 20 min): 2411 s = 40.2 min.**
MDSG_WORKERS=63 MDSG_TIMEOUT=2411 MDSG_CEILING=11870
ADVANCE_GO_CONDITION projection_h=2.198 under_8h=TRUE

- **W = 63.** It minimizes projected wall: 42.95 replicates per minute against 31.65 at 32 and 18.34 at 16. **What limits it:** the template's worker cap, physical cores − 1 = 63 (`sim_fs_maxeffCons_mr_field_md_template.qmd:115–116`), i.e. cores; scaling falloff is visible (per-replicate `fit_mr_secs` 1.25× the 16-worker mean) but throughput still rises to the cap; memory is not binding (peak summed RSS 56.6 GB of 251 GB).
- **Projection for Stage 2 at W = 63: 7,913 s (131.9 min, 2.20 h)**; ceiling 1.5 × projection = **11,870 s (197.8 min, 3.30 h)**; per-render timeout 2 × the longest projected batch (20.1 min) = **2,411 s (40.2 min)**, above the 20-min floor. Runner environment: `MDSG_WORKERS=63 MDSG_TIMEOUT=2411 MDSG_CEILING=11870`.
- Every calibration replicate declared (48, 96, 189) and none carries a warning.

## §1.8 Gate 1 — PASS; the advance go applies

Every Stage 1 gate is green (S1, S2, §1.6(a), §1.6(b), the S5 gate) and the §1.7 projection (2.20 h) is under 8 hours. **The advance go, as the GRF task states it (Dispositions):** "if every Stage 1 gate is green and the §1.7 projection for Stage 2 is under 8 hours, do not stop at Gate 1. Note the advance go and its condition in the Gate 1 record, run Stage 2 at the worker count chosen in §1.7 with the record's ceiling and per-render timeout, and run Stage 3 once Stage 2 is green." Its condition is met. Stage 2 launches after this record's commit with `MDSG_WORKERS=63 MDSG_TIMEOUT=2411 MDSG_CEILING=11870`, knob `grf`, `FS_MD_DMIN_GRF=30`, the rule knobs, cells in `mdsgnb20`'s order; no reduction of cells, replicates, gates or knobs.

## Findings

- **R1.** After P1, MR's kept family on the GRF path is the full enumerated pool (`n_family` 1,185–1,372 on the smoke; 1,257 on sim_id 1), and GRF's admitted set is larger (`admitted_n` 256–1,172; sim_id 1: 474 against 234 before the fix). The first Stage 1's F3 (family = scorable continuous-only candidates) no longer holds.
- **R2.** The GRF path costs more after the fix: `fit_mr_secs` 32.5 s per replicate at 20 workers against 22.8 s in the first smoke; the field pass 18.9 s against 13.8 s — the family doubled. The calibration and projection use the post-fix cost.
- **R3.** The S5 NA-membership gate is computed by a read-only re-identification (MR off), because the recorder has no such column and S2 freezes the template; its selections and `admitted_n` match the smoke bundle on 20 of 20, which ties the count to the recorded replicates.
- **R4.** S3 deleted the first smoke's driver with its outputs; `logs_mdgrf/smoke.sh` was re-created with the same knob set and render lines.
- **R5.** Stage 1 compute in this run: smoke 172 s of renders, checks, calibration 603 s — under the GRF task's 2-h ceiling.
- **R6.** The projection allows 90 s per combine render (the `mdsgnb20` allowance); `mdsgnb20`'s combines took 40–41 s.
- **R7.** Stage 3 preparation, done in the waiting time: `summary_continuous_field_mdgrf.qmd` (the §3.1 transplant, untracked until Stage 3) dry-rendered against the md40 n500 GRF smoke bundle (`MDSG_SUMMARY_TAG=mdgrfsmoke MDSG_SUMMARY_GLOB=res_1_20`, exit 0; 277 GRF extract rows and 4 FS comparator rows) under `logs_mdgrf/dryrun/` (untracked).

## Untracked outputs of this stage (for §3.5 deletion)

`mr_md_harm/fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdgrfsmokefs_d5000/`, `mr_md_harm/grf_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdgrfsmoke_d5000/`, `mr_md_harm/grf_effMaxSG_mr_field_md40_knoise0_n700_nb20_mdgrfcal{16,32,63}_d5000/`, and in `logs_mdgrf/`: `smoke.sh`, `smoke_driver.log`, `fs.log`, `grf.log`, `fs.peak_mb`, `grf.peak_mb`, `smoke_fs.html`, `smoke_grf.html`, `smoke_identity_fs.txt`, `smoke_identity_grf.txt`, `na_membership_check.R`, `na_membership_check.txt`, `calib.sh`, `calib_driver.log`, `calib_summary.R`, `calib_summary.txt`, `mdgrfcal{16,32,63}.{log,peak_mb,html}`, `dryrun/`.

## git log --oneline for this stage

```
b62ae1aa Add TASK_md_grf_resume_2026-09-16 as received
<this record: the next commit>
```
