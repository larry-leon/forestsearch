# REPORT — GRF on the ACTG175 continuous (MD) design: campaign `mdgrf`

Date: 2026-09-17 (UTC). Machine: `pop-os` (AMD Ryzen Threadripper PRO 5995WX, 64 physical cores, 251 GB; R 4.6.1, reference BLAS). Branch `feature/glm-extension`. Tasks: `dev/tasks/TASK_md_grf_2026-09-16.md` ("the GRF task", `e5ee1008`) and `dev/tasks/TASK_md_grf_resume_2026-09-16.md` (`b62ae1aa`), which executes the GRF task from §1.6 after the membership fix. Records of this campaign: the stopped Stage 1 `REPORT_md_grf_stage1_2026-09-16.md` (`a4c063bf`), the fix `REPORT_grf_dina_fixes_2026-09-16.md` (`baaf0bc6`), the resumed Stage 1 `REPORT_md_grf_stage1_resume_2026-09-16.md` (`f5cc256d`), Gate 2 `REPORT_md_grf_gate2_2026-09-16.md` (per cell, committed by the runner), this record. Installed forestsearch 0.3.5, `Built: R 4.6.1; ; 2026-09-17 04:47:31 UTC; unix` (the fix task's install); no `R/` change in the campaign. Every render ran with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`.

GRF's candidate family is generated from a fitted surface, so the fixed-family condition does not hold: every coverage figure below is coverage of the estimand conditional on the proposed family. Comparisons with FS are descriptive, not a contest: the identifier and the family construction differ, and each summary conditions on a different set of detected replicates.

## Dispositions (Larry, 2026-09-16), as executed

- **Identifier:** GRF only (`FS_MD_METHOD=grf`: `subgroup_method = "grf"`, `grf_selection = "frontier"`, `grf_select_statistic = "effect"`, `grf_depth = 2`). DINA was not run.
- **Rule:** `effMaxSG`, `effect_neighborhood = 0.20`, `selection_rule = "neighborhood"` (`mdsgnb20`'s `meta`), through `FS_MD_FOCUS` / `FS_MD_NBHD`.
- **GRF admission floor:** `dmin.grf = 30` (`FS_MD_DMIN_GRF`), recorded in every bundle's `meta` (`dmin_grf 30`); it floors the DR-score pre-filter, while the harm-oriented MD admission floor 30 is `effect.threshold`'s (first Stage 1 record, F1).
- **Other GRF arguments:** the survival `grfmr` campaign's (`frontier`, `effect`, depth 2), with `frac.tau` omitted as survival-only (first Stage 1 record §1.4).
- **MR:** `ci_method = "field"`, `draws = 5000`, `include_complement`, `field_complement`, `field_scale_complement = "selected"`, `return_reselection = TRUE`, `ij_residual = "two_term"`, `confirm_rule = "point"`, `t_confirm` near-null — the MD template's.
- **Re-selection alignment:** aligned (Stage 0 S0.4; first Stage 1 record §1.3).
- **Labelling:** the conditional-family sentence opens this record and the summary, and ends every summary table and figure caption.
- **Gate 1 advance go:** its condition (every Stage 1 gate green, projection under 8 hours) held — projection 2.20 h — so Stages 2 and 3 ran (resumed Stage 1 record §1.8).
- **Resume substitutions S1–S7** (the resume task): preconditions, no repeat of §1.3–§1.5, the stopped smoke cleaned, FS regression re-run, the (c) gate replaced by zero factor warnings and zero NA memberships, the resume record name, and the DINA open-work line worded as the fix record words P2.

## Gate 1 and Gate 2 in brief

- **Gate 1.** First Stage 1 (`a4c063bf`): template edits `894da993`, scripts `f0b9c844`, FS regression bit-identical, GRF smoke green, **stop at §1.6(c)** on the factor-membership defect (645 of 1,257 candidates dropped on sim_id 1). The fix P1 landed (`0cd33f7b`). Resumed Stage 1 (`f5cc256d`): FS regression bit-identical to `mdsgnb20` on 20 of 20 again; GRF smoke green with **zero factor-comparison warnings and zero NA-membership candidates on 20 of 20** (selections and admitted counts re-derived and matched); sim_id 1 selects `{preanti <= 792.8} & {cd40 > 364}`, n 156, admitted 474, family 1,257; calibration W = 63 (42.95 replicates per minute; peak 56.6 GB), projection 7,913 s (2.20 h), ceiling 11,870 s, timeout 2,411 s.
- **Gate 2** (`REPORT_md_grf_gate2_2026-09-16.md`, `scripts_mdgrf/gate2.R`): every cell passed **66 of 66** checks — 2,000 rows with `sim_id` exactly 1–2000 and the batch files matching the combined bundle on all 172 columns; `meta` carrying `grf`, `dmin_grf 30`, `frontier` / depth 2 / `effect`, `effMaxSG` / 0.20 / `neighborhood`, `ci_method field`, `field_scale_complement selected`, `pkg_version 0.3.5`, `hostname pop-os`, 63 workers; the same draws as `mdsgnb20` in both directions (`n_true` identical on 2,000 rows; oracle columns equal, max relative difference 0; complement only in the null cell); field-s and recorder checks on every declared replicate; `admitted_n` ≥ 1 on every declared replicate; **0 MR failures** and **0 rows with a captured warning** in every cell. Declared: 2,000 of 2,000 in every cell. Cell commits `c4c49572`, `f84261e3`, `15bee7cd`, `432f69dc`; campaign complete `1a65b28a`. No halt.
- **Render walls** (the progress log, UTC):
```
2026-09-17T05:07:47Z	campaign	start	HEAD=f5cc256d workers=63 timeout_s=2411 ceiling_s=11870 built=R 4.6.1; ; 2026-09-17 04:47:31 UTC; unix
2026-09-17T05:07:47Z	md40_n500	start	stem=grf_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdgrf
2026-09-17T05:21:58Z	md40_n500	batch_1_1000	wall_s=851 peak_mb=57326 rc=0 cumulative_s=851
2026-09-17T05:36:14Z	md40_n500	batch_1001_2000	wall_s=856 peak_mb=57837 rc=0 cumulative_s=1707
2026-09-17T05:36:55Z	md40_n500	combine_1_2000	wall_s=41 peak_mb=1178 rc=0 cumulative_s=1748
2026-09-17T05:36:55Z	md40_n500	done	cell_wall_s=1748 GATE_COUNTS run=66 passed=66 failed=0
2026-09-17T05:36:55Z	md120_n500	start	stem=grf_effMaxSG_mr_field_md120_knoise0_n500_nb20_mdgrf
2026-09-17T05:52:39Z	md120_n500	batch_1_1000	wall_s=944 peak_mb=57554 rc=0 cumulative_s=2692
2026-09-17T06:08:33Z	md120_n500	batch_1001_2000	wall_s=954 peak_mb=57621 rc=0 cumulative_s=3646
2026-09-17T06:09:09Z	md120_n500	combine_1_2000	wall_s=36 peak_mb=1155 rc=0 cumulative_s=3682
2026-09-17T06:09:09Z	md120_n500	done	cell_wall_s=1934 GATE_COUNTS run=66 passed=66 failed=0
2026-09-17T06:09:09Z	null_n500	start	stem=grf_effMaxSG_mr_field_mdnull_knoise0_n500_nb20_mdgrf
2026-09-17T06:23:11Z	null_n500	batch_1_1000	wall_s=842 peak_mb=57536 rc=0 cumulative_s=4524
2026-09-17T06:37:07Z	null_n500	batch_1001_2000	wall_s=836 peak_mb=58198 rc=0 cumulative_s=5360
2026-09-17T06:37:47Z	null_n500	combine_1_2000	wall_s=40 peak_mb=1109 rc=0 cumulative_s=5400
2026-09-17T06:37:47Z	null_n500	done	cell_wall_s=1718 GATE_COUNTS run=66 passed=66 failed=0
2026-09-17T06:37:48Z	md40_n700	start	stem=grf_effMaxSG_mr_field_md40_knoise0_n700_nb20_mdgrf
2026-09-17T06:54:22Z	md40_n700	batch_1_1000	wall_s=994 peak_mb=65133 rc=0 cumulative_s=6394
2026-09-17T07:10:51Z	md40_n700	batch_1001_2000	wall_s=989 peak_mb=64611 rc=0 cumulative_s=7383
2026-09-17T07:11:31Z	md40_n700	combine_1_2000	wall_s=40 peak_mb=1184 rc=0 cumulative_s=7423
2026-09-17T07:11:31Z	md40_n700	done	cell_wall_s=2023 GATE_COUNTS run=66 passed=66 failed=0
2026-09-17T07:11:32Z	campaign	complete	cumulative_s=7423
```
  **7,423 s (123.7 min) in total**, against the 7,913-s projection and the 11,870-s ceiling; peak summed RSS 57.3–58.2 GB at n = 500 and 64.6–65.1 GB at n = 700 (63 workers); `fit_mr_secs` mean 45.1 / 51.2 / 44.1 / 53.5 s (md40 n500 / md120 / null / md40 n700), of which the field pass 27.2 / 32.3 / 26.3 / 29.3 s and the GRF fit 6.0 / 6.5 / 6.0 / 7.2 s.

## Tables (pasted from the render of `summary_continuous_field_mdgrf.qmd`, `0333248d`; every number is in `md_grf_metrics.csv`, `0333248d`)

Conventions: harm-oriented MD scale (positive = harm); targets β(Ĥ), β(Ĥᶜ) exact per replicate, the oracle against the structural true-region effect; one-sided coverage on the exposed side (Ĥ: LOWER bound; Ĥᶜ: UPPER bound); Wilson 95% intervals in parentheses; Monte Carlo SEs in parentheses in the ladder tables. The null cell is labelled by its truth (no subgroup; homogeneous +26); its oracle row on Ĥ is blank (Q is empty). **Every table is conditional on the proposed family.**

### Declaration (conditional on the proposed family)

|cell|replicates|declared|rate|mc_se|fs_rate|
|---|---|---|---|---|---|
|md40 n500|2000|2000|1|0|0.9990|
|md120 n500|2000|2000|1|0|1.0000|
|null n500 (no subgroup; homogeneous +26)|2000|2000|1|0|0.9965|
|md40 n700|2000|2000|1|0|0.9995|

- GRF declares a subgroup on 2,000 of 2,000 replicates in every cell, the null included; FS (`mdsgnb20`) declares on 0.9990 / 1.0000 / 0.9965 / 0.9995 on the same seeds.
- In the null cell every declared Ĥ has β(Ĥ) = +26.26, so its rows below are read as coverage of a homogeneous harm, not as a false-claim rate.

### Ĥ — unadjusted, oracle, IJ two-term, field (conditional on the proposed family)

|cell|estimator|n|declaration rate|bias (MD)|bias (SD units)|SE/SD|one-sided LOWER coverage (Wilson)|two-sided coverage (Wilson)|
|---|---|---|---|---|---|---|---|---|
|md40 n500|naive|2000|1|58.624|3.293|1.450|0.158 (0.143, 0.175)|0.320 (0.300, 0.341)|
|md40 n500|oracle|2000|1|-0.684|-0.034|0.991|0.953 (0.942, 0.961)|0.943 (0.932, 0.952)|
|md40 n500|MR (IJ)|2000|1|15.422|0.745|1.628|0.969 (0.960, 0.975)|0.992 (0.986, 0.995)|
|md40 n500|MR (field)|2000|1|7.400|0.333|1.133|0.939 (0.928, 0.949)|0.965 (0.955, 0.972)|
|md120 n500|naive|2000|1|32.906|1.682|1.165|0.599 (0.577, 0.620)|0.713 (0.693, 0.732)|
|md120 n500|oracle|2000|1|-0.684|-0.034|0.991|0.953 (0.942, 0.961)|0.943 (0.932, 0.952)|
|md120 n500|MR (IJ)|2000|1|-0.455|-0.019|1.382|0.989 (0.983, 0.992)|0.983 (0.976, 0.988)|
|md120 n500|MR (field)|2000|1|-4.088|-0.161|1.004|0.956 (0.946, 0.964)|0.926 (0.914, 0.937)|
|null n500 (no subgroup; homogeneous +26)|naive|2000|1|60.180|3.383|1.473|0.143 (0.128, 0.159)|0.300 (0.280, 0.320)|
|null n500 (no subgroup; homogeneous +26)|oracle|0|1|NaN|NaN|NA|NaN (NA, NA)|NaN (NA, NA)|
|null n500 (no subgroup; homogeneous +26)|MR (IJ)|2000|1|16.271|0.784|1.639|0.972 (0.963, 0.978)|0.991 (0.986, 0.994)|
|null n500 (no subgroup; homogeneous +26)|MR (field)|2000|1|7.966|0.357|1.136|0.941 (0.929, 0.950)|0.969 (0.960, 0.975)|
|md40 n700|naive|2000|1|60.124|3.450|1.447|0.100 (0.087, 0.113)|0.235 (0.217, 0.255)|
|md40 n700|oracle|2000|1|0.060|0.003|0.967|0.941 (0.930, 0.951)|0.943 (0.932, 0.952)|
|md40 n700|MR (IJ)|2000|1|15.787|0.790|1.587|0.967 (0.958, 0.974)|0.989 (0.983, 0.992)|
|md40 n700|MR (field)|2000|1|7.596|0.355|1.149|0.943 (0.932, 0.952)|0.967 (0.958, 0.974)|

- The field's one-sided lower bound on β(Ĥ) covers at 0.939 / 0.956 / 0.941 / 0.943 (md40 n500 / md120 n500 / null / md40 n700); the Wilson interval contains 0.95 at md120 (0.946–0.964), the null (0.929–0.950) and n = 700 (0.932–0.952), and lies below it at md40 n500 (0.928–0.949). The oracle's is 0.941–0.953.
- Retained bias in SD units: unadjusted +3.293 / +1.682 / +3.383 / +3.450, IJ +0.745 / −0.019 / +0.784 / +0.790, field +0.333 / −0.161 / +0.357 / +0.355. The field's SE/SD is 1.004–1.149; IJ's is 1.382–1.639, with one-sided coverage 0.967–0.989.
- The field's two-sided interval covers at 0.965 / 0.926 / 0.969 / 0.967; at md120 it is below nominal while the one-sided lower bound is above it.

### Ĥᶜ — unadjusted, oracle, IJ two-term, field (beside), field-s (evaluated) (conditional on the proposed family)

|cell|estimator|n|declaration rate|bias (MD)|bias (SD units)|SE/SD|one-sided UPPER coverage (Wilson)|two-sided coverage (Wilson)|
|---|---|---|---|---|---|---|---|---|
|md40 n500|naive|2000|1|-16.315|-1.334|1.078|0.677 (0.656, 0.697)|0.799 (0.781, 0.816)|
|md40 n500|oracle|2000|1|-0.301|-0.021|1.002|0.956 (0.946, 0.964)|0.953 (0.942, 0.961)|
|md40 n500|MR (IJ)|2000|1|-5.810|-0.461|1.885|0.995 (0.990, 0.997)|0.999 (0.996, 1.000)|
|md40 n500|MR (field)|2000|1|-3.666|-0.286|1.006|0.916 (0.903, 0.927)|0.942 (0.930, 0.951)|
|md40 n500|MR (field-s)|2000|1|-3.677|-0.287|1.011|0.920 (0.907, 0.931)|0.938 (0.926, 0.947)|
|md120 n500|naive|2000|1|-12.443|-0.868|0.980|0.790 (0.772, 0.807)|0.866 (0.850, 0.880)|
|md120 n500|oracle|2000|1|-0.301|-0.021|1.002|0.956 (0.946, 0.964)|0.953 (0.942, 0.961)|
|md120 n500|MR (IJ)|2000|1|-2.768|-0.185|1.671|0.997 (0.993, 0.999)|1.000 (0.997, 1.000)|
|md120 n500|MR (field)|2000|1|-1.212|-0.079|0.890|0.926 (0.914, 0.937)|0.929 (0.917, 0.940)|
|md120 n500|MR (field-s)|2000|1|-1.228|-0.080|0.914|0.931 (0.920, 0.942)|0.934 (0.923, 0.945)|
|null n500 (no subgroup; homogeneous +26)|naive|2000|1|-16.094|-1.329|1.084|0.688 (0.667, 0.708)|0.798 (0.780, 0.815)|
|null n500 (no subgroup; homogeneous +26)|oracle|2000|1|-0.431|-0.037|1.012|0.950 (0.940, 0.959)|0.954 (0.943, 0.962)|
|null n500 (no subgroup; homogeneous +26)|MR (IJ)|2000|1|-5.687|-0.455|1.900|0.994 (0.989, 0.996)|0.999 (0.996, 1.000)|
|null n500 (no subgroup; homogeneous +26)|MR (field)|2000|1|-3.561|-0.280|1.013|0.915 (0.902, 0.927)|0.942 (0.931, 0.951)|
|null n500 (no subgroup; homogeneous +26)|MR (field-s)|2000|1|-3.572|-0.282|1.017|0.920 (0.907, 0.931)|0.944 (0.933, 0.953)|
|md40 n700|naive|2000|1|-11.224|-1.089|1.046|0.745 (0.725, 0.763)|0.831 (0.813, 0.846)|
|md40 n700|oracle|2000|1|0.166|0.014|1.004|0.955 (0.944, 0.963)|0.950 (0.939, 0.958)|
|md40 n700|MR (IJ)|2000|1|-3.585|-0.339|1.889|0.996 (0.991, 0.998)|0.999 (0.996, 1.000)|
|md40 n700|MR (field)|2000|1|-1.987|-0.185|0.990|0.929 (0.917, 0.939)|0.945 (0.934, 0.954)|
|md40 n700|MR (field-s)|2000|1|-1.997|-0.187|0.994|0.931 (0.919, 0.941)|0.947 (0.936, 0.956)|

- The field-s one-sided upper bound on β(Ĥᶜ) covers at 0.920 / 0.931 / 0.920 / 0.931; its Wilson upper limits (0.931 / 0.942 / 0.931 / 0.941) are below 0.95 in every cell. The unstudentized field is 0.916 / 0.926 / 0.915 / 0.929, so field-s adds 0.002–0.005.
- Field-s retained bias is −0.287 / −0.080 / −0.282 / −0.187 SD units with SE/SD 0.914–1.017; IJ covers at 0.994–0.997 with SE/SD 1.671–1.900; the oracle at 0.950–0.956.

### Bound location on Ĥ and Ĥᶜ — the D3 ladder (conditional on the proposed family)

Ĥ: one-sided 95% LOWER bound (field | oracle; IJ and unadjusted for orientation).

|cell|estimator|n|mean|q05|q25|median|q75|q95|P(L>=0)|P(L>=10)|P(L>=20)|P(L>=30)|P(L>=40)|P(L>=60)|P(L>=80)|P(L>=100)|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|md40 n500|MR (field)|2000|-4.2|-37.9|-19.2|-5.3|9.2|33.5|0.397 (0.011)|0.236 (0.009)|0.126 (0.007)|0.067 (0.006)|0.032 (0.004)|0.007 (0.002)|0.001 (0.001)|0.000 (0.000)|
|md40 n500|oracle|2000|6.7|-27.0|-5.9|7.3|20.1|39.8|0.638 (0.011)|0.442 (0.011)|0.252 (0.010)|0.118 (0.007)|0.048 (0.005)|0.003 (0.001)|0.000 (0.000)|0.000 (0.000)|
|md40 n500|MR (IJ)|2000|-8.4|-39.7|-21.4|-9.1|4.0|24.8|0.318 (0.010)|0.167 (0.008)|0.075 (0.006)|0.034 (0.004)|0.013 (0.003)|0.001 (0.001)|0.000 (0.000)|0.000 (0.000)|
|md40 n500|naive|2000|47.8|22.0|36.7|47.0|58.0|76.0|0.999 (0.001)|0.994 (0.002)|0.963 (0.004)|0.867 (0.008)|0.671 (0.011)|0.212 (0.009)|0.034 (0.004)|0.004 (0.001)|
|md120 n500|MR (field)|2000|45.2|6.0|28.1|44.3|60.5|87.4|0.964 (0.004)|0.931 (0.006)|0.848 (0.008)|0.719 (0.010)|0.570 (0.011)|0.258 (0.010)|0.086 (0.006)|0.023 (0.003)|
|md120 n500|oracle|2000|86.7|53.0|74.1|87.3|100.1|119.8|1.000 (0.000)|1.000 (0.000)|0.999 (0.001)|0.997 (0.001)|0.986 (0.003)|0.907 (0.006)|0.638 (0.011)|0.252 (0.010)|
|md120 n500|MR (IJ)|2000|38.5|3.2|24.4|38.3|51.8|74.2|0.961 (0.004)|0.909 (0.006)|0.818 (0.009)|0.661 (0.011)|0.467 (0.011)|0.149 (0.008)|0.035 (0.004)|0.005 (0.002)|
|md120 n500|naive|2000|87.9|59.0|76.4|87.4|99.1|119.0|1.000 (0.000)|1.000 (0.000)|1.000 (0.000)|1.000 (0.000)|0.997 (0.001)|0.945 (0.005)|0.668 (0.011)|0.235 (0.009)|
|null n500 (no subgroup; homogeneous +26)|MR (field)|2000|-9.2|-42.5|-24.6|-10.4|4.6|28.5|0.326 (0.010)|0.180 (0.009)|0.092 (0.006)|0.045 (0.005)|0.022 (0.003)|0.004 (0.001)|0.001 (0.000)|0.000 (0.000)|
|null n500 (no subgroup; homogeneous +26)|oracle|0|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|
|null n500 (no subgroup; homogeneous +26)|MR (IJ)|2000|-13.4|-44.6|-26.9|-13.9|-0.7|20.4|0.239 (0.010)|0.116 (0.007)|0.051 (0.005)|0.018 (0.003)|0.008 (0.002)|0.001 (0.001)|0.000 (0.000)|0.000 (0.000)|
|null n500 (no subgroup; homogeneous +26)|naive|2000|43.3|17.2|32.1|42.7|53.6|71.5|0.999 (0.001)|0.987 (0.003)|0.931 (0.006)|0.798 (0.009)|0.559 (0.011)|0.150 (0.008)|0.021 (0.003)|0.002 (0.001)|
|md40 n700|MR (field)|2000|-3.3|-35.4|-16.9|-5.0|9.3|34.0|0.390 (0.011)|0.241 (0.010)|0.123 (0.007)|0.065 (0.005)|0.032 (0.004)|0.010 (0.002)|0.002 (0.001)|0.001 (0.000)|
|md40 n700|oracle|2000|12.5|-15.5|0.8|12.0|24.0|41.5|0.763 (0.010)|0.552 (0.011)|0.322 (0.010)|0.165 (0.008)|0.059 (0.005)|0.004 (0.001)|0.000 (0.000)|0.000 (0.000)|
|md40 n700|MR (IJ)|2000|-4.8|-34.0|-16.9|-5.4|6.4|27.6|0.378 (0.011)|0.190 (0.009)|0.090 (0.006)|0.040 (0.004)|0.017 (0.003)|0.002 (0.001)|0.001 (0.000)|0.000 (0.000)|
|md40 n700|naive|2000|50.2|26.2|39.9|49.4|59.1|77.5|1.000 (0.000)|0.997 (0.001)|0.981 (0.003)|0.918 (0.006)|0.748 (0.010)|0.232 (0.009)|0.038 (0.004)|0.007 (0.002)|

Ĥᶜ: one-sided 95% UPPER bound (field-s | oracle; field, IJ and unadjusted for orientation).

|cell|estimator|n|mean|q05|q25|median|q75|q95|P(U<=0)|P(U<=10)|P(U<=20)|P(U<=30)|P(U<=40)|P(U<=60)|P(U<=80)|P(U<=100)|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|md40 n500|MR (field-s)|2000|48.4|27.3|40.3|48.4|56.8|69.7|0.000 (0.000)|0.003 (0.001)|0.018 (0.003)|0.070 (0.006)|0.243 (0.010)|0.818 (0.009)|0.995 (0.002)|1.000 (0.000)|
|md40 n500|oracle|2000|49.7|26.8|39.8|49.9|59.5|72.5|0.001 (0.001)|0.006 (0.002)|0.021 (0.003)|0.081 (0.006)|0.253 (0.010)|0.758 (0.010)|0.985 (0.003)|1.000 (0.000)|
|md40 n500|MR (field)|2000|48.4|27.0|40.1|48.3|57.1|69.4|0.000 (0.000)|0.002 (0.001)|0.020 (0.003)|0.075 (0.006)|0.246 (0.010)|0.811 (0.009)|0.994 (0.002)|1.000 (0.000)|
|md40 n500|MR (IJ)|2000|64.1|43.5|56.1|64.1|72.4|84.8|0.000 (0.000)|0.000 (0.000)|0.001 (0.001)|0.004 (0.001)|0.032 (0.004)|0.369 (0.011)|0.895 (0.007)|0.998 (0.001)|
|md40 n500|naive|2000|36.2|16.4|28.6|36.1|44.2|56.4|0.001 (0.001)|0.021 (0.003)|0.090 (0.006)|0.294 (0.010)|0.635 (0.011)|0.976 (0.003)|1.000 (0.000)|1.000 (0.000)|
|md120 n500|MR (field-s)|2000|65.7|41.1|55.8|66.3|75.8|89.4|0.000 (0.000)|0.001 (0.001)|0.003 (0.001)|0.011 (0.002)|0.044 (0.005)|0.339 (0.011)|0.839 (0.008)|0.992 (0.002)|
|md120 n500|oracle|2000|49.7|26.8|39.8|49.9|59.5|72.5|0.001 (0.001)|0.006 (0.002)|0.021 (0.003)|0.081 (0.006)|0.253 (0.010)|0.758 (0.010)|0.985 (0.003)|1.000 (0.000)|
|md120 n500|MR (field)|2000|65.1|38.9|55.1|66.0|75.5|89.4|0.000 (0.000)|0.001 (0.001)|0.003 (0.001)|0.014 (0.003)|0.054 (0.005)|0.360 (0.011)|0.844 (0.008)|0.991 (0.002)|
|md120 n500|MR (IJ)|2000|82.3|58.2|72.8|83.0|92.1|105.5|0.000 (0.000)|0.000 (0.000)|0.000 (0.000)|0.001 (0.001)|0.003 (0.001)|0.062 (0.005)|0.424 (0.011)|0.901 (0.007)|
|md120 n500|naive|2000|54.6|31.8|45.7|55.3|64.0|77.5|0.001 (0.000)|0.003 (0.001)|0.011 (0.002)|0.041 (0.004)|0.147 (0.008)|0.649 (0.011)|0.970 (0.004)|1.000 (0.000)|
|null n500 (no subgroup; homogeneous +26)|MR (field-s)|2000|43.9|23.2|35.9|43.7|52.2|64.8|0.001 (0.001)|0.006 (0.002)|0.030 (0.004)|0.135 (0.008)|0.376 (0.011)|0.894 (0.007)|0.998 (0.001)|1.000 (0.000)|
|null n500 (no subgroup; homogeneous +26)|oracle|2000|45.0|26.3|37.4|45.1|52.6|64.0|0.000 (0.000)|0.003 (0.001)|0.018 (0.003)|0.093 (0.007)|0.334 (0.011)|0.897 (0.007)|1.000 (0.000)|1.000 (0.000)|
|null n500 (no subgroup; homogeneous +26)|MR (field)|2000|43.9|22.5|35.7|43.7|52.3|65.0|0.001 (0.001)|0.008 (0.002)|0.030 (0.004)|0.135 (0.008)|0.376 (0.011)|0.896 (0.007)|0.998 (0.001)|1.000 (0.000)|
|null n500 (no subgroup; homogeneous +26)|MR (IJ)|2000|59.6|39.4|51.6|59.5|67.8|80.4|0.000 (0.000)|0.000 (0.000)|0.001 (0.001)|0.012 (0.002)|0.056 (0.005)|0.517 (0.011)|0.946 (0.005)|1.000 (0.000)|
|null n500 (no subgroup; homogeneous +26)|naive|2000|31.8|12.2|24.1|31.6|39.7|51.5|0.007 (0.002)|0.033 (0.004)|0.155 (0.008)|0.446 (0.011)|0.756 (0.010)|0.992 (0.002)|1.000 (0.000)|1.000 (0.000)|
|md40 n700|MR (field-s)|2000|46.3|29.3|39.1|46.1|53.7|63.8|0.000 (0.000)|0.001 (0.000)|0.007 (0.002)|0.056 (0.005)|0.279 (0.010)|0.895 (0.007)|0.999 (0.001)|1.000 (0.000)|
|md40 n700|oracle|2000|46.5|27.2|38.4|46.5|54.6|66.5|0.000 (0.000)|0.002 (0.001)|0.012 (0.002)|0.086 (0.006)|0.295 (0.010)|0.872 (0.007)|0.997 (0.001)|1.000 (0.000)|
|md40 n700|MR (field)|2000|46.3|28.8|39.0|46.1|53.8|64.1|0.000 (0.000)|0.001 (0.001)|0.007 (0.002)|0.055 (0.005)|0.272 (0.010)|0.894 (0.007)|0.999 (0.001)|1.000 (0.000)|
|md40 n700|MR (IJ)|2000|60.1|43.3|53.2|60.0|67.3|77.6|0.000 (0.000)|0.000 (0.000)|0.000 (0.000)|0.004 (0.001)|0.028 (0.004)|0.501 (0.011)|0.969 (0.004)|1.000 (0.000)|
|md40 n700|naive|2000|37.4|20.7|30.8|37.1|44.4|54.4|0.001 (0.000)|0.006 (0.002)|0.043 (0.005)|0.231 (0.009)|0.609 (0.011)|0.988 (0.002)|1.000 (0.000)|1.000 (0.000)|

- On Ĥ the field lower bound sits at or above 0 on 0.397 / 0.964 / 0.326 / 0.390 of declared replicates and at or above 40 on 0.032 / 0.570 / 0.022 / 0.032; the oracle's sits at or above 40 on 0.048 / 0.986 / – / 0.059. At md120 the field bound reaches 100 on 0.023 of replicates against the oracle's 0.252.
- In the null cell (truth +26 everywhere) the field lower bound sits at or above 30 on 0.045 of replicates; the unadjusted bound on 0.798.
- On Ĥᶜ the field-s upper bound sits at or below 30 on 0.070 / 0.011 / 0.135 / 0.056 and at or below 60 on 0.818 / 0.339 / 0.894 / 0.895, beside the oracle's 0.081 / 0.081 / 0.093 / 0.086 and 0.758 / 0.758 / 0.897 / 0.872.

### Joint pair (Ĥ lower, Ĥᶜ upper) (conditional on the proposed family)

|cell|pair|declared|both_bounds|share_both|joint|joint_mc_se|cov_H|cov_Hc|margin_H|margin_Hc|
|---|---|---|---|---|---|---|---|---|---|---|
|md40 n500|Bonferroni field-s (gamma = 0.025)|2000|2000|1|0.919 (0.907, 0.931)|0.006|0.969|0.950|60.754|27.477|
|md40 n500|Bonferroni unstudentized (gamma = 0.025)|2000|2000|1|0.921 (0.909, 0.932)|0.006|0.969|0.952|60.754|27.500|
|md40 n500|separate 95% bounds: field lower, field-s upper|2000|2000|1|0.864 (0.848, 0.878)|0.008|0.939|0.920|51.168|23.405|
|md120 n500|Bonferroni field-s (gamma = 0.025)|2000|2000|1|0.943 (0.931, 0.952)|0.005|0.980|0.962|56.364|28.866|
|md120 n500|Bonferroni unstudentized (gamma = 0.025)|2000|2000|1|0.942 (0.930, 0.951)|0.005|0.980|0.961|56.364|28.246|
|md120 n500|separate 95% bounds: field lower, field-s upper|2000|2000|1|0.893 (0.879, 0.906)|0.007|0.956|0.931|46.862|24.474|
|null n500 (no subgroup; homogeneous +26)|Bonferroni field-s (gamma = 0.025)|2000|2000|1|0.930 (0.918, 0.940)|0.006|0.972|0.957|61.404|27.342|
|null n500 (no subgroup; homogeneous +26)|Bonferroni unstudentized (gamma = 0.025)|2000|2000|1|0.927 (0.915, 0.938)|0.006|0.972|0.955|61.404|27.412|
|null n500 (no subgroup; homogeneous +26)|separate 95% bounds: field lower, field-s upper|2000|2000|1|0.865 (0.849, 0.879)|0.008|0.941|0.920|51.767|23.289|
|md40 n700|Bonferroni field-s (gamma = 0.025)|2000|2000|1|0.934 (0.923, 0.945)|0.006|0.971|0.963|60.374|22.373|
|md40 n700|Bonferroni unstudentized (gamma = 0.025)|2000|2000|1|0.933 (0.921, 0.943)|0.006|0.971|0.961|60.374|22.381|
|md40 n700|separate 95% bounds: field lower, field-s upper|2000|2000|1|0.877 (0.862, 0.891)|0.007|0.943|0.931|50.634|19.061|

|cell|metric|value|mc_se|n|
|---|---|---|---|---|
|md40 n500|gamma_mean_s|0.0252|0.0000|2000|
|md40 n500|corr_s|0.0923|0.0014|2000|
|md40 n500|gamma_mean|0.0251|0.0000|2000|
|md40 n500|corr|0.0970|0.0014|2000|
|md120 n500|gamma_mean_s|0.0252|0.0000|2000|
|md120 n500|corr_s|0.0256|0.0016|2000|
|md120 n500|gamma_mean|0.0251|0.0000|2000|
|md120 n500|corr|0.0300|0.0017|2000|
|null n500 (no subgroup; homogeneous +26)|gamma_mean_s|0.0252|0.0000|2000|
|null n500 (no subgroup; homogeneous +26)|corr_s|0.0932|0.0014|2000|
|null n500 (no subgroup; homogeneous +26)|gamma_mean|0.0251|0.0000|2000|
|null n500 (no subgroup; homogeneous +26)|corr|0.0976|0.0014|2000|
|md40 n700|gamma_mean_s|0.0251|0.0000|2000|
|md40 n700|corr_s|0.0926|0.0013|2000|
|md40 n700|gamma_mean|0.0251|0.0000|2000|
|md40 n700|corr|0.0965|0.0013|2000|

- The field-s Bonferroni pair covers (β(Ĥ), β(Ĥᶜ)) jointly on 0.919 / 0.943 / 0.930 / 0.934 of declared replicates (Wilson upper limits 0.931 / 0.952 / 0.940 / 0.945); both bounds exist on every declared replicate. The unstudentized pair is within 0.003 of it; the separate 95% pair covers at 0.864–0.893.
- The calibrated split returns the Bonferroni floor (mean γ 0.0251–0.0252) with corr(Λ*, Λ*ᶜ) 0.0256–0.0932 (field-s).

### Identification, GRF beside FS (conditional on the proposed family)

|cell|quantity|GRF|FS_mdsgnb20|
|---|---|---|---|
|md40 n500|declaration rate|1.0000|0.9990|
|md40 n500|mean size of Hhat (n_sel)|109.46|111.62|
|md40 n500|mean true positives (sens x n_true)|43.02|-|
|md40 n500|sensitivity (mean over declared; NA where n_true = 0)|0.2497|0.2679|
|md40 n500|PPV (mean over declared)|0.3909|0.4122|
|md40 n500|mean family size K (n_family, MR’s kept family)|1222.4|-|
|md40 n500|mean admitted set (admitted_n)|624.6|-|
|md40 n500|size of Hhat GRF larger / equal / smaller than FS (paired by sim_id, both declared)|865 / 78 / 1055|-|
|md120 n500|declaration rate|1.0000|1.0000|
|md120 n500|mean size of Hhat (n_sel)|143.61|152.74|
|md120 n500|mean true positives (sens x n_true)|103.13|-|
|md120 n500|sensitivity (mean over declared; NA where n_true = 0)|0.5989|0.7116|
|md120 n500|PPV (mean over declared)|0.7149|0.7980|
|md120 n500|mean family size K (n_family, MR’s kept family)|1222.4|-|
|md120 n500|mean admitted set (admitted_n)|1059.8|-|
|md120 n500|size of Hhat GRF larger / equal / smaller than FS (paired by sim_id, both declared)|752 / 41 / 1207|-|
|null n500 (no subgroup; homogeneous +26)|declaration rate|1.0000|0.9965|
|null n500 (no subgroup; homogeneous +26)|mean size of Hhat (n_sel)|106.14|107.49|
|null n500 (no subgroup; homogeneous +26)|mean true positives (sens x n_true)|NA|-|
|null n500 (no subgroup; homogeneous +26)|sensitivity (mean over declared; NA where n_true = 0)|NA|NA|
|null n500 (no subgroup; homogeneous +26)|PPV (mean over declared)|0.0000|0.0000|
|null n500 (no subgroup; homogeneous +26)|mean family size K (n_family, MR’s kept family)|1222.4|-|
|null n500 (no subgroup; homogeneous +26)|mean admitted set (admitted_n)|508.5|-|
|null n500 (no subgroup; homogeneous +26)|size of Hhat GRF larger / equal / smaller than FS (paired by sim_id, both declared)|865 / 87 / 1041|-|
|md40 n700|declaration rate|1.0000|0.9995|
|md40 n700|mean size of Hhat (n_sel)|116.33|117.99|
|md40 n700|mean true positives (sens x n_true)|45.16|-|
|md40 n700|sensitivity (mean over declared; NA where n_true = 0)|0.1873|0.2041|
|md40 n700|PPV (mean over declared)|0.3863|0.4101|
|md40 n700|mean family size K (n_family, MR’s kept family)|1354.5|-|
|md40 n700|mean admitted set (admitted_n)|711.4|-|
|md40 n700|size of Hhat GRF larger / equal / smaller than FS (paired by sim_id, both declared)|877 / 57 / 1065|-|

- Mean |Ĥ| is 109.46 / 143.61 / 106.14 / 116.33 patients for GRF against FS's 111.62 / 152.74 / 107.49 / 117.99; paired by `sim_id`, GRF's Ĥ is larger on 865 / 752 / 865 / 877, equal on 78 / 41 / 87 / 57 and smaller on 1,055 / 1,207 / 1,041 / 1,065 replicates.
- GRF's Ĥ holds 43.02 / 103.13 / – / 45.16 truly harmed patients on average (sensitivity × `n_true`; never |Ĥ|), with sensitivity 0.2497 / 0.5989 / – / 0.1873 and PPV 0.3909 / 0.7149 / 0 / 0.3863, against FS's 0.2679 / 0.7116 / – / 0.2041 and 0.4122 / 0.7980 / 0 / 0.4101. In the null cell PPV is 0 by construction and sensitivity undefined.
- The family MR re-selects over averages 1,222.4 candidates in the three n = 500 cells and 1,354.5 at n = 700 (the enumerated pool, which depends on the covariates only, and the three n = 500 cells share every covariate draw); the admitted set averages 624.6 / 1,059.8 / 508.5 / 711.4.

### Regime diagnostics (conditional on the proposed family)

|cell|n_det|p_hat_mean|p_hat_lt05|sd_btc_naive|lamc_naive|lamc_s_naive|ij_sd_H|ij_sd_Hc|fit_secs|field_secs|comp_secs|
|---|---|---|---|---|---|---|---|---|---|---|---|
|md40 n500|2000|0.099|0.998|0.957|0.978|0.982|1.628|1.885|45.148|27.202|1.719|
|md120 n500|2000|0.128|0.996|1.064|0.969|0.994|1.382|1.671|51.187|32.274|1.862|
|null n500 (no subgroup; homogeneous +26)|2000|0.103|0.998|0.952|0.980|0.982|1.639|1.900|44.115|26.252|1.708|
|md40 n700|2000|0.099|1.000|0.980|0.984|0.986|1.587|1.889|53.470|29.309|2.389|

- p̂(Ĥ) is 0.099 / 0.128 / 0.103 / 0.099, with 0.996–1.000 of replicates below 0.5: the tie regime.
- SD(β̃ᶜ) / mean naive SEᶜ is 0.952–1.064 and λ-SDᶜ / naive SEᶜ 0.969–0.984 (field), 0.982–0.994 (field-s): the complement is not in the moved regime.
- IJ SE / empirical SD is 1.382–1.639 on Ĥ and 1.671–1.900 on Ĥᶜ; seconds per replicate 44.115–53.470 (field 26.252–32.274; complement 1.708–2.389).

### The display (identity scale) (conditional on the proposed family)

|cell|block|estimator|n|bias (MD)|SD|mean SE|b|r|1-sided cov|1-sided ref|2-sided cov|2-sided ref|
|---|---|---|---|---|---|---|---|---|---|---|---|---|
|md40 n500|H|naive|2000|58.624|17.801|25.810|3.293|1.450|0.158|0.182|0.320|0.326|
|md40 n500|H|mr|2000|15.422|20.704|33.703|0.745|1.628|0.969|0.973|0.992|0.993|
|md40 n500|H|fld|2000|7.400|22.230|25.196|0.333|1.133|0.939|0.937|0.965|0.965|
|md40 n500|Hc|naive|2000|-16.315|12.234|13.185|-1.334|1.078|0.677|0.670|0.799|0.782|
|md40 n500|Hc|mr|2000|-5.810|12.616|23.783|-0.461|1.885|0.995|0.996|0.999|0.999|
|md40 n500|Hc|fld|2000|-3.666|12.826|12.902|-0.286|1.006|0.916|0.914|0.942|0.942|
|md40 n500|Hc|fld_s|2000|-3.677|12.797|12.943|-0.287|1.011|0.920|0.916|0.938|0.943|
|md120 n500|H|naive|2000|32.906|19.567|22.794|1.682|1.165|0.599|0.593|0.713|0.726|
|md120 n500|H|mr|2000|-0.455|23.552|32.555|-0.019|1.382|0.989|0.989|0.983|0.993|
|md120 n500|H|fld|2000|-4.088|25.431|25.544|-0.161|1.004|0.956|0.965|0.926|0.948|
|md120 n500|Hc|naive|2000|-12.443|14.332|14.046|-0.868|0.980|0.790|0.772|0.866|0.851|
|md120 n500|Hc|mr|2000|-2.768|14.940|24.969|-0.185|1.671|0.997|0.995|1.000|0.999|
|md120 n500|Hc|fld|2000|-1.212|15.293|13.607|-0.079|0.890|0.926|0.917|0.929|0.918|
|md120 n500|Hc|fld_s|2000|-1.228|15.276|13.959|-0.080|0.914|0.931|0.923|0.934|0.926|
|null n500 (no subgroup; homogeneous +26)|H|naive|2000|60.180|17.789|26.205|3.383|1.473|0.143|0.169|0.300|0.310|
|null n500 (no subgroup; homogeneous +26)|H|mr|2000|16.271|20.758|34.023|0.784|1.639|0.972|0.972|0.991|0.992|
|null n500 (no subgroup; homogeneous +26)|H|fld|2000|7.966|22.337|25.371|0.357|1.136|0.941|0.935|0.969|0.964|
|null n500 (no subgroup; homogeneous +26)|Hc|naive|2000|-16.094|12.106|13.127|-1.329|1.084|0.688|0.675|0.798|0.787|
|null n500 (no subgroup; homogeneous +26)|Hc|mr|2000|-5.687|12.491|23.737|-0.455|1.900|0.994|0.996|0.999|0.999|
|null n500 (no subgroup; homogeneous +26)|Hc|fld|2000|-3.561|12.700|12.868|-0.280|1.013|0.915|0.917|0.942|0.944|
|null n500 (no subgroup; homogeneous +26)|Hc|fld_s|2000|-3.572|12.672|12.886|-0.282|1.017|0.920|0.918|0.944|0.945|
|md40 n700|H|naive|2000|60.124|17.429|25.229|3.450|1.447|0.100|0.143|0.236|0.270|
|md40 n700|H|mr|2000|15.787|19.979|31.700|0.790|1.587|0.967|0.966|0.989|0.990|
|md40 n700|H|fld|2000|7.596|21.414|24.601|0.355|1.149|0.943|0.938|0.967|0.967|
|md40 n700|Hc|naive|2000|-11.224|10.306|10.783|-1.089|1.046|0.745|0.736|0.831|0.831|
|md40 n700|Hc|mr|2000|-3.585|10.563|19.950|-0.339|1.889|0.996|0.997|0.999|1.000|
|md40 n700|Hc|fld|2000|-1.987|10.714|10.605|-0.185|0.990|0.929|0.925|0.945|0.944|
|md40 n700|Hc|fld_s|2000|-1.997|10.702|10.637|-0.187|0.994|0.931|0.926|0.947|0.945|

- Field points on Ĥ against the Gaussian reference: 0.939 vs 0.937, 0.956 vs 0.965, 0.941 vs 0.935, 0.943 vs 0.938 (b −0.161 to +0.357, r 1.004–1.149); the md120 two-sided point (0.926 vs 0.948) is the largest departure on this block.
- Field-s points on Ĥᶜ: 0.920 vs 0.916, 0.931 vs 0.923, 0.920 vs 0.918, 0.931 vs 0.926 (b −0.287 to −0.080, r 0.914–1.017).
- Figures: `fig_mdgrf_bias_coverage_display_H.png`, `fig_mdgrf_bias_coverage_display_Hc.png` (also embedded in the summary HTML), captioned with the conditional-family sentence.

## The confound

GRF's candidate family is generated from a fitted surface, so the fixed-family condition does not hold: every coverage figure below is coverage of the estimand conditional on the proposed family. Comparisons with FS are descriptive, not a contest: the identifier and the family construction differ, and each summary conditions on a different set of detected replicates.

## Scope

These are operating characteristics of the GRF identifier on one continuous design, conditional on the proposed family. They do not verify condition (A3), no construction is promoted on this design's performance, and DINA was not run.

## Findings

1. **Stop and resume.** The first Stage 1 stopped on the factor-membership defect; P1 (`0cd33f7b`) fixed the evaluator, and the resumed Stage 1 found zero factor warnings and zero NA memberships on every smoke replicate. Stage 2 ran under the advance go without a halt.
2. **Cost.** 7,423 s for 8,000 replicates at 63 workers (123.7 min against the 131.9-min projection); per replicate 44.1–53.5 s, the field pass 26.3–32.3 s of it — the post-fix family (about 1,200–1,350 candidates) roughly doubles the pre-fix GRF cost measured in the first smoke.
3. **Family size depends on the covariates only.** MR's kept family equals GRF's enumerated pool (min 1,170, median 1,199 at n = 500; median 1,337 at n = 700) and is identical across the three n = 500 cells, which share every covariate draw; the admitted set does depend on the outcome.
4. **Declaration is total.** GRF declares on every replicate of every cell, the null included, where FS left 2 / 0 / 7 / 1 undeclared; the admitted set is never empty (minimum 15 / 327 / 7 / 17).
5. **Complement one-sided coverage is below nominal in every cell** for field-s (Wilson upper limits 0.931–0.942), conditional on the proposed family. Reported, not interpreted.
6. **`dmin.grf = 30`** is recorded and applied; under the effect re-selection it floors only the DR-score pre-filter (first Stage 1 record, F1), so the admission is `effect.threshold`'s 30 on the harm-oriented MD.
7. **Convention 9.** The recorder's `n_harm` is |Ĥ| (template `:772`), and the FS extract's `mean_n_harm` equals its `mean_n_sel`; the GRF extract reports |Ĥ| as `mean_size_hhat` and the true-positive count as sensitivity × `n_true`, and does not copy `mean_n_harm`. The FS extract was not re-run.
8. **Attribution in the runner's commits.** The committed runner (`f0b9c844`, frozen by S2) writes the trailer `Co-Authored-By: Claude Fable 5.1`; the four cell commits and the campaign-complete commit carry it.
9. **Catalog inventory.** `scripts_mdsgnb20/status_inventory.R` has no `mdgrf` rows, so the campaign's bundles, renders and scripts fall under its generic rows in `current_status.md` §3; the script was not changed.
10. **Rounding.** The md40 n700 unadjusted two-sided coverage prints 0.235 in the Ĥ table (Wilson formatting) and 0.236 in the display table (kable rounding) — the same value, 0.2355, by two routes.
11. **Extract, two rounds.** The first render's extract lacked the regime table's timing columns; the summary was amended to export them and re-rendered (`0333248d`), its tables unchanged, so every number the record's tables print is in the CSV.
12. **Raw campaign logs** stay untracked under `logs_mdgrf/` (runner log, per-render logs and peaks, `GATE2_*.txt`, batch HTML); the closeout deletes the smoke, calibration, check and dry-run items only.

## Commits of this campaign (`git log --oneline`, oldest last; the catalog commits follow this record)

```
0333248d mdgrf summary and extract: export the regime table's timing columns (secs_fit_mr_mean, secs_field_mean, secs_complement_mean) so ev
31360161 mdgrf extract (TASK_md_grf_2026-09-16 §3.2): md_grf_metrics.csv, 1,110 rows (1,094 grf + 16 FS comparator rows copied from md_field
a441fa86 mdgrf summary (TASK_md_grf_2026-09-16 §3.1): summary_continuous_field_mdgrf.qmd, the transplant of summary_continuous_field_mdsgnb2
1a65b28a mdgrf: campaign complete (cumulative render wall 7423 s)
432f69dc mdgrf md40_n700: md 40 n 700, 2000 replicates + combine; Gate 2 PASS (GATE_COUNTS run=66 passed=66 failed=0)
15bee7cd mdgrf null_n500: md null n 500, 2000 replicates + combine; Gate 2 PASS (GATE_COUNTS run=66 passed=66 failed=0)
f84261e3 mdgrf md120_n500: md 120 n 500, 2000 replicates + combine; Gate 2 PASS (GATE_COUNTS run=66 passed=66 failed=0)
c4c49572 mdgrf md40_n500: md 40 n 500, 2000 replicates + combine; Gate 2 PASS (GATE_COUNTS run=66 passed=66 failed=0)
f5cc256d GRF on the MD design, campaign mdgrf, Stage 1 resumed after the membership fix (TASK_md_grf_resume_2026-09-16): S1, S2 pass; stoppe
b62ae1aa Add TASK_md_grf_resume_2026-09-16 as received
…   (the fix task, 0a697abb..019be60f, between the two stages)
a4c063bf GRF on the ACTG175 continuous (MD) design, campaign mdgrf, Stage 1 record (TASK_md_grf_2026-09-16): gates 1.1, 1.2, 1.3, 1.6(a) (FS
f0b9c844 scripts_mdgrf (TASK_md_grf_2026-09-16 §1.6, transplants of scripts_mdsgnb20): mem_sampler.sh unchanged; run_mdgrf.sh (the FS runner
894da993 MD template E1-E4 (TASK_md_grf_2026-09-16 §1.5, transplanted from the survival m1 template): FS_MD_METHOD identifier knob (default 
e5ee1008 Add TASK_md_grf_2026-09-16 as received
<this record; then status_curated.md; then current_status.md alone>
```
