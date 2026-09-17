# REPORT — DINA on the ACTG175 continuous (MD) design: campaign `mddina`

Date: 2026-09-17 (UTC). Machine: `pop-os` (64 physical cores, 251 GB; R 4.6.1, reference BLAS). Branch `feature/glm-extension`. Task: `dev/tasks/TASK_md_dina_campaign_2026-09-17.md` (`189d4ec2`), which transplants the mechanics of `dev/tasks/TASK_md_grf_2026-09-16.md` and `TASK_md_grf_resume_2026-09-16.md`. Records of this campaign: Gate 1 `REPORT_md_dina_stage1_2026-09-17.md` (`37395393`), Gate 2 `REPORT_md_dina_gate2_2026-09-17.md` (per cell, committed by the runner), this record. The fix it depends on: `REPORT_grf_dina_fixes_2026-09-16.md` (P2 `064fce91`). Installed forestsearch 0.3.5, `Built: R 4.6.1; ; 2026-09-17 04:47:31 UTC; unix`; no `R/` change, no install. Every render ran with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`.

DINA's and GRF's candidate families are generated from fitted surfaces, so the fixed-family condition does not hold for them: every DINA and GRF coverage figure below is coverage of the estimand conditional on the proposed family. FS's family is the prespecified cut grid. The three identifiers are compared descriptively, not as a contest: the identifier, the family construction and the set of detected replicates each summary conditions on all differ.

## Dispositions (Larry, 2026-09-17), as executed

- **Identifier:** DINA only, through the MD template's identifier knob from `894da993`, `FS_MD_METHOD`, at its DINA value `dina` (`subgroup_method = "dina"`). This runs DINA's identifier path (`R/forestsearch_main.R:2238`; `.forestsearch_dina_select()`, the lines P2 changed), not the `use_dina` screening path under FS (Gate 1 record §1.3).
- **Rule:** `effMaxSG`, `effect_neighborhood = 0.20`, `selection_rule = "neighborhood"` (`mdsgnb20`'s `meta`), through `FS_MD_FOCUS` / `FS_MD_NBHD`.
- **DINA floors:** as the package applies them after P2 — the proposal floor `m_diff = 30` on −tau-hat and the admission floor 30 on the harm-oriented MD, both at the effect threshold. Both were read as applied on every smoke replicate (Gate 1 record §1.6(c)), and every declared campaign replicate has its smallest proposed oriented tau-hat at or above 30 (Gate 2).
- **Other DINA arguments:** `dinamr`'s — `dina_args = list()`, `dina_select_statistic = "effect"` — recorded in every bundle's `meta`. No survival-only argument was passed, so none was omitted (Gate 1 record §1.4).
- **MR:** `ci_method = "field"` passed explicitly, `draws = 5000`, `include_complement`, `field_complement`, `field_scale_complement = "selected"`, `return_reselection = TRUE`, `ij_residual = "two_term"`, `confirm_rule = "point"`, `t_confirm` near-null. None is hard-coded on the DINA path, and the re-selection is aligned (Gate 1 record §1.3).
- **Labelling:** the conditional-family sentence opens this record and the summary, and ends every summary table and figure caption.
- **Carried fixes:** the runner writes no hard-coded model or co-author trailer (`d65f3b9b`; the four cell commits and the campaign-complete commit carry none). The catalog generator gained the `mdgrf` and `mddina` rows (`7e26ad42`).
- **Gate 1 advance go:** its condition held (every Stage 1 gate green; projection 6.10 h, under 8 h), so Stages 2 and 3 ran.

## Gate 1 and Gate 2 in brief

- **Gate 1** (`37395393`):
  - Template edits `0927669c` (DINA recorder fields and `meta` keys; the argument block was already present). Scripts `d65f3b9b`.
  - FS and GRF regressions identical to `mdsgnb20` and `mdgrf` on 20 of 20 replicates (max relative difference 0, zero selection flips).
  - DINA smoke: 20 of 20 declared; sim_id 1 = the fix record's F5 after P2 (`{cd40 >= 400} & {cd80 >= 1040}`, n 78; 8,324 searched, 2,690 proposed, 1,768 admitted); both floors at 30 on the harm-oriented scale; no proposed candidate below oriented tau-hat 30; zero factor-comparison warnings.
  - Calibration: W = 63 (14.26 replicates per minute; peak 101.7 GB); projection 21,959 s (6.10 h); ceiling 32,939 s; timeout 6,864 s.
- **Gate 2** (`REPORT_md_dina_gate2_2026-09-17.md`, `scripts_mddina/gate2.R`): every cell passed **65 of 65** checks.
  - 2,000 rows with `sim_id` exactly 1–2000, and the batch files match the combined bundle on all 175 columns.
  - `meta` carries `dina`, `dina_select_statistic effect`, `dina_args list()`, `effMaxSG` / 0.20 / `neighborhood`, `ci_method field`, `field_scale_complement selected`, `pkg_version 0.3.5`, `hostname pop-os`, 63 workers.
  - The same draws as `mdsgnb20` in both directions: `n_true` identical on 2,000 rows; oracle columns equal, max relative difference 0 (complement only in the null cell).
  - The E2 fields are filled on every declared replicate, and every proposed candidate sits at oriented tau-hat ≥ 30 (minimum 30.0000 in every cell).
  - Field-s and recorder checks pass; **0 MR failures** and **0 rows with a captured warning** in every cell.
  - Declared: 1,996 / 2,000 / 1,987 / 1,997.
  - Cell commits `1856c107`, `c9f1b948`, `9051b702`, `21e48ee5`; campaign complete `30b282eb`. No halt.
- **Render walls** (the progress log, UTC):
```
2026-09-17T08:53:08Z	campaign	start	HEAD=37395393 workers=63 timeout_s=6864 ceiling_s=32939 built=R 4.6.1; ; 2026-09-17 04:47:31 UTC; unix
2026-09-17T08:53:08Z	md40_n500	start	stem=dina_effMaxSG_mr_field_md40_knoise0_n500_nb20_mddina
2026-09-17T09:36:17Z	md40_n500	batch_1_1000	wall_s=2589 peak_mb=109791 rc=0 cumulative_s=2589
2026-09-17T10:19:06Z	md40_n500	batch_1001_2000	wall_s=2569 peak_mb=110447 rc=0 cumulative_s=5158
2026-09-17T10:19:47Z	md40_n500	combine_1_2000	wall_s=41 peak_mb=1185 rc=0 cumulative_s=5199
2026-09-17T10:19:47Z	md40_n500	done	cell_wall_s=5199 GATE_COUNTS run=65 passed=65 failed=0
2026-09-17T10:19:47Z	md120_n500	start	stem=dina_effMaxSG_mr_field_md120_knoise0_n500_nb20_mddina
2026-09-17T11:45:53Z	md120_n500	batch_1_1000	wall_s=5166 peak_mb=121791 rc=0 cumulative_s=10365
2026-09-17T13:11:25Z	md120_n500	batch_1001_2000	wall_s=5132 peak_mb=120815 rc=0 cumulative_s=15497
2026-09-17T13:12:00Z	md120_n500	combine_1_2000	wall_s=35 peak_mb=1112 rc=0 cumulative_s=15532
2026-09-17T13:12:01Z	md120_n500	done	cell_wall_s=10334 GATE_COUNTS run=65 passed=65 failed=0
2026-09-17T13:12:01Z	null_n500	start	stem=dina_effMaxSG_mr_field_mdnull_knoise0_n500_nb20_mddina
2026-09-17T13:47:11Z	null_n500	batch_1_1000	wall_s=2110 peak_mb=102441 rc=0 cumulative_s=17642
2026-09-17T14:20:00Z	null_n500	batch_1001_2000	wall_s=1969 peak_mb=103052 rc=0 cumulative_s=19611
2026-09-17T14:20:40Z	null_n500	combine_1_2000	wall_s=40 peak_mb=1106 rc=0 cumulative_s=19651
2026-09-17T14:20:41Z	null_n500	done	cell_wall_s=4120 GATE_COUNTS run=65 passed=65 failed=0
2026-09-17T14:20:41Z	md40_n700	start	stem=dina_effMaxSG_mr_field_md40_knoise0_n700_nb20_mddina
2026-09-17T15:29:21Z	md40_n700	batch_1_1000	wall_s=4120 peak_mb=121470 rc=0 cumulative_s=23771
2026-09-17T16:34:37Z	md40_n700	batch_1001_2000	wall_s=3916 peak_mb=123321 rc=0 cumulative_s=27687
2026-09-17T16:35:17Z	md40_n700	combine_1_2000	wall_s=40 peak_mb=1199 rc=0 cumulative_s=27727
2026-09-17T16:35:18Z	md40_n700	done	cell_wall_s=8077 GATE_COUNTS run=65 passed=65 failed=0
2026-09-17T16:35:18Z	campaign	complete	cumulative_s=27727
```
  **27,727 s (462.1 min, 7.70 h) in total** — 26% over the 21,959-s projection and under the 32,939-s ceiling. The md120 cell took 10,334 s against its projected 5,664 s (finding 2).
  - Peak summed RSS (63 workers): 109.8–110.4 GB (md40 n500), 120.8–121.8 GB (md120), 102.4–103.1 GB (null), 121.5–123.3 GB (md40 n700).
  - `fit_mr_secs` mean: 134.2 / 304.9 / 100.3 / 199.7 s (md40 n500 / md120 / null / md40 n700).
  - Of which the field pass: 69.8 / 149.2 / 54.4 / 88.0 s; the complement field: 7.00 / 17.65 / 4.95 / 12.56 s.
  - DINA fit (`id_secs`): 12.05 / 21.22 / 9.32 / 14.55 s.

## Tables (pasted from the render of `summary_continuous_field_mddina.qmd`, `1ffef474`; every number is in `md_dina_metrics.csv`, `34067b7d`)

Conventions:
- Scale: harm-oriented MD (positive = harm).
- Targets: β(Ĥ) and β(Ĥᶜ), exact per replicate; the oracle is scored against the structural true-region effect.
- Coverage is one-sided on the exposed side (Ĥ: the LOWER bound; Ĥᶜ: the UPPER bound).
- Parentheses hold Wilson 95% intervals in the coverage tables and Monte Carlo SEs in the ladder tables.
- The null cell is labelled by its truth (no subgroup; homogeneous +26); its oracle row on Ĥ is blank because Q is empty.
- **Every table is conditional on the proposed family.**

### Declaration (conditional on the proposed family)

|cell|replicates|declared|rate|mc_se|fs_rate|
|---|---|---|---|---|---|
|md40 n500|2000|1996|0.9980|0.0010|0.9990|
|md120 n500|2000|2000|1.0000|0.0000|1.0000|
|null n500 (no subgroup; homogeneous +26)|2000|1987|0.9935|0.0018|0.9965|
|md40 n700|2000|1997|0.9985|0.0009|0.9995|

- DINA declares on 1,996 / 2,000 / 1,987 / 1,997 of 2,000 replicates (rates 0.9980 / 1.0000 / 0.9935 / 0.9985); FS (`mdsgnb20`) on 0.9990 / 1.0000 / 0.9965 / 0.9995 on the same seeds; GRF (`mdgrf`) on every replicate.
- In the null cell every declared Ĥ has β(Ĥ) = +26.26, so its rows below are read as coverage of a homogeneous harm, not as a false-claim rate.

### Ĥ — unadjusted, oracle, IJ two-term, field (conditional on the proposed family)

|cell|estimator|n|declaration rate|bias (MD)|bias (SD units)|SE/SD|one-sided LOWER coverage (Wilson)|two-sided coverage (Wilson)|
|---|---|---|---|---|---|---|---|---|
|md40 n500|naive|1996|0.998|61.464|3.591|1.500|0.111 (0.098, 0.125)|0.244 (0.226, 0.264)|
|md40 n500|oracle|1996|0.998|-0.662|-0.033|0.991|0.952 (0.942, 0.961)|0.943 (0.932, 0.952)|
|md40 n500|MR (IJ)|1996|0.998|19.566|1.084|1.919|0.976 (0.968, 0.982)|0.993 (0.988, 0.996)|
|md40 n500|MR (field)|1996|0.998|11.615|0.600|1.278|0.937 (0.926, 0.947)|0.971 (0.963, 0.978)|
|md120 n500|naive|2000|1.000|31.462|1.727|1.225|0.596 (0.574, 0.617)|0.705 (0.685, 0.725)|
|md120 n500|oracle|2000|1.000|-0.684|-0.034|0.991|0.953 (0.942, 0.961)|0.943 (0.932, 0.952)|
|md120 n500|MR (IJ)|2000|1.000|-2.764|-0.122|1.419|0.991 (0.985, 0.994)|0.984 (0.977, 0.988)|
|md120 n500|MR (field)|2000|1.000|-5.940|-0.238|1.010|0.965 (0.955, 0.972)|0.909 (0.896, 0.921)|
|null n500 (no subgroup; homogeneous +26)|naive|1987|0.994|62.973|3.602|1.493|0.115 (0.101, 0.130)|0.233 (0.215, 0.252)|
|null n500 (no subgroup; homogeneous +26)|oracle|0|0.994|NaN|NaN|NA|NaN (NA, NA)|NaN (NA, NA)|
|null n500 (no subgroup; homogeneous +26)|MR (IJ)|1987|0.994|21.827|1.209|1.963|0.977 (0.969, 0.983)|0.992 (0.987, 0.995)|
|null n500 (no subgroup; homogeneous +26)|MR (field)|1987|0.994|13.800|0.712|1.291|0.931 (0.919, 0.941)|0.967 (0.959, 0.974)|
|md40 n700|naive|1997|0.999|63.496|3.705|1.472|0.071 (0.061, 0.083)|0.168 (0.152, 0.185)|
|md40 n700|oracle|1997|0.999|0.100|0.006|0.968|0.941 (0.930, 0.950)|0.943 (0.932, 0.952)|
|md40 n700|MR (IJ)|1997|0.999|19.906|1.092|1.789|0.960 (0.951, 0.968)|0.985 (0.979, 0.989)|
|md40 n700|MR (field)|1997|0.999|11.646|0.596|1.243|0.930 (0.918, 0.941)|0.965 (0.957, 0.973)|

- **Field lower-bound coverage on β(Ĥ):** 0.937 / 0.965 / 0.931 / 0.930 (md40 n500 / md120 n500 / null / md40 n700).
  - The Wilson interval lies below 0.95 at md40 n500 (0.926–0.947), in the null cell (0.919–0.941) and at n = 700 (0.918–0.941).
  - It lies above 0.95 at md120 (0.955–0.972).
  - The oracle's coverage is 0.941–0.953.
- **Retained bias in SD units:** unadjusted +3.591 / +1.727 / +3.602 / +3.705; IJ +1.084 / −0.122 / +1.209 / +1.092; field +0.600 / −0.238 / +0.712 / +0.596.
- **SE/SD:** the field's is 1.010–1.291; IJ's is 1.419–1.963, with one-sided coverage 0.960–0.991.
- **Two-sided field interval:** covers at 0.971 / 0.909 / 0.967 / 0.965. At md120 it is below nominal while the one-sided lower bound is above it.

### Ĥᶜ — unadjusted, oracle, IJ two-term, field (beside), field-s (evaluated) (conditional on the proposed family)

|cell|estimator|n|declaration rate|bias (MD)|bias (SD units)|SE/SD|one-sided UPPER coverage (Wilson)|two-sided coverage (Wilson)|
|---|---|---|---|---|---|---|---|---|
|md40 n500|naive|1996|0.998|-16.964|-1.425|1.105|0.659 (0.638, 0.679)|0.778 (0.759, 0.795)|
|md40 n500|oracle|1996|0.998|-0.234|-0.016|1.007|0.958 (0.948, 0.966)|0.953 (0.943, 0.962)|
|md40 n500|MR (IJ)|1996|0.998|-6.728|-0.519|1.844|0.991 (0.986, 0.994)|1.000 (0.998, 1.000)|
|md40 n500|MR (field)|1996|0.998|-4.630|-0.348|0.970|0.894 (0.880, 0.907)|0.930 (0.918, 0.940)|
|md40 n500|MR (field-s)|1996|0.998|-4.638|-0.350|0.975|0.897 (0.883, 0.909)|0.933 (0.922, 0.943)|
|md120 n500|naive|2000|1.000|-12.253|-0.853|0.978|0.788 (0.770, 0.806)|0.868 (0.852, 0.882)|
|md120 n500|oracle|2000|1.000|-0.301|-0.021|1.002|0.956 (0.946, 0.964)|0.953 (0.942, 0.961)|
|md120 n500|MR (IJ)|2000|1.000|-2.197|-0.145|1.654|0.995 (0.990, 0.997)|0.999 (0.996, 1.000)|
|md120 n500|MR (field)|2000|1.000|-0.674|-0.043|0.878|0.934 (0.923, 0.945)|0.927 (0.915, 0.938)|
|md120 n500|MR (field-s)|2000|1.000|-0.704|-0.045|0.902|0.937 (0.925, 0.947)|0.936 (0.925, 0.946)|
|null n500 (no subgroup; homogeneous +26)|naive|1987|0.994|-16.580|-1.420|1.122|0.671 (0.650, 0.692)|0.788 (0.770, 0.806)|
|null n500 (no subgroup; homogeneous +26)|oracle|1987|0.994|-0.262|-0.023|1.028|0.955 (0.945, 0.963)|0.958 (0.948, 0.966)|
|null n500 (no subgroup; homogeneous +26)|MR (IJ)|1987|0.994|-6.824|-0.528|1.853|0.992 (0.988, 0.995)|0.999 (0.997, 1.000)|
|null n500 (no subgroup; homogeneous +26)|MR (field)|1987|0.994|-4.809|-0.362|0.969|0.890 (0.875, 0.903)|0.931 (0.919, 0.941)|
|null n500 (no subgroup; homogeneous +26)|MR (field-s)|1987|0.994|-4.817|-0.364|0.972|0.893 (0.879, 0.906)|0.933 (0.921, 0.943)|
|md40 n700|naive|1997|0.999|-11.642|-1.127|1.041|0.720 (0.700, 0.739)|0.822 (0.805, 0.838)|
|md40 n700|oracle|1997|0.999|0.201|0.017|1.006|0.955 (0.945, 0.964)|0.950 (0.940, 0.959)|
|md40 n700|MR (IJ)|1997|0.999|-4.207|-0.382|1.816|0.995 (0.991, 0.997)|0.998 (0.996, 0.999)|
|md40 n700|MR (field)|1997|0.999|-2.659|-0.237|0.943|0.910 (0.897, 0.922)|0.926 (0.914, 0.937)|
|md40 n700|MR (field-s)|1997|0.999|-2.664|-0.237|0.946|0.911 (0.898, 0.923)|0.928 (0.916, 0.938)|

- **Field-s upper-bound coverage on β(Ĥᶜ):** 0.897 / 0.937 / 0.893 / 0.911.
  - Its Wilson upper limits (0.909 / 0.947 / 0.906 / 0.923) are below 0.95 in every cell.
  - The unstudentized field covers at 0.894 / 0.934 / 0.890 / 0.910, so field-s adds 0.001–0.003.
- **Field-s retained bias:** −0.350 / −0.045 / −0.364 / −0.237 SD units, with SE/SD 0.902–0.975.
- **References:** IJ covers at 0.991–0.995 with SE/SD 1.654–1.853; the oracle at 0.955–0.958.

### Bound location on Ĥ and Ĥᶜ — the D3 ladder (conditional on the proposed family)

Ĥ: one-sided 95% LOWER bound (field | oracle; IJ and unadjusted for orientation).

|cell|estimator|n|mean|q05|q25|median|q75|q95|P(L>=0)|P(L>=10)|P(L>=20)|P(L>=30)|P(L>=40)|P(L>=60)|P(L>=80)|P(L>=100)|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|md40 n500|MR (field)|1996|1.2|-27.3|-12.2|0.4|12.8|33.3|0.508 (0.011)|0.301 (0.010)|0.157 (0.008)|0.070 (0.006)|0.029 (0.004)|0.004 (0.001)|0.001 (0.001)|0.000 (0.000)|
|md40 n500|oracle|1996|6.8|-27.0|-5.9|7.3|20.1|39.8|0.638 (0.011)|0.443 (0.011)|0.252 (0.010)|0.119 (0.007)|0.048 (0.005)|0.003 (0.001)|0.000 (0.000)|0.000 (0.000)|
|md40 n500|MR (IJ)|1996|-5.8|-34.4|-18.0|-6.2|5.8|23.6|0.370 (0.011)|0.184 (0.009)|0.076 (0.006)|0.030 (0.004)|0.011 (0.002)|0.001 (0.001)|0.000 (0.000)|0.000 (0.000)|
|md40 n500|naive|1996|50.8|24.5|39.8|50.7|61.3|77.5|0.996 (0.001)|0.990 (0.002)|0.971 (0.004)|0.905 (0.007)|0.747 (0.010)|0.279 (0.010)|0.039 (0.004)|0.004 (0.001)|
|md120 n500|MR (field)|2000|49.1|10.9|32.3|47.4|64.5|89.9|0.986 (0.003)|0.956 (0.005)|0.892 (0.007)|0.782 (0.009)|0.627 (0.011)|0.311 (0.010)|0.107 (0.007)|0.025 (0.003)|
|md120 n500|oracle|2000|86.7|53.0|74.1|87.3|100.1|119.8|1.000 (0.000)|1.000 (0.000)|0.999 (0.001)|0.997 (0.001)|0.986 (0.003)|0.907 (0.006)|0.638 (0.011)|0.252 (0.010)|
|md120 n500|MR (IJ)|2000|41.8|7.0|28.1|41.2|55.4|77.7|0.979 (0.003)|0.931 (0.006)|0.852 (0.008)|0.718 (0.010)|0.521 (0.011)|0.193 (0.009)|0.038 (0.004)|0.003 (0.001)|
|md120 n500|naive|2000|92.3|64.6|80.5|91.7|103.7|121.6|1.000 (0.000)|1.000 (0.000)|1.000 (0.000)|1.000 (0.000)|1.000 (0.000)|0.972 (0.004)|0.760 (0.010)|0.318 (0.010)|
|null n500 (no subgroup; homogeneous +26)|MR (field)|1987|-2.4|-32.1|-15.4|-3.2|9.7|29.8|0.431 (0.011)|0.243 (0.010)|0.120 (0.007)|0.050 (0.005)|0.017 (0.003)|0.003 (0.001)|0.001 (0.001)|0.000 (0.000)|
|null n500 (no subgroup; homogeneous +26)|oracle|0|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|NA|
|null n500 (no subgroup; homogeneous +26)|MR (IJ)|1987|-10.2|-38.6|-22.8|-10.3|1.5|19.8|0.278 (0.010)|0.133 (0.008)|0.050 (0.005)|0.018 (0.003)|0.005 (0.002)|0.001 (0.001)|0.000 (0.000)|0.000 (0.000)|
|null n500 (no subgroup; homogeneous +26)|naive|1987|46.3|18.7|35.5|46.3|57.4|73.5|0.993 (0.002)|0.979 (0.003)|0.943 (0.005)|0.836 (0.008)|0.656 (0.011)|0.200 (0.009)|0.021 (0.003)|0.002 (0.001)|
|md40 n700|MR (field)|1997|1.6|-26.3|-11.0|-0.4|13.0|34.6|0.494 (0.011)|0.290 (0.010)|0.162 (0.008)|0.075 (0.006)|0.035 (0.004)|0.008 (0.002)|0.001 (0.001)|0.000 (0.000)|
|md40 n700|oracle|1997|12.5|-15.4|0.8|12.1|24.0|41.5|0.764 (0.009)|0.553 (0.011)|0.322 (0.010)|0.165 (0.008)|0.059 (0.005)|0.004 (0.001)|0.000 (0.000)|0.000 (0.000)|
|md40 n700|MR (IJ)|1997|-2.2|-28.8|-13.8|-3.3|8.1|27.8|0.421 (0.011)|0.222 (0.009)|0.107 (0.007)|0.044 (0.005)|0.017 (0.003)|0.002 (0.001)|0.000 (0.000)|0.000 (0.000)|
|md40 n700|naive|1997|53.6|29.3|43.7|53.1|62.7|79.9|0.996 (0.001)|0.995 (0.002)|0.985 (0.003)|0.945 (0.005)|0.825 (0.009)|0.316 (0.010)|0.050 (0.005)|0.007 (0.002)|

Ĥᶜ: one-sided 95% UPPER bound (field-s | oracle; field, IJ and unadjusted for orientation).

|cell|estimator|n|mean|q05|q25|median|q75|q95|P(U<=0)|P(U<=10)|P(U<=20)|P(U<=30)|P(U<=40)|P(U<=60)|P(U<=80)|P(U<=100)|
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
|md40 n500|MR (field-s)|1996|47.4|25.8|39.0|47.7|56.2|69.2|0.000 (0.000)|0.006 (0.002)|0.024 (0.003)|0.094 (0.007)|0.271 (0.010)|0.830 (0.008)|0.996 (0.001)|1.000 (0.000)|
|md40 n500|oracle|1996|49.7|27.1|39.9|50.0|59.5|72.5|0.001 (0.001)|0.006 (0.002)|0.020 (0.003)|0.079 (0.006)|0.252 (0.010)|0.758 (0.010)|0.984 (0.003)|1.000 (0.000)|
|md40 n500|MR (field)|1996|47.4|25.1|38.9|47.7|56.6|69.5|0.000 (0.000)|0.007 (0.002)|0.026 (0.004)|0.096 (0.007)|0.275 (0.010)|0.830 (0.008)|0.996 (0.001)|1.000 (0.000)|
|md40 n500|MR (IJ)|1996|63.4|42.8|55.1|63.7|71.5|84.6|0.000 (0.000)|0.000 (0.000)|0.000 (0.000)|0.007 (0.002)|0.033 (0.004)|0.392 (0.011)|0.900 (0.007)|0.999 (0.001)|
|md40 n500|naive|1996|35.5|16.1|27.8|35.6|43.3|55.8|0.002 (0.001)|0.020 (0.003)|0.092 (0.006)|0.318 (0.010)|0.658 (0.011)|0.984 (0.003)|1.000 (0.000)|1.000 (0.000)|
|md120 n500|MR (field-s)|2000|63.6|39.7|53.3|64.0|74.0|88.0|0.000 (0.000)|0.001 (0.000)|0.006 (0.002)|0.015 (0.003)|0.053 (0.005)|0.398 (0.011)|0.868 (0.008)|0.993 (0.002)|
|md120 n500|oracle|2000|49.7|26.8|39.8|49.9|59.5|72.5|0.001 (0.001)|0.006 (0.002)|0.021 (0.003)|0.081 (0.006)|0.253 (0.010)|0.758 (0.010)|0.985 (0.003)|1.000 (0.000)|
|md120 n500|MR (field)|2000|63.1|38.4|52.8|63.5|73.9|88.0|0.000 (0.000)|0.001 (0.000)|0.008 (0.002)|0.018 (0.003)|0.065 (0.005)|0.413 (0.011)|0.871 (0.007)|0.993 (0.002)|
|md120 n500|MR (IJ)|2000|80.3|57.0|70.5|80.5|90.1|104.1|0.000 (0.000)|0.000 (0.000)|0.000 (0.000)|0.001 (0.000)|0.008 (0.002)|0.082 (0.006)|0.486 (0.011)|0.918 (0.006)|
|md120 n500|naive|2000|52.2|29.6|43.0|52.4|61.7|74.3|0.000 (0.000)|0.005 (0.002)|0.015 (0.003)|0.054 (0.005)|0.186 (0.009)|0.711 (0.010)|0.981 (0.003)|1.000 (0.000)|
|null n500 (no subgroup; homogeneous +26)|MR (field-s)|1987|42.6|21.2|34.1|42.6|51.3|64.7|0.001 (0.001)|0.010 (0.002)|0.045 (0.005)|0.170 (0.008)|0.423 (0.011)|0.903 (0.007)|0.999 (0.001)|1.000 (0.000)|
|null n500 (no subgroup; homogeneous +26)|oracle|1987|45.2|26.8|37.5|45.2|52.6|64.0|0.000 (0.000)|0.002 (0.001)|0.014 (0.003)|0.088 (0.006)|0.330 (0.011)|0.896 (0.007)|0.999 (0.001)|1.000 (0.000)|
|null n500 (no subgroup; homogeneous +26)|MR (field)|1987|42.6|21.1|34.0|42.6|51.6|64.9|0.001 (0.001)|0.011 (0.002)|0.048 (0.005)|0.170 (0.008)|0.426 (0.011)|0.900 (0.007)|0.999 (0.001)|1.000 (0.000)|
|null n500 (no subgroup; homogeneous +26)|MR (IJ)|1987|58.8|38.5|50.4|58.8|67.0|80.5|0.000 (0.000)|0.000 (0.000)|0.002 (0.001)|0.014 (0.003)|0.063 (0.005)|0.543 (0.011)|0.947 (0.005)|1.000 (0.000)|
|null n500 (no subgroup; homogeneous +26)|naive|1987|31.2|12.5|23.6|31.1|38.9|51.4|0.007 (0.002)|0.030 (0.004)|0.159 (0.008)|0.458 (0.011)|0.783 (0.009)|0.994 (0.002)|1.000 (0.000)|1.000 (0.000)|
|md40 n700|MR (field-s)|1997|45.6|27.0|38.0|45.8|53.3|64.2|0.000 (0.000)|0.002 (0.001)|0.013 (0.002)|0.079 (0.006)|0.303 (0.010)|0.904 (0.007)|0.998 (0.001)|1.000 (0.000)|
|md40 n700|oracle|1997|46.5|27.3|38.5|46.6|54.6|66.5|0.000 (0.000)|0.002 (0.001)|0.012 (0.002)|0.085 (0.006)|0.294 (0.010)|0.872 (0.007)|0.996 (0.001)|1.000 (0.000)|
|md40 n700|MR (field)|1997|45.6|27.0|37.9|45.7|53.3|64.1|0.000 (0.000)|0.002 (0.001)|0.013 (0.003)|0.079 (0.006)|0.307 (0.010)|0.904 (0.007)|0.998 (0.001)|1.000 (0.000)|
|md40 n700|MR (IJ)|1997|59.5|41.9|52.2|59.4|66.9|77.9|0.000 (0.000)|0.000 (0.000)|0.001 (0.001)|0.005 (0.001)|0.036 (0.004)|0.520 (0.011)|0.970 (0.004)|1.000 (0.000)|
|md40 n700|naive|1997|36.9|20.5|30.0|36.8|43.8|54.1|0.001 (0.001)|0.005 (0.002)|0.047 (0.005)|0.251 (0.010)|0.624 (0.011)|0.989 (0.002)|1.000 (0.000)|1.000 (0.000)|

- **Field lower bound on Ĥ, share of declared replicates:**
  - at or above 0: 0.508 / 0.986 / 0.431 / 0.494;
  - at or above 40: 0.029 / 0.627 / 0.017 / 0.035.
  - The oracle's bound sits at or above 40 on 0.048 / 0.986 / – / 0.059.
  - At md120 the field bound reaches 100 on 0.025 of replicates, against the oracle's 0.252.
- **Null cell** (truth +26 everywhere): the field lower bound sits at or above 30 on 0.050 of replicates; the unadjusted bound on 0.836.
- **Field-s upper bound on Ĥᶜ, share of declared replicates:**
  - at or below 30: 0.094 / 0.015 / 0.170 / 0.079 (oracle 0.079 / 0.081 / 0.088 / 0.085);
  - at or below 60: 0.830 / 0.398 / 0.903 / 0.904 (oracle 0.758 / 0.758 / 0.896 / 0.872).

### Joint pair (Ĥ lower, Ĥᶜ upper) (conditional on the proposed family)

|cell|pair|declared|both_bounds|share_both|joint|joint_mc_se|cov_H|cov_Hc|margin_H|margin_Hc|
|---|---|---|---|---|---|---|---|---|---|---|
|md40 n500|Bonferroni field-s (gamma = 0.025)|1996|1996|1|0.914 (0.901, 0.926)|0.006|0.972|0.941|59.103|27.392|
|md40 n500|Bonferroni unstudentized (gamma = 0.025)|1996|1996|1|0.911 (0.898, 0.923)|0.006|0.972|0.937|59.103|27.390|
|md40 n500|separate 95% bounds: field lower, field-s upper|1996|1996|1|0.840 (0.823, 0.856)|0.008|0.937|0.897|49.926|23.346|
|md120 n500|Bonferroni field-s (gamma = 0.025)|2000|2000|1|0.950 (0.940, 0.959)|0.005|0.983|0.967|54.754|28.896|
|md120 n500|Bonferroni unstudentized (gamma = 0.025)|2000|2000|1|0.945 (0.934, 0.954)|0.005|0.983|0.962|54.754|28.267|
|md120 n500|separate 95% bounds: field lower, field-s upper|2000|2000|1|0.903 (0.890, 0.916)|0.007|0.965|0.937|45.698|24.523|
|null n500 (no subgroup; homogeneous +26)|Bonferroni field-s (gamma = 0.025)|1987|1987|1|0.910 (0.897, 0.922)|0.006|0.969|0.941|59.719|27.218|
|null n500 (no subgroup; homogeneous +26)|Bonferroni unstudentized (gamma = 0.025)|1987|1987|1|0.909 (0.895, 0.921)|0.006|0.969|0.939|59.719|27.251|
|null n500 (no subgroup; homogeneous +26)|separate 95% bounds: field lower, field-s upper|1987|1987|1|0.832 (0.815, 0.848)|0.008|0.931|0.893|50.475|23.165|
|md40 n700|Bonferroni field-s (gamma = 0.025)|1997|1997|1|0.911 (0.898, 0.923)|0.006|0.966|0.944|59.135|22.295|
|md40 n700|Bonferroni unstudentized (gamma = 0.025)|1997|1997|1|0.909 (0.896, 0.921)|0.006|0.966|0.943|59.135|22.296|
|md40 n700|separate 95% bounds: field lower, field-s upper|1997|1997|1|0.846 (0.829, 0.861)|0.008|0.930|0.911|49.855|18.989|

|cell|metric|value|mc_se|n|
|---|---|---|---|---|
|md40 n500|gamma_mean_s|0.0252|0.0000|1996|
|md40 n500|corr_s|0.0879|0.0013|1996|
|md40 n500|gamma_mean|0.0251|0.0000|1996|
|md40 n500|corr|0.0905|0.0014|1996|
|md120 n500|gamma_mean_s|0.0252|0.0000|2000|
|md120 n500|corr_s|0.0111|0.0017|2000|
|md120 n500|gamma_mean|0.0252|0.0000|2000|
|md120 n500|corr|0.0148|0.0017|2000|
|null n500 (no subgroup; homogeneous +26)|gamma_mean_s|0.0252|0.0000|1987|
|null n500 (no subgroup; homogeneous +26)|corr_s|0.0855|0.0013|1987|
|null n500 (no subgroup; homogeneous +26)|gamma_mean|0.0251|0.0000|1987|
|null n500 (no subgroup; homogeneous +26)|corr|0.0878|0.0013|1987|
|md40 n700|gamma_mean_s|0.0252|0.0000|1997|
|md40 n700|corr_s|0.0879|0.0012|1997|
|md40 n700|gamma_mean|0.0251|0.0000|1997|
|md40 n700|corr|0.0901|0.0013|1997|

- **Field-s Bonferroni pair:** covers (β(Ĥ), β(Ĥᶜ)) jointly on 0.914 / 0.950 / 0.910 / 0.911 of declared replicates.
  - Its Wilson upper limits are 0.926 / 0.959 / 0.922 / 0.923.
  - Both bounds exist on every declared replicate.
- **Beside it:** the unstudentized pair is within 0.005 of it; the separate 95% pair covers at 0.832–0.903.
- **Calibrated split:** returns the Bonferroni floor (mean γ 0.0251–0.0252), with corr(Λ*, Λ*ᶜ) 0.0111–0.0879 (field-s).

### Identification, DINA beside GRF and FS (conditional on the proposed family)

|cell|quantity|DINA|GRF_mdgrf|FS_mdsgnb20|
|---|---|---|---|---|
|md40 n500|declaration rate|0.9980|1.0000|0.9990|
|md40 n500|mean size of Hhat (n_sel)|108.53|109.46|111.62|
|md40 n500|mean true positives (sens x n_true)|42.39|-|-|
|md40 n500|sensitivity (mean over declared; NA where n_true = 0)|0.2462|0.2497|0.2679|
|md40 n500|PPV (mean over declared)|0.3888|0.3909|0.4122|
|md40 n500|mean proposed family (dina_proposed_n)|3341.7|-|-|
|md40 n500|mean family size K (n_family, MR’s kept family)|3341.7|1222.4|-|
|md40 n500|mean admitted set (admitted_n)|2631.0|624.6|-|
|md40 n500|size of Hhat DINA larger / equal / smaller than FS (paired by sim_id, both declared)|852 / 93 / 1050|-|-|
|md120 n500|declaration rate|1.0000|1.0000|1.0000|
|md120 n500|mean size of Hhat (n_sel)|146.59|143.61|152.74|
|md120 n500|mean true positives (sens x n_true)|113.50|-|-|
|md120 n500|sensitivity (mean over declared; NA where n_true = 0)|0.6596|0.5989|0.7116|
|md120 n500|PPV (mean over declared)|0.7666|0.7149|0.7980|
|md120 n500|mean proposed family (dina_proposed_n)|5940.0|-|-|
|md120 n500|mean family size K (n_family, MR’s kept family)|5940.0|1222.4|-|
|md120 n500|mean admitted set (admitted_n)|5599.1|1059.8|-|
|md120 n500|size of Hhat DINA larger / equal / smaller than FS (paired by sim_id, both declared)|757 / 78 / 1165|-|-|
|null n500 (no subgroup; homogeneous +26)|declaration rate|0.9935|1.0000|0.9965|
|null n500 (no subgroup; homogeneous +26)|mean size of Hhat (n_sel)|105.15|106.14|107.49|
|null n500 (no subgroup; homogeneous +26)|mean true positives (sens x n_true)|NA|-|-|
|null n500 (no subgroup; homogeneous +26)|sensitivity (mean over declared; NA where n_true = 0)|NA|NA|NA|
|null n500 (no subgroup; homogeneous +26)|PPV (mean over declared)|0.0000|0.0000|0.0000|
|null n500 (no subgroup; homogeneous +26)|mean proposed family (dina_proposed_n)|2644.1|-|-|
|null n500 (no subgroup; homogeneous +26)|mean family size K (n_family, MR’s kept family)|2644.1|1222.4|-|
|null n500 (no subgroup; homogeneous +26)|mean admitted set (admitted_n)|1959.0|508.5|-|
|null n500 (no subgroup; homogeneous +26)|size of Hhat DINA larger / equal / smaller than FS (paired by sim_id, both declared)|867 / 104 / 1010|-|-|
|md40 n700|declaration rate|0.9985|1.0000|0.9995|
|md40 n700|mean size of Hhat (n_sel)|113.84|116.33|117.99|
|md40 n700|mean true positives (sens x n_true)|44.92|-|-|
|md40 n700|sensitivity (mean over declared; NA where n_true = 0)|0.1862|0.1873|0.2041|
|md40 n700|PPV (mean over declared)|0.3895|0.3863|0.4101|
|md40 n700|mean proposed family (dina_proposed_n)|4007.0|-|-|
|md40 n700|mean family size K (n_family, MR’s kept family)|4007.0|1354.5|-|
|md40 n700|mean admitted set (admitted_n)|3202.8|711.4|-|
|md40 n700|size of Hhat DINA larger / equal / smaller than FS (paired by sim_id, both declared)|873 / 80 / 1043|-|-|

- **Mean |Ĥ|:** 108.53 / 146.59 / 105.15 / 113.84 patients for DINA, against GRF's 109.46 / 143.61 / 106.14 / 116.33 and FS's 111.62 / 152.74 / 107.49 / 117.99.
- **Paired by `sim_id` with FS:** DINA's Ĥ is larger on 852 / 757 / 867 / 873, equal on 93 / 78 / 104 / 80, and smaller on 1,050 / 1,165 / 1,010 / 1,043 replicates.
- **Truly harmed patients in Ĥ:** 42.39 / 113.50 / – / 44.92 on average (sensitivity × `n_true`; never |Ĥ|).
- **Sensitivity:** 0.2462 / 0.6596 / – / 0.1862 (GRF 0.2497 / 0.5989 / – / 0.1873; FS 0.2679 / 0.7116 / – / 0.2041).
- **PPV:** 0.3888 / 0.7666 / 0 / 0.3895 (GRF 0.3909 / 0.7149 / 0 / 0.3863; FS 0.4122 / 0.7980 / 0 / 0.4101). In the null cell PPV is 0 by construction and sensitivity is undefined.
- **Family sizes:**
  - DINA's proposed family averages 3,341.7 / 5,940.0 / 2,644.1 / 4,007.0 candidates, and MR re-selects over all of it (`n_family` equals the proposed count).
  - The admitted set averages 2,631.0 / 5,599.1 / 1,959.0 / 3,202.8.
  - GRF's family is 1,222.4 / 1,222.4 / 1,222.4 / 1,354.5, and its admitted set 624.6 / 1,059.8 / 508.5 / 711.4.
  - DINA's family depends on the outcome through the fitted surface, and it is largest in the md120 cell.

### The three identifiers (conditional on the proposed family)

Rows are the cells, one line per identifier. FS and GRF values are copied from `md_field_metrics.csv` and `md_grf_metrics.csv`. FS's family is the prespecified cut grid; the GRF and DINA figures are conditional on their proposed families.

|cell|identifier|declaration rate|field lower coverage, Hhat|field-s upper coverage, Hhat^c|Bonferroni joint (field-s)|mean size of Hhat|sensitivity|PPV|
|---|---|---|---|---|---|---|---|---|
|md40 n500|FS (mdsgnb20)|0.9990|0.9479|0.9139|0.9279|111.6241|0.2679|0.4122|
|md40 n500|GRF (mdgrf)|1.0000|0.9390|0.9200|0.9195|109.4590|0.2497|0.3909|
|md40 n500|DINA (mddina)|0.9980|0.9374|0.8968|0.9143|108.5326|0.2462|0.3888|
|md120 n500|FS (mdsgnb20)|1.0000|0.9670|0.9405|0.9490|152.7385|0.7116|0.7980|
|md120 n500|GRF (mdgrf)|1.0000|0.9560|0.9315|0.9425|143.6125|0.5989|0.7149|
|md120 n500|DINA (mddina)|1.0000|0.9645|0.9370|0.9500|146.5870|0.6596|0.7666|
|null n500 (no subgroup; homogeneous +26)|FS (mdsgnb20)|0.9965|0.9473|0.9172|0.9318|107.4922|NA|0.0000|
|null n500 (no subgroup; homogeneous +26)|GRF (mdgrf)|1.0000|0.9405|0.9200|0.9300|106.1385|NA|0.0000|
|null n500 (no subgroup; homogeneous +26)|DINA (mddina)|0.9935|0.9311|0.8933|0.9104|105.1500|NA|0.0000|
|md40 n700|FS (mdsgnb20)|0.9995|0.9395|0.9280|0.9275|117.9905|0.2041|0.4101|
|md40 n700|GRF (mdgrf)|1.0000|0.9430|0.9310|0.9345|116.3330|0.1873|0.3863|
|md40 n700|DINA (mddina)|0.9985|0.9304|0.9114|0.9109|113.8353|0.1862|0.3895|

**Reading.**
- The confound applies throughout: each identifier's figures condition on its own declared replicates, and on a family built differently (a prespecified grid for FS, fitted surfaces for GRF and DINA). The rows are therefore read side by side, not ranked.
- **Declaration:** the three identifiers agree to within 0.007 in every cell (0.9935–1.0000).
- **Field lower-bound coverage on Ĥ:** DINA 0.9374 / 0.9645 / 0.9311 / 0.9304, GRF 0.9390 / 0.9560 / 0.9405 / 0.9430, FS 0.9479 / 0.9670 / 0.9473 / 0.9395. DINA's is within 0.017 of GRF's and FS's in every cell, and the lowest of the three in the md40 n500, null and n = 700 cells.
- **Field-s upper-bound coverage on Ĥᶜ:** DINA 0.8968 / 0.9370 / 0.8933 / 0.9114, the lowest of the three in the md40 n500, null and n = 700 cells (GRF 0.9200 / 0.9315 / 0.9200 / 0.9310; FS 0.9139 / 0.9405 / 0.9172 / 0.9280). At md120, DINA's 0.9370 sits between GRF's and FS's.
- **Field-s Bonferroni joint coverage:** DINA 0.9143 / 0.9500 / 0.9104 / 0.9109; GRF 0.9195 / 0.9425 / 0.9300 / 0.9345; FS 0.9279 / 0.9490 / 0.9318 / 0.9275.
- **Mean |Ĥ|:** DINA's is the smallest of the three in every cell except md120, where it lies between GRF's and FS's. At md120, DINA's sensitivity (0.6596) and PPV (0.7666) lie between GRF's and FS's; in the two MD 40 cells all three are within 0.022 on sensitivity and 0.024 on PPV.

### Regime diagnostics (conditional on the proposed family)

|cell|n_det|p_hat_mean|p_hat_lt05|sd_btc_naive|lamc_naive|lamc_s_naive|ij_sd_H|ij_sd_Hc|fit_secs|field_secs|comp_secs|
|---|---|---|---|---|---|---|---|---|---|---|---|
|md40 n500|1996|0.080|0.996|0.985|0.979|0.982|1.919|1.844|134.235|69.798|6.999|
|md120 n500|2000|0.083|0.999|1.076|0.972|0.998|1.419|1.654|304.939|149.201|17.650|
|null n500 (no subgroup; homogeneous +26)|1987|0.092|0.990|0.986|0.981|0.983|1.963|1.853|100.276|54.362|4.951|
|md40 n700|1997|0.076|0.997|1.024|0.985|0.987|1.789|1.816|199.655|88.036|12.555|

- **p̂(Ĥ):** 0.080 / 0.083 / 0.092 / 0.076, with 0.990–0.999 of replicates below 0.5 — the tie regime.
- **Complement scale:** SD(β̃ᶜ) / mean naive SEᶜ is 0.985–1.076, and λ-SDᶜ / naive SEᶜ is 0.972–0.985 (field) and 0.982–0.998 (field-s); the complement is not in the moved regime.
- **IJ SE / empirical SD:** 1.419–1.963 on Ĥ, 1.654–1.853 on Ĥᶜ.
- **Seconds per replicate:** 100.276–304.939 (field 54.362–149.201; complement 4.951–17.650).

### The display (identity scale) (conditional on the proposed family)

|cell|block|estimator|n|bias (MD)|SD|mean SE|b|r|1-sided cov|1-sided ref|2-sided cov|2-sided ref|
|---|---|---|---|---|---|---|---|---|---|---|---|---|
|md40 n500|H|naive|1996|61.464|17.117|25.674|3.591|1.500|0.111|0.131|0.244|0.257|
|md40 n500|H|mr|1996|19.566|18.048|34.641|1.084|1.919|0.976|0.981|0.993|0.996|
|md40 n500|H|fld|1996|11.615|19.362|24.750|0.600|1.278|0.937|0.934|0.971|0.971|
|md40 n500|Hc|naive|1996|-16.964|11.907|13.160|-1.425|1.105|0.659|0.653|0.778|0.771|
|md40 n500|Hc|mr|1996|-6.728|12.960|23.904|-0.519|1.844|0.991|0.994|1.000|0.999|
|md40 n500|Hc|fld|1996|-4.630|13.289|12.889|-0.348|0.970|0.894|0.894|0.930|0.927|
|md40 n500|Hc|fld_s|1996|-4.638|13.263|12.927|-0.350|0.975|0.897|0.895|0.933|0.929|
|md120 n500|H|naive|2000|31.462|18.216|22.308|1.727|1.225|0.596|0.613|0.705|0.750|
|md120 n500|H|mr|2000|-2.764|22.680|32.192|-0.122|1.419|0.991|0.993|0.984|0.994|
|md120 n500|H|fld|2000|-5.940|24.997|25.239|-0.238|1.010|0.965|0.971|0.909|0.946|
|md120 n500|Hc|naive|2000|-12.253|14.367|14.051|-0.853|0.978|0.788|0.775|0.868|0.854|
|md120 n500|Hc|mr|2000|-2.197|15.122|25.016|-0.145|1.654|0.995|0.995|0.999|0.999|
|md120 n500|Hc|fld|2000|-0.674|15.555|13.657|-0.043|0.878|0.934|0.919|0.927|0.914|
|md120 n500|Hc|fld_s|2000|-0.704|15.549|14.019|-0.045|0.902|0.937|0.925|0.936|0.922|
|null n500 (no subgroup; homogeneous +26)|H|naive|1987|62.973|17.485|26.102|3.602|1.493|0.115|0.126|0.233|0.250|
|null n500 (no subgroup; homogeneous +26)|H|mr|1987|21.827|18.048|35.421|1.209|1.963|0.977|0.978|0.992|0.996|
|null n500 (no subgroup; homogeneous +26)|H|fld|1987|13.800|19.374|25.012|0.712|1.291|0.931|0.921|0.967|0.965|
|null n500 (no subgroup; homogeneous +26)|Hc|naive|1987|-16.580|11.679|13.101|-1.420|1.122|0.671|0.665|0.788|0.782|
|null n500 (no subgroup; homogeneous +26)|Hc|mr|1987|-6.824|12.914|23.930|-0.528|1.853|0.992|0.994|0.999|0.999|
|null n500 (no subgroup; homogeneous +26)|Hc|fld|1987|-4.809|13.267|12.858|-0.362|0.969|0.890|0.891|0.931|0.926|
|null n500 (no subgroup; homogeneous +26)|Hc|fld_s|1987|-4.817|13.244|12.873|-0.364|0.972|0.893|0.892|0.933|0.927|
|md40 n700|H|naive|1997|63.496|17.137|25.225|3.705|1.472|0.071|0.100|0.168|0.206|
|md40 n700|H|mr|1997|19.906|18.235|32.616|1.092|1.789|0.960|0.968|0.985|0.992|
|md40 n700|H|fld|1997|11.646|19.546|24.288|0.596|1.243|0.930|0.926|0.965|0.966|
|md40 n700|Hc|naive|1997|-11.642|10.330|10.753|-1.127|1.041|0.720|0.721|0.822|0.819|
|md40 n700|Hc|mr|1997|-4.207|11.010|19.997|-0.382|1.816|0.995|0.995|0.998|0.999|
|md40 n700|Hc|fld|1997|-2.659|11.229|10.586|-0.237|0.943|0.910|0.906|0.926|0.928|
|md40 n700|Hc|fld_s|1997|-2.664|11.219|10.617|-0.237|0.946|0.911|0.906|0.928|0.929|

- **Field on Ĥ against the Gaussian reference:** 0.937 vs 0.934, 0.965 vs 0.971, 0.931 vs 0.921, 0.930 vs 0.926 (b −0.238 to +0.712, r 1.010–1.291). The md120 two-sided point (0.909 vs 0.946) is the largest departure on this block.
- **Field-s on Ĥᶜ:** 0.897 vs 0.895, 0.937 vs 0.925, 0.893 vs 0.892, 0.911 vs 0.906 (b −0.364 to −0.045, r 0.902–0.975).
- **Figures:** `fig_mddina_bias_coverage_display_H.png` and `fig_mddina_bias_coverage_display_Hc.png`, also embedded in the summary HTML, captioned with the conditional-family sentence.

## The confound

DINA's and GRF's candidate families are generated from fitted surfaces, so the fixed-family condition does not hold for them: every DINA and GRF coverage figure below is coverage of the estimand conditional on the proposed family. FS's family is the prespecified cut grid. The three identifiers are compared descriptively, not as a contest: the identifier, the family construction and the set of detected replicates each summary conditions on all differ.

## Scope

These are operating characteristics of the DINA identifier on one continuous design, conditional on the proposed family. They do not verify condition (A3), no construction is promoted on this design's performance, and DINA ran with the proposal-floor orientation fix landed in the commits recorded in REPORT_grf_dina_fixes_2026-09-16.md.

## Findings

1. **P2 in operation.**
   - The oriented proposal floor holds on every declared replicate of every cell (smallest proposed oriented tau-hat 30.0000).
   - The admitted set is never empty on a declared replicate; its minimum is 1 in three cells and 397 at md120.
   - DINA declares on 99.35–100% of replicates.
   - Before P2, DINA on this design proposed benefit candidates and selected nothing (Stage 0, S0.3).
2. **Cost, and the projection.**
   - 27,727 s for 8,000 replicates at 63 workers, 26% over the 21,959-s projection and within the ceiling.
   - The projection scaled per-cell cost by FS's cell ratios, but DINA's cost follows its family size. That size depends on the outcome and is largest at md120 (mean 5,940 candidates), where the cell took 10,334 s against the projected 5,664 s.
   - Per replicate, DINA with MR costs 100–305 s, against GRF's 44–53 s and FS's 54–85 s. Most of it is the field pass over the large family (54–149 s).
   - Peak memory reached 123.3 GB.
3. **MR's family is DINA's whole proposed family.** `n_family` equals `dina_proposed_n` on every declared replicate. Its minimum is 2 / 415 / 1 / 1 and its maximum 6,683 / 6,712 / 6,617 / 7,706 candidates, so on some replicates the family has a single member.
4. **Non-detections.**
   - 4 / 0 / 13 / 3 replicates are undeclared (FS: 2 / 0 / 7 / 1).
   - The E2 fields are NA on them by construction (DINA's selection object is dropped on a non-detection), so the bundle cannot say whether the proposal or the admission came up empty (Gate 1 record F3).
5. **Complement one-sided coverage is below nominal** for field-s in every cell (Wilson upper limits 0.906–0.947), conditional on the proposed family, and below GRF's and FS's in three of the four cells. Reported, not interpreted.
6. **Convention 9.**
   - The recorder's `n_harm` is |Ĥ| (template `:794`).
   - The FS extract's `mean_n_harm` equals its `mean_n_sel`.
   - The DINA extract reports |Ĥ| as `mean_size_hhat` and the true-positive count as sensitivity × `n_true`. It does not copy `mean_n_harm`, and COLUMNS states the equality.
7. **The use_dina screening path** under FS still applies the unoriented floor (fix record F3). It was not used here and is out of scope.
8. **Recorder, authored fields.** `dina_searched_n`, `dina_proposed_n`, `dina_tau_min` and the DINA fill of `admitted_n` have no survival precedent; they carry the P2 checks per replicate (Gate 1 record F3).
9. **Record tables.** The tables were taken from the rendered HTML (xml2, in order, with captions), and labels are pipe-free ("size of Hhat"). kable's row-name column is dropped.
10. **Raw campaign logs** stay untracked under `logs_mddina/`: the runner log, per-render logs and peaks, `GATE2_*.txt`, and batch HTML. The closeout deletes the smoke, calibration, check and dry-run items only.
11. **Process check (Gate 1 record F1).** An orphaned memory-sampling loop from another session, running since 2026-09-06, matched the task's `ps` pattern; it runs no R.

## Commits of this campaign (`git log --oneline`, oldest last; the catalog commits follow this record)

```
34067b7d mddina extract (TASK_md_dina_campaign_2026-09-17 §3.2): md_dina_metrics.csv, 1,174 rows (1,110 dina + 28 FS rows copied from md_fi
1ffef474 mddina summary (TASK_md_dina_campaign_2026-09-17 §3.1): summary_continuous_field_mddina.qmd, the transplant of summary_continuous_
30b282eb mddina: campaign complete (cumulative render wall 27727 s)
21e48ee5 mddina md40_n700: md 40 n 700, 2000 replicates + combine; Gate 2 PASS (GATE_COUNTS run=65 passed=65 failed=0)
9051b702 mddina null_n500: md null n 500, 2000 replicates + combine; Gate 2 PASS (GATE_COUNTS run=65 passed=65 failed=0)
c9f1b948 mddina md120_n500: md 120 n 500, 2000 replicates + combine; Gate 2 PASS (GATE_COUNTS run=65 passed=65 failed=0)
1856c107 mddina md40_n500: md 40 n 500, 2000 replicates + combine; Gate 2 PASS (GATE_COUNTS run=65 passed=65 failed=0)
37395393 DINA on the ACTG175 continuous (MD) design, campaign mddina, Stage 1 record (TASK_md_dina_campaign_2026-09-17 §1.8): gates 1.1-1.3
d65f3b9b scripts_mddina (TASK_md_dina_campaign_2026-09-17 §1.6, transplants of scripts_mdgrf): mem_sampler.sh unchanged; run_mddina.sh with
7e26ad42 actg175/continuous catalog generator (TASK_md_dina_campaign_2026-09-17 §1.5, carried fix): current_status_regen.R gains campaign r
0927669c MD template E1-E3 (TASK_md_dina_campaign_2026-09-17 §1.5): E1 no edit -- the DINA argument block (dina_args = list(), dina_select_
189d4ec2 dev/tasks: TASK_md_dina_campaign_2026-09-17.md as received (campaign mddina, Stages 1-3)
<this record; then status_curated.md; then current_status.md alone>
```
