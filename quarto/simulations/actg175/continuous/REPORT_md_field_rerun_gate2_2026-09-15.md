# REPORT — mdsgnb20 Gate 2 (per cell)

Task: `dev/tasks/TASK_md_field_rerun_2026-09-15.md` §2.3. Runner: `quarto/simulations/actg175/continuous/scripts_mdsgnb20/run_mdsgnb20.sh`; checker: `quarto/simulations/actg175/continuous/scripts_mdsgnb20/gate2.R`.

## md40_n500 — md 40, n 500

- Stem: `fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdsgnb20`; HEAD before this cell's commit: 34580fb0; workers 63; threads 1.
- Knobs: `FS_MD_FOCUS=effMaxSG FS_MD_NBHD=0.20 FS_MD_FIELD_SCALEC=selected FS_MD_CI=field FS_MD_CAMPAIGN=mdsgnb20 FS_MD_FB=none FS_MD_MD=40 FS_MD_N=500 FS_MD_WORKERS=63`.
- Seeds: 8316951 + sim_id; batches sim_id 1-1000 and 1001-2000, then combine.
- Cell wall (renders, from the progress log): 1045 + 1035 + 41 = 2121 s; cumulative render wall 2121 s (ceiling 19071 s).
- Render: WALL_SECONDS=1045 RC=0 PEAK_MB=72407 OUT=fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdsgnb20_batch_1_1000.html
- Render: WALL_SECONDS=41 RC=0 PEAK_MB=1207 OUT=fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdsgnb20_combine_1_2000.html
- Render: WALL_SECONDS=1035 RC=0 PEAK_MB=73284 OUT=fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdsgnb20_batch_1001_2000.html
- **Finding (checker, not data):** the runner's first Gate 2 pass on this cell halted at 07:03:42Z with `GATE_COUNTS run=53 passed=52 failed=1`: every data check had passed and the one failure was the checker itself erroring (`too few arguments`: the same-draws oracle label's `sprintf` lacked its tolerance argument, `gate2.R:47`). The checker was fixed and committed by explicit path (`34580fb0`; the check is unchanged), re-run on the completed cell output without re-rendering, and passed 58 of 58 (below). The cell's commit and the halt file's removal were then done by hand with the runner's own paths, and the runner relaunched (it skips this cell as tracked and clean). Halt commit: `414467b3`.
- Gate counts: `GATE_COUNTS run=58 passed=58 failed=0`.

```

########## GATE 2 (mdsgnb20): md40_n500 -- md 40, n 500 ##########
  bundle: ../mr_md_harm/fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdsgnb20_d5000/fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdsgnb20_combined_1_2000.rds
  combined payload on disk                                       PASS     
  --- combine assertions ---
  exactly 2 batch files, res_1_1000 and res_1001_2000            PASS     (2)
  2,000 rows                                                     PASS     (2000)
  combined sim_id == 1:2000                                      PASS     
  batch sim_id sets 1:1000 and 1001:2000                         PASS     (batch1 1000, batch2 1000)
  batch files match the combined bundle on every column          PASS     (156 columns)
  no CONFIG-ERROR replicate                                      PASS     (0)
  --- meta (both batches) ---
  meta: subgroup_method == consistency                           PASS     (consistency / consistency)
  meta: sg_focus == effMaxSG                                     PASS     (effMaxSG / effMaxSG)
  meta: effect_neighborhood == 0.2                               PASS     (0.2 / 0.2)
  meta: selection_rule == neighborhood                           PASS     (neighborhood / neighborhood)
  meta: consistency_method == resample                           PASS     (resample / resample)
  meta: ci_method == field                                       PASS     (field / field)
  meta: mr_draws == 5000                                         PASS     (5000 / 5000)
  meta: field_uniform == FALSE                                   PASS     (FALSE / FALSE)
  meta: field_complement == TRUE                                 PASS     (TRUE / TRUE)
  meta: field_scale_complement == selected                       PASS     (selected / selected)
  meta: ij_residual == two_term                                  PASS     (two_term / two_term)
  meta: return_reselection == TRUE                               PASS     (TRUE / TRUE)
  meta: fb_mode == none                                          PASS     (none / none)
  meta: seed_base == 8316951                                     PASS     (8316951 / 8316951)
  meta: campaign_tag == mdsgnb20                                 PASS     (mdsgnb20 / mdsgnb20)
  meta: n_sample == 500                                          PASS     (500 / 500)
  meta: null_cell == FALSE                                       PASS     (FALSE / FALSE)
  meta: effect_threshold == 30                                   PASS     (30 / 30)
  meta: consistency_threshold == 10                              PASS     (10 / 10)
  meta: pkg_version == 0.3.5                                     PASS     (0.3.5 / 0.3.5)
  meta: hostname == pop-os                                       PASS     (pop-os / pop-os)
  meta: n_workers == 63                                          PASS     (63 / 63)
  meta seed_base 8316951 | host pop-os | R 4.6.1 | built_at 2026-09-16 00:02:42 / 2026-09-15 23:45:27
  >> DECLARATION RATE         : 0.9990 (1998 / 2000)
  MR failures on declared replicates <= max(20, 2 x mdf1's 0)    PASS     (0)
  every finiteness column present                                PASS     
  harm products finite on declared replicates with a field block PASS     (1998 rows)
  nine fld_Hc_*_s and nine fld_joint_s_* columns present         PASS     
  fld_Hc_*_s finite on all 1998 filled replicates                PASS     
  fld_joint_s_* finite on all 1998 filled replicates             PASS     
  complement field filled on 1998 of 1998 declared replicates (0 notes)
  invariant harm  : fld_H_lo1s <= fld_H_est2                     PASS     
  invariant compl : fld_Hc_est2 <= fld_Hc_up1s                   PASS     
  invariant _s    : fld_Hc_lo1s_s <= fld_Hc_up1s_s               PASS     
  invariant _s    : fld_Hc_lo2s_s <= fld_Hc_hi2s_s               PASS     
  invariant _s    : fld_Hc_est2_s <= fld_Hc_up1s_s               PASS     
  invariant joint : bonf_loH <= fld_H_est2                       PASS     
  invariant jointS: fld_Hc_est2_s <= bonf_upHc_s                 PASS     
  invariant IJ    : mr_H_lo <= est <= mr_H_hi                    PASS     
  gamma (joint)   in [0.025, 0.05]                               PASS     [0.02500, 0.02800]
  gamma (joint-s) in [0.025, 0.05]                               PASS     [0.02500, 0.02800]
  identity: field-s inverted around the same beta-tilde^c        PASS     max |diff| = 1.42e-14
  identity: Bonferroni harm bound joint == joint_s where draw counts agree PASS     (1998 of 1998 agree; max |diff| 0)
  identity: lo1s = beta-tilde - q95; up1s = beta-tilde^c - q05   PASS     max |diff| = 0
  p-hat in [0, 1]                                                PASS     (mean 0.079, share < 0.5: 0.999)
  >> CLASSIFICATION           : sens 0.2679 ppv 0.4122 | mean |Hhat| (n_harm) 111.6 | mdf1 72.1
  --- same-draws: mdsgnb20 -> mdf1 ---
  [mdsgnb20 -> mdf1] same sim_id set (2000 rows)                 PASS     (comparator 2000 rows)
  [mdsgnb20 -> mdf1] n_true identical() on all rows              PASS     
  [mdsgnb20 -> mdf1] oracle columns (H and Hc) <= 1e-08 relative on all rows PASS     (max 2.78e-11)
  --- same-draws: mdf1 -> mdsgnb20 ---
  [mdf1 -> mdsgnb20] same sim_id set (2000 rows)                 PASS     (comparator 2000 rows)
  [mdf1 -> mdsgnb20] n_true identical() on all rows              PASS     
  [mdf1 -> mdsgnb20] oracle columns (H and Hc) <= 1e-08 relative on all rows PASS     (max 2.78e-11)
  size <= 100 MB: fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdsgnb20_res_1_1000.rds PASS     (744612 B)
  size <= 100 MB: fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdsgnb20_res_1001_2000.rds PASS     (743883 B)
  size <= 100 MB: fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdsgnb20_combined_1_2000.rds PASS     (1467414 B)
  timing: fit_mr_secs mean 56.2 median 56.8 p90 65.8 max 77.0 | fld_H_secs mean 32.7 | fld_Hc_secs mean 3.33

GATE_COUNTS run=58 passed=58 failed=0
```

## md120_n500 — md 120, n 500

- Stem: `fs_effMaxSG_mr_field_md120_knoise0_n500_nb20_mdsgnb20`; HEAD before this cell's commit: 76d6a8fd; workers 63; threads 1.
- Knobs: `FS_MD_FOCUS=effMaxSG FS_MD_NBHD=0.20 FS_MD_FIELD_SCALEC=selected FS_MD_CI=field FS_MD_CAMPAIGN=mdsgnb20 FS_MD_FB=none FS_MD_MD=120 FS_MD_N=500 FS_MD_WORKERS=63`.
- Seeds: 8316951 + sim_id; batches sim_id 1-1000 and 1001-2000, then combine.
- Cell wall: 2521 s; cumulative render wall 2521 s (ceiling 19071 s).
- Render: WALL_SECONDS=1245 RC=0 PEAK_MB=72511 OUT=fs_effMaxSG_mr_field_md120_knoise0_n500_nb20_mdsgnb20_batch_1001_2000.html
- Render: WALL_SECONDS=1235 RC=0 PEAK_MB=72875 OUT=fs_effMaxSG_mr_field_md120_knoise0_n500_nb20_mdsgnb20_batch_1_1000.html
- Render: WALL_SECONDS=41 RC=0 PEAK_MB=1168 OUT=fs_effMaxSG_mr_field_md120_knoise0_n500_nb20_mdsgnb20_combine_1_2000.html
- Gate counts: `GATE_COUNTS run=58 passed=58 failed=0`.

```

########## GATE 2 (mdsgnb20): md120_n500 -- md 120, n 500 ##########
  bundle: ../mr_md_harm/fs_effMaxSG_mr_field_md120_knoise0_n500_nb20_mdsgnb20_d5000/fs_effMaxSG_mr_field_md120_knoise0_n500_nb20_mdsgnb20_combined_1_2000.rds
  combined payload on disk                                       PASS     
  --- combine assertions ---
  exactly 2 batch files, res_1_1000 and res_1001_2000            PASS     (2)
  2,000 rows                                                     PASS     (2000)
  combined sim_id == 1:2000                                      PASS     
  batch sim_id sets 1:1000 and 1001:2000                         PASS     (batch1 1000, batch2 1000)
  batch files match the combined bundle on every column          PASS     (156 columns)
  no CONFIG-ERROR replicate                                      PASS     (0)
  --- meta (both batches) ---
  meta: subgroup_method == consistency                           PASS     (consistency / consistency)
  meta: sg_focus == effMaxSG                                     PASS     (effMaxSG / effMaxSG)
  meta: effect_neighborhood == 0.2                               PASS     (0.2 / 0.2)
  meta: selection_rule == neighborhood                           PASS     (neighborhood / neighborhood)
  meta: consistency_method == resample                           PASS     (resample / resample)
  meta: ci_method == field                                       PASS     (field / field)
  meta: mr_draws == 5000                                         PASS     (5000 / 5000)
  meta: field_uniform == FALSE                                   PASS     (FALSE / FALSE)
  meta: field_complement == TRUE                                 PASS     (TRUE / TRUE)
  meta: field_scale_complement == selected                       PASS     (selected / selected)
  meta: ij_residual == two_term                                  PASS     (two_term / two_term)
  meta: return_reselection == TRUE                               PASS     (TRUE / TRUE)
  meta: fb_mode == none                                          PASS     (none / none)
  meta: seed_base == 8316951                                     PASS     (8316951 / 8316951)
  meta: campaign_tag == mdsgnb20                                 PASS     (mdsgnb20 / mdsgnb20)
  meta: n_sample == 500                                          PASS     (500 / 500)
  meta: null_cell == FALSE                                       PASS     (FALSE / FALSE)
  meta: effect_threshold == 30                                   PASS     (30 / 30)
  meta: consistency_threshold == 10                              PASS     (10 / 10)
  meta: pkg_version == 0.3.5                                     PASS     (0.3.5 / 0.3.5)
  meta: hostname == pop-os                                       PASS     (pop-os / pop-os)
  meta: n_workers == 63                                          PASS     (63 / 63)
  meta seed_base 8316951 | host pop-os | R 4.6.1 | built_at 2026-09-16 00:46:24 / 2026-09-16 00:25:40
  >> DECLARATION RATE         : 1.0000 (2000 / 2000)
  MR failures on declared replicates <= max(20, 2 x mdf1's 0)    PASS     (0)
  every finiteness column present                                PASS     
  harm products finite on declared replicates with a field block PASS     (2000 rows)
  nine fld_Hc_*_s and nine fld_joint_s_* columns present         PASS     
  fld_Hc_*_s finite on all 2000 filled replicates                PASS     
  fld_joint_s_* finite on all 2000 filled replicates             PASS     
  complement field filled on 2000 of 2000 declared replicates (0 notes)
  invariant harm  : fld_H_lo1s <= fld_H_est2                     PASS     
  invariant compl : fld_Hc_est2 <= fld_Hc_up1s                   PASS     
  invariant _s    : fld_Hc_lo1s_s <= fld_Hc_up1s_s               PASS     
  invariant _s    : fld_Hc_lo2s_s <= fld_Hc_hi2s_s               PASS     
  invariant _s    : fld_Hc_est2_s <= fld_Hc_up1s_s               PASS     
  invariant joint : bonf_loH <= fld_H_est2                       PASS     
  invariant jointS: fld_Hc_est2_s <= bonf_upHc_s                 PASS     
  invariant IJ    : mr_H_lo <= est <= mr_H_hi                    PASS     
  gamma (joint)   in [0.025, 0.05]                               PASS     [0.02500, 0.02700]
  gamma (joint-s) in [0.025, 0.05]                               PASS     [0.02500, 0.02700]
  identity: field-s inverted around the same beta-tilde^c        PASS     max |diff| = 1.42e-14
  identity: Bonferroni harm bound joint == joint_s where draw counts agree PASS     (2000 of 2000 agree; max |diff| 0)
  identity: lo1s = beta-tilde - q95; up1s = beta-tilde^c - q05   PASS     max |diff| = 0
  p-hat in [0, 1]                                                PASS     (mean 0.089, share < 0.5: 1.000)
  >> CLASSIFICATION           : sens 0.7116 ppv 0.7980 | mean |Hhat| (n_harm) 152.7 | mdf1 74.8
  --- same-draws: mdsgnb20 -> mdf1 ---
  [mdsgnb20 -> mdf1] same sim_id set (2000 rows)                 PASS     (comparator 2000 rows)
  [mdsgnb20 -> mdf1] n_true identical() on all rows              PASS     
  [mdsgnb20 -> mdf1] oracle columns (H and Hc) <= 1e-08 relative on all rows PASS     (max 2.78e-11)
  --- same-draws: mdf1 -> mdsgnb20 ---
  [mdf1 -> mdsgnb20] same sim_id set (2000 rows)                 PASS     (comparator 2000 rows)
  [mdf1 -> mdsgnb20] n_true identical() on all rows              PASS     
  [mdf1 -> mdsgnb20] oracle columns (H and Hc) <= 1e-08 relative on all rows PASS     (max 2.78e-11)
  size <= 100 MB: fs_effMaxSG_mr_field_md120_knoise0_n500_nb20_mdsgnb20_res_1_1000.rds PASS     (731187 B)
  size <= 100 MB: fs_effMaxSG_mr_field_md120_knoise0_n500_nb20_mdsgnb20_res_1001_2000.rds PASS     (731111 B)
  size <= 100 MB: fs_effMaxSG_mr_field_md120_knoise0_n500_nb20_mdsgnb20_combined_1_2000.rds PASS     (1440770 B)
  timing: fit_mr_secs mean 68.0 median 69.1 p90 77.6 max 87.5 | fld_H_secs mean 41.3 | fld_Hc_secs mean 3.51

GATE_COUNTS run=58 passed=58 failed=0
```

## null_n500 — md null, n 500

- Stem: `fs_effMaxSG_mr_field_mdnull_knoise0_n500_nb20_mdsgnb20`; HEAD before this cell's commit: b12f5983; workers 63; threads 1.
- Knobs: `FS_MD_FOCUS=effMaxSG FS_MD_NBHD=0.20 FS_MD_FIELD_SCALEC=selected FS_MD_CI=field FS_MD_CAMPAIGN=mdsgnb20 FS_MD_FB=none FS_MD_MD=null FS_MD_N=500 FS_MD_WORKERS=63`.
- Seeds: 8316951 + sim_id; batches sim_id 1-1000 and 1001-2000, then combine.
- Cell wall: 2065 s; cumulative render wall 4586 s (ceiling 19071 s).
- Render: WALL_SECONDS=1010 RC=0 PEAK_MB=73020 OUT=fs_effMaxSG_mr_field_mdnull_knoise0_n500_nb20_mdsgnb20_batch_1001_2000.html
- Render: WALL_SECONDS=1014 RC=0 PEAK_MB=72582 OUT=fs_effMaxSG_mr_field_mdnull_knoise0_n500_nb20_mdsgnb20_batch_1_1000.html
- Render: WALL_SECONDS=41 RC=0 PEAK_MB=1184 OUT=fs_effMaxSG_mr_field_mdnull_knoise0_n500_nb20_mdsgnb20_combine_1_2000.html
- Gate counts: `GATE_COUNTS run=58 passed=58 failed=0`.

```

########## GATE 2 (mdsgnb20): null_n500 -- md null, n 500 ##########
  bundle: ../mr_md_harm/fs_effMaxSG_mr_field_mdnull_knoise0_n500_nb20_mdsgnb20_d5000/fs_effMaxSG_mr_field_mdnull_knoise0_n500_nb20_mdsgnb20_combined_1_2000.rds
  combined payload on disk                                       PASS     
  --- combine assertions ---
  exactly 2 batch files, res_1_1000 and res_1001_2000            PASS     (2)
  2,000 rows                                                     PASS     (2000)
  combined sim_id == 1:2000                                      PASS     
  batch sim_id sets 1:1000 and 1001:2000                         PASS     (batch1 1000, batch2 1000)
  batch files match the combined bundle on every column          PASS     (156 columns)
  no CONFIG-ERROR replicate                                      PASS     (0)
  --- meta (both batches) ---
  meta: subgroup_method == consistency                           PASS     (consistency / consistency)
  meta: sg_focus == effMaxSG                                     PASS     (effMaxSG / effMaxSG)
  meta: effect_neighborhood == 0.2                               PASS     (0.2 / 0.2)
  meta: selection_rule == neighborhood                           PASS     (neighborhood / neighborhood)
  meta: consistency_method == resample                           PASS     (resample / resample)
  meta: ci_method == field                                       PASS     (field / field)
  meta: mr_draws == 5000                                         PASS     (5000 / 5000)
  meta: field_uniform == FALSE                                   PASS     (FALSE / FALSE)
  meta: field_complement == TRUE                                 PASS     (TRUE / TRUE)
  meta: field_scale_complement == selected                       PASS     (selected / selected)
  meta: ij_residual == two_term                                  PASS     (two_term / two_term)
  meta: return_reselection == TRUE                               PASS     (TRUE / TRUE)
  meta: fb_mode == none                                          PASS     (none / none)
  meta: seed_base == 8316951                                     PASS     (8316951 / 8316951)
  meta: campaign_tag == mdsgnb20                                 PASS     (mdsgnb20 / mdsgnb20)
  meta: n_sample == 500                                          PASS     (500 / 500)
  meta: null_cell == TRUE                                        PASS     (TRUE / TRUE)
  meta: effect_threshold == 30                                   PASS     (30 / 30)
  meta: consistency_threshold == 10                              PASS     (10 / 10)
  meta: pkg_version == 0.3.5                                     PASS     (0.3.5 / 0.3.5)
  meta: hostname == pop-os                                       PASS     (pop-os / pop-os)
  meta: n_workers == 63                                          PASS     (63 / 63)
  meta seed_base 8316951 | host pop-os | R 4.6.1 | built_at 2026-09-16 01:20:50 / 2026-09-16 01:03:57
  >> DECLARATION RATE         : 0.9965 (1993 / 2000)
  MR failures on declared replicates <= max(20, 2 x mdf1's 0)    PASS     (0)
  every finiteness column present                                PASS     
  harm products finite on declared replicates with a field block PASS     (1993 rows)
  nine fld_Hc_*_s and nine fld_joint_s_* columns present         PASS     
  fld_Hc_*_s finite on all 1993 filled replicates                PASS     
  fld_joint_s_* finite on all 1993 filled replicates             PASS     
  complement field filled on 1993 of 1993 declared replicates (0 notes)
  invariant harm  : fld_H_lo1s <= fld_H_est2                     PASS     
  invariant compl : fld_Hc_est2 <= fld_Hc_up1s                   PASS     
  invariant _s    : fld_Hc_lo1s_s <= fld_Hc_up1s_s               PASS     
  invariant _s    : fld_Hc_lo2s_s <= fld_Hc_hi2s_s               PASS     
  invariant _s    : fld_Hc_est2_s <= fld_Hc_up1s_s               PASS     
  invariant joint : bonf_loH <= fld_H_est2                       PASS     
  invariant jointS: fld_Hc_est2_s <= bonf_upHc_s                 PASS     
  invariant IJ    : mr_H_lo <= est <= mr_H_hi                    PASS     
  gamma (joint)   in [0.025, 0.05]                               PASS     [0.02500, 0.02800]
  gamma (joint-s) in [0.025, 0.05]                               PASS     [0.02500, 0.02800]
  identity: field-s inverted around the same beta-tilde^c        PASS     max |diff| = 7.11e-15
  identity: Bonferroni harm bound joint == joint_s where draw counts agree PASS     (1993 of 1993 agree; max |diff| 0)
  identity: lo1s = beta-tilde - q95; up1s = beta-tilde^c - q05   PASS     max |diff| = 0
  p-hat in [0, 1]                                                PASS     (mean 0.082, share < 0.5: 0.999)
  >> CLASSIFICATION           : sens NaN ppv 0.0000 | mean |Hhat| (n_harm) 107.5 | mdf1 72.1
  --- same-draws: mdsgnb20 -> mdf1 ---
  [mdsgnb20 -> mdf1] same sim_id set (2000 rows)                 PASS     (comparator 2000 rows)
  [mdsgnb20 -> mdf1] n_true identical() on all rows              PASS     
  [mdsgnb20 -> mdf1] oracle columns (complement only) <= 1e-08 relative on all rows PASS     (max 2.54e-10)
  --- same-draws: mdf1 -> mdsgnb20 ---
  [mdf1 -> mdsgnb20] same sim_id set (2000 rows)                 PASS     (comparator 2000 rows)
  [mdf1 -> mdsgnb20] n_true identical() on all rows              PASS     
  [mdf1 -> mdsgnb20] oracle columns (complement only) <= 1e-08 relative on all rows PASS     (max 2.54e-10)
  size <= 100 MB: fs_effMaxSG_mr_field_mdnull_knoise0_n500_nb20_mdsgnb20_res_1_1000.rds PASS     (684263 B)
  size <= 100 MB: fs_effMaxSG_mr_field_mdnull_knoise0_n500_nb20_mdsgnb20_res_1001_2000.rds PASS     (682752 B)
  size <= 100 MB: fs_effMaxSG_mr_field_mdnull_knoise0_n500_nb20_mdsgnb20_combined_1_2000.rds PASS     (1349124 B)
  timing: fit_mr_secs mean 53.9 median 54.4 p90 63.5 max 75.3 | fld_H_secs mean 31.3 | fld_Hc_secs mean 3.28

GATE_COUNTS run=58 passed=58 failed=0
```

## md40_n700 — md 40, n 700

- Stem: `fs_effMaxSG_mr_field_md40_knoise0_n700_nb20_mdsgnb20`; HEAD before this cell's commit: dbeee517; workers 63; threads 1.
- Knobs: `FS_MD_FOCUS=effMaxSG FS_MD_NBHD=0.20 FS_MD_FIELD_SCALEC=selected FS_MD_CI=field FS_MD_CAMPAIGN=mdsgnb20 FS_MD_FB=none FS_MD_MD=40 FS_MD_N=700 FS_MD_WORKERS=63`.
- Seeds: 8316951 + sim_id; batches sim_id 1-1000 and 1001-2000, then combine.
- Cell wall: 3076 s; cumulative render wall 7662 s (ceiling 19071 s).
- Render: WALL_SECONDS=1503 RC=0 PEAK_MB=82506 OUT=fs_effMaxSG_mr_field_md40_knoise0_n700_nb20_mdsgnb20_batch_1001_2000.html
- Render: WALL_SECONDS=1533 RC=0 PEAK_MB=81838 OUT=fs_effMaxSG_mr_field_md40_knoise0_n700_nb20_mdsgnb20_batch_1_1000.html
- Render: WALL_SECONDS=40 RC=0 PEAK_MB=1199 OUT=fs_effMaxSG_mr_field_md40_knoise0_n700_nb20_mdsgnb20_combine_1_2000.html
- Gate counts: `GATE_COUNTS run=58 passed=58 failed=0`.

```

########## GATE 2 (mdsgnb20): md40_n700 -- md 40, n 700 ##########
  bundle: ../mr_md_harm/fs_effMaxSG_mr_field_md40_knoise0_n700_nb20_mdsgnb20_d5000/fs_effMaxSG_mr_field_md40_knoise0_n700_nb20_mdsgnb20_combined_1_2000.rds
  combined payload on disk                                       PASS     
  --- combine assertions ---
  exactly 2 batch files, res_1_1000 and res_1001_2000            PASS     (2)
  2,000 rows                                                     PASS     (2000)
  combined sim_id == 1:2000                                      PASS     
  batch sim_id sets 1:1000 and 1001:2000                         PASS     (batch1 1000, batch2 1000)
  batch files match the combined bundle on every column          PASS     (156 columns)
  no CONFIG-ERROR replicate                                      PASS     (0)
  --- meta (both batches) ---
  meta: subgroup_method == consistency                           PASS     (consistency / consistency)
  meta: sg_focus == effMaxSG                                     PASS     (effMaxSG / effMaxSG)
  meta: effect_neighborhood == 0.2                               PASS     (0.2 / 0.2)
  meta: selection_rule == neighborhood                           PASS     (neighborhood / neighborhood)
  meta: consistency_method == resample                           PASS     (resample / resample)
  meta: ci_method == field                                       PASS     (field / field)
  meta: mr_draws == 5000                                         PASS     (5000 / 5000)
  meta: field_uniform == FALSE                                   PASS     (FALSE / FALSE)
  meta: field_complement == TRUE                                 PASS     (TRUE / TRUE)
  meta: field_scale_complement == selected                       PASS     (selected / selected)
  meta: ij_residual == two_term                                  PASS     (two_term / two_term)
  meta: return_reselection == TRUE                               PASS     (TRUE / TRUE)
  meta: fb_mode == none                                          PASS     (none / none)
  meta: seed_base == 8316951                                     PASS     (8316951 / 8316951)
  meta: campaign_tag == mdsgnb20                                 PASS     (mdsgnb20 / mdsgnb20)
  meta: n_sample == 700                                          PASS     (700 / 700)
  meta: null_cell == FALSE                                       PASS     (FALSE / FALSE)
  meta: effect_threshold == 30                                   PASS     (30 / 30)
  meta: consistency_threshold == 10                              PASS     (10 / 10)
  meta: pkg_version == 0.3.5                                     PASS     (0.3.5 / 0.3.5)
  meta: hostname == pop-os                                       PASS     (pop-os / pop-os)
  meta: n_workers == 63                                          PASS     (63 / 63)
  meta seed_base 8316951 | host pop-os | R 4.6.1 | built_at 2026-09-16 02:12:06 / 2026-09-16 01:47:01
  >> DECLARATION RATE         : 0.9995 (1999 / 2000)
  MR failures on declared replicates <= max(20, 2 x mdf1's 0)    PASS     (0)
  every finiteness column present                                PASS     
  harm products finite on declared replicates with a field block PASS     (1999 rows)
  nine fld_Hc_*_s and nine fld_joint_s_* columns present         PASS     
  fld_Hc_*_s finite on all 1999 filled replicates                PASS     
  fld_joint_s_* finite on all 1999 filled replicates             PASS     
  complement field filled on 1999 of 1999 declared replicates (0 notes)
  invariant harm  : fld_H_lo1s <= fld_H_est2                     PASS     
  invariant compl : fld_Hc_est2 <= fld_Hc_up1s                   PASS     
  invariant _s    : fld_Hc_lo1s_s <= fld_Hc_up1s_s               PASS     
  invariant _s    : fld_Hc_lo2s_s <= fld_Hc_hi2s_s               PASS     
  invariant _s    : fld_Hc_est2_s <= fld_Hc_up1s_s               PASS     
  invariant joint : bonf_loH <= fld_H_est2                       PASS     
  invariant jointS: fld_Hc_est2_s <= bonf_upHc_s                 PASS     
  invariant IJ    : mr_H_lo <= est <= mr_H_hi                    PASS     
  gamma (joint)   in [0.025, 0.05]                               PASS     [0.02500, 0.02700]
  gamma (joint-s) in [0.025, 0.05]                               PASS     [0.02500, 0.02700]
  identity: field-s inverted around the same beta-tilde^c        PASS     max |diff| = 1.42e-14
  identity: Bonferroni harm bound joint == joint_s where draw counts agree PASS     (1999 of 1999 agree; max |diff| 0)
  identity: lo1s = beta-tilde - q95; up1s = beta-tilde^c - q05   PASS     max |diff| = 0
  p-hat in [0, 1]                                                PASS     (mean 0.080, share < 0.5: 0.999)
  >> CLASSIFICATION           : sens 0.2041 ppv 0.4101 | mean |Hhat| (n_harm) 118.0 | mdf1 73.3
  --- same-draws: mdsgnb20 -> mdf1 ---
  [mdsgnb20 -> mdf1] same sim_id set (2000 rows)                 PASS     (comparator 2000 rows)
  [mdsgnb20 -> mdf1] n_true identical() on all rows              PASS     
  [mdsgnb20 -> mdf1] oracle columns (H and Hc) <= 1e-08 relative on all rows PASS     (max 4.01e-11)
  --- same-draws: mdf1 -> mdsgnb20 ---
  [mdf1 -> mdsgnb20] same sim_id set (2000 rows)                 PASS     (comparator 2000 rows)
  [mdf1 -> mdsgnb20] n_true identical() on all rows              PASS     
  [mdf1 -> mdsgnb20] oracle columns (H and Hc) <= 1e-08 relative on all rows PASS     (max 4.01e-11)
  size <= 100 MB: fs_effMaxSG_mr_field_md40_knoise0_n700_nb20_mdsgnb20_res_1_1000.rds PASS     (740620 B)
  size <= 100 MB: fs_effMaxSG_mr_field_md40_knoise0_n700_nb20_mdsgnb20_res_1001_2000.rds PASS     (740802 B)
  size <= 100 MB: fs_effMaxSG_mr_field_md40_knoise0_n700_nb20_mdsgnb20_combined_1_2000.rds PASS     (1461232 B)
  timing: fit_mr_secs mean 84.9 median 86.9 p90 99.8 max 113.7 | fld_H_secs mean 42.0 | fld_Hc_secs mean 5.75

GATE_COUNTS run=58 passed=58 failed=0
```

