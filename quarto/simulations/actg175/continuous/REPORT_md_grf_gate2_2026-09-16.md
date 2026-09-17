## md40_n500 — md 40, n 500

- Stem: `grf_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdgrf`; HEAD before this cell's commit: f5cc256d; workers 63; threads 1.
- Knobs: `FS_MD_METHOD=grf FS_MD_DMIN_GRF=30 FS_MD_FOCUS=effMaxSG FS_MD_NBHD=0.20 FS_MD_FIELD_SCALEC=selected FS_MD_CI=field FS_MD_CAMPAIGN=mdgrf FS_MD_FB=none FS_MD_MD=40 FS_MD_N=500 FS_MD_WORKERS=63`.
- Seeds: 8316951 + sim_id; batches sim_id 1-1000 and 1001-2000, then combine.
- Cell wall: 1748 s; cumulative render wall 1748 s (ceiling 11870 s).
- Render: WALL_SECONDS=856 RC=0 PEAK_MB=57837 OUT=grf_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdgrf_batch_1001_2000.html
- Render: WALL_SECONDS=851 RC=0 PEAK_MB=57326 OUT=grf_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdgrf_batch_1_1000.html
- Render: WALL_SECONDS=41 RC=0 PEAK_MB=1178 OUT=grf_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdgrf_combine_1_2000.html
- Gate counts: `GATE_COUNTS run=66 passed=66 failed=0`.

```

########## GATE 2 (mdgrf): md40_n500 -- md 40, n 500 ##########
  Every coverage figure of this campaign is conditional on the proposed family.
  bundle: ../mr_md_harm/grf_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdgrf_d5000/grf_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdgrf_combined_1_2000.rds
  combined payload on disk                                       PASS     
  --- combine assertions ---
  exactly 2 batch files, res_1_1000 and res_1001_2000            PASS     (2)
  2,000 rows                                                     PASS     (2000)
  combined sim_id == 1:2000                                      PASS     
  batch sim_id sets 1:1000 and 1001:2000                         PASS     (batch1 1000, batch2 1000)
  batch files match the combined bundle on every column          PASS     (172 columns)
  no CONFIG-ERROR replicate                                      PASS     (0)
  --- meta (both batches) ---
  meta: subgroup_method == grf                                   PASS     (grf / grf)
  meta: dmin_grf == 30                                           PASS     (30 / 30)
  meta: grf_selection == frontier                                PASS     (frontier / frontier)
  meta: grf_depth == 2                                           PASS     (2 / 2)
  meta: grf_select_statistic == effect                           PASS     (effect / effect)
  meta: sg_focus == effMaxSG                                     PASS     (effMaxSG / effMaxSG)
  meta: effect_neighborhood == 0.2                               PASS     (0.2 / 0.2)
  meta: selection_rule == neighborhood                           PASS     (neighborhood / neighborhood)
  meta: ci_method == field                                       PASS     (field / field)
  meta: mr_draws == 5000                                         PASS     (5000 / 5000)
  meta: field_uniform == FALSE                                   PASS     (FALSE / FALSE)
  meta: field_complement == TRUE                                 PASS     (TRUE / TRUE)
  meta: field_scale_complement == selected                       PASS     (selected / selected)
  meta: ij_residual == two_term                                  PASS     (two_term / two_term)
  meta: return_reselection == TRUE                               PASS     (TRUE / TRUE)
  meta: fb_mode == none                                          PASS     (none / none)
  meta: seed_base == 8316951                                     PASS     (8316951 / 8316951)
  meta: campaign_tag == mdgrf                                    PASS     (mdgrf / mdgrf)
  meta: n_sample == 500                                          PASS     (500 / 500)
  meta: null_cell == FALSE                                       PASS     (FALSE / FALSE)
  meta: effect_threshold == 30                                   PASS     (30 / 30)
  meta: consistency_threshold == 10                              PASS     (10 / 10)
  meta: pkg_version == 0.3.5                                     PASS     (0.3.5 / 0.3.5)
  meta: hostname == pop-os                                       PASS     (pop-os / pop-os)
  meta: n_workers == 63                                          PASS     (63 / 63)
  meta seed_base 8316951 | host pop-os | R 4.6.1 | built_at 2026-09-16 22:35:52 / 2026-09-16 22:21:37
  >> DECLARATION RATE         : 1.0000 (2000 / 2000)
  MR failures on declared replicates <= 40                       PASS     (0)
  >> ADMITTED_N (qualified)   : min 15  q25 383  MED 633  q75 865.25  p90 1020  max 1277  (mean 624.6, n 2000); == 0 on 0 row(s)
  >> N_FAMILY (ENUMERATED POOL, not the qualified set): min 1170  med 1199  p90 1283  max 1526
  admitted_n recorded on every row                               PASS     (2000 of 2000)
  admitted_n >= 1 on every declared replicate                    PASS     
  n_family finite on every declared replicate with a gate        PASS     
  p_hat_H and p_hat_sum recorded on declared replicates with a gate PASS     
  p-hat validity (0<=p<=1, p_H<=sum)                             PASS     
  STRUCTURAL n_cons_qual  present, all-NA -- STRUCTURAL on GRF (no consistency screen), NOT a failure
  STRUCTURAL band_n       present, all-NA -- STRUCTURAL on GRF (no consistency screen), NOT a failure
  >> WARNINGS                 : 0 of 2000 rows carry warn_msg; distinct messages: none
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
  >> p-hat(H): mean 0.099, share < 0.5: 0.998 | CLASSIFICATION: sens 0.2497 ppv 0.3909 | mean |Hhat| (n_harm) 109.5
  mdsgnb20 comparator on disk                                    PASS     fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdsgnb20_combined_1_2000.rds
  --- same-draws: mdgrf -> mdsgnb20 ---
  [mdgrf -> mdsgnb20] same sim_id set (2000 rows)                PASS     (comparator 2000 rows)
  [mdgrf -> mdsgnb20] n_true identical() on all rows             PASS     
  [mdgrf -> mdsgnb20] oracle columns (H and Hc) <= 1e-08 relative on all rows PASS     (max 0)
  --- same-draws: mdsgnb20 -> mdgrf ---
  [mdsgnb20 -> mdgrf] same sim_id set (2000 rows)                PASS     (comparator 2000 rows)
  [mdsgnb20 -> mdgrf] n_true identical() on all rows             PASS     
  [mdsgnb20 -> mdgrf] oracle columns (H and Hc) <= 1e-08 relative on all rows PASS     (max 0)
  FS comparator (same cell): declaration 0.9990 | mean |Hhat| 111.6
  size <= 100 MB: grf_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdgrf_res_1_1000.rds PASS     (766849 B)
  size <= 100 MB: grf_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdgrf_res_1001_2000.rds PASS     (767454 B)
  size <= 100 MB: grf_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdgrf_combined_1_2000.rds PASS     (1508903 B)
  timing: fit_mr_secs mean 45.1 median 46.0 p90 49.3 max 53.9 | id_secs (GRF fit) mean 6.00 median 6.14 max 8.16 | fld_H_secs mean 27.2 | fld_Hc_secs mean 1.72

GATE_COUNTS run=66 passed=66 failed=0
```

## md120_n500 — md 120, n 500

- Stem: `grf_effMaxSG_mr_field_md120_knoise0_n500_nb20_mdgrf`; HEAD before this cell's commit: c4c49572; workers 63; threads 1.
- Knobs: `FS_MD_METHOD=grf FS_MD_DMIN_GRF=30 FS_MD_FOCUS=effMaxSG FS_MD_NBHD=0.20 FS_MD_FIELD_SCALEC=selected FS_MD_CI=field FS_MD_CAMPAIGN=mdgrf FS_MD_FB=none FS_MD_MD=120 FS_MD_N=500 FS_MD_WORKERS=63`.
- Seeds: 8316951 + sim_id; batches sim_id 1-1000 and 1001-2000, then combine.
- Cell wall: 1934 s; cumulative render wall 3682 s (ceiling 11870 s).
- Render: WALL_SECONDS=954 RC=0 PEAK_MB=57621 OUT=grf_effMaxSG_mr_field_md120_knoise0_n500_nb20_mdgrf_batch_1001_2000.html
- Render: WALL_SECONDS=944 RC=0 PEAK_MB=57554 OUT=grf_effMaxSG_mr_field_md120_knoise0_n500_nb20_mdgrf_batch_1_1000.html
- Render: WALL_SECONDS=36 RC=0 PEAK_MB=1155 OUT=grf_effMaxSG_mr_field_md120_knoise0_n500_nb20_mdgrf_combine_1_2000.html
- Gate counts: `GATE_COUNTS run=66 passed=66 failed=0`.

```

########## GATE 2 (mdgrf): md120_n500 -- md 120, n 500 ##########
  Every coverage figure of this campaign is conditional on the proposed family.
  bundle: ../mr_md_harm/grf_effMaxSG_mr_field_md120_knoise0_n500_nb20_mdgrf_d5000/grf_effMaxSG_mr_field_md120_knoise0_n500_nb20_mdgrf_combined_1_2000.rds
  combined payload on disk                                       PASS     
  --- combine assertions ---
  exactly 2 batch files, res_1_1000 and res_1001_2000            PASS     (2)
  2,000 rows                                                     PASS     (2000)
  combined sim_id == 1:2000                                      PASS     
  batch sim_id sets 1:1000 and 1001:2000                         PASS     (batch1 1000, batch2 1000)
  batch files match the combined bundle on every column          PASS     (172 columns)
  no CONFIG-ERROR replicate                                      PASS     (0)
  --- meta (both batches) ---
  meta: subgroup_method == grf                                   PASS     (grf / grf)
  meta: dmin_grf == 30                                           PASS     (30 / 30)
  meta: grf_selection == frontier                                PASS     (frontier / frontier)
  meta: grf_depth == 2                                           PASS     (2 / 2)
  meta: grf_select_statistic == effect                           PASS     (effect / effect)
  meta: sg_focus == effMaxSG                                     PASS     (effMaxSG / effMaxSG)
  meta: effect_neighborhood == 0.2                               PASS     (0.2 / 0.2)
  meta: selection_rule == neighborhood                           PASS     (neighborhood / neighborhood)
  meta: ci_method == field                                       PASS     (field / field)
  meta: mr_draws == 5000                                         PASS     (5000 / 5000)
  meta: field_uniform == FALSE                                   PASS     (FALSE / FALSE)
  meta: field_complement == TRUE                                 PASS     (TRUE / TRUE)
  meta: field_scale_complement == selected                       PASS     (selected / selected)
  meta: ij_residual == two_term                                  PASS     (two_term / two_term)
  meta: return_reselection == TRUE                               PASS     (TRUE / TRUE)
  meta: fb_mode == none                                          PASS     (none / none)
  meta: seed_base == 8316951                                     PASS     (8316951 / 8316951)
  meta: campaign_tag == mdgrf                                    PASS     (mdgrf / mdgrf)
  meta: n_sample == 500                                          PASS     (500 / 500)
  meta: null_cell == FALSE                                       PASS     (FALSE / FALSE)
  meta: effect_threshold == 30                                   PASS     (30 / 30)
  meta: consistency_threshold == 10                              PASS     (10 / 10)
  meta: pkg_version == 0.3.5                                     PASS     (0.3.5 / 0.3.5)
  meta: hostname == pop-os                                       PASS     (pop-os / pop-os)
  meta: n_workers == 63                                          PASS     (63 / 63)
  meta seed_base 8316951 | host pop-os | R 4.6.1 | built_at 2026-09-16 23:08:15 / 2026-09-16 22:52:23
  >> DECLARATION RATE         : 1.0000 (2000 / 2000)
  MR failures on declared replicates <= 40                       PASS     (0)
  >> ADMITTED_N (qualified)   : min 327  q25 1003  MED 1085  q75 1146  p90 1186.1  max 1420  (mean 1059.8, n 2000); == 0 on 0 row(s)
  >> N_FAMILY (ENUMERATED POOL, not the qualified set): min 1170  med 1199  p90 1283  max 1526
  admitted_n recorded on every row                               PASS     (2000 of 2000)
  admitted_n >= 1 on every declared replicate                    PASS     
  n_family finite on every declared replicate with a gate        PASS     
  p_hat_H and p_hat_sum recorded on declared replicates with a gate PASS     
  p-hat validity (0<=p<=1, p_H<=sum)                             PASS     
  STRUCTURAL n_cons_qual  present, all-NA -- STRUCTURAL on GRF (no consistency screen), NOT a failure
  STRUCTURAL band_n       present, all-NA -- STRUCTURAL on GRF (no consistency screen), NOT a failure
  >> WARNINGS                 : 0 of 2000 rows carry warn_msg; distinct messages: none
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
  >> p-hat(H): mean 0.128, share < 0.5: 0.996 | CLASSIFICATION: sens 0.5989 ppv 0.7149 | mean |Hhat| (n_harm) 143.6
  mdsgnb20 comparator on disk                                    PASS     fs_effMaxSG_mr_field_md120_knoise0_n500_nb20_mdsgnb20_combined_1_2000.rds
  --- same-draws: mdgrf -> mdsgnb20 ---
  [mdgrf -> mdsgnb20] same sim_id set (2000 rows)                PASS     (comparator 2000 rows)
  [mdgrf -> mdsgnb20] n_true identical() on all rows             PASS     
  [mdgrf -> mdsgnb20] oracle columns (H and Hc) <= 1e-08 relative on all rows PASS     (max 0)
  --- same-draws: mdsgnb20 -> mdgrf ---
  [mdsgnb20 -> mdgrf] same sim_id set (2000 rows)                PASS     (comparator 2000 rows)
  [mdsgnb20 -> mdgrf] n_true identical() on all rows             PASS     
  [mdsgnb20 -> mdgrf] oracle columns (H and Hc) <= 1e-08 relative on all rows PASS     (max 0)
  FS comparator (same cell): declaration 1.0000 | mean |Hhat| 152.7
  size <= 100 MB: grf_effMaxSG_mr_field_md120_knoise0_n500_nb20_mdgrf_res_1_1000.rds PASS     (756462 B)
  size <= 100 MB: grf_effMaxSG_mr_field_md120_knoise0_n500_nb20_mdgrf_res_1001_2000.rds PASS     (755524 B)
  size <= 100 MB: grf_effMaxSG_mr_field_md120_knoise0_n500_nb20_mdgrf_combined_1_2000.rds PASS     (1487382 B)
  timing: fit_mr_secs mean 51.2 median 52.3 p90 55.3 max 61.1 | id_secs (GRF fit) mean 6.50 median 6.62 max 9.22 | fld_H_secs mean 32.3 | fld_Hc_secs mean 1.86

GATE_COUNTS run=66 passed=66 failed=0
```

