## md40_n500 — md 40, n 500

- Stem: `dina_effMaxSG_mr_field_md40_knoise0_n500_nb20_mddina`; HEAD before this cell's commit: 37395393; workers 63; threads 1.
- Knobs: `FS_MD_METHOD=dina FS_MD_FOCUS=effMaxSG FS_MD_NBHD=0.20 FS_MD_FIELD_SCALEC=selected FS_MD_CI=field FS_MD_CAMPAIGN=mddina FS_MD_FB=none FS_MD_MD=40 FS_MD_N=500 FS_MD_WORKERS=63`.
- Seeds: 8316951 + sim_id; batches sim_id 1-1000 and 1001-2000, then combine.
- Cell wall: 5199 s; cumulative render wall 5199 s (ceiling 32939 s).
- Render: WALL_SECONDS=2569 RC=0 PEAK_MB=110447 OUT=dina_effMaxSG_mr_field_md40_knoise0_n500_nb20_mddina_batch_1001_2000.html
- Render: WALL_SECONDS=2589 RC=0 PEAK_MB=109791 OUT=dina_effMaxSG_mr_field_md40_knoise0_n500_nb20_mddina_batch_1_1000.html
- Render: WALL_SECONDS=41 RC=0 PEAK_MB=1185 OUT=dina_effMaxSG_mr_field_md40_knoise0_n500_nb20_mddina_combine_1_2000.html
- Gate counts: `GATE_COUNTS run=65 passed=65 failed=0`.

```

########## GATE 2 (mddina): md40_n500 -- md 40, n 500 ##########
  Every coverage figure of this campaign is conditional on the proposed family.
  bundle: ../mr_md_harm/dina_effMaxSG_mr_field_md40_knoise0_n500_nb20_mddina_d5000/dina_effMaxSG_mr_field_md40_knoise0_n500_nb20_mddina_combined_1_2000.rds
  combined payload on disk                                       PASS     
  --- combine assertions ---
  exactly 2 batch files, res_1_1000 and res_1001_2000            PASS     (2)
  2,000 rows                                                     PASS     (2000)
  combined sim_id == 1:2000                                      PASS     
  batch sim_id sets 1:1000 and 1001:2000                         PASS     (batch1 1000, batch2 1000)
  batch files match the combined bundle on every column          PASS     (175 columns)
  no CONFIG-ERROR replicate                                      PASS     (0)
  --- meta (both batches) ---
  meta: subgroup_method == dina                                  PASS     (dina / dina)
  meta: dina_select_statistic == effect                          PASS     (effect / effect)
  meta: dina_args == list()                                      PASS     (list() / list())
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
  meta: campaign_tag == mddina                                   PASS     (mddina / mddina)
  meta: n_sample == 500                                          PASS     (500 / 500)
  meta: null_cell == FALSE                                       PASS     (FALSE / FALSE)
  meta: effect_threshold == 30                                   PASS     (30 / 30)
  meta: consistency_threshold == 10                              PASS     (10 / 10)
  meta: pkg_version == 0.3.5                                     PASS     (0.3.5 / 0.3.5)
  meta: hostname == pop-os                                       PASS     (pop-os / pop-os)
  meta: n_workers == 63                                          PASS     (63 / 63)
  meta seed_base 8316951 | host pop-os | R 4.6.1 | built_at 2026-09-17 03:18:44 / 2026-09-17 02:35:58
  >> DECLARATION RATE         : 0.9980 (1996 / 2000)
  MR failures on declared replicates <= 40                       PASS     (0)
  >> DINA_PROPOSED_N (>= floor): min 2  q25 1284  MED 3507.5  q75 5386.5  p90 6088  max 6683  (mean 3341.7, n 1996); searched 8234-8916
  >> ADMITTED_N (qualified)   : min 1  q25 847.75  MED 2227.5  q75 4397.5  p90 5357.5  max 6494  (mean 2631.0, n 1996); == 0 on 0 row(s)
  >> N_FAMILY (MR's family: proposed candidates with >= n.min members): min 2  med 3507.5  p90 6088  max 6683
  E2 fields filled on every declared replicate                   PASS     (dina_searched_n 1996/1996, dina_proposed_n 1996/1996, dina_tau_min 1996/1996, admitted_n 1996/1996)
  STRUCTURAL E2 fields on non-detections: all NA TRUE (DINA's selection object is absent when nothing is selected)
  every proposed candidate at oriented tau-hat >= 30 (dina_tau_min, declared replicates) PASS     (min 30.0000)
  1 <= admitted_n <= dina_proposed_n on every declared replicate PASS     
  n_family finite on every declared replicate with a gate        PASS     
  p_hat_H and p_hat_sum recorded on declared replicates with a gate PASS     
  p-hat validity (0<=p<=1, p_H<=sum)                             PASS     
  STRUCTURAL n_cons_qual  present, all-NA -- STRUCTURAL on DINA (no consistency screen), NOT a failure
  STRUCTURAL band_n       present, all-NA -- STRUCTURAL on DINA (no consistency screen), NOT a failure
  >> WARNINGS                 : 0 of 2000 rows carry warn_msg; distinct messages: none
  every finiteness column present                                PASS     
  harm products finite on declared replicates with a field block PASS     (1996 rows)
  nine fld_Hc_*_s and nine fld_joint_s_* columns present         PASS     
  fld_Hc_*_s finite on all 1996 filled replicates                PASS     
  fld_joint_s_* finite on all 1996 filled replicates             PASS     
  complement field filled on 1996 of 1996 declared replicates (0 notes)
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
  identity: field-s inverted around the same beta-tilde^c        PASS     max |diff| = 7.11e-15
  identity: Bonferroni harm bound joint == joint_s where draw counts agree PASS     (1996 of 1996 agree; max |diff| 0)
  identity: lo1s = beta-tilde - q95; up1s = beta-tilde^c - q05   PASS     max |diff| = 0
  >> p-hat(H): mean 0.080, share < 0.5: 0.996 | CLASSIFICATION: sens 0.2462 ppv 0.3888 | mean |Hhat| (n_harm) 108.5
  mdsgnb20 comparator on disk                                    PASS     fs_effMaxSG_mr_field_md40_knoise0_n500_nb20_mdsgnb20_combined_1_2000.rds
  --- same-draws: mddina -> mdsgnb20 ---
  [mddina -> mdsgnb20] same sim_id set (2000 rows)               PASS     (comparator 2000 rows)
  [mddina -> mdsgnb20] n_true identical() on all rows            PASS     
  [mddina -> mdsgnb20] oracle columns (H and Hc) <= 1e-08 relative on all rows PASS     (max 0)
  --- same-draws: mdsgnb20 -> mddina ---
  [mdsgnb20 -> mddina] same sim_id set (2000 rows)               PASS     (comparator 2000 rows)
  [mdsgnb20 -> mddina] n_true identical() on all rows            PASS     
  [mdsgnb20 -> mddina] oracle columns (H and Hc) <= 1e-08 relative on all rows PASS     (max 0)
  FS comparator (same cell): declaration 0.9990 | mean |Hhat| 111.6
  size <= 100 MB: dina_effMaxSG_mr_field_md40_knoise0_n500_nb20_mddina_res_1_1000.rds PASS     (789909 B)
  size <= 100 MB: dina_effMaxSG_mr_field_md40_knoise0_n500_nb20_mddina_res_1001_2000.rds PASS     (789677 B)
  size <= 100 MB: dina_effMaxSG_mr_field_md40_knoise0_n500_nb20_mddina_combined_1_2000.rds PASS     (1549947 B)
  timing: fit_mr_secs mean 134.2 median 109.0 p90 272.5 max 343.8 | id_secs (DINA fit) mean 12.05 median 11.11 max 30.66 | fld_H_secs mean 69.8 | fld_Hc_secs mean 7.00

GATE_COUNTS run=65 passed=65 failed=0
```

## md120_n500 — md 120, n 500

- Stem: `dina_effMaxSG_mr_field_md120_knoise0_n500_nb20_mddina`; HEAD before this cell's commit: 1856c107; workers 63; threads 1.
- Knobs: `FS_MD_METHOD=dina FS_MD_FOCUS=effMaxSG FS_MD_NBHD=0.20 FS_MD_FIELD_SCALEC=selected FS_MD_CI=field FS_MD_CAMPAIGN=mddina FS_MD_FB=none FS_MD_MD=120 FS_MD_N=500 FS_MD_WORKERS=63`.
- Seeds: 8316951 + sim_id; batches sim_id 1-1000 and 1001-2000, then combine.
- Cell wall: 10334 s; cumulative render wall 15532 s (ceiling 32939 s).
- Render: WALL_SECONDS=5132 RC=0 PEAK_MB=120815 OUT=dina_effMaxSG_mr_field_md120_knoise0_n500_nb20_mddina_batch_1001_2000.html
- Render: WALL_SECONDS=5166 RC=0 PEAK_MB=121791 OUT=dina_effMaxSG_mr_field_md120_knoise0_n500_nb20_mddina_batch_1_1000.html
- Render: WALL_SECONDS=35 RC=0 PEAK_MB=1112 OUT=dina_effMaxSG_mr_field_md120_knoise0_n500_nb20_mddina_combine_1_2000.html
- Gate counts: `GATE_COUNTS run=65 passed=65 failed=0`.

```

########## GATE 2 (mddina): md120_n500 -- md 120, n 500 ##########
  Every coverage figure of this campaign is conditional on the proposed family.
  bundle: ../mr_md_harm/dina_effMaxSG_mr_field_md120_knoise0_n500_nb20_mddina_d5000/dina_effMaxSG_mr_field_md120_knoise0_n500_nb20_mddina_combined_1_2000.rds
  combined payload on disk                                       PASS     
  --- combine assertions ---
  exactly 2 batch files, res_1_1000 and res_1001_2000            PASS     (2)
  2,000 rows                                                     PASS     (2000)
  combined sim_id == 1:2000                                      PASS     
  batch sim_id sets 1:1000 and 1001:2000                         PASS     (batch1 1000, batch2 1000)
  batch files match the combined bundle on every column          PASS     (175 columns)
  no CONFIG-ERROR replicate                                      PASS     (0)
  --- meta (both batches) ---
  meta: subgroup_method == dina                                  PASS     (dina / dina)
  meta: dina_select_statistic == effect                          PASS     (effect / effect)
  meta: dina_args == list()                                      PASS     (list() / list())
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
  meta: campaign_tag == mddina                                   PASS     (mddina / mddina)
  meta: n_sample == 500                                          PASS     (500 / 500)
  meta: null_cell == FALSE                                       PASS     (FALSE / FALSE)
  meta: effect_threshold == 30                                   PASS     (30 / 30)
  meta: consistency_threshold == 10                              PASS     (10 / 10)
  meta: pkg_version == 0.3.5                                     PASS     (0.3.5 / 0.3.5)
  meta: hostname == pop-os                                       PASS     (pop-os / pop-os)
  meta: n_workers == 63                                          PASS     (63 / 63)
  meta seed_base 8316951 | host pop-os | R 4.6.1 | built_at 2026-09-17 06:11:08 / 2026-09-17 04:45:32
  >> DECLARATION RATE         : 1.0000 (2000 / 2000)
  MR failures on declared replicates <= 40                       PASS     (0)
  >> DINA_PROPOSED_N (>= floor): min 415  q25 5768.75  MED 6161  q75 6369.25  p90 6468.1  max 6712  (mean 5940.0, n 2000); searched 8234-8916
  >> ADMITTED_N (qualified)   : min 397  q25 5341  MED 5809.5  q75 6116.25  p90 6290  max 6645  (mean 5599.1, n 2000); == 0 on 0 row(s)
  >> N_FAMILY (MR's family: proposed candidates with >= n.min members): min 415  med 6161  p90 6468.1  max 6712
  E2 fields filled on every declared replicate                   PASS     (dina_searched_n 2000/2000, dina_proposed_n 2000/2000, dina_tau_min 2000/2000, admitted_n 2000/2000)
  STRUCTURAL E2 fields on non-detections: all NA TRUE (DINA's selection object is absent when nothing is selected)
  every proposed candidate at oriented tau-hat >= 30 (dina_tau_min, declared replicates) PASS     (min 30.0000)
  1 <= admitted_n <= dina_proposed_n on every declared replicate PASS     
  n_family finite on every declared replicate with a gate        PASS     
  p_hat_H and p_hat_sum recorded on declared replicates with a gate PASS     
  p-hat validity (0<=p<=1, p_H<=sum)                             PASS     
  STRUCTURAL n_cons_qual  present, all-NA -- STRUCTURAL on DINA (no consistency screen), NOT a failure
  STRUCTURAL band_n       present, all-NA -- STRUCTURAL on DINA (no consistency screen), NOT a failure
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
  >> p-hat(H): mean 0.083, share < 0.5: 0.999 | CLASSIFICATION: sens 0.6596 ppv 0.7666 | mean |Hhat| (n_harm) 146.6
  mdsgnb20 comparator on disk                                    PASS     fs_effMaxSG_mr_field_md120_knoise0_n500_nb20_mdsgnb20_combined_1_2000.rds
  --- same-draws: mddina -> mdsgnb20 ---
  [mddina -> mdsgnb20] same sim_id set (2000 rows)               PASS     (comparator 2000 rows)
  [mddina -> mdsgnb20] n_true identical() on all rows            PASS     
  [mddina -> mdsgnb20] oracle columns (H and Hc) <= 1e-08 relative on all rows PASS     (max 0)
  --- same-draws: mdsgnb20 -> mddina ---
  [mdsgnb20 -> mddina] same sim_id set (2000 rows)               PASS     (comparator 2000 rows)
  [mdsgnb20 -> mddina] n_true identical() on all rows            PASS     
  [mdsgnb20 -> mddina] oracle columns (H and Hc) <= 1e-08 relative on all rows PASS     (max 0)
  FS comparator (same cell): declaration 1.0000 | mean |Hhat| 152.7
  size <= 100 MB: dina_effMaxSG_mr_field_md120_knoise0_n500_nb20_mddina_res_1_1000.rds PASS     (775501 B)
  size <= 100 MB: dina_effMaxSG_mr_field_md120_knoise0_n500_nb20_mddina_res_1001_2000.rds PASS     (777226 B)
  size <= 100 MB: dina_effMaxSG_mr_field_md120_knoise0_n500_nb20_mddina_combined_1_2000.rds PASS     (1526571 B)
  timing: fit_mr_secs mean 304.9 median 324.0 p90 360.4 max 389.4 | id_secs (DINA fit) mean 21.22 median 22.20 max 30.07 | fld_H_secs mean 149.2 | fld_Hc_secs mean 17.65

GATE_COUNTS run=65 passed=65 failed=0
```

## null_n500 — md null, n 500

- Stem: `dina_effMaxSG_mr_field_mdnull_knoise0_n500_nb20_mddina`; HEAD before this cell's commit: c9f1b948; workers 63; threads 1.
- Knobs: `FS_MD_METHOD=dina FS_MD_FOCUS=effMaxSG FS_MD_NBHD=0.20 FS_MD_FIELD_SCALEC=selected FS_MD_CI=field FS_MD_CAMPAIGN=mddina FS_MD_FB=none FS_MD_MD=null FS_MD_N=500 FS_MD_WORKERS=63`.
- Seeds: 8316951 + sim_id; batches sim_id 1-1000 and 1001-2000, then combine.
- Cell wall: 4120 s; cumulative render wall 19651 s (ceiling 32939 s).
- Render: WALL_SECONDS=1969 RC=0 PEAK_MB=103052 OUT=dina_effMaxSG_mr_field_mdnull_knoise0_n500_nb20_mddina_batch_1001_2000.html
- Render: WALL_SECONDS=2110 RC=0 PEAK_MB=102441 OUT=dina_effMaxSG_mr_field_mdnull_knoise0_n500_nb20_mddina_batch_1_1000.html
- Render: WALL_SECONDS=40 RC=0 PEAK_MB=1106 OUT=dina_effMaxSG_mr_field_mdnull_knoise0_n500_nb20_mddina_combine_1_2000.html
- Gate counts: `GATE_COUNTS run=65 passed=65 failed=0`.

```

########## GATE 2 (mddina): null_n500 -- md null, n 500 ##########
  Every coverage figure of this campaign is conditional on the proposed family.
  bundle: ../mr_md_harm/dina_effMaxSG_mr_field_mdnull_knoise0_n500_nb20_mddina_d5000/dina_effMaxSG_mr_field_mdnull_knoise0_n500_nb20_mddina_combined_1_2000.rds
  combined payload on disk                                       PASS     
  --- combine assertions ---
  exactly 2 batch files, res_1_1000 and res_1001_2000            PASS     (2)
  2,000 rows                                                     PASS     (2000)
  combined sim_id == 1:2000                                      PASS     
  batch sim_id sets 1:1000 and 1001:2000                         PASS     (batch1 1000, batch2 1000)
  batch files match the combined bundle on every column          PASS     (175 columns)
  no CONFIG-ERROR replicate                                      PASS     (0)
  --- meta (both batches) ---
  meta: subgroup_method == dina                                  PASS     (dina / dina)
  meta: dina_select_statistic == effect                          PASS     (effect / effect)
  meta: dina_args == list()                                      PASS     (list() / list())
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
  meta: campaign_tag == mddina                                   PASS     (mddina / mddina)
  meta: n_sample == 500                                          PASS     (500 / 500)
  meta: null_cell == TRUE                                        PASS     (TRUE / TRUE)
  meta: effect_threshold == 30                                   PASS     (30 / 30)
  meta: consistency_threshold == 10                              PASS     (10 / 10)
  meta: pkg_version == 0.3.5                                     PASS     (0.3.5 / 0.3.5)
  meta: hostname == pop-os                                       PASS     (pop-os / pop-os)
  meta: n_workers == 63                                          PASS     (63 / 63)
  meta seed_base 8316951 | host pop-os | R 4.6.1 | built_at 2026-09-17 07:19:39 / 2026-09-17 06:46:50
  >> DECLARATION RATE         : 0.9935 (1987 / 2000)
  MR failures on declared replicates <= 40                       PASS     (0)
  >> DINA_PROPOSED_N (>= floor): min 1  q25 741.5  MED 2046  q75 4608.5  p90 5726.4  max 6617  (mean 2644.1, n 1987); searched 8234-8916
  >> ADMITTED_N (qualified)   : min 1  q25 447  MED 1306  q75 3373.5  p90 4778.2  max 6422  (mean 1959.0, n 1987); == 0 on 0 row(s)
  >> N_FAMILY (MR's family: proposed candidates with >= n.min members): min 1  med 2046  p90 5726.4  max 6617
  E2 fields filled on every declared replicate                   PASS     (dina_searched_n 1987/1987, dina_proposed_n 1987/1987, dina_tau_min 1987/1987, admitted_n 1987/1987)
  STRUCTURAL E2 fields on non-detections: all NA TRUE (DINA's selection object is absent when nothing is selected)
  every proposed candidate at oriented tau-hat >= 30 (dina_tau_min, declared replicates) PASS     (min 30.0000)
  1 <= admitted_n <= dina_proposed_n on every declared replicate PASS     
  n_family finite on every declared replicate with a gate        PASS     
  p_hat_H and p_hat_sum recorded on declared replicates with a gate PASS     
  p-hat validity (0<=p<=1, p_H<=sum)                             PASS     
  STRUCTURAL n_cons_qual  present, all-NA -- STRUCTURAL on DINA (no consistency screen), NOT a failure
  STRUCTURAL band_n       present, all-NA -- STRUCTURAL on DINA (no consistency screen), NOT a failure
  >> WARNINGS                 : 0 of 2000 rows carry warn_msg; distinct messages: none
  every finiteness column present                                PASS     
  harm products finite on declared replicates with a field block PASS     (1987 rows)
  nine fld_Hc_*_s and nine fld_joint_s_* columns present         PASS     
  fld_Hc_*_s finite on all 1987 filled replicates                PASS     
  fld_joint_s_* finite on all 1987 filled replicates             PASS     
  complement field filled on 1987 of 1987 declared replicates (0 notes)
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
  identity: field-s inverted around the same beta-tilde^c        PASS     max |diff| = 7.11e-15
  identity: Bonferroni harm bound joint == joint_s where draw counts agree PASS     (1987 of 1987 agree; max |diff| 0)
  identity: lo1s = beta-tilde - q95; up1s = beta-tilde^c - q05   PASS     max |diff| = 0
  >> p-hat(H): mean 0.092, share < 0.5: 0.990 | CLASSIFICATION: sens NaN ppv 0.0000 | mean |Hhat| (n_harm) 105.1
  mdsgnb20 comparator on disk                                    PASS     fs_effMaxSG_mr_field_mdnull_knoise0_n500_nb20_mdsgnb20_combined_1_2000.rds
  --- same-draws: mddina -> mdsgnb20 ---
  [mddina -> mdsgnb20] same sim_id set (2000 rows)               PASS     (comparator 2000 rows)
  [mddina -> mdsgnb20] n_true identical() on all rows            PASS     
  [mddina -> mdsgnb20] oracle columns (complement only) <= 1e-08 relative on all rows PASS     (max 0)
  --- same-draws: mdsgnb20 -> mddina ---
  [mdsgnb20 -> mddina] same sim_id set (2000 rows)               PASS     (comparator 2000 rows)
  [mdsgnb20 -> mddina] n_true identical() on all rows            PASS     
  [mdsgnb20 -> mddina] oracle columns (complement only) <= 1e-08 relative on all rows PASS     (max 0)
  FS comparator (same cell): declaration 0.9965 | mean |Hhat| 107.5
  size <= 100 MB: dina_effMaxSG_mr_field_mdnull_knoise0_n500_nb20_mddina_res_1_1000.rds PASS     (727316 B)
  size <= 100 MB: dina_effMaxSG_mr_field_mdnull_knoise0_n500_nb20_mddina_res_1001_2000.rds PASS     (726613 B)
  size <= 100 MB: dina_effMaxSG_mr_field_mdnull_knoise0_n500_nb20_mddina_combined_1_2000.rds PASS     (1433600 B)
  timing: fit_mr_secs mean 100.3 median 69.5 p90 228.3 max 306.2 | id_secs (DINA fit) mean 9.32 median 7.10 max 31.53 | fld_H_secs mean 54.4 | fld_Hc_secs mean 4.95

GATE_COUNTS run=65 passed=65 failed=0
```

