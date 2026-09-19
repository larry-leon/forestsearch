# HALT — ACTG175 binary/OR Stage 2

- Halted at (UTC): 2026-09-19T04:32:27Z
- Cell: orfs_or150_n500
- Assertion: gate2.R failed: GATE_COUNTS run=73 passed=72 failed=1
- HEAD at halt: 81e671a9
- Cumulative render wall: 1242 s (ceiling 30000 s)
- Completed cells stay committed; nothing is rolled back.

## Observed values

```

########## GATE 2 (orfs): orfs_or150_n500 -- target_or_h 1.5, n 500, identifier consistency ##########
  FS's candidate family is the prespecified cut grid: the fixed-family condition holds.
  bundle: ../mr_or_harm/fs_effMaxSG_mr_field_or150_n500_nb20_orfs_d5000/fs_effMaxSG_mr_field_or150_n500_nb20_orfs_combined_1_1000.rds
  combined payload on disk                                               PASS     
  --- combine assertions ---
  exactly 1 batch file, res_1_1000                                       PASS     (1)
  1,000 rows                                                             PASS     (1000)
  combined sim_id == 1:1000 exactly                                      PASS     
  batch sim_id set 1:1000                                                PASS     (batch1 1000)
  batch files match the combined bundle on every column                  PASS     (180 columns)
  no CONFIG-ERROR replicate                                              PASS     (0)
  --- meta (the batch) ---
  meta: subgroup_method == consistency                                   PASS     (consistency)
  meta: sg_focus == effMaxSG                                             PASS     (effMaxSG)
  meta: effect_neighborhood == 0.2                                       PASS     (0.2)
  meta: selection_rule == neighborhood                                   PASS     (neighborhood)
  meta: effect_threshold == 0.9                                          PASS     (0.9)
  meta: consistency_threshold == 0.8                                     PASS     (0.8)
  meta: pconsistency == 0.9                                              PASS     (0.9)
  meta: adverse_outcome == TRUE                                          PASS     (TRUE)
  meta: outcome_type == binary                                           PASS     (binary)
  meta: effect_measure == OR                                             PASS     (OR)
  meta: target_or_h == 1.5                                               PASS     (1.5)
  meta: design_tag == or150                                              PASS     (or150)
  meta: dgm_model == alt                                                 PASS     (alt)
  meta: sg_quantile == 0.6285                                            PASS     (0.6285)
  meta: n_super == 100000                                                PASS     (100000)
  meta: eval_seed == 20260628                                            PASS     (20260628)
  meta: ci_method == field                                               PASS     (field)
  meta: mr_draws == 5000                                                 PASS     (5000)
  meta: field_uniform == FALSE                                           PASS     (FALSE)
  meta: field_complement == TRUE                                         PASS     (TRUE)
  meta: field_scale_complement == selected                               PASS     (selected)
  meta: ij_residual == two_term                                          PASS     (two_term)
  meta: return_reselection == TRUE                                       PASS     (TRUE)
  meta: fb_mode == none                                                  PASS     (none)
  meta: seed_base == 8316951                                             PASS     (8316951)
  meta: campaign_tag == orfs                                             PASS     (orfs)
  meta: n_sample == 500                                                  PASS     (500)
  meta: k_random_noise == 0                                              PASS     (0)
  meta: consistency_method == resample                                   PASS     (resample)
  meta: pkg_version == 0.3.5.9000                                        PASS     (0.3.5.9000)
  meta: hostname == pop-os                                               PASS     (pop-os)
  meta: n_workers == 63                                                  PASS     (63)
  meta carries the truths                                                PASS     (marg_H 1.5000000002 | marg_Hc 0.6564077470 | cde_H 1.5414142152 | cde_Hc 0.6313905111 | prev 0.149170)
  meta seed_base 8316951 | seed_scheme pre-generated table indexed by global sim_id | host pop-os | R 4.6.1 | pkg_commit 81e671a9 | built_at 2026-09-18 21:31:51
  --- launch record (combined meta) ---
  meta: n_sims == 1000                                                   PASS     (1000)
  meta: Stage 0 feasibility gate green, not overridden                   PASS     (feasible TRUE, override FALSE, tol 0.05, shares n500=0.015|n750=0.000|n1000=0.000|n2000=0.000)
  meta: the oracle helper assertion passed in the render                 PASS     (TRUE)
  meta: wall clock recorded                                              PASS     (1071 s total; by batch 1071)
  meta: the non-estimable / NA-oracle counts recorded                    PASS     (n_na_oracle_H=0, n_na_oracle_Hc=0, n_nonestimable_H=0, n_nonestimable_Hc=0)
  >> NA-ORACLE (true region)  : H 0 / 1000, Hc 0 / 1000
  >> NON-ESTIMABLE (selected) : H 0, Hc 0 (of the declared replicates)
  >> PREVALENCE(H) / sg_quantile: 0.149170 / 0.62850 | pkg 0.3.5.9000 built 2026-09-18 21:31:51 | workers 63
  >> DECLARATION RATE         : 0.9310 (931 / 1000)
  MR failures on declared replicates <= 40                               PASS     (0)
  >> N_FAMILY (MR's kept family K): min 1850  med 2226  p90 2260  max 2293
  STRUCTURAL n_cons_qual  POPULATED (expected on FS)
  STRUCTURAL band_n       POPULATED (expected on FS)
  n_family finite on every declared replicate with a gate                PASS     
  p_hat_H and p_hat_sum recorded on declared replicates with a gate      PASS     
  p-hat validity (0<=p<=1, p_H<=sum)                                     PASS     
  >> WARNINGS                 : 0 of 1000 rows carry warn_msg; distinct: none
  zero factor-comparison warnings                                        PASS     
  every finiteness column present (incl. C_dagger_* / C_ddagger_*)       PASS     
  harm products finite on declared replicates with a field block         PASS     (931 rows)
  nine fld_Hc_*_s and nine fld_joint_s_* columns present                 PASS     
  fld_Hc_*_s finite on all 931 filled replicates                         PASS     
  fld_joint_s_* finite on all 931 filled replicates                      PASS     
  complement field filled on 931 of 931 declared replicates (0 notes)
  invariant harm  : fld_H_lo1s <= fld_H_est2                             PASS     
  invariant compl : fld_Hc_est2 <= fld_Hc_up1s                           PASS     
  invariant _s    : fld_Hc_lo1s_s <= fld_Hc_up1s_s                       PASS     
  invariant _s    : fld_Hc_lo2s_s <= fld_Hc_hi2s_s                       PASS     
  invariant _s    : fld_Hc_est2_s <= fld_Hc_up1s_s                       PASS     
  invariant joint : bonf_loH <= fld_H_est2                               PASS     
  invariant jointS: fld_Hc_est2_s <= bonf_upHc_s                         PASS     
  invariant IJ    : mr_H_lo <= est <= mr_H_hi                            PASS     
  every ESTIMATE is a positive OR (15 columns)                           PASS     
  no NEGATIVE bound (33 columns)                                         PASS     
  >> DEGENERATE BOUNDS (separation) : none | oracle rows affected: H 0, Hc 0 of 1000
  gamma (joint)   in [0.025, 0.05]                                       PASS     [0.02500, 0.02700]
  gamma (joint-s) in [0.025, 0.05]                                       PASS     [0.02500, 0.02700]
  identity: field-s inverted around the same beta-tilde^c (log scale)    PASS     max |diff| = 2.22e-16
  identity: log(est2) + lambda_mean = log(beta-tilde^c), field and field-s PASS     max |diff| = 2.22e-16
  identity: Bonferroni harm bound joint == joint_s where draw counts agree PASS     (931 of 931 agree; max |diff| 0)
  identity: log(lo1s) = log(beta-tilde) - q95; log(up1s) = log(beta-tilde^c) - q05 PASS     max |diff| = 2.22e-16
  C_dagger_* / C_ddagger_* equal the truth table on every row            PASS     
  >> p-hat(H): mean 0.086, share < 0.5: 0.996 | CLASSIFICATION: sens 0.4068 ppv 0.3256 | mean |Hhat| 99.9
  >> BOUNDS: mean fld_H_lo1s 0.4108 | share >= 1.0 0.020 | mean fld_Hc_up1s_s 0.9821 | share <= 1.0 0.556
  same-draws: this IS the orfs campaign (the reference for the other two).
  checker ran without error                                              **FAIL** object 'f1' not found

GATE_COUNTS run=73 passed=72 failed=1
FAILED: checker ran without error 
```
