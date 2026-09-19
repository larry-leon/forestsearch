# REPORT — ACTG175 binary/OR Gate 2 (per cell)

Task: `dev/tasks/TASK_actg175_binary_campaign_2026-09-17.md` Stage 2. Runner:
`quarto/simulations/actg175/binary_020/scripts_or/run_or.sh`; checker:
`quarto/simulations/actg175/binary_020/scripts_or/gate2.R`.

GRF's and DINA's candidate families are generated from fitted surfaces, so every coverage figure
of the `orgrf` and `ordina` campaigns is coverage of the estimand **conditional on the proposed
family**; FS's family is the prespecified cut grid. Comparisons across the three identifiers are
descriptive.

*(This header was written by hand. The runner emits it only when the record file does not exist,
but its `[ -f ... ]` test runs inside the `>> record` redirection, which has already created the
file — so the header is always skipped. The committed `mdgrf` and `mddina` Gate 2 records have the
same gap; see the findings.)*

## The cell-3 halt, its cause, and the convention that follows from it

**The halt.** Stage 2 halted at 2026-09-18T00:23:34Z on cell 3, `orfs_or150_n500`, batch
`1001_2000`, with `OR positivity violated -- non-positive finite values in: or_H_lo (1)`. Cells 1
and 2 were already committed and are untouched.

**The cause, reproduced exactly.** `sim_id` **1179**, design `or150`, n = 500. The **true** harm
region holds 43 subjects and its **treated arm is all-event** — 17 events, 0 non-events: complete
separation. The study's `.logit_or_ci()` guards only the *overall* ≥ 5 events and ≥ 5 non-events
(`maxeffCons_mr_coverage_sweep_or075.qmd:366–368`), and this subset has 25 and 18, so it fits. The
logistic MLE then diverges — **β̂ = 20.38 with SE = 2608** — so
`exp(20.38 − 1.96 × 2608) = exp(−5092)` **underflows to 0** and the upper bound overflows to `Inf`.
Mathematically the bound is positive; the double cannot hold it. The committed study records
exactly these values and drops such rows from coverage through its `is.finite(lo) & is.finite(hi)`
masks. **The fault was in my guard, not in the data recipe.**

**The patch commits.**

- `1f99dfe4` — the guard splits in two: **estimates** must be strictly positive and stay fatal; a
  **bound** is fatal only when *negative*, while a bound that underflows to 0 or overflows to Inf is
  counted and reported. Applied in the template, `gate2.R` and `smoke_identity.R`.
- `6f80f292` — clears `HALT_or.md` so Stage 2 can resume from cell 3.

**Check 1 — the patch is assertion-only, verified.** `sim_id` 1–20 of the committed cell
`orfs_or075_n500` were re-rendered with the patched template under a separate campaign tag
(`orassert`, so no committed bundle could be written) and compared against the committed rows:
**174 columns compared** (the 6 `*_secs` columns excluded), identical column sets, and a **maximum
relative difference of exactly 0** — NA matching NA and character columns exact. The patch changes
no recorded value. Regression tests alongside it: all four smoke modes PASS, and `gate2.R` passes
**70 / 70** on both committed cells (70 rather than 69 because the estimate/bound split adds a
check).

**The convention for non-convergent fits — one rule, every estimator.** A replicate is
**non-convergent for an estimator** when that estimator's point estimate or either of its two-sided
bounds is non-finite, or its point estimate is ≤ 0. Such rows are excluded from **that** estimator's
coverage, location, spread and bound-location statistics — exactly as the committed study's
finiteness masks exclude them from coverage — and the count is reported per cell × block ×
estimator as `n_nonconvergent_fits` in `or_metrics.csv`, and as a `non-convergent` column in every
table that carries an affected row. **No row is dropped silently**, and no estimator's rows are
dropped on another estimator's account. Coverage was already protected by the finiteness masks; a
**mean** and an **SD** were not, and one separated fit carries an odds ratio of order 10⁸.

**The data recipe is the committed study's, verbatim, and is not touched: every recorded value is
exactly what the recorder wrote. The convention above is a consumer-side convention, applied in
`summary_actg175_or.qmd`.**

**One aborted launch, recorded.** A relaunch was started at 2026-09-18T16:30:51Z and stopped within
minutes, before any cell completed, so that check 1 could be run first. It correctly skipped the two
committed cells; its partial cell-3 output was discarded. Its three heartbeat lines are in
`LOG_or_progress.txt`.

---

## orfs_or075_n500 — identifier consistency, target_or_h 0.75, n 500

- Stem: `fs_effMaxSG_mr_field_or075_n500_nb20_orfs`; HEAD before this cell's commit: 6233870a; workers 63; threads 1.
- Knobs: `FS_OR_METHOD=consistency FS_OR_FOCUS=effMaxSG FS_OR_NBHD=0.20 FS_OR_RULE=neighborhood FS_OR_CI=field FS_OR_FIELD_SCALEC=selected FS_OR_CAMPAIGN=orfs FS_OR_TARGET=0.75 FS_OR_N=500 FS_OR_WORKERS=63`.
- Seeds: the study's pre-generated table indexed by global sim_id (seed_base 8316951); batches sim_id 1-1000 and 1001-2000, then combine.
- Cell wall: 2248 s; cumulative render wall 2248 s (ceiling 105165 s).
- Render: WALL_SECONDS=1129 RC=0 PEAK_MB=78870 OUT=fs_effMaxSG_mr_field_or075_n500_nb20_orfs_batch_1001_2000.html
- Render: WALL_SECONDS=1093 RC=0 PEAK_MB=77454 OUT=fs_effMaxSG_mr_field_or075_n500_nb20_orfs_batch_1_1000.html
- Render: WALL_SECONDS=26 RC=0 PEAK_MB=1233 OUT=fs_effMaxSG_mr_field_or075_n500_nb20_orfs_combine_1_2000.html
- Gate counts: `GATE_COUNTS run=69 passed=69 failed=0`.

```

########## GATE 2 (orfs): orfs_or075_n500 -- target_or_h 0.75, n 500, identifier consistency ##########
  FS's candidate family is the prespecified cut grid: the fixed-family condition holds.
  bundle: ../mr_or_harm/fs_effMaxSG_mr_field_or075_n500_nb20_orfs_d5000/fs_effMaxSG_mr_field_or075_n500_nb20_orfs_combined_1_2000.rds
  combined payload on disk                                               PASS     
  --- combine assertions ---
  exactly 2 batch files, res_1_1000 and res_1001_2000                    PASS     (2)
  2,000 rows                                                             PASS     (2000)
  combined sim_id == 1:2000 exactly                                      PASS     
  batch sim_id sets 1:1000 and 1001:2000                                 PASS     (batch1 1000, batch2 1000)
  batch files match the combined bundle on every column                  PASS     (180 columns)
  no CONFIG-ERROR replicate                                              PASS     (0)
  --- meta (both batches) ---
  meta: subgroup_method == consistency                                   PASS     (consistency / consistency)
  meta: sg_focus == effMaxSG                                             PASS     (effMaxSG / effMaxSG)
  meta: effect_neighborhood == 0.2                                       PASS     (0.2 / 0.2)
  meta: selection_rule == neighborhood                                   PASS     (neighborhood / neighborhood)
  meta: effect_threshold == 0.9                                          PASS     (0.9 / 0.9)
  meta: consistency_threshold == 0.8                                     PASS     (0.8 / 0.8)
  meta: pconsistency == 0.9                                              PASS     (0.9 / 0.9)
  meta: adverse_outcome == TRUE                                          PASS     (TRUE / TRUE)
  meta: outcome_type == binary                                           PASS     (binary / binary)
  meta: effect_measure == OR                                             PASS     (OR / OR)
  meta: target_or_h == 0.75                                              PASS     (0.75 / 0.75)
  meta: design_tag == or075                                              PASS     (or075 / or075)
  meta: dgm_model == alt                                                 PASS     (alt / alt)
  meta: sg_quantile == 0.7                                               PASS     (0.7 / 0.7)
  meta: n_super == 100000                                                PASS     (100000 / 100000)
  meta: eval_seed == 20260628                                            PASS     (20260628 / 20260628)
  meta: ci_method == field                                               PASS     (field / field)
  meta: mr_draws == 5000                                                 PASS     (5000 / 5000)
  meta: field_uniform == FALSE                                           PASS     (FALSE / FALSE)
  meta: field_complement == TRUE                                         PASS     (TRUE / TRUE)
  meta: field_scale_complement == selected                               PASS     (selected / selected)
  meta: ij_residual == two_term                                          PASS     (two_term / two_term)
  meta: return_reselection == TRUE                                       PASS     (TRUE / TRUE)
  meta: fb_mode == none                                                  PASS     (none / none)
  meta: seed_base == 8316951                                             PASS     (8316951 / 8316951)
  meta: campaign_tag == orfs                                             PASS     (orfs / orfs)
  meta: n_sample == 500                                                  PASS     (500 / 500)
  meta: k_random_noise == 0                                              PASS     (0 / 0)
  meta: consistency_method == resample                                   PASS     (resample / resample)
  meta: pkg_version == 0.3.5                                             PASS     (0.3.5 / 0.3.5)
  meta: hostname == pop-os                                               PASS     (pop-os / pop-os)
  meta: n_workers == 63                                                  PASS     (63 / 63)
  meta carries the truths                                                PASS     (marg_H 0.7499999955 | marg_Hc 0.6560116720 | cde_H 0.7321189340 | cde_Hc 0.6313905111 | prev 0.096320)
  meta seed_base 8316951 | seed_scheme pre-generated table indexed by global sim_id | host pop-os | R 4.6.1 | pkg_commit 6233870a | built_at 2026-09-17 13:59:37 / 2026-09-17 13:40:48
  >> DECLARATION RATE         : 0.7615 (1523 / 2000)
  MR failures on declared replicates <= 40                               PASS     (0)
  >> N_FAMILY (MR's kept family K): min 1860  med 2224  p90 2259  max 2293
  STRUCTURAL n_cons_qual  POPULATED (expected on FS)
  STRUCTURAL band_n       POPULATED (expected on FS)
  n_family finite on every declared replicate with a gate                PASS     
  p_hat_H and p_hat_sum recorded on declared replicates with a gate      PASS     
  p-hat validity (0<=p<=1, p_H<=sum)                                     PASS     
  >> WARNINGS                 : 0 of 2000 rows carry warn_msg; distinct: none
  zero factor-comparison warnings                                        PASS     
  every finiteness column present (incl. C_dagger_* / C_ddagger_*)       PASS     
  harm products finite on declared replicates with a field block         PASS     (1523 rows)
  nine fld_Hc_*_s and nine fld_joint_s_* columns present                 PASS     
  fld_Hc_*_s finite on all 1523 filled replicates                        PASS     
  fld_joint_s_* finite on all 1523 filled replicates                     PASS     
  complement field filled on 1523 of 1523 declared replicates (0 notes)
  invariant harm  : fld_H_lo1s <= fld_H_est2                             PASS     
  invariant compl : fld_Hc_est2 <= fld_Hc_up1s                           PASS     
  invariant _s    : fld_Hc_lo1s_s <= fld_Hc_up1s_s                       PASS     
  invariant _s    : fld_Hc_lo2s_s <= fld_Hc_hi2s_s                       PASS     
  invariant _s    : fld_Hc_est2_s <= fld_Hc_up1s_s                       PASS     
  invariant joint : bonf_loH <= fld_H_est2                               PASS     
  invariant jointS: fld_Hc_est2_s <= bonf_upHc_s                         PASS     
  invariant IJ    : mr_H_lo <= est <= mr_H_hi                            PASS     
  every bound is an OR, so positive (46 columns)                         PASS     
  gamma (joint)   in [0.025, 0.05]                                       PASS     [0.02500, 0.02700]
  gamma (joint-s) in [0.025, 0.05]                                       PASS     [0.02500, 0.02700]
  identity: field-s inverted around the same beta-tilde^c (log scale)    PASS     max |diff| = 3.33e-16
  identity: log(est2) + lambda_mean = log(beta-tilde^c), field and field-s PASS     max |diff| = 2.22e-16
  identity: Bonferroni harm bound joint == joint_s where draw counts agree PASS     (1523 of 1523 agree; max |diff| 0)
  identity: log(lo1s) = log(beta-tilde) - q95; log(up1s) = log(beta-tilde^c) - q05 PASS     max |diff| = 2.22e-16
  C_dagger_* / C_ddagger_* equal the truth table on every row            PASS     
  >> p-hat(H): mean 0.077, share < 0.5: 0.998 | CLASSIFICATION: sens 0.2451 ppv 0.1261 | mean |Hhat| 100.6
  >> BOUNDS: mean fld_H_lo1s 0.3408 | share >= 1.0 0.003 | mean fld_Hc_up1s_s 0.9181 | share <= 1.0 0.710
  same-draws: this IS the orfs campaign (the reference for the other two).
  size <= 100 MB: fs_effMaxSG_mr_field_or075_n500_nb20_orfs_res_1_1000.rds PASS     (644173 B)
  size <= 100 MB: fs_effMaxSG_mr_field_or075_n500_nb20_orfs_res_1001_2000.rds PASS     (638339 B)
  size <= 100 MB: fs_effMaxSG_mr_field_or075_n500_nb20_orfs_combined_1_2000.rds PASS     (1258996 B)
  timing: fit_mr_secs mean 45.0 median 54.7 p90 64.7 max 76.1 | id_secs mean 8.39 median 8.26 max 16.85 | fld_H_secs mean 25.7 | fld_Hc_secs mean 3.69

GATE_COUNTS run=69 passed=69 failed=0
```

## orfs_or075_n2000 — identifier consistency, target_or_h 0.75, n 2000

- Stem: `fs_effMaxSG_mr_field_or075_n2000_nb20_orfs`; HEAD before this cell's commit: 59b144ea; workers 63; threads 1.
- Knobs: `FS_OR_METHOD=consistency FS_OR_FOCUS=effMaxSG FS_OR_NBHD=0.20 FS_OR_RULE=neighborhood FS_OR_CI=field FS_OR_FIELD_SCALEC=selected FS_OR_CAMPAIGN=orfs FS_OR_TARGET=0.75 FS_OR_N=2000 FS_OR_WORKERS=63`.
- Seeds: the study's pre-generated table indexed by global sim_id (seed_base 8316951); batches sim_id 1-1000 and 1001-2000, then combine.
- Cell wall: 9986 s; cumulative render wall 12233 s (ceiling 105165 s).
- Render: WALL_SECONDS=4939 RC=0 PEAK_MB=122180 OUT=fs_effMaxSG_mr_field_or075_n2000_nb20_orfs_batch_1001_2000.html
- Render: WALL_SECONDS=5021 RC=0 PEAK_MB=123052 OUT=fs_effMaxSG_mr_field_or075_n2000_nb20_orfs_batch_1_1000.html
- Render: WALL_SECONDS=25 RC=0 PEAK_MB=1289 OUT=fs_effMaxSG_mr_field_or075_n2000_nb20_orfs_combine_1_2000.html
- Gate counts: `GATE_COUNTS run=69 passed=69 failed=0`.

```

########## GATE 2 (orfs): orfs_or075_n2000 -- target_or_h 0.75, n 2000, identifier consistency ##########
  FS's candidate family is the prespecified cut grid: the fixed-family condition holds.
  bundle: ../mr_or_harm/fs_effMaxSG_mr_field_or075_n2000_nb20_orfs_d5000/fs_effMaxSG_mr_field_or075_n2000_nb20_orfs_combined_1_2000.rds
  combined payload on disk                                               PASS     
  --- combine assertions ---
  exactly 2 batch files, res_1_1000 and res_1001_2000                    PASS     (2)
  2,000 rows                                                             PASS     (2000)
  combined sim_id == 1:2000 exactly                                      PASS     
  batch sim_id sets 1:1000 and 1001:2000                                 PASS     (batch1 1000, batch2 1000)
  batch files match the combined bundle on every column                  PASS     (180 columns)
  no CONFIG-ERROR replicate                                              PASS     (0)
  --- meta (both batches) ---
  meta: subgroup_method == consistency                                   PASS     (consistency / consistency)
  meta: sg_focus == effMaxSG                                             PASS     (effMaxSG / effMaxSG)
  meta: effect_neighborhood == 0.2                                       PASS     (0.2 / 0.2)
  meta: selection_rule == neighborhood                                   PASS     (neighborhood / neighborhood)
  meta: effect_threshold == 0.9                                          PASS     (0.9 / 0.9)
  meta: consistency_threshold == 0.8                                     PASS     (0.8 / 0.8)
  meta: pconsistency == 0.9                                              PASS     (0.9 / 0.9)
  meta: adverse_outcome == TRUE                                          PASS     (TRUE / TRUE)
  meta: outcome_type == binary                                           PASS     (binary / binary)
  meta: effect_measure == OR                                             PASS     (OR / OR)
  meta: target_or_h == 0.75                                              PASS     (0.75 / 0.75)
  meta: design_tag == or075                                              PASS     (or075 / or075)
  meta: dgm_model == alt                                                 PASS     (alt / alt)
  meta: sg_quantile == 0.7                                               PASS     (0.7 / 0.7)
  meta: n_super == 100000                                                PASS     (100000 / 100000)
  meta: eval_seed == 20260628                                            PASS     (20260628 / 20260628)
  meta: ci_method == field                                               PASS     (field / field)
  meta: mr_draws == 5000                                                 PASS     (5000 / 5000)
  meta: field_uniform == FALSE                                           PASS     (FALSE / FALSE)
  meta: field_complement == TRUE                                         PASS     (TRUE / TRUE)
  meta: field_scale_complement == selected                               PASS     (selected / selected)
  meta: ij_residual == two_term                                          PASS     (two_term / two_term)
  meta: return_reselection == TRUE                                       PASS     (TRUE / TRUE)
  meta: fb_mode == none                                                  PASS     (none / none)
  meta: seed_base == 8316951                                             PASS     (8316951 / 8316951)
  meta: campaign_tag == orfs                                             PASS     (orfs / orfs)
  meta: n_sample == 2000                                                 PASS     (2000 / 2000)
  meta: k_random_noise == 0                                              PASS     (0 / 0)
  meta: consistency_method == resample                                   PASS     (resample / resample)
  meta: pkg_version == 0.3.5                                             PASS     (0.3.5 / 0.3.5)
  meta: hostname == pop-os                                               PASS     (pop-os / pop-os)
  meta: n_workers == 63                                                  PASS     (63 / 63)
  meta carries the truths                                                PASS     (marg_H 0.7499999955 | marg_Hc 0.6560116720 | cde_H 0.7321189340 | cde_Hc 0.6313905111 | prev 0.096320)
  meta seed_base 8316951 | seed_scheme pre-generated table indexed by global sim_id | host pop-os | R 4.6.1 | pkg_commit 59b144ea | built_at 2026-09-17 16:46:00 / 2026-09-17 15:23:40
  >> DECLARATION RATE         : 0.8950 (1790 / 2000)
  MR failures on declared replicates <= 40                               PASS     (0)
  >> N_FAMILY (MR's kept family K): min 2525  med 2985  p90 2991  max 3001
  STRUCTURAL n_cons_qual  POPULATED (expected on FS)
  STRUCTURAL band_n       POPULATED (expected on FS)
  n_family finite on every declared replicate with a gate                PASS     
  p_hat_H and p_hat_sum recorded on declared replicates with a gate      PASS     
  p-hat validity (0<=p<=1, p_H<=sum)                                     PASS     
  >> WARNINGS                 : 0 of 2000 rows carry warn_msg; distinct: none
  zero factor-comparison warnings                                        PASS     
  every finiteness column present (incl. C_dagger_* / C_ddagger_*)       PASS     
  harm products finite on declared replicates with a field block         PASS     (1790 rows)
  nine fld_Hc_*_s and nine fld_joint_s_* columns present                 PASS     
  fld_Hc_*_s finite on all 1790 filled replicates                        PASS     
  fld_joint_s_* finite on all 1790 filled replicates                     PASS     
  complement field filled on 1790 of 1790 declared replicates (0 notes)
  invariant harm  : fld_H_lo1s <= fld_H_est2                             PASS     
  invariant compl : fld_Hc_est2 <= fld_Hc_up1s                           PASS     
  invariant _s    : fld_Hc_lo1s_s <= fld_Hc_up1s_s                       PASS     
  invariant _s    : fld_Hc_lo2s_s <= fld_Hc_hi2s_s                       PASS     
  invariant _s    : fld_Hc_est2_s <= fld_Hc_up1s_s                       PASS     
  invariant joint : bonf_loH <= fld_H_est2                               PASS     
  invariant jointS: fld_Hc_est2_s <= bonf_upHc_s                         PASS     
  invariant IJ    : mr_H_lo <= est <= mr_H_hi                            PASS     
  every bound is an OR, so positive (46 columns)                         PASS     
  gamma (joint)   in [0.025, 0.05]                                       PASS     [0.02500, 0.02700]
  gamma (joint-s) in [0.025, 0.05]                                       PASS     [0.02500, 0.02700]
  identity: field-s inverted around the same beta-tilde^c (log scale)    PASS     max |diff| = 2.22e-16
  identity: log(est2) + lambda_mean = log(beta-tilde^c), field and field-s PASS     max |diff| = 2.22e-16
  identity: Bonferroni harm bound joint == joint_s where draw counts agree PASS     (1790 of 1790 agree; max |diff| 0)
  identity: log(lo1s) = log(beta-tilde) - q95; log(up1s) = log(beta-tilde^c) - q05 PASS     max |diff| = 2.22e-16
  C_dagger_* / C_ddagger_* equal the truth table on every row            PASS     
  >> p-hat(H): mean 0.103, share < 0.5: 0.996 | CLASSIFICATION: sens 0.1057 ppv 0.1616 | mean |Hhat| 141.5
  >> BOUNDS: mean fld_H_lo1s 0.3377 | share >= 1.0 0.003 | mean fld_Hc_up1s_s 0.7698 | share <= 1.0 0.998
  same-draws: this IS the orfs campaign (the reference for the other two).
  size <= 100 MB: fs_effMaxSG_mr_field_or075_n2000_nb20_orfs_res_1_1000.rds PASS     (711595 B)
  size <= 100 MB: fs_effMaxSG_mr_field_or075_n2000_nb20_orfs_res_1001_2000.rds PASS     (702304 B)
  size <= 100 MB: fs_effMaxSG_mr_field_or075_n2000_nb20_orfs_combined_1_2000.rds PASS     (1392267 B)
  timing: fit_mr_secs mean 278.1 median 319.7 p90 340.4 max 372.6 | id_secs mean 26.07 median 26.76 max 42.93 | fld_H_secs mean 70.3 | fld_Hc_secs mean 28.14

GATE_COUNTS run=69 passed=69 failed=0
```

## orfs_or150_n500 — identifier consistency, target_or_h 1.5, n 500

- Stem: `fs_effMaxSG_mr_field_or150_n500_nb20_orfs`; HEAD before this cell's commit: 236ef82a; workers 63; threads 1.
- Knobs: `FS_OR_METHOD=consistency FS_OR_FOCUS=effMaxSG FS_OR_NBHD=0.20 FS_OR_RULE=neighborhood FS_OR_CI=field FS_OR_FIELD_SCALEC=selected FS_OR_CAMPAIGN=orfs FS_OR_TARGET=1.5 FS_OR_N=500 FS_OR_WORKERS=63`.
- Seeds: the study's pre-generated table indexed by global sim_id (seed_base 8316951); one batch, sim_id 1-1000, then combine.
- Replicates: 1000 (TASK_binary_launch_v2_2026-09-18; the superseded 2 x 1,000 layout is not used).
- Package: R 4.6.1; ; 2026-09-19 03:54:02 UTC; unix.
- Cell wall: 1237 s; cumulative render wall 1237 s (ceiling 30000 s).
- Render: WALL_SECONDS=1054 RC=1 PEAK_MB=80219 OUT=fs_effMaxSG_mr_field_or150_n500_nb20_orfs_batch_1001_2000.html
- Render: WALL_SECONDS=1216 RC=0 PEAK_MB=79041 OUT=fs_effMaxSG_mr_field_or150_n500_nb20_orfs_batch_1_1000.html
- Render: WALL_SECONDS=21 RC=0 PEAK_MB=1255 OUT=fs_effMaxSG_mr_field_or150_n500_nb20_orfs_combine_1_1000.html
- Gate counts: `GATE_COUNTS run=74 passed=74 failed=0`.

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
  meta seed_base 8316951 | seed_scheme pre-generated table indexed by global sim_id | host pop-os | R 4.6.1 | pkg_commit 236ef82a | built_at 2026-09-18 21:53:45
  --- launch record (combined meta) ---
  meta: n_sims == 1000                                                   PASS     (1000)
  meta: Stage 0 feasibility gate green, not overridden                   PASS     (feasible TRUE, override FALSE, tol 0.05, shares n500=0.015|n750=0.000|n1000=0.000|n2000=0.000)
  meta: the oracle helper assertion passed in the render                 PASS     (TRUE)
  meta: wall clock recorded                                              PASS     (1070 s total; by batch 1070)
  meta: the non-estimable / NA-oracle counts recorded                    PASS     (n_na_oracle_H=0, n_na_oracle_Hc=0, n_nonestimable_H=0, n_nonestimable_Hc=0)
  >> NA-ORACLE (true region)  : H 0 / 1000, Hc 0 / 1000
  >> NON-ESTIMABLE (selected) : H 0, Hc 0 (of the declared replicates)
  >> PREVALENCE(H) / sg_quantile: 0.149170 / 0.62850 | pkg 0.3.5.9000 built 2026-09-18 21:53:45 | workers 63
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
  size <= 100 MB: fs_effMaxSG_mr_field_or150_n500_nb20_orfs_res_1_1000.rds PASS     (734652 B)
  size <= 100 MB: fs_effMaxSG_mr_field_or150_n500_nb20_orfs_combined_1_1000.rds PASS     (734706 B)
  timing: fit_mr_secs mean 57.0 median 60.8 p90 69.6 max 83.0 | id_secs mean 9.95 median 9.72 max 19.23 | fld_H_secs mean 27.6 | fld_Hc_secs mean 3.76

GATE_COUNTS run=74 passed=74 failed=0
```

## orfs_or075_n500 — identifier consistency, target_or_h 0.75, n 500

- Stem: `fs_effMaxSG_mr_field_or075_n500_nb20_orfs`; HEAD before this cell's commit: 13d5b6c4; workers 63; threads 1.
- Knobs: `FS_OR_METHOD=consistency FS_OR_FOCUS=effMaxSG FS_OR_NBHD=0.20 FS_OR_RULE=neighborhood FS_OR_CI=field FS_OR_FIELD_SCALEC=selected FS_OR_CAMPAIGN=orfs FS_OR_TARGET=0.75 FS_OR_N=500 FS_OR_WORKERS=63`.
- Seeds: the study's pre-generated table indexed by global sim_id (seed_base 8316951); one batch, sim_id 1-1000, then combine.
- Replicates: 1000 (TASK_binary_launch_v2_2026-09-18; the superseded 2 x 1,000 layout is not used).
- Package: R 4.6.1; ; 2026-09-19 03:54:02 UTC; unix.
- Cell wall: 1105 s; cumulative render wall 2341 s (ceiling 30000 s).
- Render: WALL_SECONDS=1129 RC=0 PEAK_MB=78870 OUT=fs_effMaxSG_mr_field_or075_n500_nb20_orfs_batch_1001_2000.html
- Render: WALL_SECONDS=1079 RC=0 PEAK_MB=78865 OUT=fs_effMaxSG_mr_field_or075_n500_nb20_orfs_batch_1_1000.html
- Render: WALL_SECONDS=25 RC=0 PEAK_MB=1217 OUT=fs_effMaxSG_mr_field_or075_n500_nb20_orfs_combine_1_1000.html
- Render: WALL_SECONDS=26 RC=0 PEAK_MB=1233 OUT=fs_effMaxSG_mr_field_or075_n500_nb20_orfs_combine_1_2000.html
- Gate counts: `GATE_COUNTS run=74 passed=74 failed=0`.

```

########## GATE 2 (orfs): orfs_or075_n500 -- target_or_h 0.75, n 500, identifier consistency ##########
  FS's candidate family is the prespecified cut grid: the fixed-family condition holds.
  bundle: ../mr_or_harm/fs_effMaxSG_mr_field_or075_n500_nb20_orfs_d5000/fs_effMaxSG_mr_field_or075_n500_nb20_orfs_combined_1_1000.rds
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
  meta: target_or_h == 0.75                                              PASS     (0.75)
  meta: design_tag == or075                                              PASS     (or075)
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
  meta carries the truths                                                PASS     (marg_H 0.7499999940 | marg_Hc 0.6564077470 | cde_H 0.7345268290 | cde_Hc 0.6313905111 | prev 0.149170)
  meta seed_base 8316951 | seed_scheme pre-generated table indexed by global sim_id | host pop-os | R 4.6.1 | pkg_commit 13d5b6c4 | built_at 2026-09-18 22:12:02
  --- launch record (combined meta) ---
  meta: n_sims == 1000                                                   PASS     (1000)
  meta: Stage 0 feasibility gate green, not overridden                   PASS     (feasible TRUE, override FALSE, tol 0.05, shares n500=0.015|n750=0.000|n1000=0.000|n2000=0.000)
  meta: the oracle helper assertion passed in the render                 PASS     (TRUE)
  meta: wall clock recorded                                              PASS     (947 s total; by batch 947)
  meta: the non-estimable / NA-oracle counts recorded                    PASS     (n_na_oracle_H=0, n_na_oracle_Hc=0, n_nonestimable_H=0, n_nonestimable_Hc=0)
  >> NA-ORACLE (true region)  : H 0 / 1000, Hc 0 / 1000
  >> NON-ESTIMABLE (selected) : H 0, Hc 0 (of the declared replicates)
  >> PREVALENCE(H) / sg_quantile: 0.149170 / 0.62850 | pkg 0.3.5.9000 built 2026-09-18 22:12:02 | workers 63
  >> DECLARATION RATE         : 0.7810 (781 / 1000)
  MR failures on declared replicates <= 40                               PASS     (0)
  >> N_FAMILY (MR's kept family K): min 1851  med 2225  p90 2259  max 2293
  STRUCTURAL n_cons_qual  POPULATED (expected on FS)
  STRUCTURAL band_n       POPULATED (expected on FS)
  n_family finite on every declared replicate with a gate                PASS     
  p_hat_H and p_hat_sum recorded on declared replicates with a gate      PASS     
  p-hat validity (0<=p<=1, p_H<=sum)                                     PASS     
  >> WARNINGS                 : 0 of 1000 rows carry warn_msg; distinct: none
  zero factor-comparison warnings                                        PASS     
  every finiteness column present (incl. C_dagger_* / C_ddagger_*)       PASS     
  harm products finite on declared replicates with a field block         PASS     (781 rows)
  nine fld_Hc_*_s and nine fld_joint_s_* columns present                 PASS     
  fld_Hc_*_s finite on all 781 filled replicates                         PASS     
  fld_joint_s_* finite on all 781 filled replicates                      PASS     
  complement field filled on 781 of 781 declared replicates (0 notes)
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
  identity: Bonferroni harm bound joint == joint_s where draw counts agree PASS     (781 of 781 agree; max |diff| 0)
  identity: log(lo1s) = log(beta-tilde) - q95; log(up1s) = log(beta-tilde^c) - q05 PASS     max |diff| = 2.22e-16
  C_dagger_* / C_ddagger_* equal the truth table on every row            PASS     
  >> p-hat(H): mean 0.076, share < 0.5: 0.997 | CLASSIFICATION: sens 0.2401 ppv 0.1861 | mean |Hhat| 101.7
  >> BOUNDS: mean fld_H_lo1s 0.3434 | share >= 1.0 0.004 | mean fld_Hc_up1s_s 0.9213 | share <= 1.0 0.703
  same-draws: this IS the orfs campaign (the reference for the other two).
  size <= 100 MB: fs_effMaxSG_mr_field_or075_n500_nb20_orfs_res_1_1000.rds PASS     (653065 B)
  size <= 100 MB: fs_effMaxSG_mr_field_or075_n500_nb20_orfs_combined_1_1000.rds PASS     (653120 B)
  timing: fit_mr_secs mean 46.6 median 56.2 p90 65.3 max 77.0 | id_secs mean 8.77 median 8.70 max 16.86 | fld_H_secs mean 25.8 | fld_Hc_secs mean 3.75

GATE_COUNTS run=74 passed=74 failed=0
```

