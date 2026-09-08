# REPORT — Continuous/MD field campaign `mdf1`: Gate 2 records (per cell)

Date: 2026-09-07. Branch `feature/glm-extension-mac`; forestsearch 0.3.5; template `sim_fs_maxeffCons_mr_field_md_template.qmd`; 13 workers; `ci_method = "field"`, `include_complement = TRUE`, `field_complement = TRUE`, `return_reselection = TRUE`, `ij_residual = "two_term"`, 5,000 draws, field R_out/R_in = 1000/500 (package defaults). Each cell: batch 1 = `sim_id` 1–1000 on the committed seeds (`seed_base + sim_id`, L'Ecuyer-CMRG), batch 2 = 1001–2000 under the same scheme, then a combine render. Gate 2 checker: the session's `gate2_check.R` (its verbatim output is quoted per cell). Pre-authorization (M-5): 4 h cumulative wall ceiling, 1.5 h hard timeout per render, stop if a render exceeds 1.5× its Gate 1 projection.

**Pairing-proof convention.** Pre-existing columns (the twin's 59 minus timings, messages and the FB columns) compared by `sim_id` on batch 1 against the committed bundle at ≤ 1e-8 relative, rule strings as term sets. Every row that is not identical is enumerated and classified: *selection flip* (the sample membership differs: n_harm or naive/oracle estimates differ), *same selection, MR numerics differ*, *label tie with target move* (identical sample membership, the cut label differs, and the exact super-population target β(Ĥ) attached to the label moves), *pure label tie* (identical membership, numerics identical, label only). Enumerated rows are excluded from the pairing proof, never from the results.

## Cell 1 — md40 n = 500 (priority 1; FB joined on batch 1)

Wall: batch 1 **1,189 s** (19.8 min; projection 21 min, stop threshold 31.5 min), batch 2 1,209 s, combine 24 s; cumulative 2,422 s. **Peak memory on the first render: 17.5 GB** (summed RSS of all R/Quarto processes, sampled every 5 s; 36 GB machine). Per replicate: fit + MR + field 14.9 s mean (max 19.6), batch 2 15.0 s. Detection 1,998 of 2,000. FB joined on `sim_id` 1–100 from the committed 0.2.0 FB bundle under the 1e-12 naive-identity licence (no skips).

Gate 2 output (verbatim):

```
== GATE 2: md=40 n=500 ==
  [PASS] exists: fs_maxeffCons_mr_field_md40_knoise0_n500_mdf1_res_1_1000.rds
  [PASS] exists: fs_maxeffCons_mr_field_md40_knoise0_n500_mdf1_res_1001_2000.rds
  [PASS] exists: fs_maxeffCons_mr_field_md40_knoise0_n500_mdf1_combined_1_2000.rds
  [PASS] save guard: all three bundle paths untracked at save time (guard cannot have been bypassed)
  [PASS] completeness: sim_id 1-1000, 1001-2000, combined 1-2000
  [PASS] no CONFIG-ERROR rows (status table: DETECTED 1998, NO-DETECTION 2)
  pairing proof: 986 of 1000 sim_ids identical on 38 pre-existing numeric cols (max rel diff among identical rows 4.32e-10); rule-string-order-only differences on 12 sim_ids
  ENUMERATED FLIPS (14): 8, 9, 92, 161, 206, 261, 267, 296, 412, 448, 809, 817, 850, 946
    sim 8 -> label tie, no numeric consequence
    sim 9 -> label tie, super-population target moves
    sim 92 -> label tie, super-population target moves
    sim 161 -> label tie, super-population target moves
    sim 206 -> label tie, no numeric consequence
    sim 261 -> label tie, no numeric consequence
    sim 267 -> same selection, MR numerics differ
    sim 296 -> label tie, no numeric consequence
    sim 412 -> label tie, super-population target moves
    sim 448 -> label tie, super-population target moves
    sim 809 -> label tie, no numeric consequence
    sim 817 -> label tie, no numeric consequence
    sim 850 -> label tie, no numeric consequence
    sim 946 -> label tie, no numeric consequence
  classification: 0 selection flips; 1 MR-numerics; 5 label ties with target move; 8 pure label ties
  [PASS] pairing proof: all 986 non-flip rows identical (<= 1e-8); 14 enumerated and excluded (0 selection flips)
  [PASS] new columns present
  finite share on 1998 detected reps: min 1.0000 (fld_H_est2); field notes set on 0 (H) / 0 (Hc)
  [PASS] fields finite on >= 99% of detected replicates
  [PASS] interval invariants lo <= hi
  [PASS] bound identities (max abs 0.0e+00)
  [PASS] gamma in [0.025, 0.05] (range 0.025-0.027; mean 0.0251)
  [PASS] p_hat(Hhat) finite in [0,1] (mean 0.155, share < 0.5: 0.990)
  FB joined on 100 sim_ids (fs_maxeffCons_fb_mr_md40_knoise0_n500_res_1_100.rds)
  meta: pkg 0.3.5 | workers 13/13 | ci field | complement TRUE | ij_residual two_term | built 2026-09-07 13:27:23
  timing: batch1 fit+MR mean 14.9 s (max 19.6), batch2 15.0 s; detection 0.999
GATE 2 md=40 n=500: PASS
```

Enumerated rows, with the committed and current rule strings (the full list is in `gate2_flips.txt` beside the bundles): sims 8, 206, 261, 296, 809, 817, 850, 946 are pure label ties (`{str2}` ≡ `!{preanti <= 0}` on the sample; `karnof <= 90` ≡ `karnof <= 95` since karnof takes values 70/80/90/100; identical estimates to 1e-13); sims 9, 92, 161, 412, 448 are the `{preanti <= 0}` / `!{str2}` label tie, identical on the sample but `str2 == 1` and `preanti > 0` disagree on 7 of 1,083 patients, so β(Ĥ)/β(Ĥᶜ) attached to the label on the super-population move by ≤ 0.07 MD units (e.g. sim 92: −32.069 vs −31.997); sim 267 has the identical rule and naive estimate with the MR columns differing at 2.4e-4 relative (cross-vintage; one replicate). **No selection flip in the survival sense (membership identical on all 1,000).** The tie-break between duplicate-membership labels is the 0.2.2 → 0.3.5 difference.

**Gate 2 cell 1: PASS.** Committed: `mr_md_harm/fs_maxeffCons_mr_field_md40_knoise0_n500_mdf1_d5000/` (two batch bundles, the combined bundle, `gate2_flips.txt`) and the rendered combine document `fs_maxeffCons_mr_field_md40_knoise0_n500_mdf1_combine_1_2000.html`.

## Cell 2 — md120 n = 500 (priority 2; no FB)

Wall: batch 1 **1,371 s** (22.9 min; projection 21 min, stop threshold 31.5 min), batch 2 1,344 s, combine 20 s; cumulative 5,157 s (86 min). Per replicate: fit + MR + field 17.3 s mean (max 20.5); the md120 cell's field costs more than md40's (more candidates admitted at the larger effect). Detection 2,000 of 2,000.

Gate 2 output (verbatim):

```
== GATE 2: md=120 n=500 ==
  [PASS] exists: fs_maxeffCons_mr_field_md120_knoise0_n500_mdf1_res_1_1000.rds
  [PASS] exists: fs_maxeffCons_mr_field_md120_knoise0_n500_mdf1_res_1001_2000.rds
  [PASS] exists: fs_maxeffCons_mr_field_md120_knoise0_n500_mdf1_combined_1_2000.rds
  [PASS] save guard: all three bundle paths untracked at save time (guard cannot have been bypassed)
  [PASS] completeness: sim_id 1-1000, 1001-2000, combined 1-2000
  [PASS] no CONFIG-ERROR rows (status table: DETECTED 2000)
  pairing proof: 992 of 1000 sim_ids identical on 38 pre-existing numeric cols (max rel diff among identical rows 1.47e-11); rule-string-order-only differences on 9 sim_ids
  ENUMERATED FLIPS (8): 178, 267, 412, 567, 650, 668, 851, 880
    sim 178 -> label tie, no numeric consequence
    sim 267 -> same selection, MR numerics differ
    sim 412 -> label tie, super-population target moves
    sim 567 -> label tie, no numeric consequence
    sim 650 -> label tie, super-population target moves
    sim 668 -> same selection, MR numerics differ
    sim 851 -> label tie, super-population target moves
    sim 880 -> label tie, no numeric consequence
  classification: 0 selection flips; 2 MR-numerics; 3 label ties with target move; 3 pure label ties
  [PASS] pairing proof: all 992 non-flip rows identical (<= 1e-8); 8 enumerated and excluded (0 selection flips)
  [PASS] new columns present
  finite share on 2000 detected reps: min 1.0000 (fld_H_est2); field notes set on 0 (H) / 0 (Hc)
  [PASS] fields finite on >= 99% of detected replicates
  [PASS] interval invariants lo <= hi
  [PASS] bound identities (max abs 0.0e+00)
  [PASS] gamma in [0.025, 0.05] (range 0.025-0.028; mean 0.0252)
  [PASS] p_hat(Hhat) finite in [0,1] (mean 0.189, share < 0.5: 0.974)
  meta: pkg 0.3.5 | workers 13/13 | ci field | complement TRUE | ij_residual two_term | built 2026-09-07 14:14:02
  timing: batch1 fit+MR mean 17.3 s (max 20.5), batch2 17.0 s; detection 1.000
GATE 2 md=120 n=500: PASS
```

Enumerated rows: the same two mechanisms as cell 1 (label ties on duplicate-membership cuts, five of them the `{preanti <= 0}` / `!{str2}` pair whose super-population target moves; sims 267 and 668 with identical selection and MR columns differing at the 1e-4 level). **No selection flip.** Committed: the three bundles, `gate2_flips.txt`, and `fs_maxeffCons_mr_field_md120_knoise0_n500_mdf1_combine_1_2000.html`.

## Cell 3 — null n = 500 (priority 3; no FB)

Wall: batch 1 **1,200 s** (20.0 min; projection 20 min, stop threshold 30 min), batch 2 1,163 s, combine 23 s; cumulative 7,543 s (126 min). Per replicate: fit + MR + field 15.0 s mean (max 20.9). Detection 1,993 of 2,000 (the null DGM still yields a screened winner on almost every replicate: the consistency screen is non-discriminating at these thresholds, as the twin's D2 note records). Under the null, β(Ĥ) is constant (−26.26 raw, the ITT effect) for every rule, so a label tie can never move the target here.

Gate 2 output (verbatim; the checker was made NA-safe for the NO-DETECTION rows before this run):

```
== GATE 2: md=null n=500 ==
  [PASS] exists: fs_maxeffCons_mr_field_mdnull_knoise0_n500_mdf1_res_1_1000.rds
  [PASS] exists: fs_maxeffCons_mr_field_mdnull_knoise0_n500_mdf1_res_1001_2000.rds
  [PASS] exists: fs_maxeffCons_mr_field_mdnull_knoise0_n500_mdf1_combined_1_2000.rds
  [PASS] save guard: all three bundle paths untracked at save time (guard cannot have been bypassed)
  [PASS] completeness: sim_id 1-1000, 1001-2000, combined 1-2000
  [PASS] no CONFIG-ERROR rows (status table: DETECTED 1993, NO-DETECTION 7)
  pairing proof: 984 of 1000 sim_ids identical on 38 pre-existing numeric cols (max rel diff among identical rows 2.54e-10); rule-string-order-only differences on 12 sim_ids
  ENUMERATED FLIPS (16): 8, 9, 38, 95, 98, 226, 267, 331, 403, 412, 491, 549, 560, 615, 825, 868
    sim 8 -> label tie, no numeric consequence
    sim 9 -> label tie, no numeric consequence
    sim 38 -> label tie, no numeric consequence
    sim 95 -> same selection, MR numerics differ
    sim 98 -> label tie, no numeric consequence
    sim 226 -> label tie, no numeric consequence
    sim 267 -> same selection, MR numerics differ
    sim 331 -> label tie, no numeric consequence
    sim 403 -> label tie, no numeric consequence
    sim 412 -> label tie, no numeric consequence
    sim 491 -> label tie, no numeric consequence
    sim 549 -> label tie, no numeric consequence
    sim 560 -> label tie, no numeric consequence
    sim 615 -> same selection, MR numerics differ
    sim 825 -> label tie, no numeric consequence
    sim 868 -> label tie, no numeric consequence
  classification: 0 selection flips; 3 MR-numerics; 0 label ties with target move; 13 pure label ties
  [PASS] pairing proof: all 984 non-flip rows identical (<= 1e-8); 16 enumerated and excluded (0 selection flips)
  [PASS] new columns present
  finite share on 1993 detected reps: min 1.0000 (fld_H_est2); field notes set on 0 (H) / 0 (Hc)
  [PASS] fields finite on >= 99% of detected replicates
  [PASS] interval invariants lo <= hi
  [PASS] bound identities (max abs 0.0e+00)
  [PASS] gamma in [0.025, 0.05] (range 0.025-0.027; mean 0.0251)
  [PASS] p_hat(Hhat) finite in [0,1] (mean 0.155, share < 0.5: 0.987)
  meta: pkg 0.3.5 | workers 13/13 | ci field | complement TRUE | ij_residual two_term | built 2026-09-07 14:53:58
  timing: batch1 fit+MR mean 15.0 s (max 20.9), batch2 14.5 s; detection 0.997
GATE 2 md=null n=500: PASS
```

Enumerated rows: 13 pure label ties; sims 95, 267, 615 identical rule and naive estimate with the MR columns differing at ≤ 4.5e-5 relative. **No selection or detection flip.** Committed: the three bundles, `gate2_flips.txt`, and `fs_maxeffCons_mr_field_mdnull_knoise0_n500_mdf1_combine_1_2000.html`.

## Cell 4 — md40 n = 700 (priority 4; no FB)

Wall: batch 1 **1,391 s** (23.2 min; projection 24 min, stop threshold 36 min), batch 2 1,415 s, combine 22 s; cumulative **10,371 s = 172.8 min** of the 240 min ceiling. Per replicate: fit + MR + field 17.5 s mean (max 22.1). Detection 1,999 of 2,000.

Gate 2 output (verbatim):

```
== GATE 2: md=40 n=700 ==
  [PASS] exists: fs_maxeffCons_mr_field_md40_knoise0_n700_mdf1_res_1_1000.rds
  [PASS] exists: fs_maxeffCons_mr_field_md40_knoise0_n700_mdf1_res_1001_2000.rds
  [PASS] exists: fs_maxeffCons_mr_field_md40_knoise0_n700_mdf1_combined_1_2000.rds
  [PASS] save guard: all three bundle paths untracked at save time (guard cannot have been bypassed)
  [PASS] completeness: sim_id 1-1000, 1001-2000, combined 1-2000
  [PASS] no CONFIG-ERROR rows (status table: DETECTED 1999, NO-DETECTION 1)
  pairing proof: 983 of 1000 sim_ids identical on 38 pre-existing numeric cols (max rel diff among identical rows 1.56e-10); rule-string-order-only differences on 9 sim_ids
  ENUMERATED FLIPS (17): 82, 176, 184, 249, 303, 486, 499, 527, 571, 654, 663, 716, 758, 791, 795, 943, 956
    sim 82 -> same selection, MR numerics differ
    sim 176 -> same selection, MR numerics differ
    sim 184 -> same selection, MR numerics differ
    sim 249 -> same selection, MR numerics differ
    sim 303 -> label tie, no numeric consequence
    sim 486 -> label tie, no numeric consequence
    sim 499 -> label tie, no numeric consequence
    sim 527 -> label tie, super-population target moves
    sim 571 -> label tie, no numeric consequence
    sim 654 -> same selection, MR numerics differ
    sim 663 -> same selection, MR numerics differ
    sim 716 -> label tie, super-population target moves
    sim 758 -> label tie, super-population target moves
    sim 791 -> label tie, super-population target moves
    sim 795 -> label tie, no numeric consequence
    sim 943 -> label tie, no numeric consequence
    sim 956 -> label tie, super-population target moves
  classification: 0 selection flips; 6 MR-numerics; 5 label ties with target move; 6 pure label ties
  [PASS] pairing proof: all 983 non-flip rows identical (<= 1e-8); 17 enumerated and excluded (0 selection flips)
  [PASS] new columns present
  finite share on 1999 detected reps: min 1.0000 (fld_H_est2); field notes set on 0 (H) / 0 (Hc)
  [PASS] fields finite on >= 99% of detected replicates
  [PASS] interval invariants lo <= hi
  [PASS] bound identities (max abs 0.0e+00)
  [PASS] gamma in [0.025, 0.05] (range 0.025-0.027; mean 0.0251)
  [PASS] p_hat(Hhat) finite in [0,1] (mean 0.163, share < 0.5: 0.986)
  meta: pkg 0.3.5 | workers 13/13 | ci field | complement TRUE | ij_residual two_term | built 2026-09-07 15:41:22
  timing: batch1 fit+MR mean 17.5 s (max 22.1), batch2 17.8 s; detection 1.000
GATE 2 md=40 n=700: PASS
```

Enumerated rows: 6 pure label ties, 5 label ties with a super-population target move (the `{preanti <= 0}` / `!{str2}` pair), 6 with identical selection and MR columns differing at the 1e-5 to 1e-4 level. **No selection flip.** Committed: the three bundles, `gate2_flips.txt`, and `fs_maxeffCons_mr_field_md40_knoise0_n700_mdf1_combine_1_2000.html`.

## Campaign close

| cell | batch 1 wall | batch 2 wall | combine | per-replicate mean | Gate 2 | identical / enumerated (selection flips) |
|---|---|---|---|---|---|---|
| md40 n500 | 1,189 s | 1,209 s | 24 s | 14.9 s | PASS | 986 / 14 (0) |
| md120 n500 | 1,371 s | 1,344 s | 20 s | 17.3 s | PASS | 992 / 8 (0) |
| null n500 | 1,200 s | 1,163 s | 23 s | 15.0 s | PASS | 984 / 16 (0) |
| md40 n700 | 1,391 s | 1,415 s | 22 s | 17.5 s | PASS | 983 / 17 (0) |

8,000 replicates in 172.8 min cumulative wall on 13 workers; peak memory 17.5 GB (first render); no render exceeded 1.5× its projection; no cell deferred; every save went through `.refuse_if_tracked()` into untracked campaign paths. Across the four cells the pairing proof holds on 3,945 of 4,000 anchored replicates at ≤ 4.3e-10 relative; the 55 enumerated rows are 35 pure label ties, 10 label ties whose super-population target moves (`str2` vs `preanti > 0`, 7 of 1,083 patients), and 12 replicates with identical selection and MR columns differing at ≤ 2.4e-4 (0.2.2/0.3.1 → 0.3.5); **zero selection flips**.
