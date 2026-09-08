# REPORT — Studentized complement field, E0 (instrumented smoke on nb20-A HR 1.75 n500, 200 replicates; report-and-wait)

**Task:** `dev/tasks/TASK_field_studentize_stage1_e0_2026-09-08.md` (1bb8bdf1), E0. Proposal: `dev/tasks/PROPOSAL_complement_field_scale_2026-09-08_v2.md` §3, §6. Stage 1 record: `REPORT_field_studentize_stage1_2026-09-08.md` (7245e898; Gate 1 PASS, G1a and G1b). Decisions D-1–D-4 at defaults (200 replicates; the four recorder columns; nb20-A HR 1.75 n500; ceiling 30 min / timeout 1 h).
**Date:** 2026-09-08. Executor: Claude Code (Linux), unattended. Compute: this run only. Winner-only and winner-floor excluded from every table and line. **No R1 / R2 / R0 recommendation is made here — the record reports; the Linux chat decides.**

---

## GATE E0-a: PASS — every pre-existing column of the 200 e0stud rows `identical()` to the committed nb20-A HR 1.75 bundle rows sim_id 1–200 (131 / 131 non-timing columns; `truth` identical). The instrumentation is inert on everything old under production settings; the new columns are finite on all 200 detected rows.

## Run

`sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` (committed, 7245e898) driven by env only: `FS_S7_FOCUS=effMaxSG FS_S7_Z1Q=0.60 FS_S7_NBHD=0.20 FS_S7_N=500 FS_S7_HR=1.75 FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_FIELD_DECOMP=TRUE FS_S7_CAMPAIGN=e0stud FS_S7_NSIMS=200 FS_S7_START=1 FS_S7_WORKERS=100` (J = 10 default; `return_reselection = TRUE` in the template's `mr_inference_args`), seeds `8316951 + sim_id`, sim_id 1–200, one batch. Wall **240 s** (4.0 min; ceiling 30 min, timeout 1 h; mean fit+MR per replicate 73.5 s under the 100-worker load). Meta: `campaign_tag = e0stud`, `target_hr_harm = 1.75`, `n_sample = 500`, `effect_neighborhood = 0.20`, `er_jcuts = 10`, `harm_z1_quantile = 0.60`, `field_complement = TRUE`, `ij_residual = two_term`, `fb_mode = none`, `n_workers = 100`, forestsearch 0.3.5 (installed from 7245e898). Outputs beside the campaign's: `results/fs_effMaxSG_fb_mr_field_m1_h175_knoise0_n500_z1q60_nb20_e0stud_res_1_200.rds`, `..._e0stud_batch_1_200.html`; analysis script `e0stud_analysis.R` (this directory; every number below is its verbatim output). 200 / 200 detected; 200 analysed (every input finite). |H| 127–191, |Ĥ| 61–233, p̂(Ĥ) 0.002–0.753.

## Stage 0 quotes (from HEAD 1bb8bdf1 before the edit; full set with line numbers in the Stage 1 record)

```r
.fs_mr_field_complement <- function(df, spec, kept, Bc, bh_c, winset, sel, bdc,      # 989
                                    G_out, W_in, Xo, Xi_f, to_eff, z975,
                                    lam_H = NULL, beta_deb = NA_real_,
                                    alpha = 0.05) {
  ...
  fit_ok <- is.finite(bh_c)                                                          # 1016
  Zo_c <- crossprod(Bc, Xo)                  # Ncol x R_out : zeta^c (outer)         # 1017
  Zi_c <- crossprod(Bc, Xi_f)                # Ncol x R_in  : zeta'^c (inner)        # 1018
  ...
    lam_c[r]  <- Zo_c[G, r] - mean(Zi_c[cbind(wi[ok_in], ok_in)])                    # 1029
  ...
  sd_c <- stats::sd(lf)                                                              # 1045
  ...
    se_field = sd_c,                                                                 # 1062
```
Objects used by the instrumentation (Stage 1, add-only, after `lam_c` is final): `s[g] = sqrt(colSums(Bc * Bc))[g]`; `s_sel = s[sel]`; `s_win = s[G_out[ok_c]]`; ρᶜ = `s_sel / mean(s_win)`; `scale_win_cv = sd(s_win) / mean(s_win)`. Recorded as `fld_Hc_scale_sel`, `fld_Hc_scale_win`, `fld_Hc_scale_cv`, `fld_Hc_scale_ratio`.

**One identity to read the tables by.** `fld_Hc_scale_sel` equals `nv_Hc_se` on every replicate (corr 1.000, mean ratio 1.000): the selected complement's influence norm √Σ dfbeta² *is* the naive (robust) SE. Hence ρᶜ · `fld_Hc_se` / `nv_Hc_se` = `fld_Hc_se` / `fld_Hc_scale_win` exactly — "the corrected ratio ≈ 1" and "λ-SDᶜ ≈ the mean winner scale s̄_G" are the same statement.

## Table 1 — corr(ρᶜ, `fld_Hc_se`/`nv_Hc_se`) and mean ρᶜ by tertile (tertiles within the 200; n = 67 / 66 / 67)

| statistic | value | 95% CI / p |
|---|---|---|
| Pearson corr(ρᶜ, `fld_Hc_se`/`nv_Hc_se`) | −0.934 | [−0.950, −0.914] |
| Spearman corr(ρᶜ, `fld_Hc_se`/`nv_Hc_se`) | −0.930 | p ≈ 0 |
| Pearson / Spearman corr(ρᶜ, p̂(Ĥ)) | −0.522 / −0.634 | |
| Pearson / Spearman corr(ρᶜ, \|Ĥ\|/\|H\|) | 0.862 / 0.870 | |

(The deficit `fld_Hc_se`/`nv_Hc_se` falls as ρᶜ rises, so the sign is negative; the prediction "ρᶜ tracks the deficit" is the magnitude.)

**By p̂(Ĥ) tertile:**

| p̂ tertile | n | range | mean ρᶜ (SE) | median ρᶜ | share ρᶜ > 1 [Wilson] | mean `fld`/`nv` (SE) | mean ρᶜ·`fld`/`nv` (SE) | √(mean λ²)/√(mean nSE²) | corrected RMS form |
|---|---|---|---|---|---|---|---|---|---|
| T1 | 67 | [0.002, 0.048] | 1.112 (0.008) | 1.114 | 0.97 [0.90, 0.99] | 0.902 (0.008) | 0.999 (0.003) | 0.896 | 0.999 |
| T2 | 66 | [0.048, 0.133] | 1.047 (0.007) | 1.035 | 0.82 [0.71, 0.89] | 0.957 (0.007) | 0.999 (0.003) | 0.953 | 0.999 |
| T3 | 67 | [0.134, 0.753] | 1.009 (0.005) | 1.010 | 0.57 [0.45, 0.68] | 0.995 (0.005) | 1.003 (0.003) | 0.993 | 1.003 |

**By |Ĥ|/|H| tertile:**

| \|Ĥ\|/\|H\| tertile | n | range | mean ρᶜ (SE) | median ρᶜ | share ρᶜ > 1 [Wilson] | mean `fld`/`nv` (SE) | mean ρᶜ·`fld`/`nv` (SE) | √(mean λ²)/√(mean nSE²) | corrected RMS form |
|---|---|---|---|---|---|---|---|---|---|
| T1 | 67 | [0.350, 0.669] | 0.998 (0.004) | 0.999 | 0.48 [0.36, 0.60] | 1.006 (0.005) | 1.003 (0.003) | 1.005 | 1.003 |
| T2 | 66 | [0.669, 0.966] | 1.037 (0.004) | 1.038 | 0.88 [0.78, 0.94] | 0.967 (0.004) | 1.003 (0.003) | 0.966 | 1.003 |
| T3 | 67 | [0.968, 1.493] | 1.132 (0.006) | 1.125 | 1.00 [0.95, 1.00] | 0.881 (0.006) | 0.996 (0.003) | 0.880 | 0.996 |

**Against the prediction on record** (ρᶜ > 1 concentrated in low-p̂ / large-Ĥ, ≈ 1 at the other end): mean ρᶜ 1.11 at p̂ < 0.05 (97% of replicates above 1) falling to 1.01 at p̂ > 0.13 (57% above 1 — i.e. centred on 1); 1.13 in the top |Ĥ|/|H| tertile (every replicate above 1) against 1.00 in the bottom (48% above 1). Direction as predicted at both ends, with the top-|Ĥ| tertile the sharpest stratum.

## Table 2 — (ρᶜ · `fld_Hc_se`)/`nv_Hc_se` beside the uncorrected `fld_Hc_se`/`nv_Hc_se`

**By |Ĥ|/|H| tertile (A2's stratification; A2's committed pattern on the full cell: 1.00 / 0.92 / 0.79 of the naive SE²-scale ratio, i.e. ≈ 1.00 / 0.96 / 0.89 on the SE scale):**

| tertile | n | uncorrected mean `fld`/`nv` (SE) | corrected mean ρᶜ·`fld`/`nv` (SE) | uncorrected RMS | corrected RMS | mean `nv_Hc_se` | mean `fld_Hc_se` | mean ρᶜ·`fld_Hc_se` |
|---|---|---|---|---|---|---|---|---|
| T1 | 67 | 1.006 (0.005) | 1.003 (0.003) | 1.005 | 1.003 | 0.134 | 0.134 | 0.134 |
| T2 | 66 | 0.967 (0.004) | 1.003 (0.003) | 0.966 | 1.003 | 0.140 | 0.135 | 0.140 |
| T3 | 67 | 0.881 (0.006) | 0.996 (0.003) | 0.880 | 0.996 | 0.155 | 0.136 | 0.154 |
| all | 200 | 0.951 (0.005) | 1.000 (0.002) | 0.945 | 1.000 | 0.143 | 0.135 | 0.143 |

**By p̂(Ĥ) tertile:**

| tertile | n | uncorrected mean `fld`/`nv` (SE) | corrected mean ρᶜ·`fld`/`nv` (SE) | uncorrected RMS | corrected RMS | mean `nv_Hc_se` | mean `fld_Hc_se` | mean ρᶜ·`fld_Hc_se` |
|---|---|---|---|---|---|---|---|---|
| T1 | 67 | 0.902 (0.008) | 0.999 (0.003) | 0.896 | 0.999 | 0.150 | 0.134 | 0.149 |
| T2 | 66 | 0.957 (0.007) | 0.999 (0.003) | 0.953 | 0.999 | 0.142 | 0.135 | 0.141 |
| T3 | 67 | 0.995 (0.005) | 1.003 (0.003) | 0.993 | 1.003 | 0.137 | 0.136 | 0.138 |
| all | 200 | 0.951 (0.005) | 1.000 (0.002) | 0.945 | 1.000 | 0.143 | 0.135 | 0.143 |

Per-replicate: slope of `fld`/`nv` on |Ĥ|/|H| **−0.209** (SE 0.010); slope of ρᶜ·`fld`/`nv` on |Ĥ|/|H| **−0.010** (SE 0.007). SD across replicates: `fld_Hc_se` 0.006, ρᶜ·`fld_Hc_se` 0.011, `nv_Hc_se` 0.011; corr(`fld_Hc_se`, `nv_Hc_se`) 0.356 (A2 reported 0.34–0.47 on the full cells), corr(ρᶜ·`fld_Hc_se`, `nv_Hc_se`) **0.950**.

**Against the prediction on record** (the corrected ratio ≈ 1 in every tertile): 0.996–1.003 in all six strata under both stratifications, with SEs 0.003; the uncorrected ratio reproduces A2's tertile pattern on these 200 (1.006 / 0.967 / 0.881 ≈ the committed cell's 1.00 / 0.96 / 0.89). The rescaled field SD tracks the naive SE replicate by replicate (corr 0.95 against 0.36; the same SD across replicates, 0.011).

## Table 3 — `fld_Hc_scale_cv` (per-replicate CV of the used outer winners' complement scales s_G; the R1-vs-R2 per-draw stability input)

| mean | q10 | q25 | q50 | q75 | q90 | q99 | max |
|---|---|---|---|---|---|---|---|
| 0.054 | 0.037 | 0.045 | 0.054 | 0.064 | 0.072 | 0.087 | 0.089 |

| p̂ tertile | mean CV | q50 CV | q90 CV | mean `scale_sel` | mean `scale_win` |
|---|---|---|---|---|---|
| T1 | 0.058 | 0.056 | 0.081 | 0.150 | 0.134 |
| T2 | 0.056 | 0.055 | 0.072 | 0.142 | 0.135 |
| T3 | 0.049 | 0.047 | 0.065 | 0.137 | 0.136 |

Reading, not a recommendation: the winners' scales within a replicate vary by 4–9% (CV), against a selected-vs-mean-winner gap of 0–15% (ρᶜ); `scale_win` is flat across p̂ tertiles (0.134–0.136) while `scale_sel` moves (0.150 → 0.137) — the field's SD sits at the family-average scale in every stratum, as diagnosed.

## Table 4 — Context: ρᶜ against the cell's committed average deficit

| quantity | value |
|---|---|
| mean ρᶜ (SE) | 1.056 (0.005) |
| median ρᶜ | 1.038 |
| q10 / q90 ρᶜ | 0.978 / 1.146 |
| share ρᶜ > 1 [Wilson] | 0.79 [0.72, 0.84] (157 of 200) |
| these 200: λ² / naive SE² (mean of squares) | 0.893 |
| these 200: λ / naive SE | 0.945 |
| these 200: implied mean ρᶜ = naive SE / λ | 1.058 |
| these 200: corrected (ρᶜ·λ)² / naive SE² | 1.001 |
| these 200: mean per-replicate `fld`/`nv` | 0.951 |
| these 200: mean per-replicate ρᶜ·`fld`/`nv` | 1.000 |
| committed cell (A1, nb20-A HR 1.75 n500, n = 1999): λ² / naive SE² | 0.894 |
| committed cell: λ / naive SE = √0.894 | 0.945 |
| committed cell: implied mean ρᶜ = 1 / 0.945 | 1.058 |

The 200-replicate subset reproduces the committed cell's average deficit (0.893 vs 0.894), and the observed mean ρᶜ (1.056, SE 0.005) equals the deficit's implied value (1.058).

## What this record does and does not say

Reported: (i) ρᶜ is concentrated above 1 in the low-p̂ / large-Ĥ strata and centred on 1 at the other end; (ii) ρᶜ · λ-SDᶜ restores tracking of the naive SE in every tertile (0.996–1.003) and per replicate (corr 0.95); (iii) the per-draw winner-scale CV is 0.05 at the median, 0.07 at q90. Not said: which of R1 / R2 / R0 follows, any coverage number (E0 computes no bound with scaling; that is E1's shape test), and anything about the 2–3-point residual at the concentrated-pick ends (out of scope, handoff §5). Report and wait.
