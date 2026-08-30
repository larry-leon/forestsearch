# REPORT — breadth, stage 1: the interaction ladder on the MD40 DGM, and a locked forecast

**Task:** `dev/tasks/cc_task_oc_breadth_ladder_2026-08-30.md` (arrived in `~/Downloads` as `cc_task_oc_breadth_ladder_20260830.md`, copied as `cc_task_oc_breadth_ladder_2026-08-30.md`, committed alone as `f5388486`).
**Script:** `dev/glm-continuous-sims/oc_breadth_ladder_2026-08-30.R` (parts `gate`, `rung <q>`, `forecast <q>`); logs `oc_breadth_ladder_2026-08-30_{gate,rung60,…,rung160,forecast120}.log`; outputs `oc_breadth_ladder_2026-08-30_{gate,rung60,…,rung160,forecast120}.rds` (each 130–350 KB; the per-threshold `results` lists are not stored).
**No `R/` change. No replicate was drawn. No push, no install, no render.**

## 0. Provenance and the build quotations

```
pop-os
/home/larryleon/Documents/GitHub/forestsearch
feature/glm-extension
f5388486          (after the task-doc commit; 2276ee10 was HEAD, clean tree, before it)
                  (git status --porcelain: empty)
f5388486 dev/tasks: add cc_task_oc_breadth_ladder_2026-08-30
2276ee10 report — record commit hash
fbd564de feat: consistency screen runs the plain loop for plan = "sequential" ... (0.3.1)
911d73ba task doc — parallel dispatch, narrowed: plan == "sequential" only
acb49b82 report — record commit hash
2a4c0126 residual quantiles — ...
packageVersion("forestsearch") = 0.3.1
```

**How the driver built MD40** (`quarto/simulations/actg175/continuous/sim_fs_maxeffCons_mr_md40_knoise0_n500_batch_1_1000.qmd`, chunks `setup` L123–269 and `build-dgm` L282–336; the n = 700 driver differs from it in one line only, `n_sample <- 700L` (L133, verified by `diff`); the null cell is the same document with `null_cell = TRUE`, which routes L321–326):

```r
seed_base <- 8316951L; target_md_harm <- -40; adverse_outcome <- FALSE
actg_arms <- c(1L, 3L); actg_treat_arm <- 1L; actg_age_cut <- 34; actg_preanti_cut <- 744.5
dgm_n_super <- 5000L; cal_k_grid_range <- c(0, 120); cal_grid_step <- 2
dgm_factor_vars <- paste0("z", 1:12); dgm_subgroup_vars <- c("z1","z2"); dgm_subgroup_cuts <- list(z1 = 1L, z2 = 1L)
actg_df <- subset(speff2trial::ACTG175, arms %in% actg_arms)
actg_df$id <- seq_len(nrow(actg_df))
actg_df$treat_orig <- ifelse(actg_df$arms == actg_treat_arm, 1L, 0L)
actg_df$treat      <- 1L - actg_df$treat_orig               # ddI = 1
actg_df$cd4_change <- actg_df$cd420 - actg_df$cd40
actg_df <- actg_df[!is.na(actg_df$cd420), ]
actg_df$z1 <- as.factor(ifelse(actg_df$age > actg_age_cut, 1L, 0L))
actg_df$z2 <- as.factor(ifelse(actg_df$preanti <= actg_preanti_cut, 1L, 0L))
actg_df$z3 <- as.factor(ifelse(actg_df$wtkg <= 75, 1L, 0L))
actg_df$z4..z6 <- karnof/cd40/cd80 <= median ; z7..z12 <- as.factor(hemo, homo, drugs, race, gender, symptom)
for (v in analysis_binary_vars) actg_df[[v]] <- as.factor(actg_df[[v]])   # hemo homo drugs race gender symptom str2
dgm <- calibrate_glm_interaction(data = actg_df, factor_vars = dgm_factor_vars,
    outcome_var = "cd4_change", treatment_var = "treat", target_effect = target_md_harm,
    outcome_type = "continuous", effect_measure = "MD", subgroup_vars = dgm_subgroup_vars,
    subgroup_cuts = dgm_subgroup_cuts, k_inter_range = cal_k_grid_range, grid_step = cal_grid_step,
    n_super = dgm_n_super, seed = seed_base, verbose = FALSE)
# null cell: generate_glm_dgm(<same>, model = "null", n_super = dgm_n_super, seed = seed_base)
```
`continuous_vars` is not passed (default `NULL`); `k_treat` is not passed (default 1); `adverse_outcome` is not passed to the calibrator (default `FALSE`). The three drivers agree on the build.

**The `forestsearch_args` and settings the corrected run used** (`oc_wrapper_grid_corrected_2026-08-30.R` and the `.rds`'s `fs_args` / `settings`): `confounders.name` = the driver's 13 (`age preanti wtkg karnof cd40 cd80 hemo homo drugs race gender symptom str2`), `conf.cont_jcuts = list(age = 10, preanti = 10)`, `n.min = 60`, `maxk = 2`, `sg_focus = "maxeffCons"`, `effect.threshold c1 = 30`, `consistency.threshold c2 = 10`, `pconsistency.threshold = 0.9`; `draws = 2e5`, `block = Inf`, `seed = 20260825`. Used unchanged for every rung. One difference of the corrected script from the driver is on the record: its `actg_frame()` does not convert the seven analysis binaries to factor and carries no `treat_orig`; it built the DGM directly with `generate_glm_dgm(..., k_inter = truth$beta_inter)`. §2 shows this makes no difference to anything the wrapper consumes.

## 1. THE FORECAST — locked before any replicate is drawn

**Forecast rung: q = 120 (`|m_tau[Q]| = 120`, `k_inter = −93.7447641240`, `bint = 93.7448`), n = 500, driver `c2 = 10`, `pconsistency = 0.9`.**
Chosen by §5's rule: power at `c1_05` is 0.6125 at q = 100 and 0.8372 at q = 120; linear interpolation gives q = 100 + 20·(0.80 − 0.6125)/(0.8372 − 0.6125) = 116.69, rounded up to the nearest multiple of 5 = **120**; the rung was rebuilt under §3's gates in the `forecast` part (all TRUE, §5 below).

**`c1*` = 135.7412** — the largest `c1` with predicted power ≥ 0.80 under `"resample"` (`fs_oc_invert`, target 0.80, the exact order statistic on the 2e5-draw set): achieved 0.80000 (MC SE 0.00089), ceiling 1.00000, attainable. Cross-check: the blocked grid's `det_rate` at the nearest integer `c1 = 136` is 0.79630 (at 135: above 0.80, monotone). Under `"split"` the same inversion gives 135.68.

| quantity at (n = 500, c1* = 135.741, c2 = 10) | `"resample"` | `"split"` |
|---|---:|---:|
| **predicted power** (declaration rate) ± MC SE | **0.8000 ± 0.0009** | 0.7992 ± 0.0009 |
| **predicted null type-I** at c1* ± MC SE (stored null family, resample) | **0.03868 ± 0.00043** | — (not run; §4 note) |
| `EnH` (E\|Ĥ\| given declaration) | 72.84 | 72.81 |
| `Esens` | 0.3326 | 0.3327 |
| `Espec` | 0.9526 | 0.9527 |
| `Eppv` | 0.7752 | 0.7755 |
| `Enpv` | 0.7320 | 0.7320 |
| `EbetaH` (oriented E[β(Ĥ)] given declaration; truth in Q is 120) | 98.93 | 98.96 |
| `Enaive_bias` | 62.26 | 62.24 |
| `mass_below` | 1.0000 | 1.0000 |
| `sel_c` mass on Q itself | 0 (Q is not a candidate; see below) | 0 |
| `sel_c` mass on candidates with `PQg ≥ 0.95` (15 candidates) | 0.3448 | 0.3449 |

Only `det_rate` carries an MC SE in the `fs_oc_predict` object (`det_rate_se`; `P1_se`, `p_sel_se` are per-candidate); the other columns are draw-set expectations without a stored SE.

**The null side at c1*:** `fs_oc_predict(family = <re-enumerated null family, gated identical to the stored one on lab/Pg/PQg/beta_g/se_g/sens_g/spec_g/M>, n = 500, c1 = 135.7412, c2 = 10, "resample", pconsistency 0.9, draws 2e5, seed 20260825)` → false declaration **0.03868 (SE 0.00043)**, `EnH` 66.55, `EbetaH` 26.26 (= the common effect), naive bias 118.41. It is below 0.05 because c1* = 135.74 > c1_05 = 133.235.

**Selection given declaration (resample; split is the same to 3 decimals).** Q itself — `age > 34 & preanti ≤ 744.5` — is **not a member of the family**: the population-quantile cuts on `age` are {25, 28, 30, 31, 33, 35, 37, 39, 42, 47} and on `preanti` {0, 39, 188, 402, 672, 875, 1055}, so no candidate reproduces Q's boundaries; the 15 candidates with `PQg = 1` are strict subsets of Q, the largest with `sens_g = 0.860`. The top 15 by `sel_c`:

| rank | candidate | `sel_c` | `Pg` | `PQg` | `beta_g` |
|---:|---|---:|---:|---:|---:|
| 1 | age > 39 & preanti ≤ 188 | 0.0795 | 0.1252 | 1.0000 | 120.00 |
| 2 | age > 42 & preanti ≤ 875 | 0.0596 | 0.1268 | 0.9259 | 113.05 |
| 3 | age > 39 & preanti ≤ 402 | 0.0382 | 0.1498 | 1.0000 | 120.00 |
| 4 | age > 37 & preanti ≤ 0 | 0.0363 | 0.1384 | 1.0000 | 120.00 |
| 5 | age > 37 & preanti ≤ 39 | 0.0320 | 0.1458 | 1.0000 | 120.00 |
| 6 | age > 35 & homo != 1 | 0.0287 | 0.1260 | 0.7381 | 95.45 |
| 7 | age > 37 & str2 != 1 | 0.0277 | 0.1390 | 1.0000 | 120.00 |
| 8 | age > 33 & cd80 ≤ 645 | 0.0270 | 0.1236 | 0.7152 | 93.30 |
| 9 | age > 35 & cd40 ≤ 260 | 0.0269 | 0.1236 | 0.7233 | 94.06 |
| 10 | age > 35 & wtkg > 83 | 0.0260 | 0.1202 | 0.7238 | 94.11 |
| 11 | age > 39 & preanti ≤ 672 | 0.0233 | 0.1750 | 1.0000 | 120.00 |
| 12 | age > 33 & cd40 > 420 | 0.0210 | 0.1224 | 0.6912 | 91.05 |
| 13 | age > 35 & preanti ≤ 0 | 0.0208 | 0.1696 | 1.0000 | 120.00 |
| 14 | age > 35 & age ≤ 39 | 0.0194 | 0.1708 | 0.7436 | 95.96 |
| 15 | age > 33 & race | 0.0187 | 0.1328 | 0.6913 | 91.06 |

**The same rows at the driver's own c1 = 30, this harm** (the saturated screen, for the record beside the crossover):

| (n = 500, c1 = 30, c2 = 10) | det_rate ± SE | EnH | Esens | Espec | Eppv | Enpv | EbetaH | naive | mass_below |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| resample | 1.0000 ± 0.0000 | 73.58 | 0.3277 | 0.9478 | 0.7556 | 0.7296 | 97.09 | 56.92 | 0.0040 |
| split | 1.0000 ± 0.0000 | 73.56 | 0.3277 | 0.9478 | 0.7558 | 0.7296 | 97.11 | 56.86 | 0.0038 |

**The three predictions Stage 2 will be scored against.**
(a) At `c1 = 135.741` (c2 = 10, pconsistency 0.9, n = 500, resample) the declaration rate is **0.800** (MC SE 0.0009 on the prediction; a 1000-replicate cell has binomial SE ≈ 0.013).
(b) The null false-declaration rate at that same `c1*` is **0.0387** (SE 0.0004); Stage 2 checks this against the **existing** null payload (`fs_maxeffCons_mr_mdnull_knoise0_n500_res_1_1000.rds`) by re-thresholding, with no new null run.
(c) Given declaration, the selected rule has `E[β(Ĥ)]` = **98.9** (oriented; truth in Q = 120), specificity **0.953**, PPV **0.775**, sensitivity **0.333**, `E|Ĥ|` = **72.8** — **carrying the documented limitation**: `E|Ĥ|`, PPV, sensitivity and the naive bias inherit a between-rule gap of +2.11 / +0.61 / +1.65 subjects at MD40 that five mechanisms failed to explain (`REPORT_residual_2026-08-30.md`, closed as a documented limitation), and the population-versus-sample offset on PPV and sensitivity is expected; whether that gap moves with the harm is a finding of Stage 2, not a prediction of Stage 1.

## 2. The ladder (§3) and power at fixed type-I (§4)

Rungs are `q = |m_tau[Q]|`; `s = sign(m_tau[Q]) = −1` (MD40 scale table: `m_tau` Q = −40.000000, Qc = −26.2552358760, S = −30.9916815932); `k_inter(q) = k40 − (q − 40)` with `k40 = −13.7447641240`. Every rung: `|m_tau[Q]| = q` exactly (printed to 12 decimals), `|m_tau[Qc]| = 26.255235876036` unchanged, family `lab, Pg, PQg, se_g, sens_g, spec_g, M` `identical()` to MD40's, `max|beta_g − beta_g_MD40 − (q − 40)·PQg|` ≤ 2.8e-14. `bint = q − 26.2552`.

**At the driver's own c1 = 30, c2 = 10 (resample; the q = 40 row is the §2 reproduction, identical to the stored object):**

| q | k_inter | bint | det_rate ± SE | EnH | Esens | Espec | Eppv | Enpv | EbetaH | naive bias | mass_below |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 40 | −13.7448 | 13.7448 | 0.9991 ± 0.0001 | 70.79 | 0.1638 | 0.8701 | 0.4000 | 0.6643 | 31.75 | 75.53 | 0.2647 |
| 60 | −33.7448 | 33.7448 | 0.9999 ± 0.0000 | 70.78 | 0.2029 | 0.8907 | 0.4938 | 0.6802 | 42.92 | 73.27 | 0.0907 |
| 80 | −53.7448 | 53.7448 | 1.0000 ± 0.0000 | 71.12 | 0.2441 | 0.9113 | 0.5879 | 0.6968 | 57.85 | 69.15 | 0.0397 |
| 100 | −73.7448 | 73.7448 | 1.0000 ± 0.0000 | 72.04 | 0.2863 | 0.9307 | 0.6768 | 0.7134 | 76.17 | 63.49 | 0.0130 |
| 120 | −93.7448 | 93.7448 | 1.0000 ± 0.0000 | 73.58 | 0.3277 | 0.9478 | 0.7556 | 0.7296 | 97.09 | 56.92 | 0.0040 |
| 140 | −113.7448 | 113.7448 | 1.0000 ± 0.0000 | 75.65 | 0.3668 | 0.9620 | 0.8215 | 0.7446 | 119.70 | 50.10 | 0.0010 |
| 160 | −133.7448 | 133.7448 | 1.0000 ± 0.0000 | 78.08 | 0.4021 | 0.9731 | 0.8736 | 0.7580 | 143.10 | 43.68 | 0.0002 |

**Power at the null's thresholds — the central column.** From the corrected run's stored null inversions (n = 500, resample): `c1_05 = 133.235` (target 0.05, achieved 0.05000, ceiling 0.99694) and `c1_10 = 125.700` (target 0.10; **stored**, so no new inversion was run). Read at the nearest integer grid point: 133.235 → 133, 125.700 → 126.

| q | power at c1_05 (grid 133) ± SE | power at c1_10 (grid 126) ± SE | EnH @133 | Eppv @133 | Esens @133 | EbetaH @133 |
|---:|---:|---:|---:|---:|---:|---:|
| 60 | 0.1768 ± 0.0009 | 0.2833 ± 0.0010 | 67.68 | 0.5507 | 0.2165 | 44.84 |
| 80 | 0.3584 ± 0.0011 | 0.5011 ± 0.0011 | 68.82 | 0.6402 | 0.2573 | 60.66 |
| 100 | 0.6125 ± 0.0011 | 0.7415 ± 0.0010 | 70.73 | 0.7118 | 0.2956 | 78.75 |
| **120** | **0.8372 ± 0.0008** | 0.9118 ± 0.0006 | 72.96 | 0.7720 | 0.3319 | 98.63 |
| 140 | 0.9580 ± 0.0004 | 0.9820 ± 0.0003 | 75.44 | 0.8264 | 0.3679 | 120.26 |
| 160 | 0.9938 ± 0.0002 | 0.9979 ± 0.0001 | 78.04 | 0.8744 | 0.4022 | 143.21 |

**The inverted c1 per rung (`fs_oc_invert`, resample, seed 20260825, one-block draw set) with the blocked grid's `det_rate` at the nearest integer as cross-check** (agreement to MC precision, not bit-for-bit, as expected):

| q | c1 @ 80% | grid det @ nearest | c1 @ 90% | grid det | c1 @ 95% | grid det | ceiling |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 60 | 100.541 | 0.7922 @101 | 93.266 | 0.9029 @93 | 87.368 | 0.9524 @87 | 0.99987 |
| 80 | 110.642 | 0.7940 @111 | 103.054 | 0.9006 @103 | 96.949 | 0.9497 @97 | 0.99999 |
| 100 | 122.340 | 0.8050 @122 | 114.383 | 0.9037 @114 | 108.065 | 0.9504 @108 | 1.00000 |
| 120 | 135.741 | 0.7963 @136 | 127.304 | 0.9029 @127 | 120.622 | 0.9477 @121 | 1.00000 |
| 140 | 150.589 | 0.7942 @151 | 141.682 | 0.8973 @142 | 134.589 | 0.9476 @135 | 1.00000 |
| 160 | 166.807 | 0.7975 @167 | 157.419 | 0.9034 @157 | 149.938 | 0.9497 @150 | 1.00000 |

Every inversion's achieved rate is exactly the target (0.80000 / 0.90000 / 0.95000; SE 0.00089 / 0.00067 / 0.00049).

**The power-against-c1 figure (resample `det_rate`, per harm; type-I from the null: 0.05 at 133.2, 0.10 at 125.7):**

| c1 | q60 | q80 | q100 | q120 | q140 | q160 |
|---:|---:|---:|---:|---:|---:|---:|
| 0–60 | ≥0.9996 | 1.0000 | 1.0000 | 1.0000 | 1.0000 | 1.0000 |
| 80 | 0.9824 | 0.9965 | 0.9996 | 1.0000 | 1.0000 | 1.0000 |
| 100 | 0.8091 | 0.9284 | 0.9830 | 0.9977 | 0.9998 | 1.0000 |
| 120 | 0.3990 | 0.6271 | 0.8339 | 0.9535 | 0.9925 | 0.9993 |
| 140 | 0.1024 | 0.2391 | 0.4722 | 0.7337 | 0.9147 | 0.9842 |
| 160 | 0.0146 | 0.0519 | 0.1579 | 0.3708 | 0.6539 | 0.8770 |
| 180 | 0.0013 | 0.0069 | 0.0319 | 0.1168 | 0.3113 | 0.5979 |
| 200 | 0.0001 | 0.0006 | 0.0044 | 0.0233 | 0.0955 | 0.2766 |
| 220 | 0.0000 | 0.0000 | 0.0004 | 0.0032 | 0.0193 | 0.0843 |
| 240–300 | 0 | 0 | 0 | ≤0.0003 | ≤0.0027 | ≤0.0174 |

The full 0..300 tables are in each rung's `.rds` (`$grid$table`). The rung-120 resample grid table and the forecast part's resample grid table are `identical()` (same call, same seed).

## 3. The `se_g` diagnostic per rung (§4) — reported, not acted on

`ratio = se_g / se_direct` with `se_direct = sqrt(V_eff[g] / (n·Pg))` from `fs_dgm_scale(dgm, regions = <each family$memb column>)` (P_g gate `all.equal` to `fam$Pg` passed on every rung), n = 500, over all M = 1696 candidates. At MD40 the band was [0.992, 1.015] (`REPORT_residual_2026-08-30.md` §2).

| q | min | 5% | median | 95% | max | cor(ratio, Pg) | outside [0.992, 1.015] | outside ±2% |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 40 | 0.99226 | 0.99543 | 0.99870 | 1.00677 | 1.01062 | −0.232 | 0 | 0 |
| 60 | 0.99170 | 0.99427 | 0.99737 | 1.00503 | 1.00935 | −0.329 | 4 | 0 |
| 80 | 0.98636 | 0.98977 | 0.99443 | 1.00152 | 1.00768 | −0.439 | 473 | 0 |
| 100 | 0.97833 | 0.98177 | 0.98872 | 0.99935 | 1.00717 | −0.453 | 1184 | 11 |
| 120 | 0.96745 | 0.97093 | 0.98034 | 0.99887 | 1.00643 | −0.417 | 1442 | 813 |
| 140 | 0.95394 | 0.95778 | 0.96951 | 0.99876 | 1.00600 | −0.382 | 1486 | 1363 |
| 160 | 0.93853 | 0.94286 | 0.95647 | 0.99872 | 1.00600 | −0.358 | 1505 | 1425 |

By purity band (median [min, max]; n = 436 / 1245 / 15):

| q | PQg < 0.25 | 0.25 ≤ PQg < 0.95 | PQg ≥ 0.95 |
|---:|---|---|---|
| 40 | 0.9977 [0.9923, 1.0080] | 0.9990 [0.9930, 1.0106] | 1.0002 [0.9995, 1.0007] |
| 60 | 0.9972 [0.9922, 1.0080] | 0.9974 [0.9917, 1.0094] | 1.0002 [0.9995, 1.0007] |
| 80 | 0.9965 [0.9899, 1.0077] | 0.9934 [0.9864, 1.0058] | 1.0002 [0.9995, 1.0007] |
| 100 | 0.9950 [0.9849, 1.0072] | 0.9868 [0.9783, 1.0005] | 1.0002 [0.9995, 1.0007] |
| 120 | 0.9926 [0.9772, 1.0064] | 0.9777 [0.9675, 0.9961] | 1.0002 [0.9995, 1.0007] |
| 140 | 0.9894 [0.9674, 1.0060] | 0.9664 [0.9539, 0.9930] | 1.0002 [0.9995, 1.0007] |
| 160 | 0.9859 [0.9560, 1.0060] | 0.9530 [0.9385, 0.9891] | 1.0002 [0.9995, 1.0007] |

**Plain statement.** The band stays inside MD40's [0.992, 1.015] at q = 40 and (to 4 candidates, min 0.9917) at q = 60; it **leaves that band at q = 80** (473 candidates, min 0.986) and leaves ±2% at q = 100 (11 candidates) and decisively at q = 120 (813 candidates, min 0.967, median 0.980). The drift is one-sided (`se_g` too small: the prevalence scaling omits the `½·k_inter²·PQg(1−PQg)` bracket gain) and lives entirely in the mixed-purity band; the pure-Q band (PQg ≥ 0.95) is unchanged at every rung, as the algebra says it must be. At the forecast rung the median understatement of `se_g` is 2.0% overall and 2.2% in the mixed band. Reported, not acted on; used only in §4 below as a sensitivity.

## 4. The sensitivity (§6) beside the forecast — forecast rung only

`fs_oc_predict(family = <rung-120 family with se_g replaced by se_direct>, n = 500, c1 = c1* = 135.741, c2 = 10, "resample", pconsistency 0.9, draws 2e5, seed 20260825)`:

| at c1* (resample) | forecast (wrapper as it is) | se_direct family |
|---|---:|---:|
| power | 0.8000 ± 0.0009 | 0.8101 ± 0.0009 |
| EnH | 72.84 | 72.63 |
| Eppv | 0.7752 | 0.7633 |
| Esens | 0.3326 | 0.3265 |
| EbetaH | 98.93 | 97.81 |
| sel_c mass on Q (Q not a candidate) | 0 | 0 |
| sel_c mass on PQg ≥ 0.95 | 0.3448 | 0.3246 |

This is a sensitivity showing what the prevalence-scaling approximation costs at this harm. It is **not** a proposed constructor, it is **not** adopted, and the forecast of §5 stands as issued with the wrapper as it is. Whether `se_g` should ever come from the direct per-candidate `V_eff` is Larry's decision, taken on these numbers.

## 5. Every gate's arithmetic

**§1.** hostname pop-os; branch `feature/glm-extension`; tree clean; `2276ee10` in `git log -6`; version 0.3.1. PASS.

**§2 — reproduce MD40.**
- `calibrate_glm_interaction(<driver call>)` → `model_params$k_inter = −13.7447641240 = model_params$beta_inter`; payload `truth$beta_inter = −13.7447641240`; |difference| = 0.000e+00 < 1e-9. PASS.
- `fs_dgm_scale`: m_tau Q = −40.0000000000, Qc = −26.2552358760, S = −30.9916815932 (P_g 0.3446 / 0.6554 / 1; V_eff 67453.760 / 68034.437 / 67891.079).
- `identical(df_super[calibrate], df_super[generate_glm_dgm(k_inter = k40)])` = TRUE (every column, `mu1` included). PASS.
- Family at n = 500, `fs_oc_family_enumerate(dgm, fs_args, n = 500, max_M = 5000)`: cut columns 74, enumerated 2775, dropped empty 127 / rmin 107 / size 656, kept 1885, duplicates 189, **M = 1696**; `identical()` to the stored corrected family on `lab, Pg, PQg, beta_g, se_g, sens_g, spec_g, M`: all TRUE. **`ovl` and `memb` are not in the stored record** (`family_record()` in the corrected script keeps neither), so they could not be gated against the store; they are deterministic functions of `memb`, which is a deterministic function of `df_super` (identical) and the cuts (identical `lab`). Stated, not worked around.
- Draw-level: the stored n = 500 alternative point was `fs_oc_predict()` (one block; the corrected script's `runs`), not a grid cell. Re-run as `fs_oc_predict(family = fam, n = 500, c1 = 30, c2 = 10, consistency_method = "resample", pconsistency = 0.9, draws = 2e5, seed = 20260825)`: `identical()` to the stored object = **TRUE** (366 s; det_rate 0.99906). A first attempt with `n = 500L` was identical on every computed field and differed only in `settings$n` (`int` vs the stored `num`, the corrected script having called `run_alt(500)`); the call form was corrected to `n = 500` and re-run — the failing log line is preserved in this report for the record, the committed `gate.log` is the passing run. PASS. (This is the standing guard that 0.3.1 changed nothing the wrapper touches.)

**§3 — per rung (60, 80, 100, 120, 140, 160, and the forecast rebuild of 120):** `||m_tau[Q]| − q|` = 0 to 12 printed decimals; `|m_tau[Qc]|` = 26.255235876036 at every rung = MD40's; family `lab, Pg, PQg, se_g, sens_g, spec_g, M` identical to MD40's; `max|beta_g − beta_g_MD40 − (q − 40)·PQg|` = 1.07e-14 (60), 1.42e-14 (80, 100, 120, 140), 2.84e-14 (160). All PASS.

**§4 — null.** `c1_05 = 133.235` and `c1_10 = 125.700` both read from the stored null inversion table (resample; ceiling 0.99694); no new inversion, nothing else re-derived. For the §5 null rate at c1* the null family was re-enumerated (`generate_glm_dgm(model = "null", <driver args>)`) and gated `identical()` to the stored null family on `lab, Pg, PQg, beta_g, se_g, sens_g, spec_g, M`: all TRUE. PASS.

**§5 — forecast rung.** §3's gates at the rebuilt q = 120: all TRUE. `c1*` = 135.7412, attainable TRUE (ceiling 1.00000). PASS.

**Compute.** gate 378 s; rungs 2237–2255 s each (six concurrent processes: draws ≈ 330 s, 301-point sweep ≈ 1220 s, inversion ≈ 680 s); forecast 5338 s (two grids 1155 + 1444 s, two inversions, five one-block predicts). Wall-clock ≈ 2 h 15 min.

## 6. Ten-line summary

1. Task committed alone as `f5388486` (arrived hyphen-stripped as `cc_task_oc_breadth_ladder_20260830.md`); §1 provenance gate passed on 0.3.1.
2. MD40 rebuilt bit-for-bit as the driver did (calibrator, k40 = −13.7447641240 = payload truth), direct route identical on `df_super`, family identical (M = 1696), `fs_oc_predict` identical to the stored object — 0.3.1 touched nothing the wrapper calls.
3. Six rungs q = 60…160 built by `k_inter = k40 − (q − 40)`; every DGM and family gate exact (beta shift to 1e-14), as source said.
4. At the driver's c1 = 30 the screen saturates at every harm (det ≥ 0.9991): no `c1 = 30` cell separates power from type-I anywhere on the ladder.
5. At the null's 5% threshold (c1_05 = 133.235) power runs 0.18 / 0.36 / 0.61 / 0.84 / 0.96 / 0.99 — the crossover regime the wrapper predicts.
6. **Forecast rung q = 120** (interpolated 116.7, rounded up to 120); **c1\* = 135.741**; predicted power **0.800 ± 0.001**, predicted null type-I at c1\* **0.0387 ± 0.0004**.
7. At c1\*: E|Ĥ| 72.8, sensitivity 0.333, specificity 0.953, PPV 0.775, E[β(Ĥ)] 98.9 (truth 120), naive bias 62.3; split gate within 0.001 of resample on every column.
8. Q itself is not a family member (cuts at age 33/35, preanti 672/875); `sel_c` mass on PQg ≥ 0.95 candidates is 0.345, led by `age > 39 & preanti ≤ 188`.
9. `se_g`/`se_direct` leaves MD40's 2% band from q = 80 and reaches a 2.0% median understatement at q = 120 (mixed-purity candidates only); the se_direct sensitivity moves power at c1\* from 0.800 to 0.810 and PPV from 0.775 to 0.763 — reported, not adopted.
10. Lock: the report was committed as `c9cb0ca2`; this line was added in the immediately following commit. Stage 2 is a separate go/no-go.
