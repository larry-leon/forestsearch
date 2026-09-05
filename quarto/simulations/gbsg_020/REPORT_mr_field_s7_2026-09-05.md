# REPORT — MR (field) on the Section 7 forestsearch cells: Stage 3 (final)

**Task:** `dev/tasks/TASK_mr_field_section7_2026-09-05.md` (bb516569). Records: Stage 0 ca2996b3 · Gate 1 STOP 8293a3ce · fix a6702fd8 · Gate 1 PASS 7f94fdc9 · Stage 2 STOP 9d39673d · Gate 2 amendment 98f608c5 · Gate 2 record beside this file.
**Date:** 2026-09-05. Data: campaign `s7` combined bundles (2 cells × 2,000 replicates, forestsearch 0.3.5, seeds `8316951 + sim_id`), FB joined from the committed `fb_mr` bundle on sim_id 1–500 (minus the enumerated flip sim 393), never recomputed. Rendered combined documents: `sim_fs_maxeffCons_fb_mr_field_m1_h100_knoise0_n500_s7_combine_1_2000.html`, `..._h175_..._s7_combine_1_2000.html`.

Conventions (the committed documents' own): estimates, bounds and lengths on the HR scale; SD/SE on the log-HR (β) scale; relative bias = 100·mean((est−target)/target); two-sided 95% coverage with Wilson intervals is the primary convention (F4), one-sided lower supplementary; MR (field) = `est2` with the two-sided Λ-quantile interval, Ĥ block only (F3 — Ĥᶜ stays on MR (IJ)); λ-SD (`se_field`) is the field's SE column.

## Estimation and coverage (Table S4 layout)

### Cell h100 (target_hr_harm = 1): 2000 replicates, 1361 detections (68.0%)

**Harm Ĥ block** (mean |Ĥ| = 79, 15.8% of n = 500; targets: θ† = 1.000, θ‡ = 1.000)

| Estimator | n | mean est | b_or % | b_β(Ĥ) % | b_θ‡ % | b_θ† % | SD_β | mean SE_β | SE/SD | Cov_or (Wilson) | Cov_β(Ĥ) (Wilson) | Cov_θ‡ (Wilson) | Cov_θ† (Wilson) | len mean | len med |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1361 | 1.629 | +56.5 | +130.7 | +62.9 | +62.8 | 0.163 | 0.324 | 1.99 | 0.791 (0.769, 0.812) | 0.157 (0.138, 0.177) | 0.869 (0.850, 0.886) | 0.869 (0.850, 0.886) | 2.25 | 2.17 |
| oracle | 1361 | 1.132 | +0.0 | +57.2 | +13.2 | +13.2 | 0.312 | 0.325 | 1.04 | 1.000 (0.997, 1.000) | 0.801 (0.779, 0.821) | 0.958 (0.946, 0.968) | 0.958 (0.946, 0.968) | 1.55 | 1.42 |
| FB | 328 | 1.223 | +16.5 | +72.0 | +22.3 | +22.2 | 0.186 | 0.356 | 1.91 | 0.973 (0.949, 0.985) | 0.814 (0.768, 0.852) | 0.994 (0.978, 0.998) | 0.994 (0.978, 0.998) | 1.89 | 1.80 |
| MR (IJ) | 1361 | 0.933 | -10.9 | +32.0 | -6.7 | -6.7 | 0.222 | 0.399 | 1.79 | 0.966 (0.955, 0.975) | 0.986 (0.978, 0.991) | 0.999 (0.995, 1.000) | 0.999 (0.995, 1.000) | 1.66 | 1.54 |
| MR (field) | 1361 | 0.848 | -19.2 | +19.8 | -15.2 | -15.2 | 0.253 | 0.357 | 1.41 | 0.813 (0.791, 0.832) | 0.988 (0.981, 0.993) | 0.940 (0.927, 0.952) | 0.940 (0.926, 0.951) | 1.19 | 1.07 |

One-sided lower 95% coverage (supplementary; field uses stored `lower_1s`, others exp(log est − 1.645·SE)):

| Estimator | vs oracle | vs β(Ĥ) | vs θ‡ | vs θ† |
|---|---|---|---|---|
| naive | 0.684 | 0.051 | 0.702 | 0.703 |
| oracle | 1.000 | 0.675 | 0.929 | 0.930 |
| FB | 0.945 | 0.655 | 0.976 | 0.976 |
| MR (IJ) | 0.993 | 0.966 | 0.999 | 0.999 |
| MR (field) | 0.993 | 0.963 | 0.997 | 0.997 |

**Complement Ĥᶜ block** (mean |Ĥ| = 421, 84.2% of n = 500; targets: θ† = 0.657, θ‡ = 0.585)

| Estimator | n | mean est | b_or % | b_β(Ĥ) % | b_θ‡ % | b_θ† % | SD_β | mean SE_β | SE/SD | Cov_or (Wilson) | Cov_β(Ĥ) (Wilson) | Cov_θ‡ (Wilson) | Cov_θ† (Wilson) | len mean | len med |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1361 | 0.609 | -7.3 | -7.9 | +4.1 | -7.3 | 0.131 | 0.139 | 1.06 | 0.997 (0.992, 0.999) | 0.911 (0.895, 0.925) | 0.960 (0.949, 0.969) | 0.919 (0.903, 0.933) | 0.34 | 0.33 |
| oracle | 1361 | 0.658 | +0.0 | -0.5 | +12.5 | +0.2 | 0.135 | 0.139 | 1.03 | 1.000 (0.997, 1.000) | 0.965 (0.954, 0.974) | 0.888 (0.870, 0.903) | 0.958 (0.946, 0.968) | 0.36 | 0.36 |
| FB | 328 | 0.632 | -3.7 | -4.4 | +8.0 | -3.8 | 0.135 | 0.193 | 1.43 | 1.000 (0.988, 1.000) | 0.985 (0.965, 0.993) | 0.982 (0.961, 0.992) | 0.982 (0.961, 0.992) | 0.49 | 0.48 |
| MR (IJ) | 1361 | 0.653 | -0.6 | -1.2 | +11.7 | -0.6 | 0.137 | 0.255 | 1.85 | 1.000 (0.997, 1.000) | 1.000 (0.997, 1.000) | 0.999 (0.995, 1.000) | 1.000 (0.997, 1.000) | 0.68 | 0.68 |

One-sided lower 95% coverage (supplementary; field uses stored `lower_1s`, others exp(log est − 1.645·SE)):

| Estimator | vs oracle | vs β(Ĥ) | vs θ‡ | vs θ† |
|---|---|---|---|---|
| naive | 1.000 | 0.996 | 0.936 | 0.993 |
| oracle | 1.000 | 0.977 | 0.807 | 0.961 |
| FB | 1.000 | 0.991 | 0.970 | 0.991 |
| MR (IJ) | 1.000 | 1.000 | 0.993 | 0.999 |

Retained bias vs β(Ĥ), harm block (log-HR / SD units): naive +0.8134 / +5.00; MR (IJ) β̃ +0.2442 / +1.10; MR (field) est2 +0.1405 / +0.56.
SE calibration, harm block: IJ SE/SD = 1.79; λ-SD/SD (field) = 1.41. Field draw usage: mean n_out 976/1000 (min 697), mean n_in 497/500; field secs mean 14.7.

### Cell h175 (target_hr_harm = 1.75): 2000 replicates, 1900 detections (95.0%)

**Harm Ĥ block** (mean |Ĥ| = 72, 14.4% of n = 500; targets: θ† = 1.769, θ‡ = 2.036)

| Estimator | n | mean est | b_or % | b_β(Ĥ) % | b_θ‡ % | b_θ† % | SD_β | mean SE_β | SE/SD | Cov_or (Wilson) | Cov_β(Ĥ) (Wilson) | Cov_θ‡ (Wilson) | Cov_θ† (Wilson) | len mean | len med |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1900 | 1.968 | +4.4 | +72.0 | -3.3 | +11.2 | 0.222 | 0.310 | 1.40 | 0.978 (0.971, 0.984) | 0.645 (0.623, 0.666) | 0.991 (0.985, 0.994) | 0.980 (0.973, 0.985) | 2.56 | 2.42 |
| oracle | 1900 | 1.970 | +0.0 | +68.3 | -3.2 | +11.4 | 0.297 | 0.311 | 1.05 | 1.000 (0.998, 1.000) | 0.682 (0.661, 0.703) | 0.961 (0.951, 0.968) | 0.960 (0.950, 0.968) | 2.59 | 2.36 |
| MR (IJ) | 1900 | 1.237 | -35.7 | +6.4 | -39.2 | -30.1 | 0.312 | 0.393 | 1.26 | 0.862 (0.845, 0.876) | 0.961 (0.951, 0.969) | 0.719 (0.699, 0.739) | 0.847 (0.830, 0.863) | 2.15 | 1.92 |
| MR (field) | 1900 | 1.180 | -39.3 | +0.5 | -42.0 | -33.3 | 0.358 | 0.363 | 1.01 | 0.610 (0.588, 0.632) | 0.887 (0.872, 0.901) | 0.458 (0.436, 0.481) | 0.610 (0.588, 0.632) | 1.73 | 1.47 |

One-sided lower 95% coverage (supplementary; field uses stored `lower_1s`, others exp(log est − 1.645·SE)):

| Estimator | vs oracle | vs β(Ĥ) | vs θ‡ | vs θ† |
|---|---|---|---|---|
| naive | 0.984 | 0.539 | 0.987 | 0.960 |
| oracle | 1.000 | 0.568 | 0.979 | 0.943 |
| MR (IJ) | 1.000 | 0.983 | 0.999 | 0.999 |
| MR (field) | 1.000 | 0.966 | 0.999 | 0.996 |

**Complement Ĥᶜ block** (mean |Ĥ| = 428, 85.6% of n = 500; targets: θ† = 0.657, θ‡ = 0.585)

| Estimator | n | mean est | b_or % | b_β(Ĥ) % | b_θ‡ % | b_θ† % | SD_β | mean SE_β | SE/SD | Cov_or (Wilson) | Cov_β(Ĥ) (Wilson) | Cov_θ‡ (Wilson) | Cov_θ† (Wilson) | len mean | len med |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1900 | 0.628 | -1.3 | -5.9 | +7.5 | -4.3 | 0.138 | 0.139 | 1.01 | 1.000 (0.998, 1.000) | 0.921 (0.907, 0.932) | 0.927 (0.914, 0.938) | 0.935 (0.923, 0.945) | 0.35 | 0.34 |
| oracle | 1900 | 0.638 | +0.0 | -4.6 | +9.1 | -2.9 | 0.140 | 0.139 | 0.99 | 1.000 (0.998, 1.000) | 0.939 (0.928, 0.949) | 0.913 (0.899, 0.925) | 0.943 (0.931, 0.952) | 0.35 | 0.35 |
| MR (IJ) | 1900 | 0.669 | +5.0 | +0.2 | +14.4 | +1.9 | 0.143 | 0.257 | 1.80 | 1.000 (0.998, 1.000) | 1.000 (0.998, 1.000) | 0.999 (0.996, 1.000) | 1.000 (0.998, 1.000) | 0.70 | 0.70 |

One-sided lower 95% coverage (supplementary; field uses stored `lower_1s`, others exp(log est − 1.645·SE)):

| Estimator | vs oracle | vs β(Ĥ) | vs θ‡ | vs θ† |
|---|---|---|---|---|
| naive | 1.000 | 0.985 | 0.883 | 0.978 |
| oracle | 1.000 | 0.981 | 0.858 | 0.972 |
| MR (IJ) | 1.000 | 0.999 | 0.983 | 0.999 |

Retained bias vs β(Ĥ), harm block (log-HR / SD units): naive +0.4809 / +2.17; MR (IJ) β̃ -0.0101 / -0.03; MR (field) est2 -0.0748 / -0.21.
SE calibration, harm block: IJ SE/SD = 1.26; λ-SD/SD (field) = 1.01. Field draw usage: mean n_out 993/1000 (min 766), mean n_in 499/500; field secs mean 16.3.

Paired vs FB (key cell, harm block, 328 joined detections): MR (IJ) β̃ − FB log-scale mean Δ -0.2793, SD Δ 0.0761, corr 0.9623.
Paired vs FB (key cell, harm block, 328 joined detections): MR (field) est2 − FB log-scale mean Δ -0.3811, SD Δ 0.1190, corr 0.9314.

## The two vintage flips (h100) and two cut-label findings (h175)

The Stage 2 STOP diagnosis (9d39673d) established, and the Gate 2 record confirms on the full runs: between the committed 0.2.0/0.2.1 bundles and the current 0.3.5 package, exactly **2 of 2,000** committed replicates change their *selection* — both in the borderline-null cell, both knife-edge near-ties at the consistency screen, with the **data identical and the oracle exact** on each:

- **sim 393:** `{er <= 121} & {pgr <= 0}` (n=62, naive HR 1.397) → `{pgr <= 0}` (n=64, naive 1.342; winner Pcons 0.92). The two-term rule's higher HR won under 0.2.0; under 0.3.5 it evidently falls just below the 0.90 consistency screen.
- **sim 780:** `{pgr <= 38} & {er <= 0}` → `{pgr <= 0} & {er <= 17}` (both n=62, naive 1.426 → 1.383).

Both are excluded from the FB join per the Gate 2 amendment (`fb_join_skip`, recorded in meta and printed with both selections in the render) and enter every other summary as ordinary fresh replicates. Additionally, two h175 replicates (**38**, **721**) carry a one-grid-step different *cut value* in the realized rule (`{age <= 53}` vs `{age <= 52}`; `{age <= 52}` vs `{age <= 51}`) with the identical trial subgroup — visible only in the β(Ĥ) truth attachment (rel diffs 0.041 / 0.012), every trial-level column exact. At 2 selection flips + 2 cut-label flips per 2,000 replicates, no cell summary moves beyond its Monte Carlo error.

## Reading criteria (Larry's, from the task — factual mapping, not gates)

- *"Conditional coverage within MC error of nominal or above (≥ 0.93 at ~1,300 detections), both cells and both blocks."* Ĥ block, Cov_β(Ĥ) two-sided: key cell — MR (field) **0.988** (0.981, 0.993), MR (IJ) 0.986 — met with room; h175 — MR (field) **0.887** (0.872, 0.901) — **below 0.93 beyond the Wilson interval** (MR (IJ): 0.961, met). Ĥᶜ block (on IJ per F3): 1.000 / 1.000 — met.
- *"SE/SD in [0.9, 1.2]."* MR (field) λ-SD/SD: h175 **1.01** — met; key cell **1.41** — above the band (conservative side; MR (IJ) is 1.79 there).
- *"est2's retained bias below β̃'s."* Key cell: **+0.141 vs +0.244 log-HR** (+0.56 vs +1.10 SD units) — met. h175: est2 **−0.075** vs β̃ −0.010 — |bias| larger with the sign flipped to overshoot: **not met**.
- *"Marginal coverage on the h175 cell not below MR (IJ)'s 0.85 by more than MC error."* MR (field) Cov_θ† = **0.610** (0.588, 0.632) vs MR (IJ) 0.847 — **not met, decisively**.

Summary of what the tables show (no interpretation beyond them): on the borderline-null key cell the field interval is the best-behaved of the five rows against the conditional target — coverage at/above nominal with intervals ~28% shorter than IJ (mean length 1.19 vs 1.66) and the smallest retained bias. On the genuine-harm h175 cell the field's tighter interval pays: conditional two-sided coverage drops below the criterion (0.887, one-sided still 0.966) and marginal coverage falls well below MR (IJ)'s, with est2 overshooting β(Ĥ) downward — the over-shrinkage risk this cell was chosen to expose. The complement block (IJ per F3) is unremarkable in both cells. FB↔MR agreement on the key cell: β̃ − FB mean Δ −0.279 (corr 0.962), est2 − FB mean Δ −0.381 (corr 0.931), log scale, 328 joined detections.

## Proposed follow-up (NOT started): vintage-flip attribution task

Per Larry's request, a proposal only. **Goal:** name the mechanism behind the 4 enumerated flips and decide whether a NEWS entry is warranted. **Sketch:** (1) `git worktree` at the 0.2.1 tree (the FB bundle's vintage), `R CMD INSTALL` into a temporary library (`R_LIBS_USER` override); (2) under each of 0.2.1 and 0.3.5, re-run h10 sims 393/780 and h175 sims 38/721 with the template's exact replicate settings under `RNGkind("L'Ecuyer-CMRG")` (the worker regime — a plain Mersenne-Twister session does NOT reproduce these replicates), instrumenting the search to dump the full candidate table (every candidate's effect, Pcons, size, and the cut grid for `age`/`er`/`pgr`) before the maxeffCons argmax; (3) diff the tables: the flip mechanism is then read off directly (a Pcons value crossing 0.90, a cut-grid change, or an RNG-stream realignment from the 0.3.4 effect-screen reordering), and the NEWS decision follows from whether the mechanism is behavioral (grid/screen semantics) or numerical (stream alignment). Estimated cost: ~1 h, no changes to `R/`, nothing re-run from the committed documents. Requires Larry's go.

## Done-means checklist (task §Done means)

- Stage 3 report (this file) and the two rendered combined documents: committed. ✔
- Gate 2 record beside them (`REPORT_mr_field_s7_gate2_2026-09-05.md`): committed. ✔
- Template committed under its own stem (7f94fdc9 + amendment 98f608c5); committed simulation documents untouched throughout. ✔
- Every Stage 1 identity listed with its measured value (`REPORT_mr_field_s7_stage1_2026-09-05.md`, 48/48). ✔
- Branch left unpushed for Larry. ✔
- F6 (template adoption as the standard document) remains Larry's call on this report.
