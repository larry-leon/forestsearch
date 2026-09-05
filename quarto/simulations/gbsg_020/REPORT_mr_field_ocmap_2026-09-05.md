# REPORT — OC map: MR (IJ) and MR (field) on cells not yet covered (Stage 3, final)

**Task:** `dev/tasks/TASK_mr_field_ocmap_2026-09-05.md` (0ef90878). Records: Stage 0 26718386 · Stage 1 1533b99d · Gate 2 record beside this file. Unattended run, H-Q1–H-Q6 defaults.
**Data:** campaign `map1`, five cells × 2,000 replicates (seeds `8316951 + sim_id`, two seed-disjoint batches each), committed template at b1316d22 driven by `FS_S7_*` env only, `ci_method = "field"` (both MR rows + the one-sided bound from one gate call), FB none, uniform-κ knob off. No `R/`, template, or committed-document changes. Rendered combined documents: `fs_maxeffCons_fb_mr_field_m1_{h150_knoise0_n500,h075_knoise0_n500,h100_knoise0_n1000,h175_knoise3_n500,h150_knoise0_n1500}_map1_combine_1_2000.html`.

Conventions: log-HR scale for bias/SD/SE/half-widths (z·SE = 1.96·mean SE, not HR-scale lengths); two-sided coverage of β(Ĥ) with Wilson intervals primary; one-sided lower supplementary (field: stored `lower_1s`; others exp(log est − 1.645·SE)); θ†/θ‡ coverage reported, not scored; MR (field) = est2 with the two-sided Λ-quantile interval, Ĥ block only (F3 convention); M_eff/p̂_Ĥ are not recorded by the committed template (noted per task — the uniform task's recorder extension carries them via M/mass when its knob is on).

## Per-cell tables

### Cell c1 — h1.50, n=500, knoise0: 2000 replicates, 1822 detections (91.1%)

**Harm Ĥ** (θ† = 1.509, θ‡ = 1.671)

| Estimator | n | bias vs β(Ĥ) log-HR | in SD units | SD_β | mean SE_β | SE/SD | Cov_β(Ĥ) 2s (Wilson) | Cov_β(Ĥ) 1s | z·SE (log) | Cov_θ† | Cov_θ‡ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1822 | +0.5802 | +2.92 | 0.199 | 0.315 | 1.59 | 0.529 (0.506, 0.552) | 0.417 | 0.618 | 0.973 | 0.991 |
| oracle | 1822 | +0.4813 | +1.65 | 0.291 | 0.313 | 1.07 | 0.682 (0.660, 0.703) | 0.558 | 0.613 | 0.963 | 0.968 |
| MR (IJ) | 1822 | +0.0586 | +0.21 | 0.276 | 0.395 | 1.43 | 0.977 (0.969, 0.983) | 0.978 | 0.773 | 0.902 | 0.823 |
| MR (field) | 1822 | -0.0209 | -0.07 | 0.316 | 0.359 | 1.14 | 0.919 (0.905, 0.930) | 0.965 | 0.704 | 0.654 | 0.537 |

**Complement Ĥᶜ** (θ† = 0.657, θ‡ = 0.585)

| Estimator | n | bias vs β(Ĥ) log-HR | in SD units | SD_β | mean SE_β | SE/SD | Cov_β(Ĥ) 2s (Wilson) | Cov_β(Ĥ) 1s | z·SE (log) | Cov_θ† | Cov_θ‡ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1822 | -0.0806 | -0.59 | 0.136 | 0.139 | 1.02 | 0.915 (0.902, 0.927) | 0.990 | 0.272 | 0.935 | 0.937 |
| oracle | 1822 | -0.0555 | -0.40 | 0.139 | 0.139 | 1.00 | 0.943 (0.932, 0.953) | 0.985 | 0.273 | 0.946 | 0.911 |
| MR (IJ) | 1822 | -0.0148 | -0.10 | 0.141 | 0.256 | 1.81 | 1.000 (0.998, 1.000) | 0.999 | 0.502 | 1.000 | 0.997 |

Field: draw usage n_out mean 989 (min 817), n_in mean 499; field secs mean 15.9; fit+MR+field secs mean 39.1.

### Cell c2 — h0.75 (protective), n=500: 2000 replicates, 1042 detections (52.1%)

**Harm Ĥ** (θ† = 0.753, θ‡ = 0.701)

| Estimator | n | bias vs β(Ĥ) log-HR | in SD units | SD_β | mean SE_β | SE/SD | Cov_β(Ĥ) 2s (Wilson) | Cov_β(Ĥ) 1s | z·SE (log) | Cov_θ† | Cov_θ‡ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1042 | +0.9051 | +6.08 | 0.149 | 0.329 | 2.21 | 0.001 (0.000, 0.005) | 0.000 | 0.646 | 0.259 | 0.004 |
| oracle | 1042 | +0.2151 | +0.63 | 0.339 | 0.341 | 1.00 | 0.916 (0.897, 0.931) | 0.845 | 0.668 | 0.955 | 0.938 |
| MR (IJ) | 1042 | +0.3224 | +1.57 | 0.205 | 0.402 | 1.96 | 0.984 (0.974, 0.990) | 0.964 | 0.787 | 0.999 | 0.993 |
| MR (field) | 1042 | +0.2123 | +0.91 | 0.234 | 0.359 | 1.54 | 0.988 (0.980, 0.993) | 0.953 | 0.705 | 0.996 | 0.994 |

**Complement Ĥᶜ** (θ† = 0.657, θ‡ = 0.585)

| Estimator | n | bias vs β(Ĥ) log-HR | in SD units | SD_β | mean SE_β | SE/SD | Cov_β(Ĥ) 2s (Wilson) | Cov_β(Ĥ) 1s | z·SE (log) | Cov_θ† | Cov_θ‡ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1042 | -0.0848 | -0.65 | 0.131 | 0.140 | 1.07 | 0.923 (0.905, 0.938) | 0.992 | 0.275 | 0.893 | 0.968 |
| oracle | 1042 | +0.0410 | +0.31 | 0.130 | 0.139 | 1.06 | 0.960 (0.946, 0.970) | 0.929 | 0.272 | 0.964 | 0.858 |
| MR (IJ) | 1042 | -0.0170 | -0.12 | 0.138 | 0.255 | 1.84 | 1.000 (0.996, 1.000) | 0.999 | 0.499 | 0.998 | 0.999 |

Field: draw usage n_out mean 966 (min 653), n_in mean 497; field secs mean 13.9; fit+MR+field secs mean 21.7.

### Cell c3 — h1.00, n=1000: 2000 replicates, 1319 detections (66.0%)

**Harm Ĥ** (θ† = 1.000, θ‡ = 1.000)

| Estimator | n | bias vs β(Ĥ) log-HR | in SD units | SD_β | mean SE_β | SE/SD | Cov_β(Ĥ) 2s (Wilson) | Cov_β(Ĥ) 1s | z·SE (log) | Cov_θ† | Cov_θ‡ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1319 | +0.6043 | +4.26 | 0.142 | 0.252 | 1.78 | 0.269 (0.246, 0.294) | 0.161 | 0.494 | 0.910 | 0.909 |
| oracle | 1319 | +0.3502 | +1.60 | 0.219 | 0.226 | 1.03 | 0.698 (0.673, 0.722) | 0.560 | 0.443 | 0.949 | 0.949 |
| MR (IJ) | 1319 | +0.1806 | +0.94 | 0.191 | 0.298 | 1.56 | 0.980 (0.971, 0.987) | 0.936 | 0.585 | 0.992 | 0.992 |
| MR (field) | 1319 | +0.1056 | +0.49 | 0.216 | 0.264 | 1.23 | 0.958 (0.945, 0.967) | 0.920 | 0.518 | 0.901 | 0.901 |

**Complement Ĥᶜ** (θ† = 0.657, θ‡ = 0.585)

| Estimator | n | bias vs β(Ĥ) log-HR | in SD units | SD_β | mean SE_β | SE/SD | Cov_β(Ĥ) 2s (Wilson) | Cov_β(Ĥ) 1s | z·SE (log) | Cov_θ† | Cov_θ‡ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1319 | -0.0441 | -0.48 | 0.092 | 0.097 | 1.05 | 0.943 (0.929, 0.954) | 0.989 | 0.189 | 0.941 | 0.892 |
| oracle | 1319 | -0.0172 | -0.18 | 0.095 | 0.098 | 1.03 | 0.966 (0.955, 0.974) | 0.981 | 0.192 | 0.951 | 0.836 |
| MR (IJ) | 1319 | +0.0004 | +0.00 | 0.097 | 0.179 | 1.85 | 1.000 (0.997, 1.000) | 1.000 | 0.351 | 1.000 | 0.995 |

Field: draw usage n_out mean 972 (min 626), n_in mean 497; field secs mean 20.9; fit+MR+field secs mean 52.1.

### Cell c4 — h1.75, knoise3, n=500: 2000 replicates, 1945 detections (97.2%)

**Harm Ĥ** (θ† = 1.769, θ‡ = 2.036)

| Estimator | n | bias vs β(Ĥ) log-HR | in SD units | SD_β | mean SE_β | SE/SD | Cov_β(Ĥ) 2s (Wilson) | Cov_β(Ĥ) 1s | z·SE (log) | Cov_θ† | Cov_θ‡ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1945 | +0.4911 | +2.30 | 0.214 | 0.319 | 1.49 | 0.644 (0.617, 0.669) | 0.539 | 0.626 | 0.979 | 0.992 |
| oracle | 1945 | +0.4701 | +1.55 | 0.304 | 0.311 | 1.02 | 0.685 (0.659, 0.710) | 0.563 | 0.609 | 0.956 | 0.953 |
| MR (IJ) | 1945 | -0.0654 | -0.22 | 0.296 | 0.395 | 1.34 | 0.950 (0.937, 0.961) | 0.982 | 0.774 | 0.857 | 0.698 |
| MR (field) | 1945 | -0.1422 | -0.42 | 0.338 | 0.362 | 1.07 | 0.866 (0.846, 0.883) | 0.981 | 0.709 | 0.549 | 0.401 |

**Complement Ĥᶜ** (θ† = 0.657, θ‡ = 0.585)

| Estimator | n | bias vs β(Ĥ) log-HR | in SD units | SD_β | mean SE_β | SE/SD | Cov_β(Ĥ) 2s (Wilson) | Cov_β(Ĥ) 1s | z·SE (log) | Cov_θ† | Cov_θ‡ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1945 | -0.0823 | -0.60 | 0.137 | 0.138 | 1.01 | 0.899 (0.882, 0.914) | 0.989 | 0.270 | 0.931 | 0.930 |
| oracle | 1945 | -0.0691 | -0.49 | 0.140 | 0.139 | 0.99 | 0.926 (0.911, 0.939) | 0.984 | 0.273 | 0.943 | 0.915 |
| MR (IJ) | 1945 | -0.0110 | -0.08 | 0.140 | 0.256 | 1.83 | 1.000 (0.997, 1.000) | 0.999 | 0.502 | 1.000 | 0.998 |

Field: draw usage n_out mean 998 (min 874), n_in mean 500; field secs mean 24.9; fit+MR+field secs mean 83.8.

### Cell c5 — h1.50, n=1500: 2000 replicates, 1976 detections (98.8%)

**Harm Ĥ** (θ† = 1.509, θ‡ = 1.671)

| Estimator | n | bias vs β(Ĥ) log-HR | in SD units | SD_β | mean SE_β | SE/SD | Cov_β(Ĥ) 2s (Wilson) | Cov_β(Ĥ) 1s | z·SE (log) | Cov_θ† | Cov_θ‡ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1976 | +0.1447 | +0.91 | 0.158 | 0.183 | 1.16 | 0.866 (0.851, 0.881) | 0.808 | 0.359 | 0.974 | 0.949 |
| oracle | 1976 | +0.1560 | +0.91 | 0.171 | 0.176 | 1.03 | 0.830 (0.813, 0.846) | 0.778 | 0.346 | 0.960 | 0.930 |
| MR (IJ) | 1976 | -0.0842 | -0.37 | 0.227 | 0.246 | 1.08 | 0.901 (0.887, 0.913) | 0.981 | 0.482 | 0.830 | 0.712 |
| MR (field) | 1976 | -0.0863 | -0.33 | 0.259 | 0.230 | 0.89 | 0.848 (0.832, 0.863) | 0.956 | 0.451 | 0.737 | 0.610 |

**Complement Ĥᶜ** (θ† = 0.657, θ‡ = 0.585)

| Estimator | n | bias vs β(Ĥ) log-HR | in SD units | SD_β | mean SE_β | SE/SD | Cov_β(Ĥ) 2s (Wilson) | Cov_β(Ĥ) 1s | z·SE (log) | Cov_θ† | Cov_θ‡ |
|---|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1976 | -0.0196 | -0.25 | 0.079 | 0.079 | 1.01 | 0.953 (0.943, 0.962) | 0.977 | 0.156 | 0.935 | 0.831 |
| oracle | 1976 | -0.0264 | -0.34 | 0.077 | 0.080 | 1.04 | 0.937 (0.925, 0.947) | 0.979 | 0.157 | 0.927 | 0.861 |
| MR (IJ) | 1976 | +0.0066 | +0.08 | 0.082 | 0.151 | 1.84 | 0.999 (0.996, 1.000) | 0.998 | 0.296 | 0.999 | 0.988 |

Field: draw usage n_out mean 994 (min 747), n_in mean 499; field secs mean 36.0; fit+MR+field secs mean 145.7.

### Cross-cell summary (harm block, coverage of β(Ĥ))

| Cell | Detection | IJ 2-sided (Wilson) | IJ SE/SD | Field 1-sided | Field 2-sided (Wilson) | λ-SD/SD |
|---|---|---|---|---|---|---|
| h1.50, n=500, knoise0 | 91.1% | 0.977 (0.969, 0.983) | 1.43 | 0.965 | 0.919 (0.905, 0.930) | 1.14 |
| h0.75 (protective), n=500 | 52.1% | 0.984 (0.974, 0.990) | 1.96 | 0.953 | 0.988 (0.980, 0.993) | 1.54 |
| h1.00, n=1000 | 66.0% | 0.980 (0.971, 0.987) | 1.56 | 0.920 | 0.958 (0.945, 0.967) | 1.23 |
| h1.75, knoise3, n=500 | 97.2% | 0.950 (0.937, 0.961) | 1.34 | 0.981 | 0.866 (0.846, 0.883) | 1.07 |
| h1.50, n=1500 | 98.8% | 0.901 (0.887, 0.913) | 1.08 | 0.956 | 0.848 (0.832, 0.863) | 0.89 |

## Findings (in the record; nothing blocks, no task proposed)

1. **The field's one-sided bound holds across the whole map**: 0.920–0.981 against β(Ĥ) in all five cells (0.953 on the protective cell, 0.956 at n=1500) — consistent with the uniform-validity claim it carries.
2. **Moderate harm at large n (c5) is the first regime where MR (IJ)'s own two-sided interval under-covers**: 0.901 (0.887, 0.913) with SE/SD 1.08 — the conservatism that propped up IJ's two-sided coverage elsewhere (SE/SD 1.3–2.0) is gone at n=1500, and the field's plain two-sided interval is lower still (0.848). At n=500 (c1) IJ still covers (0.977) but only via SE/SD 1.43.
3. **The protective cell (c2) is a naive catastrophe and an MR success**: naive conditional coverage 0.001 with retained bias +0.91 log-HR (+6.1 SD); MR (IJ) 0.984 and MR (field) 0.988 two-sided, 0.964/0.953 one-sided. The sweep-era concern about IJ's conditional coverage under a protective truth does not reproduce under the current package on the conditional target.
4. **A larger family (knoise3, c4) behaves like its knoise0 counterpart**: IJ 0.950 vs 0.961 at h175 knoise0 (s7); field two-sided 0.866 vs 0.887 — the two-sided dip persists, slightly deeper, with detection 97.2% and field draw usage healthy (min n_out 874).
5. **Retained bias ordering (field est2 below IJ β̃) holds in every harm block except c5**, where the two are equal (−0.086 vs −0.084) and both overshoot slightly; the c2 protective cell shows the largest remaining MR retained bias (+0.32 IJ / +0.21 field, log-HR).
6. One cross-vintage selection flip surfaced (c3 sim 309, oracle exact) — enumerated in the Gate 2 record, same class as the s7 arc's; cumulative enumerated flips across all anchored comparisons to date: h10 393/780, h175 cut-labels 38/721, h100-n1000 309.
7. Cost: fit+MR+field per replicate at 100 workers — 21.7 s (c2) to 145.7 s (c5, n=1500); field block 14–36 s. Cumulative Stage 2 wall 2 h 25 m of the 5 h ceiling.

## Done means (task §Done)

Stage 3 report + five rendered combined documents committed; Gate 2 record beside them; every run under campaign `map1`; committed documents and bundles untouched; branch unpushed. **Cells completed: all five (1–5). Deferred: none. Dropped: none.**
