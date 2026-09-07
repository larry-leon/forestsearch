# REPORT — The field method on the continuous/MD path (campaign `mdf1`, Stage 3)

Date: 2026-09-07. Machine: Mac Studio (M4 Max, 13 workers). Branch `feature/glm-extension-mac`; forestsearch 0.3.5 (with the add-only `scale` argument on `fs_sim_bias_coverage()`, commit `2f118042`). Task: `dev/tasks/TASK_continuous_field_mac_2026-09-07.md`. Records: Stage 0 `REPORT_continuous_field_stage0_2026-09-07.md`, Stage 1 `REPORT_continuous_field_stage1_2026-09-07.md`, Gate 2 `REPORT_continuous_field_gate2_2026-09-07.md` (all beside this file). Rendered documents: one combine document per cell (`fs_maxeffCons_mr_field_<cell>_mdf1_combine_1_2000.html`) and the cross-cell `summary_continuous_field_mdf1.html` (self-contained; reads the four committed combined bundles beside it). Every table below is computed from those bundles (the session's aggregator reproduces the summary document's tables).

## Design and conventions

- Cells (M-1, the continuous twin's committed cells): md40 n = 500, md120 n = 500, null n = 500, md40 n = 700; 2,000 replicates each (M-2), `sim_id` 1–2,000, seeds `8316951 + sim_id` under L'Ecuyer-CMRG; `sim_id` 1–1,000 are the committed seeds and carry the pairing proof (Gate 2 record). ACTG175 CD4-change DGM, true region `age > 34 & preanti <= 744.5` (prevalence 0.345), structural harm-region MD −40 (md40) or −120 (md120), complement −26.26 raw; the null cell has no treatment-by-subgroup interaction (a homogeneous −26.26). Identifier `maxeffCons` on the consistency engine, thresholds 30 / 10 / 0.90, J = 10 cut grids on `age` and `preanti`, `str2` in the pool, 5,000 multiplier draws, field R_out/R_in = 1000/500, complement field on, two-term IJ, FB joined on md40 n = 500 `sim_id` 1–100 only.
- **Scale and sides.** All estimator columns are oriented (positive = harm; the gate works on −`cd4_change`); β(Ĥ), β(Ĥᶜ) are the exact super-population targets at the realized rule, oriented with −1. The harm block's directional product is the one-sided 95% **lower** bound ("harm ≥ L"); the complement's is the one-sided 95% **upper** bound ("harm ≤ U", i.e. benefit ≥ −U raw). The field rows use their stored `fld_H_lo1s` / `fld_Hc_up1s` (Λ* q95 / q05 inversion about β̃); the normal-based rows use est ∓ 1.645 SE. The oracle is scored against the structural true-region effect; under the null the harm oracle is undefined (Q is empty) and its row is blank.
- **Reading convention.** Bounds are read by location against MD thresholds (0 = no harm; 10 = the design's consistency threshold; 20; 30 = the search's effect threshold; 40 = the planted MD), never as significance at the null. Coverage carries Wilson 95% intervals; bias in MD units and in SD units (bias / empirical SD). Rows: naive, oracle, MR (IJ two-term), MR (field); the joint pair separate / Bonferroni / calibrated. The winner-only and winner-floor variants are excluded from every table (recorder columns exist).
- Survival ranges quoted for comparison are the handoff's §3 (12.5%-prevalence cells s7/map1/s7c; 31% prevalence effMaxSG p30sg), all at 2,000 replicates.

## Table-2 layout, both blocks (per cell; rows naive / oracle / MR (IJ) / MR (field))


| cell | block | estimator | n | bias | bias_sd | SD | SE | SE_SD | cov2_w | cov1_w | side | halfwidth | margin1 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| md40 n500 | H | naive | 1998.000 | 74.656 | 4.134 | 18.061 | 31.086 | 1.721 | 0.232 (0.214, 0.251) | 0.101 (0.088, 0.115) | lower | 60.928 | 51.132 |
| md40 n500 | H | oracle | 1998.000 | -0.660 | -0.033 | 19.985 | 19.806 | 0.991 | 0.943 (0.932, 0.952) | 0.952 (0.942, 0.961) | lower | 38.820 | 32.578 |
| md40 n500 | H | MR (IJ) | 1998.000 | 20.711 | 0.924 | 22.424 | 38.904 | 1.735 | 0.984 (0.977, 0.989) | 0.964 (0.955, 0.972) | lower | 76.251 | 63.992 |
| md40 n500 | H | MR (field) | 1998.000 | 9.629 | 0.392 | 24.578 | 29.547 | 1.202 | 0.978 (0.971, 0.984) | 0.950 (0.940, 0.959) | lower | 57.554 | 62.126 |
| md40 n500 | Hc | naive | 1998.000 | -12.740 | -1.069 | 11.918 | 12.556 | 1.054 | 0.845 (0.829, 0.861) | 0.760 (0.741, 0.778) | upper | 24.610 | 20.653 |
| md40 n500 | Hc | oracle | 1998.000 | -0.262 | -0.018 | 14.343 | 14.418 | 1.005 | 0.953 (0.943, 0.962) | 0.957 (0.947, 0.965) | upper | 28.260 | 23.716 |
| md40 n500 | Hc | MR (IJ) | 1998.000 | -4.385 | -0.363 | 12.087 | 23.521 | 1.946 | 0.999 (0.997, 1.000) | 0.994 (0.990, 0.997) | upper | 46.099 | 38.688 |
| md40 n500 | Hc | MR (field) | 1998.000 | -2.680 | -0.219 | 12.224 | 12.392 | 1.014 | 0.943 (0.932, 0.952) | 0.924 (0.911, 0.935) | upper | 24.186 | 22.101 |
| md120 n500 | H | naive | 2000.000 | 57.606 | 2.685 | 21.451 | 30.927 | 1.442 | 0.564 (0.542, 0.585) | 0.408 (0.387, 0.430) | lower | 60.615 | 50.870 |
| md120 n500 | H | oracle | 2000.000 | -0.684 | -0.034 | 19.989 | 19.807 | 0.991 | 0.943 (0.932, 0.952) | 0.953 (0.942, 0.961) | lower | 38.821 | 32.579 |
| md120 n500 | H | MR (IJ) | 2000.000 | 9.898 | 0.360 | 27.485 | 40.828 | 1.485 | 0.990 (0.984, 0.993) | 0.970 (0.961, 0.976) | lower | 80.022 | 67.156 |
| md120 n500 | H | MR (field) | 2000.000 | 1.840 | 0.061 | 30.387 | 31.412 | 1.034 | 0.942 (0.930, 0.951) | 0.948 (0.937, 0.957) | lower | 61.074 | 62.113 |
| md120 n500 | Hc | naive | 2000.000 | -10.045 | -0.781 | 12.869 | 12.823 | 0.996 | 0.892 (0.878, 0.905) | 0.829 (0.811, 0.844) | upper | 25.132 | 21.092 |
| md120 n500 | Hc | oracle | 2000.000 | -0.301 | -0.021 | 14.389 | 14.418 | 1.002 | 0.953 (0.942, 0.961) | 0.956 (0.946, 0.964) | upper | 28.259 | 23.715 |
| md120 n500 | Hc | MR (IJ) | 2000.000 | -2.657 | -0.202 | 13.136 | 24.078 | 1.833 | 1.000 (0.998, 1.000) | 0.997 (0.993, 0.998) | upper | 47.193 | 39.606 |
| md120 n500 | Hc | MR (field) | 2000.000 | -1.416 | -0.106 | 13.345 | 12.739 | 0.955 | 0.946 (0.935, 0.955) | 0.934 (0.923, 0.945) | upper | 24.855 | 22.202 |
| null n500 | H | naive | 1993.000 | 75.217 | 4.215 | 17.845 | 31.112 | 1.744 | 0.224 (0.206, 0.243) | 0.093 (0.081, 0.106) | lower | 60.979 | 51.175 |
| null n500 | H | oracle | 0.000 |  NaN |   NA |   NA |   NA |   NA | NaN (NA, NA) | NaN (NA, NA) | lower |   NA |   NA |
| null n500 | H | MR (IJ) | 1993.000 | 21.068 | 0.949 | 22.200 | 38.871 | 1.751 | 0.984 (0.978, 0.989) | 0.965 (0.956, 0.972) | lower | 76.186 | 63.937 |
| null n500 | H | MR (field) | 1993.000 | 9.897 | 0.406 | 24.368 | 29.527 | 1.212 | 0.978 (0.971, 0.984) | 0.950 (0.940, 0.959) | lower | 57.511 | 62.167 |
| null n500 | Hc | naive | 1993.000 | -12.777 | -1.084 | 11.782 | 12.557 | 1.066 | 0.846 (0.829, 0.861) | 0.759 (0.740, 0.777) | upper | 24.611 | 20.655 |
| null n500 | Hc | oracle | 1993.000 | -0.317 | -0.028 | 11.363 | 11.653 | 1.026 | 0.956 (0.946, 0.964) | 0.953 (0.943, 0.962) | upper | 22.839 | 19.167 |
| null n500 | Hc | MR (IJ) | 1993.000 | -4.406 | -0.368 | 11.966 | 23.514 | 1.965 | 0.999 (0.997, 1.000) | 0.996 (0.992, 0.998) | upper | 46.087 | 38.677 |
| null n500 | Hc | MR (field) | 1993.000 | -2.692 | -0.222 | 12.107 | 12.393 | 1.024 | 0.945 (0.934, 0.954) | 0.926 (0.913, 0.936) | upper | 24.181 | 22.105 |
| md40 n700 | H | naive | 1999.000 | 76.442 | 4.359 | 17.536 | 30.959 | 1.765 | 0.189 (0.172, 0.206) | 0.067 (0.056, 0.078) | lower | 60.679 | 50.924 |
| md40 n700 | H | oracle | 1999.000 | 0.074 | 0.004 | 17.339 | 16.771 | 0.967 | 0.943 (0.932, 0.952) | 0.941 (0.930, 0.950) | lower | 32.872 | 27.586 |
| md40 n700 | H | MR (IJ) | 1999.000 | 21.350 | 0.961 | 22.224 | 36.904 | 1.661 | 0.974 (0.967, 0.981) | 0.942 (0.931, 0.952) | lower | 72.330 | 60.702 |
| md40 n700 | H | MR (field) | 1999.000 | 10.069 | 0.410 | 24.564 | 29.205 | 1.189 | 0.970 (0.962, 0.977) | 0.947 (0.936, 0.956) | lower | 56.766 | 61.905 |
| md40 n700 | Hc | naive | 1999.000 | -8.590 | -0.829 | 10.359 | 10.384 | 1.002 | 0.870 (0.855, 0.884) | 0.797 (0.779, 0.814) | upper | 20.351 | 17.080 |
| md40 n700 | Hc | oracle | 1999.000 | 0.186 | 0.015 | 12.100 | 12.179 | 1.007 | 0.950 (0.940, 0.959) | 0.955 (0.945, 0.963) | upper | 23.871 | 20.033 |
| md40 n700 | Hc | MR (IJ) | 1999.000 | -2.628 | -0.252 | 10.424 | 19.812 | 1.901 | 0.999 (0.997, 1.000) | 0.996 (0.993, 0.998) | upper | 38.830 | 32.587 |
| md40 n700 | Hc | MR (field) | 1999.000 | -1.380 | -0.131 | 10.500 | 10.296 | 0.981 | 0.939 (0.928, 0.949) | 0.936 (0.925, 0.946) | upper | 20.092 | 18.175 | 


Compact format (bias / SD / SE / SE-to-SD / two-sided / one-sided):


| cell | block | estimator | bias (MD | SD) | SD | SE | SE/SD | two-sided | one-sided |
|---|---|---|---|---|---|---|---|---|
| md40 n500 | H | naive | +74.66 | +4.134 | 18.06 | 31.09 | 1.721 | 0.232 (0.214, 0.251) | 0.101 (0.088, 0.115) lower |
| md40 n500 | H | oracle | -0.66 | -0.033 | 19.98 | 19.81 | 0.991 | 0.943 (0.932, 0.952) | 0.952 (0.942, 0.961) lower |
| md40 n500 | H | MR (IJ) | +20.71 | +0.924 | 22.42 | 38.90 | 1.735 | 0.984 (0.977, 0.989) | 0.964 (0.955, 0.972) lower |
| md40 n500 | H | MR (field) | +9.63 | +0.392 | 24.58 | 29.55 | 1.202 | 0.978 (0.971, 0.984) | 0.950 (0.940, 0.959) lower |
| md40 n500 | Hc | naive | -12.74 | -1.069 | 11.92 | 12.56 | 1.054 | 0.845 (0.829, 0.861) | 0.760 (0.741, 0.778) upper |
| md40 n500 | Hc | oracle | -0.26 | -0.018 | 14.34 | 14.42 | 1.005 | 0.953 (0.943, 0.962) | 0.957 (0.947, 0.965) upper |
| md40 n500 | Hc | MR (IJ) | -4.38 | -0.363 | 12.09 | 23.52 | 1.946 | 0.999 (0.997, 1.000) | 0.994 (0.990, 0.997) upper |
| md40 n500 | Hc | MR (field) | -2.68 | -0.219 | 12.22 | 12.39 | 1.014 | 0.943 (0.932, 0.952) | 0.924 (0.911, 0.935) upper |
| md120 n500 | H | naive | +57.61 | +2.685 | 21.45 | 30.93 | 1.442 | 0.564 (0.542, 0.585) | 0.408 (0.387, 0.430) lower |
| md120 n500 | H | oracle | -0.68 | -0.034 | 19.99 | 19.81 | 0.991 | 0.943 (0.932, 0.952) | 0.953 (0.942, 0.961) lower |
| md120 n500 | H | MR (IJ) | +9.90 | +0.360 | 27.49 | 40.83 | 1.485 | 0.990 (0.984, 0.993) | 0.970 (0.961, 0.976) lower |
| md120 n500 | H | MR (field) | +1.84 | +0.061 | 30.39 | 31.41 | 1.034 | 0.942 (0.930, 0.951) | 0.948 (0.937, 0.957) lower |
| md120 n500 | Hc | naive | -10.04 | -0.781 | 12.87 | 12.82 | 0.996 | 0.892 (0.878, 0.905) | 0.829 (0.811, 0.844) upper |
| md120 n500 | Hc | oracle | -0.30 | -0.021 | 14.39 | 14.42 | 1.002 | 0.953 (0.942, 0.961) | 0.956 (0.946, 0.964) upper |
| md120 n500 | Hc | MR (IJ) | -2.66 | -0.202 | 13.14 | 24.08 | 1.833 | 1.000 (0.998, 1.000) | 0.997 (0.993, 0.998) upper |
| md120 n500 | Hc | MR (field) | -1.42 | -0.106 | 13.35 | 12.74 | 0.955 | 0.946 (0.935, 0.955) | 0.934 (0.923, 0.945) upper |
| null n500 | H | naive | +75.22 | +4.215 | 17.84 | 31.11 | 1.744 | 0.224 (0.206, 0.243) | 0.093 (0.081, 0.106) lower |
| null n500 | H | oracle | NaN | NA | NA | NA | NA | NaN (NA, NA) | NaN (NA, NA) lower |
| null n500 | H | MR (IJ) | +21.07 | +0.949 | 22.20 | 38.87 | 1.751 | 0.984 (0.978, 0.989) | 0.965 (0.956, 0.972) lower |
| null n500 | H | MR (field) | +9.90 | +0.406 | 24.37 | 29.53 | 1.212 | 0.978 (0.971, 0.984) | 0.950 (0.940, 0.959) lower |
| null n500 | Hc | naive | -12.78 | -1.084 | 11.78 | 12.56 | 1.066 | 0.846 (0.829, 0.861) | 0.759 (0.740, 0.777) upper |
| null n500 | Hc | oracle | -0.32 | -0.028 | 11.36 | 11.65 | 1.026 | 0.956 (0.946, 0.964) | 0.953 (0.943, 0.962) upper |
| null n500 | Hc | MR (IJ) | -4.41 | -0.368 | 11.97 | 23.51 | 1.965 | 0.999 (0.997, 1.000) | 0.996 (0.992, 0.998) upper |
| null n500 | Hc | MR (field) | -2.69 | -0.222 | 12.11 | 12.39 | 1.024 | 0.945 (0.934, 0.954) | 0.926 (0.913, 0.936) upper |
| md40 n700 | H | naive | +76.44 | +4.359 | 17.54 | 30.96 | 1.765 | 0.189 (0.172, 0.206) | 0.067 (0.056, 0.078) lower |
| md40 n700 | H | oracle | +0.07 | +0.004 | 17.34 | 16.77 | 0.967 | 0.943 (0.932, 0.952) | 0.941 (0.930, 0.950) lower |
| md40 n700 | H | MR (IJ) | +21.35 | +0.961 | 22.22 | 36.90 | 1.661 | 0.974 (0.967, 0.981) | 0.942 (0.931, 0.952) lower |
| md40 n700 | H | MR (field) | +10.07 | +0.410 | 24.56 | 29.21 | 1.189 | 0.970 (0.962, 0.977) | 0.947 (0.936, 0.956) lower |
| md40 n700 | Hc | naive | -8.59 | -0.829 | 10.36 | 10.38 | 1.002 | 0.870 (0.855, 0.884) | 0.797 (0.779, 0.814) upper |
| md40 n700 | Hc | oracle | +0.19 | +0.015 | 12.10 | 12.18 | 1.007 | 0.950 (0.940, 0.959) | 0.955 (0.945, 0.963) upper |
| md40 n700 | Hc | MR (IJ) | -2.63 | -0.252 | 10.42 | 19.81 | 1.901 | 0.999 (0.997, 1.000) | 0.996 (0.993, 0.998) upper |
| md40 n700 | Hc | MR (field) | -1.38 | -0.131 | 10.50 | 10.30 | 0.981 | 0.939 (0.928, 0.949) | 0.936 (0.925, 0.946) upper | 


## Bound location against MD thresholds (reading aids: 0 / 10 / 20 / 30 / 40 CD4 cells/mm³, oriented)

Harm block, one-sided 95% LOWER bound: location summary and the share of replicates whose bound is at or above each threshold ("harm at least τ" supported). For orientation, the mean oriented β(Ĥ) is +31.7 in the md40 cells (median 31.4; Ĥ overlaps Q weakly: sensitivity 0.17 / 0.12, PPV 0.40), +26.3 under the null (no subgroup; homogeneous +26 on the harm-oriented scale, so the null cell's bound locations are read against +26, not as a false-claim rate) and +96.0 in md120 (median 94.4, q10–q90 64–120; sensitivity 0.33, PPV 0.76).


| cell | block | estimator | n | mean | median | q10 | q90 | P(L>=0) | P(L>=10) | P(L>=20) | P(L>=30) | P(L>=40) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| md40 n500 | H | naive | 1998.000 | 55.193 | 54.555 | 31.810 | 78.850 | 1.000 | 0.999 | 0.981 | 0.918 | 0.797 |
| md40 n500 | H | MR (IJ) | 1998.000 | -11.612 | -13.145 | -39.874 | 18.203 | 0.282 | 0.170 | 0.088 | 0.042 | 0.020 |
| md40 n500 | H | MR (field) | 1998.000 | -9.746 | -11.398 | -37.017 | 20.761 | 0.307 | 0.185 | 0.107 | 0.054 | 0.025 |
| md120 n500 | H | naive | 2000.000 | 102.738 | 101.856 | 75.021 | 131.026 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| md120 n500 | H | MR (IJ) | 2000.000 | 38.743 | 37.405 | 5.624 | 73.244 | 0.942 | 0.864 | 0.753 | 0.603 | 0.468 |
| md120 n500 | H | MR (field) | 2000.000 | 43.787 | 41.676 | 9.360 | 81.077 | 0.956 | 0.896 | 0.799 | 0.661 | 0.524 |
| null n500 | H | naive | 1993.000 | 50.298 | 49.769 | 27.281 | 73.881 | 1.000 | 0.999 | 0.960 | 0.867 | 0.709 |
| null n500 | H | MR (IJ) | 1993.000 | -16.614 | -18.294 | -44.618 | 12.997 | 0.211 | 0.120 | 0.065 | 0.030 | 0.015 |
| null n500 | H | MR (field) | 1993.000 | -14.844 | -16.708 | -42.663 | 15.641 | 0.240 | 0.138 | 0.076 | 0.037 | 0.017 |
| md40 n700 | H | naive | 1999.000 | 57.191 | 56.353 | 35.339 | 79.101 | 1.000 | 0.999 | 0.986 | 0.946 | 0.848 |
| md40 n700 | H | MR (IJ) | 1999.000 | -7.679 | -9.481 | -34.745 | 20.732 | 0.337 | 0.194 | 0.103 | 0.056 | 0.035 |
| md40 n700 | H | MR (field) | 1999.000 | -8.882 | -10.810 | -35.174 | 19.716 | 0.309 | 0.185 | 0.097 | 0.059 | 0.035 | 


Complement block, one-sided 95% UPPER bound: share at or below each threshold ("harm at most τ" supported). The complement's oriented β(Ĥᶜ) is about +31 in the md40 cells, +26 under the null (no subgroup; homogeneous +26 on the harm-oriented scale) and +52 in md120 (the treatment lowers CD4 change everywhere on this DGM, and in md120 two thirds of Q sit in Ĥᶜ), so a benefit claim on Ĥᶜ is not supported at any threshold, and the shares below 30 are correctly small.


| cell | block | estimator | n | mean | median | q10 | q90 | P(U<=0) | P(U<=10) | P(U<=20) | P(U<=30) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| md40 n500 | Hc | naive | 1998.000 | 38.796 | 38.894 | 23.288 | 53.976 | 0.001 | 0.011 | 0.058 | 0.216 |
| md40 n500 | Hc | MR (IJ) | 1998.000 | 65.185 | 65.389 | 49.673 | 80.818 | 0.000 | 0.000 | 0.000 | 0.004 |
| md40 n500 | Hc | MR (field) | 1998.000 | 48.598 | 48.832 | 32.887 | 64.423 | 0.000 | 0.002 | 0.015 | 0.068 |
| md120 n500 | Hc | naive | 2000.000 | 62.625 | 62.896 | 46.388 | 79.074 | 0.000 | 0.001 | 0.001 | 0.007 |
| md120 n500 | Hc | MR (IJ) | 2000.000 | 88.527 | 88.840 | 71.690 | 105.217 | 0.000 | 0.000 | 0.000 | 0.000 |
| md120 n500 | Hc | MR (field) | 2000.000 | 71.123 | 71.514 | 54.277 | 88.340 | 0.000 | 0.000 | 0.001 | 0.002 |
| null n500 | Hc | naive | 1993.000 | 34.133 | 34.081 | 19.056 | 49.490 | 0.003 | 0.022 | 0.114 | 0.353 |
| null n500 | Hc | MR (IJ) | 1993.000 | 60.526 | 60.590 | 44.881 | 76.266 | 0.000 | 0.000 | 0.001 | 0.009 |
| null n500 | Hc | MR (field) | 1993.000 | 43.955 | 44.153 | 28.250 | 59.483 | 0.000 | 0.006 | 0.029 | 0.120 |
| md40 n700 | Hc | naive | 1999.000 | 39.398 | 39.138 | 26.301 | 52.775 | 0.000 | 0.003 | 0.030 | 0.182 |
| md40 n700 | Hc | MR (IJ) | 1999.000 | 60.868 | 60.534 | 48.170 | 74.475 | 0.000 | 0.000 | 0.000 | 0.003 |
| md40 n700 | Hc | MR (field) | 1999.000 | 46.456 | 46.228 | 33.073 | 60.076 | 0.000 | 0.001 | 0.006 | 0.054 | 


## Joint pair (Ĥ lower, Ĥᶜ upper)


| cell | pair | n | joint | cov_H | cov_Hc | margin_H | margin_Hc |
|---|---|---|---|---|---|---|---|
| md40 n500 | separate 95% field bounds | 1998.000 | 0.879 (0.864, 0.892) | 0.950 | 0.924 | 62.126 | 22.101 |
| md40 n500 | Bonferroni (gamma = 0.025) | 1998.000 | 0.940 (0.929, 0.950) | 0.981 | 0.958 | 73.422 | 26.002 |
| md40 n500 | calibrated gamma | 1998.000 | 0.940 (0.929, 0.950) | 0.981 | 0.957 | 73.339 | 25.975 |
| md120 n500 | separate 95% field bounds | 2000.000 | 0.886 (0.872, 0.900) | 0.948 | 0.934 | 62.113 | 22.202 |
| md120 n500 | Bonferroni (gamma = 0.025) | 2000.000 | 0.943 (0.932, 0.952) | 0.977 | 0.965 | 73.654 | 26.202 |
| md120 n500 | calibrated gamma | 2000.000 | 0.943 (0.931, 0.952) | 0.977 | 0.965 | 73.550 | 26.162 |
| null n500 | separate 95% field bounds | 1993.000 | 0.880 (0.865, 0.894) | 0.950 | 0.926 | 62.167 | 22.105 |
| null n500 | Bonferroni (gamma = 0.025) | 1993.000 | 0.940 (0.928, 0.949) | 0.979 | 0.959 | 73.457 | 26.015 |
| null n500 | calibrated gamma | 1993.000 | 0.939 (0.928, 0.949) | 0.979 | 0.959 | 73.374 | 25.985 |
| md40 n700 | separate 95% field bounds | 1999.000 | 0.884 (0.870, 0.898) | 0.947 | 0.936 | 61.905 | 18.175 |
| md40 n700 | Bonferroni (gamma = 0.025) | 1999.000 | 0.936 (0.924, 0.946) | 0.972 | 0.963 | 73.071 | 21.425 |
| md40 n700 | calibrated gamma | 1999.000 | 0.935 (0.924, 0.945) | 0.971 | 0.963 | 72.984 | 21.399 | 


## Regime diagnostics (per cell)

`p_hat_*`: p̂(Ĥ), the winner's re-selection frequency (mean, median, share < 0.5 = tie regime, share ≥ 0.9 = settled); `gamma_mean`: the calibrated joint γ; `corr_lam`: corr(Λ*, Λ*ᶜ); `nfit_mean` / `share_newfit`: complement fits per replicate and the share of draw-winner readings needing a lazy fit; `sd_btc_naive`: SD(β̃ᶜ)/naive SE; `lamc_naive`: λ-SDᶜ/naive SE; `ij_sd_*`: IJ SE / empirical SD per block; seconds per replicate (fit + MR + field; field; complement) at 13 workers.


| cell | n_det | mean_pH | p_hat_mean | p_hat_med | p_hat_lt05 | p_hat_ge09 | gamma_mean | corr_lam | nfit_mean | share_newfit | sd_btc_naive | lamc_naive | ij_sd_H | ij_sd_Hc | fit_secs | field_secs | comp_secs |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| md40 n500 | 1998.000 | 72.078 | 0.155 | 0.134 | 0.990 | 0.000 | 0.025 | 0.091 | 581.506 | 0.057 | 0.963 | 0.987 | 1.735 | 1.946 | 14.953 | 10.982 | 0.399 |
| md120 n500 | 2000.000 | 74.843 | 0.189 | 0.166 | 0.974 | 0.000 | 0.025 | 0.045 | 522.751 | 0.062 | 1.024 | 0.993 | 1.485 | 1.833 | 17.151 | 12.502 | 0.379 |
| null n500 | 1993.000 | 72.077 | 0.155 | 0.133 | 0.987 | 0.000 | 0.025 | 0.091 | 581.906 | 0.057 | 0.953 | 0.987 | 1.751 | 1.965 | 14.744 | 10.876 | 0.410 |
| md40 n700 | 1999.000 | 73.257 | 0.163 | 0.146 | 0.986 | 0.000 | 0.025 | 0.084 | 545.793 | 0.054 | 1.004 | 0.992 | 1.661 | 1.901 | 17.648 | 12.816 | 0.420 | 


## The display (identity scale, both blocks)

`fs_sim_bias_coverage(scale = "identity")`: b = retained bias / empirical SD, r = mean SE / empirical SD; `cov1_ref` = Φ(1.645·r − b) (harm block, lower bound) or Φ(1.645·r + b) (complement, upper bound); `cov2_ref` = Φ(1.96·r − b) − Φ(−1.96·r − b). Figures: `fig_mdf1_bias_coverage_display_H.png`, `fig_mdf1_bias_coverage_display_Hc.png` (also embedded in `summary_continuous_field_mdf1.html`).


| cell | block | estimator | n | bias_log | sd_emp | se_mean | b | r | cov1 | cov1_ref | cov2 | cov2_ref |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| md40 n500 | H | mr | 1998.000 | 20.711 | 22.424 | 38.904 | 0.924 | 1.735 | 0.964 | 0.973 | 0.984 | 0.993 |
| md40 n500 | H | fld | 1998.000 | 9.629 | 24.578 | 29.547 | 0.392 | 1.202 | 0.950 | 0.944 | 0.978 | 0.972 |
| md40 n500 | Hc | mr | 1998.000 | -4.385 | 12.087 | 23.521 | -0.363 | 1.946 | 0.994 | 0.998 | 0.999 | 1.000 |
| md40 n500 | Hc | fld | 1998.000 | -2.680 | 12.224 | 12.392 | -0.219 | 1.014 | 0.924 | 0.926 | 0.943 | 0.948 |
| md120 n500 | H | mr | 2000.000 | 9.898 | 27.485 | 40.828 | 0.360 | 1.485 | 0.970 | 0.981 | 0.990 | 0.994 |
| md120 n500 | H | fld | 2000.000 | 1.840 | 30.387 | 31.412 | 0.061 | 1.034 | 0.948 | 0.949 | 0.942 | 0.957 |
| md120 n500 | Hc | mr | 2000.000 | -2.657 | 13.136 | 24.078 | -0.202 | 1.833 | 0.997 | 0.998 | 1.000 | 1.000 |
| md120 n500 | Hc | fld | 2000.000 | -1.416 | 13.345 | 12.739 | -0.106 | 0.955 | 0.934 | 0.928 | 0.946 | 0.937 |
| null n500 | H | mr | 1993.000 | 21.068 | 22.200 | 38.871 | 0.949 | 1.751 | 0.965 | 0.973 | 0.984 | 0.993 |
| null n500 | H | fld | 1993.000 | 9.897 | 24.368 | 29.527 | 0.406 | 1.212 | 0.950 | 0.944 | 0.978 | 0.973 |
| null n500 | Hc | mr | 1993.000 | -4.406 | 11.966 | 23.514 | -0.368 | 1.965 | 0.996 | 0.998 | 0.999 | 1.000 |
| null n500 | Hc | fld | 1993.000 | -2.692 | 12.107 | 12.393 | -0.222 | 1.024 | 0.926 | 0.928 | 0.945 | 0.950 |
| md40 n700 | H | mr | 1999.000 | 21.350 | 22.224 | 36.904 | 0.961 | 1.661 | 0.942 | 0.962 | 0.974 | 0.989 |
| md40 n700 | H | fld | 1999.000 | 10.069 | 24.564 | 29.205 | 0.410 | 1.189 | 0.947 | 0.939 | 0.970 | 0.970 |
| md40 n700 | Hc | mr | 1999.000 | -2.628 | 10.424 | 19.812 | -0.252 | 1.901 | 0.996 | 0.998 | 0.999 | 1.000 |
| md40 n700 | Hc | fld | 1999.000 | -1.380 | 10.500 | 10.296 | -0.131 | 0.981 | 0.936 | 0.931 | 0.939 | 0.943 | 

## Findings (reading criteria per the task; Larry's criteria, not gates)

**1. The field's one-sided bounds are calibrated on both blocks, and the continuous path is the cleanest yet.** Harm block, field one-sided LOWER coverage of β(Ĥ): **0.950, 0.948, 0.950, 0.947** (md40 n500, md120 n500, null n500, md40 n700; Wilson half-widths ≈ 0.010). Every cell sits on the nominal level, inside the survival ranges (0.920–0.981 at 12.5% prevalence; 0.946–0.971 at 31%) and without the survival tie cell's 0.920 shortfall — even though **every continuous cell is a tie regime** (finding 5). Complement, field one-sided UPPER coverage of β(Ĥᶜ): **0.924, 0.934, 0.926, 0.936** (Wilson 0.911–0.946), matching the survival n = 500 complement (0.930–0.937) and the p30sg range (0.911–0.932), slightly under nominal at n = 500 and closer at n = 700 (0.936), the direction the survival n ≥ 1,000 cells took (0.953–0.956). λ-SDᶜ/naive SE 0.987–0.993 (survival 0.97–0.99); SD(β̃ᶜ)/naive SE 0.95–1.02 — the complement is in the survival n = 500 regime, not the p30sg moved regime (1.08–1.12).

**2. Two-sided intervals.** MR (IJ two-term) two-sided coverage 0.984 / 0.990 / 0.984 / 0.974 with SE/SD **1.74 / 1.49 / 1.75 / 1.66** on Ĥ (survival 1.08–1.96) and **1.95 / 1.83 / 1.97 / 1.90** on Ĥᶜ (survival 1.80–1.85): conservative by construction, as on survival, with the complement's IJ upper bound uninformative (coverage 0.994–0.997; mean upper bound 61–89 oriented against a target of 26–52). The field's two-sided interval covers at **0.978 / 0.942 / 0.978 / 0.970** — no under-coverage on this path (survival harm cells: 0.85–0.92), the md120 cell nearest nominal. Field half-widths are 25% narrower than IJ's on Ĥ (57 vs 76 MD units at n = 500) and half IJ's on Ĥᶜ (24 vs 46); the field's one-sided margin on Ĥ (62) is close to IJ's (64) because it inverts about β̃ with an asymmetric Λ*.

**3. Retained bias (SD units), Ĥ: naive → IJ → field = +4.1 → +0.92 → +0.39** (md40 n500), **+2.7 → +0.36 → +0.06** (md120), **+4.2 → +0.95 → +0.41** (null), **+4.4 → +0.96 → +0.41** (n700); survival: +2–6 → +0.1–1.6 → −0.4–0.9. The field retains a positive residual (0.4 SD) in the three tie cells and is essentially unbiased at md120; it does not over-correct anywhere on this path. In MD units the naive optimism is 58–76 CD4 cells/mm³ on Ĥ; IJ removes about 72%, the field about 87% (md40), 97% (md120). On Ĥᶜ the retained bias is −1.1 → −0.36 → −0.22 SD (md40); the two-term correction moves the complement the right way and the field halves what is left.

**4. The display: all points on the Gaussian reference — the smallest departures of any path.** Field points: |observed − reference| ≤ **0.008** for the one-sided coverage on both blocks in every cell (H: 0.950 vs 0.944, 0.948 vs 0.949, 0.950 vs 0.944, 0.947 vs 0.939; Ĥᶜ: 0.924 vs 0.926, 0.934 vs 0.928, 0.926 vs 0.928, 0.936 vs 0.931) and ≤ 0.015 for the two-sided. IJ points sit within 0.02 (the largest gap 0.942 vs 0.962, one-sided, n700). The field's r is 1.03–1.21 on Ĥ and 0.96–1.02 on Ĥᶜ: λ-SD tracks the empirical SD of est2 to within 20% (harm) and 5% (complement). No tail, on either block, in any cell.

**5. Regime: every continuous cell is a tie regime, and the field's calibration survives it.** p̂(Ĥ) mean **0.155 / 0.189 / 0.155 / 0.163**, median 0.13–0.17, below 0.5 on **97–99%** of replicates, above 0.9 on none. The family has 1,842 candidates (36 J-quantile and default cuts × two directions, ≤ 2 conjunctions) with many duplicate-membership labels (`str2` ≡ `preanti > 0`; `karnof ≤ 90` ≡ `≤ 95`), and the winner's own re-selection frequency is small in every cell — the tie regime is the design's regime, not a corner. The enumerated ties of the Gate 2 record are ties between **labels with identical membership on the analysis sample**, so p̂ computed on labels understates the settledness of the membership: draws that move the winner between such labels change nothing in the selected subgroup or its estimate, and that is why this tie regime is benign for the one-sided bound. On survival the one tie cell gave 0.920 one-sided; here the tie regime gives 0.947–0.950. A membership-based p̂ (re-selection frequency of the selected *set*, pooling labels with identical membership) would be the more informative diagnostic — a finding only, no task. The joint γ calibrates to Bonferroni in every cell (mean γ 0.0251–0.0252; corr(Λ*, Λ*ᶜ) 0.05–0.09), as on survival.

**6. Joint pair.** Bonferroni joint coverage **0.940 / 0.943 / 0.940 / 0.936** (Wilson ≈ ±0.011), the calibrated pair identical (γ = 0.025 to the grid), the separate 95% pair 0.88 as expected of two independent 95% claims; survival 0.940–0.960 and 0.938–0.942. Rule for a claim on both subgroups: Bonferroni, as before.

**7. Bound location (reading aids, not gates).** md40 cells: the field's lower bound on Ĥ lies at or above 0 on 31% of replicates, above 10 on 19%, above 30 on 5–6%, against a β(Ĥ) of about +32 — the design's harm region is only weakly recovered (PPV 0.40), and the selection-adjusted bounds say so: a harm claim of any size is supported on fewer than a third of the trials, the price of selection on a 1,842-candidate family at n = 500–700 (naive: "harm ≥ 30" on 92–95% of trials, of which most are optimism). md120: lower bound above 0 on 96%, above 30 on 66%, above 40 on 52% against a β(Ĥ) of +96 — the bounds locate the harm at the design's larger effect. Null (truth: no subgroup; homogeneous +26 on the harm-oriented scale): the field's lower bound sits above 0 on 24% of trials — a bound location against +26, not a false-claim rate, since every patient is harmed by 26 and the bound is read against that; above 10 on 14%, above 20 on 8%; the naive lower bound sits above 30 on 87% — the null cell is where the correction matters most and the field's coverage there (0.950) is the cleanest evidence that the construction, not the Cox setting, carries the calibration. Complement: the field's upper bound lies below 30 on 5–12% of md40/null trials and below 20 on 0.6–3%; against β(Ĥᶜ) of +26 to +31 a benefit claim is (correctly) not available, and the IJ upper bound (below 30 on ≤ 0.9%) is uninformative.

**8. Identification.** Detection 0.997–1.000 in every cell including the null (the consistency screen is non-discriminating at these thresholds, the twin's D2 note); mean |Ĥ| 72–75 of 500 (73 of 700); sensitivity 0.17 / 0.33 / – / 0.12, PPV 0.40 / 0.76 / 0.00 / 0.40; the identification-structure tables in each combine document give the anchor/partner/proxy classification.

**9. Cost.** 14.9–17.6 s per replicate at 13 workers (field 11–13 s, complement 0.4 s); 8,000 replicates in 172.8 min cumulative wall; peak memory 17.5 GB.

## What this establishes

The field's one-sided calibration, the complement block's one-sided-upper product and the joint pair's Bonferroni rule behave on the mean-difference scale exactly as on the Cox scale, with the smallest Gaussian-reference departures of any path so far, in a design where every cell is a tie regime and the null cell carries a homogeneous treatment effect: they are properties of the construction. Nothing blocks; no task is proposed.
