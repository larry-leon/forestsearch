# REPORT — Field block for the complement Ĥᶜ: Stage 3 (final)

**Task:** `dev/tasks/TASK_mr_field_complement_2026-09-06.md` (2c519337). Records: Stage 0 92d9ed20 · Stage 1 / Gate 1 PASS ef3e609a, 5e540de2 · Gate 2 record beside this file (`REPORT_mr_field_complement_gate2_2026-09-06.md`). H-C1–H-C4 at defaults; H-C5 unattended pre-authorization.
**Date:** 2026-09-06. Data: campaigns `s7c` (h100 n500, h175 n500) and `map1c` (h150 n500, h150 n1500, h075, h100 n1000, h175 knoise3), 2,000 replicates each, `field_complement = TRUE`, FB none, forestsearch 0.3.5, seeds `8316951 + sim_id` — **the same replicates as the committed `s7`/`map1` bundles, every non-complement-field column identical (Gate 2)**. Rendered documents: the seven `..._{s7c,map1c}_combine_1_2000.html` (each with the complement rows in its tables, the complement diagnostics section and the complement display) and the cross-cell `summary_mr_field_complement.html` (+ `.qmd`).

**Conventions.** Estimates and bounds on the HR scale; SD/SE/bias on the log-HR (β) scale; the target is the conditional β(Ĥᶜ) (θ† and θ‡ coverage reported, not scored). MR (field) for the complement = `est2` (β̃ᶜ − λ̄ᶜ) with the two-sided Λᶜ-quantile interval; its **primary product is the one-sided 95% upper bound `up1s` = exp(β̃ᶜ − q₀.₀₅(Λ*ᶜ))** — the benefit-claim orientation. One-sided coverage in the complement block is `β(Ĥᶜ) ≤ U` for every estimator (naive/oracle/IJ: `exp(log est + 1.645·SE)`). "1s upper margin" = mean of `log U − log(centre)` (centre = β̃ᶜ for the field, which inverts about it; the estimator's own point estimate otherwise); "log half-width" = mean of `(log hi − log lo)/2`. Wilson 95% intervals on coverage. Cell 7 (knoise3) scores β(Ĥᶜ) on the 1,311 of 1,945 detections where the rule evaluates on the super-population (`betaHhat_status = "ok"`; the 634 `unresolved` are the same replicates as in the committed `map1` bundle).

## Per-cell complement block (naive / oracle / MR (IJ) / MR (field))

### Cell h100 n=500 (campaign s7c): 2000 replicates, 1361 detections (68.0%); complement mean |Ĥᶜ| = 421; targets θ† = 0.657, θ‡ = 0.585

| Estimator | n | mean est | bias vs β(Ĥᶜ) log-HR | bias SD units | SD_β | mean SE_β | SE/SD | Cov2s β(Ĥᶜ) (Wilson) | Cov1s upper β(Ĥᶜ) (Wilson) | Cov2s θ† | Cov2s θ‡ | Cov1s-up θ† | Cov1s-up θ‡ | log half-width | 1s upper margin |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1361 | 0.609 | -0.0904 | -0.69 | 0.131 | 0.139 | 1.06 | 0.911 (0.895, 0.925) | 0.859 (0.839, 0.876) | 0.919 | 0.960 | 0.866 | 0.977 | 0.273 | 0.229 |
| oracle | 1361 | 0.658 | -0.0133 | -0.10 | 0.135 | 0.139 | 1.03 | 0.965 (0.954, 0.974) | 0.949 (0.936, 0.960) | 0.958 | 0.888 | 0.942 | 0.993 | 0.272 | 0.228 |
| MR (IJ) | 1361 | 0.653 | -0.0214 | -0.16 | 0.137 | 0.255 | 1.85 | 1.000 (0.997, 1.000) | 0.997 (0.992, 0.999) | 1.000 | 0.999 | 0.997 | 1.000 | 0.499 | 0.419 |
| MR (field) | 1361 | 0.661 | -0.0096 | -0.07 | 0.141 | 0.136 | 0.97 | 0.950 (0.937, 0.960) | 0.930 (0.915, 0.943) | 0.945 | 0.851 | 0.935 | 0.988 | 0.265 | 0.236 |

Display (block Hc, side upper): naive b=-0.688 r=1.060 cov1=0.859 (ref 0.855) cov2=0.911 (ref 0.915); mr b=-0.156 r=1.854 cov1=0.997 (ref 0.998) cov2=1.000 (ref 1.000); fld b=-0.069 r=0.966 cov1=0.930 (ref 0.936) cov2=0.950 (ref 0.941)
Complement λ-SD/SD = 0.97; IJ SE/SD = 1.85; est2 bias − β̃ᶜ bias = +0.0117 (MCSE of a mean log error ≈ 0.0037). Fits: n_complement_fits mean 414.7 (max 572), share_newfit mean 0.048 (max 0.579); complement secs mean 2.90 (q90 3.89) vs harm field 14.6; fit+MR 41.3 s/rep; NA est2 0, notes 0; n_out mean 976 (min 697).

### Cell h175 n=500 (campaign s7c): 2000 replicates, 1900 detections (95.0%); complement mean |Ĥᶜ| = 428; targets θ† = 0.657, θ‡ = 0.585

| Estimator | n | mean est | bias vs β(Ĥᶜ) log-HR | bias SD units | SD_β | mean SE_β | SE/SD | Cov2s β(Ĥᶜ) (Wilson) | Cov1s upper β(Ĥᶜ) (Wilson) | Cov2s θ† | Cov2s θ‡ | Cov1s-up θ† | Cov1s-up θ‡ | log half-width | 1s upper margin |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1900 | 0.628 | -0.0703 | -0.51 | 0.138 | 0.139 | 1.01 | 0.921 (0.907, 0.932) | 0.876 (0.861, 0.890) | 0.935 | 0.927 | 0.895 | 0.983 | 0.272 | 0.228 |
| oracle | 1900 | 0.638 | -0.0559 | -0.40 | 0.140 | 0.139 | 0.99 | 0.939 (0.928, 0.949) | 0.900 (0.886, 0.913) | 0.943 | 0.913 | 0.911 | 0.984 | 0.273 | 0.229 |
| MR (IJ) | 1900 | 0.669 | -0.0082 | -0.06 | 0.143 | 0.257 | 1.80 | 1.000 (0.998, 1.000) | 0.999 (0.996, 1.000) | 1.000 | 0.999 | 0.999 | 1.000 | 0.504 | 0.423 |
| MR (field) | 1900 | 0.675 | -0.0006 | -0.00 | 0.147 | 0.136 | 0.93 | 0.934 (0.922, 0.944) | 0.936 (0.924, 0.946) | 0.927 | 0.795 | 0.945 | 0.991 | 0.266 | 0.232 |

Display (block Hc, side upper): naive b=-0.511 r=1.007 cov1=0.876 (ref 0.874) cov2=0.921 (ref 0.922); mr b=-0.058 r=1.800 cov1=0.999 (ref 0.998) cov2=1.000 (ref 1.000); fld b=-0.004 r=0.927 cov1=0.936 (ref 0.936) cov2=0.934 (ref 0.931)
Complement λ-SD/SD = 0.93; IJ SE/SD = 1.80; est2 bias − β̃ᶜ bias = +0.0076 (MCSE of a mean log error ≈ 0.0033). Fits: n_complement_fits mean 382.0 (max 539), share_newfit mean 0.104 (max 0.663); complement secs mean 3.01 (q90 3.85) vs harm field 16.0; fit+MR 45.9 s/rep; NA est2 0, notes 0; n_out mean 993 (min 766).

### Cell h150 n=500 (campaign map1c): 2000 replicates, 1822 detections (91.1%); complement mean |Ĥᶜ| = 427; targets θ† = 0.657, θ‡ = 0.585

| Estimator | n | mean est | bias vs β(Ĥᶜ) log-HR | bias SD units | SD_β | mean SE_β | SE/SD | Cov2s β(Ĥᶜ) (Wilson) | Cov1s upper β(Ĥᶜ) (Wilson) | Cov2s θ† | Cov2s θ‡ | Cov1s-up θ† | Cov1s-up θ‡ | log half-width | 1s upper margin |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1822 | 0.624 | -0.0806 | -0.59 | 0.136 | 0.139 | 1.02 | 0.915 (0.902, 0.927) | 0.864 (0.848, 0.879) | 0.935 | 0.937 | 0.886 | 0.981 | 0.272 | 0.228 |
| oracle | 1822 | 0.640 | -0.0555 | -0.40 | 0.139 | 0.139 | 1.00 | 0.943 (0.932, 0.953) | 0.900 (0.885, 0.913) | 0.946 | 0.911 | 0.917 | 0.987 | 0.273 | 0.229 |
| MR (IJ) | 1822 | 0.667 | -0.0148 | -0.10 | 0.141 | 0.256 | 1.81 | 1.000 (0.998, 1.000) | 0.997 (0.993, 0.998) | 1.000 | 0.997 | 0.999 | 1.000 | 0.502 | 0.421 |
| MR (field) | 1822 | 0.674 | -0.0052 | -0.04 | 0.145 | 0.136 | 0.94 | 0.935 (0.923, 0.946) | 0.936 (0.924, 0.947) | 0.935 | 0.802 | 0.947 | 0.993 | 0.265 | 0.233 |

Display (block Hc, side upper): naive b=-0.592 r=1.018 cov1=0.864 (ref 0.861) cov2=0.915 (ref 0.915); mr b=-0.105 r=1.813 cov1=0.997 (ref 0.998) cov2=1.000 (ref 1.000); fld b=-0.036 r=0.938 cov1=0.936 (ref 0.934) cov2=0.935 (ref 0.934)
Complement λ-SD/SD = 0.94; IJ SE/SD = 1.81; est2 bias − β̃ᶜ bias = +0.0096 (MCSE of a mean log error ≈ 0.0033). Fits: n_complement_fits mean 394.6 (max 565), share_newfit mean 0.080 (max 0.803); complement secs mean 3.04 (q90 3.89) vs harm field 15.7; fit+MR 44.8 s/rep; NA est2 0, notes 0; n_out mean 989 (min 817).

### Cell h150 n=1500 (campaign map1c): 2000 replicates, 1976 detections (98.8%); complement mean |Ĥᶜ| = 1315; targets θ† = 0.657, θ‡ = 0.585

| Estimator | n | mean est | bias vs β(Ĥᶜ) log-HR | bias SD units | SD_β | mean SE_β | SE/SD | Cov2s β(Ĥᶜ) (Wilson) | Cov1s upper β(Ĥᶜ) (Wilson) | Cov2s θ† | Cov2s θ‡ | Cov1s-up θ† | Cov1s-up θ‡ | log half-width | 1s upper margin |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1976 | 0.635 | -0.0196 | -0.25 | 0.079 | 0.079 | 1.01 | 0.953 (0.943, 0.962) | 0.923 (0.910, 0.934) | 0.935 | 0.831 | 0.891 | 0.996 | 0.156 | 0.131 |
| oracle | 1976 | 0.630 | -0.0264 | -0.34 | 0.077 | 0.080 | 1.04 | 0.937 (0.925, 0.947) | 0.911 (0.898, 0.923) | 0.927 | 0.861 | 0.876 | 0.995 | 0.157 | 0.132 |
| MR (IJ) | 1976 | 0.652 | +0.0066 | +0.08 | 0.082 | 0.151 | 1.84 | 0.999 (0.996, 1.000) | 0.999 (0.996, 1.000) | 0.999 | 0.988 | 0.997 | 1.000 | 0.296 | 0.248 |
| MR (field) | 1976 | 0.652 | +0.0068 | +0.08 | 0.084 | 0.079 | 0.94 | 0.946 (0.936, 0.955) | 0.956 (0.947, 0.965) | 0.937 | 0.707 | 0.929 | 0.997 | 0.155 | 0.130 |

Display (block Hc, side upper): naive b=-0.249 r=1.011 cov1=0.923 (ref 0.921) cov2=0.953 (ref 0.945); mr b=+0.081 r=1.837 cov1=0.999 (ref 0.999) cov2=0.999 (ref 1.000); fld b=+0.081 r=0.941 cov1=0.956 (ref 0.948) cov2=0.946 (ref 0.934)
Complement λ-SD/SD = 0.94; IJ SE/SD = 1.84; est2 bias − β̃ᶜ bias = +0.0001 (MCSE of a mean log error ≈ 0.0018). Fits: n_complement_fits mean 347.1 (max 572), share_newfit mean 0.202 (max 0.932); complement secs mean 14.14 (q90 19.56) vs harm field 35.9; fit+MR 161.6 s/rep; NA est2 0, notes 0; n_out mean 994 (min 747).

### Cell h075 n=500 (campaign map1c): 2000 replicates, 1042 detections (52.1%); complement mean |Ĥᶜ| = 419; targets θ† = 0.657, θ‡ = 0.585

| Estimator | n | mean est | bias vs β(Ĥᶜ) log-HR | bias SD units | SD_β | mean SE_β | SE/SD | Cov2s β(Ĥᶜ) (Wilson) | Cov1s upper β(Ĥᶜ) (Wilson) | Cov2s θ† | Cov2s θ‡ | Cov1s-up θ† | Cov1s-up θ‡ | log half-width | 1s upper margin |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1042 | 0.594 | -0.0848 | -0.65 | 0.131 | 0.140 | 1.07 | 0.923 (0.905, 0.938) | 0.872 (0.851, 0.891) | 0.893 | 0.968 | 0.834 | 0.964 | 0.275 | 0.230 |
| oracle | 1042 | 0.674 | +0.0410 | +0.31 | 0.130 | 0.139 | 1.06 | 0.960 (0.946, 0.970) | 0.979 (0.968, 0.986) | 0.964 | 0.858 | 0.962 | 0.996 | 0.272 | 0.228 |
| MR (IJ) | 1042 | 0.636 | -0.0170 | -0.12 | 0.138 | 0.255 | 1.84 | 1.000 (0.996, 1.000) | 0.995 (0.989, 0.998) | 0.998 | 0.999 | 0.994 | 1.000 | 0.499 | 0.419 |
| MR (field) | 1042 | 0.644 | -0.0056 | -0.04 | 0.142 | 0.137 | 0.96 | 0.943 (0.928, 0.956) | 0.936 (0.919, 0.949) | 0.936 | 0.883 | 0.910 | 0.982 | 0.266 | 0.236 |

Display (block Hc, side upper): naive b=-0.649 r=1.071 cov1=0.872 (ref 0.867) cov2=0.923 (ref 0.923); mr b=-0.123 r=1.843 cov1=0.995 (ref 0.998) cov2=1.000 (ref 1.000); fld b=-0.039 r=0.960 cov1=0.936 (ref 0.938) cov2=0.943 (ref 0.940)
Complement λ-SD/SD = 0.96; IJ SE/SD = 1.84; est2 bias − β̃ᶜ bias = +0.0114 (MCSE of a mean log error ≈ 0.0043). Fits: n_complement_fits mean 419.9 (max 562), share_newfit mean 0.042 (max 0.382); complement secs mean 2.81 (q90 3.79) vs harm field 13.9; fit+MR 39.3 s/rep; NA est2 0, notes 0; n_out mean 966 (min 653).

### Cell h100 n=1000 (campaign map1c): 2000 replicates, 1319 detections (66.0%); complement mean |Ĥᶜ| = 871; targets θ† = 0.657, θ‡ = 0.585

| Estimator | n | mean est | bias vs β(Ĥᶜ) log-HR | bias SD units | SD_β | mean SE_β | SE/SD | Cov2s β(Ĥᶜ) (Wilson) | Cov1s upper β(Ĥᶜ) (Wilson) | Cov2s θ† | Cov2s θ‡ | Cov1s-up θ† | Cov1s-up θ‡ | log half-width | 1s upper margin |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1319 | 0.631 | -0.0441 | -0.48 | 0.092 | 0.097 | 1.05 | 0.943 (0.929, 0.954) | 0.908 (0.891, 0.922) | 0.941 | 0.892 | 0.894 | 0.991 | 0.189 | 0.159 |
| oracle | 1319 | 0.648 | -0.0172 | -0.18 | 0.095 | 0.098 | 1.03 | 0.966 (0.955, 0.974) | 0.945 (0.931, 0.956) | 0.951 | 0.836 | 0.928 | 0.995 | 0.192 | 0.161 |
| MR (IJ) | 1319 | 0.660 | +0.0004 | +0.00 | 0.097 | 0.179 | 1.85 | 1.000 (0.997, 1.000) | 0.998 (0.993, 0.999) | 1.000 | 0.995 | 0.998 | 1.000 | 0.351 | 0.295 |
| MR (field) | 1319 | 0.665 | +0.0076 | +0.08 | 0.099 | 0.095 | 0.96 | 0.945 (0.932, 0.956) | 0.953 (0.940, 0.963) | 0.940 | 0.711 | 0.948 | 0.997 | 0.185 | 0.164 |

Display (block Hc, side upper): naive b=-0.479 r=1.050 cov1=0.908 (ref 0.894) cov2=0.943 (ref 0.937); mr b=+0.004 r=1.854 cov1=0.998 (ref 0.999) cov2=1.000 (ref 1.000); fld b=+0.077 r=0.958 cov1=0.953 (ref 0.951) cov2=0.945 (ref 0.939)
Complement λ-SD/SD = 0.96; IJ SE/SD = 1.85; est2 bias − β̃ᶜ bias = +0.0072 (MCSE of a mean log error ≈ 0.0027). Fits: n_complement_fits mean 416.8 (max 555), share_newfit mean 0.069 (max 0.692); complement secs mean 6.73 (q90 10.04) vs harm field 20.7; fit+MR 79.9 s/rep; NA est2 0, notes 0; n_out mean 972 (min 626).

### Cell h175 knoise3 n=500 (campaign map1c): 2000 replicates, 1945 detections (97.2%); complement mean |Ĥᶜ| = 429; targets θ† = 0.657, θ‡ = 0.585

| Estimator | n | mean est | bias vs β(Ĥᶜ) log-HR | bias SD units | SD_β | mean SE_β | SE/SD | Cov2s β(Ĥᶜ) (Wilson) | Cov1s upper β(Ĥᶜ) (Wilson) | Cov2s θ† | Cov2s θ‡ | Cov1s-up θ† | Cov1s-up θ‡ | log half-width | 1s upper margin |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1311 | 0.625 | -0.0823 | -0.60 | 0.137 | 0.138 | 1.01 | 0.899 (0.882, 0.914) | 0.853 (0.833, 0.871) | 0.919 | 0.933 | 0.864 | 0.977 | 0.270 | 0.227 |
| oracle | 1311 | 0.637 | -0.0691 | -0.49 | 0.140 | 0.139 | 0.99 | 0.926 (0.911, 0.939) | 0.873 (0.853, 0.890) | 0.930 | 0.926 | 0.884 | 0.979 | 0.273 | 0.229 |
| MR (IJ) | 1311 | 0.672 | -0.0110 | -0.08 | 0.140 | 0.256 | 1.83 | 1.000 (0.997, 1.000) | 0.998 (0.994, 1.000) | 1.000 | 0.998 | 0.999 | 1.000 | 0.502 | 0.422 |
| MR (field) | 1311 | 0.680 | -0.0017 | -0.01 | 0.143 | 0.135 | 0.94 | 0.937 (0.922, 0.949) | 0.937 (0.923, 0.949) | 0.933 | 0.818 | 0.946 | 0.993 | 0.263 | 0.233 |

Display (block Hc, side upper): naive b=-0.601 r=1.008 cov1=0.853 (ref 0.855) cov2=0.899 (ref 0.910); mr b=-0.079 r=1.833 cov1=0.998 (ref 0.998) cov2=1.000 (ref 1.000); fld b=-0.012 r=0.945 cov1=0.937 (ref 0.938) cov2=0.937 (ref 0.936)
Complement λ-SD/SD = 0.94; IJ SE/SD = 1.83; est2 bias − β̃ᶜ bias = +0.0093 (MCSE of a mean log error ≈ 0.0032). Fits: n_complement_fits mean 603.2 (max 839), share_newfit mean 0.092 (max 0.915); complement secs mean 6.45 (q90 8.69) vs harm field 24.7; fit+MR 90.9 s/rep; NA est2 0, notes 0; n_out mean 998 (min 874).

## Across cells

### The complement display (block Ĥᶜ, side upper) — `summary_mr_field_complement.html`

| Cell | Estimator | b | r | 1s upper cov (Wilson) | ref Φ(1.645r+b) | 2s cov (Wilson) | ref |
|---|---|---|---|---|---|---|---|
| null 1.0, n=500 | IJ | −0.156 | 1.854 | 0.997 (0.992, 0.999) | 0.998 | 1.000 | 1.000 |
| null 1.0, n=500 | field | −0.069 | 0.966 | 0.930 (0.915, 0.943) | 0.936 | 0.950 (0.937, 0.960) | 0.941 |
| harm 1.75, n=500 | IJ | −0.058 | 1.800 | 0.999 | 0.998 | 1.000 | 1.000 |
| harm 1.75, n=500 | field | −0.004 | 0.927 | 0.936 (0.924, 0.946) | 0.936 | 0.934 (0.922, 0.944) | 0.931 |
| harm 1.5, n=500 | IJ | −0.105 | 1.813 | 0.997 | 0.998 | 1.000 | 1.000 |
| harm 1.5, n=500 | field | −0.036 | 0.938 | 0.936 (0.924, 0.947) | 0.934 | 0.935 (0.923, 0.946) | 0.934 |
| harm 1.5, n=1500 | IJ | +0.081 | 1.837 | 0.999 | 0.999 | 0.999 | 1.000 |
| harm 1.5, n=1500 | field | +0.081 | 0.941 | 0.956 (0.947, 0.965) | 0.948 | 0.946 (0.936, 0.955) | 0.934 |
| protective 0.75 | IJ | −0.123 | 1.843 | 0.995 | 0.998 | 1.000 | 1.000 |
| protective 0.75 | field | −0.039 | 0.960 | 0.936 (0.919, 0.949) | 0.938 | 0.943 (0.928, 0.956) | 0.940 |
| null 1.0, n=1000 | IJ | +0.004 | 1.854 | 0.998 | 0.999 | 1.000 | 1.000 |
| null 1.0, n=1000 | field | +0.077 | 0.958 | 0.953 (0.940, 0.963) | 0.951 | 0.945 (0.932, 0.956) | 0.939 |
| harm 1.75, 3 noise | IJ | −0.079 | 1.833 | 0.998 | 0.998 | 1.000 | 1.000 |
| harm 1.75, 3 noise | field | −0.012 | 0.945 | 0.937 (0.923, 0.949) | 0.938 | 0.937 (0.922, 0.949) | 0.936 |

Every one of the 14 points sits on its Gaussian reference (max |gap| 0.008 one-sided, 0.012 two-sided): in the complement block the coverage patterns are entirely the arithmetic of retained bias and SE calibration — there is no over-correction tail here, unlike the harm block's two-sided points.

### Reading criteria (Larry's, not gates)

| Criterion | Result across the seven cells |
|---|---|
| Field one-sided **upper** coverage of β(Ĥᶜ) within Monte Carlo error of nominal | 0.930, 0.936, 0.936, 0.956, 0.936, 0.953, 0.937 — Wilson intervals include 0.95 in 2 of 7 (n=1500, n=1000); the five n = 500 cells sit at 0.930–0.937 with upper Wilson limits 0.943–0.949, i.e. **1.3–2.0 points under nominal, just outside MC error**, tracking λ-SD/SD = 0.93–0.97 (the field's SE runs 3–7% under the empirical SD at n = 500) plus a residual b of −0.01 to −0.07 SD. At n ≥ 1000 (r = 0.94–0.96, b ≈ +0.08) coverage is nominal. |
| λ-SD/SD in [0.9, 1.2] | 0.97, 0.93, 0.94, 0.94, 0.96, 0.96, 0.94 — **all inside**, at the low end. |
| Margin materially below IJ's (expected roughly half) | Field one-sided upper margin 0.232–0.236 (n=500), 0.164 (n=1000), 0.130 (n=1500) vs IJ 0.419–0.423 / 0.295 / 0.248: **0.53–0.56 of IJ's** (n=500), 0.56 (n=1000), 0.52 (n=1500). Two-sided log half-width likewise 0.263–0.266 vs 0.499–0.504. IJ's r = 1.80–1.85 is the SE/SD 1.8–1.9 the task named; its coverage of 0.997–0.999 is the price. |
| est₂ᶜ bias not worse than β̃ᶜ's beyond Monte Carlo error | est₂ᶜ − β̃ᶜ retained bias = +0.012, +0.008, +0.010, +0.000, +0.011, +0.007, +0.009 log-HR (MCSE of a mean 0.002–0.004): est₂ᶜ is **closer to β(Ĥᶜ)** in six cells (|bias| ≤ 0.010 vs IJ's 0.008–0.021) and equal at n = 1500; its bias sits at −0.07 to +0.08 SD everywhere. The complement's λ̄ᶜ is small and *negative* (−0.005 to −0.011 log-HR): the harm field's second-order correction on the complement side moves the estimate slightly *away from benefit*, the right direction. |

### Bound locations (the Reading's inputs)

| Cell | mean β(Ĥᶜ) | naive est / 1s-upper | IJ β̃ᶜ / 1s-upper | field est₂ᶜ / 1s-upper | share of field U < 0.85 / < 0.80 | share of IJ U < 0.85 / < 0.80 |
|---|---|---|---|---|---|---|
| null 1.0, n=500 | 0.661 | 0.609 / 0.765 | 0.653 / 0.992 | 0.661 / 0.827 | 0.59 / 0.41 | 0.15 / 0.08 |
| harm 1.75, n=500 | 0.669 | 0.628 / 0.789 | 0.669 / 1.021 | 0.675 / 0.844 | 0.54 / 0.37 | 0.12 / 0.05 |
| harm 1.5, n=500 | 0.671 | 0.624 / 0.784 | 0.667 / 1.016 | 0.674 / 0.843 | 0.54 / 0.37 | 0.13 / 0.05 |
| harm 1.5, n=1500 | 0.646 | 0.635 / 0.723 | 0.652 / 0.835 | 0.652 / 0.743 | 0.95 / 0.83 | 0.60 / 0.31 |
| protective 0.75 | 0.641 | 0.594 / 0.748 | 0.636 / 0.967 | 0.644 / 0.806 | 0.67 / 0.49 | 0.19 / 0.11 |
| null 1.0, n=1000 | 0.657 | 0.631 / 0.739 | 0.660 / 0.886 | 0.665 / 0.777 | 0.83 / 0.63 | 0.34 / 0.15 |
| harm 1.75, 3 noise | 0.664 | 0.625 / 0.783 | 0.672 / 1.024 | 0.680 / 0.848 | 0.52 / 0.36 | 0.11 / 0.04 |

## Reading

Read against effect size rather than against HR = 1.0, with 0.80 / 0.85 named as reading aids for a clinically meaningful benefit, not as decision rules. In every cell the complement's true conditional effect is a substantial benefit, β(Ĥᶜ) ≈ 0.64–0.67, and the naive estimate overstates it by 0.02–0.09 log-HR (0.61–0.63 at n = 500) — the beneficial optimism the task named. The naive one-sided upper bound sits at 0.74–0.79 and covers β(Ĥᶜ) only 85–92% of the time: the benefit it would claim is too strong too often. The two-term IJ estimate β̃ᶜ removes that optimism (bias −0.02 to +0.01 log-HR) but its upper bound sits at **0.97–1.02 at n = 500** (0.89 at n = 1000, 0.84 at n = 1500): with the IJ SE at 1.8× the empirical SD, the adjusted interval can rule out essentially nothing at n = 500 — only 12–15% of replicates place U below 0.85 — although the underlying benefit is 0.65. The complement field keeps the same centre (est₂ᶜ within 0.01 of β̃ᶜ, bias ≤ 0.08 SD) and places its upper bound at **0.83–0.85 at n = 500, 0.78 at n = 1000, 0.74 at n = 1500**: it locates the benefit the data can support — a true HR of 0.65 read as "at most ≈ 0.84" — with U below 0.85 in 52–59% of n = 500 replicates and 83–95% at n ≥ 1000, while covering β(Ĥᶜ) 93–96% of the time. The spread between the adjusted upper bounds — IJ ≈ 1.00 vs field ≈ 0.84 at n = 500, 0.19 log-HR — is what the two constructions pay for selection: the IJ pays it in width (a 0.42 one-sided margin, coverage 0.997–0.999), the field pays a residual 1–2 points of one-sided coverage at n = 500 (0.930–0.937, λ-SD/SD 0.93–0.97) that closes at n ≥ 1000. Unlike the harm block, the complement's coverages are fully explained by the Gaussian (b, r) bookkeeping — the display's 14 points all lie on their reference curves — so the small n = 500 shortfall is an SE-calibration effect of a few percent, not a tail phenomenon. The complement fits cost 2.8–3.0 s per replicate at n = 500 (≈ 7% of the gate) and 14 s at n = 1,500; every draw winner's complement was fit in every cell.

No task proposed; nothing blocked. Findings in the record. (H-C6 — the complement rows in the GBSG frozen-family interval table — is queued as its own task.)
