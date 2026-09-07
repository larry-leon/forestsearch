# REPORT — Complement inference refinements: Stage 3 (final)

**Task:** `dev/tasks/TASK_complement_refinements_2026-09-06.md` (95e6c8c5). Records: Stage 0 291df669 · Stage 1 / Gate 1 PASS 8fbb3bdc · Gate 2 record beside this file (`REPORT_complement_refinements_gate2_2026-09-06.md`, including the cell-7 campaign-tag incident and its correction). K-1–K-3 at defaults; K-4 unattended pre-authorization.
**Date:** 2026-09-06. Data: campaigns `s7w` (h100 n500, h175 n500) and `map1w` (h150 n500, h150 n1500, h075, h100 n1000, h175 knoise3), 2,000 replicates each, `ci_method = "field"`, `include_complement = TRUE`, `field_complement = TRUE`, `ij_residual = "two_term"` (reported rows unchanged; variants recorded), FB none — **the same replicates as `s7c`/`map1c`, every existing column identical (Gate 2)**. Rendered documents: the seven `..._{s7w,map1w}_combine_1_2000.html` (both blocks' tables carry the two winner rows; the joint-pair table; the displays with the winner rows) and the cross-cell `summary_complement_refinements.html` (+ `.qmd`).

**Conventions.** Bias/SD/SE on the log-HR (β) scale against the conditional targets β(Ĥ), β(Ĥᶜ); one-sided coverage on each block's exposed side — Ĥ: `β(Ĥ) ≥ L` with `L = exp(log est − 1.645·SE)` (field: stored `lo1s`); Ĥᶜ: `β(Ĥᶜ) ≤ U` with `U = exp(log est + 1.645·SE)` (field: stored `up1s`); two-sided from each estimator's stored bounds; Wilson 95% intervals. The three IJ rows share the two-term point estimate β̃ and differ only in SE: two-term (the paper's), winner-only, winner-floor (= max(winner, naive SE)). "1s margin" = |log bound − log centre| (centre = β̃ for the IJ rows and the field). Joint table: coverage of `β(Ĥ) ≥ lower_H and β(Ĥᶜ) ≤ upper_Hc` for the two separate 95% field bounds, the Bonferroni pair (γ = 0.025) and the calibrated pair; margins on the log scale from β̃ / β̃ᶜ. Cell 7 scores on the 1,311 of 1,945 detections with `betaHhat_status = "ok"` (as in `map1`/`map1c`).

## Per-cell tables

### Cell h100 n=500 (campaign s7w): 2000 replicates, 1361 detections (68.0%)

**Harm Ĥ block** (one-sided = lower bound)

| Estimator | n | bias vs β log-HR | bias SD units | SD_β | mean SE_β | SE/SD | Cov2s β (Wilson) | Cov1s lower β (Wilson) | log half-width | 1s margin |
|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1361 | +0.8134 | +5.00 | 0.163 | 0.324 | 1.99 | 0.157 (0.138, 0.177) | 0.051 (0.040, 0.064) | 0.635 | 0.533 |
| MR (IJ, two-term) | 1361 | +0.2442 | +1.10 | 0.222 | 0.399 | 1.79 | 0.986 (0.978, 0.991) | 0.966 (0.955, 0.975) | 0.782 | 0.656 |
| MR (IJ, winner) | 1361 | +0.2442 | +1.10 | 0.222 | 0.175 | 0.79 | 0.642 (0.616, 0.667) | 0.580 (0.553, 0.606) | 0.343 | 0.288 |
| MR (IJ, winner-floor) | 1361 | +0.2442 | +1.10 | 0.222 | 0.324 | 1.46 | 0.950 (0.937, 0.960) | 0.892 (0.874, 0.907) | 0.635 | 0.533 |
| MR (field) | 1361 | +0.1405 | +0.56 | 0.253 | 0.357 | 1.41 | 0.988 (0.981, 0.993) | 0.963 (0.951, 0.971) | 0.693 | 0.734 |

**Complement Ĥᶜ block** (one-sided = upper bound)

| Estimator | n | bias vs β log-HR | bias SD units | SD_β | mean SE_β | SE/SD | Cov2s β (Wilson) | Cov1s upper β (Wilson) | log half-width | 1s margin |
|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1361 | -0.0904 | -0.69 | 0.131 | 0.139 | 1.06 | 0.911 (0.895, 0.925) | 0.859 (0.839, 0.876) | 0.273 | 0.229 |
| MR (IJ, two-term) | 1361 | -0.0214 | -0.16 | 0.137 | 0.255 | 1.85 | 1.000 (0.997, 1.000) | 0.997 (0.992, 0.999) | 0.499 | 0.419 |
| MR (IJ, winner) | 1361 | -0.0214 | -0.16 | 0.137 | 0.123 | 0.89 | 0.917 (0.901, 0.930) | 0.896 (0.878, 0.911) | 0.240 | 0.202 |
| MR (IJ, winner-floor) | 1361 | -0.0214 | -0.16 | 0.137 | 0.139 | 1.01 | 0.954 (0.941, 0.964) | 0.932 (0.917, 0.944) | 0.273 | 0.229 |
| MR (field) | 1361 | -0.0096 | -0.07 | 0.141 | 0.136 | 0.97 | 0.950 (0.937, 0.960) | 0.930 (0.915, 0.943) | 0.265 | 0.236 |

**Joint (Ĥ lower, Ĥᶜ upper)** — n = 1361; mean γ = 0.0250 (share at 0.025: 0.968); achieved joint prob. mean 0.9499 (Bonferroni 0.9500); mean corr(Λ*, Λ*ᶜ) = +0.047

| Pair | joint cov (Wilson) | cov Ĥ | cov Ĥᶜ | margin Ĥ (log) | margin Ĥᶜ (log) |
|---|---|---|---|---|---|
| separate 95% field bounds | 0.894 (0.877, 0.909) | 0.963 | 0.930 | 0.734 | 0.236 |
| Bonferroni (gamma = 0.025) | 0.955 (0.943, 0.965) | 0.988 | 0.967 | 0.883 | 0.279 |
| calibrated gamma | 0.955 (0.943, 0.965) | 0.988 | 0.967 | 0.883 | 0.279 |

### Cell h175 n=500 (campaign s7w): 2000 replicates, 1900 detections (95.0%)

**Harm Ĥ block** (one-sided = lower bound)

| Estimator | n | bias vs β log-HR | bias SD units | SD_β | mean SE_β | SE/SD | Cov2s β (Wilson) | Cov1s lower β (Wilson) | log half-width | 1s margin |
|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1900 | +0.4809 | +2.17 | 0.222 | 0.310 | 1.40 | 0.645 (0.623, 0.666) | 0.539 (0.516, 0.561) | 0.607 | 0.510 |
| MR (IJ, two-term) | 1900 | -0.0101 | -0.03 | 0.312 | 0.393 | 1.26 | 0.961 (0.951, 0.969) | 0.983 (0.976, 0.988) | 0.770 | 0.646 |
| MR (IJ, winner) | 1900 | -0.0101 | -0.03 | 0.312 | 0.179 | 0.57 | 0.607 (0.585, 0.629) | 0.767 (0.748, 0.786) | 0.351 | 0.295 |
| MR (IJ, winner-floor) | 1900 | -0.0101 | -0.03 | 0.312 | 0.310 | 0.99 | 0.889 (0.874, 0.902) | 0.929 (0.917, 0.940) | 0.607 | 0.510 |
| MR (field) | 1900 | -0.0748 | -0.21 | 0.358 | 0.363 | 1.01 | 0.887 (0.872, 0.901) | 0.966 (0.957, 0.974) | 0.707 | 0.699 |

**Complement Ĥᶜ block** (one-sided = upper bound)

| Estimator | n | bias vs β log-HR | bias SD units | SD_β | mean SE_β | SE/SD | Cov2s β (Wilson) | Cov1s upper β (Wilson) | log half-width | 1s margin |
|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1900 | -0.0703 | -0.51 | 0.138 | 0.139 | 1.01 | 0.921 (0.907, 0.932) | 0.876 (0.861, 0.890) | 0.272 | 0.228 |
| MR (IJ, two-term) | 1900 | -0.0082 | -0.06 | 0.143 | 0.257 | 1.80 | 1.000 (0.998, 1.000) | 0.999 (0.996, 1.000) | 0.504 | 0.423 |
| MR (IJ, winner) | 1900 | -0.0082 | -0.06 | 0.143 | 0.124 | 0.87 | 0.914 (0.901, 0.926) | 0.913 (0.899, 0.925) | 0.243 | 0.204 |
| MR (IJ, winner-floor) | 1900 | -0.0082 | -0.06 | 0.143 | 0.139 | 0.97 | 0.947 (0.936, 0.957) | 0.938 (0.927, 0.948) | 0.272 | 0.228 |
| MR (field) | 1900 | -0.0006 | -0.00 | 0.147 | 0.136 | 0.93 | 0.934 (0.922, 0.944) | 0.936 (0.924, 0.946) | 0.266 | 0.232 |

**Joint (Ĥ lower, Ĥᶜ upper)** — n = 1900; mean γ = 0.0251 (share at 0.025: 0.924); achieved joint prob. mean 0.9502 (Bonferroni 0.9503); mean corr(Λ*, Λ*ᶜ) = +0.014

| Pair | joint cov (Wilson) | cov Ĥ | cov Ĥᶜ | margin Ĥ (log) | margin Ĥᶜ (log) |
|---|---|---|---|---|---|
| separate 95% field bounds | 0.905 (0.891, 0.918) | 0.966 | 0.936 | 0.699 | 0.232 |
| Bonferroni (gamma = 0.025) | 0.956 (0.946, 0.964) | 0.990 | 0.965 | 0.847 | 0.275 |
| calibrated gamma | 0.956 (0.946, 0.964) | 0.990 | 0.965 | 0.846 | 0.275 |

### Cell h150 n=500 (campaign map1w): 2000 replicates, 1822 detections (91.1%)

**Harm Ĥ block** (one-sided = lower bound)

| Estimator | n | bias vs β log-HR | bias SD units | SD_β | mean SE_β | SE/SD | Cov2s β (Wilson) | Cov1s lower β (Wilson) | log half-width | 1s margin |
|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1822 | +0.5802 | +2.92 | 0.199 | 0.315 | 1.59 | 0.529 (0.506, 0.552) | 0.417 (0.394, 0.439) | 0.618 | 0.518 |
| MR (IJ, two-term) | 1822 | +0.0586 | +0.21 | 0.276 | 0.395 | 1.43 | 0.977 (0.969, 0.983) | 0.978 (0.970, 0.984) | 0.773 | 0.649 |
| MR (IJ, winner) | 1822 | +0.0586 | +0.21 | 0.276 | 0.177 | 0.64 | 0.626 (0.604, 0.648) | 0.727 (0.706, 0.747) | 0.348 | 0.292 |
| MR (IJ, winner-floor) | 1822 | +0.0586 | +0.21 | 0.276 | 0.315 | 1.14 | 0.923 (0.910, 0.935) | 0.926 (0.914, 0.938) | 0.618 | 0.518 |
| MR (field) | 1822 | -0.0209 | -0.07 | 0.316 | 0.359 | 1.14 | 0.919 (0.905, 0.930) | 0.965 (0.956, 0.973) | 0.699 | 0.710 |

**Complement Ĥᶜ block** (one-sided = upper bound)

| Estimator | n | bias vs β log-HR | bias SD units | SD_β | mean SE_β | SE/SD | Cov2s β (Wilson) | Cov1s upper β (Wilson) | log half-width | 1s margin |
|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1822 | -0.0806 | -0.59 | 0.136 | 0.139 | 1.02 | 0.915 (0.902, 0.927) | 0.864 (0.848, 0.879) | 0.272 | 0.228 |
| MR (IJ, two-term) | 1822 | -0.0148 | -0.10 | 0.141 | 0.256 | 1.81 | 1.000 (0.998, 1.000) | 0.997 (0.993, 0.998) | 0.502 | 0.421 |
| MR (IJ, winner) | 1822 | -0.0148 | -0.10 | 0.141 | 0.123 | 0.87 | 0.917 (0.903, 0.928) | 0.905 (0.891, 0.918) | 0.242 | 0.203 |
| MR (IJ, winner-floor) | 1822 | -0.0148 | -0.10 | 0.141 | 0.139 | 0.98 | 0.951 (0.940, 0.960) | 0.938 (0.926, 0.948) | 0.272 | 0.228 |
| MR (field) | 1822 | -0.0052 | -0.04 | 0.145 | 0.136 | 0.94 | 0.935 (0.923, 0.946) | 0.936 (0.924, 0.947) | 0.265 | 0.233 |

**Joint (Ĥ lower, Ĥᶜ upper)** — n = 1822; mean γ = 0.0251 (share at 0.025: 0.943); achieved joint prob. mean 0.9501 (Bonferroni 0.9502); mean corr(Λ*, Λ*ᶜ) = +0.030

| Pair | joint cov (Wilson) | cov Ĥ | cov Ĥᶜ | margin Ĥ (log) | margin Ĥᶜ (log) |
|---|---|---|---|---|---|
| separate 95% field bounds | 0.905 (0.891, 0.918) | 0.965 | 0.936 | 0.710 | 0.233 |
| Bonferroni (gamma = 0.025) | 0.953 (0.942, 0.962) | 0.990 | 0.963 | 0.858 | 0.276 |
| calibrated gamma | 0.953 (0.942, 0.962) | 0.990 | 0.963 | 0.858 | 0.276 |

### Cell h150 n=1500 (campaign map1w): 2000 replicates, 1976 detections (98.8%)

**Harm Ĥ block** (one-sided = lower bound)

| Estimator | n | bias vs β log-HR | bias SD units | SD_β | mean SE_β | SE/SD | Cov2s β (Wilson) | Cov1s lower β (Wilson) | log half-width | 1s margin |
|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1976 | +0.1447 | +0.91 | 0.158 | 0.183 | 1.16 | 0.866 (0.851, 0.881) | 0.808 (0.790, 0.825) | 0.359 | 0.302 |
| MR (IJ, two-term) | 1976 | -0.0842 | -0.37 | 0.227 | 0.246 | 1.08 | 0.901 (0.887, 0.913) | 0.981 (0.974, 0.986) | 0.482 | 0.404 |
| MR (IJ, winner) | 1976 | -0.0842 | -0.37 | 0.227 | 0.110 | 0.49 | 0.573 (0.551, 0.595) | 0.863 (0.848, 0.878) | 0.217 | 0.182 |
| MR (IJ, winner-floor) | 1976 | -0.0842 | -0.37 | 0.227 | 0.183 | 0.81 | 0.810 (0.792, 0.826) | 0.941 (0.930, 0.951) | 0.359 | 0.302 |
| MR (field) | 1976 | -0.0863 | -0.33 | 0.259 | 0.230 | 0.89 | 0.848 (0.832, 0.863) | 0.956 (0.947, 0.965) | 0.447 | 0.398 |

**Complement Ĥᶜ block** (one-sided = upper bound)

| Estimator | n | bias vs β log-HR | bias SD units | SD_β | mean SE_β | SE/SD | Cov2s β (Wilson) | Cov1s upper β (Wilson) | log half-width | 1s margin |
|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1976 | -0.0196 | -0.25 | 0.079 | 0.079 | 1.01 | 0.953 (0.943, 0.962) | 0.923 (0.910, 0.934) | 0.156 | 0.131 |
| MR (IJ, two-term) | 1976 | +0.0066 | +0.08 | 0.082 | 0.151 | 1.84 | 0.999 (0.996, 1.000) | 0.999 (0.996, 1.000) | 0.296 | 0.248 |
| MR (IJ, winner) | 1976 | +0.0066 | +0.08 | 0.082 | 0.073 | 0.89 | 0.936 (0.924, 0.946) | 0.945 (0.934, 0.954) | 0.143 | 0.120 |
| MR (IJ, winner-floor) | 1976 | +0.0066 | +0.08 | 0.082 | 0.080 | 0.97 | 0.953 (0.943, 0.962) | 0.960 (0.950, 0.968) | 0.156 | 0.131 |
| MR (field) | 1976 | +0.0068 | +0.08 | 0.084 | 0.079 | 0.94 | 0.946 (0.936, 0.955) | 0.956 (0.947, 0.965) | 0.155 | 0.130 |

**Joint (Ĥ lower, Ĥᶜ upper)** — n = 1976; mean γ = 0.0251 (share at 0.025: 0.901); achieved joint prob. mean 0.9503 (Bonferroni 0.9505); mean corr(Λ*, Λ*ᶜ) = -0.049

| Pair | joint cov (Wilson) | cov Ĥ | cov Ĥᶜ | margin Ĥ (log) | margin Ĥᶜ (log) |
|---|---|---|---|---|---|
| separate 95% field bounds | 0.916 (0.903, 0.927) | 0.956 | 0.956 | 0.398 | 0.130 |
| Bonferroni (gamma = 0.025) | 0.957 (0.947, 0.965) | 0.974 | 0.982 | 0.484 | 0.155 |
| calibrated gamma | 0.957 (0.947, 0.965) | 0.974 | 0.982 | 0.483 | 0.155 |

### Cell h075 n=500 (campaign map1w): 2000 replicates, 1042 detections (52.1%)

**Harm Ĥ block** (one-sided = lower bound)

| Estimator | n | bias vs β log-HR | bias SD units | SD_β | mean SE_β | SE/SD | Cov2s β (Wilson) | Cov1s lower β (Wilson) | log half-width | 1s margin |
|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1042 | +0.9051 | +6.08 | 0.149 | 0.329 | 2.21 | 0.001 (0.000, 0.005) | 0.000 (0.000, 0.004) | 0.646 | 0.542 |
| MR (IJ, two-term) | 1042 | +0.3224 | +1.57 | 0.205 | 0.402 | 1.96 | 0.984 (0.974, 0.990) | 0.964 (0.950, 0.973) | 0.787 | 0.661 |
| MR (IJ, winner) | 1042 | +0.3224 | +1.57 | 0.205 | 0.175 | 0.86 | 0.559 (0.528, 0.588) | 0.465 (0.435, 0.496) | 0.343 | 0.288 |
| MR (IJ, winner-floor) | 1042 | +0.3224 | +1.57 | 0.205 | 0.329 | 1.61 | 0.944 (0.929, 0.957) | 0.876 (0.855, 0.895) | 0.646 | 0.542 |
| MR (field) | 1042 | +0.2123 | +0.91 | 0.234 | 0.359 | 1.54 | 0.988 (0.980, 0.993) | 0.953 (0.938, 0.964) | 0.697 | 0.747 |

**Complement Ĥᶜ block** (one-sided = upper bound)

| Estimator | n | bias vs β log-HR | bias SD units | SD_β | mean SE_β | SE/SD | Cov2s β (Wilson) | Cov1s upper β (Wilson) | log half-width | 1s margin |
|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1042 | -0.0848 | -0.65 | 0.131 | 0.140 | 1.07 | 0.923 (0.905, 0.938) | 0.872 (0.851, 0.891) | 0.275 | 0.230 |
| MR (IJ, two-term) | 1042 | -0.0170 | -0.12 | 0.138 | 0.255 | 1.84 | 1.000 (0.996, 1.000) | 0.995 (0.989, 0.998) | 0.499 | 0.419 |
| MR (IJ, winner) | 1042 | -0.0170 | -0.12 | 0.138 | 0.123 | 0.89 | 0.915 (0.896, 0.930) | 0.899 (0.879, 0.916) | 0.240 | 0.202 |
| MR (IJ, winner-floor) | 1042 | -0.0170 | -0.12 | 0.138 | 0.140 | 1.01 | 0.955 (0.941, 0.966) | 0.934 (0.917, 0.947) | 0.275 | 0.230 |
| MR (field) | 1042 | -0.0056 | -0.04 | 0.142 | 0.137 | 0.96 | 0.943 (0.928, 0.956) | 0.936 (0.919, 0.949) | 0.266 | 0.236 |

**Joint (Ĥ lower, Ĥᶜ upper)** — n = 1042; mean γ = 0.0250 (share at 0.025: 0.971); achieved joint prob. mean 0.9498 (Bonferroni 0.9498); mean corr(Λ*, Λ*ᶜ) = +0.044

| Pair | joint cov (Wilson) | cov Ĥ | cov Ĥᶜ | margin Ĥ (log) | margin Ĥᶜ (log) |
|---|---|---|---|---|---|
| separate 95% field bounds | 0.890 (0.869, 0.907) | 0.953 | 0.936 | 0.747 | 0.236 |
| Bonferroni (gamma = 0.025) | 0.950 (0.935, 0.962) | 0.988 | 0.962 | 0.896 | 0.280 |
| calibrated gamma | 0.950 (0.935, 0.962) | 0.988 | 0.962 | 0.896 | 0.280 |

### Cell h100 n=1000 (campaign map1w): 2000 replicates, 1319 detections (66.0%)

**Harm Ĥ block** (one-sided = lower bound)

| Estimator | n | bias vs β log-HR | bias SD units | SD_β | mean SE_β | SE/SD | Cov2s β (Wilson) | Cov1s lower β (Wilson) | log half-width | 1s margin |
|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1319 | +0.6043 | +4.26 | 0.142 | 0.252 | 1.78 | 0.269 (0.246, 0.294) | 0.161 (0.142, 0.182) | 0.494 | 0.414 |
| MR (IJ, two-term) | 1319 | +0.1806 | +0.94 | 0.191 | 0.298 | 1.56 | 0.980 (0.971, 0.987) | 0.936 (0.922, 0.948) | 0.585 | 0.491 |
| MR (IJ, winner) | 1319 | +0.1806 | +0.94 | 0.191 | 0.125 | 0.65 | 0.585 (0.558, 0.611) | 0.549 (0.522, 0.576) | 0.245 | 0.205 |
| MR (IJ, winner-floor) | 1319 | +0.1806 | +0.94 | 0.191 | 0.252 | 1.32 | 0.938 (0.923, 0.950) | 0.882 (0.863, 0.898) | 0.494 | 0.414 |
| MR (field) | 1319 | +0.1056 | +0.49 | 0.216 | 0.264 | 1.23 | 0.958 (0.945, 0.967) | 0.920 (0.905, 0.934) | 0.514 | 0.542 |

**Complement Ĥᶜ block** (one-sided = upper bound)

| Estimator | n | bias vs β log-HR | bias SD units | SD_β | mean SE_β | SE/SD | Cov2s β (Wilson) | Cov1s upper β (Wilson) | log half-width | 1s margin |
|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1319 | -0.0441 | -0.48 | 0.092 | 0.097 | 1.05 | 0.943 (0.929, 0.954) | 0.908 (0.891, 0.922) | 0.189 | 0.159 |
| MR (IJ, two-term) | 1319 | +0.0004 | +0.00 | 0.097 | 0.179 | 1.85 | 1.000 (0.997, 1.000) | 0.998 (0.993, 0.999) | 0.351 | 0.295 |
| MR (IJ, winner) | 1319 | +0.0004 | +0.00 | 0.097 | 0.087 | 0.90 | 0.927 (0.912, 0.940) | 0.939 (0.924, 0.950) | 0.170 | 0.143 |
| MR (IJ, winner-floor) | 1319 | +0.0004 | +0.00 | 0.097 | 0.097 | 1.00 | 0.958 (0.945, 0.967) | 0.953 (0.940, 0.963) | 0.189 | 0.159 |
| MR (field) | 1319 | +0.0076 | +0.08 | 0.099 | 0.095 | 0.96 | 0.945 (0.932, 0.956) | 0.953 (0.940, 0.963) | 0.185 | 0.164 |

**Joint (Ĥ lower, Ĥᶜ upper)** — n = 1319; mean γ = 0.0250 (share at 0.025: 0.968); achieved joint prob. mean 0.9499 (Bonferroni 0.9500); mean corr(Λ*, Λ*ᶜ) = +0.038

| Pair | joint cov (Wilson) | cov Ĥ | cov Ĥᶜ | margin Ĥ (log) | margin Ĥᶜ (log) |
|---|---|---|---|---|---|
| separate 95% field bounds | 0.876 (0.858, 0.893) | 0.920 | 0.953 | 0.542 | 0.164 |
| Bonferroni (gamma = 0.025) | 0.940 (0.926, 0.952) | 0.961 | 0.978 | 0.650 | 0.194 |
| calibrated gamma | 0.940 (0.926, 0.952) | 0.961 | 0.978 | 0.650 | 0.194 |

### Cell h175 knoise3 n=500 (campaign map1w): 2000 replicates, 1945 detections (97.2%)

**Harm Ĥ block** (one-sided = lower bound)

| Estimator | n | bias vs β log-HR | bias SD units | SD_β | mean SE_β | SE/SD | Cov2s β (Wilson) | Cov1s lower β (Wilson) | log half-width | 1s margin |
|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1311 | +0.4911 | +2.30 | 0.214 | 0.319 | 1.49 | 0.644 (0.617, 0.669) | 0.539 (0.511, 0.565) | 0.626 | 0.525 |
| MR (IJ, two-term) | 1311 | -0.0654 | -0.22 | 0.296 | 0.395 | 1.34 | 0.950 (0.937, 0.961) | 0.982 (0.973, 0.988) | 0.774 | 0.650 |
| MR (IJ, winner) | 1311 | -0.0654 | -0.22 | 0.296 | 0.176 | 0.60 | 0.603 (0.577, 0.629) | 0.798 (0.775, 0.819) | 0.345 | 0.290 |
| MR (IJ, winner-floor) | 1311 | -0.0654 | -0.22 | 0.296 | 0.319 | 1.08 | 0.880 (0.862, 0.897) | 0.950 (0.936, 0.960) | 0.626 | 0.525 |
| MR (field) | 1311 | -0.1422 | -0.42 | 0.338 | 0.362 | 1.07 | 0.866 (0.846, 0.883) | 0.981 (0.972, 0.987) | 0.705 | 0.717 |

**Complement Ĥᶜ block** (one-sided = upper bound)

| Estimator | n | bias vs β log-HR | bias SD units | SD_β | mean SE_β | SE/SD | Cov2s β (Wilson) | Cov1s upper β (Wilson) | log half-width | 1s margin |
|---|---|---|---|---|---|---|---|---|---|---|
| naive | 1311 | -0.0823 | -0.60 | 0.137 | 0.138 | 1.01 | 0.899 (0.882, 0.914) | 0.853 (0.833, 0.871) | 0.270 | 0.227 |
| MR (IJ, two-term) | 1311 | -0.0110 | -0.08 | 0.140 | 0.256 | 1.83 | 1.000 (0.997, 1.000) | 0.998 (0.994, 1.000) | 0.502 | 0.422 |
| MR (IJ, winner) | 1311 | -0.0110 | -0.08 | 0.140 | 0.124 | 0.88 | 0.912 (0.896, 0.926) | 0.908 (0.891, 0.922) | 0.242 | 0.203 |
| MR (IJ, winner-floor) | 1311 | -0.0110 | -0.08 | 0.140 | 0.138 | 0.99 | 0.947 (0.934, 0.958) | 0.938 (0.924, 0.950) | 0.270 | 0.227 |
| MR (field) | 1311 | -0.0017 | -0.01 | 0.143 | 0.135 | 0.94 | 0.937 (0.922, 0.949) | 0.937 (0.923, 0.949) | 0.263 | 0.233 |

**Joint (Ĥ lower, Ĥᶜ upper)** — n = 1311; mean γ = 0.0251 (share at 0.025: 0.899); achieved joint prob. mean 0.9503 (Bonferroni 0.9505); mean corr(Λ*, Λ*ᶜ) = +0.027

| Pair | joint cov (Wilson) | cov Ĥ | cov Ĥᶜ | margin Ĥ (log) | margin Ĥᶜ (log) |
|---|---|---|---|---|---|
| separate 95% field bounds | 0.919 (0.903, 0.933) | 0.981 | 0.937 | 0.711 | 0.232 |
| Bonferroni (gamma = 0.025) | 0.960 (0.948, 0.970) | 0.993 | 0.967 | 0.858 | 0.275 |
| calibrated gamma | 0.960 (0.947, 0.969) | 0.993 | 0.966 | 0.857 | 0.274 |

## Across cells

### Complement Ĥᶜ (one-sided **upper** coverage of β(Ĥᶜ); r = mean SE / empirical SD)

| Cell | two-term cov (r) | winner cov (r) | winner-floor cov (r) | field cov (r) | winner-floor two-sided | field two-sided | 1s margin: two-term / winner / floor / field |
|---|---|---|---|---|---|---|---|
| null 1.0, n=500 | 0.997 (1.85) | 0.896 (0.89) | 0.932 (1.01) | 0.930 (0.97) | 0.954 | 0.950 | 0.419 / 0.202 / 0.229 / 0.236 |
| harm 1.75, n=500 | 0.999 (1.80) | 0.913 (0.87) | 0.938 (0.97) | 0.936 (0.93) | 0.947 | 0.934 | 0.423 / 0.204 / 0.228 / 0.232 |
| harm 1.5, n=500 | 0.997 (1.81) | 0.905 (0.87) | 0.938 (0.98) | 0.936 (0.94) | 0.951 | 0.935 | 0.421 / 0.203 / 0.228 / 0.233 |
| harm 1.5, n=1500 | 0.999 (1.84) | 0.945 (0.89) | 0.960 (0.97) | 0.956 (0.94) | 0.953 | 0.946 | 0.248 / 0.120 / 0.131 / 0.130 |
| protective 0.75 | 0.995 (1.84) | 0.899 (0.89) | 0.934 (1.01) | 0.936 (0.96) | 0.955 | 0.943 | 0.419 / 0.202 / 0.230 / 0.236 |
| null 1.0, n=1000 | 0.998 (1.85) | 0.939 (0.90) | 0.953 (1.00) | 0.953 (0.96) | 0.958 | 0.945 | 0.295 / 0.143 / 0.159 / 0.164 |
| harm 1.75, 3 noise | 0.998 (1.83) | 0.908 (0.88) | 0.938 (0.99) | 0.937 (0.94) | 0.947 | 0.937 | 0.422 / 0.203 / 0.227 / 0.233 |

### Harm Ĥ (one-sided **lower** coverage of β(Ĥ))

| Cell | two-term cov (r) | winner cov (r) | winner-floor cov (r) | field cov (r) | two-sided: two-term / floor / field | retained bias β̃ (SD units) | 1s margin: two-term / winner / floor / field |
|---|---|---|---|---|---|---|---|
| null 1.0, n=500 | 0.966 (1.79) | 0.580 (0.79) | 0.892 (1.46) | 0.963 (1.41) | 0.986 / 0.950 / 0.988 | +1.10 | 0.656 / 0.288 / 0.533 / 0.734 |
| harm 1.75, n=500 | 0.983 (1.26) | 0.767 (0.57) | 0.929 (0.99) | 0.966 (1.01) | 0.961 / 0.889 / 0.887 | −0.03 | 0.646 / 0.295 / 0.510 / 0.699 |
| harm 1.5, n=500 | 0.978 (1.43) | 0.727 (0.64) | 0.926 (1.14) | 0.965 (1.14) | 0.977 / 0.923 / 0.919 | +0.21 | 0.649 / 0.292 / 0.518 / 0.710 |
| harm 1.5, n=1500 | 0.981 (1.08) | 0.863 (0.49) | 0.941 (0.81) | 0.956 (0.89) | 0.901 / 0.810 / 0.848 | −0.37 | 0.404 / 0.182 / 0.302 / 0.398 |
| protective 0.75 | 0.964 (1.96) | 0.465 (0.86) | 0.876 (1.61) | 0.953 (1.54) | 0.984 / 0.944 / 0.988 | +1.57 | 0.661 / 0.288 / 0.542 / 0.747 |
| null 1.0, n=1000 | 0.936 (1.56) | 0.549 (0.65) | 0.882 (1.32) | 0.920 (1.23) | 0.980 / 0.938 / 0.958 | +0.94 | 0.491 / 0.205 / 0.414 / 0.542 |
| harm 1.75, 3 noise | 0.982 (1.34) | 0.798 (0.60) | 0.950 (1.08) | 0.981 (1.07) | 0.950 / 0.880 / 0.866 | −0.22 | 0.650 / 0.290 / 0.525 / 0.717 |

### The simultaneous pair

| Cell | separate 95% bounds: joint cov | Bonferroni: joint cov | calibrated: joint cov | mean γ (share at 0.025) | mean corr | margins Ĥ / Ĥᶜ: separate → calibrated |
|---|---|---|---|---|---|---|
| null 1.0, n=500 | 0.894 (0.877, 0.909) | 0.955 | 0.955 (0.943, 0.965) | 0.0250 (0.968) | +0.047 | 0.734 / 0.236 → 0.883 / 0.279 |
| harm 1.75, n=500 | 0.905 | 0.956 | 0.956 (0.946, 0.964) | 0.0251 (0.924) | +0.014 | 0.699 / 0.232 → 0.846 / 0.275 |
| harm 1.5, n=500 | 0.905 | 0.953 | 0.953 (0.942, 0.962) | 0.0251 (0.943) | +0.030 | 0.710 / 0.233 → 0.858 / 0.276 |
| harm 1.5, n=1500 | 0.916 | 0.957 | 0.957 (0.947, 0.965) | 0.0251 (0.901) | −0.049 | 0.398 / 0.130 → 0.483 / 0.155 |
| protective 0.75 | 0.890 | 0.950 | 0.950 (0.935, 0.962) | 0.0250 (0.971) | +0.044 | 0.747 / 0.236 → 0.896 / 0.280 |
| null 1.0, n=1000 | 0.876 | 0.940 | 0.940 (0.926, 0.952) | 0.0250 (0.968) | +0.038 | 0.542 / 0.164 → 0.650 / 0.194 |
| harm 1.75, 3 noise | 0.919 | 0.960 | 0.960 (0.947, 0.969) | 0.0251 (0.899) | +0.027 | 0.711 / 0.232 → 0.857 / 0.274 |

### Reading criteria (Larry's, not gates)

| Criterion | Result |
|---|---|
| Complement `winner` SE/SD in [0.9, 1.1] with one-sided upper coverage within MC error of nominal in every cell | **Not met by the bare winner-only SE**: r = 0.87–0.90 (at or just under the band's floor) and upper coverage 0.896–0.945 — 1–5 points under nominal, outside the Wilson intervals in six of seven cells. The complement is *nearly* in the one-dominant-candidate regime, not exactly (a 10–13% SE shortfall). **Met by `winner_floor`** in the sense the task intends: r = 0.97–1.01, upper coverage 0.932–0.960 — the Wilson interval contains 0.95 at n ≥ 1000 and sits 1.2–1.8 points under at n = 500 (upper Wilson limits 0.944–0.950), the same footprint as the field's (0.930–0.956). Two-sided 0.947–0.958, nominal everywhere. |
| Harm-side `winner` under-covering at ties as expected; `winner_floor` at or above nominal | `winner` under-covers as the PoC predicted (r = 0.49–0.86; one-sided 0.47–0.86). **`winner_floor` is not at nominal on the harm side**: one-sided 0.876–0.950 (nominal only in the knoise3 cell), two-sided 0.810–0.950. Two regimes: in the null/protective cells it is the retained bias of β̃ (+0.9 to +1.6 SD) that the naive-SE floor does not absorb where the two-term's r ≈ 1.8–2.0 does; in the harm cells (r ≈ 0.8–1.1, bias ≤ 0.4 SD) it is the moderately separated regime the PoC flagged — β̃'s SD exceeds the naive SE by 10–25%, most at n = 1500 (r = 0.81, two-sided 0.810). The field's one-sided lower coverage (0.920–0.981) remains the better-calibrated harm-side bound; the floor is 2–8 points under it at 30% narrower margin. |
| Joint calibrated coverage within MC error of 0.95 with margins below Bonferroni's | Coverage **met**: 0.940–0.960, Wilson intervals containing 0.95 in six of seven (0.940 at h100 n1000, upper limit 0.952). Margins **not below Bonferroni's — they equal them**: γ lands at 0.025 on 90–97% of replicates and never above 0.026, because the harm and complement field draws are nearly independent (mean corr −0.05 to +0.05: a candidate's harm and complement influences have disjoint supports, and the winner's noise couples to the complement only through the second-order inner mean). Under independence the exact equal-tail solution is 1 − √0.95 = 0.0253, which the 0.001 grid rounds to 0.025. There is nothing for calibration to recover here; the value of `field$joint` is that it *shows* this. Meanwhile the two separate 95% bounds jointly cover only 0.876–0.919 — 3–7 points short — which is the quantity a development claim stating both bounds would actually have. |

### Bound locations (the Reading's inputs; means over detected replicates)

| Cell | mean β(Ĥ) | Ĥ 95% lower: naive / two-term / winner / floor / field | share floor-L ≥ 0.85 / ≥ 0.95 | mean β(Ĥᶜ) | Ĥᶜ 95% upper: naive / two-term / winner / floor / field | share floor-U < 0.85 / < 0.80 | share field-U < 0.85 / < 0.80 |
|---|---|---|---|---|---|---|---|
| null 1.0, n=500 | 0.719 | 0.95 / 0.48 / 0.70 / 0.55 / 0.45 | 0.03 / 0.01 | 0.661 | 0.77 / 0.99 / 0.80 / 0.82 / 0.83 | 0.61 / 0.44 | 0.59 / 0.41 |
| harm 1.75, n=500 | 1.253 | 1.18 / 0.64 / 0.92 / 0.74 / 0.62 | 0.26 / 0.16 | 0.669 | 0.79 / 1.02 / 0.82 / 0.84 / 0.84 | 0.55 / 0.38 | 0.54 / 0.37 |
| harm 1.5, n=500 | 1.041 | 1.08 / 0.57 / 0.82 / 0.66 / 0.55 | 0.14 / 0.08 | 0.671 | 0.78 / 1.02 / 0.82 / 0.84 / 0.84 | 0.55 / 0.39 | 0.54 / 0.37 |
| harm 1.5, n=1500 | 1.348 | 1.15 / 0.83 / 1.04 / 0.92 / 0.85 | 0.59 / 0.40 | 0.646 | 0.72 / 0.84 / 0.74 / 0.74 / 0.74 | 0.95 / 0.83 | 0.95 / 0.83 |
| protective 0.75 | 0.636 | 0.92 / 0.46 / 0.67 / 0.52 / 0.43 | 0.01 / 0.01 | 0.641 | 0.75 / 0.97 / 0.78 / 0.80 / 0.81 | 0.69 / 0.51 | 0.67 / 0.49 |
| null 1.0, n=1000 | 0.766 | 0.92 / 0.56 / 0.75 / 0.61 / 0.54 | 0.05 / 0.02 | 0.657 | 0.74 / 0.89 / 0.76 / 0.77 / 0.78 | 0.85 / 0.64 | 0.83 / 0.63 |
| harm 1.75, 3 noise | 1.306 | 1.21 / 0.62 / 0.88 / 0.70 / 0.58 | 0.19 / 0.12 | 0.664 | 0.78 / 1.02 / 0.82 / 0.84 / 0.85 | 0.54 / 0.37 | 0.52 / 0.36 |

## Reading

Read by bound location against clinically meaningful effect sizes (0.80 / 0.85 for a benefit in Ĥᶜ; the harm the lower bound can rule out in Ĥ), not as significance at HR = 1.0. **Complement.** The true conditional benefit is β(Ĥᶜ) ≈ 0.64–0.67 in every cell. The naive upper bound (0.72–0.79) over-claims — it covers only 85–92% of the time — and the two-term IJ bound sits at 0.97–1.02 at n = 500, ruling out no benefit threshold at all (U < 0.85 on 11–19% of replicates). The winner-floor and the field place the upper bound at **0.82–0.84 at n = 500, 0.77–0.78 at n = 1000, 0.74 at n = 1500** — within 0.01 of each other everywhere, both with U < 0.85 on 54–69% of n = 500 replicates and 83–95% at n ≥ 1000, and both at 93–96% upper coverage. So the complement's benefit claim has two constructions in agreement, one of which (`winner_floor`) costs no field simulation at all: the naive SE *is* the right SE for β̃ᶜ here (r = 0.97–1.01), and the two-term IJ's factor of 1.8 is the doubled same-draws term the task named. The bare winner-only SE (upper bound 0.80–0.82, coverage 0.90–0.94) under-covers by a few points and is not the one to report. **Harm.** With a naive estimate at 1.3–1.7 and β(Ĥ) at 0.72–1.35, the adjusted lower bounds rule out little meaningful harm anywhere: the winner-floor puts L at 0.52–0.74 at n = 500 (0.92 at n = 1500), the two-term at 0.46–0.64, the field at 0.43–0.62 — L ≥ 0.85 on 1–26% of n = 500 replicates (floor) and 59% at n = 1500. The spread between the constructions (0.09–0.19 log-HR between floor and field at n = 500) is what each pays for selection: the floor buys a 30% shorter margin at the price of 1–8 points of one-sided coverage in the moderately separated harm cells (0.876–0.950), where β̃'s dispersion exceeds the naive SE by 10–25% and only the field's λ-SD tracks it. **The pair.** A claim stating both bounds from two separate 95% constructions holds jointly only 88–92% of the time; the Bonferroni pair holds 94–96% at a cost of 0.10–0.15 (Ĥ) and 0.025–0.045 (Ĥᶜ) log-HR of margin, and the calibrated pair *is* the Bonferroni pair in these cells because the two draws are nearly independent (corr −0.05 to +0.05) — the joint construction's finding is that there is no dependence to exploit, which is itself worth knowing before promising a tighter pair.

No task proposed; nothing blocked. Findings in the record. K-5 (whether `winner_floor` becomes the complement's reported SE and whether the pair enters the standard output) is Larry's call on this record: on the evidence, the complement's `winner_floor` matches the field at zero cost; the harm side should stay on the two-term or the field.
