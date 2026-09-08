# REPORT — The identifier band dial on the J = 10 harm cells: minSG, effMaxSG at ε = 0.10 / 0.20 / 0.30, maxSG (final)

**Task:** `dev/tasks/TASK_complement_variance_banddial_2026-09-07.md` (4b245516), Part B as instructed for this run (`maxSG` and `minSG` at their own definitions, `effMaxSG` at ε = 0.30; the two J = 10 harm cells). Records: Stage 0 3dd238fc (Gate 0 PASS, all three arms) · Stage 1 / Gate 1 PASS c133902e (knob-inert identity exact; smoke and nesting; alignment 4/4; projection 2.7 h) · Gate 2 record beside this file (`REPORT_banddial_gate2_2026-09-07.md`: six cells PASS, 2 h 36 m, none deferred). Parts A (fe5eeaf5) and C (b49be1a6) precede. Decisions T-1–T-5 at their defaults or as instructed (no ε = 0.40 setting: `maxSG` and ε = 0.30 differ materially, but a fourth band was not the instruction). **Winner-only and winner-floor excluded** from every table, figure and line below.
**Date:** 2026-09-07. Data: five settings on identical replicates (`n_true` and `truth` identical, the same family K on every replicate; 1,999 of 2,000 detected under every setting — the same undetected replicate): `minSG` (`banddial`), `effMaxSG` ε 0.10 (`p30sg`), ε 0.20 (nb20 arm A `p30sgnb20`), ε 0.30 (`banddial`), `maxSG` (`banddial`); `ci_method = "field"`, `include_complement = TRUE`, `field_complement = TRUE`, `ij_residual = "two_term"`, FB none. Rendered: the six `…_combine_1_2000.html` and the cross-setting `summary_banddial.html` (+ `.qmd`).

**Conventions.** As the nb20 report: log-HR scale against β(Ĥ), β(Ĥᶜ); one-sided coverage on the exposed side (Ĥ lower, Ĥᶜ upper); Wilson 95%; bound locations against 0.85 / 0.95 (Ĥ lower) and 0.85 / 0.80 (Ĥᶜ upper). "SD units" and "SD(β̃ᶜ)/naive SE" are the reports' marginal-SD convention; the error-scale ratio of Part A (`REPORT_complement_variance_2026-09-07.md`, A0) is given beside it wherever the two differ. The band is chosen by Larry on the capture / specificity trade-off; this record does not recommend one.

## 1. The dial table (paired by replicate; n = 1,999 per cell)

### HR 1.50 n500 (true harm 153 of 500; θ†(H) = 1.499, θ†(Hᶜ) = 0.721)

| Quantity | minSG | ε 0.10 | ε 0.20 | ε 0.30 | maxSG |
|---|---|---|---|---|---|
| detection | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| mean \|Ĥ\| | 62.5 | 97.3 | 127.1 | 168.9 | 367.3 |
| \|Ĥ\|/\|H\| median (q10, q90) | 0.41 (0.37, 0.45) | 0.60 (0.43, 0.93) | 0.79 (0.53, 1.18) | 1.09 (0.70, 1.48) | 2.52 (1.54, 3.07) |
| share \|Ĥ\| ≥ \|H\| (at or above the planted size) | 0.001 | 0.064 | 0.264 | 0.634 | 0.989 |
| share \|Ĥ\| ≥ 0.9\|H\| | 0.001 | 0.113 | 0.391 | 0.743 | 0.993 |
| sensitivity | 0.248 | 0.450 | 0.590 | 0.737 | 0.915 |
| specificity | 0.929 | 0.917 | 0.893 | 0.838 | 0.344 |
| PPV | 0.603 | 0.697 | 0.706 | 0.682 | 0.402 |
| NPV | 0.737 | 0.794 | 0.837 | 0.883 | 0.915 |
| naive optimism Ĥ: log-HR (marginal SD units; error-SD units) | +0.384 (+2.03; +1.53) | +0.548 (+2.48; +1.98) | +0.430 (+1.81; +1.53) | +0.299 (+1.33; +1.23) | +0.082 (+1.47; +0.90) |
| naive optimism Ĥᶜ: log-HR (marginal SD units) | −0.056 (−0.40) | −0.112 (−0.76) | −0.123 (−0.79) | −0.137 (−0.84) | −0.225 (−0.67) |
| mean β(Ĥ) vs θ†(H) = 1.499 | 1.163 | 1.252 | 1.256 | 1.228 | 0.977 |
| mean β(Ĥᶜ) vs θ†(Hᶜ) = 0.721 | 0.874 | 0.835 | 0.808 | 0.779 | 0.762 |
| complement SD(β̃ᶜ)/naive SE (marginal; error ratio SD(eᶜ)/naive SE) | 1.06 (1.04) | 1.08 (1.01) | 1.10 (1.01) | 1.07 (1.00) | 1.00 (1.02) |
| complement λ-SDᶜ/naive SE | 0.99 | 0.97 | 0.95 | 0.92 | 1.02 |
| harm SD(β̃)/naive SE (marginal; error) | 0.73 (0.90) | 0.94 (1.06) | 1.10 (1.19) | 1.18 (1.19) | 0.46 (0.87) |
| p̂(Ĥ): mean (median); share = argmax; top-3 mass | 0.092 (0.043); 0.10; 0.74 | not recorded | 0.103 (0.069); 0.25; 0.36 | 0.064 (0.042); 0.15; 0.30 | 0.291 (0.127); 0.39; 0.56 |
| band on the observed effects (mean) | 3.4 (ε 0.10, informational) | not recorded | 10.3 | 29.1 | 3.4 (informational) |
| consistency-qualifying; family K | 202; 1228 | — | 202; 1228 | 202; 1228 | 202; 1228 |
| complement fits/rep | 133 | 484 | 615 | 749 | 586 |

### HR 1.75 n500 (true harm 153; θ†(H) = 1.746, θ†(Hᶜ) = 0.721)

| Quantity | minSG | ε 0.10 | ε 0.20 | ε 0.30 | maxSG |
|---|---|---|---|---|---|
| detection | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| mean \|Ĥ\| | 61.8 | 99.3 | 128.9 | 166.8 | 401.6 |
| \|Ĥ\|/\|H\| median (q10, q90) | 0.40 (0.37, 0.45) | 0.62 (0.43, 0.95) | 0.82 (0.54, 1.16) | 1.09 (0.70, 1.43) | 2.75 (1.93, 3.13) |
| share \|Ĥ\| ≥ \|H\| | 0.001 | 0.069 | 0.302 | 0.647 | 0.999 |
| share \|Ĥ\| ≥ 0.9\|H\| | 0.001 | 0.135 | 0.429 | 0.768 | 0.999 |
| sensitivity | 0.249 | 0.506 | 0.660 | 0.792 | 0.948 |
| specificity | 0.931 | 0.936 | 0.919 | 0.868 | 0.260 |
| PPV | 0.613 | 0.767 | 0.779 | 0.740 | 0.372 |
| NPV | 0.738 | 0.815 | 0.866 | 0.909 | 0.936 |
| naive optimism Ĥ: log-HR (marginal SD units; error-SD units) | +0.326 (+1.50; +1.20) | +0.483 (+2.12; +1.67) | +0.365 (+1.48; +1.24) | +0.254 (+1.07; +1.00) | +0.046 (+0.74; +0.49) |
| naive optimism Ĥᶜ: log-HR (marginal SD units) | −0.049 (−0.35) | −0.098 (−0.65) | −0.102 (−0.63) | −0.113 (−0.69) | −0.191 (−0.49) |
| mean β(Ĥ) vs θ†(H) = 1.746 | 1.304 | 1.512 | 1.515 | 1.454 | 1.014 |
| mean β(Ĥᶜ) vs θ†(Hᶜ) = 0.721 | 0.915 | 0.848 | 0.808 | 0.774 | 0.758 |
| complement SD(β̃ᶜ)/naive SE (marginal; error ratio) | 1.06 (1.03) | 1.12 (1.01) | 1.14 (1.02) | 1.09 (1.02) | 1.04 (1.05) |
| complement λ-SDᶜ/naive SE | 0.99 | 0.97 | 0.95 | 0.93 | 0.99 |
| harm SD(β̃)/naive SE (marginal; error) | 0.80 (0.95) | 1.00 (1.13) | 1.18 (1.28) | 1.27 (1.26) | 0.59 (0.92) |
| p̂(Ĥ): mean (median); share = argmax; top-3 mass | 0.086 (0.033); 0.09; 0.78 | not recorded | 0.119 (0.083); 0.28; 0.40 | 0.080 (0.055); 0.18; 0.33 | 0.415 (0.299); 0.54; 0.66 |
| band on the observed effects (mean) | 3.3 (informational) | not recorded | 9.9 | 27.8 | 3.3 (informational) |
| consistency-qualifying; family K | 277; 1228 | — | 277; 1228 | 277; 1228 | 277; 1228 |
| complement fits/rep | 98 | 457 | 580 | 708 | 465 |

**Reading of the dial.** The five settings are nested on every replicate (minSG ≤ ε 0.10 ≤ ε 0.20 ≤ ε 0.30 ≤ maxSG; the J = 10 family and the qualifying set are identical, only the pick moves), so the dial is one monotone path through the same candidates. **Capture:** the median |Ĥ|/|H| runs 0.40 → 0.60 → 0.79 → 1.09 → 2.5–2.8; the share at or above the planted size 0.00 → 0.06 → 0.26 → 0.63 → 0.99; sensitivity 0.25 → 0.45 / 0.51 → 0.59 / 0.66 → 0.74 / 0.79 → 0.92 / 0.95. **Its price:** specificity 0.93 → 0.92–0.94 → 0.89–0.92 → 0.84–0.87 → 0.26–0.34; PPV peaks at ε 0.20 (0.71 / 0.78) and falls to 0.40 / 0.37 under `maxSG`. ε = 0.30 is the first setting whose median region reaches the planted size — and it overshoots on the top third (q90 1.43–1.48, |Ĥ| up to 368). `maxSG` is not a subgroup rule at this design: it returns 73–80% of the sample (the broadest candidate clearing the harm screen and the consistency floor, with a de-biased β(Ĥ) of 0.98 / 1.01 — the overall population's effect), and `minSG` sits on the n.min floor (61–64 on 99% of replicates), a two-factor rule of 60 patients with sensitivity 0.25. **β(Ĥ)** is flat at 1.25–1.26 / 1.51–1.52 from ε 0.10 to 0.20, dips to 1.23 / 1.45 at ε 0.30 (the overshoot dilutes the region) and collapses to ≈ 1.0 under `maxSG`; **β(Ĥᶜ)** moves toward θ†(Hᶜ) = 0.721 monotonically (0.87 / 0.92 → 0.84 / 0.85 → 0.81 → 0.78 / 0.77 → 0.76). **Naive optimism on Ĥ** falls along the dial in log-HR (+0.55 → +0.43 → +0.30 → +0.08 at HR 1.50); in marginal-SD units it is +2.5 → +1.8 → +1.3, then +1.5 under `maxSG` because the marginal SD collapses (0.055–0.063: the pick is nearly the same broad rule every time); in error-SD units +2.0 → +1.5 → +1.2 → +0.9. `minSG`'s optimism is +2.0 / +1.5 marginal SD (a 60-patient region carries a large winner's curse in absolute terms: +0.38 / +0.33 log-HR). **Regime:** p̂(Ĥ) is lowest at ε 0.30 (mean 0.06–0.08, argmax on 15–18% of replicates) and highest under `maxSG` (0.29–0.42; argmax 39–54%; top-3 mass 0.56–0.66) — the size rule re-selects the same broad candidate; `minSG` is competitive (p̂ 0.09, top-3 mass 0.74–0.78 concentrated on a few smallest rules). **Complement regime:** the marginal diagnostic SD(β̃ᶜ)/naive SE is 1.06–1.14 at ε 0.10–0.30 and 1.00–1.04 at the two ends; the error ratio is 1.00–1.05 everywhere (Part A); λ-SDᶜ/naive SE falls 0.99 → 0.97 → 0.95 → 0.92–0.93 along the band and returns to 0.99–1.02 under `maxSG`.

## 2. Constructions: one-sided coverage on both blocks at each setting (Wilson 95%; b = bias / marginal SD, r = mean SE / marginal SD)

### HR 1.50 n500

| Setting | Ĥ lower: naive | Ĥ lower: IJ two-term (r) | Ĥ lower: field (r) | Ĥᶜ upper: naive | Ĥᶜ upper: IJ two-term (r) | Ĥᶜ upper: field (r) |
|---|---|---|---|---|---|---|
| minSG | 0.752 (0.733, 0.771) | 0.974 (0.967, 0.981) (1.62) | 0.942 (0.931, 0.951) (1.18) | 0.867 (0.852, 0.882) | 0.997 (0.993, 0.999) (1.82) | 0.923 (0.911, 0.934) (0.94) |
| ε 0.10 | 0.381 (0.360, 0.402) | 0.978 (0.971, 0.984) (1.41) | 0.970 (0.962, 0.977) (1.12) | 0.793 (0.775, 0.810) | 0.996 (0.993, 0.998) (1.70) | 0.913 (0.900, 0.925) (0.89) |
| ε 0.20 | 0.484 (0.462, 0.506) | 0.985 (0.979, 0.990) (1.29) | 0.974 (0.967, 0.981) (1.08) | 0.781 (0.762, 0.798) | 0.995 (0.991, 0.997) (1.64) | 0.897 (0.883, 0.910) (0.85) |
| ε 0.30 | 0.633 (0.612, 0.654) | 0.995 (0.991, 0.997) (1.31) | 0.980 (0.973, 0.986) (1.10) | 0.774 (0.755, 0.792) | 0.996 (0.992, 0.998) (1.63) | 0.904 (0.890, 0.916) (0.85) |
| maxSG | 0.943 (0.932, 0.952) | 0.999 (0.997, 1.000) (3.33) | 0.949 (0.938, 0.958) (1.69) | 0.726 (0.706, 0.745) | 0.995 (0.991, 0.997) (1.59) | 0.928 (0.916, 0.939) (0.97) |

### HR 1.75 n500

| Setting | Ĥ lower: naive | Ĥ lower: IJ two-term (r) | Ĥ lower: field (r) | Ĥᶜ upper: naive | Ĥᶜ upper: IJ two-term (r) | Ĥᶜ upper: field (r) |
|---|---|---|---|---|---|---|
| minSG | 0.804 (0.786, 0.821) | 0.971 (0.963, 0.977) (1.51) | 0.936 (0.924, 0.946) (1.11) | 0.879 (0.864, 0.893) | 0.997 (0.994, 0.999) (1.83) | 0.929 (0.917, 0.940) (0.94) |
| ε 0.10 | 0.474 (0.452, 0.496) | 0.982 (0.975, 0.987) (1.34) | 0.964 (0.955, 0.972) (1.08) | 0.819 (0.801, 0.835) | 0.996 (0.993, 0.998) (1.65) | 0.919 (0.906, 0.930) (0.86) |
| ε 0.20 | 0.573 (0.551, 0.595) | 0.985 (0.979, 0.989) (1.23) | 0.970 (0.962, 0.977) (1.04) | 0.826 (0.809, 0.842) | 0.995 (0.991, 0.998) (1.59) | 0.912 (0.899, 0.924) (0.82) |
| ε 0.30 | 0.680 (0.660, 0.700) | 0.994 (0.990, 0.997) (1.23) | 0.976 (0.968, 0.982) (1.04) | 0.807 (0.790, 0.824) | 0.996 (0.993, 0.998) (1.60) | 0.903 (0.890, 0.916) (0.83) |
| maxSG | 0.947 (0.936, 0.956) | 0.999 (0.997, 1.000) (2.77) | 0.947 (0.936, 0.956) (1.46) | 0.776 (0.757, 0.794) | 0.995 (0.991, 0.997) (1.60) | 0.925 (0.913, 0.936) (0.92) |

Field two-sided coverage: Ĥ 0.94 / 0.92 (minSG), 0.94 / 0.93 (ε 0.10), 0.93 / 0.91 (ε 0.20), 0.94 / 0.93 (ε 0.30), 0.95 / 0.96 (maxSG); Ĥᶜ 0.94 / 0.94, 0.94 / 0.94, 0.92 / 0.92, 0.92 / 0.91, 0.94 / 0.94. Field bias on Ĥ (log-HR): +0.03 / +0.02 (minSG), −0.00 / −0.02 (ε 0.10), −0.06 / −0.08 (ε 0.20), −0.09 / −0.10 (ε 0.30), +0.00 / −0.01 (maxSG); on Ĥᶜ −0.01 to −0.04 everywhere.

**Bound locations** (mean; shares against the reading thresholds)

| Cell | Setting | Ĥ field lower: mean (share ≥ 0.85 / ≥ 0.95) | Ĥ IJ two-term lower | Ĥᶜ field upper: mean (share < 0.85 / < 0.80) | Ĥᶜ IJ two-term upper | mean β(Ĥ) / β(Ĥᶜ) |
|---|---|---|---|---|---|---|
| HR 1.50 | minSG | 0.71 (0.21 / 0.12) | 0.68 (0.13 / 0.05) | 1.08 (0.05 / 0.02) | 1.31 (0.00 / 0.00) | 1.16 / 0.87 |
| HR 1.50 | ε 0.10 | 0.71 (0.22 / 0.13) | 0.73 (0.22 / 0.12) | 1.02 (0.13 / 0.06) | 1.22 (0.01 / 0.00) | 1.25 / 0.84 |
| HR 1.50 | ε 0.20 | 0.69 (0.20 / 0.12) | 0.71 (0.18 / 0.10) | 0.99 (0.19 / 0.11) | 1.19 (0.02 / 0.01) | 1.26 / 0.81 |
| HR 1.50 | ε 0.30 | 0.68 (0.17 / 0.09) | 0.67 (0.13 / 0.06) | 0.96 (0.25 / 0.16) | 1.16 (0.03 / 0.01) | 1.23 / 0.78 |
| HR 1.50 | maxSG | 0.80 (0.16 / 0.03) | 0.69 (0.01 / 0.00) | 1.30 (0.14 / 0.10) | 1.97 (0.02 / 0.01) | 0.98 / 0.76 |
| HR 1.75 | minSG | 0.79 (0.32 / 0.21) | 0.74 (0.21 / 0.12) | 1.13 (0.03 / 0.01) | 1.37 (0.00 / 0.00) | 1.30 / 0.92 |
| HR 1.75 | ε 0.10 | 0.84 (0.40 / 0.28) | 0.85 (0.43 / 0.29) | 1.04 (0.11 / 0.06) | 1.25 (0.01 / 0.00) | 1.51 / 0.85 |
| HR 1.75 | ε 0.20 | 0.82 (0.37 / 0.26) | 0.83 (0.38 / 0.24) | 1.00 (0.18 / 0.10) | 1.21 (0.02 / 0.01) | 1.52 / 0.81 |
| HR 1.75 | ε 0.30 | 0.81 (0.34 / 0.22) | 0.79 (0.31 / 0.18) | 0.97 (0.25 / 0.15) | 1.17 (0.03 / 0.01) | 1.45 / 0.77 |
| HR 1.75 | maxSG | 0.83 (0.27 / 0.08) | 0.71 (0.02 / 0.00) | 1.39 (0.14 / 0.10) | 2.32 (0.02 / 0.01) | 1.01 / 0.76 |

**Joint (Ĥ lower, Ĥᶜ upper), field:** Bonferroni / calibrated / separate — HR 1.50: minSG 0.939 / 0.938 / 0.869; ε 0.10 0.940 / 0.939 / 0.884; ε 0.20 0.932 / 0.932 / 0.872; ε 0.30 0.935 / 0.935 / 0.886; maxSG 0.937 / 0.936 / 0.881. HR 1.75: 0.930 / 0.929 / 0.869; 0.940 / 0.940 / 0.886; 0.933 / 0.933 / 0.885; 0.935 / 0.935 / 0.881; 0.942 / 0.941 / 0.876. Mean γ 0.025 (0.026 under maxSG); corr(Λ*, Λ*ᶜ) +0.01 to +0.05. Margins (log): Ĥ 0.60 / 0.66 / 0.62 / 0.56 / 0.21 (Bonferroni 0.72 / 0.80 / 0.76 / 0.69 / 0.25); Ĥᶜ 0.22 / 0.23 / 0.24 / 0.25 / 0.54 at HR 1.50.

**Reading of the constructions.** **Harm block:** the field's one-sided lower bound holds at 0.94–0.98 across the whole dial (0.94 under minSG and maxSG, 0.96–0.98 in the band settings; two-sided 0.91–0.96); the IJ two-term at 0.97–0.999 with r 1.2–1.6 (3.3 / 2.8 under maxSG, whose marginal SD is tiny). The field's retained bias on Ĥ grows with the band (−0.00 → −0.06 → −0.09 log-HR; −0.3 marginal SD at ε 0.30) — the correction over-shoots as the region overshoots. The lower bound's location is what changes along the dial: under `maxSG` the field lower bound is 0.80 / 0.83 against β(Ĥ) ≈ 1.0 (a bound on the whole population), under `minSG` 0.71 / 0.79 against 1.16 / 1.30; the band settings sit at 0.68–0.71 / 0.81–0.84 against 1.23–1.26 / 1.45–1.52, with the share ≥ 0.85 highest at ε 0.10 (0.22 / 0.40) and lowest at ε 0.30 (0.17 / 0.34). **Complement block:** the field's upper coverage is 0.90–0.92 in the three band settings (r 0.82–0.89) and 0.92–0.93 at the two ends (minSG r 0.94, maxSG r 0.92–0.97), never reaching 0.95; the IJ two-term 0.995–0.997 with an upper bound of 1.2–1.4 (2.0–2.3 under maxSG) that rules out nothing. The field's upper bound moves toward the design location as the band widens: 1.02 / 1.04 (ε 0.10) → 0.99 / 1.00 → 0.96 / 0.97 (ε 0.30; below 0.85 on 25%), against β(Ĥᶜ) 0.84 → 0.81 → 0.78; under `minSG` it is 1.08 / 1.13 (β(Ĥᶜ) 0.87 / 0.92: the complement is 88% of the sample and carries most of the planted harm) and under `maxSG` 1.30 / 1.39 (a 100–130-patient complement, naive SE 0.31–0.36). The joint Bonferroni pair holds at 0.93–0.94 at every setting (Wilson upper limits 0.94–0.95), calibrated = Bonferroni. **Against Part A:** λ-SDᶜ/naive SE falls 0.99 → 0.92 along the band as the spread of |Ĥ|/|H| widens (q10–q90 0.37–0.45 → 0.70–1.48), and the field's complement under-coverage (0.90–0.92) tracks it; under `maxSG`, where the pick is nearly the same broad rule every time (p̂ 0.29–0.42) and the complement's identity barely varies, λ-SDᶜ is back at the naive SE (1.02 / 0.99) and the field's upper coverage at 0.93.

## 3. Across settings (one row per cell and setting)

| Cell | Setting | det | mean \|Ĥ\| (true) | \|Ĥ\|/\|H\| median | share ≥ \|H\| | sens | spec | PPV | β(Ĥ) / β(Ĥᶜ) | Ĥ field 1s-lower cov (r) | Ĥ IJ 1s-lower cov (r) | Ĥᶜ field 1s-upper cov (r) | Ĥᶜ IJ 1s-upper cov (r) | Ĥᶜ SD ratio marginal (error) | λ-SDᶜ/naive SE | joint Bonferroni (Wilson) | p̂(Ĥ) mean |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HR 1.50 | minSG | 1.000 | 63 (153) | 0.41 | 0.001 | 0.248 | 0.929 | 0.603 | 1.16 / 0.87 | 0.942 (1.18) | 0.974 (1.62) | 0.923 (0.94) | 0.997 (1.82) | 1.06 (1.04) | 0.99 | 0.939 (0.928, 0.949) | 0.092 |
| HR 1.50 | ε 0.10 | 1.000 | 97 (153) | 0.60 | 0.064 | 0.450 | 0.917 | 0.697 | 1.25 / 0.84 | 0.970 (1.12) | 0.978 (1.41) | 0.913 (0.89) | 0.996 (1.70) | 1.08 (1.01) | 0.97 | 0.940 (0.929, 0.950) | — |
| HR 1.50 | ε 0.20 | 1.000 | 127 (153) | 0.79 | 0.264 | 0.590 | 0.893 | 0.706 | 1.26 / 0.81 | 0.974 (1.08) | 0.985 (1.29) | 0.897 (0.85) | 0.995 (1.64) | 1.10 (1.01) | 0.95 | 0.932 (0.920, 0.942) | 0.103 |
| HR 1.50 | ε 0.30 | 1.000 | 169 (153) | 1.09 | 0.634 | 0.737 | 0.838 | 0.682 | 1.23 / 0.78 | 0.980 (1.10) | 0.995 (1.31) | 0.904 (0.85) | 0.996 (1.63) | 1.07 (1.00) | 0.92 | 0.935 (0.923, 0.945) | 0.064 |
| HR 1.50 | maxSG | 1.000 | 367 (153) | 2.52 | 0.989 | 0.915 | 0.344 | 0.402 | 0.98 / 0.76 | 0.949 (1.69) | 0.999 (3.33) | 0.928 (0.97) | 0.995 (1.59) | 1.00 (1.02) | 1.02 | 0.937 (0.925, 0.947) | 0.291 |
| HR 1.75 | minSG | 1.000 | 62 (153) | 0.40 | 0.001 | 0.249 | 0.931 | 0.613 | 1.30 / 0.92 | 0.936 (1.11) | 0.971 (1.51) | 0.929 (0.94) | 0.997 (1.83) | 1.06 (1.03) | 0.99 | 0.930 (0.918, 0.940) | 0.086 |
| HR 1.75 | ε 0.10 | 1.000 | 99 (153) | 0.62 | 0.069 | 0.506 | 0.936 | 0.767 | 1.51 / 0.85 | 0.964 (1.08) | 0.982 (1.34) | 0.919 (0.86) | 0.996 (1.65) | 1.12 (1.01) | 0.97 | 0.940 (0.929, 0.950) | — |
| HR 1.75 | ε 0.20 | 1.000 | 129 (153) | 0.82 | 0.302 | 0.660 | 0.919 | 0.779 | 1.52 / 0.81 | 0.970 (1.04) | 0.985 (1.23) | 0.912 (0.82) | 0.995 (1.59) | 1.14 (1.02) | 0.95 | 0.933 (0.922, 0.944) | 0.119 |
| HR 1.75 | ε 0.30 | 1.000 | 167 (153) | 1.09 | 0.647 | 0.792 | 0.868 | 0.740 | 1.45 / 0.77 | 0.976 (1.04) | 0.994 (1.23) | 0.903 (0.83) | 0.996 (1.60) | 1.09 (1.02) | 0.93 | 0.935 (0.923, 0.945) | 0.080 |
| HR 1.75 | maxSG | 1.000 | 402 (153) | 2.75 | 0.999 | 0.948 | 0.260 | 0.372 | 1.01 / 0.76 | 0.947 (1.46) | 0.999 (2.77) | 0.925 (0.92) | 0.995 (1.60) | 1.04 (1.05) | 0.99 | 0.942 (0.931, 0.952) | 0.415 |

## Reading

By location and coverage, never as significance at 1.0; the band choice is Larry's. **The dial.** The three bands and the two size rules are one nested path through the same J = 10 candidates: each step outward adds members and sensitivity (0.25 → 0.45 / 0.51 → 0.59 / 0.66 → 0.74 / 0.79 → 0.92 / 0.95) and takes specificity (0.93 → 0.92–0.94 → 0.89–0.92 → 0.84–0.87 → 0.26–0.34). ε = 0.30 is where the median region first reaches the planted size (|Ĥ|/|H| 1.09; at or above |H| on 63–65% of replicates) and it does so by overshooting on a third of them (q90 1.4–1.5), which costs 5 points of specificity against ε = 0.20 and 2–4 points of PPV, and moves β(Ĥ) down (1.23 / 1.45 against 1.26 / 1.52) for the first time on the dial; β(Ĥᶜ) reaches 0.78 / 0.77 (θ†(Hᶜ) 0.721). The two ends are not subgroup rules at this design: `maxSG` returns three-quarters of the sample with β(Ĥ) ≈ 1.0 and a field lower bound of 0.80–0.83 that says nothing about a harmed region, and `minSG` returns the 60-patient floor with a winner's curse of +0.33–0.38 log-HR. **Constructions along the dial.** The harm-side field lower bound is inside 0.94–0.98 at every setting — robust to the pick as it was to the band and the grid at nb20 — with its location moving from 0.68–0.71 / 0.81–0.84 (the bands) to 0.80 / 0.83 under `maxSG`; its retained bias grows with the band (to −0.09 / −0.10 at ε 0.30). The complement-side field upper bound covers at 0.90–0.92 in every band setting and 0.92–0.93 at the ends, tracking λ-SDᶜ/naive SE (0.92–0.97 in the bands, 0.99–1.02 at the ends) as Part A predicted: the shortfall follows the spread of |Ĥ|/|H| across replicates, not its level, and vanishes where the pick stops varying. The joint Bonferroni pair holds at 0.93–0.94 everywhere. **What this says for the choice.** Between ε 0.20 and 0.30 the trade is 14 points of sensitivity and a doubling of the share reaching the planted size against 5 points of specificity, a slightly diluted β(Ĥ), and 1–2 more points of complement under-coverage; nothing on the constructions side breaks at either. Larry's call.

## Side issues (flagged, not fixed)

1. **γ off the Bonferroni floor under `maxSG`** on 44–50% of replicates (max 0.030) where the complement is 100–130 patients; the calibrated pair still equals Bonferroni to three decimals. Informational.
2. **A leftover process from a session of 2026-09-05** (`free -m` memory monitor loop in a bash from the `uburst` renders, pid 2086198) is still alive on the machine; it takes no CPU (load 4 on 128 cores at the time of Stage 1) and was not touched.
3. Part A's finding stands over this report too: every "SD units" and "SD(β̃ᶜ)/naive SE" in the tables above is on the marginal SD; the error-scale values are printed beside them where they differ. Whether the tables should switch conventions is Larry's call.

No task proposed; nothing blocked. Findings in the record.
