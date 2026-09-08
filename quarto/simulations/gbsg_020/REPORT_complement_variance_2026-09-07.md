# REPORT — Decomposition of the complement's error variance across the regime sequence (Part A, no compute)

**Task:** `dev/tasks/TASK_complement_variance_banddial_2026-09-07.md` (4b245516), Part A. Decision T-1 at default (cells as listed). Standing conventions: winner-only and winner-floor excluded; bounds by location; Wilson intervals; bias in SD units; verify from source.
**Date:** 2026-09-07. Executor: Claude Code, unattended. No compute, no `R/` change. Data: the committed pooled bundles of `s7c` (HR 1.00 / 1.75 n500; M1 prevalence 12.4%, `maxeffCons`), `p30` (HR 1.00 / 1.50 / 1.75 n500, HR 1.00 n1000; 31%, `maxeffCons`), `p30sg` (the same four cells; `effMaxSG` ε 0.10, J 10), nb20 arm A `p30sgnb20` (HR 1.50 / 1.75 n500; ε 0.20, J 10) and arm B `p30sgnb20j20` (HR 1.00 / 1.50 / 1.75 n500, HR 1.50 / 1.00 n1000; ε 0.20, J 20) — 17 cells. Document: `summary_complement_variance.qmd` (rendered `.html` beside it; reads the bundles by campaign glob). Detected replicates with every input finite (n per cell = the detections; none dropped).

**Objects, per detected replicate, log scale:** a = log(`nv_Hc_est`) − log(`betaHhat_Hc`); c = log(`nv_Hc_est`) − log(`mr_Hc_est`); e = a − c = log(`mr_Hc_est`) − log(`betaHhat_Hc`). Reported variances: mean `nv_Hc_se`² (naive SE²), mean `fld_Hc_se`² (λ-SDᶜ², the SD of Λ*ᶜ over the field's outer draws), mean `mr_Hc_se_ij`² (IJ two-term). **The complement's `selection_bias` / `fixed_bias` split of c is not recorded in any bundle** (the gate returns the two components for the harm block only; the template records neither for the complement), so c is read whole.

---

## A0 — A definitional finding that changes the reading of every "regime diagnostic" so far

The constructions tables of every report (`fs_sim_bias_coverage()`: `sd_emp = sd(log est)`; the template's `.sd_log()`) take the **marginal SD of the log estimate across replicates**, not the SD of its error against the per-replicate target. So "SD_β", "SE/SD", "bias in SD units" and the regime diagnostic **"SD(β̃ᶜ)/naive SE = 1.08–1.20"** all carry the spread of the target β(Ĥᶜ) across replicates: Var(log β̃ᶜ) = Var(e) + Var(log β(Ĥᶜ)) + 2Cov(e, log β(Ĥᶜ)).

| Cell (HR 1.75 n500) | SD(log β̃ᶜ) (reported) | SD(log β(Ĥᶜ)) | corr(e, target) | SD(e) | naive SE | reported ratio SD(β̃ᶜ)/naive SE | error ratio SD(e)/naive SE | λ-SDᶜ/SD(e) | λ-SDᶜ/SD(β̃ᶜ) (reported) |
|---|---|---|---|---|---|---|---|---|---|
| s7c (12%, maxeffCons) | 0.143 | 0.058 | −0.17 | 0.141 | 0.139 | 1.03 | 1.02 | 0.97 | 0.95 |
| p30 (31%, maxeffCons) | 0.140 | 0.056 | −0.10 | 0.134 | 0.133 | 1.06 | 1.01 | 0.97 | 0.93 |
| p30sg (ε 0.10) | 0.152 | 0.081 | −0.11 | 0.138 | 0.137 | 1.12 | 1.01 | 0.95 | 0.87 |
| nb20 A (ε 0.20, J 10) | 0.163 | 0.099 | −0.16 | 0.147 | 0.143 | 1.14 | 1.02 | 0.92 | 0.83 |
| nb20 B (ε 0.20, J 20) | 0.164 | 0.101 | −0.16 | 0.146 | 0.143 | 1.15 | 1.02 | 0.92 | 0.82 |

HR 1.50 n500 reads the same (reported ratio 1.03 → 1.08 → 1.10 → 1.10; error ratio 1.00 → 1.01 → 1.01 → 1.01; SD(target) 0.046 → 0.066 → 0.082 → 0.084); harm 1.5 n1000: reported 1.20, error ratio 1.00, SD(target) 0.084. The null cells: reported 0.98–1.01, error ratio 0.95–0.99, SD(target) 0.02–0.04.

**Reading.** The rise of the regime diagnostic across the sequence (1.03 → 1.15 at n = 500, 1.20 at n = 1000) is the growth of the **target's** spread as the region grows (SD(log β(Ĥᶜ)) 0.056 → 0.101; the complement's identity varies with the winner and so does the complement's true effect), partly offset by a negative correlation between the error and the target (−0.10 → −0.16). The **error** variance of the de-biased complement estimate is 1.00–1.02 × the naive SE² in every harm cell of every campaign. What has been called "the complement's variability growing with |Ĥ|" is the target moving, not the estimator's error. This does not change any coverage number (coverage is of β(Ĥᶜ) by the bound, an error-scale event), but it changes what the diagnostic measures and what a repair would have to match: **the field's λ-SDᶜ must be read against SD(e), not against SD(log β̃ᶜ)** — 0.92–0.97 of it, not 0.82–0.95.

## A1 — Per cell: the three components, the reported variances, the ratios

Variances on the (log-HR)² scale; the identity Var(e) = Var(a) + Var(c) − 2Cov(a, c) holds to ≤ 4e-18 in every cell (column `check` of the document).

| Cell | regime | n | Var(a) | Var(c) | Cov(a,c) [corr] | Var(c) − 2Cov | Var(e) | naive SE² | λ-SDᶜ² | IJ se² | Var(a)/nSE² | Var(e)/nSE² [r] | λ²/Var(e) | λ²/Var(a) | λ²/nSE² |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HR 1.00 n500 | s7c | 1361 | 0.01672 | 0.00021 | −0.00069 [−0.37] | 0.00159 | 0.01831 | 0.01946 | 0.01851 | 0.0649 | 0.859 | 0.941 [0.970] | 1.011 | 1.107 | 0.951 |
| HR 1.00 n500 | p30 | 1841 | 0.01625 | 0.00019 | −0.00023 [−0.13] | 0.00066 | 0.01691 | 0.01822 | 0.01751 | 0.0624 | 0.892 | 0.928 [0.963] | 1.035 | 1.077 | 0.961 |
| HR 1.00 n500 | p30sg | 1841 | 0.01704 | 0.00023 | −0.00022 [−0.11] | 0.00067 | 0.01770 | 0.01889 | 0.01786 | 0.0631 | 0.902 | 0.937 [0.968] | 1.009 | 1.048 | 0.945 |
| HR 1.00 n500 | nb20 B | 1879 | 0.01780 | 0.00027 | −0.00012 [−0.06] | 0.00052 | 0.01832 | 0.02019 | 0.01842 | 0.0650 | 0.882 | 0.907 [0.953] | 1.005 | 1.035 | 0.912 |
| HR 1.00 n1000 | p30 | 1909 | 0.00821 | 0.00007 | −0.00011 [−0.15] | 0.00030 | 0.00851 | 0.00876 | 0.00853 | 0.0310 | 0.937 | 0.971 [0.985] | 1.003 | 1.039 | 0.974 |
| HR 1.00 n1000 | p30sg | 1909 | 0.00839 | 0.00009 | −0.00010 [−0.11] | 0.00028 | 0.00867 | 0.00924 | 0.00878 | 0.0316 | 0.908 | 0.939 [0.969] | 1.012 | 1.046 | 0.950 |
| HR 1.00 n1000 | nb20 B | 1923 | 0.00919 | 0.00012 | −0.00006 [−0.06] | 0.00023 | 0.00942 | 0.01027 | 0.00923 | 0.0331 | 0.895 | 0.918 [0.958] | 0.979 | 1.004 | 0.899 |
| HR 1.50 n500 | p30 | 1999 | 0.01745 | 0.00021 | −0.00003 [−0.02] | 0.00028 | 0.01772 | 0.01764 | 0.01701 | 0.0617 | 0.989 | 1.005 [1.002] | 0.960 | 0.975 | 0.964 |
| HR 1.50 n500 | p30sg | 1999 | 0.01841 | 0.00024 | −0.00005 [−0.03] | 0.00035 | 0.01876 | 0.01859 | 0.01746 | 0.0631 | 0.990 | 1.009 [1.005] | 0.930 | 0.948 | 0.939 |
| HR 1.50 n500 | nb20 A | 1999 | 0.02047 | 0.00030 | −0.00007 [−0.03] | 0.00043 | 0.02091 | 0.02045 | 0.01833 | 0.0661 | 1.001 | 1.023 [1.011] | 0.877 | 0.895 | 0.897 |
| HR 1.50 n500 | nb20 B | 1999 | 0.02009 | 0.00027 | −0.00009 [−0.04] | 0.00045 | 0.02054 | 0.02025 | 0.01817 | 0.0658 | 0.992 | 1.014 [1.007] | 0.884 | 0.904 | 0.897 |
| HR 1.50 n1000 | nb20 B | 2000 | 0.00998 | 0.00013 | −0.00021 [−0.19] | 0.00056 | 0.01053 | 0.01051 | 0.00941 | 0.0345 | 0.949 | 1.002 [1.001] | 0.894 | 0.944 | 0.896 |
| HR 1.75 n500 | s7c | 1900 | 0.01882 | 0.00026 | −0.00040 [−0.18] | 0.00106 | 0.01988 | 0.01925 | 0.01857 | 0.0663 | 0.978 | 1.033 [1.016] | 0.934 | 0.987 | 0.965 |
| HR 1.75 n500 | p30 | 1999 | 0.01772 | 0.00021 | −0.00006 [−0.03] | 0.00032 | 0.01804 | 0.01756 | 0.01694 | 0.0616 | 1.009 | 1.028 [1.014] | 0.939 | 0.956 | 0.965 |
| HR 1.75 n500 | p30sg | 1999 | 0.01868 | 0.00024 | −0.00011 [−0.05] | 0.00047 | 0.01915 | 0.01864 | 0.01745 | 0.0634 | 1.002 | 1.027 [1.013] | 0.911 | 0.934 | 0.936 |
| HR 1.75 n500 | nb20 A | 1999 | 0.02082 | 0.00029 | −0.00020 [−0.08] | 0.00069 | 0.02151 | 0.02057 | 0.01840 | 0.0669 | 1.012 | 1.045 [1.022] | 0.856 | 0.884 | 0.894 |
| HR 1.75 n500 | nb20 B | 1999 | 0.02070 | 0.00028 | −0.00020 [−0.08] | 0.00067 | 0.02137 | 0.02042 | 0.01825 | 0.0666 | 1.014 | 1.046 [1.023] | 0.854 | 0.882 | 0.894 |

Means: a = −0.07 to −0.15 (the naive optimism on the complement), c = −0.05 to −0.09 (the correction removes 60–85% of it), e = −0.01 to −0.06. The IJ two-term se² is 3.1–3.7 × Var(e) everywhere.

**Reading of A1.** (1) **Var(c) is small**: 0.0001–0.0003, i.e. SD(c) 0.011–0.017 log-HR — the correction is nearly constant across datasets, a tenth of SD(a). (2) **Cov(a, c) is small and negative** (corr −0.02 to −0.08 at 31% prevalence; −0.18 / −0.37 in s7c, whose complement is 88% of the sample), so the correction adds Var(c) − 2Cov = 0.0003–0.0007 to the error variance: 1.5–3.3% of Var(e) at 31%, 5–9% in s7c. (3) **Var(a) ≈ naive SE²** in the harm cells at n = 500 (0.98–1.01) and below it in the null cells (0.86–0.90) and at n = 1000 (0.90–0.95): the naive complement estimate's realized error variance is what its SE says. (4) So **Var(e) = 1.00–1.05 × naive SE²** in the harm cells (r 1.00–1.02) and 0.91–0.97 in the null cells — the de-biased complement estimate is not appreciably noisier than the naive one. (5) **λ-SDᶜ² sits below the naive SE²**, by a margin that grows along the sequence: 0.965 (s7c, p30) → 0.936–0.939 (p30sg) → 0.894–0.897 (every nb20 cell, n = 500 and 1000 alike); against Var(e): 0.93–0.96 → 0.91–0.93 → 0.85–0.89.

## A2 — Stratified

**By p̂(Ĥ) tertile (nb20 cells; T1 = p̂ < 0.03–0.05, T3 = p̂ > 0.07–0.13):**

| Cell | tertile | n | \|Ĥ\|/\|H\| | Var(a) | Var(c) − 2Cov | Var(e) | naive SE² | λ² | Var(a)/nSE² | r | λ²/Var(e) | λ²/Var(a) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HR 1.50 n500 A | T1 / T2 / T3 | 669 / 664 / 666 | 1.01 / 0.81 / 0.68 | 0.0218 / 0.0194 / 0.0200 | 0.0000 / 0.0003 / 0.0006 | 0.0218 / 0.0197 / 0.0205 | 0.0221 / 0.0202 / 0.0190 | 0.0180 / 0.0183 / 0.0187 | 0.99 / 0.96 / 1.05 | 0.99 / 0.99 / 1.04 | **0.82** / 0.93 / 0.91 | 0.82 / 0.94 / 0.94 |
| HR 1.75 n500 A | T1 / T2 / T3 | 667 / 666 / 666 | 1.01 / 0.83 / 0.70 | 0.0216 / 0.0197 / 0.0202 | 0.0004 / 0.0006 / 0.0004 | 0.0220 / 0.0203 / 0.0206 | 0.0222 / 0.0203 / 0.0193 | 0.0180 / 0.0184 / 0.0189 | 0.98 / 0.97 / 1.05 | 1.00 / 1.00 / 1.03 | **0.82** / 0.90 / 0.92 | 0.83 / 0.93 / 0.94 |
| HR 1.50 n500 B | T1 / T2 / T3 | 673 / 660 / 666 | 0.98 / 0.80 / 0.66 | 0.0207 / 0.0198 / 0.0192 | 0.0004 / −0.0001 / 0.0006 | 0.0211 / 0.0197 / 0.0199 | 0.0219 / 0.0200 / 0.0188 | 0.0180 / 0.0181 / 0.0184 | 0.94 / 0.99 / 1.02 | 0.98 / 0.99 / 1.03 | **0.85** / 0.92 / 0.93 | 0.87 / 0.92 / 0.96 |
| HR 1.75 n500 B | T1 / T2 / T3 | 669 / 664 / 666 | 0.99 / 0.82 / 0.68 | 0.0218 / 0.0192 / 0.0203 | 0.0004 / 0.0004 / 0.0005 | 0.0223 / 0.0196 / 0.0208 | 0.0221 / 0.0202 / 0.0190 | 0.0179 / 0.0182 / 0.0186 | 0.99 / 0.95 / 1.07 | 1.00 / 0.98 / 1.05 | **0.81** / 0.93 / 0.90 | 0.82 / 0.95 / 0.92 |
| HR 1.50 n1000 B | T1 / T2 / T3 | 667 / 666 / 667 | 1.04 / 0.85 / 0.69 | 0.0114 / 0.0095 / 0.0085 | 0.0003 / 0.0004 / 0.0003 | 0.0117 / 0.0099 / 0.0088 | 0.0115 / 0.0104 / 0.0096 | 0.0094 / 0.0094 / 0.0095 | 0.99 / 0.91 / 0.88 | 1.01 / 0.97 / 0.96 | **0.80** / 0.95 / 1.08 | 0.82 / 0.99 / 1.12 |
| HR 1.00 n500 B | T1 / T2 / T3 | 628 / 626 / 625 | 0.93 / 0.72 / 0.62 | 0.0191 / 0.0160 / 0.0175 | 0.0002 / 0.0004 / 0.0008 | 0.0193 / 0.0163 / 0.0183 | 0.0218 / 0.0199 / 0.0189 | 0.0183 / 0.0184 / 0.0185 | 0.88 / 0.80 / 0.93 | 0.94 / 0.91 / 0.98 | 0.95 / 1.13 / 1.01 | 0.96 / 1.15 / 1.06 |
| HR 1.00 n1000 B | T1 / T2 / T3 | 643 / 639 / 641 | 0.99 / 0.74 / 0.60 | 0.0097 / 0.0087 / 0.0089 | 0.0002 / 0.0002 / 0.0000 | 0.0100 / 0.0089 / 0.0089 | 0.0115 / 0.0100 / 0.0093 | 0.0093 / 0.0092 / 0.0092 | 0.85 / 0.87 / 0.96 | 0.93 / 0.94 / 0.98 | 0.93 / 1.04 / 1.03 | 0.95 / 1.06 / 1.03 |

Low p̂ goes with large Ĥ (|Ĥ|/|H| ≈ 1.0 in T1 vs 0.7 in T3: the largest in-band candidate is the least re-selected), so the two stratifications are the same axis.

**By |Ĥ|/|H| tertile (the harm cells; the document has all 17 cells):**

| Cell | tertile bounds | Var(e) | naive SE² | λ² | r | λ²/nSE² | λ²/Var(e) |
|---|---|---|---|---|---|---|---|
| HR 1.75 n500 s7c | [0.66, 1.04] / [1.04, 1.20] / [1.21, 3.66] | 0.0207 / 0.0213 / 0.0174 | 0.0191 / 0.0191 / 0.0195 | 0.0187 / 0.0186 / 0.0185 | 1.04 / 1.06 / 0.94 | 0.98 / 0.97 / 0.95 | 0.90 / 0.87 / 1.07 |
| HR 1.75 n500 p30 | [0.35, 0.45] / [0.45, 0.56] / [0.56, 1.40] | 0.0183 / 0.0184 / 0.0175 | 0.0169 / 0.0173 / 0.0186 | 0.0169 / 0.0170 / 0.0169 | 1.04 / 1.03 / 0.97 | 1.00 / 0.99 / 0.91 | 0.92 / 0.93 / 0.97 |
| HR 1.75 n500 p30sg | [0.35, 0.54] / [0.54, 0.71] / [0.71, 1.40] | 0.0187 / 0.0189 / 0.0199 | 0.0171 / 0.0182 / 0.0206 | 0.0174 / 0.0174 / 0.0175 | 1.05 / 1.02 / 0.98 | 1.01 / 0.96 / 0.85 | 0.93 / 0.92 / 0.88 |
| HR 1.75 n500 nb20 A | [0.35, 0.71] / [0.71, 0.97] / [0.97, 1.93] | 0.0181 / 0.0219 / 0.0244 | 0.0181 / 0.0200 / 0.0236 | 0.0182 / 0.0184 / 0.0186 | 1.00 / 1.05 / 1.02 | 1.00 / 0.92 / 0.79 | 1.00 / 0.84 / **0.77** |
| HR 1.75 n500 nb20 B | [0.35, 0.69] / [0.69, 0.97] / [0.97, 1.73] | 0.0187 / 0.0213 / 0.0241 | 0.0180 / 0.0198 / 0.0235 | 0.0180 / 0.0183 / 0.0185 | 1.02 / 1.04 / 1.01 | 1.00 / 0.92 / 0.79 | 0.96 / 0.86 / **0.77** |
| HR 1.50 n500 nb20 A | [0.37, 0.69] / [0.69, 0.95] / [0.95, 2.09] | 0.0190 / 0.0206 / 0.0232 | 0.0181 / 0.0198 / 0.0235 | 0.0181 / 0.0183 / 0.0186 | 1.02 / 1.02 / 0.99 | 1.00 / 0.93 / 0.79 | 0.96 / 0.89 / **0.80** |
| HR 1.50 n1000 nb20 B | [0.31, 0.70] / [0.70, 1.02] / [1.02, 1.66] | 0.0091 / 0.0112 / 0.0113 | 0.0089 / 0.0105 / 0.0121 | 0.0091 / 0.0094 / 0.0098 | 1.01 / 1.03 / 0.97 | 1.02 / 0.89 / 0.81 | 1.00 / 0.84 / 0.86 |

**Reading of A2.** Within a cell the error ratio r stays at 1.0 in every tertile; the naive SE² and Var(e) both grow with |Ĥ| (the complement shrinks: from 0.018 to 0.024 across the nb20 tertiles), and **λ-SDᶜ² does not follow** (0.0182 → 0.0186, flat). Per replicate, `fld_Hc_se` tracks `nv_Hc_se` with slope 0.18–0.26 and correlation 0.34–0.47 in the p30sg / nb20 cells (SD of `fld_Hc_se` 0.005–0.006 vs 0.008–0.010 for `nv_Hc_se`; the ratio `fld_Hc_se`/`nv_Hc_se` is 1.00 in the lowest |Ĥ| tertile and 0.89 (nb20) / 0.93 (p30sg) in the highest), against slope 0.36 / 0.61 in p30 / s7c, whose Ĥ varies less. The shortfall is concentrated in the replicates whose Ĥ reaches or exceeds the planted size — exactly the ones the wider band was meant to produce.

## A3 — Across the regime sequence

HR 1.75 n500 (excess_a = Var(a) − naive SE²; excess_c = Var(c) − 2Cov(a,c); shortfall = Var(e) − λ²):

| regime | \|Ĥ\| | p̂ | Var(a) | Var(c) | Cov | Var(e) | naive SE² | λ² | excess_a | excess_c | Var(e) − nSE² | λ² − nSE² | shortfall | r | λ²/Var(e) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| s7c (12%, maxeffCons) | 72 | — | 0.01882 | 0.00026 | −0.00040 | 0.01988 | 0.01925 | 0.01857 | −0.00043 | 0.00106 | 0.00064 | −0.00067 | 0.00131 | 1.016 | 0.934 |
| p30 (31%, maxeffCons) | 80 | — | 0.01772 | 0.00021 | −0.00006 | 0.01804 | 0.01756 | 0.01694 | 0.00016 | 0.00032 | 0.00049 | −0.00062 | 0.00110 | 1.014 | 0.939 |
| p30sg (ε 0.10, J 10) | 99 | — | 0.01868 | 0.00024 | −0.00011 | 0.01915 | 0.01864 | 0.01745 | 0.00004 | 0.00047 | 0.00050 | −0.00119 | 0.00170 | 1.013 | 0.911 |
| nb20 A (ε 0.20, J 10) | 129 | 0.119 | 0.02082 | 0.00029 | −0.00020 | 0.02151 | 0.02057 | 0.01840 | 0.00025 | 0.00069 | 0.00093 | −0.00217 | 0.00311 | 1.022 | 0.856 |
| nb20 B (ε 0.20, J 20) | 126 | 0.097 | 0.02070 | 0.00028 | −0.00020 | 0.02137 | 0.02042 | 0.01825 | 0.00028 | 0.00067 | 0.00095 | −0.00217 | 0.00312 | 1.023 | 0.854 |

HR 1.50 n500: shortfall 0.00072 (p30) → 0.00131 (p30sg) → 0.00257 / 0.00238 (nb20 A / B); λ²/Var(e) 0.960 → 0.930 → 0.877 / 0.884; r 1.002 → 1.005 → 1.011 / 1.007; excess_c 0.00028 → 0.00035 → 0.00043 / 0.00045; λ² − naive SE² −0.00064 → −0.00113 → −0.00211 / −0.00208. HR 1.00 n500: Var(e) < naive SE² throughout (r 0.95–0.97) and the shortfall is ≤ 0.0006 in magnitude (λ²/Var(e) 1.01–1.04); n = 1000: harm 1.5 shortfall 0.00112 (λ²/Var(e) 0.894), null 0.98–1.01.

**Which component grows.** Along the sequence at HR 1.75 the shortfall grows ×2.4 (0.0013 → 0.0031) while (a) Var(a) − naive SE² stays within ±0.0003, (b) Var(c) − 2Cov grows 0.0003 → 0.0007 (29% → 28% → 22% of the shortfall at 31% prevalence — a minor share that shrinks; 81% in s7c, where corr(a, c) is −0.18), and (c) **λ² − naive SE² grows −0.0006 → −0.0022** and accounts for 55–70% of the shortfall at 31% prevalence (p30 56%, p30sg 70%, nb20 70%; s7c 51%). What grows is the gap between the field's SD and the naive SE, not any component of the data's error variance.

## A4 — Reading (for the Linux chat's decision)

**(i) Var(a) ≈ naive SE² holds in every harm cell (0.98–1.01).** The excess of Var(e) over Var(a) is Var(c) − 2Cov(a, c) = 0.0003–0.0007 at 31% prevalence — the correction's data-to-data variability is real but small (SD(c) 0.015–0.017 against SD(a) 0.14) and it is not what the field misses: the field's λ-SDᶜ² falls short of Var(e) by 0.0026–0.0031 in the nb20 harm cells, of which the correction's term is 0.0004–0.0007 and the remaining 0.0021–0.0022 is λ-SDᶜ² lying **below the naive SE² itself** (0.894–0.897 of it in every nb20 harm cell, 0.899–0.912 in its null cells; 0.936–0.939 in p30sg's harm cells; 0.964–0.965 in p30 / s7c). So the premise of (i) holds but its conclusion does not: the field under-simulates neither the correction's variability (that would be a 2–3% effect) nor the naive noise as such; it under-scales the complement's noise for the replicates whose Ĥ is large. A2 shows where: `fld_Hc_se` is flat across replicates (0.134–0.136 across the |Ĥ| tertiles) while `nv_Hc_se` and SD(e) grow with |Ĥ| (0.134 → 0.154). The complement field draws ζᶜ = Bcᵀξ at the **re-selected** winner G_r of each outer draw, so Λ*ᶜ carries the complement noise of whichever candidates the field re-selects — a family-average complement scale — not the selected complement's own; when the selected Ĥ is larger than the field's typical re-selection (the low-p̂ replicates: p̂ < 0.04, |Ĥ|/|H| ≈ 1.0, λ²/Var(e) 0.80–0.85), the field's SD is the smaller one. **The analysis-time objects a diagnosis or repair would need**, all inside `.fs_mr_field_complement()` with no extra draws: the two pieces of Λ*ᶜ_r recorded separately — ζᶜ_{r,G_r} (`Zo_c[G, r]`) and the inner mean m̂ᶜ(v_r) (`mean(Zi_c[cbind(wi, ok_in)])`) — giving the field's own Var(ζᶜ_G), Var(m̂ᶜ) and Cov(ζᶜ_G, m̂ᶜ) as the analogues of Var(a), Var(c), Cov(a, c); and the per-candidate complement noise scale diag(Bcᵀ Bc)[g] for the outer winners against that of the selected complement (`Bc[, sel]`), which is what an analysis can compare to its own naive SE². A repair that scales Λ*ᶜ by the selected complement's SE relative to the field's mean winner SE would be a document-level or `R/` change for the Linux chat to propose; this record does not.

**(ii) Var(a) does not exceed naive SE²** in any cell (it is below it in the null cells and at n = 1000). The complement's identity varying with the winner shows up **in the target, not the error**: SD(log β(Ĥᶜ)) grows 0.056 → 0.081 → 0.099–0.101 along the sequence and is what lifts the reported SD(β̃ᶜ)/naive SE from 1.03 to 1.15 (A0). Against Var(a), λ-SDᶜ² falls short by 1.3–4.4% (s7c, p30), 5–7% (p30sg) and 10–12% (nb20 n = 500; 6% at harm 1.5 n1000); against SD(a), 0.7–2.2%, 2.6–3.4%, 5–6%.

**(iii) The shortfall is not a stable fraction of the regime diagnostic.** In error terms the diagnostic barely moves (r² 1.03 → 1.03 → 1.03 → 1.045 at HR 1.75; 1.005 → 1.009 → 1.02 at HR 1.50) while the shortfall 1 − λ²/Var(e) grows 0.066 → 0.061 → 0.089 → 0.144–0.146 and 0.040 → 0.070 → 0.116–0.123; in the reported (marginal) terms the diagnostic grows 1.03 → 1.15 and the shortfall tracks it closely — because both are driven by the same thing, the widening spread of Ĥ across replicates, which moves the target and starves the field of the selected complement's scale on the large-Ĥ replicates. The rule as documented ("field bound with its shortfall stated against SD(β̃ᶜ)/naive SE and p̂") therefore states the shortfall against a quantity that mixes the target's spread with the error; the same rule stated against SD(e)/naive SE would read "error ratio 1.00–1.02 in every harm cell, field SD 0.92–0.97 of the error SD, falling with the spread of |Ĥ|/|H| (tertile ratios 1.00 / 0.92 / 0.79 of the naive SE at ε = 0.20)". The IJ two-term at 3.1–3.7 × Var(e) remains the conservative option whether the diagnostic is read either way.

## Side observations (not tasks)

1. **Every constructions table so far reports the marginal SD.** `fs_sim_bias_coverage()` (`sd_emp = sd(log est)`) and the template's `.sd_log()` are both marginal; "bias in SD units" and "SE/SD" in the p30 / p30sg / nb20 reports are on that scale, and the reports' "SD(β̃ᶜ)/naive SE" is therefore Var(e) + Var(target) + 2Cov, not the error. On the harm block the same applies (β(Ĥ) varies with Ĥ too). Coverage numbers are unaffected. Whether the tables should carry the error SD alongside is Larry's call; the per-replicate columns support both.
2. The complement's `selection_bias` / `fixed_bias` split is returned by the gate for the harm block only and recorded for neither block; a template-level recorder addition would need the gate to expose the complement's two components first (an `R/` change).
3. Field draw usage: `fld_Hc_lam_mean` (the field's own estimate of the complement correction) is −0.007 to −0.014, while the two-term correction c it inverts around is −0.05 to −0.08; corr(c, `fld_Hc_lam_mean`) = 0.60–0.68 across replicates. Recorded here as a reference for the Linux chat; not analysed further.

Nothing blocked; no task proposed by this record. The decision on a repair proposal is the Linux chat's, from A4.
