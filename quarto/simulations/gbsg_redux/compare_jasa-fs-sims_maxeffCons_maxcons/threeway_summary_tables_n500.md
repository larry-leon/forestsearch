---
title: "Three-way summary: t1_t2_legacy vs maxeffCons vs maxcons (M1, h10, n = 500)"
bibliography: []
---

All three bundles: `n_sample` 500, `n_sims` 500, `nb_boots` 300, draws 5000,
`seed_base` 8316951, `consistency_method = "resample"`, `forestsearch 0.2.0`,
identical `truth` (HR 0.685, marg_H 1.000, marg_Hc 0.657, cde_H 1.000, cde_Hc 0.585).
`sg_focus`: `"eff"` / `"maxeffCons"` / `"maxcons"`.

Subgroup and estimation rows are over the detected replicates.

# A. Detection and subgroup identification

| | t1_t2_legacy | maxeffCons | maxcons |
|---|---|---|---|
| Detections / 500 | 346 | 346 | 346 |
| Detection rate | 69.2% | 69.2% | 69.2% |
| mean \|Ĥ\| | 82.3 | 78.1 | 87.5 |
| SD \|Ĥ\| | 21.8 | 17.8 | 28.1 |
| median \|Ĥ\| | 76 | 73 | 79 |
| min \|Ĥ\| | 61 | 61 | 61 |
| max \|Ĥ\| | 184 | 167 | 220 |
| mean \|Ĥ\| as % ITT | 16.5% | 15.6% | 17.5% |
| mean true Q | 61.9 | 61.9 | 61.9 |
| Sensitivity | 0.390 | 0.368 | 0.405 |
| Specificity | 0.867 | 0.874 | 0.858 |
| PPV | 0.315 | 0.309 | 0.309 |
| NPV | 0.909 | 0.907 | 0.910 |
| MR harm flags | 145 | 109 | 134 |
| MR returned estimate | 346 / 346 | 346 / 346 | 346 / 346 |
| mean β(Ĥ) | 0.728 | 0.726 | 0.725 |

# B. Estimation and coverage — harm block Ĥ

Avg = mean point estimate; Len = median interval length; Bias% and Cov are against
β(Ĥ).

| Estimator | Metric | t1_t2_legacy | maxeffCons | maxcons |
|---|---|---|---|---|
| oracle | Avg | 1.135 | 1.135 | 1.135 |
| | Len | 1.441 | 1.441 | 1.441 |
| | Bias% | 54.9 | 55.2 | 55.6 |
| | Cov | 0.835 | 0.838 | 0.827 |
| naive | Avg | 1.629 | 1.643 | 1.611 |
| | Len | 2.082 | 2.157 | 1.964 |
| | Bias% | 128.0 | 130.9 | 125.9 |
| | Cov | 0.145 | 0.153 | 0.130 |
| FB | Avg | 1.222 | 1.229 | 1.223 |
| | Len | 1.659 | 1.749 | 1.648 |
| | Bias% | 71.0 | 72.6 | 71.3 |
| | Cov | 0.801 | 0.806 | 0.783 |
| MR | Avg | 1.003 | 0.946 | 0.992 |
| | Len | 1.473 | 1.536 | 1.395 |
| | Bias% | 40.0 | 32.8 | 38.7 |
| | Cov | 0.968 | 0.986 | 0.968 |

The oracle refits on the true region H, so its Bias% and Cov columns are computed
against a different region than it estimates.

# C. Estimation and coverage — complement block Ĥᶜ

| Estimator | Metric | t1_t2_legacy | maxeffCons | maxcons |
|---|---|---|---|---|
| oracle | Avg | 0.660 | 0.660 | 0.660 |
| | Len | 0.360 | 0.360 | 0.360 |
| | Bias% | 0.1 | −0.1 | 0.1 |
| | Cov | 0.960 | 0.960 | 0.960 |
| naive | Avg | 0.605 | 0.609 | 0.598 |
| | Len | 0.335 | 0.336 | 0.336 |
| | Bias% | −8.2 | −7.6 | −9.2 |
| | Cov | 0.916 | 0.928 | 0.908 |
| FB | Avg | 0.631 | 0.636 | 0.626 |
| | Len | 0.490 | 0.487 | 0.485 |
| | Bias% | −4.2 | −3.6 | −5.0 |
| | Cov | 0.983 | 0.988 | 0.986 |
| MR | Avg | 0.665 | 0.655 | 0.658 |
| | Len | 0.672 | 0.679 | 0.667 |
| | Bias% | 0.9 | −0.7 | −0.2 |
| | Cov | 1.000 | 1.000 | 1.000 |

# D. Timings (seconds per replicate)

Legacy columns `t1_secs`/`t2_secs`; new arms `fb_secs`/`fit_mr_secs` — same measured
quantities. FB rows over detected replicates, fit+MR rows over all 500.

| | t1_t2_legacy | maxeffCons | maxcons |
|---|---|---|---|
| FB mean | 166.4 | 136.2 | 167.5 |
| FB median | 117.2 | 135.2 | 185.4 |
| FB min | 100.8 | 105.2 | 113.3 |
| FB max | 397.6 | 165.3 | 192.1 |
| fit+MR mean | 12.78 | 5.01 | 12.66 |
| fit+MR median | 11.38 | 5.41 | 10.99 |
| fit+MR min | 3.37 | 2.14 | 3.33 |
| fit+MR max | 23.64 | 6.97 | 23.47 |
| batches pooled | 4 | 3 | 3 |

Timings are wall-clock under 115 workers and vary by batch within each arm (legacy FB
batch means range 113–372 s; maxcons 137–188 s; maxeffCons 132–143 s), so
cross-arm differences reflect render conditions as well as any code effect.
