# REPORT — ACTG175 binary (OR) campaigns `orfs`, `orgrf`, `ordina`: Stage 3 summarization of the 18 committed cells (2026-09-28)

Task: `dev/tasks/TASK_actg175_or_stage3_2026-09-28.md` (committed `70d8a488`). Branch `feature/glm-extension`, HEAD at the
start `70d8a488` (task doc) on `70ef2d21`. Machine `pop-os`, R 4.6.1.

**Summarization only.** No simulation was run; no bundle was written, moved or deleted. R read the 18 committed combined
bundles and wrote `or_metrics.csv`, `COLUMNS_or.md` and the render `summary_actg175_or.html`. The md5 of every bundle read
was taken before and after (`§6`); all 18 unchanged. Untracked files present at the start
(`mr_or_harm/…_redes_d5000/`, `…_relaunch_d5000/`, `smoke_redes.html`, `smoke_relaunch.html`,
`gbsg_020/scripts_dinamr/logs/nullmr_findings.err`) were left alone. The closeout deletion of the smoke and calibration
bundles (catalog §3.5) was not done (not in this task).

**Every GRF and DINA coverage figure below is coverage of the estimand conditional on the proposed family**
(`current_status.md` §5). FS's family is the prespecified cut grid. The three identifiers are compared descriptively,
never ranked.

---

## 1. Design and settings (cited, not re-derived)

Template of record `sim_fs_mr_field_or_template.qmd`; catalog `current_status.md` §2 (pin `d356bd8a`); inventory
`dev/reports/REPORT_glm_campaign_inventory_2026-09-28.md` §2.

- **Data and outcome.** ACTG175 arms 1 (ZDV+ddI) vs 3 (ddI); `y_neg = 1 − 1{cd420 > cd40}`; `outcome_type "binary"`,
  `effect_measure "OR"` (`:316–317`); `adverse_outcome = TRUE` in the analysis, so **OR > 1 is harm** (`:293`); DGM
  constructed with `adverse_outcome = FALSE` (`:401`, intentional).
- **Planted region and design points.** H = {wtkg > q} ∩ {cd40 > q}, `sg_quantile 0.62850` (`:307`, `:311–313`),
  prevalence(H) 14.917 % (design of record since 2026-09-18, `REPORT_binary_redesign_2026-09-18.md`); `n_super 100000`
  (`:296`); `target_effect = FS_OR_TARGET` (`:207`, `:395`) = 0.75 (protective), 1.0 (borderline null), 1.5 (harm); n = 500
  and 2000 (`:210`). Six cells per identifier, 18 in all.
- **Replicates and seeds.** 1,000 per cell, one batch over `sim_id` 1–1000 (`:146`); seeds from a pre-generated table,
  `seed_base 8316951` (`:164–166`), indexed by global `sim_id`, so identical across identifiers and machines (verified per
  replicate, `REPORT_binary_stage2_2026-09-19.md` §4).
- **Rule and thresholds.** `sg_focus "effMaxSG"` (`:179`), `effect_neighborhood 0.20` (`:186`), `selection_rule
  "neighborhood"` (`:268`); c1 = 0.90, c2 = 0.80, p⋆ = 0.90 on the OR scale (`:277–279`); `fs.splits 500`, `maxk 2`,
  `n.min 60`, `d0/d1.min 10` (`:280`); `consistency_method "resample"` (`:267`). Identifiers: `consistency` (`orfs`),
  `grf` (frontier / effect / depth 2 / `dmin.grf 0`, `:287–290`), `dina` (`dina_args = list()`, `effect`, `:291–292`).
- **MR and bounds.** `ci_method "field"` (`:331`), 5,000 draws (`:149`), complement field on (`:336`),
  `field_scale_complement "selected"` (`:339`), `ij_residual "two_term"` (`:341`), `return_reselection TRUE` (`:343`);
  FB never run (`:148`, `:354`). Recorded: unadjusted (`nv_*`), oracle on the true region (`or_*`), IJ two-term (`mr_*`),
  field on Ĥ (`fld_H_*`), field and field-s on Ĥᶜ (`fld_Hc_*`, `fld_Hc_*_s`), the joint pairs (`fld_joint_*`,
  `fld_joint_s_*`). Every bundle's `meta` confirms: `ci_method field`, `field_scale_complement selected`,
  `effect_neighborhood 0.2`, `effect_threshold 0.9`, `consistency_threshold 0.8`, `pconsistency 0.9`, `fb_mode none`,
  `feas_feasible TRUE`, `feas_override FALSE`, `sg_quantile 0.6285`, `n_sims 1000`.

## 2. The design's truths (Step 3; read from `truth` / `meta` of each design point's bundles, identical across identifiers and n)

| design point | `target_or_h` | prevalence(H) `truth$prevalence_Q` | marginal OR in H `truth$marg_H` (θ†) | marginal OR in Hᶜ `truth$marg_Hc` | whole-trial marginal OR `truth$or_causal` (= `truth$effect_ITT`) | CDE in H `truth$cde_H` (θ‡) | CDE in Hᶜ `truth$cde_Hc` | calibrated `beta_inter` |
|---|---|---|---|---|---|---|---|---|
| OR 0.75 | 0.75 | 0.14917 | 0.7499999940 | 0.6564077470 | 0.6722232577 | 0.7345268290 | 0.6313905111 | 0.1513019741 |
| OR 1.0 | 1.00 | 0.14917 | 0.9999999999 | 0.6564077470 | 0.7043635080 | 0.9999999999 | 0.6313905111 | 0.4598307312 |
| OR 1.5 | 1.50 | 0.14917 | 1.5000000002 | 0.6564077470 | 0.7501835133 | 1.5414142152 | 0.6313905111 | 0.8925310479 |

Every field was present in every bundle; nothing was computed. The complement's marginal OR and CDE are the same at
all three design points (the complement inherits the fitted ACTG175 effect; catalog §2). The committed study's
9.632 %-prevalence truths (θ†(Hᶜ) 0.6560, overall 0.6666) differ from these and are not a comparator (catalog §2).

## 3. What was done (Steps 1, 2, 4)

**Step 1 — `summary_actg175_or.qmd`, six edits, no logic change** (`git diff` in the commit):
1. `:64` `glob_new` default `"combined_1_2000"` → `"combined_1_1000"` (the committed layout).
2. `:441` caption "2,000 replicates per cell x identifier" → "1,000".
3. `:541` and `:550` (study-comparison prose and caption) "2,000 replicates" → "1,000".
4. `:597` COLUMNS header gains "Stage 3 run under TASK_actg175_or_stage3_2026-09-28".
5. `:602` COLUMNS `design` line: "prevalence 9.632 % … θ†(Ĥᶜ) 0.6560" → "the design of record's planted prevalence
   14.917 % (`sg_quantile` 0.62850) … θ†(Ĥᶜ) 0.6564 / θ‡(Ĥᶜ) 0.6314" (the values the bundles carry, §2).
The document's computations, conventions (non-convergent-fit masking per estimator, log-OR working scale, the OR ladder)
and output paths are unchanged.

**Step 2 — extract.** `or_metrics.csv`: 5,814 rows; columns `campaign, identifier, design, cell, block, estimator,
metric, tau, value, mc_se, wilson_lo, wilson_hi, n, commit` — `md_dina_metrics.csv`'s schema (itself
`md_field_metrics.csv` + `identifier`) **plus `design`**, as `COLUMNS_or.md` states. Every one of the 18 identifier × cell
combinations is present (1,932 rows per campaign; 18 quoted rows from the committed study grid under campaign
`maxeffCons_actg175_or075_seedtab_s1000`). Bias is `bias_log` = mean(log est − log target) with the target β(Ĥ) / β(Ĥᶜ)
per replicate (θ† for the oracle row), exactly as the continuous extracts define `bias_md` on their scale; `bias_or` and
`bias_sd_units` beside it. `n_nonconvergent_fits` is **0 in every block × estimator of every cell**. `COLUMNS_or.md`
defines every column and metric.

**Step 4 — render.** `quarto render summary_actg175_or.qmd` (quarto 1.9.38, 24 s): "read 18 of 18 cell x identifier
bundles; 5796 extract rows" (+ 18 study rows appended in the extract chunk). The render also wrote
`fig_or_bias_coverage_display_{H,Hc}.png` (the summary's `ggsave` lines `:575`, `:582`); the same figures are embedded in
the html, and the task's post-conditions list no figure files, so **the two PNGs were removed after the render** (they
are this task's own outputs, not bundles). `current_status.md` edited in place (§7).

## 4. Comparative tables (one per block; from `or_metrics.csv`; rows identifier × cell; declared replicates)

Conventions: `declaration` over 1,000 replicates; `n` = declared replicates (every estimator converged on all of them);
bias on the log-OR scale (Monte Carlo SE in parentheses) against β(Ĥ) / β(Ĥᶜ); "SD units" = `bias_log` / SD(log est);
one-sided coverage on the exposed side with Wilson 95 % interval; two-sided coverage of each estimator's own interval, the
oracle's against θ†. GRF and DINA rows conditional on the proposed family.

### 4.1 Ĥ (harm block; one-sided = the field LOWER bound)

| identifier | cell | declaration | n | bias log-OR: naive | IJ two-term | field | field bias, SD units | field one-sided LOWER coverage (Wilson) | two-sided: naive | oracle (θ†) | IJ two-term | field |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| FS | OR 0.75 n500 | 0.781 | 781 | 1.158 (0.009) | 0.287 (0.011) | 0.130 (0.012) | 0.40 | 0.974 (0.961, 0.983) | 0.010 | 0.948 | 0.992 | 0.988 |
| FS | OR 0.75 n2000 | 0.910 | 910 | 1.055 (0.009) | 0.216 (0.010) | 0.066 (0.011) | 0.20 | 0.978 (0.966, 0.986) | 0.022 | 0.931 | 0.989 | 0.992 |
| FS | OR 1.0 n500 | 0.852 | 852 | 1.114 (0.009) | 0.248 (0.011) | 0.092 (0.012) | 0.26 | 0.978 (0.965, 0.986) | 0.110 | 0.964 | 0.995 | 0.991 |
| FS | OR 1.0 n2000 | 0.961 | 961 | 1.018 (0.010) | 0.189 (0.011) | 0.042 (0.012) | 0.11 | 0.974 (0.962, 0.982) | 0.117 | 0.946 | 0.984 | 0.975 |
| FS | OR 1.5 n500 | 0.931 | 931 | 1.032 (0.010) | 0.181 (0.012) | 0.030 (0.013) | 0.07 | 0.977 (0.966, 0.985) | 0.284 | 0.954 | 0.995 | 0.983 |
| FS | OR 1.5 n2000 | 0.998 | 998 | 0.871 (0.013) | 0.082 (0.015) | −0.048 (0.016) | −0.10 | 0.964 (0.950, 0.974) | 0.439 | 0.950 | 0.951 | 0.909 |
| GRF | OR 0.75 n500 | 0.999 | 999 | 1.135 (0.012) | 0.297 (0.014) | 0.136 (0.014) | 0.31 | 0.931 (0.913, 0.945) | 0.272 | 0.950 | 0.992 | 0.962 |
| GRF | OR 0.75 n2000 | 0.999 | 999 | 1.009 (0.010) | 0.252 (0.012) | 0.105 (0.012) | 0.28 | 0.945 (0.929, 0.957) | 0.179 | 0.931 | 0.980 | 0.969 |
| GRF | OR 1.0 n500 | 0.999 | 999 | 1.120 (0.012) | 0.289 (0.014) | 0.130 (0.014) | 0.29 | 0.936 (0.919, 0.950) | 0.302 | 0.962 | 0.991 | 0.958 |
| GRF | OR 1.0 n2000 | 1.000 | 1000 | 0.980 (0.011) | 0.234 (0.013) | 0.092 (0.013) | 0.24 | 0.945 (0.929, 0.958) | 0.212 | 0.943 | 0.979 | 0.963 |
| GRF | OR 1.5 n500 | 1.000 | 1000 | 1.076 (0.013) | 0.264 (0.014) | 0.109 (0.015) | 0.24 | 0.941 (0.925, 0.954) | 0.343 | 0.952 | 0.995 | 0.958 |
| GRF | OR 1.5 n2000 | 1.000 | 1000 | 0.809 (0.014) | 0.107 (0.015) | −0.017 (0.015) | −0.04 | 0.948 (0.932, 0.960) | 0.448 | 0.948 | 0.947 | 0.907 |
| DINA | OR 0.75 n500 | 0.944 | 944 | 1.196 (0.014) | 0.504 (0.013) | 0.358 (0.014) | 0.84 | 0.870 (0.847, 0.890) | 0.244 | 0.951 | 0.988 | 0.932 |
| DINA | OR 0.75 n2000 | 0.870 | 870 | 1.060 (0.014) | 0.461 (0.014) | 0.333 (0.015) | 0.76 | 0.869 (0.845, 0.890) | 0.298 | 0.933 | 0.975 | 0.918 |
| DINA | OR 1.0 n500 | 0.968 | 968 | 1.179 (0.014) | 0.475 (0.013) | 0.328 (0.014) | 0.77 | 0.881 (0.859, 0.900) | 0.272 | 0.963 | 0.994 | 0.940 |
| DINA | OR 1.0 n2000 | 0.934 | 934 | 1.069 (0.014) | 0.429 (0.014) | 0.298 (0.015) | 0.69 | 0.878 (0.855, 0.897) | 0.273 | 0.945 | 0.974 | 0.931 |
| DINA | OR 1.5 n500 | 0.979 | 979 | 1.140 (0.014) | 0.413 (0.014) | 0.264 (0.014) | 0.60 | 0.903 (0.883, 0.920) | 0.325 | 0.953 | 0.989 | 0.940 |
| DINA | OR 1.5 n2000 | 0.987 | 987 | 1.002 (0.016) | 0.307 (0.016) | 0.172 (0.016) | 0.37 | 0.903 (0.883, 0.920) | 0.402 | 0.949 | 0.969 | 0.938 |

### 4.2 Ĥᶜ (complement block; one-sided = the field-s UPPER bound, the evaluated construction; the unstudentized field beside it)

| identifier | cell | declaration | bias log-OR: naive | IJ two-term | field-s | field-s bias, SD units | field-s one-sided UPPER coverage (Wilson) | field (unstud.) one-sided UPPER | two-sided: naive | oracle (θ†) | IJ two-term | field-s |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| FS | OR 0.75 n500 | 0.781 | −0.199 (0.006) | −0.061 (0.006) | −0.034 (0.006) | −0.19 | 0.950 (0.932, 0.963) | 0.941 | 0.899 | 0.971 | 1.000 | 0.969 |
| FS | OR 0.75 n2000 | 0.910 | −0.057 (0.003) | −0.020 (0.003) | −0.012 (0.003) | −0.13 | 0.942 (0.925, 0.955) | 0.938 | 0.913 | 0.955 | 1.000 | 0.948 |
| FS | OR 1.0 n500 | 0.852 | −0.210 (0.006) | −0.072 (0.006) | −0.044 (0.006) | −0.24 | 0.948 (0.931, 0.961) | 0.928 | 0.876 | 0.969 | 1.000 | 0.968 |
| FS | OR 1.0 n2000 | 0.961 | −0.061 (0.003) | −0.025 (0.003) | −0.018 (0.003) | −0.19 | 0.938 (0.920, 0.951) | 0.939 | 0.924 | 0.950 | 1.000 | 0.941 |
| FS | OR 1.5 n500 | 0.931 | −0.214 (0.006) | −0.078 (0.006) | −0.051 (0.006) | −0.27 | 0.928 (0.910, 0.943) | 0.921 | 0.861 | 0.963 | 1.000 | 0.959 |
| FS | OR 1.5 n2000 | 0.998 | −0.048 (0.003) | −0.015 (0.003) | −0.008 (0.003) | −0.08 | 0.945 (0.929, 0.957) | 0.943 | 0.924 | 0.947 | 1.000 | 0.949 |
| GRF | OR 0.75 n500 | 0.999 | −0.212 (0.006) | −0.068 (0.006) | −0.036 (0.006) | −0.18 | 0.934 (0.917, 0.948) | 0.931 | 0.845 | 0.958 | 1.000 | 0.957 |
| GRF | OR 0.75 n2000 | 0.999 | −0.060 (0.003) | −0.021 (0.003) | −0.011 (0.003) | −0.11 | 0.934 (0.917, 0.948) | 0.931 | 0.895 | 0.947 | 0.999 | 0.944 |
| GRF | OR 1.0 n500 | 0.999 | −0.209 (0.006) | −0.067 (0.006) | −0.034 (0.006) | −0.17 | 0.939 (0.922, 0.952) | 0.934 | 0.847 | 0.958 | 0.999 | 0.958 |
| GRF | OR 1.0 n2000 | 1.000 | −0.062 (0.003) | −0.023 (0.003) | −0.013 (0.003) | −0.14 | 0.931 (0.914, 0.945) | 0.931 | 0.897 | 0.946 | 1.000 | 0.939 |
| GRF | OR 1.5 n500 | 1.000 | −0.204 (0.006) | −0.065 (0.006) | −0.034 (0.006) | −0.17 | 0.941 (0.925, 0.954) | 0.942 | 0.857 | 0.957 | 1.000 | 0.947 |
| GRF | OR 1.5 n2000 | 1.000 | −0.050 (0.003) | −0.013 (0.003) | −0.004 (0.003) | −0.04 | 0.944 (0.928, 0.957) | 0.943 | 0.918 | 0.946 | 1.000 | 0.942 |
| DINA | OR 0.75 n500 | 0.944 | −0.196 (0.006) | −0.086 (0.006) | −0.062 (0.007) | −0.31 | 0.916 (0.897, 0.932) | 0.913 | 0.886 | 0.962 | 1.000 | 0.936 |
| DINA | OR 0.75 n2000 | 0.870 | −0.038 (0.003) | −0.015 (0.003) | −0.009 (0.003) | −0.10 | 0.940 (0.922, 0.954) | 0.939 | 0.939 | 0.952 | 1.000 | 0.945 |
| DINA | OR 1.0 n500 | 0.968 | −0.202 (0.006) | −0.089 (0.006) | −0.063 (0.007) | −0.31 | 0.908 (0.888, 0.925) | 0.905 | 0.862 | 0.962 | 1.000 | 0.940 |
| DINA | OR 1.0 n2000 | 0.934 | −0.048 (0.003) | −0.021 (0.003) | −0.015 (0.003) | −0.16 | 0.928 (0.910, 0.943) | 0.926 | 0.927 | 0.946 | 1.000 | 0.940 |
| DINA | OR 1.5 n500 | 0.979 | −0.198 (0.006) | −0.080 (0.006) | −0.054 (0.006) | −0.27 | 0.919 (0.901, 0.935) | 0.917 | 0.880 | 0.961 | 1.000 | 0.947 |
| DINA | OR 1.5 n2000 | 0.987 | −0.043 (0.003) | −0.015 (0.003) | −0.008 (0.003) | −0.09 | 0.939 (0.923, 0.952) | 0.936 | 0.933 | 0.947 | 1.000 | 0.945 |

### 4.3 Joint pair (Ĥ lower, Ĥᶜ upper) on (β(Ĥ), β(Ĥᶜ)); every declared replicate carries both bounds (`share_both_bounds` = 1.000 in all 18 cells)

| identifier | cell | declaration | Bonferroni field-s joint coverage (Wilson) | n | its Ĥ side | its Ĥᶜ side | separate 95 % pair (field lower, field-s upper) | unstudentized Bonferroni |
|---|---|---|---|---|---|---|---|---|
| FS | OR 0.75 n500 | 0.781 | 0.960 (0.944, 0.972) | 781 | 0.988 | 0.972 | 0.924 | 0.953 |
| FS | OR 0.75 n2000 | 0.910 | 0.963 (0.948, 0.973) | 910 | 0.993 | 0.969 | 0.921 | 0.960 |
| FS | OR 1.0 n500 | 0.852 | 0.964 (0.949, 0.974) | 852 | 0.991 | 0.973 | 0.926 | 0.957 |
| FS | OR 1.0 n2000 | 0.961 | 0.955 (0.940, 0.967) | 961 | 0.990 | 0.966 | 0.914 | 0.954 |
| FS | OR 1.5 n500 | 0.931 | 0.957 (0.942, 0.968) | 931 | 0.990 | 0.967 | 0.905 | 0.946 |
| FS | OR 1.5 n2000 | 0.998 | 0.959 (0.945, 0.970) | 998 | 0.987 | 0.972 | 0.910 | 0.957 |
| GRF | OR 0.75 n500 | 0.999 | 0.936 (0.919, 0.950) | 999 | 0.969 | 0.967 | 0.866 | 0.933 |
| GRF | OR 0.75 n2000 | 0.999 | 0.942 (0.926, 0.955) | 999 | 0.974 | 0.968 | 0.881 | 0.938 |
| GRF | OR 1.0 n500 | 0.999 | 0.934 (0.917, 0.948) | 999 | 0.966 | 0.968 | 0.879 | 0.930 |
| GRF | OR 1.0 n2000 | 1.000 | 0.939 (0.922, 0.952) | 1000 | 0.976 | 0.963 | 0.879 | 0.940 |
| GRF | OR 1.5 n500 | 1.000 | 0.933 (0.916, 0.947) | 1000 | 0.970 | 0.963 | 0.884 | 0.932 |
| GRF | OR 1.5 n2000 | 1.000 | 0.949 (0.934, 0.961) | 1000 | 0.980 | 0.969 | 0.894 | 0.949 |
| DINA | OR 0.75 n500 | 0.944 | 0.881 (0.859, 0.900) | 944 | 0.938 | 0.943 | 0.793 | 0.878 |
| DINA | OR 0.75 n2000 | 0.870 | 0.906 (0.885, 0.923) | 870 | 0.936 | 0.969 | 0.814 | 0.906 |
| DINA | OR 1.0 n500 | 0.968 | 0.893 (0.871, 0.911) | 968 | 0.943 | 0.947 | 0.803 | 0.890 |
| DINA | OR 1.0 n2000 | 0.934 | 0.898 (0.877, 0.916) | 934 | 0.937 | 0.961 | 0.812 | 0.896 |
| DINA | OR 1.5 n500 | 0.979 | 0.904 (0.884, 0.921) | 979 | 0.945 | 0.957 | 0.830 | 0.895 |
| DINA | OR 1.5 n2000 | 0.987 | 0.923 (0.905, 0.938) | 987 | 0.953 | 0.969 | 0.847 | 0.922 |

### 4.4 Identification and bound location (for reading §4.1–4.3; `or_metrics.csv` rows `mean_size_hhat`, `sensitivity_mean`, `ppv_mean`, `mean_n_family`, `share_lower_ge_tau` / `share_upper_le_tau` at τ = 1.0, and the θ† rows)

| identifier | cell | mean \|Ĥ\| | sensitivity | PPV | MR family K | field lower ≥ 1.0 on Ĥ (share) | field-s upper ≤ 1.0 on Ĥᶜ (share) | IJ two-sided coverage of θ† | field one-sided coverage of θ† |
|---|---|---|---|---|---|---|---|---|---|
| FS | OR 0.75 n500 | 101.7 | 0.240 | 0.186 | 2138 | 0.004 | 0.703 | 0.999 | 0.985 |
| FS | OR 0.75 n2000 | 138.9 | 0.089 | 0.194 | 2963 | 0.002 | 0.997 | 0.990 | 0.987 |
| FS | OR 1.0 n500 | 101.3 | 0.303 | 0.237 | 2140 | 0.007 | 0.641 | 0.999 | 0.993 |
| FS | OR 1.0 n2000 | 139.9 | 0.146 | 0.319 | 2962 | 0.010 | 0.991 | 0.960 | 0.990 |
| FS | OR 1.5 n500 | 99.9 | 0.407 | 0.326 | 2138 | 0.020 | 0.556 | 0.984 | 0.997 |
| FS | OR 1.5 n2000 | 137.5 | 0.280 | 0.588 | 2962 | 0.044 | 0.961 | 0.894 | 0.987 |
| GRF | OR 0.75 n500 | 91.1 | 0.198 | 0.160 | 1093 | 0.010 | 0.708 | 0.998 | 0.959 |
| GRF | OR 0.75 n2000 | 134.9 | 0.078 | 0.160 | 1367 | 0.009 | 0.997 | 0.990 | 0.968 |
| GRF | OR 1.0 n500 | 90.0 | 0.245 | 0.200 | 1093 | 0.017 | 0.614 | 0.999 | 0.983 |
| GRF | OR 1.0 n2000 | 140.3 | 0.131 | 0.254 | 1367 | 0.015 | 0.984 | 0.961 | 0.985 |
| GRF | OR 1.5 n500 | 90.5 | 0.332 | 0.272 | 1093 | 0.024 | 0.527 | 0.966 | 0.996 |
| GRF | OR 1.5 n2000 | 154.3 | 0.286 | 0.496 | 1367 | 0.029 | 0.956 | 0.863 | 0.992 |
| DINA | OR 0.75 n500 | 86.7 | 0.207 | 0.177 | 852 | 0.017 | 0.750 | 0.996 | 0.923 |
| DINA | OR 0.75 n2000 | 109.0 | 0.068 | 0.186 | 253 | 0.018 | 0.997 | 0.986 | 0.911 |
| DINA | OR 1.0 n500 | 86.1 | 0.243 | 0.210 | 1053 | 0.022 | 0.676 | 0.999 | 0.978 |
| DINA | OR 1.0 n2000 | 112.8 | 0.109 | 0.291 | 423 | 0.034 | 0.990 | 0.996 | 0.966 |
| DINA | OR 1.5 n500 | 87.3 | 0.333 | 0.288 | 1433 | 0.041 | 0.575 | 0.997 | 0.993 |
| DINA | OR 1.5 n2000 | 115.1 | 0.207 | 0.512 | 947 | 0.070 | 0.956 | 0.972 | 0.983 |

## 5. Reading (facts from §4; bounds read by location against the OR ladder, never as significance at OR = 1)

1. On Ĥ, the field one-sided lower bound covers β(Ĥ) at 0.964–0.978 for FS (above nominal, retained bias +0.07 to +0.40 SD, −0.10 at OR 1.5 n2000), at 0.931–0.948 for GRF, and at 0.869–0.903 for DINA (below nominal, retained bias +0.37 to +0.84 SD; conditional on the proposed family).
2. The unadjusted two-sided interval covers β(Ĥ) on 1–45 % of declared replicates; the IJ two-term interval on 0.947–0.995 (SE/SD in the extract); the oracle on 0.931–0.964.
3. On Ĥᶜ, the field-s upper bound covers β(Ĥᶜ) at 0.928–0.950 (FS), 0.931–0.944 (GRF), 0.908–0.940 (DINA), within 0.01 of the unstudentized field in every cell; the IJ two-sided interval on Ĥᶜ covers at 0.999–1.000.
4. The field-s Bonferroni pair covers (β(Ĥ), β(Ĥᶜ)) jointly at 0.955–0.964 (FS), 0.933–0.949 (GRF), 0.881–0.923 (DINA); the separate 95 % pair at 0.79–0.93.
5. Declaration: FS 0.781–0.998, rising with n and with the design point; GRF 0.999–1.000 in every cell, the protective and null points included; DINA 0.870–0.987, lowest at OR 0.75 n2000.
6. The planted region is weakly recovered in every cell (sensitivity 0.07–0.41, PPV 0.16–0.59), and the field lower bound on Ĥ sits at or above 1.0 on at most 7 % of declared replicates, including the OR 1.5 cells (FS 0.020 / 0.044); the field-s upper bound on Ĥᶜ sits at or below 1.0 on 53–75 % at n = 500 and 96–100 % at n = 2000.
7. Every estimator converged on every declared replicate in every cell (`n_nonconvergent_fits` = 0 throughout); every declared replicate carries both joint bounds.
8. Against the committed study at OR 0.75 (quoted rows, `maxeffCons` ε 0.10, IJ only, 9.632 % prevalence): FS detection 0.768 / 0.903 (n 500 / 2000) vs 0.781 / 0.910 here; the study's IJ two-sided coverage of θ† 0.997 / 0.997 vs 0.999 / 0.990 here (FS, `cov2_theta_dagger`); different rule, prevalence and constructions, so a placement, not a contrast.

## 6. Provenance: the package build that produced the bundles, and the render's

- **Bundles.** Every one of the 18 combined bundles carries `meta$pkg_version 0.3.5.9000`. The eight pop-os cells
  (`orfs` × 6, `orgrf_or150_n500`, `ordina_or150_n500`; `meta$hostname pop-os`, 63 workers; `meta$pkg_commit`
  `236ef82a`, `13d5b6c4`, `86f543b1`, `3ee0db52`, `27170d87`, `ed607c82`, `29689584`, `5ef6b748`; bundles written
  2026-09-18 21:54 to 2026-09-19 03:37 local) ran on the build of **2026-09-19 03:54:02 UTC** (catalog §2); the ten Mac
  cells (`Mac-Studio-3.local`, 13 workers; `pkg_commit` `b3d4921d` … `699254af`; written 2026-09-19 13:50 to 17:39 local)
  on the build of **2026-09-19 20:23:40 UTC** from the same package source (`REPORT_binary_stage2_2026-09-19.md` §1).
- **Every bundle predates the 2026-09-23 MR admission alignment** (`7713942e`, `96f84ad8`, `06ac5391`;
  `dev/reports/REPORT_mr_admission_alignment_2026-09-23.md`); no ACTG175 binary cell has been re-run under it.
- **The render** ran on the installed forestsearch 0.3.5.9000 built 2026-09-24 01:50:39 UTC (post-alignment); the
  summary calls the package only for `fs_sim_bias_coverage()` / `fs_plot_bias_coverage()` (display bookkeeping on the
  stored columns) and computes no MR, so the alignment cannot enter any number here.
- **Bundle md5s**, taken before the render and after every write of this task (18 of 18 identical):

```
733cb8adbcfa5c781b7117f52c737dd5  fs_effMaxSG_mr_field_or075_n2000_nb20_orfs_combined_1_1000.rds
8be2748d176b7d7f00365b688cfe53ef  fs_effMaxSG_mr_field_or075_n500_nb20_orfs_combined_1_1000.rds
321d82122cef232748289ef3c2da7569  fs_effMaxSG_mr_field_or100_n2000_nb20_orfs_combined_1_1000.rds
29adca66675561d92dc164bf021caf6b  fs_effMaxSG_mr_field_or100_n500_nb20_orfs_combined_1_1000.rds
6d2368b60153d3a13e060e3f7352e8fb  fs_effMaxSG_mr_field_or150_n2000_nb20_orfs_combined_1_1000.rds
d7c454204cc772c8e38b9a8ab8f02765  fs_effMaxSG_mr_field_or150_n500_nb20_orfs_combined_1_1000.rds
e3f2b2cdccc7ccb5bc51f74e638b5255  grf_effMaxSG_mr_field_or075_n2000_nb20_orgrf_combined_1_1000.rds
0184b338d379debbdabb23bbe2f867e4  grf_effMaxSG_mr_field_or075_n500_nb20_orgrf_combined_1_1000.rds
c24f87e6c9a38c9bf63284304d2dafa6  grf_effMaxSG_mr_field_or100_n2000_nb20_orgrf_combined_1_1000.rds
ceebd206b38cbc234cfa560165e83eb9  grf_effMaxSG_mr_field_or100_n500_nb20_orgrf_combined_1_1000.rds
7c5369ed72614a4d26671ea60cb53dd6  grf_effMaxSG_mr_field_or150_n2000_nb20_orgrf_combined_1_1000.rds
3c29d466b9389deda218f31b89ae1ab6  grf_effMaxSG_mr_field_or150_n500_nb20_orgrf_combined_1_1000.rds
2b35d27817fffdec020a63402543587e  dina_effMaxSG_mr_field_or075_n2000_nb20_ordina_combined_1_1000.rds
9eddc985b79a68d2f6edd769dfb71022  dina_effMaxSG_mr_field_or075_n500_nb20_ordina_combined_1_1000.rds
d524fd56f722267387214d77b229abae  dina_effMaxSG_mr_field_or100_n2000_nb20_ordina_combined_1_1000.rds
045d823f1cba1cdcf6c9630fa6f376c0  dina_effMaxSG_mr_field_or100_n500_nb20_ordina_combined_1_1000.rds
d0c2b24abf8f8714c01b44fbf52b1036  dina_effMaxSG_mr_field_or150_n2000_nb20_ordina_combined_1_1000.rds
4738d678192ce32b6ad86a302a607964  dina_effMaxSG_mr_field_or150_n500_nb20_ordina_combined_1_1000.rds
```

## 7. Files, catalog, post-conditions

- Added: `or_metrics.csv`, `COLUMNS_or.md`, `summary_actg175_or.html`, this report. Modified: `summary_actg175_or.qmd`
  (§3), `current_status.md` (in place: a Stage 3 line in its preamble and the four Stage 3 documents in §3.1's "campaign
  documents" column for each campaign). Nothing else. The two render-side PNGs were removed (§3).
- `status_curated.md` was **not** edited (the task's post-conditions do not list it), so the directory's regenerator
  (`scripts_or/current_status_regen.R`, which includes `status_curated.md` verbatim) would not carry the Stage 3 preamble
  line on its next run; §3.1's document column is generated and would pick the new files up by its `docre` pattern.
  The catalog's pin (`d356bd8a`) is unchanged; `check_current_status.sh` was not run (its pin rule is for the regenerator's
  closeout commit).
- Copies in `~/Downloads/actg175_or_stage3_2026-09-28/`: this report, `or_metrics.csv`, `COLUMNS_or.md`.
- Scope, verbatim from the campaign task: "These are operating characteristics on one binary design at three effect
  sizes, with GRF's and DINA's figures conditional on the proposed family. They do not verify condition (A3), no
  construction is promoted on this design's performance, and they supersede nothing in the committed study, which ran a
  different rule and different constructions." No recommendations.
