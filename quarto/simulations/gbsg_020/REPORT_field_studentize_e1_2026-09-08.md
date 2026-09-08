# REPORT — Studentized complement field E1: the six-cell campaign `e1stud` (field-s, variant R1) — Gate 2 records, the standard tables and the five pre-registered findings (report-and-wait)

**Task:** `dev/tasks/TASK_field_studentize_e1_2026-09-08.md` (9a6e21cc), Stages 2–3. Stage 0/1 record: `REPORT_field_studentize_e1_stage1_2026-09-08.md` (3880022e; G2a/G2b/G2c PASS; Gate 1 compute go). Proposal §5–§6 (R1). Predecessors: E0 record 9d71f736, Part A `REPORT_complement_variance_2026-09-07.md`, `REPORT_banddial_2026-09-07.md`.
**Date:** 2026-09-08. Executor: Claude Code (Linux), unattended. Winner-only and winner-floor excluded from every table and line. Bounds read by location; Wilson intervals; the marginal-vs-error SD caveat (Part A, A0) applies to every SD-unit column — both SDs are shown. **No adoption recommendation; no change to any documented rule; the record reports, the chat and Larry decide.** The ε 0.30 cells are a mapped stress comparator only (ε ≤ 0.25 cap).

---

## GATE 2: PASS on all six cells — 2,000 rows each, no duplicates, meta as set (`field_decompose = TRUE`, `field_scale_complement = selected`, seed 8316951); **131 / 131 pre-existing non-timing columns `identical()` to the committed comparator on every cell** (the harm block, the unscaled complement field, the joint pair and the identification are untouched: the add-beside proof), `truth` identical; every `_s` and `joint_s` column finite on all 1,999 detected rows per cell; interval invariants hold on the `_s`, unscaled and harm blocks; `joint_s` γ ∈ [0.025, 0.05] with achieved probability ≥ 0.95 off the α/2 fallback. Campaign wall **2 h 53 m** at 100 workers (31–32 min per band cell, 23–24 min per size-rule cell; ceiling 3.5 h, timeout 5 h); **six of six cells completed, none deferred, none dropped.**

## Run

Committed template (3880022e) driven by env only: `FS_S7_Z1Q=0.60 FS_S7_N=500 FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_FIELD_DECOMP=TRUE FS_S7_FIELD_SCALEC=selected FS_S7_CAMPAIGN=e1stud FS_S7_WORKERS=100` (J = 10 default; `return_reselection = TRUE`), per cell `FS_S7_FOCUS` / `FS_S7_HR` / `FS_S7_NBHD`, two seed-disjoint batches (`FS_S7_START` 1 and 1001, 1,000 each) then `FS_S7_MODE=combine`; seeds `8316951 + sim_id`. Driver `e1_campaign.sh` (session scratchpad; re-projects each cell from the realized wall of its kind and defers any cell crossing the ceiling — none did). Stems: `fs_effMaxSG_…_z1q60_nb20_e1stud` (HR 1.50 / 1.75; comparator `…_nb20_p30sgnb20`), `fs_effMaxSG_…_z1q60_nb30_e1stud` (comparator `…_nb30_banddial`), `fs_maxSG_…_z1q60_e1stud`, `fs_minSG_…_z1q60_e1stud` (comparators `…_banddial`); per cell two batch renders and one combine render committed beside the bundles. Cross-cell document `summary_e1stud.qmd` (transplanted from the committed `summary_banddial.qmd`: bundle list, field-s rows and the `joint_s` pairs added by column substitution — no fresh authorship of the table machinery) rendered to `summary_e1stud.html`; the tables below are the verbatim output of `e1stud_findings.R` (this directory), which implements the same definitions.

| cell | wall | Gate 2 | pre-existing identical | `_s` finite | mean ρᶜ | `se_field_s / (ρᶜ · se_field)` mean (min–max) | joint_s at α/2 fallback (joint) |
|---|---|---|---|---|---|---|---|
| effMaxSG ε 0.20, HR 1.50 | 30 m 58 s | PASS | 131 / 131 | 1999 / 1999 | 1.054 | 0.998 (0.983–1.013) | 0.864 (0.917) |
| effMaxSG ε 0.20, HR 1.75 | 31 m 41 s | PASS | 131 / 131 | 1999 / 1999 | 1.058 | 0.999 (0.983–1.011) | 0.859 (0.910) |
| effMaxSG ε 0.30, HR 1.50 | 31 m 31 s | PASS | 131 / 131 | 1999 / 1999 | 1.087 | 0.996 (0.963–1.015) | 0.862 (0.938) |
| effMaxSG ε 0.30, HR 1.75 | 31 m 58 s | PASS | 131 / 131 | 1999 / 1999 | 1.083 | 0.998 (0.970–1.014) | 0.858 (0.937) |
| maxSG, HR 1.75 | 23 m 29 s | PASS | 131 / 131 | 1999 / 1999 | 0.997 | 0.969 (0.813–1.031) | 0.873 (0.561) |
| minSG, HR 1.75 | 23 m 45 s | PASS | 131 / 131 | 1999 / 1999 | 1.006 | 1.000 (0.995–1.004) | 0.847 (0.860) |

## The five pre-registered findings, stated as findings

**1. Coverage (band cells; target ≥ 0.92 with Wilson support; committed field 0.897–0.912).** Field-s one-sided upper coverage of β(Ĥᶜ): **0.912 [0.899, 0.924], 0.919 [0.907, 0.931], 0.921 [0.909, 0.932], 0.919 [0.906, 0.930]** against the field's 0.897, 0.912, 0.904, 0.903 — gains of +1.6, +0.7, +1.8, +1.6 points, each from replicates that flip to cover (1.0–1.8%) against 0.1–0.3% that flip to miss. The point estimate reaches 0.92 in one band cell (ε 0.30 HR 1.50, 0.921); **no band cell reaches 0.92 with Wilson support** (lower limits 0.899–0.909; the Wilson intervals of field and field-s overlap in every cell). IJ two-term sits at 0.995–0.996. The harm block and the unscaled complement are unchanged (structural, by pairing). On the error scale the studentized SD restores the calibration the field lacked: mean SE / error SD 0.96–0.97 (field 0.89–0.92; naive 0.99–1.00), λ-SDᶜ/naive SE 0.997–1.002 (field 0.917–0.947); the Gaussian reference coverage for field-s is 0.900–0.910 against the observed 0.912–0.921.

**2. Shape (upper coverage by |Ĥ|/|H| tertile).** The field's falling pattern is **flattened** on the |Ĥ|/|H| axis in all four band cells: 0.918 / 0.892 / 0.881 → 0.916 / 0.904 / 0.917 (ε 0.20 HR 1.50); 0.922 / 0.905 / 0.910 → 0.915 / 0.908 / 0.935 (ε 0.20 HR 1.75); 0.909 / 0.916 / 0.887 → 0.910 / 0.931 / 0.923 (ε 0.30 HR 1.50); 0.913 / 0.919 / 0.878 → 0.913 / 0.929 / 0.914 (ε 0.30 HR 1.75) — the top tertile gains 2.5–3.6 points, the bottom is unchanged (ρᶜ ≈ 1.00–1.02 there), and the studentized per-replicate ratio `fld_Hc_se_s`/`nv_Hc_se` is 0.993–1.003 in every tertile where the unscaled one fell from 1.00 to 0.86–0.89. **On the p̂ axis the falling pattern remains**: field-s 0.936 / 0.917 / 0.884, 0.954 / 0.928 / 0.877, 0.945 / 0.932 / 0.887, 0.936 / 0.937 / 0.884 — the top-p̂ tertile (p̂ > 0.07–0.13) sits at 0.877–0.887 under both constructions (field 0.872–0.878), where ρᶜ is already 1.01–1.04 and the unscaled ratio 0.97–0.99: the shortfall in the concentrated-pick tertile is not a scale effect, consistent with the proposal's stated out-of-scope residual (§2 (iii)); the low-p̂ tertile gains 1.1–2.5 points.

**3. No-regression ends.** `maxSG`: field 0.925 [0.913, 0.936] → field-s 0.926 [0.914, 0.937] (mean ρᶜ 0.997); `minSG`: 0.929 [0.917, 0.940] → 0.929 [0.917, 0.939] (mean ρᶜ 1.006) — within Wilson overlap, differences +0.002 / −0.001. One structure inside `maxSG` worth the record: with the pick spread over candidates of very different complement sizes (winner-scale CV 0.26 per replicate, ρᶜ 0.83 / 1.03 / 1.13 across the |Ĥ|/|H| tertiles), field-s redistributes coverage across tertiles — 0.916 / 0.925 / 0.934 → 0.888 / 0.937 / 0.956 (1.3% flip to cover, 1.3% flip to miss) — with the cell total unchanged; the bottom tertile (|Ĥ|/|H| < 2.5, the replicates where `maxSG` picked a smaller candidate than the field typically re-selects, so ρᶜ < 1 and the studentized SD is *smaller*) loses 2.8 points. `minSG` is inert (ρᶜ 1.004–1.007 in every tertile).

**4. Bound locations (HR scale, complement upper bound).** Band cells: field-s upper mean 0.997 / 1.010 / 0.977 / 0.981 against field 0.987 / 0.999 / 0.961 / 0.965 and IJ 1.188 / 1.207 / 1.159 / 1.171; share < 0.85: field-s 0.164 / 0.156 / 0.214 / 0.210, field 0.192 / 0.178 / 0.252 / 0.249, IJ 0.019–0.028; share < 0.80: 0.092 / 0.084 / 0.129 / 0.128 vs 0.109 / 0.104 / 0.160 / 0.151 vs 0.006–0.011; mean margin over β̃ᶜ 0.246–0.266 log units (field 0.234–0.248, IJ 0.42–0.44). Joint pair: `joint_s` Bonferroni / calibrated 0.942 / 0.939 / 0.949 / 0.951–0.952 against `joint` 0.932 / 0.933 / 0.935 / 0.935 in the band cells (the calibrated γ sits at the α/2 fallback on 86% of replicates for `joint_s`, 91–94% for `joint`, so calibrated ≈ Bonferroni for both); ends 0.939 vs 0.941–0.942 (`maxSG`), 0.930 vs 0.929–0.930 (`minSG`). Separate 95% pairs: 0.888–0.904 (field-s) vs 0.872–0.886 (field) in the band cells. Harm-side locations identical between field and field-s by construction.

**5. R1 vs the global rescale (informational).** `se_field_s / (ρᶜ · se_field)`: mean 0.996–0.999 with SD 0.0035–0.0057 in the band cells (q01–q99 0.978–1.008; winner-scale CV 0.05–0.07 mean, 0.07–0.09 q90), 0.9999 (SD 0.0009) at `minSG` (CV 0.015), **0.969 (SD 0.032, q01 0.88, min 0.81) at `maxSG`** (CV 0.26 mean, 0.375 q90) — the per-draw and global forms coincide to within 1% wherever the winners' scales are homogeneous and separate only at `maxSG`, where R1's candidate-wise weights are on average smaller than the global ρᶜ (the correction estimate moves with the studentization there; the E2-side question, not this record's). Per replicate the studentized SD tracks the naive SE (corr 0.94–0.96 in the band cells vs 0.34–0.42 for the unscaled; SD across replicates 0.0105–0.0149 vs 0.0054–0.0070, matching the naive SE's 0.0102–0.0148).

## What this record does and does not say

Reported: the add-beside construction is inert on every committed column in all six cells; field-s raises the complement upper coverage in the four band cells by 0.7–1.8 points to 0.912–0.921 without reaching 0.92 with Wilson support; it flattens the |Ĥ|/|H| shape and leaves the p̂ shape (the concentrated-pick shortfall) in place; the ends are unchanged in total, with a within-cell redistribution at `maxSG`; the studentized SD is calibrated to the naive SE on both the average and the per-replicate scale; R1 = global to 1% except at `maxSG`. Not said: whether field-s is adopted, whether any documented rule changes, or what E2 should be. Report and wait.

---

# Tables (verbatim output of `e1stud_findings.R`)

## Gate 2 per cell (pairing identity to the committed comparator)

| cell | rows | detected | pre-existing non-timing columns identical | truth identical | _s finite on detected | mean rho^c | joint_s at the alpha/2 fallback (joint) |
|---|---|---|---|---|---|---|---|
| effMaxSG eps 0.20, HR 1.50 | 2000 | 1999 | 131 / 131 | TRUE | TRUE | 1.054 | 0.864 (0.917) |
| effMaxSG eps 0.20, HR 1.75 | 2000 | 1999 | 131 / 131 | TRUE | TRUE | 1.058 | 0.859 (0.910) |
| effMaxSG eps 0.30 (stress), HR 1.50 | 2000 | 1999 | 131 / 131 | TRUE | TRUE | 1.087 | 0.862 (0.938) |
| effMaxSG eps 0.30 (stress), HR 1.75 | 2000 | 1999 | 131 / 131 | TRUE | TRUE | 1.083 | 0.858 (0.937) |
| maxSG, HR 1.75 | 2000 | 1999 | 131 / 131 | TRUE | TRUE | 0.997 | 0.873 (0.561) |
| minSG, HR 1.75 | 2000 | 1999 | 131 / 131 | TRUE | TRUE | 1.006 | 0.847 (0.860) |

## Constructions per cell (identical replicates): naive / field / field-s / IJ two-term

| cell | block | construction | n | bias (log) | marginal SD | error SD | mean SE | b | r = SE/marg SD | SE/error SD | one-sided cov [Wilson] | two-sided cov | Gaussian ref |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| effMaxSG eps 0.20, HR 1.50 | Hhat (lower) | naive | 1999 | 0.430 | 0.238 | 0.281 | 0.250 | 1.809 | 1.049 | 0.887 | 0.484 [0.462, 0.506] | 0.619 | 0.467 |
| effMaxSG eps 0.20, HR 1.50 | Hhat (lower) | field | 1999 | -0.056 | 0.299 | 0.314 | 0.324 | -0.185 | 1.081 | 1.029 | 0.974 [0.967, 0.981] | 0.925 | 0.975 |
| effMaxSG eps 0.20, HR 1.50 | Hhat (lower) | IJ two-term | 1999 | -0.003 | 0.279 | 0.303 | 0.360 | -0.011 | 1.289 | 1.186 | 0.985 [0.979, 0.990] | 0.981 | 0.983 |
| effMaxSG eps 0.20, HR 1.50 | Hhat^c (upper) | naive | 1999 | -0.123 | 0.156 | 0.143 | 0.143 | -0.789 | 0.912 | 0.997 | 0.781 [0.762, 0.798] | 0.849 | 0.762 |
| effMaxSG eps 0.20, HR 1.50 | Hhat^c (upper) | field | 1999 | -0.031 | 0.159 | 0.147 | 0.135 | -0.193 | 0.852 | 0.923 | 0.897 [0.883, 0.910] | 0.921 | 0.886 |
| effMaxSG eps 0.20, HR 1.50 | Hhat^c (upper) | field-s | 1999 | -0.031 | 0.159 | 0.147 | 0.142 | -0.194 | 0.896 | 0.970 | 0.912 [0.899, 0.924] | 0.938 | 0.900 |
| effMaxSG eps 0.20, HR 1.50 | Hhat^c (upper) | IJ two-term | 1999 | -0.044 | 0.157 | 0.145 | 0.257 | -0.282 | 1.636 | 1.776 | 0.995 [0.991, 0.997] | 1.000 | 0.992 |
| effMaxSG eps 0.20, HR 1.75 | Hhat (lower) | naive | 1999 | 0.365 | 0.246 | 0.294 | 0.243 | 1.481 | 0.985 | 0.825 | 0.573 [0.551, 0.595] | 0.694 | 0.555 |
| effMaxSG eps 0.20, HR 1.75 | Hhat (lower) | field | 1999 | -0.078 | 0.312 | 0.329 | 0.323 | -0.251 | 1.035 | 0.983 | 0.970 [0.962, 0.977] | 0.911 | 0.975 |
| effMaxSG eps 0.20, HR 1.75 | Hhat (lower) | IJ two-term | 1999 | -0.037 | 0.291 | 0.317 | 0.356 | -0.129 | 1.227 | 1.125 | 0.985 [0.979, 0.989] | 0.972 | 0.984 |
| effMaxSG eps 0.20, HR 1.75 | Hhat^c (upper) | naive | 1999 | -0.102 | 0.162 | 0.144 | 0.143 | -0.631 | 0.884 | 0.991 | 0.826 [0.809, 0.842] | 0.880 | 0.795 |
| effMaxSG eps 0.20, HR 1.75 | Hhat^c (upper) | field | 1999 | -0.019 | 0.165 | 0.149 | 0.136 | -0.116 | 0.822 | 0.911 | 0.912 [0.899, 0.924] | 0.918 | 0.892 |
| effMaxSG eps 0.20, HR 1.75 | Hhat^c (upper) | field-s | 1999 | -0.019 | 0.165 | 0.149 | 0.143 | -0.117 | 0.868 | 0.961 | 0.919 [0.907, 0.931] | 0.933 | 0.905 |
| effMaxSG eps 0.20, HR 1.75 | Hhat^c (upper) | IJ two-term | 1999 | -0.030 | 0.163 | 0.147 | 0.258 | -0.183 | 1.586 | 1.761 | 0.995 [0.991, 0.998] | 0.999 | 0.992 |
| effMaxSG eps 0.30 (stress), HR 1.50 | Hhat (lower) | naive | 1999 | 0.299 | 0.224 | 0.244 | 0.213 | 1.332 | 0.949 | 0.873 | 0.633 [0.612, 0.654] | 0.735 | 0.590 |
| effMaxSG eps 0.30 (stress), HR 1.50 | Hhat (lower) | field | 1999 | -0.093 | 0.270 | 0.267 | 0.298 | -0.344 | 1.103 | 1.116 | 0.980 [0.973, 0.986] | 0.937 | 0.985 |
| effMaxSG eps 0.30 (stress), HR 1.50 | Hhat (lower) | IJ two-term | 1999 | -0.068 | 0.254 | 0.258 | 0.334 | -0.266 | 1.311 | 1.294 | 0.995 [0.991, 0.997] | 0.987 | 0.992 |
| effMaxSG eps 0.30 (stress), HR 1.50 | Hhat^c (upper) | naive | 1999 | -0.137 | 0.162 | 0.153 | 0.153 | -0.844 | 0.944 | 1.001 | 0.774 [0.755, 0.792] | 0.854 | 0.761 |
| effMaxSG eps 0.30 (stress), HR 1.50 | Hhat^c (upper) | field | 1999 | -0.031 | 0.167 | 0.156 | 0.141 | -0.187 | 0.845 | 0.900 | 0.904 [0.890, 0.916] | 0.916 | 0.886 |
| effMaxSG eps 0.30 (stress), HR 1.50 | Hhat^c (upper) | field-s | 1999 | -0.032 | 0.167 | 0.157 | 0.152 | -0.190 | 0.915 | 0.974 | 0.921 [0.909, 0.932] | 0.936 | 0.906 |
| effMaxSG eps 0.30 (stress), HR 1.50 | Hhat^c (upper) | IJ two-term | 1999 | -0.047 | 0.164 | 0.154 | 0.266 | -0.289 | 1.625 | 1.726 | 0.996 [0.992, 0.998] | 1.000 | 0.991 |
| effMaxSG eps 0.30 (stress), HR 1.75 | Hhat (lower) | naive | 1999 | 0.254 | 0.239 | 0.254 | 0.210 | 1.065 | 0.879 | 0.825 | 0.680 [0.660, 0.700] | 0.777 | 0.648 |
| effMaxSG eps 0.30 (stress), HR 1.75 | Hhat (lower) | field | 1999 | -0.097 | 0.286 | 0.277 | 0.297 | -0.339 | 1.038 | 1.074 | 0.976 [0.968, 0.982] | 0.925 | 0.980 |
| effMaxSG eps 0.30 (stress), HR 1.75 | Hhat (lower) | IJ two-term | 1999 | -0.082 | 0.271 | 0.268 | 0.333 | -0.303 | 1.231 | 1.243 | 0.994 [0.990, 0.997] | 0.978 | 0.990 |
| effMaxSG eps 0.30 (stress), HR 1.75 | Hhat^c (upper) | naive | 1999 | -0.113 | 0.164 | 0.153 | 0.152 | -0.691 | 0.928 | 0.996 | 0.807 [0.790, 0.824] | 0.875 | 0.798 |
| effMaxSG eps 0.30 (stress), HR 1.75 | Hhat^c (upper) | field | 1999 | -0.021 | 0.171 | 0.158 | 0.141 | -0.123 | 0.828 | 0.894 | 0.903 [0.890, 0.916] | 0.914 | 0.892 |
| effMaxSG eps 0.30 (stress), HR 1.75 | Hhat^c (upper) | field-s | 1999 | -0.022 | 0.171 | 0.158 | 0.153 | -0.127 | 0.894 | 0.965 | 0.919 [0.906, 0.930] | 0.941 | 0.910 |
| effMaxSG eps 0.30 (stress), HR 1.75 | Hhat^c (upper) | IJ two-term | 1999 | -0.033 | 0.167 | 0.155 | 0.268 | -0.198 | 1.600 | 1.723 | 0.996 [0.993, 0.998] | 0.999 | 0.993 |
| maxSG, HR 1.75 | Hhat (lower) | naive | 1999 | 0.046 | 0.063 | 0.095 | 0.133 | 0.739 | 2.118 | 1.400 | 0.947 [0.936, 0.956] | 0.974 | 0.997 |
| maxSG, HR 1.75 | Hhat (lower) | field | 1999 | -0.005 | 0.087 | 0.131 | 0.126 | -0.056 | 1.455 | 0.961 | 0.947 [0.936, 0.956] | 0.956 | 0.993 |
| maxSG, HR 1.75 | Hhat (lower) | IJ two-term | 1999 | 0.002 | 0.080 | 0.123 | 0.221 | 0.025 | 2.773 | 1.799 | 0.999 [0.997, 1.000] | 0.996 | 1.000 |
| maxSG, HR 1.75 | Hhat^c (upper) | naive | 1999 | -0.191 | 0.393 | 0.406 | 0.355 | -0.487 | 0.903 | 0.876 | 0.776 [0.757, 0.794] | 0.844 | 0.841 |
| maxSG, HR 1.75 | Hhat^c (upper) | field | 1999 | -0.041 | 0.392 | 0.395 | 0.361 | -0.104 | 0.921 | 0.915 | 0.925 [0.913, 0.936] | 0.940 | 0.921 |
| maxSG, HR 1.75 | Hhat^c (upper) | field-s | 1999 | -0.039 | 0.393 | 0.395 | 0.358 | -0.100 | 0.912 | 0.906 | 0.926 [0.914, 0.937] | 0.940 | 0.919 |
| maxSG, HR 1.75 | Hhat^c (upper) | IJ two-term | 1999 | -0.064 | 0.389 | 0.394 | 0.621 | -0.164 | 1.598 | 1.576 | 0.995 [0.991, 0.997] | 0.999 | 0.993 |
| minSG, HR 1.75 | Hhat (lower) | naive | 1999 | 0.326 | 0.217 | 0.272 | 0.320 | 1.504 | 1.476 | 1.177 | 0.804 [0.786, 0.821] | 0.881 | 0.822 |
| minSG, HR 1.75 | Hhat (lower) | field | 1999 | 0.022 | 0.284 | 0.326 | 0.315 | 0.078 | 1.110 | 0.967 | 0.936 [0.924, 0.946] | 0.922 | 0.960 |
| minSG, HR 1.75 | Hhat (lower) | IJ two-term | 1999 | 0.074 | 0.258 | 0.305 | 0.391 | 0.288 | 1.514 | 1.280 | 0.971 [0.963, 0.977] | 0.985 | 0.986 |
| minSG, HR 1.75 | Hhat^c (upper) | naive | 1999 | -0.049 | 0.139 | 0.136 | 0.130 | -0.350 | 0.932 | 0.958 | 0.879 [0.864, 0.893] | 0.924 | 0.882 |
| minSG, HR 1.75 | Hhat^c (upper) | field | 1999 | -0.013 | 0.137 | 0.134 | 0.129 | -0.094 | 0.940 | 0.967 | 0.929 [0.917, 0.940] | 0.940 | 0.927 |
| minSG, HR 1.75 | Hhat^c (upper) | field-s | 1999 | -0.013 | 0.137 | 0.134 | 0.130 | -0.094 | 0.945 | 0.972 | 0.929 [0.917, 0.939] | 0.943 | 0.928 |
| minSG, HR 1.75 | Hhat^c (upper) | IJ two-term | 1999 | -0.019 | 0.137 | 0.134 | 0.251 | -0.139 | 1.830 | 1.881 | 0.997 [0.994, 0.999] | 1.000 | 0.998 |

**Across cells** (one-sided coverage: min / mean / max; mean two-sided; mean b; mean r; mean SE/error SD):

| block | construction | cells | cov1 min | cov1 mean | cov1 max | cov2 mean | b mean | r mean | SE/error SD mean |
|---|---|---|---|---|---|---|---|---|---|
| Hhat (lower) | naive | 6 | 0.484 | 0.687 | 0.947 | 0.780 | 1.322 | 1.243 | 0.998 |
| Hhat (lower) | field | 6 | 0.936 | 0.964 | 0.980 | 0.929 | -0.183 | 1.137 | 1.022 |
| Hhat (lower) | IJ two-term | 6 | 0.971 | 0.988 | 0.999 | 0.983 | -0.066 | 1.557 | 1.321 |
| Hhat^c (upper) | naive | 6 | 0.774 | 0.807 | 0.879 | 0.871 | -0.632 | 0.917 | 0.970 |
| Hhat^c (upper) | field | 6 | 0.897 | 0.912 | 0.929 | 0.925 | -0.136 | 0.868 | 0.918 |
| Hhat^c (upper) | field-s | 6 | 0.912 | 0.921 | 0.929 | 0.939 | -0.137 | 0.905 | 0.958 |
| Hhat^c (upper) | IJ two-term | 6 | 0.995 | 0.996 | 0.997 | 1.000 | -0.209 | 1.646 | 1.741 |

## Finding 1 -- complement one-sided upper coverage: field / field-s / IJ two-term (Wilson)

| cell | n | field | field-s | IJ two-term | field-s - field | flips to cover / to miss | mean rho^c | lambda-SD^c/nSE unscaled -> studentized | field-s Wilson lower >= 0.92 | point >= 0.92 |
|---|---|---|---|---|---|---|---|---|---|---|
| effMaxSG eps 0.20, HR 1.50 | 1999 | 0.897 [0.883, 0.910] | 0.912 [0.899, 0.924] | 0.995 [0.991, 0.997] | +0.016 | 0.018 / 0.002 | 1.054 | 0.947 -> 0.998 | FALSE | FALSE |
| effMaxSG eps 0.20, HR 1.75 | 1999 | 0.912 [0.899, 0.924] | 0.919 [0.907, 0.931] | 0.995 [0.991, 0.998] | +0.007 | 0.010 / 0.003 | 1.058 | 0.946 -> 1.001 | FALSE | FALSE |
| effMaxSG eps 0.30 (stress), HR 1.50 | 1999 | 0.904 [0.890, 0.916] | 0.921 [0.909, 0.932] | 0.996 [0.992, 0.998] | +0.018 | 0.018 / 0.001 | 1.087 | 0.917 -> 0.997 | FALSE | TRUE |
| effMaxSG eps 0.30 (stress), HR 1.75 | 1999 | 0.903 [0.890, 0.916] | 0.919 [0.906, 0.930] | 0.996 [0.993, 0.998] | +0.016 | 0.017 / 0.002 | 1.083 | 0.925 -> 1.002 | FALSE | FALSE |
| maxSG, HR 1.75 | 1999 | 0.925 [0.913, 0.936] | 0.926 [0.914, 0.937] | 0.995 [0.991, 0.997] | +0.002 | 0.015 / 0.013 | 0.997 | 0.988 -> 1.014 | FALSE | TRUE |
| minSG, HR 1.75 | 1999 | 0.929 [0.917, 0.940] | 0.929 [0.917, 0.939] | 0.997 [0.994, 0.999] | -0.001 | 0.002 / 0.002 | 1.006 | 0.994 -> 1.000 | FALSE | TRUE |

## Finding 2 -- shape: complement upper coverage by |Hhat|/|H| tertile and by p-hat tertile (field -> field-s)

| cell | stratification | T1 range / T2 / T3 | n | mean rho^c T1 / T2 / T3 | fld/nv T1 / T2 / T3 | flds/nv T1 / T2 / T3 | field T1 / T2 / T3 | field-s T1 / T2 / T3 | field-s Wilson T1 / T2 / T3 |
|---|---|---|---|---|---|---|---|---|---|
| effMaxSG eps 0.20, HR 1.50 | |Hhat|/|H| | [0.37, 0.69] / [0.69, 0.95] / [0.95, 2.09] | 667 / 666 / 666 | 1.001 / 1.041 / 1.120 | 1.002 / 0.964 / 0.892 | 0.999 / 1.000 / 0.994 | 0.918 / 0.892 / 0.881 | 0.916 / 0.904 / 0.917 | [0.893, 0.935] / [0.879, 0.924] / [0.894, 0.936] |
| effMaxSG eps 0.20, HR 1.50 | p-hat | [0.00, 0.04] / [0.04, 0.11] / [0.11, 0.82] | 669 / 664 / 666 | 1.103 / 1.047 / 1.012 | 0.908 / 0.956 / 0.994 | 0.994 / 0.996 / 1.003 | 0.912 / 0.907 / 0.872 | 0.936 / 0.917 / 0.884 | [0.915, 0.952] / [0.894, 0.936] / [0.858, 0.906] |
| effMaxSG eps 0.20, HR 1.75 | |Hhat|/|H| | [0.35, 0.71] / [0.71, 0.97] / [0.97, 1.93] | 667 / 666 / 666 | 1.002 / 1.047 / 1.125 | 1.002 / 0.960 / 0.890 | 1.002 / 1.003 / 0.997 | 0.922 / 0.905 / 0.910 | 0.915 / 0.908 / 0.935 | [0.891, 0.933] / [0.884, 0.928] / [0.914, 0.952] |
| effMaxSG eps 0.20, HR 1.75 | p-hat | [0.00, 0.05] / [0.05, 0.13] / [0.13, 0.93] | 667 / 666 / 666 | 1.107 / 1.049 / 1.017 | 0.906 / 0.955 / 0.992 | 0.997 / 0.998 / 1.007 | 0.943 / 0.916 / 0.878 | 0.954 / 0.928 / 0.877 | [0.935, 0.967] / [0.906, 0.945] / [0.850, 0.900] |
| effMaxSG eps 0.30 (stress), HR 1.50 | |Hhat|/|H| | [0.37, 0.97] / [0.97, 1.24] / [1.24, 2.41] | 667 / 666 / 666 | 1.017 / 1.082 / 1.162 | 0.987 / 0.928 / 0.864 | 0.999 / 1.000 / 0.993 | 0.909 / 0.916 / 0.887 | 0.910 / 0.931 / 0.923 | [0.886, 0.929] / [0.909, 0.948] / [0.901, 0.941] |
| effMaxSG eps 0.30 (stress), HR 1.50 | p-hat | [0.00, 0.03] / [0.03, 0.07] / [0.07, 0.69] | 671 / 662 / 666 | 1.151 / 1.073 / 1.037 | 0.874 / 0.935 / 0.971 | 0.993 / 0.996 / 1.003 | 0.920 / 0.917 / 0.875 | 0.945 / 0.932 / 0.887 | [0.925, 0.960] / [0.910, 0.949] / [0.861, 0.909] |
| effMaxSG eps 0.30 (stress), HR 1.75 | |Hhat|/|H| | [0.39, 0.98] / [0.98, 1.22] / [1.22, 2.35] | 667 / 666 / 666 | 1.020 / 1.080 / 1.150 | 0.987 / 0.932 / 0.876 | 1.003 / 1.003 / 0.999 | 0.913 / 0.919 / 0.878 | 0.913 / 0.929 / 0.914 | [0.889, 0.932] / [0.907, 0.947] / [0.891, 0.933] |
| effMaxSG eps 0.30 (stress), HR 1.75 | p-hat | [0.00, 0.03] / [0.04, 0.08] / [0.08, 0.87] | 667 / 666 / 666 | 1.143 / 1.069 / 1.037 | 0.880 / 0.940 / 0.975 | 0.996 / 1.000 / 1.009 | 0.918 / 0.920 / 0.872 | 0.936 / 0.937 / 0.884 | [0.914, 0.952] / [0.916, 0.953] / [0.858, 0.906] |
| maxSG, HR 1.75 | |Hhat|/|H| | [0.45, 2.52] / [2.52, 2.89] / [2.89, 3.64] | 667 / 678 / 654 | 0.828 / 1.033 / 1.131 | 1.262 / 1.029 / 0.924 | 0.971 / 1.010 / 1.021 | 0.916 / 0.925 / 0.934 | 0.888 / 0.937 / 0.956 | [0.861, 0.909] / [0.916, 0.953] / [0.937, 0.969] |
| maxSG, HR 1.75 | p-hat | [0.01, 0.12] / [0.12, 0.64] / [0.64, 1.00] | 667 / 666 / 666 | 0.824 / 1.086 / 1.081 | 1.271 / 0.993 / 0.953 | 0.973 / 1.003 / 1.026 | 0.948 / 0.874 / 0.953 | 0.922 / 0.893 / 0.964 | [0.899, 0.940] / [0.868, 0.915] / [0.947, 0.976] |
| minSG, HR 1.75 | |Hhat|/|H| | [0.33, 0.39] / [0.39, 0.42] / [0.42, 1.04] | 679 / 654 / 666 | 1.004 / 1.007 / 1.007 | 0.995 / 0.994 / 0.995 | 0.999 / 1.000 / 1.001 | 0.946 / 0.927 / 0.916 | 0.944 / 0.931 / 0.911 | [0.924, 0.959] / [0.909, 0.948] / [0.887, 0.931] |
| minSG, HR 1.75 | p-hat | [0.00, 0.01] / [0.01, 0.08] / [0.08, 0.87] | 669 / 664 / 666 | 1.005 / 1.006 / 1.007 | 0.998 / 0.997 / 0.989 | 1.002 / 1.002 / 0.996 | 0.955 / 0.935 / 0.898 | 0.957 / 0.932 / 0.898 | [0.938, 0.970] / [0.911, 0.949] / [0.873, 0.919] |

## Finding 4 -- bound locations (HR scale) and the joint pair

| cell | construction | mean beta(Hhat) | H lower mean | share >= 0.85 | share >= 0.95 | mean beta(Hhat^c) | Hc upper mean | share < 0.85 | share < 0.80 | margin_Hc (log) |
|---|---|---|---|---|---|---|---|---|---|---|
| effMaxSG eps 0.20, HR 1.50 | naive | 1.256 | 1.289 | 0.985 | 0.927 | 0.808 | 0.909 | 0.349 | 0.216 | 0.156 |
| effMaxSG eps 0.20, HR 1.50 | field | 1.256 | 0.693 | 0.196 | 0.115 | 0.808 | 0.987 | 0.192 | 0.109 | 0.237 |
| effMaxSG eps 0.20, HR 1.50 | field-s | 1.256 | 0.693 | 0.196 | 0.115 | 0.808 | 0.997 | 0.164 | 0.092 | 0.247 |
| effMaxSG eps 0.20, HR 1.50 | IJ two-term | 1.256 | 0.705 | 0.181 | 0.095 | 0.808 | 1.188 | 0.020 | 0.007 | 0.422 |
| effMaxSG eps 0.20, HR 1.75 | naive | 1.515 | 1.473 | 0.998 | 0.986 | 0.808 | 0.929 | 0.299 | 0.181 | 0.163 |
| effMaxSG eps 0.20, HR 1.75 | field | 1.515 | 0.821 | 0.371 | 0.260 | 0.808 | 0.999 | 0.178 | 0.104 | 0.234 |
| effMaxSG eps 0.20, HR 1.75 | field-s | 1.515 | 0.821 | 0.371 | 0.260 | 0.808 | 1.010 | 0.156 | 0.084 | 0.246 |
| effMaxSG eps 0.20, HR 1.75 | IJ two-term | 1.515 | 0.827 | 0.383 | 0.243 | 0.808 | 1.207 | 0.019 | 0.006 | 0.425 |
| effMaxSG eps 0.30 (stress), HR 1.50 | naive | 1.228 | 1.173 | 0.950 | 0.835 | 0.779 | 0.880 | 0.428 | 0.284 | 0.162 |
| effMaxSG eps 0.30 (stress), HR 1.50 | field | 1.228 | 0.678 | 0.165 | 0.092 | 0.779 | 0.961 | 0.252 | 0.160 | 0.248 |
| effMaxSG eps 0.30 (stress), HR 1.50 | field-s | 1.228 | 0.678 | 0.165 | 0.092 | 0.779 | 0.977 | 0.214 | 0.129 | 0.266 |
| effMaxSG eps 0.30 (stress), HR 1.50 | IJ two-term | 1.228 | 0.672 | 0.127 | 0.056 | 0.779 | 1.159 | 0.028 | 0.011 | 0.438 |
| effMaxSG eps 0.30 (stress), HR 1.75 | naive | 1.454 | 1.337 | 0.990 | 0.950 | 0.774 | 0.894 | 0.402 | 0.261 | 0.170 |
| effMaxSG eps 0.30 (stress), HR 1.75 | field | 1.454 | 0.805 | 0.340 | 0.224 | 0.774 | 0.965 | 0.249 | 0.151 | 0.245 |
| effMaxSG eps 0.30 (stress), HR 1.75 | field-s | 1.454 | 0.805 | 0.340 | 0.224 | 0.774 | 0.981 | 0.210 | 0.128 | 0.262 |
| effMaxSG eps 0.30 (stress), HR 1.75 | IJ two-term | 1.454 | 0.786 | 0.311 | 0.181 | 0.774 | 1.171 | 0.028 | 0.008 | 0.440 |
| maxSG, HR 1.75 | naive | 1.014 | 0.852 | 0.351 | 0.082 | 0.758 | 1.258 | 0.330 | 0.257 | 0.457 |
| maxSG, HR 1.75 | field | 1.014 | 0.828 | 0.271 | 0.079 | 0.758 | 1.389 | 0.143 | 0.102 | 0.584 |
| maxSG, HR 1.75 | field-s | 1.014 | 0.828 | 0.271 | 0.079 | 0.758 | 1.451 | 0.144 | 0.103 | 0.611 |
| maxSG, HR 1.75 | IJ two-term | 1.014 | 0.706 | 0.015 | 0.001 | 0.758 | 2.321 | 0.019 | 0.009 | 1.022 |
| minSG, HR 1.75 | naive | 1.304 | 1.062 | 0.846 | 0.584 | 0.915 | 1.088 | 0.039 | 0.015 | 0.184 |
| minSG, HR 1.75 | field | 1.304 | 0.790 | 0.318 | 0.206 | 0.915 | 1.125 | 0.025 | 0.010 | 0.218 |
| minSG, HR 1.75 | field-s | 1.304 | 0.790 | 0.318 | 0.206 | 0.915 | 1.127 | 0.023 | 0.010 | 0.219 |
| minSG, HR 1.75 | IJ two-term | 1.304 | 0.738 | 0.209 | 0.118 | 0.915 | 1.368 | 0.001 | 0.000 | 0.413 |

**Joint pair** (coverage of (beta(Hhat) >= lower_H, beta(Hhat^c) <= upper_Hc), Wilson; margins in log units):

| cell | pair | n | joint cov [Wilson] | cov H | cov Hc | margin H | margin Hc | mean gamma | mean corr |
|---|---|---|---|---|---|---|---|---|---|
| effMaxSG eps 0.20, HR 1.50 | separate 95% (field) | 1999 | 0.872 [0.857, 0.886] | 0.974 | 0.897 | 0.624 | 0.237 | - | - |
| effMaxSG eps 0.20, HR 1.50 | separate 95% (field-s) | 1999 | 0.888 [0.873, 0.901] | 0.974 | 0.912 | 0.624 | 0.247 | - | - |
| effMaxSG eps 0.20, HR 1.50 | Bonferroni (joint) | 1999 | 0.932 [0.920, 0.942] | 0.989 | 0.942 | 0.762 | 0.280 | 0.025 | 0.041 |
| effMaxSG eps 0.20, HR 1.50 | Bonferroni (joint_s) | 1999 | 0.942 [0.931, 0.951] | 0.989 | 0.952 | 0.762 | 0.292 | 0.025 | 0.031 |
| effMaxSG eps 0.20, HR 1.50 | calibrated (joint) | 1999 | 0.932 [0.920, 0.942] | 0.989 | 0.942 | 0.761 | 0.279 | 0.025 | 0.041 |
| effMaxSG eps 0.20, HR 1.50 | calibrated (joint_s) | 1999 | 0.942 [0.931, 0.951] | 0.989 | 0.952 | 0.761 | 0.291 | 0.025 | 0.031 |
| effMaxSG eps 0.20, HR 1.75 | separate 95% (field) | 1999 | 0.885 [0.870, 0.898] | 0.970 | 0.912 | 0.611 | 0.234 | - | - |
| effMaxSG eps 0.20, HR 1.75 | separate 95% (field-s) | 1999 | 0.892 [0.878, 0.905] | 0.970 | 0.919 | 0.611 | 0.246 | - | - |
| effMaxSG eps 0.20, HR 1.75 | Bonferroni (joint) | 1999 | 0.933 [0.922, 0.944] | 0.986 | 0.947 | 0.747 | 0.276 | 0.025 | 0.019 |
| effMaxSG eps 0.20, HR 1.75 | Bonferroni (joint_s) | 1999 | 0.939 [0.928, 0.949] | 0.986 | 0.953 | 0.747 | 0.290 | 0.025 | 0.009 |
| effMaxSG eps 0.20, HR 1.75 | calibrated (joint) | 1999 | 0.933 [0.921, 0.943] | 0.986 | 0.946 | 0.746 | 0.276 | 0.025 | 0.019 |
| effMaxSG eps 0.20, HR 1.75 | calibrated (joint_s) | 1999 | 0.939 [0.928, 0.949] | 0.986 | 0.953 | 0.746 | 0.290 | 0.025 | 0.009 |
| effMaxSG eps 0.30 (stress), HR 1.50 | separate 95% (field) | 1999 | 0.886 [0.872, 0.900] | 0.980 | 0.904 | 0.555 | 0.248 | - | - |
| effMaxSG eps 0.30 (stress), HR 1.50 | separate 95% (field-s) | 1999 | 0.904 [0.890, 0.916] | 0.980 | 0.921 | 0.555 | 0.266 | - | - |
| effMaxSG eps 0.30 (stress), HR 1.50 | Bonferroni (joint) | 1999 | 0.935 [0.923, 0.945] | 0.993 | 0.941 | 0.691 | 0.294 | 0.025 | 0.039 |
| effMaxSG eps 0.30 (stress), HR 1.50 | Bonferroni (joint_s) | 1999 | 0.949 [0.939, 0.958] | 0.993 | 0.955 | 0.691 | 0.314 | 0.025 | 0.025 |
| effMaxSG eps 0.30 (stress), HR 1.50 | calibrated (joint) | 1999 | 0.935 [0.923, 0.945] | 0.993 | 0.941 | 0.690 | 0.294 | 0.025 | 0.039 |
| effMaxSG eps 0.30 (stress), HR 1.50 | calibrated (joint_s) | 1999 | 0.949 [0.939, 0.958] | 0.993 | 0.955 | 0.690 | 0.314 | 0.025 | 0.025 |
| effMaxSG eps 0.30 (stress), HR 1.75 | separate 95% (field) | 1999 | 0.881 [0.867, 0.895] | 0.976 | 0.903 | 0.542 | 0.245 | - | - |
| effMaxSG eps 0.30 (stress), HR 1.75 | separate 95% (field-s) | 1999 | 0.897 [0.883, 0.910] | 0.976 | 0.919 | 0.542 | 0.262 | - | - |
| effMaxSG eps 0.30 (stress), HR 1.75 | Bonferroni (joint) | 1999 | 0.935 [0.923, 0.945] | 0.992 | 0.943 | 0.676 | 0.290 | 0.025 | 0.012 |
| effMaxSG eps 0.30 (stress), HR 1.75 | Bonferroni (joint_s) | 1999 | 0.952 [0.942, 0.961] | 0.992 | 0.960 | 0.676 | 0.310 | 0.025 | -0.000 |
| effMaxSG eps 0.30 (stress), HR 1.75 | calibrated (joint) | 1999 | 0.935 [0.923, 0.945] | 0.992 | 0.943 | 0.675 | 0.290 | 0.025 | 0.012 |
| effMaxSG eps 0.30 (stress), HR 1.75 | calibrated (joint_s) | 1999 | 0.951 [0.941, 0.960] | 0.992 | 0.959 | 0.675 | 0.310 | 0.025 | -0.000 |
| maxSG, HR 1.75 | separate 95% (field) | 1999 | 0.876 [0.861, 0.890] | 0.947 | 0.925 | 0.204 | 0.584 | - | - |
| maxSG, HR 1.75 | separate 95% (field-s) | 1999 | 0.878 [0.863, 0.892] | 0.947 | 0.926 | 0.204 | 0.611 | - | - |
| maxSG, HR 1.75 | Bonferroni (joint) | 1999 | 0.942 [0.931, 0.952] | 0.975 | 0.966 | 0.243 | 0.710 | 0.025 | 0.035 |
| maxSG, HR 1.75 | Bonferroni (joint_s) | 1999 | 0.939 [0.928, 0.949] | 0.975 | 0.963 | 0.243 | 0.722 | 0.025 | 0.073 |
| maxSG, HR 1.75 | calibrated (joint) | 1999 | 0.941 [0.930, 0.951] | 0.975 | 0.965 | 0.241 | 0.706 | 0.026 | 0.035 |
| maxSG, HR 1.75 | calibrated (joint_s) | 1999 | 0.939 [0.928, 0.949] | 0.975 | 0.963 | 0.242 | 0.721 | 0.025 | 0.073 |
| minSG, HR 1.75 | separate 95% (field) | 1999 | 0.869 [0.854, 0.883] | 0.936 | 0.929 | 0.589 | 0.218 | - | - |
| minSG, HR 1.75 | separate 95% (field-s) | 1999 | 0.869 [0.853, 0.883] | 0.936 | 0.929 | 0.589 | 0.219 | - | - |
| minSG, HR 1.75 | Bonferroni (joint) | 1999 | 0.930 [0.918, 0.940] | 0.965 | 0.963 | 0.702 | 0.258 | 0.025 | 0.020 |
| minSG, HR 1.75 | Bonferroni (joint_s) | 1999 | 0.930 [0.918, 0.941] | 0.965 | 0.963 | 0.702 | 0.259 | 0.025 | 0.019 |
| minSG, HR 1.75 | calibrated (joint) | 1999 | 0.929 [0.917, 0.940] | 0.964 | 0.963 | 0.701 | 0.258 | 0.025 | 0.020 |
| minSG, HR 1.75 | calibrated (joint_s) | 1999 | 0.930 [0.918, 0.940] | 0.964 | 0.963 | 0.701 | 0.259 | 0.025 | 0.019 |

## Finding 5 -- se_field_s / (rho^c x se_field) per cell (R1 vs the global rescale; informational)

| cell | n | mean | sd | q01 | q10 | q50 | q90 | q99 | min | max | scale CV mean / q90 | mean rho^c | share rho^c > 1 | corr(se_s, nSE) / corr(se, nSE) | SD across reps: se_s / se / nSE |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| effMaxSG eps 0.20, HR 1.50 | 1999 | 0.9979 | 0.0036 | 0.9874 | 0.9933 | 0.9980 | 1.0021 | 1.0061 | 0.9831 | 1.0128 | 0.052 / 0.069 | 1.054 | 0.807 | 0.941 / 0.343 | 0.0105 / 0.0054 / 0.0104 |
| effMaxSG eps 0.20, HR 1.75 | 1999 | 0.9987 | 0.0035 | 0.9896 | 0.9943 | 0.9989 | 1.0029 | 1.0064 | 0.9833 | 1.0105 | 0.054 / 0.071 | 1.058 | 0.827 | 0.942 / 0.362 | 0.0106 / 0.0058 / 0.0102 |
| effMaxSG eps 0.30 (stress), HR 1.50 | 1999 | 0.9960 | 0.0057 | 0.9784 | 0.9888 | 0.9967 | 1.0022 | 1.0068 | 0.9626 | 1.0150 | 0.071 / 0.094 | 1.087 | 0.865 | 0.962 / 0.418 | 0.0149 / 0.0068 / 0.0148 |
| effMaxSG eps 0.30 (stress), HR 1.75 | 1999 | 0.9975 | 0.0050 | 0.9821 | 0.9910 | 0.9980 | 1.0031 | 1.0075 | 0.9700 | 1.0135 | 0.070 / 0.092 | 1.083 | 0.874 | 0.958 / 0.414 | 0.0136 / 0.0070 / 0.0132 |
| maxSG, HR 1.75 | 1999 | 0.9688 | 0.0319 | 0.8816 | 0.9252 | 0.9723 | 1.0055 | 1.0178 | 0.8125 | 1.0308 | 0.260 / 0.375 | 0.997 | 0.542 | 0.996 / 0.883 | 0.1271 / 0.0808 / 0.1192 |
| minSG, HR 1.75 | 1999 | 0.9999 | 0.0009 | 0.9974 | 0.9988 | 1.0000 | 1.0009 | 1.0020 | 0.9945 | 1.0040 | 0.015 / 0.022 | 1.006 | 0.651 | 0.741 / 0.511 | 0.0047 / 0.0044 / 0.0035 |

## Identification context per cell (unchanged from the comparators by pairing)

| cell | detected | mean |Hhat| | |Hhat|/|H| median | sens | spec | p-hat mean | complement fits/rep | fit+MR s/rep (e1stud) |
|---|---|---|---|---|---|---|---|---|
| effMaxSG eps 0.20, HR 1.50 | 1999 | 127.1 | 0.793 | 0.590 | 0.893 | 0.103 | 615 | 72.7 |
| effMaxSG eps 0.20, HR 1.75 | 1999 | 128.9 | 0.821 | 0.660 | 0.919 | 0.119 | 580 | 76.0 |
| effMaxSG eps 0.30 (stress), HR 1.50 | 1999 | 168.9 | 1.093 | 0.737 | 0.838 | 0.064 | 749 | 73.4 |
| effMaxSG eps 0.30 (stress), HR 1.75 | 1999 | 166.8 | 1.087 | 0.792 | 0.868 | 0.080 | 708 | 76.7 |
| maxSG, HR 1.75 | 1999 | 401.6 | 2.752 | 0.948 | 0.260 | 0.415 | 465 | 55.7 |
| minSG, HR 1.75 | 1999 | 61.8 | 0.404 | 0.249 | 0.931 | 0.086 | 98 | 53.5 |
