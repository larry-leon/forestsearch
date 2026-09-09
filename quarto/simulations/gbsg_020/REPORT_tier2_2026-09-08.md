# REPORT — Tier 2 (Mac Studio): the finalized constructions at 12.5% prevalence, the dominated-complement regime (campaign `tier2`) — Gate 2 records, the standard tables and the pre-registered criteria (report-and-wait)

**Task:** `dev/tasks/TASK_tier2_mac_2026-09-08.md` (2a84f6c4), Stages 2–3. Stage 0/1 record: `REPORT_tier2_stage1_2026-09-08.md` (b9b4b5ed; smoke PASS, Gate 1 compute go). Predecessors: `REPORT_mr_field_complement_2026-09-06.md` and `REPORT_complement_refinements_2026-09-06.md` (the `s7c` / `map1c` unscaled-field records at this prevalence), `REPORT_field_studentize_e1_2026-09-08.md`, `dev/notes/HANDOFF_mr_field_linux_2026-09-08.md`.
**Date:** 2026-09-08. Executor: Claude Code on **Mac-Studio-3** (`Mac-Studio-3.local`, M4 Max, 14 physical cores, 36 GiB), unattended; ran in parallel with the Linux `cert20` campaign, outputs disjoint. **No `R/`, template or shared-document change of any kind.** Winner-only and winner-floor excluded from every table and line. Bounds read by location; Wilson intervals; marginal and error SDs side by side. **No recommendation change; report and wait.**

---

## GATE 2: PASS on all four cells — 2,000 rows each, sim_id 1–2000, no duplicates, no CONFIG-ERROR, meta knobs as set (`field_complement TRUE`, `field_decompose TRUE`, `field_scale_complement "selected"`, `ij_residual "two_term"`, J = 10, `harm_z1_quantile 0.25`, seed 8316951, `forestsearch 0.3.5`, `Mac-Studio-3.local` on every batch); every harm / complement / `_s` / joint / β(Ĥ) / β(Ĥᶜ) / p̂ / K / ρᶜ quantity finite on detected replicates; interval invariants hold on all three blocks; γ ∈ [0.025, 0.028] on both the joint and `joint_s` pairs; **zero replicates below the achieved-probability bar 0.95 − 2/n_joint** (tested per replicate with that replicate's own `n_joint`); bound↔quantile identities **0** to machine precision; realized prevalence 0.1236–0.1240 against the M1 super-population 0.1242. **Campaign wall 2 h 04 m (7,442 s) at 12 workers; four of four cells completed, none deferred, none dropped** (ceiling 9 h, hard timeout 11 h). Peak memory 21.6 GB total RSS over 13 R processes, inside the 24 GB rule.

## Run

Committed template (unchanged, HEAD `2a84f6c4`; last `R/` commit `c0f48a7c`) driven by environment only, every knob explicit:
`FS_S7_FOCUS=maxeffCons FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_CAMPAIGN=tier2 FS_S7_WORKERS=12`, per cell `FS_S7_HR` / `FS_S7_N`, two seed-disjoint batches (`FS_S7_START` 1 and 1001, `FS_S7_NSIMS=1000`) then `FS_S7_MODE=combine`; seeds 8316951 + sim_id; J = 10 default; `return_reselection = TRUE`. **`FS_S7_Z1Q` was never set**, so the M1 default 0.25 (12.4 % prevalence) stands and the committed `s7c` / `map1c` configuration is reproduced. Single-threaded BLAS (`VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`). Driver `tier2_campaign.sh` (session scratchpad; re-projects each cell before starting it and defers any that would cross the hard timeout — none did). `.refuse_if_tracked()` live throughout. Stems `fs_maxeffCons_fb_mr_field_m1_h{175,175,175,100}_knoise0_n{500,1000,1500,500}_tier2`; per cell two batch renders and one combine render committed beside the bundles. Cross-cell document `summary_tier2.qmd`, **transplanted from the committed `summary_e1stud.qmd`** (bundle globs, cell labels, comparator names, the absent-comparator guard, the DGM-draw pairing table and the 0.93 threshold constant — no fresh authorship of the table machinery), rendered to `summary_tier2.html`. The tables below are the verbatim output of `tier2_findings.R` (this directory), which implements the same definitions.

| # | cell | batch 1 | batch 2 | combine | cell wall | comparator |
|---|---|---|---|---|---|---|
| 1 | HR 1.75, n = 500 | 920 s | 893 s | 10 s | **1,823 s (30 m 23 s)** | `s7c` h175 n500 |
| 2 | HR 1.75, n = 1000 | 1,014 s | 1,037 s | 9 s | **2,060 s (34 m 20 s)** | none at this cell |
| 3 | HR 1.75, n = 1500 | 1,045 s | 1,048 s | 9 s | **2,102 s (35 m 02 s)** | none exact (`map1c` h150 n1500 is a different HR) |
| 4 | HR 1.00, n = 500 (null) | 690 s | 757 s | 10 s | **1,457 s (24 m 17 s)** | `s7c` null |
| | **campaign** | | | | **7,442 s (2 h 04 m 02 s)** | Gate 1 projection was ≈ 2 h 45 m |

## The cross-machine pairing, stated plainly

The Gate 2 pairing to a comparator is by **DGM draws**: `n_true` is `identical()` on all 2,000 rows at both cells that have one (HR 1.75 n500 against `s7c`, and the null against `s7c` null). **`truth` is *not* `identical()`** — it agrees to a maximum relative difference of **6.55e−16** (HR 1.75 n500) and **3.80e−16** (null), the cross-machine BLAS-precision difference the task's own note anticipates; per that note the exact-identity assertion is not applied to it, and the ~1e−8 tolerance is met by twelve orders of magnitude. Informationally (never asserted): over the paired detected rows the Mac reproduces the committed comparator's fitted values to **max |Δ `fld_Hc_up1s`| = 6.66e−15 / 5.88e−15** and **max |Δ β(Ĥᶜ)| = 1.11e−16 / 2.22e−16**, and the unscaled field's coverage matches the comparator's to four decimals (**0.9358 vs 0.9358**; **0.9302 vs 0.9302**).

## The findings, stated as findings

**1. The dominated regime holds, cell by cell.** ρᶜ is essentially 1 everywhere: mean **1.0166 / 1.0105 / 1.0095 / 1.0196** with q10–q90 **[0.991, 1.044] / [0.989, 1.031] / [0.992, 1.026] / [0.988, 1.053]** and share > 1 of 0.76–0.80 — the complement is dominated by the selected candidate's own scale, exactly as pre-registered, and the spread narrows with n. λ-SDᶜ/naive SE is **0.982 / 0.995 / 1.005 / 0.975** and `se_field_s`/naive SE **0.998 / 1.006 / 1.015 / 0.993**, i.e. both field SDs already sit at the naive scale here (contrast the 31 % band cells, where the unscaled ratio was 0.917–0.947). **field-s ≈ field in every cell**: +0.0047, +0.0015, +0.0010, +0.0096, with the Wilson intervals overlapping at all four and flips almost one-directional (to cover 0.005 / 0.002 / 0.001 / 0.011 against to miss 0.000 / 0.0005 / 0.000 / 0.0015). The pre-registered expectation for this regime is met.

**2. Against the committed records (0.930–0.956 at this prevalence).** The field's complement upper coverage is **0.936 / 0.955 / 0.960 / 0.930** across the four cells — inside or at the top of the committed band, and rising monotonically with n at HR 1.75. Where a comparator exists the reproduction is exact to four decimals (above). The two larger-n cells have no comparator on record at this HR; they are new.

**3. The stable-pick (high-p̂) caveat is present but shrinks with n, and it moves between blocks.** By p̂ tertile the complement field falls **0.954 → 0.946 → 0.906** at n = 500, **0.968 → 0.958 → 0.938** at n = 1000 and **0.979 → 0.958 → 0.941** at n = 1500: the T3 shortfall against T1 narrows from −4.8 to −3.0 to −3.8 points, and the T3 level itself rises from 0.906 to 0.941. The harm block shows the same shape and does **not** shrink the same way — **0.983 / 0.986 / 0.930** at n = 500, **0.983 / 0.959 / 0.888** at n = 1000, **0.982 / 0.964 / 0.934** at n = 1500 — the n = 1000 harm T3 (0.888 [0.862, 0.910]) is the weakest stratum in the campaign. Gaussian-implied tracks observed closely on both blocks: over the 24 harm/complement tertile pairs, 16 agree within 0.010 and 22 within 0.020, the maximum gap being 0.037 at the null cell's harm T3 (0.916 observed vs 0.879 implied); the complement pairs agree within 0.014 everywhere. So the shortfall is a location/scale effect the Gaussian reference already predicts, not an artefact of the tail.

**4. Identification improves sharply with n, and the complement's target is stable.** Detection **0.950 / 0.995 / 0.999** at HR 1.75 and 0.680 at the null; |Ĥ|/|H| median **1.112 / 1.000 / 1.028** (q90 1.507 / 1.148 / 1.132); sens **0.659 / 0.787 / 0.883**, spec 0.930 / 0.972 / 0.983, PPV 0.590 / 0.794 / 0.877, NPV 0.951 / 0.970 / 0.984. β(Ĥ) climbs toward the planted marginal 1.769 (**1.253 / 1.523 / 1.629**; planted CDE 2.036) while β(Ĥᶜ) sits on the planted 0.657 throughout (**0.669 / 0.653 / 0.642**) — the harm block's target is the one that moves with identification quality, the complement's is not. p̂ mean rises 0.234 → 0.374 → 0.498; complement fits per replicate 382 / 346 / 308 with 10 % / 19 % / 29 % of draws needing a new fit.

**5. The null cell (recorded without criterion).** Detection 0.680 on 2,000 replicates; the naive harm interval covers 0.051 of the time (b = 5.00) where the field covers 0.963 and IJ 0.966; the complement field 0.930 and field-s 0.940 (the largest field-s gain in the campaign, +0.0096); joint Bonferroni 0.955.

## The pre-registered acceptance criteria

| cell | harm field one-sided lower (≥ 0.94 **with Wilson support**) | field-s complement upper (≥ 0.93) | IJ two-sided (≥ 0.93) | joint Bonferroni (≥ 0.93) |
|---|---|---|---|---|
| HR 1.75, n = 500 | 0.9663 [0.9572, 0.9741] — **MET** | 0.9405 [0.9290, 0.9503] — **MET** | H 0.9611 / Hᶜ 1.0000 — **MET** | 0.9558 (joint_s 0.9568) — **MET** |
| HR 1.75, n = 1000 | 0.9437 [0.9327, 0.9531] — **NOT MET** (point ≥ 0.94, Wilson lower 0.9327 < 0.94) | 0.9563 [0.9464, 0.9644] — **MET** | H **0.9166** / Hᶜ 0.9995 — **NOT MET** | 0.9412 (joint_s 0.9412) — **MET** |
| HR 1.75, n = 1500 | 0.9600 [0.9504, 0.9678] — **MET** | 0.9605 [0.9510, 0.9682] — **MET** | H **0.9129** / Hᶜ 0.9995 — **NOT MET** | 0.9635 (joint_s 0.9640) — **MET** |
| HR 1.00, n = 500 (null) — **recorded, no criterion** | 0.9625 [0.9511, 0.9713] | 0.9398 [0.9258, 0.9512] | H 0.9860 / Hᶜ 1.0000 | 0.9552 (joint_s 0.9566) |

**Two criteria are not met, both on the harm block and both at HR 1.75 with n ≥ 1000.** (i) The harm field's one-sided lower coverage at n = 1000 is 0.9437 — the point estimate clears 0.94, the Wilson lower limit (0.9327) does not, so the criterion as written ("with Wilson support") fails by 0.007; at n = 500 and n = 1500 it clears comfortably. (ii) The IJ two-term **two-sided** coverage on the harm block is 0.9166 (n = 1000) and 0.9129 (n = 1500) against the 0.93 bar, while at n = 500 it is 0.9611 and at the null 0.9860; on the complement block IJ two-sided is 0.9995–1.0000 everywhere. The IJ shortfall tracks a residual harm-block bias that does not vanish with n (field bias −0.075 / −0.091 / −0.075 log units; IJ −0.010 / −0.081 / −0.093), so the two-sided interval — which is the only construction penalised on both sides — loses coverage as n grows even as the one-sided lower bound holds. **Nothing here is a repair proposal.** The field-s complement criterion (≥ 0.93) is **met at every cell of this prevalence**, which is the finalized product's own bar, and the joint Bonferroni criterion is met at every cell.

## What this record does and does not say

Reported: Gate 2 passes on four of four cells with the DGM-draw pairing exact where a comparator exists; the dominated regime is confirmed (ρᶜ ≈ 1.01–1.02, both field SDs at the naive scale) and field-s is inert against field there, gaining 0.1–1.0 points with overlapping Wilson intervals; the field's complement coverage 0.930–0.960 reproduces and extends the committed 0.930–0.956 band, rising with n; the stable-pick high-p̂ shortfall is present on both blocks, shrinks with n on the complement and does not on the harm block; field-s ≥ 0.93 and joint Bonferroni ≥ 0.93 hold at every cell; the harm field's ≥ 0.94-with-Wilson-support criterion fails at n = 1000 only, and IJ two-sided on the harm block falls below 0.93 at n = 1000 and n = 1500. Not said: whether any construction is adopted, whether any documented rule changes, whether the two unmet criteria warrant a repair, or what the next tier should be. **Report and wait.**

---

# Tables (verbatim output of `tier2_findings.R`)

## Gate 2 per cell

| cell | rows | sim_id | dups | config_err | detected | prevalence | all_finite | invariants | gamma_range | gamma_s_range | joint_n_range | joint_below_bar | joint_s_below_bar | bound_id_max | n_true_identical | truth_max_reldiff | machine | version | knobs |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HR 1.75, n = 500 | 2000.000000 | 1-2000 | 0.000000 | 0.000000 | 1900.000000 | 0.123628 | TRUE | TRUE | [0.0250, 0.0270] | [0.0250, 0.0270] | 766-1000 | 0.000000 | 0.000000 | 0.000000 | TRUE | 6.55e-16 | Mac-Studio-3.local | 0.3.5 | TRUE/TRUE/selected/two_term/J=10/z1q=0.25 |
| HR 1.75, n = 1000 | 2000.000000 | 1-2000 | 0.000000 | 0.000000 | 1990.000000 | 0.123876 | TRUE | TRUE | [0.0250, 0.0280] | [0.0250, 0.0280] | 707-1000 | 0.000000 | 0.000000 | 0.000000 | NA | NA | Mac-Studio-3.local | 0.3.5 | TRUE/TRUE/selected/two_term/J=10/z1q=0.25 |
| HR 1.75, n = 1500 | 2000.000000 | 1-2000 | 0.000000 | 0.000000 | 1998.000000 | 0.124012 | TRUE | TRUE | [0.0250, 0.0270] | [0.0250, 0.0270] | 888-1000 | 0.000000 | 0.000000 | 0.000000 | NA | NA | Mac-Studio-3.local | 0.3.5 | TRUE/TRUE/selected/two_term/J=10/z1q=0.25 |
| HR 1.00, n = 500 (null) | 2000.000000 | 1-2000 | 0.000000 | 0.000000 | 1361.000000 | 0.123628 | TRUE | TRUE | [0.0250, 0.0270] | [0.0250, 0.0270] | 697-1000 | 0.000000 | 0.000000 | 0.000000 | TRUE | 3.80e-16 | Mac-Studio-3.local | 0.3.5 | TRUE/TRUE/selected/two_term/J=10/z1q=0.25 |


## 1. Constructions per cell (identical replicates): naive / field / field-s / IJ two-term

| cell | block | construction | n | bias_log | sd_emp | sd_err | se_mean | b | r | se_over_sd_err | cov1_wilson | cov2 | cov1_ref |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HR 1.75, n = 500 | Hhat (lower) | naive | 1900.000 | 0.481 | 0.222 | 0.349 | 0.310 | 2.165 | 1.395 | 0.888 | 0.539 [0.516, 0.561] | 0.645 | 0.551 |
| HR 1.75, n = 500 | Hhat (lower) | field | 1900.000 | -0.075 | 0.358 | 0.399 | 0.363 | -0.209 | 1.014 | 0.910 | 0.966 [0.957, 0.974] | 0.887 | 0.970 |
| HR 1.75, n = 500 | Hhat (lower) | IJ two-term | 1900.000 | -0.010 | 0.312 | 0.380 | 0.393 | -0.032 | 1.260 | 1.034 | 0.983 [0.976, 0.988] | 0.961 | 0.982 |
| HR 1.75, n = 500 | Hhat^c (upper) | naive | 1900.000 | -0.070 | 0.138 | 0.137 | 0.139 | -0.511 | 1.007 | 1.011 | 0.876 [0.861, 0.890] | 0.921 | 0.874 |
| HR 1.75, n = 500 | Hhat^c (upper) | field | 1900.000 | -0.001 | 0.147 | 0.144 | 0.136 | -0.004 | 0.927 | 0.948 | 0.936 [0.924, 0.946] | 0.934 | 0.936 |
| HR 1.75, n = 500 | Hhat^c (upper) | field-s | 1900.000 | -0.001 | 0.147 | 0.144 | 0.138 | -0.005 | 0.942 | 0.962 | 0.941 [0.929, 0.950] | 0.938 | 0.939 |
| HR 1.75, n = 500 | Hhat^c (upper) | IJ two-term | 1900.000 | -0.008 | 0.143 | 0.141 | 0.257 | -0.058 | 1.800 | 1.825 | 0.999 [0.996, 1.000] | 1.000 | 0.998 |
| HR 1.75, n = 1000 | Hhat (lower) | naive | 1990.000 | 0.213 | 0.202 | 0.291 | 0.228 | 1.056 | 1.134 | 0.786 | 0.757 [0.737, 0.775] | 0.824 | 0.791 |
| HR 1.75, n = 1000 | Hhat (lower) | field | 1990.000 | -0.091 | 0.323 | 0.339 | 0.280 | -0.282 | 0.867 | 0.826 | 0.944 [0.933, 0.953] | 0.850 | 0.956 |
| HR 1.75, n = 1000 | Hhat (lower) | IJ two-term | 1990.000 | -0.081 | 0.285 | 0.323 | 0.302 | -0.283 | 1.059 | 0.933 | 0.975 [0.968, 0.981] | 0.917 | 0.979 |
| HR 1.75, n = 1000 | Hhat^c (upper) | naive | 1990.000 | -0.026 | 0.099 | 0.096 | 0.097 | -0.258 | 0.979 | 1.006 | 0.926 [0.914, 0.937] | 0.949 | 0.912 |
| HR 1.75, n = 1000 | Hhat^c (upper) | field | 1990.000 | 0.008 | 0.106 | 0.101 | 0.096 | 0.080 | 0.913 | 0.958 | 0.955 [0.945, 0.963] | 0.936 | 0.943 |
| HR 1.75, n = 1000 | Hhat^c (upper) | field-s | 1990.000 | 0.008 | 0.106 | 0.101 | 0.098 | 0.078 | 0.922 | 0.968 | 0.956 [0.946, 0.964] | 0.940 | 0.945 |
| HR 1.75, n = 1000 | Hhat^c (upper) | IJ two-term | 1990.000 | 0.007 | 0.103 | 0.099 | 0.184 | 0.071 | 1.784 | 1.855 | 0.998 [0.996, 0.999] | 0.999 | 0.999 |
| HR 1.75, n = 1500 | Hhat (lower) | naive | 1998.000 | 0.083 | 0.170 | 0.208 | 0.179 | 0.488 | 1.051 | 0.860 | 0.881 [0.866, 0.894] | 0.910 | 0.893 |
| HR 1.75, n = 1500 | Hhat (lower) | field | 1998.000 | -0.075 | 0.262 | 0.261 | 0.231 | -0.286 | 0.882 | 0.887 | 0.960 [0.950, 0.968] | 0.878 | 0.959 |
| HR 1.75, n = 1500 | Hhat (lower) | IJ two-term | 1998.000 | -0.093 | 0.236 | 0.247 | 0.258 | -0.393 | 1.090 | 1.045 | 0.989 [0.983, 0.993] | 0.913 | 0.986 |
| HR 1.75, n = 1500 | Hhat^c (upper) | naive | 1998.000 | -0.013 | 0.079 | 0.078 | 0.080 | -0.162 | 1.010 | 1.027 | 0.935 [0.923, 0.945] | 0.957 | 0.933 |
| HR 1.75, n = 1500 | Hhat^c (upper) | field | 1998.000 | 0.005 | 0.084 | 0.081 | 0.080 | 0.061 | 0.955 | 0.988 | 0.959 [0.950, 0.967] | 0.950 | 0.949 |
| HR 1.75, n = 1500 | Hhat^c (upper) | field-s | 1998.000 | 0.005 | 0.084 | 0.081 | 0.081 | 0.059 | 0.964 | 0.998 | 0.960 [0.951, 0.968] | 0.951 | 0.950 |
| HR 1.75, n = 1500 | Hhat^c (upper) | IJ two-term | 1998.000 | 0.007 | 0.082 | 0.080 | 0.153 | 0.089 | 1.860 | 1.910 | 0.999 [0.996, 1.000] | 0.999 | 0.999 |
| HR 1.00, n = 500 (null) | Hhat (lower) | naive | 1361.000 | 0.813 | 0.163 | 0.211 | 0.324 | 4.997 | 1.989 | 1.538 | 0.051 [0.040, 0.064] | 0.157 | 0.042 |
| HR 1.00, n = 500 (null) | Hhat (lower) | field | 1361.000 | 0.140 | 0.253 | 0.279 | 0.357 | 0.556 | 1.410 | 1.279 | 0.963 [0.951, 0.971] | 0.988 | 0.961 |
| HR 1.00, n = 500 (null) | Hhat (lower) | IJ two-term | 1361.000 | 0.244 | 0.222 | 0.255 | 0.399 | 1.098 | 1.794 | 1.564 | 0.966 [0.955, 0.975] | 0.986 | 0.968 |
| HR 1.00, n = 500 (null) | Hhat^c (upper) | naive | 1361.000 | -0.090 | 0.131 | 0.129 | 0.139 | -0.688 | 1.060 | 1.078 | 0.859 [0.839, 0.876] | 0.911 | 0.855 |
| HR 1.00, n = 500 (null) | Hhat^c (upper) | field | 1361.000 | -0.010 | 0.141 | 0.139 | 0.136 | -0.069 | 0.966 | 0.980 | 0.930 [0.915, 0.943] | 0.950 | 0.936 |
| HR 1.00, n = 500 (null) | Hhat^c (upper) | field-s | 1361.000 | -0.010 | 0.141 | 0.139 | 0.138 | -0.069 | 0.983 | 0.998 | 0.940 [0.926, 0.951] | 0.952 | 0.939 |
| HR 1.00, n = 500 (null) | Hhat^c (upper) | IJ two-term | 1361.000 | -0.021 | 0.137 | 0.135 | 0.255 | -0.156 | 1.854 | 1.882 | 0.997 [0.992, 0.999] | 1.000 | 0.998 |


**Across cells** (one-sided coverage min / mean / max; mean two-sided; mean b; mean r; mean SE/error SD):

| block | construction | cells | cov1_min | cov1_mean | cov1_max | cov2_mean | b_mean | r_mean | se_over_sd_err_mean |
|---|---|---|---|---|---|---|---|---|---|
| Hhat (lower) | field | 4.000 | 0.944 | 0.958 | 0.966 | 0.901 | -0.055 | 1.043 | 0.975 |
| Hhat (lower) | IJ two-term | 4.000 | 0.966 | 0.978 | 0.989 | 0.944 | 0.097 | 1.301 | 1.144 |
| Hhat (lower) | naive | 4.000 | 0.051 | 0.557 | 0.881 | 0.634 | 2.177 | 1.392 | 1.018 |
| Hhat^c (upper) | field | 4.000 | 0.930 | 0.945 | 0.959 | 0.942 | 0.017 | 0.940 | 0.969 |
| Hhat^c (upper) | field-s | 4.000 | 0.940 | 0.949 | 0.960 | 0.945 | 0.016 | 0.953 | 0.981 |
| Hhat^c (upper) | IJ two-term | 4.000 | 0.997 | 0.998 | 0.999 | 1.000 | -0.013 | 1.825 | 1.868 |
| Hhat^c (upper) | naive | 4.000 | 0.859 | 0.899 | 0.935 | 0.934 | -0.405 | 1.014 | 1.030 |


## 2. The dominated-regime check: rho^c, the SD ratios, field vs field-s

| cell | n | rho_c_mean | rho_c_q10 | rho_c_q90 | rho_c_share_gt1 | lamSD_over_nSE | lamSDs_over_nSE | field | field_lo | field_hi | field_s | field_s_lo | field_s_hi | diff | wilson_overlap | flip_to_cover | flip_to_miss |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HR 1.75, n = 500 | 1900.0000 | 1.0166 | 0.9912 | 1.0438 | 0.7895 | 0.9824 | 0.9984 | 0.9358 | 0.9239 | 0.9460 | 0.9405 | 0.9290 | 0.9503 | 0.0047 | TRUE | 0.0047 | 0.0000 |
| HR 1.75, n = 1000 | 1990.0000 | 1.0105 | 0.9894 | 1.0310 | 0.7598 | 0.9949 | 1.0057 | 0.9548 | 0.9447 | 0.9631 | 0.9563 | 0.9464 | 0.9644 | 0.0015 | TRUE | 0.0020 | 0.0005 |
| HR 1.75, n = 1500 | 1998.0000 | 1.0095 | 0.9921 | 1.0264 | 0.8013 | 1.0048 | 1.0150 | 0.9595 | 0.9499 | 0.9673 | 0.9605 | 0.9510 | 0.9682 | 0.0010 | TRUE | 0.0010 | 0.0000 |
| HR 1.00, n = 500 (null) | 1361.0000 | 1.0196 | 0.9884 | 1.0532 | 0.7627 | 0.9751 | 0.9931 | 0.9302 | 0.9154 | 0.9426 | 0.9398 | 0.9258 | 0.9512 | 0.0096 | TRUE | 0.0110 | 0.0015 |


## 3. Against the committed unscaled-field comparators (pairing by DGM draws, not fitted values)

| cell | comparator | n | tier2_field | tier2_lo | tier2_hi | comp_field | comp_lo | comp_hi | diff | wilson_overlap | max_absdiff_up1s | max_absdiff_betaHc |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HR 1.75, n = 500 | s7c | 1900.0000 | 0.9358 | 0.9239 | 0.9460 | 0.9358 | 0.9239 | 0.9460 | 0.0000 | TRUE | 6.66e-15 | 1.11e-16 |
| HR 1.75, n = 1000 | none at this cell | 1990.0000 | 0.9548 | 0.9447 | 0.9631 |    NA |    NA |    NA |    NA | NA | NA | NA |
| HR 1.75, n = 1500 | none exact (map1c h150 n1500 is a different HR) | 1998.0000 | 0.9595 | 0.9499 | 0.9673 |    NA |    NA |    NA |    NA | NA | NA | NA |
| HR 1.00, n = 500 (null) | s7c null | 1361.0000 | 0.9302 | 0.9154 | 0.9426 | 0.9302 | 0.9154 | 0.9426 | 0.0000 | TRUE | 5.88e-15 | 2.22e-16 |


## 4. By p-hat tertile, both blocks (observed with Wilson; Gaussian-implied beside)

| cell | tertile | n | p_hat_mean | rho_c | harm_field | harm_lo | harm_hi | harm_gauss | comp_field | comp_lo | comp_hi | comp_gauss | comp_field_s | comp_s_lo | comp_s_hi | comp_s_gauss |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HR 1.75, n = 500 | T1 [0.000, 0.138] | 634.000 | 0.080 | 1.028 | 0.983 | 0.969 | 0.990 | 0.962 | 0.954 | 0.935 | 0.968 | 0.961 | 0.964 | 0.946 | 0.976 | 0.965 |
| HR 1.75, n = 500 | T2 [0.138, 0.268] | 635.000 | 0.196 | 1.014 | 0.986 | 0.973 | 0.993 | 0.980 | 0.946 | 0.926 | 0.961 | 0.955 | 0.951 | 0.932 | 0.965 | 0.957 |
| HR 1.75, n = 500 | T3 [0.268, 0.991] | 631.000 | 0.426 | 1.008 | 0.930 | 0.908 | 0.948 | 0.927 | 0.906 | 0.881 | 0.927 | 0.903 | 0.906 | 0.881 | 0.927 | 0.905 |
| HR 1.75, n = 1000 | T1 [0.002, 0.261] | 665.000 | 0.173 | 1.018 | 0.983 | 0.971 | 0.991 | 0.983 | 0.968 | 0.952 | 0.979 | 0.967 | 0.970 | 0.954 | 0.980 | 0.969 |
| HR 1.75, n = 1000 | T2 [0.262, 0.431] | 662.000 | 0.341 | 1.009 | 0.959 | 0.941 | 0.972 | 0.977 | 0.958 | 0.940 | 0.971 | 0.960 | 0.958 | 0.940 | 0.971 | 0.961 |
| HR 1.75, n = 1000 | T3 [0.431, 0.968] | 663.000 | 0.610 | 1.005 | 0.888 | 0.862 | 0.910 | 0.902 | 0.938 | 0.917 | 0.954 | 0.930 | 0.941 | 0.921 | 0.957 | 0.931 |
| HR 1.75, n = 1500 | T1 [0.007, 0.377] | 666.000 | 0.258 | 1.016 | 0.982 | 0.969 | 0.990 | 0.989 | 0.979 | 0.965 | 0.987 | 0.972 | 0.979 | 0.965 | 0.987 | 0.974 |
| HR 1.75, n = 1500 | T2 [0.377, 0.592] | 666.000 | 0.483 | 1.008 | 0.964 | 0.947 | 0.976 | 0.978 | 0.958 | 0.940 | 0.971 | 0.957 | 0.959 | 0.942 | 0.972 | 0.958 |
| HR 1.75, n = 1500 | T3 [0.594, 0.996] | 666.000 | 0.753 | 1.005 | 0.934 | 0.912 | 0.950 | 0.947 | 0.941 | 0.921 | 0.957 | 0.940 | 0.943 | 0.923 | 0.958 | 0.941 |
| HR 1.00, n = 500 (null) | T1 [0.000, 0.090] | 456.000 | 0.046 | 1.040 | 0.978 | 0.960 | 0.988 | 0.962 | 0.947 | 0.923 | 0.964 | 0.961 | 0.965 | 0.944 | 0.978 | 0.966 |
| HR 1.00, n = 500 (null) | T2 [0.091, 0.189] | 451.000 | 0.136 | 1.016 | 0.993 | 0.981 | 0.998 | 0.988 | 0.951 | 0.927 | 0.968 | 0.960 | 0.958 | 0.935 | 0.973 | 0.962 |
| HR 1.00, n = 500 (null) | T3 [0.190, 0.797] | 454.000 | 0.293 | 1.003 | 0.916 | 0.887 | 0.938 | 0.879 | 0.892 | 0.860 | 0.917 | 0.892 | 0.896 | 0.865 | 0.921 | 0.892 |


## 5. Identification and n

| cell | n | rows | detection | n_true_mean | n_sel_mean | relsize_median | relsize_q90 | sens | spec | ppv | npv | betaHhat_H_mean | planted_marg_H | planted_cde_H | betaHhat_Hc_mean | planted_marg_Hc | planted_cde_Hc | p_hat_mean | n_family_mean | comp_fits_mean | share_new_fit |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HR 1.75, n = 500 | 500.000 | 2000.000 | 0.950 | 61.814 | 72.030 | 1.112 | 1.507 | 0.659 | 0.930 | 0.590 | 0.951 | 1.253 | 1.769 | 2.036 | 0.669 | 0.657 | 0.585 | 0.234 | 1228.865 | 382.024 | 0.104 |
| HR 1.75, n = 1000 | 1000.000 | 2000.000 | 0.995 | 123.876 | 122.507 | 1.000 | 1.148 | 0.787 | 0.972 | 0.794 | 0.970 | 1.523 | 1.769 | 2.036 | 0.653 | 0.657 | 0.585 | 0.374 | 1323.515 | 346.443 | 0.190 |
| HR 1.75, n = 1500 | 1500.000 | 2000.000 | 0.999 | 186.018 | 186.437 | 1.028 | 1.132 | 0.883 | 0.983 | 0.877 | 0.984 | 1.629 | 1.769 | 2.036 | 0.642 | 0.657 | 0.585 | 0.498 | 1324.148 | 307.572 | 0.288 |
| HR 1.00, n = 500 (null) | 500.000 | 2000.000 | 0.680 | 61.814 | 78.904 | 1.194 | 1.741 | 0.352 | 0.870 | 0.292 | 0.904 | 0.719 | 1.000 | 1.000 | 0.661 | 0.657 | 0.585 | 0.158 | 1228.161 | 414.744 | 0.048 |


## 6. The pre-registered acceptance criteria, evaluated as findings

| cell | n | harm_field | harm_wilson_lo | harm_met | field_s | field_s_wilson_lo | field_s_met | ij2_H | ij2_Hc | ij2_met | joint_bonf | joint_s_bonf | joint_met |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| HR 1.75, n = 500 | 1900.0000 | 0.9663 | 0.9572 | MET | 0.9405 | 0.9290 | MET | 0.9611 | 1.0000 | MET | 0.9558 | 0.9568 | MET |
| HR 1.75, n = 1000 | 1990.0000 | 0.9437 | 0.9327 | NOT MET | 0.9563 | 0.9464 | MET | 0.9166 | 0.9995 | NOT MET | 0.9412 | 0.9412 | MET |
| HR 1.75, n = 1500 | 1998.0000 | 0.9600 | 0.9504 | MET | 0.9605 | 0.9510 | MET | 0.9129 | 0.9995 | NOT MET | 0.9635 | 0.9640 | MET |
| HR 1.00, n = 500 (null) | 1361.0000 | 0.9625 | 0.9511 | recorded (null) | 0.9398 | 0.9258 | recorded (null) | 0.9860 | 1.0000 | recorded (null) | 0.9552 | 0.9566 | recorded (null) |


## Reading lines

- HR 1.75, n = 500 (n 1900): rho^c mean 1.0166 [q10 0.9912, q90 1.0438], share>1 0.789; lambda-SD^c/naive 0.9824, se_field_s/naive 0.9984; field 0.936 [0.924, 0.946] vs field-s 0.941 [0.929, 0.950] (diff +0.0047, Wilson overlap TRUE)
- HR 1.75, n = 1000 (n 1990): rho^c mean 1.0105 [q10 0.9894, q90 1.0310], share>1 0.760; lambda-SD^c/naive 0.9949, se_field_s/naive 1.0057; field 0.955 [0.945, 0.963] vs field-s 0.956 [0.946, 0.964] (diff +0.0015, Wilson overlap TRUE)
- HR 1.75, n = 1500 (n 1998): rho^c mean 1.0095 [q10 0.9921, q90 1.0264], share>1 0.801; lambda-SD^c/naive 1.0048, se_field_s/naive 1.0150; field 0.959 [0.950, 0.967] vs field-s 0.960 [0.951, 0.968] (diff +0.0010, Wilson overlap TRUE)
- HR 1.00, n = 500 (null) (n 1361): rho^c mean 1.0196 [q10 0.9884, q90 1.0532], share>1 0.763; lambda-SD^c/naive 0.9751, se_field_s/naive 0.9931; field 0.930 [0.915, 0.943] vs field-s 0.940 [0.926, 0.951] (diff +0.0096, Wilson overlap TRUE)

- HR 1.75, n = 500: harm field 0.966 [0.957, 0.974] [MET]; field-s 0.941 [0.929, 0.950] [MET]; IJ two-sided H 0.961 / Hc 1.000 [MET]; joint Bonferroni 0.956 (joint_s 0.957) [MET]
- HR 1.75, n = 1000: harm field 0.944 [0.933, 0.953] [NOT MET]; field-s 0.956 [0.946, 0.964] [MET]; IJ two-sided H 0.917 / Hc 0.999 [NOT MET]; joint Bonferroni 0.941 (joint_s 0.941) [MET]
- HR 1.75, n = 1500: harm field 0.960 [0.950, 0.968] [MET]; field-s 0.960 [0.951, 0.968] [MET]; IJ two-sided H 0.913 / Hc 0.999 [NOT MET]; joint Bonferroni 0.963 (joint_s 0.964) [MET]
- HR 1.00, n = 500 (null): harm field 0.963 [0.951, 0.971] [recorded (null)]; field-s 0.940 [0.926, 0.951] [recorded (null)]; IJ two-sided H 0.986 / Hc 1.000 [recorded (null)]; joint Bonferroni 0.955 (joint_s 0.957) [recorded (null)]
