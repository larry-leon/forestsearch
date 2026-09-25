# REPORT — declcal c0 summary: the S1.8 payload (2026-09-25)

**Task:** `dev/tasks/TASK_declcal_c0_summary_2026-09-25.md`. **Branch:** `feature/glm-extension`. **Machine:** pop-os.
**Outputs:**
- script `dev/analysis/declcal_c0_summary/declcal_c0_summary.R`;
- `quarto/simulations/gbsg_020/results/declcalc0_summary_tables.csv` (401 rows);
- `quarto/simulations/gbsg_020/results/declcalc0_summary_tables.provenance.md` (definitions, cell sets, checksums);
- the typed input `~/Downloads/S18_typed_tables_2026-09-25.csv` (124 values), which is not committed.

## Gates and post-conditions

- **G0.1.** No tracked file was modified at start. These untracked files were already in the tree and were left alone:
  - `dev/analysis/`: this task's script, written in the first session;
  - `quarto/simulations/actg175/binary_020/mr_or_harm/fs_effMaxSG_mr_field_or075_n500_nb20_redes_d5000/`;
  - `quarto/simulations/actg175/binary_020/mr_or_harm/fs_effMaxSG_mr_field_or075_n500_nb20_relaunch_d5000/`;
  - `quarto/simulations/actg175/binary_020/smoke_redes.html`;
  - `quarto/simulations/actg175/binary_020/smoke_relaunch.html`;
  - `quarto/simulations/gbsg_020/results/declcalc0_summary_tables.csv`: this task's CSV, first written in the first session;
  - `quarto/simulations/gbsg_020/scripts_dinamr/logs/nullmr_findings.err`.
- **G0.2.** The branch is `feature/glm-extension`. After `git fetch origin`, the checkout is 1 ahead (the unpushed task commit `051f5c3d`) and 0 behind `origin/feature/glm-extension`.
- **G0.3 / P2.** All 23 payloads are present and tracked. The SHA-256 of each is unchanged after the run, and none shows a `git diff`. The checksums are in the provenance file.
- **The task ran in two sessions.** The first wrote the script and a 393-row CSV (Steps 1–2), then stopped because the typed-values file was missing. This session added the typed file, then ran Steps 3–6.
  - The script gained one block, which adds the comparators' `implied_pstar` rows. The typed S4 prints an implied p\* for the conventional screen, the fixed cutoff and the claim threshold, and P3 needs a computed row for each.
  - Re-running the script reproduced the 393 earlier rows unchanged and added 8.

## What was computed

- **(a)** The conventional screen's rate in B1–B6 (S3).
- **(b)** For each α ∈ {0.10 (S4), 0.05 (S5)} and each c0:
  - the calibrated rate in every B cell, the maximum over B with ties, and the rate in every C cell;
  - comparators: the conventional screen, the fixed cutoff `max_T_post >= 2.0` and the claim threshold (unshifted `declcal` rule);
  - the implied p\* rows.
- **(c)** κ̂ by n, under four candidate definitions, with the implied level `2Φ(κ̂)−1` of each.
- **(d)** The re-selection footprint.
- Cell identity, design HR and n were read from `meta`. The C cells from `meta` are C1 = HR 1.5 / n 1000, C2 = HR 1.5 / n 1500, C3 = HR 2.0 / n 1000 and C4 = HR 2.0 / n 1500. This agrees with the label mapping given for the typed columns.

**How typed rows were keyed to computed rows:**
- typed `alpha` selects S4 or S5;
- typed `c0` `c070`…`c085` → `calibrated_<c0>`;
- `ccons` → `claim_threshold`, the `declcal` payloads' unshifted rule (this is an identification, and every `ccons` value matches under it);
- blank c0 → the conventional or fixed row;
- `worst_false_declaration` → `max_B`;
- `power_HR<h>_n=<n>` → the C cell with that `meta` design and n;
- `implied_pstar_n=<n>` → S6 `level_<def>` at that n;
- S6 `kappa_n=<n>` → S6 `kappa_<def>` at that n.

**Precision.** A typed value matches if |computed − typed| ≤ half a unit in its last printed digit.

## Step 3 — comparison against the typed tables

**Every typed value has a computed row (P3). There is 1 mismatch in 124 values.** For S6 and `implied_pstar_n=` the table below uses the winning definition (item c). For the claim threshold's all-n `implied_pstar` it uses the median of per-cell medians; the pooled median gives the same printed value.

### Mismatches

| # | table | typed row | typed col | α | c0 | typed | computed | source column | cell set |
|---|---|---|---|---|---|---|---|---|---|
| 100 | S6 | log0.70 | kappa_n=1500 | 0.05 | c070 | 1.76 | 1.7545 (median of per-cell medians, all cells) | `kappa_hat_05_c070` (`declcalc0`) | B3;B6;C2;C4 |

- **No definition reproduces 1.76 here.** The four candidates at this position are:
  - pooled, all cells: 1.7657;
  - median of cell medians, all cells: 1.7545;
  - pooled, B only: 1.7805;
  - median of cell medians, B only: 1.7807.
- **The typed implied p\* in the same position (S5 `implied_pstar_n=1500`, c070, 0.921) matches the computed κ̂, not the typed one.** 2Φ(1.7545)−1 = 0.92066 prints as 0.921, while 2Φ(1.76)−1 = 0.92159 would print as 0.922.
- **The same typed value also appears in `status_curated.md` §2.11.** Its "c0 headline" reads κ̂ "2.45 / 2.07 / 1.76" for c0 0.70. That text was not edited (Step 5 permits one added line only).

### Half-unit boundary cases

Of the matches, 27 sit exactly at half a unit. These are rates of *k*/2000 whose third decimal is 5, such as 0.7905 typed as 0.791. Each is typed rounded half up. They count as matches under the stated precision and are marked "match (half-unit boundary)" in the full table.

### Full comparison, all 124 typed values

`#` is the typed file's data-row number. Computed values are rounded to 6 decimals here; the CSV has full precision.

| # | table | typed row | typed col | α | c0 | typed | computed | count/den | cell(s) | source column | result |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | S3 | 0.657 | n=500 |  |  | 0.0950 | 0.095 | 190/2000 | B1 | declared_conv | match |
| 2 | S3 | 0.657 | n=1000 |  |  | 0.0325 | 0.0325 | 65/2000 | B2 | declared_conv | match |
| 3 | S3 | 0.657 | n=1500 |  |  | 0.0070 | 0.007 | 14/2000 | B3 | declared_conv | match |
| 4 | S3 | 0.721 | n=500 |  |  | 0.2325 | 0.2325 | 465/2000 | B4 | declared_conv | match |
| 5 | S3 | 0.721 | n=1000 |  |  | 0.1110 | 0.111 | 222/2000 | B5 | declared_conv | match |
| 6 | S3 | 0.721 | n=1500 |  |  | 0.0325 | 0.0325 | 65/2000 | B6 | declared_conv | match |
| 7 | S4 | Conventional pstar=0.90 | worst_false_declaration | 0.10 |  | 0.2325 | 0.2325 | 465/2000 | B4 | declared_conv | match |
| 8 | S4 | Conventional pstar=0.90 | power_HR1.5_n=1000 | 0.10 |  | 0.699 | 0.6985 | 1397/2000 | C1 | declared_conv | match (half-unit boundary) |
| 9 | S4 | Conventional pstar=0.90 | power_HR1.5_n=1500 | 0.10 |  | 0.791 | 0.7905 | 1581/2000 | C2 | declared_conv | match (half-unit boundary) |
| 10 | S4 | Conventional pstar=0.90 | power_HR2.0_n=1000 | 0.10 |  | 0.964 | 0.9635 | 1927/2000 | C3 | declared_conv | match (half-unit boundary) |
| 11 | S4 | Conventional pstar=0.90 | power_HR2.0_n=1500 | 0.10 |  | 0.988 | 0.988 | 1976/2000 | C4 | declared_conv | match |
| 12 | S4 | Conventional pstar=0.90 | implied_pstar | 0.10 |  | 0.90 | 0.9 |  | B1;B2;B3;B4;B5;B6;C1;C2;C3;C4 | meta$p_star | match |
| 13 | S4 | Tuned fixed pstar=0.9545 | worst_false_declaration | 0.10 |  | 0.0725 | 0.0725 | 145/2000 | B4 | max_T_post>=2.0 | match |
| 14 | S4 | Tuned fixed pstar=0.9545 | power_HR1.5_n=1000 | 0.10 |  | 0.526 | 0.5255 | 1051/2000 | C1 | max_T_post>=2.0 | match (half-unit boundary) |
| 15 | S4 | Tuned fixed pstar=0.9545 | power_HR1.5_n=1500 | 0.10 |  | 0.648 | 0.6475 | 1295/2000 | C2 | max_T_post>=2.0 | match (half-unit boundary) |
| 16 | S4 | Tuned fixed pstar=0.9545 | power_HR2.0_n=1000 | 0.10 |  | 0.910 | 0.9095 | 1819/2000 | C3 | max_T_post>=2.0 | match (half-unit boundary) |
| 17 | S4 | Tuned fixed pstar=0.9545 | power_HR2.0_n=1500 | 0.10 |  | 0.971 | 0.9705 | 1941/2000 | C4 | max_T_post>=2.0 | match (half-unit boundary) |
| 18 | S4 | Tuned fixed pstar=0.9545 | implied_pstar | 0.10 |  | 0.9545 | 0.9545 |  |  | 2*pnorm(2.0)-1 | match |
| 19 | S4 | Calibrated c0=log0.70 | worst_false_declaration | 0.10 | c070 | 0.0705 | 0.0705 | 141/2000 | B4 | declared_cal10_c070 | match |
| 20 | S4 | Calibrated c0=log0.70 | power_HR1.5_n=1000 | 0.10 | c070 | 0.617 | 0.6165 | 1233/2000 | C1 | declared_cal10_c070 | match (half-unit boundary) |
| 21 | S4 | Calibrated c0=log0.70 | power_HR1.5_n=1500 | 0.10 | c070 | 0.824 | 0.8235 | 1647/2000 | C2 | declared_cal10_c070 | match (half-unit boundary) |
| 22 | S4 | Calibrated c0=log0.70 | power_HR2.0_n=1000 | 0.10 | c070 | 0.943 | 0.9425 | 1885/2000 | C3 | declared_cal10_c070 | match (half-unit boundary) |
| 23 | S4 | Calibrated c0=log0.70 | power_HR2.0_n=1500 | 0.10 | c070 | 0.992 | 0.9915 | 1983/2000 | C4 | declared_cal10_c070 | match (half-unit boundary) |
| 24 | S4 | Calibrated c0=log0.70 | implied_pstar_n=500 | 0.10 | c070 | 0.972 | 0.971629 |  | B1;B4 | 2*pnorm(kappa_hat_10_c070)-1 | match |
| 25 | S4 | Calibrated c0=log0.70 | implied_pstar_n=1000 | 0.10 | c070 | 0.932 | 0.931654 |  | B2;B5;C1;C3 | 2*pnorm(kappa_hat_10_c070)-1 | match |
| 26 | S4 | Calibrated c0=log0.70 | implied_pstar_n=1500 | 0.10 | c070 | 0.868 | 0.867861 |  | B3;B6;C2;C4 | 2*pnorm(kappa_hat_10_c070)-1 | match |
| 27 | S4 | Calibrated c0=log0.75 | worst_false_declaration | 0.10 | c075 | 0.0350 | 0.035 | 70/2000 | B4 | declared_cal10_c075 | match |
| 28 | S4 | Calibrated c0=log0.75 | power_HR1.5_n=1000 | 0.10 | c075 | 0.501 | 0.501 | 1002/2000 | C1 | declared_cal10_c075 | match |
| 29 | S4 | Calibrated c0=log0.75 | power_HR1.5_n=1500 | 0.10 | c075 | 0.723 | 0.7225 | 1445/2000 | C2 | declared_cal10_c075 | match (half-unit boundary) |
| 30 | S4 | Calibrated c0=log0.75 | power_HR2.0_n=1000 | 0.10 | c075 | 0.896 | 0.896 | 1792/2000 | C3 | declared_cal10_c075 | match |
| 31 | S4 | Calibrated c0=log0.75 | power_HR2.0_n=1500 | 0.10 | c075 | 0.982 | 0.982 | 1964/2000 | C4 | declared_cal10_c075 | match |
| 32 | S4 | Calibrated c0=log0.75 | implied_pstar_n=500 | 0.10 | c075 | 0.983 | 0.983326 |  | B1;B4 | 2*pnorm(kappa_hat_10_c075)-1 | match |
| 33 | S4 | Calibrated c0=log0.75 | implied_pstar_n=1000 | 0.10 | c075 | 0.963 | 0.962982 |  | B2;B5;C1;C3 | 2*pnorm(kappa_hat_10_c075)-1 | match |
| 34 | S4 | Calibrated c0=log0.75 | implied_pstar_n=1500 | 0.10 | c075 | 0.931 | 0.930915 |  | B3;B6;C2;C4 | 2*pnorm(kappa_hat_10_c075)-1 | match |
| 35 | S4 | Calibrated c0=log0.80 | worst_false_declaration | 0.10 | c080 | 0.0165 | 0.0165 | 33/2000 | B4 | declared_cal10_c080 | match |
| 36 | S4 | Calibrated c0=log0.80 | power_HR1.5_n=1000 | 0.10 | c080 | 0.370 | 0.3695 | 739/2000 | C1 | declared_cal10_c080 | match (half-unit boundary) |
| 37 | S4 | Calibrated c0=log0.80 | power_HR1.5_n=1500 | 0.10 | c080 | 0.599 | 0.5985 | 1197/2000 | C2 | declared_cal10_c080 | match (half-unit boundary) |
| 38 | S4 | Calibrated c0=log0.80 | power_HR2.0_n=1000 | 0.10 | c080 | 0.833 | 0.833 | 1666/2000 | C3 | declared_cal10_c080 | match |
| 39 | S4 | Calibrated c0=log0.80 | power_HR2.0_n=1500 | 0.10 | c080 | 0.964 | 0.9635 | 1927/2000 | C4 | declared_cal10_c080 | match (half-unit boundary) |
| 40 | S4 | Calibrated c0=log0.80 | implied_pstar_n=500 | 0.10 | c080 | 0.990 | 0.990402 |  | B1;B4 | 2*pnorm(kappa_hat_10_c080)-1 | match |
| 41 | S4 | Calibrated c0=log0.80 | implied_pstar_n=1000 | 0.10 | c080 | 0.981 | 0.980807 |  | B2;B5;C1;C3 | 2*pnorm(kappa_hat_10_c080)-1 | match |
| 42 | S4 | Calibrated c0=log0.80 | implied_pstar_n=1500 | 0.10 | c080 | 0.966 | 0.966227 |  | B3;B6;C2;C4 | 2*pnorm(kappa_hat_10_c080)-1 | match |
| 43 | S4 | Calibrated c0=log0.85 | worst_false_declaration | 0.10 | c085 | 0.0055 | 0.0055 | 11/2000 | B4 | declared_cal10_c085 | match |
| 44 | S4 | Calibrated c0=log0.85 | power_HR1.5_n=1000 | 0.10 | c085 | 0.265 | 0.265 | 530/2000 | C1 | declared_cal10_c085 | match |
| 45 | S4 | Calibrated c0=log0.85 | power_HR1.5_n=1500 | 0.10 | c085 | 0.459 | 0.4585 | 917/2000 | C2 | declared_cal10_c085 | match (half-unit boundary) |
| 46 | S4 | Calibrated c0=log0.85 | power_HR2.0_n=1000 | 0.10 | c085 | 0.756 | 0.7555 | 1511/2000 | C3 | declared_cal10_c085 | match (half-unit boundary) |
| 47 | S4 | Calibrated c0=log0.85 | power_HR2.0_n=1500 | 0.10 | c085 | 0.921 | 0.921 | 1842/2000 | C4 | declared_cal10_c085 | match |
| 48 | S4 | Calibrated c0=log0.85 | implied_pstar_n=500 | 0.10 | c085 | 0.995 | 0.994605 |  | B1;B4 | 2*pnorm(kappa_hat_10_c085)-1 | match |
| 49 | S4 | Calibrated c0=log0.85 | implied_pstar_n=1000 | 0.10 | c085 | 0.990 | 0.990453 |  | B2;B5;C1;C3 | 2*pnorm(kappa_hat_10_c085)-1 | match |
| 50 | S4 | Calibrated c0=log0.85 | implied_pstar_n=1500 | 0.10 | c085 | 0.985 | 0.984653 |  | B3;B6;C2;C4 | 2*pnorm(kappa_hat_10_c085)-1 | match |
| 51 | S4 | Calibrated c0=ccons | worst_false_declaration | 0.10 | ccons | 0.0005 | 0.0005 | 1/2000 | B1;B4 | declared_cal10 | match |
| 52 | S4 | Calibrated c0=ccons | power_HR1.5_n=1000 | 0.10 | ccons | 0.083 | 0.0825 | 165/2000 | C1 | declared_cal10 | match (half-unit boundary) |
| 53 | S4 | Calibrated c0=ccons | power_HR1.5_n=1500 | 0.10 | ccons | 0.145 | 0.1445 | 289/2000 | C2 | declared_cal10 | match (half-unit boundary) |
| 54 | S4 | Calibrated c0=ccons | power_HR2.0_n=1000 | 0.10 | ccons | 0.448 | 0.4475 | 895/2000 | C3 | declared_cal10 | match (half-unit boundary) |
| 55 | S4 | Calibrated c0=ccons | power_HR2.0_n=1500 | 0.10 | ccons | 0.726 | 0.7255 | 1451/2000 | C4 | declared_cal10 | match (half-unit boundary) |
| 56 | S4 | Calibrated c0=ccons | implied_pstar | 0.10 | ccons | 0.9992 | 0.99916 |  | B1;B2;B3;B4;B5;B6;C1;C2;C3;C4 | 2*pnorm(kappa_hat_10)-1 | match |
| 57 | S5 | Calibrated c0=log0.70 | worst_false_declaration | 0.05 | c070 | 0.0315 | 0.0315 | 63/2000 | B4 | declared_cal05_c070 | match |
| 58 | S5 | Calibrated c0=log0.70 | power_HR1.5_n=1000 | 0.05 | c070 | 0.503 | 0.5025 | 1005/2000 | C1 | declared_cal05_c070 | match (half-unit boundary) |
| 59 | S5 | Calibrated c0=log0.70 | power_HR1.5_n=1500 | 0.05 | c070 | 0.747 | 0.747 | 1494/2000 | C2 | declared_cal05_c070 | match |
| 60 | S5 | Calibrated c0=log0.70 | power_HR2.0_n=1000 | 0.05 | c070 | 0.896 | 0.896 | 1792/2000 | C3 | declared_cal05_c070 | match |
| 61 | S5 | Calibrated c0=log0.70 | power_HR2.0_n=1500 | 0.05 | c070 | 0.983 | 0.983 | 1966/2000 | C4 | declared_cal05_c070 | match |
| 62 | S5 | Calibrated c0=log0.70 | implied_pstar_n=500 | 0.05 | c070 | 0.986 | 0.985756 |  | B1;B4 | 2*pnorm(kappa_hat_05_c070)-1 | match |
| 63 | S5 | Calibrated c0=log0.70 | implied_pstar_n=1000 | 0.05 | c070 | 0.962 | 0.961708 |  | B2;B5;C1;C3 | 2*pnorm(kappa_hat_05_c070)-1 | match |
| 64 | S5 | Calibrated c0=log0.70 | implied_pstar_n=1500 | 0.05 | c070 | 0.921 | 0.920659 |  | B3;B6;C2;C4 | 2*pnorm(kappa_hat_05_c070)-1 | match |
| 65 | S5 | Calibrated c0=log0.75 | worst_false_declaration | 0.05 | c075 | 0.0130 | 0.013 | 26/2000 | B4 | declared_cal05_c075 | match |
| 66 | S5 | Calibrated c0=log0.75 | power_HR1.5_n=1000 | 0.05 | c075 | 0.373 | 0.373 | 746/2000 | C1 | declared_cal05_c075 | match |
| 67 | S5 | Calibrated c0=log0.75 | power_HR1.5_n=1500 | 0.05 | c075 | 0.620 | 0.62 | 1240/2000 | C2 | declared_cal05_c075 | match |
| 68 | S5 | Calibrated c0=log0.75 | power_HR2.0_n=1000 | 0.05 | c075 | 0.834 | 0.834 | 1668/2000 | C3 | declared_cal05_c075 | match |
| 69 | S5 | Calibrated c0=log0.75 | power_HR2.0_n=1500 | 0.05 | c075 | 0.966 | 0.966 | 1932/2000 | C4 | declared_cal05_c075 | match |
| 70 | S5 | Calibrated c0=log0.75 | implied_pstar_n=500 | 0.05 | c075 | 0.992 | 0.992018 |  | B1;B4 | 2*pnorm(kappa_hat_05_c075)-1 | match |
| 71 | S5 | Calibrated c0=log0.75 | implied_pstar_n=1000 | 0.05 | c075 | 0.980 | 0.980378 |  | B2;B5;C1;C3 | 2*pnorm(kappa_hat_05_c075)-1 | match |
| 72 | S5 | Calibrated c0=log0.75 | implied_pstar_n=1500 | 0.05 | c075 | 0.961 | 0.96112 |  | B3;B6;C2;C4 | 2*pnorm(kappa_hat_05_c075)-1 | match |
| 73 | S5 | Calibrated c0=log0.80 | worst_false_declaration | 0.05 | c080 | 0.0040 | 0.004 | 8/2000 | B4;B5 | declared_cal05_c080 | match |
| 74 | S5 | Calibrated c0=log0.80 | power_HR1.5_n=1000 | 0.05 | c080 | 0.271 | 0.2705 | 541/2000 | C1 | declared_cal05_c080 | match (half-unit boundary) |
| 75 | S5 | Calibrated c0=log0.80 | power_HR1.5_n=1500 | 0.05 | c080 | 0.491 | 0.491 | 982/2000 | C2 | declared_cal05_c080 | match |
| 76 | S5 | Calibrated c0=log0.80 | power_HR2.0_n=1000 | 0.05 | c080 | 0.753 | 0.753 | 1506/2000 | C3 | declared_cal05_c080 | match |
| 77 | S5 | Calibrated c0=log0.80 | power_HR2.0_n=1500 | 0.05 | c080 | 0.928 | 0.9275 | 1855/2000 | C4 | declared_cal05_c080 | match (half-unit boundary) |
| 78 | S5 | Calibrated c0=log0.80 | implied_pstar_n=500 | 0.05 | c080 | 0.996 | 0.995607 |  | B1;B4 | 2*pnorm(kappa_hat_05_c080)-1 | match |
| 79 | S5 | Calibrated c0=log0.80 | implied_pstar_n=1000 | 0.05 | c080 | 0.990 | 0.990375 |  | B2;B5;C1;C3 | 2*pnorm(kappa_hat_05_c080)-1 | match |
| 80 | S5 | Calibrated c0=log0.80 | implied_pstar_n=1500 | 0.05 | c080 | 0.982 | 0.98208 |  | B3;B6;C2;C4 | 2*pnorm(kappa_hat_05_c080)-1 | match |
| 81 | S5 | Calibrated c0=log0.85 | worst_false_declaration | 0.05 | c085 | 0.0010 | 0.001 | 2/2000 | B4;B5;B6 | declared_cal05_c085 | match |
| 82 | S5 | Calibrated c0=log0.85 | power_HR1.5_n=1000 | 0.05 | c085 | 0.188 | 0.188 | 376/2000 | C1 | declared_cal05_c085 | match |
| 83 | S5 | Calibrated c0=log0.85 | power_HR1.5_n=1500 | 0.05 | c085 | 0.363 | 0.3625 | 725/2000 | C2 | declared_cal05_c085 | match (half-unit boundary) |
| 84 | S5 | Calibrated c0=log0.85 | power_HR2.0_n=1000 | 0.05 | c085 | 0.676 | 0.676 | 1352/2000 | C3 | declared_cal05_c085 | match |
| 85 | S5 | Calibrated c0=log0.85 | power_HR2.0_n=1500 | 0.05 | c085 | 0.882 | 0.882 | 1764/2000 | C4 | declared_cal05_c085 | match |
| 86 | S5 | Calibrated c0=log0.85 | implied_pstar_n=500 | 0.05 | c085 | 0.998 | 0.997638 |  | B1;B4 | 2*pnorm(kappa_hat_05_c085)-1 | match |
| 87 | S5 | Calibrated c0=log0.85 | implied_pstar_n=1000 | 0.05 | c085 | 0.995 | 0.995465 |  | B2;B5;C1;C3 | 2*pnorm(kappa_hat_05_c085)-1 | match |
| 88 | S5 | Calibrated c0=log0.85 | implied_pstar_n=1500 | 0.05 | c085 | 0.992 | 0.992353 |  | B3;B6;C2;C4 | 2*pnorm(kappa_hat_05_c085)-1 | match |
| 89 | S5 | Calibrated c0=ccons | worst_false_declaration | 0.05 | ccons | 0.0005 | 0.0005 | 1/2000 | B1;B4 | declared_cal05 | match |
| 90 | S5 | Calibrated c0=ccons | power_HR1.5_n=1000 | 0.05 | ccons | 0.050 | 0.0495 | 99/2000 | C1 | declared_cal05 | match (half-unit boundary) |
| 91 | S5 | Calibrated c0=ccons | power_HR1.5_n=1500 | 0.05 | ccons | 0.094 | 0.0935 | 187/2000 | C2 | declared_cal05 | match (half-unit boundary) |
| 92 | S5 | Calibrated c0=ccons | power_HR2.0_n=1000 | 0.05 | ccons | 0.355 | 0.355 | 710/2000 | C3 | declared_cal05 | match |
| 93 | S5 | Calibrated c0=ccons | power_HR2.0_n=1500 | 0.05 | ccons | 0.641 | 0.641 | 1282/2000 | C4 | declared_cal05 | match |
| 94 | S5 | Calibrated c0=ccons | implied_pstar | 0.05 | ccons | 0.9997 | 0.99967 |  | B1;B2;B3;B4;B5;B6;C1;C2;C3;C4 | 2*pnorm(kappa_hat_05)-1 | match |
| 95 | S6 | log0.70 | kappa_n=500 | 0.10 | c070 | 2.19 | 2.19212 |  | B1;B4 | kappa_hat_10_c070 | match |
| 96 | S6 | log0.70 | kappa_n=1000 | 0.10 | c070 | 1.82 | 1.822721 |  | B2;B5;C1;C3 | kappa_hat_10_c070 | match |
| 97 | S6 | log0.70 | kappa_n=1500 | 0.10 | c070 | 1.51 | 1.505719 |  | B3;B6;C2;C4 | kappa_hat_10_c070 | match |
| 98 | S6 | log0.70 | kappa_n=500 | 0.05 | c070 | 2.45 | 2.451045 |  | B1;B4 | kappa_hat_05_c070 | match |
| 99 | S6 | log0.70 | kappa_n=1000 | 0.05 | c070 | 2.07 | 2.071718 |  | B2;B5;C1;C3 | kappa_hat_05_c070 | match |
| 100 | S6 | log0.70 | kappa_n=1500 | 0.05 | c070 | 1.76 | 1.754522 |  | B3;B6;C2;C4 | kappa_hat_05_c070 | **MISMATCH** |
| 101 | S6 | log0.75 | kappa_n=500 | 0.10 | c075 | 2.39 | 2.39381 |  | B1;B4 | kappa_hat_10_c075 | match |
| 102 | S6 | log0.75 | kappa_n=1000 | 0.10 | c075 | 2.09 | 2.08556 |  | B2;B5;C1;C3 | kappa_hat_10_c075 | match |
| 103 | S6 | log0.75 | kappa_n=1500 | 0.10 | c075 | 1.82 | 1.817863 |  | B3;B6;C2;C4 | kappa_hat_10_c075 | match |
| 104 | S6 | log0.75 | kappa_n=500 | 0.05 | c075 | 2.65 | 2.652838 |  | B1;B4 | kappa_hat_05_c075 | match |
| 105 | S6 | log0.75 | kappa_n=1000 | 0.05 | c075 | 2.33 | 2.33349 |  | B2;B5;C1;C3 | kappa_hat_05_c075 | match |
| 106 | S6 | log0.75 | kappa_n=1500 | 0.05 | c075 | 2.07 | 2.065451 |  | B3;B6;C2;C4 | kappa_hat_05_c075 | match |
| 107 | S6 | log0.80 | kappa_n=500 | 0.10 | c080 | 2.59 | 2.589987 |  | B1;B4 | kappa_hat_10_c080 | match |
| 108 | S6 | log0.80 | kappa_n=1000 | 0.10 | c080 | 2.34 | 2.341758 |  | B2;B5;C1;C3 | kappa_hat_10_c080 | match |
| 109 | S6 | log0.80 | kappa_n=1500 | 0.10 | c080 | 2.12 | 2.122767 |  | B3;B6;C2;C4 | kappa_hat_10_c080 | match |
| 110 | S6 | log0.80 | kappa_n=500 | 0.05 | c080 | 2.85 | 2.848463 |  | B1;B4 | kappa_hat_05_c080 | match |
| 111 | S6 | log0.80 | kappa_n=1000 | 0.05 | c080 | 2.59 | 2.589029 |  | B2;B5;C1;C3 | kappa_hat_05_c080 | match |
| 112 | S6 | log0.80 | kappa_n=1500 | 0.05 | c080 | 2.37 | 2.36726 |  | B3;B6;C2;C4 | kappa_hat_05_c080 | match |
| 113 | S6 | log0.85 | kappa_n=500 | 0.10 | c085 | 2.78 | 2.782476 |  | B1;B4 | kappa_hat_10_c085 | match |
| 114 | S6 | log0.85 | kappa_n=1000 | 0.10 | c085 | 2.59 | 2.591811 |  | B2;B5;C1;C3 | kappa_hat_10_c085 | match |
| 115 | S6 | log0.85 | kappa_n=1500 | 0.10 | c085 | 2.42 | 2.424082 |  | B3;B6;C2;C4 | kappa_hat_10_c085 | match |
| 116 | S6 | log0.85 | kappa_n=500 | 0.05 | c085 | 3.04 | 3.040473 |  | B1;B4 | kappa_hat_05_c085 | match |
| 117 | S6 | log0.85 | kappa_n=1000 | 0.05 | c085 | 2.84 | 2.838356 |  | B2;B5;C1;C3 | kappa_hat_05_c085 | match |
| 118 | S6 | log0.85 | kappa_n=1500 | 0.05 | c085 | 2.67 | 2.667259 |  | B3;B6;C2;C4 | kappa_hat_05_c085 | match |
| 119 | S6 | ccons | kappa_n=500 | 0.10 | ccons | 3.34 | 3.337849 |  | B1;B4 | kappa_hat_10 | match |
| 120 | S6 | ccons | kappa_n=1000 | 0.10 | ccons | 3.34 | 3.342749 |  | B2;B5;C1;C3 | kappa_hat_10 | match |
| 121 | S6 | ccons | kappa_n=1500 | 0.10 | ccons | 3.33 | 3.331844 |  | B3;B6;C2;C4 | kappa_hat_10 | match |
| 122 | S6 | ccons | kappa_n=500 | 0.05 | ccons | 3.60 | 3.597673 |  | B1;B4 | kappa_hat_05 | match |
| 123 | S6 | ccons | kappa_n=1000 | 0.05 | ccons | 3.59 | 3.591427 |  | B2;B5;C1;C3 | kappa_hat_05 | match |
| 124 | S6 | ccons | kappa_n=1500 | 0.05 | ccons | 3.58 | 3.58012 |  | B3;B6;C2;C4 | kappa_hat_05 | match |

## Item (c) — the definition of the S6 cutoff

**Winner: `median_of_cell_medians_allcells`.** This is the median of the per-cell medians of `kappa_hat`, over every cell (B and C) at that n. Recorded as the definition.

- **κ̂ values reproduced, of 30 typed:**
  - median of cell medians, all cells: **29**, the one miss being #100 above;
  - pooled median, all cells: 18;
  - pooled median, B only: 10;
  - median of cell medians, B only: 10.
- **Implied-p\* values reproduced (S4/S5 `implied_pstar_n=`), of 24 typed:**
  - median of cell medians, all cells: **24**;
  - pooled median, all cells: 10;
  - median of cell medians, B only: 9;
  - pooled median, B only: 8.
- **The manuscript's wording** is "medians across cells at each sample size". That is the median of per-cell medians, and the cell set is every cell at that n, including the power cells.
- **At n 500 the four candidates reduce to two**, because the all-cells set is B1 and B4 only. With two cells, the median of cell medians is their mean.
- **The claim threshold's single implied p\*** (0.9992 at α 0.10, 0.9997 at α 0.05) is reproduced by both all-n summaries (0.99916 / 0.99916; 0.99967 / 0.99967).

κ̂ under all four definitions, by typed position:

| # | typed row | typed col | α | typed | pooled, all cells | median of cell medians, all cells | pooled, B only | median of cell medians, B only |
|---|---|---|---|---|---|---|---|---|
|  95 | log0.70 | kappa_n=500 | 0.10 | 2.19 | 2.1902 | 2.1921 | 2.1902 | 2.1921 |
|  96 | log0.70 | kappa_n=1000 | 0.10 | 1.82 | 1.8306 | 1.8227 | 1.8462 | 1.8471 |
|  97 | log0.70 | kappa_n=1500 | 0.10 | 1.51 | 1.5157 | 1.5057 | 1.5317 | 1.5325 |
|  98 | log0.70 | kappa_n=500 | 0.05 | 2.45 | 2.4505 | 2.4510 | 2.4505 | 2.4510 |
|  99 | log0.70 | kappa_n=1000 | 0.05 | 2.07 | 2.0813 | 2.0717 | 2.0978 | 2.0984 |
| 100 | log0.70 | kappa_n=1500 | 0.05 | 1.76 | 1.7657 | 1.7545 | 1.7805 | 1.7807 |
| 101 | log0.75 | kappa_n=500 | 0.10 | 2.39 | 2.3942 | 2.3938 | 2.3942 | 2.3938 |
| 102 | log0.75 | kappa_n=1000 | 0.10 | 2.09 | 2.0931 | 2.0856 | 2.1068 | 2.1085 |
| 103 | log0.75 | kappa_n=1500 | 0.10 | 1.82 | 1.8276 | 1.8179 | 1.8411 | 1.8416 |
| 104 | log0.75 | kappa_n=500 | 0.05 | 2.65 | 2.6535 | 2.6528 | 2.6535 | 2.6528 |
| 105 | log0.75 | kappa_n=1000 | 0.05 | 2.33 | 2.3423 | 2.3335 | 2.3574 | 2.3571 |
| 106 | log0.75 | kappa_n=1500 | 0.05 | 2.07 | 2.0741 | 2.0655 | 2.0892 | 2.0892 |
| 107 | log0.80 | kappa_n=500 | 0.10 | 2.59 | 2.5901 | 2.5900 | 2.5901 | 2.5900 |
| 108 | log0.80 | kappa_n=1000 | 0.10 | 2.34 | 2.3480 | 2.3418 | 2.3607 | 2.3606 |
| 109 | log0.80 | kappa_n=1500 | 0.10 | 2.12 | 2.1305 | 2.1228 | 2.1435 | 2.1434 |
| 110 | log0.80 | kappa_n=500 | 0.05 | 2.85 | 2.8497 | 2.8485 | 2.8497 | 2.8485 |
| 111 | log0.80 | kappa_n=1000 | 0.05 | 2.59 | 2.5951 | 2.5890 | 2.6101 | 2.6094 |
| 112 | log0.80 | kappa_n=1500 | 0.05 | 2.37 | 2.3756 | 2.3673 | 2.3898 | 2.3885 |
| 113 | log0.85 | kappa_n=500 | 0.10 | 2.78 | 2.7825 | 2.7825 | 2.7825 | 2.7825 |
| 114 | log0.85 | kappa_n=1000 | 0.10 | 2.59 | 2.5981 | 2.5918 | 2.6090 | 2.6081 |
| 115 | log0.85 | kappa_n=1500 | 0.10 | 2.42 | 2.4296 | 2.4241 | 2.4417 | 2.4410 |
| 116 | log0.85 | kappa_n=500 | 0.05 | 3.04 | 3.0397 | 3.0405 | 3.0397 | 3.0405 |
| 117 | log0.85 | kappa_n=1000 | 0.05 | 2.84 | 2.8450 | 2.8384 | 2.8567 | 2.8556 |
| 118 | log0.85 | kappa_n=1500 | 0.05 | 2.67 | 2.6738 | 2.6673 | 2.6853 | 2.6852 |
| 119 | ccons | kappa_n=500 | 0.10 | 3.34 | 3.3374 | 3.3378 | 3.3374 | 3.3378 |
| 120 | ccons | kappa_n=1000 | 0.10 | 3.34 | 3.3448 | 3.3427 | 3.3478 | 3.3480 |
| 121 | ccons | kappa_n=1500 | 0.10 | 3.33 | 3.3333 | 3.3318 | 3.3368 | 3.3364 |
| 122 | ccons | kappa_n=500 | 0.05 | 3.60 | 3.5975 | 3.5977 | 3.5975 | 3.5977 |
| 123 | ccons | kappa_n=1000 | 0.05 | 3.59 | 3.5943 | 3.5914 | 3.5982 | 3.5976 |
| 124 | ccons | kappa_n=1500 | 0.05 | 3.58 | 3.5820 | 3.5801 | 3.5865 | 3.5860 |

## Item (d) — the re-selection footprint

This is over the C cells at α 0.05: among replicates with `declared_cal05_<c0> == 1`, the share with `n_admitted_cal05_<c0> > 1`.

| c0 | share | count / declared |
|---|---|---|
| 0.70 | 0.9682 | 6058 / 6257 |
| 0.75 | 0.9656 | 5394 / 5586 |
| 0.80 | 0.9652 | 4714 / 4884 |
| 0.85 | 0.9666 | 4076 / 4217 |
| pooled | **0.9665** | 20242 / 20944 |

The pooled figure of 96.6% is the "~97%" of `status_curated.md` §2.11. The column exists at α 0.05 only.

## Ties in the B maximum (reported, not resolved)

- **Most rows have one worst cell, B4 (design 0.721, n 500).** This covers every calibrated row at α 0.10, c070 and c075 at α 0.05, the conventional screen and the fixed cutoff.
- **Ties:**
  - S5 c080: B4 and B5 (0.0040);
  - S5 c085: B4, B5 and B6 (0.0010);
  - the claim threshold at both α: B1 and B4 (0.0005).

## Not established from source

- **The fixed cutoff's implied p\* (0.9545)** is `2Φ(2.0)−1` of the task-specified constant `k_fix = 2.0`, not a payload column. It matches the typed value. No payload records a "tuned" p\* or how 2.0 was tuned: `not established from source`.
- **Item (e)**, the plug-in and free-check figures, is out of scope by the task and was not computed.
