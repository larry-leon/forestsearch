# REPORT — approximate c0 table from 120 field captures

Task: `dev/tasks/TASK_declcal_c0_approx_2026-09-22.md`. Branch `feature/glm-extension`, pop-os. No `R/` change. No push.

**These are plug-in fixed-cutoff rates.** Each rate is `mean(max_T_pre >= median κ̂(c0))`, where the median is taken over the captures at that n and applied to every committed replicate at that n. They are not the per-replicate calibrated rule, which uses each replicate's own κ̂(c0). The campaign task (`TASK_declcal_c0_campaign_2026-09-22.md`) supplies the exact version.

## Pins and commits

- Step 0: **`80abee8d`**, `docs(tasks): add c0 approximate-table task (2026-09-22)`.
- Installed `forestsearch` 0.3.5.9000, built 2026-09-23 02:32 UTC. `declaration_c0` is a formal of the installed `forestsearch:::fs_mr_inference` (internal, not exported) with default `NULL`, so no reinstall was needed.
- Comparators: the committed `results/declcal_{inull,power}_*_res_1_2000.rds` at `81752681` (B_cal 500). They are read-only, and the script asserts they are unchanged since `81752681`.
- Results: the commit that adds this file.

## Replicate count

- Section 1 of the task lists "replicates 1–20 of each of B1, B4, B2, B5, B3, B6". That is 6 × 20 = **120** captures, and all 120 were run.
- The "60 in all" in the task was an arithmetic error. Larry confirmed that the 1–20 range is the instruction.
- The committed table therefore uses **40 captures per n** (20 from each of the two B cells at that n).
- The 60-capture reading (reps 1–10 per cell) is kept as a sensitivity note below. It needed no extra compute: each replicate is seeded by `8316951 + rep`, so reps 1–10 are exactly what a 60-search run would have produced.

## Captures

- **Scripts:** `scripts_dinamr/declcal_c0approx_run.R` is a transplant of `declcal_run.R` at `80abee8d`. The lines changed are listed in its header, using `declcal_run.R` line numbers:
  - :57, the replicate-range default (20);
  - :73, the campaign tag;
  - :326, `declaration_c0 = c(0.70, 0.75, 0.80, 0.85)` added to the `fs_mr_inference()` call;
  - after :213, the per-c0 columns added to the record;
  - after :354, those columns filled from `fld$Mstar_c0` with the same type-1 quantile rule as the unshifted κ̂;
  - :419, the meta.
- B = 2000 comes from `DECLCAL_BCAL`, the engine's own default.
- **Driver:** `scripts_dinamr/declcal_c0approx.sh` uses the `declcal.sh` per-cell environment:
  - `FS_S7_C1 = FS_S7_C2 = 1.0`, `effMaxSG`, neighbourhood 0.20, template floors, centred Poisson multipliers;
  - the cells from `declcal_inull.cells`;
  - `FS_S7_NSIMS=20`, `FS_S7_WORKERS=10`, with the six cells run concurrently (60 workers).
- **Run:** 86 s wall. 120/120 status `ok`. In all six cells the structural-null design-point check passed and no replicate disagreed with the search's own declaration.
- **Post-processing:** `scripts_dinamr/declcal_c0approx_findings.R` writes `results/declcal_c0approx_res.rds` (120 rows × 63 columns, plus aux, the approx vectors and per-cell meta) and `scripts_dinamr/logs/declcal_c0approx.txt`, which holds the full tables with Wilson intervals, per-cell medians and α 0.10.

## Identity gate — PASS

- `max_T_pre` and `G_pre` are `identical()` to the committed declcal payloads on all 120 replicates.
- Difference in unshifted κ̂₀.₀₅ (B 2000 − B 500): mean +0.0083, SD 0.0572, max |diff| 0.1722.
- 119 of 120 are within 0.15. The one exception is B6 rep 20 at 0.172. The task gave 0.15 as an expectation, not a stop condition, and 0.172 is about 3 SD of the B 2000 vs B 500 Monte-Carlo difference.

## Table 1 — median κ̂(c0) per n (40 captures, B = 2000)

| n | c0 | median κ̂₀.₀₅ (min–max) | IQR | implied p\* | median κ̂₀.₁₀ | mean fw_1645 |
|---|---|---|---|---|---|---|
| 500 | 0.70 | 2.448 (2.284–2.687) | 2.419–2.521 | 0.9856 | 2.200 | 0.323 |
| 500 | 0.75 | 2.657 (2.482–2.865) | 2.621–2.703 | 0.9921 | 2.402 | 0.447 |
| 500 | 0.80 | 2.848 (2.703–3.030) | 2.819–2.897 | 0.9956 | 2.601 | 0.578 |
| 500 | 0.85 | 3.046 (2.915–3.204) | 3.015–3.087 | 0.9977 | 2.793 | 0.696 |
| 1000 | 0.70 | 2.080 (1.942–2.240) | 2.048–2.120 | 0.9625 | 1.841 | 0.163 |
| 1000 | 0.75 | 2.343 (2.205–2.477) | 2.313–2.370 | 0.9809 | 2.095 | 0.284 |
| 1000 | 0.80 | 2.586 (2.476–2.733) | 2.563–2.616 | 0.9903 | 2.348 | 0.439 |
| 1000 | 0.85 | 2.835 (2.748–2.974) | 2.817–2.868 | 0.9954 | 2.595 | 0.608 |
| 1500 | 0.70 | 1.758 (1.624–1.956) | 1.729–1.817 | 0.9212 | 1.514 | 0.072 |
| 1500 | 0.75 | 2.077 (1.954–2.246) | 2.035–2.117 | 0.9622 | 1.827 | 0.160 |
| 1500 | 0.80 | 2.391 (2.261–2.513) | 2.352–2.429 | 0.9832 | 2.129 | 0.306 |
| 1500 | 0.85 | 2.690 (2.570–2.796) | 2.651–2.720 | 0.9929 | 2.434 | 0.502 |

- **Direction:** κ̂(c0) rises with c0 at every n, and falls with n at fixed c0. As n grows, the shift (c2 − c0)/σ_D(g) grows, because σ_D shrinks.
- **The approximation's error bar:** the min–max spread is about ±0.15 around the median.
- **The unshifted reference (c0 = c2):** the median κ̂₀.₀₅ is about 3.59 at every n (see the log).
- **fw_1645:** the family-wise size of the p\* = 0.90 rule if the true null sat at c0.

## Table 2 — declaration rate at α 0.05 (x/2000; Wilson 95% in the log)

| cell | n | p\* 0.90 as executed | c0 0.70 | c0 0.75 | c0 0.80 | c0 0.85 | c0 = c2 (committed) |
|---|---|---|---|---|---|---|---|
| B1 | 500 | 0.0950 | 0.0065 | 0.0025 | 0.0015 | 0.0005 | 0.0005 |
| B2 | 1000 | 0.0325 | 0.0065 | 0.0030 | 0.0010 | 0.0005 | 0.0000 |
| B3 | 1500 | 0.0070 | 0.0040 | 0.0020 | 0.0015 | 0.0000 | 0.0000 |
| B4 | 500 | 0.2325 | 0.0270 | 0.0080 | 0.0030 | 0.0010 | 0.0005 |
| B5 | 1000 | 0.1110 | 0.0255 | 0.0100 | 0.0035 | 0.0005 | 0.0000 |
| B6 | 1500 | 0.0325 | 0.0185 | 0.0040 | 0.0020 | 0.0010 | 0.0000 |
| C1 | 1000 | 0.6985 | 0.5025 | 0.3650 | 0.2685 | 0.1925 | 0.0495 |
| C2 | 1500 | 0.7905 | 0.7505 | 0.6190 | 0.4785 | 0.3515 | 0.0935 |
| C3 | 1000 | 0.9635 | 0.8945 | 0.8325 | 0.7570 | 0.6740 | 0.3550 |
| C4 | 1500 | 0.9880 | 0.9850 | 0.9650 | 0.9295 | 0.8780 | 0.6410 |

- **"p\* 0.90 as executed"** is `declared_conv`: the rounded rule on the post-reduction family.
- **"c0 = c2"** is `declared_cal05`: the per-replicate κ̂ at B = 500.
- **The c0 columns** apply the Table 1 median to `max_T_pre`, the pre-reduction family.
- **α 0.10:** the table is in the log. There the worst B rate is 0.057 at c0 0.70 (B4), 0.030 at 0.75, 0.012 at 0.80 and 0.003 at 0.85.

## Table 3 — worst uniform-benefit rate vs power

**α 0.05**

| rule | worst B (cell) | HR 1.5 n 1000 | HR 1.5 n 1500 | HR 2.0 n 1000 | HR 2.0 n 1500 |
|---|---|---|---|---|---|
| fixed p\* 0.9545 (k 2.0, committed, max_T_pre) | 0.0955 (B4) | 0.5305 | 0.6500 | 0.9120 | 0.9710 |
| c0 0.70 | 0.0270 (B4) | 0.5025 | 0.7505 | 0.8945 | 0.9850 |
| c0 0.75 | 0.0100 (B5) | 0.3650 | 0.6190 | 0.8325 | 0.9650 |
| c0 0.80 | 0.0035 (B5) | 0.2685 | 0.4785 | 0.7570 | 0.9295 |
| c0 0.85 | 0.0010 (B4) | 0.1925 | 0.3515 | 0.6740 | 0.8780 |
| c0 = c2 (committed, per-replicate) | 0.0005 (B1) | 0.0495 | 0.0935 | 0.3550 | 0.6410 |

**α 0.10**

| rule | worst B (cell) | HR 1.5 n 1000 | HR 1.5 n 1500 | HR 2.0 n 1000 | HR 2.0 n 1500 |
|---|---|---|---|---|---|
| fixed p\* 0.9545 (k 2.0, committed, max_T_pre) | 0.0955 (B4) | 0.5305 | 0.6500 | 0.9120 | 0.9710 |
| c0 0.70 | 0.0570 (B4) | 0.6055 | 0.8225 | 0.9380 | 0.9915 |
| c0 0.75 | 0.0300 (B4) | 0.4980 | 0.7185 | 0.8915 | 0.9810 |
| c0 0.80 | 0.0120 (B4) | 0.3635 | 0.6005 | 0.8310 | 0.9600 |
| c0 0.85 | 0.0030 (B4) | 0.2655 | 0.4615 | 0.7560 | 0.9210 |
| c0 = c2 (committed, per-replicate) | 0.0005 (B1) | 0.0825 | 0.1445 | 0.4475 | 0.7255 |

- The k 2.0 row reproduces the `max_T_pre` column of the committed `logs/declcal_fixedk_practical.txt` exactly. On the post-reduction family the committed values are 0.0725 (B4), 0.5255, 0.6475, 0.9095 and 0.9705.
- **Compared with fixed k 2.0 at α 0.05:**
  - c0 0.70 holds the worst uniform-benefit rate to 0.027, against 0.0955.
  - It gives up about 0.03 of power at n 1000 (C1 0.50 vs 0.53; C3 0.89 vs 0.91).
  - It gains at n 1500 (C2 0.75 vs 0.65), because κ̂(c0) falls with n while k 2.0 is fixed.
- **Rising c0** trades power steadily for size. The largest effect is at HR 1.5, n 1000, where power falls from 0.50 to 0.19 between c0 0.70 and 0.85.
- **Every c0 in the grid** recovers most of the power that the c0 = c2 calibration gives up.

## Sensitivity — the 60-capture reading (reps 1–10 per cell, 20 per n)

- **Identity gate:** PASS on 60/60. The largest |κ̂₀.₀₅ diff| is 0.131, and all 60 are under 0.15.
- **Medians:** they move by at most 0.035 against the 40-capture medians. For example, at n 1500 and c0 0.70 the median is 1.738 against 1.758.
- **Table 3 at α 0.05:** every entry is within 0.012 of the 40-capture value.
  - The worst B rates are identical at c0 0.70, 0.75 and 0.85. At c0 0.80 it is 0.0030 (B4) against 0.0035 (B5).
  - The largest power difference is C3 at c0 0.85, 0.6630 against 0.6740.
- **Conclusion:** with this many captures, the per-n median is not the limiting error. The replicate-to-replicate spread of κ̂(c0) in Table 1 is.

## OPEN ITEMS

- The payload is `results/declcal_c0approx_res.rds` (the path the task named) and holds 120 rows. The task's "60-replicate payload" wording follows its own arithmetic slip.
- `current_status.md` was not regenerated, as the task says; the campaign's closeout covers it.
- The untracked declcalc0 files already in the tree (`results/declcalc0_inull_B{1,2,3}_res_1_2000.rds`, their logs, `declcalc0_findings.R`) were not touched and are not in this commit.
