# TASK — Approximate c0 table from 60 field captures (runs between the R change and the campaign)

**Repo / branch:** `larry-leon/forestsearch`, `feature/glm-extension`. **Machine:** pop-os, 64 workers.
**Prerequisite:** `TASK_declcal_c0_rchange_2026-09-22.md` green and installed (`devtools::install()`).
**Kind:** 60 replicate searches (≈1 min wall on 64 workers; ≈25 min serial) plus read-only post-processing of
the committed `declcal` payloads. **No `R/` change.** Small compute, stated: 60 searches at B = 2000.

## 0. First action
Copy this file to `dev/tasks/TASK_declcal_c0_approx_2026-09-22.md`; `git add` that path; commit
`docs(tasks): add c0 approximate-table task (2026-09-22)`. Assert `declaration_c0` is a formal of the **installed**
`fs_mr_inference` (default `NULL`); if not, `devtools::install()` from HEAD and re-assert; STOP only if still failing.

## 1. The captures
- Replicates **1–20** of each of B1, B4 (n 500), B2, B5 (n 1000), B3, B6 (n 1500): 60 in all, at the campaign's
  own seeds (`8316951 + rep`), c1 = c2 = 1.0, floors as the template, **B = 2000**, centred Poisson,
  `declaration_c0 = c(0.70, 0.75, 0.80, 0.85)`, `keep_field_matrix = FALSE`.
- Transplant the per-replicate engine from the committed `declcal_run.R`; change only the replicate range, `B`,
  and the `declaration_c0` argument. Record the lines changed.
- Per replicate record: `cell, n, rep, max_T_pre, G_pre, kappa_hat_05, kappa_hat_10` (unshifted) and, per `c0`,
  `kappa_hat_05_c0, kappa_hat_10_c0, pstar_implied_05_c0, Mstar_c0_q90/q95/q99, fw_1645_c0`.

## 2. Identity gate (STOP on failure)
For all 60 replicates, `max_T_pre` and `G_pre` must be **identical** to the committed
`results/declcal_inull_<cell>_res_1_2000.rds` values at `81752681`, and the unshifted `kappa_hat_05` within the
B = 2000 vs B = 500 Monte-Carlo difference (report the differences; expect |diff| < 0.15). Disagreement on
`max_T_pre` or `G_pre` is a STOP.

## 3. The approximation
- For each n and each `c0`: **median** `kappa_hat_05_c0` and `kappa_hat_10_c0` over the 20 replicates, with
  min / IQR / max reported beside it (the spread is the approximation's error bar).
- Apply each median as a **fixed cutoff** to the stored `max_T_pre` of all 2,000 replicates in every B and C cell
  at that n (committed payloads, read-only): rate = `mean(max_T_pre >= median_kappa)`, Wilson 95%.
- Also report, per n and `c0`, the mean `fw_1645_c0` over the 20 replicates.

## 4. Output — `scripts_dinamr/logs/declcal_c0_approx.txt` and a short record
`dev/reports/REPORT_declcal_c0_approx_2026-09-22.md`
- Table 1: n × c0 → median κ̂₀.₀₅ (min–max), median κ̂₀.₁₀, implied p\*, mean fw_1645.
- Table 2: cells × {p\* 0.90 as executed (committed), c0 0.70, 0.75, 0.80, 0.85, c0 = c2 (committed)} at α 0.05,
  Wilson intervals; B rows then C rows; no Block A. Same at α 0.10.
- Table 3: c0 as rows → worst uniform-benefit rate (which cell), HR 1.5 power n 1000 / 1500, HR 2.0 power
  n 1000 / 1500 — with the committed fixed-p\* (k 2.0) and c0 = c2 rows beside them.
- State plainly that these are plug-in fixed-cutoff rates using a per-n median κ̂, not the per-replicate rule; the
  campaign task supplies the exact version.
- Commit the table, the record, the 60-replicate payload (`results/declcal_c0approx_res.rds`) and the script,
  explicit paths. No push. No `current_status.md` regeneration (the campaign's closeout covers it).
