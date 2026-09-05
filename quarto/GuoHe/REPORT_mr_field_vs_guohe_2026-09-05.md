# REPORT — MR (field) vs MR (IJ) vs Guo & He (2021) on the replicated Guo–He simulations

Date: 2026-09-05. Task: `dev/tasks/TASK_mr_field_vs_guohe_2026-09-05.md`; run unattended under Larry's pre-authorization.
Full tables, figures and provenance: `mr_field_vs_guohe.qmd` (rendered `mr_field_vs_guohe.html`, untracked per the `GuoHe/*.html` convention; reproduces from the 16 committed `mr_field_vs_guohe_*.rds` bundles). Gate 2 record: `REPORT_mr_field_gate2_2026-09-05.md` (GREEN; includes the mid-run watchdog kill and identical-settings resume of one cell). Stage 1 identities with measured values: `REPORT_mr_field_stage1_2026-09-05.md`. PoC provenance: `dev/tasks/POC_mr_interval_alternatives_2026-09-05.md` (first PoC) and `dev/tasks/poc2_results_2026-09-05.csv` (authoritative reference for the shrunk-field / second-order design and the 1b values).

Design in one line: the 2026-09-04 paired comparison re-run with one gate call per replicate at `ci_method = "field"` — MR (IJ) provably unchanged (`identical()` to the 2026-09-04 bundles on all 32,000 replicates) — adding the MR (field) row: `est2 = β̃ − λ̄` (E1/E3), primary one-sided bound `β̃ − q95(Λ*)` (E4), two-sided quantile interval supplementary, R_out/R_in = 1000/500 (E2), diagnostics joined per E5.

## Headline numbers (log-HR; one-sided 95%; G&H at r = 1/30; per-cell Wilson intervals in the qmd)

| Cell | bias: IJ / field / G&H | cover: IJ / field / G&H | margin: IJ / field / G&H |
|---|---|---|---|
| t35 β₂=0.0 | +0.028 / **+0.010** / +0.001 | 0.999 / **0.945** / 0.949 | 0.50 / **0.30** / 0.27 |
| t35 β₂=0.1 | +0.020 / +0.004 / −0.008 | 0.999 / 0.949 / 0.956 | 0.50 / 0.30 / 0.26 |
| t35 β₂=0.2 | +0.016 / +0.003 / −0.020 | 0.999 / 0.945 / 0.950 | 0.51 / 0.30 / 0.26 |
| t35 β₂=0.3 | +0.005 / −0.003 / −0.038 | 1.000 / 0.936 / 0.947 | 0.52 / 0.31 / 0.26 |
| t35 β₂=0.4 | −0.009 / −0.012 / −0.061 | 0.997 / 0.936 / 0.953 | 0.54 / 0.31 / 0.26 |
| t35 β₂=0.5 | −0.013 / −0.011 / −0.072 | 0.997 / 0.942 / 0.954 | 0.56 / 0.31 / 0.26 |
| t6 k=2 | +0.026 / +0.008 / +0.000 | 0.998 / 0.944 / 0.950 | 0.50 / 0.30 / 0.27 |
| t6 k=6 | +0.072 / +0.030 / +0.004 | 0.987 / 0.941 / 0.948 | 0.40 / 0.29 / 0.22 |
| t6 k=10 | +0.088 / +0.038 / +0.004 | 0.982 / 0.944 / 0.953 | 0.37 / 0.29 / 0.20 |
| t6 k=12 | +0.093 / +0.041 / +0.002 | 0.972 / 0.942 / 0.949 | 0.36 / 0.29 / 0.20 |
| t7 β₂=0.0 | +0.045 / +0.030 / +0.011 | 0.998 / 0.940 / 0.954 | 0.58 / 0.31 / 0.31 |
| t7 β₂=0.1 | +0.040 / +0.026 / +0.005 | 0.998 / 0.933 / 0.954 | 0.58 / 0.31 / 0.31 |
| t7 β₂=0.2 | +0.039 / +0.025 / −0.000 | 0.999 / 0.938 / 0.958 | 0.59 / 0.31 / 0.31 |
| t7 β₂=0.3 | +0.033 / +0.021 / −0.010 | 0.999 / 0.937 / 0.962 | 0.60 / 0.31 / 0.30 |
| t7 β₂=0.4 | +0.028 / +0.017 / −0.021 | 0.999 / 0.941 / 0.967 | 0.60 / 0.32 / 0.30 |
| t7 β₂=0.5 | +0.015 / +0.007 / −0.041 | 0.999 / 0.951 / 0.973 | 0.61 / 0.32 / 0.31 |

(Naive and the other three r columns in the qmd; coverage MCSE ≈ 0.005. "Margin" = est − one-sided lower; for the field row est = est2, lower = β̃ − q95.)

## What the numbers show

1. **Coverage at half the width.** MR (field)'s one-sided coverage sits at 0.933–0.951 across all 16 cells — near-nominal everywhere the current interval over-covers at 0.97–1.00 — with margins of 0.29–0.32 against MR (IJ)'s 0.36–0.61: the field calibration removes roughly half the interval width, landing within 0.00–0.10 of the replicated G&H's width (and equal to it on the nested family) while keeping, unlike G&H at small r, a point estimate that never over-corrects (field bias is ≥ −0.012 everywhere; G&H reaches −0.072).
2. **Second-order bias removal.** In tie scenarios `est2` cuts the retained naive optimism from β̃'s 0.26–0.31 to **0.08–0.13** on the disjoint families and from 0.37 to **0.24** on the nested family (est2 bias +0.010 to +0.041 vs β̃'s +0.026 to +0.093; per-cell MCSE ≈ 0.004). The graded cells shrink toward zero without sign flip. The nested family remains the hardest case, consistent with the first PoC's reading that the excess is a property of unequal candidate variances under nesting.
3. **Calibration of the field SD.** `lambda_sd`/SD(est2) runs 0.86–1.09 — near 1 in tie/null and nested cells, dipping to 0.86–0.93 in the moderately separated t35 cells, which is also where the one-sided coverage shows its mild low points (0.936) and the supplementary two-sided quantile interval dips to 0.895–0.902 (elsewhere 0.92–0.97). This is the first PoC's S3 signature, much reduced by the shrunk field but not eliminated; the one-sided bound (the primary comparison, E4) stays ≥ 0.933 everywhere.
4. **Nothing upstream moved.** MR (IJ) columns are `identical()` to the 2026-09-04 bundles on every replicate, the naive column is `identical()` to the replication bundles, G&H is stored arithmetic — every difference in this report is attributable to the field construction alone.
5. **Cost.** The field pass adds ~0.03 s (k = 2) to ~0.24–0.36 s (nested) per replicate inside a gate call of ~0.1–2.4 s — a few percent to ~15% on top of the MR pass, and no M_eff Monte Carlo ran (E5). Whole grid: ~40 min wall in total across the interrupted run and the one-cell resume (per-cell walls in the qmd cost table).

Interpretation stops at the tables. No task is proposed: nothing blocked, and the method disposition (whether "field" should become a default) is Larry's call on this record.

## Provenance and decisions

- Engine: `forestsearch:::fs_mr_inference` at `87880a24` (`ci_method = "field"`, add-only; default path identical component-for-component to `736088e3`, wall-clock timing field excluded; suite totals unchanged at 0 fail / 4941 pass).
- Bundles: `mr_field_vs_guohe_*.rds` × 16 (driver `mr_field_vs_guohe_run.R`); read-only inputs `guohe_repro_*.rds` and `mr_vs_guohe_*.rds`, both unmodified.
- Decisions: E1–E5 defaults as specified; E6 resolved by the pre-authorization (Gate 1 green, projection 47 core-h / ~23 min ≪ 12 h).
