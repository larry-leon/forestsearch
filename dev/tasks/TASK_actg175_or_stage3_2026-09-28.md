# TASK — ACTG 175 binary (odds ratio): Stage 3 summarization of the 18 committed cells (2026-09-28)

**Repo:** `forestsearch`, branch `feature/glm-extension`.
**Kind:** summarization of committed results. **No simulation is run, no bundle is written, moved or deleted.**
R reads the committed combined bundles and writes summary files only.
**Source of truth for the design and settings:** `quarto/simulations/actg175/binary_020/current_status.md`,
`sim_fs_mr_field_or_template.qmd`, and `dev/reports/REPORT_glm_campaign_inventory_2026-09-28.md` §2.

**Why:** the three binary campaigns (`orfs`, `orgrf`, `ordina`; 18 cells × 1,000 replicates) are compute-complete and
committed, but their Stage 3 synthesis was never done: no coverage or bias table exists for any cell. The
manuscript needs one, in the same form as the continuous campaigns' extracts.

**First action:** copy this file into `dev/tasks/` and commit it.

## Gates (stop on failure)
- No tracked file modified before any step; untracked files present at the start are listed and left alone.
- The 18 committed bundles `…_orfs_d5000/…_combined_1_1000.rds`, `…_orgrf_…`, `…_ordina_…` exist and are tracked.
- md5 of every bundle read is taken before and after; all unchanged at the end.

## Step 1 — the summary source
`binary_020/summary_actg175_or.qmd` reads the superseded `*_combined_1_2000.rds` layout. Point it at the committed
`*_combined_1_1000.rds` bundles of the three campaigns (18 cells); change nothing else in its logic unless needed to
read the current layout, and record each change.

## Step 2 — the extract
Write `binary_020/or_metrics.csv` in **the same schema as** `continuous/md_field_metrics.csv` (and
`md_grf_metrics.csv`'s `identifier` column): one row per identifier × cell × block (region Ĥ, complement Ĥᶜ, joint) ×
estimator × metric, with Monte Carlo standard errors. Metrics at least: declaration rate; bias (on the log-odds
scale, against the planted target as the continuous extracts define it); one-sided coverage of the field lower bound
on β(Ĥ) and of the field-s upper bound on β(Ĥᶜ); joint coverage of the Bonferroni pair; two-sided coverage of the
unadjusted, oracle and IJ two-term intervals as references. Bias and coverage over declaring replicates, as in the
continuous extracts. Write `binary_020/COLUMNS_or.md` defining every column, modelled on `COLUMNS_md_field.md`.

## Step 3 — the design's truths
From each design point's bundles (`meta` / `truth`), read the stored marginal odds ratios in the planted region, in
its complement, and over the whole trial, and the planted prevalence. Report them per design point (OR 0.75 / 1.0 /
1.5). If a field is absent, write `not established from source`; do not compute one.

## Step 4 — render and report
- Render `summary_actg175_or.qmd` to its html.
- Write `binary_020/REPORT_actg175_or_2026-09-28.md`: the design and settings (cite the template and catalog); the
  Step 3 truths; **one comparative table per block** (rows identifier × cell; columns declaration, bias, and each
  coverage), then a reading of at most eight lines; the package build that produced the bundles (from `meta`), and a
  one-line note that it predates the 2026-09-23 admission alignment. No recommendations.
- Update `binary_020/current_status.md` in place: Stage 3 done, with the extract and report paths.
- Copy the report, `or_metrics.csv` and `COLUMNS_or.md` to `~/Downloads/actg175_or_stage3_2026-09-28/`.

## Not in this task
The closeout deletion of the smoke and calibration bundles (catalog §3.5), and the two untracked `redes` /
`relaunch` directories: leave all of them exactly as they are.

## Post-conditions (stop on failure)
- Bundle md5s unchanged; no bundle added, moved or deleted.
- Added or modified: this task file, `summary_actg175_or.qmd` (and its html), `or_metrics.csv`, `COLUMNS_or.md`, the
  report, `current_status.md`. Nothing else.
- Every cell of the 18 appears in `or_metrics.csv` for every identifier it belongs to.
- Commit. **Do not push.**
