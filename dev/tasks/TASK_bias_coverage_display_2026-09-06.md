# TASK — Bias-vs-coverage display as a standard output: package functions and a cross-cell summary document

Date: 2026-09-06. Author: chat (spec). Executor: Claude Code. Approver: Larry.
Reference: `bias_vs_coverage.png/.pdf` and `bias_coverage_points.csv` (chat outputs, 2026-09-06; the CSV is the identity fixture — 14 rows: cell, method, b, r, c1, c1_pred, c2, c2_pred). Source records: `REPORT_mr_field_s7_2026-09-05.md`, `REPORT_mr_field_ocmap_2026-09-05.md`; bundles under `quarto/simulations/gbsg_020/results/` for campaigns `s7` and `map1`.

## Protocol

- First action: copy this file and the CSV to `dev/tasks/` and commit. Do not push.
- R/ change: adds two functions; no existing behaviour changes; suite unchanged. Classification: adds code.
- Verify from source: the recorder's column names in the committed bundles (`nv_H_*`, `mr_H_*`, `fld_H_*`, `betaHhat_H`, `detected`), and the template's existing diagnostics section (its one-sided convention and Wilson helper), before writing the functions.
- Committed bundles and documents are read-only inputs.

## The display

Three panels, harm block by default: (1) one-sided 95% lower-bound coverage against residual bias in SD units, with Gaussian-reference curves Φ(1.645·r − b) for r ∈ {1, 1.25, 1.5, 2}; (2) two-sided 95% coverage against bias, curves Φ(1.96·r − b) − Φ(−1.96·r − b); (3) observed coverage against the Gaussian reference at each point's own (b, r), with the diagonal. Points: one per (cell, estimator); IJ and field distinguished by marker; nominal line at 0.95; cells labelled. The reference is bookkeeping, not theory: b = mean log error / empirical SD, r = mean reported SE / empirical SD, both from the cell's own replicates.

## Functions (R/, add-only)

1. `fs_sim_bias_coverage(results, block = c("H", "Hc"), estimators = c("naive", "mr", "fld"), level = 0.95, target = "betaHhat")`: from a combined bundle's `results` (detected replicates), returns one row per estimator with `n`, `bias_log`, `sd_emp`, `se_mean`, `b` (= bias_log / sd_emp), `r` (= se_mean / sd_emp), `cov1`, `cov1_wilson_lo/hi`, `cov2`, `cov2_wilson_lo/hi`, `cov1_ref`, `cov2_ref`. One-sided bound per estimator: field from `fld_H_lo1s`; naive and IJ from exp(log est − z₀.₉₅·SE) (the template's existing convention). Truth from `betaHhat_<block>` (or the oracle/CDE/marginal columns when `target` says so).
2. `fs_plot_bias_coverage(tbl, labels = TRUE, curves = c(1, 1.25, 1.5, 2))`: draws the three panels (ggplot2 + patchwork or base grid; whichever the package already imports) from a stacked table with a `cell` column; returns the plot object.

## Documents

3. New `quarto/simulations/gbsg_020/summary_bias_coverage.qmd`: knobs — a campaign glob (default `*_{s7,map1}_combine_*.rds`) and block; reads every matching combined bundle, calls (1) per bundle, stacks, draws (2), and prints the table. Self-contained: no reads from `dev/`.
4. Combined template: a new chunk in the MR (field) diagnostics section that calls (1) on its own cell and (2) with the two points, so every future render carries the display for that cell.

## Stages

- **Stage 0** (no compute): quote the recorder columns and the template's one-sided/Wilson code; confirm the s7 and map1 combined bundles present (7 cells: h100/h175 s7; h150, h075, h100-n1000, h175-knoise3, h150-n1500 map1); confirm ggplot2 availability in the package imports.
- **Stage 1**: implement (1)–(2); Rd; suite. **Identity:** running (1) on the seven committed bundles reproduces `bias_coverage_points.csv` — b, r, cov1, cov2 within 0.005 and cov1_ref/cov2_ref within 0.001 — for all 14 rows; the plot renders. Render (3); add (4) and re-render one s7 cell from its committed bundle in combine mode (no new replicates) to confirm the chunk works. Gate 1: identity passes; renders present.
- **Stage 2**: commit functions, documents, rendered HTML per the directory's convention, and `REPORT_bias_coverage_display_<date>.md` (identity table, the rendered figure, one paragraph on the reading). No compute beyond renders.

## Decisions (defaults in brackets)

- J1 Default block: harm Ĥ [default]; complement available by argument.
- J2 Curves: r ∈ {1, 1.25, 1.5, 2} [default].
- J3 Guo–He bundles: out of scope here (different recorder); a later extension if wanted [default: out].
- J4 Function names as above [default].

## Done means

Functions landed add-only with the 14-point identity recorded; summary document and template chunk committed; report beside the s7 records; branch left unpushed.
