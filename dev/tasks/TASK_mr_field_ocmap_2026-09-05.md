# TASK — OC map: MR (IJ) and MR (field) on cells not yet covered (committed template, no code changes)

Date: 2026-09-05. Author: chat (spec). Executor: Claude Code, unattended. Approver: Larry.
Sequencing: start only after `TASK_mr_field_uniform_2026-09-05.md` has stopped at Gate 1 (its Stage 1 projection must be measured on an unloaded machine). Do not touch that task's state.
Predecessors: `TASK_mr_field_section7_2026-09-05.md` (template `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd`, campaign `s7` records, all in `quarto/simulations/gbsg_020/`); the committed coverage-sweep documents and bundles for the consistency detector (`maxeffCons_mr_coverage_sweep_h*_knoise*.qmd`, results under their `mr_sweep`/results directories).

## Protocol

- First action: copy this file to `dev/tasks/` and commit. Do not push.
- Gates stop on failure. Gate 1 is pre-authorized under the ceiling in H-Q6.
- No change under `R/`. No change to the committed template or documents; every run is an environment-variable render of the template under a new campaign tag.
- Verify from source; reports beside the s7 records.

## Purpose

Extend the operating-characteristics map of MR (IJ) and MR (field) — the same gate call, `ci_method = "field"`, yielding both rows plus the one-sided field bound — to regimes not yet covered: where IJ's conditional coverage is expected near-nominal without being over-wide (moderate harm), the disclosed failure regime (protective truth), larger n, and a larger candidate family (noise covariates). All cells use the M1 design of the committed documents, the maxeffCons identifier, and the template's defaults; FB is `none` throughout (no committed FB for these cells).

Cells, in priority order (decision H-Q1), 2,000 replicates each (`sim_id` 1–2,000, two batches, seeds `8316951 + sim_id`), campaign tag `map1`:

| # | target_hr_harm | n | knoise | env overrides | why |
|---|---|---|---|---|---|
| 1 | 1.50 | 500 | 0 | `FS_S7_HR=1.5` | moderate harm: the regime where IJ should be near-nominal and the field's two-sided interval is most stressed |
| 2 | 0.75 | 500 | 0 | `FS_S7_HR=0.75` | protective truth: conditional coverage of IJ reported well below nominal in the sweep records |
| 3 | 1.00 | 1000 | 0 | `FS_S7_HR=1.0 FS_S7_N=1000` | does the tie-regime picture at n = 500 persist at n = 1000 |
| 4 | 1.75 | 500 | 3 | `FS_S7_HR=1.75 FS_S7_KNOISE=3` | larger family (noise covariates): M_eff and the field under a bigger competition |
| 5 | 1.50 | 1500 | 0 | `FS_S7_HR=1.5 FS_S7_N=1500` | moderate harm at large n |

## Stage 0 — Discovery (no compute)

0a. For each cell, identify whether a committed bundle with the same identifier, DGM settings and seed scheme exists (the consistency-detector sweep documents at h075/h10/h15/h20 and their per-cell bundles; the knoise3 documents). Quote paths, `n_sims`, `seed_base`, `sg_focus`, `subgroup_method`, thresholds, and confirm whether the seed scheme matches the template's (`seed_base + sim_id`, noise `+ 1000000L`). Where it matches, that bundle is the identity anchor for the MR (IJ) and naive columns on shared `sim_id`s (tolerance 1e-12, memberships not strings, enumerated cross-vintage flips reported individually as on 2026-09-05).
0b. Confirm the template's DGM calibration handles `target_hr_harm = 0.75` (protective) and `n = 1000/1500` without edits (`calibrate_k_inter`, `n_super`), and that the stem for each cell is new (`h150`, `h075`, `h100_n1000`, `h175_knoise3`, `h150_n1500`, campaign `map1`).
0c. Cost anchors from the s7 cost table (h100/h175 at 100 workers) and, for n = 1000/1500 and knoise 3, from the committed documents' own timing lines if present.

Gate 0: 0a–0c quoted; a cell whose DGM calibration fails or whose stem would collide is dropped from the list with a note, not fixed.

Output: `REPORT_mr_field_ocmap_stage0_<date>.md`.

## Stage 1 — Smoke and projection

1a. 5 replicates per cell (committed seeds, `campaign=map1smoke`): render completes; MR (IJ) and naive columns identical (≤ 1e-12) to the anchor bundle where one exists; field finite; interval invariants hold.
1b. Projection at 100 workers from the smoke timings and the s7 anchors, per cell and cumulative in priority order.

Gate 1 (pre-authorized): if 1a passes for a cell and the cumulative projection through that cell is within the H-Q6 ceiling, that cell runs in Stage 2; cells beyond the ceiling are deferred and listed. If 1a fails for a cell, that cell is dropped with the failure quoted and the rest proceed.

Output: `REPORT_mr_field_ocmap_stage1_<date>.md`.

## Stage 2 — Runs

Cells in priority order, two batches each, then combine; 100 workers; hard timeout per H-Q6; a failed batch stops that cell and reports, the next cell proceeds. Gate 2 per cell: completeness 2,000/2,000; identity to the anchor bundle on shared `sim_id`s where one exists; field NA count zero or documented; interval invariants.

## Stage 3 — Report

`REPORT_mr_field_ocmap_<date>.md` plus the rendered combined documents. Per cell, harm and complement blocks, in the s7 report layout: rows naive / oracle / MR (IJ) / MR (field); bias against β(Ĥ) in log-HR and SD units, empirical SD, mean SE, SE/SD, two-sided coverage of β(Ĥ) with Wilson intervals, one-sided coverage, z·SE on the log scale (not HR-scale lengths), coverage of θ† and θ‡ reported but not scored; detection rate; M_eff and p̂_Ĥ where recorded. Then one cross-cell table: for each cell, IJ two-sided coverage and SE/SD beside field one-sided coverage, field two-sided coverage and λ-SD/SD. Findings in the record; no task proposed unless something blocks.

## Decisions (defaults in brackets; Larry is offline — defaults apply)

- H-Q1 Cells and order: as tabled [default].
- H-Q2 Replicates: 2,000 per cell [default].
- H-Q3 Workers: 100 [default].
- H-Q4 FB: none [default].
- H-Q5 Anchors: sweep bundles where seed-compatible; otherwise no identity gate for that cell beyond completeness and invariants [default].
- H-Q6 Compute ceiling: cumulative projection ≤ 5 h wall; hard timeout 6 h; cells beyond the ceiling deferred [default].

## Done means

Stage 3 report and rendered documents committed, Gate 2 records beside them, every run under campaign `map1`, committed documents and bundles untouched, branch left unpushed; end with a one-paragraph summary listing cells completed, deferred, or dropped, and the commit range.
