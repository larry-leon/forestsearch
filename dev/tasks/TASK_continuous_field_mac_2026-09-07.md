# TASK (Mac Studio) — Continuous/MD path: port the field additions to the continuous twin and evaluate the current constructions

Date: 2026-09-07. Author: chat (spec). Executor: Claude Code on the Mac Studio. Approver: Larry.
Machine and branch: fresh clone or pull of `origin/feature/glm-extension`, then work on `feature/glm-extension-mac` (create from the pulled tip; push that branch; never push to `feature/glm-extension`). The Linux box is concurrently running `TASK_p30sg_nb20_2026-09-07.md` on `feature/glm-extension`; this task must not edit any file that task touches (the survival template `quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_field_m1_template.qmd`, anything under `quarto/simulations/gbsg_020/`, `R/fs_mr_inference.R`, `R/forestsearch_main.R`).
Predecessors: the continuous/MD workstream (`HANDOFF_continuous_2026-08-27_v5`; the continuous twin template and its committed bundles), the survival-side additions of 2026-09-04→07 (field block, complement field, joint pair, `return_reselection`, `fs_sim_bias_coverage()` / `fs_plot_bias_coverage()`, the campaign-tag stem, the `.refuse_if_tracked()` save guard, `FS_*` env knobs) — all in the survival template and in `R/`.

## Standing conventions this session has not yet learned

- Protocol: task documents travel via ~/Downloads → `dev/tasks/` and are committed as the first action; gates stop on failure; Gate 1 is the compute go/no-go; verify from source; records beside the results, not in `dev/tasks/`; no push to `feature/glm-extension`.
- Check ~/Downloads for stale variants of a task document before reading it.
- Simulation replicates are seeded `seed_base + sim_id`; workers run under `RNGkind("L'Ecuyer-CMRG")` (`future` with `seed = TRUE`), so a standalone reproduction of a replicate must set the same regime.
- Never write to a git-tracked path from a driver; the survival template's `.refuse_if_tracked()` guard is to be ported here before any batch save.
- Interpretation convention: bounds are read by location against clinically meaningful effect sizes, never as significance at the null; coverage tables carry Wilson intervals; bias is reported in SD units alongside the natural scale.
- Rows: naive, oracle, MR (IJ two-term), MR (field); the joint pair (separate / Bonferroni / calibrated). The IJ winner-only and winner-floor variants are excluded from every table, figure and report line (recorder columns may exist).
- Cross-machine numerics: identities against Linux-committed bundles are tolerance-based (≤ 1e-8 relative on numerics; memberships compared as sets), never `identical()`; seeds reproduce exactly.

## Purpose

Every field-method evaluation to date is survival. The gate is outcome-agnostic, and the continuous path is where the theory is cleanest (the influence linearization closest to exact for a mean difference, Λ* closest to Gaussian), so it is the right place to establish that the field's one-sided calibration, the complement block, and the joint pair are properties of the construction rather than of the Cox setting. Deliverables: the continuous twin brought to parity with the survival template (same knobs, recorder, tables, display, guard), and the current constructions evaluated on the continuous design's cells with identity anchors on everything that already existed.

## Stage 0 — Discovery (no compute, no R/ edits)

0a. Locate the continuous twin template and its committed bundles (paths, cells, `n_sims`, seeds, `sg_focus`, thresholds, the DGM and its harm orientation on the MD scale, the truth attachment for β(Ĥ)/β(Ĥᶜ), θ†, θ‡); quote the gate call and its `mr_inference_args`; list which columns the committed bundles carry (naive / oracle / FB / MR (IJ) and the coverage indicators).
0b. Confirm from source that `fs_mr_inference()` under `ci_method = "field"`, `include_complement = TRUE`, `field_complement = TRUE`, `return_reselection = TRUE` runs on the continuous specification (effect scale identity, `to_eff` the identity, harm orientation as the gate's sign convention), and what `fs_sim_bias_coverage()` assumes about the scale (it takes logs of HR-scale columns; on the MD scale a `scale = c("log", "identity")` argument is needed — the one permitted R/ change, add-only, in `R/fs_bias_coverage.R` only, default `"log"` preserving the 14-point identity fixture).
0c. Enumerate the survival-template additions to port (recorder columns `fld_H_*`, `fld_Hc_*`, joint, `p_hat_H` / top-3; the Table-2-layout rows for MR (field); the complement block's one-sided-upper convention; the field diagnostics section; the display chunk; `campaign_tag` stem; `FS_*` knobs; `.refuse_if_tracked()`; combine-mode poolability keys) and map each to its place in the continuous twin. Quote the twin's stem construction.
0d. Mac cost anchors: cores, a single-replicate timing of the twin at its committed settings; worker count = physical cores − 1 unless memory says otherwise.

Gate 0: 0a–0d resolved; STOP if the gate's field path needs any change beyond the `scale` argument, or if the continuous twin's bundles lack the columns needed for the identity anchor.

Output: `REPORT_continuous_field_stage0_<date>.md` beside the continuous twin.

## Stage 1 — Port, identities, smoke, projection

1a. Port the additions to the continuous twin (document-level); add `scale` to `fs_sim_bias_coverage()` (add-only; Rd; suite; the 14-point survival fixture still passing under the default).
1b. Identities: with all new knobs at default and `ci_method = "ij"`, 5 replicates per cell reproduce the committed continuous bundles on every pre-existing column to ≤ 1e-8 relative (cross-machine tolerance), memberships as sets; with `ci_method = "field"`: MR (IJ) columns unchanged to the same tolerance, field and complement-field finite, bound identities (`lower_1s = β̃ − q95`, complement `upper_1s = β̃ᶜ − q05`, etc.) at ≤ 1e-12, γ in range, p̂ finite; the observed Ĥ reproduced by the gate's re-selection map on the unperturbed effects on one replicate; `.refuse_if_tracked()` verified to block a tracked path and allow an untracked one.
1c. Projection on the Mac at the calibrated worker count, measured under load.

Gate 1: identities pass; projection reported. Compute go per M-5.

Output: `REPORT_continuous_field_stage1_<date>.md`.

## Stage 2 — Runs

The continuous design's committed cells (M-1), `sim_id` 1–2,000 (or the committed count if smaller — state it), seeds as committed, `ci_method = "field"`, `include_complement = TRUE`, `field_complement = TRUE`, `return_reselection = TRUE`, `ij_residual = "two_term"`, FB per M-3, campaign `cont_field_mac`. Gate 2 per cell: completeness; pre-existing columns within tolerance of the committed bundles (the pairing proof); fields finite; invariants; γ in range.

## Stage 3 — Report

`REPORT_continuous_field_<date>.md` plus rendered documents, in the survival reports' layout: per cell, both blocks, rows naive / oracle / MR (IJ) / MR (field): bias on the MD scale and in SD units, empirical SD, mean SE, SE/SD, two-sided coverage of β(·) with Wilson intervals, one-sided coverage on the exposed side (harm: lower bound in the harm direction of the MD scale; complement: the benefit-direction bound), half-widths; bound-location tables against clinically meaningful MD thresholds (state the thresholds used as reading aids); the joint table; the display for both blocks (`scale = "identity"`); regime diagnostics (p̂(Ĥ), complement fits, corr(Λ*, Λ*ᶜ), SD(β̃ᶜ)/naive SE, λ-SDᶜ/naive SE); the two summary-table formats (bias / SD / SE / SE-to-SD / two-sided / one-sided) for Ĥ and Ĥᶜ. Reading criteria (Larry's, not gates): the field's one-sided bounds on both blocks against the survival ranges (harm 0.92–0.98; complement 0.91–0.96); IJ two-term's SE/SD; whether all display points lie on the Gaussian reference (the continuous path should show the smallest departures). Findings in the record; no task proposed unless something blocks.

## Decisions (defaults in brackets)

- M-1 Cells: the continuous twin's committed cells [default]; Stage 0 lists them for confirmation.
- M-2 Replicates: the committed count [default]; 2,000 if the committed count is smaller and cost allows.
- M-3 FB: joined from committed bundles where they exist, otherwise none; never re-run [default].
- M-4 Workers: physical cores − 1 on the Mac, measured [default].
- M-5 Compute: go at Gate 1 on the measured Mac projection, or a pre-authorization with a ceiling set from it.
- M-6 Branch: `feature/glm-extension-mac`, pushed as commits land; merge into `feature/glm-extension` on Linux when both sides are quiet — Larry's step.

## Done means

Continuous twin at parity with the survival template; `scale` argument landed add-only with the survival fixture intact; Stage 3 report and rendered documents committed on `feature/glm-extension-mac` with Gate 2 records beside them; the branch pushed; no file touched that the Linux task edits.
