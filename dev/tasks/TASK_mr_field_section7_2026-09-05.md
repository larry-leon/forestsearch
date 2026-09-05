# TASK — MR (field) on the Section 7 forestsearch cells: key test and the combined FB / MR (IJ) / MR (field) template

Date: 2026-09-05. Author: chat (spec). Executor: Claude Code. Approver: Larry.
Predecessors: `TASK_mr_field_vs_guohe_2026-09-05.md` (ci_method = "field" landed add-only at 87880a24; record `quarto/GuoHe/REPORT_mr_field_vs_guohe_2026-09-05.md`); the paper's Section 7 / Table S4 simulation documents under the simulations directory (`sim_fs_maxeffCons_mr_m1_*.qmd`, `sim_fs_maxeffCons_fb_mr_m1_*.qmd`).

## Protocol

- First action: copy this file to `dev/tasks/` and commit. Do not push.
- Gates stop on failure; Gate 1 is the compute go/no-go (may be pre-authorized).
- Verify from source; quote paths, signatures and line numbers.
- Reports beside the simulation documents, next to the existing results bundles' directory.
- The committed simulation documents and their results bundles are read-only inputs. No re-run of anything committed: FB is joined, never recomputed.
- No change to `R/` is expected. If Stage 0 finds the wrapper does not pass `ci_method` through, or drops the `field` element, the fix is add-only (pass-through / return-field) and is called out as such in the Stage 1 record; anything touching the gate's arithmetic is a STOP.

## Purpose

Test the field-calibrated interval where it has not yet been tested: under the paper's own identifier — enumerated overlapping family, effect screen, consistency requirement, size tie-break — on the two Section 7 cells that bracket the regimes. Deliver, as the same artifact, the template that runs FB, MR (IJ) and MR (field) together going forward.

Cells (decision F1):
- **Key cell:** borderline null, `target_hr_harm = 1.0`, n = 500, knoise 0, FS `maxeffCons` — the Table S4 cell; MR (IJ) there: conditional coverage 0.98, SE/SD ≈ 1.6, retained bias ≈ +0.26 log-HR.
- **Separation check:** genuine harm, `target_hr_harm = 1.75`, same settings — MR (IJ) 0.94 conditional / 0.85 marginal; the over-shrinkage risk case for "field".

Comparators per replicate, all on the same data and the same selected subgroup: naive, oracle, FB (key cell only, joined from the committed `fb_mr` bundle), MR (IJ), MR (field). Targets: β(Ĥ) conditional (primary), θ† marginal, θ‡ CDE, oracle. Blocks: Ĥ and Ĥᶜ.

## Stage 0 — Discovery (no compute, no R/ edits)

0a. `forestsearch()` → gate: quote where `mr_inference_args` is forwarded to `fs_mr_inference()`, confirm `ci_method = "field"` passes through untouched, and where the gate's return lands in the fit object (the `debias_gate` element). Confirm the `field` element survives into the fit object, and whether a `field` block exists for the complement when `include_complement = TRUE` (decision F3 depends on this).
0b. Committed bundles for the two cells: paths of `results/<stem>_res_1_1000.rds` for the MR-only h10 and h175 documents and of the h10 `fb_mr` bundle; their `n_sims`, `nb_boots`, `mr_draws`, `seed_base`; confirm the seed scheme (`seed_base + sim_id` for DGM, search and noise) is the same in the MR-only and FB documents, and that on common `sim_id`s the naive and oracle columns are identical across the two committed bundles (this licenses the FB join).
0c. Recorder: quote the per-replicate columns for MR (est, SE, bounds, `mr_ok`, `mr_harm_flag`, complement columns), the four coverage targets and how β(Ĥ) truth is attached (`fs_attach_betaHhat` or equivalent), and the detection/selection columns.
0d. Cost: per-replicate `fit_mr_secs` distribution from the committed bundles; the FB-side cost is irrelevant (joined).

Gate 0: 0a–0d quoted; STOP if the wrapper needs an arithmetic change, if the seed schemes differ between the MR-only and FB documents, or if the naive/oracle cross-bundle identity fails.

Output: `REPORT_mr_field_s7_stage0_<date>.md`.

## Stage 1 — Template, identities, smoke, projection

1a. Template. Derive `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` from the committed `fb_mr` document with knob deltas only: `ci_method = "field"`; recorder extended to capture the field block for Ĥ (and Ĥᶜ per F3): `est2`, `lambda_mean`, `lambda_sd`, the seven quantiles, one-sided bound, two-sided quantile interval, SE-type interval, draw-usage counts; an `fb_mode` knob with values `run` (nb_boots ≥ 1), `join` (read the committed FB bundle and join by `sim_id`), `none`; estimation/coverage tables in the Table S4 layout gaining MR (field) rows for every target and block; the FB↔MR agreement panel gaining `est2`. The committed documents stay untouched; the template carries a `campaign_tag` so it never writes into an existing stem.
1b. Identities (machine-checked):
- On 5 replicates of the key cell at the committed seeds: naive, oracle, detection flag and selected-subgroup definition identical to the committed MR-only bundle; MR (IJ) est / SE / bounds identical; FB columns joined for the same `sim_id`s with naive identical across bundles; MR (field) finite.
- Same on 5 replicates of the h175 cell (no FB).
- Two-sided and one-sided field coverage indicators recomputed from the stored bounds agree with the recorder's.
1c. Projection on a loaded machine at ~60 workers (the host is memory-bandwidth-bound for these runs; 120 workers gained nothing on 2026-09-05): per-replicate cost including the search, and core-hours / wall for 2 cells × 2,000 replicates.

Gate 1: 1b passes; projection reported. Compute go per F7.

Output: `REPORT_mr_field_s7_stage1_<date>.md`.

## Stage 2 — Run

Both cells, `sim_id` 1–2,000 as two seed-disjoint batches each (1–1,000 and 1,001–2,000), ~60 workers, then `run_mode = "combine"`. FB joined for the key cell on `sim_id` 1–1,000 (the committed range); FB absent for 1,001–2,000 and for h175, recorded as such. No mid-run changes; a failed batch stops and reports.

Gate 2: completeness per batch; on `sim_id` 1–1,000 naive / oracle / detection / MR (IJ) identical to the committed bundles; FB join count equals the committed FB detections; MR (field) NA count zero or documented.

## Stage 3 — Report

`REPORT_mr_field_s7_<date>.md` plus the rendered combined document. Per cell, Table S4 layout: rows naive / oracle / FB / MR (IJ) / MR (field), columns mean estimate, relative bias against each target, empirical SD, mean SE, SE/SD, two-sided 95% coverage with Wilson intervals (paper convention, primary), one-sided lower coverage (supplementary), mean and median length; both blocks. Add: detection rate; retained bias of β̃ and est2 against β(Ĥ) in log-HR and in SD units; λ-sd/SD; and the paired panel est2 vs FB on the key cell. Interpretation limited to the tables; the reading criteria below are for Larry, not gates.

Reading criteria (proposed, not gates): "field" holds up if, on both cells and both blocks, conditional coverage lies within Monte Carlo error of nominal or above it (≥ 0.93 at ~1,300 detections) with SE/SD in [0.9, 1.2], est2's retained bias is below β̃'s, and marginal coverage on the h175 cell is not below MR (IJ)'s 0.85 by more than Monte Carlo error. A failure on the complement alone is a scoping result (F3), not a method result.

## Decisions required from Larry (defaults in brackets)

- F1 Cells: h10 (key) + h175 (separation) [default]; h075 protective and the n-sweeps deferred to a second stage.
- F2 Replicates: 2,000 per cell in two batches [default]; 1,000 if cost demands.
- F3 Complement: MR (field) for Ĥᶜ if the gate already produces it under `include_complement = TRUE` or it is add-only cheap; otherwise Ĥᶜ stays on IJ with a note [default].
- F4 Primary convention: two-sided 95% (paper's tables), one-sided supplementary [default].
- F5 FB: joined from the committed h10 bundle only; no FB run anywhere [default].
- F6 Template adoption as the standard document: Larry's call after the Stage 3 report [—].
- F7 Compute: go at Gate 1, or an unattended pre-authorization with a wall ceiling.

## Done means

Stage 3 report and rendered document committed, Gate 2 record beside them, template committed under its own stem with the committed documents untouched, every Stage 1 identity listed with its measured value, branch left unpushed for Larry.
