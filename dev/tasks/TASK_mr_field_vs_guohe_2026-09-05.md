# TASK — Field-calibrated MR interval ("MR (field)") and the three-way Guo–He comparison

Date: 2026-09-05. Author: chat (spec). Executor: Claude Code. Approver: Larry.
Predecessor: `dev/tasks/TASK_mr_vs_guohe_2026-09-04.md` + addendum A; Stage 3 record `quarto/GuoHe/REPORT_mr_vs_guohe_2026-09-04.md`; bundles `mr_vs_guohe_*.rds` (16, tracked).

## Protocol

- First action: copy this file to `dev/tasks/` and commit. Do not push.
- Gates stop on failure; Gate 1 is the compute go/no-go and stops to ask.
- Verify from source; quote paths, signatures and line numbers. Nothing in this document about internals is authoritative.
- Reports beside the existing `REPORT_mr_vs_guohe_*` files in `quarto/GuoHe/`.
- Do not touch the Guo–He replication code or bundles, and do not modify the existing `mr_vs_guohe_*.rds` bundles; the new run writes new bundles.

## Status of the method change (called out separately)

This task adds a new interval method to the MR gate. Classification: **changes the method** — it is a method proposal, requested by Larry for evaluation on 2026-09-05 after the proof of concept (`POC_mr_interval_alternatives_2026-09-05.md`, `poc2_results_2026-09-05.csv`). Implementation constraint: add-only. A new `ci_method` value with new return fields; the default `ci_method = "ij"` output must remain byte-identical (exact comparison of the return object on fixed seeds); the existing test suite must pass unchanged. No change to the two-term point estimate β̃, to the multiplier draws, to re-selection, or to the IJ variance.

## The method: MR (field)

Inputs the gate already holds after its per-candidate fits: the candidate effect vector β̂ (length K, the kept family), the influence matrix B_eff (n × K), the observed winner Ĥ, the two-term estimate β̃(Ĥ), and the re-selection routine S (here `reselection = "maxeff"` with unrestricted admission; in general the configured rule with thresholds held fixed).

1. Shrunk field: `w = β̂; w[Ĥ] = β̃(Ĥ)` (decision E1).
2. Gaussian-multiplier perturbations, so no Cholesky and no Σ̂ formed explicitly: `ζ = B_effᵀ ξ` with `ξ ~ N(0, I_n)`. Draw R_out outer columns ζ* and R_in inner columns ζ' (decision E2: 1,000 / 500), from a seed derived from the gate's seed.
3. For each outer draw r: `v_r = w + ζ*_r`; `G_r = S(v_r)`; `m̂(v_r) = mean_j ζ'_{j, S(v_r + ζ'_j)}` over inner draws that produce a winner (draws with no winner are skipped, as in the gate's own bias_sel convention); `Λ*_r = ζ*_{r, G_r} − m̂(v_r)`. Outer draws with no winner are skipped and their count recorded.
4. Outputs, added to the return list under `field`: `lambda_mean`, `lambda_sd`, quantiles q05/q25/q50/q75/q95 and q025/q975 of Λ*, `n_out_used`, `n_in_used_mean`; the second-order point estimate `est2 = β̃ − lambda_mean`; the one-sided 95% lower bound `lower_1s = β̃ − q95` (primary); the two-sided 95% interval `[β̃ − q975, β̃ − q025]`; and `se_field = lambda_sd` with the SE-type interval `est2 ± 1.96·se_field` (supplementary).
5. The point estimate reported for MR (field) is `est2` (decision E3); the current MR row keeps β̃ and the IJ interval unchanged.

Cost: (R_out + R_in) products of B_effᵀ with an n-vector, plus R_out × R_in applications of S over K candidates — vectorized per outer draw as an R_in × K matrix argmax. Expected ≤ 0.1 s for K ≤ 12 and ~1–3 s for K ≈ 151 single-core. The addendum-A M_eff Monte Carlo is not run in this task (decision E5); its values are already in the 2026-09-04 bundles and are joined by seed.

## Stage 0 — Discovery (no compute, no R/ edits)

0a. In `R/fs_mr_inference.R` (and whatever it calls), locate and quote: where β̂ over the kept family and B_eff are available after the fits; the re-selection routine applied per draw and whether it is callable on an arbitrary K-vector with the configured admission rule; the existing `multiplier = "gaussian"` path; how the gate seeds its draws; the return-list assembly (the D6 pattern from 736088e3 is the model).
0b. Confirm the Stage 1 harness `quarto/GuoHe/mr_vs_guohe_sim.R` and driver `mr_vs_guohe_run.R` regenerate replicates from the stored seeds bit-identically (they did on 2026-09-04) and identify the minimal change to record the new fields.
0c. State the K, n and orientation handling per lineage as before (treat-flip for §5.1; raw for §5.2).

Gate 0: all of 0a–0c quoted; any need to change existing arithmetic to reach the field quantities → STOP and report.

Output: `REPORT_mr_field_stage0_<date>.md`.

## Stage 1 — Implementation, identities, smoke test, projection

1a. Implement `ci_method = "field"` as specified, add-only; regenerate Rd; run the suite. Byte-identity of the default return object on fixed seeds is a Gate 1 check.

1b. Identities (machine-checked, tolerances stated):
- K = 1: `lambda_mean` within Monte Carlo error of 0, `lambda_sd`/σ̂_{D,Ĥ} in [0.95, 1.05], `est2` = β̃ within Monte Carlo error, the one-sided margin q95 within 5% of 1.645·σ̂_{D,Ĥ}.
- Covariance of the Gaussian-multiplier perturbations equals B_effᵀB_eff to within Monte Carlo error on a small family (report the max relative deviation over entries at R = 20,000).
- Seed reproducibility: two calls with the same seed give identical `field` output.
- Default-path byte-identity: `ci_method = "ij"` return object identical to the pre-change object on three fixed-seed cases (one per lineage plus K = 1).
- Sign of the correction: on an exchangeable disjoint tie family (K = 10, ng = 200, simulated), `lambda_mean` > 0 and `est2` closer to θ_{k̂} than β̃ on average over 200 replicates (reference from the proof of concept: retained bias 0.063 → 0.025 at K = 10; 0.018 → 0.004 at K = 2).

1c. Smoke test: 5 paired replicates each of one t35 cell, one t6 cell and one t7 cell; naive column `identical()` to stored; MR (current) est/se_ij identical to the 2026-09-04 bundle values on the same seeds; MR (field) finite everywhere.

1d. Projection, measured on a loaded machine: time the per-replicate cost with ≥ 30 concurrent workers, not single-core (the 2026-09-04 projection missed by 9× for this reason); report per-cell and full-grid core-hours and wall at the intended worker count.

Gate 1: 1b identities pass; 1c holds; suite green; projection reported. STOP for Larry's compute go.

Output: `REPORT_mr_field_stage1_<date>.md`.

## Stage 2 — Full run (after the compute go)

16 cells × 2,000 replicates, identical seeds and settings to 2026-09-04 (B = 5,000 centred Poisson for the current MR; R_out/R_in per E2 for the field). Per-replicate record: everything the 2026-09-04 record carried except the M_eff columns (joined from the old bundles by seed) plus the `field` fields, the two coverage indicators for MR (field) (one-sided at q95; two-sided quantile), and timings. New bundle names `mr_field_vs_guohe_<id>.rds`.

Gate 2: completeness 2,000/2,000 per cell; naive `identical()` to stored throughout; MR (current) est/se_ij identical to the 2026-09-04 bundles on every replicate (the pairing proof for the new run); MR (field) NA count zero or documented.

## Stage 3 — Report

`REPORT_mr_field_vs_guohe_<date>.md` plus `mr_field_vs_guohe.qmd` rendered. Same structure as the 2026-09-04 report, rows naive / G&H (four r) / MR (IJ, current) / MR (field, new); one-sided coverage and margin primary, two-sided for both MR rows supplementary; the tie table gains "retained after second-order" for `est2`; the diagnostics section reuses the joined M_eff/p̂/(A6) columns; the cost table adds the field pass. Findings in the record; no task proposed unless something blocks.

## Decisions required from Larry (defaults in brackets)

- E1 Field shrinkage: winner's entry → β̃ [default]; plug-in β̂ unshrunk as a second field [default: no].
- E2 Draws: R_out = 1,000, R_in = 500, Gaussian multipliers [default].
- E3 MR (field) point estimate: `est2 = β̃ − lambda_mean` [default]; alternative: keep β̃ and change only the interval.
- E4 Primary comparison: one-sided 95% lower bound, G&H convention [default]; two-sided quantile interval supplementary.
- E5 M_eff pass: off, joined from the 2026-09-04 bundles [default].
- E6 Compute: go/no-go on the Stage 1 loaded-machine projection.

## Done means

Stage 3 report and rendered qmd committed, Gate 2 record beside them, every Stage 1 identity listed with its measured value, branch left unpushed for Larry.
