# TASK — Complement inference refinements: the winner-only IJ variance and simultaneous (Ĥ lower, Ĥᶜ upper) bounds

Date: 2026-09-06. Author: chat (spec). Executor: Claude Code. Approver: Larry.
Sequencing: after `TASK_mr_field_complement_2026-09-06.md` has reached Stage 3 (campaigns `s7c`/`map1c` are this task's identity anchors). Fresh session.
Predecessors: the complement field block (`field_complement`), the field block, `fs_sim_bias_coverage(side=)`, the combined template. Proof-of-concept reference for the winner-only IJ: `POC_mr_interval_alternatives_2026-09-05.md` and `poc_ci_results_2026-09-05.csv` in `dev/tasks/` ("bag" and "bagfloor" constructions: exact when one candidate dominates, under-covering at ties on the harm side).

## Protocol

- First action: copy this file to `dev/tasks/` and commit. Do not push.
- Gates stop on failure; Gate 1 is the compute go/no-go (may be pre-authorized).
- Verify from source: the IJ residual construction for both blocks (`r` and `r_c`, the same-draws terms `P[sel, ]` / `Pc[sel, ]`), the finite-B correction, the complement field block's per-draw `Λ*ᶜ` and its alignment with the harm field's `Λ*` by outer-draw index.
- Committed bundles and documents are read-only; new runs write new campaign stems.
- No change to the working model, the selection map, the point estimates, or any existing default output.

## Purpose

Two additions for the complement, where the benefit claim lives, both free of new draws:

1. **Winner-only IJ variance** (item 1). The current IJ residual `r_b = (bias_sel + bias_fix) − D_{Ĥ*_b}(b) − D_Ĥ(b)` counts the same-draws term twice in the variance, which tends to 4σ² when one candidate dominates. The complement is effectively always in that regime (its identity barely changes across re-selection), so its IJ SE is 1.8–1.9× the truth in every cell. The winner-only residual `r_b^w = bias_sel − D_{Ĥ*_b}(b)` (complement: `sbc − selb_c`) is the total-derivative IJ and is exact when one candidate dominates. It costs nothing — same draws, one term dropped — and applies to the plain `"ij"` path, so it benefits users who never enable the field. On the harm side the proof of concept showed it under-covers at ties unless floored at the naive variance; both variants are recorded for both blocks so the table shows where each holds.
2. **Simultaneous bounds for the pair** (item 3). A development claim states both a harm lower bound on Ĥ and a benefit upper bound on Ĥᶜ. The field simulation yields the joint draws `(Λ*_r, Λ*ᶜ_r)` on the same multipliers and winners, so a joint 95% guarantee for the pair is a quantile of that joint distribution rather than two separate 95% bounds (whose joint coverage is below nominal) or a Bonferroni split (conservative).

## Method (add-only)

**A. `ij_residual = c("two_term", "winner", "winner_floor")`** on `fs_mr_inference()`, default `"two_term"` (current). Under `"winner"` the IJ variance for each block is built from the winner-only residual; under `"winner_floor"` it is `max(V_winner, σ̂²_naive)` for that block. The de-biased point estimates are unchanged. To keep the default output byte-identical and to evaluate all variants in one run, the gate always returns, beside the reported `se_ij`, the alternative SEs as `se_ij_winner` and `se_ij_winner_floor` (harm and complement) with their lower/upper bounds — additional list elements only; `ij_residual` selects which one populates the reported `debiased$lower/upper/se_ij`. Roxygen states the regime result: exact when one candidate dominates, under-covering at ties on the harm side unless floored.

**B. `field$joint`** (requires `ci_method = "field"`, `include_complement = TRUE`, `field_complement = TRUE`): from the aligned outer draws `(Λ*_r, Λ*ᶜ_r)` over the draws where both exist, find the equal-tail level γ such that `P*(Λ*_r ≤ q_{1−γ}(Λ*) and Λ*ᶜ_r ≥ q_γ(Λ*ᶜ)) ≥ 1 − α` (grid on γ from α down to α/2, finest step 0.001; report the achieved joint probability). Return `gamma`, the joint pair `lower_H = β̃ − q_{1−γ}(Λ*)` and `upper_Hc = β̃ᶜ − q_γ(Λ*ᶜ)` on the effect scale, the Bonferroni pair (γ = α/2) for reference, the correlation of `(Λ*, Λ*ᶜ)`, and `n_joint_draws`. No new draws; the marginal one-sided bounds are unchanged.

Classification: A changes the method (an alternative variance; a proposal for evaluation), add-only with the default byte-identical; B adds a derived quantity from existing draws, add-only. Forwarding lines in `forestsearch()` for `ij_residual`; template knob `FS_S7_IJ_RESIDUAL` (default `two_term`); recorder columns `mr_H_se_w`, `mr_H_lo_w/hi_w`, `mr_H_se_wf`, `mr_H_lo_wf/hi_wf` and the `Hc` counterparts; `fld_joint_gamma`, `fld_joint_loH`, `fld_joint_upHc`, `fld_joint_bonf_loH`, `fld_joint_bonf_upHc`, `fld_joint_corr`, `fld_joint_n`; coverage indicators: per-block one-sided (lower for Ĥ, upper for Ĥᶜ) and two-sided for the two IJ variants; joint coverage `β(Ĥ) ≥ lower_H and β(Ĥᶜ) ≤ upper_Hc` for the calibrated pair, the Bonferroni pair, and the two separate 95% field bounds.

## Stage 0 — Discovery (no compute, no R/ edits)

0a. Quote the IJ residual lines for both blocks, the finite-B correction and the fallbacks (`ij_source`), and where `se_ij`/`sec` feed the reported bounds.
0b. Quote the harm field's `lam` vector and the complement field's `Λ*ᶜ` vector and confirm they are indexed by the same outer draw `r` (the complement re-reads `G_out`), so the joint draws are aligned without re-drawing; quote the dropped-draw handling on each side.
0c. Cost anchors from the s7c/map1c Gate 2 records.

Gate 0: 0a–0b quoted; STOP if the joint draws cannot be aligned without changing either block.

## Stage 1 — Implementation, identities, smoke, projection

1a. Implement A and B, forwarding, template, recorder, tables (the complement block gains rows "MR (IJ, winner)" and "MR (IJ, winner-floor)"; a joint-bounds table with the three pairs' joint coverage), Rd, tests.
1b. Identities (machine-checked):
- Default-path byte-identity: `"ij"` and `"field"` (with and without complement) identical to pre-change objects on the three Guo–He fixed-seed cases and 5 replicates of h100/h175 at the s7c seeds — the new list elements excluded, everything else exact.
- K = 1: `se_ij_winner` equals the naive SE within Monte Carlo error on both blocks, and `se_ij` (two-term) ≈ 2× it (the 4σ² identity); `winner_floor` equals the naive SE exactly.
- Exchangeable K = 10 tie (simulated): harm-side `se_ij_winner` below the empirical SD (the known under-coverage regime), `winner_floor` ≥ naive SE.
- Joint: with the complement's field switched on, `gamma ∈ [α/2, α]`; the achieved joint probability ≥ 1 − α; when `(Λ*, Λ*ᶜ)` are independent in a constructed case, γ ≈ 1 − √(1 − α) ≈ 0.0253.
- Smoke: 5 replicates per cell, all existing columns identical to the s7c/map1c bundles.
1c. Projection at 100 workers (gate cost unchanged; the additions are arithmetic on existing draws).

Gate 1: identities pass; projection reported. Compute go per K-4.

## Stage 2 — Runs

Seven cells, `sim_id` 1–2,000, seeds as s7/map1, `ci_method = "field"`, `include_complement = TRUE`, `field_complement = TRUE`, `ij_residual = "two_term"` (reported rows unchanged; variants recorded), campaigns `s7w` / `map1w`, priority order as the complement task. Gate 2 per cell: completeness; all existing columns identical to the s7c/map1c bundles; new columns finite; interval invariants; `gamma` in range.

## Stage 3 — Report

`REPORT_complement_refinements_<date>.md` plus rendered documents. Per cell, complement block: rows naive / MR (IJ, two-term) / MR (IJ, winner) / MR (IJ, winner-floor) / MR (field): SE/SD, one-sided **upper** coverage with Wilson intervals, two-sided coverage, margin — and the same rows for the harm block with the lower bound. Joint table per cell: joint coverage and the two margins for (i) separate 95% field bounds, (ii) Bonferroni, (iii) calibrated γ; the mean γ and the mean correlation. The complement display (`block = "Hc"`, `side = "upper"`) gains the winner-only rows. Reading criteria (Larry's, not gates): complement `winner` SE/SD in [0.9, 1.1] with one-sided upper coverage within Monte Carlo error of nominal in every cell; harm-side `winner` under-covering at ties as expected with `winner_floor` at or above nominal; joint calibrated coverage within Monte Carlo error of 0.95 with margins below Bonferroni's. Findings in the record; no task proposed unless something blocks.

## Decisions (defaults in brackets)

- K-1 Variants recorded: `winner` and `winner_floor` for both blocks [default]; reported row stays two-term for this evaluation.
- K-2 Joint convention: equal-tail γ from the joint draws, Bonferroni as reference [default]; alternative — unequal tails weighted toward the complement.
- K-3 Cells: all seven [default].
- K-4 Compute: pre-authorization with a ceiling (expected ~3 h at 100 workers, as the complement run).
- K-5 Defaults after Stage 3: whether `ij_residual = "winner_floor"` becomes the complement's reported SE, and whether the joint pair is added to the standard output — Larry's call on the record.

## Done means

Stage 3 report and rendered documents committed; Gate 2 records beside them; `ij_residual` and `field$joint` landed add-only with default-path identities recorded; template updated under its own stems; branch left unpushed for Larry.
