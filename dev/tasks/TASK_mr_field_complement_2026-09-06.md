# TASK — Field block for the complement Ĥᶜ: the one-sided upper bound for a benefit claim

Date: 2026-09-06. Author: chat (spec). Executor: Claude Code. Approver: Larry.
Predecessors: `TASK_mr_field_section7_2026-09-05.md` (field block for Ĥ; s7 bundles; combined template), `TASK_mr_field_ocmap_2026-09-05.md` (map1 bundles), `TASK_bias_coverage_display_2026-09-06.md` (`fs_sim_bias_coverage()`, `fs_plot_bias_coverage()`). Complement-interval fix a6702fd8 (complement IJ interval under `"field"`).

## Protocol

- First action: copy this file to `dev/tasks/` and commit. Do not push.
- Gates stop on failure; Gate 1 is the compute go/no-go (may be pre-authorized).
- Verify from source: the gate's existing complement path (`include_complement`, the complement fits for draw winners, `complement$debiased`, `bias_sel_c`/`bias_fix_c`), the field block's outer/inner draw structure and seeds, and the template's complement recorder columns, before writing anything.
- Committed bundles and documents are read-only; new runs write new campaign stems.

## Purpose

Everything the field delivers for the harm subgroup — a calibrated one-sided bound, est₂, the plain two-sided interval, diagnostics — is missing for the complement, which carries only the IJ interval (SE/SD 1.8–1.9, coverage 1.000 in every cell) and is, for development purposes, the subgroup that matters: it is where treatment is expected to work. This task adds the complement's field block, evaluated on the same replicates as everything else.

Orientation. Selection acts on the harm-subgroup effects; the complement is its by-product, and its naive estimate is pushed toward *benefit* (−0.6 SD across the cells; the corrected estimate is within 0.1 SD). For a benefit claim ("HR in Ĥᶜ is at most U") the exposed limit is therefore the **upper** bound on the complement's log-HR: under-correction of the beneficial optimism pushes U too low and over-claims. The mirror of the harm-side result applies — an upper bound is vulnerable to under-correction and immune to over-correction. The primary complement product is the one-sided 95% upper bound; the lower bound and the two-sided interval are also returned.

## Method (add-only)

New argument `field_complement = FALSE` on `fs_mr_inference()`; active only when `ci_method = "field"` and `include_complement = TRUE`. Default FALSE reproduces every existing return byte-for-byte (the s7/map1 bundles were produced with `"field"` + `include_complement = TRUE`, so the default must not add computation to that path).

When TRUE, after the harm field block:

1. Selection is unchanged: each outer draw r uses the harm field's `v_r = w + ζ*_r` and winner `G_r = S(v_r)`; each inner draw likewise. The complement enters only through its own effects and influences.
2. Complement effects and influences `β̂ᶜ(g)`, `B_effᶜ[, g]` for every candidate g that wins on any multiplier, outer, or inner draw: fit lazily with a cache keyed by candidate (the gate already fits complements for the multiplier-draw winners and the selected candidate; extend the same routine). Record the number of distinct complement fits and the share of draws that required a new fit.
3. Complement perturbations from the **same** multipliers as the harm field: `ζᶜ_r = B_effᶜᵀ ξ_r` with the ξ_r already drawn (so the complement's noise carries its correct correlation with every candidate's noise); same for the inner draws.
4. `Λ*ᶜ_r = ζᶜ_{r, G_r} − m̂ᶜ(v_r)`, with `m̂ᶜ(v) = mean over inner draws of ζᶜ'_{G(v+ζ')}` (winners from the harm field, complement noise read for those winners); draws without a winner dropped, as in the harm block.
5. Returned under `field$complement`: `est2` = β̃ᶜ − mean(Λ*ᶜ); `upper_1s` = β̃ᶜ − q₀.₀₅(Λ*ᶜ) (primary); `lower_1s` = β̃ᶜ − q₀.₉₅(Λ*ᶜ); `lower_2s`/`upper_2s` = β̃ᶜ − q₀.₉₇₅ / β̃ᶜ − q₀.₀₂₅; `se_field` = sd(Λ*ᶜ); `lambda_mean`; the seven quantiles; `n_out_used`, `n_in_used_mean`, `n_complement_fits`, `share_draws_new_fit`, `timing_seconds`; all effect-scale quantities exponentiated as the harm block does. β̃ᶜ itself is the existing two-term complement estimate, unchanged.

Classification: changes the method (a complement field block; a proposal requested for evaluation); add-only; default byte-identical; suite unchanged; the harm field block's draws and outputs untouched. The forwarding line in `forestsearch()` for `field_complement` is add-only pass-through.

Display: `fs_sim_bias_coverage()` gains `side = c("lower", "upper")` (default `"lower"`, current behaviour); under `"upper"` the one-sided coverage is `truth ≤ upper bound` and the Gaussian reference is Φ(1.645·r + b). The `block = "Hc"` path stops dropping the field estimator when the `fld_Hc_*` columns are present. Add-only.

## Stage 0 — Discovery (no compute, no R/ edits)

0a. Quote the gate's complement path: where complement fits happen for draw winners, the complement influence matrix, `complement$debiased` and its bias terms; the field block's outer/inner structure (`ξ` draws, winners, `m̂`) and seed offsets; the return assembly.
0b. Quote the template's complement recorder columns (`mr_Hc_est/lo/hi/se_ij`, `nv_Hc_*`, `or_Hc_*`, `betaHhat_Hc`, `cde_Hc`, `marg_Hc`) and the coverage helpers, and the `fs_sim_bias_coverage()` `block = "Hc"` branch.
0c. Cost anchors: complement fit time per candidate; the s7 per-replicate gate cost (~28–40 s under load); the field block cost (~15 s).

Gate 0: 0a–0c quoted; STOP if the complement field cannot reuse the harm field's ξ draws and winners without changing them.

Output: `REPORT_mr_field_complement_stage0_<date>.md`.

## Stage 1 — Implementation, identities, smoke, projection

1a. Implement per the method; roxygen (state the orientation and the guarantee mirror); Rd; suite. Template: knob `FS_S7_FIELD_COMPLEMENT` (default FALSE), recorder columns `fld_Hc_est2`, `fld_Hc_up1s` (primary), `fld_Hc_lo1s`, `fld_Hc_lo2s`, `fld_Hc_hi2s`, `fld_Hc_lo_se`, `fld_Hc_hi_se`, `fld_Hc_se`, `fld_Hc_lam_mean`, the seven quantiles, `fld_Hc_nout`, `fld_Hc_nin_mean`, `fld_Hc_nfit`, `fld_Hc_share_newfit`, `fld_Hc_secs`, `fld_Hc_note`; coverage indicators for Ĥᶜ: two-sided, and one-sided **upper** (β(Ĥᶜ) ≤ U); tables gain the MR (field) row in the complement block; the diagnostics section adds the complement's λ-SD/SD and retained bias; the display chunk calls `fs_sim_bias_coverage(block = "Hc", side = "upper")`.
1b. Identities (machine-checked):
- Default-path byte-identity: `"ij"`, and `"field"` with `field_complement = FALSE` (with and without `include_complement`), identical to pre-change objects on the three fixed-seed Guo–He cases and 5 s7 replicates per cell; the harm `field` block identical with `field_complement = TRUE` (its draws are shared, not redrawn).
- K = 1 with complement: the complement's Λ*ᶜ has mean 0 within Monte Carlo error and sd equal to the complement's naive SE within 5%; `upper_1s` equals the naive one-sided upper bound within Monte Carlo error.
- Coupling check on one replicate: the correlation across outer draws between the winner's harm perturbation and the complement's perturbation equals the corresponding entry of `B_effᵀB_effᶜ` within Monte Carlo error.
- Bound identities: `upper_1s = exp(β̃ᶜ − q05)`, `lower_2s = exp(β̃ᶜ − q975)`, `hi_se = exp(log est2 + 1.96·se)` at ≤ 1e-12.
- Smoke: 5 replicates per cell (h100, h175) at the committed seeds — all non-complement-field columns identical to the s7 bundles; complement field finite; `n_complement_fits` and `share_draws_new_fit` reported.
1c. Projection at 100 workers under load for the H-C1 cells.

Gate 1: identities pass; projection reported. Compute go per H-C5.

Output: `REPORT_mr_field_complement_stage1_<date>.md`.

## Stage 2 — Runs (after the compute go)

Cells in priority order (H-C1), `sim_id` 1–2,000 each, same seeds as s7/map1, `field_complement = TRUE`, FB none, campaign `s7c` (s7 cells) / `map1c` (map1 cells): h100 n500, h175 n500, h150 n500, h150 n1500, then h075, h100 n1000, h175 knoise3 if the ceiling allows. Gate 2 per cell: completeness; all non-complement-field columns identical (≤ 1e-12) to the s7/map1 bundles (the pairing proof); complement field NA count zero or documented; interval invariants.

## Stage 3 — Report

`REPORT_mr_field_complement_<date>.md` plus rendered documents. Per cell, complement block, rows naive / oracle / MR (IJ) / MR (field): bias vs β(Ĥᶜ) in log-HR and SD units, empirical SD, mean SE, SE/SD, two-sided coverage of β(Ĥᶜ) with Wilson intervals, one-sided **upper** coverage, log half-width and one-sided margin; coverage of θ†(Ĥᶜ) and θ‡(Ĥᶜ) reported, not scored; the complement display (block Hc, side upper) across cells; complement-fit counts and cost. Reading criteria (Larry's, not gates): field one-sided upper coverage of β(Ĥᶜ) within Monte Carlo error of nominal in every cell; λ-SD/SD in [0.9, 1.2]; margin materially below IJ's (expected roughly half); est₂ᶜ bias not worse than β̃ᶜ's beyond Monte Carlo error. Findings in the record; no task proposed unless something blocks.

## Decisions (defaults in brackets)

- H-C1 Cells and order: as listed [default]; the two s7 cells only if cost demands.
- H-C2 Complement fits: lazy with cache over all draw winners [default]; alternative — precompute for the mass-carrying set only, with draws outside it dropped and counted.
- H-C3 Primary complement bound: one-sided upper (benefit) [default]; lower and two-sided returned.
- H-C4 Workers: 100 [default].
- H-C5 Compute: go at Gate 1, or an unattended pre-authorization with a wall ceiling (expected 35–50 min per cell at 100 workers).
- H-C6 Analysis document: after Stage 3, the GBSG frozen-family document gains the complement rows in its interval table — a separate small task.

## Done means

Stage 3 report and rendered documents committed; Gate 2 records beside them; `field_complement` landed add-only with default-path identities recorded; `fs_sim_bias_coverage(side=)` landed add-only; template updated under its own stems; branch left unpushed for Larry.
