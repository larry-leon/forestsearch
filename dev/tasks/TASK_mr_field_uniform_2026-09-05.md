# TASK — Field interval, uniform calibration: the κ(Σ̂) two-sided interval and the one-sided bound as the standard product

Date: 2026-09-05. Author: chat (spec). Executor: Claude Code. Approver: Larry.
Supersedes `TASK_mr_field_adaptive_2026-09-05.md` (drafted, never kicked off; withdrawn — the limit-experiment analysis showed adaptive shrinkage moves the two-sided worst case from 0.896 to 0.912 and is not the mechanism).
Predecessors: `TASK_mr_field_section7_2026-09-05.md` (s7 campaign bundles, template `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` in `quarto/simulations/gbsg_020/`), `TASK_mr_field_vs_guohe_2026-09-05.md` (harness `quarto/GuoHe/mr_vs_guohe_sim.R`, field bundles). Reference computations: `limit_sweep.py`, `limit_sweep_K.py`, `limit_sweep_K2_rules_2026-09-05.csv`, `limit_sweep_K10_kappa_2026-09-05.csv` (chat outputs, 2026-09-05; to travel with this document into `dev/tasks/`).

## Protocol

- First action: copy this file and the four reference files to `dev/tasks/` and commit. Do not push.
- Gates stop on failure; Gate 1 is the compute go/no-go (may be pre-authorized).
- Verify from source; quote paths, signatures, line numbers.
- Reports beside the s7 records; the Guo–He check writes beside `quarto/GuoHe/REPORT_mr_field_vs_guohe_2026-09-05.md`.
- Committed bundles and documents are read-only inputs; new runs write new campaign stems.

## Purpose and method

The field interval's one-sided bound is uniformly at nominal over the local parameter (limit experiment: 0.93–0.957 across separations 0–4 SE, K = 2 and K = 10, disjoint and overlapping; empirically 0.933–0.966 in every cell of both arcs). Its two-sided quantile interval is not: worst case ≈ 0.90 at 2–2.5 SE separation for small overlapping families, down to 0.77 for many disjoint competitors, and no shrinkage rule repairs it (structural: a separated truth can look tied in the data). This task (i) makes the one-sided bound the field method's standard product and (ii) adds a two-sided interval whose validity is uniform over a stated protection family, obtained by widening the field quantiles by a factor κ computed at analysis time from the trial's own Σ̂ — nothing tuned to any simulation design.

**The κ computation (per analysis).** Inputs already inside the gate after the fits and the multiplier pass: the candidate effects β̂ (K), the influence matrix B_eff (n × K, so Σ̂ = B_effᵀB_eff), the winner Ĥ, β̃(Ĥ), the re-selection frequencies p̂, and the selection routine S with its thresholds and sizes fixed.

1. Mass-carrying set: order candidates by p̂, keep the smallest M whose cumulative p̂ ≥ 0.99, with M capped (H3) and Ĥ always included; restrict β̂, B_eff, S to it. Record M and the mass covered.
2. Protection family (H1, "winner-profile"): for δ on a grid (H2), the hypothetical true field μ_δ = w (the gate's shrunk field) with the winner's entry set to max_{g ≠ Ĥ} w_g + δ·σ̂_Ĥ.
3. For each δ, R_rep hypothetical trials: ζ_r = B_effᵀξ_r, ξ_r ~ N(0, I_n) (Gaussian multipliers, covariance Σ̂ exactly); W_r = μ_δ + ζ_r; G_r = S(W_r), trials with no winner dropped (conditioning on detection, as in the real analysis); Λ_r = ζ_{r,G_r} − m̂(W_r), m̂ by R_in inner draws. Then the procedure as it runs on real data: the shrunk field from W_r, R_out field draws, Λ*_r, its quantiles (q₀.₀₂₅, q₀.₉₅, q₀.₉₇₅) and mean c_r.
4. Coverage profiles: C₁(δ) = P(Λ_r ≤ q₀.₉₅,r); C₂(δ; κ) = P(c_r + κ(q₀.₀₂₅,r − c_r) ≤ Λ_r ≤ c_r + κ(q₀.₉₇₅,r − c_r)).
5. κ* = the smallest κ ∈ [1, 2] (grid or root-find) with min_δ C₂(δ; κ) ≥ 1 − α; report κ*, min_δ C₁, the profiles, and the Monte Carlo SE of κ* (from the replication count).
6. Return the uniform two-sided interval [β̃ − c − κ*(q₀.₉₇₅ − c), β̃ − c − κ*(q₀.₀₂₅ − c)] for the real analysis (c = the real field's λ̄), alongside the plain quantile interval and the one-sided bound.

Cost: R_rep × R_out × R_in × M argmaxes per δ; at H4 defaults and M ≤ 12, seconds per δ on a prespecified family, minutes per analysis on an enumerated family after the top-M reduction. Stage 1 measures it.

## Method change (called out)

New internal function (working name `fs_mr_field_uniform()`) and an add-only `field_uniform = FALSE` argument on `fs_mr_inference()` (returns the κ block under `field$uniform` when TRUE). Changes the method: a proposal, requested for evaluation. Default output byte-identical (both `"ij"` and `"field"` with `field_uniform = FALSE`); the suite unchanged; the field's own draws and seeds untouched (the κ sweep uses a separately derived seed). Roxygen documents the protection family, the mass-carrying reduction, and the guarantee: one-sided uniformly valid; two-sided uniformly valid over the winner-profile family at the reported κ*; plain quantile interval approximate. Anything touching the IJ path, β̃, or S is a STOP.

## Stage 0 — Discovery (no compute, no R/ edits)

0a. Quote where B_eff, β̂, Ĥ, β̃, p̂ and the re-selection routine are in scope after the field block in `R/fs_mr_inference.R`; the field's seed derivation; the return-list assembly (the D6 pattern).
0b. Quote the template's `mr_inference_args` and recorder lines and the Guo–He harness call site for the pass-through of `field_uniform` and the new columns.
0c. Confirm the s7 bundles' seed scheme (`seed_base + sim_id`) regenerates replicates bit-identically (it did on 2026-09-05), so κ can be evaluated on the same replicates.
0d. Cost anchors: gate call per replicate under load from the s7 cost table (~39 s incl. search and field at 100 workers).

Gate 0: 0a–0d quoted; STOP if κ cannot be computed without altering the field's own draws or the re-selection routine.

Output: `REPORT_mr_field_uniform_stage0_<date>.md`.

## Stage 1 — Implementation, identities, smoke, projection

1a. Implement per the method; template and harness gain the pass-through and columns `fld_H_kappa`, `fld_H_M`, `fld_H_mass`, `fld_H_minC1`, `fld_H_lo2u`, `fld_H_hi2u`, `fld_H_kappa_mcse`, `fld_H_uniform_secs`.
1b. Identities (machine-checked, tolerances stated):
- Default-path byte-identity: `"ij"` and `"field"` (with `field_uniform = FALSE`) identical to the pre-change objects on the three fixed-seed Guo–He cases and 5 s7 replicates per cell.
- **K = 2, Σ = I, pure argmax (the closed-form reference, `limit_sweep_K2_rules_2026-09-05.csv`, rule A):** with R_rep ≥ 2,000, two-sided C₂(δ; 1) within 0.015 of 0.945/0.93/0.905/0.90/0.935 at δ = 0/1/2/2.5/4; one-sided C₁ in [0.93, 0.96] at every δ; κ* within 0.05 of 1.20; C₂(2.5; 1.2) within 0.015 of 0.947.
- **K = 10 exchangeable (`limit_sweep_K10_kappa_2026-09-05.csv`, R_rep = 800):** disjoint: C₂(δ; 1) within 0.025 of 0.966/0.935/0.840/0.766 at δ = 0/1/2/3, C₁ in [0.92, 0.96]; ρ = 0.7: C₂(δ; 1) within 0.025 of 0.948/0.931/0.919/0.936 and C₂(δ; 1.2) of 0.991/0.970/0.959/0.974.
- K = 1: C₁ and C₂ at nominal for every δ, κ* = 1 within Monte Carlo error.
- Top-M reduction: on a K = 12 family where two candidates carry ≥ 0.99 of p̂, κ* from the full family and from the reduced family agree within its Monte Carlo SE.
- Seed reproducibility: two calls with the same seed give identical uniform blocks; the field block itself is unchanged by `field_uniform = TRUE` (identical to `FALSE`).
1c. Smoke: 5 replicates per cell (h100, h175) at the committed seeds with `field_uniform = TRUE`: all non-uniform columns identical to the s7 bundles; κ*, M, mass, minC₁ finite; timing per replicate.
1d. Projection at 100 workers, measured under load, for Stage 2 as scoped by H5.

Gate 1: identities pass; projection reported. Compute go per H6.

Output: `REPORT_mr_field_uniform_stage1_<date>.md`.

## Stage 2 — Runs (after the compute go)

2a. **κ on real Σ̂:** h100 and h175, `sim_id` 1–2,000, campaign `s7u`, `field_uniform = TRUE` — one gate call per replicate yields MR (IJ), MR (field) and the uniform block. Gate 2a: completeness; all non-uniform columns identical to the s7 bundles; κ* finite on every detected replicate. Report the distribution of κ*, M, mass covered, minC₁ per cell.
2b. **Guo–He check set (H5):** t35 β₂ = 0.3 and 0.5 (the two-sided dip), t6 k = 10 (many disjoint), t7 β₂ = 0.0 (nested, top-M in action), campaign tag `uniform`. Gate 2b as 2a against the 2026-09-05 field bundles.

## Stage 3 — Report

`REPORT_mr_field_uniform_<date>.md` plus rendered documents. Per cell, harm block, rows MR (IJ) / MR (field, one-sided) / MR (field, plain two-sided) / MR (field, uniform two-sided): coverage of β(Ĥ) with Wilson intervals, retained bias, SE/SD, z·SE or half-width on the log scale; the κ* distribution per cell; runtime. Reading criteria (Larry's, not gates): uniform two-sided coverage of β(Ĥ) ≥ 0.93 in every cell including the t35 dip cells and h175; one-sided unchanged; κ* on forestsearch families reported as found. Findings in the record; no task proposed unless something blocks.

## Decisions required from Larry (defaults in brackets)

- H1 Protection family: winner-profile [default]; the two-parameter family (winner separation × number of tied competitors) as an option [default: no].
- H2 δ grid: 0 to 4 SE by 0.5 [default].
- H3 Mass-carrying set: cumulative p̂ ≥ 0.99, M ≤ 12 [default].
- H4 Draws: R_rep = 300, R_out = 300, R_in = 150 per δ [default]; Stage 1 reports the Monte Carlo SE of κ* at these counts.
- H5 Cells: h100, h175, and the four Guo–He check cells [default].
- H6 Compute: go at Gate 1, or an unattended pre-authorization with a wall ceiling.
- H7 Policy (after Stage 3): one-sided bound as the standard reported quantity; two-sided by κ(Σ̂) or by IJ — Larry's call on the record.

## Done means

Stage 3 report and rendered documents committed, Gate 2 records beside them, `field_uniform` landed add-only with default-path identities recorded, the four reference files in `dev/tasks/`, template and harness updated under their own stems, branch left unpushed for Larry.
