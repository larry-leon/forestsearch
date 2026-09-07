# TASK — Prevalence ~30%: the current constructions under a larger harm subgroup

Date: 2026-09-06. Author: chat (spec). Executor: Claude Code. Approver: Larry.
Predecessors: the s7/map1 campaigns (12.5% prevalence, M1 design), the complement field (`field_complement`), the joint pair (`field$joint`), the bias-coverage display; combined template `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` in `quarto/simulations/gbsg_020/`.

## Protocol

- First action: copy this file to `dev/tasks/` and commit. Do not push.
- Gates stop on failure; Gate 1 is the compute go/no-go (may be pre-authorized).
- No change under `R/`. The prevalence change is a document-level DGM knob on the template; the committed template's defaults must reproduce the s7 columns exactly (identity in Stage 1).
- Verify from source: how the template's DGM defines the harm subgroup and its prevalence, how `calibrate_k_inter`, the CDE/marginal truths and `fs_attach_betaHhat` depend on that definition, and what the harm rule's form is (so the ~30% rule keeps the same form).
- Committed bundles and documents read-only; new runs under campaign `p30`.
- Before Stage 2, the batch save must refuse to write to a git-tracked path (the guard from the cell-7 incident): add it to the template if absent — add-only, document-level.

## Purpose

Every complement result so far comes from designs where the harm subgroup is ~12.5% of the trial and near-tied candidates have near-identical complements, so the complement's selection variability is negligible. That is a regime, not a property of any method: it is why the naive SE happened to be adequate for the complement, and it is the reason an unadjusted SE cannot be promoted. At ~30% prevalence the complement is smaller, competing candidates' complements differ more, and the complement's selection variability should surface. The question is whether the selection-adjusted constructions hold there and the unadjusted one does not.

Constructions evaluated (the current set; nothing else): naive and oracle as references; MR (IJ, two-term) — point estimate β̃, two-sided interval, one-sided bound; MR (field) — est₂, plain two-sided, one-sided bound (lower for Ĥ, upper for Ĥᶜ); the joint pair (separate / Bonferroni / calibrated). The winner-only and winner-floor IJ variants are excluded from every table, figure and report line (decision H-P3); their recorder columns may exist but are not summarised.

## Stage 0 — Discovery (no compute)

0a. Quote the DGM section of the template: the harm-subgroup rule, its covariates and cut(s), how prevalence ≈ 12.5% arises, the `target_hr_harm` calibration, and how the truth targets (θ†, θ‡, per-replicate β(Ĥ)/β(Ĥᶜ)) are attached — confirm they are computed from whatever rule the document defines, not hard-coded to M1.
0b. Propose a harm rule of the same form as M1's (single cut or conjunction on the same covariates) with super-population prevalence in [0.28, 0.32]; quote its realized prevalence on the 100k super-population, its calibrated `k_inter` at the H-P1 targets, and the corresponding θ†/θ‡ for both blocks. If the form must change, say so and stop for Larry's choice (decision H-P2).
0c. Knob: a `FS_S7_PREV` (or rule-selecting) knob whose default reproduces M1 exactly; quote the lines.
0d. Cost anchors from s7c/map1c.

Gate 0: 0a–0c quoted; STOP if the truth attachment does not follow the rule, or if a same-form rule cannot reach [0.28, 0.32].

Output: `REPORT_p30_stage0_<date>.md`.

## Stage 1 — Knob, identities, smoke, projection

1a. Add the knob and the save guard to the template (document-level, add-only).
1b. Identities: with the knob at default, 5 replicates of h100 at the s7 seeds identical (≤ 1e-12) to the s7c bundle on every column — the knob is inert at default; with the knob at ~30%, 5 replicates per H-P1 cell: realized prevalence per replicate in [0.24, 0.36], truth targets finite, detection non-zero, field and complement field finite, joint γ in range; the bound identities as before.
1c. Projection at 100 workers (expect the s7 pace; the complement fits are more numerous when complements differ).

Gate 1: identities pass; projection reported. Compute go per H-P5.

Output: `REPORT_p30_stage1_<date>.md`.

## Stage 2 — Runs

Cells (H-P1), `sim_id` 1–2,000, two batches, campaign `p30`, `ci_method = "field"`, `include_complement = TRUE`, `field_complement = TRUE`, `ij_residual = "two_term"`, FB none, 100 workers, fail-fast per cell:
1. prevalence 30%, HR 1.00, n = 500 (null; the tie regime, where the harm-side one-sided bound is most exposed)
2. prevalence 30%, HR 1.50, n = 500 (moderate harm)
3. prevalence 30%, HR 1.75, n = 500 (harm)
4. prevalence 30%, HR 1.00, n = 1,000 (the null at the n where the screen effect appeared)
Gate 2 per cell: completeness; realized prevalence; finite fields; invariants; γ in range. No identity anchor exists for these cells beyond the knob-inert check of Stage 1.

## Stage 3 — Report

`REPORT_p30_<date>.md` plus rendered documents. Per cell, both blocks, rows naive / oracle / MR (IJ) / MR (field): bias vs β(·) in log-HR and SD units, empirical SD, mean SE, SE/SD, two-sided coverage of β(·) with Wilson intervals, one-sided coverage on the exposed side (Ĥ lower, Ĥᶜ upper), log half-width and one-sided margin; the bound-location tables (share of Ĥ lower ≥ 0.85/0.95; share of Ĥᶜ upper < 0.85/0.80); the joint table; the display for both blocks; regime diagnostics: p̂(Ĥ), M_eff if recorded, and for the complement the ratios SD(β̃ᶜ)/naive SE and λ-SDᶜ/naive SE — the two numbers that say whether the complement's selection variability is negligible. Findings in the record; no task proposed unless something blocks.

Reading criteria (Larry's, not gates), read by location and coverage, never as significance at 1.0:
- Complement: SD(β̃ᶜ)/naive SE materially above 1 (the regime has changed); the field's one-sided upper coverage of β(Ĥᶜ) within Monte Carlo error of nominal with λ-SDᶜ tracking the SD; IJ two-term's SE/SD reported as found.
- Harm: the field's one-sided lower coverage as at 12.5% (0.92–0.98); the two-sided dip reported as found.
- Joint: Bonferroni within Monte Carlo error of 0.95; the correlation of the draws reported (it need not be near zero here).

## Decisions (defaults in brackets)

- H-P1 Cells: the four above [default].
- H-P2 Harm rule: same form as M1, prevalence in [0.28, 0.32] [default]; a changed form only on Larry's choice at Gate 0.
- H-P3 Winner variants: excluded from all tables, figures and reports [default, per Larry].
- H-P4 Workers: 100 [default].
- H-P5 Compute: go at Gate 1, or an unattended pre-authorization with a wall ceiling (expected ~2.5 h at the s7 pace).

## Done means

Stage 3 report and rendered documents committed; Gate 2 records beside them; the knob and save guard landed on the template with the knob-inert identity recorded; branch left unpushed for Larry.
