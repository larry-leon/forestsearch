# TASK — Identifier variant `effMaxSG` on the prevalence-30% cells: capture of the planted region, and the constructions under it

Date: 2026-09-06. Author: chat (spec). Executor: Claude Code. Approver: Larry.
Predecessors: `TASK_prevalence30_2026-09-06.md` (campaign `p30`, four cells, `sg_focus = "maxeffCons"`; records `REPORT_p30_*` beside the s7 records; combined template with `FS_S7_Z1Q`, the `.refuse_if_tracked()` guard); all constructions as evaluated there.

## Protocol

- First action: copy this file to `dev/tasks/` and commit. Do not push.
- Gates stop on failure; Gate 1 is the compute go/no-go (may be pre-authorized).
- No change under `R/`. The identifier switch is a document-level knob on the template whose default reproduces `p30` exactly (identity in Stage 1).
- Verify from source: what `sg_focus = "effMaxSG"` does inside `forestsearch()` (the selection functional: effect ordering, the size neighbourhood and its width, the consistency screen's role, tie-breaks), and what re-selection map `fs_mr_inference()` applies under it on every multiplier draw and in the field's outer/inner draws. Quote both.
- Committed bundles and documents read-only; new runs under campaign `p30sg`.

## Purpose

At 31% planted prevalence the `maxeffCons` identifier returned about half the region (|Ĥ| 80 of 500 against 153 true), so the complement stayed at ~85% of the trial and the regime did not move. The effect-argmax rewards the most extreme sub-region. `effMaxSG` — the largest subgroup within a neighbourhood of the maximal effect — should return a larger Ĥ, closer to the planted region, at some cost in the estimated effect and with less optimism to remove. Two questions, on the same replicates as `p30`:

1. **Identification.** How much of the planted region does `effMaxSG` capture (sensitivity, specificity, |Ĥ|, detection), against `maxeffCons` on identical data.
2. **The constructions under this identifier.** Naive / oracle / MR (IJ two-term) / MR (IJ winner-only) / MR (field), both blocks, the joint pair, and the regime diagnostics — with the complement now smaller if question 1 comes out as expected. The winner-floor variant is excluded from every table, figure and report line.

## Alignment (the stop condition)

The correction is faithful only if the rule re-applied on each draw is the rule that selected Ĥ. Stage 0 must establish that the gate's re-selection under `effMaxSG` (and the field's `sel_one`) implements the same functional — effect ordering, the same neighbourhood width, the same size tie-break, the consistency thresholds held fixed — as the identifier. If it does not, STOP at Gate 0 and report the gap; do not run a misaligned correction. If alignment requires an R/ change, that is a separate proposal for Larry, not part of this task.

## Stage 0 — Discovery (no compute)

0a. Quote the `effMaxSG` selection code in `forestsearch()` and the neighbourhood parameter (name, default, where it is set in the template), and the gate's re-selection code path for this focus (`reselection`, `selection_rule`, `nbhd_pct` or equivalent) including the field's outer/inner `sel_one`. State whether they are the same functional.
0b. Quote the template's `sg_focus` / `focus_tag` lines; propose `FS_S7_FOCUS` (default `maxeffCons`) with the stem following `focus_tag` as it already does.
0c. Cost anchors: the p30 per-replicate costs; whether the neighbourhood rule forces the field's slow selection path (it already does under `maxeffCons` thresholds).

Gate 0: alignment established and quoted; knob lines quoted; STOP on misalignment.

Output: `REPORT_p30sg_stage0_<date>.md`.

## Stage 1 — Knob, identities, smoke, projection

1a. Add `FS_S7_FOCUS` (document-level, add-only).
1b. Identities: knob at default — 5 replicates of the p30 HR 1.00 n500 cell at the p30 seeds identical (≤ 1e-12) to the `p30` bundle on every column; knob at `effMaxSG` — 5 replicates per cell: selection differs from `maxeffCons` on at least some replicates (otherwise the knob is not reaching the identifier), fields finite, the gate's re-selection frequencies sum correctly, bound identities, γ in range; on one replicate, verify that the observed Ĥ is reproduced by applying the gate's re-selection map to the unperturbed effect vector (the alignment check in numbers).
1c. Projection at 100 workers.

Gate 1: identities pass; projection reported. Compute go per Q5.

Output: `REPORT_p30sg_stage1_<date>.md`.

## Stage 2 — Runs

The four p30 cells (Q1), `sim_id` 1–2,000, seeds `8316951 + sim_id`, `FS_S7_Z1Q=0.60`, `FS_S7_FOCUS=effMaxSG`, `ci_method = "field"`, `include_complement = TRUE`, `field_complement = TRUE`, `ij_residual = "two_term"`, FB none, 100 workers, fail-fast per cell, campaign `p30sg`. Gate 2 per cell: completeness; realized prevalence (same draws as p30 — identical by seed; assert); all fields finite; invariants; γ in range; sens/spec/|Ĥ| recorded.

## Stage 3 — Report

`REPORT_p30sg_<date>.md` plus rendered documents.

1. **Identification, paired by replicate (p30 vs p30sg):** detection, mean |Ĥ|, sensitivity, specificity, PPV, NPV; the distribution of |Ĥ|/|H|; the share of replicates where Ĥ under `effMaxSG` strictly contains Ĥ under `maxeffCons`; β(Ĥ) and β(Ĥᶜ) means against θ†(H) and θ†(Hᶜ) (how much closer the found region's effect is to the planted one).
2. **Constructions under `effMaxSG`:** per cell, both blocks, rows naive / oracle / MR (IJ two-term) / MR (IJ winner-only) / MR (field): bias in log-HR and SD units, empirical SD, mean SE, SE/SD, two-sided coverage of β(·) with Wilson intervals, one-sided coverage on the exposed side, half-widths; the bound-location tables (Ĥ lower ≥ 0.85/0.95; Ĥᶜ upper < 0.85/0.80); the joint table; the display for both blocks.
3. **Regime diagnostics:** p̂(Ĥ), complement fits per replicate, corr(Λ*, Λ*ᶜ), and for the complement SD(β̃ᶜ)/naive SE and λ-SDᶜ/naive SE — beside the same numbers from p30.
4. **Naive optimism under each identifier:** naive bias in SD units, both blocks, p30 vs p30sg — whether the larger region carries less winner's curse.

Reading criteria (Larry's, not gates; by location and coverage, never significance at 1.0): sensitivity materially above p30's (≈ 0.5); |Ĥ| approaching 153 / 306; the complement's SD(β̃ᶜ)/naive SE reported as found — if it rises above ~1.05 the regime has moved and the winner-only and field rows on the complement become the test they were meant to be; the field's one-sided bounds within 0.92–0.98 on both blocks; whether the complement's upper bounds in the harm cells fall back below 0.85/0.80 once more of the harm region is captured.

## Decisions (defaults in brackets)

- Q1 Cells: the four p30 cells [default]; the 12.5% s7 cells under `effMaxSG` as a follow-on [default: no].
- Q2 Neighbourhood width: the identifier's own default as found in source, reported [default]; a different width only on Larry's choice at Gate 0.
- Q3 Rows: winner-only evaluated; winner-floor excluded [per Larry].
- Q4 Workers: 100 [default].
- Q5 Compute: go at Gate 1, or an unattended pre-authorization with a wall ceiling (expected ~2 h).

## Done means

Stage 3 report and rendered documents committed; Gate 2 records beside them; the knob landed with the knob-inert identity recorded; alignment quoted in the Stage 0 record; branch left unpushed for Larry.
