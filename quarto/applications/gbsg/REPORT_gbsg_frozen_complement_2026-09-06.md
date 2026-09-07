# REPORT — GBSG frozen-family illustration: the complement (benefit-claim) rows and the joint pair (final, reframed)

**Task:** `dev/tasks/TASK_gbsg_frozen_complement_2026-09-06.md` (eac66f3b); Stage 0 record 791448e0; first render and record 2d295d7d. **Reframed per Larry's dispositions (2026-09-06):** the reported complement construction is MR (field); the winner-only / winner-floor IJ variants are rejected on theoretical grounds, removed from Table 3 and the diagnostics, and not evaluated going forward; the reading is by bound location only (no "below 1" clauses). L1–L4 at defaults. No `R/` changes.
**Date:** 2026-09-06. Render: `analysis_gbsg_survival_frozen_family.qmd` with the committed LOO cache (`LOO_CACHE` → `_payloads_2026-09-01_complete/…/cv_out_…_maxeff_neighborhood.rds`), forestsearch 0.3.5 at 8fbb3bdc; the render is the compute. Payload `_payloads/analysis_gbsg_survival_frozen_family_complement/…_payload.rds` (the document's `dirout` set so the committed frozen-family payload is not overwritten); tracked HTML refreshed.

## Identities against the committed payload: ALL PASS (12/12), harm side unchanged

Against `_payloads/analysis_gbsg_survival_frozen_family/…_payload.rds` (ca2a2f93): selected subgroup identical (`{pgr <= 32.5} {er <= 0}`); MR naive est/lower/upper, MR (IJ) `debiased` (every pre-existing element, est 1.753114), the harm field block (21 elements; est2 1.7851, lower_1s 0.8820), the uniform κ block (κ* 1.590), `p_hat`/`n_family`/`selection_bias`, FB `H_estimates` **and** `Hc_estimates`, and G&H (2.2218 / 1.3605 / 0.9617) — all at **0.0 relative difference** (FB via its fixed seed), on every render (three: first, reading-sentence fix, reframe). Complement blocks finite (no degenerate note); bound identities (`up1s = exp(β̃ᶜ − q05)`, `lo2s`, `hi_se`) ≤ 1e-12; `gamma = 0.025 ∈ [0.025, 0.05]`; `ij_residual = "two_term"` recorded. Check script `gbsg_identity.R` (session scratchpad).

## Table 3 — Complement Ĥᶜ (benefit claim), as rendered

Complement of `{pgr <= 32.5} ∩ {er <= 0}` (HR scale; SE log-HR). Primary column: the one-sided 95% **upper** bound. **Reported construction: MR (field)**; naive, FB and MR (IJ, two-term) are comparators.

| Method | Point est. | Two-sided 95% | One-sided 95% upper | SE (log-HR) |
|---|---|---|---|---|
| Naive | 0.61 | (0.47, 0.79) | 0.76 | 0.134 |
| Full bootstrap | 0.62 | (0.41, 0.94) | 0.88 | 0.210 |
| MR (IJ, two-term) | 0.62 | (0.38, 1.02) | 0.94 | 0.249 |
| **MR (field)** | 0.62 | (0.48, 0.82) | **0.79** | 0.139 |

Conventions (the table's footnote): naive and IJ rows `exp(log est + 1.645·SE)`; FB `exp(log H2 + 1.645·sdH2/H2)` (delta method, as Table 1); MR (field) the gate's stored `upper_1s = exp(β̃ᶜ − q₀.₀₅(Λ*ᶜ))`, built on β̃ᶜ, not est2. Complement field: λ̄ᶜ = −0.0002, λ-SDᶜ 0.139, 1000/1000 outer draws, 92 complement fits (1.2% of draw-winner readings needed a new fit), 0.3 s.

## Table 4 — Joint pair (Ĥ lower bound, Ĥᶜ upper bound), as rendered

| Pair | Ĥ lower | Ĥᶜ upper | Joint prob. (field draws) |
|---|---|---|---|
| Separate one-sided 95% field bounds | 0.88 | 0.79 | < 0.95 by construction |
| Bonferroni (γ = 0.025) | 0.76 | 0.82 | 0.950 |
| Calibrated (γ = 0.025) | 0.76 | 0.82 | 0.950 |

γ = 0.025 on the 0.025–0.050 grid; corr(Λ*, Λ*ᶜ) = −0.104 over 1,000 aligned draws — too weak a dependence to lift γ off the Bonferroni floor, so the calibrated pair coincides with Bonferroni's.

## Diagnostics beside the harm side

| Quantity | Ĥ (harm) | Ĥᶜ (complement) |
|---|---|---|
| p̂(Ĥ) | 0.678 | — |
| naive SE (log) | 0.327 | 0.134 |
| IJ two-term SE / naive | 0.484 / 1.48× | 0.249 / **1.85×** |
| field λ-SD / naive; / se_ij | 0.417 / 1.27×; 0.86× | 0.139 / **1.03×**; 0.56× |
| one-sided bound: naive → IJ → field | 1.30 → 0.79 → 0.88 (lower) | 0.76 → 0.94 → 0.79 (upper) |

The field's complement SE coincides with the naive SE (1.03×) because this analysis sits in the dominated-selection regime (p̂(Ĥ) = 0.68: the complement's identity barely changes across re-selection, so the second-order term is negligible, λ̄ᶜ ≈ 0) — a diagnosis of this analysis, not a rule. The two-term IJ's 1.85× is the doubled same-draws term in that regime.

## Timing

Fit + gate 79 s · FB 145 s · MR gate 0.6 s (harm field 0.6 s; **complement field 0.3 s**; κ sweep 75 s) · G&H 28 s · 10-fold CV 13 s · LOO from cache. Whole render ≈ 6 min.

## What the illustration now shows for both subgroups (three sentences)

In the selected subgroup the naive HR of 2.22 carries a one-sided lower bound of 1.30, and the selection-adjusted lower bounds sit at 0.79–0.96 (MR (IJ) 0.79, field 0.88, FB 0.93, G&H 0.96): no clinically meaningful harm is established once selection is accounted for, and the spread across the adjusted bounds is what the constructions pay for selection. In the complement, where treatment is expected to work, the naive HR 0.61 and the selection-adjusted point estimates all sit at 0.62, and the reported field upper bound is 0.79 — below a 0.80 or 0.85 benefit threshold read as an aid — against FB's 0.88 and the two-term IJ's 0.94, the latter the price of its doubled SE (1.85× naive) where the field's λ-SD sits at 1.03× naive because this analysis is in the dominated-selection regime (p̂(Ĥ) = 0.68). Stating both bounds together, the joint 95% pair is (0.76 for Ĥ, 0.82 for Ĥᶜ), the calibrated level coinciding with Bonferroni's because the two field draws are nearly uncorrelated (−0.10).

## Committed with this report

Document (`analysis_gbsg_survival_frozen_family.qmd`: gate-call knobs `include_complement`/`field_complement`, `dirout`, Tables 3–4 and the generated complement reading in `@sec-intervals`, payload schema note and blocks, complement-field timing row), refreshed tracked HTML, the payload directory, this report and the Stage 0 record. Committed payloads untouched. No task proposed; nothing blocked.
