# REPORT — GBSG frozen-family illustration: the complement (benefit-claim) rows and the joint pair (final)

**Task:** `dev/tasks/TASK_gbsg_frozen_complement_2026-09-06.md` (eac66f3b); Stage 0 record 791448e0. L1–L4 at defaults. No `R/` changes.
**Date:** 2026-09-06. Render: `analysis_gbsg_survival_frozen_family.qmd` with the committed LOO cache (`LOO_CACHE` → `_payloads_2026-09-01_complete/…/cv_out_…_maxeff_neighborhood.rds`), forestsearch 0.3.5 at 8fbb3bdc; the render is the compute. New payload `_payloads/analysis_gbsg_survival_frozen_family_complement/analysis_gbsg_survival_frozen_family_complement_payload.rds` (the document's `dirout` set so the committed frozen-family payload is not overwritten); tracked HTML refreshed.

## Stage 1 — Identities against the committed payload: ALL PASS (12/12), harm side unchanged

Against `_payloads/analysis_gbsg_survival_frozen_family/…_payload.rds` (ca2a2f93): selected subgroup identical (`{pgr <= 32.5} {er <= 0}`); MR naive est/lower/upper, MR (IJ) `debiased` (every pre-existing element, est 1.753114), the harm field block (21 elements; est2 1.7851, lower_1s 0.8820), the uniform κ block (κ* 1.590), `p_hat`/`n_family`/`selection_bias`, FB `H_estimates` **and** `Hc_estimates`, and G&H (2.2218 / 1.3605 / 0.9617) — all at **0.0 relative difference** (FB via its fixed seed). Complement blocks finite (no degenerate note); bound identities (`up1s = exp(β̃ᶜ − q05)`, `lo2s`, `hi_se`, `upper_1s_wf`) ≤ 1e-12; `gamma = 0.025 ∈ [0.025, 0.05]`; `ij_residual = "two_term"` recorded. Check script `gbsg_identity.R` (session scratchpad; console record in the session log).

## Table 3 — Complement Ĥᶜ (benefit claim), as rendered

Complement of `{pgr <= 32.5} ∩ {er <= 0}` (HR scale; SE log-HR). Primary column: the one-sided 95% **upper** bound.

| Method | Point est. | Two-sided 95% | One-sided 95% upper | SE (log-HR) |
|---|---|---|---|---|
| Naive | 0.61 | (0.47, 0.79) | 0.76 | 0.134 |
| Full bootstrap | 0.62 | (0.41, 0.94) | 0.88 | 0.210 |
| MR (IJ, two-term) | 0.62 | (0.38, 1.02) | 0.94 | 0.249 |
| MR (IJ, winner-floor) | 0.62 | (0.48, 0.81) | 0.78 | 0.134 |
| MR (field) | 0.62 | (0.48, 0.82) | 0.79 | 0.139 |

Conventions (the table's footnote): naive and IJ rows `exp(log est + 1.645·SE)`; FB `exp(log H2 + 1.645·sdH2/H2)` (delta method, as Table 1); MR (field) the gate's stored `upper_1s = exp(β̃ᶜ − q₀.₀₅(Λ*ᶜ))`, built on β̃ᶜ, not est2; winner-floor = the winner-only IJ residual floored at the naive variance (the unfloored winner-only SE is 0.121; the floor binds at 0.134). Complement field diagnostics: λ̄ᶜ = −0.0002, λ-SDᶜ 0.139, 1000/1000 outer draws, 92 complement fits (1.2% of draw-winner readings needed a new fit), 0.29 s.

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
| IJ two-term SE / naive | 0.484 / 1.48× | 0.249 / 1.85× |
| IJ winner-only SE (unfloored) / naive | 0.219 / 0.67× | 0.121 / 0.90× |
| IJ winner-floor SE / naive | 0.327 / 1.00× (floor binds) | 0.134 / 1.00× (floor binds) |
| field λ-SD / naive; / se_ij | 0.417 / 1.27×; 0.86× | 0.139 / 1.03×; 0.56× |
| one-sided bound: naive → IJ → floor → field | 1.30 → 0.79 → 1.02* → 0.88 (lower) | 0.76 → 0.94 → 0.78 → 0.79 (upper) |

*Harm-side winner-floor lower bound `exp(log 1.753 − 1.645·0.327)` = 1.02 — shown for the comparison only and not read: the harm side stays reported on the two-term and the field, because the refinements evaluation found the floor 1–8 points under nominal in the moderately separated harm regime (p̂ ≈ 0.7 is exactly that regime), where β̃'s dispersion exceeds the naive SE.

## Timing

Fit + gate 78.5 s · FB 146.7 s · MR gate 0.61 s (harm field 0.62 s; **complement field 0.29 s**; κ sweep 75.0 s) · G&H 28.7 s · 10-fold CV 12.9 s · LOO from cache. Whole render ≈ 6 min (rendered twice: the second render only re-worded one generated sentence of the reading; identities re-verified on the second payload).

## What the illustration now shows for both subgroups (three sentences)

In the selected subgroup the naive HR of 2.22 carries a one-sided lower bound of 1.30, and every selection-adjusted construction moves that bound below 1 (FB 0.93, MR (IJ) 0.79, field 0.88, G&H 0.96), so no clinically meaningful harm is established once selection is accounted for — the spread 0.79–0.96 across the adjusted bounds being what the constructions pay for selection. In the complement, where treatment is expected to work, the naive HR 0.61 and the four selection-adjusted point estimates all sit at 0.62, and the one-sided 95% upper bounds are 0.94 (IJ two-term), 0.88 (FB), 0.79 (field) and 0.78 (winner-floor): the benefit claim survives every selection-adjusted construction (all upper bounds below 1), and the field and winner-floor constructions — which agree to 0.01 because the complement is in the dominant-selection regime (p̂(Ĥ) = 0.68) where the naive SE is the right SE for β̃ᶜ — place the benefit at "at most ≈ 0.79", below a 0.85 or 0.80 threshold, while the two-term IJ's doubled SE (1.85× naive) leaves only "at most 0.94". Stating both bounds together, the joint 95% pair is (0.76 for Ĥ, 0.82 for Ĥᶜ) — the calibrated level coincides with Bonferroni's because the two field draws are nearly uncorrelated (−0.10) — so the development claim that can be made jointly is: no meaningful harm established in Ĥ, and a benefit of at most 0.82 in Ĥᶜ.

## Committed with this report

Document (`analysis_gbsg_survival_frozen_family.qmd`: gate-call knobs `include_complement`/`field_complement`, `dirout`, Tables 3–4 and the generated complement reading in `@sec-intervals`, payload schema note and blocks, complement-field timing row), refreshed tracked HTML, the new payload directory, this report and the Stage 0 record. Committed payloads untouched. No task proposed; nothing blocked.
