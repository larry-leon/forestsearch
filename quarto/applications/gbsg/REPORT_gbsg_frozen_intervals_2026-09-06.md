# REPORT — GBSG frozen-family interval constructions (final)

**Task:** `dev/tasks/TASK_gbsg_frozen_intervals_2026-09-05.md` (084ec961). Gate 0 STOP (878d13c2 — `return_reselection` not forwarded) resolved by Larry: add-only pass-throughs committed at **c1b957b2** (`return_reselection` + `field_M_cap` on `forestsearch()`, defaults FALSE/NULL, no behaviour change; Rd regenerated; suite 0/4941/3/32; fixed-seed spot byte-identity ij + field both exact). **I1 revised:** `field_uniform = TRUE` with the mass-carrying cap lifted to the full family (`field_M_cap = 1000L ≥ K = 133`); I2–I6 at defaults.
**Date:** 2026-09-06. Render: LOO from the committed cache (`LOO_CACHE` env → `_payloads_2026-09-01_complete/...`); the render is the compute.

## Stage 1 — Identities: ALL PASS (12/12)

Against the committed payload `_payloads_2026-09-01_complete/analysis_gbsg_survival_frozen_family/analysis_gbsg_survival_frozen_family_payload.rds`:

- Selected subgroup identical: `{pgr <= 32.5} {er <= 0}`.
- Naive est/lower/upper ≤ 1e-12. **MR (IJ)** est/lower/upper/lower_1s/se/se_ij/se_wald/var_ij ≤ 1e-12 (STOP gate — est 1.753114 reproduced exactly under `ci_method = "field"`); n_family (133) and selection_bias identical. **G&H** naive/debiased/bound ≤ 1e-12 (STOP gate — 2.22 / 1.3605 / 0.9617). **FB** H_estimates identical at 1e-12 (seed-reproducible via the formal default `seed = 8316951L`, as recorded at Gate 0).
- 1c: field block finite, no degenerate note, `lower_2s ≤ upper_2s`, `lower_1s ≤ est` (log scale); p̂ sums to the selection rate (1.000 — every multiplier draw produced a winner under the plain maxeff argmax); uniform block finite with the κ-widened interval containing the plain quantile interval.

## Table 1 — Interval constructions for the selected subgroup (HR scale; SE log-HR)

| Method | Point est. | Two-sided 95% | One-sided 95% lower | SE (log) |
|---|---|---|---|---|
| Naive | 2.22 | (1.17, 4.22) | 1.30 | 0.327 |
| Full bootstrap | 1.96 | (0.81, 4.74) | 0.93 | 0.451 |
| MR (IJ) | 1.75 | (0.68, 4.52) | 0.79 | 0.484 |
| MR (field) | 1.79 | (0.76, 3.78) | 0.88 | 0.417 |
| MR (field, uniform 2s) | 1.79 | (0.46, 5.89) | 0.88 | 0.417 |
| Guo & He (2021) | 1.36 | — (one-sided by design) | 0.96 | — |

One-sided conventions: naive and MR (IJ) `exp(log est − 1.645·SE)` (IJ's is the gate's stored selection-adjusted bound); FB `exp(log H2 − 1.645·sdH2/H2)` (I2 middle option — the object carries an SE, not draws); the field's bound is built on β̃, not est2; the uniform row is the conservative uniformly-valid option (κ* widening of the field quantiles about their mean), computed on this trial's own Σ̂ with the cap lifted to the full family.

## Table 2 — Re-selection and field diagnostics

| Quantity | Value |
|---|---|
| p̂(Ĥ)  [selected: q6.1 & q8.1] | **0.678** |
| Top re-selected: q6.1 & q8.1 | 0.678 |
| 2nd: q8.1 | 0.106 |
| 3rd: q1.0 & q2.1 | 0.092 |
| Selection rate | 1.000 |
| Field λ̄ (second-order term, log) | −0.0181 |
| λ-SD vs se_ij | 0.417 / 0.484 = 0.86 |
| λ-SD vs naive SE | 0.417 / 0.327 = 1.27 |
| Draw usage | 1000/1000 outer; 500/500 inner |
| Uniform κ* (mcse) | 1.59 (0.054) |
| Uniform M / family; mass; minC₁ | 17 / 133; 0.991; 0.930 |

(The full-family cap was inert in the best way: the ≥ 0.99 mass target was reached at M = 17 — the cap never bound.)

## Timing

Fit + gate 85 s · FB 148 s · MR gate 0.6 s (field block 0.43 s; **κ sweep 82.1 s**) · G&H 28 s · LOO from cache. Whole render ≈ 6 min.

## What the illustration shows (three sentences)

Every selection-adjusted construction pulls the naive HR of 2.22 down (FB 1.96, MR (IJ) 1.75, field est2 1.79, G&H 1.36), and all four selection-adjusted one-sided 95% lower bounds fall below 1 (0.93, 0.79, 0.88, 0.96) while the naive bound sits at 1.30 — so the one-sided harm claim survives only the analysis that ignores selection. The re-selection diagnostics place this analysis in the dominant-selection regime (p̂(Ĥ) = 0.68 against a runner-up at 0.11), where the second-order field term is small (λ̄ = −0.018) and the constructions differ mostly in width, not location. The κ(Σ̂)-widened interval (κ* = 1.59, mass 0.991 at M = 17 of 133) prices the uniformly-valid two-sided statement at (0.46, 5.89) against the field's (0.76, 3.78) — width bought as a guarantee, changing no conclusion here.

## Committed with this report

Document (`analysis_gbsg_survival_frozen_family.qmd` — gate-call knobs, the new `@sec-intervals` section with both tables and the generated reading, payload schema bump, timing rows), refreshed tracked HTML, and the new payload under `_payloads/analysis_gbsg_survival_frozen_family/` (committed payloads under `_payloads_2026-09-01_complete/` untouched). No task proposed; nothing blocked after the Gate 0 resolution.
