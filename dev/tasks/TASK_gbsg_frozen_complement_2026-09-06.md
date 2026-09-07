# TASK — GBSG frozen-family illustration: the complement (benefit-claim) rows and the joint pair

Date: 2026-09-06. Author: chat (spec). Executor: Claude Code. Approver: Larry.
Sequencing: after `TASK_mr_field_complement_2026-09-06.md` and `TASK_complement_refinements_2026-09-06.md` have both reached Stage 3 (this document uses `field_complement`, `ij_residual`, and `field$joint`).
Document: `quarto/applications/gbsg/analysis_gbsg_survival_frozen_family.qmd` (section `@sec-intervals` added by `TASK_gbsg_frozen_intervals_2026-09-05.md`; payload under `_payloads/`).

## Protocol

- First action: copy this file to `dev/tasks/` and commit. Do not push.
- No change under `R/`. The document stays self-contained (no reads from `dev/` or simulation directories; wrapper calls visible and few; interpretation inside).
- Verify from source: the element names of `fs$mr_inference$complement`, `fs$mr_inference$field$complement`, `field$joint`, and the winner-only SE elements; what `fs_bc` carries for the complement (FB `Hc` estimates, if any); the committed payload's harm-side fields (the identity anchor).
- Committed payloads read-only; the render writes a new payload.

## Changes

1. **Gate call**: `mr_inference_args` gains `include_complement = TRUE, field_complement = TRUE`; `ij_residual` left at its default (`"two_term"`) so the reported harm-side numbers are unchanged; the winner-only SEs are read from the additional elements. Everything else as committed.
2. **Section `@sec-intervals` gains**:
   - **Table 3 — Complement Ĥᶜ (benefit claim)**, rows Naive / Full bootstrap (if `fs_bc` carries the complement; otherwise omitted with a note) / MR (IJ, two-term) / MR (IJ, winner-floor) / MR (field): point estimate (HR); two-sided 95% interval; **one-sided 95% upper bound** (primary; naive and IJ rows `exp(log est + 1.645·SE)`, field from `field$complement$upper_1s`); SE on the log-HR scale. Footnote: orientation and why the upper bound is the exposed limit for a benefit claim.
   - **Table 4 — Joint pair**: (Ĥ lower bound, Ĥᶜ upper bound) for (i) the two separate 95% field bounds, (ii) Bonferroni, (iii) the calibrated pair from `field$joint`, with γ and the (Λ*, Λ*ᶜ) correlation; footnote on what a joint 95% guarantee means for a claim that states both.
   - **Reading** (inside the document): whether the benefit claim in the complement survives every selection-adjusted construction (upper bounds below 1), set beside the harm-claim reading already there; p̂_Ĥ and the complement's λ-SD against se_ij and the naive SE; the width of the winner-only IJ against the two-term IJ on the complement. Conclusions as the numbers allow.
3. **Payload**: the complement, complement-field, winner-only and joint blocks added.
4. **Timing table**: complement-field seconds.

## Stages

- **Stage 0** (no compute): quote the element names above and the payload's harm-side fields; confirm the LOO cache.
- **Stage 1**: edit, render with the LOO cache (~6–10 min). Identities against the committed payload: naive, FB, MR (IJ), MR (field, harm side), κ, G&H — all unchanged (≤ 1e-12; FB via its fixed seed); the selected subgroup identical. STOP on any harm-side difference. Complement blocks finite; bound identities; `gamma ∈ [0.025, 0.05]`.
- **Stage 2**: `REPORT_gbsg_frozen_complement_<date>.md` beside the document: Tables 3–4 as rendered, identities, timing, and a three-sentence statement of what the illustration now shows for both subgroups. Commit document, payload, HTML per the directory's convention, and the report.

## Decisions (defaults in brackets)

- L1 Complement rows: as listed [default]; FB complement only if already computed by the committed FB block.
- L2 Joint pair: included [default].
- L3 κ for the complement: not built; not shown [default].
- L4 Placement: extend `@sec-intervals` in this document [default].

## Done means

Document, payload and report committed; harm-side identities recorded; branch left unpushed for Larry.
