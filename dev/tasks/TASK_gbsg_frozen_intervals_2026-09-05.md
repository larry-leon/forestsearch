# TASK — Interval constructions illustrated on the GBSG frozen-family analysis

Date: 2026-09-05. Author: chat (spec). Executor: Claude Code. Approver: Larry.
Document: `quarto/applications/gbsg/analysis_gbsg_survival_frozen_family.qmd` (frozen 116-candidate family, all cuts fixed; naive / FB / MR (IJ) / Guo–He / CV / LOO already present; payload under `_payloads/`).
Related records: `quarto/simulations/gbsg_020/REPORT_mr_field_*` (field interval evaluations), `SUMMARY_ij_vs_field_2026-09-05.md`.
Sequencing: after the uniform task's Stage 2/3 finish in the current session (or in a fresh session at reduced workers — decision I6).

## Protocol

- First action: copy this file to `dev/tasks/` and commit. Do not push.
- No change under `R/`. The document stays self-contained: no reads from `dev/` or from simulation directories; wrapper calls visible and few; interpretation inside the document.
- Verify from source: the field block's element names in `fs$mr_inference$field`, the reselection block, and what `fs_bc$H_estimates` carries for FB, before writing any table.
- Committed payloads are read-only; the new render writes a new payload and the identity check compares against the old one.

## Purpose

Add to the frozen-family document a section that presents, for the selected subgroup, every interval construction side by side — two-sided 95% interval, one-sided 95% lower bound on the harm effect, and the SE each is built from — together with the re-selection diagnostics that say which regime the analysis sits in, and a plain-language reading of what the harm claim rests on. The document's existing purpose (the FB–MR frozen-family correspondence) is untouched; the new material is additive and reads the same fit.

## Changes to the document

1. **Gate call** (`{r fit}`): `mr_inference_args = list(draws = mr_draws, ci_method = "field", return_reselection = TRUE, field_uniform = FALSE)`. Under `"field"` the IJ block is unchanged, so every existing MR number reproduces; `field_uniform` is a documented knob left off (decision I1) with a one-line note that the κ-calibrated two-sided interval is under evaluation and will be switched on when adopted.
2. **New section "Interval constructions for the selected subgroup"**, placed after "Corrections in the frozen family":
   - Table 1 (HR scale unless stated), rows Naive / Full bootstrap / MR (IJ) / MR (field) / Guo & He: point estimate; two-sided 95% interval; one-sided 95% lower bound; SE on the log-HR scale; one-sided convention per row in a footnote — naive and IJ `exp(log est − 1.645·SE)`; FB from `fs_bc$H_estimates` (the bootstrap distribution's 5th percentile if the object carries the draws or an SE, decision I2); field `exp(β̃ − q₀.₉₅(Λ*))` from `field$lower_1s`; G&H its own bound. MR (field)'s point estimate is `est2` with the two-sided Λ-quantile interval; the row's footnote says its one-sided bound is built on β̃, not est2.
   - Table 2, re-selection diagnostics from `fs$mr_inference$reselection`: p̂_Ĥ; the top three re-selected candidates with their frequencies and definitions; selection rate; the field's λ̄ (the second-order term, log scale), λ-SD against se_ij and the naive SE; draw usage.
   - Interpretation (prose, inside the document): the target is the effect in the region found with the family held fixed; what each interval's guarantee is — IJ conservative by construction; field one-sided uniformly at nominal in the fixed-family limit experiment and across the simulation cells except the tie regime at larger n under the screen; field two-sided approximate with a known worst case at moderate separation; G&H one-sided by design; κ pending — and then the reading for this subgroup: whether any selection-adjusted lower bound exceeds 1, and what p̂_Ĥ says about how competitive the selection was. State conclusions as the numbers allow; no claim beyond them.
3. **Payload**: add the field block and the reselection block to the reproducibility payload; bump the payload's schema note.
4. **Timing table**: the field block's seconds as a row.

## Stage 0 — Discovery (no compute)

0a. Quote the element names of `fs$mr_inference$field` and `$reselection` under the wrapper at the current package, and the contents of `fs_bc$H_estimates` (whether an SE or the bootstrap draws for H2 are available for a one-sided bound).
0b. Quote the committed payload's fields for naive, FB, MR (IJ), G&H (the identity anchor) and whether FB is seed-reproducible (seed in `forestsearch_bootstrap_dofuture` / plan).
0c. LOO cache presence (`LOO_CACHE` / `gh_dir` path) so the re-render is minutes, not hours.

Gate 0: 0a–0c quoted; STOP if the field block is not present under the wrapper (it is, per the s7 records) or the LOO cache is absent (then decision I4 applies).

## Stage 1 — Edit, render, identity

1a. Apply the changes; render with the LOO cache.
1b. Identity against the committed payload: naive, MR (IJ) est/lower/upper, G&H naive/debiased/bound identical (≤ 1e-12); FB identical if seed-reproducible, otherwise within its Monte Carlo error with the comparison shown; selected subgroup identical. STOP on any MR (IJ) difference.
1c. Field block finite; interval invariants; `lower_1s ≤ est` on the log scale; p̂ sums to the selection rate.

Gate 1: identities pass; rendered HTML present. (No separate compute stage: the render is the compute, ~5 min.)

## Stage 2 — Record

`REPORT_gbsg_frozen_intervals_<date>.md` beside the document: the two tables as rendered, the identity results, the timing row, and a three-sentence statement of what the illustration shows. Commit the document, the new payload, the HTML if the directory's convention tracks it (check `.gitignore`), and the report. No task proposed unless something blocks.

## Decisions (defaults in brackets)

- I1 κ interval now: off, knob present [default]; on, if the uniform task's Stage 3 has been accepted by then.
- I2 FB one-sided bound: bootstrap 5th percentile if available, else `exp(log H2 − 1.645·SE_boot)`, else omitted with a note [default in that order].
- I3 Scope: this document only [default]; the headline maxeff analysis (regenerating family, ~1,340 candidates) as a later, separate task.
- I4 LOO: use the cache; if absent, `run_loo = FALSE` for this render with a note [default].
- I5 Placement: a new section in this document [default]; alternative — a sibling document that re-fits the frozen family for the interval illustration alone.
- I6 Session: queue behind the uniform task in the current session [default]; alternative — a fresh session now with FB workers reduced to avoid contending with the κ runs.

## Done means

Document, payload and report committed; identities recorded; branch left unpushed for Larry.
