# TASK — Field re-selection recovery: membership agreement against the observed subgroup (add-only)

Date: 2026-09-09. Author: chat (spec). Executor: Claude Code (Linux). Approver: Larry (R-1 approved, R-2 deferred, 2026-09-09). Reviewer: the Linux MR-field chat.
Governing proposal: `PROPOSAL_field_recovery_2026-09-09.md` (in `~/Downloads`; commit it to `dev/tasks/` alongside this file). Predecessors: `dev/notes/NOTE_survival_products_2026-09-09.md`, `REPORT_print_vignette_2026-09-09.md`, `R/forestsearch_cross_validation.R` (the CV vocabulary this transplants).

**Decisions fixed by Larry.** **R-1 yes** — membership agreement is the deliverable. **R-2 deferred** — the rule-name family (`Exact`, `At least 1`, `Cov1`, …) is **out of scope**: covariate-name comparison is the fragile part and membership agreement answers the question better. **R-3** `field_recovery`, prefix `fld_recov_`. **R-4** now. **R-5** validation against the FB/CV metrics on a committed cell.

**Why.** The field currently reports one recovery number, p̂(Ĥ) — exact-match frequency. It cannot distinguish "the field re-selected near-twins of Ĥ" from "the field re-selected unrelated regions." The GBSG fit's p̂ = 0.006 over a family of 1,744 makes this concrete: standing alone that number invites misreading. Membership agreement answers it directly and is computable from draws already made.

## Protocol

- First action: archive stale variants in `~/Downloads` to `~/Downloads/cc_archive/` — **do not archive `HANDOFF_guohe_comparison_2026-09-09.md`** (current, belongs to another chat). Copy this file **and** `PROPOSAL_field_recovery_2026-09-09.md` to `dev/tasks/` and commit.
- **Commit only; do not push. Larry reviews before this work is treated as final.** (GitHub Desktop on this machine pushes when driven from its UI; that is Larry's business, not a failure of this task, and it needs no comment in the record.)
- One `R/` change, classified **adds code; byte-identical defaults**, confined to `R/fs_mr_inference.R` plus the one-line `.g_mr` pass-through in `R/forestsearch_main.R` if the knob must reach the gate (verify at Stage 0 — the same pattern as `field_decompose`). Template recorder columns are a document-level change to the survival template. No change to any construction, bound, or default.
- Compute: verification renders (≤ 5 replicates) and one validation read of committed bundles. No campaign.
- Gates stop on failure; **on failure revert the touched files to HEAD, re-install, record the failure, and stop.**
- Standing conventions: transplant-first; verify from source; winner-only and winner-floor excluded; NPV reported alongside sensitivity/specificity/PPV where classification metrics appear. Leave the seven pre-existing untracked files alone.

## Stage 0 — Discovery (quote from HEAD)

1. Quote `kept <- candidates[asm$keep]` with its enclosing `if (isTRUE(include_complement)) {` and enough context to show that `candidates` and `asm$keep` are both in scope **before** that branch. **The hoist below is byte-identical only if `kept` is a pure subset of `candidates` with no fitting, drawing or side effects — confirm this from source and STOP if not.**
2. Quote `G_out` (its construction and its meaning per outer draw), `sel`, `ok_c`, and the field block's gate `if (ci_method == "field")`.
3. Quote the CV membership cross-tab that this transplants — `table(treat.recommend, treat.recommend.original)` and the `sens_H` / `sens_Hc` / `ppv_H` / `ppv_Hc` definitions in `R/forestsearch_cross_validation.R` — so the field's metrics reuse those definitions rather than parallel ones.
4. Quote the template's recorder rows for `fld_Hc_scale_*` (the insertion pattern) and the `field_decompose` knob and its `.g_mr` pass-through.
5. Name a committed cell carrying FB or CV metrics on the same replicates for the R-5 validation; if none exists, say so and report what the closest available comparison is.

## Part R — The construction

**R1. Hoist.** Move `kept <- candidates[asm$keep]` (and `Nall`, `Ncol` if they are equally pure) above the `include_complement` branch, leaving the branch otherwise untouched. Comment that the hoist is byte-identical and why. Gate R2a proves it.

**R2. New argument** `field_recovery = FALSE` on `fs_mr_inference()`, forwarded into the field block; `.g_mr` pass-through in `forestsearch_main.R` if Stage 0 shows the knob cannot otherwise reach the gate. Roxygen: one paragraph stating that these are **descriptive diagnostics computed from draws already made — no new fits, no new randomness, and no construction reads them**, and that they answer a narrower question than FB/CV (re-selection within the fixed kept family under perturbation, not re-discovery from scratch).

**R3. The metrics.** Over the used outer draws `r ∈ ok_c` with re-selected winner `g_r = G_out[r]`, comparing membership `kept[[g_r]]` to the observed `kept[[sel]]` (all set operations on patient indices; `n = Nall`):

- per draw: `a_r = |kept[[g_r]] ∩ kept[[sel]]|`, `b_r = |kept[[g_r]]|`, `c_r = |kept[[sel]]|`;
- `sens_H` = mean over draws of `a_r / c_r` — **the share of the identified patients the re-selections retain** (the primary quantity);
- `ppv_H` = mean of `a_r / b_r`; and the complement pair `sens_Hc`, `ppv_Hc` from the same 2×2, matching the CV definitions quoted at Stage 0. **Also report `npv_Hc`** so the four classification metrics follow the standing convention.
- the distribution of the per-draw containment `a_r / c_r`: q10, q50, q90, and the share equal to 1 (a re-selection containing all of Ĥ).
- guards: draws with `b_r = 0` or an unfit candidate are excluded and counted (`n_used`, `n_skipped`); when `ok_c` is empty every metric is `NA_real_`.

Attach as `field$recovery` (a named list) when `field_recovery = TRUE`; absent when FALSE. Cost: set intersections over ≤ 1,000 draws — **measure it against the standing per-replicate reference and report it.**

**R4. Recorder.** Template gains `FS_S7_FIELD_RECOV` (default FALSE, default-inert) and columns `fld_recov_sens_H`, `fld_recov_ppv_H`, `fld_recov_sens_Hc`, `fld_recov_npv_Hc`, `fld_recov_q10`, `fld_recov_q50`, `fld_recov_q90`, `fld_recov_share1`, `fld_recov_n_used`, filled `%||% NA_real_`, beside the `fld_Hc_scale_*` block.

**R5. `NEWS.md`:** one bullet under the development header.

**Gate R** (5 replicates, ε 0.20 HR 1.50 n500 config, seeds `8316951 + sim_id`, sim_id 1–5, tags `recov_off` / `recov_on`; timing columns excluded):
- **R2a — the hoist is inert.** With `field_recovery = FALSE`, every non-timing column `identical()` to the committed `e1stud` rows 1–5 and `truth` `identical()`; the new columns present and all `NA`.
- **R2b — on.** With `field_recovery = TRUE`, every pre-existing column still `identical()` to the same comparator; the new columns finite; invariants `0 ≤ sens_H, ppv_H, sens_Hc, npv_Hc ≤ 1`, `q10 ≤ q50 ≤ q90`, `n_used + n_skipped = |ok_c|`.
- **R2c — arithmetic check, one replicate, by hand.** Recompute `sens_H` and `ppv_H` for a single replicate directly from `kept`, `G_out` and `sel` in a scratch script, independent of the new code path, and assert agreement to ≤ 1e-12. Quote both numbers.

## Part V2 — Validation and reporting (R-5)

- On the cell named at Stage 0: report the field's exact-match rate (p̂) and `sens_H` **beside** the FB or CV metrics available for the same replicates. They should be **related but not equal** — different resampling schemes, and the field's is family-conditional. Quantify the gap; state plainly that neither is a substitute for the other.
- Add `sens_H` to `summary.forestsearch()`'s post-selection block (one line, beside p̂), guarded on presence so output is unchanged when the diagnostics are absent — re-run the absent-MR invariance check (`Pa` from the print/vignette task) to confirm. **Do not add it to `print()`** in this task; the one-line summary of recovery belongs in `summary()` until the validation is read.
- Update the vignette's p̂ section with one short paragraph and the GBSG fit's `sens_H`, so the p̂ = 0.006 example is no longer presented without it. Rebuild and report the build time.

## Done means

Stage 0 quotes (or a STOP with the reason); Part R committed with R2a/R2b/R2c concrete values and the measured cost; Part V2 committed with the validation table, the `summary()` line, the rebuilt vignette and the re-run invariance check; full test suite re-run and its tally reported (target FAIL 0); `devtools::check()` tally reported against the previous task's; `REPORT_field_recovery_2026-09-09.md` beside the other records; one-paragraph closing summary with the gate results and the commit range. **Out of scope:** the rule-name family (R-2, deferred), `.fs_apply_mr()`'s `ci_method`, any construction or default change, any campaign.
