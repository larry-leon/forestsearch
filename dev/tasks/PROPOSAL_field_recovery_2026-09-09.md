# PROPOSAL — Field re-selection recovery diagnostics: membership agreement and the rule-name family, in the CV/FB vocabulary

**Date:** 2026-09-09. Author: the Linux MR-field chat. **Status: proposal for Larry's decision (§6); not a task. No code changes until approval.**
**Provenance (verified from the `R_09Sep2026` source Larry uploaded; CC re-verifies against HEAD at Stage 0):** `forestsearch_cross_validation.R` — the membership cross-tab `table(treat.recommend, treat.recommend.original)` giving `sens_H`, `sens_Hc`, `ppv_H`, `ppv_Hc` (l. 1409–1424) and `find_metrics` = `Any`, `Exact`, `At least 1`, `Cov1`, `Cov2`, `Cov 1 & 2`, `Cov1 exact`, `Cov2 exact` (l. 1396–1407); `fs_mr_inference.R` — `kept <- candidates[asm$keep]` (l. 703), `kept[[w]]` (l. 710), `G_out` (l. 857, 864), `sel` (l. 556). Applied precedent: `quarto/applications/gbsg/analysis_gbsg_survival_effMaxSG.qmd` (`cv-original-agreement`, `cv_out$find_metrics[["Exact"]]`, `sens_H`/`ppv_H` reported for GBSG).

## 1. Correction of an earlier claim, for the record

An earlier chat message said the CV/FB recovery diagnostics "cannot run on a real dataset." That is wrong. The CV family's comparator is the **original full-data subgroup** (`treat.recommend.original`), so it runs on any real analysis and does so in the GBSG document. Only `fs_identification_summary.R`'s anchor/partner/proxy classification is simulation-only, because it requires a planted triple. The membership-weighted quantity Larry singled out **already exists** as CV's `sens_H` / `ppv_H`.

## 2. What is proposed

Field analogues of the existing diagnostics, computed from objects the field already holds, in the **same vocabulary**:

- **Membership agreement (primary).** For each outer draw *r* with a re-selected winner `G_out[r]`, compare that candidate's membership `kept[[G_out[r]]]` to the observed pick's `kept[[sel]]`, then average over used draws:
  - `field_sens_H` = mean over draws of |G_r ∩ Ĥ| / |Ĥ| — **the share of the identified patients the field's re-selections retain** (the "membership-weighted containment" Larry asked for);
  - `field_ppv_H` = mean of |G_r ∩ Ĥ| / |G_r| — how concentrated the re-selection is on Ĥ;
  - the complement pair `field_sens_Hc`, `field_ppv_Hc` from the same cross-tab, matching CV's four-metric shape exactly;
  - the distribution, not only the mean: q10 / q50 / q90 of the per-draw containment, since a mean of 0.84 from a tight cloud and from a bimodal one mean different things.
- **Rule-name family (secondary), reusing CV's names.** Over the same draws: `Exact` (= p̂, already recorded), `At least 1`, `Cov1`, `Cov2`, `Cov 1 & 2`, `Cov1 exact`, `Cov2 exact`, computed against the **observed** rule via `.fs_rule_columns()` on each candidate's `sg_def`. Names and definitions transplanted from `forestsearch_cross_validation.R`; no new terms invented.
- Both are deterministic given the multiplier draws already made: **no new randomness, no new fits, no change to any bound.**

## 3. What it is not

- Not an input to any interval. Nothing in the construction reads these; adding them changes no bound, and no p̂- or agreement-dependent interval is proposed (that would be the closed κ line under another name).
- Not a replacement for FB or CV. **FB and CV re-run the search**; the field **re-selects within the fixed kept family** under perturbation. The field's version answers the narrower question — is the pick reproducible given this family — and is free in every analysis; FB/CV answer the stronger question at real compute cost. The record and any printed output must say which is which.
- Not a claim about coverage. These are descriptive; the calibrated meaning of p̂ (high p̂ ⇒ stable-pick regime ⇒ bounds modestly optimistic) comes from the certification record and C-4, not from this proposal.

## 4. Design (add-only, defaults off)

- **Stage 0 (verification, no code):** confirm from HEAD that `kept`, `G_out`, `sel` are in scope at a single insertion point for **both** the `field_complement = TRUE` and `FALSE` paths (`kept` is currently built inside the complement branch at l. 703 — if it is not available when the complement block is off, either hoist its construction (add-only, byte-identical) or gate the diagnostics on the complement path and say so). Quote the lines. STOP and report if neither is clean.
- **Stage 1:** new argument `field_recovery = FALSE` on `fs_mr_inference()`, forwarded to the field block; when TRUE, `field$recovery` gains the membership metrics, their quantiles, and the rule-name family. Template recorder columns `fld_recov_*`. Classification: **adds code; byte-identical defaults**. Identity gate: with `field_recovery = FALSE`, all outputs byte-identical on the standing identity cells; with TRUE, every existing column unchanged, new columns finite.
- **Cost:** set intersections over ≤ 1,000 draws against a membership vector; expected negligible, but measured at Stage 1 against the standing per-replicate reference (fit+MR ≈ 36 s light load) and reported before any campaign.
- **Validation (the point of the exercise):** on one committed simulation cell where FB or CV metrics exist for the same replicates, report the field's `Exact` beside FB/CV's `Exact`, and `field_sens_H` beside CV's `sens_H`. They should be *related but not equal* — different resampling schemes, and the field's is family-conditional. Quantifying that gap is what makes the field's version interpretable to a reader who knows the FB numbers.

## 5. Where it appears in output

- Main results block, beside |Ĥ| and the bounds: `p̂(Ĥ)` (re-selection frequency) and `field_sens_H` (mean share of Ĥ retained under re-selection), one line each with plain-language labels.
- Diagnostics block: the containment quantiles, the top-3 re-selection mass (already recorded), and the rule-name family.
- Vignette: one paragraph distinguishing the field's family-conditional reproducibility from FB/CV re-discovery, with the GBSG numbers as the worked example.

## 6. Decisions (defaults in brackets)

- **R-1** Proceed with the membership agreement as the primary new diagnostic [yes].
- **R-2** Include the rule-name family in the same pass, names transplanted from CV [yes] — or defer it and ship membership only.
- **R-3** Argument name `field_recovery`; column prefix `fld_recov_` [default]. Alternative if Larry prefers CV's own naming: `field_find_*`.
- **R-4** Sequencing: after C-4 reports, folded into the documentation package (NOTE wording, print/summary method, vignette) [default] — or ahead of C-4 if Larry wants it in the GBSG headline analysis sooner.
- **R-5** Validation cell for §4: [one committed cell with CV or FB metrics on the same replicates; CC names it at Stage 0].

## 7. If declined

p̂ alone remains the analysis-time flag, as certified; the FB and CV diagnostics continue to carry the recovery question in analyses that run them. Nothing else changes.
