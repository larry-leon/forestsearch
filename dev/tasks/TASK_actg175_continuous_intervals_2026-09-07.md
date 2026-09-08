# TASK — ACTG175 continuous: the interval constructions in the headline analysis document

**File:** `dev/tasks/TASK_actg175_continuous_intervals_2026-09-07.md` · **Issued:** 2026-09-07 by chat, commissioned by Larry
**Machine:** Mac-Studio-3 · **Repo:** `~/Documents/GitHub/forestsearch` · **Branch:** `feature/glm-extension-mac` (continue on it; HEAD ≥ 8cbba125)
**Transport:** `~/Downloads` → `dev/tasks/`, committed first (§1).
**Precedent to mirror:** the intervals section of `quarto/applications/gbsg/analysis_gbsg_survival_frozen_family.qmd` (added 2026-09-05/06: one flag on the gate call, a table of two-sided intervals and one-sided bounds, re-selection diagnostics, the complement rows and the Bonferroni pair, interpretation by bound location inside the document). Read it as the pattern; do not edit it.

**Purpose.** Campaign `mdf1` (REPORT_continuous_field_2026-09-07.md) established the field constructions on the continuous/MD path. This task makes them visible where an analyst would look: the headline ACTG175 continuous analysis document gains the current interval set for the found subgroup Ĥ and its complement — naive, MR (IJ two-term), MR (field) with its one-sided bounds as the directional products, the Bonferroni pair for a claim on both — with p̂(Ĥ), read by bound location against mean-difference thresholds inside the document. It is the continuous mirror of the GBSG frozen-family intervals section and gives the first real-data reading of the complement's selection-adjusted upper bound on this path.

---

## Standing conventions (govern this session)

1. **Unattended.** Gates are stop-on-failure: on a failed gate, commit what is green, write the record with the failure stated, stop. Never ask, never work around.
2. **Verify from source.** Records quote path, signature, chunk label and line numbers from the current tree.
3. **No `R/` change of any kind.** No edit to any package file. If the document needs something the package does not provide on this path, that is a finding and a stop.
4. **Off-limits (Linux nb20 campaign in flight on `feature/glm-extension`):** `quarto/simulations/gbsg_020/*`, `R/fs_mr_inference.R`, `R/forestsearch_main.R`. Never touch `feature/glm-extension`; never push.
5. **Self-contained analysis document.** No reads from `dev/` or from any simulation directory; wrapper calls visible and few; the interpretation is a section of the document, not a separate memo. Follow the directory's existing conventions for the `_payloads/` object and for tracking rendered HTML.
6. **Rows:** naive, MR (IJ two-term), MR (field), Bonferroni pair; FB only if a committed FB result for this document already exists in its payload (never run a new bootstrap); no Guo–He row (the adapters are survival-only); **no winner-only, winner-floor, κ / `field_uniform`, hybrid κ, covariate adjustment.**
7. **Interpretation convention:** every bound is read by its location against the MD thresholds (N-3), never as significance at the null; no "below/above the null", no "significant" anywhere in the new text. The spread across the adjusted bounds is reported as the price of selection. The complement's upper bound (benefit claim) is the practically important product and is read first.
8. **Do not re-run committed work.** The document's existing results (the found Ĥ, its naive estimate, any committed MR (IJ) interval) must reproduce, not be replaced; identity against the committed payload is the gate.
9. The paper is not a topic.
10. Report raw output, not summaries; every number in the record is computed, never typed.

---

## ⚠ CATEGORY

**No `R/` change.** **One existing application document edited** (flags on its gate call, one new section, one new payload element) plus its re-rendered HTML and payload per the directory's convention. **Compute:** one or two renders of that document (a GBSG-class render was ≈ 5 min; the ACTG175 continuous search with GRF may take longer — measured in §2). **Pre-authorized: ≤ 2 h cumulative wall, 1 h hard timeout per render** (use the Perl alarm wrapper committed under `scripts_mdf1/`). If §2's projection exceeds the ceiling, stop and report.

---

## Decisions (all at defaults; unattended)

- **N-1 — the document.** The headline ACTG175 continuous analysis document is the one whose `forestsearch()` call on the continuous ACTG175 outcome yields the anchor subgroup **Ĥ = {age ≤ 37 & cd40 > 507}, n = 66** (the anchor recorded by the applied-OC Stage 0 of 2026-08-31; T̂ ≈ 87.92; N = 1,083; ITT ≈ −27.6 on the raw scale). §2 identifies it from source; if the anchor does not reproduce in any candidate, stop.
- **N-2 — rows.** As in convention 6.
- **N-3 — thresholds.** Oriented scale: positive = harm, the gate's working scale (−cd4_change), as in `summary_continuous_field_mdf1.qmd`. Harm block, one-sided lower bound read against 0, 10, 20, 30, 40 CD4 cells/mm³ (0 = no harm; 10 = the design's consistency threshold; 30 = the search's effect threshold; 40 = the planted MD of the simulation design). Complement block, one-sided upper bound read against 0, 10, 20, 30. Every quantity is also shown on the raw `cd4_change` scale in the same table (sign flipped), since that is how a reader of the ACTG175 analysis thinks. Larry may replace these thresholds later; state in the document that they are the simulation design's reference points.
- **N-4 — compute.** Pre-authorized as in CATEGORY.

---

## 1. Provenance and first commit — GATE

```bash
cd ~/Documents/GitHub/forestsearch
hostname; pwd; git branch --show-current; git rev-parse --short HEAD; git status --porcelain; git log --oneline -3
Rscript -e 'cat(as.character(packageVersion("forestsearch")), "\n")'
```

*GATE:* branch `feature/glm-extension-mac`, HEAD at or after `8cbba125`, installed 0.3.5. The 16 known untracked entries (the mdf1 identity/calibration bundles and the pre-existing payload dirt) are left alone; **stop only if** dirt touches this task's own paths (the N-1 document, its `_payloads/` subdirectory, `dev/tasks/`). Copy this document into `dev/tasks/` and commit it alone.

**First repo action after that commit (carried over from the mdf1 close-out):** two wording edits to `quarto/simulations/actg175/continuous/REPORT_continuous_field_2026-09-07.md`, no re-runs, one commit: (i) in the p̂ / tie-regime reading, state that the enumerated ties are between labels with identical membership on the analysis sample, so p̂ computed on labels understates the settledness of the membership, which is why this tie regime is benign for the one-sided bound; note that a membership-based p̂ would be the more informative diagnostic (finding only, no task); add that in the md40 and null cells the field's nominal one-sided coverage combines λ-SD/SD ≈ 1.2 with +0.4 SD retained bias, whereas at md120 both are near 1 and 0. (ii) Wherever the null cell's bound locations are quoted, label its truth explicitly as "no subgroup; homogeneous +26 on the harm-oriented scale" so the 24% is read as bound location against +26, not as a false-claim rate.

## 2. Stage 0 — locate and anchor — GATE 0

Record: `REPORT_actg175_continuous_intervals_stage0_2026-09-07.md` beside the document.

2.1 List every `.qmd` under `quarto/applications/actg175/` that calls `forestsearch()` with a continuous outcome; for each quote the call (outcome, effect measure, `sg_focus`, thresholds, the covariate list), whether it calls the MR gate today (and with which arguments), and its committed payload path and contents. Identify the N-1 document by the anchor: read its committed payload (not the HTML) and confirm Ĥ = {age ≤ 37 & cd40 > 507}, n = 66. *GATE:* exactly one document reproduces the anchor from its committed payload. Several or none: stop.

2.2 From the committed payload, list the numbers the document already reports that must be unchanged after this task (the found rule, its membership count, the naive estimate and interval, the MR (IJ) interval if present, any FB result). These are the identity anchors for §3.

2.3 Quote the gate call as it stands and state, from `?forestsearch` and the 0.3.5 signature, that `ci_method = "field"`, `field_complement = TRUE`, `return_reselection = TRUE` reach the gate on the continuous path (the mdf1 Stage 0 record established this for the template; confirm the same entry point is used here). Quote where the returned object carries the field block, the complement block, `field$joint` and p̂.

2.4 Render cost: from the committed HTML's timing (if the document records it) or from a one-off timing of the search chunk alone, project the full render. *GATE:* projection ≤ 1 h per render; else stop.

## 3. Stage 1 — the flags and the section — GATE 1

3.1 On the document's gate call add `ci_method = "field"`, `field_complement = TRUE`, `return_reselection = TRUE`; nothing else on the call. Store the returned interval objects in the payload under a new element (name per the directory's convention), leaving every existing payload element byte-identical.

3.2 New section, placed after the document's existing subgroup/estimation results and titled for what it is (e.g. "Post-selection intervals for the found subgroup and its complement"), containing — every number inline R from the payload, nothing typed:

- **Table 1 — Ĥ:** rows naive, MR (IJ), MR (field) [FB if committed]: point estimate (IJ: β̃; field: est₂ and β̃), SE (IJ SE; field λ-SD), two-sided 95% interval, **one-sided 95% lower bound**; on the oriented scale and on the raw `cd4_change` scale. Plus a line: the retained-optimism spread naive → IJ → field on the oriented scale.
- **Table 2 — Ĥᶜ:** the same rows with the **one-sided 95% upper bound**; plus λ-SDᶜ / naive SE as the regime diagnostic.
- **Table 3 — joint:** the Bonferroni pair (Ĥ lower at γ = 0.025, Ĥᶜ upper at γ = 0.025) and the calibrated pair from `field$joint`, with the realised γ.
- **Diagnostics:** p̂(Ĥ) with the two-line reading (near 1 = settled; below ~0.5 = tie regime), the top re-selected candidates and their re-selection frequencies if the object carries them.
- **Reading (in the document):** for each adjusted bound, its location against the N-3 thresholds in plain language — the complement's field upper bound first ("harm at most τ is supported / not supported on the complement"), then Ĥ's field lower bound, then the IJ bounds as the conservative reference, then the Bonferroni pair for a claim on both. The spread across adjusted bounds stated as the price of selection. No significance language.
- **Provenance line:** package version, `ci_method`, draws, seed, render time.

3.3 Render with the timeout wrapper. *GATE 1:* (a) every §2.2 anchor reproduces from the new payload exactly (rule, n, naive estimate; MR (IJ) interval identical if it existed — the IJ block is untouched by `ci_method = "field"`, as mdf1 showed bit-identical); (b) the new section renders with all numbers finite; (c) the field's bound identities hold (one-sided lower = β̃ − q₀.₉₅(Λ*), upper = β̃ − q₀.₀₅(Λ*ᶜ)); (d) no other section's rendered numbers changed (diff the rendered text of the old and new HTML outside the new section; list any difference).

3.4 Commit by explicit path: the `.qmd`, the payload, the rendered HTML per the directory's convention.

## 4. Stage 2 — the record

`REPORT_actg175_continuous_intervals_2026-09-07.md` beside the document: provenance; §2's classification and anchors; the three tables as rendered (pasted from the payload, not retyped); the diagnostics; the Gate 1 diff result; render wall; commits (`git log --oneline 8cbba125..HEAD`); a ten-line reading ending with the two numbers Larry will read first — the complement's field one-sided upper bound and Ĥ's field one-sided lower bound, both on the raw `cd4_change` scale with their oriented values in parentheses. Findings go in the record; a task is proposed only if something blocks.

## 5. Out of scope

No `R/` change; no new bootstrap; no Guo–He row; no edit to any other application document (the binary and survival ACTG175/GBSG documents, the OC-evaluation document); no simulation; no threshold other than N-3's; no touch of the off-limits paths; no push.
