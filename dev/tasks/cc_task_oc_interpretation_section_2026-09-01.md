# CC task — fold the interpretive memo into the self-contained OC document as its interpretation section

**File:** `dev/tasks/cc_task_oc_interpretation_section_2026-09-01.md` · **Issued:** 2026-09-01 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first (§1).
**Amends:** `quarto/applications/actg175/analysis_actg175_continuous_oc.qmd` (commit `c93890db`).

**Larry's instruction, 2026-09-01:** the interpretive memo should live **inside** the analysis document — as a memo/appendix section — rather than as a separate file, so the analysis and its interpretation are one source. The memo's text is in `claude/memo_oc_applied_takeaways_v2_2026-09-01.md` in the project KB; **the KB is not on this machine**, so §3 below carries the content to be written. **Every number must be inline R read from this document's own computed objects** — the document's standing rule. The few external facts are cited as such.

---

## ⚠ CATEGORY

**No `R/` change.** Edits: the one document above, its rendered HTML, its payload. Plus this task document. The archived `analysis_actg175_continuous_oc_evaluation.qmd` is **not touched**.

**Compute:** a re-render (measured 37.1 min at `n_workers = 14`, ≈ 86 GB peak), plus §5's optional caching verification (a warm re-render, expected ~1–2 min). *GATE:* if the cold render projects above **1.5 h**, stop and report.

**Unattended.** Gates stop the task; never ask, never work around. **Provenance gate is dirt-tolerant**; stop only if dirt touches this task's own paths.

---

## 1. Provenance and first commit — GATE

Standard block plus the installed version (expect **0.3.2**). *GATE:* branch `feature/glm-extension`; `c79f1b51`, `c93890db`, `3c86d5e8` in the log. Copy this document into `dev/tasks/` and commit it alone; report both filenames.

## 2. Two computed quantities to add first

1. **The unadjusted ITT mean difference.** In §1's `data-prep` chunk (or a two-line chunk beside it): `itt <- coef(summary(lm(y_decline ~ treat, data = actg_df)))["treat", ]`, giving the estimate and SE. This is *not* the document's `beta_treat` (the covariate-adjusted GLM coefficient); the interpretation section refers to the unadjusted ITT, so both must exist and the prose must not conflate them. Print it.
2. **The no-subgroup null's tail probability.** `hom_at_Tobs <- at(och, T_obs, params$c2_ratio * 10)` in the §7.1 chunk — the value already sits in the homogeneous grid and is currently unread.

## 3. The new section

Add **`## 12. Interpretation: what the analysis may claim`** after §11, before the payload chunk. It restates nothing computational — every figure is inline R. Content, in this order:

**Opening.** One short paragraph: the analysis found Ĥ (inline: definition, n(Ĥ) of N, T̂_obs, consistency) against an unadjusted ITT of `itt` — a benefit on average — and the sections above ask what the procedure as configured does under truths anchored to this data; this section says what may therefore be claimed. State that the caveats governing all of it are §11's and are not repeated here.

**12.1 The declaration event is close to uninformative.** At `(c1 = 10, c2 = 8)` the search declares some harm subgroup at the no-subgroup rate (inline) and the sub-threshold rate (inline); at the fixed `(10, 5)` gate the latter is (inline). Two consequences: "ForestSearch identified a subgroup" carries almost no evidential weight at this configuration; and because the two nulls differ by only ~(inline: the difference) the rate is driven by multiplicity over M = (inline) correlated candidates and a liberal threshold, not by a lurking weak subgroup. The thresholds are clinical-relevance floors, not error controls.

**12.2 The observed magnitude — a range, not a p-value.** A small table, built from `hom_at_Tobs$det_rate` and `det0` (the per-variant sub-threshold rates already computed in the `null-not-shared` chunk), with rows: no subgroup (homogeneous); Q = Ĥ with zero effect; Q = wider; Q = broad — each with its MC SE where available. Prose: this is the selection-honest tail probability of the observed statistic, computed against the distribution of the search's own maximum so multiplicity is already inside it; the answer is a **range** that depends on how broad the truly affected population may be, and under the broad reading the observation is not surprising even at zero harm. Cite the archived deep run's 0.041 at 2×10⁵ draws as the external cross-check, naming it as such. Conclusion: suggestive, not decisive under any null in the family.

**12.3 The magnitude of the harm is weakly identified.** The 0.05 consonance crossings per variant (inline from `knob$q05`), against the naive T̂_obs (inline) and the analysis's own `c1 = 10`; note that under the supersets the lower bound collapses toward the first rung, so if the harmed population is not exactly the rule the search returned the data are consistent with arbitrarily small harm. State plainly: **do not quote T̂_obs as an effect size**; the de-biased (MR/bootstrap) estimates are the estimation-side response and this curve is the design-side statement of why they are needed.

**12.4 What error control costs.** The `c2 = 5` versus scaled-`c2` contrast in one paragraph, with `c1_05_hom`, `c1_05_diag` and `c1_05_fix` inline, and the observation that no configuration on this ladder is a test — strict control costs most of the power at every planted effect examined, so the finding is hypothesis-generating and needs prespecified confirmation, which this same machinery can size.

**12.5 Breadth versus severity.** Two or three sentences with `knob$q50` and `knob$q05` inline and the power/sensitivity movement across variants (`knob$ptop`, `knob$sens`): a narrow true subgroup must be severely harmed to make T̂_obs typical, a broad one need only be mildly harmed. Reference §10's table rather than repeating it.

**12.6 A draft results paragraph.** Set in a blockquote and explicitly labelled *a draft for the analyst to adapt, not a conclusion of this document*: three sentences carrying the declaration rates, the tail-probability range, the identifiability statement, and the confirmation recommendation — all inline R, no typed figures.

**External facts, cited as external** wherever they appear (they are not computed here): the archived deep run's values; the ~2-subject `E|Ĥ|` offset from the verification fixtures; the `se_direct` sensitivity's ≤ 1.4-point movement. Name the source document for each.

*GATE:* `grep` the new section for typed numerals that duplicate computed quantities — the only bare numbers permitted are the cited external constants, the thresholds `10` / `5` / `0.8`, and section cross-references. Report what the grep found.

## 4. Payload

Add `extras$interpretation`: the tail-probability table as a data frame (variant, rate, se), the ITT vector, and the crossing summary already in `knob`. Nothing else changes.

## 5. Optional — chunk caching, hard-gated

Every prose edit currently costs a full render. If it can be done cleanly in one pass, enable knitr caching on the heavy chunks (`anchor`, `family`, `fundamental-call`, `evaluation-loop`, `homogeneous-null`, `null-not-shared`) with explicit `dependson`, and `cache.extra = list(params, q_rungs, q_shared, c1_ladder, c2_vec, Q_variants)` so any structural or parameter change invalidates the cache.

*GATE, all three required:* (a) the cold render passes every existing in-document gate; (b) an immediate warm re-render completes in **under 5 minutes**; (c) the warm render's payload is `identical()` to the cold render's on `table`, `extras$type1`, `extras$declared`, `extras$calibration` and `extras$q_variants$table`. **If any fails, revert the caching entirely** (leave the section), say so in the report, and commit the uncached document. The interpretation section is the deliverable; caching is a convenience.

## 6. Render and close-out

Re-render, recording wall-clock and peak memory. *GATES:* every existing gate; §3's grep; HTML free of `NA ±` and error text in the body; payload written with the new element. Confirm the computed results are unchanged from `c93890db` — the two additions in §2 are new reads, not new truths — by `identical()` on `extras$q_variants$table` and `extras$type1$diagonal` against the payload recovered from that commit.

Commit by explicit path (qmd, HTML, payload, report). **No push. No install. No `R/` change.** Report `REPORT_oc_interpretation_section_2026-09-01.md`: provenance · the section as written (quote its rendered prose) · the grep result · §5's three gates and whether caching survived · the unchanged-results check · render wall-clock and memory · commits · ten-line summary.

## 7. Out of scope

- No `R/` change; no new rungs, variants or draw-count change; no change to the archived evaluation; no gbsg; no push.
- No restatement of §11's limitations inside §12 — reference them.
