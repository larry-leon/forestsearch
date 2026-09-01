# CC task — finalise the self-contained applied OC document: ladder refinement and three text corrections

**File:** `dev/tasks/cc_task_oc_finalise_2026-09-01.md` · **Issued:** 2026-09-01 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first (§1).
**Amends:** `quarto/applications/actg175/analysis_actg175_continuous_oc.qmd` (the interpretation-section commit).

**This is the closing pass on this document.** Larry's instruction, 2026-09-01: finish it comprehensively now rather than leave a pending-edits list. Two things: refine the effect ladder so the consonance crossings are properly resolved, and three text corrections from chat's review.

---

## ⚠ CATEGORY

**No `R/` change.** Edits: the one document above, its rendered HTML, its payload. Plus this task document. The archived `analysis_actg175_continuous_oc_evaluation.qmd` is **not touched**.

**Compute:** the ladder grows from 15 grid evaluations to 21. The 15-evaluation render measured **37.1 min at `n_workers = 14`, ≈ 86 GB peak**, so 21 projects to **≈ 55 min at the same width** (two waves). *GATE:* recompute the projection once the first wave's timing is visible; if the total projects above **2 h**, stop and report. Keep `n_workers = 14` — the host absorbed 86 GB comfortably but a wider wave is not needed.

**Unattended.** Gates stop the task; never ask, never work around. **Provenance gate is dirt-tolerant**; stop only if dirt touches this task's own paths.

---

## 1. Provenance and first commit — GATE

Standard block plus the installed version (expect **0.3.2**). *GATE:* branch `feature/glm-extension`, and the interpretation-section commits in the log (read their hashes from `REPORT_oc_interpretation_section_2026-09-01.md`). Copy this document into `dev/tasks/` and commit it alone; report both filenames.

## 2. Ladder refinement — why, and what

The consonance crossings reported in §10 and read in §12.3 are interpolated between adjacent rungs. The primary's 0.05 crossing currently falls in the `q ∈ (0.01, 10)` gap and the supersets' in `(0.01, 20)`, so those numbers are limited by **ladder resolution, not Monte-Carlo noise** — raising `draws` would not improve them. Bracket every 0.05 crossing within five units:

    q_rungs   <- c(0.01, 5, 10, 15, 20, 40, 60, T_obs)          # primary: 8 rungs
    q_shared  <- c(0.01, 5, 10, 20, 40, T_obs)                  # each superset: 6 rungs

(as literals in the setup chunk with `87 + 11/12` in the last slot, replaced by `T_obs` in the anchor chunk exactly as now — the existing `q_rungs[length(q_rungs)] <- T_obs` and `q_shared[length(q_shared)] <- T_obs` assignments still work unchanged.)

That is 8 + 6 + 6 = **20 (variant, rung) jobs plus the homogeneous grid = 21 evaluations.** Nothing else about the loop changes.

**Unaffected by construction, and to be confirmed:** `c1_05_diag`, `c1_10_diag`, `c1_05_hom`, `c1_10_hom` and the type-I surface all read the `q = 0.01` rung and the homogeneous grid only; the knob table's power columns are pinned at `c(20, 40, T_obs)`. *GATE:* these are `identical()` to the previous render's payload values (recover it with `git show <interpretation commit>:<payload path>`), since their inputs have not moved. Report the comparison.

**One new report, not a gate.** With more rungs at rates near 0.05 and MC SE ≈ 0.004, a variant's calibration column could be non-monotone in `q` by chance, which would make `cross()`'s first-crossing read noisy. For each variant, check whether `pT` is non-decreasing in `q` within 2 MC SEs and **print the result**; where it is not, note it in one sentence beside the knob table so the crossing is read with that in mind. Do not gate on it and do not "fix" a non-monotone column.

## 3. Three text corrections in §12

1. **§12.2** — the prose reads "from `min(tail_tab$rate)` with no subgroup at all to `max(tail_tab$rate)` under the broad reading". The homogeneous ≤ primary ordering is gated, but primary ≤ wider ≤ broad is not, so the labels could drift from the values at other settings. Index the rows explicitly instead: `tail_tab$rate[1]` (no subgroup) and `tail_tab$rate[4]` (broad).
2. **§12's closing paragraph** — it says "Two external facts" and then lists three (the archived calibration values, the ~2-subject `E|Ĥ|` offset, the `se_direct` sensitivity). Change to "Three".
3. **§12.6** — "the harm magnitude … ranges from `qfmt(knob$q05[3])` to `qfmt(knob$q05[1])` CD4 units" renders awkwardly when the broad entry is a censored string. Reword to name the endpoints: *"…ranges from essentially zero if the affected population is the broad variant to `r qfmt(knob$q05[1])` CD4 units if it is Ĥ itself"* — keeping the inline R for the primary value and describing the censored end in words.

## 4. Render and gates

Re-render, recording wall-clock and peak memory. *GATES:* every existing in-document gate (anchor, `M = 4508`, per-rung `|m_tau_Q| = q`, orientation, Q-nesting, null monotonicity pointwise, ladder monotonicity, the between-variant difference gate); §2's `identical()` check on the threshold constants; HTML free of `NA ±` and error text in the body; payload written.

Report the **before/after crossings**: the four crossing columns (`q05`, `q05c`, `q50`, `q50c`) per variant beside the previous render's values, with one sentence on how much the refinement moved them.

## 5. Close-out

Commit by explicit path (qmd, HTML, payload, report). **No push. No install. No `R/` change.** Report `REPORT_oc_finalise_2026-09-01.md`: provenance · the rung change and the 21-evaluation timing · §2's identical-constants check · the before/after crossing table · the monotonicity report · §3's three corrections quoted as rendered · render wall-clock and memory · commits · ten-line summary stating plainly that this closes the document.

## 6. Out of scope

- No `R/` change; no `draws` change (the crossings are resolution-bound, and the archived companion holds the high-draw analysis); no caching (reverted in the previous pass and not revisited); no new Q-variants; no change to the archived evaluation; no gbsg; no push.
