# CC task — amendment to the self-contained applied OC document: both nulls, and the calibration at the analyst's own gate

**File:** `dev/tasks/cc_task_oc_selfcontained_amend_2026-09-01.md` · **Issued:** 2026-09-01 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first (§1).
**Amends:** `quarto/applications/actg175/analysis_actg175_continuous_oc.qmd` (commit `c73a6b06`) per `REPORT_oc_applied_selfcontained_2026-09-01.md` (`b964db2e`). Two corrections, both chat's errors in the prior spec, both cheap.

---

## ⚠ CATEGORY

**No `R/` change.** Edits: the one document above, its rendered HTML, its payload. Plus this task document. The archived `analysis_actg175_continuous_oc_evaluation.qmd` is **not touched** (its pointer note stands; not re-rendered).

**Compute:** one additional job (the homogeneous-null family + one `fs_oc_grid()`, ≈ 10 min) inside a re-render of the document (≈ 27 min at twelve-wide, ≈ 71 GB peak measured — keep `n_workers = 12` unless memory forces lower; say so if you lower it). **Projected ≈ 30–40 min.** *GATE:* if the re-render projects above **1.5 h**, stop and report.

**Unattended.** Gates stop the task; never ask, never work around. **Provenance gate is dirt-tolerant**; stop only if dirt touches this task's own paths.

---

## 1. Provenance and first commit — GATE

Standard block plus the installed version (expect **0.3.2**). *GATE:* branch `feature/glm-extension`; `cf4e85ba`, `c73a6b06`, `b964db2e` in the log. Copy this document into `dev/tasks/` and commit it alone; report both filenames.

## 2. What is wrong, stated

**(A) The zero-plus rung is mislabelled.** With `k_inter(q) = q − beta_treat`, the within-Q effect is `q` and the outside-Q effect stays `beta_treat = −26.978725`. So `q = 0.01` plants a subgroup whose effect is **zero against a background that benefits by 27** — a DGM with 27 units of genuine effect heterogeneity. That is **the sub-threshold-subgroup null** ("a true subgroup exists but sits far below the thresholds"), not the no-subgroup null. The prior document calls it "the structural null"; that label is wrong and the reading built on it is wrong.

*Confirmation from the archived run, quote it in the report:* at that rung `E[β(Ĥ)] = −22.3` with `E[PPV] = 0.173`, and `−26.978725 + (0.01 − (−26.978725)) × 0.173 = −22.3` — the arithmetic identifies exactly what was computed.

**The no-subgroup null cannot be reached by `k_inter = 0`** under 0.3.2: both region effects are then equal and negative, so `s = sign(m_tau[Q]) = −1` and every candidate's oriented mean becomes `+26.98`, clearing `c1 = 10` trivially (declaration ≈ 1 — the same pathology stage 0 measured on the null branch). §3 builds it correctly instead.

**(B) The calibration is read at the wrong `c2`.** `c2` sets candidate *eligibility*; the calibration asks what truths are consistent with the statistic an analysis **already conducted at the analyst's thresholds** produced. Reading `P(T ≥ T̂_obs | q)` at the diagonal point `(T̂_obs, 0.8·T̂_obs) = (87.92, 70.33)` imposes an eligibility screen (`Bhat ≥ 70.33 + z_p·se_g`) that no analysis ran — which is why the 0.05 crossing moved from the archived ≈ 7.5 to ≈ 61. **Type-I and power stay on the diagonal** (they are about threshold policy); **the calibration moves to the analyst's operating `c2`.**

## 3. Fix A — the homogeneous null, built by construction

Add one visible chunk, immediately before the evaluation loop's null discussion, and one job to the loop. **No new DGM**: under a homogeneous truth every candidate's true mean *is* the common effect, so take the primary family object and set its means directly — the same object-level substitution the archived document's `se_direct` sensitivity used:

```r
fam_hom <- fam                       # the primary family, already built
fam_hom$beta_g <- rep(beta_treat, fam$M)   # homogeneous truth: no subgroup
```

with two sentences of prose saying exactly that, and why `k_inter = 0` is not the route (§2's orientation reason, condensed). Evaluate one `fs_oc_grid()` on `fam_hom` over the same `c1_ladder × c2_vec` crossing, same draws/block/seed, and add it to the results as `variant = "homogeneous"`, `q = NA`.

Report for it, in the type-I section (§7 of the document), **beside** the sub-threshold rung:

- the diagonal false-declaration curve and `c1_05^hom` / `c1_10^hom` at ladder resolution;
- the rate at the analyst's `(10, 0.8×10)` and at the comparator `(10, 5)`;
- **`det_rate` only** — `PQg`, `sens_g`, `spec_g` and every classification functional are undefined with no Q, and the document must say so rather than printing them.

*GATE:* `c1_05^hom ≤ c1_05^sub-threshold` and the homogeneous rate at any `(c1, c2)` is ≤ the sub-threshold rate there (the sub-threshold truth has Q-enriched candidates with elevated means, so it can only declare more often). A violation means the construction is wrong — stop and report.

## 4. Fix B — the calibration at the analyst's gate

Read every `P(T ≥ T̂_obs | q, Q)` at **`c1 = T̂_obs`, `c2 = params$c2_ratio × 10`** (the analyst's operating consistency floor under the policy, = 8) — a point already present in the computed crossing, so **no new draws**. Add a second column at the comparator `c2 = 5`, which is what the archived deep run used, and state in one sentence that the two should sit close and that the archived crossing was ≈ 7.5.

Apply the same change to the knob table's `q at 0.05` / `q at 0.5` columns and to the three-curve calibration figure. Replace the diagonal read-off entirely; do not keep it as a third column.

One sentence of prose, in the calibration section, giving the reason: **`c2` sets eligibility, and the calibration interrogates the analysis as actually conducted, so it is read at the analyst's own consistency floor while type-I and power are read on the policy diagonal.**

## 5. Relabelling and the reading

- Every occurrence of "the zero-plus structural null" / "the structural null" for the `q = 0.01` rung becomes **"the sub-threshold-subgroup null (a true Q with zero effect against the beneficial background)"**, with the homogeneous row named **"the no-subgroup null"**.
- The document must state plainly that **both of the analyst's nulls are now present and that they bracket the question**: no subgroup at all, and a subgroup real but below the thresholds; the second is the more adversarial of the two and is where the larger false-declaration rate comes from.
- The reading section: update every sentence whose number moved (the type-I headline now has two values; the calibration crossings return to their analyst-gate values). **Do not re-assert the prior text's conclusions where the numbers no longer support them** — rewrite from the recomputed values.
- The payload gains the homogeneous row and the corrected calibration; the `c2_policy` element records that the policy applies to type-I and power and **not** to the calibration read-off, with the reason.

## 6. Render and gates

Re-render (repo quarto convention), recording wall-clock and peak memory. *GATES:* §3's monotonicity gate; the anchor, M = 4508 and Q-nesting gates as before; the calibration column at `c2 = 8` and at `c2 = 5` agree within ~2 MC SEs of each other's crossings and the `c2 = 5` crossing is near the archived ≈ 7.5 (report the numbers either way; a large discrepancy is a finding, not a stop); HTML free of `NA ±` and error text; payload written.

## 7. Close-out

Commit by explicit path (qmd, HTML, payload, report). **No push. No install. No `R/` change.** Report `REPORT_oc_selfcontained_amend_2026-09-01.md`: provenance · §2 quoted with the archived arithmetic check · the two nulls' numbers side by side (diagonal type-I at the analyst's gate, `c1_05` for each) · the corrected calibration crossings beside the pre-amendment ≈ 61 and the archived ≈ 7.5 · what changed in the reading · render wall-clock and memory · commits · ten-line summary.

## 8. Out of scope

- No `R/` change — the orientation's behaviour at a homogeneous truth is a recorded limitation, documented in prose here, not a task.
- No change to the archived evaluation document or its HTML; no re-run of its precompute.
- No new Q-variants, no new rungs, no draw-count change, no gbsg, no push.
