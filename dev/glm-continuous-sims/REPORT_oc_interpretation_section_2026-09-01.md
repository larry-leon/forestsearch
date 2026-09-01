# REPORT — the interpretation section folded into the self-contained OC document (2026-09-01)

**Task:** `dev/tasks/cc_task_oc_interpretation_section_2026-09-01.md` (commit
`e7d2c740`)
**Verdict up front: the section is in and every figure in it is inline R —
all gates pass; the optional §5 caching was attempted, failed its own gate
(knitr's cache cannot serialize the homogeneous grid object: lazyLoadDB
"long vectors not supported"), and was reverted entirely per the gate's
instruction; the committed document is uncached.  Render 37.1 min,
computed results bit-identical to `c93890db`.**

## Provenance

- Host `pop-os` · `/home/larryleon/Documents/GitHub/forestsearch` · branch
  `feature/glm-extension` · tree clean at start · installed version
  **0.3.2** ✓ · commits `c79f1b51`, `c93890db`, `3c86d5e8` in the log ✓.
- **Transport names.** Arrived hyphen-stripped as
  `cc_task_oc_interpretation_section_20260901.md` in `~/Downloads`
  (identified by content; expected stem
  `cc_task_oc_interpretation_section_2026-09-01.md`).  Committed in-repo
  under the canonical name
  `dev/tasks/cc_task_oc_interpretation_section_2026-09-01.md`, alone, as
  `e7d2c740` — the first action.

## §2's two computed quantities

- **Unadjusted ITT** (in `data-prep`): −27.591 (SE 7.889) CD4 units —
  printed, kept distinct from the covariate-adjusted `beta_treat`
  (−26.979) throughout; §12's opening names both and says they are not
  conflated.
- **`hom_at_Tobs`** (in the §7.1 chunk): the no-subgroup null's tail
  probability at (c1 = T̂_obs, c2 = 8) = **0.0237 ± 0.0011** — read from
  the homogeneous grid where it already sat unread.

## The section as rendered (§12, after §11, before the payload)

Quoted from the rendered HTML (all figures produced by inline R):

> **12. Interpretation: what the analysis may claim.**  "The analysis
> found Ĥ = {age <= 37} & !{cd40 <= 507} — n(Ĥ) = 66 of 1083, fitted
> mean difference T̂_obs = 87.9 CD4 units, consistency 0.95 — against an
> unadjusted ITT mean difference of −27.6 CD4 units (SE 7.9), a benefit
> on average. […] The caveats governing all of it are section 11's and
> are not repeated here."
>
> **12.1** "At (c1 = 10, c2 = 8) the search declares some harm subgroup
> with probability 0.423 under the no-subgroup null and 0.504 under the
> sub-threshold-subgroup null; at the fixed (10, 5) gate the latter is
> 0.588. […] because the two nulls differ by only 0.081, the rate is
> driven by multiplicity over M = 4508 correlated candidates and a
> liberal threshold, not by a lurking weak subgroup.  The thresholds are
> clinical-relevance floors, not error controls."
>
> **12.2** [table: no subgroup 0.0237 ± 0.0011 · Q = Ĥ 0.0382 ± 0.0014 ·
> wider 0.0483 ± 0.0015 · broad 0.0568 ± 0.0016]  "This is the
> selection-honest tail probability of the observed statistic, computed
> against the distribution of the search's own maximum, so multiplicity
> is already inside it.  The answer is a range — from 0.0237 with no
> subgroup at all to 0.0568 under the broad reading — […] As an external
> cross-check (computed in the archived
> analysis_actg175_continuous_oc_evaluation deep run, not here), the
> primary row's value at 2×10⁵ draws was 0.041.  The observed magnitude
> is suggestive, not decisive, under every null in the family."
>
> **12.3** "The 0.05 consonance crossings of section 10 put the smallest
> planted harm consistent with T̂_obs at q ≈ 8.9 CD4 units if the harmed
> population is exactly Ĥ, but q ≈ 0.5 under the wider reading and
> ≤ 0.01 under the broad one […] **Do not quote T̂_obs as an effect
> size**: the de-biased (multiplier-resampling / bootstrap) estimates are
> the estimation-side response, and this curve is the design-side
> statement of why they are needed."
>
> **12.4** "Holding the family-wise false-declaration rate at 0.05
> requires c1 = 50 on the scaled diagonal (c1 = 40 against the
> no-subgroup null alone, and 87 on the fixed c2 = 5 row) — thresholds
> far above the clinical floor — and at that gate power is 0.58 even at
> the observed-statistic rung […] the finding is hypothesis-generating
> and needs prespecified confirmation — which this same machinery can
> size."
>
> **12.5** "For T̂_obs to be typical (the 0.5 crossing), a subgroup
> exactly Ĥ must be harmed by q ≈ 76.3 CD4 units, the wider population
> by q ≈ 60.1, the broad one by q ≈ 54.0; meanwhile top-rung power rises
> 0.58 → 0.95 across the variants as the declared rule's sensitivity to
> the true Q falls 0.78 → 0.34."
>
> **12.6** (blockquote, labelled *a draft for the analyst to adapt, not a
> conclusion of this document*): "ForestSearch, configured with
> clinical-relevance floors rather than error controls, identified
> Ĥ = {age <= 37} & !{cd40 <= 507} (n = 66 of 1083) with fitted mean
> difference 87.9 CD4 units; since the procedure declares some subgroup
> with probability 0.42–0.50 under nulls anchored to these data, the
> identification is hypothesis-generating.  The observed statistic's
> selection-honest tail probability ranges from 0.0237 (no subgroup) to
> 0.057 (a broad mildly-affected population), and the harm magnitude
> consistent with it at the 0.05 consonance level ranges from ≤ 0.01 to
> 8.9 CD4 units depending on the breadth of the truly affected
> population, so the fitted magnitude should not be quoted as an effect
> estimate.  A prespecified confirmation in independent data — sized with
> the same operating-characteristic machinery — is required before any
> clinical claim."

External facts appear only as cited externals, each named to its source:
the archived deep run's 0.041 at 2×10⁵ draws
(`analysis_actg175_continuous_oc_evaluation`), the ~2-subject E|Ĥ|
machinery offset (verification fixtures, via section 11), and the
archived `se_direct` sensitivity's ≤ 1.4-point movement.

## §3's grep gate — PASS

Scanning the new section with inline-R spans and the one code chunk
stripped, the distinct bare numerals are exactly: the cited external
constants (0.041, 2×10⁵, ~2, 1.4); the thresholds 10 and 5 and the
0.05 / 0.5 consonance levels; section cross-references (7.2, 10, 11,
12–12.6); and letter-adjacent fragments the pattern picks out of "c1",
"c2" and "CD4".  **No typed numeral duplicates a computed quantity** —
every rate, crossing, count and estimate in the section is inline R.

## §5 caching — attempted, failed its gate, reverted

- Caching was enabled in one pass on the six named chunks with explicit
  `dependson` chains and the prescribed `cache.extra` (the `anchor`
  chunk's key necessarily omitted `Q_variants`, which does not exist yet
  at that point in the document).
- **Gate (a) failed:** the cold render died at the `homogeneous-null`
  chunk after 2251 s with
  `Error in lazyLoadDBinsertVariable(...): long vectors not supported
  yet: connections.c:6604` — knitr's lazy-load cache cannot serialize
  the chunk's grid object (`g_hom` carries the full draw machinery,
  beyond the cache database's object-size limit).  Gates (b) and (c)
  were therefore never reached.
- **Per the gate's instruction the caching was reverted entirely** — all
  six option blocks removed, the 114 MB partial cache directory deleted —
  and the uncached document (with the section, which is the deliverable)
  was rendered and committed.  No warm-render timing exists to report.
- For the record: `*_cache/` is not in `.gitignore`, so even a working
  cache would have sat as untracked dirt; nothing was added to
  `.gitignore` (outside this task's edit paths).

## The unchanged-results check — PASS

Against the payload recovered from `c93890db` (`git show`):
`identical()` TRUE on both `extras$q_variants$table` and
`extras$type1$diagonal`.  The §2 additions are new reads of existing
grids, not new truths, as the task stated.

## Payload

`extras$interpretation` added: `tail_table` (truth, rate, se — the four
rows above), `itt` (the full coefficient row: estimate −27.5909, SE
7.8889, t −3.4974, p 0.0005), and `crossings`
(variant, P, q05, q05c, q50, q50c from `knob`).  Nothing else changed.

## Render and gates

- Final (uncached) render: exit 0, **wall-clock 2225 s = 37.1 min**
  (matching the measured 37.1 min baseline; far under 1.5 h), ≈ 85.2 GB
  peak attributable (baseline 33.1 GB → 120.4 GB), `n_workers = 14`.
  The failed caching attempt's render cost an additional 2251 s before
  the revert.
- All existing in-document gates pass (`stopifnot`, exit 0): anchor,
  M = 4508 / s = +1 / m_tau_Q at all 14 jobs, Q-nesting, between-variant
  null gate, null monotonicity pointwise, ladder monotonicity.  HTML: 0
  hits for `NA ±`, 0 for error text.  Payload written with the new
  element.

## Commits

- `e7d2c740` tasks: the task document, alone (first action).
- `a0f732a0` applications: the amended qmd + rendered HTML + payload, by
  explicit path (no cache directory exists to exclude).
- (this file) reports: committed alone, after.

No push.  No install.  No `R/` change.

## Ten-line summary

1. Task arrived hyphen-stripped
   (`cc_task_oc_interpretation_section_20260901.md`), committed
   canonically alone as `e7d2c740`; provenance clean on 0.3.2.
2. Two new computed quantities: unadjusted ITT −27.591 (SE 7.889),
   distinct from `beta_treat` −26.979; and the no-subgroup tail
   probability at (T̂_obs, 8) = 0.0237 ± 0.0011, previously unread.
3. §12 added after §11 with subsections 12.1–12.6 exactly as specified;
   every figure is inline R; §11's caveats referenced, not restated.
4. The tail-probability table brackets the observed statistic:
   0.0237 (no subgroup) → 0.0382 (Ĥ) → 0.0483 (wider) → 0.0568 (broad),
   with the archived 0.041 cited as the external cross-check.
5. The grep gate passes: no typed numeral duplicates a computed
   quantity; bare numbers are only cited externals, thresholds/levels,
   and cross-references.
6. §5 caching failed gate (a) — knitr's lazyLoadDB cannot store the
   homogeneous grid object ("long vectors not supported") — and was
   reverted entirely; the committed document is uncached, the section
   retained.
7. The unchanged-results check passes: `extras$q_variants$table` and
   `extras$type1$diagonal` are `identical()` to the `c93890db` payload.
8. Payload gains `extras$interpretation` (tail table, ITT row, crossing
   summary); nothing else changed.
9. Final render 37.1 min, ≈ 85.2 GB peak attributable, exit 0, HTML
   clean; the failed caching attempt cost one extra 37.5-min render.
10. Commits `e7d2c740` (task), `a0f732a0` (document set), this report
    alone; no push, no install, no `R/` change; archived evaluation
    untouched.

## Deviations

- **Caching reverted** (detailed above) — the task's own fallback path,
  taken without iteration: the failure is mechanical (cache object-size
  limit at `homogeneous-null`), gates (b)/(c) were never evaluated, and
  the task's "one pass" framing forecloses debugging the cache; the
  section is the deliverable and is in.
- **`anchor`'s `cache.extra`** (while caching was in place) omitted
  `Q_variants`, which is only defined two chunks later — the prescribed
  list was applied verbatim to the other five chunks.  Moot after the
  revert.
- **Plot-label shorthand carried over**: §12 uses the established short
  forms ("no-subgroup null", "sub-threshold-subgroup null") from the
  amend-1 relabelling.
- The document's pre-existing double section numbering in the rendered
  HTML ("0.12 12. Interpretation…") affects every section (theme-level
  quarto numbering) and predates this task; not touched.
- None otherwise: no `R/` change, archived evaluation untouched, no new
  rungs/variants/draws, no gbsg, no push.
