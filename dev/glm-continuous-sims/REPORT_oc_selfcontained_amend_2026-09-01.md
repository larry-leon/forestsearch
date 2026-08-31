# REPORT — amendment to the self-contained applied OC document: both nulls, and the calibration at the analyst's gate (2026-09-01)

**Task:** `dev/tasks/cc_task_oc_selfcontained_amend_2026-09-01.md` (commit
`59166fc2`)
**Verdict up front: all gates PASS; the document now carries both of the
analyst's nulls (no-subgroup type-I 0.423 at the analyst's gate against
sub-threshold 0.504) and the calibration crossings return to their
analyst-gate values (0.05 at q ≈ 8.9, against the pre-amendment ≈ 61 and
the archived ≈ 7.5); re-render 36.1 min at twelve-wide, exit 0.**

## Provenance

- Host `pop-os` · `/home/larryleon/Documents/GitHub/forestsearch` · branch
  `feature/glm-extension` · tree clean at start (no dirt at all, so the
  dirt-tolerant gate had nothing to tolerate) · installed version
  **0.3.2** ✓ · commits `cf4e85ba`, `c73a6b06`, `b964db2e` in the log ✓.
- **Transport names.** The task arrived hyphen-stripped as
  `cc_task_oc_selfcontained_amend_20260901.md` in `~/Downloads`
  (identified by content; the expected stem was
  `cc_task_oc_selfcontained_amend_2026-09-01.md`).  Committed in-repo
  under the canonical name
  `dev/tasks/cc_task_oc_selfcontained_amend_2026-09-01.md`, alone, as
  `59166fc2` — the first action.

## §2 of the task, quoted against the archived arithmetic

The task's diagnosis (A): with `k_inter(q) = q − beta_treat`, the
q = 0.01 rung plants a subgroup whose within-Q effect is ≈ 0 against a
background that benefits by `beta_treat = −26.978725` — the
**sub-threshold-subgroup null**, not the no-subgroup null.  The archived
run's own numbers identify exactly what was computed, as the task
instructed be quoted: at that rung **E[β(Ĥ)] = −22.3 with
E[PPV] = 0.173**, and

> −26.978725 + (0.01 − (−26.978725)) × 0.173 = −22.31 ≈ −22.3 ✓

— the declared rule's expected true effect is the PPV-weighted mixture of
the planted zero and the beneficial background, i.e. a truth with ≈ 27
CD4 units of genuine effect heterogeneity.  (This document's own §7.3
sub-threshold row reads E[β(Ĥ)] = −22.29 at E[PPV] = 0.173 — same
arithmetic at 2e4 draws.)

And (B): the prior document read `P(T ≥ T̂_obs | q)` at the diagonal
point (87.92, 70.33), imposing an eligibility screen no analysis ran,
which is why its 0.05 crossing sat at ≈ 61 instead of the archived ≈ 7.5.

## What was done (all in `analysis_actg175_continuous_oc.qmd`; no `R/` change)

- **Fix A** — one visible chunk (`homogeneous-null`) immediately before
  the evaluation loop's null discussion: `fam_hom <- fam;
  fam_hom$beta_g <- rep(beta_treat, fam$M)` (the same object-level
  substitution the archived `se_direct` sensitivity used), one
  `fs_oc_grid()` over the same `c1_ladder × c2_vec` crossing at the same
  draws/block/seed, appended to the results as `variant = "homogeneous"`,
  `q = NA`.  Prose states the construction and why `k_inter = 0` is not
  the route (orientation flips every mean to +26.98; declaration ≈ 1).
- **Fix B** — every calibration read-off moved from the diagonal
  (T̂_obs, 0.8·T̂_obs) to **c1 = T̂_obs, c2 = 8** (the analyst's
  operating floor, `params$c2_ratio × 10`), with a second column at the
  comparator **c2 = 5**; applied to the §7.2 ladder table, the §10 knob
  table's crossing columns (now four), and the calibration figure (solid
  c2 = 8, dashed c2 = 5).  The diagonal read-off is gone, not kept as a
  third column.  One-sentence reason in §7.2: c2 sets eligibility; the
  calibration interrogates the analysis as conducted.
- **Relabelling** — every "zero-plus structural null" / "structural null"
  became "the sub-threshold-subgroup null (a true Q with zero effect
  against the beneficial background)" (shortened to "sub-threshold null"
  inside plot labels); the homogeneous row is "the no-subgroup null"; a
  new §6 paragraph states plainly that both nulls are present and bracket
  the question, the sub-threshold being the more adversarial; §11 records
  the `k_inter = 0` orientation behaviour as a limitation in prose.
- **Payload** — gains `extras$type1$homogeneous` (diagonal, c1_05/c1_10,
  rates at (10, 8) and (10, 5), full surface), the corrected
  `calibration` (pT/pT5 columns), and a `c2_policy` element recording
  that the policy applies to type-I and power and not to the calibration
  read-off, with the reason.
- The archived `analysis_actg175_continuous_oc_evaluation.qmd` untouched,
  not re-rendered.

## The two nulls, side by side (diagonal c2 = 0.8·c1, 2e4 draws)

| quantity | sub-threshold-subgroup null (q = 0.01) | no-subgroup null (homogeneous) |
|---|---|---|
| rate at the analyst's (10, 8) | **0.5042 ± 0.0035** | **0.4229 ± 0.0035** |
| rate at the comparator (10, 5) | 0.5884 ± 0.0035 | 0.5078 ± 0.0035 |
| c1_05 (ladder resolution) | 50 (rate 0.0273) | **40** (rate 0.0442) |
| c1_10 (ladder resolution) | 40 | 40 |
| classification functionals | full set | `det_rate` only — undefined with no Q, and the document says so |

*GATE (§3, enforced as `stopifnot` in the render):* c1_05^hom = 40 ≤
c1_05^sub-threshold = 50 ✓, and the homogeneous rate ≤ the sub-threshold
rate **pointwise over all 156 (c1, c2) points** ✓ (same seed and se_g, so
common random numbers make the elevated-means ordering exact).  The
sub-threshold null is confirmed as the more adversarial of the two
everywhere; the no-subgroup null is far from innocuous — 0.42 at the
analyst's own gate.

## The corrected calibration, beside the pre-amendment and archived values

Primary Q, P(T ≥ T̂_obs | q) crossings at c1 = T̂_obs:

| read-off | 0.05 crossing | 0.5 crossing |
|---|---|---|
| **amended, c2 = 8 (analyst's floor)** | **q ≈ 8.90** | **q ≈ 76.3** |
| **amended, c2 = 5 (comparator)** | **q ≈ 8.90** | **q ≈ 76.3** |
| pre-amendment (diagonal, c2 ≈ 70.3) | q ≈ 61.1 | beyond the ladder |
| archived deep run (c2 = 5, 2e5 draws, 11 rungs) | q ≈ 7.54 | q ≈ 75.95 |

- **Finding, not a stop (per §6):** the c2 = 8 and c2 = 5 columns do not
  merely agree within 2 MC SEs — they are **identical to four decimals at
  every rung** (e.g. 0.6305 ± 0.0034 at the top rung under both).  At
  c1 = T̂_obs = 87.92 the eligibility floor of 8 vs 5 excludes no
  candidate that the declaration screen would have admitted, on the same
  draw set.  The two crossings are therefore trivially within tolerance.
- The c2 = 5 crossing q ≈ 8.90 is near the archived ≈ 7.54; the gap is
  ladder resolution (linear interpolation between the q = 0.01 rung at
  0.0382 and the q = 10 rung at 0.0515 — a shallow segment where a
  1-MC-SE wobble moves the crossing by several q units) plus 2e4 vs 2e5
  draws, not configuration.  The 0.5 crossing 76.3 vs archived 75.95
  agrees closely.
- Knob crossings (0.05, at the analyst's c2): primary 8.90, wider 3.19,
  broad 2.28 — against the pre-amendment 61.1 / 43.8 / 42.2.

## What changed in the reading

- The type-I headline has **two values** (0.504 sub-threshold, 0.423
  no-subgroup at (10, 8)), and §6 states that the nulls bracket the
  question — no subgroup at all, and a subgroup real but below the
  thresholds — with the sub-threshold the more adversarial.
- The prior text's claim that the q = 0.01 rung is "the same truth as
  'no subgroup' to a hundredth of a CD4 unit" is gone — it was the wrong
  reading built on the wrong label; that rung carries 27 CD4 units of
  genuine heterogeneity.
- The calibration is no longer "far more conservative under the c2
  policy": that conservatism was an artifact of reading at an eligibility
  screen no analysis ran.  At the analyst's own gate the crossings sit
  where the archived run put them, and the breadth–severity reading of
  §10 survives on the corrected numbers (0.5 crossing 76.3 → 60.2 → 54.0
  across primary → wider → broad).
- §11 now states the policy's scope explicitly: the c2 policy governs
  type-I and power (threshold-policy quantities, read on the diagonal),
  not the calibration read-off.

## Render and gates

- Re-render (RStudio-bundled quarto 1.9.38): exit 0, **wall-clock
  2163 s = 36.1 min** — inside the projected 30–40 min, far under the
  1.5 h stop; `n_workers = 12` kept; evaluation loop 1107.5 s = 18.5 min
  over the twelve jobs, the homogeneous grid sequential after it.
- **Memory:** host used sampled every 15 s — baseline 32.8 GB, peak
  97.1 GB → ≈ 62.8 GB attributable; no pressure on 251 GB.
- Anchor gate (66/66, T̂_obs = 87.916667), M = 4508 with s = +1 at all
  twelve jobs, Q-nesting gate, §3 monotonicity gate — all pass (render
  exit 0 with every gate a `stopifnot`).  HTML: 0 hits for `NA ±`, 0 for
  error text.  Payload written with the homogeneous row and corrected
  calibration.

## Commits

- `59166fc2` tasks: the task document, alone (first action).
- `5ba546b3` applications: the amended qmd + its rendered HTML + its
  payload, by explicit path.
- (this file) reports: committed alone, after.

No push.  No install.  No `R/` change.

## Ten-line summary

1. Task arrived hyphen-stripped (`cc_task_oc_selfcontained_amend_20260901.md`),
   committed canonically alone as `59166fc2`; provenance clean on 0.3.2.
2. The q = 0.01 rung is relabelled throughout as the
   sub-threshold-subgroup null; the archived arithmetic
   (−26.978725 + 26.988725 × 0.173 = −22.3) is quoted as the identification.
3. The no-subgroup null is built by construction — `fam_hom$beta_g`
   overwritten with `beta_treat` — in one visible chunk, one extra
   `fs_oc_grid()` job, same crossing/draws/block/seed.
4. `k_inter = 0` is documented (prose, §6 and §11) as unusable: the 0.3.2
   orientation flips every mean to +26.98, declaration ≈ 1.
5. Two-null headline at the analyst's (10, 8): 0.5042 sub-threshold vs
   0.4229 no-subgroup; c1_05 = 50 vs 40; monotonicity gate (pointwise
   over all 156 grid points) passes.
6. Every calibration read-off moved to (c1 = T̂_obs, c2 = 8) with a
   c2 = 5 comparator column; the diagonal read-off is gone.
7. The corrected 0.05 crossing is q ≈ 8.90 (both c2 columns identical to
   four decimals — the eligibility floor does not bind at c1 = T̂_obs),
   against the pre-amendment ≈ 61 and near the archived ≈ 7.54; the 0.5
   crossing 76.3 vs archived 75.95.
8. The reading was rewritten from the recomputed values: two-value
   type-I, the bracketing statement, no re-assertion of the diagonal
   read-off's conservatism.
9. Re-render 36.1 min at twelve-wide, ≈ 62.8 GB peak attributable, exit
   0; HTML clean; payload gains the homogeneous row, the pT/pT5
   calibration, and the c2_policy scope note.
10. Commits `59166fc2` (task), `5ba546b3` (document set), this report
    alone; no push, no install, no `R/` change; the archived evaluation
    document untouched.

## Deviations

- **Label shortening in plots:** the full relabel "the
  sub-threshold-subgroup null (a true Q with zero effect against the
  beneficial background)" is used in prose and first mentions; axis
  labels, figure legends, and the surface title use "sub-threshold null"
  / "sub-threshold-subgroup null" without the parenthetical for space.
- **The homogeneous job runs sequentially after the loop**, in the
  visible chunk itself, rather than as a thirteenth `mclapply` job — the
  chunk needs `fam` and the loop's results object, and one sequential
  ~9-minute grid is what the task's §CATEGORY costed (render landed at
  36.1 min, inside 30–40).
- **The c2 = 8 vs c2 = 5 calibration columns are identical**, not merely
  close (reported above as a finding); the knob table therefore shows
  duplicate values in its paired crossing columns, by construction, and
  the document's required "should sit close" sentence stands beside them.
- Peak attributable memory ≈ 62.8 GB, below the prior run's ≈ 70.6 GB
  (the 15 s sampler may straddle peaks differently); `n_workers = 12`
  was never lowered.
- None otherwise: no `R/` change, no touch of the archived evaluation
  document, no new DGM, no draw-count change, no push.
