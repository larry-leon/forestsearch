# REPORT — finalise the self-contained applied OC document: ladder refinement and three text corrections (2026-09-01)

**Task:** `dev/tasks/cc_task_oc_finalise_2026-09-01.md` (commit `e3cdae3c`)
**Verdict up front: all gates PASS — this closes the document.  The
refined 21-evaluation ladder brackets every 0.05 consonance crossing
within five units (primary 8.90 → 9.04, wider 0.54 → 0.89, broad stays
censored at ≤ 0.01 by its own endpoint), the threshold constants are
`identical()` to the previous render, every calibration column is
monotone within 2 MC SEs, and the three §12 corrections are in.  Render
52.0 min at `n_workers = 14`.**

## Provenance

- Host `pop-os` · `/home/larryleon/Documents/GitHub/forestsearch` · branch
  `feature/glm-extension` · tree clean at start · installed version
  **0.3.2** ✓ · the interpretation-section commits — `e7d2c740` (task),
  `a0f732a0` (document set), `eb6a3a7d` (report), hashes read from
  `REPORT_oc_interpretation_section_2026-09-01.md` — in the log ✓.
  (One commit has landed since that report: `258b6228`, Larry's KB memo
  file `claude/memo_oc_applied_takeaways_v2_20260901.md`; it touches
  nothing of this task's paths.)
- **Transport names.** Arrived hyphen-stripped as
  `cc_task_oc_finalise_20260901.md` in `~/Downloads` (identified by
  content; expected stem `cc_task_oc_finalise_2026-09-01.md`).
  Committed in-repo under the canonical name
  `dev/tasks/cc_task_oc_finalise_2026-09-01.md`, alone, as `e3cdae3c` —
  the first action.

## The rung change and the 21-evaluation timing

- Setup-chunk literals exactly as specified:
  `q_rungs <- c(0.01, 5, 10, 15, 20, 40, 60, 87 + 11/12)` (8 primary
  rungs), `q_shared <- c(0.01, 5, 10, 20, 40, 87 + 11/12)` (6 per
  superset); the existing last-slot → `T_obs` replacements in the anchor
  chunk work unchanged.  8 + 6 + 6 = 20 (variant, rung) jobs + the
  homogeneous grid = **21 evaluations**; §6's prose count updated
  (fourteen → twenty) with a clause naming the bracketing purpose.
  Nothing else in the loop changed; no draws change; no caching.
- **Timing:** evaluation loop 2056.4 s = **34.3 min** over 20 jobs at
  `n_workers = 14` (two waves, as projected), homogeneous grid
  sequential after it; total render **3122 s = 52.0 min** — right at the
  task's ≈ 55 min projection.  The §CATEGORY recompute-at-first-wave
  check was performed live at the 30-minute mark: still in the loop,
  host at ≈ 76.6 GB, projection ≈ 55 min ≪ 2 h → continued.
- **Memory:** baseline 33.2 GB, peak 124.2 GB → **≈ 88.9 GB
  attributable**; `n_workers = 14` kept throughout.

## §2's identical-constants check — GATE PASS

Against the payload recovered from `a0f732a0`:

| quantity | identical() |
|---|---|
| `c1_05_diag` (= 50) | TRUE |
| `c1_10_diag` (= 40) | TRUE |
| `c1_05_hom` (= 40) | TRUE |
| `c1_10_hom` (= 40) | TRUE |
| type-I surface (`extras$type1$surface`) | TRUE |
| primary type-I diagonal (`extras$type1$diagonal`) | TRUE |
| knob power columns p20 / p40 / ptop | TRUE |

One adjacent object is worth reporting precisely: the *homogeneous
diagonal* data frame fails a strict `identical()` **on attributes only** —
its values are bit-identical (`all.equal` at tolerance 0 with attributes
ignored is TRUE, as are the homogeneous surface and both point reads);
the data-frame row names differ because the homogeneous rows now sit
after 20 × 156 rather than 14 × 156 rows of `oc`.  No computed value
moved anywhere the inputs did not.

## The before/after crossing table

At c1 = T̂_obs (c2 = 8 and c2 = 5 columns remain identical to each
other, as in every render since the eligibility floor does not bind
there):

| variant | q05 before | q05 after | q50 before | q50 after |
|---|---|---|---|---|
| primary | 8.903 | **9.040** | 76.301 | 76.301 |
| wider | 0.542 | **0.889** | 60.148 | 60.148 |
| broad | ≤ 0.01 | ≤ 0.01 | 54.009 | 54.009 |

(q05c / q50c identical to q05 / q50 throughout.)  The refinement moved
the 0.05 crossings by +0.14 (primary, now interpolated inside (5, 10):
rates 0.0439 → 0.0515) and +0.35 (wider, now inside (0.01, 5): rates
0.0483 → 0.0580) — both bracketed within five units as intended; the
broad crossing remains censored at the first rung because its own
sub-threshold endpoint (0.0568) already exceeds 0.05, and the 0.5
crossings did not move because their bracketing rungs (40, 60, T̂_obs)
were untouched.  §12.3's inline reads follow automatically (9.0 / 0.9 /
≤ 0.01).

## The monotonicity report (report, not a gate)

Computed in the knob chunk and printed in the rendered document:

    calibration column non-decreasing in q within 2 MC SEs:
    primary   wider   broad
       TRUE    TRUE    TRUE

Every variant's calibration column is monotone within tolerance — the
richer ladder introduced no chance dips — so the inline sentence beside
the knob table renders as its clean branch: "Each variant's calibration
column is non-decreasing in q within 2 MC SEs, so the crossing reads are
clean."  (Had any column failed, the sentence names the variants and
flags the first-crossing read; the branch exists in the source.)

## §3's three corrections, as rendered

1. **§12.2 indexed range** — now "from 0.0237 with no subgroup at all to
   0.0568 under the broad reading", read from `tail_tab$rate[1]` and
   `tail_tab$rate[4]` rather than min/max, so the labels can no longer
   drift from the values.
2. **Closing paragraph** — "Three external facts cited above and in
   section 11, named as external: …" (was "Two", listing three).
3. **§12.6 endpoints named** — "…the harm magnitude consistent with it
   at the 0.05 consonance level ranges from essentially zero if the
   affected population is the broad variant to 9.0 CD4 units if it is Ĥ
   itself, so the fitted magnitude should not be quoted as an effect
   estimate." — the censored end described in words, the primary value
   still inline R.

## Render and gates

- All existing in-document gates pass (`stopifnot`, exit 0): anchor
  (66/66, T̂_obs = 87.916667), M = 4508 / s = +1 / |m_tau_Q| = q at all
  20 jobs, Q-nesting, null monotonicity pointwise over all 156 (c1, c2)
  points, ladder monotonicity, the between-variant difference gate
  (max dev 16.35 / 19.75, unchanged).
- HTML: 0 hits for `NA ±`, 0 for error text in the body.
- Payload written; `extras$calibration` now carries 8 primary rungs
  (0.0382, 0.0439, 0.0515, 0.0614, 0.0734, 0.1549, 0.3168, 0.6305) and
  the per-variant curves 6–8 points each.

## Commits

- `e3cdae3c` tasks: the task document, alone (first action).
- `9aae6ce1` applications: the finalised qmd + rendered HTML + payload,
  by explicit path.
- (this file) reports: committed alone, after.

No push.  No install.  No `R/` change.

## Ten-line summary

1. Task arrived hyphen-stripped (`cc_task_oc_finalise_20260901.md`),
   committed canonically alone as `e3cdae3c`; provenance clean on 0.3.2
   with the interpretation commits verified from their report.
2. The ladder grew to 8 + 6 + 6 = 20 jobs plus the homogeneous grid —
   21 evaluations — with the specified literals; nothing else in the
   loop changed, no draws change, no caching.
3. Render 52.0 min at `n_workers = 14` (loop 34.3 min in two waves),
   ≈ 88.9 GB peak attributable; the 30-minute projection recheck kept
   the total far under the 2 h gate.
4. The threshold constants and type-I surface are `identical()` to the
   `a0f732a0` payload; the homogeneous diagonal's values are
   bit-identical with only row-name attributes shifted by the larger
   `oc`.
5. The 0.05 crossings are now bracketed within five units: primary
   8.90 → 9.04 (inside (5, 10)), wider 0.54 → 0.89 (inside (0.01, 5));
   broad stays ≤ 0.01 by its own endpoint; every 0.5 crossing unchanged.
6. The calibration columns are monotone in q within 2 MC SEs for all
   three variants (printed in-document); the knob-table sentence renders
   its clean branch.
7. §12.2 reads its range by explicit row index, the closing paragraph
   says "Three external facts", and §12.6 names its endpoints in words
   with the primary value inline.
8. All existing gates pass; HTML clean; payload written with the refined
   curves.
9. Commits `e3cdae3c` (task), `9aae6ce1` (document set), this report
   alone; no push, no install, no `R/` change; archived evaluation
   untouched.
10. **This closes the document**: the ladder resolves what it reports,
    the interpretation section carries the memo with every figure
    computed in place, and no pending-edits list remains.

## Deviations

- **The homogeneous-diagonal attribute note** (above): reported rather
  than silently passed, since a strict `identical()` on that adjacent
  object is FALSE for row-name reasons while every gated quantity and
  every value is unchanged.
- **The monotonicity sentence is a computed conditional** (inline R with
  a clean branch and a naming branch) rather than static prose, so the
  document stays correct if a future re-render at other settings turns a
  column non-monotone.
- The midpoint projection recheck used elapsed time and live memory
  (the loop prints its timing only at completion); at 30 min the render
  was mid-loop at ≈ 76.6 GB and the ≈ 55 min projection stood.
- None otherwise: no `R/` change, no draws change, no caching revisit,
  no new variants, archived evaluation untouched, no gbsg, no push.
