# REPORT — amendment 2 to the self-contained applied OC document: each Q-variant gets its own sub-threshold null (2026-09-01)

**Task:** `dev/tasks/cc_task_oc_selfcontained_amend2_2026-09-01.md` (commit
`c79f1b51`)
**Verdict up front: all gates PASS.  The sub-threshold null is measurably
not shared — between-variant |Δbeta_g| reaches 16.3 (wider) and 19.8
(broad) CD4 units — and the corrected wider/broad 0.05 crossings drop to
0.54 and ≤ 0.01 (from 3.19 and 2.28), while the primary column is
bit-identical to the previous render.  Re-render 37.1 min at
`n_workers = 14`, exit 0.**

## Provenance

- Host `pop-os` · `/home/larryleon/Documents/GitHub/forestsearch` · branch
  `feature/glm-extension` · tree clean at start (nothing for the
  dirt-tolerant gate to tolerate) · installed version **0.3.2** ✓ ·
  commits `59166fc2`, `5ba546b3`, `8b4e31c7` in the log ✓.
- **Transport names.** Arrived hyphen-stripped as
  `cc_task_oc_selfcontained_amend2_20260901.md` in `~/Downloads`
  (identified by content; expected stem
  `cc_task_oc_selfcontained_amend2_2026-09-01.md`).  Committed in-repo
  under the canonical name
  `dev/tasks/cc_task_oc_selfcontained_amend2_2026-09-01.md`, alone, as
  `c79f1b51` — the first action.

## §2 restated, with the measured between-variant differences

The previous document asserted the q = 0.01 sub-threshold null is shared
across Q-variants and used the primary's grid as every variant's left
calibration endpoint.  It is not shared: at q → 0 a candidate's mean is
`beta_treat`·(1 − PQg^v), Q-dependent through PQg^v.  The old
`null-shared` check compared each family to **its own** Q-dependent
background — an algebraic identity (deviation = q·PQg ≤ 0.01 by
construction) that passes for every variant and cannot see the
between-variant difference.  The measured facts, from the new
`null-not-shared` chunk (member-aligned families, alignment asserted on
identical `lab` vectors):

| vs primary, at q = 0.01 | max \|Δbeta_g\| | median \|Δbeta_g\| |
|---|---|---|
| wider | **16.35 CD4 units** | 1.56 |
| broad | **19.75 CD4 units** | 3.34 |

*GATE (§3):* between-variant max > 1 CD4 unit for at least one superset —
passes by an order of magnitude; the premise of the amendment is
confirmed (`stopifnot` inside the render; exit 0).

Each variant's own declaration rate at (T̂_obs, 8) on its own
sub-threshold null: primary 0.0382, wider 0.0483, broad 0.0568 — against
the single 0.0382 the shared shortcut had been lending to all three.

## What was done (all in `analysis_actg175_continuous_oc.qmd`; no `R/` change)

- `q_shared <- c(0.01, 20, 40, T_obs)` — each superset runs its own
  sub-threshold rung; `jobs` = 6 + 4 + 4 = **14**; YAML `n_workers`
  default 12 → **14**.
- `knob-table`: the shared-null branch
  (`d <- if (qv == q_rungs[1]) ocp else …`) is deleted — every variant
  reads its own rows for every q, including 0.01.  The power-column rungs
  are pinned to `c(20, 40, T_obs)` explicitly (previously indexed off
  `q_shared`, which would have silently shifted with the new element).
- The `null-shared` chunk is **replaced** by `null-not-shared`: per-superset
  max/median |beta_g^v − beta_g^primary| over the M members, each
  variant's own rate at (T̂_obs, 8), the > 1-CD4-unit gate, and two
  sentences of prose stating that the null is a property of the planted Q
  and that the shortcut was invalid and removed.
- `knob-figure` caption and §10 prose lose the "sharing the sub-threshold
  null point" claim (now "each anchored by its own sub-threshold null
  point").
- §7.1 gains one sentence: each superset's own sub-threshold null carries
  its own type-I — at (10, 8): wider 0.5675 ± 0.0035, broad
  0.6172 ± 0.0034 — with full diagonal curves in the payload
  (`extras$type1$subthreshold_variants`); the headline stays the
  primary's (0.5042).

## §4's three hygiene items

1. **Payload trap** — `g_hom$table[, c("Eppv", "Esens", "Espec", "Enpv")]
   <- NA_real_` before the homogeneous rows enter `oc`, with the
   one-line comment (fam_hom inherits PQg/sens_g/spec_g from the primary
   family; those functionals were being computed against a Q that does
   not exist under that truth).  `EnH`, `EbetaH`, `mass_below` stay as
   computed.  Verified in the payload: all four columns NA on every
   homogeneous row; the HTML body has zero `NA ±` hits (the homogeneous
   rows are only ever printed as `det_rate`).
2. **Per-rung truth gate** — `eval_job` now asserts
   `abs(abs(fam_q$orientation$m_tau_Q) - jobs$q[i]) < 1e-9` alongside the
   `M` and `s` gates (read from the family's stored `orientation$m_tau_Q`,
   which `fs_oc_family_enumerate` records from the DGM's region scale).
   Held at all 14 jobs (render exit 0).
3. **Inert argument** — `adverse_outcome = TRUE` in `dgm_at()` was
   **kept, with a comment** ("inert for a continuous outcome: the
   negation branch is binary-only") rather than dropped — functionally
   identical and it documents the inertness at the call site.

## The corrected knob table, beside the previous render's

Crossings of P(T ≥ T̂_obs | q, Q) at c1 = T̂_obs (c2 = 8 and c2 = 5
columns remain identical to each other, as in the previous render — the
eligibility floor does not bind at c1 = T̂_obs):

| variant | q at 0.05 (corrected) | q at 0.05 (previous) | q at 0.5 (corrected) | q at 0.5 (previous) |
|---|---|---|---|---|
| primary | **8.90** | 8.90 — unchanged | **76.30** | 76.30 — unchanged |
| wider | **0.54** | 3.19 | **60.15** | 60.15 |
| broad | **≤ 0.01** | 2.28 | **54.01** | 54.01 |

- The corrections run in the direction §2 predicted: the shared shortcut
  understated the supersets' low-q rates, pushing their 0.05 crossings
  conservatively high.  With its own null the broad variant's left
  endpoint (0.0568) already exceeds 0.05, so its crossing clamps to the
  first rung and the document prints "≤ 0.01" by its existing `qfmt`
  convention.
- The 0.5 crossings are unchanged — they interpolate between q = 40 and
  T̂_obs, where no endpoint enters — so the breadth–severity reading's
  headline numbers (76.3 → 60.2 → 54.0) stand.

## The primary-unchanged check (bit-exact, vs the `5ba546b3` payload)

- Primary knob row (q05, q50, p20, p40, ptop, sens, ppv): identical at
  tolerance 0 ✓.
- Primary calibration column (all six rungs): `identical()` TRUE ✓.
- Primary type-I diagonal and homogeneous diagonal: `identical()` TRUE ✓
  (0.5042 at (10, 8), c1_05 = 50 / c1_10 = 40; homogeneous 0.4229,
  c1_05^hom = 40).
- Wider/broad power columns p20/p40/ptop: identical ✓ (only their
  calibration endpoints — and hence 0.05 crossings — moved, as scoped).

## Render and gates

- Re-render (RStudio-bundled quarto 1.9.38): exit 0, **wall-clock
  2225 s = 37.1 min** (within the ≈ 40 min projection, far under 1.5 h);
  evaluation loop 1169.6 s = 19.5 min over 14 jobs at `n_workers = 14`.
- **Memory:** baseline 33.3 GB, peak 121.4 GB → **≈ 86.1 GB
  attributable** — above the task's ≈ 75 GB projection (14 concurrent
  jobs plus the parent), with wide margin on the 251 GB host;
  `n_workers = 14` never lowered.
- Gates, all enforced as `stopifnot` in the render: §3 between-variant
  gate ✓ · anchor (66/66, T̂_obs = 87.916667) ✓ · M = 4508 and s = +1 at
  all 14 jobs ✓ · per-rung m_tau_Q truth gate at all 14 jobs ✓ ·
  Q-nesting ✓ · null monotonicity pointwise over all 156 (c1, c2) points
  ✓ · ladder monotonicity ✓ · family alignment (`lab` identical) in
  `null-not-shared` ✓.  HTML: 0 hits for `NA ±`, 0 for error text.
  Payload written with `subthreshold_variants` and `null_not_shared`
  (replacing the invalid `null_dev`).

## Commits

- `c79f1b51` tasks: the task document, alone (first action).
- `c93890db` applications: the amended qmd + rendered HTML + payload, by
  explicit path.
- (this file) reports: committed alone, after.

No push.  No install.  No `R/` change.

## Ten-line summary

1. Task arrived hyphen-stripped
   (`cc_task_oc_selfcontained_amend2_20260901.md`), committed canonically
   alone as `c79f1b51`; provenance clean on 0.3.2.
2. The shared-null shortcut was removed: each superset now runs its own
   q = 0.01 rung — 14 (variant, rung) jobs at `n_workers = 14`.
3. The invalid identity check was replaced by a between-variant
   difference report with its own gate: max |Δbeta_g| vs primary 16.35
   (wider) / 19.75 (broad) CD4 units — gate > 1 passes decisively.
4. Per-variant own-null rates at (T̂_obs, 8): 0.0382 / 0.0483 / 0.0568;
   corrected 0.05 crossings 8.90 / 0.54 / ≤ 0.01 against the previous
   8.90 / 3.19 / 2.28 — the shortcut's bias was conservative, as scoped.
5. The 0.5 crossings (76.30 / 60.15 / 54.01) and every power column are
   unchanged; the primary column is bit-identical to the `5ba546b3`
   payload.
6. §7.1 now reports each superset's own sub-threshold type-I (0.5675
   wider, 0.6172 broad at (10, 8)); the headline remains the primary's
   0.5042.
7. Hygiene: homogeneous rows' Eppv/Esens/Espec/Enpv NA'd before entering
   `oc` (verified in payload); per-rung
   |`m_tau_Q`| = q gate added in `eval_job`; `adverse_outcome = TRUE`
   kept with an inert-here comment.
8. The knob's power rungs were pinned to `c(20, 40, T_obs)` explicitly so
   the enlarged `q_shared` could not silently shift the p20/p40 columns.
9. Re-render 37.1 min, ≈ 86.1 GB peak attributable (above the ≈ 75 GB
   projection, ample host margin), exit 0; HTML clean; payload written.
10. Commits `c79f1b51` (task), `c93890db` (document set), this report
    alone; no push, no install, no `R/` change; archived evaluation
    untouched.

## Deviations

- **Member-alignment assertion added** to `null-not-shared`
  (`identical(lab)` across the three q = 0.01 families): the task's
  max/median of member-wise differences presumes aligned families; a
  silent misalignment would have made the report garbage, so it is
  asserted rather than assumed.
- **Power-rung pinning** (`sapply(c(20, 40, T_obs), …)` replacing
  `sapply(q_shared, …)` in the knob chunk): not named in the task, but
  required by §3's `q_shared` change — otherwise p20/p40/ptop would have
  silently become p0.01/p20/p40.  Verified bit-identical to the previous
  render's power columns.
- **Peak memory ≈ 86.1 GB** vs the task's ≈ 75 GB projection (the
  15 s-sampled peak caught the full 14-wide wave); no pressure at
  251 GB, `n_workers = 14` kept.
- **§4.3 choice**: kept `adverse_outcome = TRUE` with the inert-here
  comment (the task offered drop-or-keep; keeping changes nothing and
  documents the inertness).
- None otherwise: no `R/` change, archived evaluation untouched, no new
  variants/rungs beyond the two added jobs, no draw-count change, no
  push.
