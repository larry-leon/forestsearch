# REPORT — applied OC evaluation: type-I, the effect ladder, and the calibration curve (2026-08-31)

**Task:** `dev/tasks/cc_task_oc_applied_evaluation_2026-08-31.md` (commit `fa3a7d34`)
**Verdict up front: all gates PASS; the evaluation ran in full, the triggered
`se_g` sensitivity included; the document renders in 8 s.**

## Provenance

- Host `pop-os` · `/home/larryleon/Documents/GitHub/forestsearch` · branch
  `feature/glm-extension` · start HEAD `15ff5640` · installed version
  **0.3.2** ✓ · stage-0 report commit `209a8e85` and signed-orientation
  commits `eb136a35`, `15ff5640` in the log ✓.
- §1 dirt check: `git status --porcelain` **empty** — no untracked or
  modified paths at all, so the dirt-tolerant clause had nothing to list.
- Task document committed alone as `fa3a7d34`; the `~/Downloads` stem
  arrived hyphen-stripped as `cc_task_oc_applied_evaluation_20260831.md`,
  committed as `cc_task_oc_applied_evaluation_2026-08-31.md`.
- Scripts (all under `dev/glm-continuous-sims/`, committed with their logs
  and objects): `oc_applied_eval_build_2026-08-31.R` (§2–§3),
  `oc_applied_eval_rung_2026-08-31.R` (§4, one process per rung),
  `oc_applied_eval_seg_band_2026-08-31.R` (§5),
  `oc_applied_eval_sens_seg_2026-08-31.R` (§5 sensitivity),
  `oc_applied_eval_summary_2026-08-31.R` (§4 post-processing). R 4.6.1.

## §2 The anchored context — GATE: PASS

From `stage0_oc_applied_2026-08-31.rds`, not memory: N = 1083;
`beta_treat = −26.978724502131`; `k_inter(q) = q + 26.978725`;
`T̂_obs = 87.916666667` (= 87 + 11/12, the stored full-precision value, used
throughout; the task's literal 87.916667 is its 6-dp print); clause mapping
`age = 37`, `cd40 = list(type = "greater", value = 507)`. The reconstructed
flag on the replicated analysis frame matches the anchor membership
**66/66, exact** (`setequal` on ids and count).

## §3 The rungs — GATES: PASS at all eleven

q ∈ {0.01, 2.5, 5, 7.5, 10, 15, 20, 30, 40, 60, 87.916667}. Per rung
(build log, committed):

- `|m_tau[Q]| − q` ≤ 1.8e-15 (bound 1e-9) ✓; `m_tau[Qc] − beta_treat`
  ≤ 3.6e-15 ✓ (`fs_dgm_scale` readback).
- `orientation$s = +1` at every rung ✓; **M = 4508** at every rung ✓.
- `lab`, `Pg`, `PQg`, `se_g`, `sens_g`, `spec_g`, `ovl`, `memb` all
  **`identical()`** (bit-identical) to the first rung's ✓.
- `max |beta_g(q) − beta_g(0.01) − (q − 0.01)·PQg|` ≤ 1.4e-14 (bound 1e-9) ✓.

## §4 The evaluation — one 2e5-draw set per rung

Per rung: `fs_oc_grid(family, n = 1083, c1 = sort(c(0:200, T̂_obs)),
c2 = 5, "resample", pconsistency = 0.90, draws = 2e5, block = 5e4,
seed = 8316951)`. Eleven concurrent detached processes (all 11 at once,
≤ 12 allowed); the zero-plus rung derived `c1_05`/`c1_10` from its own grid
and published them; the others picked them up on completion. Checkpoints:
grid `table` + full `fs_oc_predict` objects at the named points only
(c1 = 10, c1 = T̂_obs, c1_05 = 86, c1_10 = 79); no per-threshold `results`
lists saved.

**Compute gate (6 h wall-clock), recomputed at first completion as
required:** first rung measured **15 237 s ≈ 4.23 h** — the task's ~67-min
estimate missed the blocked sweep's per-threshold cost (draw ≈ 5 760 s,
sweep ≈ 9 880 s over 202 points). With all rungs concurrent and identically
sized, the projection at first completion was max-of-stragglers ≈ 4.4 h
plus the sensitivity (launched 15 min later) ≈ **4.7 h < 6 h → continue.**
Realized: all eleven done at 4.40 h; sensitivity at ≈ 4.6 h after launch.

### §4.1 Type-I (zero-plus structural null, q = 0.01)

- **Headline: type-I at the analyst's (c1 = 10, c2 = 5) = 0.58500
  (MC SE 0.00110).**
- At 1-unit resolution (the document's stated resolution):
  **c1_05 = 86** (rate 0.04961 ± 0.00049), **c1_10 = 79**
  (rate 0.09818 ± 0.00067).

### §4.2 The ladder

| q | rate at c1 = 10 ± SE | power at c1_05 = 86 ± SE | power at c1_10 = 79 ± SE | largest int c1 with ≥ 0.80 / 0.90 / 0.95 |
|---|---|---|---|---|
| 0.01 | 0.5850 ± 0.0011 | 0.0496 ± 0.0005 | 0.0982 ± 0.0007 | — / — / — |
| 2.5 | 0.5945 ± 0.0011 | 0.0528 ± 0.0005 | 0.1030 ± 0.0007 | — / — / — |
| 5 | 0.6048 ± 0.0011 | 0.0566 ± 0.0005 | 0.1083 ± 0.0007 | — / — / — |
| 7.5 | 0.6156 ± 0.0011 | 0.0604 ± 0.0005 | 0.1147 ± 0.0007 | — / — / — |
| 10 | 0.6267 ± 0.0011 | 0.0650 ± 0.0006 | 0.1213 ± 0.0007 | — / — / — |
| 15 | 0.6497 ± 0.0011 | 0.0756 ± 0.0006 | 0.1377 ± 0.0008 | — / — / — |
| 20 | 0.6733 ± 0.0010 | 0.0891 ± 0.0006 | 0.1572 ± 0.0008 | — / — / — |
| 30 | 0.7248 ± 0.0010 | 0.1268 ± 0.0007 | 0.2087 ± 0.0009 | — / — / — |
| 40 | 0.7789 ± 0.0009 | 0.1819 ± 0.0009 | 0.2779 ± 0.0010 | — / — / — |
| 60 | 0.8784 ± 0.0007 | 0.3500 ± 0.0011 | 0.4639 ± 0.0011 | 59 / — / — |
| 87.9167 | 0.9660 ± 0.0004 | 0.6604 ± 0.0011 | 0.7534 ± 0.0010 | 75 / 65 / 57 |

(The task's "smallest integer c1 with rate ≥ 0.80/0.90/0.95" is trivially 0
on a nonincreasing curve; read as the **largest** such c1 — the breadth
ladder's c1\* convention. "—" = not attained at any integer c1 ≥ 0.)

**No rung on the ladder reaches 80% power at c1_05** — the
observed-statistic rung attains 0.660. Power columns monotone in q ✓.

### §4.3 What is declared at (c1 = 10, c2 = 5)

| q | E\|Ĥ\| ± SE | E[PPV] ± SE | E[sens] | E[spec] | E[β(Ĥ)] ± SE | naive bias ± ~SE | mass below c1 ± SE |
|---|---|---|---|---|---|---|---|
| 0.01 | 75.77 ± 0.06 | 0.174 ± 0.0009 | 0.180 | 0.938 | −22.29 ± 0.021 | 89.93 ± ~0.09 | 1.0000 ± 0 |
| 2.5 | 75.60 ± 0.05 | 0.186 ± 0.0009 | 0.192 | 0.939 | −21.49 ± 0.024 | 89.38 ± ~0.09 | 1.0000 ± 0 |
| 5 | 75.42 ± 0.05 | 0.200 ± 0.0010 | 0.206 | 0.940 | −20.59 ± 0.026 | 88.74 ± ~0.09 | 1.0000 ± 0 |
| 7.5 | 75.20 ± 0.05 | 0.214 ± 0.0011 | 0.219 | 0.941 | −19.59 ± 0.029 | 88.04 ± ~0.09 | 1.0000 ± 0 |
| 10 | 74.99 ± 0.05 | 0.229 ± 0.0011 | 0.234 | 0.942 | −18.51 ± 0.032 | 87.27 ± ~0.09 | 1.0000 ± 0 |
| 15 | 74.52 ± 0.05 | 0.262 ± 0.0012 | 0.265 | 0.945 | −15.99 ± 0.039 | 85.51 ± ~0.09 | 0.9109 ± 0.0009 |
| 20 | 74.00 ± 0.05 | 0.297 ± 0.0013 | 0.299 | 0.947 | −13.03 ± 0.045 | 83.45 ± ~0.09 | 0.8498 ± 0.0011 |
| 30 | 72.87 ± 0.04 | 0.374 ± 0.0014 | 0.373 | 0.953 | −5.65 ± 0.057 | 78.31 ± ~0.09 | 0.6425 ± 0.0014 |
| 40 | 71.67 ± 0.04 | 0.458 ± 0.0015 | 0.452 | 0.960 | 3.72 ± 0.067 | 71.86 ± ~0.09 | 0.5216 ± 0.0015 |
| 60 | 69.37 ± 0.03 | 0.624 ± 0.0014 | 0.607 | 0.973 | 27.30 ± 0.080 | 56.50 ± ~0.09 | 0.3323 ± 0.0014 |
| 87.9167 | 67.07 ± 0.02 | 0.809 ± 0.0011 | 0.778 | 0.987 | 65.97 ± 0.078 | 34.88 ± ~0.09 | 0.1191 ± 0.0010 |

MC SEs for the selection-weighted functionals are exact delta-method SEs
from the multinomial selection counts,
√((Σ sel_c·x² − E²)/(draws·det_rate)); the naive-bias SE (marked ~) treats
the winner's per-draw noise scale as `se_g` and is approximate (formula in
the summary script header).

**Coherence, reported:** power monotone in q at all three thresholds ✓.
The zero-plus rung's declared rules carry **E[β(Ĥ)] = −22.29** — a benefit,
not a harm (the honest false-positive reading), though pulled up from the
pure background −26.98 by the selection mass on partially-Q-overlapping
rules (E[PPV] = 0.174); the task's "≈ −27" holds in sign and order, not to
the unit. All zero-plus selection mass is on rules with true mean below c1
(mass_below = 1.0000) ✓.

### §4.4 The calibration curve (grid read-off at c1 = T̂_obs)

| q | P(T ≥ T̂_obs \| q) ± SE |
|---|---|
| 0.01 | 0.0408 ± 0.0004 |
| 2.5 | 0.0432 ± 0.0005 |
| 5 | 0.0463 ± 0.0005 |
| 7.5 | 0.0499 ± 0.0005 |
| 10 | 0.0539 ± 0.0005 |
| 15 | 0.0635 ± 0.0005 |
| 20 | 0.0754 ± 0.0006 |
| 30 | 0.1101 ± 0.0007 |
| 40 | 0.1614 ± 0.0008 |
| 60 | 0.3222 ± 0.0010 |
| 87.9167 | 0.6333 ± 0.0011 |

Crossings by linear interpolation between adjacent rungs (stated as
interpolation): **0.05 at q ≈ 7.54; 0.5 at q ≈ 75.95; 0.95 not reached on
the ladder** (the top rung sits at 0.633). The observed statistic is
entirely consistent with planted effects from ≈ 7.5 CD4 units upward — the
curve is flat and low over the small-q half of the ladder because the
selection maximum, not the planted effect, carries most of T's mass there.

## §5 The `se_g` band — and the triggered sensitivity

Ratio `se_direct/se_g` over the 4508 members (`se_direct =
√(V_eff[g]/(n·Pg))`, per-member `fs_dgm_scale` regions):

| q | range | median | median PQg ≥ 0.95 (3 members) | median mid (118) | median PQg < 0.25 (4387) |
|---|---|---|---|---|---|
| 0.01 | [0.9659, 1.0413] | 1.0059 | 1.0014 | 1.0025 | 1.0061 |
| 20 | [0.9659, 1.0492] | 1.0099 | 1.0015 | 1.0119 | 1.0099 |
| 87.9167 | [0.9659, 1.0914] | 1.0290 | 1.0021 | 1.0619 | 1.0282 |

Top rung max |ratio − 1| = **9.14% > ~5% → the sensitivity ran** (top-rung
family with `se_g` replaced by `se_direct`, one 2e5-draw set, same seed —
labeled as the breadth ladder labeled it: **a sensitivity, not a
constructor, not adopted**):

| run | rate at c1 = 10 | power at 86 | power at 79 | P(T ≥ T̂_obs) | E\|Ĥ\| | E[β(Ĥ)] | E[PPV] |
|---|---|---|---|---|---|---|---|
| primary (`se_g`) | 0.9660 ± 0.0004 | 0.6604 ± 0.0011 | 0.7534 ± 0.0010 | 0.6333 ± 0.0011 | 67.07 | 65.97 | 0.809 |
| sensitivity (`se_direct`) | 0.9652 ± 0.0004 | 0.6742 ± 0.0010 | 0.7669 ± 0.0009 | 0.6473 ± 0.0011 | 67.30 | 63.48 | 0.787 |

Movement: rates ≤ 1.4 percentage points, E\|Ĥ\| +0.2 subjects, E[PPV]
−0.022, E[β(Ĥ)] −2.49 CD4 units (−3.8%) — small beside every §4 contrast
and immaterial to the conclusions.

## §6 Wall-clocks and concurrency

Eleven concurrent rung processes (launched together) + the sensitivity
(launched 15 min later after the band computed); host 128 cores, peak
memory ≈ 162 GB of 251 during the accumulation phase.

| rung (q) | grid secs | | rung (q) | grid secs |
|---|---|---|---|---|
| 1 (0.01) | 15 237 | | 7 (20) | 15 782 |
| 2 (2.5) | 15 304 | | 8 (30) | 15 548 |
| 3 (5) | 15 524 | | 9 (40) | 15 453 |
| 4 (7.5) | 15 382 | | 10 (60) | 15 195 |
| 5 (10) | 15 584 | | 11 (87.92) | 15 715 |
| 6 (15) | 15 732 | | sensitivity | 15 644 |

Total ≈ 51.7 CPU-h; realized wall-clock ≈ **4.7 h** (< the 6 h gate;
> the task's ~2 h estimate — see Deviations). §2–§3 build: 11 × ~44 s
enumerations, all gates in one process (~9 min). `se_g` band: ~6 s/rung.

## §7 The application document and the render check

`quarto/applications/actg175/analysis_actg175_continuous_oc_evaluation.qmd`
— reads `oc_applied_eval_summary_2026-08-31.rds` only, computes nothing
heavy, types no number (inline R throughout; the two gt/base-figure
sections per task §6, the orientation paragraph citing stage 0's 1.000,
the verbatim model-based caveat, the four adapted limitations). Payload
per the conventions: `table labels meta extras est_scale built_at
forestsearch_version`, `est_scale = "md"`, `results_dir`/`dirout` NULL →
`_payloads/<stem>/<stem>_payload.rds` (written, committed).

**Render check:** RStudio-bundled quarto 1.9.38,
`quarto render analysis_actg175_continuous_oc_evaluation.qmd` — **8.0 s
wall-clock** (< 2 min required ✓); `Output created:
analysis_actg175_continuous_oc_evaluation.html`; payload sentinel written
at render time (path echoed in the rendered HTML); no chunk recomputes a
draw set. HTML committed per the directory's convention (rendered HTML
tracked beside its qmd).

## Commits

- `fa3a7d34` tasks: the task document, alone (first action).
- `3ab3a6b0` sims: precompute scripts + per-rung/sensitivity/band/build/
  summary objects and logs (35 files; every `.rds` ≤ 744 KB, all tracked).
- `3118660b` applications: the qmd, its rendered HTML, its payload.
- (this file) reports: committed alone, after.

No push. No install. No `R/` change.

## Ten-line summary

1. All §1–§3 gates pass: clean tree, 0.3.2, anchor rebuilt 66/66, eleven
   rung families bit-identical off the planted effect, `beta_g` exactly
   linear in q through purity.
2. Eleven 2e5-draw grids ran concurrently; measured 4.2–4.4 h per rung
   (not 67 min — blocked-sweep cost); recomputed projection 4.7 h < 6 h,
   so the run continued and finished.
3. The analyst's gate (c1 = 10, c2 = 5) has **type-I 0.585** under the
   structural null; holding 5%/10% needs c1 = 86/79.
4. Power at c1_05 = 86 rises monotonically but reaches only 0.660 at the
   observed-statistic rung: **no rung on the ladder attains 80% power at
   fixed 5% type-I.**
5. What the analyst's gate declares is ~75 subjects with PPV 0.17–0.23 and
   a *beneficial* true effect (≈ −22 to −18) until q ≈ 35; PPV 0.81 and
   E[β(Ĥ)] ≈ 66 at the top rung; naive bias 90 → 35 across the ladder.
6. The calibration curve puts the observed T̂_obs = 87.92 at consonance
   0.041 under the structural null, crossing 0.05 near q ≈ 7.5 and 0.5
   near q ≈ 76; 0.95 is beyond the ladder.
7. The `se_g` band widens with q (max 9.14% at the top rung), so the
   triggered sensitivity ran: rates move ≤ 1.4 pp and E[β(Ĥ)] by −2.5
   units — immaterial to the conclusions, not adopted.
8. The application document reads the stored objects, renders in 8 s, and
   writes its payload under the conventions.
9. Committed: task doc alone, precompute + objects, document + payload;
   report alone follows. No push, no install, no `R/` change.
10. **The three numbers: type-I at the analyst's (10, 5) = 0.585; smallest
    q with 80% power at c1_05 = none on the ladder (0.660 at
    q = T̂_obs = 87.92); P(T ≥ 87.92 | q) at the rung nearest the 0.5
    crossing (q ≈ 76 → nearest rung 87.92) = 0.633.**

## Deviations

- **Wall-clock vs estimate:** the task's ≈ 67 min/rung (stage 0's 2e4-draw
  single-point cost × 10) omitted the blocked grid's per-threshold sweep
  (202 points × 4 blocks); measured ≈ 4.3 h/rung. The §CATEGORY gate was
  recomputed at first completion as instructed: projected ≈ 4.7 h < 6 h →
  continued. Total ≈ 51.7 CPU-h vs the projected 12.3.
- **Checkpoint size:** the saved named-point `fs_oc_predict` objects carry
  their family stripped of `memb`, `ovl`, `scale`, `cuts` (the ~250 MB
  carriers, reproducible from the committed build script); every numeric
  OC quantity, `lab`, and the family vectors are retained. Checkpoints are
  ~690 KB and tracked.
- **"Smallest integer c1 with rate ≥ 0.80/0.90/0.95"** read as the
  *largest* such c1 (nonincreasing curve; the breadth ladder's c1\*
  convention), stated in the document and above.
- **MC SEs for the declared functionals** are not package outputs; exact
  delta-method SEs from the multinomial selection counts were derived (and
  an approximate one for `Enaive_bias`), formulas in the summary script.
- **Coherence note:** zero-plus E[β(Ĥ)] = −22.29, not −27: the sign and
  reading the task expects (a benefit, the honest false positive), offset
  by selection tilt onto partially-Q-overlapping rules; reported, no gate.
- `T̂_obs` used at stored full precision (87.9166667) everywhere,
  including the c1 grid point, rather than the task's 6-dp literal.
- The eleven rungs were launched together and the zero-plus rung published
  `c1_05`/`c1_10` to the waiting rungs on completion, so the named-point
  predicts come from each rung's own in-memory grid (no re-draws).
- None otherwise: no `R/` change, no edit to any existing document, no
  push, no install, no `fs_oc_invert()` runs.
