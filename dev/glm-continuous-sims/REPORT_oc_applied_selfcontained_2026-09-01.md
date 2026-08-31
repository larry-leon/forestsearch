# REPORT — the self-contained applied OC document: wrapper-forward, c2 = 0.8·c1, the Q_variants knob (2026-09-01)

**Task:** `dev/tasks/cc_task_oc_applied_selfcontained_2026-09-01.md` (amended
version, commit `cf4e85ba`, superseding `d3f11f2a`)
**Verdict up front: all gates PASS; the document renders self-contained in
26.7 min at twelve-wide and writes its payload; the archived evaluation got
its two-sentence pointer and was not re-rendered.**

## Provenance

- Host `pop-os` · `/home/larryleon/Documents/GitHub/forestsearch` · branch
  `feature/glm-extension` · installed version **0.3.2** ✓ · commits
  `209a8e85` (stage 0), `eb136a35` (0.3.2), `d597f289` (the evaluation
  report) in the log ✓.
- **Transport names.** The original task arrived hyphen-stripped as
  `cc_task_oc_applied_selfcontained_20260901.md` (committed `d3f11f2a`);
  the amended version arrived as
  `cc_task_oc_applied_selfcontained_20260901 (1).md` — identified by its
  `Q_variants` content and used; both copies remain in `~/Downloads`.
  Committed in-repo as `cc_task_oc_applied_selfcontained_2026-09-01.md`,
  alone, as `cf4e85ba` (first action of the amended run).
- §1 dirt check at the amended run's start: exactly two paths, both this
  task's own and both explicitly carried over from the prior session by
  instruction — the drafted
  `quarto/applications/actg175/analysis_actg175_continuous_oc.qmd`
  (untracked) and the pointer note in
  `analysis_actg175_continuous_oc_evaluation.qmd` (modified).  Nothing
  else.  Per instruction, kept, updated, not reverted.

## §2 gates — all PASS (enforced as `stopifnot` inside the rendered document; render exit 0)

- **Anchor** (in-document `forestsearch()`, sequential, ~8 s):
  Ĥ = `{age <= 37} & !{cd40 <= 507}` (clause set, order-insensitive),
  n(Ĥ) = 66, T̂_obs = 87.916667 (= 87 + 11/12 within 5e-7),
  p.consistency = 0.95 ✓.
- **Family**: M = 4508 with `orientation$s = +1` at the anchor build and at
  every one of the twelve (variant, rung) jobs ✓.  Q representability:
  nearest member `age <= 37 & cd40 > 506`, purity 0.9969, not forced ✓.
- **Q_variants nesting** (strict, on `df_super`): primary
  `age ≤ 37 & cd40 > 507` (P = 0.0634) ⊂ wider `age ≤ 39 & cd40 > 449`
  (P = 0.1240, **1.96×**) ⊂ broad `age ≤ 39 & cd40 > 412` (P = 0.1906,
  **3.01×**) ✓ — both thresholds relaxed vs the primary, chosen from the
  family's own population quantile grid (age grid contains 37 and 39;
  cd40 grid 204.5…506) nearest 2×/3× the primary prevalence.
- **Shared structural null**: per-variant
  max |beta_g(0.01) − beta_treat·(1 − PQg)| = 0.009969 / 0.010000 /
  0.010000 ≤ 0.01 ✓ (wider/broad sit at exactly 0.01·1: each has an
  on-grid pure member).  The primary's zero-plus grid stands as the null
  for all three; it was not re-run.
- **Self-containment**: `grep 'dev/'` on the qmd — zero hits; no `readRDS`
  of anything the document did not itself write; every number in prose is
  inline R ✓.
- **§3 render gates**: zero-plus diagonal type-I at c1 = 10 finite
  (0.5042); its declared-rule row oriented-negative
  (E[β(Ĥ)] = −22.29) ✓; ladder power columns monotone in q within 2 MC
  SEs ✓; payload written ✓; rendered HTML free of `NA ±` and error text
  (0 grep hits each) ✓.

## Compute gate and render

- **Probe** (scratch call of the job chunk code, one job, 156-point grid at
  2e4 draws): 604.5 s ≈ 10.1 min — enumerate 45.6 s + draw 410.5 s +
  sweep 147.4 s — peak ≈ 5.1 GB.  Twelve-wide projection ≈ 22–25 min
  **≪ 2.5 h → GO.**  (The task's ~6 s/grid-point estimate was for the
  blocked path; at draws = 2e4 < block = 5e4 the exact one-block
  reduction runs ≈ 0.94 s/point.)
- **Render** (RStudio-bundled quarto 1.9.38): exit 0,
  **wall-clock 1599 s = 26.7 min** (inside the expected 35–55 min band's
  floor); evaluation loop alone 1091.9 s = 18.2 min over the twelve jobs
  at `n_workers = 12`.
- **Memory**: host used-memory sampled every 15 s — baseline 32.7 GB, peak
  103.3 GB → **≈ 70.6 GB attributable to the render** (twelve concurrent
  jobs at ~5 GB each plus the parent; the task's ≈ 10 GB estimate appears
  to have carried the earlier six-rung arithmetic).  Host: 128 cores,
  251 GB — no pressure.

## The diagonal headline, beside the archived fixed-c2 run

| quantity | this document (diagonal c2 = 0.8·c1, 2e4 draws) | archived (fixed c2 = 5, 2e5 draws) |
|---|---|---|
| type-I at the analyst's gate | **0.5042 ± 0.0035** at (10, 8) | 0.5850 ± 0.0011 at (10, 5) |
| fixed-type-I 0.05 threshold | **c1_05^diag = 50** (rate 0.0273 ± 0.0012, ladder resolution) | c1_05 = 86 (integer resolution) |
| fixed-type-I 0.10 threshold | c1_10^diag = 40 (rate 0.0657 ± 0.0018) | c1_10 = 79 |
| power at the 0.05 threshold, top rung | 0.5766 ± 0.0035 at (50, 40) | 0.6604 ± 0.0011 at (86, 5) |
| calibration read-off point | (c1 = T̂_obs, c2 = 0.8·T̂_obs ≈ 70.3) | (c1 = T̂_obs, c2 = 5) |
| P(T ≥ T̂_obs) crosses 0.05 at | q ≈ 61.1 | q ≈ 7.54 |
| P(T ≥ T̂_obs) crosses 0.5 at | beyond the ladder (0.189 at the top rung) | q ≈ 75.95 |
| comparator: zero-plus rate at (10, 5) | 0.5884 ± 0.0035 | 0.5850 ± 0.0011 — agrees to MC precision |

**What the c2 policy changed.**  Scaling the consistency floor with the
screening floor moves the binding screen: for noisy candidates
0.8·c1 + z_p·se_g overtakes c1, so the consistency gate does the excluding
that the fixed c2 = 5 left entirely to the effect screen.  Five-percent
type-I control arrives at c1 = 50 on the diagonal versus 86 on the fixed
row — a usable threshold rather than one sitting at the observed statistic
itself — at the cost that any read-off *at* the observed statistic becomes
far more conservative: with c2 ≈ 70.3 accompanying c1 = T̂_obs, reaching
the observed 87.92 requires clearing a consistency floor the archived run
never imposed, so the calibration crossings shift far upward (0.05 at
q ≈ 61 vs ≈ 7.5) and the 0.5 crossing leaves the ladder.  The declaration
event and the declared-rule quality at the analyst's own gate barely move
(E|Ĥ| ≈ 76, E[PPV] 0.17, E[β(Ĥ)] ≈ −22 at the zero-plus rung under either
c2) — the policy's effect is concentrated where it should be, on the
thresholds that control error rates, not on what a liberal gate declares.

**The knob** (new §10 of the document): power at (50, 40) at the top rung
rises 0.577 → 0.856 → 0.951 across primary → wider → broad while
sensitivity of the declared rule against the variant's own Q falls
0.779 → 0.499 → 0.343 (PPV 0.811 → 0.861 → 0.870); the 0.05 calibration
crossing falls q ≈ 61.1 → 43.8 → 42.2 — the observed statistic constrains
a trade between breadth and severity, as the amended task's reading
predicted the table would have to say on its own numbers.

## Commits

- `cf4e85ba` tasks: the amended task document, alone (first action;
  supersedes `d3f11f2a`).
- `c73a6b06` applications: the new qmd + its rendered HTML + its payload
  (`_payloads/analysis_actg175_continuous_oc/…_payload.rds`) + the
  two-sentence pointer edit to
  `analysis_actg175_continuous_oc_evaluation.qmd` (its HTML **not**
  re-rendered).
- (this file) reports: committed alone, after.

No push.  No install.  No `R/` change.

## Ten-line summary

1. Amended task (identified by `Q_variants` content among two `~/Downloads`
   copies) committed alone as `cf4e85ba`, superseding `d3f11f2a`; all
   provenance gates pass on 0.3.2.
2. The carried-over draft qmd was updated in place to the amended spec; the
   pointer note stood as-is; neither was reverted.
3. Compute probe: one (variant, rung) job = 10.1 min at 2e4 draws —
   twelve-wide projection ≈ 25 min ≪ 2.5 h → rendered.
4. The document is self-contained: anchor `forestsearch()` (8 s, gated
   66/66, T̂_obs = 87.916667), `generate_glm_dgm()`/`fs_dgm_scale()`,
   `fs_oc_family_enumerate()` (M = 4508), one visible `fs_oc_predict()`,
   and a ~15-line `mclapply` loop of twelve `fs_oc_grid()` jobs; zero
   `dev/` references.
5. Render: 26.7 min wall-clock, ≈ 70.6 GB peak, exit 0; payload written;
   HTML free of `NA ±` and error text.
6. Diagonal headline: type-I 0.5042 at (10, 8); c1_05^diag = 50 and
   c1_10^diag = 40 at ladder resolution, against the archived 86/79 at
   fixed c2 = 5.
7. The comparator row reconciles: zero-plus 0.5884 ± 0.0035 at (10, 5) vs
   the archived 0.5850 ± 0.0011.
8. Calibration at the diagonal read-off is far more conservative: 0.05
   crossing at q ≈ 61 (archived ≈ 7.5), 0.5 beyond the ladder (archived
   ≈ 76).
9. The Q_variants knob (1×/1.96×/3.01× prevalence, strictly nested,
   grid-chosen; shared null verified at ≤ 0.01) shows the breadth–severity
   trade: top-rung power at the 0.05 gate 0.58 → 0.95, sensitivity
   0.78 → 0.34, 0.05-crossing q 61 → 42.
10. Commits `cf4e85ba` (task), `c73a6b06` (document set), this report
    alone; no push, no install, no `R/` change.

## Deviations

- **Knob placement**: the `Q_variants` list is defined in the document's
  §3 knob chunk, not the setup chunk — the superset thresholds need the
  population frame (`df_super`) that exists only after the first DGM
  build; the rung and ladder vectors are setup-chunk defaults as
  specified.
- **Knob table scoring rung**: the task fixed no rung for the PPV/sens
  columns at (10, 8); the document scores them at the observed-statistic
  rung and its table note says so.
- **Compute estimates vs measured**: sweep ≈ 0.94 s/point (one-block exact
  path at draws < block), not ~6 s; job 10.1 min, not 22–25; render
  26.7 min (below the 35–55 band); peak memory ≈ 70.6 GB, not ≈ 10 GB.
  All on the favorable side except memory, which fit the host with wide
  margin.
- **Calibration 0.5 crossings** for the primary and wider variants lie
  beyond the ladder at the diagonal read-off point (top-rung values 0.189
  and 0.384; the broad variant crosses at q ≈ 86.5, top-rung 0.515); the
  document prints "beyond the ladder" rather than a number, by
  construction.
- **`null_dev` for wider/broad = 0.01 exactly** (each superset has a pure
  on-grid family member, so the bound 0.01·max(PQg) is attained); the gate
  compares with a 1e-9 floating-point guard and passes.
- The probe ran as a scratch call in the session scratchpad (task-allowed
  alternative to a probe render); a first probe launch from the
  pre-amendment session was killed mid-grid on the user's stop and re-run
  to completion under the amended task.
- None otherwise: no `R/` change, no re-render of the archived evaluation,
  no other document touched, no 2e5-draw run, no gbsg, no push.
