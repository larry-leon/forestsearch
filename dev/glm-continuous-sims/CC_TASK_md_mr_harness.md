# CC task: GLM/continuous MR harness — harm direction

Repo `forestsearch`, branch `feature/mr-in-replicates`. This is the priority
workstream. If the regeneration campaign is running, it is background-only and
never blocks or interleaves with this task.

**Decisions, made and final — do not reopen:**
- `focus = "harm"` everywhere.
- DGM: the verified harm calibration — `calibrate_glm_interaction(outcome_type
  = "continuous", effect_measure = "MD")` with `cal_target_md = -40`, all else
  as the continuous template. Convergence at machine precision is already
  established; do **not** recalibrate, and do not "fix" the selection rate —
  rate ≈ 1.0 at this cell is the phenomenon, on the record.
- Search orientation: `adverse_outcome = FALSE` with thresholds per the
  committed `dev/glm-continuous-sims/NOTE_threshold_sign_md.md` — the same
  corrected configuration as the 100-replicate harm sweep (the one that
  produced mean selection effect −61.87). Reuse that configuration verbatim.
- Aim: unbiasedness and interval coverage of the MR-corrected estimate against
  the **exact** β(Ĥ), with θ† reference columns. One harness, both readouts.
- Identifier scope v1: the FS/consistency configuration from that sweep only.
  DINA/GRF cells are a follow-on (GRF is now safe via the resolver, but out of
  scope here).

## The document

`quarto/simulations/actg175/continuous/mr_coverage_sweep_md_harm.qmd`,
mirroring the structure of the binary MR coverage/sweep documents (same
`run_cell()` / bundle / meta conventions, same `B` for the multiplier
resampling as the committed binary documents — read theirs, don't invent).
Swap: the harm-calibrated continuous DGM; scoring frame
`fs_build_eval_frame(dgm, "continuous")` (returns `df_super`; the target is
exact); `fs_attach_betaHhat(..., focus = "harm", outcome_type = "continuous",
effect_measure = "MD")`; θ† via `fs_betaHhat_theta_dagger_check()`.

Seeds pre-generated and indexed by **global replicate id**, not within-batch
position — results must be invariant to worker count. `devtools::install()`
before any parallel run. Uniform bundle meta signature plus: harm direction,
`cal_target_md`, threshold orientation, package commit.

## Per-cell readouts

For each cell, from the same replicates:
1. **Bias**: mean and SE of (estimate − β(Ĥ)) for the naive in-trial MD and
   for the MR-corrected estimate, side by side — the comparison is the
   unbiasedness story, and every deviation is the method, since the target
   carries zero Monte Carlo error.
2. **Coverage**: MR interval coverage of β(Ĥ), with `n_eff` from the parity
   guard (runs automatically; a parity failure stops the run).
3. **Selection**: `selection_rate`, and the selection-effect distribution
   (in-trial MD − β(Ĥ)): mean, SD, quartiles.
4. **Identification** vs Q: sensitivity, specificity, PPV, **NPV** (project
   convention — NPV, not Jaccard).
5. Partition from the bundle: `nH_eval + nHc_eval == N` on every row,
   including undetected rows (which now carry `0 + N` with the ITT
   complement — those rows are part of the Ĥᶜ coverage denominator, not
   holes).
6. θ† reference columns beside β(Ĥ), clearly labeled as distinct estimands.

## Phasing

**STOP A — pilot.** One cell: `n = 1000`, 50 replicates, MR on, full scoring
path. Report every readout above verbatim, the per-replicate wall-clock
including MR, and a proposed full-grid schedule. Grid to propose:
`n ∈ {500, 1000, 2000, 4000}`, `s = 1000` replicates per cell — adjust only
if the measured cost says otherwise, and say so.

**STOP B — full grid**, on approval. All cells, PID-watched (never `pgrep -f`
a string the watcher's own command line contains). Deliver the four-cell
tables for readouts 1–4, the bias-vs-n and coverage-vs-n trends, and the
document's own rendered tables. Hold for review; commit only on approval —
document plus bundles, provenance meta complete.

## Constraints

- No `R/` changes. No DGM recalibration. No changes to the benefit vignette,
  the acceptance document, shims, or frozen directories.
- New files only: the document, its results directory, and (if the binary
  pattern uses one) a small runner.
- A guard failure (parity, partition) stops the phase and is reported
  verbatim — not reconciled in place.
- Disagreements with this document's decisions are reported, not resolved:
  the decisions are closed.

## After STOP B (for the record, not for action now)

Binary is next and is this same document with `outcome_type = "binary"`,
`effect_measure = "OR"`, and its own committed DGM calibration — the port is
a parameter change, which is the point of the consolidated infrastructure.
