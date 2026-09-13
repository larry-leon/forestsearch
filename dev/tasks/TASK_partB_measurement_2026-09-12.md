# TASK — Part B, initial measurement: what MR costs, and what the six criteria resolve to

- **Date:** 2026-09-12
- **Machine:** Mac Studio (the committed comparison bundles are Mac-produced, so walls must be measured against them there).
- **Branch:** `feature/glm-extension`, the current tip.
- **Purpose:** decide Part B's replicate count and machine allocation from measured numbers rather than assumption. **This task does not run the sweep.** It ends with a report and stops.
- **Commit only, do not push. No `R/` change.** A template change may be needed (Stage 0c) — if so, report it and stop; do not make it.

## Background, so nothing is re-derived

- **Part A** (not this task): FS at `effMaxSG` ε 0.20, the nine missing 12.4% cells, full MR. Takes FS to 18/18 and matches DINA and GRF.
- **Part B** (this task measures for it): identification only, no MR, across six `sg_focus` criteria × 18 cells × 3 identifiers.
- `mr_inference = FALSE` is the default in `forestsearch()` (`R/forestsearch_main.R:1298`) and is gated on every engine — `.mr_eligible` (`:3320`), `.mr_dina_ok` (`:2299`), `.mr_grf_ok` (`:2490`). With `FALSE`, `out$mr_inference <- NULL` and `.fs_apply_mr()` is never reached. The documentation states that `TRUE` and `FALSE` give identical `sg.harm`, `df.est` and `max_sg_est` (`:1041–1044`).
- The band arguments `selection_rule` and `effect_neighborhood` are consulted **only** by `effMaxSG` / `effMinSG` (`:602–606`); supplying a non-default `selection_rule` with a non-band focus is an error, not a no-op.

## Stage 0 — verify from source, no compute

Quote every finding with file and line numbers.

### 0a. The resolution map at HEAD

Confirm `fs_focus_tag()` (`R/fs_focus_tag.R:64–88`) resolves as follows, and report any departure:

| Spelling | consistency | dina / grf |
|---|---|---|
| `eff`, `hr`, `maxcons` | `maxcons` | `eff` |
| `maxeff` | `maxeff` | `eff` |
| `maxeffCons` | `maxeffCons` | `eff` |
| `effMaxSG`, `hrMaxSG` | `effMaxSG` | `effMaxSG` |
| `effMinSG`, `hrMinSG` | `effMinSG` | `effMinSG` |
| `maxSG`, `minSG` | pass through | pass through |

- State plainly whether, on DINA and GRF, `maxeffCons` and `maxeff` are the same run. If they are, Part B is 6 + 5 + 5 = **288** cell-runs rather than 324.
- Confirm that the tag is what drives the run, not merely the output filename — i.e. that two spellings sharing a tag genuinely execute the same selection.

### 0b. ε scope

- Confirm `effect_neighborhood` is consulted only by `effMaxSG` / `effMinSG` on **all three** engines (`forestsearch_main.R:602–606`, `dina_subgroup.R:214–222`, `grf_main.R:47–52`).
- Report GRF's default (`grf_main.R:51` states 0.10), and confirm 0.20 must be set explicitly for GRF's band rules.
- **Check against the committed bundles:** the 12.4% FS comparators (`tier2`, `p12ext`) ran `maxeffCons` with `effect_neighborhood = 0.10`. If `maxeffCons` does not consult the band, that setting was **inert** on those runs, and those nine cells differ from the 31% ones in the *rule*, not in ε. Confirm or refute from the bundle metas. Every report in this directory currently says "different ε"; if that is wrong it is a records correction, and you should say so without editing the records in this task.

### 0c. Is `mr_inference` reachable from the template?

- Determine whether the simulation template exposes `mr_inference` as an `FS_S7_*` knob, or hard-codes it.
- If a knob exists, name it and its default. If not, say what an add-only, default-inert template change would touch — **do not make it**.

## Stage M — the measurement

**Coordinate:** 31% prevalence (`FS_S7_Z1Q=0.60`), HR 1.50, n 500. All three engines have a committed bundle there at `effMaxSG` ε 0.20, so each measurement has an exact comparator.

**Runs:** three, one per engine (`consistency`, `dina`, `grf`), **200 replicates each**, `mr_inference = FALSE`, everything else exactly as the committed bundle at that coordinate — same seeds (8316951 + sim_id), same knobs, 12 workers with the three thread variables at 1. Tag them so they cannot collide with a campaign glob.

**If Stage 0c finds no knob:** stop and report. Do not proceed to Stage M by other means.

**Compare per-replicate compute, not wall.** The committed bundles record `fit_mr_secs` per replicate, so the comparison does not require matching replicate counts. Report per engine:

- median, q25, q75 and p90 of `fit_mr_secs`, MR off, over the 200 replicates;
- the same quantiles from the committed bundle at that coordinate, MR on, over its 2,000;
- **the MR share of per-replicate cost**: 1 − (median MR-off / median MR-on);
- whether the identification output is unchanged — `sg_def`, `n_sel`, `detected` and the classification columns identical to the committed bundle on the shared `sim_id` range, which the documented guarantee predicts. **Report any difference; do not fix it.**
- bundle size on disk, MR off against MR on, per replicate.

## Gate P — project Part B, then stop

From the measured MR-off per-replicate cost, project the full sweep. Assume the cost is roughly flat across `sg_focus` within an engine unless the measurement suggests otherwise, and say so if it does.

Report a table of total compute for **288 cell-runs** (108 FS, 90 DINA, 90 GRF) at:

| Replicates per cell | Projected total |
|---|---|
| 2,000 | |
| 1,000 | |
| 500 | |

- Break it down per engine, since the MR share will differ between them.
- State the Monte Carlo resolution at each count: Wilson half-width on a rate near 0.50 and near 0.90, and the detected-replicate count at the thinnest cell in the grid (DINA, 12.4%, HR 1.00, n 1500, detection 0.344).
- **Then stop and report.** Do not run any sweep cell. The replicate count and the machine allocation are Larry's decisions and depend on these numbers.

## Closeout

- Commit the measurement bundles, the driver and the report per the standing rule in `scripts_dinamr/README.md`; report path, tracked status and size.
- Write the findings to `REPORT_partB_measurement_2026-09-12.md`.
- Last action: regenerate `current_status.md` with the pin post-condition (`check_current_status.sh`).
- Give the commit range to push.

## Report

One paragraph: the three Stage 0 determinations; the MR share per engine; the projected totals at 2,000 / 1,000 / 500; and the commit range. Bullet form for the detail, one item per bullet.
