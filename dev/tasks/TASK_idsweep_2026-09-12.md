# TASK — Part B: the identification sweep (`idsweep`)

- **Date:** 2026-09-12
- **Machine:** Mac Studio. **Branch:** `feature/glm-extension`, current tip.
- **Campaign tag:** `idsweep`. Report-and-wait; proceed unattended.
- **Commit only, do not push. No `R/` change and no template change** — the enabling change landed in `0ab5d1c5`. If anything appears to need a code change, STOP and report.

## What this is

- **288 cell-runs**: three engines × six selection criteria × 18 cells, **500 replicates each**, **MR off**.
- It measures **identification and classification only**. No MR, no field constructions, no coverage, no bias, no bounds. The MR-derived columns are NA by design and are not reported.
- **A fresh, self-contained campaign.** Do not merge, reuse or cite any earlier bundle as part of this sweep, even where a committed bundle covers the same engine, criterion and cell. Everything in `idsweep` is run here, under one tag and one knob set.

## The grid

**18 cells**, both prevalences × three hazard ratios × three sample sizes:

| Prevalence | HR | n |
|---|---|---|
| 12.4% (`FS_S7_Z1Q` unset) | 1.50, 1.75, 1.00 | 500, 1000, 1500 |
| 31% (`FS_S7_Z1Q=0.60`) | 1.50, 1.75, 1.00 | 500, 1000, 1500 |

**16 runs per cell** — the engine × criterion structure verified in the smoke:

| Engine | Criteria | Runs |
|---|---|---|
| consistency | `effMaxSG` (ε 0.20), `effMinSG` (ε 0.20), `maxeffCons`, `maxeff`, `maxSG`, `minSG` | 6 |
| DINA | `effMaxSG` (ε 0.20), `effMinSG` (ε 0.20), `maxSG`, `minSG`, `eff` | 5 |
| GRF | `effMaxSG` (ε 0.20), `effMinSG` (ε 0.20), `maxSG`, `minSG`, `eff` | 5 |

- On DINA and GRF, `eff` stands for `eff` = `maxeff` = `maxeffCons` — run once, and say so in every table.
- ε is set only for `effMaxSG` and `effMinSG`; supplying it elsewhere is an error under the widened guard.
- On GRF the rule is carried by `frontier_rule` with `grf_selection = "frontier"`, and its ε default is 0.10 — set 0.20 explicitly.

## Knobs

As `campaign.sh` pins them, plus:

- **MR off** — the knob added in `0ab5d1c5`. Confirm it resolves to off on every batch; `_nomr` should appear in every output filename.
- `FS_S7_CAMPAIGN=idsweep`; `FS_S7_ER_JCUTS` unset.
- Seeds 8316951 + sim_id, sim_id 1–500, **one batch per run** (500 replicates does not need splitting; say so in the record).
- 12 workers with the three thread variables at 1.

## Run order — defer from the tail

All 16 runs complete within a cell before the next cell starts, so a deferral costs cells rather than criteria.

1. 12.4%, HR 1.50, n 500 / 1000 / 1500
2. 31%, HR 1.50, n 500 / 1000 / 1500
3. 12.4%, HR 1.75, n 500 / 1000 / 1500
4. 31%, HR 1.75, n 500 / 1000 / 1500
5. 12.4%, HR 1.00, n 500 / 1000 / 1500
6. 31%, HR 1.00, n 500 / 1000 / 1500

## Gate 1 — compute go/no-go

- **Ceiling 13 h wall; hard timeout 16 h.** Record both.
- **Reconcile the projection before trusting it.** `REPORT_partB_enabling_2026-09-12.md` gives an adjusted total of 11.1 h at 500 replicates, but applying its own per-engine multipliers (×2.08 consistency, ×2.33 DINA, ×1.26 GRF) to its per-engine flat figures (3.2 / 1.4 / 2.9 h) gives about 13.6 h. State which is right and why, and project from the corrected basis.
- Project per criterion from the smoke's measured per-replicate costs — `consistency maxeff` ran +63% over the other criteria, so do not assume flatness there.
- Render overhead was measured at 16.4 s on 30-replicate renders only. A 500-replicate MR-off render is unmeasured; measure it on the first cell and re-project.
- If the 288 do not fit, defer from the tail and record it. **Do not reduce replicates** — 500 is Larry's decision.

## Gate A — alignment, per run, stop-on-failure

- MR resolves to off; `_nomr` in the filename.
- `sg_focus` resolves to the intended value; ε is 0.20 for the two band foci and unset otherwise.
- On GRF: `grf_selection = "frontier"`, `frontier_rule` as intended, ε 0.20.
- Report all resolved values.

## Gate I — integrity, per cell, stop-on-failure

- The identification and classification columns are present and populated on every detected replicate: `detected`, `n_sel`, `label`, `sens`, `spec`, `ppv`, `npv`, and `admitted_n` on GRF.
- No errors; non-detections recorded as such and counted.
- Realized trial prevalence against the super-population value.
- `n_family`, and `n_cons_qual` / `band_n` on DINA and GRF, are NA by design — report as structural, never as failures.
- **Same-draws check:** where a committed bundle exists at the same coordinate, `n_true` `identical()` on the shared sim_id range 1–500. A mismatch is a finding about the DGM path, not a cell failure. This is a check on the draws only; **no committed bundle is used as data in this sweep.**

## Stage 3 — the summary

- Build `summary_idsweep.qmd` from the smoke's OC table script, extended to the full grid. **Do not transplant `summary_dinamr.qmd` or `summary_grfmr.qmd`** — they are built around MR products that are NA here.
- Guard per-cell chunks so absent cells skip rather than render empty.
- **Per cell**, the OC table exactly as the smoke produced it: engine, `sg_focus`, ε, detection or selection rate with a Wilson interval, sensitivity, specificity, PPV, NPV, mean |Ĥ| and mean |H| and their ratio, `admitted_n` median and range (GRF), realized prevalence.
- **Across cells**, one table or figure each:
  - detection or selection rate against n, by engine and criterion, both prevalences;
  - mean |Ĥ| / mean |H| against n, by engine and criterion;
  - sensitivity, specificity, PPV and NPV against n, by engine and criterion;
  - criterion agreement — for each engine and cell, the share of detected replicates on which each pair of criteria selected the identical subgroup. The smoke found `maxeffCons`/`maxeff` agreeing on 27 of 29 and `effMinSG`/`minSG` on 26 of 29 on the consistency engine, far less on DINA and GRF; the sweep should show whether that holds across the grid.
- At the HR 1.00 cells: call the rate a **selection rate**, use no harm vocabulary — "differentially null against a benefiting complement" — and note that the region clears the sub-null log(0.90) floor so returning it is admissible.
- Wilson intervals on every rate. No coverage, no bias, no bounds, no MR quantity anywhere.
- No acceptance criterion and no recommendation.

## Closeout

- Commit everything per the standing rule in `scripts_dinamr/README.md`; report path, tracked status and size; flag anything over 50 MB, hard stop at 100 MB.
- Write the findings to `REPORT_idsweep_2026-09-12.md`.
- Last action: regenerate `current_status.md` with the pin post-condition (`check_current_status.sh`), adding `idsweep` to §2 and updating §7.
- Give the commit range to push.

## Report

One paragraph: cells completed, deferred and dropped with walls; the headline on how detection, subgroup size and classification move with the criterion and with n; the Gate 1 projection against realized; and the commit range. Bullet form for the detail, one item per bullet.
