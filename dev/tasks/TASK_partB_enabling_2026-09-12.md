# TASK — Part B enabling change, and the OC smoke

- **Date:** 2026-09-12
- **Machine:** Mac Studio. **Branch:** `feature/glm-extension`, current tip.
- **Purpose:** make the identification-only sweep runnable, then run a 30-replicate smoke across every engine × criterion and produce the operating-characteristics table, so Larry can confirm what Part B will produce before it runs.
- **Commit only, do not push. No `R/` change** — everything here is the simulation template and its recorder. If anything appears to need an `R/` change, STOP and report.

## Background — settled, do not re-determine

- Part B records **identification and classification only**. MR is not run, and no MR quantity is wanted.
- `mr_inference = FALSE` is already the `forestsearch()` default (`R/forestsearch_main.R:1298`) and is gated on all three engines (`:3320`, `:2299`, `:2490`). The template hard-codes `TRUE` at line 1017.
- Part B is **288 cell-runs** — FS 108, DINA 90, GRF 90 — across six criteria × 18 cells. On DINA and GRF, `maxeff` and `maxeffCons` are the same run.
- The six criteria: `effMaxSG` (ε 0.20), `maxeffCons`, `effMinSG` (ε 0.20), `maxSG`, `minSG`, `maxeff`. ε is consulted only by the two band foci; supplying it elsewhere is an error.

## Part 1 — three changes

### 1a. Turn MR off

- Make `mr_inference` settable from the template, defaulting to the current `TRUE` so nothing changes when it is unset. Implement it however fits the template's existing knob pattern — an `FS_S7_*` variable is the obvious route, but the route is yours.
- Echo it in the knob audit line and the settings readout; record it in batch and combined meta.

### 1b. Widen the focus guard

- The guard at template line 327 admits four criteria. Widen it to all six: `effMaxSG`, `maxeffCons`, `effMinSG`, `maxSG`, `minSG`, `maxeff`. Both additions are valid `sg_focus` values (`forestsearch_main.R:566–578`; `fs_focus_tag.R:64–88`).
- Keep it a guard — reject anything outside that set.
- Enforce what the source already requires: `effect_neighborhood` and a non-default `selection_rule` are accepted only for `effMaxSG` and `effMinSG`, and are an error otherwise.

### 1c. The recorder

- With MR off, `n_sel`, `label`, `n_family` and `band_n` come back NA because the recorder fills them from the MR object.
- Route them from the `forestsearch()` result instead, which carries the identification regardless of MR. **Take whichever route is available**; name the fields you read and quote the lines.
- MR-derived columns stay exactly as they are — NA when MR is off, by design.
- If one of the four genuinely cannot be recovered from the result without an `R/` change, record which and why, leave it NA, and carry on. Do not make an `R/` change and do not work around it.

### Gate T3 — stop-on-failure

- **Unset:** 5 replicates at the standing identity cell with the new knob unset — every non-timing column and truth `identical()` to a pre-change run on this machine. Revert and STOP on failure.
- **MR off:** 5 replicates at the same cell with MR off — `n_sel`, `label`, `n_family`, `band_n`, `detected` and the classification columns present and populated, and **identical to the unset run on the same seeds**. MR cannot change which subgroup is identified (`forestsearch_main.R:1041–1044`); this gate tests that. Report any difference and STOP.

## Part 2 — the OC smoke

Runs only if Gate T3 passes. **Cap 90 minutes wall**; if the projection exceeds it, cut the replicate count and say so.

- **Coordinate:** 12.4% (`FS_S7_Z1Q` unset), HR 1.50, n 500.
- **30 replicates**, MR off, 12 workers with the three thread variables at 1. Tag so nothing can collide with a campaign glob.
- **Sixteen runs** — the full engine × criterion structure at one cell:

| Engine | Criteria | Runs |
|---|---|---|
| consistency (FS) | all six | 6 |
| DINA | `effMaxSG`, `effMinSG`, `maxSG`, `minSG`, `eff` | 5 |
| GRF | `effMaxSG`, `effMinSG`, `maxSG`, `minSG`, `eff` | 5 |

- ε 0.20 for `effMaxSG` and `effMinSG` only. On GRF the rule is carried by `frontier_rule` with `grf_selection = "frontier"`, and its ε default is 0.10, so set 0.20 explicitly.
- On DINA and GRF, `maxeff` and `maxeffCons` both resolve to `eff` — run that once, and say in the table that it stands for all three spellings there.

### The OC table

One row per run. This is what Larry is confirming, so print it in full:

- engine, `sg_focus`, ε (or "inert")
- detection or selection rate, with a Wilson interval
- sensitivity, specificity, PPV, NPV
- mean |Ĥ| and mean |H|, and their ratio
- `n_family` median and range
- `admitted_n` median and range (GRF)
- realized trial prevalence
- median per-replicate wall, and the run's total wall

Then, below the table:

- any column that came back NA and why;
- whether any criterion selected the same subgroup as another on the same seeds — if two rules never differ at this cell, say so, since it bears on what the sweep will show;
- the projected cost of the full 288 cell-runs at 2,000, 1,000 and 500 replicates per cell, from the measured per-replicate wall, broken out by engine.

**Then stop.** Run no sweep cell. The replicate count is Larry's decision.

## Closeout

- Commit everything per the standing rule in `scripts_dinamr/README.md`; report path, tracked status and size.
- Write the findings to `REPORT_partB_enabling_2026-09-12.md`.
- Last action: regenerate `current_status.md` with the pin post-condition (`check_current_status.sh`).
- Give the commit range to push.

## Report

One paragraph: whether Gate T3 passed both halves; what the OC table shows; the projected totals; and the commit range. Bullet form for the detail, one item per bullet. Print the OC table in full.
