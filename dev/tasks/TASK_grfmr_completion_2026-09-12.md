# TASK — `grfmr` completion: the two deferred harm cells and the six HR 1.00 cells

- **Date:** 2026-09-12
- **Campaign:** `grfmr` (continuation). Report-and-wait; proceed unattended.
- **Machine:** Mac Studio. Fresh CC session.
- **Predecessors, committed:** `dev/tasks/TASK_grfmr_campaign_2026-09-11.md`, `dev/tasks/TASK_dinamr_blockC_grfprobe_2026-09-11.md`, `dev/tasks/TASK_dinamr_campaign_2026-09-10.md`. Read them for the DGM, cell definitions, Gate 2 content and the structurally-NA column list. This document governs where they differ.
- **Tooling:** reuse `quarto/simulations/gbsg_020/scripts_dinamr/` — `campaign.sh`, `render.sh`, `gate2G.R`, `projectC.R`, `grfmr_numbers.R`, `grfmr_tables.R`, `check_current_status.sh`. See that directory's `README.md` for the two standing rules.
- **Commit only, do not push. No `R/` change and no template change.** `admitted_n` is already recorded (Part T2, `23648d7d`), so no Part T is needed. If anything appears to need a code change, STOP and report.

## What this completes

- Eight cells take `grfmr` to **18 of 18**, matching the DINA grid cell for cell.
- No acceptance criterion, no recommendation, anywhere.

## Settled already — do not re-determine

- **HR 1.00 is not a global null.** The planted region carries HR 1.00 against a benefiting complement (HR 0.657 at 12.4%, 0.721 at 31%), established from source in the `dinamr` Block C Gate 0 on this same template. Cite that determination; do not repeat it.
- **Naming.** Call it `selection_rate`. Never a false-selection rate: the planted region clears the deliberately sub-null log(0.90) floor, so returning it is an admissible selection. Use **no harm vocabulary** for these cells — "differentially null against a benefiting complement".
- **Nothing in the classification set is structurally undefined.** The planted region exists, so sensitivity, specificity, PPV, NPV and the oracle target (log 1.00 = 0) are all defined. The structural columns remain `n_cons_qual`, `band_n` and `p_star`, structural on every GRF cell and not on account of HR 1.00.
- **Captions** state the products, the cells, and that rates are over detected replicates. No conditional-versus-unconditional qualifier.
- **Stratify on `admitted_n`, never on `n_family`.** On GRF `n_family` is the outcome-independent enumerated pool.

## The question these cells answer

- At the null cells, FS's selection rate falls (0.62–0.68 at 12.4%, 0.92–0.96 at 31%) and DINA's falls hard (0.714 → 0.344 at 12.4%). GRF's null probe selected on 35 of 36.
- **If GRF selects at or near 1 even where the planted region is null, its selection rate carries no information and bound location is the entire diagnostic.** Report the selection rate and the location shares as two separate quantities and let the shares carry the question, exactly as the DINA Block C tables do.

## Cells — eight at 2,000 replicates

Knobs as `campaign.sh` pins them with the GRF engine: `FS_S7_METHOD=grf`, `FS_S7_FOCUS=effMaxSG`, `FS_S7_NBHD=0.20`, the four field knobs, `FS_S7_IJ_RESIDUAL=two_term`, `FS_S7_FB=none`, `FS_S7_CAMPAIGN=grfmr`, `FS_S7_ER_JCUTS` unset, `dmin.grf = 0.0`, seeds 8316951 + sim_id, two batches of 1,000 combined, 12 workers with the three thread variables at 1.

Run order, highest priority first; defer from the tail:

1. 31% (`FS_S7_Z1Q=0.60`), HR 1.75, n 1500 — deferred at the `grfmr` checkpoint.
2. 31%, HR 1.75, n 1000 — same.
3. 12.4% (`FS_S7_Z1Q` unset), HR 1.00 at n 500, 1000, 1500.
4. 31%, HR 1.00 at n 500, 1000, 1500.

This follows the established defer order: the harm n-trajectory ranks above the null block, and within the null block 31% defers before 12.4%.

## Gates

### Stage 1 — one 5-replicate smoke at (12.4%, HR 1.00, n 500), stop-on-failure

- `admitted_n` present and populated on the null path; p̂, ρᶜ and the nine recovery columns present and populated; all constructions finite; invariants and γ in range; realized prevalence checked.
- Record the selection rate and the `admitted_n` distribution.
- The `grfmr` probes produced one non-detection at the null corner that the mechanism diagnostic could not classify — `NO-DETECTION` rather than `CONFIG-ERROR`, so nothing errored. Record how many non-detections this smoke produces and whether `admitted_n` is NA or 0 on them. **A non-detection is not a failure; do not chase it, and do not add a recorder change.**

### Gate 1 — compute go/no-go

- **Ceiling 9 h wall; hard timeout 10 h.** Record both.
- Project per cell from the probe cost distribution, calibrated on `grfmr`'s realized-over-projected ratio of 0.896, charging overhead per batch render as `projectC.R` does.
- Reference: ≈ 1.9 h for the two harm cells and ≈ 4.9 h for the six null cells, ≈ 6.8 h in total. The null corner probe measured 13.50 s per replicate, the cheapest of the five.
- If the eight do not fit, defer from the tail and record it. **Do not reduce replicates** — 2,000 per cell is what makes the grid comparable to the FS and DINA grids.

### Gate 3 — alignment, per batch, stop-on-failure

`grf_select_statistic` → `"effect"`, `grf_selection` → `"frontier"`, `dmin.grf` → 0.0. Report all three resolved values.

### Gate 2 — per cell, via `gate2G.R`

- The corrected identity `log(est2_s) + lam_mean_s == log(est2) + lam_mean`, not the `fld_joint_bonf_*` vs `fld_joint_*` pair.
- Amendment 3's same-draws assertion against the committed FS comparator: `n_true` `identical()` on all 2,000 rows, truth `all.equal()` at 1e-8 with `identical()` and the maximum absolute and relative discrepancy beside it. A mismatch is a DGM-path finding, not a cell failure.
- Confirm the designated comparator resolves for each new cell before running the assertion; `e1stud` has no HR 1.00 bundle, and `map1`/`s7` bundles must not be substituted.
- Record per cell, prominently: selection rate; the `admitted_n` distribution (median, quartiles, p90, min, max); `n_family` beside it, stated to be the enumerated pool and not the qualified set.

## Stage 3 — extend, do not restart

- Extend `summary_grfmr.qmd` to the full 18 cells. The absent-cell guard should now skip nothing if all eight land.
- Carry through for the new cells: the `admitted_n` and p̂ strata with their joint count table; the error-SD Gaussian reference beside the marginal one with both formulas stated; Wilson intervals on every rate; marginal SD and error SD side by side; absolute coverage levels for every product.
- **Add the location table for the null cells**, the same columns the DINA null tables use: median field lower bound, median realized θ(Ĥ), the gap in difference, ratio and paired-ratio form, the planted target, and the shares of field lower bounds at or above 1.00 and at or above 1.25, with Wilson intervals.
- FS beside GRF with the comparator's `sg_focus` and ε per cell and the confound sentence with every comparison — criterion-matched at 31% (`effMaxSG` ε 0.20), not at 12.4% (`maxeffCons` ε 0.10).
- Quote FS product ranges only with the cell set and the construction named, per `REPORT_fs_products_reconciliation_2026-09-12.md`.

## Closeout

- Every artifact carrying payload or analysis content is committed, per the standing rule in `scripts_dinamr/README.md`. Report path, tracked status and size; flag anything over 50 MB, hard stop at 100 MB.
- **Last action: regenerate `current_status.md`**, with the post-condition that its stated pin equals HEAD at commit time (`check_current_status.sh`). Update §2.3 for the completed grid and §7 for what remains.

## Report

One-paragraph summary: cells completed, deferred and dropped with walls; the selection-rate and location-share headline at the null cells; the `admitted_n` distribution; Gate 1 projection against realized; and the commit range. Bullet form for the detail, one item per bullet.
