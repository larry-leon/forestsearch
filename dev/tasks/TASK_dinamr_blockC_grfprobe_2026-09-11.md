# TASK — `dinamr` Block C and the deferred Block B cell, then GRF cost probes

- **Date:** 2026-09-11
- **Campaign:** `dinamr` (continuation). Report-and-wait; proceed unattended.
- **Machine:** Mac Studio. Fresh CC session.
- **Predecessor:** `dev/tasks/TASK_dinamr_campaign_2026-09-10.md` (committed). This document does not restate it — read it for the DGM, cell definitions, knob set, Gate 2 content and the structurally-NA column list. Where the two differ, this document governs.
- **Commit only, do not push.** No `R/` change. If anything appears to need one, STOP and report.

## Framing for the record

- This does not certify DINA. The fixed-family condition does not hold, so every coverage number is coverage of the conditional-on-proposed-family estimand and every table must say so.
- Part B is cost measurement for GRF, not a GRF campaign and not an adopted GRF configuration. No coverage claim, no acceptance criterion, no recommendation comes out of Part B.

## Tooling — reuse what is committed

- The campaign's drivers and checkers are committed and execution-verified in `quarto/simulations/gbsg_020/scripts_dinamr/`. Reuse them rather than rewriting: `campaign.sh` (knob set, 12-worker pinning), `render.sh` (the three thread variables), `gate2.R` (Gate 2, corrected identity, `TOL_TRUTH <- 1e-8`), `project.R` (Gate 1 projection), `checkpoint.R`, `stage1_checks.R`, `probe.sh`.
- Add a `blockC.cells` list in that directory beside `blockA.cells` / `blockB.cells`, and commit it.
- `stage1_checks.R` still carries the superseded bonf-vs-raw comparison, labelled as such. Use `gate2.R`'s corrected identity for any gate.
- Anything new this task needs goes in the same directory and is committed, not left in the session scratchpad.

## Part A — seven cells at 2,000 replicates

Knobs, seeds and criterion exactly as the predecessor's Part C, as pinned in `campaign.sh`: `FS_S7_METHOD=dina`, `FS_S7_FOCUS=effMaxSG`, `FS_S7_NBHD=0.20`, the four field knobs, `FS_S7_IJ_RESIDUAL=two_term`, `FS_S7_FB=none`, `FS_S7_CAMPAIGN=dinamr`, `FS_S7_ER_JCUTS` unset, seeds 8316951 + sim_id, two batches of 1,000 combined, 12 workers with the three thread variables at 1. Replicate count is 2,000 per cell and is not traded for cell count.

Run order, highest priority first; defer from the tail:

1. Block B, HR 1.75, n 1500, 31% (`FS_S7_Z1Q=0.60`) — the cell deferred at the Block A checkpoint.
2. Block C at 12.4% (`FS_S7_Z1Q` unset): HR 1.00 at n 500, 1000, 1500.
3. Block C at 31% (`FS_S7_Z1Q=0.60`): HR 1.00 at n 500, 1000, 1500.

This order follows the predecessor's stated defer order (Block C defers before B's n1500; within Block C, 31% defers before 12.4%).

### Gate 0 — Stage 0, verify from source

- Quote with file and line numbers: Block C's cell definition in the predecessor task document, and whatever in the harness makes a planted-harm quantity undefined at HR 1.00.
- State which quantities are structurally undefined at a null (the oracle target, and any classification metric defined against a planted H) and how each will be reported. Report them as structurally undefined, never as failures — as with `n_cons_qual`, `band_n`, `p_star`.
- STOP if Block C's definition in the predecessor document is ambiguous about whether HR 1.00 is a global null or a planted region with HR 1.00 inside it. Do not choose.

### Gate 1 — compute go/no-go

- **Ceiling 9 h wall for Part A; hard timeout 12 h.** Record the ceiling and its source in the Gate 1 section.
- Null-cell family size is unmeasured. State which of the original Gate 1 probes covered a Block C cell; if none did, run 36-replicate probes at (12.4%, HR 1.00, n 500) and (31%, HR 1.00, n 1500) before projecting.
- Project from the family-size distribution, never from a mean and never from the FS walls. Prior realized-over-projected ratios: Block A 1.27, Block B 0.92.
- Reference projections to beat or correct: the deferred B cell 2.225 h (the checkpoint's by-n re-projection), Block C 2.17 h (the original Gate 1).
- If the seven cells do not fit, defer from the tail of the run order and record the decision. Do not reduce replicates.

### Gate 2 — per cell, as the predecessor defines it

- Use the corrected bound↔quantile identity — field-s inverted around the same bdc, `log(est2_s) + lam_mean_s == log(est2) + lam_mean` — not the `fld_joint_bonf_*` vs `fld_joint_*` pair, which coincide only when γ sits at its 0.025 floor.
- Same-draws assertion where a committed FS comparator shares the DGM draws: `n_true` `identical()` on all 2,000 rows, truth by `all.equal()` at tolerance 1e-8 with `identical()` and the maximum absolute and relative discrepancy reported beside it. A mismatch is a finding about the DGM path, not a cell failure.
- Record detection and the family-size distribution prominently. At a null, detection is a false-selection rate — label it as such, not as a detection rate.

### Stage 3 — extend the existing summary

- Extend `summary_dinamr.qmd`; do not start a new document. The absent-cell guard should now skip nothing at 12.4%/31% harm and null if all seven land.
- Every coverage column stays labelled as the conditional-on-proposed-family estimand.
- Carry through, for the new cells: the `n_family` and p̂ strata with their joint count table, the error-SD Gaussian reference beside the marginal one with both formulas stated, and the FS comparator's `sg_focus` and ε beside every FS number with the confound sentence extended to the criterion wherever it differs.
- For the 31% cells state that the FS comparator is criterion-matched (`effMaxSG` ε 0.20); for 12.4%, that it is not (`maxeffCons` ε 0.10).
- Report absolute coverage levels for every product, not only differences from FS.

## Part B — GRF cost probes

Runs only after Part A completes or defers. **Hard cap 1.5 h wall.** If Part A leaves less than that under the 12 h timeout, skip Part B and say so.

- `FS_S7_METHOD=grf`, otherwise the Part A knob set. `dmin.grf = 0.0` — decided by Larry on 2026-09-11; record the decision and its rationale (below) in the probe record.
- **Gate 3, alignment — STOP on failure:** `grf_select_statistic` must resolve to `"effect"` and `grf_selection` to `"frontier"`. Report both resolved values.
- Five 36-replicate probes: (12.4%, HR 1.50, n 500), (12.4%, HR 1.50, n 1500), (31%, HR 1.50, n 500), (31%, HR 1.50, n 1500), (12.4%, HR 1.00, n 500).
- Record per probe:
  - wall per replicate — median, 90th percentile, maximum — and wall against family size;
  - proposed-family size quantiles;
  - detection, with a Wilson interval;
  - **how often the frontier band comes back empty, and what a replicate records when it does** — GRF's `effMaxSG` applies the band as a frontier-only filter with no empty-band fallback, and the frequency of that path is the evidence for Larry's remaining open decision;
  - peak memory;
  - whether p̂, ρᶜ and the nine recovery columns are present and populated on the GRF path.
- If the empty-band path errors or produces a malformed record, STOP and report. Do not fix it.
- **No coverage table, no FS comparison, no acceptance criterion, no recommendation.** Cost and mechanism only.

### Rationale to record for `dmin.grf = 0.0`

- Larry, 2026-09-11: GRF's DR-scores target RMST for survival outcomes, whereas FS and DINA both target the Cox hazard ratio. FS's and DINA's floors are alignable with each other; GRF's is not alignable with either, so there is no GRF value that reproduces DINA's sub-null log(0.90) floor. 0.0 is the null point on GRF's own scale and reproduces the setting used in the manuscript's own GRF runs.
- Consequence for every GRF record: GRF must not be called "FS-analogous". A GRF-to-FS or GRF-to-DINA comparison differs in identifier, family construction, detection set, selection criterion **and the scale of the selection criterion**.
- Unchanged: the inference products are computed on β(Ĥ) via the Cox model on the identified region whichever identifier proposed it, so the estimand is the same kind of object across all three.

## Report

End with a one-paragraph summary: Part A cells completed, deferred and dropped with walls; the null-cell detection and family-size headline; the Gate 1 projection against realized; and for Part B the cost surface, the empty-band frequency and the commit range. Bullet form for the detail, one item per bullet.
