# TASK — GRF MR campaign (`grfmr`), with the family determination and the FS extraction

- **Date:** 2026-09-11
- **Campaign:** `grfmr`. Report-and-wait; proceed unattended.
- **Machine:** Mac Studio. Fresh CC session.
- **Predecessors, committed:** `dev/tasks/TASK_dinamr_campaign_2026-09-10.md` and `dev/tasks/TASK_dinamr_blockC_grfprobe_2026-09-11.md`. Read them for the DGM, cell definitions, Gate 2 content and the structurally-NA column list. This document governs where they differ.
- **Tooling:** reuse `quarto/simulations/gbsg_020/scripts_dinamr/` (`campaign.sh`, `render.sh`, `gate2.R`, `project.R`, `grfprobe.sh`, `grf_mechanism.R`), all execution-verified with the `SP` export fixed. Anything new goes in the same directory and is committed.
- **Commit only, do not push.** Part T2 below is the **only** authorized code change, it is to the simulation template and not to `R/`, and it is conditional. If anything else appears to need a change, STOP and report.

## Framing for the record

- This campaign does not certify GRF. Every coverage number is coverage of β(Ĥ) conditional on the proposed family, over selected replicates, and every table must say so. Gate 0 may find the underlying picture more favourable than that; the caption does not change on account of it — see §Gate 0.
- GRF is **not** "FS-analogous". A GRF-to-FS or GRF-to-DINA comparison differs in identifier, family construction, detection set, and — at the DR pre-filter only, not at admission — the scale of the selection criterion.
- No acceptance criterion, no recommendation, anywhere.

---

## GATE 0 — three source determinations, no compute

Quote every finding with file and line numbers. Record all three in the report whatever they say.

### 0a. Which candidate set MR resampling re-evaluates — the determination that labels the campaign

- `.grf_dr_candidates()` (`R/grf_subgroup_labels.R:255–277`) enumerates from quantiles of X subject to `n_min`, so the enumerated pool does not depend on the outcome; the DR scores enter only the effect column, and the frontier select then admits on them.
- **Determine from source which of those two sets the MR resampling re-evaluates on each draw:** the full enumerated pool, or the forest-admitted subset. Trace the object that reaches the MR layer and say which it is.
- Record the mechanism and the deciding lines. **Do not draw any conclusion about whether GRF satisfies a fixed-family condition, and do not change any caption on account of this finding** — every table stays labelled conditional-on-proposed-family. The determination is for Larry.
- STOP only if the answer cannot be read from source. Ambiguity is a STOP; an inconvenient answer is not.

### 0b. Is `admitted_n` a recorded per-replicate column?

- GRF's `n_family` is the enumerated pool: near-constant, and identical across prevalences replicate by replicate. Stratifying on it would stratify on something outcome-independent, so `admitted_n` — the forest-qualified count — is the analogue of DINA's family-size stratifier and the campaign's core diagnostic.
- Determine whether `admitted_n` reaches the bundle as a column. If it does, note it and skip Part T2.
- If it does not, run **Part T2** below before Gate 1.

### 0c. What verification standard applies

- GRF fits are reproducible within a session (3 of 3) but not across contexts (0 of 179 rows matched between the probe diagnostic and the pipeline).
- State what that does and does not affect. Specifically: Amendment 3's same-draws assertion uses `n_true` and truth, which come from the DGM and not from GRF, so it is unaffected — confirm that from source rather than assuming it.
- State plainly that any GRF diagnostic re-run outside the pipeline measures statistically equivalent but not identical fits.

## PART T2 — conditional: record `admitted_n`

Run **only** if Gate 0b finds `admitted_n` absent from the bundle.

- Add `admitted_n` to the recorder in the simulation template. **Add-only**, written `NA_integer_` on every non-GRF path so that DINA and consistency bundles are byte-identical to what they are today.
- **Gate T2, stop-on-failure:** 5 replicates at the standing identity cell with `FS_S7_METHOD` unset — every non-timing column and truth `identical()` to the committed rows, on this machine, pre-change against post-change. Revert and STOP on failure.
- No `R/` change. If `admitted_n` cannot be recorded without one, STOP and report.

## GATE 1 — compute go/no-go

- **Ceiling 9 h wall; hard timeout 10 h.** Larry's available window is 9–10 h. Record both.
- Project per cell from the probes' per-replicate cost distribution, not a mean, calibrated on the `dinamr` realized-over-projected record (Block A 1.27 before the overhead correction, 0.949 after; Block B 0.92). Charge overhead as `projectC.R` does: per batch render, not 3 × uniform.
- Reference: the probes measured 13.5–19.8 s median per replicate, rising with n and with prevalence, and **flat in family size** (|ρ| ≤ 0.171) — the opposite of DINA's profile. Chat's rough projection is ≈ 9.3 h compute plus ≈ 1 h overhead for the 12 cells, so expect to defer.
- If the twelve do not fit, defer from the tail of the run order. **Do not reduce replicates** — 2,000 per cell is what makes the grid comparable to the FS and DINA grids and is not traded for cell count.

## PART A — twelve harm cells at 2,000 replicates

Knobs as `campaign.sh` pins them, with the GRF engine: `FS_S7_METHOD=grf`, `FS_S7_FOCUS=effMaxSG`, `FS_S7_NBHD=0.20`, the four field knobs, `FS_S7_IJ_RESIDUAL=two_term`, `FS_S7_FB=none`, `FS_S7_CAMPAIGN=grfmr`, `FS_S7_ER_JCUTS` unset, `dmin.grf = 0.0`, seeds 8316951 + sim_id, two batches of 1,000 combined, 12 workers with the three thread variables at 1.

Run order, highest priority first; defer from the tail:

1. 12.4% (`FS_S7_Z1Q` unset): HR 1.50 at n 500, 1000, 1500.
2. 31% (`FS_S7_Z1Q=0.60`): HR 1.50 at n 500, 1000, 1500.
3. 12.4%: HR 1.75 at n 500, 1000, 1500.
4. 31%: HR 1.75 at n 500, 1000, 1500.

This keeps the HR 1.50 n-trajectory complete at both prevalences first, mirroring the `dinamr` checkpoint decision. The HR 1.00 null cells are **not** in this task; they go to a follow-up session.

- **Gate 3, alignment, per batch, stop-on-failure:** `grf_select_statistic` must resolve to `"effect"`, `grf_selection` to `"frontier"`, `dmin.grf` to 0.0. Report all three resolved values.
- **Gate 2 per cell**, via `gate2.R`, exactly as for `dinamr`: the corrected identity `log(est2_s) + lam_mean_s == log(est2) + lam_mean`; Amendment 3's same-draws assertion against the committed FS comparator (`n_true` `identical()` on all 2,000 rows, truth `all.equal()` at 1e-8 with `identical()` and the maximum absolute and relative discrepancy beside it; a mismatch is a DGM-path finding, not a cell failure); structurally-NA columns reported as such.
- Record per cell, prominently: detection; **`admitted_n` distribution** (median, quartiles, 90th percentile, min, max); `n_family` beside it, stated to be the enumerated pool and not the qualified set.

## STAGE 3 — a new summary, transplanted

- Transplant `summary_dinamr.qmd` to `summary_grfmr.qmd`: globs and labels only, plus the GRF-specific substitutions below. Guard per-cell chunks so absent cells skip rather than render empty.
- Every coverage column labelled as the conditional-on-proposed-family estimand.
- **Substitute `admitted_n` for `n_family` as the stratifier** in the strata section, and state in the caption why: `n_family` on GRF is the outcome-independent enumerated pool. Keep the `n_family` column in the descriptive tables.
- Carry through unchanged: the p̂ strata and the joint count table; the error-SD Gaussian reference beside the marginal one with both formulas stated; Wilson intervals on every rate; marginal SD and error SD side by side; absolute coverage levels for every product, not only differences from FS.
- FS beside GRF with the FS comparator's `sg_focus` and ε stated per cell, and the confound sentence with every comparison. Criterion is matched at 31% (`effMaxSG` ε 0.20) and not at 12.4% (`maxeffCons` ε 0.10).
- GRF's detection was 1.0000 at all four harm probes. If that holds on the full cells, say so and note that it removes the detection-conditioning confound that qualifies every DINA n-trend.

## PART C — FS extraction for the manuscript (no compute)

Runs **only** while cells are computing, and **never** ahead of a gate, a projection or a launch. If it would delay any of those, defer it.

From the committed FS comparator bundles (tier2, p12ext at 12.4%; cert20, e1stud at 31%), for every cell those bundles cover — harm and null:

- Classification against the planted region: sensitivity, specificity, PPV, NPV, and mean |Ĥ| against |H|, with Wilson intervals on the rates.
- Bound location on the HR scale, the same columns as the DINA location tables: median field lower bound, median realized θ(Ĥ), the gap in difference, ratio and paired-ratio form, and the shares of field lower bounds at or above 1.00 and at or above 1.25.
- The FS comparator's `sg_focus` and ε beside every row.

Write these to `REPORT_fs_extraction_2026-09-11.md` beside the campaign report. **Reading committed bundles only — no re-run, and no new simulation.** If a quantity is not derivable from the committed columns, say so and do not add a recorder change for it.

## Report

End with a one-paragraph summary: the three Gate 0 determinations; whether Part T2 ran; Part A cells completed, deferred and dropped with walls; the detection and `admitted_n` headline; Gate 1 projection against realized; whether Part C completed; and the commit range. Bullet form for the detail, one item per bullet.
