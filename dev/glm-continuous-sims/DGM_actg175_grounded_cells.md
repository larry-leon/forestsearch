# ACTG175-grounded DGM cells for the continuous/MD simulations

Sources: `actg175_literature_review.qmd` (cited by line) and
`analysis_actg175_continuous_compare_all.qmd` (the analysis actually run on
the dataset). Everything here is a configuration block for the existing
machinery — `calibrate_glm_interaction()` knobs inside the committed
twin/batch structure. No machinery changes; the committed −40 cell remains
the verified anchor and these run alongside it as arms.

## Cell L1 — literature-grounded subgroup

The published modifier set (litreview :140–143): age, Karnofsky, baseline
CD4 are the three most important treatment-effect modifiers, with the
median age split at 34; prior-ART status is the trial's ~43% naive stratum
(:49). Grounded Q, two candidate forms for sign-off:

```r
# L1a: published modifiers, harm on the OLDER side (current direction)
dgm_subgroup_vars <- c("age", "karnof")
dgm_subgroup_cuts <- list(age    = list(type = "greater", value = 34),
                          karnof = list(type = "less_eq", value = 90))

# L1b: ART-naive stratum in place of the quantile preanti cut
dgm_subgroup_vars <- c("age", "preanti")
dgm_subgroup_cuts <- list(age     = list(type = "greater", value = 34),
                          preanti = list(type = "less_eq", value = 0))  # naive
```

DECISION (yours): the age SIDE. The literature's interesting effects are
for *younger* patients (≤ 34); the current DGM plants harm on > 34. Either
is defensible — mirroring the literature's side (`less_eq 34`) makes the
citation direct; keeping > 34 preserves continuity with the committed cell.
State the choice in the cell's header either way.

## Cell L2 — analysis-mirroring thresholds

The earlier analysis's own configuration (compare_all :107–109), CD4
units, adverse scale:

```r
md_threshold   <- 10   # screening: MD >= 10 indicates harm
md_consistency <- 5    # consistency: MD >= 5
```

This resolves alignment decision D2 without touching the DGM's truth: the
simulation evaluates the procedure as it was configured on the real data —
externally justified floors, no truth-derived tuning. Note for
interpretation: under the −40/−26.26 truth both floors are sub-CATE
(permissive, manuscript-stress posture), so detection is expected near 1
and the tie dynamics of the committed grid apply; that is the point — this
cell mirrors practice, the committed 30/10 cell mirrors a discriminating
regime, and the pair brackets the design space.

## Cell L3 — dual-scale outcome framing (structural, all cells)

Adopt compare_all's dual-scale convention (:92–95) in the simulation
documents: construct `y_decline` (adverse scale, FS-native) and
`y_improve` (beneficial scale) once at data build, and feed each engine
its native column. This retires the orientation bridge as a per-readout
sign map — the fragility class behind CHECK C's catch and the FB
orientation fix — at the cost of one data-prep chunk. Ledger entry, not a
rebuild.

## Magnitude — to confirm, not assumed

Whether −40 CD4 cells/mm^3 is a literature-plausible harm magnitude should
be anchored to the observed effects in compare_all's results sections (and
the litreview's cd420 analyses) before L1 calibrates. If the observed ITT
and modifier-conditional effects suggest a different plausible scale, the
calibrator's `target_effect` is the one knob to move — with the SE-vs-floor
arithmetic (clearance probabilities) restated for the chosen value.

## What does not change

The identifier and inference configuration (maxeffCons, resample, MR 5000,
FB 300), the four-estimator schema, the batch/combine structure, the exact
β(Ĥ) scoring, and the committed −40 cell as the verified anchor. Each L
cell is one config block + one calibration run + the standard batch.
