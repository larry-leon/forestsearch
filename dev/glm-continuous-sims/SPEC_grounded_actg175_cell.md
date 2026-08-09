# Grounded ACTG175 continuous cell — extraction (G1) then design (G2)

Goal: the companion paper's continuous simulation section, under a DGM that
mimics the realistic ACTG175 analysis in
`analysis_actg175_continuous_compare_all.qmd`. Machinery unchanged: the
committed twin (`sim_fs_maxeffCons_fb_mr_*`), the calibrator, exact-β(Ĥ)
scoring. Only DGM knobs and the outcome construction change.

## G1 — extract the data anchors (no decisions; verbatim numbers)

From the compare_all document (re-render the needed chunks, or read its
written .rds bundles if present in the repo/local dir), report verbatim:

1. Analysis cohort N (arms 1 vs 3, non-missing cd420).
2. ITT on the adverse scale: the `lm(y_decline ~ treat)` treat coefficient
   — estimate, SE.
3. Residual SD of `y_decline` (overall, and by arm if one line).
4. The tidy selected-subgroups table for the `maxeffCons` arm:
   definition, n, MD estimate, Pcons, flags — every row, verbatim.
5. If the analysis produced any de-biased (MR/FB) estimate for a selected
   region, report it beside the naive.

Deliverable: `dev/glm-continuous-sims/G1_actg175_anchors.md`, committed and
pushed. Stop there — G2's choices are made from this table.

## G2 — the grounded cell (config only; runs after the picks)

Two decisions from the G1 table, then one config block:

- **Q** = the chosen real-data region (or its literature-cleaned form —
  the age-34 modifier and ART-naive stratum per the design note).
- **target_effect** = the chosen planted truth on the y_decline scale.
  The real selected-region MD is selection-inflated; anchor the truth
  BELOW it (the de-biased value if G1 item 5 exists, else a stated
  fraction of the naive), so the simulation measures exactly that
  inflation against a truth chosen with eyes open.

Fixed by the data and the analysis itself:

```r
# Outcome: the analysis's own adverse-scale construction (positive = harm)
actg_df$y_decline <- actg_df$cd40 - actg_df$cd420   # compare_all :208
outcome_name    <- "y_decline"
adverse_outcome <- TRUE          # y_decline is already adverse-oriented
n_sample        <- <G1 item 1>   # the real cohort size
md_threshold    <- 10            # the analysis's own screening floor
md_consistency  <- 5             # the analysis's own consistency floor
```

Everything else — maxeffCons, resample, MR 5000, FB 300 on the anchor
runs, the batch/combine structure, four-estimator schema — is the
committed twin verbatim. The −40 planted cell remains as the machinery-
validation and comparison arm; this cell is the paper's scenario.

G2 adopts the same J = 10 quantile grids on all six continuous covariates as
Step 1's real-data arm (`age`, `preanti`, `wtkg`, `karnof`, `cd40`, `cd80`),
with the FB-arm cost priced from a measured smoke before any batch schedule is
set — and with `max_subgroups_search` reviewed alongside them, since at J = 10
the committed cap of 30 discarded 98 of 128 candidate subgroups before
evaluation (G1 § 8, finding 3).
