# COLUMNS — md_grf_metrics.csv (campaign mdgrf; TASK_md_grf_2026-09-16 §3.2)

Written by `summary_continuous_field_mdgrf.qmd` from the same objects its tables print. The schema of `md_field_metrics.csv` (see `COLUMNS_md_field.md`) plus an `identifier` column. One row per cell × block × estimator × metric (× tau for bound-location rows).

## Columns

- `campaign`: `mdgrf`, or `mdsgnb20` for the FS comparator rows.
- `identifier`: `grf` (this campaign), or `fs` for the FS comparator rows, copied from the committed `md_field_metrics.csv` with their own `commit`, not recomputed.
- `cell`: `md40 n500`, `md120 n500`, `null n500 (no subgroup; homogeneous +26)`, `md40 n700`.
- `block`: `H` (the selected subgroup Ĥ), `Hc` (its complement Ĥᶜ), `joint` (the pair).
- `estimator`, `metric`, `tau`, `value`, `mc_se`, `wilson_lo`, `wilson_hi`, `n`: as in `COLUMNS_md_field.md`.
- `commit`: the commit that added the cell's combined bundle (for `fs` rows, the value in the committed FS extract).

## Estimator codes

`naive`, `oracle`, `mr`, `fld`, `fld_s`, `bonf_s`, `bonf`, `separate_s`, `calibrated`, `all`: as in `COLUMNS_md_field.md`, applied to the GRF-selected Ĥ.

## Metrics

As in `COLUMNS_md_field.md`, except the identification rows:

- `mean_size_hhat`: mean |Ĥ| over declared replicates (the recorder's `n_sel`; its `n_harm` is the same count, the size of Ĥ, not a true-positive count).
- `mean_true_positives`: mean number of truly harmed patients in Ĥ, sensitivity × `n_true`, over declared replicates with `n_true` > 0 (undefined in the null cell).
- `sensitivity_mean`, `ppv_mean`: the recorder's `sens`, `ppv` (`.classify`), mean over declared replicates.
- `mean_n_family`: mean of MR's kept family size K; `mean_admitted_n`: mean of GRF's admitted set (candidates clearing the harm-oriented MD floor 30).
- `share_size_larger_than_fs`, `share_size_equal_than_fs`, `share_size_smaller_than_fs` and the matching `count_size_*_than_fs`: GRF's |Ĥ| against FS's (`mdsgnb20`) on the same `sim_id`, over replicates both identifiers declared.
- FS rows (`identifier = fs`): `declaration_rate`, `mean_n_sel`, `sensitivity_mean`, `ppv_mean`. The FS extract's `mean_n_harm` equals its `mean_n_sel` (a labelling defect: `n_harm` is |Ĥ|) and is not copied.
- The `mdf1` rule-contrast metrics of the FS extract (`*_vs_mdf1`) have no counterpart here.

## Scale convention

Every estimate, bound, bias and threshold is on the **harm-oriented mean-difference scale** (positive = harm; adverse_outcome = FALSE); `betaHhat_*` in the bundles are raw cd4_change and are oriented with -1 here. The null cell's truth is a homogeneous +26.255 on this scale for every Ĥ and Ĥᶜ. **Every coverage figure is coverage of the estimand conditional on the proposed family**: GRF's candidate family is generated from a fitted surface, so the fixed-family condition does not hold. Comparisons with FS are descriptive. No significance language: bounds are read by location against the ladder.
