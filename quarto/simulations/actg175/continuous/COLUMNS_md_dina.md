# COLUMNS — md_dina_metrics.csv (campaign mddina; TASK_md_dina_campaign_2026-09-17 §3.2)

Written by `summary_continuous_field_mddina.qmd` from the same objects its tables print. The schema of `md_grf_metrics.csv` (see `COLUMNS_md_grf.md`, itself `md_field_metrics.csv`'s plus `identifier`). One row per cell × block × estimator × metric (× tau for bound-location rows).

## Columns

- `campaign`: `mddina`; `mdsgnb20` for the FS comparator rows; `mdgrf` for the GRF comparator rows.
- `identifier`: `dina` (this campaign); `fs` for the FS comparator rows, copied from the committed `md_field_metrics.csv`; `grf` for the GRF comparator rows, copied from the committed `md_grf_metrics.csv`. Comparator rows keep their own `commit` and are not recomputed.
- `cell`: `md40 n500`, `md120 n500`, `null n500 (no subgroup; homogeneous +26)`, `md40 n700`.
- `block`: `H` (the selected subgroup Ĥ), `Hc` (its complement Ĥᶜ), `joint` (the pair), `all` (timing).
- `estimator`, `metric`, `tau`, `value`, `mc_se`, `wilson_lo`, `wilson_hi`, `n`: as in `COLUMNS_md_field.md`.
- `commit`: the commit that added the cell's combined bundle (for `fs` and `grf` rows, the value in the committed extract).

## Estimator codes

`naive`, `oracle`, `mr`, `fld`, `fld_s`, `bonf_s`, `bonf`, `separate_s`, `calibrated`, `all`: as in `COLUMNS_md_field.md`, applied to the DINA-selected Ĥ.

## Metrics

As in `COLUMNS_md_grf.md` (itself `COLUMNS_md_field.md`), with these identification rows:

- `mean_size_hhat`: mean |Ĥ| over declared replicates (the recorder's `n_sel`; its `n_harm` is the same count, the size of Ĥ, not a true-positive count).
- `mean_true_positives`: mean number of truly harmed patients in Ĥ, sensitivity × `n_true`, over declared replicates with `n_true` > 0 (undefined in the null cell).
- `sensitivity_mean`, `ppv_mean`: the recorder's `sens`, `ppv` (`.classify`), mean over declared replicates.
- `mean_proposed_n`: mean of DINA's proposed family (`dina_proposed_n`: candidates with oriented tau-hat at or above the floor 30 and at least `n.min` members); `mean_n_family`: mean of MR's kept family size K; `mean_admitted_n`: mean of DINA's admitted set (proposed candidates whose harm-oriented MD clears the admission floor 30). All three over declared replicates: DINA's selection object, which carries the proposed and admitted counts, exists only when a subgroup is selected.
- `share_size_larger_than_fs`, `share_size_equal_than_fs`, `share_size_smaller_than_fs` and the matching `count_size_*_than_fs`: DINA's |Ĥ| against FS's (`mdsgnb20`) on the same `sim_id`, over replicates both identifiers declared.
- FS rows (`identifier = fs`): `declaration_rate`, `mean_n_sel`, `sensitivity_mean`, `ppv_mean`, and the three coverage rows of the three-identifier table (`H`/`fld`/`cov1_lower`, `Hc`/`fld_s`/`cov1_upper`, `joint`/`bonf_s`/`joint_coverage`). **FS's `mean_n_harm` equals its `mean_n_sel`** in `md_field_metrics.csv` (a labelling defect: `n_harm` is |Ĥ|, not a true-positive count); it is not copied.
- GRF rows (`identifier = grf`): `declaration_rate`, `mean_size_hhat`, `sensitivity_mean`, `ppv_mean`, `mean_n_family`, `mean_admitted_n`, and the same three coverage rows.
- `secs_fit_mr_mean`, `secs_field_mean`, `secs_complement_mean` (block `all`): mean seconds per replicate for the fit with MR (all replicates), the field pass and the complement field (declared replicates), as the regime table prints them; `fit_mr_secs` contains the other two, which are never summed.

## Scale convention

Every estimate, bound, bias and threshold is on the **harm-oriented mean-difference scale** (positive = harm; adverse_outcome = FALSE); `betaHhat_*` in the bundles are raw cd4_change and are oriented with -1 here. The null cell's truth is a homogeneous +26.255 on this scale for every Ĥ and Ĥᶜ. **Every DINA and GRF coverage figure is coverage of the estimand conditional on the proposed family**: their candidate families are generated from fitted surfaces, so the fixed-family condition does not hold. FS's family is the prespecified cut grid. Comparisons across the three identifiers are descriptive: the identifier, the family construction and the detected set all differ. No significance language: bounds are read by location against the ladder.
