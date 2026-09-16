# COLUMNS — md_field_metrics.csv (campaign mdsgnb20; TASK_md_field_rerun_2026-09-15 §3.2)

Written by `summary_continuous_field_mdsgnb20.qmd` from the same objects its tables print. One row per cell × block × estimator × metric (× tau for bound-location rows).

## Columns

- `campaign`: `mdsgnb20`, or `mdf1` for the paired comparator rows of the rule contrast (declaration rate, identification, field one-sided coverage on the same seeds).
- `cell`: `md40 n500`, `md120 n500`, `null n500 (no subgroup; homogeneous +26)`, `md40 n700`.
- `block`: `H` (the selected subgroup Ĥ), `Hc` (its complement Ĥᶜ), `joint` (the pair).
- `estimator`: see the code map below.
- `metric`: see the metric list below.
- `tau`: the D3 ladder threshold for bound-location rows (0, 10, 20, 30, 40, 60, 80, 100 on the harm-oriented MD scale); empty otherwise.
- `value`: the metric.
- `mc_se`: Monte Carlo standard error -- proportions sqrt(p(1-p)/n); bias SD/sqrt(n); empirical SD SD/sqrt(2(n-1)); mean SE SD(SE)/sqrt(n); SE/SD by the delta method (r * sqrt((mc_se(mean SE)/mean SE)^2 + (mc_se(SD)/SD)^2)); empty where none is defined (quantiles, display bookkeeping).
- `wilson_lo`, `wilson_hi`: Wilson 95% interval for proportions; empty otherwise.
- `n`: the denominator (declared replicates with the quantity defined; total replicates for the declaration rate).
- `commit`: the commit that added the cell's combined bundle.

## Estimator codes

- `naive`: the unadjusted plug-in estimate on the selected region (one-sided bound est -/+ 1.645 SE, robust SE).
- `oracle`: the refit on the true region Q / Q^c (its target is the structural effect; empty on Ĥ in the null cell).
- `mr`: MR (IJ two-term): the de-biased estimate with the infinitesimal-jackknife SE (one-sided bound est -/+ 1.645 SE_IJ).
- `fld`: MR (field): the shrunk-field estimate `est2` with the Λ*-quantile intervals; on Ĥ the one-sided lower bound `lo1s`; on Ĥᶜ the unstudentized complement field (`up1s`), shown as the paired before/after.
- `fld_s`: MR (field-s): the studentized complement field on Ĥᶜ (`est2_s`, `up1s_s`, `se_field_s`) -- the evaluated complement construction.
- `bonf_s` / `bonf` / `separate_s`: the field-s Bonferroni pair (gamma = 0.025 each side) / the unstudentized Bonferroni pair / the separate 95% pair (field lower on Ĥ, field-s upper on Ĥᶜ).
- `calibrated`: diagnostics of the calibrated joint split (mean gamma, corr(Λ*, Λ*ᶜ)), `_s` = studentized.
- `all`: cell-level quantities not tied to an estimator (declaration rate, identification, p-hat).

## Metrics

- `declaration_rate`: declared replicates / replicates.
- `bias_md`, `bias_sd_units`: mean(estimate - target), in MD units and divided by the empirical SD; the target is beta(Ĥ) / beta(Ĥᶜ) per replicate (the structural effect for the oracle).
- `sd_emp`, `se_mean`, `se_over_sd`: empirical SD of the estimate, mean reported SE, their ratio.
- `cov1_lower` (H) / `cov1_upper` (Hc): one-sided 95% coverage on the exposed side; `cov2`: two-sided 95% coverage.
- `halfwidth_md`, `margin1_md`: mean two-sided half-width; mean one-sided margin from beta-tilde (field rows) or from the estimate.
- `share_lower_ge_tau` (H): share of replicates whose one-sided lower bound is >= tau ('harm of at least tau supported'); `share_upper_le_tau` (Hc): share whose one-sided upper bound is <= tau ('harm of at most tau supported').
- `bound_mean`, `bound_q05` ... `bound_q95`: location of the one-sided bound.
- `joint_coverage`, `share_both_bounds`, `cov_H`, `cov_Hc`, `margin_H_md`, `margin_Hc_md`: joint pair rows on declared replicates carrying both bounds; margins in MD units from beta-tilde.
- `gamma_mean`, `corr`, `gamma_mean_s`, `corr_s`: calibrated-split diagnostics.
- `mean_n_harm`, `mean_n_sel`, `sensitivity_mean`, `ppv_mean`, `share_nharm_{grew,stayed,shrank}_vs_mdf1`, `count_nharm_{grew,stayed,shrank}_vs_mdf1`: identification, paired with mdf1 by sim_id (shares and counts of the paired replicates).
- `p_hat_mean`, `p_hat_share_lt_05`: re-selection frequency of the winner and the tie-regime share.
- `sd_btc_over_naive_se` (mr), `lamsd_over_naive_se` (fld, fld_s): regime diagnostics on Ĥᶜ.
- `display_b`, `display_r`, `display_cov1`, `display_cov1_ref`, `display_cov2`, `display_cov2_ref`: fs_sim_bias_coverage(scale = 'identity') bookkeeping.

## Scale convention

Every estimate, bound, bias and threshold is on the **harm-oriented mean-difference scale** (positive = harm; the gate works on -cd4_change because adverse_outcome = FALSE); `betaHhat_*` in the bundles are raw cd4_change and are oriented with -1 here. The null cell's truth is a homogeneous +26.255 on this scale for every Ĥ and Ĥᶜ. No significance language: bounds are read by location against the ladder.
