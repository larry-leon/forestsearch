# COLUMNS — or_metrics.csv (campaigns orfs, orgrf, ordina; TASK_actg175_binary_campaign_2026-09-17 §3.2; Stage 3 run under TASK_actg175_or_stage3_2026-09-28)

Written by `summary_actg175_or.qmd` from the same objects its tables print. The schema of `../continuous/md_dina_metrics.csv` (see `COLUMNS_md_dina.md`) plus a `design` column. One row per design × cell × identifier × block × estimator × metric (× tau for bound-location rows).

## Columns

- `campaign`: `orfs`, `orgrf` or `ordina`; `maxeffCons_actg175_or075_seedtab_s1000` for the rows quoted from the committed study grid.
- `identifier`: `fs`, `grf` or `dina` — the identifier that produced the row. Every `grf` and `dina` coverage row is coverage of the estimand **conditional on the proposed family**.
- `design`: `or075` (the supplement's protective design), `or150` (harm) or `or100` (borderline null: the planted region at the null against a protective complement). All three share the design of record's planted prevalence 14.917% (`sg_quantile` 0.62850, `REPORT_binary_redesign_2026-09-18.md`) and the complement references θ†(Ĥᶜ) 0.6564 / θ‡(Ĥᶜ) 0.6314; they differ only in `target_effect`.
- `cell`: `OR 0.75 n500`, `OR 0.75 n2000`, `OR 1.5 n500`, `OR 1.5 n2000`, `OR 1.0 n500`, `OR 1.0 n2000`.
- `block`: `H` (the selected subgroup Ĥ), `Hc` (its complement Ĥᶜ), `joint` (the pair), `all` (timing).
- `estimator`: `naive` (unadjusted), `oracle`, `mr` (IJ two-term), `fld` (field), `fld_s` (field-s), `bonf_s`, `bonf`, `separate_s`, `calibrated`, `all`.
- `metric`: see below. `tau`: the ladder point, on bound-location rows only.
- `value`, `mc_se`: the figure and its Monte Carlo SE. For a proportion, `mc_se` is sqrt(p(1-p)/n) and `wilson_lo` / `wilson_hi` are the Wilson 95% interval; for a mean, `mc_se` is SD/sqrt(n) and the Wilson columns are empty.
- `n`: the denominator the figure was computed on.
- `commit`: the commit that added the cell's combined bundle; `committed grid` for the quoted study rows.

## Metrics

- Estimation, per block × estimator: `bias_or` (mean estimate minus target, OR difference), `bias_log` (mean log estimate minus log target), `bias_sd_units` (`bias_log` / SD(log est)), `sd_emp_or`, `sd_emp_log`, `se_mean` (the stored SE, log-OR scale), `se_over_sd` (`se_mean` / SD(log est)), `halfwidth_or`, `margin1_or`.
- Coverage, per block × estimator: `cov2` (two-sided, the estimator's own interval) and `cov1_lower` / `cov1_upper` (one-sided on the block's exposed side), each against the row's target — β(Ĥ) for the conditional rows and θ† for the `oracle` row. The harm block additionally carries `cov2_theta_dagger`, `cov1_lower_theta_dagger`, `cov2_theta_ddagger` and `cov1_lower_theta_ddagger`.
- Bound location: `share_lower_ge_tau` (Ĥ) and `share_upper_le_tau` (Ĥᶜ) at each `tau`; `bound_mean` and `bound_q05` … `bound_q95`.
- Joint: `joint_coverage`, `share_both_bounds`, `cov_H`, `cov_Hc`, `margin_H_or`, `margin_Hc_or`; `gamma_mean`, `gamma_mean_s`, `corr`, `corr_s` under estimator `calibrated`.
- `n_nonconvergent_fits` (every block × estimator): declared replicates that are NON-CONVERGENT FOR THAT ESTIMATOR — its point estimate or either two-sided bound is non-finite, or its point estimate is <= 0. Such rows are excluded from that estimator's coverage, location, spread and bound-location statistics, exactly as the committed study's `is.finite(lo) & is.finite(hi)` masks exclude them from coverage, and the count is carried here and as a `non-convergent` column in every table that has an affected row. No row is dropped silently, and no estimator's rows are dropped on another estimator's account. The motivating case is complete separation in an arm of the true region (the study's `.logit_or_ci()` guards only the OVERALL >= 5 events / >= 5 non-events), where the logistic MLE diverges and the Wald interval degenerates to (0, Inf) with a point estimate of order 1e8. **The data recipe is the committed study's, verbatim: every recorded value is exactly what the recorder wrote. This is a consumer-side convention.**
- Identification (block `H`, estimator `all`): `declaration_rate`, `mean_size_hhat` (the recorder's `n_sel`), `mean_true_positives` (sensitivity × `n_true`, never |Ĥ|), `sensitivity_mean`, `ppv_mean`, `spec_mean`, `mean_n_true`, `mean_n_family` (MR's kept family K), `mean_admitted_n` (GRF's forest-qualified / DINA's admitted count; NA on FS), `mean_proposed_n` (DINA only), `mean_n_cons_qual` and `mean_band_n` (FS only).
- Regime: `p_hat_mean`, `p_hat_share_lt_05`; `sd_btc_over_naive_se` (block `Hc`, estimator `mr`) and `lamsd_over_naive_se` (block `Hc`, estimators `fld` and `fld_s`), all on the log-OR scale.
- Display: `display_b`, `display_r`, `display_cov1`, `display_cov1_ref`, `display_cov2`, `display_cov2_ref` from `fs_sim_bias_coverage(scale = "log")`.
- Timing (block `all`, estimator `all`): `secs_fit_mr_mean` (the fit with MR, all replicates), `secs_id_mean` (identification inside forestsearch()), `secs_field_mean` and `secs_complement_mean` (declared replicates). `fit_mr_secs` CONTAINS the other three, which are never summed.
- Quoted study rows: `study_declaration_rate` and `study_cov2_theta_dagger` — the committed study's detection rate and two-sided coverage of θ† under ITS rule (`maxeffCons`, ε 0.10), IJ intervals only, 1,000 replicates per cell. Read from its grid, never recomputed.

## Scale convention

Every estimate, bound and threshold is an **odds ratio**; `adverse_outcome = TRUE`, so **OR > 1 is harm** and no orientation flip is applied anywhere. Two scales coexist and are never mixed (`R/fs_mr_inference.R:480–488`): the **effect** scale (the OR) carries every estimate and bound, and the **working** scale (the log-OR) carries the stored Wald and IJ SEs, the field's `lambda_mean` / `se_field` and the Λ* quantiles. So every bias in SD units, every SE/SD and every regime ratio here is computed on the log-OR scale, while `bias_or`, `halfwidth_or`, `margin*_or` and the bound columns are ORs.

## Reading convention

Bounds are read by **location** against the ladder τ = 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 2.0, never as significance at OR = 1: on Ĥ, a one-sided 95% lower bound at or above τ reads "harm of at least τ"; on Ĥᶜ, a one-sided 95% upper bound at or below τ reads "harm of at most τ".

## Conditional-family reading

GRF's and DINA's candidate families are generated from fitted surfaces, so the fixed-family condition does not hold for them: **every `grf` and `dina` coverage row is coverage of the estimand conditional on the proposed family.** FS's family is the prespecified cut grid. Comparisons across the three identifiers are descriptive: the identifier, the family construction and the set of detected replicates all differ, so rows are read side by side and never ranked.
