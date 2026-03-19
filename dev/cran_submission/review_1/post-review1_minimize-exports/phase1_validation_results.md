# Phase 1 Validation Results: Vignette-Driven Export Classification

## Methodology

Every `.Rmd` and `.qmd` file was scanned programmatically.
Only function calls inside **executable R code blocks** (`` ```{r ... ``` ``)
and **inline R expressions** (`` `r func()` ``) were counted.
Functions appearing only in prose, comments, or `gt` table strings were excluded.
Functions from `weightedsurv` (`df_counting`, `plot_weighted_km`) and
functions defined inline within vignettes were excluded.

---

## Vignette-by-Vignette Findings

### forestsearch.Rmd (Getting Started)

18 forestsearch functions called:
`forestsearch()`, `forestsearch_bootstrap_dofuture()`,
`summarize_bootstrap_results()`, `summarize_bootstrap_events()`,
`forestsearch_Kfold()`, `forestsearch_tenfold()`, `forestsearch_KfoldOut()`,
`cv_summary_tables()`, `cv_metrics_tables()`, `grf.subg.harm.survival()`,
`sg_tables()`, `create_summary_table()`, `plot_subgroup_results_forestplot()`,
`render_forestplot()`, `create_forest_theme()`, `plot_km_band_forestsearch()`,
`plot_sg_weighted_km()`, `figure_note()`

### extreme_subgroups.Rmd

5 forestsearch functions called:
`gg_forest()`, `generate_aft_dgm_flex()`, `simulate_from_dgm()`,
`calibrate_cens_adjust()`, `check_censoring_dgm()`

### paper_simulations.Rmd

17 forestsearch functions called:
`create_summary_table()`, `plot_detection_curve()`,
`generate_aft_dgm_flex()`, `setup_gbsg_dgm()`, `compute_dgm_cde()`,
`simulate_from_dgm()`, `calibrate_k_inter()`, `validate_k_inter_effect()`,
`run_simulation_analysis()`, `summarize_simulation_results()`,
`format_oc_results()`, `build_classification_table()`,
`build_estimation_table()`, `interpret_estimation_table()`,
`render_reference_table()`, `compute_detection_probability()`,
`generate_detection_curve()`

### treatment_effect_definitions.Rmd

5 forestsearch functions called:
`cox_ahr_cde_analysis()`, `setup_gbsg_dgm()`, `compute_dgm_cde()`,
`simulate_from_dgm()`, `build_estimation_table()`

### biomarker_comparison.qmd

1 forestsearch function called:
`cox_cs_fit()`

### cox_causal_review.qmd

1 forestsearch function called:
`setup_gbsg_dgm()`

### biomarker_effects.Rmd — NO code calls (prose only)
### causal_effects_brief_review.Rmd — NO code calls (prose only)
### methodology.Rmd — NO code calls (prose only)

---

## Validated Tier 1 List (48 functions)

### Directly Called in Vignette Code (37 unique)

| # | Function | Vignette(s) |
|---|----------|-------------|
| 1 | `forestsearch()` | forestsearch, extreme_subgroups |
| 2 | `grf.subg.harm.survival()` | forestsearch |
| 3 | `forestsearch_bootstrap_dofuture()` | forestsearch |
| 4 | `summarize_bootstrap_results()` | forestsearch |
| 5 | `summarize_bootstrap_events()` | forestsearch |
| 6 | `forestsearch_Kfold()` | forestsearch |
| 7 | `forestsearch_tenfold()` | forestsearch |
| 8 | `forestsearch_KfoldOut()` | forestsearch |
| 9 | `cv_summary_tables()` | forestsearch |
| 10 | `cv_metrics_tables()` | forestsearch |
| 11 | `cox_ahr_cde_analysis()` | treatment_effect_definitions |
| 12 | `cox_cs_fit()` | biomarker_comparison |
| 13 | `sg_tables()` | forestsearch |
| 14 | `create_summary_table()` | forestsearch, paper_simulations |
| 15 | `gg_forest()` | extreme_subgroups |
| 16 | `plot_subgroup_results_forestplot()` | forestsearch |
| 17 | `render_forestplot()` | forestsearch |
| 18 | `create_forest_theme()` | forestsearch |
| 19 | `plot_km_band_forestsearch()` | forestsearch |
| 20 | `plot_sg_weighted_km()` | forestsearch |
| 21 | `plot_detection_curve()` | paper_simulations |
| 22 | `figure_note()` | forestsearch |
| 23 | `generate_aft_dgm_flex()` | extreme_subgroups, paper_simulations |
| 24 | `setup_gbsg_dgm()` | paper_simulations, treatment_effect_defs, cox_causal |
| 25 | `compute_dgm_cde()` | paper_simulations, treatment_effect_defs |
| 26 | `simulate_from_dgm()` | extreme_subgroups, paper_simulations, treatment_effect_defs |
| 27 | `calibrate_cens_adjust()` | extreme_subgroups |
| 28 | `check_censoring_dgm()` | extreme_subgroups |
| 29 | `calibrate_k_inter()` | paper_simulations |
| 30 | `validate_k_inter_effect()` | paper_simulations |
| 31 | `run_simulation_analysis()` | paper_simulations |
| 32 | `summarize_simulation_results()` | paper_simulations |
| 33 | `format_oc_results()` | paper_simulations |
| 34 | `build_classification_table()` | paper_simulations |
| 35 | `build_estimation_table()` | paper_simulations, treatment_effect_defs |
| 36 | `interpret_estimation_table()` | paper_simulations |
| 37 | `render_reference_table()` | paper_simulations |
|    | `compute_detection_probability()` | paper_simulations |
|    | `generate_detection_curve()` | paper_simulations |

### S3 Methods Required for Dispatch (8)

| Function | Class |
|----------|-------|
| `print.forestsearch()` | forestsearch |
| `summary.forestsearch()` | forestsearch |
| `plot.forestsearch()` | forestsearch |
| `print.fs_kfold()` | fs_kfold |
| `print.fs_tenfold()` | fs_tenfold |
| `print.cox_ahr_cde()` | cox_ahr_cde |
| `summary.cox_ahr_cde()` | cox_ahr_cde |
| `print.gbsg_dgm()` | gbsg_dgm |

### README / Configuration Entry Points (3)

| Function | Rationale |
|----------|-----------|
| `select_best_subgroup()` | README Quick Start workflow |
| `default_fs_params()` | Configuration entry point |
| `default_grf_params()` | Configuration entry point |

**Total Tier 1: 48 functions**

---

## Tier 2: Internal (~70 functions)

All remaining currently-exported functions become `@keywords internal`
with `@export` removed and `@examples` deleted entirely.

Major categories:

- **Subgroup evaluation internals** (10): `analyze_subgroup`, `assign_subgroup_membership`, `evaluate_subgroup_consistency`, `extract_subgroup`, `get_subgroup_membership`, `prepare_subgroup_data`, `sort_subgroups`, `sort_subgroups_preview`, `remove_near_duplicate_subgroups`, `remove_redundant_subgroups`
- **Consistency internals** (7): `subgroup.consistency`, `evaluate_consistency_twostage`, `run_single_consistency_split`, `setup_parallel_SGcons`, `sg_consistency_out`, `wilson_ci`, `early_stop_decision`
- **Bootstrap internals** (7): `bootstrap_results`, `bootstrap_ystar`, `count_boot_id`, `generate_bootstrap_synthetic`, `generate_bootstrap_with_noise`, `generate_gbsg_bootstrap_general`, `get_dfRes`
- **Bootstrap summary internals** (7): `summarize_bootstrap_subgroups`, `summarize_factor_presence_robust`, `format_bootstrap_table`, `format_bootstrap_diagnostics_table`, `format_bootstrap_timing_table`, `format_subgroup_summary_tables`, `create_factor_summary_tables`
- **CV internals** (3): `CV_sgs`, `cv_summary_text`, `cv_compare_results`
- **GRF internals** (8): `grf.subg.eval`, `fit_causal_forest`, `fit_policy_trees`, `create_grf_config`, `validate_grf_data`, `print_grf_details`, `compute_node_metrics`, `find_leaf_split`
- **Cox internals** (8): `cox_summary`, `cox_summary_batch`, `cox_summary_legacy`, `cox_summary_vectorized`, `build_cox_formula`, `fit_cox_models`, `get_Cox_sg`, `get_split_hr_fast`, `rmst_calculation`
- **Subgroup table internals** (7): `SG_tab_estimates`, `SGplot_estimates`, `FS_labels`, `create_summary_table_compact`, `create_summary_table_minimal`, `create_summary_table_presentation`, `create_summary_table_publication`, `create_sample_size_table`
- **Visualization internals** (14): `save_forestplot`, `print.fs_forest_theme`, `print.fs_forestplot`, `plot.fs_forestplot`, `create_subgroup_summary_df`, `quick_km_band_plot`, `plot_sg_results`, `plot_sg_distribution`, `print.fs_weighted_km`, `plot_subgroup`, `plot.fs_sg_plot`, `print.fs_sg_plot`, `plot_subgroup_effects`, `plot_spline_treatment_effect`, `compare_detection_curves`, `sens_text`, `km_summary`
- **Data prep internals** (13): `get_dfpred`, `dummy_encode`, `add_id_column`, `evaluate_comparison`, `evaluate_cuts_once`, `detect_variable_types`, `is_flag_continuous`, `is_flag_drop`, `is.continuous`, `get_cut_name`, `cut_var`, `lasso_selection`, `filter_by_lassokeep`
- **Tree cut internals** (4): `extract_all_tree_cuts`, `extract_selected_tree_cuts`, `extract_tree_cuts`, `extract_idx_flagredundancy`
- **DGM internals** (6): `create_gbsg_dgm`, `create_dgm_for_mrct`, `simulate_from_gbsg_dgm`, `get_dgm_with_output`, `find_k_inter_for_target_hr`, `sensitivity_analysis_k_inter`
- **Simulation internals** (3): `find_required_sample_size`, `create_null_result`, `create_success_result`
- **Formatting internals** (17): `format_CI`, `format_results`, `format_fs_details`, `hrCI_format`, `n_pcnt`, `ci_est`, `calc_cov`, `get_targetEst`, `calculate_counts`, `calculate_potential_hr`, `density_threshold_both`, `find_quantile_for_proportion`, `qlow`, `qhigh`, `get_best_survreg`, `compare_multiple_survreg`, `print.multi_survreg_comparison`
- **Core internals** (2): `subgroup.search`, `get_FSdata`
- **Internal utilities** (6): `filter_call_args`, `get_combinations_info`, `get_conf_force`, `get_covs_in`, `process_conf_force_expr`, `safe_eval_expr`
- **MRCT internals** (3): `mrct_region_sims`, `summaryout_mrct`, `validate_mrct_data`

---

## Decisions Requiring Larry's Input

| # | Function | Current | Recommendation | Question |
|---|----------|---------|----------------|----------|
| 1 | `subgroup.consistency()` | Export | Internal | Only called by `forestsearch()`. OK to internalize? |
| 2 | `create_gbsg_dgm()` | Export | Internal | Superseded by `setup_gbsg_dgm()`. OK to internalize? |
| 3 | `mrct_region_sims()` | Export | Internal | Not in any vignette. Keep for future MRCT vignette? |
| 4 | `create_dgm_for_mrct()` | Export | Internal | Same as above. |
| 5 | `plot_spline_treatment_effect()` | Export | Internal | Not in vignette code; `cox_cs_fit()` is the entry point. |
| 6 | `cox_summary()` | Export | Internal | Generally useful but not in vignettes. |
| 7 | `get_FSdata()` | Export | Internal | Data prep helper, not called directly. |
| 8 | `grf.subg.eval()` | Export | Internal | Companion to `grf.subg.harm.survival()`. |
| 9 | `format_CI()` | Export | Internal | Convenient but not in vignette code. |
| 10 | `simulate_from_gbsg_dgm()` | Export | Internal | Not in any vignette. `simulate_from_dgm()` used instead. |

---

## Next Steps

1. Larry reviews the 10 borderline decisions above
2. Proceed to Phase 2: Remove `@export` and delete `@examples` for all Tier 2 functions
3. Proceed to Phase 3: Fix examples for Tier 1 functions
4. Comprehensive `\dontrun{}` scan
5. `devtools::document()` → `devtools::check()`
