# Phase 2 Change Summary

## Statistics
- **81 functions** internalized (removed `@export`, added `@keywords internal`, deleted `@examples`)
- **32 \dontrun{}** blocks converted to `\donttest{}` in Tier 1 functions
- **53 functions** remain exported (Tier 1)
- **0 \dontrun{}** blocks remain
- **39 files** modified out of 48 total

## Tier 1 Exports (53 functions)

### Core Algorithm
`forestsearch()`, `print.forestsearch()`, `summary.forestsearch()`, `plot.forestsearch()`

### Subgroup Evaluation
`select_best_subgroup()`, `subgroup.consistency()`

### Bootstrap
`forestsearch_bootstrap_dofuture()`, `summarize_bootstrap_results()`, `summarize_bootstrap_events()`

### Cross-Validation
`forestsearch_Kfold()`, `forestsearch_tenfold()`, `forestsearch_KfoldOut()`, `cv_summary_tables()`, `cv_metrics_tables()`, `print.fs_kfold()`, `print.fs_tenfold()`

### GRF Integration
`grf.subg.harm.survival()`

### Cox Model Utilities
`cox_summary()`, `cox_ahr_cde_analysis()`, `print.cox_ahr_cde()`, `summary.cox_ahr_cde()`, `cox_cs_fit()`

### Subgroup Tables
`sg_tables()`, `create_summary_table()`

### Visualization
`gg_forest()`, `plot_subgroup_results_forestplot()`, `render_forestplot()`, `create_forest_theme()`, `plot_km_band_forestsearch()`, `plot_sg_weighted_km()`, `plot_detection_curve()`, `plot_spline_treatment_effect()`, `figure_note()`

### DGM / Simulation
`generate_aft_dgm_flex()`, `setup_gbsg_dgm()`, `print.gbsg_dgm()`, `compute_dgm_cde()`, `create_dgm_for_mrct()`, `simulate_from_dgm()`, `calibrate_cens_adjust()`, `check_censoring_dgm()`, `calibrate_k_inter()`, `validate_k_inter_effect()`, `run_simulation_analysis()`, `default_fs_params()`, `default_grf_params()`, `summarize_simulation_results()`, `format_oc_results()`, `build_classification_table()`, `build_estimation_table()`, `interpret_estimation_table()`, `render_reference_table()`, `compute_detection_probability()`, `generate_detection_curve()`

### MRCT
`mrct_region_sims()`

## Modified Files (39)

### Cox_estimation_helpers.R
- 3 function(s) internalized
- 3 \`dontrun{}\` → \`donttest{}\`

### bootstrap_analysis_dofuture.R
- 1 function(s) internalized
- 1 \`dontrun{}\` → \`donttest{}\`

### bootstrap_calculations_helpers.R
- 4 function(s) internalized
- 2 \`dontrun{}\` → \`donttest{}\`

### bootstrap_summaries_helpers.R
- 1 function(s) internalized
- 2 \`dontrun{}\` → \`donttest{}\`

### cox_spline_fit.R
- 1 function(s) internalized
- 2 \`dontrun{}\` → \`donttest{}\`

### create_summary_table.R
- 1 \`dontrun{}\` → \`donttest{}\`

### cv_summary_tables.R
- 2 function(s) internalized
- 4 \`dontrun{}\` → \`donttest{}\`

### find_k_inter_main.R
- 3 function(s) internalized

### forestsearch_cross_validation.R
- 1 function(s) internalized

### forestsearch_helpers.R
- 2 function(s) internalized
- 2 \`dontrun{}\` → \`donttest{}\`

### forestsearch_main.R
- 1 \`dontrun{}\` → \`donttest{}\`

### forestsearch_methods.R
- 3 \`dontrun{}\` → \`donttest{}\`

### format_fs_details.R
- 1 function(s) internalized
- 1 \`dontrun{}\` → \`donttest{}\`

### generate_aft_dgm_helpers.R
- 1 function(s) internalized
- 2 \`dontrun{}\` → \`donttest{}\`

### get_FSdata_helpers.R
- 6 function(s) internalized

### get_fsdata.R
- 1 function(s) internalized

### gg_forest.R
- 1 \`dontrun{}\` → \`donttest{}\`

### grf_helpers.R
- 8 function(s) internalized

### grf_main.R
- 1 function(s) internalized
- 2 \`dontrun{}\` → \`donttest{}\`

### mrct_simulation.R
- 3 function(s) internalized
- 5 \`dontrun{}\` → \`donttest{}\`

### oc_analyses_gbsg.R
- 4 \`dontrun{}\` → \`donttest{}\`

### plot_km_band_forestsearch.R
- 1 function(s) internalized
- 2 \`dontrun{}\` → \`donttest{}\`

### plot_sg_distribution.R
- 1 function(s) internalized
- 1 \`dontrun{}\` → \`donttest{}\`

### plot_sg_results.R
- 3 function(s) internalized

### plot_sg_weighted_km.R
- 1 function(s) internalized

### plot_subgroup_results_forestplot.R
- 1 \`dontrun{}\` → \`donttest{}\`

### render_forestplot.R
- 4 function(s) internalized

### run_simulation_analysis.R
- 1 \`dontrun{}\` → \`donttest{}\`

### setup_gbsg_dgm.R
- 1 \`dontrun{}\` → \`donttest{}\`

### sim_aft_gbsg.R
- 2 function(s) internalized

### simulate_from_dgm.R
- 3 \`dontrun{}\` → \`donttest{}\`

### simulation_tables.R
- 4 \`dontrun{}\` → \`donttest{}\`

### subgroup_consistency_helpers.R
- 10 function(s) internalized

### subgroup_consistency_main.R
- 1 \`dontrun{}\` → \`donttest{}\`

### subgroup_search.R
- 4 function(s) internalized
- 1 \`dontrun{}\` → \`donttest{}\`

### summarize_bootstrap_results.R
- 1 \`dontrun{}\` → \`donttest{}\`

### summary_utility_functions.R
- 7 function(s) internalized
- 5 \`dontrun{}\` → \`donttest{}\`

### synthetic_data_perturbation.R
- 3 function(s) internalized

### truefind_asymptotic.R
- 4 function(s) internalized

## Next Steps

1. Copy all 39 modified `.R` files to `R/` in your local repo
2. Run `devtools::document()` to regenerate NAMESPACE
3. Run `devtools::load_all()`
4. Run `devtools::check(args = "--as-cran")`
5. Update `_pkgdown.yml` reference sections (Phase 3)
6. Win-builder testing
