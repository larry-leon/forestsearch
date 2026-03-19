# Package index

## Core Algorithm

Main forestsearch function, subgroup selection, and consistency
evaluation.

- [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  : ForestSearch: Exploratory Subgroup Identification
- [`print(`*`<forestsearch>`*`)`](https://larry-leon.github.io/forestsearch/reference/print.forestsearch.md)
  : Print Method for forestsearch Objects
- [`summary(`*`<forestsearch>`*`)`](https://larry-leon.github.io/forestsearch/reference/summary.forestsearch.md)
  : Summary Method for forestsearch Objects
- [`plot(`*`<forestsearch>`*`)`](https://larry-leon.github.io/forestsearch/reference/plot.forestsearch.md)
  : Plot ForestSearch Results
- [`select_best_subgroup()`](https://larry-leon.github.io/forestsearch/reference/select_best_subgroup.md)
  : Select best subgroup based on criterion
- [`subgroup.consistency()`](https://larry-leon.github.io/forestsearch/reference/subgroup.consistency.md)
  : Evaluate Subgroup Consistency

## Bootstrap Bias Correction

Bootstrap methods for bias-corrected hazard ratio estimation using
infinitesimal jackknife variance estimation (Leon et al., 2024).

- [`forestsearch_bootstrap_dofuture()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md)
  : ForestSearch Bootstrap with doFuture Parallelization
- [`summarize_bootstrap_results()`](https://larry-leon.github.io/forestsearch/reference/summarize_bootstrap_results.md)
  : Enhanced Bootstrap Results Summary
- [`summarize_bootstrap_events()`](https://larry-leon.github.io/forestsearch/reference/summarize_bootstrap_events.md)
  : Summarize Bootstrap Event Counts

## Cross-Validation

K-fold and repeated cross-validation for assessing subgroup
identification stability and agreement metrics.

- [`forestsearch_Kfold()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_Kfold.md)
  : ForestSearch K-Fold Cross-Validation
- [`forestsearch_tenfold()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_tenfold.md)
  : ForestSearch Repeated K-Fold Cross-Validation
- [`forestsearch_KfoldOut()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_KfoldOut.md)
  : ForestSearch K-Fold Cross-Validation Output Summary
- [`cv_summary_tables()`](https://larry-leon.github.io/forestsearch/reference/cv_summary_tables.md)
  : Create Summary Tables from forestsearch_KfoldOut Results
- [`cv_metrics_tables()`](https://larry-leon.github.io/forestsearch/reference/cv_metrics_tables.md)
  : Create Metrics Tables for Cross-Validation Results
- [`print(`*`<fs_kfold>`*`)`](https://larry-leon.github.io/forestsearch/reference/print.fs_kfold.md)
  : Print Method for K-Fold CV Results
- [`print(`*`<fs_tenfold>`*`)`](https://larry-leon.github.io/forestsearch/reference/print.fs_tenfold.md)
  : Print Method for Repeated K-Fold CV Results

## GRF Integration

Generalized Random Forest methods for heterogeneous treatment effect
estimation and variable importance screening.

- [`grf.subg.harm.survival()`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.survival.md)
  : GRF Subgroup Identification for Survival Data

## Cox Model Utilities

Cox proportional hazards model wrappers with robust standard errors,
spline fitting, and average hazard ratio calculations.

- [`cox_summary()`](https://larry-leon.github.io/forestsearch/reference/cox_summary.md)
  : Cox model summary for subgroup (OPTIMIZED)
- [`cox_ahr_cde_analysis()`](https://larry-leon.github.io/forestsearch/reference/cox_ahr_cde_analysis.md)
  : Comprehensive Wrapper for Cox Spline Analysis with AHR and CDE
  Plotting
- [`print(`*`<cox_ahr_cde>`*`)`](https://larry-leon.github.io/forestsearch/reference/print.cox_ahr_cde.md)
  : Print method for cox_ahr_cde objects
- [`summary(`*`<cox_ahr_cde>`*`)`](https://larry-leon.github.io/forestsearch/reference/summary.cox_ahr_cde.md)
  : Summary method for cox_ahr_cde objects
- [`cox_cs_fit()`](https://larry-leon.github.io/forestsearch/reference/cox_cs_fit.md)
  : Fit Cox Model with Cubic Spline for Treatment Effect Heterogeneity

## Subgroup Tables & Estimates

Formatted tables of subgroup estimates.

- [`sg_tables()`](https://larry-leon.github.io/forestsearch/reference/sg_tables.md)
  : Enhanced Subgroup Summary Tables (gt output)
- [`create_summary_table()`](https://larry-leon.github.io/forestsearch/reference/create_summary_table.md)
  : Create Enhanced Summary Table for Baseline Characteristics

## Visualization

Publication-ready plotting functions for forest plots, survival curves,
subgroup characteristics, and detection probability curves.

- [`gg_forest()`](https://larry-leon.github.io/forestsearch/reference/gg_forest.md)
  : ggplot2 / patchwork forest plot
- [`plot_subgroup_results_forestplot()`](https://larry-leon.github.io/forestsearch/reference/plot_subgroup_results_forestplot.md)
  : Plot Subgroup Results Forest Plot
- [`render_forestplot()`](https://larry-leon.github.io/forestsearch/reference/render_forestplot.md)
  : Render ForestSearch Forest Plot
- [`create_forest_theme()`](https://larry-leon.github.io/forestsearch/reference/create_forest_theme.md)
  : Create Forest Plot Theme with Size Controls
- [`plot_km_band_forestsearch()`](https://larry-leon.github.io/forestsearch/reference/plot_km_band_forestsearch.md)
  : Plot Kaplan-Meier Survival Difference Bands for ForestSearch
  Subgroups
- [`plot_sg_weighted_km()`](https://larry-leon.github.io/forestsearch/reference/plot_sg_weighted_km.md)
  : Plot Weighted Kaplan-Meier Curves for ForestSearch Subgroups
- [`figure_note()`](https://larry-leon.github.io/forestsearch/reference/figure_note.md)
  : Generate Figure Note for Quarto/RMarkdown
- [`plot_spline_treatment_effect()`](https://larry-leon.github.io/forestsearch/reference/plot_spline_treatment_effect.md)
  : Plot Spline Treatment Effect Function
- [`plot_detection_curve()`](https://larry-leon.github.io/forestsearch/reference/plot_detection_curve.md)
  : Plot Detection Probability Curve

## Data Generating Mechanisms

Functions for creating data generating mechanisms (DGMs) for simulation
studies with configurable treatment effect heterogeneity.

- [`generate_aft_dgm_flex()`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md)
  : Generate Synthetic Survival Data using AFT Model with Flexible
  Subgroups
- [`setup_gbsg_dgm()`](https://larry-leon.github.io/forestsearch/reference/setup_gbsg_dgm.md)
  : Set Up a GBSG-Based AFT Data Generating Mechanism
- [`print(`*`<gbsg_dgm>`*`)`](https://larry-leon.github.io/forestsearch/reference/print.gbsg_dgm.md)
  : Print Method for gbsg_dgm Objects
- [`compute_dgm_cde()`](https://larry-leon.github.io/forestsearch/reference/compute_dgm_cde.md)
  : Compute and Attach CDE Values to a DGM Object
- [`simulate_from_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md)
  : Simulate Survival Data from AFT Data Generating Mechanism
- [`calibrate_cens_adjust()`](https://larry-leon.github.io/forestsearch/reference/calibrate_cens_adjust.md)
  : Calibrate Censoring Adjustment to Match DGM Reference Distribution
- [`check_censoring_dgm()`](https://larry-leon.github.io/forestsearch/reference/check_censoring_dgm.md)
  : Diagnose Censoring Consistency Between DGM Source Data and Simulated
  Data
- [`calibrate_k_inter()`](https://larry-leon.github.io/forestsearch/reference/calibrate_k_inter.md)
  : Calibrate k_inter for Target Subgroup Hazard Ratio
- [`validate_k_inter_effect()`](https://larry-leon.github.io/forestsearch/reference/validate_k_inter_effect.md)
  : Validate k_inter Effect on HR Heterogeneity

## Simulation Studies

Running and summarizing simulation studies for operating
characteristics.

- [`run_simulation_analysis()`](https://larry-leon.github.io/forestsearch/reference/run_simulation_analysis.md)
  : Run One Simulation Replicate
- [`summarize_simulation_results()`](https://larry-leon.github.io/forestsearch/reference/summarize_simulation_results.md)
  : Summarize Simulation Results
- [`format_oc_results()`](https://larry-leon.github.io/forestsearch/reference/format_oc_results.md)
  : Format Operating Characteristics Results as GT Table
- [`build_classification_table()`](https://larry-leon.github.io/forestsearch/reference/build_classification_table.md)
  : Build Classification Rate Table from Simulation Results
- [`build_estimation_table()`](https://larry-leon.github.io/forestsearch/reference/build_estimation_table.md)
  : Build Estimation Properties Table from Simulation Results
- [`interpret_estimation_table()`](https://larry-leon.github.io/forestsearch/reference/interpret_estimation_table.md)
  : Generate Narrative Interpretation of Estimation Properties
- [`render_reference_table()`](https://larry-leon.github.io/forestsearch/reference/render_reference_table.md)
  : Render Reference Simulation Table as gt
- [`compute_detection_probability()`](https://larry-leon.github.io/forestsearch/reference/compute_detection_probability.md)
  : Compute Probability of Detecting True Subgroup
- [`generate_detection_curve()`](https://larry-leon.github.io/forestsearch/reference/generate_detection_curve.md)
  : Generate Detection Probability Curve

## Multi-Regional Clinical Trials

Specialized functions for multi-regional clinical trial (MRCT) subgroup
analysis with regional consistency evaluation.

- [`mrct_region_sims()`](https://larry-leon.github.io/forestsearch/reference/mrct_region_sims.md)
  : MRCT Regional Subgroup Simulation
- [`create_dgm_for_mrct()`](https://larry-leon.github.io/forestsearch/reference/create_dgm_for_mrct.md)
  : Create Data Generating Mechanism for MRCT Simulations

## Internal

Internal implementation details. These functions are not intended to be
called directly by users.

- [`CV_sgs()`](https://larry-leon.github.io/forestsearch/reference/CV_sgs.md)
  : Cross-Validation Subgroup Match Summary
- [`FS_labels()`](https://larry-leon.github.io/forestsearch/reference/FS_labels.md)
  : Convert Factor Code to Label
- [`SG_tab_estimates()`](https://larry-leon.github.io/forestsearch/reference/SG_tab_estimates.md)
  : Subgroup summary table estimates
- [`SGplot_estimates()`](https://larry-leon.github.io/forestsearch/reference/SGplot_estimates.md)
  : Violin/Boxplot Visualization of HR Estimates
- [`acm.disjctif()`](https://larry-leon.github.io/forestsearch/reference/acm.disjctif.md)
  : Disjunctive (dummy) coding for factor columns
- [`add_id_column()`](https://larry-leon.github.io/forestsearch/reference/add_id_column.md)
  : Add ID Column to Data Frame
- [`add_unprocessed_vars()`](https://larry-leon.github.io/forestsearch/reference/add_unprocessed_vars.md)
  : Add Unprocessed Variables from Original Data
- [`analyze_subgroup()`](https://larry-leon.github.io/forestsearch/reference/analyze_subgroup.md)
  : Analyze subgroup for summary table (OPTIMIZED)
- [`apply_spline_constraint()`](https://larry-leon.github.io/forestsearch/reference/apply_spline_constraint.md)
  : Apply Spline Constraint to Treatment Effect Coefficients
- [`assemble_results()`](https://larry-leon.github.io/forestsearch/reference/assemble_results.md)
  : Assemble Final Results Object
- [`assign_subgroup_membership()`](https://larry-leon.github.io/forestsearch/reference/assign_subgroup_membership.md)
  : Assign data to subgroups based on selected node
- [`bootstrap_results()`](https://larry-leon.github.io/forestsearch/reference/bootstrap_results.md)
  : Bootstrap Results for ForestSearch with Bias Correction
- [`bootstrap_ystar()`](https://larry-leon.github.io/forestsearch/reference/bootstrap_ystar.md)
  : Bootstrap Ystar Matrix
- [`build_cox_formula()`](https://larry-leon.github.io/forestsearch/reference/build_cox_formula.md)
  : Build Cox Model Formula
- [`calc_cov()`](https://larry-leon.github.io/forestsearch/reference/calc_cov.md)
  : Calculate Covariance for Bootstrap Estimates
- [`calculate_counts()`](https://larry-leon.github.io/forestsearch/reference/calculate_counts.md)
  : Calculate counts for subgroup summary
- [`calculate_event_counts()`](https://larry-leon.github.io/forestsearch/reference/calculate_event_counts.md)
  : Calculate Event Counts by Treatment Arm
- [`calculate_hazard_ratios()`](https://larry-leon.github.io/forestsearch/reference/calculate_hazard_ratios.md)
  : Calculate Hazard Ratios from Potential Outcomes
- [`calculate_linear_predictors()`](https://larry-leon.github.io/forestsearch/reference/calculate_linear_predictors.md)
  : Calculate Linear Predictors for Potential Outcomes
- [`calculate_max_combinations()`](https://larry-leon.github.io/forestsearch/reference/calculate_max_combinations.md)
  : Calculate Maximum Combinations
- [`calculate_potential_hr()`](https://larry-leon.github.io/forestsearch/reference/calculate_potential_hr.md)
  : Calculate potential outcome hazard ratio
- [`calculate_skewness()`](https://larry-leon.github.io/forestsearch/reference/calculate_skewness.md)
  : Calculate Skewness
- [`ci_est()`](https://larry-leon.github.io/forestsearch/reference/ci_est.md)
  : Confidence Interval for Estimate
- [`compare_detection_curves()`](https://larry-leon.github.io/forestsearch/reference/compare_detection_curves.md)
  : Compare Detection Curves Across Sample Sizes
- [`compare_multiple_survreg()`](https://larry-leon.github.io/forestsearch/reference/compare_multiple_survreg.md)
  : Compare Multiple Survival Regression Models
- [`compute_ahr()`](https://larry-leon.github.io/forestsearch/reference/compute_ahr.md)
  : Compute AHR from loghr_po
- [`compute_cde()`](https://larry-leon.github.io/forestsearch/reference/compute_cde.md)
  : Compute CDE from theta_0 and theta_1
- [`compute_detection_probability_single()`](https://larry-leon.github.io/forestsearch/reference/compute_detection_probability_single.md)
  : Compute Detection Probability for Single Theta (Internal)
- [`compute_node_metrics()`](https://larry-leon.github.io/forestsearch/reference/compute_node_metrics.md)
  : Compute node metrics for a policy tree
- [`compute_sg_hr()`](https://larry-leon.github.io/forestsearch/reference/compute_sg_hr.md)
  : Compute Hazard Ratio for a Single Subgroup
- [`compute_sg_hr_estimates()`](https://larry-leon.github.io/forestsearch/reference/compute_sg_hr_estimates.md)
  : Compute Hazard Ratio Estimates for Subgroups
- [`compute_sg_summary()`](https://larry-leon.github.io/forestsearch/reference/compute_sg_summary.md)
  : Compute Summary Statistics for Subgroups
- [`count_boot_id()`](https://larry-leon.github.io/forestsearch/reference/count_boot_id.md)
  : Count ID Occurrences in Bootstrap Sample
- [`cox_summary_batch()`](https://larry-leon.github.io/forestsearch/reference/cox_summary_batch.md)
  : Batch Cox summaries with caching
- [`cox_summary_vectorized()`](https://larry-leon.github.io/forestsearch/reference/cox_summary_vectorized.md)
  : Cox model summary for subgroup - vectorized version
- [`create_bootstrap_caption()`](https://larry-leon.github.io/forestsearch/reference/create_bootstrap_caption.md)
  : Calculate Bootstrap Table Caption
- [`create_bootstrap_diagnostic_plots()`](https://larry-leon.github.io/forestsearch/reference/create_bootstrap_diagnostic_plots.md)
  : Create Bootstrap Diagnostic Plots
- [`create_factor_summary_tables()`](https://larry-leon.github.io/forestsearch/reference/create_factor_summary_tables.md)
  : Create Factor Summary Tables from Bootstrap Results
- [`create_fs_subgroup_indicators()`](https://larry-leon.github.io/forestsearch/reference/create_fs_subgroup_indicators.md)
  : Create Subgroup Indicator Columns from ForestSearch
- [`create_gbsg_dgm()`](https://larry-leon.github.io/forestsearch/reference/create_gbsg_dgm.md)
  : Create GBSG-Based AFT Data Generating Mechanism
- [`create_grf_config()`](https://larry-leon.github.io/forestsearch/reference/create_grf_config.md)
  : Helper Functions for GRF Subgroup Analysis
- [`create_null_result()`](https://larry-leon.github.io/forestsearch/reference/create_null_result.md)
  : Create result object when no subgroup is found
- [`create_reference_subgroup_columns()`](https://larry-leon.github.io/forestsearch/reference/create_reference_subgroup_columns.md)
  : Create Reference Subgroup Indicator Columns
- [`create_result_row()`](https://larry-leon.github.io/forestsearch/reference/create_result_row.md)
  : Create Result Row
- [`create_sample_size_table()`](https://larry-leon.github.io/forestsearch/reference/create_sample_size_table.md)
  : Create Sample Size Table for Multiple Scenarios
- [`create_spline_variables()`](https://larry-leon.github.io/forestsearch/reference/create_spline_variables.md)
  : Create Spline Variables
- [`create_subgroup_indicator()`](https://larry-leon.github.io/forestsearch/reference/create_subgroup_indicator.md)
  : Create Subgroup Indicator from Factor Definitions
- [`create_subgroup_summary_df()`](https://larry-leon.github.io/forestsearch/reference/create_subgroup_summary_df.md)
  : Create Subgroup Summary Data Frame for Forest Plot
- [`create_success_result()`](https://larry-leon.github.io/forestsearch/reference/create_success_result.md)
  : Create result object for successful subgroup identification
- [`create_summary_table_compact()`](https://larry-leon.github.io/forestsearch/reference/create_summary_table_compact.md)
  : Preset: Compact Table
- [`create_summary_table_minimal()`](https://larry-leon.github.io/forestsearch/reference/create_summary_table_minimal.md)
  : Preset: Minimal Table (No Highlighting, No Alternating)
- [`create_summary_table_presentation()`](https://larry-leon.github.io/forestsearch/reference/create_summary_table_presentation.md)
  : Preset: Presentation Table (Large Fonts)
- [`create_summary_table_publication()`](https://larry-leon.github.io/forestsearch/reference/create_summary_table_publication.md)
  : Preset: Publication-Ready Table
- [`create_timing_summary_table()`](https://larry-leon.github.io/forestsearch/reference/create_timing_summary_table.md)
  : Create Timing Summary Table
- [`cut_numeric()`](https://larry-leon.github.io/forestsearch/reference/cut_numeric.md)
  : Discretize Continuous Variable into Quantile-Based Categories
- [`cut_size()`](https://larry-leon.github.io/forestsearch/reference/cut_size.md)
  : Discretize Continuous Variable by Size Categories
- [`cut_var()`](https://larry-leon.github.io/forestsearch/reference/cut_var.md)
  : Generate cut expressions for a variable
- [`cv_compare_results()`](https://larry-leon.github.io/forestsearch/reference/cv_compare_results.md)
  : Compare Multiple CV Results
- [`cv_summary_text()`](https://larry-leon.github.io/forestsearch/reference/cv_summary_text.md)
  : Create Compact CV Summary Text
- [`default_fs_params()`](https://larry-leon.github.io/forestsearch/reference/default_fs_params.md)
  : Default ForestSearch Parameters for GBSG Simulations
- [`default_grf_params()`](https://larry-leon.github.io/forestsearch/reference/default_grf_params.md)
  : Default GRF Parameters for GBSG Simulations
- [`default_grf_params_gen()`](https://larry-leon.github.io/forestsearch/reference/default_grf_params_gen.md)
  : Default GRF parameters (general)
- [`default_sim_params()`](https://larry-leon.github.io/forestsearch/reference/default_sim_params.md)
  : Default ForestSearch parameters (general)
- [`define_subgroups()`](https://larry-leon.github.io/forestsearch/reference/define_subgroups.md)
  : Define Subgroups with Flexible Cutpoints
- [`density_threshold_both()`](https://larry-leon.github.io/forestsearch/reference/density_threshold_both.md)
  : Bivariate Density for Split-Sample HR Threshold Detection
- [`density_threshold_integrand()`](https://larry-leon.github.io/forestsearch/reference/density_threshold_integrand.md)
  : Vectorized Density for Integration
- [`detect_variable_types()`](https://larry-leon.github.io/forestsearch/reference/detect_variable_types.md)
  : Automatically Detect Variable Types in a Dataset
- [`dummy_encode()`](https://larry-leon.github.io/forestsearch/reference/dummy_encode.md)
  : Dummy-code a data frame (numeric pass-through, factors expanded)
- [`early_stop_decision()`](https://larry-leon.github.io/forestsearch/reference/early_stop_decision.md)
  : Early Stopping Decision
- [`evaluate_combination_with_status()`](https://larry-leon.github.io/forestsearch/reference/evaluate_combination_with_status.md)
  : Evaluate a Single Factor Combination with Status Tracking
- [`evaluate_comparison()`](https://larry-leon.github.io/forestsearch/reference/evaluate_comparison.md)
  : Evaluate a Comparison Expression Without eval(parse())
- [`evaluate_consistency_twostage()`](https://larry-leon.github.io/forestsearch/reference/evaluate_consistency_twostage.md)
  : Evaluate Consistency (Two-Stage Algorithm)
- [`evaluate_cuts_once()`](https://larry-leon.github.io/forestsearch/reference/evaluate_cuts_once.md)
  : Cache and validate cut expressions efficiently
- [`evaluate_subgroup_consistency()`](https://larry-leon.github.io/forestsearch/reference/evaluate_subgroup_consistency.md)
  : Evaluate Single Subgroup for Consistency (Fixed-Sample)
- [`extract_all_tree_cuts()`](https://larry-leon.github.io/forestsearch/reference/extract_all_tree_cuts.md)
  : Extract all cuts from fitted trees
- [`extract_fs_estimates()`](https://larry-leon.github.io/forestsearch/reference/extract_fs_estimates.md)
  : Extract Estimates from ForestSearch Results
- [`extract_fs_subgroup_definition()`](https://larry-leon.github.io/forestsearch/reference/extract_fs_subgroup_definition.md)
  : Extract Subgroup Definition from ForestSearch Object
- [`extract_grf_estimates()`](https://larry-leon.github.io/forestsearch/reference/extract_grf_estimates.md)
  : Extract Estimates from GRF Results
- [`extract_idx_flagredundancy()`](https://larry-leon.github.io/forestsearch/reference/extract_idx_flagredundancy.md)
  : Extract redundancy flag for subgroup combinations
- [`extract_selected_tree_cuts()`](https://larry-leon.github.io/forestsearch/reference/extract_selected_tree_cuts.md)
  : Extract cuts from selected tree only
- [`extract_subgroup()`](https://larry-leon.github.io/forestsearch/reference/extract_subgroup.md)
  : Extract Subgroup Information
- [`extract_tree_cuts()`](https://larry-leon.github.io/forestsearch/reference/extract_tree_cuts.md)
  : Extract cut information from a policy tree
- [`filter_by_lassokeep()`](https://larry-leon.github.io/forestsearch/reference/filter_by_lassokeep.md)
  : Filter a vector by LASSO-selected variables
- [`filter_call_args()`](https://larry-leon.github.io/forestsearch/reference/filter_call_args.md)
  : Filter and merge arguments for function calls
- [`find_covariate_any_match()`](https://larry-leon.github.io/forestsearch/reference/find_covariate_any_match.md)
  : Find Covariate Any Match
- [`find_k_inter_for_target_hr()`](https://larry-leon.github.io/forestsearch/reference/find_k_inter_for_target_hr.md)
  : Find k_inter Value to Achieve Target Harm Subgroup Hazard Ratio
- [`find_leaf_split()`](https://larry-leon.github.io/forestsearch/reference/find_leaf_split.md)
  : Find the split that leads to a specific leaf node
- [`find_quantile_for_proportion()`](https://larry-leon.github.io/forestsearch/reference/find_quantile_for_proportion.md)
  : Find Quantile for Target Subgroup Proportion
- [`find_required_sample_size()`](https://larry-leon.github.io/forestsearch/reference/find_required_sample_size.md)
  : Find Minimum Sample Size for Target Detection Power
- [`fit_aft_model()`](https://larry-leon.github.io/forestsearch/reference/fit_aft_model.md)
  : Fit AFT Model with Optional Spline Treatment Effect
- [`fit_aft_model_spline()`](https://larry-leon.github.io/forestsearch/reference/fit_aft_model_spline.md)
  : Fit AFT Model with Spline Treatment Effect
- [`fit_aft_model_standard()`](https://larry-leon.github.io/forestsearch/reference/fit_aft_model_standard.md)
  : Fit Standard AFT Model (Non-Spline)
- [`fit_causal_forest()`](https://larry-leon.github.io/forestsearch/reference/fit_causal_forest.md)
  : Fit causal survival forest
- [`fit_cox_for_subgroup()`](https://larry-leon.github.io/forestsearch/reference/fit_cox_for_subgroup.md)
  : Fit Cox Model for Subgroup
- [`fit_cox_models()`](https://larry-leon.github.io/forestsearch/reference/fit_cox_models.md)
  : Fit Cox Models for Subgroups
- [`fit_policy_trees()`](https://larry-leon.github.io/forestsearch/reference/fit_policy_trees.md)
  : Fit policy trees up to specified depth
- [`forestsearch-package`](https://larry-leon.github.io/forestsearch/reference/forestsearch-package.md)
  [`forestsearch-imports`](https://larry-leon.github.io/forestsearch/reference/forestsearch-package.md)
  : forestsearch: Exploratory Subgroup Identification in Clinical Trials
  with Survival Endpoints
- [`format_CI()`](https://larry-leon.github.io/forestsearch/reference/format_CI.md)
  : Format Confidence Interval for Estimates
- [`format_bootstrap_diagnostics_table()`](https://larry-leon.github.io/forestsearch/reference/format_bootstrap_diagnostics_table.md)
  : Format Bootstrap Diagnostics Table with gt
- [`format_bootstrap_table()`](https://larry-leon.github.io/forestsearch/reference/format_bootstrap_table.md)
  : Format Bootstrap Results Table with gt
- [`format_bootstrap_timing_table()`](https://larry-leon.github.io/forestsearch/reference/format_bootstrap_timing_table.md)
  : Format Bootstrap Timing Table with gt
- [`format_continuous_definition()`](https://larry-leon.github.io/forestsearch/reference/format_continuous_definition.md)
  : Format Continuous Variable Definition for Display
- [`format_fs_details()`](https://larry-leon.github.io/forestsearch/reference/format_fs_details.md)
  : Format ForestSearch Details Output for Beamer Two-Column Display
- [`format_results()`](https://larry-leon.github.io/forestsearch/reference/format_results.md)
  : Format results for subgroup summary
- [`format_search_results()`](https://larry-leon.github.io/forestsearch/reference/format_search_results.md)
  : Format Search Results
- [`format_subgroup_summary_tables()`](https://larry-leon.github.io/forestsearch/reference/format_subgroup_summary_tables.md)
  : Format Subgroup Summary Tables with gt
- [`generate_bootstrap_synthetic()`](https://larry-leon.github.io/forestsearch/reference/generate_bootstrap_synthetic.md)
  : Generate Synthetic Data using Bootstrap with Perturbation
- [`generate_bootstrap_with_noise()`](https://larry-leon.github.io/forestsearch/reference/generate_bootstrap_with_noise.md)
  : Generate Bootstrap Sample with Added Noise
- [`generate_combination_indices()`](https://larry-leon.github.io/forestsearch/reference/generate_combination_indices.md)
  : Generate Combination Indices
- [`generate_complement_expression()`](https://larry-leon.github.io/forestsearch/reference/generate_complement_expression.md)
  : Generate Complement Expression
- [`generate_gbsg_bootstrap_general()`](https://larry-leon.github.io/forestsearch/reference/generate_gbsg_bootstrap_general.md)
  : Generate Synthetic GBSG Data using Generalized Bootstrap
- [`generate_readable_sg_labels()`](https://larry-leon.github.io/forestsearch/reference/generate_readable_sg_labels.md)
  : Generate Readable Subgroup Labels from ForestSearch Object
- [`generate_super_population()`](https://larry-leon.github.io/forestsearch/reference/generate_super_population.md)
  : Generate Super Population and Calculate Linear Predictors
- [`get_Cox_sg()`](https://larry-leon.github.io/forestsearch/reference/get_Cox_sg.md)
  : Fit Cox Model for Subgroup
- [`get_FSdata()`](https://larry-leon.github.io/forestsearch/reference/get_FSdata.md)
  : ForestSearch Data Preparation and Feature Selection
- [`get_best_survreg()`](https://larry-leon.github.io/forestsearch/reference/get_best_survreg.md)
  : Get Best Model from Comparison
- [`get_bootstrap_exports()`](https://larry-leon.github.io/forestsearch/reference/get_bootstrap_exports.md)
  : Get all exported functions from ForestSearch namespace
- [`get_combinations_info()`](https://larry-leon.github.io/forestsearch/reference/get_combinations_info.md)
  : Get all combinations of subgroup factors up to maxk
- [`get_conf_force()`](https://larry-leon.github.io/forestsearch/reference/get_conf_force.md)
  : Get forced cut expressions for variables
- [`get_covs_in()`](https://larry-leon.github.io/forestsearch/reference/get_covs_in.md)
  : Get indicator vector for selected subgroup factors
- [`get_cut_name()`](https://larry-leon.github.io/forestsearch/reference/get_cut_name.md)
  : Get variable name from cut expression
- [`get_dfRes()`](https://larry-leon.github.io/forestsearch/reference/get_dfRes.md)
  : Bootstrap Confidence Interval and Bias Correction Results
- [`get_dfpred()`](https://larry-leon.github.io/forestsearch/reference/get_dfpred.md)
  : Generate Prediction Dataset with Subgroup Treatment Recommendation
- [`get_dgm_hr()`](https://larry-leon.github.io/forestsearch/reference/get_dgm_hr.md)
  : Extract HR from DGM (Backward Compatible)
- [`get_dgm_with_output()`](https://larry-leon.github.io/forestsearch/reference/get_dgm_with_output.md)
  : Create DGM with Output File Path
- [`get_param()`](https://larry-leon.github.io/forestsearch/reference/get_param.md)
  : Get Parameter with Default Fallback
- [`get_split_hr_fast()`](https://larry-leon.github.io/forestsearch/reference/get_split_hr_fast.md)
  : Fast Cox Model HR Estimation
- [`get_subgroup_membership()`](https://larry-leon.github.io/forestsearch/reference/get_subgroup_membership.md)
  : Get subgroup membership vector
- [`get_targetEst()`](https://larry-leon.github.io/forestsearch/reference/get_targetEst.md)
  : Target Estimate and Standard Error for Bootstrap
- [`grf.subg.eval()`](https://larry-leon.github.io/forestsearch/reference/grf.subg.eval.md)
  : GRF Subgroup Evaluation and Performance Metrics
- [`has_positive_variance()`](https://larry-leon.github.io/forestsearch/reference/has_positive_variance.md)
  : Check if Matrix Has Positive Variance
- [`hrCI_format()`](https://larry-leon.github.io/forestsearch/reference/hrCI_format.md)
  : Format Hazard Ratio and Confidence Interval
- [`is.continuous()`](https://larry-leon.github.io/forestsearch/reference/is.continuous.md)
  : Check if a variable is continuous
- [`is_flag_continuous()`](https://larry-leon.github.io/forestsearch/reference/is_flag_continuous.md)
  : Check if cut expression is for a continuous variable (OPTIMIZED)
- [`is_flag_drop()`](https://larry-leon.github.io/forestsearch/reference/is_flag_drop.md)
  : Check if cut expression should be dropped
- [`km_summary()`](https://larry-leon.github.io/forestsearch/reference/km_summary.md)
  : KM median summary for subgroup
- [`lasso_selection()`](https://larry-leon.github.io/forestsearch/reference/lasso_selection.md)
  : LASSO selection for Cox model
- [`meets_event_criteria()`](https://larry-leon.github.io/forestsearch/reference/meets_event_criteria.md)
  : Check Event Count Criteria
- [`meets_prevalence_threshold()`](https://larry-leon.github.io/forestsearch/reference/meets_prevalence_threshold.md)
  : Check Prevalence Threshold
- [`n_pcnt()`](https://larry-leon.github.io/forestsearch/reference/n_pcnt.md)
  : Calculate n and percent
- [`parse_sg_harm_to_expression()`](https://larry-leon.github.io/forestsearch/reference/parse_sg_harm_to_expression.md)
  : Parse sg.harm Factor Names to Expression
- [`plot(`*`<fs_forestplot>`*`)`](https://larry-leon.github.io/forestsearch/reference/plot.fs_forestplot.md)
  : Plot Method for ForestSearch Forest Plot
- [`plot(`*`<fs_sg_plot>`*`)`](https://larry-leon.github.io/forestsearch/reference/plot.fs_sg_plot.md)
  : Plot Method for fs_sg_plot Objects
- [`plot_sg_distribution()`](https://larry-leon.github.io/forestsearch/reference/plot_sg_distribution.md)
  : Plot Distribution of Identified Subgroups
- [`plot_sg_forest()`](https://larry-leon.github.io/forestsearch/reference/plot_sg_forest.md)
  : Plot Forest Plot of Hazard Ratios
- [`plot_sg_km()`](https://larry-leon.github.io/forestsearch/reference/plot_sg_km.md)
  : Plot Kaplan-Meier Survival Curves for Subgroups
- [`plot_sg_results()`](https://larry-leon.github.io/forestsearch/reference/plot_sg_results.md)
  : Plot ForestSearch Subgroup Results
- [`plot_sg_summary_panel()`](https://larry-leon.github.io/forestsearch/reference/plot_sg_summary_panel.md)
  : Plot Summary Statistics Panel
- [`plot_subgroup()`](https://larry-leon.github.io/forestsearch/reference/plot_subgroup.md)
  : Plot Subgroup Survival Curves
- [`plot_subgroup_effects()`](https://larry-leon.github.io/forestsearch/reference/plot_subgroup_effects.md)
  : Plot Subgroup Analysis Results
- [`prepare_censoring_model()`](https://larry-leon.github.io/forestsearch/reference/prepare_censoring_model.md)
  : Prepare Censoring Model Parameters
- [`prepare_search_data()`](https://larry-leon.github.io/forestsearch/reference/prepare_search_data.md)
  : Prepare Data for Subgroup Search
- [`prepare_subgroup_data()`](https://larry-leon.github.io/forestsearch/reference/prepare_subgroup_data.md)
  : Prepare subgroup data for analysis
- [`prepare_working_dataset()`](https://larry-leon.github.io/forestsearch/reference/prepare_working_dataset.md)
  : Prepare Working Dataset with Processed Covariates
- [`print(`*`<fs_forest_theme>`*`)`](https://larry-leon.github.io/forestsearch/reference/print.fs_forest_theme.md)
  : Print Method for ForestSearch Forest Theme
- [`print(`*`<fs_forestplot>`*`)`](https://larry-leon.github.io/forestsearch/reference/print.fs_forestplot.md)
  : Print Method for ForestSearch Forest Plot
- [`print(`*`<fs_sg_plot>`*`)`](https://larry-leon.github.io/forestsearch/reference/print.fs_sg_plot.md)
  : Print Method for fs_sg_plot Objects
- [`print(`*`<fs_weighted_km>`*`)`](https://larry-leon.github.io/forestsearch/reference/print.fs_weighted_km.md)
  : Print Method for fs_weighted_km Objects
- [`print(`*`<multi_survreg_comparison>`*`)`](https://larry-leon.github.io/forestsearch/reference/print.multi_survreg_comparison.md)
  : Print method for survreg_comparison objects
- [`print_cv_params()`](https://larry-leon.github.io/forestsearch/reference/print_cv_params.md)
  : Print CV ForestSearch Parameters
- [`print_grf_details()`](https://larry-leon.github.io/forestsearch/reference/print_grf_details.md)
  : Print detailed output for debugging
- [`process_conf_force_expr()`](https://larry-leon.github.io/forestsearch/reference/process_conf_force_expr.md)
  : Process forced cut expression for a variable
- [`process_continuous_subgroup()`](https://larry-leon.github.io/forestsearch/reference/process_continuous_subgroup.md)
  : Process Continuous Variable for Subgroup Definition
- [`process_continuous_vars()`](https://larry-leon.github.io/forestsearch/reference/process_continuous_vars.md)
  : Process Continuous Variables
- [`process_cutpoint()`](https://larry-leon.github.io/forestsearch/reference/process_cutpoint.md)
  : Process Cutpoint Specification for Subgroup Definition
- [`process_factor_subgroup()`](https://larry-leon.github.io/forestsearch/reference/process_factor_subgroup.md)
  : Process Factor Variable for Subgroup Definition
- [`process_factor_vars()`](https://larry-leon.github.io/forestsearch/reference/process_factor_vars.md)
  : Process Factor Variables with Largest Value as Reference
- [`qhigh()`](https://larry-leon.github.io/forestsearch/reference/qhigh.md)
  : 75th Percentile (Quantile High)
- [`qlow()`](https://larry-leon.github.io/forestsearch/reference/qlow.md)
  : 25th Percentile (Quantile Low)
- [`quick_km_band_plot()`](https://larry-leon.github.io/forestsearch/reference/quick_km_band_plot.md)
  : Quick Plot KM Bands from ForestSearch
- [`remove_near_duplicate_subgroups()`](https://larry-leon.github.io/forestsearch/reference/remove_near_duplicate_subgroups.md)
  : Remove Near-Duplicate Subgroups
- [`remove_redundant_subgroups()`](https://larry-leon.github.io/forestsearch/reference/remove_redundant_subgroups.md)
  : Remove Redundant Subgroups
- [`resolve_bootstrap_parallel_args()`](https://larry-leon.github.io/forestsearch/reference/resolve_bootstrap_parallel_args.md)
  : Resolve parallel processing arguments for bootstrap
- [`resolve_cv_parallel_args()`](https://larry-leon.github.io/forestsearch/reference/resolve_cv_parallel_args.md)
  : Resolve Parallel Arguments for Cross-Validation
- [`rmst_calculation()`](https://larry-leon.github.io/forestsearch/reference/rmst_calculation.md)
  : RMST calculation for subgroup
- [`run_forestsearch_analysis()`](https://larry-leon.github.io/forestsearch/reference/run_forestsearch_analysis.md)
  : Run ForestSearch Analysis
- [`run_grf_analysis()`](https://larry-leon.github.io/forestsearch/reference/run_grf_analysis.md)
  : Run GRF Analysis
- [`run_single_consistency_split()`](https://larry-leon.github.io/forestsearch/reference/run_single_consistency_split.md)
  : Run Single Consistency Split
- [`safe_eval_expr()`](https://larry-leon.github.io/forestsearch/reference/safe_eval_expr.md)
  : Evaluate an expression string in a data-frame scope
- [`safe_subset()`](https://larry-leon.github.io/forestsearch/reference/safe_subset.md)
  : Subset a data frame using an expression string
- [`save_forestplot()`](https://larry-leon.github.io/forestsearch/reference/save_forestplot.md)
  : Save ForestSearch Forest Plot to File
- [`sens_text()`](https://larry-leon.github.io/forestsearch/reference/sens_text.md)
  : Generate Cross-Validation Sensitivity Text
- [`sensitivity_analysis_k_inter()`](https://larry-leon.github.io/forestsearch/reference/sensitivity_analysis_k_inter.md)
  : Sensitivity Analysis of Hazard Ratios to k_inter
- [`setup_parallel_SGcons()`](https://larry-leon.github.io/forestsearch/reference/setup_parallel_SGcons.md)
  : Set up parallel processing for subgroup consistency
- [`sg_consistency_out()`](https://larry-leon.github.io/forestsearch/reference/sg_consistency_out.md)
  : Output Subgroup Consistency Results
- [`simulate_from_gbsg_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_gbsg_dgm.md)
  : Simulate Trial Data from GBSG DGM
- [`sort_subgroups()`](https://larry-leon.github.io/forestsearch/reference/sort_subgroups.md)
  : Sort Subgroups by Focus
- [`sort_subgroups_preview()`](https://larry-leon.github.io/forestsearch/reference/sort_subgroups_preview.md)
  : Sort Subgroups by Focus at consistency stage (consistency not
  available at this point)
- [`subgroup.search()`](https://larry-leon.github.io/forestsearch/reference/subgroup.search.md)
  : Subgroup Search for Treatment Effect Heterogeneity (Improved,
  Parallelized)
- [`summarize_bootstrap_subgroups()`](https://larry-leon.github.io/forestsearch/reference/summarize_bootstrap_subgroups.md)
  : Summarize Bootstrap Subgroup Analysis Results
- [`summarize_factor_presence_robust()`](https://larry-leon.github.io/forestsearch/reference/summarize_factor_presence_robust.md)
  : Summarize Factor Presence Across Bootstrap Subgroups
- [`summarize_single_analysis()`](https://larry-leon.github.io/forestsearch/reference/summarize_single_analysis.md)
  : Summarize Single Analysis Results
- [`summaryout_mrct()`](https://larry-leon.github.io/forestsearch/reference/summaryout_mrct.md)
  : Summary Tables for MRCT Simulation Results
- [`validate_grf_data()`](https://larry-leon.github.io/forestsearch/reference/validate_grf_data.md)
  : Validate input data for GRF analysis
- [`validate_inputs()`](https://larry-leon.github.io/forestsearch/reference/validate_inputs.md)
  : Validate Input Parameters
- [`validate_mrct_data()`](https://larry-leon.github.io/forestsearch/reference/validate_mrct_data.md)
  : Validate Dataset for MRCT Simulations
- [`validate_spline_spec()`](https://larry-leon.github.io/forestsearch/reference/validate_spline_spec.md)
  : Validate Spline Specification
- [`wilson_ci()`](https://larry-leon.github.io/forestsearch/reference/wilson_ci.md)
  : Wilson Score Confidence Interval
