# Package index

## Core Algorithm

Main forestsearch function and subgroup search engine.

- [`forestsearch()`](https://larry-leon.github.io/forestsearch/reference/forestsearch.md)
  : ForestSearch: Exploratory Subgroup Identification
- [`print(`*`<forestsearch>`*`)`](https://larry-leon.github.io/forestsearch/reference/print.forestsearch.md)
  : Print Method for forestsearch Objects
- [`summary(`*`<forestsearch>`*`)`](https://larry-leon.github.io/forestsearch/reference/summary.forestsearch.md)
  : Summary Method for forestsearch Objects
- [`subgroup.search()`](https://larry-leon.github.io/forestsearch/reference/subgroup.search.md)
  : Subgroup Search for Treatment Effect Heterogeneity (Improved,
  Parallelized)
- [`get_FSdata()`](https://larry-leon.github.io/forestsearch/reference/get_FSdata.md)
  : ForestSearch Data Preparation and Feature Selection

## Subgroup Evaluation & Selection

Functions for evaluating, sorting, and selecting candidate subgroups.

- [`analyze_subgroup()`](https://larry-leon.github.io/forestsearch/reference/analyze_subgroup.md)
  : Analyze subgroup for summary table (OPTIMIZED)
- [`assign_subgroup_membership()`](https://larry-leon.github.io/forestsearch/reference/assign_subgroup_membership.md)
  : Assign data to subgroups based on selected node
- [`evaluate_subgroup_consistency()`](https://larry-leon.github.io/forestsearch/reference/evaluate_subgroup_consistency.md)
  : Evaluate Single Subgroup for Consistency (Fixed-Sample)
- [`extract_subgroup()`](https://larry-leon.github.io/forestsearch/reference/extract_subgroup.md)
  : Extract Subgroup Information
- [`get_subgroup_membership()`](https://larry-leon.github.io/forestsearch/reference/get_subgroup_membership.md)
  : Get subgroup membership vector
- [`prepare_subgroup_data()`](https://larry-leon.github.io/forestsearch/reference/prepare_subgroup_data.md)
  : Prepare subgroup data for analysis
- [`select_best_subgroup()`](https://larry-leon.github.io/forestsearch/reference/select_best_subgroup.md)
  : Select best subgroup based on criterion
- [`sort_subgroups()`](https://larry-leon.github.io/forestsearch/reference/sort_subgroups.md)
  : Sort Subgroups by Focus
- [`sort_subgroups_preview()`](https://larry-leon.github.io/forestsearch/reference/sort_subgroups_preview.md)
  : Sort Subgroups by Focus at consistency stage (consistency not
  available at this point)
- [`remove_near_duplicate_subgroups()`](https://larry-leon.github.io/forestsearch/reference/remove_near_duplicate_subgroups.md)
  : Remove Near-Duplicate Subgroups
- [`remove_redundant_subgroups()`](https://larry-leon.github.io/forestsearch/reference/remove_redundant_subgroups.md)
  : Remove Redundant Subgroups

## Consistency Evaluation

Split-sample consistency evaluation including the two-stage sequential
method.

- [`subgroup.consistency()`](https://larry-leon.github.io/forestsearch/reference/subgroup.consistency.md)
  : Evaluate Subgroup Consistency
- [`evaluate_consistency_twostage()`](https://larry-leon.github.io/forestsearch/reference/evaluate_consistency_twostage.md)
  : Evaluate Consistency (Two-Stage Algorithm)
- [`run_single_consistency_split()`](https://larry-leon.github.io/forestsearch/reference/run_single_consistency_split.md)
  : Run Single Consistency Split
- [`setup_parallel_SGcons()`](https://larry-leon.github.io/forestsearch/reference/setup_parallel_SGcons.md)
  : Set up parallel processing for subgroup consistency
- [`sg_consistency_out()`](https://larry-leon.github.io/forestsearch/reference/sg_consistency_out.md)
  : Output Subgroup Consistency Results
- [`wilson_ci()`](https://larry-leon.github.io/forestsearch/reference/wilson_ci.md)
  : Wilson Score Confidence Interval
- [`early_stop_decision()`](https://larry-leon.github.io/forestsearch/reference/early_stop_decision.md)
  : Early Stopping Decision

## Bootstrap Bias Correction

Bootstrap methods for bias-corrected hazard ratio estimation using
infinitesimal jackknife variance estimation (Leon et al., 2024).

- [`forestsearch_bootstrap_dofuture()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_bootstrap_dofuture.md)
  : ForestSearch Bootstrap with doFuture Parallelization
- [`bootstrap_results()`](https://larry-leon.github.io/forestsearch/reference/bootstrap_results.md)
  : Bootstrap Results for ForestSearch with Bias Correction
- [`bootstrap_ystar()`](https://larry-leon.github.io/forestsearch/reference/bootstrap_ystar.md)
  : Bootstrap Ystar Matrix
- [`count_boot_id()`](https://larry-leon.github.io/forestsearch/reference/count_boot_id.md)
  : Count ID Occurrences in Bootstrap Sample
- [`generate_bootstrap_synthetic()`](https://larry-leon.github.io/forestsearch/reference/generate_bootstrap_synthetic.md)
  : Generate Synthetic Data using Bootstrap with Perturbation
- [`generate_bootstrap_with_noise()`](https://larry-leon.github.io/forestsearch/reference/generate_bootstrap_with_noise.md)
  : Generate Bootstrap Sample with Added Noise
- [`generate_gbsg_bootstrap_general()`](https://larry-leon.github.io/forestsearch/reference/generate_gbsg_bootstrap_general.md)
  : Generate Synthetic GBSG Data using Generalized Bootstrap
- [`get_dfRes()`](https://larry-leon.github.io/forestsearch/reference/get_dfRes.md)
  : Bootstrap Confidence Interval and Bias Correction Results

## Bootstrap Summaries

Summarizing and formatting bootstrap results into publication-ready
tables and diagnostics.

- [`summarize_bootstrap_results()`](https://larry-leon.github.io/forestsearch/reference/summarize_bootstrap_results.md)
  : Enhanced Bootstrap Results Summary
- [`summarize_bootstrap_subgroups()`](https://larry-leon.github.io/forestsearch/reference/summarize_bootstrap_subgroups.md)
  : Summarize Bootstrap Subgroup Analysis Results
- [`summarize_bootstrap_events()`](https://larry-leon.github.io/forestsearch/reference/summarize_bootstrap_events.md)
  : Summarize Bootstrap Event Counts
- [`summarize_factor_presence_robust()`](https://larry-leon.github.io/forestsearch/reference/summarize_factor_presence_robust.md)
  : Summarize Factor Presence Across Bootstrap Subgroups
- [`format_bootstrap_table()`](https://larry-leon.github.io/forestsearch/reference/format_bootstrap_table.md)
  : Format Bootstrap Results Table with gt
- [`format_bootstrap_diagnostics_table()`](https://larry-leon.github.io/forestsearch/reference/format_bootstrap_diagnostics_table.md)
  : Format Bootstrap Diagnostics Table with gt
- [`format_bootstrap_timing_table()`](https://larry-leon.github.io/forestsearch/reference/format_bootstrap_timing_table.md)
  : Format Bootstrap Timing Table with gt
- [`format_subgroup_summary_tables()`](https://larry-leon.github.io/forestsearch/reference/format_subgroup_summary_tables.md)
  : Format Subgroup Summary Tables with gt
- [`create_factor_summary_tables()`](https://larry-leon.github.io/forestsearch/reference/create_factor_summary_tables.md)
  : Create Factor Summary Tables from Bootstrap Results

## Cross-Validation

K-fold and repeated cross-validation for assessing subgroup
identification stability and agreement metrics.

- [`forestsearch_Kfold()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_Kfold.md)
  : ForestSearch K-Fold Cross-Validation
- [`forestsearch_tenfold()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_tenfold.md)
  : ForestSearch Repeated K-Fold Cross-Validation
- [`forestsearch_KfoldOut()`](https://larry-leon.github.io/forestsearch/reference/forestsearch_KfoldOut.md)
  : ForestSearch K-Fold Cross-Validation Output Summary
- [`CV_sgs()`](https://larry-leon.github.io/forestsearch/reference/CV_sgs.md)
  : Cross-Validation Subgroup Match Summary
- [`cv_summary_tables()`](https://larry-leon.github.io/forestsearch/reference/cv_summary_tables.md)
  : Create Summary Tables from forestsearch_KfoldOut Results
- [`cv_metrics_tables()`](https://larry-leon.github.io/forestsearch/reference/cv_metrics_tables.md)
  : Create Metrics Tables for Cross-Validation Results
- [`cv_summary_text()`](https://larry-leon.github.io/forestsearch/reference/cv_summary_text.md)
  : Create Compact CV Summary Text
- [`cv_compare_results()`](https://larry-leon.github.io/forestsearch/reference/cv_compare_results.md)
  : Compare Multiple CV Results
- [`print(`*`<fs_kfold>`*`)`](https://larry-leon.github.io/forestsearch/reference/print.fs_kfold.md)
  : Print Method for K-Fold CV Results
- [`print(`*`<fs_tenfold>`*`)`](https://larry-leon.github.io/forestsearch/reference/print.fs_tenfold.md)
  : Print Method for Repeated K-Fold CV Results

## GRF Integration

Generalized Random Forest methods for heterogeneous treatment effect
estimation and variable importance screening.

- [`grf.subg.harm.survival()`](https://larry-leon.github.io/forestsearch/reference/grf.subg.harm.survival.md)
  : GRF Subgroup Identification for Survival Data
- [`grf.subg.eval()`](https://larry-leon.github.io/forestsearch/reference/grf.subg.eval.md)
  : GRF Subgroup Evaluation and Performance Metrics
- [`fit_causal_forest()`](https://larry-leon.github.io/forestsearch/reference/fit_causal_forest.md)
  : Fit causal survival forest
- [`fit_policy_trees()`](https://larry-leon.github.io/forestsearch/reference/fit_policy_trees.md)
  : Fit policy trees up to specified depth
- [`create_grf_config()`](https://larry-leon.github.io/forestsearch/reference/create_grf_config.md)
  : Helper Functions for GRF Subgroup Analysis
- [`validate_grf_data()`](https://larry-leon.github.io/forestsearch/reference/validate_grf_data.md)
  : Validate input data for GRF analysis
- [`print_grf_details()`](https://larry-leon.github.io/forestsearch/reference/print_grf_details.md)
  : Print detailed output for debugging
- [`compute_node_metrics()`](https://larry-leon.github.io/forestsearch/reference/compute_node_metrics.md)
  : Compute node metrics for a policy tree
- [`find_leaf_split()`](https://larry-leon.github.io/forestsearch/reference/find_leaf_split.md)
  : Find the split that leads to a specific leaf node

## Cox Model Utilities

Cox proportional hazards model wrappers with robust standard errors,
spline fitting, and average hazard ratio calculations.

- [`cox_summary()`](https://larry-leon.github.io/forestsearch/reference/cox_summary.md)
  : Cox model summary for subgroup (OPTIMIZED)
- [`cox_summary_batch()`](https://larry-leon.github.io/forestsearch/reference/cox_summary_batch.md)
  : Batch Cox summaries with caching
- [`cox_summary_legacy()`](https://larry-leon.github.io/forestsearch/reference/cox_summary_legacy.md)
  : Cox model summary for subgroup
- [`cox_summary_vectorized()`](https://larry-leon.github.io/forestsearch/reference/cox_summary_vectorized.md)
  : Cox model summary for subgroup - vectorized version
- [`cox_ahr_cde_analysis()`](https://larry-leon.github.io/forestsearch/reference/cox_ahr_cde_analysis.md)
  : Comprehensive Wrapper for Cox Spline Analysis with AHR and CDE
  Plotting
- [`print(`*`<cox_ahr_cde>`*`)`](https://larry-leon.github.io/forestsearch/reference/print.cox_ahr_cde.md)
  : Print method for cox_ahr_cde objects
- [`summary(`*`<cox_ahr_cde>`*`)`](https://larry-leon.github.io/forestsearch/reference/summary.cox_ahr_cde.md)
  : Summary method for cox_ahr_cde objects
- [`cox_cs_fit()`](https://larry-leon.github.io/forestsearch/reference/cox_cs_fit.md)
  : Fit Cox Model with Cubic Spline for Treatment Effect Heterogeneity
- [`build_cox_formula()`](https://larry-leon.github.io/forestsearch/reference/build_cox_formula.md)
  : Build Cox Model Formula
- [`fit_cox_models()`](https://larry-leon.github.io/forestsearch/reference/fit_cox_models.md)
  : Fit Cox Models for Subgroups
- [`get_Cox_sg()`](https://larry-leon.github.io/forestsearch/reference/get_Cox_sg.md)
  : Fit Cox Model for Subgroup
- [`get_split_hr_fast()`](https://larry-leon.github.io/forestsearch/reference/get_split_hr_fast.md)
  : Fast Cox Model HR Estimation
- [`rmst_calculation()`](https://larry-leon.github.io/forestsearch/reference/rmst_calculation.md)
  : RMST calculation for subgroup

## Subgroup Tables & Estimates

Formatted tables of subgroup estimates and labels.

- [`sg_tables()`](https://larry-leon.github.io/forestsearch/reference/sg_tables.md)
  : Enhanced Subgroup Summary Tables (gt output)
- [`SG_tab_estimates()`](https://larry-leon.github.io/forestsearch/reference/SG_tab_estimates.md)
  : Subgroup summary table estimates
- [`SGplot_estimates()`](https://larry-leon.github.io/forestsearch/reference/SGplot_estimates.md)
  : Violin/Boxplot Visualization of HR Estimates
- [`FS_labels()`](https://larry-leon.github.io/forestsearch/reference/FS_labels.md)
  : Convert Factor Code to Label
- [`create_summary_table()`](https://larry-leon.github.io/forestsearch/reference/create_summary_table.md)
  : Create Enhanced Summary Table for Baseline Characteristics
- [`create_summary_table_compact()`](https://larry-leon.github.io/forestsearch/reference/create_summary_table_compact.md)
  : Preset: Compact Table
- [`create_summary_table_minimal()`](https://larry-leon.github.io/forestsearch/reference/create_summary_table_minimal.md)
  : Preset: Minimal Table (No Highlighting, No Alternating)
- [`create_summary_table_presentation()`](https://larry-leon.github.io/forestsearch/reference/create_summary_table_presentation.md)
  : Preset: Presentation Table (Large Fonts)
- [`create_summary_table_publication()`](https://larry-leon.github.io/forestsearch/reference/create_summary_table_publication.md)
  : Preset: Publication-Ready Table
- [`create_sample_size_table()`](https://larry-leon.github.io/forestsearch/reference/create_sample_size_table.md)
  : Create Sample Size Table for Multiple Scenarios

## Visualization

Publication-ready plotting functions for forest plots, survival curves,
and subgroup characteristics.

- [`gg_forest()`](https://larry-leon.github.io/forestsearch/reference/gg_forest.md)
  : ggplot2 / patchwork forest plot
- [`plot_subgroup_results_forestplot()`](https://larry-leon.github.io/forestsearch/reference/plot_subgroup_results_forestplot.md)
  : Plot Subgroup Results Forest Plot
- [`render_forestplot()`](https://larry-leon.github.io/forestsearch/reference/render_forestplot.md)
  : Render ForestSearch Forest Plot
- [`save_forestplot()`](https://larry-leon.github.io/forestsearch/reference/save_forestplot.md)
  : Save ForestSearch Forest Plot to File
- [`create_forest_theme()`](https://larry-leon.github.io/forestsearch/reference/create_forest_theme.md)
  : Create Forest Plot Theme with Size Controls
- [`print(`*`<fs_forest_theme>`*`)`](https://larry-leon.github.io/forestsearch/reference/print.fs_forest_theme.md)
  : Print Method for ForestSearch Forest Theme
- [`print(`*`<fs_forestplot>`*`)`](https://larry-leon.github.io/forestsearch/reference/print.fs_forestplot.md)
  : Print Method for ForestSearch Forest Plot
- [`plot(`*`<fs_forestplot>`*`)`](https://larry-leon.github.io/forestsearch/reference/plot.fs_forestplot.md)
  : Plot Method for ForestSearch Forest Plot
- [`create_subgroup_summary_df()`](https://larry-leon.github.io/forestsearch/reference/create_subgroup_summary_df.md)
  : Create Subgroup Summary Data Frame for Forest Plot
- [`plot_km_band_forestsearch()`](https://larry-leon.github.io/forestsearch/reference/plot_km_band_forestsearch.md)
  : Plot Kaplan-Meier Survival Difference Bands for ForestSearch
  Subgroups
- [`quick_km_band_plot()`](https://larry-leon.github.io/forestsearch/reference/quick_km_band_plot.md)
  : Quick Plot KM Bands from ForestSearch
- [`plot_sg_results()`](https://larry-leon.github.io/forestsearch/reference/plot_sg_results.md)
  : Plot ForestSearch Subgroup Results
- [`plot_sg_distribution()`](https://larry-leon.github.io/forestsearch/reference/plot_sg_distribution.md)
  : Plot Distribution of Identified Subgroups
- [`plot_sg_weighted_km()`](https://larry-leon.github.io/forestsearch/reference/plot_sg_weighted_km.md)
  : Plot Weighted Kaplan-Meier Curves for ForestSearch Subgroups
- [`print(`*`<fs_weighted_km>`*`)`](https://larry-leon.github.io/forestsearch/reference/print.fs_weighted_km.md)
  : Print Method for fs_weighted_km Objects
- [`plot_subgroup()`](https://larry-leon.github.io/forestsearch/reference/plot_subgroup.md)
  : Plot Subgroup Survival Curves
- [`plot(`*`<fs_sg_plot>`*`)`](https://larry-leon.github.io/forestsearch/reference/plot.fs_sg_plot.md)
  : Plot Method for fs_sg_plot Objects
- [`print(`*`<fs_sg_plot>`*`)`](https://larry-leon.github.io/forestsearch/reference/print.fs_sg_plot.md)
  : Print Method for fs_sg_plot Objects
- [`plot_subgroup_effects()`](https://larry-leon.github.io/forestsearch/reference/plot_subgroup_effects.md)
  : Plot Subgroup Analysis Results
- [`plot(`*`<forestsearch>`*`)`](https://larry-leon.github.io/forestsearch/reference/plot.forestsearch.md)
  : Plot ForestSearch Results
- [`plot_spline_treatment_effect()`](https://larry-leon.github.io/forestsearch/reference/plot_spline_treatment_effect.md)
  : Plot Spline Treatment Effect Function
- [`plot_detection_curve()`](https://larry-leon.github.io/forestsearch/reference/plot_detection_curve.md)
  : Plot Detection Probability Curve
- [`compare_detection_curves()`](https://larry-leon.github.io/forestsearch/reference/compare_detection_curves.md)
  : Compare Detection Curves Across Sample Sizes
- [`sens_text()`](https://larry-leon.github.io/forestsearch/reference/sens_text.md)
  : Generate Cross-Validation Sensitivity Text
- [`figure_note()`](https://larry-leon.github.io/forestsearch/reference/figure_note.md)
  : Generate Figure Note for Quarto/RMarkdown
- [`km_summary()`](https://larry-leon.github.io/forestsearch/reference/km_summary.md)
  : KM median summary for subgroup

## Data Preparation & Encoding

Cut point generation, dummy variable creation, and variable selection.

- [`get_dfpred()`](https://larry-leon.github.io/forestsearch/reference/get_dfpred.md)
  : Generate Prediction Dataset with Subgroup Treatment Recommendation
- [`dummy_encode()`](https://larry-leon.github.io/forestsearch/reference/dummy_encode.md)
  : Dummy-code a data frame (numeric pass-through, factors expanded)
- [`add_id_column()`](https://larry-leon.github.io/forestsearch/reference/add_id_column.md)
  : Add ID Column to Data Frame
- [`evaluate_comparison()`](https://larry-leon.github.io/forestsearch/reference/evaluate_comparison.md)
  : Evaluate a Comparison Expression Without eval(parse())
- [`evaluate_cuts_once()`](https://larry-leon.github.io/forestsearch/reference/evaluate_cuts_once.md)
  : Cache and validate cut expressions efficiently
- [`detect_variable_types()`](https://larry-leon.github.io/forestsearch/reference/detect_variable_types.md)
  : Automatically Detect Variable Types in a Dataset
- [`is_flag_continuous()`](https://larry-leon.github.io/forestsearch/reference/is_flag_continuous.md)
  : Check if cut expression is for a continuous variable (OPTIMIZED)
- [`is_flag_drop()`](https://larry-leon.github.io/forestsearch/reference/is_flag_drop.md)
  : Check if cut expression should be dropped
- [`is.continuous()`](https://larry-leon.github.io/forestsearch/reference/is.continuous.md)
  : Check if a variable is continuous
- [`get_cut_name()`](https://larry-leon.github.io/forestsearch/reference/get_cut_name.md)
  : Get variable name from cut expression
- [`cut_var()`](https://larry-leon.github.io/forestsearch/reference/cut_var.md)
  : Generate cut expressions for a variable
- [`lasso_selection()`](https://larry-leon.github.io/forestsearch/reference/lasso_selection.md)
  : LASSO selection for Cox model
- [`filter_by_lassokeep()`](https://larry-leon.github.io/forestsearch/reference/filter_by_lassokeep.md)
  : Filter a vector by LASSO-selected variables

## Tree Cut Extraction

Extract and process cut points from tree-based models.

- [`extract_all_tree_cuts()`](https://larry-leon.github.io/forestsearch/reference/extract_all_tree_cuts.md)
  : Extract all cuts from fitted trees
- [`extract_selected_tree_cuts()`](https://larry-leon.github.io/forestsearch/reference/extract_selected_tree_cuts.md)
  : Extract cuts from selected tree only
- [`extract_tree_cuts()`](https://larry-leon.github.io/forestsearch/reference/extract_tree_cuts.md)
  : Extract cut information from a policy tree
- [`extract_idx_flagredundancy()`](https://larry-leon.github.io/forestsearch/reference/extract_idx_flagredundancy.md)
  : Extract redundancy flag for subgroup combinations

## Data Generating Mechanisms

Functions for creating data generating mechanisms (DGMs) for simulation
studies with configurable treatment effect heterogeneity.

- [`generate_aft_dgm_flex()`](https://larry-leon.github.io/forestsearch/reference/generate_aft_dgm_flex.md)
  : Generate Synthetic Survival Data using AFT Model with Flexible
  Subgroups
- [`setup_gbsg_dgm()`](https://larry-leon.github.io/forestsearch/reference/setup_gbsg_dgm.md)
  : Set Up a GBSG-Based AFT Data Generating Mechanism
- [`create_gbsg_dgm()`](https://larry-leon.github.io/forestsearch/reference/create_gbsg_dgm.md)
  : Create GBSG-Based AFT Data Generating Mechanism
- [`print(`*`<gbsg_dgm>`*`)`](https://larry-leon.github.io/forestsearch/reference/print.gbsg_dgm.md)
  : Print Method for gbsg_dgm Objects
- [`compute_dgm_cde()`](https://larry-leon.github.io/forestsearch/reference/compute_dgm_cde.md)
  : Compute and Attach CDE Values to a DGM Object
- [`create_dgm_for_mrct()`](https://larry-leon.github.io/forestsearch/reference/create_dgm_for_mrct.md)
  : Create Data Generating Mechanism for MRCT Simulations
- [`simulate_from_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_dgm.md)
  : Simulate Survival Data from AFT Data Generating Mechanism
- [`simulate_from_gbsg_dgm()`](https://larry-leon.github.io/forestsearch/reference/simulate_from_gbsg_dgm.md)
  : Simulate Trial Data from GBSG DGM
- [`get_dgm_with_output()`](https://larry-leon.github.io/forestsearch/reference/get_dgm_with_output.md)
  : Create DGM with Output File Path
- [`calibrate_cens_adjust()`](https://larry-leon.github.io/forestsearch/reference/calibrate_cens_adjust.md)
  : Calibrate Censoring Adjustment to Match DGM Reference Distribution
- [`check_censoring_dgm()`](https://larry-leon.github.io/forestsearch/reference/check_censoring_dgm.md)
  : Diagnose Censoring Consistency Between DGM Source Data and Simulated
  Data
- [`calibrate_k_inter()`](https://larry-leon.github.io/forestsearch/reference/calibrate_k_inter.md)
  : Calibrate k_inter for Target Subgroup Hazard Ratio
- [`find_k_inter_for_target_hr()`](https://larry-leon.github.io/forestsearch/reference/find_k_inter_for_target_hr.md)
  : Find k_inter Value to Achieve Target Harm Subgroup Hazard Ratio
- [`validate_k_inter_effect()`](https://larry-leon.github.io/forestsearch/reference/validate_k_inter_effect.md)
  : Validate k_inter Effect on HR Heterogeneity
- [`sensitivity_analysis_k_inter()`](https://larry-leon.github.io/forestsearch/reference/sensitivity_analysis_k_inter.md)
  : Sensitivity Analysis of Hazard Ratios to k_inter

## Simulation Studies

Running and summarizing simulation studies for operating
characteristics.

- [`run_simulation_analysis()`](https://larry-leon.github.io/forestsearch/reference/run_simulation_analysis.md)
  : Run One Simulation Replicate
- [`default_fs_params()`](https://larry-leon.github.io/forestsearch/reference/default_fs_params.md)
  : Default ForestSearch Parameters for GBSG Simulations
- [`default_grf_params()`](https://larry-leon.github.io/forestsearch/reference/default_grf_params.md)
  : Default GRF Parameters for GBSG Simulations
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
- [`find_required_sample_size()`](https://larry-leon.github.io/forestsearch/reference/find_required_sample_size.md)
  : Find Minimum Sample Size for Target Detection Power
- [`create_null_result()`](https://larry-leon.github.io/forestsearch/reference/create_null_result.md)
  : Create result object when no subgroup is found
- [`create_success_result()`](https://larry-leon.github.io/forestsearch/reference/create_success_result.md)
  : Create result object for successful subgroup identification

## Multi-Regional Clinical Trials

Specialized functions for multi-regional clinical trial (MRCT) subgroup
analysis with regional consistency evaluation.

- [`mrct_region_sims()`](https://larry-leon.github.io/forestsearch/reference/mrct_region_sims.md)
  : MRCT Regional Subgroup Simulation
- [`summaryout_mrct()`](https://larry-leon.github.io/forestsearch/reference/summaryout_mrct.md)
  : Summary Tables for MRCT Simulation Results
- [`validate_mrct_data()`](https://larry-leon.github.io/forestsearch/reference/validate_mrct_data.md)
  : Validate Dataset for MRCT Simulations

## Formatting & Estimation Helpers

Formatting utilities, confidence intervals, and estimation helpers.

- [`format_CI()`](https://larry-leon.github.io/forestsearch/reference/format_CI.md)
  : Format Confidence Interval for Estimates
- [`format_results()`](https://larry-leon.github.io/forestsearch/reference/format_results.md)
  : Format results for subgroup summary
- [`format_fs_details()`](https://larry-leon.github.io/forestsearch/reference/format_fs_details.md)
  : Format ForestSearch Details Output for Beamer Two-Column Display
- [`hrCI_format()`](https://larry-leon.github.io/forestsearch/reference/hrCI_format.md)
  : Format Hazard Ratio and Confidence Interval
- [`n_pcnt()`](https://larry-leon.github.io/forestsearch/reference/n_pcnt.md)
  : Calculate n and percent
- [`ci_est()`](https://larry-leon.github.io/forestsearch/reference/ci_est.md)
  : Confidence Interval for Estimate
- [`calc_cov()`](https://larry-leon.github.io/forestsearch/reference/calc_cov.md)
  : Calculate Covariance for Bootstrap Estimates
- [`get_targetEst()`](https://larry-leon.github.io/forestsearch/reference/get_targetEst.md)
  : Target Estimate and Standard Error for Bootstrap
- [`calculate_counts()`](https://larry-leon.github.io/forestsearch/reference/calculate_counts.md)
  : Calculate counts for subgroup summary
- [`calculate_potential_hr()`](https://larry-leon.github.io/forestsearch/reference/calculate_potential_hr.md)
  : Calculate potential outcome hazard ratio
- [`density_threshold_both()`](https://larry-leon.github.io/forestsearch/reference/density_threshold_both.md)
  : Bivariate Density for Split-Sample HR Threshold Detection
- [`find_quantile_for_proportion()`](https://larry-leon.github.io/forestsearch/reference/find_quantile_for_proportion.md)
  : Find Quantile for Target Subgroup Proportion
- [`qlow()`](https://larry-leon.github.io/forestsearch/reference/qlow.md)
  : 25th Percentile (Quantile Low)
- [`qhigh()`](https://larry-leon.github.io/forestsearch/reference/qhigh.md)
  : 75th Percentile (Quantile High)
- [`get_best_survreg()`](https://larry-leon.github.io/forestsearch/reference/get_best_survreg.md)
  : Get Best Model from Comparison
- [`compare_multiple_survreg()`](https://larry-leon.github.io/forestsearch/reference/compare_multiple_survreg.md)
  : Compare Multiple Survival Regression Models
- [`print(`*`<multi_survreg_comparison>`*`)`](https://larry-leon.github.io/forestsearch/reference/print.multi_survreg_comparison.md)
  : Print method for survreg_comparison objects

## Internal Utilities

Lower-level helpers. Most users will not call these directly.

- [`filter_call_args()`](https://larry-leon.github.io/forestsearch/reference/filter_call_args.md)
  : Filter and merge arguments for function calls
- [`get_combinations_info()`](https://larry-leon.github.io/forestsearch/reference/get_combinations_info.md)
  : Get all combinations of subgroup factors up to maxk
- [`get_conf_force()`](https://larry-leon.github.io/forestsearch/reference/get_conf_force.md)
  : Get forced cut expressions for variables
- [`get_covs_in()`](https://larry-leon.github.io/forestsearch/reference/get_covs_in.md)
  : Get indicator vector for selected subgroup factors
- [`process_conf_force_expr()`](https://larry-leon.github.io/forestsearch/reference/process_conf_force_expr.md)
  : Process forced cut expression for a variable
