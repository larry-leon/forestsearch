# Phase 1 Validation: Vignette-Driven Export Classification
# =========================================================
#
# Methodology: Every .Rmd and .qmd file was scanned for forestsearch
# function calls inside executable R code blocks (```{r ... ```).
# Inline R expressions (`r func()`) were also checked.
# Functions appearing only in prose, comments, or gt table strings
# were excluded.
#
# Functions from weightedsurv (df_counting, plot_weighted_km) and
# functions defined inline within vignettes were excluded.
# =========================================================


# ── DEFINITIVE: forestsearch functions called in executable code ─────────────

# forestsearch.Rmd (Getting Started vignette)
#   forestsearch()
#   forestsearch_bootstrap_dofuture()
#   summarize_bootstrap_results()
#   summarize_bootstrap_events()
#   forestsearch_Kfold()
#   forestsearch_tenfold()
#   forestsearch_KfoldOut()
#   cv_summary_tables()
#   cv_metrics_tables()
#   grf.subg.harm.survival()
#   sg_tables()
#   create_summary_table()
#   plot_subgroup_results_forestplot()
#   render_forestplot()
#   create_forest_theme()
#   plot_km_band_forestsearch()
#   plot_sg_weighted_km()
#   figure_note()                     [inline R expression]

# extreme_subgroups.Rmd
#   gg_forest()
#   generate_aft_dgm_flex()
#   simulate_from_dgm()
#   calibrate_cens_adjust()
#   check_censoring_dgm()

# paper_simulations.Rmd
#   create_summary_table()            [also in forestsearch.Rmd]
#   plot_detection_curve()
#   generate_aft_dgm_flex()           [also in extreme_subgroups.Rmd]
#   setup_gbsg_dgm()
#   compute_dgm_cde()
#   simulate_from_dgm()               [also in extreme_subgroups.Rmd]
#   calibrate_k_inter()
#   validate_k_inter_effect()
#   run_simulation_analysis()
#   summarize_simulation_results()
#   format_oc_results()
#   build_classification_table()
#   build_estimation_table()
#   interpret_estimation_table()
#   render_reference_table()
#   compute_detection_probability()
#   generate_detection_curve()

# treatment_effect_definitions.Rmd
#   cox_ahr_cde_analysis()
#   setup_gbsg_dgm()                  [also in paper_simulations.Rmd]
#   compute_dgm_cde()                 [also in paper_simulations.Rmd]
#   simulate_from_dgm()               [also above]
#   build_estimation_table()           [also in paper_simulations.Rmd]

# biomarker_comparison.qmd
#   cox_cs_fit()

# cox_causal_review.qmd
#   setup_gbsg_dgm()                  [also above]

# biomarker_effects.Rmd              — NO forestsearch code calls (prose only)
# causal_effects_brief_review.Rmd    — NO forestsearch code calls (prose only)
# methodology.Rmd                    — NO forestsearch code calls (prose only)


# ── UNIQUE FUNCTIONS CALLED IN EXECUTABLE CODE (deduplicated) ────────────────
# Total: 37 unique forestsearch functions

vignette_called <- c(
  # Core algorithm
  "forestsearch",
  "grf.subg.harm.survival",

  # Bootstrap
  "forestsearch_bootstrap_dofuture",
  "summarize_bootstrap_results",
  "summarize_bootstrap_events",


  # Cross-validation
  "forestsearch_Kfold",
  "forestsearch_tenfold",
  "forestsearch_KfoldOut",
  "cv_summary_tables",
  "cv_metrics_tables",

  # Cox utilities
  "cox_ahr_cde_analysis",
  "cox_cs_fit",

  # Subgroup tables
  "sg_tables",
  "create_summary_table",

  # Visualization
  "gg_forest",
  "plot_subgroup_results_forestplot",
  "render_forestplot",
  "create_forest_theme",
  "plot_km_band_forestsearch",
  "plot_sg_weighted_km",
  "plot_detection_curve",
  "figure_note",

  # DGM / Simulation
  "generate_aft_dgm_flex",
  "setup_gbsg_dgm",
  "compute_dgm_cde",
  "simulate_from_dgm",
  "calibrate_cens_adjust",
  "check_censoring_dgm",
  "calibrate_k_inter",
  "validate_k_inter_effect",
  "run_simulation_analysis",
  "summarize_simulation_results",
  "format_oc_results",
  "build_classification_table",
  "build_estimation_table",
  "interpret_estimation_table",
  "render_reference_table",
  "compute_detection_probability",
  "generate_detection_curve"
)


# ── S3 METHODS THAT MUST REMAIN EXPORTED (dispatch requires it) ──────────────
# Even though users don't call these explicitly, R's S3 dispatch
# requires them in NAMESPACE when the corresponding class is user-facing.

s3_methods_keep <- c(
  "print.forestsearch",
  "summary.forestsearch",
  "plot.forestsearch",
  "print.fs_kfold",
  "print.fs_tenfold",
  "print.cox_ahr_cde",
  "summary.cox_ahr_cde",
  "print.gbsg_dgm"
)


# ── ADDITIONAL EXPORTS: README Quick Start (not in vignettes) ────────────────
# The README shows select_best_subgroup() and plot_sg_weighted_km() (already
# counted). No additional unique functions beyond the vignette list.
# However, default_fs_params() and default_grf_params() are important
# configuration entry points even if not shown in current vignette code.

readme_additions <- c(
  "select_best_subgroup",      # README workflow step
  "default_fs_params",         # configuration entry point
  "default_grf_params"         # configuration entry point
)


# ── BORDERLINE: functions mentioned in prose but possibly needed ─────────────
# These are mentioned in explanatory text within vignettes. They are NOT
# called in code. Decision: keep as INTERNAL unless Larry promotes them.

prose_only_not_exported <- c(
  "run_fs_analysis",           # prose in treatment_effect_definitions.Rmd
  "run_grf_analysis",          # prose in treatment_effect_definitions.Rmd
  "plot_subgroup_effects",     # prose in treatment_effect_definitions.Rmd
  "get_dgm_hr",               # prose in paper_simulations.Rmd
  "format_CI",                 # not called in any vignette code
  "format_results"             # not called in any vignette code
)


# ══════════════════════════════════════════════════════════════════════════════
# FINAL TIER 1 (USER-FACING) LIST: 37 + 8 S3 + 3 README = 48 functions
# ══════════════════════════════════════════════════════════════════════════════

tier1_final <- sort(unique(c(vignette_called, s3_methods_keep, readme_additions)))
cat("TIER 1 count:", length(tier1_final), "\n")

# ══════════════════════════════════════════════════════════════════════════════
# TIER 2 (INTERNAL): Everything currently exported that is NOT in tier1_final
# Estimated: 118 - 48 = ~70 functions
# ══════════════════════════════════════════════════════════════════════════════


# ── TIER 2 FUNCTIONS BY CATEGORY (from _pkgdown.yml, minus Tier 1) ──────────

tier2_internal <- list(

  core_algorithm = c(
    "subgroup.search",         # called only by forestsearch() internally
    "get_FSdata"               # data preparation helper
  ),

  subgroup_evaluation = c(
    "analyze_subgroup",
    "assign_subgroup_membership",
    "evaluate_subgroup_consistency",
    "extract_subgroup",
    "get_subgroup_membership",
    "prepare_subgroup_data",
    # select_best_subgroup → TIER 1
    "sort_subgroups",
    "sort_subgroups_preview",
    "remove_near_duplicate_subgroups",
    "remove_redundant_subgroups"
  ),

  consistency = c(
    "subgroup.consistency",
    "evaluate_consistency_twostage",
    "run_single_consistency_split",
    "setup_parallel_SGcons",
    "sg_consistency_out",
    "wilson_ci",
    "early_stop_decision"
  ),

  bootstrap_internals = c(
    "bootstrap_results",
    "bootstrap_ystar",
    "count_boot_id",
    "generate_bootstrap_synthetic",
    "generate_bootstrap_with_noise",
    "generate_gbsg_bootstrap_general",
    "get_dfRes"
  ),

  bootstrap_summary_internals = c(
    # summarize_bootstrap_results → TIER 1
    # summarize_bootstrap_events → TIER 1
    "summarize_bootstrap_subgroups",
    "summarize_factor_presence_robust",
    "format_bootstrap_table",
    "format_bootstrap_diagnostics_table",
    "format_bootstrap_timing_table",
    "format_subgroup_summary_tables",
    "create_factor_summary_tables"
  ),

  cv_internals = c(
    # forestsearch_Kfold, _tenfold, _KfoldOut → TIER 1
    "CV_sgs",
    # cv_summary_tables, cv_metrics_tables → TIER 1
    "cv_summary_text",
    "cv_compare_results"
  ),

  grf_internals = c(
    # grf.subg.harm.survival → TIER 1
    "grf.subg.eval",
    "fit_causal_forest",
    "fit_policy_trees",
    "create_grf_config",
    "validate_grf_data",
    "print_grf_details",
    "compute_node_metrics",
    "find_leaf_split"
  ),

  cox_internals = c(
    "cox_summary",
    "cox_summary_batch",
    "cox_summary_legacy",
    "cox_summary_vectorized",
    # cox_ahr_cde_analysis → TIER 1
    # cox_cs_fit → TIER 1
    "build_cox_formula",
    "fit_cox_models",
    "get_Cox_sg",
    "get_split_hr_fast",
    "rmst_calculation"
  ),

  subgroup_table_internals = c(
    # sg_tables → TIER 1
    "SG_tab_estimates",
    "SGplot_estimates",
    "FS_labels",
    # create_summary_table → TIER 1
    "create_summary_table_compact",
    "create_summary_table_minimal",
    "create_summary_table_presentation",
    "create_summary_table_publication",
    "create_sample_size_table"
  ),

  viz_internals = c(
    # gg_forest → TIER 1
    # plot_subgroup_results_forestplot → TIER 1
    # render_forestplot → TIER 1
    "save_forestplot",
    # create_forest_theme → TIER 1
    "print.fs_forest_theme",
    "print.fs_forestplot",
    "plot.fs_forestplot",
    "create_subgroup_summary_df",
    # plot_km_band_forestsearch → TIER 1
    "quick_km_band_plot",
    "plot_sg_results",
    "plot_sg_distribution",
    # plot_sg_weighted_km → TIER 1
    "print.fs_weighted_km",
    "plot_subgroup",
    "plot.fs_sg_plot",
    "print.fs_sg_plot",
    "plot_subgroup_effects",
    # plot.forestsearch → TIER 1 (S3)
    "plot_spline_treatment_effect",
    # plot_detection_curve → TIER 1
    "compare_detection_curves",
    "sens_text",
    # figure_note → TIER 1
    "km_summary"
  ),

  data_prep = c(
    "get_dfpred",
    "dummy_encode",
    "add_id_column",
    "evaluate_comparison",
    "evaluate_cuts_once",
    "detect_variable_types",
    "is_flag_continuous",
    "is_flag_drop",
    "is.continuous",
    "get_cut_name",
    "cut_var",
    "lasso_selection",
    "filter_by_lassokeep"
  ),

  tree_cuts = c(
    "extract_all_tree_cuts",
    "extract_selected_tree_cuts",
    "extract_tree_cuts",
    "extract_idx_flagredundancy"
  ),

  dgm_internals = c(
    # generate_aft_dgm_flex → TIER 1
    # setup_gbsg_dgm → TIER 1
    "create_gbsg_dgm",         # NOTE: TIER 1 has print.gbsg_dgm but
    #                            create_gbsg_dgm is superseded by setup_gbsg_dgm
    # compute_dgm_cde → TIER 1
    "create_dgm_for_mrct",
    # simulate_from_dgm → TIER 1
    "simulate_from_gbsg_dgm",
    "get_dgm_with_output",
    # calibrate_cens_adjust → TIER 1
    # check_censoring_dgm → TIER 1
    # calibrate_k_inter → TIER 1
    "find_k_inter_for_target_hr",
    # validate_k_inter_effect → TIER 1
    "sensitivity_analysis_k_inter"
  ),

  simulation_internals = c(
    # run_simulation_analysis → TIER 1
    # default_fs_params, default_grf_params → TIER 1
    # summarize_simulation_results → TIER 1
    # format_oc_results → TIER 1
    # build_classification_table → TIER 1
    # build_estimation_table → TIER 1
    # interpret_estimation_table → TIER 1
    # render_reference_table → TIER 1
    # compute_detection_probability → TIER 1
    # generate_detection_curve → TIER 1
    "find_required_sample_size",
    "create_null_result",
    "create_success_result"
  ),

  formatting_internals = c(
    "format_CI",
    "format_results",
    "format_fs_details",
    "hrCI_format",
    "n_pcnt",
    "ci_est",
    "calc_cov",
    "get_targetEst",
    "calculate_counts",
    "calculate_potential_hr",
    "density_threshold_both",
    "find_quantile_for_proportion",
    "qlow",
    "qhigh",
    "get_best_survreg",
    "compare_multiple_survreg",
    "print.multi_survreg_comparison"
  ),

  internal_utilities = c(
    "filter_call_args",
    "get_combinations_info",
    "get_conf_force",
    "get_covs_in",
    "process_conf_force_expr",
    "safe_eval_expr"
  ),

  mrct_internals = c(
    "mrct_region_sims",
    "summaryout_mrct",
    "validate_mrct_data"
  )
)

tier2_count <- sum(sapply(tier2_internal, length))
cat("TIER 2 count:", tier2_count, "\n")
cat("TOTAL:", length(tier1_final) + tier2_count, "\n")


# ══════════════════════════════════════════════════════════════════════════════
# SPECIAL DECISIONS REQUIRING LARRY'S INPUT
# ══════════════════════════════════════════════════════════════════════════════

# 1. subgroup.consistency() — currently Tier 2.
#    It's listed prominently in README but not called directly in vignettes.
#    forestsearch() calls it internally. Keep internal?
#    RECOMMENDATION: Keep internal — forestsearch() is the user entry point.

# 2. create_gbsg_dgm() — superseded by setup_gbsg_dgm() per memory.
#    setup_gbsg_dgm() is Tier 1. create_gbsg_dgm() stays Tier 2?
#    RECOMMENDATION: Keep internal — but document "superseded" status.

# 3. simulate_from_gbsg_dgm() — not called in any vignette.
#    RECOMMENDATION: Internal.

# 4. mrct_region_sims() — listed in README as key feature but NOT called
#    in any current vignette code.
#    RECOMMENDATION: Keep Tier 2 for now. Promote if/when MRCT vignette added.

# 5. create_dgm_for_mrct() — same situation as mrct_region_sims().
#    RECOMMENDATION: Internal.

# 6. plot_spline_treatment_effect() — in README features table but not in
#    vignette code (cox_cs_fit() is the entry point in biomarker_comparison.qmd).
#    RECOMMENDATION: Internal (unless Larry wants it user-facing).

# 7. cox_summary() — a generally useful standalone function even if not
#    in current vignettes. However, the principle is "export only what
#    vignettes use."
#    RECOMMENDATION: Internal.

# 8. format_CI(), format_results() — convenient formatting utilities.
#    Not in vignette code.
#    RECOMMENDATION: Internal.

# 9. get_FSdata(), subgroup.search() — listed as "Core Algorithm" in pkgdown
#    but not called directly by users. forestsearch() is the entry point.
#    RECOMMENDATION: Internal.

# 10. grf.subg.eval() — companion to grf.subg.harm.survival() but not
#     called in vignettes.
#     RECOMMENDATION: Internal.
