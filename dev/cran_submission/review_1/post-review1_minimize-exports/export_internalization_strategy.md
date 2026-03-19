# ForestSearch CRAN Review: Export Internalization Strategy

## Reviewer Concerns Summary

| ID  | Concern | Core Issue |
|-----|---------|------------|
| (A) | Comment-only `@examples` (descriptive stubs) | CRAN reviewer will flag as insufficient |
| (B) | Complex `\dontrun{}` examples for functions users never call | Wasted documentation effort |
| (C) | 118 exports, only ~42 used in vignettes | Too many user-facing functions |
| (D) | `\dontrun{}` still present in `R/oc_analyses_gbsg.R` | Comprehensive scan needed |

---

## Guiding Principle

> **Export only what an end user needs to call directly.**
> Everything else becomes internal (`@keywords internal`, no `@export`).
> Internal functions need only a title + `@description` + `@param` tags — no `@examples`.

---

## Step-by-Step Plan

### Step 1: Classify All 118 Exported Functions

Divide every currently exported function into one of two tiers:

**Tier 1 — User-facing (KEEP `@export`)**
Functions directly called in vignettes/articles, the README Quick Start,
or that an end user would reasonably need for a complete analysis.
These receive thorough documentation with real, runnable `@examples`.

**Tier 2 — Internal (REMOVE `@export`)**
Infrastructure called only by other package functions.
These receive minimal roxygen (title + description + params + `@keywords internal`).
All `@examples` sections are **completely removed**.

The classification is driven by two evidence sources:

1. **Vignettes/articles** — functions called in the 7 pkgdown articles
   (forestsearch, methodology, treatment_effect_definitions,
   biomarker_effects, causal_effects_brief_review, paper_simulations,
   extreme_subgroups)
2. **Quarto documents** — functions called in `biomarker_comparison.qmd`
   and `cox_causal_review.qmd`
3. **README Quick Start** — the canonical 4-step workflow

---

### Step 2: Proposed Tier 1 (User-Facing) Functions (~45-50)

Based on the README, vignettes, and Quarto documents, the following
functions should remain exported. This is the **proposed** list —
to be validated against actual vignette source code.

#### Core Algorithm (5)
- `forestsearch()`
- `print.forestsearch()`
- `summary.forestsearch()`
- `plot.forestsearch()`
- `get_FSdata()`

#### Subgroup Evaluation & Selection (3)
- `select_best_subgroup()`
- `subgroup.consistency()`
- `extract_subgroup()`

#### Bootstrap Inference (2)
- `forestsearch_bootstrap_dofuture()`
- `summarize_bootstrap_results()`

#### Cross-Validation (5)
- `forestsearch_Kfold()`
- `forestsearch_tenfold()`
- `cv_metrics_tables()`
- `cv_summary_tables()`
- `print.fs_kfold()` / `print.fs_tenfold()`

#### GRF Integration (2)
- `grf.subg.harm.survival()`
- `grf.subg.eval()`

#### Cox Model Utilities (4)
- `cox_summary()`
- `cox_cs_fit()`
- `cox_ahr_cde_analysis()`
- `print.cox_ahr_cde()` / `summary.cox_ahr_cde()`

#### Data Generating Mechanisms (5)
- `generate_aft_dgm_flex()`
- `setup_gbsg_dgm()`
- `create_gbsg_dgm()` / `print.gbsg_dgm()`
- `simulate_from_dgm()`
- `simulate_from_gbsg_dgm()`

#### Simulation Studies (3)
- `run_simulation_analysis()`
- `summarize_simulation_results()`
- `default_fs_params()` / `default_grf_params()`

#### MRCT (2)
- `mrct_region_sims()`
- `create_dgm_for_mrct()`

#### Visualization (8-10)
- `plot_subgroup_results_forestplot()`
- `plot_sg_weighted_km()`
- `plot_km_band_forestsearch()`
- `plot_spline_treatment_effect()`
- `plot_detection_curve()`
- `gg_forest()`
- `plot_sg_results()`
- `plot_subgroup_effects()`

#### Formatting (2-3)
- `format_CI()`
- `format_results()`
- `sg_tables()` / `SG_tab_estimates()`

---

### Step 3: Proposed Tier 2 (Internal) Functions (~70-75)

Everything NOT in the Tier 1 list above becomes internal. Major
categories include:

#### Subgroup Evaluation Internals
- `analyze_subgroup()`, `assign_subgroup_membership()`
- `evaluate_subgroup_consistency()`, `get_subgroup_membership()`
- `prepare_subgroup_data()`, `sort_subgroups()`, `sort_subgroups_preview()`
- `remove_near_duplicate_subgroups()`, `remove_redundant_subgroups()`

#### Consistency Evaluation Internals
- `evaluate_consistency_twostage()`, `run_single_consistency_split()`
- `setup_parallel_SGcons()`, `sg_consistency_out()`
- `wilson_ci()`, `early_stop_decision()`

#### Bootstrap Internals
- `bootstrap_results()`, `bootstrap_ystar()`, `count_boot_id()`
- `generate_bootstrap_synthetic()`, `generate_bootstrap_with_noise()`
- `generate_gbsg_bootstrap_general()`, `get_dfRes()`

#### Bootstrap Summary Internals
- `summarize_bootstrap_subgroups()`, `summarize_bootstrap_events()`
- `summarize_factor_presence_robust()`
- `format_bootstrap_table()`, `format_bootstrap_diagnostics_table()`
- `format_bootstrap_timing_table()`, `format_subgroup_summary_tables()`
- `create_factor_summary_tables()`

#### Cross-Validation Internals
- `forestsearch_KfoldOut()`, `CV_sgs()`
- `cv_summary_text()`, `cv_compare_results()`

#### GRF Internals
- `fit_causal_forest()`, `fit_policy_trees()`
- `create_grf_config()`, `validate_grf_data()`
- `print_grf_details()`, `compute_node_metrics()`, `find_leaf_split()`

#### Cox Model Internals
- `cox_summary_batch()`, `cox_summary_legacy()`, `cox_summary_vectorized()`
- `build_cox_formula()`, `fit_cox_models()`
- `get_Cox_sg()`, `get_split_hr_fast()`, `rmst_calculation()`

#### Subgroup Table Internals
- `SGplot_estimates()`, `FS_labels()`
- `create_summary_table()` (and _compact, _minimal, _presentation, _publication)
- `create_sample_size_table()`

#### Visualization Internals
- `render_forestplot()`, `save_forestplot()`
- `create_forest_theme()`, `print.fs_forest_theme()`
- `print.fs_forestplot()`, `plot.fs_forestplot()`
- `create_subgroup_summary_df()`
- `quick_km_band_plot()`, `plot_sg_distribution()`
- `print.fs_weighted_km()`, `plot_subgroup()`
- `plot.fs_sg_plot()`, `print.fs_sg_plot()`
- `compare_detection_curves()`, `sens_text()`, `figure_note()`, `km_summary()`

#### Data Preparation Internals
- `get_dfpred()`, `dummy_encode()`, `add_id_column()`
- `evaluate_comparison()`, `evaluate_cuts_once()`
- `detect_variable_types()`, `is_flag_continuous()`, `is_flag_drop()`
- `is.continuous()`, `get_cut_name()`, `cut_var()`
- `lasso_selection()`, `filter_by_lassokeep()`

#### Tree Cut Extraction Internals
- `extract_all_tree_cuts()`, `extract_selected_tree_cuts()`
- `extract_tree_cuts()`, `extract_idx_flagredundancy()`

#### DGM Internals
- `compute_dgm_cde()`, `get_dgm_with_output()`
- `calibrate_cens_adjust()`, `check_censoring_dgm()`
- `calibrate_k_inter()`, `find_k_inter_for_target_hr()`
- `validate_k_inter_effect()`, `sensitivity_analysis_k_inter()`

#### Simulation Internals
- `format_oc_results()`, `build_classification_table()`
- `build_estimation_table()`, `interpret_estimation_table()`
- `render_reference_table()`, `compute_detection_probability()`
- `generate_detection_curve()`, `find_required_sample_size()`
- `create_null_result()`, `create_success_result()`

#### Formatting Internals
- `format_fs_details()`, `hrCI_format()`, `n_pcnt()`
- `ci_est()`, `calc_cov()`, `get_targetEst()`
- `calculate_counts()`, `calculate_potential_hr()`
- `density_threshold_both()`, `find_quantile_for_proportion()`
- `qlow()`, `qhigh()`, `get_best_survreg()`
- `compare_multiple_survreg()`, `print.multi_survreg_comparison()`

#### Already-Identified Internal Utilities
- `filter_call_args()`, `get_combinations_info()`
- `get_conf_force()`, `get_covs_in()`, `process_conf_force_expr()`

#### MRCT Internals
- `summaryout_mrct()`, `validate_mrct_data()`

---

### Step 4: For Each Tier 2 Function — Specific Edits

For every function moving to Tier 2:

1. **Remove `@export`** from the roxygen block
2. **Add `@keywords internal`** if not already present
3. **Delete the entire `@examples` section** — no stubs, no `\dontrun{}`,
   no `\donttest{}`, nothing
4. **Keep** title, `@description`, `@param`, `@return`, and `@seealso`
   (but these can be minimal)
5. Optionally add `@noRd` if the function truly needs no man page at all
   (e.g., trivial one-liner helpers) — but `@keywords internal` is
   generally safer because it still generates a man page that `R CMD check`
   can validate

---

### Step 5: For Each Tier 1 Function — Fix Examples

For every function remaining in Tier 1:

1. **Remove all comment-only example stubs** (concern A)
2. **Replace `\dontrun{}` with real examples or `\donttest{}`** (concern D)
3. Apply the decision tree:
   - If the example is short and runs in < 5 seconds → **keep as runnable**
   - If the example is moderate and runs in 5-30 seconds → **wrap in `\donttest{}`**
   - If the example requires parallel workers, large data, or external
     resources → **wrap in `\donttest{}`** with a brief comment
   - **Never use `\dontrun{}`** — this is the hard rule from concern (D)

---

### Step 6: Comprehensive `\dontrun{}` Scan (Concern D)

After Steps 4-5, perform a global scan:

```r
grep -rn "\\\\dontrun" R/*.R
```

The result must be **zero matches**. Any remaining `\dontrun{}` blocks
must be converted to either:
- Real runnable examples (preferred)
- `\donttest{}` (for slow/resource-intensive examples)
- Removed entirely (for internal functions)

**Known location:** `R/oc_analyses_gbsg.R` has 4 `\dontrun{}` blocks.
These functions are likely Tier 2 (internal), so the examples should be
removed entirely.

---

### Step 7: Update NAMESPACE and _pkgdown.yml

1. Run `devtools::document()` to regenerate NAMESPACE from the updated
   `@export` tags
2. Update `_pkgdown.yml` to reorganize the reference:
   - Tier 1 functions get prominent placement in topic sections
   - Remove Tier 2 functions from explicit reference listings
     (they'll still appear under "Internal" if they have man pages)
3. Verify: `devtools::check()` must pass with 0 errors, 0 warnings

---

### Step 8: Final Verification

1. `grep -rn "\\\\dontrun" R/*.R` → 0 matches
2. `grep -rn "# Example:" R/*.R` (or similar stub patterns) → 0 matches
   in Tier 2 functions
3. Confirm NAMESPACE has ~45-50 exports (down from 118)
4. Confirm `man/` directory has ~45-50 user-facing pages + some
   `@keywords internal` pages
5. `devtools::check(args = "--as-cran")` → 0 errors, 0 warnings
6. All vignettes still build successfully (no broken references to
   now-internal functions — vignettes use `:::` or the functions are
   called indirectly via exported wrappers)

---

## Execution Order

| Phase | Action | Files Touched |
|-------|--------|---------------|
| 1 | Validate Tier 1 list against actual vignette `.Rmd` source | vignettes/*.Rmd |
| 2 | Remove `@export` + delete `@examples` for all Tier 2 functions | ~30-35 R/*.R files |
| 3 | Fix Tier 1 examples (remove stubs, convert `\dontrun{}`) | ~15-20 R/*.R files |
| 4 | Comprehensive `\dontrun{}` scan | All R/*.R files |
| 5 | `devtools::document()` → regenerate NAMESPACE | NAMESPACE |
| 6 | Update `_pkgdown.yml` reference sections | _pkgdown.yml |
| 7 | `devtools::check(args = "--as-cran")` | Full package |
| 8 | Win-builder testing | Tarball |

---

## Risk Mitigation

**Risk: Vignettes call now-internal functions directly**
- Mitigation: Vignettes should only call Tier 1 functions. If a vignette
  calls a Tier 2 function, either promote it to Tier 1 or refactor the
  vignette to use the exported wrapper.

**Risk: S3 methods for non-exported generics**
- S3 methods (`print.forestsearch`, `plot.forestsearch`, etc.) should
  remain exported even if the user doesn't call them directly — R's
  method dispatch requires they be in the NAMESPACE. Keep all `print.*`,
  `summary.*`, and `plot.*` methods that correspond to Tier 1 classes.

**Risk: Over-pruning breaks downstream code**
- The package is pre-CRAN (no reverse dependencies). The only consumers
  are the vignettes and Quarto documents. Both must be checked.

---

## Expected Outcome

| Metric | Before | After |
|--------|--------|-------|
| Exported functions | 118 | ~45-50 |
| Man pages | ~255 | ~60-70 (incl. internal) |
| `\dontrun{}` blocks | 4+ | 0 |
| Comment-only example stubs | Many | 0 |
| CRAN reviewer friction | High | Low |
