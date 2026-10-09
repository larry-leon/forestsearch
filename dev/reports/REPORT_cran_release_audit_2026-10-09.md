# REPORT — CRAN release audit, Stage 0 (2026-10-09)

- Task: `dev/tasks/TASK_cran_release_audit_2026-10-09.md` (committed alone at `541fef30`).
- Audit run at HEAD `541fef30c92d272c84178ad09925fa9259e58828` on `feature/glm-extension`.
- Type: read-only audit. No file under `R/`, `DESCRIPTION`, `NAMESPACE`, `man/`, `tests/`, `vignettes/`, `README*`, `NEWS.md` or `.Rbuildignore` was edited. Nothing was installed and nothing was uploaded.
- All build and check work ran in the session scratchpad, outside the repo. Those directories were deleted after extraction.

## 1. Summary

### Blockers (CRAN would reject or a reviewer would return the submission)

- **WARNING `checking Rd \usage sections`:** 4 exported functions in `R/mrct_inject_dgm.R` have no `@param` (undocumented arguments) and no `@return`:
  - `analysis_covariates` (:425)
  - `mrct_truth_by_region` (:443)
  - `plot_region_modifier` (:474)
  - `sim_region_metrics` (:502)
  - These are new since the 2026-09-18 reference (module added 2026-10-02).
- **Missing `\value`:** 11 exported Rd files.
  - The 4 above, plus `print_grf_details` (`R/grf_helpers.R:625`).
  - 6 S3 method Rd files:
    - `plot.fs_sg_plot` (`R/plot_sg_results.R:914`)
    - `print.fs_sg_plot` (`R/plot_sg_results.R:871`)
    - `print.forestsearch_comparison` (`R/compare_selection_rules.R:811`)
    - `print.fs_kfold` (`R/forestsearch_cross_validation.R:1984`)
    - `print.fs_tenfold` (`R/forestsearch_cross_validation.R:2011`)
    - `print.gbsg_dgm` (`R/sim_aft_gbsg.R:1155`)
- **PDF manual not certified (ERROR + WARNING, environment-caused):** `pdflatex is not available`.
  - Rd→LaTeX conversion completed, but the PDF was never compiled, so the LaTeX safety of the Rd files is **unverified** on this run.
  - `~/Library/TinyTeX` exists but is not on `PATH`. I did not fix this, per the task's "do not attempt to fix" rule for PATH tools.
- **Possible more than 2 cores under check:** the `\donttest` example of `run_null_calibration` used 74.8 s user CPU against 27.3 s elapsed. That is a ratio of **2.74**, above CRAN's 2.5 CPU-to-elapsed NOTE threshold.
  - This happened even though the example sets `plan = "sequential"`, `workers = 1`.
  - Cause not verified. Implicit threading inside a dependency (the example sets `use_grf = TRUE`) is a hypothesis only.
- **Defaults based on `parallel::detectCores()`, uncapped under `_R_CHECK_LIMIT_CORES_`, in exported entry points.** None of these was exercised by the check's examples or tests:
  - `subgroup.search(parallel_workers = parallel::detectCores())` (`R/subgroup_search.R:101`). It also calls `future::plan(multisession, ...)` at :164 without restoring the previous plan.
  - `mrct_region_sims()`: `workers <- max(1, detectCores() - 1)` when unset (`R/mrct_simulation.R:266`).
  - `fs_fdr_report()`: `n_workers` defaults to `0.90 * (detectCores(logical = FALSE) - 1)` (`R/fs_fdr_report.R:219`).
  - The internal `resolve_cv_parallel_args()` defaults to `detectCores() - 1` (`R/forestsearch_cross_validation.R:1928`).
- **Default writes outside `tempdir()`:**
  - `save_leff_calibration(output_dir = NULL)` resolves to `"_data"` relative to the working directory, then runs `dir.create()` and `saveRDS()` there (`R/leff_calibration_io.R:66-71, 99`).
  - `cox_ahr_cde_analysis()`'s `output_dir` default changed from `tempdir()` (0.1.0) to `"plots"` (`R/cox_ahr_cde_wrapper.R:122`). It only writes when `save_plots = TRUE`, which is off by default (`:373-374`).
- **`Title` is 66 characters,** one over the 65-character limit.
- **No `SystemRequirements: Quarto command line tools`,** although `vignettes/survival-post-selection.qmd` uses `%\VignetteEngine{quarto::html}` and `VignetteBuilder` lists `quarto`.
- **Not a blocker (checked):**
  - No `Remotes:` field is present.
  - All 3 `weightedsurv` functions forestsearch uses exist in CRAN weightedsurv 0.1.0.

### NOTE-level items

- CRAN incoming NOTE: `Version contains large components (0.3.5.9002)`. Also "VignetteBuilder field but no prebuilt vignette index", which comes from the `--no-build-vignettes` build.
- `checking top-level files` NOTE: README/NEWS cannot be checked without pandoc (environment).
- `checking HTML version of manual` NOTE: `tidy` is too old (environment; same as 2026-09-18).
- `non-standard things in the check directory` NOTE: `forestsearch-manual.tex`, left over from the failed PDF step.
- Suggests `latentcor` is not installed locally. The first check run stopped with `ERROR` at `checking package dependencies`; I reran with `_R_CHECK_FORCE_SUGGESTS_=false` (see §8).
- 145 Rd files use `\dontrun{}`, and every one of them wraps the whole example in it. 141 of these document exported symbols.
- 39 exported functions and 22 exported S3 methods have no `\examples` at all.
- `bootstrap_dofuture_main.R:347` sets `options(warn = -1)` for the parallel section. It is restored via `on.exit`.
- `plot_subgroup()` (internal) at `R/subgroup_consistency_helpers.R:968` sets `par(mfrow = c(1, 2))` and resets it to `c(1, 1)` instead of restoring the old `par`.
- `set.seed()` appears at 44 sites in package code. All are seeded from a user argument or an object field; 16 are guarded by `if (!is.null(seed))` and 28 are unconditional.
- `inst/patches/2026-05-21_extract_grf_args.md` (a dev note) ships in the tarball.
- `run_null_calibration` example takes 27.3 s elapsed and `suggest_thresholds` 6.0 s; both are `\donttest`.
- Positional-argument risk: new formals were inserted mid-signature in `forestsearch()` and `subgroup.consistency()` (see §3). Not checked further.

## 2. Repo state

- P1 — `git status --porcelain` before copying the task file: *(empty; clean)*. After copying, before the commit: `?? dev/tasks/TASK_cran_release_audit_2026-10-09.md`.
- P2 — Branch `feature/glm-extension`.
  - HEAD before the task commit: `c8d99093`.
  - Task-doc commit and **audit HEAD: `541fef30c92d272c84178ad09925fa9259e58828`**.
  - `DESCRIPTION` `Version: 0.3.5.9002`.
- P3 — `pandoc` and `quarto` are **not on `PATH`** (`command not found` for both).
  - Bundled copies exist at `/Applications/RStudio.app/Contents/Resources/app/quarto/bin/quarto` and `.../quarto/bin/tools/{aarch64,x86_64}/pandoc`. They were not used.
  - Per P3, the vignette build was skipped: built with `--no-build-vignettes`, checked with `--ignore-vignettes`.
  - `pdflatex` is also not on `PATH`; `~/Library/TinyTeX` exists.
- Environment: R 4.5.2 (aarch64-apple-darwin20), macOS 27.0.1, roxygen2 8.1.0, 14 logical cores.
- The CLAUDE.md pandoc hint (`/usr/lib/rstudio/...x86_64`) is a Linux path and does not apply on this machine.

## 3. Release delta since 0.1.0

### How 0.1.0 was identified

- Lightweight tag `v0.1.0` → `29890c7c` (2026-03-17, "CRAN check (finally) passed"). `DESCRIPTION` there reads `Version: 0.1.0`.
- **The tag is not an ancestor of HEAD.**
  - It lives on `remotes/origin/cran-review-1`.
  - The merge-base with HEAD is `18e13d0d`.
- All deltas below are direct tree comparisons (`v0.1.0..HEAD`), not branch history.

### Breaking changes first

- **Removed exports (4).** All 4 functions still exist in `R/` as internal (unexported) functions. There are no renames.
  - `create_sample_size_table` (`R/truefind_asymptotic.R:712`)
  - `detect_variable_types` (`R/synthetic_data_perturbation.R:285`)
  - `extract_selected_tree_cuts` (`R/grf_helpers.R:422`)
  - `find_quantile_for_proportion` (`R/find_k_inter_main.R:383`)
  - Their `man/*.Rd` pages still exist, with `\dontrun` examples.
- **Removed argument:** `showten_subgroups` was removed from `forestsearch()` and `subgroup.consistency()`. Neither has `...`, so old calls passing it now error with "unused argument".
- **Changed defaults of exported functions:**
  - `forestsearch()`:
    - `use_lasso`: TRUE → FALSE
    - `use_grf`: TRUE → FALSE
    - `dmin.grf`: 12 → 0
    - `frac.tau`: 0.6 → 0.8
    - `d0.min`: 12 → 10
    - `d1.min`: 12 → 10
    - `max_subgroups_search`: 10 → Inf
    - `vi.grf.min`: -0.2 → NULL
    - `stop_threshold`: 0.95 → `pconsistency.threshold`
    - `parallel_args`: `list(plan = "callr", workers = 6, ...)` → `list(plan = "multisession", workers = .default_parallel_workers(), ...)`
  - `subgroup.consistency()`: `stop_Kgroups`: 10 → 200.
  - `forestsearch_Kfold()`: `Kfolds`: `nrow(fs.est$df.est)` → `.fs_cv_nobs(fs.est)`.
  - `cox_ahr_cde_analysis()`: `output_dir`: `tempdir()` → `"plots"`.
- **New arguments on common exports:** 141 in total. Most come from the GLM, DINA, MR and threshold work; for example `forestsearch()` gained 34.
  - The full list was produced by parsing formals from both trees. It is reproducible with the script described in §6's method note and is not repeated here.

### Additive

- Exports: 118 → 226. **+112 added**, 4 removed.
- S3 registrations: 16 → 60. +44 added, 0 removed.
- `git diff --stat v0.1.0..HEAD -- R/ DESCRIPTION NAMESPACE vignettes/ tests/ inst/ data/`: **174 files, +62,779 / −4,713**.
  - `R/`: 107 files (65 added, 41 modified, 1 deleted: `oc_analyses_gbsg.R`), +48,519 / −3,962.
  - `NAMESPACE`: +406 / −186.
  - `DESCRIPTION`: +33 / −19.
  - `tests/`: 59 files, all added, +12,190 (0.1.0 had no tests).
  - `vignettes/`: 5 files. Added `forestsearch.Rmd` and `survival-post-selection.qmd`; deleted `forestsearch.R`; modified 2 articles.
  - `inst/`: 1 file added (`inst/patches/...md`).
  - `data/`: no change.

## 4. DESCRIPTION

- `Title: Exploratory Subgroup Identification with Survival and GLM Outcomes`
  - 66 characters (**over the 65 limit**).
  - Title case: OK (`tools::toTitleCase` is idempotent).
  - No "package" or "R" and no trailing period.
  - It no longer says "with Survival Endpoints" (that was the 0.1.0 Title). The current Title does describe the GLM extension.
- `Description:`
  - Begins "Implements advanced statistical methods…" (not "This package").
  - The DOI is in `<doi:10.1002/sim.10163>` form.
  - Single quotes are applied to `'GRF'` and `'MRCT'`, which are a method name and an acronym, not software. CRAN's quoting rule is for software and package names.
  - `LASSO`, `Cox` and `GLMs` are unquoted, which is correct for those.
  - "advanced" is promotional wording that CRAN reviewers often query.
- `Authors@R: person("Larry", "Leon", , "larry.leon.05@post.harvard.edu", role = c("aut", "cre"))`. The Maintainer is derived from it. No ORCID.
- `License: MIT + file LICENSE`. `LICENSE` reads `YEAR: 2025-2026` / `COPYRIGHT HOLDER: Larry Leon` (CRAN template form). OK.
- `URL`: GitHub repo and pkgdown site. `BugReports`: GitHub issues. `Depends: R (>= 4.1.0)`.
- `Imports` (23): data.table, doFuture, dplyr, foreach, future, future.apply, future.callr, ggplot2, glmnet (>= 5.0), grf, gt, Matrix, mvtnorm, parallel, patchwork, policytree, progressr, randomForest, rlang, stats, survival, utils, weightedsurv.
  - Since 0.1.0: `stringr` was dropped; `Matrix`, `mvtnorm`, `parallel`, `stats` and `utils` were added.
- `Suggests` (15): DiagrammeR, doRNG, htmltools, MASS, sandwich, tidyr, forestploter, cubature, svglite, knitr, rmarkdown, katex, latentcor, quarto, testthat (>= 3.0.0).
- `LinkingTo`: none. `VignetteBuilder: knitr, quarto`. **`SystemRequirements`: absent.** It is needed for the quarto vignette.
- `RoxygenNote`: **absent.** It is replaced by `Config/roxygen2/version: 8.1.0`, which matches the installed roxygen2 8.1.0. `Roxygen: list(markdown = TRUE)`: present.
- **`Remotes:`: absent** (removed before 0.1.0, in `33f652a1`, 2026-02-09). CLAUDE.md's statement that `Remotes: larry-leon/weightedsurv` is present is stale.
- **weightedsurv:**
  - No version requirement is stated (bare `weightedsurv`).
  - Imported via `importFrom(weightedsurv, df_counting, plot_weighted_km)` (roxygen tags at `R/subgroup_consistency_helpers.R:956` and `R/forestsearch_main.R:1508`).
  - `weightedsurv::` calls:
    - `plotKM.band_subgroups` (`R/plot_km_band_forestsearch.R:349`)
    - `df_counting` (`R/plot_sg_weighted_km.R:345, 354`; `R/subgroup_consistency_helpers.R:966, 967`; `R/summary_utility_functions.R:430`)
    - `plot_weighted_km` (`R/plot_sg_weighted_km.R:369, 386`; `R/subgroup_consistency_helpers.R:969, 970`)
  - Installed weightedsurv is `0.1.0` but is a **GitHub install** (RemoteSha `edd61a2c`, packaged 2026-06-16), not the CRAN binary. All 3 functions are in its namespace.
  - I also checked CRAN's own page: weightedsurv 0.1.0 was published 2025-12-22. Its reference manual documents all three: `df_counting`, `plot_weighted_km` and `plotKM.band_subgroups`. Its 46 function topics match the installed export list one-for-one.
  - Argument-level compatibility with the CRAN build was **not** verified.
  - Two error messages still tell users to `devtools::install_github('larry-leon/weightedsurv')` (`R/plot_km_band_forestsearch.R:342`, `R/plot_sg_weighted_km.R:334`), and roxygen at `R/plot_km_band_forestsearch.R:171` says the same.
- **Declared vs used:**
  - Every `Imports` entry is referenced, through `pkg::` or a `NAMESPACE` import.
  - No undeclared non-base package is used in code. The only extra hit was `devtools`, which appears inside message strings only.
  - Suggests used in code: sandwich, MASS, latentcor, forestploter, svglite, cubature. `latentcor` is guarded by `requireNamespace` (`R/mrct_inject_dgm.R:95`).
  - Suggests not referenced in `R/`: DiagrammeR, doRNG, htmltools, tidyr, plus the build and test tools. They may be used in vignettes or tests; not checked further.
  - Base-priority packages `grDevices`, `graphics`, `grid`, `splines` and `tools` are used without being listed in `Imports`. R CMD check did not flag this (dependencies and R-code sections OK).

## 5. Build contents and size

- `R CMD build --no-build-vignettes --no-manual` produced `forestsearch_0.3.5.9002.tar.gz`, **1,322,197 bytes (1.26 MB)**, 616 entries.
- Top-level contents: `DESCRIPTION`, `LICENSE`, `NAMESPACE`, `NEWS.md`, `README.md`, `R/` (113), `man/` (428), `tests/` (62), `vignettes/` (2 sources), `inst/patches/` (1), `build/partial.rdb`.
- None of these are in the tarball: `dev/`, `quarto/`, `_payloads/`, `manuscript/`, `claude/`, `.claude/`, `CLAUDE.md`, `docs/`, `pkgdown/`, `REPORT_*`, `TASK_*`, `HANDOFF_*`, `.html` files.
- Flagged contents:
  - `inst/patches/2026-05-21_extract_grf_args.md`, a dev note that installs into the user library. Candidate for `.Rbuildignore`.
  - `tests/testthat/fixtures/grf_disjunction.rds`, a test fixture. Legitimate; listed for completeness.
- Largest 20 files (bytes):
  - `R/forestsearch_main.R` 201,012
  - `R/forestsearch_helpers.R` 115,619
  - `R/forestsearch_cross_validation.R` 85,696
  - `NEWS.md` 85,161
  - `R/subgroup_consistency_helpers.R` 82,313
  - `R/fs_mr_inference.R` 81,038
  - `man/forestsearch.Rd` 69,357
  - `R/run_simulation_analysis.R` 65,389
  - `R/generate_aft_dgm_helpers.R` 65,041
  - `R/mrct_simulation.R` 64,749
  - `R/dina_subgroup.R` 60,338
  - `R/simulation_tables.R` 59,378
  - `R/subgroup_consistency_main.R` 54,867
  - `R/plot_subgroup_results_forestplot.R` 51,836
  - `R/summary_utility_functions.R` 48,900
  - `R/glm_effect_estimators.R` 48,396
  - `R/bootstrap_analysis_dofuture.R` 47,484
  - `R/dina_subgroup_bootstrap.R` 46,867
  - `R/summarize_cv_results.R` 46,719
  - `R/run_subgroup_sims.R` 45,870
- The temp tarball was deleted.

## 6. Documentation and code policy

- **Missing `\value`:** the 11 listed in §1.
- **No `\examples`:**
  - 39 exports: `analysis_covariates`, `compute_aor`, `compute_cde_glm`, `compute_frontier_cis`, `consistency_resample_compare`, `create_glm_row`, `expand_covariates`, `fs_attach_betaHhat`, `fs_betaHhat_counts`, `fs_betaHhat_neff_parity`, `fs_betaHhat_neff_report`, `fs_betaHhat_one`, `fs_betaHhat_table`, `fs_betaHhat_theta_dagger_check`, `fs_build_eval_frame`, `fs_identification_structure`, `fs_identification_table`, `fs_mr_inference`, `fs_plot_bias_coverage`, `fs_rule_covariate_pairs`, `fs_scale_se`, `fs_sim_bias_coverage`, `get_leff`, `glm_effect_profile`, `inject_mrct_structure`, `interpret_search_config`, `make_region_model`, `mr_estimates_table`, `mrct_truth_by_region`, `plot_fs_identification_involvement`, `plot_fs_identification_structure`, `plot_fs_rule_pairs`, `plot_region_modifier`, `sim_region_metrics`, `test_constant_cate`, `test_hte`, `test_hte_from_forestsearch`, `test_zero_cate`, `validate_subgroups`.
  - 22 S3 methods: `confint.dina`, `predict.dina`, `print.dina_bagged`, `print.dina_subgroup`, `print.dina_subgroup_bootstrap`, `print.dina_subgroup_refit`, `print.forestsearch_comparison`, `print.fpr_calibration`, `print.fs_cv_summary`, `print.fs_dgm_scale`, `print.fs_fdr_report`, `print.fs_mr_oc`, `print.fs_weighted_km`, `print.glm_dgm`, `print.glm_effect_profile`, `print.guohe_a3`, `print.guohe_ar`, `print.hte_test`, `print.leff_calibration`, `summary.dina`, `summary.fpr_calibration`, `summary.subgroup_sims`.
  - 16 registered S3 methods have no Rd alias. This is allowed for S3 methods; listed for completeness: `coef.dina`, `vcov.dina`, `print.dina`, `print.summary.dina`, `print.detection_curve`, `print.fpr_approximation`, `print.fs_identification`, `print.fs_oc_family`, `print.fs_oc_grid`, `print.fs_oc_invert`, `print.fs_oc_invert_list`, `print.fs_oc_predict`, `print.subgroup_sims`, `print.subgroup_sims_summary`, `print.summary.fs_oc_grid`, `summary.fs_oc_grid`.
- **`\dontrun{}`:** 145 Rd files (out of 427) and 146 roxygen `\dontrun{` lines in `R/`. 18 Rd files use `\donttest{}`.
  - Every `\dontrun` example has no runnable code outside the wrapper.
  - The 4 that do not document an export are the Rd pages of the 4 de-exported functions.
  - Rd files: `analyze_subgroup`, `bootstrap_results`, `bootstrap_ystar`, `build_classification_table`, `build_estimation_table`, `calibrate_cens_adjust`, `calibrate_glm_interaction`, `calibrate_k_inter`, `calibrate_k_treat`, `check_censoring_dgm`, `compare_detection_curves`, `compare_fs_subgroups`, `compare_multiple_survreg`, `compare_selection_rules`, `compare_subgroup_sims`, `compute_detection_probability`, `compute_dgm_cde`, `compute_node_metrics`, `consistency_resample`, `cox_ahr_cde_analysis`, `cox_cs_fit`, `cox_summary`, `create_dgm_for_mrct`, `create_forest_theme`, `create_gbsg_dgm`, `create_sample_size_table`, `create_summary_table`, `cv_compare_results`, `cv_metrics_tables`, `CV_sgs`, `cv_summary_tables`, `cv_summary_text`, `detect_variable_types`, `dina_subgroup_bootstrap`, `dina_subgroup_refit`, `estimate_propensity_scores`, `evaluate_comparison`, `evaluate_consistency_twostage`, `evaluate_subgroup_consistency`, `explain_pareto_selection`, `extract_grf_args`, `extract_selected_subgroups`, `extract_selected_tree_cuts`, `figure_note`, `filter_call_args`, `find_k_inter_for_target_hr`, `find_leaf_split`, `find_quantile_for_proportion`, `find_required_sample_size`, `fit_cox_models`, `fit_policy_trees`, `forest_height`, `forestsearch_bootstrap_dofuture`, `forestsearch_Kfold`, `forestsearch_KfoldOut`, `forestsearch_tenfold`, `forestsearch`, `format_bootstrap_table`, `format_fs_details`, `format_oc_results`, `fpr_calibration_transfer`, `fpr_calibration`, `frontier_member_flags`, `fs_declaration_calibration`, `fs_dgm_scale`, `fs_fdr_report`, `fs_identification_summary`, `FS_labels`, `fs_mr_oc_summary`, `generate_aft_dgm_flex`, and 75 more.
  - The full 145-file list can be regenerated with `grep -l '\\dontrun' man/*.Rd`.
- **Parallelism:**
  - Examples: 34 roxygen example lines in 14 `R/` files set `workers =` or a multisession, multicore or callr plan. All of them sit inside `\dontrun`.
  - Tests: every forestsearch run uses `plan = "sequential", workers = 1`, except two:
    - `test-search-reproducibility.R:50` (multisession; guarded by `skip_on_cran()` at :66)
    - `test-cv-no-subgroup-edges.R:247` (multisession, `workers = 2L`; that file is guarded by `skip_on_cran()`)
  - Shipped vignette `forestsearch.Rmd`: `n_cores <- 2` is hard-coded (:264) and used for multisession at :303, :476, :604 and :640. `detectCores()` at :265 is for display only.
  - Package-code defaults: see §1.
  - `run_subgroup_sims()` uses `future::availableCores()` (`R/run_subgroup_sims.R:821, 825`), which respects `_R_CHECK_LIMIT_CORES_`.
  - `.default_parallel_workers()` (`R/forestsearch_main.R:181`) caps workers at 2 under the check limit.
- **Side effects in package code:**
  - `setwd`: none.
  - `Sys.setenv`: none.
  - `installed.packages()`: none.
  - `:::`: none.
  - `T`/`F` literals: none (all hits were inside strings).
  - `options(`: `run_subgroup_sims.R:880` (restored via `on.exit`); `bootstrap_dofuture_main.R:347` (`warn = -1`, restored via `on.exit` :349-350 and :494).
  - `par(` without restore: `R/subgroup_consistency_helpers.R:968` only (internal `plot_subgroup()`). Every other `par(` call has an `on.exit(par(old))`: `plot_sg_weighted_km.R:365`, `plot_sg_results.R:563/693/806`, `find_k_inter_main.R:284`, `cox_ahr_cde_wrapper.R:378`, `cox_spline_fit.R:129/599`, `glm_spline_fit.R:155`, `mrct_inject_dgm.R:480`.
  - File writes:
    - `leff_calibration_io.R:70, 99` (default `"_data"`; see §1)
    - `cox_ahr_cde_wrapper.R:374` (gated by `save_plots`)
    - `render_forestplot.R:336-344` and `mrct_inject_dgm.R:478` (to a user-supplied filename only)
    - `latest_rds(output_dir = "_data")` (`R/leff_calibration_io.R:219`) reads from the working directory by default.
  - `set.seed`: 44 sites (see §1). The unconditional ones are at:
    - `sim_aft_gbsg.R:478, 604, 856, 884, 900, 908`
    - `setup_gbsg_dgm.R:130`
    - `fpr_calibration.R:301`
    - `forestsearch_cross_validation.R:369, 958`
    - `subgroup_consistency_main.R:58, 494`
    - `dina.R:433`
    - `grf_subg_harm_glm.R:460`
    - `generate_glm_dgm.R:402, 418`
    - `dina_bagged.R:436`
    - `run_simulation_analysis.R:514`
    - `run_subgroup_sims.R:512`
    - `mrct_inject_dgm.R:79, 392, 448`
    - `generate_aft_dgm_main.R:248`
    - `bootstrap_analysis_dofuture.R:334`
    - `get_FSdata_helpers.R:259`
    - `synthetic_data_perturbation.R:98, 506`
    - `truefind_asymptotic_glm.R:338`
  - `cat(`: 1,115 non-comment lines in `R/`; `print(`: 46. I did not classify which of these are inside print methods; the volume is noted only.
- **Non-ASCII in `R/`:** 169 lines in 27 files, **all in comments**: 162 full-comment lines and 7 trailing comments (`simulation_tables.R:700-701`, `run_simulation_analysis.R:1031`, `summary_utility_functions.R:779-784`). The check did not flag them. `R/fs_bias_coverage.R`, the 2026-09-18 offender, now has 0.
- **`globalVariables`:** the 2026-09-18 NOTE variables (9 in `fs_plot_bias_coverage`) are covered (`R/globals.R:264`). `checking R code for possible problems ... OK`.
- **Roxygen:**
  - `Roxygen: list(markdown = TRUE)` is present.
  - There is no `RoxygenNote`; `Config/roxygen2/version: 8.1.0` matches the installed version.
  - **`man/` and `NAMESPACE` are not stale.** I ran `roxygen2::roxygenise()` on a temp copy, both incrementally and as a full regeneration from an empty `man/`. All 427 Rd files and `NAMESPACE` came out byte-identical to the repo.
  - The full regeneration printed "Could not resolve link" messages. These are artifacts of regenerating from an empty `man/`; the check's Rd cross-reference section is OK.
  - The temp copy was deleted.
- Method note: exports and formals were compared by parsing `NAMESPACE` and top-level `name <- function(...)` definitions from `git archive v0.1.0` against the working tree. Nothing was installed.

## 7. Tests and vignettes

- Tests:
  - 55 `test-*.R` files and 57 `.R` files in `tests/testthat/` (including helpers).
  - 459 `test_that(` blocks.
  - `skip_on_cran()`: 41 calls in 13 files.
- Under check: **`[ FAIL 0 | WARN 27 | SKIP 77 | PASS 5452 ]`, elapsed 132.7 s** (user 144.5 s).
  - Skips: 58 "On CRAN".
  - 10 are probe skips: "resolver evaluated to an error: could not find function `.probe_record`" in `test-directive-c.R:662`, `test-threshold-pair-directive-a.R:45/66/87/127` and `test-threshold-sync.R:106/123/141/152/165`. These files now skip under check instead of failing as they did on 2026-09-18.
  - Also: baseline CSV not in tree (1), campaign template not in tree (1), package source not available (1), and others.
- Tests that run heavy entry points without `skip_on_cran()`:
  - `test-mr-inference.R` runs `forestsearch_bootstrap_dofuture` (:367, :394, :483, :490; `nb_boots = 2L` or `nb`), `forestsearch_Kfold` (:380, `Kfolds = 2L`) and `forestsearch_tenfold` (:514, :529, :563, :569).
  - These all run under a sequential plan (`.mri_seq`, :347).
  - 21 other files without `skip_on_cran()` call `forestsearch()`, also sequentially.
  - Total test time under check is 133 s.
- Vignettes shipped:
  - `vignettes/forestsearch.Rmd`: `knitr::rmarkdown`, 41 KB. Evaluated live, not precomputed; bootstrap with `NB <- 2` and CV with `Ksims <- 1`. `oob-cv` chunk has `eval = FALSE`.
  - `vignettes/survival-post-selection.qmd`: `quarto::html`, 12 KB, `parallel_args = list(plan = "sequential")`.
  - `vignettes/articles/` (552 KB) is `.Rbuildignore`d.
- `inst/doc`: not present in the source tree.
- Installed size with vignettes ignored: **3.8 MB** (`R` 2.0 MB, `help` 1.6 MB, `NEWS.md` 84 KB). There is no `doc/`, so this is not comparable to CRAN 0.1.0's 6.1 MB / doc 4.6 MB. The vignette build was skipped (P3).
- `NEWS.md`: present, 1,533 lines. Last heading is `# forestsearch (development version)`; the previous one is `# forestsearch 0.3.5`. There is no `0.1.0`-relative entry for the CRAN release.
- `cran-comments.md`: present (`.Rbuildignore`d).
- Vignette build elapsed: **not run** (P3).

## 8. CRAN check

### Commands

- Build: `R CMD build --no-build-vignettes .`, run in a scratch directory. The P3 fallback applied because pandoc and quarto are not on PATH.
- Run 1 (12:38):
  - `_R_CHECK_CRAN_INCOMING_REMOTE_=FALSE _R_CHECK_LIMIT_CORES_=TRUE R CMD check --as-cran --ignore-vignettes`
  - Stopped at dependencies.
- Run 2 (12:39–12:43, 265 s):
  - Same command plus `_R_CHECK_FORCE_SUGGESTS_=false`.
  - Its `00check.log` is copied to `dev/reports/cran_audit_00check_2026-10-09.log`.
- Hard timeout: 4,500 s via `perl alarm` (macOS has no `timeout`). Not reached.

### Run 1 status: `Status: 1 ERROR, 1 NOTE`

- `checking package dependencies ... ERROR`
  ```
  Package suggested but not available: 'latentcor'

  The suggested packages are required for a complete check.
  Checking can be attempted without them by setting the environment
  variable _R_CHECK_FORCE_SUGGESTS_ to a false value.
  ```

### Run 2 status: `Status: 1 ERROR, 2 WARNINGs, 4 NOTEs`

- `checking CRAN incoming feasibility ... NOTE`
  ```
  Maintainer: 'Larry Leon <larry.leon.05@post.harvard.edu>'

  Version contains large components (0.3.5.9002)

  Package has a VignetteBuilder field but no prebuilt vignette index.
  ```
- `checking package dependencies ... INFO`: `Package suggested but not available for checking: 'latentcor'`
- `checking top-level files ... NOTE`: `Files 'README.md' or 'NEWS.md' cannot be checked without 'pandoc' being installed.`
- `checking Rd \usage sections ... WARNING`
  ```
  Undocumented arguments in Rd file 'analysis_covariates.Rd'
    'dgm' 'case' 'drop'

  Undocumented arguments in Rd file 'mrct_truth_by_region.Rd'
    'dgm' 'seed'

  Undocumented arguments in Rd file 'plot_region_modifier.Rd'
    'dgm' 'file' 'main'

  Undocumented arguments in Rd file 'sim_region_metrics.Rd'
    'sim' 'region_col' 'subgroup_index' 'pi_preserve'

  Functions with \usage entries need to have the appropriate \alias
  entries, and all their arguments documented.
  The \usage entries must correspond to syntactically valid R code.
  See chapter 'Writing R documentation files' in the 'Writing R
  Extensions' manual.
  ```
- `checking files in 'vignettes' ... SKIPPED`
- `checking PDF version of manual ... WARNING`
  ```
  LaTeX errors when creating PDF version.
  This typically indicates Rd problems.
  ```
- `checking PDF version of manual without index ... ERROR`
  - The log shows only `Re-running with no redirection of stdout/stderr.`
  - The console shows `Converting parsed Rd's to LaTeX ......` then `Error in texi2dvi(...) : pdflatex is not available`.
  - So this is the environment, not an Rd defect, but it also means the Rd files are not certified.
- `checking HTML version of manual ... NOTE`
  ```
  Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.
  ```
- `checking for non-standard things in the check directory ... NOTE`: `'forestsearch-manual.tex'` (left over from the failed PDF step).
- Everything else is OK. That includes:
  - R code for possible problems (19 s)
  - Rd cross-references
  - Examples (14 s)
  - Examples with `--run-donttest` (100 s user / 51 s elapsed)
  - Tests (152 s / 133 s)

### Comparison to 2026-09-18 (`0 errors | 1 warning | 2 notes`; certification run `1 error | 2 warnings | 2 notes`)

- Non-ASCII WARNING in `R/fs_bias_coverage.R`: **gone.**
- 9 "no visible binding" NOTE in `fs_plot_bias_coverage`: **gone.**
- HTML Tidy NOTE: **same.**
- Certification-run ERROR (`FAIL 14` probe tests): **gone.** Now FAIL 0; those cells skip.
- Certification-run WARNING (Rd cross-reference `subgroup_search` in `evaluate_combination_with_status.Rd`): **gone.**
- **New, package-caused:** the Rd `\usage` WARNING (`R/mrct_inject_dgm.R`, added 2026-10-02).
- **New, environment-caused on this machine:** the PDF manual ERROR and WARNING (no pdflatex on PATH), the README/NEWS NOTE (no pandoc), the `manual.tex` leftover NOTE, and the latentcor dependency ERROR/INFO.
  - The 2026-09-18 runs had LaTeX available ("checking PDF version of manual ... OK").
- The incoming NOTE (version components) was not recorded in the 2026-09-18 summary.

### Examples timing (`forestsearch-Ex.timings`)

- Examples over 5 s elapsed:
  - `run_null_calibration`: 74.828 user / 27.296 elapsed. **CPU/elapsed 2.74.**
  - `suggest_thresholds`: 5.635 user / 5.987 elapsed.
- The timings file covers all 208 timed topics, totalling 83.6 s user / 37.2 s elapsed. The run producing it included `\donttest` code.
- Tests total elapsed: 132.7 s.
- Vignette build elapsed: not run (P3).

### Other tools

- `urlchecker::url_check()`: installed, but **could not run**: `pandoc is not installed and on the PATH.`
- `spelling::spell_check_package()`: ran.
  - 531 distinct words flagged; no `inst/WORDLIST` and no `Language:` field (defaults to en-US).
  - The top hits are domain terms and abbreviations: ForestSearch, GRF, DGM, IRR, et al, Hc, ITT, AHR, GBSG, IJ, MRCT, estimand, RMST, IPTW, PPV/NPV.
  - Also British spellings: behaviour, randomisation.
  - No triage of real typos was done.
- Temp build and check directories were deleted after extraction.

## 9. Open questions for Larry

1. **TeX for the certification run.** The PDF manual cannot be certified on this machine because `pdflatex` is not on PATH, even though `~/Library/TinyTeX` exists. Should Stage 1's verification put TinyTeX (and the RStudio pandoc/quarto) on `PATH` for the check? The 2026-09-18 certification run evidently had all three.
2. **Which tree is "0.1.0"?** The `v0.1.0` tag (`29890c7c`, on `cran-review-1`) is not an ancestor of `feature/glm-extension`, and it contains no `vignettes/forestsearch.Rmd`. Yet CRAN 0.1.0 shipped a 4.6 MB `doc/`. Is the tag tree exactly the uploaded tarball, or should Stage 1 diff against the CRAN tarball itself?
3. **De-exports and the removed `showten_subgroups`:** are these intentional API breaks to document in NEWS, or should Stage 1 restore them, for example as deprecated shims?
4. **`run_null_calibration` CPU ratio (2.74):** should Stage 1 diagnose the threading source, or simply shrink or wrap the example? The cause is not verified from source.
