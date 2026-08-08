# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

`forestsearch` is an R package (CRAN-targeted, lifecycle: experimental) implementing the Forest Search methodology from Leon et al. (2024, *Statistics in Medicine*, doi:10.1002/sim.10163). It performs exploratory post-hoc subgroup identification in clinical trials: it enumerates candidate subgroups from baseline factors, screens them against an effect threshold indicative of harm/benefit, selects the subgroup "maximally consistent with harm" via repeated split-sample evaluation, and produces bias-corrected treatment-effect estimates with infinitesimal-jackknife confidence intervals.

Originally survival-only (Cox HR), the package is being extended on the current branch (`feature/glm-extension`) to GLM outcomes: binary (OR/RR/RD), continuous (mean difference), and count/rate (IRR via Poisson/quasi-Poisson/negative-binomial). Much code branches on `outcome_type` (`"survival"`/`"binary"`/`"continuous"`/`"count"`) and `effect_measure`.

## Development commands

This is a standard R package built with `devtools`/`roxygen2`. Run from the package root (an R session, e.g. RStudio, or `Rscript`):

```r
devtools::load_all()          # load all R/ functions for interactive dev (use this, not library())
devtools::document()          # regenerate man/*.Rd + NAMESPACE from roxygen (run after editing @importFrom/@export or roxygen)
devtools::test()              # run the full testthat suite
devtools::test(filter = "glm-pipeline")   # run one test file (tests/testthat/test-glm-pipeline-integration.R)
devtools::check()             # R CMD check
devtools::check(args = c("--as-cran", "--no-build-vignettes"))  # CRAN-style check, faster
devtools::install(dependencies = FALSE)   # install locally
```

Notes:
- After editing any roxygen block, `@importFrom`, or `@export`, run `devtools::document()` — `NAMESPACE` and `man/*.Rd` are generated, never edit by hand.
- `dev/devtools_package.R` is a scratch log of every dev/CRAN/pkgdown command the author has used (not a script to run top-to-bottom); consult it for the exact incantation of a rare task (e.g. corrupt-install recovery, pkgdown deploy).
- Vignettes/articles are Quarto (`.qmd`) under `quarto/` and `vignettes/`; render with `quarto render <file>.qmd`. Preserve `embed-resources: true`. The package `vignettes/forestsearch.Rmd` is the knitr-built one.
- If package install/load errors with a corrupt `.rdb`, remove the installed copy (`rm -rf ~/R/x86_64-pc-linux-gnu-library/<ver>/forestsearch`) and restart R before reinstalling.

## Architecture

The pipeline is a four-stage flow; each stage is a separate exported entry point that consumes the previous stage's result object.

**1. Search — `forestsearch()`** (`R/forestsearch_main.R`, the largest file at ~2700 lines; helpers in `forestsearch_helpers.R`).
Constructs candidate binary factors (LASSO via `glmnet`, GRF via `grf`, quartile splits), enumerates one- and two-factor subgroup combinations meeting min-sample/event criteria, screens against the effect threshold, and runs split-sample consistency to select the estimated subgroup. Returns a `forestsearch` S3 object. Data preparation (factor processing, dummy encoding, candidate construction) lives in `get_fsdata.R` + `get_FSdata_helpers.R`. Subgroup enumeration/scoring is in `subgroup_search.R`.

**2. Consistency — `subgroup.consistency()`** (`subgroup_consistency_main.R` + `subgroup_consistency_helpers.R`).
Repeated random 50/50 splits (default R=400), optionally with two-stage sequential early stopping. `select_best_subgroup()` applies the selection rule (`"hr"`, `"maxSG"`, `"minSG"`, `"hrMaxSG"`, `"hrMinSG"`).

**3. Bootstrap bias correction — `forestsearch_bootstrap_dofuture()`** (`bootstrap_dofuture_main.R`, `bootstrap_analysis_dofuture.R`, `bootstrap_calculations_helpers.R`, summaries in `bootstrap_summaries_helpers.R` / `summarize_bootstrap_results.R`).
Re-runs the full search inside each bootstrap replicate, applies two-source bias correction on the log-effect scale, and computes IJ variance. Parallelised via `doFuture`/`future`/`foreach`.

**4. Cross-validation — `forestsearch_Kfold()` / `forestsearch_tenfold()`** (`forestsearch_cross_validation.R`, `summarize_cv_results.R`).
N-fold (leave-one-out) and repeated K-fold CV for stability; `cv_metrics_tables()` reports sensCV/ppvCV/exact-match.

**Effect estimation backends.** Cox-based estimation: `Cox_estimation_helpers.R`, `cox_ahr_cde_wrapper.R`, `cox_spline_fit.R`. GLM-based estimation (the extension): `glm_effect_estimators.R`, `glm_effect_profile.R`, `glm_spline_fit.R`, `generate_glm_dgm.R`, `grf_subg_harm_glm.R`, `calibrate_glm_interaction.R`. `make_effect_estimator()` dispatches on outcome type.

**GRF layer.** `grf_main.R`, `grf_helpers.R`, `grf_args.R`, `grf_subgroup_labels.R` wrap Generalized Random Forests for candidate/policy-tree subgroup discovery.

**DINA.** The `dina*.R` family (`dina.R`, `dina_subgroup.R`, `dina_bagged.R`, `dina_subgroup_bootstrap.R`, `dina_subgroup_refit.R`) is a distinct "difference-in-nested-... " bagged-subgroup method with its own S3 classes and print/summary/plot methods.

**Simulation & operating characteristics.** `generate_aft_dgm_*.R` (flexible AFT data-generating mechanism), `setup_gbsg_dgm.R`/`sim_aft_gbsg.R` (GBSG-based DGM for reproducing paper sims), `simulate_from_dgm.R`, `run_simulation_analysis.R`, `oc_analyses.R`, `simulation_tables.R`, `truefind_asymptotic*.R`. MRCT regional-consistency sims in `mrct_simulation.R`.

**Visualization & tables.** `plot_subgroup_results_forestplot.R` (main forest plot), `plot_sg_weighted_km.R`, `plot_km_band_forestsearch.R`, `plot_spline_treatment_effect.R`, `gg_forest.R`, `render_forestplot.R`; `gt`-based tables in `create_summary_table.R`, `format_subgroup_summary_tables.R`, `summary_utility_functions.R`.

**S3 dispatch** is the primary polymorphism mechanism — `print`/`summary`/`plot`/`coef`/`confint`/`predict` methods across `forestsearch`, `dina*`, `fs_forestplot`, `fs_kfold`, `fs_tenfold`, `fpr_calibration`, etc. `forestsearch_methods.R` and `fs_debias_gate_methods.R` hold method definitions. Check `NAMESPACE` (generated) for the full S3 method registry.

## Conventions specific to this repo

- **File organization:** one major exported function per `R/<name>.R` file, named after the function. Single-use helpers live in the same file, prefixed with `.` and marked `@noRd`. Helpers shared across functions get their own `*_helpers.R` file.
- **Global variables:** all `utils::globalVariables()` declarations are consolidated in `R/globals.R` (to suppress NSE NOTEs from data.table/ggplot2/dplyr). Do not scatter `globalVariables()` calls elsewhere — add new NSE column names to `globals.R`.
- **Dependency `weightedsurv`** is a GitHub-only package (`Remotes: larry-leon/weightedsurv`), installed automatically from source, not CRAN.
- **Parallelism:** entry points take a `parallel_args = list(plan = "multisession", workers = N)`. Default worker count is `floor(0.80 * physical_cores)`, auto-capped to 2 under `R CMD check` (`_R_CHECK_LIMIT_CORES_`) for CRAN compliance — see `.default_parallel_workers()`.
- **Argument capture:** `forestsearch()` mirrors resolved formals into `args_call_all` via `.sync_args_call_all()` so bootstrap replicates can reconstruct the exact call; when adding a formal that must survive into bootstrap, add it to that sync list.
- **DGM wrappers** around `generate_aft_dgm_flex()` take a `base_args` named list (not `...`) and override with `utils::modifyList(base_args, list(...))`, never `c()`.
- **CRAN hygiene:** use `pkg::fn()` with a matching `@importFrom`; ASCII only in R source (use `>=` not `≥`, or `≥` in strings); wrap slow examples in `\dontrun{}`/`\donttest{}`.

## Testing conventions

- Tests live in `tests/testthat/test-<name>.R` (edition 3). Current suite focuses on cross-cutting contracts: input validation, return-shape contracts, cross-outcome parity (survival vs GLM giving consistent shapes), get_fsdata diagnostics, GLM pipeline integration, and downstream consumers. `helper-synthetic-dgm.R` provides small synthetic data.
- When a test fails, confirm the fault is in the test before editing it. Do not modify tests to pass without an actionable reason — the author has explicitly flagged this as a recurring problem.

## Working style (author preferences)

These reflect standing instructions the author keeps in `quarto/claude code notes/CLAUDE.md`:

- **Scope discipline:** modify only what is explicitly requested. Do not opportunistically fix, refactor, or restyle unrelated code. Flag side issues separately at the end of a turn rather than silently fixing them.
- **Source of truth:** the file at the path cited, or the most recent version shared in conversation. After editing a file, re-read it before further edits rather than trusting a stale view.
- **Don't commit or push** without explicit approval; propose commit messages instead. Never rewrite history (`rebase`/`amend`/`reset --hard`) without asking. Don't bump `DESCRIPTION` version/deps, run `install.packages()` for nontrivial chains, or delete files without asking first.
- Present tradeoffs as a short numbered list with a clear recommendation. Honest disagreement is welcome; explain reasoning for substantive statistical claims rather than asserting.
- **Statistical vocabulary** (used precisely throughout vignettes/code): distinguish *marginal Cox HR* (`dgm$hazard_ratios$overall`, attenuated by non-collapsibility — never call it the "true HR"), *average hazard ratio (AHR)*, and *individual/patient-level HR*. Calibration targets the AHR by default.
- **Certification surface** is bare `R CMD check --as-cran` (via `rcmdcheck::rcmdcheck(args = "--as-cran")`) — it builds the PDF manual, so it catches LaTeX-unsafe Rd content that other surfaces miss.
- **`devtools::check()` is the dev loop, not certification:** its `manual = FALSE` default passes `--no-manual`, skipping the PDF manual build, so it can report `0 errors | 0 warnings | 0 notes` on a tree that the certification surface flags.
