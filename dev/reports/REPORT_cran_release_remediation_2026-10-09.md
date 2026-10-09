# REPORT — CRAN release remediation, Stage 1 (2026-10-09)

- Task: `dev/tasks/TASK_cran_release_remediation_2026-10-09.md` (committed alone at `dc3ef230`).
- Branch `feature/glm-extension`. No branch created or switched. No fetch, pull, push, tag or upload.

## 1. Outcome

- **Final check status: `Status: 1 NOTE`** (0 ERRORs, 0 WARNINGs). The NOTE is the HTML Tidy environment NOTE.
- P1 / P2 / P3: **PASS**.
- Gate 2 (documentation): **PASS**.
- Gate 3 (DESCRIPTION): **PASS**.
- Gate 4 (weightedsurv messages): **PASS**. One extra site was edited; see §4 and PC1.
- Step 5 (`.Rbuildignore`): done.
- Step 6 (`save_leff_calibration()`): **PASS** (new test green).
- Gate 7a (diagnosis): **PASS**. grf is the source, so 7b–7c ran.
- Gate 7 (thread cap): **PASS**. Ratio 1.43; new test green.
- Step 8 (NEWS, cran-comments): done.
- Gate 9 (release check): **PASS**.
- Post-conditions: PC1–PC6 **PASS**. PC1 passes under Larry's amendment, which also allows `R/summary_utility_functions.R` (§9). Detail in §9.
- Disposition (Larry, 2026-10-09): deviations 1–6 reported at the end of the run are **accepted**. They are: the extra Step 4 site; base R in place of `withr` in the test; latentcor dependencies installed as binaries into the temp library; the `save_leff_calibration` doc and example edits; the NEWS layout choices; and the separate `forestsearch-package.Rd` commit.

## 2. Repo state

- P1: `git status --porcelain` before the task-doc commit returned *(empty; clean)*.
- P2: HEAD at start was `00ec78a4`, a descendant of `541fef30` (`merge-base --is-ancestor` exit 0). The task-doc commit is `dc3ef230ef3d1b7c409043fde52f4af41824e8b1`.
- P3: all three tools are on `PATH`:
  - `pandoc 3.8.3`
  - `1.9.38` (quarto)
  - `pdfTeX 3.141592653-2.6-1.40.29 (TeX Live 2026)`
- Environment: R 4.5.2 (aarch64-apple-darwin20), macOS 27.0.1.
- User-library weightedsurv at start: version 0.1.0, `RemoteSha` `edd61a2cecf4b08352719af0a24b10fa1e0b041c` (a GitHub install).
- HEAD at end: the report commit. Commits in order:
  - `dc3ef230` task doc
  - `19f3dc2b` docs (Step 2)
  - `48674761` DESCRIPTION, messages, `.Rbuildignore` (Steps 3–5)
  - `e282d914` `forestsearch-package.Rd` regenerated from the new DESCRIPTION
  - `77e35a45` `save_leff_calibration` (Step 6)
  - `fdd3299e` thread cap (Step 7)
  - `e48a17f1` NEWS
  - `a1951f8f` cran-comments
  - report and check log

## 3. Step 1 — CRAN 0.1.0 tarball vs the `v0.1.0` tag

- Source: `https://cran.r-project.org/src/contrib/forestsearch_0.1.0.tar.gz`. The main URL worked, so the Archive fallback was not needed.
  - `Packaged: 2026-03-19 04:27:29 UTC`; `Date/Publication: 2026-03-23 17:20:14 UTC`.
  - It ships `inst/doc/forestsearch.{Rmd,R,html}` and `vignettes/forestsearch.Rmd`.
- Method: same as the audit. Parse `NAMESPACE` and the top-level `name <- function(...)` formals from the extracted tarball, then compare them with the working tree.
- **Exports: the tarball has 45; HEAD has 226.** 181 were added and **0 removed**. S3 registrations went from 8 to 60, with 0 removed.
- **This does not match the audit's tag-based §3** (118 exports at the tag, 4 removed).
  - The audit's 4 "removed exports" were never exported in the CRAN 0.1.0 tarball: `create_sample_size_table`, `detect_variable_types`, `extract_selected_tree_cuts` and `find_quantile_for_proportion`.
  - So they are **not** API breaks relative to CRAN, and NEWS does not list them.
  - The `v0.1.0` tag tree (`29890c7c`) therefore differs from the uploaded tarball, at least in `NAMESPACE`.
- **Removed argument:** `showten_subgroups` from `forestsearch()` and `subgroup.consistency()`. This matches the audit.
- **Changed defaults:** identical to the audit's list.
  - `forestsearch()`:
    - `use_lasso` TRUE→FALSE
    - `use_grf` TRUE→FALSE
    - `dmin.grf` 12→0
    - `frac.tau` 0.6→0.8
    - `d0.min` 12→10
    - `d1.min` 12→10
    - `max_subgroups_search` 10→Inf
    - `vi.grf.min` -0.2→NULL
    - `stop_threshold` 0.95→`pconsistency.threshold`
    - `parallel_args` callr/6 → multisession/`.default_parallel_workers()`
  - `subgroup.consistency()`: `stop_Kgroups` 10→200.
  - `forestsearch_Kfold()`: `Kfolds` `nrow(fs.est$df.est)`→`.fs_cv_nobs(fs.est)`.
  - `cox_ahr_cde_analysis()`: `output_dir` `tempdir()`→`"plots"`.
- The tarball was extracted in the scratchpad (read-only) and deleted at the end.

## 4. Steps 2–6

### Step 2 — documentation

- All cited sites were re-read and were at the audit's lines: `R/mrct_inject_dgm.R:425/443/474/502`, `R/grf_helpers.R:625`, `R/plot_sg_results.R:871/914`, `R/compare_selection_rules.R:811`, `R/forestsearch_cross_validation.R:1984/2011`, `R/sim_aft_gbsg.R:1155`.
- `@param` was added to:
  - `analysis_covariates` (`dgm`, `case`, `drop`)
  - `mrct_truth_by_region` (`dgm`, `seed`)
  - `plot_region_modifier` (`dgm`, `file`, `main`)
  - `sim_region_metrics` (`sim`, `region_col`, `subgroup_index`, `pi_preserve`)
  - Each description was written from the function body.
- `@return` was added to all 11 topics. Each was verified from the body:
  - The 4 MRCT functions return, respectively: a character vector; a data frame with attribute `smd_x_pred`; `NULL` invisibly; a list.
  - `print_grf_details`: incidental, invisible — `tree_cuts$all` or `NULL`.
  - `print.fs_sg_plot`, `print.forestsearch_comparison`, `print.fs_kfold`, `print.fs_tenfold`, `print.gbsg_dgm`: `x`, invisibly.
  - `plot.fs_sg_plot`: `x`, invisibly. The body draws nothing; it only messages which `plot_sg_results()` call to make, and errors on an unknown character `which`. The `@return` says so.
- Gate 2:
  - (a) `NAMESPACE` is byte-identical to `541fef30`.
  - (b) Exactly the 11 Rd files changed.
  - (c) In the 6 R/ files, every added or removed line matches `^[+-]\s*#'`.

### Step 3 — DESCRIPTION

- `Version: 0.4.0`.
- `Title: Exploratory Subgroup Identification in Survival and GLM Outcomes`: 64 characters, and `toTitleCase` is idempotent.
- `Description`: dropped "advanced"; `('GRF')`→`(GRF)`; `('MRCT')`→`(MRCT)`. No other wording changed.
- Added `SystemRequirements: Quarto command line tools (https://github.com/quarto-dev/quarto-cli)`.
- `weightedsurv` → `weightedsurv (>= 0.1.0)`.
- Gate 3 was confirmed via `read.dcf`.
- Consequence: `man/forestsearch-package.Rd` (Title and Description) was regenerated by roxygen and committed separately (`e282d914`).

### Step 4 — weightedsurv messages

- All advice now says `install.packages("weightedsurv")`:
  - `R/plot_km_band_forestsearch.R:171` (roxygen, which says "installed from CRAN"; regenerates `man/plot_km_band_forestsearch.Rd`)
  - `R/plot_km_band_forestsearch.R:342`
  - `R/plot_sg_weighted_km.R:334`
- **Extra site, not cited by the audit:** `R/summary_utility_functions.R:428` (`rmst_calculation()`) also said `install_github('larry-leon/weightedsurv')`.
  - Gate 4 cannot pass without editing it, so its message string was changed the same way.
  - This file was not in PC1's original allowed list. **PC1 amended by Larry** to allow it (Step 4, messaging-only).
- Gate 4: `grep -rn "install_github" R/ man/` returns nothing.

### Step 5 — build exclusions

- `^inst/patches$` was appended to `.Rbuildignore`. The file was not moved.

### Step 6 — `save_leff_calibration()`

- `R/leff_calibration_io.R:68-73`: when `output_dir` is `NULL`, the function now stops. The message says `output_dir` must be supplied and names the former default `"_data"`.
  - `latest_rds()` and `load_leff_calibration()` are untouched.
- Documentation edits beyond the code line:
  - The `@param output_dir` text no longer claims a `"_data"` default.
  - The `\dontrun` example now passes `output_dir = tempdir()`; without that it would error.
- Callers of `save_leff_calibration(` that omit `output_dir` (searched `R/`, `tests/`, `vignettes/`, `quarto/`, `dev/`, `inst/`):
  - `quarto/guides/leff_workflow_guide.md:52` (a code block in a guide). **Not edited.**
  - `quarto/legacy calibration files/leff_calibration_io.R:36` is the roxygen example of a legacy copy of the function itself, not a call to the package function. Listed for completeness; not edited.
- All other code calls pass `output_dir`:
  - `quarto/calibration/calibration_survival_fourway.qmd:1173`
  - `calibration_binary_leff.qmd:1121`
  - `selection_survival_frontier.qmd:1508`
  - `calibration_survival_leff_grid.qmd:743`
  - the 5 calls in `quarto/legacy calibration files/`
- No callers under `tests/`, `vignettes/`, `inst/` or `R/`.
- New test `tests/testthat/test-save-leff-calibration.R`:
  - Calling without `output_dir` errors, and nothing is created in the working directory.
  - With `output_dir = tempfile()`, the bundle is written there.
  - Result: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 6 ]`.

## 5. Step 7 — `run_null_calibration` CPU ratio

- The example was extracted from `man/run_null_calibration.Rd` with `tools::Rd2ex()`, with `\donttest` code included.
- It was timed with `system.time()` around `source()` after `devtools::load_all()`, under `_R_CHECK_LIMIT_CORES_=TRUE OMP_THREAD_LIMIT=2`.

| Run | user | system | elapsed | (user+system)/elapsed |
|---|---|---|---|---|
| 7a as written (`use_grf = TRUE`) | 74.724 | 2.718 | 27.073 | **2.86** |
| 7a `use_grf = FALSE` | 8.503 | 0.478 | 9.110 | **0.99** |
| 7c as written, after the cap | 55.864 | 1.634 | 40.083 | **1.43** |
| Step 9 check (`forestsearch-Ex.timings`) | 56.348 | 1.671 | 40.298 | **1.44** |

- Gate 7a: 0.99 ≤ 2.0 and 2.86 > 2.5, so grf is the source.
  - grf trains with its own `std::thread` pool, which neither the future plan nor `OMP_THREAD_LIMIT` controls.
- 7b: new file `R/utils_check_threads.R` with internal `.fs_grf_num_threads(user = NULL)`, as specified.
- None of the existing grf calls forwarded a `num.threads`, so every site now passes `num.threads = .fs_grf_num_threads()`. Without `_R_CHECK_LIMIT_CORES_` this is `NULL`, which is grf's own default, so arguments are identical.
- Sites changed (9):
  - `R/forestsearch_main.R:3256`, `:3263` (`causal_survival_forest`)
  - `R/grf_helpers.R:74`, `:84` (`causal_survival_forest`, in `fit_causal_forest()`)
  - `R/grf_helpers.R:123`, `:132` (`causal_forest`, in `fit_causal_forest_glm()`)
  - `R/grf_subg_harm_glm.R:463`, `:467` (`causal_forest`)
  - `R/glm_effect_estimators.R:1354` (`regression_forest`)
- No other `grf::*forest` training call exists in `R/`.
- Not capped (outside the spec, recorded only): `predict()` and `grf::variable_importance()` on fitted forests. The measured ratio shows these do not push the example over 2.0.
- The cap costs elapsed time: the example went from 27.1 s to 40.1 s under check settings.
- 7c test `tests/testthat/test-check-thread-cap.R` checks:
  - `NULL` when the variable is unset
  - `NULL` for `""`, `"false"` and `"FALSE"`
  - `2L` for `"TRUE"`
  - the user value when supplied
  - Result: `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 7 ]`.
  - **Deviation:** the test uses base `Sys.setenv()` / `Sys.unsetenv()` with `on.exit()` restore, not `withr::local_envvar`. `withr` is not declared in `DESCRIPTION`, and using it in tests would trigger "unstated dependencies in tests".
- Gate 7: ratio 1.43 ≤ 2.0 and the test passes.

## 6. Step 8 — NEWS and cran-comments

- `NEWS.md`:
  - `# forestsearch (development version)` → `# forestsearch 0.4.0`.
  - Directly under it: `## Changes since the last CRAN release (0.1.0)`.
  - `### Breaking changes` comes from the tarball-derived list (§3), so it has **no de-export bullets**: none of the 4 was exported on CRAN. It covers:
    - `showten_subgroups` (with its replacement `show_candidate_summary`, which is verified to be a formal of both functions)
    - every changed default, old → new
    - `save_leff_calibration()` now requiring `output_dir`
  - `### New features`, one bullet per area, each supported by existing NEWS text or source:
    - GLM outcomes (OR / MD / IRR)
    - DINA and GRF identifiers (`subgroup_method`, `dina_subgroup()`, `grf.subg.harm.glm()`)
    - MR with IJ SE and field bounds (`fs_mr_inference()`)
    - `fs_oc_family_enumerate()`, `fs_oc_predict()`, `fs_oc_grid()` and `fs_oc_invert()`
    - the c2 = 0.8·c1 rule and the c2 > c1 error (survival HR and binary OR only)
    - `fs_declaration_calibration()`
    - `pconsistency.digits`
    - the Step 7 thread cap
  - Placement deviation: the Step 6 change is listed under Breaking changes rather than New features, because it removes a default write.
  - The 14 headings `# forestsearch 0.2.0` … `0.3.5` became `## forestsearch <v> (development; not released on CRAN)`. Their bodies are unchanged. `# forestsearch 0.1.0` stays a top-level heading.
  - Added heading (deviation): `## Development changes after 0.3.5 (not released on CRAN)` above the former development-version body. Without it, those bullets would read as part of `### New features`.
  - Known cosmetic effect: after demotion, the intermediate sections' own `##` subheadings sit at the same level as the version headings. They were left unchanged because the task requires unchanged bodies.
- `cran-comments.md` was rewritten:
  - Submission: 0.1.0 → 0.4.0; intermediate versions were development-only; API changes are in NEWS.
  - Test environments: the local check, plus win-builder and mac-builder TBD.
  - Check results: verbatim, with a one-line explanation of the NOTE.
  - Reverse dependencies: `tools::package_dependencies("forestsearch", db = available.packages(repos = "https://cloud.r-project.org"), reverse = TRUE, which = "all")` returned `$forestsearch: chr(0)`, i.e. none.

## 7. Step 9 — release check

- Temp library `<scratch>/rlib`:
  - `weightedsurv` 0.1.0, `Repository: CRAN`, built from CRAN source.
  - `latentcor` 2.0.2, `Repository: CRAN`, built from CRAN source.
  - Confirmed: under `R_LIBS=<scratch>/rlib`, `find.package("weightedsurv")` resolves to the temp library (`Repository: CRAN`).
  - **Deviation:** latentcor's dependencies failed to compile from source because there is no gfortran on this machine (`gss`, `qap`, `pcaPP`, `mnormt`, …). The 9 missing dependencies (`pcaPP fMultivar mnormt heatmaply fBasics sn seriation gss qap`, plus their own dependencies) were installed as **CRAN macOS binaries into the same temp library**. Nothing went into the user library.
- Build: `R CMD build <repo>` with vignettes, in a scratch directory, `R_LIBS=<scratch>/rlib`.
  - Exit 0, 50 s total (including the vignette build; the build log does not time vignettes separately).
  - Tarball `forestsearch_0.4.0.tar.gz`: **3,346,014 bytes**.
- Check: `R CMD check --as-cran` with `R_LIBS=<scratch>/rlib _R_CHECK_CRAN_INCOMING_REMOTE_=FALSE _R_CHECK_LIMIT_CORES_=TRUE OMP_THREAD_LIMIT=2 _R_CHECK_EXAMPLE_TIMING_CPU_TO_ELAPSED_THRESHOLD_=2.5 _R_CHECK_TEST_TIMING_CPU_TO_ELAPSED_THRESHOLD_=2.5`.
  - No `_R_CHECK_FORCE_SUGGESTS_` override.
  - 331 s wall time.
- **Status line: `Status: 1 NOTE`**
- Findings, verbatim:
  - `checking installed package size ... INFO`
    ```
    installed size is  8.6Mb
    sub-directories of 1Mb or more:
      R      2.0Mb
      doc    4.8Mb
      help   1.6Mb
    ```
  - `checking HTML version of manual ... NOTE`
    ```
    Skipping checking HTML validation: 'tidy' doesn't look like recent enough HTML Tidy.
    Please obtain a recent version of HTML Tidy by downloading a binary
    release or compiling the source code from <https://www.html-tidy.org/>.
    ```
- Everything else is OK. That includes:
  - CRAN incoming feasibility: no version-components NOTE now
  - top-level files: pandoc now available
  - Rd `\usage` sections
  - **PDF version of manual: OK, 8 s / 16 s, so the Rd files are certified**
  - files in `vignettes`
  - re-building of vignette outputs: 22 s / 27 s
- Examples: 13 s / 13 s. With `--run-donttest`: 80 s / 64 s.
- Examples over 5 s elapsed:
  - `run_null_calibration`: user 56.348, system 1.671, elapsed 40.298, ratio 1.44
  - `suggest_thresholds`: user 5.580, system 0.305, elapsed 5.885, ratio 1.00
  - Total over 207 timed topics: user 64.8, system 2.2, elapsed 50.3.
- Tests: `[142s/135s]`, `[ FAIL 0 | WARN 27 | SKIP 77 | PASS 5465 ]`.
  - No test failed under the Step 7 cap.
  - No failure was traceable to CRAN weightedsurv 0.1.0.
- Vignette re-build elapsed: 27 s.
- Installed size: 8.6 MB.
- `00check.log` was copied to `dev/reports/cran_release_00check_2026-10-09.log`.
- The tarball was kept out of the repo. The temp library, build and check directories, and the 0.1.0 extraction were deleted.

## 8. Recorded, not actioned (S1-F)

- `\dontrun` → `\donttest` conversion: 145 Rd files.
- Missing `\examples`: 39 exports and 22 S3 methods.
- `detectCores()`-based defaults not exercised under check: `subgroup.search()`, `mrct_region_sims()`, `fs_fdr_report()`, `resolve_cv_parallel_args()`.
- `subgroup.search()` sets `future::plan(multisession)` without restoring the previous plan.
- `cox_ahr_cde_analysis()` `output_dir = "plots"`. Opt-in via `save_plots`.
- Unconditional `set.seed()` sites (28).
- Spelling: no `inst/WORDLIST`.
- Internal `plot_subgroup()` resets `par(mfrow = c(1, 1))` instead of restoring the old `par`.
- Deferred API items: `sg_objective`, `max.minutes`.

## 9. Post-conditions

- **PC1: PASS (amended).** Every file changed since `dc3ef230^` is in the allowed set.
  - **Amendment (Larry, 2026-10-09):** PC1 also allows `R/summary_utility_functions.R` as a Step 4 messaging-only edit. Its single changed line is the `rmst_calculation()` error string (§4).
  - The other changed R/ files are all from Steps 2, 4, 6 and 7: the 6 Step-2 files; `plot_km_band_forestsearch.R` and `plot_sg_weighted_km.R`; `leff_calibration_io.R`; `utils_check_threads.R`; and the 4 files with grf forest calls.
- **PC2:** `NAMESPACE` is byte-identical to `541fef30`.
- **PC3:** `git status --porcelain` is clean. No tarball, `.Rcheck` directory or temp library is in the repo.
- **PC4:** `Version: 0.4.0`. The first `NEWS.md` heading is `# forestsearch 0.4.0`.
- **PC5:** the committed log's status line is `Status: 1 NOTE`.
- **PC6:** the user-library weightedsurv `RemoteSha` is `edd61a2cecf4b08352719af0a24b10fa1e0b041c` at both start and end.
