# REPORT — forestsearch code audit (read-only), 2026-11-09

Task: `dev/tasks/TASK_code_audit_2026-11-09_v2.md`. Branch `fs-code-audit`. The audit finds and ranks; it changed no code.

## 1. Header

| Item | Value |
|---|---|
| SHA audited | `5b99417776d1c59e4e60510824d2af6c828d4745` (`master` after Step 0.3; the task-doc commit; `R/`, `tests/`, `man/`, `DESCRIPTION`, `NAMESPACE` are identical to its parent `71a7fc1e`) |
| Version (DESCRIPTION) | 0.4.1.9002 |
| Date | 2026-11-09 |
| Files | `R/`: 125 `.R` files, 793 function definitions (785 distinct names; `%||%` is defined in 9 files). `tests/testthat/`: 84 entries (78 `test-*.R`, 4 `helper-*.R`, `fixtures/`, `Rplots.pdf`) |
| `fs-articles-wrappers` | Merged: `master` was fast-forwarded to `46a5ddd5` (reflog `refs/heads/master@{16}: merge fs-articles-wrappers: Fast-forward`); `46a5ddd5` and the A4 commit `ae0af1bb` are ancestors of `master`; no local or remote branch of that name remains |
| Installed copy | forestsearch 0.4.1.9002 built 2026-11-06; all 787 namespace functions deparse identically to the audited tree, so multisession workers in the reproducers ran the audited code |
| Reproducer runs | 201 runs of 196 scripts; total wall clock 713.0 s (11.9 min) against the 20 min cap; longest single run 34.7 s against the 120 s cap |

What ran and what did not:

- Step 0 on `master`: status check, task doc copied and committed (`5b994177`), branch `fs-code-audit` created from it.
- Inventory: every file in `R/` parsed (793 function definitions, all inventoried); exported status from `NAMESPACE`; test references by whole-word grep under `tests/testthat/`; static call graph from `codetools::findGlobals()` over each parsed body (199 Tier 1 functions in 26 files).
- Static sweep (main session): grep for the dimension A and F patterns; `NAMESPACE` imports against use; `DESCRIPTION` Suggests against use and `requireNamespace()` guards; worker counts in examples and tests; roxygen-stated defaults against formals in every file; `codetools::checkUsageEnv()` on the namespace loaded with `pkgload::load_all()` from `tempdir()` (296 messages: 240 NSE bindings, 54 unused locals, 2 NSE function names; no partial-argument-match or possible-error message). `lintr` is not installed and was not used.
- Deep read and sweep: 29 reader agents (13 Tier 1 line by line, 16 Tier 2) returned 545 candidate findings. Every candidate was re-traced at source by an independent verifier; each S1 candidate by two (a refutation pass, then a configuration-and-consequence pass). 6 were refuted, 2 set out of scope, 507 kept; 108 duplicates were merged into 29 groups; 3 defects inside the closed `field_uniform` option were moved to section 10.
- Each invariant was adjudicated by a dedicated verifier over all readers' trails, and invariant 5 also by reproducer R03.
- Reproducers: see the header row above and the appendix. Each ran from `tempdir()` with `pkgload::load_all()` and wrote nothing into the repo; `git status --porcelain` of the forestsearch tree was empty after every batch.
- Not run: the test suite, `devtools::check()`, `R CMD check`, `document()`, `build()`, pkgdown, vignette or article renders, covr, package installs, profiling.
- Incident outside this repository: one sweep script ran with the session's working directory set to the `mrct-consistency` checkout and wrote two scratch CSV files there. They were moved to the scratch area at once; that checkout's status is clean. Nothing was written into the forestsearch tree.

## 2. Gate checklist

| Gate | Result | Evidence |
|---|---|---|
| Step 0.2: `git status --porcelain` empty on `master` | PASS | 0 lines |
| Step 0.3: task doc copied to `dev/tasks/` and committed on `master` | PASS | `5b994177 dev: add code-audit task v2`; v1 was not in `dev/tasks/` and nothing was done about it |
| Step 0.6: no branch `fs-code-audit`; no `REPORT_code_audit_2026-11-09.md` on any local branch | PASS | branch list and `git cat-file -e` over all local branches |
| Step 0.7: `fs-code-audit` created from `master` | PASS | from `5b994177` |
| Inventory: functions inventoried = functions parsed | PASS | 793 = 793 |
| Step 3: each of the five invariants has a verdict, evidence class and trail | PASS | section 4 |
| Read-only: no change under `R/`, `tests/`, `man/`, `vignettes/`, `inst/`, `DESCRIPTION`, `NAMESPACE` | PASS | staged diff against `master` before the commit lists only the two report files |
| Every S1 is reproduced or traced (Step 4) | PASS | S1 rows: A-001 reproduced, A-002 reproduced, A-003 reproduced, A-004 reproduced |
| Each reproducer within 120 s; total within 20 min | PASS | longest 34.7 s; total 713.0 s; 0 prepared scripts not run |
| Step 6.4: CSV rows = ranked-table rows; no S1 with evidence = suspected | PASS | both 414; checked when this file was written and again after the commit (Step 6.7 print) |
| Step 6.2: `git diff --name-status master..fs-code-audit` is exactly two `A` lines | see note | checked on the staged tree before the commit and re-run after it; the post-commit result is printed at close (Step 6.7) |
| Step 6.3: `git status --porcelain` empty after the commit | see note | as above |

## 3. Summary counts

Severity by conducted-analysis path:

| Severity | on path | off path | total |
|---|---|---|---|
| S1 | 2 | 2 | 4 |
| S2 | 80 | 85 | 165 |
| S3 | 69 | 115 | 184 |
| S4 | 19 | 42 | 61 |
| total | 170 | 244 | 414 |

Dimension by severity:

| Dimension | S1 | S2 | S3 | S4 | total |
|---|---|---|---|---|---|
| A | 4 | 159 | 86 | 13 | 262 |
| B | 0 | 0 | 11 | 13 | 24 |
| C | 0 | 1 | 53 | 8 | 62 |
| D | 0 | 1 | 18 | 12 | 31 |
| E | 0 | 0 | 4 | 0 | 4 |
| F | 0 | 0 | 10 | 5 | 15 |
| G | 0 | 0 | 0 | 10 | 10 |
| INV | 0 | 4 | 2 | 0 | 6 |

Evidence class:

| Evidence | S1 | S2 | S3 | S4 |
|---|---|---|---|---|
| reproduced | 4 | 151 | 0 | 0 |
| traced | 0 | 12 | 173 | 59 |
| suspected | 0 | 2 | 11 | 2 |

Fix class (alters_numbers = yes in brackets):

| Fix class | S1 | S2 | S3 | S4 |
|---|---|---|---|---|
| add-only | 0 [0] | 38 [0] | 98 [0] | 24 [0] |
| moves code | 0 [0] | 0 [0] | 13 [0] | 20 [0] |
| changes behaviour | 4 [3] | 123 [11] | 68 [3] | 17 [0] |
| method-adjacent | 0 [0] | 4 [1] | 5 [0] | 0 [0] |

## 4. Invariants

| # | Invariant | Verdict | Evidence | Finding | Trail (short) |
|---|---|---|---|---|---|
| 1 | With use_grf = FALSE and use_lasso = FALSE no model-based pre-processing occurs | holds | traced | none; related ungated input D-001 | GRF fits only under use_grf (forestsearch_main.R:3071) or subgroup_method = "grf"; glmnet only in lasso_selection() under use_lasso (get_fsdata.R:326); otherwise cuts come from covariate quantiles and caller-supplied cuts. |
| 2 | future::plan, options() and RNG state restored on exit, including on error | violated | traced | INV-2; also INV-2a, INV-2b | subgroup.search() replaces the caller's plan on every forestsearch() run (subgroup_search.R:174, :191); subgroup.consistency() calls set.seed(8316951) with no restore (subgroup_consistency_main.R:497-499); further sites in the bootstrap, MR and simulation code. |
| 3 | Every documented default equals the formal default | violated | traced | INV-3; also A-153 | 12 sites in 11 functions plus three nested list-element defaults, headed by forestsearch(use_twostage) documented FALSE with formal TRUE (forestsearch_main.R:928 vs :1700) and mr_inference_args$field_complement documented FALSE, resolved TRUE (:1001 vs :3885). |
| 4 | Deprecated arguments are inert | holds | traced | none | vi.grf.min, max_n_confounders and max.minutes are read only by three warning() guards (forestsearch_main.R:1744-1768); max.minutes is dropped before subgroup.search() (:3493); none reaches get_FSdata(), the GRF/DINA builders or subgroup.consistency(). |
| 5 | Sequential and multisession give the same selection for the same seed | violated | reproduced (V110b; R03) | INV-5; also INV-5a | Holds on the default path (closed-form resample rate, no RNG drawn; R03: same selection and same Pcons). Violated wherever literal splitting runs (consistency_method = "split", effect_measure = "IRD", NA closed-form fallback): one seeded stream sequentially versus per-candidate L'Ecuyer streams in parallel (subgroup_consistency_main.R:717-723 vs :970-999). V110b: under split, per-candidate Pcons differed between plans in 3 of 3 datasets and the selected subgroup differed in 1 (data seed 203: {b2} & {b3} sequential, {b1} & {b3} multisession); under the default resample all comparisons were identical. R03 (one dataset): Pcons differs under split, same selection. |

Full trails:

- **Invariant 1 (holds).**
  - Re-read at source: forestsearch() body 1722-2780 and 3066-3500 line by line, get_FSdata() 84-748 line by line, plus a package-wide grep of every model-fitting call (grf::, glmnet::, policytree::, dina/dina_fit, coxph/glm/lm) mapped to its enclosing function and tier.
  - Formals: forestsearch_main.R:1638 use_lasso = FALSE, :1639 use_grf = FALSE, :1642 use_dina = FALSE, :1647 subgroup_method default "consistency".
  - GRF: the only GRF fits on the forestsearch() graph are forestsearch_main.R:3093 / :3113, inside `if (use_grf && ...)` at :3071, and forestsearch_helpers.R:1877 / :1901 inside .forestsearch_grf_select(), reached only from :2945 under `if (subgroup_method == "grf")` (:2892, returns at :3063). grf::variable_importance appears only at grf_subg_harm_glm.R:476 (inside grf.subg.harm.glm, display). randomForest is imported in roxygen (:1618) but never called in R/.
  - GRF cuts in get_FSdata(): get_fsdata.R:252 turns use_grf off when grf_cuts is NULL; grf_cuts is consumed only under use_grf (:437, :440, :447, :495-507); grf_cuts_keep stays character(0) (:492) so the protect vector (:570) is unaffected. A caller-supplied grf_cuts / grf_res with use_grf = FALSE is therefore inert (apart from the type check at :263-268).
  - Lasso: the only glmnet calls on the path are get_FSdata_helpers.R:319 and :326 inside lasso_selection(), whose single caller is get_fsdata.R:328 inside `if(use_lasso){` (:326). lassokeep/lassoomit are NULL otherwise (:324-325); the filter at :355-356 is inside that block; the remaining lasso branches are guarded at :362 and :386; :395 and :423 are the non-lasso branches.
  - Non-model pool construction with both flags FALSE: is.continuous() unique-value count (get_FSdata_helpers.R:11), mean/median/quartile/J-quantile cuts (get_conf_force :129, get_conf_force_jq :183, median at get_fsdata.R:480), conf_force / exclude_cuts supplied by the caller, and collapse_redundant_cuts() using sd(x)/sqrt(n) of the covariate only (get_FSdata_helpers.R:692-693). No coxph/glm/lm/grf/glmnet call in get_fsdata.R or get_FSdata_helpers.R other than lasso_selection (grep).
  - After get_FSdata(): forestsearch_main.R:3374 conf.screen <- FSconfounders.name (no importance screening or ordering), :3386 df.confounders <- df[, conf.screen], :3496 subgroup.search. The Cox/GLM fits in subgroup_search.R and the consistency files are the per-candidate estimators of the method, not pre-processing.
  - DINA (third front end, not named in the invariant's premise): fit at forestsearch_main.R:3206 only inside `if (use_dina && is.null(dina_cuts))` (:3180), or via .forestsearch_dina_select() under `if (subgroup_method == "dina")` (:2761, returns at :2878). With use_dina = TRUE and use_grf = use_lasso = FALSE a fitted model does shape the candidates (get_fsdata.R:517-526), by explicit request; the code's own statement of the contract names all three switches (forestsearch_main.R:1741-1745, :3370-3373; .fs_family_status forestsearch_helpers.R:2182-2185).
  - Propensity model: is.RCT = FALSE on a GLM outcome resolves ps_method to "grf" (:2598) and fits estimate_propensity_scores() (:2688; glm_effect_estimators.R:1479 / :1496 / :1519). Its outputs are the columns ps_hat, sw, ips_covar (glm_effect_estimators.R:1450-1452); confounders.name is unchanged and get_FSdata() builds cuts from confounders.name only, so cuts and candidate factors are untouched. On a survival outcome this configuration stops at :2619.
  - Reader reports: 7 'holds' and 2 'not established' (fs-helpers, oc); both 'not established' were scope limits of those readers (they had not read forestsearch() / get_FSdata()), not counter-evidence. No reader claimed a violation. Every cited line was re-opened and matches.
- **Invariant 2 (violated).**
  - Site list from grep of all 125 files in R/ for plan(, options(, set.seed(, .Random.seed, RNGkind(, par(, handlers(, setDTthreads, setwd, Sys.setenv, sink, registerDo*, withr (the last six: no hits), each site then opened.
  - PLAN - NOT restored:
  - subgroup_search.R:174 future::plan(multisession, workers = parallel_workers), :191 future::plan(sequential); subgroup.search() has no on.exit and never reads the previous plan. Reached on every forestsearch() run: forestsearch_main.R:3455 (parallel_workers = 1) -> :3496 do.call(subgroup.search, ...). forestsearch_main.R contains no plan(, on.exit, set.seed or options( call. Normal exit: caller's plan is sequential. Error inside search_combinations_parallel (:176-188): plan left multisession, and forestsearch() swallows the error at :3495-3500 and carries on.
  - bootstrap_dofuture_main.R:366 setup_parallel_SGcons(parallel_args); exit handler :349-356 does plan(.Last.future.plan) if that object exists, else plan("sequential"). `.Last.future.plan` is assigned nowhere in R/ (only declared in globals.R:30), so the entry plan is never captured; exit state is sequential.
  - mrct_simulation.R:358/365/368/371/375/377 set the plan in mrct_region_sims(); :704 plan(sequential) on normal exit only; no on.exit in the function.
  - fs_fdr_report.R:398-399 (.fs_fdr_run_null) plan(sequential) then plan(multisession); :411 plan(sequential) on normal exit only; no on.exit.
  - fpr_calibration.R:272 plan(multisession) when n_workers > 1; on.exit :273-276 sets plan("sequential"), not the previous plan (error-safe, wrong target).
  - PLAN - restored (capture + on.exit registered before the change, so also on error): subgroup_consistency_main.R:883-884 (parallel branch only; sequential branch :717-722 never touches the plan); forestsearch_cross_validation.R:332-333 and :840-841; run_subgroup_sims.R:817-819 and :829-834.
  - PLAN - setters by documented purpose (no restore expected): setup_parallel_SGcons() subgroup_consistency_helpers.R:57/60/66/69; reset_workers() forestsearch_helpers.R:806/834/836.
  - PLAN - inside a %dofuture% body (effect on the calling session under a sequential plan not established): fpr_calibration.R:296, fs_fdr_report.R:406.
  - OPTIONS: bootstrap_dofuture_main.R:347 options(warn = -1), restored by on.exit :349-350 (and inline :494) - holds. run_subgroup_sims.R:880-881 saved and restored by on.exit - holds. mrct_simulation.R:392 progressr::handlers("progress") changes the session's progressr handler setting and is not restored (mechanism inside progressr not read at source). No other options() call in R/.
  - RNG - restored by on.exit (also on error): get_FSdata_helpers.R:248-259 (lasso_selection); dina.R:420-433 (.dina_make_folds); dina_bagged.R:423-436; subgroup_consistency_main.R:44-58 (.make_candidate_rng_seeds, kind and state); mrct_case_study.R:263-271 with RNGkind at :280 (kind and state). Partial: fs_dgm_feasibility.R:162-166 restores only when a seed existed at entry (set.seed at :181 otherwise leaves a new .Random.seed). Inline restore, not on.exit: generate_glm_dgm.R:416-422, setup_gbsg_dgm.R:128-134 (and both functions reseed unconditionally before it).
  - RNG - NOT restored, conducted-analysis path:
  - subgroup_consistency_main.R:497-499 set.seed(seed), formal seed = 8316951 (:397); forestsearch() passes no seed (overrides :3567-3592) -> every forestsearch() run that reaches the consistency stage (:3655) leaves the session RNG reseeded. The only on.exit in subgroup.consistency() is the plan one (:884).
  - fs_mr_inference.R:615 and :912 (reached from forestsearch_main.R:3844 with seed = seedit at :3912 when mr_inference = TRUE); fs_mr_field_uniform.R:93 (from fs_mr_inference.R:1020).
  - grf_subg_harm_glm.R:460 set.seed(seedit) (from forestsearch_main.R:3113 and forestsearch_helpers.R:1901, GLM outcomes with GRF).
  - bootstrap_calculations_helpers.R:27 (.make_boot_index_matrix, from bootstrap_dofuture_main.R:383 and from exported bootstrap_ystar :103-104); bootstrap_analysis_dofuture.R:334 (bootstrap_results, main process).
  - consistency_resample.R:459 (only with a user seed and method mc/both; the forestsearch callers at subgroup_consistency_helpers.R:1487/1508/1870/1874 pass no seed and method = "closed").
  - RNG - NOT restored, off path (Tier 2): consistency_resample.R:505; forestsearch_cross_validation.R:369; frontier_cis.R:164; fs_oc_predict.R:233; fs_oc_grid.R:311, :492; guohe_algorithm3.R:395; guohe_adaptive_r.R:205; generate_aft_dgm_main.R:248; generate_glm_dgm.R:402; simulate_from_dgm.R:208; simulate_from_glm_dgm.R:97; sim_aft_gbsg.R:478, :604, :856, :884, :900, :908; synthetic_data_perturbation.R:98, :506; mrct_inject_dgm.R:97, :277, :572, :664; run_simulation_analysis.R:509; truefind_asymptotic_glm.R:338. None of these files has an on.exit or .Random.seed handling for the seed.
  - RNG - inside a %dofuture% body or per-trial worker (effect on the calling session under a sequential plan not established): fpr_calibration.R:301; forestsearch_cross_validation.R:958; run_subgroup_sims.R:512.
  - PAR (in the task's site list, not in the invariant as worded): restored by on.exit at cox_ahr_cde_wrapper.R:377-378, cox_spline_fit.R:128-129 and :598-599, glm_spline_fit.R:154-155, find_k_inter_main.R:283-284, plot_sg_results.R:562-563/692-693/805-806, plot_sg_weighted_km.R:364-365, mrct_inject_dgm.R:706-707, mrct_trial_helpers.R:236 and :259-260. NOT restored: subgroup_consistency_helpers.R:976 par(mfrow = c(1, 2)) in plot_subgroup(), called at :1062 when plot.sg = TRUE.
  - DATA.TABLE THREADS: no setDTthreads call in R/.
- **Invariant 3 (violated).**
  - All paths under R/. Every reader-claimed mismatch was re-opened; two package-wide sweeps were added (roxygen @param blocks vs formals: 2789 params / 481 functions; man/*.Rd \arguments vs \usage: 402 pages / 2718 args, which also covers @inheritParams and @rdname).
  - A. Formal has a literal default, roxygen states a different value (6):
  - forestsearch(use_twostage): doc forestsearch_main.R:928 'Default FALSE for backward compatibility'; formal :1700 use_twostage = TRUE.
  - get_Cox_sg(cox_initial): doc Cox_estimation_helpers.R:8-9 'pre-fitted Cox model object ... Default NULL'; formal :33 cox_initial = log(1), used as coxph init at :53.
  - calibrate_glm_interaction(k_inter_range): doc calibrate_glm_interaction.R:48 'Default c(0, 10)'; formal :103 c(-10, 10).
  - compute_detection_probability_glm(c2): doc truefind_asymptotic_glm.R:180 'Default: 1.0 for ratio, 0.0 for difference'; formal :241 c2 = 1.0, no scale-dependent branch (:249-292; difference scale takes k2 <- c2 at :279). All four in-package callers pass c2 explicitly (calibrate_null_correction.R:180, :373; fpr_calibration.R:250, :527).
  - plot_subgroup_effects(hrz_crit): doc cox_spline_fit.R:537 'Default is 1 (HR=1 on log scale is 0)'; formal :590 hrz_crit = 0, compared on the log scale at :640.
  - dina_subgroup(n_min.frac): doc dina_subgroup.R:166 'Default 0.10' followed by a stray :168 'Default 60L'; formal :321 n_min.frac = 0.10.
  - B. Formal is a NULL sentinel and the body-resolved default differs from the documented one (3):
  - fpr_approx(thresholds): doc fpr_approximation.R:57 'Default c(1.0, 1.5, 2.0)'; formal :103 NULL; :165 sort(unique(c(c1, thresholds))) with no NULL substitution, so only c1 is evaluated.
  - gg_forest(widths): doc gg_forest.R:92 'Default: c(3.5, 5, rep(1, n_annot))'; formal :141 NULL; resolved :434-435 to c(3.5, 5, rep(1.1, n_annot)).
  - interpret_estimation_table(scenario): doc simulation_tables.R:1043 '"null" or "alt" (default)'; formal :1092 NULL; inferred from dgm at :1158-1169 (can resolve to "null").
  - C. Documented default, formal has none (1 exported + inherited text):
  - bootstrap_results(show_three): doc bootstrap_analysis_dofuture.R:78-79 'Default: FALSE'; formal :301 has no default; evaluated unconditionally at :428.
  - Inherited via @inheritParams onto internal helpers with default-less formals: compute_detection_probability_single (truefind_asymptotic.R:262 inherit, formals :270-279; prop_cens 0.3, method "cubature", n_mc 100000, tol 1e-4 from :127-137) and density_threshold_integrand (:98 inherit, formal :103; prop_cens).
  - Orphaned duplicate roxygen block subgroup_consistency_helpers.R:1170-1223 (states defaults NULL / "split") fused onto .make_eval_subgroup_consistency (:1239-1245, no defaults, @noRd :1238); not rendered.
  - D. Nested list-element defaults documented in roxygen vs the code that resolves them (3):
  - forestsearch mr_inference_args$field_complement: doc forestsearch_main.R:1001 'default FALSE ... reproduces prior output exactly'; resolved :3885 to TRUE; fs_mr_inference formal fs_mr_inference.R:601 TRUE (its own doc :352 says TRUE).
  - forestsearch mr_inference_args$ci_method: doc forestsearch_main.R:978 '"field" (default)'; true on the consistency branch (:3870) but "ij" on the DINA/GRF branches (fs_mr_inference_methods.R:171, called from forestsearch_main.R:2860, :3045).
  - subgroup.consistency parallel_args$batch_size: doc subgroup_consistency_main.R:276-279 says 1 for sg_focus "hr"/"minSG"; code :908-913 gives 1 for "maxeffCons" only.
  - E. Cross-reference: run_subgroup_sims.R:200-201 (subgroup_glm) says make_effect_estimator()'s binary default is "RD"; glm_effect_estimators.R:120-123 resolves binary to "OR" (its own doc :47 is right).
  - Literal-only, not counted (formal NULL, body resolves to exactly the documented value): mrct_region_sims(conf_force) mrct_simulation.R:54 vs :271, :326-328; load_leff_calibration(output_dir) leff_calibration_io.R:123 vs :138, :141-143. Not counted: forestsearch_Kfold(Kfolds) forestsearch_cross_validation.R:190 vs :289 (.fs_cv_nobs equals nrow(df.est) when a subgroup was found, :148-154); fs_mr_inference.R:332/:339 kappa in [1,2] and M <= 12 vs fs_mr_field_uniform.R:87/:90 (40L, 1 to 3) is prose on the closed field_uniform option.
  - Checked and matching: every other forestsearch() @param default (formals :1622-1720), including ps_method (:2596-2598), collapse_cuts_args (get_fsdata.R:551-552), twostage_args keys (subgroup_consistency_main.R:518-524), and the remaining mr_inference_args keys (:3853-3912).
- **Invariant 4 (holds).**
  - Re-checked at source (tree 5b994177, branch fs-code-audit, clean). All 10 reader reports under key 4 were read: 8 'holds', 2 'not established' (boot-core, cv) only because they had not traced subgroup_search.R:181/274 and the code downstream of forestsearch_main.R:1768; both gaps are closed below. No reader claimed a violation.
  - Complete non-roxygen site list for the three names (grep over all of R/, also for alias spellings minutes / n_confounders / vi.grf / vi_grf):
  - forestsearch_main.R:1648, 1689, 1698 -- formals (max_n_confounders = 1000, max.minutes = 3, vi.grf.min = NULL).
  - forestsearch_main.R:1744-1751 -- if (!is.null(vi.grf.min)) warning(...) only.
  - forestsearch_main.R:1752-1758 -- if (!missing(max_n_confounders) && !identical(as.numeric(max_n_confounders), 1000)) warning(...) only.
  - forestsearch_main.R:1761-1768 -- if (!missing(max.minutes) && !identical(as.numeric(max.minutes), 3)) warning(...) only.
  - forestsearch_main.R:1880-1881 -- args_call_all <- mget(names(formals())): the values are recorded (by design, comment :1742-1743) and returned in $args_call_all (:1927, :2807, :2996, :3330, :3504, :3964).
  - forestsearch_main.R:3493 -- search_args$max.minutes <- NULL before do.call(subgroup.search, search_args) at :3496.
  - subgroup_search.R:100 (formal, default 30), :112-118 (warning when supplied), :181 (forwarded), :274 (formal of search_combinations_parallel); body :282-350 never reads max.minutes or t.start (t.start is used only at :171/:201 for the elapsed-time print).
  - fs_family_report.R:142-143, 164, 266-269, 353-354 -- reporting only; status is the constant INE ('inert') whatever the value; fs_family_report() is not called by forestsearch().
  - forestsearch_main.R:1736-1741, 1760, 3372, 3491 are comments.
  - Why the recorded values cannot reach anything else in forestsearch()'s call graph:
  - Whole-list consumers of args_call_all inside forestsearch(): .validate_mr_configuration (forestsearch_helpers.R:2070-2133) and .fs_family_status (:2174-2186) read by exact name with config[[nm]] (subgroup_method, dina_select_statistic, grf_selection, grf_select_statistic, use_lasso, use_grf, use_dina) -- none of the three.
  - filter_call_args() (summary_utility_functions.R:679-688) keeps only names that are formals of the target; get_FSdata (get_fsdata.R:84-91) and subgroup.consistency (subgroup_consistency_main.R:376-408) have no such formals and no '...'. Call sites forestsearch_main.R:3311-3320 and :3629-3655.
  - Search: :3487-3490 keeps only names(formals(subgroup.search)), then :3493 drops max.minutes, so subgroup.search() runs with its unforced default and does not warn.
  - GRF and DINA calls are built from explicit named arguments (:3078-3113 via .build_grf_survival_args/.build_grf_glm_args, grf_args.R:25-44, 86-102; :3198-3220 from da$fit / da$frontier) -- none of the three.
  - .sync_args_call_all (:21-30) copies only explicitly named variables (:2149, 2461, 2535, 2632, 2716, 3307).
  - No dynamic access that could pick them up: no match.call/sys.call/parent.frame/dynGet/as.list(environment()) in forestsearch_main.R or its callees other than the mget at :1881; no partial-matching accessor ($max..., $vi..., [[..., exact = FALSE]]) on an argument list anywhere in R/.
  - Replays (outside the invariant's scope, checked because the values are replayed): bootstrap_analysis_dofuture.R:614 try(do.call(forestsearch, ...)) inside suppressWarnings (:408); forestsearch_cross_validation.R:479, :1001 suppressWarnings(try(...)); compare_selection_rules.R:560-572, calibrate_null_correction.R:357-360, fpr_calibration.R:307, run_simulation_analysis.R:818-821 catch errors only; mrct_trial_helpers.R:46-52 and mrct_simulation.R:570-578 use calling handlers that muffle or pass; fs_fdr_report.R:361-376 promotes a warning to a failure only when its text matches 'subgroup.consistency|subscript', which none of the three deprecation texts (forestsearch_main.R:1745-1750, 1754-1757, 1763-1767) contains. So a stored non-default value cannot fail or alter a replicate.
- **Invariant 5 (violated).**
  - WHAT IS PLAN-INDEPENDENT (holds by reading)
  - forestsearch_main.R:2059-2092: parallel_args is only validated in forestsearch(); its single consumer is subgroup.consistency(), reached through args_call_all (:1880-1881) and filter_call_args (:3629-3633; summary_utility_functions.R:679-688).
  - Search stage: forestsearch_main.R:3455 forces parallel_workers = 1 for every call; subgroup_search.R:174 sets multisession(workers = 1), :285 future_lapply without future.seed, :191 resets to sequential. No RNG call exists in subgroup_search.R (grep for sample/runif/rnorm/rbinom/set.seed: none). The candidate table is therefore the same under either parallel_args.
  - GRF / lasso / DINA pre-processing and MR run in the calling process with explicit seeds and never read parallel_args.
  - Default consistency path: evaluate_subgroup_consistency() subgroup_consistency_helpers.R:1485-1491 and evaluate_consistency_twostage() :1863-1913 call consistency_resample(method = "closed"); consistency_resample.R:455-457 computes 2*pnorm(delta/sigma_D)-1 with no draw (the only draws are in the mc branch :458-468, not requested). Both the plain loop and the worker closures (.make_eval_* helpers:1242-1294) call the same evaluators, so Pcons, admission and sort are identical.
  - sg_focus = "maxeff": subgroup_consistency_main.R:758-786 evaluates in-process on both plans and the key (-hr, K) has no Pcons term.
  - Early stopping through forestsearch(): stop_threshold is reset to NULL for hr/maxSG/minSG/hrMaxSG/hrMinSG (forestsearch_main.R:2116-2148) and for maxeff (:2027); for maxeffCons the key (-hr, K) equals the enumeration order, so extra candidates in a parallel batch cannot outrank the first qualifier.
  - WHERE THE TWO PLANS DIVERGE
  - subgroup_consistency_main.R:397 seed = 8316951 (forestsearch() has no `seed` formal and consistency_overrides :3568-3592 passes none, so seedit never reaches this stage); :497-499 set.seed(seed) on the session RNG.
  - :702 use_parallel; :717-723 plan == "sequential" sets use_parallel <- FALSE -> plain loop :793-871: all candidates consume the single stream from :498 in evaluation order.
  - Any other plan -> :878-1000: :970-974 .make_candidate_rng_seeds(seed, n_candidates) (:40-66, set.seed under L'Ecuyer-CMRG then parallel::nextRNGStream per candidate), :991-999 future_lapply(..., future.seed = candidate_seeds[batch_indices]).
  - Literal splits draw at subgroup_consistency_helpers.R:272 (run_single_consistency_split, via .consistency_via_splits :355-375) and :1733 (two-stage embedded copy, Stage 1 :1919-1925, Stage 2 :1982-1988). Different generator and stream -> different splits -> different Pcons for the same candidate.
  - Effect on selection: admission helpers:1541 / :2021 (Pcons < pconsistency.threshold -> NULL), two-stage screen :1939 and Wilson stop :1965-1980; sort_subgroups() :586-587 orders the default sg_focus = "hr" by (-Pcons, -hr, K), so Pcons is the PRIMARY key; sg_consistency_out() :1045-1049 takes row 1.
  - CONFIGURATIONS UNDER WHICH THE SELECTION CAN DIFFER (forestsearch(), same seed, parallel_args plan sequential vs multisession)
  - 1. consistency_method = "split", any outcome type, any sg_focus except "maxeff" (helpers:1529-1534; two-stage falls through :1863-1868 to :1919). use_twostage = TRUE is the forestsearch() default (:1700).
  - 2. Default consistency_method = "resample" with a GLM effect_measure = "IRD": .glm_resample_supported() helpers:343-346 excludes IRD, so helpers:1502-1504 fails and :1529-1534 splits (spec always supplied at forestsearch_main.R:3616-3625).
  - 3. Default "resample" when rate_closed is NA for a candidate: helpers:1492-1500 (Cox), :1520-1527 (GLM), two-stage :1890 -> :1919. NA arises when the pieces are NULL: consistency_resample.R:62, :74, :76 (any NA coefficient, e.g. an adjust_covariates term aliased inside the subgroup), :82, :87 -> :438.
  - In 2 and 3 the statement 'Results under consistency_method = "resample" are identical on every path (the closed form consumes no RNG)' (roxygen subgroup_consistency_main.R:263-264; comment :714-715) does not hold.
  - Not executed; the size of the difference is Monte-Carlo noise in Pcons at the stated fs.splits.

## 5. Ranked findings

Order: S1 to S4; on-path before off-path; reproduced, traced, suspected. A location list ("Also at") marks one finding that repeats at several sites. Line numbers are at the audited SHA.

| id | tier | dimension | file | line | function | severity | on_path | evidence | summary | alters_numbers | fix_summary | fix_class |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| A-001 | 1 | A | forestsearch_main.R | 2766 | forestsearch | S1 | on | reproduced | For a continuous (gaussian/MD) outcome with both threshold spellings left unset, both DINA routes take the harm floor from the raw formal hr.threshold = 1.25 instead of the resolved MD screening threshold 0.0, so the fit builds its candidate family under a 1.25-outcome-unit floor while the returned object reports a floor of 0 and every bootstrap/CV replay of args_call_all runs DINA with floor 0. Also at: forestsearch_main.R:3237 (use_dina screening, selected_only = TRUE: m_diff_sel <- hr.threshold for gaussian); forestsearch_helpers.R:1525 (consumer: m_diff <- hr.threshold for gaussian); ... (full list in the CSV) | yes | Derive the DINA floor at both sites (the subgroup_method = "dina" call and the use_dina selected-cut screening) from the resolved natural-scale screening threshold that is recorded in threshold_config and written to args_call_all$effect.threshold, not from the raw hr.threshold formal. | changes behaviour |
| A-002 | 1 | A | forestsearch_main.R | 3722 | forestsearch | S1 | on | reproduced | On the default consistency path forestsearch() attaches treat.recommend to df.est and df.predict by a merge on the id column, so inside every full-bootstrap replicate (resample keeps the original ids, repeated k_i times) df.est has k_i^2 rows per subject and df.predict is the in-bag subjects only (k_i copies each, out-of-bag subjects NA), which makes the replicate's 'new subgroup on original data' estimate the bootstrap-sample estimate and its 'new subgroup on bootstrap data' estimate a k_i^2-weighted one. Also at: forestsearch_main.R:3733 (df.predict merge); flag table subgroup_consistency_helpers.R:939-942; consumers bootstrap_analysis_dofuture.R:629-630, 651, 667, 688-697, 700-711, 726-728, 741-769, 872-875, 889; ... (full list in the CSV) | yes | On the consistency path attach treat.recommend to df.est by row position and to df.predict by evaluating the selected subgroup definition on its rows (as df.test and the DINA/GRF branches already do) instead of joining on id. | changes behaviour |
| A-003 | 2 | A | run_simulation_analysis.R | 1025 | .extract_fs_estimates_gen | S1 | off | reproduced | For a replicate under an alternative DGM in which nothing is identified, the simulation extractors store spec = \|Hc\|/N (the value of NPV for an empty identified set) where specificity is 1, so the all-replicate specificity shown by the exported OC summaries is too low by about (1 - detection rate) x \|H\|/N; npv is stored as 1, which matches the table-fill convention but not the found-case formula. Also at: run_simulation_analysis.R:1026 (npv = 1), :1066, :1392, :1418 (same assignments in the three other nothing-found branches; :1392 and :1418 are in .extract_grf_estimates_gen); ... (full list in the CSV) | yes | Store spec = 1 in the four nothing-found branches under an alternative, and choose one package-wide value for npv in that case (1, as the table fill and fs_mr_oc_summary use, or \|Hc\|/N, as the found-case formula gives). | changes behaviour |
| A-004 | 2 | A | synthetic_data_perturbation.R | 366 | generate_gbsg_bootstrap_general | S1 | off | reproduced | The exported generate_gbsg_bootstrap_general() loads its source data with data("gbsg", package = "survival"), which loads nothing with survival 3.8-12 (gbsg ships inside the cancer data file), so gbsg stays NULL and every call, including the default call, stops with "'data' must be a data frame". | no | Take the data from survival::gbsg, as the package's other examples do. | changes behaviour |
| A-005 | 2 | A | Cox_estimation_helpers.R | 54 | get_Cox_sg | S2 | on | reproduced | get_Cox_sg() fits coxph() under suppressWarnings() with no per-arm event check, so a slice with a zero-event arm is expected to return a large finite log-HR that passes the bootstrap's is.finite() guard, whereas the search scorer and .make_cox_estimator() return non-estimable for the same slice. Also at: bootstrap_analysis_dofuture.R:478, :507 (the two calls that can meet a zero-event arm: original H / Hc evaluated on a resample); bootstrap_analysis_dofuture.R:724, :767 (guards test is.finite only); ... (full list in the CSV) | no | Apply the existing per-arm existence check (.fs_existence_reason for HR) inside get_Cox_sg(), or in its bootstrap callers, and return NA for a non-estimable slice so the draw is dropped like a non-finite one. | changes behaviour |
| A-006 | 2 | A | bootstrap_analysis_dofuture.R | 436 | bootstrap_results | S2 | on | reproduced | When the exported bootstrap_results() is called directly with its default boot_index_mat = NULL, it draws the resampling indices from per-iteration doFuture RNG streams and does not return them, while bootstrap_ystar() with the same seed draws a different Mersenne-Twister index set, so a direct caller who pairs the two in get_dfRes() (the route the @seealso and the wrapper's seed documentation point to) gets an IJ variance computed against unrelated multiplicities, with no error. Also at: bootstrap_dofuture_main.R:30-36 (stale @param seed: the wrapper no longer calls bootstrap_ystar()); bootstrap_analysis_dofuture.R:250-251 (@seealso), 344-346 (comment), 101-111 (@param boot_index_mat); ... (full list in the CSV) | no | Return the indices actually used (or require boot_index_mat for direct callers) and correct the @seealso, the wrapper's @param seed and the in-code comment that promise seed alignment with bootstrap_ystar(). | add-only |
| A-007 | 2 | A | bootstrap_analysis_dofuture.R | 534 | bootstrap_results | S2 | on | reproduced | The four event counts for the original subgroup on each resample (events_H_0, events_H_1, events_Hc_0, events_Hc_1) are computed and then overwritten with NA before the result row is built, so these documented columns are NA on every replicate and summarize_bootstrap_events() always reports zero low-event resamples for the original H and Hc. Also at: bootstrap_analysis_dofuture.R:535-537 (the other three overwrites); :916-919 (written to the row); :133-136 (roxygen documents integer counts); bootstrap_summaries_helpers.R:1091-1103 (consumer) | no | Remove the four re-initialisations at :534-537, or move the NA initialisation block above the computation at :464. | changes behaviour |
| A-008 | 2 | A | bootstrap_analysis_dofuture.R | 614 | bootstrap_results | S2 | on | reproduced | A bootstrap replicate whose inner forestsearch() call throws is recorded exactly like a replicate that found no subgroup (NA estimates, any_found = 0), its warning is suppressed twice, and if every replicate throws the user is told the selection criteria were not reached. Also at: bootstrap_analysis_dofuture.R:408 (suppressWarnings around the foreach); bootstrap_analysis_dofuture.R:619-621 (the warning that is never seen); bootstrap_dofuture_main.R:346-347 (options(warn = -1)); ... (full list in the CSV) | no | Record each replicate's caught error (a column or attribute) and report the count and first message after the loop, outside the suppressed section. | add-only |
| A-009 | 2 | A | bootstrap_dofuture_main.R | 302 | forestsearch_bootstrap_dofuture | S2 | on | reproduced | For GLM outcomes the bootstrap rebuilds its effect estimator from the stored call without adverse_outcome and without ps_adjust_method, so its observed and bias-corrected estimates are on the unflipped, unweighted scale even when the search that selected the subgroup used adverse_outcome = FALSE (the resolved default for continuous outcomes) or propensity adjustment. Also at: mr_estimates_table.R:204-208 (same omission, but the effect column is then overwritten from fs$mr_inference at :283-286); forestsearch_cross_validation.R:1452-1458, 1508-1515 (reader's citation, not re-verified); ... (full list in the CSV) | yes | Pass the recorded adverse_outcome and ps_adjust_method when rebuilding the bootstrap estimator, or document that the bootstrap reports on the natural (unflipped, unweighted) scale. | changes behaviour |
| A-010 | 2 | A | bootstrap_dofuture_main.R | 533 | forestsearch_bootstrap_dofuture | S2 | on | reproduced | The success guard covers only zero re-identifying replicates; with exactly one, the infinitesimal-jackknife standard error is exactly 0 and a zero-width bias-corrected interval is reported with no warning. Also at: bootstrap_calculations_helpers.R:345-346, 351, 360-371 (get_targetEst returns sehat_new = 0); bootstrap_calculations_helpers.R:156 (ci_est accepts sd = 0) | no | Treat fewer than two usable replicates as having no variance estimate (NA interval with a warning) instead of returning a zero standard error. | changes behaviour |
| A-011 | 2 | A | bootstrap_dofuture_main.R | 736 | forestsearch_bootstrap_dofuture | S2 | on | reproduced | For a survival fit made with est.scale = "1/hr", each row of the bootstrap summary table FSsg_tab pairs one subgroup's label, counts, medians and unadjusted estimate with the other subgroup's bias-corrected estimate (the Qstnbl row, which holds H = treat.recommend == 0, shows Hc_bc, and the Recmnd row shows H_bc). Also at: bootstrap_dofuture_main.R:737 (paired assignment); ... (full list in the CSV) | no | On est.scale = "1/hr" pass each bias-corrected string to its own subgroup's row (H_bc with treat.recommend == 0, Hc_bc with treat.recommend == 1) and, if the reversed reading is wanted, swap the two row labels as the cross-validation and KM-plot code do. | changes behaviour |
| A-012 | 1 | A | consistency_resample.R | 81 | .consistency_cox_pieces | S2 | on | reproduced | For a Cox model with more than one coefficient (a linear adjust_covariates term such as "x2") the treatment dfbeta column is looked up by column name, but survival's dfbeta residual matrix has no column names, so the function returns NULL for every candidate and the resample approximation and MR pieces are silently unavailable. Also at: subgroup_consistency_helpers.R:1491-1500 and :1885-1890 (silent fallback to literal splitting); fs_mr_inference.R:96-97 and :637-641 (every candidate dropped, MR returns a note with harm_flag NA); ... (full list in the CSV) | no | Select the treatment column of the dfbeta matrix by position matched against names(coef(fit)) rather than by the matrix's column names. | changes behaviour |
| A-013 | 1 | A | dina.R | 270 | dina_fit | S2 | on | reproduced | Cross-fit coefficients are combined with colMeans(beta_mat, na.rm = TRUE), so a step-2 coefficient that is aliased (NA) in one fold is silently replaced by the other fold's estimate alone, and one that is NA in every fold becomes NaN, is returned without error, and later makes dina_subgroup() fail with 'missing value where TRUE/FALSE needed'. Also at: dina.R:852 (glm.fit NA coefficient passed on unchecked); dina.R:911-923 (coxph NA coefficient; singular-information branch returns an NA sigma); dina_subgroup.R:749 and :888 (if (mean_tau < m_diff) errors on NaN) | no | Detect NA step-2 coefficients per fold and stop or warn naming the covariate instead of averaging over the remaining folds, and have dina_subgroup() reject a fit with non-finite coefficients with a clear message. | changes behaviour |
| A-014 | 1 | A | dina.R | 501 | .dina_fit_baselines_glm | S2 | on | reproduced | dina_fit(family = "gaussian", baseline_method = list(eta0 = , eta1 = )) is documented and passes validation but always stops with 'could not find function "nu_fn"', because the eta0/eta1 list branch returns nu_fn = NULL and the Gaussian nuisance combiner evaluates only nu_fn. Also at: dina.R:766 (the failing call nu_fn(X_new) with nu_fn bound to NULL); dina.R:86-89 (roxygen documenting the eta0/eta1 list with no family restriction); ... (full list in the CSV) | no | For the Gaussian family with an eta0/eta1 list, either build nu(x) from the supplied functions and the propensity (the a*eta1 + (1 - a)*eta0 combination already used for the other families, with a = e) or reject that list up front with an explicit message and correct the roxygen and the stop text. | add-only |
| A-015 | 1 | A | dina.R | 660 | .dina_fit_censoring_cox | S2 | on | reproduced | cens_type = "unknown" fits one pooled logistic model of status on X and returns the same function for both arms, so the two not-censoring probabilities cancel in a(x) and the option gives the same estimator as cens_type = NULL up to the eps = 1e-6 stabiliser. Also at: dina.R:648-659 (single pooled fit, no arm term); dina.R:801-805 (identical probabilities cancel in a(x)); dina.R:629 (helper header states a per-arm contract); roxygen dina.R:95-98 | no | Either make the 'unknown' branch produce arm-specific probabilities as the helper's stated contract says, or document that the option is currently equivalent to no censoring adjustment; which of the two is the author's method decision. | method-adjacent |
| A-016 | 1 | A | dina_df.R | 244 | .dina_extract_from_df | S2 | on | reproduced | dina() and dina_bagged() document `df` as 'data frame (or tibble / data.table)', but the shared extraction helper selects columns with `df[, required_cols]`, which in this data.table-importing namespace is evaluated with data.table semantics and stops for every data.table input. Also at: dina_df.R:256 (same idiom; not reached for a data.table because 244 stops first); ... (full list in the CSV) | no | Coerce `df` to a plain data frame at the top of both extraction helpers so the documented tibble / data.table inputs work (or, alternatively, drop the data.table claim from the dina() roxygen and reject data.table input with the package's own message). | add-only |
| A-017 | 1 | A | dina_subgroup.R | 425 | dina_subgroup | S2 | on | reproduced | dina_subgroup() and dina_frontier() check only that the fit has length(covariates) + 1 coefficients and then pair coefficients with columns by position, so a `covariates` vector in a different order (or naming different columns) than the fit was trained on silently gives a wrong tau-hat and a wrong subgroup. Also at: dina_subgroup.R:435 (positional tau_hat); dina_subgroup.R:1236 and dina_subgroup.R:1240 (same length-only check and positional product in dina_frontier); ... (full list in the CSV) | no | Compare names(coef(fit))[-1] with `covariates` in both functions and stop (or reorder by name) on a mismatch. | add-only |
| A-018 | 1 | A | dina_subgroup.R | 583 | dina_subgroup | S2 | on | reproduced | The returned `labels`/`label` render the threshold with paste0() on the raw double (15 significant digits), so evaluating the documented label with get_dfpred() can drop the observations sitting exactly on the threshold and disagree with the returned `mask` and `n_subgroup`. Also at: dina_subgroup.R:584-587 (`label`, same rendering); dina_subgroup.R:249-251 (roxygen: labels are 'in the form consumed by get_dfpred()'); ... (full list in the CSV) | no | Render the threshold in labels/label with a round-trip format (17 significant digits), as .forestsearch_dina_select() already does. | changes behaviour |
| A-019 | 1 | A | dina_subgroup.R | 738 | .dina_collect_candidates | S2 | on | reproduced | The depth-1 enumeration includes the non-cuts `x <= max(x)` and `x >= min(x)` (the whole sample, empty complement) for every covariate, so dina_subgroup() can return the entire data set as the selected 'subgroup' and dina_frontier() carries full-sample 'cuts' on every wide frontier. Also at: dina_subgroup.R:743-749 (admitted on size and floor only); dina_subgroup.R:536 (maxSG ranks size first); ... (full list in the CSV) | no | Decide explicitly how a candidate whose mask covers the whole sample (empty complement) is to be treated by the collectors; at present it is admitted without any check or flag. | method-adjacent |
| A-020 | 1 | A | dina_subgroup.R | 1190 | dina_frontier | S2 | on | reproduced | The harm orientation (tau_sign) that forestsearch applies to DINA for continuous/binary outcomes with adverse_outcome = FALSE exists only on the subgroup_method = "dina" path: dina_frontier() has no such argument and the use_dina screening call to dina_subgroup() omits it, so on those outcomes the DINA screening cuts are ranked and floored in the opposite (benefit) direction, contrary to the documented 'same cut subgroup_method = "dina" would select'. Also at: dina_subgroup.R:1240 (tau_hat without orientation); forestsearch_main.R:3217 (use_dina screening calls dina_frontier with da$frontier only); ... (full list in the CSV) | no | Give dina_frontier() the same tau_sign orientation as dina_subgroup() and pass .dina_tau_sign() from all forestsearch call sites, including the use_dina screening call to dina_subgroup(). | changes behaviour |
| A-021 | 1 | A | forestsearch_helpers.R | 122 | get_dfpred | S2 | on | reproduced | The exported get_dfpred() called with sg.harm = NULL or character(0), which is what a no-subgroup fit carries, silently flags every row as a harm-subgroup member (treat.recommend = 0) instead of raising an error or returning nobody. Also at: forestsearch_helpers.R:137 | no | Reject a NULL or zero-length sg.harm with an informative error (or return an all-complement flag) instead of the vacuous all-member result. | changes behaviour |
| A-022 | 1 | A | forestsearch_helpers.R | 202 | evaluate_comparison | S2 | on | reproduced | evaluate_comparison() coerces numeric-level factors but not numeric-valued character columns, so a cut label on such a column is compared as strings, unlike the DINA and GRF coding that selected the cut on the numeric values and whose comments say the coercion is the same. Also at: forestsearch_helpers.R:1344-1347, 1373-1374, 1510-1512 (roxygen / comments stating the coercion is identical to evaluate_comparison()) | no | Apply the same all-numeric coercion to character columns as to factors in evaluate_comparison(), so the label evaluator matches the DINA / GRF coding. | add-only |
| A-023 | 1 | A | forestsearch_helpers.R | 1311 | .dina_reselect_on_effect | S2 | on | reproduced | Effect-based re-selection overwrites only some winner fields, so the stored identifier object (grp.consistency$out_sg / fs$grf_res) mixes the new winner's definition with the native winner's label, SE, CI (DINA) or sg.harm.id, data$treat.recommend, grf.gsub, selected_depth/best_depth, cate_sg (GRF). Also at: forestsearch_helpers.R:1793-1794 (.grf_reselect_on_effect replaces only sg_def); forestsearch_helpers.R:1953-1957 (.forestsearch_grf_select details line prints the native winner's depth) | no | After re-selection, recompute or clear every winner-specific field of the stored identifier object (DINA: labels, label, se_mean_tau_hat, ci; GRF: sg.harm.id, data$treat.recommend, grf.gsub, selected_depth/best_depth, cate_sg). | add-only |
| A-024 | 1 | A | forestsearch_helpers.R | 1531 | .forestsearch_dina_select | S2 | on | reproduced | For count outcomes forestsearch() requires offset.name for IRR, but the DINA fit is called without exposure and dina()/dina_fit() have no offset argument, so DINA's tau-hat and the m_diff = log(threshold) proposal floor are on the mean-count ratio while forestsearch's text calls that floor log(IRR). Also at: forestsearch_main.R:3200-3206 (use_dina screening fit, same call shape); dina_df.R:84-90 and dina.R:157-167 (no offset/exposure formal); dina.R:537-544 and 842-845 (Poisson fits of the raw count); ... (full list in the CSV) | unknown | Either carry the exposure offset through dina()/dina_fit() for the Poisson family, or state in the documentation and at run time that DINA's proposal ignores offset.name and correct the 'log(IRR)' wording. | method-adjacent |
| A-025 | 1 | A | forestsearch_helpers.R | 1759 | .grf_reselect_on_effect | S2 | on | reproduced | When effect-based re-selection scores no candidate or errors, the native-statistic winner is returned silently with no effect floor applied; this happens deterministically for subgroup_method = "grf" with effect_measure = "IRD", because the per-candidate estimator returns NULL for IRD. Also at: forestsearch_helpers.R:1252 (DINA nothing-scorable return); forestsearch_helpers.R:1620 and 1924 (error = function(e) sg / grf_res); ... (full list in the CSV) | no | Raise a warning (an error under mr_inference) when re-selection fails or scores no candidate, record on the result which statistic actually selected, and refuse or handle estimands the per-candidate scorer cannot represent (IRD) instead of falling back. | changes behaviour |
| C-001 | 1 | C | forestsearch_main.R | 382 | forestsearch | S2 | on | reproduced | forestsearch() documents est.scale as "hr" or "rmst", never reads or validates it, and no code implements "rmst" (every consumer of the stored value knows only "hr" and "1/hr"), so the documented "rmst" is silently treated as "hr" by the search, tables and cross-validation and makes the full bootstrap run all replicates and then return no estimates. Also at: forestsearch_main.R:1637 (formal est.scale = "hr", never read in the body), 1880-1881 (captured unvalidated into args_call_all); man/forestsearch.Rd:150; bootstrap_dofuture_main.R:267, 559, 568, 574-616, 618-639; ... (full list in the CSV) | no | Correct the @param to the values the code supports ("hr", "1/hr") and validate est.scale at entry to forestsearch() so any other value fails before the search. | add-only |
| INV-3 | 1 | INV | forestsearch_main.R | 928 | forestsearch | S2 | on | reproduced | violated: roxygen states a default that the formal (or the code's effective default) contradicts at 12 confirmed sites in 11 functions, headed by forestsearch(use_twostage) documented FALSE with formal TRUE; three nested list-element defaults and one cross-reference are also wrong, and no reader claim was refuted at source. Also at: Cox_estimation_helpers.R:8; calibrate_glm_interaction.R:48; truefind_asymptotic_glm.R:180; cox_spline_fit.R:537; dina_subgroup.R:168; fpr_approximation.R:57; gg_forest.R:92; simulation_tables.R:1043; ... (full list in the CSV) | no | Correct the roxygen text at each listed site to the formal or effective default (or, where the documented value was the intended one, decide per site whether the formal should change instead). | add-only |
| A-026 | 1 | A | forestsearch_main.R | 1711 | forestsearch | S2 | on | reproduced | The documented tune_grf argument of forestsearch() is never read in the function body or forwarded to any GRF fit, so tune_grf = TRUE silently changes nothing on both the use_grf screening path and the subgroup_method = "grf" path. Also at: forestsearch_main.R:1077-1083 (roxygen claim), 3078-3092 and 3095-3112 (use_grf call sites), 2945-2965 (subgroup_method = "grf"); grf_args.R:25-44, 86-102 (builders have no tune_grf formal); ... (full list in the CSV) | no | Forward tune_grf through the two GRF argument builders and .forestsearch_grf_select() to the GRF fits, or remove the formal and its documentation. | changes behaviour |
| A-027 | 1 | A | forestsearch_main.R | 2220 | forestsearch | S2 | on | reproduced | The overdispersion argument is forwarded only to the GRF routines and never reaches the effect estimator, so on the default subgroup_method = "consistency" path with a count outcome and use_grf = FALSE, "quasi" and "negbin" change nothing. Also at: forestsearch_main.R:1709 (formal), 1067-1071 (roxygen), 2704-2711 (PS rebuild, also without it), 2952 and 3110 (the only consumers: GRF paths); ... (full list in the CSV) | no | Either carry overdispersion into the count effect estimator or document (and warn) that it applies to the GRF paths only. | changes behaviour |
| A-028 | 1 | A | forestsearch_main.R | 2291 | forestsearch | S2 | on | reproduced | A user-supplied effect or consistency threshold above 1 is replaced by the package default for effect_measure = "IRD" as well as "RD" on the ground that both are bounded in [-1, 1], which is not true of an incidence-rate difference, and the warning names 0.05 as the applied value when 0.01 is applied for IRD. Also at: forestsearch_main.R:296-310 (validation warning 'A risk difference cannot exceed 1.0' for IRD), 2303 (message hard-codes 0.05 while 2305-2306 applies 0.01 for IRD), 2312-2324 (consistency threshold > 1 reset to 0.0), 127; ... (full list in the CSV) | no | Restrict the greater-than-1 remap to RD, leave supplied IRD thresholds as given, and make the warning report the value actually applied. | changes behaviour |
| A-029 | 1 | A | forestsearch_main.R | 2656 | forestsearch | S2 | on | reproduced | A user-supplied ps_hat is silently discarded, unvalidated, whenever ps_method resolves to "none", which is the default for is.RCT = TRUE, so a call supplying propensity scores with ps_adjust_method = "iptw" runs unadjusted with no message. Also at: forestsearch_main.R:2596-2599 (ps_method resolves to "none" for is.RCT = TRUE), 2639-2641 (ps_adjust_method forced to "none"), 1131-1142 (ps_hat roxygen) | no | Warn or stop when ps_hat is supplied (or ps_adjust_method is not "none") while ps_method resolves to "none". | add-only |
| A-030 | 1 | A | forestsearch_main.R | 2705 | forestsearch | S2 | on | reproduced | When ps_method resolves to anything other than "none" on a GLM outcome, the estimator closure is rebuilt without adverse_outcome and reverts to make_effect_estimator()'s default TRUE, so the candidate search screens and ranks on the un-flipped outcome whenever the resolved adverse_outcome is FALSE (the continuous default, or binary with adverse_outcome = FALSE). Also at: forestsearch_main.R:2227 (first build passes adverse_outcome), forestsearch_main.R:3623 and :3841 (resample-consistency spec and MR spec carry the resolved adverse_outcome), glm_effect_estimators.R:107 (formal default ad; ... (full list in the CSV) | no | Pass the resolved adverse_outcome to the propensity-branch rebuild of the estimator closure so it matches the first construction. | changes behaviour |
| A-031 | 1 | A | forestsearch_main.R | 2875 | forestsearch | S2 | on | reproduced | On the subgroup_method = "dina" and "grf" branches mr_harm_confirmed is isTRUE(harm_flag), so an MR result whose harm_flag is NA (selected subgroup could not be fitted) is reported as FALSE, whereas the consistency branch and the documentation return NA for that case. Also at: forestsearch_main.R:3060-3061 (GRF branch, same expression); contrast forestsearch_main.R:3956-3957 (consistency branch keeps NA); forestsearch_main.R:3919-3927 (no .mr_skip message for the short variant on any branch); ... (full list in the CSV) | no | Return harm_flag unchanged on the DINA and GRF branches as the consistency branch does, and route the short variant through the existing not-performed message on all three branches. | changes behaviour |
| A-032 | 1 | A | forestsearch_main.R | 3071 | forestsearch | S2 | on | reproduced | With use_grf = TRUE a supplied grf_res that already has tree.cuts is never read for its cuts (no GRF cut enters the candidate family, with no message), and a supplied grf_cuts without grf_res is overwritten by a fresh GRF fit, so the documented reuse works only when both are supplied together. Also at: forestsearch_main.R:3124 and :3164 (user-supplied grf_cuts overwritten after a refit); forestsearch_main.R:3316 (grf_cuts passed on); get_fsdata.R:252 (use_grf switched off silently when grf_cuts is NULL); ... (full list in the CSV) | no | Take grf_cuts from a supplied grf_res when grf_cuts is NULL, and skip the refit when grf_cuts is supplied. | changes behaviour |
| A-033 | 1 | A | forestsearch_main.R | 3350 | forestsearch | S2 | on | reproduced | df.analysis is never coerced to a base data.frame, and Y, Event, Treat, id and the dummy matrix are built with df[, name], which does not drop to a vector for a tibble (and is not column selection by a character variable for a data.table), so such inputs are expected to fail with an internal error or a 'no subgroup' early return instead of being coerced or rejected at entry. Also at: forestsearch_main.R:3352, :3358, :3390 (single-column `[` relied on to drop to a vector); :3386-3387 with get_FSdata_helpers.R:747 (as.factor(df[, i])); :3391 (data.frame(Y, Event, Treat, id, ...)); ... (full list in the CSV) | no | Coerce df.analysis (and df.predict / df.test) to a base data.frame at entry, or reject non-base classes there with the message get_FSdata() already carries. | add-only |
| A-034 | 1 | A | forestsearch_main.R | 3391 | forestsearch | S2 | on | reproduced | With propensity adjustment requested on a GLM outcome (ps_method not "none" and ps_adjust_method "iptw" or "dr_gcomp") the weights are used only in the candidate search: the consistency frame df.fs is built without sw / ips_covar / ps_hat, so the consistency stage (split and resample) and MR silently fit unadjusted models, and the documented automatic fallback to "split" for propensity-adjusted effects does not exist. Also at: forestsearch_main.R:3398-3409 (only outcome / treat / offset re-attached to df.fs); :3616-3625 (glm_resample_spec has no PS field); :3837-3841 (MR spec has no PS field); ... (full list in the CSV) | no | Carry the propensity columns onto df.fs and either route propensity-adjusted fits to a consistency/MR path that uses them or refuse the combination, and correct the fallback claim in the roxygen. | changes behaviour |
| A-035 | 1 | A | forestsearch_main.R | 3399 | forestsearch | S2 | on | reproduced | On GLM outcomes the user's outcome/treatment/offset column is attached to the consistency frame only if its name is not already taken, so a GLM outcome column literally named "Event" is silently replaced by the all-ones placeholder and the consistency stage (and GLM MR) evaluate a constant outcome while the search scored the real one. Also at: forestsearch_main.R:3356 (GLM placeholder Event <- rep(1L, ...)), :3391 (df.fs built with internal Y/Event/Treat/id), :3402, :3405-3407 (same not-already-present guard for treat.name and offset.name), :3419-3425 (reserve; ... (full list in the CSV) | no | Stop (as is already done for adjust_covariates) when outcome.name, treat.name or offset.name names a different internal column of the consistency frame, or build that frame without reusing user-visible names. | add-only |
| A-036 | 1 | A | forestsearch_main.R | 3545 | forestsearch | S2 | on | reproduced | Under sg_focus = "maxeff" (documented as the effect argmax with no effect or consistency threshold, and recorded with an unrestricted admission set) the consistency-stage entry gate still requires some candidate effect to exceed c2, so when the largest candidate effect is <= c2 the stage is skipped and no subgroup is returned instead of the argmax. Also at: forestsearch_main.R:3540-3544 (check_threshold), forestsearch_main.R:3552 and :3745 (stage skipped, no else branch), forestsearch_main.R:134-135 (c2 derived as 0.80 * c1 for every focus, so hr.threshold alone moves the g; ... (full list in the CSV) | no | Either make the Section 8 gate read the resolved admission decision (.admit_applies) so the c2 comparison is not applied when no floor applies, or keep the gate and state it in the maxeff documentation, the override warning and the recorded admission set. | changes behaviour |
| A-037 | 1 | A | forestsearch_main.R | 3568 | forestsearch | S2 | on | reproduced | forestsearch() never passes its seedit to subgroup.consistency(seed = ), so the random consistency splits are always seeded with the callee's own default 8316951 and changing seedit has no effect on them. Also at: forestsearch_main.R:3629-3633 (filter_call_args call), :381 (roxygen 'Random seed'), :1636 (formal seedit); ... (full list in the CSV) | no | Pass seedit as the seed of the subgroup.consistency() call (or document that the split seed is fixed and independent of seedit). | changes behaviour |
| A-038 | 1 | A | forestsearch_main.R | 3722 | forestsearch | S2 | on | reproduced | merge() is called with its default sort = TRUE, so on the consistency path df.est and df.predict come back ordered by the id key while grp.consistency$sg.harm.id and df_flag stay in input row order, and positional consumers misalign whenever df.analysis is not already in ascending id order. Also at: forestsearch_main.R:3733; ... (full list in the CSV) | no | Attach the flag positionally (or merge with sort = FALSE and restore input order) so df.est and df.predict rows stay aligned with df.analysis and sg.harm.id on every path. | changes behaviour |
| A-039 | 1 | A | forestsearch_main.R | 3830 | forestsearch | S2 | on | reproduced | On the survival MR path adjust_covariates is intersected with the column names of the scoring frame, so any adjustment term that is not a bare column name (the documented "strata(x1)", or pspline()) is silently dropped and MR is computed from an unadjusted or partially adjusted Cox model while search and consistency used the adjusted one. Also at: fs_mr_inference_methods.R:207-208 (.fs_mr_spec: the same intersect() for the DINA and GRF MR branches, called at forestsearch_main.R:2845 and :3033); forestsearch_main.R:1088-1094 (roxygen documenting "strata(x1)") | no | Pass the adjustment terms to the MR spec exactly as the identifier used them (checking the underlying variables rather than the term strings), or stop/skip MR with a message when a term cannot be honoured. | changes behaviour |
| A-040 | 1 | A | forestsearch_main.R | 3844 | forestsearch | S2 | on | reproduced | On the consistency branch (the default subgroup_method) forestsearch() calls fs_mr_inference() without field_R_out / field_R_in, so those two mr_inference_args entries are silently ignored and the field block always runs at 1000 x 500 outer/inner draws. Also at: fs_mr_inference.R:597-598 (formal defaults 1000L / 500L), :919-920, :947, :985 (use and echo); fs_mr_inference_methods.R:178-179 (sibling DINA/GRF wrapper that does forward them); ... (full list in the CSV) | no | Forward mr_inference_args$field_R_out and $field_R_in in the consistency-branch call (ideally by routing that branch through the same wrapper as DINA/GRF so the forwarded sets cannot diverge). | changes behaviour |
| A-041 | 1 | A | fs_mr_inference_methods.R | 157 | .fs_apply_mr | S2 | on | reproduced | .fs_apply_mr() forwards 25 of fs_mr_inference()'s 29 arguments and drops keep_declaration_field, keep_field_matrix, declaration_c0 and pconsistency.digits, so on the DINA and GRF branches those mr_inference_args entries are silently ignored although the roxygen and NEWS say the full set is forwarded. Also at: fs_mr_inference_methods.R:120 (roxygen: 'forwards the full fs_mr_inference() argument set'); NEWS.md:584-594 (same claim); forestsearch_main.R:2860, :3045 (DINA / GRF call sites); ... (full list in the CSV) | no | Forward the four remaining arguments from mr_inference_args (and pconsistency.digits from the caller), or state in the roxygen/NEWS that the declaration-field arguments are not supported on the DINA and GRF branches. | add-only |
| A-042 | 1 | A | get_FSdata_helpers.R | 337 | lasso_selection | S2 | on | reproduced | For binary, continuous and count outcomes lasso_selection() returns "(Intercept)" in `selected`, so get_FSdata()'s 'if any selected' guard is always true and, when LASSO keeps no covariate, every confounder is filtered out (error / NULL fit) whereas the survival path keeps the whole pool. Also at: get_fsdata.R:354 (guard that the intercept defeats); forestsearch_main.R:3342-3344 (count includes the intercept and is labelled 'Cox-LASSO' for every outcome type) | unknown | Drop the intercept from the selected/omitted sets returned by lasso_selection() and apply one empty-selection rule to all outcome types. | changes behaviour |
| A-043 | 1 | A | get_FSdata_helpers.R | 692 | collapse_redundant_cuts | S2 | on | reproduced | A factor confounder with at least cont.cutoff (default 4) observed values and numeric-parsable levels makes get_FSdata() stop under the default collapse_cuts = TRUE, because collapse_redundant_cuts() calls stats::sd() on the raw factor column, so forestsearch() returns its early-exit list (sg.harm = NULL, error_log stage 'get_FSdata') with a base-R message about var() on a factor. Also at: get_FSdata_helpers.R:703-705 (membership safety check applies <= / >= to the raw factor; unreachable today because sd() stops first, but it would give NA and zero flips if only sd() were patched); ... (full list in the CSV) | no | Coerce numeric-level factor columns to their numeric values (as process_conf_force_expr() and evaluate_comparison() already do) before sd(), the membership comparisons of the safety check and the median computation, or reject such columns at entry with a clear message. | add-only |
| A-044 | 1 | A | get_FSdata_helpers.R | 859 | evaluate_cuts_once | S2 | on | reproduced | Missing confounder values are counted as a distinct level both when a variable is classed continuous and when a cut is tested for having more than one level, so a cut that is constant on every observed row is kept and returned as a one-level candidate factor. Also at: get_FSdata_helpers.R:11 | no | Count distinct non-missing values in both tests and decide explicitly at entry (error, or a documented rule) how missing confounder values are handled. | changes behaviour |
| A-045 | 1 | A | get_fsdata.R | 291 | get_FSdata | S2 | on | reproduced | Supplying conf.cont_medians or conf.cont_medians_force silently removes every cut on all other continuous confounders, and under the default cut_type the variables named in conf.cont_medians receive the four default cuts instead of a median cut, both contrary to the roxygen. Also at: get_fsdata.R:399-414; get_fsdata.R:31-32; forestsearch_main.R:581-588 | no | Make the code and the documentation agree on whether these arguments add median cuts to the default grid or restrict which continuous confounders are cut, and state the chosen rule in the roxygen of both functions. | changes behaviour |
| A-046 | 1 | A | get_fsdata.R | 486 | get_FSdata | S2 | on | reproduced | Every confounder with fewer than cont.cutoff (default 4) unique values enters the candidate pool as its bare name, which is evaluated as `x == 1` with no check that the variable is a 0/1 indicator, so a categorical confounder with no level equal to 1 (factor F/M, numeric 0/2) is removed from the search with no warning and a 3-level one contributes only 'level 1 versus the rest' under a label that reads as a binary indicator. Also at: get_fsdata.R:185-195 (type check accepts any numeric or factor); get_fsdata.R:272-279,286-287 (fewer than cont.cutoff unique values -> categorical, bare name); ... (full list in the CSV) | no | Where categorical confounders enter the pool, either expand each one to one indicator per level or stop (or warn unconditionally) when a categorical confounder is not a 0/1 indicator, and state the rule in the confounders.name documentation and in the methodology article. | changes behaviour |
| A-047 | 1 | A | get_fsdata.R | 499 | get_FSdata | S2 | on | reproduced | Existing candidates are compared with GRF and DINA cuts by unanchored grepl() regex matching, so a GRF or DINA cut is discarded as already represented whenever any categorical confounder name (or median-cut string) is merely a substring of it. Also at: get_fsdata.R:520; get_fsdata.R:451 | unknown | Compare on the parsed variable name (and operator and value where a full-expression match is intended) with exact or word-anchored matching instead of substring regex. | changes behaviour |
| D-001 | 1 | D | get_fsdata.R | 517 | get_FSdata | S2 | on | reproduced | Caller-supplied dina_cuts enter the candidate pool whenever non-empty, with no use_dina gate, unlike grf_cuts (inert unless use_grf = TRUE), while the family status, the family enumerators and the bootstrap/CV replays all treat the run as having no DINA cuts. Also at: forestsearch_main.R:3317 | no | Either gate the DINA merge on use_dina the way the GRF merge is gated on use_grf, or document that dina_cuts alone activates them and make the family status, enumerators and replays carry the supplied cuts consistently. | changes behaviour |
| A-048 | 1 | A | get_fsdata.R | 527 | get_FSdata | S2 | on | reproduced | conf_force expressions are not validated: one whose right-hand side is not a number (compound expression, unsupported function, typo) is compared as a character string and silently becomes a mislabelled candidate or is silently dropped. Also at: forestsearch_helpers.R:207-208; get_FSdata_helpers.R:414-417; get_FSdata_helpers.R:435-438 | no | Validate each conf_force element as a single 'variable operator number' comparison (or a supported deferred form) and stop on anything else, and refuse ordered comparisons whose right-hand side is not numeric. | add-only |
| A-049 | 1 | A | get_fsdata.R | 535 | get_FSdata | S2 | on | reproduced | exclude_cuts is applied only before collapse_redundant_cuts(), so the default coarsening (collapse_cuts = TRUE, digits = 0) can rewrite a surviving cut into exactly the excluded expression, which is then returned and searched. Also at: get_fsdata.R:550-572; get_FSdata_helpers.R:702-713 | no | Re-apply the exclude_cuts filter after the coarsening step (keeping the existing first pass so excluded thresholds do not enter cluster centroids). | add-only |
| A-050 | 1 | A | get_fsdata.R | 634 | get_FSdata | S2 | on | reproduced | Cuts that evaluate to a single value are dropped silently (the only message is under details = TRUE, which forestsearch() hard-codes to FALSE) and the zero-candidate diagnostic is tested only before that filter, so get_FSdata() can return an empty candidate set and forestsearch() then stops in dummy() with an unrelated message. Also at: get_fsdata.R:574; get_fsdata.R:648; get_FSdata_helpers.R:428-434; forestsearch_main.R:3318; forestsearch_main.R:3386-3387; get_FSdata_helpers.R:776 | no | Repeat the zero-candidate diagnostic after the validity filter and emit a warning (independent of details) naming the cuts or confounders dropped for having one level. | add-only |
| A-051 | 1 | A | get_fsdata.R | 678 | get_FSdata | S2 | on | reproduced | Candidate factor columns are written as q1, q2, ... with no collision check, so an input column carrying one of those names (for example a covariate called q1) is silently overwritten by a 0/1 factor in the returned data frame. Also at: get_fsdata.R:666-667; get_fsdata.R:685; get_fsdata.R:700; get_fsdata.R:708; bootstrap_analysis_dofuture.R:393-396 | no | Stop with a clear message when a generated factor name collides with a column named in confounders.name or the outcome/event columns (or generate names that cannot collide). | add-only |
| A-052 | 1 | A | glm_effect_estimators.R | 163 | make_effect_estimator | S2 | on | reproduced | For a binary outcome with a rate measure (IRR) and adverse_outcome = FALSE, the search estimator fits the rate ratio of Y (the flag is not forwarded to the Poisson rate estimator) while the default resample consistency fit flips Y to 1 - Y first, so the consistency rate describes a different quantity from the one the search ranked. Also at: consistency_resample.R:255-260; forestsearch_main.R:3616-3624 | no | Apply one convention on both sides for binary rate measures: either reject adverse_outcome = FALSE with IRR/IRD or make the resample fit skip the flip exactly where the estimator does. | changes behaviour |
| A-053 | 1 | A | glm_effect_estimators.R | 187 | make_effect_estimator | S2 | on | reproduced | effect_measure is validated against the outcome type for binary and count outcomes but not for continuous, so a ratio or risk measure supplied with outcome_type = 'continuous' is accepted: the estimate is a mean difference while forestsearch() log-transforms the thresholds, labels the effect with the supplied measure and exponentiates it for display. Also at: forestsearch_main.R:1789-1795; forestsearch_main.R:2340-2368; forestsearch_main.R:3607-3609; forestsearch_main.R:3616-3625; forestsearch_main.R:1789 | no | Validate effect_measure in the continuous branch (only 'MD') and stop on a mismatch, in make_effect_estimator() and at forestsearch() entry. | add-only |
| A-054 | 1 | A | glm_effect_estimators.R | 429 | .make_glm_binary_estimator | S2 | on | reproduced | A requested ps_adjust_method has no implementation for several measure/method pairs (RR with iptw or dr_gcomp, RD with dr_gcomp, MD with dr_gcomp), yet forestsearch() accepts every method for every GLM measure, fits those cases unadjusted and still reports the requested method; the literal-split consistency fits and the bootstrap/table estimators are also unadjusted because their data or closure carry no propensity columns. Also at: glm_effect_estimators.R:414-419; glm_effect_estimators.R:465; glm_effect_estimators.R:952-966; forestsearch_main.R:2639-2641; forestsearch_main.R:2705-2712; forestsearch_main.R:3391-3409; ... (full list in the CSV) | no | Reject at forestsearch() entry (or in the factory) a ps_adjust_method that has no implementation for the resolved effect measure, and either carry the propensity columns onto the consistency frame or state that split consistency is unadjusted. | changes behaviour |
| A-055 | 1 | A | glm_effect_estimators.R | 535 | .estimate_rd | S2 | on | reproduced | With effect_measure = "RD" and no covariates, a slice with an arm at 0% or 100% events gets its standard error from the logistic tier as (logit treatment-coefficient SE) x mean(p1(1-p1)), which is near zero when the treated arm is degenerate and in the hundreds to thousands when the control arm is, while the point estimate is correct and converged = TRUE is returned. Also at: glm_effect_estimators.R:536 (se_rd), :541 (converged = TRUE), :545 (method_used), :548 (error-only handler, so tier 3 at :553-556 is skipped); ... (full list in the CSV) | unknown | On unadjusted slices with a degenerate arm, return a standard error that is valid at the boundary (the raw-proportion form already written in tier 3) instead of accepting the quasi-separated logistic fit in tier 2. | changes behaviour |
| A-056 | 1 | A | glm_effect_estimators.R | 671 | .estimate_rr | S2 | on | reproduced | With the non-default effect_measure = "RR", the log-binomial attempt in .estimate_rr() handles warnings in an exiting tryCatch() handler that calls invokeRestart("muffleWarning"), so a glm "step size truncated" warning becomes the uncaught error "no 'restart' 'muffleWarning' found" and the modified-Poisson fallback at :680-705 is never reached. Also at: glm_effect_estimators.R:505-510 (.estimate_rd tier 1, identical handler; ... (full list in the CSV) | no | Handle non-"converge" warnings from the log-binomial attempt (and the identity-link tier in .estimate_rd) without invoking a restart from an exiting handler, either by muffling them in a calling handler or by routing them to the next fallback tier. | changes behaviour |
| A-057 | 1 | A | glm_effect_estimators.R | 776 | .make_poisson_rate_estimator | S2 | on | reproduced | The IPTW branch of the Poisson rate estimator never looks at effect_measure, so with effect_measure = "IRD" and ps_adjust_method = "iptw" it returns the log incidence rate ratio while every caller treats the value as a rate difference. Also at: glm_effect_estimators.R:786-794 (returns coefficient, measure = "log-IRR"); glm_effect_estimators.R:749 (IRR existence check skipped because effect_measure is IRD); ... (full list in the CSV) | no | Make the IPTW branch honour effect_measure = "IRD" (compute the weighted rate difference as the other branches do) or reject the IRD + iptw combination with an error. | changes behaviour |
| A-058 | 1 | A | grf_helpers.R | 684 | .grf_code_column | S2 | on | reproduced | .grf_code_column() derives the coding of a character column from the values present in the frame it is handed (unique(x), NA included), so a GRF cut can be evaluated on a different numeric scale on the unfiltered df, df.predict or df.test than the scale the forest split on. Also at: grf_helpers.R:690 (codes from the values present); grf_subgroup_labels.R:225 (evaluator); forestsearch_helpers.R:1987, 1990, 1993 (evaluator applied to df, df.predict, df.test); ... (full list in the CSV) | no | Derive the coding from a level set fixed at fit time and carried with the subgroup definition, and ignore NA when testing whether the labels are all numeric. | changes behaviour |
| A-059 | 1 | A | grf_main.R | 226 | grf.subg.harm.survival | S2 | on | reproduced | The survival GRF path builds its covariate matrix with apply(data[, confounders], 2, as.numeric), which passes through a character matrix whenever any confounder is a factor or character: non-numeric labels become NA, numeric columns are rounded to 7 significant digits, and the coding differs from the GLM path and from the subgroup evaluator. Also at: grf_main.R:308 and forestsearch_helpers.R:1987 (membership re-derived on the raw columns); grf_helpers.R:427-436 (tree-path membership from predict(tree, X)); grf_subg_harm_glm.R:868-894 (sibling .build_grf_X coding) | unknown | Build the survival covariate matrix with the same per-column coding the GLM path and the evaluator use (.build_grf_X / .grf_code_column). | changes behaviour |
| A-060 | 1 | A | grf_main.R | 235 | grf.subg.harm.survival | S2 | on | reproduced | When validate_grf_data() rejects the data (one treatment level, one event level such as fully uncensored data, or an arm smaller than n.min), the early exit calls create_null_result() with trees = list(), which stops with 'subscript out of bounds' instead of returning the null result. Also at: grf_helpers.R:535 (trees[[1]] on an empty list), grf_helpers.R:536, grf_helpers.R:559; grf_helpers.R:535 | no | Pass NULL (as the frontier no-selection exit does) or make create_null_result() tolerate an empty tree list so the early exit returns the null result. | add-only |
| A-061 | 1 | A | grf_main.R | 295 | grf.subg.harm.survival | S2 | on | reproduced | Both GRF identifiers test the frontier selection with is.null() only, but .grf_frontier_select() returns a zero-row data frame when the band is empty (effMaxSG/effMinSG with a negative maximum eligible effect), and the next step stops with 'argument is of length zero'. Also at: grf_subg_harm_glm.R:528 and 547 (same is.null-only test); grf_subgroup_labels.R:387-391 (zero-row return), grf_subgroup_labels.R:405 (where it stops) | no | Treat a zero-row selection as no subgroup in both identifiers, as .grf_reselect_on_effect() already does. | add-only |
| A-062 | 1 | A | grf_main.R | 297 | grf.subg.harm.survival | S2 | on | reproduced | When no DR candidate clears dmin.grf, the frontier no-selection return omits the candidates table, so the effect-based re-selection (default grf_select_statistic = 'effect') is skipped and dmin.grf acts as an all-or-nothing gate on the DR scale, although it is not applied to individual candidates once any one clears it. Also at: grf_subg_harm_glm.R:532-545 (same omission on the GLM identifier); forestsearch_helpers.R:1910-1913 (consumer that skips the re-selection); forestsearch_helpers.R:1913 | unknown | Either attach the candidate family to the frontier no-selection returns of both identifiers so the effect floor is the only admission rule, or document that dmin.grf is applied first as an existence gate on the DR scale. | changes behaviour |
| A-063 | 1 | A | grf_subg_harm_glm.R | 303 | grf.subg.harm.glm | S2 | on | reproduced | return_selected_cuts_only is accepted by grf.subg.harm.glm() but never acted on (cuts always come from the single best-depth tree), and frac.tau and conf.level there, and id.name in grf.subg.harm.survival(), are likewise accepted and unused. Also at: grf_subg_harm_glm.R:297 (frac.tau), 301 (conf.level), 900 (n and n.min unused in .select_best_sg_glm); grf_main.R:155 (id.name unused); grf_helpers.R:35-36 (config fields nothing reads) | no | Honour return_selected_cuts_only on the GLM identifier as the survival one does, or document it as having no effect there, and remove or document the other unused arguments. | changes behaviour |
| A-064 | 1 | A | grf_subg_harm_glm.R | 443 | grf.subg.harm.glm | S2 | on | reproduced | For outcome_type = "count" the adverse_outcome orientation is never applied to the GRF outcome, so the GRF 'harm' subgroup is the region where treatment LOWERS the count, the opposite of the documented meaning of adverse_outcome = TRUE and of the package's own count convention (IRR > 1 = harm). Also at: grf_subg_harm_glm.R:451 (comment 'flip not needed'); grf_subg_harm_glm.R:100-103 (roxygen: adverse event count => harm = treatment increases outcome); ... (full list in the CSV) | yes | Orient the count GRF outcome so that a higher count is the worse outcome when adverse_outcome is TRUE (as already done for binary and continuous), and correct the comment at line 451. | changes behaviour |
| A-065 | 1 | A | grf_subg_harm_glm.R | 914 | .select_best_sg_glm | S2 | on | reproduced | With sg.criterion = 'Nsg' the GLM selector returns the largest subgroup without comparing its score with dmin.grf, whereas the 'mDiff' branch and the survival sibling both apply the threshold. | no | Apply the dmin.grf threshold in the 'Nsg' branch as select_best_subgroup() does. | changes behaviour |
| A-066 | 1 | A | subgroup_consistency_main.R | 403 | subgroup.consistency | S2 | on | reproduced | consistency_threshold = NULL is documented as falling back to hr.consistency, but on the split path it is never resolved: with estimator_fn supplied and the threshold NULL the splitter skips the GLM estimator and runs the survival branch, so a direct GLM caller silently gets no subgroup (or a Cox comparison if the frame happens to carry Y/Event/Treat). Also at: subgroup_consistency_main.R:209-210 (roxygen: "Default NULL (uses hr.consistency)"); subgroup_consistency_helpers.R:241-242 (same statement for run_single_consistency_split); ... (full list in the CSV) | no | Resolve a NULL consistency_threshold to hr.consistency when estimator_fn is supplied (as the roxygen states), or make the GLM branch depend on estimator_fn alone. | changes behaviour |
| A-067 | 1 | A | subgroup_consistency_main.R | 550 | subgroup.consistency | S2 | on | reproduced | Under sg_focus = "maxeff" a finite m1.threshold is accepted and silently not applied, because the branch that skips the effect floor also skips the m1 filter, and neither the maxeff override warning nor the documentation nor fs_family_report() says so. Also at: forestsearch_main.R:2017-2041 (maxeff override block: m1.threshold neither reset nor named in the warning); forestsearch_main.R:699-701 (maxeff roxygen lists what is disabled, m1.threshold not among them); ... (full list in the CSV) | no | Have forestsearch() announce (and fs_family_report() record) that m1.threshold is not applied under maxeff, or apply the m1 filter on the maxeff branch if that is the intent. | add-only |
| INV-5 | 1 | INV | subgroup_consistency_main.R | 717 | subgroup.consistency | S2 | on | reproduced | violated: the selection is plan-independent on the default path (consistency_method = "resample" with an estimable closed form, where no RNG is drawn), but whenever the consistency rate is computed by literal splitting - consistency_method = "split", effect_measure = "IRD", or the per-candidate fallback when the closed form is NA - plan = "sequential" draws every candidate's splits from one stream seeded once while plan = "multisession" gives each candidate its own L'Ecuyer-CMRG stream, so Pcons, admission and the selected subgroup can differ for the same seed. Also at: subgroup_consistency_main.R:498, :793, :970, :994; subgroup_consistency_helpers.R:272, :343, :1492, :1520, :1529, :1733, :1890; consistency_resample.R:76; forestsearch_main.R:3568; subgroup_consistency_main.R:717; ... (full list in the CSV) | no | Drive the plain sequential loop from the same per-candidate streams the batched path uses (or the reverse) so both plans draw identical splits, or else state the plan dependence of every literal-splitting path in forestsearch()'s own documentation and correct the 'identical on every path' wording for the IRD and NA-fallback cases. | changes behaviour |
| INV-2 | 1 | INV | subgroup_search.R | 174 | subgroup.search | S2 | on | reproduced | violated: on the default forestsearch() path subgroup.search() replaces the caller's future plan (multisession, then sequential; no capture, no on.exit) and subgroup.consistency() calls set.seed(8316951) with no restore of .Random.seed; further unrestored plan / RNG sites exist in the bootstrap, MR inference and Tier-2 simulation code. Also at: betaHhat_truth.R:849; bootstrap_analysis_dofuture.R:334; bootstrap_dofuture_main.R:352; forestsearch_cross_validation.R:369; fpr_calibration.R:272; frontier_cis.R:164; fs_fdr_report.R:398; fs_mr_inference.R:615; ... (full list in the CSV) | no | Capture the previous plan / .Random.seed (and kind) at entry of each listed function and restore them with on.exit, as subgroup.consistency()'s plan handling and .make_candidate_rng_seeds() already do. | changes behaviour |
| A-068 | 1 | A | subgroup_search.R | 370 | prepare_search_data | S2 | on | reproduced | With one candidate column the Z matrix collapses to a vector (no drop = FALSE) so the column count is NULL, and with fewer candidate columns than maxk combn() is called unguarded; both stop with the bare error "n < m" instead of evaluating the available candidates or giving an informative error. Also at: subgroup_search.R:132; subgroup_search.R:389; subgroup_search.R:439 | no | Keep Z two-dimensional when subsetting and either cap the enumerated combination size at the number of candidate columns or stop with an explicit message. | changes behaviour |
| A-069 | 1 | A | subgroup_search.R | 628 | evaluate_combination_with_status | S2 | on | reproduced | When subgroup.search() is given estimator_fn (or adjust_covariates) without a row-aligned df_analysis, it silently scores candidates with the unadjusted treatment-only Cox fit on (Y, Event) and returns those numbers in the usual columns instead of stopping. Also at: subgroup_search.R:138-147; subgroup_search.R:991-999 | no | Stop with an error when estimator_fn or adjust_covariates is supplied without a row-aligned df_analysis rather than falling through to the treatment-only Cox path. | changes behaviour |
| A-070 | 1 | A | subgroup_search.R | 771 | has_positive_variance | S2 | on | reproduced | subgroup.search() documents Z as 'Matrix or data frame' and its own example passes a data frame, but Z is never coerced, so a data-frame Z reaches `t(x) %*% x` and a direct call stops with a base matrix-product error on the first combination; only a numeric matrix works. Also at: subgroup_search.R:15 (roxygen '@param Z Matrix or data frame of candidate subgroup factors (binary indicators).'; same text in man/subgroup.search.Rd:41); ... (full list in the CSV) | no | Coerce Z to a numeric matrix on entry to subgroup.search() (or state matrix-only in the roxygen) and replace the example with one that passes 0/1 indicator columns. | add-only |
| A-071 | 2 | A | summary_utility_functions.R | 131 | format_results | S2 | on | reproduced | format_results() builds a table row with c() and includes the optional fields by per-row is.na() tests, so the exported SG_tab_estimates() fails with its own default NULL subgroup names, and rows of different widths arise when only one adjusted HR or one potential-outcome HR is NA. Also at: summary_utility_functions.R:128-141 (all four c() rows), :606 (NULL name defaults), :617-626 (rbind + column names keyed on hr_1a and potentialOutcome.name only), :518-522 (empty-arm early return always passes hr_po = NA; ... (full list in the CSV) | no | Return fixed-width rows (name coerced to a string, explicit NA placeholders) and decide the column set from the arguments rather than from per-row NA tests. | changes behaviour |
| A-072 | 2 | A | summary_utility_functions.R | 366 | SG_tab_estimates_glm | S2 | on | reproduced | For GLM outcomes est.scale = "1/hr" reverses only the descriptive arm columns (n1, Rate/Mean(C), Rate/Mean(T), Diff) while the effect estimate still comes from the estimator closure bound to the original treatment column, so a row shows reversed arms beside an un-reversed effect. Also at: summary_utility_functions.R:52-55 (prepare_subgroup_data adds treat2), :195-204 (arm counts and rates from the reversed column), :219-233 (effect from the closure on the original column); ... (full list in the CSV) | no | Apply one coding consistently in the GLM table: either invert the effect estimate along with the arm columns or leave the arm columns un-reversed. | changes behaviour |
| A-073 | 2 | A | summary_utility_functions.R | 537 | analyze_subgroup | S2 | on | reproduced | A subgroup that contains patients from one treatment arm only is not guarded: cox_summary() returns an NA string for it, but the following km_summary() and rmst_calculation() calls are unguarded two-arm computations, so the whole table is expected to abort instead of giving an NA row. Also at: summary_utility_functions.R:73-76 (km_summary), :540 and :430-431 (rmst_calculation), :518-522 (the guard that covers only the empty case) | no | Extend the degenerate-subgroup guard to the single-arm case and return an NA row of the normal width. | add-only |
| A-074 | 2 | A | summary_utility_functions.R | 629 | SG_tab_estimates | S2 | on | reproduced | Supplying the documented optional hr_1a/hr_0a arguments to sg_tables() (or calling SG_tab_estimates() with SG_flag = "ITT" and a non-NA hr_1a) always stops with a base-R error, because the ITT row is built without an adjusted-estimate field while its names, and the later rbind, are sized as if it had one. Also at: summary_utility_functions.R:630-636 (ITT names sized by hr_1a; error at 634, or 635 with a potential-outcome column); 1024-1037 (sg_tables survival path forwards hr_1a/hr_0a to the ITT call); ... (full list in the CSV) | no | Give the ITT row the same width as the subgroup rows when adjusted estimates are supplied (an empty adjusted cell), in both the survival and the GLM table builders. | add-only |
| A-075 | 2 | A | summary_utility_functions.R | 629 | SG_tab_estimates | S2 | on | reproduced | With est.scale = "1/hr" the ITT branch of SG_tab_estimates() never reads est.scale while the subgroup branch reverses the treatment arms, so the sg_tables() estimates table shows HR and original-arm n1/m1/m0 in the ITT row and 1/HR and swapped-arm n1/m1/m0 in the subgroup rows under the same column headers. Also at: summary_utility_functions.R:609-616 with 52-55 (subgroup branch reverses the arms); 1036-1043 (sg_tables stacks ITT and subgroup rows under one header); ... (full list in the CSV) | no | Apply the est.scale arm reversal to the ITT row as well, so every row of the estimates table is on one coding (and make the Table 2 source note follow its column). | changes behaviour |
| A-076 | 1 | A | forestsearch_main.R | 1920 | forestsearch | S2 | on | traced | A df.analysis that is not a data frame is routed into the zero-row early exit, so it is reported as 'df.analysis has zero rows' and returned as a no-subgroup result, and with quiet = TRUE nothing is signalled at all. Also at: forestsearch_main.R:1921-1924; get_fsdata.R:99-106 | no | Split the two conditions: stop with a type error for a non-data.frame input and keep the quiet early exit for a zero-row data frame only. | changes behaviour |
| A-077 | 1 | A | glm_effect_estimators.R | 529 | .estimate_rd | S2 | on | traced | On a slice with an empty treatment arm the risk-difference estimator returns estimate = 0 with se = NA and converged = TRUE instead of a non-estimable result, unlike MD (NA) and OR/RR (NA with a reason). Also at: glm_effect_estimators.R:399; glm_effect_estimators.R:533 | no | Return the non-estimable result for RD when either arm of the slice is empty. | changes behaviour |
| INV-5a | 1 | INV | subgroup_consistency_main.R | 897 | subgroup.consistency | S2 | on | traced | In a direct subgroup.consistency() call with a non-NULL stop_threshold and a Pcons-keyed focus ("hr", "maxSG", "minSG", "hrMaxSG", "hrMinSG"), the sequential loop stops at the first candidate reaching stop_threshold while a parallel plan evaluates and keeps the whole batch containing it, so the selected subgroup can differ between plans and between worker counts with no randomness involved (also under the default "resample" method). Also at: subgroup_consistency_main.R:853, :913, :1005, :1014; roxygen subgroup_consistency_main.R:156-161, :276-279 | no | Make subgroup.consistency() itself ignore (or reject with a warning) a stop_threshold for every focus whose key contains Pcons, as forestsearch() already does, and correct the roxygen that calls early stopping valid for those foci. | changes behaviour |
| A-078 | 1 | A | subgroup_consistency_main.R | 448 | subgroup.consistency | S2 | on | suspected | subgroup.consistency() neither documents nor checks that df must contain an id column; without one it returns a selected subgroup together with an all-zero sg.harm.id and a df_flag that has no id column and recommends treatment for every subject, with no error or warning. Also at: subgroup_consistency_main.R:75-76 (roxygen lists only Y, Event, Treat); subgroup_consistency_helpers.R:934-937; subgroup_consistency_helpers.R:939-942 | no | Require and document the id column in subgroup.consistency(), or derive membership from row position instead of id matching. | add-only |
| A-079 | 2 | A | ahr_analytic.R | 267 | .ahr_set_index | S2 | off | reproduced | A standard set name that has no patient in the population (for example "Hstar" on a DGM built without an x3 subgroup) is not rejected, and ahr_limit() with b and h given returns a finite hr (the uniroot bracket end, exp(-10)) with NaN ahr/cde and no error, while the same empty set passed as a logical vector stops with "The set is empty.". Also at: ahr_analytic.R:420 (.ahr_fast adds no part when n == 0); ahr_analytic.R:328 (.ahr_core: uniroot on an identically zero score); ahr_analytic.R:621-622 (ahr_design computes requested standard sets the same way); ... (full list in the CSV) | no | Apply the existing empty-set stop to named sets as well, before either route is entered. | add-only |
| A-080 | 2 | A | ahr_analytic.R | 622 | ahr_design | S2 | off | reproduced | ahr_design() does not validate b and h: with b = NULL (alone or with h = NULL) every set is computed on its own direct-route grid from the DGM's stored effects, a supplied h is ignored yet echoed in the result, and the nested-set correlations and trigger probabilities are formed from node arrays that do not share a grid. Also at: ahr_analytic.R:612-618 (validation block covers pop, n, sets only); ahr_analytic.R:456-458 (.ahr_exact: `fast <- two && ...`, so force_fast is ignored when b or h is NULL); ... (full list in the CSV) | no | Reject NULL or non-scalar b and h at the top of ahr_design(), as ahr_design_table() already does. | add-only |
| A-081 | 2 | A | bootstrap_summaries_helpers.R | 722 | create_bootstrap_diagnostic_plots | S2 | off | reproduced | For a fit with est.scale = "1/hr" the 'observed' and 'bias-corrected' reference lines of the bootstrap diagnostic histograms are drawn at the negative of their log-scale positions, because log() is applied to estimates already inverted to 1/HR while the histogram is of uninverted log estimates. Also at: bootstrap_summaries_helpers.R:752-753 (Hc plot); summarize_bootstrap_results.R:304-310 (call site that has est.scale but does not pass it) | no | Pass the object's est.scale to the plot helper and place the reference lines on the same uninverted log scale as the histogram (or invert the histogram consistently). | changes behaviour |
| A-082 | 2 | A | bootstrap_summaries_helpers.R | 887 | summarize_factor_presence_robust | S2 | off | reproduced | The helper re-filters its input to rows with non-missing Pcons, so for DINA-selection bootstraps (Pcons always NA although any_found = 1) factor_presence and factor_presence_specific come back empty even though subgroups were identified. Also at: summarize_bootstrap_subgroups.R:107-108 (caller selects found rows by any_found), 469-474 (call), 628-637 (returned as factor_presence / factor_presence_specific) | no | Do not re-filter on Pcons inside the helper (the caller already passes the found rows), or apply the same any_found-first rule there. | changes behaviour |
| A-083 | 2 | A | calibrate_null_correction.R | 357 | run_null_calibration | S2 | off | reproduced | Every forestsearch() error inside the calibration loop is discarded and counted as 'no subgroup found', so a fs_args list that makes every call fail returns a fitted leff_calibration built only from the clamp value, with no warning or failure count. Also at: calibrate_null_correction.R:362-365 (NULL scored as not found); calibrate_null_correction.R:392-393 (zero FPR clamped to 0.5/n_sims) | no | Count errored replicates separately from non-detections and stop or warn (with the first condition message) when replicates fail. | changes behaviour |
| A-084 | 2 | A | compare_fs_subgroups.R | 123 | compare_fs_subgroups | S2 | off | reproduced | A GRF result in which no subgroup was found is rejected with an input-type error, although the roxygen says non-detection is handled automatically by classifying all subjects to the complement. Also at: compare_fs_subgroups.R:133-137 (the stop reached); grf_helpers.R:540-543 (null GRF result returns the input data without treat.recommend); roxygen compare_fs_subgroups.R:17-19, 58-61 | no | Treat a GRF result list whose data lacks treat.recommend (sg.harm.id NULL) as non-detection in the same way as a forestsearch object with no df.est, or correct the documentation. | add-only |
| A-085 | 2 | A | compare_selection_rules.R | 582 | compare_selection_rules | S2 | off | reproduced | Result slots are filled with `list[[k]] <- value` where value can be NULL, which deletes element k rather than storing NULL, so the returned $fs, $errors, $ci_tab and $plots lists come back shorter than n_combos with names shifted onto the wrong combinations. Also at: compare_selection_rules.R:585 (error_list[[k]] <- err_k, NULL on every successful fit); :603 (ci_list[[k]] <- ci_k, NULL when compute_frontier_cis errors); ... (full list in the CSV) | no | Store NULL results without deleting the slot (assign only when non-NULL, or assign through single-bracket with list(NULL)) for all four lists. | changes behaviour |
| A-086 | 2 | A | cox_ahr_cde_wrapper.R | 692 | cox_ahr_cde_analysis | S2 | off | reproduced | With plot_style = "grid", a threshold and CDE columns, the function stops in density() after all computation whenever one of the two subgroups has fewer than two subjects, which includes its own 'no subgroup found' fallback, so no result is returned. Also at: cox_ahr_cde_wrapper.R:693; cox_ahr_cde_wrapper.R:224-225 and :270-275 (NaN AHR/CDE for an empty subgroup, no message) | no | Skip or guard the subgroup density overlays when a subgroup has fewer than two subjects, and flag the empty-subgroup NaN summaries. | add-only |
| A-087 | 2 | A | cox_ahr_cde_wrapper.R | 849 | summary.cox_ahr_cde | S2 | off | reproduced | summary.cox_ahr_cde() applies `$knots` to the splines::ns() basis (an atomic matrix whose knots are an attribute), so the exported summary method stops with a `$`-on-atomic error for every object returned by cox_ahr_cde_analysis(). | no | Take the knot count from the basis attribute (or report the number of basis columns) instead of using `$knots` on the ns matrix. | changes behaviour |
| A-088 | 2 | A | create_summary_table.R | 116 | create_summary_table | S2 | off | reproduced | Records with a missing treatment value are not excluded before the arm split, so each one becomes an all-NA row in BOTH arm data frames, inflating the header n of both columns and producing spurious 'Missing' rows for every variable. Also at: create_summary_table.R:117 (second arm); consequences at 148/153, 211/214, 325-326, 333/336, 438-439, 446/449, 548-549, 556/559, 601-602, 874-881 | no | Drop (or stop on / report) records with missing treatment before splitting the data by arm. | changes behaviour |
| A-089 | 2 | A | create_summary_table.R | 116 | create_summary_table | S2 | off | reproduced | The first sorted level of treat_var is always labelled "Control" and the second "Treatment", with the level values never shown, so a two-level coding whose control level does not sort first (e.g. character "active"/"placebo") is displayed with the arms swapped. Also at: create_summary_table.R:107 (sort order decides the arm), 117, 123-124, 601-602, 873-881 (hard-coded 'Control'/'Treatment' labels); roxygen :15 | no | Either label the two columns with the actual level values or validate and document that treat_var must be coded 0 = control, 1 = treatment. | changes behaviour |
| A-090 | 2 | A | create_summary_table.R | 747 | create_summary_table | S2 | off | reproduced | With the default use_alternating_rows = TRUE the even-row index is built as seq(2, nrow, by = 2), which errors ('wrong sign in by argument') when the table has one row (or none), so a summary of a single continuous or binary variable with no missing values cannot be produced. | no | Build the even-row index with a construction that is empty when the table has fewer than two rows, before styling. | add-only |
| A-091 | 2 | A | cv_summary_tables.R | 601 | cv_metrics_tables | S2 | off | reproduced | With include_raw = TRUE and table_style = "combined" the is.list(result) test is TRUE for a data.frame and for a gt table, so the raw matrices are assigned into the table object instead of being wrapped in a list: with use_gt = FALSE this errors for a forestsearch_tenfold result and returns a data.frame with a junk list column for a forestsearch_Kfold result. Also at: cv_summary_tables.R:602-608 (assignments into the table object); 609-620 (list-wrapping branch unreachable for 'combined'); roxygen 311-312 | no | Choose between appending and wrapping by table_style (or by whether the result is a bare table object) rather than by is.list(result). | changes behaviour |
| A-092 | 2 | A | dina_bagged.R | 194 | dina_fit_bagged | S2 | off | reproduced | dina_fit_bagged(seed = ) seeds only the bootstrap indices; the per-bag cross-fitting folds (on by default) are drawn unseeded from the caller's RNG stream, so the same seed does not reproduce the fit and the documented 'global RNG state is restored on exit' does not hold. Also at: dina_bagged.R:357 (per-bag dina_fit(seed = NULL)); dina_bagged.R:204-207 (parallel branch: seed = TRUE streams from the ambient RNG state); dina_bagged.R:43 with dina.R:109-110 (inherited @param seed text; ... (full list in the CSV) | no | Derive the per-bag fold seeds from `seed` in both the sequential and the parallel branch, and give dina_fit_bagged its own @param seed text that states what the seed controls. | changes behaviour |
| A-093 | 2 | A | dina_bagged.R | 365 | .dina_subset_Y | S2 | off | reproduced | dina_fit_bagged() with family = "cox" and Y supplied as a two-column matrix (a form its inherited @param Y documents and dina_fit() accepts) subsets Y by linear indexing, so every bag fails and the call stops. Also at: dina_bagged.R:385 (call site in .dina_safe_one_bag); dina_bagged.R:149-158 (dina_fit_bagged validates/coerces X only, never Y); ... (full list in the CSV) | no | Coerce a matrix Y to a data frame once at the top of dina_fit_bagged() for the Cox family, as dina_fit() does (or make the subsetting helper row-subset matrices). | add-only |
| A-094 | 2 | A | dina_subgroup_bootstrap.R | 643 | dina_subgroup_bootstrap | S2 | off | reproduced | A supplied `sg` computed with tau_sign = -1 is accepted, but the bootstrap has no tau_sign: the fixed-subgroup effect CI is computed on the fit's own scale (opposite sign to sg$mean_tau_hat and sg$ci) and every resample applies the m_diff floor to the un-negated tau-hat, i.e. the opposite tail. Also at: dina_subgroup_bootstrap.R:805-809 (per-resample search with tau_sign = 1); dina_subgroup_bootstrap.R:412-459 (no tau_sign check on sg); dina_subgroup_bootstrap.R:940 vs 946-949 (printed side by side); ... (full list in the CSV) | no | Carry the orientation of a supplied sg (or a tau_sign argument) into the effect functional and every per-resample dina_subgroup() call, or reject an sg whose call records tau_sign = -1. | changes behaviour |
| A-095 | 2 | A | dina_subgroup_bootstrap.R | 805 | .dina_safe_one_bootstrap | S2 | off | reproduced | dina_subgroup_bootstrap() has no max_depth, grid_probs or selection_rule argument and does not compare them with a supplied `sg`, so an `sg` found with max_depth = 1 or selection_rule = "pareto"/"both" is accepted while every resample is searched with max_depth = 2 and the neighborhood rule. Also at: dina_subgroup_bootstrap.R:470-475, 483-488 (point-estimate searches, same omission); dina_subgroup_bootstrap.R:273-296 (formals); dina_subgroup_bootstrap.R:412-459 (sg validation block); ... (full list in the CSV) | no | Expose max_depth, grid_probs and selection_rule on dina_subgroup_bootstrap(), forward them to all three dina_subgroup() calls, and reject a supplied sg whose stored max_depth / selection_rule differ. | changes behaviour |
| A-096 | 2 | A | dina_subgroup_refit.R | 435 | .dina_signature_string | S2 | off | reproduced | The within-subgroup refit helpers and print.dina_subgroup_bootstrap() still assume a single-cut subgroup (scalar `if (direction == "left")`, `if (covariate %in% ...)`), so for a two-cut subgroup, which dina_subgroup() can return under its default max_depth = 2L, dina_subgroup_refit() stops with 'the condition has length > 1', dina_subgroup_bootstrap() drops its default refit with that message as a warning (or aborts if refit_confounders is a character vector), and print() of the bootstrap result errors; no reported number is wrong. Also at: dina_subgroup_refit.R:274 (.dina_resolve_confounders, only with a character-vector confounders; called outside the tryCatch at dina_subgroup_bootstrap.R:511-515, so there it aborts the bootstrap); ... (full list in the CSV) | no | Restore depth awareness in the three places that lost or never had it: make the bootstrap print method and the signature string/mask and confounder check work over the cut set (or, if conjunction refits stay unsupported, give dina_subgroup_refit() and the bootstrap an explicit depth-2 guard with a clear message, as the bootstrap had before c6be5ca3). | changes behaviour |
| A-097 | 2 | A | dina_subgroup_refit.R | 498 | print.dina_subgroup_refit | S2 | off | reproduced | The ratio-scale line of print.dina_subgroup_refit() is hard-coded to read "95% CI" although the limits printed are the exponentiated 1 - alpha limits, so any alpha other than 0.05 is mislabelled. | no | Use the computed ci_pct in the ratio-scale label instead of the literal 95. | changes behaviour |
| A-098 | 2 | A | explain_pareto_selection.R | 240 | .epx_classify_loss | S2 | off | reproduced | The loser classification and its wording hard-code the hrMaxSG/hrMinSG key order (band, N, Pcons, effect, K), so under sg_focus = "maxeff"/"maxeffCons" a loser is reported as losing on Pcons although the rule ranks on effect only, and under the default "hr" Pcons is described as "the secondary tiebreaker after N" although it is the primary key and N is not in the key. Also at: explain_pareto_selection.R:260 (Pcons step applied for maxeff / maxeffCons); explain_pareto_selection.R:267 (effect step applied for maxSG / minSG); ... (full list in the CSV) | no | Derive the loss-axis order and the tiebreaker wording from the sg_focus-specific sort key, with explicit branches for hr, maxeff and maxeffCons. | changes behaviour |
| A-099 | 2 | A | explain_pareto_selection.R | 341 | .epx_selected | S2 | off | reproduced | When the selected subgroup is not on the (effect, N) frontier, explain_pareto_selection() attributes this to `hrMinSG` whatever the sg_focus, promises frontier rows "below" that are never emitted, and then states that the frontier has a single candidate. Also at: explain_pareto_selection.R:130-137 (winner NULL -> loser list emptied); explain_pareto_selection.R:467-469 (single-candidate sentence keyed on zero losers) | no | In the off-frontier branch name the actual sg_focus, either list the frontier rows or drop the promise, and reserve the single-candidate sentence for a frontier of one row. | changes behaviour |
| A-100 | 2 | A | forestsearch_cross_validation.R | 501 | forestsearch_Kfold | S2 | off | reproduced | A fold refit that fails inside forestsearch() at a stage forestsearch() itself converts to a warning plus an sg.harm = NULL return (get_FSdata, subgroup.search, subgroup.consistency) is recorded by both CV loops as an ordinary no-subgroup (ITT) fold with the warning suppressed, so the CV's stated fatal-on-refit-error policy covers only errors that escape forestsearch(). Also at: forestsearch_cross_validation.R:485, 1008, 1025 (only try-error is tested; error_log never read); warning suppression at :464, :479, :950, :1001; ... (full list in the CSV) | no | Have the fold bodies distinguish a returned failure object (error_log present, or a consistency-stage failure flag) from a clean no-subgroup result and re-raise genuine failures with fold context, as is already done for a try-error. | changes behaviour |
| A-101 | 2 | A | forestsearch_cross_validation.R | 1417 | forestsearch_KfoldOut | S2 | off | reproduced | When the full-data analysis has a subgroup but no held-out subject of a CV run is assigned to a CV subgroup, the one-row agreement table is indexed as if its first row were the CV-H row, so sens_H is returned as 1 (the complement of its definition, which gives 0) and ppv_H as the prevalence of the original H (undefined by definition). Also at: forestsearch_cross_validation.R:1419 (ppv_H, same row-1 indexing); :1418 (sens_Hc forced NA in the mirror single-row case); consumers: forestsearch_cross_validation.R:616, 1171, 1201, 1229; ... (full list in the CSV) | yes | Build the agreement table on fixed 0/1 levels for both indicators so that a run with no CV-H subject cannot report sens_H = 1 (the author must choose 0, per the a/c definition, or NA, per the no-original-subgroup guard and the MR recovery's skip rule) and ppv_H is NA, then update the pinned expectation in test-cv-mixed-and-fatal.R:142. | changes behaviour |
| A-102 | 2 | A | forestsearch_cross_validation.R | 1672 | forestsearch_KfoldOut | S2 | off | reproduced | For a survival fit with est.scale = "1/hr", the two subgroup-table calls in forestsearch_KfoldOut(outall = TRUE) omit est.scale, so the CV subgroup rows are hazard ratios on the uninverted scale under labels swapped for the inverted scale, and the treat2 column / local treat.name that forestsearch_Kfold builds for the flip are never used. Also at: forestsearch_cross_validation.R:1688 (second subgroup-table call without est.scale); :574-576 (treat2 column and local treat.name written, never read); :577-578 (labels overwritten); ... (full list in the CSV) | no | Pass est.scale to the two survival subgroup-table calls in forestsearch_KfoldOut, and drop or actually use the treat2 / local treat.name assignment in forestsearch_Kfold. | changes behaviour |
| A-103 | 2 | A | forestsearch_cross_validation.R | 1798 | CV_sgs | S2 | off | reproduced | When the original subgroup has one factor, the CV 'Exact' find metric scores a fold as an exact match whenever either fold factor equals that factor, so a fold that selected a two-factor subgroup containing it (a different, smaller subgroup) counts as exact, and the package's other exact-match measure scores the same fold as not exact. Also at: forestsearch_cross_validation.R:1776 (depth-2 branch uses the same containment test); label at cv_summary_tables.R:418; the stricter measure at summarize_cv_results.R:940-943 | yes | In the depth-1 branch require the fold's subgroup to consist of the original factor alone for exact_match, keeping the present containment test for one_match and cov1_exact. | changes behaviour |
| A-104 | 2 | A | forestsearch_cross_validation.R | 1887 | find_covariate_any_match | S2 | off | reproduced | The covariate-any match misattributes covariates when one confounder name is a prefix or substring of another: the disambiguating substr() is one character short, and the final test is an unanchored grepl() of the bare name against the fold labels. Also at: forestsearch_cross_validation.R:1889 (negated-label variant of the off-by-one); :1902-1903 (unanchored grepl of the bare name, also treated as a regular expression) | unknown | Extract the variable token from the target and from each fold label and compare variable names for equality instead of using prefix charmatch, substr and substring grepl. | changes behaviour |
| A-105 | 2 | A | format_fs_details.R | 147 | format_fs_details | S2 | off | reproduced | The line wrapper inside format_fs_details() never terminates when, after an indent, the next token has no space or comma and is longer than the remaining width: the last-space fallback picks a space inside the indent, emits a blank piece and rebuilds the same remainder, so the loop runs forever while its result vector grows. Also at: format_fs_details.R:129-158 (the while loop), :144-151 (last-space fallback), :154-157 (blank piece appended, continuation prefix re-added) | no | When choosing the fallback break ignore spaces inside the leading or continuation indent, and hard-break at the width when no later break point exists so every pass consumes characters. | changes behaviour |
| A-106 | 2 | A | fpr_approximation.R | 112 | fpr_approx | S2 | off | reproduced | fpr_approx() requires effect_true > 0 and c1 > 0 before it looks at effect_type, so effect_type = "difference" stops on a zero or negative true difference (and on a zero or negative c1), although the difference branch itself would compute with those values. Also at: fpr_approximation.R:115 | no | Apply the positivity checks on effect_true and c1 only when effect_type is "ratio". | changes behaviour |
| A-107 | 2 | A | fpr_calibration.R | 236 | fpr_calibration | S2 | off | reproduced | With outcome_type_approx left at its NULL default and no outcome_type in fs_params (the documented minimum for a survival analysis), the P1 approximation silently assumes a binary outcome while the simulated forestsearch() runs are survival, so d_eff, P1 and L_eff_implied use n_min*p_event*(1-p_event) instead of n_min*(1-prop_cens) and a supplied prop_cens is ignored. Also at: fpr_calibration.R:245 (default arm of switch(): an unrecognised type string, e.g. a misspelt outcome_type_approx, also falls silently to d_eff_binary; otype is never validated); ... (full list in the CSV) | yes | Resolve the approximation outcome type as forestsearch() does (survival when fs_params carries none), and stop on a value outside the four supported types instead of falling back to binary. | changes behaviour |
| A-108 | 2 | A | fpr_calibration.R | 253 | fpr_calibration | S2 | off | reproduced | fpr_calibration() always computes the per-subgroup approximation P1 on the ratio scale and requires c1 > 0 and c2 > 0, so for identity-scale measures (RD, IRD, MD) the reported P1 and L_eff_implied are computed from log(c1), log(c2) and the identity-scale consistency threshold 0 is rejected outright. Also at: fpr_calibration.R:210; fpr_calibration.R:211; fpr_calibration.R:530 | no | Choose the approximation scale (and the null value of theta) from the effect measure and apply the positivity checks only to ratio-scale measures. | changes behaviour |
| A-109 | 2 | A | fpr_calibration.R | 259 | fpr_calibration | S2 | off | reproduced | fpr_calibration() injects c1/c2 as hr.threshold/hr.consistency without removing effect.threshold, consistency.threshold or seedit from fs_params, so an fs_params list carrying the preferred spellings silently overrides c1/c2 on identity-scale measures or makes every replicate error on ratio measures, and seedit in fs_params makes every replicate error. Also at: fpr_calibration.R:306 | no | Remove or reject the reserved names (both threshold spellings and seedit) from fs_params before building the call and report the first replicate error when all replicates fail. | changes behaviour |
| A-110 | 2 | A | fpr_calibration.R | 543 | fpr_calibration_transfer | S2 | off | reproduced | The validation simulation in fpr_calibration_transfer() (on by default, nsims_check = 200) cannot run with the calibrated configuration: the inner call is built with df.analysis = NULL and fs_params = list(), supplying them through ... duplicates those arguments, and the resulting error is swallowed, so fpr_sim, n_sim and in_ci come back NA / 0 / NA with no message. Also at: fpr_calibration.R:545; fpr_calibration.R:547; fpr_calibration.R:556 | no | Take the data, null_fn and fs_params for the check as explicit arguments (or store them in the calibration object), merge ... without duplicating names, and report a failed validation run instead of swallowing the error. | changes behaviour |
| A-111 | 2 | A | fs_fdr_report.R | 377 | .fs_fdr_one_rep | S2 | off | reproduced | The replicate harness counts a failure only on a thrown error, a NULL return or a 'subgroup.consistency\|subscript' warning, so a forestsearch() call that returns its failure-shaped list (sg.harm = NULL plus error_log, after get_FSdata or subgroup.search fails) is scored as a valid non-declaration, and a run in which every replicate failed reports A = 0 with n_valid = nsims and verdict CLEAN. Also at: fs_fdr_report.R:368-376 (warning handler promotes only messages matching 'subgroup.consistency\|subscript', muffles all others); ... (full list in the CSV) | no | Treat a forestsearch() return that carries error_log as a failed replicate (flag -1 with the error_log text) in .fs_fdr_one_rep and as n_failed in fpr_calibration's batch loop, and have fs_fdr_report() stop (or fall back to the GBSG set) once up front when the inherited confounders or conf.cont_jcuts names do not survive the intersection with the DGM columns. | changes behaviour |
| A-112 | 2 | A | fs_identification_summary.R | 321 | plot_fs_identification_involvement | S2 | off | reproduced | With a non-default `roles`, plot_fs_identification_involvement() pairs the i-th requested covariate with the i-th of the fixed anchor / true partner / proxy labels, so the legend attributes one covariate's involvement curve to another. Also at: fs_identification_summary.R:314-318 (labels and default palette built in fixed anchor/partner/proxy order) | no | Build each requested covariate's legend label from its own role, looked up by name, instead of by position. | changes behaviour |
| A-113 | 2 | A | generate_aft_dgm_helpers.R | 698 | process_factor_subgroup | S2 | off | reproduced | For a non-numeric subgroup variable, any cut specification other than a single value or list(type = "multiple") silently falls back to 'variable == first level', so the DGM's true subgroup is not the one requested. Also at: generate_aft_dgm_helpers.R:705 (any other non-matching spec, e.g. a vector of levels); ... (full list in the CSV) | no | Replace the two silent first-level fallbacks for unrecognised cut specifications with an error, keeping the first-level default only when cut_spec is NULL. | changes behaviour |
| A-114 | 2 | A | generate_aft_dgm_main.R | 226 | generate_aft_dgm_flex | S2 | off | reproduced | The exported generate_aft_dgm_flex() has the formal default treatment_var = NULL and its roxygen says NULL means treatment is simulated 50/50, but no code implements that, so a call that omits treatment_var stops with an obscure subscript error instead of building a DGM. Also at: generate_aft_dgm_main.R:19-20 (roxygen sentence); man/generate_aft_dgm_flex.Rd:57-58; generate_aft_dgm_helpers.R:431 (validation lets NULL through); generate_aft_dgm_helpers.R:459 (point of failure) | no | Either make treatment_var a required argument with a clear validation error and remove the roxygen sentence, or implement the documented random 50/50 assignment when it is NULL. | changes behaviour |
| A-115 | 2 | A | gg_forest.R | 226 | gg_forest | S2 | off | reproduced | The label panel decides bold/plain by comparing each row colour with cat_colours["ITT"], so any non-NULL cat_colours without an element named "ITT" gives geom_text an all-NA fontface and the returned plot fails when it is drawn. Also at: plot_subgroup_sims.R:274, plot_subgroup_sims.R:376 (forward x$cat_cols unchecked, do.call at :438); ... (full list in the CSV) | no | Decide which labels are bold from the category vector (or an explicit argument) instead of colour equality with an "ITT" entry, so that fontface can never be NA. | changes behaviour |
| A-116 | 2 | A | glm_effect_profile.R | 377 | .gep_resolve_family | S2 | off | reproduced | For effect_measure = "log_RR", overdispersion = "quasi" or "negbin" fits a plain Poisson model and also turns off the sandwich variance, so the profile's SEs and CIs are model-based Poisson while the result reports the requested overdispersion. Also at: glm_effect_profile.R:137 (validation admits log_RR for overdispersion); glm_effect_profile.R:370-379 (log_RR branch has no quasi/negbin arm); ... (full list in the CSV) | no | In the log_RR branch either honour the overdispersion choice (quasipoisson / glm.nb) or warn and ignore it while keeping the sandwich variance. | changes behaviour |
| A-117 | 2 | A | glm_spline_fit.R | 320 | glm_cs_fit | S2 | off | reproduced | With family = "gaussian", glm_cs_fit() compares the raw deviance difference (a residual-sum-of-squares difference in squared outcome units) to a chi-square without dividing by the dispersion, so the returned lrt_pvalue changes with the units of the outcome and is not a valid interaction p-value. Also at: glm_spline_fit.R:322 (pchisq on the unscaled statistic); glm_spline_fit.R:326-327 (printed as 'Chi-sq' when verbose); glm_spline_fit.R:461 (returned as lrt_pvalue); glm_spline_fit.R:489 (print.glm_cs_fit) | no | For the gaussian family, divide the deviance difference by the dispersion estimate of the larger model (or use the F form of the nested-model test) before computing the p-value, leaving binomial and poisson as they are. | changes behaviour |
| A-118 | 2 | A | grf_args.R | 189 | extract_grf_args | S2 | off | reproduced | For a fit made with subgroup_method = "grf", extract_grf_args() omits grf_selection, frontier_rule, effect_neighborhood and selection_rule and reports the stored return_selected_cuts_only, so the list it calls exact reproduces a policy-tree run with default band settings, not the frontier run that was made. Also at: grf_args.R:206-223 (GLM branch, same omission); grf_args.R:130-137 (roxygen claim of the exact list) | no | When the fit used subgroup_method = "grf", have extract_grf_args() rebuild the selection call (resolved grf_selection, frontier_rule mapped from sg_focus, selection_rule, effect_neighborhood, return_selected_cuts_only = TRUE), or document that it reproduces only the use_grf screening call. | changes behaviour |
| A-119 | 2 | A | guohe_algorithm3.R | 79 | .g3_coef | S2 | off | reproduced | For a binary outcome the closed-form scorer treats a candidate with an empty 2x2 cell as non-estimable, but the glm scorer used when fast = FALSE or adjust_covariates is supplied only checks arm sizes and that y takes two values, so it admits the separated fit's finite divergent coefficient and the two paths documented as the same estimator can select different candidates and return different numbers. Also at: guohe_algorithm3.R:84, guohe_algorithm3.R:89 (the two glm calls), guohe_algorithm3.R:61 (suppressWarnings), guohe_algorithm3.R:310-311 (closed-form zero-cell rule the model path lacks), guohe_algorithm3.R:183-187 and :27; ... (full list in the CSV) | no | Give the binary model path (and .ar_coef_se) the same empty-cell rule the closed form already applies, so both paths declare the same candidates non-estimable. | changes behaviour |
| A-120 | 2 | A | guohe_algorithm3.R | 112 | .g3_cox_lean | S2 | off | reproduced | The Cox scorers of the Guo-He comparator accept any subgroup with at least min_events total events, so a bootstrap resample in which one arm of a candidate has no events yields a monotone-likelihood fit whose finite divergent coefficient (warning suppressed) wins that resample's supremum and is averaged into the de-biased estimate and the bootstrap SD, although the package's own Cox estimator treats a zero-event arm as non-estimable. Also at: guohe_algorithm3.R:110-111 and :113-122 (.g3_cox_lean guards and suppressed fit), guohe_algorithm3.R:57, :64, :68-69, :75 (.g3_coef survival branch, same rule), guohe_algorithm3.R:399-401 (finite value enters the supremu; ... (full list in the CSV) | unknown | Apply the package's existing HR existence condition (at least one event in each arm) in .g3_cox_lean, the survival branch of .g3_coef and .ar_coef_se, so a resampled candidate with a zero-event arm is non-estimable and counted in drop_rate. | changes behaviour |
| A-121 | 2 | A | mrct_case_study.R | 282 | mrct_case_study | S2 | off | reproduced | mrct_case_study() restores the caller's RNG kind before running the identifier while mrct_region_sims() runs the identifier under the workers' L'Ecuyer-CMRG kind, so the identifier's set.seed()-seeded random steps (consistency splits with consistency_method = "split", the lasso folds with use_lasso = TRUE, the split fallback of the resample path) draw different streams and the documented claim that the case study returns the simulation's H-hat "under any generator" does not hold for those configurations. Also at: mrct_case_study.R:149-152 (the documented guarantee); mrct_trial_helpers.R:41 and mrct_simulation.R:534 (both force plan = "sequential", i.e. the plain-loop split path); ... (full list in the CSV) | no | Either keep the workers' generator kind in force for the identifier call too (restoring the caller's kind only on exit), or narrow the roxygen statement to the configurations in which the identifier consumes no random numbers. | changes behaviour |
| A-122 | 2 | A | mrct_inject_dgm.R | 414 | inject_mrct_structure | S2 | off | reproduced | In the planted-region mode region$loghr is not validated, so region = list(loghr = NULL) deletes the default 0 and, together with x3 or region_treat_loghr, the interaction log-HR is silently recycled into the region's prognostic coefficient. Also at: mrct_inject_dgm.R:489-490 (set_var keeps z_<region>, beta_var is NULL); mrct_inject_dgm.R:494 and :498 (the interaction value is appended); generate_aft_dgm_helpers_outcome.R:154-155 (silent recycling); ... (full list in the CSV) | no | Validate region$loghr in the planted mode as one finite number, as is already done for the observed mode. | add-only |
| A-123 | 2 | A | mrct_inject_dgm.R | 558 | inject_mrct_structure | S2 | off | reproduced | The region-dependent effect (region_loghr_shift / region_harm_loghr_shift) is added to dgm$df_super only, while dgm$df_source, the frame simulate_from_dgm(baseline = "fixed") draws from, keeps the unshifted potential outcomes, so a fixed-baseline simulation from a shifted DGM silently carries no region effect. Also at: mrct_inject_dgm.R:568 (only df_super is replaced); generate_aft_dgm_main.R:388-402, 424 (df_source built from the unshifted fit); simulate_from_dgm.R:241-245 (baseline = "fixed" reads dgm$df_source) | no | Apply the same shift to dgm$df_source using its own region and flag_harm columns, or refuse baseline = "fixed" for a DGM carrying $mrct$region_effect. | changes behaviour |
| A-124 | 2 | A | mrct_inject_dgm.R | 635 | analysis_covariates | S2 | off | reproduced | The pattern that removes the spline helper columns, "_treat$\|_k$", is not tied to the spline variable, so any real covariate whose z_ name ends in _treat or _k (e.g. z_prior_treat) is silently dropped from the covariate set handed to forestsearch(). | no | Remove exactly the three helper columns built from dgm$model_params$spline_info$var instead of matching on name endings. | changes behaviour |
| A-125 | 2 | A | mrct_simulation.R | 319 | mrct_region_sims | S2 | off | reproduced | With confounders.name left at its formal default NULL, mrct_region_sims() hands forestsearch() every z_ column of dgm$df_super, which for any spline DGM (always the case for create_dgm_for_mrct()) includes the DGM-internal columns z_bm_treat, z_bm_k and z_bm_k_treat (and the region indicator), so the search evaluates and can declare cuts on columns that are not baseline covariates of the simulated trial. Also at: mrct_case_study.R:288-291 (mrct_case_study() rebuilds the same default set); roxygen examples that rely on the default: mrct_simulation.R:232-239, mrct_case_study.R:65-67 | yes | Remove the region variable and the spline helper columns (names built from dgm$model_params$spline_info$var, not from a name-suffix pattern) from the default confounder set in mrct_region_sims() and in mrct_case_study(). | changes behaviour |
| A-126 | 2 | A | mrct_simulation.R | 422 | mrct_region_sims | S2 | off | reproduced | A user-supplied sim_args$entry_var is silently removed whenever dgm$df_super has no column named 'entrytime', because the override list carries entry_var = NULL and modifyList() deletes a key whose new value is NULL, so entry times are drawn Uniform(0, max_entry) instead of read from the named column. Also at: mrct_case_study.R:277 | no | Add entry_var = "entrytime" to the override list only when that column exists, so a user-supplied sim_args$entry_var is otherwise left in place (or reject it with an error if it is not meant to be supported). | changes behaviour |
| A-127 | 2 | A | mrct_simulation.R | 570 | mrct_region_sims | S2 | off | reproduced | An error thrown by forestsearch() inside a replicate is caught and turned into NULL, and the replicate is returned as an ordinary 'no subgroup found' row (any_found = 0, sg_found = 'none') with no warning, count or flag, so a call that fails in every replicate (for example a misspelled name in fs_args) yields a clean table with a 0% identification rate. Also at: mrct_simulation.R:508 (forestsearch skipped when the region Cox fit gives hr_test = NA; the replicate is likewise returned with any_found = 0) | no | Record the identifier's failure per replicate (an error column plus one warning giving the count and the first message, or stop when every replicate fails) instead of coding it as any_found = 0. | changes behaviour |
| A-128 | 2 | A | mrct_trial_helpers.R | 60 | .mrct_analyze_trial | S2 | off | reproduced | mrct_trial_analysis() does not check the declared rule's covariates for missing values, so a region patient with NA in a variable of the declared rule gets in_Hhat = NA, is silently left out of the region(G-hat) fit although not counted in H-hat either, and makes the 2 x 2 counts NA, whereas mrct_four_sets() rejects an in_Hhat containing NA. Also at: mrct_four_sets.R:126 (.mrct_check_trial checks NA only on the outcome, event, treatment, region and truth columns); mrct_trial_helpers.R:91 and :108-109 (set index and 2 x 2 counts built from in_Hhat); ... (full list in the CSV) | no | Stop (or state the handling) when the declared rule evaluates to NA on any trial row, mirroring mrct_four_sets()'s in_Hhat check, or check the confounder columns for missing values up front. | add-only |
| A-129 | 2 | A | pareto_frontier_table.R | 288 | pareto_frontier_table | S2 | off | reproduced | compute_frontier_cis() accepts ci_level but returns no record of it, and every consumer labels the intervals '95%', so a table or explanation built from a ci_table computed at another level displays those intervals as 95% intervals. Also at: pareto_frontier_table.R:296; pareto_frontier_table.R:314; pareto_frontier_table.R:510; explain_pareto_selection.R:365; frontier_cis.R:328 (returned table carries no ci_level) | no | Carry ci_level on the returned table (attribute or column) and build the column, spanner and sentence labels from it. | changes behaviour |
| A-130 | 2 | A | plot_km_band_forestsearch.R | 351 | plot_km_band_forestsearch | S2 | off | reproduced | sg_labels are handed to weightedsurv as row-subsetting expressions while the validated 0/1 indicator columns are never used, so the default bare-column label, an automatic label for a single binary-indicator subgroup such as '{meno}', and a reference subgroup's display label are evaluated as code, giving a numeric 0/1 row index (row 1 repeated) or a parse/eval error instead of the intended subgroup. Also at: plot_km_band_forestsearch.R:292 (default sg_labels <- sg_cols); :459-461 (Mode 1 label from sg.harm_label); :666, :269, :273 (ref_subgroups display label appended to sg_labels); roxygen :31 and :68 | no | Pass expressions built from the indicator columns (the documented "<col> == 1" form) to weightedsurv and keep the readable labels for the legend only. | changes behaviour |
| A-131 | 2 | A | plot_km_band_forestsearch.R | 487 | generate_readable_sg_labels | S2 | off | reproduced | When a fit has no grp.consistency$out_sg$sg.harm_label (subgroup_method = "dina" or "grf") a multi-factor sg.harm is pasted into one string before generate_complement_expression(), whose single-string branch flips only the first matching operator type, so the label handed to weightedsurv as the 'complement' selects a strict subset of the true complement and the second curve and risk row are drawn for the wrong patients. Also at: plot_km_band_forestsearch.R:522-549 (single-string operator flip), 613-617 (conjunction collapsed to one string) | no | Negate a compound expression as a whole (wrap it in !(...)) and reserve the operator flip for a single simple comparison, or subset from the Qrecommend/Brecommend indicators already created. | changes behaviour |
| A-132 | 2 | A | plot_sg_results.R | 214 | plot_sg_results | S2 | off | reproduced | est.scale = "1/hr" only reverses the arm labels and swaps the two subgroup names: the reversed treatment column df$treat_plot is created and never used, so every hazard ratio, KM curve and summary is computed exactly as for "hr" and the figure shows un-inverted HRs beside reversed labels. Also at: plot_sg_results.R:215, 221 (treat_plot assigned, never read), 245-254, 278-296, 598-613, 629-632, 764-765 | no | Either use the reversed treatment column for the Cox fits and curves when est.scale = "1/hr" so labels and estimates agree, or remove the dead treat_plot assignment and document that only labels change. | changes behaviour |
| A-133 | 2 | A | plot_sg_results.R | 504 | compute_sg_summary | S2 | off | reproduced | compute_sg_summary() extracts the event totals with vapply(..., integer(1), "events"), so plot_sg_results() and plot() on a forestsearch fit stop with a vapply type error, for every plot type, whenever the 0/1 event column is stored as double. Also at: plot_sg_results.R:505 (same integer(1) template on "events" inside event_rate); double value produced at plot_sg_results.R:472, 481, 490; entry point forestsearch_methods.R:621 | no | Extract the event counts with a numeric template (or coerce the event sums to integer) so a double-stored 0/1 event column is accepted. | changes behaviour |
| A-134 | 2 | A | plot_subgroup_results_forestplot.R | 597 | plot_subgroup_results_forestplot | S2 | off | reproduced | conf.level changes the column header and the unadjusted rows but not the bias-corrected rows, which copy the bootstrap object's fixed 95% limits, so with conf.level != 0.95 the forest plot mixes interval levels under one "HR (xx% CI)" header. Also at: plot_subgroup_results_forestplot.R:598-600; 1346-1349 (create_subgroup_summary_df, no caller); plot_sg_glm_outcomes.R:317-323 | no | Recompute the bias-corrected limits at conf.level from the stored estimate and standard error, or state the bootstrap's own level in the row label or footnote (or reject conf.level != 0.95 when bias-corrected rows are present). | changes behaviour |
| A-135 | 2 | A | run_simulation_analysis.R | 410 | run_simulation_analysis | S2 | off | reproduced | The documented argument cox_formula ('Optional Cox formula for unadjusted ITT') is threaded through four internal functions and never read, so hr.itt is always the fit of Surv(outcome, event) ~ treat whatever formula is supplied; debug is likewise accepted, documented and never read. Also at: run_simulation_analysis.R:263 (roxygen); pass-through only at :594, :619, :644, :681, :775, :847, :923, :1173, :1199, :1264, :1307; debug: :271, :416, :647, :1175, :1202, :1267, :1309; ... (full list in the CSV) | no | Either use cox_formula for the ITT fit when it is supplied, or remove or deprecate it (and debug) with a message. | changes behaviour |
| A-136 | 2 | A | run_simulation_analysis.R | 587 | run_simulation_analysis | S2 | off | reproduced | On the legacy boolean path (methods = NULL, the default) an explicit fs_params = list(n.min = NULL) is deleted by utils::modifyList() before it reaches forestsearch(), so the sample-size-adaptive rule never fires and n.min silently stays at the formal default 60. Also at: run_simulation_analysis.R:612-615 (same merge for the FSlg arm); misleading comment at :811-814 | no | Use the NULL-preserving merge on the legacy boolean path as well. | changes behaviour |
| A-137 | 2 | A | run_simulation_analysis.R | 818 | .run_fs_analysis_gen | S2 | off | reproduced | An error raised by forestsearch() (or the standalone GRF call) inside a simulation replicate is downgraded to a warning and the replicate is returned as an ordinary 'no subgroup found' row (any.H = 0, size.H = 0, sens = 0, ppv = 0, spec/npv filled) with no column that marks it as failed. Also at: run_simulation_analysis.R:1257-1260 (standalone GRF call, same warning-then-NULL pattern); run_simulation_analysis.R:1195-1211 (GRF function not found -> warning + non-detection row) | unknown | Record a per-row completion status (and the error message) and leave the detection and classification columns NA for failed replicates so that summaries can count or exclude them. | add-only |
| A-138 | 2 | A | run_simulation_analysis.R | 953 | .extract_fs_estimates_gen | S2 | off | reproduced | In GLM simulation mode the closure that fills hr.itt / hr.H.hat / hr.Hc.hat / hr.H.true / hr.Hc.true returns a log odds ratio for effect_measure = "RD", an odds ratio for "RR", a log rate ratio for "IRD", and NA for every estimate when effect_measure is left NULL (with ahr.* then computed by the odds formula even for continuous or count outcomes), instead of the effect on the named or default scale. Also at: run_simulation_analysis.R:961 (log_scale); run_simulation_analysis.R:977-978 (return value); run_simulation_analysis.R:1108 and :1137 (NULL measure falls back to "OR" for compute_aor() whatever the outcome type); ... (full list in the CSV) | no | Resolve a NULL effect_measure once from outcome_type as forestsearch() does, compute each supported measure on its own scale (or delegate to the package's estimator layer) in both extractors, and stop on an unsupported measure instead of returning NA. | changes behaviour |
| A-139 | 2 | A | run_subgroup_sims.R | 464 | validate_subgroups | S2 | off | reproduced | validate_subgroups() evaluates each subgroup expression with enclos = parent.frame(), so a free name that exists only in the user's workspace (not in the data or in `cutpoints`) passes validation, while the per-trial evaluation on a multisession worker cannot see it and swallows the error, giving a silent all-NA column that the validation step exists to prevent. Also at: run_subgroup_sims.R:538-541 (per-trial evaluation; tryCatch turns the error into NULL and the subgroup is skipped); run_subgroup_sims.R:799 (call site) | no | Evaluate the validation expressions with only the visibility the per-trial evaluation has on a worker (data columns, cut-points, benchmark columns, base functions), and surface rather than swallow per-trial expression errors. | changes behaviour |
| A-140 | 2 | A | setup_gbsg_dgm.R | 155 | setup_gbsg_dgm | S2 | off | reproduced | setup_gbsg_dgm() does not put a lin_pred_obs column on df_super, so a direct call simulate_from_dgm(<GBSG dgm>) with n left at its default NULL stops with an opaque data-frame error ('replacement has 0 rows, data has <n_super>') instead of returning the super-population simulation the roxygen describes. Also at: simulate_from_dgm.R:229-234 (n = NULL branch sets only follow_up and treat_sim), simulate_from_dgm.R:310 (reads lin_pred_obs), simulate_from_dgm.R:350 (where the error surfaces) | no | Give the GBSG object the column the n = NULL path reads (lin_pred_obs equal to the stored lin.conf.true), or have the n = NULL branch of simulate_from_dgm() derive it from treat and lin_pred_1 / lin_pred_0 when it is absent. | add-only |
| A-141 | 2 | A | setup_gbsg_dgm.R | 175 | setup_gbsg_dgm | S2 | off | reproduced | With cens_type = "uniform", setup_gbsg_dgm() builds model_params$censoring with no min/max, so simulate_from_dgm() calls runif(n, min = NULL, max = NULL) and stops with an uninformative error. Also at: simulate_from_dgm.R:330-333 | no | Either populate censoring$min and censoring$max for the uniform case in setup_gbsg_dgm() or reject cens_type = "uniform" there with an explicit error. | changes behaviour |
| A-142 | 2 | A | simulate_from_dgm.R | 274 | simulate_from_dgm | S2 | off | reproduced | A strata_rand value that is not a column of the sampled data is silently ignored and marginal allocation is used in place of the requested stratified allocation. | no | Stop (or at least warn) when strata_rand is supplied but is not a column of the sampled frame. | changes behaviour |
| A-143 | 2 | A | simulation_tables.R | 132 | invert_hr_columns | S2 | off | reproduced | The benefit-notation inversion list omits cde.H.true and cde.Hc.true, so format_oc_results(subgroup_notation = "benefit") shows the true-subgroup CDE means on the switched scale beside identified-subgroup CDE means (and all HR/AHR columns) on the original scale. Also at: oc_analyses.R:658-660 (call), oc_analyses.R:732-737 (only consumer of the un-inverted columns), oc_analyses.R:768-771 (output columns), oc_analyses.R:919-924 and :992-1002 (G / Gc labels and common spanner) | no | Add the two true-subgroup CDE columns to the inversion list so every effect column is transformed together. | changes behaviour |
| A-144 | 2 | A | simulation_tables.R | 645 | build_estimation_table | S2 | off | reproduced | With subgroup_notation = "benefit" and a DGM that has no stored CDE, the CDE truth is recomputed from the raw super-population after the DGM has already been inverted, so the b-ddagger bias column, the CDE rows and the group header compare inverted estimates with an un-inverted CDE truth. Also at: simulation_tables.R:1367 (interpret_estimation_table: same fallback after the inversion at :1104-1105; un-inverted cde_H then used at :1384-1386 and in the bias-corrected CDE sentence) | no | Compute the fallback CDE before inverting the DGM (or apply the same inversion to the fallback values) in both functions. | changes behaviour |
| A-145 | 2 | A | simulation_tables.R | 1279 | interpret_estimation_table | S2 | off | reproduced | The generated narrative contains fixed evaluative sentences (low detection rate confirms type-I error control; bias is upward; bias correction improves accuracy) that are not conditioned on the detection rate or on the sign and size of the biases computed in the same call, so the text can contradict the figures it prints. Also at: simulation_tables.R:1311-1317 ("relative bias above the true value ... This upward bias ... worst-looking patients ... modestly better", printed whatever the sign); ... (full list in the CSV) | no | Condition each evaluative sentence on the computed detection rate and on the sign and relative size of the computed biases, or remove the evaluative wording from the template. | changes behaviour |
| A-146 | 2 | A | summarize_bootstrap_subgroups.R | 387 | summarize_bootstrap_subgroups | S2 | off | reproduced | For maxk >= 3 the agreement table labels each bootstrap subgroup from M.1 and M.2 only, so three-factor subgroups lose their third factor and distinct subgroups sharing the first two factors are counted as one row. Also at: summarize_bootstrap_subgroups.R:382 (branch maxk >= 2), :394 (group by Subgroup, K_sg); summarize_bootstrap_results.R:387-391 (maxk passed is always 5) | no | Build the agreement label from all non-empty M.1..M.maxk columns, as the original_agreement section already does. | changes behaviour |
| A-147 | 2 | A | summarize_bootstrap_subgroups.R | 417 | summarize_bootstrap_subgroups | S2 | off | reproduced | Both cut() calls use R's default right-closed intervals while the bin labels describe left-closed bins, so a value exactly on a break is counted one row below its label (Pcons 0.95 under "0.85-0.95" instead of ">=0.95"; N_sg 50 under "<50", 100 under "50-99", 300 under "200-299"). Also at: summarize_bootstrap_subgroups.R:444 (N_sg bins, same pattern); correct sibling for contrast: summarize_cv_results.R:764-765 (right = FALSE) | yes | Bin both variables left-closed, as the CV sibling already does, so that each count matches its label. | changes behaviour |
| A-148 | 2 | A | synthetic_data_perturbation.R | 168 | generate_bootstrap_synthetic | S2 | off | reproduced | The ordinal adjacent-value step calls sample(possible_new, 1) on what is a length-1 numeric at either end of the scale, so base sample() draws from 1:possible_new instead of returning the single adjacent value. Also at: synthetic_data_perturbation.R:163-165 (possible_new has one element at either end of the scale) | no | Select the replacement by index, or guard the length-1 case, so that a single candidate is returned unchanged. | changes behaviour |
| A-149 | 2 | A | synthetic_data_perturbation.R | 180 | generate_bootstrap_synthetic | S2 | off | reproduced | A numeric two-valued categorical variable is flipped with 1 - x, which is only right for 0/1 coding; any other numeric coding (for example 1/2) silently yields values that are not in the data. | no | Flip a numeric binary to the other observed value, as the factor/character branch already does, instead of computing 1 - x. | changes behaviour |
| A-150 | 2 | A | synthetic_data_perturbation.R | 528 | generate_bootstrap_with_noise | S2 | off | reproduced | With continuous_vars and cat_vars both left NULL (the default), every numeric column is classed continuous and then removed from the auto-detected categorical set, so numeric categorical variables with 10 or fewer values are never perturbed as categorical although the roxygen says they are. Also at: synthetic_data_perturbation.R:514-517; synthetic_data_perturbation.R:522-525 | yes | In auto-detection, classify low-cardinality numeric columns as categorical first and exclude them from the auto-detected continuous set. | changes behaviour |
| A-151 | 2 | A | test_hte_crump.R | 399 | test_hte_from_forestsearch | S2 | off | reproduced | test_hte_from_forestsearch() uses fit$confounders.evaluated as column names of the user's data, but a default forestsearch() fit stores cut expressions such as "age <= 50" there for every continuous covariate, so those covariates are silently left out of the HTE test (only bare-name categorical covariates are tested) and the function returns NULL with no message when all covariates are continuous. Also at: test_hte_crump.R:410-411 (labels that are not columns of df skipped with no message); test_hte_crump.R:416 (silent NULL when nothing matched, unlike the warning at :406); ... (full list in the CSV) | yes | Resolve the evaluated cut labels to their underlying data columns (or take the original confounder names from the stored call arguments) and warn or stop when a label cannot be resolved or no covariate remains. | changes behaviour |
| A-152 | 2 | A | truefind_asymptotic.R | 652 | find_required_sample_size | S2 | off | reproduced | With tol below 1 the bisection loop never terminates, because once the integer bounds are one apart the rounded midpoint equals one of the bounds and neither bound moves. | no | Validate tol >= 1 on entry, or stop the loop when the midpoint coincides with a bound. | add-only |
| A-153 | 2 | A | truefind_asymptotic_glm.R | 241 | compute_detection_probability_glm | S2 | off | reproduced | The roxygen for c2 says the default is 1.0 on the ratio scale and 0.0 on the difference scale, but the formal default is a fixed 1.0, so a direct call with effect_scale = "difference" and c2 omitted silently uses a consistency threshold of 1.0 and returns a much smaller probability than the documented default would give. Also at: truefind_asymptotic_glm.R:179-180 (roxygen statement of the default); man/compute_detection_probability_glm.Rd:37-38 (same text; the Usage block at line 11 shows c2 = 1) | no | Either make the c2 default depend on effect_scale as the documentation says (0.0 on the difference scale) or correct the roxygen to say the default is 1.0 on both scales and that c2 must be supplied for difference measures. | changes behaviour |
| A-154 | 2 | A | forestsearch_methods.R | 277 | print.forestsearch | S2 | off | traced | For GLM fits print() and summary() show the selected subgroup's effect under the label 'HR' / 'Hazard ratio', and for OR, RR and IRR the number shown is the un-exponentiated log-scale estimate. Also at: forestsearch_methods.R:499 ('Hazard ratio:'); summarize_bootstrap_subgroups.R:227 ('Hazard Ratio (hr_sg)' heading for every outcome type); ... (full list in the CSV) | no | Label the line with the object's effect_measure and exponentiate the value for OR/RR/IRR, as the Search Space block of summary() already does. | changes behaviour |
| A-155 | 2 | A | fs_family_report.R | 293 | fs_family_report | S2 | off | traced | For an argument list with outcome_type = "binary" and no effect_measure, fs_family_report() resolves the measure to "RD" while forestsearch() resolves it to "OR", so the estimability-boundary row is reported as inert ('no existence condition') although the engine applies the four-cell OR condition, and status_counts and the verdict sentence count one data-dependent stage too few. Also at: same binary = "RD" fallback, not traced here: frontier_cis.R:138; plot_sg_glm_outcomes.R:151; fs_dgm_feasibility.R:153; forestsearch_cross_validation.R:1449, :1498 | no | Resolve a missing effect_measure for binary to "OR", the same default forestsearch() uses. | changes behaviour |
| A-156 | 2 | A | fs_fdr_report.R | 258 | fs_fdr_report | S2 | off | traced | The harm null's k_inter is calibrated on calibrate_k_inter()'s default 5000-row super-population but the DGM is then built with fs_fdr_report()'s n_super (default 25000), so the DGM labelled 'target HR_H = 1.0' is not the population the root was found on and its realised HR_H is not the calibrated value. Also at: fs_fdr_report.R:261-262 (DGM built at n_super); sim_aft_gbsg.R:1024-1027 (calibration builds at the 5000 default) | yes | Forward n_super (and the seed) to calibrate_k_inter() so that the root is found on the same super-population the null is simulated from. | changes behaviour |
| A-157 | 2 | A | generate_aft_dgm_helpers.R | 1370 | compare_multiple_survreg | S2 | off | traced | The exported argument criteria is documented as choosing the comparison criteria but is never read: AIC and BIC are always both computed and the best model is always the AIC winner. Also at: generate_aft_dgm_helpers.R:1346 (roxygen); :1593, :1602, :1609, :1616 (best_model always the AIC winner) | no | Either honour criteria in the ranking and best-model choice or remove the argument and its roxygen entry. | changes behaviour |
| A-158 | 2 | A | mr_method_tables.R | 91 | mr_method_comparison | S2 | off | traced | mr_method_comparison() decides whether the bootstrap point is on the 1/hr scale from the fit's args_call_all$est.scale instead of the scale the bootstrap used, so for a GLM fit on an identity measure made with est.scale = "1/hr" the bootstrap column shows the reciprocal of the difference. Also at: bootstrap_dofuture_main.R:559-568, 784; forestsearch_main.R:382, 1637, 1881 | no | Read the scale from the bootstrap result (boot$est.scale) rather than from the fit's call. | changes behaviour |
| A-159 | 2 | A | mrct_inject_dgm.R | 226 | make_region_model | S2 | off | traced | make_region_model() called with its own formal defaults for both term arguments (x_pred = NULL, x1_vars = character(0), each documented as a valid 'none' value) stops with an unrelated base-R error, because the odds-ratio vector becomes NULL. Also at: mrct_inject_dgm.R:228, :250, :256, :270 | no | Build ors as a zero-length numeric vector when there are no terms, so the linear predictor reduces to the intercept plus the optional band term. | add-only |
| A-160 | 2 | A | oc_analyses.R | 648 | format_oc_results | S2 | off | traced | With subgroup_notation = "benefit", a dgm supplied and no subtitle, the auto-generated subtitle prints the DGM's subgroup effects on the switched-treatment scale under the G/Gc labels while the table body is inverted to the original scale. Also at: oc_analyses.R:641-655 (auto-subtitle block), :657-660 (only `results` is inverted); siblings that invert the DGM: simulation_tables.R:297-301, :599-602, :1105 | no | Build the subtitle from the DGM after the same benefit inversion the sibling tables apply (via get_dgm_hr()), with a fallback when the subgroup effects are absent. | changes behaviour |
| A-161 | 2 | A | plot_subgroup_results_forestplot.R | 917 | plot_subgroup_results_forestplot | S2 | off | traced | On difference scales (MD by default, RD/IRD with xlim_method = "data") the axis-limit outlier test is a ratio test that only works for positive values: when the two lowest values are both negative it always fires and sets the lower limit to second/extreme_ci_cap, inside the second-lowest value, and when the second value is 0 it sets the limit to 0, so CI bounds that are not outliers are cut off, printed as "<x"/">x" and overwritten in the returned data. Also at: plot_subgroup_results_forestplot.R:908, 973-986, 1174-1191 | no | Apply the outlier comparison to distances from the reference value on each side separately, so it fires only when the most extreme bound is more than twice as far from zero as the next one. | changes behaviour |
| A-162 | 2 | A | synthetic_data_perturbation.R | 145 | generate_bootstrap_synthetic | S2 | off | traced | ordinal_vars is honoured only for variables that are also listed in cat_vars, so in the documented call pattern and in generate_gbsg_bootstrap_general() the ordinal variable is validated but never perturbed. Also at: synthetic_data_perturbation.R:152 (ordinal branch inside the cat_vars loop), :217 (remaining-variables note ignores ordinal_vars), :20-29 and :42-50 (roxygen examples), :368-378 (generate_gbsg_bootstrap_general call) | yes | Run the categorical loop over the union of cat_vars and ordinal_vars and count ordinal_vars as specified in the remaining-variables note. | changes behaviour |
| A-163 | 2 | A | sim_aft_gbsg.R | 638 | .create_gbsg_dgm_ | S2 | off | suspected | The censoring linear-predictor coefficients are back-transformed with the event model's scale (sigma) rather than the censoring model's (sigma_cens), so the stored censoring predictors equal z'gamma_cens times sigma/sigma_cens and do not match the fitted censoring AFT coefficients that the same object stores in cens_params$model$gamma. Also at: sim_aft_gbsg.R:646-647 (linC1.conf / linC0.conf built from b_cens); sim_aft_gbsg.R:652-657 (cens_model$gamma stores the unscaled gamma_cens) | yes | Decide which censoring model the GBSG DGM is meant to use and make the stored predictors and the stored censoring coefficients agree with it. | method-adjacent |
| E-003 | 2 | E | Cox_estimation_helpers.R | 33 | get_Cox_sg and 27 others | S3 | on | traced | 28 on-path Tier 2 functions (bootstrap, bias correction and MR helpers) are never referenced under tests/testthat/ (static whole-word grep), for example get_Cox_sg(), build_cox_formula(), fit_cox_models(), fit_subgroup_effect(); the full list is in section 8. | no | Add tests that call these functions directly, starting with the on-path ones that carry selection or estimation logic. | add-only |
| E-002 | 1 | E | Cox_estimation_helpers.R | 137 | .fs_adjust_terms and 109 others | S3 | on | traced | 110 of the 158 internal Tier 1 functions are never referenced under tests/testthat/ (static whole-word grep), for example .fs_adjust_terms(), .consistency_cox_pieces(), .consistency_adj_terms(), .dina_one_pass(); the full list is in section 8. | no | Add tests that call these functions directly, starting with the on-path ones that carry selection or estimation logic. | add-only |
| A-164 | 2 | A | bootstrap_calculations_helpers.R | 221 | format_CI | S3 | on | traced | The exported format_CI() documents its input as 'Data frame or data.table' but selects columns with the data.table-only '..' prefix, so a plain data frame input errors. Also at: bootstrap_calculations_helpers.R:206 | no | Select the columns in a way that works for both classes, or narrow the documentation to data.table. | add-only |
| C-002 | 2 | C | bootstrap_dofuture_main.R | 102 | forestsearch_bootstrap_dofuture | S3 | on | traced | The roxygen of forestsearch_bootstrap_dofuture(), bootstrap_results() and get_dfRes() describes arguments, return elements, seeding and warning behaviour that the current code does not have. Also at: bootstrap_dofuture_main.R:180-181, 186; bootstrap_analysis_dofuture.R:87-89, 144-145, 200, 204, 211-212, 230; bootstrap_calculations_helpers.R:233-234, 247-249; bootstrap_dofuture_main.R:30 | no | Bring the roxygen of the three functions (arguments, return list, seeding, warning behaviour, examples) in line with the current code. | add-only |
| D-002 | 2 | D | bootstrap_dofuture_main.R | 540 | forestsearch_bootstrap_dofuture | S3 | on | traced | The three return paths of forestsearch_bootstrap_dofuture() build fs_bootstrap objects with different element sets. Also at: bootstrap_dofuture_main.R:621-636; bootstrap_dofuture_main.R:722-730; bootstrap_dofuture_main.R:769-789; bootstrap_dofuture_main.R:905-906 | no | Assemble the return list in one place so every exit carries the same elements, with NULL where a value is not available, and word the print fallback to cover the zero-success case. | add-only |
| A-165 | 2 | A | bootstrap_dofuture_setup_helpers.R | 22 | resolve_bootstrap_parallel_args | S3 | on | traced | With the bootstrap's default parallel_args = list() (documented as inheriting the fit's settings), the bootstrap stops with an error when the fit itself was run with parallel_args = list(), which forestsearch() accepts and runs sequentially. Also at: bootstrap_dofuture_main.R:366; subgroup_consistency_helpers.R:41 | no | Fall back to a defined plan (for example sequential, as forestsearch() effectively does for an empty list) when neither the bootstrap call nor the fit supplies one. | changes behaviour |
| A-166 | 1 | A | consistency_resample.R | 466 | consistency_resample | S3 | on | traced | With the default method = "both" (or "mc") the multiplier matrix has nrow(df) columns while dfbeta has one entry per row the fit actually used, so any NA-dropped row makes the exported function stop with a non-conformable error, and an NA event value fails the Cox guard. Also at: consistency_resample.R:62 (event guard not NA-safe); consistency_resample.R:90 and 322 (n = nrow(df)); consistency_resample.R:460-465 | no | Size the multiplier draws by length(pieces$dfbeta) (and report that as n) and make the event guard NA-safe. | add-only |
| E-001 | 1 | E | dina.R | 157 | dina_fit and 24 others | S3 | on | traced | 25 exported Tier 1 functions are never referenced under tests/testthat/ (static whole-word grep), for example dina_fit(), evaluate_comparison(), get_conf_force(), get_conf_force_jq(); the full list is in section 8. | no | Add tests that call these functions directly, starting with the on-path ones that carry selection or estimation logic. | add-only |
| C-003 | 1 | C | dina_subgroup.R | 241 | dina_subgroup | S3 | on | traced | The candidates data frame returned by dina_subgroup() (the qualifying family that forestsearch's effect re-ranking and MR consume) is not described in @return, and it is present only in the found = TRUE branch. Also at: dina_subgroup.R:631-635 (candidates built, found = TRUE only); dina_subgroup.R:481-506 (found = FALSE list has no candidates element); dina_subgroup.R:238 (tau_sign text refers to candidates$tau_hat) | no | Add a candidates item to @return, stating that it is absent (NULL) when found = FALSE. | add-only |
| C-004 | 1 | C | dina_subgroup.R | 1095 | dina_frontier | S3 | on | traced | The roxygen 'Cut form' paragraph of dina_frontier() says cut_expr is always '<covariate> <= <threshold>' with direction kept for reference only, but the code writes '>=' for right-direction rows. Also at: dina_subgroup.R:1286-1291 (code and its comment); dina_subgroup.R:1278-1279 (stale inline comment 'render the canonical "var <= t" cut'); dina_subgroup.R:1156 (@return 'the canonical forestsearch cut expression') | no | Rewrite the 'Cut form' paragraph (and the stale inline comment at 1278-1279) to describe direction-specific operators. | add-only |
| C-005 | 1 | C | forestsearch_main.R | 1394 | forestsearch | S3 | on | traced | Several roxygen passages describe superseded behaviour: the Details algorithm overview, the 'default split-consistency search' wording, the unqualified 'both spellings disagreeing is an error' claim, and the d0.min scope. Also at: forestsearch_main.R:1398, 1403-1421, 951-953, 661-662, 668-669, 887 | no | Update the listed roxygen passages to the current defaults and scopes. | add-only |
| F-001 | 1 | F | forestsearch_main.R | 1619 | forestsearch | S3 | on | traced | randomForest is a hard dependency (DESCRIPTION Imports plus an @importFrom on forestsearch()) that no code in R/ calls. Also at: DESCRIPTION:50; NAMESPACE:489; plot_pareto_frontier.R:59 | no | Remove the unused randomForest import and the Imports entry. | add-only |
| A-167 | 1 | A | forestsearch_main.R | 1628 | forestsearch | S3 | on | traced | The documented formal flag_harm.name is accepted by forestsearch() but nothing in the function body or in any consumer of args_call_all reads it. Also at: forestsearch_main.R:349 (roxygen) | no | Document flag_harm.name as a recorded-only label or remove the formal. | add-only |
| A-168 | 1 | A | forestsearch_main.R | 1789 | forestsearch | S3 | on | traced | For outcome_type = "survival" a supplied effect_measure is documented as ignored but is never reset, so it is returned in the result, can trigger the identity-scale threshold warning, and makes the DINA path stop. Also at: forestsearch_main.R:2186 and 298-310 (stray measure reaches the bounded-identity warning), 3236 and forestsearch_helpers.R:1446-1459, 1524 (DINA refusal), forestsearch_main.R:3969 (returned unreset); ... (full list in the CSV) | no | Reset effect_measure to NULL (with a message) when outcome_type is "survival", as part of the same entry validation. | add-only |
| A-169 | 1 | A | forestsearch_main.R | 2017 | forestsearch | S3 | on | traced | The sg_focus = "maxeff" override block and its warning fire regardless of subgroup_method, so on the DINA and GRF paths the call warns that the effect floor is disabled and that consistency settings were overridden although those paths keep the effect floor and never read the overridden values. Also at: forestsearch_main.R:2032-2041 (warning text); 2043-2049 (split-cost warning); 2446, 2495 (threshold_config$pconsistency) | no | Condition the maxeff override block and its warning on the consistency engine, using the same predicate as .fs_admission_applies(). | changes behaviour |
| D-003 | 1 | D | forestsearch_main.R | 2813 | forestsearch | S3 | on | traced | The object forestsearch() returns differs in element set and meaning across the three identifier paths, and three early exits return an unclassed three-element list that the @return section does not describe. Also at: forestsearch_main.R:3002 (subgroup_method element, GRF return); :2821, :3012 (mr_inference element absent unless MR ran); :3368 with :3733 vs forestsearch_helpers.R:1666, 1989-1991 (df.predict); ... (full list in the CSV) | no | Build every exit from one constructor with a fixed element set and the forestsearch class, and document the failure shape in @return. | changes behaviour |
| A-170 | 1 | A | forestsearch_main.R | 3272 | forestsearch | S3 | on | traced | With use_dina = TRUE an invalid dina_args (for example an unrecognised key) gives a warning and a run without DINA cuts when details = FALSE but a hard error when details = TRUE, because the details printout repeats .resolve_dina_args() outside the tryCatch. Also at: forestsearch_main.R:3182-3185 (the same call inside the screening tryCatch), :3263-3268 (its handler); forestsearch_helpers.R:1101-1103 | no | Resolve dina_args once before the screening block and reuse the result for the details printout. | moves code |
| A-171 | 1 | A | forestsearch_main.R | 3554 | forestsearch | S3 | on | traced | The by.risk fallback computed when plot.sg = TRUE and by.risk = NULL is assigned to a local that nothing reads, so subgroup.consistency() still receives by.risk = NULL from args_call_all. Also at: forestsearch_main.R:3568-3592 (consistency_overrides omit by.risk), :3629-3633; subgroup_consistency_main.R:384, 1149; subgroup_consistency_helpers.R:1065, 974-975 | no | Pass the resolved by.risk in consistency_overrides (or sync it into args_call_all). | add-only |
| D-004 | 1 | D | forestsearch_main.R | 3654 | forestsearch | S3 | on | traced | An error inside subgroup.consistency() is turned into a warning and a normal classed return with sg.harm = NULL and no error_log, so unlike the get_FSdata and subgroup.search failures it cannot be told apart in the object from a legitimate 'no subgroup found'. Also at: forestsearch_main.R:3678-3683, :3942-3982 (normal return, no error marker); contrast :3327-3333 and :3501-3507 (error_log on the other two stage failures); bootstrap_analysis_dofuture.R:408, :626 | no | Attach the same error_log marker (stage = subgroup.consistency) to the returned object when the consistency stage errors. | add-only |
| A-172 | 1 | A | forestsearch_main.R | 3715 | forestsearch | S3 | on | traced | id.name is never validated up front: the 'unrelated id column' stop runs only after the full search and only when a subgroup was found, a missing id column fails late with a base-R subsetting error, and the DINA/GRF paths silently substitute other ids for the same input. Also at: forestsearch_main.R:3726-3732 (same guard for df.predict); :3390 (first use of id.name; a missing column fails here with base R's 'undefined columns selected'); ... (full list in the CSV) | no | Validate id.name (presence, and the 'id' collision rule) once in the input-validation section, identically for all three identifier paths. | moves code |
| A-173 | 1 | A | forestsearch_main.R | 3821 | forestsearch | S3 | on | traced | The MR candidate family is enumerated with only a non-strict size floor (size >= n.min), whereas the search rejects size <= n.min and also drops candidates by zero variance, prevalence (minp) and redundancy (rmin), so the family differs from the searched one in more ways than the two gaps the code comment and NEWS acknowledge. Also at: forestsearch_main.R:3802-3807 (comment: 'identical ... with two known gaps'); subgroup_search.R:610-612 (variance), :615-617 (minp), :620-623 (rmin redundancy), :656 and :706 (nx <= n.min rejects size exactly n.min); ... (full list in the CSV) | unknown | Either use the search's own strict size rule and pre-fit filters when enumerating the MR family, or extend the comment and the NEWS known-limitation entry to list these additional gaps. | method-adjacent |
| C-006 | 1 | C | fs_mr_inference.R | 590 | fs_mr_inference | S3 | on | traced | The exported fs_mr_inference() has a selection_rule argument with no @param entry. Also at: man/fs_mr_inference.Rd:17 | no | Add a @param selection_rule entry describing the three band rules and that only the effMaxSG / effMinSG re-selection rules consult it. | add-only |
| A-174 | 1 | A | fs_mr_inference.R | 665 | fs_mr_inference | S3 | on | traced | An admission list with a consistency floor but no effect floor falls into the 'unrestricted' branch, so the supplied consistency floor is not applied as an admission filter and no condition is raised. Also at: fs_mr_inference.R:678; fs_mr_inference.R:681-685 | no | Either apply the consistency floor on its own as a fourth branch or stop with an error for the consistency-only pattern, and state the choice in the @param. | changes behaviour |
| A-175 | 1 | A | get_FSdata_helpers.R | 519 | collapse_cuts_round | S3 | on | traced | A threshold in (-0.5, 0) rounds to negative zero and is formatted as "-0", so the coarsened cut is labelled 'x <= -0' and is not de-duplicated against 'x <= 0', leaving two candidate cuts with identical membership. Also at: get_FSdata_helpers.R:534; get_FSdata_helpers.R:581 | unknown | Normalise negative zero to zero in collapse_cuts_round() (or before formatting) so the label is 'x <= 0' and string de-duplication works. | add-only |
| C-007 | 1 | C | get_fsdata.R | 48 | get_FSdata | S3 | on | traced | Several roxygen statements for get_FSdata's cut arguments and return value, and the lasso_selection title, disagree with the code. Also at: get_fsdata.R:52; get_fsdata.R:62; get_fsdata.R:13-21; get_FSdata_helpers.R:214-216; forestsearch_main.R:578; forestsearch_main.R:600-608 | no | Bring the roxygen for exclude_cuts, outcome_type, the confs return entry, the collapse and rounding behaviour and the lasso_selection title into line with the code. | add-only |
| A-176 | 1 | A | get_fsdata.R | 328 | get_FSdata | S3 | on | traced | The LASSO cross-validation seed is fixed at lasso_selection()'s own default because get_FSdata() has no seed argument, so forestsearch(seedit = ...) does not reach the LASSO folds. Also at: get_FSdata_helpers.R:242 | no | Thread the caller's seed through get_FSdata() to lasso_selection(), or state in the seedit documentation that the LASSO fold seed is fixed. | add-only |
| A-177 | 1 | A | get_fsdata.R | 470 | get_FSdata | S3 | on | traced | conf.cont_medians_force is appended to conf.cont_medians twice, so with collapse_cuts = FALSE and no other de-duplicating branch the same median cut is returned as two identical candidate factors. Also at: get_fsdata.R:434-436 | no | Re-introduce the forced names only when they were removed (or de-duplicate confs unconditionally before evaluation). | changes behaviour |
| C-008 | 1 | C | glm_effect_estimators.R | 59 | make_effect_estimator | S3 | on | traced | The make_effect_estimator roxygen says adverse_outcome is ignored for continuous outcomes although the code negates Y there, omits the count measures from effect_measure, and omits the reason field from the return value. Also at: glm_effect_estimators.R:47; glm_effect_estimators.R:79; glm_effect_estimators.R:4; glm_effect_estimators.R:369; glm_effect_estimators.R:413 | no | Bring the @param adverse_outcome, @param effect_measure and @return text into line with the dispatcher. | add-only |
| A-178 | 1 | A | glm_effect_estimators.R | 655 | .estimate_rr | S3 | on | traced | The log-binomial start vector is hard-coded to length 2, so with adjust_covariates the log-binomial tier always errors and the adjusted RR is always the modified-Poisson fit. Also at: glm_effect_estimators.R:484 | no | Size the start vector from the model matrix so the log-binomial tier is attempted for adjusted fits. | changes behaviour |
| D-005 | 1 | D | glm_effect_estimators.R | 879 | .make_poisson_rate_estimator | S3 | on | traced | The converged field means different things across the sibling estimator closures, and the consumers disagree on whether to honour it. Also at: glm_effect_estimators.R:333; glm_effect_estimators.R:441; glm_effect_estimators.R:560; glm_effect_estimators.R:791; glm_effect_estimators.R:833; glm_effect_estimators.R:846; glm_effect_estimators.R:897; ... (full list in the CSV) | no | Define one meaning for converged across the estimator factories and make the consumers apply the same rule. | changes behaviour |
| A-179 | 1 | A | glm_effect_estimators.R | 888 | .make_poisson_rate_estimator | S3 | on | traced | The value returned as IRD is not one quantity: the unadjusted branch returns (rate1 - rate0) multiplied by the slice's geometric-mean follow-up (covariates, if any, set to zero), while the dr_gcomp branch and the package's own DGM truth use the arithmetic mean of expected counts at each subject's own follow-up, and none of them is a rate per unit time as the roxygen says. Also at: glm_effect_estimators.R:889-893 (unadjusted delta chain), :823-824 and :840 (dr_gcomp), :719 (roxygen); generate_glm_dgm.R:475-487, :676 (package truth for IRD); ... (full list in the CSV) | unknown | State the IRD definition in one place and make the unadjusted branch, the dr_gcomp branch and the documentation refer to that same quantity. | method-adjacent |
| A-180 | 1 | A | glm_effect_estimators.R | 1031 | .fs_lm_two_group | S3 | on | traced | The fast unadjusted closures are not equivalent to the slow ones when the outcome or treatment column of the slice has a missing value: the matrix fits error and the closure returns NA with a failure label, where lm()/glm() drop the incomplete rows and return an estimate, contrary to the roxygen statement that failure behaviour is identical. Also at: glm_effect_estimators.R:1149 (glm.fit core), :1061-1068 and :1164-1171 (roxygen claiming identical failure behaviour), :1090-1108, :1203-1213 (error handlers) | no | Drop incomplete rows in the two-group cores before fitting, or state and check a complete-data precondition and correct the roxygen claim. | changes behaviour |
| A-181 | 1 | A | glm_effect_estimators.R | 1250 | .robust_or_model_se | S3 | on | traced | With robust_se = TRUE (the default) any failure of the sandwich computation, including the Suggests package not being installed, silently yields the model-based SE with no message and no change to method_used, so the reported SE of the IRR and modified-Poisson RR fits depends on an optional package without trace. Also at: glm_effect_estimators.R:61-62 (roxygen), :688, :787, :875 (the only call sites); run_subgroup_sims.R:207-213 (documents the no-op scope) | no | Record which SE was used (for example in method_used or a message), restrict the handler to the missing-package case, and state in make_effect_estimator's roxygen which measures the flag applies to. | add-only |
| C-009 | 1 | C | glm_effect_estimators.R | 1346 | estimate_propensity_scores | S3 | on | traced | The roxygen describes trim as quantile bounds, ps_hat as the raw score and trimmed as observations trimmed, but the code clips the scores to the absolute bounds trim[1] and trim[2], stores the clipped vector, removes nothing, and the method = "none" return has no sw element although @return lists one. Also at: glm_effect_estimators.R:1334 (ps_hat called raw), :1357 (trimmed), :1395-1402 ("none" return), :1416-1422 (clipping code), :1450 | no | Reword the trim, ps_hat and trimmed documentation to describe clipping at absolute probability bounds and document, or align, the element set returned for method = "none". | add-only |
| A-182 | 1 | A | glm_effect_estimators.R | 1496 | .ps_lasso | S3 | on | traced | With ps_method = "lasso" the propensity model's cross-validation folds are drawn from the ambient RNG: the seed forestsearch() passes (seedit) reaches only the GRF branch, so the propensity scores are not reproducible from seedit and the caller's RNG stream is advanced. Also at: glm_effect_estimators.R:1408-1412 (seed passed to .ps_grf only), :1345 (roxygen: seed for GRF); forestsearch_main.R:2688-2694 (call site passing seed = seedit) | no | Derive the cv.glmnet fold assignment deterministically from the supplied seed without leaving the caller's RNG state changed. | changes behaviour |
| C-010 | 1 | C | grf_main.R | 57 | grf.subg.harm.survival | S3 | on | traced | The roxygen of grf.subg.harm.survival() misstates the length of sg.harm.id, the tree.cuts behaviour when no subgroup is found, and the fields of the frontier-mode return, and forestsearch() documents dmin.grf as 'Minimum events for GRF'. Also at: grf_main.R:62-64, 81, 99-112; forestsearch_main.R:555 | no | Correct the @return entries, the no-subgroup sentence and the dmin.grf description in forestsearch() to match the code. | add-only |
| A-183 | 1 | A | grf_main.R | 183 | grf.subg.harm.survival | S3 | on | traced | maxdepth is checked only for '> 3' on the survival path and not at all on the GLM path, and the survival helpers loop with 1:maxdepth, so 0 or negative values fail with an indexing error (survival) or return 'no subgroup' silently (GLM), and depths above 3 are accepted by the GLM sibling only. Also at: grf_helpers.R:163, 327, 384, 558; grf_subg_harm_glm.R:577-581; grf_helpers.R:36 | no | Validate maxdepth once as an integer in 1..3 in both identifiers and iterate with seq_len(). | add-only |
| C-011 | 1 | C | grf_subg_harm_glm.R | 22 | grf.subg.harm.glm | S3 | on | traced | The roxygen of grf.subg.harm.glm() describes a subgroup-enumeration and GLM effect-estimation stage via create_glm_row() that the function does not perform, documents inert arguments as active, and its @return does not match the returned object. Also at: grf_subg_harm_glm.R:30-32, 84-85, 92-99, 134-160, 170, 204-209, 215-219, 339-344 | no | Rewrite the description, details and @return to describe what the function does (causal forest, policy-tree or frontier selection, returned fields) and drop the estimation claims and the inert-argument descriptions. | add-only |
| A-184 | 1 | A | grf_subg_harm_glm.R | 582 | grf.subg.harm.glm | S3 | on | traced | Several GRF steps discard errors or warnings without a message, so a failed policy tree, candidate enumeration or cut extraction is indistinguishable from 'no subgroup' or 'no candidates'. Also at: grf_main.R:269-276; grf_subg_harm_glm.R:499-506, 929-939; grf_helpers.R:73, 83; forestsearch_helpers.R:1924 | no | Surface the caught conditions as a warning or a recorded message in the result instead of returning NULL silently. | add-only |
| D-006 | 1 | D | grf_subg_harm_glm.R | 685 | grf.subg.harm.glm | S3 | on | traced | grf.subg.harm.survival() and grf.subg.harm.glm() return differently shaped objects for the same roles (tree.cuts type and rounding, tree vs trees, selected_depth vs best_depth, null-result contents, class, rows of $data), so callers convert per type and plot.grf = TRUE draws nothing for GLM outcomes without a message. Also at: grf_helpers.R:246, 469, 480, 491, 540-555; grf_subg_harm_glm.R:623-643, 927-939; forestsearch_main.R:3130-3134, 3155; get_fsdata.R:258-261; forestsearch_helpers.R:1954-1957 | unknown | Give the two identifiers one return contract (same field names and types for cuts, trees, depth and the null result) and convert at the source instead of in each caller. | changes behaviour |
| C-012 | 1 | C | interpret_search_config.R | 147 | interpret_search_config | S3 | on | traced | Three passages of the printed search diagnostic describe behaviour the executing code does not have: 'GRF sees 1-Y' for a continuous adverse outcome, an unconditional 'GRF sees log(Y + 0.5)' for counts, and an RD/RR search direction that ignores the estimator's Y flip. Also at: interpret_search_config.R:226-227; interpret_search_config.R:172-176; interpret_search_config.R:112-115, 119-120; interpret_search_config.R:62-65 | no | Correct the three message texts (passing grf_count_transform through for the count note) and validate outcome_type at entry. | add-only |
| C-013 | 2 | C | mr_values.R | 42 | mr_values | S3 | on | traced | The roxygen for ij_lower / ij_upper calls them the two-term IJ interval and lists only ci_method = "wald" and ij_residual != "two_term" as exceptions, but when the IJ variance is unusable (ij_source = "wald_fallback") the returned limits are robust-SE limits and no returned column identifies that case. Also at: mr_values.R:136-144 (two_term_ci gate and its comment), :163-164; fs_mr_inference.R:220-225, :755, :770, :1070-1071 (H side); :846, :854, :861-862 (Hc side) | no | Document the wald_fallback case for ij_lower / ij_upper as is done for ij_se, or expose ij_source for H and Hc (or return NA limits in that case). | add-only |
| B-001 | 2 | B | mr_values.R | 241 | .mr_get | S3 | on | traced | The exact-matching path walker and the scalar-or-NA coercion for an mr_inference list are implemented three times (.mr_get/.mr_num, and local closures in .mrct_mr_record() and in mr_estimates_table.R) with only incidental differences. Also at: mr_values.R:252 (.mr_num); mrct_simulation.R:773-775 (num), :776-783 (g) in .mrct_mr_record; mr_estimates_table.R:234 (num) | no | Have .mrct_mr_record() and the mr_estimates_table.R reader call .mr_get() / .mr_num() instead of their local copies. | moves code |
| D-007 | 1 | D | subgroup_consistency_helpers.R | 571 | sort_subgroups | S3 | on | traced | Sibling sg_focus branches of the exported sorter have different contracts: most reorder the caller's data.table in place with setorder (NA keys first), while hrMaxSG/hrMinSG return a reordered copy built with order() (NA keys last). Also at: subgroup_consistency_helpers.R:583, :587, :591, :595, :605 (setorder) versus :620-626 (order + copy); sort_subgroups_preview: :665, :672, :676, :680, :684, :694 versus :722-728 | no | Give every branch the same contract: one sorting primitive, one NA placement, and either always copy or always sort in place, documented. | changes behaviour |
| F-002 | 1 | F | subgroup_consistency_helpers.R | 976 | plot_subgroup | S3 | on | traced | plot_subgroup() sets par(mfrow = c(1, 2)) without saving and restoring the previous graphical parameters, so forestsearch(plot.sg = TRUE) on a survival outcome leaves the device in a two-panel layout. | no | Capture the previous par() value and restore it with on.exit() in plot_subgroup(). | add-only |
| C-014 | 1 | C | subgroup_consistency_helpers.R | 1170 | .make_eval_subgroup_consistency | S3 | on | traced | Stale documentation in the file: an orphaned 'Evaluate Single Subgroup for Consistency (Fixed-Sample)' roxygen block is fused onto a @noRd factory, the sorters' sg_focus lists omit maxeff/maxeffCons and the aliases, a comment calls the internal preview sorter exported, and the pareto summary differs from the rule the code applies. Also at: subgroup_consistency_helpers.R:472-473, :494-508, :641-642 (sg_focus lists); :658-659 (comment); :461-463 (pareto wording) | no | Delete the orphan block and bring the sg_focus lists, the pareto sentence and the comment in line with the dispatch. | add-only |
| C-015 | 1 | C | subgroup_consistency_helpers.R | 1354 | evaluate_subgroup_consistency | S3 | on | traced | Both exported evaluators document a 'Named numeric vector' return but return a named character vector, because the numeric fields are combined with the character factor labels in one c() call. Also at: subgroup_consistency_helpers.R:1639 (two-stage @return), :1216 (orphan block); construction at :1433-1434, :1561-1570, :1905-1906, :2040-2049 | no | Document the character return, or return a typed one-row structure so numeric fields are not coerced through text. | changes behaviour |
| B-002 | 1 | B | subgroup_consistency_helpers.R | 1683 | evaluate_consistency_twostage | S3 | on | traced | evaluate_consistency_twostage() defines local copies of four package-level helpers (Wilson interval, early-stop decision, fast split HR, single split) that duplicate the exported versions and have already diverged in small ways, while a third Wilson helper with a different signature exists under the same dotted name in fpr_calibration.R. Also at: subgroup_consistency_helpers.R:97, :134, :192, :266, :1695, :1704, :1730; fpr_calibration.R:597 | no | Have the two-stage evaluator call the package-level helpers and keep a single Wilson / early-stop implementation. | moves code |
| C-016 | 1 | C | subgroup_consistency_main.R | 158 | subgroup.consistency | S3 | on | traced | The roxygen of subgroup.consistency() says early stopping is valid for sg_focus 'hr', 'maxSG' and 'minSG' and that batch_size defaults to 1 for 'hr'/'minSG', whereas the code treats early stopping as sound, and uses batch size 1, for 'maxeffCons' only, a value the sg_focus documentation does not list. Also at: subgroup_consistency_main.R:276-279 (batch_size default); subgroup_consistency_main.R:103-125 (sg_focus values omit maxeffCons) | no | Rewrite the stop_threshold and batch_size passages to match the maxeffCons-only logic and add maxeffCons to the documented sg_focus values. | add-only |
| C-017 | 1 | C | subgroup_consistency_main.R | 172 | subgroup.consistency | S3 | on | traced | show_candidate_summary is documented as printing two diagnostic tables when TRUE, but every print site also requires details = TRUE, a dependence that is not documented and that the roxygen examples do not satisfy. Also at: subgroup_consistency_main.R:683; subgroup_consistency_main.R:1065; subgroup_consistency_main.R:1177; subgroup_consistency_main.R:351; subgroup_consistency_main.R:361; forestsearch_main.R:871-875 | no | Document that the tables are printed only when details = TRUE (or remove the details condition). | add-only |
| A-185 | 1 | A | subgroup_consistency_main.R | 398 | subgroup.consistency | S3 | on | traced | The arguments checking and Lsg are accepted and documented but are never referenced in the body of subgroup.consistency(). Also at: subgroup_consistency_main.R:387 (Lsg formal); :101 and :205 (roxygen); forestsearch_main.R:3571 | no | Document both arguments as unused and retained for compatibility, or remove them through a deprecation step. | add-only |
| D-008 | 1 | D | subgroup_consistency_main.R | 780 | subgroup.consistency | S3 | on | traced | For sg_focus = 'maxeff' the behaviour described in this function's roxygen holds only when forestsearch() has applied its overrides: on a direct call the winner's consistency is reported as NA whenever it falls below pconsistency.threshold, the pool is still truncated by stop_Kgroups, and algorithm reports 'twostage' although the two-stage evaluator is never used. Also at: subgroup_consistency_main.R:648-664 (truncation and comment); subgroup_consistency_main.R:1200 (algorithm label; also :464, :542, :1084); roxygen :107-112 | no | Apply the maxeff overrides (threshold 0 for the winner's evaluation, no truncation, correct algorithm label) inside subgroup.consistency(), or document that they hold only through forestsearch(). | changes behaviour |
| A-186 | 1 | A | subgroup_consistency_main.R | 990 | subgroup.consistency | S3 | on | traced | The parallel branch wraps the whole per-batch candidate evaluation in suppressWarnings(), so every warning raised while evaluating candidates is discarded there, whereas the sequential branch lets the same warnings through. Also at: subgroup_consistency_main.R:887 | no | Muffle only the intended 'was built under R version' warning by message, as mrct_simulation.R:690-697 already does, instead of all warnings. | changes behaviour |
| A-187 | 1 | A | subgroup_consistency_main.R | 1140 | subgroup.consistency | S3 | on | traced | Any error in the selection step sg_consistency_out() is turned into a warning and the function returns an object with out_sg, sg.harm, sg.harm.id and df_flag all NULL while n_passed is positive, which downstream code reads as no subgroup identified. Also at: forestsearch_main.R:3654-3676 (same pattern one level up) | no | Let the error propagate, or return an explicit failure marker distinct from the no-subgroup result. | changes behaviour |
| A-188 | 1 | A | subgroup_search.R | 98 | subgroup.search | S3 | on | traced | The ID argument of subgroup.search() is documented as an optional vector of subject IDs but is never used, search_combinations_parallel() still takes max.minutes and t.start formals that its body never reads, and prop_max_count is a hard-coded 1 that forestsearch() still stores as prop_maxk. Also at: subgroup_search.R:14 (roxygen); subgroup_search.R:181; subgroup_search.R:274; subgroup_search.R:1238; subgroup_search.R:1266; forestsearch_main.R:3521 | no | Deprecate or document the unused ID argument and drop the dead internal pass-through parameters and the constant. | changes behaviour |
| A-189 | 1 | A | subgroup_search.R | 218 | subgroup.search | S3 | on | traced | The details filtering summary prints the effect threshold from hr.threshold rather than the threshold actually applied and exponentiates it for any ratio measure, so a direct survival call with effect_measure = 'HR' prints exp(1.25) = 3.490, and it prints a per-arm 'subjects' criterion for continuous and count outcomes where that filter is skipped. Also at: subgroup_search.R:211-213; subgroup_search.R:235-236 | no | Print the threshold that was applied (screen_threshold) and back-transform only when it is on the link scale. | changes behaviour |
| C-018 | 1 | C | subgroup_search.R | 470 | get_subgroup_membership | S3 | on | traced | get_subgroup_membership() is documented as returning a numeric 1/0 vector but returns a logical vector, except when no factor is selected, where it returns numeric ones, so the return type depends on the input. Also at: subgroup_search.R:454, 458 (roxygen); subgroup_search.R:468 (numeric branch) | no | Return the same type from both branches and state that type in @return. | changes behaviour |
| A-190 | 1 | A | subgroup_search.R | 549 | extract_idx_flagredundancy | S3 | on | traced | The redundancy screen compares each factor only with the running intersection of the factors before it in column order, so a nested pair whose larger set comes first (B inside A, A first) is kept as a two-factor candidate with exactly the membership of B alone, while the same pair in the other order is dropped. Also at: subgroup_search.R:540-554; fs_oc_family.R:505-515 (population mirror of the same ordered walk) | unknown | State in the function's roxygen that the screen is an ordered column walk and is therefore not symmetric in the factors; any symmetric variant would change the candidate family and is a method decision. | method-adjacent |
| C-019 | 1 | C | subgroup_search.R | 656 | evaluate_combination_with_status | S3 | on | traced | The size, redundancy and effect screens are strict (n > n.min, reduction > rmin, effect > threshold) while the roxygen and the printed filter summary call them minimums and print '>=', so a subgroup of exactly n.min subjects is rejected although the text says it is admissible. Also at: subgroup_search.R:706 (survival path, same test); labels and roxygen for size: subgroup_search.R:16, 214, 237, 566; forestsearch_main.R:614; redundancy: subgroup_search.R:549 vs :46, :533, :576; ... (full list in the CSV) | no | Reword the roxygen and the printed filter labels to the strict comparisons the code applies. | add-only |
| A-191 | 1 | A | subgroup_search.R | 1060 | fit_glm_for_subgroup | S3 | on | traced | Every per-candidate fit is wrapped in try() or tryCatch(error = NULL), so an error that is the same for every candidate (for example an adjustment column absent from df_analysis in a direct subgroup.search() call) is counted as 'failed model fit' for all of them and the search returns no candidates with no error or warning. Also at: subgroup_search.R:960-982, 1004-1010 (Cox try()); subgroup_search.R:1132-1135, 1174-1175 (fast GLM paths); subgroup_search.R:992-999 (missing adjustment column dropped from the frame but kept in the formula) | no | Keep the first fit error message, or warn once when every candidate that reached the fit stage failed, so a systematic failure can be told apart from an empty candidate set. | add-only |
| A-192 | 2 | A | summary_utility_functions.R | 454 | rmst_calculation | S3 | on | traced | rmst_calculation() returns its RMST values through signif(x, 2), so the RMST column of every survival estimates table (and m1/m0 when return_medians = FALSE) carries only two significant digits and the precision cannot be recovered by the caller. Also at: summary_utility_functions.R:556-561 (analyze_subgroup places the rounded values in the table row); summary_utility_functions.R:75 (medians rounded to one decimal in the same row) | yes | Return unrounded RMST values from rmst_calculation() and round to a fixed number of decimals where the table row is formatted. | changes behaviour |
| C-020 | 2 | C | summary_utility_functions.R | 587 | SG_tab_estimates | S3 | on | traced | Several roxygen statements for SG_tab_estimates() and sg_tables() disagree with the code: the return class, two arguments described as doing something they do not do, the shape of the sg_tables() return, and the reach of ndecimals. Also at: summary_utility_functions.R:583-584 (draws, details); 726-727 and 736-737 (sg_tables include_search_info and @return); 719-725 (ndecimals) | no | Correct the @return and @param text to describe the objects and behaviour the code produces. | add-only |
| D-009 | 2 | D | summary_utility_functions.R | 604 | SG_tab_estimates | S3 | on | traced | The sibling table builders SG_tab_estimates() and SG_tab_estimates_glm() differ in return class, in the argument names for the bias-corrected estimates, in the subgroup-label defaults and in whether precision is an argument. Also at: summary_utility_functions.R:328-337 (SG_tab_estimates_glm formals); 349-362, 391-400 (data-frame return); 617-637 (matrix / vector return); 764-765 (sg_tables exposes only hr_1a / hr_0a) | no | Align the two builders on one return class, one naming scheme for the bias-corrected estimates and the same label defaults. | changes behaviour |
| A-193 | 1 | A | grf_helpers.R | 77 | fit_causal_forest | S3 | on | suspected | The forests are trained with num.threads = NULL (grf's default, all hardware threads) outside R CMD check and with 2 threads under it, and no user-facing argument fixes the thread count, so for a fixed seedit the fitted forest and the GRF-derived cuts may differ between machines and between check and interactive runs. Also at: grf_helpers.R:87; grf_subg_harm_glm.R:465, :469; glm_effect_estimators.R:1480; utils_check_threads.R:22-27 | unknown | Expose and record a thread count for the GRF fits (the helper already accepts a user value) so that a seed identifies one forest regardless of machine or check environment. | changes behaviour |
| A-194 | 1 | A | subgroup_consistency_helpers.R | 49 | setup_parallel_SGcons | S3 | on | suspected | An NA from parallel::detectCores() is not guarded where the worker count is capped or defaulted, so on a platform where the core count is unknown the workers value passed to future::plan() becomes NA. Also at: forestsearch_main.R:2063-2076 (same min(), result not used afterwards); forestsearch_cross_validation.R:1928 and :1950-1951 (max(1, NA - 1) is NA; 'if (workers > NA)' errors); mrct_simulation.R:361; ... (full list in the CSV) | no | Route every detectCores() use through one helper that maps NA to a safe count, as .default_parallel_workers() already does. | add-only |
| A-195 | 1 | A | subgroup_consistency_helpers.R | 213 | get_split_hr_fast | S3 | on | suspected | On the literal-split path a half-sample with no events in one arm is fitted with warnings suppressed and its divergent exp(beta) is scored as a valid split, whereas the search treats a zero-event arm as a fit failure. Also at: subgroup_consistency_helpers.R:312 (only a total-events check per half); embedded copy :1720-1727, :1763, :1768-1769 | no | Treat a half without events in both arms (or a non-converged fit) as an invalid split rather than a scored one; this changes the split statistic, so it is a method decision. | method-adjacent |
| INV-2b | 1 | INV | subgroup_consistency_main.R | 883 | subgroup.consistency | S3 | on | suspected | The plan is saved with future::plan() called without arguments, which is understood to return only the top strategy, so a caller's nested plan stack would be restored as a single-level plan. Also at: forestsearch_cross_validation.R:332-333; forestsearch_cross_validation.R:840-841; run_subgroup_sims.R:817-818; run_subgroup_sims.R:829-830 | no | Restore from the value returned by the plan-setting call (or from plan('list')) so the whole stack is reinstated. | changes behaviour |
| A-196 | 1 | A | subgroup_search.R | 734 | evaluate_combination_with_status | S3 | on | suspected | On the survival path the fitted HR is compared with the threshold with no NA guard, so if the adjusted Cox fit returns an NA treatment coefficient for one candidate, `if (NA)` aborts the whole search, whereas the GLM paths turn an NA estimate into status 5. Also at: subgroup_search.R:1004-1031 (adjusted Cox arm returns the row without a finiteness check); contrast subgroup_search.R:1071-1076, 1137, 1181 (GLM paths map an NA estimate to a failed fit) | no | Treat a non-finite fitted HR as a failed fit (status 5) before the threshold comparison. | add-only |
| A-197 | 1 | A | subgroup_search.R | 988 | fit_cox_for_subgroup | S3 | on | suspected | The adjusted search-stage frame uses the internal column names Y, E, Treat and id.x, but forestsearch() reserves only Y, Event, Treat and id, so an adjustment covariate named E (or id.x) is bound beside the internal column of the same name and the formula is expected to resolve to the internal one. Also at: forestsearch_main.R:3419 (reserved-name list) | no | Extend the reserved-name check to the names the search-stage frame actually uses (E, id.x), or build that frame with collision-proof internal names. | add-only |
| C-021 | 2 | C | betaHhat_truth.R | 574 | fs_attach_betaHhat | S3 | off | traced | The roxygen says an undetected replicate always receives nH_eval = 0, nHc_eval = nrow(frame), betaHhat_H = NA and betaHhat_Hc = the ITT effect, but under focus = "benefit" the code assigns the reverse record. Also at: betaHhat_truth.R:310; betaHhat_truth.R:623 | no | State the benefit-focus no-subgroup record in the roxygen, or build the no-subgroup record without the focus inversion, whichever the author intends. | add-only |
| E-004 | 2 | E | bootstrap_calculations_helpers.R | 103 | bootstrap_ystar and 116 others | S3 | off | traced | 117 exported off-path functions are never referenced under tests/testthat/ (static whole-word grep), for example bootstrap_ystar(), count_boot_id(), format_bootstrap_table(), summarize_bootstrap_events(); the full list is in section 8. | no | Add tests that call these functions directly, starting with the on-path ones that carry selection or estimation logic. | add-only |
| B-003 | 2 | B | bootstrap_dofuture_setup_helpers.R | 42 | get_bootstrap_exports | S3 | off | traced | get_bootstrap_exports() is unexported and has no caller anywhere; the bootstrap loop relies on automatic globals detection, yet the roxygen of bootstrap_results() still tells readers the workers need the functions it lists. Also at: bootstrap_analysis_dofuture.R:238; bootstrap_analysis_dofuture.R:439; bootstrap_dofuture_main.R:285; bootstrap_calculations_helpers.R:259 | no | Remove the unused helper, its man page and the roxygen note that refers to it. | moves code |
| C-022 | 2 | C | bootstrap_summaries_helpers.R | 33 | format_bootstrap_table | S3 | off | traced | The format_bootstrap_table() example passes summaries$table_data, an element summarize_bootstrap_results() never returns, and omits the required nb_boots; four further roxygen statements in the same files disagree with the code. Also at: bootstrap_summaries_helpers.R:548, 692 (effect_label formal, no @param; absent from \arguments in man/format_bootstrap_diagnostics_table.Rd and man/create_bootstrap_diagnostic_plots.Rd); ... (full list in the CSV) | no | Correct the example (pass fs_bc$FSsg_tab and nb_boots, or show summaries$table) and bring the four roxygen statements in line with the code. | add-only |
| A-198 | 2 | A | bootstrap_summaries_helpers.R | 38 | format_bootstrap_table | S3 | off | traced | The documented est.scale argument of format_bootstrap_table() is accepted but never read, and summarize_bootstrap_results() only forwards it, so supplying "1/hr" changes nothing in either function. Also at: summarize_bootstrap_results.R:65 (formal), :172 (only use: forwarded); callers passing it: mr_estimates_table.R:156, 306 | no | Either use est.scale to label the effect columns and footnotes, or remove (deprecate) the argument in both functions and their callers. | changes behaviour |
| C-023 | 2 | C | calibrate_glm_interaction.R | 49 | calibrate_glm_interaction | S3 | off | traced | grid_step is documented as the grid resolution but is never used, and the description and @details still describe a grid search although the code solves by root-finding. Also at: calibrate_glm_interaction.R:12-13; calibrate_glm_interaction.R:61-63; calibrate_glm_interaction.R:104 | no | Rewrite the description and @details for root-finding with the tol_rel stop, and document grid_step as ignored (kept for backward compatibility). | add-only |
| A-199 | 2 | A | calibrate_k_treat.R | 120 | calibrate_k_treat | S3 | off | traced | The 'uniroot objective may be noisy' warning and the matching @details sentence apply to every call whose base_args simply omits seed, although generate_aft_dgm_flex() then runs with its fixed default seed and the objective is deterministic. Also at: calibrate_k_treat.R:44-50 (@details) | no | Warn only when base_args carries an explicit NULL seed, and correct the @details sentence. | changes behaviour |
| A-200 | 2 | A | compare_fs_subgroups.R | 113 | compare_fs_subgroups | S3 | off | traced | When id.name is supplied but the column is missing from a result object, the function silently pairs rows by position instead of reporting that the requested ID merge was not done. Also at: compare_fs_subgroups.R:125 (GRF branch); compare_fs_subgroups.R:160 and 172-181 (silent positional alignment) | no | Warn or stop when id.name is given but is not found in a result object. | add-only |
| A-201 | 2 | A | cox_ahr_cde_wrapper.R | 163 | cox_ahr_cde_analysis | S3 | off | traced | cox_ahr_cde_analysis() cannot set cox_cs_fit()'s grid step, so with the fixed z_by = 1 a biomarker spanning less than one unit up to its 90th percentile gives a one-point profile and the spline-estimate line and band are not drawn. Also at: cox_spline_fit.R:116, 208-210 | no | Expose the grid arguments in the wrapper (or derive the step from the biomarker range) and forward them to cox_cs_fit(). | add-only |
| A-202 | 2 | A | cox_ahr_cde_wrapper.R | 384 | cox_ahr_cde_analysis | S3 | off | traced | The profile panel's legend is hard-coded to "95% CI" although the band is drawn at level 1 - alpha, which is 80% under the function's default alpha = 0.20. Also at: cox_ahr_cde_wrapper.R:472, 581 (legend drawn) | no | Build the legend label from alpha, as cox_cs_fit() does. | changes behaviour |
| C-024 | 2 | C | cox_ahr_cde_wrapper.R | 776 | print.cox_ahr_cde | S3 | off | traced | The \dontrun examples of print.cox_ahr_cde and summary.cox_ahr_cde call cox_ahr_cde_analysis() with arguments it does not have (outcome.name, event.name, treat.name, confounders.name) on survival::gbsg, which also lacks the required loghr_po column. Also at: cox_ahr_cde_wrapper.R:835-837 | no | Replace both examples with a call using the function's real argument names on data that carries the required columns (as the main function's own example at 72-99 does). | add-only |
| C-025 | 2 | C | cox_spline_fit.R | 61 | cox_cs_fit | S3 | off | traced | cox_cs_fit() documents the plot parameter show_target (target effect line) but the only code that reads it is commented out, so it has no effect; the @return lists of cox_cs_fit() and glm_cs_fit() also omit elements the functions return. Also at: cox_spline_fit.R:375-376 (defaults set); cox_spline_fit.R:441-443, 475-480 (commented-out readers); cox_spline_fit.R:63-74 vs 514-518 (@return omits alpha, ci_level, class); ... (full list in the CSV) | no | Remove show_target from the Plot Parameters documentation (or restore the line) and complete the two @return lists. | add-only |
| F-003 | 2 | F | cox_spline_fit.R | 128 | cox_cs_fit | S3 | off | traced | cox_cs_fit() and glm_cs_fit() call par(no.readonly = TRUE) as their first statement, before show_plot is consulted, so a call that draws nothing still opens a graphics device (Rplots.pdf in the working directory in a non-interactive session). Also at: glm_spline_fit.R:154-155; glm_spline_fit.R:154 | no | Capture and restore par only inside the show_plot branch of both functions. | moves code |
| D-010 | 2 | D | cox_spline_fit.R | 643 | plot_subgroup_effects | S3 | off | traced | When no subject meets the threshold, plot_subgroup_effects() falls back to the cutpoint min(z) (everyone in the 'optimal' subgroup) while the sibling cox_ahr_cde_analysis() falls back to max(z) (only the top value in 'Recommend'), for the same subgroup definition z >= cutpoint. Also at: cox_ahr_cde_wrapper.R:198 | no | Use one fallback convention in both functions, preferably an explicit 'no subgroup' state. | changes behaviour |
| A-203 | 2 | A | create_summary_table.R | 365 | create_summary_table | S3 | off | traced | For a categorical variable with three or more levels and any cell count below 5, the p-value comes from fisher.test(simulate.p.value = TRUE), a Monte Carlo test drawn from the session RNG with no seed handling, so the printed p-value varies between calls, the caller's RNG stream is advanced, and the footnote still calls the test exact. Also at: create_summary_table.R:830-833 (footnote says 'Fisher's exact') | yes | Compute the exact test where feasible, or make the simulated p-value reproducible and say in the footnote that it is simulated. | changes behaviour |
| C-026 | 2 | C | dina_subgroup_bootstrap.R | 21 | dina_subgroup_bootstrap | S3 | off | traced | The roxygen title and opening sentence call the function a 'Selection-adjusted bootstrap CI' that 'performs selection-adjusted inference', while the same roxygen block, the file header and the print method state that neither CI adjusts for selection. Also at: dina_subgroup_bootstrap.R:23; man/dina_subgroup_bootstrap.Rd:5,191 | no | Reword the title and first sentence to describe conditional-on-signature bootstrap CIs plus selection-stability diagnostics. | add-only |
| C-027 | 2 | C | explain_pareto_selection.R | 62 | explain_pareto_selection | S3 | off | traced | The roxygen says format = "character" returns elements with no markdown markers and that verbosity = "detailed" adds resolved indicator values and a downstream-tie note, but the character output keeps the ** and '- ' markers and the closing paragraph holds only the lex-key text and candidate count. Also at: explain_pareto_selection.R:20-25; explain_pareto_selection.R:26-28 | no | Correct the @return and Verbosity text to describe what is produced (or strip the markers for format = "character"). | add-only |
| C-028 | 2 | C | find_k_inter_main.R | 35 | find_k_inter_for_target_hr | S3 | off | traced | @return and @details of find_k_inter_for_target_hr() describe a classed object with a `method` element and boundary diagnostics on failure, none of which the code produces. Also at: find_k_inter_main.R:43; find_k_inter_main.R:50-52; calibrate_k_treat.R:24-27; calibrate_glm_interaction.R:61-63 | no | Align the roxygen with the returned unclassed list and the stop-on-failure behaviour of the shared root-finder. | add-only |
| D-011 | 2 | D | find_k_inter_main.R | 88 | find_k_inter_for_target_hr | S3 | off | traced | The four effect calibrators that share .calibrate_by_root() expose different defaults for the search interval, tolerance and verbosity and return three different kinds of object. Also at: sim_aft_gbsg.R:1004-1045; calibrate_k_treat.R:96-153; calibrate_glm_interaction.R:89-218 | no | Agree one convention for the interval, tol and verbose defaults and for the returned object across the calibrator family. | changes behaviour |
| A-204 | 2 | A | find_k_inter_main.R | 405 | find_quantile_for_proportion | S3 | off | traced | The returned actual_proportion is the requested target_prop, not the proportion achieved at the returned cutpoint. | no | Return the proportion recomputed at the selected cutpoint, or remove the function. | changes behaviour |
| F-004 | 2 | F | forestsearch_cross_validation.R | 377 | forestsearch_Kfold | S3 | off | traced | Progress and diagnostic output in the cross-validation functions is written with cat()/print() rather than message(), and forestsearch_tenfold() emits it on a default call because details defaults to TRUE. Also at: forestsearch_cross_validation.R:586-588, 849-852, 1236-1239, 1386-1393, 1427, 1467, 1537, 1558, 1580, 1638, 1722, 2060-2103 (print_cv_params, called at 451 and 944); bootstrap_dofuture_main.R:391; ... (full list in the CSV) | no | Route non-print-method progress output through message() so it can be suppressed uniformly. | changes behaviour |
| A-205 | 2 | A | forestsearch_cross_validation.R | 606 | forestsearch_Kfold | S3 | off | traced | An error raised by forestsearch_KfoldOut() inside forestsearch_Kfold() is swallowed by tryCatch and replaced by all-NA sens_summary/find_summary, with a warning only when details = TRUE. | no | Signal the warning unconditionally, or let the error propagate as forestsearch_tenfold() does. | changes behaviour |
| D-012 | 2 | D | forestsearch_cross_validation.R | 832 | forestsearch_tenfold | S3 | off | traced | forestsearch_tenfold() does not check Kfolds against the number of observations as forestsearch_Kfold() does, so Kfolds > n fails late inside the simulation loop with a data-frame error that does not mention Kfolds. | no | Add the same Kfolds <= number-of-observations check to forestsearch_tenfold()'s input validation. | add-only |
| B-004 | 2 | B | forestsearch_cross_validation.R | 1191 | forestsearch_tenfold | S3 | off | traced | The partial-failure handling in forestsearch_tenfold() is unreachable now that any refit error is fatal, and two requireNamespace() fallbacks in summarize_cv_results.R guard packages that are hard Imports. Also at: forestsearch_cross_validation.R:1089, 1099, 1193-1199, 1238, 1254; summarize_cv_results.R:341-348, 1271-1274 | no | Remove the unreachable partial-failure handling and redundant try-error re-tests (or restore a deliberate partial-failure policy), and drop the Imports guards. | moves code |
| C-029 | 2 | C | forestsearch_cross_validation.R | 1278 | forestsearch_KfoldOut | S3 | off | traced | Several roxygen statements in the CV files disagree with the code: the digits headroom claim for cv_summary_tables(), elements and defaults documented as required that are optional or computed differently, and error/dropped-simulation outcomes that can no longer occur. Also at: forestsearch_cross_validation.R:188-190, 710-713, 739-741, 754-755, 769-770, 1272-1274; cv_summary_tables.R:33, 311-312, 637 | no | Correct the listed roxygen statements to match the code. | add-only |
| A-206 | 2 | A | forestsearch_cross_validation.R | 1402 | forestsearch_KfoldOut | S3 | off | traced | For a single-factor original subgroup the 'Cov 1 & 2' find metric evaluates to 0 or NaN depending on the data, instead of being consistently missing like 'Cov2'. Also at: forestsearch_cross_validation.R:1391 | yes | Return NA for 'Cov 1 & 2' whenever the original subgroup has fewer than two factors. | changes behaviour |
| C-030 | 2 | C | forestsearch_helpers.R | 753 | reset_workers | S3 | off | traced | reset_workers() roxygen says a failed strategy leaves the previous plan in place and that workers = NULL keeps the current plan's worker count, but the code tears the plan down to sequential first, leaves it sequential on failure, and starts the new plan with future's own default worker count. Also at: forestsearch_helpers.R:757-760 | no | Correct the two roxygen statements to describe what the code does (or save the previous plan and worker count before teardown and reuse them). | add-only |
| C-031 | 2 | C | forestsearch_methods.R | 222 | print.forestsearch | S3 | off | traced | The print.forestsearch roxygen says that with ci_method = "ij" only the two-sided interval is shown, but the code prints no post-selection block at all when the MR result has no field element. Also at: forestsearch_methods.R:85-86 (code), 310-311, 521-522 (callers); man/print.forestsearch.Rd:35 | no | Either print the two-sided IJ interval when the field block is absent or change the roxygen sentence to say nothing is printed. | add-only |
| A-207 | 2 | A | forestsearch_methods.R | 384 | summary.forestsearch | S3 | off | traced | summary.forestsearch() (and the bootstrap/CV parameter echoes) print the raw call values from args_call_all as 'Analysis Parameters', which differ from the values the fit used after the RD/IRD/MD threshold remap or the sg_focus = 'maxeff' overrides. Also at: forestsearch_methods.R:385-390, 394 (hr.consistency, pconsistency.threshold, use_twostage echoed from args_call_all); bootstrap_dofuture_main.R:415-436; forestsearch_cross_validation.R:2064-2072; ... (full list in the CSV) | no | Have summary()/the bootstrap and CV echoes print the resolved values (threshold_config, admission, resolved gating settings) or label the args_call_all values as the call as supplied. | add-only |
| A-208 | 2 | A | forestsearch_methods.R | 471 | summary.forestsearch | S3 | off | traced | Under a default call summary() prints 'Stage 1 screening splits:' and 'Stage 2 batch size:' with no values, because it reads the user-supplied twostage_args (an empty list by default) rather than the resolved two-stage settings. Also at: forestsearch_methods.R:472-475 (two cat() lines); forestsearch_methods.R:431-452 ('Search Space' block prints NA for DINA/GRF fits: prop_maxk / max_sg_est are NA_real_ at forestsearch_main.R:2804-2805, 2993-2994) | no | Print the resolved two-stage settings (or skip a line whose value is absent) and skip the search-space lines when the values are NA. | changes behaviour |
| D-013 | 2 | D | forestsearch_methods.R | 601 | plot.forestsearch | S3 | off | traced | plot.forestsearch() defaults the column names to "Y"/"Event"/"Treat" while forestsearch() defaults to "tte"/"event"/"treat" and stores the names it used on the object, so plot(fs) stops with a missing-columns error unless the caller repeats the names. Also at: forestsearch_methods.R:602-603; plot_sg_results.R:134-136, 174-179; forestsearch_main.R:1623-1625 | no | Resolve the three names from x$args_call_all when the caller does not supply them. | changes behaviour |
| B-005 | 2 | B | format_subgroup_summary_tables.R | 12 | format_subgroup_summary_tables | S3 | off | traced | Six unexported functions have no caller in R/, tests/, vignettes/ or inst/: format_subgroup_summary_tables(), its only callee create_factor_summary_tables(), and the four create_summary_table_* presets. Also at: format_subgroup_summary_tables.R:419; create_summary_table.R:904, 920, 939, 955; _pkgdown.yml:159-160, 245-248; Cox_estimation_helpers.R:265; ahr_analytic.R:735; bootstrap_summaries_helpers.R:821; ... (full list in the CSV) | no | Remove the six functions together with their man pages and pkgdown entries, or wire the formatter into the bootstrap summary path, updating the quarto development document that calls it through :::. | moves code |
| C-032 | 2 | C | fpr_calibration.R | 118 | fpr_calibration | S3 | off | traced | Several @return / @param statements disagree with the code: fpr_calibration's nsims is documented as replicates attempted but holds replicates completed, five returned fields are undocumented, and there is similar drift in fpr_calibration_transfer(), fpr_approx() and generate_detection_curve(). Also at: fpr_calibration.R:112-129 (vs :372-377); fpr_calibration.R:463-465 (vs :497, :509-512); fpr_approximation.R:61-77 (vs :201); truefind_asymptotic.R:336-337 (vs :379); truefind_asymptotic.R:341-346 (vs :406) | no | Bring the roxygen for these items in line with the returned objects and the argument handling. | add-only |
| C-033 | 2 | C | frontier_cis.R | 151 | compute_frontier_cis | S3 | off | traced | compute_frontier_cis() always refits treatment-only, unweighted models regardless of the adjust_covariates / propensity settings stored on the fit, and the roxygen does not say that the intervals are unadjusted. Also at: frontier_cis.R:228-232; roxygen frontier_cis.R:11-14 | no | State in the roxygen that the frontier intervals are unadjusted, and optionally warn when the fit used covariate or propensity adjustment. | add-only |
| A-209 | 2 | A | frontier_member_flags.R | 91 | frontier_member_flags | S3 | off | traced | When a frontier member's cut cannot be evaluated on fs$df.est the member's column is returned as all zeros with no warning or marker, whereas the sibling compute_frontier_cis() warns and returns NA for the same condition. Also at: frontier_member_flags.R:82-87; sibling frontier_cis.R:204-209 | no | Warn as compute_frontier_cis() does and mark the failed column (NA or a reported index), or at least document the all-zero convention in @return. | changes behaviour |
| C-034 | 2 | C | fs_bias_coverage.R | 37 | fs_sim_bias_coverage | S3 | off | traced | The roxygen tells the caller to pass a bundle's `readRDS(...)$results` and offers scale = "identity" for MD/RD, but never states that `betaHhat_<block>` must already be on the estimator columns' oriented scale, although the function compares the two columns as stored and the sibling fs_mr_oc_summary() applies the orientation bridge itself. Also at: fs_bias_coverage.R:107 (target read as stored); :155 (bias); :159 (two-sided coverage); :163-164 (one-sided coverage); contrast fs_mr_oc_summary.R:114,130,134 (sibling applies the orientation bridge itself) | no | State in the roxygen that betaHhat_<block> must be on the estimators' oriented scale (multiply by -1 for adverse_outcome = FALSE bundles), or add an orientation argument applying the same bridge as fs_mr_oc_summary(). | add-only |
| A-210 | 2 | A | fs_bias_coverage.R | 215 | fs_plot_bias_coverage | S3 | off | traced | A single NA/NaN in tbl$b, which fs_sim_bias_coverage() itself returns for a cell with one detected replicate, makes fs_plot_bias_coverage() stop inside seq() with an unrelated message instead of dropping that point. Also at: fs_bias_coverage.R:156, 167 | no | Build the curve grid from the finite values of tbl$b (and let the point layers drop the non-finite rows). | add-only |
| INV-2a | 2 | INV | fs_dgm_feasibility.R | 165 | fs_dgm_feasibility | S3 | off | traced | The RNG restore on exit runs only when .Random.seed existed at entry, so in a session with no seed yet the function leaves its own seeded stream in the global environment although the roxygen says the caller's RNG state is restored. Also at: fs_dgm_feasibility.R:162-166 | no | When no .Random.seed existed at entry, remove it from the global environment on exit. | add-only |
| A-211 | 2 | A | fs_dgm_feasibility.R | 205 | fs_dgm_feasibility | S3 | off | traced | For a count DGM with effect measure IRR, share_nonestimable is always 0 and the print method says no existence condition applies, although the roxygen states that IRR needs an event in each arm under the estimator boundary's condition and the count estimator does apply that condition. Also at: fs_dgm_feasibility.R:37-41 (roxygen); fs_dgm_feasibility.R:159, 193-196; fs_dgm_feasibility.R:260-264 (print note) | no | Either apply the IRR existence condition on count draws (subjects with outcome > 0 per arm) or correct the roxygen and print note to say it is not evaluated for count outcomes. | changes behaviour |
| A-212 | 2 | A | fs_dgm_scale.R | 372 | .fs_jensen_factor | S3 | off | traced | The Jensen factor documented as exact and conditional on both arms being non-empty keeps the probability mass of a region of size 1 in its normalisation while giving it a zero contribution, and conditions the arm split separately within each region size, so it is not the stated conditional expectation. Also at: fs_dgm_scale.R:374-379; roxygen fs_dgm_scale.R:212-216; internal doc fs_dgm_scale.R:356-364 | no | Make the conditioning event the same in the normalisation and in the contributions, or reword the documentation to describe the quantity actually computed. | method-adjacent |
| A-213 | 2 | A | fs_family_report.R | 151 | fs_family_report | S3 | off | traced | On the argument-list route the effect-screen and consistency-screen rows print hr.threshold / hr.consistency (survival) or effect.threshold / consistency.threshold (other outcomes) as read independently from the list or the formal default, without forestsearch()'s alias merge or its c2 = 0.80 * c1 derivation, so the values column can show thresholds the engine would not use. Also at: fs_family_report.R:152; fs_family_report.R:312; fs_family_report.R:336 | no | Apply the same alias merge and c2 derivation as forestsearch() before formatting the two threshold rows on the list route. | changes behaviour |
| A-214 | 2 | A | fs_fdr_report.R | 218 | fs_fdr_report | S3 | off | traced | The default n_workers is computed from parallel::detectCores(logical = FALSE) with no guard for an NA return and no check-environment cap, unlike the package's own .default_parallel_workers(). Also at: fs_fdr_report.R:399 | no | Derive the default from .default_parallel_workers() or apply the same NA guard and check-time cap. | changes behaviour |
| A-215 | 2 | A | fs_fdr_report.R | 239 | fs_fdr_report | S3 | off | traced | fs_fdr_report() forces mr_inference = TRUE on every replicate while inheriting use_lasso / use_grf / use_dina and the DINA / GRF selector settings verbatim, so for a fit made with a model-based front end (or a misaligned selector) every replicate stops at the MR configuration guard and even the raw FDR A, which needs no MR, is returned as NA. Also at: forestsearch_main.R:1894; forestsearch_helpers.R:2121 | no | Check the inherited configuration once before the replicate loop and either stop with a clear message or compute A without MR when the configuration cannot support MR. | changes behaviour |
| B-006 | 2 | B | fs_fdr_report.R | 240 | fs_fdr_report | S3 | off | traced | The fallback MR settings (ci_method = 'ij', draws = 1000, include_complement = TRUE) are guarded by is.null(fs_params$mr_inference_args), but a fitted object stores the formal default list(), so the fallback never applies to an ordinary fit and the replicates run with forestsearch()'s own MR defaults instead. Also at: fs_fdr_report.R:241-242; fs_fdr_report.R:323 | unknown | Either treat an empty list like NULL so the intended fallback applies, or delete the fallback and document that forestsearch()'s MR defaults are used. | changes behaviour |
| C-035 | 2 | C | fs_identification_summary.R | 542 | plot_fs_rule_pairs | S3 | off | traced | plot_fs_rule_pairs() is documented as accepting an fs_identification object and taking true_covariates from it, but the code unconditionally stops for that class, and the true_covariates assignment placed before the stop can never be used. Also at: fs_identification_summary.R:544-545 (roxygen); fs_identification_summary.R:556-561 (code); fs_identification_summary.R:557-558 (dead assignment) | no | Make the roxygen say an fs_identification object is rejected and drop the unreachable assignment, or retain sg_def in the object so the documented behaviour is real. | add-only |
| B-007 | 2 | B | fs_mr_inference.R | 42 | .fs_mr_confirm_null | S3 | off | traced | .fs_mr_confirm_null() has no caller in R/ and is not exported; the t_confirm default it describes is computed inline by a separate rule, so the unit test that pins the helper does not test the default the package uses. Also at: fs_mr_inference.R:630; tests/testthat/test-mr-inference.R:130-138; man/dot-fs_mr_confirm_null.Rd | no | Remove the helper with its test, or make fs_mr_inference() call it so the default has one definition. | moves code |
| A-216 | 2 | A | fs_mr_oc_summary.R | 131 | fs_mr_oc_summary | S3 | off | traced | The H-block fields are read with `$` (truth$effect_Q, truth$marg_H, rd$betaHhat_H), and each name is a strict prefix of its Hc sibling, so a payload that lacks the H field but carries the Hc one is scored against the complement's value by partial matching instead of yielding NA. Also at: fs_mr_oc_summary.R:130; fs_mr_oc_summary.R:132; fs_mr_oc_summary.R:185; fs_mr_oc_summary.R:187; fs_mr_oc_summary.R:298-300 (rd$n_harm, rd$sens, rd$spec, rd$ppv, rd$npv with $ against exact res[[col]] at :287) | no | Read payload fields with exact matching ([[ ]]) so an absent H field yields NA rather than the Hc value. | changes behaviour |
| A-217 | 2 | A | generate_aft_dgm_helpers.R | 391 | prepare_censoring_model | S3 | off | traced | With select_censoring = TRUE and cens_type = "uniform", cens_params is neither ignored as documented nor honoured consistently: min and max are used when both are supplied, and both are silently replaced by the data range when only one is supplied. Also at: generate_aft_dgm_main.R:49-50 (roxygen: cens_params 'Ignored' when select_censoring = TRUE); generate_aft_dgm_helpers.R:132-141 (select_censoring = FALSE branch defaults each bound separately with a message); ... (full list in the CSV) | no | Decide one rule for uniform censoring under select_censoring = TRUE (ignore cens_params as documented, or default each bound independently with a message as the FALSE branch does) and make code and roxygen state it. | changes behaviour |
| A-218 | 2 | A | generate_aft_dgm_helpers.R | 459 | prepare_working_dataset | S3 | off | traced | The treatment and event columns are used without any coding check: a treatment not coded 0/1 is combined with potential outcomes fixed at 0 and 1, and any event code other than 1 is silently treated as censored. Also at: generate_aft_dgm_helpers.R:460 (event recode); generate_aft_dgm_helpers.R:415-445 (validate_inputs checks existence only); generate_glm_dgm.R:214-220 (treatment documented '0/1 coded' at :42-43, not checked) | no | Add input checks that the treatment column takes only the values 0 and 1 and the event column only 0/1 (or logical). | add-only |
| A-219 | 2 | A | generate_aft_dgm_helpers.R | 551 | process_factor_vars | S3 | off | traced | Dummy columns for factors with three or more levels are named from the raw level labels and then pasted unquoted into model formulas, so labels that are not syntactic names make the DGM build stop with a parser or 'object not found' error. Also at: generate_aft_dgm_helpers_outcome.R:46 and :129 (formula pasted and parsed with as.formula); generate_aft_dgm_helpers.R:222 and :310 (reformulate on zcens_ names, same exposure expected) | no | Sanitise the generated dummy names or back-quote them when the formulas are built. | changes behaviour |
| A-220 | 2 | A | generate_aft_dgm_helpers.R | 934 | calculate_linear_predictors | S3 | off | traced | A spline_spec$var that is not a z_-prefixed covariate column passes validation but makes the build stop with the uninformative error 'Column mismatch in X_super'. Also at: generate_aft_dgm_helpers_outcome.R:208 (validate_spline_spec accepts any df_work column); generate_aft_dgm_main.R:88-90 (roxygen: 'var (variable name)') | no | Check in validate_spline_spec() that var is one of the z_ covariate columns (or map a raw name to it) and state the requirement in the roxygen. | add-only |
| A-221 | 2 | A | generate_aft_dgm_helpers.R | 1437 | compare_multiple_survreg | S3 | off | traced | The tryCatch error handler in the convergence check assigns to its own local copy of conv_status, so a model whose vcov() throws is still reported as converged with no issue recorded. | no | Have the handler return a flag that the enclosing function applies to conv_status. | changes behaviour |
| C-036 | 2 | C | generate_aft_dgm_helpers_outcome.R | 309 | plot_spline_treatment_effect | S3 | off | traced | Two \dontrun examples in the DGM files call functions with arguments that are not formals, and three roxygen sentences describe behaviour the code does not have. Also at: generate_aft_dgm_helpers.R:1769; generate_aft_dgm_helpers.R:1759; generate_aft_dgm_helpers.R:1707; simulate_from_dgm.R:719; generate_glm_dgm.R:51 | no | Rewrite the two examples to the current signatures and correct the three stale sentences and the class name in the print-method title. | add-only |
| C-037 | 2 | C | generate_aft_dgm_main.R | 54 | generate_aft_dgm_flex | S3 | off | traced | Several @param statements of generate_aft_dgm_flex contradict the executing code: the censoring modes, the stored censoring linear predictor, the factor reference level, inheritance of the *_cens vectors and the required spline fields. Also at: generate_aft_dgm_main.R:14; generate_aft_dgm_main.R:35; generate_aft_dgm_main.R:41; generate_aft_dgm_main.R:49; generate_aft_dgm_main.R:82; generate_aft_dgm_main.R:88; generate_aft_dgm_main.R:172 | no | Rewrite the affected @param entries to describe the three censoring modes, the zero covariate contribution, the two-level reference rule, the actual inheritance rule and the spline defaults, and replace the author placeholder. | add-only |
| A-222 | 2 | A | glm_spline_fit.R | 179 | glm_cs_fit | S3 | off | traced | A family object passed as family is reduced to its name and rebuilt with the canonical link, so a non-canonical link such as binomial(link = "log") is discarded without a message. Also at: glm_spline_fit.R:184-188, :36-42 (roxygen) | no | Reject a family object whose link is not the canonical one, or use the supplied object as given and label the profile accordingly. | add-only |
| F-005 | 2 | F | globals.R | 13 | (top-level utils::globalVariables call) | S3 | off | traced | A large share of the 238 names declared in globals.R are no longer used as free (non-standard-evaluation) symbols anywhere in R/, so the declarations are stale and would hide a genuine undefined-global NOTE for those names, and the comment that all declarations are consolidated in globals.R is no longer true. Also at: globals.R:13-274 (entries); forestsearch-package.R:19 (comment); ... (full list in the CSV) | no | Remove the stale names from the globalVariables() vector, confirming each removal with R CMD check, and either consolidate the six per-file declarations or correct the comment in forestsearch-package.R. | add-only |
| B-008 | 2 | B | grf_args.R | 231 | %\|\|% | S3 | off | traced | The null-coalescing operator %\|\|% is defined at top level nine times in R/ (plus two guarded top-level definitions and one function-local one) with two different meanings, and because they share one namespace only the last collated definition is in effect, so the variants that also treat a zero-length value as missing are overridden. Also at: bootstrap_dofuture_main.R:981; calibrate_glm_interaction.R:222; explain_pareto_selection.R:177; fpr_calibration.R:594; frontier_cis.R:362; oc_analyses.R:27; pareto_frontier_table.R:4; plot_pareto_frontier.R:238; ... (full list in the CSV) | no | Keep a single definition of %\|\|% in one utilities file with the intended meaning and delete the others. | moves code |
| C-038 | 2 | C | grf_main.R | 454 | grf.subg.eval | S3 | off | traced | grf.subg.eval() is exported and documented as returning a data frame of evaluation metrics, but its body is a placeholder that ignores all seven arguments and returns NULL. Also at: grf_main.R:433-453 (roxygen); NAMESPACE:223; _pkgdown.yml:190; man/grf.subg.eval.Rd | no | Remove the stub from the exports and the reference index (or restore an implementation) and drop its man page. | changes behaviour |
| A-223 | 2 | A | grf_subg_harm_glm.R | 971 | .cr_get_vcov | S3 | off | traced | create_glm_row() falls back without any message from the sandwich variance to the model-based Poisson variance for effect_measure = "log_RR" when the suggested package sandwich is not installed, so its SE and CI depend on what is installed. Also at: grf_subg_harm_glm.R:823-825, 830-833, 843-855 | no | Warn (or stop) when the robust variance is called for by the effect measure but sandwich is unavailable, and let the missing-MASS condition surface. | add-only |
| A-224 | 2 | A | guohe_adaptive_r.R | 251 | guohe_adaptive_r | S3 | off | traced | Every error from the training-fold guohe_algorithm3() call is swallowed by try(silent = TRUE) and the fold skipped, so a configuration error surfaces only as the generic 'Cross-validation objective undefined for every r.' and a partial failure is not reported at all. Also at: guohe_adaptive_r.R:57-95, 197-202, 260, 267-270 | no | Validate the column arguments up front and keep the first engine error message to report when all folds fail (and warn when some do). | add-only |
| A-225 | 2 | A | leff_calibration_io.R | 107 | save_leff_calibration | S3 | off | traced | n_min is validated only as a positive numeric under the label 'must be a positive integer' and is then formatted with %d after the file has been written, so a non-integer n_min saves the bundle and then errors without returning the path. Also at: leff_calibration_io.R:62-63, 158-159; calibrate_null_correction.R:119-120, 126-127 | no | Check integer-valuedness in the validation (or format with %s / %g) so the message cannot fail after the write. | add-only |
| C-039 | 2 | C | leff_calibration_io.R | 255 | latest_rds | S3 | off | traced | latest_rds() orders by an 8-digit substring of the filename whenever every match has one, although its roxygen says it uses modification time, and it does not break ties between files with the same stamp. Also at: leff_calibration_io.R:200-204, 212, 233, 248-253 | no | Make the roxygen describe the date-stamp preference (or make the code use modification time) and break stamp ties by modification time. | changes behaviour |
| C-040 | 2 | C | mr_estimates_table.R | 318 | .mr_table_from_fs | S3 | off | traced | The template-free table always heads the de-biased column 'De-biased (95% IJ CI)' and footnotes it as an infinitesimal-jackknife interval, although under ci_method = "wald" the limits it prints are the robust-SE interval. Also at: mr_estimates_table.R:326-328 (footnote); values read at mr_estimates_table.R:253-254, 290-292; template route at mr_estimates_table.R:143-147 | no | Make the header and footnote of the template-free route depend on ci_method, as the template route's footnote already does. | add-only |
| A-226 | 2 | A | mrct_four_sets.R | 129 | .mrct_check_trial | S3 | off | traced | .mrct_check_trial() enforces 0/1 coding on the treatment, region and truth columns but not on the event column, so an event column in another coding that Surv() accepts (for example 1 = censored, 2 = event) passes and yields wrong events_exp / events_ctl counts beside a correct Cox fit. Also at: mrct_simulation.R:825-826 (event sums); mrct_trial_helpers.R:95-101 (hand-off), :147-148 and :275-276 (counts printed and drawn); callers mrct_four_sets.R:49, mrct_trial_analysis.R:92 | no | Add event.name to the columns that must be 0/1 (or logical) in .mrct_check_trial(). | add-only |
| C-041 | 2 | C | mrct_simulation.R | 31 | mrct_region_sims | S3 | off | traced | n_sample is documented as 'If NULL (default), uses the entire super-population from dgm', but the code sets the trial size to nrow(dgm$df_super) and draws that many rows with replacement, then applies entry times, re-randomisation and censoring as for any other size. Also at: mrct_simulation.R:309-311, :421; simulate_from_dgm.R:221, :251 | no | Reword the @param to say that NULL sets the trial size to nrow(dgm$df_super), sampled with replacement like any other size. | add-only |
| C-042 | 2 | C | mrct_simulation.R | 53 | mrct_region_sims | S3 | off | traced | conf_force is documented with the default c("z_age <= 65", "z_bm <= 0", ...) while its formal default is NULL, which the body silently replaces by those five case-study cuts, so 'no forced cuts' requires character(0) (undocumented) whereas in forestsearch() conf_force = NULL means none. Also at: mrct_simulation.R:271 (formal), :326-328 (NULL replaced), :529 (passed on); mrct_case_study.R:292; forestsearch_main.R:1655 (conf_force = NULL means none) | no | Either put the vector in the signature or document that NULL selects the case-study cuts and that character(0) requests none. | add-only |
| C-043 | 2 | C | mrct_simulation.R | 110 | mrct_region_sims | S3 | off | traced | The @return items call n_sg, hr_sg, POhr_sg, prev_sg and the _train versions the 'Subgroup' quantities directly after sg_found ('identified subgroup'), but the code computes them on treat.recommend == 1, the complement of the identified subgroup (the whole set when nothing is found), and nothing in the documentation says so. Also at: mrct_simulation.R:111-116, :139-140 (details steps 5-6), :1348 (label_sg), :1893-1896 (table labels); code at :592-597, :605-634; forestsearch_helpers.R:69, :137 | no | State in the @return items (and the table and plot labels) that these columns describe the patients outside the identified subgroup sg_found, equal to the whole set when none is declared. | add-only |
| F-006 | 2 | F | mrct_simulation.R | 392 | mrct_region_sims | S3 | off | traced | mrct_region_sims() calls progressr::handlers("progress") unconditionally, discarding the previous handlers it returns and never restoring them, and the 'progress' handler depends on the 'progress' package, which forestsearch does not declare or guard. Also at: mrct_simulation.R:692 (with_progress); DESCRIPTION Imports/Suggests (no 'progress'); mrct_simulation.R:392 | no | Leave the handler choice to the caller, or save the previous handlers and restore them with on.exit and guard or declare the 'progress' dependency. | changes behaviour |
| F-007 | 2 | F | mrct_simulation.R | 1029 | validate_mrct_data | S3 | off | traced | The base package 'tools' is used through tools:: at four places but is not listed in DESCRIPTION Imports (stats, utils and parallel are). Also at: pareto_frontier_table.R:186; oc_analyses.R:1012; render_forestplot.R:332; globals.R:290; oc_analyses.R:1012 | no | Add tools to Imports in DESCRIPTION. | add-only |
| A-227 | 2 | A | mrct_simulation.R | 1197 | create_dgm_for_mrct | S3 | off | traced | confounder_var is not validated: given without confounder_effect it is silently ignored, and a name that is not one of the model's single z_ columns (a typo, a variable outside the fixed covariate set, or a factor with more than two levels) is appended as a new coefficient and ends in the unrelated error 'Column mismatch in X_super'. Also at: mrct_simulation.R:1199-1208, :1090-1093 and :1113-1117 (roxygen); generate_aft_dgm_helpers_outcome.R:154-156; generate_aft_dgm_helpers.R:1032-1035, :545-550 | no | Check that confounder_effect accompanies confounder_var and that z_<confounder_var> is one of the model's coefficient names, with an error that names the argument. | add-only |
| D-014 | 2 | D | mrct_simulation.R | 1642 | summaryout_mrct | S3 | off | traced | The optional pop_summary (default NULL) is the first formal and the required mrct_sims the second, so a positional call summaryout_mrct(results) fails, and the element names pop_summary must carry are neither documented nor produced by any function in the package. Also at: mrct_simulation.R:1562-1563 (@param pop_summary), 1729-1737 (process_pop_summary reads ps$AHR, ps$ITT_unadj, ps$ITT_sR, ps$ITT_sRw, ps$W_1) | no | Document the pop_summary elements and either move mrct_sims to the first position or state that the arguments must be named. | changes behaviour |
| A-228 | 2 | A | mrct_simulation.R | 1787 | summaryout_mrct | S3 | off | traced | The 'Subgroup Classification' rows are counted by unanchored, case-insensitive substring matches on the sg_found label, so a subgroup is counted in a row whenever any covariate name merely contains the fragment (e.g. 'age' in 'stage', 'bm' in 'bmi', 'reg' in 'regimen'). Also at: mrct_simulation.R:1792, 1800, 1803, 1806, 1810, 1814, 1818, 1822 | no | Classify from the variable names parsed out of each cut label against a user-supplied (or documented) name-to-category map, or document that the rows are name-pattern matches tied to the original case-study variable names. | changes behaviour |
| A-229 | 2 | A | oc_analyses.R | 204 | create_subgroup_indicator | S3 | off | traced | A factor definition that does not parse as '<column>.<level>' with an existing column is skipped without warning, so the internal helper returns membership 1 for every row when no definition resolves. Also at: oc_analyses.R:199-200 (split on '.', definitions without a dot skipped); run_simulation_analysis.R:1048-1055 (only caller, a fallback branch) | no | Make the helper fail or warn when a definition cannot be resolved to a column and level, or remove the fallback branch if it cannot be reached. | changes behaviour |
| D-015 | 2 | D | oc_analyses.R | 457 | summarize_single_analysis | S3 | off | traced | summarize_simulation_results() works only on a data.table because its helper indexes with `with = FALSE`, whereas the sibling table functions coerce data.frame input first. Also at: oc_analyses.R:430 (subset() keeps the input class), :404 (@param says data.table); siblings that coerce: oc_analyses.R:663-665 (format_oc_results), simulation_tables.R:290 (build_classification_table) | no | Coerce the input to data.table at the top of summarize_simulation_results(), as the sibling functions do. | add-only |
| A-230 | 2 | A | pareto_frontier_table.R | 425 | .render_pareto_gt | S3 | off | traced | With include_dominated = TRUE the gt header still reads 'Pareto frontier on (...)' with the subtitle '<n> non-dominated candidate(s)', where n is the number of all passing candidates, dominated ones included. Also at: pareto_frontier_table.R:432 (title fixed as 'Pareto frontier on (...)') | no | Pass the mode to the renderer and word the title and subtitle for the passing-set view (number of passing candidates, of which k on the frontier). | changes behaviour |
| D-016 | 2 | D | plot_km_band_forestsearch.R | 182 | plot_km_band_forestsearch | S3 | off | traced | The five subgroup-plot entry points use different defaults and spellings for the same arguments (outcome/event/treatment column defaults, risk-table interval, CI and log-rank switches, accepted fs.est types), and none takes the column names from the fit. Also at: plot_sg_results.R:134-149, :159; plot_sg_weighted_km.R:97-105, :271-276; plot_sg_glm_outcomes.R:95-96, :126; plot_subgroup_results_forestplot.R:199-204, :241 | no | Document one naming and default convention for the plot family and, when names are not supplied, resolve the outcome, event and treatment columns from the fit's stored call arguments. | changes behaviour |
| C-044 | 2 | C | plot_pareto_frontier.R | 39 | plot_pareto_frontier | S3 | off | traced | @details says that when the selected subgroup is not on the frontier 'the selected marker appears off the polyline', but only frontier rows are plotted, so in that case no selected marker is drawn at all. Also at: plot_pareto_frontier.R:98, 107-108, 157-165 | no | Either add the selected row to the plotted points when it is off the frontier, or correct the @details paragraph. | add-only |
| D-017 | 2 | D | plot_pareto_frontier.R | 223 | plot_pareto_frontier | S3 | off | traced | The two frontier plots decide whether to shade the effect-neighborhood band from sg_focus alone, while pareto_frontier_table() and explain_pareto_selection() also require selection_rule to be "neighborhood" or "both", so for a selection_rule = "pareto" fit the plots shade a band that the table and explanation treat as not used. Also at: plot_pareto_combined.R:282-283; plot_pareto_combined.R:375-378; pareto_frontier_table.R:204-205; explain_pareto_selection.R:115-116 | no | Decide band applicability in one shared helper (sg_focus and selection_rule) used by the plots, the table and the explanation, and align the sibling defaults. | changes behaviour |
| C-045 | 2 | C | plot_sg_distribution.R | 79 | plot_sg_distribution | S3 | off | traced | The first roxygen example passes fs$grp.consistency$result_new, an element that does not exist on a forestsearch result, so the example as written stops with "Column 'any_found' not found in 'results'". | no | Replace the first example with one built from mrct_region_sims() output. | add-only |
| A-231 | 2 | A | plot_sg_glm_outcomes.R | 385 | plot_sg_glm_outcomes | S3 | off | traced | Count outcomes share the binary bar-chart branch: the y-axis is capped at 1 whenever no effect annotation is added, and the axis is labelled 'Event Rate (per person-year)' even when no offset is supplied and arm means are what is plotted. Also at: plot_sg_glm_outcomes.R:198, 343-347 | no | Apply the upper cap of 1 to binary outcomes only and derive the count axis label from whether an offset is in use. | changes behaviour |
| A-232 | 2 | A | plot_sg_glm_outcomes.R | 496 | plot_sg_glm_outcomes | S3 | off | traced | The generated figure note always says '95 percent CI' or '+/- 1.96 SE' although the error bars and effect intervals are drawn at the user's conf.level. Also at: plot_sg_glm_outcomes.R:499 | no | Build the note text from conf.level and z_alpha. | add-only |
| A-233 | 2 | A | plot_sg_results.R | 138 | plot_sg_results | S3 | off | traced | The documented arguments by.risk and show_events have no effect: by.risk is auto-computed and passed to plot_sg_km(), which never reads it, and show_events is never referenced. Also at: plot_sg_results.R:146, 208-211, 284, 547, 451-464 | no | Implement the risk-table interval and event-count display, or remove the two arguments and the unused median computation from the signature and documentation. | changes behaviour |
| A-234 | 2 | A | plot_sg_results.R | 278 | plot_sg_results | S3 | off | traced | The returned plots element is always an empty list because each plotting helper returns invisible(NULL) and assigning NULL with $<- adds no element, so print() lists no available plots and plot.fs_sg_plot() always reports 'No plots available', contrary to the documented return value. Also at: plot_sg_results.R:300, 311, 661, 775, 849, 894, 920-935 | no | Record which plots were drawn (a marker or a recorded plot per type), or correct the documented return value and the plot method. | changes behaviour |
| C-046 | 2 | C | plot_sg_results.R | 869 | print.fs_sg_plot | S3 | off | traced | The \dontrun examples of print.fs_sg_plot() and plot.fs_sg_plot() build their object with plot_subgroup(fs), an unexported helper with a six-argument signature that never returns an fs_sg_plot object. Also at: plot_sg_results.R:915 | no | Change both examples to build the object with plot_sg_results(). | add-only |
| C-047 | 2 | C | plot_sg_weighted_km.R | 29 | plot_sg_weighted_km | S3 | off | traced | The roxygen for sg0_name/sg1_name promises automatic titles "H: {definition}" and "Hc: NOT {definition}", but the code uses the bare definition and "NOT <definition>", and rounds thresholds only in the H title. Also at: plot_sg_weighted_km.R:244, 255, 260; plot_km_band_forestsearch.R:117-118 | no | Make the roxygen describe the titles actually generated and build the Hc title from the tidied label. | add-only |
| A-235 | 2 | A | plot_sg_weighted_km.R | 303 | plot_sg_weighted_km | S3 | off | traced | The H/Hc split uses `df[df$treat.recommend == 0, ]` with no NA handling, so rows with a missing treat.recommend become all-NA rows in both subgroups and are counted in the n shown in the panel titles, unlike plot_sg_results() which excludes them. Also at: plot_sg_weighted_km.R:304; plot_sg_glm_outcomes.R:181-182 | no | Drop or reject rows whose treat.recommend is NA before splitting, as plot_sg_results() does. | add-only |
| A-236 | 2 | A | plot_subgroup_results_forestplot.R | 645 | plot_subgroup_results_forestplot | S3 | off | traced | subgroup_list elements are not validated, and an element with no `type` makes Filter() pick the wrong elements because `x$type == "reference"` is logical(0) for it and drops out of the index vector, shifting the later TRUE/FALSE values by one. Also at: plot_subgroup_results_forestplot.R:667 | no | Check at entry that every subgroup_list element has subset_expr, name and a type of "reference" or "posthoc", and stop otherwise. | add-only |
| A-237 | 2 | A | plot_subgroup_results_forestplot.R | 747 | plot_subgroup_results_forestplot | S3 | off | traced | est.scale is never validated and its two uses test different values, so any string other than "hr" or "1/hr" (e.g. "HR") leaves treatment un-reversed but swaps the benefit/questionable rows, their bias-corrected estimates and the CV agreement order, with no error. Also at: plot_subgroup_results_forestplot.R:204, 611, 713 | no | Validate est.scale against c("hr", "1/hr") at entry. | add-only |
| A-238 | 2 | A | plot_subgroup_results_forestplot.R | 1150 | plot_subgroup_results_forestplot | S3 | off | traced | For ratio-scale GLM measures (OR/RR/IRR) a finite point estimate outside the axis limits is set to NA in the plot and in the returned data and printed as "NE" with the footnote "NE = not estimable", a behaviour the roxygen does not mention. Also at: plot_subgroup_results_forestplot.R:1151-1156, 1197 | no | Document the rule in the roxygen and label such rows as outside the axis range rather than not estimable, keeping the numeric estimate in the returned data. | changes behaviour |
| F-008 | 2 | F | render_forestplot.R | 209 | print.fs_forest_theme | S3 | off | traced | The print.fs_forest_theme example runs create_forest_theme() outside \dontrun with no requireNamespace guard, and create_forest_theme() stops when the Suggests-only package forestploter is missing although it never calls forestploter. Also at: render_forestplot.R:141 | no | Drop the forestploter requirement from the list constructor, or guard the example with requireNamespace(). | changes behaviour |
| F-009 | 2 | F | render_forestplot.R | 335 | save_forestplot | S3 | off | traced | save_forestplot() opens a graphics device and closes it only on the success path, so an error while drawing leaves the file device open and current. Also at: render_forestplot.R:358 | no | Register on.exit(grDevices::dev.off()) immediately after the device is opened. | add-only |
| C-048 | 2 | C | run_simulation_analysis.R | 341 | run_simulation_analysis | S3 | off | traced | The 'GLM Parameters' roxygen section says outcome_type, effect_measure and offset.name go only to the estimation step and not to forestsearch(), 'which uses Cox PH for subgroup identification in v0.1.x', but the code forwards all three to forestsearch(). | no | Rewrite the section to say these parameters are forwarded to forestsearch() as well as to the estimation step. | add-only |
| B-009 | 2 | B | run_simulation_analysis.R | 1464 | .extract_grf_estimates_gen | S3 | off | traced | Three identically bodied merge closures are defined inside .extract_grf_estimates_gen(), and the non-detection fill block and the GLM effect closure are copied between the FS and GRF extractors. Also at: run_simulation_analysis.R:1491-1498, :1513-1520 (same merge closure body); non-detection fill repeated at :1009-1038, :1056-1073, :1381-1399, :1408-1425; ... (full list in the CSV) | no | Factor the merge, the non-detection fill and the GLM effect closure into single internal helpers shared by both extractors. | moves code |
| C-049 | 2 | C | run_subgroup_sims.R | 200 | subgroup_glm | S3 | off | traced | The effect_measure roxygen and an inline comment of subgroup_glm() state that make_effect_estimator()'s own binary default is "RD", but make_effect_estimator() resolves a NULL binary measure to "OR". Also at: run_subgroup_sims.R:293-296 (inline comment repeating the claim and pointing at the Rd) | no | Remove the stale divergence note from the roxygen and from the inline comment. | add-only |
| C-050 | 2 | C | run_subgroup_sims.R | 597 | run_subgroup_sims | S3 | off | traced | The @param entries for baseline, dgm and n say the fixed baseline is survival-only, that a glm_dgm always resamples and that n is required for a glm_dgm, whereas the code runs baseline = "fixed" for a glm_dgm and then requires n to be NULL. Also at: run_subgroup_sims.R:585-586 (@param dgm: 'baseline = "fixed" (survival only)'); run_subgroup_sims.R:601-603 (@param n: 'Required for "resample" and for a "glm_dgm"'); ... (full list in the CSV) | no | Update the dgm, baseline and n @param entries (and the 'n is required' clause in the details) to describe the GLM fixed-baseline mode. | add-only |
| D-018 | 2 | D | run_subgroup_sims.R | 926 | run_subgroup_sims | S3 | off | traced | For a glm_dgm the returned object stores the literal "glm" in both design and sim_config$baseline, so nothing in the result records whether the run was fixed-baseline or resampled, unlike survival results. Also at: run_subgroup_sims.R:919 (design = "glm"); run_subgroup_sims.R:862-866 (verbose message prints the real baseline); compare_subgroup_sims.R:117-124 (expect_designs guard compares $design) | no | Record the actual baseline in sim_config$baseline for GLM results, keeping design = "glm" if the documented payload format requires it. | changes behaviour |
| C-051 | 2 | C | sim_aft_gbsg.R | 216 | create_gbsg_dgm | S3 | off | traced | The deprecation message of create_gbsg_dgm() recommends a generate_aft_dgm_flex() call with argument names that function does not have, and then says to 'wrap with setup_gbsg_dgm()', which takes no DGM argument. Also at: sim_aft_gbsg.R:223 ('Then wrap with setup_gbsg_dgm()'); sim_aft_gbsg.R:805-807 (simulate_from_gbsg_dgm message) | no | Rewrite the deprecation messages to point at setup_gbsg_dgm() with its real arguments, and say that simulate_from_dgm() needs a setup_gbsg_dgm() object. | add-only |
| C-052 | 2 | C | sim_aft_gbsg.R | 758 | simulate_from_gbsg_dgm | S3 | off | traced | The roxygen says a positive muC_adj increases censoring although the code lengthens censoring times (less censoring), and it calls min_cens / max_cens required for uniform censoring although the code silently defaults them. Also at: sim_aft_gbsg.R:760-763 (min_cens / max_cens 'Required') vs sim_aft_gbsg.R:919-922 (silent defaults) | no | Correct the muC_adj description and document the uniform-bound defaults (or make the bounds genuinely required). | add-only |
| C-053 | 2 | C | simulate_from_dgm.R | 65 | simulate_from_dgm | S3 | off | traced | The roxygen calls hrz_crit a log-HR threshold and hrz_flag a harm-zone indicator, but the flag is computed, exactly as the same @param states, from lin_pred_1 - lin_pred_0, which is on the AFT log-time scale with the opposite sign to the log hazard ratio (and a factor tau), so with hrz_crit = 0 it marks subjects who benefit. Also at: simulate_from_dgm.R:108 (return item 'Individual harm-zone indicator'); simulate_from_dgm.R:360-362 (the computation) | no | Make the description and the computation agree: either describe hrz_crit as a threshold on the log-time-ratio difference and drop the 'log-HR'/'harm-zone' wording, or compute the flag from the stored individual log hazard ratio. | changes behaviour |
| A-239 | 2 | A | simulate_from_dgm.R | 330 | simulate_from_dgm | S3 | off | traced | cens_adjust is silently ignored when the DGM's censoring type is "uniform", and calibrate_cens_adjust() on such a DGM then fails with a message that blames the target. Also at: simulate_from_dgm.R:851-867 (calibrate_cens_adjust objective); simulate_from_dgm.R:55-57 and :135-139 (documentation with no exception) | no | Reject (or document as ignored) a non-zero cens_adjust for uniform censoring and stop early with a specific message in calibrate_cens_adjust() for uniform DGMs. | add-only |
| A-240 | 2 | A | simulate_from_dgm.R | 483 | check_censoring_dgm | S3 | off | traced | The by-arm reference censoring rates and KM medians pair the observed y / event of df_super with df_super$treat, which is an arbitrary first-half / second-half relabelling when the DGM was built with draw_treatment = TRUE. Also at: simulate_from_dgm.R:538-539 (ref_km_0 / ref_km_1) | no | Take the reference arm from a column that preserves the observed assignment, or drop the by-arm reference rows when the DGM re-drew treatment. | changes behaviour |
| D-019 | 2 | D | simulate_from_glm_dgm.R | 106 | simulate_from_glm_dgm | S3 | off | traced | The AFT and GLM sibling simulators and builders differ in what n = NULL does with draw_treatment, in argument defaults and validation style, and in the fields of the returned DGM object. Also at: simulate_from_dgm.R:221-234 (n = NULL keeps super-population treatment, draw_treatment not applied); simulate_from_glm_dgm.R:121-123 (Bernoulli allocation); ... (full list in the CSV) | no | Document the differences side by side, or align the two families on one n = NULL rule and one set of defaults and return fields. | changes behaviour |
| B-010 | 2 | B | subgroup_consistency_helpers.R | 1155 | remove_redundant_subgroups | S3 | off | traced | remove_redundant_subgroups() has no caller and is not exported, and its filter does not do what its documentation says (it effectively de-duplicates on HR alone). Also at: _pkgdown.yml:112 and man/remove_redundant_subgroups.Rd (reference entries only) | no | Remove the function together with its man page and pkgdown reference entry, or, if it is to be kept, compute a row-wise any-difference before indexing. | moves code |
| A-241 | 2 | A | summarize_bootstrap_results.R | 92 | summarize_bootstrap_results | S3 | off | traced | The success rate and the 'Failed to find subgroup' count are computed from non-missing H_biasadj_2 and not from any_found, so resamples that identified a subgroup but had a non-finite fit are reported as not finding one, and the figure can disagree with n_found in the same returned object. Also at: summarize_bootstrap_results.R:275-277, 337-340; bootstrap_summaries_helpers.R:267-271; bootstrap_dofuture_main.R:675-678; contrast summarize_bootstrap_subgroups.R:107-117 | no | Report identification from any_found and report the number of usable bias-corrected estimates as a separate figure. | changes behaviour |
| A-242 | 2 | A | summarize_cv_results.R | 707 | .cv_decompose_no_sg | S3 | off | traced | When GRF was not used at all (use_grf = FALSE, the forestsearch() default), the no-subgroup decomposition assigns every no-subgroup fold to 'GRF returned no cut' and none to 'consistency rejected', and the table footnote then advises relaxing dmin.grf. Also at: summarize_cv_results.R:245, 254, 269-291, 528-532, 1136-1141; forestsearch_cross_validation.R:1099-1105, 1136 | no | Detect that GRF was not used (all grf_cuts missing, or a flag carried from the fit) and suppress or relabel the GRF-specific tables and guidance in that case. | changes behaviour |
| C-054 | 2 | C | summary_subgroup_sims.R | 83 | summary.subgroup_sims | S3 | off | traced | The @return section lists the fields and says legacy calls return 'exactly this set', but the returned list also contains n_per_trial, single_cats and combo_cats. Also at: summary_subgroup_sims.R:293-311 (returned list) | no | Add the three fields to @return and drop or correct the word 'exactly'. | add-only |
| A-243 | 2 | A | summary_utility_functions.R | 763 | sg_tables | S3 | off | traced | Four documented formals are accepted and never read: sg_tables(potentialOutcome.name=), SG_tab_estimates(draws=, details=) and format_results(return_medians=). Also at: summary_utility_functions.R:606 (draws), summary_utility_functions.R:607 (details), summary_utility_functions.R:128 (return_medians) | no | Either pass sg_tables' potentialOutcome.name into the override list of the SG_tab_estimates call, or remove/deprecate the inert formals and their roxygen lines. | changes behaviour |
| A-244 | 2 | A | summary_utility_functions.R | 934 | sg_tables | S3 | off | traced | which_df is not validated, so any value other than "est" or "testing" leaves the local df unassigned and the is.null(df) guard then tests whatever object named df is found outside the function (normally the stats::df function), ending in an unrelated error. Also at: summary_utility_functions.R:935, summary_utility_functions.R:956 | no | Validate which_df against the two allowed values (match.arg) so an unknown value stops with a clear message. | add-only |
| B-011 | 2 | B | synthetic_data_perturbation.R | 481 | generate_bootstrap_with_noise | S3 | off | traced | generate_bootstrap_synthetic() and generate_bootstrap_with_noise() are two exported near-duplicate implementations of one bootstrap-plus-perturbation routine, with different argument order, defaults and categorical handling. Also at: synthetic_data_perturbation.R:53-231 | no | Keep one implementation and make the other a thin wrapper with one argument convention. | moves code |
| F-010 | 2 | F | test_hte_crump.R | 51 | sandwich_var_ols | S3 | off | traced | MASS::ginv() is called at six sites with no requireNamespace() guard although MASS is listed only under Suggests. Also at: test_hte_crump.R:66, test_hte_crump.R:107, test_hte_crump.R:111, test_hte_crump.R:176, test_hte_crump.R:222 | no | Guard the generalized-inverse fallback with requireNamespace("MASS", quietly = TRUE) and a clear error, or move MASS to Imports. | add-only |
| A-245 | 2 | A | compare_selection_rules.R | 737 | compare_selection_rules | S3 | off | suspected | When a single signature group covers all valid fits but plot_pareto_combined() returns NULL, `plot_combined_subsets[[1L]]` indexes an empty list and aborts the call after all fits are done; no real input reaching that state was found. Also at: compare_selection_rules.R:716-731 (the NULL group_plot state handled two lines earlier) | no | Set plot_combined only when the single covering group actually produced a plot, otherwise leave it NULL and report combined_skip_reason. | add-only |
| A-246 | 2 | A | forestsearch_helpers.R | 802 | reset_workers | S3 | off | suspected | reset_workers() labels the plan with class(future::plan())[1L], which under the installed future is expected to be the generic "FutureStrategy" rather than the strategy name, so the returned previous/current values and the messages would not identify the plan. Also at: forestsearch_helpers.R:838; mrct_simulation.R:383 | no | Derive the plan label from the first class that is not FutureStrategy/tweaked (or test inherits() against the known strategy names). | add-only |
| A-247 | 2 | A | generate_aft_dgm_helpers.R | 1445 | compare_multiple_survreg | S3 | off | suspected | The 'max iterations reached' convergence check tests fit$maxiter, a component survreg objects do not appear to carry, so the check can never fire. | no | Compare fit$iter with the iteration limit actually used for the fit instead of a component the fit does not have. | changes behaviour |
| A-248 | 2 | A | mrct_simulation.R | 407 | mrct_region_sims | S3 | off | suspected | The replicate loop combines .errorhandling = "pass" with .combine = "rbind" and never separates condition objects from result rows, so an error raised outside the inner tryCatch() (for example an invalid sim_args name) is not reported as itself. Also at: mrct_simulation.R:408, :428, :706-707 | no | Collect the replicates as a list, separate the condition objects, report their count and first message (stop when all fail), and rbind the successful rows. | changes behaviour |
| A-249 | 2 | A | run_simulation_analysis.R | 476 | run_simulation_analysis | S3 | off | suspected | Replicate data are generated by set.seed(seed_base + sim_id) under whatever RNG kind is active, so the same sim_id is expected to give different data in a top-level call than inside a %dofuture% loop run with seed = TRUE, the idiom the simulation article shows. Also at: run_simulation_analysis.R:482, 509; simulate_from_dgm.R:208; simulate_from_glm_dgm.R:97 | unknown | Fix the RNG kind where the replicate seed is set, or document that the loop must use seed = NULL as run_subgroup_sims() does. | changes behaviour |
| A-250 | 1 | A | dina_df.R | 103 | dina | S4 | on | traced | dina() stores match.call() as fit$call, and forestsearch invokes it through do.call() with the evaluated data frame, so the stored call embeds the whole data set and print(summary(fs$dina_res)) deparses it to the console. Also at: forestsearch_helpers.R:1531-1534; forestsearch_main.R:3200-3206 (do.call(dina, list(df = <data>, ...))); forestsearch_main.R:2802, 3959 (fit stored as dina_res); ... (full list in the CSV) | no | Have the two forestsearch call sites build the call with symbols, or replace the stored call's df argument with a placeholder. | changes behaviour |
| A-251 | 1 | A | forestsearch_helpers.R | 90 | get_dfpred | S4 | on | traced | get_dfpred()'s version argument is accepted and never used, and the roxygen does not say it has no effect. Also at: forestsearch_helpers.R:65-66 | no | State in the roxygen that version is ignored (or deprecate the argument). | add-only |
| D-020 | 1 | D | forestsearch_helpers.R | 1144 | .resolve_dina_args | S4 | on | traced | The internal DINA/GRF selector helpers carry fallback defaults opposite to forestsearch()'s public defaults, and the two sibling selectors disagree on the adverse_outcome default. Also at: forestsearch_helpers.R:1487, 1838, 1843, 1848; grf_args.R (grf_selection = "tree" formals) | no | Make the internal fallbacks match the public defaults, or drop them so callers must supply the values. | add-only |
| C-055 | 1 | C | forestsearch_helpers.R | 1981 | .forestsearch_grf_select | S4 | on | traced | Internal documentation of how GRF membership and the disjunctive sg.harm string are handled contradicts the code, and get_dfpred()'s roxygen does not mention the disjunction input form it evaluates. Also at: forestsearch_helpers.R:1823-1824; :43-70 (get_dfpred roxygen); :2254 | no | Update the stale comment and internal roxygen to the structured-definition route and document the disjunction input form in get_dfpred(). | add-only |
| F-011 | 1 | F | forestsearch_main.R | 185 | .default_parallel_workers | S4 | on | traced | The default-worker cap is applied only when _R_CHECK_LIMIT_CORES_ parses as logical TRUE, so values such as "warn" leave the default at 80% of the physical cores, while the sibling helper .fs_grf_num_threads() treats any value other than "", "false", "FALSE" as an active limit. Also at: utils_check_threads.R:24-25; forestsearch_main.R:1631 | no | Parse the environment variable once, in a shared helper, treating any non-empty value other than false as an active limit. | changes behaviour |
| B-012 | 1 | B | forestsearch_main.R | 2062 | forestsearch | S4 | on | traced | forestsearch() carries a set of statements and guards that cannot have any effect (unused worker clamp, repeated match.arg/default resolution, always-true exists() tests, try-error tests on tryCatch results, NULL assignment to an absent list element). Also at: forestsearch_main.R:2063, 2073-2077 (n_workers/max_cores never read); 2172-2174 (repeat of match.arg at :1730-1732); 2186 (exists("effect_measure") on a formal); 2198-2200 (repeat of :1797-1799); ... (full list in the CSV) | no | Delete the no-op statements and either use or remove the worker clamp. | moves code |
| B-013 | 1 | B | fs_mr_inference_methods.R | 96 | .fs_mr_reselection_from_focus | S4 | on | traced | The switch branches for the alias spellings maxcons, effMaxSG and effMinSG, and the error handler around the normalizer, cannot be reached for a character-scalar sg_focus because the focus is normalized to its canonical spelling first. Also at: fs_mr_inference_methods.R:104; fs_mr_inference_methods.R:105; fs_mr_inference_methods.R:89 (tryCatch handler) | no | Delete the three alias branches and the tryCatch, leaving the canonical names only. | moves code |
| D-021 | 1 | D | get_fsdata.R | 89 | get_FSdata | S4 | on | traced | replace_med_grf defaults to TRUE in get_FSdata() but FALSE in forestsearch(), and the median cut is rounded to two decimals on one path and one decimal on the other. Also at: forestsearch_main.R:1659; get_fsdata.R:480; get_FSdata_helpers.R:430; get_fsdata.R:324-325 | no | Give replace_med_grf the same default in both functions and use one rounding convention for the median cut. | changes behaviour |
| A-252 | 1 | A | get_fsdata.R | 429 | get_FSdata | S4 | on | traced | x[1:toprint] with toprint = 0 prints 'NA' in the details output when the vector is character(0). Also at: get_fsdata.R:419 | no | Index with seq_len(toprint) or use head(). | add-only |
| D-022 | 1 | D | subgroup_consistency_helpers.R | 306 | run_single_consistency_split | S4 | on | traced | The per-split success test is non-strict (>=) on the GLM branch and strict (>) on the survival branch of the same function. Also at: subgroup_consistency_helpers.R:320, :1758, :1769 | no | Use one comparison convention for both outcome families and state it in the threshold documentation. | changes behaviour |
| D-023 | 1 | D | subgroup_consistency_helpers.R | 1385 | evaluate_subgroup_consistency | S4 | on | traced | The two exported evaluators default consistency_method to "split" while subgroup.consistency() and forestsearch() default to "resample". Also at: subgroup_consistency_helpers.R:1675 | no | Align the evaluators' default with the engine's, or remove the default so the method must be stated. | changes behaviour |
| C-056 | 1 | C | subgroup_consistency_helpers.R | 1958 | evaluate_consistency_twostage | S4 | on | traced | Stage 2 always runs whole batches, so the total number of splits can exceed n.splits.max (documented as 'Maximum total splits') by up to batch.size - 1, and the variable all_flags is assigned and never read. Also at: subgroup_consistency_helpers.R:1613 (roxygen), :1954 (unused all_flags); forestsearch_main.R:847-848 | no | State the bound as approximate (rounded up to a whole batch) in the roxygen, or cap the last batch at the remaining budget, and drop the unused variable. | add-only |
| B-014 | 1 | B | subgroup_consistency_main.R | 429 | subgroup.consistency | S4 | on | traced | Several guards cannot fire inside the installed package: exists() checks for functions defined in the same namespace, a zero-row check repeated after an earlier return, an emptiness check after de-duplication, and an always-true any.found > 0 test. Also at: consistency_resample.R:184 (fallback :187-189); consistency_resample.R:495; subgroup_consistency_main.R:537 (nrow half of the condition only); subgroup_consistency_main.R:625; subgroup_consistency_main.R:1123 | no | Remove the unreachable guards and call .fs_adjust_terms() directly. | moves code |
| D-024 | 1 | D | subgroup_consistency_main.R | 653 | subgroup.consistency | S4 | on | traced | The pool-truncation warning names max_subgroups_search, which is forestsearch()'s argument, although the argument of this exported function is stop_Kgroups. Also at: subgroup_consistency_main.R:1174 (comment) | no | Name stop_Kgroups in the warning (optionally adding that forestsearch() calls it max_subgroups_search). | add-only |
| G-010 | 2 | G | bootstrap_analysis_dofuture.R | 300 | bootstrap_results | S4 | on | traced | 737 lines, 12 arguments, brace depth 6. Replicate setup, the inner forestsearch() call, four subgroup fits, bias corrections, MR capture and timing in one foreach body. | no | Split into named stage functions with the current behaviour pinned by tests before the split. | moves code |
| G-002 | 1 | G | forestsearch_main.R | 1622 | forestsearch | S4 | on | traced | 2362 lines, 86 arguments, brace depth 6. Validation, alias and threshold resolution, propensity scores, the GRF and DINA selectors, cut construction, search, consistency, re-selection, MR and return assembly in one function. | no | Split into named stage functions with the current behaviour pinned by tests before the split. | moves code |
| G-009 | 1 | G | fs_mr_inference.R | 584 | fs_mr_inference | S4 | on | traced | 557 lines, 29 arguments, brace depth 4. Argument resolution, family enumeration, several variance constructions, field calibration, admission and output shaping in one function. | no | Split into named stage functions with the current behaviour pinned by tests before the split. | moves code |
| G-006 | 1 | G | subgroup_consistency_helpers.R | 1654 | evaluate_consistency_twostage | S4 | on | traced | 406 lines, 22 arguments, brace depth 6. Two-stage screening, Wilson stopping, the closed-form path and an embedded copy of the split routine. | no | Split into named stage functions with the current behaviour pinned by tests before the split. | moves code |
| G-001 | 1 | G | subgroup_consistency_main.R | 376 | subgroup.consistency | S4 | on | traced | 834 lines, 33 arguments, brace depth 8. Argument resolution, seeding, a sequential loop and a separate batched parallel loop with early stopping, the maxeff path, and result assembly in one body; the two evaluation loops duplicate logic. | no | Split into named stage functions with the current behaviour pinned by tests before the split. | moves code |
| F-012 | 2 | F | DESCRIPTION | 57 | (DESCRIPTION Suggests) | S4 | off | traced | Suggests lists DiagrammeR, doRNG and tidyr, which nothing in R/, tests/, vignettes/ or inst/ uses; they appear only under quarto/, which .Rbuildignore excludes. Also at: DESCRIPTION:58, 62 | no | Remove the three packages from Suggests, or keep them deliberately for the quarto material and say so. | add-only |
| D-025 | 2 | D | betaHhat_truth.R | 437 | fs_betaHhat_one | S4 | off | traced | fs_betaHhat_one(), fs_betaHhat_table() and fs_attach_betaHhat() default effect_measure to NULL and do not check it for continuous or count outcomes, so the NULL reaches switch() in compute_aor() and fails with a generic error, whereas the sibling fs_betaHhat_theta_dagger_check() stops with an explicit message. Also at: betaHhat_truth.R:498; betaHhat_truth.R:605; betaHhat_truth.R:905 | no | Apply the same explicit effect_measure requirement in the three functions as in fs_betaHhat_theta_dagger_check(). | add-only |
| B-015 | 2 | B | bootstrap_analysis_dofuture.R | 337 | bootstrap_results and others | S4 | off | traced | 54 local variables in 28 files are assigned and never read (codetools::checkUsageEnv), most in summarize_bootstrap_subgroups.R (9) and fs_oc_predict.R (5). Also at: bootstrap_analysis_dofuture.R:447 outcome_var; bootstrap_summaries_helpers.R:474 perf_color; cox_spline_fit.R:183 stratum_reference; cv_summary_tables.R:389 suffix; dina.R:830 d; dina_bagged.R:484 p; ... (full list in the CSV) | no | Remove the dead assignments after checking each is not a forgotten use. | moves code |
| D-026 | 2 | D | compare_fs_subgroups.R | 91 | compare_fs_subgroups | S4 | off | traced | For subgroup_notation = "benefit" this function labels the subgroups Q/Qc while get_sg_labels(), used by the other table functions, gives G/Gc, and its source note says 'H' under both notations. Also at: compare_fs_subgroups.R:96 (Qc label); compare_fs_subgroups.R:320, 322 (source note hard-codes 'H:'); simulation_tables.R:47-55 (get_sg_labels gives G / Gc for benefit) | no | Take the default labels from get_sg_labels() and use the notation-aware symbol in the source note. | changes behaviour |
| F-013 | 2 | F | cox_ahr_cde_wrapper.R | 374 | cox_ahr_cde_analysis | S4 | off | traced | With save_plots = TRUE the default output_dir = "plots" creates a directory and writes a PDF under the current working directory rather than a user-specified or temporary location. Also at: cox_ahr_cde_wrapper.R:708-716 | no | Require an explicit output_dir when save_plots = TRUE, or default it to tempdir(). | changes behaviour |
| B-016 | 2 | B | cox_spline_fit.R | 183 | cox_cs_fit | S4 | off | traced | Three values are assigned and never read: stratum_reference in cox_cs_fit(), the is_poisson_type field of every .gep_resolve_family() descriptor, and null_ref in glm_cs_fit() whose two branches are both 0. Also at: cox_spline_fit.R:186; glm_effect_profile.R:355, 361, 367, 376, 387, 394, 401; glm_spline_fit.R:420 | no | Remove the unused assignments and field and replace the two-branch null_ref with a constant. | moves code |
| C-057 | 2 | C | dina_bagged.R | 83 | dina_fit_bagged | S4 | off | traced | @return says vcov is 'the raw IJ if project_psd = FALSE', but in that case the code returns the finite-sample-corrected, unprojected matrix; 'raw' is defined two lines later as the uncorrected IJ variance. Also at: dina_bagged.R:87-88 | no | Reword the vcov and vcov_ij_corrected entries to state what each holds under the two flags. | add-only |
| F-014 | 2 | F | forestsearch_helpers.R | 601 | collect_results | S4 | off | traced | The \donttest examples of collect_results() and reset_workers() change the future plan and do not restore the caller's previous plan. Also at: forestsearch_helpers.R:778-782 | no | Capture the previous plan in the two examples and restore it at the end. | add-only |
| A-253 | 2 | A | forestsearch_methods.R | 172 | .fs_print_mr_products | S4 | off | traced | summary() labels the complement SE line 'field-s (Hc) se_field_s' unconditionally, even when the value printed is the unscaled se_field. Also at: forestsearch_methods.R:111 (value chosen), 88-89 and 102 with 163-164 (the bound line carries the correct 'field (unscaled; field-s absent)' tag) | no | Derive the SE line's label from the same scaled flag that sets comp_tag. | changes behaviour |
| A-254 | 2 | A | format_subgroup_summary_tables.R | 222 | format_subgroup_summary_tables | S4 | off | traced | The green highlight rule matches the label ">= 0.95" (with a space) while the producer emits ">=0.95", so the top consistency bin is never coloured although the footnote says green covers >= 0.85. | no | Use the producer's exact label string in the highlight rule. | add-only |
| B-017 | 2 | B | fpr_calibration.R | 425 | print.fpr_calibration | S4 | off | traced | in_ci is computed and never used, and as written it tests fpr_verify against its own confidence interval rather than the corrected FPR against the verification interval. | no | Remove the unused assignment, or compute and print the intended corrected-versus-verification comparison. | add-only |
| A-255 | 2 | A | fs_identification_summary.R | 287 | print.fs_identification | S4 | off | traced | A user-supplied confounders vector that omits the anchor is accepted by fs_identification_summary(), after which print.fs_identification() and fs_identification_table() fail with a subscript error because they index the involvement vector by the anchor name. Also at: fs_identification_summary.R:408; no validation at fs_identification_summary.R:144 and :235 | no | Check at construction that the anchor (or the whole triple) is in confounders, or union it in. | add-only |
| D-027 | 2 | D | fs_oc_family.R | 385 | fs_oc_family_enumerate | S4 | off | traced | The enumeration's size floor keeps a candidate with Pg exactly equal to n.min/n, whereas the search it follows rejects a subgroup of size exactly n.min, so the two use different boundary conventions. Also at: subgroup_search.R:656, :706 (search rejects nx <= n.min); fs_oc_family.R:59 (roxygen states Pg >= n.min / n) | no | Choose one boundary convention for the size floor and state it identically in the enumeration code, its roxygen and the search. | changes behaviour |
| A-256 | 2 | A | fs_oc_grid.R | 147 | fs_oc_grid | S4 | off | traced | For a hand-built family without args_used (the function's own example) the stored floor is numeric(0), and print.fs_oc_grid() silently omits the per-n 'M = ... (floor ...)' line. Also at: fs_oc_grid.R:422 (print.fs_oc_grid, where the line is lost); fs_oc_grid.R:82-93 (the @examples family has no args_used and prints the grid); ... (full list in the CSV) | no | Store NA when the family records no n.min and print the floor only when it is available. | add-only |
| B-018 | 2 | B | fs_oc_predict.R | 179 | fs_oc_predict | S4 | off | traced | fs_oc_predict() repeats inline the pconsistency resolution and check that .fs_oc_resolve_pcons() provides, and the EnH/Esens/Espec/Eppv/Enpv/EbetaH/mass_below formulas exist in three copies. Also at: fs_oc_grid.R:479-484 (.fs_oc_resolve_pcons, same statements); ... (full list in the CSV) | no | Call .fs_oc_resolve_pcons() from fs_oc_predict() and compute the selection-weighted functionals in one shared helper. | moves code |
| C-058 | 2 | C | gg_forest.R | 92 | gg_forest | S4 | off | traced | Three small roxygen statements contradict the code: the documented default widths in gg_forest, the footnote font-size formula in render_forestplot, and the description of the simulation column names in plot_effect_distribution. Also at: render_forestplot.R:36; plot_effect_distribution.R:33 | no | Correct the three roxygen statements. | add-only |
| D-028 | 2 | D | glm_spline_fit.R | 243 | glm_cs_fit | S4 | off | traced | The three sibling profile functions use different conventions for the same arguments: the lower end of the grid, the defaults of z_by, alpha and z_quantile, the argument naming style, and the names of the returned estimate vectors. Also at: cox_spline_fit.R:208, :113-118; glm_effect_profile.R:184, :80-90, :297-302; glm_spline_fit.R:139-142, :449-453 | no | Align the grid convention, defaults and names across the three functions, or document the differences in each help page. | changes behaviour |
| F-015 | 2 | F | globals.R | 287 | (roxygen @importFrom) | S4 | off | traced | Ten imported symbols are never used in R/: ggplot2 annotation_logticks and facet_wrap, graphics barplot, stats integrate, optim, qt and rgamma, utils read.table, sessionInfo and write.table. Also at: gg_forest.R:111; plot_sg_glm_outcomes.R:90; forestsearch-package.R:9; globals.R:288, 290, 292, 293; truefind_asymptotic.R:19; find_k_inter_main.R:85 | no | Drop the unused symbols from the @importFrom tags and regenerate NAMESPACE. | add-only |
| B-019 | 2 | B | guohe_adaptive_r.R | 51 | .ar_coef_se | S4 | off | traced | .ar_coef_se() is a statement-for-statement copy of .g3_coef() that additionally returns the standard error, so the slow-path estimability guards and model calls exist in two files. Also at: guohe_algorithm3.R:51-95 | no | Have one internal slow-path scorer return both the coefficient and its standard error and use it from both files, or extend the header note to say the slow-path copy is also deliberate. | moves code |
| D-029 | 2 | D | mrct_simulation.R | 926 | validate_mrct_data | S4 | off | traced | Sibling functions name the same inputs differently: the case-study data frame is df.case in validate_mrct_data() but df_case in create_dgm_for_mrct(), and the trial size is n_sample in mrct_region_sims() but n in mrct_case_study(). Also at: mrct_simulation.R:1160 (df_case), :263 (n_sample); mrct_case_study.R:241 (n), :162-164 | no | Align the argument names, keeping the old ones as deprecated aliases if needed. | changes behaviour |
| C-059 | 2 | C | mrct_trial_analysis.R | 41 | mrct_trial_analysis | S4 | off | traced | The roxygen for mr_inference_args says to give `seed` "for a reproducible result", which implies the result is otherwise not reproducible, but without a seed the MR step is seeded with forestsearch()'s fixed seedit default, so it is reproducible either way. Also at: forestsearch_main.R:3912 and :1636 (seed falls back to seedit = 8316951); fs_mr_inference.R:615 (set.seed(seed)); mrct_case_study.R:294-297 (the sibling injects seed + k when none is given) | no | Reword the sentence to say that the multiplier draws default to forestsearch()'s seedit and that `seed` overrides it. | add-only |
| D-030 | 2 | D | mrct_trial_analysis.R | 85 | mrct_trial_analysis | S4 | off | traced | Sibling wrappers name the same thing differently: a forestsearch fit is the argument x, fs or fit depending on the function, the region column is region.name in mrct_trial_analysis()/mrct_four_sets() but region_var in mrct_case_study(), and fs_declaration_table() defaults alpha to 0.10 while the function it wraps defaults to 0.05. Also at: mrct_four_sets.R:47 (region.name, truth.name) vs mrct_case_study.R:241 (region_var); ... (full list in the CSV) | no | Document the naming differences (or converge on one first-argument name for a fit in a future API pass) and cross-reference the two alpha defaults in each help page. | add-only |
| B-020 | 2 | B | oc_analyses.R | 332 | compute_cde_glm | S4 | off | traced | The same CDE quantities are computed by near-duplicate helpers: two GLM versions and three places that write out the survival CDE ratio. Also at: generate_glm_dgm.R:633 (.compute_cde_glm, same mean-odds ratio); ... (full list in the CSV) | no | Route the DGM constructors and the simulation extraction through one CDE helper per outcome family. | moves code |
| C-060 | 2 | C | oc_analyses.R | 414 | summarize_simulation_results | S4 | off | traced | Two roxygen examples, one @param description and one emitted sentence describe behaviour the code no longer has. Also at: oc_analyses.R:614 (example comment says format_oc_results() is called by summarize_simulation_results(); :429-432 calls only summarize_single_analysis()); ... (full list in the CSV) | no | Update the two examples, the n_boots description and the trim-disclosure sentence to match the current code. | add-only |
| C-061 | 2 | C | plot_pareto_combined.R | 55 | plot_pareto_combined | S4 | off | traced | The ci_table_list roxygen says CI bars are drawn 'when provided and the tables match (within tolerance)', but the CI tables are never compared; the first table is drawn unconditionally. Also at: plot_pareto_combined.R:325-327 | no | Reword the @param to say the first table is used as given (or add the comparison the sentence promises). | add-only |
| A-257 | 2 | A | plot_pareto_combined.R | 315 | plot_pareto_combined | S4 | off | traced | The combined plot draws every passing candidate, dominated ones included, and its legend labels all non-selected points 'frontier'. Also at: plot_pareto_combined.R:296-297, 310-312 | no | Label the FALSE level as passing candidates, or flag frontier membership separately. | changes behaviour |
| B-021 | 2 | B | plot_sg_distribution.R | 112 | plot_sg_distribution | S4 | off | traced | requireNamespace() guards for packages listed in Imports (ggplot2, data.table, patchwork, gt) can never fail once the package is loaded, so their stop/warning/fallback branches are unreachable. Also at: plot_effect_distribution.R:154; plot_effect_distribution.R:157; gg_forest.R:154; plot_pareto_combined.R:246; pareto_frontier_table.R:389; plot_sg_glm_outcomes.R:119; create_summary_table.R:865; mrct_simulation.R:1962 | no | Remove the guards for imported packages. | add-only |
| B-022 | 2 | B | run_simulation_analysis.R | 1195 | .run_grf_analysis_gen | S4 | off | traced | The 'GRF function not found' fallback branch is unreachable because both names looked up in the package namespace are unconditional top-level functions of the package. Also at: run_simulation_analysis.R:1189-1193 (the lookup that feeds the dead branch) | no | Call the two GRF functions directly and remove the namespace lookup and its fallback branch. | moves code |
| A-258 | 2 | A | sim_aft_gbsg.R | 869 | .simulate_from_gbsg_dgm_ | S4 | off | traced | When the control-arm size rounds to zero, (n1 + 1):n counts backwards, so the simulated data gain an all-NA row and a duplicated subject labelled as control. | no | Stop when either arm would be empty, or index the control rows with seq_len-style logic. | add-only |
| A-259 | 2 | A | simulation_tables.R | 392 | build_classification_table | S4 | off | traced | Analysis labels become data-frame column names with name checking on, so the supported hyphenated method labels (e.g. "fs-effMaxSG") appear in the table header as "fs.effMaxSG". | no | Build the row data frame with name checking off so labels are kept verbatim. | add-only |
| A-260 | 2 | A | simulation_tables.R | 436 | build_classification_table | S4 | off | traced | The bolding documented as flagging an inflated Type I error is applied to every any(H) row, so in alternative-hypothesis sections the detection rate (power) above 0.05 is bolded as well. Also at: simulation_tables.R:1541-1546 (render_reference_table, same rule; there the roxygen at :1485-1486 describes the rule literally) | no | Restrict the bolding to rows of null-hypothesis scenarios, or reword the bold_threshold documentation to say it applies to all any(H) rows. | changes behaviour |
| A-261 | 2 | A | simulation_tables.R | 1131 | interpret_estimation_table | S4 | off | traced | The zero-detection sentence is completed with paste0() using a sprintf-style '%%' escape, so the returned and printed text contains a literal '0%% detection rate'. | no | Use a single percent sign in the pasted literal. | changes behaviour |
| B-023 | 2 | B | subgroup_search.R | 434 | get_combinations_info | S4 | off | traced | The exported get_combinations_info() repeats the enumeration done by the internal generate_combination_indices() that the search actually uses, has no caller anywhere in the package, and differs from it in limits (no maxk <= 3 stop) and in return shape. Also at: subgroup_search.R:379-407 (generate_combination_indices); subgroup_search.R:412-417 (calculate_max_combinations) | no | Keep one enumeration helper and make the exported function a thin wrapper around it (or retire the exported one through the usual deprecation route). | moves code |
| C-062 | 2 | C | test_hte_crump.R | 198 | test_constant_cate | S4 | off | traced | The @return text calls K the number of covariates and says NULL is returned when fewer than two covariates are retained, but K counts basis columns including the intercept and NULL is returned only when no non-intercept column remains. Also at: test_hte_crump.R:201-202 | no | Describe K as the number of basis terms including the intercept and state that NULL is returned when no covariate term is retained. | add-only |
| D-031 | 2 | D | truefind_asymptotic.R | 629 | find_required_sample_size | S4 | off | traced | find_required_sample_size() returns three different list shapes depending on the exit taken (the documented four fields, four plus note, or six fields) and fixes the integration method to "cubature" with no method argument, unlike its sibling compute_detection_probability(). Also at: truefind_asymptotic.R:639, truefind_asymptotic.R:677, truefind_asymptotic.R:619 | no | Return one fixed set of fields from every exit and expose the integration method as an argument. | changes behaviour |
| G-003 | 2 | G | create_summary_table.R | 75 | create_summary_table | S4 | off | traced | 819 lines, 27 arguments, brace depth 8. Data assembly, statistics and gt formatting in one function. | no | Split into named stage functions with the current behaviour pinned by tests before the split. | moves code |
| G-005 | 2 | G | mrct_simulation.R | 260 | mrct_region_sims | S4 | off | traced | 484 lines, 20 arguments, brace depth 8. Plan setup, DGM simulation, analysis, MR and per-replicate record assembly in one function. | no | Split into named stage functions with the current behaviour pinned by tests before the split. | moves code |
| G-008 | 2 | G | plot_sg_weighted_km.R | 94 | plot_sg_weighted_km | S4 | off | traced | 372 lines, 19 arguments, brace depth 5. Subgroup extraction, KM computation, layout and annotation in one function. | no | Split into named stage functions with the current behaviour pinned by tests before the split. | moves code |
| G-004 | 2 | G | plot_subgroup_results_forestplot.R | 195 | plot_subgroup_results_forestplot | S4 | off | traced | 1083 lines, 28 arguments, brace depth 6. Row construction, effect fitting, CV and bootstrap text, theme and forestploter drawing in one function. | no | Split into named stage functions with the current behaviour pinned by tests before the split. | moves code |
| G-007 | 2 | G | run_simulation_analysis.R | 389 | run_simulation_analysis | S4 | off | traced | 377 lines, 28 arguments, brace depth 5. Simulation, ForestSearch and GRF runs, and estimate extraction for survival and GLM outcomes in one function. | no | Split into named stage functions with the current behaviour pinned by tests before the split. | moves code |
| A-262 | 2 | A | plot_pareto_frontier.R | 168 | plot_pareto_frontier | S4 | off | suspected | The colour scale pairs named values with unnamed positional labels c("frontier", "selected"), so when the only plotted level is TRUE (a one-row frontier that is the selected subgroup) the legend is expected to label the selected point 'frontier'. | no | Name the labels by level as plot_pareto_combined() does. | changes behaviour |
| B-024 | 2 | B | summarize_bootstrap_results.R | 387 | summarize_bootstrap_results | S4 | off | suspected | Two defensive branches cannot fire for objects the package itself produces: the 'maxk' attribute lookup and its 'maxk == 0' fallback in summarize_bootstrap_results(), and the 'Factor names' line in summary.forestsearch(). Also at: summarize_bootstrap_results.R:390; forestsearch_methods.R:489-491 | no | Either set the maxk attribute where the bootstrap result is built and store the technical factor names separately from the labels, or delete the two branches. | changes behaviour |

## 6. Notes per finding (S1 and S2)

- **A-001** (S1, on path, reproduced) `forestsearch_main.R:2766` forestsearch
  - Reproducer: Appendix V010, Appendix V002 (verdict demonstrated; wall clock in the appendix).
  - Merged duplicate reports: 1.
  - Configuration (all documented, none the package default): outcome_type = "continuous" (formal :1677-1678, default "survival"), subgroup_method = "dina" (formal :1647, default "consistency"; roxygen :489-512) or use_dina = TRUE (formal :1642, default FALSE); effect.threshold = NULL and hr.threshold = 1.25 left at their formals (:1668, :1670).
  - forestsearch_main.R:1805 user_set_threshold = FALSE; :1811 is the only assignment to hr.threshold in the function and does not fire (effect.threshold is NULL), so hr.threshold stays 1.25.
  - :2275 effect_threshold <- hr.threshold; :2332-2334 the MD branch remaps only the LOCAL effect_threshold <- 0.0 (comment :2333: 1.25 'is meaningless for MD').
  - :2436-2439 threshold_config$screening = 0; :2532-2537 writes effect.threshold = 0 into args_call_all (hr.threshold stays 1.25 there, captured at :1881); :2548-2552 admission_resolved$effect_floor = 0 (forestsearch_helpers.R:2482).
  - :2761-2766 .forestsearch_dina_select(..., hr.threshold = hr.threshold, ...) passes the raw 1.25, not effect_threshold.
  - forestsearch_helpers.R:1492 family = .map_dina_family("continuous") = "gaussian" (:1013); :1525 m_diff <- hr.threshold = 1.25; :1539 tau_sign = -1 under the continuous default adverse_outcome = FALSE (main :2198-2200), +1 when TRUE; :1590-1602 dina_subgroup(m_diff = 1.25, tau_sign); dina_subgroup.R:749 and :888 drop every candidate whose oriented mean tau-hat < 1.25.
  - Effect in the fit: the qualifying family is {mean oriented tau-hat >= 1.25 outcome units}. If no candidate reaches 1.25 the fit returns sg.harm = NULL (helpers :1643-1651); otherwise the default dina_select_statistic = "effect" re-ranking (helpers :1259-1261) applies the admission floor 0 only inside the family that already survived 1.25. The documented floor for this call is 0.0 (resolved-defaults table main :1268; main :446-449 and dina_subgroup.R:144-146 'from the resolved effect.threshold').
  - What the user sees: the threshold banner is suppressed on this route (main :2410, subgroup_method != "dina"), and the 'Harm floor: m_diff = 1.2500' line is printed only with details = TRUE (helpers :1541-1553; default details = FALSE). The returned object carries threshold_config$screening = 0 and admission$effect_floor = 0, so nothing in a default run shows the 1.25.
  - Effect in a replay: bootstrap_analysis_dofuture.R:406, :558, :614 and forestsearch_cross_validation.R:345/:479, :859/:1001 do.call(forestsearch, <args_call_all>) with effect.threshold = 0 -> main :1811 hr.threshold <- 0 -> m_diff = 0 in every replicate / fold. .fs_resolve_threshold_pair returns early for dina (main :97-100), so the disagreeing spellings (0 vs 1.25) raise nothing. The bootstrap therefore bias-corrects a selection made under floor 1.25 with replicates selected under floor 0 (it runs only when the parent found a subgroup, bootstrap_dofuture_main.R:236-238); CV folds use floor 0 even when the parent fit, under 1.25, found nothing.
  - Screening route: main :3237 m_diff_sel <- hr.threshold for gaussian, so with use_dina = TRUE (selected_only = TRUE default, helpers :1113) the single DINA cut is contributed only if mean tau-hat >= 1.25 in the fit, floor 0 in replays.
  - Refutation attempts that failed: no reassignment of hr.threshold between :1811 and :2766 / :3237; .validate_outcome_threshold_config (main :206-310) has no MD check; .dina_assert_ratio_estimand cannot fire on gaussian (helpers :1448-1451); no test exercises the unset-threshold case (tests/testthat/test-directive-c.R:165-173 hr.threshold = 0.30, :463-466 0.10, :694-700 effect.threshold = 0.30; test-dina-proposal-orientation.R:97 hr.threshold = 30; test-directive-c.R:113 pins only the text of the derivation line); NEWS has no entry making 1.25 the intended MD floor under DINA, and NEWS.md:520-541 promises that replicates resolve the thresholds the fit resolved.
  - Not affected: ratio families (log(hr.threshold) equals the resolved value); RD/IRD (refused); any MD run that supplies effect.threshold or hr.threshold.
  - Note: Second verifier, configuration-and-consequence lens. I could not refute it, and I agree with the first verifier on the defect, the location (2766 / 3237 / helpers 1525 all re-checked), the trace and the fix. All three members are one defect. Configuration. Supported and documented: yes. subgroup_method = "dina" is a documented value (roxygen :489-512), continuous is a documented outcome_type, the roxygen states the gaussian floor is derived 'from the resolved effect.threshold' (:446-449), and the resolved-defaults table gives MD c1 = 0.0 (:1268). Default or commonly used: no. It needs a ...
- **A-002** (S1, on path, reproduced) `forestsearch_main.R:3722` forestsearch
  - Reproducer: Appendix V001, Appendix V003, Appendix V004, Appendix R01_bootstrap_id_merge, Appendix R02_bootstrap_end_to_end (verdict demonstrated; wall clock in the appendix).
  - Merged duplicate reports: 2.
  - Configuration: forestsearch() with its default subgroup_method (forestsearch_main.R:1648 formal c("consistency","dina","grf"), match.arg at :1733, roxygen :489 '"consistency" (default)') followed by the exported forestsearch_bootstrap_dofuture(fs.est, nb_boots) (bootstrap_dofuture_main.R:211). No non-default argument is needed.
  - bootstrap_dofuture_main.R:478-480: bootstrap_results(fs.est, df_boot_analysis = fs.est$df.est, ...).
  - bootstrap_analysis_dofuture.R:393-396: only fs.est$confounders.candidate (the q-cut columns, forestsearch_main.R:3945) and treat.recommend are dropped; the id column stays in dfnew_predict (original rows, unique ids).
  - bootstrap_analysis_dofuture.R:438-439: df_boot <- df_boot_analysis[in_boot, ] keeps the original id values (subject i appears k_i times); id_boot is created but id.name is never changed (no assignment to id.name in bootstrap_analysis_dofuture.R or bootstrap_dofuture_main.R).
  - bootstrap_analysis_dofuture.R:550, 558-560, 614: do.call(forestsearch, args) with df.analysis = resample, df.predict = original data, id.name inherited from fs.est$args_call_all.
  - forestsearch_main.R:2725 and :3338: df is the resample (plus cut columns from get_FSdata); id.name is used in forestsearch_main.R only at :3083/:3100 (GRF args), :3390 and the merge block; no uniqueness check, no id reset (get_fsdata.R does not touch the id).
  - forestsearch_main.R:3390-3391: id <- df[, id.name] (with repeats) -> df.fs; :3569 df = df.fs into subgroup.consistency; subgroup_consistency_main.R:1140-1162 -> sg_consistency_out -> extract_subgroup (subgroup_consistency_helpers.R:939-942): df_flag = data.frame(id = df$id, treat.recommend), n rows with repeated ids.
  - forestsearch_main.R:3708, 3722: merge(df, temp, by id, all.x = TRUE): both sides carry k_i rows per id, base merge returns every match, so k_i^2 rows per subject, sum(k_i^2) rows (about 2n) instead of n.
  - forestsearch_main.R:3733: merge(df.predict, temp, by id, all.x = TRUE): an in-bag subject gets k_i rows, an out-of-bag subject one row with treat.recommend = NA.
  - bootstrap_analysis_dofuture.R:629-630 read both frames; :651 / :667 subset(df_PredBoot, treat.recommend == 0 / 1) drops the NA rows, so df_Hstar / df_Hcstar are exactly the bootstrap sample restricted to H* / Hc*; Hstar_obs (:688-697) and Hcstar_obs (:741-750) are therefore bootstrap-sample estimates, not the original-data estimates stated in the comments (:650, :686) and roxygen (:183, :193-194; man/forestsearch_bootstrap_dofuture.Rd:133).
  - bootstrap_analysis_dofuture.R:700-711 / :753-764 fit Hstar_star / Hcstar_star on the k_i^2-row frame; get_Cox_sg (Cox_estimation_helpers.R:33-82) is a plain coxph on the rows given (no weights, no de-duplication), so these are k_i^2-weighted estimates.
  - bootstrap_analysis_dofuture.R:726-728, 768-769: H_biasadj_1/2 and Hc_biasadj_1/2 use (Hstar_star - Hstar_obs) = (k^2-weighted) minus (k-weighted) estimate on the same in-bag subjects instead of (bootstrap) minus (original).
  - bootstrap_dofuture_main.R:573-601: get_dfRes() turns results$H_biasadj_1/2 and Hc_biasadj_1/2 into the reported bias-corrected estimates and their intervals.
  - bootstrap_analysis_dofuture.R:872-875, 889: the replicate's N_sg / E_sg are counted on the k_i^2-row frame and hr_sg is exp of the k^2-weighted estimate.
  - Not affected: the top-level (non-bootstrap) forestsearch() result, because df.predict defaults to df (:3368) and ids are unique there; the DINA and GRF selection branches, which evaluate the rule row-wise with get_dfpred (forestsearch_helpers.R:1664-1666, 1987-1994); df.test (forestsearch_main.R:3737-3742).
  - Note: Second verifier (configuration and consequence lens): I traced it independently and AGREE with the first verifier on verdict, line (3722 first, 3733 second), severity and scope. Configuration: this is the package default, not an edge. subgroup_method defaults to "consistency" (formal :1648, roxygen :489) and the exported forestsearch_bootstrap_dofuture() needs only fs.est and nb_boots; its only call to bootstrap_results() passes fs.est$df.est, and bootstrap_results() is the only place in R/ that calls forestsearch() with df.predict different from df.analysis (grep of df.predict call sites: ...
- **A-003** (S1, off path, reproduced) `run_simulation_analysis.R:1025` .extract_fs_estimates_gen
  - Reproducer: Appendix V005 (verdict demonstrated; wall clock in the appendix).
  - Configuration: run_simulation_analysis() is exported (NAMESPACE:278); formals run_fs = TRUE, run_fs_grf = TRUE, run_grf = TRUE, harm_col = "flag_harm" (run_simulation_analysis.R:403-406); @return documents 'classification metrics against the DGM-stored subgroup (sens, spec, ppv, npv)' (:302-303). Default arguments, any alternative DGM.
  - Reach: :591 / :616 / :678 call .run_fs_analysis_gen(); forestsearch() with no subgroup returns out_sg = NULL and sg.harm = NULL (subgroup_consistency_main.R:1079-1082, forestsearch_main.R:3514-3518), an error returns NULL (:818-821); has_result is FALSE (:823-831); :844-845 pass fs_res = NULL, fs_full = NULL.
  - Defect: .has_treat_rec FALSE (:1005-1007) -> branch :1009; flag_harm has members under an alternative (generate_aft_dgm_helpers.R:601-638; zero only under model = 'null', :597) -> :1020 TRUE; :1023-1026 store sens = 0, ppv = 0, spec = sum(!true_H)/nrow(df), npv = 1; return at :1037. No later assignment to spec or npv in the file other than :1031-1032, :1066-1068, :1127-1129, :1392-1394, :1418-1420, :1531-1533.
  - Correct value: same function, :1123-1129, spec = tn/(tn+fp), npv = tn/(tn+fn); an empty identified set has tn = |Hc|, fp = 0, fn = |H|, so spec = 1 and npv = |Hc|/N. Specificity of an empty set needs no convention (tn + fp = |Hc| > 0).
  - Same assignments at :1066 (rows but no M.* factor), :1392 (grf_est NULL), :1418 (GRF returns no subgroup).
  - Effect on output: build_classification_table() (exported) uses all rows of the method (simulation_tables.R:358), fills only NA values (:369-372, so the fill spec = 1 no longer applies) and reports sens(Hc) = mean(r$spec) (:379); summarize_single_analysis() means over all rows (oc_analyses.R:455-457); format_oc_results() likewise (oc_analyses.R:688). Reported specificity = det x E[spec | detected] + (1 - det) x |Hc|/N in place of + (1 - det) x 1.
  - Null DGMs are not affected (:1030-1032 store 1 and 1). forestsearch(), bootstrap and MR results are not affected. .fs_oc_identification() overwrites non-detected rows with spec = 1 (fs_mr_oc_summary.R:286-290, :308) and reads a different result schema.
  - Note: Second verifier (configuration and consequence). I agree with the first verifier: confirmed, S1, traced; both members are one defect. Configuration: default arguments of the exported driver (formals run_fs = TRUE, run_fs_grf = TRUE, run_grf = TRUE, harm_col = "flag_harm"); the workflow run_simulation_analysis() -> format_oc_results() / build_classification_table() is the one shown in vignettes/articles/paper_simulations.qmd (Steps 3-4, lines 205-260) and in README. Every non-detected replicate under any alternative DGM reaches the branch, on the legacy boolean path and on the methods path ...
- **A-004** (S1, off path, reproduced) `synthetic_data_perturbation.R:366` generate_gbsg_bootstrap_general
  - Reproducer: Appendix R04_gbsg_data_call (verdict demonstrated; wall clock in the appendix).
  - synthetic_data_perturbation.R:365 gbsg <- NULL; :366 data("gbsg", package = "survival", envir = environment()) warns 'data set gbsg not found' and assigns nothing.
  - :368-369 generate_bootstrap_synthetic(data = gbsg, ...) receives NULL and stops with "'data' must be a data frame".
  - The function is exported (NAMESPACE:207) with a roxygen example (in dontrun), so the default call is a documented configuration; it is not on the conducted-analysis path.
  - Not raised by the readers (static reading assumed data() works); surfaced by verifier reproducer V185, whose script for another finding stopped here, then confirmed by R04.
  - Note: Found in the main session from a reproducer failure; not through the two-verifier gate. The other in-package examples use survival::gbsg.
- **A-005** (S2, on path, reproduced) `Cox_estimation_helpers.R:54` get_Cox_sg
  - Reproducer: Appendix V043 (verdict demonstrated; wall clock in the appendix).
  - Cox_estimation_helpers.R:54 `fit <- suppressWarnings(do.call(coxph, fit_args))`; :57-82 extract coef and robust SE unconditionally; there is no event-count or convergence check anywhere in :33-83.
  - Sibling guard that get_Cox_sg lacks: glm_effect_estimators.R:303-308 (.make_cox_estimator calls .fs_existence_reason("HR", ...) and returns non-estimable); glm_effect_estimators.R:269-282 (HR needs e0 and e1 >= 1); subgroup_search.R:715.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-006** (S2, on path, reproduced) `bootstrap_analysis_dofuture.R:436` bootstrap_results
  - Reproducer: Appendix V006 (verdict demonstrated; wall clock in the appendix).
  - Trigger (non-default route): direct call of exported bootstrap_results() with formal boot_index_mat = NULL (bootstrap_analysis_dofuture.R:305); the roxygen example (:276-286) calls it this way.
  - :334 set.seed(seed); :408-413 foreach with .options.future seed = TRUE; :436 in_boot <- sample.int(NN, NN, replace = TRUE) inside the body draws from per-iteration L'Ecuyer streams. Indices are not in the returned table (:905-960) or its attributes.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-007** (S2, on path, reproduced) `bootstrap_analysis_dofuture.R:534` bootstrap_results
  - Reproducer: Appendix V045 (verdict demonstrated; wall clock in the appendix).
  - bootstrap_analysis_dofuture.R:464-465 (survival) / :472-473 (GLM) compute events_H_0, events_H_1; :494-495 / :502-503 compute events_Hc_0, events_Hc_1; their only use before the overwrite is the message at :467-470 and :497-500.
  - bootstrap_analysis_dofuture.R:534-537: `events_H_0 <- NA; events_H_1 <- NA; events_Hc_0 <- NA; events_Hc_1 <- NA`, inside the 'initialise as NA' block that sits after the computation; grep shows no later assignment to these four names.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-008** (S2, on path, reproduced) `bootstrap_analysis_dofuture.R:614` bootstrap_results
  - Reproducer: Appendix V046 (verdict demonstrated; wall clock in the appendix).
  - bootstrap_analysis_dofuture.R:614 run_bootstrap <- try(do.call(forestsearch, args_FS_boot), TRUE): the error is caught silently.
  - :619-621 warning("Bootstrap ", boot, " failed: ...") is raised inside the foreach body, which is wrapped in suppressWarnings({...}) opened at :408; the caller also runs the whole section under options(warn = -1) (bootstrap_dofuture_main.R:346-347, restored only at :492).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-009** (S2, on path, reproduced) `bootstrap_dofuture_main.R:302` forestsearch_bootstrap_dofuture
  - Reproducer: Appendix V047 (verdict demonstrated; wall clock in the appendix).
  - forestsearch_main.R:1797-1798 resolves adverse_outcome NULL -> FALSE for outcome_type = "continuous" before the capture at :1880-1881, so args_call_all$adverse_outcome = FALSE; :2220-2229 (and the fast closures at :2244-2258) build the search estimator with adverse_outcome = adverse_outcome; :2641 stores args_call_all$ps_adjust_method; :2705-2712 rebuilds the search estimator with ps_adjust_method when ps_method != "none".
  - bootstrap_dofuture_main.R:302-310 calls make_effect_estimator() with outcome_type, treat.name, outcome.name, event.name, offset.name, effect_measure, adjust_covariates only.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-010** (S2, on path, reproduced) `bootstrap_dofuture_main.R:533` forestsearch_bootstrap_dofuture
  - Reproducer: Appendix V049 (verdict demonstrated; wall clock in the appendix).
  - Trigger: exactly one replicate has a non-NA H_biasadj_* (very small nb_boots, or a subgroup that is almost never re-identified).
  - bootstrap_dofuture_main.R:532-533: n_success == 0L is the only guard, so n_success == 1 proceeds to get_dfRes() at :574-587 (and :597-609 for Hc).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-011** (S2, on path, reproduced) `bootstrap_dofuture_main.R:736` forestsearch_bootstrap_dofuture
  - Reproducer: Appendix V007 (verdict demonstrated; wall clock in the appendix).
  - Configuration: forestsearch(est.scale = "1/hr"). forestsearch_main.R:1637 formal est.scale = "hr" (default); the value is never read by the search (only occurrences in the file: :382 roxygen, :1487 example, :1637 formal) and is not validated; it is captured at forestsearch_main.R:1880-1881 (mget(names(formals()))) into args_call_all. Selection and the coding treat.recommend == 0 = sg.harm are therefore identical on both scales.
  - Entry: forestsearch_bootstrap_dofuture() has no est.scale argument (formals :211-218); it reads the fit's value at bootstrap_dofuture_main.R:267 and :559. No normalisation for survival (:564-571 touches GLM identity measures only).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-012** (S2, on path, reproduced) `consistency_resample.R:81` .consistency_cox_pieces
  - Reproducer: Appendix V050 (verdict demonstrated; wall clock in the appendix).
  - consistency_resample.R:64-66: adjust_covariates = "x2" -> formula Surv(Y, Event) ~ Treat + x2; :68-69 fitted with robust = FALSE.
  - :79 dfbeta <- residuals(fit, type = "dfbeta"): n x 2 matrix. :80-82 j <- match(treat.name, colnames(dfbeta)); if (is.na(j)) return(NULL).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-013** (S2, on path, reproduced) `dina.R:270` dina_fit
  - Reproducer: Appendix V051 (verdict demonstrated; wall clock in the appendix).
  - Trigger A (one fold): a sparse 0/1 covariate whose non-zero rows all fall in one test fold (folds are an unstratified permutation, dina.R:435-437). In the other pass X_test[, j] is all zero, so X_resW[, j] = 0 (dina.R:833 / :895).
  - GLM families: glm.fit returns NA for the aliased coefficient with no warning; dina.R:848 only checks convergence; :852 passes the NA on. Cox: coxph returns NA and solve(fit$var) fails -> warning + NA sigma (dina.R:916-923).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-014** (S2, on path, reproduced) `dina.R:501` .dina_fit_baselines_glm
  - Reproducer: Appendix V008 (verdict demonstrated; wall clock in the appendix).
  - Configuration: dina.R:160 formal `baseline_method = NULL` (default); dina.R:218-220 replaces NULL by "glm"/"cox" only, so the default path never reaches the list branch. The trigger is a NON-default, user-written list of closures.
  - Documented: dina.R:86-89 roxygen, 'a list with components `eta0` and `eta1`, each a function ... For Gaussian, a list with a single component `nu` may be supplied instead' (no family restriction on eta0/eta1; `nu` is described as an optional alternative). dina.R:508-509 stop text says the same ('either `nu` (Gaussian only) or both `eta0` and `eta1`').
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-015** (S2, on path, reproduced) `dina.R:660` .dina_fit_censoring_cox
  - Reproducer: Appendix V052 (verdict demonstrated; wall clock in the appendix).
  - dina.R:648-660: for cens_type 'unknown' a single glm(.S ~ X) is fitted on the pooled training fold (no treatment term, no split by arm) and list(censor0_fn = pred_fn, censor1_fn = pred_fn) is returned; pred_fn ignores eta_new.
  - dina.R:801-805: p0_nc = p1_nc = p, so num = e*p, den = e*p + (1-e)*p + eps = p + eps, a_x = e*p/(p + eps).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-016** (S2, on path, reproduced) `dina_df.R:244` .dina_extract_from_df
  - Reproducer: Appendix V009 (verdict demonstrated; wall clock in the appendix).
  - CONFIGURATION
  - dina_df.R:35-36: `@param df data frame (or tibble / `data.table`) containing the outcome, treatment, and covariate columns.` This one roxygen line is the only place the package promises data.table input for DINA; it dates from the commit that added DINA (81bffbf7). dina_bagged() inherits it (dina_df.R:129).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-017** (S2, on path, reproduced) `dina_subgroup.R:425` dina_subgroup
  - Reproducer: Appendix V053 (verdict demonstrated; wall clock in the appendix).
  - dina.R:289-293: dina_fit() names the coefficients '(Intercept)' + colnames(X); dina_df.R:261 builds X from df[, covariates], so the names are the training covariate names in training order.
  - dina_subgroup.R:408: X <- .dina_extract_X_from_df(df, covariates) takes columns in the order of the CALL's `covariates` (dina_subgroup.R:655-683 validates existence/NA/numeric only).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-018** (S2, on path, reproduced) `dina_subgroup.R:583` dina_subgroup
  - Reproducer: Appendix V054 (verdict demonstrated; wall clock in the appendix).
  - dina_subgroup.R:738: depth-1 thresholds are data values (sort(unique(x_j))), so at least one row equals the winning threshold exactly; depth-2 thresholds are type-7 quantiles (dina_subgroup.R:839-842), which are data values or interpolations.
  - dina_subgroup.R:575-579: the mask uses the exact double (`x_t <= win_q[t]` / `x_t >= win_q[t]`).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-019** (S2, on path, reproduced) `dina_subgroup.R:738` .dina_collect_candidates
  - Reproducer: Appendix V055 (verdict demonstrated; wall clock in the appendix).
  - dina_subgroup.R:738 `qs <- sort(unique(x_j))` includes max(x_j) and min(x_j); dina_subgroup.R:743-744: dir 'left' with q = max, or 'right' with q = min, gives a mask that is all TRUE and n_S = n.
  - dina_subgroup.R:746-749: the only tests are n_S >= n_min and mean_tau >= m_diff; there is no upper size bound or complement check anywhere in the collector or its two callers (grep for a full-sample guard in dina_subgroup.R finds none).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-020** (S2, on path, reproduced) `dina_subgroup.R:1190` dina_frontier
  - Reproducer: Appendix V056, Appendix V072 (verdict demonstrated; wall clock in the appendix).
  - Merged duplicate reports: 1.
  - dina_subgroup.R:329, 435: dina_subgroup() takes tau_sign and negates tau_hat before the floor and the ranking.
  - dina_subgroup.R:1190-1197: dina_frontier() has no tau_sign formal; dina_subgroup.R:1240 computes raw tau_hat; dina_subgroup.R:1243-1245 applies the scope = 'harm' floor to it and dina_subgroup.R:1262, 1282-1284 build the (effect, N) frontier and effect ranking on it (so the direction matters for scope = 'wide' too, not only 'harm').
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-021** (S2, on path, reproduced) `forestsearch_helpers.R:122` get_dfpred
  - Reproducer: Appendix V057 (verdict demonstrated; wall clock in the appendix).
  - forestsearch_helpers.R:93 labels is NULL / character(0); :101 length(labels) == 1L is FALSE, so the disjunction branch is skipped.
  - :122 in_harm <- rep(TRUE, nrow(df.pred)); :124 the for loop over labels has zero iterations; :137 treat.recommend <- ifelse(in_harm, 0L, 1L) gives 0 for every row.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-022** (S2, on path, reproduced) `forestsearch_helpers.R:202` evaluate_comparison
  - Reproducer: Appendix V058 (verdict demonstrated; wall clock in the appendix).
  - forestsearch_helpers.R:199-206: only is.factor(col) triggers as.numeric(as.character(col)); a character column is left unchanged.
  - :207-218: val is numeric; 'col >= val' with a character col makes R compare strings ('10' >= '9' is FALSE).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-023** (S2, on path, reproduced) `forestsearch_helpers.R:1311` .dina_reselect_on_effect
  - Reproducer: Appendix V059 (verdict demonstrated; wall clock in the appendix).
  - forestsearch_main.R:1646-1651: dina_select_statistic and grf_select_statistic default to "effect", grf_selection to "frontier"; injected at :1787.
  - dina_subgroup.R:604-616: native winner fields include labels, label, se_mean_tau_hat, ci.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-024** (S2, on path, reproduced) `forestsearch_helpers.R:1531` .forestsearch_dina_select
  - Reproducer: Appendix V060 (verdict demonstrated; wall clock in the appendix).
  - forestsearch_main.R:2211-2216: effect_measure IRR/IRD stops without offset.name.
  - forestsearch_helpers.R:1011-1017: count maps to DINA family "poisson".
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-025** (S2, on path, reproduced) `forestsearch_helpers.R:1759` .grf_reselect_on_effect
  - Reproducer: Appendix V061 (verdict demonstrated; wall clock in the appendix).
  - Trigger: forestsearch(subgroup_method = "grf", effect_measure = "IRD", offset.name = ...) with the default grf_selection = "frontier", grf_select_statistic = "effect" (forestsearch_main.R:1650-1651); no guard refuses IRD on the GRF path, and .dina_assert_ratio_estimand()'s message (forestsearch_helpers.R:1455-1457) recommends "grf" for IRD.
  - forestsearch_helpers.R:1752-1754: each candidate is scored by .fs_mr_pieces(); fs_mr_inference.R:67-71 dispatches to .consistency_glm_pieces(); consistency_resample.R:249 returns NULL for effect_measure == "IRD"; helpers:1754 skips the candidate.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **C-001** (S2, on path, reproduced) `forestsearch_main.R:382` forestsearch
  - Reproducer: Appendix V011 (verdict demonstrated; wall clock in the appendix).
  - Configuration: forestsearch_main.R:382 roxygen '@param est.scale Character. Estimation scale ("hr" or "rmst"). Default "hr".' (same text man/forestsearch.Rd:150); formal at :1637 is est.scale = "hr". So "rmst" is documented but non-default; no vignette, example or test calls forestsearch() with it (every call in vignettes/ and tests/ uses "hr", one test uses "1/hr").
  - Not implemented: grep of "rmst" as a quoted value over R/ returns only line 382. Introduced by commit f8ca1eae; no NEWS entry.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **INV-3** (S2, on path, reproduced) `forestsearch_main.R:928` forestsearch
  - Reproducer: Appendix V063 (verdict demonstrated; wall clock in the appendix).
  - Merged duplicate reports: 15.
  - All paths under R/. Every reader-claimed mismatch was re-opened; two package-wide sweeps were added (roxygen @param blocks vs formals: 2789 params / 481 functions; man/*.Rd \arguments vs \usage: 402 pages / 2718 args, which also covers @inheritParams and @rdname).
  - A. Formal has a literal default, roxygen states a different value (6):
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-026** (S2, on path, reproduced) `forestsearch_main.R:1711` forestsearch
  - Reproducer: Appendix V064, Appendix V065 (verdict demonstrated; wall clock in the appendix).
  - Merged duplicate reports: 1.
  - Input: forestsearch(use_grf = TRUE, tune_grf = TRUE) or forestsearch(subgroup_method = "grf", tune_grf = TRUE).
  - forestsearch_main.R:1711 formal tune_grf = FALSE; grep finds the name only at :1077 (roxygen) and :1711 in the whole file; it is captured into args_call_all only.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-027** (S2, on path, reproduced) `forestsearch_main.R:2220` forestsearch
  - Reproducer: Appendix V067 (verdict demonstrated; wall clock in the appendix).
  - Input: forestsearch(outcome_type = "count", offset.name = ..., overdispersion = "negbin" or "quasi") with the default use_grf = FALSE (:1639) and subgroup_method = "consistency".
  - forestsearch_main.R:1709 formal; :1731 and :2173 match.arg; the only other uses in the body are :2952 (.forestsearch_grf_select, subgroup_method = "grf") and :3110 (.build_grf_glm_args, use_grf = TRUE).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-028** (S2, on path, reproduced) `forestsearch_main.R:2291` forestsearch
  - Reproducer: Appendix V068 (verdict demonstrated; wall clock in the appendix).
  - Input: outcome_type = "count" (or binary), effect_measure = "IRD", offset.name set, effect.threshold = 1.5 (a rate difference per unit of follow-up time).
  - forestsearch_main.R:1805 user_set_threshold = TRUE; :2180-2190 -> :298-310 warns 'A risk difference cannot exceed 1.0' for both RD and IRD.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-029** (S2, on path, reproduced) `forestsearch_main.R:2656` forestsearch
  - Reproducer: Appendix V069 (verdict demonstrated; wall clock in the appendix).
  - Input: GLM outcome, default is.RCT = TRUE (:1635), ps_hat = <vector>, ps_adjust_method = "iptw", ps_method not supplied.
  - forestsearch_main.R:2596-2599 ps_method <- "none" because is.RCT is TRUE.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-030** (S2, on path, reproduced) `forestsearch_main.R:2705` forestsearch
  - Reproducer: Appendix V012 (verdict demonstrated; wall clock in the appendix).
  - Configuration: NOT a default. Formals forestsearch_main.R:1635 is.RCT = TRUE, :1704 outcome_type = c("survival", ...), :1715 ps_method = NULL; :2597-2599 ps_method resolves to "none" when is.RCT = TRUE, so the branch at :2656 is skipped in every RCT run. Trigger needs all three: GLM outcome, is.RCT = FALSE (ps_method then resolves to "grf") or an explicit ps_method, and a resolved adverse_outcome = FALSE.
  - Documentation of the trigger: roxygen :1064-1066 (adverse_outcome FALSE = default for continuous); :1102-1109 ps_method 'implemented but has not been tested against the manuscript analyses ... Default: "none" for RCT, "grf" for observational'; section :1423-1431 'They are untested, not known to be wrong. Users who need them should check them on their own data'.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-031** (S2, on path, reproduced) `forestsearch_main.R:2875` forestsearch
  - Reproducer: Appendix V070 (verdict demonstrated; wall clock in the appendix).
  - Trigger: mr_inference = TRUE, subgroup_method = "dina" or "grf", and the selected subgroup is dropped by .fs_mr_assemble() (fs_mr_inference.R:95-97: fewer than 6 rows, .fs_mr_pieces() NULL/error, or dfbeta length mismatch; pieces return NULL at consistency_resample.R:62 (<2 events), :249 (effect_measure IRD), :265 (<3 per arm), :305 (non-convergence)).
  - fs_mr_inference.R:631 sel <- match(sel_lab, asm$names) is NA; fs_mr_inference.R:637-641 returns the short list with harm_flag = NA and a note.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-032** (S2, on path, reproduced) `forestsearch_main.R:3071` forestsearch
  - Reproducer: Appendix V071 (verdict demonstrated; wall clock in the appendix).
  - Roxygen :389-393 says cuts are 'taken from grf_res / grf_cuts when supplied'; :399 'grf_res GRF results object (optional, for reuse)'.
  - Case 1, forestsearch(use_grf = TRUE, grf_res = <result with tree.cuts>): :3071 condition use_grf && (is.null(grf_res) || is.null(grf_res$tree.cuts)) is FALSE, so the block holding the only assignment grf_cuts <- grf_res$tree.cuts (:3124) is skipped and there is no else branch; grf_cuts keeps its formal NULL (:1641).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-033** (S2, on path, reproduced) `forestsearch_main.R:3350` forestsearch
  - Reproducer: Appendix V073, Appendix V084 (verdict demonstrated; wall clock in the appendix).
  - Merged duplicate reports: 1.
  - forestsearch_main.R:1920 and get_fsdata.R:99 test only is.data.frame(), TRUE for tbl_df and data.table; grep finds no as.data.frame() in forestsearch_main.R or on the get_FSdata() path; get_FSdata() uses [[ ]] throughout and returns df.FS with its class intact (get_fsdata.R:107, :746).
  - forestsearch_main.R:3338 df <- FSdata$df; :3350 Y <- df[, outcome.name], :3352 Event <- df[, event.name], :3358 Treat <- df[, treat.name], :3390 id <- df[, c(id.name)]: for a tibble each is a one-column tibble, not a vector.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-034** (S2, on path, reproduced) `forestsearch_main.R:3391` forestsearch
  - Reproducer: Appendix V021, Appendix V074 (verdict demonstrated; wall clock in the appendix).
  - Merged duplicate reports: 2.
  - Input: outcome_type binary/continuous/count, ps_method != "none", ps_adjust_method = "iptw" (or "dr_gcomp").
  - forestsearch_main.R:2677-2685 (user ps_hat) or :2688-2695 (estimated) append ps_hat, sw, ips_covar to df.analysis; :2704-2712 rebuild estimator_fn with ps_adjust_method.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-035** (S2, on path, reproduced) `forestsearch_main.R:3399` forestsearch
  - Reproducer: Appendix V075 (verdict demonstrated; wall clock in the appendix).
  - Trigger: outcome_type = "binary" (or continuous/count), outcome.name = "Event".
  - forestsearch_main.R:3350 Y <- df[, outcome.name] (real outcome); :3353-3356 Event <- rep(1L, nrow(df)) for every non-survival outcome; :3391 df.fs <- data.frame(Y, Event, Treat, id, df.confounders) -> df.fs$Event is constant 1.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-036** (S2, on path, reproduced) `forestsearch_main.R:3545` forestsearch
  - Reproducer: Appendix V013 (verdict demonstrated; wall clock in the appendix).
  - Configuration: sg_focus = "maxeff" is a supported, documented, NON-default focus (formal forestsearch_main.R:1672 sg_focus = "hr"; roxygen :699-701 'The effect maximiser over the enumerated family, with no consistency threshold, effect threshold, early stopping ...', :720 '"maxeff" disables both', :745 'effect argmax, no floors').
  - Trigger: forestsearch(sg_focus = "maxeff", subgroup_method = "consistency") on data where every estimable candidate has effect <= c2. With default thresholds (hr.consistency = 1.0, :1671) that needs every candidate HR <= 1 (strong uniform benefit or a small family); with hr.threshold supplied alone, c2 is derived as 0.80 * c1 (:134-135, called at :1827-1833, not focus-conditional), e.g. hr.threshold = 1.5 puts the gate at 1.2.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-037** (S2, on path, reproduced) `forestsearch_main.R:3568` forestsearch
  - Reproducer: Appendix V076 (verdict demonstrated; wall clock in the appendix).
  - Trigger: forestsearch(..., seedit = <anything but 8316951>) on a path that draws splits: consistency_method = "split", or the automatic fallback to splitting (IRD, propensity-adjusted effects, non-convergent resample fit).
  - forestsearch_main.R:3568-3592 consistency_overrides (plus GLM additions :3598-3626) contains no `seed` entry.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-038** (S2, on path, reproduced) `forestsearch_main.R:3722` forestsearch
  - Reproducer: Appendix V014 (verdict demonstrated; wall clock in the appendix).
  - Merged duplicate reports: 1.
  - Trigger (non-default data layout, not a setting): df.analysis whose rows are not in ascending order of id.name (shuffled rows, or character ids such as P1, P2, ..., P10, which sort lexicographically).
  - subgroup_consistency_helpers.R:934-942 sg.harm.id and df_flag are built in df.fs row order (= df.analysis order, forestsearch_main.R:3390-3391).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-039** (S2, on path, reproduced) `forestsearch_main.R:3830` forestsearch
  - Reproducer: Appendix V077 (verdict demonstrated; wall clock in the appendix).
  - Trigger: forestsearch(outcome_type = "survival", adjust_covariates = "strata(x1)", mr_inference = TRUE); the term form is documented at forestsearch_main.R:1088-1094.
  - :3416-3434 only the raw variables (.fs_adjust_vars -> "x1", Cox_estimation_helpers.R:158-166) are attached to df.fs; the full terms go verbatim to the search (:3481-3484) and the consistency engine (:3591).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-040** (S2, on path, reproduced) `forestsearch_main.R:3844` forestsearch
  - Reproducer: Appendix V078 (verdict demonstrated; wall clock in the appendix).
  - Merged duplicate reports: 1.
  - Trigger: forestsearch(..., mr_inference = TRUE, mr_inference_args = list(field_R_out = 50, field_R_in = 25)) with subgroup_method = "consistency"; ci_method defaults to "field" on this branch (forestsearch_main.R:3870).
  - :3844-3912 the direct fs_mr_inference() call lists every other control but not field_R_out / field_R_in (grep: the strings do not occur in forestsearch_main.R).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-041** (S2, on path, reproduced) `fs_mr_inference_methods.R:157` .fs_apply_mr
  - Reproducer: Appendix V081 (verdict demonstrated; wall clock in the appendix).
  - Merged duplicate reports: 1.
  - fs_mr_inference.R:584-609: 29 formals; the last four are keep_declaration_field, keep_field_matrix, declaration_c0, pconsistency.digits.
  - fs_mr_inference_methods.R:157-190: the wrapper's call ends at field_recovery; none of the four is named, and there is no ... or do.call pass-through.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-042** (S2, on path, reproduced) `get_FSdata_helpers.R:337` lasso_selection
  - Reproducer: Appendix V082 (verdict demonstrated; wall clock in the appendix).
  - get_FSdata_helpers.R:333-338: coefs come from coef(fit) with its row names; selected <- names(coefs)[coefs != 0]. For glmnet families other than 'cox' (set at :283-303) coef() has an '(Intercept)' row, which is non-zero in general, so it is always in `selected`.
  - get_fsdata.R:336 lassokeep <- get_lasso$selected; :354 `if (length(lassokeep) > 0)` is therefore always TRUE for GLM outcome types.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-043** (S2, on path, reproduced) `get_FSdata_helpers.R:692` collapse_redundant_cuts
  - Reproducer: Appendix V015 (verdict demonstrated; wall clock in the appendix).
  - Configuration: every argument at its formal default (forestsearch_main.R:1638-1665: use_lasso = FALSE, use_grf = FALSE, cut_type = "default", cont.cutoff = 4, collapse_cuts = TRUE; get_fsdata.R:84-91 has the same defaults). The only non-standard element is the INPUT TYPE: a confounder stored as a factor with >= 4 distinct values whose levels parse as numbers (e.g. factor(grade) with grades 1-4, factor(ecog) 0-3).
  - Call chain from the exported entry point: forestsearch() does not coerce df.analysis before data preparation (the only coercion, .coerce_covariates_numeric at forestsearch_main.R:3195, writes to a DINA-local copy df.dina); forestsearch_main.R:3310-3319 calls get_FSdata via filter_call_args(args_call_all, ...) with df.analysis unchanged and collapse_cuts = TRUE carried from the formals.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-044** (S2, on path, reproduced) `get_FSdata_helpers.R:859` evaluate_cuts_once
  - Reproducer: Appendix V083 (verdict demonstrated; wall clock in the appendix).
  - No NA check on confounders: get_fsdata.R:99-195 validates only class/type; forestsearch_main.R has no NA guard before the call at :3310-3319.
  - get_FSdata_helpers.R:11: is.continuous() uses length(unique(x)) >= cutoff, and unique() keeps NA, so x in {1,2,3} with some NA has 4 unique values and is classed continuous at get_fsdata.R:272-279 (cont.cutoff = 4).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-045** (S2, on path, reproduced) `get_fsdata.R:291` get_FSdata
  - Reproducer: Appendix V085 (verdict demonstrated; wall clock in the appendix).
  - forestsearch_main.R:581-588 and get_fsdata.R:31-32 document conf.cont_medians as confounders 'to cut at their median' and conf.cont_medians_force as 'additional' confounders to force a median cut on.
  - get_fsdata.R:291-293: all continuous confounders are queued for cutting only when BOTH arguments are NULL.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-046** (S2, on path, reproduced) `get_fsdata.R:486` get_FSdata
  - Reproducer: Appendix V016 (verdict demonstrated; wall clock in the appendix).
  - Configuration: forestsearch() formals use_lasso = FALSE, use_grf = FALSE, cont.cutoff = 4 (forestsearch_main.R:1638-1639, 1660); roxygen for confounders.name is only 'Character vector. Names of candidate subgroup-defining variables.' (forestsearch_main.R:350) and for cont.cutoff 'Cutoff for continuous vs categorical. Default 4.' (forestsearch_main.R:580). No coding requirement is stated anywhere in roxygen.
  - Call chain: forestsearch() does not touch the confounder columns before Section 4 (the only coercion, .coerce_covariates_numeric, is on the DINA paths, forestsearch_main.R:3195 and forestsearch_helpers.R:1517); it calls get_FSdata through filter_call_args with df.analysis, grf_cuts, dina_cuts and details = FALSE overridden (forestsearch_main.R:3310-3320; filter_call_args at summary_utility_functions.R:679-688 passes confounders.name and cont.cutoff through unchanged).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-047** (S2, on path, reproduced) `get_fsdata.R:499` get_FSdata
  - Reproducer: Appendix V086 (verdict demonstrated; wall clock in the appendix).
  - Input: use_grf = TRUE (non-default), binary confounder 'bm', continuous 'bmi', GRF cut 'bmi <= 25.3' (format from grf_helpers.R:249; same for z1/z10, er/her2).
  - get_fsdata.R:286-287, :486: confs = 'bm' (under cut_type = 'default' no median strings exist at this point).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **D-001** (S2, on path, reproduced) `get_fsdata.R:517` get_FSdata
  - Reproducer: Appendix V087b (verdict demonstrated; wall clock in the appendix).
  - Other reproducer runs: Appendix V087 (no verdict line (exit 1)).
  - Merged duplicate reports: 1.
  - Not a violation of invariant 1 as worded (no model is fitted by the call), but the one place where the front-end flag does not gate the candidate pool.
  - Trigger: forestsearch(..., use_dina = FALSE (default), dina_cuts = <character vector>), e.g. the $dina_cuts returned by an earlier fit (returned at forestsearch_main.R:3960).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-048** (S2, on path, reproduced) `get_fsdata.R:527` get_FSdata
  - Reproducer: Appendix V088 (verdict demonstrated; wall clock in the appendix).
  - No validation of conf_force in forestsearch() (forestsearch_main.R:2728-2729 only prints it) or at get_FSdata() entry.
  - Input conf_force = 'age <= 50 & er > 0'.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-049** (S2, on path, reproduced) `get_fsdata.R:535` get_FSdata
  - Reproducer: Appendix V089 (verdict demonstrated; wall clock in the appendix).
  - get_fsdata.R:535-541: exclusion is a grepl() filter on the pre-coarsening strings (e.g. 'x <= 49.6' does not match the pattern 'x <= 50').
  - get_fsdata.R:550-571: collapse_redundant_cuts() runs afterwards (collapse_cuts = TRUE is the formal default, :88); only GRF/DINA/user conf_force cuts are protected (:570).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-050** (S2, on path, reproduced) `get_fsdata.R:634` get_FSdata
  - Reproducer: Appendix V090 (verdict demonstrated; wall clock in the appendix).
  - get_fsdata.R:573-603: the 'no candidate cut expressions' stop() tests n_confs before evaluation only.
  - get_fsdata.R:618-638: evaluate_cuts_once() marks single-valued cuts invalid (get_FSdata_helpers.R:859); they are dropped at :634-645; the only report is cat() under details (:636-638).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-051** (S2, on path, reproduced) `get_fsdata.R:678` get_FSdata
  - Reproducer: Appendix V091 (verdict demonstrated; wall clock in the appendix).
  - get_fsdata.R:666-667 generates names q1..qK; :678 and :700 assign df.FS[[names_new[i]]] with no test against names(df.FS).
  - Cuts are evaluated from the cache built at :618 before the overwrite, so the candidate columns themselves are correct; the damage is to the user's column in the returned frame (:746).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-052** (S2, on path, reproduced) `glm_effect_estimators.R:163` make_effect_estimator
  - Reproducer: Appendix V092 (verdict demonstrated; wall clock in the appendix).
  - Config: forestsearch(outcome_type = 'binary', effect_measure = 'IRR', offset.name = ..., adverse_outcome = FALSE); consistency_method = 'resample' is the default (forestsearch_main.R:1702).
  - forestsearch_main.R:2220-2229 passes adverse_outcome to make_effect_estimator(); glm_effect_estimators.R:162-172 routes IRR/IRD to .make_poisson_rate_estimator() without it, and that factory has no such formal (:726-735) and no flip in its closure (:744-915).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-053** (S2, on path, reproduced) `glm_effect_estimators.R:187` make_effect_estimator
  - Reproducer: Appendix V066, Appendix V093 (verdict demonstrated; wall clock in the appendix).
  - Merged duplicate reports: 1.
  - glm_effect_estimators.R:158-161 (binary) and :200-203 (count) run match.arg(); the continuous branch :187-195 never reads effect_measure and always builds the OLS estimator (measure 'MD', :975).
  - forestsearch() resolves a default only when the argument is NULL (forestsearch_main.R:1789-1795) and has no outcome-type cross-check (grep for a validation of effect_measure in R/ finds none outside the DGM and simulation helpers).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-054** (S2, on path, reproduced) `glm_effect_estimators.R:429` .make_glm_binary_estimator
  - Reproducer: Appendix V094 (verdict demonstrated; wall clock in the appendix).
  - Unimplemented pairs, unadjusted fit returned with an unadjusted method_used: RR at glm_effect_estimators.R:429-433 (.estimate_rr at :645 has no propensity argument); RD with dr_gcomp at :414-419 (.estimate_rd handles only iptw, :465); MD with dr_gcomp at :952-966 (only iptw, :956).
  - forestsearch_main.R:2639-2641 validates ps_adjust_method by match.arg() only and stores the requested value in args_call_all; :2705-2712 rebuilds the closure for any GLM measure. So forestsearch(outcome_type = 'binary', effect_measure = 'RR', is.RCT = FALSE, ps_adjust_method = 'iptw') estimates the score, attaches sw, and runs an unadjusted search.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-055** (S2, on path, reproduced) `glm_effect_estimators.R:535` .estimate_rd
  - Reproducer: Appendix V017 (verdict demonstrated; wall clock in the appendix).
  - Configuration: forestsearch(outcome_type = "binary", effect_measure = "RD"). Documented and supported but NOT the default: roxygen forestsearch_main.R:1049-1050 'For binary: "OR" (default), "RD", "RR", ...'; formal effect_measure = NULL (:1706) resolves to OR for binary (:1789-1790). RD was the binary default until the change announced in NEWS.md:552-561.
  - Call chain: forestsearch_main.R:2220-2228 make_effect_estimator(effect_measure = effect_measure, adjust_covariates = ...) -> glm_effect_estimators.R:174 .make_glm_binary_estimator -> closure :414-419 -> .estimate_rd. The fast-closure swap at forestsearch_main.R:2250 is OR only, so RD always uses this code. The existence guard at glm_effect_estimators.R:399 is OR/RR only.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-056** (S2, on path, reproduced) `glm_effect_estimators.R:671` .estimate_rr
  - Reproducer: Appendix V018 (verdict demonstrated; wall clock in the appendix).
  - Configuration: forestsearch() formals are outcome_type = "survival" (roxygen forestsearch_main.R:1047-1048) and effect_measure = NULL (:1706), resolved to "OR" for binary (:1789-1794). Roxygen :1049-1051 lists "RR" as a supported binary measure and NEWS.md:1457-1460 announces it, so RR is documented and supported but needs two non-default arguments (outcome_type = "binary", effect_measure = "RR"). No test or vignette runs forestsearch() with RR (the only tests/ hit, test-fs-dgm-feasibility.R:98, does not call the estimator); function_index gives .estimate_rr n_test_files = 0.
  - Reach: forestsearch_main.R:2220-2229 make_effect_estimator() -> glm_effect_estimators.R:174 .make_glm_binary_estimator -> RR branch :429-434 -> .estimate_rr. The fast closure replaces it only for continuous and binary OR (forestsearch_main.R:2241-2260), so RR always uses this code. With adjust_covariates the start vector (:655, length 2) mismatches the coefficient count, glm errors, :670 returns NULL and the fallback is taken, so only the unadjusted RR fit is exposed.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-057** (S2, on path, reproduced) `glm_effect_estimators.R:776` .make_poisson_rate_estimator
  - Reproducer: Appendix V095 (verdict demonstrated; wall clock in the appendix).
  - Config: outcome_type = "count" (or binary with IRR/IRD), effect_measure = "IRD", ps_method != "none", ps_adjust_method = "iptw".
  - forestsearch_main.R:2640 resolves ps_adjust_resolved; :2680 / glm_effect_estimators.R:1450 attach the `sw` column to df.analysis; :2705-2712 rebuild estimator_fn via make_effect_estimator(effect_measure = effect_measure, ps_adjust_method = ps_adjust_resolved).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-058** (S2, on path, reproduced) `grf_helpers.R:684` .grf_code_column
  - Reproducer: Appendix V096 (verdict demonstrated; wall clock in the appendix).
  - Trigger: a CHARACTER (not factor, not numeric) confounder on the GLM GRF path (grf.subg.harm.glm via use_grf or subgroup_method = "grf"). forestsearch() does no type validation of confounders before the GRF dispatch (forestsearch_main.R:1874-2892 touch df.analysis only for the zero-row check :1920, the outcome/threshold check :2181 and the PS columns), and tests/testthat/test-grf-membership-coding.R:59 treats a character column as a supported input.
  - Case A (numeric strings with NA, e.g. c("0","1",NA)): fit time -- grf_subg_harm_glm.R:373-380 drops incomplete rows, then :399 .build_grf_X -> grf_helpers.R:684 lvls = unique(x) has no NA -> :687 numeric 0/1, so the cut is on the 0/1 scale. Evaluation -- forestsearch_helpers.R:1987 calls .grf_evaluate_subgroup(sg_def, df) on the UNFILTERED frame -> grf_subgroup_labels.R:225 -> grf_helpers.R:684 unique(x) now contains NA -> :685 anyNA(as.numeric(lvls)) is TRUE -> :690 integer codes 1/2. A cut "z <= 0" then selects nobody and "z > 0" selects everybody with a non-missing z. A factor with the same values is unaffected (levels() excludes NA).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-059** (S2, on path, reproduced) `grf_main.R:226` grf.subg.harm.survival
  - Reproducer: Appendix V097 (verdict demonstrated; wall clock in the appendix).
  - Trigger: outcome_type = "survival" with use_grf = TRUE or subgroup_method = "grf" (both non-default: use_grf formal is FALSE at forestsearch_main.R:1639) and at least one factor/character confounder. forestsearch passes the raw frame (forestsearch_main.R:3078-3093; forestsearch_helpers.R:1858-1877).
  - grf_main.R:226: apply() on a data frame calls as.matrix(); with any factor/character column base R returns a character matrix in which numeric columns are format()-ed (7 significant digits by default) and logical columns become "TRUE"/"FALSE"; as.numeric then yields (a) rounded values for continuous columns and (b) NA, with a coercion warning, for factor/character columns with non-numeric labels and for logical columns. With all-numeric confounders none of this happens.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-060** (S2, on path, reproduced) `grf_main.R:235` grf.subg.harm.survival
  - Reproducer: Appendix V019, Appendix V098 (verdict demonstrated; wall clock in the appendix).
  - Merged duplicate reports: 1.
  - grf_helpers.R:632-658: validate_grf_data() warns and returns FALSE for length(unique(W)) < 2, length(unique(D)) < 2, or n_treat/n_control < n.min.
  - grf_main.R:234-236: if (!validate_grf_data(...)) return(create_null_result(data, NULL, list(), config)).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-061** (S2, on path, reproduced) `grf_main.R:295` grf.subg.harm.survival
  - Reproducer: Appendix V099 (verdict demonstrated; wall clock in the appendix).
  - Trigger: grf_selection = 'frontier', frontier_rule effMaxSG (the default rule) or effMinSG, dmin.grf < 0, every eligible candidate effect negative. dmin.grf is not range-checked anywhere (grep of forestsearch_main.R, grf_args.R, grf_main.R, grf_subg_harm_glm.R).
  - grf_subgroup_labels.R:361-362: elig is non-empty because effect >= dmin with dmin negative.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-062** (S2, on path, reproduced) `grf_main.R:297` grf.subg.harm.survival
  - Reproducer: Appendix V062, Appendix V100 (verdict demonstrated; wall clock in the appendix).
  - Merged duplicate reports: 1.
  - grf_main.R:291-303: .grf_frontier_select() returns NULL when no candidate has effect >= dmin.grf (grf_subgroup_labels.R:361-362); the return is create_null_result() plus sg_def and selection, with no $candidates. Same at grf_subg_harm_glm.R:528-545.
  - The success returns attach the unfiltered DR family: grf_main.R:269-276 and 331; grf_subg_harm_glm.R:499-506 and 573.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-063** (S2, on path, reproduced) `grf_subg_harm_glm.R:303` grf.subg.harm.glm
  - Reproducer: Appendix V101 (verdict demonstrated; wall clock in the appendix).
  - return_selected_cuts_only: forestsearch_main.R:3107 -> grf_args.R:64 -> formal at grf_subg_harm_glm.R:303. In the body (312-714) it appears only as an echoed element of the return lists (543, 566, 639, 702).
  - grf_subg_harm_glm.R:657-660: tree_cuts <- .extract_tree_cuts(best_tree, ...) with best_tree = trees[[best$depth]]; 927-940 reads that one tree only.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-064** (S2, on path, reproduced) `grf_subg_harm_glm.R:443` grf.subg.harm.glm
  - Reproducer: Appendix V020 (verdict demonstrated; wall clock in the appendix).
  - Configuration: forestsearch() formals forestsearch_main.R:1639 use_grf = FALSE, :1647 subgroup_method = c("consistency", ...), :1704 outcome_type = c("survival", ...): the default call never reaches grf.subg.harm.glm(). Reached only with outcome_type = 'count' AND (subgroup_method = 'grf' at :2892-2961, or use_grf = TRUE at :3071-3113), or by a direct call / the simulation GRF arm.
  - forestsearch_main.R:1797-1798 and 2198-2199: adverse_outcome resolves to TRUE for count; passed unchanged at :2951 -> forestsearch_helpers.R:1879-1901 -> grf_args.R:65, and at :3108 -> grf_args.R:65. No data-level flip of the count outcome anywhere upstream.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-065** (S2, on path, reproduced) `grf_subg_harm_glm.R:914` .select_best_sg_glm
  - Reproducer: Appendix V103 (verdict demonstrated; wall clock in the appendix).
  - grf_subg_harm_glm.R:616: best <- .select_best_sg_glm(values, sg.criterion, dmin.grf, n, n.min).
  - 904-912 ('mDiff'): returns NULL when max(scores) < dmin.grf.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-066** (S2, on path, reproduced) `subgroup_consistency_main.R:403` subgroup.consistency
  - Reproducer: Appendix V104 (verdict demonstrated; wall clock in the appendix).
  - Merged duplicate reports: 1.
  - subgroup_consistency_main.R:403 formal consistency_threshold = NULL; the body only forwards it (:776, :817, :837, :942, :954), never resolves it.
  - :427-436: a GLM caller (estimator_fn) without glm_resample_spec is switched to consistency_method = "split" with a warning.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-067** (S2, on path, reproduced) `subgroup_consistency_main.R:550` subgroup.consistency
  - Reproducer: Appendix V108 (verdict demonstrated; wall clock in the appendix).
  - subgroup_consistency_main.R:550-553 `if (identical(sg_focus, "maxeff")) found.hrs <- hr.subgroups`; the whole else block :554-574 is bypassed, and it holds both the HR >= hr.threshold floor and the m1 filter (:555-570).
  - m1.threshold is consumed nowhere else in the selection path (grep of R/: only forestsearch formal/doc at forestsearch_main.R:849/:1676, fs_family_report.R, and a read-back in fs_declaration_calibration.R:305).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **INV-5** (S2, on path, reproduced) `subgroup_consistency_main.R:717` subgroup.consistency
  - Reproducer: Appendix V109, Appendix V110b (verdict demonstrated; wall clock in the appendix).
  - Other reproducer runs: Appendix V110 (inconclusive (a forestsearch run failed; see FAILED lines)), Appendix R03_inv5_seq_vs_multisession (selection identical in both; Pcons differs between plans und).
  - Merged duplicate reports: 3.
  - WHAT IS PLAN-INDEPENDENT (holds by reading)
  - forestsearch_main.R:2059-2092: parallel_args is only validated in forestsearch(); its single consumer is subgroup.consistency(), reached through args_call_all (:1880-1881) and filter_call_args (:3629-3633; summary_utility_functions.R:679-688).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **INV-2** (S2, on path, reproduced) `subgroup_search.R:174` subgroup.search
  - Reproducer: Appendix V120, Appendix V044, Appendix V048, Appendix V138, Appendix V147, Appendix V149, Appendix V152, Appendix V079, Appendix V154, Appendix V156, Appendix V102, Appendix V159, Appendix V161, Appendix V165, Appendix V170, Appendix V176, Appendix V106, Appendix V112, Appendix V107, Appendix V190 (verdict demonstrated; wall clock in the appendix).
  - Other reproducer runs: Appendix V111 (NOT_REPRODUCED).
  - Merged duplicate reports: 25.
  - Site list from grep of all 125 files in R/ for plan(, options(, set.seed(, .Random.seed, RNGkind(, par(, handlers(, setDTthreads, setwd, Sys.setenv, sink, registerDo*, withr (the last six: no hits), each site then opened.
  - PLAN - NOT restored:
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-068** (S2, on path, reproduced) `subgroup_search.R:370` prepare_search_data
  - Reproducer: Appendix V113 (verdict demonstrated; wall clock in the appendix).
  - One column: subgroup_search.R:363 cbind(Y, Event, Treat, Z) is a 4-column matrix; :370 Z = temp[, -c(1, 2, 3)] drops to a vector; :132 L <- ncol(zz) is NULL; :160 generate_combination_indices(NULL, maxk) -> :389 combn(NULL, 1) stops with 'n < m'.
  - L < maxk: e.g. two columns and maxk = 3: :388-389 at k = 3 calls combn(2, 3) -> 'n < m' (no check of L against maxk anywhere in :98-160 or :379-407).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-069** (S2, on path, reproduced) `subgroup_search.R:628` evaluate_combination_with_status
  - Reproducer: Appendix V114 (verdict demonstrated; wall clock in the appendix).
  - subgroup_search.R:137-147: df_clean is built only when df_analysis is non-NULL; with df_analysis = NULL nothing is signalled; on a row mismatch :143-144 warns 'Using vector-based path' and sets df_clean <- NULL.
  - :628 the GLM branch needs !is.null(estimator_fn) && !is.null(df_clean); with df_clean NULL control falls to the survival path at :696 onward.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-070** (S2, on path, reproduced) `subgroup_search.R:771` has_positive_variance
  - Reproducer: Appendix V022 (verdict demonstrated; wall clock in the appendix).
  - Configuration: subgroup_search.R:15 roxygen reads 'Matrix or data frame of candidate subgroup factors (binary indicators)'; the only example (:87-95, inside \dontrun, so never run by R CMD check) builds Z <- df[, c("age", "meno", "size", "grade3", "nodes", "pgr", "er")], a data frame. Z has no default (formal at :98). So the type is documented, but it is reached only by a direct call to this exported lower-level engine.
  - subgroup_search.R:112-125: no type check or coercion of Z before prepare_search_data(Y, Event, Treat, Z).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-071** (S2, on path, reproduced) `summary_utility_functions.R:131` format_results
  - Reproducer: Appendix V115 (verdict demonstrated; wall clock in the appendix).
  - Default trigger: SG_tab_estimates(df, SG_flag = <flag column>) with the formal defaults sg1_name = NULL, sg0_name = NULL (summary_utility_functions.R:604-606).
  - :615-616 pass subgroup_name = NULL to analyze_subgroup -> format_results (:560-561); c(NULL, n, n_treat, d, m1, m0, drmst, hr) at :131 drops the name, giving 7 elements.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-072** (S2, on path, reproduced) `summary_utility_functions.R:366` SG_tab_estimates_glm
  - Reproducer: Appendix V116 (verdict demonstrated; wall clock in the appendix).
  - Trigger: binary or count outcome with a ratio measure (OR, RR, IRR) and est.scale = "1/hr"; the bootstrap keeps "1/hr" for these (bootstrap_dofuture_main.R:559-568) and passes it at :719; mr_estimates_table accepts it for ratio measures (mr_estimates_table.R:185-189) and passes it at :213.
  - summary_utility_functions.R:366 calls prepare_subgroup_data, which adds treat2 = 1 - treat and returns treat.name = "treat2" (:52-55); :367 treat_name_use = "treat2".
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-073** (S2, on path, reproduced) `summary_utility_functions.R:537` analyze_subgroup
  - Reproducer: Appendix V117 (verdict demonstrated; wall clock in the appendix).
  - Trigger: SG_tab_estimates on a flag whose 0- or 1-group has one arm only (reachable from the unwrapped CV calls forestsearch_cross_validation.R:1672-1699 when a cross-validated group is tiny).
  - analyze_subgroup guards only the empty or all-NA case (summary_utility_functions.R:518-522).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-074** (S2, on path, reproduced) `summary_utility_functions.R:629` SG_tab_estimates
  - Reproducer: Appendix V023 (verdict demonstrated; wall clock in the appendix).
  - Configuration: sg_tables() formals hr_1a = NA, hr_0a = NA (summary_utility_functions.R:764-765); roxygen :717-718 "Adjusted HR for subgroup 1 (optional)" / "subgroup 0 (optional)". Non-default; the only vignette use passes NA (vignettes/articles/forestsearch-2024.Rmd:423-424); no test passes a non-NA value.
  - Survival path: :1024-1034 filter_call_args() (:679-688) copies hr_1a/hr_0a into the argument list; :1036-1037 sets SG_flag = "ITT" and calls SG_tab_estimates with that list.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-075** (S2, on path, reproduced) `summary_utility_functions.R:629` SG_tab_estimates
  - Reproducer: Appendix V024 (verdict demonstrated; wall clock in the appendix).
  - Configuration: forestsearch() formal est.scale = "hr" (forestsearch_main.R:1637); its roxygen (:382) lists only "hr" or "rmst", but "1/hr" is the documented alternative of the downstream functions (SG_tab_estimates roxygen summary_utility_functions.R:586; bootstrap_dofuture_main.R:119; plot_sg_results.R:43, 140; plot_subgroup_results_forestplot.R:37, 127) and a fit made with est.scale = "1/hr" is exercised by tests/testthat/test-cv-no-subgroup-edges.R:62-67. Non-default; every vignette and the README use "hr".
  - forestsearch() only stores est.scale in args_call_all (no other executing use in forestsearch_main.R); sg_tables() has no est.scale argument and takes it from fs$args_call_all via filter_call_args (:1024-1034), then calls SG_tab_estimates twice with the same list (:1036-1040).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-076** (S2, on path, traced) `forestsearch_main.R:1920` forestsearch
  - Input: forestsearch(df.analysis = <matrix or other non-data.frame>, ...); n.min defaults to 60 (:1666) so the n.min block at :1848 is skipped, and no line between the formals and :1920 inspects df.analysis otherwise.
  - :1920 if (!is.data.frame(df.analysis) || nrow(df.analysis) == 0L) is TRUE for the type error.
  - :1921-1924 the only signal is warning('df.analysis has zero rows ...'), gated on !quiet.
  - :1925-1932 returns list(sg.harm = NULL, args_call_all, error_log = list(stage = "input_validation", message = 'df.analysis has zero rows; no subgroup search performed.')).
  - get_fsdata.R:99-106 holds the intended actionable stop() for a non-data.frame input, which forestsearch() can no longer reach.
  - Note: Description confirmed as written; severity raised from the reader's S3 to S2 because this is the brief's 'silent failure where an error is warranted' case (fully silent under quiet = TRUE, mislabelled otherwise), and the explicit stop() in get_FSdata shows an error is the intended behaviour. Invalid-input path only: tibbles and data.tables inherit from data.frame and are unaffected. Tests pin the three-field early-return shape for genuine failures (test-return-shape-contract.R:231) but none exercises a non-data.frame input.
- **A-077** (S2, on path, traced) `glm_effect_estimators.R:529` .estimate_rd
  - Input: a slice with n1 = 0 (or n0 = 0), effect_measure = 'RD', ps_adjust_method = 'none'. The existence check is scoped to OR/RR (:399), so RD falls through to .estimate_rd() (:414-419).
  - :461-462: the empty arm's proportion is NaN; :484 puts it in the start vector, glm() rejects the starting values, the error handler returns NULL (:504), tier 1 is skipped (:512).
  - :516-521: the logistic fit has the treatment column aliased (coefficient NA) and still converges; :527-528 predictions under treat = 0 and treat = 1 use only the intercept, raising a rank-deficient warning that the error-only tryCatch (:548) does not catch; :529 gives exactly 0.
  - :533: vcov() carries NA for the aliased coefficient, so se is NA; :538-546 return estimate 0, se NA, converged TRUE, method_used 'logistic_margins'.
  - Consumers test only is.na(estimate) (subgroup_search.R:1071; Cox_estimation_helpers.R:220), so the 0 is taken as an effect: in the search it is compared with the threshold (subgroup_search.R:681), and in the H/Hc fits it is reported with se NA.
  - Reach: the binary search requires d0.min/d1.min events per arm (subgroup_search.R:637-643; defaults 10, forestsearch_main.R:1687-1688), which rules it out at the defaults; it is reachable with a floor of 0, in H/Hc fits and through the exported closure.
  - Note: Raised from the reader's S3: returning a number for an estimand that does not exist is wrong behaviour in an edge configuration, which is S2 under the brief; it is not S1 because the default floors exclude it. The tier-2 steps rest on base-R handling of an aliased coefficient (predict warns and drops it; vcov returns NA for it), which the package's own comment at :1016-1018 states for the lm case; the reader's reproducer (rid glm#9, in the batch file) checks the returned list directly and should be run centrally. Evidence of intent: the comment at :235-238, NEWS (the estimability-boundary ...
- **INV-5a** (S2, on path, traced) `subgroup_consistency_main.R:897` subgroup.consistency
  - Not reachable through forestsearch(): stop_threshold is reset to NULL for hr/maxSG/minSG/hrMaxSG/hrMinSG (forestsearch_main.R:2116-2148) and for maxeff (:2027); subgroup.consistency() itself has no such reset (stop_threshold is only validated at subgroup_consistency_main.R:472-483). subgroup.consistency() is exported and has no other caller in R/.
  - Sequential plan: subgroup_consistency_main.R:717-723 -> loop :793-871; :853-869 breaks immediately after the first candidate m* with Pcons >= stop_threshold, so only candidates 1..m* can enter the result.
  - Parallel plan: :897-914; for any focus other than maxeffCons the batch size is min(n_workers, max(1, n_candidates %/% 4), n_candidates) (:913), i.e. 2 or more once workers >= 2 and n_candidates >= 8. The whole batch is evaluated (:990-1000) and every non-NULL row stored (:1005-1009) BEFORE the stop check (:1014-1038), so candidates after m* in the same batch enter the result.
  - Selection: sg_consistency_out() helpers:1045-1049 takes row 1 of sort_subgroups(); for "hr" the key is (-Pcons, -hr, K) (helpers:586-587), for maxSG/minSG Pcons is the tie-break (:590-597), for hrMaxSG/hrMinSG the band is computed over the rows present (:599-626). A later candidate in the batch with higher Pcons therefore wins on a parallel plan but was never evaluated on the sequential plan.
  - The batch size depends on n_workers (:913), so the same holds between workers = 1 and workers = 2 on one plan.
  - The source comment :899-907 acknowledges that the truncation is not sound for these foci; the roxygen at :156-161 nevertheless states that early stopping 'is valid' for "hr", "maxSG", "minSG".
  - Note: Distinct mechanism from INV-5: deterministic, no RNG, and present under the default consistency_method. It is outside what forestsearch() can pass, so it does not affect any forestsearch(), bootstrap or cross-validation analysis; it concerns direct callers of the exported function only. The documentation side (early stopping described as valid for hr/maxSG/minSG, batch size 1 documented for hr/minSG) is already cluster C422; the plan- and worker-dependent selection is the consequence that cluster does not state. No reproducer: a direct call needs the internal df.fs frame (Y/Event/Treat/id ...
- **A-078** (S2, on path, suspected) `subgroup_consistency_main.R:448` subgroup.consistency
  - Other reproducer runs: Appendix V105 (NOT_REPRODUCED).
  - Roxygen :75-76 names only Y, Event, Treat as required; validation :447-453 checks those three (and nothing when estimator_fn is supplied).
  - Selection proceeds; :1140-1158 sg_consistency_out() -> subgroup_consistency_helpers.R:1050 extract_subgroup().
  - helpers:934-937 `sg.harm.id <- rep(0, nrow(df)); if (nrow(df_temp) > 0 && "id" %in% names(df_temp)) sg.harm.id[df$id %in% df_temp$id] <- 1` -> without id the vector stays all 0.
  - helpers:939-942 `data.frame(id = df$id, treat.recommend = ifelse(sg.harm.id == 1, 0, 1))` -> df$id is NULL and is dropped; treat.recommend = 1 for everyone.
  - subgroup_consistency_main.R:1160-1164 and :1193-1198: returned as df_flag / sg.harm.id next to a non-NULL sg.harm.
  - Guard for the main entry: forestsearch() always builds df.fs with an id column (forestsearch_main.R:3390-3391), so only direct callers of the exported function are affected. No test calls subgroup.consistency() without id; no NEWS entry.
  - Note: Reproducer V105 did not reproduce the claim: without an id column the run ended in a caught error inside sg_consistency_out() (a warning), not in a silent all-zero result. The defect (no check for the id column) stands; its described symptom is not established. Direct-call path only (exported, tier 1 function, but forestsearch() always supplies id). S2 as a silent wrong membership vector where an input error is warranted; no conducted analysis is affected.
- **A-079** (S2, off path, reproduced) `ahr_analytic.R:267` .ahr_set_index
  - Reproducer: Appendix V118 (verdict demonstrated; wall clock in the appendix).
  - ahr_analytic.R:263-267: a character set returns `.ahr_cell(pop) %in% .ahr_set_cells[[set]]` directly; the emptiness stop at :272 is reached only by logical sets.
  - ahr_analytic.R:456-459: with b, h supplied and inside [0.2, 30] the cached route .ahr_fast() is taken.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-080** (S2, off path, reproduced) `ahr_analytic.R:622` ahr_design
  - Reproducer: Appendix V119 (verdict demonstrated; wall clock in the appendix).
  - ahr_analytic.R:610-618: only pop, n and sets are checked; :622 passes b, h (and b_R = b, h_R = h) to .ahr_exact(..., force_fast = TRUE).
  - ahr_analytic.R:455 -> :289: .ahr_psi() returns pop$psi as soon as either b or h is NULL, before the positivity checks at :290-295; ahr_limit() rejects the one-NULL case at :534, ahr_design() has no such line.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-081** (S2, off path, reproduced) `bootstrap_summaries_helpers.R:722` create_bootstrap_diagnostic_plots
  - Reproducer: Appendix V121 (verdict demonstrated; wall clock in the appendix).
  - forestsearch() stores est.scale unvalidated in args_call_all (formal at forestsearch_main.R:1637, capture at :1880-1881); the bootstrap reads it at bootstrap_dofuture_main.R:559.
  - get_dfRes(est.scale = "1/hr") -> ci_est(scale = "1/hr") returns H0 = exp(-Hobs) and H2 = exp(-q) (bootstrap_calculations_helpers.R:183-187, 286-312); these become boot_results$H_estimates / $Hc_estimates.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-082** (S2, off path, reproduced) `bootstrap_summaries_helpers.R:887` summarize_factor_presence_robust
  - Reproducer: Appendix V122 (verdict demonstrated; wall clock in the appendix).
  - DINA-mode replicate rows: Pcons stays NA_real_ (initialised at bootstrap_analysis_dofuture.R:787; the consistency block :797-828 does not run) while M.1.. are filled from sg.harm (:841-868) and any_found = 1 (:939-943).
  - summarize_bootstrap_subgroups.R:107-108 selects found rows by any_found (the comment at :101-106 says keying on Pcons mis-reports DINA runs) and passes sg_found to summarize_factor_presence_robust() at :469-474.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-083** (S2, off path, reproduced) `calibrate_null_correction.R:357` run_null_calibration
  - Reproducer: Appendix V123 (verdict demonstrated; wall clock in the appendix).
  - calibrate_null_correction.R:346-354: fs_call_args = c(list(df.analysis), fs_args, list(hr.threshold = c1, hr.consistency = c2, seedit = seed_s)); if fs_args also carries one of these names (roxygen :208-211 only says they are 'excluded'), or any unknown name (forestsearch() has no '...', forestsearch_main.R:1622-1712), do.call(forestsearch, ...) errors on every replicate.
  - :357-360: tryCatch(do.call(forestsearch, fs_call_args), error = function(e) NULL); no counter, message or warning.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-084** (S2, off path, reproduced) `compare_fs_subgroups.R:123` compare_fs_subgroups
  - Reproducer: Appendix V124 (verdict demonstrated; wall clock in the appendix).
  - grf_main.R:235, :297, :363: no subgroup -> create_null_result(data, ...); grf_helpers.R:540-543 returns list(data = data, sg.harm.id = NULL, ...) with the caller's data unchanged; treat.recommend is only added on the success path (grf_main.R:309).
  - compare_fs_subgroups.R:102 not numeric/logical; :107 not class 'forestsearch'; :123-124 requires 'treat.recommend' %in% names(obj$data) -> FALSE.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-085** (S2, off path, reproduced) `compare_selection_rules.R:582` compare_selection_rules
  - Reproducer: Appendix V125 (verdict demonstrated; wall clock in the appendix).
  - compare_selection_rules.R:533-538: six lists pre-allocated with vector('list', n_combos) and named by combo_labels.
  - :579-585: fs_k is NULL when forestsearch() errored, err_k is NULL when it succeeded; `fs_list[[k]] <- fs_k` and `error_list[[k]] <- err_k`. In R, `x[[k]] <- NULL` on a list removes element k (and is a no-op when k > length(x)).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-086** (S2, off path, reproduced) `cox_ahr_cde_wrapper.R:692` cox_ahr_cde_analysis
  - Reproducer: Appendix V126 (verdict demonstrated; wall clock in the appendix).
  - cox_ahr_cde_wrapper.R:193-196: if every subject has loghr_po < log(hr_threshold), optimal_cutpoint = min(z) and (:206-213) question_idx is all FALSE; :197-203: if no subject qualifies, the fallback sets optimal_cutpoint = max(z) with a warning, so recommend_idx holds one subject when the maximum is unique.
  - :224-225 and :270-275: the empty-side AHR/CDE are NaN (mean of numeric(0)), reported in subgroup_stats without a message.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-087** (S2, off path, reproduced) `cox_ahr_cde_wrapper.R:849` summary.cox_ahr_cde
  - Reproducer: Appendix V025 (verdict demonstrated; wall clock in the appendix).
  - Configuration: unconditional. cox_ahr_cde_wrapper.R:116 formal `spline_df = 3`; no branch on spline_df, hr_threshold, plot_select, plot_style or the theta columns changes what is stored.
  - cox_ahr_cde_wrapper.R:163-173 cox_ahr_cde_analysis() always calls cox_cs_fit(); :723-724 stores the result unchanged as `cox_fit`; :736 sets class c('cox_ahr_cde','list'). Single return at :761.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-088** (S2, off path, reproduced) `create_summary_table.R:116` create_summary_table
  - Reproducer: Appendix V127 (verdict demonstrated; wall clock in the appendix).
  - Trigger: data[[treat_var]] contains NA (plain data.frame or tibble input).
  - create_summary_table.R:107 sort(unique(.)) drops NA, so the 2-level check at :110 passes.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-089** (S2, off path, reproduced) `create_summary_table.R:116` create_summary_table
  - Reproducer: Appendix V128 (verdict demonstrated; wall clock in the appendix).
  - create_summary_table.R:15 documents only 'must have 2 levels'; no 0/1 requirement is stated or validated (:110 checks the count of levels only).
  - :107 treat_levels <- sort(unique(data[[treat_var]])); :116 data_0 = first sorted level, :117 data_1 = second.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-090** (S2, off path, reproduced) `create_summary_table.R:747` create_summary_table
  - Reproducer: Appendix V129 (verdict demonstrated; wall clock in the appendix).
  - Trigger: create_summary_table(data, treat_var, vars_continuous = "age") on complete data (same for a single vars_binary, or a single 0/1-coded vars_categorical).
  - create_summary_table.R:193-203 adds exactly one row per continuous variable; :206 adds no Missing row when miss_0 = miss_1 = 0.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-091** (S2, off path, reproduced) `cv_summary_tables.R:601` cv_metrics_tables
  - Reproducer: Appendix V130 (verdict demonstrated; wall clock in the appendix).
  - Trigger: cv_metrics_tables(cv, use_gt = FALSE, include_raw = TRUE) (table_style default 'combined', :347/:355).
  - cv_summary_tables.R:447-456 result <- combined_df, a 12-row data.frame (4 agreement + 8 finding rows).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-092** (S2, off path, reproduced) `dina_bagged.R:194` dina_fit_bagged
  - Reproducer: Appendix V131 (verdict demonstrated; wall clock in the appendix).
  - Merged duplicate reports: 1.
  - dina_bagged.R:194 passes seed only to .dina_make_bootstrap_indices(); :423-437 set.seed(seed), draw the indices, and restore the previous .Random.seed on exit, so the stream is back at the caller's state before any bag is fitted.
  - Sequential branch :222-237 -> .dina_safe_one_bag :383 -> .dina_fit_one_bag :348-357 calls dina_fit(..., seed = NULL) for every bag.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-093** (S2, off path, reproduced) `dina_bagged.R:365` .dina_subset_Y
  - Reproducer: Appendix V132 (verdict demonstrated; wall clock in the appendix).
  - dina_bagged.R:43 `@inheritParams dina_fit` -> dina.R:63-66 documents Y for cox as 'a data frame (or matrix coercible to one)'; dina.R:183 `if (is.matrix(Y)) Y <- as.data.frame(Y)` makes that work in dina_fit().
  - dina_bagged.R:149-158: dina_fit_bagged() coerces X only; Y is passed untouched to .dina_safe_one_bag (dina_bagged.R:209-219 parallel, 226-236 sequential).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-094** (S2, off path, reproduced) `dina_subgroup_bootstrap.R:643` dina_subgroup_bootstrap
  - Reproducer: Appendix V133 (verdict demonstrated; wall clock in the appendix).
  - dina_subgroup.R:435: tau_hat <- tau_sign * (beta0 + X beta); with tau_sign = -1 the floor, ranking, mean_tau_hat and ci (:596-598, :610-612) are on the negated scale; tau_sign is kept only in sg$call (:372-373), not as a list element.
  - dina_subgroup_bootstrap.R:412-459: the sg validation never inspects sg$call$tau_sign; sg is accepted (sg_point <- sg at :460).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-095** (S2, off path, reproduced) `dina_subgroup_bootstrap.R:805` .dina_safe_one_bootstrap
  - Reproducer: Appendix V134 (verdict demonstrated; wall clock in the appendix).
  - dina_subgroup_bootstrap.R:273-296: formals hold m_diff, n_min, direction, sg_focus, effect_neighborhood, alpha; grep for max_depth|grid_probs|selection_rule|tau_sign in the file returns nothing; `...` goes to dina() only (:391, :799-803).
  - dina_subgroup_bootstrap.R:412-459: a supplied sg is checked for class, n_total, mask length, m_diff, n_min, alpha, family, sg_focus and (band foci only) effect_neighborhood; sg$max_depth and sg$selection_rule (stored at dina_subgroup.R:618, 624) are never read. The comment at :443-445 states the intent that point and per-iteration search use the same rule.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-096** (S2, off path, reproduced) `dina_subgroup_refit.R:435` .dina_signature_string
  - Reproducer: Appendix V026 (verdict demonstrated; wall clock in the appendix).
  - CONFIGURATION
  - dina_subgroup.R:323 formal `max_depth = 2L`; roxygen :172-177 '`1L` or `2L` (default `2L` ...)'; return shape documented :246-248 'length 2 (one entry per covariate) when `max_depth = 2L` selects a conjunction'. :458-476 appends pair candidates; :565-572, :606-608 return length-2 covariate/direction/threshold when a pair ranks first.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-097** (S2, off path, reproduced) `dina_subgroup_refit.R:498` print.dina_subgroup_refit
  - Reproducer: Appendix V135 (verdict demonstrated; wall clock in the appendix).
  - dina_subgroup_refit.R:156, 179-181: alpha is a user argument in (0, 1).
  - dina_subgroup_refit.R:207-209: ci = effect -/+ qnorm(1 - alpha/2) * se; alpha stored at :232.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-098** (S2, off path, reproduced) `explain_pareto_selection.R:240` .epx_classify_loss
  - Reproducer: Appendix V136 (verdict demonstrated; wall clock in the appendix).
  - Actual keys, subgroup_consistency_helpers.R:571 (maxeff: -hr, K), :583 (maxeffCons: -hr, K), :587 (hr: -Pcons, -hr, K), :591 (maxSG: -N, -Pcons, K), :595 (minSG: N, -Pcons, K), :621-624 (hrMaxSG/hrMinSG: -in_band, -/+N, -Pcons, -hr, K).
  - explain_pareto_selection.R:240-245: N_loses switch has branches only for maxSG/hrMaxSG/minSG/hrMinSG (default FALSE) - correct for hr/maxeff.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-099** (S2, off path, reproduced) `explain_pareto_selection.R:341` .epx_selected
  - Reproducer: Appendix V137 (verdict demonstrated; wall clock in the appendix).
  - subgroup_consistency_helpers.R:1078-1079: pareto_frontier is computed for every sg_focus from the full sorted result table; sort keys at :571-595 (hr = (-Pcons, -hr, K), minSG = (N, -Pcons, K)) make a dominated top row ordinary for the default sg_focus = "hr" (forestsearch_main.R:1672) and for minSG / hrMinSG.
  - explain_pareto_selection.R:107, 112: selected_m = result[1]$m; selected_on_frontier <- selected_m %in% ft$m is FALSE.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-100** (S2, off path, reproduced) `forestsearch_cross_validation.R:501` forestsearch_Kfold
  - Reproducer: Appendix V139 (verdict demonstrated; wall clock in the appendix).
  - Policy: forestsearch_cross_validation.R:15-23 and :481-484 say any refit error is fatal and only a clean no-subgroup fold falls back to ITT.
  - forestsearch_main.R:3310-3325 tryCatch around get_FSdata turns an error into warning(); :3327-3332 returns list(sg.harm = NULL, args_call_all, error_log). Same shape at :3495-3506 (subgroup.search). :3654-3683 turns a subgroup.consistency error into a warning and continues with grp.consistency = NULL (sg.harm stays NULL, no error_log).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-101** (S2, off path, reproduced) `forestsearch_cross_validation.R:1417` forestsearch_KfoldOut
  - Reproducer: Appendix V027 (verdict demonstrated; wall clock in the appendix).
  - Configuration: default arguments. forestsearch_tenfold(fs.est, sims, Kfolds = 10, ...) (formals :803-811) and forestsearch_Kfold(fs.est, Kfolds = .fs_cv_nobs(fs.est), ...) (formals :287-296) are exported and documented, and both appear in the case-study article vignettes/articles/forestsearch-2024.Rmd:606, 644, 660. The trigger is a DATA outcome, not a setting: every held-out row of one K-fold run (or of one tenfold simulation) has treat.recommend = 1.
  - Coding: get_dfpred() sets treat.recommend 0 = in H, 1 = complement (forestsearch_helpers.R:117, :137); treat.recommend.original is df.est's flag renamed (:133-134), so it has both levels whenever fs.est$sg.harm is non-NULL.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-102** (S2, off path, reproduced) `forestsearch_cross_validation.R:1672` forestsearch_KfoldOut
  - Reproducer: Appendix V140 (verdict demonstrated; wall clock in the appendix).
  - Trigger: fs.est fitted with est.scale = '1/hr' (survival), forestsearch_Kfold(), then forestsearch_KfoldOut(res, outall = TRUE).
  - forestsearch_cross_validation.R:574-579: resCV$treat2 <- 1 - treat; local treat.name <- 'treat2'; sg0.name <- 'Recommend'; sg1.name <- 'Not recommend'. Lines :580-648 never read treat.name again; cv_args$treat.name keeps the original name; grep finds no reader of resCV$treat2 in R/.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-103** (S2, off path, reproduced) `forestsearch_cross_validation.R:1798` CV_sgs
  - Reproducer: Appendix V141 (verdict demonstrated; wall clock in the appendix).
  - Trigger: original sg.harm of length 1 (sg_depth == 1, forestsearch_cross_validation.R:1794) under maxk = 2 (the default), and a fold whose selected subgroup is {A} & {B} with {A} the original factor.
  - :1798 exact_match <- ifelse((sg1 == sg1a | sg2 == sg1a), 1, 0) gives 1 for (sg1 = A, sg2 = B) and for (sg1 = B, sg2 = A); nothing requires the fold's other factor to be absent. (sg1 = A, sg2 = NA gives TRUE | NA = TRUE, the genuine exact case.)
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-104** (S2, off path, reproduced) `forestsearch_cross_validation.R:1887` find_covariate_any_match
  - Reproducer: Appendix V142 (verdict demonstrated; wall clock in the appendix).
  - Called from CV_sgs (forestsearch_cross_validation.R:1791-1792, :1806) for the Cov1 / Cov2 / 'Cov 1 & 2' find metrics (:1400-1402).
  - (a) Substring match, single-candidate path: confs = c('age','stage','sex'), sg_target = '{age <= 50}'. :1868 charmatch of '{age}' fails, :1873-1875 only '{age' is a prefix, :1894 cfs = 'age'. :1902-1903 grepl('age', '{stage <= 2}') is TRUE, so cov_any = 1 for a fold that used stage, not age. Same for x1 against x10, x12.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-105** (S2, off path, reproduced) `format_fs_details.R:147` format_fs_details
  - Reproducer: Appendix V143 (verdict demonstrated; wall clock in the appendix).
  - Trigger: a non-table line longer than max_width (default 48) containing an unbroken token (no space, no comma) longer than max_width minus the continuation indent, e.g. fs_output = c("Candidate factors: 1", paste0("  ", strrep("x", 60))).
  - format_fs_details.R:199-215: the line is not blank, not a separator (:204 needs a leading run of = or -), not columnar (:213), so wrap_line(ln, max_width) is called (:220).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-106** (S2, off path, reproduced) `fpr_approximation.R:112` fpr_approx
  - Reproducer: Appendix V144 (verdict demonstrated; wall clock in the appendix).
  - fpr_approximation.R:44-46 roxygen documents effect_type = "difference" (RD, MD, identity-scale SE); :99 it is a formal choice.
  - :111-119 stopifnot(..., effect_true > 0 (:112), ..., c1 > 0 (:115), ...) runs unconditionally, right after match.arg at :107.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-107** (S2, off path, reproduced) `fpr_calibration.R:236` fpr_calibration
  - Reproducer: Appendix V028 (verdict demonstrated; wall clock in the appendix).
  - Configuration is documented and uses default formals: fpr_calibration.R:199 outcome_type_approx = NULL; roxygen :69-71 requires outcome_type in fs_params only 'for binary/continuous/count'; roxygen :89-92 says the NULL default is 'inferred from fs_params$outcome_type' and states no binary fallback.
  - fpr_calibration.R:236: otype <- outcome_type_approx %||% fs_params$outcome_type %||% "binary" (%||% at :594 returns the right side only on NULL).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-108** (S2, off path, reproduced) `fpr_calibration.R:253` fpr_calibration
  - Reproducer: Appendix V145 (verdict demonstrated; wall clock in the appendix).
  - Trigger: fs_params with effect_measure = "RD" / "MD"; the file header (:24-25) lists RD and MD as supported, :244 has a continuous d_eff branch, and the roxygen (:49-58) says c1/c2 are on the natural scale.
  - fpr_calibration.R:210-211 validation requires c1 > 0 and c2 > 0; forestsearch_main.R:2284-2285 (RD/IRD) and :2336-2337 (MD) resolve the default identity-scale consistency threshold to 0.0, which this wrapper therefore cannot calibrate ('c2 must be a positive scalar').
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-109** (S2, off path, reproduced) `fpr_calibration.R:259` fpr_calibration
  - Reproducer: Appendix V146 (verdict demonstrated; wall clock in the appendix).
  - fpr_calibration.R:66-68 roxygen: fs_params holds 'All forestsearch() arguments except hr.threshold and hr.consistency ... and df.analysis'; effect.threshold, consistency.threshold and seedit are not excluded.
  - :259-265 modifyList(fs_params, list(hr.threshold = c1, hr.consistency = c2, ...)) leaves fs_params$effect.threshold / $consistency.threshold in place.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-110** (S2, off path, reproduced) `fpr_calibration.R:543` fpr_calibration_transfer
  - Reproducer: Appendix V148 (verdict demonstrated; wall clock in the appendix).
  - fpr_calibration.R:383-390 the fpr_calibration object stores call_args = nsims, nsims_verify, n_min, otype, n_workers, sim_seed only; no data, no fs_params, no null_fn.
  - :545-546 df.analysis = cal_result$call_args$df.analysis %||% eval(cal_result$call_args$null_fn_expr): both fields are absent, so NULL; :547 fs_params = list() is hard-coded; :553 list(...) is appended after them.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-111** (S2, off path, reproduced) `fs_fdr_report.R:377` .fs_fdr_one_rep
  - Reproducer: Appendix V029 (verdict demonstrated; wall clock in the appendix).
  - Configuration: fs_fdr_report(fit, df) with the default confounders = NULL (formal at fs_fdr_report.R:187) on a survival fit whose covariates are not GBSG columns. Two triggers traced: (a) confounders.name shares no name with the simulated frame; (b) partial overlap, and the inherited conf.cont_jcuts names a covariate that the intersection drops.
  - fs_fdr_report.R:223-230: conf_use = inherited confounders.name whenever it is non-empty; the GBSG fallback (:229) is taken only when the fit carried none, before any intersection.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-112** (S2, off path, reproduced) `fs_identification_summary.R:321` plot_fs_identification_involvement
  - Reproducer: Appendix V153 (verdict demonstrated; wall clock in the appendix).
  - fs_identification_summary.R:313 roles defaults to c(anchor, partner, proxy); :314-316 labs_role is always the three labels in anchor/partner/proxy order, whatever `roles` is.
  - :320 keeps the involvement rows whose covariate is in `roles`; :321-322 factor(d$covariate, levels = roles, labels = labs_role[seq_along(roles)]) pairs levels and labels by position.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-113** (S2, off path, reproduced) `generate_aft_dgm_helpers.R:698` process_factor_subgroup
  - Reproducer: Appendix V155 (verdict demonstrated; wall clock in the appendix).
  - generate_aft_dgm_helpers.R:613-622: define_subgroups() sends a variable that is neither in continuous_vars nor numeric to process_factor_subgroup().
  - :687-690 handles only a length-1 numeric/character/factor value; :692-697 handles only list(type = "multiple").
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-114** (S2, off path, reproduced) `generate_aft_dgm_main.R:226` generate_aft_dgm_flex
  - Reproducer: Appendix V030 (verdict demonstrated; wall clock in the appendix).
  - Configuration: generate_aft_dgm_main.R:226 formal 'treatment_var = NULL'; roxygen :19-20 'If NULL, treatment will be randomly simulated with 50/50 allocation' (same text in man/generate_aft_dgm_flex.Rd:57-58). So the trigger is literally the documented default of an exported function.
  - Reach: treatment_var is used in the body only at :296 and :302 (git grep treatment_var -- R/generate_aft_dgm_main.R); there is no is.null(treatment_var) branch anywhere in R/.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-115** (S2, off path, reproduced) `gg_forest.R:226` gg_forest
  - Reproducer: Appendix V031 (verdict demonstrated; wall clock in the appendix).
  - Configuration: non-default. Formal is cat_colours = NULL (gg_forest.R:121); roxygen gg_forest.R:66-67 says only 'Optional named character vector mapping category labels to colours. Defaults to grey for all rows.' No ITT requirement is stated anywhere (man example is an empty \dontrun, gg_forest.R:105-109).
  - gg_forest.R:146-155: validation checks lengths only; names(cat_colours) is never inspected.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-116** (S2, off path, reproduced) `glm_effect_profile.R:377` .gep_resolve_family
  - Reproducer: Appendix V157 (verdict demonstrated; wall clock in the appendix).
  - Trigger: glm_effect_profile(..., effect_measure = "log_RR", overdispersion = "quasi" or "negbin").
  - glm_effect_profile.R:137-144: the applicability check accepts log_RR, so overdispersion keeps its value with no warning.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-117** (S2, off path, reproduced) `glm_spline_fit.R:320` glm_cs_fit
  - Reproducer: Appendix V032 (verdict demonstrated; wall clock in the appendix).
  - Configuration: glm_spline_fit.R:137 formal is family = "binomial"; roxygen :36-42 documents 'One of "binomial" (default), "poisson", or "gaussian"' and :100-101 gives the gaussian model, so gaussian is supported and documented but NOT the default.
  - glm_spline_fit.R:179-188: family_name "gaussian" passes the check and maps to stats::gaussian().
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-118** (S2, off path, reproduced) `grf_args.R:189` extract_grf_args
  - Reproducer: Appendix V158 (verdict demonstrated; wall clock in the appendix).
  - grf_args.R:189-203 (survival) and :206-223 (GLM): the builders are called without grf_selection / frontier_rule / effect_neighborhood / selection_rule, so the builder defaults apply ("tree", "effMaxSG", 0.10, "neighborhood"; grf_args.R:41-44, :99-102), and return_selected_cuts_only = a$return_selected_cuts_only.
  - The actual subgroup_method = "grf" call: forestsearch_main.R:2893 resolves grf_selection (formal default c("frontier", "tree") at :1650, i.e. "frontier"); :2935-2944 maps sg_focus to frontier_rule; :2945-2965 passes them with effect_neighborhood and selection_rule to .forestsearch_grf_select(); forestsearch_helpers.R:1858-1901 builds the args with return_selected_cuts_only = TRUE (forced), grf_selection, frontier_rule, effect_neighborhood, selection_rule and data = df.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-119** (S2, off path, reproduced) `guohe_algorithm3.R:79` .g3_coef
  - Reproducer: Appendix V033 (verdict demonstrated; wall clock in the appendix).
  - Configuration: formal fast = NULL (guohe_algorithm3.R:248) resolves to the closed form for unadjusted binary/continuous (:277-279). The glm path (:355-365 -> .g3_coef) is reached only with an explicit fast = FALSE or a non-NULL adjust_covariates; both are documented arguments (:183-189; :165-174, the latter stated as untested). No caller inside R/ passes fast = FALSE: guohe_adaptive_r() forwards its own fast (default NULL, guohe_adaptive_r.R:194, :257, :287) and fs_to_guohe() is survival only (fs_to_guohe.R:357-362).
  - Closed form: guohe_algorithm3.R:310-311 bad when an arm < min_events OR any of a, b, cc, dd <= 0; :317 NA; full sample :370-377 dropped with a warning; resamples :400-401 excluded and counted in drop_rate (:402, :424).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-120** (S2, off path, reproduced) `guohe_algorithm3.R:112` .g3_cox_lean
  - Reproducer: Appendix V034 (verdict demonstrated; wall clock in the appendix).
  - Configuration is the default one: outcome = "survival" is the first choice (guohe_algorithm3.R:233), min_events = 5L (:244), fast = NULL -> FALSE for survival (:277-279), lean_cox TRUE (:326-327). The exported bridge fs_to_guohe() calls exactly this with orient = +1 and no adjustment (fs_to_guohe.R:357-362; min_events default 5L at :341).
  - Guards: guohe_algorithm3.R:110 two arms, :111 per-arm SUBJECT count, :112 TOTAL events >= min_events (as documented at :162-163). No per-arm event requirement.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-121** (S2, off path, reproduced) `mrct_case_study.R:282` mrct_case_study
  - Reproducer: Appendix V160 (verdict demonstrated; wall clock in the appendix).
  - mrct_simulation.R:405 .options.future = list(seed = TRUE): each replicate's future gets an L'Ecuyer-CMRG .Random.seed, so the active kind is L'Ecuyer for simulate_from_dgm() and for forestsearch() (do.call at :572). The package states this mechanism itself at run_subgroup_sims.R:869-879.
  - mrct_case_study.R:280-282: RNGkind("L'Ecuyer-CMRG") only around simulate_from_dgm(); :282 do.call(RNGkind, as.list(old_kind)) puts the caller's kind (Mersenne-Twister by default) back before .mrct_analyze_trial() (:301) calls forestsearch() (mrct_trial_helpers.R:47).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-122** (S2, off path, reproduced) `mrct_inject_dgm.R:414` inject_mrct_structure
  - Reproducer: Appendix V162 (verdict demonstrated; wall clock in the appendix).
  - Input: inject_mrct_structure(..., region = list(..., loghr = NULL), x3 = list(vars, cuts, loghr = L)) without region$observed. The roxygen describes loghr = NULL for the observed mode only (:304-314).
  - mrct_inject_dgm.R:414-416: utils::modifyList(list(..., loghr = 0, ...), region) removes loghr (a NULL value deletes the element); the validation at :409-412 is inside the observed branch only.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-123** (S2, off path, reproduced) `mrct_inject_dgm.R:558` inject_mrct_structure
  - Reproducer: Appendix V163 (verdict demonstrated; wall clock in the appendix).
  - Input: inject_mrct_structure(..., region_loghr_shift != 0 or region_harm_loghr_shift != 0), then simulate_from_dgm(dgm, baseline = "fixed").
  - mrct_inject_dgm.R:506-528: generate_aft_dgm_flex() on the seed rows; generate_aft_dgm_main.R:388-402 builds df_source with lin_pred_1 / theta_1 / loghr_po from the fitted coefficients and :424 stores it in the DGM.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-124** (S2, off path, reproduced) `mrct_inject_dgm.R:635` analysis_covariates
  - Reproducer: Appendix V164 (verdict demonstrated; wall clock in the appendix).
  - Input: a seed trial with a covariate such as prior_treat or risk_k in continuous_vars / factor_vars of inject_mrct_structure() (prior_treat is the name the package's own MRCT example data uses, mrct_simulation.R:933, 1187).
  - generate_aft_dgm_helpers.R (process_continuous_vars / process_factor_vars, marker = "z_") names its column z_prior_treat in df_super; the only helper columns are z_<x_pred>_treat, z_<x_pred>_k, z_<x_pred>_k_treat (generate_aft_dgm_helpers_outcome.R:241-254).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-125** (S2, off path, reproduced) `mrct_simulation.R:319` mrct_region_sims
  - Reproducer: Appendix V035 (verdict demonstrated; wall clock in the appendix).
  - CONFIGURATION
  - mrct_simulation.R:270 formal `confounders.name = NULL`; roxygen :51-52 'If NULL, automatically extracted from dgm'; the function's own example (:232-239, \dontrun) omits confounders.name on a create_dgm_for_mrct() DGM. So the trigger is the formal default plus the package's own DGM builder.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-126** (S2, off path, reproduced) `mrct_simulation.R:422` mrct_region_sims
  - Reproducer: Appendix V166 (verdict demonstrated; wall clock in the appendix).
  - Trigger: sim_args = list(..., entry_var = "enroll") and dgm$df_super without a column 'entrytime'.
  - mrct_simulation.R:419-426: modifyList(sim_args, list(dgm=, n=, entry_var = if ("entrytime" %in% names(dgm$df_super)) "entrytime" else NULL, ...)); list() keeps the NULL element and utils::modifyList (keep.null = FALSE) assigns x[["entry_var"]] <- NULL, which deletes the user's entry.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-127** (S2, off path, reproduced) `mrct_simulation.R:570` mrct_region_sims
  - Reproducer: Appendix V167 (verdict demonstrated; wall clock in the appendix).
  - Trigger (corrected): any input for which forestsearch() itself stops, e.g. a misspelled fs_args name (forestsearch has no '...' formal, forestsearch_main.R:1626-1720, so do.call raises 'unused argument'), or use_twostage not a single logical (forestsearch_main.R:2154-2156), or the other stop() calls in the body (e.g. :1779, :1854, :2617).
  - mrct_simulation.R:569-584: tryCatch(do.call(forestsearch, fs_call_args), error = function(e) { if (verbose_fs) cat(...); NULL }); verbose_fs = details && sim <= verbose_n_sims (line 503), FALSE by default, so nothing is printed.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-128** (S2, off path, reproduced) `mrct_trial_helpers.R:60` .mrct_analyze_trial
  - Reproducer: Appendix V168 (verdict demonstrated; wall clock in the appendix).
  - Trigger: an observed trial in which some region patients have NA in a covariate that the declared rule uses (the covariates named in fs_args$confounders.name).
  - mrct_trial_analysis.R:93 -> mrct_four_sets.R:120-128: .mrct_check_trial() tests anyNA only on the named analysis columns; confounders are not examined.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-129** (S2, off path, reproduced) `pareto_frontier_table.R:288` pareto_frontier_table
  - Reproducer: Appendix V169 (verdict demonstrated; wall clock in the appendix).
  - Trigger: ci <- compute_frontier_cis(fs, ci_level = 0.90); pareto_frontier_table(fs, ci_table = ci) or explain_pareto_selection(..., ci_table = ci).
  - frontier_cis.R:161: z_crit <- qnorm(0.5 + ci_level / 2) drives all three intervals (:256-257, :299-304).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-130** (S2, off path, reproduced) `plot_km_band_forestsearch.R:351` plot_km_band_forestsearch
  - Reproducer: Appendix V171 (verdict demonstrated; wall clock in the appendix).
  - plot_km_band_forestsearch.R:349-367 passes only sg_labels; sg_cols is used for existence checks (:284) and messages (:329-331).
  - weightedsurv 0.1.0 R/kmplotting_helpers.R:499: sg_flags <- lapply(sg_labels, function(expr) with(df, eval(parse(text = expr)))); :508 df_sg <- df[sg_flags[[i]], ].
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-131** (S2, off path, reproduced) `plot_km_band_forestsearch.R:487` generate_readable_sg_labels
  - Reproducer: Appendix V172 (verdict demonstrated; wall clock in the appendix).
  - Trigger: plot_km_band_forestsearch(df, fs.est = fs) with fs from subgroup_method = "dina" (max_depth 2 default) or "grf" and a 2-cut subgroup; sg.harm = c("{x1 <= 3.5}", "{x2 >= 10}") (forestsearch_helpers.R:1660-1662) or the GRF label vector / union string (forestsearch_helpers.R:1983). out_sg is the selector object (forestsearch_helpers.R:1672-1678, 1998-2004); sg.harm_label is produced only in subgroup_consistency_helpers.R:926/946, so Method 1 at plot_km_band_forestsearch.R:458 is skipped; subgroup_definition (Method 2, :475) is set by no fit.
  - plot_km_band_forestsearch.R:484-485 parse_sg_harm_to_expression(): each element matches [<>=!] and is kept (:592-594); :613-617 returns ONE string "{x1 <= 3.5} & {x2 >= 10}".
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-132** (S2, off path, reproduced) `plot_sg_results.R:214` plot_sg_results
  - Reproducer: Appendix V173 (verdict demonstrated; wall clock in the appendix).
  - Trigger: plot_sg_results(fs, est.scale = "1/hr") (roxygen :42-43 documents it as the inverse effect scale; also reachable through plot.forestsearch(...) dots, forestsearch_methods.R:621-628).
  - plot_sg_results.R:214-219: df$treat_plot <- 1 - df[[treat.name]]; treat_labels <- rev(treat_labels); sg0_name/sg1_name swapped. grep over R/ and tests/: treat_plot occurs only at :215 and :221.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-133** (S2, off path, reproduced) `plot_sg_results.R:504` compute_sg_summary
  - Reproducer: Appendix V036 (verdict demonstrated; wall clock in the appendix).
  - Configuration: event indicator column numeric 0/1 with storage type double. Supported: get_fsdata.R:177-180 only requires is.numeric() ('must be numeric (0/1)'); roxygen forestsearch_main.R:345 and plot_sg_results.R:29 say only 'event indicator (1=event, 0=censored)'. No formal or doc requires integer storage.
  - Pipeline keeps the type: forestsearch_main.R:2725 df <- df.analysis; :3310-3318 get_FSdata(df.analysis = df.analysis); get_fsdata.R:107 df.FS <- df.analysis and :746 return(list(df = df.FS, ...)) with the event column never reassigned; forestsearch_main.R:3338 df <- FSdata$df; :3722 df.est_out <- merge(df, temp, by.x = id.name, by.y = "id", all.x = TRUE) (temp = df_flag with id + treat.recommend only, forestsearch_helpers.R:1672-1673); :3947 df.est = df.est_out. No as.integer()/storage.mode on the event column anywhere on that path (grep of forestsearch_main.R, get_fsdata.R, get_FSdata_helpers.R).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-134** (S2, off path, reproduced) `plot_subgroup_results_forestplot.R:597` plot_subgroup_results_forestplot
  - Reproducer: Appendix V174 (verdict demonstrated; wall clock in the appendix).
  - Trigger: plot_subgroup_results_forestplot(fs_results = list(fs.est, fs_bc), ..., conf.level = 0.90).
  - plot_subgroup_results_forestplot.R:282 z_alpha from conf.level; :380 (GLM :320-321) ci_label "HR (90% CI)"; :535 summary(fit, conf.int = conf.level) for ITT / reference / Full-analysis rows.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-135** (S2, off path, reproduced) `run_simulation_analysis.R:410` run_simulation_analysis
  - Reproducer: Appendix V175 (verdict demonstrated; wall clock in the appendix).
  - Trigger: run_simulation_analysis(..., cox_formula = <any formula that is not the plain treatment-only model, e.g. with strata()>).
  - grep of run_simulation_analysis.R: cox_formula (as opposed to cox_formula_adj) appears only in the roxygen (:263), formals (:410, :775, :923, :1173, :1307) and pass-through arguments (:594, :619, :644, :681, :847, :1199, :1264); no expression evaluates it. No other file in R/ passes it.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-136** (S2, off path, reproduced) `run_simulation_analysis.R:587` run_simulation_analysis
  - Reproducer: Appendix V177 (verdict demonstrated; wall clock in the appendix).
  - Trigger: run_simulation_analysis(..., methods = NULL, run_fs or run_fs_grf = TRUE, fs_params = list(n.min = NULL)).
  - run_simulation_analysis.R:587-590 (FS) and :612-615 (FSlg): fs_p <- modifyList(modifyList(fs_defaults, ...), fs_params); fs_defaults has n.min = 60 (:59) and modifyList assigns x[["n.min"]] <- NULL, which removes the element.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-137** (S2, off path, reproduced) `run_simulation_analysis.R:818` .run_fs_analysis_gen
  - Reproducer: Appendix V178 (verdict demonstrated; wall clock in the appendix).
  - Trigger: any error from forestsearch() in a replicate, e.g. an fs_params value that fails match.arg (forestsearch_main.R:1733).
  - run_simulation_analysis.R:818-821: tryCatch(do.call(forestsearch, fs_args), error = function(e) { warning(label, " failed: ", e$message); NULL }).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-138** (S2, off path, reproduced) `run_simulation_analysis.R:953` .extract_fs_estimates_gen
  - Reproducer: Appendix V037 (verdict demonstrated; wall clock in the appendix).
  - Entry point: run_simulation_analysis() (exported, run_simulation_analysis.R:389) with dgm of class glm_dgm (roxygen :204-206 lists it as accepted) and fs_params = list(outcome_type = "binary"|"count"|"continuous", effect_measure = "RD"|"RR"|"IRD" or omitted); roxygen :341-343 says these GLM parameters must be passed inside fs_params.
  - No validation or resolution of effect_measure in run_simulation_analysis(): fs_params is merged into the defaults at :587-590 / :612-615 / :674-677 (default_sim_params() :47-70 holds neither outcome_type nor effect_measure) and handed to .run_fs_analysis_gen().
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-139** (S2, off path, reproduced) `run_subgroup_sims.R:464` validate_subgroups
  - Reproducer: Appendix V179 (verdict demonstrated; wall clock in the appendix).
  - Trigger: subgroup id such as "age <= age_med" with age_med defined in the global environment and not passed in `cutpoints`; default workers = NULL (multisession, :816-821) or workers > 1.
  - run_subgroup_sims.R:799 validate_subgroups(subgroups, probe, cutpoints, benchmarks).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-140** (S2, off path, reproduced) `setup_gbsg_dgm.R:155` setup_gbsg_dgm
  - Reproducer: Appendix V038 (verdict demonstrated; wall clock in the appendix).
  - Configuration, as documented: simulate_from_dgm.R:173 formal `n = NULL`; simulate_from_dgm.R:23-24 '@param n ... If NULL (default), uses the entire super population without sampling'; simulate_from_dgm.R:126-131 '## n = NULL path'; setup_gbsg_dgm.R:25-26 'fully compatible with simulate_from_dgm'. Against that, simulate_from_dgm.R:21-22 documents `dgm` as an object 'created by generate_aft_dgm_flex', and for those objects the path works (generate_aft_dgm_helpers.R:1051 creates lin_pred_obs).
  - sim_aft_gbsg.R:592 and :598: the GBSG builder stores the assigned-arm predictor as `lin.conf.true`; sim_aft_gbsg.R:677 returns the frame as df_super_rand. grep of R/: lin_pred_obs is created only at generate_aft_dgm_helpers.R:1051 and simulate_from_dgm.R:298, and adjusted at mrct_inject_dgm.R:564.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-141** (S2, off path, reproduced) `setup_gbsg_dgm.R:175` setup_gbsg_dgm
  - Reproducer: Appendix V180 (verdict demonstrated; wall clock in the appendix).
  - setup_gbsg_dgm.R:95 and :103: cens_type = c("weibull", "uniform"), match.arg; "uniform" is an accepted, documented value (sim_aft_gbsg.R:108-109 via @inheritParams).
  - sim_aft_gbsg.R:620-625: under "uniform" no censoring model is fitted (mu_cens, sigma_cens stay NULL); :703-708 cens_params = list(type = "uniform", mu = NULL, sigma = NULL, model = NULL).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-142** (S2, off path, reproduced) `simulate_from_dgm.R:274` simulate_from_dgm
  - Reproducer: Appendix V181 (verdict demonstrated; wall clock in the appendix).
  - simulate_from_dgm.R:62-64 roxygen: strata_rand names a column for within-stratum balanced allocation; NULL means marginal allocation.
  - :274 if (!is.null(strata_rand) && strata_rand %in% names(df_sim)): a misspelt or absent column makes the condition FALSE.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-143** (S2, off path, reproduced) `simulation_tables.R:132` invert_hr_columns
  - Reproducer: Appendix V039 (verdict demonstrated; wall clock in the appendix).
  - Configuration: oc_analyses.R:629 formal `subgroup_notation = c("harm", "benefit")`, match.arg at :632, so the default is "harm"; roxygen :586-588 documents "benefit" ("labels as G/Gc for benefit-search analyses (treatment switching)"). Supported and documented, but non-default.
  - Reaching the code: format_oc_results() is exported (oc_analyses.R:617-618) and has no caller in R/ (grep finds only the definition and a header comment at :23); it is reached only by a user call on run_simulation_analysis() output. Not on the forestsearch() / bootstrap / MR path (function_index: invert_hr_columns tier 2, on_path FALSE; format_oc_results tier 2, on_path FALSE, 0 test files).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-144** (S2, off path, reproduced) `simulation_tables.R:645` build_estimation_table
  - Reproducer: Appendix V182 (verdict demonstrated; wall clock in the appendix).
  - Trigger: subgroup_notation = "benefit" + DGM without stored CDE but with theta_0/theta_1 in df_super_rand; create_gbsg_dgm() output is such a DGM (sim_aft_gbsg.R:666-673 hazard_ratios has no CDE* entries; df_super_rand = df_big carries theta_0/theta_1 from :488-489).
  - simulation_tables.R:599-602: results <- invert_hr_columns(); dgm <- invert_dgm_hrs() (1/x on hazard_ratios and top-level fields, :169-190; super-population untouched).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-145** (S2, off path, reproduced) `simulation_tables.R:1279` interpret_estimation_table
  - Reproducer: Appendix V183 (verdict demonstrated; wall clock in the appendix).
  - :1215 detect_rate computed; :1275-1285 null paragraph is a fixed template: "%d of %d simulations (%s) identified a subgroup ... This low detection rate confirms controlled type-I error." with no test on detect_rate.
  - :1226-1228 hr_H_bias may be negative; :1308-1322 prints "representing %s relative bias above the true value ... This upward bias reflects selection" regardless of sign.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-146** (S2, off path, reproduced) `summarize_bootstrap_subgroups.R:387` summarize_bootstrap_subgroups
  - Reproducer: Appendix V184 (verdict demonstrated; wall clock in the appendix).
  - maxk = 3 is a supported search depth (subgroup_search.R:412-416 accepts 1, 2 or 3).
  - The bootstrap worker stores K_sg and M.1..M.7 per iteration (bootstrap_analysis_dofuture.R:788 initialises all seven; :822-828 fill M.1..M.K; DINA path :862-868 the same).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-147** (S2, off path, reproduced) `summarize_bootstrap_subgroups.R:417` summarize_bootstrap_subgroups
  - Reproducer: Appendix V040 (verdict demonstrated; wall clock in the appendix).
  - Configuration is the default one: forestsearch_main.R:1677 pconsistency.threshold = 0.90 and :1684 pconsistency.digits = 2 (roxygen :858-860: the consistency proportion is rounded to that many decimals before comparison).
  - subgroup_consistency_helpers.R:1491, :1519, :1885, :2012: every consistency path stores round(rate, pconsistency.digits).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-148** (S2, off path, reproduced) `synthetic_data_perturbation.R:168` generate_bootstrap_synthetic
  - Reproducer: Appendix V186 (verdict demonstrated; wall clock in the appendix).
  - Trigger: a numeric ordinal variable listed in both cat_vars and ordinal_vars (the only way the branch at synthetic_data_perturbation.R:152 is reached, see C467), with a flipped row at the lowest or highest level.
  - :163-165 for the lowest level only the upper neighbour is added, for the highest only the lower neighbour, so possible_new has length 1.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-149** (S2, off path, reproduced) `synthetic_data_perturbation.R:180` generate_bootstrap_synthetic
  - Reproducer: Appendix V187 (verdict demonstrated; wall clock in the appendix).
  - synthetic_data_perturbation.R:145-149: for each cat_vars entry, flip_indices are drawn with probability cat_flip_prob (default noise_level / 2, :93-95).
  - :152 / :172-176: a variable not listed in ordinal_vars with exactly two unique values takes the binary branch.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-150** (S2, off path, reproduced) `synthetic_data_perturbation.R:528` generate_bootstrap_with_noise
  - Reproducer: Appendix V188 (verdict demonstrated; wall clock in the appendix).
  - synthetic_data_perturbation.R:514-517: continuous_vars NULL -> all numeric columns except id_var.
  - :520-525: cat_vars NULL -> factors, logicals and numeric columns with <= 10 unique values (as @param cat_vars at :395-397 states).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-151** (S2, off path, reproduced) `test_hte_crump.R:399` test_hte_from_forestsearch
  - Reproducer: Appendix V041 (verdict demonstrated; wall clock in the appendix).
  - Configuration: forestsearch() formals use_lasso = FALSE, use_grf = FALSE, cut_type = "default", subgroup_method = c("consistency", ...) (forestsearch_main.R:1638-1657); the helper's roxygen says fs_result is 'A forestsearch result object' and df is 'Data frame used in the forestsearch() call' (test_hte_crump.R:379-380). So the trigger is a default fit with at least one continuous covariate, passed as documented.
  - get_fsdata.R:286-287: non-continuous confounders stay as bare names; continuous ones become cut expressions (default cuts via get_conf_force at get_fsdata.R:401-410 and process_conf_force_expr at :518-524; median cuts 'name <= median' at :482).
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-152** (S2, off path, reproduced) `truefind_asymptotic.R:652` find_required_sample_size
  - Reproducer: Appendix V189 (verdict demonstrated; wall clock in the appendix).
  - Trigger: find_required_sample_size(theta, tol = 0.5); tol is a documented argument (truefind_asymptotic.R:575, formal default 1 at :603) and is not validated anywhere in :596-650.
  - :624-646: early returns leave the loop only with power(n_low) <= 0 and power(n_high) >= 0.
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-153** (S2, off path, reproduced) `truefind_asymptotic_glm.R:241` compute_detection_probability_glm
  - Reproducer: Appendix V042 (verdict demonstrated; wall clock in the appendix).
  - truefind_asymptotic_glm.R:179-180: '@param c2 Numeric. Consistency threshold on the natural scale. Default: 1.0 for ratio, 0.0 for difference.' (identical text at man/compute_detection_probability_glm.Rd:37-38)
  - truefind_asymptotic_glm.R:240-242: formals 'c1 = 1.25, c2 = 1.0, effect_scale = c("ratio", "difference")'; the default does not depend on effect_scale
  - (further trace steps omitted here: the reproducer demonstrates the effect)
- **A-154** (S2, off path, traced) `forestsearch_methods.R:277` print.forestsearch
  - Input: forestsearch(outcome_type = "binary") (default effect_measure OR), subgroup_method = "consistency", a subgroup found.
  - glm_effect_estimators.R:12 the estimator returns `estimate` on the log scale for OR/RR/IRR; subgroup_search.R:1083-1084 hr = res$estimate; the search table's HR column therefore holds log values (stated at forestsearch_main.R:3603-3611, which passes effect_log_scale = TRUE for exactly these measures; subgroup_consistency_main.R:242-245).
  - subgroup_consistency_helpers.R:1433-1435, 1572, 1905-1907, 2051: the result row's 'hr' is found.hrs$HR[m], unchanged; exp() is applied only to local copies for ranking/display (subgroup_consistency_helpers.R:602, 691, 875); grep finds no exponentiation of out_sg$result$hr anywhere in R/.
  - forestsearch_methods.R:276-277 prints round(top$hr, 3) as 'HR:'; :498-499 as 'Hazard ratio:'. A true OR of 1.5 is displayed as 'HR: 0.405'.
  - The same summary() exponentiates max_sg_est for these measures and labels it with effect_measure (forestsearch_methods.R:438-450), so the two lines of one summary disagree in scale.
  - Note: Corrected and UPGRADED from the reader's S3: the reader reported a wrong label only; the trace shows the value is also on the wrong scale for the ratio measures, so a harmful OR between 1 and e prints as an 'HR' below 1. Stored results are correct; only the printed line is wrong, and GLM outcomes are a non-default configuration, hence S2 rather than S1. Not run (no R allowed): a central check is cheap - fit any binary-outcome forestsearch() that finds a subgroup and compare print(fs) with exp(fs$grp.consistency$out_sg$result$hr[1]). No test pins either printed line. Identity measures (RD, MD, ...
- **A-155** (S2, off path, traced) `fs_family_report.R:293` fs_family_report
  - Other reproducer runs: Appendix V150 (no verdict line (exit 1)).
  - Trigger: fs_family_report(list(...), outcome_type = "binary") with no effect_measure; the roxygen (:31-32) says the report mirrors forestsearch()'s own argument resolution.
  - fs_family_report.R:291 em <- .arg("effect_measure") -> NULL (forestsearch formal default NULL, forestsearch_main.R:1706); :292-294 switch(outcome_type, survival = "HR", binary = "RD", ...).
  - Engine: forestsearch_main.R:1789-1794 binary = "OR"; roxygen forestsearch_main.R:1049-1050 'For binary: "OR" (default)'. The OR estimator applies the existence condition: glm_effect_estimators.R:1196-1201 (.fs_existence_reason("OR", ...) -> .fs_nonestimable).
  - fs_family_report.R:296-301 ratio_cond NULL for RD; :302-311 row 'estimability boundary' gets status 'inert' and 'RD exists on any slice' where OR would give 'all four cells ... >= 1' and 'data-dependent'.
  - :358-367 status_counts and verdict are built from the statuses, so they are off by one stage.
  - Evidence of intent against the RD fallback: NEWS.md:543-561 announces that binary now resolves to "OR" (previously "RD"); the drift-guard test (tests/testthat/test-fs-family-report.R:149) only exercises binary with an explicit effect_measure = "OR", so nothing pins the RD fallback.
  - Fitted-object input is unaffected (args_call_all carries the resolved measure).
  - Note: Kept at S2: a reporter that changes no result, but its output (a status and a count) is wrong for binary with the default effect measure, and it explicitly claims to mirror the engine's resolution. Tier 2, not on the conducted-analysis path. The other files with the same binary = "RD" fallback were listed by the reader but not traced; they may be the same leftover from before the NEWS.md:543 change and deserve one look each (I did not verify them).
- **A-156** (S2, off path, traced) `fs_fdr_report.R:258` fs_fdr_report
  - Other reproducer runs: Appendix V151 (INCONCLUSIVE).
  - fs_fdr_report.R:258-259 calibrate_k_inter(target_hr_harm = 1.0, model = "alt", use_ahr = FALSE): n_super is not forwarded (formal n_super = 25000L at :188).
  - sim_aft_gbsg.R:1023-1028 eff_at() calls .create_gbsg_dgm_(model, k_treat, k_inter, cens_type, verbose, ...) and reads dgm$hr_H_true; :242 n_super = DEFAULT_N_SUPER, :35 DEFAULT_N_SUPER <- 5000L.
  - sim_aft_gbsg.R:478-481 set.seed(seed); sample(..., size = n_super): the super-population is a different sample for a different n_super; :498 epsilon drawn with n_super; :522-530 hr_H_true is the empirical Cox HR on that sample.
  - fs_fdr_report.R:261-262 setup_gbsg_dgm(model = "alt", k_inter = k_harm, n_super = n_super) -> setup_gbsg_dgm.R:106-117 builds at n_super = 25000 (same seed default, different size).
  - Effect: hr_H_true of the simulated DGM differs from the calibrated 1.0 by the sampling difference between the two populations. It is shown in hr_structure (fs_fdr_report.R:293-295, .fs_fdr_hr_row :489) but the DGM label (:261), the roxygen (:67-69, :123-124 'sits exactly at HR = 1.0') and quarto/calibration/gbsg_aft_dgm_appendix.qmd say exactly 1.0.
  - Not established by reading: the size of the gap (the reproducer measures it).
  - Note: The mismatch itself is certain from the source; its magnitude is not, and may be small (Monte Carlo difference between a 5000- and a 25000-row resample of the same fitted model). alters_numbers = yes only for fs_fdr_report()'s harm-null rows (tier 2, not on the conducted-analysis path): k_harm would change, so the roxygen's 'published reference (n = 700): A ~ 0.138 (hom), 0.278 (harm)' (fs_fdr_report.R:163) was produced under the mismatch and the harm figure would move after a fix. The same calibrate-at-default / build-at-another-n_super pattern is used in ...
- **A-157** (S2, off path, traced) `generate_aft_dgm_helpers.R:1370` compare_multiple_survreg
  - generate_aft_dgm_helpers.R:1346 roxygen: 'Character vector of criteria to use ("AIC", "BIC", or both)'; :1370 formal criteria = c("AIC", "BIC").
  - grep 'criteria' over the file: only roxygen prose (:1341, :1346, :1353) and the formal (:1370); the body :1372-1694 never refers to it.
  - :1542 the table is sorted by AIC; :1593, :1602, :1609, :1616 set results$best_model <- best_aic in every branch.
  - Trigger: compare_multiple_survreg(fit1, fit2, criteria = "BIC"). Effect: same table, ranking and best_model as the default, i.e. an AIC-based selection although BIC was requested.
  - The package's own calls (:320, :1193) do not pass criteria, so they are unaffected.
  - Note: Facts as the reader reported; severity raised S3 -> S2 because in the non-default call criteria = "BIC" an exported function silently returns the AIC choice (misleading behaviour in a non-default configuration, brief S2). If the author regards the argument as vestigial and simply deletes it, this is S3 hygiene. Tier 2, exported, no test file mentions the function.
- **A-158** (S2, off path, traced) `mr_method_tables.R:91` mr_method_comparison
  - forestsearch_main.R:1637 formal est.scale = 'hr'; :1881 args_call_all <- mget(...) stores whatever was passed; there is no check of est.scale in forestsearch() (the three grep hits are the roxygen, the formal and an example).
  - bootstrap_dofuture_main.R:559: est.scale <- args_forestsearch_call$est.scale; :564-568: for RD / IRD / MD it is reset to 'hr'; H_estimates / Hc_estimates are computed on that scale (:574-600) and the result stores est.scale = 'hr' (:784).
  - mr_method_tables.R:89: point <- c(H_estimates$H2, Hc_estimates$H2) via .mr_method_point() (:160-165).
  - mr_method_tables.R:91-92: stored <- fit$args_call_all$est.scale; if (identical(stored, '1/hr')) point <- 1 / point -> for an RD/MD/IRD fit called with est.scale = '1/hr' the bootstrap difference is replaced by its reciprocal and printed in the bootstrap column (:99), next to naive and multiplier columns that are correct.
  - The guard at :81-83 concerns only this function's own est.scale argument and does not catch the case.
  - Note: Raised from the reader's S3 to S2 by the brief's definition: the chain is complete and the output is a wrong number, in an edge configuration. It is an unusual one: est.scale = '1/hr' has no meaning for a difference measure and forestsearch()'s roxygen lists only 'hr' or 'rmst' for est.scale (forestsearch_main.R:382, itself out of step with the '1/hr' value the rest of the package uses), but forestsearch() accepts it and the bootstrap explicitly neutralises it (:568), which is what this function fails to follow. The behaviour for ratio measures is deliberate and pinned ...
- **A-159** (S2, off path, traced) `mrct_inject_dgm.R:226` make_region_model
  - Input: make_region_model(pop, prevalence = 0.2), optionally with band = list(var, cut, or); formals at :213-214, documented as valid at :199-202.
  - :225 terms <- c(NULL, character(0)) = character(0); :226-227 ors <- c(NULL, NULL) = NULL.
  - :228 names(ors) <- terms sets an attribute on NULL (base R error); independently, :256 as.vector(scale_fn(df) %*% log(ors)) and :270 coefs = log(ors) evaluate log(NULL), also an error, so the call fails whichever is reached first.
  - :250 if (!length(terms)) return(matrix(0, nrow(df), 0)) shows the empty-terms case was meant to work; it is unreachable.
  - Only direct calls are affected: the one internal caller (:474) always passes x_pred, guarded by stopifnot(x_pred %in% continuous_vars) at :432.
  - Note: Raised from the reader's S3 to S2: an exported function fails at its documented defaults (a loud error, no wrong numbers), but only on a direct call of a rarely used helper in a degenerate configuration, hence S2 at most and not an S1 candidate. No test in tests/testthat calls make_region_model() directly. The reader supplied a two-line reproducer (formal defaults; band only) that can be run centrally.
- **A-160** (S2, off path, traced) `oc_analyses.R:648` format_oc_results
  - oc_analyses.R:641-654: when subtitle is NULL and dgm is given, hr_H <- dgm$hazard_ratios$harm_subgroup, hr_Hc <- dgm$hazard_ratios$no_harm_subgroup, formatted as '%d sims, %s(%s) = %.2f, %s(%s) = %.2f' with L$sg / L$sg_c (G / Gc under benefit notation).
  - :657-660: if (subgroup_notation == "benefit") results <- invert_hr_columns(results, effect_measure) -- only the results are inverted (1/x, or -x on identity scale; simulation_tables.R:131-143). invert_dgm_hrs() is never called in format_oc_results.
  - The sibling table functions do invert the DGM truths in the same mode (simulation_tables.R:297-301 build_classification_table, :599-602 build_estimation_table, :1105 interpret_estimation_table), which shows the intended convention.
  - Effect: body shows original-scale effects, subtitle shows e.g. 'HR(G) = <switched-scale value>'.
  - Second, unestablished effect: for a DGM whose hazard_ratios has no harm_subgroup / no_harm_subgroup (null-model aft_dgm_flex: generate_aft_dgm_helpers.R:1134-1141 adds them only when model == "alt"), sprintf() receives NULL and returns character(0), handed to gt::tab_header(subtitle = character(0)) at :815-819; what gt does with it was not checked.
  - Note: Reader's description verified; severity raised from S3 to S2 because a displayed truth value is on the wrong scale relative to its label and to the table body in a documented, non-default configuration (benefit notation + dgm + default subtitle). No computation is affected and the user can override subtitle. The roxygen (:569-573) promises 'true effect sizes using notation-aware labels', which the code does not deliver in benefit mode. No test or in-repo vignette call exercises this combination (the benefit-notation calls are in quarto/archive). The zero-length-subtitle part is suspected only ...
- **A-161** (S2, off path, traced) `plot_subgroup_results_forestplot.R:917` plot_subgroup_results_forestplot
  - Reached for effect_measure = "MD" under the default xlim_method = "clinical" (:346 xlim_deferred <- effect_measure == "MD") and for RD/IRD with xlim_method = "data" (:349); user-supplied xlim bypasses it (:352).
  - :981-983 lower_lim <- robust_lower(c(low_valid, est_valid, 0), extreme_ci_cap).
  - :912-917 robust_lower: bottom <- vals[1]; second <- vals[2]; `if (bottom < second / 2) second / cap_frac else bottom`. For bottom < second < 0, second/2 > second > bottom, so the test is always TRUE and the result is second/1.5, nearer zero than second. Example lows -3 and -2, largest value 2: lower_lim = -1.33; :992-994 pad = 0.33, xlim[1] = -1.67; :1000 rounds to -2.0, so the bound -3 lies outside the axis although it is not more than twice the second.
  - Reader's case (one negative value, second = 0): -2.5 < 0 is TRUE, lower_lim = 0/1.5 = 0, xlim[1] = -pad. Mirror on the upper side (:903-908): a lone positive value with second = 0 gives top > 0 TRUE and upper_lim = 0.
  - :1183-1186 (and :1176-1179): a bound beyond xlim is replaced by xlim in dt and printed as "<%.2f"/">%.2f"; :1269 returns the mutated dt as result$data.
  - Documented rule (roxygen :88-92): cap only if the most extreme bound exceeds twice the second-most extreme, at extreme_ci_cap times that bound.
  - Note: Corrected, not just confirmed: the reader described only the single-value-across-zero case; the lower-side test fires for every pair of negative lowest values, which is the normal situation for an MD or RD forest plot, so the scope is wider. Raised S3 -> S2 for that reason (misleading output in a non-default but supported configuration: GLM outcomes with MD, or RD/IRD with xlim_method = "data"). Not S1: the estimates are untouched, the cut bound is marked with "<"/">" and a footnote (:1196-1197), and whether a bound is actually cut depends on padding and rounding (e.g. lows -0.9/-0.8 on the ...
- **A-162** (S2, off path, traced) `synthetic_data_perturbation.R:145` generate_bootstrap_synthetic
  - Other reproducer runs: Appendix V185 (no verdict line (exit 1)).
  - synthetic_data_perturbation.R:79-82 validates ordinal_vars against the data, so the argument is accepted.
  - :145 the only categorical perturbation loop iterates for (var in cat_vars).
  - :152 the ordinal branch (var %in% ordinal_vars) sits inside that loop; a variable listed only in ordinal_vars is never visited.
  - Both roxygen examples (:20-29 grade, :42-50 education) and the exported wrapper generate_gbsg_bootstrap_general (:368-378: cat_vars = meno, hormon, status; ordinal_vars = grade) pass the ordinal variable in ordinal_vars and not in cat_vars.
  - Effect: grade is returned as the plain bootstrap resample; :217-228 only emits a note listing it among variables 'kept as-is', while @param ordinal_vars (:12-13) says these will be perturbed to adjacent values.
  - Note: Tier 2, off the conducted-analysis path; both functions are exported with no caller in R/, tests/ or vignettes and no NEWS entry. alters_numbers = yes only in the sense that synthetic datasets produced by generate_gbsg_bootstrap_general() or by the documented generate_bootstrap_synthetic() pattern would change; no forestsearch analysis result is affected. Fixing this alone makes C468 (sample() on a length-1 numeric) live for numeric ordinals such as gbsg grade, so the two should be fixed together.
- **A-163** (S2, off path, suspected) `sim_aft_gbsg.R:638` .create_gbsg_dgm_
  - sim_aft_gbsg.R:383 sigma <- fit_aft$scale (event model); not reassigned before :638 (grep for 'sigma <-' finds only :383 and :898).
  - :628-636 Weibull censoring fit; sigma_cens, mu_cens, gamma_cens extracted.
  - :637 b0_cens <- -gamma_cens / sigma_cens; :638 b_cens <- -b0_cens * sigma = gamma_cens * sigma / sigma_cens. The event-model analogue keeps the AFT coefficients unchanged (:440 b_true <- gamma; :456-458).
  - :646-647 linC1.conf / linC0.conf = z_cens %*% b_cens; used as the AFT-scale term in log C = mu_cens + sigma_cens * eps + linC.conf at :915 and, after the rename at setup_gbsg_dgm.R:157-158, at simulate_from_dgm.R:326-327, whose stated convention for lin_pred_cens_* is gamma_c'X (generate_aft_dgm_helpers.R:9-20).
  - :652-657 cens_model$gamma = gamma_cens, so the object reports coefficients that differ from those used to build the predictors.
  - Default cens_type = "weibull": every GBSG DGM from create_gbsg_dgm()/setup_gbsg_dgm() carries it. True HR/AHR values depend on the event model only and are unaffected; only the covariate dependence of simulated censoring times is scaled.
  - Missing for 'traced': the author's intent. The line dates from the original 'gbsg sims' commit (git log -S), so it may be a legacy construction carried over deliberately rather than a slip; no comment, test or NEWS entry states which.
  - Note: Tier 2, not on the conducted-analysis path (simulation DGM constructor). The arithmetic is confirmed at source; whether it is a defect depends on intent, hence evidence = suspected and no S1 candidacy: the DGM is internally usable (construction and simulation both use the same stored predictors with mu_cens and sigma_cens), no analysis estimate or stored truth is wrong, and the only established inconsistency is between the predictors and cens_params$model$gamma inside the returned object. Any change alters every simulated GBSG data set and therefore breaks the exact-reproduction property that ...

## 7. Call graph and path definition (Tier 1)

- Tier 1: forestsearch() and the 198 package functions reachable from it in the static call graph (codetools::findGlobals() over each parsed body, restricted to package functions).
- Conducted-analysis path (on_path = on): Tier 1, plus the static call graph of forestsearch_bootstrap_dofuture() (bootstrap and bias correction; 30 further functions), plus fs_mr_inference() and mr_values() with their helpers (MR inference on a forestsearch() result; 5 further functions): 234 functions. Every other function is off the path.
- The graph is static. Calls made by string or S3 dispatch count only where the symbol appears in the body; print, summary and plot methods are off the path.
- On-path rows include non-default branches of on-path functions (for example the count outcome, the DINA and GRF routes, the split fallback); the function and summary columns say which configuration triggers each.

Direct callees of `forestsearch()` (static, from the parsed source):

| Callee of forestsearch() | File | Stage | Tier-1 functions in its subtree |
|---|---|---|---|
| .announce_sg_focus | forestsearch_helpers.R | setup / validation / helpers | 2 |
| .build_grf_glm_args | grf_args.R | setup / validation / helpers | 1 |
| .build_grf_survival_args | grf_args.R | setup / validation / helpers | 1 |
| .coerce_covariates_numeric | forestsearch_helpers.R | setup / validation / helpers | 1 |
| .default_parallel_workers | forestsearch_main.R | setup / validation / helpers | 1 |
| .dina_assert_ratio_estimand | forestsearch_helpers.R | setup / validation / helpers | 1 |
| .forestsearch_dina_select | forestsearch_helpers.R | DINA identifier | 47 |
| .forestsearch_grf_select | forestsearch_helpers.R | GRF identifier | 45 |
| .fs_adjust_terms | Cox_estimation_helpers.R | setup / validation / helpers | 1 |
| .fs_adjust_vars | Cox_estimation_helpers.R | setup / validation / helpers | 2 |
| .fs_admission_applies | forestsearch_helpers.R | setup / validation / helpers | 2 |
| .fs_apply_mr | fs_mr_inference_methods.R | MR inference | 52 |
| .fs_c2_inert_note | forestsearch_helpers.R | setup / validation / helpers | 1 |
| .fs_family_status | forestsearch_helpers.R | setup / validation / helpers | 1 |
| .fs_mr_family_from_table | fs_mr_inference_methods.R | setup / validation / helpers | 4 |
| .fs_mr_reselection_from_focus | fs_mr_inference_methods.R | setup / validation / helpers | 2 |
| .fs_mr_spec | fs_mr_inference_methods.R | setup / validation / helpers | 1 |
| .fs_resolve_admission | forestsearch_helpers.R | setup / validation / helpers | 3 |
| .fs_resolve_threshold_pair | forestsearch_main.R | setup / validation / helpers | 1 |
| .make_glm_or_estimator_fast | glm_effect_estimators.R | outcome-type estimators | 5 |
| .make_lm_estimator_fast | glm_effect_estimators.R | outcome-type estimators | 2 |
| .mr_announce | fs_mr_inference_methods.R | setup / validation / helpers | 2 |
| .mr_msg | fs_mr_inference_methods.R | setup / validation / helpers | 1 |
| .mr_skip | fs_mr_inference_methods.R | setup / validation / helpers | 2 |
| .normalize_sg_focus | forestsearch_helpers.R | setup / validation / helpers | 1 |
| .resolve_dina_args | forestsearch_helpers.R | setup / validation / helpers | 2 |
| .sync_args_call_all | forestsearch_main.R | setup / validation / helpers | 1 |
| .validate_effect_neighborhood | subgroup_consistency_helpers.R | setup / validation / helpers | 1 |
| .validate_mr_configuration | forestsearch_helpers.R | setup / validation / helpers | 1 |
| .validate_outcome_threshold_config | forestsearch_main.R | setup / validation / helpers | 1 |
| .validate_selection_rule | subgroup_consistency_helpers.R | setup / validation / helpers | 1 |
| calculate_max_combinations | subgroup_search.R | setup / validation / helpers | 1 |
| dina | dina_df.R | DINA identifier | 16 |
| dina_frontier | dina_subgroup.R | DINA identifier | 6 |
| dina_subgroup | dina_subgroup.R | DINA identifier | 9 |
| estimate_propensity_scores | glm_effect_estimators.R | propensity adjustment | 5 |
| filter_call_args | summary_utility_functions.R | setup / validation / helpers | 1 |
| fs_mr_inference | fs_mr_inference.R | MR inference | 51 |
| generate_combination_indices | subgroup_search.R | setup / validation / helpers | 2 |
| get_FSdata | get_fsdata.R | candidate cuts | 15 |
| get_covs_in | subgroup_search.R | setup / validation / helpers | 1 |
| get_dfpred | forestsearch_helpers.R | setup / validation / helpers | 2 |
| get_subgroup_membership | subgroup_search.R | setup / validation / helpers | 1 |
| grf.subg.harm.glm | grf_subg_harm_glm.R | GRF identifier (GLM) | 20 |
| grf.subg.harm.survival | grf_main.R | GRF identifier (survival) | 30 |
| interpret_search_config | interpret_search_config.R | setup / validation / helpers | 2 |
| make_effect_estimator | glm_effect_estimators.R | outcome-type estimators | 13 |
| subgroup.consistency | subgroup_consistency_main.R | consistency screen and selection | 34 |
| subgroup.search | subgroup_search.R | subgroup search | 27 |
| tidy_cut_display | get_FSdata_helpers.R | setup / validation / helpers | 3 |

Tier 1 functions by file (199 functions in 26 files):

| File | Tier-1 functions |
|---|---|
| Cox_estimation_helpers.R | .fs_adjust_terms, .fs_adjust_vars |
| consistency_resample.R | .consistency_adj_terms, .consistency_cox_pieces, .consistency_glm_pieces, .dfbeta_glm, consistency_resample |
| dina.R | .dina_base_cumhaz, .dina_clip01, .dina_combine_nuisance_cox, .dina_combine_nuisance_glm, .dina_fit_baselines_cox, .dina_fit_baselines_glm, .dina_fit_censoring_cox, .dina_fit_propensity, .dina_link_and_variance, .dina_make_folds, .dina_one_pass, .dina_step2_cox, .dina_step2_glm, dina_fit |
| dina_df.R | .dina_extract_from_df, dina |
| dina_subgroup.R | .dina_check_cap, .dina_collect_candidates, .dina_collect_candidates_depth2, .dina_extract_X_from_df, .dina_warn_cap_trim, dina_frontier, dina_subgroup |
| forestsearch_helpers.R | .announce_sg_focus, .coerce_covariates_numeric, .dina_assert_ratio_estimand, .dina_reselect_on_effect, .dina_tau_sign, .forestsearch_dina_select, .forestsearch_grf_select, .fs_admission_applies, .fs_c2_inert_note, .fs_family_status, .fs_resolve_admission, .grf_reselect_on_effect, .map_dina_family, .normalize_sg_focus, .resolve_dina_args, .sg_focus_aliases_of, .validate_mr_configuration, evaluate_comparison, get_dfpred |
| forestsearch_main.R | .default_parallel_workers, .fs_resolve_threshold_pair, .sync_args_call_all, .validate_outcome_threshold_config, forestsearch |
| fs_declaration_calibration.R | .fs_decl_c0_cmp, .fs_decl_c0_shift, .fs_decl_field, .fs_decl_shifted_max, .fs_pcons_eff |
| fs_mr_field_uniform.R | fs_mr_field_uniform |
| fs_mr_inference.R | .fs_mr_assemble, .fs_mr_field_complement, .fs_mr_field_joint, .fs_mr_field_recovery, .fs_mr_ij_floor, .fs_mr_ij_var, .fs_mr_multipliers, .fs_mr_pieces, .fs_mr_se_from_ij, .fs_mr_select, fs_mr_inference |
| fs_mr_inference_methods.R | .fs_apply_mr, .fs_mr_family_from_table, .fs_mr_members_from_conj, .fs_mr_reselection_from_focus, .fs_mr_spec, .mr_announce, .mr_msg, .mr_skip |
| get_FSdata_helpers.R | collapse_cuts_fmt, collapse_cuts_round, collapse_redundant_cuts, cut_var, cut_var_jq, evaluate_cuts_once, filter_by_lassokeep, get_conf_force, get_conf_force_jq, is.continuous, is_flag_continuous, lasso_selection, process_conf_force_expr, tidy_cut_display |
| get_fsdata.R | get_FSdata |
| glm_effect_estimators.R | .build_adjusted_formula, .estimate_or, .estimate_rd, .estimate_rr, .fs_binary_cells, .fs_existence_reason, .fs_lm_two_group, .fs_logit_two_group, .fs_nonestimable, .make_cox_estimator, .make_glm_binary_estimator, .make_glm_or_estimator_fast, .make_lm_estimator, .make_lm_estimator_fast, .make_poisson_rate_estimator, .ps_grf, .ps_lasso, .ps_logistic, .robust_or_model_se, estimate_propensity_scores, make_effect_estimator |
| grf_args.R | .build_grf_glm_args, .build_grf_survival_args |
| grf_helpers.R | .grf_code_column, assign_subgroup_membership, compute_node_metrics, create_grf_config, create_null_result, create_success_result, extract_all_tree_cuts, extract_selected_tree_cuts, extract_tree_cuts, fit_causal_forest, fit_policy_trees, print_grf_details, select_best_subgroup, validate_grf_data |
| grf_main.R | grf.subg.harm.survival |
| grf_subg_harm_glm.R | .build_grf_X, .extract_tree_cuts, .select_best_sg_glm, grf.subg.harm.glm |
| grf_subgroup_labels.R | .grf_build_subgroup_definition, .grf_conjunction_labels, .grf_dr_candidates, .grf_dr_candidates_d2, .grf_evaluate_subgroup, .grf_frontier_select, .grf_leaf_conjunction, .grf_mark_frontier, .grf_parent_map, .grf_sg_def_from_candidate, .grf_simplify_conjunction |
| interpret_search_config.R | interpret_search_config |
| print_candidate_summary.R | print_candidate_preview, print_candidate_summary |
| subgroup_consistency_helpers.R | .compute_inclusion_band, .consistency_via_splits, .glm_resample_supported, .make_eval_consistency_twostage, .make_eval_subgroup_consistency, .maxeff_membership_dedup, .pareto_dominated_xy, .validate_effect_neighborhood, .validate_selection_rule, FS_labels, compute_pareto_frontier, evaluate_consistency_twostage, evaluate_subgroup_consistency, extract_subgroup, get_split_hr_fast, pareto_dominated_flags, plot_subgroup, remove_near_duplicate_subgroups, run_single_consistency_split, setup_parallel_SGcons, sg_consistency_out, sort_subgroups, sort_subgroups_preview |
| subgroup_consistency_main.R | .make_candidate_rng_seeds, subgroup.consistency |
| subgroup_search.R | .fit_or_for_subgroup_vector, .fs_cox_lean_conventions, .fs_cox_lean_ok, .fs_cox_two_group, calculate_event_counts, calculate_max_combinations, create_result_row, evaluate_combination_with_status, extract_idx_flagredundancy, fit_cox_for_subgroup, fit_glm_for_subgroup, fit_glm_for_subgroup_fast, format_search_results, generate_combination_indices, get_covs_in, get_subgroup_membership, has_positive_variance, meets_event_criteria, meets_prevalence_threshold, prepare_search_data, search_combinations_parallel, subgroup.search |
| summary_utility_functions.R | filter_call_args |
| utils_check_threads.R | .fs_grf_num_threads |

## 8. Test-gap list (dimension E)

Static reference only: a function counts as referenced when its name appears as a whole word in any file under `tests/testthat/`. S3 methods are not listed.

- **Exported Tier 1 functions with no reference** (25; finding E-001):
  - dina.R: dina_fit
  - forestsearch_helpers.R: evaluate_comparison
  - get_FSdata_helpers.R: evaluate_cuts_once, filter_by_lassokeep, get_conf_force, get_conf_force_jq, lasso_selection, process_conf_force_expr
  - glm_effect_estimators.R: estimate_propensity_scores
  - grf_helpers.R: compute_node_metrics, create_grf_config, fit_causal_forest, fit_policy_trees, print_grf_details, select_best_subgroup, validate_grf_data
  - subgroup_consistency_helpers.R: FS_labels, evaluate_consistency_twostage, get_split_hr_fast, run_single_consistency_split, setup_parallel_SGcons, sg_consistency_out
  - subgroup_search.R: get_covs_in, get_subgroup_membership
  - summary_utility_functions.R: filter_call_args
- **Tier 1 internal functions with no reference** (110; finding E-002):
  - Cox_estimation_helpers.R: .fs_adjust_terms
  - consistency_resample.R: .consistency_adj_terms, .consistency_cox_pieces
  - dina.R: .dina_base_cumhaz, .dina_clip01, .dina_combine_nuisance_cox, .dina_combine_nuisance_glm, .dina_fit_baselines_cox, .dina_fit_baselines_glm, .dina_fit_censoring_cox, .dina_fit_propensity, .dina_link_and_variance, .dina_make_folds, .dina_one_pass, .dina_step2_cox, .dina_step2_glm
  - dina_df.R: .dina_extract_from_df
  - dina_subgroup.R: .dina_check_cap, .dina_collect_candidates, .dina_collect_candidates_depth2, .dina_extract_X_from_df, .dina_warn_cap_trim
  - forestsearch_helpers.R: .coerce_covariates_numeric, .dina_reselect_on_effect, .fs_family_status, .map_dina_family, .validate_mr_configuration
  - forestsearch_main.R: .default_parallel_workers, .sync_args_call_all, .validate_outcome_threshold_config
  - fs_declaration_calibration.R: .fs_decl_shifted_max
  - fs_mr_field_uniform.R: fs_mr_field_uniform
  - fs_mr_inference.R: .fs_mr_field_complement, .fs_mr_field_joint, .fs_mr_field_recovery, .fs_mr_ij_floor, .fs_mr_pieces, .fs_mr_se_from_ij
  - fs_mr_inference_methods.R: .fs_apply_mr, .fs_mr_family_from_table, .fs_mr_members_from_conj, .fs_mr_spec, .mr_announce, .mr_msg, .mr_skip
  - get_FSdata_helpers.R: collapse_cuts_fmt, collapse_cuts_round, collapse_redundant_cuts, cut_var, cut_var_jq, is.continuous, is_flag_continuous, tidy_cut_display
  - glm_effect_estimators.R: .build_adjusted_formula, .estimate_or, .estimate_rr, .fs_lm_two_group, .fs_logit_two_group, .fs_nonestimable, .make_cox_estimator, .make_lm_estimator_fast, .make_poisson_rate_estimator, .ps_grf, .ps_lasso, .ps_logistic, .robust_or_model_se
  - grf_args.R: .build_grf_glm_args, .build_grf_survival_args
  - grf_helpers.R: assign_subgroup_membership, create_null_result, create_success_result, extract_all_tree_cuts, extract_selected_tree_cuts, extract_tree_cuts
  - grf_subg_harm_glm.R: .extract_tree_cuts, .select_best_sg_glm
  - grf_subgroup_labels.R: .grf_build_subgroup_definition, .grf_conjunction_labels, .grf_dr_candidates, .grf_dr_candidates_d2, .grf_leaf_conjunction, .grf_parent_map, .grf_sg_def_from_candidate, .grf_simplify_conjunction
  - print_candidate_summary.R: print_candidate_preview, print_candidate_summary
  - subgroup_consistency_helpers.R: .consistency_via_splits, .glm_resample_supported, .make_eval_consistency_twostage, .make_eval_subgroup_consistency, .maxeff_membership_dedup, .validate_effect_neighborhood, .validate_selection_rule, compute_pareto_frontier, extract_subgroup, pareto_dominated_flags, plot_subgroup, remove_near_duplicate_subgroups
  - subgroup_search.R: .fit_or_for_subgroup_vector, calculate_event_counts, calculate_max_combinations, create_result_row, evaluate_combination_with_status, extract_idx_flagredundancy, generate_combination_indices, has_positive_variance, meets_event_criteria, meets_prevalence_threshold, prepare_search_data, search_combinations_parallel
- **Exported on-path Tier 2 functions (bootstrap and MR path) with no reference** (9; finding E-003):
  - Cox_estimation_helpers.R: build_cox_formula, cox_summary, fit_cox_models, get_Cox_sg
  - bootstrap_calculations_helpers.R: format_CI, get_dfRes
  - summary_utility_functions.R: hrCI_format, n_pcnt, rmst_calculation
- **Internal on-path Tier 2 functions with no reference** (19; finding E-003):
  - Cox_estimation_helpers.R: fit_effect_models, fit_subgroup_effect
  - bootstrap_calculations_helpers.R: .boot_indices_to_ystar, .make_boot_index_matrix, calc_cov, ci_est, get_targetEst
  - bootstrap_dofuture_setup_helpers.R: ensure_packages, resolve_bootstrap_parallel_args
  - mr_values.R: .mr_get, .mr_int, .mr_num, .mr_values_input
  - summary_utility_functions.R: SG_tab_estimates_glm, analyze_subgroup_glm, calculate_counts, calculate_potential_hr, format_results, km_summary
- **Exported off-path functions with no reference** (117; finding E-004):
  - bootstrap_calculations_helpers.R: bootstrap_ystar, count_boot_id
  - bootstrap_summaries_helpers.R: format_bootstrap_table, summarize_bootstrap_events
  - calibrate_glm_interaction.R: calibrate_glm_interaction
  - calibrate_k_treat.R: calibrate_k_treat
  - calibrate_null_correction.R: calibrate_L_eff, predict_fpr_corrected, run_null_calibration
  - compare_fs_subgroups.R: compare_fs_subgroups
  - compare_selection_rules.R: compare_selection_rules, extract_selected_subgroups
  - cox_ahr_cde_wrapper.R: cox_ahr_cde_analysis
  - cox_spline_fit.R: cox_cs_fit, plot_subgroup_effects
  - create_summary_table.R: create_summary_table
  - cv_summary_tables.R: cv_compare_results, cv_metrics_tables, cv_summary_tables, cv_summary_text
  - dina_bagged.R: dina_fit_bagged
  - dina_df.R: dina_bagged
  - dina_subgroup_refit.R: dina_subgroup_refit
  - explain_pareto_selection.R: explain_pareto_selection
  - find_k_inter_main.R: find_k_inter_for_target_hr, sensitivity_analysis_k_inter
  - forestsearch_cross_validation.R: CV_sgs
  - forestsearch_helpers.R: collect_results, reset_workers
  - format_fs_details.R: format_fs_details
  - fpr_approximation.R: fpr_approx
  - fpr_calibration.R: fpr_calibration, fpr_calibration_transfer
  - frontier_cis.R: compute_frontier_cis
  - frontier_member_flags.R: frontier_member_flags
  - fs_bias_coverage.R: fs_plot_bias_coverage
  - fs_fdr_report.R: fs_fdr_report
  - fs_identification_summary.R: fs_identification_structure, fs_identification_summary, fs_identification_table, fs_rule_covariate_pairs, plot_fs_identification_involvement, plot_fs_identification_structure, plot_fs_rule_pairs
  - fs_to_guohe.R: fs_assert_membership, fs_export_maxeff_family, fs_to_guohe
  - generate_aft_dgm_helpers.R: compare_multiple_survreg
  - generate_aft_dgm_helpers_outcome.R: plot_spline_treatment_effect
  - get_FSdata_helpers.R: dummy_encode
  - glm_effect_profile.R: glm_effect_profile
  - glm_spline_fit.R: glm_cs_fit
  - grf_args.R: extract_grf_args
  - grf_helpers.R: find_leaf_split
  - grf_main.R: grf.subg.eval
  - grf_subg_harm_glm.R: create_glm_row
  - guohe_adaptive_r.R: guohe_adaptive_r
  - guohe_algorithm3.R: guohe_algorithm3
  - leff_calibration_io.R: get_leff, latest_rds, load_leff_calibration
  - mrct_inject_dgm.R: make_region_model, plot_region_modifier
  - mrct_simulation.R: SGplot_estimates, create_dgm_for_mrct, summaryout_mrct, validate_mrct_data
  - oc_analyses.R: compute_aor, compute_cde_glm, format_oc_results, summarize_simulation_results
  - pareto_frontier_table.R: pareto_frontier_table
  - plot_effect_distribution.R: plot_effect_distribution
  - plot_km_band_forestsearch.R: plot_km_band_forestsearch, quick_km_band_plot
  - plot_pareto_combined.R: plot_pareto_combined
  - plot_pareto_frontier.R: plot_pareto_frontier
  - plot_sg_distribution.R: plot_sg_distribution, sgdist_fig_height
  - plot_sg_weighted_km.R: figure_note
  - render_forestplot.R: create_forest_theme, render_forestplot, save_forestplot
  - run_simulation_analysis.R: run_simulation_analysis
  - run_subgroup_sims.R: subgroup_glm
  - sim_aft_gbsg.R: create_gbsg_dgm, simulate_from_gbsg_dgm, validate_k_inter_effect
  - simulate_from_dgm.R: calibrate_cens_adjust, check_censoring_dgm
  - simulation_tables.R: build_classification_table, build_estimation_table, interpret_estimation_table, render_reference_table
  - subgroup_consistency_helpers.R: early_stop_decision, wilson_ci
  - subgroup_search.R: get_combinations_info
  - suggest_thresholds.R: suggest_thresholds
  - summarize_bootstrap_results.R: summarize_bootstrap_results
  - summarize_bootstrap_subgroups.R: summarize_bootstrap_subgroups
  - synthetic_data_perturbation.R: generate_bootstrap_synthetic, generate_bootstrap_with_noise, generate_gbsg_bootstrap_general
  - test_hte_crump.R: test_constant_cate, test_hte, test_hte_from_forestsearch, test_zero_cate
  - truefind_asymptotic.R: compare_detection_curves, compute_detection_probability, find_required_sample_size, generate_detection_curve, plot_detection_curve
  - truefind_asymptotic_glm.R: compute_detection_probability_glm, d_eff_binary, d_eff_continuous, d_eff_count, d_eff_survival

## 9. Maintainability: ten worst

Ranked by the sum of ranks of length, argument count and brace-nesting depth across all 793 functions. One row each; the rows also appear in the ranked table as G-001 to G-010 (S4 by convention: a structural cost, no defect).

| ID | Function | File:line | Lines | Arguments | Max brace depth | Mixed responsibilities |
|---|---|---|---|---|---|---|
| G-001 | subgroup.consistency() | subgroup_consistency_main.R:376 | 834 | 33 | 8 | Argument resolution, seeding, a sequential loop and a separate batched parallel loop with early stopping, the maxeff path, and result assembly in one body; the two evaluation loops duplicate logic. |
| G-002 | forestsearch() | forestsearch_main.R:1622 | 2362 | 86 | 6 | Validation, alias and threshold resolution, propensity scores, the GRF and DINA selectors, cut construction, search, consistency, re-selection, MR and return assembly in one function. |
| G-003 | create_summary_table() | create_summary_table.R:75 | 819 | 27 | 8 | Data assembly, statistics and gt formatting in one function. |
| G-004 | plot_subgroup_results_forestplot() | plot_subgroup_results_forestplot.R:195 | 1083 | 28 | 6 | Row construction, effect fitting, CV and bootstrap text, theme and forestploter drawing in one function. |
| G-005 | mrct_region_sims() | mrct_simulation.R:260 | 484 | 20 | 8 | Plan setup, DGM simulation, analysis, MR and per-replicate record assembly in one function. |
| G-006 | evaluate_consistency_twostage() | subgroup_consistency_helpers.R:1654 | 406 | 22 | 6 | Two-stage screening, Wilson stopping, the closed-form path and an embedded copy of the split routine. |
| G-007 | run_simulation_analysis() | run_simulation_analysis.R:389 | 377 | 28 | 5 | Simulation, ForestSearch and GRF runs, and estimate extraction for survival and GLM outcomes in one function. |
| G-008 | plot_sg_weighted_km() | plot_sg_weighted_km.R:94 | 372 | 19 | 5 | Subgroup extraction, KM computation, layout and annotation in one function. |
| G-009 | fs_mr_inference() | fs_mr_inference.R:584 | 557 | 29 | 4 | Argument resolution, family enumeration, several variance constructions, field calibration, admission and output shaping in one function. |
| G-010 | bootstrap_results() | bootstrap_analysis_dofuture.R:300 | 737 | 12 | 6 | Replicate setup, the inner forestsearch() call, four subgroup fits, bias corrections, MR capture and timing in one foreach body. |

## 10. Closed items encountered

Listed once each with a location; no severity. None was entered in the ranked table.

- Existence of the deprecated arguments vi.grf.min, max_n_confounders, max.minutes: forestsearch_main.R:1648, 1689, 1698 (warnings :1744-1768; max.minutes dropped at :3493); subgroup_search.R:100; reported as inert by fs_family_report.R:266, 353. Their inertness is invariant 4 (holds).
- Winner-only and winner-floor IJ variance paths: fs_mr_inference.R:763-767, 849-853, 869-881, 1078-1088, 1409-1413; forestsearch_main.R:1008-1013, 3901; fs_mr_inference_methods.R:188; mr_values.R:138-144; fs_bias_coverage.R:96-103, 130-141.
- field_uniform research option: fs_mr_inference.R:599, 1016-1045; fs_mr_field_uniform.R:79-254; forestsearch_main.R:997-1000, 3873; fs_mr_inference_methods.R:180. Three defects inside it were found and verified and are recorded here, not ranked: a length-1 t_g indexed as length K on the effect-floor-only admission branch (fs_mr_inference.R:1022, verifier S2); kappa = 1 returned with no note when no hypothetical trial is usable (fs_mr_field_uniform.R:228, S3); roxygen says kappa in [1, 2] and M <= 12 while the code uses a 1 to 3 grid and a cap of 40 (fs_mr_inference.R:332, S3).
- Unscaled and `_s` complement bounds and the joint / joint_s pair (kept on purpose as research options): fs_mr_inference.R:1308-1452; forestsearch_methods.R:88-106; mrct_simulation.R:798-801; mr_values.R:132-133, 187-203; mrct_case_study.R:246-248.
- mr_estimates_table()'s "mr_values" attribute reading stored limits: mr_estimates_table.R:246-255, 336.
- The case study's MR print block reading stored limits: mrct_case_study.R:347, mrct_trial_helpers.R:197-209.
- No labels argument on the KM plot: plot_sg_weighted_km.R:94; mrct_case_study.R:369-377, mrct_trial_helpers.R:252-279.
- mrct_record() reads CSV only: mrct_minor_helpers.R:27.
- Not encountered: the removal of fit_causal_forest_glm(); articles with untracked HTML (outside R/).

Known-and-recorded cross-check: none of the closed items above appears as a finding. Two verifier items were classed out of scope (methodology or efficiency) and six were refuted at source; they are not in the table.

## 11. Coverage statement

- Tier 1, read line by line (26 files, by owner): Cox_estimation_helpers.R, consistency_resample.R, dina.R, dina_df.R, dina_subgroup.R, forestsearch_helpers.R, forestsearch_main.R, fs_declaration_calibration.R, fs_mr_field_uniform.R, fs_mr_inference.R, fs_mr_inference_methods.R, get_FSdata_helpers.R, get_fsdata.R, glm_effect_estimators.R, grf_args.R, grf_helpers.R, grf_main.R, grf_subg_harm_glm.R, grf_subgroup_labels.R, interpret_search_config.R, print_candidate_summary.R, subgroup_consistency_helpers.R, subgroup_consistency_main.R, subgroup_search.R, summary_utility_functions.R, utils_check_threads.R. forestsearch_main.R was read three times over (roxygen and helpers; body 1622-2800; body 2800-3983), with each body reader also reading the other range in full.
- Tier 2, swept (the other 97 files): every file in R/ was assigned to exactly one reader; the files on the bootstrap and MR path (bootstrap_dofuture_main.R, bootstrap_analysis_dofuture.R, bootstrap_calculations_helpers.R, bootstrap_dofuture_setup_helpers.R, mr_values.R, and the on-path functions of summary_utility_functions.R and Cox_estimation_helpers.R) were read line by line.
- Cross-cutting passes over all 125 files: invariants 2, 4 and 5 (grep for every plan/options/RNG/par/thread site, each enclosing function read), invariant 3 (all roxygen @param blocks against formals, and all man/*.Rd \arguments against \usage), and the main-session static sweep.
- Nothing in R/ was left unread. Limits: the readers did not run code, so draw-time behaviour of ggplot2, grid, gt and forestploter was judged from source only (two plotting findings say so); the internals of future, doFuture and progressr were not read; man/ pages were compared only for defaults and arguments.
- Reproducers: every S1 and S1-candidate finding and every S2 finding with a script were run until the budget allowed (see the header and appendix); findings without a run keep their traced or suspected class.

## 12. Appendix: reproducers

Each script ran with `Rscript --vanilla` from `tempdir()` after `pkgload::load_all()` on the checked-out tree, sequential plan unless the script states otherwise, under `timeout 120`. Output is the tail of the run (at most 25 lines). Scripts V001 to V190 were written by the verifiers and are pasted unchanged; V087b and V110b are V087 and V110 with the missing package load added; R01 to R04 were written in the main session. A verdict line in capitals (CONFIRMED, REPRODUCED, DEFECT_PRESENT) is the verifier's own wording for demonstrated; every verdict is computed from the observed values in the script.

### R01_bootstrap_id_merge

- Finding: A-002 — Bootstrap replicate built by hand as bootstrap_results() builds it: frame sizes and the two H* estimates
- Runs: exit 0, 3.0 s. Verdict: demonstrated

```r
suppressMessages(pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE))
future::plan(future::sequential)
set.seed(1); n <- 600
d <- data.frame(id = seq_len(n), x1 = rbinom(n, 1, 0.5), x2 = rnorm(n), x3 = rnorm(n), treat = rbinom(n, 1, 0.5))
t0 <- rexp(n, 0.1 * exp(ifelse(d$x1 == 1, log(3), log(0.7)) * d$treat)); cc <- runif(n, 5, 30)
d$tte <- pmin(t0, cc); d$event <- as.integer(t0 <= cc)
idx <- sample.int(n, n, replace = TRUE); k <- tabulate(idx, n)
boot <- d[idx, ]   # what bootstrap_results() builds at bootstrap_analysis_dofuture.R:438 (ids kept, repeated)
fs <- suppressWarnings(suppressMessages(forestsearch(df.analysis = boot, df.predict = d,
  confounders.name = c("x1", "x2", "x3"), outcome.name = "tte", event.name = "event", treat.name = "treat", id.name = "id",
  use_lasso = FALSE, use_grf = FALSE, maxk = 1,
  parallel_args = list(plan = "sequential", workers = 1, show_message = FALSE), details = FALSE)))
stopifnot(!is.null(fs$sg.harm))
cat("sg.harm:", paste(fs$sg.harm, collapse = " & "), "\n")
cat("df.est rows:", nrow(fs$df.est), "| bootstrap n:", nrow(boot), "| sum(k^2):", sum(k^2), "\n")
cat("df.predict rows:", nrow(fs$df.predict), "| original n:", nrow(d), "| sum(k) + out-of-bag:", sum(k) + sum(k == 0), "\n")
cat("NA treat.recommend in df.predict:", sum(is.na(fs$df.predict$treat.recommend)), "| out-of-bag ids:", sum(k == 0), "\n")
ref_orig <- get_dfpred(d, fs$sg.harm, version = 2)
ref_boot <- get_dfpred(boot, fs$sg.harm, version = 2)
hr <- function(df) unname(exp(coef(survival::coxph(survival::Surv(tte, event) ~ treat, data = df))))
H_orig_rule  <- hr(subset(ref_orig, treat.recommend == 0))             # H* on ORIGINAL data (rule evaluated row-wise)
H_orig_used  <- hr(subset(fs$df.predict, treat.recommend == 0))        # what bootstrap_results() uses as 'H* on original'
H_boot_rule  <- hr(subset(ref_boot, treat.recommend == 0))             # H* on the BOOTSTRAP sample
H_boot_used  <- hr(subset(fs$df.est, treat.recommend == 0))            # what bootstrap_results() uses as 'H* on bootstrap'
cat(sprintf("HR(H*) on original data : by rule %.4f | from fs$df.predict %.4f\n", H_orig_rule, H_orig_used))
cat(sprintf("HR(H*) on bootstrap data: by rule %.4f | from fs$df.est     %.4f\n", H_boot_rule, H_boot_used))
cat(sprintf("fs$df.predict-based 'original' HR equals the bootstrap-sample HR: %s\n", isTRUE(all.equal(H_orig_used, H_boot_rule))))
demo <- nrow(fs$df.est) == sum(k^2) && nrow(fs$df.est) != nrow(boot) && sum(is.na(fs$df.predict$treat.recommend)) == sum(k == 0) &&
  !isTRUE(all.equal(H_orig_rule, H_orig_used)) && !isTRUE(all.equal(H_boot_rule, H_boot_used))
cat("REPRO_VERDICT:", if (demo) "demonstrated" else "not demonstrated", "\n")
```

```text
sg.harm: {x1} 
df.est rows: 1176 | bootstrap n: 600 | sum(k^2): 1176 
df.predict rows: 817 | original n: 600 | sum(k) + out-of-bag: 817 
NA treat.recommend in df.predict: 217 | out-of-bag ids: 217 
HR(H*) on original data : by rule 2.9590 | from fs$df.predict 3.1306
HR(H*) on bootstrap data: by rule 3.1306 | from fs$df.est     3.1221
fs$df.predict-based 'original' HR equals the bootstrap-sample HR: TRUE
REPRO_VERDICT: demonstrated
```

### R02_bootstrap_end_to_end

- Finding: A-002 — End to end: forestsearch_bootstrap_dofuture(nb_boots = 4) with forestsearch() wrapped to record inputs and outputs (first attempts failed on the script; every attempt is counted)
- Runs: exit 0, 8.1 s; exit 1, 7.8 s; exit 1, 0.1 s; exit 0, 7.7 s; exit 0, 7.7 s; exit 0, 7.8 s. Verdict: demonstrated

```r
# End-to-end: what the real bootstrap hands to, and reads back from, forestsearch() in each replicate.
suppressMessages(pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE))
future::plan(future::sequential)
set.seed(1); n <- 600
d <- data.frame(id = seq_len(n), x1 = rbinom(n, 1, 0.5), x2 = rnorm(n), x3 = rnorm(n), treat = rbinom(n, 1, 0.5))
t0 <- rexp(n, 0.1 * exp(ifelse(d$x1 == 1, log(3), log(0.7)) * d$treat)); cc <- runif(n, 5, 30)
d$tte <- pmin(t0, cc); d$event <- as.integer(t0 <= cc)
fs <- suppressWarnings(suppressMessages(forestsearch(df.analysis = d,
  confounders.name = c("x1", "x2", "x3"), outcome.name = "tte", event.name = "event", treat.name = "treat", id.name = "id",
  use_lasso = FALSE, use_grf = FALSE, maxk = 1,
  parallel_args = list(plan = "sequential", workers = 1, show_message = FALSE), details = FALSE)))
stopifnot(!is.null(fs$sg.harm))
cat("original fit: sg.harm =", paste(fs$sg.harm, collapse = " & "), "| nrow(df.est) =", nrow(fs$df.est), "\n")
rec <- list()
ns <- asNamespace("forestsearch")
orig_fs <- get("forestsearch", envir = ns)
wrapper <- function(...) {          # records inputs and outputs, changes nothing
  a <- list(...)
  rv <- orig_fs(...)
  idn <- if (is.null(a$id.name)) "id" else a$id.name
  ids <- a$df.analysis[[idn]]
  k <- as.integer(table(ids))
  r <- get("rec", envir = globalenv())
  r[[length(r) + 1L]] <- c(
    n_in = nrow(a$df.analysis), dup_ids = sum(duplicated(ids)), sum_k2 = sum(k^2),
    found = !is.null(rv$sg.harm),
    n_df_est = if (is.null(rv$df.est)) NA else nrow(rv$df.est),
    n_predict_in = if (is.null(a$df.predict)) NA else nrow(a$df.predict),
    n_df_predict = if (is.null(rv$df.predict)) NA else nrow(rv$df.predict),
    na_flag_predict = if (is.null(rv$df.predict)) NA else sum(is.na(rv$df.predict$treat.recommend)),
    oob = if (is.null(a$df.predict)) NA else sum(!(a$df.predict[[idn]] %in% ids)))
  assign("rec", r, envir = globalenv())
  rv
}
pe <- as.environment("package:forestsearch")   # futures resolve package globals through the search path
unlockBinding("forestsearch", ns); assign("forestsearch", wrapper, envir = ns)
if (bindingIsLocked("forestsearch", pe)) unlockBinding("forestsearch", pe); assign("forestsearch", wrapper, envir = pe)
wl <- character(); bt <- withCallingHandlers(suppressMessages(forestsearch_bootstrap_dofuture(fs.est = fs, nb_boots = 4, seed = 11L,
  parallel_args = list(plan = "sequential", workers = 1, show_message = FALSE))), warning = function(w) { wl <<- c(wl, conditionMessage(w)); invokeRestart("muffleWarning") })
cat("warnings:", length(wl), "\n"); cat(substr(unique(wl), 1, 300), sep = "\n"); print(as.data.frame(bt$results)[, 1:min(8, ncol(bt$results))])
assign("forestsearch", orig_fs, envir = ns); assign("forestsearch", orig_fs, envir = pe); lockBinding("forestsearch", ns)
cat("calls recorded:", length(rec), "\n")
tab <- as.data.frame(do.call(rbind, rec)); print(tab)
ok <- tab[tab$found == 1 & !is.na(tab$n_df_est), ]
demo <- nrow(ok) > 0 && all(ok$dup_ids > 0) && all(ok$n_df_est == ok$sum_k2) && all(ok$n_df_est > ok$n_in) &&
        all(ok$na_flag_predict == ok$oob) && all(ok$n_df_predict > ok$n_predict_in)
cat("replicates with a subgroup:", nrow(ok), "\n")
cat("REPRO_VERDICT:", if (demo) "demonstrated" else "not demonstrated", "\n")
```

```text
warnings: 1 
nb_boots < 100 may produce unreliable confidence intervals. Recommend at least 500 for production use.
  boot_id H_biasadj_1 H_biasadj_2 Hc_biasadj_1 Hc_biasadj_2 max_sg_est  L
1       1   1.1211469    1.189380   -0.2670052 -0.400418338   2.763831 16
2       2   0.9773275    1.047181   -0.2848833 -0.131660405   2.759358 14
3       3   1.1273308    1.202030   -0.1621497 -0.002727427   2.746018 16
4       4   1.1226248    1.275705   -0.2389456 -0.059283083   2.539002 14
  max_count
1        16
2        14
3        16
4        14
calls recorded: 4 
  n_in dup_ids sum_k2 found n_df_est n_predict_in n_df_predict na_flag_predict
1  600     214   1176     1     1176          600          814             214
2  600     217   1172     1     1172          600          817             217
3  600     215   1160     1     1160          600          815             215
4  600     226   1214     1     1214          600          826             226
  oob
1 214
2 217
3 215
4 226
replicates with a subgroup: 4 
REPRO_VERDICT: demonstrated
```

### R03_inv5_seq_vs_multisession

- Finding: INV-5 — Invariant 5: forestsearch() under sequential and multisession (2 workers), resample and split, one dataset
- Runs: exit 0, 34.7 s. Verdict: selection identical in both; Pcons differs between plans under split

```r
# Invariant 5: same input and seed, sequential vs multisession (2 workers).
suppressMessages(pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE))
old_plan <- future::plan(); on.exit(future::plan(old_plan), add = TRUE)
set.seed(1); n <- 600
d <- data.frame(id = seq_len(n), x1 = rbinom(n, 1, 0.5), x2 = rnorm(n), x3 = rnorm(n), treat = rbinom(n, 1, 0.5))
t0 <- rexp(n, 0.1 * exp(ifelse(d$x1 == 1, log(3), log(0.7)) * d$treat)); cc <- runif(n, 5, 30)
d$tte <- pmin(t0, cc); d$event <- as.integer(t0 <= cc)
run <- function(plan, workers, method) {
  fs <- suppressWarnings(suppressMessages(forestsearch(df.analysis = d,
    confounders.name = c("x1", "x2", "x3"), outcome.name = "tte", event.name = "event", treat.name = "treat", id.name = "id",
    use_lasso = FALSE, use_grf = FALSE, maxk = 2, consistency_method = method, fs.splits = 200,
    parallel_args = list(plan = plan, workers = workers, show_message = FALSE), details = FALSE)))
  g <- fs$grp.consistency$out_sg$result
  if (is.null(g)) g <- fs$grp.consistency$result
  list(sg = paste(fs$sg.harm, collapse = " & "), tab = as.data.frame(g))
}
cmp <- function(method) {
  a <- run("sequential", 1, method); b <- run("multisession", 2, method)
  cat("\n== consistency_method =", method, "\n")
  cat("sequential  : sg.harm =", a$sg, "\n"); cat("multisession: sg.harm =", b$sg, "\n")
  pc <- intersect(c("Pcons", "pcons", "Pconsistency"), names(a$tab))
  if (length(pc) && nrow(a$tab) && nrow(b$tab)) {
    k <- min(nrow(a$tab), nrow(b$tab), 6)
    cat("Pcons (first candidates) sequential  :", format(a$tab[[pc[1]]][seq_len(k)], digits = 6), "\n")
    cat("Pcons (first candidates) multisession:", format(b$tab[[pc[1]]][seq_len(k)], digits = 6), "\n")
    same_p <- isTRUE(all.equal(a$tab[[pc[1]]], b$tab[[pc[1]]]))
  } else { cat("columns:", paste(names(a$tab), collapse = ", "), "\n"); same_p <- NA }
  cat("same selection:", identical(a$sg, b$sg), "| same Pcons:", same_p, "\n")
  c(sel = identical(a$sg, b$sg), pcons = same_p)
}
r1 <- cmp("resample"); r2 <- cmp("split")
future::plan(future::sequential)
cat("\nDefault (resample / closed form): selection identical =", r1[["sel"]], "; Pcons identical =", r1[["pcons"]], "\n")
cat("Split: selection identical =", r2[["sel"]], "; Pcons identical =", r2[["pcons"]], "\n")
cat("REPRO_VERDICT:", if (!isTRUE(r1[["sel"]]) || !isTRUE(r2[["sel"]])) "demonstrated (selection differs)" else
    if (identical(r2[["pcons"]], FALSE)) "selection identical in both; Pcons differs between plans under split" else "not demonstrated", "\n")
```

```text

== consistency_method = resample 
sequential  : sg.harm = {x3 <= -0.7} & {x1} 
multisession: sg.harm = {x3 <= -0.7} & {x1} 
Pcons (first candidates) sequential  : 1 1 1 1 1 1 
Pcons (first candidates) multisession: 1 1 1 1 1 1 
same selection: TRUE | same Pcons: TRUE 

== consistency_method = split 
sequential  : sg.harm = {x3 <= -0.7} & {x1} 
multisession: sg.harm = {x3 <= -0.7} & {x1} 
Pcons (first candidates) sequential  : 1 1 1 1 1 1 
Pcons (first candidates) multisession: 1 1 1 1 1 1 
same selection: TRUE | same Pcons: FALSE 

Default (resample / closed form): selection identical = TRUE ; Pcons identical = TRUE 
Split: selection identical = TRUE ; Pcons identical = FALSE 
REPRO_VERDICT: selection identical in both; Pcons differs between plans under split
```

### R04_gbsg_data_call

- Finding: A-004 — generate_gbsg_bootstrap_general(): the data() call loads nothing and every call fails
- Runs: exit 0, 2.3 s. Verdict: demonstrated

```r
# generate_gbsg_bootstrap_general(): data("gbsg", package = "survival") at synthetic_data_perturbation.R:366
suppressMessages(pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE))
cat("survival version:", as.character(packageVersion("survival")), "\n")
e <- new.env()
suppressWarnings(data("gbsg", package = "survival", envir = e))
cat("data(\"gbsg\", package = \"survival\") loads gbsg:", exists("gbsg", envir = e, inherits = FALSE), "\n")
cat("survival::gbsg is available:", is.data.frame(survival::gbsg), "with", nrow(survival::gbsg), "rows\n")
r <- tryCatch(suppressWarnings(generate_gbsg_bootstrap_general(n = 50, seed = 1)), error = function(err) err)
cat("generate_gbsg_bootstrap_general(n = 50, seed = 1):", if (inherits(r, "error")) paste("ERROR -", conditionMessage(r)) else paste("returned", nrow(r), "rows"), "\n")
r0 <- tryCatch(suppressWarnings(generate_gbsg_bootstrap_general()), error = function(err) err)
cat("generate_gbsg_bootstrap_general() with defaults:", if (inherits(r0, "error")) paste("ERROR -", conditionMessage(r0)) else paste("returned", nrow(r0), "rows"), "\n")
cat("REPRO_VERDICT:", if (inherits(r, "error") && inherits(r0, "error") && !exists("gbsg", envir = e, inherits = FALSE)) "demonstrated" else "not demonstrated", "\n")
```

```text
survival version: 3.8.12 
data("gbsg", package = "survival") loads gbsg: FALSE 
survival::gbsg is available: TRUE with 686 rows
generate_gbsg_bootstrap_general(n = 50, seed = 1): ERROR - 'data' must be a data frame 
generate_gbsg_bootstrap_general() with defaults: ERROR - 'data' must be a data frame 
REPRO_VERDICT: demonstrated
```

### V087b

- Finding: D-001 — V087 rerun with the missing pkgload::load_all() line added (V087 failed on the script)
- Runs: exit 0, 2.6 s. Verdict: demonstrated

```r
suppressMessages(pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE))  # added: the verifier script omitted the load
## INV-1-a: dina_cuts enter the pool with every front-end flag FALSE; grf_cuts do not.
## Run from tempdir() after pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1)
n <- 300
d <- data.frame(id = seq_len(n), tte = rexp(n), event = rbinom(n, 1, 0.7),
                treat = rbinom(n, 1, 0.5),
                age = round(rnorm(n, 55, 10), 1),
                bm = round(rnorm(n, 20, 5), 1),
                sexf = rbinom(n, 1, 0.5))
conf <- c("age", "bm", "sexf")
cut <- "age <= 51.37"
base <- get_FSdata(d, use_lasso = FALSE, use_grf = FALSE, confounders.name = conf, details = FALSE)
g    <- get_FSdata(d, use_lasso = FALSE, use_grf = FALSE, grf_cuts = cut, confounders.name = conf, details = FALSE)
dn   <- get_FSdata(d, use_lasso = FALSE, use_grf = FALSE, dina_cuts = cut, confounders.name = conf, details = FALSE)
cat("baseline candidates:\n"); print(base$confs)
cat("with grf_cuts, use_grf = FALSE:\n"); print(g$confs)
cat("with dina_cuts, no flag:\n"); print(dn$confs)
in_dina <- any(grepl("51.37", dn$confs, fixed = TRUE))
in_grf  <- any(grepl("51.37", g$confs, fixed = TRUE))
status <- .fs_family_status(list(use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE,
                                 subgroup_method = "consistency", dina_cuts = cut))
cat("dina cut in pool:", in_dina, "| grf cut in pool:", in_grf,
    "| grf pool identical to baseline:", identical(g$confs, base$confs),
    "| family_status:", status, "\n")
demonstrated <- isTRUE(in_dina) && !isTRUE(in_grf) && identical(g$confs, base$confs) &&
  identical(status, "no-front-end")
cat("REPRO_VERDICT:", if (demonstrated) "demonstrated" else "not demonstrated", "\n")
```

```text
baseline candidates:
[1] "age <= 54" "age <= 46" "age <= 62" "bm <= 20"  "bm <= 17"  "bm <= 23" 
[7] "sexf"     
with grf_cuts, use_grf = FALSE:
[1] "age <= 54" "age <= 46" "age <= 62" "bm <= 20"  "bm <= 17"  "bm <= 23" 
[7] "sexf"     
with dina_cuts, no flag:
[1] "age <= 51.37" "age <= 54"    "age <= 46"    "age <= 62"    "bm <= 20"    
[6] "bm <= 17"     "bm <= 23"     "sexf"        
dina cut in pool: TRUE | grf cut in pool: FALSE | grf pool identical to baseline: TRUE | family_status: no-front-end 
REPRO_VERDICT: demonstrated
```

### V110b

- Finding: INV-5 — V110 rerun with the missing pkgload::load_all() line added (V110 failed on the script): invariant 5 over three datasets
- Runs: exit 0, 28.7 s. Verdict: demonstrated (sequential and multisession select different subgroups for the same seed in: B: consistency_method='split', data seed 203 )

```r
suppressMessages(pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE))  # added: the verifier script omitted the load
## INV-5 reproducer: forestsearch() with parallel_args plan = "sequential" vs "multisession", same seed.
## Run from tempdir() after pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE).
## A = default consistency_method ("resample"): expected identical.
## B = consistency_method = "split": Pcons expected to differ; selection may differ.
## REPRO_VERDICT 'demonstrated' means the two plans SELECT DIFFERENTLY in at least one comparison.
t_start <- proc.time()[["elapsed"]]
old_plan <- future::plan()
had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
old_seed <- if (had_seed) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL

inst <- suppressWarnings(tryCatch(utils::packageDescription("forestsearch", lib.loc = .libPaths()),
                                  error = function(e) NA))
if (is.list(inst)) {
  cat("installed forestsearch (the copy multisession workers load): Version", inst$Version,
      "| Packaged", inst$Packaged, "\n")
} else {
  cat("no installed forestsearch found in .libPaths(); multisession workers may fail (verdict would be 'inconclusive')\n")
}

make_df <- function(seed, n, hr_harm, two_planted) {
  set.seed(seed)
  b <- matrix(stats::rbinom(5L * n, 1L, 0.5), ncol = 5L)
  treat <- stats::rbinom(n, 1L, 0.5)
  harm <- b[, 1L] == 1L & b[, 2L] == 1L
  if (two_planted) harm <- harm | (b[, 3L] == 1L & b[, 4L] == 1L)
  rate <- 0.08 * exp(treat * ifelse(harm, log(hr_harm), log(0.8)))
  t_ev <- stats::rexp(n, rate)
  t_c <- stats::runif(n, 8, 30)
  data.frame(id = seq_len(n), tte = pmin(t_ev, t_c), event = as.integer(t_ev <= t_c),
             treat = treat, b1 = b[, 1L], b2 = b[, 2L], b3 = b[, 3L], b4 = b[, 4L],
             b5 = b[, 5L])
}

run_fs <- function(df, plan, method, pthr, splits, twostage) {
  w <- character(0)
  fs <- NULL
  junk <- utils::capture.output(
    fs <- withCallingHandlers(
      tryCatch(
        suppressMessages(forestsearch(
          df.analysis = df, outcome.name = "tte", event.name = "event",
          treat.name = "treat", id.name = "id",
          confounders.name = c("b1", "b2", "b3", "b4", "b5"),
          consistency_method = method, use_twostage = twostage,
          fs.splits = splits, pconsistency.threshold = pthr,
          n.min = 40, d0.min = 5, d1.min = 5, maxk = 2,
          details = FALSE, quiet = TRUE, seedit = 8316951,
          parallel_args = list(plan = plan, workers = 2, show_message = FALSE))),
        error = function(e) structure(list(error = conditionMessage(e)), class = "fs_error")),
      warning = function(cond) {
        w <<- c(w, conditionMessage(cond))
        invokeRestart("muffleWarning")
      }))
  list(fs = fs, warn = w)
}

failed <- function(r) {
  inherits(r$fs, "fs_error") ||
    any(grepl("Error in subgroup\\.consistency|Error in subgroup\\.search|Subgroup search failed|FSdata failure",
              r$warn))
}
sel_label <- function(r) {
  if (failed(r)) return("<FAILED>")
  if (is.null(r$fs$sg.harm)) return("<none>")
  paste(r$fs$sg.harm, collapse = " & ")
}
res_tab <- function(r) {
  if (inherits(r$fs, "fs_error")) return(NULL)
  x <- r$fs$grp.consistency$out_sg$result
  if (is.null(x)) return(NULL)
  x <- as.data.frame(x)
  rownames(x) <- NULL
  x
}
top_vals <- function(r) {
  x <- res_tab(r)
  if (is.null(x) || nrow(x) == 0L) return(NULL)
  unlist(x[1L, c("Pcons", "hr", "N", "E", "K")])
}
fmt <- function(v) {
  if (is.null(v)) return("<no row>")
  paste(sprintf("%s=%s", names(v), signif(v, 6)), collapse = " ")
}
pcons_by_cand <- function(r) {
  x <- res_tab(r)
  if (is.null(x)) return(numeric(0))
  stats::setNames(as.numeric(x$Pcons), as.character(x$m))
}
pcons_diff <- function(ra, rb) {
  pa <- pcons_by_cand(ra)
  pb <- pcons_by_cand(rb)
  ids <- union(names(pa), names(pb))
  sum(vapply(ids, function(i) !isTRUE(all.equal(unname(pa[i]), unname(pb[i]))), logical(1)))
}

compare <- function(label, df, method, pthr, splits, twostage) {
  rs <- run_fs(df, "sequential", method, pthr, splits, twostage)
  rp <- run_fs(df, "multisession", method, pthr, splits, twostage)
  fail <- failed(rs) || failed(rp)
  s_seq <- sel_label(rs)
  s_par <- sel_label(rp)
  cand <- function(r) if (inherits(r$fs, "fs_error")) NULL else r$fs$find.grps$out.found$hr.subgroups
  search_same <- isTRUE(all.equal(as.data.frame(cand(rs)), as.data.frame(cand(rp)),
                                  check.attributes = FALSE))
  est_same <- isTRUE(all.equal(top_vals(rs), top_vals(rp)))
  n_pc <- if (fail) NA_integer_ else pcons_diff(rs, rp)
  cat(sprintf("[%s] sequential  : selection = %s | top row: %s\n", label, s_seq, fmt(top_vals(rs))))
  cat(sprintf("[%s] multisession: selection = %s | top row: %s\n", label, s_par, fmt(top_vals(rp))))
  cat(sprintf("[%s] search table identical = %s | candidates whose Pcons/admission differ = %s | selection identical = %s | top-row estimates identical = %s\n",
              label, search_same, n_pc, identical(s_seq, s_par), est_same))
  if (fail) {
    msgs <- unique(c(rs$warn, rp$warn,
                     if (inherits(rs$fs, "fs_error")) rs$fs[["error"]],
                     if (inherits(rp$fs, "fs_error")) rp$fs[["error"]]))
    cat(sprintf("[%s] FAILED run(s): %s\n", label, paste(msgs, collapse = " || ")))
  }
  list(label = label, fail = fail,
       sel_differs = !fail && !identical(s_seq, s_par),
       est_differs = !fail && !est_same,
       n_pc = n_pc, search_same = search_same)
}

res <- tryCatch({
  A <- compare("A: default consistency_method='resample'",
               make_df(101L, 400L, 2.5, FALSE), "resample", 0.90, 100L, TRUE)
  B <- list()
  for (ds in 201:208) {
    if (proc.time()[["elapsed"]] - t_start > 65) break
    B[[length(B) + 1L]] <- compare(sprintf("B: consistency_method='split', data seed %d", ds),
                                   make_df(ds, 300L, 2.0, TRUE), "split", 0.60, 40L, FALSE)
    last <- B[[length(B)]]
    if (isTRUE(last$fail) || isTRUE(last$sel_differs)) break
  }
  list(A = A, B = B)
}, finally = {
  future::plan(old_plan)
  if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
})

flag <- function(lst, nm) vapply(lst, function(z) isTRUE(z[[nm]]), logical(1))
any_fail <- isTRUE(res$A$fail) || any(flag(res$B, "fail"))
b_sel <- flag(res$B, "sel_differs")
b_pc <- sum(vapply(res$B, function(z) isTRUE(z$n_pc > 0L), logical(1)))
cat("plan restored:", identical(class(future::plan()), class(old_plan)), "\n")
cat("RESAMPLE_DEFAULT_IDENTICAL (selection and top-row estimates):",
    !isTRUE(res$A$fail) && !isTRUE(res$A$sel_differs) && !isTRUE(res$A$est_differs), "\n")
cat("SPLIT comparisons run:", length(res$B),
    "| with different per-candidate Pcons/admission:", b_pc,
    "| with different selection:", sum(b_sel), "\n")
if (any_fail) {
  cat("REPRO_VERDICT: inconclusive (a forestsearch run failed; see FAILED lines)\n")
} else if (isTRUE(res$A$sel_differs) || any(b_sel)) {
  where <- c(if (isTRUE(res$A$sel_differs)) res$A$label,
             vapply(res$B[b_sel], function(z) z$label, character(1)))
  cat("REPRO_VERDICT: demonstrated (sequential and multisession select different subgroups for the same seed in:",
      paste(where, collapse = "; "), ")\n")
} else {
  cat("REPRO_VERDICT: not demonstrated (selection identical in all", 1L + length(res$B),
      "comparisons; split-path Pcons differed in", b_pc, "of", length(res$B), ")\n")
}
```

```text
installed forestsearch (the copy multisession workers load): Version 0.4.1.9002 | Packaged 2026-11-06 23:42:45 UTC; larryleon 
[A: default consistency_method='resample'] sequential  : selection = {b1} & {b2} | top row: Pcons=1 hr=2.71324 N=95 E=77 K=2
[A: default consistency_method='resample'] multisession: selection = {b1} & {b2} | top row: Pcons=1 hr=2.71324 N=95 E=77 K=2
[A: default consistency_method='resample'] search table identical = TRUE | candidates whose Pcons/admission differ = 0 | selection identical = TRUE | top-row estimates identical = TRUE
[B: consistency_method='split', data seed 201] sequential  : selection = {b2} & {b5} | top row: Pcons=1 hr=2.95317 N=75 E=61 K=2
[B: consistency_method='split', data seed 201] multisession: selection = {b2} & {b5} | top row: Pcons=1 hr=2.95317 N=75 E=61 K=2
[B: consistency_method='split', data seed 201] search table identical = TRUE | candidates whose Pcons/admission differ = 27 | selection identical = TRUE | top-row estimates identical = TRUE
[B: consistency_method='split', data seed 202] sequential  : selection = {b1} & {b2} | top row: Pcons=1 hr=2.64039 N=76 E=67 K=2
[B: consistency_method='split', data seed 202] multisession: selection = {b1} & {b2} | top row: Pcons=1 hr=2.64039 N=76 E=67 K=2
[B: consistency_method='split', data seed 202] search table identical = TRUE | candidates whose Pcons/admission differ = 12 | selection identical = TRUE | top-row estimates identical = TRUE
[B: consistency_method='split', data seed 203] sequential  : selection = {b2} & {b3} | top row: Pcons=1 hr=2.52574 N=74 E=61 K=2
[B: consistency_method='split', data seed 203] multisession: selection = {b1} & {b3} | top row: Pcons=1 hr=2.27233 N=75 E=65 K=2
[B: consistency_method='split', data seed 203] search table identical = TRUE | candidates whose Pcons/admission differ = 15 | selection identical = FALSE | top-row estimates identical = FALSE
plan restored: TRUE 
RESAMPLE_DEFAULT_IDENTICAL (selection and top-row estimates): TRUE 
SPLIT comparisons run: 3 | with different per-candidate Pcons/admission: 3 | with different selection: 1 
REPRO_VERDICT: demonstrated (sequential and multisession select different subgroups for the same seed in: B: consistency_method='split', data seed 203 )
```

### V001

- Finding: A-002 — On the default consistency path each bootstrap replicate passes forestsearch() a resample with repeated subject ids, and forestsearch() attaches the subgroup fl
- Runs: exit 0, 7.5 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
future::plan(future::sequential)
verdict <- "not demonstrated"
mk <- function(seed, N = 300L) {
  set.seed(seed)
  d <- data.frame(id = seq_len(N), treat = rbinom(N, 1L, 0.5),
                  age = round(rnorm(N, 55, 12)),
                  z1 = factor(rbinom(N, 1L, 0.5), levels = c(0L, 1L)),
                  z2 = factor(rbinom(N, 1L, 0.5), levels = c(0L, 1L)))
  lam <- 0.05 * exp(ifelse(d$z1 == "1", log(3), log(0.6)) * d$treat)
  tt <- rexp(N, lam)
  d$time <- pmin(tt, 60)
  d$event <- as.integer(tt <= 60)
  d
}
fsargs <- list(outcome.name = "time", event.name = "event", treat.name = "treat", id.name = "id",
  confounders.name = c("age", "z1", "z2"), use_grf = FALSE, use_lasso = FALSE,
  hr.threshold = 1.25, hr.consistency = 1.0, pconsistency.threshold = 0.80,
  fs.splits = 100L, n.min = 40L, d0.min = 5L, d1.min = 5L, maxk = 2L, max_subgroups_search = 5L,
  details = FALSE, quiet = TRUE, plot.sg = FALSE, plot.grf = FALSE,
  parallel_args = list(plan = "sequential", workers = 1L, show_message = FALSE))
f <- build_cox_formula("time", "event", "treat")
est <- function(d) unname(get_Cox_sg(df_sg = d[which(d$treat.recommend == 0), , drop = FALSE],
                                     cox.formula = f, est.loghr = TRUE, treat.name = "treat")$est_obs)
fs <- NULL
for (s in 1:5) {
  fs <- suppressWarnings(suppressMessages(do.call(forestsearch, c(list(df.analysis = mk(s)), fsargs))))
  if (!is.null(fs$sg.harm)) break
}
if (is.null(fs) || is.null(fs$sg.harm)) {
  cat("Original fit identified no subgroup; the bootstrap cannot be exercised.\n")
} else {
  ## Part A: one replicate built as bootstrap_results() builds it (lines 393-396, 438, 550, 558-607, 614)
  N <- nrow(fs$df.est)
  keep <- !(names(fs$df.est) %in% c(fs$confounders.candidate, "treat.recommend"))
  a <- fs$args_call_all
  a$df.predict <- fs$df.est[, keep, drop = FALSE]
  a$details <- FALSE; a$show_candidate_summary <- FALSE; a$plot.sg <- FALSE
  a$plot.grf <- FALSE; a$quiet <- TRUE
  a$grf_res <- NULL; a$grf_cuts <- NULL; a$dina_res <- NULL; a$dina_cuts <- NULL; a$ps_hat <- NULL
  a$parallel_args$plan <- "sequential"; a$parallel_args$workers <- 1L
  a$parallel_args$show_message <- FALSE
  a$mr_inference <- FALSE
  rb <- NULL; idx <- NULL
  for (b in 1:20) {
    set.seed(1000 + b)
    idx <- sample.int(N, size = N, replace = TRUE)
    a$df.analysis <- fs$df.est[idx, keep, drop = FALSE]
    rb <- suppressWarnings(suppressMessages(try(do.call(forestsearch, a), silent = TRUE)))
    if (!inherits(rb, "try-error") && !is.null(rb$sg.harm)) break
  }
  if (inherits(rb, "try-error") || is.null(rb$sg.harm)) {
    cat("No replicate identified a subgroup; nothing to compare.\n")
  } else {
    k <- tabulate(idx, nbins = N)
    n_oob <- sum(k == 0L)
    n_est <- nrow(rb$df.est); n_pred <- nrow(rb$df.predict)
    n_na <- sum(is.na(rb$df.predict$treat.recommend))
    ## correct 'H* on the bootstrap sample': the N resampled rows with the replicate's own flag
    boot_ok <- a$df.analysis
    boot_ok$treat.recommend <- rb$df.est$treat.recommend[match(boot_ok$id, rb$df.est$id)]
    code_star <- est(rb$df.est); code_obs <- est(rb$df.predict); ok_star <- est(boot_ok)
    ## correct 'H* on the original data': the definition applied to all N original subjects
    ok_obs <- tryCatch({
      o <- suppressWarnings(get_dfpred(df.predict = a$df.predict, sg.harm = rb$sg.harm, version = 2))
      same <- all(o$treat.recommend[match(boot_ok$id, o$id)] == boot_ok$treat.recommend)
      if (isTRUE(same)) est(o) else NA_real_
    }, error = function(e) NA_real_)
    cat(sprintf("N = %d; sum(k) = %d; sum(k^2) = %d; out-of-bag subjects = %d\n", N, sum(k), sum(k^2), n_oob))
    cat(sprintf("rows of replicate $df.est (frame for Hstar_star): observed %d, correct %d\n", n_est, N))
    cat(sprintf("rows of replicate $df.predict (frame for Hstar_obs): observed %d with %d NA flags, correct %d with 0 NA flags\n", n_pred, n_na, N))
    cat(sprintf("Hstar_star (log HR of H* on bootstrap data): code %.6f, correct %.6f\n", code_star, ok_star))
    cat(sprintf("Hstar_obs (log HR of H* on ORIGINAL data): code %.6f, correct %s\n", code_obs, format(ok_obs, digits = 7)))
    cat(sprintf("code's Hstar_obs equals H* on the bootstrap sample (%.6f): %s\n", ok_star, isTRUE(abs(code_obs - ok_star) < 1e-6)))
    cat(sprintf("bias term Hstar_star - Hstar_obs: code %.6f, correct %s\n", code_star - code_obs, format(ok_star - ok_obs, digits = 7)))
    demoA <- isTRUE(n_est == sum(k^2)) && n_est != N && isTRUE(n_pred == N + n_oob) && n_oob > 0L &&
      n_na == n_oob && isTRUE(abs(code_obs - ok_star) < 1e-6) && isTRUE(abs(code_star - ok_star) > 1e-6)
    if (demoA) verdict <- "demonstrated"
  }
  ## Part B: end to end. If 'original data' were the original data, events(H*) + events(Hc*) would equal
  ## the trial's event count in every replicate; the code returns the resample's event count instead.
  B <- 5L
  bc <- tryCatch(suppressWarnings(suppressMessages(forestsearch_bootstrap_dofuture(fs.est = fs, nb_boots = B,
          parallel_args = list(plan = "sequential", workers = 1L, show_message = FALSE)))),
        error = function(e) NULL)
  if (!is.null(bc)) {
    r <- as.data.frame(bc$results); ymat <- bc$Ystar_mat
  } else {
    future::plan(future::sequential)
    set.seed(8316951L)
    im <- matrix(sample.int(N, size = N * B, replace = TRUE), nrow = B, ncol = N, byrow = TRUE)
    cf <- fit_cox_models(fs$df.est, f, treat.name = "treat")
    r <- as.data.frame(suppressWarnings(suppressMessages(bootstrap_results(fs.est = fs,
           df_boot_analysis = fs$df.est, cox.formula.boot = f, nb_boots = B, show_three = FALSE,
           H_obs = cf$H_obs, Hc_obs = cf$Hc_obs, boot_index_mat = im))))
    ymat <- t(apply(im, 1, tabulate, nbins = N))[r$boot_id, , drop = FALSE]
  }
  pos <- which(r$any_found == 1 & !is.na(r$events_Hstar_0))
  if (length(pos) > 0L && nrow(ymat) == nrow(r)) {
    ev_code <- with(r[pos, ], events_Hstar_0 + events_Hstar_1 + events_Hcstar_0 + events_Hcstar_1)
    ev_boot <- as.vector(ymat[pos, , drop = FALSE] %*% fs$df.est$event)
    print(data.frame(boot_id = r$boot_id[pos], events_code_original_data = ev_code,
                     events_in_bootstrap_sample = ev_boot, correct_events_in_original_data = sum(fs$df.est$event)))
    cat("Part B (end to end) shows the same defect:",
        all(ev_code == ev_boot) && any(ev_code != sum(fs$df.est$event)), "\n")
  } else {
    cat("Part B: no replicate with an identified subgroup to tabulate.\n")
  }
}
cat(paste0("REPRO_VERDICT: ", verdict, "\n"))
```

```text
N = 300; sum(k) = 300; sum(k^2) = 578; out-of-bag subjects = 108
rows of replicate $df.est (frame for Hstar_star): observed 578, correct 300
rows of replicate $df.predict (frame for Hstar_obs): observed 408 with 108 NA flags, correct 300 with 0 NA flags
Hstar_star (log HR of H* on bootstrap data): code 0.809118, correct 0.960551
Hstar_obs (log HR of H* on ORIGINAL data): code 0.960551, correct 1.132732
code's Hstar_obs equals H* on the bootstrap sample (0.960551): TRUE
bias term Hstar_star - Hstar_obs: code -0.151432, correct -0.1721818
  boot_id events_code_original_data events_in_bootstrap_sample
1       1                       270                        270
2       2                       278                        278
3       3                       274                        274
4       4                       275                        275
5       5                       276                        276
  correct_events_in_original_data
1                             277
2                             277
3                             277
4                             277
5                             277
Part B (end to end) shows the same defect: TRUE 
REPRO_VERDICT: demonstrated
```

### V002

- Finding: A-001 — For a continuous (gaussian/MD) outcome with both threshold spellings left unset, both DINA routes take the harm floor from the raw formal hr.threshold = 1.25 in
- Runs: exit 0, 4.3 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(20261109)
n <- 400
# covariates rounded to one decimal so the DINA candidate family stays small
d <- data.frame(id = seq_len(n),
                x1 = round(rnorm(n), 1), x2 = round(rnorm(n), 1),
                treat = rbinom(n, 1, 0.5))
# adverse_outcome = TRUE below: higher y is worse.
# Harm of 0.6 outcome units (below 1.25, above 0) where x1 > 0.
d$y <- 0.6 * d$treat * (d$x1 > 0) + rnorm(n, sd = 0.5)

run <- function(args) {
  msgs <- character(0)
  err  <- NA_character_
  fs   <- NULL
  invisible(utils::capture.output(
    fs <- tryCatch(
      withCallingHandlers(
        do.call(forestsearch, args),
        message = function(m) {
          msgs <<- c(msgs, conditionMessage(m))
          invokeRestart("muffleMessage")
        },
        warning = function(w) invokeRestart("muffleWarning")),
      error = function(e) { err <<- conditionMessage(e); NULL })))
  # floor actually handed to dina_subgroup(): the details-time line
  # '  Harm floor:          m_diff = %.4f' (forestsearch_helpers.R:1548)
  ln <- grep("Harm floor", unlist(strsplit(msgs, "\n", fixed = TRUE)), value = TRUE)
  fl <- if (length(ln)) suppressWarnings(as.numeric(
          sub(".*m_diff = *([-0-9.eE+]+).*", "\\1", ln[1]))) else NA_real_
  # fallback when a subgroup was found: dina_subgroup() records its m_diff
  if (is.na(fl) && !is.null(fs$grp.consistency$out_sg$m_diff))
    fl <- as.numeric(fs$grp.consistency$out_sg$m_diff)
  list(fs = fs, floor = fl, err = err)
}

base <- list(df.analysis = d, outcome.name = "y", treat.name = "treat", id.name = "id",
             confounders.name = c("x1", "x2"), outcome_type = "continuous",
             adverse_outcome = TRUE, subgroup_method = "dina",
             details = TRUE, quiet = TRUE,
             parallel_args = list(plan = "sequential", workers = 1L, show_message = FALSE))

a <- run(base)                                   # thresholds unset: documented resolved MD c1 = 0.0
b <- run(c(base, list(effect.threshold = 0)))    # the same resolved value, supplied explicitly
r <- if (!is.null(a$fs$args_call_all)) run(a$fs$args_call_all) else
       list(fs = NULL, floor = NA_real_, err = "no args_call_all to replay")  # what a bootstrap replicate / CV fold replays

for (nm in c("a", "b", "r")) {
  e <- get(nm)$err
  if (!is.na(e)) cat("run", nm, "raised an error:", e, "\n")
}

sgl <- function(x) if (is.null(x$fs)) "<error>" else if (is.null(x$fs$sg.harm)) "<none>" else paste(x$fs$sg.harm, collapse = " & ")
exp_floor <- a$fs$threshold_config$screening
if (is.null(exp_floor)) exp_floor <- NA_real_
cat("default fit reports: threshold_config$screening =", format(exp_floor),
    "| admission$effect_floor =", format(a$fs$admission$effect_floor), "\n")
cat("default fit args_call_all: effect.threshold =", format(a$fs$args_call_all$effect.threshold),
    "| hr.threshold =", format(a$fs$args_call_all$hr.threshold), "\n")
cat(sprintf("DINA harm floor m_diff | default fit: %s | explicit effect.threshold = 0: %s | replay of default fit: %s | correct (all three): %s\n",
            format(a$floor), format(b$floor), format(r$floor), format(exp_floor)))
cat(sprintf("sg.harm                | default fit: %s | explicit effect.threshold = 0: %s | replay of default fit: %s | correct: all three identical\n",
            sgl(a), sgl(b), sgl(r)))

fit_floor_wrong   <- isTRUE(abs(a$floor - exp_floor) > 1e-8)       # fit does not use the floor it reports
replay_floor_diff <- isTRUE(abs(a$floor - r$floor) > 1e-8)         # replicate uses a different floor from its parent
selection_differs <- !is.null(a$fs) && !is.null(b$fs) && !is.null(r$fs) &&
  (!identical(a$fs$sg.harm, b$fs$sg.harm) || !identical(a$fs$sg.harm, r$fs$sg.harm))
cat("fit floor differs from reported resolved floor:", fit_floor_wrong,
    "| replay floor differs from fit floor:", replay_floor_diff,
    "| selection differs across the three runs:", selection_differs, "\n")

demonstrated <- isTRUE(exp_floor == 0) &&
  (fit_floor_wrong || replay_floor_diff || selection_differs)
cat("REPRO_VERDICT:", if (demonstrated) "demonstrated" else "not demonstrated", "\n")
```

```text
default fit reports: threshold_config$screening = 0 | admission$effect_floor = 0 
default fit args_call_all: effect.threshold = 0 | hr.threshold = 1.25 
DINA harm floor m_diff | default fit: 1.25 | explicit effect.threshold = 0: 0 | replay of default fit: 0 | correct (all three): 0
sg.harm                | default fit: <none> | explicit effect.threshold = 0: {x1 >= 0.29999999999999999} & {x2 >= 0.33999999999999775} | replay of default fit: {x1 >= 0.29999999999999999} & {x2 >= 0.33999999999999775} | correct: all three identical
fit floor differs from reported resolved floor: TRUE | replay floor differs from fit floor: TRUE | selection differs across the three runs: TRUE 
REPRO_VERDICT: demonstrated
```

### V003

- Finding: A-002 — On the default consistency path df.est and df.predict are built by merge() on the id column against a flag table that has one row per analysis row, so inside ev
- Runs: exit 0, 6.5 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
library(survival)
set.seed(101); n <- 500
d <- data.frame(id = seq_len(n), x1 = rbinom(n, 1, 0.5), x2 = rbinom(n, 1, 0.5),
                x3 = rnorm(n), treat = rbinom(n, 1, 0.5))
lhr <- ifelse(d$x1 == 1, 1.2, -0.6) * d$treat
t0 <- rexp(n, 0.10 * exp(lhr)); cc <- rexp(n, 0.03)
d$tte <- pmin(t0, cc); d$event <- as.numeric(t0 <= cc)
fs0 <- suppressWarnings(suppressMessages(forestsearch(
  df.analysis = d, outcome.name = "tte", event.name = "event", treat.name = "treat",
  id.name = "id", confounders.name = c("x1", "x2", "x3"), use_lasso = FALSE, use_grf = FALSE,
  n.min = 40, d0.min = 5, d1.min = 5, fs.splits = 100, maxk = 1, details = FALSE,
  quiet = TRUE, plot.sg = FALSE,
  parallel_args = list(plan = "sequential", workers = 1L, show_message = FALSE))))
demonstrated <- FALSE
if (is.null(fs0$sg.harm) || is.null(fs0$df.est)) {
  cat("original fit found no subgroup; nothing to resample\n")
} else {
  # Build the replicate call exactly as bootstrap_results() does
  # (bootstrap_analysis_dofuture.R:393-396, 438, 550, 558-607).
  keep <- !(names(fs0$df.est) %in% c(fs0$confounders.candidate, "treat.recommend"))
  pred <- fs0$df.est[, keep, drop = FALSE]
  N <- nrow(pred)
  # Index rows are fixed BEFORE any replicate runs, so seeding inside forestsearch()
  # cannot make the ten draws identical.
  set.seed(4242)
  IDX <- t(replicate(10, sample.int(N, N, replace = TRUE)))
  fb <- NULL; idx <- NULL; boot <- NULL
  for (b in seq_len(nrow(IDX))) {
    idx_b  <- IDX[b, ]
    boot_b <- fs0$df.est[idx_b, keep, drop = FALSE]
    a <- fs0$args_call_all
    a$df.analysis <- boot_b; a$df.predict <- pred
    a$df.test <- pred   # extra: the package's own definition-based scoring of the original rows
    a$details <- FALSE; a$show_candidate_summary <- FALSE; a$plot.sg <- FALSE
    a$plot.grf <- FALSE; a$quiet <- TRUE; a$mr_inference <- FALSE
    a$grf_res <- NULL; a$grf_cuts <- NULL; a$dina_res <- NULL; a$dina_cuts <- NULL; a$ps_hat <- NULL
    a$parallel_args$plan <- "sequential"; a$parallel_args$workers <- 1L
    a$parallel_args$show_message <- FALSE
    r <- try(suppressWarnings(suppressMessages(do.call(forestsearch, a))), silent = TRUE)
    if (!inherits(r, "try-error") && !is.null(r$sg.harm) && !is.null(r$df.est) && !is.null(r$df.predict)) {
      fb <- r; idx <- idx_b; boot <- boot_b; break
    }
  }
  if (is.null(fb)) {
    cat("no replicate found a subgroup in 10 draws\n")
  } else {
    k <- as.vector(table(idx)); n_oob <- N - length(k)
    lh <- function(dd) tryCatch(unname(coef(coxph(Surv(tte, event) ~ treat, data = dd))[1]),
                                error = function(e) NA_real_)
    Hstar_obs_code  <- lh(subset(fb$df.predict, treat.recommend == 0))  # used by the bootstrap as 'original data'
    Hstar_star_code <- lh(subset(fb$df.est, treat.recommend == 0))      # used by the bootstrap as 'bootstrap data'
    ids_H <- unique(fb$df.est$id[fb$df.est$treat.recommend == 0])
    Hstar_star_true <- lh(boot[boot$id %in% ids_H, , drop = FALSE])     # H* on the resample, multiplicity k
    Hstar_obs_true  <- if (!is.null(fb$df.test)) lh(subset(fb$df.test, treat.recommend == 0)) else NA_real_  # H* on all original rows
    scoring_ok <- NA
    if (!is.null(fb$df.test)) {
      fl  <- tapply(fb$df.est$treat.recommend, fb$df.est$id, function(z) z[1])
      chk <- fb$df.test$treat.recommend[match(as.integer(names(fl)), fb$df.test$id)]
      scoring_ok <- !anyNA(chk) && all(chk == as.vector(fl))
    }
    cat("replicate subgroup:", paste(fb$sg.harm, collapse = " & "), "\n")
    cat("nrow(df.est)     observed", nrow(fb$df.est), " correct", N, " sum(k^2) =", sum(k^2), "\n")
    cat("nrow(df.predict) observed", nrow(fb$df.predict), " correct", N, " N + out-of-bag =", N + n_oob, "\n")
    cat("NA flags in df.predict observed", sum(is.na(fb$df.predict$treat.recommend)), " correct 0 ; out-of-bag subjects =", n_oob, "\n")
    cat("max copies of one original subject in df.predict:", max(table(fb$df.predict$id)), " correct 1\n")
    cat("definition-based scoring of original rows agrees with in-bag flags:", scoring_ok, "\n")
    cat("Hstar_obs  (log HR) code:", Hstar_obs_code, " correct (all original rows):", Hstar_obs_true,
        " resample estimate:", Hstar_star_true, "\n")
    cat("Hstar_star (log HR) code:", Hstar_star_code, " correct (resample, multiplicity k):", Hstar_star_true, "\n")
    cat("bias term Hstar_star - Hstar_obs  code:", Hstar_star_code - Hstar_obs_code,
        " correct:", Hstar_star_true - Hstar_obs_true, "\n")
    # Informational end-to-end check: the package's own worker on the same index row.
    e2e <- try({
      cf <- build_cox_formula(outcome.name = "tte", event.name = "event", treat.name = "treat")
      ef <- fit_cox_models(fs0$df.est, cf, treat.name = "treat")
      br <- suppressWarnings(suppressMessages(bootstrap_results(
        fs.est = fs0, df_boot_analysis = fs0$df.est, cox.formula.boot = cf,
        nb_boots = 1L, show_three = FALSE, H_obs = ef$H_obs, Hc_obs = ef$Hc_obs,
        boot_index_mat = matrix(idx, nrow = 1L))))
      c(H_obs = unname(ef$H_obs), pkg = as.numeric(br$H_biasadj_1[1]))
    }, silent = TRUE)
    if (!inherits(e2e, "try-error") && length(e2e) == 2L && all(is.finite(e2e))) {
      cat("bootstrap_results() H_biasadj_1 for this index row: observed", e2e[["pkg"]],
          " ; H_obs - (code term) =", e2e[["H_obs"]] - (Hstar_star_code - Hstar_obs_code),
          " ; correct H_obs - (documented term) =", e2e[["H_obs"]] - (Hstar_star_true - Hstar_obs_true), "\n")
    } else {
      cat("end-to-end bootstrap_results() check not available (informational only)\n")
    }
    c1 <- nrow(fb$df.est) == sum(k^2) && nrow(fb$df.est) > N
    c2 <- n_oob > 0 && nrow(fb$df.predict) == N + n_oob &&
          sum(is.na(fb$df.predict$treat.recommend)) == n_oob
    c3 <- is.finite(Hstar_obs_code) && is.finite(Hstar_star_true) &&
          isTRUE(all.equal(Hstar_obs_code, Hstar_star_true, tolerance = 1e-6))
    c4 <- is.finite(Hstar_star_code) && is.finite(Hstar_star_true) &&
          abs(Hstar_star_code - Hstar_star_true) > 1e-6
    demonstrated <- c1 && c2 && c3 && c4
  }
}
cat(paste0("REPRO_VERDICT: ", if (demonstrated) "demonstrated" else "not demonstrated", "\n"))
```

```text
replicate subgroup: {x1} 
nrow(df.est)     observed 962  correct 500  sum(k^2) = 962 
nrow(df.predict) observed 671  correct 500  N + out-of-bag = 671 
NA flags in df.predict observed 171  correct 0 ; out-of-bag subjects = 171 
max copies of one original subject in df.predict: 5  correct 1
definition-based scoring of original rows agrees with in-bag flags: TRUE 
Hstar_obs  (log HR) code: 1.058859  correct (all original rows): 1.058029  resample estimate: 1.058859 
Hstar_star (log HR) code: 1.081682  correct (resample, multiplicity k): 1.058859 
bias term Hstar_star - Hstar_obs  code: 0.02282311  correct: 0.0008303906 
bootstrap_results() H_biasadj_1 for this index row: observed 1.035206  ; H_obs - (code term) = 1.035206  ; correct H_obs - (documented term) = 1.057198 
REPRO_VERDICT: demonstrated
```

### V004

- Finding: A-002 — On the default consistency path forestsearch() attaches treat.recommend to df.est and df.predict by a merge on the id column, so inside every full-bootstrap rep
- Runs: exit 0, 4.1 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
suppressPackageStartupMessages(library(survival))
qt <- function(expr) suppressWarnings(suppressMessages(expr))
set.seed(1); n <- 500
d <- data.frame(id = seq_len(n), treat = rbinom(n, 1, 0.5), x1 = rnorm(n), x2 = rnorm(n))
d$tte <- rexp(n, 0.1 * ifelse(d$treat == 1, ifelse(d$x1 < -0.6, 3, 0.45), 1))
d$event <- as.integer(d$tte < 30); d$tte <- pmin(d$tte, 30)
# default consistency path (subgroup_method not supplied)
fs0 <- qt(forestsearch(df.analysis = d, confounders.name = c("x1", "x2"), use_grf = FALSE, use_lasso = FALSE,
  maxk = 1, fs.splits = 200, parallel_args = list(plan = "sequential", workers = 1, show_message = FALSE),
  quiet = TRUE, details = FALSE, plot.sg = FALSE))
verdict <- "not demonstrated"
if (is.null(fs0$sg.harm)) {
  cat("no subgroup found on the original data; nothing to replay\n")
} else {
  cat("original fit: df.est rows", nrow(fs0$df.est), "of n =", n, "(unique ids, join is harmless here)\n")
  # one bootstrap replicate, built exactly as bootstrap_results() builds it
  # (bootstrap_analysis_dofuture.R:393-396, 438-439, 550, 558-560, 566-600, 614)
  base <- fs0$df.est
  keep <- !(names(base) %in% c(fs0$confounders.candidate, "treat.recommend"))
  orig <- base[, keep, drop = FALSE]
  rb <- NULL; a <- NULL
  for (s in 11:22) {
    set.seed(s); idx <- sample.int(n, n, replace = TRUE)
    df_boot <- base[idx, ]
    df_boot$id_boot <- seq_len(nrow(df_boot))
    a <- fs0$args_call_all
    a$df.analysis <- df_boot[, keep, drop = FALSE]
    a$df.predict <- orig
    a$details <- FALSE; a$quiet <- TRUE; a$plot.sg <- FALSE; a$plot.grf <- FALSE
    a$show_candidate_summary <- FALSE
    a$grf_res <- NULL; a$grf_cuts <- NULL; a$dina_res <- NULL; a$dina_cuts <- NULL; a$ps_hat <- NULL
    a$parallel_args$plan <- "sequential"; a$parallel_args$workers <- 1L; a$parallel_args$show_message <- FALSE
    a$mr_inference <- FALSE
    r <- try(qt(do.call(forestsearch, a)), silent = TRUE)
    if (!inherits(r, "try-error") && !is.null(r$sg.harm) && !is.null(r$df.est) && !is.null(r$df.predict)) { rb <- r; break }
  }
  if (is.null(rb)) {
    cat("no replicate found a subgroup\n")
  } else {
    k <- tabulate(match(a$df.analysis$id, orig$id), nbins = n)   # multiplicity of each original subject
    cf <- Surv(tte, event) ~ treat
    est <- function(df) unname(get_Cox_sg(df_sg = df, cox.formula = cf, est.loghr = TRUE, treat.name = "treat")$est_obs)
    rule_boot <- get_dfpred(a$df.analysis, rb$sg.harm, version = 2)$treat.recommend
    flag_boot <- rb$grp.consistency$df_flag$treat.recommend       # engine flag, one per resample row
    if (length(flag_boot) != nrow(a$df.analysis)) flag_boot <- rule_boot
    rule_orig <- get_dfpred(orig, rb$sg.harm, version = 2)        # the rule on the ORIGINAL data
    coded_obs  <- est(subset(rb$df.predict, treat.recommend == 0))  # Hstar_obs as bootstrap_results() computes it (:651, :688)
    coded_star <- est(subset(rb$df.est, treat.recommend == 0))      # Hstar_star as bootstrap_results() computes it (:702)
    true_obs   <- est(rule_orig[rule_orig$treat.recommend == 0, , drop = FALSE])
    true_star  <- est(a$df.analysis[flag_boot == 0, , drop = FALSE])
    n_na <- sum(is.na(rb$df.predict$treat.recommend))
    cat("subgroup selected in the replicate:", paste(rb$sg.harm, collapse = " & "), "\n")
    cat("rule vs engine flag on the resample, agreement:", mean(rule_boot == flag_boot), "\n")
    cat("df.est rows observed:", nrow(rb$df.est), "| correct:", n, "| sum(k^2):", sum(k^2), "\n")
    cat("df.predict rows observed:", nrow(rb$df.predict), "| correct:", n, "| NA flags:", n_na, "| out-of-bag subjects:", sum(k == 0), "\n")
    cat("H* size 'on original data' observed:", sum(rb$df.predict$treat.recommend == 0, na.rm = TRUE),
        "| correct (rule on original data):", sum(rule_orig$treat.recommend == 0),
        "| H* rows in the bootstrap sample:", sum(flag_boot == 0), "\n")
    cat("Hstar_obs  (log HR) observed:", coded_obs, "| correct (rule on original data):", true_obs,
        "| bootstrap-sample value:", true_star, "\n")
    cat("Hstar_obs equals the bootstrap-sample value:", isTRUE(abs(coded_obs - true_star) < 1e-6), "\n")
    cat("Hstar_star (log HR) observed:", coded_star, "| correct (bootstrap sample):", true_star, "\n")
    cat("optimism term Hstar_star - Hstar_obs observed:", coded_star - coded_obs, "| correct:", true_star - true_obs, "\n")
    ok <- nrow(rb$df.est) == sum(k^2) && nrow(rb$df.est) > n &&
      n_na == sum(k == 0) && sum(k == 0) > 0 &&
      abs(coded_obs - true_obs) > 1e-6 && abs(coded_star - true_star) > 1e-6
    if (isTRUE(ok)) verdict <- "demonstrated"
  }
}
cat(paste0("REPRO_VERDICT: ", verdict, "\n"))
```

```text
original fit: df.est rows 500 of n = 500 (unique ids, join is harmless here)
subgroup selected in the replicate: {x1 <= -1} 
rule vs engine flag on the resample, agreement: 1 
df.est rows observed: 990 | correct: 500 | sum(k^2): 990 
df.predict rows observed: 683 | correct: 500 | NA flags: 183 | out-of-bag subjects: 183 
H* size 'on original data' observed: 101 | correct (rule on original data): 91 | H* rows in the bootstrap sample: 101 
Hstar_obs  (log HR) observed: 1.179357 | correct (rule on original data): 1.360192 | bootstrap-sample value: 1.179357 
Hstar_obs equals the bootstrap-sample value: TRUE 
Hstar_star (log HR) observed: 1.218404 | correct (bootstrap sample): 1.179357 
optimism term Hstar_star - Hstar_obs observed: 0.0390461 | correct: -0.1808342 
REPRO_VERDICT: demonstrated
```

### V005

- Finding: A-003 — For a replicate under an alternative DGM in which nothing is identified, the simulation extractors store spec = \|Hc\|/N (the value of NPV for an empty identified
- Runs: exit 0, 2.5 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1)
n <- 200
df <- data.frame(id = seq_len(n), y_sim = rexp(n), event_sim = rbinom(n, 1, 0.7),
                 treat_sim = rep(0:1, n / 2),
                 flag_harm = c(rep(1L, 40), rep(0L, 160)))
ns <- asNamespace("forestsearch")
f <- get(".extract_fs_estimates_gen", envir = ns)
args0 <- list(df = df, dgm = NULL, cox_formula = NULL, cox_formula_adj = NULL,
              outcome_name = "y_sim", event_name = "event_sim", treat_name = "treat_sim",
              harm_col = "flag_harm", analysis = "FS", verbose = FALSE)
# (a) nothing identified: what .run_fs_analysis_gen() passes when forestsearch() finds nothing
none <- do.call(f, c(args0, list(fs_res = NULL, fs_full = NULL)))
# confusion-matrix values for an empty identified set (the function's own formulas, lines 1123-1129)
true_H <- df$flag_harm == 1L; hat_H <- rep(FALSE, n)
tn <- sum(!true_H & !hat_H); fp <- sum(!true_H & hat_H); fn <- sum(true_H & !hat_H)
formula_spec <- tn / (tn + fp); formula_npv <- tn / (tn + fn)
cat("nothing-found branch : spec =", none$spec, " npv =", none$npv, "\n")
cat("formula, empty set   : spec =", formula_spec, " npv =", formula_npv, "\n")
# (b) a detected replicate whose identified set is exactly the true subgroup (spec = 1)
hit <- tryCatch(
  do.call(f, c(args0, list(fs_res = NULL,
    fs_full = list(df.est = data.frame(treat.recommend = ifelse(df$flag_harm == 1L, 0L, 1L)))))),
  error = function(e) { cat("found-path call failed:", conditionMessage(e), "\n"); NULL })
if (!is.null(hit)) cat("exact detection      : spec =", hit$spec, " npv =", hit$npv, "\n")
# (c) consequence in the exported summary: 2 exact detections + 2 non-detections.
# No subject outside H is ever placed in the identified set, so mean specificity must be 1.
summ_spec <- NA_real_
if (!is.null(hit)) {
  res <- data.table::rbindlist(list(hit, hit, none, none), fill = TRUE)
  summ <- tryCatch(summarize_simulation_results(res),
                   error = function(e) { cat("summary failed:", conditionMessage(e), "\n"); NULL })
  if (!is.null(summ)) {
    summ_spec <- as.numeric(summ["spec", "FS"])
    cat("summarize_simulation_results(): mean spec =", summ_spec,
        " (correct = 1; any.H =", as.numeric(summ["any.H", "FS"]), ")\n")
  }
}
cat("observed stored spec =", none$spec, "; correct spec for an empty identified set = 1\n")
demonstrated <- isTRUE(all.equal(as.numeric(none$spec), sum(!true_H) / n)) &&
  isTRUE(all.equal(formula_spec, 1)) &&
  abs(as.numeric(none$spec) - 1) > 1e-8 &&
  isTRUE(all.equal(as.numeric(none$npv), 1)) &&
  (is.na(summ_spec) || abs(summ_spec - 1) > 1e-8)
cat(if (demonstrated) "REPRO_VERDICT: demonstrated\n" else "REPRO_VERDICT: not demonstrated\n")
```

```text
nothing-found branch : spec = 0.8  npv = 1 
formula, empty set   : spec = 1  npv = 0.8 
exact detection      : spec = 1  npv = 1 
summarize_simulation_results(): mean spec = 0.9  (correct = 1; any.H = 0.5 )
observed stored spec = 0.8 ; correct spec for an empty identified set = 1
REPRO_VERDICT: demonstrated
```

### V006

- Finding: A-006 — When bootstrap_results() is called directly with its default boot_index_mat = NULL, the resampling indices are drawn inside the foreach workers and are not retu
- Runs: exit 0, 7.8 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
future::plan(future::sequential)
NN <- 40L; B <- 5L; seed <- 8316951L
df <- data.frame(id = seq_len(NN), rowpos = seq_len(NN), treat = rep(0:1, length.out = NN),
                 y = as.numeric(seq_len(NN)), treat.recommend = 0L)
## What the documentation tells a direct caller to pair with bootstrap_results(): same df, nb_boots, seed
Y_doc <- bootstrap_ystar(df, nb_boots = B, seed = seed)
## The index matrix the wrapper would pass (bootstrap_dofuture_main.R:383-387), used as a positive control
set.seed(seed)
im <- matrix(sample.int(NN, size = NN * B, replace = TRUE), nrow = B, ncol = NN, byrow = TRUE)
## fs.est stub: the inner forestsearch() call fails at once (invalid sg_focus) and is caught by the
## try() at line 614, so each iteration only draws its resample and fits H on it (line 485).
stub <- list(confounders.candidate = character(0),
             args_call_all = list(outcome.name = "y", event.name = "y", treat.name = "treat",
                                  id.name = "id", sg_focus = "not_a_focus"))
## All rows have treat.recommend == 0, so the estimator closure sees the whole resample and records it.
rec <- new.env(); rec$rows <- list()
recorder <- function(d) {
  if (nrow(d) > 0L) {
    rec$rows[[length(rec$rows) + 1L]] <- d$rowpos
    message("IDX:", paste(d$rowpos, collapse = ","))
  }
  list(estimate = 0, se = 1)
}
run <- function(index_mat) {
  rec$rows <- list(); msgs <- character(0)
  out <- withCallingHandlers(
    suppressWarnings(bootstrap_results(fs.est = stub, df_boot_analysis = df, cox.formula.boot = NULL,
      nb_boots = B, show_three = FALSE, H_obs = 0, Hc_obs = 0, seed = seed,
      estimator_fn = recorder, boot_index_mat = index_mat)),
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  rows <- rec$rows
  if (length(rows) != B) {
    mm <- sub("\n$", "", msgs[startsWith(msgs, "IDX:")])
    rows <- lapply(strsplit(sub("^IDX:", "", mm), ","), as.integer)
  }
  if (length(rows) == B) t(vapply(rows, tabulate, integer(NN), nbins = NN)) else NULL
}
same_rows <- function(a, b) if (is.null(a)) NA_integer_ else sum(rowSums(a != b) == 0)
Y_ctl <- run(im)     # boot_index_mat supplied (the wrapper's route)
Y_null <- run(NULL)  # default route of a direct caller
cat("replicates:", B, "\n")
cat("index matrix supplied: rows of the multiplicity matrix actually used that equal bootstrap_ystar():",
    same_rows(Y_ctl, Y_doc), "(correct:", B, ")\n")
cat("boot_index_mat = NULL: rows of the multiplicity matrix actually used that equal bootstrap_ystar():",
    same_rows(Y_null, Y_doc), "(correct:", B, ")\n")
ctl_ok <- !is.null(Y_ctl) && isTRUE(same_rows(Y_ctl, Y_doc) == B)
demo <- ctl_ok && !is.null(Y_null) && isTRUE(same_rows(Y_null, Y_doc) < B)
cat(paste0("REPRO_VERDICT: ", if (demo) "demonstrated" else "not demonstrated", "\n"))
```

```text
replicates: 5 
index matrix supplied: rows of the multiplicity matrix actually used that equal bootstrap_ystar(): 5 (correct: 5 )
boot_index_mat = NULL: rows of the multiplicity matrix actually used that equal bootstrap_ystar(): 0 (correct: 5 )
REPRO_VERDICT: demonstrated
```

### V007

- Finding: A-011 — For a survival fit made with est.scale = "1/hr", each row of the bootstrap summary table FSsg_tab pairs one subgroup's label, counts, medians and unadjusted est
- Runs: exit 0, 10.7 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
future::plan(future::sequential)
# Data and search settings of tests/testthat/test-mr-method-tables.R:113-125 (without MR), est.scale added.
set.seed(1); n <- 500
d <- data.frame(id = seq_len(n), treat = rbinom(n, 1, 0.5), x1 = rnorm(n), x2 = rnorm(n))
d$time <- rexp(n, 0.1 * ifelse(d$treat == 1, ifelse(d$x1 < -0.6, 3, 0.45), 1))
d$event <- as.integer(d$time < 30); d$time <- pmin(d$time, 30)
run <- function(scale) {
  fs <- suppressWarnings(suppressMessages(forestsearch(
    df.analysis = d, outcome.name = "time", event.name = "event", treat.name = "treat", id.name = "id",
    confounders.name = c("x1", "x2"), use_grf = FALSE, use_lasso = FALSE, maxk = 1, fs.splits = 200,
    sg_focus = "hrMaxSG", est.scale = scale, parallel_args = list(plan = "sequential"),
    quiet = TRUE, details = FALSE, plot.sg = FALSE, show_candidate_summary = FALSE)))
  stopifnot(!is.null(fs$sg.harm), identical(fs$args_call_all$est.scale, scale))
  bc <- suppressWarnings(suppressMessages(forestsearch_bootstrap_dofuture(
    fs.est = fs, nb_boots = 6L, seed = 8316951L, parallel_args = list(plan = "sequential"))))
  list(fs = fs, bc = bc)
}
lead <- function(s) suppressWarnings(as.numeric(sub("^\\s*([-0-9.eE]+).*$", "\\1", s)))
report <- function(r, scale) {
  fs <- r$fs; bc <- r$bc
  cat("\n==== est.scale =", scale, "====\n")
  cat("sg.harm:", paste(fs$sg.harm, collapse = " & "), "\n")
  if (is.null(bc$FSsg_tab) || is.null(bc$SG_CIs) || is.null(bc$H_estimates) || is.null(bc$Hc_estimates)) {
    cat("no FSsg_tab / estimates returned\n"); return(NA)
  }
  m <- as.matrix(bc$FSsg_tab)
  print(m[, c("Subgroup", "n", "HR (95% CI)", "HR*")])
  dfe <- fs$df.est
  n_H <- sum(dfe$treat.recommend == 0); n_Hc <- sum(dfe$treat.recommend == 1)
  hr_of <- function(flag) {
    h <- unname(exp(stats::coef(survival::coxph(survival::Surv(time, event) ~ treat,
                                                data = dfe[dfe$treat.recommend == flag, ]))))
    if (scale == "1/hr") 1 / h else h
  }
  H_bc <- bc$SG_CIs$H_bc; Hc_bc <- bc$SG_CIs$Hc_bc
  H2 <- as.numeric(bc$H_estimates$H2); Hc2 <- as.numeric(bc$Hc_estimates$H2)
  row1_bc <- unname(m[1, "HR*"]); row2_bc <- unname(m[2, "HR*"])
  cat("row 1 label:", m[1, "Subgroup"], "| row 1 n:", lead(m[1, "n"]),
      "| n with treat.recommend == 0 (H):", n_H, "| n complement:", n_Hc, "\n")
  cat("row 1 unadjusted estimate:", lead(m[1, "HR (95% CI)"]),
      "| independent Cox on H:", round(hr_of(0), 2), "| on complement:", round(hr_of(1), 2), "\n")
  cat("SG_CIs$H_bc  (H, bias-corrected):          ", H_bc, " [H_estimates$H2 =", round(H2, 4), "]\n")
  cat("SG_CIs$Hc_bc (complement, bias-corrected): ", Hc_bc, " [Hc_estimates$H2 =", round(Hc2, 4), "]\n")
  cat("row 1 HR* observed:", row1_bc, "| correct (H's own):", H_bc, "\n")
  cat("row 2 HR* observed:", row2_bc, "| correct (complement's own):", Hc_bc, "\n")
  usable <- n_H != n_Hc && isTRUE(lead(m[1, "n"]) == n_H) && !grepl("NA", H_bc) && !grepl("NA", Hc_bc) &&
    H_bc != Hc_bc && is.finite(H2) && is.finite(Hc2) && abs(H2 - Hc2) > 1e-3
  if (!usable) { cat("run not usable for the comparison\n"); return(NA) }
  swapped <- row1_bc == Hc_bc && row2_bc == H_bc &&
    abs(lead(row1_bc) - Hc2) < 1e-3 && abs(lead(row2_bc) - H2) < 1e-3
  cat("row 1 holds H's data but the complement's bias-corrected estimate:", swapped, "\n")
  swapped
}
sw_inv <- report(run("1/hr"), "1/hr")
sw_hr <- tryCatch(report(run("hr"), "hr"), error = function(e) NA)
cat("\nswapped on 1/hr:", sw_inv, "| swapped on hr (control, expected FALSE):", sw_hr, "\n")
cat("REPRO_VERDICT:", if (isTRUE(sw_inv) && !isTRUE(sw_hr)) "demonstrated" else "not demonstrated", "\n")
```

```text
res_0 "Qstnbl" "134 (26.8%)" "0.28 (0.19, 0.41)" "1.7867 (0.1327,24.059)"
res_1 "Recmnd" "366 (73.2%)" "1.87 (1.49, 2.33)" "0.3576 (0.0818,1.5632)"
row 1 label: Qstnbl | row 1 n: 134 | n with treat.recommend == 0 (H): 134 | n complement: 366 
row 1 unadjusted estimate: 0.28 | independent Cox on H: 0.28 | on complement: 1.87 
SG_CIs$H_bc  (H, bias-corrected):           0.3576 (0.0818,1.5632)  [H_estimates$H2 = 0.3576 ]
SG_CIs$Hc_bc (complement, bias-corrected):  1.7867 (0.1327,24.059)  [Hc_estimates$H2 = 1.7867 ]
row 1 HR* observed: 1.7867 (0.1327,24.059) | correct (H's own): 0.3576 (0.0818,1.5632) 
row 2 HR* observed: 0.3576 (0.0818,1.5632) | correct (complement's own): 1.7867 (0.1327,24.059) 
row 1 holds H's data but the complement's bias-corrected estimate: TRUE 

==== est.scale = hr ====
sg.harm: {x1 <= -0.7} 
      Subgroup n             HR (95% CI)         HR*                      
res_0 "Qstnbl" "134 (26.8%)" "3.56 (2.45, 5.16)" "2.7961 (0.6397,12.2214)"
res_1 "Recmnd" "366 (73.2%)" "0.54 (0.43, 0.67)" "0.5597 (0.0416,7.5365)" 
row 1 label: Qstnbl | row 1 n: 134 | n with treat.recommend == 0 (H): 134 | n complement: 366 
row 1 unadjusted estimate: 3.56 | independent Cox on H: 3.56 | on complement: 0.54 
SG_CIs$H_bc  (H, bias-corrected):           2.7961 (0.6397,12.2214)  [H_estimates$H2 = 2.7961 ]
SG_CIs$Hc_bc (complement, bias-corrected):  0.5597 (0.0416,7.5365)  [Hc_estimates$H2 = 0.5597 ]
row 1 HR* observed: 2.7961 (0.6397,12.2214) | correct (H's own): 2.7961 (0.6397,12.2214) 
row 2 HR* observed: 0.5597 (0.0416,7.5365) | correct (complement's own): 0.5597 (0.0416,7.5365) 
row 1 holds H's data but the complement's bias-corrected estimate: FALSE 

swapped on 1/hr: TRUE | swapped on hr (control, expected FALSE): FALSE 
REPRO_VERDICT: demonstrated
```

### V008

- Finding: A-014 — dina_fit(family = "gaussian", baseline_method = list(eta0 = , eta1 = )) is documented and passes validation but always stops with 'could not find function "nu_f
- Runs: exit 0, 2.5 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1)
n <- 300
X <- matrix(stats::runif(n * 2, -1, 1), n, 2)
W <- stats::rbinom(n, 1, 0.5)
Y <- X[, 1] + W * (0.5 + X[, 1]) + stats::rnorm(n)
eta0 <- function(X) X[, 1]
eta1 <- function(X) 0.5 + 2 * X[, 1]

run <- function(...) tryCatch(dina_fit(...), error = function(e) e)
show <- function(r) if (inherits(r, "error")) paste("ERROR:", conditionMessage(r)) else
  paste(format(r$coefficients, digits = 5), collapse = " ")

# Documented configuration: Gaussian with list(eta0, eta1); with and without cross-fitting
r_cf   <- run(X, Y, W, family = "gaussian", propensity_method = 0.5,
              baseline_method = list(eta0 = eta0, eta1 = eta1), seed = 1)
r_nocf <- run(X, Y, W, family = "gaussian", propensity_method = 0.5,
              baseline_method = list(eta0 = eta0, eta1 = eta1),
              cross_fitting = FALSE, seed = 1)

# What a working implementation would return: for Gaussian a(x) = e(x), so
# nu(x) = e * eta1 + (1 - e) * eta0; supplied through the working `nu` route.
ref <- run(X, Y, W, family = "gaussian", propensity_method = 0.5,
           baseline_method = list(nu = function(X) 0.5 * eta1(X) + 0.5 * eta0(X)),
           seed = 1)

# Control (informational only): the same kind of list is accepted for binomial.
Yb <- stats::rbinom(n, 1, stats::plogis(X[, 1] + W * 0.5))
ctl <- run(X, Yb, W, family = "binomial", propensity_method = 0.5,
           baseline_method = list(eta0 = function(X) X[, 1],
                                  eta1 = function(X) 0.5 + X[, 1]), seed = 1)

cat("observed, gaussian list(eta0, eta1), cross-fitting   :", show(r_cf), "\n")
cat("observed, gaussian list(eta0, eta1), no cross-fitting:", show(r_nocf), "\n")
cat("expected, equivalent list(nu) route                  :", show(ref), "\n")
cat("control,  binomial list(eta0, eta1)                  :", show(ctl), "\n")

is_nu_err <- function(r) inherits(r, "error") && grepl("nu_fn", conditionMessage(r), fixed = TRUE)
demonstrated <- is_nu_err(r_cf) && is_nu_err(r_nocf) &&
  !inherits(ref, "error") && all(is.finite(ref$coefficients))
cat("REPRO_VERDICT:", if (demonstrated) "demonstrated" else "not demonstrated", "\n")
```

```text
observed, gaussian list(eta0, eta1), cross-fitting   : ERROR: could not find function "nu_fn" 
observed, gaussian list(eta0, eta1), no cross-fitting: ERROR: could not find function "nu_fn" 
expected, equivalent list(nu) route                  : 0.62209 0.61874 0.27792 
control,  binomial list(eta0, eta1)                  :  0.38412 -0.15872 -0.46970 
REPRO_VERDICT: demonstrated
```

### V009

- Finding: A-016 — dina() and dina_bagged() document `df` as 'data frame (or tibble / data.table)', but the shared extraction helper selects columns with `df[, required_cols]`, wh
- Runs: exit 0, 3.6 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)

msg_of <- function(x) if (inherits(x, "try-error")) conditionMessage(attr(x, "condition")) else "<no error>"
cls_of <- function(x) paste(class(x), collapse = "/")
# the data.table error for a bare-symbol j that is not a column
is_j_err <- function(x) inherits(x, "try-error") &&
  grepl("single symbol|\\.\\.required_cols|\\.\\.covariates", msg_of(x))

cat("data.table version:", as.character(utils::packageVersion("data.table")), "\n")
cat("forestsearch namespace imports data.table:",
    "data.table" %in% names(getNamespaceImports(asNamespace("forestsearch"))), "\n\n")

set.seed(1)
n <- 200
df0 <- data.frame(w = rbinom(n, 1, 0.5), x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
df0$y <- 1 + df0$x1 + df0$w * (0.5 + 0.8 * df0$x1) + rnorm(n)
dt0 <- data.table::as.data.table(df0)
covs <- c("x1", "x2")

## ---- Part 1: dina() / dina_bagged(), the documented input classes ----
r_df <- try(dina(df = df0, outcome = "y", treatment = "w", covariates = covs,
                 family = "gaussian", seed = 1L), silent = TRUE)
r_dt <- try(dina(df = dt0, outcome = "y", treatment = "w", covariates = covs,
                 family = "gaussian", seed = 1L), silent = TRUE)
r_wk <- try(dina(df = as.data.frame(dt0), outcome = "y", treatment = "w", covariates = covs,
                 family = "gaussian", seed = 1L), silent = TRUE)
r_bag <- try(dina_bagged(df = dt0, outcome = "y", treatment = "w", covariates = covs,
                         family = "gaussian", n_bags = 5L, seed = 1L), silent = TRUE)
cat("dina(data.frame)                -> ", cls_of(r_df), "\n")
cat("dina(data.table)                -> ", cls_of(r_dt), "\n   message: ", msg_of(r_dt), "\n")
cat("dina(as.data.frame(data.table)) -> ", cls_of(r_wk), "\n")
cat("dina_bagged(data.table)         -> ", cls_of(r_bag), "\n   message: ", msg_of(r_bag), "\n")
if (requireNamespace("tibble", quietly = TRUE)) {
  r_tb <- try(dina(df = tibble::as_tibble(df0), outcome = "y", treatment = "w", covariates = covs,
                   family = "gaussian", seed = 1L), silent = TRUE)
  cat("dina(tibble)                    -> ", cls_of(r_tb), " (expected to work)\n")
}

## ---- Part 2 (informational): sibling helper .dina_extract_X_from_df ----
if (inherits(r_df, "dina")) {
  f_df <- try(dina_frontier(fit = r_df, df = df0, covariates = covs, n_min = 20L), silent = TRUE)
  f_dt <- try(dina_frontier(fit = r_df, df = dt0, covariates = covs, n_min = 20L), silent = TRUE)
  s_df <- try(dina_subgroup(fit = r_df, df = df0, covariates = covs, m_diff = 0, n_min = 20L), silent = TRUE)
  s_dt <- try(dina_subgroup(fit = r_df, df = dt0, covariates = covs, m_diff = 0, n_min = 20L), silent = TRUE)
  cat("\ndina_frontier(data.frame) -> ", cls_of(f_df), "\n")
  cat("dina_frontier(data.table) -> ", cls_of(f_dt), " | data.table j error: ", is_j_err(f_dt), "\n")
  cat("dina_subgroup(data.frame) -> ", cls_of(s_df), "\n")
  cat("dina_subgroup(data.table) -> ", cls_of(s_dt), " | data.table j error: ", is_j_err(s_dt), "\n")
}

## ---- Part 3 (informational): the pipeline entry point, subgroup_method = "dina" ----
## forestsearch() documents df.analysis as 'Data frame' only; this shows whether a
## data.table reaches dina() through forestsearch_helpers.R:1534. Not part of the verdict.
set.seed(7)
N <- 160
sv <- data.frame(id = seq_len(N), treat = rbinom(N, 1L, 0.5),
                 age = round(rnorm(N, 55, 12)), bm = rnorm(N), sexn = rbinom(N, 1L, 0.5))
tt <- rexp(N, rate = 0.05)
sv$time  <- pmin(tt, 60)
sv$event <- as.integer(tt <= 60)
fs_args <- list(
  outcome_type = "survival", outcome.name = "time", event.name = "event",
  treat.name = "treat", id.name = "id", confounders.name = c("age", "bm", "sexn"),
  subgroup_method = "dina", sg_focus = "maxSG", hr.threshold = 10, hr.consistency = 1,
  pconsistency.threshold = 0.80, max_subgroups_search = 5L, use_grf = TRUE, use_lasso = FALSE,
  is.RCT = TRUE, cont.cutoff = 4L, cut_type = "default", seedit = 42L, maxk = 2L,
  n.min = 40L, d0.min = 5L, d1.min = 5L, fs.splits = 100L, details = FALSE, quiet = TRUE,
  plot.sg = FALSE, plot.grf = FALSE,
  parallel_args = list(plan = "sequential", workers = 1L, show_message = FALSE))
fs_df <- try(suppressWarnings(do.call(forestsearch, c(list(df.analysis = sv), fs_args))), silent = TRUE)
fs_dt <- try(suppressWarnings(do.call(forestsearch,
             c(list(df.analysis = data.table::as.data.table(sv)), fs_args))), silent = TRUE)
cat("\nforestsearch(data.frame, subgroup_method = 'dina') -> ", cls_of(fs_df),
    if (inherits(fs_df, "try-error")) paste0("\n   message: ", msg_of(fs_df)) else "", "\n")
cat("forestsearch(data.table, subgroup_method = 'dina') -> ", cls_of(fs_dt),
    "\n   message: ", msg_of(fs_dt), "\n")
cat("   pipeline reaches the defect (control ran, data.table stopped with the j error): ",
    !inherits(fs_df, "try-error") && is_j_err(fs_dt), "\n")

## ---- Verdict (Part 1 only) ----
cat("\nCorrect implementation: dina() on the data.table returns a 'dina' object with the same\n",
    "coefficients as on the data.frame (the roxygen names data.table as an accepted input).\n", sep = "")
df_ok   <- inherits(r_df, "dina")
dt_jerr <- is_j_err(r_dt)
wk_same <- df_ok && inherits(r_wk, "dina") &&
  isTRUE(all.equal(unname(coef(r_wk)), unname(coef(r_df))))
dt_same <- df_ok && inherits(r_dt, "dina") &&
  isTRUE(all.equal(unname(coef(r_dt)), unname(coef(r_df))))
cat("Observed: data.frame fit ok = ", df_ok,
    " | data.table call stopped with the data.table j-symbol error = ", dt_jerr,
    " | same data as plain data.frame gives identical coefficients = ", wk_same,
    " | data.table fit equals data.frame fit = ", dt_same, "\n", sep = "")
cat(if (df_ok && dt_jerr && wk_same) "REPRO_VERDICT: demonstrated\n" else "REPRO_VERDICT: not demonstrated\n")
```

```text
data.table version: 1.18.6.1 
forestsearch namespace imports data.table: TRUE 

dina(data.frame)                ->  dina 
dina(data.table)                ->  try-error 
   message:  j (the 2nd argument inside [...]) is a single symbol but column name 'required_cols' is not found. If you intended to select columns using a variable in calling scope, please try DT[, ..required_cols]. The .. prefix conveys one-level-up similar to a file system path. 
dina(as.data.frame(data.table)) ->  dina 
dina_bagged(data.table)         ->  try-error 
   message:  j (the 2nd argument inside [...]) is a single symbol but column name 'required_cols' is not found. If you intended to select columns using a variable in calling scope, please try DT[, ..required_cols]. The .. prefix conveys one-level-up similar to a file system path. 
dina(tibble)                    ->  dina  (expected to work)

dina_frontier(data.frame) ->  data.frame 
dina_frontier(data.table) ->  try-error  | data.table j error:  TRUE 
dina_subgroup(data.frame) ->  dina_subgroup 
dina_subgroup(data.table) ->  try-error  | data.table j error:  TRUE 

forestsearch(data.frame, subgroup_method = 'dina') ->  forestsearch/list  
forestsearch(data.table, subgroup_method = 'dina') ->  try-error 
   message:  j (the 2nd argument inside [...]) is a single symbol but column name 'required_cols' is not found. If you intended to select columns using a variable in calling scope, please try DT[, ..required_cols]. The .. prefix conveys one-level-up similar to a file system path. 
   pipeline reaches the defect (control ran, data.table stopped with the j error):  TRUE 

Correct implementation: dina() on the data.table returns a 'dina' object with the same
coefficients as on the data.frame (the roxygen names data.table as an accepted input).
Observed: data.frame fit ok = TRUE | data.table call stopped with the data.table j-symbol error = TRUE | same data as plain data.frame gives identical coefficients = TRUE | data.table fit equals data.frame fit = FALSE
REPRO_VERDICT: demonstrated
```

### V010

- Finding: A-001 — With outcome_type = "continuous" (MD), subgroup_method = "dina" and neither threshold supplied, DINA's candidacy floor m_diff is the raw survival default hr.thr
- Runs: exit 0, 7.9 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)

# Continuous outcome, true harm of 0.8 units in x1 > 0 (above 0, below 1.25).
set.seed(11)
n <- 400
d <- data.frame(id = seq_len(n), treat = rbinom(n, 1, 0.5),
                x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
d$y <- 0.5 * d$x2 + d$treat * 0.8 * (d$x1 > 0) + rnorm(n, sd = 0.5)

# Run forestsearch() and read the floor DINA really applied from the
# 'Harm floor: m_diff = ...' details message of .forestsearch_dina_select().
run <- function(args) {
  msgs <- character(0)
  fit <- NULL
  invisible(utils::capture.output(
    fit <- withCallingHandlers(
      do.call(forestsearch, args),
      message = function(m) {
        msgs <<- c(msgs, conditionMessage(m))
        invokeRestart("muffleMessage")
      },
      warning = function(w) invokeRestart("muffleWarning"))))
  txt <- paste(msgs, collapse = "\n")
  hit <- regmatches(txt, regexpr("m_diff = -?[0-9.]+", txt))
  list(fit = fit,
       m_diff = if (length(hit)) as.numeric(sub("m_diff = ", "", hit)) else NA_real_)
}
lab <- function(fit) if (is.null(fit$sg.harm)) "<none>" else paste(fit$sg.harm, collapse = " & ")
num <- function(x) if (is.null(x) || !length(x)) NA_real_ else as.numeric(x)[1L]

base <- list(df.analysis = d, outcome.name = "y", treat.name = "treat", id.name = "id",
             confounders.name = c("x1", "x2"), outcome_type = "continuous",
             adverse_outcome = TRUE, subgroup_method = "dina", n.min = 60,
             mr_inference = FALSE, details = TRUE, quiet = TRUE,
             parallel_args = list(plan = "sequential", workers = 1L, show_message = FALSE))

r_def  <- run(base)                                 # thresholds left at their defaults
r_zero <- run(c(base, list(effect.threshold = 0)))  # the documented resolved MD default, made explicit
# What a bootstrap replicate runs: the default fit's own args_call_all (here on the same data).
replay_args <- r_def$fit$args_call_all
replay_args$df.analysis <- d
r_replay <- run(replay_args)

scr <- num(r_def$fit$threshold_config$screening)
flo <- num(r_def$fit$admission$effect_floor)
rep_thr <- num(r_def$fit$args_call_all$effect.threshold)
cat("Default call (continuous / MD, subgroup_method = 'dina', no threshold supplied):\n")
cat("  reported threshold_config$screening       :", format(scr), "  (documented MD default: 0)\n")
cat("  reported admission$effect_floor           :", format(flo), "\n")
cat("  args_call_all$effect.threshold (replayed) :", format(rep_thr), "\n")
cat("  DINA m_diff actually applied              :", format(r_def$m_diff), "  (correct: equal to the screening value above)\n")
cat("  selected subgroup                         :", lab(r_def$fit), "\n")
cat("Explicit effect.threshold = 0:\n")
cat("  DINA m_diff actually applied              :", format(r_zero$m_diff), "\n")
cat("  selected subgroup                         :", lab(r_zero$fit), "  (correct: same as the default call)\n")
cat("Replay of the default fit's own args_call_all:\n")
cat("  DINA m_diff actually applied              :", format(r_replay$m_diff), "  (correct: same as the original fit)\n")
cat("  selected subgroup                         :", lab(r_replay$fit), "  (correct: same as the original fit)\n")

floor_differs_from_reported <- is.finite(r_def$m_diff) && is.finite(scr) &&
  abs(r_def$m_diff - scr) > 1e-8
replay_differs_from_parent <- is.finite(r_def$m_diff) && is.finite(r_replay$m_diff) &&
  abs(r_def$m_diff - r_replay$m_diff) > 1e-8
cat("floor differs from reported screening threshold:", floor_differs_from_reported, "\n")
cat("replay floor differs from parent floor         :", replay_differs_from_parent, "\n")
cat("selection differs (default vs explicit 0)      :", !identical(lab(r_def$fit), lab(r_zero$fit)), "\n")
cat("selection differs (default vs replay)          :", !identical(lab(r_def$fit), lab(r_replay$fit)), "\n")
cat("REPRO_VERDICT:",
    if (floor_differs_from_reported && replay_differs_from_parent) "demonstrated" else "not demonstrated",
    "\n")
```

```text
Default call (continuous / MD, subgroup_method = 'dina', no threshold supplied):
  reported threshold_config$screening       : 0   (documented MD default: 0)
  reported admission$effect_floor           : 0 
  args_call_all$effect.threshold (replayed) : 0 
  DINA m_diff actually applied              : 1.25   (correct: equal to the screening value above)
  selected subgroup                         : <none> 
Explicit effect.threshold = 0:
  DINA m_diff actually applied              : 0 
  selected subgroup                         : {x1 >= 0.4160062456969173} & {x2 <= 0.25437322240322824}   (correct: same as the default call)
Replay of the default fit's own args_call_all:
  DINA m_diff actually applied              : 0   (correct: same as the original fit)
  selected subgroup                         : {x1 >= 0.4160062456969173} & {x2 <= 0.25437322240322824}   (correct: same as the original fit)
floor differs from reported screening threshold: TRUE 
replay floor differs from parent floor         : TRUE 
selection differs (default vs explicit 0)      : TRUE 
selection differs (default vs replay)          : TRUE 
REPRO_VERDICT: demonstrated
```

### V011

- Finding: C-001 — forestsearch() documents est.scale as "hr" or "rmst", never reads or validates it, and no code implements "rmst" (every consumer of the stored value knows only 
- Runs: exit 0, 2.4 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
pkg <- "/home/larryleon/Documents/GitHub/forestsearch"

## 1. The values documented for forestsearch(est.scale = )
rd <- readLines(file.path(pkg, "man", "forestsearch.Rd"))
doc_line <- grep("\\item{est.scale}", rd, fixed = TRUE, value = TRUE)
cat("Rd text:", doc_line, "\n")
doc_says_rmst <- length(doc_line) == 1L && grepl("rmst", doc_line, fixed = TRUE)
doc_says_inv  <- length(doc_line) == 1L && grepl("1/hr", doc_line, fixed = TRUE)
cat("doc mentions 'rmst':", doc_says_rmst, "  (correct doc: FALSE)\n")
cat("doc mentions '1/hr':", doc_says_inv,  "  (correct doc: TRUE)\n")

## 2. forestsearch() neither reads nor validates est.scale
b <- body(forestsearch)
body_reads <- ("est.scale" %in% all.names(b)) ||
  any(grepl("est.scale", deparse(b), fixed = TRUE))
cat("formal default:", deparse(formals(forestsearch)$est.scale), "\n")
cat("forestsearch() body references est.scale:", body_reads,
    "  (a validating implementation: TRUE)\n")

## 3. The bootstrap consumer of args_call_all$est.scale: get_dfRes(), called inside
##    try() at bootstrap_dofuture_main.R:574-609, with the documented value
set.seed(1)
B <- 8L; N <- 20L
ystar <- t(stats::rmultinom(B, N, rep(1 / N, N)))   # B x N multiplicities
h_adj <- stats::rnorm(B, 0.3, 0.1)
call_dfRes <- function(scale) try(
  get_dfRes(Hobs = log(1.5), seHobs = 0.3, H1_adj = h_adj, H2_adj = h_adj,
            ystar = ystar, est.scale = scale, est.loghr = TRUE), silent = TRUE)
r_rmst <- call_dfRes("rmst"); r_hr <- call_dfRes("hr"); r_inv <- call_dfRes("1/hr")
show <- function(r) if (inherits(r, "try-error")) trimws(as.character(r)) else "returned estimates"
cat("get_dfRes(est.scale = 'rmst'):", show(r_rmst), "  (if 'rmst' were supported: returned estimates)\n")
cat("get_dfRes(est.scale = 'hr')  :", show(r_hr), "\n")
cat("get_dfRes(est.scale = '1/hr'):", show(r_inv), "\n")
rmst_fails <- inherits(r_rmst, "try-error") && grepl("est.scale", as.character(r_rmst), fixed = TRUE)
hr_ok <- !inherits(r_hr, "try-error") && !inherits(r_inv, "try-error")

## 4. The table consumer (summary_utility_functions.R:52) treats 'rmst' as 'hr' with no message
d <- data.frame(treat = c(0, 1, 0, 1, 0, 1), sg_flag = c(1, 1, 0, 0, 1, 0))
p_rmst <- prepare_subgroup_data(d, SG_flag = "sg_flag", est.scale = "rmst", treat.name = "treat")
p_hr   <- prepare_subgroup_data(d, SG_flag = "sg_flag", est.scale = "hr",   treat.name = "treat")
tables_silent <- identical(p_rmst, p_hr)
cat("prepare_subgroup_data('rmst') identical to ('hr'):", tables_silent,
    "  (a validating implementation: error)\n")

demonstrated <- isTRUE(doc_says_rmst) && !isTRUE(body_reads) && isTRUE(rmst_fails) &&
  isTRUE(hr_ok) && isTRUE(tables_silent)
cat("REPRO_VERDICT:", if (demonstrated) "demonstrated" else "not demonstrated", "\n")
```

```text
Rd text: \item{est.scale}{Character. Estimation scale ("hr" or "rmst"). Default "hr".} 
doc mentions 'rmst': TRUE   (correct doc: FALSE)
doc mentions '1/hr': FALSE   (correct doc: TRUE)
formal default: "hr" 
forestsearch() body references est.scale: FALSE   (a validating implementation: TRUE)
get_dfRes(est.scale = 'rmst'): Error in get_dfRes(Hobs = log(1.5), seHobs = 0.3, H1_adj = h_adj, H2_adj = h_adj,  : 
  'est.scale' must be either 'hr' or '1/hr'.   (if 'rmst' were supported: returned estimates)
get_dfRes(est.scale = 'hr')  : returned estimates 
get_dfRes(est.scale = '1/hr'): returned estimates 
prepare_subgroup_data('rmst') identical to ('hr'): TRUE   (a validating implementation: error)
REPRO_VERDICT: demonstrated
```

### V012

- Finding: A-030 — When ps_method resolves to anything other than "none" on a GLM outcome, the estimator closure is rebuilt without adverse_outcome and reverts to make_effect_esti
- Runs: exit 0, 4.4 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
future::plan(future::sequential)
set.seed(42); N <- 400L
df <- data.frame(id = seq_len(N), treat = rbinom(N, 1L, 0.5),
                 age = round(rnorm(N, 55, 12)), biomarker = rnorm(N),
                 biomarker_hi = rbinom(N, 1L, 0.5), sex = rbinom(N, 1L, 0.5))
# Higher y = better (continuous default: adverse_outcome resolves to FALSE).
# Treatment LOWERS y by 1.5 in biomarker_hi == 1 (true harm) and RAISES y by 1.0 in biomarker_hi == 0 (benefit).
df$y <- rnorm(N, mean = df$treat * ifelse(df$biomarker_hi == 1, -1.5, 1.0), sd = 1)
harm <- df[df$biomarker_hi == 1, ]
ben  <- df[df$biomarker_hi == 0, ]
raw_md <- function(d) mean(d$y[d$treat == 1]) - mean(d$y[d$treat == 0])

# --- 1. Closure level: the two constructions forestsearch() makes (forestsearch_main.R:2220 and :2705)
e_first   <- make_effect_estimator(outcome_type = "continuous", treat.name = "treat", outcome.name = "y",
                                   effect_measure = "MD", adverse_outcome = FALSE)
e_rebuilt <- make_effect_estimator(outcome_type = "continuous", treat.name = "treat", outcome.name = "y",
                                   offset.name = NULL, effect_measure = "MD", ps_adjust_method = "none")
est_first   <- e_first(harm)$estimate
est_rebuilt <- e_rebuilt(harm)$estimate
cat("true-harm subgroup, harm-scale effect: first build =", round(est_first, 4),
    "| PS rebuild =", round(est_rebuilt, 4), "| correct rebuild =", round(-raw_md(harm), 4), "\n")

# --- 2. End to end. ps_adjust_method is left at its default ("none"): no adjustment is requested in any run.
run <- function(rct, ps, cm = "resample") {
  res <- NULL
  invisible(utils::capture.output(res <- tryCatch(suppressWarnings(suppressMessages(forestsearch(
    df.analysis = df, confounders.name = c("age", "biomarker", "biomarker_hi", "sex"),
    outcome.name = "y", event.name = "y", treat.name = "treat", id.name = "id",
    outcome_type = "continuous", effect_measure = "MD",
    effect.threshold = 0.5, consistency.threshold = 0.0, pconsistency.threshold = 0.5,
    fs.splits = 20, n.min = 30, maxk = 1, use_lasso = FALSE, use_grf = FALSE,
    is.RCT = rct, ps_method = ps, consistency_method = cm,
    details = FALSE, quiet = TRUE, seedit = 1L, parallel_args = list(plan = "sequential")))),
    error = function(e) list(sg.harm = paste("ERROR:", conditionMessage(e)), args_call_all = list()))))
  res
}
tab  <- function(fs) { t <- fs$find.grps$out.found$hr.subgroups; if (is.null(t)) NULL else as.data.frame(t) }
top  <- function(fs) { t <- tab(fs); if (is.null(t) || nrow(t) == 0L) NA_real_ else max(t$HR) }
topn <- function(fs) { t <- tab(fs); if (is.null(t) || nrow(t) == 0L) NA_real_ else t$n[which.max(t$HR)] }
sg   <- function(fs) if (is.null(fs$sg.harm)) "NULL" else paste(fs$sg.harm, collapse = " & ")
a  <- run(TRUE,  NULL)                 # RCT default: ps_method resolves to "none", branch at :2656 skipped
b  <- run(FALSE, "logistic")           # explicit ps_method (mirrors tests/testthat/test-iptw-continuous.R)
g  <- run(FALSE, NULL)                 # is.RCT = FALSE only: ps_method resolves to "grf" (:2598)
b2 <- run(FALSE, "logistic", "split")  # split consistency uses the same closure as the search
topA <- top(a); topB <- top(b); topG <- top(g)
cat("resolved adverse_outcome: RCT =", a$args_call_all$adverse_outcome, "| logistic =", b$args_call_all$adverse_outcome,
    "| is.RCT=FALSE only =", g$args_call_all$adverse_outcome, "\n")
cat("resolved ps_method / ps_adjust_method: RCT =", a$args_call_all$ps_method, "/", a$args_call_all$ps_adjust_method,
    "| logistic =", b$args_call_all$ps_method, "/", b$args_call_all$ps_adjust_method,
    "| is.RCT=FALSE only =", g$args_call_all$ps_method, "/", g$args_call_all$ps_adjust_method, "\n")
cat("hand-computed: harm-scale effect in biomarker_hi == 1 (-raw MD) =", round(-raw_md(harm), 4), "(n =", nrow(harm), ")",
    "| RAW MD in biomarker_hi == 0, a benefit =", round(raw_md(ben), 4), "(n =", nrow(ben), ")\n")
cat("top screening effect (n of that candidate): RCT =", round(topA, 4), "(", topn(a), ")",
    "| logistic =", round(topB, 4), "(", topn(b), ")",
    "| is.RCT=FALSE only =", round(topG, 4), "(", topn(g), ")",
    "| correct for all three (no adjustment requested) =", round(topA, 4), "(", topn(a), ")\n")
cat("candidates: RCT =", NROW(tab(a)), "| logistic =", NROW(tab(b)), "| is.RCT=FALSE only =", NROW(tab(g)), "\n")
cat("sg.harm: RCT =", sg(a), "| logistic (resample) =", sg(b), "| is.RCT=FALSE only (resample) =", sg(g),
    "| logistic (split) =", sg(b2), "| correct = same as RCT run\n")

closure_flip   <- is.finite(est_first) && is.finite(est_rebuilt) && est_first > 0 && est_rebuilt < 0 &&
                  isTRUE(all.equal(est_rebuilt, -est_first))
same_config    <- isFALSE(a$args_call_all$adverse_outcome) && isFALSE(b$args_call_all$adverse_outcome) &&
                  identical(b$args_call_all$ps_adjust_method, "none")
screen_differs <- is.finite(topA) && !isTRUE(all.equal(topA, topB))
default_trigger <- isFALSE(g$args_call_all$adverse_outcome) && identical(g$args_call_all$ps_method, "grf") &&
                   is.finite(topA) && !isTRUE(all.equal(topA, topG))
cat("closure_flip =", closure_flip, "| same_config =", same_config, "| screen_differs =", screen_differs,
    "| selection_differs =", !identical(a$sg.harm, b$sg.harm),
    "| reached_by_is.RCT_FALSE_alone =", default_trigger, "\n")
cat("REPRO_VERDICT:", if (closure_flip && same_config && screen_differs) "demonstrated" else "not demonstrated", "\n")
```

```text
true-harm subgroup, harm-scale effect: first build = 1.6117 | PS rebuild = -1.6117 | correct rebuild = 1.6117 
resolved adverse_outcome: RCT = FALSE | logistic = FALSE | is.RCT=FALSE only = FALSE 
resolved ps_method / ps_adjust_method: RCT = none / none | logistic = logistic / none | is.RCT=FALSE only = grf / none 
hand-computed: harm-scale effect in biomarker_hi == 1 (-raw MD) = 1.6117 (n = 210 ) | RAW MD in biomarker_hi == 0, a benefit = 0.7296 (n = 190 )
top screening effect (n of that candidate): RCT = 1.6117 ( 210 ) | logistic = 0.7296 ( 190 ) | is.RCT=FALSE only = 0.7296 ( 190 ) | correct for all three (no adjustment requested) = 1.6117 ( 210 )
candidates: RCT = 3 | logistic = 1 | is.RCT=FALSE only = 1 
sg.harm: RCT = {biomarker_hi} | logistic (resample) = NULL | is.RCT=FALSE only (resample) = NULL | logistic (split) = !{biomarker_hi} | correct = same as RCT run
closure_flip = TRUE | same_config = TRUE | screen_differs = TRUE | selection_differs = TRUE | reached_by_is.RCT_FALSE_alone = TRUE 
REPRO_VERDICT: demonstrated
```

### V013

- Finding: A-036 — Under sg_focus = "maxeff" (documented as the effect argmax with no effect or consistency threshold, and recorded with an unrestricted admission set) the consist
- Runs: exit 0, 4.1 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
future::plan(future::sequential)
set.seed(20261109)
n <- 600
d <- data.frame(id = seq_len(n),
                x1 = rnorm(n), x2 = rnorm(n), x3 = rbinom(n, 1, 0.5),
                treat = rep(0:1, n / 2))
# uniform strong benefit: true HR = 0.25 in every subgroup, so every candidate HR is well below 1
t0 <- rexp(n, rate = 0.10 * exp(log(0.25) * d$treat))
cc <- runif(n, 5, 40)
d$tte <- pmin(t0, cc)
d$event <- as.integer(t0 <= cc)

run <- function(...) suppressWarnings(suppressMessages(forestsearch(
  df.analysis = d, outcome.name = "tte", event.name = "event",
  treat.name = "treat", id.name = "id",
  confounders.name = c("x1", "x2", "x3"), sg_focus = "maxeff",
  use_lasso = FALSE, use_grf = FALSE, n.min = 60, maxk = 1, fs.splits = 50,
  details = FALSE, quiet = TRUE, plot.sg = FALSE,
  parallel_args = list(plan = "sequential", workers = 1L, show_message = FALSE), ...)))
lab <- function(x) if (is.null(x$sg.harm)) "NULL (no subgroup)" else paste(x$sg.harm, collapse = " & ")

# Run a: maxeff with default thresholds (hr.threshold = 1.25, hr.consistency = c2 = 1.0)
a  <- run()
hs <- a$find.grps$out.found$hr.subgroups
n_cand <- if (is.null(hs)) 0L else nrow(hs)
max_hr <- if (n_cand > 0) max(hs$HR, na.rm = TRUE) else NA_real_
expected <- NA_character_
if (n_cand > 0) {
  hsd <- as.data.frame(hs)
  zc  <- setdiff(names(hsd), c("grp", "K", "n", "E", "d1", "m1", "m0", "HR", "L(HR)", "U(HR)"))
  top <- hsd[which.max(hsd$HR), , drop = FALSE]
  sel <- zc[unlist(top[zc]) == 1]
  # map cut-column names to their labels when the fit carries the mapping
  cand <- a$confounders.candidate; ev <- a$confounders.evaluated
  if (!is.null(cand) && !is.null(ev) && length(cand) == length(ev) && all(sel %in% cand)) {
    sel <- ev[match(sel, cand)]
  }
  expected <- paste(sel, collapse = " & ")
}
cat("candidates in enumerated family      :", n_cand, "\n")
cat("largest candidate HR                 :", round(max_hr, 4), "\n")
cat("recorded admission set               :", forestsearch:::.fs_format_admission(a$admission), "\n")
cat("documented maxeff result (argmax cut):", expected, " (a non-NULL subgroup)\n")
cat("observed, c2 = 1.0 (default)         :", lab(a), "\n")

# Run b: same data, same focus; only c2 lowered below the largest candidate HR
b <- NULL
if (is.finite(max_hr) && max_hr > 0) {
  c2_low <- signif(0.9 * max_hr, 3)
  b <- tryCatch(run(hr.consistency = c2_low), error = function(e) NULL)
  cat("observed, c2 =", c2_low, "            :", if (is.null(b)) "error" else lab(b), "\n")
}

# Run c (informational): only hr.threshold changed -- the maxeff warning calls it disabled,
# but c2 is derived as 0.80 * hr.threshold, which moves the gate
if (is.finite(max_hr) && max_hr > 0) {
  c1_low <- signif(1.125 * max_hr, 3)
  cc_run <- tryCatch(run(hr.threshold = c1_low), error = function(e) e)
  cat("observed, hr.threshold =", c1_low, "only :",
      if (inherits(cc_run, "error")) paste("error:", conditionMessage(cc_run)) else lab(cc_run), "\n")
}

demonstrated <- n_cand > 0 && is.finite(max_hr) && max_hr <= 1 &&
  is.null(a$sg.harm) &&
  is.null(a$admission$effect_floor) && is.null(a$admission$consistency) &&
  !is.null(b) && !is.null(b$sg.harm)
cat("REPRO_VERDICT:", if (isTRUE(demonstrated)) "demonstrated" else "not demonstrated", "\n")
```

```text
candidates in enumerated family      : 14 
largest candidate HR                 : 0.3624 
recorded admission set               : unrestricted (no effect floor, no consistency floor) 
documented maxeff result (argmax cut): q7.0  (a non-NULL subgroup)
observed, c2 = 1.0 (default)         : NULL (no subgroup) 
observed, c2 = 0.326             : !{x3} 
observed, hr.threshold = 0.408 only : !{x3} 
REPRO_VERDICT: demonstrated
```

### V014

- Finding: A-038 — merge() is called with its default sort = TRUE, so on the consistency path df.est and df.predict come back ordered by the id key while grp.consistency$sg.harm.i
- Runs: exit 0, 2.6 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(101); n <- 500
d <- data.frame(id = seq_len(n), x1 = rbinom(n, 1, 0.5), x2 = rbinom(n, 1, 0.5),
                x3 = rnorm(n), treat = rbinom(n, 1, 0.5))
lhr <- ifelse(d$x1 == 1, 1.2, -0.6) * d$treat
t0 <- rexp(n, 0.10 * exp(lhr)); cc <- rexp(n, 0.03)
d$tte <- pmin(t0, cc); d$event <- as.numeric(t0 <= cc)
d <- d[sample.int(n), ]; rownames(d) <- NULL          # rows not in id order
fs <- suppressWarnings(suppressMessages(forestsearch(
  df.analysis = d, outcome.name = "tte", event.name = "event", treat.name = "treat",
  id.name = "id", confounders.name = c("x1", "x2", "x3"), use_lasso = FALSE, use_grf = FALSE,
  n.min = 40, d0.min = 5, d1.min = 5, fs.splits = 100, maxk = 1, details = FALSE,
  quiet = TRUE, plot.sg = FALSE,
  parallel_args = list(plan = "sequential", workers = 1L, show_message = FALSE))))
demonstrated <- FALSE
if (is.null(fs$sg.harm) || is.null(fs$df.est)) {
  cat("no subgroup found; nothing to compare\n")
} else {
  own <- as.integer(fs$grp.consistency$sg.harm.id)            # input row order
  est <- as.integer(fs$df.est$treat.recommend == 0)           # df.est row order
  same_order <- identical(as.vector(fs$df.est$id), as.vector(d$id))
  mism <- if (length(own) == length(est)) sum(own != est) else NA_integer_
  est_re <- est[match(d$id, fs$df.est$id)]                    # df.est put back in input order
  mism_re <- sum(own != est_re)
  cat("df.est in input row order: observed", same_order, " correct TRUE\n")
  cat("df.est sorted by id:", !is.unsorted(fs$df.est$id), "\n")
  cat("rows where sg.harm.id[i] disagrees with df.est$treat.recommend[i] == 0: observed", mism,
      "of", n, " correct 0\n")
  cat("same comparison after re-aligning df.est by id:", mism_re, "(0 shows the flags are right and only the order differs)\n")
  demonstrated <- !same_order && isTRUE(mism > 0) && isTRUE(mism_re == 0)
}
cat(paste0("REPRO_VERDICT: ", if (demonstrated) "demonstrated" else "not demonstrated", "\n"))
```

```text
df.est in input row order: observed FALSE  correct TRUE
df.est sorted by id: TRUE 
rows where sg.harm.id[i] disagrees with df.est$treat.recommend[i] == 0: observed 246 of 500  correct 0
same comparison after re-aligning df.est by id: 0 (0 shows the flags are right and only the order differs)
REPRO_VERDICT: demonstrated
```

### V015

- Finding: A-043 — A factor confounder with at least cont.cutoff (default 4) observed values and numeric-parsable levels makes get_FSdata() stop under the default collapse_cuts = 
- Runs: exit 0, 2.7 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1)
n <- 200
stage_num <- sample(1:5, n, TRUE)
df_fac <- data.frame(id = seq_len(n), treat = rbinom(n, 1, 0.5),
                     tte = rexp(n), event = rbinom(n, 1, 0.7),
                     stage = factor(stage_num), z = rnorm(n))
df_num <- df_fac
df_num$stage <- stage_num

as_err <- function(e) structure(paste("ERROR:", conditionMessage(e)), class = "repro_error")
is_err <- function(x) inherits(x, "repro_error")
run <- function(d, ...) {
  tryCatch(
    suppressWarnings(get_FSdata(df.analysis = d, confounders.name = c("stage", "z"),
                                outcome.name = "tte", event.name = "event",
                                details = FALSE, ...)$confs),
    error = as_err)
}

r_fac_default <- run(df_fac)                         # formal default collapse_cuts = TRUE
r_num_default <- run(df_num)                         # same values, numeric column
r_fac_nocoll  <- run(df_fac, collapse_cuts = FALSE)  # collapse switched off
r_fac_median  <- run(df_fac, cut_type = "median", collapse_cuts = FALSE)
r_direct <- tryCatch(collapse_redundant_cuts(cuts = "stage <= 3", df = df_fac),
                     error = as_err)

cat("factor stage, defaults (observed):\n"); print(unclass(r_fac_default))
cat("numeric stage, defaults (what a correct implementation would also give for the factor):\n"); print(unclass(r_num_default))
cat("factor stage, collapse_cuts = FALSE (observed; runs):\n"); print(unclass(r_fac_nocoll))
cat("factor stage, cut_type = 'median', collapse_cuts = FALSE (observed):\n"); print(unclass(r_fac_median))
cat("direct collapse_redundant_cuts('stage <= 3', df_fac) (observed):\n"); print(unclass(r_direct))

# Informational only (not part of the verdict): what the exported entry point returns.
wmsg <- character(0)
fs <- tryCatch(
  withCallingHandlers(
    forestsearch(df.analysis = df_fac, outcome.name = "tte", event.name = "event",
                 treat.name = "treat", id.name = "id",
                 confounders.name = c("stage", "z"),
                 parallel_args = list(plan = "sequential", workers = 1L, show_message = FALSE),
                 details = FALSE, quiet = TRUE),
    warning = function(w) { wmsg <<- c(wmsg, conditionMessage(w)); invokeRestart("muffleWarning") }),
  error = as_err)
cat("forestsearch() on the factor data, all defaults (observed):\n")
if (is_err(fs)) print(unclass(fs)) else {
  cat("  sg.harm is NULL:", is.null(fs$sg.harm), "\n")
  cat("  error_log stage:", if (is.null(fs$error_log$stage)) "<none>" else fs$error_log$stage, "\n")
}
cat("  warnings:\n"); print(wmsg)

factor_msg <- is_err(r_fac_default) && grepl("factor", r_fac_default, fixed = TRUE)
demonstrated <- factor_msg && !is_err(r_num_default) &&
  !is_err(r_fac_nocoll) && is_err(r_direct)
cat("default-path error names a factor:", factor_msg, "\n")
cat("median path also fails on the factor:", is_err(r_fac_median), "\n")
cat("REPRO_VERDICT:", if (demonstrated) "demonstrated" else "not demonstrated", "\n")
```

```text
factor stage, defaults (observed):
[1] "ERROR: Calling var(x) on a factor x is defunct.\n  Use something like 'all(duplicated(x)[-1L])' to test for a constant vector."
numeric stage, defaults (what a correct implementation would also give for the factor):
[1] "stage <= 2.9" "stage <= 3"   "stage <= 2"   "stage <= 4"   "z <= -0"     
[6] "z <= -0.7"    "z <= 0.6"    
factor stage, collapse_cuts = FALSE (observed; runs):
[1] "stage <= 2.9" "stage <= 3"   "stage <= 2"   "stage <= 4"   "z <= -0.1"   
[6] "z <= -0.7"    "z <= 0.6"    
factor stage, cut_type = 'median', collapse_cuts = FALSE (observed):
[1] "ERROR: need numeric data"
direct collapse_redundant_cuts('stage <= 3', df_fac) (observed):
[1] "ERROR: Calling var(x) on a factor x is defunct.\n  Use something like 'all(duplicated(x)[-1L])' to test for a constant vector."
forestsearch() on the factor data, all defaults (observed):
  sg.harm is NULL: TRUE 
  error_log stage: get_FSdata 
  warnings:
[1] "Error in get_FSdata: Calling var(x) on a factor x is defunct.\n  Use something like 'all(duplicated(x)[-1L])' to test for a constant vector."
[2] "FSdata failure - returning NULL result"                                                                                                      
default-path error names a factor: TRUE 
median path also fails on the factor: TRUE 
REPRO_VERDICT: demonstrated
```

### V016

- Finding: A-046 — Every confounder with fewer than cont.cutoff (default 4) unique values enters the candidate pool as its bare name, which is evaluated as `x == 1` with no check 
- Runs: exit 0, 2.3 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1)
n <- 300
sex01 <- rbinom(n, 1, 0.5)
df <- data.frame(
  tte   = rexp(n),
  event = rbinom(n, 1, 0.7),
  z     = rnorm(n, 50, 10),
  sex01 = sex01,                                   # reference coding 0/1
  sexFM = factor(ifelse(sex01 == 1, "M", "F")),    # the same variable as a factor F/M
  sex02 = 2 * sex01,                               # the same variable coded 0/2
  grade = sample(1:3, n, replace = TRUE)           # three levels coded 1/2/3 (as survival::gbsg$grade)
)
stopifnot(length(unique(df$sex01)) == 2L, length(unique(df$grade)) == 3L)
msgs <- character(0)
printed <- utils::capture.output(
  o <- withCallingHandlers(
    get_FSdata(df.analysis = df,
               confounders.name = c("z", "sex01", "sexFM", "sex02", "grade"),
               outcome.name = "tte", event.name = "event",
               use_lasso = FALSE, use_grf = FALSE,
               details = FALSE),   # the value forestsearch() hard-codes (forestsearch_main.R:3318)
    warning = function(w) { msgs <<- c(msgs, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") }
  )
)
cat("Candidate pool returned by get_FSdata():\n"); print(unname(o$confs))
vars <- c("sex01", "sexFM", "sex02", "grade")
represented <- vapply(vars, function(v) any(grepl(v, o$confs, fixed = TRUE)), logical(1))
cat("\nVariable represented in the candidate pool (observed):\n"); print(represented)
cat("Correct implementation: sexFM and sex02 carry exactly the information of sex01, so TRUE for all four (or an error/warning naming sexFM and sex02)\n")
cand <- lapply(o$df[, o$confs_names, drop = FALSE], function(x) as.numeric(as.character(x)))
level_has_candidate <- vapply(1:3, function(k) {
  any(vapply(cand, function(x) all(x == as.numeric(df$grade == k)) ||
                               all(x == as.numeric(df$grade != k)), logical(1)))
}, logical(1))
names(level_has_candidate) <- paste0("grade==", 1:3)
cat("\nLevel of 'grade' that has its own candidate indicator (observed):\n"); print(level_has_candidate)
cat("Correct implementation (one indicator per level, as vignettes/articles/methodology.qmd Step 1 states): TRUE TRUE TRUE\n")
cat("Label under which the grade candidate is reported (observed):", o$confs[grepl("grade", o$confs, fixed = TRUE)], " -- it means grade == 1\n")
silent <- length(printed) == 0L && length(msgs) == 0L
cat("\nPrinted lines:", length(printed), " warnings/messages:", length(msgs), " silent about the drop:", silent, "\n")
demonstrated <- isTRUE(represented[["sex01"]]) && !represented[["sexFM"]] && !represented[["sex02"]] &&
  identical(unname(level_has_candidate), c(TRUE, FALSE, FALSE)) && silent
cat("REPRO_VERDICT:", if (demonstrated) "demonstrated" else "not demonstrated", "\n")
```

```text
Candidate pool returned by get_FSdata():
[1] "z <= 50" "z <= 44" "z <= 57" "sex01"   "grade"  

Variable represented in the candidate pool (observed):
sex01 sexFM sex02 grade 
 TRUE FALSE FALSE  TRUE 
Correct implementation: sexFM and sex02 carry exactly the information of sex01, so TRUE for all four (or an error/warning naming sexFM and sex02)

Level of 'grade' that has its own candidate indicator (observed):
grade==1 grade==2 grade==3 
    TRUE    FALSE    FALSE 
Correct implementation (one indicator per level, as vignettes/articles/methodology.qmd Step 1 states): TRUE TRUE TRUE
Label under which the grade candidate is reported (observed): grade  -- it means grade == 1

Printed lines: 0  warnings/messages: 0  silent about the drop: TRUE 
REPRO_VERDICT: demonstrated
```

### V017

- Finding: A-055 — With effect_measure = "RD" and no covariates, a slice with an arm at 0% or 100% events gets its standard error from the logistic tier as (logit treatment-coeffi
- Runs: exit 0, 2.4 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
fn <- make_effect_estimator(outcome_type = "binary", treat.name = "treat",
                            outcome.name = "y", effect_measure = "RD",
                            adverse_outcome = TRUE)
n <- 20L
mk <- function(e0, e1) data.frame(
  treat = rep(0:1, each = n),
  y = c(rep(1:0, c(e0, n - e0)), rep(1:0, c(e1, n - e1))))
# row 1: interior control case; rows 2-5: one arm at 0% or 100% events
cfgs <- list(c(10, 14), c(10, 20), c(10, 0), c(20, 10), c(0, 10))
one <- function(cfg) {
  r <- tryCatch(suppressWarnings(fn(mk(cfg[1], cfg[2]))),
                error = function(e) list(estimate = NA_real_, se = NA_real_,
                  method_used = paste("ERROR:", conditionMessage(e)),
                  converged = NA))
  p0 <- cfg[1] / n; p1 <- cfg[2] / n
  data.frame(p0 = p0, p1 = p1, est = r$estimate, est_expected = p1 - p0,
             se = r$se,
             se_expected = sqrt(p0 * (1 - p0) / n + p1 * (1 - p1) / n),
             method = r$method_used, converged = r$converged,
             stringsAsFactors = FALSE)
}
out <- do.call(rbind, lapply(cfgs, one))
out$se_ratio <- out$se / out$se_expected
print(out, digits = 4)
# the same SE as it reaches the candidate table (Wald limits in the search)
g <- suppressWarnings(fit_glm_for_subgroup(mk(10, 20), rep(1L, 2L * n), fn))
cat(sprintf("candidate-table limits, p0=0.5 p1=1: est=%.3f lower=%.6f upper=%.6f | expected about %.3f to %.3f\n",
            g$hr, g$lower, g$upper,
            0.5 - 1.96 * sqrt(0.25 / n), 0.5 + 1.96 * sqrt(0.25 / n)))
# the same SE as it reaches the reporting path (se_obs -> unadjusted CI)
fe <- suppressWarnings(fit_subgroup_effect(mk(20, 10), fn))
cat(sprintf("reporting path, p0=1 p1=0.5: est_obs=%.3f se_obs=%.4g | expected se about %.3f\n",
            fe$est_obs, fe$se_obs, sqrt(0.25 / n)))
deg <- out[-1, ]
bad <- deg$method %in% "logistic_margins" & deg$converged %in% TRUE &
  abs(deg$est - deg$est_expected) < 1e-6 &
  (deg$se_ratio < 0.01 | deg$se_ratio > 100)
control_ok <- isTRUE(abs(out$se_ratio[1] - 1) < 0.05) &&
  identical(out$method[1], "identity_binomial")
cat("degenerate rows with correct estimate, converged = TRUE and SE off by more than 100x:",
    sum(bad), "of", nrow(deg), "| interior control row SE ratio:",
    format(out$se_ratio[1], digits = 4), "method:", out$method[1], "\n")
cat("REPRO_VERDICT:", if (any(bad) && control_ok) "demonstrated" else "not demonstrated", "\n")
```

```text
   p0  p1  est est_expected        se se_expected            method converged
1 0.5 0.7  0.2          0.2 1.517e-01      0.1517 identity_binomial      TRUE
2 0.5 1.0  0.5          0.5 7.649e-06      0.1118  logistic_margins      TRUE
3 0.5 0.0 -0.5         -0.5 7.649e-06      0.1118  logistic_margins      TRUE
4 1.0 0.5 -0.5         -0.5 6.012e+02      0.1118  logistic_margins      TRUE
5 0.0 0.5  0.5          0.5 6.012e+02      0.1118  logistic_margins      TRUE
   se_ratio
1 1.000e+00
2 6.842e-05
3 6.842e-05
4 5.377e+03
5 5.377e+03
candidate-table limits, p0=0.5 p1=1: est=0.500 lower=0.499985 upper=0.500015 | expected about 0.281 to 0.719
reporting path, p0=1 p1=0.5: est_obs=-0.500 se_obs=601.2 | expected se about 0.112
degenerate rows with correct estimate, converged = TRUE and SE off by more than 100x: 4 of 4 | interior control row SE ratio: 1 method: identity_binomial 
REPRO_VERDICT: demonstrated
```

### V018

- Finding: A-056 — With the non-default effect_measure = "RR", the log-binomial attempt in .estimate_rr() handles warnings in an exiting tryCatch() handler that calls invokeRestar
- Runs: exit 0, 2.4 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
fn <- make_effect_estimator(outcome_type = "binary", treat.name = "treat",
                            outcome.name = "y", effect_measure = "RR")
# same estimator for a beneficial (response) outcome: Y is flipped to 1 - Y inside
fn_ben <- make_effect_estimator(outcome_type = "binary", treat.name = "treat",
                                outcome.name = "y", effect_measure = "RR",
                                adverse_outcome = FALSE)
run <- function(d, f) tryCatch({
  r <- f(d)
  list(ok = TRUE, estimate = r$estimate, method = r$method_used, msg = "")
}, error = function(e) list(ok = FALSE, estimate = NA_real_,
                            method = NA_character_, msg = conditionMessage(e)))
# (a) control 10/20 events, treated 20/20 events: RR = 2, events in both arms
da <- data.frame(treat = rep(0:1, each = 20),
                 y = c(rep(1:0, each = 10), rep(1L, 20)))
# (b) no boundary proportion: control 12/40, treated 8/10: RR = 2.667
db <- data.frame(treat = rep(0:1, c(40, 10)),
                 y = c(rep(1:0, c(12, 28)), rep(1:0, c(8, 2))))
# (c) control case that does not overshoot: control 10/20, treated 14/20: RR = 1.4
dc <- data.frame(treat = rep(0:1, each = 20),
                 y = c(rep(1:0, each = 10), rep(1:0, c(14, 6))))
# (d) response outcome, adverse_outcome = FALSE: control 6/20 responders, treated 0/20;
#     on the failure scale control 14/20, treated 20/20: RR = 1/0.7
dd <- data.frame(treat = rep(0:1, each = 20),
                 y = c(rep(1:0, c(6, 14)), rep(0L, 20)))
cases <- list(a = list(d = da, f = fn,     flip = FALSE, correct = log(2)),
              b = list(d = db, f = fn,     flip = FALSE, correct = log(0.8 / 0.3)),
              c = list(d = dc, f = fn,     flip = FALSE, correct = log(1.4)),
              d = list(d = dd, f = fn_ben, flip = TRUE,  correct = log(1 / 0.7)))
# what the intended modified-Poisson fallback returns on the same slice
pois <- function(d, flip) {
  if (flip) d$y <- 1L - d$y
  unname(coef(glm(y ~ treat, family = poisson(link = "log"), data = d))["treat"])
}
res <- lapply(cases, function(k) run(k$d, k$f))
for (nm in names(res)) {
  r <- res[[nm]]; k <- cases[[nm]]
  cat(sprintf("case %s: ok=%s observed log-RR=%s method=%s | correct log-RR=%.6f (modified Poisson gives %.6f)%s\n",
              nm, r$ok, format(r$estimate, digits = 7), r$method,
              k$correct, pois(k$d, k$flip),
              if (nzchar(r$msg)) paste0(" | ERROR: ", r$msg) else ""))
}
# what the search loop does with case (a): the error is swallowed, candidate dropped
sg <- tryCatch(fit_glm_for_subgroup(da, rep(1L, nrow(da)), fn),
               error = function(e) paste("ERROR:", conditionMessage(e)))
cat("fit_glm_for_subgroup on case (a) returns NULL (candidate dropped, status 5):", is.null(sg),
    "| a correct implementation returns a fit with RR = 2\n")
hit <- function(r) !r$ok && grepl("muffleWarning", r$msg, fixed = TRUE)
demonstrated <- hit(res$a) && hit(res$b) &&
  res$c$ok && isTRUE(abs(res$c$estimate - cases$c$correct) < 1e-4) &&
  is.null(sg)
cat("case d (response outcome, zero responders on treatment) also errors:", hit(res$d), "\n")
cat("REPRO_VERDICT:", if (demonstrated) "demonstrated" else "not demonstrated", "\n")
```

```text
case a: ok=FALSE observed log-RR=NA method=NA | correct log-RR=0.693147 (modified Poisson gives 0.693147) | ERROR: no 'restart' 'muffleWarning' found
case b: ok=FALSE observed log-RR=NA method=NA | correct log-RR=0.980829 (modified Poisson gives 0.980829) | ERROR: no 'restart' 'muffleWarning' found
case c: ok=TRUE observed log-RR=0.3364722 method=log_binomial | correct log-RR=0.336472 (modified Poisson gives 0.336472)
case d: ok=FALSE observed log-RR=NA method=NA | correct log-RR=0.356675 (modified Poisson gives 0.356675) | ERROR: no 'restart' 'muffleWarning' found
fit_glm_for_subgroup on case (a) returns NULL (candidate dropped, status 5): TRUE | a correct implementation returns a fit with RR = 2
case d (response outcome, zero responders on treatment) also errors: TRUE 
REPRO_VERDICT: demonstrated
```

### V019

- Finding: A-060 — When validate_grf_data() rejects the data, grf.subg.harm.survival() calls create_null_result() with trees = list(), and trees[[1]] on an empty list throws "subs
- Runs: exit 0, 3.3 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1)
run <- function(expr) tryCatch(suppressWarnings(expr), error = function(e) e)
run_w <- function(expr) {
  w <- character(0)
  v <- tryCatch(withCallingHandlers(expr, warning = function(c) {
    w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning")
  }), error = function(e) e)
  list(value = v, warnings = w)
}
describe <- function(r) if (inherits(r, "error")) paste0("ERROR: ", conditionMessage(r)) else paste0("returned ", class(r)[1], " with sg.harm.id = ", deparse(r$sg.harm.id))
is_oob <- function(r) inherits(r, "error") && grepl("subscript out of bounds", conditionMessage(r))

# Case 1: no censoring (event column has one unique value), arms 100/100
n <- 200
d1 <- data.frame(id = 1:n, y = rexp(n), e = 1L, tr = rep(0:1, each = n / 2),
                 x1 = rnorm(n), x2 = rnorm(n))
r1 <- run(grf.subg.harm.survival(data = d1, confounders.name = c("x1", "x2"),
                                 outcome.name = "y", event.name = "e",
                                 id.name = "id", treat.name = "tr"))

# Case 2: censored data, arms 50/50, default n.min = 60 per arm
n <- 100
d2 <- data.frame(id = 1:n, y = rexp(n), e = rep(c(1L, 1L, 0L, 1L), n / 4), tr = rep(0:1, each = n / 2),
                 x1 = rnorm(n), x2 = rnorm(n))
r2 <- run(grf.subg.harm.survival(data = d2, confounders.name = c("x1", "x2"),
                                 outcome.name = "y", event.name = "e",
                                 id.name = "id", treat.name = "tr"))

# What the same helper returns when the tree list is NULL (the convention of the call at grf_main.R:297)
cfg <- create_grf_config(frac.tau = 1, n.min = 60, dmin.grf = 0, RCT = TRUE,
                         sg.criterion = "mDiff", maxdepth = 2, seedit = 1)
expected <- run(create_null_result(d2, NULL, NULL, cfg))

cat("Case 1 (no censoring)      observed:", describe(r1), "\n")
cat("Case 2 (arms 50/50 < 60)   observed:", describe(r2), "\n")
cat("Correct implementation     expected:", describe(expected), "(a null-result list, no error)\n")

# Informational (not part of the verdict): the same input through forestsearch()
seq_args <- list(plan = "sequential", workers = 1, show_message = FALSE)
r3 <- run(forestsearch(df.analysis = d2, outcome.name = "y", event.name = "e",
                       treat.name = "tr", id.name = "id",
                       confounders.name = c("x1", "x2"),
                       subgroup_method = "grf", quiet = TRUE,
                       parallel_args = seq_args))
cat("forestsearch(subgroup_method = 'grf') on case 2 observed:",
    if (inherits(r3, "error")) paste0("ERROR: ", conditionMessage(r3)) else paste0("returned object, sg.harm = ", deparse(r3$sg.harm)),
    "| expected: a result with sg.harm = NULL | aborts with out-of-bounds:", is_oob(r3), "\n")
r4 <- run_w(forestsearch(df.analysis = d2, outcome.name = "y", event.name = "e",
                         treat.name = "tr", id.name = "id",
                         confounders.name = c("x1", "x2"),
                         use_grf = TRUE, maxk = 1, quiet = TRUE,
                         parallel_args = seq_args))
cat("forestsearch(use_grf = TRUE) on case 2: run",
    if (inherits(r4$value, "error")) paste0("ERRORED: ", conditionMessage(r4$value)) else "completed",
    "| warning 'GRF analysis failed: subscript out of bounds' raised:",
    any(grepl("GRF analysis failed: subscript out of bounds", r4$warnings, fixed = TRUE)),
    "| expected: only the 'Insufficient sample size' warning\n")

demonstrated <- is_oob(r1) && is_oob(r2) && !inherits(expected, "error") && is.null(expected$sg.harm.id)
cat("REPRO_VERDICT:", if (demonstrated) "demonstrated" else "not demonstrated", "\n")
```

```text
Case 1 (no censoring)      observed: ERROR: subscript out of bounds 
Case 2 (arms 50/50 < 60)   observed: ERROR: subscript out of bounds 
Correct implementation     expected: returned list with sg.harm.id = NULL (a null-result list, no error)
forestsearch(subgroup_method = 'grf') on case 2 observed: ERROR: subscript out of bounds | expected: a result with sg.harm = NULL | aborts with out-of-bounds: TRUE 
forestsearch(use_grf = TRUE) on case 2: run completed | warning 'GRF analysis failed: subscript out of bounds' raised: TRUE | expected: only the 'Insufficient sample size' warning
REPRO_VERDICT: demonstrated
```

### V020

- Finding: A-064 — For outcome_type = "count" the adverse_outcome orientation is never applied to the GRF outcome, so the GRF 'harm' subgroup is the region where treatment LOWERS 
- Runs: exit 0, 4.5 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
future::plan("sequential")

irr_in <- function(dat, flag) {
  if (sum(flag) < 10L) return(NA_real_)
  fit <- stats::glm(y ~ treat + offset(log(pt)), family = stats::poisson(),
                    data = dat[flag, , drop = FALSE])
  unname(exp(stats::coef(fit)["treat"]))
}
par_seq <- list(plan = "sequential", workers = 1, show_message = FALSE)

# ---------------- Part A (verdict): tree path, direction reversal ----------------
set.seed(101)
n <- 800
z <- rbinom(n, 1, 0.5); x <- rnorm(n); treat <- rbinom(n, 1, 0.5); pt <- runif(n, 0.8, 1.2)
# Adverse count: treatment RAISES the rate in z == 1 (true harm) and lowers it in z == 0 (benefit)
y <- rpois(n, exp(0.5 + 0.9 * treat * (z == 1) - 0.9 * treat * (z == 0)) * pt)
d <- data.frame(id = seq_len(n), treat = treat, z = z, x = x, pt = pt, y = y)
d$ylog <- log(d$y + 0.5)

# Count path under test (adverse_outcome = TRUE is what forestsearch() passes for count)
cnt <- grf.subg.harm.glm(data = d, confounders.name = c("z", "x"), outcome.name = "y",
                         treat.name = "treat", id.name = "id", outcome_type = "count",
                         offset.name = "pt", adverse_outcome = TRUE, n.min = 60L,
                         maxdepth = 1L, seedit = 42L)
# Reference: the same transformed outcome through the continuous branch, where the
# adverse_outcome flip IS applied (-log(Y + 0.5)); what a correctly oriented count branch gives
ref <- grf.subg.harm.glm(data = d, confounders.name = c("z", "x"), outcome.name = "ylog",
                         treat.name = "treat", id.name = "id", outcome_type = "continuous",
                         adverse_outcome = TRUE, n.min = 60L, maxdepth = 1L, seedit = 42L)

flag_cnt <- cnt$data$treat.recommend == 0L
flag_ref <- ref$data$treat.recommend == 0L
irr_cnt <- irr_in(d, flag_cnt); irr_ref <- irr_in(d, flag_ref)
cat("A. count path : sg.harm.id =", paste(cnt$sg.harm.id, collapse = " & "), "| n flagged =", sum(flag_cnt),
    "| share z==1 among flagged =", round(mean(d$z[flag_cnt] == 1), 3),
    "| IRR(treat) in flagged 'harm' subgroup =", round(irr_cnt, 3), "\n")
cat("A. reference  : sg.harm.id =", paste(ref$sg.harm.id, collapse = " & "), "| n flagged =", sum(flag_ref),
    "| share z==1 among flagged =", round(mean(d$z[flag_ref] == 1), 3),
    "| IRR(treat) in flagged 'harm' subgroup =", round(irr_ref, 3), "\n")
cat("A. Correct behaviour: with adverse_outcome = TRUE the flagged subgroup is z == 1 with IRR > 1 (about exp(0.9) = 2.46).\n")
print(table(count_flagged_harm = flag_cnt, z = d$z))
print(table(reference_flagged_harm = flag_ref, z = d$z))

# Supplementary (not used for the verdict): same thing through forestsearch(), GRF selection, tree mode
fsA <- tryCatch(
  forestsearch(df.analysis = d, outcome.name = "y", event.name = "y", treat.name = "treat",
               id.name = "id", confounders.name = c("z", "x"), outcome_type = "count",
               effect_measure = "IRR", offset.name = "pt", subgroup_method = "grf",
               grf_selection = "tree", grf_select_statistic = "dr", grf_depth = 1, n.min = 60,
               details = FALSE, quiet = TRUE, plot.sg = FALSE, parallel_args = par_seq),
  error = function(e) { cat("A. forestsearch() tree call errored:", conditionMessage(e), "\n"); NULL })
if (!is.null(fsA)) {
  cat("A. forestsearch(count, subgroup_method = 'grf', grf_selection = 'tree') sg.harm =",
      paste(fsA$sg.harm, collapse = " & "), "\n")
  if (!is.null(fsA$df.est)) {
    fl <- fsA$df.est$treat.recommend == 0L
    cat("   share z==1 in that sg.harm =", round(mean(d$z[fl] == 1), 3),
        "| IRR(treat) in it =", round(irr_in(d, fl), 3), " (correct: z == 1, IRR > 1)\n")
  }
}

# ------- Part B (supplementary): default frontier + effect, treatment raises the count everywhere -------
set.seed(202)
y2 <- rpois(n, exp(0.3 + 0.8 * treat + 0.6 * treat * (z == 1)) * pt)
d2 <- data.frame(id = seq_len(n), treat = treat, z = z, x = x, pt = pt, y = y2)
d2$ylog <- log(d2$y + 0.5)
cntB <- grf.subg.harm.glm(data = d2, confounders.name = c("z", "x"), outcome.name = "y",
                          treat.name = "treat", id.name = "id", outcome_type = "count",
                          offset.name = "pt", adverse_outcome = TRUE, n.min = 60L, maxdepth = 2L,
                          seedit = 42L, grf_selection = "frontier", frontier_rule = "eff")
refB <- grf.subg.harm.glm(data = d2, confounders.name = c("z", "x"), outcome.name = "ylog",
                          treat.name = "treat", id.name = "id", outcome_type = "continuous",
                          adverse_outcome = TRUE, n.min = 60L, maxdepth = 2L,
                          seedit = 42L, grf_selection = "frontier", frontier_rule = "eff")
cat("B. True IRR(treat): z==0 about", round(exp(0.8), 2), ", z==1 about", round(exp(1.4), 2),
    "; observed IRR in z==1 =", round(irr_in(d2, d2$z == 1), 3), "\n")
cat("B. count frontier  : sg.harm.id =", if (is.null(cntB$sg.harm.id)) "NULL (no subgroup)" else paste(cntB$sg.harm.id, collapse = " & "),
    "| $candidates present =", !is.null(cntB$candidates), "\n")
cat("B. reference       : sg.harm.id =", if (is.null(refB$sg.harm.id)) "NULL (no subgroup)" else paste(refB$sg.harm.id, collapse = " & "),
    "| $candidates present =", !is.null(refB$candidates), "\n")
fsB <- tryCatch(
  forestsearch(df.analysis = d2, outcome.name = "y", event.name = "y", treat.name = "treat",
               id.name = "id", confounders.name = c("z", "x"), outcome_type = "count",
               effect_measure = "IRR", offset.name = "pt", subgroup_method = "grf", n.min = 60,
               details = FALSE, quiet = TRUE, plot.sg = FALSE, parallel_args = par_seq),
  error = function(e) { cat("B. forestsearch() default-GRF call errored:", conditionMessage(e), "\n"); NULL })
if (!is.null(fsB)) {
  cat("B. forestsearch(count, subgroup_method = 'grf') [default frontier + effect] sg.harm =",
      if (is.null(fsB$sg.harm)) "NULL (no subgroup)" else paste(fsB$sg.harm, collapse = " & "), "\n")
}
cat("B. Correct behaviour: a harm subgroup is reported (every candidate has IRR well above the 1.25 floor); per the trace the count path reports none.\n")

demonstrated <- isTRUE(irr_cnt < 1) && isTRUE(irr_ref > 1)
cat("REPRO_VERDICT:", if (demonstrated) "demonstrated" else "not demonstrated", "\n")
```

```text
A. count path : sg.harm.id = z <= 0 | n flagged = 410 | share z==1 among flagged = 0 | IRR(treat) in flagged 'harm' subgroup = 0.445 
A. reference  : sg.harm.id = z > 0 | n flagged = 390 | share z==1 among flagged = 1 | IRR(treat) in flagged 'harm' subgroup = 2.14 
A. Correct behaviour: with adverse_outcome = TRUE the flagged subgroup is z == 1 with IRR > 1 (about exp(0.9) = 2.46).
                  z
count_flagged_harm   0   1
             FALSE   0 390
             TRUE  410   0
                      z
reference_flagged_harm   0   1
                 FALSE 410   0
                 TRUE    0 390
A. forestsearch(count, subgroup_method = 'grf', grf_selection = 'tree') sg.harm = {z <= 0} 
   share z==1 in that sg.harm = 0 | IRR(treat) in it = 0.445  (correct: z == 1, IRR > 1)
B. True IRR(treat): z==0 about 2.23 , z==1 about 4.06 ; observed IRR in z==1 = 4.054 
B. count frontier  : sg.harm.id = NULL (no subgroup) | $candidates present = FALSE 
B. reference       : sg.harm.id = z > 0 & x <= -0.87109222988074286 | $candidates present = TRUE 
B. forestsearch(count, subgroup_method = 'grf') [default frontier + effect] sg.harm = NULL (no subgroup) 
B. Correct behaviour: a harm subgroup is reported (every candidate has IRR well above the 1.25 floor); per the trace the count path reports none.
REPRO_VERDICT: demonstrated
```

### V021

- Finding: A-034 — In a propensity-adjusted GLM run (ps_method != "none" with ps_adjust_method = "iptw" or "dr_gcomp") the search ranks candidates on the adjusted effect, but the 
- Runs: exit 0, 5.3 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)

## ---- Part A (decides the verdict): closure and spec built as forestsearch() builds them ----
## closure: forestsearch_main.R:2704-2712; spec: forestsearch_main.R:3616-3625
set.seed(20261109)
n  <- 2000L
x  <- rbinom(n, 1, 0.5)
ps <- plogis(-1.5 + 3 * x)                                  # strong confounding by x
Treat <- rbinom(n, 1, ps)
Y  <- rbinom(n, 1, plogis(-1.5 - 0.5 * Treat + 2.5 * x))    # treatment protective given x
pt <- mean(Treat)
df <- data.frame(id = seq_len(n), Y = Y, Treat = Treat, q1.1 = 1L,
                 ps_hat = ps,
                 sw = ifelse(Treat == 1, pt / ps, (1 - pt) / (1 - ps)),
                 ips_covar = ifelse(Treat == 1, 1 / ps, 1 / (1 - ps)))
est <- make_effect_estimator(outcome_type = "binary", treat.name = "Treat",
                             outcome.name = "Y", effect_measure = "OR",
                             ps_adjust_method = "iptw")
spec <- list(outcome_type = "binary", effect_measure = "OR", treat.name = "Treat",
             outcome.name = "Y", offset.name = NULL, adjust_covariates = NULL,
             adverse_outcome = TRUE, comparison_threshold = 0)
gate <- .glm_resample_supported(spec)

b_iptw  <- suppressWarnings(est(df))$estimate               # what the search ranks on
b_unadj <- unname(coef(glm(Y ~ Treat, family = binomial(), data = df))[["Treat"]])
## frame shaped like df.fs (forestsearch_main.R:3391-3409): no sw / ps_hat / ips_covar
df_fs_like <- df[, c("Y", "Treat", "id", "q1.1")]
b_nosw  <- suppressWarnings(est(df_fs_like))$estimate

found <- data.frame(HR = b_iptw, n = n, E = sum(Y), grp = 1)
run <- function(data, method, nsplit) {
  r <- suppressWarnings(evaluate_subgroup_consistency(
    m = 1, index.Z = data.frame(q1.1 = 1), names.Z = "q1.1", df = data,
    found.hrs = found, n.splits = nsplit, hr.consistency = 0,
    pconsistency.threshold = 0, pconsistency.digits = 4, maxk = 1,
    confs_labels = "x", estimator_fn = est, consistency_threshold = 0,
    consistency_method = method, glm_resample_spec = spec))
  if (is.null(r)) NA_real_ else as.numeric(r[["Pcons"]])
}
## n.splits = 0: a real split fallback would return NULL (NA here); a number means the closed form ran
p_resample <- run(df, "resample", 0L)
cr <- consistency_resample(df, method = "closed", outcome_type = "binary",
                           effect_measure = "OR", comparison_threshold = 0,
                           treat.name = "Treat", outcome.name = "Y")
p_closed_unadj <- round(cr$rate_closed, 4)
set.seed(1); p_split_iptw    <- run(df, "split", 200L)          # splitting WITH the IPTW estimator (what the docs promise)
set.seed(1); p_split_fsframe <- run(df_fs_like, "split", 200L)  # splitting on a df.fs-like frame (what a fallback would really do)

cat("gate .glm_resample_supported(PS-adjusted OR spec):", gate, "  (documented: FALSE / fall back to split)\n")
cat("IPTW log-OR (search ranking effect):   ", round(b_iptw, 4), "\n")
cat("unadjusted log-OR:                     ", round(b_unadj, 4), "\n")
cat("closed-form beta_hat behind Pcons:     ", round(cr$beta_hat, 4), "\n")
cat("estimator_fn on frame without sw:      ", round(b_nosw, 4), "  (silently unadjusted)\n")
cat("Pcons returned, method = 'resample':   ", p_resample, "\n")
cat("unadjusted closed-form rate:           ", p_closed_unadj, "\n")
cat("Pcons by splitting the IPTW estimator: ", p_split_iptw, "  (rate of the estimator actually ranked)\n")
cat("Pcons by splitting on df.fs-like frame:", p_split_fsframe, "\n")

ok <- isTRUE(gate) &&
  !is.na(p_resample) && isTRUE(abs(p_resample - p_closed_unadj) < 1e-8) &&
  isTRUE(abs(b_iptw - b_unadj) > 0.5) &&
  isTRUE(abs(cr$beta_hat - b_unadj) < 1e-6) &&
  isTRUE(abs(b_nosw - b_unadj) < 1e-6) &&
  !is.na(p_split_iptw) && isTRUE(abs(p_resample - p_split_iptw) > 0.5)

## ---- Part B (supplementary, cannot change the verdict): what forestsearch() hands the engine ----
ns <- asNamespace("forestsearch")
seen <- new.env()
e2e <- tryCatch({
  suppressMessages(trace("subgroup.consistency", where = ns, print = FALSE,
    tracer = bquote({
      assign("cols", names(df), envir = .(seen))
      assign("spec", glm_resample_spec, envir = .(seen))
      assign("method", consistency_method, envir = .(seen))
    })))
  set.seed(42); N <- 300L
  d2 <- data.frame(id = seq_len(N), treat = rbinom(N, 1L, 0.5),
                   age = round(rnorm(N, 55, 12)), biomarker = rnorm(N),
                   biomarker_hi = rbinom(N, 1L, 0.5), sex = rbinom(N, 1L, 0.5))
  in_harm <- d2$age <= 50 & d2$biomarker_hi == 1
  d2$y <- rnorm(N, mean = d2$treat * ifelse(in_harm, 1.5, 0), sd = 1)
  fit <- suppressWarnings(suppressMessages(forestsearch(
    df.analysis = d2, confounders.name = c("age", "biomarker", "biomarker_hi", "sex"),
    outcome.name = "y", event.name = "y", treat.name = "treat", id.name = "id",
    outcome_type = "continuous", effect_measure = "MD", adverse_outcome = TRUE,
    effect.threshold = 0.5, consistency.threshold = 0.0, pconsistency.threshold = 0.5,
    fs.splits = 20, n.min = 30, maxk = 2, use_lasso = FALSE, use_grf = FALSE,
    sg_focus = "eff", is.RCT = FALSE, ps_method = "logistic", ps_adjust_method = "iptw",
    details = FALSE, quiet = TRUE, seedit = 1L, parallel_args = list(plan = "sequential"))))
  if (is.null(seen$cols)) {
    "subgroup.consistency() was not reached in the end-to-end run"
  } else {
    paste0("consistency frame has PS columns: ",
           any(c("sw", "ps_hat", "ips_covar") %in% seen$cols),
           " (adjusted consistency needs TRUE); consistency_method = ", seen$method[1],
           "; gate on forestsearch()'s own spec = ", .glm_resample_supported(seen$spec),
           " (documented fallback needs FALSE); spec fields: ",
           paste(names(seen$spec), collapse = ", "))
  }
}, error = function(e) paste("supplementary end-to-end check did not run:", conditionMessage(e)),
   finally = try(suppressMessages(untrace("subgroup.consistency", where = ns)), silent = TRUE))
cat("END-TO-END (supplementary):", e2e, "\n")

cat("REPRO_VERDICT:", if (ok) "demonstrated" else "not demonstrated", "\n")
```

```text
gate .glm_resample_supported(PS-adjusted OR spec): TRUE   (documented: FALSE / fall back to split)
IPTW log-OR (search ranking effect):    -0.4281 
unadjusted log-OR:                      0.9279 
closed-form beta_hat behind Pcons:      0.9279 
estimator_fn on frame without sw:       0.9279   (silently unadjusted)
Pcons returned, method = 'resample':    1 
unadjusted closed-form rate:            1 
Pcons by splitting the IPTW estimator:  0   (rate of the estimator actually ranked)
Pcons by splitting on df.fs-like frame: 1 
END-TO-END (supplementary): consistency frame has PS columns: FALSE (adjusted consistency needs TRUE); consistency_method = resample; gate on forestsearch()'s own spec = TRUE (documented fallback needs FALSE); spec fields: outcome_type, effect_measure, treat.name, outcome.name, offset.name, adjust_covariates, adverse_outcome, comparison_threshold 
REPRO_VERDICT: demonstrated
```

### V022

- Finding: A-070 — subgroup.search() documents Z as 'Matrix or data frame' and its own example passes a data frame, but Z is never coerced, so a data-frame Z reaches `t(x) %*% x` 
- Runs: exit 0, 2.7 s. Verdict: demonstrated

```r
Sys.setenv(LANGUAGE = "en")
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
old_plan <- future::plan()
set.seed(1); n <- 300
Y  <- rexp(n); E <- rbinom(n, 1, 0.8); Tr <- rbinom(n, 1, 0.5)
Zm <- cbind(z1 = rbinom(n, 1, 0.5), z2 = rbinom(n, 1, 0.5), z3 = rbinom(n, 1, 0.4))
Zd <- as.data.frame(Zm)
run <- function(Z) {
  out <- tryCatch(
    subgroup.search(Y = Y, Event = E, Treat = Tr, Z = Z,
                    n.min = 20, d0.min = 5, d1.min = 5,
                    hr.threshold = 0.5, maxk = 2, parallel_workers = 1),
    error = function(e) e)
  future::plan(future::sequential)
  out
}
res_mat <- run(Zm)
res_df  <- run(Zd)
mat_ok <- !inherits(res_mat, "error")
df_err <- inherits(res_df, "error")
# localise the failure to the cited line (has_positive_variance, subgroup_search.R:771)
hpv_mat <- tryCatch(has_positive_variance(Zm[, 1, drop = FALSE]), error = function(e) e)
hpv_df  <- tryCatch(has_positive_variance(Zd[, 1, drop = FALSE]), error = function(e) e)
cat("class(Zm):", class(Zm)[1], " class(Zd):", class(Zd)[1], "\n")
cat("matrix Z     : ", if (mat_ok) paste0("ran; combinations evaluated = ", res_mat$filter_counts$n_evaluated, ", candidates = ", NROW(res_mat$out.found$hr.subgroups)) else paste("ERROR:", conditionMessage(res_mat)), "\n")
cat("data.frame Z : ", if (df_err) paste("ERROR:", conditionMessage(res_df)) else paste0("ran; combinations evaluated = ", res_df$filter_counts$n_evaluated, ", candidates = ", NROW(res_df$out.found$hr.subgroups)), "\n")
cat("has_positive_variance(matrix column)    : ", if (inherits(hpv_mat, "error")) paste("ERROR:", conditionMessage(hpv_mat)) else hpv_mat, "\n")
cat("has_positive_variance(data.frame column): ", if (inherits(hpv_df, "error")) paste("ERROR:", conditionMessage(hpv_df)) else hpv_df, "\n")
cat("expected (roxygen: 'Matrix or data frame'): the data.frame call runs and returns the same candidates as the matrix call\n")
future::plan(old_plan)
demonstrated <- isTRUE(mat_ok) && isTRUE(df_err) && inherits(hpv_df, "error") && !inherits(hpv_mat, "error")
cat("REPRO_VERDICT:", if (demonstrated) "demonstrated" else "not demonstrated", "\n")
```

```text
Caught simpleError. Canceling all iterations ...
class(Zm): matrix  class(Zd): data.frame 
matrix Z     :  ran; combinations evaluated = 6, candidates = 6 
data.frame Z :  ERROR: requires numeric/complex matrix/vector arguments 
has_positive_variance(matrix column)    :  TRUE 
has_positive_variance(data.frame column):  ERROR: requires numeric/complex matrix/vector arguments 
expected (roxygen: 'Matrix or data frame'): the data.frame call runs and returns the same candidates as the matrix call
REPRO_VERDICT: demonstrated
```

### V023

- Finding: A-074 — Supplying the documented optional hr_1a/hr_0a arguments to sg_tables() (or calling SG_tab_estimates() with SG_flag = "ITT" and a non-NA hr_1a) always stops with
- Runs: exit 0, 2.7 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 200
df <- data.frame(treat = rbinom(n, 1, 0.5), treat.recommend = rbinom(n, 1, 0.6))
df$tte <- rexp(n, 0.05 * exp(-0.5 * df$treat))
df$event <- rbinom(n, 1, 0.7)
df$y <- rbinom(n, 1, 0.4)
run <- function(expr) tryCatch({ suppressWarnings(expr); "OK" }, error = function(e) paste("ERROR:", conditionMessage(e)))
tab <- function(flag, h1 = NA, h0 = NA) SG_tab_estimates(df, SG_flag = flag, outcome.name = "tte", event.name = "event", treat.name = "treat", sg1_name = "Recommend", sg0_name = "Questionable", hr_1a = h1, hr_0a = h0)
h1 <- "0.90 (0.70, 1.20)"; h0 <- "1.50 (1.00, 2.20)"
# (a) exported SG_tab_estimates() directly
a1 <- run(tab("ITT"))                       # no adjusted HR: works
a2 <- run(tab("treat.recommend", h1, h0))   # subgroup rows with adjusted HR: works (the bootstrap's call)
a3 <- run(tab("ITT", h1, h0))               # ITT row with adjusted HR: errors
# (b) sg_tables() survival path with a minimal fs-like object (only the fields Table 1 reads)
fs <- list(df.est = df, sg.harm = "{x<=1}", args_call_all = list(outcome.name = "tte", event.name = "event", treat.name = "treat", est.scale = "hr", maxk = 1))
b1 <- run(sg_tables(fs))                          # documented default: works
b2 <- run(sg_tables(fs, hr_1a = h1, hr_0a = h0))  # documented optional arguments: errors
# (c) GLM table builder: ITT frame (7 columns) stacked on the subgroup frame with adjusted effects (8 columns), as sg_tables() line 1043 does
glm_itt <- SG_tab_estimates_glm(df, SG_flag = "ITT", outcome.name = "y", treat.name = "treat")
glm_sg  <- SG_tab_estimates_glm(df, SG_flag = "treat.recommend", outcome.name = "y", treat.name = "treat", effect_a_1 = "0.01 (0.00, 0.02)", effect_a_0 = "0.03 (0.01, 0.05)")
c1 <- run(as.data.frame(rbind(glm_itt, glm_sg)))
cat("SG_tab_estimates ITT, hr_1a = NA         :", a1, " | expected: OK\n")
cat("SG_tab_estimates subgroup, hr_1a supplied:", a2, " | expected: OK\n")
cat("SG_tab_estimates ITT, hr_1a supplied     :", a3, " | expected (correct implementation): OK\n")
cat("sg_tables(fs)                            :", b1, " | expected: OK\n")
cat("sg_tables(fs, hr_1a=, hr_0a=)            :", b2, " | expected (correct implementation): OK\n")
cat("GLM: ncol ITT =", ncol(glm_itt), ", ncol subgroup =", ncol(glm_sg), "; rbind:", c1, " | expected (correct implementation): equal widths, OK\n")
direct <- identical(a1, "OK") && identical(a2, "OK") && grepl("^ERROR", a3) && grepl("'names' attribute", a3, fixed = TRUE)
viasg <- identical(b1, "OK") && grepl("^ERROR", b2) && grepl("'names' attribute", b2, fixed = TRUE)
cat("REPRO_VERDICT:", if (direct && viasg) "demonstrated" else "not demonstrated", "\n")
```

```text
SG_tab_estimates ITT, hr_1a = NA         : OK  | expected: OK
SG_tab_estimates subgroup, hr_1a supplied: OK  | expected: OK
SG_tab_estimates ITT, hr_1a supplied     : ERROR: 'names' attribute [9] must be the same length as the vector [8]  | expected (correct implementation): OK
sg_tables(fs)                            : OK  | expected: OK
sg_tables(fs, hr_1a=, hr_0a=)            : ERROR: 'names' attribute [9] must be the same length as the vector [8]  | expected (correct implementation): OK
GLM: ncol ITT = 7 , ncol subgroup = 8 ; rbind: ERROR: numbers of columns of arguments do not match  | expected (correct implementation): equal widths, OK
REPRO_VERDICT: demonstrated
```

### V024

- Finding: A-075 — With est.scale = "1/hr" the ITT branch of SG_tab_estimates() never reads est.scale while the subgroup branch reverses the treatment arms, so the sg_tables() est
- Runs: exit 0, 2.8 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 400
df <- data.frame(treat = rbinom(n, 1, 0.5), treat.recommend = rbinom(n, 1, 0.6))
df$tte <- rexp(n, 0.05 * exp(-0.7 * df$treat))
df$event <- rbinom(n, 1, 0.8)
f <- function(flag, sc) suppressWarnings(SG_tab_estimates(df, SG_flag = flag, outcome.name = "tte", event.name = "event", treat.name = "treat", sg1_name = "Recommend", sg0_name = "Questionable", est.scale = sc))
itt_hr <- f("ITT", "hr"); itt_inv <- f("ITT", "1/hr")
sg_hr <- f("treat.recommend", "hr"); sg_inv <- f("treat.recommend", "1/hr")
# the ITT row on the same reversed coding the subgroup rows use under "1/hr"
df2 <- df; df2$treat2 <- 1 - df2$treat
expected <- unname(suppressWarnings(analyze_subgroup(df2, "tte", "event", "treat2", NULL, "ITT", NA, NULL, TRUE, nrow(df2))))
names(expected) <- names(itt_inv)
cols <- c("n1", "m1", "m0", "HR (95% CI)")
cat("Rows as stacked by sg_tables() line 1043 under est.scale = '1/hr':\n")
print(rbind(itt_inv, sg_inv)[, c("Subgroup", cols)])
cat("\nITT row, observed vs correct (same coding as the subgroup rows):\n")
print(rbind(observed = itt_inv[cols], expected = expected[cols]))
# the same table through sg_tables() itself, minimal fs-like object (printed evidence only)
fs <- list(df.est = df, sg.harm = "{x<=1}", args_call_all = list(outcome.name = "tte", event.name = "event", treat.name = "treat", est.scale = "1/hr", maxk = 1))
via <- tryCatch(as.data.frame(suppressWarnings(sg_tables(fs))$tab_estimates[["_data"]])[, c("Subgroup", cols)], error = function(e) paste("could not extract gt data:", conditionMessage(e)))
cat("\nsg_tables(fs)$tab_estimates with fs$args_call_all$est.scale = '1/hr':\n"); print(via)
num <- function(x) as.numeric(sub(" .*$", "", x))
hr_sg_hr <- num(sg_hr[, "HR (95% CI)"]); hr_sg_inv <- num(sg_inv[, "HR (95% CI)"])
cat("\nsubgroup HR under 'hr':", hr_sg_hr, " under '1/hr':", hr_sg_inv, "\n")
cat("ITT HR under 'hr':", num(itt_hr["HR (95% CI)"]), " observed under '1/hr':", num(itt_inv["HR (95% CI)"]), " correct under '1/hr':", num(expected["HR (95% CI)"]), "\n")
itt_ignores <- identical(itt_hr, itt_inv)
sg_reversed <- all(abs(log(hr_sg_inv) + log(hr_sg_hr)) < 0.05) && !identical(sg_hr, sg_inv)
itt_differs <- !identical(unname(itt_inv), unname(expected))
cat("ITT row ignores est.scale:", itt_ignores, "| subgroup rows inverted:", sg_reversed, "| ITT row differs from reversed-coding row:", itt_differs, "\n")
cat("REPRO_VERDICT:", if (itt_ignores && sg_reversed && itt_differs) "demonstrated" else "not demonstrated", "\n")
```

```text
Rows as stacked by sg_tables() line 1043 under est.scale = '1/hr':
        Subgroup       n1            m1     m0     HR (95% CI)        
itt_inv "ITT"          "183 (45.8%)" "41.2" "20.1" "0.49 (0.39, 0.62)"
res_0   "Questionable" "91 (54.5%)"  "23.3" "33.2" "1.64 (1.15, 2.35)"
res_1   "Recommend"    "126 (54.1%)" "16.3" "44"   "2.42 (1.77, 3.30)"

ITT row, observed vs correct (same coding as the subgroup rows):
         n1            m1     m0     HR (95% CI)        
observed "183 (45.8%)" "41.2" "20.1" "0.49 (0.39, 0.62)"
expected "217 (54.2%)" "20.1" "41.2" "2.04 (1.61, 2.58)"

sg_tables(fs)$tab_estimates with fs$args_call_all$est.scale = '1/hr':
      Subgroup          n1   m1   m0       HR (95% CI)
1          ITT 183 (45.8%) 41.2 20.1 0.49 (0.39, 0.62)
2 Questionable  91 (54.5%) 23.3 33.2 1.64 (1.15, 2.35)
3    Recommend 126 (54.1%) 16.3   44 2.42 (1.77, 3.30)

subgroup HR under 'hr': 0.61 0.41  under '1/hr': 1.64 2.42 
ITT HR under 'hr': 0.49  observed under '1/hr': 0.49  correct under '1/hr': 2.04 
ITT row ignores est.scale: TRUE | subgroup rows inverted: TRUE | ITT row differs from reversed-coding row: TRUE 
REPRO_VERDICT: demonstrated
```

### V025

- Finding: A-087 — summary.cox_ahr_cde() applies `$knots` to the splines::ns() basis (an atomic matrix whose knots are an attribute), so the exported summary method stops with a `
- Runs: exit 0, 2.4 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1)
n <- 200
df <- data.frame(
  os_time   = rexp(n, 0.05),
  os_event  = rbinom(n, 1, 0.7),
  treat     = rbinom(n, 1, 0.5),
  biomarker = rnorm(n, 50, 10)
)
df$loghr_po <- -0.05 * (df$biomarker - 50)

# Default configuration of the wrapper (spline_df = 3, hr_threshold = 0.7,
# plot_select = "all", plot_style = "combined", save_plots = FALSE).
# Null device: the wrapper draws base plots; nothing is written.
grDevices::pdf(NULL)
res <- tryCatch(
  suppressWarnings(cox_ahr_cde_analysis(df, verbose = FALSE)),
  error = function(e) e
)
invisible(grDevices::dev.off())

if (inherits(res, "error")) {
  cat("cox_ahr_cde_analysis() itself failed:", conditionMessage(res), "\n")
  cat("REPRO_VERDICT: not demonstrated\n")
} else {
  basis <- res$cox_fit$spline_basis
  cat("class(res)          :", paste(class(res), collapse = ", "), "\n")
  cat("class(spline_basis) :", paste(class(basis), collapse = ", "),
      "| is.atomic:", is.atomic(basis), "\n")

  # the expression used at cox_ahr_cde_wrapper.R:849, in isolation
  dollar_err <- tryCatch({ basis$knots; NA_character_ },
                         error = function(e) conditionMessage(e))

  # print method (shows the object itself is usable)
  print_err <- tryCatch({ invisible(utils::capture.output(print(res))); NA_character_ },
                        error = function(e) conditionMessage(e))

  # summary method: keep what was printed before the failure
  summ_err <- NA_character_
  out <- utils::capture.output(
    summ_err <- tryCatch({ summary(res); NA_character_ },
                         error = function(e) conditionMessage(e))
  )
  printed_primary <- any(grepl("Primary effect", out, fixed = TRUE))
  printed_df      <- any(grepl("Spline df", out, fixed = TRUE))
  printed_diag    <- any(grepl("Diagnostics", out, fixed = TRUE))

  intended <- length(attr(basis, "knots")) + 1L

  cat("basis$knots error   :", if (is.na(dollar_err)) "<none>" else dollar_err, "\n")
  cat("print(res) error    :", if (is.na(print_err)) "<none>" else print_err, "\n")
  cat("summary(res) error  :", if (is.na(summ_err)) "<none>" else summ_err, "\n")
  cat("summary printed 'Primary effect':", printed_primary,
      "| 'Spline df':", printed_df, "| 'Diagnostics':", printed_diag, "\n")
  cat("expected (correct implementation): no error; a 'Spline df:' line showing",
      intended, "(spline_df passed = 3), then the print block and Diagnostics\n")

  demonstrated <- is.atomic(basis) &&
    !is.na(dollar_err) && !is.na(summ_err) && identical(summ_err, dollar_err) &&
    is.na(print_err) &&
    printed_primary && !printed_df && !printed_diag &&
    isTRUE(intended == 3L)
  cat(if (demonstrated) "REPRO_VERDICT: demonstrated\n" else "REPRO_VERDICT: not demonstrated\n")
}
```

```text
class(res)          : cox_ahr_cde, list 
class(spline_basis) : ns, basis, matrix | is.atomic: TRUE 
basis$knots error   : $ operator is invalid for atomic vectors 
print(res) error    : <none> 
summary(res) error  : $ operator is invalid for atomic vectors 
summary printed 'Primary effect': TRUE | 'Spline df': FALSE | 'Diagnostics': FALSE 
expected (correct implementation): no error; a 'Spline df:' line showing 3 (spline_df passed = 3), then the print block and Diagnostics
REPRO_VERDICT: demonstrated
```

### V026

- Finding: A-096 — The within-subgroup refit helpers and print.dina_subgroup_bootstrap() still assume a single-cut subgroup (scalar `if (direction == "left")`, `if (covariate %in%
- Runs: exit 0, 3.0 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(11)
n <- 600
df <- data.frame(w = rbinom(n, 1, 0.5), x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
df$y <- df$w * (df$x1 + df$x2) + rnorm(n, sd = 0.5)
covs <- c("x1", "x2")
fit <- dina(df, outcome = "y", treatment = "w", covariates = covs, family = "gaussian", seed = 1L)
# max_depth is left at its default (2L)
sg <- dina_subgroup(fit, df, covariates = covs, m_diff = 0, n_min = 60L, sg_focus = "hr")
cat("R version:", as.character(getRversion()), "\n")
cat("found:", sg$found, " depth:", sg$depth, " label:", sg$label, "\n")
cat("covariate:", sg$covariate, " direction:", sg$direction, " threshold:", format(sg$threshold), "\n")
is_depth2 <- isTRUE(sg$found) && identical(as.integer(sg$depth), 2L)

# Reference values a depth-aware refit would report
expected_effect <- NA_real_
expected_sig <- NA_character_
if (isTRUE(sg$found)) {
  expected_effect <- unname(stats::coef(stats::glm(y ~ w, data = df[sg$mask, , drop = FALSE], family = "gaussian"))["w"])
  expected_sig <- paste(paste(sg$covariate, ifelse(sg$direction == "left", "<=", ">="), format(sg$threshold)), collapse = " & ")
}

# (1) exported dina_subgroup_refit() on the default-depth result
r1 <- tryCatch(dina_subgroup_refit(sg, df, treatment = "w", outcome = "y", covariates = covs),
               error = function(e) e)
refit_errored <- inherits(r1, "error")
refit_sig_wrong <- FALSE
if (refit_errored) {
  cat("(1) dina_subgroup_refit(): ERROR:", conditionMessage(r1), "\n")
  cat("    expected: effect =", format(expected_effect), " signature =", expected_sig, "\n")
} else {
  refit_sig_wrong <- length(r1$signature) != 1L ||
    !(grepl("x1", r1$signature[1L], fixed = TRUE) && grepl("x2", r1$signature[1L], fixed = TRUE))
  cat("(1) dina_subgroup_refit(): effect =", format(r1$effect), " expected =", format(expected_effect),
      "\n    signature =", paste(r1$signature, collapse = " | "), " expected =", expected_sig, "\n")
}

# (2) bootstrap with its defaults (refit = TRUE, refit_confounders = "none")
warns <- character(0)
bs <- tryCatch(
  withCallingHandlers(
    dina_subgroup_bootstrap(df = df, outcome = "y", treatment = "w", covariates = covs,
                            family = "gaussian", m_diff = 0, n_min = 60L, sg_focus = "hr",
                            n_boot = 10L, seed = 1L, sg = sg),
    warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") }),
  error = function(e) e)
boot_errored <- inherits(bs, "error")
refit_dropped <- FALSE
print_errored <- FALSE
if (boot_errored) {
  cat("(2) dina_subgroup_bootstrap(): ERROR:", conditionMessage(bs), "\n")
} else {
  cat("(2) bootstrap warnings:\n"); for (w in warns) cat("    -", w, "\n")
  refit_dropped <- is.null(bs$refit) && all(is.na(bs$refit_effect_ci))
  cat("    bs$refit is NULL:", is.null(bs$refit), " refit_effect_ci:", format(bs$refit_effect_ci),
      " (expected: a refit object with effect", format(expected_effect), "and a finite CI, or an explicit 'depth-2 not supported' message)\n")
  cat("    DINA effect_ci (unaffected):", format(bs$effect_ci), "\n")
  p <- tryCatch({ utils::capture.output(print(bs)); "ok" }, error = function(e) paste("ERROR:", conditionMessage(e)))
  print_errored <- !identical(p, "ok")
  cat("    print(bs):", p, " (expected: ok)\n")
}

# (3) bootstrap with a character refit_confounders: the scalar %in% test sits outside the tryCatch
bs3 <- tryCatch(
  suppressWarnings(dina_subgroup_bootstrap(df = df, outcome = "y", treatment = "w", covariates = covs,
                            family = "gaussian", m_diff = 0, n_min = 60L, sg_focus = "hr",
                            n_boot = 4L, seed = 1L, sg = sg, refit_confounders = "x2")),
  error = function(e) e)
boot3_errored <- inherits(bs3, "error")
cat("(3) bootstrap with refit_confounders = \"x2\":",
    if (boot3_errored) paste("ERROR:", conditionMessage(bs3)) else "completed",
    " (expected: completed)\n")

demonstrated <- is_depth2 && (refit_errored || refit_sig_wrong) &&
  (boot_errored || refit_dropped || print_errored || boot3_errored)
if (!is_depth2) cat("NOTE: the search did not select a depth-2 subgroup on this data, so the defect was not exercised.\n")
cat("REPRO_VERDICT:", if (demonstrated) "demonstrated" else "not demonstrated", "\n")
```

```text
R version: 4.6.2 
found: TRUE  depth: 2  label: {x1 >= 0.232128357328474 & x2 >= 0.462484237272292} 
covariate: x1 x2  direction: right right  threshold: 0.2321284 0.4624842 
(1) dina_subgroup_refit(): ERROR: the condition has length > 1 
    expected: effect = 1.117932  signature = x1 >= 0.2321284 & x2 >= 0.4624842 
(2) bootstrap warnings:
    - `n_boot` < 100; bootstrap percentile CI will be unstable. Consider increasing `n_boot`. 
    - Within-subgroup standard-model point estimate failed (the condition has length > 1); skipping the within-subgroup refit. 
    bs$refit is NULL: TRUE  refit_effect_ci: NA NA  (expected: a refit object with effect 1.117932 and a finite CI, or an explicit 'depth-2 not supported' message)
    DINA effect_ci (unaffected): 1.142458 1.280723 
    print(bs): ERROR: the condition has length > 1  (expected: ok)
(3) bootstrap with refit_confounders = "x2": ERROR: the condition has length > 1  (expected: completed)
REPRO_VERDICT: demonstrated
```

### V027

- Finding: A-101 — When the full-data analysis has a subgroup but no held-out subject of a CV run is assigned to a CV subgroup, the one-row agreement table is indexed as if its fi
- Runs: exit 0, 2.4 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)

# Shape that forestsearch_Kfold() (R/forestsearch_cross_validation.R:517, 596-607) and
# forestsearch_tenfold() (:1082, 1149-1168) hand to forestsearch_KfoldOut() when the
# full-data analysis found a subgroup but NO fold found one: every held-out row gets
# treat.recommend = 1 (not in H); coding is 0 = in H (get_dfpred).
cvidx <- rep(1:4, each = 25L)
orig  <- c(rep(0, 35), rep(1, 65))   # full-data analysis: 35 of 100 subjects in H
base  <- list(
  cv_args = list(confounders.name = c("age", "sex"), outcome.name = "time",
                 event.name = "event", treat.name = "treat"),
  sg0.name = "Not recommend", sg1.name = "Recommend", Kfolds = 4L,
  sg_analysis = c("{age <= 50}", "{sex >= 1}"))

res_zero <- c(base, list(resCV = data.frame(
  cvindex = cvidx, sg1 = NA_character_, sg2 = NA_character_,
  treat.recommend = 1, treat.recommend.original = orig,
  stringsAsFactors = FALSE)))
out_zero <- forestsearch_KfoldOut(res_zero, outall = FALSE)

rec <- res_zero$resCV$treat.recommend
cat("Agreement table (rows = CV flag, cols = full-data flag; 0 = in H):\n")
print(table(treat.recommend = rec, treat.recommend.original = orig))

# Definitions (those of the 2x2 branch, :1411-1414): sens_H = P(CV-H | orig-H),
# ppv_H = P(orig-H | CV-H)
sens_H_def <- sum(rec == 0 & orig == 0) / sum(orig == 0)
ppv_H_def  <- if (sum(rec == 0) > 0) sum(rec == 0 & orig == 0) / sum(rec == 0) else NA_real_

obs_any  <- out_zero$find_metrics[["Any"]]
obs_sens <- out_zero$sens_metrics_original[["sens_H"]]
obs_ppv  <- out_zero$sens_metrics_original[["ppv_H"]]
cat("find_metrics['Any'] (share of folds identifying):", obs_any, "\n")
cat("sens_H observed:", obs_sens, " | by definition:", sens_H_def, "\n")
cat("ppv_H  observed:", obs_ppv,  " | by definition:", ppv_H_def, "\n")

# What a user sees for a single K-fold run of this shape (forestsearch_Kfold()
# returns exactly these two vectors as sens_summary / find_summary, :616-617)
cat("cv_summary_text():",
    cv_summary_text(list(find_summary = out_zero$find_metrics,
                         sens_summary = out_zero$sens_metrics_original)), "\n")

# Effect on the forestsearch_tenfold() summary: same rbind + apply(median, na.rm = TRUE)
# as :1201 and :1229. Two zero-identification simulations plus one simulation with
# real agreement 30/35.
rec2  <- c(rep(0, 40), rep(1, 60))
orig2 <- c(rep(0, 30), rep(1, 10), rep(0, 5), rep(1, 55))
fold_sg1 <- c("{age <= 50}", "{age <= 60}", NA, "{sex >= 1}")
fold_sg2 <- c("{sex >= 1}", NA, NA, NA)
res_ident <- c(base, list(resCV = data.frame(
  cvindex = cvidx, sg1 = fold_sg1[cvidx], sg2 = fold_sg2[cvidx],
  treat.recommend = rec2, treat.recommend.original = orig2,
  stringsAsFactors = FALSE)))
out_ident <- forestsearch_KfoldOut(res_ident, outall = FALSE)
sens_out <- rbind(out_zero$sens_metrics_original, out_zero$sens_metrics_original,
                  out_ident$sens_metrics_original)
med_obs <- apply(sens_out, 2, median, na.rm = TRUE)[["sens_H"]]
med_def0  <- median(c(sens_H_def, sens_H_def, 30 / 35))      # zero-identification sims counted as 0
med_defNA <- median(c(NA, NA, 30 / 35), na.rm = TRUE)        # ... or skipped as NA
cat("tenfold-style median sens_H observed:", med_obs,
    " | zero-id sims as 0:", med_def0, " | zero-id sims as NA:", med_defNA, "\n")

demonstrated <- isTRUE(obs_any == 0) &&
  isTRUE(all.equal(obs_sens, 1)) &&
  isTRUE(all.equal(sens_H_def, 0)) &&
  !is.na(obs_ppv) && is.na(ppv_H_def) &&
  isTRUE(all.equal(med_obs, 1)) && !isTRUE(all.equal(med_obs, med_def0)) &&
  !isTRUE(all.equal(med_obs, med_defNA))
cat("REPRO_VERDICT:", if (demonstrated) "demonstrated" else "not demonstrated", "\n")
```

```text
Agreement table (rows = CV flag, cols = full-data flag; 0 = in H):
               treat.recommend.original
treat.recommend  0  1
              1 35 65
find_metrics['Any'] (share of folds identifying): 0 
sens_H observed: 1  | by definition: 0 
ppv_H  observed: 0.35  | by definition: NA 
cv_summary_text(): CV found = 0%, Agree(+,-) = NA%, 100% 
tenfold-style median sens_H observed: 1  | zero-id sims as 0: 0  | zero-id sims as NA: 0.8571429 
REPRO_VERDICT: demonstrated
```

### V028

- Finding: A-107 — With outcome_type_approx left at its NULL default and no outcome_type in fs_params (the documented minimum for a survival analysis), the P1 approximation silent
- Runs: exit 0, 3.1 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
future::plan(future::sequential)
set.seed(1)
n <- 200
df <- data.frame(id = seq_len(n), treat = rbinom(n, 1, 0.5), tte = rexp(n),
                 event = rbinom(n, 1, 0.7), x1 = rnorm(n))
# Survival analysis. 'zz_no_such_arg' makes each forestsearch() replicate fail at argument
# matching (forestsearch has no '...'), so the script is fast; d_eff and P1 are computed
# before the replicates and do not depend on them.
cat("forestsearch() has '...':", "..." %in% names(formals(forestsearch)), "\n")
base <- list(outcome.name = "tte", event.name = "event", treat.name = "treat",
             id.name = "id", confounders.name = "x1", zz_no_such_arg = 1)
run <- function(fs_params, prop_cens = 0.3) tryCatch(
  suppressWarnings(fpr_calibration(c1 = 1.25, c2 = 1.0, df.analysis = df,
                                   fs_params = fs_params, nsims = 1L, n_min = 60L,
                                   prop_cens = prop_cens, n_workers = 1L, quiet = TRUE)),
  error = function(e) e)
r_omit  <- run(base)                                        # outcome_type left at forestsearch() default
r_omit2 <- run(base, prop_cens = 0.6)                       # same, different censoring proportion
r_surv  <- run(c(base, list(outcome_type = "survival")))   # same analysis, default spelled out
fs_default <- eval(formals(forestsearch)$outcome_type)[1]
cat("forestsearch() default outcome_type:", fs_default, "\n")
errs <- Filter(function(r) inherits(r, "error"), list(r_omit, r_omit2, r_surv))
if (length(errs) > 0) {
  for (e in errs) cat("call error:", conditionMessage(e), "\n")
  demonstrated <- FALSE
} else {
  expected_d <- d_eff_survival(60, 0.3)                     # 42 expected events
  cat("otype used : omitted =", r_omit$call_args$otype, "| explicit survival =", r_surv$call_args$otype, "| correct = survival\n")
  cat("d_eff      : omitted =", r_omit$d_eff, "| explicit survival =", r_surv$d_eff, "| correct =", expected_d, "\n")
  cat("d_eff with prop_cens 0.3 vs 0.6, type omitted:", r_omit$d_eff, "vs", r_omit2$d_eff,
      "| correct =", d_eff_survival(60, 0.3), "vs", d_eff_survival(60, 0.6), "\n")
  cat("P1         : omitted =", r_omit$P1, "| explicit survival =", r_surv$P1, "| correct =", r_surv$P1, "\n")
  cat("L_eff at a hypothetical FPR of 0.10: omitted =", log(0.9) / log(1 - r_omit$P1),
      "| correct =", log(0.9) / log(1 - r_surv$P1), "\n")
  demonstrated <- identical(fs_default, "survival") &&
    identical(r_omit$call_args$otype, "binary") &&
    isTRUE(all.equal(r_surv$d_eff, expected_d)) &&
    !isTRUE(all.equal(r_omit$d_eff, expected_d)) &&
    isTRUE(all.equal(r_omit$d_eff, r_omit2$d_eff)) &&
    is.finite(r_omit$P1) && is.finite(r_surv$P1) &&
    !isTRUE(all.equal(r_omit$P1, r_surv$P1))
}
cat("REPRO_VERDICT:", if (isTRUE(demonstrated)) "demonstrated" else "not demonstrated", "\n")
```

```text
forestsearch() has '...': FALSE 
forestsearch() default outcome_type: survival 
otype used : omitted = binary | explicit survival = survival | correct = survival
d_eff      : omitted = 15 | explicit survival = 42 | correct = 42 
d_eff with prop_cens 0.3 vs 0.6, type omitted: 15 vs 15 | correct = 42 vs 24 
P1         : omitted = 0.22222 | explicit survival = 0.180038 | correct = 0.180038 
L_eff at a hypothetical FPR of 0.10: omitted = 0.4192426 | correct = 0.5307907 
REPRO_VERDICT: demonstrated
```

### V029

- Finding: A-111 — The replicate harness counts a failure only on a thrown error, a NULL return or a 'subgroup.consistency\|subscript' warning, so a forestsearch() call that return
- Runs: exit 0, 3.7 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
future::plan(future::sequential)
set.seed(1)
n <- 200
d <- data.frame(id = seq_len(n), y_sim = rexp(n), event_sim = rbinom(n, 1, 0.7),
                treat_sim = rbinom(n, 1, 0.5), flag_harm = 0L, er = rnorm(n))
dn <- list(outcome = "y_sim", event = "event_sim", treat = "treat_sim",
           id = "id", harm = "flag_harm")
fsp <- list(mr_inference = TRUE,
            mr_inference_args = list(ci_method = "ij", draws = 200L),
            quiet = TRUE, details = FALSE, plot.sg = FALSE)

# The same forestsearch() call the harness builds, made directly, with warnings collected.
direct <- function(confs, extra = list()) {
  w <- character(0)
  args <- c(list(df.analysis = d, outcome.name = "y_sim", event.name = "event_sim",
                 treat.name = "treat_sim", id.name = "id", flag_harm.name = "flag_harm",
                 confounders.name = confs, hr.threshold = 1.25, hr.consistency = 1.0,
                 seedit = 1L, parallel_args = list(plan = "sequential")), fsp, extra)
  fs <- withCallingHandlers(
    tryCatch(do.call(forestsearch, args), error = function(e) e),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  threw <- inherits(fs, "error")
  list(threw = threw,
       failed_shape = !threw && is.null(fs$sg.harm) && !is.null(fs$error_log),
       stage = if (!threw && !is.null(fs$error_log)) fs$error_log$stage else NA_character_,
       warnings = w)
}

# Trigger (a): analyst confounders absent from the simulated frame -> empty intersection
out_a <- .fs_fdr_one_rep(d, 1L, fsp, c1 = 1.25, c2 = 1.0, dgm_names = dn,
                         conf_use = c("bmi", "ecog"))
dir_a <- direct(character(0))
cat("(a) direct call: threw =", dir_a$threw, " failure-shaped return =", dir_a$failed_shape,
    " error_log stage =", dir_a$stage, "\n")
cat("(a) direct call warnings:\n"); print(dir_a$warnings)
cat("(a) harness flag observed:", out_a$flag, " | correct for a failed fit: -1\n")
cat("(a) harness msg observed :", out_a$msg, " | correct: the failure text\n")

# Trigger (b): partial overlap, inherited conf.cont_jcuts names a covariate the intersection drops
jc <- list(conf.cont_jcuts = list(bmi = 4L))
out_b <- .fs_fdr_one_rep(d, 1L, c(fsp, jc), c1 = 1.25, c2 = 1.0, dgm_names = dn,
                         conf_use = c("er", "bmi"))
dir_b <- direct("er", jc)
cat("(b) direct call: threw =", dir_b$threw, " failure-shaped return =", dir_b$failed_shape,
    " error_log stage =", dir_b$stage, "\n")
cat("(b) direct call warnings:\n"); print(dir_b$warnings)
cat("(b) harness flag observed:", out_b$flag, " | correct for a failed fit: -1\n")
demo_b <- isTRUE(dir_b$failed_shape) && identical(out_b$flag, 0L) && is.na(out_b$msg)
cat("(b) second trigger demonstrated:", demo_b, "\n")

# Supplementary: aggregate level (does not enter the verdict)
agg <- tryCatch(
  .fs_fdr_run_null(gen_df = function(seed) d, label = "all-fail", fs_params = fsp,
                   c1 = 1.25, c2 = 1.0, dgm_names = dn, conf_use = c("bmi", "ecog"),
                   nsims = 3L, n_workers = 1L, seed_base = 100L, c_confirm = 1.0,
                   n = n, event_min = 5L, quiet = TRUE),
  error = function(e) { cat("aggregate step not run:", conditionMessage(e), "\n"); NULL })
future::plan(future::sequential)
if (!is.null(agg))
  cat("aggregate observed: n_valid =", agg$n_valid, " fs_errors =", agg$fs_errors,
      " unhandled =", agg$unhandled, " A =", agg$rows$A[1],
      " | correct: n_valid = 0, fs_errors = 3, A = NA\n")

demonstrated <- isTRUE(dir_a$failed_shape) && identical(out_a$flag, 0L) && is.na(out_a$msg)
cat("REPRO_VERDICT:", if (demonstrated) "demonstrated" else "not demonstrated", "\n")
```

```text
(a) direct call: threw = FALSE  failure-shaped return = TRUE  error_log stage = get_FSdata 
(a) direct call warnings:
[1] "Error in get_FSdata: 'confounders.name' must be a non-empty character vector."
[2] "FSdata failure - returning NULL result"                                       
(a) harness flag observed: 0  | correct for a failed fit: -1
(a) harness msg observed : NA  | correct: the failure text
(b) direct call: threw = FALSE  failure-shaped return = TRUE  error_log stage = get_FSdata 
(b) direct call warnings:
[1] "Error in get_FSdata: conf.cont_jcuts names not in confounders.name: 'bmi'."
[2] "FSdata failure - returning NULL result"                                    
(b) harness flag observed: 0  | correct for a failed fit: -1
(b) second trigger demonstrated: TRUE 
aggregate observed: n_valid = 3  fs_errors = 0  unhandled = 0  A = 0  | correct: n_valid = 0, fs_errors = 3, A = NA
REPRO_VERDICT: demonstrated
```

### V030

- Finding: A-114 — The exported generate_aft_dgm_flex() has the formal default treatment_var = NULL and its roxygen says NULL means treatment is simulated 50/50, but no code imple
- Runs: exit 0, 2.6 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
df <- survival::gbsg
build <- function(...) {
  tryCatch(
    suppressMessages(suppressWarnings(generate_aft_dgm_flex(
      data            = df,
      continuous_vars = c("age", "size", "nodes", "pgr", "er"),
      factor_vars     = "meno",
      outcome_var     = "rfstime",
      event_var       = "status",
      model           = "null",
      n_super         = 500,
      verbose         = FALSE,
      ...
    ))),
    error = function(e) e
  )
}
# Control: same call with the treatment column named (the roxygen example; expected to work)
ctrl <- build(treatment_var = "hormon")
# Documented default: treatment_var omitted (= NULL)
res  <- build()
cat("formal default of treatment_var:", deparse(formals(generate_aft_dgm_flex)$treatment_var), "\n")
cat("control (treatment_var = 'hormon') class:", class(ctrl)[1], "\n")
if (inherits(ctrl, "error")) cat("  control error:", conditionMessage(ctrl), "\n")
cat("default (treatment_var = NULL) class:", class(res)[1], "\n")
if (inherits(res, "error")) {
  cat("  observed: error ->", conditionMessage(res), "\n")
} else {
  cat("  observed: DGM returned; mean(df_super$treat) =", mean(res$df_super$treat), "\n")
}
cat("  documented/correct: an aft_dgm_flex object with treatment simulated 50/50 (mean treat near 0.5)\n")
demonstrated <- inherits(ctrl, "aft_dgm_flex") && inherits(res, "error")
cat("REPRO_VERDICT:", if (demonstrated) "demonstrated" else "not demonstrated", "\n")
```

```text
formal default of treatment_var: NULL 
control (treatment_var = 'hormon') class: aft_dgm_flex 
default (treatment_var = NULL) class: simpleError 
  observed: error -> attempt to select less than one element in get1index 
  documented/correct: an aft_dgm_flex object with treatment simulated 50/50 (mean treat near 0.5)
REPRO_VERDICT: demonstrated
```

### V031

- Finding: A-115 — The label panel decides bold/plain by comparing each row colour with cat_colours["ITT"], so any non-NULL cat_colours without an element named "ITT" gives geom_t
- Runs: exit 0, 3.2 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
grDevices::pdf(NULL)   # null device: nothing is written

mk <- function(cc, cv = c("x", "y")) {
  gg_forest(subgroups = c("A", "B"), est = c(1.0, 1.2), lo = c(0.8, 0.9),
            hi = c(1.3, 1.6), cat_vec = cv, cat_colours = cc,
            xlim = c(0.5, 2), ticks_at = c(0.5, 1, 2))
}
build <- function(cc, cv = c("x", "y")) tryCatch(mk(cc, cv), error = function(e) e)
draw <- function(p) {
  if (inherits(p, "error")) return(paste("ERROR at construction:", conditionMessage(p)))
  tryCatch({
    suppressWarnings(print(p))
    "drawn without error"
  }, error = function(e) {
    msg <- conditionMessage(e)
    if (!is.null(e$parent)) msg <- paste(msg, "|", conditionMessage(e$parent))
    paste("ERROR at draw:", gsub("[[:space:]]+", " ", msg))
  })
}
# fontface constant stored on the label panel (first patch of the patchwork)
ff <- function(p) tryCatch(
  paste(format(p[[1]]$layers[[1]]$aes_params$fontface), collapse = ", "),
  error = function(e) "<could not extract>")

p_no  <- build(c(x = "blue", y = "red"))                 # documented form, no ITT entry
p_itt <- build(c(x = "blue", y = "red", ITT = "black"))  # control: same call plus an ITT entry
p_nul <- build(NULL, cv = NULL)                          # default: no categories at all

r_no  <- draw(p_no)
r_itt <- draw(p_itt)
r_nul <- draw(p_nul)

cat("cat_colours without ITT : fontface =", ff(p_no),  "| ", r_no,  "\n")
cat("   correct implementation: fontface = plain, plain | drawn without error\n")
cat("cat_colours with ITT    : fontface =", ff(p_itt), "| ", r_itt, "\n")
cat("   correct implementation: fontface = plain, plain | drawn without error\n")
cat("cat_colours = NULL      : fontface =", ff(p_nul), "| ", r_nul, "(side observation: every label bold)\n")
cat("error message mentions fontface:", grepl("fontface", r_no, fixed = TRUE), "\n")
invisible(grDevices::dev.off())

# demonstrated = the no-ITT map fails while the otherwise identical control draws
demonstrated <- grepl("^ERROR", r_no) && identical(r_itt, "drawn without error")
cat("REPRO_VERDICT:", if (demonstrated) "demonstrated" else "not demonstrated", "\n")
```

```text
cat_colours without ITT : fontface = NA, NA |  ERROR at draw: Problem while converting geom to grob. ℹ Error occurred in the 1st layer. Caused by error in `FUN()`: ! invalid fontface NA | invalid fontface NA 
   correct implementation: fontface = plain, plain | drawn without error
cat_colours with ITT    : fontface = plain, plain |  drawn without error 
   correct implementation: fontface = plain, plain | drawn without error
cat_colours = NULL      : fontface = bold, bold |  drawn without error (side observation: every label bold)
error message mentions fontface: TRUE 
REPRO_VERDICT: demonstrated
```

### V032

- Finding: A-117 — With family = "gaussian", glm_cs_fit() compares the raw deviance difference (a residual-sum-of-squares difference in squared outcome units) to a chi-square with
- Runs: exit 0, 2.4 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1)
n <- 400
df <- data.frame(treat = rbinom(n, 1, 0.5), bm = rnorm(n))
df$y   <- 1 + 0.5 * df$treat + rnorm(n)   # no treat x bm interaction, residual SD 1
df$y10 <- 10 * df$y                        # the same data in other units

# glm_cs_fit() calls par() on entry (line 154), which opens a graphics device;
# a null pdf device is opened first so that no file is written.
grDevices::pdf(NULL)
fits <- tryCatch({
  list(
    f1 = glm_cs_fit(df, outcome_name = "y",   treat_name = "treat", z_name = "bm",
                    family = "gaussian", verbose = FALSE),
    f2 = glm_cs_fit(df, outcome_name = "y10", treat_name = "treat", z_name = "bm",
                    family = "gaussian", verbose = FALSE)
  )
}, finally = grDevices::dev.off())
f1 <- fits$f1
f2 <- fits$f2

# Reference: F test for the same nested gaussian models, from the stored fits
# (model_main = main effects only, model_fit = with treat x spline interaction)
ref_p <- function(f) {
  rss0 <- f$model_main$deviance
  rss1 <- f$model_fit$deviance
  q    <- f$model_main$df.residual - f$model_fit$df.residual
  Fst  <- ((rss0 - rss1) / q) / (rss1 / f$model_fit$df.residual)
  stats::pf(Fst, q, f$model_fit$df.residual, lower.tail = FALSE)
}
stat <- function(f) f$model_main$deviance - f$model_fit$deviance

out <- data.frame(
  outcome       = c("y", "10*y"),
  pkg_statistic = c(stat(f1), stat(f2)),
  pkg_lrt_p     = c(f1$lrt_pvalue, f2$lrt_pvalue),
  correct_F_p   = c(ref_p(f1), ref_p(f2))
)
print(out, digits = 6)
cat("ratio of package statistics (10*y over y):", stat(f2) / stat(f1),
    " (a valid test gives 1)\n")

ref_invariant <- isTRUE(all.equal(ref_p(f1), ref_p(f2), tolerance = 1e-8))
pkg_differs   <- abs(f1$lrt_pvalue - f2$lrt_pvalue) > 1e-6
ratio_100     <- isTRUE(all.equal(stat(f2) / stat(f1), 100, tolerance = 1e-6))
cat("reference p invariant to units:", ref_invariant,
    "| package p changes with units:", pkg_differs,
    "| statistic scales by 100:", ratio_100, "\n")
cat(if (ref_invariant && pkg_differs && ratio_100) "REPRO_VERDICT: demonstrated" else "REPRO_VERDICT: not demonstrated", "\n")
```

```text
  outcome pkg_statistic   pkg_lrt_p correct_F_p
1       y       0.67908 8.78113e-01     0.90163
2    10*y      67.90799 1.19709e-14     0.90163
ratio of package statistics (10*y over y): 100  (a valid test gives 1)
reference p invariant to units: TRUE | package p changes with units: TRUE | statistic scales by 100: TRUE 
REPRO_VERDICT: demonstrated 
Warning message:
In par(oldpar) : calling par(new=TRUE) with no plot
```

### V033

- Finding: A-119 — For a binary outcome the closed-form scorer treats a candidate with an empty 2x2 cell as non-estimable, but the glm scorer used when fast = FALSE or adjust_cova
- Runs: exit 0, 2.4 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
# Binary outcome; candidate g1 has NO events among its 40 controls (one empty 2x2 cell),
# candidate g2 has all four cells populated.
n <- 200
treat <- rep(c(0L, 1L), each = 100)
g1 <- rep(rep(c(1L, 0L), times = c(40, 60)), 2)
g2 <- as.integer(seq_len(n) %% 2 == 0)
y <- integer(n)
y[101:112] <- 1L   # g1, treated: 12 events of 40
y[41:55]   <- 1L   # not g1, control: 15 events of 60
y[141:158] <- 1L   # not g1, treated: 18 events of 60
df <- data.frame(treat = treat, y = y, g1 = g1, g2 = g2)
cells <- function(s) table(treat = df$treat[s], y = factor(df$y[s], levels = 0:1))
cat("2x2 cells, g1:\n"); print(cells(df$g1 == 1))
cat("2x2 cells, g2:\n"); print(cells(df$g2 == 1))

# the model-path scorer on the g1 subgroup (the closed form marks it NA: a cell is 0)
slow_coef <- .g3_coef(df[df$g1 == 1, ], "binary", "treat", NULL, NULL, "y", 5L)
cat("model-path coefficient for g1 (.g3_coef):", format(slow_coef),
    "  | closed form: NA (cell count 0)\n")

run <- function(fast) {
  w <- character(0)
  fit <- withCallingHandlers(
    guohe_algorithm3(df, outcome = "binary", treatment = "treat",
                     candidates = c("g1", "g2"), y = "y", orient = +1,
                     B = 40, seed = 2, fast = fast),
    warning = function(cond) {
      w <<- c(w, conditionMessage(cond))
      invokeRestart("muffleWarning")
    })
  list(fit = fit, warnings = w)
}
a <- run(TRUE)    # closed form (what the default fast = NULL resolves to here)
b <- run(FALSE)   # glm path
row_of <- function(z) c(selected = z$fit$selected,
                        n_estimable = z$fit$n_estimable,
                        naive_OR = format(z$fit$naive, digits = 6),
                        debiased_OR = format(z$fit$debiased, digits = 6),
                        bound = format(z$fit$bound_one_sided, digits = 6))
cat("\nobserved (documented as the SAME estimator; a correct implementation gives identical rows):\n")
print(rbind(fast_TRUE = row_of(a), fast_FALSE = row_of(b)))
cat("warnings fast=TRUE :", if (length(a$warnings)) unique(a$warnings) else "<none>", "\n")
cat("warnings fast=FALSE:", if (length(b$warnings)) unique(b$warnings) else "<none>", "\n")

differs <- !identical(a$fit$selected, b$fit$selected) ||
  a$fit$n_estimable != b$fit$n_estimable ||
  !isTRUE(all.equal(a$fit$naive, b$fit$naive))
demonstrated <- isTRUE(is.finite(slow_coef)) && isTRUE(differs)
cat("REPRO_VERDICT:", if (demonstrated) "demonstrated" else "not demonstrated", "\n")
```

```text
2x2 cells, g1:
     y
treat  0  1
    0 40  0
    1 28 12
2x2 cells, g2:
     y
treat  0  1
    0 43  7
    1 35 15
model-path coefficient for g1 (.g3_coef): 18.71877   | closed form: NA (cell count 0)

observed (documented as the SAME estimator; a correct implementation gives identical rows):
           selected n_estimable naive_OR    debiased_OR bound     
fast_TRUE  "g2"     "1"         "2.63265"   "2.28537"   "0.799029"
fast_FALSE "g1"     "2"         "134728292" "57380465"  "27833183"
warnings fast=TRUE : 1 of 2 candidates not estimable on the full sample; dropped. 1 of 40 resamples had no estimable candidate; dropped. 
warnings fast=FALSE: <none> 
REPRO_VERDICT: demonstrated
```

### V034

- Finding: A-120 — The Cox scorers of the Guo-He comparator accept any subgroup with at least min_events total events, so a bootstrap resample in which one arm of a candidate has 
- Runs: exit 0, 5.2 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
# Survival outcome, guohe_algorithm3() defaults (lean Cox path, min_events = 5, B = 2000,
# r = 0.03) with orient = +1 as fs_to_guohe() passes. On the FULL sample every candidate
# has >= 5 events in EACH arm; only bootstrap resamples can leave an arm of g1 without events.
set.seed(11)
n <- 200
treat <- rep(c(0L, 1L), times = n / 2)
g1 <- as.integer(seq_len(n) <= 60)     # 30 controls, 30 treated
g2 <- as.integer(seq_len(n) > 40)      # 160 subjects, ordinary
t0 <- rexp(n, 0.05)
time <- pmin(t0, 30)
event <- as.integer(t0 <= 30)
c1 <- which(g1 == 1 & treat == 0)      # the 30 controls of g1
time[c1] <- 30; event[c1] <- 0L
time[c1[1:5]] <- c(4, 9, 14, 19, 24); event[c1[1:5]] <- 1L   # exactly 5 control events
df <- data.frame(time = time, event = event, treat = treat, g1 = g1, g2 = g2)
cat("full-sample events by arm (rows: control, treated):\n")
print(sapply(c("g1", "g2"), function(k) {
  s <- df[[k]] == 1
  tapply(df$event[s], df$treat[s], sum)
}))

B <- 2000L; SEED <- 3L
w <- character(0)
fit <- withCallingHandlers(
  guohe_algorithm3(df, outcome = "survival", treatment = "treat",
                   candidates = c("g1", "g2"), time = "time", event = "event",
                   orient = +1, B = B, seed = SEED),
  warning = function(cond) {
    w <<- c(w, conditionMessage(cond))
    invokeRestart("muffleWarning")
  })

# Independent recomputation of the same resamples with the package's own scorer, once as
# the package does it and once treating a resampled candidate with a zero-event arm as
# non-estimable (the condition forestsearch's own Cox estimator applies).
ctl <- survival::coxph.control()
M <- cbind(g1 = df$g1, g2 = df$g2)
arm_events <- function(idx, k) {
  s <- M[idx, k] == 1L
  tb <- df$treat[idx][s]; eb <- df$event[idx][s]
  c(control = sum(eb[tb == 0L]), treated = sum(eb[tb == 1L]))
}
score <- function(idx, per_arm) vapply(1:2, function(k) {
  s <- M[idx, k] == 1L
  if (per_arm && any(arm_events(idx, k) == 0)) return(NA_real_)
  .g3_cox_lean(df$time[idx][s], df$event[idx][s], df$treat[idx][s], 5L, ctl)
}, numeric(1))
obs <- score(seq_len(n), FALSE)
gamma_max <- max(obs)
offset <- (1 - n^(0.03 - 0.5)) * (gamma_max - obs)
set.seed(SEED)
idx_list <- lapply(seq_len(B), function(b) sample.int(n, n, replace = TRUE))
u_of <- function(idx, per_arm) {
  v <- score(idx, per_arm) + offset
  f <- is.finite(v)
  if (any(f)) max(v[f]) - gamma_max else NA_real_
}
U_obs <- vapply(idx_list, u_of, numeric(1), per_arm = FALSE)
U_ref <- vapply(idx_list, u_of, numeric(1), per_arm = TRUE)
U_obs <- U_obs[is.finite(U_obs)]; U_ref <- U_ref[is.finite(U_ref)]
n_mono <- sum(vapply(idx_list, function(idx) {
  e <- arm_events(idx, 1L)
  e[["control"]] == 0 && sum(e) >= 5
}, logical(1)))
z <- stats::qnorm(0.975)
summ <- function(U) {
  d <- gamma_max - mean(U)
  c(debiased_HR = exp(d),
    bound_one_sided = exp(gamma_max - as.numeric(stats::quantile(U, 0.95, names = FALSE))),
    wald_lo = exp(d - z * stats::sd(U)), wald_hi = exp(d + z * stats::sd(U)),
    sd_U = stats::sd(U), max_U = max(U))
}
man <- summ(U_obs); ref <- summ(U_ref)
replicated <- isTRUE(all.equal(fit$debiased, man[["debiased_HR"]])) &&
  isTRUE(all.equal(fit$ci_wald, c(man[["wald_lo"]], man[["wald_hi"]])))

cat("\nselected:", fit$selected, " naive HR:", format(fit$naive, digits = 5),
    " n_estimable:", fit$n_estimable, "\n")
cat("resamples (of", B, ") in which g1 has NO control events yet is scored:", n_mono, "\n")
cat("drop_rate reported for g1:", format(fit$drop_rate[["g1"]]),
    " | warnings:", if (length(w)) paste(unique(w), collapse = " | ") else "<none>", "\n")
cat("manual recomputation reproduces the package output:", replicated, "\n\n")
cat("observed (package):   debiased HR =", format(fit$debiased, digits = 5),
    " bound =", format(fit$bound_one_sided, digits = 5),
    " Wald CI = (", format(fit$ci_wald[1], digits = 5), ",", format(fit$ci_wald[2], digits = 5), ")",
    " sd(U) =", format(man[["sd_U"]], digits = 4), " max(U) =", format(man[["max_U"]], digits = 4), "\n")
cat("reference (zero-event arm = non-estimable, same resamples): debiased HR =",
    format(ref[["debiased_HR"]], digits = 5),
    " bound =", format(ref[["bound_one_sided"]], digits = 5),
    " Wald CI = (", format(ref[["wald_lo"]], digits = 5), ",", format(ref[["wald_hi"]], digits = 5), ")",
    " sd(U) =", format(ref[["sd_U"]], digits = 4), " max(U) =", format(ref[["max_U"]], digits = 4), "\n")

demonstrated <- isTRUE(n_mono > 0) && isTRUE(replicated) &&
  isTRUE(man[["max_U"]] > 8) &&
  isTRUE(abs(log(fit$debiased / ref[["debiased_HR"]])) > 0.02)
cat("REPRO_VERDICT:", if (demonstrated) "demonstrated" else "not demonstrated", "\n")
```

```text
full-sample events by arm (rows: control, treated):
  g1 g2
0  5 51
1 22 61

selected: g1  naive HR: 6.74  n_estimable: 2 
resamples (of 2000 ) in which g1 has NO control events yet is scored: 12 
drop_rate reported for g1: 0  | warnings: <none> 
manual recomputation reproduces the package output: TRUE 

observed (package):   debiased HR = 4.7892  bound = 2.0407  Wald CI = ( 0.24396 , 94.019 )  sd(U) = 1.519  max(U) = 19.88 
reference (zero-event arm = non-estimable, same resamples): debiased HR = 5.3701  bound = 2.1535  Wald CI = ( 2.1884 , 13.177 )  sd(U) = 0.458  max(U) = 2.181 
REPRO_VERDICT: demonstrated
```

### V035

- Finding: A-125 — With confounders.name left at its formal default NULL, mrct_region_sims() hands forestsearch() every z_ column of dgm$df_super, which for any spline DGM (always
- Runs: exit 0, 34.5 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
future::plan(future::sequential)
t0 <- proc.time()[["elapsed"]]
set.seed(11); n <- 600
d <- data.frame(tte = rexp(n, 0.05), event = rbinom(n, 1, 0.75), treat = rep(0:1, n / 2),
                age = round(rnorm(n, 62, 8)), bm = round(runif(n, 0, 12), 1),
                male = rbinom(n, 1, 0.6), histology = rbinom(n, 1, 0.4),
                prior_treat = rbinom(n, 1, 0.3), regA = rbinom(n, 1, 0.3))
dgms <- list(
  alt  = create_dgm_for_mrct(d, model_type = "alt", log_hrs = log(c(3, 1.25, 0.5)), verbose = FALSE),
  null = create_dgm_for_mrct(d, model_type = "null", verbose = FALSE))
dgm <- dgms$alt

# (1) the default set exactly as built at mrct_simulation.R:314-323
v <- names(dgm$df_super)
default_set <- setdiff(v[grepl("^z_|^ecog|^strat", v)],
  c("id", "treat", "tte", "event", "y_sim", "treat_sim", "event_sim", "t_true", "c_time",
    "lin_pred_0", "lin_pred_1", "loghr_po", "entrytime", "flag_harm"))
helpers  <- paste0("z_bm", c("_treat", "_k", "_k_treat"))
intended <- setdiff(default_set, c("z_regA", helpers))
cat("default set  (observed):", default_set, "\n")
cat("intended set (correct) :", intended, "(baseline covariates only)\n")

# (2) one trial drawn as mrct_region_sims() draws it: z_bm_treat is z_bm x SOURCE arm
s <- simulate_from_dgm(dgm, n = 500, rand_ratio = 1, draw_treatment = TRUE,
                       analysis_time = 60, cens_adjust = 0, seed = 101)
cat("z_bm_treat == z_bm * source arm :", all(s$z_bm_treat == s$z_bm * s$treat), "\n")
cat("z_bm_treat == z_bm * trial arm  :", all(s$z_bm_treat == s$z_bm * s$treat_sim), "\n")

# (3) candidate pool forestsearch() builds from the default set on the training (non-region) patients
cf <- c("z_age <= 65", "z_bm <= 0", "z_bm <= 1", "z_bm <= 2", "z_bm <= 5")
fsd <- tryCatch(get_FSdata(df.analysis = as.data.frame(s[s$z_regA == 0, ]),
                           confounders.name = default_set, conf_force = cf,
                           outcome.name = "y_sim", event.name = "event_sim", details = FALSE),
                error = function(e) { cat("get_FSdata error:", conditionMessage(e), "\n"); NULL })
helper_cuts <- if (is.null(fsd)) character(0) else grep("z_bm_treat|z_bm_k", fsd$confs, value = TRUE)
cat("candidate cuts on helper columns (correct: none):\n"); print(helper_cuts)

# (4) same seeds, default set vs intended set (trial data identical; only the candidate set differs).
#     Chunks of 4 replicates, alternating the alt and the null DGM, until about 50 s have been used.
neq <- function(x, y) !((is.na(x) & is.na(y)) | (!is.na(x) & !is.na(y) & abs(x - y) < 1e-8))
run <- function(g, cn, seed) tryCatch(
  mrct_region_sims(g, n_sims = 4, n_sample = 500, region_var = "z_regA",
                   confounders.name = cn, fs_args = list(fs.splits = 100),
                   parallel_args = list(plan = "sequential", show_message = FALSE), seed = seed),
  error = function(e) { cat("mrct_region_sims error:", conditionMessage(e), "\n"); NULL })
tab <- NULL
for (chunk in 1:6) {
  lab <- if (chunk %% 2 == 1) "alt" else "null"
  a <- run(dgms[[lab]], NULL, 100 * chunk)
  b <- run(dgms[[lab]], intended, 100 * chunk)
  if (!is.null(a) && !is.null(b) && nrow(a) == nrow(b)) {
    cat("chunk", chunk, lab, ": same simulated trials in both runs:", isTRUE(all.equal(a$hr_itt, b$hr_itt)), "\n")
    tab <- rbind(tab, data.frame(model = lab, seed = 100 * chunk + a$sim,
      sg_default = a$sg_found, sg_intended = b$sg_found,
      n_sg_default = a$n_sg, n_sg_intended = b$n_sg,
      hr_sg_default = round(a$hr_sg, 3), hr_sg_intended = round(b$hr_sg, 3),
      differs = (a$any_found != b$any_found) | neq(a$n_sg, b$n_sg) | neq(a$hr_sg, b$hr_sg)))
  }
  if (proc.time()[["elapsed"]] - t0 > 50) break
}
if (!is.null(tab)) {
  print(tab)
  cat("S1_EVIDENCE: declared subgroup uses a source-arm column (z_bm_treat / z_bm_k_treat) in",
      sum(grepl("z_bm_treat|z_bm_k_treat", tab$sg_default)), "of", nrow(tab),
      "replicates; uses z_bm_k in", sum(grepl("z_bm_k(?!_treat)", tab$sg_default, perl = TRUE)),
      "; replicates where any_found / n_sg / hr_sg differ between the default and the intended set:",
      sum(tab$differs), "(correct: 0)\n")
}

demonstrated <- all(c("z_regA", helpers) %in% default_set) && length(helper_cuts) > 0
cat("REPRO_VERDICT:", if (demonstrated) "demonstrated" else "not demonstrated", "\n")
```

```text
2          1.081          1.081   FALSE
3          1.515          1.515   FALSE
4          0.623          0.623   FALSE
5          0.753          0.753   FALSE
6          0.659          0.659   FALSE
7          0.524          0.524   FALSE
8          0.794          0.794   FALSE
9          1.310          1.310   FALSE
10         1.007          1.007   FALSE
11         1.197          1.197   FALSE
12         0.820          0.820   FALSE
13         0.713          0.713   FALSE
14         0.792          0.792   FALSE
15         0.737          0.737   FALSE
16         0.807          0.807   FALSE
17         0.821          0.821   FALSE
18         0.969          0.969   FALSE
19         1.203          1.203   FALSE
20         1.265          1.265   FALSE
21         0.679          0.679   FALSE
22         0.755          0.755   FALSE
23         0.588          0.588   FALSE
24         0.613          0.613   FALSE
S1_EVIDENCE: declared subgroup uses a source-arm column (z_bm_treat / z_bm_k_treat) in 0 of 24 replicates; uses z_bm_k in 0 ; replicates where any_found / n_sg / hr_sg differ between the default and the intended set: 0 (correct: 0)
REPRO_VERDICT: demonstrated
```

### V036

- Finding: A-133 — compute_sg_summary() extracts the event totals with vapply(..., integer(1), "events"), so plot_sg_results() and plot() on a forestsearch fit stop with a vapply 
- Runs: exit 0, 2.6 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 200L
d_int <- data.frame(Y = rexp(n), Event = rbinom(n, 1L, 0.7), Treat = rbinom(n, 1L, 0.5),
                    treat.recommend = rbinom(n, 1L, 0.7))
d_dbl <- d_int; d_dbl$Event <- as.numeric(d_dbl$Event)   # same 0/1 values, storage type double
stopifnot(is.integer(d_int$Event), is.double(d_dbl$Event), all(d_int$Event == d_dbl$Event))

# 1. exported function, list input
run <- function(d) {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  tryCatch(plot_sg_results(list(df.est = d, sg.harm = "{x <= 1}"), plot_type = "summary"),
           error = function(e) e)
}
r_int <- run(d_int)
r_dbl <- run(d_dbl)

# 2. S3 entry point plot(<forestsearch>), every type, double event
fs_like <- structure(list(df.est = d_dbl, sg.harm = "{x <= 1}"), class = c("forestsearch", "list"))
run_s3 <- function(type) {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  tryCatch(plot(fs_like, type = type, outcome.name = "Y", event.name = "Event", treat.name = "Treat"),
           error = function(e) e)
}
types <- c("combined", "km", "forest", "summary")
r_s3 <- lapply(types, run_s3)

# 3. localise to the internal helper
s_dbl <- tryCatch(
  compute_sg_summary(df = d_dbl,
                     df_H  = d_dbl[d_dbl$treat.recommend == 0, ],
                     df_Hc = d_dbl[d_dbl$treat.recommend == 1, ],
                     outcome.name = "Y", event.name = "Event", treat.name = "Treat",
                     sg0_name = "H", sg1_name = "Hc"),
  error = function(e) e)

expected_events <- c(sum(d_dbl$Event), sum(d_dbl$Event[d_dbl$treat.recommend == 0]),
                     sum(d_dbl$Event[d_dbl$treat.recommend == 1]))
describe <- function(r) if (inherits(r, "error")) paste("ERROR:", gsub("\\s+", " ", conditionMessage(r))) else paste("OK, class", class(r)[1])
is_vapply_type_error <- function(r) inherits(r, "error") &&
  grepl("integer", conditionMessage(r)) && grepl("double", conditionMessage(r))

cat("integer Event -> plot_sg_results():", describe(r_int), "\n")
if (inherits(r_int, "fs_sg_plot")) cat("  observed events (ITT, H, Hc), integer run:", r_int$summary$events, "\n")
cat("double  Event -> plot_sg_results():", describe(r_dbl), "\n")
for (i in seq_along(types)) cat("double  Event -> plot(<forestsearch>, type = '", types[i], "'): ", describe(r_s3[[i]]), "\n", sep = "")
cat("double  Event -> compute_sg_summary():", describe(s_dbl), "\n")
cat("  correct result for both storage types, events (ITT, H, Hc):", expected_events, "\n")

ok_int   <- inherits(r_int, "fs_sg_plot") && isTRUE(all(r_int$summary$events == expected_events))
fail_dbl <- is_vapply_type_error(r_dbl) && is_vapply_type_error(s_dbl) &&
  all(vapply(r_s3, is_vapply_type_error, logical(1)))
cat("REPRO_VERDICT:", if (ok_int && fail_dbl) "demonstrated" else "not demonstrated", "\n")
```

```text
integer Event -> plot_sg_results(): OK, class fs_sg_plot 
  observed events (ITT, H, Hc), integer run: 138 48 90 
double  Event -> plot_sg_results(): ERROR: values must be type 'integer', but FUN(X[[1]]) result is type 'double' 
double  Event -> plot(<forestsearch>, type = 'combined'): ERROR: values must be type 'integer', but FUN(X[[1]]) result is type 'double'
double  Event -> plot(<forestsearch>, type = 'km'): ERROR: values must be type 'integer', but FUN(X[[1]]) result is type 'double'
double  Event -> plot(<forestsearch>, type = 'forest'): ERROR: values must be type 'integer', but FUN(X[[1]]) result is type 'double'
double  Event -> plot(<forestsearch>, type = 'summary'): ERROR: values must be type 'integer', but FUN(X[[1]]) result is type 'double'
double  Event -> compute_sg_summary(): ERROR: values must be type 'integer', but FUN(X[[1]]) result is type 'double' 
  correct result for both storage types, events (ITT, H, Hc): 138 48 90 
REPRO_VERDICT: demonstrated
```

### V037

- Finding: A-138 — In GLM simulation mode the closure that fills hr.itt / hr.H.hat / hr.Hc.hat / hr.H.true / hr.Hc.true returns a log odds ratio for effect_measure = "RD", an odds
- Runs: exit 0, 2.5 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(2)
n  <- 400
tr <- rep(0:1, n / 2)
y  <- rbinom(n, 1, ifelse(tr == 1, 0.5, 0.3))
df <- data.frame(id = seq_len(n), y_sim = y, event_sim = y,
                 treat_sim = tr, flag_harm = 0L)
extract_fn <- get(".extract_fs_estimates_gen", envir = asNamespace("forestsearch"))
run1 <- function(em, data = df, fs_full = NULL, ot = "binary") {
  extract_fn(
    df = data, fs_res = NULL, fs_full = fs_full, dgm = NULL,
    cox_formula = NULL, cox_formula_adj = NULL,
    outcome_name = "y_sim", event_name = "event_sim", treat_name = "treat_sim",
    harm_col = "flag_harm", analysis = "FS", verbose = FALSE,
    outcome_type = ot, effect_measure = em, offset_name = NULL
  )
}
itt <- function(em) as.numeric(run1(em)$hr.itt)

# Part 1: binary outcome, ITT effect per measure
p1 <- mean(y[tr == 1]); p0 <- mean(y[tr == 0])
or <- (p1 / (1 - p1)) / (p0 / (1 - p0))
rep_OR   <- itt("OR")
rep_RD   <- itt("RD")
rep_RR   <- itt("RR")
rep_NULL <- itt(NULL)
cat(sprintf("OR   : reported %.6f | correct (sample OR) %.6f\n", rep_OR, or))
cat(sprintf("RD   : reported %.6f | correct (p1 - p0) %.6f | log(OR) %.6f\n", rep_RD, p1 - p0, log(or)))
cat(sprintf("RR   : reported %.6f | correct (p1 / p0) %.6f | OR %.6f\n", rep_RR, p1 / p0, or))
cat(sprintf("NULL : reported %s | correct (NULL resolves to OR in forestsearch()) %.6f\n", format(rep_NULL), or))
tol <- 1e-4
control_ok <- isTRUE(abs(rep_OR - or) < tol)
rd_wrong   <- isTRUE(abs(rep_RD - log(or)) < tol) && isTRUE(abs(rep_RD - (p1 - p0)) > 1e-2)
rr_wrong   <- isTRUE(abs(rep_RR - or) < tol) && isTRUE(abs(rep_RR - p1 / p0) > 1e-2)
null_na    <- is.na(rep_NULL)
cat("control_ok =", control_ok, " rd_is_logOR =", rd_wrong,
    " rr_is_OR =", rr_wrong, " null_is_NA =", null_na, "\n")

# Part 2 (extra, informative): continuous outcome, effect_measure omitted, subgroup found
set.seed(3)
g   <- rep(0:1, each = n / 2)
mu0 <- rep(10, n)
mu1 <- ifelse(g == 1, 13, 11)
yc  <- ifelse(tr == 1, mu1, mu0) + rnorm(n)
dfc <- data.frame(id = seq_len(n), y_sim = yc, event_sim = 1L, treat_sim = tr,
                  flag_harm = g, mu0 = mu0, mu1 = mu1)
fsf <- list(df.est = data.frame(treat.recommend = 1L - g))
rc  <- run1(NULL, data = dfc, fs_full = fsf, ot = "continuous")
odds_on_means <- (13 / (1 - 13)) / (10 / (1 - 10))
cat(sprintf("continuous, measure NULL: hr.H.hat reported %s | correct (MD in subgroup) about 3\n",
            format(as.numeric(rc$hr.H.hat))))
cat(sprintf("continuous, measure NULL: ahr.H.hat reported %.6f | correct (mean(mu1 - mu0) in subgroup) %.6f | odds formula on means %.6f\n",
            as.numeric(rc$ahr.H.hat), 3, odds_on_means))
extra <- is.na(as.numeric(rc$hr.H.hat)) &&
  isTRUE(abs(as.numeric(rc$ahr.H.hat) - odds_on_means) < 1e-8)
cat("extra_continuous_null_ahr_is_odds_formula =", extra, "\n")

cat("REPRO_VERDICT:", if (control_ok && rd_wrong && rr_wrong && null_na) "demonstrated" else "not demonstrated", "\n")
```

```text
OR   : reported 1.595238 | correct (sample OR) 1.595238
RD   : reported 0.467023 | correct (p1 - p0) 0.110000 | log(OR) 0.467023
RR   : reported 1.595238 | correct (p1 / p0) 1.333333 | OR 1.595238
NULL : reported NA | correct (NULL resolves to OR in forestsearch()) 1.595238
control_ok = TRUE  rd_is_logOR = TRUE  rr_is_OR = TRUE  null_is_NA = TRUE 
continuous, measure NULL: hr.H.hat reported NA | correct (MD in subgroup) about 3
continuous, measure NULL: ahr.H.hat reported 0.975000 | correct (mean(mu1 - mu0) in subgroup) 3.000000 | odds formula on means 0.975000
extra_continuous_null_ahr_is_odds_formula = TRUE 
REPRO_VERDICT: demonstrated
```

### V038

- Finding: A-140 — setup_gbsg_dgm() does not put a lin_pred_obs column on df_super, so a direct call simulate_from_dgm(<GBSG dgm>) with n left at its default NULL stops with an op
- Runs: exit 0, 2.4 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
n_super <- 1000L   # same size as the fixture in tests/testthat/test-knoise-population.R
dgm <- suppressWarnings(setup_gbsg_dgm(model = "alt", n_super = n_super, verbose = FALSE))
stopifnot(inherits(dgm, "aft_dgm_flex"), nrow(dgm$df_super) == n_super)
has_obs  <- "lin_pred_obs" %in% names(dgm$df_super)
has_true <- "lin.conf.true" %in% names(dgm$df_super)
cat("lin_pred_obs in df_super :", has_obs, " (the n = NULL path reads this column)\n")
cat("lin.conf.true in df_super:", has_true, "\n")

# 1. default n = NULL on the object as built
r_null <- tryCatch(simulate_from_dgm(dgm, seed = 1), error = function(e) e)
null_failed <- inherits(r_null, "error")
msg <- if (null_failed) conditionMessage(r_null) else ""
opaque <- null_failed && grepl("replacement has 0 row", msg, fixed = TRUE)
cat("observed, n = NULL (default):", if (null_failed) paste("ERROR:", msg) else paste("rows =", nrow(r_null)), "\n")
cat("correct                     : a data frame with", n_super, "rows (roxygen: 'uses the entire super population')\n")
cat("error is the unguarded zero-length assignment:", opaque, "\n")

# 2. controls: explicit n on the same object works
r_n <- tryCatch(simulate_from_dgm(dgm, n = 100, seed = 1), error = function(e) e)
n_ok <- is.data.frame(r_n) && nrow(r_n) == 100
cat("n = 100                     :", if (n_ok) paste("rows =", nrow(r_n)) else paste("ERROR:", conditionMessage(r_n)), "\n")
r_all <- tryCatch(simulate_from_dgm(dgm, n = n_super, replace = FALSE, seed = 1), error = function(e) e)
all_ok <- is.data.frame(r_all) && nrow(r_all) == n_super
cat("n = n_super, replace = FALSE:", if (all_ok) paste("rows =", nrow(r_all)) else paste("ERROR:", conditionMessage(r_all)), "\n")

# 3. same object with the missing column supplied in memory (what a fix would give)
fix_ok <- FALSE
if (has_true) {
  dgm2 <- dgm
  dgm2$df_super$lin_pred_obs <- dgm2$df_super$lin.conf.true
  r_fix <- tryCatch(simulate_from_dgm(dgm2, seed = 1), error = function(e) e)
  fix_ok <- is.data.frame(r_fix) && nrow(r_fix) == n_super &&
    all(is.finite(r_fix$y_sim)) && all(r_fix$event_sim %in% c(0L, 1L))
  cat("n = NULL with lin_pred_obs added:", if (is.data.frame(r_fix)) paste("rows =", nrow(r_fix), " events =", sum(r_fix$event_sim)) else paste("ERROR:", conditionMessage(r_fix)), "\n")
}

demonstrated <- (!has_obs) && null_failed && n_ok && all_ok && fix_ok
cat("REPRO_VERDICT:", if (demonstrated) "demonstrated" else "not demonstrated", "\n")
```

```text
lin_pred_obs in df_super : FALSE  (the n = NULL path reads this column)
lin.conf.true in df_super: TRUE 
observed, n = NULL (default): ERROR: replacement has 0 rows, data has 1000 
correct                     : a data frame with 1000 rows (roxygen: 'uses the entire super population')
error is the unguarded zero-length assignment: TRUE 
n = 100                     : rows = 100 
n = n_super, replace = FALSE: rows = 1000 
n = NULL with lin_pred_obs added: rows = 1000  events = 447 
REPRO_VERDICT: demonstrated
```

### V039

- Finding: A-143 — The benefit-notation inversion list omits cde.H.true and cde.Hc.true, so format_oc_results(subgroup_notation = "benefit") shows the true-subgroup CDE means on t
- Runs: exit 0, 2.5 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
res <- data.frame(
  sim = 1:4, analysis = "FS", any.H = 1L,
  sens = 0.8, spec = 0.9, ppv = 0.8, npv = 0.9, size.H = 50, hr.itt = 1,
  hr.H.hat = 2, hr.Hc.hat = 0.5, hr.H.true = 2, hr.Hc.true = 0.5,
  ahr.H.hat = 2, ahr.Hc.hat = 0.5, ahr.H.true = 2, ahr.Hc.true = 0.5,
  cde.H.true = 2, cde.Hc.true = 0.5, cde.H.hat = 2, cde.Hc.hat = 0.5,
  stringsAsFactors = FALSE
)
h <- as.data.frame(format_oc_results(res, subgroup_notation = "harm", use_gt = FALSE))
b <- as.data.frame(format_oc_results(res, subgroup_notation = "benefit", use_gt = FALSE))
cols <- c("HR_H_true", "HR_H_hat", "AHR_H_true", "AHR_H_hat",
          "CDE_H_true", "CDE_H_hat", "CDE_Hc_true", "CDE_Hc_hat")
cat("harm notation (input scale):\n"); print(h[, cols])
cat("benefit notation (observed):\n"); print(b[, cols])
expected <- 1 / h[, cols]
cat("benefit notation (correct = 1/harm for every effect column):\n"); print(expected)
same <- function(x, y) isTRUE(all.equal(as.numeric(x), as.numeric(y)))
hr_inverted   <- same(b$HR_H_true, 1 / h$HR_H_true) && same(b$HR_H_hat, 1 / h$HR_H_hat)
hat_inverted  <- same(b$CDE_H_hat, 1 / h$CDE_H_hat) && same(b$CDE_Hc_hat, 1 / h$CDE_Hc_hat)
true_inverted <- same(b$CDE_H_true, 1 / h$CDE_H_true) && same(b$CDE_Hc_true, 1 / h$CDE_Hc_true)
cat(sprintf("CDE_H_true observed = %g, correct = %g; CDE_H_hat observed = %g, correct = %g\n",
            b$CDE_H_true, 1 / h$CDE_H_true, b$CDE_H_hat, 1 / h$CDE_H_hat))
cat(sprintf("CDE_Hc_true observed = %g, correct = %g; CDE_Hc_hat observed = %g, correct = %g\n",
            b$CDE_Hc_true, 1 / h$CDE_Hc_true, b$CDE_Hc_hat, 1 / h$CDE_Hc_hat))
cat("REPRO_VERDICT:", if (hr_inverted && hat_inverted && !true_inverted) "demonstrated" else "not demonstrated", "\n")
```

```text
harm notation (input scale):
  HR_H_true HR_H_hat AHR_H_true AHR_H_hat CDE_H_true CDE_H_hat CDE_Hc_true
1         2        2          2         2          2         2         0.5
  CDE_Hc_hat
1        0.5
benefit notation (observed):
  HR_H_true HR_H_hat AHR_H_true AHR_H_hat CDE_H_true CDE_H_hat CDE_Hc_true
1       0.5      0.5        0.5       0.5          2       0.5         0.5
  CDE_Hc_hat
1          2
benefit notation (correct = 1/harm for every effect column):
  HR_H_true HR_H_hat AHR_H_true AHR_H_hat CDE_H_true CDE_H_hat CDE_Hc_true
1       0.5      0.5        0.5       0.5        0.5       0.5           2
  CDE_Hc_hat
1          2
CDE_H_true observed = 2, correct = 0.5; CDE_H_hat observed = 0.5, correct = 0.5
CDE_Hc_true observed = 0.5, correct = 2; CDE_Hc_hat observed = 2, correct = 2
REPRO_VERDICT: demonstrated
```

### V040

- Finding: A-147 — Both cut() calls use R's default right-closed intervals while the bin labels describe left-closed bins, so a value exactly on a break is counted one row below i
- Runs: exit 0, 2.3 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)

# Pcons as the pipeline produces it: a rate rounded to 2 digits, passed through
# character and back (subgroup_consistency_helpers.R:1560-1572,
# subgroup_consistency_main.R:1103-1111)
rates <- c(0.9512, 0.9549, 0.9451, 0.96, 0.90, 1.00)
pcons <- as.numeric(as.character(round(rates, 2)))   # 0.95 0.95 0.95 0.96 0.90 1.00
nsg   <- as.numeric(c(50, 100, 150, 200, 300, 99))   # N is numeric in the bootstrap results
cat("Pcons values:", format(pcons, nsmall = 2), "\n")
cat("Pcons values identical to the break 0.95:", sum(pcons == 0.95), "\n")

res <- data.frame(
  boot_id   = 1:6,
  any_found = 1L,
  Pcons     = pcons,
  hr_sg     = 2,
  N_sg      = nsg,
  K_sg      = 1,
  M.1       = "{x <= 1}",
  M.2       = NA_character_,
  stringsAsFactors = FALSE
)

out <- suppressWarnings(
  summarize_bootstrap_subgroups(results = res, nb_boots = 6, maxk = 2)
)

cd     <- as.data.frame(out$consistency_dist)
sd_tab <- as.data.frame(out$size_dist)
print(cd)
print(sd_tab)

get_count <- function(tab, label) {
  v <- tab$Count[tab[[1]] == label]
  if (length(v) == 0L) NA_integer_ else v
}

obs_p_top <- get_count(cd, ">=0.95");     exp_p_top <- sum(pcons >= 0.95)
obs_p_mid <- get_count(cd, "0.85-0.95");  exp_p_mid <- sum(pcons >= 0.85 & pcons < 0.95)
obs_n_top <- get_count(sd_tab, ">=300");  exp_n_top <- sum(nsg >= 300)
obs_n_low <- get_count(sd_tab, "<50");    exp_n_low <- sum(nsg < 50)
obs_n_50  <- get_count(sd_tab, "50-99");  exp_n_50  <- sum(nsg >= 50 & nsg <= 99)

cat(sprintf("Pcons '>=0.95'   : observed %s, label implies %d\n", obs_p_top, exp_p_top))
cat(sprintf("Pcons '0.85-0.95': observed %s, label implies %d\n", obs_p_mid, exp_p_mid))
cat(sprintf("N_sg  '>=300'    : observed %s, label implies %d\n", obs_n_top, exp_n_top))
cat(sprintf("N_sg  '<50'      : observed %s, label implies %d\n", obs_n_low, exp_n_low))
cat(sprintf("N_sg  '50-99'    : observed %s, label implies %d (the two observed are N = 100 and N = 99; N = 50 is under '<50')\n",
            obs_n_50, exp_n_50))

# What the CV sibling's convention (right = FALSE) gives for the same Pcons values
print(table(cut(pcons, breaks = c(0, 0.5, 0.7, 0.8, 0.85, 0.95, 1.0),
                labels = c("<0.5", "0.5-0.7", "0.7-0.8", "0.8-0.85", "0.85-0.95", ">=0.95"),
                include.lowest = TRUE, right = FALSE)))

mismatch <- !isTRUE(obs_p_top == exp_p_top) || !isTRUE(obs_p_mid == exp_p_mid) ||
            !isTRUE(obs_n_top == exp_n_top) || !isTRUE(obs_n_low == exp_n_low)
cat("REPRO_VERDICT:", if (mismatch) "demonstrated" else "not demonstrated", "\n")
```

```text
Pcons values: 0.95 0.95 0.95 0.96 0.90 1.00 
Pcons values identical to the break 0.95: 3 
  Consistency Range Count  Percent
1              <0.5     0  0.00000
2           0.5-0.7     0  0.00000
3           0.7-0.8     0  0.00000
4          0.8-0.85     0  0.00000
5         0.85-0.95     4 66.66667
6            >=0.95     2 33.33333
  Size Range Count  Percent
1        <50     1 16.66667
2      50-99     2 33.33333
3    100-149     1 16.66667
4    150-199     1 16.66667
5    200-299     1 16.66667
6      >=300     0  0.00000
Pcons '>=0.95'   : observed 2, label implies 5
Pcons '0.85-0.95': observed 4, label implies 1
N_sg  '>=300'    : observed 0, label implies 1
N_sg  '<50'      : observed 1, label implies 0
N_sg  '50-99'    : observed 2, label implies 2 (the two observed are N = 100 and N = 99; N = 50 is under '<50')

     <0.5   0.5-0.7   0.7-0.8  0.8-0.85 0.85-0.95    >=0.95 
        0         0         0         0         1         5 
REPRO_VERDICT: demonstrated
```

### V041

- Finding: A-151 — test_hte_from_forestsearch() uses fit$confounders.evaluated as column names of the user's data, but a default forestsearch() fit stores cut expressions such as 
- Runs: exit 0, 4.5 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(11)
n <- 300
df <- data.frame(id = seq_len(n), treat = rbinom(n, 1, 0.5),
                 age = round(rnorm(n, 55, 12)),
                 bmi = round(rnorm(n, 27, 4), 1),
                 bm  = rbinom(n, 1, 0.5))
harm <- df$age <= 50 & df$bm == 1
tt <- rexp(n, 0.05 * exp(log(2.5) * df$treat * harm))
df$time <- pmin(tt, 60)
df$event <- as.integer(tt <= 60)

## Default pre-processing (use_grf = FALSE, use_lasso = FALSE, cut_type = "default",
## subgroup_method = "consistency"); only search-size settings are reduced for speed.
run_fs <- function(confs) {
  tryCatch(suppressWarnings(suppressMessages(forestsearch(
    df.analysis = df, outcome.name = "time", event.name = "event",
    treat.name = "treat", id.name = "id", confounders.name = confs,
    sg_focus = "hr", use_grf = FALSE, use_lasso = FALSE, fs.splits = 50,
    n.min = 30, d0.min = 5, d1.min = 5, hr.threshold = 1.1,
    hr.consistency = 1.0, pconsistency.threshold = 0.5, maxk = 2,
    details = FALSE, quiet = TRUE, plot.sg = FALSE,
    parallel_args = list(plan = "sequential", workers = 1L,
                         show_message = FALSE)))),
    error = function(e) { cat("forestsearch error:", conditionMessage(e), "\n"); NULL })
}
run_hte <- function(fit) {
  w <- character(0)
  res <- withCallingHandlers(
    test_hte_from_forestsearch(fit, df, outcome.name = "time", treat.name = "treat"),
    warning = function(cond) { w <<- c(w, conditionMessage(cond)); invokeRestart("muffleWarning") })
  list(res = res, warnings = w)
}

## Case 1: two continuous + one binary covariate
confs1 <- c("age", "bmi", "bm")
fit1 <- run_fs(confs1)
ev1 <- fit1$confounders.evaluated
cat("Case 1 class of fit:", paste(class(fit1), collapse = "/"), "\n")
cat("Case 1 confounders.evaluated:\n"); print(ev1)
cat("Case 1 confounders.candidate:\n"); print(fit1$confounders.candidate)
labels_not_columns <- length(ev1) > 0 && any(!ev1 %in% names(df))
cat("Some evaluated labels are not columns of df:", labels_not_columns, "\n")
h1 <- if (length(ev1) > 0) run_hte(fit1) else list(res = NULL, warnings = character(0))
used1 <- if (is.null(h1$res)) character(0) else h1$res$covariates_used
ref1 <- test_hte(df$time, df$treat, as.matrix(df[, confs1]))
cat("Observed covariates_used :", if (length(used1)) used1 else "<none / NULL result>", "\n")
cat("Expected covariates      :", confs1, "\n")
if (!is.null(h1$res)) {
  cc <- h1$res$hte_test$constant_cate
  cat("Observed constant-CATE chi-sq / df:",
      if (is.null(cc)) "NULL" else paste(signif(cc$chi_sq, 6), "/", cc$df),
      "| source:", h1$res$source, "\n")
}
cat("Expected constant-CATE chi-sq / df (all three covariates):",
    signif(ref1$constant_cate$chi_sq, 6), "/", ref1$constant_cate$df, "\n")
cat("Warnings raised in case 1:", length(h1$warnings), "\n")
case1_dropped <- labels_not_columns && !all(c("age", "bmi") %in% used1)

## Case 2: continuous covariates only -> silent NULL
confs2 <- c("age", "bmi")
fit2 <- run_fs(confs2)
ev2 <- fit2$confounders.evaluated
cat("Case 2 confounders.evaluated:\n"); print(ev2)
h2 <- if (length(ev2) > 0) run_hte(fit2) else list(res = "not run", warnings = character(0))
cat("Case 2 result is NULL:", is.null(h2$res), "| warnings raised:", length(h2$warnings),
    "| expected: a test on age and bmi, or an error/warning\n")
case2_silent_null <- length(ev2) > 0 && is.null(h2$res) && length(h2$warnings) == 0

demonstrated <- isTRUE(case1_dropped) || isTRUE(case2_silent_null)
cat("case1_dropped:", case1_dropped, " case2_silent_null:", case2_silent_null, "\n")
cat("REPRO_VERDICT:", if (demonstrated) "demonstrated" else "not demonstrated", "\n")
```

```text
Case 1 class of fit: forestsearch/list 
Case 1 confounders.evaluated:
[1] "age <= 55" "age <= 47" "age <= 63" "bmi <= 27" "bmi <= 25" "bmi <= 30"
[7] "bm"       
Case 1 confounders.candidate:
[1] "q1" "q2" "q3" "q4" "q5" "q6" "q7"
Some evaluated labels are not columns of df: TRUE 
Observed covariates_used : bm 
Expected covariates      : age bmi bm 
Observed constant-CATE chi-sq / df: 0.113008 / 1 | source: confounders.evaluated (GRF/LASSO selected) 
Expected constant-CATE chi-sq / df (all three covariates): 1.30938 / 3 
Warnings raised in case 1: 0 
Case 2 confounders.evaluated:
[1] "age <= 55" "age <= 47" "age <= 63" "bmi <= 27" "bmi <= 25" "bmi <= 30"
Case 2 result is NULL: TRUE | warnings raised: 0 | expected: a test on age and bmi, or an error/warning
case1_dropped: TRUE  case2_silent_null: TRUE 
REPRO_VERDICT: demonstrated
```

### V042

- Finding: A-153 — The roxygen for c2 says the default is 1.0 on the ratio scale and 0.0 on the difference scale, but the formal default is a fixed 1.0, so a direct call with effe
- Runs: exit 0, 2.4 s. Verdict: demonstrated

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)

# Documented default, read from the roxygen source (read-only)
src <- readLines("/home/larryleon/Documents/GitHub/forestsearch/R/truefind_asymptotic_glm.R")
doc_says_zero <- any(grepl("Default: 1.0 for ratio, 0.0 for difference", src, fixed = TRUE))
formal_c2 <- formals(compute_detection_probability_glm)[["c2"]]

d <- d_eff_continuous(n_sg = 100, sigma_y = 1.5)   # 44.44
sd_s <- sqrt(8 / d)

# c2 omitted on the difference scale (relies on the default)
p_dflt <- compute_detection_probability_glm(
  theta = 0.5, d_eff = d, c1 = 0.2,
  effect_scale = "difference", method = "monte_carlo"
)
# documented default supplied explicitly
p_doc <- compute_detection_probability_glm(
  theta = 0.5, d_eff = d, c1 = 0.2, c2 = 0.0,
  effect_scale = "difference", method = "monte_carlo"
)
# what the omitted-c2 call actually used: c2 = 1.0 supplied explicitly (same seed, so identical)
p_one <- compute_detection_probability_glm(
  theta = 0.5, d_eff = d, c1 = 0.2, c2 = 1.0,
  effect_scale = "difference", method = "monte_carlo"
)
# closed form when the consistency threshold is 1.0 (screening is then implied)
p_closed <- stats::pnorm((1 - 0.5) / sd_s, lower.tail = FALSE)^2

cat("roxygen states c2 default 0.0 for difference:", doc_says_zero, "\n")
cat("formal default of c2:", format(formal_c2), "\n")
cat(sprintf("observed, c2 omitted               : %.4f\n", p_dflt))
cat(sprintf("explicit c2 = 1.0                  : %.4f\n", p_one))
cat(sprintf("closed form with c2 = 1.0          : %.4f\n", p_closed))
cat(sprintf("correct per documentation (c2 = 0) : %.4f\n", p_doc))

demonstrated <- isTRUE(doc_says_zero) &&
  isTRUE(all.equal(as.numeric(formal_c2), 1)) &&
  isTRUE(all.equal(p_dflt, p_one)) &&
  abs(p_dflt - p_closed) < 0.005 &&
  (p_doc - p_dflt) > 0.25
cat("REPRO_VERDICT:", if (demonstrated) "demonstrated" else "not demonstrated", "\n")
```

```text
roxygen states c2 default 0.0 for difference: TRUE 
formal default of c2: 1 
observed, c2 omitted               : 0.0143
explicit c2 = 1.0                  : 0.0143
closed form with c2 = 1.0          : 0.0142
correct per documentation (c2 = 0) : 0.7431
REPRO_VERDICT: demonstrated
```

### V043

- Finding: A-005 — get_Cox_sg() fits coxph() under suppressWarnings() with no per-arm event check, so a slice with a zero-event arm is expected to return a large finite log-HR tha
- Runs: exit 0, 2.3 s. Verdict: CONFIRMED finite divergent log-HR returned with no warning for a zero-event arm

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(3); n <- 60
d <- data.frame(tte = rexp(n), treat = rep(0:1, each = n / 2))
d$event <- ifelse(d$treat == 1, 0L, 1L)   # zero events in the treated arm
f <- build_cox_formula("tte", "event", "treat")
w <- character(0)
r <- tryCatch(withCallingHandlers(get_Cox_sg(d, cox.formula = f, treat.name = "treat"),
  warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") }), error = function(e) e)
print(r); cat("warnings surfaced:", length(w), "\n")
ok <- !inherits(r, "error") && length(r$est_obs) == 1L && is.finite(r$est_obs) && abs(r$est_obs) > 5 && length(w) == 0L
cat("REPRO_VERDICT:", if (ok) "CONFIRMED finite divergent log-HR returned with no warning for a zero-event arm" else "NOT REPRODUCED", "\n")
```

```text
$est_obs
[1] -21.21119

$se_obs
[1] 0.2967228

warnings surfaced: 0 
REPRO_VERDICT: CONFIRMED finite divergent log-HR returned with no warning for a zero-event arm
```

### V044

- Finding: INV-2 — The bootstrap reseeds the session RNG with set.seed(seed) (default 8316951L) in bootstrap_results() and .make_boot_index_matrix() and nothing restores the calle
- Runs: exit 0, 2.3 s. Verdict: CONFIRMED caller's RNG stream after the call is fixed by the bootstrap seed

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
df <- data.frame(id = 1:50)
set.seed(1); invisible(bootstrap_ystar(df, nb_boots = 3)); a <- runif(3)
set.seed(2); invisible(bootstrap_ystar(df, nb_boots = 3)); b <- runif(3)
print(rbind(a, b))
cat("REPRO_VERDICT:", if (identical(a, b)) "CONFIRMED caller's RNG stream after the call is fixed by the bootstrap seed" else "NOT REPRODUCED", "\n")
```

```text
       [,1]     [,2]     [,3]
a 0.5278842 0.469567 0.977497
b 0.5278842 0.469567 0.977497
REPRO_VERDICT: CONFIRMED caller's RNG stream after the call is fixed by the bootstrap seed
```

### V045

- Finding: A-007 — The four event counts for the original subgroup on each resample (events_H_0, events_H_1, events_Hc_0, events_Hc_1) are computed and then overwritten with NA be
- Runs: exit 0, 4.8 s. Verdict: CONFIRMED the four original-subgroup event-count columns are NA although the true counts are positive

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 200L; B <- 3L
d <- data.frame(id = seq_len(n), treat = rep(0:1, n / 2), time = rexp(n), event = rbinom(n, 1L, 0.8), x = rnorm(n))
d$treat.recommend <- as.integer(d$x > 0)
idx <- matrix(sample.int(n, n * B, replace = TRUE), B, n, byrow = TRUE)
# the inner forestsearch() call fails by construction (unused argument); the row is still built and the counts are computed before that call
fake <- list(confounders.candidate = character(0), args_call_all = list(not_a_forestsearch_argument = 1))
res <- as.data.frame(suppressWarnings(suppressMessages(bootstrap_results(fake, d, build_cox_formula("time", "event", "treat"),
  nb_boots = B, show_three = FALSE, H_obs = 0, Hc_obs = 0, boot_index_mat = idx))))
cols <- c("events_H_0", "events_H_1", "events_Hc_0", "events_Hc_1")
truth <- vapply(seq_len(B), function(b) { db <- d[idx[b, ], ]; sum(db$event[db$treat.recommend == 0 & db$treat == 0]) }, 0)
print(res[, c("boot_id", cols)]); cat("true events_H_0 per replicate:", truth, "\n")
ok <- nrow(res) == B && all(is.na(unlist(res[, cols]))) && all(truth > 0)
cat("REPRO_VERDICT:", if (ok) "CONFIRMED the four original-subgroup event-count columns are NA although the true counts are positive" else "NOT REPRODUCED", "\n")
```

```text
  boot_id events_H_0 events_H_1 events_Hc_0 events_Hc_1
1       1         NA         NA          NA          NA
2       2         NA         NA          NA          NA
3       3         NA         NA          NA          NA
true events_H_0 per replicate: 32 35 39 
REPRO_VERDICT: CONFIRMED the four original-subgroup event-count columns are NA although the true counts are positive
```

### V046

- Finding: A-008 — A bootstrap replicate whose inner forestsearch() call throws is recorded exactly like a replicate that found no subgroup (NA estimates, any_found = 0), its warn
- Runs: exit 0, 4.9 s. Verdict: DEFECT_PRESENT (replicate errors hidden, reported as no subgroup)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
future::plan(future::sequential)
set.seed(1); n <- 200L
d <- data.frame(id = seq_len(n), treat = rep(0:1, n / 2), time = rexp(n),
                event = rbinom(n, 1L, 0.8), x = rnorm(n))
d$treat.recommend <- as.integer(d$x > 0)
# stand-in fit: the unknown argument makes every replicate's forestsearch() call error inside try()
fake <- list(df.est = d, sg.harm = "{x <= 0}", confounders.candidate = character(0),
  args_call_all = list(outcome.name = "time", event.name = "event", treat.name = "treat",
    id.name = "id", est.scale = "hr", not_a_forestsearch_argument = 1,
    parallel_args = list(plan = "sequential", workers = 1L, show_message = FALSE)))
w <- character(0)
bc <- withCallingHandlers(
  suppressMessages(forestsearch_bootstrap_dofuture(fake, nb_boots = 3)),
  warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
inner <- tryCatch(do.call(forestsearch, c(list(df.analysis = d), fake$args_call_all)),
                  error = function(e) conditionMessage(e))
cat("inner error each replicate hit:", substr(as.character(inner)[1], 1, 120), "\n")
cat("warnings the caller saw:\n"); print(substr(w, 1, 160))
r <- as.data.frame(bc$results)
print(r[, c("boot_id", "any_found", "H_biasadj_2")])
ok <- is.character(inner) && nrow(r) == 3 && all(r$any_found == 0) &&
  !any(grepl("failed", w)) && any(grepl("re-identified", w))
cat("REPRO_VERDICT:", if (ok) "DEFECT_PRESENT (replicate errors hidden, reported as no subgroup)" else "NOT_REPRODUCED", "\n")
```

```text
inner error each replicate hit: unused argument (not_a_forestsearch_argument = 1) 
warnings the caller saw:
[1] "nb_boots < 100 may produce unreliable confidence intervals. Recommend at least 500 for production use."                                                          
[2] "No bootstrap resample re-identified a subgroup (0/3 succeeded), so no bias-corrected estimates could be computed. This usually means the selection criteria (e.g"
  boot_id any_found H_biasadj_2
1       1         0          NA
2       2         0          NA
3       3         0          NA
REPRO_VERDICT: DEFECT_PRESENT (replicate errors hidden, reported as no subgroup)
```

### V047

- Finding: A-009 — For GLM outcomes the bootstrap rebuilds its effect estimator from the stored call without adverse_outcome and without ps_adjust_method, so its observed and bias
- Runs: exit 0, 2.3 s. Verdict: DEFECT_PRESENT (bootstrap estimator omits adverse_outcome / ps_adjust_method; sign opposite to search scale)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 200L
d <- data.frame(treat = rep(0:1, n / 2)); d$y <- rnorm(n) + 0.5 * d$treat
# the search estimator for a default continuous fit (adverse_outcome resolved to FALSE)
e_search <- make_effect_estimator("continuous", "treat", "y", effect_measure = "MD", adverse_outcome = FALSE)(d)$estimate
# the estimator as rebuilt by the bootstrap (argument omitted -> TRUE)
e_boot <- make_effect_estimator("continuous", "treat", "y", effect_measure = "MD")(d)$estimate
src <- deparse(body(forestsearch_bootstrap_dofuture))
i <- grep("make_effect_estimator", src)[1]
call_txt <- paste(src[i:(i + 9)], collapse = " ")
cat("search-scale estimate:", e_search, " bootstrap-closure estimate:", e_boot, "\n")
cat("bootstrap call passes adverse_outcome:", grepl("adverse_outcome", call_txt),
    "| passes ps_adjust_method:", grepl("ps_adjust_method", call_txt), "\n")
ok <- isTRUE(all.equal(e_search, -e_boot)) && !grepl("adverse_outcome", call_txt) && !grepl("ps_adjust_method", call_txt)
cat("REPRO_VERDICT:", if (ok) "DEFECT_PRESENT (bootstrap estimator omits adverse_outcome / ps_adjust_method; sign opposite to search scale)" else "NOT_REPRODUCED", "\n")
```

```text
search-scale estimate: -0.4274624  bootstrap-closure estimate: 0.4274624 
bootstrap call passes adverse_outcome: FALSE | passes ps_adjust_method: FALSE 
REPRO_VERDICT: DEFECT_PRESENT (bootstrap estimator omits adverse_outcome / ps_adjust_method; sign opposite to search scale)
```

### V048

- Finding: INV-2 — The bootstrap entry point does not restore the caller's future plan on exit: its handler looks for an object .Last.future.plan that nothing in the package creat
- Runs: exit 0, 5.8 s. Verdict: DEFECT_PRESENT (caller's plan replaced by sequential)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 200L
d <- data.frame(id = seq_len(n), treat = rep(0:1, n / 2), time = rexp(n),
                event = rbinom(n, 1L, 0.8), x = rnorm(n))
d$treat.recommend <- as.integer(d$x > 0)
# stand-in fit, enough to pass the entry point; replicates fail fast inside try()
fake <- list(df.est = d, sg.harm = "{x <= 0}", confounders.candidate = character(0),
  args_call_all = list(outcome.name = "time", event.name = "event", treat.name = "treat",
    id.name = "id", est.scale = "hr", not_a_forestsearch_argument = 1,
    parallel_args = list(plan = "sequential", workers = 1L, show_message = FALSE)))
has_obj <- exists(".Last.future.plan")
future::plan(future::multisession, workers = 2)
before <- inherits(future::plan(), "multisession")
bc <- suppressWarnings(suppressMessages(forestsearch_bootstrap_dofuture(fake, nb_boots = 2)))
after_multi <- inherits(future::plan(), "multisession")
after_seq <- inherits(future::plan(), "sequential")
cat(".Last.future.plan exists:", has_obj, "| multisession before:", before,
    "| multisession after:", after_multi, "| sequential after:", after_seq, "\n")
future::plan(future::sequential)
cat("REPRO_VERDICT:", if (before && !after_multi && after_seq) "DEFECT_PRESENT (caller's plan replaced by sequential)" else "NOT_REPRODUCED", "\n")
```

```text
.Last.future.plan exists: FALSE | multisession before: TRUE | multisession after: FALSE | sequential after: TRUE 
REPRO_VERDICT: DEFECT_PRESENT (caller's plan replaced by sequential)
```

### V049

- Finding: A-010 — The success guard covers only zero re-identifying replicates; with exactly one, the infinitesimal-jackknife standard error is exactly 0 and a zero-width bias-co
- Runs: exit 0, 2.3 s. Verdict: DEFECT_PRESENT (zero SE, zero-width interval from one replicate)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
ystar <- bootstrap_ystar(data.frame(a = 1:60), nb_boots = 25, seed = 1L)
x <- rep(NA_real_, 25); x[7] <- 0.40   # exactly one replicate re-identified a subgroup
r <- as.data.frame(get_dfRes(Hobs = 0.5, seHobs = 0.2, H1_adj = x, H2_adj = x, ystar = ystar))
print(r)
ok <- isTRUE(r$sdH2 == 0) && isTRUE(r$H2_lower == r$H2_upper)
cat("REPRO_VERDICT:", if (ok) "DEFECT_PRESENT (zero SE, zero-width interval from one replicate)" else "NOT_REPRODUCED", "\n")
```

```text
        H0      sdH0 H0_lower H0_upper       H1 sdH1 H1_lower H1_upper       H2
1 1.648721 0.3297443 1.114056 2.439987 1.491825    0 1.491825 1.491825 1.491825
  sdH2 H2_lower H2_upper
1    0 1.491825 1.491825
REPRO_VERDICT: DEFECT_PRESENT (zero SE, zero-width interval from one replicate)
```

### V050

- Finding: A-012 — For a Cox model with more than one coefficient (a linear adjust_covariates term such as "x2") the treatment dfbeta column is looked up by column name, but survi
- Runs: exit 0, 2.3 s. Verdict: CONFIRMED (linear-adjusted Cox pieces are NULL; closed-form rate NA)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 200
d <- data.frame(Y = rexp(n), Event = rbinom(n, 1, 0.7), Treat = rep(0:1, n / 2), x2 = rnorm(n))
fit <- survival::coxph(survival::Surv(Y, Event) ~ Treat + x2, data = d)
cn <- colnames(residuals(fit, type = "dfbeta"))
cat("colnames of dfbeta residual matrix: "); print(cn)
u <- is.null(.consistency_cox_pieces(d)); a <- is.null(.consistency_cox_pieces(d, adjust_covariates = "x2"))
s <- is.null(.consistency_cox_pieces(d, adjust_covariates = "strata(I(x2 > 0))"))
rc <- consistency_resample(d, adjust_covariates = "x2", method = "closed")$rate_closed
cat("pieces NULL - unadjusted:", u, " linear-adjusted:", a, " strata-adjusted:", s, " closed rate (linear-adjusted):", rc, "\n")
cat("REPRO_VERDICT:", if (is.null(cn) && !u && a && !s && is.na(rc)) "CONFIRMED (linear-adjusted Cox pieces are NULL; closed-form rate NA)" else "NOT_CONFIRMED", "\n")
```

```text
colnames of dfbeta residual matrix: NULL
pieces NULL - unadjusted: FALSE  linear-adjusted: TRUE  strata-adjusted: FALSE  closed rate (linear-adjusted): NA 
REPRO_VERDICT: CONFIRMED (linear-adjusted Cox pieces are NULL; closed-form rate NA)
```

### V051

- Finding: A-013 — Cross-fit coefficients are combined with colMeans(beta_mat, na.rm = TRUE), so a step-2 coefficient that is aliased (NA) in one fold is silently replaced by the 
- Runs: exit 0, 2.4 s. Verdict: CONFIRMED (one-fold NA averaged away silently; all-fold NA -> NaN coefficient -> cryptic dina_subgroup error)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(2); n <- 300
folds <- .dina_make_folds(n, 2L, seed = 1L)
z <- numeric(n); z[folds$test[[1]][1:6]] <- 1          # non-zero only inside test fold 1
X <- cbind(X1 = runif(n, -1, 1), X2 = z)
W <- rbinom(n, 1, 0.5); Y <- X[, 1] + W * (0.5 + X[, 1]) + rnorm(n)
tr <- folds$train[[2]]; te <- folds$test[[2]]
r2 <- suppressWarnings(.dina_one_pass(X[tr, , drop = FALSE], Y[tr], W[tr], X[te, , drop = FALSE], Y[te], W[te],
        family = "gaussian", propensity_method = "logistic", baseline_method = "glm",
        cens_type = NULL, cens_params = list(), n_grid = 500L, eps = 1e-6))
fit <- suppressWarnings(dina_fit(X, Y, W, family = "gaussian", seed = 1L))
cat("fold-2 beta:", format(r2$beta), "\n"); print(coef(fit))
df <- data.frame(y = Y, w = W, x1 = X[, 1], x0 = 0)     # constant covariate: NA in every fold
fit0 <- suppressWarnings(dina(df, "y", "w", c("x1", "x0"), family = "gaussian", seed = 1L))
print(coef(fit0))
e <- try(dina_subgroup(fit0, df, c("x1", "x0"), m_diff = 0.5, n_min = 30L), silent = TRUE); print(e)
cat("REPRO_VERDICT:", if (is.na(r2$beta[3]) && is.finite(coef(fit)[3]) && is.nan(coef(fit0)[3]) && inherits(e, "try-error")) "CONFIRMED (one-fold NA averaged away silently; all-fold NA -> NaN coefficient -> cryptic dina_subgroup error)" else "NOT REPRODUCED", "\n")
```

```text
fold-2 beta: 0.394113 1.022762       NA 
(Intercept)          X1          X2 
  0.4613209   1.0982547   0.7756229 
(Intercept)          x1          x0 
  0.4770136   1.0991329         NaN 
[1] "Error in if (mean_tau < m_diff) next : \n  missing value where TRUE/FALSE needed\n"
attr(,"class")
[1] "try-error"
attr(,"condition")
<simpleError in if (mean_tau < m_diff) next: missing value where TRUE/FALSE needed>
REPRO_VERDICT: CONFIRMED (one-fold NA averaged away silently; all-fold NA -> NaN coefficient -> cryptic dina_subgroup error)
```

### V052

- Finding: A-015 — cens_type = "unknown" fits one pooled logistic model of status on X and returns the same function for both arms, so the two not-censoring probabilities cancel i
- Runs: exit 0, 2.4 s. Verdict: CONFIRMED (cens_type = 'unknown' numerically identical to NULL)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 400
X <- matrix(runif(n * 2, -1, 1), n, 2); W <- rbinom(n, 1, 0.5)
Tt <- rexp(n, exp(0.3 * X[, 1] + W * (0.2 + 0.5 * X[, 1])))
Cc <- rexp(n, 0.5 * exp(W))                      # arm-dependent censoring
Y <- data.frame(time = pmin(Tt, Cc), status = as.integer(Tt <= Cc))
f0 <- dina_fit(X, Y, W, family = "cox", cens_type = NULL, seed = 1)
f1 <- dina_fit(X, Y, W, family = "cox", cens_type = "unknown", seed = 1)
print(rbind(null = coef(f0), unknown = coef(f1)))
dd <- max(abs(coef(f0) - coef(f1))); cat("max |diff| =", dd, "\n")
cat("REPRO_VERDICT:", if (is.finite(dd) && dd < 1e-4) "CONFIRMED (cens_type = 'unknown' numerically identical to NULL)" else "NOT REPRODUCED", "\n")
```

```text
        (Intercept)        X1          X2
null     0.06412882 0.7126826 -0.06987324
unknown  0.06412916 0.7126826 -0.06987314
max |diff| = 3.442481e-07 
REPRO_VERDICT: CONFIRMED (cens_type = 'unknown' numerically identical to NULL)
```

### V053

- Finding: A-017 — dina_subgroup() and dina_frontier() check only that the fit has length(covariates) + 1 coefficients and then pair coefficients with columns by position, so a `c
- Runs: exit 0, 2.4 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 300
df <- data.frame(x1 = runif(n, -1, 1), x2 = runif(n, -1, 1), x3 = runif(n, -1, 1))
fit <- structure(list(coefficients = c("(Intercept)" = 0, x1 = 1, x2 = 0, x3 = 0),
                      vcov = diag(4) * 0.01, family = "gaussian"), class = "dina")
a <- dina_subgroup(fit, df, covariates = c("x1", "x2", "x3"), m_diff = 0.3, n_min = 30L, max_depth = 1L)
b <- tryCatch(dina_subgroup(fit, df, covariates = c("x3", "x2", "x1"), m_diff = 0.3, n_min = 30L, max_depth = 1L),
              error = function(e) e)
cat("coef names:", names(coef(fit)), "\n")
cat("training order:", a$label, a$n_subgroup, "\n")
cat("permuted order:", if (inherits(b, "error")) paste("ERROR", conditionMessage(b)) else paste(b$label, b$n_subgroup), "\n")
cat("REPRO_VERDICT:", if (!inherits(b, "error") && isTRUE(b$found) && !identical(a$mask, b$mask)) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
coef names: (Intercept) x1 x2 x3 
training order: {x1 >= -0.321854124311358} 196 
permuted order: {x3 >= -0.371596340090036} 212 
REPRO_VERDICT: CONFIRMED
```

### V054

- Finding: A-018 — The returned `labels`/`label` render the threshold with paste0() on the raw double (15 significant digits), so evaluating the documented label with get_dfpred()
- Runs: exit 0, 2.2 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 300
df <- data.frame(x1 = (3 * sample.int(200, n, replace = TRUE) - 1) / 3, x2 = rnorm(n))  # values k - 1/3
fit <- structure(list(coefficients = c("(Intercept)" = 0, x1 = 1, x2 = 0),
                      vcov = diag(3) * 0.01, family = "gaussian"), class = "dina")
sg <- dina_subgroup(fit, df, covariates = c("x1", "x2"), m_diff = 150, n_min = 30L,
                    direction = "right", max_depth = 1L, sg_focus = "maxSG")
cat("threshold (17g):", formatC(sg$threshold, digits = 17, format = "g"), " label:", sg$labels, "\n")
n_lab <- sum(get_dfpred(df, sg$labels)$treat.recommend == 0L)
cat("n from mask:", sum(sg$mask), " n_subgroup:", sg$n_subgroup, " n from get_dfpred(labels):", n_lab, "\n")
cat("REPRO_VERDICT:", if (isTRUE(sg$found) && n_lab < sum(sg$mask)) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
threshold (17g): 110.66666666666667  label: {x1 >= 110.666666666667} 
n from mask: 113  n_subgroup: 113  n from get_dfpred(labels): 110 
REPRO_VERDICT: CONFIRMED
```

### V055

- Finding: A-019 — The depth-1 enumeration includes the non-cuts `x <= max(x)` and `x >= min(x)` (the whole sample, empty complement) for every covariate, so dina_subgroup() can r
- Runs: exit 0, 2.3 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
fit <- structure(list(coefficients = c("(Intercept)" = 1, x = 0.01),
                      vcov = diag(2) * 0.01, family = "gaussian"), class = "dina")
df <- data.frame(x = as.numeric(1:100))
sg <- dina_subgroup(fit, df, "x", m_diff = 0.5, n_min = 10L, max_depth = 1L, sg_focus = "maxSG")
cat(sg$label, " n_subgroup =", sg$n_subgroup, " n_total =", sg$n_total, " all in mask:", all(sg$mask), "\n")
fr <- dina_frontier(fit, df, "x", n_min = 10L)
full <- fr[fr$n_subgroup == nrow(df), c("direction", "threshold", "n_subgroup", "cut_expr")]
print(full)
cat("REPRO_VERDICT:", if (isTRUE(sg$found) && all(sg$mask) && sg$n_subgroup == sg$n_total && nrow(full) >= 1L) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
{x <= 100}  n_subgroup = 100  n_total = 100  all in mask: TRUE 
   direction threshold n_subgroup cut_expr
91      left       100        100 x <= 100
92     right         1        100   x >= 1
REPRO_VERDICT: CONFIRMED
```

### V056

- Finding: A-020 — The harm orientation (tau_sign) that forestsearch applies to DINA for continuous/binary outcomes with adverse_outcome = FALSE exists only on the subgroup_method
- Runs: exit 0, 2.4 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
has_sg <- "tau_sign" %in% names(formals(dina_subgroup)); has_fr <- "tau_sign" %in% names(formals(dina_frontier))
in_main <- any(grepl("tau_sign", deparse(body(forestsearch)), fixed = TRUE))
cat("tau_sign formal: dina_subgroup", has_sg, " dina_frontier", has_fr, " mentioned in forestsearch() body:", in_main, "\n")
fit <- structure(list(coefficients = c("(Intercept)" = 0, x = -5), vcov = diag(c(1, 0.01)), family = "gaussian"), class = "dina")
df <- data.frame(x = as.numeric(1:20))
sg <- dina_subgroup(fit, df, "x", m_diff = 30, n_min = 5L, max_depth = 1L, tau_sign = -1)
sg0 <- dina_subgroup(fit, df, "x", m_diff = 30, n_min = 5L, max_depth = 1L)   # what the use_dina screening call does
fr <- dina_frontier(fit, df, "x", scope = "harm", m_diff = 30, n_min = 5L)
cat("harm-oriented candidates:", nrow(sg$candidates), " un-oriented found:", sg0$found, " harm-scope frontier rows:", nrow(fr), "\n")
cat("REPRO_VERDICT:", if (has_sg && !has_fr && !in_main && isTRUE(sg$found) && !isTRUE(sg0$found) && nrow(fr) == 0L) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
tau_sign formal: dina_subgroup TRUE  dina_frontier FALSE  mentioned in forestsearch() body: FALSE 
harm-oriented candidates: 26  un-oriented found: FALSE  harm-scope frontier rows: 0 
REPRO_VERDICT: CONFIRMED
```

### V057

- Finding: A-021 — The exported get_dfpred() called with sg.harm = NULL or character(0), which is what a no-subgroup fit carries, silently flags every row as a harm-subgroup membe
- Runs: exit 0, 2.3 s. Verdict: CONFIRMED (no subgroup definition -> every row flagged as harm-subgroup member, no error)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
df <- data.frame(x = 1:5)
r1 <- tryCatch(get_dfpred(df, NULL)$treat.recommend, error = function(e) e)
r2 <- tryCatch(get_dfpred(df, character(0))$treat.recommend, error = function(e) e)
print(r1); print(r2)
cat("REPRO_VERDICT:", if (identical(r1, rep(0L, 5)) && identical(r2, rep(0L, 5)))
  "CONFIRMED (no subgroup definition -> every row flagged as harm-subgroup member, no error)" else "REFUTED", "\n")
```

```text
[1] 0 0 0 0 0
[1] 0 0 0 0 0
REPRO_VERDICT: CONFIRMED (no subgroup definition -> every row flagged as harm-subgroup member, no error)
```

### V058

- Finding: A-022 — evaluate_comparison() coerces numeric-level factors but not numeric-valued character columns, so a cut label on such a column is compared as strings, unlike the
- Runs: exit 0, 2.3 s. Verdict: CONFIRMED (character digits compared as strings; factor and DINA coding compare numerically)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
d_chr <- data.frame(x = c("2", "10", "9"), stringsAsFactors = FALSE)
d_fac <- data.frame(x = factor(c("2", "10", "9")))
a <- evaluate_comparison("x >= 9", d_chr)                 # character column
b <- evaluate_comparison("x >= 9", d_fac)                 # factor column, coerced
k <- .coerce_covariates_numeric(d_chr, "x")$x >= 9        # what DINA selects on
print(a); print(b); print(k)
print(get_dfpred(d_chr, "{x >= 9}")$treat.recommend)
cat("REPRO_VERDICT:", if (!identical(a, k) && identical(b, k))
  "CONFIRMED (character digits compared as strings; factor and DINA coding compare numerically)" else "REFUTED", "\n")
```

```text
[1] FALSE FALSE  TRUE
[1] FALSE  TRUE  TRUE
[1] FALSE  TRUE  TRUE
[1] 1 1 0
REPRO_VERDICT: CONFIRMED (character digits compared as strings; factor and DINA coding compare numerically)
```

### V059

- Finding: A-023 — Effect-based re-selection overwrites only some winner fields, so the stored identifier object (grp.consistency$out_sg / fs$grf_res) mixes the new winner's defin
- Runs: exit 0, 4.6 s. Verdict: DEFECT_REPRODUCED (out_sg label / CI belong to a different subgroup than its signature and mean)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
future::plan("sequential")
shown <- FALSE; tried <- 0L
for (s in 1:6) {
  set.seed(s); n <- 400
  d <- data.frame(id = seq_len(n), treat = rbinom(n, 1, 0.5), x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
  lp <- d$treat * (0.9 * (d$x1 > 0) - 0.3 * (d$x1 <= 0)) + 0.3 * d$x2
  t0 <- rexp(n, 0.1 * exp(lp)); cc <- runif(n, 5, 30)
  d$tte <- pmin(t0, cc); d$event <- as.integer(t0 <= cc)
  fs <- suppressMessages(suppressWarnings(forestsearch(df.analysis = d, confounders.name = c("x1", "x2", "x3"),
    subgroup_method = "dina", sg_focus = "hr", n.min = 60, details = FALSE, quiet = TRUE,
    parallel_args = list(plan = "sequential", workers = 1L, show_message = FALSE))))
  sg <- fs$grp.consistency$out_sg
  if (is.null(sg)) next
  tried <- tried + 1L
  ops <- ifelse(sg$direction == "left", "<=", ">=")
  lab_now <- paste0("{", paste0(sg$covariate, " ", ops, " ", sg$threshold, collapse = " & "), "}")
  cat("seed", s, "| fs$sg.harm:", paste(fs$sg.harm, collapse = " & "), "| out_sg$label:", sg$label,
      "| mean_tau_hat:", sg$mean_tau_hat, "| CI midpoint:", mean(sg$ci), "\n")
  if (!identical(lab_now, sg$label) || abs(mean(sg$ci) - sg$mean_tau_hat) > 1e-8) { shown <- TRUE; print(sg); break }
}
cat("REPRO_VERDICT:", if (shown) "DEFECT_REPRODUCED (out_sg label / CI belong to a different subgroup than its signature and mean)" else sprintf("NOT_TRIGGERED (effect winner equalled the native winner in all %d fits with a subgroup)", tried), "\n")
```

```text
seed 1 | fs$sg.harm: {x1 >= 0.89234145446392354} & {x3 >= -0.89946754333462919} | out_sg$label: {x1 >= 1.08342991008328} | mean_tau_hat: 1.133083 | CI midpoint: 1.135233 
DINA subgroup search
  Family:      cox
  Focus:       eff (most extreme effect)
  m_diff:      0.2231
  n_min:       60
  Candidates:  3372 searched, 766 satisfied m_diff and size

  Signature:   x1 >=  0.8923 & x3 >= -0.8995
  Subgroup n:  66 / 400 (16.5%)
  mean tau-hat:1.133 (SE 0.2004)
  Wald 95% CI:  [0.7425, 1.528]

  Note: CI is conditional on the selected (covariate, threshold)
  and does NOT adjust for the selection across candidates.
REPRO_VERDICT: DEFECT_REPRODUCED (out_sg label / CI belong to a different subgroup than its signature and mean)
```

### V060

- Finding: A-024 — For count outcomes forestsearch() requires offset.name for IRR, but the DINA fit is called without exposure and dina()/dina_fit() have no offset argument, so DI
- Runs: exit 0, 2.5 s. Verdict: DEFECT_REPRODUCED (tau-hat is the count ratio; a true IRR of 2 fails the log(1.25) floor)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 4000
d <- data.frame(treat = rbinom(n, 1, 0.5), x1 = rnorm(n))
d$fu <- ifelse(d$treat == 1, 0.5, 1)                      # treated arm followed half as long
d$y  <- rpois(n, d$fu * exp(-0.2 + log(2) * d$treat))     # true rate ratio = 2 everywhere
fit  <- dina(df = d, outcome = "y", treatment = "treat", covariates = "x1", family = "poisson", seed = 1L)
tau_bar <- mean(coef(fit)[1] + coef(fit)[2] * d$x1)
log_irr <- unname(coef(glm(y ~ treat + offset(log(fu)), family = poisson, data = d))["treat"])
has_off <- any(c("offset", "offset.name", "exposure") %in% c(names(formals(dina)), names(formals(dina_fit))))
cat("DINA mean tau-hat:", tau_bar, "| offset-adjusted log IRR:", log_irr, "| floor log(1.25):", log(1.25), "| dina has an offset formal:", has_off, "\n")
cat("REPRO_VERDICT:", if (!has_off && abs(tau_bar) < 0.15 && log_irr > 0.5) "DEFECT_REPRODUCED (tau-hat is the count ratio; a true IRR of 2 fails the log(1.25) floor)" else "NOT_REPRODUCED", "\n")
```

```text
DINA mean tau-hat: 0.00816948 | offset-adjusted log IRR: 0.6883667 | floor log(1.25): 0.2231436 | dina has an offset formal: FALSE 
REPRO_VERDICT: DEFECT_REPRODUCED (tau-hat is the count ratio; a true IRR of 2 fails the log(1.25) floor)
```

### V061

- Finding: A-025 — When effect-based re-selection scores no candidate or errors, the native-statistic winner is returned silently with no effect floor applied; this happens determ
- Runs: exit 0, 4.0 s. Verdict: DEFECT_REPRODUCED (subgroup returned below the effect floor; no candidate was scored on the effect)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
future::plan("sequential")
set.seed(5); n <- 500
d <- data.frame(id = seq_len(n), treat = rbinom(n, 1, 0.5), x1 = rnorm(n), x2 = rnorm(n), fu = round(runif(n, 1, 5), 2))
d$y <- rpois(n, d$fu * 0.2 * exp(d$treat * ifelse(d$x1 > 0, log(2), 0)))   # rates 0.2-0.4 per unit time
run <- function(em, thr) tryCatch(suppressMessages(suppressWarnings(forestsearch(
  df.analysis = d, outcome.name = "y", treat.name = "treat", id.name = "id", confounders.name = c("x1", "x2"),
  outcome_type = "count", effect_measure = em, offset.name = "fu", effect.threshold = thr,
  subgroup_method = "grf", n.min = 60, details = FALSE, quiet = TRUE,
  parallel_args = list(plan = "sequential", workers = 1L, show_message = FALSE)))),
  error = function(e) { cat("ERROR (", em, "):", conditionMessage(e), "\n"); NULL })
a <- run("IRD", 0.99)   # unreachable floor: no subgroup can have a rate difference >= 0.99
b <- run("IRR", 50)     # control: unreachable ratio floor on a scorable estimand
ok <- FALSE
if (!is.null(a) && !is.null(a$sg.harm)) {
  e <- a$df.est; m <- e$treat.recommend == 0
  ird <- sum(e$y[m & e$treat == 1]) / sum(e$fu[m & e$treat == 1]) - sum(e$y[m & e$treat == 0]) / sum(e$fu[m & e$treat == 0])
  se <- a$grf_res$candidates$sel_effect
  cat("IRD run: sg.harm =", paste(a$sg.harm, collapse = " & "), "| crude IRD in it =", ird, "| effect floor =", a$admission$effect_floor,
      "| candidates scored:", sum(is.finite(se)), "of", length(se), "| grf_res$select_statistic:", deparse(a$grf_res$select_statistic), "\n")
  ok <- !any(is.finite(se)) && ird < a$admission$effect_floor
} else cat("IRD run returned no subgroup\n")
cat("IRR control run: sg.harm =", if (is.null(b) || is.null(b$sg.harm)) "NULL" else paste(b$sg.harm, collapse = " & "), "\n")
cat("REPRO_VERDICT:", if (ok) "DEFECT_REPRODUCED (subgroup returned below the effect floor; no candidate was scored on the effect)" else "NOT_REPRODUCED", "\n")
```

```text
IRD run: sg.harm = {x1 <= -0.51493997018469095} | crude IRD in it = -0.01730942 | effect floor = 0.99 | candidates scored: 0 of 81 | grf_res$select_statistic: NULL 
IRR control run: sg.harm = NULL 
REPRO_VERDICT: DEFECT_REPRODUCED (subgroup returned below the effect floor; no candidate was scored on the effect)
```

### V062

- Finding: A-062 — Under the default grf_select_statistic = "effect", dmin.grf does not restrict the candidate family that is ranked on the effect, but it still vetoes the whole s
- Runs: exit 0, 4.7 s. Verdict: DEFECT_REPRODUCED (dmin.grf vetoes effect re-selection)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
future::plan("sequential")
set.seed(3); n <- 500
d <- data.frame(id = seq_len(n), treat = rbinom(n, 1, 0.5), x1 = rnorm(n), x2 = rnorm(n), x3 = rnorm(n))
lp <- d$treat * (0.9 * (d$x1 > 0) - 0.3 * (d$x1 <= 0))
t0 <- rexp(n, 0.1 * exp(lp)); cc <- runif(n, 5, 30)
d$tte <- pmin(t0, cc); d$event <- as.integer(t0 <= cc)
run <- function(dm) suppressMessages(suppressWarnings(forestsearch(df.analysis = d, confounders.name = c("x1", "x2", "x3"),
  subgroup_method = "grf", dmin.grf = dm, n.min = 60, details = FALSE, quiet = TRUE,
  parallel_args = list(plan = "sequential", workers = 1L, show_message = FALSE))))
A <- run(0)
cand <- A$grf_res$candidates
gate <- FALSE; ignored <- NA
if (!is.null(A$sg.harm) && !is.null(cand)) {
  dr_max <- max(cand$effect); fin <- which(is.finite(cand$sel_effect))
  w <- fin[which.max(cand$sel_effect[fin])]; dr_w <- cand$effect[w]
  B <- run(dr_max + 1)                      # no DR candidate clears dmin.grf
  cat("A (dmin.grf = 0): sg.harm =", paste(A$sg.harm, collapse = " & "), "| candidates above the effect floor:", A$grf_res$admitted_n, "\n")
  cat("B (dmin.grf = max DR effect + 1): sg.harm =", if (is.null(B$sg.harm)) "NULL" else paste(B$sg.harm, collapse = " & "), "\n")
  gate <- is.null(B$sg.harm) && isTRUE(A$grf_res$admitted_n > 0)
  if (dr_w < dr_max) {
    dmid <- (dr_w + dr_max) / 2; C <- run(dmid)   # winner's own DR effect is below dmin.grf
    ignored <- identical(C$sg.harm, A$sg.harm)
    cat("C (dmin.grf =", dmid, "> winner's DR effect", dr_w, "): sg.harm =", paste(C$sg.harm, collapse = " & "), "| same as A:", ignored, "\n")
  }
} else cat("A found no subgroup; change the seed\n")
cat("REPRO_VERDICT:", if (gate) paste0("DEFECT_REPRODUCED (dmin.grf vetoes effect re-selection", if (isTRUE(ignored)) "; and is ignored for the ranked family when the veto does not fire)" else ")") else "NOT_REPRODUCED", "\n")
```

```text
A (dmin.grf = 0): sg.harm = {x1 > 0.19540177435839479} & {x3 > 0.2712408755590815} | candidates above the effect floor: 66 
B (dmin.grf = max DR effect + 1): sg.harm = NULL 
REPRO_VERDICT: DEFECT_REPRODUCED (dmin.grf vetoes effect re-selection)
```

### V063

- Finding: INV-3 — Invariant 3 is violated: the formal default is use_twostage = TRUE while the roxygen says 'Default FALSE for backward compatibility', and two documented mr_infe
- Runs: exit 0, 2.2 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
rd <- paste(readLines("/home/larryleon/Documents/GitHub/forestsearch/man/forestsearch.Rd"), collapse = "\n")
f_two <- formals(forestsearch)$use_twostage
doc_two_false <- grepl("Default FALSE for backward compatibility", rd, fixed = TRUE)
f_fc <- formals(fs_mr_inference)$field_complement
doc_fc_false <- grepl("field_complement}}{Logical, default \\code{FALSE}", rd, fixed = TRUE)
ij <- grepl('ci_method, "ij"', paste(deparse(body(forestsearch:::.fs_apply_mr)), collapse = " "), fixed = TRUE)
cat("formal use_twostage:", deparse(f_two), "| Rd says Default FALSE:", doc_two_false, "\n")
cat("fs_mr_inference field_complement formal:", deparse(f_fc), "| Rd says default FALSE:", doc_fc_false, "\n")
cat("DINA/GRF wrapper ci_method fallback is ij:", ij, "\n")
cat("REPRO_VERDICT:", if (isTRUE(f_two) && doc_two_false) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
formal use_twostage: TRUE | Rd says Default FALSE: TRUE 
fs_mr_inference field_complement formal: TRUE | Rd says default FALSE: TRUE 
DINA/GRF wrapper ci_method fallback is ij: FALSE 
REPRO_VERDICT: CONFIRMED
```

### V064

- Finding: A-026 — The documented tune_grf argument of forestsearch() is never read in the function body or forwarded to any GRF fit, so tune_grf = TRUE silently changes nothing o
- Runs: exit 0, 2.2 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
in_formals <- "tune_grf" %in% names(formals(forestsearch))
in_body <- "tune_grf" %in% all.names(body(forestsearch))
b1 <- "tune_grf" %in% names(formals(forestsearch:::.build_grf_survival_args))
b2 <- "tune_grf" %in% names(formals(forestsearch:::.build_grf_glm_args))
b3 <- "tune_grf" %in% names(formals(forestsearch:::.forestsearch_grf_select))
cat("formal of forestsearch():", in_formals, "\n")
cat("referenced in forestsearch() body:", in_body, "\n")
cat("formal of survival builder / glm builder / grf_select:", b1, b2, b3, "\n")
cat("accepted by grf.subg.harm.survival:", "tune_grf" %in% names(formals(grf.subg.harm.survival)), "\n")
cat("REPRO_VERDICT:", if (in_formals && !in_body && !b1 && !b2 && !b3) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
formal of forestsearch(): TRUE 
referenced in forestsearch() body: FALSE 
formal of survival builder / glm builder / grf_select: FALSE FALSE FALSE 
accepted by grf.subg.harm.survival: TRUE 
REPRO_VERDICT: CONFIRMED
```

### V065

- Finding: A-026 — forestsearch(tune_grf = TRUE) is accepted and silently ignored: the formal is never read and neither GRF argument builder carries it, so the inner GRF fit alway
- Runs: exit 0, 2.3 s. Verdict: CONFIRMED (tune_grf formal is never forwarded)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
is_formal <- "tune_grf" %in% names(formals(forestsearch))
body_refs <- sum(grepl("tune_grf", deparse(body(forestsearch)), fixed = TRUE))
helper_refs <- sum(grepl("tune_grf", deparse(body(.forestsearch_grf_select)), fixed = TRUE))
in_surv <- "tune_grf" %in% names(formals(.build_grf_survival_args))
in_glm  <- "tune_grf" %in% names(formals(.build_grf_glm_args))
print(c(is_formal = is_formal, body_refs = body_refs, grf_select_refs = helper_refs, in_survival_builder = in_surv, in_glm_builder = in_glm))
ok <- is_formal && body_refs == 0 && helper_refs == 0 && !in_surv && !in_glm
cat("REPRO_VERDICT:", if (ok) "CONFIRMED (tune_grf formal is never forwarded)" else "NOT_REPRODUCED", "\n")
```

```text
          is_formal           body_refs     grf_select_refs in_survival_builder 
                  1                   0                   0                   0 
     in_glm_builder 
                  0 
REPRO_VERDICT: CONFIRMED (tune_grf formal is never forwarded)
```

### V066

- Finding: A-053 — effect_measure is not validated against outcome_type for continuous outcomes, so outcome_type = "continuous" with e.g. effect_measure = "OR" or "RD" silently ru
- Runs: exit 0, 2.8 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
future::plan(future::sequential)
est <- try(make_effect_estimator(outcome_type = "continuous", treat.name = "treat", outcome.name = "y", effect_measure = "OR"), silent = TRUE)
cat("estimator built for continuous + OR without error:", is.function(est), "\n")
set.seed(6); n <- 400
d <- data.frame(id = seq_len(n), x1 = rbinom(n, 1, 0.5), x2 = rbinom(n, 1, 0.5), treat = rbinom(n, 1, 0.5))
d$y <- 1.5 * d$treat * d$x1 + rnorm(n)
fs <- try(suppressWarnings(suppressMessages(forestsearch(df.analysis = d, outcome.name = "y", treat.name = "treat", id.name = "id",
  confounders.name = c("x1", "x2"), outcome_type = "continuous", effect_measure = "OR", adverse_outcome = TRUE, maxk = 1,
  quiet = TRUE, parallel_args = list(plan = "sequential", workers = 1L, show_message = FALSE)))), silent = TRUE)
ok <- !inherits(fs, "try-error")
cat("forestsearch() errored:", !ok, if (!ok) conditionMessage(attr(fs, "condition")) else "", "\n")
tc <- if (ok) fs$threshold_config else NULL
if (ok) str(tc[c("outcome_type", "effect_measure", "scale", "screening")])
cat("REPRO_VERDICT:", if (is.function(est) && ok && identical(tc$effect_measure, "OR") && identical(tc$scale, "log")) "CONFIRMED" else if (is.function(est)) "PARTIAL (estimator accepts the mismatch; full run did not complete)" else "NOT_REPRODUCED", "\n")
```

```text
estimator built for continuous + OR without error: TRUE 
forestsearch() errored: FALSE  
List of 4
 $ outcome_type  : chr "continuous"
 $ effect_measure: chr "OR"
 $ scale         : chr "log"
 $ screening     : num 0.223
REPRO_VERDICT: CONFIRMED
```

### V067

- Finding: A-027 — The overdispersion argument is forwarded only to the GRF routines and never reaches the effect estimator, so on the default subgroup_method = "consistency" path
- Runs: exit 0, 2.3 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
in_formals <- "overdispersion" %in% names(formals(forestsearch))
in_est <- "overdispersion" %in% names(formals(make_effect_estimator))
in_pois <- "overdispersion" %in% names(formals(forestsearch:::.make_poisson_rate_estimator))
src <- deparse(body(forestsearch))
i <- grep("make_effect_estimator(", src, fixed = TRUE)
passed <- vapply(i, function(k) any(grepl("overdispersion", src[k:min(length(src), k + 10L)], fixed = TRUE)), logical(1))
cat("forestsearch() formal:", in_formals, "\n")
cat("make_effect_estimator formal:", in_est, "| Poisson estimator formal:", in_pois, "\n")
cat("make_effect_estimator() calls in forestsearch():", length(i), "| any passes overdispersion:", any(passed), "\n")
cat("REPRO_VERDICT:", if (in_formals && !in_est && !in_pois && length(i) > 0 && !any(passed)) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
forestsearch() formal: TRUE 
make_effect_estimator formal: FALSE | Poisson estimator formal: FALSE 
make_effect_estimator() calls in forestsearch(): 2 | any passes overdispersion: FALSE 
REPRO_VERDICT: CONFIRMED
```

### V068

- Finding: A-028 — A user-supplied effect or consistency threshold above 1 is replaced by the package default for effect_measure = "IRD" as well as "RD" on the ground that both ar
- Runs: exit 0, 3.5 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
future::plan(future::sequential)
set.seed(1); n <- 300
d <- data.frame(id = seq_len(n), x1 = rbinom(n, 1, 0.5), x2 = rbinom(n, 1, 0.5), treat = rbinom(n, 1, 0.5), fu = runif(n, 0.5, 2))
d$y <- rpois(n, d$fu * (1 + 2 * d$treat * d$x1))
w <- character(0)
fs <- try(withCallingHandlers(suppressMessages(forestsearch(df.analysis = d, outcome.name = "y", treat.name = "treat", id.name = "id",
  confounders.name = c("x1", "x2"), outcome_type = "count", effect_measure = "IRD", offset.name = "fu", effect.threshold = 1.5,
  maxk = 1, quiet = TRUE, parallel_args = list(plan = "sequential", workers = 1L, show_message = FALSE))),
  warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") }), silent = TRUE)
remap <- grep("Remapping to default", w, value = TRUE)
cat(remap, sep = "\n")
used <- if (!inherits(fs, "try-error")) fs$threshold_config$screening else NA
cat("supplied effect.threshold = 1.5; screening threshold used:", used, "\n")
cat("REPRO_VERDICT:", if (length(remap) > 0 && any(grepl("'IRD'", remap, fixed = TRUE)) && any(grepl("Remapping to default: 0.05", remap, fixed = TRUE))) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
effect.threshold = 1.50 appears to be on a ratio scale, but effect_measure = 'IRD' uses the identity scale (range [-1, 1]).
  Did you mean effect.threshold = 1.50 on the OR scale?
  Remapping to default: 0.05.
  To use a custom IRD threshold, pass a value in [-1, 1] (e.g., effect.threshold = 0.07 for a 7 pct-point threshold).
supplied effect.threshold = 1.5; screening threshold used: 0.01 
REPRO_VERDICT: CONFIRMED
```

### V069

- Finding: A-029 — A user-supplied ps_hat is silently discarded, unvalidated, whenever ps_method resolves to "none", which is the default for is.RCT = TRUE, so a call supplying pr
- Runs: exit 0, 2.8 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
future::plan(future::sequential)
set.seed(505); n <- 400
df <- data.frame(id = seq_len(n), x1 = rbinom(n, 1, 0.5), x2 = rbinom(n, 1, 0.5), treat = rbinom(n, 1, 0.5))
df$y <- rbinom(n, 1, plogis(-1 + 1.2 * df$treat * df$x1))
cond <- character(0)
a <- try(withCallingHandlers(forestsearch(df.analysis = df, outcome.name = "y", treat.name = "treat", id.name = "id",
  confounders.name = c("x1", "x2"), outcome_type = "binary", n.min = 40, d0.min = 5, d1.min = 5, maxk = 1, quiet = TRUE,
  ps_hat = runif(n, 0.05, 0.95), ps_adjust_method = "iptw",
  parallel_args = list(plan = "sequential", workers = 1L, show_message = FALSE)),
  warning = function(cnd) { cond <<- c(cond, conditionMessage(cnd)); invokeRestart("muffleWarning") },
  message = function(cnd) { cond <<- c(cond, conditionMessage(cnd)); invokeRestart("muffleMessage") }), silent = TRUE)
ok <- !inherits(a, "try-error")
n_ps <- sum(grepl("ps_hat|ps_adjust|propensity", cond, ignore.case = TRUE))
cat("run completed:", ok, "| warnings/messages mentioning the PS inputs:", n_ps, "\n")
if (ok) cat("recorded ps_method / ps_adjust_method:", a$args_call_all$ps_method, "/", a$args_call_all$ps_adjust_method, "\n")
cat("REPRO_VERDICT:", if (ok && n_ps == 0 && identical(a$args_call_all$ps_method, "none") && identical(a$args_call_all$ps_adjust_method, "none")) "CONFIRMED" else if (!ok) "INCONCLUSIVE (run failed)" else "NOT_REPRODUCED", "\n")
```

```text
run completed: TRUE | warnings/messages mentioning the PS inputs: 0 
recorded ps_method / ps_adjust_method: none / none 
REPRO_VERDICT: CONFIRMED
```

### V070

- Finding: A-031 — On the subgroup_method = "dina" and "grf" branches mr_harm_confirmed is isTRUE(harm_flag), so an MR result whose harm_flag is NA (selected subgroup could not be
- Runs: exit 0, 2.3 s. Verdict: confirmed

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 200
d <- data.frame(Y = rexp(n), Event = rbinom(n, 1, 0.7), Treat = rbinom(n, 1, 0.5), x = rnorm(n))
spec <- list(outcome_type = "survival", effect_measure = "HR", treat.name = "Treat", outcome.name = "Y",
             event.name = "Event", offset.name = NULL, adjust_covariates = NULL, adverse_outcome = TRUE)
# selected subgroup of 4 rows cannot be fitted (fs_mr_inference.R:95): the short variant
r <- .fs_apply_mr(d, list(a = which(d$x > 0)), 1:4, spec,
                  admission = list(effect_floor = log(0.5), consistency = NULL),
                  effect_neighborhood = 0.1, reselection_default = "maxeff",
                  mr_inference_args = list(draws = 100L), seedit = 1L)
str(r[c("selected_index", "harm_flag", "note")])
dina_grf <- if (!is.null(r)) isTRUE(r$harm_flag) else NA   # forestsearch_main.R:2875 / 3060
cons     <- if (!is.null(r)) r$harm_flag else NA           # forestsearch_main.R:3956
cat("DINA/GRF expression:", dina_grf, "| consistency expression:", cons, "\n")
cat("REPRO_VERDICT:", if (!is.null(r) && length(r$harm_flag) == 1L && is.na(r$harm_flag) && identical(dina_grf, FALSE) && is.na(cons)) "confirmed" else "not reproduced", "\n")
```

```text
List of 3
 $ selected_index: int NA
 $ harm_flag     : logi NA
 $ note          : chr "selected subgroup could not be fit in the reconstructed family"
DINA/GRF expression: FALSE | consistency expression: NA 
REPRO_VERDICT: confirmed
```

### V071

- Finding: A-032 — With use_grf = TRUE a supplied grf_res that already has tree.cuts is never read for its cuts (no GRF cut enters the candidate family, with no message), and a su
- Runs: exit 0, 4.5 s. Verdict: confirmed | supplied grf_cuts overwritten: TRUE

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
future::plan(future::sequential)
set.seed(3); n <- 500
d <- data.frame(id = seq_len(n), x1 = rnorm(n), x2 = rnorm(n), treat = rbinom(n, 1, 0.5))
t0 <- rexp(n, 0.1 * exp(ifelse(d$x1 > 0.37, log(2.5), 0) * d$treat)); cc <- runif(n, 5, 30)
d$tte <- pmin(t0, cc); d$event <- as.integer(t0 <= cc)
prior <- list(tree.cuts = "x1 <= 0.37")   # shape of an earlier grf.subg.harm.survival() result
run <- function(...) suppressWarnings(suppressMessages(forestsearch(d, confounders.name = c("x1", "x2"),
  use_grf = TRUE, use_lasso = FALSE, collapse_cuts = FALSE, maxk = 1, quiet = TRUE,
  parallel_args = list(plan = "sequential", workers = 1, show_message = FALSE), ...)))
a <- run(grf_res = prior)                             # documented reuse
b <- run(grf_res = prior, grf_cuts = "x1 <= 0.37")   # both supplied
cc2 <- run(grf_cuts = "x2 <= 9.99")                  # cuts only: refit overwrites them
cat("grf_res only : class", class(a)[1], "| grf_cuts =", deparse(a$grf_cuts), "| n factors:", length(a$confounders.evaluated), "\n")
cat("grf_res+cuts : class", class(b)[1], "| grf_cuts =", deparse(b$grf_cuts), "| n factors:", length(b$confounders.evaluated), "\n")
cat("grf_cuts only: grf_cuts =", deparse(cc2$grf_cuts), "\n")
ok <- inherits(a, "forestsearch") && inherits(b, "forestsearch") && is.null(a$grf_cuts) && identical(b$grf_cuts, "x1 <= 0.37")
cat("REPRO_VERDICT:", if (ok) "confirmed" else "not reproduced", "| supplied grf_cuts overwritten:", !identical(cc2$grf_cuts, "x2 <= 9.99"), "\n")
```

```text
grf_res only : class forestsearch | grf_cuts = NULL | n factors: 6 
grf_res+cuts : class forestsearch | grf_cuts = "x1 <= 0.37" | n factors: 7 
grf_cuts only: grf_cuts = c("x1 <= 0.2", "x2 <= -0.27", "x1 <= 0.27") 
REPRO_VERDICT: confirmed | supplied grf_cuts overwritten: TRUE
```

### V072

- Finding: A-020 — The use_dina selected-cut screening call to dina_subgroup() omits tau_sign, so for continuous outcomes (default adverse_outcome = FALSE) and binary outcomes wit
- Runs: exit 0, 2.5 s. Verdict: confirmed

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(7); n <- 400
d <- data.frame(id = seq_len(n), x1 = rnorm(n), x2 = rnorm(n), treat = rbinom(n, 1, 0.5))
d$y <- 10 + d$treat * ifelse(d$x1 > 0, 3, -3) + rnorm(n)   # treatment LOWERS y on x1 <= 0: the harm side when adverse_outcome = FALSE
fit <- dina(d, outcome = "y", treatment = "treat", covariates = c("x1", "x2"), family = "gaussian")
a <- list(fit = fit, df = d, covariates = c("x1", "x2"), m_diff = 1.25, n_min = 60, max_depth = 1L, sg_focus = "hr")
scr <- do.call(dina_subgroup, a)                                                          # argument set of forestsearch_main.R:3239 (no tau_sign)
sel <- do.call(dina_subgroup, c(a, list(tau_sign = .dina_tau_sign("continuous", FALSE))))  # argument set of forestsearch_helpers.R:1590
no_ts <- !any(grepl("tau_sign", deparse(body(forestsearch)), fixed = TRUE))
cat("forestsearch() body never mentions tau_sign:", no_ts, "\n")
cat("screening-style call :", scr$covariate, scr$direction, signif(scr$threshold, 3), "\n")
cat("selection-style call :", sel$covariate, sel$direction, signif(sel$threshold, 3), "\n")
cat("REPRO_VERDICT:", if (no_ts && isTRUE(scr$found) && isTRUE(sel$found) && !identical(scr$direction, sel$direction)) "confirmed" else "not reproduced", "\n")
```

```text
forestsearch() body never mentions tau_sign: TRUE 
screening-style call : x1 right 1.05 
selection-style call : x1 left -0.983 
REPRO_VERDICT: confirmed
```

### V073

- Finding: A-033 — df.analysis is never coerced to a base data.frame, and Y, Event, Treat, id and the dummy matrix are built with df[, name], which does not drop to a vector for a
- Runs: exit 0, 3.4 s. Verdict: confirmed

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(101); n <- 500
d <- data.frame(id = seq_len(n), x1 = rbinom(n, 1, 0.5), x2 = rbinom(n, 1, 0.5), x3 = rnorm(n), treat = rbinom(n, 1, 0.5))
t0 <- rexp(n, 0.10 * exp(ifelse(d$x1 == 1, 1.2, -0.6) * d$treat)); cc <- rexp(n, 0.03)
d$tte <- pmin(t0, cc); d$event <- as.numeric(t0 <= cc)
run <- function(x) { w <- character()
  r <- tryCatch(withCallingHandlers(forestsearch(df.analysis = x, confounders.name = c("x1", "x2", "x3"), n.min = 40, d0.min = 5, d1.min = 5,
      maxk = 1, details = FALSE, quiet = TRUE, parallel_args = list(plan = "sequential", workers = 1L, show_message = FALSE)),
    warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") },
    message = function(c) invokeRestart("muffleMessage")), error = function(e) structure(list(msg = conditionMessage(e)), class = "fs_err"))
  s <- if (inherits(r, "fs_err")) paste("ERROR:", r$msg) else if (!inherits(r, "forestsearch")) paste("early return at stage:", r$error_log$stage)
       else paste("ok; sg.harm:", paste(r$sg.harm, collapse = " & "))
  cat(s, "\n"); if (length(w)) cat(paste0("  warning: ", unique(w)), sep = "\n"); s }
cat("data.frame -> "); s_df <- run(d)
cat("data.table -> "); s_dt <- run(data.table::as.data.table(d))
s_tb <- if (requireNamespace("tibble", quietly = TRUE)) { cat("tibble     -> "); run(tibble::as_tibble(d)) } else s_df
cat("REPRO_VERDICT:", if (startsWith(s_df, "ok") && (!identical(s_dt, s_df) || !identical(s_tb, s_df))) "confirmed" else "not reproduced", "\n")
```

```text
data.frame -> ok; sg.harm: {x1} 
data.table -> ERROR: j (the 2nd argument inside [...]) is a single symbol but column name 'outcome.name' is not found. If you intended to select columns using a variable in calling scope, please try DT[, ..outcome.name]. The .. prefix conveys one-level-up similar to a file system path. 
tibble     -> ERROR: cannot xtfrm data frames 
REPRO_VERDICT: confirmed
```

### V074

- Finding: A-034 — With propensity adjustment requested on a GLM outcome (ps_method not "none" and ps_adjust_method "iptw" or "dr_gcomp") the weights are used only in the candidat
- Runs: exit 0, 3.0 s. Verdict: confirmed

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(11); n <- 600
d <- data.frame(id = seq_len(n), x1 = rbinom(n, 1, 0.5), x2 = rbinom(n, 1, 0.5), x3 = rnorm(n), treat = rbinom(n, 1, 0.5))
d$resp <- rbinom(n, 1, plogis(-1 + d$treat * ifelse(d$x1 == 1, 1.5, -0.3)))
seen <- NULL; ns <- asNamespace("forestsearch")
suppressMessages(trace("subgroup.consistency", where = ns, print = FALSE,
  tracer = quote(assign("seen", list(cols = names(df), spec = names(glm_resample_spec)), envir = globalenv()))))
fs <- suppressWarnings(suppressMessages(forestsearch(df.analysis = d, outcome.name = "resp", event.name = "resp", treat.name = "treat", id.name = "id",
  confounders.name = c("x1", "x2", "x3"), outcome_type = "binary", ps_method = "logistic", ps_adjust_method = "iptw",
  n.min = 60, maxk = 1, details = FALSE, quiet = TRUE, parallel_args = list(plan = "sequential", workers = 1L, show_message = FALSE))))
suppressMessages(untrace("subgroup.consistency", where = ns))
cat("reported ps_adjust_method:", fs$args_call_all$ps_adjust_method, "\n")
cat("consistency frame columns:", paste(seen$cols, collapse = ", "), "\n")
cat("resample spec fields     :", paste(seen$spec, collapse = ", "), "\n")
cat("REPRO_VERDICT:", if (is.null(seen)) "inconclusive (consistency stage not reached)" else if (identical(fs$args_call_all$ps_adjust_method, "iptw") && !any(c("sw", "ips_covar", "ps_hat") %in% seen$cols) && !any(grepl("ps|sw", seen$spec))) "confirmed" else "not reproduced", "\n")
```

```text
reported ps_adjust_method: iptw 
consistency frame columns: Y, Event, Treat, id, q1.0, q1.1, q2.0, q2.1, q3.0, q3.1, q4.0, q4.1, q5.0, q5.1, resp, treat 
resample spec fields     : outcome_type, effect_measure, treat.name, outcome.name, offset.name, adjust_covariates, adverse_outcome, comparison_threshold 
REPRO_VERDICT: confirmed
```

### V075

- Finding: A-035 — On GLM outcomes the user's outcome/treatment/offset column is attached to the consistency frame only if its name is not already taken, so a GLM outcome column l
- Runs: exit 0, 4.0 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(11); n <- 600
d <- data.frame(id = seq_len(n), x1 = rbinom(n, 1, 0.5), x2 = rbinom(n, 1, 0.5), treat = rbinom(n, 1, 0.5))
d$resp <- rbinom(n, 1, plogis(-1 + d$treat * ifelse(d$x1 == 1, 1.5, -0.3)))
d$Event <- d$resp   # identical outcome under a different column name
run <- function(nm) suppressWarnings(suppressMessages(forestsearch(df.analysis = d[, c("id", "x1", "x2", "treat", nm)],
  outcome.name = nm, treat.name = "treat", id.name = "id", confounders.name = c("x1", "x2"), outcome_type = "binary",
  n.min = 60, d0.min = 5, d1.min = 5, maxk = 1, details = FALSE, quiet = TRUE, plot.sg = FALSE,
  parallel_args = list(plan = "sequential", workers = 1L, show_message = FALSE))))
a <- run("resp"); b <- run("Event")
lab <- function(f) if (is.null(f$sg.harm)) "NULL" else paste(f$sg.harm, collapse = " & ")
cat("outcome column 'resp'  -> sg.harm:", lab(a), "\noutcome column 'Event' -> sg.harm:", lab(b), "\n")
cat("REPRO_VERDICT:", if (identical(lab(a), "NULL")) "INCONCLUSIVE (control run found nothing)" else if (!identical(lab(a), lab(b))) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
outcome column 'resp'  -> sg.harm: {x1} 
outcome column 'Event' -> sg.harm: NULL 
REPRO_VERDICT: CONFIRMED
```

### V076

- Finding: A-037 — forestsearch() never passes its seedit to subgroup.consistency(seed = ), so the random consistency splits are always seeded with the callee's own default 831695
- Runs: exit 0, 4.1 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(101); n <- 500
d <- data.frame(id = seq_len(n), x1 = rbinom(n, 1, 0.5), x2 = rbinom(n, 1, 0.5), treat = rbinom(n, 1, 0.5))
t0 <- rexp(n, 0.10 * exp(ifelse(d$x1 == 1, 1.2, -0.6) * d$treat)); cc <- rexp(n, 0.03)
d$tte <- pmin(t0, cc); d$event <- as.numeric(t0 <= cc)
run <- function(s) suppressWarnings(suppressMessages(forestsearch(df.analysis = d, outcome.name = "tte", event.name = "event",
  treat.name = "treat", id.name = "id", confounders.name = c("x1", "x2"), seedit = s,
  consistency_method = "split", use_twostage = FALSE, fs.splits = 40, n.min = 40, d0.min = 5, d1.min = 5, maxk = 1,
  details = FALSE, quiet = TRUE, plot.sg = FALSE, parallel_args = list(plan = "sequential", workers = 1L, show_message = FALSE))))
a <- run(123L); b <- run(456L)
sa <- a$grp.consistency$seed; sb <- b$grp.consistency$seed
cat("'seed' is a forestsearch formal:", "seed" %in% names(formals(forestsearch)), "| subgroup.consistency default seed:", formals(subgroup.consistency)$seed, "\n")
cat("seedit = 123 -> consistency seed:", format(sa), "| seedit = 456 -> consistency seed:", format(sb), "\n")
cat("REPRO_VERDICT:", if (is.null(sa) || is.null(sb)) "INCONCLUSIVE (consistency stage did not run)" else if (identical(sa, sb) && sa == 8316951) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
'seed' is a forestsearch formal: FALSE | subgroup.consistency default seed: 8316951 
seedit = 123 -> consistency seed: 8316951 | seedit = 456 -> consistency seed: 8316951 
REPRO_VERDICT: CONFIRMED
```

### V077

- Finding: A-039 — On the survival MR path adjust_covariates is intersected with the column names of the scoring frame, so any adjustment term that is not a bare column name (the 
- Runs: exit 0, 2.6 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
strata <- survival::strata
set.seed(7); n <- 800
d <- data.frame(id = seq_len(n), treat = rbinom(n, 1, 0.5), x1 = rbinom(n, 1, 0.5), x2 = rbinom(n, 1, 0.5), s = rbinom(n, 1, 0.5))
t0 <- rexp(n, 0.05 * exp(1.6 * d$s + ifelse(d$x1 == 1, 0.9, -0.3) * d$treat)); cc <- runif(n, 5, 40)
d$tte <- pmin(t0, cc); d$event <- as.integer(t0 <= cc)
fit <- try(suppressWarnings(suppressMessages(forestsearch(d, confounders.name = c("x1", "x2"), outcome.name = "tte", event.name = "event",
  treat.name = "treat", id.name = "id", maxk = 1, adjust_covariates = "strata(s)", mr_inference = TRUE,
  mr_inference_args = list(draws = 200, ci_method = "ij"), details = FALSE, quiet = TRUE, plot.sg = FALSE,
  parallel_args = list(plan = "sequential", workers = 1L, show_message = FALSE)))), silent = TRUE)
mr <- if (inherits(fit, "try-error")) NULL else fit$mr_inference
if (is.null(mr$naive$est)) { cat("REPRO_VERDICT: INCONCLUSIVE (no MR result)\n") } else {
  sel <- d[fit$grp.consistency$sg.harm.id == 1, ]
  hr_u <- exp(coef(survival::coxph(survival::Surv(tte, event) ~ treat, data = sel)))[["treat"]]
  hr_a <- exp(coef(survival::coxph(survival::Surv(tte, event) ~ treat + strata(s), data = sel)))[["treat"]]
  cat("MR naive HR:", mr$naive$est, "| unadjusted Cox:", hr_u, "| strata(s)-adjusted Cox:", hr_a, "\n")
  mu <- abs(mr$naive$est - hr_u) < 1e-4; ma <- abs(mr$naive$est - hr_a) < 1e-4
  cat("REPRO_VERDICT:", if (mu && !ma && abs(hr_u - hr_a) > 1e-3) "CONFIRMED" else if (ma && !mu) "NOT_REPRODUCED" else "INCONCLUSIVE", "\n")
}
```

```text
MR naive HR: 1.693294 | unadjusted Cox: 1.693294 | strata(s)-adjusted Cox: 2.247317 
REPRO_VERDICT: CONFIRMED
```

### V078

- Finding: A-040 — On the consistency branch (the default subgroup_method) forestsearch() calls fs_mr_inference() without field_R_out / field_R_in, so those two mr_inference_args 
- Runs: exit 0, 11.0 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
static <- !grepl("field_R_out", paste(deparse(body(forestsearch)), collapse = "\n"), fixed = TRUE)
set.seed(1); n <- 500
d <- data.frame(id = seq_len(n), treat = rbinom(n, 1, 0.5), x1 = rnorm(n), x2 = rnorm(n))
d$time <- rexp(n, 0.1 * ifelse(d$treat == 1, ifelse(d$x1 < -0.6, 3, 0.45), 1)); d$event <- as.integer(d$time < 30); d$time <- pmin(d$time, 30)
fit <- try(suppressWarnings(suppressMessages(forestsearch(df.analysis = d, outcome.name = "time", event.name = "event", treat.name = "treat",
  id.name = "id", confounders.name = c("x1", "x2"), maxk = 1, sg_focus = "hrMaxSG", mr_inference = TRUE,
  mr_inference_args = list(draws = 200, field_R_out = 50, field_R_in = 25), quiet = TRUE, details = FALSE, plot.sg = FALSE,
  parallel_args = list(plan = "sequential", workers = 1L, show_message = FALSE)))), silent = TRUE)
f <- if (inherits(fit, "try-error")) NULL else fit$mr_inference$field
cat("forestsearch() body never mentions field_R_out:", static, "\n")
if (!is.null(f$R_out)) cat("requested R_out / R_in: 50 / 25; used:", f$R_out, "/", f$R_in, "\n")
cat("REPRO_VERDICT:", if (!is.null(f$R_out)) (if (f$R_out != 50L) "CONFIRMED" else "NOT_REPRODUCED") else if (static) "CONFIRMED (static only: no field block was produced)" else "NOT_REPRODUCED", "\n")
```

```text
forestsearch() body never mentions field_R_out: TRUE 
requested R_out / R_in: 50 / 25; used: 1000 / 500 
REPRO_VERDICT: CONFIRMED
```

### V079

- Finding: INV-2 — fs_mr_inference() (and fs_mr_field_uniform()) call set.seed() on the global RNG and never restore .Random.seed, so a direct call with a seed, forestsearch(mr_in
- Runs: exit 0, 2.6 s. Verdict: CONFIRMED (caller's stream replaced by the MR seed)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 300
d <- data.frame(Y = rexp(n), Event = rbinom(n, 1, 0.7), Treat = rbinom(n, 1, 0.5), x = rnorm(n), z = rnorm(n))
cands <- list(a = which(d$x > 0), b = which(d$z > 0))
spec <- list(outcome_type = "survival", effect_measure = "HR", treat.name = "Treat", outcome.name = "Y",
             event.name = "Event", offset.name = NULL, adjust_covariates = NULL, adverse_outcome = TRUE)
run <- function() try(fs_mr_inference(d, cands, spec, cands$a, admission = list(effect_floor = NULL, consistency = NULL),
                                      reselection = "maxeff", draws = 100L, ci_method = "ij", seed = 7L), silent = TRUE)
set.seed(42); before <- .Random.seed; r <- run(); restored <- identical(before, .Random.seed); a <- runif(1)
set.seed(43); invisible(run()); b <- runif(1)
cat("call completed without error:", !inherits(r, "try-error"), "\n")
cat("RNG state restored after fs_mr_inference(seed = 7):", restored, "\n")
cat("next draw identical after caller seeds 42 and 43:", identical(a, b), "\n")
cat("REPRO_VERDICT:", if (!restored && identical(a, b)) "CONFIRMED (caller's stream replaced by the MR seed)" else "NOT_REPRODUCED", "\n")
```

```text
call completed without error: TRUE 
RNG state restored after fs_mr_inference(seed = 7): FALSE 
next draw identical after caller seeds 42 and 43: TRUE 
REPRO_VERDICT: CONFIRMED (caller's stream replaced by the MR seed)
```

### V080

- Finding: (not in the table) — On the effect-floor-only admission branch fs_mr_inference() passes a length-1 t_g to fs_mr_field_uniform(), which indexes it as a length-K vector (t_g[keep]), s
- Runs: exit 0, 2.5 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
# Part 1: the sweep, same effect floor as a length-K vector (documented contract) and as a scalar.
K <- 5; Sigma <- diag(K) * 0.04
w <- c(0.10, 0.50, 0.30, 0.20, 0.00); p_hat <- c(0.05, 0.60, 0.20, 0.10, 0.05)
run <- function(tg) suppressWarnings(fs_mr_field_uniform(
  Sigma = Sigma, w = w, sel = 2L, sigma_sel = 0.2, p_hat = p_hat, t_g = tg,
  reselection = "maxeff", delta_grid = c(0, 2), R_rep = 40L, R_out = 40L,
  R_in = 20L, kappa_grid = seq(1, 3, by = 0.1), seed = 11L))
a <- run(rep(5, K)); b <- run(5)            # a floor of 5 that no candidate can clear
cat("vector floor n_kept:", a$n_kept, "| scalar floor n_kept:", b$n_kept, "\n")
part1 <- all(a$n_kept == 0L) && any(b$n_kept > 0L)
# Part 2: what fs_mr_inference() actually passes on the effect-floor-only branch (sweep stubbed in this session only).
cap <- new.env()
ok2 <- tryCatch({
  set.seed(1); n <- 300
  d <- data.frame(Y = rexp(n), Event = rbinom(n, 1, 0.7), Treat = rbinom(n, 1, 0.5), x = rnorm(n), z = rnorm(n))
  cands <- list(a = which(d$x > 0), b = which(d$z > 0), c = which(d$x > 0 & d$z > 0))
  spec <- list(outcome_type = "survival", effect_measure = "HR", treat.name = "Treat", outcome.name = "Y",
               event.name = "Event", offset.name = NULL, adjust_covariates = NULL, adverse_outcome = TRUE)
  utils::assignInNamespace("fs_mr_field_uniform", function(..., w, t_g = NULL) {
    cap$K <- length(w); cap$len_tg <- length(t_g); stop("captured") }, ns = "forestsearch")
  invisible(fs_mr_inference(d, cands, spec, cands$a, admission = list(effect_floor = log(0.5), consistency = NULL),
    reselection = "maxeff", draws = 200L, ci_method = "field", field_uniform = TRUE, field_R_out = 30L,
    field_R_in = 10L, field_complement = FALSE, include_complement = FALSE, seed = 1L))
  TRUE }, error = function(e) { cat("part 2 did not run:", conditionMessage(e), "\n"); FALSE })
cat("caller passed length(t_g) =", cap$len_tg, "for K =", cap$K, "\n")
part2 <- isTRUE(ok2) && identical(cap$len_tg, 1L) && isTRUE(cap$K > 1L)
cat("REPRO_VERDICT:", if (part1 && part2) "CONFIRMED" else if (part1) "PARTIAL (sweep mis-indexes a scalar floor; caller capture did not run)" else "NOT_REPRODUCED", "\n")
```

```text
vector floor n_kept: 0 0 | scalar floor n_kept: 40 40 
caller passed length(t_g) = 1 for K = 3 
REPRO_VERDICT: CONFIRMED
```

### V081

- Finding: A-041 — .fs_apply_mr() forwards 25 of fs_mr_inference()'s 29 arguments and drops keep_declaration_field, keep_field_matrix, declaration_c0 and pconsistency.digits, so o
- Runs: exit 0, 2.6 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 300
d <- data.frame(Y = rexp(n), Event = rbinom(n, 1, 0.7), Treat = rbinom(n, 1, 0.5), x = rnorm(n), z = rnorm(n))
cands <- list(a = which(d$x > 0), b = which(d$z > 0), c = which(d$x > 0 & d$z > 0))
spec <- list(outcome_type = "survival", effect_measure = "HR", treat.name = "Treat", outcome.name = "Y",
             event.name = "Event", offset.name = NULL, adjust_covariates = NULL, adverse_outcome = TRUE)
adm <- list(effect_floor = log(0.5), consistency = NULL)   # the DINA/GRF admission shape
r_wrap <- .fs_apply_mr(d, cands, cands$a, spec, admission = adm, effect_neighborhood = 0.1,
  reselection_default = "maxeff",
  mr_inference_args = list(draws = 200L, ci_method = "ij", keep_declaration_field = TRUE), seedit = 1L)
r_direct <- fs_mr_inference(d, cands, spec, cands$a, admission = adm, reselection = "maxeff", draws = 200L,
  ci_method = "ij", include_complement = TRUE, seed = 1L, keep_declaration_field = TRUE)
src <- paste(deparse(body(.fs_apply_mr)), collapse = " ")
fm <- names(formals(fs_mr_inference))
not_fwd <- fm[!vapply(fm, function(a) grepl(paste0(a, " *="), src), logical(1))]
cat("formals not named in the wrapper's call:", not_fwd, "\n")
cat("wrapper keeps declaration_field:", !is.null(r_wrap$declaration_field), "| direct call keeps it:", !is.null(r_direct$declaration_field), "\n")
cat("REPRO_VERDICT:", if (!is.null(r_wrap) && is.null(r_wrap$declaration_field) && !is.null(r_direct$declaration_field)) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
formals not named in the wrapper's call: keep_declaration_field keep_field_matrix declaration_c0 pconsistency.digits 
wrapper keeps declaration_field: FALSE | direct call keeps it: TRUE 
REPRO_VERDICT: CONFIRMED
```

### V082

- Finding: A-042 — For binary, continuous and count outcomes lasso_selection() returns "(Intercept)" in `selected`, so get_FSdata()'s 'if any selected' guard is always true and, w
- Runs: exit 0, 2.4 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
cn <- c("z1", "z2", "z3")
mk <- function(s) { set.seed(s); n <- 300
  data.frame(y = 5 + rnorm(n), z1 = rnorm(n), z2 = rnorm(n), z3 = rnorm(n), ev = 1) }
hit <- NULL; has_int <- NA
for (s in 1:40) {                      # look for a draw where LASSO keeps no covariate
  df <- mk(s)
  ls <- lasso_selection(df, confounders.name = cn, outcome.name = "y", event.name = "ev", outcome_type = "continuous")
  has_int <- "(Intercept)" %in% ls$selected
  if (!length(intersect(ls$selected, cn))) { hit <- s; break }
}
cat("'(Intercept)' in selected:", has_int, "| selected:", ls$selected, "\n")
res <- "no all-zero LASSO draw found in 40 seeds"
if (!is.null(hit)) res <- tryCatch({
  fd <- get_FSdata(df, use_lasso = TRUE, confounders.name = cn, outcome.name = "y", event.name = "ev",
                   details = FALSE, outcome_type = "continuous")
  paste("returned", length(fd$confs), "cuts")
}, error = function(e) paste("ERROR:", substr(conditionMessage(e), 1, 70)))
cat("continuous outcome, LASSO kept no covariate (data seed", hit, "):", res, "\n")
cat("REPRO_VERDICT:", if (isTRUE(has_int) && grepl("^ERROR: get_FSdata: no candidate", res)) "CONFIRMED" else
  if (isTRUE(has_int)) "PARTIAL (intercept is in the selected set; empty-selection case not exercised)" else "NOT_REPRODUCED", "\n")
```

```text
'(Intercept)' in selected: TRUE | selected: (Intercept) 
continuous outcome, LASSO kept no covariate (data seed 2 ): ERROR: get_FSdata: no candidate cut expressions could be constructed from the 
REPRO_VERDICT: CONFIRMED
```

### V083

- Finding: A-044 — Missing confounder values are counted as a distinct level both when a variable is classed continuous and when a cut is tested for having more than one level, so
- Runs: exit 0, 2.3 s. Verdict: REPRODUCED (one-level candidate factor returned)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 300
x <- rep(c(1, 2, 3), 100); x[1:10] <- NA
df <- data.frame(tte = rexp(n), event = 1, x = x, z = rnorm(n))
o <- get_FSdata(df, confounders.name = c("x", "z"), outcome.name = "tte",
                event.name = "event", details = FALSE)
nl <- vapply(o$df[o$confs_names], nlevels, integer(1)); names(nl) <- o$confs
print(nl)
cat("3-level x given continuous cuts:", any(grepl("^x <=", o$confs)), "\n")
cat("REPRO_VERDICT:", if (any(nl == 1L)) "REPRODUCED (one-level candidate factor returned)" else "NOT_REPRODUCED", "\n")
```

```text
   x <= 2    x <= 1    x <= 3    z <= 0 z <= -0.7  z <= 0.7 
        2         2         1         2         2         2 
3-level x given continuous cuts: TRUE 
REPRO_VERDICT: REPRODUCED (one-level candidate factor returned)
```

### V084

- Finding: A-033 — A tibble or data.table passes the entry check uncoerced (the stop message that tells such users to coerce can never fire for them), while downstream code extrac
- Runs: exit 0, 2.7 s. Verdict: REPRODUCED (tibble/data.table input does not behave like a data.frame)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 200
d <- data.frame(id = seq_len(n), tte = rexp(n), event = rbinom(n, 1, 0.7),
                treat = rbinom(n, 1, 0.5), x1 = rnorm(n), x2 = rbinom(n, 1, 0.5))
run <- function(dd) tryCatch({
  o <- get_FSdata(dd, confounders.name = c("x1", "x2"), outcome.name = "tte", event.name = "event", details = FALSE)
  nm <- "tte"; Y <- o$df[, nm]                                  # as forestsearch_main.R:3350
  Z <- dummy_encode(o$df[, o$confs_names, drop = FALSE])       # as forestsearch_main.R:3386-3387
  list(Y_is_vector = is.atomic(Y), dimZ = dim(Z), namesZ = names(Z))
}, error = function(e) paste("ERROR:", conditionMessage(e)))
base <- run(d); dt <- run(data.table::as.data.table(d))
tb <- if (requireNamespace("tibble", quietly = TRUE)) run(tibble::as_tibble(d)) else base
str(list(base = base, data.table = dt, tibble = tb))
cat("REPRO_VERDICT:", if (!identical(dt, base) || !identical(tb, base)) "REPRODUCED (tibble/data.table input does not behave like a data.frame)" else "NOT_REPRODUCED", "\n")
```

```text
List of 3
 $ base      :List of 3
  ..$ Y_is_vector: logi TRUE
  ..$ dimZ       : int [1:2] 200 8
  ..$ namesZ     : chr [1:8] "q1.0" "q1.1" "q2.0" "q2.1" ...
 $ data.table: chr "ERROR: j (the 2nd argument inside [...]) is a single symbol but column name 'nm' is not found. If you intended "| __truncated__
 $ tibble    : chr "ERROR: cannot xtfrm data frames"
REPRO_VERDICT: REPRODUCED (tibble/data.table input does not behave like a data.frame)
```

### V085

- Finding: A-045 — Supplying conf.cont_medians or conf.cont_medians_force silently removes every cut on all other continuous confounders, and under the default cut_type the variab
- Runs: exit 0, 2.7 s. Verdict: REPRODUCED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
df <- survival::gbsg; df$grade3 <- as.integer(df$grade == 3)
cf <- c("age", "meno", "size", "grade3", "nodes", "pgr", "er")
f <- function(...) get_FSdata(df, confounders.name = cf, outcome.name = "rfstime",
                              event.name = "status", details = FALSE, ...)$confs
vars <- function(s) unique(sub(" .*$", "", s))
a <- f(); b <- f(conf.cont_medians_force = "age"); m <- f(conf.cont_medians = "age")
print(a); print(b); print(m)
lost_b <- setdiff(vars(a), vars(b)); lost_m <- setdiff(vars(a), vars(m))
cat("lost with _force:", lost_b, "| lost with conf.cont_medians:", lost_m, "| age cuts under conf.cont_medians:", sum(grepl("^age ", m)), "\n")
cat("REPRO_VERDICT:", if (length(lost_b) > 0 && length(lost_m) > 0 && sum(grepl("^age ", m)) > 1) "REPRODUCED" else "NOT_REPRODUCED", "\n")
```

```text
 [1] "age <= 53"  "age <= 46"  "age <= 61"  "size <= 29" "size <= 25"
 [6] "size <= 20" "size <= 35" "nodes <= 5" "nodes <= 3" "nodes <= 1"
[11] "nodes <= 7" "pgr <= 110" "pgr <= 33"  "pgr <= 7"   "pgr <= 132"
[16] "er <= 96"   "er <= 36"   "er <= 8"    "er <= 114"  "meno"      
[21] "grade3"    
[1] "age <= 53" "meno"      "grade3"   
[1] "age <= 53" "age <= 46" "age <= 61" "meno"      "grade3"   
lost with _force: size nodes pgr er | lost with conf.cont_medians: size nodes pgr er | age cuts under conf.cont_medians: 3 
REPRO_VERDICT: REPRODUCED
```

### V086

- Finding: A-047 — Existing candidates are compared with GRF and DINA cuts by unanchored grepl() regex matching, so a GRF or DINA cut is discarded as already represented whenever 
- Runs: exit 0, 2.6 s. Verdict: REPRODUCED (GRF cut on bmi dropped only because 'bm' is a confounder name)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 200
df <- data.frame(tte = rexp(n), event = rbinom(n, 1, 0.7), bm = rbinom(n, 1, 0.5),
                 bmi = rnorm(n, 27, 4), age = rnorm(n, 60, 10))
g <- c("bmi <= 25.3", "age <= 57.2")
f <- function(d, cf) get_FSdata(d, confounders.name = cf, use_grf = TRUE, grf_cuts = g,
                                outcome.name = "tte", event.name = "event", details = FALSE)$confs
bug <- f(df, c("bm", "bmi", "age"))
df2 <- df; names(df2)[names(df2) == "bm"] <- "zz"          # control: same data, no name overlap
ctl <- f(df2, c("zz", "bmi", "age"))
print(bug); print(ctl)
cat("REPRO_VERDICT:", if (!("bmi <= 25.3" %in% bug) && ("bmi <= 25.3" %in% ctl) && ("age <= 57.2" %in% bug)) "REPRODUCED (GRF cut on bmi dropped only because 'bm' is a confounder name)" else "NOT_REPRODUCED", "\n")
```

```text
[1] "age <= 57.2" "bmi <= 27"   "bmi <= 24"   "bmi <= 30"   "age <= 59"  
[6] "age <= 50"   "age <= 67"   "bm"         
[1] "bmi <= 25.3" "age <= 57.2" "bmi <= 27"   "bmi <= 24"   "bmi <= 30"  
[6] "age <= 59"   "age <= 50"   "age <= 67"   "zz"         
REPRO_VERDICT: REPRODUCED (GRF cut on bmi dropped only because 'bm' is a confounder name)
```

### V087

- Finding: D-001 — Caller-supplied dina_cuts enter the candidate pool whenever non-empty, with no use_dina gate, unlike grf_cuts (inert unless use_grf = TRUE), while the family st
- Runs: exit 1, 0.1 s. Verdict: no verdict line (exit 1)

```r
## INV-1-a: dina_cuts enter the pool with every front-end flag FALSE; grf_cuts do not.
## Run from tempdir() after pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1)
n <- 300
d <- data.frame(id = seq_len(n), tte = rexp(n), event = rbinom(n, 1, 0.7),
                treat = rbinom(n, 1, 0.5),
                age = round(rnorm(n, 55, 10), 1),
                bm = round(rnorm(n, 20, 5), 1),
                sexf = rbinom(n, 1, 0.5))
conf <- c("age", "bm", "sexf")
cut <- "age <= 51.37"
base <- get_FSdata(d, use_lasso = FALSE, use_grf = FALSE, confounders.name = conf, details = FALSE)
g    <- get_FSdata(d, use_lasso = FALSE, use_grf = FALSE, grf_cuts = cut, confounders.name = conf, details = FALSE)
dn   <- get_FSdata(d, use_lasso = FALSE, use_grf = FALSE, dina_cuts = cut, confounders.name = conf, details = FALSE)
cat("baseline candidates:\n"); print(base$confs)
cat("with grf_cuts, use_grf = FALSE:\n"); print(g$confs)
cat("with dina_cuts, no flag:\n"); print(dn$confs)
in_dina <- any(grepl("51.37", dn$confs, fixed = TRUE))
in_grf  <- any(grepl("51.37", g$confs, fixed = TRUE))
status <- .fs_family_status(list(use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE,
                                 subgroup_method = "consistency", dina_cuts = cut))
cat("dina cut in pool:", in_dina, "| grf cut in pool:", in_grf,
    "| grf pool identical to baseline:", identical(g$confs, base$confs),
    "| family_status:", status, "\n")
demonstrated <- isTRUE(in_dina) && !isTRUE(in_grf) && identical(g$confs, base$confs) &&
  identical(status, "no-front-end")
cat("REPRO_VERDICT:", if (demonstrated) "demonstrated" else "not demonstrated", "\n")
```

```text
Error in get_FSdata(d, use_lasso = FALSE, use_grf = FALSE, confounders.name = conf,  : 
  could not find function "get_FSdata"
Calls: source -> withVisible -> eval -> eval
Execution halted
```

### V088

- Finding: A-048 — conf_force expressions are not validated: one whose right-hand side is not a number (compound expression, unsupported function, typo) is compared as a character
- Runs: exit 0, 2.3 s. Verdict: REPRODUCED (compound conf_force accepted and mis-evaluated under its own label)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 300
df <- data.frame(tte = rexp(n), event = 1, age = round(runif(n, 30, 110)), er = rnorm(n))
lab <- "age <= 50 & er > 0"
o <- get_FSdata(df, confounders.name = c("age", "er"), conf_force = lab,
                outcome.name = "tte", event.name = "event", details = FALSE)
i <- match(lab, o$confs)
flag <- if (is.na(i)) rep(NA, n) else as.character(o$df[[o$confs_names[i]]]) == "1"
intended <- df$age <= 50 & df$er > 0
print(table(flag = flag, intended = intended, useNA = "ifany"))
print(table(flag = flag, age_100_plus = df$age >= 100, useNA = "ifany"))
cat("REPRO_VERDICT:", if (!is.na(i) && any(flag != intended)) "REPRODUCED (compound conf_force accepted and mis-evaluated under its own label)" else "NOT_REPRODUCED", "\n")
```

```text
       intended
flag    FALSE TRUE
  FALSE   192    0
  TRUE     79   29
       age_100_plus
flag    FALSE TRUE
  FALSE   192    0
  TRUE     69   39
REPRO_VERDICT: REPRODUCED (compound conf_force accepted and mis-evaluated under its own label)
```

### V089

- Finding: A-049 — exclude_cuts is applied only before collapse_redundant_cuts(), so the default coarsening (collapse_cuts = TRUE, digits = 0) can rewrite a surviving cut into exa
- Runs: exit 0, 2.5 s. Verdict: DEFECT_PRESENT (excluded cut returned after coarsening)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
x <- rep(c(30, 45, 49.6, 49.6, 55, 68.4), each = 50); n <- length(x)
set.seed(1)
df <- data.frame(tte = rexp(n), event = 1, x = x)
o  <- get_FSdata(df, confounders.name = "x", exclude_cuts = "x <= 50",
                 outcome.name = "tte", event.name = "event", details = FALSE)
o0 <- get_FSdata(df, confounders.name = "x", exclude_cuts = "x <= 50", collapse_cuts = FALSE,
                 outcome.name = "tte", event.name = "event", details = FALSE)
print(o$confs); print(o0$confs)
cat("REPRO_VERDICT:", if ("x <= 50" %in% o$confs) "DEFECT_PRESENT (excluded cut returned after coarsening)" else "NOT_REPRODUCED", "\n")
```

```text
[1] "x <= 50" "x <= 45" "x <= 55"
[1] "x <= 49.6" "x <= 45"   "x <= 55"  
REPRO_VERDICT: DEFECT_PRESENT (excluded cut returned after coarsening)
```

### V090

- Finding: A-050 — Cuts that evaluate to a single value are dropped silently (the only message is under details = TRUE, which forestsearch() hard-codes to FALSE) and the zero-cand
- Runs: exit 0, 2.2 s. Verdict: DEFECT_PRESENT (empty candidate set, no error, no warning)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 100
df <- data.frame(tte = rexp(n), event = 1, b = runif(n, 0.001, 0.04),
                 sex = factor(sample(c("F", "M"), n, TRUE)))
w <- character(0)
o <- withCallingHandlers(
  tryCatch(get_FSdata(df, confounders.name = c("b", "sex"), outcome.name = "tte",
                      event.name = "event", details = FALSE), error = function(e) e),
  warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
if (inherits(o, "error")) print(conditionMessage(o)) else { str(o$confs); str(o$confs_names) }
print(w)
cat("REPRO_VERDICT:", if (!inherits(o, "error") && length(o$confs) == 0L && length(w) == 0L)
  "DEFECT_PRESENT (empty candidate set, no error, no warning)" else "NOT_REPRODUCED", "\n")
```

```text
 chr(0) 
 NULL
character(0)
REPRO_VERDICT: DEFECT_PRESENT (empty candidate set, no error, no warning)
```

### V091

- Finding: A-051 — Candidate factor columns are written as q1, q2, ... with no collision check, so an input column carrying one of those names (for example a covariate called q1) 
- Runs: exit 0, 2.3 s. Verdict: DEFECT_PRESENT (numeric covariate q1 replaced by a 0/1 factor)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 200
df <- data.frame(tte = rexp(n), event = 1, q1 = rnorm(n, 50, 10), z = rnorm(n))
o <- get_FSdata(df, confounders.name = c("z", "q1"), outcome.name = "tte",
                event.name = "event", details = FALSE)
str(df$q1); str(o$df$q1); print(o$confs)
cat("REPRO_VERDICT:", if (is.numeric(df$q1) && is.factor(o$df$q1))
  "DEFECT_PRESENT (numeric covariate q1 replaced by a 0/1 factor)" else "NOT_REPRODUCED", "\n")
```

```text
 num [1:200] 62.1 37.7 59.8 52.2 35.3 ...
 Factor w/ 2 levels "0","1": 1 1 1 2 1 1 2 2 1 2 ...
[1] "z <= 0"    "z <= -0.7" "z <= 0.7"  "q1 <= 50"  "q1 <= 45"  "q1 <= 57" 
REPRO_VERDICT: DEFECT_PRESENT (numeric covariate q1 replaced by a 0/1 factor)
```

### V092

- Finding: A-052 — For a binary outcome with a rate measure (IRR) and adverse_outcome = FALSE, the search estimator fits the rate ratio of Y (the flag is not forwarded to the Pois
- Runs: exit 0, 2.4 s. Verdict: DEFECT_PRESENT (estimator and resample fit disagree only under adverse_outcome = FALSE)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 400
d <- data.frame(treat = rep(0:1, n / 2), y = rbinom(n, 1, 0.3), time = runif(n, 0.5, 2))
f <- function(adv) c(
  estimator = make_effect_estimator("binary", "treat", "y", offset.name = "time",
                                    effect_measure = "IRR", adverse_outcome = adv)(d)$estimate,
  resample  = .consistency_glm_pieces(d, "binary", "IRR", treat.name = "treat", outcome.name = "y",
                                      offset.name = "time", adverse_outcome = adv)$beta_hat)
a <- f(TRUE); b <- f(FALSE); print(rbind(adverse_TRUE = a, adverse_FALSE = b))
cat("REPRO_VERDICT:", if (isTRUE(all.equal(a[[1]], a[[2]])) && !isTRUE(all.equal(b[[1]], b[[2]])))
  "DEFECT_PRESENT (estimator and resample fit disagree only under adverse_outcome = FALSE)" else "NOT_REPRODUCED", "\n")
```

```text
                estimator    resample
adverse_TRUE  -0.04064802 -0.04064802
adverse_FALSE -0.04064802  0.03336849
REPRO_VERDICT: DEFECT_PRESENT (estimator and resample fit disagree only under adverse_outcome = FALSE)
```

### V093

- Finding: A-053 — effect_measure is validated against the outcome type for binary and count outcomes but not for continuous, so a ratio or risk measure supplied with outcome_type
- Runs: exit 0, 2.3 s. Verdict: DEFECT_PRESENT (continuous + OR accepted, returns MD; sibling branch rejects)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1)
d <- data.frame(treat = rep(0:1, 50), y = rnorm(100))
r <- tryCatch(make_effect_estimator("continuous", "treat", "y", effect_measure = "OR")(d), error = function(e) e)
b <- tryCatch(make_effect_estimator("binary", "treat", "y", effect_measure = "MD"), error = function(e) e)
if (!inherits(r, "error")) print(r[c("estimate", "measure", "method_used")])
cat("binary + MD rejected:", inherits(b, "error"), "\n")
cat("REPRO_VERDICT:", if (!inherits(r, "error") && identical(r$measure, "MD") && inherits(b, "error"))
  "DEFECT_PRESENT (continuous + OR accepted, returns MD; sibling branch rejects)" else "NOT_REPRODUCED", "\n")
```

```text
$estimate
[1] -0.1732765

$measure
[1] "MD"

$method_used
[1] "ols"

binary + MD rejected: TRUE 
REPRO_VERDICT: DEFECT_PRESENT (continuous + OR accepted, returns MD; sibling branch rejects)
```

### V094

- Finding: A-054 — A requested ps_adjust_method has no implementation for several measure/method pairs (RR with iptw or dr_gcomp, RD with dr_gcomp, MD with dr_gcomp), yet forestse
- Runs: exit 0, 2.3 s. Verdict: DEFECT_PRESENT (iptw changes OR but is silently ignored for RR)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(2); n <- 200
d <- data.frame(treat = rep(0:1, n / 2), y = rbinom(n, 1, .4), sw = runif(n, .5, 2))
mk <- function(m, ps) suppressWarnings(make_effect_estimator("binary", "treat", "y", effect_measure = m, ps_adjust_method = ps)(d))
rr_w <- mk("RR", "iptw"); rr_0 <- mk("RR", "none"); or_w <- mk("OR", "iptw"); or_0 <- mk("OR", "none")
print(c(RR_iptw = rr_w$estimate, RR_none = rr_0$estimate, OR_iptw = or_w$estimate, OR_none = or_0$estimate))
print(c(rr_w$method_used, or_w$method_used))
cat("REPRO_VERDICT:", if (identical(rr_w$estimate, rr_0$estimate) && !identical(or_w$estimate, or_0$estimate))
  "DEFECT_PRESENT (iptw changes OR but is silently ignored for RR)" else "NOT_REPRODUCED", "\n")
```

```text
    RR_iptw     RR_none     OR_iptw     OR_none 
-0.08004271 -0.08004271 -0.17629397 -0.12805193 
[1] "log_binomial"  "logistic_iptw"
REPRO_VERDICT: DEFECT_PRESENT (iptw changes OR but is silently ignored for RR)
```

### V095

- Finding: A-057 — The IPTW branch of the Poisson rate estimator never looks at effect_measure, so with effect_measure = "IRD" and ps_adjust_method = "iptw" it returns the log inc
- Runs: exit 0, 2.3 s. Verdict: CONFIRMED (IRD requested, log-IRR returned under iptw)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 400
d <- data.frame(treat = rep(0:1, n / 2), time = runif(n, 1, 3))
d$y <- rpois(n, 0.2 * d$time * ifelse(d$treat == 1, 1.5, 1)); d$sw <- 1
mk <- function(em, ps = "none") make_effect_estimator("count", "treat", "y", offset.name = "time", effect_measure = em, ps_adjust_method = ps)
u <- mk("IRD")(d); p <- mk("IRD", "iptw")(d); r <- mk("IRR")(d)
print(rbind(IRD_unadj = c(u$estimate, u$measure, u$method_used), IRD_iptw_unit_weights = c(p$estimate, p$measure, p$method_used), IRR_unadj = c(r$estimate, r$measure, r$method_used)))
ok <- identical(p$measure, "log-IRR") && isTRUE(all.equal(p$estimate, r$estimate)) && !isTRUE(all.equal(p$estimate, u$estimate))
cat("REPRO_VERDICT:", if (ok) "CONFIRMED (IRD requested, log-IRR returned under iptw)" else "NOT_REPRODUCED", "\n")
```

```text
                      [,1]                [,2]      [,3]                  
IRD_unadj             "0.213321486486412" "IRD"     "poisson_offset_delta"
IRD_iptw_unit_weights "0.446190326008219" "log-IRR" "poisson_iptw"        
IRR_unadj             "0.446190326008219" "log-IRR" "poisson_offset"      
REPRO_VERDICT: CONFIRMED (IRD requested, log-IRR returned under iptw)
```

### V096

- Finding: A-058 — .grf_code_column() derives the coding of a character column from the values present in the frame it is handed (unique(x), NA included), so a GRF cut can be eval
- Runs: exit 0, 2.3 s. Verdict: CONFIRMED (character coding depends on the frame)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
c_na  <- .grf_code_column(c("0", "1", NA, "1"))            # character with NA
c_cc  <- .grf_code_column(c("0", "1", "1"))                # same column after complete-case filtering (fit time)
c_fac <- .grf_code_column(factor(c("0", "1", NA, "1")))    # factor: stable
print(list(character_with_NA = c_na, complete_cases = c_cc, factor_with_NA = c_fac))
tr <- c("a", "b", "c", "b"); nw <- c("b", "c", "c")
print(list(train = .grf_code_column(tr), new = .grf_code_column(nw)))
def <- list(conjunctions = list(data.frame(variable = "g", op = "<=", value = 1, stringsAsFactors = FALSE)))
m_tr <- .grf_evaluate_subgroup(def, data.frame(g = tr, stringsAsFactors = FALSE))
m_nw <- .grf_evaluate_subgroup(def, data.frame(g = nw, stringsAsFactors = FALSE))
print(list(members_train = tr[m_tr == 0L], members_new = nw[m_nw == 0L]))   # 'a' vs 'b'
ok <- identical(c_cc, c(0, 1, 1)) && identical(c_na, c(1L, 2L, NA, 2L)) && identical(c_fac, c(0, 1, NA, 1)) && identical(tr[m_tr == 0L], "a") && identical(nw[m_nw == 0L], "b")
cat("REPRO_VERDICT:", if (ok) "CONFIRMED (character coding depends on the frame)" else "NOT_REPRODUCED", "\n")
```

```text
$character_with_NA
[1]  1  2 NA  2

$complete_cases
[1] 0 1 1

$factor_with_NA
[1]  0  1 NA  1

$train
[1] 1 2 3 2

$new
[1] 1 2 2

$members_train
[1] "a"

$members_new
[1] "b"

REPRO_VERDICT: CONFIRMED (character coding depends on the frame)
```

### V097

- Finding: A-059 — The survival GRF path builds its covariate matrix with apply(data[, confounders], 2, as.numeric), which passes through a character matrix whenever any confounde
- Runs: exit 0, 2.3 s. Verdict: CONFIRMED (7-digit rounding and NA coding in the survival matrix build)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
d <- data.frame(x = c(61.234561234, 2.5, 70.000000123), f = factor(c("0", "1", "0")))
X <- apply(d[, c("x", "f"), drop = FALSE], 2, as.numeric)   # grf_main.R:226 verbatim
print(cbind(raw = d$x, X = X[, "x"], diff = d$x - X[, "x"]), digits = 15)
Xn <- apply(d[, "x", drop = FALSE], 2, as.numeric)           # no factor present: no rounding
d2 <- data.frame(x = c(1.5, 2.5, 3.5), sex = factor(c("m", "f", "m")))
Xs <- suppressWarnings(apply(d2[, c("x", "sex"), drop = FALSE], 2, as.numeric))  # survival coding
Xg <- .build_grf_X(d2, c("x", "sex"))                                              # GLM coding
print(Xs); print(Xg)
ok <- any(d$x != X[, "x"]) && all(d$x == Xn[, "x"]) && all(is.na(Xs[, "sex"])) && !anyNA(Xg)
cat("REPRO_VERDICT:", if (ok) "CONFIRMED (7-digit rounding and NA coding in the survival matrix build)" else "NOT_REPRODUCED", "\n")
```

```text
              raw        X                 diff
[1,] 61.234561234 61.23456 1.23399999552021e-06
[2,]  2.500000000  2.50000 0.00000000000000e+00
[3,] 70.000000123 70.00000 1.23000006624352e-07
       x sex
[1,] 1.5  NA
[2,] 2.5  NA
[3,] 3.5  NA
       x sex
[1,] 1.5   2
[2,] 2.5   1
[3,] 3.5   2
REPRO_VERDICT: CONFIRMED (7-digit rounding and NA coding in the survival matrix build)
```

### V098

- Finding: A-060 — When validate_grf_data() rejects the data (one treatment level, one event level such as fully uncensored data, or an arm smaller than n.min), the early exit cal
- Runs: exit 0, 2.4 s. Verdict: REPRODUCED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 200
d <- data.frame(id = 1:n, tte = rexp(n), event = 1L, treat = rbinom(n, 1, 0.5),
                x1 = rnorm(n), x2 = rnorm(n))   # no censoring -> validate_grf_data() FALSE
r1 <- tryCatch(suppressWarnings(grf.subg.harm.survival(d, c("x1", "x2"), "tte", "event", "id", "treat", n.min = 20)),
               error = function(e) e)
d2 <- d; d2$event <- rbinom(n, 1, 0.7); d2$treat <- rep(c(1, 0), c(15, n - 15))  # arm < n.min
r2 <- tryCatch(suppressWarnings(grf.subg.harm.survival(d2, c("x1", "x2"), "tte", "event", "id", "treat", n.min = 20)),
               error = function(e) e)
print(r1); print(r2)
hit <- function(r) inherits(r, "error") && grepl("subscript out of bounds", conditionMessage(r))
cat("REPRO_VERDICT:", if (hit(r1) && hit(r2)) "REPRODUCED" else "NOT_REPRODUCED", "\n")
```

```text
<subscriptOutOfBoundsError in trees[[1]]: subscript out of bounds>
<subscriptOutOfBoundsError in trees[[1]]: subscript out of bounds>
REPRO_VERDICT: REPRODUCED
```

### V099

- Finding: A-061 — Both GRF identifiers test the frontier selection with is.null() only, but .grf_frontier_select() returns a zero-row data frame when the band is empty (effMaxSG/
- Runs: exit 0, 2.3 s. Verdict: REPRODUCED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
cand <- data.frame(v1 = "x", d1 = c("left", "right"), c1 = 0, v2 = NA_character_, d2 = NA_character_,
                   c2 = NA_real_, effect = c(-0.5, -0.2), size = c(100L, 80L), stringsAsFactors = FALSE)
sel <- .grf_frontier_select(cand, dmin = -1, rule = "effMaxSG", nbhd = 0.10)
cat("is.null(sel):", is.null(sel), " nrow(sel):", NROW(sel), "\n")
r <- tryCatch(.grf_sg_def_from_candidate(sel), error = function(e) e)
print(r)
cat("REPRO_VERDICT:", if (!is.null(sel) && NROW(sel) == 0L && inherits(r, "error")) "REPRODUCED" else "NOT_REPRODUCED", "\n")
```

```text
is.null(sel): FALSE  nrow(sel): 0 
<simpleError in if (cand_row$d1 == "left") "<=" else ">": argument is of length zero>
REPRO_VERDICT: REPRODUCED
```

### V100

- Finding: A-062 — When no DR candidate clears dmin.grf, the frontier no-selection return omits the candidates table, so the effect-based re-selection (default grf_select_statisti
- Runs: exit 0, 3.0 s. Verdict: REPRODUCED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(3); n <- 400
dd <- data.frame(id = 1:n, treat = rbinom(n, 1, 0.5), x1 = rnorm(n), x2 = rnorm(n))
t0 <- rexp(n, 0.05 * exp(0.7 * dd$treat * (dd$x1 > 0))); cens <- runif(n, 5, 40)
dd$tte <- pmin(t0, cens); dd$event <- as.integer(t0 <= cens)
run <- function(dmin) suppressWarnings(grf.subg.harm.survival(dd, c("x1", "x2"), "tte", "event", "id", "treat",
  n.min = 40, frac.tau = 0.8, grf_selection = "frontier", dmin.grf = dmin, frontier_rule = "eff"))
ok <- run(-1e6); nul <- run(1e6)
cat("winner found  -> candidates rows:", NROW(ok$candidates), "\n")
cat("no DR winner  -> has candidates:", !is.null(nul$candidates), "\n")
cat("REPRO_VERDICT:", if (NROW(ok$candidates) > 0L && is.null(nul$candidates)) "REPRODUCED" else "NOT_REPRODUCED", "\n")
```

```text
winner found  -> candidates rows: 88 
no DR winner  -> has candidates: FALSE 
REPRO_VERDICT: REPRODUCED
```

### V101

- Finding: A-063 — return_selected_cuts_only is accepted by grf.subg.harm.glm() but never acted on (cuts always come from the single best-depth tree), and frac.tau and conf.level 
- Runs: exit 0, 2.3 s. Verdict: REPRODUCED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
b <- deparse(body(grf.subg.harm.glm))
cnt <- function(a, bb = b) sum(grepl(a, bb, fixed = TRUE))
for (a in c("frac.tau", "conf.level", "return_selected_cuts_only"))
  cat(a, ":", cnt(a), "reference(s):", paste(trimws(grep(a, b, fixed = TRUE, value = TRUE)), collapse = " || "), "\n")
echo_only <- all(grepl("return_selected_cuts_only = return_selected_cuts_only", grep("return_selected_cuts_only", b, fixed = TRUE, value = TRUE), fixed = TRUE))
id_surv <- cnt("id.name", deparse(body(grf.subg.harm.survival)))
cat("id.name references in grf.subg.harm.survival body:", id_surv, "\n")
cat("REPRO_VERDICT:", if (cnt("frac.tau") == 0L && cnt("conf.level") == 0L && echo_only && id_surv == 0L) "REPRODUCED" else "NOT_REPRODUCED", "\n")
```

```text
frac.tau : 0 reference(s):  
conf.level : 0 reference(s):  
return_selected_cuts_only : 4 reference(s): return_selected_cuts_only = return_selected_cuts_only, || return_selected_cuts_only = return_selected_cuts_only, || return_selected_cuts_only = return_selected_cuts_only, || return_selected_cuts_only = return_selected_cuts_only, 
id.name references in grf.subg.harm.survival body: 0 
REPRO_VERDICT: REPRODUCED
```

### V102

- Finding: INV-2 — grf.subg.harm.glm() calls set.seed(seedit) on the global RNG and never restores the previous state, so every GLM-outcome GRF run, including inside forestsearch(
- Runs: exit 0, 2.5 s. Verdict: REPRODUCED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 300
d <- data.frame(id = 1:n, treat = rbinom(n, 1, 0.5), x1 = rnorm(n), x2 = rnorm(n))
d$y <- rbinom(n, 1, plogis(-0.5 + 0.8 * d$treat * (d$x1 > 0)))
set.seed(999); before <- .Random.seed
r <- tryCatch(suppressWarnings(grf.subg.harm.glm(d, c("x1", "x2"), "y", n.min = 30, seedit = 42)), error = function(e) e)
after <- .Random.seed
set.seed(42); same_as_seedit <- identical(after, .Random.seed)
cat("call errored:", inherits(r, "error"), "\n")
cat("RNG state unchanged by the call:", identical(before, after), "\n")
cat("state after the call equals set.seed(seedit):", same_as_seedit, "\n")
cat("REPRO_VERDICT:", if (!identical(before, after)) "REPRODUCED" else "NOT_REPRODUCED", "\n")
```

```text
call errored: FALSE 
RNG state unchanged by the call: FALSE 
state after the call equals set.seed(seedit): TRUE 
REPRO_VERDICT: REPRODUCED
```

### V103

- Finding: A-065 — With sg.criterion = 'Nsg' the GLM selector returns the largest subgroup without comparing its score with dmin.grf, whereas the 'mDiff' branch and the survival s
- Runs: exit 0, 2.3 s. Verdict: REPRODUCED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
vals <- list(list(depth = 1L, n_sg = 150, score = 0.01), list(depth = 2L, n_sg = 90, score = 0.50))
g <- .select_best_sg_glm(vals, "Nsg", dmin.grf = 0.2, n = 400, n.min = 60)
m <- .select_best_sg_glm(vals, "mDiff", dmin.grf = 0.2, n = 400, n.min = 60)
s <- select_best_subgroup(data.frame(diff = c(0.01, 0.50), Nsg = c(150, 90)), "Nsg", dmin.grf = 0.2, n.max = 400)
cat("GLM Nsg pick: depth", g$depth, "score", g$score, "\n")
cat("GLM mDiff pick: depth", m$depth, "score", m$score, "\n")
cat("survival Nsg pick: diff", s$diff, "Nsg", s$Nsg, "\n")
cat("REPRO_VERDICT:", if (!is.null(g) && g$score < 0.2 && isTRUE(s$diff >= 0.2)) "REPRODUCED" else "NOT_REPRODUCED", "\n")
```

```text
GLM Nsg pick: depth 1 score 0.01 
GLM mDiff pick: depth 2 score 0.5 
survival Nsg pick: diff 0.5 Nsg 90 
REPRO_VERDICT: REPRODUCED
```

### V104

- Finding: A-066 — consistency_threshold = NULL is documented as falling back to hr.consistency, but on the split path it is never resolved: with estimator_fn supplied and the thr
- Runs: exit 0, 2.6 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(7); n <- 300
q1 <- rbinom(n, 1, 0.5); q2 <- rbinom(n, 1, 0.5); Treat <- rep(0:1, n / 2)
dfc <- data.frame(id = seq_len(n), Y = 2 * Treat * q1 + rnorm(n), Treat = Treat, q1 = q1, q2 = q2)
est <- function(d) list(estimate = unname(stats::coef(stats::lm(Y ~ Treat, data = d))[2]), se = NA_real_)
hsc <- data.table::data.table(grp = 1, K = 1, n = sum(q1), E = sum(q1), d1 = sum(q1 == 1 & Treat == 1),
  m1 = NA_real_, m0 = NA_real_, HR = 2, `L(HR)` = 1.5, `U(HR)` = 2.5, q1 = 1, q2 = 0)
f <- function(ct) suppressWarnings(subgroup.consistency(dfc, hsc, confs_labels = c("a", "b"),
  hr.threshold = 0, hr.consistency = 0, estimator_fn = est, consistency_threshold = ct,
  consistency_method = "split", n.splits = 50, parallel_args = list()))
p_null <- f(NULL)$n_passed; p_zero <- f(0)$n_passed
cat("consistency_threshold = NULL (documented: uses hr.consistency = 0) -> n_passed:", p_null, "\n")
cat("consistency_threshold = 0 explicitly -> n_passed:", p_zero, "\n")
cat("REPRO_VERDICT:", if (isTRUE(p_null == 0) && isTRUE(p_zero == 1)) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
consistency_threshold = NULL (documented: uses hr.consistency = 0) -> n_passed: 0 
consistency_threshold = 0 explicitly -> n_passed: 1 
REPRO_VERDICT: CONFIRMED
```

### V105

- Finding: A-078 — subgroup.consistency() neither documents nor checks that df must contain an id column; without one it returns a selected subgroup together with an all-zero sg.h
- Runs: exit 0, 2.3 s. Verdict: NOT_REPRODUCED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(11); n <- 400
q1 <- rbinom(n, 1, 0.5); q2 <- rbinom(n, 1, 0.5); Treat <- rep(0:1, n / 2)
Y <- rexp(n, rate = ifelse(q1 == 1 & Treat == 1, 5, 1)); Event <- rep(1L, n)
df <- data.frame(Y = Y, Event = Event, Treat = Treat, q1 = q1, q2 = q2)   # no 'id'
hs <- data.table::data.table(grp = 1, K = 1, n = sum(q1), E = sum(q1),
  d1 = sum(q1 == 1 & Treat == 1), m1 = 0.14, m0 = 0.69, HR = 4,
  `L(HR)` = 3, `U(HR)` = 5, q1 = 1, q2 = 0)
res <- subgroup.consistency(df, hs, confs_labels = c("q1 lab", "q2 lab"), parallel_args = list())
cat("selected:", res$sg.harm, " true subgroup size:", sum(q1), "\n")
cat("sum(sg.harm.id):", sum(res$sg.harm.id), " df_flag columns:", names(res$df_flag), "\n")
print(table(res$df_flag$treat.recommend))
ok <- !is.null(res$sg.harm) && sum(res$sg.harm.id) == 0 && !("id" %in% names(res$df_flag)) && all(res$df_flag$treat.recommend == 1)
cat("REPRO_VERDICT:", if (isTRUE(ok)) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
selected:  true subgroup size: 191 
sum(sg.harm.id): 0  df_flag columns: 
< table of extent 0 >
REPRO_VERDICT: NOT_REPRODUCED 
Warning message:
In value[[3L]](cond) :
  Error in sg_consistency_out for 'hr': arguments imply differing number of rows: 0, 400
```

### V106

- Finding: INV-2 — subgroup.consistency() calls set.seed(seed) on the global RNG (default 8316951) and never restores the caller's .Random.seed, so every forestsearch() call that 
- Runs: exit 0, 2.6 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(11); n <- 400
q1 <- rbinom(n, 1, 0.5); q2 <- rbinom(n, 1, 0.5); Treat <- rep(0:1, n / 2)
Y <- rexp(n, rate = ifelse(q1 == 1 & Treat == 1, 5, 1)); Event <- rep(1L, n)
df <- data.frame(id = seq_len(n), Y = Y, Event = Event, Treat = Treat, q1 = q1, q2 = q2)
hs <- data.table::data.table(grp = 1, K = 1, n = sum(q1), E = sum(q1),
  d1 = sum(q1 == 1 & Treat == 1), m1 = 0.14, m0 = 0.69, HR = 4,
  `L(HR)` = 3, `U(HR)` = 5, q1 = 1, q2 = 0)
run <- function() invisible(subgroup.consistency(df, hs,
  confs_labels = c("q1 lab", "q2 lab"), parallel_args = list()))
set.seed(1); run(); a <- runif(1)
set.seed(2); run(); b <- runif(1)
set.seed(8316951); ref <- runif(1)
print(c(a = a, b = b, ref = ref))
cat("equals first draw after set.seed(8316951):", identical(a, ref), "\n")
cat("REPRO_VERDICT:", if (identical(a, b)) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
        a         b       ref 
0.3564042 0.3564042 0.3564042 
equals first draw after set.seed(8316951): TRUE 
REPRO_VERDICT: CONFIRMED
```

### V107

- Finding: INV-2 — subgroup.consistency() calls set.seed(seed) on the session RNG and never restores the caller's RNG state, so every forestsearch() call that reaches the consiste
- Runs: exit 0, 2.4 s. Verdict: CONFIRMED (caller's RNG state not restored)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(11); n <- 400
q1 <- rbinom(n, 1, 0.5); q2 <- rbinom(n, 1, 0.5); Treat <- rep(0:1, n / 2)
Y <- rexp(n, rate = ifelse(Treat == 1, 0.8, 1))
df2 <- data.frame(id = seq_len(n), Y = Y, Event = 1L, Treat = Treat, q1 = q1, q2 = q2)
nn <- c(sum(q1), sum(q2), sum(q1 & q2))
hs3 <- data.table::data.table(grp = 1:3, K = c(1, 1, 2), n = nn, E = nn, d1 = c(90, 80, 40),
  m1 = c(0.5, 0.6, 0.7), m0 = c(0.6, 0.7, 0.8), HR = c(1.3, 1.2, 1.1),
  `L(HR)` = c(1, 0.9, 0.8), `U(HR)` = c(1.7, 1.6, 1.5), q1 = c(1, 0, 1), q2 = c(0, 1, 1))
set.seed(99); before <- .Random.seed
r <- try(suppressWarnings(subgroup.consistency(df2, hs3, confs_labels = c("a", "b"),
  pconsistency.threshold = 0, parallel_args = list(plan = "sequential"))), silent = TRUE)
after <- .Random.seed
set.seed(8316951); fixed <- .Random.seed
cat("call returned without error:", !inherits(r, "try-error"), "\n")
cat("state after call equals the state right after set.seed(8316951):", identical(after, fixed), "\n")
cat("REPRO_VERDICT:", if (!identical(before, after)) "CONFIRMED (caller's RNG state not restored)" else "NOT_REPRODUCED", "\n")
```

```text
call returned without error: TRUE 
state after call equals the state right after set.seed(8316951): TRUE 
REPRO_VERDICT: CONFIRMED (caller's RNG state not restored)
```

### V108

- Finding: A-067 — Under sg_focus = "maxeff" a finite m1.threshold is accepted and silently not applied, because the branch that skips the effect floor also skips the m1 filter, a
- Runs: exit 0, 2.6 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(11); n <- 400
q1 <- rbinom(n, 1, 0.5); q2 <- rbinom(n, 1, 0.5); Treat <- rep(0:1, n / 2)
Y <- rexp(n, rate = ifelse(q1 == 1 & Treat == 1, 5, 1)); Event <- rep(1L, n)
df <- data.frame(id = seq_len(n), Y = Y, Event = Event, Treat = Treat, q1 = q1, q2 = q2)
hs2 <- data.table::data.table(grp = 1:2, K = c(1, 1), n = c(sum(q1), sum(q2)),
  E = c(sum(q1), sum(q2)), d1 = c(100, 95), m1 = c(5, 0.1), m0 = c(0.7, 0.7),
  HR = c(4, 2), `L(HR)` = c(3, 1.5), `U(HR)` = c(5, 2.6), q1 = c(1, 0), q2 = c(0, 1))
g <- function(focus) suppressWarnings(subgroup.consistency(df, hs2, confs_labels = c("q1 lab", "q2 lab"),
  sg_focus = focus, m1.threshold = 1, parallel_args = list()))
r_me <- g("maxeff"); r_hr <- g("hr")
me <- paste(r_me$sg.harm, collapse = " & "); hr <- if (is.null(r_hr$sg.harm)) "none" else paste(r_hr$sg.harm, collapse = " & ")
cat("maxeff, m1.threshold = 1 -> selected:", me, "(candidate q1 has m1 = 5)\n")
cat("hr,     m1.threshold = 1 -> selected:", hr, "\n")
cat("REPRO_VERDICT:", if (grepl("q1", me, fixed = TRUE) && !grepl("q1", hr, fixed = TRUE)) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
maxeff, m1.threshold = 1 -> selected: {q1 lab} (candidate q1 has m1 = 5)
hr,     m1.threshold = 1 -> selected: {q2 lab} 
REPRO_VERDICT: CONFIRMED
```

### V109

- Finding: INV-5 — Under the default consistency_method = "resample", any candidate that is evaluated by literal splitting (closed form NA, or an unsupported measure such as IRD) 
- Runs: exit 0, 3.4 s. Verdict: CONFIRMED (Pcons differs between plans under consistency_method = 'resample')

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(11); n <- 400
q1 <- rbinom(n, 1, 0.5); q2 <- rbinom(n, 1, 0.5); Treat <- rep(0:1, n / 2)
Y <- rexp(n, rate = ifelse(Treat == 1, 1.25, 1))
df2 <- data.frame(id = seq_len(n), Y = Y, Event = 1L, Treat = Treat, q1 = q1, q2 = q2)
nn <- c(sum(q1), sum(q2), sum(q1 & q2))
hs3 <- data.table::data.table(grp = 1:3, K = c(1, 1, 2), n = nn, E = nn, d1 = c(90, 80, 40),
  m1 = c(0.5, 0.6, 0.7), m0 = c(0.6, 0.7, 0.8), HR = c(1.3, 1.2, 1.1),
  `L(HR)` = c(1, 0.9, 0.8), `U(HR)` = c(1.7, 1.6, 1.5), q1 = c(1, 0, 1), q2 = c(0, 1, 1))
est <- function(d) list(estimate = mean(d$Y[d$Treat == 1]) - mean(d$Y[d$Treat == 0]), se = 1)
# default method "resample" with an unsupported measure (IRD) -> literal splitting on both plans
go <- function(pa) {
  r <- suppressWarnings(subgroup.consistency(df2, hs3, confs_labels = c("a", "b"), sg_focus = "hr",
    pconsistency.threshold = 0, pconsistency.digits = 4, n.splits = 200,
    consistency_method = "resample", estimator_fn = est, consistency_threshold = -0.2,
    glm_resample_spec = list(effect_measure = "IRD"), parallel_args = pa))
  if (is.null(r$out_sg)) return(NULL)
  x <- as.data.frame(r$out_sg$result); x <- x[order(as.integer(x$m)), c("m", "Pcons")]; rownames(x) <- NULL; x
}
a <- go(list(plan = "sequential", workers = 1, show_message = FALSE))
b <- go(list(plan = "multicore", workers = 2, show_message = FALSE))
future::plan(future::sequential)
cat("sequential plan:\n"); print(a); cat("multicore plan, same seed:\n"); print(b)
cat("REPRO_VERDICT:", if (is.null(a) || is.null(b)) "INCONCLUSIVE (no result table)" else if (!isTRUE(all.equal(a, b))) "CONFIRMED (Pcons differs between plans under consistency_method = 'resample')" else "NOT_REPRODUCED", "\n")
```

```text
sequential plan:
  m Pcons
1 1 0.035
2 2 0.000
3 3 0.075
multicore plan, same seed:
  m Pcons
1 1  0.03
2 2  0.00
3 3  0.06
REPRO_VERDICT: CONFIRMED (Pcons differs between plans under consistency_method = 'resample')
```

### V110

- Finding: INV-5 — violated: the selection is plan-independent on the default path (consistency_method = "resample" with an estimable closed form, where no RNG is drawn), but when
- Runs: exit 0, 0.3 s. Verdict: inconclusive (a forestsearch run failed; see FAILED lines)

```r
## INV-5 reproducer: forestsearch() with parallel_args plan = "sequential" vs "multisession", same seed.
## Run from tempdir() after pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE).
## A = default consistency_method ("resample"): expected identical.
## B = consistency_method = "split": Pcons expected to differ; selection may differ.
## REPRO_VERDICT 'demonstrated' means the two plans SELECT DIFFERENTLY in at least one comparison.
t_start <- proc.time()[["elapsed"]]
old_plan <- future::plan()
had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
old_seed <- if (had_seed) get(".Random.seed", envir = globalenv(), inherits = FALSE) else NULL

inst <- suppressWarnings(tryCatch(utils::packageDescription("forestsearch", lib.loc = .libPaths()),
                                  error = function(e) NA))
if (is.list(inst)) {
  cat("installed forestsearch (the copy multisession workers load): Version", inst$Version,
      "| Packaged", inst$Packaged, "\n")
} else {
  cat("no installed forestsearch found in .libPaths(); multisession workers may fail (verdict would be 'inconclusive')\n")
}

make_df <- function(seed, n, hr_harm, two_planted) {
  set.seed(seed)
  b <- matrix(stats::rbinom(5L * n, 1L, 0.5), ncol = 5L)
  treat <- stats::rbinom(n, 1L, 0.5)
  harm <- b[, 1L] == 1L & b[, 2L] == 1L
  if (two_planted) harm <- harm | (b[, 3L] == 1L & b[, 4L] == 1L)
  rate <- 0.08 * exp(treat * ifelse(harm, log(hr_harm), log(0.8)))
  t_ev <- stats::rexp(n, rate)
  t_c <- stats::runif(n, 8, 30)
  data.frame(id = seq_len(n), tte = pmin(t_ev, t_c), event = as.integer(t_ev <= t_c),
             treat = treat, b1 = b[, 1L], b2 = b[, 2L], b3 = b[, 3L], b4 = b[, 4L],
             b5 = b[, 5L])
}

run_fs <- function(df, plan, method, pthr, splits, twostage) {
  w <- character(0)
  fs <- NULL
  junk <- utils::capture.output(
    fs <- withCallingHandlers(
      tryCatch(
        suppressMessages(forestsearch(
          df.analysis = df, outcome.name = "tte", event.name = "event",
          treat.name = "treat", id.name = "id",
          confounders.name = c("b1", "b2", "b3", "b4", "b5"),
          consistency_method = method, use_twostage = twostage,
          fs.splits = splits, pconsistency.threshold = pthr,
          n.min = 40, d0.min = 5, d1.min = 5, maxk = 2,
          details = FALSE, quiet = TRUE, seedit = 8316951,
          parallel_args = list(plan = plan, workers = 2, show_message = FALSE))),
        error = function(e) structure(list(error = conditionMessage(e)), class = "fs_error")),
      warning = function(cond) {
        w <<- c(w, conditionMessage(cond))
        invokeRestart("muffleWarning")
      }))
  list(fs = fs, warn = w)
}

failed <- function(r) {
  inherits(r$fs, "fs_error") ||
    any(grepl("Error in subgroup\\.consistency|Error in subgroup\\.search|Subgroup search failed|FSdata failure",
              r$warn))
}
sel_label <- function(r) {
  if (failed(r)) return("<FAILED>")
  if (is.null(r$fs$sg.harm)) return("<none>")
  paste(r$fs$sg.harm, collapse = " & ")
}
res_tab <- function(r) {
  if (inherits(r$fs, "fs_error")) return(NULL)
  x <- r$fs$grp.consistency$out_sg$result
  if (is.null(x)) return(NULL)
  x <- as.data.frame(x)
  rownames(x) <- NULL
  x
}
top_vals <- function(r) {
  x <- res_tab(r)
  if (is.null(x) || nrow(x) == 0L) return(NULL)
  unlist(x[1L, c("Pcons", "hr", "N", "E", "K")])
}
fmt <- function(v) {
  if (is.null(v)) return("<no row>")
  paste(sprintf("%s=%s", names(v), signif(v, 6)), collapse = " ")
}
pcons_by_cand <- function(r) {
  x <- res_tab(r)
  if (is.null(x)) return(numeric(0))
  stats::setNames(as.numeric(x$Pcons), as.character(x$m))
}
pcons_diff <- function(ra, rb) {
  pa <- pcons_by_cand(ra)
  pb <- pcons_by_cand(rb)
  ids <- union(names(pa), names(pb))
  sum(vapply(ids, function(i) !isTRUE(all.equal(unname(pa[i]), unname(pb[i]))), logical(1)))
}

compare <- function(label, df, method, pthr, splits, twostage) {
  rs <- run_fs(df, "sequential", method, pthr, splits, twostage)
  rp <- run_fs(df, "multisession", method, pthr, splits, twostage)
  fail <- failed(rs) || failed(rp)
  s_seq <- sel_label(rs)
  s_par <- sel_label(rp)
  cand <- function(r) if (inherits(r$fs, "fs_error")) NULL else r$fs$find.grps$out.found$hr.subgroups
  search_same <- isTRUE(all.equal(as.data.frame(cand(rs)), as.data.frame(cand(rp)),
                                  check.attributes = FALSE))
  est_same <- isTRUE(all.equal(top_vals(rs), top_vals(rp)))
  n_pc <- if (fail) NA_integer_ else pcons_diff(rs, rp)
  cat(sprintf("[%s] sequential  : selection = %s | top row: %s\n", label, s_seq, fmt(top_vals(rs))))
  cat(sprintf("[%s] multisession: selection = %s | top row: %s\n", label, s_par, fmt(top_vals(rp))))
  cat(sprintf("[%s] search table identical = %s | candidates whose Pcons/admission differ = %s | selection identical = %s | top-row estimates identical = %s\n",
              label, search_same, n_pc, identical(s_seq, s_par), est_same))
  if (fail) {
    msgs <- unique(c(rs$warn, rp$warn,
                     if (inherits(rs$fs, "fs_error")) rs$fs[["error"]],
                     if (inherits(rp$fs, "fs_error")) rp$fs[["error"]]))
    cat(sprintf("[%s] FAILED run(s): %s\n", label, paste(msgs, collapse = " || ")))
  }
  list(label = label, fail = fail,
       sel_differs = !fail && !identical(s_seq, s_par),
       est_differs = !fail && !est_same,
       n_pc = n_pc, search_same = search_same)
}

res <- tryCatch({
  A <- compare("A: default consistency_method='resample'",
               make_df(101L, 400L, 2.5, FALSE), "resample", 0.90, 100L, TRUE)
  B <- list()
  for (ds in 201:208) {
    if (proc.time()[["elapsed"]] - t_start > 65) break
    B[[length(B) + 1L]] <- compare(sprintf("B: consistency_method='split', data seed %d", ds),
                                   make_df(ds, 300L, 2.0, TRUE), "split", 0.60, 40L, FALSE)
    last <- B[[length(B)]]
    if (isTRUE(last$fail) || isTRUE(last$sel_differs)) break
  }
  list(A = A, B = B)
}, finally = {
  future::plan(old_plan)
  if (had_seed) assign(".Random.seed", old_seed, envir = globalenv())
})

flag <- function(lst, nm) vapply(lst, function(z) isTRUE(z[[nm]]), logical(1))
any_fail <- isTRUE(res$A$fail) || any(flag(res$B, "fail"))
b_sel <- flag(res$B, "sel_differs")
b_pc <- sum(vapply(res$B, function(z) isTRUE(z$n_pc > 0L), logical(1)))
cat("plan restored:", identical(class(future::plan()), class(old_plan)), "\n")
cat("RESAMPLE_DEFAULT_IDENTICAL (selection and top-row estimates):",
    !isTRUE(res$A$fail) && !isTRUE(res$A$sel_differs) && !isTRUE(res$A$est_differs), "\n")
cat("SPLIT comparisons run:", length(res$B),
    "| with different per-candidate Pcons/admission:", b_pc,
    "| with different selection:", sum(b_sel), "\n")
if (any_fail) {
  cat("REPRO_VERDICT: inconclusive (a forestsearch run failed; see FAILED lines)\n")
} else if (isTRUE(res$A$sel_differs) || any(b_sel)) {
  where <- c(if (isTRUE(res$A$sel_differs)) res$A$label,
             vapply(res$B[b_sel], function(z) z$label, character(1)))
  cat("REPRO_VERDICT: demonstrated (sequential and multisession select different subgroups for the same seed in:",
      paste(where, collapse = "; "), ")\n")
} else {
  cat("REPRO_VERDICT: not demonstrated (selection identical in all", 1L + length(res$B),
      "comparisons; split-path Pcons differed in", b_pc, "of", length(res$B), ")\n")
}
```

```text
installed forestsearch (the copy multisession workers load): Version 0.4.1.9002 | Packaged 2026-11-06 23:42:45 UTC; larryleon 
[A: default consistency_method='resample'] sequential  : selection = <FAILED> | top row: <no row>
[A: default consistency_method='resample'] multisession: selection = <FAILED> | top row: <no row>
[A: default consistency_method='resample'] search table identical = TRUE | candidates whose Pcons/admission differ = NA | selection identical = TRUE | top-row estimates identical = TRUE
[A: default consistency_method='resample'] FAILED run(s): could not find function "forestsearch"
[B: consistency_method='split', data seed 201] sequential  : selection = <FAILED> | top row: <no row>
[B: consistency_method='split', data seed 201] multisession: selection = <FAILED> | top row: <no row>
[B: consistency_method='split', data seed 201] search table identical = TRUE | candidates whose Pcons/admission differ = NA | selection identical = TRUE | top-row estimates identical = TRUE
[B: consistency_method='split', data seed 201] FAILED run(s): could not find function "forestsearch"
plan restored: TRUE 
RESAMPLE_DEFAULT_IDENTICAL (selection and top-row estimates): FALSE 
SPLIT comparisons run: 1 | with different per-candidate Pcons/admission: 0 | with different selection: 0 
REPRO_VERDICT: inconclusive (a forestsearch run failed; see FAILED lines)
```

### V111

- Finding: INV-2 — subgroup.search() unconditionally sets future::plan(multisession) and then future::plan(sequential) without saving the caller's plan or using on.exit, so every 
- Runs: exit 0, 3.4 s. Verdict: NOT_REPRODUCED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
old <- future::plan()
tryCatch({
  set.seed(3); n <- 200
  Z <- cbind(q1 = rbinom(n, 1, 0.5), q2 = rbinom(n, 1, 0.5), q3 = rbinom(n, 1, 0.5))
  Treat <- rep(0:1, n / 2); Y <- rexp(n); Event <- rbinom(n, 1, 0.8)
  future::plan(future::multisession, workers = 2)
  before <- class(future::plan())[1]
  r <- try(suppressWarnings(subgroup.search(Y = Y, Event = Event, Treat = Treat, Z = Z,
    n.min = 20, d0.min = 5, d1.min = 5, maxk = 1, parallel_workers = 1)), silent = TRUE)
  after <- class(future::plan())[1]
  cat("plan before:", before, " plan after:", after, " search errored:", inherits(r, "try-error"), "\n")
  cat("REPRO_VERDICT:", if (!identical(before, after)) "CONFIRMED" else "NOT_REPRODUCED", "\n")
}, finally = future::plan(old))
```

```text
plan before: FutureStrategy  plan after: FutureStrategy  search errored: FALSE 
REPRO_VERDICT: NOT_REPRODUCED
```

### V112

- Finding: INV-2 — subgroup.search() sets the future plan to multisession and then to sequential without saving the caller's plan and without on.exit, so every forestsearch() call
- Runs: exit 0, 3.3 s. Verdict: CONFIRMED (caller's multisession plan replaced by sequential)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 300
Y <- rexp(n); E <- rbinom(n, 1, 0.8); Tr <- rbinom(n, 1, 0.5)
Z <- cbind(z1 = rbinom(n, 1, 0.5), z2 = rbinom(n, 1, 0.5), z3 = rbinom(n, 1, 0.4))
future::plan(future::multisession, workers = 2)
wb <- future::nbrOfWorkers(); mb <- inherits(future::plan(), "multisession")
r <- subgroup.search(Y = Y, Event = E, Treat = Tr, Z = Z, n.min = 20, d0.min = 5, d1.min = 5,
                     hr.threshold = 0.5, maxk = 2, parallel_workers = 1)
wa <- future::nbrOfWorkers(); sa <- inherits(future::plan(), "sequential")
cat("workers before:", wb, " after:", wa, " plan sequential after:", sa, "\n")
cat("REPRO_VERDICT:", if (mb && wb == 2 && wa == 1) "CONFIRMED (caller's multisession plan replaced by sequential)" else "NOT_REPRODUCED", "\n")
future::plan(future::sequential)
```

```text
workers before: 2  after: 1  plan sequential after: TRUE 
REPRO_VERDICT: CONFIRMED (caller's multisession plan replaced by sequential)
```

### V113

- Finding: A-068 — With one candidate column the Z matrix collapses to a vector (no drop = FALSE) so the column count is NULL, and with fewer candidate columns than maxk combn() i
- Runs: exit 0, 2.5 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 300
Y <- rexp(n); E <- rbinom(n, 1, 0.8); Tr <- rbinom(n, 1, 0.5)
Zm <- cbind(z1 = rbinom(n, 1, 0.5), z2 = rbinom(n, 1, 0.5))
base <- list(Y = Y, Event = E, Treat = Tr, n.min = 20, d0.min = 5, d1.min = 5,
             hr.threshold = 0.5, parallel_workers = 1)
msg <- function(x) if (inherits(x, "try-error")) conditionMessage(attr(x, "condition")) else "ran"
one <- try(do.call(subgroup.search, c(base, list(Z = Zm[, 1, drop = FALSE], maxk = 1))), silent = TRUE)
two <- try(do.call(subgroup.search, c(base, list(Z = Zm, maxk = 3))), silent = TRUE)
ctl <- try(do.call(subgroup.search, c(base, list(Z = Zm, maxk = 2))), silent = TRUE)
cat("one column, maxk = 1 :", msg(one), "\n"); cat("two columns, maxk = 3:", msg(two), "\n")
cat("two columns, maxk = 2 (control):", msg(ctl), "\n")
cat("REPRO_VERDICT:", if (identical(msg(one), "n < m") && identical(msg(two), "n < m") && identical(msg(ctl), "ran")) "CONFIRMED" else "NOT_REPRODUCED", "\n")
future::plan(future::sequential)
```

```text
one column, maxk = 1 : n < m 
two columns, maxk = 3: n < m 
two columns, maxk = 2 (control): ran 
REPRO_VERDICT: CONFIRMED
```

### V114

- Finding: A-069 — When subgroup.search() is given estimator_fn (or adjust_covariates) without a row-aligned df_analysis, it silently scores candidates with the unadjusted treatme
- Runs: exit 0, 2.5 s. Verdict: CONFIRMED (estimator_fn ignored without error or warning; Cox HR returned)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 300
Tr <- rbinom(n, 1, 0.5); y <- rexp(n)
Zm <- cbind(z1 = rbinom(n, 1, 0.5), z2 = rbinom(n, 1, 0.5))
est <- function(d) list(estimate = 9, se = 1)   # sentinel: a row scored by it would carry 9
w <- character(0)
res <- withCallingHandlers(
  try(subgroup.search(Y = y, Event = rep(1L, n), Treat = Tr, Z = Zm, n.min = 20, d0.min = 5, d1.min = 5,
        hr.threshold = 0, maxk = 1, parallel_workers = 1, estimator_fn = est, outcome_type = "continuous"),
      silent = TRUE),
  warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
future::plan(future::sequential)
tab <- if (inherits(res, "try-error") || is.null(res$out.found)) NULL else as.data.frame(res$out.found$hr.subgroups)
if (!is.null(tab)) print(tab[, c("K", "n", "HR", "m1", "m0")])
cat("errored:", inherits(res, "try-error"), " warnings:", length(w), "\n")
ok <- !is.null(tab) && nrow(tab) > 0 && !any(abs(tab$HR - 9) < 1e-8) && !any(grepl("df_analysis|estimator", w))
cat("REPRO_VERDICT:", if (ok) "CONFIRMED (estimator_fn ignored without error or warning; Cox HR returned)" else "NOT_REPRODUCED", "\n")
```

```text
  K   n        HR        m1        m0
1 1 141 1.0202917 0.6748331 0.6510394
2 1 136 0.9475531 0.7377696 0.5284578
errored: FALSE  warnings: 0 
REPRO_VERDICT: CONFIRMED (estimator_fn ignored without error or warning; Cox HR returned)
```

### V115

- Finding: A-071 — format_results() builds a table row with c() and includes the optional fields by per-row is.na() tests, so the exported SG_tab_estimates() fails with its own de
- Runs: exit 0, 2.5 s. Verdict: CONFIRMED (default NULL names fail, named call works)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 200
df <- data.frame(tte = rexp(n, 0.05), event = rbinom(n, 1, 0.7), treat = rbinom(n, 1, 0.5), treat.recommend = rbinom(n, 1, 0.6))
f <- function(...) tryCatch(suppressWarnings(SG_tab_estimates(df, SG_flag = "treat.recommend", outcome.name = "tte", event.name = "event", treat.name = "treat", ...)), error = function(e) e)
r_named   <- f(sg1_name = "R", sg0_name = "Q")            # control: should work
r_default <- f()                                          # exported defaults, names NULL
r_onesided <- f(sg1_name = "R", sg0_name = "Q", hr_0a = "1.5 (1.0, 2.2)")
show <- function(x) if (inherits(x, "error")) paste("ERROR:", conditionMessage(x)) else paste("ok,", paste(dim(x), collapse = " x "))
cat("named   :", show(r_named), "\n"); cat("default :", show(r_default), "\n"); cat("onesided:", show(r_onesided), "\n")
cat("REPRO_VERDICT:", if (!inherits(r_named, "error") && inherits(r_default, "error")) "CONFIRMED (default NULL names fail, named call works)" else "NOT REPRODUCED", "\n")
```

```text
named   : ok, 2 x 8 
default : ERROR: length of 'dimnames' [2] not equal to array extent 
onesided: ERROR: length of 'dimnames' [2] not equal to array extent 
REPRO_VERDICT: CONFIRMED (default NULL names fail, named call works)
```

### V116

- Finding: A-072 — For GLM outcomes est.scale = "1/hr" reverses only the descriptive arm columns (n1, Rate/Mean(C), Rate/Mean(T), Diff) while the effect estimate still comes from 
- Runs: exit 0, 2.4 s. Verdict: CONFIRMED (arm columns reversed, OR column unchanged)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 400
df <- data.frame(treat = rbinom(n, 1, 0.5), treat.recommend = rbinom(n, 1, 0.6))
df$y <- rbinom(n, 1, 0.2 + 0.25 * df$treat)
est <- make_effect_estimator(outcome_type = "binary", treat.name = "treat", outcome.name = "y", effect_measure = "OR")
g <- function(sc) SG_tab_estimates_glm(df, SG_flag = "treat.recommend", outcome.name = "y", treat.name = "treat", estimator_fn = est, effect_measure = "OR", outcome_type = "binary", est.scale = sc, digits = 3)
a <- g("hr"); b <- g("1/hr"); print(a); print(b)
arms_swapped <- identical(a[["Rate(C)"]], b[["Rate(T)"]]) && !identical(a[["Rate(C)"]], b[["Rate(C)"]])
eff_same <- identical(a[["OR (95% CI)"]], b[["OR (95% CI)"]]) && !any(a[["OR (95% CI)"]] == "NA")
cat("REPRO_VERDICT:", if (arms_swapped && eff_same) "CONFIRMED (arm columns reversed, OR column unchanged)" else "NOT REPRODUCED", "\n")
```

```text
      Subgroup           n          n1 Rate(C) Rate(T)     Diff
1 Questionable 167 (41.8%)  76 (45.5%) 18.681% 51.316% +32.634%
2    Recommend 233 (58.2%) 107 (45.9%) 18.254% 40.187% +21.933%
           OR (95% CI)
1 4.588 (2.295, 9.175)
2 3.009 (1.660, 5.454)
      Subgroup           n          n1 Rate(C) Rate(T)     Diff
1 Questionable 167 (41.8%)  91 (54.5%) 51.316% 18.681% -32.634%
2    Recommend 233 (58.2%) 126 (54.1%) 40.187% 18.254% -21.933%
           OR (95% CI)
1 4.588 (2.295, 9.175)
2 3.009 (1.660, 5.454)
REPRO_VERDICT: CONFIRMED (arm columns reversed, OR column unchanged)
```

### V117

- Finding: A-073 — A subgroup that contains patients from one treatment arm only is not guarded: cox_summary() returns an NA string for it, but the following km_summary() and rmst
- Runs: exit 0, 2.3 s. Verdict: CONFIRMED (single-arm subgroup aborts the table)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(2); n <- 120
df <- data.frame(tte = rexp(n, 0.05), event = rbinom(n, 1, 0.7), treat = rbinom(n, 1, 0.5), treat.recommend = 1)
df$treat.recommend[which(df$treat == 1)[1:6]] <- 0   # sg0 = six treated patients, no controls
r <- tryCatch(suppressWarnings(SG_tab_estimates(df, SG_flag = "treat.recommend", outcome.name = "tte", event.name = "event", treat.name = "treat", sg1_name = "Recommend", sg0_name = "Questionable")), error = function(e) e)
if (inherits(r, "error")) cat("ERROR:", conditionMessage(r), "\n") else print(r)
cat("REPRO_VERDICT:", if (inherits(r, "error")) "CONFIRMED (single-arm subgroup aborts the table)" else "NOT REPRODUCED (table returned; inspect the single-arm row above)", "\n")
```

```text
ERROR: incorrect number of dimensions 
REPRO_VERDICT: CONFIRMED (single-arm subgroup aborts the table)
```

### V118

- Finding: A-079 — A standard set name that has no patient in the population (for example "Hstar" on a DGM built without an x3 subgroup) is not rejected, and ahr_limit() with b an
- Runs: exit 0, 2.5 s. Verdict: CONFIRMED empty named set returns a finite hr with no error; the same set as a logical vector errors

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); N <- 400
band <- rep(FALSE, N); region <- rbinom(N, 1, 0.3) == 1      # a population with no H* patients
pop <- structure(list(theta0 = rnorm(N, -3, 0.3), psi = rep(log(0.55), N), band = band, region = region,
  delay = rep(0, N), delays = 0, group = rep(1L, N), mu = 0, tau = 1,
  cell = paste0(ifelse(band, "band", "nonband"), "_", ifelse(region, "region", "nonregion")),
  x = 30, R = 6, N = N, cache = new.env()), class = "ahr_population")
e <- tryCatch(ahr_limit(pop, "Hstar", b = 0.55, h = 3), error = function(err) err)
lg <- tryCatch(ahr_limit(pop, pop$band, b = 0.55, h = 3), error = function(err) conditionMessage(err))
if (!inherits(e, "error")) print(c(hr = e$hr, exp_minus_10 = exp(-10), marghr_inf = e$marghr_inf, ahr = e$ahr, n_pop = e$n_pop)) else print(conditionMessage(e))
print(lg)
ok <- !inherits(e, "error") && isTRUE(e$n_pop == 0) && is.finite(e$hr) && identical(lg, "The set is empty.")
cat("REPRO_VERDICT:", if (ok) "CONFIRMED empty named set returns a finite hr with no error; the same set as a logical vector errors" else "NOT REPRODUCED", "\n")
```

```text
          hr exp_minus_10   marghr_inf          ahr        n_pop 
4.539993e-05 4.539993e-05 4.539993e-05          NaN 0.000000e+00 
[1] "The set is empty."
REPRO_VERDICT: CONFIRMED empty named set returns a finite hr with no error; the same set as a logical vector errors
```

### V119

- Finding: A-080 — ahr_design() does not validate b and h: with b = NULL (alone or with h = NULL) every set is computed on its own direct-route grid from the DGM's stored effects,
- Runs: exit 0, 2.6 s. Verdict: CONFIRMED: NULL b accepted, h = 7 ignored but echoed; max |corr diff| vs supported call = 0.441

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); N <- 600
band <- rbinom(N, 1, 0.2) == 1; region <- rbinom(N, 1, 0.3) == 1
pop <- structure(list(theta0 = rnorm(N, -3, 0.3), psi = log(0.55) + log(3 / 0.55) * band, band = band, region = region,
  delay = rep(0, N), delays = 0, group = rep(1L, N), mu = 0, tau = 1,
  cell = paste0(ifelse(band, "band", "nonband"), "_", ifelse(region, "region", "nonregion")),
  x = 30, R = 6, N = N, cache = new.env()), class = "ahr_population")
nw <- 0L
run <- function(...) tryCatch(withCallingHandlers(ahr_design(pop, ..., n = 800),
  warning = function(w) { nw <<- nw + 1L; invokeRestart("muffleWarning") }), error = function(e) e)
a <- run(b = 0.55, h = 3)            # supported call; pop$psi encodes the same b, h
z <- run(b = NULL, h = NULL); z2 <- run(b = NULL, h = 7)
if (is.list(z) && !inherits(z, "error")) { print(cbind(a$corr[, 1:2], supported = a$corr$exact, null_bh = z$corr$exact)); print(rbind(supported = a$prob, null_bh = z$prob)) } else print(conditionMessage(z))
cat("recycling or other warnings muffled:", nw, "\n")
v <- if (inherits(z, "error")) paste("PARTIAL: NULL b/h not validated, fails downstream with:", conditionMessage(z)) else if (!inherits(z2, "error") && isTRUE(all.equal(z2$prob, z$prob)) && identical(z2$h, 7)) sprintf("CONFIRMED: NULL b accepted, h = 7 ignored but echoed; max |corr diff| vs supported call = %.3g", max(abs(a$corr$exact - z$corr$exact))) else "NOT REPRODUCED"
cat("REPRO_VERDICT:", v, "\n")
```

```text
     big        small supported   null_bh
1    all       region 0.5502303 0.5415686
2    all region_Gstar 0.4873107 0.2514524
3 region region_Gstar 0.8853833 0.4448649
          p_itt_ub_lt1 p_region_itt_gt1 p_region_Gstar_lt1  p_trigger
supported    0.8945623       0.03778545          0.9991754 0.01987151
null_bh      0.8945623       0.03778545          0.9991754 0.02019224
          p_trigger_and_Gstar_lt1
supported              0.01961178
null_bh                0.02004399
recycling or other warnings muffled: 32 
REPRO_VERDICT: CONFIRMED: NULL b accepted, h = 7 ignored but echoed; max |corr diff| vs supported call = 0.441
```

### V120

- Finding: INV-2 — fs_build_eval_frame() reseeds the caller's global RNG through its non-NULL default eval_seed (set.seed inside the simulators) and never restores it, so the call
- Runs: exit 0, 2.4 s. Verdict: CONFIRMED caller's RNG stream after the call is the same for caller seeds 1 and 2

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
dgm <- suppressWarnings(setup_gbsg_dgm(model = "alt", k_inter = 1, n_super = 500L, seed = 1L))
r <- tryCatch({
  set.seed(1); invisible(fs_build_eval_frame(dgm, "survival")); a <- runif(3)
  set.seed(2); invisible(fs_build_eval_frame(dgm, "survival")); b <- runif(3)
  print(rbind(a, b)); identical(a, b)
}, error = function(e) { print(conditionMessage(e)); NA })
cat("REPRO_VERDICT:", if (isTRUE(r)) "CONFIRMED caller's RNG stream after the call is the same for caller seeds 1 and 2" else "NOT REPRODUCED", "\n")
```

```text
       [,1]      [,2]      [,3]
a 0.2739099 0.6535292 0.3595728
b 0.2739099 0.6535292 0.3595728
REPRO_VERDICT: CONFIRMED caller's RNG stream after the call is the same for caller seeds 1 and 2
```

### V121

- Finding: A-081 — For a fit with est.scale = "1/hr" the 'observed' and 'bias-corrected' reference lines of the bootstrap diagnostic histograms are drawn at the negative of their 
- Runs: exit 0, 2.5 s. Verdict: DEFECT_PRESENT (reference lines mirrored about zero under 1/hr)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1)
res <- data.frame(H_biasadj_2 = rnorm(200, log(2), 0.15),
                  Hc_biasadj_2 = rnorm(200, log(0.7), 0.10))
mk <- function(b0, b2) data.frame(H0 = ci_est(b0, 0.2, scale = "1/hr")$est,
                                  H2 = ci_est(b2, 0.2, scale = "1/hr")$est)
p <- create_bootstrap_diagnostic_plots(res, mk(log(2.2), log(2)), mk(log(0.65), log(0.7)),
                                       effect_label = "HR")
m <- mean(res$H_biasadj_2)
xi <- ggplot2::layer_data(p$H_distribution, 2)$xintercept
cat("H histogram mean (log HR):", round(m, 3), "| reference lines at:", round(xi, 3), "\n")
ok <- m > 0.5 && all(xi < -0.5)
cat("REPRO_VERDICT:", if (ok) "DEFECT_PRESENT (reference lines mirrored about zero under 1/hr)" else "NOT_REPRODUCED", "\n")
```

```text
H histogram mean (log HR): 0.698 | reference lines at: -0.693 -0.788 
REPRO_VERDICT: DEFECT_PRESENT (reference lines mirrored about zero under 1/hr)
```

### V122

- Finding: A-082 — The helper re-filters its input to rows with non-missing Pcons, so for DINA-selection bootstraps (Pcons always NA although any_found = 1) factor_presence and fa
- Runs: exit 0, 2.2 s. Verdict: DEFECT_PRESENT (found rows discarded by the Pcons re-filter)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
res <- data.frame(boot_id = 1:4, any_found = 1L, Pcons = NA_real_, hr_sg = 2,
  N_sg = 80L, K_sg = 1L,
  M.1 = c("{age <= 50}", "{age <= 50}", "{er <= 0}", "{age <= 50}"),
  M.2 = NA_character_, stringsAsFactors = FALSE)
sg_found <- res[res$any_found == 1L & !is.na(res$any_found), ]   # the caller's own filter
out <- summarize_factor_presence_robust(sg_found, maxk = 2, threshold = 10, as_gt = FALSE)
cat("rows passed as found:", nrow(sg_found), "| base_factors rows:", nrow(out$base_factors),
    "| specific_factors rows:", nrow(out$specific_factors), "\n")
ok <- nrow(sg_found) == 4 && nrow(out$base_factors) == 0 && nrow(out$specific_factors) == 0
cat("REPRO_VERDICT:", if (ok) "DEFECT_PRESENT (found rows discarded by the Pcons re-filter)" else "NOT_REPRODUCED", "\n")
```

```text
rows passed as found: 4 | base_factors rows: 0 | specific_factors rows: 0 
REPRO_VERDICT: DEFECT_PRESENT (found rows discarded by the Pcons re-filter)
```

### V123

- Finding: A-083 — Every forestsearch() error inside the calibration loop is discarded and counted as 'no subgroup found', so a fs_args list that makes every call fail returns a f
- Runs: exit 0, 2.4 s. Verdict: CONFIRMED (all fits errored, calibration object returned from clamp values, no failure reported)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
sim_null <- function(n, seed) { set.seed(seed)
  data.frame(id = seq_len(n), treat = rbinom(n, 1, 0.5), x1 = rnorm(n), tte = rexp(n), event = rbinom(n, 1, 0.7)) }
# fs_args repeats hr.threshold, which the wrapper also sets -> every forestsearch() call errors at argument matching
fs_args <- list(confounders.name = "x1", outcome.name = "tte", event.name = "event", treat.name = "treat", id.name = "id", hr.threshold = 1.25)
direct <- tryCatch(do.call(forestsearch, c(list(df.analysis = sim_null(100, 1)), fs_args, list(hr.threshold = 1.25, quiet = TRUE))), error = function(e) conditionMessage(e))
wmsg <- character(0)
cal <- withCallingHandlers(tryCatch(run_null_calibration(sim_null, fs_args, d_eff_fn = function(n) 0.7 * n, N_values = c(100, 200), n_sims = 10, verbose = FALSE), error = function(e) e), warning = function(w) { wmsg <<- c(wmsg, conditionMessage(w)); invokeRestart("muffleWarning") })
cat("direct forestsearch call:", if (is.character(direct)) direct else "no error", "\n")
cat("warnings raised by run_null_calibration:", if (length(wmsg)) paste(unique(wmsg), collapse = " | ") else "none", "\n")
if (inherits(cal, "leff_calibration")) print(cal$data[, c("N", "sim_fpr")]) else print(cal)
ok <- is.character(direct) && inherits(cal, "leff_calibration") && all(cal$data$sim_fpr == 0.5 / 10) && !any(grepl("forestsearch|fail|error", wmsg, ignore.case = TRUE))
cat("REPRO_VERDICT:", if (ok) "CONFIRMED (all fits errored, calibration object returned from clamp values, no failure reported)" else "NOT_CONFIRMED", "\n")
```

```text
direct forestsearch call: formal argument "hr.threshold" matched by multiple actual arguments 
warnings raised by run_null_calibration: none 
    N sim_fpr
1 100    0.05
2 200    0.05
REPRO_VERDICT: CONFIRMED (all fits errored, calibration object returned from clamp values, no failure reported)
```

### V124

- Finding: A-084 — A GRF result in which no subgroup was found is rejected with an input-type error, although the roxygen says non-detection is handled automatically by classifyin
- Runs: exit 0, 2.3 s. Verdict: CONFIRMED (GRF non-detection rejected with input-type error)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
d <- data.frame(id = 1:100, x = rnorm(100))
cfg <- list(maxdepth = 2, tau.rmst = 1, dmin.grf = 0, frac.tau = 1, n.min = 60, return_selected_cuts_only = FALSE)
grf_null <- create_null_result(d, NULL, NULL, cfg)   # the object grf.subg.harm.survival() returns when nothing is found
has_tr <- "treat.recommend" %in% names(grf_null$data)
r <- tryCatch(compare_fs_subgroups(c(rep(0L, 30), rep(1L, 70)), grf_null, label1 = "A", label2 = "GRF"), error = function(e) conditionMessage(e))
cat("treat.recommend in null GRF data:", has_tr, "\n"); print(if (is.character(r)) r else names(r))
cat("REPRO_VERDICT:", if (!has_tr && is.character(r) && grepl("must be a forestsearch object", r)) "CONFIRMED (GRF non-detection rejected with input-type error)" else "NOT_CONFIRMED", "\n")
```

```text
treat.recommend in null GRF data: FALSE 
[1] "fs2 must be a forestsearch object, a GRF result with data$treat.recommend, or a numeric/logical vector."
REPRO_VERDICT: CONFIRMED (GRF non-detection rejected with input-type error)
```

### V125

- Finding: A-085 — Result slots are filled with `list[[k]] <- value` where value can be NULL, which deletes element k rather than storing NULL, so the returned $fs, $errors, $ci_t
- Runs: exit 0, 2.3 s. Verdict: CONFIRMED ($fs shorter than n_combos; names shifted)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
# forestsearch() has no '...', so an unknown argument makes every fit fail at once (no fitting time)
out <- compare_selection_rules(df.analysis = data.frame(), sg_focus = c("hr", "hr", "hr"),
  selection_rule = c("pareto", "both", "neighborhood"), combo_labels = c("a", "b", "c"),
  compute_cis = FALSE, verbose = FALSE, not_a_forestsearch_arg = 1)
cat("n_combos:", nrow(out$combos), " length(out$fs):", length(out$fs), " names:", paste0("'", names(out$fs), "'", collapse = ", "), "\n")
# the same assignment pattern on the errors slot when every fit succeeds (err_k is NULL each time)
e <- vector("list", 3); names(e) <- c("a", "b", "c"); for (k in 1:3) e[[k]] <- NULL
cat("errors slot after 3 successes: length", length(e), "\n")
# and on the fs slot when only combo a fails
x <- vector("list", 3); names(x) <- c("a", "b", "c"); x[[1]] <- NULL; x[[2]] <- "fit_b"; x[[3]] <- "fit_c"
cat("x[['c']] holds:", x[["c"]], "\n")
ok <- length(out$fs) != nrow(out$combos) && length(e) != 3L && identical(x[["c"]], "fit_b")
cat("REPRO_VERDICT:", if (ok) "CONFIRMED ($fs shorter than n_combos; names shifted)" else "NOT_CONFIRMED", "\n")
```

```text
n_combos: 3  length(out$fs): 1  names: 'b' 
errors slot after 3 successes: length 1 
x[['c']] holds: fit_b 
REPRO_VERDICT: CONFIRMED ($fs shorter than n_combos; names shifted)
```

### V126

- Finding: A-086 — With plot_style = "grid", a threshold and CDE columns, the function stops in density() after all computation whenever one of the two subgroups has fewer than tw
- Runs: exit 0, 2.7 s. Verdict: CONFIRMED (grid style aborts in density(); combined style returns NaN rows)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(3); n <- 200
df <- data.frame(os_time = rexp(n, 0.05), os_event = rbinom(n, 1, 0.7), treat = rbinom(n, 1, 0.5), biomarker = rnorm(n, 50, 10))
df$theta_0 <- 0; df$theta_1 <- log(0.5); df$loghr_po <- df$theta_1 - df$theta_0   # every HR 0.5 < 0.7: Questionable is empty
grDevices::pdf(NULL)
run <- function(d, style) tryCatch(withCallingHandlers(cox_ahr_cde_analysis(d, hr_threshold = 0.7, plot_style = style, verbose = FALSE), warning = function(w) invokeRestart("muffleWarning")), error = function(e) paste("ERROR:", conditionMessage(e)))
r1 <- run(df, "grid"); r0 <- run(df, "combined")
df$theta_1 <- 0; df$loghr_po <- 0                                                 # nobody below 0.7: the max-biomarker fallback
r2 <- run(df, "grid")
grDevices::dev.off()
print(if (is.character(r1)) r1 else "grid/all-beneficial: returned"); print(if (is.character(r2)) r2 else "grid/none-beneficial: returned")
print(if (is.character(r0)) r0 else r0$subgroup_stats)
hit <- function(r) is.character(r) && grepl("at least 2 points", r)
cat("REPRO_VERDICT:", if (hit(r1) && hit(r2) && !is.character(r0)) "CONFIRMED (grid style aborts in density(); combined style returns NaN rows)" else "NOT_CONFIRMED", "\n")
```

```text
[1] "ERROR: need at least 2 points to select a bandwidth automatically"
[1] "ERROR: need at least 2 points to select a bandwidth automatically"
      Subgroup   n n_pct events AHR_po CDE
1 Questionable   0     0      0    NaN NaN
2    Recommend 200   100    134    0.5 0.5
3      Overall 200   100    134    0.5 0.5
REPRO_VERDICT: CONFIRMED (grid style aborts in density(); combined style returns NaN rows)
```

### V127

- Finding: A-088 — Records with a missing treatment value are not excluded before the arm split, so each one becomes an all-NA row in BOTH arm data frames, inflating the header n 
- Runs: exit 0, 2.4 s. Verdict: CONFIRMED (header n = 23 per arm, true 20; Missing rows although no covariate is missing)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1)
d <- data.frame(treat = c(rep(0:1, each = 20), NA, NA, NA), age = rnorm(43, 60, 10), bmi = rnorm(43, 27, 4))
out <- suppressMessages(create_summary_table(d, treat_var = "treat", vars_continuous = c("age", "bmi")))
is_gt <- inherits(out, "gt_tbl")
tab  <- if (is_gt) as.data.frame(out[["_data"]]) else out
labs <- if (is_gt) unlist(out[["_boxhead"]]$column_label) else names(out)
print(tab); print(labs)
cat("REPRO_VERDICT:", if (any(grepl("n=23", labs)) && sum(tab$Level == "Missing") == 2) "CONFIRMED (header n = 23 per arm, true 20; Missing rows although no covariate is missing)" else "NOT REPRODUCED", "\n")
```

```text
  Variable     Level    Control  Treatment P_value_formatted SMD_formatted
1      age Mean (SD) 61.9 (9.1) 59.9 (8.7)             0.489          0.22
2            Missing  3 (13.0%)  3 (13.0%)                                
3      bmi Mean (SD) 28.1 (3.8) 26.9 (3.7)             0.312          0.32
4            Missing  3 (13.0%)  3 (13.0%)                                
[1] "Characteristic"    ""                  "Control\n(n=23)"  
[4] "Treatment\n(n=23)" "P-value"           "SMD"              
REPRO_VERDICT: CONFIRMED (header n = 23 per arm, true 20; Missing rows although no covariate is missing)
```

### V128

- Finding: A-089 — The first sorted level of treat_var is always labelled "Control" and the second "Treatment", with the level values never shown, so a two-level coding whose cont
- Runs: exit 0, 2.3 s. Verdict: CONFIRMED ('Control' column shows the active arm, age mean ~70)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1)
d <- data.frame(arm = rep(c("placebo", "active"), each = 30), age = c(rnorm(30, 50, 2), rnorm(30, 70, 2)), bmi = rnorm(60, 27, 4))
out <- suppressMessages(create_summary_table(d, treat_var = "arm", vars_continuous = c("age", "bmi")))
tab <- if (inherits(out, "gt_tbl")) as.data.frame(out[["_data"]]) else out
print(tab)
ctrl_age <- as.numeric(sub(" .*", "", tab[[3]][1]))   # mean of age shown under 'Control'
cat("REPRO_VERDICT:", if (ctrl_age > 60) "CONFIRMED ('Control' column shows the active arm, age mean ~70)" else "NOT REPRODUCED", "\n")
```

```text
  Variable     Level    Control  Treatment P_value_formatted SMD_formatted
1      age Mean (SD) 70.3 (1.6) 50.2 (1.8)            <0.001         11.66
2      bmi Mean (SD) 27.5 (3.5) 27.4 (3.8)             0.990          0.00
REPRO_VERDICT: CONFIRMED ('Control' column shows the active arm, age mean ~70)
```

### V129

- Finding: A-090 — With the default use_alternating_rows = TRUE the even-row index is built as seq(2, nrow, by = 2), which errors ('wrong sign in by argument') when the table has 
- Runs: exit 0, 2.3 s. Verdict: CONFIRMED (one-row table errors in seq(2, 1, by = 2))

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1)
d <- data.frame(treat = rep(0:1, each = 20), age = rnorm(40, 60, 10))
r <- tryCatch(class(create_summary_table(d, treat_var = "treat", vars_continuous = "age")), error = function(e) paste("ERROR:", conditionMessage(e)))
print(r)
cat("REPRO_VERDICT:", if (any(grepl("wrong sign", r))) "CONFIRMED (one-row table errors in seq(2, 1, by = 2))" else "NOT REPRODUCED", "\n")
```

```text
[1] "ERROR: wrong sign in 'by' argument"
REPRO_VERDICT: CONFIRMED (one-row table errors in seq(2, 1, by = 2))
```

### V130

- Finding: A-091 — With include_raw = TRUE and table_style = "combined" the is.list(result) test is TRUE for a data.frame and for a gt table, so the raw matrices are assigned into
- Runs: exit 0, 2.3 s. Verdict: CONFIRMED (tenfold-shaped input errors; Kfold-shaped input returns a data.frame with a metadata list column)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
cv <- list(
  sens_summary = c(sens_H = .8, sens_Hc = .9, ppv_H = .7, ppv_Hc = .95),
  find_summary = c(Any = 1, Exact = .5, `At least 1` = .6, Cov1 = .6, Cov2 = .5,
                   `Cov 1 & 2` = .5, `Cov1 exact` = .5, `Cov2 exact` = .5),
  sens_out = matrix(0.5, 5, 4), find_out = matrix(0.5, 5, 8), sims = 5, Kfolds = 10)
r1 <- try(cv_metrics_tables(cv, use_gt = FALSE, include_raw = TRUE), silent = TRUE); print(r1)
cv1 <- cv; cv1$sens_out <- cv1$find_out <- NULL; cv1$sims <- NULL   # Kfold-shaped input
r2 <- cv_metrics_tables(cv1, use_gt = FALSE, include_raw = TRUE); str(r2, max.level = 1)
cat("REPRO_VERDICT:", if (inherits(r1, "try-error") && is.data.frame(r2) && "metadata" %in% names(r2)) "CONFIRMED (tenfold-shaped input errors; Kfold-shaped input returns a data.frame with a metadata list column)" else "NOT REPRODUCED", "\n")
```

```text
[1] "Error in `$<-.data.frame`(`*tmp*`, sens_out, value = structure(c(0.5,  : \n  replacement has 5 rows, data has 12\n"
attr(,"class")
[1] "try-error"
attr(,"condition")
<simpleError in `$<-.data.frame`(`*tmp*`, sens_out, value = structure(c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5), dim = 5:4)): replacement has 5 rows, data has 12>
'data.frame':	12 obs. of  5 variables:
 $ Category   : chr  "Agreement" "Agreement" "Agreement" "Agreement" ...
 $ Metric     : chr  "Sensitivity (H)" "Sensitivity (Hc)" "PPV (H)" "PPV (Hc)" ...
 $ Description: chr  "Agreement rate for subgroup H" "Agreement rate for complement Hc" "Positive predictive value for H" "Positive predictive value for Hc" ...
 $ Value      : num  80 90 70 95 100 50 60 60 50 50 ...
 $ metadata   :List of 12
REPRO_VERDICT: CONFIRMED (tenfold-shaped input errors; Kfold-shaped input returns a data.frame with a metadata list column)
```

### V131

- Finding: A-092 — dina_fit_bagged(seed = ) seeds only the bootstrap indices; the per-bag cross-fitting folds (on by default) are drawn unseeded from the caller's RNG stream, so t
- Runs: exit 0, 2.7 s. Verdict: CONFIRMED (seed fixes the bags but not the per-bag folds; RNG state advanced)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 200
X <- matrix(runif(n * 2, -1, 1), n, 2); colnames(X) <- c("x1", "x2")
W <- rbinom(n, 1, 0.5); Y <- X[, 1] + W * (0.5 + X[, 1]) + rnorm(n)
f <- function() suppressWarnings(dina_fit_bagged(X, Y, W, family = "gaussian", n_bags = 10L, seed = 1L))
f1 <- f(); f2 <- f()
same_idx  <- identical(f1$bag_inclusion, f2$bag_inclusion)
same_coef <- isTRUE(all.equal(coef(f1), coef(f2)))
s0 <- .Random.seed; invisible(f()); restored <- identical(s0, .Random.seed)
print(rbind(coef(f1), coef(f2)))
cat("same bootstrap indices:", same_idx, " same coefficients:", same_coef, " RNG restored:", restored, "\n")
cat("REPRO_VERDICT:", if (same_idx && !same_coef && !restored) "CONFIRMED (seed fixes the bags but not the per-bag folds; RNG state advanced)" else "NOT REPRODUCED", "\n")
```

```text
     (Intercept)       x1         x2
[1,]   0.5641430 1.079716 -0.3850158
[2,]   0.5690567 1.103241 -0.3028680
same bootstrap indices: TRUE  same coefficients: FALSE  RNG restored: FALSE 
REPRO_VERDICT: CONFIRMED (seed fixes the bags but not the per-bag folds; RNG state advanced)
```

### V132

- Finding: A-093 — dina_fit_bagged() with family = "cox" and Y supplied as a two-column matrix (a form its inherited @param Y documents and dina_fit() accepts) subsets Y by linear
- Runs: exit 0, 2.7 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 200
X <- matrix(runif(n * 2, -1, 1), n, 2, dimnames = list(NULL, c("x1", "x2")))
W <- rbinom(n, 1, 0.5)
Y <- cbind(time = rexp(n, exp(0.3 * W)), status = 1)
single_ok <- !inherits(try(dina_fit(X, Y, W, family = "cox", seed = 1L), silent = TRUE), "try-error")
run <- function(Yarg) suppressWarnings(try(dina_fit_bagged(X, Yarg, W, family = "cox", n_bags = 5L, seed = 1L), silent = TRUE))
r_mat <- run(Y); r_df <- run(as.data.frame(Y))
cat("dina_fit, matrix Y ok:", single_ok, "\n")
cat("dina_fit_bagged, matrix Y:", if (inherits(r_mat, "try-error")) conditionMessage(attr(r_mat, "condition")) else "ok", "\n")
cat("dina_fit_bagged, data-frame Y ok:", !inherits(r_df, "try-error"), "\n")
cat("REPRO_VERDICT:", if (single_ok && inherits(r_mat, "try-error") && !inherits(r_df, "try-error")) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
dina_fit, matrix Y ok: TRUE 
dina_fit_bagged, matrix Y: Fewer than 2 bags succeeded; cannot compute IJ variance. 
dina_fit_bagged, data-frame Y ok: TRUE 
REPRO_VERDICT: CONFIRMED
```

### V133

- Finding: A-094 — A supplied `sg` computed with tau_sign = -1 is accepted, but the bootstrap has no tau_sign: the fixed-subgroup effect CI is computed on the fit's own scale (opp
- Runs: exit 0, 2.8 s. Verdict: CONFIRMED (effect CI has the opposite sign from the point estimate)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(3); n <- 500
df <- data.frame(w = rbinom(n, 1, 0.5), x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
df$y <- df$w * (-0.6 - 0.8 * df$x1) + rnorm(n, sd = 0.5)
covs <- c("x1", "x2")
fit <- dina(df, outcome = "y", treatment = "w", covariates = covs, family = "gaussian", seed = 1L)
sg <- dina_subgroup(fit, df, covariates = covs, m_diff = 0.3, n_min = 60L, max_depth = 1L, tau_sign = -1)
stopifnot(isTRUE(sg$found))
bs <- suppressWarnings(dina_subgroup_bootstrap(df, outcome = "y", treatment = "w", covariates = covs,
  family = "gaussian", m_diff = 0.3, n_min = 60L, n_boot = 10L, seed = 1L, sg = sg, refit = FALSE))
cat("point mean_tau_hat:", sg$mean_tau_hat, " point Wald CI:", sg$ci, "\n")
cat("bootstrap effect_ci:", bs$effect_ci, "\n")
cat("resamples found:", bs$n_boot_found, "of", bs$n_boot, "\n")
opp <- all(is.finite(bs$effect_ci)) && sg$mean_tau_hat > 0 && all(bs$effect_ci < 0)
cat("REPRO_VERDICT:", if (opp) "CONFIRMED (effect CI has the opposite sign from the point estimate)" else "NOT_REPRODUCED", "\n")
```

```text
point mean_tau_hat: 0.6751932  point Wald CI: 0.5891678 0.7612187 
bootstrap effect_ci: -0.7197923 -0.6040229 
resamples found: 0 of 10 
REPRO_VERDICT: CONFIRMED (effect CI has the opposite sign from the point estimate)
```

### V134

- Finding: A-095 — dina_subgroup_bootstrap() has no max_depth, grid_probs or selection_rule argument and does not compare them with a supplied `sg`, so an `sg` found with max_dept
- Runs: exit 0, 2.8 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(3); n <- 500
df <- data.frame(w = rbinom(n, 1, 0.5), x1 = runif(n, -1, 1), x2 = runif(n, -1, 1))
df$y <- df$w * (0.2 + 0.9 * df$x1 + 0.9 * df$x2) + rnorm(n, sd = 0.5)
covs <- c("x1", "x2")
fit <- dina(df, outcome = "y", treatment = "w", covariates = covs, family = "gaussian", seed = 1L)
sg <- dina_subgroup(fit, df, covariates = covs, m_diff = 0.9, n_min = 60L, max_depth = 1L)
stopifnot(isTRUE(sg$found), sg$max_depth == 1L)
bs <- tryCatch(suppressWarnings(dina_subgroup_bootstrap(df, outcome = "y", treatment = "w", covariates = covs,
  family = "gaussian", m_diff = 0.9, n_min = 60L, n_boot = 10L, seed = 1L, sg = sg, refit = FALSE)), error = function(e) e)
no_arg <- !any(c("max_depth", "grid_probs", "selection_rule") %in% names(formals(dina_subgroup_bootstrap)))
accepted <- !inherits(bs, "error")
depth2 <- accepted && any(grepl(" & ", bs$boot_results$covariate, fixed = TRUE), na.rm = TRUE)
if (accepted) print(bs$boot_results[, c("found", "covariate", "direction", "n_subgroup")])
cat("point depth:", sg$depth, "| sg$max_depth:", sg$max_depth, "| accepted:", accepted, "| depth-2 resample selections:", depth2, "\n")
cat("REPRO_VERDICT:", if (no_arg && accepted && depth2) "CONFIRMED" else if (no_arg && accepted) "PARTIAL (accepted, no depth-2 resample seen)" else "NOT_REPRODUCED", "\n")
```

```text
   found covariate     direction n_subgroup
1   TRUE   x1 & x2 right & right        149
2   TRUE   x1 & x2 right & right        127
3   TRUE   x1 & x2 right & right        117
4   TRUE   x1 & x2 right & right        146
5   TRUE   x1 & x2 right & right        123
6   TRUE   x1 & x2 right & right        155
7   TRUE   x1 & x2 right & right        154
8   TRUE   x1 & x2 right & right        137
9   TRUE   x1 & x2 right & right        175
10  TRUE   x1 & x2 right & right        130
point depth: 1 | sg$max_depth: 1 | accepted: TRUE | depth-2 resample selections: TRUE 
REPRO_VERDICT: CONFIRMED
```

### V135

- Finding: A-097 — The ratio-scale line of print.dina_subgroup_refit() is hard-coded to read "95% CI" although the limits printed are the exponentiated 1 - alpha limits, so any al
- Runs: exit 0, 2.3 s. Verdict: CONFIRMED (90% limits printed under a 95% label)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
z <- qnorm(0.95)
x <- structure(list(family = "cox", signature = "{x1 >= 0}", n_subgroup = 100L, n_treated = 50L, n_control = 50L,
  n_events = 40L, n_events_treated = 22L, n_events_control = 18L, confounders_used = character(0),
  strata_used = character(0), alpha = 0.10, effect = 0.5, se = 0.2,
  ci = c(lower = 0.5 - z * 0.2, upper = 0.5 + z * 0.2), effect_scale = "log-HR", ratio_scale = TRUE),
  class = "dina_subgroup_refit")
out <- capture.output(print(x))
cat(out, sep = "\n")
hit <- any(grepl("Wald 90% CI", out, fixed = TRUE)) && any(grepl("95% CI (HR)", out, fixed = TRUE))
cat("REPRO_VERDICT:", if (hit) "CONFIRMED (90% limits printed under a 95% label)" else "NOT_REPRODUCED", "\n")
```

```text
Standard within-subgroup model (DINA-identified subgroup)
  Family:       cox
  Signature:    {x1 >= 0}
  Subgroup n:   100 (treated 50, control 50)
  Events:       40 (treated 22, control 18)
  Adjustment:   unadjusted

  Treatment effect (log-HR):  0.5  (SE 0.2)
  Wald 90% CI (log-HR):  [0.171, 0.829]
  HR:  1.649   95% CI (HR):  [1.187, 2.291]

  Note: CI is conditional on the discovered signature treated as
  pre-specified; it does NOT adjust for signature selection.
REPRO_VERDICT: CONFIRMED (90% limits printed under a 95% label)
```

### V136

- Finding: A-098 — The loser classification and its wording hard-code the hrMaxSG/hrMinSG key order (band, N, Pcons, effect, K), so under sg_focus = "maxeff"/"maxeffCons" a loser 
- Runs: exit 0, 2.4 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
library(data.table)
res <- data.table(Pcons = c(0.99, 0.93), hr = c(2.5, 2.0), N = c(100, 200), E = c(50, 90),
                  g = "1", m = c("1", "2"), K = c(1, 1), M.1 = c("{b <= 1}", "{c <= 1}"), M.2 = "")
mk <- function(f) structure(list(grp.consistency = list(out_sg = list(result = res, pareto_frontier = res)),
  args_call_all = list(sg_focus = f, selection_rule = "neighborhood", effect_neighborhood = 0.10)),
  class = c("forestsearch", "list"))
o1 <- explain_pareto_selection(mk("maxeffCons"), format = "list")
o2 <- explain_pareto_selection(mk("hr"), format = "list")
cat("maxeffCons loss axis:", o1$losers[[1]]$loss_axis, "\n ", o1$losers[[1]]$plain, "\n")
cat("hr loss axis:", o2$losers[[1]]$loss_axis, "\n ", o2$losers[[1]]$plain, "\n")
hit <- identical(o1$losers[[1]]$loss_axis, "Pcons") &&
  grepl("secondary tiebreaker after N", o2$losers[[1]]$plain, fixed = TRUE)
cat("REPRO_VERDICT:", if (hit) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text

Attaching package: ‘data.table’

The following object is masked from ‘package:base’:

    %notin%

maxeffCons loss axis: Pcons 
  - **{c <= 1}** (hazard ratio = 2.000, N = 200, Pcons = 0.930, K = 1). Lower consistency (Pcons=0.930 vs 0.990); Pcons is the secondary tiebreaker after N. 
hr loss axis: Pcons 
  - **{c <= 1}** (hazard ratio = 2.000, N = 200, Pcons = 0.930, K = 1). Lower consistency (Pcons=0.930 vs 0.990); Pcons is the secondary tiebreaker after N. 
REPRO_VERDICT: CONFIRMED
```

### V137

- Finding: A-099 — When the selected subgroup is not on the (effect, N) frontier, explain_pareto_selection() attributes this to `hrMinSG` whatever the sg_focus, promises frontier 
- Runs: exit 0, 2.4 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
library(data.table)
res <- data.table(Pcons = c(0.99, 0.95, 0.93), hr = c(1.5, 2.5, 2.0), N = c(80, 100, 200), E = c(40, 50, 90),
                  g = "1", m = c("3", "1", "2"), K = c(1, 1, 1),
                  M.1 = c("{a <= 1}", "{b <= 1}", "{c <= 1}"), M.2 = "")
# top row (highest Pcons, as sg_focus = "hr" sorts) is dominated on (hr, N); frontier = rows 2 and 3
fs <- structure(list(grp.consistency = list(out_sg = list(result = res, pareto_frontier = res[2:3])),
                     args_call_all = list(sg_focus = "hr", selection_rule = "neighborhood", effect_neighborhood = 0.10)),
                class = c("forestsearch", "list"))
o <- explain_pareto_selection(fs, format = "list")
txt <- explain_pareto_selection(fs)
cat(txt, "\n\n")
cat("n_frontier:", o$meta$n_frontier, "| losers listed:", length(o$losers), "\n")
hit <- o$meta$n_frontier == 2L && length(o$losers) == 0L && grepl("allowed under `hrMinSG`", txt, fixed = TRUE) &&
  grepl("frontier has a single candidate", txt, fixed = TRUE)
cat("REPRO_VERDICT:", if (hit) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text

Attaching package: ‘data.table’

The following object is masked from ‘package:base’:

    %notin%

**Selection criterion.** The configured criterion is highest-consistency (sg_focus = "hr") under the effect-size neighborhood rule.

**Selected subgroup.** The selected subgroup does not appear on the Pareto frontier shown above. This is allowed under `hrMinSG`, which may pick a small-N candidate that is itself dominated by a larger-N candidate with similar effect. The frontier rows below describe the non-dominated alternatives that were evaluated but not chosen.

**Non-selected candidates.** The frontier has a single candidate, so no comparison applies. 

n_frontier: 2 | losers listed: 0 
REPRO_VERDICT: CONFIRMED
```

### V138

- Finding: INV-2 — forestsearch_Kfold() calls set.seed(seedit) in the caller's session when Kfolds < n and never restores the caller's .Random.seed, violating invariant 2.
- Runs: exit 0, 5.7 s. Verdict: CONFIRMED (caller RNG state after the call does not depend on the caller's seed)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(42); N <- 200L
age <- round(rnorm(N, 55, 12)); bmk <- rnorm(N); sex <- rbinom(N, 1L, 0.5); treat <- rbinom(N, 1L, 0.5)
tt <- rexp(N, 0.05 * exp(ifelse(age <= 50, log(2.5), 0) * treat))
df <- data.frame(id = seq_len(N), treat = treat, time = pmin(tt, 60), event = as.integer(tt <= 60), age = age, bmk = bmk, sex = sex)
pa <- list(plan = "sequential", workers = 1L, show_message = FALSE)
fs <- suppressWarnings(forestsearch(df, outcome.name = "time", event.name = "event", treat.name = "treat", id.name = "id",
  confounders.name = c("age", "bmk", "sex"), use_grf = FALSE, use_lasso = FALSE, hr.threshold = 1.25, hr.consistency = 1.0,
  pconsistency.threshold = 0.80, maxk = 1L, n.min = 40L, d0.min = 5L, d1.min = 5L, fs.splits = 50L, seedit = 42L,
  details = FALSE, quiet = TRUE, parallel_args = pa))
set.seed(1);   invisible(suppressWarnings(forestsearch_Kfold(fs, Kfolds = 3L, parallel_args = pa))); s1 <- .Random.seed; r1 <- runif(3)
set.seed(999); invisible(suppressWarnings(forestsearch_Kfold(fs, Kfolds = 3L, parallel_args = pa))); s2 <- .Random.seed; r2 <- runif(3)
print(rbind(r1, r2))
same <- identical(s1, s2) && identical(r1, r2)
cat("REPRO_VERDICT:", if (same) "CONFIRMED (caller RNG state after the call does not depend on the caller's seed)" else "NOT_REPRODUCED", "\n")
```

```text
        [,1]      [,2]      [,3]
r1 0.4463279 0.8574682 0.6407764
r2 0.4463279 0.8574682 0.6407764
REPRO_VERDICT: CONFIRMED (caller RNG state after the call does not depend on the caller's seed)
```

### V139

- Finding: A-100 — A fold refit that fails inside forestsearch() at a stage forestsearch() itself converts to a warning plus an sg.harm = NULL return (get_FSdata, subgroup.search,
- Runs: exit 0, 3.8 s. Verdict: CONFIRMED (all 3 failed refits recorded as ITT no-subgroup folds, no error)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(42); N <- 200L
age <- round(rnorm(N, 55, 12)); bmk <- rnorm(N); sex <- rbinom(N, 1L, 0.5); treat <- rbinom(N, 1L, 0.5)
tt <- rexp(N, 0.05 * exp(ifelse(age <= 50, log(2.5), 0) * treat))
df <- data.frame(id = seq_len(N), treat = treat, time = pmin(tt, 60),
                 event = as.integer(tt <= 60), age = age, bmk = bmk, sex = sex)
pa <- list(plan = "sequential", workers = 1L, show_message = FALSE)
fs <- suppressWarnings(forestsearch(df, outcome.name = "time", event.name = "event",
  treat.name = "treat", id.name = "id", confounders.name = c("age", "bmk", "sex"),
  use_grf = FALSE, use_lasso = FALSE, hr.threshold = 1.25, hr.consistency = 1.0,
  pconsistency.threshold = 0.80, maxk = 1L, n.min = 40L, d0.min = 5L, d1.min = 5L,
  fs.splits = 50L, seedit = 42L, details = FALSE, quiet = TRUE, parallel_args = pa))
# Every training refit now fails in data preparation; by the stated policy CV should halt.
kf <- tryCatch(testthat::with_mocked_bindings(
  suppressWarnings(forestsearch_Kfold(fs, Kfolds = 3L, parallel_args = pa)),
  get_FSdata = function(...) stop("simulated fold-level data-prep failure"),
  .package = "forestsearch"), error = function(e) e)
if (inherits(kf, "error")) {
  cat("Kfold raised:", conditionMessage(kf), "\n")
  cat("REPRO_VERDICT: REFUTED (CV halted on the failed refits)\n")
} else {
  cat("prop_SG_found =", kf$prop_SG_found, "\n"); print(table(kf$resCV$treat.recommend, useNA = "ifany"))
  ok <- all(kf$resCV$treat.recommend == 1) && isTRUE(kf$prop_SG_found == 0)
  cat("REPRO_VERDICT:", if (ok) "CONFIRMED (all 3 failed refits recorded as ITT no-subgroup folds, no error)" else "INCONCLUSIVE", "\n")
}
```

```text
prop_SG_found = 0 

  1 
200 
REPRO_VERDICT: CONFIRMED (all 3 failed refits recorded as ITT no-subgroup folds, no error)
```

### V140

- Finding: A-102 — For a survival fit with est.scale = "1/hr", the two subgroup-table calls in forestsearch_KfoldOut(outall = TRUE) omit est.scale, so the CV subgroup rows are haz
- Runs: exit 0, 2.6 s. Verdict: CONFIRMED (est.scale = 1/hr requested, CV subgroup rows are on the hr scale)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 240
d <- data.frame(id = 1:n, treat = rbinom(n, 1, .5))
d$time <- rexp(n, 0.1 * exp(0.5 * d$treat)); d$event <- rbinom(n, 1, .8)
d$treat.recommend.original <- rbinom(n, 1, .5); d$treat.recommend <- d$treat.recommend.original
d$cvindex <- rep(1:4, length.out = n); d$sg1 <- "{age <= 50}"; d$sg2 <- NA_character_
res <- list(cv_args = list(confounders.name = "age", outcome.name = "time", event.name = "event",
                           treat.name = "treat", est.scale = "1/hr"),
            sg0.name = "Recommend", sg1.name = "Not recommend", Kfolds = 4L,
            resCV = d, sg_analysis = "{age <= 50}")
out <- forestsearch_KfoldOut(res, outall = TRUE)
ref <- function(sc) as.character(unlist(SG_tab_estimates(d, SG_flag = "treat.recommend", outcome.name = "time",
  event.name = "event", treat.name = "treat", sg0_name = "Recommend", sg1_name = "Not recommend",
  est.scale = sc)[, "HR (95% CI)"]))
got <- as.character(unlist(out$SG_tab_Kfold[, "HR (95% CI)"]))
cat("CV table   :", got, "\n"); cat("hr scale   :", ref("hr"), "\n"); cat("1/hr scale :", ref("1/hr"), "\n")
cat("REPRO_VERDICT:", if (identical(got, ref("hr")) && !identical(got, ref("1/hr")))
  "CONFIRMED (est.scale = 1/hr requested, CV subgroup rows are on the hr scale)" else "REFUTED", "\n")
```

```text
CV table   : 1.98 (1.30, 3.01) 1.66 (1.10, 2.51) 
hr scale   : 1.98 (1.30, 3.01) 1.66 (1.10, 2.51) 
1/hr scale : 0.50 (0.33, 0.77) 0.60 (0.40, 0.91) 
REPRO_VERDICT: CONFIRMED (est.scale = 1/hr requested, CV subgroup rows are on the hr scale)
```

### V141

- Finding: A-103 — When the original subgroup has one factor, the CV 'Exact' find metric scores a fold as an exact match whenever either fold factor equals that factor, so a fold 
- Runs: exit 0, 2.3 s. Verdict: CONFIRMED (two-factor folds scored as exact matches of a one-factor subgroup)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
out <- CV_sgs(sg1 = c("{age <= 50}", "{age <= 50}", "{sex >= 1}"), sg2 = c(NA, "{sex >= 1}", "{age <= 50}"),
              confs = c("age", "sex"), sg_analysis = "{age <= 50}")
print(out$exact_match)   # folds 2 and 3 selected {age <= 50} & {sex >= 1}, not the original subgroup
strict <- .cv_canonicalize_sg(c("{age <= 50}", "{age <= 50}", "{sex >= 1}"), c(NA, "{sex >= 1}", "{age <= 50}")) == "{age <= 50}"
print(as.numeric(strict))
cat("REPRO_VERDICT:", if (identical(as.numeric(out$exact_match), c(1, 1, 1)) && identical(as.numeric(strict), c(1, 0, 0)))
  "CONFIRMED (two-factor folds scored as exact matches of a one-factor subgroup)" else "REFUTED", "\n")
```

```text
[1] 1 1 1
[1] 1 0 0
REPRO_VERDICT: CONFIRMED (two-factor folds scored as exact matches of a one-factor subgroup)
```

### V142

- Finding: A-104 — The covariate-any match misattributes covariates when one confounder name is a prefix or substring of another: the disambiguating substr() is one character shor
- Runs: exit 0, 2.2 s. Verdict: CONFIRMED (substring and off-by-one misattribution in covariate-any match)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
# (a) original subgroup cuts on age; the folds cut on stage only
a <- CV_sgs(sg1 = c("{stage <= 2}", "{stage <= 1}", NA), sg2 = c(NA, NA, NA),
            confs = c("age", "stage", "sex"), sg_analysis = "{age <= 50}")$cov1_any
print(a)   # expected 0 0 0
# (b) target on x11, fold cut on x1 only
b <- find_covariate_any_match("{x11 <= 2}", sg1 = "{x1 <= 0}", sg2 = NA, confs = c("x1", "x11"))
print(b)   # expected 0
# (c) negated target, fold reproduces the identical label
c3 <- find_covariate_any_match("!{x11 <= 2}", sg1 = "!{x11 <= 2}", sg2 = NA, confs = c("x1", "x11"))
print(c3)  # expected 1
cat("REPRO_VERDICT:", if (identical(as.numeric(a), c(1, 1, 0)) && b == 1 && c3 == 0)
  "CONFIRMED (substring and off-by-one misattribution in covariate-any match)" else "PARTIAL/REFUTED - inspect printed values", "\n")
```

```text
[1] 1 1 0
[1] 1
[1] 0
REPRO_VERDICT: CONFIRMED (substring and off-by-one misattribution in covariate-any match)
```

### V143

- Finding: A-105 — The line wrapper inside format_fs_details() never terminates when, after an indent, the next token has no space or comma and is longer than the remaining width:
- Runs: exit 0, 7.3 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
ok  <- c("Candidate factors: 1", "  a short line")
bad <- c("Candidate factors: 1", paste0("  ", strrep("x", 60)))
r0 <- tryCatch({ invisible(utils::capture.output(format_fs_details(ok))); "returned" }, error = function(e) conditionMessage(e))
r1 <- tryCatch({ setTimeLimit(elapsed = 5, transient = TRUE)
                 invisible(utils::capture.output(format_fs_details(bad))); "returned" },
               error = function(e) conditionMessage(e), finally = setTimeLimit(elapsed = Inf))
cat("control input:", r0, "| long unbroken token:", r1, "\n")
cat("REPRO_VERDICT:", if (identical(r0, "returned") && grepl("time limit", r1)) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
control input: returned | long unbroken token: reached elapsed time limit 
REPRO_VERDICT: CONFIRMED
```

### V144

- Finding: A-106 — fpr_approx() requires effect_true > 0 and c1 > 0 before it looks at effect_type, so effect_type = "difference" stops on a zero or negative true difference (and 
- Runs: exit 0, 2.3 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
f <- function(...) tryCatch(fpr_approx(..., effect_type = "difference", verbose = FALSE), error = function(e) conditionMessage(e))
r0 <- f(effect_true = 0, c1 = 0.05)       # null difference
r1 <- f(effect_true = -0.05, c1 = 0.05)   # protective difference
r2 <- f(effect_true = 0.05, c1 = 0.10)    # control: positive difference is accepted
print(r0); print(r1); print(class(r2))
cat("REPRO_VERDICT:", if (is.character(r0) && is.character(r1) && inherits(r2, "fpr_approximation")) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
[1] "effect_true > 0 is not TRUE"
[1] "effect_true > 0 is not TRUE"
[1] "fpr_approximation"
REPRO_VERDICT: CONFIRMED
```

### V145

- Finding: A-108 — fpr_calibration() always computes the per-subgroup approximation P1 on the ratio scale and requires c1 > 0 and c2 > 0, so for identity-scale measures (RD, IRD, 
- Runs: exit 0, 2.8 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
future::plan(future::sequential)
set.seed(1)
df <- data.frame(id = 1:200, treat = rbinom(200, 1, 0.5), x1 = rnorm(200), y = rbinom(200, 1, 0.3))
fsp <- list(confounders.name = "not_a_column", outcome.name = "y", event.name = "y", treat.name = "treat",
            id.name = "id", outcome_type = "binary", effect_measure = "RD")
e0 <- tryCatch(fpr_calibration(c1 = 0.05, c2 = 0, df.analysis = df, fs_params = fsp, nsims = 1, quiet = TRUE),
               error = function(e) conditionMessage(e))
res <- suppressWarnings(fpr_calibration(c1 = 0.05, c2 = 0.01, df.analysis = df, fs_params = fsp, nsims = 1, quiet = TRUE))
p_id <- compute_detection_probability_glm(theta = 0, d_eff = res$d_eff, c1 = 0.05, c2 = 0.01,
                                          effect_scale = "difference", method = "monte_carlo")
print(e0); print(c(P1_reported = res$P1, P1_identity_scale = p_id, L_eff_implied = res$L_eff_implied))
cat("REPRO_VERDICT:", if (is.character(e0) && isTRUE(abs(res$P1 - p_id) > 0.1)) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
[1] "c2 must be a positive scalar"
      P1_reported P1_identity_scale     L_eff_implied 
          1.00000           0.24378                NA 
REPRO_VERDICT: CONFIRMED
```

### V146

- Finding: A-109 — fpr_calibration() injects c1/c2 as hr.threshold/hr.consistency without removing effect.threshold, consistency.threshold or seedit from fs_params, so an fs_param
- Runs: exit 0, 4.4 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
future::plan(future::sequential)
set.seed(1)
df <- data.frame(id = 1:200, treat = rbinom(200, 1, 0.5), x1 = rnorm(200), y = rbinom(200, 1, 0.3))
base <- list(confounders.name = "x1", outcome.name = "y", event.name = "y", treat.name = "treat", id.name = "id",
             outcome_type = "binary", effect_measure = "OR",
             parallel_args = list(plan = "sequential", workers = 1, show_message = FALSE))
run <- function(extra) suppressWarnings(fpr_calibration(c1 = 1.7, c2 = 1.7, df.analysis = df,
                         fs_params = c(base, extra), nsims = 2, quiet = TRUE))
r0 <- run(list()); r1 <- run(list(effect.threshold = 1.2)); r2 <- run(list(seedit = 123))
print(rbind(control = c(failed = r0$n_failed, done = r0$nsims), new_spelling = c(r1$n_failed, r1$nsims),
            seedit = c(r2$n_failed, r2$nsims)))
cat("REPRO_VERDICT:", if (r0$n_failed != 0) "INCONCLUSIVE_CONTROL_FAILED" else if (r1$n_failed == 2 && r2$n_failed == 2) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
             failed done
control           0    2
new_spelling      2    0
seedit            2    0
REPRO_VERDICT: CONFIRMED
```

### V147

- Finding: INV-2 — With n_workers > 1, fpr_calibration() sets a multisession future plan and on exit sets plan("sequential") instead of restoring the plan the caller had (invarian
- Runs: exit 0, 5.0 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
future::plan(future::multisession, workers = 2)
before <- inherits(future::plan(), "multisession")
set.seed(1)
df <- data.frame(id = 1:200, treat = rbinom(200, 1, 0.5), x1 = rnorm(200), y = rbinom(200, 1, 0.3))
fsp <- list(confounders.name = "x1", outcome.name = "y", event.name = "y", treat.name = "treat", id.name = "id",
            outcome_type = "binary", effect_measure = "OR")
invisible(try(suppressWarnings(fpr_calibration(c1 = 1.7, c2 = 1.7, df.analysis = df, fs_params = fsp,
                                                nsims = 1, n_workers = 2, quiet = TRUE)), silent = TRUE))
after <- inherits(future::plan(), "multisession")
future::plan(future::sequential)
print(c(caller_plan_multisession_before = before, after = after))
cat("REPRO_VERDICT:", if (before && !after) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
caller_plan_multisession_before                           after 
                           TRUE                           FALSE 
REPRO_VERDICT: CONFIRMED
```

### V148

- Finding: A-110 — The validation simulation in fpr_calibration_transfer() (on by default, nsims_check = 200) cannot run with the calibrated configuration: the inner call is built
- Runs: exit 0, 2.4 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
future::plan(future::sequential)
cal <- structure(list(fpr = 0.10, P1 = 0.02, d_eff = 15, L_eff_implied = 5.2, c1 = 1.5, c2 = 1.25,
  call_args = list(nsims = 200L, nsims_verify = 0L, n_min = 60L, otype = "binary", n_workers = 1L, sim_seed = 42L)),
  class = "fpr_calibration")
set.seed(1)
df <- data.frame(id = 1:300, treat = rbinom(300, 1, 0.5), x1 = rnorm(300), y = rbinom(300, 1, 0.3))
fsp <- list(confounders.name = "x1", outcome.name = "y", event.name = "y", treat.name = "treat", id.name = "id",
            outcome_type = "binary", effect_measure = "OR",
            parallel_args = list(plan = "sequential", workers = 1, show_message = FALSE))
a <- fpr_calibration_transfer(cal, delta = 0.05, nsims_check = 2L, quiet = TRUE)   # as in the roxygen example
b <- fpr_calibration_transfer(cal, delta = 0.05, nsims_check = 2L, df.analysis = df, fs_params = fsp, quiet = TRUE)  # as the source comment instructs
print(a[, c("fpr_pred", "fpr_sim", "n_sim", "in_ci")]); print(b[, c("fpr_pred", "fpr_sim", "n_sim", "in_ci")])
dup <- tryCatch(do.call(fpr_calibration, c(list(c1 = 1.55, c2 = 1.30, df.analysis = NULL, fs_params = list(), nsims = 2L),
                                            list(df.analysis = df, fs_params = fsp))), error = conditionMessage)
print(dup)
cat("REPRO_VERDICT:", if (is.na(a$fpr_sim) && a$n_sim == 0 && is.na(b$fpr_sim) && b$n_sim == 0 && is.character(dup)) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
   fpr_pred fpr_sim n_sim in_ci
1 0.4737228      NA     0    NA
   fpr_pred fpr_sim n_sim in_ci
1 0.4737228      NA     0    NA
[1] "formal argument \"df.analysis\" matched by multiple actual arguments"
REPRO_VERDICT: CONFIRMED
```

### V149

- Finding: INV-2 — compute_frontier_cis(seed = ) calls set.seed(seed) on the global RNG and never restores the caller's .Random.seed (invariant 2).
- Runs: exit 0, 2.4 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(99)
df <- data.frame(y = rexp(60), ev = 1L, tr = rep(0:1, 30), x = 0)
fr <- data.table::data.table(Pcons = 0.95, hr = 1.1, N = 60, E = 60, g = "1", m = "1", K = 1, M.1 = "{x <= 0}", M.2 = "")
fs <- structure(list(grp.consistency = list(out_sg = list(result = fr, pareto_frontier = fr)), df.est = df,
                     args_call_all = list(treat.name = "tr", outcome.name = "y", event.name = "ev"),
                     outcome_type = "survival"), class = c("forestsearch", "list"))
set.seed(1); r1 <- try(compute_frontier_cis(fs, n_splits = 5, seed = 123), silent = TRUE); a <- runif(1)
set.seed(2); r2 <- try(compute_frontier_cis(fs, n_splits = 5, seed = 123), silent = TRUE); b <- runif(1)
print(class(r1)); print(c(a = a, b = b))
cat("REPRO_VERDICT:", if (identical(a, b)) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
[1] "data.table" "data.frame"
        a         b 
0.7845753 0.7845753 
REPRO_VERDICT: CONFIRMED
```

### V150

- Finding: A-155 — For an argument list with outcome_type = "binary" and no effect_measure, fs_family_report() resolves the measure to "RD" while forestsearch() resolves it to "OR
- Runs: exit 1, 2.4 s. Verdict: no verdict line (exit 1)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
a <- list(confounders.name = c("a", "b"))
row <- function(r) r[r$stage == "estimability boundary", c("values", "status")]
r0 <- fs_family_report(a, outcome_type = "binary")                                   # no effect_measure
r1 <- fs_family_report(c(a, list(effect_measure = "OR")), outcome_type = "binary")   # what the engine resolves to
print(row(r0)); print(row(r1))
engine_or <- any(grepl('binary = "OR"', deparse(body(forestsearch)), fixed = TRUE))
cat("forestsearch() body defaults binary to OR:", engine_or, "\n")
cat("verdict default :", attr(r0, "verdict"), "\n"); cat("verdict OR      :", attr(r1, "verdict"), "\n")
cat("REPRO_VERDICT:", if (engine_or && !identical(row(r0)$status, row(r1)$status)) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text


Error in data.frame(stage = x$stage, status = x$status, values = ifelse(nchar(x$values) >  : 
  arguments imply differing number of rows: 0, 1
Calls: source ... eval -> print -> print.fs_family_report -> data.frame
Execution halted
```

### V151

- Finding: A-156 — The harm null's k_inter is calibrated on calibrate_k_inter()'s default 5000-row super-population but the DGM is then built with fs_fdr_report()'s n_super (defau
- Runs: exit 0, 5.3 s. Verdict: INCONCLUSIVE

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
k   <- suppressWarnings(calibrate_k_inter(target_hr_harm = 1.0, model = "alt", use_ahr = FALSE))
h5  <- unname(suppressWarnings(setup_gbsg_dgm(model = "alt", k_inter = k, n_super = 5000L))$hr_H_true)
h25 <- unname(suppressWarnings(setup_gbsg_dgm(model = "alt", k_inter = k, n_super = 25000L))$hr_H_true)
print(c(k_inter = k, HR_H_at_calibration_n5000 = h5, HR_H_as_simulated_n25000 = h25, gap = h25 - 1))
# 0.0005 = half a unit of the 3-decimal rounding used in the report's hr_structure table
cat("REPRO_VERDICT:", if (abs(h5 - 1) < 0.0005 && abs(h25 - 1) >= 0.0005) sprintf("CONFIRMED (simulated HR_H = %.4f, not 1.000)", h25) else if (abs(h25 - 1) < 0.0005) "NEGLIGIBLE (gap below table rounding)" else "INCONCLUSIVE", "\n")
```

```text
                  k_inter HR_H_at_calibration_n5000  HR_H_as_simulated_n25000 
              0.567916949               1.006568218               1.001587485 
                      gap 
              0.001587485 
REPRO_VERDICT: INCONCLUSIVE
```

### V152

- Finding: INV-2 — future::plan is replaced with no saved previous plan and no on.exit, then reset to "sequential" rather than to the caller's plan, so the caller's plan is lost o
- Runs: exit 0, 3.6 s. Verdict: CONFIRMED (caller's multisession plan replaced)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
future::plan(future::multisession, workers = 2)          # the caller's plan
was_ms <- inherits(future::plan(), "multisession")
r <- try(.fs_fdr_run_null(gen_df = function(s) NULL, label = "x", fs_params = list(), c1 = 1.25, c2 = 1,
                          dgm_names = list(), conf_use = "a", nsims = 0L, n_workers = 1L, seed_base = 1L,
                          c_confirm = 1, n = 10L, event_min = 5L, quiet = TRUE), silent = TRUE)
cat("call completed without error:", !inherits(r, "try-error"), "\n")
still_ms <- inherits(future::plan(), "multisession")
print(class(future::plan()))
future::plan(future::sequential)
cat("REPRO_VERDICT:", if (was_ms && !still_ms) "CONFIRMED (caller's multisession plan replaced)" else "NOT_REPRODUCED", "\n")
```

```text
call completed without error: TRUE 
[1] "FutureStrategy" "sequential"     "uniprocess"     "future"        
[5] "function"      
REPRO_VERDICT: CONFIRMED (caller's multisession plan replaced)
```

### V153

- Finding: A-112 — With a non-default `roles`, plot_fs_identification_involvement() pairs the i-th requested covariate with the i-th of the fixed anchor / true partner / proxy lab
- Runs: exit 0, 2.4 s. Verdict: CONFIRMED (meno and age both labelled 'er (anchor)')

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
res <- data.frame(detected = 1L,
  sg_def = c("{er <= 0} & {meno == 0}", "{er <= 0} & {age <= 50}", "{er <= 0}", "{meno == 0}"),
  stringsAsFactors = FALSE)
id <- fs_identification_summary(list(`500` = res), anchor = "er", partner = "meno", proxy = "age")
print(id$involvement)                                   # er 0.75, meno 0.50, age 0.25
d <- plot_fs_identification_involvement(id, roles = c("meno", "er"))$data
print(d[, c("covariate", "involvement", "role")])
d1 <- plot_fs_identification_involvement(id, roles = "age")$data
print(d1[, c("covariate", "involvement", "role")])
bad <- identical(as.character(d$role[d$covariate == "meno"]), "er (anchor)") &&
       identical(as.character(d1$role[d1$covariate == "age"]), "er (anchor)")
cat("REPRO_VERDICT:", if (bad) "CONFIRMED (meno and age both labelled 'er (anchor)')" else "NOT_REPRODUCED", "\n")
```

```text
  point covariate involvement
1   500        er        0.75
2   500      meno        0.50
3   500       age        0.25
  covariate involvement                role
1        er        0.75 meno (true partner)
2      meno        0.50         er (anchor)
  covariate involvement        role
3       age        0.25 er (anchor)
REPRO_VERDICT: CONFIRMED (meno and age both labelled 'er (anchor)')
```

### V154

- Finding: INV-2 — fs_oc_predict(), fs_oc_invert() and fs_oc_grid() call set.seed(seed) on the global RNG when seed is non-NULL and never restore the caller's .Random.seed (invari
- Runs: exit 0, 2.3 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
piQ <- 0.34
fam <- structure(list(lab = c("Q", "P", "D"), Pg = c(piQ, 0.45, 0.31), PQg = c(1, 0.28 / 0.45, 1),
  sens_g = c(1, 0.28 / piQ, 0.31 / piQ), spec_g = c(1, 1 - 0.17 / (1 - piQ), 1),
  ovl = matrix(c(piQ, 0.28, 0.31, 0.28, 0.45, 0.28 * 0.31 / piQ, 0.31, 0.28 * 0.31 / piQ, 0.31), 3, 3),
  M = 3L, PQ = piQ, n = 500), class = c("fs_oc_family", "list"))
fam$beta_g <- 26 + 14 * fam$PQg
fam$se_g   <- 13.7 * sqrt(2) * sqrt(piQ / fam$Pg)
chk <- function(f) { set.seed(99); invisible(runif(1)); s0 <- get(".Random.seed", envir = globalenv())
  tryCatch({ f(); identical(s0, get(".Random.seed", envir = globalenv())) }, error = function(e) NA) }
res <- c(
  fs_oc_predict = chk(function() fs_oc_predict(family = fam, n = 500, c1 = 30, c2 = 10, draws = 1000, seed = 1)),
  fs_oc_grid    = chk(function() fs_oc_grid(family = fam, n = 500, c1 = 30, c2 = 10, consistency_method = "resample", draws = 1000, seed = 1)),
  fs_oc_invert  = chk(function() fs_oc_invert(family = fam, n = 500, target = 0.5, solve_for = "c1", c2 = 10, draws = 1000, seed = 1)))
print(res)   # TRUE = caller's RNG state restored; NA = call errored
cat("REPRO_VERDICT:", if (identical(res[["fs_oc_predict"]], FALSE)) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
fs_oc_predict    fs_oc_grid  fs_oc_invert 
        FALSE         FALSE         FALSE 
REPRO_VERDICT: CONFIRMED
```

### V155

- Finding: A-113 — For a non-numeric subgroup variable, any cut specification other than a single value or list(type = "multiple") silently falls back to 'variable == first level'
- Runs: exit 0, 2.2 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
g <- factor(c(1, 2, 3, 2, 3, 1, 2, 3))
r1 <- process_factor_subgroup(g, c("2", "3"), "grade", FALSE)                           # vector of levels
r2 <- process_factor_subgroup(g, list(type = "quantile", value = 0.5), "grade", FALSE)   # list type other than "multiple"
r3 <- process_factor_subgroup(g, list(type = "multiple", values = c("2", "3")), "grade", FALSE)  # supported form
print(c(vector_spec = r1$definition, quantile_spec = r2$definition, multiple_spec = r3$definition))
d <- data.frame(grade = g, treat = rep(0:1, 4))
sg <- define_subgroups(d, d, "grade", list(grade = c("2", "3")), character(0), "alt", FALSE)
print(table(grade = g, flag_harm = sg$flag_harm))    # requested grades 2 and 3
ok <- identical(r1$definition, "grade == 1") && identical(r2$definition, "grade == 1") &&
  identical(as.numeric(sg$flag_harm), as.numeric(g == "1"))
cat("REPRO_VERDICT:", if (ok) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
    vector_spec   quantile_spec   multiple_spec 
   "grade == 1"    "grade == 1" "grade in 2, 3" 
     flag_harm
grade 0 1
    1 0 2
    2 3 0
    3 3 0
REPRO_VERDICT: CONFIRMED
```

### V156

- Finding: INV-2 — Twenty-nine set.seed() calls in functions off the forestsearch() path (DGM builders, simulators, OC, Guo-He, cross-validation and calibration helpers) reseed th
- Runs: exit 0, 2.7 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
df <- survival::gbsg
build <- function() generate_aft_dgm_flex(data = df, outcome_var = "rfstime", event_var = "status",
  treatment_var = "hormon", continuous_vars = c("age", "size", "nodes", "pgr", "er"), factor_vars = "meno",
  model = "null", n_super = 500, verbose = FALSE)
set.seed(1); invisible(build()); a <- runif(1)
set.seed(2); invisible(build()); b <- runif(1)
set.seed(1); a0 <- runif(1); set.seed(2); b0 <- runif(1)
cat("next caller draw after the build, caller seeds 1 and 2:", a, b, "\n")
cat("same draws without the build:", a0, b0, "\n")
cat("REPRO_VERDICT:", if (identical(a, b) && !identical(a0, b0)) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
next caller draw after the build, caller seeds 1 and 2: 0.5400066 0.5400066 
same draws without the build: 0.2655087 0.1848823 
REPRO_VERDICT: CONFIRMED
```

### V157

- Finding: A-116 — For effect_measure = "log_RR", overdispersion = "quasi" or "negbin" fits a plain Poisson model and also turns off the sandwich variance, so the profile's SEs an
- Runs: exit 0, 2.4 s. Verdict: CONFIRMED (quasi/negbin ignored for log_RR, model-based Poisson SEs)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(2); n <- 600
df <- data.frame(treat = rbinom(n, 1, 0.5), bm = rnorm(n))
df$y <- rbinom(n, 1, plogis(-0.5 + 0.3 * df$treat))
a  <- glm_effect_profile(df, "y", effect_measure = "log_RR", overdispersion = "none")
b  <- glm_effect_profile(df, "y", effect_measure = "log_RR", overdispersion = "quasi")
nb <- glm_effect_profile(df, "y", effect_measure = "log_RR", overdispersion = "negbin")
print(c(family_used = b$family_used, overdispersion = b$overdispersion, dispersion = b$dispersion))
print(c(quasi_fit_family = family(b$model_fit)$family, negbin_fit_class = class(nb$model_fit)[1]))
print(summary(b$se / a$se))  # ratio > 1 when 'sandwich' is installed: the robust SE was dropped
ok <- identical(family(b$model_fit)$family, "poisson") && !inherits(nb$model_fit, "negbin") && identical(b$se, nb$se) && identical(b$overdispersion, "quasi") && is.na(b$dispersion)
cat("REPRO_VERDICT:", if (ok) "CONFIRMED (quasi/negbin ignored for log_RR, model-based Poisson SEs)" else "NOT_REPRODUCED", "\n")
```

```text
                             family_used 
"poisson(log) [modified Poisson for RR]" 
                          overdispersion 
                                 "quasi" 
                              dispersion 
                                      NA 
quasi_fit_family negbin_fit_class 
       "poisson"            "glm" 
   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
  1.247   1.290   1.313   1.333   1.379   1.448 
REPRO_VERDICT: CONFIRMED (quasi/negbin ignored for log_RR, model-based Poisson SEs)
```

### V158

- Finding: A-118 — For a fit made with subgroup_method = "grf", extract_grf_args() omits grf_selection, frontier_rule, effect_neighborhood and selection_rule and reports the store
- Runs: exit 0, 2.2 s. Verdict: CONFIRMED (extracted args are the tree defaults, not the frontier call)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
# extract_grf_args() reads only fs$outcome_type and fs$args_call_all, so a stub with the stored formals of a subgroup_method = 'grf' run suffices
fs <- list(outcome_type = "survival", args_call_all = list(df.analysis = data.frame(), confounders.name = "x",
  outcome.name = "tte", event.name = "event", id.name = "id", treat.name = "treat", frac.tau = 0.8, n.min = 60,
  dmin.grf = 0, is.RCT = TRUE, grf_depth = 2, seedit = 1, return_selected_cuts_only = FALSE,
  subgroup_method = "grf", grf_selection = c("frontier", "tree"), sg_focus = "hr", selection_rule = "pareto",
  effect_neighborhood = 0.20))
a <- extract_grf_args(fs)
print(a[c("grf_selection", "frontier_rule", "selection_rule", "effect_neighborhood", "return_selected_cuts_only")])
# the run itself used: frontier / eff (sg_focus = 'hr') / pareto / 0.20 / TRUE
ok <- identical(a$grf_selection, "tree") && identical(a$frontier_rule, "effMaxSG") && identical(a$selection_rule, "neighborhood") && identical(a$effect_neighborhood, 0.10) && identical(a$return_selected_cuts_only, FALSE)
cat("REPRO_VERDICT:", if (ok) "CONFIRMED (extracted args are the tree defaults, not the frontier call)" else "NOT_REPRODUCED", "\n")
```

```text
$grf_selection
[1] "tree"

$frontier_rule
[1] "effMaxSG"

$selection_rule
[1] "neighborhood"

$effect_neighborhood
[1] 0.1

$return_selected_cuts_only
[1] FALSE

REPRO_VERDICT: CONFIRMED (extracted args are the tree defaults, not the frontier call)
```

### V159

- Finding: INV-2 — guohe_algorithm3() and guohe_adaptive_r() call set.seed(seed) on the caller's global RNG when a seed is supplied and never restore the previous state (invariant
- Runs: exit 0, 2.4 s. Verdict: REPRODUCED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 160
df <- data.frame(treat = rbinom(n, 1, 0.5), y = rbinom(n, 1, 0.3), g1 = rbinom(n, 1, 0.5), g2 = rbinom(n, 1, 0.5))
g <- function() invisible(tryCatch(suppressWarnings(guohe_algorithm3(df, outcome = "binary", treatment = "treat",
  candidates = c("g1", "g2"), y = "y", orient = +1, B = 20, seed = 2)), error = function(e) e))
set.seed(99);  g(); s99  <- .Random.seed; a <- runif(1)
set.seed(100); g(); s100 <- .Random.seed; b <- runif(1)
cat("caller RNG state identical after calls made under seeds 99 and 100:", identical(s99, s100), "\n")
cat("next draw identical:", identical(a, b), "\n")
cat("REPRO_VERDICT:", if (identical(s99, s100) && identical(a, b)) "REPRODUCED" else "NOT_REPRODUCED", "\n")
```

```text
caller RNG state identical after calls made under seeds 99 and 100: TRUE 
next draw identical: TRUE 
REPRO_VERDICT: REPRODUCED
```

### V160

- Finding: A-121 — mrct_case_study() restores the caller's RNG kind before running the identifier while mrct_region_sims() runs the identifier under the workers' L'Ecuyer-CMRG kin
- Runs: exit 0, 4.9 s. Verdict: CONFIRMED (same seed, different split outcomes under the two kinds)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
# The identifier's split step under the plain loop: set.seed(seed) then sample() in run_single_consistency_split().
# Same seed, the two generator kinds (case study = caller's Mersenne-Twister; simulation worker = L'Ecuyer-CMRG).
old_kind <- RNGkind()
set.seed(1)
df <- data.table::data.table(Y = rexp(200), Event = rbinom(200, 1, 0.7), Treat = rep(0:1, 100))
splits <- function(kind) {
  RNGkind(kind); set.seed(8316951)
  vapply(seq_len(400), function(i) run_single_consistency_split(df, N.x = 200, hr.consistency = 1.0), numeric(1))
}
a <- splits("Mersenne-Twister"); b <- splits("L'Ecuyer-CMRG")
do.call(RNGkind, as.list(old_kind))
cat("Pcons, Mersenne-Twister:", mean(a, na.rm = TRUE), " Pcons, L'Ecuyer-CMRG:", mean(b, na.rm = TRUE), "\n")
cat("REPRO_VERDICT:", if (!identical(a, b)) "CONFIRMED (same seed, different split outcomes under the two kinds)" else "NOT_REPRODUCED", "\n")
```

```text
Pcons, Mersenne-Twister: 0.01  Pcons, L'Ecuyer-CMRG: 0.0125 
REPRO_VERDICT: CONFIRMED (same seed, different split outcomes under the two kinds)
```

### V161

- Finding: INV-2 — The MRCT DGM builders (expand_covariates, the make_region_model draw closure, inject_mrct_structure, mrct_truth_by_region) call set.seed() on the session genera
- Runs: exit 0, 2.2 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
d <- data.frame(a = rnorm(50), b = rbinom(50, 1, 0.5))
set.seed(42); invisible(runif(1)); before <- .Random.seed
x <- expand_covariates(d, c("a", "b"), 20, method = "bootstrap")
after <- .Random.seed
set.seed(20260903L); invisible(sample.int(50, 20, replace = TRUE)); fixed <- .Random.seed
cat("caller's stream preserved:", identical(before, after), "\n")
cat("stream left at the function's default-seed state:", identical(after, fixed), "\n")
cat("REPRO_VERDICT:", if (!identical(before, after) && identical(after, fixed)) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
caller's stream preserved: FALSE 
stream left at the function's default-seed state: TRUE 
REPRO_VERDICT: CONFIRMED
```

### V162

- Finding: A-122 — In the planted-region mode region$loghr is not validated, so region = list(loghr = NULL) deletes the default 0 and, together with x3 or region_treat_loghr, the 
- Runs: exit 0, 2.9 s. Verdict: CONFIRMED (x3 log-HR recycled into the region coefficient)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
c1 <- survival::cgd; c1 <- c1[c1$enum == 1, ]
d <- data.frame(tte = c1$tstop - c1$tstart, event = c1$status, treat = as.numeric(c1$treat == "rIFN-g"),
  age = c1$age, height = c1$height, weight = c1$weight, female = as.numeric(c1$sex == "female"),
  autosom = as.numeric(c1$inherit == "autosomal"), steroids = c1$steroids, propylac = c1$propylac)
mk <- function(region) suppressWarnings(inject_mrct_structure(seed_data = d, outcome_var = "tte",
  event_var = "event", treatment_var = "treat", continuous_vars = c("age", "height", "weight"),
  factor_vars = c("female", "autosom", "steroids", "propylac"), x_pred = "age",
  spline_spec = list(knot = 12, zeta = 25, log_hrs = log(c(0.35, 0.80, 1.30))),
  region = region, x3 = list(vars = "weight", cuts = list(weight = 40), loghr = log(2)),
  n_super = 2000, expand = "bootstrap", seed = 1))
base <- list(prevalence = 0.20, or_pred = 10, x1_vars = "female", or_x1 = 3)
g0 <- mk(base); g2 <- mk(c(base, list(loghr = NULL)))
v0 <- unname(g0$model_params$b0["z_region"]); v2 <- unname(g2$model_params$b0["z_region"])
print(c(default = v0, loghr_NULL = v2, x3_loghr = log(2), reported = g2$mrct$region$prognostic_loghr))
cat("REPRO_VERDICT:", if (isTRUE(all.equal(v0, 0)) && isTRUE(all.equal(v2, log(2)))) "CONFIRMED (x3 log-HR recycled into the region coefficient)" else "NOT_REPRODUCED", "\n")
```

```text
   default loghr_NULL   x3_loghr   reported 
 0.0000000  0.6931472  0.6931472  0.6931472 
REPRO_VERDICT: CONFIRMED (x3 log-HR recycled into the region coefficient)
```

### V163

- Finding: A-123 — The region-dependent effect (region_loghr_shift / region_harm_loghr_shift) is added to dgm$df_super only, while dgm$df_source, the frame simulate_from_dgm(basel
- Runs: exit 0, 2.9 s. Verdict: CONFIRMED (shift in df_super, absent from df_source)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
c1 <- survival::cgd; c1 <- c1[c1$enum == 1, ]
d <- data.frame(tte = c1$tstop - c1$tstart, event = c1$status, treat = as.numeric(c1$treat == "rIFN-g"),
  age = c1$age, height = c1$height, weight = c1$weight, female = as.numeric(c1$sex == "female"),
  autosom = as.numeric(c1$inherit == "autosomal"), steroids = c1$steroids, propylac = c1$propylac)
mk <- function(...) suppressWarnings(inject_mrct_structure(seed_data = d, outcome_var = "tte",
  event_var = "event", treatment_var = "treat", continuous_vars = c("age", "height", "weight"),
  factor_vars = c("female", "autosom", "steroids", "propylac"), x_pred = "age",
  spline_spec = list(knot = 12, zeta = 25, log_hrs = log(c(0.35, 0.80, 1.30))),
  region = list(prevalence = 0.20, or_pred = 10, x1_vars = "female", or_x1 = 3, loghr = 0),
  n_super = 2000, expand = "bootstrap", seed = 1, ...))
g0 <- mk(); g1 <- mk(region_loghr_shift = log(2))
super_shifted <- max(abs(g1$df_super$lin_pred_1 - g0$df_super$lin_pred_1)) > 0
source_same <- !is.null(g1$df_source) && identical(g0$df_source$lin_pred_1, g1$df_source$lin_pred_1) &&
  identical(g0$df_source$loghr_po, g1$df_source$loghr_po)
print(g1$mrct$region_effect); print(c(super_shifted = super_shifted, source_unshifted = source_same))
cat("REPRO_VERDICT:", if (super_shifted && source_same) "CONFIRMED (shift in df_super, absent from df_source)" else "NOT_REPRODUCED", "\n")
```

```text
$region_loghr_shift
[1] 0.6931472

$region_harm_loghr_shift
[1] 0

   super_shifted source_unshifted 
            TRUE             TRUE 
REPRO_VERDICT: CONFIRMED (shift in df_super, absent from df_source)
```

### V164

- Finding: A-124 — The pattern that removes the spline helper columns, "_treat$\|_k$", is not tied to the spline variable, so any real covariate whose z_ name ends in _treat or _k 
- Runs: exit 0, 2.2 s. Verdict: CONFIRMED (real covariates dropped by the name pattern)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
dgm_f <- list(
  df_super = data.frame(z_age = 1, z_prior_treat = 0, z_risk_k = 2, z_x = 1,
                        z_x_treat = 0, z_x_k = 0, z_x_k_treat = 0, z_region = 0),
  mrct = list(region = list(z_name = "z_region"), classes = list(X2 = "x")))
out <- analysis_covariates(dgm_f)
print(out)   # documented: z_age, z_prior_treat, z_risk_k, z_x
cat("REPRO_VERDICT:", if (!any(c("z_prior_treat", "z_risk_k") %in% out) && all(c("z_age", "z_x") %in% out)) "CONFIRMED (real covariates dropped by the name pattern)" else "NOT_REPRODUCED", "\n")
```

```text
[1] "z_age" "z_x"  
REPRO_VERDICT: CONFIRMED (real covariates dropped by the name pattern)
```

### V165

- Finding: INV-2 — mrct_region_sims(), .fs_fdr_run_null() and fpr_calibration() replace the caller's future plan without capturing it: on normal exit the session is left on sequen
- Runs: exit 0, 5.1 s. Verdict: CONFIRMED (caller's multisession plan replaced and not restored)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
# A stub DGM passes the input checks and reaches the plan set-up (mrct_simulation.R:357-379);
# whatever happens afterwards, the caller's plan must still be in force when the function returns or stops.
dgm_stub <- structure(list(df_super = data.frame(z_regA = c(0, 1))), class = "aft_dgm_flex")
future::plan(future::multisession, workers = 2)          # the caller's plan
before <- class(future::plan())[1]
r <- try(suppressWarnings(mrct_region_sims(dgm_stub, n_sims = 1,
  parallel_args = list(plan = "sequential", show_message = FALSE), seed = 1)), silent = TRUE)
after <- class(future::plan())[1]
cat("plan before:", before, " plan after:", after, " call errored:", inherits(r, "try-error"), "\n")
cat("REPRO_VERDICT:", if (inherits(future::plan(), "multisession")) "NOT_REPRODUCED" else "CONFIRMED (caller's multisession plan replaced and not restored)", "\n")
future::plan(future::sequential)
```

```text
plan before: FutureStrategy  plan after: FutureStrategy  call errored: TRUE 
REPRO_VERDICT: CONFIRMED (caller's multisession plan replaced and not restored)
```

### V166

- Finding: A-126 — A user-supplied sim_args$entry_var is silently removed whenever dgm$df_super has no column named 'entrytime', because the override list carries entry_var = NULL
- Runs: exit 0, 3.1 s. Verdict: CONFIRMED (supplied entry_var was dropped)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(11); n <- 400
d <- data.frame(tte = rexp(n, 0.05), event = rbinom(n, 1, 0.75), treat = rep(0:1, n / 2),
                age = round(rnorm(n, 62, 8)), bm = round(runif(n, 0, 12), 1),
                male = rbinom(n, 1, 0.6), histology = rbinom(n, 1, 0.4),
                prior_treat = rbinom(n, 1, 0.3), regA = rbinom(n, 1, 0.3))
dgm <- create_dgm_for_mrct(d, model_type = "alt", verbose = FALSE)
dgm$df_super$enroll <- 0     # every patient enters at time 0; no 'entrytime' column
cs <- mrct_case_study(dgm, seed = 1, k = 1, n = 500, region_var = "z_regA",
  fs_args = list(confounders.name = c("z_age", "z_bm"), conf_force = character(0), fs.splits = 100),
  sim_args = list(rand_ratio = 1, draw_treatment = TRUE, entry_var = "enroll"),
  analysis_time = 60, mr_inference = FALSE)
passed <- names(cs$call_args$sim_args)      # the arguments actually given to simulate_from_dgm()
print(passed)
print(utils::modifyList(list(entry_var = "enroll"), list(entry_var = NULL)))   # the mechanism: list()
cat("REPRO_VERDICT:", if (!"entry_var" %in% passed) "CONFIRMED (supplied entry_var was dropped)" else "NOT REPRODUCED", "\n")
```

```text
[1] "rand_ratio"     "draw_treatment" "n"              "analysis_time" 
[5] "cens_adjust"    "seed"          
named list()
REPRO_VERDICT: CONFIRMED (supplied entry_var was dropped)
```

### V167

- Finding: A-127 — An error thrown by forestsearch() inside a replicate is caught and turned into NULL, and the replicate is returned as an ordinary 'no subgroup found' row (any_f
- Runs: exit 0, 4.7 s. Verdict: CONFIRMED (every replicate failed in forestsearch() and was returned as 'none' with no error or warning)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(11); n <- 400
d <- data.frame(tte = rexp(n, 0.05), event = rbinom(n, 1, 0.75), treat = rep(0:1, n / 2),
                age = round(rnorm(n, 62, 8)), bm = round(runif(n, 0, 12), 1),
                male = rbinom(n, 1, 0.6), histology = rbinom(n, 1, 0.4),
                prior_treat = rbinom(n, 1, 0.3), regA = rbinom(n, 1, 0.3))
dgm <- create_dgm_for_mrct(d, model_type = "alt", verbose = FALSE)
w <- character(0)
res <- withCallingHandlers(
  mrct_region_sims(dgm, n_sims = 3, n_sample = 500, region_var = "z_regA",
    confounders.name = c("z_age", "z_bm"), conf_force = character(0),
    fs_args = list(use_grff = FALSE, fs.splits = 100),     # misspelled argument: forestsearch() errors
    parallel_args = list(plan = "sequential", show_message = FALSE), seed = 1),
  warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
r <- as.data.frame(res)
print(r[, c("sim", "any_found", "sg_found", "hr_test", "hr_sg")]); print(w)
ok <- nrow(r) == 3 && all(r$any_found == 0) && all(r$sg_found == "none") && all(is.na(r$hr_sg[!is.na(r$hr_test)])) && !any(grepl("unused argument|ForestSearch|forestsearch", w))
cat("REPRO_VERDICT:", if (ok) "CONFIRMED (every replicate failed in forestsearch() and was returned as 'none' with no error or warning)" else "NOT REPRODUCED", "\n")
```

```text
  sim any_found sg_found   hr_test hr_sg
1   1         0     none 0.7322864    NA
2   2         0     none 0.9897811    NA
3   3         0     none 0.8366032    NA
character(0)
REPRO_VERDICT: CONFIRMED (every replicate failed in forestsearch() and was returned as 'none' with no error or warning)
```

### V168

- Finding: A-128 — mrct_trial_analysis() does not check the declared rule's covariates for missing values, so a region patient with NA in a variable of the declared rule gets in_H
- Runs: exit 0, 2.7 s. Verdict: CONFIRMED (NA membership returned; patients silently dropped from region_Ghat)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 800
d <- data.frame(id = seq_len(n), treat = rbinom(n, 1, 0.5), x1 = rnorm(n), x2 = rnorm(n))
d$region <- rbinom(n, 1, plogis(-1.5 - 1.5 * d$x1))
d$time <- rexp(n, 0.1 * ifelse(d$treat == 1, ifelse(d$x1 < -0.6, 3, 0.45), 1))
d$event <- as.integer(d$time < 30); d$time <- pmin(d$time, 30)
miss <- which(d$region == 1)[1:10]; d$x1[miss] <- NA; d$x2[miss] <- NA   # 10 region patients, covariates missing
fit <- tryCatch(mrct_trial_analysis(d, "time", "event", "treat", "id", region.name = "region",
  fs_args = list(sg_focus = "hrMaxSG", confounders.name = c("x1", "x2"), use_grf = FALSE, use_lasso = FALSE, fs.splits = 200),
  mr_inference = FALSE, proceed_without_trigger = TRUE), error = function(e) conditionMessage(e))
if (is.character(fit)) { print(fit); cat("REPRO_VERDICT: INCONCLUSIVE (the call stopped)\n") } else {
  print(fit$stage); print(fit$Hhat); print(table(in_Hhat = fit$in_Hhat, useNA = "always"))
  nG <- fit$sets$n[fit$sets$set == "region_Ghat"]
  expect <- sum(d$region == 1) - sum(fit$in_Hhat & d$region == 1, na.rm = TRUE)
  cat("region_Ghat n:", nG, " region minus H-hat members:", expect, "\n")
  cat("REPRO_VERDICT:", if (fit$stage != "complete") "INCONCLUSIVE (no subgroup declared)" else if (anyNA(fit$in_Hhat) && nG < expect) "CONFIRMED (NA membership returned; patients silently dropped from region_Ghat)" else "NOT REPRODUCED", "\n")
}
```

```text
[1] "complete"
[1] "{x1 <= -0.4}"
in_Hhat
FALSE  TRUE  <NA> 
  500   290    10 
region_Ghat n: 56  region minus H-hat members: 66 
REPRO_VERDICT: CONFIRMED (NA membership returned; patients silently dropped from region_Ghat)
```

### V169

- Finding: A-129 — compute_frontier_cis() accepts ci_level but returns no record of it, and every consumer labels the intervals '95%', so a table or explanation built from a ci_ta
- Runs: exit 0, 2.6 s. Verdict: CONFIRMED (50% bounds shown under 'Naive 95% CI'; level not carried on the table)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
library(data.table)
set.seed(99)
df <- data.frame(y = rexp(60), ev = 1L, tr = rep(0:1, 30), x = 0)
fr <- data.table(Pcons = 0.95, hr = 1.1, N = 60, E = 60, g = "1", m = "1", K = 1, M.1 = "{x <= 0}", M.2 = "")
fs <- structure(list(grp.consistency = list(out_sg = list(result = fr, pareto_frontier = fr)), df.est = df,
                     args_call_all = list(treat.name = "tr", outcome.name = "y", event.name = "ev", sg_focus = "hr"),
                     outcome_type = "survival"), class = c("forestsearch", "list"))
ci50 <- compute_frontier_cis(fs, n_splits = 5, ci_level = 0.50, seed = 1)
ci95 <- compute_frontier_cis(fs, n_splits = 5, ci_level = 0.95, seed = 1)
print(rbind(level50 = unlist(ci50[, .(naive_lcl, naive_ucl)]), level95 = unlist(ci95[, .(naive_lcl, naive_ucl)])))
carried <- "ci_level" %in% c(names(ci50), names(attributes(ci50)))
tab <- tryCatch(as.data.frame(pareto_frontier_table(fs, ci_table = ci50, format = "data.table")), error = function(e) conditionMessage(e))
print(tab)
narrower <- isTRUE((ci50$naive_ucl - ci50$naive_lcl) < (ci95$naive_ucl - ci95$naive_lcl))
labelled95 <- is.data.frame(tab) && "Naive 95% CI" %in% names(tab)
cat("REPRO_VERDICT:", if (narrower && !carried && labelled95) "CONFIRMED (50% bounds shown under 'Naive 95% CI'; level not carried on the table)" else "NOT REPRODUCED", "\n")
```

```text

Attaching package: ‘data.table’

The following object is masked from ‘package:base’:

    %notin%

        naive_lcl naive_ucl
level50 1.0651962  1.499862
level95 0.7687841  2.078149
  is_selected      M.1 M.2  N  E  HR Pcons K Naive 95% CI Split ~ 95% CI
1        TRUE {x <= 0}     60 60 1.1  0.95 1 (1.07, 1.50)   (1.04, 1.54)
      FSBC ~ 95% CI
1 1.31 (1.07, 1.60)
REPRO_VERDICT: CONFIRMED (50% bounds shown under 'Naive 95% CI'; level not carried on the table)
```

### V170

- Finding: INV-2 — Every call leaves the caller's global RNG at the state following a fixed seed, because weightedsurv::plotKM.band_subgroups() is invoked without saving and resto
- Runs: exit 0, 2.4 s. Verdict: CONFIRMED (the caller's RNG state is overwritten with a fixed seed)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 200
d <- data.frame(tte = rexp(n), event = rbinom(n, 1, 0.7), treat = rbinom(n, 1, 0.5), g = rbinom(n, 1, 0.5))
grDevices::pdf(NULL)
run <- function() try(plot_km_band_forestsearch(d, sg_cols = "g", sg_labels = "g == 1", sg_colors = "blue", draws_band = 20), silent = TRUE)
set.seed(42); r1 <- suppressWarnings(run()); a <- runif(1)
set.seed(99); r2 <- suppressWarnings(run()); b <- runif(1)
set.seed(42); a0 <- runif(1)              # what the user's stream would give without the call
grDevices::dev.off()
cat("plot errored:", inherits(r1, "try-error"), "\n"); print(c(after42 = a, after99 = b, plain42 = a0))
cat("REPRO_VERDICT:", if (identical(a, b) && !identical(a, a0)) "CONFIRMED (the caller's RNG state is overwritten with a fixed seed)" else "NOT REPRODUCED", "\n")
```

```text
plot errored: FALSE 
  after42   after99   plain42 
0.1836983 0.1836983 0.9148060 
REPRO_VERDICT: CONFIRMED (the caller's RNG state is overwritten with a fixed seed)
```

### V171

- Finding: A-130 — sg_labels are handed to weightedsurv as row-subsetting expressions while the validated 0/1 indicator columns are never used, so the default bare-column label, a
- Runs: exit 0, 2.3 s. Verdict: CONFIRMED (auto label '{meno}' indexes row 1 repeatedly; display label is not a parseable expression)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
d <- data.frame(id = 1:6, meno = c(0, 1, 1, 0, 1, 0), age = c(40, 50, 60, 70, 80, 45))
# (ii) one-factor subgroup on a binary indicator, consistency-path label shape
fs <- list(grp.consistency = list(out_sg = list(sg.harm_label = "{meno}")), sg.harm = "{meno}")
lab <- generate_readable_sg_labels(fs); print(lab)
idx <- with(d, eval(parse(text = lab[1])))   # what weightedsurv::plotKM.band_subgroups does with a label
print(class(idx)); print(d[idx, ])           # numeric index: row 1 three times, not the meno == 1 rows
# (iii) documented display label for a reference subgroup
r <- create_reference_subgroup_columns(d, list(young = list(subset_expr = "age < 65", label = "Young patients")))
print(r$labels); p <- tryCatch(parse(text = r$labels), error = function(e) "PARSE_ERROR"); print(p)
ok <- is.numeric(idx) && nrow(d[idx, ]) == 3 && all(d[idx, "id"] == 1) && identical(r$labels, "Young patients") && identical(p, "PARSE_ERROR")
cat("REPRO_VERDICT:", if (ok) "CONFIRMED (auto label '{meno}' indexes row 1 repeatedly; display label is not a parseable expression)" else "NOT REPRODUCED", "\n")
```

```text
[1] "{meno}"    "!({meno})"
[1] "numeric"
    id meno age
1    1    0  40
1.1  1    0  40
1.2  1    0  40
[1] "Young patients"
[1] "PARSE_ERROR"
REPRO_VERDICT: CONFIRMED (auto label '{meno}' indexes row 1 repeatedly; display label is not a parseable expression)
```

### V172

- Finding: A-131 — When a fit has no grp.consistency$out_sg$sg.harm_label (subgroup_method = "dina" or "grf") a multi-factor sg.harm is pasted into one string before generate_comp
- Runs: exit 0, 2.3 s. Verdict: REPRODUCED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
fs <- list(sg.harm = c("{x1 <= 3.5}", "{x2 >= 10}"))   # shape of a DINA depth-2 sg.harm
lab <- generate_readable_sg_labels(fs)
print(lab)
d <- expand.grid(x1 = c(1, 5), x2 = c(5, 15))
h  <- with(d, eval(parse(text = lab[1])))
hc <- with(d, eval(parse(text = lab[2])))
print(cbind(d, harm = h, plotted_complement = hc, true_complement = !h))
cat("REPRO_VERDICT:", if (!identical(hc, !h)) "REPRODUCED" else "NOT_REPRODUCED", "\n")
```

```text
[1] "{x1 <= 3.5} & {x2 >= 10}" "{x1 > 3.5} & {x2 >= 10}" 
  x1 x2  harm plotted_complement true_complement
1  1  5 FALSE              FALSE            TRUE
2  5  5 FALSE              FALSE            TRUE
3  1 15  TRUE              FALSE           FALSE
4  5 15 FALSE               TRUE            TRUE
REPRO_VERDICT: REPRODUCED
```

### V173

- Finding: A-132 — est.scale = "1/hr" only reverses the arm labels and swaps the two subgroup names: the reversed treatment column df$treat_plot is created and never used, so ever
- Runs: exit 0, 2.5 s. Verdict: REPRODUCED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); n <- 300
d <- data.frame(Y = rexp(n), Event = rbinom(n, 1, 0.7), Treat = rbinom(n, 1, 0.5),
                treat.recommend = rbinom(n, 1, 0.7))
fs <- list(df.est = d, sg.harm = "{x <= 1}")
grDevices::pdf(NULL)
a <- plot_sg_results(fs, plot_type = "km", est.scale = "hr")
b <- plot_sg_results(fs, plot_type = "km", est.scale = "1/hr")
grDevices::dev.off()
print(a$hr_estimates); print(b$hr_estimates)
cat("REPRO_VERDICT:", if (isTRUE(all.equal(a$hr_estimates, b$hr_estimates)) && !"treat_plot" %in% all.names(body(plot_sg_km))) "REPRODUCED" else "NOT_REPRODUCED", "\n")
```

```text
  Subgroup   n events       HR  HR_lower HR_upper
1      ITT 300    205 1.496714 1.1313013 1.980157
2        H  86     64 1.259944 0.7671687 2.069244
3       Hc 214    141 1.561116 1.1109672 2.193659
  Subgroup   n events       HR  HR_lower HR_upper
1      ITT 300    205 1.496714 1.1313013 1.980157
2        H  86     64 1.259944 0.7671687 2.069244
3       Hc 214    141 1.561116 1.1109672 2.193659
REPRO_VERDICT: REPRODUCED
```

### V174

- Finding: A-134 — conf.level changes the column header and the unadjusted rows but not the bias-corrected rows, which copy the bootstrap object's fixed 95% limits, so with conf.l
- Runs: exit 0, 2.8 s. Verdict: REPRODUCED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
if (requireNamespace("forestploter", quietly = TRUE)) {
  set.seed(1); n <- 300
  d <- data.frame(id = 1:n, Y = rexp(n), Event = rbinom(n, 1, 0.7), Treat = rbinom(n, 1, 0.5),
                  treat.recommend = rbinom(n, 1, 0.6))
  fs <- list(df.est = d, sg.harm = "{x <= 1}")
  bc <- list(H_estimates  = data.frame(H2 = 1.50, H2_lower = 1.00, H2_upper = 2.25),
             Hc_estimates = data.frame(H2 = 0.80, H2_lower = 0.60, H2_upper = 1.07))
  run <- function(cl) tryCatch(plot_subgroup_results_forestplot(list(fs.est = fs, fs_bc = bc), d,
      outcome.name = "Y", event.name = "Event", treat.name = "Treat", conf.level = cl,
      xlim = c(0.1, 10), ticks_at = c(0.1, 1, 10))$data, error = function(e) { message("ERR: ", conditionMessage(e)); NULL })
  a <- run(0.95); b <- run(0.50)
  ok <- !is.null(a) && !is.null(b)
  if (ok) {
    ib <- grepl("Bias-corrected", a$Subgroup); iu <- grepl("Full-analysis", a$Subgroup)
    print(names(a)); print(names(b)); print(a[ib | iu, c("Subgroup", "low", "hi")]); print(b[ib | iu, c("Subgroup", "low", "hi")])
    ok <- any(ib) && isTRUE(all.equal(a$low[ib], b$low[ib])) && isTRUE(all.equal(a$hi[ib], b$hi[ib])) &&
          !isTRUE(all.equal(a$low[iu], b$low[iu])) && !identical(names(a), names(b))
  }
  cat("REPRO_VERDICT:", if (ok) "REPRODUCED" else "NOT_REPRODUCED", "\n")
} else cat("REPRO_VERDICT: SKIPPED (forestploter not installed)\n")
```

```text
[1] "Subgroup"     "Experimental" "Control"      "est"          "low"         
[6] "hi"           "se"           " "            "HR (95% CI)" 
[1] "Subgroup"     "Experimental" "Control"      "est"          "low"         
[6] "hi"           "se"           " "            "HR (50% CI)" 
          Subgroup       low       hi
3    Full-analysis 1.2179095 2.484251
4   Bias-corrected 0.6000000 1.070000
6    Full-analysis 0.7197888 1.730049
7   Bias-corrected 1.0000000 2.250000
          Subgroup       low       hi
3    Full-analysis 1.5386389 1.966409
4   Bias-corrected 0.6000000 1.070000
6    Full-analysis 0.9596205 1.297669
7   Bias-corrected 1.0000000 2.250000
REPRO_VERDICT: REPRODUCED
```

### V175

- Finding: A-135 — The documented argument cox_formula ('Optional Cox formula for unadjusted ITT') is threaded through four internal functions and never read, so hr.itt is always 
- Runs: exit 0, 2.4 s. Verdict: REPRODUCED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
uses <- function(f, nm) nm %in% all.names(body(f))
r <- c(fs_extractor_reads_cox_formula  = uses(.extract_fs_estimates_gen,  "cox_formula"),
       grf_extractor_reads_cox_formula = uses(.extract_grf_estimates_gen, "cox_formula"),
       fs_extractor_reads_cox_formula_adj = uses(.extract_fs_estimates_gen, "cox_formula_adj"),
       grf_extractor_reads_debug = uses(.extract_grf_estimates_gen, "debug"))
print(r)   # the two leaf functions that compute hr.itt never reference cox_formula
cat("REPRO_VERDICT:", if (!r[[1]] && !r[[2]] && r[[3]] && !r[[4]]) "REPRODUCED" else "NOT_REPRODUCED", "\n")
```

```text
    fs_extractor_reads_cox_formula    grf_extractor_reads_cox_formula 
                             FALSE                              FALSE 
fs_extractor_reads_cox_formula_adj          grf_extractor_reads_debug 
                              TRUE                              FALSE 
REPRO_VERDICT: REPRODUCED
```

### V176

- Finding: INV-2 — Invariant 2: run_simulation_analysis() reseeds the global RNG on every call (through the simulator's seed argument, and directly when n_add_noise > 0) and never
- Runs: exit 0, 2.3 s. Verdict: REPRODUCED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
dgm <- structure(list(
  df_super = data.frame(v1 = rnorm(500), mu0 = 0, mu1 = 0.2, flag_harm = 0L),
  outcome_type = "continuous", model_params = list(sigma = 1)),
  class = c("glm_dgm", "list"))
set.seed(42); s0 <- get(".Random.seed", envir = globalenv())
r <- try(suppressWarnings(run_simulation_analysis(
  sim_id = 1, dgm = dgm, n_sample = 50, confounders_base = "v1", n_add_noise = 1L,
  run_fs = FALSE, run_fs_grf = FALSE, run_grf = FALSE)), silent = TRUE)   # no analysis is run
if (inherits(r, "try-error")) message("call ended in error (state compared anyway): ", conditionMessage(attr(r, "condition")))
changed <- !identical(s0, get(".Random.seed", envir = globalenv()))
cat("REPRO_VERDICT:", if (changed) "REPRODUCED" else "NOT_REPRODUCED", "\n")
```

```text
REPRO_VERDICT: REPRODUCED
```

### V177

- Finding: A-136 — On the legacy boolean path (methods = NULL, the default) an explicit fs_params = list(n.min = NULL) is deleted by utils::modifyList() before it reaches forestse
- Runs: exit 0, 2.3 s. Verdict: REPRODUCED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
fs_defaults <- default_sim_params()
user <- list(n.min = NULL)
legacy <- utils::modifyList(utils::modifyList(fs_defaults, list(use_lasso = TRUE, use_grf = FALSE)), user)   # lines 587-590
meth   <- .modify_keep_null(.modify_keep_null(fs_defaults, user), list())                                   # lines 674-677
u <- filter_valid_args(legacy, forestsearch,
       exclude = c("df.analysis", "confounders.name", "details", "plot.sg", "quiet"))
cat("legacy path : n.min entry present =", "n.min" %in% names(legacy), "; passed to forestsearch() =", "n.min" %in% names(u), "\n")
cat("methods path: n.min entry present =", "n.min" %in% names(meth), "; is NULL =", is.null(meth[["n.min"]]), "\n")
cat("forestsearch formal n.min =", formals(forestsearch)$n.min, "\n")
cat("REPRO_VERDICT:", if (!("n.min" %in% names(u)) && "n.min" %in% names(meth) && is.null(meth[["n.min"]])) "REPRODUCED" else "NOT_REPRODUCED", "\n")
```

```text
legacy path : n.min entry present = FALSE ; passed to forestsearch() = FALSE 
methods path: n.min entry present = TRUE ; is NULL = TRUE 
forestsearch formal n.min = 60 
REPRO_VERDICT: REPRODUCED
```

### V178

- Finding: A-137 — An error raised by forestsearch() (or the standalone GRF call) inside a simulation replicate is downgraded to a warning and the replicate is returned as an ordi
- Runs: exit 0, 2.3 s. Verdict: CONFIRMED (failed forestsearch() call returned as an unmarked non-detection row)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(3)
n <- 200
df <- data.frame(id = 1:n, y_sim = rexp(n), event_sim = rbinom(n, 1, 0.7),
                 treat_sim = rep(0:1, n / 2), v1 = rnorm(n),
                 flag_harm = c(rep(1L, 40), rep(0L, 160)))
w <- character(0)
out <- tryCatch(withCallingHandlers(
  forestsearch:::.run_fs_analysis_gen(
    data = df, confounders_name = "v1",
    params = list(outcome.name = "y_sim", event.name = "event_sim",
                  treat.name = "treat_sim", subgroup_method = "bogus"),  # forestsearch() errors in match.arg
    dgm = NULL, cox_formula = NULL, cox_formula_adj = NULL,
    outcome_name = "y_sim", event_name = "event_sim", treat_name = "treat_sim",
    harm_col = "flag_harm", label = "FS", verbose = FALSE),
  warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") }),
  error = function(e) e)
print(w)
if (inherits(out, "error")) {
  cat("REPRO_VERDICT: NOT_REPRODUCED (helper itself errored:", conditionMessage(out), ")\n")
} else {
  o <- as.data.frame(out)
  print(o[, intersect(c("analysis", "any.H", "size.H", "sens", "ppv", "spec", "npv", "sg.def", "fs_minutes"), names(o))])
  failed  <- any(grepl("FS failed", w, fixed = TRUE))
  looks_ok <- nrow(o) == 1 && isTRUE(o$any.H == 0) && isTRUE(o$size.H == 0) && isTRUE(o$sens == 0)
  marked  <- any(grepl("status|error|fail", names(o), ignore.case = TRUE))
  cat("REPRO_VERDICT:", if (failed && looks_ok && !marked) "CONFIRMED (failed forestsearch() call returned as an unmarked non-detection row)" else "NOT_REPRODUCED", "\n")
}
```

```text
[1] "FS failed: 'arg' should be one of “consistency”, “dina”, “grf”"
  analysis any.H size.H sens ppv spec npv sg.def fs_minutes
1       FS     0      0    0   0  0.8   1                NA
REPRO_VERDICT: CONFIRMED (failed forestsearch() call returned as an unmarked non-detection row)
```

### V179

- Finding: A-139 — validate_subgroups() evaluates each subgroup expression with enclos = parent.frame(), so a free name that exists only in the user's workspace (not in the data o
- Runs: exit 0, 2.3 s. Verdict: CONFIRMED (validation accepts a name that exists only in the calling session)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
d  <- data.frame(age = c(40, 50, 60, 70))
sg <- list(list(id = "age <= age_med", name = "Young", grp = "Clinical"))
assign("age_med", 55, envir = globalenv())   # workspace only; NOT passed through `cutpoints`
runner <- function() validate_subgroups(sg, d, cutpoints = list(), benchmarks = NULL)
passed <- tryCatch(isTRUE(runner()), error = function(e) { cat("validation error:", conditionMessage(e), "\n"); FALSE })
# what a fresh worker can resolve: data columns and base only
worker_sees <- tryCatch({ eval(parse(text = sg[[1]]$id), envir = d, enclos = baseenv()); TRUE },
                        error = function(e) { cat("worker-side error:", conditionMessage(e), "\n"); FALSE })
cat("validation passed:", passed, "; resolvable without the workspace:", worker_sees, "\n")
cat("REPRO_VERDICT:", if (passed && !worker_sees) "CONFIRMED (validation accepts a name that exists only in the calling session)" else "NOT_REPRODUCED", "\n")
rm("age_med", envir = globalenv())
```

```text
worker-side error: object 'age_med' not found 
validation passed: TRUE ; resolvable without the workspace: FALSE 
REPRO_VERDICT: CONFIRMED (validation accepts a name that exists only in the calling session)
```

### V180

- Finding: A-141 — With cens_type = "uniform", setup_gbsg_dgm() builds model_params$censoring with no min/max, so simulate_from_dgm() calls runif(n, min = NULL, max = NULL) and st
- Runs: exit 0, 2.3 s. Verdict: CONFIRMED (uniform-censoring GBSG DGM cannot be simulated)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
dgm_u <- tryCatch(setup_gbsg_dgm(model = "null", cens_type = "uniform", verbose = FALSE),
                  error = function(e) e)
if (inherits(dgm_u, "error")) {
  cat("setup_gbsg_dgm itself failed:", conditionMessage(dgm_u), "\n")
  cat("REPRO_VERDICT: NOT_REPRODUCED (failed earlier than claimed)\n")
} else {
  str(dgm_u$model_params$censoring)
  r <- tryCatch(simulate_from_dgm(dgm_u, n = 100, seed = 1), error = function(e) e)
  if (inherits(r, "error")) cat("simulate_from_dgm error:", conditionMessage(r), "\n")
  no_range <- is.null(dgm_u$model_params$censoring$min) && is.null(dgm_u$model_params$censoring$max)
  cat("REPRO_VERDICT:", if (inherits(r, "error") && no_range) "CONFIRMED (uniform-censoring GBSG DGM cannot be simulated)" else "NOT_REPRODUCED", "\n")
}
```

```text
List of 3
 $ type: chr "uniform"
 $ mu  : NULL
 $ tau : NULL
simulate_from_dgm error: invalid arguments 
REPRO_VERDICT: CONFIRMED (uniform-censoring GBSG DGM cannot be simulated)
```

### V181

- Finding: A-142 — A strata_rand value that is not a column of the sampled data is silently ignored and marginal allocation is used in place of the requested stratified allocation
- Runs: exit 0, 2.3 s. Verdict: CONFIRMED (unknown strata_rand silently gives marginal allocation)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
N <- 200
dgm <- structure(list(
  df_super = data.frame(id = 1:N, treat = rep(0:1, N / 2), stratum = rep(c("a", "b"), each = N / 2),
                        lin_pred_1 = 0.2, lin_pred_0 = 0, lin_pred_cens_1 = 0, lin_pred_cens_0 = 0),
  model_params = list(mu = 3, tau = 1, censoring = list(type = "weibull", mu = 4, tau = 1))),
  class = c("aft_dgm_flex", "list"))
w <- character(0)
res <- tryCatch({
  a <- withCallingHandlers(
    simulate_from_dgm(dgm, n = 100, seed = 1, strata_rand = "no_such_column"),
    warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
  b <- simulate_from_dgm(dgm, n = 100, seed = 1)                    # marginal allocation
  s <- simulate_from_dgm(dgm, n = 100, seed = 1, strata_rand = "stratum")  # genuine stratified allocation
  list(same_as_marginal = identical(a$treat_sim, b$treat_sim),
       differs_from_stratified = !identical(a$treat_sim, s$treat_sim))
}, error = function(e) e)
if (inherits(res, "error")) {
  cat("error raised:", conditionMessage(res), "\nREPRO_VERDICT: NOT_REPRODUCED (an error is raised)\n")
} else {
  print(res); cat("warnings:", length(w), "\n")
  cat("REPRO_VERDICT:", if (isTRUE(res$same_as_marginal) && length(w) == 0) "CONFIRMED (unknown strata_rand silently gives marginal allocation)" else "NOT_REPRODUCED", "\n")
}
```

```text
$same_as_marginal
[1] TRUE

$differs_from_stratified
[1] TRUE

warnings: 0 
REPRO_VERDICT: CONFIRMED (unknown strata_rand silently gives marginal allocation)
```

### V182

- Finding: A-144 — With subgroup_notation = "benefit" and a DGM that has no stored CDE, the CDE truth is recomputed from the raw super-population after the DGM has already been in
- Runs: exit 0, 2.6 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(1); N <- 2000
sup <- data.frame(flag.harm = rep(c(1L, 0L), c(400, 1600)))
sup$theta_0 <- rnorm(N, 0, 0.3)
sup$theta_1 <- sup$theta_0 + ifelse(sup$flag.harm == 1, log(2), log(0.6))  # CDE(H) = 2 exactly
dgm <- list(df_super_rand = sup, hr_H_true = 2, hr_Hc_true = 0.6, hr_causal = 0.75,
            hazard_ratios = list(overall = 0.75, harm_subgroup = 2, no_harm_subgroup = 0.6),
            model_type = "alt")
res <- data.frame(analysis = "FSlg", any.H = 1L, size.H = 100, size.Hc = 400,
                  hr.H.hat = c(1.9, 2.0, 2.1, 2.0, 2.0, 2.0), hr.Hc.hat = 0.6)
dh <- as.data.frame(build_estimation_table(res, dgm, subgroup_notation = "harm")[["_data"]])
db <- as.data.frame(build_estimation_table(res, dgm, subgroup_notation = "benefit")[["_data"]])
print(dh); print(db)
col <- grep("‡", names(db), value = TRUE)[1]
bh <- as.numeric(dh[1, col]); bb <- as.numeric(db[1, col])
cat("b-ddagger, H row: harm notation =", bh, "; benefit notation =", bb, "\n")
cat("REPRO_VERDICT:", if (isTRUE(abs(bh) < 5) && isTRUE(bb < -50)) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
  Estimator Avg   SD Min Max b‡ (%) b† (%)
1      θ̂(Ĥ) 2.0 0.06 1.9 2.1      0      0
2     θ̂(Ĥᶜ) 0.6 0.00 0.6 0.6      0      0
                                             Subgroup
1 Ĥ: 6 estimable, avg |Ĥ| = 100, θ†(H) = 2, θ‡(H) = 2
2      Ĥᶜ: avg |Ĥᶜ| = 400, θ†(Hᶜ) = 0.6, θ‡(Hᶜ) = 0.6
  Estimator  Avg   SD  Min  Max b‡ (%) b† (%)
1      θ̂(Ĝ) 0.50 0.02 0.48 0.53 -74.98   0.08
2     θ̂(Ĝᶜ) 1.67 0.00 1.67 1.67 177.78   0.00
                                               Subgroup
1 Ĝ: 6 estimable, avg |Ĝ| = 100, θ†(G) = 0.5, θ‡(G) = 2
2       Ĝᶜ: avg |Ĝᶜ| = 400, θ†(Gᶜ) = 1.67, θ‡(Gᶜ) = 0.6
b-ddagger, H row: harm notation = 0 ; benefit notation = -74.98 
REPRO_VERDICT: CONFIRMED
```

### V183

- Finding: A-145 — The generated narrative contains fixed evaluative sentences (low detection rate confirms type-I error control; bias is upward; bias correction improves accuracy
- Runs: exit 0, 2.3 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
res <- data.frame(analysis = "GRF", any.H = rep(c(1L, 0L), c(40, 60)),
                  size.H = 80, size.Hc = 420, hr.H.hat = 0.7, hr.Hc.hat = 1.1)
dgm <- list(hr_H_true = NA_real_, hr_Hc_true = 1, hr_causal = 1, model_type = "null")
txt <- interpret_estimation_table(res, dgm, analysis_method = "GRF", scenario = "null", cat = FALSE)
cat(txt, "\n")
ok <- grepl("40.0%", txt, fixed = TRUE) &&
  grepl("low detection rate confirms controlled type-I error", txt, fixed = TRUE) &&
  grepl("-30.0% relative bias above", txt, fixed = TRUE) && grepl("upward bias", txt, fixed = TRUE)
cat("REPRO_VERDICT:", if (isTRUE(ok)) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
Under the null hypothesis (true Cox HR = 1 uniformly), 40 of 100 simulations (40.0%) identified a subgroup using GRF. This low detection rate confirms controlled type-I error. Among those 40 false detections, the identified subgroup averaged 80 patients.

The naive Cox HR in the identified subgroup averaged 0.7 (SD = 0), representing -30.0% relative bias above the true value of 1. This upward bias reflects selection: the algorithm identified whichever patients happened to look most like a harm subgroup by chance. In the complement, the Cox HR averaged 1.1 (10.0% bias), showing the expected mirror effect where removing the worst-looking patients makes the remainder appear modestly better.

These results underscore that under the null, the few false detections produce highly biased estimates, reinforcing the need for bootstrap bias correction for any subgroup identified by a data-driven search. 
REPRO_VERDICT: CONFIRMED
```

### V184

- Finding: A-146 — For maxk >= 3 the agreement table labels each bootstrap subgroup from M.1 and M.2 only, so three-factor subgroups lose their third factor and distinct subgroups
- Runs: exit 0, 2.3 s. Verdict: CONFIRMED (4 distinct 3-factor subgroups shown as 1 two-factor row)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
res <- data.frame(boot_id = 1:4, any_found = 1L, Pcons = 0.92, hr_sg = 2, N_sg = 80, K_sg = 3,
  M.1 = "{a <= 1}", M.2 = "{b <= 1}",
  M.3 = c("{c <= 1}", "{d <= 1}", "{e <= 1}", "{f <= 1}"), stringsAsFactors = FALSE)
out <- tryCatch(suppressWarnings(summarize_bootstrap_subgroups(res, nb_boots = 4, maxk = 3)), error = function(e) e)
if (inherits(out, "error")) {
  cat("REPRO_VERDICT: INCONCLUSIVE -", conditionMessage(out), "\n")
} else {
  ag <- as.data.frame(out$agreement); print(ag)
  merged <- nrow(ag) == 1L && ag$N[1] == 4L && !grepl("c <= 1", ag$Subgroup[1], fixed = TRUE)
  cat("REPRO_VERDICT:", if (merged) "CONFIRMED (4 distinct 3-factor subgroups shown as 1 two-factor row)" else "NOT REPRODUCED", "\n")
}
```

```text
             Subgroup K_sg N Percent_of_successful Rank
1 {a <= 1} & {b <= 1}    3 4                   100    1
REPRO_VERDICT: CONFIRMED (4 distinct 3-factor subgroups shown as 1 two-factor row)
```

### V185

- Finding: A-162 — ordinal_vars is honoured only for variables that are also listed in cat_vars, so in the documented call pattern and in generate_gbsg_bootstrap_general() the ord
- Runs: exit 1, 2.3 s. Verdict: no verdict line (exit 1)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
# noise_level = 1 => cat_flip_prob = 0.5: about half of the grade values should move if grade were perturbed
d <- suppressMessages(generate_gbsg_bootstrap_general(n = 2000, seed = 1, noise_level = 1))
g <- survival::gbsg
set.seed(1); idx <- sample(nrow(g), 2000, replace = TRUE)
same_grade <- identical(as.numeric(d$grade), as.numeric(g$grade[idx]))
meno_moved <- mean(d$meno != g$meno[idx])
cat("grade identical to plain bootstrap resample:", same_grade, "\n")
cat("share of meno changed (a cat_var, for contrast):", round(meno_moved, 3), "\n")
cat("REPRO_VERDICT:", if (same_grade && meno_moved > 0.3) "CONFIRMED (ordinal variable not perturbed, categorical one is)" else "NOT REPRODUCED", "\n")
```

```text
Error in generate_bootstrap_synthetic(data = gbsg, continuous_vars = c("age",  : 
  'data' must be a data frame
Calls: source ... generate_gbsg_bootstrap_general -> generate_bootstrap_synthetic
In addition: Warning message:
In data("gbsg", package = "survival", envir = environment()) :
  data set ‘gbsg’ not found
Execution halted
```

### V186

- Finding: A-148 — The ordinal adjacent-value step calls sample(possible_new, 1) on what is a length-1 numeric at either end of the scale, so base sample() draws from 1:possible_n
- Runs: exit 0, 2.2 s. Verdict: CONFIRMED (non-level values written by sample() on a length-1 numeric)

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
d <- data.frame(g = rep(c(10, 20, 30), 100), x = rnorm(300))
s <- suppressMessages(generate_bootstrap_synthetic(d, continuous_vars = "x", cat_vars = "g",
        ordinal_vars = "g", seed = 1, noise_level = 1))
bad <- sort(setdiff(unique(s$g), c(10, 20, 30)))
cat("values outside the original level set {10,20,30}:", if (length(bad)) paste(bad, collapse = " ") else "none", "\n")
cat("REPRO_VERDICT:", if (length(bad) > 0) "CONFIRMED (non-level values written by sample() on a length-1 numeric)" else "NOT REPRODUCED", "\n")
```

```text
values outside the original level set {10,20,30}: 1 2 3 4 5 6 7 8 9 11 12 13 14 15 16 17 18 19 
REPRO_VERDICT: CONFIRMED (non-level values written by sample() on a length-1 numeric)
```

### V187

- Finding: A-149 — A numeric two-valued categorical variable is flipped with 1 - x, which is only right for 0/1 coding; any other numeric coding (for example 1/2) silently yields 
- Runs: exit 0, 2.3 s. Verdict: CONFIRMED - values outside the original 1/2 coding: -1, 0

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
d <- data.frame(sex = rep(c(1, 2), 50), x = seq_len(100) / 10)
s <- generate_bootstrap_synthetic(d, continuous_vars = "x", cat_vars = "sex", seed = 1, noise_level = 1)
print(table(s$sex))
bad <- setdiff(unique(s$sex), unique(d$sex))
cat("REPRO_VERDICT:", if (length(bad) > 0) "CONFIRMED" else "NOT_REPRODUCED",
    "- values outside the original 1/2 coding:", paste(sort(bad), collapse = ", "), "\n")
```

```text

-1  0  1  2 
25 27 24 24 
REPRO_VERDICT: CONFIRMED - values outside the original 1/2 coding: -1, 0
```

### V188

- Finding: A-150 — With continuous_vars and cat_vars both left NULL (the default), every numeric column is classed continuous and then removed from the auto-detected categorical s
- Runs: exit 0, 2.4 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
g <- survival::gbsg
vars <- c("meno", "hormon", "status", "grade")
s_auto <- generate_bootstrap_with_noise(g, seed = 1, noise_level = 0.1)
s_cont <- generate_bootstrap_with_noise(g, seed = 1, noise_level = 0.1,
  continuous_vars = c("age", "size", "nodes", "pgr", "er", "rfstime"))
set.seed(1); idx <- sample(nrow(g), nrow(g), replace = TRUE)
chg_auto <- sapply(vars, function(v) mean(s_auto[[v]] != g[[v]][idx]))
chg_cont <- sapply(vars, function(v) mean(s_cont[[v]] != g[[v]][idx]))
print(rbind(both_NULL = chg_auto, continuous_given = chg_cont))  # documented rate: 0.05
cat("REPRO_VERDICT:", if (all(chg_auto == 0) && all(chg_cont > 0)) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
                       meno     hormon     status      grade
both_NULL        0.00000000 0.00000000 0.00000000 0.00000000
continuous_given 0.04810496 0.06268222 0.04810496 0.05102041
REPRO_VERDICT: CONFIRMED
```

### V189

- Finding: A-152 — With tol below 1 the bisection loop never terminates, because once the integer bounds are one apart the rounded midpoint equals one of the bounds and neither bo
- Runs: exit 0, 4.8 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
stopifnot(requireNamespace("cubature", quietly = TRUE))
seen <- character()
setTimeLimit(elapsed = 90, transient = TRUE)   # fallback only; the message guard should end the call first
res <- tryCatch(withCallingHandlers(
  find_required_sample_size(theta = 2, tol = 0.5, verbose = TRUE),
  message = function(m) {
    seen <<- c(seen, trimws(conditionMessage(m)))
    if (length(seen) >= 30) stop("LOOP_GUARD")   # bisection on 20..500 needs at most 9 steps
    invokeRestart("muffleMessage")
  }), error = function(e) conditionMessage(e))
setTimeLimit(elapsed = Inf)
print(tail(seen, 5))
cat("REPRO_VERDICT:", if (identical(res, "LOOP_GUARD") && length(unique(tail(seen, 10))) == 1) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
[1] "n = 40, power = 0.796" "n = 40, power = 0.796" "n = 40, power = 0.796"
[4] "n = 40, power = 0.796" "n = 40, power = 0.796"
REPRO_VERDICT: CONFIRMED
```

### V190

- Finding: INV-2 — The Monte Carlo branch calls set.seed(seed) (default 42) and never restores the caller's RNG state, so compute_detection_probability_glm(method = "monte_carlo")
- Runs: exit 0, 2.4 s. Verdict: CONFIRMED

```r
pkgload::load_all("/home/larryleon/Documents/GitHub/forestsearch", quiet = TRUE)
set.seed(123); x1 <- runif(1)
set.seed(123); invisible(compute_detection_probability_glm(theta = 2, d_eff = 30, method = "monte_carlo")); x2 <- runif(1)
print(c(expected = x1, after_call = x2))   # differ: the caller's stream was replaced
set.seed(1); invisible(compute_detection_probability_glm(theta = 2, d_eff = 30, method = "monte_carlo")); u1 <- runif(1)
set.seed(2); invisible(compute_detection_probability_glm(theta = 1.5, d_eff = 30, method = "monte_carlo")); u2 <- runif(1)
print(c(u1 = u1, u2 = u2))                 # identical: every call leaves the same fixed state
cat("REPRO_VERDICT:", if (!identical(x1, x2) && identical(u1, u2)) "CONFIRMED" else "NOT_REPRODUCED", "\n")
```

```text
  expected after_call 
 0.2875775  0.1021029 
       u1        u2 
0.1021029 0.1021029 
REPRO_VERDICT: CONFIRMED
```

