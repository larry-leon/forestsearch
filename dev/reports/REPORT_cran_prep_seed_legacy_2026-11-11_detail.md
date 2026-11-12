# CRAN preparation, part 2: detail of Parts 1.1, 2 and 3, 2026-11-11

Companion to `dev/reports/REPORT_cran_prep_seed_legacy_2026-11-11.md`. Every row below is the result returned by an independent adversarial verifier, which re-derived the investigator's result from the source at `c40d5309` (Parts 2, 3), whose `R/`, `vignettes/` and `quarto/` are identical to the start SHA `3fbbc55c`, or at `3fbbc55c` itself (Part 1.1). The verifier's corrections to the investigator follow each item. Nothing was changed for Parts 2 and 3.

## Part 1.1 — every set.seed() site under R/ (start SHA 3fbbc55c)


The rows below give one entry per set.seed() site in package code under R/ of the forestsearch R package, 44 sites in all, at the start SHA 3fbbc55c. They were classified by four agents, each independently verified. Every field is reproduced verbatim from the verified rows, in the inventory's order; the verifiers' corrections follow at the end.

### R/bootstrap_analysis_dofuture.R:334 — bootstrap_results

- **site:** R/bootstrap_analysis_dofuture.R:334
- **fun:** bootstrap_results
- **seed_source:** Formal `seed = 8316951L` of bootstrap_results (R/bootstrap_analysis_dofuture.R:302). On a direct call the user sets it. From forestsearch_bootstrap_dofuture it is that function's own formal `seed = 8316951L` (R/bootstrap_dofuture_main.R:213), passed through `boot_seed <- seed` (:388) and `seed = boot_seed` (:492). The user can change it from both exported entry points.
- **classification:** argument
- **restored:** no. set.seed(seed) runs unconditionally on the main process, and the caller's .Random.seed is never saved. The on.exit in forestsearch_bootstrap_dofuture (R/bootstrap_dofuture_main.R:359-370) restores only options(warn) and future::plan (plus a gc).
- **exported_reaching:** bootstrap_results (direct); forestsearch_bootstrap_dofuture -> bootstrap_results(seed = boot_seed) at R/bootstrap_dofuture_main.R:484-497. No other internal caller. The static list is correct.
- **notes:** There is no is.null guard. seed = NULL is documented as allowed (R/bootstrap_dofuture_main.R:35-36), and set.seed(NULL) then re-initialises the RNG from entropy, whereas .make_boot_index_matrix skips the call. forestsearch_bootstrap_dofuture always passes boot_index_mat, so this seed does not drive the bootstrap indices (doc :86-89; branch at :430-434 area). It does fix the RNG state from which `%dofuture%` with .options.future seed = TRUE (:408-413) derives the per-iteration L'Ecuyer streams. On a direct call with boot_index_mat = NULL, those streams drive the per-iteration sample.int. Each replicate's inner forestsearch() re-seeds its own steps with fs.est$args_call_all$seedit (for example get_FSdata, restored at R/get_FSdata_helpers.R:249-259).

### R/bootstrap_calculations_helpers.R:27 — .make_boot_index_matrix

- **site:** R/bootstrap_calculations_helpers.R:27
- **fun:** .make_boot_index_matrix
- **seed_source:** Formal `seed = NULL` of the internal helper (:26), guarded by `if (!is.null(seed))`. From bootstrap_ystar it is that function's formal `seed = 8316951L` (:103), passed positionally at :104. From forestsearch_bootstrap_dofuture it is that function's formal `seed = 8316951L` (R/bootstrap_dofuture_main.R:213), via `boot_seed <- seed` and `seed = boot_seed` (:388-393). The user can change it from both exported entry points.
- **classification:** argument
- **restored:** no. Neither the helper nor its two callers saves or restores .Random.seed. In forestsearch_bootstrap_dofuture the state is then re-seeded with the same value at R/bootstrap_analysis_dofuture.R:334.
- **exported_reaching:** bootstrap_ystar -> .make_boot_index_matrix (R/bootstrap_calculations_helpers.R:104); forestsearch_bootstrap_dofuture -> .make_boot_index_matrix (R/bootstrap_dofuture_main.R:389-393). The static list is correct.
- **notes:** This matrix is the single source of the bootstrap row indices for forestsearch_bootstrap_dofuture. It is drawn on the main process with the current-kind (normally Mersenne-Twister) sample.int. With seed = NULL nothing is seeded and the running stream is consumed.

### R/consistency_resample.R:459 — consistency_resample

- **site:** R/consistency_resample.R:459
- **fun:** consistency_resample
- **seed_source:** Formal `seed = NULL` (R/consistency_resample.R:417), guarded by `!is.null(seed)`. It runs only inside `if (method %in% c("mc","both"))` (:458). The direct-call default method is "both" via match.arg (:410, :419), and the line is not reached when pieces is NULL (early return :437). On a direct call the user sets the seed. From consistency_resample_compare it is that function's formal `seed = NULL` (:494), passed as `seed = seed` with method = "both" (:514-518). Every internal package caller passes method = "closed" and no seed: R/subgroup_consistency_helpers.R:1487-1490 and :1508-1518 (evaluate_subgroup_consistency), and :1870-1883 (evaluate_consistency_twostage). From forestsearch() and everything above it, this set.seed never executes.
- **classification:** argument
- **restored:** no. Nothing saves or restores .Random.seed.
- **exported_reaching:** Effectively only consistency_resample (direct) and consistency_resample_compare (R/consistency_resample.R:514-518). The following reach the function only through the closed-form calls, never the set.seed line: evaluate_subgroup_consistency, evaluate_consistency_twostage, subgroup.consistency (via forestsearch's do.call at R/forestsearch_main.R:3645), forestsearch, and the forestsearch callers (bootstrap, CV, fpr_calibration, compare_selection_rules, run_null_calibration, run_simulation_analysis, fs_fdr_report, mrct_*). Static over-reach: fs_declaration_calibration, fs_declaration_table, fs_mr_inference, fs_family_report and fs_oc_* do not reach consistency_resample at all. They only read formals() of subgroup.consistency or forestsearch.
- **notes:** The closed-form screen used by forestsearch(consistency_method = "resample") (the default, R/forestsearch_main.R:1706) draws no random numbers here. The MC multiplier draws (rbinom/rnorm/rpois, :461-465) happen only for method mc or both, which only a direct call or a _compare call requests.

### R/consistency_resample.R:505 — consistency_resample_compare

- **site:** R/consistency_resample.R:505
- **fun:** consistency_resample_compare
- **seed_source:** Formal `seed = NULL` of consistency_resample_compare (R/consistency_resample.R:494), guarded by `!is.null(seed)`. The user can change it from the exported function. No package code calls this function.
- **classification:** argument
- **restored:** no. Nothing saves or restores .Random.seed.
- **exported_reaching:** consistency_resample_compare only (no internal callers). The static list is correct.
- **notes:** The seed governs the R_true literal splits (run_single_consistency_split's sample(), R/subgroup_consistency_helpers.R:272) at :506-510. The same value is then passed to consistency_resample(seed = seed, method = "both") (:514-518), which re-seeds at :459. The MC multiplier draws therefore restart from the RNG state the literal splits started from. With seed = NULL both stages use the running stream.

### R/dina_bagged.R:436 — .dina_make_bootstrap_indices

- **site:** R/dina_bagged.R:436
- **fun:** .dina_make_bootstrap_indices
- **seed_source:** Formal `seed = NULL` of the internal helper (:422), guarded by `if (!is.null(seed))` (:423). From dina_fit_bagged it is that function's formal `seed = NULL` (R/dina_bagged.R:141), passed at :194. From dina_bagged it is forwarded through `...` to dina_fit_bagged (R/dina_df.R:179-183). From dina_subgroup_bootstrap it is that function's formal `seed = NULL` (R/dina_subgroup_bootstrap.R:289), passed unconditionally at :535. The user can change it from all three exported entry points.
- **classification:** argument
- **restored:** yes, when seed is non-NULL. on.exit at R/dina_bagged.R:424-434 saves globalenv .Random.seed and assigns it back, or rm()s it when none existed at entry. This is the same pattern as R/dina.R:421-433.
- **exported_reaching:** dina_fit_bagged (:194); dina_bagged -> dina_fit_bagged via `...`; dina_subgroup_bootstrap (R/dina_subgroup_bootstrap.R:535). The static list is correct.
- **notes:** With the default seed = NULL nothing is seeded and the current stream is consumed without restoration, because the restore sits inside the non-NULL branch. The per-bag and per-replicate parallel foreach loops use .options.future seed = TRUE (R/dina_bagged.R:206, R/dina_subgroup_bootstrap.R:544). Those streams derive from the restored caller state, not from `seed`. Per-bag dina_fit calls hard-code seed = NULL (R/dina_bagged.R:357).

### R/dina.R:433 — .dina_make_folds

- **site:** R/dina.R:433
- **fun:** .dina_make_folds
- **seed_source:** Formal `seed = NULL` of the helper (R/dina.R:419), guarded by !is.null(seed) (:420). It is passed from dina_fit's formal `seed = NULL` (:167) at :241. That call is made only when cross_fitting is TRUE and n_folds != 1, which holds by default (:161-162, :226). By entry point: (a) dina_fit: user formal. (b) dina: user passes seed through `...` (R/dina_df.R:98-102). (c) forestsearch: seed = dina_args$seed if supplied, else forestsearch's `seedit = 8316951` (R/forestsearch_main.R:1640), resolved in .resolve_dina_args at R/forestsearch_helpers.R:1154-1157 (seed key :1156, get_arg :1105). This covers use_dina screening (R/forestsearch_main.R:3189-3215, fit only when dina_cuts and dina_res are NULL) and subgroup_method = "dina" (R/forestsearch_main.R:2766-2774 -> R/forestsearch_helpers.R:1490-1537, fit only when dina_res is NULL). The user can set it. (d) dina_subgroup_bootstrap: its formal `seed`, used only for the point refit, and only when both sg and fit are NULL (R/dina_subgroup_bootstrap.R:415, :461, :476-482). Per-replicate refits strip the seed (:397, :799-803), so they pass NULL. (e) bootstrap_results / forestsearch_bootstrap_dofuture / forestsearch_Kfold / forestsearch_tenfold re-run forestsearch with args_FS_template or cv_args = fs.est$args_call_all, with dina_res and dina_cuts nulled (R/bootstrap_analysis_dofuture.R:581-582; R/forestsearch_cross_validation.R:415-416, :905-906). The value is the original call's seedit or dina_args$seed, stored in the object. (f) fpr_calibration passes seedit = seed_i = sim_seed + offset + i (R/fpr_calibration.R:298, :307), argument-derived; DINA is enabled only via fs_params. Other forestsearch callers (compare_selection_rules dots, run_simulation_analysis fs_params, mrct_* fs_args/fs_spec) let the user set seedit or dina_args. run_null_calibration and fs_fdr_report set seedit from their own seed arguments (R/calibrate_null_correction.R:352, R/fs_fdr_report.R:358).
- **classification:** argument
- **restored:** yes, when seed is non-NULL. on.exit at R/dina.R:421-432 saves globalenv .Random.seed and assigns it back, or rm()s it when none existed at entry.
- **exported_reaching:** dina_fit; dina; dina_subgroup_bootstrap (point refit only, when sg and fit are NULL); forestsearch (subgroup_method = "dina" or use_dina = TRUE); and, through actual do.call(forestsearch) sites: bootstrap_results, forestsearch_bootstrap_dofuture, forestsearch_Kfold, forestsearch_tenfold, fpr_calibration, compare_selection_rules, run_null_calibration, run_simulation_analysis, mrct_case_study, mrct_trial_analysis, mrct_region_sims, fs_fdr_report. Removed from the static list: dina_bagged and dina_fit_bagged (.dina_fit_one_bag hard-codes seed = NULL, R/dina_bagged.R:357). Also removed: fpr_calibration_transfer (fs_params = list() is hard-coded at R/fpr_calibration.R:548, fpr_calibration has no `...`, and forestsearch defaults use_dina = FALSE / subgroup_method = "consistency", R/forestsearch_main.R:1646,1651). Also removed: fs_family_report and fs_oc_family_enumerate / fs_oc_grid / fs_oc_invert / fs_oc_predict, which only read formals(forestsearch) and never call it.
- **notes:** One row with mixed provenance. At the site the value is a formal, so it is classed as argument. Through bootstrap and CV it is stored in the object (fs.est$args_call_all$seedit or $dina_args$seed). Through the calibration loops it is derived per replicate. Within one bootstrap or CV run, every replicate's DINA refit re-seeds with the same value, so the cross-fitting fold assignment (sample.int(n)) is identical across replicates of equal n.

### R/forestsearch_cross_validation.R:370 — forestsearch_Kfold

- **site:** R/forestsearch_cross_validation.R:370
- **fun:** forestsearch_Kfold
- **seed_source:** Formal `seedit = 8316951L` of forestsearch_Kfold (R/forestsearch_cross_validation.R:291). The user can change it. No package code calls forestsearch_Kfold.
- **classification:** argument
- **restored:** no. The on.exit at :334 restores only future::plan, and .Random.seed is never saved. The call runs on the main process.
- **exported_reaching:** forestsearch_Kfold only (no internal callers). The static list is correct.
- **notes:** It executes only when Kfolds < number of observations (:369). The default Kfolds = .fs_cv_nobs(fs.est) (:290) is leave-one-out, so with defaults the set.seed never runs. The seed controls only the scramble and fold assignment (:371-375). The per-fold forestsearch() refits use cv_args$seedit from fs.est$args_call_all (cv_args <- fs_args, :388), a separate stored value that this formal does not override.

### R/forestsearch_cross_validation.R:959 — forestsearch_tenfold (inside the %dofuture% foreach body over ksim)

- **site:** R/forestsearch_cross_validation.R:959
- **fun:** forestsearch_tenfold (inside the %dofuture% foreach body over ksim)
- **seed_source:** seed + 1000 * ksim. `seed` is the formal `seed = 8316951L` of forestsearch_tenfold (R/forestsearch_cross_validation.R:809), which the user can change. ksim is the foreach iterator seq_len(sims) (:952).
- **classification:** argument-derived
- **restored:** no, from package code. The on.exit at :842 restores only future::plan. The call runs inside the %dofuture% task (.options.future seed = TRUE, :951-955). On the default multisession plan (:810) it changes only the worker's RNG. Under a sequential plan, the future framework's own RNG reset (future:::evalFutureInternal) restores the caller's state, not this function.
- **exported_reaching:** forestsearch_tenfold only (no internal callers). The static list is correct.
- **notes:** It controls the per-simulation scramble and fold assignment (:958-964). Replicate ksim is reproducible regardless of worker scheduling. The inner per-fold forestsearch() refits use args_call_all$seedit, which is stored in the object.

### R/fpr_calibration.R:302 — fpr_calibration, nested closure .run_batch (inside the %dofuture% body)

- **site:** R/fpr_calibration.R:302
- **fun:** fpr_calibration, nested closure .run_batch (inside the %dofuture% body)
- **seed_source:** seed_i <- sim_seed + seed_offset + i (R/fpr_calibration.R:298). sim_seed is the formal `sim_seed = 42L` (:204), which the user can change. seed_offset is 0L for the calibration batch (:321) or nsims + 10000L for the verification batch (:350). i is the foreach index (:293). The line runs only in permutation mode, when null_fn is NULL (:218, :299-302). In dgm mode seed_i goes to null_fn(seed_i) instead. seed_i is also passed to forestsearch as seedit (:307).
- **classification:** argument-derived
- **restored:** no, from package code. Nothing saves .Random.seed. The on.exit at :274-277 restores the plan, and only when n_workers > 1. The call runs inside the %dofuture% task (.options.future seed = TRUE, :292-296). With the default n_workers = 1L it runs under the caller's current plan, where any RNG reset comes from the future framework.
- **exported_reaching:** fpr_calibration (direct). fpr_calibration_transfer calls fpr_calibration via do.call (:557) with sim_seed = cal_result$call_args$sim_seed + 20000L + round(dc1[k]*1000) (:552-553). That value is stored in the object, plus a literal, plus a term derived from the delta / delta_c1 argument. In practice it cannot reach this line. It passes df.analysis = cal_result$call_args$df.analysis %||% eval(cal_result$call_args$null_fn_expr) (:546-547), and fpr_calibration never stores either field (call_args :384-391), so df.analysis is NULL and permutation mode stops at :221-222. The stop is swallowed by tryCatch at :557. A df.analysis supplied via `...` duplicates a named argument and errors. A null_fn supplied via `...` selects dgm mode.
- **notes:** Static over-reach: from fpr_calibration_transfer this line is reachable only if the user hand-edits cal_result$call_args$df.analysis (or null_fn_expr to yield a data frame).

### R/frontier_cis.R:164 — compute_frontier_cis

- **site:** R/frontier_cis.R:164
- **fun:** compute_frontier_cis
- **seed_source:** Formal `seed = NULL` (R/frontier_cis.R:109), guarded by `!is.null(seed)`. The line is reached only when the frontier is non-empty (:115-118) and df.est and args_call_all exist (:121-130). On a direct call the user sets it. From compare_selection_rules it is that function's formal `ci_seed = 1L` (R/compare_selection_rules.R:488), passed as `seed = ci_seed` at :596 when compute_cis = TRUE (the default, :486, :594). The user can change it from both exported entry points.
- **classification:** argument
- **restored:** no. Neither compute_frontier_cis nor compare_selection_rules saves or restores .Random.seed.
- **exported_reaching:** compute_frontier_cis (direct); compare_selection_rules -> compute_frontier_cis(seed = ci_seed) at R/compare_selection_rules.R:596, once per sg_focus x selection_rule combination. No other internal caller. The static list is correct.
- **notes:** The seed governs the per-split halving over the frontier members. compare_selection_rules re-seeds with the same ci_seed for every combination k, so every combination uses the same split stream. Direct calls default to NULL, meaning no seeding.

### R/fs_dgm_feasibility.R:181 — fs_dgm_feasibility (inner replicate loop)

- **site:** R/fs_dgm_feasibility.R:181
- **fun:** fs_dgm_feasibility (inner replicate loop)
- **seed_source:** s_i <- as.integer(seed) + i (R/fs_dgm_feasibility.R:180). `seed` is the formal `seed = NULL` (:123), which the user can change, and i is the replicate index (:179). The set.seed is skipped when seed is NULL.
- **classification:** argument-derived
- **restored:** partly. :162-166 saves globalenv .Random.seed unconditionally and restores it on exit, but only when it existed at entry. There is no rm branch, so a caller with no .Random.seed is left with the one that set.seed and the draws created. Compare R/dina.R:427-432.
- **exported_reaching:** fs_dgm_feasibility only (no internal callers). The static list is correct.
- **notes:** The same s_i is reused for every sample size n[j] (outer loop :174, inner :179), as documented at :73-74. s_i is also passed to simulate_from_glm_dgm(seed = s_i), which re-seeds with the same value (R/simulate_from_glm_dgm.R:97). The RNG kind is never changed (:161). Because the save is not gated on seed, a seed = NULL call with an existing .Random.seed also rolls back its consumed draws on exit, so repeated unseeded calls replay the same stream.

### R/fs_mr_field_uniform.R:93 — fs_mr_field_uniform (internal, not exported)

- **site:** R/fs_mr_field_uniform.R:93
- **fun:** fs_mr_field_uniform (internal, not exported)
- **seed_source:** Comes from the formal `seed` (default NULL, fs_mr_field_uniform.R:91), via `if (!is.null(seed)) set.seed(as.integer(seed))`. The only caller is fs_mr_inference.R:1020-1029, which passes `seed = if (!is.null(seed)) as.integer(seed) + 910000L else NULL` (line 1029). In a direct fs_mr_inference(seed=) call the user sets this seed (default NULL at fs_mr_inference.R:595, and NULL means no set.seed). From forestsearch(), the consistency branch passes `.g_mr(mr_inference_args$seed, seedit)` (forestsearch_main.R:3930), and the DINA and GRF branches pass `.g(mr_inference_args$seed, seedit)` (fs_mr_inference_methods.R:172). seedit defaults to 8316951 (forestsearch_main.R:1640), so the default is 9226951; the user can change it through seedit= or mr_inference_args$seed. Wrapper sources: mrct_case_study uses seed + k (mrct_case_study.R:297); fpr_calibration uses sim_seed + offset + i (fpr_calibration.R:298, 307); run_null_calibration uses seed_base + N*1000 + s (calibrate_null_correction.R:343, 352); fs_fdr_report uses seed_base + 100000*j + i (fs_fdr_report.R:287, 358, 409-410); bootstrap and CV use the value stored in fs.est$args_call_all. Gates: field_uniform = TRUE (formal FALSE at fs_mr_inference.R:599; FALSE at forestsearch_main.R:3891; no wrapper sets it), ci_method == "field" (fs_mr_inference.R:910), and at least 2 finite outer draws (fs_mr_inference.R:968).
- **classification:** argument-derived
- **restored:** no (neither fs_mr_field_uniform.R nor fs_mr_inference.R touches .Random.seed or uses on.exit)
- **exported_reaching:** fs_mr_inference (direct, field_uniform = TRUE); forestsearch (mr_inference = TRUE, mr_inference_args$field_uniform = TRUE; on the consistency branch, or on the DINA/GRF branches with mr_inference_args$ci_method = "field"); through forestsearch: compare_selection_rules (...), fpr_calibration (fs_params), fs_fdr_report (only if the fit stored field_uniform = TRUE), mrct_case_study, mrct_trial_analysis, mrct_region_sims (fs_args), run_null_calibration (fs_args), run_simulation_analysis (fs_params; front ends must be off on FS arms); bootstrap_results, forestsearch_bootstrap_dofuture, forestsearch_Kfold, forestsearch_tenfold only with mr_in_replicates = TRUE. Removed: fs_family_report, fs_oc_family_enumerate, fs_oc_grid, fs_oc_invert, fs_oc_predict (they only read formals(forestsearch): fs_family_report.R:109,123; fs_oc_family.R:217; fs_oc_predict.R:181; fs_oc_grid.R:480), and fpr_calibration_transfer (fs_params = list() at fpr_calibration.R:548; a user fs_params in ... would be a duplicate formal and an extra forestsearch arg would be an unused argument, so either do.call errors). Main path: forestsearch -> fs_mr_inference (forestsearch_main.R:3862) -> fs_mr_field_uniform (fs_mr_inference.R:1020).
- **notes:** The +910000L offset mirrors the field block's +900000L, and the sweep draws only after the field stream is used up. On overflow, set.seed(NA) errors and the tryCatch at fs_mr_inference.R:1019 turns that into field$uniform$note, but only for seeds in (integer.max - 910000, integer.max - 900000]. Larger seeds fail earlier at line 912.

### R/fs_mr_inference.R:615 — fs_mr_inference

- **site:** R/fs_mr_inference.R:615
- **fun:** fs_mr_inference
- **seed_source:** Comes from the formal `seed` (default NULL, fs_mr_inference.R:595), via `if (!is.null(seed)) set.seed(seed)`. A direct call is an argument, and NULL means unseeded. From forestsearch(): mr_inference_args$seed, else seedit (formal default 8316951, forestsearch_main.R:1640), passed at forestsearch_main.R:3930 (consistency branch) or through .fs_apply_mr(seedit = seedit) (forestsearch_main.R:2881, 3066 -> fs_mr_inference_methods.R:172) on the DINA and GRF branches. Per wrapper: compare_selection_rules passes it through ...; mrct_region_sims and run_simulation_analysis use fs_args$seedit / fs_params$seedit, else 8316951, constant across replicates; mrct_trial_analysis uses mr_inference_args$seed, else fs_args$seedit, else 8316951; mrct_case_study uses mr_inference_args$seed, else seed + k (mrct_case_study.R:293-297); fpr_calibration uses sim_seed + seed_offset + i (fpr_calibration.R:298, 307); run_null_calibration uses seed_base + N*1000 + s (calibrate_null_correction.R:343, 352); fs_fdr_report uses seed_base + 100000*j + i (fs_fdr_report.R:287, 358, 409-410) unless the fit stored mr_inference_args$seed; bootstrap and CV use the seedit / mr_inference_args$seed stored in fs.est$args_call_all, identical for every replicate.
- **classification:** argument
- **restored:** no. fs_mr_inference saves nothing. The only package-level restore on a reaching path is mrct_case_study's wrapper-level on.exit (mrct_case_study.R:263-271).
- **exported_reaching:** fs_mr_inference; forestsearch (mr_inference = TRUE); compare_selection_rules; fpr_calibration; fs_fdr_report (forces mr_inference = TRUE, fs_fdr_report.R:239); mrct_case_study and mrct_trial_analysis (mr_inference = TRUE by default); mrct_region_sims; run_null_calibration; run_simulation_analysis (FS arms need use_lasso, use_grf and use_dina set to FALSE in fs_params, or the dina-/grf- labels); bootstrap_results, forestsearch_bootstrap_dofuture, forestsearch_Kfold, forestsearch_tenfold only with mr_in_replicates = TRUE (otherwise forced FALSE at bootstrap_analysis_dofuture.R:607 and forestsearch_cross_validation.R:423, 914). Removed: fs_family_report, fs_oc_family_enumerate, fs_oc_grid, fs_oc_invert, fs_oc_predict (formals only), fpr_calibration_transfer (fs_params = list()). Main path: forestsearch(mr_inference = TRUE) -> fs_mr_inference (forestsearch_main.R:3862).
- **notes:** This is the first statement after the match.arg calls, so it reseeds the global stream before the multiplier draws. A direct call is unseeded by default. From forestsearch it is seeded by default (8316951). With the consistency method, MR together with any front end (use_lasso, use_grf or use_dina = TRUE) is rejected by the guard at forestsearch_main.R:1898-1901, before this site.

### R/fs_mr_inference.R:912 — fs_mr_inference (field block)

- **site:** R/fs_mr_inference.R:912
- **fun:** fs_mr_inference (field block)
- **seed_source:** Computed as as.integer(seed) + 900000L from the same `seed` formal as line 615 (default NULL, fs_mr_inference.R:595), so every source listed for line 615 applies, offset by 900000. The forestsearch default is 9216951. A user can change it through fs_mr_inference(seed=), forestsearch(seedit=) or mr_inference_args$seed. The site runs only inside `if (ci_method == "field")` (line 910). "field" is the default for a direct call (fs_mr_inference.R:594), for forestsearch's consistency branch (forestsearch_main.R:3888) and for the default mr_inference_args of mrct_case_study and mrct_trial_analysis. The DINA and GRF branches default to "ij" (fs_mr_inference_methods.R:171).
- **classification:** argument-derived
- **restored:** no
- **exported_reaching:** Same corrected list as line 615 (the DINA/GRF branches need ci_method = "field"). Removed: fs_family_report, the four fs_oc_* functions and fpr_calibration_transfer. Main path: forestsearch consistency branch, where ci_method defaults to "field".
- **notes:** fs_fdr_report's "ij" fallback (fs_fdr_report.R:240-242) fires only if the fit's mr_inference_args is NULL. A default fit stores list() (formal at forestsearch_main.R:1724, captured by mget at 1885), so this site is reached. An integer overflow here (seed > integer.max - 900000) is not caught locally; it propagates to the caller's tryCatch.

### R/fs_oc_grid.R:311 — fs_oc_invert

- **site:** R/fs_oc_grid.R:311
- **fun:** fs_oc_invert
- **seed_source:** Comes from fs_oc_invert's own formal `seed` (default NULL, fs_oc_grid.R:270), via `if (!is.null(seed)) set.seed(seed)`. The user sets it directly, and NULL means no set.seed. No internal caller.
- **classification:** argument
- **restored:** no
- **exported_reaching:** fs_oc_invert only (confirmed; nothing in R/ calls fs_oc_invert).
- **notes:** fs_oc_family_enumerate (fs_oc_grid.R:305) runs before the set.seed, unless `family` is supplied. One .fs_oc_draw follows at line 312, and the bisection reuses those fixed draws. The seed is recorded in settings$seed (line 316).

### R/fs_oc_grid.R:492 — .fs_oc_sweep (internal)

- **site:** R/fs_oc_grid.R:492
- **fun:** .fs_oc_sweep (internal)
- **seed_source:** Comes from the formal `seed` (no default, fs_oc_grid.R:489). The sole caller, fs_oc_grid.R:152, forwards fs_oc_grid's formal `seed = NULL` (fs_oc_grid.R:100) unchanged. The user can set it from fs_oc_grid, and NULL means no set.seed.
- **classification:** argument
- **restored:** no
- **exported_reaching:** fs_oc_grid only (the sole call is fs_oc_grid.R:152).
- **notes:** fs_oc_grid reseeds with the same seed for every (n, gate) pair (loops at fs_oc_grid.R:140-153), which gives common random numbers. In blocked mode (lines 525-546), successive blocks continue one stream from that single set.seed.

### R/fs_oc_predict.R:234 — fs_oc_predict

- **site:** R/fs_oc_predict.R:234
- **fun:** fs_oc_predict
- **seed_source:** Comes from fs_oc_predict's formal `seed` (default NULL, fs_oc_predict.R:176), via `if (!is.null(seed)) set.seed(seed)`. The user sets it directly. No internal caller.
- **classification:** argument
- **restored:** no
- **exported_reaching:** fs_oc_predict only (confirmed).
- **notes:** fs_oc_family_enumerate (fs_oc_predict.R:213) runs before the set.seed. The value is echoed into the result via .fs_oc_result (line 241).

### R/generate_aft_dgm_main.R:248 — generate_aft_dgm_flex

- **site:** R/generate_aft_dgm_main.R:248
- **fun:** generate_aft_dgm_flex
- **seed_source:** Comes from the formal `seed` (default 8316951, generate_aft_dgm_main.R:238) and is set unconditionally, so seed = NULL reseeds randomly. A direct call is an argument. inject_mrct_structure passes its own `seed` formal (default 8316951L, mrct_inject_dgm.R:386) at mrct_inject_dgm.R:525, so it is an argument there. mrct_build_dgm reaches the same parameter through ... -> inject_mrct_structure (mrct_minor_helpers.R:125), also an argument. In calibrate_k_treat it is base_args$seed, via modifyList and do.call (calibrate_k_treat.R:134-135): absent means the formal default 8316951, and an explicit NULL means a random reseed at each evaluation, which the warning at 119-123 flags. In sensitivity_analysis_k_inter it is an argument through ... (find_k_inter_main.R:262-266). From create_dgm_for_mrct it is fixed at 8316951 (no seed formal and no seed passed, mrct_simulation.R:1162-1170, 1222-1245). From find_k_inter_for_target_hr it is also fixed at 8316951 (no seed formal or ..., find_k_inter_main.R:88-103; calls at 106-121 and 139-153).
- **classification:** argument
- **restored:** no
- **exported_reaching:** generate_aft_dgm_flex, inject_mrct_structure, mrct_build_dgm, calibrate_k_treat and sensitivity_analysis_k_inter (user-settable); create_dgm_for_mrct and find_k_inter_for_target_hr (fixed 8316951). The static list is confirmed and nothing is missing (create_gbsg_dgm only mentions it in a .Deprecated message, sim_aft_gbsg.R:212).
- **notes:** In calibrate_k_treat and find_k_inter_for_target_hr the build repeats at every root-finding evaluation and reseeds to the same value each time, which keeps the objective smooth. Every call leaves the global stream positioned after the seed's draws.

### R/generate_glm_dgm.R:402 — generate_glm_dgm

- **site:** R/generate_glm_dgm.R:402
- **fun:** generate_glm_dgm
- **seed_source:** Comes from the formal `seed` (default 8316951L, generate_glm_dgm.R:195) and is set unconditionally. A direct call is an argument. calibrate_glm_interaction passes its own formal `seed` (default 8316951L, calibrate_glm_interaction.R:107) at calibrate_glm_interaction.R:165 (root objective) and 199 (final build), so it is an argument there too. calibrate_glm_interaction has no ... argument.
- **classification:** argument
- **restored:** no. The save at generate_glm_dgm.R:416-417 comes after this set.seed and after the sample() at 403, and only when k_random_noise > 0. The caller's stream is replaced.
- **exported_reaching:** generate_glm_dgm; calibrate_glm_interaction (confirmed; no other caller in R/).
- **notes:** Bootstrap super-population sampling happens at line 403. The same seed is reset at every root evaluation inside calibrate_glm_interaction.

### R/generate_glm_dgm.R:418 — generate_glm_dgm (k_random_noise block)

- **site:** R/generate_glm_dgm.R:418
- **fun:** generate_glm_dgm (k_random_noise block)
- **seed_source:** Comes from the formal `noise_seed` (default 20260807L, generate_glm_dgm.R:198). The site runs only when k_random_noise > 0 (default 0L, line 197; branch at 414). Only a direct generate_glm_dgm() call can set either value. calibrate_glm_interaction passes neither (calibrate_glm_interaction.R:149-167, 183-201) and has no ..., so the branch never runs from there.
- **classification:** argument
- **restored:** yes. Explicit save at generate_glm_dgm.R:416-417 and assign/rm at 421-422 (not on.exit). The restored state is the stream after the sample() at line 403, so the noise draws leave the downstream stream unchanged.
- **exported_reaching:** generate_glm_dgm only. calibrate_glm_interaction is removed because it can never set k_random_noise > 0.
- **notes:** The `else rm(.Random.seed)` branch (line 422) is dead in practice, because set.seed(seed) at 402 always creates .Random.seed. Because the restore is not done with on.exit, an error in the rnorm loop (419-420) would skip it. setup_gbsg_dgm.R:126-130 uses the same scheme.

### R/get_FSdata_helpers.R:259 — lasso_selection

- **site:** R/get_FSdata_helpers.R:259
- **fun:** lasso_selection
- **seed_source:** Comes from the formal `seedit` (default 8316951, get_FSdata_helpers.R:242). The only internal caller, get_FSdata (get_fsdata.R:328-335, inside if(use_lasso) at 326), does not pass seedit, and get_FSdata has no seed formal (get_fsdata.R:84-90). forestsearch calls get_FSdata through filter_call_args(args_call_all, get_FSdata, ...) (forestsearch_main.R:3324-3333), which drops seedit. So from get_FSdata, forestsearch and every wrapper the value is the literal 8316951; per-replicate seedit values in fpr_calibration and run_null_calibration have no effect here. It is an argument only in a direct lasso_selection(seedit=) call, and NULL skips seeding (line 248).
- **classification:** fixed
- **restored:** yes: the old state is captured at get_FSdata_helpers.R:249-251 and restored, or removed if it did not exist, by the on.exit at 252-258, registered before the set.seed at 259.
- **exported_reaching:** lasso_selection (direct, argument); get_FSdata (use_lasso = TRUE); forestsearch (use_lasso = TRUE, default FALSE at forestsearch_main.R:1642, and only with mr_inference = FALSE, because the consistency-path MR guard at forestsearch_main.R:1898-1901 / forestsearch_helpers.R:2124-2129 stops the call first and the dina/grf branches return at 2887/3072 before get_FSdata); run_simulation_analysis (legacy FS arms force use_lasso = TRUE, run_simulation_analysis.R:53, 588, 613); compare_selection_rules, fpr_calibration, mrct_region_sims and run_null_calibration (through fs args); mrct_case_study and mrct_trial_analysis (only with the non-default mr_inference = FALSE); bootstrap_results, forestsearch_bootstrap_dofuture, forestsearch_Kfold and forestsearch_tenfold (stored use_lasso = TRUE with no replicate MR). Removed: fs_fdr_report (forces mr_inference = TRUE at fs_fdr_report.R:239, so a use_lasso fit errors at the guard), fs_family_report (use_lasso = FALSE hard-coded, fs_family_report.R:214), fs_oc_family_enumerate (use_lasso = FALSE, fs_oc_family.R:316), fs_oc_grid, fs_oc_invert and fs_oc_predict (they reach get_FSdata only through fs_oc_family_enumerate), and fpr_calibration_transfer (fs_params = list(), fpr_calibration.R:548). Main path: forestsearch(use_lasso = TRUE) -> get_FSdata -> lasso_selection.
- **notes:** Same pattern as subgroup.consistency: from every exported entry point except lasso_selection itself, the LASSO cross-validation fold seed is the literal 8316951, whatever seedit the user passes to forestsearch. The caller's RNG stream is preserved.

### R/grf_subg_harm_glm.R:460 — grf.subg.harm.glm

- **site:** R/grf_subg_harm_glm.R:460
- **fun:** grf.subg.harm.glm
- **seed_source:** Comes from the formal `seedit` (default 8316951L, grf_subg_harm_glm.R:302), set unconditionally and also passed to grf::causal_forest(seed = seedit) at 463/467. A direct call is an argument. From forestsearch it is the seedit formal (default 8316951, forestsearch_main.R:1640), passed to .build_grf_glm_args(seedit) (grf_args.R:63) and then do.call at forestsearch_helpers.R:1904 (subgroup_method = "grf", via .forestsearch_grf_select at forestsearch_main.R:2954-2963, only when no grf_res is supplied) and at forestsearch_main.R:3122 (use_grf screening). Both apply only when outcome_type != "survival". In run_simulation_analysis the standalone GRF arm uses grf_params$seedit, else the literal 8316951L from default_grf_params_gen (run_simulation_analysis.R:86), merged at 552 and filtered at 1236-1241; the FS arms use fs_params$seedit. Elsewhere: fpr_calibration uses seed_i (argument-derived); run_null_calibration uses seed_s (argument-derived); bootstrap and CV use fs.est$args_call_all$seedit (stored in object; grf_res is nulled per replicate at bootstrap_analysis_dofuture.R:573 and in CV); compare_selection_rules passes it through ....
- **classification:** argument
- **restored:** no
- **exported_reaching:** grf.subg.harm.glm; forestsearch (GLM outcome_type with subgroup_method = "grf", or use_grf = TRUE with mr_inference = FALSE); run_simulation_analysis (GLM outcome_type with run_grf, run_fs_grf or methods); compare_selection_rules; fpr_calibration; run_null_calibration; bootstrap_results, forestsearch_bootstrap_dofuture, forestsearch_Kfold, forestsearch_tenfold (GLM + GRF original fit); mrct_trial_analysis (fs_args with a GLM outcome_type plus subgroup_method = "grf", or use_grf with mr_inference = FALSE); mrct_region_sims and mrct_case_study only technically, by forcing a GLM outcome_type onto the survival y_sim through fs_args. Removed: fs_family_report, fs_oc_family_enumerate, fs_oc_grid, fs_oc_invert, fs_oc_predict (no forestsearch call, and get_FSdata never fits GRF), fs_fdr_report (stops for a GLM outcome_type, fs_fdr_report.R:207-212), fpr_calibration_transfer (fs_params = list(), so survival with no GRF). Main path: forestsearch(outcome_type = binary/continuous/count, subgroup_method = "grf") -> .forestsearch_grf_select -> grf.subg.harm.glm.
- **notes:** The set.seed seeds R's global stream, while the forest's own RNG comes from causal_forest(seed = seedit). The two defaults differ only in type: 8316951 (double) in forestsearch, 8316951L in grf.subg.harm.glm. The function does no on.exit or .Random.seed handling.

### R/guohe_adaptive_r.R:205 — guohe_adaptive_r

- **site:** R/guohe_adaptive_r.R:205
- **fun:** guohe_adaptive_r
- **seed_source:** Formal `seed = NULL` (R/guohe_adaptive_r.R:190), guarded by `if (!is.null(seed))`. guohe_adaptive_r is exported (NAMESPACE:226) and no function in R/ calls it, so only a direct caller sets the value, and any value can be passed. With the default NULL there is no set.seed.
- **classification:** argument
- **restored:** no. R/guohe_adaptive_r.R has no save/restore of .Random.seed and no on.exit. The fold draw `sample(rep_len(seq_len(v), n))` (:215) and the fold-level guohe_algorithm3 bootstrap indices (:251-259, no seed passed) consume this stream, and the caller's state is left changed.
- **exported_reaching:** guohe_adaptive_r (direct only; no internal callers in R/)
- **notes:** With a non-NULL seed, the same value is applied again in the full-data refit `guohe_algorithm3(..., seed = seed, ...)` (R/guohe_adaptive_r.R:283-288, seed on :286) when refit = TRUE. The refit's bootstrap indices therefore match a standalone guohe_algorithm3(seed = seed) call.

### R/guohe_algorithm3.R:395 — guohe_algorithm3

- **site:** R/guohe_algorithm3.R:395
- **fun:** guohe_algorithm3
- **seed_source:** Formal `seed = NULL` (R/guohe_algorithm3.R:243), guarded by `if (!is.null(seed))`. (1) Direct call (exported, NAMESPACE:227): set by the user. (2) fs_to_guohe: its formal `seed = NULL` (R/fs_to_guohe.R:340) is passed through as `seed = seed` (:360), so the user sets it. (3) guohe_adaptive_r: the fold fits (R/guohe_adaptive_r.R:251-259) pass no seed, so it is NULL and this set.seed is skipped. The refit (:283-288) passes `seed = seed` (:286), which comes from guohe_adaptive_r's formal (:190).
- **classification:** argument
- **restored:** no. R/guohe_algorithm3.R has no save/restore and no on.exit. The B index vectors are drawn serially right after (:396), and the caller's RNG is left advanced from seed.
- **exported_reaching:** guohe_algorithm3 (direct); fs_to_guohe (R/fs_to_guohe.R:357-362, seed at :360); guohe_adaptive_r (refit only, R/guohe_adaptive_r.R:283-288; the fold fits at :251-259 never run the set.seed)
- **notes:** The default NULL means no seeding. All indices are drawn before the fitting step, and future_lapply runs with `future.seed = FALSE` (:412-417, flag at :414), so parallel = TRUE reproduces the serial result for the same seed.

### R/mrct_inject_dgm.R:97 — expand_covariates

- **site:** R/mrct_inject_dgm.R:97
- **fun:** expand_covariates
- **seed_source:** Formal `seed = 20260903L` (R/mrct_inject_dgm.R:90), set unconditionally by `set.seed(seed)`. (1) Direct expand_covariates() call (NAMESPACE:143): the user sets it, default 20260903L. (2) inject_mrct_structure always passes `seed = seed` (:466) from its own formal `seed = 8316951L` (:386), so on that path the 20260903L default is never used. (3) mrct_build_dgm forwards `...` to inject_mrct_structure (R/mrct_minor_helpers.R:120, :125), so mrct_build_dgm(data, ..., seed = X) sets it.
- **classification:** argument
- **restored:** no. There is no save/restore and no on.exit. The bootstrap branch (`sample.int` at :100) and the copula branch (`mvtnorm::rmvnorm` at :143) both leave the caller's RNG advanced.
- **exported_reaching:** expand_covariates (direct); inject_mrct_structure (:466); mrct_build_dgm (via ... at R/mrct_minor_helpers.R:125)
- **notes:** The effective default depends on the entry point: 20260903L for a direct call, 8316951L through inject_mrct_structure or mrct_build_dgm. The set.seed at :97 runs before the latentcor fit (:120), so any RNG use inside latentcor is also covered by this seed.

### R/mrct_inject_dgm.R:277 — make_region_model -> nested closure `draw` (R/mrct_inject_dgm.R:276-279), returned as list element $draw

- **site:** R/mrct_inject_dgm.R:277
- **fun:** make_region_model -> nested closure `draw` (R/mrct_inject_dgm.R:276-279), returned as list element $draw
- **seed_source:** `seed` is a formal of the closure `draw = function(df, seed = NULL)` (:276), guarded by `if (!is.null(seed))`. make_region_model never calls draw. On the main path, inject_mrct_structure calls `rm_$draw(df_pop, seed = seed + 1L)` (:478) and `rm_$draw(seed_df, seed = seed + 2L)` (:479), where seed is inject_mrct_structure's formal `seed = 8316951L` (:386). The user sets that formal directly or through mrct_build_dgm's ... (R/mrct_minor_helpers.R:125). A user who calls make_region_model() directly (NAMESPACE:236) can call the returned $draw(df, seed = X) with any value, or with NULL for no seeding. The closure is also stored as dgm$mrct$region_model (:595). Nothing in R/ calls it again.
- **classification:** argument-derived
- **restored:** no. The closure does not save or restore the state, and inject_mrct_structure does not wrap the two draws.
- **exported_reaching:** inject_mrct_structure (:478-479, seed+1L / seed+2L); mrct_build_dgm (via ... at R/mrct_minor_helpers.R:125); make_region_model (only when the user calls the returned $draw closure; make_region_model itself never runs this set.seed). Not reached in observed-region mode (`if (observed)` at :470, `rm_ <- NULL` at :472).
- **notes:** The static row's 'seed = local' refers to the closure's own formal (default NULL). On the inject path the value is seed + 1L / seed + 2L (argument-derived). For a direct $draw() caller it is a plain argument. The seed-row guard `sample.int` at :481-482 draws from the stream left by the seed+2L draw.

### R/mrct_inject_dgm.R:572 — inject_mrct_structure

- **site:** R/mrct_inject_dgm.R:572
- **fun:** inject_mrct_structure
- **seed_source:** `seed + 3L`, where seed is the formal `seed = 8316951L` (R/mrct_inject_dgm.R:386). It is set unconditionally, and seed is not reassigned in the body. A direct inject_mrct_structure(seed = X) call (NAMESPACE:229) changes it, and mrct_build_dgm(data, ..., seed = X) forwards it via ... (R/mrct_minor_helpers.R:120, :125).
- **classification:** argument-derived
- **restored:** no. inject_mrct_structure has no save/restore and no on.exit. calculate_hazard_ratios (:573) draws `rexp` (R/generate_aft_dgm_helpers.R:1080) from this stream, and the caller's state is left changed.
- **exported_reaching:** inject_mrct_structure (direct); mrct_build_dgm (R/mrct_minor_helpers.R:125)
- **notes:** The same seed family is used at seed (expand_covariates :466; generate_aft_dgm_flex :525) and at seed+1L / seed+2L (region draws :478-479). Seeding immediately before calculate_hazard_ratios makes dgm$hazard_ratios reproducible from seed.

### R/mrct_inject_dgm.R:664 — mrct_truth_by_region

- **site:** R/mrct_inject_dgm.R:664
- **fun:** mrct_truth_by_region
- **seed_source:** Formal `seed = 1L` (R/mrct_inject_dgm.R:659), set unconditionally before the common error draw `rexp` (:665). It is exported (NAMESPACE:249), and nothing in R/ calls it, so a direct caller sets the value.
- **classification:** argument
- **restored:** no. There is no save/restore and no on.exit. The caller's RNG is left at set.seed(seed) followed by nrow(ds) rexp draws.
- **exported_reaching:** mrct_truth_by_region (direct only)
- **notes:** With the default seed = 1L, every call resets the global stream to set.seed(1) without warning.

### R/run_simulation_analysis.R:509 — run_simulation_analysis

- **site:** R/run_simulation_analysis.R:509
- **fun:** run_simulation_analysis
- **seed_source:** `seed_base + 1000L * sim_id`. seed_base is the formal `seed_base = 8316951L` (R/run_simulation_analysis.R:413), and sim_id is a formal with no default (:390). Both are set by the direct caller (exported, NAMESPACE:278), and nothing in R/ calls run_simulation_analysis. The call runs only `if (n_add_noise > 0L)` (:508), and the formal is `n_add_noise = 0L` (:399), so it is off by default.
- **classification:** argument-derived
- **restored:** no. There is no save/restore and no on.exit. The rnorm noise draws (:511) leave the caller's RNG advanced.
- **exported_reaching:** run_simulation_analysis (direct only)
- **notes:** The same replicate's data draw uses seed_base + sim_id (:476, :482). The noise seed for sim_id = k equals the data seed for sim_id = 1000k, so the two collide once n_sims >= 1000. The roxygen (:266-267) documents only 'replicate seed = seed_base + sim_id', and the n_add_noise doc (:221-222) does not mention the 1000L * sim_id noise offset.

### R/run_subgroup_sims.R:512 — .run_subgroup_sims_one (internal; called inside the foreach %dofuture% body of run_subgroup_sims)

- **site:** R/run_subgroup_sims.R:512
- **fun:** .run_subgroup_sims_one (internal; called inside the foreach %dofuture% body of run_subgroup_sims)
- **seed_source:** `seed_base + ss + rand_seed_offset`, all formals of the internal function (R/run_subgroup_sims.R:494-495). run_subgroup_sims passes them at :894-897: seed_base from its formal `seed_base = 0L` (:675), rand_seed_offset from its formal `rand_seed_offset = 1e6L` (:676), and ss from the foreach index `ss = seq_len(n_sims)` (:888), where n_sims is a user formal. None is reassigned in :683-893. The call runs only `if (!is.null(benchmarks))` (:511), and the default `benchmarks = benchmark_spec()` (:672) is non-NULL, so it runs by default.
- **classification:** argument-derived
- **restored:** no (package code). Neither function saves or restores .Random.seed. The on.exit calls at :818, :830 and :881 restore only the future plan and the future.rng.onMisuse option. Under multisession the set.seed runs in worker processes. Under plan(sequential), a re-run check (future 1.80.0) confirmed the caller's .Random.seed is identical after the %dofuture% loop. future provides that restore, not this package.
- **exported_reaching:** run_subgroup_sims (only caller, R/run_subgroup_sims.R:894)
- **notes:** The same iteration's DGM draw uses seed_base + ss (:503). The benchmark seed of trial ss equals the draw seed of trial ss + rand_seed_offset, so the two collide only if n_sims > 1e6 at the default offset. The defaults (0L, 1e6L) reproduce the vignettes (roxygen :625-626).

### R/setup_gbsg_dgm.R:130 — setup_gbsg_dgm

- **site:** R/setup_gbsg_dgm.R:130
- **fun:** setup_gbsg_dgm
- **seed_source:** Formal `noise_seed = 20260807L` (R/setup_gbsg_dgm.R:100). The call runs only `if (k_random_noise > 0)` (:126), and the formal is `k_random_noise = 0L` (:99). A direct setup_gbsg_dgm(k_random_noise > 0, noise_seed = X) call (NAMESPACE:286) sets it. fs_fdr_report calls setup_gbsg_dgm at R/fs_fdr_report.R:256 and :262 without k_random_noise or noise_seed, so this set.seed never runs from fs_fdr_report. Even if it did, the value would be fixed at 20260807L.
- **classification:** argument
- **restored:** yes, by explicit save/restore, not on.exit. It saves with `has_seed <- exists(".Random.seed", ...)` and `old_seed <- get(...)` (:128-129), and after the rnorm loop (:131-132) restores with `assign(".Random.seed", old_seed, ...)` (:133) or `rm(...)` (:134).
- **exported_reaching:** setup_gbsg_dgm (direct, only when k_random_noise > 0). fs_fdr_report is in the static list but never runs this set.seed (k_random_noise stays 0L at R/fs_fdr_report.R:256, :262).
- **notes:** The saved state is not the caller's pre-call state. .create_gbsg_dgm_ (called at :106-117) has already run set.seed(seed) and set.seed(seed + 1) (R/sim_aft_gbsg.R:478, :604) without restoring. As a result .Random.seed always exists at :128, and the rm() branch at :134 cannot run from here. The restore is not in on.exit, so an error or interrupt in the loop would leave the RNG on the noise_seed stream. Afterwards the caller's stream matches a k_random_noise = 0 build, as the comment at :121-124 states. The value is recorded as dgm$noise_seed (:187).

### R/sim_aft_gbsg.R:478 — .create_gbsg_dgm_ (internal)

- **site:** R/sim_aft_gbsg.R:478
- **fun:** .create_gbsg_dgm_ (internal)
- **seed_source:** Formal `seed = SEED_BASE` (R/sim_aft_gbsg.R:245), where SEED_BASE <- 8316951L is a package constant (:29). It is set unconditionally and not reassigned. The user can set it through: setup_gbsg_dgm's formal `seed = 8316951L` (R/setup_gbsg_dgm.R:97), passed at :115; create_gbsg_dgm's formal `seed = SEED_BASE` (R/sim_aft_gbsg.R:208), passed at :229; calibrate_k_inter's `...` (:1014), forwarded at :1024-1027; validate_k_inter_effect's `...` (:1232), forwarded at :1254. fs_fdr_report passes no seed: setup_gbsg_dgm at R/fs_fdr_report.R:256 and :262, and calibrate_k_inter at :258-259 with no `...`. From fs_fdr_report the value is therefore fixed at 8316951L. Its seed_base (:193) feeds only the per-replicate seeds (:287).
- **classification:** argument
- **restored:** no. Neither .create_gbsg_dgm_ nor any wrapper saves or restores the state. setup_gbsg_dgm's save/restore (:128-134) covers only its later noise draw.
- **exported_reaching:** setup_gbsg_dgm (:115); create_gbsg_dgm (deprecated, :229); calibrate_k_inter (via ..., :1026); validate_k_inter_effect (via ..., :1254); fs_fdr_report (fixed 8316951L, via setup_gbsg_dgm and calibrate_k_inter with no seed). get_dgm_with_output (:1108) is internal, not exported, and has no callers.
- **notes:** The value is an argument from setup_gbsg_dgm, create_gbsg_dgm, calibrate_k_inter and validate_k_inter_effect, and fixed (8316951L) from fs_fdr_report. It is stored as `seed = seed` (:732) and read later by the legacy simulator as dgm$seed + sim_id (:856, :884), + 1000 (:900) and + 2000 (:908). The calibrate_k_inter root search (eff_at, :1023-1029) re-seeds identically on every evaluation.

### R/sim_aft_gbsg.R:604 — .create_gbsg_dgm_ (internal)

- **site:** R/sim_aft_gbsg.R:604
- **fun:** .create_gbsg_dgm_ (internal)
- **seed_source:** `seed + 1`, where seed is the same formal `seed = SEED_BASE` (= 8316951L; R/sim_aft_gbsg.R:245, :29). It is set unconditionally. The user can set it through the same pass-throughs: setup_gbsg_dgm (R/setup_gbsg_dgm.R:97, :115), create_gbsg_dgm (R/sim_aft_gbsg.R:208, :229), calibrate_k_inter ... (:1026) and validate_k_inter_effect ... (:1254). From fs_fdr_report it is fixed at 8316952 (R/fs_fdr_report.R:256, :258-259 and :262 pass no seed).
- **classification:** argument-derived
- **restored:** no. There is no save/restore and no on.exit. The rexp draw (:605) leaves the caller's RNG on the seed+1 stream (setup_gbsg_dgm's later save/restore preserves this state rather than undoing it).
- **exported_reaching:** setup_gbsg_dgm; create_gbsg_dgm (deprecated); calibrate_k_inter (via ...); validate_k_inter_effect (via ...); fs_fdr_report (fixed 8316952, no pass-through)
- **notes:** The value is argument-derived from setup_gbsg_dgm, create_gbsg_dgm, calibrate_k_inter and validate_k_inter_effect, and fixed (8316952) from fs_fdr_report. `seed + 1` is a double. set.seed coerces it to integer, so the stream is unaffected.

### R/sim_aft_gbsg.R:856 — .simulate_from_gbsg_dgm_ (n non-NULL, draw_treatment = TRUE branch)

- **site:** R/sim_aft_gbsg.R:856
- **fun:** .simulate_from_gbsg_dgm_ (n non-NULL, draw_treatment = TRUE branch)
- **seed_source:** dgm$seed + sim_id. dgm$seed is stored at R/sim_aft_gbsg.R:732 from the seed formal of .create_gbsg_dgm_ (:245, default SEED_BASE = 8316951L, :29). It gets there through create_gbsg_dgm(seed = SEED_BASE) (exported, deprecated; formal :208, passed :229) or setup_gbsg_dgm(seed = 8316951L) (exported; R/setup_gbsg_dgm.R:97, passed :115; class c(aft_dgm_flex, gbsg_dgm, list) at :182, so it passes the inherits check at :836). The internal get_dgm_with_output (:1108, passes seed through ...) is not exported. sim_id is the formal of .simulate_from_gbsg_dgm_ (:827, default 1), passed from simulate_from_gbsg_dgm(sim_id = 1) (:793 -> :814). The user can change both: seed when building the DGM, sim_id at simulation time.
- **classification:** stored in object (from an argument)
- **restored:** no (R/sim_aft_gbsg.R has no on.exit or .Random.seed handling, so the caller's RNG stream is overwritten)
- **exported_reaching:** simulate_from_gbsg_dgm (only caller: R/sim_aft_gbsg.R:813, deprecated wrapper). The dgm comes from create_gbsg_dgm or setup_gbsg_dgm. Static list confirmed.
- **notes:** Runs only when n is non-NULL and draw_treatment = TRUE. Offsets overlap across replicates: seed+sim_id here equals seed+(sim_id-1000)+1000 at :900, so replicate k's sampling stream is replicate k-1000's event-time stream (likewise :908 for k-2000). With sim_id = 0 the seed equals the super-population sampling seed at :478. With the default sim_id = 1 it equals the seed+1 used at :604 for the super-population's observed event times.

### R/sim_aft_gbsg.R:884 — .simulate_from_gbsg_dgm_ (n non-NULL, draw_treatment = FALSE branch)

- **site:** R/sim_aft_gbsg.R:884
- **fun:** .simulate_from_gbsg_dgm_ (n non-NULL, draw_treatment = FALSE branch)
- **seed_source:** Same as :856: dgm$seed (stored at R/sim_aft_gbsg.R:732 from create_gbsg_dgm(seed=) :208 -> :229 or setup_gbsg_dgm(seed=) R/setup_gbsg_dgm.R:97 -> :115; default 8316951L) + sim_id (formal :827, default 1, from simulate_from_gbsg_dgm(sim_id=) :793 -> :814). The user can change it.
- **classification:** stored in object (from an argument)
- **restored:** no
- **exported_reaching:** simulate_from_gbsg_dgm (via R/sim_aft_gbsg.R:813). Static list confirmed.
- **notes:** Mutually exclusive with :856 and uses the same seed value, so both branches start from the same stream. Not reached when n = NULL (the whole super-population is used). Same :478/:604 collisions as :856 at sim_id = 0 and 1.

### R/sim_aft_gbsg.R:900 — .simulate_from_gbsg_dgm_ (event-time draw)

- **site:** R/sim_aft_gbsg.R:900
- **fun:** .simulate_from_gbsg_dgm_ (event-time draw)
- **seed_source:** dgm$seed + sim_id + 1000. The +1000 offset is a literal. dgm$seed comes from the stored object (R/sim_aft_gbsg.R:732, set through the create_gbsg_dgm or setup_gbsg_dgm seed formal, default 8316951L), and sim_id from simulate_from_gbsg_dgm(sim_id = 1) :793 -> :814. The user can change both.
- **classification:** stored in object (from an argument)
- **restored:** no
- **exported_reaching:** simulate_from_gbsg_dgm. Static list confirmed.
- **notes:** Always executed. The +1000 offset collides with the sampling seed of replicate sim_id+1000 and with the censoring seed of replicate sim_id-1000, so runs of more than 1000 replicates reuse streams across replicates.

### R/sim_aft_gbsg.R:908 — .simulate_from_gbsg_dgm_ (censoring draw)

- **site:** R/sim_aft_gbsg.R:908
- **fun:** .simulate_from_gbsg_dgm_ (censoring draw)
- **seed_source:** dgm$seed + sim_id + 2000. The +2000 offset is a literal. dgm$seed comes from the object (R/sim_aft_gbsg.R:732, from the create_gbsg_dgm or setup_gbsg_dgm seed formal, default 8316951L), and sim_id from simulate_from_gbsg_dgm (:793 -> :814). The user can change both.
- **classification:** stored in object (from an argument)
- **restored:** no
- **exported_reaching:** simulate_from_gbsg_dgm. Static list confirmed.
- **notes:** Always executed: the set.seed comes before the Weibull/uniform if. Same cross-replicate overlap as :900 (offsets 2000, 1000 and 0).

### R/simulate_from_dgm.R:208 — simulate_from_dgm

- **site:** R/simulate_from_dgm.R:208
- **fun:** simulate_from_dgm
- **seed_source:** Formal seed = NULL (R/simulate_from_dgm.R:180); set.seed runs only when it is non-NULL. Exported callers: calibrate_cens_adjust(seed = 42) (:795), passed verbatim at :859 (sim_metric) and :893; fs_build_eval_frame(eval_seed = 20260628L) (R/betaHhat_truth.R:823), passed at :851; fs_fdr_report(seed_base = 8316951L) (R/fs_fdr_report.R:193) -> seed_base + 100000L*j (:287) + i (:409), passed at :270; mrct_case_study(seed, k) (required formals, R/mrct_case_study.R:241) -> seed + k (:279), do.call at :281; mrct_region_sims(seed = NULL) (R/mrct_simulation.R:280) -> seed + sim (:428), do.call at :431 (a NULL value removes the element in modifyList, so there is no seeding, and any sim_args$seed is discarded too); run_simulation_analysis(seed_base = 8316951L, sim_id) (R/run_simulation_analysis.R:413) -> seed_base + sim_id (:482); run_subgroup_sims(seed_base = 0L) (R/run_subgroup_sims.R:675) -> seed_base + ss (:503, sim_fun set at :773/:787). The user can change it on every path.
- **classification:** argument
- **restored:** no in simulate_from_dgm itself. Among the callers only mrct_case_study restores: it saves the state at R/mrct_case_study.R:263-265 and its on.exit (:266-271) resets RNGkind and then reassigns or rm()s .Random.seed. It forces L'Ecuyer-CMRG at :280 and resets the kind at :282. The other callers do not restore; fs_fdr_report, mrct_region_sims and run_subgroup_sims draw inside %dofuture% futures.
- **exported_reaching:** simulate_from_dgm (direct); calibrate_cens_adjust; fs_build_eval_frame; fs_fdr_report; mrct_case_study; mrct_region_sims; run_simulation_analysis; run_subgroup_sims. Static list confirmed. Main path: mrct_region_sims or run_simulation_analysis, with a per-replicate seed of seed/seed_base plus the replicate index.
- **notes:** An argument on direct calls, and argument-derived from every wrapper except calibrate_cens_adjust and fs_build_eval_frame, which pass their formal verbatim. calibrate_cens_adjust deliberately reuses one seed for every root-finding evaluation (common random numbers, comment at :849). Under mrct_region_sims the futures use seed = TRUE (:408), so the worker RNG kind is L'Ecuyer-CMRG and set.seed(seed+sim) seeds that kind. In run_simulation_analysis the simulation seed seed_base+sim_id (:482) overlaps the noise seed seed_base+1000*sim_id (:509) of another replicate (sim_id 1000 vs 1).

### R/simulate_from_glm_dgm.R:97 — simulate_from_glm_dgm

- **site:** R/simulate_from_glm_dgm.R:97
- **fun:** simulate_from_glm_dgm
- **seed_source:** Formal seed = NULL (R/simulate_from_glm_dgm.R:72); set.seed runs only when it is non-NULL. Callers: fs_build_eval_frame(eval_seed = 20260628L) (R/betaHhat_truth.R:823), passed at :855 (binary outcome only); fs_dgm_feasibility(seed = NULL) (R/fs_dgm_feasibility.R:123) -> s_i = as.integer(seed) + i (:180), passed at :182-183 (also set.seed(s_i) at :181); run_simulation_analysis(seed_base = 8316951L) -> seed_base + sim_id (R/run_simulation_analysis.R:476); run_subgroup_sims(seed_base = 0L) -> seed_base + ss (R/run_subgroup_sims.R:503, sim_fun set at :747/:757). The user can change it on every path.
- **classification:** argument
- **restored:** no in simulate_from_glm_dgm. fs_dgm_feasibility restores partly: its on.exit (R/fs_dgm_feasibility.R:162-166) reassigns the old .Random.seed only if one existed and does not rm a newly created one. The other callers do not restore.
- **exported_reaching:** simulate_from_glm_dgm (direct); fs_build_eval_frame; fs_dgm_feasibility; run_simulation_analysis; run_subgroup_sims. Static list confirmed.
- **notes:** Argument-derived (seed + index) from fs_dgm_feasibility, run_simulation_analysis and run_subgroup_sims; verbatim from fs_build_eval_frame. fs_dgm_feasibility sets the same s_i twice (:181 and inside the callee), which is redundant but harmless. as.integer() truncates a non-integer seed.

### R/subgroup_consistency_main.R:58 — .make_candidate_rng_seeds

- **site:** R/subgroup_consistency_main.R:58
- **fun:** .make_candidate_rng_seeds
- **seed_source:** Formal seed (no default, :40). Its only caller is subgroup.consistency at R/subgroup_consistency_main.R:971, on the parallel batch branch (else at :878) and only when seed is non-NULL (:970). subgroup.consistency forwards its own formal seed = 8316951 (:397). From forestsearch() and every wrapper the value is fixed at 8316951: forestsearch has no seed formal and no '...' (only seedit, R/forestsearch_main.R:1640); args_call_all holds only forestsearch formals (:1884-1885); consistency_overrides (:3582-3606) adds no seed; filter_call_args (:3643-3647) followed by do.call (:3669) leaves the default. Only a direct subgroup.consistency(seed = ) call can change it.
- **classification:** argument
- **restored:** yes: the on.exit at R/subgroup_consistency_main.R:49-56 restores RNGkind(oldkind[1]) and then the .Random.seed saved at :44-47, or rm()s it when none existed. The kind was switched at :48.
- **exported_reaching:** Only callers that end up on the parallel branch: forestsearch (default parallel_args = multisession, R/forestsearch_main.R:1634-1636; fixed 8316951); compare_selection_rules (passes ... into forestsearch, R/compare_selection_rules.R:561; fixed); fpr_calibration and run_simulation_analysis, only when fs_params supplies a non-sequential parallel_args (defaults are sequential: R/fpr_calibration.R:266-268, R/run_simulation_analysis.R:68; fixed); subgroup.consistency direct, only with a non-sequential parallel_args that has plan and workers (the default list() at :401 does not reach it; seed is an argument there). Removed from the static list because they force plan = "sequential": bootstrap_results and forestsearch_bootstrap_dofuture (R/bootstrap_analysis_dofuture.R:595), forestsearch_Kfold (R/forestsearch_cross_validation.R:389), forestsearch_tenfold (:879), fpr_calibration_transfer (fs_params = list() at R/fpr_calibration.R:548), run_null_calibration (R/calibrate_null_correction.R:324-326), fs_fdr_report (R/fs_fdr_report.R:359), mrct_region_sims (R/mrct_simulation.R:537), mrct_trial_analysis and mrct_case_study (R/mrct_trial_helpers.R:41). Also removed because they only reference formals(): fs_declaration_calibration, fs_declaration_table, fs_family_report, fs_mr_inference, fs_oc_family_enumerate, fs_oc_grid, fs_oc_invert, fs_oc_predict.
- **notes:** Effectively fixed at 8316951 on every forestsearch path; it is an argument only through a direct subgroup.consistency call. Requires sg_focus != "maxeff" (fast path at :757), use_parallel (:702-742) and seed non-NULL. Bootstrap, CV and simulation replicates run sequentially, so they never use these per-candidate streams; they use the single global stream re-seeded at :498.

### R/subgroup_consistency_main.R:498 — subgroup.consistency

- **site:** R/subgroup_consistency_main.R:498
- **fun:** subgroup.consistency
- **seed_source:** Formal seed = 8316951 (R/subgroup_consistency_main.R:397); set.seed runs only when it is non-NULL (:497). A direct subgroup.consistency(seed = ) call can change it. From forestsearch() it is fixed at 8316951: forestsearch has no seed formal and no '...' (seedit, R/forestsearch_main.R:1640, is not forwarded); args_call_all holds forestsearch formals only (:1884-1885), and no .sync_args_call_all call adds a seed; consistency_overrides (:3582-3606) has no seed; filter_call_args + do.call at :3643-3647 / :3669. Every other exported caller reaches it only through do.call(forestsearch, ...), so a seed cannot be injected.
- **classification:** argument
- **restored:** no. The body has no .Random.seed or on.exit RNG handling, so it overwrites the caller's global RNG state whenever seed is non-NULL, including on the parallel, maxeff and resample paths that do not consume this stream.
- **exported_reaching:** subgroup.consistency (direct: argument); forestsearch (main path: fixed 8316951); bootstrap_results (R/bootstrap_analysis_dofuture.R:614); forestsearch_bootstrap_dofuture (via bootstrap_results, R/bootstrap_dofuture_main.R:484); compare_selection_rules (R/compare_selection_rules.R:561); forestsearch_Kfold (R/forestsearch_cross_validation.R:480); forestsearch_tenfold (:1002); fpr_calibration (R/fpr_calibration.R:308); fpr_calibration_transfer (via do.call(fpr_calibration) at :557); run_null_calibration (R/calibrate_null_correction.R:358); fs_fdr_report (R/fs_fdr_report.R:364); mrct_region_sims (R/mrct_simulation.R:575); mrct_trial_analysis (R/mrct_trial_analysis.R:120) and mrct_case_study (R/mrct_case_study.R:301), both via .mrct_analyze_trial (R/mrct_trial_helpers.R:47); run_simulation_analysis (via the internal .run_fs_analysis_gen, R/run_simulation_analysis.R:819, defined at :773). All are fixed at 8316951. Removed (formals() or string references only): fs_declaration_calibration, fs_declaration_table, fs_family_report, fs_mr_inference, fs_oc_family_enumerate, fs_oc_grid, fs_oc_invert, fs_oc_predict.
- **notes:** Fixed in practice: every forestsearch-based run, including all bootstrap, CV and simulation replicates, re-seeds to 8316951 before the consistency screen, whatever seedit or the replicate seed is. It is user-settable only by calling subgroup.consistency directly. Not reached on the early return at :461 (no candidates). This stream drives the sequential split loop (:788), which every sequential-plan wrapper uses; the parallel branch uses the :58 streams. Under the default consistency_method = "resample" (R/forestsearch_main.R:1706), the closed-form path (R/subgroup_consistency_helpers.R:1484-1490) draws no random numbers except the split fallback (:1491-1499), so the reseed mainly clobbers the caller's RNG state.

### R/synthetic_data_perturbation.R:98 — generate_bootstrap_synthetic

- **site:** R/synthetic_data_perturbation.R:98
- **fun:** generate_bootstrap_synthetic
- **seed_source:** Formal seed = 123 (R/synthetic_data_perturbation.R:57), used verbatim. Also reached from generate_gbsg_bootstrap_general(seed = 123) (:362), passed verbatim at :375 in the call at :368. The user can change it from both.
- **classification:** argument
- **restored:** no
- **exported_reaching:** generate_bootstrap_synthetic (direct); generate_gbsg_bootstrap_general (R/synthetic_data_perturbation.R:368). Static list confirmed.
- **notes:** Unconditional set.seed with a non-NULL default, so every call overwrites the caller's RNG stream. seed is not validated, and passing seed = NULL still calls set.seed(NULL), which re-initialises the RNG randomly.

### R/synthetic_data_perturbation.R:506 — generate_bootstrap_with_noise

- **site:** R/synthetic_data_perturbation.R:506
- **fun:** generate_bootstrap_with_noise
- **seed_source:** Formal seed = 123 (R/synthetic_data_perturbation.R:486), checked as numeric at :502-504 and used verbatim. The user can change it.
- **classification:** argument
- **restored:** no
- **exported_reaching:** generate_bootstrap_with_noise only (no internal callers). Static list confirmed.
- **notes:** Unconditional set.seed with a non-NULL default; NULL is rejected by the numeric check, so seeding cannot be turned off and every call overwrites the caller's RNG state.

### R/truefind_asymptotic_glm.R:338 — .detect_prob_glm_single (method = "monte_carlo" branch)

- **site:** R/truefind_asymptotic_glm.R:338
- **fun:** .detect_prob_glm_single (method = "monte_carlo" branch)
- **seed_source:** Formal seed (no default, :308-309), passed verbatim from compute_detection_probability_glm(seed = 42L) (R/truefind_asymptotic_glm.R:245) at :285 (vector theta) and :291. Only a direct compute_detection_probability_glm(seed = ) call can change it. Every other exported caller omits seed, so the value is fixed at 42L: predict_fpr_corrected (R/calibrate_null_correction.R:180-183), run_null_calibration (:373-376), fpr_calibration (R/fpr_calibration.R:249-256), fpr_calibration_transfer (:527-531, and again through do.call(fpr_calibration) at :557), suggest_thresholds (R/suggest_thresholds.R:144-148).
- **classification:** argument
- **restored:** no: none of .detect_prob_glm_single, compute_detection_probability_glm or these callers saves or restores .Random.seed. fpr_calibration's on.exit at :274 restores only the future plan.
- **exported_reaching:** compute_detection_probability_glm (direct: argument); predict_fpr_corrected, run_null_calibration, fpr_calibration, fpr_calibration_transfer, suggest_thresholds (fixed at 42L). Static list confirmed.
- **notes:** Reached only when method = "monte_carlo". A direct call defaults to method = "cubature" (:243), which never seeds. All five wrappers force method = "monte_carlo", so each P1 evaluation re-seeds the global RNG to 42; in suggest_thresholds this repeats for every (c1, c2) pair, giving common random numbers. Effectively fixed on every wrapper path.

### Verifier corrections

- R/dina.R:433 (.dina_make_folds), seed_source (c): the forestsearch seed default is resolved at R/forestsearch_helpers.R:1154-1157. The seed key is :1156, `seed = get_arg("seed", seed_default)`, and get_arg is defined at :1105. The cited :1149-1153 is the end of the select_statistic check plus a comment. Also, the use_dina fit runs only when dina_cuts is NULL (R/forestsearch_main.R:3189) and dina_res is NULL (:3206). The subgroup_method = "dina" fit runs only when dina_res is NULL (R/forestsearch_helpers.R:1532).
- R/dina.R:433, exported_reaching: drop fpr_calibration_transfer. It hard-codes fs_params = list() (R/fpr_calibration.R:548). fpr_calibration has no `...`, so use_dina or subgroup_method cannot be forwarded: an unused-argument error is caught at :557. A second fs_params errors as a duplicate name. forestsearch() defaults are use_dina = FALSE (R/forestsearch_main.R:1646) and subgroup_method = "consistency" (:1651), so the DINA fit is never requested from that entry point.
- R/dina.R:433, exported_reaching: drop fs_family_report, fs_oc_family_enumerate, fs_oc_grid, fs_oc_invert and fs_oc_predict. None of them calls forestsearch(). They only read formals(forestsearch): R/fs_family_report.R:109,123; R/fs_oc_family.R:217; R/fs_oc_grid.R:480; R/fs_oc_predict.R:181. The actual call sites of forestsearch() are R/bootstrap_analysis_dofuture.R:614, R/compare_selection_rules.R:561, R/calibrate_null_correction.R:358, R/fs_fdr_report.R:364, R/fpr_calibration.R:308, R/mrct_trial_helpers.R:47 (mrct_trial_analysis, mrct_case_study), R/forestsearch_cross_validation.R:480,1002, R/run_simulation_analysis.R:819 and R/mrct_simulation.R:575.
- R/consistency_resample.R:459, exported_reaching: the conclusion stands (no internal caller executes the set.seed). However, the claim that the other listed entry points reach the function through the closed-form calls is wrong for several of them. fs_declaration_calibration, fs_declaration_table, fs_mr_inference, fs_family_report and fs_oc_* do not reach consistency_resample() at all. They only read formals(subgroup.consistency) (R/fs_declaration_calibration.R:564, R/fs_mr_inference.R:673) or formals(forestsearch). subgroup.consistency() is called only from forestsearch() (R/forestsearch_main.R:3645).
- R/dina.R:433, seed_source (d): the dina_subgroup_bootstrap point refit, which forwards `seed`, runs only when both `sg` and `fit` are NULL (R/dina_subgroup_bootstrap.R:415, :461, :476-482). When the caller supplies either one, `seed` reaches only .dina_make_bootstrap_indices (:535).
- R/fpr_calibration.R:302: the permutation-mode stop is at :221-222 (`if (is.null(df.analysis))` / stop), not :220-221. From fpr_calibration_transfer the sim_seed is more precisely the object-stored cal_result$call_args$sim_seed plus the literal 20000L plus round(dc1[k]*1000). That last term is derived from the transfer's delta / delta_c1 argument (:552-553).
- R/forestsearch_cross_validation.R:370: note that the default Kfolds = .fs_cv_nobs(fs.est) (:290) is leave-one-out. With default arguments the set.seed(seedit) never executes. It runs only when the user passes Kfolds < n (:369).
- R/fs_dgm_feasibility.R:181: the replicate loops start at :174 (j) and :179 (i), not :173. The save at :162-166 is unconditional, not gated on seed. With seed = NULL and an existing .Random.seed, the draws consumed by simulate_from_glm_dgm are also rolled back on exit, so the next unseeded call replays the same stream.
- R/get_FSdata_helpers.R:259 (lasso_selection), exported_reaching: fs_fdr_report removed. It forces mr_inference = TRUE (fs_fdr_report.R:239). On the consistency path with use_lasso = TRUE, the MR alignment guard (forestsearch_main.R:1898-1901, which calls .validate_mr_configuration at forestsearch_helpers.R:2124-2129) stops forestsearch() before get_FSdata (forestsearch_main.R:3324). The dina and grf branches return at forestsearch_main.R:2887 and 3072, also before get_FSdata. So fs_fdr_report can never reach lasso_selection.
- R/get_FSdata_helpers.R:259, qualification: forestsearch() reaches lasso_selection only when mr_inference = FALSE (or the consistency-path guard stops the call). mrct_case_study and mrct_trial_analysis default to mr_inference = TRUE (mrct_case_study.R:245, mrct_trial_analysis.R:87), so they reach it only when the user passes mr_inference = FALSE together with fs_args$use_lasso = TRUE.
- R/generate_aft_dgm_main.R:248, calibrate_k_treat: the formal default 8316951 applies only when base_args has no seed element. An explicit base_args$seed = NULL survives modifyList and is forwarded by do.call (calibrate_k_treat.R:134-135), so set.seed(NULL) reseeds randomly at every root-finding evaluation. The warning at 119-123 covers exactly this case.
- R/grf_subg_harm_glm.R:460: the call sites are do.call(grf.subg.harm.glm, ...) at forestsearch_helpers.R:1904 and forestsearch_main.R:3122. Lines 1893 and 3115 are only the seedit = seedit argument lines. The use_grf screening route on the consistency path cannot run together with mr_inference = TRUE (alignment guard, forestsearch_helpers.R:2124-2129). mrct_trial_analysis (default mr_inference = TRUE) therefore reaches the site only through subgroup_method = "grf" or with mr_inference = FALSE. Under conservative "may call", mrct_region_sims and mrct_case_study also reach it if fs_args forces a GLM outcome_type, so they stay on the list as technically reachable.
- R/fs_mr_field_uniform.R:93, exported_reaching: the site is also reachable through forestsearch's DINA and GRF branches. .fs_apply_mr forwards field_uniform (fs_mr_inference_methods.R:180) and ci_method (line 171, default "ij"), so the site runs there when mr_inference_args sets both ci_method = "field" and field_uniform = TRUE. The seed is the same .g(mr_inference_args$seed, seedit) + 910000L (fs_mr_inference_methods.R:172).
- R/fs_mr_field_uniform.R:93, overflow note refined: when seed > .Machine$integer.max - 900000, the field block at fs_mr_inference.R:912 overflows first. It has no local tryCatch, so fs_mr_inference fails and the caller's tryCatch catches it (forestsearch_main.R:3861 or fs_mr_inference_methods.R:156). The overflow that the tryCatch at fs_mr_inference.R:1019 absorbs therefore applies only to seeds in (integer.max - 910000, integer.max - 900000].
- R/fs_mr_inference.R:615 / 912 and R/fs_mr_field_uniform.R:93, qualification for run_simulation_analysis: the legacy FS arms force use_lasso = TRUE before fs_params is merged (run_simulation_analysis.R:588, 613), and the 'fs-*' method labels set use_grf = TRUE and use_dina = TRUE (run_simulation_analysis.R:133-134). MR runs on those arms only if fs_params also sets those front ends to FALSE. Otherwise the alignment guard stops forestsearch, and the arm's tryCatch (run_simulation_analysis.R:818-821) records it as a failure. The dina-* and grf-* labels are not affected. The site is still reachable, so the function stays on the list.
- R/generate_glm_dgm.R:418, restored wording: the saved and restored state is the one after the sample() at line 403, not the state right after line 402.
- No classification, restored, or exported-caller finding changes. All 11 rows were re-derived from the start_tree source and the substance holds. Only line citations are tightened and one note is added.
- R/guohe_adaptive_r.R:205 (and the :395 row): the fold-level guohe_algorithm3 call runs from R/guohe_adaptive_r.R:251 to :259 (`), silent = TRUE)` is at :259), not :251-258. In the refit at :283-288, `seed = seed` is on :286.
- R/guohe_algorithm3.R:395: the parallel call is `future.apply::future_lapply(...)` at :412-417, with `future.seed = FALSE` at :414. :409-411 is the comment. The cited range :410-415 is slightly off.
- R/setup_gbsg_dgm.R:130 (added note): .create_gbsg_dgm_ has already called set.seed at R/sim_aft_gbsg.R:478 and :604 earlier in the same call, so `.Random.seed` always exists at :128. has_seed is therefore always TRUE in practice, and the `rm(".Random.seed")` branch at :134 cannot run from setup_gbsg_dgm. The restore is still explicit assign/rm and not on.exit, as the investigator said.
- R/mrct_inject_dgm.R:277: in observed-region mode, `if (observed)` is at :470 and `rm_ <- NULL` is at :472 (the investigator cited :472-473). The finding is unchanged.
- R/run_subgroup_sims.R:512: I re-ran the investigator's future check (future 1.80.0, plan(sequential), `%dofuture%` with .options.future seed = NULL) in /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2p1/agents/v_seed. The caller's .Random.seed was identical afterwards, which confirms the restore comes from future and not from package code.
- R/subgroup_consistency_main.R:58 exported_reaching: the static list ('may call') is too broad. .make_candidate_rng_seeds is called only on the parallel batch branch: the `else` at :878, preceded by `if (identical(sg_focus, "maxeff"))` at :757 and `else if (!use_parallel)` at :788; the call itself is at :970-971. use_parallel is set to FALSE when parallel_args$plan is "sequential" (:717-723). These wrappers always force plan = "sequential" on the inner forestsearch call, so they reach :498 but never reach :58: bootstrap_results and forestsearch_bootstrap_dofuture (R/bootstrap_analysis_dofuture.R:595), forestsearch_Kfold (R/forestsearch_cross_validation.R:389), forestsearch_tenfold (:879), run_null_calibration (R/calibrate_null_correction.R:324-326), fs_fdr_report (R/fs_fdr_report.R:359), mrct_region_sims (R/mrct_simulation.R:537; modifyList at :542 keeps plan = "sequential"), and mrct_trial_analysis and mrct_case_study (R/mrct_trial_helpers.R:41). fpr_calibration_transfer is also out: it always passes fs_params = list() (R/fpr_calibration.R:548), and a second fs_params in ... makes do.call fail inside tryCatch, so fpr_calibration then defaults to sequential (:266-268). The site is reached from: forestsearch (default parallel_args = multisession, R/forestsearch_main.R:1634-1636); compare_selection_rules (passes ... into forestsearch, R/compare_selection_rules.R:561); fpr_calibration and run_simulation_analysis, but only when fs_params supplies a non-sequential parallel_args (their defaults are sequential: R/fpr_calibration.R:266-268, R/run_simulation_analysis.R:68); and a direct subgroup.consistency call only with a non-sequential parallel_args that has plan and workers (the default parallel_args = list() at :401 never reaches it).
- R/subgroup_consistency_main.R:58 notes: the claim 'every bootstrap, CV and simulation replicate uses the same master seed for the per-candidate L'Ecuyer streams' is wrong. Bootstrap, CV, mrct, fs_fdr_report and run_null_calibration replicates run the consistency screen with plan = "sequential", so they never build these streams; they use the single global stream re-seeded to 8316951 at :498.
- R/subgroup_consistency_main.R:498 exported_reaching line refs: the forestsearch call cited at R/run_simulation_analysis.R:819 is inside the internal .run_fs_analysis_gen (defined at :773), which run_simulation_analysis calls. bootstrap_results' own call site is R/bootstrap_analysis_dofuture.R:614 (bootstrap_dofuture_main.R:484 is where forestsearch_bootstrap_dofuture calls bootstrap_results). mrct_trial_analysis enters .mrct_analyze_trial at R/mrct_trial_analysis.R:120, and mrct_case_study at R/mrct_case_study.R:301.
- R/subgroup_consistency_main.R:498 notes refinement: under the default consistency_method = "resample" (forestsearch formal, R/forestsearch_main.R:1706), the closed-form consistency (R/subgroup_consistency_helpers.R:1484-1490) draws no random numbers. The :498 stream is used only for consistency_method = "split" or for the split fallback when the closed form returns NA (:1491-1499). The reseed still overwrites the caller's RNG state on every call that reaches :498.
- R/sim_aft_gbsg.R:856/:884 notes addition: with the default sim_id = 1, the sampling seed dgm$seed + 1 is the same seed .create_gbsg_dgm_ set at R/sim_aft_gbsg.R:604 (set.seed(seed + 1)) to draw the super-population's observed event times. This adds to the sim_id = 0 collision with :478.


## Part 2 — vignettes and articles (per file)

### vignettes/forestsearch.qmd

- **Evaluation:** Every chunk is evaluated. The YAML header (:1-18) has no execute: block. Neither vignettes/_quarto.yml nor a root _quarto.yml exists or is tracked. The setup chunk (:20-26) has message: false and warning: false and sets only comment and fig.align through opts_chunk. There are 17 r chunks (:20, :136, :224, :251, :272, :294, :326, :364, :395, :416, :438, :457, :469, :498, :556, :619, :727); the {=html} block at :28 is raw HTML, not R. None of them has eval, include, purl or cache options. The only other option in any chunk is echo: false at :499, which does not affect evaluation. Nothing is read from saved files (no readRDS or load). Every result is computed live from survival::gbsg. The inline `r` expressions are at :155, :259, :355, :356, :448, :450, :493, :505-529, :590-594, :610 and :646-651. They only format values already computed, using sprintf, round, pnorm, lengths, min, max and paste, and draw nothing.
- **Calls reaching the fixed seed:**
  - vignettes/forestsearch.qmd:228 (chunk fit-fs): forestsearch(..., subgroup_method = "consistency", consistency_method = "resample", seedit = 8316951, mr_inference = TRUE). Path: has_subgroups is set at R/forestsearch_main.R:3559 and tested by if (has_subgroups) at :3566. Then sc_filtered_args <- filter_call_args(args_call_all, subgroup.consistency, consistency_overrides) at :3643-3647. args_call_all holds forestsearch's formals only (:1884-1885), and those have seedit (:1640) but no seed. consistency_overrides (:3582-3640) has no seed either, and filter_call_args (R/summary_utility_functions.R:679-687) keeps only target formals. So seed is not forwarded, and do.call(subgroup.consistency, sc_filtered_args) at :3669 runs with the default seed = 8316951 (R/subgroup_consistency_main.R:397), reaching if (!is.null(seed)) set.seed(seed) at :497-498. The vignette shows the stage ran: it reads fit_fs$grp.consistency$n_passed at :281 and fit_fs$grp.consistency$out_sg$result at :295.
- **set.seed() calls and seed arguments:**
  - vignettes/forestsearch.qmd has no set.seed() call.
  - vignettes/forestsearch.qmd:232 seed argument seedit = 8316951 (forestsearch, chunk fit-fs)
  - vignettes/forestsearch.qmd:403 seed argument seedit = 8316951 (forestsearch, chunk fit-dina)
  - vignettes/forestsearch.qmd:424 seed argument seedit = 8316951 (forestsearch, chunk fit-grf)
- **Rows:** 0
- **Notes:** No evaluated random draw follows the fixed-seed call at :228 without its own seed first, so rows = [].

The calls in evaluation order that draw, or could draw:

(1) :228 fit_fs. This is the only call that reaches the fixed seed. After the consistency stage, the multiplier stage reseeds: forestsearch_main.R:3862 calls fs_mr_inference(..., seed = .g_mr(mr_inference_args$seed, seedit)), with the seed argument at :3930. That reaches if (!is.null(seed)) set.seed(seed) at R/fs_mr_inference.R:615, before the multiplier draw at :702 (.fs_mr_multipliers, :112-114). The field block reseeds at :912 with seed + 900000L before its rnorm draws at :919-920. Nothing draws before the fixed seed: use_lasso, use_grf and use_dina are FALSE (forestsearch_main.R:1642-1646), and the setup helpers and subgroup_search.R draw nothing.

(2) :399 fit_dina, forestsearch(subgroup_method = "dina", seedit = 8316951). It does not reach the fixed seed: the DINA branch starts at forestsearch_main.R:2766 and returns at :2887, before Section 9. Its draws are all seeded:
- Cross-fitting folds. .resolve_dina_args sets fit$seed = seedit (R/forestsearch_helpers.R:1154-1157), and .forestsearch_dina_select calls do.call(dina, fit_call). From there: dina (R/dina_df.R:84) -> dina_fit (R/dina.R:157; cross_fitting = TRUE, n_folds = 2L) -> .dina_make_folds(n, n_folds, seed) at R/dina.R:241. That function saves the caller's .Random.seed and restores it on exit (or removes it if absent) (:421-431), then calls set.seed(seed) (:433) and sample.int (:435). It neither reads nor disturbs the caller's stream.
- Everything else on the DINA path is deterministic. dina_subgroup.R, the DINA nuisance fits (glm/coxph), and .dina_reselect_on_effect (forestsearch_helpers.R:1204) with .fs_mr_pieces (fs_mr_inference.R:61) contain no draws.
- Multiplier stage. .fs_apply_mr (R/fs_mr_inference_methods.R:141, seed = .g(mr_inference_args$seed, seedit) at :172) calls fs_mr_inference, which runs set.seed(seed) at :615 before any draw.

(3) :420 fit_grf, forestsearch(subgroup_method = "grf", seedit = 8316951). It does not reach the fixed seed: the GRF branch runs from :2901 and returns at :3072. Its draws are all seeded:
- grf::causal_survival_forest(..., W.hat = 0.5, seed = seedit) at R/grf_helpers.R:74-78. In the installed grf 2.6.1, the nuisance regression_forest is skipped because W.hat is supplied. The survival_forest nuisance fits get seed = seed through args.nuisance, and the main forest gets seed = seed. Since seed is supplied, the default seed = runif(...) is never evaluated, and R's RNG is untouched.
- On the frontier path, grf.subg.harm.survival returns at R/grf_main.R:336, before fit_policy_trees (:340). Only policytree::double_robust_scores (:264) runs, and it is deterministic. .grf_dr_candidates, .grf_dr_candidates_d2, .grf_frontier_select (R/grf_subgroup_labels.R:258, :284, :358) and .grf_reselect_on_effect (forestsearch_helpers.R:1727) contain no draws.
- Multiplier stage: seeded exactly as for DINA.
- Scanning R/grf*.R, R/dina_subgroup.R and R/forestsearch_helpers.R for sample, runif, rnorm, rbinom, rpois and set.seed finds only grf_subg_harm_glm.R:460 set.seed(seedit), which is on the GLM path and not used here.

(4) :558 and :559 fs_declaration_calibration(fit_fs, ...). The function is on the 23-name list only through eval(formals(subgroup.consistency)$pconsistency.digits) at R/fs_declaration_calibration.R:564. That is a lookup of a default, not a call, and it is guarded by if (is.null(digits)) at :563. Here digits = fit_fs$args_call_all$pconsistency.digits = 2, so the lookup is not even evaluated. The body (:514-679), .fs_decl_reduction (:280) with remove_near_duplicate_subgroups, and .fs_decl_c0_block (:680-727) take only quantiles and means of the stored Mstar and Mstar_c0 (the capture is used because declaration_c0 was set at :246). They draw nothing. fs_mr_inference is on the list through the same kind of lookup (fs_mr_inference.R:673), which is likewise not evaluated: fit_fs passes pconsistency.digits at forestsearch_main.R:3929.

The vignette itself makes no base-R draws (no sample, rnorm, runif, rbinom or jitter), never calls simulate_from_dgm, and makes no plots. The other chunks use kable, aggregate, sprintf, pnorm, qnorm and sessionInfo, all deterministic. No path on these calls uses future, so no future seed is taken from the caller's stream. The future.seed and .options.future = list(seed = TRUE) sites in the package (bootstrap, CV, fdr, fpr, dina_bagged, MRCT simulation) are not reached. subgroup.consistency's per-candidate streams (.make_candidate_rng_seeds, R/subgroup_consistency_main.R:40-66) restore the RNG kind and state, and they come after the fixed seed in any case.

Effect of the planned save/restore change: none on this vignette's RNG-dependent outputs. Nothing after :228 reads the caller's stream. fit_dina and fit_grf draw only under their own seedit-derived seeds (DINA's fold draw already saves and restores the caller's state), and the declaration chunks draw nothing. These would therefore stay unchanged and stable from render to render: the identified regions, family sizes, Table 1 (chunk table1, :469-482; also table1-helper :326 and fs-supporting :364), the DINA/GRF selections (:438), the agreement partition (:457), the declaration tables (:556, :619), and every inline value. The elapsed-time lines printed in fs-selected (:251) and dina-grf-selected (:438) already differ on every render, independent of the RNG.

Build: the file is a CRAN-built vignette and also the pkgdown Introduction article. Evidence:
- The vignette-related .Rbuildignore entries are ^vignettes/articles$, ^vignettes/results$, ^vignettes/.*_cache$ and ^vignettes/.*_files$. None of them, nor \.quarto or ^quarto$, matches vignettes/forestsearch.qmd.
- DESCRIPTION:30 has VignetteBuilder: knitr, quarto, with quarto in Suggests (:70).
- The header declares %\VignetteEngine{quarto::html} (:16).
- _pkgdown.yml:32,40 link articles/forestsearch.html.
- **Verifier (corrected) corrections:**
  - The finding stands: rows = [] is right. I found no evaluated draw after :228 that takes from the caller's stream without its own seed. Every correction below is a line number or a clarification.
  - forestsearch_main.R path, line numbers. The filter_call_args( call starts at :3643; :3645 is its subgroup.consistency argument. consistency_overrides runs :3582-3605, with GLM-only additions at :3612-3640, not :3576-3640. The fs_mr_inference( call starts at :3862, but its seed = .g_mr(mr_inference_args$seed, seedit) argument is at :3930.
  - .fs_apply_mr is defined at R/fs_mr_inference_methods.R:141 and calls fs_mr_inference at :157. The seed = .g(mr_inference_args$seed, seedit) argument is at :172, not :173.
  - DINA call chain. dina() is in R/dina_df.R:84 and calls dina_fit at R/dina_df.R:98. dina_fit is defined at R/dina.R:157, with defaults cross_fitting = TRUE and n_folds = 2L at :161-162. It calls .dina_make_folds at R/dina.R:241. .dina_make_folds (:419) saves .Random.seed or registers its removal at :421-431, calls set.seed at :433 and sample.int at :435.
  - GRF. The RCT branch of fit_causal_forest is R/grf_helpers.R:74-78, with seed = seedit at :75. W.hat = 0.5 is supplied, so the installed grf 2.6.1 skips the nuisance regression_forest (the body tests is.null(W.hat)). Only the survival_forest nuisance fits (sf.survival, sf.censor) run, and they get seed = seed through args.nuisance. With grf_selection = "frontier", grf.subg.harm.survival returns at R/grf_main.R:336, before fit_policy_trees at :340. So policytree::policy_tree is never called; only double_robust_scores at :264 runs, and it is deterministic.
  - The formals() lookups are not even evaluated here. In fs_declaration_calibration, R/fs_declaration_calibration.R:564 sits inside if (is.null(digits)) at :563, and digits = fit_fs$args_call_all$pconsistency.digits = 2. In fs_mr_inference, :673 sits inside if (is.null(digits)), and fit_fs passes pconsistency.digits = pconsistency.digits (= 2) at forestsearch_main.R:3929. On the DINA and GRF calls, :673 is never reached, because the admission set has no consistency floor.
  - Displayed outputs. The planned change leaves the selections, family sizes, Table 1, the agreement partition and both declaration tables unchanged and stable. But the fs-selected (:251) and dina-grf-selected (:438) chunks also print elapsed fit times (fit_sec_fs, g_fs$timing_seconds, fit_sec_dina, fit_sec_grf). Those differ on every render anyway, for reasons that have nothing to do with the RNG.
  - Inline expressions. The ones near :448-450 are on lines :448 and :450 only; :449 has none.
  - Added: fit_fs does not read the document's own stream before its fixed-seed call either. use_lasso, use_grf and use_dina are FALSE (forestsearch_main.R:1642-1646). The setup helpers (.validate_mr_configuration, .fs_family_status, .fs_resolve_admission in R/forestsearch_helpers.R) and subgroup_search.R draw nothing. lasso_selection, which is off here, would scope its own seed at R/get_FSdata_helpers.R:245-259. So no output depends on the render's unseeded initial state.

### vignettes/articles/actg175-glm-outcomes.qmd

- **Evaluation:** Every chunk is evaluated. The YAML header (:1-18) has no execute: or eval: option. The setup chunk (:20-30) sets only comment and fig.align (:23), with message/warning false. No chunk has eval: false, eval: !expr, include: false, cache or purl. Two chunks are echo: false (s16-values :304, s3-values :562). The five fit chunks are results: hide (:153, :181, :419, :439, :458), and all of them still run. No result comes from a saved file: there is no readRDS(), load() or cache. A tracked, pre-rendered actg175-glm-outcomes.html sits beside the .qmd, but nothing reads it. quarto/_quarto.yml does not cover vignettes/articles. There are 21 inline `r` expressions (:310, :311, :312, :314, :315, :317, :318, :320-:324, :326, :395 x2, :397 x2, :571, :581, :582, :604). They only format values that already exist. The working-tree file is byte-identical to git show c40d5309:vignettes/articles/actg175-glm-outcomes.qmd, and it is tracked.
- **Calls reaching the fixed seed:**
  - actg175-glm-outcomes.qmd:155 forestsearch() (chunk fit-bin), subgroup_method = 'consistency' (:165), resolved by match.arg at forestsearch_main.R:1737. A region is found (stopifnot :210), so Section 9 runs do.call(subgroup.consistency, sc_filtered_args) (forestsearch_main.R:3669). sc_filtered_args comes from filter_call_args (summary_utility_functions.R:679), and args_call_all holds only forestsearch's formals (mget at forestsearch_main.R:1884-1885), so no 'seed' is passed. set.seed(8316951) then runs at subgroup_consistency_main.R:497-499 (formal :397). The MR stage follows: set.seed(seedit) at fs_mr_inference.R:615 (seed passed at forestsearch_main.R:3930), then the field stage's set.seed(seedit + 900000L) at fs_mr_inference.R:912 (ci_method = 'field', :171).
  - actg175-glm-outcomes.qmd:183 forestsearch() (chunk fit-cont): default subgroup_method 'consistency', region found (stopifnot :211), same path to subgroup_consistency_main.R:498. Then fs_mr_inference.R:615 and :912 (ci_method = 'field', :198).
  - actg175-glm-outcomes.qmd:421 forestsearch() (chunk fit-s3-fs): subgroup_method = 'consistency', region = Part 1's (stopifnot :477), same path to subgroup_consistency_main.R:498. Then fs_mr_inference.R:615, and :912 through the default ci_method 'field' (forestsearch_main.R:3888).
  - actg175-glm-outcomes.qmd:522 forestsearch_bootstrap_dofuture() (chunk bootstrap), loop iteration m = 'FS'. Each replicate replays forestsearch(subgroup_method = 'consistency') through do.call(forestsearch, args_FS_boot) at bootstrap_analysis_dofuture.R:614, inside foreach %dofuture% (:408, seed = TRUE :411), so each replicate reaches subgroup_consistency_main.R:498. With n_workers = 2 these run in multisession workers. The call also sets the main process's RNG through its explicit seed = 8316951L: set.seed at bootstrap_calculations_helpers.R:27 (via bootstrap_dofuture_main.R:389-393) and at bootstrap_analysis_dofuture.R:334. The DINA and GRF iterations replay the dina and grf branches, which never reach subgroup.consistency.
  - Not fixed-seed calls in the subgroup.consistency sense: :441 forestsearch(subgroup_method = 'dina') returns at forestsearch_main.R:2887, and :460 forestsearch(subgroup_method = 'grf') returns at :3072, both before Section 9. Each still seeds itself from the document's explicit seeds. DINA: set.seed(dina_args seed) at dina.R:433, with the caller's state saved and restored at dina.R:420-431 (via dina.R:241). GRF: set.seed(seedit) at grf_subg_harm_glm.R:460, and causal_forest(seed = seedit) at :463, the RCT branch. MR for both: set.seed(seedit) at fs_mr_inference.R:615 (seed at fs_mr_inference_methods.R:172), with ci_method 'ij' (:171), so no field reseed.
  - On the static reach list but not reaching the seed: :346 fs_declaration_table(), called from s17_column() at :366 and :367. It calls fs_declaration_calibration() twice (fs_declaration_table.R:68-69). That function's only code link to subgroup.consistency is eval(formals(subgroup.consistency)$pconsistency.digits) at fs_declaration_calibration.R:564, which is not a call. Neither file has a draw site; its remove_near_duplicate_subgroups (:315) does not draw either.
- **set.seed() calls and seed arguments:**
  - No set.seed() call anywhere in the document.
  - actg175-glm-outcomes.qmd:159 forestsearch(seedit = 8316951) (fit-bin; reaches fs_mr_inference.R:615 via forestsearch_main.R:3930, and :912)
  - actg175-glm-outcomes.qmd:187 forestsearch(seedit = 8316951) (fit-cont; same as :159)
  - actg175-glm-outcomes.qmd:425 forestsearch(seedit = 8316951) (fit-s3-fs; same as :159)
  - actg175-glm-outcomes.qmd:445 forestsearch(seedit = 8316951) (fit-s3-dina; MR seed via fs_mr_inference_methods.R:172 -> fs_mr_inference.R:615)
  - actg175-glm-outcomes.qmd:448 dina_args = list(seed = 8316951L, n_folds = 2L) (fit-s3-dina; .resolve_dina_args forestsearch_helpers.R:1155-1158 -> dina_fit -> .dina_make_folds, set.seed at dina.R:433 with save/restore)
  - actg175-glm-outcomes.qmd:464 forestsearch(seedit = 8316951) (fit-s3-grf; set.seed(seedit) at grf_subg_harm_glm.R:460, causal_forest(seed = seedit) at :463; MR at fs_mr_inference.R:615)
  - actg175-glm-outcomes.qmd:523 forestsearch_bootstrap_dofuture(seed = 8316951L) (chunk bootstrap, used in each of the 3 loop iterations; bootstrap_calculations_helpers.R:27, bootstrap_analysis_dofuture.R:334)
- **Rows:** 0
- **Notes:** No such place exists in this file, so rows = []. The first fixed-seed call is :155 (fit-bin). Every later call that uses the RNG sets its own seed before its first draw, and that seed comes from an argument the document passes explicitly.
(1) :183 fit-cont and :421 fit-s3-fs: nothing draws between entry and subgroup_consistency_main.R:498.
- forestsearch_main.R has no draw site.
- use_lasso, use_grf and use_dina are FALSE (:1642, :1643, :1646); lasso_selection (get_FSdata_helpers.R:241-259) is gated off at get_fsdata.R:326 and saves/restores anyway.
- ps_method resolves to 'none' because is.RCT = TRUE (forestsearch_main.R:2602).
- get_fsdata.R and subgroup_search.R have no draw sites. subgroup.search's future_lapply (subgroup_search.R:287) has no future.seed and runs with parallel_workers = 1 (forestsearch_main.R:3469).
- subgroup.consistency lines :404-:497 only validate.
- The split/resample draw sites (subgroup_consistency_helpers.R:272, :1733; consistency_resample.R:459-464) are reached only inside subgroup.consistency, after :498.
(2) :346, :366, :367 declaration chunk: arithmetic on the captured field only. The document asserts source == 'capture' at :347.
(3) :441 fit-s3-dina: the DINA path makes one draw, the fold permutation sample.int at dina.R:435, under its own seed with save/restore. dina_subgroup.R, dina_df.R and forestsearch_helpers.R have no draw sites. dina_fit_bagged and dina_subgroup_bootstrap are not called. MR then reseeds.
(4) :460 fit-s3-grf: nothing draws before set.seed(seedit) at grf_subg_harm_glm.R:460.
(5) :522 bootstrap:
- Before its set.seed, the call runs only validation, make_effect_estimator, fit_effect_models (GLM fits) and setup_parallel_SGcons (bootstrap_dofuture_main.R:372 -> future::plan(multisession) at subgroup_consistency_helpers.R:57).
- I re-ran the check (Rscript --vanilla, future 1.80.0, parallelly 1.48.0): .Random.seed was identical after plan(multisession, workers = 2), after a future and after resetting to sequential.
- In any case nothing displayed depends on draws made before the explicit set.seed.
- doFuture seed = TRUE builds its L'Ecuyer streams after set.seed(seed) (bootstrap_analysis_dofuture.R:334). Replicate rows come from the pre-generated boot_index_mat.
(6) Later chunks make no draws: mr_method_comparison and mr_cost_table (mr_method_tables.R), mr_values (mr_values.R), fs_in_subgroup (mrct_minor_helpers.R:77 -> get_dfpred), kable and sessionInfo. Every rnorm/rbinom/set.seed in mr_method_tables.R and mr_values.R is in roxygen examples.
Effect of the planned save/restore change: none of the displayed outputs change, either from one render to the next or once. That covers:
- the data chunk's counts and whole-trial effects, computed before any fixed-seed call;
- the Table S16 and S17 kables and the 21 inline values;
- the fit summaries printed by the regions and s3-selected chunks, whose timings vary anyway;
- the s3-columns table, including the bootstrap point estimates;
- the cost table and round(boot_secs, 1), whose timings vary anyway;
- the three static published tables.
Only the RNG state left after each call would differ, and nothing later reads it. The DINA fold helper (dina.R:419-433) and lasso_selection already use the save/restore pattern.
Build: pkgdown only. .Rbuildignore:20 '^vignettes/articles$' keeps the article out of the CRAN tarball, even though DESCRIPTION:30 declares VignetteBuilder: knitr, quarto. _pkgdown.yml lists it at :49 (navbar href) and :584 (articles index).
- **Verifier (corrected) corrections:**
  - The file has 21 inline `r` expressions, not 19. The list left out :604 `r plan_label` in the cost prose, which formats n_workers from the bootstrap chunk and is deterministic. In full: :310, :311, :312, :314, :315, :317, :318, :320, :321, :322, :323, :324, :326, :395 (x2), :397 (x2), :571, :581, :582, :604.
  - subgroup_method's match.arg is at forestsearch_main.R:1737, not :1736. Line :1736 is grf_count_transform's match.arg. The formal is at :1651, as stated.
  - One internal seed was left out. The three consistency-path fits (:155, :183, :421) run the field stage, which reseeds with set.seed(as.integer(seed) + 900000L) at fs_mr_inference.R:912, after :615. fit-bin and fit-cont pass ci_method = 'field' explicitly (:171, :198). fit-s3-fs gets 'field' as the consistency-path default at forestsearch_main.R:3888. The DINA and GRF fits run MR with ci_method defaulting to 'ij' (fs_mr_inference_methods.R:171), so only :615 runs there. Rows are unchanged.
  - GRF: only grf_subg_harm_glm.R:463, the RCT branch with W.hat = 0.5, runs. Line :467 is the is.RCT = FALSE branch and does not run, because is.RCT = TRUE is forwarded at grf_args.R:60. The investigator's ':463/:467' should be read as two alternatives, not two calls.
  - Location of the bootstrap replicate call: each replicate runs do.call(forestsearch, args_FS_boot) at bootstrap_analysis_dofuture.R:614, inside the foreach at :408 with seed = TRUE at :411. Replicates are forced to plan 'sequential' (:595-597) and mr_inference = FALSE (:607). The bootstrap row indices come from boot_index_mat, built once in the main process under set.seed(seed) (bootstrap_calculations_helpers.R:27-29, called from bootstrap_dofuture_main.R:389-393). So the per-iteration future streams decide nothing that is displayed.
  - Added context, with no change to the evaluation status: a pre-rendered vignettes/articles/actg175-glm-outcomes.html is tracked at c40d5309 (git ls-tree). The .qmd never reads it, and every chunk is still computed live.

### vignettes/articles/biomarker_effects.qmd

- **Evaluation:** No document-level execute or eval option is set. The YAML front matter (lines 1-19) has only format: html options (embed-resources, toc, toc-depth, number-sections, code-overflow, fig-width/height) and the vignette block (VignetteEngine quarto::html at :17). No _quarto.yml or _metadata.yml in vignettes/articles/, vignettes/ or the package root could add one. The only _quarto.yml in the tree is quarto/_quarto.yml, which applies to the separate quarto/ project and not to vignettes/articles. The root .Rprofile only sets options(pkgdown.internet = FALSE). The file has five chunks, all evaluated: none has eval: false, an eval: !expr option, purl or include, and there is no inline `r ` code. (1) The {css, echo=FALSE} chunk at :21-34 is a CSS style block with no R code. (2) The {r param-table, echo=FALSE} chunk at :171-189 is a knitr::kable() of a literal data.frame. (3) The {r working-example, echo=FALSE} chunk at :307-324 is a knitr::kable() of a literal data.frame. (4) The {r pipeline-table, echo=FALSE} chunk at :561-585 is a knitr::kable() of a literal data.frame. (5) The {r correspondence-table, echo=FALSE} chunk at :608-630 is a knitr::kable() of a literal data.frame. No chunk calls library(), readRDS(), load() or source(), so no results are read from saved files. Every evaluated R chunk only formats a hard-coded table, and none draws random numbers.
- **Calls reaching the fixed seed:**
- **set.seed() calls and seed arguments:**
- **Rows:** 0
- **Notes:** No such place exists, so rows = []. No evaluated call reaches the fixed seed: no chunk calls forestsearch(), subgroup.consistency() or any of the other 21 names in reach_fixed_seed.txt, and no chunk loads forestsearch. The document draws no random numbers: there is no sample, r* distribution call, simulate_from_dgm, future/doFuture or jitter. It also has no set.seed() call and no seed argument. Searching the file for seed|sample|rnorm|runif|rbinom|rexp|forestsearch(|consistency|simulate|jitter|readRDS|load(|source(|execute|eval|library|require finds only prose lines (64, 90, 148, 235, 240, 267, 550), none inside a chunk. Package functions are named only in prose, never called in code: create_gbsg_dgm() and generate_aft_dgm_flex() at :66 and :587, cox_ahr_cde_analysis() and calculate_hazard_ratios() at :590-591. The four R chunks call only knitr::kable() on literal data.frames (:172, :308, :562, :609). The planned save/restore of .Random.seed therefore changes no output of this article. All four kable tables (param-table, working-example, pipeline-table, correspondence-table) stay byte-stable across renders before and after the change. Build: the file is in vignettes/articles/, which .Rbuildignore:20 (^vignettes/articles$) excludes from the source tarball. Even though DESCRIPTION:30 has VignetteBuilder: knitr, quarto and the YAML declares VignetteEngine quarto::html, CRAN never builds it. pkgdown alone builds it: _pkgdown.yml:53 (href articles/biomarker_effects.html) and :592 ("articles/biomarker_effects"). For reference, the fixed seed is confirmed at R/subgroup_consistency_main.R:397 (seed = 8316951) and :498 (set.seed(seed), inside if (!is.null(seed)) at :497). The do.call(subgroup.consistency, sc_filtered_args) call is at R/forestsearch_main.R:3669. Neither is reached from this document.
- **Verifier (confirmed) corrections:**

### vignettes/articles/causal_effects_brief_review.qmd

- **Evaluation:** The YAML header (:1-19) has no execute/eval option. There is no _quarto.yml, _metadata.yml or _freeze in vignettes/articles, vignettes/ or the package root. knitr defaults therefore apply and every chunk runs. The setup chunk at :21-28 (include=FALSE, so its code runs but shows nothing) calls knitr::opts_chunk$set(echo = FALSE, message = FALSE, warning = FALSE, fig.align = "center") at :22-27 and does not turn off eval. There are exactly four chunks, all R and all evaluated: setup :21-28; estimand-taxonomy-table :108-195 (library(gt) :109, literal data.frame :111-143, gt(tax) pipeline :145-194); reconciliation-table :432-533 (literal data.frame :433-500, gt(df) :502-532); estimand-mapping-table :789-913 (literal data.frame :790-882, gt(df2) :884-912). No chunk sets eval: false, #| options, purl or a non-R engine. There are no inline `r ` expressions. The backticked `hr.H.hat` (:955, :970) and `dgm$hr_H_true` (:955, :973) are plain inline code and do not run. Nothing loads precomputed results: there is no readRDS, load() or source(). forestsearch is never loaded, and no forestsearch function is called.
- **Calls reaching the fixed seed:**
- **set.seed() calls and seed arguments:**
- **Rows:** 0
- **Notes:** No place in this file has a random draw following a call that reaches the fixed seed, so rows = []. No evaluated code reaches the fixed seed: grepping the file for 'name(' for each of the 23 names in reach_fixed_seed.txt finds zero hits. The file has no set.seed() and no seed argument. Correction: the file does make random draws. Each gt table is printed without an id. gt's as.tags.gt_tbl then calls random_id(), which is paste(sample(letters, n, replace = TRUE), collapse = "") with n = 10 (installed gt 1.4.0), so three sample() draws per render come from the global stream, one at each of :145, :502 and :884. They set each table's container div id and the selector of its scoped CSS. The tracked HTML shows id="duxkegtsbt", "hjqsfyvuca" and "rksvpglvjr", with matching '#<id> table' rules. No fixed-seed call precedes them, and the document seeds nothing, so they draw from R's time- and process-id-seeded stream and already differ from one render to the next. They do not change the visible content of the tables. The table contents are literal strings and deterministic. The planned save/restore of .Random.seed in forestsearch's exported entry points cannot affect this file. For reference, the fixed seed is set at R/subgroup_consistency_main.R:498, set.seed(seed) inside if (!is.null(seed)) at :497, with default seed = 8316951 at :397. forestsearch reaches it through R/forestsearch_main.R:3669, do.call(subgroup.consistency, sc_filtered_args). Build: pkgdown only, not CRAN. .Rbuildignore:20 is ^vignettes/articles$, so R CMD build drops the articles directory, even though DESCRIPTION:30 declares VignetteBuilder: knitr, quarto and the YAML at :14-18 declares the %\VignetteEngine{quarto::html} vignette (engine line :17). _pkgdown.yml lists the article in the navbar at :54-55 (href articles/causal_effects_brief_review.html) and in the articles index at :593. A pre-rendered .html is also tracked.
- **Verifier (corrected) corrections:**
  - Line numbers for the inline code: `hr.H.hat` appears at :955 and :970, not :971. `dgm$hr_H_true` appears at :955 and :973. Plain inline code, not evaluated (grep -n 'hr.H.hat\|hr_H_true').
  - The notes say the file 'contains no random draws'. That is wrong. Each of the three gt tables (gt(tax) at :145, gt(df) at :502, gt(df2) at :884) is passed no id. When knitr prints it, gt's as.tags.gt_tbl calls random_id() because table_id is NA, and random_id is paste(sample(letters, n, replace = TRUE), collapse = "") with n = 10 (checked in the installed gt 1.4.0). That is three draws from the global RNG stream. The draws set the HTML id of each table's container div and the selector of its scoped CSS. The tracked pre-rendered vignettes/articles/causal_effects_brief_review.html shows them as id="duxkegtsbt", "hjqsfyvuca", "rksvpglvjr" with matching '#<id> table' CSS rules. Rows stays [] because no fixed-seed call comes before these draws. They come from the session's own unseeded stream: there is no set.seed in the file, so R seeds from time and process id. These ids therefore already differ on every render and do not change the visible table content. The planned save/restore cannot affect them.

### vignettes/articles/extreme_subgroups.qmd

- **Evaluation:** Every R chunk is evaluated. The YAML execute block (:7-10) sets only echo: true, warning: false and message: false. No chunk sets eval:, purl or eval: !expr; a grep for 'eval' finds only the gbsg_eval object and eval(parse()) in code (:1057-1071, :1237). Three chunk options hide output but still run the code: '#| include: false' on ub-forest-highrisk-setup (:1673), '#| echo: false' on timing-summary (:1837) and '#| message: false' on parallel-setup (:1144). The fences at :218-220, :224-229 and :259-261 are plain ```r display blocks, not {r} chunks. Nothing configures execution from outside the file: there is no _quarto.yml or _metadata.yml in vignettes/articles/ or at the repository root, and _pkgdown.yml has no execute or eval settings. Nothing is read back: there is no readRDS, load() or file read, so every value, including the 150-trial %dofuture% loop in uniform-sims (:1251-1288), is computed live. The file is tracked (git ls-files) and is byte-identical to c40d5309:vignettes/articles/extreme_subgroups.qmd (cmp).
- **Calls reaching the fixed seed:**
- **set.seed() calls and seed arguments:**
  - vignettes/articles/extreme_subgroups.qmd:559 seed = 99 in base_args_uniform. calibrate_k_treat() (:572) carries it into every objective evaluation (R/calibrate_k_treat.R:133-136, modifyList keeps seed), and do.call(generate_aft_dgm_flex, ...) (:582-585) uses it too. Each call runs set.seed(seed) at R/generate_aft_dgm_main.R:248 before drawing (generate_super_population at :338, calculate_hazard_ratios at :345)
  - vignettes/articles/extreme_subgroups.qmd:708 simulate_from_dgm(seed = 1), call at :699, chunk pilot-sim. It runs if (!is.null(seed)) set.seed(seed) at R/simulate_from_dgm.R:208, and every draw in the function (:251-331) comes after that line
  - vignettes/articles/extreme_subgroups.qmd:778 calibrate_cens_adjust(seed = 42), call at :771, chunk calibrate. The formal default is also 42 (R/simulate_from_dgm.R:795). The seed is passed to simulate_from_dgm at :859 (objective) and :893 (final)
  - vignettes/articles/extreme_subgroups.qmd:813 simulate_from_dgm(seed = 123), call at :804, chunk verify-cal
  - vignettes/articles/extreme_subgroups.qmd:858 simulate_from_dgm(seed = 100), call at :849, chunk single-trial
  - vignettes/articles/extreme_subgroups.qmd:959 set.seed(8316951), chunk define-uniform-subgroups. This is the document's own explicit seed, not a forestsearch call. Nothing in the main session draws afterwards before :1905, so no output depends on it
  - vignettes/articles/extreme_subgroups.qmd:1209 simulate_from_dgm(seed = sim_config$seed_base + ss) inside run_one_sim() (call at :1203). It runs in each worker, once per iteration dispatched from :1268. seed_base = 0L is set at :674
  - vignettes/articles/extreme_subgroups.qmd:1220 set.seed(sim_config$seed_base + ss + sim_config$rand_seed_offset) inside run_one_sim(), before sample.int at :1222. rand_seed_offset = 1e6L is set at :675
  - vignettes/articles/extreme_subgroups.qmd:1264 .options.future = list(seed = NULL) on foreach %dofuture%. doFuture 1.4.0 maps NULL to future's default seed = FALSE (doFuture2 deparsed lines 147-155), so make_rng_seeds() returns NULL and no streams are derived from the caller's RNG; each iteration seeds itself at :1209 and :1220
  - vignettes/articles/extreme_subgroups.qmd:1905 set.seed(42), chunk hte-build-dgm. It is redundant: the next draw comes from generate_aft_dgm_flex, which reseeds itself
  - vignettes/articles/extreme_subgroups.qmd:1935 generate_aft_dgm_flex(seed = 42), call at :1907 (set.seed at R/generate_aft_dgm_main.R:248)
  - vignettes/articles/extreme_subgroups.qmd:1981 calibrate_cens_adjust(seed = 42), call at :1974, chunk hte-calibrate
  - vignettes/articles/extreme_subgroups.qmd:2003 simulate_from_dgm(seed = 100), call at :1994, chunk hte-single-trial
- **Rows:** 0
- **Notes:** There are no rows. No evaluated code in this file reaches the fixed seed, so the planned save/restore of .Random.seed changes nothing here. All displayed outputs are already identical from one render to the next, because every draw is seeded, and the change will not alter them.

1. None of the 23 names in reach_fixed_seed.txt is called in an evaluated chunk. The word 'forestsearch' appears at :26-29 (version check and stop() strings), :89, :211, :213, :219, :226, :255, :260 (prose and display blocks), :271 (library), :600, :604 (comment and stop() string), :908, :1265 (packages option) and :2076 (prose naming `forestsearch()`). None of these is a call.

2. The package functions the document does call are:
- generate_aft_dgm_flex: :582-583 via do.call, and :1907.
- calibrate_k_treat: :572.
- simulate_from_dgm: :699, :804, :849, :1203 and :1994.
- check_censoring_dgm: :720 and :820.
- calibrate_cens_adjust: :771 and :1974.
- create_summary_table: :313.
- gg_forest: :1467, :1513, :1571, :1620, :1705 and :1762.

They are defined at R/generate_aft_dgm_main.R:218, R/calibrate_k_treat.R:96, R/simulate_from_dgm.R:172, :440 and :789, R/create_summary_table.R:75 and R/gg_forest.R:115. None is in the reach list.

I ran an independent static check (scratch script /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2ro/extreme_subgroups_verify/cg.R). It parses every R/*.R file, builds a graph from all.names() plus string constants (so do.call with a string is caught), and runs BFS. It finds that none of the seven functions reaches subgroup.consistency. As sanity checks, forestsearch and fs_oc_grid do reach it. The fixed seed itself is at R/subgroup_consistency_main.R:397 (seed = 8316951) and :497-499 (if (!is.null(seed)) set.seed(seed)).

3. The RNG primitives inside those call graphs are:
- generate_aft_dgm_flex: set.seed (:248), then sample (R/generate_aft_dgm_helpers.R:858) and rexp (:1080). All three run after the set.seed.
- simulate_from_dgm: set.seed (:208), then sample, runif, rbinom, rexp and rnorm (:251-331). All run after the set.seed.
- create_summary_table: fisher.test(simulate.p.value = TRUE) at R/create_summary_table.R:365. This branch is not taken for this input: the smallest grade-by-hormon cell is 33, so chisq.test at :359 runs instead.
- check_censoring_dgm and gg_forest: none.

Every draw is therefore seeded by an argument given in the document.

4. The document's set.seed(8316951) at :959 uses the same value as subgroup.consistency's default and as generate_aft_dgm_flex's default (R/generate_aft_dgm_main.R:238). It is an explicit document seed, not a package call. Nothing between :959 and the next seed draws in the main session: median and quantile are deterministic; plan, foreach and %dofuture% with seed = NULL do not consume the caller's RNG; coxph, ggplot2/gg_forest and kable do not draw. set.seed(42) at :1905 is likewise followed only by self-seeded package calls.

5. The %dofuture% loop works as follows. doFuture2 replaces seed = NULL with FALSE, so no per-element seeds are generated, and the caller's RNG is neither read nor advanced (the oseed save/restore at deparsed lines 161-162 sits only in the seeds-non-NULL branch). Each worker iteration calls set.seed itself (:1209 via simulate_from_dgm; :1220). With seed = FALSE, future's RNG-misuse check may raise a warning because the futures change the RNG. I did not verify this. Either way it is hidden by warning: false and does not affect values.

6. Build: .Rbuildignore:20 has ^vignettes/articles$, so the article is not in the CRAN tarball. DESCRIPTION:30 VignetteBuilder: knitr, quarto does not apply to articles/. The article is built by pkgdown only (_pkgdown.yml:58-59 navbar href articles/extreme_subgroups.html; :599 "articles/extreme_subgroups").
- **Verifier (corrected) corrections:**
  - The conclusion stands: nothing evaluated in this file reaches the fixed seed, there are no rows, and the planned save/restore of .Random.seed changes nothing here. All corrections below are line-number or completeness fixes.
  - check_censoring_dgm's second call is at :820 (diag_cal <- check_censoring_dgm(), chunk verify-cal), not :818. Line :818 is the closing parenthesis of the simulate_from_dgm() call above it.
  - calibrate_cens_adjust's first call is at :771 (cal_uniform <- calibrate_cens_adjust(), chunk calibrate), not :772. Line :772 is the dgm = dgm_uniform argument.
  - The HTE call to generate_aft_dgm_flex is at :1907 (dgm_hte <- generate_aft_dgm_flex(), chunk hte-build-dgm), not :1906. Line :1906 is blank. The do.call at :582 names the function on :583.
  - The simulate_from_dgm locations were given as chunk spans (:696-717, :803-828, :848-872, :1993-2018). The calls themselves start at :699, :804, :849, :1203 and :1994. The gg_forest calls were not given lines; they are at :1467, :1513, :1571, :1620, :1705 and :1762.
  - The grep list of 'forestsearch' hits misses :28 and :29. Both are string literals inside stop(), which does not run when the version check passes. :1265 appeared both in the prose list and as the packages option. None of these lines is a call.
  - Added mechanism for :1264 seed = NULL. In the installed doFuture 1.4.0, doFuture2 replaces a NULL seed with formals(future)$seed, which is FALSE (deparsed doFuture2 lines 147-155). make_rng_seeds(nX, seed = FALSE) then returns NULL, so the caller's .Random.seed is never read or advanced. The conclusion is unchanged. This evidence comes from the installed package, not from repository source.
  - Added completeness for create_summary_table (:313). It can draw from the caller's RNG through fisher.test(tab, simulate.p.value = TRUE) at R/create_summary_table.R:365. That branch runs only for a multi-level categorical variable with a cell count below 5. For grade by hormon in gbsg the smallest cell is 33 (cells 48/33, 281/163, 111/50), so the chisq.test branch at :359 runs and nothing is drawn. The binary meno uses fisher.test or chisq.test without simulation (:274/:279). The investigator's statement that every random step has its own seed therefore holds, but this function had not been checked.
  - Added context. generate_aft_dgm_flex's own default seed is also 8316951 (R/generate_aft_dgm_main.R:238), but every call in the document overrides it: seed 99 through base_args, and seed 42. calibrate_k_treat's objective eff_at (R/calibrate_k_treat.R:133-136) uses modifyList, which keeps base_args$seed = 99, so set.seed(99) runs at R/generate_aft_dgm_main.R:248 on every root-finding evaluation. generate_aft_dgm_flex does not reach subgroup.consistency, so it is not a fixed-seed call.

### vignettes/articles/forestsearch-2024.Rmd

- **Evaluation:** There is no document-level eval or execute option. The YAML is :1-20, and the setup chunk (:33, include=FALSE, still evaluated) sets only echo, warning, message and fig.* via opts_chunk at :34-44, plus a code_fold knit hook at :48-58 and timings at :61-62, none of which draw. Chunk eval settings: bootstrap :470 eval=TRUE, bootstrap-plots :585 eval=TRUE, kfold-cv :598 eval=TRUE, oob-cv :640 eval=FALSE (so forestsearch_Kfold :644, forestsearch_KfoldOut :660 and cv_summary_tables :668 are NOT evaluated). timing-summary :900 has echo=FALSE, code_fold=TRUE and is evaluated. Every other chunk is evaluated. Inline code at :5 (Sys.Date) and :541 (figure_note) is evaluated and draws nothing. Nothing is read from saved files: there is no readRDS, load, cache or child, and every result is computed live.
- **Calls reaching the fixed seed:**
  - forestsearch-2024.Rmd:301 forestsearch(), main process. args_call_all holds formals only (R/forestsearch_main.R:1884-1885), and forestsearch has no 'seed' formal (it has seedit), so filter_call_args (R/summary_utility_functions.R:679-688) passes no seed. R/forestsearch_main.R:3669 do.call(subgroup.consistency, sc_filtered_args) -> R/subgroup_consistency_main.R:397 seed = 8316951 -> :498 set.seed(seed), in the main process. Defaults: subgroup_method 'consistency' (:1651), mr_inference = FALSE (:1723). A subgroup is found: bootstrap :478 would stop() otherwise (R/bootstrap_dofuture_main.R:244-251), and :872 reads fs$grp.consistency$out_sg.
  - forestsearch-2024.Rmd:478 forestsearch_bootstrap_dofuture(). Path: R/bootstrap_dofuture_main.R:484 bootstrap_results() -> R/bootstrap_analysis_dofuture.R:408-415 foreach(.options.future = list(seed = TRUE, globals = TRUE)) %dofuture% -> :614 do.call(forestsearch, args_FS_boot) -> forestsearch_main.R:3669 -> subgroup_consistency_main.R:498, inside the 2 multisession workers. In the main process the call's own seed = 8316951L (:481) resets the stream: set.seed at R/bootstrap_calculations_helpers.R:27 (via .make_boot_index_matrix, bootstrap_dofuture_main.R:389), then sample.int :29, then set.seed(seed) at bootstrap_analysis_dofuture.R:334, then the %dofuture% seed derivation moves it forward 2 sample.int(1) steps.
  - forestsearch-2024.Rmd:606 forestsearch_tenfold(). Path: R/forestsearch_cross_validation.R:951-955 foreach(.options.future = list(seed = TRUE)) %dofuture% -> :1002 do.call(forestsearch, cv_args) -> subgroup_consistency_main.R:498, inside the workers. The main process has no set.seed. future::make_rng_seeds (via as_lecyer_cmrg_seed(TRUE)) derives L'Ecuyer streams from the caller's .Random.seed. The caller's stream then moves forward 2 sample.int(1) steps: one next_random_seed() in future:::make_rng_seeds and one in doFuture:::doFuture2 (future 1.80.0, doFuture 1.4.0; checked empirically).
  - NOT evaluated: forestsearch-2024.Rmd:644 forestsearch_Kfold() (chunk oob-cv, eval = FALSE at :640)
- **set.seed() calls and seed arguments:**
  - The file contains no set.seed() call (grep for 'seed' finds only the 5 argument lines below).
  - forestsearch-2024.Rmd:245 seedit = 8316951 -> grf.subg.harm.survival. Passed to grf as seed = seedit (R/grf_helpers.R:75, :85). There is no R-level set.seed, and grf::causal_survival_forest does not draw from R's stream when seed is supplied (its runif default is never evaluated).
  - forestsearch-2024.Rmd:322 seedit = 8316951 -> forestsearch. Used internally, e.g. set.seed(seedit) at R/get_FSdata_helpers.R:259 inside lasso_selection's save/restore of .Random.seed (:248-258).
  - forestsearch-2024.Rmd:481 seed = 8316951L -> forestsearch_bootstrap_dofuture (set.seed at R/bootstrap_calculations_helpers.R:27 and R/bootstrap_analysis_dofuture.R:334).
  - forestsearch-2024.Rmd:611 seed = 8316951L -> forestsearch_tenfold. Used only inside the foreach body: set.seed(seed + 1000 * ksim) at R/forestsearch_cross_validation.R:959.
  - NOT evaluated: forestsearch-2024.Rmd:647 seedit = 8316951L -> forestsearch_Kfold (oob-cv eval = FALSE).
  - Implicit, not written in the file: forestsearch-2024.Rmd:812 plot_km_band_forestsearch -> weightedsurv::plotKM.band_subgroups (R/plot_km_band_forestsearch.R:349; no seedstart passed, so the default seedstart = 8316951 applies) -> KM_diff, which runs set.seed(seedstart) before rnorm when draws.band = 20 > 0 (installed weightedsurv 0.1.0).
  - Build tool, not the file: pkgdown 2.2.1 build_articles(seed = 1014L) -> build_article -> build_rmarkdown_article -> callr::r_safe(rmarkdown_render_with_seed), which runs set.seed(1014) and htmlwidgets::setWidgetIdSeed(1014) before rmarkdown::render.
- **Rows:** 5
  - Row 1:
    - file_line: forestsearch-2024.Rmd:432
    - later_call: knitr auto-print of res_tabs$sg10_out, a gt_tbl built in sg_tables' maxk == 2 branch (R/summary_utility_functions.R:1193) at :1196-1197 (sg10 <- sg10[, c('M.1','M.2','N','E','d1','hr','Pcons')] at :1194). Call chain: gt:::knit_print.gt_tbl -> as.tags -> gt:::as.tags.gt_tbl -> random_id() = paste(sample(letters, 10, replace = TRUE)), installed gt 1.4.0. It draws because table_id is NA: R/ sets no gt table id (no table.id or id = anywhere). gt construction does not draw: create_unique_id_vals draws only when na_is_index = FALSE, and no caller uses that.
    - fixed_seed_call: forestsearch() at forestsearch-2024.Rmd:301 -> R/forestsearch_main.R:3669 do.call(subgroup.consistency, ...) with no seed -> R/subgroup_consistency_main.R:498 set.seed(8316951) (formal at :397), in the main process.
    - seeded: No. gt does not seed itself, the file has no set.seed, and sg_tables (:415), called between the two, neither seeds nor draws: its static closure (30 functions) has no RNG primitive. Only pkgdown's set.seed(1014) before render fixes the stream, and that comes from the build tool.
    - displayed_outputs: Chunk fs-results (:412-433), the 'Identified Subgroups' gt table. The draw sets only the table's random 10-letter HTML container id (<div id=...>) and the CSS scoped to #id. The visible cells (M.1, M.2, N, E, d1, hr, Pcons and the footnotes) do not depend on it, so no shown number or figure changes. After the change, by the document alone (no set.seed anywhere) the id would vary from one render to the next under a plain rmarkdown::render or Knit, as the :165 baseline table's id already does. Under pkgdown (set.seed(1014) before render) it changes once relative to the old render and then stays stable. Today it is already stable under both routes.
    - built_by: pkgdown only. .Rbuildignore:20 '^vignettes/articles$' keeps it out of the tarball, so it is not a CRAN vignette despite DESCRIPTION:30 'VignetteBuilder: knitr, quarto'. _pkgdown.yml:67-68 (navbar) and :611-613 (articles index) list it, and .github/workflows/pkgdown.yaml:40 runs pkgdown::build_site_github_pages(new_process = FALSE, install = FALSE) against the installed forestsearch (cran::forestsearch at :36).
  - Row 2:
    - file_line: forestsearch-2024.Rmd:439
    - later_call: knitr auto-print of res_tabs$tab_estimates, a gt_tbl built at R/summary_utility_functions.R:1045 (gt::gt(tab_est, auto_align = TRUE)). Same chain: gt:::as.tags.gt_tbl -> random_id() -> sample(letters, 10, replace = TRUE).
    - fixed_seed_call: forestsearch() at forestsearch-2024.Rmd:301 -> forestsearch_main.R:3669 -> subgroup_consistency_main.R:498 set.seed(8316951).
    - seeded: No. There is no set.seed or seed argument between :301 and :439. The only call in between is the gt print at :432, which draws but does not seed.
    - displayed_outputs: Chunk fs-estimates (:437-440), the 'Treatment Effect Estimates' gt table. Only its HTML container id and scoped CSS change; the visible HR, CI and count cells do not. After the change it would vary from render to render by the document alone. Under pkgdown (set.seed(1014)) it changes once and then stays stable.
    - built_by: pkgdown only (.Rbuildignore:20 ^vignettes/articles$; _pkgdown.yml:67-68, :611-613; .github/workflows/pkgdown.yaml:40).
  - Row 3:
    - file_line: forestsearch-2024.Rmd:510
    - later_call: knitr auto-print of summaries$table, a gt_tbl from format_bootstrap_table (R/bootstrap_summaries_helpers.R:157-159, or the placeholder gt at :88 when there are no estimates), called at R/summarize_bootstrap_results.R:169. Chain: gt:::as.tags.gt_tbl -> random_id() -> sample().
    - fixed_seed_call: forestsearch_bootstrap_dofuture() at forestsearch-2024.Rmd:478. Path to the fixed seed: bootstrap_dofuture_main.R:484 -> bootstrap_analysis_dofuture.R:408-415 %dofuture% -> :614 do.call(forestsearch, ...) -> subgroup_consistency_main.R:498, in the workers. Its main-process stream is reset by set.seed(8316951L) at bootstrap_calculations_helpers.R:27 and bootstrap_analysis_dofuture.R:334, then moved forward 2 steps by the %dofuture% seed derivation. The earlier fixed-seed call forestsearch() :301 also precedes it.
    - seeded: No. The seed = 8316951L at :481 belongs to the fixed-seed call, not to the later call, and after the change the bootstrap's restore on exit undoes its effect on the main stream. summarize_bootstrap_results (:501) neither seeds nor draws: its closure has no RNG primitive, its plots use only geom_histogram and geom_vline, and factor presence is built with as_gt = FALSE.
    - displayed_outputs: Chunk bootstrap-summary (:499-511), the bias-corrected estimates gt table. Only the HTML container id and scoped CSS change; the visible bias-corrected HRs do not. Old: the id comes from the 8316951L stream the bootstrap leaves. New (wrapper restores on exit): it comes from the document's own stream, so it varies render to render by the document alone. Under pkgdown (set.seed(1014)) it changes once and then stays stable.
    - built_by: pkgdown only (.Rbuildignore:20 ^vignettes/articles$; _pkgdown.yml:67-68, :611-613; .github/workflows/pkgdown.yaml:40).
  - Row 4:
    - file_line: forestsearch-2024.Rmd:564
    - later_call: knitr auto-print of summaries$diagnostics_table_gt, a gt_tbl from format_bootstrap_diagnostics_table (R/bootstrap_summaries_helpers.R:630-632), called at R/summarize_bootstrap_results.R:283. Chain: gt:::as.tags.gt_tbl -> random_id() -> sample().
    - fixed_seed_call: forestsearch_bootstrap_dofuture() at forestsearch-2024.Rmd:478 (same path as the :510 row; main-process reset at bootstrap_calculations_helpers.R:27 and bootstrap_analysis_dofuture.R:334).
    - seeded: No. Between :478 and :564 run plot_sg_weighted_km (:518), figure_note (:541 inline) and summarize_bootstrap_events (:554). plot_sg_weighted_km calls weightedsurv::df_counting (R/plot_sg_weighted_km.R:345, :354) with default draws = 0, so it takes no set.seed branch, and weightedsurv::plot_weighted_km (:369, :386) has no RNG. The other two have no RNG and no seed; summarize_bootstrap_events only cat()s.
    - displayed_outputs: Chunk bootstrap-diagnostics (:562-565), the 'Bootstrap Diagnostics Summary' gt table. Only the HTML container id and scoped CSS change; the visible diagnostics do not. After the change it would vary render to render by the document alone. Under pkgdown (set.seed(1014)) it changes once and then stays stable.
    - built_by: pkgdown only (.Rbuildignore:20 ^vignettes/articles$; _pkgdown.yml:67-68, :611-613; .github/workflows/pkgdown.yaml:40).
  - Row 5:
    - file_line: forestsearch-2024.Rmd:632
    - later_call: knitr auto-print of metrics_tables, a gt_tbl from cv_metrics_tables with table_style 'combined' and use_gt TRUE (R/cv_summary_tables.R:459-460). Chain: gt:::as.tags.gt_tbl -> random_id() -> sample().
    - fixed_seed_call: forestsearch_tenfold() at forestsearch-2024.Rmd:606. Path to the fixed seed: R/forestsearch_cross_validation.R:951-955 %dofuture% with seed = TRUE -> :1002 do.call(forestsearch, cv_args) -> subgroup_consistency_main.R:498, in the workers. In the main process it only reads the caller's stream to derive the worker streams and moves it forward 2 sample.int(1) steps (future:::make_rng_seeds and doFuture:::doFuture2). So in the old code the main stream here is still the one fixed by forestsearch_bootstrap_dofuture :478 (set.seed(8316951L)), plus the gt draws at :510 and :564 and tenfold's 2-step advance.
    - seeded: No. The seed = 8316951L at :611 is the fixed-seed call's own argument, used only inside the foreach body (:959). cv_metrics_tables (:622) neither seeds nor draws.
    - displayed_outputs: Chunk kfold-cv (:598-634), the 'Cross-Validation Metrics' gt table. Only the HTML container id and scoped CSS change; the visible CV metric values do not, because they come from set.seed(seed + 1000*ksim) inside tenfold, applied before the body's first draw. After the change it would vary render to render by the document alone. Under pkgdown (set.seed(1014)) it changes once and then stays stable.
    - built_by: pkgdown only (.Rbuildignore:20 ^vignettes/articles$; _pkgdown.yml:67-68, :611-613; .github/workflows/pkgdown.yaml:40).
- **Notes:** The only draws this change can affect are five gt table auto-prints (rows above). They set only the hidden random 10-letter HTML id and its scoped CSS. No visible number, table cell, figure or printed value in the article depends on them. The investigator's analysis is confirmed apart from the line-number, column-list and step-count corrections. Verifier scratch work is in /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2ro/fs2024v/ (closure.R, closure2.R, dofut.R, plan.R). No repository files were touched.

**Not rows:**
- **plot_km_band_forestsearch (:812)** draws: KM_diff calls rnorm(draws.band*n0) and rnorm(draws.band*n1) with draws_band = 20. It seeds itself first: R/plot_km_band_forestsearch.R:349 passes no seedstart, so set.seed(8316951) runs in KM_diff right before rnorm. The draw branch of calculate_risk_event_counts is not taken, because the weights are all 1 (weight.name NULL), and that branch also seeds first. It is not in the 23-name list.
- **timing gt (:914)** draws an id, but the set.seed(8316951) inside KM_diff at :812 runs between :606 and :914. So its id is the same in old and new code, unless the planned wrapper also covers plot_km_band_forestsearch. Its visible timings differ on every render anyway.
- **forestsearch_bootstrap_dofuture (:478)** does not draw from the caller's stream before its own set.seed(8316951L) (bootstrap_calculations_helpers.R:27; first draw sample.int at :29). setup_parallel_SGcons (bootstrap_dofuture_main.R:372) runs earlier but only sets the future plan, and worker launch uses parallelly's stealth_sample, which preserves .Random.seed.
- **forestsearch_tenfold (:606)** does draw from the caller's stream in the main process, with no seed in between: the %dofuture% seed = TRUE derivation moves it forward 2 steps. Its seed argument covers only the workers. Even so, set.seed(seed + 1000*ksim) at forestsearch_cross_validation.R:959 comes before the body's first draw at :960, so fs_kfold and the CV metrics do not depend on the caller's stream. After the change, tenfold's own restore on exit undoes the advance.
- **Other later calls** (sg_tables, summarize_bootstrap_results, plot_sg_weighted_km, figure_note, summarize_bootstrap_events, cv_metrics_tables, create_forest_theme, plot_subgroup_results_forestplot, render_forestplot) have no RNG primitive anywhere in their static call closure over R/. Their external calls (gt construction, ggplot2, forestploter, patchwork, survival::coxph, MASS::glm.nb, grid, weightedsurv with draws = 0) do not draw.
- **Prints that draw nothing:** :574 (original_agreement, a data.table) and :579 (factor_presence, a data.frame, because summarize_factor_presence_robust is called with as_gt = FALSE). The ggplot/patchwork print at :587 and plan('sequential') at :398/:490/:618 also draw nothing; plan switching was checked empirically.
- **gt** draws only when rendering HTML (as.tags.gt_tbl), never when building a table.

**Earlier draws:**
- Before :301, the only draw from the document's stream is the gt print at :165 (create_summary_table returns a gt_tbl).
- The DiagrammeR htmlwidgets at :257-258 save and restore the global seed (htmlwidgets 1.6.4 createWidgetId).
- grf.subg.harm.survival (:231) passes seedit to grf.
- create_summary_table's chisq.test and fisher.test run without simulate.p.value, so they do not draw.

**Render-to-render behaviour:**
- The task's rule (vary if the document never calls set.seed) gives "vary every render" for all five rows, and that holds for a plain rmarkdown::render or RStudio Knit.
- The real build route is pkgdown only. pkgdown 2.2.1 calls set.seed(1014) before rendering: build_articles -> build_article(seed = 1014L, new_process default TRUE) -> build_rmarkdown_article -> callr::r_safe(rmarkdown_render_with_seed). So on the site the five ids change once relative to the old render and then stay stable.
- Today the post-:301 ids are already stable under both routes, because of the 8316951 / 8316951L resets.
- Caveat: the CI job installs cran::forestsearch (.github/workflows/pkgdown.yaml:36) and builds with install = FALSE (:40), and the article render loads the installed package without load_all. So the deployed page changes only once a released forestsearch that contains the wrapper is installed, not at the c40d5309 source. A local pkgdown build with the dev version installed shows the change immediately.
- **Verifier (corrected) corrections:**
  - Row :432: the investigator cited the wrong gt::gt branch for sg10_out. The article calls forestsearch with maxk = 2 (:351), so sg_tables (args_fs <- fs$args_call_all at R/summary_utility_functions.R:937) takes the maxk == 2 branch at :1193 and builds sg10_out at R/summary_utility_functions.R:1196-1197. The cited :1165-1166 is the maxk == 1 branch. The visible columns are M.1, M.2, N, E, d1, hr, Pcons (:1194), not M.1, N, E, d1, hr, Pcons.
  - Fixed-seed call :478: the foreach ... %dofuture% in bootstrap_results is at R/bootstrap_analysis_dofuture.R:408-415 (foreach at :408, seed = TRUE at :411, %dofuture% at :415), not :409-416.
  - Fixed-seed call :301: the guard in forestsearch_bootstrap_dofuture that stops when no subgroup was found is at R/bootstrap_dofuture_main.R:244-251 (if at :244, stop at :246), not :245-251.
  - Fixed-seed call :606 and row :632: foreach(.options.future = list(seed = TRUE)) %dofuture% moves the caller's stream forward TWO sample.int(1) steps, not one. future 1.80.0 make_rng_seeds calls next_random_seed() once after as_lecyer_cmrg_seed(TRUE); that function's RNGkind switch is undone by its own on.exit. doFuture 1.4.0 doFuture2 then calls next_random_seed() again and keeps the advanced state. An Rscript check under both plan(sequential) and plan(multisession, workers = 2) showed exactly 2 steps (scratch dofut.R). The same 2-step advance follows set.seed(8316951L) at R/bootstrap_analysis_dofuture.R:334 in the bootstrap's main process. All conclusions are unchanged, because those positions are fixed in the old code.
  - Not-rows refinement for forestsearch_tenfold (:606): its seed argument does not cover the main process. There it derives the L'Ecuyer worker streams from the caller's stream and moves it forward 2 steps, with no seed in between. In the old code that stream was left by :478's set.seed(8316951L); in the new code it is the document's own stream. So it is an unseeded draw after the fixed-seed calls :301 and :478. It still does not qualify as a row with a displayed effect: set.seed(seed + 1000 * ksim) at R/forestsearch_cross_validation.R:959 overwrites the worker stream before the body's first draw at :960, and after the change tenfold's own restore on exit undoes the main-stream advance.
  - Build-route caveat (added to notes): .github/workflows/pkgdown.yaml:36 installs cran::forestsearch, and :40 runs build_site_github_pages(new_process = FALSE, install = FALSE). Articles render through build_articles -> build_article(seed = 1014L, new_process left at its TRUE default) -> build_rmarkdown_article -> callr::r_safe(rmarkdown_render_with_seed), which loads the INSTALLED forestsearch with no load_all. So on the CI-deployed site the 'change once, then stable' happens when the released package that contains the wrapper is installed, not at the c40d5309 source. A local pkgdown build with the dev version installed shows it immediately. _pkgdown.yml also lists the article in the articles index at :611-613 as well as in the navbar at :67-68.

### vignettes/articles/guohe-nested-continuum.qmd

- **Evaluation:** Every chunk is evaluated.
- Configuration: the YAML header (:1-18) has no execute or eval option. There is no _quarto.yml or _metadata.yml in vignettes/articles/, vignettes/ or the package root. No chunk option sets eval, include or purl; the only #| options are message, warning, echo, code-fold, code-summary and fig-*. knitr::opts_chunk$set at :48 sets only comment and fig.align. The ```r block at :895-901 is plain markdown, so it never runs.
- Nothing is precomputed or read from files. reference-data (:527-617) is hard-coded data.frames.
- Inline r at :99-100, :203, :381-385, :389-392, :398 and :783-784 makes no draws.
- Evaluation order:
  - settings :23
  - setup :45, with RNGkind("Mersenne-Twister", "Inversion", "Rejection") at :53
  - design :127 (definitions only)
  - design-one :179: set.seed(1) at :182, then sim_data draws at :183
  - truth :207: truth_curve per cell at :232-233, which calls set.seed(seed) at :211
  - bounds :284 (definitions only)
  - one-replicate :344: bounds_one at :350. The truth_curve at :348 is skipped, because "0.3" is in names(truth) (cells = c(0, 0.3) at :28).
  - one-reading :376
  - montecarlo :404: plan at :416; future_lapply with future.seed = TRUE at :409-411, once per cell from :417
  - summaries :424, the own-* tables, figure-lower :619, figure-complement :660, timing :693, the ref-* tables, session :937
- No chunk after montecarlo draws random numbers: quantile, median, sd and the base-graphics calls are deterministic.
- **Calls reaching the fixed seed:**
- **set.seed() calls and seed arguments:**
  - :182 set.seed(1) (chunk design-one)
  - :211 set.seed(seed) inside truth_curve(). It runs for each cell from :232-233 with seed = 20260721L + round(100*b2), set at :233.
  - :297 set.seed(base + m) inside bounds_one(). It runs at :350 (one-replicate, main session, Mersenne-Twister) and in every future at :410 (montecarlo, where the kind is L'Ecuyer-CMRG under future.seed = TRUE).
  - :310 seed = base + m + 500000L, passed to guohe_algorithm3(). The function seeds at R/guohe_algorithm3.R:395 before its first draw (sample.int at :396); the full-sample fits at :369 draw nothing.
  - :315 seed = base + m + 700000L, passed to fs_mr_inference(). The function seeds at R/fs_mr_inference.R:615 before any draw (.fs_mr_multipliers rpois, called at :702). It seeds again with seed + 900000L at :912 before the field rnorm draws at :919-920.
  - :348 seed = 20260721L + 30L, passed to truth_curve(). It is in the file but does not run at these settings: the if at :347 takes the truth[["0.3"]] branch.
  - Not a fixed seed: :411 future.seed = TRUE. The future streams come from the caller's current .Random.seed: future::make_rng_seeds calls rng_config$as_seed(TRUE), which switches RNGkind to L'Ecuyer-CMRG and seeds it from the caller's state, then calls next_random_seed(), which moves the caller forward with sample.int(1L, 1L).
- **Rows:** 0
- **Notes:** 1. No place in this file matches, so rows = []. No evaluated call reaches set.seed(8316951) (R/subgroup_consistency_main.R:498, seed formal at :397).

2. The file calls two forestsearch functions, both inside bounds_one() (run at :350 and in every future at :410): guohe_algorithm3() (:307) and fs_mr_inference() (:311).
- guohe_algorithm3 is not in reach_fixed_seed.txt. Its only namespace callees are .g3_coef and .g3_cox_lean. A codetools check of t2_tree/R finds no route to subgroup.consistency, whether by calls or by symbol references.
- fs_mr_inference is in the 23-name list only because of a symbol reference, not a call: digits <- eval(formals(subgroup.consistency)$pconsistency.digits) at R/fs_mr_inference.R:673. That line reads the formal default pconsistency.digits = 2 (R/subgroup_consistency_main.R:396) and never calls the function.
- Line :673 also sits in the .has_effect && .has_cons branch (:665). That branch is not taken here, because the article passes admission = list(effect_floor = NULL, consistency = NULL) (qmd:313).
- Its callees, followed by calls only, do not reach subgroup.consistency (codetools, FALSE): .fs_mr_assemble, .fs_mr_pieces -> .consistency_cox_pieces, .fs_mr_multipliers, .fs_mr_select, .fs_mr_ij_var, .fs_mr_se_from_ij, .fs_mr_ij_floor, .fs_mr_field_complement -> .fs_mr_field_joint, plus .fs_decl_*, .fs_mr_field_recovery and fs_mr_field_uniform. The last three run only when keep_declaration_field, field_recovery or field_uniform is TRUE (:1107, :1058, :1016), and all three are FALSE here.
- Its only symbol-only reference to a namespace function is subgroup.consistency (:673). fs_mr_field_uniform is called at :1020 and also referenced in formals() at :1028.
- So the static edge from fs_mr_inference to subgroup.consistency is spurious for every argument set. This may matter for other files that call fs_mr_inference directly.

3. Even if fs_mr_inference counted as reaching, no row would follow.
- Within a replicate, every draw is explicitly seeded: :297 for the data, guohe :395 for the bootstrap indices, MR :615 for the multipliers, and MR :912 for the field. .fs_mr_field_complement draws nothing.
- After :350, the next use of the caller's state is future_lapply(future.seed = TRUE) at :409-411 (cell 0.0, from :417). plan(multisession) at :416 does not touch .Random.seed (verified).
- The per-element L'Ecuyer-CMRG streams come from the caller's state, but each element first runs set.seed(base + m) at :297 and the package functions then reseed. Inside a future the kind is L'Ecuyer-CMRG, so those set.seed calls seed L'Ecuyer-CMRG. They still fix the stream completely: in a test, inner values were identical for different caller seeds and identical under sequential and multisession plans.
- So no displayed value depends on the future streams: the walk table (:371-373), its inline text (:381-392), own-select-table :476, own-lower-table :494, the two own-complement-table kables :511, figure-lower :619 and figure-complement :660.
- set.seed(1) at :182 and the truth seeds at :211 already fix the document's stream before :350.
- Under the planned save/restore, the caller state after :350 changes at most once, depending on which functions are wrapped:
  - only fs_mr_inference wrapped: guohe_algorithm3's post-bootstrap state;
  - guohe_algorithm3 wrapped too: the state after set.seed(base+1) and sim_data;
  - fs_mr_inference left unwrapped: no change.
  Every case is fixed by explicit seeds and only moves the future streams, so no displayed output changes, either once or from render to render.

4. Package functions that seed themselves (listed as instructed, not as rows):
- guohe_algorithm3: seed = NULL by default (R/guohe_algorithm3.R:243), but the article passes it (:310). It seeds at :395 before sample.int at :396. parallel = FALSE by default (:246), so it uses lapply. The future.seed = FALSE path at :412-414 is not used.
- fs_mr_inference: seed = NULL by default (:595), but the article passes it (:315).

5. RNGkind("Mersenne-Twister", "Inversion", "Rejection") at :53 is not a seed. It comes before every seeded call.

6. Output that varies between renders anyway is wall-clock only: one_sec in the :372-373 caption, cell_sec in the timing table (:694-699), and truth_sec in the captions at :243-244 and :698-699.

7. Aside (outside the fixed-seed question): because the futures run under L'Ecuyer-CMRG, the montecarlo replicates are not the same random streams as a Mersenne-Twister sequential run with the same set.seed values. This includes the :350 replicate. This bears on the article's claim at :103-104 and :281-282 that runs reproduce the manuscript's records, not on RNG leakage.

8. Build: pkgdown only.
- .Rbuildignore:20 ^vignettes/articles$ keeps it out of the CRAN tarball; otherwise DESCRIPTION:30 VignetteBuilder: knitr, quarto would apply.
- _pkgdown.yml lists it at :47 (href articles/guohe-nested-continuum.html) and :583.
- The file is tracked: git ls-files returns it, and the c40d5309 blob is 8f9e7001, the same as t2_tree.

Scratch files: /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2ro/guohe-nested-continuum-verify/ (reach2.R, fut.R, fut2.R, fut3.R).
- **Verifier (corrected) corrections:**
  - The main finding holds. No evaluated call in the file reaches set.seed(8316951) (fixed_seed_calls = []), and no row applies (rows = []). All cited line numbers check out. Qmd: :23, :45, :53, :182, :211, :232-233, :297-298, :307-317, :347-350, :409-411, :416-418. Package source: R/subgroup_consistency_main.R:396, :397, :497-498; R/forestsearch_main.R:3669; R/fs_mr_inference.R:595, :615, :665, :673, :702, :912, :919-920; R/guohe_algorithm3.R:243, :395, :396, :414. Build lines: .Rbuildignore:20, DESCRIPTION:30, _pkgdown.yml:47 and :583. The corrections below refine wording in the notes only. None of them changes a row or a verdict.
  - Note 1 overstates how fs_mr_field_uniform is used. It is not a symbol-only reference: fs_mr_inference calls it at R/fs_mr_inference.R:1020, inside if (isTRUE(field_uniform)) at :1016. That branch is not taken here (field_uniform defaults to FALSE, :599). It is also referenced inside formals() at :1028. The only symbol-only (non-call) reference to a namespace function is subgroup.consistency at :673. A codetools re-check finds the same. fs_mr_inference's call-only callees are .fs_decl_c0_cmp, .fs_decl_c0_shift, .fs_decl_field, .fs_mr_assemble, .fs_mr_field_complement, .fs_mr_field_recovery, .fs_mr_ij_floor, .fs_mr_ij_var, .fs_mr_multipliers, .fs_mr_pieces, .fs_mr_se_from_ij, .fs_mr_select, .fs_pcons_eff and fs_mr_field_uniform. Following calls only, it does not reach subgroup.consistency. guohe_algorithm3's only namespace callees are .g3_coef and .g3_cox_lean.
  - The future.seed = TRUE mechanism needed more detail (future 1.80.0, future.apply 1.21.0). The per-element streams come from rng_config$as_seed(TRUE), which switches RNGkind to L'Ecuyer-CMRG. That switch seeds the new generator from one unif_rand of the caller's Mersenne-Twister state, and the caller's kind and state are then restored. Next, next_random_seed() moves the caller forward with sample.int(1L, 1L). Inside every future the RNG kind is therefore L'Ecuyer-CMRG. So set.seed(base + m) at qmd:297 and the package reseeds (guohe R/guohe_algorithm3.R:395; MR R/fs_mr_inference.R:615, :912) seed L'Ecuyer-CMRG there, whereas the one-replicate run at :350 seeds Mersenne-Twister in the main session. The conclusion is unchanged: within each future, every draw is fully set by those explicit seeds. I checked this by running it: inner values were identical for caller seeds 1 and 777, and identical under plan(sequential) and plan(multisession). plan(multisession) at :416 leaves the caller's .Random.seed unchanged.
  - The vary-or-change-once statement under the planned change needed a precise scope. It assumed the save/restore wraps fs_mr_inference but not guohe_algorithm3. In that case the caller state after :350 becomes fs_mr_inference's entry state: guohe_algorithm3's post-bootstrap state (set.seed(base+1+500000) at guohe :395, then B = 2000 sample.int draws at :396). If guohe_algorithm3 is wrapped too, the state is the one after set.seed(base+1) at :297 and the sim_data draws at :298. If fs_mr_inference is left unwrapped because its edge is spurious, nothing changes at all. Every case is determined by explicit seeds, so the state changes at most once, it is stable from render to render, and it only moves the future seed streams, which no displayed value reads.

### vignettes/articles/methodology.qmd

- **Evaluation:** The YAML header (:1-27) sets no execute, eval, freeze or cache option. There is no _quarto.yml or _metadata.yml in vignettes/articles/, vignettes/ or the repo root. The only _quarto.yml is quarto/_quarto.yml, a sibling directory, so it does not apply to this file. Knitr defaults apply, so every {r} chunk is evaluated.

There are seven {r} chunks:
- setup :29-52 (include=FALSE but evaluated; it only calls knitr::opts_chunk$set and knitr::knit_hooks$set, and the hook body uses knitr::knit_code$get and htmltools::htmlEscape)
- power-approximation :407-456 (echo=FALSE, code_fold=TRUE)
- cv-metrics-table :767-789
- sim-models-table :819-839
- sim-results-table :868-893
- estimation-table :917-945
- session-info :1364-1366

The four table chunks build hard-coded data.frames and pass them to knitr::kable (:788, :832, :881, :934).

Other fenced blocks are not R. The two {dot} Graphviz chunks are at :353-368 and :645-694. The two forestsearch(...) argument listings are in plain ```r fences with no braces (:1059-1159 and :1168-1268), so they are displayed and never executed.

The file has no inline `r ` expressions, never calls library(forestsearch) (the only library() call is library(ggplot2) at :437), and reads no saved results (no readRDS, load, or system.file). The package has no .onLoad or .onAttach hook in R/, so having the namespace loaded during a pkgdown build runs nothing either.
- **Calls reaching the fixed seed:**
- **set.seed() calls and seed arguments:**
  - methodology.qmd:418: set.seed(123), inside integrate_2d() within approx_power() in chunk power-approximation (:407-456). It runs once per grid point through the mapply call at :434. That is 153 calls (hr_seq = seq(0.5, 3.0, by = 0.05) gives 51 values, times n_values = c(60, 80, 100) at :431-432). Each call is followed by the rnorm draws at :420-421.
- **Rows:** 0
- **Notes:** There are no rows because no evaluated call in this file reaches the fixed seed, so no random draw can follow one.

1) Fixed-seed-reaching names. I searched for each of the 23 names in reach_fixed_seed.txt followed by '('. Eight lines match, and none is evaluated:
- forestsearch() at :1060, :1063, :1169 and :1172. :1060 and :1169 are comments; :1063 and :1172 are the calls. All four are inside the non-executed plain ```r fences that open at :1059 and :1168. The seedit = 8316951 lines at :1084 and :1193 belong to those listings.
- forestsearch() at :66, fs_declaration_calibration() at :161 and :489, and forestsearch_bootstrap_dofuture() at :636. All four are backtick code in blockquoted prose, not inline R.
None of the other 19 names appears with '('. For reference, subgroup.consistency sets set.seed(seed) at R/subgroup_consistency_main.R:498 (formal seed = 8316951 at :397), and forestsearch() reaches it through do.call(subgroup.consistency, sc_filtered_args) at R/forestsearch_main.R:3669. This file never invokes that path.

2) The only random draws are rnorm at :420 and :421 in chunk power-approximation. Each call to integrate_2d() runs set.seed(123) at :418 immediately before them, so they are seeded. Nothing from forestsearch runs before them anyway: the setup chunk calls only knitr and htmltools. These draws produce the ggplot figure with labs title 'ForestSearch Power Approximation' (:449) and fig.cap 'Approximate probability of finding H via ForestSearch' (:407). The figure is deterministic now and stays deterministic after the planned change.

3) The other evaluated code draws nothing:
- The four table chunks (:767, :819, :868, :917) only use data.frame and knitr::kable.
- The ggplot at :438-455 uses geom_line, geom_hline, geom_vline, scale_*, labs, theme and annotate, with no jitter or sampling. scales::percent at :444 draws nothing.
- sessionInfo() at :1365 draws nothing.
The RNG state that power-approximation leaves behind (set.seed(123) followed by 2 x 100000 normal draws on the 153rd call) is never used later in the document.

4) Build route: pkgdown only. .Rbuildignore:20 has ^vignettes/articles$, so R CMD build drops the file, even though its YAML declares %\VignetteEngine{quarto::html} (:25) and DESCRIPTION:30 has VignetteBuilder: knitr, quarto. pkgdown builds it through _pkgdown.yml:574-577 (title: Methodology, contents 'articles/methodology') and the navbar entry at _pkgdown.yml:42-43 (href articles/methodology.html). No pre-rendered methodology.html sits in vignettes/articles/.

5) Effect of the planned save/restore of .Random.seed on this file: none. The file never calls forestsearch code while rendering, and its only draws are reseeded locally, so its rendered outputs are unchanged.
- **Verifier (corrected) corrections:**
  - Notes item 1, hit count: the investigator says there are 7 hits but lists only 6 lines (:66, :161, :489, :636, :1063, :1172). Searching for each of the 23 names followed by '(' actually matches 8 lines. The two it missed are the comment lines :1060 and :1169, which read '# forestsearch() arguments of the 2024 workflow ...'. Both sit inside the plain ```r fences (:1059-1159 and :1168-1268), so neither is evaluated, and the conclusion stands.
  - Notes item 2, figure caption: the investigator says the figure 'ForestSearch Power Approximation' has its fig.cap at :407. That title is the ggplot labs() title at :449. The fig.cap at :407 is 'Approximate probability of finding H via ForestSearch'. This is a labelling fix only, and the determinism judgment is unchanged.

### vignettes/articles/mrct-consistency-simulations.qmd

- **Evaluation:** Engine: the file uses knitr ({r} chunks). The YAML header (:1-20) has no execute: options, and the vignette engine is quarto::html (:18). The setup chunk (:22, include=FALSE) is evaluated; its knitr::opts_chunk$set (:23-29) sets no eval or cache option, and :33 sets nsims_demo <- 100L.

Not evaluated: four chunks have eval = FALSE: own-cells :256, own-campaign :316, own-summary :346 and campaign-code :1087. The plain ```r blocks (:65, :112), the ```bash block (:118) and the raw ```{=html} style block (:36-40) are not executed. No chunk uses #| eval: !expr, cache or purl.

Evaluated: every other {r} chunk runs, from packages :42 through demo-summary :1146. No chunk after :1162 runs (the appendix is prose).

Precomputed values are read from the shipped records:
- mrct_record() at :144, :148, :215, :289, :457, :514, :540, :551, :588, :634, :659, :761, :883, :885, :908, :956-957, :1008, :1066 and :1123
- readLines(system.file(...)) at :886

Values computed live:
- the case study, mrct_trial_analysis :153
- the DGM builds at :232, :404, :429, :479, :581, :656 (via :661/:679), :675 and :686
- the analytic layer, which draws nothing
- the 100-replicate demonstration, mrct_region_sims :730
- the case-study redraw, mrct_case_study :779

Inline `r` expressions (:191, :247-251, :293-296, :489, :524-528, :605-606, :645, :668-670, :708-710, :717, :719, :811-813, :842, :870-873, :902-903, :944-945, :989-1001, :1028, :1056, :1074-1079, :1133-1137, :1141-1143, :1164-1167) only format values and draw nothing.
- **Calls reaching the fixed seed:**
  - vignettes/articles/mrct-consistency-simulations.qmd:153 mrct_trial_analysis() (chunk bridge-case). For these arguments it reaches the fixed seed on the main session. Path: R/mrct_trial_analysis.R:120 .mrct_analyze_trial(), then R/mrct_trial_helpers.R:47 do.call(forestsearch, fs_call_args), then R/forestsearch_main.R:3669 do.call(subgroup.consistency, sc_filtered_args) with no seed in consistency_overrides (:3582-3605), then R/subgroup_consistency_main.R:498 set.seed(seed) under formal seed = 8316951 (:397). The trigger is TRUE (qmd:170), so the early return at mrct_trial_analysis.R:107 is not taken, and fs$sg.harm is non-NULL (qmd:174). MR then runs with seed = .g_mr(mr_inference_args$seed, seedit) = seed_rep5 (forestsearch_main.R:3930). It calls set.seed(seed_rep5) at R/fs_mr_inference.R:615 and, with ci_method = 'field', set.seed(seed_rep5 + 900000L) at :912; qmd:174-175 confirms the field ran. field_uniform defaults to FALSE (forestsearch_main.R:3891), so the +910000L seed at fs_mr_inference.R:1029 is not reached. Nothing after the field block re-seeds or restores, and neither forestsearch_main.R nor mrct_trial_analysis.R saves .Random.seed. The call therefore leaves the caller on the fixed field stream (seed_rep5 + 900000L, plus the field's draws).
  - vignettes/articles/mrct-consistency-simulations.qmd:730 mrct_region_sims() (chunk find-sims). It is on the static reach list: forestsearch is called inside the %dofuture% body at R/mrct_simulation.R:575 (tryCatch from :573), under foreach :406-412 with .options.future seed = TRUE (:408) and plan(multisession, workers = workers_demo) (:368). For these arguments set.seed(8316951) runs only inside futures, and it does not leak to the caller. This was checked empirically with future 1.80.0 and doFuture 1.4.0, both with workers = 1 (the sequential fallback) and with workers = 2. On the main session, future:::as_lecyer_cmrg_seed(TRUE) derives the L'Ecuyer worker streams from the caller's state through an RNGkind switch and then restores that state. future:::make_rng_seeds and doFuture:::doFuture2 then each call next_random_seed() (sample.int(1L, 1L)), so the caller's MT stream moves forward by exactly 2 draws.
  - vignettes/articles/mrct-consistency-simulations.qmd:779 mrct_case_study() (chunk find-case, mr_inference = FALSE at :788). It reaches the fixed seed on the main session. Path: R/mrct_case_study.R:301 .mrct_analyze_trial(), then R/mrct_trial_helpers.R:47 forestsearch, then R/forestsearch_main.R:3669, then R/subgroup_consistency_main.R:498. qmd:802 (identical(cs$Hhat, Hhat_def)) shows a declaration, so the stage ran. On exit the function already restores the caller's RNGkind and .Random.seed (R/mrct_case_study.R:262-271), so it leaves no fixed-seed stream behind.
- **set.seed() calls and seed arguments:**
  - The file has no set.seed() call (grep: none). The evaluated seed arguments are:
  - vignettes/articles/mrct-consistency-simulations.qmd:163 seed = seed_rep5, in mr_inference_args of mrct_trial_analysis :153. It leads to R/fs_mr_inference.R:615 and :912.
  - vignettes/articles/mrct-consistency-simulations.qmd:228 seed = 7L, in mrct_build_dgm inside build_own, evaluated by build_own(h = 2) at :232 (observed-region mode, so the region draws at R/mrct_inject_dgm.R:478-479 are skipped). Chain: expand_covariates set.seed(seed) at R/mrct_inject_dgm.R:97, then generate_aft_dgm_flex set.seed(seed) at R/generate_aft_dgm_main.R:248 (via R/mrct_inject_dgm.R:525), then set.seed(seed + 3L) at R/mrct_inject_dgm.R:572, followed by rexp in calculate_hazard_ratios (call at :573, draw at R/generate_aft_dgm_helpers.R:1080).
  - vignettes/articles/mrct-consistency-simulations.qmd:416 seed = 7L, in inject_mrct_structure (chunk reader-dgm :404). It uses set.seed at R/mrct_inject_dgm.R:97, at :277 through the region draws with seed + 1L / + 2L (:478-479), at R/generate_aft_dgm_main.R:248, and at R/mrct_inject_dgm.R:572.
  - vignettes/articles/mrct-consistency-simulations.qmd:473 seed = 7L, in mrct_build_dgm inside build_cell (same chain as :416, with latent_R supplied). build_cell is evaluated at :479 (dgm_S2), :581 (the 5 cells with h != 2), :656 (every p_trig_r evaluation inside ahr_solve at :661 and :679), :675 (dgm_r2) and :686 (S2_D2). :686 is the last build before :730.
  - vignettes/articles/mrct-consistency-simulations.qmd:740 seed = seed_S2, in mrct_region_sims. Each replicate calls simulate_from_dgm(seed = seed + sim) (R/mrct_simulation.R:428), which runs set.seed at R/simulate_from_dgm.R:208 before any draw.
  - vignettes/articles/mrct-consistency-simulations.qmd:780 seed = seed_S2, in mrct_case_study: simulate_from_dgm(seed = seed + k) (R/mrct_case_study.R:279), then R/simulate_from_dgm.R:208, under RNGkind("L'Ecuyer-CMRG") at :280.
  - Not evaluated: :328 seed = 20261102 (own-campaign, eval = FALSE) and :1100 seed = seed (campaign-code, eval = FALSE).
- **Rows:** 0
- **Notes:** There is no such place, so rows = []. No evaluated random draw follows a fixed-seed call without an explicit seed in between, and no displayed output depends on the 8316951 stream.

Evaluated draws after the first fixed-seed call (:153), in evaluation order:

(0) The rest of bridge-case (:164-188) draws nothing. mr_estimates_table (:182; R/mr_estimates_table.R:73) goes through .mr_table_from_fs (:174), which calls SG_tab_estimates with draws = 0.

(1) mrct_build_dgm(seed = 7L) at :232 is seeded. In inject_mrct_structure nothing draws before expand_covariates' set.seed(seed) (R/mrct_inject_dgm.R:97, called at :466). Before that call the code only validates its inputs, runs complete.cases and converts factors (:371-463). Later in the build, generate_aft_dgm_flex re-seeds (R/generate_aft_dgm_main.R:248), and the build ends with set.seed(seed + 3L) (R/mrct_inject_dgm.R:572) followed by the rexp at R/generate_aft_dgm_helpers.R:1080. scale_dgm_hazard draws nothing (R/scale_dgm_hazard.R:42-66).

(2) inject_mrct_structure(seed = 7L) at :404 and every build_cell call (:479, :581, :656, :675, :686) are seeded the same way.

(3) mrct_region_sims at :730 derives its worker streams from the caller's state and moves the caller's MT stream forward by 2 draws (empirical, see fixed_seed_calls). This is not a row, for two reasons:
- The caller's state at :730 does not come from the 8316951 stream, either today or after the planned change. Since :153 the stream has been reset by the seed = 7L builds; the last is build_cell at :686 (S2_D2; S2_D1 reuses dgm_r2), which sets set.seed(10) at R/mrct_inject_dgm.R:572 and then draws rexp.
- No output depends on it. Each replicate seeds itself with simulate_from_dgm(seed = seed_S2 + sim) before any draw (R/mrct_simulation.R:428, R/simulate_from_dgm.R:208), and qmd:763-765 checks the demo against the shipped record to 1e-10.

(4) mrct_case_study at :779 is given seed = seed_S2. Its RNGkind("L'Ecuyer-CMRG") at R/mrct_case_study.R:280 initialises from one uniform of the caller's generator, but set.seed(seed + k) at R/simulate_from_dgm.R:208 overrides that at once, and the caller's state is restored at :262-271.

Nothing evaluated after :779 draws.

These calls draw nothing:
- ahr_population, ahr_limit, ahr_design, ahr_truths_table, ahr_design_table and ahr_solve. They use quadrature, uniroot and a deterministic scan (R/ahr_analytic.R:700-724). pmvnorm with Miwa(steps = 1024) (:660-661) leaves .Random.seed unchanged (checked with mvtnorm 1.4.2, 2-d and 3-d).
- mrct_select_case (R/mrct_case_study.R:72-128) and mrct_core_summary (R/mrct_core_summary.R:73-140).
- get_dfpred, analysis_covariates and mrct_record.
- progressr's progressor (checked: no RNG advance).
- print(cs) at :831 and plot(cs, which = 'forest') at :839, which uses points and segments with no jitter (R/mrct_trial_helpers.R:225-250).

Effect of the planned change. Today the caller is left after :153 on fs_mr_inference's field stream (seed_rep5 + 900000L). After the change it would be restored to the document's state before :153. The file has no set.seed, so .Random.seed may not exist at that point, and the restore must handle that case, as R/mrct_case_study.R:263-270 does. The next draws are seeded by seed = 7L, so no displayed output changes, either once or from render to render. Applying the restore to mrct_region_sims (:730) would only undo its 2-draw advance after the streams are derived, which also changes nothing.

Caveat on scope: if the save and restore were applied to every exported function rather than to the 23 that reach the fixed seed, the set.seed calls inside mrct_build_dgm and inject_mrct_structure would no longer persist after those calls. The stream derivation at :730 would then come from the document's unseeded stream, but for reason (3) the outputs still would not change.

Some outputs vary between renders for reasons unrelated to random numbers:
- timings: :188 seconds_mr, :243, :485, :743, and the inline values at :250, :710, :813, :1142 and :1167
- workers_demo, from detectCores (:728, shown at :743 and :1142)
- live_same at :484 and the inline text at :489, which depend on the platform's latentcor fit

Build: pkgdown only. .Rbuildignore:20 has ^vignettes/articles$, so the article is not in the CRAN tarball. DESCRIPTION:30 (VignetteBuilder: knitr, quarto) covers only vignettes/forestsearch.qmd. _pkgdown.yml:65 (navbar href) and :608 (articles index) list the article, and the rendered vignettes/articles/mrct-consistency-simulations.html is also tracked at c40d5309.

Sources outside the repo: the future and doFuture behaviour was read from the installed future 1.80.0 and doFuture 1.4.0 via Rscript --vanilla, and the RNG advances were measured there:
- as_lecyer_cmrg_seed(TRUE): no net advance
- make_rng_seeds: +1
- a full %dofuture% with seed = TRUE: +2, for workers = 1 and workers = 2
- a set.seed in the body: does not leak
- plan(): no advance
- progressr: no advance
- pmvnorm Miwa: no advance

Scratch scripts are in /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2ro/verify-sims-qmd/. rngtest.R, rngtest2.R, rngtest3.R, rngtest4.R and rngtest5.R read .Random.seed lazily (it was evaluated after the reset), so their RNG-advance numbers are invalid; rngtest6.R, rngtest7.R and rngtest8.R capture it eagerly and are the ones relied on.
- **Verifier (corrected) corrections:**
  - find-sims :730 (future seeding on the main session). The investigator said one sample.int derives the worker streams. In fact the streams come from future:::as_lecyer_cmrg_seed(TRUE): it switches to RNGkind("L'Ecuyer-CMRG"), which seeds the stream from one uniform of the caller's Mersenne-Twister state, and then restores that state. Separately, the main session's stream moves forward by two sample.int(1L, 1L) calls, one in future:::make_rng_seeds and one in doFuture:::doFuture2's own next_random_seed() with on.exit(set_random_seed(oseed)). Tested with future 1.80.0 and doFuture 1.4.0 under Rscript --vanilla: workers = 1 (sequential fallback) and workers = 2 both move the caller forward by exactly 2 MT draws, and a set.seed(8316951) inside the %dofuture% body leaks in neither case. The conclusion is unchanged: no output depends on these draws.
  - Seed chain of the DGM builds (seed = 7L at :228, :416, :473). The investigator left out generate_aft_dgm_flex()'s set.seed(seed) at R/generate_aft_dgm_main.R:248, reached through R/mrct_inject_dgm.R:525 (seed = seed). It runs between expand_covariates' set.seed at R/mrct_inject_dgm.R:97 and the final set.seed(seed + 3L) at :572. The rexp after :572 is in the calculate_hazard_ratios call at R/mrct_inject_dgm.R:573, defined at R/generate_aft_dgm_helpers.R:1076 with the draw at :1080. The 'seeded' judgment is unchanged.
  - The inline `r` list leaves out :717 and :719. Both only format nsims_demo and draw nothing.
  - scale_dgm_hazard spans R/scale_dgm_hazard.R:42-66, not :42-67.
  - The evaluation status omits the raw ```{=html} block at :36-40. It is a style block, not code.

### vignettes/articles/mrct-consistency.qmd

- **Evaluation:** Every R chunk is evaluated. The YAML header (:1-20) has no execute: key, and no _quarto.yml or _metadata.yml exists under vignettes/ (the only _quarto.yml is quarto/_quarto.yml, a separate project). The setup chunk at :22 has include=FALSE, so it runs with its output hidden; it only sets echo, message, warning and fig.* (:23-29). No chunk sets eval or #| eval:, and nothing uses purl or cache (grep for eval|execute|purl|cache finds nothing). The 19 evaluated chunks are: setup :22, packages :38, install-check :90, pre-build :279, pre-table :326, case-read :431, case-itt :450, case-region :464, case-fs :506, case-sets :549, fig-case-forest :581, case-km :590, frontier :629, mr-table :642, case-campaign :696, case-campaign-pipeline :738, one-call :787, window-limits :847 and fig-window :891. The {=html} block at :32-36 is raw HTML. The ```r blocks at :68-82 and :114-118 and the ```bash block at :120-122 have no braces, so they never run. Inline `r` expressions are evaluated: :347-348, :458, :475-476, :534-536, :597-600, :652-660, :770-779, :816 and :936-954. They only format numbers; :816 reports elapsed time. Some values are read from shipped records by mrct_record() instead of being computed: :281 and :433 (rep5_trial.csv), :437 (seed_rep5 from mr_field_S2_rep5.csv), :697 and :699 (core_full_numeric.csv, core_n1000.csv), :710 (entry_cells_core_numeric.csv) and :739 (mr_field_core.csv). The two campaign tables (case-campaign, case-campaign-pipeline) are therefore not computed in the article. These run live: forestsearch(), compute_frontier_cis(), mrct_trial_analysis(), mrct_four_sets() with plot(), mr_estimates_table(), inject_mrct_structure(), ahr_population(), ahr_design(), ahr_limit() and scale_dgm_hazard(). A rendered vignettes/articles/mrct-consistency.html is also tracked. It is a committed render; the .qmd does not read it.
- **Calls reaching the fixed seed:**
  - vignettes/articles/mrct-consistency.qmd:510 forestsearch() (chunk case-fs). Path: forestsearch_main.R:3669 do.call(subgroup.consistency, sc_filtered_args). sc_filtered_args comes from filter_call_args(args_call_all, subgroup.consistency, consistency_overrides) (:3643; summary_utility_functions.R:679-688). args_call_all holds forestsearch's formals (:1884-1885), and the only seed formal there is seedit (:1640), not seed. consistency_overrides (:3582-3640) carries no seed. So subgroup.consistency uses its default seed = 8316951 (subgroup_consistency_main.R:397) and runs set.seed(seed) at :498. The path is reached with these arguments: use_grf, use_lasso and use_dina are all FALSE; is.RCT defaults to TRUE, so ps_method is 'none' (forestsearch_main.R:2602, :2660); a subgroup is declared (stopifnot(!is.null(fs$sg.harm)) at qmd:527). hr.subgroups is therefore non-empty and the early return at subgroup_consistency_main.R:459 is not taken. Later in the same call, MR reseeds with seed = seed_rep5 (qmd:521, forestsearch_main.R:3930): fs_mr_inference.R:615 set.seed(seed), then :912 set.seed(seed + 900000L) for the field method. field_uniform defaults to FALSE (forestsearch_main.R:3891), so the +910000L seed is not used.
  - vignettes/articles/mrct-consistency.qmd:789 mrct_trial_analysis() (chunk one-call). Path: mrct_trial_analysis.R:120 .mrct_analyze_trial(), then mrct_trial_helpers.R:47 do.call(forestsearch, fs_call_args), then forestsearch_main.R:3669 and subgroup_consistency_main.R:498 set.seed(8316951). fs_args carries no seed or seedit. The trigger holds (stopifnot(trigger) at qmd:472, and .mrct_trigger at mrct_trial_helpers.R:118 uses the same Cox fits), so the early return at mrct_trial_analysis.R:107 is not taken. A subgroup is declared (stopifnot(identical(one_call$Hhat, Hhat_def)) at qmd:805). Its MR reseeds with seed = seed_rep5 (qmd:799, fs_mr_inference.R:615, :912).
- **set.seed() calls and seed arguments:**
  - No set.seed() call anywhere in the file (grep finds none).
  - vignettes/articles/mrct-consistency.qmd:296 seed = 7L, argument to inject_mrct_structure() (call opens at :287, chunk pre-build, before the first fixed-seed call). It seeds internally and does not restore: expand_covariates() runs set.seed(seed) at mrct_inject_dgm.R:97 (called at :466). The region draws seeded at :478-479 run only on the planted-region path, which is skipped because region$observed (:469-471). generate_aft_dgm_flex(seed = seed) at :506/:525, then set.seed(seed + 3L) at :572.
  - vignettes/articles/mrct-consistency.qmd:521 seed = seed_rep5 in mr_inference_args of forestsearch() (case-fs). It is used inside the same fixed-seed call, after the fixed seed (fs_mr_inference.R:615, :912).
  - vignettes/articles/mrct-consistency.qmd:632 seed = seed_rep5, argument to compute_frontier_cis() (chunk frontier). frontier_cis.R:164 runs if (!is.null(seed)) set.seed(seed) before the function's only draw, sample() at :272. The early return at :117 is not taken (stopifnot at qmd:634 needs a non-empty frontier). The function does not restore the RNG, so this seed also determines the gt table-id draw at qmd:649.
  - vignettes/articles/mrct-consistency.qmd:799 seed = seed_rep5 in mr_inference_args of mrct_trial_analysis() (one-call). It is used inside the same fixed-seed call, after the fixed seed.
- **Rows:** 0
- **Notes:** No place in this file has an unseeded random draw after a call that reaches the fixed seed, so rows = []. I checked each call after the first fixed-seed call (forestsearch() at qmd:510). A static call-graph scan from fs_in_subgroup, mrct_four_sets, plot.mrct_four_sets, print.mrct_four_sets, plot_pareto_frontier, mr_estimates_table, mr_values, mrct_record, print.mrct_trial_analysis, ahr_population, ahr_limit and scale_dgm_hazard reached 78 package functions, and none contains an RNG call.

(1) case-sets :549-575. fs_in_subgroup (mrct_minor_helpers.R:77-82, get_dfpred) and coxph fits. No draws.

(2) fig-case-forest :585-587. mrct_four_sets() (mrct_four_sets.R:47-64) calls .mrct_four_set_table (mrct_trial_helpers.R:85), which calls .mrct_record_fit (mrct_simulation.R:817-840, coxph only). plot.mrct_four_sets (mrct_four_sets.R:98) calls .mrct_plot_sets, whose forest branch (mrct_trial_helpers.R:228-250) is base graphics. No draws.

(3) case-km :594. The KM branch is mrct_trial_helpers.R:251-279: weightedsurv::df_counting at :263-265 and plot_weighted_km at :266-271 (installed weightedsurv 0.1.0). In the weightedsurv namespace, only calculate_risk_event_counts, cox_rhogamma_resample and KM_diff contain RNG code, and each seeds itself with set.seed(seedstart). df_counting reaches only calculate_risk_event_counts, and with draws = 0 (the default; the article does not pass draws) and unweighted data it takes no RNG branch. plot_weighted_km reaches none of the three. No draws.

(4) frontier :632-635. compute_frontier_cis(seed = seed_rep5) is seeded by its argument before it draws (frontier_cis.R:164 set.seed, sample at :272), so it does not qualify. plot_pareto_frontier (plot_pareto_frontier.R:158-228: geom_step, geom_point, geom_linerange, geom_text, annotate) has no jitter or repel.

(5) mr-table :643-649. The package code (mr_estimates_table, .mr_table_from_fs, format_bootstrap_table, SG_tab_estimates, mr_values) has no RNG calls. The investigator missed one draw here: auto-printing `mr_tab` at :649 goes through gt's knit_print.gt_tbl and as.tags.gt_tbl. The table has no table_id (gt::gt() at bootstrap_summaries_helpers.R:159), so gt:::random_id() runs sample(letters, 10, replace = TRUE) on the global stream. The draw comes after :510, but compute_frontier_cis's set.seed(seed_rep5) at :632 runs in between, so it is seeded, not a row. It changes only the HTML container id and scoped CSS selector (one gt div in the tracked render, id "yydwpabibo"), and the planned save/restore does not change it.

(6) case-campaign and case-campaign-pipeline :696-767. Shipped CSV records, arithmetic and knitr::kable only.

(7) one-call :789-808. mrct_trial_analysis() reaches the fixed seed itself. Before subgroup_consistency_main.R:498 it runs .mrct_check_trial and .mrct_four_set_table (coxph), then forestsearch's pre-consistency steps. Those steps make no draws with these arguments: no GRF, LASSO or DINA (lasso_selection's set.seed at get_FSdata_helpers.R:259 needs use_lasso); is.RCT = TRUE, so no propensity scores; get_FSdata has no RNG calls; and subgroup.search's future.apply::future_lapply at subgroup_search.R:287 has no future.seed, so it derives no streams. Inside the consistency stage, .make_candidate_rng_seeds (subgroup_consistency_main.R:40-66) switches to L'Ecuyer, then restores the RNG kind and state. Even if anything drew earlier, the entry state at :789 is fixed by :632 (and :649), with or without the planned change. print(one_call) (mrct_trial_analysis.R:138-150) is cat() only, and mr_values has no draws.

(8) window-limits and fig-window :847-931. ahr_population, scale_dgm_hazard (scale_dgm_hazard.R:42-66, only rewrites mu) and ahr_limit are deterministic: ahr_analytic.R uses uniroot and Gauss-Legendre via eigen (:38-50). hr_at() is closed-form, and fig-window is base graphics.

ahr_design() at :300 uses mvtnorm::pmvnorm with algorithm = Miwa (ahr_analytic.R:660-661), which is deterministic, and it runs before the fixed seed in any case.

Effect of the planned save/restore: no displayed output changes. Every draw after the fixed-seed calls is one of three kinds: inside a fixed-seed call after its own set.seed; inside compute_frontier_cis after set.seed(seed_rep5); or gt's random_id at :649, also after that seed. The document has no set.seed(), but inject_mrct_structure(seed = 7L) (:287/:296) reseeds the global stream (mrct_inject_dgm.R:97, :572) and does not restore it, so the document's own stream at :510 is fixed in a fresh session.

Today, forestsearch() at :510 leaves .Random.seed in the state MR's seed_rep5-derived draws leave (fs_mr_inference.R:615, :912). Nothing in the article uses that state before :632 reseeds.

Seeded package calls that draw under a default seed argument: none besides those above. compute_frontier_cis's seed defaults to NULL (frontier_cis.R:109) and is given here. weightedsurv::df_counting has seedstart = 8316951 by default, but with draws = 0 it never calls set.seed.

Built by pkgdown only, not CRAN: .Rbuildignore:20 ^vignettes/articles$ leaves the file out of the build tarball, even though DESCRIPTION:30 declares VignetteBuilder: knitr, quarto. _pkgdown.yml lists it at :63 (navbar href articles/mrct-consistency.html) and :607.

The background's "near :3660" for the do.call is :3669 in this tree. set.seed at subgroup_consistency_main.R:498 and seed = 8316951 at :397 are confirmed.
- **Verifier (corrected) corrections:**
  - The investigator missed one random draw, though it does not add a row. The mr-table chunk auto-prints `mr_tab` at qmd:649, a gt table. format_bootstrap_table builds it with gt::gt() and sets no table_id (bootstrap_summaries_helpers.R:159, placeholder branch :88; called from .mr_table_from_fs at mr_estimates_table.R:130). Printing it goes through gt 1.4.0's knit_print.gt_tbl and as.tags.gt_tbl, which call gt:::random_id(), i.e. sample(letters, 10, replace = TRUE), from the global stream. This happens after the fixed-seed call at :510, but compute_frontier_cis(seed = seed_rep5) at qmd:632 reseeds in between (frontier_cis.R:164), so the draw is seeded and is not a row. What it affects is the table's HTML container id and scoped CSS, not visible text; the tracked render has one gt div, id "yydwpabibo". The planned change does not alter it. So the notes' claim that the only draws after the fixed-seed calls are inside fixed-seed calls or inside compute_frontier_cis is incomplete.
  - Line references for plot.mrct_four_sets: the forest-plot code is mrct_trial_helpers.R:228-250, not :239-262. The KM code is :251-279, with weightedsurv::df_counting at :263-265 and plot_weighted_km at :266-271, not :270-281.
  - mrct_four_sets() is mrct_four_sets.R:47-64, not :47-70. plot.mrct_four_sets is at :98.

### vignettes/articles/paper_simulations.qmd

- **Evaluation:** Rendered by Quarto with the knitr engine in one R session. The YAML header (:1-20) has no execute: key. The setup chunk (:23-48, include=FALSE, still evaluated) sets only echo, message, warning and figure size through knitr::opts_chunk$set (:24-30), plus a code_fold hook. The second YAML block (:358-360) holds only a title.

Chunks NOT evaluated:
- #| eval: false: setup-gbsg-dgm-usage :371 (:372), generate-aft-equivalent :395 (:396), simulate-from-either :473 (:474).
- eval=FALSE in the chunk header: twostage-comparison :1486, noise-vars :1519, sens-analysis :1550, save-results :1586, per-sim-timing :1689, complete-example :1728, appendix-table1 :1889, appendix-table2 :1927, multi-scenario-example :2321.
- So the run_simulation_analysis() calls at :1497, :1531, :1565 and :1776 never run.
- The calls at :230 and :278 are in plain ```r display blocks (:223-232, :271-286).
- :139 is a {dot} (Graphviz) chunk, not R.

Nothing is precomputed: no readRDS, load or fread anywhere, so both simulation loops run live.

Evaluated chunks: setup, load-packages :113, timing-init :296, create-dgm-alt :500, access-hrs :530, individual-effects :559, calibrate-k-inter-cox :600, calibrate-k-inter-ahr :631, compare-calibration :662, validate-k-inter :678, create-dgm-null :695, simulate-single :720, sim-individual-effects :749, km-plot :769, summary-table :781, setup-parallel :818, define-params :831, twostage-params :891, run-alt-sims :908, run-null-sims :966, summarize-alt :1024, interpretation-alt :1069, summarize-null :1079, interpretation-null :1125, appendix-table1-earlier :1136, theoretical-detection-rates :1175, format-null-tables :1320, key-metrics :1335, timing-summary :1616, effective-per-sim :1713, ref-table4-nonoise :2405, ref-table4-noise :2492, ref-table5 :2602 (all three ref-table chunks are eval=TRUE with hard-coded data passed to gt), cleanup :2829, session-info :2837.

A mechanical extraction of evaluated chunk code (1180 lines) matched none of the 23 fixed-seed reachers except run_simulation_analysis at :925 and :983.
- **Calls reaching the fixed seed:**
  - paper_simulations.qmd:925 run_simulation_analysis(). This is the body of foreach(..., .options.future = list(seed = TRUE)) %dofuture% at :915-954 (chunk run-alt-sims :908-962).
- Path: if (is.null(method_specs) && run_fs) at R/run_simulation_analysis.R:586 (run_fs = TRUE at :939) -> .run_fs_analysis_gen (:591-599) -> do.call(forestsearch, fs_args) (:819) -> has_subgroups gate (R/forestsearch_main.R:3544-3566: some candidate HR > hr.consistency = 1.0, :3554-3559) -> do.call(subgroup.consistency, sc_filtered_args) (:3669) -> set.seed(8316951) (R/subgroup_consistency_main.R:497-499; formal seed = 8316951 at :397).
- No seed reaches subgroup.consistency: fs_params (:850-879) has none, forestsearch has only seedit (R/forestsearch_main.R:1640), and filter_call_args (R/summary_utility_functions.R:679-688) passes only shared formals.
- Each replicate is a future under plan(multisession, workers = n_workers) (:820-822). So set.seed runs in a worker process, or, when workers == 1, in a sequential future whose RNG state future and doFuture restore.
  - paper_simulations.qmd:983 run_simulation_analysis(). This is the body of foreach(..., .options.future = list(seed = TRUE)) %dofuture% at :972-1012 (chunk run-null-sims :966-1018). Same path and same evaluation context as :925, under dgm_null.
- **set.seed() calls and seed arguments:**
  - No literal set.seed() appears in any evaluated chunk.
  - paper_simulations.qmd:515 seed = 8316951 in setup_gbsg_dgm (call :507-517). It is passed at R/setup_gbsg_dgm.R:115 to .create_gbsg_dgm_, which runs set.seed(seed) at R/sim_aft_gbsg.R:478 and set.seed(seed + 1) at :604.
  - Calls with no seed argument that seed themselves from a default:
- calibrate_k_inter at :604 and :635. Each eff_at call runs .create_gbsg_dgm_ (R/sim_aft_gbsg.R:1024) with default seed = SEED_BASE = 8316951L (:29, :245).
- setup_gbsg_dgm at :614, :645 and :699. Default seed = 8316951L (R/setup_gbsg_dgm.R:97), passed at :115.
- validate_k_inter_effect at :683. It calls .create_gbsg_dgm_ for each k (R/sim_aft_gbsg.R:1254) with the default seed.
  - paper_simulations.qmd:729 seed = 1 in simulate_from_dgm (call :725-732), which runs set.seed(1) at R/simulate_from_dgm.R:208. This is the last seeding of the document session's stream; nothing after it reseeds the session.
  - paper_simulations.qmd:921 and :978 .options.future = list(seed = TRUE). doFuture2 calls make_rng_seeds(nX, seed = TRUE) to derive per-iteration L'Ecuyer-CMRG seeds from the session's current state. This forwards the session stream (next_random_seed in make_rng_seeds and again in doFuture2, deparsed lines 161-162), and the forwarded state is restored on exit.
  - paper_simulations.qmd:947 and :1005 seed_base = 8316951L. Inside each replicate this becomes simulate_from_dgm(seed = seed_base + sim_id) (R/run_simulation_analysis.R:479-485, seed at :482), which runs set.seed at R/simulate_from_dgm.R:208 in the worker.
- **Rows:** 0
- **Notes:** No row qualifies: neither fixed-seed call (:925, :983) ever sets the RNG of the R session that renders the document.

1. Where the fixed seed runs.
- With 2 or more workers it runs in separate worker processes. I verified with Rscript: plan(multisession, workers = 2), then foreach seed = TRUE whose body calls set.seed(8316951). The session's .Random.seed and the next gt:::random_id() were identical with and without that body seed. plan() itself did not move the session RNG. Script: /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2ro/paper_simulations_verify/plan_rng.R
- With workers == 1, future 1.80.0 returns SequentialFutureBackend and the futures run in-session. future's default reset includes "rng", and evalFutureInternal restores .Random.seed. doFuture2 also restores its forwarded oseed on exit, after values <- value(fs).
- So after each loop the session state is a fixed function of its state before the loop, whatever the loop body does. This matches the investigator's leak test: /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2ro/paper_simulations/leak_test.R

2. Evaluated draws after :915 from the session's state. In document order they follow :925 with no seed in between, but they are not rows because the fixed seed never reaches their stream.
- (a) The second foreach (:972, seed = TRUE): make_rng_seeds/next_random_seed.
- (b) Every gt table printed after the loops (missed by the investigator). gt 1.4.0 knit_print.gt_tbl -> as.tags.gt_tbl -> random_id() = sample(letters, 10, replace = TRUE) runs when no table_id is set. I verified that knit_print advances .Random.seed by exactly that sample.
  - No id is set by format_oc_results (R/oc_analyses.R:814), build_classification_table (R/simulation_tables.R:413, returned :462), build_estimation_table (:956, returned :1028) or render_reference_table (:1530, returned :1567), nor by gt(timing_df) :1653 or gt::gt :2728.
  - Lines: :1027, :1052, :1083, :1107, :1156, :1322, :1653, :2485, :2566, :2790.
- All of these draw from the session's own stream: fixed by set.seed(1) at :729 (R/simulate_from_dgm.R:208) and then forwarded deterministically, including twice by each foreach. The planned save/restore cannot change them. The gt draws affect only the random HTML table id and its id-scoped CSS, not any displayed value, and stay stable from render to render.
- The pre-loop gt print at :782 (create_summary_table) also draws. So may fisher.test(simulate.p.value = TRUE) at R/create_summary_table.R:365, only when some cell count is under 5. Both are before :915, outside the question.

3. Later calls checked that do not draw.
- format_oc_results, build_estimation_table, interpret_estimation_table, build_classification_table, render_reference_table: a transitive callee walk of the package R/ finds no RNG calls. The one hit, rnorm in compute_detection_probability_single, sits on the monte_carlo branch only.
- compute_detection_probability(method = "cubature") at :1202 and :1227, and generate_detection_curve (default method = "cubature", R/truefind_asymptotic.R:370), use the deterministic cubature::adaptIntegrate (R/truefind_asymptotic.R:288). rnorm appears only at :309-310.
- Also no draws: plot_detection_curve (:436-491: plot, abline, legend, grid), points/abline/legend (:1293-1316), timing code, plan("sequential") :2830, gc :2831, sessionInfo :2838.

4. Displayed simulation outputs.
- results_alt and results_null feed summarize-alt :1024, interpretation-alt :1069, summarize-null :1079, interpretation-null :1125, appendix-table1-earlier :1136, theoretical-detection-rates :1175 (empirical rates and prop_cens), format-null-tables :1320, key-metrics :1335 and effective-per-sim :1713 (timing only).
- Each replicate reseeds with simulate_from_dgm(seed = 8316951 + sim_id) before drawing. So these outputs do not depend on the session stream or on the L'Ecuyer seeds, and the planned change leaves them unchanged.
- Within a replicate, the GRF arm runs after forestsearch() but seeds grf through seedit = 8316951L (R/run_simulation_analysis.R:86 -> R/grf_helpers.R:75/:85).
- Any other draw after forestsearch() inside the replicate keeps its stream if the restore happens only at the outermost exported function (run_simulation_analysis). After that function returns, the next iteration gets a fresh seed assignment from doFuture.
- n_add_noise = 0L, so the set.seed at R/run_simulation_analysis.R:509 does not run.

5. Build. The article is built by pkgdown only.
- .Rbuildignore:20 is ^vignettes/articles$. DESCRIPTION:30 VignetteBuilder: knitr, quarto does not reach it.
- _pkgdown.yml lists it at :57 and :598.
- .github/workflows/pkgdown.yaml installs cran::forestsearch (:36) and runs build_site_github_pages(new_process = FALSE, install = FALSE) (:40). Quarto articles are rendered by quarto in a separate R process, which likely loads the installed CRAN release, not this source tree.
- pkgdown 2.2.1's build_quarto_articles sets no seed.

6. Side observation. The DGM builders seed with 8316951 (SEED_BASE) in the session. That is ordinary document seeding, not the subgroup.consistency fixed seed, and set.seed(1) at :729 supersedes it before the loops.
- **Verifier (corrected) corrections:**
  - The investigator missed later random draws. Every gt table printed after :915 draws from the document session's RNG. The path is gt 1.4.0 knit_print.gt_tbl -> as.tags.gt_tbl -> random_id(), which is paste(sample(letters, 10, replace = TRUE)) and runs whenever the table has no table_id. None of the package's gt builders set an id: R/oc_analyses.R:814, R/simulation_tables.R:413, :956, :1530 (returned at :462, :1028, :1567). Neither do the file's own gt(timing_df) at :1653 or gt::gt at :2728. Printed, evaluated tables after the loops: :1027, :1052, :1083, :1107, :1156, :1322, :1653, :2485, :2566, :2790. The note 'gt tables ... no draws' is therefore wrong, and so is the claim that the second foreach at :972 is the only later call that draws from the document's state. The conclusion does not change: these draws come from the document's own stream, not the 8316951 stream, and they affect only the random HTML table id. rows stays [].
  - Line number: setup_gbsg_dgm passes seed to .create_gbsg_dgm_ at R/setup_gbsg_dgm.R:115. The investigator gave :106, which is only the line where that call starts.
  - Path condition is more precise than 'replicates whose screen produces candidates'. forestsearch() reaches subgroup.consistency only when has_subgroups is TRUE (R/forestsearch_main.R:3544-3566), meaning some candidate HR > hr.consistency = 1.0 (:3554-3559). A seed cannot be passed through: forestsearch has no seed formal (only seedit, R/forestsearch_main.R:1640), and filter_call_args (R/summary_utility_functions.R:679-688) keeps only names that are formals of subgroup.consistency. The FS-arm block is R/run_simulation_analysis.R:586-609, with the .run_fs_analysis_gen call at :591-599 (the investigator gave :588-601).
  - Clarification of the RNG mechanics, not a reversal. make_rng_seeds forwards the session stream: as_lecyer_cmrg_seed restores the state exactly, then next_random_seed forwards it by one sample.int. doFuture2 forwards it again (deparsed body lines 161-162: oseed <- next_random_seed(); on.exit(set_random_seed(oseed), add = TRUE)). That on.exit restore of the forwarded state runs only after values <- value(fs) (deparsed line 470) inside doFuture2. So after each loop the session's .Random.seed is a fixed function of its state before the loop, for any worker count, including the workers == 1 fallback (future 1.80.0 MultisessionFutureBackend returns SequentialFutureBackend when workers == 1L).

### vignettes/articles/survival-post-selection.qmd

- **Evaluation:** Every chunk is evaluated.

- The YAML header (qmd:1-18) has no `execute:` key; the engine is quarto::html (qmd:16), running knitr.
- No chunk sets `eval: false`. The only chunk options are message/warning (qmd:24-25, 141, 162, 180, 263-264), `results: hide` (qmd:142, 163, 181) and `echo: false` (qmd:339), which still runs.
- The inline `r` expressions at qmd:356, 358, 360, 382 and 429-431 run.
- Nothing is precomputed: there is no readRDS/load. The three fits (qmd:144, 165, 183) and the three bootstraps (qmd:272) are computed at render time.
- table-s2 (qmd:321-336) is a hand-typed data.frame.
- Two chunks autoprint gt HTML tables: boot-summary (qmd:287, `summ_fs$table`) and mr-table (qmd:297, the value of `mr_estimates_table(...)`). Printing them runs gt's `knit_print.gt_tbl` -> `as.tags.gt_tbl` -> `random_id()` -> `sample()`.
- All other output is `cat`/`print` or `knitr::kable`, which do not draw.
- **Calls reaching the fixed seed:**
  - qmd:144 `forestsearch()` [fit-fs, `subgroup_method = "consistency"`]. It calls forestsearch_main.R:3669 `do.call(subgroup.consistency, sc_filtered_args)` with no `seed`: `args_call_all` is `mget(names(formals()))` at :1885, and forestsearch's formal is `seedit` (:1640). So the default `seed = 8316951` (subgroup_consistency_main.R:397) reaches `set.seed(seed)` at :498, in the main session. The same call then re-seeds MR with `set.seed(seedit)` at fs_mr_inference.R:615 (passed from forestsearch_main.R:3930) and `set.seed(seedit + 900000L)` at fs_mr_inference.R:912 (`ci_method = "field"`).
  - qmd:272 `forestsearch_bootstrap_dofuture()`, loop iteration m = "FS". Each replicate runs `do.call(forestsearch, args_FS_boot)` (bootstrap_analysis_dofuture.R:614) on the consistency path, which reaches `set.seed(8316951)` at subgroup_consistency_main.R:498. With `n_workers = min(2L, availableCores())` = 2 (qmd:265-267) this runs inside multisession workers; it runs in the main session only under the sequential fallback. The call also carries its own `seed = 8316951L` (qmd:273).
  - Not reached on these paths: qmd:165 fit-dina returns at forestsearch_main.R:2887 and qmd:183 fit-grf returns at :3072, both before :3669. The edge to `subgroup.consistency` at fs_mr_inference.R:673 only reads `formals(subgroup.consistency)`. The DINA and GRF bootstrap iterations re-run DINA/GRF fits and do not reach it.
- **set.seed() calls and seed arguments:**
  - The document has no `set.seed()` call.
  - qmd:148 `seedit = 8316951` (fit-fs): MR `set.seed(seedit)` at fs_mr_inference.R:615 (via forestsearch_main.R:3930) and `set.seed(seedit + 900000L)` at fs_mr_inference.R:912.
  - qmd:169 `seedit = 8316951` (fit-dina): `.fs_apply_mr(seedit)` at forestsearch_main.R:2881 -> `seed = .g(mr_inference_args$seed, seedit)` at fs_mr_inference_methods.R:172 -> `set.seed` at fs_mr_inference.R:615.
  - qmd:171 `dina_args = list(seed = 8316951, n_folds = 2)`: `.resolve_dina_args` puts it in `fit$seed` (forestsearch_helpers.R:1156) -> `dina()` -> `dina_fit` -> `.dina_make_folds` `set.seed(seed)` at dina.R:433, then `sample.int` at :435. The caller's .Random.seed is saved and restored around it (dina.R:420-431).
  - qmd:187 `seedit = 8316951` (fit-grf): `grf::causal_survival_forest(seed = seedit)` at grf_helpers.R:75. Grf's `runif()` default for `seed` is never evaluated, and its nuisance forests get `seed = seed` (grf 2.6.1). MR `set.seed(seedit)` at fs_mr_inference.R:615 via forestsearch_main.R:3066.
  - qmd:273 `seed = 8316951L` (`forestsearch_bootstrap_dofuture`, run 3 times in the loop): `.make_boot_index_matrix` `set.seed(seed)` at bootstrap_calculations_helpers.R:27 (called from bootstrap_dofuture_main.R:389-393), then `bootstrap_results` `set.seed(seed)` at bootstrap_analysis_dofuture.R:334, then foreach `%dofuture%` with `seed = TRUE` at :411. future 1.80's `as_lecyer_cmrg_seed` derives the streams from that state and restores it.
- **Rows:** 2
  - Row 1:
    - file_line: vignettes/articles/survival-post-selection.qmd:287
    - later_call: Autoprint of `summ_fs$table` (chunk boot-summary). The object is a gt_tbl built by `summarize_bootstrap_results()` (qmd:285) -> `format_bootstrap_table` (summarize_bootstrap_results.R:169, returned as `table` at :411) -> `gt::gt()` at bootstrap_summaries_helpers.R:159, with no table id. knitr prints it through gt 1.4.0 `knit_print.gt_tbl` -> `as.tags.gt_tbl` -> `random_id()` = `sample(letters, 10, replace = TRUE)`, which draws 10 values from the main-session RNG.
    - fixed_seed_call: qmd:144 `forestsearch()` [fit-fs] -> forestsearch_main.R:3669 `do.call(subgroup.consistency, ...)` with no `seed` -> `set.seed(8316951)` at subgroup_consistency_main.R:498, in the main session. Also qmd:272 (bootstrap, FS iteration) -> bootstrap_analysis_dofuture.R:614 `do.call(forestsearch, args_FS_boot)` -> the same `set.seed` (in workers, or the main session under the sequential fallback). In between come fit-dina (qmd:165), fit-grf (qmd:183) and the DINA/GRF bootstrap calls (qmd:272). All are exported forestsearch functions that the planned save/restore would wrap.
    - seeded: No. The draw has no seed of its own: no table id is set (gt::gt() at bootstrap_summaries_helpers.R:159, and no `table.id` anywhere in format_bootstrap_table, summarize_bootstrap_results or mr_estimates_table), and the document has no `set.seed`. Today the stream is fixed only because the last bootstrap call (GRF, `seed = 8316951L` at qmd:273) leaves its own `set.seed` (bootstrap_calculations_helpers.R:27, bootstrap_analysis_dofuture.R:334) in force when it returns. That seed belongs to the bootstrap call, and `forestsearch_bootstrap_dofuture` is exported (NAMESPACE:160) and on the reach list, so the planned save/restore would wrap it.
    - displayed_outputs: boot-summary's gt table (qmd:287). Only its HTML container id (`<div id="xxxxxxxxxx">`) and the CSS selectors in its `<style>` block (`#xxxxxxxxxx ...`) depend on these draws. The visible cells (naive and bootstrap bias-corrected HR, from boots$FS, seeded by `seed = 8316951L`) do not. Today the id is the same on every render, fixed by the GRF bootstrap's leaked `set.seed(8316951L)`. After the change the RNG is back to its state at entry to fit-fs, which no document seed fixes: it is initialised from time and PID, or absent. The id would then change on every render. Nothing on the visible page changes; the rendered HTML source does, so site diffs show it.
    - built_by: pkgdown only, not CRAN. Evidence: .Rbuildignore:20 `^vignettes/articles$` drops it from the tarball; the article note at qmd:20 says the same; _pkgdown.yml:45 and :582 list articles/survival-post-selection. DESCRIPTION:30 `VignetteBuilder: knitr, quarto` covers only the top-level vignettes that ship.
  - Row 2:
    - file_line: vignettes/articles/survival-post-selection.qmd:297
    - later_call: Autoprint of the value of `mr_estimates_table(fit_fs, boot_results = boots$FS, est.scale = "hr", digits = 2)` (chunk mr-table). It returns `format_bootstrap_table(...)` (mr_estimates_table.R:153) -> `gt::gt()` at bootstrap_summaries_helpers.R:159, with no table id. The print path is the same: gt `knit_print.gt_tbl` -> `as.tags.gt_tbl` -> `random_id()` -> `sample(letters, 10, replace = TRUE)`. It continues the main-session stream right after the 10 draws at qmd:287.
    - fixed_seed_call: Same as the qmd:287 row: qmd:144 fit-fs -> forestsearch_main.R:3669 -> `set.seed(8316951)` at subgroup_consistency_main.R:498 (main session), and qmd:272 FS iteration -> bootstrap_analysis_dofuture.R:614 -> the same. The only call in between is boot-summary's gt print (qmd:287), which draws without seeding.
    - seeded: No. The draw has no seed of its own, and the document has no `set.seed`. Today the stream is fixed only by the leaked `set.seed(8316951L)` of the last bootstrap call (bootstrap_calculations_helpers.R:27, bootstrap_analysis_dofuture.R:334), which the planned save/restore would contain.
    - displayed_outputs: mr-table's gt table (qmd:297). Only the HTML container id and its scoped CSS selectors depend on these draws. The visible cells (naive HR, and MR two-term estimates with IJ intervals from `fit_fs$mr_inference`) do not. Today the id is the same on every render; after the change it would change on every render. Nothing on the visible page changes.
    - built_by: pkgdown only, not CRAN: .Rbuildignore:20 `^vignettes/articles$`; qmd:20 article note; _pkgdown.yml:45 and :582; DESCRIPTION:30 `VignetteBuilder: knitr, quarto` covers only the shipped top-level vignettes.
- **Notes:** Two rows: the gt HTML prints at qmd:287 and qmd:297. Every model-based number in the article is unaffected by the planned change.

Evaluation order of the calls that matter:
1. fit-fs (qmd:144) reaches `set.seed(8316951)` in the main session.
2. fit-dina (qmd:165) and fit-grf (qmd:183).
3. The bootstrap loop (qmd:272, three calls, each with `seed = 8316951L`).
4. boot-summary and mr-table print gt tables (qmd:287, 297). These are the rows.
5. The kable tables and inline values follow.

Why the other later calls are not rows (checked against the source; the static draw-reach script is at /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2ro/survpostsel-verify/reach_nostr.R):
- fit-dina: the forestsearch steps before dispatch (forestsearch_main.R:1737-2757; `ps_method` is "none" at :2602) reach no draw. `.forestsearch_dina_select` reaches only `.dina_make_folds`, which seeds at dina.R:433 and restores the caller's state. Propensity is "logistic" (glm). MR is seeded at fs_mr_inference.R:615.
- fit-grf: `causal_survival_forest(seed = seedit)` (grf_helpers.R:75). The policytree and grf scoring calls do not draw from R's RNG. `grf.subg.harm.glm` (`set.seed` at :460) is not on the survival path. MR is seeded at fs_mr_inference.R:615.
- Bootstrap: nothing before bootstrap_dofuture_main.R:389 draws. parallelly's port choice uses `stealth_sample`, which saves and restores .Random.seed. future's `seed = TRUE` restores the caller's state.
- `summarize_bootstrap_results` (`create_plots = FALSE`; the plot branch is at summarize_bootstrap_results.R:302), `mr_estimates_table`, `mr_method_comparison` (through `mr_values`) and `mr_cost_table` draw nothing themselves. Static edges to forestsearch from these functions come from string literals: `inherits(x, "forestsearch")` and `asNamespace("forestsearch")`. Building a gt object does not draw: `create_unique_id_vals` is called with `na_is_index = TRUE`. The draw happens only when gt prints the table to HTML.
- `knitr::kable`, `coxph`/`confint`, `sort`, `sessionInfo`, `future::availableCores` and the inline `sprintf` calls do not draw.

Effect of the planned save/restore:
- No visible number in the article changes. The fits' numbers come from streams the calls seed themselves: `seedit`, `dina_args$seed`, and the `subgroup.consistency` default inside fit-fs. The bootstrap numbers come from `seed = 8316951L`.
- What changes is the 10-letter HTML id and scoped CSS of the two gt tables. Today they are the same on every render, because the last bootstrap call leaves its `set.seed(8316951L)` in force. Afterwards they come from the document's own RNG stream, which no `set.seed` fixes, so they change on every render. Timing values already change on every render.
- This assumes the save/restore also wraps `forestsearch_bootstrap_dofuture`, which is exported and on the reach list. If it wrapped only `forestsearch()`, the bootstrap's own seed would keep the ids stable and nothing would change.

Reach-list caveat: `forestsearch` is on the list, but for fit-dina and fit-grf the real path returns before `subgroup.consistency` (forestsearch_main.R:2887 and :3072).

How the article is built: pkgdown only (.Rbuildignore:20; _pkgdown.yml:45 and :582).
- **Verifier (corrected) corrections:**
  - The investigator missed two later random draws, so rows should not be empty. The chunks boot-summary (qmd:287, autoprint of `summ_fs$table`) and mr-table (qmd:297, autoprint of `mr_estimates_table(...)`) each print a gt table. Both tables come from `gt::gt()` with no table id (format_bootstrap_table, bootstrap_summaries_helpers.R:159, reached from summarize_bootstrap_results.R:169/:411 and mr_estimates_table.R:153). No id is set anywhere in those files. When gt 1.4.0 prints a table to HTML (`knit_print.gt_tbl` -> `as.tags.gt_tbl`), it calls `random_id()`, which runs `sample(letters, 10, replace = TRUE)` on the main-session RNG. A check with Rscript --vanilla shows the call advances .Random.seed and gives the same id under the same seed. The investigator's static walk covered only forestsearch's own R/ code, so it could not see a draw inside a dependency's print method.
  - The investigator was wrong that only timing values (proc.time) would vary between renders after the change. Today the two gt container ids are stable: the last bootstrap call (GRF, `seed = 8316951L` at qmd:273) leaves its own `set.seed` (bootstrap_calculations_helpers.R:27, bootstrap_analysis_dofuture.R:334) in force when it returns. `forestsearch_bootstrap_dofuture` is exported (NAMESPACE:160) and on the reach list, so the planned save/restore would wrap it too. Every wrapped call would then return the RNG to its state at entry to fit-fs. The document has no `set.seed`, so that state is initialised from time and PID, and the two ids would change on every render. The visible table cells would not change.
  - Line fix: the seed argument in `.fs_apply_mr` is at fs_mr_inference_methods.R:172 (`seed = .g(mr_inference_args$seed, seedit)`), not :170.
  - Line fix: the `create_plots` branch of `summarize_bootstrap_results` is at summarize_bootstrap_results.R:302. Line 241 is inside the `create_timing_summary_table()` call.

### vignettes/articles/treatment_effect_definitions.qmd

- **Evaluation:** The YAML header (:1-20) sets no document-level execute or eval option. No _quarto.yml/_quarto.yaml or _metadata.yml governs the file: there is none in vignettes/articles/, vignettes/ or the repo root. The only _quarto.yml in the tree is quarto/_quarto.yml, a separate project in a sibling directory. When pkgdown 2.2.1 builds all articles it writes a temporary vignettes/_quarto.yaml that holds only project: render: '*.qmd', plus a format-only metadata file, so no execute options come from there either. Only five chunks have braces, and four of them run. (1) setup, :22-47, include=FALSE: knitr::opts_chunk$set (:23-30), a code_fold knit hook (:34-44) and library(gt) (:46). (2) dgm-vs-sim-table, :529-637: a static data.frame and a gt() table (:570). (3) comparison-table, :755-920: a static data.frame and a gt() table (:822). (4) pipeline-table, :925-1022: a static data.frame and a gt() table (:972). The fifth, worked-example (:1030-1100), is eval=FALSE (:1030) and does not run. The fences at :459-477, :506-517, :708-711 and :716-721 are plain ```r display blocks without braces, so they are shown and never run. There is no inline `r ` code. No results are loaded from saved files. No evaluated chunk calls any forestsearch function, and library(forestsearch) appears only at :1031 inside the eval=FALSE chunk.
- **Calls reaching the fixed seed:**
- **set.seed() calls and seed arguments:**
- **Rows:** 0
- **Notes:** There is nothing to flag in this file (rows = []). No evaluated chunk calls any of the 23 functions in reach_fixed_seed.txt that may reach the fixed seed 8316951, and evaluated code has no set.seed() and no seed argument.

Unevaluated code: the eval=FALSE worked-example chunk (:1030-1100) calls four functions. They are setup_gbsg_dgm (:1034, defined at R/setup_gbsg_dgm.R:88), compute_dgm_cde (:1063, R/oc_analyses.R:131), simulate_from_dgm(dgm, n = 500, seed = 1) (:1075, R/simulate_from_dgm.R:172) and cox_ahr_cde_analysis (:1093, R/cox_ahr_cde_wrapper.R:107). None of the four is on the 23-name list, so the chunk would not reach the fixed seed even if it were turned on.

Random draws in evaluated code are hidden ones only. Each of the three gt tables (:570, :822, :972) draws a random 10-letter HTML id when it is printed. In gt 1.4.0, for HTML output, knit_print.gt_tbl calls as.tags(); as.tags.gt_tbl calls random_id() when no table_id is set, and random_id is paste(sample(letters, n, replace = TRUE), collapse = ""). The tracked rendered HTML shows three such ids: sisjzbbyee, ejpdiicgjy and abpppmmaaz. These draws come from the document's own unseeded stream. No fixed-seed call comes before any of them, so the planned save and restore of .Random.seed changes nothing here.

The ids already change from one render to the next. pkgdown 2.2.1 builds .qmd articles through quarto::quarto_render in a separate process and sets no seed; its set.seed(1014) in rmarkdown_render_with_seed is used only for Rmd articles. The repo .Rprofile only sets options(pkgdown.internet = FALSE). The ids are not visible text, and no displayed value depends on them. Evaluated code contains no jitter, sample, rnorm or runif.

Build: pkgdown only. .Rbuildignore:20 has ^vignettes/articles$, which leaves the file out of the CRAN tarball. DESCRIPTION:30 has VignetteBuilder: knitr, quarto. The article is listed in _pkgdown.yml (navbar :51, articles index :591).

A side note on the task background, not on this result: the do.call(subgroup.consistency, sc_filtered_args) in R/forestsearch_main.R is at :3669.
- **Verifier (confirmed) corrections:**

## Part 3 — the legacy GBSG simulation wrapper

### L-A

- **Method:** Read-only, using different routes from the investigator. To enumerate definitions I grepped every 'function' keyword and every line-start assignment, and ran a brace-depth awk pass that prints only top-level statements; I also searched for assign(, <<-, setMethod and lambda shorthand. I read both files in full, apart from the body of .create_gbsg_dgm_ at :340-720; for that range I listed every call token instead. I listed man/*.Rd and grepped their name, alias and keyword lines. I compared blob hashes (git hash-object against git rev-parse v0.4.0:path and c40d5309:path) for the two R files, the nine Rd pages, _pkgdown.yml, NEWS.md, calibrate_helpers.R and fs_fdr_report.R. I confirmed v0.4.0 is a lightweight tag with cat-file -t, which returns commit 83b6d708. I traced history with git log -S on '.Deprecated(', on '### Superseded functions', on 'create_gbsg_dgm' in NEWS.md, on 'scan_k_inter_effect', on 'eff_at <- function' and on 'k_random_noise', and checked ancestry with git merge-base --is-ancestor. For callers I ran word-boundary and partial-name searches ('create_gbsg', 'simulate_from_gbsg', 'gbsg_dgm_') over R/, man/ and, at v0.4.0, git grep. I also searched for indirect calls through do.call, match.fun and get(. I read .calibrate_by_root in R/calibrate_helpers.R:44-90 to confirm how eff_at is invoked and that failures stop(). I checked NAMESPACE for every S3method and for any aft_dgm_flex method. I ran no R and edited no repository files; the only file written is a JSON validation copy in my scratch folder.
- **Per function:**
  - **cut_numeric (R/sim_aft_gbsg.R:51)**
    - facts: Internal: roxygen '@keywords internal' at sim_aft_gbsg.R:50, no @export, no NAMESPACE entry. In v0.4.0: present (v0.4.0:R/sim_aft_gbsg.R is byte-identical to c40d5309) and not exported in v0.4.0:NAMESPACE. Deprecation: none (no .Deprecated/.Defunct, no deprecated/legacy/superseded wording, no NEWS.md entry). Not listed in _pkgdown.yml. Help page: man/cut_numeric.Rd (keyword{internal} at :21; blob identical at v0.4.0). Calls no wrapper function (only stats::quantile at :52). Called by .create_gbsg_dgm_ at :281 (age), :283 (pgr) and :284 (nodes). No callers elsewhere in R/.
    - count: calls 0 wrapper fns; called by 1 wrapper fn (.create_gbsg_dgm_) at 3 sites; 0 callers outside the two files
    - example: dfa$z2 <- cut_numeric(dfa$age)   # R/sim_aft_gbsg.R:281
  - **cut_size (R/sim_aft_gbsg.R:69)**
    - facts: Internal: '@keywords internal' at :68, no @export, not in NAMESPACE. In v0.4.0: present (identical) and not exported. Deprecation: none, no NEWS entry. Not listed in _pkgdown.yml. Help page: man/cut_size.Rd (keyword{internal} at :20; identical at v0.4.0). Calls no wrapper function. Called by .create_gbsg_dgm_ at :298 only. No callers elsewhere in R/.
    - count: calls 0; called by 1 wrapper fn at 1 site; 0 outside
    - example: dfa$v6 <- as.factor(cut_size(dfa$size))   # R/sim_aft_gbsg.R:298
  - **create_gbsg_dgm (R/sim_aft_gbsg.R:199)**
    - facts: Exported: '@export' at :198, giving NAMESPACE:116 export(create_gbsg_dgm). In v0.4.0: present and exported (v0.4.0:NAMESPACE:102); also exported at v0.1.0 (v0.1.0:NAMESPACE:43). Deprecation: a runtime .Deprecated(new = 'generate_aft_dgm_flex', package = 'forestsearch', msg = 'create_gbsg_dgm is deprecated ... Then wrap with setup_gbsg_dgm() ...') at :211-225. It was added in commit 834d8384 (2026-03-07), before the v0.1.0 tag 29890c7c (2026-03-17), and is present at v0.1.0:R/sim_aft_gbsg.R:211. The lifecycle package is not used. The roxygen (:82-198) and man/create_gbsg_dgm.Rd have no deprecated or superseded wording; the Rd seealso (:143-146) still points to simulate_from_gbsg_dgm and calibrate_k_inter. NEWS.md has no current entry (0 hits at c40d5309 and at v0.4.0). The 0.2.0 'Superseded functions' entry (v0.1.0:NEWS.md:31-42) named setup_gbsg_dgm as its successor; it was deleted in fdecbe7b (2026-04-16), an ancestor of v0.4.0. The replacement is named inconsistently: .Deprecated says generate_aft_dgm_flex, while the setup_gbsg_dgm roxygen (setup_gbsg_dgm.R:26-27) calls setup_gbsg_dgm 'the replacement for create_gbsg_dgm()'. Listed at _pkgdown.yml:333 in the ordinary 'Data Generating Mechanisms' section (:326). Help page: man/create_gbsg_dgm.Rd (identical at v0.4.0); setup_gbsg_dgm inherits its params via @inheritParams create_gbsg_dgm (setup_gbsg_dgm.R:44). Calls .create_gbsg_dgm_ (:226). No wrapper code calls it. It appears in roxygen examples at :174, :177, :180, :781, :990, :995 and :1153, and in doc text and comments in setup_gbsg_dgm.R (:6, :17, :27, :31, :44, :60, :75). Outside the wrapper there are no code callers in R/, only doc uses: dontrun examples at R/oc_analyses.R:123 and :413 and an @param mention at :109.
    - count: calls 1 wrapper fn at 1 site; 0 code callers; 7 roxygen example uses in sim_aft_gbsg.R plus 2 outside (oc_analyses.R:123, :413)
    - example: .create_gbsg_dgm_(model = model, k_treat = k_treat, k_inter = k_inter, ...)   # R/sim_aft_gbsg.R:226
  - **.create_gbsg_dgm_ (R/sim_aft_gbsg.R:236)**
    - facts: Internal: no @export and no NAMESPACE entry. Its roxygen is only '@keywords internal' (:235) with no title, so there is no man page (no dot-create_gbsg_dgm_.Rd at c40d5309 or v0.4.0). The plain comment at :233-234 says it is used by calibrate_k_inter(), get_dgm_with_output(), scan_k_inter_effect() and setup_gbsg_dgm(), with no deprecation warning. That comment is stale: scan_k_inter_effect() is defined nowhere, now or in history; it appears only in this comment and its copy in the Rbuildignored dev/glm/scratch/R_legacy. In v0.4.0: present (identical) and not exported. Deprecation: none. Not listed in _pkgdown.yml. It builds the legacy object with class c('gbsg_dgm','list') (:735). Calls cut_numeric (:281, :283, :284) and cut_size (:298), and uses the constants DEFAULT_N_SUPER (:242), SEED_BASE (:245) and DAYS_PER_MONTH (:273). Called by create_gbsg_dgm (:226), eff_at inside calibrate_k_inter (:1024), get_dgm_with_output (:1108), validate_k_inter_effect (:1254) and setup_gbsg_dgm (setup_gbsg_dgm.R:106). Code outside the two files never calls it; it is named only in a comment at R/fs_fdr_report.R:485.
    - count: called by 5 wrapper fns at 5 sites (the hub of the wrapper); calls 2 wrapper fns at 4 sites, plus 3 constants
    - example: dgm <- .create_gbsg_dgm_(model = model, k_treat = k_treat, ...)   # R/setup_gbsg_dgm.R:106
  - **simulate_from_gbsg_dgm (R/sim_aft_gbsg.R:789)**
    - facts: Exported: '@export' at :788, giving NAMESPACE:293. In v0.4.0: present and exported (v0.4.0:NAMESPACE:266); exported at v0.1.0 (v0.1.0:NAMESPACE:125). Deprecation: a runtime .Deprecated(new = 'simulate_from_dgm', package = 'forestsearch', msg ...) at :800-812, added in the same commit 834d8384. The message tells users to pass analysis_time = Inf and explains the column renames. The roxygen (:745-788) and man/simulate_from_gbsg_dgm.Rd have no deprecated wording; @param dgm still says it is a gbsg_dgm object from create_gbsg_dgm (:750; Rd :20). NEWS.md has no current entry; the supersession text at v0.1.0:NEWS.md:39-42 was deleted in fdecbe7b. Listed at _pkgdown.yml:338. Help page: man/simulate_from_gbsg_dgm.Rd (identical at v0.4.0). Calls .simulate_from_gbsg_dgm_ (:813). No wrapper code calls it; it appears only in the roxygen example at :782. Outside the wrapper it appears only in comments and roxygen text (R/run_simulation_analysis.R:6, :198).
    - count: calls 1 wrapper fn at 1 site; called by 0
    - example: .simulate_from_gbsg_dgm_(dgm = dgm, n = n, rand_ratio = rand_ratio, ...)   # R/sim_aft_gbsg.R:813
  - **.simulate_from_gbsg_dgm_ (R/sim_aft_gbsg.R:823)**
    - facts: Internal: the comment at :821 reads 'Internal implementation. No deprecation warning.' Its roxygen is only '@keywords internal' (:822) with no title, so there is no man page. Not in NAMESPACE. In v0.4.0: present (identical) and not exported. Deprecation: none. Not listed in _pkgdown.yml. It requires inherits(dgm, 'gbsg_dgm') (:836). It reads dgm$df_super_rand (:845), dgm$cens_params (:840, :910-912), dgm$model_params$mu and $sigma (:897-898), and dgm$seed (:856, :884, :900, :908). Calls no wrapper function. Called by simulate_from_gbsg_dgm (:813) only. No callers outside the two files.
    - count: calls 0; called by 1 wrapper fn at 1 site
    - example: stopifnot('dgm must be a gbsg_dgm object' = inherits(dgm, 'gbsg_dgm'), ...)   # R/sim_aft_gbsg.R:835-838
  - **calibrate_k_inter (R/sim_aft_gbsg.R:1004)**
    - facts: Exported: '@export' at :1003, giving NAMESPACE:92. In v0.4.0: present and exported (v0.4.0:NAMESPACE:78), with the body identical. At v0.1.0 it used a uniroot 'objective' closure (v0.1.0:R/sim_aft_gbsg.R:1019); eff_at came in 33f05195, an ancestor of v0.4.0. Deprecation: none in code or roxygen, and no NEWS entry. Stale docs: @param ... says 'Additional arguments passed to create_gbsg_dgm' (:976; Rd :43), but the code forwards them to .create_gbsg_dgm_ (:1026). The dontrun examples call the deprecated create_gbsg_dgm (:990, :995; Rd :62, :67). Listed at _pkgdown.yml:342. Help page: man/calibrate_k_inter.Rd (identical at v0.4.0). Calls .create_gbsg_dgm_ (:1024) through the nested eff_at (:1023), which it passes as effect_fun at :1041. It also calls .calibrate_by_root (:1040), which is outside the wrapper (defined at R/calibrate_helpers.R:47, @noRd). Called by get_dgm_with_output (:1093). Outside the wrapper: R/fs_fdr_report.R:258 (the same at v0.4.0).
    - count: calls 1 wrapper fn (.create_gbsg_dgm_ via eff_at) and 1 external helper; called by 1 wrapper fn at 1 site; 1 caller outside (fs_fdr_report.R:258)
    - example: k_inter <- calibrate_k_inter(target_hr_harm = target_hr_harm, model = model_harm, ...)   # R/sim_aft_gbsg.R:1093
  - **eff_at (nested in calibrate_k_inter, R/sim_aft_gbsg.R:1023)**
    - facts: A local closure inside calibrate_k_inter: not exported, no NAMESPACE entry, no man page, no roxygen. In v0.4.0: present (file identical). Deprecation: not applicable. It calls .create_gbsg_dgm_(model, k_treat, k_inter = k_val, cens_type, verbose = FALSE, ...) at :1024-1027 and returns dgm$AHR_H_true or dgm$hr_H_true (:1028). The wrapper never calls it directly; it is passed as effect_fun to .calibrate_by_root (:1041). That helper calls it at R/calibrate_helpers.R:58 (objective <- function(k) effect_fun(k) - target, used by uniroot at :61) and again at :74 (achieved <- effect_fun(k_hat)), so each calibration does one extra build.
    - count: calls 1 wrapper fn at 1 site; invoked indirectly at 2 sites in 1 external helper
    - example: dgm <- .create_gbsg_dgm_(model = model, k_treat = k_treat, k_inter = k_val, cens_type = cens_type, verbose = FALSE, ...)   # R/sim_aft_gbsg.R:1024
  - **get_dgm_with_output (R/sim_aft_gbsg.R:1076)**
    - facts: Internal: '@keywords internal' at :1075, no @export, not in NAMESPACE. In v0.4.0: present (identical) and not exported. Deprecation: none, no NEWS entry. Inconsistency: it is listed in the public 'Data Generating Mechanisms' section at _pkgdown.yml:339, and the file header calls it a key function (:13). Help page: man/get_dgm_with_output.Rd (keyword{internal} at :59; @param ... says 'passed to create_gbsg_dgm' at Rd :45 and roxygen :1066). Calls calibrate_k_inter (:1093) and .create_gbsg_dgm_ (:1108). The NA fallback at :1101-1104 is dead code, because calibrate_k_inter goes through .calibrate_by_root, which calls stop() on failure (R/calibrate_helpers.R:60-70, :77-83) and never returns NA. Called by no wrapper function. No callers in R/.
    - count: calls 2 wrapper fns at 2 sites; called by 0
    - example: dgm <- .create_gbsg_dgm_(model = model_harm, k_treat = k_treat, k_inter = k_inter, ...)   # R/sim_aft_gbsg.R:1108
  - **print.gbsg_dgm (R/sim_aft_gbsg.R:1157)**
    - facts: S3 method: '@export' at :1156, giving NAMESPACE:48 S3method(print,gbsg_dgm). In v0.4.0: present and registered (v0.4.0:NAMESPACE:43); also at v0.1.0:NAMESPACE:15. Deprecation: none; its dontrun example uses the deprecated create_gbsg_dgm() (:1153; Rd :23). Listed at _pkgdown.yml:334. Help page: man/print.gbsg_dgm.Rd (identical at v0.4.0). Calls no wrapper function. It is never called explicitly; S3 dispatch reaches it for objects of class gbsg_dgm. Those come from .create_gbsg_dgm_ (:735) and from setup_gbsg_dgm (class c('aft_dgm_flex','gbsg_dgm','list') at setup_gbsg_dgm.R:182). R/ and NAMESPACE have no print.aft_dgm_flex or any other aft_dgm_flex method, so printing a setup_gbsg_dgm object also lands here, as setup_gbsg_dgm.R:11 and :162 expect.
    - count: calls 0; 0 explicit call sites; dispatch target for objects from 2 wrapper constructors
    - example: class(result) <- c('gbsg_dgm', 'list')   # R/sim_aft_gbsg.R:735 (makes print() dispatch here)
  - **validate_k_inter_effect (R/sim_aft_gbsg.R:1229)**
    - facts: Exported: '@export' at :1228, giving NAMESPACE:312. In v0.4.0: present and exported (v0.4.0:NAMESPACE:285); also at v0.1.0:NAMESPACE:134. Deprecation: none, no NEWS entry. Stale docs: @param ... says 'Additional arguments passed to create_gbsg_dgm' (:1216; Rd :19), but the code calls .create_gbsg_dgm_. Listed at _pkgdown.yml:344. Help page: man/validate_k_inter_effect.Rd (identical at v0.4.0). Calls .create_gbsg_dgm_ (:1254) inside 'for (k in k_inter_values)' (:1253); the default has 6 values (:1230), so 6 builds. Called by no wrapper function. No callers in R/.
    - count: calls 1 wrapper fn at 1 site (inside a loop); called by 0
    - example: dgm <- .create_gbsg_dgm_(model = 'alt', k_inter = k, verbose = FALSE, ...)   # R/sim_aft_gbsg.R:1254
  - **setup_gbsg_dgm (R/setup_gbsg_dgm.R:88)**
    - facts: Exported: '@export' at setup_gbsg_dgm.R:86, giving NAMESPACE:286. In v0.4.0: present and exported (v0.4.0:NAMESPACE:259; file identical). It has been exported since v0.1.0 (v0.1.0:NAMESPACE:120). The k_random_noise and noise_seed arguments came later, in commit 83bd19de (NEWS.md:1088, under the 0.2.2 development heading at :1084). Deprecation: none. It is presented as the successor: the header at :5 calls it 'the new public API', and the roxygen at :26-27 calls it 'the replacement for create_gbsg_dgm()'. The word 'legacy' appears at :12 (legacy top-level fields), :105 (the inner legacy build) and :168 (legacy $sigma). NEWS.md mentions it at :511 and :1088. Listed at _pkgdown.yml:332. Help page: man/setup_gbsg_dgm.Rd (identical at v0.4.0). Docs and code disagree: the header (:6, :17) and roxygen :31 (Rd :90) say it calls create_gbsg_dgm(), but the code calls .create_gbsg_dgm_ (:106), so no deprecation warning fires. It hard-codes n_super = 5000L and seed = 8316951L (:94, :97) instead of using DEFAULT_N_SUPER and SEED_BASE. Calls .create_gbsg_dgm_ (:106) and the nested .rename 5 times (:155-159). Called by no wrapper function; its name appears only as a string in the .Deprecated message (sim_aft_gbsg.R:223). Outside the wrapper: R/fs_fdr_report.R:256 and :262 (the same at v0.4.0), plus a string mention at R/fs_dgm_feasibility.R:130.
    - count: calls 1 wrapper fn and 1 nested helper (5 sites); called by 0 in the wrapper; 2 call sites outside (fs_fdr_report.R)
    - example: dgm <- .create_gbsg_dgm_(model = model, k_treat = k_treat, k_inter = k_inter, ...)   # R/setup_gbsg_dgm.R:106
  - **.rename (nested in setup_gbsg_dgm, R/setup_gbsg_dgm.R:149)**
    - facts: A local closure inside setup_gbsg_dgm: not exported, no man page, no roxygen. In v0.4.0: present (file identical); also at v0.1.0:R/setup_gbsg_dgm.R:109. Deprecation: not applicable. It renames dot-notation columns to the names simulate_from_dgm() expects. Calls no wrapper function. Called by setup_gbsg_dgm at :155, :156, :157, :158 and :159.
    - count: calls 0; called 5 times by 1 wrapper fn
    - example: df_s <- .rename(df_s, 'lin1.conf', 'lin_pred_1')   # R/setup_gbsg_dgm.R:155
  - **constants SEED_BASE (:29), DAYS_PER_MONTH (:32), DEFAULT_N_SUPER (:35) in R/sim_aft_gbsg.R**
    - facts: Top-level non-function definitions: not exported (no NAMESPACE entry), no man page. In v0.4.0: present (identical). Used in the create_gbsg_dgm defaults (n_super = DEFAULT_N_SUPER at :205, seed = SEED_BASE at :208), the .create_gbsg_dgm_ defaults (:242, :245) and its body (dfa$y <- dfa$rfstime / DAYS_PER_MONTH at :273). No uses elsewhere in R/; setup_gbsg_dgm repeats the literal values instead (setup_gbsg_dgm.R:94, :97).
    - count: 3 constants; used at 5 sites, all in sim_aft_gbsg.R
    - example: dfa$y <- dfa$rfstime / DAYS_PER_MONTH   # R/sim_aft_gbsg.R:273
- **Files named:**
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/R/sim_aft_gbsg.R (current; wrapper source; blob-identical to v0.4.0 and c40d5309)
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/R/setup_gbsg_dgm.R (current; wrapper source; blob-identical to v0.4.0 and c40d5309)
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/NAMESPACE (current)
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/NEWS.md (current; differs from v0.4.0)
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/_pkgdown.yml (current; differs from v0.4.0, but the DGM section has the same entries)
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/DESCRIPTION (current; Version 0.4.1.9007)
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/man/create_gbsg_dgm.Rd (current help page)
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/man/simulate_from_gbsg_dgm.Rd (current help page)
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/man/calibrate_k_inter.Rd (current help page)
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/man/get_dgm_with_output.Rd (current help page; keyword internal)
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/man/print.gbsg_dgm.Rd (current help page)
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/man/validate_k_inter_effect.Rd (current help page)
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/man/cut_numeric.Rd (current help page; keyword internal)
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/man/cut_size.Rd (current help page; keyword internal)
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/man/setup_gbsg_dgm.Rd (current help page)
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/man/compute_dgm_cde.Rd (current; another function's page whose example and @param use create_gbsg_dgm)
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/man/summarize_simulation_results.Rd (current; another function's page whose example uses create_gbsg_dgm)
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/R/calibrate_helpers.R (current; external helper .calibrate_by_root, which calls eff_at at :58 and :74)
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/R/fs_fdr_report.R (current; calls setup_gbsg_dgm and calibrate_k_inter from outside the wrapper)
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/R/fs_dgm_feasibility.R (current; mentions setup_gbsg_dgm only in a string)
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/R/oc_analyses.R (current; roxygen-only uses of create_gbsg_dgm at :109, :123 and :413)
  - /home/larryleon/Documents/GitHub/forestsearch v0.4.0:NAMESPACE (archive; CRAN 0.4.0 release, lightweight tag on 83b6d708)
  - /home/larryleon/Documents/GitHub/forestsearch v0.4.0:R/sim_aft_gbsg.R and v0.4.0:R/setup_gbsg_dgm.R (archive; identical to current)
  - /home/larryleon/Documents/GitHub/forestsearch v0.4.0:_pkgdown.yml, v0.4.0:NEWS.md, v0.4.0:R/fs_fdr_report.R (archive)
  - /home/larryleon/Documents/GitHub/forestsearch v0.4.0:man/*.Rd for the 9 wrapper pages (archive; all identical to current)
  - /home/larryleon/Documents/GitHub/forestsearch v0.1.0:NAMESPACE, v0.1.0:R/sim_aft_gbsg.R, v0.1.0:R/setup_gbsg_dgm.R, v0.1.0:NEWS.md (archive; first CRAN release; NEWS lines 31-42 hold the removed 'Superseded functions' entry)
  - /home/larryleon/Documents/GitHub/forestsearch dev/glm/scratch/R_legacy/R/{sim_aft_gbsg.R,setup_gbsg_dgm.R}, dev/code-working/mrct_simulations/{files5/legacy_code_beforeanything,files6/legacy,legacy_code}/R/sim_aft_gbsg.R (dev copies; tracked but excluded by the .Rbuildignore pattern ^dev$, so not package code)
- **Notes:** The given list is confirmed, with additions. All 11 functions sit at the stated lines (sim_aft_gbsg.R 51, 69, 199, 236, 789, 823, 1004, 1076, 1157, 1229; setup_gbsg_dgm.R 88). The only other definitions in the two files are:
- two nested closures: eff_at in calibrate_k_inter (sim_aft_gbsg.R:1023) and .rename in setup_gbsg_dgm (setup_gbsg_dgm.R:149);
- three top-level constants: SEED_BASE :29, DAYS_PER_MONTH :32 and DEFAULT_N_SUPER :35.

Both files and all 9 wrapper Rd pages are blob-identical between v0.4.0 (83b6d708, lightweight tag) and c40d5309, so every export and internal status is unchanged since the CRAN 0.4.0 release.

**Exported** (NAMESPACE line now / at v0.4.0):
- create_gbsg_dgm (116/102)
- simulate_from_gbsg_dgm (293/266)
- calibrate_k_inter (92/78)
- validate_k_inter_effect (312/285)
- setup_gbsg_dgm (286/259)
- print.gbsg_dgm as S3method (48/43)

All six were already exported at v0.1.0.

**Internal**:
- cut_numeric, cut_size and get_dgm_with_output: @keywords internal, each with a man page;
- .create_gbsg_dgm_ and .simulate_from_gbsg_dgm_: untitled '@keywords internal' blocks, so no man page;
- the two closures and the three constants.

**Deprecation**: only create_gbsg_dgm (:211-225, which names generate_aft_dgm_flex) and simulate_from_gbsg_dgm (:800-812, which names simulate_from_dgm) are deprecated. Both use a runtime .Deprecated() added in 834d8384 before v0.1.0. Neither roxygen nor Rd says so, and the lifecycle package is not used. The NEWS supersession entry, which named setup_gbsg_dgm as the successor, was removed in fdecbe7b before v0.4.0. _pkgdown.yml lists both in the ordinary DGM section (:326-345).

**Call graph**: .create_gbsg_dgm_ is the hub, with 5 callers. The deprecated public entry points are never called by wrapper code, so no internal path triggers the warnings. Outside the wrapper, only R/fs_fdr_report.R calls into it: setup_gbsg_dgm at :256 and :262, and calibrate_k_inter at :258. This is the same at v0.4.0.

**Stale documentation and code for a cleanup**:
- the nonexistent scan_k_inter_effect in the comment at :233-234;
- setup_gbsg_dgm docs (:6, :17, :31, Rd :90) say it calls create_gbsg_dgm;
- the '... passed to create_gbsg_dgm' @param in three functions (:976, :1066, :1216);
- dontrun examples calling the deprecated create_gbsg_dgm (:990, :995, :1153, plus R/oc_analyses.R:123 and :413);
- the internal get_dgm_with_output listed in _pkgdown.yml:339;
- the hard-coded 5000L and 8316951L in setup_gbsg_dgm (:94, :97);
- newly found: the unreachable is.na(k_inter) fallback in get_dgm_with_output (:1101-1104).
- **Verifier (corrected) corrections:**
  - The substance holds up. These all match: the definition list, the six NAMESPACE lines now and at v0.4.0, the help pages and keyword lines, the _pkgdown.yml lines, the v0.4.0 and v0.1.0 checks, the .Deprecated history (834d8384), the NEWS deletion (fdecbe7b) and every call site and count in the internal call graph. The items below are small corrections and additions only.
  - eff_at is called from two places in .calibrate_by_root, not just through objective. One is R/calibrate_helpers.R:58 (objective <- function(k) effect_fun(k) - target, used by uniroot at :61). The other is :74 (achieved <- effect_fun(k_hat)). So each calibrate_k_inter call does one more .create_gbsg_dgm_ build than the root search needs.
  - setup_gbsg_dgm: the word 'legacy' appears three times in R/setup_gbsg_dgm.R, not twice: :12 ('All legacy top-level fields ... are preserved'), :105 ('Run the legacy DGM') and :168 ('legacy uses $sigma'). Only :105 refers to the inner build.
  - create_gbsg_dgm is also used in documentation outside the two wrapper files, which the result did not mention. Two dontrun roxygen examples call it: R/oc_analyses.R:123 (rendered at man/compute_dgm_cde.Rd:46) and R/oc_analyses.R:413 (man/summarize_simulation_results.Rd:32). It is also named in @param dgm at R/oc_analyses.R:109 (man/compute_dgm_cde.Rd:10). None of these is code, so 'No code callers anywhere in R/' still holds.
  - Missed stale code in get_dgm_with_output: the fallback 'if (is.na(k_inter)) { warning(Calibration failed, using k_inter = 1) ... }' at R/sim_aft_gbsg.R:1101-1104 can never run. calibrate_k_inter returns .calibrate_by_root(...)$root (:1040-1044). .calibrate_by_root calls stop() on failure instead of returning NA, both when it cannot bracket a root (R/calibrate_helpers.R:60-70) and at the tolerance gate (:77-83).
  - .simulate_from_gbsg_dgm_ also reads dgm$seed, which the list of fields it reads left out: set.seed(dgm$seed + sim_id ...) at R/sim_aft_gbsg.R:856, :884, :900 and :908.
  - Line detail for man/create_gbsg_dgm.Rd: the seealso block opens at :143, its two links are at :144-145 and it closes at :146. The result's ':144-146' is close enough. The .Deprecated message also names setup_gbsg_dgm, inside a string at R/sim_aft_gbsg.R:223, which is not a call.
  - Same calls at the release: at v0.4.0 the only callers outside the wrapper are the same three call sites. Those are setup_gbsg_dgm at v0.4.0:R/fs_fdr_report.R:256 and :262, and calibrate_k_inter at :258; there is also the string mention at v0.4.0:R/fs_dgm_feasibility.R:130. R/fs_fdr_report.R as a file does differ between v0.4.0 and c40d5309.

### L-B

- **Method:** Read-only, using routes independent of the investigator's single combined regex. (a) git hash-object of every tracked file under R/, tests/, inst/ and man/, plus _pkgdown.yml, README.md, NAMESPACE and DESCRIPTION (706 files), compared against git ls-tree c40d5309, gave 0 differences from t2_tree. (b) Definitions were re-enumerated with three greps over the two wrapper files: '<- function', 'function(' anywhere (to catch anonymous or nested functions), and '\\(' lambdas, plus top-level assignments and the roxygen tags @rdname, @describeIn, @export, @keywords and @inheritParams. (c) For each name: git grep -n -w at c40d5309 over R, tests, inst, _pkgdown.yml, README.md, README.Rmd, NAMESPACE, DESCRIPTION and .Rbuildignore; git grep -F for the dotted names (.create_gbsg_dgm_, .simulate_from_gbsg_dgm_, print.gbsg_dgm, .rename); partial and substring greps (create_gbsg, simulate_from_gbsg, sim_aft_gbsg, setup_gbsg_dgm.R, k_inter_effect, dgm_with_output, gbsg_dgm, 'gbsg_dgm(', N_SUPER, SEED_, PER_MONTH, 8316951, aft_dgm_flex, cut_, scan_k_inter); case-insensitive 'gbsg' over tests/ and inst/. (d) Every hit was read in context and classified as code, #' roxygen, # comment or string, then mapped to its enclosing function through the top-level definition lines of each file. (e) NAMESPACE S3method and export lines, the man/ file list and \keyword{internal}, the absence of any print.aft_dgm_flex, the absence of any inherits(...,'gbsg_dgm') outside :836, and print(dgm) calls in R/ and tests were all checked. (f) Test gating: skip, NOT_CRAN and return gates in each relevant test file, tests/testthat.R, the helper files and .Rbuildignore. (g) Rscript --vanilla readRDS on the 2 tracked .rds test fixtures confirmed they hold no gbsg_dgm objects. (h) git grep -w on the names in the globals.R GBSG block showed which of them are used only in the wrapper files. No files were written.
- **Per function:**
  - **create_gbsg_dgm() [R/sim_aft_gbsg.R:199; NAMESPACE:116 export; .Deprecated spans :211-225, forwards to .create_gbsg_dgm_ at :226]**
    - facts: (1) R/ outside the two wrapper files: 0 code calls, 0 comments, 3 roxygen mentions, all in R/oc_analyses.R: :109 (@param dgm of compute_dgm_cde, :131), :123 (\dontrun example at :122 of compute_dgm_cde) and :413 (\dontrun example at :412 of summarize_simulation_results, :418). R/fs_fdr_report.R and R/betaHhat_truth.R: 0. Inside the wrapper files: no code call by name. Mentions there: header comment sim_aft_gbsg.R:10; own examples :174, :177, :180; deprecation-message string :215; roxygen :750, :781, :976, :990, :995, :1066, :1153, :1216. In setup_gbsg_dgm.R: comments :6 and :17, which are stale because the code calls .create_gbsg_dgm_ at :106; roxygen :27, :31, :60, :75 (@seealso); and :44 '@inheritParams create_gbsg_dgm', which is a documentation dependency: setup_gbsg_dgm's Rd arguments are inherited from this Rd. Generated mirrors: man/compute_dgm_cde.Rd:10, :46 and man/summarize_simulation_results.Rd:32. (2) tests/: 0. (3) inst/: 0 (30 tracked files). (4) _pkgdown.yml:333, under 'Data Generating Mechanisms' (:326). (5) README.md:131, in the key-functions table next to generate_aft_dgm_flex() and run_simulation_analysis(). README.Rmd is not tracked.
    - count: R code calls 0 / roxygen 3 (oc_analyses.R) / comments 0; tests 0; inst 0; pkgdown 1; README.md 1
    - example: R/oc_analyses.R:123 (roxygen \dontrun example of compute_dgm_cde): dgm <- create_gbsg_dgm(model = "alt", k_inter = 2.0)
  - **.create_gbsg_dgm_() [R/sim_aft_gbsg.R:236; bare '#' @keywords internal' at :235 with no title, so no man/ Rd; not in NAMESPACE]**
    - facts: (1) R/ outside the wrappers: 0 code calls and 1 comment, R/fs_fdr_report.R:485 (the comment spans :484-485) inside .fs_fdr_hr_row (:482). There are 5 code callers, all inside the wrapper files: create_gbsg_dgm sim_aft_gbsg.R:226; eff_at inside calibrate_k_inter :1024; get_dgm_with_output :1108; validate_k_inter_effect :1254; setup_gbsg_dgm setup_gbsg_dgm.R:106. The comment at sim_aft_gbsg.R:233-234 also names scan_k_inter_effect(), which is defined nowhere in the tracked repo. Its only other hit is the same comment copied in dev/glm/scratch/R_legacy/R/sim_aft_gbsg.R:234, so the comment is stale. R/globals.R:196-208 declares as globalVariables the column names this function creates (lin1.conf, linC1.conf, hlin.conf.1, h1.potential, Ts, ...); that is a coupling, not a call. (2) tests: 0 direct calls. It runs indirectly through setup_gbsg_dgm (test-knoise-population.R, test-betaHhat-contract.R, and test-plan-restored.R via fs_fdr_report) and through calibrate_k_inter (test-betaHhat-contract.R:419). (3) inst: 0. (4) pkgdown: 0. (5) README: 0.
    - count: R code calls outside wrappers 0 (5 inside) / comments 1 (fs_fdr_report.R:485); tests 0 direct; inst 0; pkgdown 0; README 0
    - example: R/setup_gbsg_dgm.R:106 (setup_gbsg_dgm): dgm <- .create_gbsg_dgm_(model = model, ...)
  - **simulate_from_gbsg_dgm() [R/sim_aft_gbsg.R:789; NAMESPACE:293 export; .Deprecated spans :800-812, forwards to .simulate_from_gbsg_dgm_ at :813]**
    - facts: (1) R/ outside the wrappers: 0 code calls; 1 comment, R/run_simulation_analysis.R:6 (file header); 1 roxygen mention, R/run_simulation_analysis.R:198 (description of run_simulation_analysis, :389; mirrored in man/run_simulation_analysis.Rd:204). fs_fdr_report.R and betaHhat_truth.R: 0. Inside the wrapper: header comment :11, @seealso :193 (create_gbsg_dgm docs), deprecation-message string :804, own example :782. (2) tests: 0. (3) inst: 0. (4) _pkgdown.yml:338. (5) README: 0.
    - count: R code calls 0 / roxygen 1 / comments 1; tests 0; inst 0; pkgdown 1; README 0
    - example: R/run_simulation_analysis.R:198 (roxygen only): '#' was coupled to \code{simulate_from_gbsg_dgm()} and GBSG-specific column'
  - **.simulate_from_gbsg_dgm_() [R/sim_aft_gbsg.R:823; bare '#' @keywords internal' at :822, so no Rd; not exported]**
    - facts: The only caller anywhere is simulate_from_gbsg_dgm() at sim_aft_gbsg.R:813. The stopifnot at :836 requires inherits(dgm, 'gbsg_dgm'), and that is the package's only inherits() check on 'gbsg_dgm'. (1) R/ outside the wrappers: 0 calls, roxygen or comment mentions (fixed-string grep). fs_fdr_report.R and betaHhat_truth.R: 0. (2)-(5): 0.
    - count: none outside sim_aft_gbsg.R (1 internal caller)
    - example: R/sim_aft_gbsg.R:813 (simulate_from_gbsg_dgm): .simulate_from_gbsg_dgm_(dgm = dgm, n = n, ...)
  - **calibrate_k_inter() [R/sim_aft_gbsg.R:1004; NAMESPACE:92 export]**
    - facts: (1) R/ outside the wrappers: 1 code call, in fs_fdr_report() at R/fs_fdr_report.R:258. It sits inside if ("harm" %in% null) at :257, and null defaults to c("hom","harm") (:183) through match.arg(several.ok = TRUE) at :215, so it runs by default. Roxygen 4: fs_fdr_report.R:172 (@seealso); calibrate_helpers.R:22 (block of .calibrate_by_root :47, tagged @noRd at :46, so no Rd output); calibrate_k_treat.R:9 and :90 (calibrate_k_treat, :96). Comments 3: calibrate_helpers.R:4 (file header); calibrate_glm_interaction.R:142 (inside calibrate_glm_interaction, :89); find_k_inter_main.R:129 (inside find_k_inter_for_target_hr, :88). R/betaHhat_truth.R: 0. Inside the wrappers: code call by get_dgm_with_output at :1093; header comment :12; @seealso :194; comment :233; own examples :987, :994. Outward dependency: it calls .calibrate_by_root at :1040-1041. (2) tests: tests/testthat/test-betaHhat-contract.R:419 runs it inside the helper .small_surv_dgm (:415), cached in .surv_dgm_cache (:414). Users of the helper: T11 :432 (skip_if_not(.have_shims()) :429); T12 :450 (test_that at :449, no skip, so it runs under R CMD check and on CRAN); T13 :466 (skip_if_not :462); T14 :503 (after skip_if_not :502). The file has no skip_on_cran and no NOT_CRAN gate. quarto/ is excluded by .Rbuildignore:24, so .have_shims() is FALSE under check. test-plan-restored.R passes null = "hom" (:169), so the calibrate branch is not reached there. (3) inst: 0. (4) _pkgdown.yml:342. (5) README: 0.
    - count: R code calls 1 (fs_fdr_report) + 1 inside wrappers / roxygen 4 / comments 3; tests 1 file runs it (test-betaHhat-contract.R); inst 0; pkgdown 1; README 0
    - example: R/fs_fdr_report.R:258 (fs_fdr_report): k_harm <- calibrate_k_inter(target_hr_harm = 1.0, model = "alt", use_ahr = FALSE)
  - **eff_at() [nested closure in calibrate_k_inter, R/sim_aft_gbsg.R:1023]**
    - facts: Used once, at sim_aft_gbsg.R:1041 (effect_fun = eff_at passed to .calibrate_by_root). Three other files define independent same-named local closures that do not reference this one: R/calibrate_glm_interaction.R:148 (used :176), R/calibrate_k_treat.R:133 (used :149), R/find_k_inter_main.R:105 (used :132). tests, inst, pkgdown, README: 0.
    - count: none outside calibrate_k_inter (1 internal use)
    - example: R/sim_aft_gbsg.R:1041: effect_fun = eff_at
  - **get_dgm_with_output() [R/sim_aft_gbsg.R:1076; @keywords internal :1075; not in NAMESPACE; man/get_dgm_with_output.Rd:59 \keyword{internal}]**
    - facts: (1) R/: 0 calls, roxygen or comment mentions outside the defining file. Inside it: header comment :13, comment :233. fs_fdr_report.R and betaHhat_truth.R: 0. (2) tests: 0. (3) inst: 0. (4) _pkgdown.yml:339 lists this internal topic in the reference index. (5) README: 0. Nothing in the package calls it.
    - count: none (pkgdown entry only)
    - example: none (only listing: _pkgdown.yml:339 '- get_dgm_with_output')
  - **print.gbsg_dgm() [R/sim_aft_gbsg.R:1157; NAMESPACE:48 S3method(print,gbsg_dgm); man/print.gbsg_dgm.Rd]**
    - facts: Nothing calls it by name; it is reached only by S3 dispatch, on a .create_gbsg_dgm_ result (class set at sim_aft_gbsg.R:735) or a setup_gbsg_dgm result (class c("aft_dgm_flex","gbsg_dgm","list") at setup_gbsg_dgm.R:182). No print.aft_dgm_flex exists in R/ or NAMESPACE, so setup objects print through this method. fs_fdr_report does not keep the DGMs in its return value (:311-326), and print.fs_fdr_report (:505) prints only data frames. (1) R/ outside the wrappers: no code prints a gbsg_dgm and nothing mentions the method by name. The class string 'gbsg_dgm' appears in roxygen 3 times (oc_analyses.R:41 get_dgm_hr @param; oc_analyses.R:582 format_oc_results @param; simulation_tables.R:157 invert_dgm_hrs @param) and in comments 5 times (oc_analyses.R:133; simulation_tables.R:311, :629, :651, :1151). Inside the wrappers: \dontrun print(dgm) examples at setup_gbsg_dgm.R:82 and sim_aft_gbsg.R:1154; roxygen setup_gbsg_dgm.R:42; comments setup_gbsg_dgm.R:11, :162. (2) tests: 0 prints. test-fs-dgm-feasibility.R:26 builds structure(list(), class = "gbsg_dgm") only for expect_error. The tracked .rds fixtures contain no gbsg_dgm objects (checked with readRDS). (3) inst: 0. (4) _pkgdown.yml:334. (5) README: 0.
    - count: R code calls 0 (dispatch only) / by-name mentions outside wrappers 0 (class-string: roxygen 3, comments 5); tests 0 (1 class-string use); inst 0; pkgdown 1; README 0
    - example: none outside wrappers; dispatch example R/setup_gbsg_dgm.R:82 (\dontrun roxygen): print(dgm)
  - **validate_k_inter_effect() [R/sim_aft_gbsg.R:1229; NAMESPACE:312 export]**
    - facts: Calls .create_gbsg_dgm_ at :1254. (1) R/ outside the wrappers: 0 calls, roxygen or comment mentions; its only appearance is its own \dontrun example at :1223. fs_fdr_report.R and betaHhat_truth.R: 0. (2) tests: 0. (3) inst: 0. (4) _pkgdown.yml:344. (5) README: 0.
    - count: none (pkgdown entry + NAMESPACE export only)
    - example: none (only listing: _pkgdown.yml:344)
  - **cut_numeric() [R/sim_aft_gbsg.R:51; @keywords internal :50; not exported; man/cut_numeric.Rd:21 \keyword{internal}]**
    - facts: The only caller is .create_gbsg_dgm_, with 3 calls: :281 (age), :283 (pgr), :284 (nodes). (1) R/ outside: 0. (2) tests: 0 direct; it runs indirectly whenever setup_gbsg_dgm or calibrate_k_inter runs. (3) inst: 0. (4) pkgdown: not listed. (5) README: 0.
    - count: none outside sim_aft_gbsg.R (3 internal calls)
    - example: R/sim_aft_gbsg.R:281 (.create_gbsg_dgm_): dfa$z2 <- cut_numeric(dfa$age)
  - **cut_size() [R/sim_aft_gbsg.R:69; @keywords internal :68; not exported; man/cut_size.Rd:20 \keyword{internal}]**
    - facts: The only caller is .create_gbsg_dgm_ at :298. (1) R/ outside: 0. (2) tests: 0 direct (indirect only). (3) inst: 0. (4) pkgdown: not listed. (5) README: 0.
    - count: none outside sim_aft_gbsg.R (1 internal call)
    - example: R/sim_aft_gbsg.R:298 (.create_gbsg_dgm_): dfa$v6 <- as.factor(cut_size(dfa$size))
  - **setup_gbsg_dgm() [R/setup_gbsg_dgm.R:88; NAMESPACE:286 export]**
    - facts: (1) R/ outside the wrappers: 2 code calls, both in fs_fdr_report(): R/fs_fdr_report.R:256 (hom null, k_inter = 0) and :262 (harm null, k_inter = k_harm). fs_dgm_feasibility() names it in a stop() string at R/fs_dgm_feasibility.R:130 (the message spans :128-131); that is not a call. Roxygen 7: fs_dgm_feasibility.R:84 (@note); fs_fdr_report.R:172 (@seealso); run_simulation_analysis.R:205 (@param), :339 (@seealso), :373 (\dontrun example), all for run_simulation_analysis :389; simulate_from_dgm.R:167 (\dontrun example of simulate_from_dgm :172); :436 (\dontrun example of check_censoring_dgm :440). Comments 3: fs_fdr_report.R:32 (above .FS_GBSG_CONFOUNDERS :34), :484 (in .fs_fdr_hr_row), oc_analyses.R:12 (a mention of the file name 'setup_gbsg_dgm.R'). Inside the wrappers: deprecation string sim_aft_gbsg.R:223, comment :234. R/betaHhat_truth.R: 0 references. It reads df_super and model_params$censoring (:219, :230-233, :268, :828, :849, :854, :858), but these are the general aft_dgm_flex contract, which generate_aft_dgm_flex also produces, not a setup-specific dependency. (2) tests, 3 files run it: (a) tests/testthat/test-knoise-population.R:8, helper .make_gbsg_noise_dgm (n_super = 1000L), called at :33, :44, :59, :60, :61, :72, :73, :82 across 5 test_that blocks (:32, :43, :58, :71, :81), no skips. The name also appears in comment :2 and in the test titles :32, :43, :58, :71. (b) tests/testthat/test-betaHhat-contract.R:421, helper .small_surv_dgm (n_super = 3000L, seed = 8316951L at :422), gated like calibrate_k_inter; T12 (:449-459) is unconditional. (c) tests/testthat/test-plan-restored.R:168, transitively via fs_fdr_report(null = "hom", n_super = 2000L), with skip_on_cran() at :165 and %dofuture% mocked at :166. tests/testthat.R has no global skip. test-fs-dgm-feasibility.R:26-28 asserts 'Survival DGMs' from the stop() message that names setup_gbsg_dgm. It does not run setup_gbsg_dgm. (3) inst: 0. (4) _pkgdown.yml:332. (5) README.md: 0.
    - count: R code calls 2 (fs_fdr_report) / roxygen 7 / comments 3 / error-string 1 (fs_dgm_feasibility.R:130); tests 3 files run it (2 direct, 1 transitive; skip_on_cran only in test-plan-restored.R); inst 0; pkgdown 1; README 0
    - example: R/fs_fdr_report.R:256 (fs_fdr_report): setup_gbsg_dgm(model = "alt", k_inter = 0, n_super = n_super)
  - **.rename() [nested in setup_gbsg_dgm, R/setup_gbsg_dgm.R:149]**
    - facts: Used only at setup_gbsg_dgm.R:155-159 (5 calls: lin1.conf, lin0.conf, linC1.conf, linC0.conf, flag.harm). No other tracked R/, tests or inst file defines or calls .rename (fixed-string grep). pkgdown, README: 0.
    - count: none outside setup_gbsg_dgm (5 internal calls)
    - example: R/setup_gbsg_dgm.R:155: df_s <- .rename(df_s, "lin1.conf", "lin_pred_1")
  - **SEED_BASE (constant, R/sim_aft_gbsg.R:29 = 8316951L)**
    - facts: Used only as a default argument: sim_aft_gbsg.R:208 (create_gbsg_dgm) and :245 (.create_gbsg_dgm_). No reference to the symbol outside sim_aft_gbsg.R. setup_gbsg_dgm hard-codes the value (setup_gbsg_dgm.R:97) rather than using the constant. The literal 8316951 is the package-wide default seed and is hard-coded in about 20 R/ signatures (e.g. fs_fdr_report.R:193, bootstrap_calculations_helpers.R:103, generate_aft_dgm_main.R:238) and in tests (e.g. test-betaHhat-contract.R:422).
    - count: none outside sim_aft_gbsg.R
    - example: R/sim_aft_gbsg.R:208: seed = SEED_BASE
  - **DAYS_PER_MONTH (constant, R/sim_aft_gbsg.R:32 = 30.4375)**
    - facts: Used once, in .create_gbsg_dgm_ at :273. No other reference in R/, tests, inst, pkgdown or README.
    - count: none outside sim_aft_gbsg.R
    - example: R/sim_aft_gbsg.R:273: dfa$y <- dfa$rfstime / DAYS_PER_MONTH
  - **DEFAULT_N_SUPER (constant, R/sim_aft_gbsg.R:35 = 5000L)**
    - facts: Used only as a default argument at :205 and :242. setup_gbsg_dgm hard-codes n_super = 5000L (setup_gbsg_dgm.R:94). calibrate_k_inter's eff_at passes no n_super, so every calibration build uses this 5000 default, including the T12 test path. No other reference.
    - count: none outside sim_aft_gbsg.R
    - example: R/sim_aft_gbsg.R:205: n_super = DEFAULT_N_SUPER
  - **[coupling, not a definition] R/globals.R GBSG globalVariables block (:196-208, inside utils::globalVariables at :13)**
    - facts: The header comment at :197 reads 'GBSG DGM internal variables (from sim_aft_gbsg_refactored.R)'. That file does not exist, so the header is stale. The block declares z1-z5, zh, v1-v7, grade3, lin.conf.true, lin1.conf, lin0.conf, linC1.conf, linC0.conf, hlin.conf.1, hlin.conf.0, hlin.ratio, h1.potential, h0.potential, Ts, es, t.sim. Of these, 16 occur in R/ only in globals.R and the two wrapper files: z3, z4, z5, zh, lin.conf.true, lin1.conf, lin0.conf, linC1.conf, linC0.conf, hlin.conf.1, hlin.conf.0, hlin.ratio, h1.potential, h0.potential, Ts, t.sim. The others are shared with other R/ files. The block serves R CMD check for the wrapper's column names, and removing the wrapper would leave those 16 orphaned (harmless).
    - count: 1 block, 16 wrapper-only names
    - example: R/globals.R:203: "lin.conf.true", "lin1.conf", "lin0.conf",
- **Files named:**
  - /home/larryleon/Documents/GitHub/forestsearch/R/sim_aft_gbsg.R (current package code; wrapper definitions)
  - /home/larryleon/Documents/GitHub/forestsearch/R/setup_gbsg_dgm.R (current package code; wrapper definitions)
  - /home/larryleon/Documents/GitHub/forestsearch/R/fs_fdr_report.R (current package code; only code caller outside the wrappers: setup_gbsg_dgm x2 at :256/:262, calibrate_k_inter x1 at :258)
  - /home/larryleon/Documents/GitHub/forestsearch/R/betaHhat_truth.R (current package code; checked: 0 wrapper references; reads generic aft_dgm_flex fields only)
  - /home/larryleon/Documents/GitHub/forestsearch/R/globals.R (current package code; MISSED by investigator: GBSG globalVariables block :196-208)
  - /home/larryleon/Documents/GitHub/forestsearch/R/oc_analyses.R (current package code; roxygen :41, :109, :123, :413, :582; comments :12, :133)
  - /home/larryleon/Documents/GitHub/forestsearch/R/run_simulation_analysis.R (current package code; comment :6; roxygen :198, :205, :339, :373)
  - /home/larryleon/Documents/GitHub/forestsearch/R/simulate_from_dgm.R (current package code; roxygen \dontrun examples :167, :436)
  - /home/larryleon/Documents/GitHub/forestsearch/R/fs_dgm_feasibility.R (current package code; @note :84 + stop() string :130)
  - /home/larryleon/Documents/GitHub/forestsearch/R/calibrate_helpers.R (current package code; comment :4, @noRd roxygen :22; .calibrate_by_root callee :47)
  - /home/larryleon/Documents/GitHub/forestsearch/R/calibrate_k_treat.R (current package code; roxygen :9, :90; independent eff_at :133)
  - /home/larryleon/Documents/GitHub/forestsearch/R/calibrate_glm_interaction.R (current package code; comment :142; independent eff_at :148)
  - /home/larryleon/Documents/GitHub/forestsearch/R/find_k_inter_main.R (current package code; comment :129; independent eff_at :105)
  - /home/larryleon/Documents/GitHub/forestsearch/R/simulation_tables.R (current package code; gbsg_dgm class-string: roxygen :157, comments :311, :629, :651, :1151)
  - /home/larryleon/Documents/GitHub/forestsearch/R/generate_aft_dgm_main.R (current package code; consulted only to show df_super/censoring fields are the general aft_dgm_flex contract)
  - /home/larryleon/Documents/GitHub/forestsearch/tests/testthat/test-knoise-population.R (test; runs setup_gbsg_dgm, 5 blocks, no skips)
  - /home/larryleon/Documents/GitHub/forestsearch/tests/testthat/test-betaHhat-contract.R (test; runs calibrate_k_inter + setup_gbsg_dgm via .small_surv_dgm; T12 unconditional)
  - /home/larryleon/Documents/GitHub/forestsearch/tests/testthat/test-plan-restored.R (test; transitive setup_gbsg_dgm via fs_fdr_report(null='hom'), skip_on_cran :165)
  - /home/larryleon/Documents/GitHub/forestsearch/tests/testthat/test-fs-dgm-feasibility.R (test; class 'gbsg_dgm' stub at :26, asserts 'Survival DGMs' message; no wrapper call)
  - /home/larryleon/Documents/GitHub/forestsearch/tests/testthat.R (test runner; no global skip)
  - /home/larryleon/Documents/GitHub/forestsearch/tests/testthat/fixtures/grf_disjunction.rds (test fixture; readRDS checked, no gbsg_dgm)
  - /home/larryleon/Documents/GitHub/forestsearch/tests/testthat/fixtures/mrct_region_sims_record/point_6b5c85c0.rds (test fixture; readRDS checked, no gbsg_dgm)
  - /home/larryleon/Documents/GitHub/forestsearch/_pkgdown.yml (pkgdown config; :332-334, :338-339, :342, :344 under 'Data Generating Mechanisms' :326)
  - /home/larryleon/Documents/GitHub/forestsearch/README.md (README; :131 create_gbsg_dgm only; README.Rmd untracked)
  - /home/larryleon/Documents/GitHub/forestsearch/NAMESPACE (S3method :48; exports :92, :116, :286, :293, :312)
  - /home/larryleon/Documents/GitHub/forestsearch/.Rbuildignore (quarto excluded at :24)
  - /home/larryleon/Documents/GitHub/forestsearch/man/get_dgm_with_output.Rd (man; keyword internal :59)
  - /home/larryleon/Documents/GitHub/forestsearch/man/cut_numeric.Rd (man; keyword internal :21)
  - /home/larryleon/Documents/GitHub/forestsearch/man/cut_size.Rd (man; keyword internal :20)
  - /home/larryleon/Documents/GitHub/forestsearch/inst/ (30 tracked files: 29 in inst/extdata/mrct-consistency, 1 in inst/patches; 0 hits even case-insensitive 'gbsg')
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree (scratch export; verified identical to c40d5309 for all 706 in-scope tracked files)
- **Notes:** The definition list is confirmed. There are 11 top-level functions at the stated lines, plus nested eff_at (sim_aft_gbsg.R:1023), nested .rename (setup_gbsg_dgm.R:149) and the constants SEED_BASE (:29), DAYS_PER_MONTH (:32) and DEFAULT_N_SUPER (:35). There are no anonymous or lambda functions and no @rdname or @describeIn. .create_gbsg_dgm_ and .simulate_from_gbsg_dgm_ carry only a bare '#' @keywords internal' with no title, so neither has an Rd.

Confirmed main findings: (a) fs_fdr_report() is the only package function outside the two files that calls the wrapper: setup_gbsg_dgm at :256 and :262, calibrate_k_inter at :258. It has no callers in R/; its only caller anywhere is test-plan-restored.R:168. (b) R/betaHhat_truth.R references no wrapper name. (c) .create_gbsg_dgm_, .simulate_from_gbsg_dgm_, cut_numeric, cut_size, eff_at, .rename and the three constants are used only inside the wrapper files. (d) setup_gbsg_dgm.R:44 '@inheritParams create_gbsg_dgm' is a documentation dependency. (e) Stale comments: setup_gbsg_dgm.R:6 and :17 (say create_gbsg_dgm; the code calls .create_gbsg_dgm_ at :106) and sim_aft_gbsg.R:234 (scan_k_inter_effect is defined nowhere in the tracked repo). (f) README.md:131 lists the deprecated create_gbsg_dgm(). (g) _pkgdown.yml:339 lists the internal get_dgm_with_output. (h) On CRAN, setup_gbsg_dgm, calibrate_k_inter and, transitively, .create_gbsg_dgm_, eff_at, cut_numeric, cut_size and .rename run through test-betaHhat-contract.R T12. setup_gbsg_dgm also runs on CRAN through all 5 survival tests in test-knoise-population.R. No test calls create_gbsg_dgm, simulate_from_gbsg_dgm, get_dgm_with_output, validate_k_inter_effect or print.gbsg_dgm. (i) inst/ has 0 hits.

Corrections and additions: (1) R/globals.R:196-208 is a GBSG globalVariables block whose header comment is stale; 16 of its names are used only by the wrapper. The investigator missed it. (2) For print.gbsg_dgm, two class-string mentions classed as comments are roxygen, and oc_analyses.R:582 was missed. (3) betaHhat_truth.R's field reads are generic aft_dgm_flex and the list of them was incomplete. (4) oc_analyses.R:12 mentions a file name and its claim is inaccurate. (5) Minor mentions inside the wrapper files were missed. None of the corrections change a code-call count.
- **Verifier (corrected) corrections:**
  - MISSED COUPLING (R/ outside the wrapper files): R/globals.R:196-208, inside the utils::globalVariables(c( call that opens at R/globals.R:13. The block is headed '# GBSG DGM internal variables (from sim_aft_gbsg_refactored.R)' (:197) and declares the wrapper's NSE column names. 16 of them occur in R/ only in globals.R and the two wrapper files: z3, z4, z5, zh, lin.conf.true, lin1.conf, lin0.conf, linC1.conf, linC0.conf, hlin.conf.1, hlin.conf.0, hlin.ratio, h1.potential, h0.potential, Ts, t.sim (checked with git grep -w). The rest (z1, z2, v1-v7, grade3, es) are also used by other R/ files. The header names a file that does not exist (sim_aft_gbsg_refactored.R), so it is a stale comment. This is not a call, but removing the wrapper would leave those 16 entries orphaned. The investigator's claim that every outside reference was found therefore missed this one.
  - print.gbsg_dgm class-string classification: oc_analyses.R:41 (#' @param of get_dgm_hr, :49) and simulation_tables.R:157 (#' @param of invert_dgm_hrs, :169) are roxygen, not comments. The investigator also missed R/oc_analyses.R:582 (#' @param dgm '(\code{gbsg_dgm}, ...)' of format_oc_results, :618). Correct tally of class-string mentions outside the wrappers: roxygen 3 (oc_analyses.R:41, :582; simulation_tables.R:157) and comments 5 (oc_analyses.R:133 in compute_dgm_cde; simulation_tables.R:311 in build_classification_table, :629 and :651 in build_estimation_table, :1151 in interpret_estimation_table). Inside the wrappers, setup_gbsg_dgm.R:42 is roxygen ('#' ... \code{print.gbsg_dgm} continue to work'), not a comment. The only comments that name print.gbsg_dgm are setup_gbsg_dgm.R:11 and :162.
  - betaHhat_truth.R field reads were listed incompletely and are not specific to setup_gbsg_dgm. Besides :219 (dgm$model_params$censoring, .fs_random_cens_prob :218) and :230/:828, it reads dgm$df_super at :268 (.fs_exact_setup :246) and :849, :854, :858 (fs_build_eval_frame :820), and it requires df_super$lin_pred_cens_0/_1 at :231-233. These are the general aft_dgm_flex contract: generate_aft_dgm_flex also documents them, e.g. generate_aft_dgm_main.R:58 and generate_aft_dgm_helpers.R:80. The finding of 0 references to any wrapper name in R/betaHhat_truth.R is confirmed.
  - R/oc_analyses.R:12 ('#   - setup_gbsg_dgm.R (compute_dgm_cde)') mentions a file name, not the function. Its claim is also inaccurate: setup_gbsg_dgm.R never calls compute_dgm_cde in code, and names it only in comments and roxygen (:11, :42, :76, plus the \dontrun example :81, and :162).
  - Minor mentions inside the wrapper files were missed (they do not change any outside count). create_gbsg_dgm: header comment sim_aft_gbsg.R:10 and the deprecation-message string :215. simulate_from_gbsg_dgm: @seealso :193 and the deprecation-message string :804. calibrate_k_inter: header comment :12, @seealso :194, comment :233, and its own \dontrun examples :987 and :994. In setup_gbsg_dgm.R only :75 is @seealso. :27, :31 and :60 are description and @return text.
  - calibrate_helpers.R:22 (roxygen of .calibrate_by_root) is in a block tagged @keywords internal and @noRd (:45-46), so it produces no Rd. The classification as roxygen still stands.
  - SEED_BASE note: the hard-coded 8316951L is the package-wide default seed, not something particular to fs_fdr_report.R:193. It appears as a literal default in about 20 R/ function signatures (e.g. bootstrap_calculations_helpers.R:103, generate_aft_dgm_main.R:238) and in tests/testthat/test-betaHhat-contract.R:422 (seed passed to setup_gbsg_dgm). No file references the SEED_BASE symbol outside sim_aft_gbsg.R, which confirms the investigator's count.
  - Additional test detail for setup_gbsg_dgm: tests/testthat/test-fs-dgm-feasibility.R:26-28 expects the error 'Survival DGMs', which comes from the stop() message at R/fs_dgm_feasibility.R:128-131. That message names setup_gbsg_dgm() at :130. The test does not run setup_gbsg_dgm, and its regex does not include the name. In test-knoise-population.R the name also appears in comment :2 and in test_that description strings :32, :43, :58, :71 (not calls).

### L-C

- **Method:** Read-only. Re-derived from the c40d5309 export at /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree, after checking that every tracked vignettes/ file, both wrapper files, NAMESPACE and _pkgdown.yml are byte-identical to c40d5309 (git rev-parse of each blob vs git hash-object; 23 tracked vignettes files). Definitions: grepped `function(` plus all top-level `<-` assignments, and checked for lambdas, local() and assign(); this found the three constants. Mentions: git grep -n with partial and loose patterns (gbsg_dgm, calibrate_k, validate_k_inter, cut_numeric, cut_size, eff_at, get_dgm_with, hyphenated chunk-label variants, constant names, near names such as calibrate_k_treat, find_k_inter, get_dgm_hr), per-file line counts with git grep -c, and occurrence counts with grep -o (equal, one per line). Read every hit in context and mapped it to chunk fences and options, the YAML execute: blocks and knitr opts_chunk. Checked for inline `r` code, child includes, source() and readRDS. Indirect reach: every R/ caller of each wrapper name, string-based calls (do.call, match.fun, quoted names) and S3 methods for aft_dgm_flex and gbsg_dgm. Checked the deployed pkgdown render (read-only git show origin/gh-pages) and confirmed with git diff that its source commit differs from c40d5309 in paper_simulations.qmd only at :1340-1354. Did not run R and wrote no files.
- **Per function:**
  - **setup_gbsg_dgm() (R/setup_gbsg_dgm.R:88)**
    - facts: Exported (NAMESPACE:286). Named only in two pkgdown-only articles plus one tracked rendered HTML; vignettes/forestsearch.qmd (the only CRAN-built vignette) has 0. paper_simulations.qmd: 21 lines/21 occurrences. Evaluated calls (chunks with no eval option; document has no execute: block and opts_chunk at :24-30 sets only echo/message/warning/fig): :507 (create-dgm-alt :500), :614 (calibrate-k-inter-cox :600), :645 (calibrate-k-inter-ahr :631), :699 (create-dgm-null :695). Not run: :375 (setup-gbsg-dgm-usage :371, `#| eval: false` :372), :1756 (complete-example :1728, eval=FALSE). Plain ```r display block :271-286: call at :274. Comment in eval:false chunk simulate-from-either (:473-474): :476. Text (13): :160, :312, :362, :364, :388, :392, :442, :453, :454, :459, :469, :1821 (summary table), and :359, which sits in a mid-document YAML block (:358-360) that replaces the article's title: the deployed pkgdown page origin/gh-pages articles/paper_simulations.html has <title> (:8) and <h1> (:70) = 'setup_gbsg_dgm() and generate_aft_dgm_flex(): Relationship and Usage'. treatment_effect_definitions.qmd: 6 lines. Not run: :1034 (worked-example :1030, eval=FALSE). Comment in plain ```r block :459-477: :460. String literals in evaluated echo=FALSE, code_fold=TRUE gt-table chunks: :552 (dgm-vs-sim-table :529) and :785 (comparison-table :755); these render as table cells and as folded code, not calls. Prose: :163, :452. Document has no eval setting (opts_chunk :23-30). treatment_effect_definitions.html (tracked rendered): 8 lines, none evaluated: :2515 (prose), :2803 (prose), :2804 (plain code block comment), :3306 (gt cell), :3346 (code_fold copy), :4022 (gt cell), :4076 (code_fold copy), :4810 (eval=FALSE worked-example code shown).
    - count: paper_simulations.qmd 21 (4 evaluated calls, 2 eval=FALSE calls, 1 call in a plain display block, 1 comment, 13 text, :359 being the rendered page title); treatment_effect_definitions.qmd 6 (0 evaluated, 1 eval=FALSE call, 1 comment in a plain block, 2 string literals in table chunks, 2 prose); treatment_effect_definitions.html 8 rendered lines; forestsearch.qmd 0; all other files 0. Call expressions across files: 8 (4 evaluated, 3 eval=FALSE, 1 plain display block).
    - example: vignettes/articles/paper_simulations.qmd:507 `dgm_alt <- setup_gbsg_dgm(` in evaluated chunk create-dgm-alt (:500). eval=FALSE example: treatment_effect_definitions.qmd:1034 `dgm <- setup_gbsg_dgm(`.
  - **calibrate_k_inter() (R/sim_aft_gbsg.R:1004)**
    - facts: Exported (NAMESPACE:92). paper_simulations.qmd: 8 lines/8 occurrences. Evaluated: :604 (calibrate-k-inter-cox :600), :635 (calibrate-k-inter-ahr :631). eval=FALSE: :1749 (complete-example :1728). Commented-out call in the same chunk: :1753. Plain ```r block :271-286: :273. Text: :171, :596, :1822 (summary table). extreme_subgroups.qmd: 2 prose mentions, :202 and :485. Both are outside every block: fenced blocks are :218-220 and :224-229 (indented plain ```r), :259-261 (plain), and chunks :270-278, :280-301, :312-336, :448-456. The YAML execute block (:7-10) sets only echo, warning and message. The same file names the non-wrapper sibling calibrate_k_treat() on 8 lines (R/calibrate_k_treat.R:96), which is not counted.
    - count: paper_simulations.qmd 8 (2 evaluated, 1 eval=FALSE call, 1 commented-out, 1 in a plain display block, 3 text); extreme_subgroups.qmd 2 (prose only); all other files 0. Evaluated calls 2.
    - example: vignettes/articles/paper_simulations.qmd:604 `k_inter_cox <- calibrate_k_inter(` in evaluated chunk calibrate-k-inter-cox (:600).
  - **validate_k_inter_effect() (R/sim_aft_gbsg.R:1229)**
    - facts: Exported (NAMESPACE:312). Appears only in paper_simulations.qmd. Evaluated call :683 (validate-k-inter :678). Text: :675 (prose), :1823 (summary table). Its verbose output 'Testing k_inter effect on HR heterogeneity...' (R/sim_aft_gbsg.R:1246) appears in the deployed origin/gh-pages articles/paper_simulations.html:894.
    - count: paper_simulations.qmd 3 (1 evaluated, 2 text); all other files 0.
    - example: vignettes/articles/paper_simulations.qmd:683 `validation_results <- validate_k_inter_effect(` in evaluated chunk validate-k-inter (:678).
  - **print.gbsg_dgm() (R/sim_aft_gbsg.R:1157)**
    - facts: Registered as S3method(print,gbsg_dgm) (NAMESPACE:48). The literal name appears in no vignette or article. The method runs once, implicitly, from paper_simulations.qmd:523 `print(dgm_alt)` in evaluated chunk create-dgm-alt (:500). dgm_alt comes from setup_gbsg_dgm (:507), which sets class c('aft_dgm_flex','gbsg_dgm','list') (R/setup_gbsg_dgm.R:182). compute_dgm_cde (:520; R/oc_analyses.R:131) only sets fields. No print.aft_dgm_flex or other method for aft_dgm_flex exists in R/ or NAMESPACE, so dispatch reaches print.gbsg_dgm. The deployed origin/gh-pages articles/paper_simulations.html:523 holds its header 'GBSG-Based AFT Data Generating Mechanism (Aligned)' (R/sim_aft_gbsg.R:1158). No other print/str/summary of a dgm object and no bare dgm auto-print line exists in any .qmd/.Rmd. The class string 'gbsg_dgm' is mentioned as text at paper_simulations.qmd:383.
    - count: 0 literal mentions; 1 implicit evaluated S3 dispatch (paper_simulations.qmd:523); 1 text mention of the class (paper_simulations.qmd:383).
    - example: vignettes/articles/paper_simulations.qmd:523 `print(dgm_alt)` (evaluated; dispatches to print.gbsg_dgm).
  - **create_gbsg_dgm() (R/sim_aft_gbsg.R:199)**
    - facts: Exported (NAMESPACE:116) and deprecated (warning text at R/sim_aft_gbsg.R:215). No code call anywhere in vignettes/. The only mention is prose at biomarker_effects.qmd:66, outside all chunks (css :21-34; R chunks :171-189, :307-324, :561-585). The tracked rendered copy is biomarker_effects.html:2413 (text). Never executed by any vignette or article.
    - count: biomarker_effects.qmd 1 (text only); biomarker_effects.html 1 (rendered text); 0 code calls; 0 evaluated.
    - example: vignettes/articles/biomarker_effects.qmd:66 "`create_gbsg_dgm()` and `generate_aft_dgm_flex()` (their Supplementary" (prose).
  - **.create_gbsg_dgm_() (R/sim_aft_gbsg.R:236)**
    - facts: 0 mentions in vignettes/. Reached indirectly by all 7 evaluated article calls: the 4 setup_gbsg_dgm calls (:507, :614, :645, :699) via R/setup_gbsg_dgm.R:106; the 2 calibrate_k_inter calls (:604, :635) via eff_at at R/sim_aft_gbsg.R:1024; the validate_k_inter_effect call (:683) via R/sim_aft_gbsg.R:1254. Other R callers (R/sim_aft_gbsg.R:226 in create_gbsg_dgm, :1108 in get_dgm_with_output) are not reached from any article.
    - count: none (0 direct mentions); reached indirectly by 7 evaluated calls in paper_simulations.qmd.
    - example: none directly; indirect path: paper_simulations.qmd:507 setup_gbsg_dgm() -> R/setup_gbsg_dgm.R:106 .create_gbsg_dgm_().
  - **simulate_from_gbsg_dgm() (R/sim_aft_gbsg.R:789)**
    - facts: Exported (NAMESPACE:293) and deprecated (R/sim_aft_gbsg.R:804). 0 mentions in vignettes/. No caller in R/. The articles use simulate_from_dgm() (R/simulate_from_dgm.R:172) instead: paper_simulations.qmd:725 (evaluated, simulate-single :720), :477 (eval: false), treatment_effect_definitions.qmd:1075 (eval=FALSE). Never executed.
    - count: none
    - example: none
  - **.simulate_from_gbsg_dgm_() (R/sim_aft_gbsg.R:823)**
    - facts: 0 mentions. Its only caller is simulate_from_gbsg_dgm (R/sim_aft_gbsg.R:813), which is never called, so it is never executed.
    - count: none
    - example: none
  - **get_dgm_with_output() (R/sim_aft_gbsg.R:1076)**
    - facts: 0 mentions in vignettes/. @keywords internal (R/sim_aft_gbsg.R:1075), not exported (absent from NAMESPACE). It has man/get_dgm_with_output.Rd and a _pkgdown.yml reference entry (:339). No caller in R/; never executed. (get_dgm_hr in paper_simulations.qmd :1931, :2002-2038 is a different function.)
    - count: none
    - example: none
  - **cut_numeric() (R/sim_aft_gbsg.R:51)**
    - facts: 0 mentions in vignettes/; not exported. Runs only indirectly inside .create_gbsg_dgm_ at R/sim_aft_gbsg.R:281, :283, :284 (4-level quartile codes for age, pgr, nodes), so the 7 evaluated paper_simulations.qmd calls reach it. The eval:false 'equivalent' block at paper_simulations.qmd:395-440 does not reproduce it: it uses binary median splits (:413, :415, :416).
    - count: none
    - example: none
  - **cut_size() (R/sim_aft_gbsg.R:69)**
    - facts: 0 mentions in vignettes/; not exported. Runs only indirectly through .create_gbsg_dgm_ at R/sim_aft_gbsg.R:298 (3 levels at 20 and 50). The article's eval:false block uses a median split for size instead (paper_simulations.qmd:417).
    - count: none
    - example: none
  - **eff_at() (nested in calibrate_k_inter, R/sim_aft_gbsg.R:1023)**
    - facts: 0 mentions in vignettes/. Runs indirectly through the 2 evaluated calibrate_k_inter calls (paper_simulations.qmd:604, :635), which pass it to .calibrate_by_root() at R/sim_aft_gbsg.R:1040-1041. Same-named but separate helpers exist in R/find_k_inter_main.R:105, R/calibrate_k_treat.R:133 and R/calibrate_glm_interaction.R:148; they are not part of the wrapper.
    - count: none
    - example: none
  - **.rename() (nested in setup_gbsg_dgm, R/setup_gbsg_dgm.R:149)**
    - facts: 0 mentions in vignettes/. Runs indirectly through the 4 evaluated setup_gbsg_dgm calls; its call sites are R/setup_gbsg_dgm.R:155-159.
    - count: none
    - example: none
  - **Top-level constants SEED_BASE (R/sim_aft_gbsg.R:29), DAYS_PER_MONTH (:32), DEFAULT_N_SUPER (:35)**
    - facts: Non-function definitions in the wrapper file, missed by the investigator. 0 mentions in vignettes/ (git grep -w at c40d5309; the lowercase seed_base arguments are unrelated). Used at R/sim_aft_gbsg.R:205, :208 (create_gbsg_dgm defaults), :242, :245 (.create_gbsg_dgm_ defaults) and :273 (DAYS_PER_MONTH). Reached indirectly: DAYS_PER_MONTH by all 7 evaluated calls; SEED_BASE and DEFAULT_N_SUPER only through the calibrate_k_inter and validate_k_inter_effect calls (:604, :635, :683), since setup_gbsg_dgm passes its own seed and n_super (R/setup_gbsg_dgm.R:112, :115). The articles hard-code the same values without referencing the constants: 8316951 at paper_simulations.qmd:379, :438, :515; 5000 at :437, :513; 30.4375 at :406.
    - count: none (0 mentions); reached indirectly.
    - example: none
- **Files named:**
  - vignettes/forestsearch.qmd [current; the only CRAN-built vignette (vignettes/*.qmd, VignetteEngine quarto::html at :14-17; DESCRIPTION:30 VignetteBuilder: knitr, quarto); 0 mentions of any wrapper function or constant]
  - vignettes/articles/paper_simulations.qmd [current; pkgdown-only (.Rbuildignore:20 ^vignettes/articles$; _pkgdown.yml:56-57); setup_gbsg_dgm 21, calibrate_k_inter 8, validate_k_inter_effect 3, implicit print.gbsg_dgm 1 (:523), class string gbsg_dgm 1 (:383); all 7 evaluated wrapper calls are here; the mid-document YAML block :358-360 sets the rendered title]
  - vignettes/articles/treatment_effect_definitions.qmd [current; pkgdown-only (_pkgdown.yml:50-51); setup_gbsg_dgm 6, none evaluated as a call]
  - vignettes/articles/treatment_effect_definitions.html [tracked rendered output in the build-ignored articles dir; setup_gbsg_dgm on 8 lines, rendered text and code only]
  - vignettes/articles/extreme_subgroups.qmd [current; pkgdown-only (_pkgdown.yml:58-59); calibrate_k_inter 2, prose only; also calibrate_k_treat 8 (not a wrapper function)]
  - vignettes/articles/biomarker_effects.qmd [current; pkgdown-only (_pkgdown.yml:52-53); create_gbsg_dgm 1, prose only]
  - vignettes/articles/biomarker_effects.html [tracked rendered output; create_gbsg_dgm 1 at :2413, rendered text]
  - vignettes/articles/actg175-glm-outcomes.qmd [pkgdown-only; 0]
  - vignettes/articles/actg175-glm-outcomes.html [tracked rendered; 0]
  - vignettes/articles/causal_effects_brief_review.qmd [pkgdown-only; 0]
  - vignettes/articles/causal_effects_brief_review.html [tracked rendered; 0]
  - vignettes/articles/forestsearch-2024.Rmd [pkgdown-only; 0]
  - vignettes/articles/guohe-nested-continuum.qmd [pkgdown-only; 0]
  - vignettes/articles/methodology.qmd [pkgdown-only; 0]
  - vignettes/articles/mrct-consistency.qmd [pkgdown-only; 0]
  - vignettes/articles/mrct-consistency.html [tracked rendered; 0]
  - vignettes/articles/mrct-consistency-simulations.qmd [pkgdown-only; 0]
  - vignettes/articles/mrct-consistency-simulations.html [tracked rendered; 0]
  - vignettes/articles/survival-post-selection.qmd [pkgdown-only; 0]
  - vignettes/articles/survival-post-selection.html [tracked rendered; 0]
  - vignettes/articles/references.bib, vignettes/.gitignore, vignettes/articles/.gitignore [tracked, non-document; 0]
  - R/sim_aft_gbsg.R [source; constants :29, :32, :35; functions :51, :69, :199, :236, :789, :823, :1004 (nested eff_at :1023), :1076, :1157, :1229; cut_numeric/cut_size used at :281, :283, :284, :298]
  - R/setup_gbsg_dgm.R [source; :88, nested .rename :149; .create_gbsg_dgm_ call :106; class set :182]
  - NAMESPACE [evidence: :48 S3method(print,gbsg_dgm); exports :92 calibrate_k_inter, :116 create_gbsg_dgm, :286 setup_gbsg_dgm, :293 simulate_from_gbsg_dgm, :312 validate_k_inter_effect; get_dgm_with_output, cut_numeric, cut_size not exported; no aft_dgm_flex methods]
  - man/ [evidence: Rd pages exist for calibrate_k_inter, create_gbsg_dgm, cut_numeric, cut_size, get_dgm_with_output, print.gbsg_dgm, setup_gbsg_dgm, simulate_from_gbsg_dgm, validate_k_inter_effect]
  - .Rbuildignore [evidence: :20 ^vignettes/articles$]
  - _pkgdown.yml [evidence: article menu :36-68; reference entries :332-345]
  - DESCRIPTION [evidence: :30 VignetteBuilder: knitr, quarto]
  - R/oc_analyses.R [evidence: compute_dgm_cde :131 keeps class]
  - R/simulate_from_dgm.R [evidence: :172 simulate_from_dgm, used by the articles instead of simulate_from_gbsg_dgm]
  - R/fs_fdr_report.R [evidence: :256, :258, :262 the only R-level callers of setup_gbsg_dgm/calibrate_k_inter outside the wrapper; fs_fdr_report is not called in any vignette or article]
  - forestsearch origin/gh-pages:articles/paper_simulations.html [deployed render, d11f5514 from 19b703b (ancestor of c40d5309); evidence: <title> :8 and <h1> :70 taken from the mid-document YAML, print.gbsg_dgm output :523, validate_k_inter_effect output :894]
- **Notes:** Confirmed list: the 13 functions (11 named in the task, plus nested eff_at at R/sim_aft_gbsg.R:1023 and .rename at R/setup_gbsg_dgm.R:149). R/sim_aft_gbsg.R also defines three top-level constants: SEED_BASE (:29), DAYS_PER_MONTH (:32) and DEFAULT_N_SUPER (:35). None of the constants is mentioned in vignettes/; all are reached indirectly.

CRAN versus pkgdown: vignettes/forestsearch.qmd is the only CRAN-built vignette and mentions no wrapper function. Everything else is in vignettes/articles/, which .Rbuildignore:20 excludes, so it is pkgdown-only. The `vignette:` YAML in some articles has no effect on the CRAN build.

Evaluated use exists only in paper_simulations.qmd:
- 7 direct calls: setup_gbsg_dgm x4 (:507, :614, :645, :699), calibrate_k_inter x2 (:604, :635), validate_k_inter_effect x1 (:683).
- 1 implicit print.gbsg_dgm dispatch (:523).
- The deployed render confirms both the print output and the validate output.
- These calls run .create_gbsg_dgm_, cut_numeric, cut_size, eff_at, .rename and the constants indirectly.

Naming:
- 9 of the 13 functions are never named anywhere in vignettes/: .create_gbsg_dgm_, simulate_from_gbsg_dgm, .simulate_from_gbsg_dgm_, get_dgm_with_output, cut_numeric, cut_size, eff_at, .rename and print.gbsg_dgm.
- 4 are never executed even indirectly: create_gbsg_dgm (prose only), simulate_from_gbsg_dgm, .simulate_from_gbsg_dgm_ and get_dgm_with_output.

Other mentions:
- treatment_effect_definitions.qmd: one eval=FALSE call (:1034); the rest are prose, a plain-block comment and gt string literals.
- extreme_subgroups.qmd and biomarker_effects.qmd: prose only.
- The two tracked .html files carry the rendered copies.

Two documentation findings outside the counts:
- The mid-document YAML block at paper_simulations.qmd:358-360 replaces the article's rendered title on the deployed site.
- The article says its eval:false generate_aft_dgm_flex block (:395-440) is 'equivalent' to setup_gbsg_dgm (:392) and describes a 'binary discretisation' (:366). The wrapper actually uses quartile codes 1-4 from cut_numeric for age, pgr and nodes, and a 3-level cut_size for size (R/sim_aft_gbsg.R:281-284, :298). The article's block uses median splits instead, so the two are not equivalent.
- **Verifier (corrected) corrections:**
  - Definition list: the 13 functions are confirmed (11 top-level plus the nested eff_at at R/sim_aft_gbsg.R:1023 and .rename at R/setup_gbsg_dgm.R:149; a grep of every `function(` in both files finds exactly these 13 and no anonymous functions). But 'There are no other definitions in either file' is wrong. R/sim_aft_gbsg.R also defines three top-level constants: SEED_BASE <- 8316951L (:29), DAYS_PER_MONTH <- 30.4375 (:32) and DEFAULT_N_SUPER <- 5000L (:35). None is mentioned in any vignette or article (git grep -w at c40d5309 finds 0; the lowercase `seed_base` arguments in paper_simulations.qmd and extreme_subgroups.qmd are unrelated). All three are reached indirectly. DAYS_PER_MONTH is used in every .create_gbsg_dgm_() run (:273). SEED_BASE and DEFAULT_N_SUPER are .create_gbsg_dgm_() defaults (:245, :242), and they take effect through the evaluated calibrate_k_inter calls (paper_simulations.qmd:604, :635) and the validate_k_inter_effect call (:683): eff_at (:1024-1027) and :1254 pass no seed or n_super. They do not take effect through setup_gbsg_dgm, which passes its own defaults explicitly (R/setup_gbsg_dgm.R:94, :97, :112, :115).
  - The caveat on paper_simulations.qmd:358-360 is resolved, and :359 is not ordinary body text. The deployed pkgdown site (read-only `git show origin/gh-pages:articles/paper_simulations.html`, deploy d11f5514 built from 19b703b) shows the mid-document YAML block replacing the front-matter title from :2. The rendered page's <title> (:8), og:title (:15) and <h1> (:70) all read 'setup_gbsg_dgm() and generate_aft_dgm_flex(): Relationship and Usage', with the front-matter subtitle kept (:71). 19b703b is an ancestor of c40d5309, and paper_simulations.qmd differs between them only at :1340-1354 (a loop-variable rename). So :359 is a non-evaluated metadata mention that becomes the article's displayed title.
  - The notes say 'Seven wrapper functions have no direct mention' but list eight names: .create_gbsg_dgm_, simulate_from_gbsg_dgm, .simulate_from_gbsg_dgm_, get_dgm_with_output, cut_numeric, cut_size, eff_at and .rename. Adding print.gbsg_dgm, whose literal name also never appears (the investigator's own per-function entry says so), 9 of the 13 functions are never named in vignettes/. Only 4 are named: setup_gbsg_dgm, calibrate_k_inter, validate_k_inter_effect and create_gbsg_dgm.
  - The functions never executed by any vignette or article, directly or indirectly, are 4, not the 3 in the notes: create_gbsg_dgm (named only in prose at biomarker_effects.qmd:66), simulate_from_gbsg_dgm, .simulate_from_gbsg_dgm_ and get_dgm_with_output.
  - The setup_gbsg_dgm 'Total code calls 7, of which 4 are evaluated' does not match the investigator's own breakdown. There are 8 call expressions across files: 4 evaluated (paper_simulations.qmd:507, :614, :645, :699), 3 in eval=FALSE chunks (paper_simulations.qmd:375, :1756; treatment_effect_definitions.qmd:1034) and 1 in a plain ```r display block (paper_simulations.qmd:274). This excludes the two comment lines (paper_simulations.qmd:476; treatment_effect_definitions.qmd:460) and the two string literals (treatment_effect_definitions.qmd:552, :785). The total is 7 only if the plain-block call at :274 is dropped.
  - The cut_numeric note says paper_simulations.qmd:395-440 'writes the same discretisation out by hand'. That is wrong. The wrapper codes age, pgr and nodes into 4-level quartile codes 1-4 with cut_numeric (R/sim_aft_gbsg.R:51-58, used at :281, :283, :284) and enters them numerically in the AFT fit (:332, :403). It codes size into 3 levels at 20 and 50 with cut_size (:69-75, used at :298). The article's eval:false 'equivalent' block uses binary median splits for age, pgr, nodes and size (:413, :415, :416, :417). The article's statements that the call 'is equivalent' (:392) and that the wrapper uses 'the binary discretisation of continuous covariates' (:366) are therefore inaccurate against the source.
  - Minor completeness: extreme_subgroups.qmd also has indented plain ```r blocks at :218-220 and :224-229, which the investigator's list of nearby chunks omits. :202 and :485 are still prose outside every block. The file also names the non-wrapper sibling calibrate_k_treat() on 8 lines (defined in R/calibrate_k_treat.R:96, not one of the two wrapper files, and it does not call .create_gbsg_dgm_), so it is correctly excluded from the counts.
  - Added evidence, no change to the classification: the deployed article output shows the evaluated wrapper code ran. origin/gh-pages articles/paper_simulations.html:523 holds print.gbsg_dgm's header 'GBSG-Based AFT Data Generating Mechanism (Aligned)' (R/sim_aft_gbsg.R:1158), from print(dgm_alt) at qmd:523. :894 holds validate_k_inter_effect's verbose header (R/sim_aft_gbsg.R:1246). The only R-level callers of wrapper functions outside the two files are fs_fdr_report() (R/fs_fdr_report.R:256, :258, :262), which no vignette or article calls, and a message string in R/fs_dgm_feasibility.R:130. No vignette or article uses inline `r` code, child includes, source() or readRDS that reaches the wrapper.

### L-D

- **Method:** Independent re-derivation, read-only throughout, with scratch files only in /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2ro/LDv/. forestsearch HEAD is c40d530935e88230731ec280a9dfd07159a3fea6 on cran-prep. git status shows only one untracked directory, dev/reports/paper_sims_rerun_2026-11-11b/, which is outside quarto/. git ls-files quarto and ls-tree c40d5309 both give 6466 entries.

(1) Definitions. I grepped t2_tree for `function\s*\(`, for top-level `<-`/`=`, and for lambda, setMethod and local() forms. The definitions in R/sim_aft_gbsg.R are:
- 3 top-level constants: SEED_BASE :29, DAYS_PER_MONTH :32, DEFAULT_N_SUPER :35.
- 10 functions: cut_numeric :51, cut_size :69, create_gbsg_dgm :199, .create_gbsg_dgm_ :236, simulate_from_gbsg_dgm :789, .simulate_from_gbsg_dgm_ :823, calibrate_k_inter :1004, get_dgm_with_output :1076, print.gbsg_dgm :1157, validate_k_inter_effect :1229.
- 1 nested helper: eff_at :1023.
R/setup_gbsg_dgm.R defines setup_gbsg_dgm :88 and the nested helper .rename :149. There are no other definitions and no anonymous functions.
NAMESPACE has S3method(print,gbsg_dgm) :48 and exports calibrate_k_inter :92, create_gbsg_dgm :116, setup_gbsg_dgm :286, simulate_from_gbsg_dgm :293 and validate_k_inter_effect :312. No print, summary or format method is registered for aft_dgm_flex.

(2) Scope. I filtered git ls-files -z quarto to the .qmd, .R, .r, .Rmd, .rmarkdown and .md extensions: 1111 files (663 qmd, 232 R, 2 r, 5 Rmd, 6 rmarkdown, 203 md).

(3) Classification. I listed every directory holding in-scope files (73 directories) and read them by hand. The investigator's whole-token rule (archive[sd]?|legacy|broken|old|superseded|deprecated|obsolete) gives 189 archive and 922 current files, exactly as claimed. A looser substring rule (archiv, legacy, broken, old, supersed, deprecat, obsolet, bak, backup, retired, unused, trash, attic, stale, prev) flags no further directory. I also checked the 281 current calling files: no base name and no first-60-line header carries an archive, legacy or superseded marker. The only notable names are TEST_maxeffCons_..._standalone2.qmd and one committed .rmarkdown render intermediate next to its .qmd.

(4) Scanner. My own Python scanner (vscan.py) uses a different fence state machine from the investigator's. It also tracks YAML front matter, any-position `#| eval:` options and chunk-header eval=F. It classified 815 hits in all.

(5) Extra checks, each of which came back clean:
- No document-level eval disable: no YAML `execute: eval: false` and no knitr::opts_chunk$set(eval=...).
- All call lines sit in 12 chunk labels, e.g. {r build-dgm} with 543 call lines, and none of those chunks has an eval, include or purl option. The `#| eval: !expr run_fb` and `run_tier1` options in some calling files belong to other chunks (around lines 826 and 1261-1331).
- No call line is inside an HTML comment.
- No inline `r` code, no forestsearch::/::: prefixed calls and no line-split calls.
- No function-object references: no do.call, match.fun, get, exists or bare name passed to lapply. The only bare non-paren code-region hit is the string in gate2G.R:216.

(6) Cross-checks against the investigator's files.
- Their files_calls.txt list of 294 calling files is identical to mine (diff empty).
- All 581 per-file (function, count, first line) tuples match: 294 setup_gbsg_dgm + 287 calibrate_k_inter.

(7) print.gbsg_dgm dispatch. setup_gbsg_dgm.R:182 sets class c("aft_dgm_flex","gbsg_dgm","list"), and .create_gbsg_dgm_ sets c("gbsg_dgm","list") at :735. I collected every assigned variable: dgm 277, dgm_null 9, dgm_alt 5, dgm_gbsg_null 3, plus alt_dgms[[j]], dgm_h1, dgm_s, dgm0, d, dgm_surv_null, dgm_harm, dgm_hom. I searched for bare-line auto-print, print(), summary(), show(), str() and `(VAR <- ...)` in executed code of all 294 calling files: 0 hits.
Of the 314 setup_gbsg_dgm call lines (297 + 17), 306 are plain assignments. The other 8 are function bodies whose results are assigned or $-indexed (.hr_at, hr_at, create_dgm_h1) or are nested in compute_dgm_cde(...); none is a top-level bare call that would auto-print.

(8) Sweep of every extension. Exact-name git grep -F finds zero occurrences of get_dgm_with_output, validate_k_inter_effect (and k_inter_effect), print.gbsg, cut_numeric, cut_size, eff_at, DAYS_PER_MONTH and DEFAULT_N_SUPER. It also finds zero exact create_gbsg_dgm (see correction 1). simulate_from_gbsg_dgm appears in 1 file only. Outside the html files, the out-of-scope files that mention setup_gbsg_dgm are the 3 log .txt files.
- **Per function:**
  - **create_gbsg_dgm (R/sim_aft_gbsg.R:199; exported NAMESPACE:116)**
    - facts: No calls and no exact-name occurrences in any of the 6466 tracked quarto files, any extension. Checked with the R-identifier-boundary regex and with git grep -w. Every partial-string hit (943 lines in 163 files) is .create_gbsg_dgm_.
    - count: (a) 0 files / 0 call lines, 0 mentions. (b) 0 / 0.
    - example: none
  - **.create_gbsg_dgm_ (R/sim_aft_gbsg.R:236; internal)**
    - facts: No calls. All occurrences are text mentions in current (a) files.
- 13 prose lines in 5 files: gbsg_aft_dgm_appendix.qmd :49, :51, :94; REPORT_null_gbsg_identification_2026-09-21.md :94, :231, :293; REPORT_p30_stage0_2026-09-06.md :14, :36, :104; current_status.md :18, :191; status_curated.md :13, :186.
- 6 comment lines in 4 files: declcal_c0approx_run.R :44, declcal_run.R :27, declcalc0_run.R :43, sim_fs_maxeffCons_fb_mr_field_m1_template.qmd :380, :811, :941.
    - count: (a) calls: 0 / 0; mentions: 9 files / 19 lines. (b) 0 / 0.
    - example: (a) prose: quarto/calibration/gbsg_aft_dgm_appendix.qmd:49 "`.create_gbsg_dgm_()`, then reshapes the result into an `aft_dgm_flex` object". (a) comment: quarto/simulations/gbsg_020/scripts_dinamr/declcal_run.R:27 "#       .create_gbsg_dgm_() requires k_treat > 0 (:256).". (b): none.
  - **simulate_from_gbsg_dgm (R/sim_aft_gbsg.R:789; exported NAMESPACE:293)**
    - facts: No calls. The only file with the string in any extension is one archive file, with 2 comment or roxygen lines.
    - count: (a) 0 / 0. (b) calls 0 / 0; mentions 1 file / 2 lines.
    - example: (b) quarto/legacy calibration files/run_simulation_analysis.R:6 "# to simulate_from_gbsg_dgm() and hardcoded column names (y.sim, event.sim,"; also :76 (roxygen).
  - **.simulate_from_gbsg_dgm_ (R/sim_aft_gbsg.R:823; internal)**
    - facts: No occurrences in quarto, any extension.
    - count: (a) 0 / 0. (b) 0 / 0.
    - example: none
  - **calibrate_k_inter (R/sim_aft_gbsg.R:1004; exported NAMESPACE:92)**
    - facts: Real calls are in {r} chunks and .R scripts. All 277 current calling files also call setup_gbsg_dgm.
- (a) call files by extension: 269 .qmd, 7 .R, 1 .rmarkdown.
- (a) files/lines by directory: calibration 3/3, simulations top level 1/1, gbsg 105/105, gbsg/sims1_claude_summary 3/3, gbsg_020 74/74, gbsg_020/scripts_dinamr 6/6, gbsg_redux 85/90.
- (b) call files: _broken/gbsg_analysis.qmd :1368, archive/reproduce_table1_m1_working.qmd :225, legacy calibration files fourway :246 and frontier :344, and 6 gbsg_redux/legacy files at :249.
- (a) display-only ```r: REPORT_p30_stage0_2026-09-06.md :31, :96; _oracle_cells_kappa.md :46.
- (a) text: 33 prose lines in 18 files, plus 3 chunk-comment lines in the field template at :379, :387, :807.
- (b) text: 9 prose lines (actg175_simulations_benefit :861, :1535; ddi-vs-zdv :764, :1422; review.md :26, :296, :318 in a ```markdown fence; survival_benefit_working :915, :1993) and 2 comment lines (_broken/gbsg_analysis.qmd :190; archive/actg175_simulations_1.qmd :452).
    - count: (a) real calls: 277 files / 282 lines (+2 files / 3 lines display-only; 19 files / 36 lines mentions). (b) real calls: 10 files / 10 lines (0 display; 6 files / 11 lines mentions).
    - example: (a) call: quarto/calibration/calibration_survival_fourway.qmd:260 "k_inter_cox <- calibrate_k_inter(". (b) call: quarto/applications/gbsg/_broken/gbsg_analysis.qmd:1368 "kj   <- calibrate_k_inter(". (a) mention: quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_field_m1_template.qmd:379 (comment). (b) mention: quarto/archive/actg175_simulations_benefit.qmd:861 (prose).
  - **get_dgm_with_output (R/sim_aft_gbsg.R:1076; not exported)**
    - facts: No occurrences in quarto, any extension; the partial string dgm_with_output also has 0 hits.
    - count: (a) 0 / 0. (b) 0 / 0.
    - example: none
  - **print.gbsg_dgm (R/sim_aft_gbsg.R:1157; S3method(print,gbsg_dgm) NAMESPACE:48)**
    - facts: No explicit calls and no text mentions (the partial print.gbsg has 0 hits). Implicit dispatch would be possible: setup_gbsg_dgm.R:182 sets class c("aft_dgm_flex","gbsg_dgm","list"), and no print method is registered for aft_dgm_flex.
I checked all 294 calling files in executed code for print(VAR), summary(VAR), str(VAR), a bare VAR line, or `(VAR <- ...)`, where VAR is assigned from setup_gbsg_dgm or from create_dgm_h1: 0 hits. None of the 8 non-assignment call lines is a top-level bare call that would auto-print.
    - count: (a) 0 / 0, explicit or detectable dispatch. (b) 0 / 0.
    - example: none
  - **validate_k_inter_effect (R/sim_aft_gbsg.R:1229; exported NAMESPACE:312)**
    - facts: No occurrences in quarto, any extension; the partial string k_inter_effect also has 0 hits.
    - count: (a) 0 / 0. (b) 0 / 0.
    - example: none
  - **cut_numeric (R/sim_aft_gbsg.R:51; not exported)**
    - facts: No occurrences in quarto, any extension.
    - count: (a) 0 / 0. (b) 0 / 0.
    - example: none
  - **cut_size (R/sim_aft_gbsg.R:69; not exported)**
    - facts: No occurrences in quarto, any extension.
    - count: (a) 0 / 0. (b) 0 / 0.
    - example: none
  - **setup_gbsg_dgm (R/setup_gbsg_dgm.R:88; exported NAMESPACE:286)**
    - facts: Real calls include nested and function-body forms: nullid_designpoint.R:22 (hr_at <- function(kt) setup_gbsg_dgm(...)), nullid_designpoint.R:26 (compute_dgm_cde(setup_gbsg_dgm(...))), the .hr_at bodies at declcal_run.R:140, declcal_c0approx_run.R:157, declcalc0_run.R:156 and field template :821, and the create_dgm_h1 return value at selection_survival_frontier.qmd:380 and its legacy copy :349. 306 of the 314 call lines are plain assignments.
- (a) call files by extension: 272 .qmd, 8 .R, 1 .rmarkdown.
- (a) files/lines by directory: calibration 6/9, qc 1/1, simulations top level 1/2, gbsg 105/106, gbsg/sims1_claude_summary 3/3, gbsg_020 74/75, gbsg_020/scripts_dinamr 6/11, gbsg_redux 85/90.
- (b) files/lines: _broken 3/4, archive 1/2, legacy calibration files 3/5, gbsg_redux/legacy 6/6.
- Display-only ```r calls: (a) REPORT_p30_stage0 :32, :99 and _oracle_cells_kappa.md :48; (b) COMMIT_SUMMARY.md :91.
- (a) text: 34 prose lines in 18 files; 48 comment lines in 20 files, including cox_causal_review.qmd:1241, which is a comment inside a display ```r block; 3 strings in 3 files (gate2G.R:216, nullid_designpoint.R:16, cox_causal_review.qmd:1312).
- (b) text: 34 prose lines in 10 files, including COMMIT_SUMMARY.md:27 in a text fence; 7 comment lines in 4 files; 4 strings in 3 _broken files (gbsg_analysis.qmd :1664, :2220; 2sim-approaches :1291; working2 :1082).
    - count: (a) real calls: 281 files / 297 lines (+2 files / 3 lines display-only; 38 files / 85 lines mentions). (b) real calls: 13 files / 17 lines (+1 file / 1 line display-only; 11 files / 45 lines mentions).
    - example: (a) call: quarto/calibration/calibration_survival_fourway.qmd:269 "dgm_alt <- setup_gbsg_dgm(". (b) call: quarto/applications/gbsg/_broken/gbsg_analysis.qmd:1082 "dgm_gbsg_null <- setup_gbsg_dgm(". (a) mention: quarto/calibration/selection_survival_frontier.qmd:372 (comment). (b) mention: quarto/applications/gbsg/_broken/gbsg_analysis.qmd:94 (prose).
  - **eff_at (nested in calibrate_k_inter, R/sim_aft_gbsg.R:1023)**
    - facts: A closure that cannot be reached from outside the function. 0 occurrences in quarto (git grep -w).
    - count: (a) 0 / 0. (b) 0 / 0.
    - example: none
  - **.rename (nested in setup_gbsg_dgm, R/setup_gbsg_dgm.R:149)**
    - facts: A closure. There are no identifier-bounded occurrences. The only substring hit is quarto/_utils/gsub_in_files.R:136 file.rename(f, new_path), which is base R and unrelated.
    - count: (a) 0 / 0. (b) 0 / 0.
    - example: none
  - **SEED_BASE (top-level constant, R/sim_aft_gbsg.R:29; not exported)**
    - facts: There is no reference to the package constant. The one file with hits defines its own SEED_BASE locally, with the same value, and uses it 3 times.
    - count: (a) 0 package references; 1 file with a local redefinition (1 assignment + 3 uses). (b) 0 / 0.
    - example: quarto/simulations/gbsg_020/scripts_dinamr/grf_mechanism.R:34 "SEED_BASE      <- 8316951L"; used at :78, :83, :90.
  - **DAYS_PER_MONTH (top-level constant, R/sim_aft_gbsg.R:32; not exported)**
    - facts: No occurrences in quarto, any extension.
    - count: (a) 0 / 0. (b) 0 / 0.
    - example: none
  - **DEFAULT_N_SUPER (top-level constant, R/sim_aft_gbsg.R:35; not exported)**
    - facts: No occurrences in quarto, any extension.
    - count: (a) 0 / 0. (b) 0 / 0.
    - example: none
- **Files named:**
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/R/sim_aft_gbsg.R [wrapper source, definitions]
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/R/setup_gbsg_dgm.R [wrapper source, definitions; class set at :182]
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/NAMESPACE [exports :92 :116 :286 :293 :312; S3method :48]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/calibration/calibration_survival_fourway.qmd [current (a); calibrate @260; setup @269, :275]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/calibration/calibration_survival_leff_grid.qmd [current (a); setup @346]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/calibration/fpr_calibration_demo.qmd [current (a); setup @363]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/calibration/fs_fdr_report_smoke.qmd [current (a); setup @76]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/calibration/gbsg_aft_dgm_appendix.qmd [current (a); setup @192, :195; calibrate @193; .create_gbsg_dgm_ prose :49, :51, :94; source-file text :42-43]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/calibration/selection_survival_frontier.qmd [current (a); setup @357, :380 (create_dgm_h1 body); calibrate @375; comment :372]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/qc/qc_null_classification.R [current (a); setup @31]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/reproduce_table1_m1.qmd [current (a); calibrate @199; setup @213, :225]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg/dina_failure_diagnostic_standalone.R [current (a); calibrate @62; setup @64]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg/fs-eff_psi_t1-vs-t2_m1_h15_knoise0_n2000_batched_1_20.rmarkdown [current (a); render intermediate of the same-named .qmd; calibrate @190; setup @193]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg/sims1_claude_summary/dina-eff_psi_t2-only_m1_h15_knoise5_n500.qmd [current (a); calibrate @190; setup @193]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg/sims1_claude_summary/fs-eff_psi_t2-only_m1_h15_knoise5_n500.qmd [current (a); calibrate @190; setup @193]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg/sims1_claude_summary/grf-eff_psi_t2-only_m1_h15_knoise5_n500.qmd [current (a); calibrate @190; setup @193]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020/scripts_dinamr/declcal_c0approx_run.R [current (a); calibrate @149; setup @157, :164; .create comment :44]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020/scripts_dinamr/declcal_preflight.R [current (a); calibrate @30; setup @31]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020/scripts_dinamr/declcal_run.R [current (a); calibrate @132; setup @140, :147; .create comment :27]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020/scripts_dinamr/declcalc0_run.R [current (a); calibrate @148; setup @156, :163; .create comment :43]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020/scripts_dinamr/grf_mechanism.R [current (a); calibrate @75; setup @77; local SEED_BASE :34]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020/scripts_dinamr/nullid_designpoint.R [current (a); setup @22, :26, :62; calibrate @61; string :16; source-file text :18]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_field_m1_template.qmd [current (a); calibrate @816; setup @821, :829; comments :379-388, :807, :811, :941]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg [current (a) directory: 105 calling files]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020 [current (a) directory: 74 calling files at top level]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_redux [current (a) directory: 85 calling files at top level]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020/REPORT_p30_stage0_2026-09-06.md [current (a); display-only ```r :31, :32, :96, :99; prose]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_redux/_oracle_cells_kappa.md [current (a); display-only ```r :46, :48]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020/scripts_dinamr/gate2G.R [current (a); string :216 only]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/guides/cox_causal_review.qmd [current (a); prose :1234, display-block comment :1241, string :1312]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020/REPORT_null_gbsg_identification_2026-09-21.md [current (a); mentions only; .create prose :94, :231, :293]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020/current_status.md [current (a); mentions only; .create prose :18, :191]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020/status_curated.md [current (a); mentions only; .create prose :13, :186]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020/summary_dinamr.qmd [current (a); calling file; source-file text :942]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/claude code notes/CLAUDE.md [current (a); mention only; source-file text :59]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/HANDOFF_1_fdr_reporting_tool.md [current (a); prose :42, :54]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/applications/gbsg/_broken/gbsg_analysis.qmd [archive (b); setup @1082, :1374; calibrate @1368; prose :94; comment :190, :219; strings :1664, :2220]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/applications/gbsg/_broken/gbsg_analysis_2sim-approaches.qmd [archive (b); setup @755; string :1291]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/applications/gbsg/_broken/gbsg_analysis_working2.qmd [archive (b); setup @743; string :1082]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/archive/reproduce_table1_m1_working.qmd [archive (b); calibrate @225; setup @239, :251]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/legacy calibration files/calibration_survival_fourway.qmd [archive (b); calibrate @246; setup @255, :266]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/legacy calibration files/calibration_survival_leff_grid.qmd [archive (b); setup @129]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/legacy calibration files/selection_survival_frontier.qmd [archive (b); calibrate @344; setup @326, :349]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_redux/legacy [archive (b); 6 fs_t1_t2_m1_*_combine_*.qmd, each calibrate @249, setup @252]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/archive/COMMIT_SUMMARY.md [archive (b); text :27 in a plain ``` text fence; display-only ```r :91]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/archive/actg175_simulations_review.md [archive (b); calibrate prose :26, :296, :318 (in a ```markdown fence)]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/archive/actg175_simulations_benefit.qmd [archive (b); mention only; calibrate prose :861, :1535]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/legacy calibration files/run_simulation_analysis.R [archive (b); simulate_from_gbsg_dgm comments :6, :76; setup roxygen :83, :127, :134]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020/scripts_dinamr/logs/gate2G_B.txt [out-of-scope extension; 6 setup_gbsg_dgm log lines]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020/scripts_dinamr/logs/gate2G_C.txt [out-of-scope extension; 6 setup_gbsg_dgm log lines]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020/scripts_dinamr/logs/nullid_designpoint.txt [out-of-scope extension; 1 setup_gbsg_dgm log line]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/_utils/gsub_in_files.R [current; file.rename :136, unrelated]
  - /home/larryleon/Documents/GitHub/forestsearch/quarto/qc/gbsg_cv_qc.qmd [current; assignInNamespace :638 on run_single_consistency_split, unrelated]
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2ro/LDv/vscan.py [scratch: independent scanner]
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2ro/LDv/vhits.json [scratch: all 815 classified hits]
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2ro/LDv/callfiles.txt [scratch: 294 calling files, identical to the investigator's]
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2ro/LD/files_calls.txt [investigator scratch: 294 calling files with counts; 581 tuples verified]
- **Notes:** VERDICT: all counts, file:line examples, directory classifications and the 294-file calling list are confirmed. The two corrections concern secondary statements in the method and caveats, not any count.

HEADLINE (confirmed). In quarto/, only setup_gbsg_dgm and calibrate_k_inter are actually called.
- Current (a): 281 files call a wrapper (272 .qmd, 8 .R, 1 .rmarkdown). All 281 call setup_gbsg_dgm (297 call lines); 277 of them also call calibrate_k_inter (282 lines). The 4 that call setup only are quarto/calibration/{calibration_survival_leff_grid,fpr_calibration_demo,fs_fdr_report_smoke}.qmd and quarto/qc/qc_null_classification.R.
- Archive (b): 13 files call a wrapper, all .qmd: setup_gbsg_dgm 13 files / 17 lines, calibrate_k_inter 10 files / 10 lines. The 3 that call setup only are _broken/gbsg_analysis_2sim-approaches.qmd, _broken/gbsg_analysis_working2.qmd and legacy calibration files/calibration_survival_leff_grid.qmd.
- No calls in (a) or (b) to create_gbsg_dgm, .create_gbsg_dgm_, simulate_from_gbsg_dgm, .simulate_from_gbsg_dgm_, get_dgm_with_output, validate_k_inter_effect, cut_numeric, cut_size, print.gbsg_dgm (neither explicit nor detectable dispatch), eff_at, .rename or the three constants.

CLASSIFICATION RULE (confirmed). A file is archive (b) if any directory component matches, case-insensitively, (^|[^a-z])(archive[sd]?|legacy|broken|old|superseded|deprecated|obsolete)([^a-z]|$). This gives 922 current and 189 archive files. The (b) directories, with in-scope file counts:
- quarto/archive (13)
- quarto/applications/actg175/_archive (20)
- quarto/applications/gbsg/_archive (16)
- quarto/applications/gbsg/_broken (16)
- quarto/legacy calibration files (9)
- quarto/simulations/_legacy (16)
- quarto/simulations/actg175/_archive_pre_fsparams (3)
- quarto/simulations/actg175/binary_legacy (89, including its _fsparams 16 and _summaries 1)
- quarto/simulations/gbsg_redux/legacy (7)

CURRENT (a) CALLING FILES BY DIRECTORY (confirmed; files/lines as calibrate | setup):
- calibration: 3/3 | 6/9. Files: calibration_survival_fourway, calibration_survival_leff_grid, fpr_calibration_demo, fs_fdr_report_smoke, gbsg_aft_dgm_appendix, selection_survival_frontier.
- qc: none | 1/1 (qc_null_classification.R).
- simulations top level: 1/1 | 1/2 (reproduce_table1_m1.qmd).
- simulations/gbsg: 105/105 | 105/106. By prefix: betaHhat_sweep 2, dina 15, dina-eff_psi 8, dina_failure_diagnostic_standalone.R 1, fs 35, fs-eff_psi 14 (including the .rmarkdown), fsi 1, fs_psi 2, fs_subgroup_id_sweep 4, grf 5, grf-eff_psi 5, methods_estimation_sweep 1, mr_coverage_sweep 11, reproduce_table1 1.
- simulations/gbsg/sims1_claude_summary: 3/3 | 3/3.
- simulations/gbsg_020: 74/74 | 74/75. By prefix: gate_d2_cim_{field,ij,unset} 3, maxeffCons_mr_coverage_sweep 15, TEST_maxeffCons_mr_coverage_sweep 1, mr_coverage_sweep 1, sim_dina_* 7, sim_fs_maxeffCons_fb_mr 35 + field template 1, sim_fs_maxeffCons_mr 7, sim_grf_eff_fb_mr 4.
- simulations/gbsg_020/scripts_dinamr: 6/6 | 6/11. Files: declcal_c0approx_run.R, declcal_preflight.R, declcal_run.R, declcalc0_run.R, grf_mechanism.R, nullid_designpoint.R.
- simulations/gbsg_redux: 85/90 | 85/90. By prefix: fs_t1_t2 43, fsNEW_t1_t2 1, dina_t1_t2 4, grf_t1_t2 6, mr_coverage_sweep 8, eff{Max,Min}SG/maxcons_mr_coverage_sweep 3, sim_fs_effMaxSG_fb_mr 5, sim_fs_effMinSG_fb_mr 6, sim_fs_maxcons_fb_mr 5, sim_fs_maxeffCons_fb_mr 4.
The full list of 294 files, with per-file counts and first lines, is in /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2ro/LD/files_calls.txt. My independent copy is /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2ro/LDv/callfiles.txt; the two are identical.

ARCHIVE (b) CALLING FILES (confirmed):
- _broken/gbsg_analysis.qmd: setup :1082, :1374; calibrate :1368.
- _broken/gbsg_analysis_2sim-approaches.qmd: setup :755.
- _broken/gbsg_analysis_working2.qmd: setup :743.
- archive/reproduce_table1_m1_working.qmd: calibrate :225; setup :239, :251.
- legacy calibration files/calibration_survival_fourway.qmd: calibrate :246; setup :255, :266.
- legacy calibration files/calibration_survival_leff_grid.qmd: setup :129.
- legacy calibration files/selection_survival_frontier.qmd: calibrate :344; setup :326, :349.
- simulations/gbsg_redux/legacy: 6 fs_t1_t2_m1_*_combine_*.qmd files, each with calibrate :249 and setup :252.
- These (b) directories have no hits of any kind: applications/actg175/_archive, applications/gbsg/_archive, simulations/_legacy, simulations/actg175/_archive_pre_fsparams, simulations/actg175/binary_legacy.

MENTION-ONLY OR DISPLAY-ONLY FILES (confirmed):
- Current (a), 23 files: quarto/HANDOFF_1_fdr_reporting_tool.md, HANDOFF_2_companion_paper.md, claude code notes/{CLAUDE.md, claude_code_guide.qmd, claude_md_how_to_use.qmd}, extreme_subgroups/{gbsg_survival_fixed, gbsg_survival_resample}.qmd, guides/{cox_causal_review.qmd, extreme_subgroups.qmd, forestsearch_glm_statistical_review.md}, presentations/my use of claude/extreme_subgroups.qmd, simulations/actg175/survival/actg175_survival_benefit_simulations.qmd, simulations/gbsg_020/{REPORT_dinamr_blockC_2026-09-11, REPORT_grfmr_2026-09-11, REPORT_mr_field_ocmap_stage0_2026-09-05, REPORT_null_gbsg_identification_2026-09-21, REPORT_o1_forwarding_2026-09-10, REPORT_p30_stage0_2026-09-06, REPORT_p30_stage1_2026-09-06, current_status, status_curated}.md, simulations/gbsg_020/scripts_dinamr/gate2G.R, simulations/gbsg_redux/_oracle_cells_kappa.md.
- Archive (b), 7 files: quarto/archive/{COMMIT_SUMMARY.md, actg175_simulations_1.qmd, actg175_simulations_benefit.qmd, actg175_simulations_ddi-vs-zdv.qmd, actg175_simulations_review.md, actg175_survival_benefit_simulations_working.qmd} and legacy calibration files/run_simulation_analysis.R.
- Two fence edge cases, both classified correctly by the investigator as text: COMMIT_SUMMARY.md:27 sits inside a plain ``` commit-message text fence (lines 7-30), and actg175_simulations_review.md:318 sits inside a ```markdown fence. The only display-only ```r code calls are REPORT_p30_stage0:31, :32, :96, :99, _oracle_cells_kappa.md:46, :48 and COMMIT_SUMMARY.md:91.

CAVEATS:
1. Some (b) files are near-copies of (a) files. The three .qmd files in "legacy calibration files" duplicate quarto/calibration/, and archive/reproduce_table1_m1_working.qmd duplicates simulations/reproduce_table1_m1.qmd.
2. No quarto file redefines or sources a wrapper. Many files call devtools load_all (334 non-html files mention it), so calls resolve either to the installed export or to the dev R/ source; both are the same wrapper.
3. SEED_BASE in grf_mechanism.R:34 (used at :78, :83, :90) is a local constant with the same name and value as the package's, not a reference to the package constant.
4. quarto/qc/gbsg_cv_qc.qmd:638 uses assignInNamespace on run_single_consistency_split, which is unrelated to the wrapper.
- **Verifier (corrected) corrections:**
  - Method step 8 (sanity sweep) is wrong where it says the only exact create_gbsg_dgm hits in any extension are 3 .txt files. create_gbsg_dgm has 0 exact hits in all 6466 tracked quarto files, by both the R-identifier-boundary regex and git grep -w. All 943 partial-string hits, in 163 files (154 .html, 4 .md, 2 .qmd, 3 .R), are the internal .create_gbsg_dgm_. The 3 .txt files are quarto/simulations/gbsg_020/scripts_dinamr/logs/gate2G_B.txt (6 lines), gate2G_C.txt (6) and nullid_designpoint.txt (1). They contain only setup_gbsg_dgm (log output) and are outside the requested extensions.
  - Caveat 2 is wrong where it says the source-file names sim_aft_gbsg.R and setup_gbsg_dgm.R appear as text only in gbsg_aft_dgm_appendix.qmd:42 and REPORT_p30_stage0_2026-09-06.md:4. They appear as text on 29 lines in 13 in-scope files. Examples: quarto/claude code notes/CLAUDE.md:59, quarto/simulations/gbsg_020/REPORT_null_gbsg_identification_2026-09-21.md:94, quarto/simulations/gbsg_020/scripts_dinamr/nullid_designpoint.R:18 (a cat string), quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_field_m1_template.qmd:382 and :388 (comments), declcal_run.R:26 (comment), quarto/simulations/gbsg_020/summary_dinamr.qmd:942. The substantive conclusion stands: no quarto file source()s either file (0 hits for source\(.*(sim_aft|setup_gbsg)), and none redefines a wrapper function.

### L-E

- **Method:** Wrapper list: I listed every 'function' and '\\(' definition in t2_tree/R/sim_aft_gbsg.R (1316 lines) and R/setup_gbsg_dgm.R (190 lines). The result is the 11 named functions plus the nested eff_at (:1023) and .rename (:149), and nothing else; the list is confirmed. I cross-checked against NAMESPACE (export and S3method lines) and man/*.Rd aliases. Rd pages exist for all except .create_gbsg_dgm_ and .simulate_from_gbsg_dgm_.

fs-post-selection (main @ 5d165ac, read-only), search routes that differ from the investigator's PCRE:
(a) git grep -I -w -F per name, and partial patterns (gbsg_dgm, create_gbsg, calibrate_k, k_inter_effect, dgm_with_output, cut_num, eff_at, .rename, forestsearch:::).
(b) LaTeX-escaped forms: gbsg\\_dgm, k\\_inter. This found the 6 missed .tex files.
(c) Loose separators, for HTML- or markdown-split names. None found.
(d) Every tracked pdf, docx and zip (36 files) streamed via git cat-file into scratch, then pdftotext, unzip -p word/*.xml, and unzip -p on zip members.
(e) All 172 rds, all gzip (magic 1f8b), decompressed with zcat and grepped. No hits.

Render set: built from the include directives in main.qmd (:710, :747) and supplement.qmd (:64 to :1170; :1174-1178 are escaped {{{< >}}} and inactive), plus the documented per-file render commands (REPORT_rev1_averaging_text:97-99). Payloads were taken from the readRDS and run_tag lines in each included fragment.

Generators: cross-checked against the inventory CSV (39 render payloads: 5 AFFECTED with a generator, 34 NOT AFFECTED with none), payload meta$source_files (DINA and GRF) and the run_tag/campaign_tag formulas. In forestsearch (cran-prep @ c40d5309) I grepped each of the 37 generators at HEAD and at every historical commit of each file (git log plus git show), checked the run_mode and run_sweep guards, and confirmed the cited blobs (fded7718, 114aa11d, fc59934b, f0121563, e9ba13a8, 74b0e8eb, 8dc86a4b) and their hit counts. I checked generators for source(), load_all and forestsearch::: (only .fs_rule_columns and tidy_cut_display, neither a wrapper function). I traced package-internal reach: across R/ outside the wrapper files, the only wrapper calls are fs_fdr_report.R:256-262, and no generator calls fs_fdr_report. I re-derived the payload_manifest gbsg_redux count by expanding its generator patterns: 65 documents, all calling the wrapper. Outside manuscript_rev1 I counted both mention lines and code call lines by top directory. One Rscript --vanilla call read payload meta (under 60 s).
- **Per function:**
  - **setup_gbsg_dgm (R/setup_gbsg_dgm.R:88)**
    - facts: Exported (NAMESPACE:286). It calls .create_gbsg_dgm_ at setup_gbsg_dgm.R:106 and sets class c("aft_dgm_flex","gbsg_dgm","list") at :182. In manuscript_rev1 (tracked text) it appears in 7 files on 11 lines. (1) Six executable call lines, at :156 and :230 in each of dgm_truth_targets.qmd, dgm_truth_targets_h10.qmd and dgm_truth_targets_h15.qmd; these are standalone documents that neither the main paper nor the supplement includes. manuscript_render_manifest.html:666/:669-671 lists dgm_truth_targets.qmd and _h15 as auxiliary and does not name _h10. (2) Two rendered-code lines in dgm_truth_targets_h10.html, at :2794 and :3745. (3) Three comment lines, none of them calls: _sim_estimation_fs_h15_n1000.qmd:124 (included at supplement.qmd:1144), _sim_noise_null_sweep_standalone2.qmd:63 (included at :1148), and _sim_noise_null_sweep.qmd:63 (not included). There are no LaTeX-escaped forms, and no binary file (pdf, docx, zip or rds) contains the name. Generators: called in 34 of 37 documents, at :318/:321 (fs, grf) or :316/:319 (dina) in the sim_* documents, and at :247/:249 (h075, h15), :276/:278 (h10_knoise0 at HEAD; :250/:252 at run commit 4b6690b3) or :278/:280 (standalone2) in the sweeps. It executes in 29: the 20 batch documents (guard run_mode=="batch", :317/:315) and the 9 sweeps (guard if (run_sweep), run_sweep <- TRUE). In the 5 combine documents run_mode <- "combine" (:113), so it does not run there. Not used by the 2 application generators or the ACTG175 sweep. Outside manuscript_rev1, mention lines: analysis/ 48, manuscript/ 11, manuscript_legacy/ 9, manuscript_v40/ 8. Of these, call lines are analysis/ 18, and 6 in each manuscript* directory.
    - count: manuscript_rev1: 7 files / 11 lines (6 executable calls, all in 3 standalone dgm_truth_targets*.qmd; none in the render). Generators: 34 of 37 documents contain it; 29 execute it.
    - example: forestsearch quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_m1_h10_knoise0_n500_batch_1_50.qmd:321  dgm <- setup_gbsg_dgm(model = dgm_model, k_inter = k_inter, ... ; manuscript_rev1/dgm_truth_targets.qmd:156  dgm <- setup_gbsg_dgm(model = model, ...
  - **calibrate_k_inter (R/sim_aft_gbsg.R:1004)**
    - facts: Exported (NAMESPACE:92). Its nested eff_at (:1023) calls .create_gbsg_dgm_ at :1024, and .calibrate_by_root (:1040; effect_fun = eff_at at :1041) finds the root. In manuscript_rev1 it appears in 4 files on 16 lines. In each of dgm_truth_targets.qmd, _h10.qmd and _h15.qmd it is at :108 (prose), :139 (comment), and :150 and :227 (calls). dgm_truth_targets_h10.html has it at :2321, :2777, :2788 and :3742. No main or supplement source, render fragment, .tex or binary mentions it. Generators: one call in each of the same 34 documents (:318, :316, :247, :276 or :278), executed in 29. Outside manuscript_rev1, mention lines: analysis/ 48, manuscript/ 16, manuscript_legacy/ 16, manuscript_v40/ 16. Of these, call lines are analysis/ 18, and 6 in each manuscript* directory.
    - count: manuscript_rev1: 4 files / 16 lines (6 executable calls in 3 standalone dgm_truth_targets*.qmd; none in the render). Generators: 34 of 37 documents contain it; 29 execute it.
    - example: forestsearch quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_m1_h10_knoise0_n500_batch_1_50.qmd:318  k_inter <- calibrate_k_inter(target_hr_harm = target_hr_harm, model = dgm_model, use_ahr = FALSE) ; manuscript_rev1/dgm_truth_targets.qmd:150
  - **create_gbsg_dgm (R/sim_aft_gbsg.R:199)**
    - facts: Exported (NAMESPACE:116) and deprecated: .Deprecated() at :211, with the message 'is deprecated' at :215. It delegates to .create_gbsg_dgm_ at :226. manuscript_rev1 has no calls; the name appears in prose only, in the code-availability sentence. Plain-text files: 11 (one line each). These are the current supplement (postselection_subgroup_inference_supplement.qmd:1184, under '# Computational details and reproducibility' at :1182), _survey_2026-08/report_C_fb_mr_inventory.md:833, and 9 archived qmd files (jasa-2026 main_blinded:821, main_v5_1:1312, main_v5_6:1344, main_v5_7:1334, main_v5_7_blinded:1323, main_v5_8:832, main_v5_8_SUBMITTED:793, main_v5_8_blinded:783; pre-trim main:884). LaTeX-escaped .tex files: 6. These are the current supplement's postselection_subgroup_inference_supplement.tex:3296, and archived main_blinded.tex:1519, main_v5_6.tex:2156, main_v5_8.tex:1531, main_v5_8_blinded.tex:1487 and pre-trim main.tex:1587. Binary files (text extracted): 11. These are postselection_subgroup_inference_supplement.pdf (pdftotext line 4914), _archive/jasa-2026/ACC_form.pdf (2 lines, 66 and 168), LFL_acc_form.docx, LFL_acc_form_filled.docx, archive_submitted_v5_8.zip (member main_v5_8_SUBMITTED.qmd:793), and 6 archived main PDFs (blinded, pre-phase0, v5_6, v5_8, v5_8_blinded, pre-trim). Generators: 0 of 37. Outside manuscript_rev1: plain text manuscript/ 11 and manuscript_legacy/ 6 lines, plus escaped .tex manuscript/ 6 and manuscript_legacy/ 6; 0 calls.
    - count: manuscript_rev1: 17 text files / 17 lines (11 plain plus 6 LaTeX-escaped .tex) and 11 binary files; all prose, 0 calls. Generators: none.
    - example: manuscript_rev1/postselection_subgroup_inference_supplement.qmd:1184  "The simulation study's data-generating mechanisms are `create_gbsg_dgm()` and `generate_aft_dgm_flex()`."
  - **.create_gbsg_dgm_ (R/sim_aft_gbsg.R:236)**
    - facts: Internal. The comment at :233-234 says it is used by calibrate_k_inter, get_dgm_with_output, scan_k_inter_effect and setup_gbsg_dgm; scan_k_inter_effect is not defined anywhere in R/. Its actual callers are :226 (create_gbsg_dgm), :1024 (eff_at), :1108 (get_dgm_with_output), :1254 (validate_k_inter_effect) and setup_gbsg_dgm.R:106. No tracked file in fs-post-selection names it, text or binary, and no generator does either; there is no forestsearch::: access. Every executed generator call reaches it through setup_gbsg_dgm.R:106 and through eff_at at sim_aft_gbsg.R:1024.
    - count: none (direct); reached indirectly in the 29 generators that execute a wrapper call
    - example: none in fs-post-selection. Reached at R/setup_gbsg_dgm.R:106  dgm <- .create_gbsg_dgm_(model = model, ...
  - **simulate_from_gbsg_dgm (R/sim_aft_gbsg.R:789)**
    - facts: Exported (NAMESPACE:293) and deprecated: .Deprecated() at :800, message at :804. Zero occurrences in all tracked text of fs-post-selection, in manuscript_rev1 binaries and rds, and in all 37 generators. The generators simulate with simulate_from_dgm instead, which is not a wrapper function (e.g. sim_fs_maxeffCons_fb_mr_m1_h10_knoise0_n500_batch_1_50.qmd:468 and maxeffCons_mr_coverage_sweep_h10_knoise0_standalone2.qmd:428). In package R/ it is mentioned at run_simulation_analysis.R:6 (comment) and :198 (roxygen).
    - count: none
    - example: none
  - **.simulate_from_gbsg_dgm_ (R/sim_aft_gbsg.R:823)**
    - facts: Internal; its only caller is simulate_from_gbsg_dgm (:813). It is not named anywhere else in package R/ (the run_simulation_analysis.R:6 comment names the exported function, not this one). It is also not named in fs-post-selection or any generator, and no generator reaches it.
    - count: none
    - example: none
  - **calibrate_k_inter nested helper eff_at (R/sim_aft_gbsg.R:1023)**
    - facts: Local closure inside calibrate_k_inter; it calls .create_gbsg_dgm_ at :1024 and is passed as effect_fun at :1041. Closures with the same name in calibrate_k_treat.R:133, calibrate_glm_interaction.R:148 and find_k_inter_main.R:105 are different functions. Not named in fs-post-selection or the generators. It runs whenever calibrate_k_inter runs.
    - count: none (direct); runs in the 29 generators that execute a wrapper call
    - example: none
  - **get_dgm_with_output (R/sim_aft_gbsg.R:1076)**
    - facts: @keywords internal (:1075) and not in NAMESPACE (man/get_dgm_with_output.Rd exists). It calls calibrate_k_inter (:1093) and .create_gbsg_dgm_ (:1108). Zero occurrences in fs-post-selection (text, binary, rds) and in all generators.
    - count: none
    - example: none
  - **print.gbsg_dgm (R/sim_aft_gbsg.R:1157)**
    - facts: Registered as S3method(print,gbsg_dgm) at NAMESPACE:48. NAMESPACE has no print or summary method for aft_dgm_flex, so printing a setup_gbsg_dgm object (class set at setup_gbsg_dgm.R:182) would dispatch here. No generator or manuscript_rev1 document has an explicit print, summary or bare auto-print of a dgm object. The string gbsg_dgm (as a standalone word) occurs nowhere in fs-post-selection tracked text or rds.
    - count: none
    - example: none
  - **validate_k_inter_effect (R/sim_aft_gbsg.R:1229)**
    - facts: Exported (NAMESPACE:312); it calls .create_gbsg_dgm_ at :1254. Zero occurrences in fs-post-selection and in all generators.
    - count: none
    - example: none
  - **cut_numeric (R/sim_aft_gbsg.R:51)**
    - facts: @keywords internal (:50), not exported (man/cut_numeric.Rd exists). Not named in fs-post-selection or the generators. Its only callers are in .create_gbsg_dgm_: :281 (z2 from age), :283 (z4 from pgr) and :284 (z5 from nodes). It therefore runs in every generator that executes setup_gbsg_dgm or calibrate_k_inter.
    - count: none (direct); runs in 29 generators
    - example: none in fs-post-selection. Package: R/sim_aft_gbsg.R:281  dfa$z2 <- cut_numeric(dfa$age)
  - **cut_size (R/sim_aft_gbsg.R:69)**
    - facts: @keywords internal (:68), not exported. Not named in fs-post-selection or the generators. Its only caller is .create_gbsg_dgm_ at :298 (v6 from size).
    - count: none (direct); runs in 29 generators
    - example: none in fs-post-selection. Package: R/sim_aft_gbsg.R:298  dfa$v6 <- as.factor(cut_size(dfa$size))
  - **setup_gbsg_dgm nested helper .rename (R/setup_gbsg_dgm.R:149)**
    - facts: Local closure, used at setup_gbsg_dgm.R:155-159. No .rename( call appears in manuscript_rev1 or the generators. It runs on every setup_gbsg_dgm call.
    - count: none (direct); runs in 29 generators
    - example: none
- **Files named:**
  - [wrapper source, read-only export c40d5309] /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/R/sim_aft_gbsg.R
  - [wrapper source] /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/R/setup_gbsg_dgm.R
  - [package NAMESPACE] /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/NAMESPACE
  - [package, other caller, context] /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/R/fs_fdr_report.R
  - [package, mentions simulate_from_gbsg_dgm :6, :198] /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/R/run_simulation_analysis.R
  - [package, called by the supplement render, no wrapper] /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/R/fs_identification_summary.R
  - [current render: supplement source, prose create_gbsg_dgm :1184] /home/larryleon/Documents/GitHub/fs-post-selection/manuscript_rev1/postselection_subgroup_inference_supplement.qmd
  - [current kept .tex of the supplement, escaped create_gbsg_dgm :3296 (missed by investigator)] /home/larryleon/Documents/GitHub/fs-post-selection/manuscript_rev1/postselection_subgroup_inference_supplement.tex
  - [current rendered binary: prose create_gbsg_dgm] /home/larryleon/Documents/GitHub/fs-post-selection/manuscript_rev1/postselection_subgroup_inference_supplement.pdf
  - [current render: main source, no mentions; reads payloads :725-729] /home/larryleon/Documents/GitHub/fs-post-selection/manuscript_rev1/postselection_subgroup_inference_main.qmd
  - [current render: derivations, no mentions] /home/larryleon/Documents/GitHub/fs-post-selection/manuscript_rev1/postselection_theory_derivations.qmd
  - [project config: type default, no render list] /home/larryleon/Documents/GitHub/fs-post-selection/manuscript_rev1/_quarto.yml
  - [current render fragment: comment setup_gbsg_dgm :124] /home/larryleon/Documents/GitHub/fs-post-selection/manuscript_rev1/_sim_estimation_fs_h15_n1000.qmd
  - [current render fragment: comment setup_gbsg_dgm :63] /home/larryleon/Documents/GitHub/fs-post-selection/manuscript_rev1/_sim_noise_null_sweep_standalone2.qmd
  - [standalone, not included: comment setup_gbsg_dgm :63] /home/larryleon/Documents/GitHub/fs-post-selection/manuscript_rev1/_sim_noise_null_sweep.qmd
  - [current render fragment, library(forestsearch) :19, no wrapper] /home/larryleon/Documents/GitHub/fs-post-selection/manuscript_rev1/fs_identification_figures_h10_standalone2.qmd
  - [current render fragment, library(forestsearch) :19, no wrapper] /home/larryleon/Documents/GitHub/fs-post-selection/manuscript_rev1/fs_identification_figures_h15_standalone2.qmd
  - [auxiliary standalone, not rendered: calls :150/:156/:227/:230] /home/larryleon/Documents/GitHub/fs-post-selection/manuscript_rev1/dgm_truth_targets.qmd
  - [auxiliary standalone: calls] /home/larryleon/Documents/GitHub/fs-post-selection/manuscript_rev1/dgm_truth_targets_h10.qmd
  - [auxiliary standalone: calls] /home/larryleon/Documents/GitHub/fs-post-selection/manuscript_rev1/dgm_truth_targets_h15.qmd
  - [auxiliary rendered HTML: code lines] /home/larryleon/Documents/GitHub/fs-post-selection/manuscript_rev1/dgm_truth_targets_h10.html
  - [survey record: quote of create_gbsg_dgm prose :833] /home/larryleon/Documents/GitHub/fs-post-selection/manuscript_rev1/_survey_2026-08/report_C_fb_mr_inventory.md
  - [archive qmd x9: prose create_gbsg_dgm] /home/larryleon/Documents/GitHub/fs-post-selection/manuscript_rev1/_archive/jasa-2026/postselection_subgroup_inference_main_{blinded,v5_1,v5_6,v5_7,v5_7_blinded,v5_8,v5_8_SUBMITTED,v5_8_blinded}.qmd and /home/larryleon/Documents/GitHub/fs-post-selection/manuscript_rev1/_archive/pre-trim/postselection_subgroup_inference_main.qmd
  - [archive tex, escaped create_gbsg_dgm (missed by investigator)] /home/larryleon/Documents/GitHub/fs-post-selection/manuscript_rev1/_archive/jasa-2026/postselection_subgroup_inference_main_blinded.tex
  - [archive tex] /home/larryleon/Documents/GitHub/fs-post-selection/manuscript_rev1/_archive/jasa-2026/postselection_subgroup_inference_main_v5_6.tex
  - [archive tex] /home/larryleon/Documents/GitHub/fs-post-selection/manuscript_rev1/_archive/jasa-2026/postselection_subgroup_inference_main_v5_8.tex
  - [archive tex] /home/larryleon/Documents/GitHub/fs-post-selection/manuscript_rev1/_archive/jasa-2026/postselection_subgroup_inference_main_v5_8_blinded.tex
  - [archive tex] /home/larryleon/Documents/GitHub/fs-post-selection/manuscript_rev1/_archive/pre-trim/postselection_subgroup_inference_main.tex
  - [archive binaries: create_gbsg_dgm] /home/larryleon/Documents/GitHub/fs-post-selection/manuscript_rev1/_archive/jasa-2026/{ACC_form.pdf, LFL_acc_form.docx, LFL_acc_form_filled.docx, archive_submitted_v5_8.zip, postselection_subgroup_inference_main_blinded.pdf, _main_pre-phase0.pdf, _main_v5_6.pdf, _main_v5_8.pdf, _main_v5_8_blinded.pdf} and /home/larryleon/Documents/GitHub/fs-post-selection/manuscript_rev1/_archive/pre-trim/postselection_subgroup_inference_main.pdf
  - [provenance, stale 2026-08-17] /home/larryleon/Documents/GitHub/fs-post-selection/manuscript_rev1/payload_manifest.qmd
  - [provenance] /home/larryleon/Documents/GitHub/fs-post-selection/manuscript_rev1/REVISION_LOG.md
  - [provenance, render manifest, older; aux list :666/:669-671] /home/larryleon/Documents/GitHub/fs-post-selection/manuscript_rev1/manuscript_render_manifest.html
  - [provenance] /home/larryleon/Documents/GitHub/fs-post-selection/dev/reports/rev1_setup_inventory_2026-11-10.csv
  - [provenance] /home/larryleon/Documents/GitHub/fs-post-selection/dev/reports/REPORT_rev1_setup_inventory_2026-11-10.md
  - [provenance] /home/larryleon/Documents/GitHub/fs-post-selection/dev/reports/REPORT_rev1_reruns_2026-11-10.md
  - [provenance :11-14] /home/larryleon/Documents/GitHub/fs-post-selection/dev/reports/REPORT_a002_table_s4_2026-11-10.md
  - [provenance] /home/larryleon/Documents/GitHub/fs-post-selection/dev/reports/REPORT_rev1_mac_timing_2026-11-11.md
  - [provenance: render commands :97-99] /home/larryleon/Documents/GitHub/fs-post-selection/dev/reports/REPORT_rev1_averaging_text_2026-11-11.md
  - [current render fragments reading payloads, no wrapper mentions] /home/larryleon/Documents/GitHub/fs-post-selection/manuscript_rev1/{_gbsg_section6_1, _sim_synthesis_h10, _supp_intro_sim_summary, _supp_gbsg_table, _actg175_supp, _sim_additional_results_intro, _sim_estimation_fs_h10, _sim_estimation_dina_h10, _sim_estimation_grf_h10, _supp_sim_synthesis, _sim_estimation_fs_h10_n1000, _sim_noise_pair_n1500_standalone2, _sim_mr_coverage_{h075,h10,h15}_{conditional,marginal}_target, _sim_mr_coverage_or075_actg175_{conditional,marginal}_target, _sim_mr_coverage_h10_knoise_standalone2}.qmd
  - [generators, sim fs h10 n500: calls :318/:321; combiner guarded, not executed] /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_m1_h10_knoise0_n500_{batch_1_50,batch_51_100,batch_101_200,batch_201_300,batch_301_500,combine_1_500}.qmd
  - [generators, sim fs h10 n1000: calls :318/:321] /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_m1_h10_knoise0_n1000_{batch_1_20,batch_21_100,batch_101_110,batch_111_500,combine_1_500}.qmd
  - [generators, sim fs h15 n1000: calls :318/:321] /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_m1_h15_knoise0_n1000_{batch_1_200,batch_201_220,batch_221_300,batch_301_320,batch_321_400,combine_1_400}.qmd
  - [generators, sim dina: calls :316/:319] /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020/sim_dina_eff_fb_mr_m1_h10_knoise0_n500_{batch_1_200,batch_201_250,batch_251_500,combine_1_500}.qmd
  - [generators, sim grf: calls :318/:321] /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020/sim_grf_eff_fb_mr_m1_h10_knoise0_n500_{batch_1_200,batch_201_250,batch_251_500,combine_1_500}.qmd
  - [generator sweeps, calls :247/:249, guard :246] /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020/maxeffCons_mr_coverage_sweep_h075_knoise0.qmd and maxeffCons_mr_coverage_sweep_h15_knoise0.qmd
  - [generator sweep, calls :276/:278 at HEAD (campaign_tag standalone1 :158); :250/:252 at run commit 4b6690b3] /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020/maxeffCons_mr_coverage_sweep_h10_knoise0.qmd
  - [generator sweeps x6, calls :278/:280, guard :277, campaign_tag standalone2 :160] /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_020/maxeffCons_mr_coverage_sweep_{h10_knoise0,h10_knoise3,h10_knoise6,h15_knoise0,h15_knoise3_n1500,h20_knoise3_n1500}_standalone2.qmd
  - [generator sweep, no wrapper; GLM DGM :249/:279] /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/actg175/binary_020/maxeffCons_mr_coverage_sweep_or075.qmd
  - [generator application, no wrapper (run blob fded7718, HEAD 114aa11d)] /home/larryleon/Documents/GitHub/forestsearch/quarto/applications/gbsg/analysis_gbsg_survival_multimethod.qmd
  - [generator application, no wrapper (fc59934b = HEAD)] /home/larryleon/Documents/GitHub/forestsearch/quarto/applications/actg175/analysis_actg175_binary_multimethod_fixed_family.qmd
  - [stale-manifest generator example, calls :249/:252 and :699/:701] /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_redux/fs_t1_t2_m1_h10_knoise0_n500_combine_1_500.qmd
  - [stale-manifest generator example, calls :170/:172] /home/larryleon/Documents/GitHub/forestsearch/quarto/simulations/gbsg_redux/mr_coverage_sweep_h10.qmd
  - [verifier scratch: extracted binary text, rds copies, gbsg_redux lists, outside-rev1 call lines] /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2ro/L-E-verify/
- **Notes:** Bottom line, confirmed. Current manuscript work (manuscript_rev1 and the generators of its payloads) depends on two wrapper functions directly: setup_gbsg_dgm() and calibrate_k_inter(). Through them it also runs .create_gbsg_dgm_, cut_numeric, cut_size, eff_at and .rename. Nothing in that work calls create_gbsg_dgm, simulate_from_gbsg_dgm, .simulate_from_gbsg_dgm_, get_dgm_with_output, validate_k_inter_effect or print.gbsg_dgm. The one text dependency is the code-availability sentence naming create_gbsg_dgm(). It appears in the current supplement source (postselection_subgroup_inference_supplement.qmd:1184), its kept .tex (:3296) and its PDF (pdftotext line 4914). Removing or renaming create_gbsg_dgm would leave the published sentence pointing at a function that no longer exists.

Where the dependence sits:
- **Render (main + supplement + derivations).** No wrapper call. The main paper and the supplement read .rds from _payloads/. The supplement also runs library(forestsearch) and fs_identification_* functions (fs_identification_figures_h10/h15_standalone2.qmd:19-...), which reach no wrapper function.
- **Executable calls under manuscript_rev1.** 12 lines (6 calibrate_k_inter, 6 setup_gbsg_dgm), only in dgm_truth_targets{,_h10,_h15}.qmd at :150/:156/:227/:230. Neither the main paper nor the supplement includes them. manuscript_render_manifest.html:669-671 lists two of the three (not _h10) as auxiliary. _quarto.yml has no render list, so a whole-project quarto render would execute them.

Generators of the current render's payloads (paths relative to forestsearch/quarto/). This is the investigator's list, all re-verified:
1. **analysis_gbsg_survival_multimethod_payload.rds** (_gbsg_section6_1.qmd:21, _supp_gbsg_table.qmd:25). Generator applications/gbsg/analysis_gbsg_survival_multimethod.qmd. Blobs fded7718 (run) and 114aa11d (HEAD); no wrapper.
2. **analysis_actg175_binary_multimethod_fixed_family_payload.rds** (_actg175_supp.qmd:19). Generator applications/actg175/analysis_actg175_binary_multimethod_fixed_family.qmd. No wrapper in any version.
3. **fs_maxeffCons_fb_mr_m1_h10_knoise0_n500_combined_1_500.rds** (main.qmd:727 and others). Generator simulations/gbsg_020/sim_fs_..._n500_batch_{1_50,51_100,101_200,201_300,301_500} then _combine_1_500. The revision-1 content was re-run from scratch copies of the batch_1_50 template and the combiner (REPORT_a002_table_s4:11-14; REVISION_LOG.md:70, :616). Calls at :318/:321.
4. **dina_eff_fb_mr_m1_h10_knoise0_n500_combined_1_500.rds** (main.qmd:728). Generator sim_dina_eff_... batch_{1_200,201_250,251_500} then _combine_1_500. meta$source_files matches. Calls at :316/:319.
5. **grf_eff_..._combined_1_500.rds** (main.qmd:729). Generator sim_grf_eff_... batches then combine. Calls at :318/:321.
6. **fs_maxeffCons_..._h10_knoise0_n1000_combined_1_500.rds** (_sim_estimation_fs_h10_n1000.qmd:28). Re-run from blob f0121563 (batch_1_20 template) and e9ba13a8 (combine), both equal to HEAD. Calls at :318/:321.
7. **fs_maxeffCons_..._h15_knoise0_n1000_combined_1_400.rds**. Re-run from blobs 74b0e8eb and 8dc86a4b, equal to HEAD. Calls at :318/:321.
8. **mr_coverage_grid_m1_{h08,h10,h15}_knoise0_maxeffCons_mr_seedtab_s500.rds** (built 2026-08-18). Generators maxeffCons_mr_coverage_sweep_h075_knoise0.qmd (scenario_tag gives h08 from round(10*0.75)) and _h15_knoise0.qmd, each a single commit with calls at :247/:249. The h10 grid came from sweep_h10_knoise0.qmd at 4b6690b3 (calls :250/:252). HEAD (:276/:278) now writes a _standalone1 tag (campaign_tag :158).
9. **The six *_standalone2 folders.** Generators maxeffCons_mr_coverage_sweep_*_standalone2.qmd, with calls at :278/:280 in every historical version.
10. **mr_coverage_grid_maxeffCons_actg175_or075_seedtab_s1000.rds.** Generator simulations/actg175/binary_020/maxeffCons_mr_coverage_sweep_or075.qmd (run_tag :98). No wrapper in any of its 3 commits; it uses calibrate_glm_interaction (:249) and generate_glm_dgm (:279).

Totals: 37 generator documents. 34 contain calibrate_k_inter plus setup_gbsg_dgm. 29 execute them: 20 batch documents guarded by run_mode=="batch" (:317/:315), and 9 sweeps guarded by if (run_sweep) with run_sweep <- TRUE. The 5 combiners have run_mode "combine" at :113. The other 3 documents have no wrapper calls.

payload_manifest.qmd (generator entries dated 2026-08-17 at :26; active tribble :78-138, inactive :172-197) and manuscript_render_manifest.html are stale: no maxeffCons, eff or standalone2 names (see REPORT_rev1_setup_inventory_2026-11-10.md:61). The 65 gbsg_redux documents they name all call calibrate_k_inter and setup_gbsg_dgm, e.g. gbsg_redux/fs_t1_t2_m1_h10_knoise0_n500_combine_1_500.qmd:249/:252 and :699/:701. The 2 actg175/binary sweeps call none.

**fs-post-selection outside manuscript_rev1** (tracked text; mention lines, then code call lines):

| Directory | Mention lines (files) | calibrate_k_inter | setup_gbsg_dgm | create_gbsg_dgm | Code call lines (files) |
|---|---|---|---|---|---|
| analysis/ | 96 (33) | 48 | 48 | 0 | 36 (18): 18 calibrate_k_inter + 18 setup_gbsg_dgm |
| manuscript/ | 38 (18) | 16 | 11 | 11, plus 6 escaped .tex | 12 (3, dgm_truth_targets*.qmd) |
| manuscript_legacy/ | 31 (11) | 16 | 9 | 6, plus 6 escaped .tex | 12 (3) |
| manuscript_v40/ | 24 (4) | 16 | 8 | 0 | 12 (3) |

No other wrapper function appears in tracked text anywhere in fs-post-selection, nor in any manuscript_rev1 binary or rds. All repositories were left untouched (read-only git and file reads only). Scratch output is in /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2ro/L-E-verify/.
- **Verifier (corrected) corrections:**
  - create_gbsg_dgm count is incomplete. The search missed 6 tracked .tex files that write the name LaTeX-escaped as create\_gbsg\_dgm(), which a raw-identifier grep cannot match. One is the current supplement's kept .tex, manuscript_rev1/postselection_subgroup_inference_supplement.tex:3296. The other five are archives: _archive/jasa-2026/postselection_subgroup_inference_main_blinded.tex:1519, main_v5_6.tex:2156, main_v5_8.tex:1531, main_v5_8_blinded.tex:1487, and _archive/pre-trim/postselection_subgroup_inference_main.tex:1587. Corrected manuscript_rev1 total: 17 text files on 17 lines plus 11 binaries (28 files), all prose, 0 calls. Outside manuscript_rev1 the same escaped form adds 6 .tex files in manuscript/ and 6 in manuscript_legacy/ (prose only).
  - Wrong citation: REVISION_LOG.md:852 says nothing about the fs h10 n500 re-run. It describes the eq7 payload, which exists only in manuscript/. The re-run from scratch copies of batch_1_50 and combine_1_500 is supported by dev/reports/REPORT_a002_table_s4_2026-11-10.md:11-14 (template :11, combiner :12, scratch copies :14). In REVISION_LOG.md the matching entries are :70 (R0-002, a002 Table S4) and :616 (R1-031, the payload file).
  - The render-manifest citation is off and over-broad. manuscript_render_manifest.html has the 'D. NOT required to render' heading at :666 and the auxiliary list at :669-671. That list names only dgm_truth_targets.qmd and dgm_truth_targets_h15.qmd; dgm_truth_targets_h10.qmd appears nowhere in the manifest (0 hits). That the _h10 file is not part of the render follows from the include list, not from the manifest.
  - The claim 'the render only reads .rds files' needs refining. The supplement render runs library(forestsearch) in the included fs_identification_figures_h10_standalone2.qmd:19 and fs_identification_figures_h15_standalone2.qmd:19. It then calls fs_identification_summary, plot_fs_identification_*, and fs_identification_table (R/fs_identification_summary.R:222, :310, :347, :389). None of these reaches a wrapper function: that file contains 0 wrapper names, and the only package R/ code outside the two wrapper files that calls a wrapper is fs_fdr_report.R:256-262. So the render depends on forestsearch but on no wrapper function. The render set also includes postselection_theory_derivations.qmd, rendered to HTML (dev/reports/REPORT_rev1_averaging_text_2026-11-11.md:97-99), and it has no wrapper mentions.
  - Caveat the investigator missed: manuscript_rev1/_quarto.yml:1-2 declares project type default with no render list. A whole-project 'quarto render' in manuscript_rev1 would therefore render dgm_truth_targets*.qmd and execute their 12 wrapper call lines. The manuscript is rendered file by file (main, supplement, derivations), so this does not happen in the documented workflow.
  - Wrong attribution for .simulate_from_gbsg_dgm_: R/run_simulation_analysis.R:6 mentions the exported simulate_from_gbsg_dgm(), not .simulate_from_gbsg_dgm_. The internal function is named only at sim_aft_gbsg.R:813 (its call) and :823 (its definition). simulate_from_gbsg_dgm itself is mentioned in R/run_simulation_analysis.R at :6 (comment) and :198 (roxygen).
  - Line-number precision: the .Deprecated() calls are at sim_aft_gbsg.R:211 (create_gbsg_dgm) and :800 (simulate_from_gbsg_dgm). Lines :215 and :804 are the 'is deprecated' message strings inside those calls.
  - Missing execution guard for the sweeps: their wrapper calls sit inside if (run_sweep) at :246 (h075, h15), :275 (h10_knoise0 at HEAD; :249 at run commit 4b6690b3) and :277 (all six standalone2). run_sweep <- TRUE at :158, :168 and :170 (and :95 at 4b6690b3), so the conclusion that they execute stands. In the batch and combine documents the guard is if (identical(run_mode, "batch")) at :317 (fs and grf) and :315 (dina).
  - Incomplete example: gbsg_redux/fs_t1_t2_m1_h10_knoise0_n500_combine_1_500.qmd has a second pair of wrapper calls at :699 and :701, besides :249 and :252. The count of 65 gbsg_redux documents named by payload_manifest that call the wrapper reproduces.
  - The 'outside manuscript_rev1' table counts mention lines, not calls (line and occurrence counts are equal here). Counting only non-comment code call lines in .qmd/.R/.Rmd files: analysis/ has 36 lines in 18 files (18 calibrate_k_inter, 18 setup_gbsg_dgm). manuscript/, manuscript_legacy/ and manuscript_v40/ each have 12 lines in 3 files, the copies of dgm_truth_targets{,_h10,_h15}.qmd (6 plus 6). create_gbsg_dgm has 0 calls outside manuscript_rev1.
  - Minor: git classes 208 tracked manuscript_rev1 files as binary (172 rds, 21 pdf, 10 docx, 2 zip, 2 .DS_Store, 1 tex), not 182, and it treats 3 _archive/process PDFs as text. This changes no result: I extracted text from all 24 PDFs, 10 docx and 2 zips and decompressed all 172 rds (all gzip). The only new hits are the ones already reported (11 binaries, create_gbsg_dgm only).

### L-F

- **Method:** Independent routes; all read-only. (1) Wrapper list: grep for every 'function(' / '<- function' in the two c40d5309 files gives 10 top-level definitions in sim_aft_gbsg.R plus setup_gbsg_dgm, and the nested eff_at (:1023) and .rename (:149). Cross-checked against the roxygen tags, t2_tree NAMESPACE and the man/ Rd files. (2) In fs-glms-interpretable (main, dc94a17; ls-files gives quarto 58, manuscript 61, dev/payloads 457 = 576 tracked files), git grep -c over the three directories, with and without -I, case-insensitive, for partial names: gbsg_dgm, create_gbsg, simulate_from, k_inter, calibrate_k, validate_k, cut_numeric, cut_size, dgm_with, setup_gbsg, sim_aft, eff_at, \.rename. All are 0 except k_inter (1, SPEC:102). Also searched every 'dgm' mention, forestsearch:: / ::: usages (only plot_pareto_frontier, fs_mr_inference, .fs_decl_settable_table), library/source/include targets, and string-based do.call/get/match.fun. (3) Binaries: confirmed that all 378 .rds files have a gzip header and grepped zcat output (sanity hit confirmed); ran pdftotext on the 10 PDFs and unzip -p on all docx parts. (4) Provenance: traced which payload CSVs the manuscript reads (identifiers_metrics, idsweep_metrics, declcalc0_summary_tables, glm/*, applications/*, theory/*, guohe rds), then checked the producing forestsearch scripts with read-only git show at the pinned commits (ca9b0a62, 69ca797e, dbc9444c, 1507bfc0, b15c8b7a, a46bf9b7, 59b144ea) and in t2_tree. (5) Repo-wide git grep for wrapper names, to complete the out-of-scope context list (20 files). No files written, no R run.
- **Per function:**
  - **create_gbsg_dgm**
    - facts: Defined at R/sim_aft_gbsg.R:199, @export at :198 (NAMESPACE:116). fs-glms-interpretable is on main at dc94a17. Over the 576 tracked files in quarto/ (58), manuscript/ (61) and dev/payloads/ (457), git grep finds no match for 'create_gbsg' or 'gbsg_dgm', with or without -I, case-insensitive. Decompressed .rds and extracted PDF/docx text also give 0 hits. All 11 tracked hits in the repo are text mentions outside scope in dev/reference/post_selection/manuscript_jrssb_initialsubmit/: archive drafts _archive/jasa-2026/*_v5_1.qmd:1312, v5_6:1344, v5_7:1334, v5_7_blinded:1323, main_blinded:821, v5_8:832, v5_8_SUBMITTED:793, v5_8_blinded:783, and _archive/pre-trim/..._main.qmd:884; reference copies postselection_subgroup_inference_supplement.qmd:1182 and _survey_2026-08/report_C_fb_mr_inventory.md:833.
    - count: 0 calls, 0 text mentions (quarto/ 0, manuscript/ 0, dev/payloads/ 0)
    - example: none
  - **.create_gbsg_dgm_**
    - facts: Internal worker defined at R/sim_aft_gbsg.R:236 (@keywords internal :235). Its name does not appear in any tracked file of fs-glms-interpretable. Related output-field text, not a name mention: dev/payloads/identifiers/COLUMNS.md:39 and dev/payloads/R/identifiers_metrics.R:562 mention `dgm$hr_H_true` / `dgm$hr_Hc_true`, fields computed in this worker at sim_aft_gbsg.R:527-548.
    - count: 0 calls, 0 text mentions (quarto/ 0, manuscript/ 0, dev/payloads/ 0)
    - example: none
  - **simulate_from_gbsg_dgm**
    - facts: Defined at R/sim_aft_gbsg.R:789, @export at :788 (NAMESPACE:293). No match for 'simulate_from_gbsg' (or 'gbsg_dgm') in the three directories. The GLM simulation payloads in dev/payloads/glm use forestsearch's simulate_from_glm_dgm/generate_glm_dgm (ACTG175 templates), not this function.
    - count: 0 calls, 0 text mentions (quarto/ 0, manuscript/ 0, dev/payloads/ 0)
    - example: none
  - **.simulate_from_gbsg_dgm_**
    - facts: Internal worker defined at R/sim_aft_gbsg.R:823 (@keywords internal :822). No match in the three directories, and none anywhere in the repo.
    - count: 0 calls, 0 text mentions (quarto/ 0, manuscript/ 0, dev/payloads/ 0)
    - example: none
  - **calibrate_k_inter**
    - facts: Defined at R/sim_aft_gbsg.R:1004, @export at :1003 (NAMESPACE:92). No match for 'calibrate_k' in the three directories. A case-insensitive 'k_inter' search there has exactly 1 hit, a concept mention and not the function name: quarto/gbsg/SPEC_gbsg_grf_dina_2026-09-17.md:102 ('via `k_inter` calibration'), under '### 2.5 DGM — does not transfer' (:97). That SPEC is current: dev/manuscript_status_update.md:1045 says it specifies the GBSG trio, the 'manuscript's core analysis example'. INDIRECT, all through forestsearch-produced payloads: (a) identifiers bundles, from the gbsg_020 template per COLUMNS.md:32; the call is at template :703 (ca9b0a62/69ca797e), :740 (dbc9444c) and :816 (1507bfc0 = c40d5309); (b) idsweep bundles, from the same template at b15c8b7a (:740); (c) declcal payloads, from scripts_dinamr/declcalc0_run.R:148 and declcal_run.R:132. Out of scope: real calls in dev/reference/post_selection/manuscript_jrssb_initialsubmit/dgm_truth_targets{,_h10,_h15}.qmd:150 and :227, with text mentions at :108 and :139; the rendered _h10.html carries the same text.
    - count: 0 calls (quarto/ 0, manuscript/ 0, dev/payloads/ 0); 0 function-name mentions; 1 concept mention ('k_inter calibration') in quarto/gbsg
    - example: none (concept mention only: quarto/gbsg/SPEC_gbsg_grf_dina_2026-09-17.md:102)
  - **eff_at (nested in calibrate_k_inter)**
    - facts: Nested helper at R/sim_aft_gbsg.R:1023, inside calibrate_k_inter. Same-named nested helpers exist elsewhere in forestsearch (calibrate_glm_interaction.R:148, calibrate_k_treat.R:133, find_k_inter_main.R:105), but they are not part of the wrapper. No 'eff_at' match in the three directories.
    - count: 0 calls, 0 text mentions (quarto/ 0, manuscript/ 0, dev/payloads/ 0)
    - example: none
  - **get_dgm_with_output**
    - facts: Defined at R/sim_aft_gbsg.R:1076. It is @keywords internal at :1075 (not :1074) and absent from NAMESPACE; man/get_dgm_with_output.Rd exists. No match for 'dgm_with' in the three directories, and none anywhere in the repo.
    - count: 0 calls, 0 text mentions (quarto/ 0, manuscript/ 0, dev/payloads/ 0)
    - example: none
  - **print.gbsg_dgm**
    - facts: S3 method at R/sim_aft_gbsg.R:1157, @export at :1156, NAMESPACE:48 S3method(print,gbsg_dgm). No 'gbsg_dgm' string in the tracked text of the three directories. All 378 tracked .rds files there are gzip (header 1f8b08): 63 bundles, 18 bundles/superseded, 288 idsweep, 1 manuscript/_payloads, 3 quarto/actg175, 5 quarto/gbsg. Decompressing them and grepping for gbsg_dgm, k_inter, hr_H_true, df_super and the wrapper names gives 0 hits; a sanity grep for 'nv_H_est|or_H_est|sens' on a bundle did hit. The 10 PDFs (pdftotext) and the 1 docx (all parts) also give 0 hits. So nothing could dispatch to this method.
    - count: 0 calls (explicit or by dispatch), 0 text mentions (quarto/ 0, manuscript/ 0, dev/payloads/ 0)
    - example: none
  - **validate_k_inter_effect**
    - facts: Defined at R/sim_aft_gbsg.R:1229, @export at :1228 (NAMESPACE:312). No match for 'validate_k' in the three directories, and none anywhere in the repo.
    - count: 0 calls, 0 text mentions (quarto/ 0, manuscript/ 0, dev/payloads/ 0)
    - example: none
  - **cut_numeric**
    - facts: Internal helper at R/sim_aft_gbsg.R:51 (@keywords internal :50; man/cut_numeric.Rd). No match in the three directories, and none anywhere in the repo.
    - count: 0 calls, 0 text mentions (quarto/ 0, manuscript/ 0, dev/payloads/ 0)
    - example: none
  - **cut_size**
    - facts: Internal helper at R/sim_aft_gbsg.R:69 (@keywords internal :68; man/cut_size.Rd). No match in the three directories, and none anywhere in the repo.
    - count: 0 calls, 0 text mentions (quarto/ 0, manuscript/ 0, dev/payloads/ 0)
    - example: none
  - **setup_gbsg_dgm**
    - facts: Defined at R/setup_gbsg_dgm.R:88, @export at :86 (NAMESPACE:286). It calls .create_gbsg_dgm_ (setup_gbsg_dgm.R:106). No match for 'setup_gbsg' in the three directories. INDIRECT, through forestsearch payloads the manuscript reads: (a) identifiers_metrics.csv, read at manuscript/_05_operating_characteristics.qmd:14 and _S4_operating_characteristics_supplement.qmd:16; its bundles came from the gbsg_020 template (COLUMNS.md:32), where setup_gbsg_dgm is called at :707 (ca9b0a62/69ca797e), :744 (dbc9444c) and :821/:829 (1507bfc0 = c40d5309); (b) idsweep_metrics.csv, read at manuscript/_02_identification.qmd:4 and _S4:17, from the same template at b15c8b7a (:744; README.md:40-42; idsweep.sh:51); (c) declcalc0_summary_tables.csv, read at _09_field_constructions_supplement.qmd:824, from declcalc0_run.R:156 and :163 (a46bf9b7 and c40d5309). Output-field text: COLUMNS.md:39 and identifiers_metrics.R:562 mention `dgm$hr_H_true`. Out of scope: comment in the TODO stub sims/configs/opchar_survival_gbsg.R:2, quoted at dev/tasks/REPORT_glms_structure_review_2026-09-11.md:275; comments in dev/reference/.../_sim_estimation_fs_h15_n1000.qmd:124, _sim_noise_null_sweep.qmd:63 and _sim_noise_null_sweep_standalone2.qmd:63; real calls in dgm_truth_targets{,_h10,_h15}.qmd:156 and :230.
    - count: 0 calls, 0 text mentions (quarto/ 0, manuscript/ 0, dev/payloads/ 0)
    - example: none
  - **.rename (nested in setup_gbsg_dgm)**
    - facts: Nested helper at R/setup_gbsg_dgm.R:149, used at :155-159; for example it renames flag.harm to flag_harm at :159. No '\.rename' match in the three directories. The column it produces, flag_harm, is named in text at dev/payloads/identifiers/COLUMNS.md:34 and dev/payloads/R/identifiers_metrics.R:557 ('harm_col <- "flag_harm"'). That is a column mention, not a call.
    - count: 0 calls, 0 text mentions (quarto/ 0, manuscript/ 0, dev/payloads/ 0)
    - example: none
- **Files named:**
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/R/sim_aft_gbsg.R (forestsearch c40d5309 source: 10 top-level definitions + nested eff_at :1023)
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/R/setup_gbsg_dgm.R (forestsearch c40d5309 source: setup_gbsg_dgm :88 + nested .rename :149)
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/NAMESPACE (exports :48, :92, :116, :286, :293, :312)
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_field_m1_template.qmd (forestsearch, outside fs-glms-interpretable; upstream producer of identifiers and idsweep bundles; calls :816/:821/:829 at c40d5309 and 1507bfc0, :703/:707 at ca9b0a62 and 69ca797e, :740/:744 at dbc9444c and b15c8b7a)
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/quarto/simulations/gbsg_020/scripts_dinamr/declcalc0_run.R (forestsearch; producer of the declcalc0 payloads; calibrate_k_inter :148, setup_gbsg_dgm :156, :163)
  - /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2_tree/quarto/simulations/gbsg_020/scripts_dinamr/declcal_run.R (forestsearch; calibrate_k_inter :132, setup_gbsg_dgm :140, :147)
  - /home/larryleon/Documents/GitHub/fs-glms-interpretable/quarto/gbsg/SPEC_gbsg_grf_dina_2026-09-17.md (current work, not superseded; concept mention 'k_inter calibration' at :102 under §2.5 :97; no call)
  - /home/larryleon/Documents/GitHub/fs-glms-interpretable/dev/payloads/identifiers/COLUMNS.md (current manuscript payload doc; provenance at :32; output-field text `dgm$hr_H_true` at :39 and `flag_harm` at :34; no wrapper name)
  - /home/larryleon/Documents/GitHub/fs-glms-interpretable/dev/payloads/R/identifiers_metrics.R (current payload producer; writes identifiers_metrics.csv and COLUMNS.md (:2, :607); the same output-field text as string literals at :557 and :562; no call)
  - /home/larryleon/Documents/GitHub/fs-glms-interpretable/dev/payloads/identifiers/README.md (current; :7-8 the manuscript reads identifiers_metrics.csv; :30, :40-42 idsweep provenance from forestsearch b15c8b7a)
  - /home/larryleon/Documents/GitHub/fs-glms-interpretable/dev/payloads/R/idsweep_metrics.R (current payload producer of idsweep_metrics.csv; no call)
  - /home/larryleon/Documents/GitHub/fs-glms-interpretable/manuscript/_05_operating_characteristics.qmd (current manuscript; :14 reads identifiers_metrics.csv; no call)
  - /home/larryleon/Documents/GitHub/fs-glms-interpretable/manuscript/_02_identification.qmd (current manuscript; :4 reads idsweep_metrics.csv; no call)
  - /home/larryleon/Documents/GitHub/fs-glms-interpretable/manuscript/_S4_operating_characteristics_supplement.qmd (current supplement; :16 and :17 read both CSVs; no call)
  - /home/larryleon/Documents/GitHub/fs-glms-interpretable/manuscript/_09_field_constructions_supplement.qmd (current supplement; :823-824 reads declcalc0_summary_tables.csv; no call)
  - /home/larryleon/Documents/GitHub/fs-glms-interpretable/manuscript/_payloads/declcal/declcalc0_summary_tables.provenance.md (current; :3 producer and column semantics from forestsearch declcalc0_run.R at a46bf9b7)
  - /home/larryleon/Documents/GitHub/fs-glms-interpretable/manuscript/_R/oc_figures.R (current manuscript; :3 header comment naming identifiers_metrics.csv; no call)
  - /home/larryleon/Documents/GitHub/fs-glms-interpretable/manuscript/_06_applications.qmd (current manuscript; :4-15 reads dev/payloads/applications; no call)
  - /home/larryleon/Documents/GitHub/fs-glms-interpretable/dev/payloads/R/applications_metrics.R (current payload producer; :2-6 transforms only; :47-51 and :69 read the quarto/gbsg and quarto/actg175 payloads; no call)
  - /home/larryleon/Documents/GitHub/fs-glms-interpretable/quarto/gbsg/analysis_gbsg_mr.qmd (current; real survival::gbsg data, :3 and :62; library(forestsearch) at :150; no call)
  - /home/larryleon/Documents/GitHub/fs-glms-interpretable/dev/payloads/criteria/gbsg/run_gbsg_criteria.R (current; library(forestsearch) at :41; real survival::gbsg data; no call)
  - /home/larryleon/Documents/GitHub/fs-glms-interpretable/dev/manuscript_status_update.md (out of scope; :1045 used only to classify the SPEC)
  - /home/larryleon/Documents/GitHub/fs-glms-interpretable/sims/configs/opchar_survival_gbsg.R (out of scope; TODO stub, comment mention of setup_gbsg_dgm() at :2)
  - /home/larryleon/Documents/GitHub/fs-glms-interpretable/dev/tasks/REPORT_glms_structure_review_2026-09-11.md (out of scope; :275 quotes the sims/configs comment)
  - /home/larryleon/Documents/GitHub/fs-glms-interpretable/dev/reference/post_selection/manuscript_jrssb_initialsubmit/dgm_truth_targets.qmd (out of scope, reference copy; real calls :150, :156, :227, :230; text :108, :139; _h10.qmd, _h15.qmd and the rendered _h10.html are the same)
  - /home/larryleon/Documents/GitHub/fs-glms-interpretable/dev/reference/post_selection/manuscript_jrssb_initialsubmit/_archive/jasa-2026/postselection_subgroup_inference_main_v5_1.qmd (out of scope, archive; create_gbsg_dgm() text at :1312; also v5_6:1344, v5_7:1334, v5_7_blinded:1323, main_blinded:821, v5_8:832, v5_8_SUBMITTED:793, v5_8_blinded:783, and _archive/pre-trim/postselection_subgroup_inference_main.qmd:884)
  - /home/larryleon/Documents/GitHub/fs-glms-interpretable/dev/reference/post_selection/manuscript_jrssb_initialsubmit/postselection_subgroup_inference_supplement.qmd (out of scope, reference copy; create_gbsg_dgm text :1182; also _survey_2026-08/report_C_fb_mr_inventory.md:833)
  - /home/larryleon/Documents/GitHub/fs-glms-interpretable/dev/reference/post_selection/manuscript_jrssb_initialsubmit/_sim_noise_null_sweep.qmd (out of scope, reference copy; setup_gbsg_dgm comment :63; also _sim_noise_null_sweep_standalone2.qmd:63 and _sim_estimation_fs_h15_n1000.qmd:124)
- **Notes:** CONFIRMED CORE RESULT: no wrapper function, worker or nested helper is called or named in quarto/, manuscript/ or dev/payloads/ of fs-glms-interpretable. Every function and every directory has 0 calls and 0 name mentions, so there is no calling file to classify. The only near-mention in scope is the concept phrase 'via `k_inter` calibration' at quarto/gbsg/SPEC_gbsg_grf_dina_2026-09-17.md:102. That SPEC is current, not superseded, and the phrase describes forestsearch's gbsg_020 DGM under §2.5 'does not transfer'. Related output-field text, not function names: `dgm$hr_H_true`/`dgm$hr_Hc_true` and `flag_harm` in dev/payloads/identifiers/COLUMNS.md:34, :39 and in its producer dev/payloads/R/identifiers_metrics.R:557, :562. INDIRECT DEPENDENCE is wider than reported. The current manuscript renders from three forestsearch-produced simulation payload streams, and each was generated by code that calls calibrate_k_inter and setup_gbsg_dgm: (a) identifiers_metrics.csv, read at _05:14 and _S4:16, from gbsg_020 template bundles at ca9b0a62/69ca797e/dbc9444c/1507bfc0; (b) idsweep_metrics.csv, read at _02:4 and _S4:17, from the same template at b15c8b7a; (c) declcalc0_summary_tables.csv, read at _09:824, from scripts_dinamr/declcalc0_run.R and declcal_run.R. The GLM OC payloads (dev/payloads/glm, read by manuscript/_R/glm_oc_tables.R:31-40) come from ACTG175 GLM templates that use generate_glm_dgm/simulate_from_glm_dgm, not this wrapper. The applications payloads come from real-data analyses with no DGM. OUT OF SCOPE (not counted): 20 tracked files elsewhere in the repo carry a wrapper name. Real calls are only in dev/reference/.../dgm_truth_targets{,_h10,_h15}.qmd (plus the rendered _h10.html). The rest are text or comments in archived or reference post-selection drafts, sims/configs/opchar_survival_gbsg.R:2 and dev/tasks/REPORT_glms_structure_review_2026-09-11.md:275.
- **Verifier (corrected) corrections:**
  - Line number: get_dgm_with_output's @keywords internal tag is at R/sim_aft_gbsg.R:1075, not :1074. The other tags check out: @export at setup_gbsg_dgm.R:86 and sim_aft_gbsg.R:198, 788, 1003, 1156, 1228; @keywords internal at sim_aft_gbsg.R:50, 68, 235 (.create_gbsg_dgm_) and 822 (.simulate_from_gbsg_dgm_). NAMESPACE exports are at :48 S3method(print,gbsg_dgm), :92 calibrate_k_inter, :116 create_gbsg_dgm, :286 setup_gbsg_dgm, :293 simulate_from_gbsg_dgm and :312 validate_k_inter_effect.
  - Count: quarto/, manuscript/ and dev/payloads/ hold 10 tracked PDFs, not 11. All 10 are under manuscript/: cover_letter_jasa, fs_glms_interpretable{,_blinded}, fs_glms_interpretable_supplement{,_blinded}, and the five in submissions/jasa_2026-10-13/. pdftotext on each gives 0 hits.
  - Producing commits: the report cites only the c40d5309 template lines (:816, :821, :829). dev/payloads/identifiers/COLUMNS.md:32 names the forestsearch commits that actually produced the bundles: ca9b0a62/69ca797e (46 bundles), dbc9444c (17) and 1507bfc0 (18 s5rerun). At ca9b0a62 and 69ca797e the template calls calibrate_k_inter at :703 and setup_gbsg_dgm at :707. At dbc9444c the calls are at :740 and :744. At 1507bfc0 they are at :816, :821 and :829, the same lines as c40d5309.
  - Missed indirect dependence (idsweep): the manuscript also reads dev/payloads/identifiers/idsweep_metrics.csv, at manuscript/_02_identification.qmd:4 and manuscript/_S4_operating_characteristics_supplement.qmd:17. That CSV is derived from 288 idsweep bundles (dev/payloads/identifiers/idsweep/*.rds). README.md:30 and :40-42 say they were copied from forestsearch gbsg_020 at b15c8b7a. They come from renders of the same template: forestsearch b15c8b7a quarto/simulations/gbsg_020/scripts_dinamr/idsweep.sh:51 kills 'quarto render sim_fs_maxeffCons_fb_mr_field_m1_template.qmd'. At b15c8b7a that template calls calibrate_k_inter at :740 and setup_gbsg_dgm at :744.
  - Missed indirect dependence (declcal): manuscript/_payloads/declcal/declcalc0_summary_tables.csv is read at manuscript/_09_field_constructions_supplement.qmd:824 through manuscript/_R/declcal_tables.R (sourced at :823). It summarises forestsearch declcalc0/declcal payloads. Its provenance (declcalc0_summary_tables.provenance.md:3) cites scripts_dinamr/declcalc0_run.R at a46bf9b7, which calls calibrate_k_inter at :148 and setup_gbsg_dgm at :156 and :163 (same lines at c40d5309). The sibling declcal_run.R calls them at :132, :140 and :147.
  - Read site: oc_figures.R:3 is a header comment. The actual reads of identifiers_metrics.csv are manuscript/_05_operating_characteristics.qmd:14 and manuscript/_S4_operating_characteristics_supplement.qmd:16.
  - Missed text mentions of wrapper output fields (they are not function names, so the 0 name-mention counts stand): dev/payloads/identifiers/COLUMNS.md:39 and dev/payloads/R/identifiers_metrics.R:562 mention `dgm$hr_H_true` / `dgm$hr_Hc_true`. Those fields are computed in .create_gbsg_dgm_ (R/sim_aft_gbsg.R:527-548) and returned through setup_gbsg_dgm. COLUMNS.md:34 and identifiers_metrics.R:557 mention `harm_col <- "flag_harm"`, the column that .rename creates at R/setup_gbsg_dgm.R:159. identifiers_metrics.R is a current payload producer that writes COLUMNS.md (:2, :607), and the report did not name it.
  - Incomplete out-of-scope list. Every tracked hit for a wrapper name in the repo is in 20 files outside the three directories. The report left out these create_gbsg_dgm text mentions: _archive/jasa-2026/postselection_subgroup_inference_main_blinded.qmd:821, _v5_8.qmd:832, _v5_8_SUBMITTED.qmd:793, _v5_8_blinded.qmd:783, _archive/pre-trim/postselection_subgroup_inference_main.qmd:884, postselection_subgroup_inference_supplement.qmd:1182 and _survey_2026-08/report_C_fb_mr_inventory.md:833. It also left out the setup_gbsg_dgm comment mentions in _sim_estimation_fs_h15_n1000.qmd:124, _sim_noise_null_sweep.qmd:63 and _sim_noise_null_sweep_standalone2.qmd:63. In addition, dgm_truth_targets{,_h10,_h15}.qmd carry text mentions of calibrate_k_inter at :108 and :139 next to the real calls at :150, :156, :227 and :230. All of these files are under dev/reference/post_selection/manuscript_jrssb_initialsubmit/.

### What removal would involve (verified)

- **files to go:**
  - R/sim_aft_gbsg.R (1316 lines). Defines: constants SEED_BASE :29, DAYS_PER_MONTH :32, DEFAULT_N_SUPER :35; functions cut_numeric :51, cut_size :69, create_gbsg_dgm :199, .create_gbsg_dgm_ :236, simulate_from_gbsg_dgm :789, .simulate_from_gbsg_dgm_ :823, calibrate_k_inter :1004 (nested eff_at :1023), get_dgm_with_output :1076, print.gbsg_dgm :1157, validate_k_inter_effect :1229. Re-verified by grepping every 'function(' and top-level assignment: nothing else is defined.
  - R/setup_gbsg_dgm.R (190 lines). Defines setup_gbsg_dgm :88 and its nested .rename :149, nothing else.
  - man/create_gbsg_dgm.Rd
  - man/simulate_from_gbsg_dgm.Rd
  - man/calibrate_k_inter.Rd
  - man/validate_k_inter_effect.Rd
  - man/setup_gbsg_dgm.Rd
  - man/print.gbsg_dgm.Rd
  - man/get_dgm_with_output.Rd (keyword internal)
  - man/cut_numeric.Rd (keyword internal)
  - man/cut_size.Rd (keyword internal)
  - All 9 Rd pages are generated. .create_gbsg_dgm_ and .simulate_from_gbsg_dgm_ have no Rd page: their roxygen is an untitled '@keywords internal' at :235 and :822.
  - DESCRIPTION has no Collate field, so it needs no edit (Version 0.4.1.9007 at :3).
- **exports to go:**
  - NAMESPACE:48 S3method(print,gbsg_dgm)
  - NAMESPACE:92 export(calibrate_k_inter)
  - NAMESPACE:116 export(create_gbsg_dgm)
  - NAMESPACE:286 export(setup_gbsg_dgm)
  - NAMESPACE:293 export(simulate_from_gbsg_dgm)
  - NAMESPACE:312 export(validate_k_inter_effect)
  - History: all six are in v0.4.0:NAMESPACE (:43, :78, :102, :259, :266, :285), and all were exported at v0.1.0. NEWS.md:286 has an untagged 0.4.1 section; whether 0.4.1 reached CRAN cannot be told from the repo, but the wrapper files are identical, so the claim holds either way.
  - Deprecation: only create_gbsg_dgm (.Deprecated at R/sim_aft_gbsg.R:211-225) and simulate_from_gbsg_dgm (:800-812) are deprecated. calibrate_k_inter, validate_k_inter_effect, setup_gbsg_dgm and print.gbsg_dgm never were.
  - No import is lost. NAMESPACE:317 has import(survival), and generate_aft_dgm_main.R:189-190 declares the stats functions the wrapper used. No other R/ file registers anything for class gbsg_dgm.
- **help pages to go:**
  - Removed (9): create_gbsg_dgm.Rd, simulate_from_gbsg_dgm.Rd, calibrate_k_inter.Rd, validate_k_inter_effect.Rd, setup_gbsg_dgm.Rd, print.gbsg_dgm.Rd, get_dgm_with_output.Rd, cut_numeric.Rd, cut_size.Rd.
  - Regenerated, with links that must be fixed to avoid an Rd missing-link WARNING. This list is complete (grep of every \link{} in man/): man/calibrate_k_treat.Rd:59 and :115 (from R/calibrate_k_treat.R:9, :90); man/fs_fdr_report.Rd:182 (from R/fs_fdr_report.R:172); man/run_simulation_analysis.Rd:42 and :260 (from R/run_simulation_analysis.R:205, :339).
  - Regenerated, with stale \dontrun examples or text that name removed functions. Not run by check; every one is inside \dontrun: man/compute_dgm_cde.Rd:10 and :46 (from R/oc_analyses.R:109, :123); man/summarize_simulation_results.Rd:32 (:413); man/simulate_from_dgm.Rd:193 (R/simulate_from_dgm.R:167); man/check_censoring_dgm.Rd:72 (:436); man/run_simulation_analysis.Rd:204 and :242 (R/run_simulation_analysis.R:198, :373); man/fs_dgm_feasibility.Rd:76 (R/fs_dgm_feasibility.R:84).
  - Regenerated, with optional wording changes for the 'gbsg_dgm' class name: man/get_dgm_hr.Rd:10, man/invert_dgm_hrs.Rd:10, man/format_oc_results.Rd:54.
- **tests to go or change:**
  - tests/testthat/test-knoise-population.R: the helper .make_gbsg_noise_dgm (:7-9) and survival tests 1-5 (:32-96) have no skip and would ERROR under check. They can be ported without any new flex argument: build the flex object and draw noise onto dgm$df_super exactly as setup_gbsg_dgm.R:125-135 does (verified identical values). Test 4 (:71-79), 'perturbs nothing but the noise', then becomes a test of the caller-side helper. The GLM tests (:100-141) are unaffected.
  - tests/testthat/test-betaHhat-contract.R: the cached fixture .small_surv_dgm (:414-426) calls calibrate_k_inter (:419) and setup_gbsg_dgm (:421, seed 8316951L). T12 (:449-459) has no skip and would ERROR under check and on CRAN. T11 (:428), T13 (:461) and the tail of T14 (:502-505) are behind skip_if_not(.have_shims()); the shim files exist in quarto/, which is build-ignored (.Rbuildignore:24), so these run locally only. Rebuild with find_k_inter_for_target_hr(1.0, ...) plus generate_aft_dgm_flex on the caller-built GBSG frame: df_super lin_pred_*, theta and flag_harm then match exactly, so the fs_build_eval_frame assertions need no change. The shim quarto/simulations/gbsg_redux/betaHhat_truth.R reads no legacy-only field.
  - tests/testthat/test-plan-restored.R:164-173 calls fs_fdr_report(null = 'hom', n_super = 2000L) with skip_on_cran (:165) and a mocked %dofuture% (:166). It follows whatever fs_fdr_report becomes, and runs under devtools::test() and NOT_CRAN checks.
  - tests/testthat/test-fs-dgm-feasibility.R:22-28 needs no change if the stop() at R/fs_dgm_feasibility.R:128-131 keeps the substring 'Survival DGMs'. Its structure(list(), class = 'gbsg_dgm') stub (:26) is only a class string.
  - No test calls create_gbsg_dgm, simulate_from_gbsg_dgm, get_dgm_with_output, validate_k_inter_effect or print.gbsg_dgm. test-fs_mr_oc_summary.R:226 names run_simulation_analysis() only in a test title.
- **other package touchpoints:**
  - _pkgdown.yml: remove :332 setup_gbsg_dgm, :333 create_gbsg_dgm, :334 print.gbsg_dgm, :338 simulate_from_gbsg_dgm, :339 get_dgm_with_output, :342 calibrate_k_inter and :344 validate_k_inter_effect, all in 'Data Generating Mechanisms' (:326-345). pkgdown errors on missing topics. The file is build-ignored (.Rbuildignore:21).
  - NEWS.md: add a removal entry under :1 (0.4.1.9007). There is no current deprecation record; the old 'Superseded functions' note was deleted in fdecbe7b. Historical mentions stay: :511 (0.4.0 section) and :1088 (0.2.2 section).
  - R/globals.R:196-208: the GBSG globalVariables block; its header :197 cites a nonexistent file. 16 names are wrapper-only (re-verified): z3, z4, z5, zh, lin.conf.true, lin1.conf, lin0.conf, linC1.conf, linC0.conf, hlin.conf.1, hlin.conf.0, hlin.ratio, h1.potential, h0.potential, Ts, t.sim. Removing them is optional; leaving them raises no NOTE.
  - README.md:131 lists create_gbsg_dgm() in the key-functions table. README.Rmd is untracked and build-ignored (.Rbuildignore:14).
  - CLAUDE.md:53 (build-ignored, .Rbuildignore:33) describes setup_gbsg_dgm.R / sim_aft_gbsg.R as the GBSG DGM. Missed by the analysis; stale after removal.
  - inst/: 0 references across its 30 tracked files.
  - vignettes/forestsearch.qmd, the only CRAN-built vignette (DESCRIPTION:30): 0 references.
  - vignettes/articles/ (pkgdown-only, .Rbuildignore:20). paper_simulations.qmd has 7 evaluated wrapper calls (:507, :604, :614, :635, :645, :683, :699) and an implicit print.gbsg_dgm (:523); these make the pkgdown article build fail. After :440 it has about 49 lines that read legacy fields (:533-541, :561, :575, :622-624, :653-655, :668-670, :706-709, :1188-1203, :1596-1598 and others), some in eval: false or display blocks. Separately, the mid-document YAML at :358-360 sets the deployed page title (cosmetic). treatment_effect_definitions.qmd: eval=FALSE call :1034, plus text :163, :452, :460, :552, :785. extreme_subgroups.qmd :202, :485 and biomarker_effects.qmd :66 are prose. The tracked .html renders carry stale copies.
  - R CMD check, install: NAMESPACE must be regenerated. Stale export(...) lines or S3method(print,gbsg_dgm) for missing objects fail at install.
  - R CMD check, Rd: a missing-link WARNING unless R/calibrate_k_treat.R:9 and :90, R/fs_fdr_report.R:172 and R/run_simulation_analysis.R:205 and :339 are edited.
  - R CMD check, tests: ERROR from test-knoise-population.R survival tests 1-5 and test-betaHhat-contract.R T12 unless they are ported or removed. test-plan-restored.R:164-173 is skip_on_cran.
  - R CMD check, examples: every example that names a removed function is \dontrun (verified), so check reports nothing; the examples are only stale.
  - CRAN and API: all six names are exported in CRAN 0.4.0, so removal breaks users. Reverse dependencies were not checked. Four of the six were never deprecated. Saved gbsg_dgm objects would print with print.default. Saved setup_gbsg_dgm objects (class aft_dgm_flex) still work with simulate_from_dgm, compute_dgm_cde and fs_build_eval_frame, but have no df_source for baseline = 'fixed'.
  - quarto/ is build-ignored (.Rbuildignore:24). Still, 281 current quarto documents call setup_gbsg_dgm (281) and calibrate_k_inter (277) and would fail at run time; see current_documents.
- **dependent package functions:**
  - **fs_fdr_report()** (R/fs_fdr_report.R:181. Calls: setup_gbsg_dgm(model = 'alt', k_inter = 0, n_super) at :256; calibrate_k_inter(target_hr_harm = 1.0, model = 'alt', use_ahr = FALSE) at :258; setup_gbsg_dgm(..., k_inter = k_harm) at :262. Also: @seealso :172, .FS_GBSG_CONFOUNDERS :32-34, example reference values :163-164, analysis_time = 84 :189. Callers: tests/testthat/test-plan-restored.R:168 and one quarto document, quarto/calibration/fs_fdr_report_smoke.qmd:109 and :168. gbsg_aft_dgm_appendix.qmd only mentions it in a comment (:259).): Go together with the wrapper, or rebuild both nulls with flex. A rebuild is exact for the outcome model and truths: prepare y = rfstime/30.4375 and z1..z5 (the logic of the removed cut_numeric, moved to an internal helper); continuous_vars = z1..z5; subgroup er quantile 0.25 & meno = 0. Hom null: generate_aft_dgm_flex(model = 'alt', k_inter = 0, n_super). Harm null: find_k_inter_for_target_hr(1.0, tol = 1e-6) and then the same build. Verified: truths are identical to setup_gbsg_dgm. The censoring cannot be reproduced through flex arguments (legacy refit on the reassigned super-population, rescaled by sigma/sigma_cens, R/sim_aft_gbsg.R:604-646; flex fits on source data, helpers.R:214-336). Simulated censoring times therefore change while event times do not, and the documented A ~ 0.138 / 0.278 at n = 700 (:163-164) needs re-establishing. .FS_GBSG_CONFOUNDERS stay available through add_unprocessed_vars (helpers.R:563-589). Relink @seealso.
  - **.fs_fdr_hr_row() (internal to fs_fdr_report)** (R/fs_fdr_report.R:482-496: k_inter from dgm$effect_modifiers (:486); hr_H_true/hr_Hc_true with fallbacks hr$HR_H/HR_Hc (:491-492); AHR_H_true/AHR_Hc_true (:493-494).): Read hazard_ratios$harm_subgroup, $no_harm_subgroup, $AHR_harm and $AHR_no_harm. Flex stores these at helpers.R:1123-1141, and legacy objects already carry the same names (sim_aft_gbsg.R:666-673), so the switch works now for both. The fallbacks hr$HR_H and HR_Hc are set nowhere in R/. k_inter must come from the caller, because flex stores no k_inter or effect_modifiers (helpers.R:1320-1330). Without the change, every cell is NA on a flex DGM.
  - **run_simulation_analysis()** (R/run_simulation_analysis.R:389. The default confounders_base = c('v1',...,'v7') at :398 names columns that only .create_gbsg_dgm_ builds (R/sim_aft_gbsg.R:293-300). Also: header :6-7; roxygen :198, :205, :339; \dontrun :373.): Change the v1..v7 default, because a flex GBSG DGM has no v* columns. Retarget the docs and example. Truth reads go through get_dgm_hr (:524-527) and are used only in a verbose message (:529-538). For a flex null DGM, hr_Hc prints NA (see get_dgm_hr).
  - **simulate_from_dgm() and check_censoring_dgm()** (R/simulate_from_dgm.R:172 and :440; only their \dontrun examples (:167, :436) call setup_gbsg_dgm(model = 'null').): Docs only. They consume the generic contract: class check :189-190, df_super lin_pred_* :298-300, model_params mu/tau/censoring :310, :316-336. Flex additionally provides df_source for baseline = 'fixed' (:194-206), which setup objects lack.
  - **compute_dgm_cde()** (R/oc_analyses.R:131. Code :134 uses df_super_rand %||% df_super; flag detection :146. Doc :109-110 and :123 name create_gbsg_dgm; comment :12.): Docs only. The df_super_rand branch is then used only by DGM objects saved from older versions.
  - **summarize_simulation_results()** (R/oc_analyses.R:418; its \dontrun example :413 calls create_gbsg_dgm().): Docs only.
  - **get_dgm_hr() and format_oc_results()** (R/oc_analyses.R:49 (hazard_ratios first :52-66, top-level fallback :70-84; @param :41) and :618 (@param :582).): No code change is forced, but the fallback does not 'go dead'. Legacy objects already fill hazard_ratios (sim_aft_gbsg.R:666-673, :690), so the fallback is already redundant for them. It is still reached for a flex null DGM, which has no harm_subgroup or no_harm_subgroup and has AHR_harm = NaN (helpers.R:1104, :1134-1141). There get_dgm_hr(., 'hr_Hc') returns NA, while a legacy null returns hr_causal (sim_aft_gbsg.R:548, :672). Callers in simulation_tables already handle the null through model_type and hr_overall (:666-676). The 'gbsg_dgm' wording in the docs is optional to change.
  - **invert_dgm_hrs(), build_classification_table(), build_estimation_table(), interpret_estimation_table()** (R/simulation_tables.R:169 (:186-190), :272 (:311-330), :584 (:629-682), :1095 (:1147-1190; df_super_rand-only CDE fallback :1374-1380).): Nothing forced. They work through get_dgm_hr and df_super, and flex carries CDE in hazard_ratios (helpers.R:1107-1131).
  - **calibrate_k_treat()** (R/calibrate_k_treat.R:96; roxygen \link{calibrate_k_inter} at :9 and :90.): Docs only: retarget the two links, or check raises a missing-link WARNING. It already reproduces the quarto k_treat-on-hr_causal uniroot exactly (verified 0.6853424 = 0.6853424).
  - **.calibrate_by_root(); comments in calibrate_glm_interaction() and find_k_inter_for_target_hr()** (R/calibrate_helpers.R:47 (comment :4, @noRd text :22); R/calibrate_glm_interaction.R:142; R/find_k_inter_main.R:129.): Stay; edit the comments only. Still used by find_k_inter_for_target_hr (:131), calibrate_k_treat (:148) and calibrate_glm_interaction.
  - **fs_dgm_feasibility()** (R/fs_dgm_feasibility.R @note :82-86 and stop() :128-131.): Reword to refer to survival / aft_dgm_flex DGMs, keeping 'Survival DGMs' (asserted at test-fs-dgm-feasibility.R:28).
  - **fs_build_eval_frame() (R/betaHhat_truth.R)** (R/betaHhat_truth.R:820: dgm$df_super (:230, :268, :828-858), model_params$censoring (:219-233), simulate_from_dgm (:849). No wrapper reference.): No change.
  - **generate_aft_dgm_flex(), find_k_inter_for_target_hr(), sensitivity_analysis_k_inter(), calibrate_k_treat() (the replacements)** (R/generate_aft_dgm_main.R:218-241; R/find_k_inter_main.R:88-174, :234-277; R/calibrate_k_treat.R:96-153.): Nothing forces a change. What current documents still lack, re-ranked against the code and a run:
(1) A GBSG covariate-prep helper (months and z1..z5 as at R/sim_aft_gbsg.R:273-284, plus grade3 and v1..v7 if wanted). With it, flex reproduces the legacy outcome model and truths exactly.
(2) The legacy censoring model (:604-646), which no flex argument expresses (helpers.R:117-336).
(3) A print method for aft_dgm_flex; NAMESPACE has none.
(4) Optionally, a k_random_noise / noise_seed argument with noise_* fields. The same draw is a 3-line caller step that reproduces exactly.
(5) Optionally, an AHR target for k_inter. It is closed-form from model_params$b0.
Not needed by any current document: seed or '...' forwarding in find_k_inter_for_target_hr (every calibration runs at seed 8316951 and n_super 5000, its defaults), k_z3 (can be done in two steps with set_beta_spec) and use_rand_params (unused). To reproduce the legacy root closely, pass tol = 1e-6 to find_k_inter_for_target_hr (its default is 0.001, :99). Each of its evaluations also fits the AIC censoring set (select_censoring default TRUE); that costs time but does not affect HR(H).
- **current documents:**
  - **Group G1, forestsearch quarto/ (248 current files; list with call lines in /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2ro/LG-removal/groups.tsv): quarto/simulations/gbsg/ 104, gbsg/sims1_claude_summary/ 3, gbsg_020/ 56, gbsg_redux/ 85. The breakdown by prefix is as in the original analysis.** — step: calibrate_k_inter(target_hr_harm, model = dgm_model, use_ahr = FALSE). No seed or n_super is passed, so calibration runs at 8316951 / 5000. Then setup_gbsg_dgm(model, k_inter, n_super, seed = seed_base), compute_dgm_cde, and simulate_from_dgm(analysis_time = 84, cens_adjust = log(1.5)) and/or an eval frame. Example: quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_m1_h10_knoise0_n500_batch_1_50.qmd (build :318-331, confounders :291, simulate :468-469). Re-verified: all 248 use the raw 7 GBSG confounders. 241 of them read dgm$hr_causal / $hr_H_true / $hr_Hc_true through %||%. gbsg_redux knoise documents draw population noise caller-side onto dgm$df_super (e.g. gbsg_redux/mr_coverage_sweep_h10_knoise3.qmd:213-214). — flex + simulate_from_dgm offers the same step: no as a drop-in, but only the censoring step and the field names actually differ.

Yes, and exact, once the caller builds y (months) and z1..z5 as R/sim_aft_gbsg.R:273-284 does and passes them as continuous_vars with subgroup er quantile 0.25 & meno == 0:
- The outcome model has the same design and the same k_treat / k_inter scaling (helpers_outcome.R:45-76 vs sim_aft_gbsg.R:332-441).
- Subgroup: process_cutpoint quantile with <= (helpers.R:742-745), and meno <= 0 (:734-735).
- RNG order is the same (generate_aft_dgm_main.R:248; helpers.R:858, :1080).
- Truths come from helpers.R:1076-1159.
- Calibration: find_k_inter_for_target_hr (find_k_inter_main.R:88-174), whose defaults (seed 8316951, n_super 5000) match these calls.
- Run check: HR(H), HR(Hc), overall HR and AHR differences are 0; the calibrated k_inter differs by 2e-7 at tol = 1e-6.
- Caller-side population noise carries over unchanged.

No:
(1) Censoring. Legacy refits on the reassigned super-population and rescales by sigma/sigma_cens (R/sim_aft_gbsg.R:604-646); flex fits on source data (helpers.R:214-336). In the run, mu/tau were 4.049/0.436 vs 3.957/0.421, and only C_sim changes (t_true and treat_sim are identical per seed; mean event share at n = 700 was 0.524 vs 0.527).
(2) Field names. Truths are in hazard_ratios (helpers.R:1123-1141); the top-level names return NULL, which %||% turns into NA.
  - **Group G2, the population-noise sweeps (14 files), forestsearch quarto/simulations/gbsg_020/: maxeffCons_mr_coverage_sweep_h10_knoise0, h10_knoise3, h10_knoise3_fsonly, h10_knoise6, h15_knoise0_fsonly, h15_knoise3_fsonly_n1500, h20_knoise3_fsonly_n1500; *_standalone2 for h10_knoise0, h10_knoise3, h10_knoise6, h15_knoise0, h15_knoise3_n1500, h20_knoise3_n1500; and TEST_maxeffCons_mr_coverage_sweep_h10_knoise3_standalone2.** — step: As in G1, plus setup_gbsg_dgm(..., k_random_noise, noise_seed). The documents then assert dgm$noise_scheme, noise_names and noise_seed (e.g. maxeffCons_mr_coverage_sweep_h10_knoise3_standalone2.qmd:278-301). — flex + simulate_from_dgm offers the same step: no as a flex argument: generate_aft_dgm_main.R:218-241 has no k_random_noise or noise_seed. The step is nevertheless reproducible exactly by caller code. setup_gbsg_dgm draws after the build: set.seed(noise_seed), then rnorm(nrow) per column, with the RNG saved and restored (setup_gbsg_dgm.R:125-135). The same draw onto a flex dgm$df_super gives identical values (verified with identical(), n_super 25000), because the rows and their order are the same. This is the pre-0.2.2 driver-side scheme still used in gbsg_redux/mr_coverage_sweep_h10_knoise3.qmd:213-214. The caller would also have to set noise_names, noise_scheme and noise_seed for the asserts. Otherwise as G1, with censoring as the only simulation difference. For the 4 files with knoise 0, the noise step is a no-op.
  - **Group G3, z1_quantile variants (5 files): quarto/simulations/gbsg_020/gate_d2_cim_field.qmd (:694 / :698), gate_d2_cim_ij.qmd (:694 / :698), gate_d2_cim_unset.qmd (:698 / :702), scripts_dinamr/declcal_preflight.R (:30-31), scripts_dinamr/grf_mechanism.R (:75-77).** — step: calibrate_k_inter(target_hr_harm, model, use_ahr = FALSE, z1_quantile = harm_z1_quantile), with z1_quantile forwarded to .create_gbsg_dgm_ (R/sim_aft_gbsg.R:1024-1027). Then setup_gbsg_dgm(model, k_inter, z1_quantile, n_super, seed), plus k_treat = 1 in declcal_preflight.R, and simulate_from_dgm. — flex + simulate_from_dgm offers the same step: no as a drop-in, for the G1 censoring and field-name reasons. The z1_quantile step is a yes: the caller builds z1 at that quantile and sets subgroup_cuts er = list(type = 'quantile', value = q) (helpers.R:742-745), and find_k_inter_for_target_hr takes that same data and subgroup_cuts (find_k_inter_main.R:88-121). The truths are then exact, as in G1.
  - **Group G4, k_treat calibration on the null plus z1_quantile (5 files): quarto/simulations/gbsg_020/scripts_dinamr/declcal_run.R (:132-149), declcal_c0approx_run.R (:149-164), declcalc0_run.R (:148-163), nullid_designpoint.R (:22-26, :61-62), and quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_field_m1_template.qmd (:816-829). The template and declcalc0_run.R produce the fs-glms-interpretable payloads.** — step: k_treat comes from stats::uniroot over setup_gbsg_dgm(model = 'null', k_treat = kt, z1_quantile, n_super, seed)$hr_causal with tol = 1e-12 (declcal_run.R:139-145), or k_treat = 1e-10 for the complete null (:136-137). Alt designs use calibrate_k_inter(..., z1_quantile). Then setup_gbsg_dgm(model, k_inter, k_treat, z1_quantile, n_super, seed). The gate reads dgm$subgroup_info$fs_harm_true / grf_harm_true (declcal_run.R:165-166; template :892-893), dgm$model_params$b_hr['treat'] (declcal_run.R:176; template :911) and dgm$AHR (template :959; nullid_designpoint.R:56). Then fs_build_eval_frame and simulate_from_dgm. — flex + simulate_from_dgm offers the same step: no as a drop-in.

Yes:
- k_treat calibration through calibrate_k_treat(target, base_args with model = 'null', seed, n_super, tol = 1e-12) (R/calibrate_k_treat.R:96-153). Verified identical k_treat (0.6853424 at target 0.8).
- The null model (helpers.R:601).
- b_hr is model_params$b0 (:1288), and AHR is hazard_ratios$AHR.

No:
- Flex subgroup_info has no fs_harm_true or grf_harm_true (helpers.R:1294-1300), so the 'labels absent' gate passes vacuously.
- Censoring and field names, as in G1.
  - **Group G5, the León Table 1 reproduction (2 identical files): forestsearch quarto/simulations/reproduce_table1_m1.qmd and quarto/simulations/gbsg/reproduce_table1_m1.qmd.** — step: calibrate_k_inter(target_hr_harm = dgm_hr_h1, z1_quantile, k_treat = 1.0, n_super = dgm_n_super, seed = dgm_seed, verbose = TRUE) at :199-206, where dgm_n_super = 5000L and dgm_seed = 8316951L (:43-44). Then setup_gbsg_dgm(model = 'alt', ...) at :213 and setup_gbsg_dgm(model = 'null', ...) at :225, compute_dgm_cde, and run_simulation_analysis with confounders_base c('z1','z2','z3','z4','z5','size','grade3') (:77). — flex + simulate_from_dgm offers the same step: no overall. Calibration and truths: yes, exactly. The documents' seed and n_super equal find_k_inter_for_target_hr's fixed defaults (generate_aft_dgm_main.R:238; find_k_inter_main.R:101), so the lack of a seed argument does not matter here. z1..z5 and grade3 must be built in data, after which flex carries them into df_super (helpers.R:563-589). Simulated operating characteristics: no, because the legacy censoring model (R/sim_aft_gbsg.R:604-646) cannot be expressed in flex; only C_sim changes.
  - **Group G6, calibration notebooks and QC with the legacy covariates (5 files): quarto/calibration/calibration_survival_fourway.qmd (:260, :269, :275), selection_survival_frontier.qmd (:357, :375, :380), calibration_survival_leff_grid.qmd (:346), fpr_calibration_demo.qmd (:363; confounders :109), and quarto/qc/qc_null_classification.R (:31).** — step: calibrate_k_inter(target_hr_harm, model = 'alt', k_treat = 1.0, cens_type = 'weibull', use_ahr = FALSE, verbose = FALSE), then setup_gbsg_dgm(model = 'alt', k_treat = 1, k_inter) and/or setup_gbsg_dgm(model = 'null', k_treat = 1.0). These are followed by run_simulation_analysis or simulate_from_dgm with confounders z1..z5, size and grade3. — flex + simulate_from_dgm offers the same step: no, with the same gap as G5: covariates must be built by the caller, after which truths and calibration are exact, and censoring cannot be reproduced. The null build is a yes (helpers.R:601). cens_type does not affect HR(H): the censoring fit consumes no RNG and comes after the HR step (R/sim_aft_gbsg.R:617-658). The calibration is therefore covered by find_k_inter_for_target_hr.
  - **Group G7 (2 files): quarto/calibration/fs_fdr_report_smoke.qmd, which calls setup at :76 and calls fs_fdr_report at :109 and :168 (the latter '#| eval: true'); and quarto/calibration/gbsg_aft_dgm_appendix.qmd (setup :192, :195; calibrate :193), which does NOT call fs_fdr_report and has only a comment at :259.** — step: Smoke: setup_gbsg_dgm(model = 'alt', k_inter = 0, n_super = 25000L), simulate_from_dgm, forestsearch(), then fs_fdr_report(), which itself calls the wrapper (R/fs_fdr_report.R:256-262). Appendix: setup_gbsg_dgm(alt, k_inter = 0), then calibrate_k_inter(1.0) and setup_gbsg_dgm(alt, k_harm). It reads dgm$effect_modifiers, hr_H_true, hr_Hc_true, AHR_H_true, AHR_Hc_true and subgroup_info$proportion (:200-212), and dgm_hom$cens_params[c('type','mu','sigma')] (:227). — flex + simulate_from_dgm offers the same step: no. The smoke document follows whatever fs_fdr_report becomes. In the appendix the builds and truths are exact under flex, but a flex object has no effect_modifiers, no stored k_inter and no cens_params. Censoring is in model_params$censoring with mu/tau/gamma/type (helpers.R:1284-1291), and its values differ from the legacy model.
  - **forestsearch vignettes/articles/paper_simulations.qmd (pkgdown-only article).** — step: Evaluated calls:
- setup_gbsg_dgm(alt, k_treat 1, k_inter 2, k_z3 1, z1_quantile 0.25, n_super 5000, weibull, seed 8316951) at :507, then compute_dgm_cde and print(dgm_alt) at :523.
- calibrate_k_inter(2, use_ahr = FALSE) at :604, then setup at :614.
- calibrate_k_inter(2, use_ahr = TRUE) at :635, then setup at :645.
- validate_k_inter_effect(c(-2,-1,0,1,2,3)) at :683.
- setup(null) at :699.
- simulate_from_dgm(n = 700, analysis_time = 84, cens_adjust = log(1.5)) at :725.
The article then reads legacy fields on about 49 lines after :440. — flex + simulate_from_dgm offers the same step: no overall. Step by step:
(a) Build at :507: yes, exact, once the caller builds the covariates; k_z3 = 1 is a no-op.
(b) print at :523: no; NAMESPACE has no print method for aft_dgm_flex.
(c) :604: yes, via find_k_inter_for_target_hr (pass tol = 1e-6 for the same root).
(d) :635, use_ahr = TRUE: no packaged function. A closed form exists: loghr_po is constant in H (helpers.R:1055-1057), so k_inter = (log 2 - b0[treat]) / b0[treat_harm] from model_params$b0. Verified to equal the legacy 1.301622.
(e) :683: partly. sensitivity_analysis_k_inter(c(-2,3), n_points = 6) gives the same k grid with HR(H), HR(Hc) and overall (find_k_inter_main.R:246, :268-274), but no AHR columns or PASS checks (R/sim_aft_gbsg.R:1282-1312).
(f) :699: yes.
(g) Field reads: rename to hazard_ratios$... and df_super.
(h) :725: simulate_from_dgm already is the flex API, but censoring differs.
The article's 'equivalent' block (:395-440, eval: false) uses median splits as factors, not the wrapper's quartile codes (R/sim_aft_gbsg.R:281-284), so it is not equivalent.
  - **forestsearch vignettes/articles/treatment_effect_definitions.qmd, chunk worked-example :1030 (eval=FALSE).** — step: setup_gbsg_dgm(alt, k_treat 1, k_inter 2) at :1034. AHR and CDE from dgm$df_super_rand indexed by flag_harm (:1042-1060). compute_dgm_cde. Reads hr_causal / hr_H_true / hr_Hc_true (:1069-1071). simulate_from_dgm(dgm, n = 500, seed = 1) at :1075. — flex + simulate_from_dgm offers the same step: yes for the effect definitions, with the caller building the covariates. AHR, CDE and marginal HRs come from hazard_ratios (helpers.R:1103-1141), and df_super carries theta_0, theta_1, loghr_po and flag_harm (:1053-1057). The replicate at :1075 differs only in censoring. Existing bug, confirmed: df_super_rand has flag.harm, not flag_harm (R/sim_aft_gbsg.R:290; setup_gbsg_dgm.R renames only the df_super copy, :147-161), so :1042-1045 and :1052-1055 currently subset by NULL. Reading flex df_super would fix it.
  - **fs-post-selection (main 5d165ac) manuscript_rev1/dgm_truth_targets.qmd, dgm_truth_targets_h10.qmd and dgm_truth_targets_h15.qmd: standalone auxiliary documents, not included in the main paper or supplement. dgm_truth_targets_h10.html is a rendered copy.** — step: calibrate_k_inter(target, model = 'alt', use_ahr = FALSE, cens_type, z1_quantile) at :150 and :227, with no n_super or seed, so calibration runs at 5000 / 8316951. setup_gbsg_dgm(model, k_inter, z1_quantile, n_super = N_SUPER = 25000L, cens_type, seed = SEED = 8316951L) at :156 and :230 (constants :63-66). compute_dgm_cde. Reads hr_*_true, hr_causal, cde_*, AHR_*_true and AHR, and the flag column via intersect(flag.harm, flag_harm) (:164-183). Then model_params$gamma['zh'] and $sigma (:233-246). No simulation. — flex + simulate_from_dgm offers the same step: yes, with numerically identical output (previously 'unknown'). find_k_inter_for_target_hr with tol = 1e-6 uses defaults that match this calibration (n_super 5000, seed 8316951). generate_aft_dgm_flex(n_super = 25000, seed = 8316951) on caller-built z1..z5 then gives the same HR, AHR and CDE: in the run, the differences were 0 at n_super 25000 and the calibrated k agreed within 2e-7. compute_dgm_cde works through df_super (R/oc_analyses.R:134). Renames: hazard_ratios$harm_subgroup and related; gamma['treat_harm'] for gamma['zh']; model_params$tau for $sigma (helpers_outcome.R:88-93).
  - **Generators of the fs-post-selection manuscript_rev1 payloads (a subset of G1 and G2), in forestsearch quarto/simulations/gbsg_020/: sim_fs_maxeffCons_fb_mr_m1_h10_knoise0_n500_{batch_1_50, 51_100, 101_200, 201_300, 301_500, combine_1_500}; ..._h10_knoise0_n1000_{batch_1_20, 21_100, 101_110, 111_500, combine_1_500}; ..._h15_knoise0_n1000_{batch_1_200, 201_220, 221_300, 301_320, 321_400, combine_1_400}; sim_dina_eff_fb_mr_m1_h10_knoise0_n500_{batch_1_200, 201_250, 251_500, combine_1_500}; sim_grf_eff_..._{same}; maxeffCons_mr_coverage_sweep_{h075,h15}_knoise0 (G1); maxeffCons_mr_coverage_sweep_h10_knoise0 (G2); maxeffCons_mr_coverage_sweep_*_standalone2 x6 (G2).** — step: calibrate_k_inter(use_ahr = FALSE), then setup_gbsg_dgm. Call lines: :318/:321 (fs, grf), :316/:319 (dina), :247/:249 (h075, h15), :276/:278 (h10_knoise0) and :278/:280 (standalone2, which adds k_random_noise and noise_seed); spot-checked. Then compute_dgm_cde and simulate_from_dgm(analysis_time = 84, cens_adjust = log(1.5)). The calls run only in batch mode or when run_sweep is set; combine documents set run_mode 'combine' (e.g. :113). — flex + simulate_from_dgm offers the same step: no as a drop-in. Regenerating with flex would leave the truth targets, the calibrated k_inter (to the uniroot tolerance), the event times, the treatment draws and the population noise unchanged; the noise needs the caller-side draw. It would change the censoring times (R/sim_aft_gbsg.R:604-646 vs helpers.R:214-336), and with them the observed data and operating characteristics. In the run, the event share moved by about 0.3 percentage points at n = 700.
  - **fs-post-selection manuscript_rev1 rendered documents: postselection_subgroup_inference_main.qmd, postselection_subgroup_inference_supplement.qmd, postselection_theory_derivations.qmd and their included fragments.** — step: No wrapper call; they read _payloads/*.rds. The supplement's code-availability sentence names create_gbsg_dgm() at postselection_subgroup_inference_supplement.qmd:1184 (verified); it also appears in the kept .tex at :3296 (verified) and in the PDF. — flex + simulate_from_dgm offers the same step: not applicable. The sentence would name a removed function and would need rewording.
  - **fs-post-selection analysis/simulations/gbsg/resampling_comparisons/ (18 tracked .qmd files with code calls; missed by the analysis). Examples: m1_source/fs_t1_t2_m1_h10_knoise0_n500_combine_1_500.qmd, m1_source/dina_t1_t2_m1_h15_knoise5_n1000_batched_1_100.qmd, fs2_multiplier-resampling_vs_bootstrap_debias_n1000_simulation.qmd. Last commit touching the directory: 5f0ea31 on 2026-06-21. No manuscript_rev1 .qmd references it. Its currency is unclassified.** — step: calibrate_k_inter and setup_gbsg_dgm in the G1 / gbsg_redux pattern (calls located with grep over the 25 tracked .qmd/.R files in analysis/; lines not itemized here). — flex + simulate_from_dgm offers the same step: no as a drop-in, for the same reasons and with the same exact truths as G1. These files would fail at run time after removal.
  - **fs-glms-interpretable (main dc94a17): current files under quarto/, manuscript/ and dev/payloads/ (576 tracked files).** — step: No file calls or names a wrapper function. The dependence is indirect, through payloads produced by G4 code: identifiers_metrics.csv and idsweep_metrics.csv come from sim_fs_maxeffCons_fb_mr_field_m1_template.qmd; declcalc0_summary_tables.csv comes from scripts_dinamr/declcalc0_run.R. — flex + simulate_from_dgm offers the same step: not applicable. Regenerating those payloads would raise the G4 answer: no as a drop-in, with truths exact and censoring different.
- **Notes:** Verdict: corrected. The inventory checks out against the source: definitions, the 6 NAMESPACE lines, the 9 Rd pages and the full set of man/ links and mentions, the test files and their gating, every dependent package function, the pkgdown, README, NEWS and globals line numbers, the G1-G7 grouping (281 = 248 + 14 + 5 + 5 + 2 + 5 + 2), the 241 G1 files that read legacy truth fields, and the fact that all 248 G1 files use the raw GBSG confounders. The corrections concern the yes/no answers and a few dependency statements.

**Main correction.** Both 'unknown (not run)' answers are now settled, and several 'no' answers were too strong. A short Rscript --vanilla run sourced the c40d5309 R/ files (scripts in /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2ro/LH-adv/: equiv.R, equiv2.R, equiv3.R, equiv4.R; each ran in under 10 s). With the caller building y = rfstime/30.4375 and z1..z5 exactly as R/sim_aft_gbsg.R:273-284 does, generate_aft_dgm_flex(continuous_vars = z1..z5, subgroup er quantile 0.25 & meno == 0) reproduces setup_gbsg_dgm exactly:
- gamma, mu and tau are identical, and the flag_harm rows are identical;
- HR(H), HR(Hc), overall HR and AHR(H) differ by 0 at n_super 5000 and at 25000;
- calibrated k_inter agrees to 2e-7 at tol = 1e-6;
- k_treat calibration is identical;
- AHR-target calibration has a closed form that matches;
- the population noise drawn onto df_super after the build is identical.
For the same seed, simulate_from_dgm then gives identical event times and treatment. Only the censoring times differ, because the legacy censoring model (refit on the reassigned super-population and rescaled by sigma/sigma_cens, R/sim_aft_gbsg.R:604-646) cannot be expressed through flex arguments (helpers.R:117-336). In the run its parameters were mu 4.049 / tau 0.436, against 3.957 / 0.421 (force-fit) or 4.031 / 0.425 (AIC). The mean event share at n = 700 was 0.524 vs 0.527.

**Real gaps for current documents:**
1. A GBSG covariate-prep helper. This is the function of cut_numeric and cut_size, which would move out of the wrapper.
2. The legacy censoring model. This is the only gap that changes simulated data.
3. The field names. Documents read top-level hr_* through %||%, so the truths silently become NA. A flex null has no no_harm_subgroup, so hr_Hc must come from hazard_ratios$overall.
4. A print method for aft_dgm_flex.
5. Convenience only: a noise argument with the noise_* fields, an AHR target and fs_harm_true labels. Each is a short caller step.
Seed and '...' forwarding in find_k_inter_for_target_hr are not needed: every current calibration runs at 8316951 / 5000, its defaults. k_z3 and use_rand_params are unused; k_z3 could be done in two steps with set_beta_spec.

**Other corrections:**
- fs_fdr_report has one quarto caller (fs_fdr_report_smoke.qmd:109, :168), not two.
- get_dgm_hr's fallback does not go dead: it is still reached for flex null DGMs. Legacy objects already carry hazard_ratios.
- .fs_fdr_hr_row can switch to hazard_ratios now, for both kinds of object.
- paper_simulations.qmd has about 49 legacy-field lines after :440, not about 30. Its build fails because of the evaluated calls, not the title YAML.
- Missed document set: fs-post-selection analysis/simulations/gbsg/resampling_comparisons/ (18 calling files, last touched 2026-06-21).
- Missed touchpoint: CLAUDE.md:53.
- NEWS.md:286 has an untagged 0.4.1 section, so whether 0.4.1 reached CRAN cannot be told from the repo.

**Package removal work** is otherwise as stated:
- delete the 2 R files and their 9 Rd pages, and drop the 6 NAMESPACE lines;
- port test-knoise-population.R survival tests 1-5 (which can use the caller-side noise draw) and the test-betaHhat-contract.R fixture;
- decide on fs_fdr_report and .fs_fdr_hr_row, and change the v1..v7 default of run_simulation_analysis (:398);
- fix 5 Rd links and the stale \dontrun examples;
- remove 7 _pkgdown.yml lines, 1 README row and 16 orphan globals, and add a NEWS entry.

The G5 seed and n_super were checked and equal the flex defaults (reproduce_table1_m1.qmd:43-44).

The .Deprecated replacement call at R/sim_aft_gbsg.R:215-223 cannot run. id = and confounders = are unused arguments (no '...', main.R:218-241), and continuous_vars and factor_vars are missing; outcome, event and treat would partially match.

Everything was read-only. The only files written are the scratch scripts (equiv*.R, args.py) and three file lists (g1.txt, g1_conf.txt, psa.txt) in /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2ro/LH-adv/. The mrct-consistency repository was not searched.
- **Verifier (corrected) corrections:**
  - Numeric identity was marked 'unknown (not run)' for G1, dgm_truth_targets*.qmd and the payload generators; it is now verified. The caller builds y = rfstime/30.4375 and z1..z5 as at R/sim_aft_gbsg.R:273-284, then calls generate_aft_dgm_flex with continuous_vars = z1..z5, factor_vars = character(0) and subgroup er quantile 0.25 & meno = 0. This reproduces the legacy fit exactly. gamma, mu and tau are identical. Max |lin_pred_1 difference| is 2.2e-16. The legacy minus flex differences in HR(H), HR(Hc), overall HR and AHR(H) are exactly 0, at n_super 5000 and at 25000. Why it matches, from the code: the flex model is treat + z_* + treat_harm (generate_aft_dgm_helpers_outcome.R:45-68), the same design as the legacy treat + z1..z5 + zh (sim_aft_gbsg.R:332-379); and the RNG order set.seed, sample, rexp is the same (generate_aft_dgm_main.R:248, helpers.R:858, :1080 vs sim_aft_gbsg.R:478, :481, :498). Script: /tmp/claude-1000/-home-larryleon-Documents-GitHub-mrct-consistency/fcecc1eb-8fe3-4885-9d74-58531b66b71d/scratchpad/t2ro/LH-adv/equiv.R and equiv3.R.
  - Calibration equivalence is also verified. calibrate_k_inter(1.5) = 1.111548, and find_k_inter_for_target_hr(1.5, tol = 1e-6) on the same caller frame differs by 2.2e-7; with its default tol = 0.001 (find_k_inter_main.R:99) it gives 1.111349. The null k_treat uniroot over setup_gbsg_dgm(...)$hr_causal and calibrate_k_treat(0.8, base_args) both give 0.6853424. Script: .../t2ro/LH-adv/equiv2.R.
  - The censoring gap is confirmed and now quantified; it is the only simulated quantity that changes. Legacy Weibull censoring: mu 4.049, tau 0.436. Flex force-fit: 3.957 / 0.421. Flex AIC default: 4.031 / 0.425. lin_pred_cens differs by up to 0.216. For the same seed, simulate_from_dgm gives identical t_true and treat_sim (simulate_from_dgm.R:251-311); only C_sim differs. Over 200 seeds at n = 700, analysis_time 84 and cens_adjust log(1.5), the mean event share was 0.524 (legacy) vs 0.527 (flex force-fit).
  - G2 population noise: 'no' was stated without qualification. generate_aft_dgm_flex has no k_random_noise argument (generate_aft_dgm_main.R:218-241), but the draw is a caller-side post-build step. set.seed(noise_seed) followed by rnorm(nrow(dgm$df_super)) per column onto the flex dgm$df_super reproduces setup_gbsg_dgm's values exactly (identical() TRUE): setup draws after the build (setup_gbsg_dgm.R:125-135), and the row order is the same. Current G1 documents already do it this way (quarto/simulations/gbsg_redux/mr_coverage_sweep_h10_knoise3.qmd:213-214). The analysis's argument about source-data noise being duplicated by resampling does not apply. What flex lacks is only the argument and the noise_names / noise_scheme / noise_seed fields that the G2 drivers assert.
  - AHR target for k_inter (paper_simulations.qmd:635, and the replacement gap list): no packaged function, confirmed. But loghr_po is constant within H (generate_aft_dgm_helpers.R:1055-1057; within-H sd 8.9e-17), so AHR(H) = exp(b0[treat] + b0[treat_harm]) and k_inter = (log(target) - b0[treat]) / b0[treat_harm] can be read in closed form from a k_inter = 1 flex build's model_params$b0. This gives 1.301622, equal to calibrate_k_inter(2, use_ahr = TRUE). Script: .../t2ro/LH-adv/equiv4.R.
  - 'Seed and extra-argument forwarding in find_k_inter_for_target_hr' is not needed by any current document. Every calibrate_k_inter call in the 281 quarto files and in dgm_truth_targets*.qmd calibrates at seed 8316951 and n_super 5000, which are find_k_inter's defaults (find_k_inter_main.R:101; generate_aft_dgm_main.R:238). G1-G4, G6 and G7 pass no seed or n_super. G5 passes dgm_seed = 8316951L and dgm_n_super = 5000L (reproduce_table1_m1.qmd:43-44). z1_quantile is expressible through data and subgroup_cuts.
  - k_z3: 'no counterpart' is too strong. The fitted b0 does not depend on k (helpers_outcome.R:62-72), so a second build with set_beta_spec = list(set_var = 'z_z3', beta_var = k_z3 * b0['z_z3']) reproduces it. No current document uses a non-unit k_z3.
  - fs_fdr_report callers: the analysis said test-plan-restored.R:168 plus two quarto calibration documents. Only one quarto document calls it: quarto/calibration/fs_fdr_report_smoke.qmd:109, and :168 (the chunk is '#| eval: true' although its prose says not evaluated). gbsg_aft_dgm_appendix.qmd names it only in a commented line in an eval: false chunk (:255-259). The G7 entry is corrected to match: the appendix calls the wrapper directly and never calls fs_fdr_report.
  - get_dgm_hr: the claim that the legacy fallback 'becomes dead' is wrong. Legacy objects already fill hazard_ratios with the same names (sim_aft_gbsg.R:666-673, :690), so the fallback is already redundant for them. It stays reachable for a flex null DGM, whose hazard_ratios has no harm_subgroup / no_harm_subgroup (helpers.R:1134-1141) and has AHR_harm = NaN. So get_dgm_hr(flex_null, 'hr_Hc') returns NA where a legacy null returns hr_causal (sim_aft_gbsg.R:548, :672). build_estimation_table covers this through model_type (simulation_tables.R:666-676), and run_simulation_analysis uses the value only in a verbose message (:524-538), so nothing breaks. Documents that read hr_Hc_true on a null DGM must switch to hazard_ratios$overall.
  - .fs_fdr_hr_row: legacy objects already carry hazard_ratios$harm_subgroup, $no_harm_subgroup, $AHR_harm and $AHR_no_harm (sim_aft_gbsg.R:666-673). Switching it to hazard_ratios therefore works for both kinds of object today. Only k_inter has to come from the caller.
  - paper_simulations.qmd: there are about 49 lines after :440 that read legacy fields (some in eval: false or display blocks), not about 30. The pkgdown build fails because of the 7 evaluated wrapper calls; the mid-document YAML title (:358-360) is a separate, cosmetic issue, not the cause.
  - A document set was missing from current_documents. fs-post-selection analysis/simulations/gbsg/resampling_comparisons/ has 18 tracked .qmd files with code calls to calibrate_k_inter / setup_gbsg_dgm (e.g. m1_source/fs_t1_t2_m1_h10_knoise0_n500_combine_1_500.qmd). It was last touched 2026-06-21 (5f0ea31) and is not referenced from manuscript_rev1. Its currency is unclassified; it would be affected like G1.
  - A touchpoint was missing: CLAUDE.md:53 (build-ignored, .Rbuildignore:33) describes setup_gbsg_dgm.R / sim_aft_gbsg.R as the GBSG DGM.
  - Minor: NEWS.md:286 has a '# forestsearch 0.4.1' section, but there is no v0.4.1 tag and the CRAN-SUBMISSION record covers only 0.4.0 (83b6d708). Whether 0.4.1 reached CRAN cannot be told from the repo. The wrapper files are identical across these versions, so the export claims hold either way.

