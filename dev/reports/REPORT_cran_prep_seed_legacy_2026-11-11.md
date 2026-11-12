# REPORT — CRAN preparation, part 2: the fixed seed, and a usage inventory of the legacy GBSG wrapper, 2026-11-11

Task: `dev/tasks/TASK_cran_prep_seed_legacy_2026-11-11.md`, run unattended on 2026-11-11, on branch `cran-prep`. Nothing was pushed or merged.

## 1. Provenance

- **Host and tools.** pop-os (128 cores), R 4.6.2, devtools 2.5.2, testthat 3.3.2, pkgload 1.5.3, future 1.80.0, doFuture 1.4.0, gt 1.4.0.
- **Step 0.**
  - `git status --porcelain` was empty. **PASS.**
  - **Branch `cran-prep` did not exist,** so the first CRAN task (`dev/tasks/TASK_cran_prep_2026-11-11.md`, from `~/Downloads`) was run first, completely and with all its own gates. Its report `dev/reports/REPORT_cran_prep_2026-11-11.md` was committed on `cran-prep` as `3fbbc55c`. **PASS.**
  - That run's start was `e92d3c96`: its TASK commit on `master` (`afe96278` + the TASK file), from which `cran-prep` was created.
  - Its report is also at `~/Downloads/REPORT_cran_prep_2026-11-11.md`. Part 1 of this report summarises its outcome; the full record is that report.
- **Start SHA of this task:** `3fbbc55c`, the HEAD before this task's TASK commit. `DESCRIPTION` version at the start: **0.4.1.9007**.
- **This task's commits on `cran-prep`:**

| SHA | Commit |
|---|---|
| `a977f3d7` | dev: add CRAN preparation task, part 2 |
| `55e59792` | fix: forestsearch functions restore the caller's random-number state |
| (this commit) | dev: CRAN preparation report, part 2 |

- **Sibling repositories (read-only throughout).** Both were clean before (20:08) and after, with no file newer than the run's start marker:
  - fs-post-selection: `5d165acc` on `main`;
  - fs-glms-interpretable: `dc94a17e` on `main`.
- **How the work was done.**
  - **Part 1.** Implemented here. The `set.seed()` inventory was classified by four agents, each with an independent verifier. A separate read-only agent audited the future bodies for multisession effects.
  - **Parts 2 and 3.** One agent per vignette file and six collectors for the wrapper inventory, then a removal analysis. Each was re-derived by an independent adversarial verifier. The tables here are condensed from the verified results.
  - **Detail files** (committed with this report):
    - `dev/reports/REPORT_cran_prep_seed_legacy_2026-11-11_detail.md`: every verified row of Parts 1.1, 2 and 3, with the verifiers' corrections;
    - `dev/reports/REPORT_cran_prep_seed_legacy_2026-11-11_quarto_callers.csv`: every tracked `quarto/` file that calls the wrapper.

## 2. Status

- **First CRAN task (run as Step 0): DONE.** F1–F7 fixed (F2 narrowed to the both-unset case), the `paper_simulations` chunk fixed and the article re-rendered. Its suite and `check()` passed: 641 blocks, 8,455 expectations, 0 failures; 0/0/0. See `REPORT_cran_prep_2026-11-11.md`.
- **Part 1 — the caller's random-number state is restored: DONE (`55e59792`).**
  - The 15 exported functions that can call `subgroup.consistency()` restore `.Random.seed` at the outermost call's exit.
  - The equivalence gate passes: 6/6 calls identical, and the state is restored after each.
  - Under a multisession plan, `forestsearch_tenfold()` gives the same values as before.
  - Tests pass. Suite: 648 blocks, 8,471 expectations, 0 failures. `check()` 0/0/0.
- **Part 2 — vignettes and articles relying on the old state: DONE (read-only).**
  - 7 places, in 2 pkgdown-only articles. **None is in the vignette CRAN builds.**
  - Every one is a `gt` table's hidden random HTML id. No visible number, cell or figure depends on them.
- **Part 3 — the legacy GBSG wrapper: DONE (read-only).** The headline numbers and Larry's decision table are in §5.3.

## 3. Part 1 — the caller's random-number state is restored

Line numbers in §3.1, §4 and §5 are at the start SHA `3fbbc55c`. Commit `55e59792` inserted lines into 13 `R/` files (§3.2a), so in those files lines after an insertion are higher at the branch head. For example, `R/fs_fdr_report.R:256` at the start is `:257` at the branch head. Lines quoted "at the branch head" are at `55e59792`.

### 3.1 Inventory (Step 1): every `set.seed(` call in package code under `R/`

Method. A scan of `R/*.R` at the start SHA `3fbbc55c` for `set.seed(`, excluding roxygen and comment lines, found 44 sites in 29 files. Each was classified by reading the enclosing function and its callers: four agents took 11 sites each, and each was re-derived by an independent verifier (27 corrections, all folded in). The full rows are in the detail file.

- Seed column. "Fixed" means the user cannot change the value from any exported entry point that reaches the site. "Argument" means the value comes from a formal an exported caller exposes.
- "Exported functions that reach it" are those with a call path to the site, with the restrictions the verifiers recorded.
- Line numbers are at `3fbbc55c`; the fixed seed is row 41.

Counts:
- Seeds: 1 fixed outright (`lasso_selection`, because `seedit` is never forwarded to it), 29 argument, 10 argument-derived, 4 stored in an object.
- Seven sites are fixed on one route and an argument on another (rows 18, 21, 32, 33, 40, 41, 44). Row 41 is one: `subgroup.consistency()`'s seed is fixed at 8316951 from `forestsearch()` and every wrapper, and an argument only on a direct call.
- Restored before this change: 4 yes, 2 "yes, when seed is non-NULL", 1 partly, 37 no.

| # | Function | file:line | Seed: fixed or from an argument | Previous state restored? | Exported functions that reach it |
|---|---|---|---|---|---|
| 1 | bootstrap_results | R/bootstrap_analysis_dofuture.R:334 | argument: formal `seed = 8316951L`, exposed by forestsearch_bootstrap_dofuture(seed) via `boot_seed`; no is.null guard, so `seed = NULL` re-initialises the RNG from entropy. | no: `set.seed(seed)` runs unconditionally on the main process; the caller's .Random.seed is never saved (on.exit restores only options(warn) and future::plan). | bootstrap_results (direct); forestsearch_bootstrap_dofuture (via `bootstrap_results(seed = boot_seed)`). |
| 2 | .make_boot_index_matrix | R/bootstrap_calculations_helpers.R:27 | argument: helper formal `seed = NULL` (guarded by !is.null), set from bootstrap_ystar(`seed = 8316951L`) or forestsearch_bootstrap_dofuture(`seed = 8316951L`); the user can change it from both. | no: neither the helper nor its two callers saves or restores .Random.seed. | bootstrap_ystar; forestsearch_bootstrap_dofuture. |
| 3 | consistency_resample | R/consistency_resample.R:459 | argument: formal `seed = NULL` (guarded), used only when method is "mc" or "both"; set on a direct call or from consistency_resample_compare(seed). | no: nothing saves or restores .Random.seed. | Effectively only consistency_resample (direct) and consistency_resample_compare. forestsearch, subgroup.consistency, evaluate_subgroup_consistency, evaluate_consistency_twostage and the forestsearch callers reach the function only through method = "closed" calls with no seed, so this set.seed never executes from them; `fs_declaration_*`, fs_mr_inference, fs_family_report and `fs_oc_*` do not reach the function at all. |
| 4 | consistency_resample_compare | R/consistency_resample.R:505 | argument: formal `seed = NULL` (guarded by !is.null); no package code calls this function. | no: nothing saves or restores .Random.seed. | consistency_resample_compare only (no internal callers). |
| 5 | .dina_make_bootstrap_indices | R/dina_bagged.R:436 | argument: helper formal `seed = NULL` (guarded), set from dina_fit_bagged(`seed = NULL`), dina_bagged (via `...`) or dina_subgroup_bootstrap(`seed = NULL`). | yes, when seed is non-NULL: on.exit saves the global .Random.seed and assigns it back, or rm()s it if none existed at entry. | dina_fit_bagged; dina_bagged (via `...` to dina_fit_bagged); dina_subgroup_bootstrap. |
| 6 | .dina_make_folds | R/dina.R:433 | argument (mixed provenance): formal `seed = NULL` from dina_fit/dina; from forestsearch `dina_args$seed`, else seedit (8316951); stored in object via bootstrap/CV; argument-derived per replicate in the calibration loops. | yes, when seed is non-NULL: on.exit saves the global .Random.seed and assigns it back, or rm()s it if none existed at entry. | dina_fit; dina; dina_subgroup_bootstrap (point refit only, when sg and fit are NULL); forestsearch (subgroup_method = "dina" or use_dina = TRUE) and, through it, bootstrap_results, forestsearch_bootstrap_dofuture, forestsearch_Kfold, forestsearch_tenfold, fpr_calibration, compare_selection_rules, run_null_calibration, run_simulation_analysis, mrct_case_study, mrct_trial_analysis, mrct_region_sims, fs_fdr_report. Not dina_bagged/dina_fit_bagged (bags hard-code `seed = NULL`) or fpr_calibration_transfer. |
| 7 | forestsearch_Kfold | R/forestsearch_cross_validation.R:370 | argument: formal `seedit = 8316951L`; no package code calls forestsearch_Kfold. Runs only when Kfolds < n; the default is leave-one-out, so not by default. | no: the on.exit restores only future::plan; .Random.seed is never saved. | forestsearch_Kfold only (no internal callers). |
| 8 | forestsearch_tenfold (inside the %dofuture% foreach body over ksim) | R/forestsearch_cross_validation.R:959 | argument-derived: `seed + 1000 * ksim`, from the formal `seed = 8316951L` and the foreach iterator ksim. | no, from package code: on.exit restores only future::plan; on the default multisession plan only the worker RNG changes; under sequential, future's own reset restores it. | forestsearch_tenfold only (no internal callers). |
| 9 | fpr_calibration, nested closure .run_batch (inside the %dofuture% body) | R/fpr_calibration.R:302 | argument-derived: `seed_i = sim_seed + seed_offset + i` (formal `sim_seed = 42L`; offset 0L or nsims + 10000L; foreach index i); permutation mode only (null_fn NULL). | no, from package code: nothing saves .Random.seed (on.exit restores only the plan, when n_workers > 1); any RNG reset in the %dofuture% task comes from future. | fpr_calibration (direct). fpr_calibration_transfer calls fpr_calibration but in practice cannot reach this line: it passes a NULL df.analysis, so permutation mode stops (swallowed by tryCatch), unless the user hand-edits `cal_result$call_args`. |
| 10 | compute_frontier_cis | R/frontier_cis.R:164 | argument: formal `seed = NULL` (guarded); from compare_selection_rules its formal `ci_seed = 1L`, passed when compute_cis = TRUE (the default). | no: neither compute_frontier_cis nor compare_selection_rules saves or restores .Random.seed. | compute_frontier_cis (direct); compare_selection_rules (once per sg_focus x selection_rule combination, re-seeding with the same ci_seed). |
| 11 | fs_dgm_feasibility (inner replicate loop) | R/fs_dgm_feasibility.R:181 | argument-derived: `s_i <- as.integer(seed) + i`, from the formal `seed = NULL` and the replicate index i; skipped when seed is NULL. | partly: saves .Random.seed unconditionally and restores it on exit only if it existed at entry; no rm branch, so a newly created one remains. | fs_dgm_feasibility only (no internal callers). |
| 12 | fs_mr_field_uniform (internal, not exported) | R/fs_mr_field_uniform.R:93 | argument-derived: `seed + 910000L` from fs_mr_inference's formal seed (NULL: no set.seed); from forestsearch `mr_inference_args$seed` else seedit, so 9226951 by default; wrappers pass per-replicate or stored seeds. | no: neither fs_mr_field_uniform.R nor fs_mr_inference.R touches .Random.seed or uses on.exit. | Only with field_uniform = TRUE (default FALSE; no wrapper sets it) and ci_method = "field": fs_mr_inference (direct); forestsearch (mr_inference = TRUE); through forestsearch: compare_selection_rules, fpr_calibration, fs_fdr_report, mrct_case_study, mrct_trial_analysis, mrct_region_sims, run_null_calibration, run_simulation_analysis (FS arms need front ends off); bootstrap_results, forestsearch_bootstrap_dofuture, forestsearch_Kfold, forestsearch_tenfold only with mr_in_replicates = TRUE. Not fpr_calibration_transfer. |
| 13 | fs_mr_inference | R/fs_mr_inference.R:615 | argument: formal `seed = NULL` (direct call unseeded); from forestsearch `mr_inference_args$seed` else seedit (8316951); wrappers pass a constant, a per-replicate derived seed or the stored value. | no: fs_mr_inference saves nothing; on reaching paths only mrct_case_study restores, through a wrapper-level on.exit. | fs_mr_inference; forestsearch (mr_inference = TRUE); compare_selection_rules; fpr_calibration; fs_fdr_report; mrct_case_study; mrct_trial_analysis; mrct_region_sims; run_null_calibration; run_simulation_analysis (FS arms need use_lasso, use_grf, use_dina FALSE, or dina-/grf- labels); bootstrap_results, forestsearch_bootstrap_dofuture, forestsearch_Kfold, forestsearch_tenfold only with mr_in_replicates = TRUE. Not fpr_calibration_transfer or the formals-only fs_family_report and `fs_oc_*`. |
| 14 | fs_mr_inference (field block) | R/fs_mr_inference.R:912 | argument-derived: `seed + 900000L` from the same seed as R/fs_mr_inference.R:615 (forestsearch default 9216951); runs only when ci_method = "field". | no | Same list as R/fs_mr_inference.R:615, but only with ci_method = "field" (default for a direct call and forestsearch's consistency branch; the DINA/GRF branches default to "ij"). Not fs_family_report, the four `fs_oc_*` functions or fpr_calibration_transfer. |
| 15 | fs_oc_invert | R/fs_oc_grid.R:311 | argument: formal `seed = NULL` (NULL means no set.seed); no internal caller. | no | fs_oc_invert only (nothing in R/ calls it). |
| 16 | .fs_oc_sweep (internal) | R/fs_oc_grid.R:492 | argument: formal seed (no default), forwarded unchanged from fs_oc_grid's formal `seed = NULL`; NULL means no set.seed. | no | fs_oc_grid only (sole call at fs_oc_grid.R:152). |
| 17 | fs_oc_predict | R/fs_oc_predict.R:234 | argument: formal `seed = NULL` (NULL means no set.seed); no internal caller. | no | fs_oc_predict only. |
| 18 | generate_aft_dgm_flex | R/generate_aft_dgm_main.R:248 | argument: formal `seed = 8316951`, set unconditionally (NULL reseeds randomly); fixed at 8316951 from create_dgm_for_mrct and find_k_inter_for_target_hr. | no | generate_aft_dgm_flex, inject_mrct_structure, mrct_build_dgm, calibrate_k_treat, sensitivity_analysis_k_inter (user-settable); create_dgm_for_mrct, find_k_inter_for_target_hr (fixed 8316951). |
| 19 | generate_glm_dgm | R/generate_glm_dgm.R:402 | argument: formal `seed = 8316951L`, set unconditionally; calibrate_glm_interaction passes its own formal `seed = 8316951L`. | no: the later save (:416-417, only when k_random_noise > 0) follows this set.seed; the caller's stream is replaced. | generate_glm_dgm; calibrate_glm_interaction. |
| 20 | generate_glm_dgm (k_random_noise block) | R/generate_glm_dgm.R:418 | argument: formal `noise_seed = 20260807L`; runs only when k_random_noise > 0 (default 0L); only a direct generate_glm_dgm() call can set either. | yes: explicit save/restore (not on.exit) of the stream after the sample() at :403, so the noise draws leave the downstream stream unchanged. | generate_glm_dgm only (calibrate_glm_interaction can never set k_random_noise > 0). |
| 21 | lasso_selection | R/get_FSdata_helpers.R:259 | fixed: the literal 8316951 (formal seedit default) from get_FSdata, forestsearch and every wrapper, since seedit is dropped; an argument only on a direct lasso_selection(seedit =) call. | yes: the old state is captured and restored, or removed if absent, by an on.exit registered before the set.seed. | lasso_selection (direct, argument); get_FSdata (use_lasso = TRUE); forestsearch (use_lasso = TRUE with mr_inference = FALSE); run_simulation_analysis (legacy FS arms force use_lasso); compare_selection_rules, fpr_calibration, mrct_region_sims, run_null_calibration (via fs args); mrct_case_study, mrct_trial_analysis (only with mr_inference = FALSE); bootstrap_results, forestsearch_bootstrap_dofuture, forestsearch_Kfold, forestsearch_tenfold (stored use_lasso = TRUE, no replicate MR). Not fs_fdr_report. |
| 22 | grf.subg.harm.glm | R/grf_subg_harm_glm.R:460 | argument: formal `seedit = 8316951L`, set unconditionally; from forestsearch its seedit (default 8316951); argument-derived in fpr_calibration and run_null_calibration, stored in object in bootstrap/CV. | no | grf.subg.harm.glm (direct); with a GLM outcome_type only: forestsearch (subgroup_method = "grf", or use_grf = TRUE with mr_inference = FALSE), run_simulation_analysis (run_grf, run_fs_grf or methods), compare_selection_rules, fpr_calibration, run_null_calibration, bootstrap_results, forestsearch_bootstrap_dofuture, forestsearch_Kfold, forestsearch_tenfold, mrct_trial_analysis; mrct_region_sims and mrct_case_study only technically (forcing a GLM outcome_type onto survival data). Not fs_fdr_report or fpr_calibration_transfer. |
| 23 | guohe_adaptive_r | R/guohe_adaptive_r.R:205 | argument: formal `seed = NULL` (guarded; NULL means no set.seed); no function in R/ calls guohe_adaptive_r, so only a direct caller sets it. | no: no save/restore or on.exit; the fold draw and fold-level bootstrap indices consume the stream, leaving the caller's state changed. | guohe_adaptive_r (direct only; no internal callers). |
| 24 | guohe_algorithm3 | R/guohe_algorithm3.R:395 | argument: formal `seed = NULL` (guarded); set directly, from fs_to_guohe(`seed = NULL`), or by guohe_adaptive_r's seed on the refit (fold fits pass NULL). | no: no save/restore or on.exit; the B index vectors are drawn right after, leaving the caller's RNG advanced. | guohe_algorithm3 (direct); fs_to_guohe; guohe_adaptive_r (refit only; its fold fits never run this set.seed). |
| 25 | expand_covariates | R/mrct_inject_dgm.R:97 | argument: formal `seed = 20260903L`, set unconditionally; inject_mrct_structure always passes its own seed (default 8316951L), settable via mrct_build_dgm's `...`. | no: no save/restore or on.exit; the bootstrap and copula branches both leave the caller's RNG advanced. | expand_covariates (direct); inject_mrct_structure; mrct_build_dgm (via `...`). |
| 26 | make_region_model -> nested closure `draw` (R/mrct_inject_dgm.R:276-279), returned as list element $draw | R/mrct_inject_dgm.R:277 | argument-derived: `seed + 1L` and `seed + 2L` from inject_mrct_structure's formal `seed = 8316951L`; a plain argument (default NULL, no seeding) for a direct `$draw()` caller. | no: neither the closure nor inject_mrct_structure saves or restores the state around the two draws. | inject_mrct_structure; mrct_build_dgm (via `...`); make_region_model only if the user calls the returned `$draw` closure (make_region_model itself never runs it). Not reached in observed-region mode. |
| 27 | inject_mrct_structure | R/mrct_inject_dgm.R:572 | argument-derived: `seed + 3L` from the formal `seed = 8316951L`, set unconditionally; settable directly or via mrct_build_dgm's `...`. | no: no save/restore or on.exit; calculate_hazard_ratios draws rexp from this stream, leaving the caller's state changed. | inject_mrct_structure (direct); mrct_build_dgm. |
| 28 | mrct_truth_by_region | R/mrct_inject_dgm.R:664 | argument: formal `seed = 1L`, set unconditionally; no internal caller. | no: no save/restore or on.exit; the caller's RNG is left at set.seed(seed) plus nrow(ds) rexp draws. | mrct_truth_by_region (direct only). |
| 29 | run_simulation_analysis | R/run_simulation_analysis.R:509 | argument-derived: `seed_base + 1000L * sim_id` (formals `seed_base = 8316951L` and sim_id); runs only when n_add_noise > 0 (default 0L). | no: no save/restore or on.exit; the rnorm noise draws leave the caller's RNG advanced. | run_simulation_analysis (direct only). |
| 30 | .run_subgroup_sims_one (internal; called inside the foreach %dofuture% body of run_subgroup_sims) | R/run_subgroup_sims.R:512 | argument-derived: `seed_base + ss + rand_seed_offset` (run_subgroup_sims formals `seed_base = 0L`, `rand_seed_offset = 1e6L`; foreach index ss); runs by default (benchmarks non-NULL). | no (package code): on.exit restores only the plan and future.rng.onMisuse; under plan(sequential), future itself restored the caller's .Random.seed (checked, future 1.80.0). | run_subgroup_sims (only caller). |
| 31 | setup_gbsg_dgm | R/setup_gbsg_dgm.R:130 | argument: formal `noise_seed = 20260807L`; runs only when k_random_noise > 0 (default 0L), so never from fs_fdr_report (which would fix it at 20260807L anyway). | yes, by explicit save/restore (not on.exit); the saved state already follows .create_gbsg_dgm_'s unrestored set.seed calls, not the caller's pre-call state. | setup_gbsg_dgm (direct, only when k_random_noise > 0). fs_fdr_report is in the static list but never runs this set.seed. |
| 32 | .create_gbsg_dgm_ (internal) | R/sim_aft_gbsg.R:478 | argument: formal `seed = SEED_BASE` (8316951L), settable via setup_gbsg_dgm, create_gbsg_dgm, calibrate_k_inter (`...`) or validate_k_inter_effect (`...`); fixed at 8316951L from fs_fdr_report. | no: neither .create_gbsg_dgm_ nor any wrapper saves or restores it; setup_gbsg_dgm's save covers only its later noise draw. | setup_gbsg_dgm; create_gbsg_dgm (deprecated); calibrate_k_inter (via `...`); validate_k_inter_effect (via `...`); fs_fdr_report (fixed 8316951L). |
| 33 | .create_gbsg_dgm_ (internal) | R/sim_aft_gbsg.R:604 | argument-derived: `seed + 1` from the same formal `seed = SEED_BASE` (8316951L), via the same four wrappers; fixed at 8316952 from fs_fdr_report. | no: no save/restore or on.exit; the rexp draw leaves the caller's RNG on the seed + 1 stream. | setup_gbsg_dgm; create_gbsg_dgm (deprecated); calibrate_k_inter (via `...`); validate_k_inter_effect (via `...`); fs_fdr_report (fixed 8316952). |
| 34 | .simulate_from_gbsg_dgm_ (n non-NULL, draw_treatment = TRUE branch) | R/sim_aft_gbsg.R:856 | stored in object (from an argument): `dgm$seed + sim_id`; `dgm$seed` from the create_gbsg_dgm or setup_gbsg_dgm seed formal (default 8316951L), sim_id from simulate_from_gbsg_dgm (default 1). | no: R/sim_aft_gbsg.R has no on.exit or .Random.seed handling, so the caller's RNG stream is overwritten. | simulate_from_gbsg_dgm (deprecated wrapper, only caller), with a dgm from create_gbsg_dgm or setup_gbsg_dgm. |
| 35 | .simulate_from_gbsg_dgm_ (n non-NULL, draw_treatment = FALSE branch) | R/sim_aft_gbsg.R:884 | stored in object (from an argument): `dgm$seed + sim_id`, as at :856 (`dgm$seed` default 8316951L; sim_id default 1). | no | simulate_from_gbsg_dgm. |
| 36 | .simulate_from_gbsg_dgm_ (event-time draw) | R/sim_aft_gbsg.R:900 | stored in object (from an argument): `dgm$seed + sim_id + 1000` (literal offset); `dgm$seed` from the create_gbsg_dgm/setup_gbsg_dgm seed formal, sim_id from simulate_from_gbsg_dgm. | no | simulate_from_gbsg_dgm. |
| 37 | .simulate_from_gbsg_dgm_ (censoring draw) | R/sim_aft_gbsg.R:908 | stored in object (from an argument): `dgm$seed + sim_id + 2000` (literal offset); same sources as :900. | no | simulate_from_gbsg_dgm. |
| 38 | simulate_from_dgm | R/simulate_from_dgm.R:208 | argument: formal `seed = NULL` (guarded); passed verbatim from calibrate_cens_adjust and fs_build_eval_frame, argument-derived (seed or seed_base plus replicate index) from every other wrapper. | no in simulate_from_dgm; among callers only mrct_case_study saves the state and restores RNGkind and .Random.seed on exit. | simulate_from_dgm (direct); calibrate_cens_adjust; fs_build_eval_frame; fs_fdr_report; mrct_case_study; mrct_region_sims; run_simulation_analysis; run_subgroup_sims. |
| 39 | simulate_from_glm_dgm | R/simulate_from_glm_dgm.R:97 | argument: formal `seed = NULL` (guarded); verbatim from fs_build_eval_frame (`eval_seed = 20260628L`); argument-derived (seed + index) from fs_dgm_feasibility, run_simulation_analysis and run_subgroup_sims. | no in simulate_from_glm_dgm; fs_dgm_feasibility restores partly (only an existing .Random.seed); other callers do not restore. | simulate_from_glm_dgm (direct); fs_build_eval_frame (binary outcome only); fs_dgm_feasibility; run_simulation_analysis; run_subgroup_sims. |
| 40 | .make_candidate_rng_seeds | R/subgroup_consistency_main.R:58 | argument (formal seed, no default, from subgroup.consistency's `seed = 8316951`) only on a direct subgroup.consistency(seed =) call; fixed at 8316951 from forestsearch() and every wrapper. | yes: on.exit restores RNGkind (switched at :48) and the saved .Random.seed, or rm()s it if none existed. | Parallel batch branch only (sg_focus not "maxeff"): forestsearch (default multisession); compare_selection_rules; fpr_calibration and run_simulation_analysis only with non-sequential parallel_args in fs_params; subgroup.consistency direct only with non-sequential parallel_args. Never from bootstrap_results, forestsearch_bootstrap_dofuture, forestsearch_Kfold, forestsearch_tenfold, run_null_calibration, fs_fdr_report, mrct_region_sims, mrct_trial_analysis, mrct_case_study or fpr_calibration_transfer (sequential). |
| 41 | subgroup.consistency | R/subgroup_consistency_main.R:498 | argument (formal `seed = 8316951`) only on a direct subgroup.consistency(seed =) call; fixed at 8316951 from forestsearch() and every wrapper, whatever seedit or the replicate seed is. | no: no .Random.seed or on.exit handling; overwrites the caller's RNG whenever seed is non-NULL, even where this stream is unused. | subgroup.consistency (direct: argument); forestsearch (fixed 8316951) and the 13 exported functions that reach it through forestsearch() (fixed): bootstrap_results, forestsearch_bootstrap_dofuture, compare_selection_rules, forestsearch_Kfold, forestsearch_tenfold, fpr_calibration, fpr_calibration_transfer, run_null_calibration, fs_fdr_report, mrct_region_sims, mrct_trial_analysis, mrct_case_study, run_simulation_analysis. |
| 42 | generate_bootstrap_synthetic | R/synthetic_data_perturbation.R:98 | argument: formal `seed = 123`, used verbatim (`seed = NULL` still calls set.seed(NULL), reseeding randomly); also from generate_gbsg_bootstrap_general(`seed = 123`). | no | generate_bootstrap_synthetic (direct); generate_gbsg_bootstrap_general. |
| 43 | generate_bootstrap_with_noise | R/synthetic_data_perturbation.R:506 | argument: formal `seed = 123`, checked as numeric (NULL rejected, so seeding cannot be turned off) and used verbatim. | no | generate_bootstrap_with_noise only (no internal callers). |
| 44 | .detect_prob_glm_single (method = "monte_carlo" branch) | R/truefind_asymptotic_glm.R:338 | argument (compute_detection_probability_glm's formal `seed = 42L`) only on a direct call; fixed at 42L from predict_fpr_corrected, run_null_calibration, fpr_calibration, fpr_calibration_transfer and suggest_thresholds. | no: neither the helper, compute_detection_probability_glm nor these callers saves or restores .Random.seed. | compute_detection_probability_glm (direct: argument; only with method = "monte_carlo", the default "cubature" never seeds); predict_fpr_corrected, run_null_calibration, fpr_calibration, fpr_calibration_transfer, suggest_thresholds (fixed 42L; all force monte_carlo). |

**Scope (Step 2).**
- Only row 41 is the fixed seed this task changes, and that line is not edited (§3.2).
- The seeds that come from an argument are listed above and not changed. Exported functions outside the 15 below that set such a seed still leave the caller's state as before, e.g. `simulate_from_dgm(seed = )` and `generate_aft_dgm_flex(seed = 8316951)`.
- Inside the 15 functions, the argument seeds still run unchanged; the outermost exit restores over them.

### 3.2 The change (Rule S)

**The mechanism.** All of it is in `R/subgroup_consistency_main.R`, next to `.make_candidate_rng_seeds()`, which already saves and restores the state the same way. The fixed-seed line itself is not edited.

- **`.fs_rng`** (`:82`): a package-level environment holding `active = FALSE`.
- **`.fs_preserve_rng()`** (`:92-112`), when `.fs_rng$active` is FALSE (the outermost exported call in this R process):
  - it sets the flag;
  - it records whether `.Random.seed` exists, its value and `RNGkind()`;
  - it registers, with `on.exit(add = TRUE)` in the calling function's frame, a restore. The restore clears the flag and puts `.Random.seed` back, or, if it was absent, restores the RNG kind and removes `.Random.seed` again.
  - The restore runs on normal exit and on error.
  - When the flag is already TRUE (a nested exported call), it does nothing.
- **`.fs_rng_hold()`** (`:122-125`) sets the flag and never clears it (the multisession finding below).
- **One line, `.fs_preserve_rng()`, as the first statement of each of the 15 exported functions that can call `subgroup.consistency()`**, plus one `.fs_rng_hold()`:

| Function | Line at the branch head | Call chain to `subgroup.consistency()` |
|---|---|---|
| `forestsearch` | `R/forestsearch_main.R:1725` | direct (`do.call(subgroup.consistency, …)`) |
| `subgroup.consistency` | `R/subgroup_consistency_main.R:468` | itself (a direct call) |
| `forestsearch_bootstrap_dofuture` | `R/bootstrap_dofuture_main.R:220` | → `bootstrap_results` → `forestsearch` |
| `bootstrap_results` | `R/bootstrap_analysis_dofuture.R:307` | → `forestsearch` |
| `forestsearch_Kfold` | `R/forestsearch_cross_validation.R:298` | → `forestsearch` |
| `forestsearch_tenfold` | `R/forestsearch_cross_validation.R:815`; `.fs_rng_hold()` at `:1000` | → `forestsearch` |
| `run_simulation_analysis` | `R/run_simulation_analysis.R:420` | → `.run_fs_analysis_gen` → `forestsearch` |
| `mrct_region_sims` | `R/mrct_simulation.R:283` | → `forestsearch` |
| `mrct_case_study` | `R/mrct_case_study.R:250` | → `.mrct_analyze_trial` → `forestsearch` |
| `mrct_trial_analysis` | `R/mrct_trial_analysis.R:92` | → `.mrct_analyze_trial` → `forestsearch` |
| `compare_selection_rules` | `R/compare_selection_rules.R:494` | → `forestsearch` |
| `fpr_calibration` | `R/fpr_calibration.R:207` | → `forestsearch` |
| `fpr_calibration_transfer` | `R/fpr_calibration.R:504` | → `fpr_calibration` → `forestsearch` |
| `run_null_calibration` | `R/calibrate_null_correction.R:315` | → `forestsearch` |
| `fs_fdr_report` | `R/fs_fdr_report.R:195` | → `.fs_fdr_run_null` → `.fs_fdr_one_rep` → `forestsearch` |

**How the 15 were found.**
- A static call graph over the namespace (scratch `call_graph2.R`) counts an edge only for a call: `g(...)`, `forestsearch::g(...)`, or `g` passed as the function argument of `do.call()`, `lapply()` and the like.
- The 15 are exactly the exported functions with a call chain to `subgroup.consistency()`. The chain above is the shortest witness.
- No exported function reaches it only through a value reference.
- The verified Part 1.1 classification found the same `forestsearch()` call sites by reading (`bootstrap_analysis_dofuture.R:614`, `compare_selection_rules.R:561`, `calibrate_null_correction.R:358`, `fs_fdr_report.R:364`, `fpr_calibration.R:308`, `mrct_trial_helpers.R:47`, `forestsearch_cross_validation.R:480, :1002`, `run_simulation_analysis.R:819`, `mrct_simulation.R:575`).

**Why each exported caller, and not just `subgroup.consistency()`.**
- Rule S needs the caller's state as it was *on entry*. Several of these functions reseed before they reach `subgroup.consistency()`, from their own arguments (the bootstrap's `set.seed(seed)`, K-fold `set.seed(seedit)`, `run_simulation_analysis()`'s `set.seed(seed_base + 1000L * sim_id)`). So the state has to be captured at the entry of the outermost exported call.
- A guard in `subgroup.consistency()` alone would restore part-way through `forestsearch()`, so the `mr_inference` draws after it would come from a different stream.
- The flag makes every nested exported call a no-op: `forestsearch()` inside a bootstrap or CV loop, `subgroup.consistency()` inside `forestsearch()`.
- The argument seeds inside these functions are not edited. They are restored over, at the outermost exit, like everything else drawn inside the call.

**Multisession finding and its correction (`.fs_rng_hold()`).**
- **The question.** An independent read-only audit of the 12 future bodies in `R/` asked whether a worker draws random numbers after a nested guarded call returns, within the same future. A worker process does not share the flag, so there the first exported call is "outermost" and restores the worker's state.
- **The one case.** `forestsearch_tenfold()` (`R/forestsearch_cross_validation.R:953`) runs K `forestsearch()` calls one after another in each iteration. With a GLM outcome, `ps_method = "lasso"` and `ps_adjust_method` `"iptw"` or `"dr_gcomp"`, fold k ≥ 2 draws `cv.glmnet`'s fold ids (`R/glm_effect_estimators.R:1496`, no `foldid`) before its first `set.seed()`. Under multisession, folds 2..K would start from fold 1's entry state instead of its exit state, so the package's own values would change.
- **The other 11 bodies.** The bootstrap, K-fold, `fpr_calibration()`, `fs_fdr_report()` and `mrct_region_sims()` bodies make one guarded call per iteration and draw nothing after it; the other six bodies call no guarded function. doFuture re-seeds every element (`seed = TRUE`).
- **Measured.** Trees were installed into temporary libraries under `/tmp/seed_check_2026-11-11/`, and tenfold was run with 2 multisession workers, 2 sims and 4 folds, on continuous data with a lasso propensity score and iptw.
  - Guards without the hold: `fold_summary$training_fs_hr` differs from the start tree.
  - Guards with the hold (the final tree): identical, after dropping timing fields.
- **The correction.** One call, `.fs_rng_hold()`, in the iteration body before its fold loop (`:1000`). It sets the flag inside the worker, so the fold calls continue one stream, as they do in the main process.
  - It has to be a function call: a namespace environment named in the body would reach the worker as a copy.
  - In the main process the flag is already TRUE, so the call changes nothing.

**How the design got here (in this run).**
1. The first version guarded 23 functions, taken from a static graph that also counted `formals(forestsearch)` and `formals(subgroup.consistency)` reads as edges.
2. The verified Part 1.1 classification found that 8 of them never call `subgroup.consistency()`: `fs_declaration_calibration`, `fs_declaration_table`, `fs_family_report`, `fs_mr_inference`, `fs_oc_family_enumerate`, `fs_oc_grid`, `fs_oc_invert` and `fs_oc_predict`. Their guards were removed, so these functions behave exactly as before, their own argument seeds included.
3. The multisession audit led to `.fs_rng_hold()`.
4. A suite and `check()` run that had been started on the 23-guard version was stopped before it finished, and its results were not used. §3.5 is the single complete run on the final code.

**Also changed.** `DESCRIPTION` 0.4.1.9007 → 0.4.1.9008. `NEWS.md` has one entry in "Behaviour changes", after the part-1 entries (§3.6), and the heading now reads 0.4.1.9008.

### 3.2a Diff (`55e59792`, `git show --stat`)

```
 DESCRIPTION                              |   2 +-
 NEWS.md                                  |  21 ++++-
 R/bootstrap_analysis_dofuture.R          |   1 +
 R/bootstrap_dofuture_main.R              |   1 +
 R/calibrate_null_correction.R            |   1 +
 R/compare_selection_rules.R              |   1 +
 R/forestsearch_cross_validation.R        |   7 ++
 R/forestsearch_main.R                    |   1 +
 R/fpr_calibration.R                      |   2 +
 R/fs_fdr_report.R                        |   1 +
 R/mrct_case_study.R                      |   1 +
 R/mrct_simulation.R                      |   1 +
 R/mrct_trial_analysis.R                  |   1 +
 R/run_simulation_analysis.R              |   1 +
 R/subgroup_consistency_main.R            |  60 ++++++++++++
 tests/testthat/test-rng-state-restored.R | 152 +++++++++++++++++++++++++++++++
 16 files changed, 252 insertions(+), 2 deletions(-)
```

The helper block added to `R/subgroup_consistency_main.R` (after `.make_candidate_rng_seeds()`, lines 69-127 at the branch head):

```r
# ---------------------------------------------------------------------------
# The caller's random-number state across exported calls.
#
# subgroup.consistency() seeds the generator (set.seed(seed); seed = 8316951
# unless supplied, and forestsearch() does not supply it), and other exported
# functions seed it from their own arguments.  Every exported function that
# can reach that call starts with .fs_preserve_rng(): the OUTERMOST such call
# saves the caller's .Random.seed (or its absence) and puts it back when it
# exits, normally or by an error.  A nested exported call (forestsearch()
# inside a bootstrap or cross-validation loop, subgroup.consistency() inside
# forestsearch()) does nothing, so every draw made inside the outer call comes
# from the same stream as before and every computed value is unchanged.
# ---------------------------------------------------------------------------
.fs_rng <- new.env(parent = emptyenv())
.fs_rng$active <- FALSE

#' Restore the caller's RNG state when the outermost exported call exits
#'
#' @param envir The calling function's frame; the restore is registered there
#'   with \code{on.exit(add = TRUE)}.
#' @return Invisibly, \code{TRUE} when this call registered the restore,
#'   \code{FALSE} when an enclosing exported call already had.
#' @noRd
.fs_preserve_rng <- function(envir = parent.frame()) {
  if (isTRUE(.fs_rng$active)) return(invisible(FALSE))
  .fs_rng$active <- TRUE
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_kind <- RNGkind()
  restore <- function() {
    .fs_rng$active <- FALSE
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else {
      suppressWarnings(do.call(RNGkind, as.list(old_kind)))
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        rm(".Random.seed", envir = globalenv())
      }
    }
  }
  do.call(base::on.exit, list(as.call(list(restore)), TRUE, FALSE),
          envir = envir)
  invisible(TRUE)
}

#' Mark the rest of a future evaluation as inside an exported call
#'
#' A future worker does not share \code{.fs_rng$active} with the main
#' process, so the first exported call made there would restore the worker's
#' generator when it returns.  A future body that makes several exported calls
#' in one evaluation calls this first, so those calls leave the stream alone
#' between them, as they do in the main process.  It never clears the flag.
#' @noRd
.fs_rng_hold <- function() {
  .fs_rng$active <- TRUE
  invisible(NULL)
}


```

Every other hunk under `R/` is one added line, `  .fs_preserve_rng()`, as the first statement of the function (the table in §3.2), plus the hold in `forestsearch_tenfold()`:

```r
    # The K forestsearch() calls below run one after another in this
    # evaluation; in a worker they must not restore the generator between
    # folds (see .fs_rng_hold()).
    .fs_rng_hold()

```

`git diff -U0 3fbbc55c 55e59792 -- R/` has 79 added lines and no removed line. Every file ending is unchanged: one stray final newline that the insertion script added to `R/run_simulation_analysis.R` was removed before the commit.

### 3.3 Equivalence check (Step 4) — **PASS**

**Setup.**
- Scratch directory: `/tmp/seed_check_2026-11-11/`.
  - The start tree is `git archive 3fbbc55c` → `start/`.
  - The branch tree is the final working tree (`git ls-files -z -co --exclude-standard`) → `branch_final/`; its `R/` was `diff -r`-identical to the repository's.
- Each tree is loaded with `pkgload::load_all()` in its own `Rscript --vanilla` process (`equiv_run.R`).
- Plan: `future::plan("sequential")`, and `parallel_args = list(plan = "sequential", workers = 1)` for every call.
- Data come from the package's own test helpers, built before any `set.seed(123)` (the helpers seed themselves):
  - survival, N = 600 (`.make_survival_data(seed = 42)`; confounders age, stage, sex, noise);
  - a `df.test` of N = 300 (seed 43);
  - binary, N = 600 (`.make_binary_data(seed = 44)`).
- `set.seed(123)` immediately before each call. `.Random.seed` was recorded before and after.
- **Comparison** (`equiv_compare.R`):
  - every differing leaf is listed. A leaf whose name is a timing field is dropped: `time_search`, `minutes_all`, `timing_seconds`, `timing_minutes`, and the per-replicate columns `tmins_search` and `tmins_iteration`. A survival `time` column is data and is kept.
  - data.tables are compared as data frames, because their internal self-reference pointer is not data. Then `identical()` is applied.

| # | Call | `identical()` after dropping timing fields | Non-timing leaves that differ | Timing leaves dropped | Branch: `.Random.seed` after = before | Start: after = before | Wall, start / branch |
|---|---|---|---|---|---|---|---|
| 1 | `forestsearch()`, survival, defaults (sequential) | TRUE | 0 | `$find.grps$time_search`, `$minutes_all` | **TRUE** | FALSE | 0.6 / 0.5 s |
| 2 | same, `mr_inference = TRUE`, no seed in `mr_inference_args` | TRUE | 0 | the two above + `$mr_inference$timing_seconds`, `$mr_inference$field$timing_seconds`, `$mr_inference$field$complement$timing_seconds` | **TRUE** | FALSE | 7.8 / 7.8 s |
| 3 | same, with `df.test` | TRUE | 0 | 2 (as row 1) | **TRUE** | FALSE | 0.2 / 0.2 s |
| 4 | `forestsearch_bootstrap_dofuture()` on the fit of 1, `nb_boots = 10` | TRUE | 0 | `$results$tmins_search`, `$results$tmins_iteration` | **TRUE** | FALSE | 5.5 / 5.5 s |
| 5 | `forestsearch_Kfold()` on the fit of 1, `Kfolds = 5` | TRUE | 0 | `$timing_minutes` | **TRUE** | FALSE | 1.3 / 1.3 s |
| 6 | `forestsearch()`, binary outcome, defaults | TRUE | 0 | 2 (as row 1) | **TRUE** | FALSE | 0.6 / 0.6 s |

- **Gate: PASS.** All six are identical, and in the branch tree `.Random.seed` after each call equals `.Random.seed` before it. In the start tree it never did.
- **Wall clock.** This run took 30 s for both trees together (measured).
  - The equivalence script ran three times as the design changed: 23 guards (30 s), 23 guards with the hold (29 s), and the final 15 guards (30 s). Every run passed the gate.
  - The multisession tenfold runs took 4.5-4.7 s each (four runs; measured in R).
  - The temporary installs were not timed separately.
  - All measured Step 4 runs sum to under 2 minutes (cap 15 min).
- **Outside the gate.** The final tree also reproduced, under a real multisession plan, the tenfold result of the start tree (§3.2).

### 3.4 The consequence (Step 5, recorded, not gated)

`set.seed(123); fs <- forestsearch(...); bs <- forestsearch_bootstrap_dofuture(fs, nb_boots = 10)`, with no seed between, run twice (A, B) in each tree:

| Question | Result |
|---|---|
| Does the bootstrap output differ between the trees? | **No.** It is identical after dropping timing fields (0 differing leaves). |
| Does it repeat exactly when the three lines are run twice, in the start tree? | **Yes** (A = B). |
| … in the branch tree? | **Yes** (A = B). |

- **Why the bootstrap is unaffected.** `forestsearch_bootstrap_dofuture()` seeds itself. Its formal is `seed = 8316951L`. The row indices come from `.make_boot_index_matrix(seed = boot_seed)` (`R/bootstrap_dofuture_main.R:390-394` at the branch head), and `bootstrap_results()` calls `set.seed(seed)` before its futures derive their streams (`R/bootstrap_analysis_dofuture.R:335`). So the bootstrap never read the state `forestsearch()` left behind.
- **What does change.** An unseeded draw straight after `forestsearch()` (`after_draw.R`: `set.seed(123); forestsearch(); runif(3)`, twice per tree):
  - start tree: `0.35640416 0.17546534 0.01537917`, which continues the stream `subgroup.consistency()` seeded with 8316951; repeats exactly;
  - branch tree: `0.2875775 0.7883051 0.4089769`, the same as `runif(3)` right after `set.seed(123)`; repeats exactly.
- **NEWS wording follows the measurement.** The task's NEWS example ("a bootstrap run straight after `forestsearch()`") does not hold for this package. The entry uses `sample()` and `simulate_from_dgm(seed = NULL)` as examples instead and says the bootstrap is not affected (§3.6, OPEN ITEM 2).

### 3.4a Tests (Step 6)

`tests/testthat/test-rng-state-restored.R` (new): 7 blocks, 16 expectations, all under `skip_on_cran()` (file level), sequential, N = 400 survival, **6.4 s**. The blocks:
1. `forestsearch()` leaves `.Random.seed` as it found it, and the caller's next draw is the one after `set.seed(123)`.
2. The same with `mr_inference = TRUE`.
3. An absent `.Random.seed` is absent again.
4. A direct `subgroup.consistency()` call (arguments captured from a `forestsearch()` run) restores the state, and its result does not depend on the caller's state.
5. Two consecutive identical calls give identical results.
6. The restore runs on an error exit (the mechanism itself), and later calls restore again.
7. `.fs_rng_hold()` keeps a later exported call from restoring.

**At the start SHA.** Run against a start-SHA copy (blocks 1-6, before block 7 existed), the restore expectations failed (1, 2, 3, 4) and block 6 errored (no `.fs_preserve_rng`). That shows the file detects the old behaviour.

### 3.5 Suite and check (Step 8)

The run went once, on the final code (the working tree that became `55e59792`), by the check route of `REPORT_cran_prep_2026-11-11.md` §6.

| Run | Baseline (`REPORT_cran_prep_2026-11-11.md` §6, not re-run) | This branch | Result |
|---|---|---|---|
| Full suite (`devtools::test()`, `NOT_CRAN=true`, in the repository) | 641 blocks, 87 files, 8,455 expectations, 0 failures, 3 skips, 38 warnings; 695 s | **648 blocks, 88 files, 8,471 expectations, 0 failures, 0 errors, 3 skips, 38 warnings**; wall **686 s** (22:04:47–22:16:13; cap 30 min) | **PASS** |
| Difference | — | +7 blocks, +1 file, +16 expectations: exactly `test-rng-state-restored.R`. The skips are the same 3 (`test-cross-outcome-parity.R` ×2, `test-cv-no-subgroup-edges.R`). Warnings are unchanged. | |
| `devtools::check(document = FALSE, error_on = "never")` | 0 errors / 0 warnings / 0 notes; 1,067 s | **0 errors / 0 warnings / 0 notes**, Status OK; wall **1,121 s** (22:04:47–22:23:28; examples 36 s, `--run-donttest` 63 s, vignettes 73 s; cap 60 min) | **PASS** |

- **Check route.** An export of the working tree in the session scratchpad: `git ls-files -z -co --exclude-standard`, 9,678 files listed and 9,678 copied. The 17 files that differ from the start SHA were `cmp`-identical to the repository's.
- **Correction to Report 1 §6.** That export used `git ls-files … | tar -T -` without `-z`, so the two tracked files whose names need quoting were not copied:
  - `dev/glm/articles/An introduction to \`grf\` • grf.html`;
  - a PDF under `quarto/methodology/articles/`.
  Both are excluded by `.Rbuildignore` (`^dev$`, `^quarto$`), so that check was not affected. Its "9,660 files" was the listed count; 9,658 were copied. This run uses `-z`.
- **Contention.** The suite and `check()` ran at the same time, as in Report 1.
- **Pre-set decision.** No expectation failed, so none pinned a value of a call made after `forestsearch()` without a seed, and nothing needed updating.
- **Superseded run.** A suite and `check()` run started earlier on the 23-guard version was stopped before it finished, and its partial results were not used.

### 3.6 `NEWS.md` entry (Step 7) and commit (Step 9)

`DESCRIPTION` 0.4.1.9007 → 0.4.1.9008; `NEWS.md` heading → `# forestsearch 0.4.1.9008 (development)`. The entry, under "Behaviour changes":

> * `forestsearch()` and the other exported functions that can reach the fixed
>   seed in `subgroup.consistency()` (`set.seed(8316951)`, which `forestsearch()`
>   does not expose) no longer leave the caller's random-number state changed.
>   These are `subgroup.consistency()` itself and the functions that call
>   `forestsearch()`: `forestsearch_bootstrap_dofuture()`, `bootstrap_results()`,
>   `forestsearch_Kfold()`, `forestsearch_tenfold()`, `run_simulation_analysis()`,
>   `mrct_region_sims()`, `mrct_case_study()`, `mrct_trial_analysis()`,
>   `compare_selection_rules()`, `fpr_calibration()`, `fpr_calibration_transfer()`,
>   `run_null_calibration()` and `fs_fdr_report()`. The outermost such call saves
>   `.Random.seed` on entry and restores it on exit (absent again if it was
>   absent), whatever seeds are set inside it. Their own results are unchanged:
>   every draw inside a call comes from the same stream as before. A later call
>   that draws random numbers without its own seed, straight after
>   `forestsearch()` (for example `sample()`, or `simulate_from_dgm()` with
>   `seed = NULL`), used to continue the stream `forestsearch()` had seeded and
>   now depends on the caller's seed, so set one. `forestsearch_bootstrap_dofuture()`
>   seeds itself (`seed = 8316951L`) and is not affected. Functions that only set a
>   seed passed to them (for example `simulate_from_dgm(seed = )`) still set it.

Committed by explicit path as `55e59792` (`fix: forestsearch functions restore the caller's random-number state`): `DESCRIPTION`, `NEWS.md`, the 13 `R/` files and `tests/testthat/test-rng-state-restored.R`.

## 4. Part 2 — Random draws after a call that reaches the fixed seed (vignettes and articles)

**Method.** For each of the 13 tracked vignette and article sources, one agent traced which evaluated calls reach `set.seed(8316951)` in `subgroup.consistency()` (`R/subgroup_consistency_main.R:498`, formal `seed = 8316951` at `:397`, reached from `forestsearch()` through `do.call(subgroup.consistency, sc_filtered_args)` at `R/forestsearch_main.R:3669`) and which later draws follow with no seed in between. An independent verifier then re-derived each result from the source at `c40d5309` (identical to the start SHA `3fbbc55c` in `R/`, `vignettes/` and `quarto/`); several verifiers also ran short `Rscript --vanilla` checks of how the installed future, doFuture and gt use the RNG. Only the verified results are used below. Nothing was changed for this part. Full per-file detail is in `dev/reports/REPORT_cran_prep_seed_legacy_2026-11-11_detail.md`.

The agents worked from a provisional list of 23 names that may reach the fixed seed. The final change (Part 1) wraps 15 of them. The other 8 (`fs_declaration_calibration`, `fs_declaration_table`, `fs_family_report`, `fs_mr_inference` and the `fs_oc_*` functions) reach the seed only through `formals()` lookups and are not wrapped. **No row and no zero-row verdict below depends on any of those 8 being wrapped.** The verifiers found that they draw nothing in these documents (`fs_declaration_calibration` at `vignettes/forestsearch.qmd:558-559`; `fs_declaration_table` at `actg175-glm-outcomes.qmd:346`). For `fs_mr_inference` in `guohe-nested-continuum.qmd`, the verifier worked through the case where it is left unwrapped and found no change. All 7 rows depend only on `forestsearch`, `forestsearch_bootstrap_dofuture` and `forestsearch_tenfold`, which are among the 15 wrapped functions.

### Per-file summary

| File | Built by CRAN or pkgdown only | Evaluated calls reaching the fixed seed (count; file:line) | Rows |
|---|---|---|---|
| `vignettes/forestsearch.qmd` | CRAN-built vignette, and also the pkgdown Introduction article (`_pkgdown.yml:32,40`) | 1: `:228` `forestsearch()` | 0 |
| `vignettes/articles/actg175-glm-outcomes.qmd` | pkgdown only | 4: `:155`, `:183`, `:421` `forestsearch()`; `:522` `forestsearch_bootstrap_dofuture()`, in its FS iteration, inside the workers | 0 |
| `vignettes/articles/biomarker_effects.qmd` | pkgdown only | 0 | 0 |
| `vignettes/articles/causal_effects_brief_review.qmd` | pkgdown only | 0 | 0 |
| `vignettes/articles/extreme_subgroups.qmd` | pkgdown only | 0 | 0 |
| `vignettes/articles/forestsearch-2024.Rmd` | pkgdown only | 3: `:301` `forestsearch()`, in the main process; `:478` `forestsearch_bootstrap_dofuture()` and `:606` `forestsearch_tenfold()`, both inside the workers. `:644` `forestsearch_Kfold()` is not evaluated (`eval = FALSE` at `:640`). | 5 |
| `vignettes/articles/guohe-nested-continuum.qmd` | pkgdown only | 0 | 0 |
| `vignettes/articles/methodology.qmd` | pkgdown only | 0 | 0 |
| `vignettes/articles/mrct-consistency-simulations.qmd` | pkgdown only | 3: `:153` `mrct_trial_analysis()`, in the main session; `:730` `mrct_region_sims()`, only inside futures, with no leak to the caller; `:779` `mrct_case_study()`, in the main session, which already restores the caller's RNG on exit (`R/mrct_case_study.R:262-271`) | 0 |
| `vignettes/articles/mrct-consistency.qmd` | pkgdown only | 2: `:510` `forestsearch()`; `:789` `mrct_trial_analysis()` | 0 |
| `vignettes/articles/paper_simulations.qmd` | pkgdown only | 2: `:925` and `:983` `run_simulation_analysis()`, inside `%dofuture%` bodies. They run in workers or in sequential futures whose RNG state is restored, so they never set the rendering session's RNG. | 0 |
| `vignettes/articles/survival-post-selection.qmd` | pkgdown only | 2: `:144` `forestsearch()`, in the main session; `:272` `forestsearch_bootstrap_dofuture()`, in its FS iteration, inside the workers (or in the main session under the sequential fallback). `:165` (DINA) and `:183` (GRF) return before `:3669`. | 2 |
| `vignettes/articles/treatment_effect_definitions.qmd` | pkgdown only | 0 | 0 |

"pkgdown only" rests on `.Rbuildignore:20` `^vignettes/articles$`, which keeps `vignettes/articles/` out of the tarball even though `DESCRIPTION:30` declares `VignetteBuilder: knitr, quarto`. Each article is listed in `_pkgdown.yml`.

### Rows

| File:line | Later call | Seeded? | Displayed outputs that would vary from one render to the next | Built by CRAN, or pkgdown only |
|---|---|---|---|---|
| `forestsearch-2024.Rmd:432` | knitr auto-print of `res_tabs$sg10_out`, a gt table from `sg_tables` (maxk == 2 branch, `R/summary_utility_functions.R:1196-1197`). It comes after `forestsearch()` at `:301`. gt 1.4.0 `as.tags.gt_tbl` calls `random_id()` = `sample(letters, 10, replace = TRUE)`, because no table id is set. | No. gt does not seed, the file has no `set.seed`, and `sg_tables` (`:415`) neither seeds nor draws. Only pkgdown's `set.seed(1014)` before the render fixes the stream, and that comes from the build tool. | Chunk fs-results (`:412-433`), the 'Identified Subgroups' table: only its hidden 10-letter HTML container id and scoped CSS; no visible cell changes. After the change the id varies on every plain `rmarkdown::render` or Knit; under pkgdown it changes once and then stays stable. Today it is stable under both. | pkgdown only (`.Rbuildignore:20`; `_pkgdown.yml:67-68`, `:611-613`) |
| `forestsearch-2024.Rmd:439` | knitr auto-print of `res_tabs$tab_estimates`, a gt table built at `R/summary_utility_functions.R:1045`. It comes after `forestsearch()` at `:301`. Same chain: `as.tags.gt_tbl` -> `random_id()` -> `sample()`. | No. There is no `set.seed` or seed argument between `:301` and `:439`. The only call in between is the gt print at `:432`, which draws but does not seed. | Chunk fs-estimates (`:437-440`), the 'Treatment Effect Estimates' table: only the HTML container id and scoped CSS; the visible HR, CI and count cells do not change. After the change it varies from render to render by the document alone; under pkgdown it changes once and then stays stable. | pkgdown only (same evidence) |
| `forestsearch-2024.Rmd:510` | knitr auto-print of `summaries$table`, a gt table from `format_bootstrap_table` (`R/bootstrap_summaries_helpers.R:157-159`, via `R/summarize_bootstrap_results.R:169`). It comes after `forestsearch_bootstrap_dofuture()` at `:478` (and `:301`). | No. The `seed = 8316951L` at `:481` belongs to the fixed-seed call, and after the change the bootstrap's restore on exit undoes its effect on the main stream. `summarize_bootstrap_results` (`:501`) neither seeds nor draws. | Chunk bootstrap-summary (`:499-511`), the bias-corrected estimates table: only the HTML id and scoped CSS; the visible bias-corrected HRs do not change. Old: the id comes from the 8316951L stream the bootstrap leaves. New: it comes from the document's own stream and varies by render; under pkgdown it changes once, then stays stable. | pkgdown only (same evidence) |
| `forestsearch-2024.Rmd:564` | knitr auto-print of `summaries$diagnostics_table_gt`, a gt table from `format_bootstrap_diagnostics_table` (`R/bootstrap_summaries_helpers.R:630-632`, via `R/summarize_bootstrap_results.R:283`). It comes after `:478`. | No. Between `:478` and `:564` run `plot_sg_weighted_km` (`:518`; weightedsurv `df_counting` with draws = 0 takes no `set.seed` branch), `figure_note` (`:541`) and `summarize_bootstrap_events` (`:554`). None of them seeds; the other two have no RNG at all. | Chunk bootstrap-diagnostics (`:562-565`), the 'Bootstrap Diagnostics Summary' table: only the HTML id and scoped CSS; the visible diagnostics do not change. After the change it varies from render to render by the document alone; under pkgdown it changes once and then stays stable. | pkgdown only (same evidence) |
| `forestsearch-2024.Rmd:632` | knitr auto-print of `metrics_tables`, a gt table from `cv_metrics_tables` (`R/cv_summary_tables.R:459-460`). It comes after `forestsearch_tenfold()` at `:606`, whose `%dofuture%` seed derivation moves the main stream forward 2 steps. Same `random_id()` chain. | No. The `seed = 8316951L` at `:611` is tenfold's own argument, used only inside the foreach body (`R/forestsearch_cross_validation.R:959`). `cv_metrics_tables` (`:622`) neither seeds nor draws. | Chunk kfold-cv (`:598-634`), the 'Cross-Validation Metrics' table: only the HTML id and scoped CSS. The visible CV values do not change, because they come from `set.seed(seed + 1000*ksim)` inside tenfold. After the change the id varies by render; under pkgdown it changes once, then stays stable. | pkgdown only (same evidence) |
| `survival-post-selection.qmd:287` | Autoprint of `summ_fs$table` (chunk boot-summary), a gt table from `summarize_bootstrap_results()` (`qmd:285`) -> `format_bootstrap_table` -> `gt::gt()` (`bootstrap_summaries_helpers.R:159`), with no id. gt `knit_print` -> `random_id()` -> `sample(letters, 10)`. It comes after `qmd:144` and `qmd:272`. | No. The table has no id, and the document has no `set.seed`. Today the stream is fixed only because the last bootstrap call (GRF, `seed = 8316951L` at `qmd:273`) leaves its `set.seed` in force; the change wraps `forestsearch_bootstrap_dofuture`. | The boot-summary table: only its HTML container id and the CSS selectors in its style block. The visible naive and bias-corrected HR cells do not change. Today the id is the same on every render; after the change it changes on every render. Nothing on the visible page changes, but the HTML source does, so site diffs show it. | pkgdown only (`.Rbuildignore:20`; `qmd:20`; `_pkgdown.yml:45`, `:582`) |
| `survival-post-selection.qmd:297` | Autoprint of `mr_estimates_table(fit_fs, boot_results = boots$FS, ...)` (chunk mr-table), which returns `format_bootstrap_table(...)` (`mr_estimates_table.R:153`) -> `gt::gt()`, with no id. Same print path; it continues the stream right after the `:287` draws. | No. The table has no seed of its own, and the document has no `set.seed`. Today the stream is fixed only by the `set.seed(8316951L)` that the last bootstrap call leaves in force (`bootstrap_calculations_helpers.R:27`, `bootstrap_analysis_dofuture.R:334`), which the change contains. | The mr-table: only the HTML container id and scoped CSS. The visible cells (naive HR, and MR two-term estimates with IJ intervals) do not change. Today the id is the same on every render; after the change it changes on every render. Nothing on the visible page changes. | pkgdown only (same evidence) |

**None of the 7 rows is in a vignette that CRAN builds.** The only CRAN-built vignette, `vignettes/forestsearch.qmd`, has 0 rows. Its one fixed-seed call (`forestsearch()` at `:228`) is followed by no unseeded draw, so its RNG-dependent outputs are unchanged and stay stable from render to render. All 7 rows are in two pkgdown-only articles, `forestsearch-2024.Rmd` (5) and `survival-post-selection.qmd` (2). Every row changes only a hidden gt HTML id and its scoped CSS; no visible number, table cell or figure changes.

**Build route caveats** (from the verifiers):
- pkgdown 2.2.1 runs `set.seed(1014)` before rendering an `.Rmd` article (`build_article(seed = 1014L)` -> `rmarkdown_render_with_seed`). That is why the `forestsearch-2024.Rmd` ids change once and then stay stable on the site. The verifier for `treatment_effect_definitions.qmd` found that this seed applies only to Rmd articles: `.qmd` articles render through quarto in a separate process with no seed. That fits the survival-post-selection verdict of "changes on every render".
- The CI job installs `cran::forestsearch` (`.github/workflows/pkgdown.yaml:36`) and builds with `install = FALSE` (`:40`). So the deployed page changes only once a released forestsearch that contains the change is installed, not at the `c40d5309` source. A local pkgdown build with the dev version installed shows the change immediately.

### Notable non-rows

- **`plot_km_band_forestsearch` (`forestsearch-2024.Rmd:812`)** draws (`KM_diff` calls `rnorm` with `draws_band = 20`), but it seeds itself first. `R/plot_km_band_forestsearch.R:349` passes no `seedstart`, so weightedsurv's default `set.seed(8316951)` runs in `KM_diff` right before `rnorm`. It is not on the reach list. Because that seed runs between `:606` and `:914`, the timing gt table at `:914` gets the same id in old and new code; its visible timings differ on every render anyway.
- **`forestsearch_tenfold` (`forestsearch-2024.Rmd:606`)** does draw from the caller's stream in the main process with no seed in between: the `%dofuture%` `seed = TRUE` derivation moves it forward 2 steps. Its CV results still do not depend on that stream, because `set.seed(seed + 1000*ksim)` at `R/forestsearch_cross_validation.R:959` comes before the body's first draw at `:960`. After the change, tenfold's own restore on exit undoes the advance.
- **`forestsearch_bootstrap_dofuture` (`forestsearch-2024.Rmd:478`; `survival-post-selection.qmd:272`; `actg175-glm-outcomes.qmd:522`)** draws nothing from the caller's stream before its own `set.seed(seed)` (`R/bootstrap_calculations_helpers.R:27`). `setup_parallel_SGcons` only sets the future plan, and parallelly's port choice saves and restores `.Random.seed`.
- **Draws that seed themselves from an argument given in the document** are not rows:
  - `fs_mr_inference` reseeds at `R/fs_mr_inference.R:615` and `:912`.
  - DINA folds seed at `R/dina.R:433`, and save and restore the caller's state (`:421-431`).
  - grf gets `seed = seedit` (`R/grf_helpers.R:75`).
  - `simulate_from_dgm(seed = ...)` seeds at `R/simulate_from_dgm.R:208`.
  - `generate_aft_dgm_flex` seeds at `R/generate_aft_dgm_main.R:248`.
  - `inject_mrct_structure` / `mrct_build_dgm(seed = 7L)` seed at `R/mrct_inject_dgm.R:97` and `:572`.
  - `compute_frontier_cis(seed = seed_rep5)` at `mrct-consistency.qmd:632` seeds at `R/frontier_cis.R:164`. That seed also fixes the later gt id draw at `mrct-consistency.qmd:649`, which is therefore seeded, not a row, and unchanged.
- **Fixed-seed calls whose leftover stream nothing reads.**
  - In `forestsearch.qmd` and `actg175-glm-outcomes.qmd`, every later draw is seeded by an explicit document argument. The DINA and GRF fits return before `R/forestsearch_main.R:3669`.
  - In `mrct-consistency-simulations.qmd`, the caller's state at `mrct_region_sims` (`:730`) comes from the `seed = 7L` builds (the last is `build_cell` at `:686`), not from 8316951, and each replicate seeds itself.
- **`paper_simulations.qmd`**: the gt prints after the simulation loops (`:1027`, `:1052`, `:1083`, `:1107`, `:1156`, `:1322`, `:1653`, `:2485`, `:2566`, `:2790`) draw from the session's own stream. That stream is fixed by `set.seed(1)` at `:729` and is never reached by the fixed seed, so they are not rows, and they stay stable and unchanged.
- **gt id draws with no fixed-seed call before them**: `causal_effects_brief_review.qmd` (`:145`, `:502`, `:884`) and `treatment_effect_definitions.qmd` (`:570`, `:822`, `:972`). They already vary from one render to the next, and the change cannot affect them.
- **No reaching call at all**:
  - `extreme_subgroups.qmd`: every draw is seeded. Its `set.seed(8316951)` at `:959` is the document's own seed, not a package call.
  - `methodology.qmd`: its only draws are `rnorm` at `:420-421`, after `set.seed(123)` at `:418`.
  - `guohe-nested-continuum.qmd`: every replicate reseeds at `:297` and in the package functions.
  - `biomarker_effects.qmd` and `treatment_effect_definitions.qmd`: no forestsearch function is called in evaluated code.
- **Aside from the guohe verifier, not about RNG leakage**: under `future.seed = TRUE` the Monte Carlo replicates run under L'Ecuyer-CMRG. They are therefore not the same random streams as a sequential Mersenne-Twister run with the same `set.seed` values. This bears on the article's claim at `:103-104` and `:281-282` that runs reproduce the manuscript's records.

## 5. Part 3 — The legacy GBSG simulation wrapper: usage inventory

The scope is everything defined in `R/sim_aft_gbsg.R` and `R/setup_gbsg_dgm.R`. Six collectors (L-A to L-F), each independently verified, covered the following:
- L-A: definitions and status;
- L-B: `R/`, `tests/`, `inst/`, `_pkgdown.yml` and README;
- L-C: `vignettes/`;
- L-D: `quarto/`;
- L-E: fs-post-selection;
- L-F: fs-glms-interpretable.

A seventh agent analysed what removal would involve, and a verifier checked that analysis. Only the verified results are used. Line numbers in `R/` are at the start SHA `3fbbc55c`, which is identical to `c40d5309` in `R/`, `vignettes/` and `quarto/`. This part reports usage only; nothing was changed. The full verified results are in `dev/reports/REPORT_cran_prep_seed_legacy_2026-11-11_detail.md`. Every calling quarto file is listed in `dev/reports/REPORT_cran_prep_seed_legacy_2026-11-11_quarto_callers.csv`.

### 5.1 Per function

There are 11 top-level functions, 2 nested closures and 3 top-level constants. Both files and all 9 wrapper Rd pages are blob-identical between `v0.4.0` (a lightweight tag on `83b6d708`) and `c40d5309`. So every export and internal status is the same as in the CRAN 0.4.0 release. All six NAMESPACE entries were already exported at `v0.1.0`. `NEWS.md:286` has an untagged 0.4.1 section, and the repository cannot tell whether 0.4.1 reached CRAN. The wrapper files are identical either way.

| Function | file:line | Exported or internal | In CRAN 0.4.0 (tag v0.4.0 = 83b6d708) | Deprecation status |
|---|---|---|---|---|
| `SEED_BASE`, `DAYS_PER_MONTH`, `DEFAULT_N_SUPER` (constants) | `R/sim_aft_gbsg.R:29`, `:32`, `:35` | internal; no NAMESPACE entry, no man page | present (identical) | not applicable |
| `cut_numeric` | `R/sim_aft_gbsg.R:51` | internal (`@keywords internal` `:50`; `man/cut_numeric.Rd`) | present, not exported | none |
| `cut_size` | `R/sim_aft_gbsg.R:69` | internal (`:68`; `man/cut_size.Rd`) | present, not exported | none |
| `create_gbsg_dgm` | `R/sim_aft_gbsg.R:199` | exported (`NAMESPACE:116`) | present and exported (`v0.4.0:NAMESPACE:102`) | deprecated at run time only: `.Deprecated()` at `:211-225` names `generate_aft_dgm_flex`, added in `834d8384` before v0.1.0. The roxygen and Rd do not say so; lifecycle is not used; there is no NEWS entry (the 0.2.2 'Superseded functions' entry was deleted in `fdecbe7b` before v0.4.0). |
| `.create_gbsg_dgm_` | `R/sim_aft_gbsg.R:236` | internal (untitled `@keywords internal` `:235`; no man page) | present, not exported | none |
| `simulate_from_gbsg_dgm` | `R/sim_aft_gbsg.R:789` | exported (`NAMESPACE:293`) | present and exported (`v0.4.0:NAMESPACE:266`) | deprecated at run time only: `.Deprecated()` at `:800-812` names `simulate_from_dgm`, in the same commit `834d8384`. The roxygen and Rd are silent, and there is no NEWS entry. |
| `.simulate_from_gbsg_dgm_` | `R/sim_aft_gbsg.R:823` | internal (untitled `@keywords internal` `:822`; no man page) | present, not exported | none (comment `:821`: "No deprecation warning") |
| `calibrate_k_inter` | `R/sim_aft_gbsg.R:1004` | exported (`NAMESPACE:92`) | present and exported (`v0.4.0:NAMESPACE:78`) | none |
| `eff_at` (nested in `calibrate_k_inter`) | `R/sim_aft_gbsg.R:1023` | local closure | present | not applicable |
| `get_dgm_with_output` | `R/sim_aft_gbsg.R:1076` | internal (`@keywords internal` `:1075`; `man/get_dgm_with_output.Rd`), yet listed in the public `_pkgdown.yml` section (`:339`) | present, not exported | none |
| `print.gbsg_dgm` | `R/sim_aft_gbsg.R:1157` | registered S3 method (`NAMESPACE:48`) | present and registered (`v0.4.0:NAMESPACE:43`) | none |
| `validate_k_inter_effect` | `R/sim_aft_gbsg.R:1229` | exported (`NAMESPACE:312`) | present and exported (`v0.4.0:NAMESPACE:285`) | none |
| `setup_gbsg_dgm` | `R/setup_gbsg_dgm.R:88` | exported (`NAMESPACE:286`) | present and exported (`v0.4.0:NAMESPACE:259`) | none. It is presented as the successor: "the new public API" (`:5`) and "the replacement for create_gbsg_dgm()" (`:26-27`). |
| `.rename` (nested in `setup_gbsg_dgm`) | `R/setup_gbsg_dgm.R:149` | local closure | present | not applicable |

`calibrate_k_inter`, `validate_k_inter_effect`, `setup_gbsg_dgm` and `print.gbsg_dgm` were never deprecated. The verified notes also list stale documentation that a cleanup could fix:
- `setup_gbsg_dgm` docs (`:6`, `:17`, `:31`) say it calls `create_gbsg_dgm`; it calls `.create_gbsg_dgm_` (`:106`), so no warning fires.
- The comment at `R/sim_aft_gbsg.R:233-234` names `scan_k_inter_effect`, which is defined nowhere.
- The `is.na(k_inter)` fallback in `get_dgm_with_output` (`:1101-1104`) cannot run.
- `setup_gbsg_dgm` hard-codes `5000L` and `8316951L` (`:94`, `:97`) instead of using the constants.

### 5.2 Callers

Only `setup_gbsg_dgm` and `calibrate_k_inter` are called by name from outside the two files, apart from one article's use of `validate_k_inter_effect` and its implicit `print.gbsg_dgm`. The other functions run only indirectly, through those calls (`eff_at` through `.calibrate_by_root` in `R/calibrate_helpers.R`).

**No caller anywhere searched** (package code, tests, vignettes, quarto, fs-post-selection, fs-glms-interpretable): `create_gbsg_dgm`, `simulate_from_gbsg_dgm`, `.simulate_from_gbsg_dgm_` and `get_dgm_with_output`. `create_gbsg_dgm` is still named in docs, README and manuscript prose (see below).

#### R/

| Function | Callers outside the two files | Callers inside the two files | Example |
|---|---|---|---|
| `setup_gbsg_dgm` | 2, both in `fs_fdr_report()`: `R/fs_fdr_report.R:256` (hom null) and `:262` (harm null) | 0 | `R/fs_fdr_report.R:256` `setup_gbsg_dgm(model = "alt", k_inter = 0, n_super = n_super)` |
| `calibrate_k_inter` | 1, in `fs_fdr_report()`: `R/fs_fdr_report.R:258`, inside `if ("harm" %in% null)`, which runs by default | 1: `get_dgm_with_output` (`:1093`) | `R/fs_fdr_report.R:258` `k_harm <- calibrate_k_inter(target_hr_harm = 1.0, model = "alt", use_ahr = FALSE)` |
| `.create_gbsg_dgm_` | 0 (one comment, `R/fs_fdr_report.R:485`) | 5, the hub: `create_gbsg_dgm` `:226`, `eff_at` `:1024`, `get_dgm_with_output` `:1108`, `validate_k_inter_effect` `:1254`, `setup_gbsg_dgm` `R/setup_gbsg_dgm.R:106` | `R/setup_gbsg_dgm.R:106` |
| `cut_numeric` / `cut_size` | 0 | `.create_gbsg_dgm_` only (`:281`, `:283`, `:284` / `:298`) | `R/sim_aft_gbsg.R:281` `dfa$z2 <- cut_numeric(dfa$age)` |
| `eff_at` | invoked indirectly by `.calibrate_by_root` at `R/calibrate_helpers.R:58` and `:74` | passed as `effect_fun` at `:1041` | `R/sim_aft_gbsg.R:1041` |
| `.rename` | 0 | 5 calls in `setup_gbsg_dgm` (`:155-159`) | `R/setup_gbsg_dgm.R:155` |
| `.simulate_from_gbsg_dgm_` | 0 | 1: `simulate_from_gbsg_dgm` (`:813`) | `R/sim_aft_gbsg.R:813` |
| constants | 0 | 5 sites in `R/sim_aft_gbsg.R` (`:205`, `:208`, `:242`, `:245`, `:273`) | `R/sim_aft_gbsg.R:273` |
| `print.gbsg_dgm` | no call by name; reached by S3 dispatch on objects from `.create_gbsg_dgm_` (`:735`) and `setup_gbsg_dgm` (`:182`), since no `aft_dgm_flex` print method exists | 0 | none |

- **`fs_fdr_report()`** is the only package function outside the two files that calls the wrapper, and it was the same at v0.4.0. No `R/` function calls `fs_fdr_report()`. Its only test caller is `test-plan-restored.R:168`, and its only quarto caller is `quarto/calibration/fs_fdr_report_smoke.qmd` (`:109`, `:168`).
- **`R/betaHhat_truth.R`** references no wrapper name. It reads `df_super` and `model_params$censoring`, which are the general `aft_dgm_flex` contract that `generate_aft_dgm_flex` also produces, not a setup-specific dependency.
- The two deprecated entry points are never called by wrapper code, so no internal path triggers their warnings.
- Non-call couplings:
  - `R/globals.R:196-208`, a GBSG `globalVariables` block whose header names a nonexistent file. 16 of its names are used only by the wrapper.
  - Roxygen mentions: `create_gbsg_dgm` in `R/oc_analyses.R:109`, `:123` and `:413`, two of them `\dontrun`; `simulate_from_gbsg_dgm` in `R/run_simulation_analysis.R:6` and `:198`; `setup_gbsg_dgm` in 7 roxygen lines; `calibrate_k_inter` in 4 roxygen lines.
  - The `stop()` string at `R/fs_dgm_feasibility.R:130`.

#### tests/

- `calibrate_k_inter`: `tests/testthat/test-betaHhat-contract.R:419`, in the cached helper `.small_surv_dgm` (`:415`). T12 (`:449-459`) has no skip, so it runs under R CMD check and on CRAN. T11, T13 and the tail of T14 sit behind `skip_if_not(.have_shims())`, which is FALSE under check because `quarto/` is build-ignored.
- `setup_gbsg_dgm` is run by 3 files:
  - `test-knoise-population.R:8`: helper `.make_gbsg_noise_dgm`, used in 5 `test_that` blocks with no skips.
  - `test-betaHhat-contract.R:421`: in `.small_surv_dgm`, gated as for `calibrate_k_inter` above (T12 unconditional).
  - `test-plan-restored.R:168`: transitively, via `fs_fdr_report(null = "hom")`, with `skip_on_cran()` at `:165`.
- Indirectly, on CRAN, T12 also runs `.create_gbsg_dgm_`, `eff_at`, `cut_numeric`, `cut_size` and `.rename`. `setup_gbsg_dgm` also runs on CRAN through all 5 survival tests in `test-knoise-population.R`.
- No test calls `create_gbsg_dgm`, `simulate_from_gbsg_dgm`, `get_dgm_with_output`, `validate_k_inter_effect` or `print.gbsg_dgm`.
- `test-fs-dgm-feasibility.R:26` only builds a `structure(list(), class = "gbsg_dgm")` stub for `expect_error`. The tracked `.rds` fixtures hold no `gbsg_dgm` object.

#### vignettes/ and vignettes/articles/

| File | Build | Wrapper use |
|---|---|---|
| `vignettes/forestsearch.qmd` | CRAN | none |
| `vignettes/articles/paper_simulations.qmd` | pkgdown only | All 7 evaluated wrapper calls in the vignettes: `setup_gbsg_dgm` at `:507`, `:614`, `:645`, `:699`; `calibrate_k_inter` at `:604`, `:635`; `validate_k_inter_effect` at `:683`. Also an implicit `print.gbsg_dgm` from `print(dgm_alt)` at `:523`. Not run: `setup_gbsg_dgm` `:375` and `:1756` (eval false) and `:274` (display block); `calibrate_k_inter` `:1749` (eval false), `:1753` (commented out) and `:273` (display block). The mid-document YAML at `:358-360` makes `:359` the deployed page title. |
| `vignettes/articles/treatment_effect_definitions.qmd` | pkgdown only | `setup_gbsg_dgm` on 6 lines, none evaluated: an `eval=FALSE` call at `:1034`, a comment at `:460`, gt string literals at `:552` and `:785`, and prose at `:163` and `:452`. The tracked `.html` copy has 8 rendered lines. |
| `vignettes/articles/extreme_subgroups.qmd` | pkgdown only | `calibrate_k_inter` in prose only (`:202`, `:485`) |
| `vignettes/articles/biomarker_effects.qmd` | pkgdown only | `create_gbsg_dgm` in prose only (`:66`; the tracked `.html` has it at `:2413`) |
| `actg175-glm-outcomes.qmd`, `causal_effects_brief_review.qmd`, `forestsearch-2024.Rmd`, `guohe-nested-continuum.qmd`, `methodology.qmd`, `mrct-consistency.qmd`, `mrct-consistency-simulations.qmd`, `survival-post-selection.qmd` (all under `vignettes/articles/`) | pkgdown only | none (nor in their tracked `.html` renders) |

- 9 of the 13 functions are never named in `vignettes/`.
- 4 are never executed, even indirectly: `create_gbsg_dgm`, `simulate_from_gbsg_dgm`, `.simulate_from_gbsg_dgm_` and `get_dgm_with_output`.
- `.create_gbsg_dgm_`, `cut_numeric`, `cut_size`, `eff_at`, `.rename` and the constants run indirectly through `paper_simulations.qmd`'s 7 calls.
- Documentation finding: `paper_simulations.qmd` calls its `eval: false` `generate_aft_dgm_flex` block (`:395-440`) "equivalent" (`:392`), and it describes a "binary discretisation" (`:366`). The wrapper actually uses quartile codes 1-4 (`cut_numeric`) and a 3-level `cut_size` (`R/sim_aft_gbsg.R:281-284`, `:298`), while the block uses median splits, so the two are not equivalent.

#### quarto/ (forestsearch)

**Classification rule.** A file is archive/legacy/broken if any directory component matches, case-insensitively, `(^|[^a-z])(archive[sd]?|legacy|broken|old|superseded|deprecated|obsolete)([^a-z]|$)`. Everything else is current. The rule covers 1111 in-scope tracked files (`.qmd`, `.R`, `.r`, `.Rmd`, `.rmarkdown`, `.md`): 922 current and 189 archive. A looser substring rule flagged no further directory.

The archive directories, with their in-scope file counts:
- `quarto/archive` (13)
- `applications/actg175/_archive` (20)
- `applications/gbsg/_archive` (16)
- `applications/gbsg/_broken` (16)
- `legacy calibration files` (9)
- `simulations/_legacy` (16)
- `simulations/actg175/_archive_pre_fsparams` (3)
- `simulations/actg175/binary_legacy` (89)
- `simulations/gbsg_redux/legacy` (7)

**Only `setup_gbsg_dgm` and `calibrate_k_inter` are called in `quarto/`.**

| Class | Files calling | `setup_gbsg_dgm` files / call lines | `calibrate_k_inter` files / call lines | Example |
|---|---|---|---|---|
| current | 281 (272 `.qmd`, 8 `.R`, 1 `.rmarkdown`) | 281 / 297 | 277 / 282 | `quarto/calibration/calibration_survival_fourway.qmd:260` (calibrate), `:269` (setup) |
| archive/legacy/broken | 13 (all `.qmd`) | 13 / 17 | 10 / 10 | `quarto/applications/gbsg/_broken/gbsg_analysis.qmd:1082` (setup), `:1368` (calibrate) |

Current calling files by directory, as files/lines for calibrate | setup:
- `calibration` 3/3 \| 6/9
- `qc` none \| 1/1
- `simulations` (top level) 1/1 \| 1/2
- `simulations/gbsg` 105/105 \| 105/106
- `gbsg/sims1_claude_summary` 3/3 \| 3/3
- `simulations/gbsg_020` 74/74 \| 74/75
- `gbsg_020/scripts_dinamr` 6/6 \| 6/11
- `simulations/gbsg_redux` 85/90 \| 85/90

The 4 current files that call `setup_gbsg_dgm` only are `calibration/{calibration_survival_leff_grid,fpr_calibration_demo,fs_fdr_report_smoke}.qmd` and `qc/qc_null_classification.R`.

The archive calling files are in:
- `_broken` (3)
- `archive/reproduce_table1_m1_working.qmd` (1)
- `legacy calibration files` (3)
- `simulations/gbsg_redux/legacy` (6)

Five archive directories have no hits of any kind.

**Mention-only or display-only files:** 23 current and 7 archive.
- `simulate_from_gbsg_dgm` appears only in one archive file, `legacy calibration files/run_simulation_analysis.R:6` and `:76`, in comments.
- `.create_gbsg_dgm_` has 0 calls; it is mentioned on 19 lines in 9 current files.
- `print.gbsg_dgm` has no explicit call and no detectable dispatch: the 294 calling files contain no print, summary, str or bare auto-print of a DGM object.
- `SEED_BASE` in `gbsg_020/scripts_dinamr/grf_mechanism.R:34` is a local redefinition, not the package constant.
- Some archive files duplicate current ones: the three `.qmd` files in `legacy calibration files` duplicate `quarto/calibration/`, and `archive/reproduce_table1_m1_working.qmd` duplicates `simulations/reproduce_table1_m1.qmd`.

The CSV holds all 294 calling files: 281 current, each with its removal group from 3.4, and 13 archive.

#### inst/, _pkgdown.yml, README

- `inst/` has 0 hits across its 30 tracked files, even for a case-insensitive `gbsg`.
- `_pkgdown.yml` has 7 entries under 'Data Generating Mechanisms' (`:326`):
  - `:332` `setup_gbsg_dgm`
  - `:333` `create_gbsg_dgm`
  - `:334` `print.gbsg_dgm`
  - `:338` `simulate_from_gbsg_dgm`
  - `:339` `get_dgm_with_output` (an internal topic)
  - `:342` `calibrate_k_inter`
  - `:344` `validate_k_inter_effect`
- `README.md:131` lists the deprecated `create_gbsg_dgm()` in the key-functions table. `README.Rmd` is not tracked. No other function is in the README.

#### fs-post-selection: manuscript_rev1 and its generators

- **Rendered paper** (main, supplement, derivations): no wrapper call. These documents read `_payloads/*.rds`. The supplement also runs `library(forestsearch)` and `fs_identification_*`, which reach no wrapper function.
- **Text dependency**: the supplement's code-availability sentence names `create_gbsg_dgm()`. It is at `postselection_subgroup_inference_supplement.qmd:1184`, in the kept `.tex` at `:3296` and in the PDF. Across manuscript_rev1, `create_gbsg_dgm` appears in 17 text files (17 lines: 11 plain, 6 LaTeX-escaped `.tex`; all but the current supplement's `.qmd` and `.tex` and one survey report are archived drafts) and 11 binary files, all prose, with 0 calls.
- **Executable calls in manuscript_rev1**: 12 lines (6 `calibrate_k_inter`, 6 `setup_gbsg_dgm`), only in `dgm_truth_targets{,_h10,_h15}.qmd` at `:150`, `:156`, `:227` and `:230`.
  - These are standalone auxiliary documents. Neither the main paper nor the supplement includes them.
  - `_quarto.yml` has no render list, so a whole-project `quarto render` would execute them. The documented workflow renders the paper file by file.
- **Generators** of the current render's payloads: 37 forestsearch `quarto/` documents. 34 are in `simulations/gbsg_020/`; the others are the two application analyses and one ACTG175 sweep.
  - 34 of them contain `calibrate_k_inter` and `setup_gbsg_dgm`.
  - 29 execute them: 20 batch documents behind `run_mode == "batch"` and 9 sweeps behind `if (run_sweep)`, which is set TRUE.
  - The 5 combiners set `run_mode` to "combine", so they do not run the calls. The 3 others (the two application analyses and the ACTG175 sweep) have no wrapper call.
  - Example: `sim_fs_maxeffCons_fb_mr_m1_h10_knoise0_n500_batch_1_50.qmd:318` (calibrate) and `:321` (setup).
  - Through these calls the generators also run `.create_gbsg_dgm_`, `cut_numeric`, `cut_size`, `eff_at` and `.rename`. None of this work calls `create_gbsg_dgm`, `simulate_from_gbsg_dgm`, `.simulate_from_gbsg_dgm_`, `get_dgm_with_output`, `validate_k_inter_effect` or `print.gbsg_dgm`.
- `payload_manifest.qmd` and `manuscript_render_manifest.html` are stale. The 65 gbsg_redux documents they name all call the wrapper.
- **Outside manuscript_rev1** (code call lines):
  - `analysis/`: 36 lines in 18 files. The removal verifier located 18 calling `.qmd` files under `analysis/simulations/gbsg/resampling_comparisons/`, last touched 2026-06-21, not referenced from manuscript_rev1; their currency is unclassified.
  - `manuscript/`, `manuscript_legacy/` and `manuscript_v40/`: 12 lines in 3 files each, all copies of `dgm_truth_targets*`.

#### fs-glms-interpretable: quarto/, manuscript/, dev/payloads/

- Across the 576 tracked files, every wrapper function has 0 calls and 0 name mentions. The only near-mention is the concept phrase "via `k_inter` calibration" at `quarto/gbsg/SPEC_gbsg_grf_dina_2026-09-17.md:102`.
- **Indirect dependence** runs through three payload streams produced in forestsearch by code that calls `calibrate_k_inter` and `setup_gbsg_dgm`:
  - `identifiers_metrics.csv`: read at `manuscript/_05_operating_characteristics.qmd:14` and `_S4_operating_characteristics_supplement.qmd:16`. Its bundles come from the gbsg_020 field template at `ca9b0a62`/`69ca797e`/`dbc9444c`/`1507bfc0`.
  - `idsweep_metrics.csv`: read at `_02_identification.qmd:4` and `_S4:17`. It comes from the same template at `b15c8b7a`.
  - `declcalc0_summary_tables.csv`: read at `_09_field_constructions_supplement.qmd:824`. It comes from `scripts_dinamr/declcalc0_run.R` (calibrate `:148`, setup `:156`, `:163`).
- The GLM OC payloads use `generate_glm_dgm`/`simulate_from_glm_dgm`, and the application payloads come from real data.
- Outside these three directories, 20 tracked files name a wrapper function. Real calls among them are only in `dev/reference/.../dgm_truth_targets{,_h10,_h15}.qmd`.

### 5.3 Summary

| Function | Exported | In CRAN 0.4.0 | Called by package code | Called by a vignette or article | Called by current manuscript work | Called only by archive or legacy documents |
|---|---|---|---|---|---|---|
| `setup_gbsg_dgm` | yes | yes | yes: `fs_fdr_report()` (`:256`, `:262`) | yes: `paper_simulations.qmd`, 4 evaluated calls (pkgdown only) | yes: 29 executing generators of manuscript_rev1 payloads; 6 call lines in the standalone `dgm_truth_targets*.qmd`; fs-glms-interpretable indirectly, through payloads | no: 281 current quarto files call it |
| `calibrate_k_inter` | yes | yes | yes: `fs_fdr_report()` (`:258`); also `get_dgm_with_output` inside the wrapper | yes: `paper_simulations.qmd`, 2 evaluated calls (pkgdown only) | yes: the same 29 generators; 6 call lines in `dgm_truth_targets*.qmd`; fs-glms-interpretable indirectly | no: 277 current quarto files call it |
| `validate_k_inter_effect` | yes | yes | no | yes: `paper_simulations.qmd:683` (evaluated, pkgdown only) | no | no: its only caller is a current article |
| `print.gbsg_dgm` | S3 method registered (`NAMESPACE:48`) | yes | no: dispatch only | yes, by dispatch: `print(dgm_alt)` at `paper_simulations.qmd:523` | no | no: no other caller |
| `create_gbsg_dgm` | yes (deprecated) | yes | no: `\dontrun` examples only | no: prose only (`biomarker_effects.qmd:66`) | no call, but the manuscript_rev1 supplement names it (`:1184`) | no: no callers at all |
| `simulate_from_gbsg_dgm` | yes (deprecated) | yes | no | no | no | no: no callers at all (one archive quarto file mentions it in comments) |
| `get_dgm_with_output` | no (internal; listed in `_pkgdown.yml:339`) | yes, not exported | no | no | no | no: no callers at all |
| `.create_gbsg_dgm_` | no | yes, not exported | yes: inside the wrapper only (5 callers) | indirectly, through `paper_simulations.qmd`'s 7 calls | indirectly, through the 29 generators | no |
| `.simulate_from_gbsg_dgm_` | no | yes, not exported | only by `simulate_from_gbsg_dgm` (`:813`) | no | no | no: no callers beyond that |
| `cut_numeric`, `cut_size`; nested `eff_at`, `.rename`; constants `SEED_BASE`, `DAYS_PER_MONTH`, `DEFAULT_N_SUPER` | no | yes, not exported | only inside the wrapper (`eff_at` also through `.calibrate_by_root`) | indirectly, through `paper_simulations.qmd` | indirectly, through the 29 generators | no |

### 5.4 What removal would involve

This section is from the verified removal analysis. It describes what removal would require; nothing was removed or edited. The task asked for reading only. The verifier went further: to settle answers that reading had left open, it ran short `Rscript --vanilla` scripts that sourced the `c40d5309` `R/` files (`t2ro/LH-adv/equiv.R`, `equiv2.R`, `equiv3.R`, `equiv4.R`, each under 10 s). The answers that rest on those runs are marked **(run)** below.

**Files, exports, help pages.**
- `R/sim_aft_gbsg.R` (1316 lines) and `R/setup_gbsg_dgm.R` (190 lines) would go. `DESCRIPTION` has no Collate field, so it needs no edit.
- Six NAMESPACE lines would go: `:48` `S3method(print,gbsg_dgm)`, plus exports at `:92`, `:116`, `:286`, `:293` and `:312`. All six are in CRAN 0.4.0, and four of them were never deprecated. No import would be lost.
- Nine Rd pages would go: `create_gbsg_dgm`, `simulate_from_gbsg_dgm`, `calibrate_k_inter`, `validate_k_inter_effect`, `setup_gbsg_dgm`, `print.gbsg_dgm`, `get_dgm_with_output`, `cut_numeric` and `cut_size`.
- These regenerated pages would need link fixes to avoid an Rd missing-link WARNING:
  - `man/calibrate_k_treat.Rd:59` and `:115` (`R/calibrate_k_treat.R:9`, `:90`)
  - `man/fs_fdr_report.Rd:182` (`R/fs_fdr_report.R:172`)
  - `man/run_simulation_analysis.Rd:42` and `:260` (`R/run_simulation_analysis.R:205`, `:339`)
- These `\dontrun` examples and text would be stale but not checked:
  - `compute_dgm_cde.Rd:10`, `:46`
  - `summarize_simulation_results.Rd:32`
  - `simulate_from_dgm.Rd:193`
  - `check_censoring_dgm.Rd:72`
  - `run_simulation_analysis.Rd:204`, `:242`
  - `fs_dgm_feasibility.Rd:76`

**Tests.**
- `test-knoise-population.R` survival tests 1-5 (`:32-96`) have no skip and would ERROR. They can be ported by drawing the noise caller-side onto a flex `dgm$df_super`, as `setup_gbsg_dgm.R:125-135` does; the values are identical **(run)**. The GLM tests are unaffected.
- The `test-betaHhat-contract.R` fixture `.small_surv_dgm` (`:414-426`) feeds T12, which would ERROR under check and on CRAN. It can be rebuilt with `find_k_inter_for_target_hr(1.0, ...)` plus `generate_aft_dgm_flex` on a caller-built GBSG frame **(run)**.
- `test-plan-restored.R:164-173` (`skip_on_cran`) follows whatever `fs_fdr_report` becomes.
- `test-fs-dgm-feasibility.R:22-28` needs no change if the `stop()` keeps the substring "Survival DGMs".

**Dependent package functions.**

| Function | What it would have to do |
|---|---|
| `fs_fdr_report()` (`R/fs_fdr_report.R:181`) | Go together with the wrapper, or rebuild both nulls with flex. A rebuild is exact for the outcome model and truths **(run)**. The legacy censoring cannot be expressed through flex arguments, so simulated censoring times change and the documented A ~ 0.138 / 0.278 at n = 700 (`:163-164`) would need re-establishing. Relink `@seealso`. |
| `.fs_fdr_hr_row()` (`R/fs_fdr_report.R:482-496`) | Read `hazard_ratios$harm_subgroup`, `$no_harm_subgroup`, `$AHR_harm` and `$AHR_no_harm`, which legacy objects already carry, so the switch works today. `k_inter` must come from the caller. Without the change, every cell is NA on a flex DGM. |
| `run_simulation_analysis()` | Change the default `confounders_base = c('v1',...,'v7')` (`:398`): only `.create_gbsg_dgm_` builds those columns. Retarget the docs and example. |
| `simulate_from_dgm()`, `check_censoring_dgm()`, `compute_dgm_cde()`, `summarize_simulation_results()`, `calibrate_k_treat()` | Docs only: `\dontrun` examples, text and the two `calibrate_k_treat` links. |
| `get_dgm_hr()`, `format_oc_results()`, `invert_dgm_hrs()`, `build_classification_table()`, `build_estimation_table()`, `interpret_estimation_table()` | No code change forced. `get_dgm_hr`'s top-level fallback stays reachable for a flex null DGM, where `get_dgm_hr(., 'hr_Hc')` returns NA; a legacy null returns `hr_causal`. |
| `.calibrate_by_root()`; comments in `calibrate_glm_interaction()`, `find_k_inter_for_target_hr()` | These stay; only the comments need editing. |
| `fs_dgm_feasibility()` | Reword the `@note` and the `stop()` (`:128-131`) while keeping "Survival DGMs". |
| `fs_build_eval_frame()` (`R/betaHhat_truth.R`) | No change. |
| The replacements (`generate_aft_dgm_flex`, `find_k_inter_for_target_hr`, `sensitivity_analysis_k_inter`, `calibrate_k_treat`) | Nothing forced. Gaps for current documents: (1) a GBSG covariate-prep helper (the job of `cut_numeric`/`cut_size`); (2) the legacy censoring model, the only gap that changes simulated data; (3) a print method for `aft_dgm_flex`; optionally (4) a noise argument and (5) an AHR target for `k_inter`. Not needed by any current document: seed or `...` forwarding in `find_k_inter_for_target_hr` (every calibration runs at its defaults, 8316951 / 5000), `k_z3`, `use_rand_params`. |

**Other touchpoints.**
- `_pkgdown.yml`: 7 lines (`:332`-`:344`) would go, since pkgdown errors on missing topics.
- `NEWS.md` would need a removal entry.
- `R/globals.R:196-208`: removing the 16 orphan names is optional.
- `README.md:131`.
- `CLAUDE.md:53`.
- `inst/` and `vignettes/forestsearch.qmd`: 0 references.
- `paper_simulations.qmd`: its 7 evaluated calls and the implicit print would make the pkgdown article build fail. It has about 49 lines after `:440` that read legacy fields.
- `treatment_effect_definitions.qmd`: one `eval=FALSE` call and some text.
- R CMD check: NAMESPACE must be regenerated, the Rd links fixed, and the T12 and knoise tests ported or removed. Every example that names a removed function is `\dontrun`, so check reports nothing for them.
- API: removal breaks users of six names exported in CRAN 0.4.0. Reverse dependencies were not checked.
- The replacement call suggested in the `.Deprecated` message at `R/sim_aft_gbsg.R:215-223` cannot run as written: `id =` and `confounders =` are unused arguments, and `continuous_vars` and `factor_vars` are missing.

**Current documents.** Groups G1-G7 cover the 281 current quarto callers: 248 + 14 + 5 + 5 + 2 + 5 + 2 = 281. **(run)** marks an answer the verifier settled by running code. The "exact" statements for the other groups are the verifier's: it applied the equivalence that the run established to each group's build.

| Document(s) | Wrapper step | `generate_aft_dgm_flex()` + `simulate_from_dgm()` offers the same step? |
|---|---|---|
| G1: 248 forestsearch `quarto/` files (`simulations/gbsg/` 104, `gbsg/sims1_claude_summary/` 3, `gbsg_020/` 56, `gbsg_redux/` 85) | `calibrate_k_inter(use_ahr = FALSE)` at seed 8316951 and n_super 5000, then `setup_gbsg_dgm`, `compute_dgm_cde` and `simulate_from_dgm`. 241 files read top-level `hr_causal`/`hr_H_true`/`hr_Hc_true`. | **No as a drop-in, though only censoring and field names differ.** Once the caller builds `y` and z1..z5, the outcome model, subgroup, truths and calibration are exact: HR and AHR differences are 0, and k_inter agrees to 2e-7 **(run)**. Legacy censoring cannot be expressed in flex, so only `C_sim` changes. Truths sit in `hazard_ratios`, so top-level reads become NA. |
| G2: 14 population-noise sweeps (`gbsg_020/maxeffCons_mr_coverage_sweep_*`, the `*_standalone2` files and one TEST file) | As G1, plus `setup_gbsg_dgm(..., k_random_noise, noise_seed)`; the documents assert the `noise_*` fields. | **No as a flex argument.** The caller-side draw onto `dgm$df_super` reproduces the noise exactly **(run)**, and the caller would set the `noise_*` fields. Otherwise as G1, with censoring the only simulation difference. |
| G3: 5 `z1_quantile` variants (`gate_d2_cim_{field,ij,unset}.qmd`, `declcal_preflight.R`, `grf_mechanism.R`) | `calibrate_k_inter(..., z1_quantile)` and `setup_gbsg_dgm(..., z1_quantile)` | **No as a drop-in** (G1 reasons). The `z1_quantile` step itself is a yes: a caller-built z1 plus `subgroup_cuts`, and `find_k_inter_for_target_hr` takes them. |
| G4: 5 files (`declcal_run.R`, `declcal_c0approx_run.R`, `declcalc0_run.R`, `nullid_designpoint.R`, `sim_fs_maxeffCons_fb_mr_field_m1_template.qmd`), which produce the fs-glms-interpretable payloads | `uniroot` for `k_treat` over `setup_gbsg_dgm(model = 'null')$hr_causal`; alt designs via `calibrate_k_inter`. They read `subgroup_info$fs_harm_true`/`grf_harm_true`, `b_hr` and `AHR`. | **No as a drop-in.** `k_treat` calibration through `calibrate_k_treat` is identical (0.6853424) **(run)**, and the null model is a yes. Flex has no `fs_harm_true`/`grf_harm_true`, so the 'labels absent' gate passes vacuously. Censoring and field names differ as in G1. |
| G5: 2 identical Table 1 reproductions (`reproduce_table1_m1.qmd`, top level and `gbsg/`) | `calibrate_k_inter(..., k_treat = 1, n_super = 5000L, seed = 8316951L)`, then the alt and null `setup_gbsg_dgm` builds and `run_simulation_analysis` | **No overall.** Calibration and truths match exactly (the documents' seed and n_super equal `find_k_inter_for_target_hr`'s defaults). Simulated operating characteristics do not, because of censoring. |
| G6: 5 calibration and QC files (`calibration_survival_fourway`, `selection_survival_frontier`, `calibration_survival_leff_grid`, `fpr_calibration_demo`, `qc/qc_null_classification.R`) | `calibrate_k_inter(... cens_type = 'weibull')`, then the alt and/or null `setup_gbsg_dgm`, then simulation with z1..z5, size and grade3 | **No**, the same gap as G5. Covariates are caller-built, truths and calibration are exact, and censoring cannot be reproduced. The null build is a yes. |
| G7: `fs_fdr_report_smoke.qmd`, `gbsg_aft_dgm_appendix.qmd` | Smoke: `setup_gbsg_dgm`, then `fs_fdr_report()`. Appendix: setup, calibrate and setup; it reads `effect_modifiers`, `hr_*_true`, `AHR_*_true` and `cens_params`. | **No.** The smoke file follows `fs_fdr_report`. In the appendix the builds and truths are exact, but flex has no `effect_modifiers`, stored `k_inter` or `cens_params`, and its censoring values differ. |
| `vignettes/articles/paper_simulations.qmd` (pkgdown only) | 7 evaluated calls (`:507`, `:604`, `:614`, `:635`, `:645`, `:683`, `:699`), the print at `:523`, `simulate_from_dgm` at `:725`, and about 49 legacy-field reads | **No overall.** Step by step: build yes; print no (no `aft_dgm_flex` method); `:604` yes; `:635` (AHR target) has no packaged function, but a closed form matches 1.301622 **(run)**; `:683` partly, via `sensitivity_analysis_k_inter`; null yes; field reads need renaming; censoring differs. |
| `vignettes/articles/treatment_effect_definitions.qmd`, chunk `:1030` (`eval=FALSE`) | `setup_gbsg_dgm` at `:1034`; AHR and CDE from `df_super_rand`; `simulate_from_dgm` at `:1075` | **Yes** for the effect definitions, once the caller builds the covariates; the replicate differs only in censoring. Existing bug: `df_super_rand` has `flag.harm`, not `flag_harm`, so `:1042-1045` and `:1052-1055` subset by NULL. |
| fs-post-selection `manuscript_rev1/dgm_truth_targets{,_h10,_h15}.qmd` (standalone, not rendered with the paper) | `calibrate_k_inter` at `:150` and `:227`; `setup_gbsg_dgm(n_super = 25000L, seed = 8316951L)` at `:156` and `:230`; truth reads; no simulation | **Yes, with numerically identical output** (previously "unknown"). The differences were 0 at n_super 25000, and k agreed within 2e-7 **(run)**. Fields would need renaming. |
| Generators of the manuscript_rev1 payloads (a subset of G1 and G2) | `calibrate_k_inter` then `setup_gbsg_dgm` (e.g. `:318`/`:321`), executed in batch mode or when `run_sweep` is set | **No as a drop-in.** Truths, calibrated k_inter, event times, treatment draws and noise would be unchanged, but censoring times would change, and with them the operating characteristics: the event share moved about 0.3 percentage points at n = 700 **(run)**. |
| manuscript_rev1 rendered documents (main, supplement, derivations) | No wrapper call. The code-availability sentence names `create_gbsg_dgm()` (`supplement.qmd:1184`, `.tex:3296`, PDF). | Not applicable. The sentence would need rewording. |
| fs-post-selection `analysis/simulations/gbsg/resampling_comparisons/` (18 files; last touched 2026-06-21; currency unclassified) | `calibrate_k_inter` and `setup_gbsg_dgm` in the G1 pattern | **No as a drop-in**, for the G1 reasons; truths exact as in G1. These files would fail at run time. |
| fs-glms-interpretable `quarto/`, `manuscript/`, `dev/payloads/` | No call. The dependence is indirect, through payloads from G4 code. | Not applicable. Regenerating those payloads would meet the G4 answer: no as a drop-in, with truths exact and censoring different. |

The 13 archive quarto callers are outside groups G1-G7; the removal analysis gave no answer for them.

## 6. Part 4 — gates

Checked on the staged tree just before the report commit, and again after it (in the final printout of this run):

- **`git status --porcelain` is empty** after the commit.
- **The two sibling repositories are unchanged.** fs-post-selection is `5d165acc` on `main` and fs-glms-interpretable is `dc94a17e` on `main`, each with an empty `git status --porcelain`, as at the run's start (20:08). No tracked or untracked file in either is newer than the run's start marker.
- **`git diff 3fbbc55c --stat` lists only files this report names, and nothing under `vignettes/` or `quarto/`.** The Appendix lists every path:
  - the TASK file;
  - Part 1's 16 files (`DESCRIPTION`, `NEWS.md`, 13 `R/` files, the new test file);
  - this report and its two detail files.
- **Copied** to `~/Downloads/REPORT_cran_prep_seed_legacy_2026-11-11.md`.
- Nothing pushed or merged.

## OPEN ITEMS

1. **The NEWS example differs from the task's wording (measured).** The task asked for "a bootstrap run straight after `forestsearch()`" as the example of a call that now depends on the caller's seed. `forestsearch_bootstrap_dofuture()` seeds itself (`seed = 8316951L`), so its output is identical between the trees and repeats exactly (Step 5). The entry therefore uses `sample()` and `simulate_from_dgm(seed = NULL)` as examples, says the bootstrap is not affected, and names the 15 functions.
2. **Which exported functions restore.** Exactly the 15 that can call `subgroup.consistency()` (§3.2).
   - Other exported functions that set a seed from an argument still leave the caller's state as before. Examples: `generate_aft_dgm_flex(seed = 8316951)` by default, `simulate_from_dgm(seed = )`, `fs_oc_predict(seed = )`, `fs_mr_inference(seed = )` called directly, and `generate_bootstrap_synthetic(seed = 123)`.
   - That is within the task's scope ("seeds that come from a user-supplied argument are listed, not changed"). If Larry wants every exported function to leave the state alone, the same one-line guard applies to the functions in rows 1-44 of §3.1 whose seed is not restored; that is a wider behaviour change.
3. **Multisession and user code.** Inside a user's own `%dofuture%` loop, a guarded call (e.g. `run_simulation_analysis()`, which is documented for that use) is the outermost call in the worker. It restores the worker's state when it returns, so a draw the user's loop makes after it in the same iteration now comes from the future's own stream. This is the documented consequence (set a seed), not a change in a package value.
   - Inside the package, the only future body with several guarded calls per evaluation is `forestsearch_tenfold()`, and `.fs_rng_hold()` keeps its values unchanged (§3.2, measured).
   - A future body added later that makes more than one guarded call and draws between them would need the same one line.
4. **Installed package.** The copy in the user library (`~/R/x86_64-pc-linux-gnu-library/4.6/forestsearch`) is 0.4.1.9002. Multisession workers in Larry's own sessions load that copy until the branch is installed, so none of the changes on `cran-prep` (task 1 or this task) take effect in workers before then. No install was made into the user library.
5. **Part 2 (pkgdown only; no edit).** The 7 places are `gt` tables printed without an `id`, so `gt` draws a random 10-letter HTML id.
   - After the change, the id comes from the document's own stream instead of the 8316951 stream. The visible page is the same either way.
     - `forestsearch-2024.Rmd`: pkgdown sets `seed = 1014` before rendering an `.Rmd` article, so there the id changes once and is then stable; under a plain Knit it changes on every render.
     - `survival-post-selection.qmd`: one verifier found that pkgdown's seed does not reach Quarto articles, so the id changes on every render.
   - Options for Larry, none taken: a `set.seed()` at the top of the two articles; or a fixed `id` in the package's `gt::gt()` calls (for example `R/bootstrap_summaries_helpers.R:159`, `R/summary_utility_functions.R:1045, :1196`, `R/cv_summary_tables.R:459`).
6. **Part 3 (Larry's decision).**
   - Removing the wrapper would remove 6 functions exported in CRAN 0.4.0. Only 2 of them were ever deprecated.
   - It would break `fs_fdr_report()` and two test files, which would ERROR under check. It would also break 281 current quarto documents, 29 of the generators of fs-post-selection's `manuscript_rev1` payloads, and fs-post-selection `analysis/` documents.
   - The verifier showed by a short run that `generate_aft_dgm_flex()` reproduces the wrapper's outcome model, truths and calibration exactly once the caller builds the GBSG covariates. The legacy censoring step and some field names have no flex equivalent.
   - §5.4 has the details. Nothing was changed.
7. **Correction to Report 1 §6.** Its check export (without `-z`) did not copy two tracked files whose names need quoting, one under `dev/` and one under `quarto/`. Both are build-ignored, so that check was unaffected. This task's check used `git ls-files -z` and copied all 9,678 files (§3.5).
8. **Not rebuilt.** The pkgdown site and the articles were not rendered in this task (`vignettes/` and `quarto/` are read-only here).
9. **Scratch.** The equivalence scripts, temporary libraries and results are in `/tmp/seed_check_2026-11-11/`, the location the task names. The agents' working files are in this session's scratchpad. Neither is in the repository.

## Appendix — files changed since the start SHA

20 paths (A = added, M = modified), from `git diff --name-status 3fbbc55c` on the tree of the report commit. None is under `vignettes/` or `quarto/`.

| Status | Path |
|---|---|
| M | `DESCRIPTION` |
| M | `NEWS.md` |
| M | `R/bootstrap_analysis_dofuture.R` |
| M | `R/bootstrap_dofuture_main.R` |
| M | `R/calibrate_null_correction.R` |
| M | `R/compare_selection_rules.R` |
| M | `R/forestsearch_cross_validation.R` |
| M | `R/forestsearch_main.R` |
| M | `R/fpr_calibration.R` |
| M | `R/fs_fdr_report.R` |
| M | `R/mrct_case_study.R` |
| M | `R/mrct_simulation.R` |
| M | `R/mrct_trial_analysis.R` |
| M | `R/run_simulation_analysis.R` |
| M | `R/subgroup_consistency_main.R` |
| A | `dev/reports/REPORT_cran_prep_seed_legacy_2026-11-11.md` |
| A | `dev/reports/REPORT_cran_prep_seed_legacy_2026-11-11_detail.md` |
| A | `dev/reports/REPORT_cran_prep_seed_legacy_2026-11-11_quarto_callers.csv` |
| A | `dev/tasks/TASK_cran_prep_seed_legacy_2026-11-11.md` |
| A | `tests/testthat/test-rng-state-restored.R` |
