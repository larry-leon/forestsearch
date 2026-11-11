# REPORT — audit revisit: A-009, A-001, A-003, A-163 (read-only), 2026-11-11

Task: `dev/tasks/TASK_audit_revisit_2026-11-11.md` (committed as received, `aae7b683`).
Read-only. No fixes, diffs or recommendations appear in this report.

## 1. Provenance (Step 0)

| Item | Value |
|---|---|
| Host | `pop-os` |
| forestsearch branch / SHA at start | `master` / `02f88dac` (the handoff's expected SHA; `git status --porcelain` empty — gate passed) |
| forestsearch `DESCRIPTION` version | 0.4.1.9003 (as expected) |
| Task commit | `aae7b683 dev: add audit revisit task` (one file, `dev/tasks/TASK_audit_revisit_2026-11-11.md`) |
| `git log --oneline 02f88dac..HEAD` before the task commit | empty (HEAD was `02f88dac`) |
| Audit report | found at `master`: `dev/reports/REPORT_code_audit_2026-11-09.md` (10,299 lines); used as a pointer to sites and reproducers only |
| Audit base | `5b994177` (line numbers in the task); every site below is located by content and reported at `master` |

Sibling repositories (read only; nothing checked out, nothing changed):

| Repository | Path | Branch | SHA | Tree clean at Step 0 |
|---|---|---|---|---|
| fs-post-selection | `~/Documents/GitHub/fs-post-selection` | `main` | `5d165acc` | yes (`git status --porcelain` empty) |
| fs-glms-interpretable | `~/Documents/GitHub/fs-glms-interpretable` | `main` | `396144cd` | yes (`git status --porcelain` empty) |

`~/GitHub/fs-post-selection` and `~/GitHub/fs-glms-interpretable` do not exist; the clones above are the only ones found.

### How the work was done

- Each finding was traced by one reader from the current source, then handed to a second, independent reader whose brief was to refute it: re-open every cited `file:line` at `master`, re-derive the trace, and run its own `git grep` for exposures the first reader missed or over-counted. The sections below are the second reader's corrected text. Every correction it made is listed at the end of the section concerned (§3 "Verifier corrections").
- R executed (all `Rscript --vanilla`, sequential, nothing written to any repository):
  - A-003 reproducer `/tmp/audit_revisit_2026-11-11/A-003/repro_a003.R` (synthetic data, `pkgload::load_all()`): wall 2.22 s. This is the only reproducer; total reproducer time 2.22 s of the 10-minute cap.
  - A-003 read-only `readRDS()` scans of tracked `.rds` files under forestsearch `quarto/` (`scan_rds.R`, `scan_rds2.R`, 8.6 s each) plus three one-line structure checks (about 1 s each). These read payloads only and are not reproducers.
  - A-009, A-001, A-163: no R executed; reading settled each verdict.
- No test suite, `check()`, `document()`, `build()`, render, simulation or application re-run. No `fetch`/`pull`/`push`/`merge`.

## 2. Summary table

| ID | Site at master (file:line) | Claim (one line) | Holds at master | Exposed committed analyses | Would alter a reported number |
|---|---|---|---|---|---|
| A-009 | `R/bootstrap_dofuture_main.R:302-310` (`estimator_fn_boot <- make_effect_estimator(...)`; defaults from `R/glm_effect_estimators.R:107`, `:110`/`:113`) | On GLM outcomes, `forestsearch_bootstrap_dofuture()` rebuilds its effect estimator without the parent fit's `adverse_outcome` or `ps_adjust_method`. | **PARTLY.** Holds for the estimates the bias correction is built from (H_obs, Hc_obs, H_star, Hc_star, Hstar_obs, Hstar_star, Hcstar_*, and `FSsg_tab`): they always use `adverse_outcome = TRUE` and `ps_adjust_method = "none"`. Does not hold for the replicate's own search: `do.call(forestsearch, args_FS_boot)` carries the parent's resolved `adverse_outcome`, `ps_method` and `ps_adjust_method` via `args_call_all`. The audit's argument names are the source's. | **5**, all archived ACTG175 continuous applications, `quarto/applications/actg175/_archive/20260520_analysis_actg175_continuous{,_effMaxSG-both,_hr,_hrMaxSG-both,_hrMaxSG-pareto}.qmd` (`adverse_outcome = TRUE`, away from the continuous default FALSE; equal to the bootstrap's fallback, so no divergence). No committed non-survival analysis activates propensity adjustment. Separately, the resolved `adverse_outcome` differs at defaults (parent FALSE, bootstrap TRUE) in 8 call sites; the bootstrap runs in 2 of them. | **no** (the 2 running drivers re-orient the bootstrap output by the parent's `adverse_outcome`, and the MD is exactly sign-equivariant) |
| A-001 | `R/forestsearch_main.R:2766` (→ `R/forestsearch_helpers.R:1525-1526`, `m_diff`) and `R/forestsearch_main.R:3237-3238` (`m_diff_sel`); recorded/replayed value from `:2332-2334`, `:2436-2442`, `:2530-2537`, `:2548-2552` | Continuous outcome, DINA, `effect.threshold` and `hr.threshold` both unset: the fit's DINA floor is 1.25 while the returned object records 0 and bootstrap/CV replays resolve 0. | **YES**, on both routes (`subgroup_method = "dina"`; `use_dina = TRUE` with default `selected_only = TRUE`), in the both-unset case only. With either argument supplied, fit, record and replay agree. | **1** file, 3 calls: `quarto/dina/smoke-tests/dina_screening_smoketest.qmd:112`, `:126-127`, `:134-135` (forestsearch; `use_dina = TRUE`, both thresholds unset). Products are pass/fail checks from the fit; no bootstrap or CV. | **no** |
| A-003 | `R/run_simulation_analysis.R:1025`, `:1066` (`.extract_fs_estimates_gen`), `:1392`, `:1418` (`.extract_grf_estimates_gen`); found-case formulas `:1126-1129`, `:1530-1533` | Under an alternative, a replicate that identifies no subgroup stores specificity as \|Hc\|/N, not 1. | Per column, at all four sites: **sens AS DESIGNED** (0); **spec DEFECT** (stores \|Hc\|/N; own formula and reported definition give 1); **ppv AS DESIGNED** (stores 0; formula gives 0/0, reported convention 0); **npv DEFECT** (stores 1; own formula and reported definition give \|Hc\|/N; the stored 1 matches a non-detection convention documented at `R/fs_mr_oc_summary.R:56`). Null scenarios unaffected. | **126** forestsearch source files (125 under `quarto/`, plus `vignettes/articles/paper_simulations.qmd`), each averaging spec/npv over all replicates of an alternative via `build_classification_table()`, `format_oc_results()`, `summarize_single_analysis()`, a (Spec, NPV) figure or a direct mean; 17 have tracked renders. fs-post-selection and fs-glms-interpretable: none (different code paths). | **yes**: the H1 classification tables and the H1 (Spec, NPV) figure in the 17 tracked renders under `quarto/simulations/actg175/binary_methods/` (16) and `binary_legacy/` (1), e.g. `actg175_binary_m1_harm_methods_n700_demo.html`, table "grf-eff: Classification Under H1" (sens(Hc) 0.89 and ppv(Hc) 1.00 at any(H) 0.00). Per-cell magnitude from payloads: CANNOT DETERMINE (no per-replicate payload is tracked). |
| A-163 | `R/sim_aft_gbsg.R:638` (`b_cens <- -b0_cens * sigma`, in `.create_gbsg_dgm_`); draws at `R/sim_aft_gbsg.R:915` and `R/simulate_from_dgm.R:326-327` | Censoring-model coefficients are back-transformed with the event model's scale, not the censoring model's. | **INCONSISTENT.** Generated censoring covariate term is `z'gamma_cens * sigma / sigma_cens`; the fitted censoring model (and the coefficients stored in `cens_params$model$gamma`) is `z'gamma_cens`. No stated reason makes it correct as written. | Not searched (task: entry points only). Entry points: `create_gbsg_dgm`, `setup_gbsg_dgm`, `calibrate_k_inter`, `validate_k_inter_effect`, `fs_fdr_report` (exported), `get_dgm_with_output` (internal); draws via `simulate_from_gbsg_dgm`, `simulate_from_dgm`. | **unknown** (analyses not searched, per the task) |

## 3. Findings

### A-009 — GLM full bootstrap does not carry the parent's adverse_outcome / ps_adjust_method into its bias-correction estimator

Evidence is the source at master `aae7b683`. None of the traced files changed between the audit base `5b994177` and master, except `R/forestsearch_main.R`, whose only change (`:3719-3737`, the df.est flag attachment) is unrelated to this finding. The audit cites `bootstrap_dofuture_main.R:302`; the line is the same at master. The audit's argument names (`adverse_outcome`, `ps_adjust_method`) are the source's names. `ps_method` is the argument that turns propensity adjustment on, and the estimator never reads it.

#### The bootstrap's estimator and its defaults

`R/bootstrap_dofuture_main.R:292-310`:
```r
  outcome_type <- args_forestsearch_call$outcome_type
  if (is.null(outcome_type)) outcome_type <- "survival"
  if (length(outcome_type) > 1L) outcome_type <- outcome_type[1L]
  is_glm <- outcome_type != "survival"

  estimator_fn_boot <- NULL
  cox.formula.boot  <- NULL

  if (is_glm) {
    # GLM path: build the estimator closure from saved arguments
    estimator_fn_boot <- make_effect_estimator(
      outcome_type      = outcome_type,
      treat.name        = args_forestsearch_call$treat.name,
      outcome.name      = args_forestsearch_call$outcome.name,
      event.name        = NULL,
      offset.name       = args_forestsearch_call$offset.name,
      effect_measure    = args_forestsearch_call$effect_measure,
      adjust_covariates = args_forestsearch_call$adjust_covariates
    )
```
`args_forestsearch_call` is `fs.est$args_call_all` (`:267`). It holds `adverse_outcome`, `ps_method` and `ps_adjust_method` (see (1)), but none of them is passed here. Missing arguments fall back to the formals in `R/glm_effect_estimators.R:100-113`:
```r
make_effect_estimator <- function(
    ...
    adverse_outcome    = TRUE,
    robust_se          = TRUE,
    adjust_covariates  = NULL,
    ps_adjust_method   = c("none", "iptw", "dr_gcomp"),
    ...
) {
  ps_adjust_method <- match.arg(ps_adjust_method)
```
**Fallback values:**
- `adverse_outcome = TRUE` (`:107`).
- `ps_adjust_method = "none"`: `match.arg` returns the first choice (`:110`, `:113`).

There is no NULL path here: the argument is simply missing, so the formal default `TRUE` applies.

**What `adverse_outcome` changes in the estimate:**
- Continuous: `.make_lm_estimator`, `R/glm_effect_estimators.R:945-947`:
  ```r
  if (!adverse_outcome) {
    data_slice[[outcome.name]] <- -data_slice[[outcome.name]]
  }
  ```
  The fast closure does the same at `:1082-1084`. FALSE negates Y. The MD changes sign exactly and its SE is unchanged (`lm` + `vcov`, `:963-968`).
- Binary: `.make_glm_binary_estimator`, `:390-392`, and the fast OR closure at `:1189-1191`: `data_slice[[outcome.name]] <- 1L - data_slice[[outcome.name]]`. FALSE swaps the outcome to 1 − Y:
  - OR: the log OR changes sign exactly.
  - RD: the RD changes sign exactly.
  - RR: the result is the ratio of 1 − Y risks, which is not a sign change.
- Count, or binary with IRR/IRD: the count branch (`:198-210`) and the Poisson estimator (`:726-928`) never take `adverse_outcome` (0 occurrences in `:726-928`). It has no effect there.

**What `ps_adjust_method` changes:** with `"iptw"`, the closure fits with `weights = sw` when the slice has an `sw` column (`:956`, `:963`). With `"none"` it fits unweighted (`:965`).

#### (1) The replicate's own search: carries the parent's resolved values

The parent resolves `adverse_outcome` before it captures its arguments. `R/forestsearch_main.R:1797-1799`:
```r
  if (is.null(adverse_outcome)) {
    adverse_outcome <- (outcome_type %in% c("binary", "count"))
  }
```
`:1880-1881`:
```r
  args_names <- names(formals())
  args_call_all <- mget(args_names, envir = environment())
```
So `args_call_all$adverse_outcome` holds the resolved logical:
- FALSE for continuous when unset;
- TRUE for binary and count when unset.

`ps_method` and `ps_adjust_method` are written back after they are resolved. `:2596-2599`:
```r
  ps_method_supplied <- !is.null(ps_method)
  if (is.null(ps_method)) {
    ps_method <- if (is.RCT) "none" else "grf"
  }
```
`:2632-2641`:
```r
  args_call_all <- .sync_args_call_all(args_call_all, environment(),
                                       "ps_method")
  ...
  ps_adjust_method <- match.arg(ps_adjust_method)
  ps_adjust_resolved <- if (ps_method == "none") "none" else ps_adjust_method
  args_call_all$ps_adjust_method <- ps_adjust_resolved
```
The returned object carries this list (`:3968`, `args_call_all = args_call_all`). No later `.sync_args_call_all` call touches these names (`:2149` stop_threshold, `:2461` dmin.grf, `:2535` thresholds, `:2716` adjust_covariates, `:3307` event.name). The only other direct write is `:1984` (`sg_focus`).

The replicate clones that list. `R/bootstrap_analysis_dofuture.R:406`: `args_FS_template <- fs.est$args_call_all`. `:558-560`:
```r
    args_FS_boot <- args_FS_template
    args_FS_boot$df.analysis <- dfnew_boot
    args_FS_boot$df.predict <- dfnew_predict
```
The overrides that follow touch only output, GRF/DINA caches, `ps_hat` (`:592`, `args_FS_boot$ps_hat <- NULL`), the parallel plan and `mr_inference`. Then `:614`: `run_bootstrap <- try(do.call(forestsearch, args_FS_boot), TRUE)`.

Inside the replicate's `forestsearch()`, the non-NULL `adverse_outcome` is kept. It is used to build the search closure (`forestsearch_main.R:2220-2228` `adverse_outcome = adverse_outcome`; fast closures `:2244-2257`), which reaches the search (`:3469-3470`) and consistency (`:3598-3599`).

**Verdict for (1):**
- `adverse_outcome`: the parent's resolved value.
- `ps_method`: the parent's resolved value. The replicate re-estimates the score itself because `ps_hat` is nulled.
- `ps_adjust_method`: the parent's resolved value.

The claim does not hold for (1).

A side note on the parent itself. When `ps_method != "none"`, the parent replaces its own search closure at `forestsearch_main.R:2703-2712`, and that rebuild also omits `adverse_outcome`:
```r
    if (is_glm) {
      estimator_fn <- make_effect_estimator(
        outcome_type     = outcome_type,
        treat.name       = treat.name,
        outcome.name     = outcome.name,
        offset.name      = offset.name,
        effect_measure   = effect_measure,
        ps_adjust_method = ps_adjust_resolved
      )
```
So under active PS, the parent's search closure uses `adverse_outcome = TRUE`, and the replicate repeats this because it runs the same code. No committed GLM analysis activates PS (see the exposure section).

#### (2) The estimates the bias correction is built from: fall back to defaults

Every constituent comes from `estimator_fn_boot` through `fit_subgroup_effect()`, which only calls the closure (`R/Cox_estimation_helpers.R:214-225`; `fit_effect_models` at `:239-248`).

**Parent's subgroup on the original data (H_obs, Hc_obs).** `bootstrap_dofuture_main.R:328`:
```r
    effect_fits <- fit_effect_models(fs.est$df.est, estimator_fn_boot)
```

**Passed into the loop.** `:478-490` call `bootstrap_results(... estimator_fn = estimator_fn_boot, ...)`, and `bootstrap_analysis_dofuture.R:444` sets `is_glm_boot <- !is.null(estimator_fn)`.

**Parent's subgroup on the bootstrap sample (H_star, Hc_star).** `:485` `fitH_star <- fit_subgroup_effect(df_H, estimator_fn)` and `:514` `fitHc_star <- fit_subgroup_effect(df_Hc, estimator_fn)`.

**Bootstrap-selected subgroup on the original data (Hstar_obs).** `:651` `df_Hstar <- subset(df_PredBoot, treat.recommend == 0)` and `:695` `fitHstar_obs <- fit_subgroup_effect(df_Hstar, estimator_fn)`.

**Bootstrap-selected subgroup on the bootstrap sample (Hstar_star).** `:708-709`:
```r
        fitHstar_star <- fit_subgroup_effect(
          subset(dfboot_PredBoot, treat.recommend == 0), estimator_fn
```

**Complement versions.** `:747` and `:760-761`.

**Bias corrections.** `:724-728`:
```r
      if (all(is.finite(c(H_obs, H_star, Hstar_star, Hstar_obs)))) {
        # Method 1: Simple optimism correction
        H_biasadj_1 <- H_obs - (Hstar_star - Hstar_obs)
        # Method 2: Double correction
        H_biasadj_2 <- 2 * H_obs - (H_star + Hstar_star - Hstar_obs)
```

**Reported table.** The GLM summary table `FSsg_tab` uses the same closure (`bootstrap_dofuture_main.R:707-712`, `estimator_fn = estimator_fn_boot`).

**Verdict for (2):**
- `adverse_outcome` = TRUE regardless of the parent.
- `ps_adjust_method` = "none" regardless of the parent.

The claim holds for (2). The parent and the bootstrap differ in two cases:
- the parent's resolved `adverse_outcome` is FALSE: a continuous fit with it unset, or a binary or continuous fit with it set to FALSE;
- the parent has PS active (`ps_adjust_resolved` is "iptw" or "dr_gcomp").

For count outcomes, `adverse_outcome` makes no difference because the estimator ignores it.

The rest of the pipeline is linear and sign-symmetric on identity-scale measures:
- `ci_est` with `est.loghr = FALSE` is `x ± z·sd` (`R/bootstrap_calculations_helpers.R:151-200`);
- `get_dfRes` (`:253-330`) feeds into it the IJ target from `get_targetEst`, whose mean negates and whose `cov_i^2` variance is unchanged under negation.

So for MD (and for RD/OR on the log scale), a parent FALSE against a bootstrap TRUE produces exactly negated H0/H1/H2 and bounds with unchanged SDs.

`mr_estimates_table.R` and the cross-validation sites are not called by `forestsearch_bootstrap_dofuture()`; it only mentions `mr_estimates_table` in its roxygen (`:198`). They are skipped.

#### Exposure: every committed non-survival call of forestsearch_bootstrap_dofuture()

**Method.**
1. `git ls-files` / `git grep` over forestsearch `quarto/` and `vignettes/`, fs-post-selection `manuscript_rev1/` and `manuscript/`, and fs-glms-interpretable `quarto/` and `manuscript/`. 446 tracked `.qmd`/`.R`/`.Rmd` files and 3 tracked `.rmarkdown` render intermediates mention the function, all in forestsearch.
2. The two sibling repos have no tracked call; they only mention it in `.md` reports.
3. No package function wraps the call (`grep` of `R/` finds only comments).
4. Files were classified by `outcome_type`. For the binary simulation drivers, the sourced tracked helpers `quarto/simulations/actg175/binary_{legacy,methods}/_fsparams/fs*_params.R` were read: `outcome_type = "binary"`, `effect_measure = "OR"`, `is.RCT = TRUE`, and no adverse/ps arguments.
5. Non-flagged files that mention GLM words were checked; they are survival (for example `quarto/simulations/actg175/binary_methods/smoke_test_adjust_covariates.R:134` bootstraps a Cox fit; its `is.RCT = FALSE, ps_method = "logistic"` at `:116` is a survival guard-test fit that is not bootstrapped).
6. No non-survival file sets `is.RCT = FALSE`, `ps_method`, a non-"none" `ps_adjust_method` or `ps_hat`. `fs_is_RCT` and `is_rct` are always `TRUE`.

**Column meanings.**
- "Away from default" applies the task's literal definition.
- "Resolved values differ" compares the parent's resolved values with the bootstrap's TRUE / "none".

| # | Call site(s) (forestsearch repo) | outcome_type / measure | Parent adverse_outcome (as written → resolved) | is.RCT / ps_method / ps_adjust_method | Away from default? | Resolved values differ? |
|---|---|---|---|---|---|---|
| 1 | `quarto/applications/actg175/_archive/20260520_analysis_actg175_continuous.qmd:560`; `_continuous_hr.qmd:560`; `_continuous_hrMaxSG-both.qmd:560`; `_continuous_hrMaxSG-pareto.qmd:560` (fit `fs_result` at `:338`, `adverse_outcome = TRUE` at `:348`, `fs_outcome_type <- "continuous"` at `:148`, `fs_is_RCT <- TRUE` at `:168`); `_continuous_effMaxSG-both.qmd:584` (fit `:358`, `:368`, `:149`, `:170`) | continuous / MD (outcome `y_decline`) | TRUE → TRUE | TRUE / unset → "none" / unset → "none" | **yes** (adverse_outcome TRUE ≠ continuous default FALSE) | no (TRUE = TRUE) |
| 2 | `quarto/simulations/actg175/continuous/sim_fs_maxeffCons_fb_mr_md40_knoise0_n500_batch_1_100.qmd:670` (nb_boots 300L `:112`, guard `if (run_fb)` `:659`, fit `:550-565`, `adverse_outcome <- FALSE` `:216`, `is_rct <- TRUE` `:215`) — **runs** | continuous / MD | FALSE → FALSE | TRUE / unset / unset → "none" | no (FALSE = default) | **yes** (FALSE vs TRUE) |
| 3 | `.../sim_fs_maxeffCons_fb_mr_md40_knoise0_n500_batch_1_20.qmd:428` (nb_boots 20L `:90`, guard `:417`, `adverse_outcome = adverse_outcome` `:339`, `<- FALSE` `:150`, `is_rct` `:149`) — **runs** | continuous / MD | FALSE → FALSE | TRUE / unset / unset | no | **yes** |
| 4 | `.../sim_fs_maxeffCons_fb_mr_md40_knoise0_n500_batch_1_50.qmd:878` (nb_boots 0L `:275`, guard `:867`; `:379`, `:378`) — skipped (nb_boots 0) | continuous / MD | FALSE → FALSE | TRUE / unset / unset | no | yes, but the call is skipped |
| 5 | `.../sim_fs_maxeffCons_mr_md120_knoise0_n500_batch_1_1000.qmd:699`; `.../sim_fs_maxeffCons_mr_md40_knoise0_n500_batch_1_100.qmd:670`; `.../sim_fs_maxeffCons_mr_md40_knoise0_n500_batch_1_1000.qmd:681`; `.../sim_fs_maxeffCons_mr_md40_knoise0_n700_batch_1_1000.qmd:681` (all nb_boots NULL `:112`, guarded by `if (run_fb)`; `adverse_outcome <- FALSE`, `is_rct <- TRUE`) — skipped | continuous / MD | FALSE → FALSE | TRUE / unset / unset | no | yes, but the call is skipped |
| 6 | `quarto/resampling/test_continuous_glm_compatibility.qmd:539` (fit `fs_args` `:218-252`, no adverse_outcome) | continuous / MD | unset → FALSE | TRUE / unset / unset | no | **yes** |
| 7 | Current ACTG175 binary applications: `analysis_actg175_binary_multimethod_fixed_family.qmd:773,1282,1505` (fits `:557`/`:1221`/`:1443`, `adverse_outcome = TRUE` `:570`/`:1234`/`:1456`, `is.RCT = TRUE`); `_frontend.qmd:766,1275,1498`; `_psi_v2_2.qmd:751,1260,1483`; `_psi_v3a.qmd:751,1260,1483`; `analysis_actg175_binary_sgfocus.qmd:890` (`:662` TRUE); archive `20260730_*_psi_v2_2A{,_mac}.qmd:743,1253,1476`, `_psi_v2_2_mac.qmd:743,1253,1476`, `_psi_v2_2_mac_w2.qmd:746,1256,1479`; `20260520_analysis_actg175_binary_effMaxSG-both.qmd:613` (`:391` TRUE) | binary / OR | TRUE → TRUE | TRUE / unset / unset | no | no |
| 8 | `vignettes/articles/actg175-glm-outcomes.qmd:522` (`fits_s3` `:475` = binary fits `:421`/`:441`/`:460`, `adverse_outcome = TRUE`, `is.RCT = TRUE`) | binary / OR | TRUE → TRUE | TRUE / unset / unset | no | no |
| 9 | Binary MR/FB sweeps: `quarto/simulations/actg175/binary/{mr_coverage_sweep_or075,_s100,_s1k,_s500,_or20,_or35_s1k,effMaxSG_…_or075,maxcons_…_or075}.qmd:496`, `maxeffCons_mr_coverage_sweep_or075.qmd:490`, `mr_coverage_sweep_or10/_or15/_or15_effMaxSG/_or35_effMaxSG.qmd:467`, `binary_020/maxeffCons_mr_coverage_sweep_or075.qmd:605`, and the tracked render intermediate `binary/mr_coverage_sweep_or20.rmarkdown:496` (e.g. `_s1k`: `adverse_outcome <- TRUE` `:154`, passed `:422`, `is.RCT = TRUE` `:418`; the `adverse_outcome = FALSE` at `:227`/`:257` is in `calibrate_glm_interaction()` DGM calls, not the fit) | binary / OR | TRUE → TRUE | TRUE / unset / unset | no | no |
| 10 | 93 binary legacy/methods driver documents: 72 in `quarto/simulations/actg175/binary_legacy/*.qmd` (calls at `:323,:434`, or `:322,:460` / `:351,:574`) and 21 in `binary_methods/*.qmd` (e.g. `actg175_binary_m1_harm_methods_n1000_full.qmd:654,886`); args from `_fsparams/fs*_params.R:73-74,109` | binary / OR | unset → TRUE (a few `binary_methods` docs pass TRUE) | TRUE / unset / unset | no | no |
| 11 | `quarto/resampling/test_glm_compatibility.qmd:499`; `quarto/simulations/actg175/binary_methods/verify_glm_adjust_layers.R:114` (fit `fs_adj` from `fs_common` `:59-79` with `adjust_covariates = "S"` `:99`; no adverse_outcome, no is.RCT) | binary / OR | unset → TRUE | TRUE (unset → formal default TRUE) / unset / unset | no | no |
| 12 | Count: `quarto/applications/count_data_demo.qmd:901` (fit `:653`, `is.RCT = TRUE` `:674`); `quarto/resampling/test_count_glm_compatibility.qmd:513`; `quarto/applications/gbsg/_archive/2026-05-27_gbsg_survival_cox_poisson.qmd:1059` and its tracked render intermediate `.rmarkdown:1059` (`adverse_outcome = NULL`, `ps_method = NULL`, `ps_adjust_method = "none"`, `is.RCT = TRUE`); `quarto/applications/gbsg/_broken/gbsg_analysis_cox_poisson{,_1,_check,_checking-bootstrap,_checking-bootstrap_new,_checking-bootstrap_new-500,_maxSG}.qmd` (`fs_pois_bc` calls at `:772`, `:778`, `:772`, `:636`, `:638`, `:638`, `:638`; is.RCT unset → TRUE) | count / IRR | unset/NULL → TRUE; ignored by the count estimator | TRUE / unset / "none" | no | no (adverse_outcome has no effect on count) |
| — | Not non-survival calls: the `fs_cox_bc` survival bootstraps in the same gbsg files (`2026-05-27_gbsg_survival_cox_poisson.qmd:475` and `.rmarkdown:475`; `_broken/gbsg_analysis_cox_poisson*.qmd:375/:381/:375/:375/:377/:377/:377`, all on `fs_cox`, outcome_type unset → survival); `quarto/guides/vscode_claudecode_forestsearch.qmd:680,875` (document-level `eval: false`, `:23`); `quarto/guides/count_data_hte_summary.qmd:337`, `count_data_demo.qmd:177,893`, `test_glm_compatibility.qmd:59` (prose); `vignettes/articles/forestsearch-2024.Rmd:478` (survival) | — | — | — | — | — |

**Summary of exposure:**
- Exposed under the literal definition: row 1, 5 archived continuous application call sites. Their `adverse_outcome = TRUE` equals the bootstrap's fallback, so their bootstrap numbers are on the parent's scale.
- No committed non-survival analysis turns on propensity adjustment, so the `ps_adjust_method` half of the claim reaches no committed analysis.
- At defaults, the resolved values differ in:
  - rows 2-3, which run the bootstrap;
  - rows 4-5, where the call exists but is skipped;
  - row 6.

#### Would it alter a reported number

**No.**

Rows 2-3 re-orient the bootstrap's output by the parent's `adverse_outcome`. `sim_fs_maxeffCons_fb_mr_md40_knoise0_n500_batch_1_100.qmd:713-721`:
```r
    if (!is.null(fb) && !is.null(fb$H_estimates)) {
      he <- fb$H_estimates
      .orient <- if (isTRUE(adverse_outcome)) 1 else -1
      rec$fb_H_est <- .orient * (he$H2 %||% NA_real_)
      fb_ci <- sort(.orient * c(he$H2_lower %||% NA_real_,
                                he$H2_upper %||% NA_real_), na.last = TRUE)
      rec$fb_H_lo  <- fb_ci[1]
      rec$fb_H_hi  <- fb_ci[2]
      rec$fb_H_se  <- he$sdH2 %||% NA_real_   # scale-free, NOT negated
```
batch_1_100 does the same for the Hc arm at `:771-779`. batch_1_20 does the same for the H arm at `:471-479` and records no Hc arm (it never reads `Hc_estimates`).

Why the recorded numbers are unaffected:
- For MD, the FALSE closure gives exactly the negative of the TRUE closure's estimate with the same SE (`glm_effect_estimators.R:945-947`, `:963-968`).
- H1/H2 are linear in the constituents, and `ci_est` is symmetric on the identity scale.
- So the recorded `fb_H_*` (and, in batch_1_100, `fb_Hc_*`) values equal what a bootstrap built with the parent's FALSE would give.
- The decomposition columns `fb_src1` and `fb_src2` are documented and tabled as raw-scale (`:723-727`; table captions at `:1304`/`:1306` "RAW scale").

Row 6 reports no bootstrap estimate: its checks are structural, and `print(fs_bc)` is captured, not shown (`:777`, `:790`). Every other row has no difference between the parent's values and the bootstrap's.

#### Reproducers

None run. Reading the source settled each verdict. The audit's V047 reproducer, which shows the sign relationship for the closure alone, was not re-run.

#### Verifier corrections

- *wrong_line*: The template-clone block is at bootstrap_analysis_dofuture.R:558-560, not :557-559.
- *missed_exposure*: A tracked non-survival call site is missing from the table: quarto/simulations/actg175/binary_methods/verify_glm_adjust_layers.R:114. It bootstraps fs_adj, a binary/OR fit with adverse_outcome unset (resolves to TRUE), is.RCT unset (formal default TRUE), no ps_method or ps_adjust_method, and adjust_covariates = "S". It is not exposed and its resolved values do not differ. Two tracked .rmarkdown render intermediates also hold non-survival calls and were outside the investigator's 446 .qmd/.R/.Rmd file set: quarto/simulations/actg175/binary/mr_coverage_sweep_or20.rmarkdown:496 (binary/OR, adverse_outcome <- TRUE at :154, is.RCT = TRUE at :418) and quarto/applications/gbsg/_archive/2026-05-27_gbsg_survival_cox_poisson.rmarkdown:1059 (count; adverse_outcome = NULL, ps_method = NULL, ps_adjust_method = "none"). Neither is exposed and neither has a difference. None of this changes the verdict.
- *false_exposure*: Row 12 lists survival call sites as non-survival ones. In each gbsg cox_poisson file, the first listed line (fs_cox_bc: archive :475; _broken :375/:381/:375/:375/:377/:377/:377) bootstraps fs_cox, a survival fit, so it is not a non-survival call site. Only the fs_pois_bc lines (:1059; :772/:778/:772/:636/:638/:638/:638) are count calls. In count_data_demo.qmd, :893 is prose; the only call is :901. Neither change affects the verdict, since every row 12 site has no difference anyway.
- *other*: Row 10's count is wrong. binary_legacy holds 72 driver .qmd files that call the function and binary_methods holds 21 more, so the row covers 93 documents, not 72. All are binary/OR with no adverse_outcome FALSE and no PS setting. Not exposed, no difference.
- *trace_error*: The would-alter answer says both running drivers negate 'the reported H2, its bounds and the Hc arm'. Only batch_1_100 records an Hc arm (fb_Hc_*). batch_1_20 never reads Hc_estimates and records only fb_H_*, fb_src1/2 and fb_nres. Two quote line ranges are also off by one: the batch_1_100 H-arm block quoted as :713-720 runs to :721.


### A-001: the continuous-outcome DINA floor comes from the raw hr.threshold when both thresholds are unset

#### Did the sites change after the audit?

`git log --oneline 5b994177..HEAD -- R/forestsearch_main.R R/forestsearch_helpers.R` returns one commit, `808e284d fix: bootstrap Eq. 7 terms use row-wise subgroup flags (A-002)`. That commit touches only `R/forestsearch_main.R` around old line 3719, which is the `df.est_out` / `df.predict_out` merge. `git diff 5b994177 HEAD --stat -- R/forestsearch_helpers.R R/dina_subgroup.R R/bootstrap_analysis_dofuture.R R/forestsearch_cross_validation.R` is empty. At base, `git show 5b994177:R/forestsearch_main.R | grep -n m_diff_sel` gives 3237, and the `.forestsearch_dina_select` argument line is at 2766. Both line numbers are the same at master.

#### The two threshold arguments, as the source names them

Both are spellings of the screening threshold c1: `effect.threshold` (formal default `NULL`, `R/forestsearch_main.R:1668`) and `hr.threshold` (formal default `1.25`, `:1670`). Below, "the first argument" means `effect.threshold` and "the second argument" means `hr.threshold`. Neither DINA floor reads the consistency pair (`consistency.threshold` / `hr.consistency`).

#### How the value is resolved (R/forestsearch_main.R at master)

`R/forestsearch_main.R:1805-1811`:
```r
  user_set_threshold   <- !is.null(effect.threshold)      || !missing(hr.threshold)
  ...
  .fs_c1_legacy <- if (missing(hr.threshold))   NULL else hr.threshold
  ...
  if (!is.null(effect.threshold))      hr.threshold   <- effect.threshold
```
Line 1811 is the only assignment to the local `hr.threshold` in the function. A grep for `hr\.threshold\s*(<-|=)` matches :1670 (the formal) and :1811. It also matches argument positions at :2187, :2550, :2766 and :3452, the comparisons `hr.threshold == 1.25` at :2281 and :2332, the list element `consistency_overrides$hr.threshold <- effect_threshold` at :3604, a message string at :114 and roxygen at :1463. None of these other matches assigns the local. With both arguments unset, `hr.threshold` stays 1.25 for the whole fit.

For continuous outcomes the default effect measure is MD (`:1790-1792`: `continuous = "MD"`). The default `adverse_outcome` is FALSE (`:1797-1798`).

`R/forestsearch_main.R:2275` and `:2326-2334` (inside `if (outcome_type != "survival")`, `:2208`):
```r
    effect_threshold <- hr.threshold
    ...
    } else if (effect_measure == "MD") {
      ...
      if (!user_set_threshold && hr.threshold == 1.25) {
        # Survival default 1.25 is meaningless for MD; use 0.0
        effect_threshold <- 0.0
      }
```
Only the local `effect_threshold` is changed to 0. `hr.threshold` is not.

#### What the returned object records

`R/forestsearch_main.R:2436-2442`:
```r
    threshold_config <- list(
      outcome_type = outcome_type,
      effect_measure = effect_measure,
      screening = effect_threshold,
      ...
      screening_natural = if (is_identity) effect_threshold
```

`R/forestsearch_main.R:2530-2537`:
```r
  .fs_identity_scale <- outcome_type != "survival" &&
    effect_measure %in% c("RD", "IRD", "MD")
  effect.threshold <- if (.fs_identity_scale) effect_threshold else hr.threshold
  ...
  args_call_all <- .sync_args_call_all(args_call_all, environment(),
                                       c("effect.threshold",
                                         "consistency.threshold"))
```
`args_call_all` is first captured from the formals at `:1880-1881` (`args_call_all <- mget(args_names, envir = environment())`). This happens after :1811, so `args_call_all$hr.threshold` is 1.25 in the unset case. The sync at :2535-2537 then overwrites only `effect.threshold`, setting it to 0. No other `.sync_args_call_all` call (:2149, :2461, :2632, :2716, :3307) touches a threshold.

`R/forestsearch_main.R:2548-2552`:
```r
  admission_resolved <- .fs_resolve_admission(
    sg_focus, subgroup_method,
    hr.threshold           = threshold_config$screening,
```
`R/forestsearch_helpers.R:2482` sets `effect_floor <- as.numeric(hr.threshold)`, which here is 0. Under `subgroup_method = "dina"` the effect floor always applies (`:2419-2420`: `if (method %in% c("dina", "grf")) return(c(effect = TRUE, consistency = FALSE))`). On the consistency engine (the `use_dina` route) it applies unless `sg_focus = "maxeff"` (`:2413-2414`), in which case `effect_floor` is NULL.

Both return lists carry these values. The DINA route sets `args_call_all = args_call_all` at `R/forestsearch_main.R:2807`, and `threshold_config` and `admission = admission_resolved` at :2811-2815. The main route does the same at :3968-3982.

#### Site 1: the `subgroup_method = "dina"` route

`R/forestsearch_main.R:2761-2776`:
```r
  if (subgroup_method == "dina") {
    dsel <- .forestsearch_dina_select(
      ...
      outcome_type = outcome_type, hr.threshold = hr.threshold, n.min = n.min,
      ...
      admission = admission_resolved)
```

`R/forestsearch_helpers.R:1521-1526`:
```r
  # Harm floor on the link scale: log(hr.threshold) for ratio families
  # (cox/binomial/poisson), identity (mean difference) for gaussian.
  # Derived BEFORE the fit so the identity-scale refusal costs no model run.
  .dina_assert_ratio_estimand(da$fit$family, effect_measure)
  m_diff <- if (identical(da$fit$family, "gaussian")) hr.threshold
            else log(hr.threshold)
```
The family is gaussian for a continuous outcome (`R/forestsearch_helpers.R:1015`). `.dina_assert_ratio_estimand` refuses only RD/IRD on a non-gaussian family (`:1446-1462`), so MD passes. At `R/forestsearch_helpers.R:1590-1594`, `dina_subgroup(..., m_diff = m_diff, ...)` receives this value. `R/dina_subgroup.R:749` and `:888` then apply `if (mean_tau < m_diff) next`. The qualifying family is therefore every candidate whose oriented mean tau-hat is at least 1.25, with the orientation set by `tau_sign` at helpers :1539. The effect-statistic re-ranking then applies `admission$effect_floor` (0), but only inside that family (`R/forestsearch_helpers.R:1255-1260`: `floor_cmp <- if (is.null(admission)) NULL else admission$effect_floor`).

The multiplier-resampling family on this route is `dsel$candidates`, the family built under the fit's floor (`R/forestsearch_main.R:2858-2859`: `.mr_fam <- .fs_mr_family_from_table(.mr_df, dsel$candidates, ...)`). Multiplier resampling does not re-run `forestsearch()`.

#### Site 2: the `use_dina` selected-cut screening route

`R/forestsearch_main.R:3236-3251`:
```r
        .dina_assert_ratio_estimand(da$fit$family, effect_measure)
        m_diff_sel <- if (identical(da$fit$family, "gaussian")) hr.threshold
                      else log(hr.threshold)
        sgsel <- dina_subgroup(
          fit                 = dina_res,
          df                  = df.dina,
          covariates          = confounders.name,
          m_diff              = m_diff_sel,
          n_min               = n.min,
          max_depth           = da$select$max_depth,
          grid_probs          = da$select$grid_probs,
          sg_focus            = sg_focus,
          selection_rule      = selection_rule,
          effect_neighborhood = effect_neighborhood
        )
```
This branch runs when `use_dina && is.null(dina_cuts)` (`:3180`) and `isTRUE(da$selected_only)` (`:3225`). `selected_only` defaults to TRUE (`R/forestsearch_helpers.R:1112`: `selected_only <- get_arg("selected_only", TRUE)`). The consistency search itself screens on the local `effect_threshold`, which is 0. The search overrides use it at `R/forestsearch_main.R:3452` (`hr.threshold = if (!is.null(effect_threshold)) effect_threshold else hr.threshold,`), and the consistency stage uses it at :3604 (`consistency_overrides$hr.threshold <- effect_threshold`). So on this route only the DINA cut added to the pool uses the 1.25 floor.

#### What a replay resolves

The bootstrap copies the parent's argument list at `R/bootstrap_analysis_dofuture.R:406` (`args_FS_template <- fs.est$args_call_all`) and :558. At :559-607 it changes the data, output, re-selection and parallel settings, but no threshold. It then calls `do.call(forestsearch, args_FS_boot)` at :614. Cross-validation does the same at `R/forestsearch_cross_validation.R:345` and :387 (`cv_args <- fs_args`) and calls `do.call(forestsearch, cv_args)` at :479. The second CV function does the same at :859, :877 and :1001.

In the replay, `effect.threshold = 0` is non-NULL, so line :1811 sets `hr.threshold <- 0`. Then `m_diff = 0` (helpers :1525) and `m_diff_sel = 0` (main :3237). The replay also supplies `hr.threshold = 1.25`, so `.fs_c1_legacy` is 1.25. The mismatch check in `.fs_resolve_threshold_pair` (the `stop()` at `:113-117`) is never reached. The function returns at `:96-99` on every path except the consistency engine on a ratio estimand (`ratio_pair <- identical(..., "consistency") && (survival || binary OR)`; `if (!ratio_pair) return(c2)`). No error or warning is raised.

The bootstrap replays only when the parent found a subgroup (`R/bootstrap_dofuture_main.R:240-245`: "No identified subgroup => nothing to bias-correct").

#### Floor used by the fit and by a replay (continuous / MD, gaussian DINA family)

| Case | Route | Fit: m_diff (DINA floor) | Recorded: threshold_config$screening / admission$effect_floor / args_call_all$effect.threshold | Replay (bootstrap, CV): m_diff |
|---|---|---|---|---|
| Both unset (`effect.threshold = NULL`, `hr.threshold` missing) | `subgroup_method = "dina"` | **1.25** | 0 / 0 / 0 (`args_call_all$hr.threshold` = 1.25) | **0** |
| Both unset | `use_dina = TRUE`, `selected_only = TRUE` | **1.25** (DINA cut only; the search screen uses 0) | 0 / 0, or NULL if `sg_focus = "maxeff"` / 0 | **0** |
| Only `effect.threshold = e` | either route | e (:1811 sets `hr.threshold <- e`) | e / e / e (`args_call_all$hr.threshold` = e, captured after :1811) | e |
| Only `hr.threshold = h`, including h = 1.25 written explicitly | either route | h (`user_set_threshold` is TRUE, so :2332 does not remap) | h / h / h | h |

The two routes do not differ in floor value. Under `selected_only = FALSE` the `use_dina` route takes `fr$cut_expr` (`:3262`) and does not use `m_diff_sel`.

#### The package's own specification (pointer only)

The roxygen at `R/forestsearch_main.R:446-449` says that under `subgroup_method = "dina"`, and under `selected_only = TRUE` screening, m_diff is "computed as log(hr.threshold) -- identity for gaussian -- from the resolved effect.threshold". The resolved-defaults table at `:1268` gives `MD \tab 0.0 \tab 0.0`. `R/dina_subgroup.R:146` repeats "from the resolved `effect.threshold`". In the unset case the code takes the floor from the unresolved `hr.threshold`, which is 1.25. Which floor is correct is a method choice and is not judged here.

#### Same site as the earlier "use_dina unoriented floor" review?

Same `dina_subgroup()` call on the `use_dina` route, different defect: the earlier review is about the missing `tau_sign` (orientation), A-001 is about the floor's value, and A-001 also has a second site (`subgroup_method = "dina"`) that the earlier review does not cover.

Supporting detail: `dev/reports/REPORT_directive_C_2026-09-18.md:148-152` and `dev/reports/STATUS_thresholds_workstream_2026-09-18.md:74` placed the unoriented floor at `forestsearch_main.R:3054-3063`. At that time (c8b475b4), `m_diff_sel` is at :3051 and the `dina_subgroup()` call is at :3053-3063. At master this call is at `:3237-3251`, and it still omits `tau_sign`, which affects orientation. A-001 concerns the floor's value (the raw `hr.threshold` rather than the resolved `effect_threshold`). Its second site is at :2766 / helpers :1525. At base 5b994177, lines 3040-3075 hold the GRF route's multiplier-resampling block, not this call; the line numbers had drifted.

#### Exposure (tracked files in the scoped directories)

Search method: `git grep` in forestsearch `quarto/` and `vignettes/` (.qmd/.R/.Rmd), for files that contain both a quoted `"continuous"` and `dina`, and also for `outcome_type = "continuous"`. The `_archive` directories were included and gave no hits. Also searched: fs-post-selection `manuscript/` and `manuscript_rev1/` (including `payload_manifest.qmd`), and fs-glms-interpretable `quarto/` and `manuscript/`. No package wrapper defaults to DINA: a grep for `use_dina\s*=\s*TRUE|subgroup_method\s*=\s*"dina"` in `R/` finds only comments. Of the other `do.call(forestsearch, ...)` sites in `R/`, `fs_fdr_report()` refuses continuous outcomes (`R/fs_fdr_report.R:207-213`), and the rest do not replay `args_call_all`.

| Repo | File:line | Call as written | Thresholds | Exposed? | Reported products: replay / not replay |
|---|---|---|---|---|---|
| forestsearch | quarto/dina/smoke-tests/dina_screening_smoketest.qmd:169 → :95-100 (base list), :112 | `battery("continuous", dcont, "continuous", "gaussian", 0.2, ...)`; `run(list(use_dina = TRUE))` | none given, so `effect.threshold = NULL` and `hr.threshold` defaults to 1.25 | **YES** | Pass/fail `.tc()` checks on `f1$dina_cuts` and `f1$dina_res`, all from the fit; no bootstrap or CV |
| forestsearch | same file :126-127 | `use_dina = TRUE, dina_args = list(max_subgroups = 3L, max_per_covariate = 1L)` | none | **YES** | Pass/fail checks, from the fit; no replay |
| forestsearch | same file :134-135 | `use_dina = TRUE, dina_args = list(scope = "harm", m_diff = harm_m_diff)` (0.2; this frontier key does not set the selected-cut floor) | none | **YES** | Pass/fail checks, from the fit; no replay |
| forestsearch | same file :141 | `use_dina = TRUE, dina_cuts = c("x1 <= 0.3")` | none | no (the DINA block is skipped: main :3180 `is.null(dina_cuts)`) | — |
| forestsearch | quarto/simulations/actg175/continuous/sim_fs_maxeffCons_mr_field_md_template.qmd:690-698 (`subgroup_method` from `FS_MD_METHOD`, :139; the mddina campaign uses "dina") | `outcome_type = "continuous", effect_measure = "MD", effect.threshold = md_threshold` (:257, `md_threshold <- 30`); `use_dina <- FALSE` (:269) | effect.threshold = 30 (the meta of all four committed `*_mddina_combined_1_2000.rds` payloads records `effect_threshold 30`) | no | Feeds summary_continuous_field_mddina.qmd and md_dina_metrics.csv, which fs-glms-interpretable manuscript/_R/glm_oc_tables.R:35 and manuscript/_S4_operating_characteristics_supplement.qmd:249 onward read |
| forestsearch | quarto/dina/smoke-tests/validate_screen3.R:23 | calls the internal `.forestsearch_dina_select(..., outcome_type="continuous", hr.threshold=Fstar, ...)` directly | hr.threshold given | no (bypasses forestsearch() resolution) | — |
| forestsearch | quarto/applications/actg175/analysis_actg175_continuous_{compare_all,intervals,oc}.qmd | `use_dina = FALSE` (:434 / :118 / :119), consistency engine | effect.threshold = 10 (or `fs_md_threshold` = 10) | no (not DINA) | — |
| forestsearch | quarto/simulations/actg175/continuous/mr_coverage_sweep_md_harm.qmd:310-313; scripts_mdf1/rep1_field.R:17 | consistency engine (`methods_run <- "consistency"`, :81) | effect.threshold = md_threshold (30) | no | — |
| fs-glms-interpretable | quarto/actg175/analysis_actg175_continuous_mr.qmd:288-301, analysis_actg175_continuous_mr_pstar099.qmd:288-301 | `use_dina = FALSE` | effect.threshold = 10 | no (not DINA) | — |
| fs-post-selection | manuscript/, manuscript_rev1/ | no continuous-outcome forestsearch call; the ACTG175 sections read the binary payload (`_actg175_section6_2.qmd:19`); the DINA payloads are GBSG survival (`payload_manifest.qmd:92-94`) | — | no | — |
| forestsearch | vignettes/forestsearch.qmd:406, vignettes/articles/survival-post-selection.qmd:170, vignettes/articles/actg175-glm-outcomes.qmd:447 | `subgroup_method = "dina"` on survival / survival / binary OR | thresholds set | no (not continuous) | — |

The committed render of the exposed smoke test (`dina_screening_smoketest.html`) was last changed in d555f01c, the commit that introduced `use_dina`. That commit is an ancestor of 00c72dd1, which added `dina_args$selected_only`.

#### Would it alter a reported number?

No. No table or figure in the scoped manuscripts, vignettes or applications comes from a continuous DINA call with both thresholds unset. The one continuous DINA campaign behind reported numbers (mddina) sets `effect.threshold = 30`, so its fit and its replays use floor 30. The only exposed calls feed a smoke test's pass/fail checks, and that smoke test has no bootstrap or CV.

#### Reproducers

None run. Reading the source settled the verdict.

#### Verifier corrections

- *wrong_line*: Three line numbers in the use_dina block are off. The guard `if (use_dina && is.null(dina_cuts))` is at :3180, not :3179. The selected_only = FALSE branch `fr$cut_expr` is at :3262, not :3260. The quoted `m_diff_sel` / `dina_subgroup()` block runs to the closing parenthesis at :3251, not :3250. None of these changes the verdict.
- *misquote*: The section reports the output of a grep for `hr\.threshold\s*(<-|=)` as finding only :1670, :1811 and argument positions at :2187, :2550, :2766 and :3452. The same grep at master also matches :114 (a sprintf string), :1463 (roxygen), :2281 and :2332 (`hr.threshold == 1.25` comparisons) and :3604 (`consistency_overrides$hr.threshold <- effect_threshold`). The substantive point still holds: :1811 is the only assignment to the local `hr.threshold`. The :3604 match is a list element that carries the resolved `effect_threshold` (0) to the consistency stage, so it supports the statement that the search on the use_dina route screens at 0.
- *other*: The task asks for a ONE-line answer on whether this is the same site as the earlier 'use_dina unoriented floor' review. The section answers in a paragraph. The content is correct, so the corrected section puts a one-line answer first and keeps the supporting detail after it.


### A-003 — Specificity and NPV swapped when no subgroup is found

#### Sites at master

All four sites are at the same line numbers as at base 5b994177. The checked-out code was confirmed with `sed -n` at HEAD aae7b683.

**FS site 1**: R/run_simulation_analysis.R:1009-1033. Condition: `fs_res` is empty and there is no `treat.recommend == 0` member.
```r
1009  if ((is.null(fs_res) || nrow(fs_res) == 0) && !.has_treat_rec) {
1010    out$hr.Hc.hat <- out$hr.itt
1011    out$size.H  <- 0L
1012    out$size.Hc <- nrow(df)
1014    # Classification when nothing found: everyone is predicted Hc.
1015    # Matches fs.estimates.out() from the paper's original code.
1016    if (harm_col %in% names(df)) {
1017      true_H <- df[[harm_col]] == 1L
1018      n_true_H <- sum(true_H)
1020      if (n_true_H > 0) {
1023        out$sens <- 0
1024        out$ppv  <- 0
1025        out$spec <- sum(!true_H) / nrow(df)  # = |Hc|/N
1026        out$npv  <- 1                          # all Hc correctly in Ĥc
1027      }
1030      if (n_true_H == 0) {
1031        out$spec <- 1
1032        out$npv  <- 1
1033      }
```
**FS site 2**: R/run_simulation_analysis.R:1054-1068. `fs_res` has rows but no `M.*` factor is recovered.
```r
1054    if (length(sg_factors) > 0)
1055      df$sg_hat <- create_subgroup_indicator(df, sg_factors)
1056    else {
...
1061      if (harm_col %in% names(df)) {
1062        true_H <- df[[harm_col]] == 1L
1063        n_true_H <- sum(true_H)
1064        if (n_true_H > 0) {
1065          out$sens <- 0;  out$ppv <- 0
1066          out$spec <- sum(!true_H) / nrow(df);  out$npv <- 1
1067        } else {
1068          out$spec <- 1;  out$npv <- 1
```
**GRF site 1**: R/run_simulation_analysis.R:1381-1394. `grf_est` is NULL. `.run_grf_analysis_gen` passes `grf_est = NULL` when the GRF function is missing (:1197-1198), and its tryCatch returns NULL on error (:1257-1262).
```r
1381  if (is.null(grf_est)) {
...
1387    if (harm_col %in% names(df)) {
1388      true_H <- df[[harm_col]] == 1L
1389      n_true_H <- sum(true_H)
1390      if (n_true_H > 0) {
1391        out$sens <- 0;  out$ppv <- 0
1392        out$spec <- sum(!true_H) / nrow(df);  out$npv <- 1
1393      } else {
1394        out$spec <- 1;  out$npv <- 1
```
**GRF site 2**: R/run_simulation_analysis.R:1408-1420. There is no `sg.harm.id` or no `treat.recommend`.
```r
1408  if (!has_sg || !has_tr) {
...
1413    if (harm_col %in% names(df)) {
1414      true_H <- df[[harm_col]] == 1L
1415      n_true_H <- sum(true_H)
1416      if (n_true_H > 0) {
1417        out$sens <- 0;  out$ppv <- 0
1418        out$spec <- sum(!true_H) / nrow(df);  out$npv <- 1
1419      } else {
1420        out$spec <- 1;  out$npv <- 1
```
**Found-case formulas.** FS at R/run_simulation_analysis.R:1119-1129; GRF is identical at :1524-1533.
```r
1119  if (harm_col %in% names(df)) {
1120    true_H <- df[[harm_col]] == 1L
1121    hat_H  <- df$sg_hat == 1L
1123    tp <- sum(true_H & hat_H);   fp <- sum(!true_H & hat_H)
1124    tn <- sum(!true_H & !hat_H); fn <- sum(true_H & !hat_H)
1126    out$sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
1127    out$spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
1128    out$ppv  <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
1129    out$npv  <- if ((tn + fn) > 0) tn / (tn + fn) else NA_real_
```
The output table holds four classification columns against the truth: `ppv, npv, sens, spec` (:945 and :1330). The other per-replicate columns (`any.H`, `size.H`, and so on) are not classification against the truth. The no-find branches set `attr(out, "sg_hat") <- rep(0L, nrow(df))`, which is correct for an empty set.

#### Formula stored in each case

Notation: H is the true subgroup (`harm_col == 1`), Ĥ is the identified set, N = nrow(df). The table covers the alternative, |H| > 0.

| column | subgroup found (:1126-1129 / :1530-1533) | none found: stored (all four sites) | found-case formula applied to Ĥ = ∅ (tp = 0, fp = 0, tn = \|Hc\|, fn = \|H\|) |
|---|---|---|---|
| sens | tp/(tp+fn) | 0 | 0/\|H\| = 0 |
| spec | tn/(tn+fp) | \|Hc\|/N | \|Hc\|/\|Hc\| = 1 |
| ppv | tp/(tp+fp), NA if 0/0 | 0 | 0/0, which the code returns as NA |
| npv | tn/(tn+fn), NA if 0/0 | 1 | \|Hc\|/N |

Under the null (|H| = 0, so |Hc| = N), the branches store spec = 1 and npv = 1 and leave sens and ppv as NA. The found-case formula on an empty set gives the same: spec = N/N = 1, npv = N/N = 1, sens = 0/0 = NA, ppv = 0/0 = NA. The swap cannot be seen under the null.

#### Definitions where the columns are reported

- quarto/simulations/actg175/binary_methods/actg175_binary_m1_harm_methods.qmd:1294-1299. The same table appears in 106 tracked qmd files, and it is the one shown in the rendered binary_methods HTML.
  > `| **Spec** | $\text{sens}(\hat{Q}^c)$ | Fraction of true $Q^c$ correctly in $\hat{Q}^c$ |`
  > `| **NPV** | $\text{ppv}(\hat{Q}^c)$ | Fraction of $\hat{Q}^c$ truly in $Q^c$ |`

  The convention follows at :1301-1306: "When no subgroup is detected ($\hat{Q} = \emptyset$), the simulation contributes $\text{sens}(\hat{Q}) = 0$ and $\text{ppv}(\hat{Q}) = 0$".
- vignettes/articles/paper_simulations.qmd:1833 "**Specificity**: P(not identified | true complement)" and :1835 "**NPV**: P(true complement | not identified)". Also :1404 "`Spec` | Mean specificity" and :1406 "`NPV` | Mean negative predictive value".
- R/simulation_tables.R:239-242, the `build_classification_table` documentation: "`sens(Hc)`: Mean specificity" and "`ppv(Hc)`: Mean negative predictive value". The code comment at :360-364 gives the 0/0 convention: "ppv(Ĥ) = #{Ĥ ∩ H}/#{Ĥ} = 0/0 (convention: 0)".
- quarto/guides/review_classification_calculations.qmd:118 `spec = tn / (tn + fp)    # fraction of complement correctly excluded` and :120 `npv  = tn / (tn + fn)`. Its edge-case table at :489 reads: "| No detection (H1) | 0.00 | 1.00 | NA | $|Q^c|/N$ |" for sens | spec | ppv | npv.
- A conflicting convention is documented in the package for a different function. R/fs_mr_oc_summary.R:56 (`fs_mr_oc_summary`, unconditional rows) reads "with non-detection scored as sens = 0, ppv = 0, spec = 1, npv = 1". The NA fill in `build_classification_table` (R/simulation_tables.R:371-372) and the template's `oc_agg_from_results` (actg175_binary_methods_template.qmd:357-360) use the same values. This convention agrees with the found-case formula for sens, spec and ppv (ppv by the 0/0 convention). For npv it scores 1, which is the value the extractors store but not the value the found-case formula tn/(tn+fn) or the reported definition ppv(Q̂c) gives.

With Ĥ = ∅ we have Ĥc = everyone. Then sens(Ĥc) = #{Ĥc ∩ Hc}/#{Hc} = 1, and ppv(Ĥc) = #{Ĥc ∩ Hc}/#{Ĥc} = |Hc|/N.

#### Verdict per column (alternative, no subgroup found; the same at all four sites)

- **sens: AS DESIGNED.** The stored 0 equals the found-case formula's value and the stated convention.
- **ppv: AS DESIGNED.** The found-case code gives 0/0, which it stores as NA (:1128). The reporting convention, quoted above from the qmd at :1301-1306, simulation_tables.R:364 and R/fs_mr_oc_summary.R:56, scores it 0. The stored 0 follows that convention.
- **spec: DEFECT.** The found-case formula tn/(tn+fp) and the reported definition sens(Ĥc), "fraction of true Qc correctly in Q̂c", both give 1, and so does every documented non-detection convention. The code stores |Hc|/N. The stored value is the NPV (ppv(Ĥc)) of the empty set, which is also the overall accuracy. Rows of the `spec` column therefore mix specificity (found replicates) with NPV (no-find replicates), so the column cannot be relabelled as a single other metric.
- **npv: DEFECT.** The found-case formula tn/(tn+fn) and the reported definition ppv(Ĥc) both give |Hc|/N. The code stores 1, which is the specificity of the empty set. The code comment at :1026, "all Hc correctly in Ĥc", describes specificity. The stored 1 does match the non-detection convention documented at R/fs_mr_oc_summary.R:56 and used by the NA fill at R/simulation_tables.R:372. The task's DEFECT test compares against the column's own found-case formula, so the verdict stands, but the package documents conflicting values for this column.

The downstream fill does not change the stored values. R/simulation_tables.R:370-372 replaces only `NA` values (`r$spec[is.na(r$spec) & r$any.H == 0] <- 1`), and when `harm_col` is present the extractor never leaves NA under an alternative. The same holds for the template's `oc_agg_from_results` at actg175_binary_methods_template.qmd:357-360 (`spec_u = ifelse(is.na(spec) & any.H == 0, 1, spec)`).

#### Reproducer

Script /tmp/audit_revisit_2026-11-11/A-003/repro_a003.R; output in /tmp/audit_revisit_2026-11-11/A-003/repro_a003_out.txt. It uses synthetic data (N = 200, |H| = 50) and calls `.extract_fs_estimates_gen(fs_res = NULL, fs_full = NULL)` and `.extract_grf_estimates_gen(grf_est = NULL)` through `pkgload::load_all`. Command: `timeout 120 Rscript --vanilla repro_a003.R`. Wall clock: 2.22 s. Output, verbatim:
```
|H| = 50  |Hc| = 150  N = 200  |Hc|/N = 0.75 
  analysis any.H sens spec ppv npv
1       FS     0    0 0.75   0   1
2      GRF     0    0 0.75   0   1
found-case formula on empty set: sens = 0  spec = 1  ppv = NA  npv = 0.75 
wall 2.22 s
```
The last line of the output file (`wall 2.22 s`) was written by the shell's timing wrapper, not by the script.

#### Consumers and whether they average over all replicates

- `build_classification_table` (R/simulation_tables.R:374-381): `"sens_sgc" = mean(r$spec, na.rm = TRUE)` and `"ppv_sgc" = mean(r$npv, na.rm = TRUE)`, where `r <- res[res$analysis == a, ]` (:358). That is ALL replicates. Only `avg|H|` uses `r_found`.
- `format_oc_results` (R/oc_analyses.R:688-690): `spec <- mean(res$spec, na.rm = TRUE)` and `npv <- mean(res$npv, na.rm = TRUE)` over `res <- results[results$analysis == a, ]` (:677). That is ALL replicates. The comment at :683 says "averaged across sims with detection", but the code does not filter. Columns Spec and NPV are shown when `metrics` is "all" (the default) or includes "classification" (:787-788).
- `summarize_single_analysis` (R/oc_analyses.R:453-457): `class_cols <- c("any.H", "sens", "spec", "ppv", "npv")`, then `sapply(..., mean, na.rm = TRUE)` over ALL replicates.

#### Exposure (forestsearch, tracked files in quarto/ and vignettes/)

Each file below calls `run_simulation_analysis()` (harm_col defaults to `"flag_harm"`, :403), then summarises spec/npv over ALL replicates of an alternative scenario. The 17 tracked renders show nonzero sens(H), or in the n700_demo grf.eff cell sens(Hc) = 0.89 = 1 − p_H, which confirms `harm_col` was present in the data.

| group (count) | representative call site (file:line) and arguments as written | summary over | render tracked |
|---|---|---|---|
| vignettes/articles/paper_simulations.qmd (1) | :1027 `format_oc_results(results_alt, metrics = c("detection","classification","cde_estimates"), digits = 2, ...)`; :1156 `build_classification_table(scenario_results = scenario_list, ...)` with `alt = list(..., hypothesis = "alt")` (:1146-1151); :1346 `mean(res$spec, na.rm = TRUE)` over `results_alt` (chunk key-metrics, evaluated) | all replicates, H1 | no |
| quarto/simulations/actg175/binary_methods (18: template + 17 instances) | template :681 `format_oc_results(results_alt, metrics = c("detection", "classification", "subgroup_size"), ...)`; :732 `build_classification_table(scenario_alt, analyses = m, subgroup_notation = sg_notation, ...)` with `hypothesis = "alt"` (:729); FIGURE :751-758 `oc_agg_from_results(results_alt)` then `plot_oc_metrics(oc_disp_alt, x = "Spec", y = "NPV", ...)`; e.g. m1_harm_methods_n1000_full.qmd :1359/:1452/:1508 | all replicates, H1 | 16 html (15 with a same-named qmd, plus actg175_binary_m1_harm_methods_n700_demo.html, whose source qmd is not tracked under that name) |
| quarto/simulations/actg175/binary_legacy (72) | e.g. actg175_binary_m1_harm_maxSG_fs1.qmd :781 `format_oc_results(` (classification included); :875/:882 `build_classification_table(scenario_alt, analyses = "FS"/"GRF", ...)` with `hypothesis = "alt"` | all replicates, H1 | 1 html (effMinSG-both_fs3_recommended) |
| quarto/simulations/_legacy (14), quarto/simulations/actg175/_archive_pre_fsparams (3), quarto/archive (8), quarto/applications/gbsg/_archive (1), "quarto/legacy calibration files" (3: actg175_binary_benefit_simulations.qmd, actg175_continuous_simulations.qmd :611/:645/:663, calibration_survival_fourway.qmd :678-690) | same `format_oc_results` / `build_classification_table(... hypothesis = "alt")` pattern; calibration file `compare_metrics(results_alt_cox, results_alt_irr)` → `mean(cox$spec, na.rm = TRUE)`, `mean(cox$npv, na.rm = TRUE)` | all replicates, H1 | no |
| quarto/simulations/reproduce_table1_m1.qmd and quarto/simulations/gbsg/reproduce_table1_m1.qmd (2) | :530 `format_oc_results(results_alt, analyses = c("FSl", "FSlg", "GRF", "GRF60"), metrics = c("detection", "classification"), digits = 2, ...)`; :563 `build_classification_table(scenario_results = list(null = ..., alt = list(..., hypothesis = "alt")))` | all replicates, H1 | no |
| quarto/simulations/actg175/survival (2), quarto/simulations/actg175/continuous/actg175_continuous_simulations.qmd (1) | survival_benefit :2107/:2129/:2189/:2251/:2597; gbsg_poisson :813/:845/:862; continuous :927/:973/:991; all `format_oc_results` / `build_classification_table` on alt results | all replicates, H1 | no |
| quarto/calibration/calibration_survival_fourway.qmd (1) | :824 `compare_metrics(results_alt_cox, results_alt_irr)` → :797/:799/:807/:809 `round(mean(cox$spec, na.rm = TRUE), 3)`, `round(mean(cox$npv, na.rm = TRUE), 3)` | all replicates, H1 | no |

Total: 126 source files, with 17 tracked renders.

Possibly exposed, provenance unknown: quarto/guides/glm_extension_introduction.qmd:592-603 prints a static H1 table, "Survival Fourway — Cox PH (C) vs. Poisson+offset (A)", with Specificity 0.964 / 0.963 and NPV 0.983 / 0.974. Its rows have the same structure as `compare_metrics()` (calibration_survival_fourway.qmd:790-813, all replicates). :612-622 prints "Binary Benefit — FS vs. GRF" with Specificity 0.88 / 0.79 and NPV 0.87 / 0.96. Both tables are hand-entered and were committed in 19c83f45 (2026-04-14), one day after the no-find branch was introduced (95e2f819, 2026-04-13).

Not exposed in forestsearch:
- quarto/archive/actg175_simulations_binary_working.qmd defines and uses its own `extract_binary_estimates()` (:472). It leaves spec/npv as NA when nothing is found (:481-484, :539) and does not call `run_simulation_analysis()`.
- quarto/applications/gbsg/_broken/gbsg_analysis.qmd (:1476) and gbsg_analysis_2sim-approaches.qmd (:1176) apply `format_oc_results` to `results_null` only.
- quarto/grf/grf_scoring_binary_simulation.qmd:627/634 averages over `any.H == 1` only.
- quarto/qc/*.R, the calibration/selection frontier files and the extreme_subgroups files call the extractor without a spec/npv summary.
- The other quarto/guides/*.qmd files and vignettes/articles/treatment_effect_definitions.qmd mention the functions only in text or in `eval=FALSE` chunks.
- quarto/simulations/gbsg_redux/* and the other tracked renders that show Spec/NPV (gbsg/fs_subgroup_id_sweep*.html, gbsg_020/summary_idsweep.html, actg175/continuous/oc_wrapper_verification.html, vignettes/articles/mrct-consistency-simulations.html) are built by other code paths, not by these extractors.
- Every MR-harness summary (quarto/simulations/gbsg_020/*, gbsg/fs_subgroup_id_sweep*.qmd, actg175/continuous/md_harm_mr_simulation_summary.qmd, o1_readout.R) reads payloads with a `detected` column built by the harness's own `.classify()` (e.g. binary_020/sim_fs_mr_field_or_template.qmd:533-543), not by these extractors, and averages over `detected %in% 1L` only.

#### Sibling repositories

- **fs-post-selection** (manuscript/ and manuscript_rev1/): not exposed. No tracked file in these directories references `run_simulation_analysis`, `build_classification_table`, `format_oc_results` or the extractors. None of the generators named in manuscript/payload_manifest.qmd calls `run_simulation_analysis()`. Payloads come from the MR harness, a different code path, and every summary is detected-only. Examples: manuscript_rev1/_sim_estimation_fs_h10.qmd:125 `rH <- res[res$detected %in% 1L, , drop = FALSE]` then :141-142 `.mfin(rH$spec) ... .mfin(rH$npv)`; manuscript_rev1/_sim_synthesis_h10.qmd:59 `r <- o$results[o$results$detected == 1, , drop = FALSE]  # detected replicates only`; manuscript_rev1/_supp_sim_synthesis.qmd:85; manuscript_rev1/fs_identification_figures.qmd:42 `det <- r[r$detected == 1L, , drop = FALSE]` then :54-55. The manuscript/ copies follow the same pattern.
- **fs-glms-interpretable** (quarto/ and manuscript/): not exposed to these extractors. Table `tbl-supp-classification` (manuscript/_S4_operating_characteristics_supplement.qmd:177, `oc_tbl_classification(OC)`, manuscript/_R/oc_figures.R:377-398, rows `sens/spec/ppv/npv` × `uncond/cond`) reads identifiers_metrics.csv. That CSV is built by dev/payloads/R/identifiers_metrics.R:359-362 from MR-harness payloads, which overwrite non-declaring rows: `fill <- c(sens = 0, spec = 1, ppv = 0, npv = 1)` and `ifelse(det, r[[k]], fill[[k]])`. Its spec_uncond therefore scores non-declaration as 1, which matches the found-case formula. Its npv_uncond scores non-declaration as 1, the same value as the extractors' npv branch. That table is not produced by the four A-003 sites.
- **R/fs_mr_oc_summary.R `.fs_oc_identification`** (:280-301) sets non-detected rows to `spec = 1` and `npv = 1` through `fill()`. Its unconditional rows are stored in MR payloads (`$oc$identification`). No tracked qmd/R in the scoped directories displays them (grep for `$identification` / `fs_oc_identification`).

#### Magnitude from stored payloads

A read-only scan (/tmp/audit_revisit_2026-11-11/A-003/scan_rds2.R, 8.6 s) ran `readRDS()` on all 2,463 tracked .rds files under forestsearch quarto/, looking for any data frame with `any.H` or `spec` columns. An independent read-only re-scan of the same 2,463 files, looking for data frames holding `any.H` with `spec` or `npv`, nested up to depth 4, also found none. No file holds `run_simulation_analysis()` per-replicate classification output. The only `any.H` tables are quarto/simulations/gbsg_app_null/results/gbsg_app_null_{full,timing}.rds: a null scenario with no `spec`/`npv` columns. None of the payloads the exposed files read or write is tracked: the binary_methods bundles (`BUNDLE_DIR <- "_bundles"`, template :84), `_data/actg175_binary_simulations_200alts_200nulls.rds`, `_data/actg175_survival_benefit_simulations_50alts_50nulls.rds`, and `latest_rds("calibration_survival_fourway")`. The vignette computes at render time. So per-cell detection rate, reported mean and corrected mean cannot be computed from payloads for any exposed table (see open items).

What the code determines exactly: in each no-find replicate r under an alternative, stored spec − defined spec = −|H_r|/N, and stored npv − defined npv = +|H_r|/N. A cell mean therefore moves by (1/R)·Σ over no-find replicates of |H_r|/N, which is the no-find share times the mean of |H_r|/N over the no-find replicates. In the tracked renders the H1 tables print any(H) from 0.00 to 1.00 with header p_H = 11%.
- One cell can be read off a render without a payload, because every replicate in it is a no-find. In actg175_binary_m1_harm_methods_n700_demo.html, table "grf-eff: Classification Under H1", any(H) is 0.00, sens(H) 0.00, sens(Hc) 0.89, ppv(H) 0.00 and ppv(Hc) 1.00. Each reported value is the stored no-find value: spec = |Hc_r|/N averaged gives 0.89, and npv = 1. Scoring the no-find replicates by the definitions gives sens(Hc) ≈ 1.00 and ppv(Hc) ≈ 0.89, within rounding, since any(H) < 0.005. Which GRF no-find branch (:1392 or :1418) produced them cannot be told from the render.
- Example of a partial-detection cell: actg175_binary_m1e_harm_n1000_full_adjusted.html, dina.eff: any(H) 0.75, sens(Hc) 0.91, ppv(Hc) 0.94. This cell differs from its definition by 0.25 × (mean |H_r|/N over no-find replicates), and the exact amount needs per-replicate |H_r|, which is not stored.
- Cells printed with any(H) 1.00 (actg175_binary_m1_harm_methods.html grf.eff / grf.effMaxSG / grf.effMinSG; actg175_binary_methods_template.html grf.tree) have a no-find share below 0.005, so their shift (at most about 0.0006) cannot be seen at two decimals.

#### Verifier corrections

- *missed_exposure*: A tracked render that the investigator did not count shows an H1 cell where detection is 0. quarto/simulations/actg175/binary_methods/actg175_binary_m1_harm_methods_n700_demo.html (tracked, added in e3fb1852 on 2026-06-11; no tracked qmd has the same basename) contains ten 'Classification Under H1' tables. In the grf.eff table, every replicate is a no-find, so the stored value of each column can be read straight off the table. Reported sens(Hc) = 0.89 is the mean of |Hc_r|/N (p_H = 11%). Reported ppv(Hc) = 1.00 is the stored npv = 1. The investigator lists 15 binary_methods renders plus 1 binary_legacy render; the correct count is 16 plus 1, which is 17.
- *magnitude_error*: The would-alter statement says the rendered H1 tables show any(H) between 0.45 and 0.98 and that 'every H1 sens(Hc)/Spec cell is too low'. The H1 tables in the tracked renders actually show any(H) from 0.00 to 1.00. Where any(H) prints as 1.00, no-find replicates make up less than 0.5%, so the shift is at most about 0.005 x 0.11, about 0.0006, and cannot be seen at two decimals. The investigator's phrasing of the shift is also loose. It is (1/R) times the sum of |H_r|/N over no-find replicates, which equals the no-find share times the mean of |H_r|/N over the no-find replicates only, not over all replicates.
- *false_exposure*: quarto/archive/actg175_simulations_binary_working.qmd does not call run_simulation_analysis() or the A-003 extractors. It defines its own extractor, extract_binary_estimates() at :472, which leaves sens/spec/ppv/npv as NA when nothing is found (classification is computed only when sum(in_H) > 0, at :538-549). Its tables are not produced by the four A-003 sites. It is the only file in the investigator's 127 that does not contain run_simulation_analysis. The corrected count is 126 (quarto/archive: 8).
- *missed_exposure*: The investigator set quarto/guides/*.qmd aside as text-only or eval=FALSE. But quarto/guides/glm_extension_introduction.qmd (tracked, added in 19c83f45 on 2026-04-14) prints static H1 operating-characteristic tables with Specificity and NPV rows: :592-603 'Survival Fourway - Cox PH (C) vs. Poisson+offset (A)' (Specificity 0.964 / 0.963, NPV 0.983 / 0.974, detection 0.833 / 0.880) and :612-622 'Binary Benefit - FS vs. GRF' (Specificity 0.88 / 0.79, NPV 0.87 / 0.96, detection 0.74 / 0.72). The first table's rows match the structure of compare_metrics() in quarto/calibration/calibration_survival_fourway.qmd:790-813, which averages over all replicates of run_simulation_analysis() output. The table was committed one day after the no-find branch was introduced (95e2f819, 2026-04-13). Its provenance and the package version behind it cannot be determined.
- *other*: The investigator cites R/fs_mr_oc_summary.R:56 as a source for the ppv = 0 convention but leaves out that the same sentence scores non-detection as spec = 1 and npv = 1. The package's own build_classification_table() NA fill (R/simulation_tables.R:371-372) and the template's oc_agg_from_results() (:357-360) also use npv = 1 for non-detection. So a convention documented in the package for non-detection gives npv = 1, the value the extractors store. That convention conflicts with the found-case formula tn/(tn+fn) and with the definition used where the NPV column is reported (ppv(Q-hat^c), 'Fraction of Q-hat^c truly in Q^c'; 'P(true complement | not identified)'). Under the task's DEFECT test (the found-case formula on an empty set gives |Hc|/N, not 1), the npv verdict stays DEFECT, and the conflicting convention is now stated in the section. The spec verdict is unaffected: every one of these conventions gives spec = 1.
- *wrong_line*: Minor line errors. The format_oc_results comment 'averaged across sims with detection' is at R/oc_analyses.R:683, not :684. In the vignette, the alt scenario list starts at :1146 and hypothesis = "alt" is at :1151, not :1149. In fs-glms-interpretable dev/payloads/R/identifiers_metrics.R, the fill and ifelse lines are at :359 and :362, not :360-363.


### A-163 — GBSG generator: censoring coefficients are rescaled by the event model's scale

#### Claim (audit, "suspected")
At `R/sim_aft_gbsg.R:638` (line unchanged between base and master), the censoring-model coefficients are converted back to the AFT scale with the event model's scale `sigma` instead of the censoring model's scale `sigma_cens`.

#### Event model: fit, coefficients, draw
The event model is fitted once, and `sigma`, `mu` and `gamma` all come from `fit_aft`:

`R/sim_aft_gbsg.R:375-385`
```r
  fit_aft <- survival::survreg(
    survival::Surv(y, event) ~ z_true,
    data = dfa,
    dist = "weibull"
  )
  names(fit_aft$coefficients) <- c("(Intercept)", colnames(z_true))

  # Extract parameters
  sigma <- fit_aft$scale
  mu <- stats::coef(fit_aft)[1]
  gamma <- stats::coef(fit_aft)[-1]
```
Lines 414-420 then apply the multipliers to `gamma` (`gamma["z3"] <- k_z3 * gamma["z3"]`, `gamma["treat"] <- k_treat * gamma["treat"]`, `gamma["zh"] <- k_inter * gamma["zh"]`). The AFT coefficients are used unchanged, and the hazard-scale coefficients are derived with the event model's own scale:

`R/sim_aft_gbsg.R:438-441`
```r
  # In AFT: log(T) = mu + X*gamma + sigma*epsilon
  # Hazard: h(t|X) = h_0(t) * exp(X * b0) where b0 = -gamma/sigma
  b_true <- gamma
  b0 <- -gamma / sigma
```
`R/sim_aft_gbsg.R:456-458`
```r
  lin_conf <- as.vector(z_true %*% b_true)
  lin1_conf <- as.vector(z_true_1 %*% b_true)
  lin0_conf <- as.vector(z_true_0 %*% b_true)
```
These values are carried to `df_samp$lin1.conf` / `lin0.conf` (`:486-487`) and to `lin.conf.true` (`:592`, `:598`). Event times are drawn as follows.

`R/sim_aft_gbsg.R:606` (super-population): `log_Ts_obs <- mu + sigma * epsilon_obs + df_big$lin.conf.true`

`R/sim_aft_gbsg.R:897-902` (legacy simulator):
```r
  mu <- dgm$model_params$mu
  sigma <- dgm$model_params$sigma
  ...
  epsilon <- log(stats::rexp(n_obs))
  log_T <- mu + sigma * epsilon + df_sim$lin.conf.true
```
`R/simulate_from_dgm.R:310-311` (current simulator; `params$tau` is set to `model_params$sigma` at `R/setup_gbsg_dgm.R:173`):
```r
  epsilon  <- log(rexp(n))    # Extreme value (Gumbel) errors
  logT_sim <- params$mu + params$tau * epsilon + df_sim$lin_pred_obs
```
**Event equation.** With `eps` = log(Exp(1)), so that `eps` has the standard minimum-extreme-value distribution:

  log T = mu + sigma·eps + z′gamma

Here `mu`, `sigma` and `gamma` (after the k-multipliers) all come from `fit_aft`. The hazard-scale coefficients are `b0 = -gamma/sigma`, which is the correct Weibull AFT-to-PH conversion with the model's own scale. The event model is internally consistent.

#### Censoring model: fit, coefficients, draw
`R/sim_aft_gbsg.R:625-657`
```r
  if (cens_type == "weibull") {
    z_cens <- as.matrix(df_big[, covs_main])

    fit_cens <- survival::survreg(
      survival::Surv(y, 1 - event) ~ z_cens,
      data = df_big,
      dist = "weibull"
    )

    sigma_cens <- fit_cens$scale
    mu_cens <- stats::coef(fit_cens)[1]
    gamma_cens <- stats::coef(fit_cens)[-1]
    b0_cens <- -gamma_cens / sigma_cens
    b_cens <- -b0_cens * sigma

    # Censoring linear predictors for potential outcomes
    z_cens_1 <- z_cens
    z_cens_1[, 1] <- 1
    z_cens_0 <- z_cens
    z_cens_0[, 1] <- 0

    df_big$linC1.conf <- as.vector(z_cens_1 %*% b_cens)
    df_big$linC0.conf <- as.vector(z_cens_0 %*% b_cens)
    df_big$linC.conf <- ifelse(df_big$treat == 1,
                                df_big$linC1.conf,
                                df_big$linC0.conf)

    cens_model <- list(
      type = "weibull",
      mu = mu_cens,
      sigma = sigma_cens,
      gamma = gamma_cens
    )
```
`sigma` on line 638 is the event scale. In this file, `sigma` is assigned only at `:383` (`sigma <- fit_aft$scale`) and at `:898` (inside the separate simulator function). `covs_main` is `c("treat", "z1", "z2", "z3", "z4", "z5")` (`:332`), so column 1 of `z_cens` is `treat`.

The returned object stores the fitted censoring parameters, not the scaled ones:

`R/sim_aft_gbsg.R:703-707`
```r
    cens_params = list(
      type = cens_type,
      mu = mu_cens,
      sigma = sigma_cens,
      model = cens_model
```
**Where the censoring predictor is used.** No later step multiplies or divides it by a scale. There are two paths that draw censoring times.

(a) Legacy simulator. `linC1.conf` / `linC0.conf` are copied to `linC.conf` by assigned arm (`R/sim_aft_gbsg.R:865-867`, `:873-875`), then drawn:

`R/sim_aft_gbsg.R:910-916`
```r
  if (cens_type == "weibull" && !is.null(dgm$cens_params$mu)) {
    mu_cens <- dgm$cens_params$mu + muC_adj
    sigma_cens <- dgm$cens_params$sigma

    epsilon_cens <- log(stats::rexp(n_obs))
    log_C <- mu_cens + sigma_cens * epsilon_cens + df_sim$linC.conf
    C_sim <- exp(log_C)
```
(b) `setup_gbsg_dgm()`. The columns are renamed with no transformation (`R/setup_gbsg_dgm.R:157-158`: `df_s <- .rename(df_s, "linC1.conf", "lin_pred_cens_1")`, `df_s <- .rename(df_s, "linC0.conf", "lin_pred_cens_0")`). The scale is mapped at `R/setup_gbsg_dgm.R:175-179` (`mu = dgm$cens_params$mu`, `tau = dgm$cens_params$sigma`). The draw is:

`R/simulate_from_dgm.R:320-328`
```r
    lin_pred_cens <- ifelse(df_sim$treat_sim == 1,
                            df_sim$lin_pred_cens_1,
                            df_sim$lin_pred_cens_0)

    epsilon_cens <- if (cens_type == "weibull") log(rexp(n)) else rnorm(n)

    logC_sim <- params$censoring$mu + cens_adjust +
                params$censoring$tau * epsilon_cens + lin_pred_cens
    C_sim    <- exp(logC_sim)
```
One other consumer reads the stored predictor but draws nothing. The internal `.fs_random_cens_prob()` (`R/betaHhat_truth.R:218`, called from `.fs_exact_setup()` at `:257`) computes P(C ≤ analysis_time) from `dgm$df_super$lin_pred_cens_0/1`. It uses the same composition as path (b) and does not rescale the predictor (`R/betaHhat_truth.R:231-236`: `lp <- c(ds$lin_pred_cens_0, ds$lin_pred_cens_1)`; `z <- (log(analysis_time) - cm$mu - cens_adjust - lp) / cm$tau`). If it is given a `setup_gbsg_dgm()` object, it inherits the scaled predictor.

#### Algebra
The fitted censoring model (`fit_cens`, Weibull `survreg`) is:

  log C = mu_cens + z′gamma_cens + sigma_cens·eps

What the code builds:
- `b0_cens = -gamma_cens / sigma_cens` (`:637`). This is the correct PH-scale censoring coefficient, using the censoring model's own scale.
- `b_cens = -b0_cens · sigma = gamma_cens · sigma / sigma_cens` (`:638`).
- `linC*.conf = z_cens* %*% b_cens = z′gamma_cens · (sigma / sigma_cens)` (`:646-647`).

What is drawn (`:915`, or `simulate_from_dgm.R:326-327` with `cens_adjust = 0` / `muC_adj = 0`):

  log C = mu_cens + sigma_cens·eps + z′gamma_cens·(sigma/sigma_cens)

Covariate effects on the AFT scale:
- Generated: gamma_cens·sigma/sigma_cens
- Fitted: gamma_cens

Covariate effects on the PH scale:
- Generated censoring log-hazard coefficient: −(gamma_cens·sigma/sigma_cens)/sigma_cens = b0_cens·sigma/sigma_cens
- Fitted: b0_cens

The generated covariate effects match the fitted censoring model only when `sigma == sigma_cens`. The intercept `mu_cens` and the scale `sigma_cens` are the fitted values. Only the covariate term is rescaled.

The event model's analogous step in the same file keeps `b_true <- gamma` (`:440`). The equivalent round trip elsewhere in the package converts back with the same model's scale:

`R/generate_aft_dgm_helpers_outcome.R:58-62`
```r
  tau <- fit_aft$scale
  gamma <- coef(fit_aft)[-1]  # Coefficients (excluding intercept)

  # Weibull parameterization
  b0 <- -gamma / tau
```
`R/generate_aft_dgm_helpers_outcome.R:76`: `gamma <- -b0 * tau`

Applied to the censoring model, the same round trip would give `-b0_cens * sigma_cens = gamma_cens`. Line 638 uses `sigma` instead.

#### Is there a stated reason that makes it correct as written?
None was found:
- No comment at or near `:637-638` (the only nearby comment is `# Censoring linear predictors for potential outcomes`, `:640`).
- The roxygen for `cens_type` / `cens_params` (`R/sim_aft_gbsg.R:108-109`, `:128`) says only "Censoring distribution type" and "List with censoring model parameters".
- The `setup_gbsg_dgm` roxygen (`R/setup_gbsg_dgm.R:27-41`) describes the columns only as renamed `lin_pred_cens_*`.
- No downstream step divides the stored predictor by `sigma` or multiplies it by `sigma_cens/sigma` (paths (a) and (b) above add it directly, and `.fs_random_cens_prob()` uses it the same way).

The package's stated contract for `lin_pred_cens_*`, which the `setup_gbsg_dgm()` path feeds, is the fitted censoring coefficients:

`R/generate_aft_dgm_helpers.R:9-20`
```r
#' \code{lin_pred_cens_0} and \code{lin_pred_cens_1} store the
#' \strong{covariate contribution only} — i.e. \eqn{\gamma_c' X}, with the
#' intercept \eqn{\mu_c} excluded.  This matches the convention used for the
#' outcome model (\code{lin_pred_0}, \code{lin_pred_1} = \eqn{\gamma' X},
...
#' \code{simulate_from_dgm()} reconstructs the full log-censoring time as:
#' \deqn{\log C = \mu_c + \delta + \tau_c \epsilon + \gamma_c' X}
...
#' \eqn{\gamma_c' X} = \code{lin_pred_cens_{0|1}}.
```
The general generator meets this contract by using `predict(fit_cens, type = "linear") - mu_cens` (`R/generate_aft_dgm_helpers.R:265-266`, `:272-273`). The GBSG predictors instead equal `z′gamma_cens·sigma/sigma_cens`. They do not match the coefficients stored in the same object at `cens_params$model$gamma = gamma_cens` (`:656`), and they do not match this contract.

The line's provenance is pointer-only and is not evidence: `git log -S'b_cens <- -b0_cens * sigma' -- R` shows it was introduced in `ae8d54a2` (2026-01-20, "gbsg sims"). Copies under tracked `dev/code-working/.../legacy*/R/sim_aft_gbsg.R:601` carry the same line, also without a comment.

#### Scope of the effect, from the code only
- `cens_type` defaults to `"weibull"` (`R/sim_aft_gbsg.R:206`, `:243`, `:253`; `R/setup_gbsg_dgm.R:95`, `:103`; `calibrate_k_inter` `R/sim_aft_gbsg.R:1008`). Every default GBSG DGM therefore takes this branch. With `cens_type = "uniform"` the branch is skipped (`:625`).
- The stored truths do not use the censoring predictor: `hr_causal`, `hr_H_true`, `hr_Hc_true`, `AHR*` and `theta_*` / `loghr_po` are built from `b0` / `lin*.conf` and from stacked potential outcomes with `event = 1L` (`:456-468`, `:498-560`). The deviation affects only the covariate dependence of simulated censoring times, including treatment dependence through column 1 = `treat`.

#### Verdict
**INCONSISTENT.** The event model uses its own scale throughout (`b_true <- gamma`, `b0 <- -gamma / sigma`; drawn with `sigma`). The censoring model's coefficients are converted to the PH scale with `sigma_cens` (`:637`) but back to the AFT scale with the event scale `sigma` (`:638`). The censoring times are then drawn with `mu_cens` and `sigma_cens` plus `z′gamma_cens·sigma/sigma_cens`, not with the fitted `z′gamma_cens`. No stated reason makes this correct as written. The equality holds only if `sigma == sigma_cens`.

#### Generator entry points that reach `R/sim_aft_gbsg.R:638`
| Entry point | Export status | Path to `.create_gbsg_dgm_` (file:line at master) |
|---|---|---|
| `create_gbsg_dgm` | exported (deprecated) | `R/sim_aft_gbsg.R:199` → `:226` |
| `setup_gbsg_dgm` | exported | `R/setup_gbsg_dgm.R:106` |
| `calibrate_k_inter` | exported | `R/sim_aft_gbsg.R:1004` → `:1024` |
| `get_dgm_with_output` | internal | `R/sim_aft_gbsg.R:1076` → `:1093` (`calibrate_k_inter`), `:1108` |
| `validate_k_inter_effect` | exported | `R/sim_aft_gbsg.R:1229` → `:1254` |
| `fs_fdr_report` | exported | `R/fs_fdr_report.R:181` → `:256`, `:262` (`setup_gbsg_dgm`), `:258` (`calibrate_k_inter`) |

Censoring-time draws that consume the scaled predictor:
- `simulate_from_gbsg_dgm` (exported, deprecated): `R/sim_aft_gbsg.R:789` → `.simulate_from_gbsg_dgm_` `:823` → `:915`
- `simulate_from_dgm` (exported), when given a `setup_gbsg_dgm()` object: `R/simulate_from_dgm.R:326-327`

No exposure search of committed analyses was performed, because the task limits the INCONSISTENT case to entry points.

#### Reproducers
None. Reading only; no R was executed.

#### Verifier corrections

- *other*: The task asks to follow every use of linC1.conf / linC0.conf / lin_pred_cens_* through to where the censoring time is drawn. The section lists two paths, (a) the legacy simulator and (b) setup_gbsg_dgm() -> simulate_from_dgm(). It misses one further consumer that reads the stored predictor: the internal .fs_random_cens_prob() in R/betaHhat_truth.R, called from .fs_exact_setup() at :257. This function draws nothing. It computes P(C <= analysis_time) from dgm$df_super$lin_pred_cens_0/1 using the same composition as the draw (mu_c + cens_adjust + tau_c*eps + lp). It does not rescale lp, so the verdict and the statement that no later step multiplies or divides the predictor by a scale both still hold. If it is given a setup_gbsg_dgm() object, it inherits the sigma/sigma_cens-scaled predictor. This consumer is outside the generator-entry-point list the task asks for, so it is recorded only as a completeness note in the trace.


## 4. OPEN ITEMS

- **A-009**: none.
- **A-001**: CANNOT DETERMINE, because no reproducer was run (reading settled the verdict): whether any pass/fail check in the continuous battery of quarto/dina/smoke-tests/dina_screening_smoketest.qmd would come out differently at master with floor 0, for example the check that dina_cuts is non-empty. The committed render predates the selected-cut route (last changed in d555f01c, an ancestor of 00c72dd1, which added selected_only), and no table or figure reports these checks.
- **A-003**: CANNOT DETERMINE the per-cell magnitude (detection rate, reported mean, mean with no-find replicates scored spec = 1 / npv = |Hc|/N) for any exposed table from stored payloads, because no tracked .rds in the scoped directories holds run_simulation_analysis() per-replicate output. Two independent read-only scans of all 2,463 tracked quarto .rds files found none; the payloads the exposed files read ('_bundles', _data/actg175_*_simulations_*.rds, calibration_survival_fourway rds) are untracked; and the vignette computes its results at render time. The one cell recoverable without a payload, from the n700_demo render with any(H) = 0.00, is reported in the section.
- **A-003**: CANNOT DETERMINE which package version produced the 17 tracked renders (quarto/simulations/actg175/binary_methods/*.html, 16 of them including actg175_binary_m1_harm_methods_n700_demo.html, and binary_legacy/actg175_binary_m1_harm_effMinSG-both_fs3_recommended.html), because the HTML records no forestsearch version and the qmd files use library(forestsearch). The renders are dated 2026-06-11 and 2026-06-24, after the no-find branch was introduced (95e2f819, 2026-04-13), and the branch is present at e3fb1852, the commit that added the 2026-06-11 renders. The n700_demo grf.eff cell (sens(Hc) 0.89, ppv(Hc) 1.00 at any(H) 0.00) matches the no-find branch's output.
- **A-003**: CANNOT DETERMINE the source qmd of actg175_binary_m1_harm_methods_n700_demo.html, because no tracked qmd shares its basename and the render's title is shared with actg175_binary_m1_harm_methods.html.
- **A-003**: CANNOT DETERMINE whether the static H1 Specificity/NPV tables in quarto/guides/glm_extension_introduction.qmd:592-622 were produced with the no-find branch, because they are hand-entered numbers committed on 2026-04-14 (19c83f45), one day after the branch was introduced (95e2f819), with no recorded source run or package version.
- **A-003**: CANNOT DETERMINE whether a reported number changes in vignettes/articles/paper_simulations.qmd or the 109 other exposed source files without a tracked render, because no rendered output or per-replicate payload for them is tracked.
- **A-163**: CANNOT DETERMINE the numerical size of the factor sigma/sigma_cens at the default settings (and so how far the generated censoring covariate effects are from the fitted ones), because that needs fitting fit_aft and fit_cens, and no R was run for this reading-only step.
- **A-163**: CANNOT DETERMINE whether any reported table or figure changes, because the task limits the INCONSISTENT case to listing generator entry points; committed analyses were not searched for exposure.
- **A-163**: CANNOT DETERMINE the author's intent behind line 638, because no comment, roxygen, test or NEWS entry at master says which censoring model the GBSG DGM is meant to reproduce; the verdict rests on the code's own stored cens_params$model$gamma and the lin_pred_cens_* contract at R/generate_aft_dgm_helpers.R:9-20.

Deviation from the standing report convention: no `VERIFY_*.txt` / verification script is committed beside this report, because the task's Step 6 gate requires exactly two added files. The independent second-reader pass described in §1 stands in for it.
