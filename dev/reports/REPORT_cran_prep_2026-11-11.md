# REPORT — CRAN preparation: confirmed codebase fixes, documentation, and a check of the remaining audit findings, 2026-11-11

Task: `dev/tasks/TASK_cran_prep_2026-11-11.md`, run unattended on 2026-11-11. This task was run as Step 0 of `TASK_cran_prep_seed_legacy_2026-11-11.md` (part 2), because branch `cran-prep` did not exist. Nothing was pushed or merged.

## 1. Provenance and commits

- **Host and tools.** pop-os (128 cores), R 4.6.2, devtools 2.5.2, testthat 3.3.2, pkgload 1.5.3, roxygen2 8.1.1, future 1.80.0, doFuture 1.4.0.
- **Step 0 gates (all PASS).**
  - On `master`, `git status --porcelain` was empty.
  - `master` was at `afe96278` (= `origin/master`).
  - `DESCRIPTION` was 0.4.1.9006.
  - Branch `cran-prep` did not exist.
- **Start SHA `e92d3c96`:** the TASK commit on `master` (`dev: add CRAN preparation task`). `cran-prep` was created from it. Its R code is identical to `afe96278`.
- **Sibling repositories.** fs-post-selection (`5d165acc`, `main`) and fs-glms-interpretable (`dc94a17e`, `main`) were clean before and after. They were not read or written.
- **How the work was done.**
  - **Part 1.** One agent per fix implemented it in its own scratch `git worktree`, detached at the start SHA.
    - An independent adversarial reviewer checked each patch against the rule, the start-SHA source and the tests.
    - F1–F6 were approved in the first review. F7 went through two repair rounds; two remaining review points were closed here (§3.7).
    - The patches were applied here, one fix per commit. Each fix's own test file was re-run on the branch before its commit.
    - F2 was narrowed before it was applied (§3.2).
  - **Part 3.** One agent per audit item read a `git archive` export of the start SHA. An independent verifier then re-derived every column of each row.

**Commits on `cran-prep` (oldest first):**

| SHA | Commit | Item |
|---|---|---|
| `e92d3c96` | dev: add CRAN preparation task (on `master`; branch point) | Step 0 |
| `159225fc` | fix: GLM bootstrap uses the parent fit's adverse_outcome and ps_adjust_method (A-009) | F1 |
| `7607f113` | fix: continuous DINA floor is the resolved threshold when none is supplied (A-001) | F2 |
| `c2fbc0d6` | fix: survival GRF covariate matrix keeps exact numeric values and codes factors (A-059) | F3 |
| `e7241e09` | fix: NPV is 0 when the identified set is everyone (empty complement) | F4 |
| `67ce5b41` | fix: generate_gbsg_bootstrap_general() takes survival::gbsg (A-004) | F5 |
| `56b2429e` | fix: functions that set a future plan restore the caller's plan (INV-2) | F6 |
| `c47e26bf` | docs: stated argument defaults match the code (INV-3, C-001); version 0.4.1.9007 | F7, `DESCRIPTION`, `NEWS.md` heading |
| `c40d5309` | docs(paper_simulations): key-metrics loop variable no longer shadows the analysis column | Part 2 step 1 |
| (this commit) | dev: CRAN preparation report | Part 5 (this report, Part 2 outputs, Part 3 detail) |

## 2. Status for Larry

| Item | Fixed / BLOCKED / FAILED / not as described | What changes for a user | Outputs changed on a small example |
|---|---|---|---|
| F1 GLM bootstrap estimator (A-009) | **Fixed.** The propensity half was not blocked: every slice the bootstrap estimates on carries `sw`, `ps_hat` and `ips_covar`. | GLM bootstrap estimates are on the parent fit's scale. Fits with resolved `adverse_outcome = FALSE` (continuous at defaults) were sign-flipped. Propensity-adjusted fits are now weighted like the parent. `forestsearch()`'s own propensity-adjusted rebuild keeps `adverse_outcome`. 7 quarto drivers no longer re-orient. | Continuous, defaults: `H_estimates$H0` −1.0941 → **1.0941** (the parent's estimate); every H/Hc estimate, bound, bias-adjusted value, `SG_CIs` and `FSsg_tab` MD column flips sign, SDs unchanged. With `adverse_outcome = TRUE` and no PS, 0 of 74 outputs change. With iptw: H0 1.0941 → 1.0792. |
| F2 continuous DINA floor (A-001) | **Fixed**, narrowed to the both-unset case (§3.2). | A continuous DINA analysis with neither threshold supplied uses floor 0, the recorded and replayed value, instead of 1.25. | `subgroup_method = "dina"`: floor 1.25 → 0; qualifying candidates 67 → 1,878; `sg.harm` {x1 ≤ −0.534} (88) → {x1 ≤ 2.42} (300), equal to the replay's. `use_dina`: `dina_cuts` none → `x1 >= 0.042`. A supplied threshold gives the same floor as before. |
| F3 survival GRF matrix (A-059) | **Fixed** | Survival GRF with a factor or character confounder: no `NA` labels and no 7-digit rounding. The GRF matrix and the evaluated cut agree. | Label factor: error "Covariate matrix X contains missing values" → subgroup found (`forestsearch(subgroup_method = "grf")`; `use_grf = TRUE` keeps its GRF cuts). Numeric-level factor: `bm <= 42.603982999999999` → `42.603982585942454`; evaluated members 138 → 137 (= fitted 137). All-numeric: `identical()` matrix and return. |
| F4 NPV when everyone is identified | **Fixed** | NPV of a replicate whose identified set is everyone is 0, not `NA` or 1. | `run_simulation_analysis()` FS and GRF extractors: `npv` NA → **0**. `fs_mr_oc_summary()` unconditional NPV: one such replicate 1 → 0; 4-replicate mix 0.9125 → 0.6625. A non-empty complement is unchanged. |
| F5 `generate_gbsg_bootstrap_general()` (A-004) | **Fixed** (confirmed failing first) | The function runs. | `n = 200, seed = 42`: error "'data' must be a data frame" → a 200 × 11 data frame; the default call → 686 × 11. |
| F6 caller's `future` plan (INV-2, plan half) | **Fixed** | The caller's plan is restored after `subgroup.search()` (so `forestsearch()` and its in-session callers), `forestsearch_bootstrap_dofuture()`, `mrct_region_sims()`, `fpr_calibration()` and `fs_fdr_report()`. | The plan after each call is the caller's `multisession`, 2 workers (it was `sequential`). Fits are identical. |
| F7 documentation (INV-3, C-001) | **Fixed** (documentation only) | Help pages state the formals' defaults: 19 sites where the text named a wrong value; 20 sites where the formal is `NULL` and the text gave only what it resolves to; `est.scale` ∈ {"hr", "1/hr"}. | None. 27 `R/` files (roxygen lines only), 30 `man/` pages. |
| `paper_simulations` key-metrics chunk | **Fixed** (`c40d5309`) | §7.1 reports FS and GRF separately; they were pooled. | FS Power 0.750 → 0.886, GRF 0.750 → 0.614; FS/GRF PPV 0.691 → 0.804 / 0.578; Type I error 0.0341 → 0.0455 / 0.0227; Sens 0.692 → 0.774 / 0.610; Spec 0.991 → 0.990 / 0.992. |
| `paper_simulations` re-render | **Done**: 88 s (cap 1 h); library gate PASS; gates 5, 6, 12a, 12b PASS. | — | 1,544 numbers. 10 changed, all in §7.1 (the chunk). 21 run-dependent values changed (timings, `gc()`), counted only. Nothing changed from F1–F6. |

## 3. Part 1 — fixes

Each fix's test file was re-run on the branch after the patch was applied, before its commit, with `NOT_CRAN=true`, `pkgload::load_all()` and `testthat::test_file()`. Gate: 0 failures. All seven passed.

| Fix | Test file | Blocks | Expectations | Failed | Wall (measured) | At the start SHA |
|---|---|---|---|---|---|---|
| F1 | `tests/testthat/test-bootstrap-glm-parent-scale.R` (new) | 3 | 24 | 0 | 14.5 s | 10 of 24 fail |
| F2 | `tests/testthat/test-dina-gaussian-floor.R` (new) | 5 | 31 | 0 | 5.3 s | 4 failures + 1 error in the both-unset blocks; supplied-threshold blocks pass |
| F3 | `tests/testthat/test-grf-survival-matrix.R` (new) | 4 | 40 | 0 | 6.8 s | block 1 (identity) passes; blocks 2–4 fail/error |
| F4 | `tests/testthat/test-npv-empty-complement.R` (new) | 5 | 48 | 0 | 2.8 s | 7 of 48 fail |
| F5 | `tests/testthat/test-gbsg-bootstrap-general.R` (new) | 3 | 22 | 0 | 2.5 s | 1 failure + 3 errors |
| F6 | `tests/testthat/test-plan-restored.R` (new) | 8 | 20 | 0 | 18.5 s | the 8 plan checks fail |
| F7 | — (gate: roxygen lines only) | — | — | — | — | — |

The "At the start SHA" column was measured by the reviewers. Each ran the new file on a copy of the start-SHA tree. Every file that runs `forestsearch()` carries `skip_on_cran()` and a sequential plan. The existing test files that touch each changed function were re-run in the fix's worktree, and all passed. They are listed per fix in the scratch notes (`fixnotes/F*.md`). The full suite (§6) covers them again.

### 3.1 F1 — GLM bootstrap uses the parent fit's estimator settings (A-009)

**Site at the start SHA.**
- `R/bootstrap_dofuture_main.R:300-310`: `estimator_fn_boot <- make_effect_estimator(outcome_type, treat.name, outcome.name, event.name = NULL, offset.name, effect_measure, adjust_covariates)`.
  - It passes neither `adverse_outcome` nor `ps_adjust_method`, so they fall back to the formals: `adverse_outcome = TRUE` (`R/glm_effect_estimators.R:107`) and `ps_adjust_method = "none"` (`:110`, `:113`).
  - The closure feeds `H_obs`/`Hc_obs` (`:328`), every `bootstrap_results()` constituent (`R/bootstrap_analysis_dofuture.R:485, :514, :695, :708-709, :747, :760-761`) and `FSsg_tab` (`:707-712`).
- `R/forestsearch_main.R:2703-2712`: the parent's propensity-adjusted rebuild omits `adverse_outcome` too.
- The parent's resolved values are stored in `args_call_all`: `adverse_outcome` at `:1797-1799` (captured at `:1881`), `ps_adjust_method` at `:2641`.

**Diff** (`159225fc`):
- `R/bootstrap_dofuture_main.R`: the call now passes `adverse_outcome = boot_adverse_outcome` and `ps_adjust_method = args_forestsearch_call$ps_adjust_method`.
  - `boot_adverse_outcome` is `args_forestsearch_call$adverse_outcome`.
  - It falls back to `forestsearch()`'s own rule (`outcome_type %in% c("binary", "count")`) only when that element is `NULL`. That happens only for an object not built by `forestsearch()`; without the guard, the closure would fail and the error would be swallowed as `NA`.
- `R/forestsearch_main.R:2711`: `adverse_outcome = adverse_outcome` added to the propensity-adjusted rebuild.
- **Propensity half: not BLOCKED.** The reviewer logged all 14 slices the estimator received in a 2-replicate iptw run. `fs.est$df.est`, the resampled `df_boot`, `df_PredBoot` and the replicate's `df.est` all carry `sw`, `ps_hat` and `ips_covar`.
- **Consumers edited (source only, no re-run).** The `.orient` / `.orient_c` re-orientation of `fb$H_estimates` / `fb$Hc_estimates` was removed from the three named drivers:
  - `quarto/simulations/actg175/continuous/sim_fs_maxeffCons_fb_mr_md40_knoise0_n500_batch_1_{100,20,50}.qmd`.
- The same edit was made in the four drivers that carry the identical recorder behind `if (run_fb)` (skipped there because `nb_boots` is `NULL`). Left as they were, they would flip the sign twice if FB were turned on:
  - `sim_fs_maxeffCons_mr_md40_knoise0_n500_batch_1_100.qmd`
  - `sim_fs_maxeffCons_mr_md40_knoise0_n500_batch_1_1000.qmd`
  - `sim_fs_maxeffCons_mr_md40_knoise0_n700_batch_1_1000.qmd`
  - `sim_fs_maxeffCons_mr_md120_knoise0_n500_batch_1_1000.qmd`
- All seven purl and parse.

**Test** (three blocks).
- `adverse_outcome` FALSE and TRUE: the bootstrap's `H_estimates$H0` equals the parent's estimate for its subgroup, sign included.
- iptw parent: the bootstrap estimator carries `ps_adjust_method`, and H0 equals the parent's weighted estimate.
- The parent's propensity-adjusted rebuild keeps `adverse_outcome`.

The "unchanged with TRUE" half was checked by the implementer and the reviewer with a leaf-by-leaf comparison against the start SHA: 74 leaves, 0 differ.

**Changed outputs** (synthetic continuous, N = 300, `forestsearch(effect.threshold = 0.5, fs.splits = 20, n.min = 30, maxk = 2, sg_focus = "hr")`, then `forestsearch_bootstrap_dofuture(nb_boots = 5)`, sequential; parent subgroup {age ≤ 47} & {biomarker_hi}):

| Output | Start SHA | Branch |
|---|---|---|
| `H_estimates$H0` (default path: continuous, `adverse_outcome` unset → FALSE) | −1.0941122 | 1.0941122 (= parent `max_sg_est`) |
| `H_estimates` H0 bounds; H1; H2 bounds | (−1.6651, −0.5232); −0.6307; (−6.4402, 5.9711) | (0.5232, 1.6651); 0.6307; (−5.9711, 6.4402); `sdH0`/`sdH1`/`sdH2` unchanged |
| `Hc_estimates` H0 (bounds); H1; H2 (bounds) | −0.0958 (−0.3380, 0.1464); −0.0548; −0.0165 (−2.2124, 2.1794) | 0.0958 (−0.1464, 0.3380); 0.0548; 0.0165 (−2.1794, 2.2124) |
| `results$H_biasadj_1` (5 replicates) | −0.2911, −0.4024, −0.8897, −0.7137, −0.8565 | exact negation; likewise `H_biasadj_2`, `Hc_biasadj_*` |
| `SG_CIs` `H_raw` / `H_bc` | −1.0941 (−1.6651, −0.5232) / −0.2345 (−6.4402, 5.9711) | 1.0941 (0.5232, 1.6651) / 0.2345 (−5.9711, 6.4402) |
| `FSsg_tab` "MD (95% CI)" | −1.0941 (−1.6651, −0.5231); −0.0958 (−0.3380, 0.1464) | 1.0941 (0.5231, 1.6651); 0.0958 (−0.1464, 0.3380); `Diff`, `Mean(C)`, `Mean(T)` unchanged (raw) |
| All 74 bootstrap output leaves, `adverse_outcome = TRUE`, no PS | start values | identical |
| `H_estimates$H0` / `sdH0` (TRUE, `is.RCT = FALSE`, `ps_method = "logistic"`, iptw) | 1.0941122 (unweighted) / 0.2913 | 1.0791986 (= the parent's weighted estimate) / 0.2893 |
| `forestsearch()` `sg.harm` (FALSE, `ps_method = "logistic"`, iptw) | none found (the rebuild used TRUE) | {age ≤ 47} & {biomarker_hi}, `max_sg_est` 1.0792 |

### 3.2 F2 — continuous DINA floor comes from the resolved threshold (A-001)

**Site at the start SHA.**
- `subgroup_method = "dina"`: `R/forestsearch_main.R:2766` passes `hr.threshold` to `.forestsearch_dina_select()`, where `R/forestsearch_helpers.R:1525-1526` sets `m_diff <- if (identical(da$fit$family, "gaussian")) hr.threshold else log(hr.threshold)`.
- `use_dina = TRUE`: `R/forestsearch_main.R:3237-3238` sets `m_diff_sel` the same way.
- The recorded and replayed value comes from `:2332-2334`, where MD with both thresholds unset gives `effect_threshold <- 0.0`. Line `:2532` sets `effect.threshold`, and `:2535-2537` syncs it into `args_call_all`. The bootstrap and CV replay that list.
- Measured with both thresholds unset (MD): fit floor 1.25; record and replay 0.

**Diff** (`7607f113`):
- `.forestsearch_dina_select()` gains `effect.threshold = hr.threshold` (the default keeps direct callers unchanged), and its gaussian floor reads it.
- Both `forestsearch()` sites use `effect.threshold` when `user_set_threshold` is FALSE, and `hr.threshold` when a threshold was supplied (`R/forestsearch_main.R:2784-2785` and `:3250-3252` at the branch head).
- Non-gaussian families keep `log(hr.threshold)`.

**Deviation from the agent's patch, decided here.** The agent's version used the resolved threshold for every gaussian fit. In one case that changes the floor even though a threshold was supplied:
- The case: binary RD (or count IRD) with `dina_args = list(family = "gaussian")` and a ratio-scale `effect.threshold` such as 1.5. The package warns and remaps such a value to 0.05.
- The task says "with either threshold supplied … the floor is unchanged", and its test requires "the floor equals the supplied value". The `user_set_threshold` guard keeps that.
- In that case the fit floor (1.5) and the recorded/replayed value (0.05) still differ, as at the start SHA. This is OPEN ITEM 1.
- A test block pins it: floor 1.5 on both routes, recorded 0.05.

**Test.**
- Both thresholds unset, on both routes: the fit's floor equals `threshold_config$screening`, `args_call_all$effect.threshold` and (dina route) `admission$effect_floor`, all 0. A bootstrap-style replay of `args_call_all` gets the same floor.
- `effect.threshold = 0.4` or `hr.threshold = 0.4`: floor 0.4 at both sites.
- The binary-RD remap case above.
- A direct helper call without the new formal keeps `hr.threshold`.

**Changed outputs** (continuous MD, N = 300, `sg_focus = "maxSG"`, `n.min = 40`, `maxk = 1`, sequential; scratch `f2work/example.R`):

| Output | Start SHA | Branch |
|---|---|---|
| `subgroup_method = "dina"`, both unset: fit floor `m_diff` | 1.25 | 0 (= record = replay) |
| same: `out_sg$n_candidates_qualifying` | 67 | 1,878 |
| same: `sg.harm` (size) | {x1 ≤ −0.534} (88) | {x1 ≤ 2.42} (300), equal to the replay's selection |
| `use_dina = TRUE`, both unset: `m_diff_sel` | 1.25 | 0 |
| same: `dina_cuts`; `sg.harm` | none; {x2 ≤ 8} | `x1 >= 0.042` (= the replay's); {x2 ≤ 8} |
| binary RD, gaussian family forced, both unset: floor; `sg.harm` | 1.25; none | 0.05 (the recorded value); {x1 ≤ 2.75} |
| either threshold supplied (MD 0.4; binary RD 1.5) | supplied value | unchanged |

### 3.3 F3 — survival GRF covariate matrix (A-059)

**Site at the start SHA.**
- `R/grf_main.R:225-226`: `X <- apply(data[, confounders.name, drop = FALSE], 2, as.numeric)`.
  - With a factor or character confounder, `as.matrix.data.frame()` produces a character matrix. Numeric columns become 7-significant-digit strings, and labels such as "EU"/"US" become `NA`.
- The selected cut is evaluated on coded raw columns: `R/grf_subgroup_labels.R:225` `x <- .grf_code_column(df[[v]])`, called from `R/grf_main.R:308` and `R/forestsearch_helpers.R:1987`.
- The GLM path codes first: `R/grf_subg_harm_glm.R:868` `.build_grf_X()`, `:891-894`.

**Diff** (`c2fbc0d6`): `X <- apply(.build_grf_X(data, confounders.name), 2, as.numeric)`. It reuses the GLM path's coding unchanged, so column names and order are kept for every downstream use of `X`.

**Test.**
1. All-numeric confounders: `X` is `identical()` to the start-SHA construction, dimnames included.
2. Numeric-level factor present: the numeric columns equal the raw values exactly.
3. Fitted and evaluated membership agree.
4. A label factor no longer errors.

**Changed outputs** (`.make_survival_data(N = 600, HR_harm = 2.5, seed = 42)` + `bm`; confounders age, stage, sex, noise, bm; `n.min = 60`, `dmin.grf = 0`, `maxdepth = 2`):

| Output | Start SHA | Branch |
|---|---|---|
| all-numeric confounders (seeds 42 and 7, tree and frontier) | full return object | `identical()` |
| `sg.harm.id`, numeric-level factor present, tree | `bm <= 42.603982999999999 & noise > -0.65838929069999996` | `bm <= 42.603982585942454 & noise > -0.65838929069704344` |
| members evaluated by `.grf_evaluate_subgroup()` (fitted 137) | 138 | 137 |
| `sg.harm.id`, frontier | `bm > 70.41076910000001` | `bm > 70.410769260555512` |
| label factor `region` (EU/US/Asia), tree | error "Covariate matrix X contains missing values" | `bm <= 42.66773421489227 & noise > -0.65838929069704344`, Nsg 138, fitted = evaluated |
| `forestsearch(subgroup_method = "grf")` with the label factor | error | {noise > −0.844} & {noise ≤ 1.343}, 424 members |
| `forestsearch(use_grf = TRUE)` with the label factor: GRF warning; `grf_cuts`; `sg.harm` (size) | "GRF analysis failed: …"; none; {noise ≤ 0} & !{noise ≤ −0.7} (162) | no warning; `noise <= -0.84`, `stage <= 2`, `noise <= 1.34`; !{noise ≤ −0.84} & {noise ≤ 0} (183) |

### 3.4 F4 — identified set is everyone: NPV

**Site at the start SHA.**
- `R/run_simulation_analysis.R:1136` (FS extractor `.extract_fs_estimates_gen()`) and `:1543` (GRF extractor) both have `out$npv <- if ((tn + fn) > 0) tn / (tn + fn) else NA_real_`.
- `R/fs_mr_oc_summary.R:311` (`.fs_oc_identification()`) has `npv_u <- fill("npv", ifelse(det, 1, hc_share))`, so a detected replicate with `NA` npv is filled with 1.

**Diff** (`e7241e09`):
- Both extractors: `else 0`.
- `.fs_oc_identification()`: a detected replicate with no stored npv is filled with 0 when `n_harm` equals `meta$n_sample` (identified set = everyone), else with 1 as before.
- Comments updated. Payloads are not regenerated.

**Test.**
- Both extractors, everyone identified, under the alternative and under the null (via `df.est` and via an `M.*` label): `npv` = 0.
- `fs_mr_oc_summary()` gives 0 for such a replicate.
- A replicate with a non-empty complement, and a detected replicate with `NA` npv and a smaller set, are unchanged.

**Changed outputs:**

| Output | Start SHA | Branch |
|---|---|---|
| FS extractor `npv`, everyone identified (N = 200, 50 truly in H; alternative and null) | `NA` | 0 |
| GRF extractor `npv`, everyone identified | `NA` | 0 |
| `fs_mr_oc_summary()$identification` unconditional npv, one detected replicate with `n_harm` = `n_sample` = 200 and npv `NA` | 1 | 0 |
| same, 4 replicates (everyone `NA`; stored 0.9; detected `NA` with `n_harm` 60; non-detected \|Hc\|/N 0.75) | 0.9125 | 0.6625 |

### 3.5 F5 — `generate_gbsg_bootstrap_general()` (A-004)

**Site at the start SHA.**
- `R/synthetic_data_perturbation.R:365-366`: `gbsg <- NULL` then `data("gbsg", package = "survival", envir = environment())`.
- Confirmed by calling it on the start tree (R 4.6.2, survival 3.8.12): `generate_gbsg_bootstrap_general(n = 50, seed = 1)` warns "data set 'gbsg' not found". It then stops at `:65-67` with "'data' must be a data frame". The default call fails the same way.
- `data(package = "survival")` lists the object as "gbsg (cancer)", i.e. it ships inside the `cancer` data file.

**Diff** (`67ce5b41`): `gbsg <- survival::gbsg` (one line; survival is in Imports).

**Test.** A small call returns without error, the result is reproducible for a fixed seed, and the default call returns 686 rows.

**Changed outputs.**
- `generate_gbsg_bootstrap_general(n = 200, seed = 42)`: error → a 200 × 11 data frame (pid … status). It is `identical()` to the start-SHA `generate_bootstrap_synthetic(data = survival::gbsg, …)`.
- The default call: error → 686 × 11.

### 3.6 F6 — the caller's `future` plan is restored (INV-2, plan half)

**Sites at the start SHA, every plan-setting call in `R/`, and what was done:**

| file:line (start SHA) | Function (exported?) | Action |
|---|---|---|
| `R/subgroup_search.R:174` (`plan(multisession, workers = parallel_workers)`), `:191` (reset to sequential) | `subgroup.search` (yes); reached by `forestsearch()` (`forestsearch_main.R:3496`) and, in the session, `compare_selection_rules`, `run_simulation_analysis`, `run_null_calibration`, `mrct_trial_analysis`, `mrct_case_study` | **Restored**: `old_plan <- future::plan()` and `on.exit(future::plan(old_plan), add = TRUE)` before `:174`. The reset at `:191` is kept. |
| `R/bootstrap_dofuture_main.R:352-356` (on.exit restored `.Last.future.plan`, which nothing creates, else `"sequential"`) | `forestsearch_bootstrap_dofuture` (yes) | **Restored**: the plan is captured first, and the on.exit calls `future::plan(old_plan)`. |
| `R/fpr_calibration.R:272` (`plan("multisession")`), on.exit `:274` (`plan("sequential")`) | `fpr_calibration` (yes, when `n_workers > 1`) | **Restored**: captured, and the on.exit restores it. |
| `R/mrct_simulation.R:357-384` (set by `plan_type`), `:704` (reset) | `mrct_region_sims` (yes) | **Restored**: 2 lines before `:357`. The reset at `:704` is kept. |
| `R/fs_fdr_report.R:398, :399, :411` | `.fs_fdr_run_null` (internal; sole caller `fs_fdr_report`) | **Restored in the exported caller**: 2 lines before `runs <- lapply(...)`. One restore per call, so the caller's workers are not restarted between nulls. |
| `R/fpr_calibration.R:296`, `R/fs_fdr_report.R:406` | inside `%dofuture%` bodies | Left alone: future puts the plan back after each future (checked). |
| `R/subgroup_consistency_main.R:883-888`; `R/forestsearch_cross_validation.R:332-337`, `:840-844`; `R/run_subgroup_sims.R:817-821`, `:829-834` | `subgroup.consistency`, `forestsearch_Kfold`, `forestsearch_tenfold`, `run_subgroup_sims` | Already restore. |
| `R/subgroup_consistency_helpers.R:57, 60, 66, 69` | `setup_parallel_SGcons` (yes) | **Left alone, listed**: its documented purpose is to set the plan ("@return None. Sets up parallel backend as side effect."). |
| `R/forestsearch_helpers.R:806, 834, 836` | `reset_workers` (yes) | **Left alone, listed**: documented purpose ("Tears down the current future plan … spins up a fresh worker pool."). |

**Diff** (`56b2429e`): the 5 restores above. The fixed seed in `subgroup.consistency()` was not touched.

**Test.** The caller's plan is `multisession`, 2 workers. In each block, the plan class and worker count are the same after the call as before it. The blocks:
- `subgroup.search()`: normal exit, and error exit.
- `forestsearch()`: its fit is also identical to the fit under a sequential caller plan.
- `forestsearch_bootstrap_dofuture(nb_boots = 2)`.
- `mrct_region_sims(n_sims = 2)`: normal and error exits.
- `fpr_calibration(n_workers = 3)` and `fs_fdr_report(n_workers = 3)`, each with `%dofuture%` stubbed.
- `setup_parallel_SGcons()` still sets the plan.

**Changed outputs.**
- The plan after each of `subgroup.search()`, `forestsearch()`, `forestsearch_bootstrap_dofuture()`, `mrct_region_sims()`, `fpr_calibration()` and `fs_fdr_report()`: `sequential`, 1 worker → the caller's `multisession`, 2 workers.
- No result changes.

### 3.7 F7 — documentation matches the code (INV-3, C-001)

**What changed** (`c47e26bf`; roxygen comment lines only, then `devtools::document()` once):

**Core sites, 19.** Each text named a value that is not the formal's default.
- C-001: `forestsearch(est.scale)` "("hr" or "rmst")" → "("hr" or "1/hr")". The code accepts only these (`bootstrap_calculations_helpers.R:278`).
- The audit's INV-3 list:
  - `forestsearch(use_twostage)`: FALSE → TRUE.
  - `mr_inference_args$ci_method`: `"field"` only under `subgroup_method = "consistency"`; DINA and GRF default to `"ij"`.
  - `field_complement`: FALSE → TRUE.
  - `get_Cox_sg(cox_initial)`: a starting value, `log(1)`.
  - `calibrate_glm_interaction(k_inter_range)`: c(0, 10) → c(−10, 10).
  - `compute_detection_probability_glm(c2)`: 1.0 on both scales.
  - `plot_subgroup_effects(hrz_crit)`: 1 → 0.
  - `dina_subgroup(n_min.frac)`: stray "Default `60L`" line deleted.
  - `fpr_approx(thresholds)`: `NULL`, only `c1` evaluated.
  - `gg_forest(widths)`: `NULL`, giving `rep(1.1, …)`.
  - `interpret_estimation_table(scenario)`: `NULL`.
  - `bootstrap_results(show_three)`: no default.
  - `subgroup.consistency()` `parallel_args$batch_size`: rule names `"maxeffCons"`.
  - `subgroup_glm()` cross-reference: binary default `"OR"`.
- From the agents' sweep:
  - `forestsearch_Kfold(Kfolds)`: the no-subgroup case.
  - `plot_subgroup_results_forestplot(outcome_type)`: `NULL`, taken from the fit.
  - `forestsearch()` `@return` `mr_inference`: "the default `ci_method = "field"`" qualified as above.
- Added here, for consistency with the `@return` correction: `print.forestsearch()`'s "Post-selection inference" section, "(`ci_method = "field"`, the default)" → "…, the default under `subgroup_method = "consistency"`)".

**`NULL`-formal sites, 20.** The formal is `NULL` and the text stated only what `NULL` resolves to; each now names `NULL` and keeps the verified resolution. Under F7's literal reading these are defaults that differ from the formal.
- `load_leff_calibration(output_dir)`, `mrct_region_sims(conf_force)`, `plot_effect_distribution(panel_labels)`
- `gg_forest(tick_labels)`, `generate_bootstrap_synthetic(n, cat_flip_prob)`, `make_region_model(scale_ref)`
- `grf.subg.harm.glm(effect_measure)`, `fs_oc_predict(c1, c2, pconsistency)` (inherited by `fs_oc_grid` and `fs_oc_invert`), `fs_oc_invert(c1, c2)`
- `plot_pareto_frontier(effect_neighborhood)`, `fs_fdr_report(n)`, `forestsearch(adverse_outcome, ps_method)`
- `fs_identification_structure(confounders)`, `plot_fs_identification_involvement(roles)`, `summarize_simulation_results(analyses)`
- `SGplot_estimates(show_summary)`, added here after the second review found it missing.

**Review.**
- Round 1 found the `plot_subgroup_results_forestplot(outcome_type)` site and the `@return` inconsistency; both were fixed.
- Round 2 found `SGplot_estimates(show_summary)` and the `print.forestsearch` section; both were fixed here before the commit.
- Each corrected text was checked against `formals()` and the resolving code. The agent's table (formal file:line, old/new) is in the scratch notes `fixnotes/F7.md`.
- `NEWS.md` heading: "# forestsearch 0.4.1.9002 (development)" → "# forestsearch 0.4.1.9007 (development)". `DESCRIPTION` 0.4.1.9006 → 0.4.1.9007, once, in this last fix commit.

**Gate: PASS.**
- `git diff -U0 56b2429e c47e26bf -- R/`: 111 changed lines (51 removed, 60 added) in 27 files. All 111 start with `#'`; 0 do not.
- `NAMESPACE` unchanged; 30 `man/*.Rd` pages regenerated.
- Every line is in Appendix A.

### 3.8 `NEWS.md`

One entry per fix F1–F6 was added at the head of "### Behaviour changes" in the fix's own commit. Each says what changes and for whom.

## 4. Part 2 — `paper_simulations`

**Step 1, the chunk** (`c40d5309`).
- In `vignettes/articles/paper_simulations.qmd:1340-1355`, `results_alt` and `results_null` are data.tables. Inside `[`, the loop variable `analysis` resolved to the column, so `results_alt$analysis == analysis` was TRUE on every row, and both lines printed the pooled FS + GRF values.
- The loop variable is renamed `this_analysis` in both loops: 6 lines; nothing else in the article changed.

**Step 2, re-render** (the route of `REPORT_classification_outside_package_2026-11-11.md` §4.2):
- **Install.** `git archive HEAD` at `c40d5309` (after Part 1 and step 1) went to `/tmp/fs_src_2026-11-11b/`. It was installed with `R CMD INSTALL --no-test-load --library=/tmp/fs_lib_2026-11-11b` (16.3 s). No user or site library was touched.
- **Library gate: PASS** (`logs/library_gate.txt`). With `R_LIBS=/tmp/fs_lib_2026-11-11b`, `packageVersion("forestsearch")` = 0.4.1.9007 = `DESCRIPTION`, and `find.package()` is the temporary library. The same holds in a `plan(multisession, workers = 2)` worker, which has its own PID.
- **Scratch copy.** `/tmp/paper_sims_rerun_2026-11-11b/pkg`. The article md5 is `c0d94a6d…` there and in the repository.
- **Render command.** `cd /tmp/paper_sims_rerun_2026-11-11b/pkg; R_LIBS=/tmp/fs_lib_2026-11-11b timeout 3600 Rscript -e 'pkgdown::init_site(); pkgdown::build_article("articles/paper_simulations", lazy = FALSE, quiet = FALSE)'`.
- **Wall clock: 88 s measured** (21:19:10 to 21:20:38; cap 1 h). Exit 0; no error or warning in the log.
- **Versions** are as in the 2026-11-11 render: R 4.6.2, pkgdown 2.2.1, quarto (R) 1.5.1, Quarto CLI 1.10.18, pandoc 3.1.3, 120 workers. Settings and seeds come from the same source; only the chunk differs.
- The suite and `check()` were running at the same time, so the wall clock includes contention.
- **Outputs** (`dev/reports/paper_sims_rerun_2026-11-11b/`):
  - `paper_simulations.html` and its two figures (`paper_simulations_files/figure-html/`, added with `git add -f` because `.gitignore` ignores `*_files/`).
  - `results/`: `paper_simulations.txt`, `paper_simulations_output.txt`, and `keyvalues.txt`, which holds the precedent's gates 5, 6, 12a and 12b, all **PASS**.
  - `strip.sh`, a copy of the precedent's with `SRC` changed.
  - `logs/`: install, library gate, provenance, render log, timing.
  - `accounting/`: the precedent's two scripts with only the side labels changed, `accounting_grouped.csv` (every number), `accounting_tables.md` (one table per block), `accounting_unpaired_lines.csv` (none) and `group_lines.txt`.

**Step 3, accounting** against `dev/reports/paper_sims_rerun_2026-11-11/paper_simulations.html`:
- **Method.** The precedent's: lines aligned with numbers masked, and table cells matched by column label.
- **Coverage.** 849 lines on each side, all paired (1 with changed text); 1,544 numbers.

| Block | Metric | 2026-11-11 | 2026-11-11b | Accounted for by |
|---|---|---|---|---|
| 7.1 Key Metrics | FS Power | 0.750 | 0.886 | the chunk (= the empirical FS rate in §7, 0.886) |
| 7.1 | GRF Power | 0.750 | 0.614 | the chunk (= §7, 0.614) |
| 7.1 | FS Sens | 0.692 | 0.774 | the chunk |
| 7.1 | GRF Sens | 0.692 | 0.610 | the chunk |
| 7.1 | FS Spec | 0.991 | 0.990 | the chunk |
| 7.1 | GRF Spec | 0.991 | 0.992 | the chunk |
| 7.1 | FS PPV | 0.691 | 0.804 | the chunk (= Table 6 FS PPV 0.80) |
| 7.1 | GRF PPV | 0.691 | 0.578 | the chunk (= Table 6 GRF PPV 0.58) |
| 7.1 | FS Type I Error | 0.0341 | 0.0455 | the chunk (= §7, 0.045455) |
| 7.1 | GRF Type I Error | 0.0341 | 0.0227 | the chunk (= §7, 0.022727) |

- **10 numbers changed, all in §7.1.** Each old value is the mean of the two new ones (e.g. (0.886 + 0.614)/2 = 0.750), as expected when the pooling is removed.
- **No other computed number changed.** That includes Tables 6, 8 and 10, the Type I table and the empirical rates. F1–F6 account for no change on this page:
  - No GLM outcome, DINA or bootstrap runs in the article (F1, F2).
  - Survival GRF does run, on the confounders at `:883` (z1–z5, factors coded 0/1, plus `size` and `grade3`), and its rows came out unchanged (F3).
  - No NPV value changed (F4).
  - F5 and F6 do not change a number, and F7 is documentation.
- **Run-dependent: 21 changed, counted not listed.** Timing table §11: 11. `gc()` §17: 6. "forestsearch overall" seconds/minutes in §5.3 and §5.4: 2 + 2. The render date and session lines are unchanged: same day, same R.
- **Groups** (as defined in the precedent):
  - (i) 166 numbers, 2 changed (the two §7.1 Spec values).
  - (ii) 151 numbers, 6 changed (§7.1 Power, PPV, Type I).
  - (iii) 1,227 numbers, 23 changed: the two §7.1 Sens values and 21 run-dependent.

**Step 4.** The pkgdown site was not rebuilt.

## 5. Part 3 — remaining audit findings (read-only; nothing changed)

- **Method.** For each item, one agent read the source in a `git archive` export of the start SHA `e92d3c96`, whose `R/` is identical to `master` `afe96278`; "Site at master" refers to these lines. An independent verifier then re-derived every column; the rows below are the verified rows, condensed.
- **Full detail.** Each item's complete row, with quoted evidence, reproducer and the verifier's corrections, is in `dev/reports/REPORT_cran_prep_2026-11-11_part3_detail.md`.
- **Read-only.** No file was changed for Part 3.
- **Summary.**
  - Every item **holds** (YES). A-062 carries a qualification; A-024 is a method mismatch whose smallest fix is documentation; A-156 is exact in mechanism but negligible in effect.
  - **None changes a `forestsearch()` default-path result.**
  - Items whose fix can change the default-argument output of their own exported function (some only for some data):
    - A-103: `forestsearch_Kfold()` / `forestsearch_tenfold()` "Exact".
    - A-107: `fpr_calibration()` without an outcome type.
    - A-125: `mrct_region_sims()` with `confounders.name` omitted.
    - A-147: bootstrap summary bins, for values on a break.
    - A-150: `generate_bootstrap_with_noise()`.
    - A-151: `test_hte_from_forestsearch()`.
    - A-156: `fs_fdr_report()`, negligibly.
- **After F5.** A-162 now reaches `generate_gbsg_bootstrap_general()` on its default path: `grade` is passed only in `ordinal_vars`, so it is never perturbed. At master the wrapper failed first (A-004).

| ID | Site at master (file:line) | Holds? (YES / NO / PARTLY) | Reachable from an exported function with ordinary inputs? (how; default path or not) | Smallest fix, described | Kind (behaviour / method / documentation) | Would a default-path result change? |
|---|---|---|---|---|---|---|
| A-024 | `R/forestsearch_helpers.R:1525-1535` (DINA fit built without offset.name; floor m_diff = log(hr.threshold)); `R/forestsearch_main.R:3197-3206, 3236-3250` (use_dina fit and floor); `R/dina.R:538-542, 843-847` (Poisson fits of the raw count; no offset formal at `R/dina.R:157-167`) | YES. DINA's Poisson fit models the raw count with no exposure term and no way to pass one, so its tau-hat and log(hr.threshold) floor are on the log mean-count ratio, not the offset IRR forestsearch requires (equal only if follow-up does not depend on arm given covariates). A subgroup whose IRR clears the threshold but whose tau-hat does not is never scored or selected. | Yes: forestsearch(outcome_type = "count") with subgroup_method = "dina" or use_dina = TRUE, where offset.name is necessarily set and ignored by DINA; also dina(), dina_fit(), dina_subgroup() and dina_subgroup_refit() with family = "poisson", none of which accept an exposure. Not the default path (subgroup_method defaults to "consistency", use_dina to FALSE). | Documentation only: at dina_subgroup.R:140 and forestsearch_main.R:445 (and the two Rd files) change "log(IRR) for poisson" to "log mean-count ratio" (no exposure offset; equals log(IRR) only when follow-up does not depend on arm given covariates). Optional: same wording at dina_subgroup_refit.R:43; an add-only warning for the "poisson" family; the larger method fix, an exposure argument in dina()/dina_fit() (dina.R:541, :844). | method (estimand mismatch); the smallest fix is documentation, the optional warning is add-only | No: subgroup_method = "consistency" with use_dina = FALSE never fits DINA, and the doc fix changes nothing at run time. Carrying the offset into DINA would change results only for count outcomes with DINA, materially when mean follow-up differs by arm given covariates. |
| A-042 | `R/get_FSdata_helpers.R:333-338` (selected/omitted built from every coef(fit) row, intercept included); `R/get_fsdata.R:354-357` (guard); `R/get_fsdata.R:362-379` (pool emptied) | YES. For gaussian, binomial and poisson, coef(fit) has an almost always non-zero "(Intercept)" row, so `selected` is never empty and the guard `length(lassokeep) > 0` is always TRUE; when LASSO keeps no covariate, conf.cont_medians and conf.categorical are set to NULL. The Cox path (no intercept row) keeps the whole pool instead. | Yes, not the default path: forestsearch(outcome_type = binary/continuous/count, use_lasso = TRUE); run_simulation_analysis() with a glm_dgm and a GLM fs_params$outcome_type (use_lasso is TRUE by default there); get_FSdata()/lasso_selection() directly. Trigger: lambda.min zeroes every covariate; get_FSdata then stops (sg.harm = NULL, counted as no subgroup) or the search runs silently on forced/default/J-quantile/median-force/GRF/DINA cuts only. | In lasso_selection(), compute selected/omitted over the confounder rows only (drop "(Intercept)" before the != 0 / == 0 split, or index coefs by confounders.name); the existing empty-selection rule at get_fsdata.R:354 then applies to all outcome types. Optional, cosmetic: an outcome-dependent "Cox-LASSO" label at forestsearch_main.R:3344 and the roxygen title at get_FSdata_helpers.R:214-216. | behaviour | No: the default paths of forestsearch() and run_simulation_analysis() are survival. For GLM + LASSO: an all-zero LASSO keeps the whole pool instead of emptying it; "(Intercept)" leaves lasso_selection()$selected and FSdata$lassokeep, and the LASSO counts drop by one; cuts are unchanged when LASSO keeps a covariate. |
| A-047 | `R/get_fsdata.R:499` (GRF merge); `R/get_fsdata.R:520` (DINA merge); `R/get_fsdata.R:451` (replace_med_grf median removal, same unanchored grepl) | YES. The "already represented" test is an unanchored grepl(x, cut), so a categorical name or cut string that is a substring of another GRF/DINA cut discards it (bm vs "bmi <= 25.3"; "x <= 2" vs "x <= 25.3"; at :520 also a kept GRF cut). At :451 a median cut is removed because of a GRF cut on a variable whose name contains it. | Yes: forestsearch() with use_grf = TRUE, use_dina = TRUE or non-empty dina_cuts (not its default path); get_FSdata(); and run_simulation_analysis(), whose default path runs the GRF block (run_fs_grf = TRUE, FSlg arm) but whose default v1..v7 names do not overlap, so only user-supplied overlapping names are hit. Needs an overlapping name or cut string; :451 also needs replace_med_grf = TRUE and cut_type = "median". | At :499 and :520 parse the variable name before "<=" / ">=" and drop a cut only if that name is exactly in conf.categorical or the whole cut is exactly in the current confs, e.g. `trimws(sub('\\s*[<>]=?.*$', '', cut)) %in% conf.categorical \| cut %in% confs`. At :451 compare with %in% against the parsed GRF variable names. fixed = TRUE alone does not fix the substring case. | behaviour | No, with default inputs: forestsearch() defaults never run the blocks; run_simulation_analysis() defaults run the GRF block, but v1..v7 and cut_type "default" give no overlap. Changes only with GRF/DINA screening on and an overlapping name or cut string; the dropped cut then rejoins the pool. |
| A-055 | `R/glm_effect_estimators.R:514-548` (.estimate_rd tier 2: se_rd :536, converged = TRUE :541), reached because tier 1 `:483-512` cannot start at a boundary; correct raw-proportion SE in tier 3 `:552-565`, never reached | YES. With an arm at 0% or 100% events, tier 1 cannot start and the logistic tier silently "converges" (no warning) with an SE that drops the control-arm variance: near 0 or in the hundreds, while the estimate is p1-p0 and converged = TRUE. The formula is also wrong for interior unadjusted slices, which almost never reach tier 2. | Yes, not the default path: make_effect_estimator(outcome_type = "binary", effect_measure = "RD") and forestsearch(outcome_type = "binary", effect_measure = "RD") (defaults: survival; binary resolves to "OR"), with no adjust_covariates (the default); also subgroup_glm, compute_frontier_cis, plot_sg_glm_outcomes and the CV tables on an RD fit. The search's 10-event gate drops 0-event arms but not 100% arms; bootstrap/reporting fits have no gate. | Replace line 536 with `se_rd <- if (n_adj == 0L) sqrt(p1*(1-p1)/n1 + p0*(1-p0)/n0) else abs(avg_deriv) * se_treat`, keeping converged and method_used, so only the SE moves. Do not route these slices to tier 3: its converged = FALSE makes compute_frontier_cis and subgroup_glm drop the results. Both arms degenerate gives SE 0 (not repaired); the adjusted tier-2 SE is a separate issue. | method | No: the default path (survival, or binary resolved to "OR") never calls .estimate_rd. With RD only the SE changes, on slices with a degenerate arm (candidate Wald limits and redundancy keys, seH_obs/seHc_obs, downstream RD CIs); the point estimate and search threshold are unchanged. |
| A-062 | `R/grf_main.R:291-303` vs `:316-332` (no-selection return lacks $candidates; success return attaches .mr_candidates at :331); same in `R/grf_subg_harm_glm.R:525-546` vs `:557-573`; consumer gate `R/forestsearch_helpers.R:1910-1913` | YES, with one qualification. Under frontier + "effect", dmin.grf acts only as a DR-scale existence gate: if no candidate reaches it, no $candidates are returned and no subgroup is reported; otherwise the re-selection ranks the family without dmin.grf. On fallback paths (nothing scorable, caught error) the native DR winner, filtered by dmin.grf, stands. | Yes: forestsearch(subgroup_method = "grf") at its GRF defaults (frontier, "effect", dmin.grf = 0), plus bootstrap and CV replicates. Not forestsearch()'s overall default path (subgroup_method = "consistency"), but the GRF identifier's default configuration. With dmin.grf = 0 the gate fires when every DR harm effect is strictly negative. Direct grf.subg.harm.survival()/grf.subg.harm.glm() calls default to "tree" and do not reach it. | Attach the family in both no-selection returns: `nullres$candidates <- .mr_candidates` in grf_main.R:295-302 and `candidates = .mr_candidates` at grf_subg_harm_glm.R:533-544, leaving the admission effect floor as the only rule under "effect". Not dmin.grf = -Inf (can error via a zero-row frontier selection). Doc-only option: describe the gate in @param dmin.grf / grf_select_statistic and fix "Minimum events for GRF" (forestsearch_main.R:555). | behaviour (which admission rule is intended is a method decision; the doc-only alternative makes it documentation) | No for forestsearch()'s default. Under subgroup_method = "grf" defaults, only when every DR candidate's harm effect is negative and some Cox/GLM effect clears the admission floor: a subgroup is then returned where none is today. The fallback paths and the doc-only fix change nothing. |
| A-064 | `R/grf_subg_harm_glm.R:443-452` (adverse_outcome flip has no count branch; comment :451); `R/grf_subg_harm_glm.R:409-423` (count Y_grf); `R/grf_subgroup_labels.R:274, 307, 361-362` (frontier effect and existence gate) | YES. For count the GRF outcome (log(Y+0.5) or Y) is never negated under adverse_outcome = TRUE, so the subgroup labelled "harm" is where treatment lowers the count; adverse_outcome TRUE and FALSE give identical dr_scores and membership. Independently reproduced. | Yes, not the forestsearch() default path: grf.subg.harm.glm(outcome_type = "count", adverse_outcome = TRUE); forestsearch(outcome_type = "count") (adverse_outcome resolves to TRUE) with subgroup_method = "grf" or use_grf = TRUE. Tree mode or the "dr" statistic report the benefit region; frontier + "effect" picks correctly but the DR existence gate is inverted (harm everywhere gives no subgroup); the use_grf cut set comes from the mis-oriented tree. | Add `else if (outcome_type == "count") Y_grf <- -Y_grf` after the log or identity transform in the adverse_outcome block, and correct the comments at :437-441, :451 and :515-517. Optional clarity edit: an explicit adverse_outcome = TRUE in the roxygen example (:243-256). | behaviour | No for forestsearch() defaults and for direct calls at adverse_outcome = FALSE. Only for count with adverse_outcome = TRUE and GRF: tree/"dr" flips to the harm region; frontier + "effect" changes when the DR gate's answer differs (e.g. harm everywhere: none now, a subgroup after); GRF metadata and use_grf cuts can change. |
| A-103 | `R/forestsearch_cross_validation.R:1810-1818` (CV_sgs depth-1 branch; exact_match :1814 copied to one_match/cov1_match :1817-1818); consumers `:1378-1383, :1395-1407` | YES. With a one-factor original subgroup, :1814 scores a fold exact whenever either fold factor equals it, so a fold that chose {A} & {B} counts as an exact match of {A}, while the package's canonical-equality measure scores it not exact. The depth-2 branch matches equality unless maxk >= 3. | Yes, with ordinary inputs: CV_sgs (exported) via forestsearch_KfoldOut, on the default path of forestsearch_Kfold(fs.est) and forestsearch_tenfold(fs.est, sims). No argument triggers it, only a data outcome: a one-factor full-data sg.harm and a fold refit that selects a two-factor sg.harm containing it. Not on forestsearch()'s own path. | In the sg_depth == 1 branch, compute containment once for one_match/cov1_match, and set `exact_match <- ifelse((sg1 == sg1a & is.na(sg2)) \| (sg2 == sg1a & is.na(sg1)), 1, 0)` with NA -> 0 (verified to match canonical equality). Optional: say "Exact" means the same factor set in the roxygen and the label at cv_summary_tables.R:418. A maxk = 3 analogue is outside the smallest fix. | behaviour | Only when the full-data subgroup has one factor and some fold selects a two-factor subgroup containing it: Kfold's find_summary["Exact"] and tenfold's per-sim rows fall by that share (tenfold's median falls weakly), as do dependent tables and the details line. Other metrics and forestsearch() output are unchanged. |
| A-107 | `R/fpr_calibration.R:236` (otype fallback "binary"); `R/fpr_calibration.R:239-247` (switch, silent d_eff_binary arm at :245, run only when d_eff_approx is NULL); consumers `:249-256, :330-335` | YES. With neither outcome_type_approx nor fs_params$outcome_type, otype falls to "binary" though forestsearch() runs survival replicates, so with d_eff_approx = NULL a survival calibration gets the binary d_eff and prop_cens is ignored. An unrecognised type also falls silently to the binary formula and is stored in call_args$otype. | Yes, on fpr_calibration()'s own default-argument path (exported): a survival analysis with fs_params lacking outcome_type (the roxygen allows it) and outcome_type_approx and d_eff_approx at NULL; fs_params = list() also hits it. Off forestsearch()'s main path; the only internal caller, fpr_calibration_transfer()'s nested call, supplies d_eff_approx (metadata only). The bad d_eff flows into fpr_calibration_transfer(). | At :236 replace the literal "binary" with "survival" (or match.arg against c("survival","binary","continuous","count"), survival first). Secondary: make the switch's default arm at :245 a stop(). Optional roxygen sentence at :89-92. A documentation-only fix would not do. | behaviour | Only when outcome_type_approx = NULL and no fs_params$outcome_type (always survival replicates): with d_eff_approx NULL, d_eff 15 -> 60 and P1 0.2222 -> 0.1565, with L_eff_implied, print lines and transfer outputs following; call_args$otype changes either way. fpr_corrected and the simulated FPR do not change. |
| A-125 | `R/mrct_simulation.R:314-323` (NULL branch building the default confounders.name; setdiff :319-322); `R/mrct_simulation.R:539` (modifyList(fs_args, fs_base_args) overrides fs_args$confounders.name); mirrored `R/mrct_case_study.R:285-292` | YES. The default set contains the region indicator and all three spline helper columns. The region column is constant 0 in the training rows (cut dropped); z_bm_k is a z_bm transform that can add a new threshold (z_bm <= 7); z_bm_treat and z_bm_k_treat carry the seed-trial arm, are noise columns, and enter the candidate pool. | Yes, on the formal default path: mrct_region_sims(dgm, n_sims, ...) with confounders.name omitted, for DGMs from create_dgm_for_mrct() or inject_mrct_structure() (both always build a spline); mrct_case_study() when fs_args lacks confounders.name; three \dontrun roxygen examples. Less ordinary: a confounders.name passed via mrct_region_sims(fs_args = ...) is silently replaced. The vignette, tests and campaign pass confounders.name and avoid it. | In the NULL branch (:314-323), remove the exact names paste0(spline_info$var, c("_treat", "_k", "_k_treat")) when spline_info is non-NULL (not a suffix pattern, which would also hit z_prior_treat); removing region_var is optional (inert). Make the identical change in mrct_case_study() (:288-291) and update the roxygen at :51-52. | behaviour | Only when confounders.name is left NULL (or passed via fs_args, where the default overrides it), and only in some trials: 11 of 25 candidate cuts on the reproducer's data are removed, incl. the otherwise absent z_bm <= 7. Measured 0 differences in 12 + 24 + 15 replicates. The printed "Confounders:" line and the case study's call_args also change. |
| A-147 | `R/summarize_bootstrap_subgroups.R:417-422` (Pcons cut); `R/summarize_bootstrap_subgroups.R:441-449` (N_sg cut); contrast `R/summarize_cv_results.R:759-765` | YES. Neither cut() call sets right, so the bins are right-closed ((0.85,0.95], (0.95,1]; [0,50], (50,100], ...) while the labels ("<50", ">=0.95", ">=300") describe left-closed bins; a value exactly on a break is counted one bin below its label. | Yes, with ordinary inputs, on the default path of exported summarize_bootstrap_subgroups() and summarize_bootstrap_results() (the pkgdown articles' workflow), which calls it whenever results has a Pcons column. Shows only for values on a break: Pcons = 0.95 is common at the default 2-digit rounding and 0.90 threshold; N_sg 50, 100, 150, 200, 300 occur (also in DINA mode). Not on forestsearch()'s path (audit marks it off path). | Add right = FALSE to both cut() calls and keep include.lowest = TRUE: Pcons bins [0,0.5) ... [0.85,0.95) and [0.95,1] (1.0 still counted), N_sg bins [0,50) ... [300,Inf]. This matches the CV sibling at summarize_cv_results.R:764-765. | behaviour | Only when a Pcons or N_sg value falls exactly on a break: at the default threshold only Pcons = 0.95 (moves to ">=0.95"); N_sg breaks move up one bin. Only consistency_dist, size_dist and the gt table built from consistency_dist change; no estimates, CIs or main summary table. |
| A-150 | `R/synthetic_data_perturbation.R:514-518` (continuous auto-detect takes all numerics at :515); `R/synthetic_data_perturbation.R:520-529` (setdiff at :528); doc `:393-397` | YES. With continuous_vars and cat_vars NULL, every numeric column except id_var becomes continuous, and the setdiff then removes numeric low-cardinality columns from cat_vars, so numeric-coded categoricals (gbsg meno, hormon, status, grade) get only Gaussian noise, clamp and round, never the documented categorical change. | Yes, on the default path of exported generate_bootstrap_with_noise(), e.g. generate_bootstrap_with_noise(survival::gbsg). The first roxygen example is the same call but uses data(gbsg, package = "survival"), which loads nothing under survival 3.8.12, so it needs survival attached. Not on forestsearch()'s analysis path; nothing else in the package calls it. | At :515 keep only numeric columns with more than 10 unique values (the detect_variable_types() rule at :317-324) and update @param continuous_vars (:393-394); :522-528 stay. Equally small alternative: auto-detect cat_vars first and setdiff it out. Optional: setdiff a user-supplied cat_vars out of the auto-detected continuous_vars. | behaviour | Yes: with defaults, numeric columns with 10 or fewer values would change (gbsg, seed 123: 0.036-0.055 at noise_level 0.1 vs 0 at master); the RNG stream shifts, so later continuous columns change (age unchanged). No forestsearch analysis result changes. |
| A-162 | `R/synthetic_data_perturbation.R:145` (loop over cat_vars only); `R/synthetic_data_perturbation.R:152` (ordinal branch nested inside it); `R/synthetic_data_perturbation.R:217-228` (note ignores ordinal_vars) | YES. Apart from the formal and the existence check, ordinal_vars is used only in the branch at :152 inside the cat_vars loop, so a variable listed only in ordinal_vars comes back as the plain bootstrap resample. Both roxygen examples and generate_gbsg_bootstrap_general() use that pattern (measured: Example 1's grade identical to the resample). | Yes, with ordinary input but not the default path: generate_bootstrap_synthetic() run as roxygen Example 1 (ordinal_vars = "grade", not in cat_vars; ordinal_vars defaults to NULL). Example 1 is in \donttest, so CRAN checks exercise it silently. generate_gbsg_bootstrap_general() would hit it on its own default path but at master errors first (A-004), so it reaches the defect only once A-004 is fixed. Not on forestsearch()'s path. | Loop at :145 over union(cat_vars, ordinal_vars) and add ordinal_vars to all_specified_vars at :217. Recommended companion, not part of the minimal change: fix A-148 at :168 by index, possible_new[sample.int(length(possible_new), 1)]. Other route: documentation only, stating that ordinal_vars must be a subset of cat_vars and correcting the two examples and the wrapper. | behaviour | Only when ordinal_vars names a variable not in cat_vars (as in the examples and the wrapper): it is then perturbed to adjacent levels (measured change rate 0.077), and other columns keep their RNG draws. Default ordinal_vars = NULL is unchanged; the wrapper changes only after A-004 is fixed. No forestsearch result changes. |
| A-151 | `R/test_hte_crump.R:397-416` (covs from confounders.evaluated :399; NULL column -> next :410-411; silent return(NULL) :416); label source `R/get_fsdata.R:527-533` | YES. On the consistency path confounders.evaluated holds cut labels such as "age <= 55" for every confounder with 4 or more distinct values; looked up as df columns they are silently skipped, so only bare-name categoricals are tested, and when every confounder is cut it returns NULL without a warning. The :401-404 fallback ("q1", "q2", ...) would not help. | Yes, the default way to call exported test_hte_from_forestsearch(fit, df): any successful forestsearch() fit with default settings and at least one confounder (numeric, factor or character) with 4 or more distinct values. Off forestsearch()'s pipeline (nothing in R/ or tests/ calls it); not reached for "dina"/"grf" fits or early-return consistency fits. | Map each label to its variable: `unique(intersect(unlist(lapply(covs, function(l) all.vars(parse(text = l)))), names(df)))`, with tryCatch around parse (keeps the LASSO/GRF screening); warn on unresolved labels and when none remain instead of returning NULL silently. Second choice: args_call_all$confounders.name (drops the LASSO selection). Optional: correct the "GRF/LASSO selected" label. | behaviour | Yes, for this function's own output on any consistency fit with a continuous confounder: covariates_used, n_covariates and the Crump test include them (reproducer: chi-sq 0.113 on 1 df -> 1.309 on 3 df), and continuous-only fits get a test or warning instead of NULL. forestsearch() results do not change. |
| A-156 | `R/fs_fdr_report.R:258-259` (calibrate_k_inter(...) without n_super); `R/fs_fdr_report.R:261-262` (setup_gbsg_dgm(..., n_super = n_super), formal 25000L at :188); `R/sim_aft_gbsg.R:1023-1029` (eff_at builds at the default n_super, 5000L) | YES (the mechanism is exact; the consequence is negligible). The root is searched on 5000 rows and the null built on 25000, but the root is the structural point k* = 0.5679168436; calibrating at 25000 moves k_inter by -4.76e-7, inside uniroot's tol. The visible effect is which side of the HR_H jump at k* is hit (1.0016 vs 0.9984). | Yes, on the default path of exported fs_fdr_report(fit, df) (null = c("hom", "harm"), n_super = 25000L). It is a Tier-2 reporter, off the path of the analysis actually conducted. The internal, unexported get_dgm_with_output() passes the same ... to both builders and does not have the mismatch. | At fs_fdr_report.R:258 pass n_super = n_super to calibrate_k_inter(); its ... already reaches .create_gbsg_dgm_(), and the seed defaults already match. Consistency only; calibration becomes 4.32x slower (2.9 s -> 12.5 s). Alternative: a one-line roxygen note that the root is structural and independent of n_super. | behaviour (consistency only; the numeric effect is negligible) | No, practically: k_harm moves by -4.76e-7 (still prints 0.568) and lin_pred_1 by at most 3.2e-7; the harm-null HR_H goes 1.002 -> 0.998 (side-of-jump artefact). FDR A and B (incl. A ~ 0.278) can move by at most 0.003 (3 of 1000 y_sim orderings differ); hom null untouched; full report not re-run. |
| INV-5 | `R/subgroup_consistency_main.R:497-499` (set.seed once); `R/subgroup_consistency_main.R:788-871` (plain sequential loop); `R/subgroup_consistency_main.R:963-999` (per-candidate L'Ecuyer streams on the batched path) | YES. The plain loop draws all literal-split candidates from one stream seeded once at :498; a non-sequential plan gives each candidate its own L'Ecuyer-CMRG stream. So Pcons (and possibly admission and selection) differ between plans whenever literal splitting runs (split, IRD under "resample", NA closed form). Exception: sg_focus = "maxeff" uses the same loop under every plan. | Yes, with ordinary inputs, but not forestsearch()'s default path: forestsearch()/subgroup.consistency() with parallel_args = list(plan = "sequential") vs the default "multisession", a direct subgroup.consistency() call (default parallel_args = list()), or missing workers; combined with split, count IRD, or an NA closed-form candidate (rare on the default path). Bootstrap, CV, simulation, FDR and calibration helpers force sequential internally. | In the sequential branch only (:788-871), when seed is non-NULL build candidate_seeds <- .make_candidate_rng_seeds(seed, n_candidates) and install candidate_seeds[[m]] (L'Ecuyer-CMRG) before each candidate's evaluation, restoring kind and state on exit; the maxeff path needs no change. Documentation-only minimum: correct :263-265 and :714-715 ("identical on every path") and add a sentence to forestsearch()'s parallel_args/consistency_method docs. | behaviour | No for forestsearch()'s default call (survival, "resample", multisession): a closed form draws no RNG and multisession results stay as they are. Only when the plain loop and literal splitting both run do Pcons (and possibly admission and selection) move to the batched values; maxeff is unaffected; the doc fix changes nothing. |
| DINA tau_sign | `R/forestsearch_main.R:3239-3250` (use_dina screening dina_subgroup() call, no tau_sign); compare `R/forestsearch_helpers.R:1539` (tau_sign computed) and `R/forestsearch_helpers.R:1590-1602` (dina route passes it at :1601) | YES. The screening call passes no tau_sign, so dina_subgroup() uses its default of 1, while the subgroup_method = "dina" route passes .dina_tau_sign(outcome_type, adverse_outcome). For continuous outcomes, and binary with adverse_outcome = FALSE, the screening floor, sg_focus ranking and selected cut use the un-flipped tau-hat. | Yes, via forestsearch(), not the default path (use_dina defaults to FALSE): needs use_dina = TRUE, subgroup_method = "consistency", dina_cuts = NULL and selected_only = TRUE (default), plus either continuous with adverse_outcome NULL/FALSE or binary with a ratio measure and adverse_outcome = FALSE (RD is refused, giving no DINA cuts). Survival and count are unaffected; selected_only = FALSE never reaches this call. | Add tau_sign = .dina_tau_sign(outcome_type, adverse_outcome) to the dina_subgroup() call at :3239-3250 (adverse_outcome is already resolved). This does not cover the missing "effect" re-selection on the screening route or the unoriented selected_only = FALSE pool from dina_frontier() (A-020); both are separate. | behaviour | No: use_dina = FALSE skips the block. Under the trigger conditions tau-hat is negated before the m_diff floor and the ranking, so the screening cut(s) can change (reproducer: a direction flip on continuous data); the final subgroup changes only if the consistency stage then admits a different candidate. |
| NA members (a) | `R/run_simulation_analysis.R:1126-1136` (stop at :1133, or :1134 when every NA member is truly Hc); NA enters via `:1050` / `:1059`; GRF sibling first stops at `:1466` / `:1454` (`:1442` sum lacks na.rm) | YES. An NA in sg_hat on a truly-H row makes tp and fn NA, so `if ((tp + fn) > 0)` at :1133 stops ("missing value where TRUE/FALSE needed"); NA only on truly-Hc rows stops at :1134. It needs at least one member and the harm column; both stops were measured by direct calls to the internal extractor. | Not the default path, and narrow: run_simulation_analysis() reaches it only via a non-default methods = "grf-*" arm with a GLM outcome (glm_dgm) and NA covariates in the simulated data; the stop then aborts the whole call. Default FS/FSlg flags are 0/1; survival grf-* stops earlier and is caught; DINA rejects NA covariates. NA membership measured with forestsearch(subgroup_method = "grf", binary) (28/500), not via run_simulation_analysis(). | Add na.rm = TRUE to the four sum() calls for tp/fp/tn/fn at :1130-1131, so NA subjects drop out as they already do for size.H/size.Hc and the runner completes. For the GRF sibling, na.rm = TRUE is needed at :1442-1443 as well as :1537-1538. Alternative (a method choice): treat NA membership as not a member. | behaviour | No: the default FS/FSlg membership flag is 0/1 by construction (measured: no NA with NA covariates). Only replicates whose sg_hat contains NA change, and those stop today. |
| CV fold NA (b) | `R/forestsearch_cross_validation.R:1363` (n_sgfound counts NA); `R/forestsearch_cross_validation.R:1409-1412` (two-level branch; stop at `tabit[2, 2]`, :1412), in forestsearch_KfoldOut() | YES. With fold assignments {1, NA} or {0, NA}, length(unique()) is 2 at :1363 but table() drops NA, so tabit has one row and tabit[2, 2] stops ("subscript out of bounds"); both shapes measured. A {0, 1} fold set against a single-level full-data assignment would stop at the same line. | Yes, from exported CV functions on their default-argument path, but only with NA covariates (not with complete data): default forestsearch() fits complete, fold predictions give NA, and the stop needs a full-data subgroup with every non-NA held-out subject on one side (narrow for any K, not more likely with LOO). KfoldOut stops; Kfold catches it and sets summaries to NA; tenfold aborts. The usual {0, 1, NA} shape silently gives NA sens_Hc/ppv_Hc. | No-side-effect fix: build tabit before the branch and take the two-level branch only when `n_sgfound == 2 && nrow(tabit) == 2L` (optionally `&& ncol(tabit) == 2L`), so {1, NA} and {0, NA} fall to Rule C. One-line alternative at :1363, `length(unique(stats::na.omit(df_CV$treat.recommend)))`, gives the same values for those shapes but also changes {0, 1, NA} (NA -> values). | behaviour | No for complete covariates under either fix. With the guard variant only the shapes that stop today change; with the :1363 variant, default-argument CV on NA-covariate data also changes ({0, 1, NA} sens_Hc/ppv_Hc go from NA to values). |

**Reproducers** (cap: 2 minutes each, 15 minutes in total)
| ID | Investigator reproducer: what, wall clock | Verifier reproducer: what, wall clock |
|---|---|---|
| A-024 | none (reading settled it; cites the audit's V060, not re-run) | Synthetic count data (N = 1200; treated follow-up halved vs equal), load_all, sequential: DINA tau-hat vs offset log IRR and forestsearch(subgroup_method = "dina") selections; 19.94 s (load 2.4 s, analysis 17.4 s) |
| A-042 | none (reading) | n = 300, no signal: lasso_selection() and get_FSdata(use_lasso = TRUE) on a continuous outcome, plus the survival comparator; 3.10 s whole process (2.93 s inside R) |
| A-047 | none (cites the audit's V086, 2.6 s, not re-run) | n = 200 synthetic: GRF, DINA, median-cut and replace_med_grf overlap cases; 2.77 s outer (2.6 s in-script) |
| A-055 | make_effect_estimator(binary, RD), n = 40 per arm, five p0/p1 settings; 2.58 s | Independent re-run, adding both-arms-degenerate, adverse_outcome = FALSE and interior-formula cases; 2.61 s process (2.44 s in R) |
| A-062 | none (cites the audit's V062 and V100, not re-run) | n = 600 survival: three forestsearch(subgroup_method = "grf") runs (defaults; dmin.grf = max DR effect + 1; simulated fix); 5.26 s (load_all 2.3 s) |
| A-064 | n = 800 Poisson counts: direct grf.subg.harm.glm() tree and frontier calls; 3.26 s | n = 1000, seeds 11 and 12: direct calls, emulated fix and forestsearch(count, "grf") runs; 5.65 s (/usr/bin/time; 5.45 s in-script) |
| A-103 | none (cites the audit's V141, 2.3 s, not re-run) | Synthetic labels: CV_sgs on five fold patterns, depth-2 exhaustive check, forestsearch_KfoldOut() on a synthetic 3-fold res; 2.55 s (2.38 s inside R) |
| A-107 | none (cites the audit's V028, 3.1 s, run on the working tree, not re-run) | n = 200 survival, nsims = 1: fpr_calibration() with type omitted, survival, misspelt and d_eff_approx = 30; 4.05 s outside (3.87 s inside R) |
| A-125 | n = 600 seed trial via create_dgm_for_mrct(): default vs cleaned confounders, 12 replicates; 29.17 s (26.7 s inside R) | Same data: helper-column identities, fs_args override, case study, 15 more replicates; 34.69 s (34.5 s inside R) |
| A-147 | Base-R cut() with master's breaks and labels, package not loaded; wall clock not measured (stated only as returning at once, well under the cap) | load_all, then summarize_bootstrap_subgroups() on 8 synthetic found rows, master vs right = FALSE; 2.32 s |
| A-150 | none (reading) | One combined A-150/A-162 script, in-memory patches on gbsg; 2.79 s external (2.62 s in-script) |
| A-162 | none (reading) | Same script as A-150 (counted once, under A-150) |
| A-151 | none (cites the audit's V041, 4.5 s, not re-run) | n = 300 synthetic: default forestsearch() (fs.splits = 50), then test_hte_from_forestsearch(); 4.86 s (4.69 s inside R) |
| A-156 | Calibration at n_super 5000 vs 25000 and HR_H checks; 18.35 s (script's own "total elapsed = 18.2 s") | Two executions of one script: the first crashed on a format-string bug in the script after both calibrations, 17.83 s; the re-run took 22.82 s (40.65 s together) |
| INV-5 | none (cites the audit's V109, V110b and R03, not re-run) | n = 600 survival: forestsearch() plain loop vs batched path under split and resample; 10.5 s (load_all 2.4 s) |
| DINA tau_sign | none (reading) | n = 600 continuous: trace() on dina_subgroup() on both routes, replay with tau_sign 1 and -1; 4.86 s including R startup and load_all (2.41 s inside the script) |
| NA members (a) | none (reading) | One script shared with (b): extractor calls, na.rm-patched copy, GRF sibling, GLM and survival GRF runs; 5.5 s |
| CV fold NA (b) | none (reading) | Same script as (a) (counted once, under (a)) |

Total measured reproducer wall clock: 205.46 s (investigator 53.36 s over 4 runs + verifier 152.10 s over 17 executions of 16 scripts), counting each shared script once, both A-156 verifier executions, and the outer/whole-process figure where a row gives two (198.95 s with the in-R figures instead). Excluded: A-147's unmeasured investigator run and the audit-appendix runs the investigators cite but did not re-run. No single reproducer exceeded 120 s (longest: A-125 verifier, 34.69 s); the total is under the 15-minute (900 s) cap.

**Verifier corrections that changed a verdict**


- A-125, Holds: PARTLY -> YES. The verifier held that the claim as worded (the default set includes the region indicator and the spline helper columns) holds in full; only its consequence is narrower. Reachability stays on the formal default path (widened to inject_mrct_structure() DGMs and the fs_args override). The default-path answer stays "Only when confounders.name is left NULL", with the condition widened to a confounders.name passed through fs_args.
- A-156, Holds: PARTLY -> YES (mechanism exact, consequence negligible). Reachability (fs_fdr_report() default path) is unchanged. The default-path answer stays "practically no", but its FDR claim changed from "the published harm figure (A ~ 0.278) would not move" to "A and B can move by at most 0.003".
- A-047, Reachability: the investigator had "not on the default path" (forestsearch() defaults); the verifier found that the GRF merge block runs on the default path of exported run_simulation_analysis() (run_fs_grf = TRUE, FSlg arm), where only user-supplied overlapping names are hit. The default-path answer stays "no" but is now "with default inputs", and its reason changed from "the blocks never execute" to "the GRF block executes but the default v1..v7 names drop nothing".
- CV fold NA (b), Default-path answer: the investigator had an unqualified "no" (treating {0, 1, NA} as a non-default shape). The verifier showed {0, 1, NA} is what default-argument CV ordinarily produces on NA-covariate data, so "no" holds for complete covariates and for the guard variant, but the :1363 one-liner changes sens_Hc/ppv_Hc from NA to values; the guard variant became the recommended fix. Holds and reachability are unchanged (the LOO remark was reversed).

No verdict change in these three columns, refinements only: A-024, A-042 (adds the run_simulation_analysis() route, still not the default path), A-055, A-062 (YES kept, with a qualification), A-064, A-103, A-107 (fpr_corrected dropped from the changed outputs; still "Only when"), A-147 (default-path condition narrowed to Pcons = 0.95), A-150, A-162, A-151, INV-5 (adds the direct subgroup.consistency() default; still not forestsearch()'s default path), DINA tau_sign, NA members (a) (route narrowed to GLM grf-* arms; still not the default path).

## 6. Part 4 — tests and check

**Per fix.** Each fix's own test file was run on the branch before its commit: 0 failures for all seven (§3).

**Full suite and `devtools::check()`, once, after Parts 1 and 2:**

| Run | Baseline (`REPORT_empty_set_nofind_2026-11-11.md` §6, not re-run) | This branch (`c40d5309`) | Result |
|---|---|---|---|
| Full suite (`devtools::test()`, `NOT_CRAN=true`, in the repository) | 613 blocks, 81 files, 8,270 expectations, 0 failures, 3 skips, 38 warnings; 645 s | **641 blocks, 87 files, 8,455 expectations, 0 failures, 0 errors, 3 skips, 38 warnings**; wall **695 s** (21:19:11–21:30:46; cap 30 min) | **PASS** |
| Difference | — | +28 blocks, +6 files, +185 expectations. These are exactly the six new test files: 3 + 5 + 4 + 5 + 3 + 8 blocks and 24 + 31 + 40 + 48 + 22 + 20 expectations. The skips are the same 3 (`test-cross-outcome-parity.R` ×2, `test-cv-no-subgroup-edges.R`). Warnings are unchanged. | |
| `devtools::check(document = FALSE, error_on = "never")` | 0 errors / 0 warnings / 0 notes | **0 errors / 0 warnings / 0 notes**, Status OK; wall **1,067 s** (21:19:11–21:36:58; tests 10 m; examples 36 s; `--run-donttest` 64 s; vignettes 79 s; cap 60 min) | **PASS** |

- **Check route.** As in `REPORT_empty_set_nofind_2026-11-11.md` §6: an export of the working tree (`git ls-files -co --exclude-standard`, 9,660 files at `c40d5309`) in the session scratchpad. The 80 files the branch changed were `cmp`-identical to the repository's. `forestsearch_0.4.1.9007.tar.gz` was built there.
- **Contention.** The suite, `check()` and the Part 2 render started at the same time, so each wall clock includes contention from the others.
- **Pre-set decision.** No existing expectation failed, so none needed updating. No fix was reverted.

## 7. Gates (Part 5)

Checked on the staged tree just before the report commit, and again after it (§ the final printout of this run):

- **`git status --porcelain` is empty** after the commit (no untracked or modified file is left; the scratch worktrees were removed and `git worktree list` shows only the repository).
- **The two sibling repositories are unchanged:**
  - fs-post-selection `5d165acc` on `main` and fs-glms-interpretable `dc94a17e` on `main`, each with an empty `git status --porcelain`, both before (20:08) and after.
  - No file in either repository is newer than this run's start marker.
- **`git diff e92d3c96 --stat` lists only files this report names.** Appendix B lists every path, and each belongs to a section above:
  - F1–F7: `R/`, `man/`, `tests/testthat/test-*.R`, the 7 quarto drivers, `DESCRIPTION`, `NEWS.md`;
  - Part 2: the article and `dev/reports/paper_sims_rerun_2026-11-11b/`;
  - this report and its Part 3 detail.
- **Copied** to `~/Downloads/REPORT_cran_prep_2026-11-11.md`.
- Nothing pushed or merged.

## OPEN ITEMS

1. **F2 edge case (Larry's decision).** The setting: binary RD or count IRD, with `dina_args = list(family = "gaussian")` forced, and a supplied ratio-scale threshold (e.g. `effect.threshold = 1.5`) that the package warns about and remaps to 0.05.
   - The fit's DINA floor stays the supplied 1.5, while the object records and replays 0.05. This is the same fit/replay mismatch as A-001, unchanged from the start SHA.
   - The task's two F2 clauses conflict here ("the floor is the resolved threshold" against "with either threshold supplied … unchanged"); this branch follows the second.
   - To make fit = record = replay in this case too, drop the `user_set_threshold` guard at `R/forestsearch_main.R:2784-2785` and `:3250-3252` (the agent's first version), and flip the pinning test block (`test-dina-gaussian-floor.R:123`).
2. **F2 exposure.** The pass/fail checks in `quarto/dina/smoke-tests/dina_screening_smoketest.qmd:112, :126-127, :134-135` (continuous, `use_dina = TRUE`, both thresholds unset) may now come out differently. Not re-rendered.
3. **F1: listed, not edited** (each assumes the old raw orientation for `adverse_outcome = FALSE` fits):
   - The "RAW scale" decomposition columns `fb_src1`/`fb_src2`, now on the parent's scale. They appear in the recorder comments, in the captions at `batch_1_100.qmd:1304, :1306`, `batch_1_20.qmd:838, :840`, `batch_1_50.qmd:1586, :1588` and the matching lines in the four `mr_*` drivers, and in the prose at `:1295`, `:829`, `:1577` and `batch_1_50.qmd:129-141`.
   - The algebra check in `quarto/simulations/actg175/continuous/o2_fb_mr_batches.R:53-54` (after the fix: `(nv - fb) - (src1 + src2)`).
   - Tracked renders and payloads: `sim_fs_maxeffCons_fb_mr_md40_knoise0_n500_batch_1_{100,50}.html` and 4 `.rds` under `quarto/simulations/actg175/continuous/fb_mr_md_harm/` (not re-run). `fb_H_*`/`fb_Hc_*` are unaffected.
   - `dev/glm-continuous-sims/oc_breadth_stage2_2026-08-31_driver_purl.R:680-683, :738-741` still re-orients. It is a purl behind `run_fb`, with `nb_boots` `NULL`.
   - Package panels that build their effect with the default `adverse_outcome = TRUE` and now disagree in sign with the bias-corrected rows for FALSE parents:
     - `R/plot_sg_glm_outcomes.R:254-262` and `R/plot_subgroup_results_forestplot.R:431-439`;
     - the `FSsg_tab` `Diff` column (`R/summary_utility_functions.R:255`, raw `rate1 − rate0`);
     - `mr_estimates_table()`'s `boot_results = NULL` route (`R/mr_estimates_table.R:204-208`).
   - Scope call: the four `sim_fs_maxeffCons_mr_*.qmd` drivers outside the named glob were edited too (same recorder, currently skipped). Revert those hunks if only the named files were wanted.
4. **F4.**
   - Mirror case: when every subject is truly in H and everyone is identified, Rule F gives npv 0. The exact mirror of the empty-H PPV convention would be `NA`.
   - `fs_mr_oc_summary()`'s conditional (detected-only) row still drops a detected-everyone replicate whose npv is `NA`.
   - The 1,344 idsweep rows (NPV report OPEN ITEM 6) and other stored payloads were not regenerated.
   - The quarto MR harness `.classify()` still records `NA` for such a replicate; `fs_mr_oc_summary()` now fills it with 0.
5. **F5.**
   - `gbsg <- NULL` at `R/synthetic_data_perturbation.R:365` is now dead (kept for a one-line diff).
   - The example is still `\dontrun{}`.
   - `grade` is never perturbed (A-162, Part 3).
   - The same `data(gbsg, …)` pattern remains in `\dontrun` examples (`R/synthetic_data_perturbation.R:438`, `R/pareto_frontier_table.R:124`, `R/calibrate_k_treat.R:55`) and in `vignettes/articles/extreme_subgroups.qmd:28`.
6. **F6.**
   - INV-2b: `future::plan()` returns only the top level, so a nested plan stack comes back flattened. This is the package's existing convention and was not changed.
   - `R/globals.R:30` `".Last.future.plan"` is now unreferenced.
   - `progressr::handlers("progress")` at `R/mrct_simulation.R:394` is global and not restored.
7. **F7.** The 20 `NULL`-formal sites were edited under the literal reading of the rule. They are a separable group and can be reverted as one if Larry prefers the narrower reading.
8. **Fixed seed in `subgroup.consistency()`.** Out of scope here. It is Part 1 of `TASK_cran_prep_seed_legacy_2026-11-11.md`.
9. **Part 3.** All 18 rows hold, and none changes a `forestsearch()` default-path result. The fixes are described in §5 and nothing was changed.
10. **Scratch.** The fix agents' worktrees, patches and notes (`fixnotes/F*.md`, with per-file test tables) are in this session's scratchpad. They are not in the repository. The worktrees were removed from `git worktree list` before the final gates.

## Appendix A — F7: every line changed under `R/` (commit `c47e26bf`, generated from `git diff -U0 56b2429e c47e26bf -- R/`)

| # | old (file:line) | old text | new (file:line) | new text |
|---|---|---|---|---|
| 1 | R/Cox_estimation_helpers.R:8 | `` #' @param cox_initial Optional pre-fitted Cox model object to use instead of `` | R/Cox_estimation_helpers.R:8 | `` #' @param cox_initial Numeric. Starting value for the treatment coefficient, `` |
| 2 | R/Cox_estimation_helpers.R:9 | `` #'   fitting a new model. Default NULL `` | R/Cox_estimation_helpers.R:9 | `` #'   passed to `coxph()` as `init` when the formula has a single term. `` |
| 3 | (none) | (none) | R/Cox_estimation_helpers.R:10 | `` #'   Default `log(1)` `` |
| 4 | R/bootstrap_analysis_dofuture.R:79 | `` #'   first three bootstrap iterations for debugging purposes. Default: \code{FALSE}. `` | R/bootstrap_analysis_dofuture.R:79 | `` #'   first three bootstrap iterations for debugging purposes. No default; must be supplied. `` |
| 5 | R/calibrate_glm_interaction.R:48 | `` #'   \code{k_inter}. Default \code{c(0, 10)}. `` | R/calibrate_glm_interaction.R:48 | `` #'   \code{k_inter}. Default \code{c(-10, 10)}. `` |
| 6 | R/cox_spline_fit.R:537 | `` #'   Default is 1 (HR=1 on log scale is 0). `` | R/cox_spline_fit.R:537 | `` #'   Default is 0 (on the log scale, HR=1). `` |
| 7 | R/dina_subgroup.R:168 | `` #'   Default `60L`, matching the convention in GRF and `forestsearch()`. `` | (deleted) | (deleted) |
| 8 | R/forestsearch_cross_validation.R:190 | `` #' @param Kfolds Integer. Number of folds (default: \code{nrow(fs.est$df.est)} for LOO). `` | R/forestsearch_cross_validation.R:190 | `` #' @param Kfolds Integer. Number of folds (default: \code{nrow(fs.est$df.est)} for LOO; `` |
| 9 | (none) | (none) | R/forestsearch_cross_validation.R:191 | `` #'   when no subgroup was found, the number of rows of the original analysis data). `` |
| 10 | R/forestsearch_main.R:382 | `` #' @param est.scale Character. Estimation scale ("hr" or "rmst"). Default "hr". `` | R/forestsearch_main.R:382 | `` #' @param est.scale Character. Estimation scale ("hr" or "1/hr"). Default "hr". `` |
| 11 | R/forestsearch_main.R:928 | `` #'   improved performance. Default FALSE for backward compatibility. When TRUE, `` | R/forestsearch_main.R:928 | `` #'   improved performance. Default TRUE. When TRUE, `` |
| 12 | R/forestsearch_main.R:978 | `` #'     \item{\code{ci_method}}{\code{"field"} (default) additionally runs the `` | R/forestsearch_main.R:978 | `` #'     \item{\code{ci_method}}{\code{"field"} (default under `` |
| 13 | (none) | (none) | R/forestsearch_main.R:979 | `` #'       \code{subgroup_method = "consistency"}; the \code{"dina"} and `` |
| 14 | (none) | (none) | R/forestsearch_main.R:980 | `` #'       \code{"grf"} methods default to \code{"ij"}) additionally runs the `` |
| 15 | R/forestsearch_main.R:1001 | `` #'     \item{\code{field_complement}}{Logical, default \code{FALSE}.  Under `` | R/forestsearch_main.R:1003 | `` #'     \item{\code{field_complement}}{Logical, default \code{TRUE}.  Under `` |
| 16 | R/forestsearch_main.R:1006 | `` #'       \code{\link{fs_mr_inference}}).  Add-only pass-through: the default `` | R/forestsearch_main.R:1008 | `` #'       \code{\link{fs_mr_inference}}).  \code{FALSE} omits the complement's `` |
| 17 | R/forestsearch_main.R:1007 | `` #'       reproduces prior output exactly.} `` | R/forestsearch_main.R:1009 | `` #'       field block (the prior behaviour).} `` |
| 18 | R/forestsearch_main.R:1064 | `` #' @param adverse_outcome Logical or \code{NULL}. If \code{TRUE}, higher `` | R/forestsearch_main.R:1066 | `` #' @param adverse_outcome Logical or \code{NULL} (default, resolved by `` |
| 19 | (none) | (none) | R/forestsearch_main.R:1067 | `` #'   outcome type). If \code{TRUE}, higher `` |
| 20 | R/forestsearch_main.R:1108 | `` #'   \code{"none"}.  Default: \code{"none"} for RCT, \code{"grf"} for `` | R/forestsearch_main.R:1111 | `` #'   \code{"none"}.  Default: \code{NULL}, which resolves to \code{"none"} `` |
| 21 | R/forestsearch_main.R:1109 | `` #'   observational. `` | R/forestsearch_main.R:1112 | `` #'   for RCT and \code{"grf"} for observational. `` |
| 22 | R/forestsearch_main.R:1221 | `` #'       the infinitesimal-jackknife analogue of the FB interval under the `` | R/forestsearch_main.R:1224 | `` #'       the infinitesimal-jackknife analogue of the FB interval under `` |
| 23 | R/forestsearch_main.R:1222 | `` #'       default \code{ci_method = "field"} and under \code{"ij"}, and the `` | R/forestsearch_main.R:1225 | `` #'       \code{ci_method = "field"} (the default under `` |
| 24 | R/forestsearch_main.R:1223 | `` #'       subgroup robust-SE interval under \code{"wald"}.  Under the default `` | R/forestsearch_main.R:1226 | `` #'       \code{subgroup_method = "consistency"}) and under \code{"ij"}, and the `` |
| 25 | (none) | (none) | R/forestsearch_main.R:1227 | `` #'       subgroup robust-SE interval under \code{"wald"}.  Under \code{"field"} `` |
| 26 | R/forestsearch_methods.R:214 | `` #' (`ci_method = "field"`, the default), a further block reports: the `` | R/forestsearch_methods.R:214 | `` #' (`ci_method = "field"`, the default under `subgroup_method = "consistency"`), `` |
| 27 | (none) | (none) | R/forestsearch_methods.R:215 | `` #' a further block reports: the `` |
| 28 | R/fpr_approximation.R:57 | `` #'   besides \code{c1}.  Default \code{c(1.0, 1.5, 2.0)}. `` | R/fpr_approximation.R:57 | `` #'   besides \code{c1}.  Default \code{NULL} (only \code{c1} is evaluated). `` |
| 29 | R/fs_fdr_report.R:76 | `` #'   \code{nrow(df_analysis)}.  The published reference baseline used `` | R/fs_fdr_report.R:76 | `` #'   \code{NULL}, giving \code{nrow(df_analysis)}.  The published reference baseline used `` |
| 30 | R/fs_identification_summary.R:127 | `` #'   involvement vector. Defaults to the three named roles. `` | R/fs_identification_summary.R:127 | `` #'   involvement vector. Default `NULL`, which uses the three named roles. `` |
| 31 | R/fs_identification_summary.R:302 | `` #' @param roles Character vector. Which covariates to draw; defaults to the `` | R/fs_identification_summary.R:302 | `` #' @param roles Character vector. Which covariates to draw; default `NULL`, which uses the `` |
| 32 | R/fs_oc_grid.R:232 | `` #'   argument is ignored).  Defaults from \code{forestsearch_args}. `` | R/fs_oc_grid.R:232 | `` #'   argument is ignored).  Default \code{NULL}, which takes the value from \code{forestsearch_args}. `` |
| 33 | R/fs_oc_predict.R:98 | `` #'   \code{hr.threshold}).  Default `` | R/fs_oc_predict.R:98 | `` #'   \code{hr.threshold}).  Default \code{NULL}, which takes `` |
| 34 | R/fs_oc_predict.R:104 | `` #'   \code{forestsearch_args$consistency.threshold}; an explicit value `` | R/fs_oc_predict.R:104 | `` #'   \code{NULL}, which takes \code{forestsearch_args$consistency.threshold}; `` |
| 35 | R/fs_oc_predict.R:105 | `` #'   overrides. `` | R/fs_oc_predict.R:105 | `` #'   an explicit value overrides. `` |
| 36 | R/fs_oc_predict.R:116 | `` #'   \code{forestsearch_args$pconsistency.threshold}, itself defaulting to `` | R/fs_oc_predict.R:116 | `` #'   \code{NULL}, which takes \code{forestsearch_args$pconsistency.threshold}, `` |
| 37 | (none) | (none) | R/fs_oc_predict.R:117 | `` #'   itself defaulting to `` |
| 38 | R/gg_forest.R:80 | `` #' @param tick_labels Character vector. Custom tick labels (default: as.character(ticks_at)). `` | R/gg_forest.R:80 | `` #' @param tick_labels Character vector. Custom tick labels (default: NULL, giving as.character(ticks_at)). `` |
| 39 | R/gg_forest.R:92 | `` #'   c(label, ci, annot_1, annot_2, …). Default: c(3.5, 5, rep(1, n_annot)). `` | R/gg_forest.R:92 | `` #'   c(label, ci, annot_1, annot_2, …). Default: NULL, giving c(3.5, 5, rep(1.1, n_annot)). `` |
| 40 | R/grf_subg_harm_glm.R:55 | `` #'   Default is \code{"log_OR"} for binary, \code{"MD"} for continuous, `` | R/grf_subg_harm_glm.R:55 | `` #'   Default \code{NULL}, which resolves to \code{"log_OR"} for binary, `` |
| 41 | R/grf_subg_harm_glm.R:56 | `` #'   and \code{"log_IRR"} for count. `` | R/grf_subg_harm_glm.R:56 | `` #'   \code{"MD"} for continuous, and \code{"log_IRR"} for count. `` |
| 42 | R/leff_calibration_io.R:123 | `` #'   Default: \code{"_data"}. `` | R/leff_calibration_io.R:123 | `` #'   Default: \code{NULL}, which resolves to \code{"_data"}. `` |
| 43 | R/mrct_inject_dgm.R:205 | `` #' @param scale_ref data.frame used to define the scaling (defaults to pop) `` | R/mrct_inject_dgm.R:205 | `` #' @param scale_ref data.frame used to define the scaling (default NULL, meaning pop) `` |
| 44 | R/mrct_simulation.R:54 | `` #'   Default: c("z_age <= 65", "z_bm <= 0", "z_bm <= 1", "z_bm <= 2", "z_bm <= 5") `` | R/mrct_simulation.R:54 | `` #'   Default: NULL, which resolves to `` |
| 45 | (none) | (none) | R/mrct_simulation.R:55 | `` #'   c("z_age <= 65", "z_bm <= 0", "z_bm <= 1", "z_bm <= 2", "z_bm <= 5") `` |
| 46 | R/mrct_simulation.R:1284 | `` #'   statistics. Default: TRUE when \code{trim_fraction} is non-NULL, `` | R/mrct_simulation.R:1285 | `` #'   statistics. Default: \code{NULL}, which resolves to \code{TRUE} when `` |
| 47 | R/mrct_simulation.R:1285 | `` #'   FALSE otherwise. `` | R/mrct_simulation.R:1286 | `` #'   \code{trim_fraction} is non-NULL and \code{FALSE} otherwise. `` |
| 48 | R/oc_analyses.R:405 | `` #' @param analyses Character vector. Analysis methods to include. Default: all `` | R/oc_analyses.R:405 | `` #' @param analyses Character vector. Analysis methods to include. Default: NULL, meaning all analyses in results `` |
| 49 | R/plot_effect_distribution.R:59 | `` #'   Default: \code{c(itt = "ITT", H = "Identified H-hat", `` | R/plot_effect_distribution.R:59 | `` #'   Default: \code{NULL}, giving \code{c(itt = "ITT", H = "Identified H-hat", `` |
| 50 | R/plot_pareto_frontier.R:23 | `` #'   defaults to the value used in the original \code{forestsearch} `` | R/plot_pareto_frontier.R:23 | `` #'   default \code{NULL}, which uses the value from the original \code{forestsearch} `` |
| 51 | R/plot_subgroup_results_forestplot.R:70 | `` #' @param outcome_type Character or \code{NULL}. One of \code{"survival"} `` | R/plot_subgroup_results_forestplot.R:70 | `` #' @param outcome_type Character or \code{NULL} (default). One of `` |
| 52 | R/plot_subgroup_results_forestplot.R:71 | `` #'   (default), \code{"binary"}, or \code{"continuous"}.  If \code{NULL}, `` | R/plot_subgroup_results_forestplot.R:71 | `` #'   \code{"survival"}, \code{"binary"}, or \code{"continuous"}.  If \code{NULL}, `` |
| 53 | R/plot_subgroup_results_forestplot.R:72 | `` #'   auto-detected from \code{fs_results$fs.est$outcome_type}. `` | R/plot_subgroup_results_forestplot.R:72 | `` #'   auto-detected from \code{fs_results$fs.est$outcome_type}, or `` |
| 54 | (none) | (none) | R/plot_subgroup_results_forestplot.R:73 | `` #'   \code{"survival"} when that is absent. `` |
| 55 | R/run_subgroup_sims.R:200 | `` #'   [generate_glm_dgm()] default and the calibration target; note that `` | R/run_subgroup_sims.R:200 | `` #'   [generate_glm_dgm()] default and the calibration target; it is also `` |
| 56 | R/run_subgroup_sims.R:201 | `` #'   [make_effect_estimator()]'s own binary default is `"RD"`), or `` | R/run_subgroup_sims.R:201 | `` #'   [make_effect_estimator()]'s own binary default), or `` |
| 57 | R/simulation_tables.R:1053 | `` #' @param scenario Character. One of \code{"null"} or \code{"alt"} (default). `` | R/simulation_tables.R:1053 | `` #' @param scenario Character or \code{NULL} (default). One of \code{"null"} or \code{"alt"}. `` |
| 58 | R/subgroup_consistency_main.R:276 | `` #'       Default: When stop_threshold is set and sg_focus is "hr" or "minSG", `` | R/subgroup_consistency_main.R:276 | `` #'       Default: When stop_threshold is set and sg_focus is "maxeffCons", `` |
| 59 | R/synthetic_data_perturbation.R:6 | `` #' @param n Number of synthetic observations to generate (default: same as original) `` | R/synthetic_data_perturbation.R:6 | `` #' @param n Number of synthetic observations to generate (default: NULL, meaning same as original) `` |
| 60 | R/synthetic_data_perturbation.R:10 | `` #' @param cat_flip_prob Probability of flipping categorical values (default: noise_level/2) `` | R/synthetic_data_perturbation.R:10 | `` #' @param cat_flip_prob Probability of flipping categorical values (default: NULL, giving noise_level/2) `` |
| 61 | R/truefind_asymptotic_glm.R:180 | `` #'   Default: 1.0 for ratio, 0.0 for difference. `` | R/truefind_asymptotic_glm.R:180 | `` #'   Default: 1.0 on both scales; supply \code{c2} for difference measures. `` |

## Appendix B — files changed since the start SHA (`git diff --stat e92d3c96 HEAD`)

100 paths (A = added, M = modified), from `git diff --name-status e92d3c96` on the tree of the report commit:

| Status | Path |
|---|---|
| M | `DESCRIPTION` |
| M | `NEWS.md` |
| M | `R/Cox_estimation_helpers.R` |
| M | `R/bootstrap_analysis_dofuture.R` |
| M | `R/bootstrap_dofuture_main.R` |
| M | `R/calibrate_glm_interaction.R` |
| M | `R/cox_spline_fit.R` |
| M | `R/dina_subgroup.R` |
| M | `R/forestsearch_cross_validation.R` |
| M | `R/forestsearch_helpers.R` |
| M | `R/forestsearch_main.R` |
| M | `R/forestsearch_methods.R` |
| M | `R/fpr_approximation.R` |
| M | `R/fpr_calibration.R` |
| M | `R/fs_fdr_report.R` |
| M | `R/fs_identification_summary.R` |
| M | `R/fs_mr_oc_summary.R` |
| M | `R/fs_oc_grid.R` |
| M | `R/fs_oc_predict.R` |
| M | `R/gg_forest.R` |
| M | `R/grf_main.R` |
| M | `R/grf_subg_harm_glm.R` |
| M | `R/leff_calibration_io.R` |
| M | `R/mrct_inject_dgm.R` |
| M | `R/mrct_simulation.R` |
| M | `R/oc_analyses.R` |
| M | `R/plot_effect_distribution.R` |
| M | `R/plot_pareto_frontier.R` |
| M | `R/plot_subgroup_results_forestplot.R` |
| M | `R/run_simulation_analysis.R` |
| M | `R/run_subgroup_sims.R` |
| M | `R/simulation_tables.R` |
| M | `R/subgroup_consistency_main.R` |
| M | `R/subgroup_search.R` |
| M | `R/synthetic_data_perturbation.R` |
| M | `R/truefind_asymptotic_glm.R` |
| A | `dev/reports/REPORT_cran_prep_2026-11-11.md` |
| A | `dev/reports/REPORT_cran_prep_2026-11-11_part3_detail.md` |
| A | `dev/reports/paper_sims_rerun_2026-11-11b/accounting/accounting_grouped.csv` |
| A | `dev/reports/paper_sims_rerun_2026-11-11b/accounting/accounting_tables.md` |
| A | `dev/reports/paper_sims_rerun_2026-11-11b/accounting/accounting_unpaired_lines.csv` |
| A | `dev/reports/paper_sims_rerun_2026-11-11b/accounting/group_lines.txt` |
| A | `dev/reports/paper_sims_rerun_2026-11-11b/accounting/paper_sims_accounting.R` |
| A | `dev/reports/paper_sims_rerun_2026-11-11b/accounting/paper_sims_groups.R` |
| A | `dev/reports/paper_sims_rerun_2026-11-11b/logs/install.log` |
| A | `dev/reports/paper_sims_rerun_2026-11-11b/logs/library_gate.txt` |
| A | `dev/reports/paper_sims_rerun_2026-11-11b/logs/provenance.txt` |
| A | `dev/reports/paper_sims_rerun_2026-11-11b/logs/render_pkgdown_paper_simulations.log` |
| A | `dev/reports/paper_sims_rerun_2026-11-11b/logs/render_timing.txt` |
| A | `dev/reports/paper_sims_rerun_2026-11-11b/paper_simulations.html` |
| A | `dev/reports/paper_sims_rerun_2026-11-11b/paper_simulations_files/figure-html/km-plot-1.png` |
| A | `dev/reports/paper_sims_rerun_2026-11-11b/paper_simulations_files/figure-html/theoretical-detection-rates-1.png` |
| A | `dev/reports/paper_sims_rerun_2026-11-11b/results/keyvalues.txt` |
| A | `dev/reports/paper_sims_rerun_2026-11-11b/results/paper_simulations.txt` |
| A | `dev/reports/paper_sims_rerun_2026-11-11b/results/paper_simulations_output.txt` |
| A | `dev/reports/paper_sims_rerun_2026-11-11b/strip.sh` |
| M | `man/SGplot_estimates.Rd` |
| M | `man/bootstrap_results.Rd` |
| M | `man/calibrate_glm_interaction.Rd` |
| M | `man/compute_detection_probability_glm.Rd` |
| M | `man/dina_subgroup.Rd` |
| M | `man/forestsearch.Rd` |
| M | `man/forestsearch_Kfold.Rd` |
| M | `man/fpr_approx.Rd` |
| M | `man/fs_fdr_report.Rd` |
| M | `man/fs_identification_structure.Rd` |
| M | `man/fs_oc_grid.Rd` |
| M | `man/fs_oc_invert.Rd` |
| M | `man/fs_oc_predict.Rd` |
| M | `man/generate_bootstrap_synthetic.Rd` |
| M | `man/get_Cox_sg.Rd` |
| M | `man/gg_forest.Rd` |
| M | `man/grf.subg.harm.glm.Rd` |
| M | `man/interpret_estimation_table.Rd` |
| M | `man/load_leff_calibration.Rd` |
| M | `man/make_region_model.Rd` |
| M | `man/mrct_region_sims.Rd` |
| M | `man/plot_effect_distribution.Rd` |
| M | `man/plot_fs_identification_involvement.Rd` |
| M | `man/plot_pareto_frontier.Rd` |
| M | `man/plot_subgroup_effects.Rd` |
| M | `man/plot_subgroup_results_forestplot.Rd` |
| M | `man/print.forestsearch.Rd` |
| M | `man/subgroup.consistency.Rd` |
| M | `man/subgroup_glm.Rd` |
| M | `man/summarize_simulation_results.Rd` |
| M | `quarto/simulations/actg175/continuous/sim_fs_maxeffCons_fb_mr_md40_knoise0_n500_batch_1_100.qmd` |
| M | `quarto/simulations/actg175/continuous/sim_fs_maxeffCons_fb_mr_md40_knoise0_n500_batch_1_20.qmd` |
| M | `quarto/simulations/actg175/continuous/sim_fs_maxeffCons_fb_mr_md40_knoise0_n500_batch_1_50.qmd` |
| M | `quarto/simulations/actg175/continuous/sim_fs_maxeffCons_mr_md120_knoise0_n500_batch_1_1000.qmd` |
| M | `quarto/simulations/actg175/continuous/sim_fs_maxeffCons_mr_md40_knoise0_n500_batch_1_100.qmd` |
| M | `quarto/simulations/actg175/continuous/sim_fs_maxeffCons_mr_md40_knoise0_n500_batch_1_1000.qmd` |
| M | `quarto/simulations/actg175/continuous/sim_fs_maxeffCons_mr_md40_knoise0_n700_batch_1_1000.qmd` |
| A | `tests/testthat/test-bootstrap-glm-parent-scale.R` |
| A | `tests/testthat/test-dina-gaussian-floor.R` |
| A | `tests/testthat/test-gbsg-bootstrap-general.R` |
| A | `tests/testthat/test-grf-survival-matrix.R` |
| A | `tests/testthat/test-npv-empty-complement.R` |
| A | `tests/testthat/test-plan-restored.R` |
| M | `vignettes/articles/paper_simulations.qmd` |
