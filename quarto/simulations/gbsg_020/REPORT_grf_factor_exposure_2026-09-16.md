# REPORT — GRF factor-membership exposure on the survival path: read-only check

Date: 2026-09-16. Machine: `pop-os` (64 physical cores, 251 GB; R 4.6.1, reference BLAS). Branch `feature/glm-extension`. Task: `dev/tasks/TASK_grf_factor_exposure_2026-09-16.md` (committed as received, `e3e6f652`). Installed forestsearch 0.3.5, `Built: R 4.6.1; ; 2026-09-16 05:57:14 UTC; unix`; read-only: no `R/`, template, script, document or payload edit, no install, no campaign. The only computation: data regeneration and two GRF identification fits with `mr_inference = FALSE` from a temporary directory outside the repo (`/tmp/fs_gfx_R5js`, removed at closeout), under `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`. Line numbers at HEAD `e3e6f652`. Finding under test: `quarto/simulations/actg175/continuous/REPORT_md_grf_stage1_2026-09-16.md` §1.6(c) / F2 — `.grf_evaluate_subgroup()` (`R/grf_subgroup_labels.R:224–228`) compares raw data columns with numeric cuts without the factor coercion the forest's covariate matrix receives (`R/grf_subg_harm_glm.R:884–895`), so candidates on factor covariates get NA membership and are dropped from the effect re-selection and from MR's family. This record decides nothing.

## 1. Provenance — GATE PASS

```
pop-os
feature/glm-extension
a4c063bf
a4c063bf GRF on the ACTG175 continuous (MD) design, campaign mdgrf, Stage 1 record (TASK_md_grf_2026-09-16): ...
f0b9c844 scripts_mdgrf (TASK_md_grf_2026-09-16 §1.6, transplants of scripts_mdsgnb20): ...
894da993 MD template E1-E4 (TASK_md_grf_2026-09-16 §1.5, transplanted from the survival m1 template): ...
[tracked modifications: none]
GRF Stage 1 record in HEAD
[R / Rscript / quarto / deno processes, by process name: none]
R 4.6.1; ; 2026-09-16 05:57:14 UTC; unix
```
First commit: `e3e6f652 Add TASK_grf_factor_exposure_2026-09-16 as received`.

## 2. The survival simulation's covariates

**2.1 The campaign.** Template: `quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` — the one `render.sh` renders (`scripts_dinamr/render.sh:12`: `quarto render sim_fs_maxeffCons_fb_mr_field_m1_template.qmd --output "${OUT}.html"`). Runner: `scripts_dinamr/grfmr.sh` (`:19–21`, `KN=(FS_S7_METHOD=grf FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20 FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE FS_S7_FIELD_RECOV=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none FS_S7_CAMPAIGN=grfmr FS_S7_WORKERS=12)`; `:12–15`: "`dmin.grf = 0.0`, `grf_selection = "frontier"` and `grf_select_statistic = "effect"` are TEMPLATE LITERALS (template lines 503-506), not `FS_S7_*` knobs"; `:29`, `:33`: `env -u FS_S7_Z1Q -u FS_S7_ER_JCUTS $KN FS_S7_HR=$H FS_S7_N=$N FS_S7_MODE=$MODE FS_S7_START=$START FS_S7_NSIMS=$NS ./render.sh $OUT`), with `grfmrC.sh` / `grfmrC_smoke.sh` for the completion cells; cells in `scripts_dinamr/grfmr.cells` (e.g. `- 500 1.50 A124_h150_n500`, `0.60 500 1.50 A31_h150_n500`) and `grfmrC.cells`. Bundles: **present on this checkout** — `git ls-files quarto/simulations/gbsg_020/results | grep grfmr` lists 55 files (18 cells × {`res_1_1000`, `res_1001_2000`, `combined_1_2000`} + the `grfmrsmk` smoke bundle).

**2.2 Chunks run.** `knitr::purl()` of the template into `/tmp/fs_gfx_R5js/tpl.R`; chunks `setup` (`:211`) and `build-dgm` (`:736`) evaluated with the `grfmr` cell **A124_h150_n500** (`FS_S7_METHOD=grf FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20 FS_S7_HR=1.50 FS_S7_N=500 FS_S7_MODE=batch FS_S7_FB=none FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE FS_S7_FIELD_RECOV=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_WORKERS=1`; `FS_S7_Z1Q` and `FS_S7_ER_JCUTS` unset, as the runner unsets them for the 12.4% cells). The template has no separate data chunk: the replicate's data are drawn inside `record_replicate()` in `run-loop` (`:1003–1032`), so those lines were run verbatim for sim_id 1: `df <- simulate_from_dgm(dgm, n = n_sample, analysis_time = analysis_time, cens_adjust = cens_adjust, seed = seed_base + sim_id)` (`:1007–1008`), `df[[id_name]] <- seq_len(nrow(df))` (`:1011`), `confs <- intersect(confounders_base, names(df))` (`:1013`; `confounders_base <- c("er", "age", "meno", "pgr", "nodes", "size", "grade")`, `:704`), `k_random_noise` 0 (no noise columns). The template's echo: `Template knobs: method=grf hr=1.50 n=500 z1q=0.25 focus=effMaxSG nbhd=0.20 er_jcuts=10 ci_method=field … campaign=gfxprobe run_mode=batch`; the DGM: `Harm rule: er <= quantile(er, 0.25) & meno == 0; super-population prevalence 0.1242; k_inter 1.11155` (the bundle `meta`: `harm_prevalence_super 0.12418`, `target_hr_harm 1.5`, `harm_z1_quantile 0.25`).

**2.3 Covariate classes** (`str(df[confs])`, sim_id 1):
```
'data.frame':	500 obs. of  7 variables:
 $ er   : int  0 53 0 28 0 40 9 0 38 59 ...
 $ age  : int  64 55 48 65 48 42 43 54 59 65 ...
 $ meno : int  1 1 0 1 0 0 0 1 1 1 ...
 $ pgr  : int  227 26 0 11 0 77 23 1 88 42 ...
 $ nodes: int  1 2 1 6 1 8 1 3 2 3 ...
 $ size : int  17 25 45 33 45 9 27 30 23 10 ...
 $ grade: int  1 1 3 2 3 1 2 3 2 2 ...
```
Every covariate the identifier sees is `integer`; **no factor or character column among `confs`**. The simulated frame does carry factor columns (`v1`–`v7`, the DGM's internal dichotomizations), but they are not in `confounders_base` and are not passed. No chunk converts columns before the `forestsearch()` call (`:1039–1073` build `base_args` / `method_args` from `df` as drawn).

## 3. The applied data

**3.1 Documents that run GRF or DINA** (from `REPORT_actg175_applied_effmaxsg_stage0_2026-09-16.md` S0.1, rows with `use_dina` / `use_grf` or "DINA and GRF runs alongside FS"):
`actg175/analysis_actg175_binary_multimethod_fixed_family.qmd`, `actg175/analysis_actg175_binary_multimethod_frontend.qmd`, `actg175/analysis_actg175_binary_multimethod_psi_v2_2.qmd`, `actg175/analysis_actg175_binary_multimethod_psi_v3a.qmd` (DINA and GRF alongside FS; `subgroup_method` knob, standalone `grf.subg.harm.glm()` at `fixed_family:461` and `dina()` at `:1134`); `actg175/analysis_actg175_binary_sgfocus.qmd` (`use_grf = TRUE`, `:665`; `subgroup_method = "consistency"`, `:678`); `actg175/analysis_actg175_continuous_compare_all.qmd` (`fs_use_grf <- TRUE`, `:144`, passed to `compare_selection_rules()` at `:328`; the anchor at `:433` `use_grf = FALSE`); `gbsg/analysis_gbsg_survival_multimethod.qmd` (DINA and GRF alongside FS; `subgroup_method` knob `:212`, `use_grf = FALSE` `:710`); `count_data_demo.qmd` (`use_grf = TRUE`, `:672`; standalone `grf.subg.harm.glm()` `:381`); `validation_glm_simulation_study.qmd` (`use_grf = TRUE`, `:547`, `:638`). All paths under `quarto/applications/`.

**3.2–3.3 Data construction and classes, per document.**

| document | data lines (quoted) | covariate classes passed | any factor? |
|---|---|---|---|
| the four `analysis_actg175_binary_multimethod_*.qmd` | `actg_df <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))` (`fixed_family:369`); `cont_vars <- c("age", "preanti", "wtkg", "karnof", "cd40", "cd80")`, `bin_vars <- c("hemo", "homo", "drugs", "race", "gender", "symptom")` (`:348–349`); **`for (v in bin_vars) actg_df[[v]] <- as.numeric(actg_df[[v]])`** (`:407`; the chunk's comment `:397–406`: "Keep ALL candidate covariates numeric … binary indicators left as factors trigger hundreds of "'<=' not meaningful for factors" warnings … every engine (consistency, DINA, GRF, standalone DINA) receives numeric covariates"); the same lines at `frontend:400`, `psi_v2_2:385`, `psi_v3a:385` | six continuous: `integer`/`numeric` as in `speff2trial::ACTG175`; six binary: `numeric` after the coercion | **no** |
| `analysis_actg175_binary_sgfocus.qmd` | `:435` the same subset; **`for (v in bin_vars) actg_df[[v]] <- as.numeric(actg_df[[v]])`** (`:469`) | as above | **no** |
| `analysis_actg175_continuous_compare_all.qmd` | `:225` the same subset; **`for (v in bin_vars) actg_df[[v]] <- as.factor(actg_df[[v]])`** (`:251`, "Coerce binary indicators to factors") | six continuous numeric; **six binary `factor`** | **yes** — on the `use_grf = TRUE` screening path only (§3.4) |
| `gbsg/analysis_gbsg_survival_multimethod.qmd` | `library(survival)` (`:317`); `df.analysis <- gbsg` (`:500`, `survival::gbsg`); `grade3 <- ifelse(grade == "3", 1, 0)` (`:506`); `confounders.name <- c("age", "meno", "size", "grade3", "nodes", "pgr", "er")` (`:512`) | `survival::gbsg`: `age`, `meno`, `size`, `nodes`, `pgr`, `er` `integer`; `grade3` `numeric` | **no** |
| `count_data_demo.qmd` | `eos_low = as.factor(rbinom(n, 1, 0.35))`, `smoking = as.factor(…)`, `male = as.factor(…)`, `age = round(rnorm(n, 60, 12))` (`:241–244`); `confounders.name = c("eos_low", "smoking", "male", "age")` (`:381`, `:651`) | three `factor`, one `numeric` | **yes** — standalone `grf.subg.harm.glm()` (`:381`, count outcome) and `use_grf = TRUE` screening (`:672`) |
| `validation_glm_simulation_study.qmd` | `bm1 = as.factor(rbinom(n, 1, 0.70))`, `bm2`, `bm3` likewise, `age = round(rnorm(n, 55, 10))`, `ecog = as.factor(sample(0:1, …))` (`:169–173`; the other two generators `:198–200`, `:229–231` likewise) | three or four `factor`, one `numeric` | **yes** — `use_grf = TRUE` screening (`:547`, `:638`) |

**3.4 Which path each factor-carrying document exercises.** The evaluator under test is reached by `subgroup_method = "grf"` (frontier: `grf_subg_harm_glm.R:548` `data$treat.recommend <- .grf_evaluate_subgroup(sg_def, data)`; the effect re-selection `forestsearch_helpers.R:1612`; MR's family `fs_mr_inference_methods.R:29`) and by standalone `grf.subg.harm.glm(grf_selection = "frontier")`. The `use_grf = TRUE` screening path calls `grf.subg.harm.*` without `grf_selection` (`forestsearch_main.R:2548–2595` pass `return_selected_cuts_only` only, so `.build_grf_glm_args()`'s default `"tree"` applies) and consumes `grf_res$tree.cuts` (`:2601`, `grf_cuts <- grf_res$tree.cuts`) as `"var <= value"` strings for the consistency engine, whose `evaluate_comparison()` coerces all-numeric-level factors before comparing (`forestsearch_helpers.R:200–205`); on the tree path membership comes from `predict(trees[[d]], X)` on the coerced matrix (`grf_subg_harm_glm.R:585`). Of the three factor-carrying documents, `compare_all` and `validation_glm_simulation_study` use GRF only as the screening cut generator; `count_data_demo` additionally calls `grf.subg.harm.glm()` standalone with its default `grf_selection = "tree"` (`:381–384` pass no `grf_selection`). None of the applied documents runs `subgroup_method = "grf"` or the frontier on factor covariates: the multimethod documents that do run the GRF identifier coerce their binary covariates to numeric first (`fixed_family:407`).

## 4. Direct count on the survival replicate — the deciding fact

Two fits (the task's allowance), both `forestsearch()` with the template's `base_args` (`:1039–1056`: `is.RCT = TRUE`, `seedit = seed_base + sim_id`, `sg_focus = "effMaxSG"`, `subgroup_method = "grf"`, `hr.threshold = 0.9`, `hr.consistency = 0.8`, `pconsistency.threshold = 0.90`, `n.min = NULL`, `selection_rule = "neighborhood"`, `effect_neighborhood = 0.20`, `stop_threshold = NULL`) and the `grf` `method_args` (`:1067–1068`: `grf_selection = "frontier"`, `grf_depth = 2L`, `dmin.grf = 0`, `grf_select_statistic = "effect"`), `mr_inference = FALSE`, warnings captured with `withCallingHandlers()`; the bundle `meta` agrees (`subgroup_method grf`, `sg_focus effMaxSG`, `effect_neighborhood 0.2`, `ci_method field`, `field_scale_complement selected`; no `dmin.grf` key — a template literal).

- **Fit 1** (default RNG kind): `n_true` 56, 772 candidates, 0 with NA membership, 116 admitted, selected `{meno <= 0} & {pgr > 61.8}`, n 89, naive Cox HR 1.316244 — **not** the committed sim_id 1. Cause: the campaign's replicates ran inside `doFuture` workers (`.options.future = list(seed = TRUE)`), whose ambient generator is L'Ecuyer-CMRG, and `simulate_from_dgm()` seeds with `kind = NULL`, inheriting it; a plain session uses Mersenne-Twister (finding F2).
- **Fit 2** (`RNGkind("L'Ecuyer-CMRG")` before the draw, then the same lines): **`n_true` 72** (the bundle's row 1: 72); **warnings: none**; **candidates enumerated 779; with NA membership on the data (direct `.grf_evaluate_subgroup()` over every candidate) 0 of 779; `sel_effect` finite on 779 of 779; admitted 115**; **selected `{age > 45} & {meno <= 0}`, n = 100, naive Cox HR 1.643501** (no NA in the membership vector).
- **Committed `grfmr` sim_id 1** (`results/grf_effMaxSG_fb_mr_field_m1_h150_knoise0_n500_nb20_grfmr_combined_1_2000.rds`, present): `sg_def` `{age > 45} & {meno <= 0}`, `n_sel` 100, `nv_H_est` 1.643501, `admitted_n` 115, `n_family` 779, `n_true` 72. **They agree** on definition, n, naive effect, admitted count and family size.

The temporary directory `/tmp/fs_gfx_R5js` was removed at closeout.

## 5. DINA's evaluator — source read

DINA never compares a data-frame column: candidate membership is computed on a numeric matrix. `dina_subgroup()` builds it with `X <- .dina_extract_X_from_df(df, covariates)` (`R/dina_subgroup.R:388`), which **stops** on any non-numeric covariate — `is_num <- vapply(cov_df, is.numeric, logical(1L)); if (!all(is_num)) { stop("Covariate column(s) must be numeric: ", …) }; X <- as.matrix(cov_df); storage.mode(X) <- "double"` (`:652–658`). Candidates are evaluated on that matrix: depth 1 `x_j <- X[, j]; qs <- sort(unique(x_j)); … mask <- if (dir == "left") x_j <= q else x_j >= q` (`:714–721`), depth 2 `mask1 <- if (dir1 == "left") x1 <= q1 else x1 >= q1; … mask2 <- …; mask <- mask1 & mask2` (`:855–860`), and the selected candidate's `mask` is returned (`:593`). Before DINA is reached, every entry point coerces factors: `.forestsearch_dina_select()` — `df <- .coerce_covariates_numeric(df, confounders.name)`, likewise `df.predict`, `df.test` (`R/forestsearch_helpers.R:1396–1398`); the `use_dina` screening path — `df.dina <- .coerce_covariates_numeric(df.analysis, confounders.name)` (`R/forestsearch_main.R:2672`); `dina_subgroup_bootstrap()` (`R/dina_subgroup_bootstrap.R:309`). `.coerce_covariates_numeric()` (`helpers:1327–1345`): a factor/character column with all-numeric levels becomes `as.numeric(as.character(col))` ("matching `.build_grf_X()` and `evaluate_comparison()`"), and a column with non-numeric levels is collected as an offender and raises a `stop()` naming it. **DINA's evaluator does not have the defect:** factors are coerced (all-numeric levels) or rejected with an error (other levels) before any comparison; there is no path on which a factor is compared with `<=` / `>=`.

## The answer

The survival simulation does not carry factor covariates: the seven candidates `grfmr` searched (`er`, `age`, `meno`, `pgr`, `nodes`, `size`, `grade`) are integer columns of the simulated frame, and nothing converts them before the call. On the `grfmr` replicate A124_h150_n500 sim_id 1, regenerated with the campaign's seed scheme and RNG kind, the GRF identifier enumerated 779 candidates, 0 of which had NA membership, admitted 115, raised no warning, and selected `{age > 45} & {meno <= 0}` (n 100, naive HR 1.643501) — identical to the committed row, so that replicate dropped no candidate. Of the applied documents that run GRF or DINA, the four ACTG175 binary multimethod documents, the binary focus-sweep document and the GBSG survival multimethod document pass numeric covariates only (the ACTG175 ones coerce their binary indicators with `as.numeric()` for this reason); three documents pass factor covariates — `analysis_actg175_continuous_compare_all.qmd` (six binary factors), `count_data_demo.qmd` and `validation_glm_simulation_study.qmd` (simulated binary factors) — and all three reach GRF through the `use_grf` screening path or the standalone tree selection, not through the frontier evaluator under test, whose cuts are consumed by the consistency engine's own factor-coercing comparison. DINA's evaluator does not share the defect: every DINA entry coerces all-numeric-level factors to numeric and stops on any other factor before candidates are evaluated on a numeric matrix.

## Findings

- **F1.** The survival `grfmr` results are not exposed to the §1.6(c) defect: no factor covariate reaches the evaluator, and the direct count on sim_id 1 is 0 NA memberships of 779 candidates, with the committed selection reproduced exactly.
- **F2.** Reproducing a `grfmr` replicate outside the template needs `RNGkind("L'Ecuyer-CMRG")`: the template does not pin it (the MD template does, `sim_fs_maxeffCons_mr_field_md_template.qmd:624`), so its draws are reproducible only under the RNG kind of a `doFuture` worker. Fit 1, without the pin, gave a different sample (`n_true` 56) and a different selection.
- **F3.** The ACTG175 binary and GBSG applied documents coerce or carry numeric covariates and are not exposed; the three factor-carrying applied documents use GRF as a screening cut generator or the standalone tree, not the frontier evaluator; no applied document runs `subgroup_method = "grf"` on factor covariates.
- **F4.** DINA coerces (all-numeric levels) or rejects (other levels) factor covariates before evaluation; the MD design's binary factors would have been coerced on the DINA path.
- **F5.** The task's §2 names "setup, DGM and data chunks"; the survival template has no data chunk (the draw sits inside `record_replicate()` in `run-loop`), so those lines were run verbatim instead.
- **F6.** The `grfmr` bundle `meta` carries no `dmin.grf`, `grf_depth`, `grf_selection` or `selection_rule` key (template literals, resolved per batch by `gate3.R`); the fit used the template's literals.

## Commits

```
e3e6f652 Add TASK_grf_factor_exposure_2026-09-16 as received
<this record; then status_curated.md; then current_status.md alone>
```
