# REPORT — Two identifier fixes: GRF membership coding (P1) and DINA proposal-floor orientation (P2)

Date: 2026-09-16 (UTC 2026-09-17). Machine: `pop-os` (64 physical cores, 251 GB; R 4.6.1, reference BLAS). Branch `feature/glm-extension`. Task: `dev/tasks/TASK_grf_dina_fixes_2026-09-16.md` (committed as received, `0a697abb`). Evidence: `REPORT_md_grf_stage1_2026-09-16.md` §1.6(c), `REPORT_md_dina_grf_stage0_2026-09-16.md` S0.2–S0.3, `quarto/simulations/gbsg_020/REPORT_grf_factor_exposure_2026-09-16.md`. Every R process ran with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1` and a sequential future plan; every replicate was drawn under L'Ecuyer-CMRG with the source campaign's seed (`seed_base + sim_id`), and every DGM under R's default generator, as the template renders build it.

**Outcome.** **P1 landed** (`0cd33f7b`). **P2 landed** (`064fce91`): its correctness gate passed, and the condition held because the §2.1 list kept no call for review (F7 is empty). NEWS `323de083`. Final install `Built: R 4.6.1; ; 2026-09-17 04:47:31 UTC; unix`.

## 1. Provenance and baseline — GATE PASS

```
pop-os
feature/glm-extension
4b27f6f1
4b27f6f1 gbsg_020 closeout: regenerate current_status.md at 48a2cca2
48a2cca2 gbsg_020 status_curated.md: open-work line -- GRF factor-membership exposure check: ...
35816257 GRF factor-membership exposure on the survival path, read-only check (TASK_grf_factor_exposure_2026-09-16): ...
[tracked modifications: none]   [pending under R/ DESCRIPTION NAMESPACE tests: none]
exposure record in HEAD
R/ unchanged since the install
[R / Rscript / quarto / deno processes, by process name: none]
R 4.6.1; ; 2026-09-16 05:57:14 UTC; unix
```
First commit: `0a697abb Add TASK_grf_dina_fixes_2026-09-16 as received`.

**Baseline tests** (`devtools::test()` on HEAD `0a697abb`, `timeout 90m`; no fallback needed), 409 s:
```
[ FAIL 0 | WARN 32 | SKIP 3 | PASS 5051 ]
SUMMARY: files 41 | tests 342 | expectations passed 5051 | failed 0 | skipped 3 | warnings 32 | errors 0
TEST BLOCKS: passing 339 | failing 0 | skipped 3
BASELINE RC=0 wall_s=412
```
Failing tests: none. Skipped (3): `test-cv-no-subgroup-edges.R` "multisession and sequential tenfold agree on the no-subgroup path" (dev-load artifact), and two in `test-cross-outcome-parity.R` ("primary fit did not identify").

## 2. Stage A — exposure list and before-captures

### 2.1 Committed DINA runs

`git grep -n -E 'use_dina *= *TRUE|subgroup_method *= *"dina"' -- quarto vignettes inst tests dev ':!dev/tasks'` returns 377 lines; restricted to code (R scripts, and R chunks of `.qmd`/`.Rmd`, comments stripped) it is **95 calls in 61 files**. The `adverse_outcome` default: `forestsearch(adverse_outcome = NULL)` (`R/forestsearch_main.R:1287`), resolved by `if (is.null(adverse_outcome)) adverse_outcome <- (outcome_type %in% c("binary", "count"))` (`:1342–1343`, again `:1718–1719`): TRUE for binary and count, FALSE for continuous; DINA's survival path passes `TRUE` (`:2249–2250`). The `outcome_type` default is `"survival"` (`:1283`).

| files (calls) | outcome | `adverse_outcome` | kind |
|---|---|---|---|
| `quarto/applications/actg175/analysis_actg175_binary_multimethod_{fixed_family,frontend,psi_v2_2,psi_v3a}.qmd` (2 each: `use_dina = TRUE` `:1043`/`:1036`/`:1021`/`:1021`; `subgroup_method = "dina"` `:1236`/`:1229`/`:1214`/`:1214`) | binary | `TRUE`, explicit (7 assignments per file) | document |
| `quarto/applications/actg175/_archive/20260730_analysis_actg175_binary_multimethod_psi_v2_2{A,A_mac,_mac,_mac_w2}.qmd` (2 each) | binary | `TRUE`, explicit | archived document |
| `dev/identifier-alignment/sim_analyses/analysis_actg175_binary_multimethod.qmd` (`:1051`, `:1245`) | binary | `TRUE`, explicit | dev document |
| `quarto/applications/gbsg/analysis_gbsg_survival_multimethod.qmd` (`:1289`, `:1515`) and the six `gbsg/_archive/2026-0*` documents (2 each) | survival | n/a (DINA survival path passes `TRUE`) | document / archived |
| `dev/identifier-alignment/sim_analyses/analysis_gbsg_survival_multimethod.qmd`, `dev/review/analysis_gbsg_survival_multimethod.qmd`, `dev/replication-check/legacy_v2_2A_reconstructed.qmd`, `dev/replication-check/v2_2new_rendered_source_prerename.qmd` (2 each) | survival (`gbsg`) | n/a | dev document |
| `quarto/dina/dina_vs_forestsearch_signature_comparison{,_2factor}.qmd`, `quarto/dina/subgroup_method_signature_recovery{,_2factor}.qmd` (1 each, `fs_use_dina = TRUE`) | survival (`y_sim`/`event_sim`, `family = "cox"`) | n/a | document |
| `quarto/dina/method_equivalence_checks.qmd:107` | binary | `NULL` → `TRUE`, and `TRUE` | document |
| `quarto/dina/smoke-tests/dina_screening_smoketest.qmd` (7), `smoke_test_dina_factor.qmd` (4), `quarto/resampling/smoke_test_dina_select_statistic.R:61` | survival (explicit) | n/a | smoke document / script |
| `quarto/simulations/actg175/binary_methods/*.qmd` (19 documents, 1–3 calls each, e.g. `actg175_binary_methods_template.qmd:229`, `:307`, `:312`), `_fsparams/fs3{,b,c,d,e}_params.R:79`, `dina_grf_routing_diagnostic.R:75`, `:93`, `:185` | binary | unset → `TRUE` (three documents set `TRUE`) | simulation documents and scripts; **no committed bundles** (`git ls-files quarto/simulations/actg175/binary_methods` lists 0 `.rds`) |
| `quarto/simulations/gbsg/dina_failure_diagnostic.R:30`, `dev/terminology-work/R/11_spec_probe2.R:53` | survival | n/a | script |

**Tests** (code, not results): `tests/testthat/test-cv-no-subgroup-dina-grf.R:33`, `test-no-subgroup-bootstrap.R:32`, `test-no-subgroup-dina-grf-plot.R:37`, `test-no-subgroup-dina-grf-summary.R:32`, `test-no-subgroup-unified-contract.R:24`, `test-sg-focus-transparency.R:130`, `:178` — survival, and binary via `.fs_args_for("binary", …)`, which sets no `adverse_outcome` (`tests/testthat/helper-synthetic-dgm.R:133–140`), so `TRUE`.

**Kept for §2.3: none.** No document or template call runs DINA on a binary or continuous outcome with `adverse_outcome = FALSE`. Because the grep cannot see an identifier passed through a variable, every tracked file under `quarto/`, `vignettes/`, `inst/`, `dev/` that mentions DINA and sets `outcome_type = "continuous"` or `adverse_outcome = FALSE` was also read (97 files): the binary coverage sweeps under `quarto/simulations/actg175/binary*/` iterate `methods <- c("consistency", "dina", "grf")` but analyze with `adverse_outcome <- TRUE` (`binary/mr_coverage_sweep_or15.qmd:70`, `:148`; the `FALSE` at `:221`, `:251` is the DGM calibration); `quarto/simulations/actg175/continuous/mr_coverage_sweep_md_harm.qmd` sets `adverse_outcome <- FALSE` (`:123`) but its `forestsearch()` call (from `:309`) does not pass `subgroup_method`, which it only records in `meta` (`:418`); the MD template's `FS_MD_METHOD=dina` path has no committed bundle; the rest set `subgroup_method = "consistency"` or no identifier (finding F2).

### 2.2 P1 call sites

`.grf_evaluate_subgroup()`, before the edit: the GLM frontier selection (`R/grf_subg_harm_glm.R:548`), the survival frontier selection (`R/grf_main.R:308`), the effect re-selection (`R/forestsearch_helpers.R:1612`), the `subgroup_method = "grf"` wrapper's membership tables (`:1850`), MR's family membership (`R/fs_mr_inference_methods.R:29`), and the β(Ĥ) structured membership path (`R/betaHhat_truth.R:88`). The forest-matrix coding it now shares: `.build_grf_X()` (`R/grf_subg_harm_glm.R:866`), called once, by `grf.subg.harm.glm()` (`:399`). **None lies on the screening-cut path or the standalone tree selection** — confirmed: `use_grf = TRUE` calls `grf.subg.harm.*` without `grf_selection` (`R/forestsearch_main.R:2548–2595`; `.build_grf_glm_args()` defaults to `"tree"`) and consumes `grf_res$tree.cuts` (`:2601`); the tree path's membership is `predict(trees[[d]], X)` (`R/grf_subg_harm_glm.R:585`). The survival forest's own matrix is built separately, `X <- apply(data[, confounders.name, drop = FALSE], 2, as.numeric)` (`R/grf_main.R:226`), and was not changed (finding F4).

### 2.3 Before-captures

Harness (temporary, outside the repo): `knitr::purl()` of the survival template `quarto/simulations/gbsg_020/sim_fs_maxeffCons_fb_mr_field_m1_template.qmd` (chunks `setup`, `build-dgm`) and of the MD template `quarto/simulations/actg175/continuous/sim_fs_maxeffCons_mr_field_md_template.qmd` at `894da993`–HEAD (chunks `setup-knobs`, `build-dgm`, `machinery`), each template's **own** `base_args <- list(` … `method_args <- switch(` … `)` block read verbatim from the file (m1 `:1039–1070`; MD `:685–720`) and dispatched with `do.call(forestsearch, c(base_args, method_args))`; warnings and messages captured with `withCallingHandlers`; the returned object saved as RDS.

- **F1** — `grfmr` cell A124_h150_n500, sim_id 1–3, knobs of `scripts_dinamr/grfmr.sh:19–21` (`FS_S7_METHOD=grf FS_S7_FOCUS=effMaxSG FS_S7_NBHD=0.20 FS_S7_FIELD_COMPLEMENT=TRUE FS_S7_FIELD_SCALEC=selected FS_S7_FIELD_DECOMP=TRUE FS_S7_FIELD_RECOV=TRUE FS_S7_IJ_RESIDUAL=two_term FS_S7_FB=none`), `FS_S7_HR=1.50 FS_S7_N=500`, `FS_S7_Z1Q` unset; MR on (template `mr_inference_on` TRUE).
- **F2** — `dinamr` cell A_h150_n500 (`scripts_dinamr/blockA.cells`: `- 500 1.50 A_h150_n500`), sim_id 1–3, knobs of `scripts_dinamr/campaign.sh:9–11` (the F1 set with `FS_S7_METHOD=dina`); MR on.
- **F3** — FS on the MD design, md40 n500, sim_id 1–3, `mdsgnb20`'s knobs (`FS_MD_FOCUS=effMaxSG FS_MD_NBHD=0.20 FS_MD_FIELD_SCALEC=selected FS_MD_CI=field FS_MD_FB=none`); MR on.
- **F4** — GRF on the MD design, md40 n500, sim_id 1, the MD template's GRF block from `894da993` (`grf_selection <- "frontier"`, `grf_select_statistic <- "effect"`, `grf_depth <- 2L`, `dmin.grf <- .env_num("FS_MD_DMIN_GRF", 30)`), `FS_MD_METHOD=grf`; MR off.
- **F5** — DINA on the MD design, md40 n500, sim_id 1, `adverse_outcome = FALSE` (the template's), S0.3's settings (`FS_MD_METHOD=dina`, `dina_args = list()`, `dina_select_statistic = "effect"`, `effMaxSG`, ε 0.20, `neighborhood`), `details = TRUE`; MR off.
- **F6** — as F5 with `y_sim` negated and `adverse_outcome = TRUE`; MR off.
- **F7** — none (§2.1).

**Before-captures (installed build `2026-09-16 05:57:14 UTC`)**, pasted from the harness summary:
```
forestsearch 0.3.5 from /home/larryleon/R/x86_64-pc-linux-gnu-library/4.6/forestsearch | Built R 4.6.1; ; 2026-09-16 05:57:14 UTC; unix | R 4.6.1
F1 sim 1 (grf, MR TRUE): sg [{age > 45} & {meno <= 0}] n 100 | naive 1.643501 | warnings: none | 25.2 s
   committed sim 1: [{age > 45} & {meno <= 0}] n_sel 100 nv_H_est 1.643501
F1 sim 2 (grf, MR TRUE): sg [{er > 162}] n 98 | naive 1.502234 | warnings: none | 25.7 s
   committed sim 2: [{er > 162}] n_sel 98 nv_H_est 1.502234
F1 sim 3 (grf, MR TRUE): sg [{age <= 45} & {size <= 38.400000000000091}] n 85 | naive 1.185634 | warnings: none | 24.1 s
   committed sim 3: [{age <= 45} & {size <= 38.400000000000091}] n_sel 85 nv_H_est 1.185634
F2 sim 1 (dina, MR TRUE): sg [{age >= 47} & {meno <= 0}] n 84 | naive 1.554377 | warnings: none | 11.3 s
   committed sim 1: [{age >= 47} & {meno <= 0}] n_sel 84 nv_H_est 1.554377
F2 sim 2 (dina, MR TRUE): sg [{er >= 214}] n 78 | naive 1.851937 | warnings: none | 7.9 s
   committed sim 2: [{er >= 214}] n_sel 78 nv_H_est 1.851937
F2 sim 3 (dina, MR TRUE): sg [{meno <= 0} & {size >= 26}] n 115 | naive 0.953903 | warnings: none | 6.8 s
   committed sim 3: [{meno <= 0} & {size >= 26}] n_sel 115 nv_H_est 0.953903
F3 sim 1 (FS, MR on): sg [!{cd40 <= 415} & !{cd80 <= 948}] n 83 | naive 85.056112 | warnings: none | 34.3 s
F3 sim 2 (FS, MR on): sg [{preanti <= 0} & !{age <= 30}] n 115 | naive 80.827509 | warnings: none | 35.2 s
F3 sim 3 (FS, MR on): sg [{preanti <= 167} & !{age <= 32}] n 125 | naive 73.830528 | warnings: none | 34.5 s
   committed mdsgnb20 sim 1: [!{cd40 <= 415} & !{cd80 <= 948}] n_sel 83 nv_H_est 85.056112
   committed mdsgnb20 sim 2: [{preanti <= 0} & !{age <= 30}] n_sel 115 nv_H_est 80.827509
   committed mdsgnb20 sim 3: [{preanti <= 167} & !{age <= 32}] n_sel 125 nv_H_est 73.830528
F4 sim 1 (GRF, MR off): enumerated 1257 | NA membership 645 | sel_effect finite 612 | admitted 234 | factor-comparison warnings 730 | sg [{preanti <= 792.80000000000018} & {cd40 > 364}] n 156 | warnings: ‘<=’ not meaningful for factors [x552] | ‘>’ not meaningful for factors [x178]
   coded-column identity: helper absent in this build (factor covariates hemo,homo,drugs,race,gender,str2,symptom)
F5 sim 1 (DINA, MR off): sg [] n NA | warnings: none
   Harm floor:          m_diff = 30.0000
   Candidates searched:  8324
   Candidates qualifying (>= floor, >= n.min): 3
   SELECTED: none -- no candidate met the harm floor and size constraint.
F6 sim 1 (DINA, MR off): sg [{cd40 >= 400} & {cd80 >= 1040}] n 78 | warnings: none
   Harm floor:          m_diff = 30.0000
   Candidates searched:  8324
   Candidates qualifying (>= floor, >= n.min): 2690
   SELECTED: {cd40 >=  400} & {cd80 >= 1040}  (n = 78, mean tau-hat = 74.2003)
   proposed 2690 | tau_hat (as ranked) min 30.0064 max 111.3391 | all >= m_diff 30: TRUE | sel_effect finite 2690, >= 30: 1768 | admitted_n 1768 | cor(tau_hat, sel_effect) 0.744
F7: no kept call (the §2.1 list has no document or template DINA call on a binary or continuous outcome with adverse_outcome = FALSE)
DONE
```
F1, F2 and F3 reproduce the committed `grfmr`, `dinamr` and `mdsgnb20` rows (definition, n, naive estimate); F4 reproduces the GRF Stage 1 replicate (1,257 enumerated, 645 NA memberships, 234 admitted, `{preanti <= 792.8} & {cd40 > 364}`, n 156); F5 reproduces S0.3 (8,324 searched, 3 qualifying on the raw floor, no selection). A first capture pass was discarded: its DGMs were drawn after F1 had switched the process's generator (finding F1).

## 3. P1 — GATE PASS, landed `0cd33f7b`

**Diff** (the coding lines of `.build_grf_X()`, `R/grf_subg_harm_glm.R:885–896` at `0a697abb`, moved verbatim into `.grf_code_column()`):
```diff
diff --git a/R/grf_helpers.R b/R/grf_helpers.R
index 2d91f50a..3be1e607 100644
--- a/R/grf_helpers.R
+++ b/R/grf_helpers.R
@@ -695,3 +695,37 @@ validate_grf_data <- function(W, D, n.min) {
 
   return(TRUE)
 }
+
+
+#' Code one covariate column as the GRF forest matrix codes it
+#'
+#' The single definition of how a candidate covariate enters GRF on the GLM
+#' path.  Numeric columns pass through unchanged.  A factor or character
+#' column whose levels are all numeric strings (for example `"0"`/`"1"`)
+#' becomes those values via `as.numeric(as.character())`.  Any other factor or
+#' character column becomes integer codes via `as.integer(as.factor())`.
+#'
+#' `.build_grf_X()` applies it to build the forest's covariate matrix and
+#' `.grf_evaluate_subgroup()` applies it before comparing a column with a
+#' candidate cut, so a cut is always evaluated on the scale the forest split
+#' on.  Before this helper existed the evaluator compared the raw column, and a
+#' factor covariate returned `NA` membership for every cut on it.
+#'
+#' @param x A covariate column.
+#' @return `x` unchanged when it is neither a factor nor a character vector;
+#'   otherwise a numeric vector (all-numeric levels) or an integer vector
+#'   (other levels).
+#' @noRd
+.grf_code_column <- function(x) {
+  if (!is.factor(x) && !is.character(x)) {
+    return(x)
+  }
+  lvls <- if (is.factor(x)) levels(x) else unique(x)
+  if (!anyNA(suppressWarnings(as.numeric(lvls)))) {
+    # All-numeric levels: preserve the original values.
+    as.numeric(as.character(x))
+  } else {
+    # Non-numeric levels: integer codes.
+    as.integer(as.factor(x))
+  }
+}
diff --git a/R/grf_subg_harm_glm.R b/R/grf_subg_harm_glm.R
index 8b11b6f4..263a3c08 100644
--- a/R/grf_subg_harm_glm.R
+++ b/R/grf_subg_harm_glm.R
@@ -882,17 +882,12 @@ create_glm_row <- function(
   # with downstream cut evaluation, causing valid cuts to be
   # silently dropped (e.g., "z1 <= 1" is trivially TRUE on
   # the {0, 1} scale but meaningful on the {1, 2} scale).
+  #
+  # The per-column coding lives in .grf_code_column() (grf_helpers.R), which
+  # .grf_evaluate_subgroup() also calls, so a candidate cut is evaluated on
+  # the same scale the forest split on.
   for (v in names(X)) {
-    if (is.factor(X[[v]]) || is.character(X[[v]])) {
-      lvls <- if (is.factor(X[[v]])) levels(X[[v]]) else unique(X[[v]])
-      if (!anyNA(suppressWarnings(as.numeric(lvls)))) {
-        # All-numeric levels: preserve original values
-        X[[v]] <- as.numeric(as.character(X[[v]]))
-      } else {
-        # Non-numeric levels: integer codes
-        X[[v]] <- as.integer(as.factor(X[[v]]))
-      }
-    }
+    X[[v]] <- .grf_code_column(X[[v]])
   }
   as.matrix(X)
 }
diff --git a/R/grf_subgroup_labels.R b/R/grf_subgroup_labels.R
index 57691b82..a3bfd375 100644
--- a/R/grf_subgroup_labels.R
+++ b/R/grf_subgroup_labels.R
@@ -219,7 +219,10 @@
       v  <- cj$variable[r]; op <- cj$op[r]; val <- cj$value[r]
       if (!v %in% names(df))
         stop("Subgroup variable '", v, "' not found in data.", call. = FALSE)
-      x <- df[[v]]
+      # Code the column exactly as the forest's covariate matrix does
+      # (.grf_code_column(), shared with .build_grf_X()): a factor with levels
+      # "0"/"1" is compared as 0/1, not as a factor, which R cannot order.
+      x <- .grf_code_column(df[[v]])
       member <- switch(op,
                        "<=" = x <= val,
                        ">"  = x >  val,
```
**Tests** `tests/testthat/test-grf-membership-coding.R` (5 blocks, 18 expectations): a 0/1 factor and a numeric covariate, a cut on each and on their conjunction, membership equal to the expected logical vector with no `NA` and no warning; a numeric-only cut returning `c(0L, 0L, 0L, 1L, 1L, 1L)`; the helper's column equal to `.build_grf_X()`'s for factor, non-numeric factor, character and numeric columns; a cut on a non-numeric-level factor using the integer codes.

**Temporary install** `R CMD INSTALL --no-test-load --library=<tmp>/lib_p1 .` (`Built 2026-09-17 04:28:03 UTC`). **After-P1 captures:**
```
forestsearch 0.3.5 from /tmp/claude-1000/gdfx/tmp/lib_p1/forestsearch | Built R 4.6.1; ; 2026-09-17 04:28:03 UTC; unix | R 4.6.1
F1 sim 1 (grf, MR TRUE): sg [{age > 45} & {meno <= 0}] n 100 | naive 1.643501 | warnings: none | 25.1 s
   committed sim 1: [{age > 45} & {meno <= 0}] n_sel 100 nv_H_est 1.643501
F1 sim 2 (grf, MR TRUE): sg [{er > 162}] n 98 | naive 1.502234 | warnings: none | 25.5 s
   committed sim 2: [{er > 162}] n_sel 98 nv_H_est 1.502234
F1 sim 3 (grf, MR TRUE): sg [{age <= 45} & {size <= 38.400000000000091}] n 85 | naive 1.185634 | warnings: none | 23.9 s
   committed sim 3: [{age <= 45} & {size <= 38.400000000000091}] n_sel 85 nv_H_est 1.185634
F2 sim 1 (dina, MR TRUE): sg [{age >= 47} & {meno <= 0}] n 84 | naive 1.554377 | warnings: none | 10.9 s
   committed sim 1: [{age >= 47} & {meno <= 0}] n_sel 84 nv_H_est 1.554377
F2 sim 2 (dina, MR TRUE): sg [{er >= 214}] n 78 | naive 1.851937 | warnings: none | 7.8 s
   committed sim 2: [{er >= 214}] n_sel 78 nv_H_est 1.851937
F2 sim 3 (dina, MR TRUE): sg [{meno <= 0} & {size >= 26}] n 115 | naive 0.953903 | warnings: none | 6.8 s
   committed sim 3: [{meno <= 0} & {size >= 26}] n_sel 115 nv_H_est 0.953903
F3 sim 1 (FS, MR on): sg [!{cd40 <= 415} & !{cd80 <= 948}] n 83 | naive 85.056112 | warnings: none | 33.6 s
F3 sim 2 (FS, MR on): sg [{preanti <= 0} & !{age <= 30}] n 115 | naive 80.827509 | warnings: none | 34.9 s
F3 sim 3 (FS, MR on): sg [{preanti <= 167} & !{age <= 32}] n 125 | naive 73.830528 | warnings: none | 33.8 s
   committed mdsgnb20 sim 1: [!{cd40 <= 415} & !{cd80 <= 948}] n_sel 83 nv_H_est 85.056112
   committed mdsgnb20 sim 2: [{preanti <= 0} & !{age <= 30}] n_sel 115 nv_H_est 80.827509
   committed mdsgnb20 sim 3: [{preanti <= 167} & !{age <= 32}] n_sel 125 nv_H_est 73.830528
F4 sim 1 (GRF, MR off): enumerated 1257 | NA membership 0 | sel_effect finite 1257 | admitted 474 | factor-comparison warnings 0 | sg [{preanti <= 792.80000000000018} & {cd40 > 364}] n 156 | warnings: none
   coded-column identity (evaluator helper vs .build_grf_X vs forest X.orig), factor covariates hemo,homo,drugs,race,gender,str2,symptom: hemo=TRUE homo=TRUE drugs=TRUE race=TRUE gender=TRUE str2=TRUE symptom=TRUE
F5 sim 1 (DINA, MR off): sg [] n NA | warnings: none
   Harm floor:          m_diff = 30.0000
   Candidates searched:  8324
   Candidates qualifying (>= floor, >= n.min): 3
   SELECTED: none -- no candidate met the harm floor and size constraint.
F6 sim 1 (DINA, MR off): sg [{cd40 >= 400} & {cd80 >= 1040}] n 78 | warnings: none
   Harm floor:          m_diff = 30.0000
   Candidates searched:  8324
   Candidates qualifying (>= floor, >= n.min): 2690
   SELECTED: {cd40 >=  400} & {cd80 >= 1040}  (n = 78, mean tau-hat = 74.2003)
   proposed 2690 | tau_hat (as ranked) min 30.0064 max 111.3391 | all >= m_diff 30: TRUE | sel_effect finite 2690, >= 30: 1768 | admitted_n 1768 | cor(tau_hat, sel_effect) 0.744
F7: no kept call (the §2.1 list has no document or template DINA call on a binary or continuous outcome with adverse_outcome = FALSE)
DONE
```
**Identity, before vs after-P1** (`identical()` on the complete returned object after removing timing fields and neutralizing environments, closures and the `data.table` self-reference pointer; timing fields removed: `fit$minutes_all`, `fit$find.grps$time_search`, `fit$mr_inference$timing_seconds`, `fit$mr_inference$field$timing_seconds`, `fit$mr_inference$field$complement$timing_seconds`, `fit$mr_inference$field$recovery$timing_seconds`, each where present):
```
F1_sim1.rds  object identical: TRUE  | warnings identical: TRUE  | timing fields removed: fit$minutes_all, fit$mr_inference$field$complement$timing_seconds, fit$mr_inference$field$recovery$timing_seconds, fit$mr_inference$field$timing_seconds, fit$mr_inference$timing_seconds
F1_sim2.rds  object identical: TRUE  | warnings identical: TRUE  | timing fields removed: fit$minutes_all, fit$mr_inference$field$complement$timing_seconds, fit$mr_inference$field$recovery$timing_seconds, fit$mr_inference$field$timing_seconds, fit$mr_inference$timing_seconds
F1_sim3.rds  object identical: TRUE  | warnings identical: TRUE  | timing fields removed: fit$minutes_all, fit$mr_inference$field$complement$timing_seconds, fit$mr_inference$field$recovery$timing_seconds, fit$mr_inference$field$timing_seconds, fit$mr_inference$timing_seconds
F2_sim1.rds  object identical: TRUE  | warnings identical: TRUE  | timing fields removed: fit$minutes_all, fit$mr_inference$field$complement$timing_seconds, fit$mr_inference$field$recovery$timing_seconds, fit$mr_inference$field$timing_seconds, fit$mr_inference$timing_seconds
F2_sim2.rds  object identical: TRUE  | warnings identical: TRUE  | timing fields removed: fit$minutes_all, fit$mr_inference$field$complement$timing_seconds, fit$mr_inference$field$recovery$timing_seconds, fit$mr_inference$field$timing_seconds, fit$mr_inference$timing_seconds
F2_sim3.rds  object identical: TRUE  | warnings identical: TRUE  | timing fields removed: fit$minutes_all, fit$mr_inference$field$complement$timing_seconds, fit$mr_inference$field$recovery$timing_seconds, fit$mr_inference$field$timing_seconds, fit$mr_inference$timing_seconds
F3_sim1.rds  object identical: TRUE  | warnings identical: TRUE  | timing fields removed: fit$find.grps$time_search, fit$minutes_all, fit$mr_inference$field$complement$timing_seconds, fit$mr_inference$field$timing_seconds, fit$mr_inference$timing_seconds
F3_sim2.rds  object identical: TRUE  | warnings identical: TRUE  | timing fields removed: fit$find.grps$time_search, fit$minutes_all, fit$mr_inference$field$complement$timing_seconds, fit$mr_inference$field$timing_seconds, fit$mr_inference$timing_seconds
F3_sim3.rds  object identical: TRUE  | warnings identical: TRUE  | timing fields removed: fit$find.grps$time_search, fit$minutes_all, fit$mr_inference$field$complement$timing_seconds, fit$mr_inference$field$timing_seconds, fit$mr_inference$timing_seconds
F4_sim1.rds  object identical: FALSE | warnings identical: FALSE | timing fields removed: fit$minutes_all
   all.equal:
   Component “grp.consistency”: Component “out_sg”: Component “candidates”: Component “sel_effect”: 'is.NA' value mismatch: 0 in current 645 in target
   Component “grp.consistency”: Component “out_sg”: Component “admitted_n”: Mean relative difference: 1.025641
   Component “grf_res”: Component “candidates”: Component “sel_effect”: 'is.NA' value mismatch: 0 in current 645 in target
   Component “grf_res”: Component “admitted_n”: Mean relative difference: 1.025641
F5_sim1.rds  object identical: TRUE  | warnings identical: TRUE  | timing fields removed: fit$minutes_all
F6_sim1.rds  object identical: TRUE  | warnings identical: TRUE  | timing fields removed: fit$minutes_all
```
*GATE P1:* F1, F2, F3, F5, F6 identical (F7 empty) — **pass**; F4: 0 factor-comparison warnings and 0 NA-membership candidates — **pass**; the evaluator's coded column `identical()` to `.build_grf_X()`'s and to the forest's `X.orig` for all seven factor covariates (`hemo, homo, drugs, race, gender, str2, symptom`) — **pass**; new tests — **pass**.

| F4 (MD, sim_id 1, GRF, MR off) | before | after P1 |
|---|---|---|
| enumerated | 1,257 | 1,257 |
| NA membership | 645 | 0 |
| `sel_effect` finite | 612 | 1,257 |
| admitted (≥ 30 harm-oriented MD) | 234 | 474 |
| factor-comparison warnings | 730 (552 `‘<=’`, 178 `‘>’`) | 0 |
| selection | `{preanti <= 792.8} & {cd40 > 364}`, n 156 | the same |

The only differences in F4's object are `candidates$sel_effect` (645 NA → finite) and `admitted_n` (234 → 474), in `grf_res` and in `grp.consistency$out_sg`.

## 4. P2 — GATE PASS, condition met, landed `064fce91`

**Diff** (the orientation transplanted from the admission's effect estimator, `R/consistency_resample.R:255–261`: `if (!adverse_outcome) { if (outcome_type == "continuous") df[[outcome.name]] <- -df[[outcome.name]] else if (outcome_type == "binary") df[[outcome.name]] <- 1L - df[[outcome.name]] }` — both flips negate the link-scale treatment effect, identity and logit):
```diff
diff --git a/R/dina_subgroup.R b/R/dina_subgroup.R
index 6d18ad3b..8cbfa89f 100644
--- a/R/dina_subgroup.R
+++ b/R/dina_subgroup.R
@@ -227,6 +227,16 @@
 #'   `"maxSG"`, `"minSG"`, `"eff"`.
 #' @param alpha confidence level for the Wald interval on the
 #'   subgroup-mean tau-hat.  Default `0.05`.
+#' @param tau_sign `1` (default) or `-1`; the orientation applied to the
+#'   per-patient tau-hat before the `m_diff` floor, the `sg_focus` ranking and
+#'   the reported mean and interval.  `1` uses the fit's own scale.  `-1`
+#'   negates it, for a fit whose harm direction is the negative effect: the
+#'   floor then reads "harm of at least `m_diff`".  [forestsearch()] passes
+#'   `-1` for continuous and binary outcomes with `adverse_outcome = FALSE`,
+#'   the same outcomes its effect estimator flips (`-Y`, `1 - Y`), so DINA's
+#'   proposal floor and the admission floor point the same way.  With `-1`,
+#'   `mean_tau_hat`, `ci` and `candidates$tau_hat` are on the negated scale.
+#'   The recorded `call` carries `tau_sign` only when it is not `1`.
 #'
 #' @return An object of class `"dina_subgroup"`, a list with components:
 #'   \describe{
@@ -244,7 +254,8 @@
 #'     \item{depth}{integer; number of covariates in the selected
 #'       subgroup (`1` or `2`).}
 #'     \item{n_subgroup}{integer size of the chosen subgroup.}
-#'     \item{mean_tau_hat}{scalar subgroup-mean tau-hat.}
+#'     \item{mean_tau_hat}{scalar subgroup-mean tau-hat (in the
+#'       `tau_sign` orientation).}
 #'     \item{se_mean_tau_hat}{Wald standard error, computed as
 #'       `sqrt(a_S^T vcov(fit) a_S)` (CONDITIONAL on the chosen
 #'       subgroup -- not selection-adjusted).}
@@ -314,7 +325,8 @@ dina_subgroup <- function(fit, df, covariates,
                           sg_focus = "maxSG",
                           selection_rule = "neighborhood",
                           effect_neighborhood = 0.10,
-                          alpha = 0.05) {
+                          alpha = 0.05,
+                          tau_sign = 1) {
 
   if (!inherits(fit, "dina")) {
     stop("`fit` must be a DINA object (class \"dina\" or \"dina_bagged\").")
@@ -351,6 +363,14 @@ dina_subgroup <- function(fit, df, covariates,
       alpha <= 0 || alpha >= 1) {
     stop("`alpha` must be a single numeric in (0, 1).")
   }
+  if (length(tau_sign) != 1L || !is.numeric(tau_sign) ||
+      !tau_sign %in% c(-1, 1)) {
+    stop("`tau_sign` must be 1 or -1.")
+  }
+  # The default orientation is not recorded, so a call that leaves tau_sign
+  # at 1 reads exactly as it did before the argument existed.
+  call_rec <- match.call()
+  if (tau_sign == 1) call_rec$tau_sign <- NULL
 
   # Normalize the GLM-natural vocabulary ("eff", "effMaxSG", "effMinSG")
   # to the canonical internal form ("hr", "hrMaxSG", "hrMinSG") shared
@@ -408,8 +428,11 @@ dina_subgroup <- function(fit, df, covariates,
   }
   V <- stats::vcov(fit)
 
-  # Per-patient tau-hat
-  tau_hat <- as.numeric(beta[1L] + X %*% beta[-1L])
+  # Per-patient tau-hat, in the requested orientation (tau_sign = -1 negates
+  # it, so the floor, the ranking and the reported mean are harm-oriented for
+  # a fit whose harm direction is the negative effect).  The Wald variance
+  # below is unchanged by the sign.
+  tau_hat <- tau_sign * as.numeric(beta[1L] + X %*% beta[-1L])
 
   # ---- Collect candidates -------------------------------------------------
   # Depth-1 singletons at full resolution (all unique thresholds), shared
@@ -479,7 +502,7 @@ dina_subgroup <- function(fit, df, covariates,
       n_total                 = n,
       n_candidates_searched   = n_searched,
       n_candidates_qualifying = n_qualifying,
-      call                    = match.call()
+      call                    = call_rec
     )
     class(out) <- "dina_subgroup"
     return(out)
@@ -610,7 +633,7 @@ dina_subgroup <- function(fit, df, covariates,
       v2 = covariates[cand_j2], d2 = cand_dir2, c2 = cand_q2,
       tau_hat = cand_tau,   # native DINA ranking statistic (subgroup-mean tau-hat)
       stringsAsFactors = FALSE),
-    call                    = match.call()
+    call                    = call_rec
   )
   class(out) <- "dina_subgroup"
   out
diff --git a/R/forestsearch_helpers.R b/R/forestsearch_helpers.R
index f7611a1e..0cb1e6ba 100644
--- a/R/forestsearch_helpers.R
+++ b/R/forestsearch_helpers.R
@@ -1358,6 +1358,31 @@ reset_workers <- function(workers   = NULL,
 }
 
 
+#' Orientation of DINA's tau-hat for the proposal floor
+#'
+#' DINA is fit on the raw outcome, so its tau-hat is on the outcome's own
+#' scale.  The admission floor scores candidates with the effect estimator,
+#' which flips the outcome when `adverse_outcome = FALSE`: `-Y` for continuous
+#' and `1 - Y` for binary outcomes (`.consistency_glm_pieces()`,
+#' `R/consistency_resample.R`).  Both flips negate the link-scale treatment
+#' effect (identity and logit), so the proposal floor is applied to `-tau-hat`
+#' for exactly those outcomes.  Survival, count, and `adverse_outcome = TRUE`
+#' keep `1`.
+#'
+#' @param outcome_type Character outcome type.
+#' @param adverse_outcome Logical, as resolved by `forestsearch()`.
+#' @return `-1` or `1`.
+#' @noRd
+.dina_tau_sign <- function(outcome_type, adverse_outcome) {
+  if (!isTRUE(adverse_outcome) &&
+      outcome_type %in% c("continuous", "binary")) {
+    -1
+  } else {
+    1
+  }
+}
+
+
 #' DINA-selection mode for forestsearch (subgroup_method = "dina")
 #'
 #' Fits a DINA model and delegates subgroup selection to
@@ -1410,6 +1435,10 @@ reset_workers <- function(workers   = NULL,
   # (cox/binomial/poisson), identity (mean difference) for gaussian.
   m_diff <- if (identical(da$fit$family, "gaussian")) hr.threshold
             else log(hr.threshold)
+  # ...applied to tau-hat in the orientation the admission floor uses: DINA is
+  # fit on the raw outcome, while the admission scores candidates with the
+  # effect estimator, which flips the outcome for these outcome types.
+  tau_sign <- .dina_tau_sign(outcome_type, adverse_outcome)
 
   if (isTRUE(details)) {
     lines <- c(
@@ -1421,7 +1450,9 @@ reset_workers <- function(workers   = NULL,
       paste0("  effect_neighborhood: ", effect_neighborhood),
       paste0("  Harm floor:          ", sprintf("m_diff = %.4f", m_diff),
              if (!identical(da$fit$family, "gaussian"))
-               sprintf("  (hr.threshold = %.4g)", hr.threshold) else ""),
+               sprintf("  (hr.threshold = %.4g)", hr.threshold) else "",
+             if (tau_sign < 0)
+               "  on -tau-hat (harm-oriented; adverse_outcome = FALSE)" else ""),
       paste0("  n.min:               ", n.min)
     )
 
@@ -1461,7 +1492,8 @@ reset_workers <- function(workers   = NULL,
     grid_probs          = da$select$grid_probs,
     sg_focus            = sg_focus,
     selection_rule      = selection_rule,
-    effect_neighborhood = effect_neighborhood
+    effect_neighborhood = effect_neighborhood,
+    tau_sign            = tau_sign
   )
 
   # Effect-based re-selection (dina_args$select_statistic = "effect"): re-rank
```
Design notes: `dina_subgroup()` stores `match.call()` in its result, which `forestsearch()` embeds (`grp.consistency$out_sg$call`); the call therefore records `tau_sign` only when it is not `1`, so every existing call reads as before (the F2 identity depends on it). The `use_dina` screening path (`R/forestsearch_main.R:2694–2710`) does not pass `tau_sign` and is unchanged (finding F3).

**Tests** `tests/testthat/test-dina-proposal-orientation.R` (5 blocks, 25 expectations; constructed `"dina"` fits, no model fitted): `.dina_tau_sign()` for the seven (outcome, `adverse_outcome`) pairs; with a negative raw surface, `tau_sign = 1` qualifies nothing and `tau_sign = -1` finds a set whose every candidate has oriented tau-hat ≥ `m_diff` and raw mean ≤ −`m_diff`; with a positive surface, `tau_sign = 1` explicit is `identical()` to the default (call included); `tau_sign = 2` is an error; `.forestsearch_dina_select()` finds a subgroup under `adverse_outcome = FALSE` and none under `TRUE` on the negative surface. No existing test changed expectation (§5).

**Temporary install** `<tmp>/lib_p2` (`Built 2026-09-17 04:35:36 UTC`). **After-P2 captures:**
```
forestsearch 0.3.5 from /tmp/claude-1000/gdfx/tmp/lib_p2/forestsearch | Built R 4.6.1; ; 2026-09-17 04:35:36 UTC; unix | R 4.6.1
F1 sim 1 (grf, MR TRUE): sg [{age > 45} & {meno <= 0}] n 100 | naive 1.643501 | warnings: none | 23.7 s
   committed sim 1: [{age > 45} & {meno <= 0}] n_sel 100 nv_H_est 1.643501
F1 sim 2 (grf, MR TRUE): sg [{er > 162}] n 98 | naive 1.502234 | warnings: none | 23.9 s
   committed sim 2: [{er > 162}] n_sel 98 nv_H_est 1.502234
F1 sim 3 (grf, MR TRUE): sg [{age <= 45} & {size <= 38.400000000000091}] n 85 | naive 1.185634 | warnings: none | 23.7 s
   committed sim 3: [{age <= 45} & {size <= 38.400000000000091}] n_sel 85 nv_H_est 1.185634
F2 sim 1 (dina, MR TRUE): sg [{age >= 47} & {meno <= 0}] n 84 | naive 1.554377 | warnings: none | 10.3 s
   committed sim 1: [{age >= 47} & {meno <= 0}] n_sel 84 nv_H_est 1.554377
F2 sim 2 (dina, MR TRUE): sg [{er >= 214}] n 78 | naive 1.851937 | warnings: none | 7.4 s
   committed sim 2: [{er >= 214}] n_sel 78 nv_H_est 1.851937
F2 sim 3 (dina, MR TRUE): sg [{meno <= 0} & {size >= 26}] n 115 | naive 0.953903 | warnings: none | 6.4 s
   committed sim 3: [{meno <= 0} & {size >= 26}] n_sel 115 nv_H_est 0.953903
F3 sim 1 (FS, MR on): sg [!{cd40 <= 415} & !{cd80 <= 948}] n 83 | naive 85.056112 | warnings: none | 31.6 s
F3 sim 2 (FS, MR on): sg [{preanti <= 0} & !{age <= 30}] n 115 | naive 80.827509 | warnings: none | 32.9 s
F3 sim 3 (FS, MR on): sg [{preanti <= 167} & !{age <= 32}] n 125 | naive 73.830528 | warnings: none | 31.8 s
   committed mdsgnb20 sim 1: [!{cd40 <= 415} & !{cd80 <= 948}] n_sel 83 nv_H_est 85.056112
   committed mdsgnb20 sim 2: [{preanti <= 0} & !{age <= 30}] n_sel 115 nv_H_est 80.827509
   committed mdsgnb20 sim 3: [{preanti <= 167} & !{age <= 32}] n_sel 125 nv_H_est 73.830528
F4 sim 1 (GRF, MR off): enumerated 1257 | NA membership 0 | sel_effect finite 1257 | admitted 474 | factor-comparison warnings 0 | sg [{preanti <= 792.80000000000018} & {cd40 > 364}] n 156 | warnings: none
   coded-column identity (evaluator helper vs .build_grf_X vs forest X.orig), factor covariates hemo,homo,drugs,race,gender,str2,symptom: hemo=TRUE homo=TRUE drugs=TRUE race=TRUE gender=TRUE str2=TRUE symptom=TRUE
F5 sim 1 (DINA, MR off): sg [{cd40 >= 400} & {cd80 >= 1040}] n 78 | warnings: none
   Harm floor:          m_diff = 30.0000  on -tau-hat (harm-oriented; adverse_outcome = FALSE)
   Candidates searched:  8324
   Candidates qualifying (>= floor, >= n.min): 2690
   SELECTED: {cd40 >=  400} & {cd80 >= 1040}  (n = 78, mean tau-hat = 74.2003)
   proposed 2690 | tau_hat (as ranked) min 30.0064 max 111.3391 | all >= m_diff 30: TRUE | sel_effect finite 2690, >= 30: 1768 | admitted_n 1768 | cor(tau_hat, sel_effect) 0.744
   direction check: proposed with harm-oriented MD (sel_effect) > 0: 2675 of 2690 | tau_sign recorded in call: tau_sign
F6 sim 1 (DINA, MR off): sg [{cd40 >= 400} & {cd80 >= 1040}] n 78 | warnings: none
   Harm floor:          m_diff = 30.0000
   Candidates searched:  8324
   Candidates qualifying (>= floor, >= n.min): 2690
   SELECTED: {cd40 >=  400} & {cd80 >= 1040}  (n = 78, mean tau-hat = 74.2003)
   proposed 2690 | tau_hat (as ranked) min 30.0064 max 111.3391 | all >= m_diff 30: TRUE | sel_effect finite 2690, >= 30: 1768 | admitted_n 1768 | cor(tau_hat, sel_effect) 0.744
   direction check: proposed with harm-oriented MD (sel_effect) > 0: 2675 of 2690 | tau_sign recorded in call: absent (1)
F7: no kept call (the §2.1 list has no document or template DINA call on a binary or continuous outcome with adverse_outcome = FALSE)
DONE
```
**Identity:**
```
before vs after-P2:
F1_sim1.rds  object identical: TRUE  | warnings identical: T
F1_sim2.rds  object identical: TRUE  | warnings identical: T
F1_sim3.rds  object identical: TRUE  | warnings identical: T
F2_sim1.rds  object identical: TRUE  | warnings identical: T
F2_sim2.rds  object identical: TRUE  | warnings identical: T
F2_sim3.rds  object identical: TRUE  | warnings identical: T
F3_sim1.rds  object identical: TRUE  | warnings identical: T
F3_sim2.rds  object identical: TRUE  | warnings identical: T
F3_sim3.rds  object identical: TRUE  | warnings identical: T
F6_sim1.rds  object identical: TRUE  | warnings identical: T
after-P1 vs after-P2:
F4_sim1.rds  object identical: TRUE  | warnings identical: T
```
(F5 before vs after-P2 differs by design: `grp.consistency` and `df.est` were `NULL` before, with no subgroup found.)

*GATE P2 — correctness:* F1, F2, F3, F6 identical to before and F4 identical to after-P1 — **pass**; F5 proposes a non-empty set (2,690 candidates), every one with oriented tau-hat ≥ 30 (min 30.0064) — **pass**; the proposal and admission floors point the same way — **pass** (2,675 of the 2,690 proposed have a positive harm-oriented MD, 1,768 clear the admission floor 30, correlation of oriented tau-hat with the harm-oriented MD 0.744; before, the 3 proposed were benefit candidates and 0 were admitted); F5's selection `{cd40 >= 400} & {cd80 >= 1040}`, n 78, equals F6's — **pass** (the whole candidate summary coincides).

| F5 (MD, sim_id 1, DINA, MR off) | before | after P2 | F6 (negated outcome, `adverse_outcome = TRUE`) |
|---|---|---|---|
| proposal floor as applied | `m_diff = 30` on raw tau-hat (a benefit of 30) | `m_diff = 30` on −tau-hat (harm) | `m_diff = 30` on raw tau-hat of −Y (harm) |
| admission floor | 30 on the harm-oriented MD | the same | the same |
| searched / proposed | 8,324 / 3 | 8,324 / 2,690 | 8,324 / 2,690 |
| admitted | 0 | 1,768 | 1,768 |
| selection | none | `{cd40 >= 400} & {cd80 >= 1040}`, n 78 | the same |

**Condition (Larry's):** every F7 fit identical to its before-capture — the F7 list is empty (§2.1), so the condition holds and **P2 landed**. No patch file was written.

## 5. Tests, news, final install — GATE PASS

`devtools::test()` on the P1 + P2 tree (`timeout 90m`), 409 s:
```
[ FAIL 0 | WARN 32 | SKIP 3 | PASS 5094 ]
SUMMARY: files 43 | tests 352 | expectations passed 5094 | failed 0 | skipped 3 | warnings 32 | errors 0
TEST BLOCKS: passing 349 | failing 0 | skipped 3
FINAL RC=0 wall_s=409
```
Block by block: all **339** blocks that passed at baseline pass; the **10** new blocks pass; the same 3 are skipped; warnings 32 at both runs. No failure at baseline, none now.

`NEWS.md`: two bullets under "forestsearch (development version)", `323de083`. `devtools::document()` rewrote only `man/dina_subgroup.Rd` (the `tau_sign` parameter and the `mean_tau_hat` wording), committed with P2; its warning "Failed to parse the inline R code: `r = se_mean / sd_emp`" concerns another roxygen block and predates this task.

**Final install** `devtools::install(quick = TRUE, upgrade = FALSE)`: `Built: R 4.6.1; ; 2026-09-17 04:47:31 UTC; unix`; two doFuture multisession workers (pids 1786951, 1786952) report the same `Built` from the main library. F4 and F5 re-run against the main library:
```
F4_sim1.rds  object identical: TRUE  | warnings identical: T
F5_sim1.rds  object identical: TRUE  | warnings identical: T
```
The temporary directory and both temporary libraries were removed after this record's numbers were taken.

## 6. Findings

- **F1 (harness).** The first before-capture pass drew F2–F6's DGMs after F1's replicate had pinned the process's generator to L'Ecuyer-CMRG, so F2–F5 did not reproduce their campaigns. The harness now resets to R's default generator before each template's DGM chunks, as a template render builds them; that pass was discarded and all captures were redone. With the reset, F1–F5 reproduce the committed or recorded replicates.
- **F2 (§2.1).** The task's grep sees only literal `use_dina = TRUE` / `subgroup_method = "dina"`; several documents choose the identifier through a variable (`fs_subgroup_method`, `methods <- c(...)`, `FS_S7_METHOD`, `FS_MD_METHOD`). A second pass over DINA-mentioning files with a continuous outcome or `adverse_outcome = FALSE` found none exposed.
- **F3 (P2 scope).** The `use_dina` screening path (`R/forestsearch_main.R:2694–2710`) applies the same unoriented floor (`m_diff_sel`) when it contributes DINA's selected cut to the consistency search; it was not changed. Under `adverse_outcome = FALSE` on a continuous or binary outcome it still contributes a benefit-side cut.
- **F4 (P1 scope).** The survival GRF forest matrix is `apply(data[, confounders.name], 2, as.numeric)` (`R/grf_main.R:226`), not `.build_grf_X()`. On all-numeric data (every committed survival campaign) the evaluator is unchanged; with factor covariates the two codings agree on numeric-level factors and differ on others (the survival builder gives `NA`), and `apply()` on a mixed frame goes through a character matrix. Not changed.
- **F5.** With `details = TRUE`, the DINA wrapper's frontier display (`dina_frontier()`) still shows raw tau-hat; only the "Harm floor" line states the orientation.
- **F6.** On the MD design, P1 doubles GRF's admitted set on sim_id 1 (234 → 474) without changing that replicate's selection; the MR family on the GRF path will now include the binary-covariate candidates.
- **F7.** The P2 condition was evaluated on an empty list; it holds vacuously.
- **F8.** For a non-numeric-level factor or character column, `.grf_code_column()` codes by that frame's own levels, so a prediction frame with a different level set could code differently — a property `.build_grf_X()` already had, now shared by the evaluator.
- **F9.** `find.grps$time_search` (`R/subgroup_search.R:986`, `t.sofar`) is an elapsed time; it was added to the removed-field list after the first comparison flagged it.
- **F10.** Both new test files first failed on a test-side error (the matrix's double storage mode; `-1` recorded in the call as an expression); the tests were corrected, the code was not.
- **F11.** The P2 and NEWS commits were made after the full suite's aggregate result (0 failures) and before the block-by-block comparison, whose first attempt did not run (argument error); run again, it confirmed all 339 baseline blocks pass.
- **F12.** The first baseline launch was stopped (missing output-path variable) and restarted; the reported baseline is the restarted run.
- **F13.** F2 (`dinamr`) reproduces its committed rows on this build, so the survival DINA results are reproducible from the template with the harness recipe.

## 7. Commits

```
0a697abb Add TASK_grf_dina_fixes_2026-09-16 as received
0cd33f7b GRF membership coding (TASK_grf_dina_fixes_2026-09-16 P1): ...
064fce91 DINA proposal-floor orientation (TASK_grf_dina_fixes_2026-09-16 P2): ...
323de083 NEWS: GRF membership coding on factor covariates (P1) and DINA proposal-floor orientation (P2), TASK_grf_dina_fixes_2026-09-16
<this record; then status_curated.md; then current_status.md alone>
```
