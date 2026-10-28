# Package facts — vignettes stage 0 (2026-10-27)

Source: forestsearch at HEAD of `feature/vignettes-0.4.1` (R/ tree identical to master `174778d0`); run-time facts from the installed build (0.4.0.9000, built 2026-10-27 21:01:28 UTC, `NEWS.md` byte-identical to the repo's). Line numbers are into `R/` at that HEAD.

## 2.1 The consistency-evaluation method (GATE — passed)

**Formal.** `forestsearch(consistency_method = c("resample", "split"))` — `R/forestsearch_main.R:1686`; resolved by `match.arg` at `:1718`, so the default is **`"resample"`**. The same formal exists on `subgroup.consistency()` (`R/subgroup_consistency_main.R:407`, default `"resample"`). The value is passed through at `R/forestsearch_main.R:3573`.

**`"resample"` is the closed form.** Roxygen, `R/forestsearch_main.R:1029-1045`:

> `"resample"` (default) uses the multiplier (influence-function / `dfbeta`) approximation (`consistency_resample`), returning the rate from a single subgroup fit; `"split"` is retained as the fallback and runs the literal repeated 50/50 split-and-refit calculation, which costs `fs.splits * 2` refits per candidate evaluation. Despite the name, `"resample"` performs no numerical resampling at evaluation time -- the rate is a closed-form expression. […] configurations it cannot represent that way (IRD, propensity-adjusted effects, or a non-convergent fit on a given candidate) fall back to `"split"` automatically. When `use_twostage = TRUE`, `"resample"` bypasses the two-stage split screening entirely (the rate is computed directly).

The closed form itself, `R/consistency_resample.R:455-457`:

```r
  if (method %in% c("closed", "both")) {
    out$rate_closed <- max(0, 2 * stats::pnorm(delta / pieces$sigma_D) - 1)
  }
```
with `delta = beta_hat - log(c2)` and `sigma_D = sqrt(sum_i dfbeta[i, treat]^2)` (file header, `:8-23`). Both consistency call sites pass `method = "closed"`.

**Fallback conditions (quoted).**

1. Whole-run fallback, `R/subgroup_consistency_main.R:427-438` — only if `consistency_resample()` is missing from the session, or on the GLM `estimator_fn` path without a `glm_resample_spec`:
   ```r
   if (identical(consistency_method, "resample")) {
     if (!exists("consistency_resample", mode = "function")) {
       warning("consistency_resample() not found; using 'split'. ...")
       consistency_method <- "split"
     } else if (!is.null(estimator_fn) && is.null(glm_resample_spec)) {
       warning("consistency_method = 'resample' on the GLM (estimator_fn) path ", "requires a glm_resample_spec; using 'split'.")
       consistency_method <- "split"
   ```
   Neither applies to the survival path from `forestsearch()`.
2. Per-candidate fallback, fixed-loop evaluator, `R/subgroup_consistency_helpers.R:1485-1500`:
   ```r
   p.consistency <- round(rr$rate_closed, pconsistency.digits)
   if (is.na(p.consistency)) {
     # Inestimable closed form (e.g. non-convergent fit on a candidate):
     # fall back to literal splitting rather than dropping the subgroup.
     p.consistency <- .consistency_via_splits(df.x, N.x, n.splits, ...)
   ```
3. Per-candidate fallback, two-stage evaluator (the one the Table 1 FS fit uses: `consistency_algorithm = "twostage"`), `R/subgroup_consistency_helpers.R:1863-1907`: `cox_resample <- identical(consistency_method, "resample") && is.null(estimator_fn)`; then `resample_unsupported <- is.na(p.consistency)` (`:1890`) — if TRUE the candidate "falls through to the Stage 1/2 splitting below rather than dropping the candidate" (`:1887-1889`).

**Does `fs.splits` matter when the closed form is evaluable?** No. `fs.splits` reaches the engine only as `n.splits` (`R/forestsearch_main.R:3557`), which is used only by `.consistency_via_splits()` (`R/subgroup_consistency_helpers.R:355`) and the Stage 1/2 split loop — both reached only on fallback. RNG: "Under consistency_method = "resample" no RNG is consumed and every path gives identical results" (`R/subgroup_consistency_main.R:714-715`). So the FS search is deterministic given the data; the Monte Carlo in Table 1 is the MR multiplier draws and the field's outer/inner draws.

**What the result records.** Method: `fit$args_call_all$consistency_method` (`"resample"` in every FS rung). Algorithm: `fit$consistency_algorithm` / `fit$grp.consistency$algorithm` (`"twostage"` for FS; `"dina"` / `"grf"` otherwise). **No per-candidate closed-form/fallback flag and no fallback count is recorded** anywhere in the result (`fit$grp.consistency` names: `out_sg, sg_focus, df_flag, sg.harm, sg.harm.id, algorithm, n_candidates_evaluated, n_candidates_total, n_passed, early_stop_triggered, early_stop_candidate, stop_threshold, seed`). The ladder measured it by tracing (`ladder.R`): **S0–S3: 120 closed-form evaluations, 0 NA, 0 fallbacks**. `fs_declaration_calibration()` reports the method as `consistency_method` with `consistency_method_source = "fit$args_call_all$consistency_method"`.

**Note on the transplanted FS call.** It passes `consistency_method = "resample"` and `fs.splits = 1000` explicitly; by the above, `fs.splits` was inert for Table 1. DINA and GRF have no consistency screen: `.fs_admission_applies()` returns `consistency = FALSE` for `dina`/`grf` (`R/forestsearch_helpers.R:2416-2421`); D1/G1 traced 0 consistency evaluations.

## 2.2 Formals

- `formals(forestsearch::forestsearch)`: 86 formals with deparsed defaults in `formals_forestsearch.txt`. Defaults relevant to the vignette: `parallel_args = list(plan = "multisession", workers = .default_parallel_workers(), show_message = TRUE)`, `seedit = 8316951`, `sg_focus = "hr"`, `hr.consistency = 1`, `pconsistency.threshold = 0.9`, `fs.splits = 1000`, `consistency_method = c("resample", "split")`, `subgroup_method = c("consistency", "dina", "grf")`, `mr_inference = FALSE`, `mr_inference_args = list()`, `max.minutes = 3`, `max_n_confounders = 1000`, `vi.grf.min = NULL`, `quiet = FALSE`, `details = FALSE`, `plot.sg = FALSE`, `plot.grf = FALSE`.
- MR routine: `fs_mr_inference` (exported). `names(formals())`: `df, candidates, spec, selected_members, admission, t_confirm, confirm_rule, reselection, effect_neighborhood, selection_rule, draws (2000L), multiplier, include_complement, ci_method, seed, return_reselection, field_R_out (1000L), field_R_in (500L), field_uniform, field_M_cap, field_complement, field_decompose, field_scale_complement, ij_residual, field_recovery, keep_declaration_field, keep_field_matrix, declaration_c0, pconsistency.digits` (with defaults in `formals_forestsearch.txt`).
- Documented elements of `mr_inference_args` (`R/forestsearch_main.R:956-1028`): `confirm_rule`, `t_confirm`, `reselection`, `draws` (default 2000), `multiplier`, `ci_method`, `include_complement`, `seed` (defaults to `seedit`), `return_reselection`, `field_M_cap`, `field_complement`, `ij_residual`, `keep_declaration_field`, `keep_field_matrix`, `declaration_c0`.
- **Field outer/inner draws are not settable on the FS path.** The consistency-engine MR call (`R/forestsearch_main.R:3828-3896`) forwards none of `field_R_out` / `field_R_in` (nor `multiplier`-independent field sizes), so `fs_mr_inference()`'s defaults 1000 × 500 always apply for `subgroup_method = "consistency"`. The DINA/GRF path goes through `.fs_apply_mr()` (`R/fs_mr_inference_methods.R:141`), which **does** forward `field_R_out` / `field_R_in` from `mr_inference_args` (`:178-179`). Not documented in the `mr_inference_args` roxygen. (OPEN ITEM in the report.)

## 2.3 MR products and where they live (from rung S1, `fits/fit_S1.rds`)

`g <- fit$mr_inference`; `cc <- g$complement`; `fb <- g$field`; `fc <- fb$complement`; `js <- fb$joint_s`.

| Table 1 quantity | element path | S1 value |
|---|---|---|
| region | `fit$sg.harm` (pasted with `" & "`); MR's own label `g$selected_label` is in factor coding (`"q1.1 & q15.1"`) | `{er <= 0} & {pgr <= 26}` |
| N, events in Ĥ / Ĥᶜ | membership `fit$df.est$treat.recommend == 0` (Ĥ) matched on id; events summed from the data. `g$n_selected` (75) and `cc$n` (611) give N; **events are not stored** in the result | 75 (41) / 611 (258) |
| unadjusted HR, 95% CI, Ĥ | `g$naive$est`, `g$naive$lower`, `g$naive$upper` | 2.222 (1.169, 4.221) |
| unadjusted HR, 95% CI, Ĥᶜ | `cc$naive$est`, `$lower`, `$upper` | 0.605 (0.465, 0.788) |
| corrected (two-term), Ĥ | `g$debiased$est` | 1.277 |
| corrected (two-term), Ĥᶜ | `cc$debiased$est` | 0.652 |
| field lower bound on β(Ĥ) | `fb$lower_1s` | 0.634 |
| field-s upper bound on β(Ĥᶜ) | `fc$upper_1s_s` | 0.812 |
| Bonferroni pair | `js$bonf_lower_H`, `js$bonf_upper_Hc` (`js$bonf_joint_prob`, `js$bonf_gamma`) | 0.509 / 0.849 |
| MR family size | `g$n_family` | 1744 |
| MR draws used | `g$settings$draws`; also `g$debiased$ij_draws` | 500 |
| field outer / inner | `fb$R_out`, `fb$R_in` (used: `fb$n_out_used`, `fb$n_in_used_mean`; complement `fc$n_out_used`) | 1000 / 500 (990 / 499.09) |
| seed | **not stored** in `g$settings`; the call uses `mr_inference_args$seed` else `seedit` (`R/forestsearch_main.R:3896`); `fit$args_call_all$seedit` = 8316951; `fb$seed_offset` = 900000; `fit$grp.consistency$seed` | 8316951 |

**What prints.** `print.forestsearch` (`R/forestsearch_methods.R:244`): definition, N, HR (unadjusted), Pcons, field lower on Ĥ, field-s upper on Ĥᶜ, Bonferroni pair, two-sided IJ interval, p̂(H). `summary.forestsearch` (`:367`): the same plus settings, candidate counts, SEs (field, field-s, IJ), top-3 re-selection mass, and the family size. **Neither prints**: the corrected (two-term) estimates for Ĥ or Ĥᶜ, the unadjusted CI for Ĥ or the unadjusted HR/CI for Ĥᶜ, events in either region, the draw counts, or the seed. `fit$mr_inference` has no print method (S3 class absent; it prints as a raw list, ~3,000 lines). Both print `sg_focus: hrMaxSG` (the normalized form of `effMaxSG`). The nearest reporting helper is `mr_estimates_table(fs, boot_results, est.scale = "hr", ...)` (`R/mr_estimates_table.R:41`, exported): it renders `g$debiased$est` / `cc$debiased$est` with their IJ intervals (`:78`, `:85-87`) into the bootstrap-table format, so it **requires a `boot_results` object** and carries no field / field-s / Bonferroni bound. No exported helper assembles the nine Table 1 quantities (exports checked: `NAMESPACE` `fs_*` and `mr_*`).

**`declaration_c0`.** Reachable through `forestsearch(mr_inference_args = list(keep_declaration_field = TRUE, declaration_c0 = c(...)))` — the transplanted FS call passes both. It stores `g$declaration_field` (`Mstar`, `beta_hat`, `sigma_D`, `family_id`, `meta`, `Mstar_c0` [500 × 4]). `fs_declaration_calibration(fit, c0 = c(0.70, 0.75, 0.80, 0.85))` on S1 returns class `fs_declaration_calibration`: `kappa_hat` 3.77, implied screen `Pcons >= 0.895; T >= 1.621`, `fw_size` 0.976, admitted 9 (current) / 23 (p\* rule) / 0 (calibrated), and per-c0 rows (`kappa_hat`, `fw_size`, `pstar_set`, `z_eff`, `z_gap`, `digits_fine`, `pstar_fine`, `n_admitted`, `q95`): c0 0.70 → κ̂ 2.807, p\*-settable 1 at 2 digits; 0.75/0.80/0.85 → none at 2 digits (fine p\* 0.9969 / 0.9983 / 0.9991). These are the S1 (500-draw) values; the published diagnostic uses 5000 draws. Full output in `postprocess_log.txt`.

## 2.4 Build facts

**DESCRIPTION** (verbatim):

```
Package: forestsearch
Title: Exploratory Subgroup Identification in Survival and GLM Outcomes
Version: 0.4.0.9000
Authors@R: person("Larry", "Leon", , "larry.leon.05@post.harvard.edu", role = c("aut", "cre"))
Description: Implements statistical methods for exploratory subgroup 
    identification in clinical trials. Provides tools for identifying patient 
    subgroups with differential treatment effects using machine learning 
    approaches including Generalized Random Forests (GRF), LASSO 
    regularization, and exhaustive combinatorial search algorithms. Supports 
    survival endpoints (Cox proportional hazards), binary outcomes (log odds 
    ratio, log relative risk, risk difference), continuous outcomes (mean 
    difference), and count / rate outcomes (log incidence rate ratio via 
    Poisson, quasi-Poisson, or negative-binomial GLMs with optional 
    person-time offset). Features bootstrap bias correction using 
    infinitesimal jackknife methods to address selection bias in post-hoc 
    analyses. Designed for clinical researchers conducting exploratory 
    subgroup analyses in randomized controlled trials, particularly for 
    multi-regional clinical trials (MRCT) requiring regional consistency 
    evaluation. Methods are described in Leon et al. (2024) 
    <doi:10.1002/sim.10163>.
    Post-selection inference for the identified subgroup follows León and
    Anderson (2026) <doi:10.48550/arXiv.2609.38361>.
    The search methodology is given in León and Anderson (2026)
    <doi:10.48550/arXiv.2610.28379>.
License: MIT + file LICENSE
Depends: 
    R (>= 4.1.0)
Encoding: UTF-8
Roxygen: list(markdown = TRUE)
VignetteBuilder: knitr, quarto
SystemRequirements: Quarto command line tools (https://github.com/quarto-dev/quarto-cli)
Imports: 
    data.table,
    doFuture,
    dplyr,
    foreach,
    future,
    future.apply,
    future.callr,
    ggplot2,
    glmnet (>= 5.0),
    grf,
    gt,
    Matrix,
    mvtnorm,
    parallel,
    patchwork,
    policytree,
    progressr,
    randomForest,
    rlang,
    stats,
    survival,
    utils,
    weightedsurv (>= 0.1.0)
Suggests: 
    DiagrammeR,
    doRNG,
    htmltools,
    MASS,
    sandwich,
    tidyr,
    forestploter (>= 1.2.0),
    cubature,
    svglite,
    knitr,
    rmarkdown,
    katex,
    latentcor,
    quarto,
    speff2trial,
    testthat (>= 3.0.0)
URL: https://github.com/larry-leon/forestsearch, https://larry-leon.github.io/forestsearch/
BugReports: https://github.com/larry-leon/forestsearch/issues
Config/testthat/edition: 3
Config/roxygen2/version: 8.1.1
```

**Imports vs Suggests:** `grf` Imports; `policytree` Imports; `quarto` Suggests; `knitr` Suggests; `rmarkdown` Suggests; `speff2trial` Suggests.

**.Rbuildignore** (verbatim):

```
^vignettes/results$
\.quarto
# Development files and directories
^dev$
^.*\.Rproj$
^\.Rproj\.user$
# Git files
^\.git$
^\.gitignore$
# CI/CD
^\.github$
^\.travis\.yml$
# Documentation source
^README\.Rmd$
^cran-comments\.md$
# Other common exclusions
\.DS_Store
\.Rhistory
\.RData
^vignettes/articles$
^_pkgdown\.yml$
^docs$
^pkgdown$
^quarto$
# Knitr/vignette build artifacts
^vignettes/.*_cache$
^vignettes/.*_files$
# L_eff calibration data (consumed by quarto docs, not shipped)
^_data$
^\.positai$
^\.claude$
^_vignette_test$
^CLAUDE\.md$
^inst/patches$
^CRAN-SUBMISSION$
```

**_pkgdown.yml:** 595 lines; copied whole to `~/Downloads/vignettes_stage0_2026-10-27/_pkgdown.yml` rather than inlined here.

**vignettes/ (recursive, tracked files):**

| file | size | engine / `%\Vignette*` |
|---|---|---|
| `vignettes/forestsearch.Rmd` | 41,211 | `%\VignetteIndexEntry{Exploratory Subgroup Identification with ForestSearch}`, `%\VignetteEngine{knitr::rmarkdown}`, `%\VignetteEncoding{UTF-8}`; YAML: title "Exploratory Subgroup Identification with ForestSearch", subtitle "German Breast Cancer Study Group (GBSG) Analysis", `output: html_document` (toc 3, number_sections, code_folding hide, cosmo) |
| `vignettes/survival-post-selection.qmd` | 9,218 | `%\VignetteIndexEntry{Post-selection inference for an identified survival subgroup}`, `%\VignetteEngine{quarto::html}`, `%\VignetteEncoding{UTF-8}`; YAML: `format: html: embed-resources: true, toc, toc-depth 3, code-fold false` |
| `vignettes/.gitignore` | 28 | `/.quarto/`, `**/*.quarto_ipynb` |
| `vignettes/articles/methodology.Rmd` | 52,722 | `{ForestSearch Methodology}`, `knitr::rmarkdown`; bibliography `references.bib` |
| `vignettes/articles/paper_simulations.Rmd` | 106,952 | `{Simulation Studies for Evaluating ForestSearch Performance}`, `knitr::rmarkdown` |
| `vignettes/articles/biomarker_effects.Rmd` | 23,005 | `{Biomarker Treatment Effects in the Weibull AFT/Cox Framework}`, `knitr::rmarkdown` |
| `vignettes/articles/causal_effects_brief_review.Rmd` | 45,976 | `{Causal Interpretation of Hazard Ratios in Survival Analysis}`, `knitr::rmarkdown` |
| `vignettes/articles/treatment_effect_definitions.Rmd` | 38,425 | `{Treatment Effect Definitions in ForestSearch}`, `knitr::rmarkdown` |
| `vignettes/articles/extreme_subgroups.qmd` | 95,330 | **no `vignette:` block** (Quarto article; `execute: echo true, warning false`) |
| `vignettes/articles/references.bib` | 222,531 | — |
| `vignettes/articles/.gitignore` | 28 | as above |

Untracked/ignored on disk: `vignettes/survival-post-selection.html` (1,226,485 B, 2026-09-09; ignored by `vignettes/*.html`). `vignettes/articles` is in `.Rbuildignore`, so only `forestsearch.Rmd` and `survival-post-selection.qmd` are CRAN vignettes.

**tests referencing `vignettes/` or `inst/doc`:** none (`grep -rln 'vignettes/\|inst/doc' tests/testthat/` empty).
