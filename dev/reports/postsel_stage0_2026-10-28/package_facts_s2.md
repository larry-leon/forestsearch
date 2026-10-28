# Package facts for Table S2 — forestsearch 0.4.0.9000 at HEAD (source read 2026-10-28)

All paths are forestsearch-repo-relative; line numbers are at HEAD `51fc55cc` (R/ identical to `936b6734`).

## 1. `forestsearch_bootstrap_dofuture()` (`R/bootstrap_dofuture_main.R`)

- **Formals** (lines 211–219):
  `fs.est, nb_boots, seed = 8316951L, details = FALSE, show_three = FALSE, parallel_args = list(), digits = 4, mr_in_replicates = FALSE`.
  The Table S2 document passes no `seed`, so its bootstrap used `8316951L`.
- **What it needs from the fit** (lines 230–250): `fs.est` must hold `df.est` (a non-empty data frame) and `args_call_all`. It must also have an identified subgroup: `sg.harm` non-empty and a `treat.recommend` column in `df.est`. Otherwise it stops ("requires an identified subgroup").
- **DINA / GRF fits are accepted.** `bootstrap_results()` (`R/bootstrap_analysis_dofuture.R:384–395`) reads `fs.est$args_call_all$subgroup_method`. For `"dina"` / `"grf"` it drops only `treat.recommend` and keeps the raw confounder columns, because each replicate re-fits `dina()` or the GRF forest from those columns. For `"consistency"` it also drops `confounders.candidate`. The roxygen section "Valid for every identifier configuration" (`bootstrap_dofuture_main.R`, @section after line 127) says this explicitly.
- **What each replicate re-runs** (`bootstrap_analysis_dofuture.R:548–614`). Each replicate clones `args_call_all` and calls `do.call(forestsearch, args_FS_boot)` (line 614). That repeats the whole search plus selection, under these overrides:
  - `details = show3`, `quiet = !show3`, and `plot.sg` / `plot.grf` / `show_candidate_summary = FALSE`.
  - `grf_res`, `grf_cuts`, `dina_res`, `dina_cuts` and `ps_hat` are set to `NULL`, so the front ends re-fit.
  - The inner plan is `"sequential"` with 1 worker.
  - **MR is stripped**: `if (!isTRUE(mr_in_replicates)) args_FS_boot$mr_inference <- FALSE` (line 607).
- **How the de-biased estimate is formed:**
  - Per replicate (`bootstrap_analysis_dofuture.R:724–729`): `H_biasadj_2 <- 2 * H_obs - (H_star + Hstar_star - Hstar_obs)`, on the log-HR scale. The method-1 form is `H_obs - (Hstar_star - Hstar_obs)`. The complement analogue is at 768–769.
  - Any non-finite constituent makes that replicate `NA`.
  - Across replicates (`bootstrap_dofuture_main.R:574–616`), `get_dfRes(Hobs, seHobs, H1_adj, H2_adj, ystar = Ystar_mat, cov_method = "standard")` takes the mean of the `H_biasadj_2` values. Its interval uses the infinitesimal-jackknife SE with the finite-B correction `varhat - termc`, where `termc = N/B.eval^2 * B.eval * mean(xc^2)` (`R/bootstrap_calculations_helpers.R:360–375`).
  - If `varhat <= termc` (possible at very small B), it falls back to the **uncorrected** variance.
- **Return** (lines 769–791): an `fs_bootstrap` object, a list with these elements:
  - `results` (data.table, one row per replicate; holds `H_biasadj_1/2`, `Hc_biasadj_1/2`, event counts, `tmins_search`, `tmins_iteration`, `Pcons`, `any_found`, ...)
  - `SG_CIs` (lines 688–693: `H_raw`, `H_bc`, `Hc_raw`, `Hc_bc`, as formatted strings)
  - `FSsg_tab`
  - `Ystar_mat`
  - `H_estimates`, `Hc_estimates` (data frames carrying `H0`, `H0_lower`, `H0_upper`, `H2`, `H2_lower`, `H2_upper`, ...)
  - `nb_boots`, `original_sg`, `outcome_type`, `effect_measure`, `est.scale`
  - `mr_replicates` (`NULL` under the default).

## 2. `mr_estimates_table()` and `summarize_bootstrap_results()`

- **`mr_estimates_table(fs, boot_results, est.scale = "hr", digits = 2)`** (`R/mr_estimates_table.R:41–42`) returns a **gt** table (via `format_bootstrap_table()`, line 102).
  - It takes the full-bootstrap table `boot_results$FSsg_tab` as a template, which supplies the descriptive columns and the **naive** effect.
  - It then **replaces** the bias-corrected column (`HR*`) with MR's values:
    - Ĥ: `g$debiased$est/lower/upper` (lines 77–78)
    - Ĥᶜ: `g$complement$debiased$...` (lines 83–88)
  - So one call gives four of a Table S2 row-pair's six cells: naive Ĥ and Ĥᶜ, plus MR Ĥ and Ĥᶜ.
  - **It does not show the full-bootstrap corrected cells.** Those come from `summarize_bootstrap_results()$table`, a separate gt table with the same layout. Neither function produces Table S2's three-column (naive / full bootstrap / multiplier) layout. The manuscript built that from the payload's `$table`, whose values come from document-local code (`table_s2_calls/table_assembly.R`, chunk `export-section6-2-payload`). That code reads `g$naive`, `g$debiased`, `g$complement$naive/debiased`, `boot$H_estimates$H2/H2_lower/H2_upper` and `boot$Hc_estimates$...`.
  - The roxygen says "under the default `ci_method = "ij"`" (`mr_estimates_table.R:31`). That is stale: the default is `"field"` (`R/forestsearch_main.R:3854`). The interval it renders is the IJ one under both settings (see §3), so the stale wording does not change the output.
- **`summarize_bootstrap_results(sgharm, boot_results, create_plots = FALSE, est.scale = "hr", digits = 2)`** (`R/summarize_bootstrap_results.R:62–66`) returns a list with these elements:
  - `table` (gt, from `format_bootstrap_table`)
  - `diagnostics`
  - `diagnostics_table_gt`
  - `plots` (ggplot2, only if `create_plots = TRUE`; lines 301–313)
  - `timing`
  - `subgroup_summary`
- **Element paths in the bootstrap object:**
  - Bootstrap-corrected estimate and interval: `boot$H_estimates$H2`, `$H2_lower`, `$H2_upper` (Ĥ) and `boot$Hc_estimates$H2`, `$H2_lower`, `$H2_upper` (Ĥᶜ). These are the effect scale under `est.scale = "hr"`.
  - Naive values from the same object: `$H0`, `$H0_lower`, `$H0_upper`.
  - Formatted strings: `boot$SG_CIs$H_bc` and `$Hc_bc`.

## 3. MR intervals (`R/fs_mr_inference.R`)

- **Paths:**
  - `g$debiased$est`, `$lower`, `$upper`: effect scale, two-sided 95% (lines 1069–1071; `lower = to_eff(beta_deb - z975 * se)`).
  - `g$debiased$se_ij` and `$se_wald`.
  - Complement: `g$complement$debiased$est`, `$lower`, `$upper` (lines 858–862, built with `sec_used`).
  - Naive: `g$naive` and `g$complement$naive`, using the robust SE (lines 1067–1068).
- **IJ under both `"ij"` and `"field"`:**
  - `se <- if (ci_method == "wald") se_wald else se_ij_rep$se` (line 770).
  - `sec_used <- if (ci_method %in% c("ij", "field")) se_ijc_rep$se else sec` (line 854).
  - The roxygen agrees: "`"field"` ... the `debiased` element, its IJ SE and its IJ intervals are identical, and are returned in every call regardless of this argument" (lines 301–304). The field stream draws under `seed + 900000L` *after* the main multiplier stream (lines 322–324), so `"ij"` output is unaffected byte-for-byte.
  - Under `"field"` the extra `g$field` block adds `lower_1s`, `lower_2s`/`upper_2s`, `est2`, `lower_se`/`upper_se` and `complement$upper_1s_s` (roxygen lines 513–523).
  - The M-FS-ij rung checks all of this empirically.
- **Default change since the Table S2 render:**
  - `forestsearch()` now resolves `ci_method` to `"field"` (`R/forestsearch_main.R:3846–3854`, from TASK_cimethod_note_2026-09-09).
  - At 0.2.0 the document's own prose says the default was `"ij"` (qmd line 993: "By default the de-biased CI uses the same infinitesimal-jackknife variance as the full bootstrap (the `ci_method = 'ij'` default)").
  - That `"field"` default applies to the **consistency** branch only. The DINA and GRF branches go through `.fs_apply_mr()`, which keeps its own `"ij"` default: `ci_method = .g(mr_inference_args$ci_method, "ij")` (`R/fs_mr_inference_methods.R:171`). The roxygen at lines 137–139 states this is deliberate: "`ci_method` keeps this wrapper's own `"ij"` default, which is deliberately *not* the consistency branch's `"field"`; changing it is a separate decision."
  - The Table S2 calls pass no `ci_method`. At HEAD the FS call therefore runs the field block, and the DINA and GRF calls do not. The ladder confirms this: M-FS resolved `"field"`, M-DINA and M-GRF resolved `"ij"`.
  - `include_complement` resolves to `TRUE` (line 3844), `draws` to the call's 5000, and `ij_residual` to `"two_term"` (the `fs_mr_inference` default, line 610). "Two-term" is therefore the default estimate.

## 4. Rendering the bootstrap summaries in a quarto vignette

- **gt is the only table backend** these helpers return (`format_bootstrap_table`, `R/bootstrap_summaries_helpers.R:30` onward). It uses these gt functions:
  - `gt::gt`, `tab_header`, `cols_label`, `cols_align`, `tab_spanner`, `tab_footnote` (line 204), `tab_source_note` (line 268), `tab_options` (line 272), `tab_style`, `cell_fill`, `cell_text`
  - `gt::md` for the `n1` label (line 105).
- gt prints as inline HTML with a scoped `<style>` block, and needs no external assets. It renders inside an `embed-resources: true` / self-contained HTML vignette, as the 0.2.0 document's own render did (`summaries$table` and `mr_estimates_table(...)` printed in the committed html at `43b051b6`).
- `knitr::kable` cannot take these objects (they are `gt_tbl`). A kable rendering would have to be rebuilt from the numeric paths in §2 and §3. That is how the manuscript fragment does it (`kableExtra` over the payload's `$table`). kableExtra is in neither Imports nor Suggests of DESCRIPTION (knitr is in Suggests, line 66).
- With `create_plots = TRUE`, `summarize_bootstrap_results()` returns ggplot objects, which render in HTML. The DINA and GRF summaries in the document used `create_plots = FALSE`.
- **Nothing found that would fail in a self-contained HTML vignette.** One thing to watch: the gt footnote carries the subgroup definition string with `<=`. gt escapes it in HTML, but `gt::md()` content would be interpreted as markdown.
