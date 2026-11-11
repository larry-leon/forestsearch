# REPORT — NPV for a replicate with no subgroup: all remaining instances, 2026-11-11

Task: `dev/tasks/TASK_npv_nofind_instances_2026-11-11.md`, committed as received on `master` (`3b589565`).
Outcome: **all gates passed.** The branch `fix/npv-nofind-instances` is not pushed and not merged.

Definition applied throughout (León et al. 2024, Section 3): with no subgroup identified, the estimated complement is everyone, so sens = 0, ppv = 0 (0/0 convention), spec = sens(Hc) = 1 and npv = ppv(Hc) = |Hc|/N, which is 1 under the null.

## 1. Provenance (Step 0)

| Item | Value |
|---|---|
| Host | `pop-os`, R 4.6.2, testthat 3.3.2, devtools 2.5.2 |
| `master` at start | `0af900ab` (`dev: A-003 fix report`) |
| `git status --porcelain` on `master` at start | empty (**gate passed**) |
| A-003 gate: `grep -cE 'npv +<- +sum\(!true_H\) / nrow\(df\)' R/run_simulation_analysis.R` | 4 (**gate passed**; A-003 fix `e2f9ec18` is in `master`) |
| Task commit on `master` | `3b589565 dev: add NPV no-find instances task` (one file, copied from `~/Downloads`) |
| `DESCRIPTION` version at start | 0.4.1.9004 (as expected) |
| Branch | `fix/npv-nofind-instances` did not exist (**gate passed**); created from `3b589565` |
| Sibling repos at start | fs-post-selection `main` `5d165acc`, clean; fs-glms-interpretable `main` `396144cd`, clean |
| Fix commit / report commit | given in the closing message (the report cannot hold its own SHA) |

## 2. Part A — inventory

### 2.1 How the search was run

- **Package (`R/`, `man/`, `tests/`):** read by hand, every hit of `npv|NPV|Enpv|ppv_sgc|ppv(Hc)|ppv_Hc|npv_Hc|negative predictive` and `any.H == 0|detected == 0|NO-DETECTION|!det`. An independent read-only reviewer then repeated the search with no access to these conclusions (§7).
- **forestsearch `quarto/` and `vignettes/`:** 509 tracked source files (`.qmd .R .Rmd .md`) mention NPV. They fall into 113 clusters of files whose NPV-related lines are byte-identical (copies of one template). Eight read-only agents each read one representative per cluster and spot-checked a second member. An adversarial checker re-opened every flagged row. A completeness critic then searched with other spellings and fills (`coalesce`, `replace(is.na(npv), 1)`, prose).
- **fs-glms-interpretable:** started at `dev/payloads/R/identifiers_metrics.R:359-362`; included `manuscript/`, `quarto/`, `dev/`.
- **fs-post-selection:** `manuscript_rev1/` and its generators (per `manuscript/payload_manifest.qmd`, the generators are templates in forestsearch `quarto/simulations/gbsg_020/`).
- **Stored payloads:** every tracked `.rds` was read (read-only), and `fs_mr_oc_summary()` was recomputed on every tracked payload that stores `$oc$identification` (§2.4).
- Evidence is the current source. Line numbers for `R/` sites are on `master`, with the branch line in brackets where they differ.

### 2.2 Part A table — package sites (forestsearch `R/`, `man/`, `tests/`)

| Repo | File:line | Kind | NPV and spec for a no-subgroup replicate (on `master`) | Reachable under an alternative? | \|Hc\|/N available at the site? | Reported outputs it feeds |
|---|---|---|---|---|---|---|
| forestsearch | `R/fs_mr_oc_summary.R:48-57` (roxygen, `unconditional` item); `man/fs_mr_oc_summary.Rd:81-84` | roxygen | States npv = 1, spec = 1, "following `build_classification_table()` and Leon et al. (2024, Table 1)" | yes, it describes the code below | yes (see next row) | help page |
| forestsearch | `R/fs_mr_oc_summary.R:280-310` (`.fs_oc_identification()`, unconditional row: comment `:284-285`, `fill()` `:286-291`, `npv = .fs_oc_mean(fill("npv", 1))` `:309`) | code + comment | npv **1**; spec 1 (`fill("spec", 1)`) | **yes**: 792 of the 872 tracked MR payloads with no-detection rows are alternatives (mean \|Hc\|/N on those rows < 1; minimum 0.624) | **yes**: per-replicate `results$n_true` (true \|H\|) and `meta$n_sample` (N). Verified on 1,178 tracked payloads: the found-case NPV reconstructs from `n_true`, `n_harm`, `sens` and N = `meta$n_sample` to 1.1e-13. All 1,266 MR-shaped payloads carry `n_true`, never NA on a no-detection row | `$oc$identification` written into MR payloads by 5 batch templates (§2.3, row Q4); any caller's summary |
| forestsearch | `R/simulation_tables.R:360-372` (`build_classification_table()`: comment `:361-367`, fill `:369-372`) | code + comment | npv NA → **1** and spec NA → 1 when `any.H == 0`, alternative only. Comment says the producer stores NA for sens and ppv on a no-find row (false since `95e2f819`, and under A-003 it stores 0) | **yes, but only on rows stored with npv = NA**: output of `run_simulation_analysis()` from before `95e2f819` (2026-04-13), output with no truth column, or another producer. Current extractor output with a truth column carries the values already, so the fill does nothing there | **yes where the row carries it**: every `run_simulation_analysis()` row is `cbind(df_pop, …)` (`R/run_simulation_analysis.R:608, 632, 661, 694`) and `df_pop$propHc_true = 1 - mean(harm)` on the same simulated data (`:517-546`). **Not** when the truth column is absent (then `propHc_true` is NA too) | `ppv(Hc)` / `sens(Hc)` rows of the classification gt table (e.g. `vignettes/articles/paper_simulations.qmd:1156`, actg175 binary_methods documents) |
| forestsearch | `R/oc_analyses.R:683` (`format_oc_results()`) | comment | Comment: "averaged across sims with detection". Code (`:686-690`) is `mean(res$npv, na.rm = TRUE)` over all sims, so a no-find sim contributes its stored value (\|Hc\|/N, spec 1 since A-003) | the comment misstates it; the code is correct given the producer | yes (stored value) | NPV / Spec columns of `format_oc_results()` tables |
| forestsearch | `R/oc_analyses.R:449-456` (`summarize_single_analysis()`) | code | pass-through `mean(…, na.rm = TRUE)` of stored values | correct given the producer | yes (stored) | `summarize_simulation_results()` |
| forestsearch | `R/run_simulation_analysis.R:1014-1034, 1059-1073, 1383-1399, 1409-1425` (four no-find branches) | code | npv \|Hc\|/N, spec 1 (A-003); null: 1, 1 | correct | yes (`true_H`) | every downstream table |
| forestsearch | `tests/testthat/test-fs_mr_oc_summary.R:108` | test | pinned unconditional npv 0.85 = (0.7 + **1**)/2 | yes (alternative fixture) | (fixture had no `n_true`) | — |
| forestsearch | `tests/testthat/test-nofind-classification.R:48, 91` | test | 0.75 = \|Hc\|/N (alt); 1 (null) | correct | yes | — |

**Examined and not Part A rows** (they do not score the León NPV against the true subgroup for a replicate with no identified subgroup):

| File:line | Why not a Part A row |
|---|---|
| `R/fs_oc_predict.R:309-310`, `R/fs_oc_grid.R:556-557, 615-616` (`Enpv`) | Analytic expectation weighted by `sel_c`, the selection distribution **given declaration**: a detected-only quantity; no non-declared term exists |
| `R/fs_oc_family.R:77-78, 444` (roxygen + comment "NPV = 1 follows downstream"); `man/fs_oc_family_enumerate.Rd:119`; `tests/testthat/test-fs-oc-null.R:77` (`Enpv == 1`) | Per-candidate NPV of a **selected** rule under the null (P(Q) = 0), which is 1; not a no-subgroup replicate |
| `R/fs_identification_summary.R:149-162, 181-183` | Explicitly detected-only (`det <- r[r$detected %in% 1L, ]`); NA when nothing is detected |
| `R/fs_mr_inference.R:1159-1223` (`.fs_mr_field_recovery()`, `npv_Hc` alias) | Agreement of re-selected draws with the observed pick, not with the truth; draws with no winner are skipped and counted (`n_skipped`) |
| `R/forestsearch_cross_validation.R:1365-1376, 1409-1420` (`ppv_Hc`, `sens_Hc`); `R/cv_summary_tables.R` | Agreement of fold assignments with the full-data pick, not with the truth; the no-find value is NA, not 1. See open item 4 (the same branch inverts `sens_H`) |

### 2.3 Part A table — outside the package (listed, not changed)

Copies of one template share a row; "×n" is the number of tracked files.

| Repo | File:line | Kind | NPV and spec for a no-subgroup replicate | Reachable under an alternative? | \|Hc\|/N available at the site? | Reported outputs it feeds |
|---|---|---|---|---|---|---|
| forestsearch | Q1 `quarto/simulations/actg175/binary_methods/actg175_binary_methods_template.qmd:354-368` (`oc_agg_from_results()`, fill `:357-360`) ×1 | quarto source (code) | npv NA → **1**, spec NA → 1, no hypothesis check | only on NA rows (pre-`95e2f819` bundles or no truth column); the template's own `run_simulation_analysis()` output carries values | yes: `propHc_true` in the same rows | H1 (Spec, NPV) bubble plot and Pareto flag (`:751-759`), H0 bubble (`:807`) |
| forestsearch | Q2 same template `:119-126` ×1 (+ tracked `.html`) | quarto source (prose) | states npv = 1, spec = 1 for Ĥc "averaged unconditionally", H1 included | yes, as a statement | n/a (text) | rendered prose in the tracked `actg175_binary_methods_template.html` |
| forestsearch | Q3 `quarto/simulations/gbsg/reproduce_table1_m1.qmd:634-647` (fill `:641`) ×2 (also `quarto/simulations/reproduce_table1_m1.qmd`) | quarto source (code) | npv NA → **1** under alt; spec NA → 1; comment `:636-637` already states ppv(Hc) = \|Hc\|/N | only on NA rows; latent with the current package | yes: `ra` is `run_simulation_analysis()` output (`:332-411`), so `propHc_true` is in the row | in-document "M1 Alternative: Paper vs Reproduction" table (`ppv(Hc)`), discrepancy table; no tracked render |
| forestsearch | Q4 batch templates that store `fs_mr_oc_summary()`: `quarto/simulations/actg175/continuous/sim_fs_maxeffCons_mr_md40_knoise0_n500_batch_1_1000.qmd:875, :968` (×3 md40/md120 copies); `quarto/simulations/actg175/binary_020/sim_fs_mr_field_or_template.qmd:1220-1221, :1316` and the continuous field template (×2) | quarto source (calls package) | inherited from the package: **1** on `master`; \|Hc\|/N on this branch | yes | yes (`n_true`, `meta$n_sample` in the payload) | `$oc$identification` in 106 tracked payloads; 62 store the old value (§2.4) |
| forestsearch | Q5 `quarto/simulations/actg175/binary_legacy/_summaries/actg175_binary_m1_harm_fs1-fs6_summary.qmd:259-293` ×1 | quarto source (code) | averages stored pre-A-003 bundles: npv = **1**, spec = \|Hc\|/N for a no-find H1 replicate | yes | in the bundles (`propHc_true`), not used | OC gt table, Spec × NPV bubble, frontier table, 3D scatter (render not tracked) |
| forestsearch | Q6 17 tracked renders: `quarto/simulations/actg175/binary_methods/*.html` (16) and `binary_legacy/` (1); e.g. `actg175_binary_m1f_harm_n1000_full_adjusted_resample.html:6371` (`ppv(Hc)` 0.97 at any(H) = 0.49) | rendered output | displays ppv(Hc) and sens(Hc) from the pre-A-003 producer (npv **1**, spec \|Hc\|/N) | yes | n/a (render) | "Harm Subgroup Identification" tables (A-003 open item 8) |
| forestsearch | Q7 `dev/reports/articles_finish_2026-10-30/results/paper_simulations_output.txt:560` ×4 (also `paper_simulations.txt` and the `articles_alignment_2026-10-30` pair) | tracked text of a 2026-10-30 render | ppv(Hc) from the pre-A-003 producer (M Alt GRF: any(H) 0.61, ppv(Hc) 1.00) | yes | n/a | record of `vignettes/articles/paper_simulations.qmd` §6.2.2; the article source itself is correct at HEAD and would change on re-render |
| forestsearch | Q8 `quarto/simulations/actg175/binary_methods/actg175_binary_m1_harm_methods.qmd:1294-1308` ×105 `.qmd` (+ `.rmarkdown`/`.html` copies) | quarto source (prose) | gives sens = 0, ppv = 0 for a no-find H1 replicate and defers NPV to the package; says it matches `fs.estimates.out()` (the claim A-003 removed from the package) | through the package (correct at HEAD) | n/a | prose in those documents' renders |
| forestsearch | Q9 NA then `na.rm` (NPV/spec silently detected-only): `quarto/dina/grf_dmin_calibration.qmd:159-166, 282` ×1; `quarto/applications/validation_glm_simulation_study.qmd:560-591, 721-733` ×1; `quarto/archive/actg175_simulations_binary_working.qmd:472-552` ×1 (its `build_classification_table()` call `:728-746` then fills npv 1 under alt, since its own extractor has no `propHc_true`); `quarto/dina/dina_subgroup_signature_simulation.qmd:212-239, 405-416` ×1 (spec only) | quarto source (code) | npv NA, spec NA → dropped by `na.rm` | yes (detection < 1 in signal cells) | yes (truth vector in scope at each no-find return) | gt/kable OC tables in the tracked `.html` of each (archive render not tracked) |
| forestsearch | Q10 per-replicate NA, no summary in the template: `quarto/simulations/gbsg/mr_coverage_sweep.qmd:209-219, 236, 303-307` ×54 (+1 `.rmarkdown`); dead legacy extractor copies `quarto/legacy calibration files/run_simulation_analysis.R:435-571, 681-830` (also `dev/glm/…` ×2) | quarto source (code) | npv NA, spec NA on the stored row | stored only; dead copies are not sourced | yes (`rec$n_true`; `df[[harm_col]]`) | per-cell `*_res.rds` bundles (797 tracked) and coverage grids; nothing tracked for the legacy copies |
| forestsearch | Q11 reports quoting detected-only NPV without that label: `quarto/simulations/gbsg_020/REPORT_p30sg_2026-09-07.md:19, 79`; `REPORT_nb20_2026-09-07.md:77, 104`; `REPORT_dinamr_blockC_2026-09-11.md:359-366`; `REPORT_field_recovery_2026-09-09.md:407`; `p12x20_2026-09-12/REPORT_p12x20_A1..A9_2026-09-12.md:92` and `GATE2_*.txt` ×9 each; `REPORT_fs_extraction_2026-09-11.md:199-201` (calls non-detected classification "undefined") | report | detected-only means (non-detections excluded), unlabelled on the NPV line | yes (planted-H cells with detection < 1) | yes in the bundles | report text only |
| forestsearch | Q12 correct or not applicable (no change needed): \|Hc\|/N via `npv_none <- mean(!true_H)` in `quarto/dina/{subgroup_method_signature_recovery,dina_vs_forestsearch_signature_comparison}{,_2factor}.qmd` ×4; live package calls in `quarto/calibration/calibration_survival_fourway.qmd` ×2, `quarto/simulations/actg175/continuous/actg175_continuous_simulations.qmd` ×20, `vignettes/articles/paper_simulations.qmd:915-1168`; correct documentation in `quarto/guides/review_classification_calculations.qmd:218-346`, `quarto/simulations/classification_metrics_summary.qmd:181-189`; null-only statements (npv = 1 is right) in the `### Classification Under H0` prose ×123, `gbsg_020/current_status.md:161, 689`, `scripts_dinamr/nullthr_gateA.R:96-101` ×3; explicitly labelled detected-only tables (`sim_fs_mr_field_or_template.qmd:2021-2027` and 950 tracked renders, `summary_dinamr.qmd`, `summary_cert20.qmd`, `md_harm_mr_simulation_summary.qmd`, the 218-file cluster `quarto/simulations/gbsg/*_batched*.qmd`) | various | \|Hc\|/N and 1, or 1 under the null, or detected-only by explicit and labelled subset | — | — | — |
| fs-glms-interpretable | G1 `dev/payloads/R/identifiers_metrics.R:359-363` (`fill <- c(sens = 0, spec = 1, ppv = 0, npv = 1)`; `emit_mean(… ifelse(det, r[[k]], fill[[k]]), N)`) | code | npv **1**; spec 1 | **yes**: 211 of 369 tracked bundles have non-declaring replicates, all with `n_true` recorded | yes: `r$n_true` and `m$n_sample` | `dev/payloads/identifiers/identifiers_metrics.csv` (63 `npv_uncond` rows) → `manuscript/_R/oc_figures.R:390-395` (`oc_tbl_classification`, column "Npv") → supplement Table `tbl-supp-classification` |
| fs-glms-interpretable | G2 `dev/payloads/R/identifiers_summary.R:147-155` (doc `:11-14`) | code | npv **1** | yes | yes (`n_true`, `n_sample`) | `identifiers_summary.csv` (322 rows, `npv_uncond`, `npv_uncond_mcse`); read only by the reconcile scripts |
| fs-glms-interpretable | G3 `dev/payloads/R/idsweep_metrics.R:147-154` (doc `:31-32`) | code | npv **1** | yes (inflation up to 0.082, e.g. `dina_effMinSG` h100 n1500: 0.9666 vs 0.8844) | yes | `idsweep_metrics.csv` (288 `npv_uncond` rows); the manuscript reads this CSV but displays only sens/spec |
| fs-glms-interpretable | G4 display: `manuscript/_S4_operating_characteristics_supplement.qmd:168-177`; `manuscript/fs_glms_interpretable_supplement.tex:3018-3163` (×3 `.tex` incl. `dev/arxiv/companion_v1/{build,package}`) | quarto source / rendered | prints `npv_uncond` with the fill of 1, e.g. DINA 12.4% HR 1.50 n1500 **0.968** (\|Hc\|/N fill gives 0.943) | yes (every row is a harm cell) | no (display layer; the fix belongs in G1) | Table `tbl-supp-classification` in the supplement `.tex` and 8 tracked PDFs, including `manuscript/submissions/jasa_2026-10-13/{blinded,unblinded}_supplement.pdf` |
| fs-glms-interpretable | G5 checks and records that pin the old value: `dev/writeup/identifiers/R/reconcile_identifiers_summary.R:91-96` (asserts `npv_uncond == npv_cond*dr + 1*(1-dr)`); `dev/writeup/identifiers/R/verify_oc_section.R:772` (greps `idsweep_metrics.R` for the literal fill string) and `:488-489, 1015-1018`; `dev/writeup/identifiers/R/conv.R:4-10, 38` (console) | test / code | npv **1** | yes | reconcile: no (`identifiers_summary.csv` carries no per-replicate \|Hc\|); others via upstream | gates on the CSVs and on the supplement prose / Figure 2 caption; they fail once G1-G3 change |
| fs-glms-interpretable | G6 text stating the old convention: `dev/writeup/identifiers/identifier_evaluation_writeup.qmd:467-480, 870-871` (claims it matches León et al. 2024, Table 1); `dev/payloads/identifiers/README.md:248-249, 274`; `REQUEST_fs_classification_metrics_2026-09-12.md:131, 161, 218` and 7 other decision records; `dev/manuscript_status_update.md:965-967`; `REPORT_oc_comprehensive_2026-09-13.md:322` (quotes inflated values); `dev/writeup/identifiers/draft/_0Y_recovery_and_location.qmd:9-24, 41-45` ×2 (unconditional prose over conditional numbers; tracked `drafted_sections.pdf`); `manuscript/_05_operating_characteristics.qmd:25` (wording only; the main paper shows no NPV) | documents | npv 1 (stated) | n/a | n/a | records, drafts, one tracked draft PDF |
| fs-post-selection | P1 `manuscript_rev1/postselection_subgroup_inference_supplement.tex:2217, 2293, 2365, 2483, 2581, 2656, 3174, 3255` (and `manuscript/` frozen copy) | manuscript text | **detection-conditional** NPV (0.89-0.96), each labelled "over the … detected replicates" | no: no npv = 1 is reachable; the generators' per-replicate rows store NA for a non-detection and the summaries subset to detected | n/a | Supplement Tables S8.1-S8.3, S8.5 captions, S8.4 prose, S8.7 footnotes. The main text shows no NPV |
| fs-post-selection | P2 generators (forestsearch `quarto/simulations/gbsg_020/` templates, `.classify()` + detected-only summaries; cluster `dd2f832d` copies in `analysis/simulations/gbsg/resampling_comparisons/` ×18) | quarto source | NA stored; detected-only summaries | no | yes (`n_true`) | rev1 payloads (detected-only) |

### 2.4 Stored payloads and displayed values

- **forestsearch, unconditional NPV stored with the old value: 62 tracked payloads.**
  - 106 tracked `.rds` files store `$oc$identification` from `fs_mr_oc_summary()`. Recomputing each with the branch's `fs_mr_oc_summary()` changes the stored unconditional npv in **62** (40 in `quarto/simulations/actg175/binary_020/mr_or_harm/`, 6 in `mr_or_harm_superseded_prev09632/`, 16 in `quarto/simulations/actg175/continuous/mr_md_harm/`), by 0.00015 to 0.074 (mean 0.009). The other 44 are null cells or have no non-detected replicate, so their stored value is already right.
  - **No tracked document displays it.** Every tracked reader of `$oc$identification` (`dev/glm-continuous-sims/oc_wrapper_*.R`, `oc_breadth_stage2_score_2026-08-31.R:261`) takes the `conditional` row.
  - The 1,266 per-replicate MR payloads store `npv = NA` for a non-detected replicate (not 1), so they hold no old value.
- **forestsearch, displayed old value (pre-A-003 extractor convention, npv = 1 under H1):** the 17 actg175 renders (Q6) and the four tracked text dumps of `paper_simulations` (Q7). The sources are already correct (A-003); addressing these needs a re-render.
- **fs-glms-interpretable:** three tracked summary CSVs store `npv_uncond` built with npv = 1 (`identifiers_metrics.csv` 63 rows, `identifiers_summary.csv` 322 rows, `idsweep_metrics.csv` 288 rows). **One is displayed:** `identifiers_metrics.csv` feeds the "Npv" column of supplement Table `tbl-supp-classification`, in 3 tracked `.tex` files and 8 tracked PDFs, including the JASA submission supplement (2026-10-13).
- **fs-post-selection:** no payload stores, and no document displays, an NPV computed with npv = 1. Every NPV shown is detected-only and labelled so.

## 3. Summary for Larry

| Instance | Fixed here / BLOCKED / outside the package | Reported output affected | What addressing it needs |
|---|---|---|---|
| `fs_mr_oc_summary()` unconditional npv (`R/fs_mr_oc_summary.R`, roxygen + `.Rd`) | **Fixed here** | `$oc$identification` in 62 tracked payloads (not displayed anywhere tracked) | Done in source. The 62 payloads keep the old value until their batch documents re-run or are re-summarised (re-run, or a one-off re-summary of the stored `results`) |
| `build_classification_table()` fill (`R/simulation_tables.R`) | **Fixed here** where the row carries `propHc_true`; **BLOCKED sub-case:** a row with no truth column keeps the fill of 1 (it has no \|Hc\|/N; `propHc_true` is NA too) | `ppv(Hc)` rows, only for rows stored with npv = NA (pre-`95e2f819` output, other producers). No tracked payload has such a row with `propHc_true` < 1 (reviewer scan) | Done in source |
| `format_oc_results()` comment (`R/oc_analyses.R:683`) | **Fixed here** (comment only; the code was already correct) | none | Done |
| Q1, Q3 `oc_agg_from_results()`, `reproduce_table1_m1` fills | outside the package | bubble / Pareto plots; reproduction table (latent with current output) | source edit: take `propHc_true` instead of 1 |
| Q2, Q8 prose stating npv = 1 / "matches `fs.estimates.out()`" | outside the package | rendered prose | source edit (+ re-render) |
| Q4 batch templates calling `fs_mr_oc_summary()` | outside the package; fixed through the package | the 62 payloads above | re-run (or re-summarise) |
| Q5 actg175 binary_legacy summary over pre-A-003 bundles | outside the package | OC table, Spec × NPV plots (render not tracked) | re-run the bundles (or re-score them) |
| Q6, Q7 17 actg175 renders; `paper_simulations` text dumps | outside the package | displayed ppv(Hc) / sens(Hc) | re-render (sources correct since A-003) |
| Q9 NA then `na.rm` (grf_dmin_calibration, validation_glm, archive working, dina signature simulation) | outside the package | OC tables labelled as over all replications | source edit (\|Hc\|/N and spec 1 at the no-find return) + re-render |
| Q10 mr_coverage_sweep ×54, legacy extractor copies | outside the package | per-replicate bundles only | none needed unless summarised unconditionally |
| Q11 gbsg_020 reports with unlabelled detected-only NPV | outside the package | report text | wording edit ("over detected replicates") |
| G1-G3 fs-glms-interpretable payload builders | outside the package | `identifiers_metrics.csv` → **supplement Table tbl-supp-classification "Npv"**, in the JASA submission PDFs; two further CSVs | source edit (fill \|Hc\|/N = 1 − `n_true`/`n_sample`), re-run the payload scripts, re-render the supplement |
| G5 fs-glms-interpretable gates pinning npv = 1 | outside the package | gates | source edit, alongside G1-G3 |
| G6 fs-glms-interpretable text stating the old convention | outside the package | records, one draft PDF, supplement wording | text edit |
| P1-P2 fs-post-selection | outside the package; **no instance** (all detected-only, labelled) | none | nothing |

## 4. Part B — the change

### 4.1 What changed, by site

- **`fs_mr_oc_summary()`** (`R/fs_mr_oc_summary.R`).
  - **R1.** In the unconditional row, a non-detected replicate (`status` not `"DETECTED"`) now contributes npv = 1 − `n_true`/`meta$n_sample` = \|Hc\|/N; under the null (`n_true` = 0) that is 1. spec stays `fill("spec", 1)`. `sens`, `ppv`, `detection`, `n_used`, `mean_size_H` and the whole conditional row are unchanged. A **detected** replicate with no stored npv is still filled with 1, as before.
  - **R2.** \|Hc\|/N comes from the payload's per-replicate `results$n_true` and `meta$n_sample`, both already inputs of the exported function. No user-facing argument, payload column or input was added.
    - **Stated plainly:** the unexported, `@noRd` helper `.fs_oc_identification()` gained a parameter `n_sample = NULL`, fed from the existing `meta$n_sample` at its only call (`:188`). Read literally, R2's "Do not add arguments" forbids this. I read it as forbidding new inputs, not plumbing an existing input to an internal helper. The alternative, moving the fill out of the helper, is a refactor the task also forbids. Flagged as open item 1.
  - **R3 (where \|Hc\|/N is missing).** The unconditional npv is `NA` when \|Hc\|/N is unknown for any non-detected replicate: no numeric `n_true` column, no usable `meta$n_sample` (missing, non-numeric, not length 1, or ≤ 0), or `n_true` NA on a non-detected row (e.g. a harness `CONFIG-ERROR` row that failed before the data were simulated). No other value is substituted. No tracked MR payload hits this (all 1,266 carry numeric `n_true` and `meta$n_sample`, and none has an error row).
    - A payload with **no `status` column** marks no replicate and fills nothing, as on `master`, so its unconditional npv is unchanged (the reviewer counted 595 such tracked dev payloads).
  - Roxygen and the comment above the fill now state the definition, cite León et al. (2024), Section 3, and say how NA arises. The claim "following `build_classification_table()` … Table 1" was removed.
- **`build_classification_table()`** (`R/simulation_tables.R`).
  - **R1/R2.** Under an alternative, a row stored with npv = NA and `any.H == 0` now gets its own `propHc_true` (\|Hc\|/N, a size in the same row) instead of 1. sens, ppv and spec fills are unchanged.
  - **R3 / BLOCKED sub-case.** A row without `propHc_true`, or with it NA (no truth column: then the producer stores no truth-based value at all), keeps the fill of 1. It lacks a per-replicate true-subgroup size or flag. The comment says this, describing the code.
  - The comment's false statement about what `run_simulation_analysis()` stores was corrected.
- **`format_oc_results()`** (`R/oc_analyses.R:683`): comment only, now says the average is over every sim with a stored value and that a no-subgroup sim contributes spec = 1 and npv = \|Hc\|/N. No code change.
- `devtools::document()` ran once (after the last roxygen edit). It rewrote only `man/fs_mr_oc_summary.Rd`.
- `DESCRIPTION`: 0.4.1.9004 → 0.4.1.9005, the only line changed.
- `NEWS.md`: one entry, the first bullet under `### Behaviour changes` in the development section (whose heading still reads `0.4.1.9002`, as before; A-003 open item 6). It names both functions changed and states: a replicate with no identified subgroup now contributes NPV = \|Hc\|/N (it was 1); specificity, sensitivity, PPV and detected-only summaries are unchanged; stored payloads are not regenerated. It also states the NA rule and the BLOCKED sub-case.

### 4.2 `git diff master -- R/`, verbatim

```diff
diff --git a/R/fs_mr_oc_summary.R b/R/fs_mr_oc_summary.R
index 221a18d4..6eb406e7 100644
--- a/R/fs_mr_oc_summary.R
+++ b/R/fs_mr_oc_summary.R
@@ -51,9 +51,17 @@
 #' \describe{
 #'   \item{\code{conditional}}{Averaged over DETECTED replicates only. This is
 #'     the detected-only summary.}
-#'   \item{\code{unconditional}}{Averaged over ALL replicates, following
-#'     \code{build_classification_table()} and Leon et al. (2024, Table 1),
-#'     with non-detection scored as sens = 0, ppv = 0, spec = 1, npv = 1.}
+#'   \item{\code{unconditional}}{Averaged over ALL replicates. A replicate with
+#'     no identified subgroup has everyone in its estimated complement, so it
+#'     is scored as sens = 0, ppv = 0, spec = 1 and npv = |Hc|/N, the share of
+#'     that replicate's subjects truly in the complement (Leon et al. 2024,
+#'     Section 3), computed as \code{1 - n_true / meta$n_sample}.  Under the
+#'     null this is 1.  These are the spec and npv values
+#'     \code{run_simulation_analysis()} stores for such a replicate.  The
+#'     unconditional npv is \code{NA} when the payload has no numeric
+#'     \code{n_true} column or no usable \code{meta$n_sample}, or when
+#'     \code{n_true} is missing for any non-detected replicate (for example
+#'     a replicate that failed before its data were simulated).}
 #' }
 #' With detection at 999/1000 the two agree to three decimals; in a null cell
 #' at 0.930 they do not. Quote which one you mean.
@@ -179,7 +187,7 @@ fs_mr_oc_summary <- function(payload,
 
   out <- list(
     estimation     = est_tab,
-    identification = .fs_oc_identification(res, rd, digits),
+    identification = .fs_oc_identification(res, rd, digits, meta$n_sample),
     targets        = list(
       orient      = orient,
       struct_H    = orient * truth$effect_Q,
@@ -277,18 +285,31 @@ print.fs_mr_oc <- function(x, ...) {
 #' Classification metrics under both conventions
 #' @keywords internal
 #' @noRd
-.fs_oc_identification <- function(res, rd, digits = NULL) {
+.fs_oc_identification <- function(res, rd, digits = NULL, n_sample = NULL) {
   n_all <- nrow(res)
   det   <- res$status %in% "DETECTED"
 
-  # Unconditional: non-detection scores sens = 0, ppv = 0, spec = 1, npv = 1,
-  # matching build_classification_table() and Leon et al. (2024, Table 1).
+  # Unconditional: non-detection scores sens = 0, ppv = 0, spec = 1 and
+  # npv = |Hc|/N = 1 - n_true / n_sample, the share of that replicate's
+  # subjects truly in the complement (Leon et al. 2024, Section 3; 1 under
+  # the null), as run_simulation_analysis() stores it.  A detected replicate
+  # with no stored npv is still filled with 1.  NA when |Hc|/N is unknown
+  # for any non-detected replicate (no n_true, or no usable n_sample).
   fill <- function(col, miss) {
     v <- res[[col]]
     if (is.null(v)) return(rep(NA_real_, n_all))
-    v[!det | is.na(v)] <- miss
+    miss <- rep_len(miss, n_all)
+    v[!det | is.na(v)] <- miss[!det | is.na(v)]
     v
   }
+  hc_share <- if (!is.numeric(res$n_true) || !is.numeric(n_sample) ||
+                  length(n_sample) != 1L || !isTRUE(n_sample > 0)) {
+    rep(NA_real_, n_all)
+  } else {
+    1 - res$n_true / n_sample
+  }
+  npv_u <- fill("npv", ifelse(det, 1, hc_share))
+  npv_unknown <- any(!det & is.na(hc_share))
 
   out <- rbind(
     data.frame(
@@ -306,7 +327,8 @@ print.fs_mr_oc <- function(x, ...) {
       detection   = mean(det),
       mean_size_H = .fs_oc_mean(rd$n_harm),
       sens = .fs_oc_mean(fill("sens", 0)), spec = .fs_oc_mean(fill("spec", 1)),
-      ppv  = .fs_oc_mean(fill("ppv",  0)), npv  = .fs_oc_mean(fill("npv",  1)),
+      ppv  = .fs_oc_mean(fill("ppv",  0)),
+      npv  = if (npv_unknown) NA_real_ else .fs_oc_mean(npv_u),
       stringsAsFactors = FALSE
     )
   )
diff --git a/R/oc_analyses.R b/R/oc_analyses.R
index d8015e71..57c3c6d4 100644
--- a/R/oc_analyses.R
+++ b/R/oc_analyses.R
@@ -680,7 +680,9 @@ format_oc_results <- function(
     # Detection rate
     detection_rate <- mean(res$any.H, na.rm = TRUE)
 
-    # Classification metrics (averaged across sims with detection).
+    # Classification metrics (averaged over every sim with a stored value: a
+    # sim with no subgroup contributes spec = 1 and npv = |Hc|/N as
+    # run_simulation_analysis() stores them; NA rows drop out).
     # Under null models (León et al. 2024, Sec 3): Q = emptyset,
     # so sens = NA and ppv = 0.  NaN from mean(all-NA) is replaced
     # with NA for clean table display.
diff --git a/R/simulation_tables.R b/R/simulation_tables.R
index f4da7007..50714f84 100644
--- a/R/simulation_tables.R
+++ b/R/simulation_tables.R
@@ -362,14 +362,24 @@ build_classification_table <- function(
         # across ALL simulations (unconditional).  Under alt, when Ĥ = ∅:
         #   sens(Ĥ) = #{Ĥ ∩ H}/#{H} = 0/|H| = 0  (H exists but not found)
         #   ppv(Ĥ)  = #{Ĥ ∩ H}/#{Ĥ} = 0/0         (convention: 0)
+        #   sens(Ĥc) = 1 and ppv(Ĥc) = |Hc|/N     (Ĥc is everyone; Section 3)
         # Under null, sens is genuinely undefined (H = ∅) → keep NA.
-        # The raw values from run_simulation_analysis store NA for both
-        # sens and ppv when any.H = 0; we replace with 0 under alt only.
+        # run_simulation_analysis() already stores these values when no
+        # subgroup is found and the truth column is present, so the fills
+        # below change nothing there.  They act only on rows whose values
+        # are NA (no truth column, output of an older package version, or
+        # another producer), and only under alt.  Such a row has npv filled
+        # with its propHc_true (|Hc|/N for that replicate) when the row
+        # carries it, and with 1 when it does not (no truth column: then
+        # propHc_true is NA too).
         if (hyp != "null") {
           r$sens[is.na(r$sens) & r$any.H == 0] <- 0
           r$ppv[is.na(r$ppv) & r$any.H == 0]   <- 0
           r$spec[is.na(r$spec) & r$any.H == 0]  <- 1
-          r$npv[is.na(r$npv) & r$any.H == 0]    <- 1
+          nf <- is.na(r$npv) & r$any.H == 0
+          hc <- if ("propHc_true" %in% names(r)) r$propHc_true else rep(NA_real_, nrow(r))
+          r$npv[nf & !is.na(hc)] <- hc[nf & !is.na(hc)]
+          r$npv[nf & is.na(hc)]  <- 1
         }
 
         val <- switch(
```

### 4.3 Hunk-to-row mapping

| Hunk | Part A row |
|---|---|
| `R/fs_mr_oc_summary.R` `@@ -51,9 +51,17 @@` | roxygen `unconditional` item (`:48-57`) |
| `R/fs_mr_oc_summary.R` `@@ -179,7 +187,7 @@` | `.fs_oc_identification()` row: the call passes `meta$n_sample` (R2 plumbing) |
| `R/fs_mr_oc_summary.R` `@@ -277,18 +285,31 @@` | `.fs_oc_identification()` row: comment, `fill()`, `hc_share`, `npv_u`, `npv_unknown` |
| `R/fs_mr_oc_summary.R` `@@ -306,7 +327,8 @@` | `.fs_oc_identification()` row: the npv cell of the unconditional row (`:309` on `master`) |
| `R/oc_analyses.R` `@@ -680,7 +680,9 @@` | `format_oc_results()` comment row (`:683`) |
| `R/simulation_tables.R` `@@ -362,14 +362,24 @@` | `build_classification_table()` row (`:360-372`) |

### 4.4 Part B gates

| Gate | Result |
|---|---|
| 1. `git diff master --stat -- R/` lists only files with a Part A row | `R/fs_mr_oc_summary.R` (+31/-9), `R/oc_analyses.R` (+3/-1), `R/simulation_tables.R` (+13/-3). All three have rows. **PASS** |
| 2. Every hunk in `git diff master -- R/` maps to a row | 6 hunks, all mapped (§4.3). **PASS** |
| 3. `document()` ran: `git diff master --stat -- man/` lists only Rd files of functions whose roxygen changed | `man/fs_mr_oc_summary.Rd` only. **PASS** |

## 5. Part C — tests

### 5.1 Updated expectations

| File:line | Old value | New value |
|---|---|---|
| `tests/testthat/test-fs_mr_oc_summary.R:108` (now `:110`) | `0.85` = (0.7 + 1)/2 | `(0.7 + 459/700) / 2` = 0.677857. The fixture `make_payload()` gained `n_true = 241` (`:17`), so \|Hc\|/N = 459/700 with `meta$n_sample` = 700 |

No other expectation in the suite pins NPV = 1 for a non-detected replicate under an alternative (independent grep, §7). `test-nofind-classification.R:91` and `test-fs-oc-null.R:77` pin 1 under the null, which is correct.

### 5.2 New expectations

`tests/testthat/test-fs_mr_oc_summary.R` (function changed: `fs_mr_oc_summary()`):

| Block (line) | What it checks |
|---|---|
| `:174` alternative | Hand-built payload, 2 detected + 2 non-detected, N = 200, `n_true` 50 and 60 for the non-detected. Unconditional npv = (0.9 + 0.8 + 150/200 + 140/200)/4 and spec = (0.85 + 0.95 + 1 + 1)/4; sens and ppv unchanged. One non-detected replicate on its own: npv = 150/200, spec = 1. Conditional row unchanged |
| `:196` null | The same input with `n_true` = 0: unconditional npv = 1; non-detected replicates alone: npv = 1, spec = 1 |
| `:210` NA | npv NA when `n_true` is absent, `meta$n_sample` is absent, or one non-detected `n_true` is NA; conditional npv and unconditional spec unaffected |
| `:226` **agreement** | Four N = 200 replicates scored by the `run_simulation_analysis()` extractor `.extract_fs_estimates_gen()` (two found through `treat.recommend`, two with nothing found; `any.H` checked = 1, 1, 0, 0). The MR-style payload of the same replicates (found-case values for detected, NA for non-detected, as the harness stores them): `fs_mr_oc_summary()`'s unconditional npv equals the all-replicate mean of the extractor npv, and the same for spec. Built by hand, no search |
| `:264` | A payload with no `status` column keeps its unconditional npv (as on `master`); a non-numeric `n_sample` gives NA, not an error |

`tests/testthat/test-build_classification_table.R` (new file; function changed: `build_classification_table()`; no test file existed):

| Block | What it checks |
|---|---|
| alternative | 2 found + 2 no-find rows stored with NA (as older output) and `propHc_true` 0.75, 0.70: `ppv(Hc)` = (0.9 + 0.8 + 0.75 + 0.70)/4, `sens(Hc)` = (0.85 + 0.95 + 1 + 1)/4; sens(H), ppv(H) unchanged. One no-find row alone: 0.75 and 1. Rows already carrying values (current extractors) give the same table |
| null | The same input under the null (stored spec 1, npv 1 for no-find, as the extractors store them): `ppv(Hc)` = 1; no-find rows alone: 1 and 1 |
| BLOCKED sub-case | a row with `propHc_true` NA, or no such column: npv filled with 1 (pins the documented behaviour) |

- No search, no simulation, no parallel plan in any test.
- **Negative control.** Both files were run against `master`'s code (`git archive master` copy). `test-fs_mr_oc_summary.R`: 8 expectations fail, all on unconditional npv (blocks `:108`, alternative, NA, agreement, status/non-numeric). `test-build_classification_table.R`: 4 fail (alternative and BLOCKED blocks). The null blocks pass on both, as they should.

## 6. Part D — tests and check

| Run | Baseline (A-003 report §4, not re-run) | This branch | Result |
|---|---|---|---|
| Changed/new test files (`load_all` + `test_file`) | — | `test-fs_mr_oc_summary.R` 18 blocks, 57 expectations; `test-build_classification_table.R` 3 blocks, 13 expectations; **0 failures**, 0 warnings, 0 skips; wall 4 s together (cap 2 min each) | **PASS** |
| Full suite (`devtools::test()`) | 600 blocks, 80 files, 8,183 expectations, 0 failures, 3 skips, 38 warnings; 642 s | **608 blocks, 81 files, 8,218 expectations, 0 failures, 0 errors, 3 skips, 38 warnings**; wall **647 s** (cap 30 min) | **PASS** |
| Difference | — | +8 blocks, +1 file, +35 expectations, all from the two test files (13 new file + 22 added to `test-fs_mr_oc_summary.R`); skips (the same 3: `test-cross-outcome-parity.R` ×2, `test-cv-no-subgroup-edges.R` ×1) and warnings unchanged | |
| `devtools::check()` | 0 errors / 0 warnings / 0 notes (export route); 994 s | **0 errors / 0 warnings / 0 notes**, Status OK; wall **1,005 s** (R CMD check duration 14 m 58 s; tests 10 m, vignettes 77 s) (cap 60 min) | **PASS** |

- **Check route:** as A-003 §4: `check(document = FALSE, error_on = "never")` on an export of the working tree (`git ls-files -co --exclude-standard`, 9,624 files) in the session scratchpad; the 8 changed/new files were `cmp`-identical to the repository's. The export avoids the unrelated `legacy_vignette/` NOTE.
- **Timing:** the suite and `check()` ran at the same time (both started 15:26:33), so each wall clock includes contention from the other.
- **A first suite/check pair was stopped and discarded.** It started on an earlier version of the edits, before the independent review (§7) led to the changes in §7's table. Its results were not used.
- **Pre-set decision:** no failing expectation pinned NPV = 1 in the suite run; nothing needed updating beyond §5.1.
- **Commit:** `fix: NPV for a replicate with no subgroup is |Hc|/N in all package summaries` on `fix/npv-nofind-instances`, by explicit path: `DESCRIPTION`, `NEWS.md`, `R/fs_mr_oc_summary.R`, `R/oc_analyses.R`, `R/simulation_tables.R`, `man/fs_mr_oc_summary.Rd`, `tests/testthat/test-fs_mr_oc_summary.R`, `tests/testthat/test-build_classification_table.R`.

## 7. Independent review

Three read-only reviewers checked the uncommitted change against the task document, with no access to these conclusions. A separate read-only inventory workflow (8 batch readers, 3 sibling/payload readers, one adversarial checker per area, one completeness critic; 23 agents) produced §2.3.

| Lens | Verdict on the first version | What changed in response |
|---|---|---|
| Spec conformance | R1, Part C and gates B1-B3 pass. Two findings. (F1) The first NA rule (`anyNA` over all rows) turned a payload with no `status` column from a number into NA (595 tracked dev payloads). (F2) The BLOCKED reason for `build_classification_table()` was wrong: rows carry `propHc_true` | (F1) NA now only when \|Hc\|/N is unknown on a **non-detected** row; status-less payloads keep `master`'s value (test `:264`). (F2) The site now uses `propHc_true`; BLOCKED only where it is absent. Roxygen "these are the values `run_simulation_analysis()` stores" narrowed to spec and npv |
| Statistical correctness | PASS. 1 − `n_true`/`n_sample` is ppv(Ĥc) on an empty Ĥ. `nrow(df)` = `n_sample` in the harness (follow-up ≥ 60 at analysis time 84), confirmed on 635,141 detected rows of 920 payloads to 1.1e-13. Old → new on three alternative payloads: 0.95144 → 0.89263, 0.97060 → 0.88934, 0.80597 → 0.78330; null 1 → 1. Low findings: a character `n_sample` would error; `n_sample` = 0 would give −Inf | Guards added: non-numeric or non-positive `n_sample`, non-numeric `n_true` → NA (tested) |
| R/ inventory completeness | Inventory complete apart from: the `format_oc_results()` comment (`R/oc_analyses.R:683`), now a Part A row and corrected; the CV no-find branch, listed as not a Part A row with open item 4; the `build_classification_table()` reason (as F2) | as stated |

## OPEN ITEMS

1. **R2 literal reading.** `.fs_oc_identification()` (unexported, `@noRd`, one caller) gained `n_sample = NULL` to carry the existing `meta$n_sample`. No user-facing argument, column or input was added. If the rule is meant literally, the alternative is to compute the unconditional npv in `fs_mr_oc_summary()` itself.
2. **The 62 tracked payloads** in `quarto/simulations/actg175/{binary_020/mr_or_harm, binary_020/mr_or_harm_superseded_prev09632, continuous/mr_md_harm}/` keep the old stored `$oc$identification` unconditional npv (§2.4). Nothing tracked displays it. Re-summarising them is a one-off read of each payload's own `results`; it was not done (no payload changes in this task).
3. **fs-glms-interpretable supplement Table `tbl-supp-classification`** displays `npv_uncond` built with npv = 1 (G1, G4), including the 2026-10-13 JASA submission PDFs. The fix is upstream in `dev/payloads/R/identifiers_metrics.R:359-363` (the same file also has `identifiers_summary.R` and `idsweep_metrics.R`), with the gates in G5 updated alongside. Read-only here.
4. **CV: folds find nothing, primary subgroup exists** (`R/forestsearch_cross_validation.R:1414-1420`, `forestsearch_KfoldOut()`). The single-level branch sets `sens_Hc` and `ppv_Hc` to NA and computes `sens_H` from row 1 of `table(treat.recommend, treat.recommend.original)`. When only level `1` is present, that row is the not-H level, so `sens_H` reports 1 for zero agreement and `ppv_H` reports \|orig H\|/N. `tests/testthat/test-cv-mixed-and-fatal.R:142` pins `c(1, NA, 0.35, NA)`. This is an agreement metric against the observed pick, not the León NPV, so it is not a Part A row and was not changed. `quarto/qc/cv_code_audit.md:243-246` already recommends sens_H = 0, sens_Hc = 1, ppv_Hc = observed Hc prevalence.
5. **`fs_mr_oc_summary()` CONFIG-ERROR rows.** A non-detected replicate with `n_true` NA (harness `CONFIG-ERROR` before simulation) makes the unconditional npv NA for the whole payload, while sens, ppv and spec still score that row as a non-detection. No tracked payload has such a row.
6. **Pre-existing, unchanged:** a **detected** replicate with Ĥ = everyone (n_harm = N) has npv 0/0. `fs_mr_oc_summary()`'s unconditional row fills it with 1; `run_simulation_analysis()` stores NA, which `na.rm` then drops. This affects 1,344 rows in 12 tracked idsweep payloads. Out of scope (detected replicates were not to change).
7. **Pre-existing, unchanged:** under the null, `fs_mr_oc_summary()`'s unconditional sens is the fill 0 where `build_classification_table()` keeps NA.
8. **Parity script.** `dev/glm-continuous-sims/parity_fs_mr_oc_summary.R` compares `fs_mr_oc_summary()` with the batch documents' inline machinery; it will now disagree on the unconditional npv. Outside `R/`; not run.
9. **Outside-package instances** in §2.3 / §3 are listed only. The quickest, with the largest displayed effect, are G1-G4 (a published table) and Q6/Q7 (re-renders whose sources are already correct).
10. **NEWS heading** still reads `0.4.1.9002 (development)` while `DESCRIPTION` is 0.4.1.9005 (A-003 open item 6).
11. **No coverage previously** for `build_classification_table()`, `format_oc_results()` or `summarize_simulation_results()`; the new test file covers only the no-find npv path of the first.
