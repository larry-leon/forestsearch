# TASK: forestsearch — the MRCT consistency article, and `mr_estimates_table()` without a bootstrap template

Date: 2026-11-01
Repository: **forestsearch**, branch `mrct-analytic`, from `5db25a27`. mrct-consistency (at the commit that closed `TASK_mr-field-campaign_2026-11-01`; record it) is read only, as the source of text, summaries and seeds.
Prerequisite: `results/mr_field_campaign/summary/mr_field_core.csv` exists in mrct-consistency. Stop at step 0 if it does not.
CC does not push.

**`R/` item, called out:** §1 **changes behaviour**, additively: `mr_estimates_table()` gains `boot_results = NULL` and new columns that appear only when field quantities exist. Existing calls are unchanged. Nothing changes the method. Everything else is a pkgdown article (`vignettes/articles/`, excluded from the CRAN build) with its data files.

## 0. First action

Copy this file into `dev/tasks/` and commit it on its own, on `mrct-analytic`.

## 1. `mr_estimates_table(fs, boot_results = NULL, est.scale = "hr", digits = 2)`

- With `boot_results` supplied: exactly today's behaviour (the existing example in `R/forestsearch_main.R` and the call in `vignettes/articles/survival-post-selection.qmd` render identically; check the latter's table by rendering that chunk or by a snapshot).
- With `boot_results = NULL`: build the descriptive columns (n, events, medians as `summarize_bootstrap_results()` lays them out, or the subset computable from `fs`'s analysis data and the declared rule) and the naive effect from `fs` itself, for Ĥ and for its complement in the analysis data. The naive limits are the robust (sandwich) limits in `fs$mr_inference$naive`, and the column header says "robust". Then the MR columns: the de-biased estimate and its IJ interval, and, when `fs$mr_inference$field` exists, the field one-sided 95% bound (lower for Ĥ, upper for the complement) and the Bonferroni pair. Complement rows only when `include_complement = TRUE` was used; otherwise the cell is marked as today.
- Roxygen updated; a test on a small `forestsearch(mr_inference = TRUE)` fit (the field method at reduced `field_R_out`/`field_R_in` for speed) that the table has the expected rows and columns and that its numbers equal the elements of `fs$mr_inference` they come from (1e-12). Export list unchanged. NEWS entry under 0.4.1.9000.

## 2. The article: `vignettes/articles/mrct-consistency.qmd`

Front matter as `paper_simulations.qmd` (Quarto HTML, embedded resources, toc, numbered sections, code folded). Title: "MRCT regional consistency via subgroup identification". Add it to `_pkgdown.yml` in both places (the navbar menu and the `articles:` sections) under a new section "MRCT consistency".

**Self-contained:** the article reads only files in `vignettes/articles/mrct-consistency/` (its data directory, created by this task) and runs forestsearch. It never reads mrct-consistency at build time. Nothing in the text is typed in: every number is computed in the article or read from a shipped CSV whose provenance (mrct-consistency commit, source path) is listed in `vignettes/articles/mrct-consistency/PROVENANCE.md`.

**Shipped data (small CSVs, copied verbatim with their commits named):**
- the committed S2 build's latent correlation `latent_R` (write it as CSV from `build_full("S2")$mrct$latent_R`) and the fitted AFT parameters as a check target (`results/ident_rlog4_full/summary/build.csv`, the S2 row);
- the core OC tables: `results/ident_rlog4_S1/summary/core_n1000.csv`, and the later-enrolling cells' core table from `results/regional_entry_cells/summary/`;
- the analytic predictions: `results/ident_rlog4_full/analytic/predictions.csv`, `results/ident_rlog4_S1/analytic/predictions.csv`, `results/ident_rlog4_n700/analytic/predictions.csv`, `results/regional_entry/roots.csv`;
- the MR field summary: `results/mr_field_campaign/summary/mr_field_core.csv`;
- the replicate-5 record for the reproduction check: the four HRs and the 2 × 2 (`results/s2_example/`).

**Terminology** as the technical report: H*, Ĥ, Ĝ = region minus Ĥ, G* = region minus H*, Ĥᶜ = the non-region complement of Ĥ; the trigger is ITT ub < 1.0 and region ITT > 1.0; "joint"; cells in concrete form; built from ACTG 175; AHR, CDE, marginal hazard ratio; bounds read by location.

**Sections, in this order** (text condensed from the technical report's chapters; code written out and run):

1. **Overview: the question and the design.** The trigger and the pipeline (report Chapter 2, in a page). The DGM built from ACTG 175 with `inject_mrct_structure(..., latent_R = <shipped>)` and `scale_dgm_hazard(m = 3)`, every call visible; a check that the build matches the shipped build record (stop otherwise). The analytic layer in brief, then live: `ahr_population()`, the four truths per set and cell (`ahr_limit()`), `ahr_design()` at n = 1,000 for the six cells — the trigger predictions and the oracle joint — and the same at n = 700 (predictions only; say that the n = 700 campaign exists in the report and is not reproduced here); how r = −log 4 was chosen (`ahr_solve()` on r at h = 2, reproducing the committed root −log 4.226); the two later-enrolling-region roots reproduced with `entry_delay`. State the cells' specifications as a table.
2. **Replicate 5.** `simulate_from_dgm()` at seed 613800005 with the campaign's arguments; `forestsearch()` on the non-region patients with the campaign's arguments and `mr_inference = TRUE` (field, complement, seed 613800005; one call, about 50 s); Ĥ; the Pareto frontier (`plot_pareto_frontier(..., show_band = TRUE, effect_neighborhood = 0.20)` with `compute_frontier_cis()`); `mr_estimates_table(fs)` from §1; the four sets with their Cox HRs (one table), the forest plot and the 2 × 2 Kaplan–Meier figure; Ĥ against H* in the region and outside; the reproduction check against the shipped record (stop otherwise); the analytic layer's limits for the four sets beside the estimates.
3. **Simulation results.** The core OC table at n = 1,000 (six cells) and the later-enrolling cells' table, read from the shipped CSVs; predicted against observed for the trigger and the oracle joint; the MR field rows (the three shares) read from the shipped summary; then a short reading for each. The full campaign code (`mrct_region_sims()` with every argument, `record = "full"`, the later-enrolling cells with `entry_delay`) shown with `eval: false`, with the wall-clock of the committed runs stated from their timing records. One live demonstration: `mrct_region_sims(record = "full")` on S2 with a small replicate count set at the top of the file (as `paper_simulations.qmd` does with `nsims_alt`; default 100) and the workers from `parallel::detectCores()`, summarised into the core OC rows by the article's own code, labelled as a demonstration.
4. **Appendix: the analytic tools.** The development condensed from report Chapter 6: the counting-process setting and the least-false limit, the Struthers–Kalbfleisch root for a patient set, the clock identity and the frailty reading, the first-order and events forms, the design layer (information, √(information share) correlation, the bivariate-normal trigger probability), per-patient entry; with the formulas and a pointer to the report for derivations.

## 3. Build and checks

- `quarto render vignettes/articles/mrct-consistency.qmd` completes (record the time; the live simulation and the MR call dominate); the HTML is git-ignored as the other articles' are.
- `pkgdown::build_article("mrct-consistency")` completes; do not build the whole site.
- All `stopifnot` checks in the article hold (build record; replicate-5 HRs and Ĥ; the roots; the core tables' recount where the CSV carries counts).
- `devtools::document()` clean; `devtools::install()`; the §1 test passes; export list unchanged.

## 4. Notes file

`~/Downloads/NOTES_fs-mrct-article_2026-11-01.md`: status, commits, render time, the article's section list with any place the report's text was changed in substance, the `mr_estimates_table()` test result, and any disagreement with this task.

## 5. Commits (on `mrct-analytic`)

Three: this task; §1 with its test, docs and NEWS; the article, its data directory and `_pkgdown.yml`.
