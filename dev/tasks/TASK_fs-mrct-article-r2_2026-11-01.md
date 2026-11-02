# TASK: forestsearch — article r2: the reader's path, case-study helpers, and the Appendix equations

Date: 2026-11-01
Repository: **forestsearch**, branch `mrct-analytic`, from the second commit of `TASK_fs-exact-betaHhat-mr-record_2026-11-01_r2` (record it). mrct-consistency read only. CC does not push.

**`R/` item, called out:** §2 adds two exports, `mrct_select_case()` and `mrct_case_study()`. New functions only (additive); no existing function changes; nothing changes the method. Everything else is the article (`vignettes/articles/mrct-consistency.qmd`).

## 0. First action

Copy this file into `dev/tasks/` and commit it on its own, on `mrct-analytic`.

## 1. Purpose

Larry's review of the rendered article (2026-11-01):
- it must be self-contained for a reader who substitutes their own trial dataset, with wrapper code showing how a simulated case study is found the way replicate 5 was;
- `latent_R.csv` reads as specific to this record;
- three display equations in the Appendix spill past the column.

## 2. The helpers (`R/mrct_case_study.R`)

**`mrct_select_case(res, conditions = c("trigger", "declared", "region_Ghat_lt1"), which = c("first", "all"), itt_ub = 1.0, region_itt = 1.0)`**
- `res`: a `mrct_region_sims(record = "full")` result (its `results` table and `estimators` attribute).
- Evaluates, per trial, the trigger (ITT upper 95% limit < `itt_ub` and region HR > `region_itt`, from the `itt_all` and `itt_region` estimator rows), `declared` (`any_found == 1`), and region(Ĝ) < 1.0 (`region_Ghat` row); a missing HR never satisfies a condition. Returns the selected `sim` index or indices with the three flags and the trial's four HRs (`itt_all`, `itt_region`, `region_Ghat`, `region_Gstar` when present). `which = "first"` is the rule that chose replicate 5 (`sims/150_s2_example_km.R` in mrct-consistency: the first trial by index meeting all three).
- Stops with a clear message if `res` has no `estimators` attribute (`record = "point"`).

**`mrct_case_study(dgm, seed, k, n, region_var, fs_args, sim_args, analysis_time, cens_adjust, mr_inference = TRUE, mr_inference_args = list(ci_method = "field", include_complement = TRUE, field_complement = TRUE), entry_delay = NULL)`**
- Redraws trial k exactly as `mrct_region_sims()` draws it (`simulate_from_dgm(seed = seed + k, ...)` with the same arguments; same `entry_var`/`entrytime` convention; document that the campaign's workers drew under L'Ecuyer-CMRG and reproduce that), runs `forestsearch()` on the non-region patients with the campaign's arguments plus the MR arguments (MR seeded with `seed + k`), forms the four sets (all, region, region(G*) when `flag_harm` exists, region(Ĝ) by `get_dfpred()`), fits the unadjusted Cox model on each, and computes Ĥ against H* on both sides when `flag_harm` exists.
- Returns an object of class `mrct_case_study` holding the trial, the fit, the four-set table (n, events by arm, HR, Wald 95%), the 2 × 2 tables, and the MR object; `print()` gives the four-set table and the MR bounds for Ĥ and its non-region complement; `plot()` gives the forest plot of the four sets and, with `weightedsurv` available (Suggests; skip with a message otherwise), the 2 × 2 Kaplan–Meier figure; `mr_estimates_table(x$fs)` works on it.
- Tests (`tests/testthat/test-mrct-case-study.R`): on the S2 build with the campaign's call seed, `mrct_select_case()` on a `record = "full"` run of replicates 1–6 returns 5; `mrct_case_study(k = 5)` reproduces the four committed HRs (shipped fixture from `results/s2_example/`, 3 dp) and Ĥ; a run with no `flag_harm` omits the G* and 2 × 2 parts with a message; `print()` and `plot()` run without error. Reduced MR settings for speed where the check does not depend on them.
- Roxygen, NEWS under 0.4.1.9000; export list grows by exactly these two names.

## 3. The article

**3.1 The reader's path (§1.3 rewritten).** Build the DGM the way a reader would, with every call visible:
1. `actg <- subset(speff2trial::ACTG175, arms %in% c(1L, 3L))` (the published trial data; name the variables used);
2. `inject_mrct_structure(...)` with the campaign's arguments written out, fitting the latent correlation live (latentcor in Suggests; the fallback warning explained in one sentence);
3. `scale_dgm_hazard(dgm, m = 3)`;
4. a short list, "To use your own trial": the data frame and its outcome, event and treatment names; the covariates; the region model; H* and its (b, h); the analysis window and entry.

Then one collapsed chunk, "Reproducing this report's build exactly", that passes `latent_R = <latent_R.csv>` and feeds the build-record check and everything downstream. Say in one sentence why it exists (platform-exact reproduction of the committed record) and that a reader never needs it. If the live fit on the render machine equals the shipped matrix (the previous notes say it does), say so.

**3.2 New section 2, "Finding a case study"**, before the replicate-5 section (which becomes section 3). Four visible calls: the small `mrct_region_sims(record = "full")` run (move the live demonstration's run here, same seed and replicate count, so it is run once and reused in the results section), `mrct_select_case()`, `mrct_case_study()`, `print()`/`plot()`. State that the selection returns replicate 5 (check it) and that the rest of the article examines that trial in detail. In the results section, the demonstration subsection now refers back to this run.

**3.3 Replicate-5 section.** Keep as is, but where its code now duplicates `mrct_case_study()` (the redraw, the identifier call, the four sets), call the helper and keep the step-by-step code only where the article explains a step (the identifier's arguments; the MR table; the frontier). The reproduction check stays.

**3.4 Appendix equations.** Break into `aligned` multi-line displays: the set-level score integral; the three-form variance line; the trigger-probability line with its three z-scores. Add to the article's CSS `.math.display, mjx-container[display="true"] { overflow-x: auto; }` as a safety net. Check the render at a 1,280-pixel window with the same method as the review (every display element against the main column's right edge): none may overflow.

**3.5 Everything else unchanged**: checks, provenance file (add the two helpers' fixture if shipped), `_pkgdown.yml`.

## 4. Gates (stop on failure)

- `devtools::document()` clean; `devtools::install()`; export list grows by exactly `mrct_select_case`, `mrct_case_study`; the new test passes; `test-mr-estimates-table-nullboot.R` still passes.
- `quarto render vignettes/articles/mrct-consistency.qmd` completes; every `stopifnot` in the article holds, including the new "selection returns 5"; the overflow check of §3.4 reports none.
- `pkgdown::build_article("articles/mrct-consistency")` completes.

## 5. Notes file

`~/Downloads/NOTES_fs-mrct-article-r2_2026-11-01.md`: status, commits, render time, the overflow check's output, the helpers' test counts, and any disagreement with this task.

## 6. Commits (on `mrct-analytic`)

Three: this task; the helpers with tests, docs and NEWS; the article and its data directory.
