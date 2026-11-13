# REPORT — CRAN preparation, part 6: simulation effect scales, release 0.5.0 (2026-11-12)

- **Part 1 — PASS.** `run_simulation_analysis()` now reports RD, RR and IRD on their own scales (commit `9787bbb4`).
- **Part 2 — PASS.** None of the searched repositories has a `run_simulation_analysis()` run with RD, IRD or RR. 0 sites and 0 tracked payloads.
- **Part 3 — PARTIAL.** Done: steps 1–4 and 7, and from step 8 the cores, side-effects and spelling/links tables. Version **0.5.0** is committed as `e3479b82`.
- **Skipped at Larry's instruction** (Amendment 1, sent in chat during Part 3):
  - step 5, `R CMD check --as-cran`;
  - step 6, the full suite;
  - the timing table of step 8.
  Neither run had started, so nothing had to be stopped. `cran-comments.md` says the release check has not been run yet.
- **Ready to submit?** Not as far as this task can tell:
  - the release check has not been run;
  - two CRAN-policy items need attention (§4.1): a CRAN-run test starts up to 6 workers, and `mrct_region_sims()` changes an option without restoring it.

## 1. Provenance and commits

- Run from `~/Documents/GitHub/forestsearch` on branch `cran-prep`.
- **Step 0:** `git status --porcelain` was empty at the start. `cran-prep` = `master` = `2beb79c4`. The step's gate diff was empty and the part 5 report was present: gate passed.
- **Start SHA `9657afd4`** (the TASK commit). `DESCRIPTION` was 0.4.1.9011, as expected.
- **Start export:** package paths, written to `/tmp/cran6/start`.
- **Scratch:** under `/tmp/cran6`. Nothing was installed. Every reproducer used `pkgload::load_all()`, a sequential plan and `timeout 120`.
- **Siblings, at the start and at the end** (all three values unchanged; checked by the final gates, §5):

  | Repository | Branch | HEAD | `git status --porcelain` |
  |---|---|---|---|
  | fs-post-selection | `main` | `5d165acc` | empty |
  | fs-glms-interpretable | `main` | `dc94a17e` | empty |

- **Tools:** R 4.6.2, Pop!_OS 24.04, devtools 2.5.2, roxygen2 8.1.1, testthat 3.3.2.

| SHA | Commit |
|---|---|
| `9657afd4` | dev: add CRAN preparation task, part 6 (start) |
| `9787bbb4` | fix: run_simulation_analysis() reports RD, IRD and RR on their own scales |
| `4eb26d90` | dev: task part 6, amendment 1 (skip release check, suite and timing table) |
| `e5fad06b` | docs: NEWS bullet for run_simulation_analysis() effect scales (part 6, Part 1) |
| `e3479b82` | chore: release 0.5.0 |
| (this commit) | dev: CRAN preparation report, part 6 |

There is no `chore: release 0.5.0, check repairs` commit, because step 5 did not run.

## 2. Part 1 — simulation effect scales (`R/run_simulation_analysis.R`)

**At the start SHA** (`9657afd4`):

- FS extractor:
  - the family switch is at `:995-1001`: `"OR"` and `"RD"` give `binomial()`, `"IRR"` and `"IRD"` give `poisson()`, `"MD"` gives `gaussian()`, and the default (which includes `"RR"`) gives `binomial()`;
  - `log_scale <- effect_measure %in% c("OR", "IRR", "RR")` is at `:1003`;
  - the return is at `:1019-1020`: `b <- coef(fit)[[2]]; if (log_scale) exp(b) else b`.
- GRF extractor: the same code at `:1392-1399` and `:1414-1415`.
- What this returned: `"RD"` gave the log odds ratio, `"IRD"` the log rate ratio, and `"RR"` the odds ratio.
- Truth: `compute_aor()` returns `m1 - m0` for RD/IRD/MD and `m1 / m0` for RR/IRR (`R/oc_analyses.R:307-310`).
- Exposure: the rate ratio uses the offset `log(pmax(<offset.name>, 1e-6))`, but only when that column is in the set's data (FS `:1011-1014`, GRF `:1406-1409`).

**At head:**

- New helper `.glm_arm_effect()` (`:990`). From the same fit it computes `v0 = linkinv(b0)` and `v1 = linkinv(b0 + b1)` and returns:
  - RR: `v1 / v0`;
  - RD and IRD: `v1 - v0`;
  - a value that is not finite: `NA_real_`.
- `arm_scale <- effect_measure %in% c("RD","RR","IRD")` is tested before `log_scale`: FS `:1048` and `:1065`, GRF `:1450` and `:1466`.
- With the offset, the arm rates are per unit of exposure, so `r1 / r0` = the `"IRR"` value.
- OR, IRR, MD, survival and the NULL-measure path are unchanged.
- The `@return` text of `run_simulation_analysis()` now states the scale of every measure.
- New test file: `tests/testthat/test-sim-effect-scales.R`.

**Consumers** of `hr.itt`, `hr.H.hat`, `hr.Hc.hat`, `hr.H.true` and `hr.Hc.true` (start SHA). One research agent and one verifier found 7. None applies `log`/`exp` to these columns or compensates for the old scale, so **none is wrong for the corrected value**:

| Consumer | What it does with RD / IRD / RR |
|---|---|
| `invert_hr_columns()` `R/simulation_tables.R:131-144` | Benefit notation only: RD and IRD become `-x`, RR becomes `1/x`. |
| `build_classification_table()` `:272-463` | Inverts the columns but never shows them. |
| `summarize_simulation_results()` / `summarize_single_analysis()` `R/oc_analyses.R:428-557` | Plain means. |
| `format_oc_results()` `:628-1051` | Means. The spanner label comes from the measure (`"RD"`, `"IRD"`, `"GLM RR"`). |
| `build_estimation_table()` `R/simulation_tables.R:584-1029` | Relative bias against the DGM truth (`get_dgm_hr()`), which is on the RD/RR/IRD scale. |
| `interpret_estimation_table()` `:1095-1481` | Same as above; the narrative names the measure. |
| `plot_effect_distribution()` `R/plot_effect_distribution.R:139-404` | RD: null line at 0, flip `-x`. RR and IRD are rejected by `match.arg`, before the fix and after it. |

Before the fix these displays put "RD", "IRD" and "GLM RR" labels on a log OR, a log IRR and an OR.

**Red / green** (in the repository):

- With the test file only: 40 of 54 expectations fail (all RD, RR and IRD checks). The 10 IRR checks and 4 one-arm NA checks pass. 3.45 s.
- With the fix: 54/54 pass, 2.85 s.
- Hand working (cell totals are in the test's comments). Binary, ITT: p1 = 30/80, p0 = 25/80, so RD = 0.0625 and RR = 1.2. Count with offset, ITT: r1 = 120/160, r0 = 72/120, so IRD = 0.15 and IRR = 1.25.
- The fitted values agree with the hand values to ≤ 8.6e-10.

**Unaffected case:**

- OR, IRR (with and without an offset), MD, survival and NULL-measure extractor values for both extractors, plus `run_simulation_analysis()` with 3 replicates per outcome type (`seed_base` 8316951):
  - 26/26 `identical()` to the start export (only `fs_minutes` dropped);
  - output: `REPORT_cran_prep_release_2026-11-12_p1_unaffected.txt`.
- The existing extractor tests (5 files): 274/274 pass.
- My survival replicates used the default `confounders_base` (`v1..v7`). With those, `forestsearch()` stops in `get_FSdata`, in both trees, so those FS rows are failure paths.
- An independent verifier's own cases fill that gap, with found survival rows and found GRF count rows:
  - extractors 86/86 identical;
  - `run_simulation_analysis()` 21/21 identical;
  - `methods=` path 4/4 identical.

**Before/after example:**

- Setup: binary RD DGM, replicate 1, n = 600, GRF row.
- `hr.itt`: 0.269664 (the log OR of the drawn trial) → 0.066256 (= p1 − p0 = 0.595668 − 0.529412).
- `hr.H.true`: 1.36727 → 0.265392, against the `ahr.H.true` truth of 0.323019.

**Checks done before the amendment:**

- Two independent verifiers (one on correctness, one on the unaffected case and scope) and one comment check.
- They found no code fault.
- They found three false comment/roxygen claims, all corrected before the commit (the R code is identical when parsed):
  1. The IRD scale with an offset (see OPEN ITEMS).
  2. The claim that an empty arm gives a large RR.
  3. A stale `log_scale` comment.
- `devtools::document()` in the repository made no further change. `checkRd` is clean. The file still ends without a newline, as at the start.

## 3. Part 2 — runs with RD, IRD or RR (read-only)

| Repo | File:line | Measure | Current or archive | Outputs | Tracked payloads |
|---|---|---|---|---|---|
| — | none | — | — | — | 0 |

**Totals: 0 qualifying sites, 0 tracked payloads.** By scope:

- **forestsearch `quarto/`:** 311 `run_simulation_analysis()` run sites.
  - 62 are on the survival path, including 4 `glm_dgm` runs with no `outcome_type`.
  - 249 are GLM runs: OR 233 (45 current, 188 archive), MD 4, IRR 10, and 2 mixed scenario-validation wrappers (OR/MD/IRR/survival).
  - No GLM run leaves the measure unset.
  - Wrappers: `run_scenario_sims()`, `fit_and_score()`, `one_rep()`, plus inline `%dofuture%` loops.
  - The measure is always set by `fs_params$effect_measure`.
- **forestsearch `vignettes/`:** 8 sites, all in `articles/paper_simulations.qmd`, all survival (GBSG). `vignettes/articles` is build-ignored.
- **forestsearch `tests/`:** 0 sites.
- **forestsearch `R/`:** no code calls `run_simulation_analysis()`.
- **fs-post-selection:** 0 sites. The only mention is a comment in `analysis/simulations/actg175/binary/_template.qmd:81`. Its intended bundle is OR.
- **fs-glms-interpretable:** 0 sites in tracked files or in history. Git-ignored dumps of the 0.3.5.9000 sources mention the name; they are not runs.
- How this was checked: each scope was searched by one agent and verified by a second.
- The completeness critic was stopped under Amendment 1 before it finished, so it adds nothing.

## 4. Part 3 — release preparation

**Version (step 1):**

- `DESCRIPTION`: `Version: 0.5.0`. There is no `Date` field.
- Other version strings, listed and not edited:
  - `cran-comments.md`: rewritten in step 3; it now says 0.4.0 → 0.5.0.
  - `_pkgdown.yml:433, :435`: "Additional functions (0.4.0)".
  - `README.md`: no version.
  - `inst/CITATION`: does not exist.
  - Help text and printed messages say "Deprecated (0.4.0.9000)" or "since 0.4.0.9000":
    - `R/forestsearch_main.R:536, :937, :968` → `man/forestsearch.Rd:308, :745, :785`;
    - `R/subgroup_search.R:39`;
    - `R/plot_subgroup_results_forestplot.R:42`;
    - `R/fs_family_report.R:269, :354`.
  - 0.4.0.9000 was never released (see OPEN ITEMS).

**NEWS (step 2): gate PASS** (`REPORT_cran_prep_release_2026-11-12_news_gate.txt`; scripts `_news_restructure.R` and `_news_gate.R`).

- 80 bullets before (`e5fad06b`) and 80 after. The multisets of bullets and of lines (627) are identical, and no bullet was merged.
- Everything from `# forestsearch 0.4.0` on is byte-identical.
- There is one `# forestsearch 0.5.0` heading, with five groups:

  | Group | Bullets |
  |---|---|
  | Breaking changes | 5 |
  | Bug fixes | 34 |
  | New features | 20 |
  | Superseded and deprecated | 8 |
  | Documentation | 13 |

- The opening summary is 8 lines.
- The group rules are in the script header. Development "Behaviour changes" go to Bug fixes, except:
  - three go to Breaking changes: the bootstrap orientation, the caller's random-number state, and `fpr_calibration()`'s default and new error;
  - `mr_estimates_table()` without a template goes to New features.
- Other moves:
  - "Removed" and `survival (>= 3.5-8)` → Breaking changes;
  - "Performance" and the lean-fit helper → New features;
  - "Internal" test items → Documentation;
  - the internal `fit_causal_forest_glm()` removal → with the deprecations.
- The old development paragraph is replaced by the summary.
- The summary and the Part 1 bullet were checked by independent agents before the amendment and corrected twice.

**`cran-comments.md` (step 3)**, full text:

```
## Submission

This is an update from forestsearch 0.4.0 (the current CRAN version) to
0.5.0. NEWS.md lists the changes under "forestsearch 0.5.0".

This release removes `create_gbsg_dgm()` and `simulate_from_gbsg_dgm()`,
which were deprecated in 0.4.0 (the first CRAN release to export them,
where they already warned that they were deprecated).

## Test environments

None yet for 0.5.0. The local release check (`R CMD build`, then
`R CMD check --as-cran` on the tarball), win-builder and R-hub have not
been run for this version.

## R CMD check results

The release check has not been run yet.

## Reverse dependencies

There are no reverse dependencies: CRAN's package index (read on
2026-11-12) lists no package that depends on, imports, links to or suggests
forestsearch.
```

- `cran-comments.md` is excluded from the build (`.Rbuildignore:15`).
- The parenthesis after "deprecated in 0.4.0" is evidence from CRAN's own sources:
  - In 0.1.0 both functions call `.Deprecated()` but are not exported (`R/sim_aft_gbsg.R:191, :771`; NAMESPACE exports only `setup_gbsg_dgm`).
  - 0.4.0 exports them, still with `.Deprecated()` (`:211, :800`).

**Reverse dependencies (step 4):**

- CRAN's index (25,722 packages, read at 20:33 PST) lists forestsearch 0.4.0.
- `tools::package_dependencies(reverse = TRUE)` over Depends, Imports, LinkingTo, Suggests (and Enhances) returns **none**, with the index filtered and unfiltered. A text search of those fields also finds none.
- So there was nothing to download or search.

**Steps 5 and 6: NOT RUN** (Amendment 1). There are no NOTEs, timings, tarball size or suite counts to report.

**Default-path gate (step 7): PASS.**

- Part 4's six cases (`REPORT_cran_prep_code_fixes_2026-11-12_gate_run.R` and `_gate_compare.R`, unchanged), run on the start export and on a `git archive` of `e3479b82`.
- 6/6 are `identical()` after dropping timing fields, and `.Random.seed` is restored.
- Wall: 18.3 s (start) and 18.0 s (head). Output: `REPORT_cran_prep_release_2026-11-12_gate_compare.txt`.

### 4.1 CRAN policy audit (step 8; read-only, nothing edited)

How this was done:

- Each table comes from one research agent and one verifier, run on the start export.
- `R/` and the tests changed afterwards only by Part 1, whose new test file is sequential and writes nothing.
- `R CMD check --as-cran` sets `_R_CHECK_LIMIT_CORES_=TRUE` and runs `\donttest` examples.

**Cores:**

| Where | Workers | Runs on CRAN | Status |
|---|---|---|---|
| `tests/testthat/test-mrct-case-study.R:42-44` (`cs_runs()`, first used at `:94-96`) | `min(6, max(2, detectCores()))`: 3–6 workers whenever `detectCores() >= 3`; `_R_CHECK_LIMIT_CORES_` does not cap it | yes, when Suggests `speff2trial` and `latentcor` are installed (`mrct_skip_if_no_build()`) | **NEEDS ATTENTION** |
| `tests/testthat/helper-mrct-region-sims-record.R:26` (record and MR tests) | 2 | yes | compliant |
| `man/reset_workers.Rd:70-76` (`\donttest`) | 2 | yes | compliant |
| grf forests (`R/grf_helpers.R:84-87`, `R/grf_subg_harm_glm.R:469-475`; example `fit_causal_forest.Rd`, the vignette's `fit-grf` chunk, CRAN-run grf tests) | 2 threads under `_R_CHECK_LIMIT_CORES_` (`R/utils_check_threads.R:22-27`); all threads otherwise | yes | compliant under CRAN's settings |
| data.table (Imports; no `setDTthreads` anywhere) | `getDTthreads()` = 50% of logical cores (64 here), even with `_R_CHECK_LIMIT_CORES_`; 2 with `OMP_THREAD_LIMIT=2` | yes, implicitly | unsure: depends on the check machine's `OMP_THREAD_LIMIT` |
| `test-plan-restored.R:153-173` (`n_workers = 3`) | 3 | no (`skip_on_cran()`) | — |

- Uncapped package defaults: `forestsearch_Kfold` workers 6, `mrct_region_sims` `detectCores()-1`, `setup_parallel_SGcons` 4, `fs_fdr_report`, and others.
- No example, CRAN-run test or the vignette reaches them. The functions' own examples are `\dontrun`, and the vignette passes `plan = "sequential"`.

**Side effects:**

| Where | Effect | Runs on CRAN | Status |
|---|---|---|---|
| `R/mrct_simulation.R:414` `progressr::handlers("progress")` in `mrct_region_sims()` | sets the option `progressr.handlers`, not restored | yes: `\donttest` example `man/mrct_core_summary.Rd:74`, and the tests via `helper-mrct-region-sims-record.R:28` and `test-mrct-case-study.R:36` | **NEEDS ATTENTION** |
| `test-forestplot-xlog.R:36` and `test-no-subgroup-dina-grf-plot.R:95` | `forestploter::forest()` with no device open writes `Rplots.pdf` into `tests/testthat` | no (`skip_on_cran()`) | — |
| `R/subgroup_consistency_helpers.R:980` `par(mfrow = c(1, 2))` | not restored | no (only with `plot.sg = TRUE`) | — |
| `\dontrun` examples `save_forestplot.Rd:43`, `extract_selected_subgroups.Rd:70`, `summarize_bootstrap_subgroups.Rd:82` | write to the working directory | no | — |

- No package is installed anywhere.
- No CRAN-run test or example writes outside `tempdir()`.
- No CRAN-run test or example calls `options()` or `par()` directly.

**Spelling and links:**

- The `spelling` package is not installed, so base R's `utils::aspell` (GNU Aspell 0.60.8.1, en) was used instead.
- Distinct flagged words: Rd 343, vignette 81, NEWS 156, DESCRIPTION 8.
- They are technical terms, names, British spellings (behaviour, centred) and tokens split from code spans.
- "DIfference in NAtural" (`R/dina.R:4, :29`) spells out the DINA acronym on purpose.
- No genuine misspelling was found, so nothing was fixed. The full list is in `REPORT_cran_prep_release_2026-11-12_spelling.txt`.
- `urlchecker::url_check()` (8.2 s) found 3 moved URLs, all in the build-ignored `vignettes/articles/references.bib:901, :926, :943, :1143`. `vignettes/` is read-only, so they are listed only.
- **Time table:** skipped (Amendment 1).

## Larry's remaining steps, in order

1. Run the release check: `R CMD build`, then `R CMD check --as-cran`. Fill in `cran-comments.md` (test environments, result, NOTEs). Then run win-builder (release and devel) and R-hub.
2. Decide on the two NEEDS ATTENTION items in §4.1 before submitting.
3. Submit to CRAN.
4. Tag `v0.5.0`.
5. Merge `cran-prep` into `master`.
6. Rebuild the pkgdown site.

## OPEN ITEMS (left open by this task)

- **IRD with an offset: estimate and truth are in different units.**
  - Step 3 sets the IRD estimate per unit of exposure, the convention of the rate ratio.
  - The DGM truth (`.compute_glm_effects()`, `R/generate_glm_dgm.R:676`) and `compute_aor()` give the IRD as a difference in expected counts per patient (`mu` includes the offset).
  - With an offset, `build_estimation_table()`/`interpret_estimation_table()` bias compares unlike units. The verifiers measured the truth/estimate ratio at the mean exposure (2.25, 3.09).
  - This is not a BLOCK under step 2. No consumer was right before, and none compensates for the old scale. The mismatch comes from the convention the task specifies.
  - No committed simulation uses IRD (Part 2). The help page and the NEWS bullet state it. Larry to decide which convention.
- **Pre-existing, unchanged (code audit A-138, the other half):** with `fs_params$outcome_type` set and `effect_measure` NULL, every GLM `hr.*` is NA (`if (logical(0))` inside `tryCatch`), while `ahr.*` default to OR.
- **Pre-existing:** `.run_grf_analysis_gen()` passes `effect_measure = "IRR"` to `grf.subg.harm.glm()`, which accepts only `"log_IRR"`/`"log_RR"`, so the GRF row of a count simulation fails (seen in every count replicate the verifiers ran; `R/run_simulation_analysis.R:1294` at the start, `:1340` at head).
- **Pre-existing:** `plot_effect_distribution()` rejects RR and IRD and defaults to `"OR"`.
- **NEWS (verbatim rule):**
  - The 0.4.1 Documentation bullet says `mr_estimates_table()` "needs" a bootstrap template, while the development New features bullet says it no longer does.
  - The replaced development paragraph was the only NEWS mention of the MRCT consistency article(s).
- **Version strings:** help pages and two printed messages name the unreleased 0.4.0.9000 (§4).
- **CRAN policy:** the two NEEDS ATTENTION items (§4.1).

## Departures

- **Amendment 1:** steps 5, 6 and step 8's time table were skipped. The amendment was appended to the TASK file as its own commit (`4eb26d90`).
- **NEWS bullet for Part 1:** the task does not ask for one. It was added as its own commit (`e5fad06b`) *before* the regrouping, so the step-2 gate compares against it. The reason is that the release changes user-visible simulation results.
- **`cran-comments.md`:** the task's sentence "deprecated in 0.4.0" is kept, with a parenthesis reconciling it with the NEWS bullet's "since before 0.1.0".
- **Part 2 critic:** stopped under Amendment 1 before it returned.

## Final gates

Run after the report commit with `bash dev/reports/REPORT_cran_prep_release_2026-11-12_final_gates.sh 9657afd4` (status of the tree, siblings unchanged, files changed since the start, nothing pushed, merged or tagged). The result is in the closing message of the session, not in this file, so that the file does not change after the gates run.
