# REPORT — refs, `max.minutes` deprecation, untested-option notes, manuscript-aligned examples

**Date:** 2026-10-27. **Task doc:** `dev/tasks/cc_task_docs_refs_maxminutes_2026-10-27.md`.
**Branch:** `feature/docs-refs-maxminutes`, from `master` at `b574c618f0032f47a58774283d71cda5ca78bd7d` (verified; tree was clean).
Not pushed, merged, fetched or pulled. No R CMD check, no vignette build, no full test suite.

**Outcome:** §1.1–§1.6 applied; every gate in §2.6, §3 and §4 passed. There is one approved deviation from §1.2 (see "Deviation from §1.2" below).

## Deviation from §1.2 (approved by Larry mid-task, 2026-10-27)

The §1.2 mechanism as written did not achieve what §1.2 states:

1. `forestsearch()` builds `args_call_all <- mget(names(formals()), envir = environment())` (`R/forestsearch_main.R` L1755 at `b574c618`), so `args_call_all$max.minutes` is **always** present (default 3). The `subgroup.search()` call is assembled as `search_args <- modifyList(args_call_all, search_overrides)` and then filtered to `names(formals(subgroup.search))` (L≈3350). Removing the `max.minutes = max.minutes,` override line therefore still passes `max.minutes = 3` explicitly, and `subgroup.search()`'s `!missing()` warning would fire on **every** `forestsearch()` call.
2. Bootstrap (`do.call(forestsearch, args_FS_boot)` with `args_FS_template <- fs.est$args_call_all`, `R/bootstrap_analysis_dofuture.R` L406/L558/L614), CV (`R/forestsearch_cross_validation.R` L479, L1001), `fs_fdr_report()` (`R/fs_fdr_report.R` L362), `mrct_simulation` (`R/mrct_simulation.R` L471), `calibrate_null_correction` (`R/calibrate_null_correction.R` L355), and also `compare_selection_rules`, `fpr_calibration` and `run_simulation_analysis` all replay a full argument list via `do.call(forestsearch, ...)`. A bare `!missing(max.minutes)` trigger in `forestsearch()` would therefore warn once in every replicate or fold, even when the user never passed the argument.

Larry chose option 1:
- `forestsearch()` uses a value trigger, the same as the existing `max_n_confounders` pattern: `if (!missing(max.minutes) && !identical(as.numeric(max.minutes), 3)) warning(...)`.
- `search_args$max.minutes <- NULL` is set immediately before `do.call(subgroup.search, search_args)`, in addition to removing the override line.
- `subgroup.search()` keeps the `!missing()` trigger. It is reached with `max.minutes` only when a user calls it directly.
- The `forestsearch()` `@param max.minutes` text says "A non-default value emits one warning" where §1.2 had "Passing it emits one warning". The `subgroup.search()` text keeps "Passing it emits one warning". The NEWS bullet matches both.
- Residual behaviour, the same as `max_n_confounders` today: a user who passes a non-default `max.minutes` to `forestsearch()` gets the warning again in each replayed bootstrap/CV call. CV wraps its calls in `suppressWarnings()`; bootstrap does not.

## SHA ladder

```
a1b07bc9 dev: task doc — refs, max.minutes deprecation, untested-option notes, examples (2026-10-27)
8758718a docs: cite León & Anderson (2026a, 2026b) in forestsearch(), DESCRIPTION, README
f35399a5 deprecate: max.minutes in forestsearch() and subgroup.search(); drop call-site uses
a2d1c70a docs: untested-option notes for adjust_covariates and ps_* arguments
9fd492d4 docs: manuscript-aligned examples for forestsearch() (GBSG survival, ACTG175 binary/continuous)
2dfeea63 DESCRIPTION: Suggests speff2trial
885e0be2 docs: regenerate man pages; test-max-minutes-deprecation
(report commit: the commit that adds this file — `dev: report — refs, max.minutes deprecation, examples (2026-10-27)`)
```

Note: the NEWS bullet (§1.2 item 4) is in the `deprecate:` commit (f35399a5), not in the man/test commit as §2.8 lists it, because §1.2 item 4 is part of that change.

## §1.1 item 4 — other pages citing 2609.38361 (per-file decision)

`grep -rn "2609.38361" R/` at `b574c618` returned only `R/forestsearch_methods.R` (L138, a printed caveat string, not Rd; L225, roxygen for `print.forestsearch`). `grep -rl 2609.38361 man/` returned `forestsearch-package.Rd` and `print.forestsearch.Rd`. **`fs_mr_inference` and `fs_betaHhat_one` do not cite 2609.38361**, contrary to the 2026-10-26 handoff, so neither was changed.

| Page | Decision | Change |
|---|---|---|
| `man/forestsearch-package.Rd` | **LA-2 added.** The page is the package-level general description. | No R source edit: the page is built from `DESCRIPTION`'s `Description:` field (`R/forestsearch-package.R` is `"_PACKAGE"`), so the §1.1 item 3 sentence carries through: "The search methodology is given in León and Anderson (2026) \doi{10.48550/arXiv.2610.28379}." |
| `man/print.forestsearch.Rd` | **LA-1 only.** The page covers the MR/IJ post-selection block (field bounds, Bonferroni pair, IJ interval, p-hat). | none |
| `R/forestsearch_methods.R` L138 `.fs_mr_caveats()` | Not an Rd page; printed output for the MR block. | none |

`R/fs_mr_inference.R` L315 and `R/betaHhat_truth.R` L212 cite "León and Anderson (2026)" without the arXiv number. They are MR/IJ-specific and were left unchanged; they now read as unsuffixed against the 2026a/2026b convention (flagged under "Surprises").

## §1.2 item 3 / §1.4 — `grep -rn "max.minutes" R/ tests/ vignettes/ inst/`

### Before (`b574c618`, tracked files, `git grep`)
```
R/forestsearch_main.R:884:#' @param max.minutes Numeric. \strong{Currently inert; scheduled for
R/forestsearch_main.R:1560:                         max.minutes = 3,
R/forestsearch_main.R:3313:        max.minutes = max.minutes,
R/fs_family_report.R:164:  max.minutes <- .arg("max.minutes")
R/fs_family_report.R:353:  add("time cap", "max.minutes", sprintf("max.minutes = %s", fmt(max.minutes)), INE,
R/run_simulation_analysis.R:63:    # Inert: no path consults max.minutes (see ?forestsearch).  Kept so this
R/run_simulation_analysis.R:66:    max.minutes              = 5,
R/subgroup_search.R:39:#' @param max.minutes Numeric. \strong{Currently inert; scheduled for
R/subgroup_search.R:98:                            hr.threshold = 1.0, max.minutes = 30,
R/subgroup_search.R:171:    max.minutes = max.minutes, t.start = t.start,
R/subgroup_search.R:264:                                         minp, rmin, max.minutes, t.start,
vignettes/articles/methodology.Rmd:995:  max.minutes                = 3,                            # default 3
vignettes/articles/methodology.Rmd:1106:  max.minutes                = 3,                            # default 3
vignettes/forestsearch.Rmd:342:  max.minutes                = 3,                            # default 3; per-stage timeout
vignettes/survival-post-selection.qmd:82:  fs.splits = 200, max.minutes = 3,
```
`tests/` and `inst/` had no hits. The untracked, git-ignored `vignettes/survival-post-selection.html` (`.gitignore:15 vignettes/*.html`) also matches at L2302. It is a rendered output and was not edited (no vignette build).

### After (HEAD)
```
R/subgroup_search.R:39:#' @param max.minutes \strong{Deprecated} (0.4.0.9000): accepted and
R/subgroup_search.R:100:                            hr.threshold = 1.0, max.minutes = 30,
R/subgroup_search.R:112:  if (!missing(max.minutes)) {
R/subgroup_search.R:113:    warning("subgroup.search(): 'max.minutes' is deprecated and ignored. ",
R/subgroup_search.R:181:    max.minutes = max.minutes, t.start = t.start,
R/subgroup_search.R:274:                                         minp, rmin, max.minutes, t.start,
R/forestsearch_main.R:891:#' @param max.minutes \strong{Deprecated} (0.4.0.9000): accepted and
R/forestsearch_main.R:1673:                         max.minutes = 3,
R/forestsearch_main.R:1721:  # max.minutes are accepted so existing calls do not break, warned about
R/forestsearch_main.R:1725:  # model-based front end.  max.minutes was never enforced.  All three
R/forestsearch_main.R:1744:  # max.minutes, and bootstrap/CV replay it via do.call(forestsearch, ...).
R/forestsearch_main.R:1745:  if (!missing(max.minutes) &&
R/forestsearch_main.R:1746:      !identical(as.numeric(max.minutes), 3)) {
R/forestsearch_main.R:1747:    warning("forestsearch(): 'max.minutes' is deprecated and ignored. The ",
R/forestsearch_main.R:3475:  # max.minutes is deprecated and inert; args_call_all always carries it, so
R/forestsearch_main.R:3477:  search_args$max.minutes <- NULL
R/fs_family_report.R:164:  max.minutes <- .arg("max.minutes")
R/fs_family_report.R:353:  add("time cap", "max.minutes", sprintf("max.minutes = %s", fmt(max.minutes)), INE,
tests/testthat/test-max-minutes-deprecation.R:2:# max.minutes deprecation (0.4.0.9000): the formal is kept in forestsearch()
tests/testthat/test-max-minutes-deprecation.R:4:# warns on a non-default value (args_call_all always carries max.minutes and
tests/testthat/test-max-minutes-deprecation.R:52:.mm_passed  <- .mm_capture(.mm_run(max.minutes = 5))
tests/testthat/test-max-minutes-deprecation.R:54:test_that("a default forestsearch() call emits no max.minutes warning", {
tests/testthat/test-max-minutes-deprecation.R:55:  expect_false(any(grepl("max.minutes", .mm_default$warnings, fixed = TRUE)))
tests/testthat/test-max-minutes-deprecation.R:59:test_that("max.minutes = 5 emits exactly one deprecation warning", {
tests/testthat/test-max-minutes-deprecation.R:60:  hits <- grep("'max.minutes' is deprecated", .mm_passed$warnings, fixed = TRUE)
tests/testthat/test-max-minutes-deprecation.R:67:test_that("max.minutes is inert: the fit is identical apart from the recorded argument", {
tests/testthat/test-max-minutes-deprecation.R:70:  expect_identical(b$args_call_all$max.minutes, 5)
tests/testthat/test-max-minutes-deprecation.R:71:  a$args_call_all$max.minutes <- NULL
tests/testthat/test-max-minutes-deprecation.R:72:  b$args_call_all$max.minutes <- NULL
tests/testthat/test-max-minutes-deprecation.R:76:test_that("subgroup.search() warns once with its own prefix when passed max.minutes", {
tests/testthat/test-max-minutes-deprecation.R:86:  passed <- .mm_capture(suppressMessages(call_ss(max.minutes = 5)))
tests/testthat/test-max-minutes-deprecation.R:87:  hits <- grep("'max.minutes' is deprecated", passed$warnings, fixed = TRUE)
tests/testthat/test-max-minutes-deprecation.R:91:  expect_false(any(grepl("max.minutes", dflt$warnings, fixed = TRUE)))
tests/testthat/test-max-minutes-deprecation.R:94:test_that("the max.minutes formals are kept through this release", {
tests/testthat/test-max-minutes-deprecation.R:95:  expect_true("max.minutes" %in% names(formals(forestsearch)))
tests/testthat/test-max-minutes-deprecation.R:96:  expect_true("max.minutes" %in% names(formals(subgroup.search)))
vignettes/survival-post-selection.html:2302:<span id="cb4-17"><a href="#cb4-17" aria-hidden="true" tabindex="-1"></a>  <span class="at">fs.splits =</span> <span class="dv">200</span>, <span class="at">max.minutes =</span> <span class="dv">3</span>,</span>
```

### §1.4 removals (file : line at `b574c618`)
- `R/run_simulation_analysis.R:63–66`: the 3-line "Inert: ..." comment and `max.minutes = 5,` removed from `default_sim_params()`. That list feeds `do.call(forestsearch, fs_args)` (L823), where 5 would now warn.
- `vignettes/forestsearch.Rmd:342`: `max.minutes = 3,  # default 3; per-stage timeout` removed.
- `vignettes/survival-post-selection.qmd:82`: `fs.splits = 200, max.minutes = 3,` became `fs.splits = 200,`.
- `vignettes/articles/methodology.Rmd:995` and `:1106`: `max.minutes = 3,  # default 3` removed (both).
- `R/forestsearch_main.R:3313`: forwarding line `max.minutes = max.minutes,` removed from `search_overrides`, and `search_args$max.minutes <- NULL` added (see the deviation section).

§4.3 lists `vignettes/*.Rmd`; the `.qmd` and `articles/` files are under `vignettes/` and were edited only for this removal, as §1.4 requires ("every call in `vignettes/`").

### Hits left in place
- `R/subgroup_search.R:181, :274`: internal threading to `search_combinations_parallel()` (§1.2 item 2, out of scope).
- `R/fs_family_report.R:164, :353`: `fs_family_report()` reads `max.minutes` from a fit's `args_call_all` and reports it on an "inert" row. This is not a call site. Its note text, "subgroup.search() forwards it and never compares it", is now stale: forestsearch() no longer forwards it (flagged under "Surprises", not edited).
- **Legacy, not edited (rule §0):** `quarto/applications/gbsg/analysis_gbsg_survival_effMaxSG.qmd:580: max.minutes = 3,` (§0 gives L≈584; it is at L580). Other `.qmd` sources under `quarto/applications/` (excluding `_payloads*`/`_archive`):
```
quarto/applications/gbsg/analysis_gbsg_survival_maxeff_mrconfirm.qmd:188:  use_twostage          = TRUE, max.minutes = 3,
quarto/applications/gbsg/analysis_gbsg_survival_sgfocus.qmd:580:  max.minutes           = 3,
quarto/applications/gbsg/analysis_gbsg_survival_effMaxSG.qmd:580:  max.minutes           = 3,
quarto/applications/gbsg/analysis_gbsg_survival_frozen_family.qmd:302:  use_twostage       = TRUE, max.minutes = 3,
quarto/applications/gbsg/analysis_gbsg_survival_multimethod.qmd:775:  max.minutes               = 3,              # default 3
quarto/applications/gbsg/analysis_gbsg_survival_multimethod.qmd:1341:  max.minutes               = 3,              # default 3
```
  In total `quarto/applications/**` has 60 matching lines in 24 files, including payload and archive HTML. All of them pass the default value 3, so none of them will trigger the new `forestsearch()` warning. There are further hits outside `quarto/applications` (`quarto/GuoHe/*`, `quarto/dina/*`, `quarto/legacy calibration files/*`, `quarto/simulations/*`); none were edited.

## §1.3 — untested-option sentence: files and functions

`forestsearch()` (`R/forestsearch_main.R`): `adjust_covariates`, `ps_method`, `ps_adjust_method`, `ps_hat`, plus the new `@section Options not tested against the manuscript analyses:` placed immediately before `@examples`.

§1.3 item 3: these exported functions document `adjust_covariates`/`ps_*` and received the same first sentence:
| File : line (at `b574c618`) | Function | Param |
|---|---|---|
| `R/subgroup_search.R:65` | `subgroup.search` | adjust_covariates |
| `R/subgroup_consistency_main.R:211` | `subgroup.consistency` | adjust_covariates |
| `R/consistency_resample.R:375` | `consistency_resample` (inherited by `consistency_resample_compare` via `@inheritParams`) | adjust_covariates |
| `R/Cox_estimation_helpers.R:94` | `build_cox_formula` | adjust_covariates |
| `R/glm_effect_estimators.R:63, :65` | `make_effect_estimator` | adjust_covariates, ps_adjust_method |
| `R/guohe_algorithm3.R:165` | `guohe_algorithm3` (inherited by `guohe_adaptive_r`) | adjust_covariates |
| `R/plot_sg_glm_outcomes.R:40` | `plot_sg_glm_outcomes` | adjust_covariates |
| `R/plot_subgroup_results_forestplot.R:79` | `plot_subgroup_results_forestplot` | adjust_covariates |
| `R/subgroup_consistency_helpers.R:168` | `get_split_hr_fast` | adjust_covariates |
| `R/subgroup_consistency_helpers.R:239` | `run_single_consistency_split` | adjust_covariates |
| `R/subgroup_consistency_helpers.R:1323` | `evaluate_subgroup_consistency` | adjust_covariates |
| `R/subgroup_consistency_helpers.R:1610` | `evaluate_consistency_twostage` | adjust_covariates |

Not edited because they are not exported (absent from NAMESPACE `export()`, or `@noRd`): `.consistency_cox_pieces`, `.consistency_glm_pieces`, `.fs_adjust_terms`, `.fs_adjust_vars`, `evaluate_combination_with_status`, `fit_cox_for_subgroup`, `.make_eval_subgroup_consistency`.

## §1.5 — example argument-name check

Script: `tools::Rd2ex("man/forestsearch.Rd", commentDontrun = FALSE)` → `parse()` (7 top-level expressions) → walk every `forestsearch(...)` call → compare against `names(formals(forestsearch))`. 87 argument names across 3 calls (29 + 30 + 28); **not formals: none**. Unique names checked:
`df.analysis, outcome.name, event.name, treat.name, id.name, confounders.name, is.RCT, use_lasso, use_grf, use_dina, conf.cont_jcuts, conf_force, subgroup_method, sg_focus, selection_rule, effect_neighborhood, hr.threshold, hr.consistency, pconsistency.threshold, consistency_method, fs.splits, n.min, d0.min, d1.min, maxk, mr_inference, mr_inference_args, parallel_args, seedit, outcome_type, effect_measure, adverse_outcome, effect.threshold, consistency.threshold`.

- Note 1: formals' defaults are `use_twostage = TRUE`, `minp = 0.025`, `stop_threshold = pconsistency.threshold`, as the task doc expected.
- Note 3: all six Example 3 `mr_inference_args` elements (`ci_method`, `draws`, `include_complement`, `field_complement`, `field_scale_complement`, `return_reselection`) are formals of `fs_mr_inference()` and are each read by `forestsearch()` via `.g_mr(mr_inference_args$<name>, ...)` (`R/forestsearch_main.R` L≈3771–3803 at `a2d1c70a`). **None dropped.**
- `mr_estimates_table` is exported. `speff2trial` is installed locally. `survival::gbsg` has `pid, age, meno, size, grade, nodes, pgr, er, hormon, rfstime, status`, so every column the example uses exists.
- §1.6: there is no `data/` and no `inst/extdata/`, so ACTG175 is not bundled and `speff2trial` was added to `Suggests`.
- The examples were not run (`\dontrun{}`, minutes each; no example run was asked for).

## `git diff b574c618..HEAD -- DESCRIPTION NEWS.md`

```diff
diff --git a/DESCRIPTION b/DESCRIPTION
index 07b066a2..3fd2c1e1 100644
--- a/DESCRIPTION
+++ b/DESCRIPTION
@@ -20,6 +20,8 @@ Description: Implements statistical methods for exploratory subgroup
     <doi:10.1002/sim.10163>.
     Post-selection inference for the identified subgroup follows León and
     Anderson (2026) <doi:10.48550/arXiv.2609.38361>.
+    The search methodology is given in León and Anderson (2026)
+    <doi:10.48550/arXiv.2610.28379>.
 License: MIT + file LICENSE
 Depends: 
     R (>= 4.1.0)
@@ -66,6 +68,7 @@ Suggests:
     katex,
     latentcor,
     quarto,
+    speff2trial,
     testthat (>= 3.0.0)
 URL: https://github.com/larry-leon/forestsearch, https://larry-leon.github.io/forestsearch/
 BugReports: https://github.com/larry-leon/forestsearch/issues
diff --git a/NEWS.md b/NEWS.md
index e10dacdd..cd1a0b17 100644
--- a/NEWS.md
+++ b/NEWS.md
@@ -9,6 +9,12 @@
   cut columns is removed. The candidate family is now the cut grid alone,
   and no model is fitted in the forest-search path unless `use_lasso`,
   `use_grf` or `use_dina` is `TRUE`.
+* `max.minutes` in `forestsearch()` and `subgroup.search()` is deprecated.
+  It was documented as inert ("scheduled for deprecation in v0.3.0") and
+  that was never actioned. Both formals are kept (same positions, same
+  defaults); the argument has no effect. `forestsearch()` emits one
+  deprecation warning for a non-default value, and `subgroup.search()` one
+  warning whenever it is passed. The formals are removed in a later release.
 * What changes for callers who passed a non-NULL `vi.grf.min`: at values
   `<= 0` the screen never filtered anything, so the selected subgroup's
   membership is unchanged. The clause order of the `sg.harm` label (and the
```
`Config/roxygen2/version` was re-stamped 8.1.1 by `devtools::document()` and restored to 8.1.0 (not in the diff). `devtools::install()` did not re-stamp it.

## §2.6 gate
`NAMESPACE` unchanged (`git diff --quiet b574c618 -- NAMESPACE`). The `DESCRIPTION` diff is exactly the two fields above.

## §3 and §4.2 — test output

Package installed with `devtools::install(upgrade = FALSE)` (version 0.4.0.9000); the installed `forestsearch()` body contains the new warning.

| File | load_package | Result | Wall-clock |
|---|---|---|---|
| `tests/testthat/test-max-minutes-deprecation.R` (new) | source | `[ FAIL 0 \| WARN 0 \| SKIP 0 \| PASS 12 ]` | 4 s |
| `tests/testthat/test-vi-deprecation.R` | source | `[ FAIL 0 \| WARN 0 \| SKIP 0 \| PASS 15 ]` | 5 s |
| `tests/testthat/test-lean-fits.R` | installed | `[ FAIL 0 \| WARN 0 \| SKIP 0 \| PASS 344 ]` | 3 s |

Command: `NOT_CRAN=true Rscript -e 'testthat::test_file("<file>", package = "forestsearch", load_package = "<source|installed>")'`.

How the test file meets §3:
1. A default call (continuous fixture `.make_continuous_data(N = 300, seed = 42)`, the same as `test-vi-deprecation.R`) emits no warning containing `max.minutes` and none prefixed `subgroup.search():`.
2. `max.minutes = 5` emits exactly one `'max.minutes' is deprecated` warning, prefixed `forestsearch(): `. All warnings are captured with `withCallingHandlers`, and the test counts the hits, so a second warning fails it. There is no `subgroup.search():` warning, which confirms it is not forwarded.
3. Identity: `args_call_all$max.minutes` is set to NULL on both objects, the timing fields (`timing_seconds`, `minutes_all`, `time_search`) are stripped with the `test-declaration-c0.R` stripper, and the `saveRDS`/`readRDS` copies are compared with `expect_identical`. The fixture finds a subgroup (`{age <= 47}` & `{biomarker_hi}`), so the comparison is over a non-trivial fit.
4. A direct `subgroup.search(max.minutes = 5)` emits exactly one warning, prefixed `subgroup.search(): `; a direct default call emits none.
5. Both formals are still present.

`test-vi-deprecation.R` does not use `skip_on_cran()`, so the new file does not either.

## §4 gates
- **4.1** `tools::checkRd()` on `forestsearch.Rd` and `subgroup.search.Rd`: no errors or warnings. It does print level −1 notes, "Non-ASCII contents without declared encoding", for every `León` line (13 in `forestsearch.Rd`, 1 in `subgroup.search.Rd`). The untouched `print.forestsearch.Rd` at `b574c618` already gives the same note, because the pages carry no `\encoding{}` and the standalone call does not see `DESCRIPTION`'s `Encoding: UTF-8`. With `checkRd(..., encoding = "UTF-8", def_enc = TRUE)`, which is how R CMD check treats a package that declares an encoding, both pages report **0 problems**. I treated the gate as passed; Larry should confirm that reading. `parse(Rd2ex(...))` succeeds (7 expressions). The argument-name check passes (above). The `\description` has both `arXiv.2609.38361` and `arXiv.2610.28379`. `R CMD Rd2txt man/forestsearch.Rd | grep -c "León and Anderson"` = **8** (needed ≥ 3).
- **4.2** Green (above).
- **4.3** `git diff --stat b574c618..HEAD` touches only the permitted set: see the stat below. `NAMESPACE`, `quarto/`, `data/` and `inst/` are untouched.
- **4.4** `git status --porcelain quarto/` is empty.
- **4.5** No added line in `R/` or `DESCRIPTION` is over 80 characters (checked on the `+` lines of `git diff -U0 b574c618..HEAD -- R/ DESCRIPTION`).

```
 DESCRIPTION                                        |   3 +
 NEWS.md                                            |   6 +
 R/Cox_estimation_helpers.R                         |   6 +-
 R/consistency_resample.R                           |   6 +-
 R/forestsearch_main.R                              | 272 +++++++++++++------
 R/glm_effect_estimators.R                          |  12 +-
 R/guohe_algorithm3.R                               |   6 +-
 R/plot_sg_glm_outcomes.R                           |   6 +-
 R/plot_subgroup_results_forestplot.R               |   6 +-
 R/run_simulation_analysis.R                        |   4 -
 R/subgroup_consistency_helpers.R                   |  24 +-
 R/subgroup_consistency_main.R                      |   6 +-
 R/subgroup_search.R                                |  28 +-
 README.md                                          |  22 ++
 .../cc_task_docs_refs_maxminutes_2026-10-27.md     | 290 +++++++++++++++++++++
 man/build_cox_formula.Rd                           |   6 +-
 man/consistency_resample.Rd                        |   6 +-
 man/consistency_resample_compare.Rd                |   6 +-
 man/evaluate_consistency_twostage.Rd               |   6 +-
 man/evaluate_subgroup_consistency.Rd               |   6 +-
 man/forestsearch-package.Rd                        |   2 +-
 man/forestsearch.Rd                                | 245 ++++++++++++-----
 man/get_split_hr_fast.Rd                           |   6 +-
 man/guohe_adaptive_r.Rd                            |   6 +-
 man/guohe_algorithm3.Rd                            |   6 +-
 man/make_effect_estimator.Rd                       |  12 +-
 man/plot_sg_glm_outcomes.Rd                        |   6 +-
 man/plot_subgroup_results_forestplot.Rd            |   6 +-
 man/run_single_consistency_split.Rd                |   6 +-
 man/subgroup.consistency.Rd                        |   6 +-
 man/subgroup.search.Rd                             |  20 +-
 tests/testthat/test-max-minutes-deprecation.R      |  97 +++++++
 vignettes/articles/methodology.Rmd                 |   2 -
 vignettes/forestsearch.Rmd                         |   1 -
 vignettes/survival-post-selection.qmd              |   2 +-
 35 files changed, 959 insertions(+), 191 deletions(-)
```

## man/ MD5: before vs after

`dev/reports/docs_refs_2026-10-27/man_md5_before.txt` (at `a1b07bc9`, 427 files) vs `man_md5_after.txt` (427 files, same set). These 16 hashes changed:
```
man/build_cox_formula.Rd
man/consistency_resample_compare.Rd
man/consistency_resample.Rd
man/evaluate_consistency_twostage.Rd
man/evaluate_subgroup_consistency.Rd
man/forestsearch-package.Rd
man/forestsearch.Rd
man/get_split_hr_fast.Rd
man/guohe_adaptive_r.Rd
man/guohe_algorithm3.Rd
man/make_effect_estimator.Rd
man/plot_sg_glm_outcomes.Rd
man/plot_subgroup_results_forestplot.Rd
man/run_single_consistency_split.Rd
man/subgroup.consistency.Rd
man/subgroup.search.Rd
```
§5 expected only `forestsearch.Rd`, `subgroup.search.Rd` and the §1.1 item 4 pages. The other 13 come from §1.3 item 3:
- the 12 exported functions in the §1.3 table, including `subgroup.consistency.Rd`;
- `consistency_resample_compare.Rd` and `guohe_adaptive_r.Rd`, which inherit the edited params via `@inheritParams`.

`forestsearch-package.Rd` is the §1.1 item 4 page. `print.forestsearch.Rd` is unchanged (LA-1 only).

## Surprises
1. **§1.2 mechanism vs. `args_call_all`.** See the deviation section: the `mget` of all formals at `R/forestsearch_main.R` L1755 and the replaying callers.
2. **Handoff page list.** `fs_mr_inference` and `fs_betaHhat_one` do not cite 2609.38361. Only `forestsearch-package.Rd` (generated from `DESCRIPTION`) and `print.forestsearch.Rd` do.
3. **Unsuffixed "León and Anderson (2026)".** These citations predate the 2026a/2026b convention and were left unchanged: `R/fs_mr_inference.R:315`, `R/betaHhat_truth.R:212`, `R/forestsearch_methods.R:138` (printed caveat, ASCII "Leon & Anderson (2026)") and `:225`. The DESCRIPTION sentences also read "(2026)" twice, as §1.1 item 3 specifies verbatim. All of these are LA-1 except DESCRIPTION's second sentence.
4. **Stale note in `fs_family_report()`.** `R/fs_family_report.R:354` says "subgroup.search() forwards it and never compares it". forestsearch() no longer forwards it. Not edited (out of scope).
5. **Non-ASCII in R source.** The task's verbatim strings put `León` (é) into roxygen in 10 R files. CLAUDE.md asks for ASCII-only R source, but this repo already has `León` in roxygen (`R/forestsearch_methods.R:225`, `R/betaHhat_truth.R:212`, `R/fs_mr_inference.R:315`), and `DESCRIPTION` declares `Encoding: UTF-8`. The PDF manual (`--as-cran`) was not built here; Larry's release check covers it.
6. **README `@misc` entries** also carry a `doi = {10.48550/arXiv....}` field. That is beyond the fields §1.1 item 5 lists; remove it if unwanted.
7. **§0 line hint.** The legacy `analysis_gbsg_survival_effMaxSG.qmd` hit is at L580, not L≈584.
8. The git-ignored rendered `vignettes/survival-post-selection.html` still shows `max.minutes = 3` until the vignette is next rendered.

## Pre-check tidy (same day)

Task doc: `dev/tasks/cc_task_prechec_tidy_2026-10-27.md`. Start: `feature/docs-refs-maxminutes` at `b36da9dc0dae84e44ff774fddc33ae8978b54311`, clean tree.

**Disposition line:** A: yes  B: yes  C: yes

The task doc was first committed as received, with the disposition line blank. CC stopped there, as the doc requires. Larry supplied the dispositions in chat, and they went in as a separate commit.

### Commits

| SHA | Subject |
|---|---|
| `725d64b2` | dev: task doc — pre-check tidy (2026-10-27) |
| `c5dfe0b1` | dev: task doc — pre-check tidy dispositions (Larry, 2026-10-27) |
| `8f94d117` | DESCRIPTION: accept roxygen2 8.1.1 stamp (A) |
| `81d3dab3` | docs: fix stale note in fs_family_report.R (B) |
| `936b6734` | cleanup: remove dead fit_causal_forest_glm() (C) |
| `de298322` | docs: regenerate man pages |

### B — deviation from the task text

The task text for B says the note at L≈354 is a comment about the VI screen. In the working tree, neither part holds:

- `R/fs_family_report.R:354` is the note string on the `max.minutes` ("time cap") row. It is the one flagged in Surprises item 4 above. The VI-screen row (L266-269) already reads "deprecated and ignored since 0.4.0.9000 ...", with status `INE`, so it needed no change.
- The note is a string literal that the function prints. It is not a comment.

CC went by the doc's explicit pointer ("the one CC flagged in the report") and rewrote only that string. No logic changed:

- before: `"no path consults it (subgroup.search() forwards it and never compares it)"`
- after: `"deprecated and ignored since 0.4.0.9000: forestsearch() does not forward it and no path consults it"`

Basis for the new wording: `R/forestsearch_main.R:3477` (`search_args$max.minutes <- NULL`). No test or other file matches the old string.

### C.1 — `grep -rn "fit_causal_forest_glm" R/ tests/ vignettes/ inst/ quarto/ NAMESPACE man/` (before the edit)

```
R/grf_helpers.R:117:fit_causal_forest_glm <- function(X, Y, W, RCT, seedit,
man/fit_causal_forest_glm.Rd:3:\name{fit_causal_forest_glm}
man/fit_causal_forest_glm.Rd:4:\alias{fit_causal_forest_glm}
man/fit_causal_forest_glm.Rd:7:fit_causal_forest_glm(X, Y, W, RCT, seedit, tune_grf = FALSE)
```

### C.2 — gate: HOLDS

- (i) `NAMESPACE` has no hit, so there is no `export(fit_causal_forest_glm)`.
- (ii) The only `R/` hit is the definition (`R/grf_helpers.R:117`). Its roxygen block (L101-116, `@keywords internal`) does not name the function.
- (iii) There are no hits in `tests/`, `vignettes/`, `inst/` or `quarto/`.

C.3: the function and its roxygen block were deleted (`R/grf_helpers.R` L101-146 at `81d3dab3`, 46 lines). The block had no `@importFrom`, so there was nothing to keep. The NEWS line went under a new `### Internal` heading in the development section.

### §2.3 gates (after `devtools::document()`)

- `git diff --stat -- NAMESPACE`: empty.
- `man/`: the only change is `D man/fit_causal_forest_glm.Rd`.
- `git diff b36da9dc -- DESCRIPTION`: only `Config/roxygen2/version: 8.1.0` -> `8.1.1`. `document()` did not rewrite DESCRIPTION, so the restore step is no longer needed.

### §2.4 tests (`NOT_CRAN=true`)

| File | load_package | test_that blocks | expectations passed | failed | skipped | errors |
|---|---|---|---|---|---|---|
| `test-max-minutes-deprecation.R` | source | 5 | 12 | 0 | 0 | 0 |
| `test-vi-deprecation.R` | source | 4 | 15 | 0 | 0 | 0 |
| `test-lean-fits.R` | installed | 3 | 344 | 0 | 0 | 0 |

`devtools::install(upgrade = FALSE)` ran before `test-lean-fits.R`. `exists("fit_causal_forest_glm", envir = asNamespace("forestsearch"))` returns `FALSE` on the installed build. The first `test_file(..., load_package = "installed")` call failed because no `package =` argument was given (`'package' must be of length 1`). It was re-run with `package = "forestsearch"`. That was an invocation error; no test failed.

### §2.5

`git status --porcelain quarto/`: empty.
