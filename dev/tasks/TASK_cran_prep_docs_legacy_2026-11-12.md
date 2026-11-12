# TASK — CRAN preparation, part 3: documentation corrections, legacy GBSG wrapper, GBSG inventory (2026-11-12)

Complete on its own. Run from the forestsearch repo root (`~/Documents/GitHub/forestsearch`).
Authorized by Larry on 2026-11-12. This task runs unattended and never stops to ask.

It continues on branch `cran-prep`, after `TASK_cran_prep_2026-11-11.md` and
`TASK_cran_prep_seed_legacy_2026-11-11.md`. Larry has merged `cran-prep` into `master`.

## Purpose and kind of change

| Part | What | Kind |
|---|---|---|
| 1 | Remove four legacy GBSG functions that nothing calls | Changes the exported interface. No remaining function changes behaviour. |
| 2 | Legacy GBSG wrapper: mark superseded, document the censoring step | Documentation and source comments only |
| 3 | Synthetic-data functions: help pages match the code; tests pin what the help says | Documentation, source comments and new tests only |
| 4 | Four audit findings handled by documentation | Documentation and source comments only |
| 5 | Inventory of GBSG-specific functions | Read-only |
| 6 | Version, NEWS, suite, check, code-identity gate, report | — |

**Larry's decision that governs the whole task: no function's behaviour changes.** The only
code that changes is the deletion in Part 1. If any step would alter what an existing function
computes or returns, do not make it; report it as BLOCKED.

## Rules

- No `git fetch`, `pull`, `push` or `merge`. Every `git add` names explicit paths.
- No install into the user or site library. Temporary libraries and scratch go under `/tmp`.
- fs-post-selection and fs-glms-interpretable: read-only, no branch change.
- `vignettes/` and `quarto/`: read-only in this task. The pkgdown site is not rebuilt.
- No new package dependency (do not add `lifecycle`).
- No run-time message, warning or `.Deprecated()` call is added anywhere.
- Evidence is the current source. Quote file:line in the report.
- **The code is the truth.** Where this task gives help-page wording, check every sentence
  against the code before writing it. If the code differs from a sentence, write what the code
  does and record the difference in the report.
- Line numbers below are pointers taken at `e92d3c96` or `28e71ffa`. Locate each site by content.
- Gates are stop-on-failure for their own part, never stop-to-ask. If a part fails a gate,
  restore that part's files to their pre-part state, report the part as FAILED with the
  output, and go on to the next part. BLOCKED means a pre-set condition was not met and
  nothing was changed for that item. The report is always written.
- Optional reproducers to confirm a help-page sentence: synthetic or `survival::gbsg` data,
  sequential plan, `pkgload::load_all()` from the repo tree, `Rscript --vanilla` from
  `/tmp/cran_prep_docs_legacy_2026-11-12/`, cap 2 minutes each and 15 minutes in total.
- No timing estimates. Report measured wall clocks.

## Step 0 — branch and record

1. `git status --porcelain` is empty. **Gate.**
2. If branch `cran-prep` exists, check it out. Otherwise create it from `master`.
3. **Gate:** `git diff --quiet master cran-prep -- R man tests DESCRIPTION NAMESPACE NEWS.md vignettes _pkgdown.yml`
   exits 0 (the branch and `master` hold the same package). If not, stop: BLOCKED, the
   branches disagree and this task does not merge.
4. Copy `~/Downloads/TASK_cran_prep_docs_legacy_2026-11-12.md` to `dev/tasks/` and commit it
   on `cran-prep` by explicit path (`dev: add CRAN preparation task, part 3`). This commit is
   the **start SHA**.
5. Record: start SHA, `DESCRIPTION` version (expected 0.4.1.9008), and for each sibling
   repository its branch, SHA and `git status --porcelain`.

## Part 1 — remove four uncalled legacy GBSG functions

Functions: `create_gbsg_dgm()`, `simulate_from_gbsg_dgm()`, `.simulate_from_gbsg_dgm_()`,
`get_dgm_with_output()`.

1. **Gate, per function:** no caller at the start SHA. Search `R/`, `tests/`, `vignettes/`,
   `inst/`, `quarto/` (mark archive or legacy directories) and, read-only, both sibling
   repositories. A caller is code that calls the function; the function's own definition,
   roxygen, help page, NAMESPACE line, `_pkgdown.yml` entry and tests of that function alone
   do not count. A caller only in an archive or legacy `quarto/` directory does not block;
   list it. If a current caller exists, leave that function in place and report it BLOCKED.
   `dev/reports/REPORT_cran_prep_seed_legacy_2026-11-11.md` Part 3 is the starting inventory.
2. Remove each cleared function: its definition and roxygen, its `man/*.Rd`, its NAMESPACE
   export (through `devtools::document()`), its `_pkgdown.yml` reference entry, and any
   `test_that()` block that tests only that function.
3. An internal helper referenced only by removed functions may be removed too. List each one
   in the report. When in doubt, keep it.
4. **Keep, unchanged:** `setup_gbsg_dgm()`, `calibrate_k_inter()`, `.create_gbsg_dgm_()`,
   `validate_k_inter_effect()`, `fs_fdr_report()`, the `gbsg_dgm` class and its methods, and
   everything they reach.
5. Fix every remaining reference to a removed name in `R/` roxygen, `man/`, `tests/`,
   `_pkgdown.yml` and the forestsearch `README` (`.Rmd` and `.md` alike, text only, no
   re-knit). In prose or a list, point to `generate_aft_dgm_flex()` with
   `simulate_from_dgm()`. If the reference is in code that is evaluated, leave it and report
   it BLOCKED.
6. List, do not edit: mentions in `vignettes/`, `quarto/` and the sibling repositories
   (including the post-selection supplement near `:1184`).
7. **Gates:**
   - `git grep -nE 'create_gbsg_dgm|simulate_from_gbsg_dgm|get_dgm_with_output'` over `R/`,
     `man/`, `tests/`, `NAMESPACE`, `_pkgdown.yml`, `README*` returns nothing for the removed
     names (`.create_gbsg_dgm_` is kept and is expected to match).
   - `git diff <start> -- NAMESPACE` shows only deleted lines, each naming a removed function.
   - The test files that exercise `setup_gbsg_dgm()`, `calibrate_k_inter()` and
     `fs_fdr_report()`: 0 failures.
8. Commit by explicit path (`refactor: remove four uncalled legacy GBSG functions`).

## Part 2 — legacy GBSG wrapper: superseded, censoring step documented

No code change. Roxygen and source comments only.

1. **Confirm at source and quote in the report:**
   - (a) the censoring covariate term in `.create_gbsg_dgm_()` (`R/sim_aft_gbsg.R` near `:638`)
     is the fitted term multiplied by (event-model scale / censoring-model scale);
   - (b) `generate_aft_dgm_flex()` takes each subject's censoring linear predictor from the
     censoring fit itself, with no such factor;
   - (c) the true hazard ratios and the `k_inter` calibration do not use the censoring
     predictor.
   If (b) cannot be confirmed, omit the sentence that compares the two generators and
   describe only the wrapper's own step. If (a) or (c) does not hold, write nothing for this
   part and report it BLOCKED.
2. **`setup_gbsg_dgm()`**: add this roxygen section (plain text, no badge):

   > `@section Superseded:` This function is kept so that published simulations can be
   > reproduced. New work should use `generate_aft_dgm_flex()` with `simulate_from_dgm()`,
   > which take any survival dataset. On the GBSG data they give the same outcome model, true
   > hazard ratios and `k_inter` calibration once the caller has built the GBSG covariates;
   > there is no drop-in replacement call. The censoring step differs: here the covariate term
   > of the censoring model is the fitted term multiplied by the ratio of the event-model scale
   > to the censoring-model scale. True values are not affected; simulated censoring times are.
   > The "Extreme Subgroup Effects" article shows the general machinery on the GBSG data. This
   > function will be deprecated in a future release.

3. **`calibrate_k_inter()`**: the same section, shortened to its own role. Name a replacement
   only if an exported function exists whose documented purpose is calibrating `k_inter` for
   `generate_aft_dgm_flex()`; otherwise give no replacement name.
4. **Every other exported function whose simulated data pass through `.create_gbsg_dgm_()`**
   (expected: `fs_fdr_report()`, and `validate_k_inter_effect()` if exported; determine the
   list from the code and report it): add a `@note` saying it builds its data through the
   superseded GBSG wrapper and stating the censoring step in one sentence.
5. **Source comment** at the censoring line: two or three `#` lines saying the factor is kept
   as is so that published simulations reproduce, that true values do not depend on it, and
   pointing to the help page. No code on that line changes.
6. Run `devtools::document()`.
7. **Gate:** for every file under `R/` changed in this part, each changed line is a comment
   (`#` or `#'`):
   `git diff -U0 <pre-part> -- <file> | grep -E '^[+-]' | grep -vE '^(\+\+\+|---)' | grep -vE "^[+-][[:space:]]*#"`
   returns nothing.
8. Commit by explicit path (`docs: legacy GBSG wrapper marked superseded; censoring step documented`).

## Part 3 — synthetic-data functions (`R/synthetic_data_perturbation.R`)

**No code change. Larry's decisions:** `generate_bootstrap_with_noise()` treating every
numeric column as continuous under auto-detection is by design. `generate_bootstrap_synthetic()`
keeps its behaviour, including the loop over `cat_vars` only (near `:145`) and the
`sample(possible_new, 1)` draw (near `:168`). The help pages are corrected to say what the
code does.

### 3.1 `generate_bootstrap_with_noise()`

Replace the wording of these roxygen pieces. Everything else on the page stays.

- **Description:** "Creates a bootstrap sample from a dataset and adds controlled noise.
  Numeric columns receive Gaussian noise. Factor and logical columns, and any column named in
  `cat_vars`, are changed with a small probability. Numeric-coded categories (for example 0/1
  indicators) are treated as continuous unless they are named in `cat_vars`. This function is
  useful for generating synthetic datasets that maintain the general structure of the original
  data while introducing controlled variation."
- **`@param continuous_vars`:** "Character vector of column names to treat as continuous
  variables. If NULL (default), every numeric column except `id_var` is used, including
  numeric columns with few distinct values."
- **`@param cat_vars`:** "Character vector of column names to treat as categorical variables.
  If NULL (default), factor and logical columns are used, together with numeric columns that
  have 10 or fewer unique values and are not in `continuous_vars`. When `continuous_vars` is
  also NULL, every numeric column is already in the continuous set, so no numeric column is
  treated as categorical. When `cat_vars` is given and `continuous_vars` is NULL, a numeric
  column named here is also in the auto-detected continuous set and receives both kinds of
  change; give `continuous_vars` as well to avoid this."
- **Details, continuous section:** add "Integer-valued columns are rounded after the noise is
  added, so a 0/1 column treated as continuous is usually unchanged at a small `noise_level`."
- **Details, categorical section:** begin with "Applies to the columns in `cat_vars`, given or
  auto-detected as described for that argument." Change the ordered-factor item to "Ordered
  factors, or a column named `grade`, with more than two levels: selection is weighted toward
  adjacent levels."
- **`@note`, third item:** "For an ordered factor, or a column named `grade`, that is
  processed as categorical and has more than two levels, the perturbation favors adjacent
  levels. A numeric `grade` column is processed as categorical only when it is named in
  `cat_vars`, or when `continuous_vars` is given without it."
- **Examples:** replace `data(gbsg, package = "survival")` with `gbsg <- survival::gbsg`. On
  the first example add the comment "all gbsg columns are numeric, so all are treated as
  continuous".
- **Source comment** above the auto-detection block: "By design (L. Leon, 2026-11-12): under
  auto-detection every numeric column is continuous. Not a defect; see the help page."

### 3.2 `generate_bootstrap_synthetic()`

- **`@param ordinal_vars`:** "Optional character vector naming which of the `cat_vars` are
  ordinal. A variable listed in both `cat_vars` and `ordinal_vars` is changed to a
  neighbouring level instead of a random other level. A variable listed only in
  `ordinal_vars` is not perturbed: it is kept as resampled and is named in the 'kept as-is'
  message. Known limitation: for a numeric-coded variable at its lowest or highest level, the
  new value is not always the adjacent level."
- **Both examples:** add a comment at the `ordinal_vars` line: "listed only in ordinal_vars,
  so it is resampled but not perturbed; list it in cat_vars as well to perturb it".
- **Source comments:** at the categorical loop, "Behaviour kept as is by decision
  (2026-11-12): the loop runs over cat_vars only; see the help page for ordinal_vars." At the
  `sample(possible_new, 1)` line, "Known limitation, documented in the help page; behaviour
  kept as is by decision (2026-11-12)."

### 3.3 `generate_gbsg_bootstrap_general()`

- Add to the description: "`grade` is passed in `ordinal_vars` only, so it is resampled but
  not perturbed."
- Add: `@section Superseded:` "This GBSG-specific wrapper is kept for existing code. New work
  should call `generate_bootstrap_synthetic()` directly: the call in its first example, with
  the same `n`, `seed` and `noise_level`, returns the same data. This wrapper will be
  deprecated in a future release."

### 3.4 Same `data(gbsg, ...)` example pattern elsewhere

In the `\dontrun` examples at `R/pareto_frontier_table.R` (near `:124`) and
`R/calibrate_k_treat.R` (near `:55`), make the same replacement with `survival::gbsg`.
`vignettes/articles/extreme_subgroups.qmd:28` is not edited; list it.

### 3.5 Tests that pin what the help now says

New file `tests/testthat/test-synthetic-data-documented-behaviour.R`. Expected values come
from the resample, not from the functions: after `set.seed(seed)`, the first random call in
each function is `sample(nrow(data), n, replace = TRUE)`, so the test rebuilds those row
indices itself.

- **T1.** `generate_gbsg_bootstrap_general(n = 200, seed = 11)`: `grade` equals
  `survival::gbsg$grade[idx]`.
- **T2.** The wrapper's output is `identical()` to `generate_bootstrap_synthetic()` called
  with the lists of its first example and the same `n`, `seed`, `noise_level`.
- **T3.** `generate_bootstrap_with_noise()` under auto-detection, on a made-up data frame of
  2,000 rows with an id, a continuous column, a numeric 0/1 column and a three-level factor,
  `noise_level = 0.1`: the 0/1 column equals its resample; the factor differs from its
  resample in more than 0 and fewer than 15% of rows.
- **T4.** The same data with `cat_vars` naming the 0/1 column and `continuous_vars` naming the
  continuous column: the 0/1 column differs from its resample in more than 0 and fewer than
  15% of rows.
- Do not change a seed or a bound to make a test pass. If a test fails, the help wording is
  wrong: correct the wording to match the code and the test to match the wording, and record
  it.

### 3.6 Gates and commit

- The comment-only gate of Part 2 step 7, for every `R/` file changed in this part.
- The new test file and `tests/testthat/test-gbsg-bootstrap-general.R`: 0 failures.
- Run `devtools::document()`. Commit by explicit path
  (`docs: synthetic-data help pages match the code; tests pin the documented behaviour`).

## Part 4 — four audit findings handled by documentation

No code change, no printed label change, no run-time warning. Sites and wording come from
`dev/reports/REPORT_cran_prep_2026-11-11_part3_detail.md`; use the documentation-only option
of each row.

| ID | Correction |
|---|---|
| A-024 | Where the help says DINA's Poisson quantity is "log(IRR)" (`R/dina_subgroup.R` near `:140`, `R/forestsearch_main.R` near `:445`, and the description text in `R/dina_subgroup_refit.R` near `:43`): say it is the log mean-count ratio, that DINA takes no exposure offset, and that it equals log(IRR) only when follow-up does not depend on arm given the covariates. Do not change the printed label near `dina_subgroup_refit.R:447`. |
| A-062 | `forestsearch()` `@param dmin.grf` and `@param grf_select_statistic`: under `"effect"`, `dmin.grf` is applied first as an existence gate on the doubly robust scale and is then not applied per candidate. Correct the "Minimum events for GRF" text near `forestsearch_main.R:555`. |
| INV-5 | Correct the roxygen near `R/subgroup_consistency_main.R:263-265` and the comment near `:714-715` that say results are identical on every path. Add one sentence to `forestsearch()`'s `parallel_args` and `consistency_method` documentation: a sequential and a parallel plan can give different consistency proportions, and possibly a different selection, whenever literal splitting runs (`consistency_method = "split"`, a risk or rate difference under `"resample"`, or a candidate whose closed form is not available). |
| A-156 | `fs_fdr_report()`: one roxygen sentence saying `k_inter` is calibrated on the default super-population size and not on `n_super`, and that the calibrated value does not depend on that size. |

- Run `devtools::document()`.
- **Gate:** the comment-only gate of Part 2 step 7, for every `R/` file changed in this part.
- Commit by explicit path
  (`docs: DINA count scale, dmin.grf gate, plan-dependent splitting, fs_fdr_report calibration size`).

## Part 5 — inventory of GBSG-specific functions (read-only)

Larry wants to phase these out. A function is GBSG-specific if it loads or assumes the GBSG
data itself (reads `survival::gbsg` or `gbsg`, or relies on GBSG column names) and cannot be
given another dataset. A general function that only uses GBSG in an example is not
GBSG-specific; do not list it.

One table, every function under `R/` at the branch head, exported or internal:

| Function | File:line | Exported | In CRAN 0.4.0 | What is hard-wired | Callers: package code | Callers: tests | Callers: vignettes and articles | Callers: current quarto (count) | Callers: fs-post-selection, fs-glms-interpretable (count) | General equivalent, and what a caller must add to match it |
|---|---|---|---|---|---|---|---|---|---|---|

Then two short lists:

- general functions that carry a GBSG leftover inside them (a default such as
  `id_var = "pid"`, a special case for a column named `grade`), with file:line;
- help-page examples that load GBSG in a way that fails as written.

Change nothing for Part 5.

## Part 6 — version, NEWS, checks, report

1. `DESCRIPTION`: next development version (expected 0.4.1.9009). `NEWS.md`: entries for the
   four removed functions, the superseded status of `setup_gbsg_dgm()`, `calibrate_k_inter()`
   and `generate_gbsg_bootstrap_general()`, the censoring-step note, and the help-page
   corrections. Commit by explicit path (`chore: version and NEWS for part 3`).
2. **Code-identity gate.** Export `R/` at the start SHA and at the head with `git archive`
   into `/tmp`. With an R script, `parse(keep.source = FALSE)` every file in each export and
   deparse each top-level expression. Compare by assigned name:
   - names only at the start = exactly the functions and helpers Part 1 reports as removed;
   - names only at the head = none;
   - names in both with a different deparsed body = none;
   - any other top-level difference is listed, and is allowed only as the removal of names
     from `utils::globalVariables()`.
   Save the script and its output beside the report
   (`REPORT_cran_prep_docs_legacy_2026-11-12_code_identity.R`, `..._code_identity.txt`).
3. **Help pages:** `git diff --stat <start> -- man/` lists only the pages of functions this
   task edited or removed.
4. **Full suite** (`devtools::test()`, `NOT_CRAN=true`), cap 30 minutes. **Gate:** 0 failures,
   0 errors. Report blocks and expectations beside the baseline (648 blocks, 8,471
   expectations) and account for the difference.
5. **`devtools::check(document = FALSE, error_on = "never")`** on the check route of
   `dev/reports/REPORT_cran_prep_2026-11-11.md` §6, cap 60 minutes. **Gate:** 0 errors,
   0 warnings, 0 notes. A failure caused by a help-page or roxygen syntax error introduced by
   this task may be repaired, at most twice; anything else is FAILED and is left as it is.
6. Write `dev/reports/REPORT_cran_prep_docs_legacy_2026-11-12.md`, in bullets and tables:
   - provenance and commits;
   - a status line per part (PASS, FAILED or BLOCKED);
   - Part 1: what was removed, the caller searches, listed mentions not edited;
   - Parts 2 to 4: for each help page, the old and the new text, and each sentence that was
     changed from this task's wording because the code differed;
   - Part 3 test results;
   - the Part 5 table and lists;
   - the code-identity result, suite and check beside the baseline, with wall clocks;
   - OPEN ITEMS.
7. Commit the report and its two code-identity files by explicit path
   (`dev: CRAN preparation report, part 3`). Copy the report to `~/Downloads`.
8. **Final gates:** `git status --porcelain` is empty; both sibling repositories are
   unchanged; `git diff <start> --stat` lists only files the report names, and nothing under
   `vignettes/` or `quarto/`. Nothing is pushed or merged.
