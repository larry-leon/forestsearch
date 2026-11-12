# TASK — CRAN preparation, part 4: approved code fixes (2026-11-12)

Complete on its own. Run from the forestsearch repo root (`~/Documents/GitHub/forestsearch`).
Each fix below was authorized by Larry on 2026-11-12, one at a time, after the source was read.
This task runs unattended and never stops to ask.

It continues on branch `cran-prep`, after `TASK_cran_prep_docs_legacy_2026-11-12.md` (part 3).

## What this task changes

| Fix | Item | Kind |
|---|---|---|
| F1 | A-147 bootstrap summary bins match their labels | Changes behaviour |
| F2 | CV fold assignments {1, NA} or {0, NA} no longer stop | Changes behaviour (cases that stop today) |
| F3 | Missing-covariate messages at three sites | Adds messages only; no returned value changes |
| F4 | A-103 CV "Exact" means the same subgroup | Changes behaviour |
| F5 | A-151 `test_hte_from_forestsearch()` includes continuous covariates | Changes behaviour |
| F6 | A-042 LASSO: the intercept is not a selected variable | Changes behaviour |
| F7 | A-107 `fpr_calibration()` outcome-type fallback | Changes behaviour |
| F8 | A-047 duplicate-cut test compares names exactly | Changes behaviour |
| F9 | A-125 MRCT default confounders exclude internal spline columns | Changes behaviour |
| F10 | A-064 count outcomes with GRF are oriented | Changes behaviour |
| F11 | DINA orientation sign on the screening route | Changes behaviour |
| F12 | Three displays use the fit's own orientation | Changes behaviour |
| F13 | A-055 risk-difference standard error, unadjusted case | **Changes the method** (authorized as proposed) |
| F14 | DINA floor: fit, record and replay agree in the forced-gaussian corner | Changes behaviour |

**Not in this task, by Larry's decision:** `generate_bootstrap_with_noise()` and
`generate_bootstrap_synthetic()` keep their behaviour (A-150, A-162). The simulation
classification code still stops when subgroup membership contains NA (no `na.rm`).
The legacy GBSG censoring line is not changed.

## Rules

- No `git fetch`, `pull`, `push` or `merge`. Every `git add` names explicit paths.
- No install into the user or site library. Temporary libraries and scratch go under `/tmp`.
- fs-post-selection and fs-glms-interpretable: read-only, no branch change.
- `vignettes/` and `quarto/`: read-only. Nothing is re-rendered. The pkgdown site is not rebuilt.
- No simulation campaign is run. Reproducers use synthetic data, a sequential plan and
  `pkgload::load_all()` from the repo tree, cap 2 minutes each.
- Evidence is the current source. Line numbers below are pointers taken at `28e71ffa`
  (`cran-prep` before part 3); locate each site by content and quote file:line in the report.
  `dev/reports/REPORT_cran_prep_2026-11-11_part3_detail.md` describes each audit item.
- **Fixes are independent. One fix, one commit**, touching only that fix's files (code, its
  roxygen and `man/` page, its test file). Version and NEWS go in a separate last commit.
- **Per-fix gates, all required before the commit:**
  1. **Red, then green.** Write the test first. It must fail on the pre-fix tree and pass
     after the fix. Record both runs.
  2. **Independent expected values.** Every expected value in a test is worked out by hand or
     from the documentation, and the test says how in a comment. No expected value is copied
     from the fixed code's output.
  3. **Unaffected case.** The case named for the fix gives results identical to the start SHA.
  4. **Scope.** `git diff --stat` for the commit lists only that fix's files.
- If a fix fails a gate: restore its files to the pre-fix state, remove its test, report it
  as FAILED with the output, and go on. BLOCKED means a pre-set condition was not met and
  nothing was changed for that fix. Do not adjust a seed, a tolerance or an expected value to
  make a test pass.
- An existing expectation that fails only because it encoded the old behaviour of an
  approved fix may be updated; list each one with its old and new value. An existing
  expectation that fails for any other reason makes that fix FAILED.
- No timing estimates. Report measured wall clocks.

## Step 0 — branch and record

1. `git status --porcelain` is empty. **Gate.**
2. Check out `cran-prep`. **Gate:** `dev/reports/REPORT_cran_prep_docs_legacy_2026-11-12.md`
   exists on it. If not, stop: BLOCKED, part 3 has not finished. Record that report's status
   lines.
3. Copy `~/Downloads/TASK_cran_prep_code_fixes_2026-11-12.md` to `dev/tasks/` and commit it
   by explicit path (`dev: add CRAN preparation task, part 4`). This commit is the
   **start SHA**. Export it with `git archive` to `/tmp` for the equivalence checks.
4. Record: start SHA, `DESCRIPTION` version, and each sibling repository's branch, SHA and
   `git status --porcelain`.

## The fixes

### F1 — A-147: bootstrap summary bins (`R/summarize_bootstrap_subgroups.R`, near 417–449)

- **Change:** add `right = FALSE` to both `cut()` calls (consistency and size). Keep
  `include.lowest = TRUE`.
- **Test:** consistency values 0.95, 0.90, 0.85, 1.00, 0.949 fall in ">=0.95", "0.85-0.95",
  "0.85-0.95", ">=0.95", "0.85-0.95". Sizes 50, 49, 100, 300, 299 fall in "50-99", "<50",
  "100-149", ">=300", "200-299". These follow from the labels.
- **Unaffected case:** values not on a bin edge give identical tables.

### F2 — CV fold assignments with NA (`R/forestsearch_cross_validation.R`, near 1371, 1417–1445)

- **Change:** build the fold-by-original table before the branch. Take the two-level branch
  only when the table has two rows and two columns. Otherwise use the existing single-level
  branch. Do not change how `n_sgfound` is counted.
- **Test:** assignments {1, NA} and {0, NA} return the four agreement values the existing
  single-level rule defines, worked out by hand from the constructed data.
- **Unaffected case:** assignments {0, 1} and {0, 1, NA} give identical values.

### F3 — messages when subgroup membership is missing

No returned value changes at any of the three sites.

- **F3a. Simulation classification** (`R/run_simulation_analysis.R`, the forest search
  extractor near 1127–1138 and the GRF extractor near 1444 and 1536–1546). Where the code
  stops today because membership contains NA, stop instead with a message that gives the
  number of subjects with no subgroup assignment out of the total and says a covariate that
  defines the identified subgroup has missing values. **It must stop in exactly the cases
  that stop today.** No case that completes today may stop, and no `na.rm` is added.
- **F3b. GRF route with a binary, continuous or count outcome**, where the selected cut is
  applied to every subject (reported near `R/forestsearch_helpers.R:1987-1988`; confirm the
  site). When the resulting assignment contains NA, issue one warning with the count and the
  covariates that define the subgroup.
- **F3c. Cross-validation held-out assignments.** When they contain NA, issue one warning per
  exported call with the count. Issue it in the calling session after the folds are
  collected, not inside a worker.
- Name the covariates where the function has the subgroup definition at hand; give the count
  alone where it does not. Use `call. = FALSE`.
- **Tests:** F3a: the two NA shapes stop with the new message; a 0/1 membership completes
  with identical values. F3b and F3c: the warning appears once, its count equals the number
  of NA assignments, and the returned object is identical to the start SHA's.
- An existing test that now sees one of these warnings is changed to expect it by its
  message, not to suppress it. List each.

### F4 — A-103: CV "Exact" (`R/forestsearch_cross_validation.R`, `CV_sgs()`, near 1818–1826)

- **Change, one-factor branch only:** keep the "factor appears in the fold's subgroup"
  indicator for `one_match` and `cov1_match`. Set `exact_match` to 1 only when the fold's
  subgroup is that factor alone (the other slot is NA), with NA scored 0.
- **Test:** for full-data subgroup {A} and folds {A}, {A}&{B}, {B}&{A}, {B}, none:
  `exact_match` is 1 0 0 0 0; `one_match` and `cov1_exact` are 1 1 1 0 0.
- **Unaffected case:** a two-factor full-data subgroup, and no full-data subgroup, give
  identical output.

### F5 — A-151: `test_hte_from_forestsearch()` (`R/test_hte_crump.R`, near 397–416)

- **Change:** map each label in `confounders.evaluated` to the data columns it refers to
  (`all.vars()` of the parsed label, intersected with `names(df)`), and use the unique
  columns. Warn, naming the label, when a label cannot be resolved. Warn instead of returning
  `NULL` silently when no covariate remains. Do not change the `source` text.
- **Test:** on a fit with continuous and categorical confounders and no LASSO or GRF, the
  covariates used are the confounder names passed to `forestsearch()`, and the returned test
  equals `test_hte()` called directly on those columns. A fit with only continuous
  confounders returns a test.
- **Unaffected case:** a fit whose evaluated labels are all column names gives identical
  output.

### F6 — A-042: LASSO selection (`R/get_FSdata_helpers.R`, `lasso_selection()`, near 333–338)

- **Change:** build `selected` and `omitted` from the confounder rows only; `(Intercept)`
  stays in `coefficients`. The existing rule in `get_FSdata()` ("if LASSO selects nothing,
  keep the full pool") then applies to every outcome type. Make the printed "Cox-LASSO" label
  in `forestsearch_main.R` depend on the outcome type as `get_fsdata.R` already does, and
  correct the roxygen title that says the function is for the Cox model only.
- **Tests:** for a binary, a continuous and a count outcome, `(Intercept)` is in neither
  `selected` nor `omitted` and is in `coefficients`. With `lasso_selection()` mocked to
  select nothing, `get_FSdata()` for a non-survival outcome keeps the full candidate pool.
- **Unaffected case:** the survival path, and a non-survival fit where LASSO selects at least
  one covariate, give identical candidate cuts.

### F7 — A-107: `fpr_calibration()` (`R/fpr_calibration.R`, near 237–247)

- **Change:** when neither `outcome_type_approx` nor `fs_params$outcome_type` is given, the
  approximation type is `"survival"`, which is `forestsearch()`'s default. An unrecognised
  type stops with a clear message. Add one roxygen sentence stating the fallback.
- **Test:** with no type given, the recorded type is `"survival"` and the effective size
  equals `d_eff_survival()` evaluated in the test. A misspelt type stops.
- **Unaffected case:** with a valid type given, output is identical; with no type given, the
  simulated false-positive rate is identical.

### F8 — A-047: duplicate-cut test (`R/get_fsdata.R`, near 451, 499, 520)

- **Change:** at the GRF and DINA merges, parse the variable name from each cut (the text
  before `<=` or `>=`, trimmed). Drop the cut only if that name is exactly a categorical
  confounder or the whole cut string is exactly in the current pool. At the median-removal
  site, compare the continuous names exactly against the parsed GRF variable names.
- **Tests:** with confounders `bm` (categorical), `bmi`, `age`, the GRF cut `bmi <= 25.3` is
  kept, and likewise the DINA cut `bmi >= 25.3`. A GRF cut on a categorical confounder is
  still dropped. `x <= 25.3` is kept when the pool holds `x <= 2`. An exact duplicate is
  still dropped.
- **Unaffected case:** confounder names with no overlap give identical candidate pools.

### F9 — A-125: MRCT default confounders (`R/mrct_simulation.R` near 316–325; `R/mrct_case_study.R` near 286–293)

- **Change:** when `confounders.name` is NULL and the mechanism has spline information,
  remove the three internal columns named `<spline variable>` followed by `_treat`, `_k` and
  `_k_treat`. Build those exact names; do not match on a suffix. Make the identical change in
  both functions and state the exclusion in the roxygen. The region indicator stays.
- **Test:** for a mechanism built with a spline, both functions' default set equals the list
  written out by hand, and it contains `z_prior_treat`.
- **Unaffected case:** with `confounders.name` given, output is identical.

### F10 — A-064: count outcomes with GRF (`R/grf_subg_harm_glm.R`, near 443–452)

- **Change:** when `adverse_outcome` is TRUE and the outcome is a count, negate the
  transformed outcome before the forest is fitted. Correct the comments near 437–441, 451 and
  515–517 that say no flip is needed.
- **Test:** on made-up count data where treatment raises the count for one known region and
  lowers it elsewhere, with a higher count adverse, the identified subgroup lies mostly in
  the region where treatment raises the count. The memberships under `adverse_outcome = TRUE`
  and `FALSE` differ.
- **Unaffected case:** `adverse_outcome = FALSE` for counts, and binary and continuous
  outcomes, give identical output.

### F11 — DINA orientation sign on the screening route (`R/forestsearch_main.R`, near 3254–3265)

- **Change:** add `tau_sign = .dina_tau_sign(outcome_type, adverse_outcome)` to the
  `dina_subgroup()` call on the `use_dina = TRUE` route, as the `subgroup_method = "dina"`
  route already does (`R/forestsearch_helpers.R` near 1542, 1604). Nothing else on either
  route changes.
- **Tests:** the screening route calls `dina_subgroup()` with the sign that `.dina_tau_sign()`
  returns for the fit's outcome type and orientation. On made-up continuous data with a
  known harm region, the cut the screening route proposes points at that region.
- **Unaffected case:** survival and count outcomes, binary with `adverse_outcome = TRUE`, and
  every `subgroup_method = "dina"` fit give identical output.

### F12 — displays use the fit's own orientation

Sites reported in `dev/reports/REPORT_cran_prep_2026-11-11.md` OPEN ITEMS 3:
`R/plot_sg_glm_outcomes.R` near 254–262, `R/plot_subgroup_results_forestplot.R` near 431–439,
and `mr_estimates_table()` with `boot_results = NULL` (`R/mr_estimates_table.R` near 204–208).

- **Change:** where each display builds an effect for a `forestsearch()` fit, pass the fit's
  stored `adverse_outcome`. If it is not stored, use `forestsearch()`'s own rule, as the GLM
  bootstrap fix did.
- **Not changed:** the arm means and their raw difference in the summary table
  (`R/summary_utility_functions.R` near 255) stay raw.
- **Test:** for a continuous-outcome fit with the default orientation for that type, each
  display's estimate for the subgroup equals the fit's own estimate, sign included.
- **Unaffected case:** survival fits, and non-survival fits with `adverse_outcome = TRUE`,
  give identical display data.

### F13 — A-055: risk-difference standard error (`R/glm_effect_estimators.R`, `.estimate_rd`, near 514–548)

**Change of method, authorized as proposed.**

- **Change:** in the logistic tier, when there are no adjustment covariates and no weights,
  the standard error is `sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)`, with `p0`, `p1` the
  observed arm proportions and `n0`, `n1` the arm sizes. Otherwise the existing expression is
  kept. The estimate, `converged` and `method_used` do not change. Do not route these cases
  to the raw-proportion tier.
- Confirm how weights reach this function and report it. If the unadjusted case cannot be
  separated from a weighted one, report the fix BLOCKED.
- **Tests,** through `make_effect_estimator(outcome_type = "binary", effect_measure = "RD")`:
  with 40 per arm, control 30% and treated 100%, the standard error is sqrt(0.21 / 40) =
  0.07246; with control 0% and treated 80%, it is sqrt(0.16 / 40) = 0.06325. In both, the
  estimate is the difference of the arm proportions and `converged` is TRUE.
- **Unaffected case:** the odds-ratio and risk-ratio measures, and a risk difference with
  adjustment covariates, give identical output.
- **Record as OPEN, do not fix:** the covariate-adjusted case uses the same shortcut; with
  both arms at 0% or 100% the standard error is 0.

### F14 — DINA floor in the forced-gaussian corner (`R/forestsearch_main.R`, near 2784–2785 and 3251–3253)

Setting: binary risk difference or count rate difference, `dina_args = list(family = "gaussian")`,
and a supplied ratio-scale threshold that the package warns about and converts.

- **Change:** the DINA fit's floor is the resolved threshold at both sites, so the fit, the
  recorded value and the bootstrap and CV replays use the same floor. Update the pinning
  block in `tests/testthat/test-dina-gaussian-floor.R` (near 123) to the new value.
- **Test:** in that setting, the floor used by the fit equals the recorded value.
- **Unaffected case, as a table in the report:** for every other threshold configuration the
  existing test file covers, and for a continuous outcome with a supplied threshold, the
  floor used by the fit is unchanged from the start SHA. If any row other than the corner
  differs, the fix is FAILED.

## After the fixes

1. **Default-path gate.** Re-run the equivalence cases of
   `dev/reports/REPORT_cran_prep_seed_legacy_2026-11-11.md` (the six gate cases) at the head
   and against the start-SHA export. A default survival `forestsearch()` run must be
   identical. Any difference elsewhere must be traced to a named fix, with the mechanism
   shown; an unexplained difference is FAILED for the task and is reported, not repaired.
2. **`paper_simulations` article.** Do not re-render. State, per fix, whether any call the
   article evaluates reaches the changed lines. If one does, list it under OPEN ITEMS.
3. `DESCRIPTION`: next development version. `NEWS.md`: one entry per fix that passed, saying
   what a user sees change. Run `devtools::document()`. Commit by explicit path
   (`chore: version and NEWS for part 4`).
4. **Full suite** (`devtools::test()`, `NOT_CRAN=true`), cap 30 minutes. **Gate:** 0 failures,
   0 errors. Report blocks and expectations beside part 3's and account for the difference.
5. **`devtools::check(document = FALSE, error_on = "never")`** on the check route of
   `dev/reports/REPORT_cran_prep_2026-11-11.md` §6, cap 60 minutes. **Gate:** 0 errors,
   0 warnings, 0 notes.
6. Write `dev/reports/REPORT_cran_prep_code_fixes_2026-11-12.md`, in bullets and tables:
   - provenance and commits;
   - one status line per fix (PASS, FAILED or BLOCKED) with its commit SHA;
   - per fix: the site quoted before and after, the red and green test runs, how each
     expected value was derived, the unaffected-case result, and **one before/after example
     of what a user sees**;
   - existing expectations that were updated, with old and new values;
   - the default-path gate, the `paper_simulations` statement, suite and check with wall
     clocks;
   - OPEN ITEMS.
   Save any gate script and its output beside the report.
7. Commit the report by explicit path (`dev: CRAN preparation report, part 4`) and copy it to
   `~/Downloads`.
8. **Final gates:** `git status --porcelain` is empty; both sibling repositories are
   unchanged; `git diff <start> --stat` lists only files the report names, and nothing under
   `vignettes/` or `quarto/`. Nothing is pushed or merged.
