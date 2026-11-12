# TASK — CRAN preparation, part 5: remaining displays, adjusted standard error, documentation (2026-11-12)

Complete on its own. Run from the forestsearch repo root (`~/Documents/GitHub/forestsearch`).
Authorized by Larry on 2026-11-12, after the part 4 report. This task runs unattended and
never stops to ask.

It continues on branch `cran-prep`, after `TASK_cran_prep_code_fixes_2026-11-12.md` (part 4).

## What this task changes

| Part | Item | Kind |
|---|---|---|
| A1 | `sg_tables()` uses the fit's own orientation | Changes behaviour |
| A2 | CV summary tables use the fit's own orientation | Changes behaviour |
| A3 | `compute_frontier_cis()` fallback follows `forestsearch()`'s rule | Changes behaviour |
| A4 | Any further display of a fit found by the sweep | Changes behaviour |
| B1 | Risk-difference standard error, covariate-adjusted logistic tier | **Changes the method** (formula approved by Larry) |
| B2 | Rate-difference branch with the same shortcut | **Changes the method**, only if the code has the form stated below |
| C1–C3 | Three printed labels | Changes printed text only |
| C4 | Help pages and source comments | Documentation only |
| C5 | Examples that cannot run as written | Documentation only |
| C6 | Comments on test expectations | Comments only |
| C7 | Two article sources | Prose only, no re-render |

**Not in this task, by Larry's decisions:** everything else listed as open in the part 3 and
part 4 reports stays in the record. In particular: `generate_bootstrap_with_noise()` and
`generate_bootstrap_synthetic()` keep their behaviour; the legacy GBSG censoring line is not
changed; the simulation classification still stops on NA membership; the arm means and their
raw difference in the summary tables stay raw.

## Rules

- No `git fetch`, `pull`, `push` or `merge`. Every `git add` names explicit paths.
- No install into the user or site library. Temporary libraries and scratch go under `/tmp`.
- fs-post-selection and fs-glms-interpretable: read-only, no branch change.
- `quarto/`: read-only. `vignettes/`: only the two `.qmd` files of C7, prose only. Nothing is
  re-rendered. The pkgdown site is not rebuilt. No simulation campaign is run.
- No new package dependency.
- Evidence is the current source. Line numbers below are pointers taken at `61c15088` or
  `d93e6621`; locate each site by content and quote file:line in the report.
- **The code is the truth.** Where this task gives help-page wording or names a site, check it
  against the code first. If the code differs, write what the code does and record the
  difference.
- **One fix, one commit**, touching only that fix's files. Version and NEWS go in a separate
  last commit.
- **Gates for every code fix (A1–A4, B1, B2, C1–C3), all required before its commit:**
  1. **Red, then green.** The test is written first, fails on the pre-fix tree and passes
     after. Record both runs.
  2. **Independent expected values.** Each is worked out by hand or from the documentation,
     with the working in a comment. None is copied from the fixed code's output.
  3. **Unaffected case.** The case named for the fix is identical to the start SHA.
  4. **Scope.** `git diff --stat` for the commit lists only that fix's files.
- **Gate for C4, C5, C6:** every changed line in a file under `R/` or `tests/` is a comment
  (`#` or `#'`):
  `git diff -U0 <pre> -- <file> | grep -E '^[+-]' | grep -vE '^(\+\+\+|---)' | grep -vE "^[+-][[:space:]]*#"`
  returns nothing.
- A fix that fails a gate is restored to its pre-fix state and reported FAILED; go on to the
  next. BLOCKED means a pre-set condition was not met and nothing was changed. Do not adjust
  a seed, a tolerance or an expected value to make a test pass.
- An existing expectation that fails only because it encoded the old behaviour of a fix in
  this task may be updated; list each with its old and new value. Any other failure of an
  existing test makes that fix FAILED.
- Reproducers: synthetic data or `survival::gbsg`, sequential plan, `pkgload::load_all()`
  from the tree, cap 2 minutes each.
- No timing estimates. Report measured wall clocks.

## Step 0 — branch and record

1. `git status --porcelain` is empty. **Gate.**
2. Check out `cran-prep`. **Gate:** `dev/reports/REPORT_cran_prep_code_fixes_2026-11-12.md`
   exists on it. If not, stop: BLOCKED, part 4 is not on this branch.
3. Copy `~/Downloads/TASK_cran_prep_displays_se_docs_2026-11-12.md` to `dev/tasks/` and commit
   it by explicit path (`dev: add CRAN preparation task, part 5`). This commit is the
   **start SHA**. Export it with `git archive` to `/tmp` for the equivalence checks.
4. Record: start SHA, `DESCRIPTION` version (expected 0.4.1.9010), and each sibling
   repository's branch, SHA and `git status --porcelain`.

## Part A — displays use the fit's own orientation

**The rule, as in part 4's F12:** where a display builds an effect for a `forestsearch()` fit,
it passes the fit's stored `adverse_outcome`. If none is stored, it uses `forestsearch()`'s own
rule (`outcome_type %in% c("binary", "count")`). Arm means and their raw difference are not
touched.

- **A1. `sg_tables()`** — `R/summary_utility_functions.R`, the GRF branch near 829–839 and the
  main branch near 982–992.
- **A2. CV summary tables** — `R/forestsearch_cross_validation.R` near 1509–1517 and 1565–1575.
- **A3. `compute_frontier_cis()`** — `R/frontier_cis.R` near 150: replace the `%||% TRUE`
  fallback with the rule above. A stored value is used exactly as now.
- **A4. Sweep.** List every `make_effect_estimator(` call under `R/` with: file:line, the
  function it is in, whether it passes `adverse_outcome`, and what it is for (fitting,
  resampling, display, other). Fix a further site only if it is in a function that takes a
  `forestsearch()` fit and builds an effect for a table or plot without passing the
  orientation; give each its own commit and test. Anything else is listed, not changed.
- **Test, each of A1–A4:** for a continuous-outcome fit at the default orientation, the
  display's subgroup estimate equals the difference in arm means worked out by hand, oriented
  as the fit's own estimate (`max_sg_est`), sign included.
- **Unaffected case, each:** survival fits, and non-survival fits with
  `adverse_outcome = TRUE`, give identical display data. For A3, any object that stores the
  value gives identical output.

## Part B — standard errors by the full delta method

### B1 — covariate-adjusted risk difference (`R/glm_effect_estimators.R`, `.estimate_rd()`, logistic tier, near 534–547)

**Change of method, approved.** Today, with adjustment covariates, the standard error is
`abs(avg_deriv) * se_treat`: it uses only the treatment coefficient.

- **New, for the adjusted case only (`n_adj > 0`):** with `p1`, `p0` the fitted probabilities
  for every subject under treatment and under control, `X1`, `X0` the corresponding design
  matrices (estimable columns only), and `V` the covariance matrix of the logistic fit:

  ```
  g  = colMeans(p1 * (1 - p1) * X1) - colMeans(p0 * (1 - p0) * X0)
  SE = sqrt( g' V g )
  ```

- **Safety rule:** if that value is not a finite, non-negative number, keep the existing
  expression. No case may become an error.
- **Unchanged:** the unadjusted branch (part 4's F13), the estimate, `converged`,
  `method_used`, every other effect measure and every other tier. A small internal helper for
  the calculation is allowed.
- Correct the comment near 531–532 ("conservative for subgroup estimation"), which no longer
  describes either branch.
- **Tests:**
  - On an adjusted logistic fit, the standard error equals the delta-method value recomputed
    in the test by central finite differences of
    `mean(plogis(X1 %*% b)) - mean(plogis(X0 %*% b))` with respect to `b`, using a separately
    fitted `glm()`.
  - With no covariates the same calculation equals
    `sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)`.
  - Through `make_effect_estimator(outcome_type = "binary", effect_measure = "RD",
    adjust_covariates = ...)` on data that reaches the logistic tier: assert
    `method_used` first as the premise, then the standard error. If no deterministic dataset
    reaches that tier with covariates, test the calculation directly and say so.
- **Unaffected case:** unadjusted risk difference (the F13 values), odds ratio, risk ratio,
  and an adjusted risk difference that the identity-link tier fits, are identical.

### B2 — rate-difference branch (`.make_poisson_rate_estimator()`, `dr_gcomp`, near 856)

Reported to use `abs(avg_deriv) * se_treat` as well.

1. Write out from the code what the estimate is, as a function of which fitted coefficients,
   and how the current standard error is formed. Quote the lines.
2. **Only if** the estimate is a smooth function of the coefficients of one fitted outcome
   model, with nothing estimated elsewhere entering it (no propensity model, no augmentation
   term): replace the standard error by `sqrt(g' V g)`, with `g` the gradient of that exact
   expression with respect to every coefficient and `V` that fit's covariance matrix. Same
   safety rule, same kind of test (finite differences on the expression as documented), and
   the unaffected case is every other branch of that estimator.
3. **Otherwise:** change nothing and report B2 BLOCKED with the form found.

## Part C — printed labels, documentation, examples, comments, articles

### C1–C3 — printed labels (code, one commit each)

- **C1.** `R/format_subgroup_summary_tables.R` near 222 highlights `">= 0.95"`, with a space;
  the label the summary produces is `">=0.95"`, so the top consistency row is never
  highlighted. Use the label the summary produces. Test: every label in the highlight set is
  one of the labels of the consistency distribution.
- **C2.** The threshold remap warning (`R/forestsearch_main.R` near 2334–2350) quotes the
  risk-difference default and "OR scale" for a rate difference too, and the validator warning
  (near 301–309) says "A risk difference cannot exceed 1.0" for a rate difference. Make each
  message name the measure and the default that apply. The remapped values do not change.
  Test: the messages for a risk difference and for a rate difference carry the defaults given
  in `?forestsearch` ("Resolved defaults, per estimand").
- **C3.** The GRF footnote "(risk-difference scale)" (`R/summary_utility_functions.R` near
  917) is printed for counts as well. For a count outcome, name the scale of the outcome
  given to GRF, read from what the object records. Other outcome types are unchanged. If the
  object does not record enough to say, report C3 BLOCKED.

### C4 — help pages and source comments (one commit; comment-only gate)

| Site | Correction |
|---|---|
| `sg_tables()`, the CV summary tables, `plot_sg_glm_outcomes()`, the forest-plot display, `mr_estimates_table()` | State that the effect column follows the fit's orientation, while the arm means and their raw difference do not, so for a fit where a lower outcome is adverse the two have opposite signs. |
| `make_effect_estimator()` `@param adverse_outcome` | It says the argument is ignored for continuous outcomes; the linear-model estimators negate the outcome. Say what the code does for each outcome type. |
| `forestsearch()` roxygen near `R/forestsearch_main.R:1313-1315` | It says a supplied threshold "suppresses the corresponding remap and the supplied value is used". An out-of-range risk- or rate-difference value is remapped with a warning, and the remapped value is also the DINA floor. |
| `forestsearch()` near `R/forestsearch_main.R:440-445` | "The same cut `subgroup_method = "dina"` would select" holds for the default selection statistic only. |
| `forestsearch()` `@param consistency_method` and `R/subgroup_consistency_helpers.R` near 336–338 | They say propensity-adjusted effects fall back to `"split"`; the code decides by effect measure only. |
| `CV_sgs()` `@return` | Define "Exact" for a one-factor full-data subgroup: the fold's subgroup is that factor alone. |
| Bootstrap summary (`summarize_bootstrap_subgroups()`) | State that the bins are closed on the left. |
| `test_hte_from_forestsearch()` `@return` | Describe the mapping of cut labels to data columns and the two warnings. |
| `calibrate_k_treat()` near `R/calibrate_k_treat.R:9-11` | It calls `calibrate_k_inter()` the general-generator counterpart; that is `find_k_inter_for_target_hr()`. |
| `find_k_inter_for_target_hr()` | `@return` lists a `method` element the result lacks; `@param continuous_vars` says "standardized" where the default does not standardize. |
| `compute_dgm_cde()` `@param dgm` | It says the object must contain `df_super_rand`; the code also accepts `df_super`. |
| `generate_bootstrap_with_noise()` | "Preserves integer type" and "same structure as the input" (an integer column comes back double); "Missing values are not imputed" (a selected missing value in a categorical column is replaced by a level, and an ordered factor with more than two levels and a missing value stops); add that `continuous_vars` must be numeric. |
| `setup_gbsg_dgm()` | Add a Details paragraph carrying the model description that went with the removed `create_gbsg_dgm()` help page (subgroup definition, model, `k_inter`); recover it from `git show d93e6621:R/sim_aft_gbsg.R` and check it against the builder. |
| `#` comments | `R/dina_subgroup.R` near 402 ("log-IRR for poisson"); `R/forestsearch_main.R` near 2457–2458 ("event count threshold" for `dmin.grf`); `R/get_fsdata.R` near 327 (omits count). |

Run `devtools::document()`; each rewritten page passes `tools::checkRd()`.

### C5 — examples that cannot run as written (one commit; comment-only gate)

All are in `\dontrun`. Make each runnable as written in a fresh session after
`library(forestsearch)`: load the data with `gbsg <- survival::gbsg`, build any derived column
the example names (`grade3` is not a GBSG column), and build or drop any object the example
uses but never defines. Sites, from part 3's report §7.2:

- `R/bootstrap_summaries_helpers.R` near 28; `R/forestsearch_cross_validation.R` near 2002,
  2031; `R/plot_sg_distribution.R` near 76; `R/plot_sg_results.R` near 866, 912;
  `R/render_forestplot.R` near 372–375, 410–413; `R/run_subgroup_sims.R` near 660;
  `R/subgroup_consistency_helpers.R` near 1031, 1361, 1645;
  `R/fs_fdr_report.R` near 148–164; `R/summary_utility_functions.R` near 742–745;
  `R/pareto_frontier_table.R` near 125; and `summarize_simulation_results()`'s example, which
  passes an `nsim` argument that does not exist.
- Run each corrected example once in a fresh `Rscript --vanilla`, cap 2 minutes. Report each
  as RAN, or NOT RUN with the reason (over the cap, or it needs a long simulation). An
  example that errors is corrected until it runs or is reported.

### C6 — comments on test expectations (one commit; comment-only gate)

Add the derivation or clarifying comments the part 4 verifiers asked for (part 4 report,
OPEN ITEMS 17), in the test files of F1, F5, F6, F8, F9, F10, F12, F13 and F14. No expected
value, seed or tolerance changes. The affected test files still pass.

### C7 — two article sources (one commit; prose only)

- `vignettes/articles/biomarker_effects.qmd` near 66: the sentence names the removed
  `create_gbsg_dgm()`. Name `setup_gbsg_dgm()`.
- `vignettes/articles/paper_simulations.qmd` near 388–445 and 459–463: the prose calls an
  unevaluated `generate_aft_dgm_flex()` block "equivalent" to `setup_gbsg_dgm()`, and says
  `setup_gbsg_dgm()` "remains the correct entry point". Say instead that the block is a
  comparable construction with the general generator and is not numerically equivalent
  (state the difference the code shows), and that `setup_gbsg_dgm()` is used in this article
  to reproduce the published simulations and is superseded for new work.
- **Gate:** only these two files change under `vignettes/`; no changed line is inside a code
  chunk; no `.html` is edited. Nothing is rendered.

## After the parts

1. **Default-path gate.** Re-run the six equivalence cases with part 4's gate scripts
   (`dev/reports/REPORT_cran_prep_code_fixes_2026-11-12_gate_run.R` and `_gate_compare.R`) at
   the head against the start-SHA export. All six must be identical; an unexplained difference
   is FAILED for the task and is reported, not repaired.
2. **Articles.** Do not re-render. State, per code fix, whether any call an article evaluates
   reaches the changed lines and whether its output can change.
3. `DESCRIPTION`: next development version. `NEWS.md`: one entry per code fix that passed,
   saying what a user sees change, and one entry for the documentation pass. Commit by
   explicit path (`chore: version and NEWS for part 5`).
4. **Full suite** (`devtools::test()`, `NOT_CRAN=true`), cap 30 minutes. **Gate:** 0 failures,
   0 errors. Report blocks and expectations beside part 4's and account for the difference.
5. **`devtools::check(document = FALSE, error_on = "never")`** on the check route of
   `dev/reports/REPORT_cran_prep_2026-11-11.md` §6, cap 60 minutes. **Gate:** 0 errors,
   0 warnings, 0 notes.
6. Write `dev/reports/REPORT_cran_prep_displays_se_docs_2026-11-12.md`, in bullets and tables:
   - provenance and commits;
   - one status line per item (PASS, FAILED or BLOCKED) with its commit SHA;
   - per code fix: the site before and after, the red and green runs, how each expected value
     was derived, the unaffected case, and one before/after example of what a user sees;
   - A4: the full sweep table;
   - B2: the form found, and the reason if BLOCKED;
   - C4: old and new text per site, and any sentence changed from this task's wording;
   - C5: the RAN / NOT RUN table;
   - the default-path gate, the articles statement, suite and check with wall clocks;
   - OPEN ITEMS, limited to what this task itself leaves open.
   A detail file beside the report is allowed.
7. Commit the report by explicit path (`dev: CRAN preparation report, part 5`) and copy it,
   and any detail file, to `~/Downloads`.
8. **Final gates:** `git status --porcelain` is empty; both sibling repositories are
   unchanged; `git diff <start> --stat` lists only files the report names; under `vignettes/`
   only the two C7 files; nothing under `quarto/`. Nothing is pushed or merged.
