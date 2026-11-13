# REPORT — CRAN preparation, part 5: remaining displays, adjusted standard error, documentation, 2026-11-12

Task: `dev/tasks/TASK_cran_prep_displays_se_docs_2026-11-12.md` (committed as received at the start SHA). Run unattended from `~/Documents/GitHub/forestsearch` on branch `cran-prep`, after part 4 (`REPORT_cran_prep_code_fixes_2026-11-12.md`). Authorized by Larry on 2026-11-12; one fix, one commit.

**Result.**
- Code fixes A1, A2, A3, B1, C1, C2, C3: **7 PASS**, each in its own commit.
- A4 sweep: **PASS**, no further site qualified, so nothing changed and there is no commit.
- B2: **BLOCKED**. The rate-difference estimate depends on the estimated propensity score, so the task's pre-set condition does not hold; nothing changed.
- Documentation C4, C5, C6, C7: **4 PASS**, each in its own commit.
- 0 FAILED. The default-path gate passes (6/6 calls identical to the start SHA). Full suite: 0 failures, 0 errors. `check()`: 0 errors, 0 warnings, 0 notes. Final gates: PASS.
- Full per-item evidence: `REPORT_cran_prep_displays_se_docs_2026-11-12_detail.md` (beside this report).

## 1. Provenance

| Item | Value |
|---|---|
| Branch | `cran-prep`. The checkout was on `master` (also at `d35c2b8a`) when the task started; `git status --porcelain` was empty (gate), and Step 0 switched to `cran-prep`. `master` was not touched. |
| Part 4 gate | `dev/reports/REPORT_cran_prep_code_fixes_2026-11-12.md` exists on `cran-prep` (commit `d35c2b8a`). Its result: 14 PASS, 0 FAILED, 0 BLOCKED; suite 0 failures (703 blocks, 8,853 expectations); `check()` 0/0/0. |
| Start SHA | `b06de4e8` — `dev: add CRAN preparation task, part 5` (parent `d35c2b8a`). Exported with `git archive` (package paths: `R`, `man`, `tests`, `inst`, `DESCRIPTION`, `NAMESPACE`, `LICENSE`, `NEWS.md`, `.Rbuildignore`, `vignettes`, `_pkgdown.yml`) to `/tmp/cran_prep_displays_se_docs_2026-11-12/start/`. A first export of every path (6.1 GB with `quarto/`) was deleted and replaced by this one before anything used it. |
| `DESCRIPTION` version | 0.4.1.9010 at the start, as expected → 0.4.1.9011 |
| fs-post-selection | `main` at `5d165acce104a336e11093f1881cfdb572348462`, `git status --porcelain` empty, at the start and at the end |
| fs-glms-interpretable | `main` at `dc94a17e0e93653d5d4a9fb83bacca93bb8f87ea`, `git status --porcelain` empty, at the start and at the end |
| R / tools | R 4.6.2; devtools 2.5.2, roxygen2 8.1.1, testthat 3.3.2, pkgload 1.5.3; survival 3.8-12, grf 2.6.1, glmnet 5.1 |
| Scratch | `/tmp/cran_prep_displays_se_docs_2026-11-12/` (`$S` below): start export, per-item scratch copies `wt/`, evidence `evidence/<ID>/`, verifier copies `verify/`, gated patches `patches/`, main-session runs `main/`, gate runs `gate/`, after-the-fixes checks `cumulative/`, `articles/`, `examples_head/`, suite `suite/`, check `check/`, report drafts `report_sections/`. Nothing was installed into the user or site library. The temporary libraries are `$S/lib/C5`, `$S/lib/C5_verify` and `$S/lib/head`, used only to run examples after `library(forestsearch)`. Every scratch reproducer used `pkgload::load_all()`, a sequential plan and `timeout 120`; the main-session red/green runs used `timeout 300`. |

**Commits** (every `git add` named explicit paths; nothing fetched, pulled, pushed or merged):

| SHA | Time | Commit | Item |
|---|---|---|---|
| `b06de4e8` | 15:45 | `dev: add CRAN preparation task, part 5` | Step 0 (start SHA) |
| `2984dbee` | 16:49 | `fix: name the rate difference, IRR scale and 0.01 default in the IRD threshold warnings` | C2 |
| `f4c7dd2d` | 17:04 | `fix: sg_tables() uses the fit's own orientation` | A1 |
| `a2cd36db` | 17:05 | `fix: CV summary tables use the fit's own orientation` | A2 |
| `fb21396f` | 17:05 | `fix: compute_frontier_cis() falls back to forestsearch()'s orientation rule` | A3 |
| `8e87cf30` | 17:05 | `fix: covariate-adjusted risk-difference SE by the full delta method` | B1 |
| `9ec67b7e` | 17:05 | `fix: highlight the >=0.95 row of the bootstrap consistency table` | C1 |
| `c3323a59` | 17:05 | `fix: sg_tables() GRF CATE note names the count scale given to the forest` | C3 |
| `ba7c4ea2` | 17:06 | `docs: help pages and source comments say what the code does (part 5, C4)` | C4 |
| `5073fd67` | 17:07 | `docs: make sixteen \dontrun examples runnable as written (part 5, C5)` | C5 |
| `ffc55486` | 17:08 | `test: derivation comments on part 4's test expectations (part 5, C6)` | C6 |
| `835fac1b` | 17:08 | `docs: two article sources name setup_gbsg_dgm() and say what the flex block is (part 5, C7)` | C7 |
| `fb80dd8f` | 17:09 | `chore: version and NEWS for part 5` | After the parts, step 3 |
| `2a8be3bc` | 17:39 | `chore: NEWS wording for part 5` | NEWS review (§8) |
| (this commit) | | `dev: CRAN preparation report, part 5` | After the parts, step 7 |

C2 was committed first because it was the one code item with no comment polish pending (§1, "How the work was done"); the fixes are independent and every patch applied cleanly on top of the earlier commits.

**How the work was done.**
- **Implement (scratch, one agent per item).** Each item was built in a private copy of the start SHA (`$S/bin/mkcopy.sh`: the start export as a git repo) under `$S/wt/<ID>/`. For a code fix the agent wrote the test first, recorded red on the copy with only the test file, applied the fix, recorded green, ran the item's unaffected case (start export and fixed copy, each in its own `Rscript --vanilla` process, `identical()`), ran the existing test files that exercise the changed functions, and wrote a patch and `NOTES.md`. C4 was split by file into three agents (C4a, C4b, C4c) and committed as one.
- **Verify (one independent adversarial verifier per item).** It re-applied the patch to fresh copies, re-ran red and green, re-derived every expected value by hand, re-ran and probed the unaffected case, and checked scope and callers (for documentation: every sentence against the code, the comment-only gate, `document()` and `checkRd()`). All 15 were CONFIRMED. Three had a repair round first: A2 (three blocking text errors about count fits in the notes and proposed NEWS; no code change), A4 (two rows missing or misnamed in the sweep table) and B1 (an SE of 0 where the start gave NA on a slice where the treatment effect is not estimable; fixed so the old value is kept).
- **Polish (comments only, before commit).** The verifiers' minor findings on comments and help text were fixed in A1, A2, A3, B1, C1, C3, C4a, C4b, C4c, C6 and C7 by comment/roxygen/prose-only edits, each checked by an independent agent. Three of those checks found problems: an A1 code comment that overclaimed for counts, a B1 helper comment true only for a nonzero constant, and two problems in one C7 sentence about `simulate_from_dgm()`. The main session corrected them in the copies and regenerated the patches. The pre-polish patches are kept as `$S/patches/<ID>.verified.patch`; each polish delta is `$S/evidence/<ID>/polish.diff`.
- **Commit (main session, one item at a time).** For a code fix, `..._apply_fix.sh` applied only the patch's test file to the repository and ran it (red), then the code (green); `..._commit_fix.sh` checked that the staged diff equals the gated patch (ignoring `index` and hunk-header lines), printed `--stat` and committed. For C4 and C5 the `R/` hunks were applied, `devtools::document()` was run in the repository, and the regenerated `man/` was compared with the verified patches' `man/` (§5). The scripts are beside this report.
- **After the commits.** Read-only agents re-ran every code fix's unaffected case at the head with all fixes in (§7.2), traced the articles per code fix (§7.3), and re-ran the C5 examples against an install of the head (§7.4). Thirteen agents drafted the per-item sections from the evidence, checked against git; their NEWS checks and one further independent check led to the NEWS wording commit `2a8be3bc` (§8).
- **Measured wall clocks.** Implement and verify workflow (36 agents): 3,275 s. Polish and check (22 agents): 858 s. After-the-fixes checks (15 agents): 891 s. Report drafting (13 agents): 1,077 s. NEWS check (1 agent): 456 s. Main-session red/green per fix: 2.9–8.1 s (§3). `devtools::document()` in the repository for C4: 11.77 s (and 11.64 s for the no-op run before it); the C5 run's wall was not recorded. Default-path gate: 18.5 s (start) and 17.3 s (head). Suite 754 s; `check()` 1,103 s (§7.5). Three runs hit the 120 s cap: the `fs_fdr_report` example in the C5 implementer and verifier runs (exit 124 at 120.76 s; reported NOT RUN) and one superseded articles trace for A3 (120.05 s). Every other run finished under its cap.

## 2. Status per item

| Item | Status | Commit | One line |
|---|---|---|---|
| Step 0 | **PASS** | `b06de4e8` | Clean tree; part 4's report on `cran-prep`; task committed; start exported; version 0.4.1.9010; siblings recorded. |
| A1 | **PASS** | `f4c7dd2d` | `sg_tables()`'s two GLM closures pass the fit's `adverse_outcome` (forestsearch() fit: `args_call_all`; GRF result: its top-level value; rule if absent). |
| A2 | **PASS** | `a2cd36db` | `forestsearch_KfoldOut()`'s two GLM closures pass the fit's `adverse_outcome` from `cv_args` (rule if absent). |
| A3 | **PASS** | `fb21396f` | `compute_frontier_cis()`: the `%||% TRUE` fallback replaced by `forestsearch()`'s rule; a stored value is used as before. |
| A4 | **PASS** (no change) | — | Every `make_effect_estimator(` call listed (15 rows); no site beyond A1–A3 and part 4's F12 qualifies; other routes listed. |
| B1 | **PASS** | `8e87cf30` | Change of method: adjusted logistic-tier RD SE is `sqrt(g' V g)` over the estimable coefficients; old expression kept when that is unusable. |
| B2 | **BLOCKED** | — | The `dr_gcomp` IRD estimate uses the estimated propensity score (as a covariate and in the prediction), so the pre-set condition fails; nothing changed. |
| C1 | **PASS** | `9ec67b7e` | The top consistency row is highlighted: `">= 0.95"` → `">=0.95"`, the label the summary produces. |
| C2 | **PASS** | `2984dbee` | Both IRD warnings name the rate difference, IRR and (remap) the 0.01 it assigns; RD messages and remapped values unchanged. |
| C3 | **PASS** | `c3323a59` | The GRF CATE note names the count scale given to the forest, read from `grf_y_transform` and `offset.name`; binary and continuous unchanged. |
| C4 | **PASS** | `ba7c4ea2` | 14 help-page and comment sites rewritten from the code (18 `R/` files, 17 Rd pages); comment-only gate, `checkRd()` clean. |
| C5 | **PASS** | `5073fd67` | 16 `\dontrun{}` examples corrected; 15 RAN, `fs_fdr_report` NOT RUN (over the cap); comment-only gate, `checkRd()` clean. |
| C6 | **PASS** | `ffc55486` | 16 comment requests in 9 part-4 test files; comment-only gate; counts unchanged. |
| C7 | **PASS** | `835fac1b` | Two article sources, prose only, no changed line in a code chunk, nothing rendered. |
| Default-path gate | **PASS** | — | 6/6 calls identical to the start SHA (§7.1). |
| Version/NEWS | **PASS** | `fb80dd8f`, `2a8be3bc` | 0.4.1.9011; one NEWS entry per code fix that passed and one for the documentation pass; wording corrected after review. |
| Suite | **PASS** | — | 0 failures, 0 errors (§7.5). |
| `check()` | **PASS** | — | 0 errors, 0 warnings, 0 notes (§7.5). |
| Final gates | **PASS** | — | §10. |

## 3. Code fixes

Line numbers are at the start SHA `b06de4e8` and at head `fb80dd8f` (`R/` is identical at `835fac1b`, `fb80dd8f` and `2a8be3bc`). Each entry gives the site before and after, the red and green runs, how the expected values were derived, the unaffected case and one before/after example. The detail file has every run, every expected value and every verifier problem.

### A1 — `sg_tables()` uses the fit's own orientation

- **Status:** PASS — commit `f4c7dd2d` (parent `2984dbee`) — 2 files, +290/-0:
  - `R/summary_utility_functions.R` +22/-0;
  - `tests/testthat/test-sg-tables-orientation.R` new, +268.
  - No `man/` change. The commit equals `patches/A1.patch` (only the `index` lines differ, in abbreviation length).
- **Site, before → after:** two sites in `R/summary_utility_functions.R`.
  - GRF branch: `:828-839` at `b06de4e8` → `:850-873` at `fb80dd8f` (new argument `:869`).
  - Main branch: `:981-992` → `:1038-1059` (new argument `:1055`).
  - The task's pointers are exact at the start SHA.
  - Neither closure passed `adverse_outcome`, so `make_effect_estimator()`'s formal `TRUE` applied whatever the fit used. Each branch now resolves the value first and passes it.
  - The main branch reads `args_fs$adverse_outcome` (head `:1042-1045`).
  - The GRF branch reads `fs$adverse_outcome` (head `:856-859`). A `grf_glm_result` has no `args_call_all`; it stores the value at the top level (`R/grf_subg_harm_glm.R:551`, `:574`, `:647`, `:710`).
  - Before (start `:982-989`, main branch):
  ```r
    estimator_fn_tab <- tryCatch(
      make_effect_estimator(
        ...
        effect_measure    = effect_measure,
        adjust_covariates = args_fs$adjust_covariates
  ```
  - After (head `:1042-1056`):
  ```r
    adverse_outcome <- args_fs$adverse_outcome
    if (is.null(adverse_outcome)) {
      adverse_outcome <- outcome_type %in% c("binary", "count")
    }
    ...
        effect_measure    = effect_measure,
        adverse_outcome   = adverse_outcome,
        adjust_covariates = args_fs$adjust_covariates
  ```
- **Red → green:**
  - Main session, repository at `2984dbee` + the test file only: red, 41 failures of 155 expectations (counted from the progress line, since the reporter listed the first 10), 5.8 s.
    - First failure: `test-sg-tables-orientation.R:199:3`, "Expected MD est ITT < `5e-07 + 1e-12`. Actual comparison: 0.4573760 >= 0.0000005". 0.4573760 is 2 × 0.228688, the sign flip.
  - Green at `f4c7dd2d`: 155 passes, 0 failures, 5.5 s.
  - Implementer reproduced it: red 155/41 failed, 6.84 s; green 155/0, 5.85 s.
  - Verifier reproduced it: red 155/41, 6.83 s; green 155/0, 8.34 s.
  - Per block in every red run: 13/0, 24/11, 23/10, 46/18, 45/0, 4/2. The fixture block and the TRUE block pass on both trees.
- **Expected values:** derived by hand from the documentation; none comes from the fixed code.
  - **The effect cell.** The unadjusted OLS treatment coefficient is `mean(y | T=1) - mean(y | T=0)`. It is oriented by `?forestsearch` / `?grf.subg.harm.glm` `@param adverse_outcome` and the convention effect > 0 = harm, so FALSE gives `-(mean1 - mean0)`.
  - **The limits.** Estimate ± `qnorm(0.975)` × the pooled SE.
  - **Mean(C), Mean(T), Diff.** The raw values.
  - **Ties to the fit.**
    - `forestsearch()`: `max_sg_est` equals the hand Questionable value (1.0941122; measured difference 2.22e-16).
    - GRF result: the Questionable sign equals the sign of `-cate_sg` (0.5838).
  - **Where the comments are:** the derivation in the header `:21-35`, the tolerances in `:36-48`, `.sgo_md()` at `:102-111`, and the `max_sg_est` ties at `:166-176` and `:244-247`.
  - **Re-derivation.** The verifier re-derived every value with its own data generator. Max |table − hand| is 4.3e-7 for estimates and 1.1e-5 for limits.
- **Unaffected case:** the task's case is survival fits, and non-survival fits with `adverse_outcome = TRUE`. Display data (`_data` of every gt table, at `ndecimals` 3 and 6) is `identical()` in every run:
  - **Implementer:** 13/13, start vs fixed. The cases cover survival; continuous, binary OR/RD and count TRUE; binary and count with the value unset; count FALSE; and 5 GRF results. Walls: start 9.51 s, fixed 9.68 s, start control 9.47 s.
  - **Whole objects** are also identical in 13/13 once gt closure environments are dropped. Those environments differ between any two runs, start vs start included (192 leaves).
  - **Verifier:** 15/15 of its own cases plus 3 notation runs. Walls: 13.98 s / 13.26 s / 13.68 s.
- **Existing tests:**
  - `test-no-subgroup-dina-grf-summary.R` is the only existing caller: 21 expectations, 0 failed. Implementer 5.74 s, verifier 5.68 s.
  - Context run of part 4's `test-display-orientation.R`: 93/0 (7.07 s / 6.91 s).
  - Expectations updated: none.
- **What a user sees:** a continuous `forestsearch()` fit at the default (`max_sg_est = 1.09411`), `sg_tables(fit)$tab_estimates` (`evidence/A1/example.txt`):
  ```
                Mean(C) Mean(T) Diff     MD (95% CI) before        MD (95% CI) after
  ITT           0.003   -0.226  -0.229   -0.229 (-0.464, 0.007)    0.229 (-0.007, 0.464)
  Questionable  -0.330  -1.424  -1.094   -1.094 (-1.665, -0.523)   1.094 (0.523, 1.665)
  Recommend     0.045   -0.051  -0.096   -0.096 (-0.338, 0.146)    0.096 (-0.146, 0.338)
  ```
  - The arm columns and Diff do not change.
  - A binary `grf.subg.harm.glm()` result at its default FALSE also changes: Questionable OR 0.540 (0.288, 1.010) → 1.853 (0.990, 3.468), the OR of the non-event (`binary_false_example.txt`).
- **Verification:** CONFIRMED, no repair round. There were seven minor problems:
  - **1–2 (polish).** The reasons for the 1e-10 tolerance on `max_sg_est` and for the `+ 1e-12` slack were added as comments (test `:166-175`, `:244-246`, `:38-45`). Red and green were unchanged at 155/41 and 155/0.
  - **4 (polish, then main-session correction).** "On the fit's own scale" in the R comment became "orientation". The polish check found that wording wrong for counts, so the main session added that count and binary IRR/IRD do not change (head `:850-855`).
  - **3 (left as is).** An off-by-one line number in NOTES.
  - **5 (open item).** The read-only quarto ACTG175 continuous template (`:182`) would flip if re-rendered.
  - **6 (open item).** A pre-existing DINA stop in `sg_tables()` (head `:1175`), identical on both trees.
  - **7 (dropped).** "Diff has the opposite sign" was loose for OR/RR, and the final NEWS bullet drops it.
- **Departures:**
  - **The GRF branch reads the result's top-level `adverse_outcome`, not `args_call_all`.** A GRF result has none. The task's rule says only "the fit's stored `adverse_outcome`"; `args_call_all` is F12's pattern.
  - **Its fallback is `forestsearch()`'s rule, not `grf.subg.harm.glm()`'s FALSE default.** Only objects not built by `grf.subg.harm.glm()` reach it.
  - **The GRF test ties the sign to `-cate_sg`.** A GRF result has no `max_sg_est`.
  - **Extra tests and cases:** fallback, TRUE and mirror-twin tests, and an extended unaffected case.
  - **Binary OR/RR/RD fits that store FALSE also change.** The rule intends this.
- Full evidence: detail file, section A1.

### A2 — CV summary tables use the fit's own orientation

- **Status:** PASS — commit `a2cd36db` (parent `f4c7dd2d`) — `R/forestsearch_cross_validation.R` +12/-0, new `tests/testthat/test-cv-tables-orientation.R` +234 (2 files, +246/-0). No `man/` change. The commit equals `patches/A2.patch` (the polished patch).
- **Site, before → after:** `forestsearch_KfoldOut()`, the two `make_effect_estimator()` calls inside `if (outall)`: `R/forestsearch_cross_validation.R:1509-1517` (no-original-subgroup ITT branch) and `:1565-1575` (main GLM branch) at `b06de4e8` → `:1533-1542` (new line `:1540`) and `:1590-1601` (new line `:1597`) at `fb80dd8f`, plus a new block `:1512-1520`. The task's pointers are exact at the start SHA. Head lines are 14 lower than at the commit because of C4's roxygen (`ba7c4ea2`); the code lines are unchanged. Before, both calls used `make_effect_estimator()`'s formal `TRUE`.

Before (start `:1570-1572`; site 1 is the same):
```r
          offset.name       = offset.name,
          effect_measure    = effect_measure,
          adjust_covariates = res$cv_args$adjust_covariates
```
After (head `:1517-1520`, `:1596-1598`):
```r
    adverse_outcome <- res$cv_args$adverse_outcome
    if (is_glm && is.null(adverse_outcome)) {
      adverse_outcome <- outcome_type %in% c("binary", "count")
    }
...
          effect_measure    = effect_measure,
          adverse_outcome   = adverse_outcome,
          adjust_covariates = res$cv_args$adjust_covariates
```

- **Red → green:** Main session, repository at `f4c7dd2d` plus the test file only: 111 expectations, 25 failed, 86 passed (10 listed, then "... and 15 more"). First: `test-cv-tables-orientation.R:183:3` via `.cvo_check_row()` `:124`, "Expected MD cell ITT < `.cvo_tol`. Actual comparison: 0.457375922 >= 0.000000005" (7.7 s). Green with the code: 111 passed, 0 failed (7.3 s).
  - Implementer reproduced it: round 1 red 25 of 111 (8.34 s), green 0 of 111 (7.82 s); repair round 8.53 s / 7.90 s; polish 8.37 s / 7.79 s. Verifier: round 1 8.43 s / 7.66 s, round 2 7.92 s / 7.36 s, same counts; polish check green 0 of 111 (7.68 s).
  - Failures per block: fixture 0/6, main 11/29, no subgroup 2/11, rule fallback 12/33, `TRUE` 0/32. Every failure is a sign mismatch (difference 2|est|).
- **Expected values:** By hand; no number is a literal in the test. `.cvo_md()` (`:103-110`) gives the OLS coefficient `mean(y|T=1) - mean(y|T=0)` with pooled SE, negated under FALSE.
  - Orientation from `?forestsearch` `@param adverse_outcome` (`R/forestsearch_main.R:1104-1107` start); limits `est -/+ 1.96 se` from `analyze_subgroup_glm()` (`R/summary_utility_functions.R:223`). Working in the header (`:19-42`).
  - Fixture values: ITT 0.22868796 (-0.00700868, 0.46438460), H 1.09411222 (0.52314346, 1.66508098), Hc 0.09579320, KfA H 0.80136279, KfA Hc 0.03182661. `Diff` stays raw (`:126-128`). Hand H equals `max_sg_est` within 1e-10 (`:176`, `:230`; reason `:36-42`); observed 2.2e-16.
  - Tolerance `5e-9 + 1e-12` (`:46`, reason `:34-36`): tables built at `digits = 8`. Largest observed deviation 4.922e-09.
  - Re-derived three ways, equal to every printed digit: implementer base R; verifier r1 `t.test(var.equal = TRUE)`; verifier r2 normal equations without the package.
- **Unaffected case:** Task: survival fits, and non-survival fits with `adverse_outcome = TRUE`. Each tree ran `forestsearch()` → `forestsearch_Kfold(Kfolds = 3)` → `KfoldOut(outall = TRUE)`, `(TRUE, digits = 8)` and `(FALSE)` in its own process; only `timing_minutes` was dropped.
  - Implementer, 11 cases (survival with and without a subgroup; continuous TRUE, also adjusted and without a subgroup; binary OR rule, RD explicit TRUE, no subgroup, stored value removed; count rule): `identical(start, fixed)` TRUE for every case and overall. Walls: start 18.10 s, fixed 17.42 s, start control 17.81 s, compare 0.18 s (repair round 17.03 s / 17.31 s).
  - Verifiers added 20 cases (31.85 s / 31.96 s) and 13 cases (tenfold, Kfolds 5, `cv_summary_tables()`; 25.11 s / 25.46 s): 44 configurations, all identical. Count at explicit FALSE is identical too: the rate estimator takes no `adverse_outcome`.
- **Existing tests:** 12 files on the fixed copy and a start control, one process each: cv-exact-match 6, cross-outcome-parity 49 (2 skipped, 12 warnings), cv-no-subgroup-dina-grf 48, cv-fold-na-assignments 6, cv-no-subgroup-edges 45 (1 skipped), cv-mixed-and-fatal 58, downstream-consumers 43 (5 warnings), input-validation 19 (1 warning), rng-state-restored 16, missing-membership-messages 60, no-subgroup-cv 72, display-orientation 93. All have 0 failed and 0 errors, with the same counts on both trees. `test-mr-inference.R` (whole file, fixed only): 148 passed (181.61 s). Expectations updated: none.
- **What a user sees:** continuous fit at the default (`adverse_outcome` FALSE, `max_sg_est` 1.09411), `forestsearch_KfoldOut(forestsearch_Kfold(fit, Kfolds = 3), outall = TRUE)$tab_all`, default digits (`evidence/A2/example.txt`):
```
before (start tree)                          after (fixed tree)
Overall  -0.2287 (-0.4644, 0.0070)           Overall 0.2287 (-0.0070, 0.4644)
FA_0    -1.0941 (-1.6651, -0.5231)           FA_0     1.0941 (0.5231, 1.6651)
KfA_0   -0.8014 (-1.2420, -0.3608)           KfA_0    0.8014 (0.3608, 1.2420)
```
  The selected subgroup (`FA_0`) now reads as the fit's 1.0941. `Diff` stays -1.0941 on both trees. In a verifier probe, a binary OR fit at FALSE has its H cell move from 0.60000000 to 1.66666667, the hand OR on `1 - y`.
- **Verification:** CONFIRMED in round 2, after one repair round.
  - Round 1 PROBLEMS: three blocking problems, all text about count fits in NOTES and the proposed NEWS (count tables cannot change). Repair: text corrected; no code change; the test header's tie-tolerance comment went from 2 lines to 7.
  - Round 2 had two minor problems. (1) The test comment `:38` omitted that the search negates y under FALSE: fixed by the polish. (2) Stale round-1 line pointers in NOTES: not fixed; the detail file uses head numbers.
  - Polish: 7 comment lines rewrapped at `:36-42`; `R/` byte-identical. Polish check: CONFIRMED, no problems.
- **Departures:**
  - The fallback is resolved once, before the branch, for both sites. It is guarded by `is_glm`, so the survival path never reads it.
  - The test covers every GLM table, the no-subgroup branch (deterministic, `effect.threshold = 10`), the rule fallback with the stored value removed, and a `TRUE` guard. This is more than the task's minimum.
  - Not tested: binary OR, RD or RR with explicit FALSE, which changes. The task names only the continuous default; the verifier's probes cover these.
  - Unchanged and noted: count (IRR, IRD) and binary IRR/IRD are never oriented. Count tables print `IRR` "NA" because `resCV` lacks the offset (pre-existing).
- Full evidence: detail file, section A2.

### A3 — `compute_frontier_cis()` falls back to `forestsearch()`'s orientation rule

- **Status:** PASS — commit `fb21396f` (parent `a2cd36db`) — `R/frontier_cis.R` +9/-1, `tests/testthat/test-frontier-cis-orientation.R` (new) +194; 2 files, +203/-1. No roxygen, `man/`, `NEWS.md` or `DESCRIPTION` change. The commit equals `patches/A3.patch` (`main/A3/patch_equal.txt` empty).
- **Site, before → after:** `R/frontier_cis.R:150` at `b06de4e8` → new block `:140-147` and the call argument `:158` at `fb80dd8f`. With no stored value, the fallback is now `forestsearch()`'s rule (`b06de4e8:R/forestsearch_main.R:1840-1842`, head `:1877-1879`) on the resolved `outcome_type` (`R/frontier_cis.R:136`), not `TRUE`. A stored value is read by the same `args$adverse_outcome`.

Before
```r
      adverse_outcome   = args$adverse_outcome %||% TRUE,
```
After (without the four-line comment `:140-143`)
```r
  adverse_outcome <- args$adverse_outcome
  if (is.null(adverse_outcome)) {
    adverse_outcome <- outcome_type %in% c("binary", "count")
  }
...
      adverse_outcome   = adverse_outcome,
```

  - Only continuous output changes: the rule gives FALSE there and `.make_lm_estimator` negates `y`. Binary and count keep TRUE; survival gets FALSE, which the Cox estimator does not read.
  - Every `forestsearch()` fit stores the resolved value (resolved at start `R/forestsearch_main.R:1840-1842`, head `:1877-1879`, before the capture at start `:1923-1924`, head `:1960-1961`), so only an object built or edited by hand reaches the new branch.
- **Red → green:** main session, repository at `a2cd36db` with the test file only: 18 of 31 not passed, 13 passed (the 10 listed are all "Failure", then "... and 8 more"). First failure `test-frontier-cis-orientation.R:147:3`, "Expected `row$estimate` to equal `ref$est`." `actual` -1.1, `expected` 1.1; 3.8 s. Green with the code: 31 passed, 3.5 s.
  - Implementer reproduced it: 18 failed of 31 (per test 0/5/12/1) → 0 failed; 3.96 s → 3.91 s.
  - Verifier: the same counts; 4.17 s → 3.83 s.
  - After the polish: implementer 18 → 0 (4.15 s → 3.86 s); polish check 18 → 0 (4.17 s → 3.83 s).
  - Every failure is a sign mismatch (implementer, verifier).
- **Expected values:** by hand, from the data and the documentation.
  - MD = `-(mean(y | T=1) - mean(y | T=0))` under the continuous default `adverse_outcome = FALSE` (`?forestsearch` `@param adverse_outcome`).
  - Pooled-variance SE.
  - Naive CI `est -/+ qnorm(0.975) se`; the Split-CI midpoint is the estimate (`?compute_frontier_cis`).
  - Frontier membership is evaluated by hand from the labels.
  - The derivations are in the header `:23-38` and beside `.fco_md()` (`:86-93`) and `.fco_mask()` (`:95-110`). The tolerance `1e-10` and its reason are at `:40-42`.
  - Selected subgroup (m = 1, n = 36, 19/17): est 1.09411222 (equal to `max_sg_est`, `:144`), se 0.29131059, naive (0.52315395, 1.66507049). Rows m = 3 and m = 5: 0.78471467 and 0.58595462.
  - Checks at head: `:147-152` and `:173-180`.
  - The verifier rederived them with its own label parser in base R: hand − `max_sg_est` = -2.2e-16; fixed table within 1.8e-15. It found no value taken from `compute_frontier_cis()`.
- **Unaffected case:** task: "any object that stores the value gives identical output". `compute_frontier_cis()` table, printed verbose text and warnings, start vs fixed in separate processes, `identical()` in a third. No timing leaf; nothing dropped.
  - Implementer, 11 cases: survival; continuous FALSE/TRUE; binary OR TRUE, RD TRUE, OR FALSE; count IRR; survival, binary and count with the entry removed; continuous with the entry removed (intended). 10 of 10 non-intended cases TRUE; start vs start control 11 of 11. Walls 16.26 / 15.59 / 15.82 s.
  - Verifier, 20 cases, adding `survival::gbsg` (10/10 finite rows), odd stored values (NA, 0, 1, `"TRUE"`, `logical(0)`), RD FALSE, RR, IRR FALSE and IRD: all stored cases TRUE (45.66 / 43.29 s). Seven try-to-break objects: six identical; the seventh (continuous, no `effect_measure`, entry removed) is the intended sign change.
  - The only difference is the intended one: continuous, entry removed, every estimate negated. On the fixed tree that table is `identical()` to the stored-FALSE table.
  - After-the-fixes agent at head (`cumulative/A3/`): all non-intended cases identical to start, and whole outputs identical to the A3-only tree (implementer's script: start 16.92 s, head 16.68 s).
- **Existing tests:** none call `compute_frontier_cis()` or its caller `compare_selection_rules()` (`git grep` over `tests/` at `b06de4e8`: none).
  - Neighbour checks: `test-display-orientation.R` 93/0 (implementer 6.86 s; verifier 7.20 s) and `test-sg-focus-transparency.R` 64/0 (verifier 8.24 s).
  - Full suite at `fb80dd8f`: `✔ | 31 | frontier-cis-orientation`.
  - Expectations updated: none.
- **What a user sees:** continuous fit at the default orientation (`max_sg_est = 1.09411`, `sg.harm = {age <= 47} & {biomarker_hi}`), entry removed, `compute_frontier_cis(fit, n_splits = 200L, seed = 1L)` (`evidence/A3/example.txt`):
```
start tree (row m = 1; the Key and column-type lines are left out):
       m estimate naive_lcl naive_ucl split_lcl split_ucl fsbc_estimate
1:     1  -1.0941   -1.6651   -0.5232   -1.7340   -0.4542       -1.0925
fixed tree:
       m estimate naive_lcl naive_ucl split_lcl split_ucl fsbc_estimate
1:     1   1.0941    0.5232    1.6651    0.4542    1.7340        1.0925
```
  The selected subgroup's row was the negative of the fit's own estimate; it now equals it. Rows m = 3 and 5 flip the same way.
- **Verification:** CONFIRMED; no repair round. One polish round, comments only (+6/-2 in the test file; check CONFIRMED, no problems). Four minor problems:
  1. The test header said "blocks" (plural) for the one commented resolution block. Reworded in the polish (`:10-12`).
  2. The last expectation is relational, not a hand value. A comment saying so was added (`:188-190`).
  3. The implementer's survival case had 1 of 3 finite rows (pre-existing `.label_to_expr()` issue). Closed by the verifier's `survival::gbsg` case; the issue is an open item.
  4. Binary RR errors "no 'restart' 'muffleWarning' found". Pre-existing A-056; no change.
  - One verifier statement does not match git: it says fits saved before `bc4690b3` stored `adverse_outcome = NULL`. Every version since `99ced431` writes the resolved value back (`bc4690b3^:R/forestsearch_main.R:693`), and older objects have no frontier. So no object saved by any `forestsearch()` version reaches the new branch; the test comment `:19-20` ("saved by an older version") is inexact (detail A3.13).
- **Departures:**
  - Every `forestsearch()` fit stores the value, so a default fit would pass on both trees. The test therefore removes the stored entry, after asserting that `FALSE` was stored.
  - The fallback is a block before the call, as in F12.
  - The test checks every frontier row, the Split-CI midpoint and one relational `identical()`, besides the selected subgroup.
  - The unaffected case is extended beyond the task's wording.
  - No roxygen change: `?compute_frontier_cis` still does not say how the estimate is oriented (open item).
  - The task's site and wording match the code.
- Full evidence: detail file, section A3.

### B1 — Risk-difference standard error by the full delta method, covariate-adjusted logistic tier (change of method)

- **Status:** PASS — commit `8e87cf30` (parent `fb21396f`) — `R/glm_effect_estimators.R` +75/-4, `tests/testthat/test-rd-adjusted-se.R` +303 (new); 2 files, +378/-4.
  - No `man/` change: the new helper is `@noRd`, and `document()` wrote nothing.
  - The commit equals the final `patches/B1.patch` apart from the `index` lines.
  - After B1, C4 (`ba7c4ea2`) rewrote roxygen above the site, so B1's lines are 18 lower at head than at `8e87cf30`.
- **Site, before → after:** `.estimate_rd()`, tier 2 (logistic g-computation), `R/glm_effect_estimators.R`.
  - At `b06de4e8`: comment `:531-532`, adjusted branch `:545-547`.
  - At `fb80dd8f`: comment `:549-553`, adjusted branch `:566-578`, new helper `.fs_rd_gcomp_se()` `:611-666`, file-header helper list `:21`.
  - At the commit: `:531-535`, `:548-560`, `:593-648`.

  Before (`b06de4e8:545-547`):
  ```r
    } else {
      abs(avg_deriv) * se_treat
    }
  ```
  After (`fb80dd8f:566-578`, comment lines omitted):
  ```r
    } else {
      se_dm <- .fs_rd_gcomp_se(fit_logit, nd0, nd1, pred0, pred1)
      if (length(se_dm) == 1L && is.finite(se_dm) && se_dm >= 0) {
        se_dm
      } else {
        abs(avg_deriv) * se_treat
      }
    }
  ```
  - What the helper computes:
    - `g = colMeans(p1 (1-p1) X1) - colMeans(p0 (1-p0) X0)` and `sqrt(g' V g)` on the estimable (non-NA) coefficients, with `V = vcov(fit)`.
    - `X0`, `X1` are built as `predict.lm()` builds them.
  - When it returns `NA_real_`:
    - the value is not finite and non-negative;
    - an aliased column differs between `X1` and `X0` (`:654-655`);
    - any step errors or warns.

    The caller then keeps the old expression.
  - The old comment ("this is conservative for subgroup estimation") is replaced.
- **Red → green:** main session, on the parent `fb21396f` plus the test file only:
  - red: 7 passes, 5 failures, 4 errors, 3.0 s. First failure `test-rd-adjusted-se.R:111:3`: "Expected `res$se` to equal `.rdadj_fd_se(d, ~treat + site + x)`", `actual` 0.1291, `expected` 0.1217. The 4 errors are `could not find function ".fs_rd_gcomp_se"` (`:161`, `:243`, `:294`) and `local_mocked_bindings` (`:232`).
  - green: 35 passes, 0 failures, 2.9 s.
  - Reproduced on the final test file:
    - implementer repair round: red 12 expectations, 5 failed, 4 errors, 3.00 s; green 35/0, 2.82 s;
    - verifier round 2: red 12/5/4, 2.92 s; green 35/0, 2.87 s;
    - polish: red 12/5/4, 3.02 s; green 35/0, 3.06 s.
  - The new not-estimable block fails 4 times on the round-1 code: implementer 2.97 s, verifier 2.89 s.
  - Round 1, on the 254-line test file: implementer red 11/5/3 (3.30 s), green 27/0 (2.81 s); verifier the same counts (3.14 s, 2.88 s).
- **Expected values:** all in `tests/testthat/test-rd-adjusted-se.R`, with the working in comments. Each tolerance was set before the first run.
  - (a) `:111`, factor `site` + `x`, through `make_effect_estimator()`: the delta-method SE recomputed by central finite differences (h = 1e-5) of `mean(plogis(X1 b)) - mean(plogis(X0 b))` on a separately fitted `glm()`. Tolerance 1e-6 relative, with the truncation and rounding bounds in the comment (`:58-74`). (a2) `:124`, `:127`: the same value with an aliased `x2 = 2x`.
  - (b) `:161`, helper on an unadjusted fit, 12/40 vs 15/30: `sqrt(0.3*0.7/40 + 0.5*0.5/30) = 0.1165476` by hand, from `g'Vg` written out for the saturated model (`:146-160`). `round(., 5) = 0.11655` (`:163`).
  - (c) `:207`, `:210`, `:211`, treated arm 40/40 with binary `z`: `SE = sqrt(0.375^2*0.2*0.8/20 + 0.625^2*0.5*0.5/20) = 0.0775101` and estimate `0.6125`, by hand (`:186-206`).
  - Safety `:235`: the old expression on a separate `glm()`, with the helper mocked to NA, NaN, Inf and -0.5; tolerance 1e-12 (`:220-230`).
  - Helper failures `:243-253`: `NA_real_` by contract.
  - Not estimable `:294`, `:301`: NA by definition (`:281-289`).
  - `method_used == "logistic_gcomp_adj"` is asserted first as the premise (`:108`, `:123`, `:183`, `:234`, `:251`, `:300`).
  - Measured agreement: (a) 1.85e-11, (b) 1.87e-07, (c) 3.82e-09.
  - Both verifiers re-derived (a)-(c) with their own Newton-Raphson MLE, without `glm()` or `vcov()`.
- **Unaffected case:** task: unadjusted RD (the F13 values), OR, RR, and adjusted RD fitted by the identity-link tier. Each tree ran in its own process, compared with `identical()` in a third; timing leaves and `row.names` were dropped.
  - Implementer: 260 `make_effect_estimator()` results + 3 `forestsearch()`. Unaffected 216 TRUE (unadjusted RD 48 with F13, OR/RR 132, identity-link adjusted RD 8; also `iptw_rd`, `dr_gcomp` RD and 20 pre-existing errors). The 44 `logistic_gcomp_adj` results are identical with `se` removed. `forestsearch()` x3 TRUE. Walls 5.18 s / 5.32 s / 0.27 s; repair re-run 5.13 / 5.20 / 0.26 s.
  - Verifier round 1: 960 results, 852 TRUE, 8.42 / 8.19 / 0.22 s.
  - Verifier round 2: 1360 results + 50 displays + 4 `forestsearch()`. Unaffected 1044 TRUE (unadjusted RD 68, OR 340, RR 340, identity-link adjusted RD 132; also `iptw_rd`, `dr_gcomp` and 48 pre-existing `muffleWarning` errors). The 56 not-estimable results are identical in full. Walls 8.21 / 8.86 / 0.29 s.
  - Only an adjusted RD that reaches tier 2 changes `forestsearch()` output: the candidate `L(HR)`/`U(HR)` columns, with the same rows and the same `sg.harm`.
- **Existing tests:** run start vs fixed on 15, 17 and 28 files.
  - Final round (verifier round 2): fixed 28 files, 2390 expectations; start 27 files (no new file), 2355; the difference is the new file's 35. Per-file counts of the 27 common files are identical; 0 failed, 0 errors, 8 skipped, 16 warnings on both; 103.8 s / 104.0 s.
  - Per-file counts are identical start vs fixed. The warnings come from `cross-outcome-parity` and `get-fsdata-diagnostics` on both trees.
  - Expectations updated: none.
- **What a user sees** (`evidence/B1/example.txt`; walls 2.42 s / 2.59 s):
  ```
  (1) control 14/40, treated 40/40, adjust_covariates = "z"
  start:  se 6.79127662230129e-06   95% Wald CI: 0.6125 (0.6125, 0.6125)
  fixed:  se 0.0775100802854232     95% Wald CI: 0.6125 (0.4606, 0.7644)
  (2) control 10/30, treated 16/30, c("site", "x"), SG_tab_estimates_glm() ITT "RD (95% CI)"
  start:  0.200 (-0.053, 0.453)
  fixed:  0.200 (-0.039, 0.439)
  ```
  - With the control arm at 0% the old SE blew up instead (verifier round 2, `f13_ctl0_trt80` + `age`): `0.8044 (-830.7797, 832.3885)` → `0.8044 (0.6838, 0.9251)`.
  - The estimate is unchanged.
- **Verification:** CONFIRMED after one repair round.
  - Round 1 found one blocking problem. An aliased treatment coefficient (one arm empty, or all of its outcomes NA) was dropped like a covariate, so `g = 0` and the SE became 0 where start gave NA: cell `0.000 (0.000, 0.000)`, `subgroup_glm()` `c(0, 0)`.
  - Repair: the helper returns NA when an aliased column differs between `X1` and `X0`, and a new test block (`:256-303`) covers it. Round 2 confirmed it.
  - Polish (comments only, gate empty):
    - the helper roxygen was reworded for an aliased treatment-interaction column (verifier round 2 minor 1);
    - `.fs_rd_gcomp_se` was added to the header list `:21` (minor 2).
  - The polish check found "with age constant on the slice" false for a constant 0. The main session corrected it to "constant and nonzero" (head `:626`).
  - Minor problems left without a change:
    - test (b) has about 5× headroom (1.87e-07 vs 1e-6);
    - test (b)'s comment and redundant rounding assertion were edited after the first green run (disclosed; tolerance unchanged);
    - three NOTES sentences were not corrected: the stale §2 quote; "closure assertions pass on the start tree", which is true by the verifier's start run; and the bootstrap callers row (see departures and open items in the detail).
- **Departures:**
  - The safety rule is wider than the task's "not a finite, non-negative number". The helper also returns NA, and the old expression is kept, on any error or warning in the helper and when a column that treatment changes is aliased. Without that, a finite 0 replaced start's NA for a not-estimable effect.
  - The helper takes the fitted `pred0`/`pred1`, not `plogis(X b)`, because `binomial()` clamps |eta| > 30.
  - "Estimable columns only" means NA coefficients dropped by name.
  - Three extra test blocks (mocked safety rule, helper failures, not estimable).
  - The header-list edit at `:21` is outside the function.
  - The approved formula is the delta method conditional on the sample's covariates. Recorded, not changed.
- Full evidence: detail file, section B1.

### C1 — The consistency table highlights its `">=0.95"` row

- **Status:** PASS — commit `9ec67b7e` (parent `8e87cf30`) — `R/format_subgroup_summary_tables.R` +4/-2, new `tests/testthat/test-consistency-highlight-labels.R` +121 (2 files, +125/-2). The commit equals `patches/C1.patch`.
- **Site, before → after:** `format_subgroup_summary_tables()`, Table 4 "Consistency (P<sub>cons</sub>) Distribution", green `tab_style()`: `R/format_subgroup_summary_tables.R:222` at `b06de4e8` → `:224` at `fb80dd8f` (comment `:215` → `:215-217`). The label the summary produces is `">=0.95"` (`R/summarize_bootstrap_subgroups.R:420` at `b06de4e8`, `:432` at `fb80dd8f`), so the top row matched no highlight. The yellow and red sets already matched (`:232`, `:242` → `:234`, `:244`).

Before:
```r
            rows = `Consistency Range` %in% c("0.85-0.95", ">= 0.95")
```
After:
```r
        # Color coding: High consistency (green).  The labels in the three
        # highlight sets are the Pcons bin labels of
        # summarize_bootstrap_subgroups() (no space in ">=0.95").
...
            rows = `Consistency Range` %in% c("0.85-0.95", ">=0.95")
```

- **Red → green:** Main session, repository at `8e87cf30` plus the test file only: red `.123.45`, 2 passes (premises) and 5 failures, all on the `">=0.95"` row (`:98, :103, :106, :119, :120`). First: `:98` `actual[3:6] "#FFF3CD" "#FFF3CD" "#D4EDDA" ""` vs `expected[3:6] ... "#D4EDDA"` (3.1 s). Green with the code: 7 passes, 0 failures (3.0 s).
  - Implementer and verifier reproduced it on the verified file (lines 5 lower): red 2 passes / 5 failures (3.27 s; verifier 3.41 s), green 7 passes (3.00 s; verifier 3.34 s). After the polish: red 2/5 at the head lines (3.32 s), green 7/7 (3.26 s; polish check 3.36 s).
- **Expected values:** All by hand, none from the fixed code.
  - Row labels (`:92-95`): the documented left-closed bins; a premise.
  - Inputs: one Pcons value per bin (`:84-86`); five values >= 0.95 give Count `0 0 0 0 0 5` (`:110-112`, `:114`).
  - Fill per row, red, red, yellow, yellow, green, green (`:98-101`): read off the table's own footnote "Green (>=0.85), Yellow (0.7-0.85), Red (<0.7)".
  - Colours (`:44-46`): the formatter's `cell_fill()` literals, which the fix does not touch, upper-cased as gt stores them. Bold on green rows only (`:103`, `:120`).
  - Working in the header (`:13-42`) and beside each check (`:84-86, :91, :97, :102, :105, :110-112, :118`). The verifier re-derived every value with base R only (0.15 s).
- **Unaffected case:** The task names none for C1. The implementer's case: everything in the formatted output other than the `">=0.95"` row's highlight is identical.
  - Implementer: four inputs, 4 summaries and 30 gt tables; `ALL IDENTICAL ... TRUE`. Tables were compared with `identical(ignore.environment = TRUE)` and by rendered HTML. Two things were relaxed: gt's formatter closures, which differ even start vs start, and the fixed tree's 3 green + bold `_styles` rows on that row. The HTML differs on lines 480-482 only. Walls: start 3.17 s, fixed 3.01 s, compare 3.79 s; start-vs-start control TRUE.
  - Verifier: eight other inputs (breaks, NA, out-of-range, none found, 120 rows); `TRUE`. Walls 3.53 s / 3.32 s / 5.24 s.
- **Existing tests:** No existing test or package function calls the formatter. `test-bootstrap-summary-bins.R` (the label producer) was run anyway: 8/8 with the new file's 7/7, 0 failures (implementer 3.33 s, verifier 3.40 s; start copy 8/8, 3.12 s). Expectations updated: none.
- **What a user sees:** Seven bootstraps, Pcons 0.95, 0.96, 0.98, 0.99, 1.00, 0.90, 0.88 (`evidence/C1/example_before_after.txt`):
```
before (start tree)                               after (fixed tree)
         0.85-0.95     2   28.57 #D4EDDA bold              0.85-0.95     2   28.57 #D4EDDA bold
            >=0.95     5   71.43       (none)                 >=0.95     5   71.43 #D4EDDA bold
```
  The row that holds most bootstraps is now green and bold, as the footnote says. The four lower rows are red or yellow in both trees.
- **Verification:** CONFIRMED in round 1; no repair round. One polish, comments only: the test header no longer claims the tests prove every highlight label is a summary label. It now says what they check and what they do not (`:33-41`). The polish check: CONFIRMED, 4 mutants. Minor problems:
  - (1) The header claimed more than the tests check: a mutant with an extra `">= 0.95"` in the green set still passes 7/7. Fixed by the polish.
  - (2) NOTES says the CV Pcons table "has no highlight". It has one on summary rows (`R/summarize_cv_results.R:1179`), but no bin-label highlight. Not changed; no code effect.
  - (3) The test reads gt's internal `_styles` (gt 1.4.0). Left as an open item.
- **Departures:**
  - The task names no unaffected case; the implementer defined one.
  - The test checks that every distribution row carries its colour, the converse of the task's "every label in the highlight set is one of the labels". An extra unmatched highlight label would pass, and the header says so. The task's property holds by reading the code (the three sets name exactly the six labels), not by a test.
  - Three comment lines were added above the call.
  - The commit subject shortens the implementer's proposal.
  - Reach: the formatter is internal and called only in read-only `quarto/resampling/test_bootstrap_grf_cuts.qmd:314, 401`. No article output can change.
- Full evidence: detail file, section C1.

### C2 — Threshold warnings name the measure and the default that apply (rate difference)

- **Status:** PASS — commit `2984dbee` (parent `b06de4e8`) — `R/forestsearch_main.R` +22/-7, new `tests/testthat/test-threshold-remap-messages.R` +145 (2 files, +167/-7). No `man/` change. The commit equals `patches/C2.patch`: only the abbreviated `index` hashes differ, and the blobs are the same.
- **Site, before → after:** two sites, both in `R/forestsearch_main.R`.
  - Validator check 4, `.validate_outcome_threshold_config()`: `:301-310` at `b06de4e8` → `:301-317` at `fb80dd8f`.
  - Warned remap in `forestsearch()`: `:2334-2350` → `:2371-2395`. C4 later moved this block by +30 lines without touching it.
  - The remap message passed the literals `"OR"` and `0.05` while it assigned `switch(effect_measure, RD = 0.05, IRD = 0.01)`. It now prints the value it assigns. The validator now picks its text by measure. The RD text and the assigned values are unchanged.

  Before (start `:2346-2349`, remap):
  ```r
            effect_threshold, "OR", 0.05, effect_measure),
          call. = FALSE)
        effect_threshold <- switch(effect_measure,
          RD = 0.05, IRD = 0.01)
  ```
  After (head `:2375-2376`, `:2391`, `:2394`, remap; `:305-308`, validator):
  ```r
        .fs_remap_ratio   <- switch(effect_measure, RD = "OR", IRD = "IRR")
        .fs_remap_default <- switch(effect_measure, RD = 0.05, IRD = 0.01)
            effect_threshold, .fs_remap_ratio, .fs_remap_default,
        effect_threshold <- .fs_remap_default
    bound_line <- switch(effect_measure,
      RD  = "  A risk difference cannot exceed 1.0.\n",
      IRD = "  A rate-difference threshold above 1.0 is outside that range.\n")
    ratio_measure <- switch(effect_measure, RD = "OR", IRD = "IRR")
  ```
- **Red → green:**
  - **Main session.** Repository at `b06de4e8` plus the test file only: red, 18 pass, 9 fail, 8.1 s. The first failure is `test-threshold-remap-messages.R:84:3`, Expected `val` to match regexp "rate[- ]difference"; the actual text was "A risk difference cannot exceed 1.0. / Did you mean effect_measure = 'OR'". With the code: green, 27 pass, 0 fail, 7.7 s.
  - All nine failures are IRD text checks: validator `:84-87`, remap `:92-95` and `:97`.
  - The implementer reproduced it (red 18 pass / 9 fail, 8.18 s; green 27 pass, 7.90 s), and so did the verifier (red 18 / 9, 8.44 s; green 27, 8.12 s).
- **Expected values:**
  - **From the documentation.** The defaults RD 0.05 and IRD 0.01 come from `?forestsearch`, "Resolved defaults, per estimand" (`man/forestsearch.Rd:1207-1208` at the start, `:1222-1223` at head). The ratio counterparts OR (for RD) and IRR (for IRD) come from `@param effect_measure`.
  - **By hand.** The two full RD texts were written from the start format strings (`:302-308`, `:2335-2347`), with 1.5 printed by `%.2f` as "1.50". The in-range cases (RD 0.07, IRD 0.02, no message) follow from the condition `> 1.0`.
  - **Where the comments are.** In the test file: header `:4-17`; IRD block `:81-82`, `:89-90`, `:96`, `:99-100`; RD block `:112-116`, `:131`; in-range block `:137-138`.
  - **The verifier's re-derivation** used base R only (`v_derive.txt`): the defaults parsed from the start Rd are identical to the test's, and its hand-built RD texts equal the test's (`TRUE`, `TRUE`). No tolerance is involved.
- **Unaffected case:** the task's case is "The remapped values do not change". Each tree ran in its own process, and only the timing leaves `$find.grps$time_search` and `$minutes_all` were dropped.
  - **Implementer:** six `forestsearch()` cases, start export against the fixed copy: RD 1.5, RD 2, IRD 1.5, IRD 3, RD 0.07 and IRD 0.02, on synthetic N = 200 data.
    - Every fit and every `threshold_config` is `identical()` TRUE. Screening is 0.05 / 0.01 / 0.07 / 0.02 on both trees.
    - The RD warnings are identical. The IRD warnings differ only in the two C2 messages.
    - Walls: 23.12 s / 23.16 s, compare 0.19 s.
  - **Verifier:** 19 fit cases (legacy `hr.threshold`, unset, 1.0 boundary, DINA, IRD on binary, survival): every fit is identical, as are screening and `args_call_all`. Walls: 49.50 s / 49.51 s.
    - Validator grid, 336 cells: 328 identical; the 8 that differ are exactly {IRD, > 1, user-set}. Walls: 2.41 s / 2.61 s.
    - The 175-cell threshold-resolution probe is identical. Walls: 2.86 s / 2.83 s.
- **Existing tests:** 0 failures in every run.
  - Implementer, on the fixed copy:
    - Batch 1, 25.08 s; the start control has the same counts, 24.74 s: binary-default-or-entry-points, binary-default-or, dina-gaussian-floor, threshold-pair-directive-a and threshold-sync.
    - Batch 2, 48.29 s: declaration-c0, directive-c, estimability-boundary, fs_dgm_scale, fs-dgm-feasibility, lean-fits and rd-standard-error.
    - Three skips, for `dev/` and quarto files missing from the export. With the baselines placed, directive-c gives 182/0/0 and threshold-pair-directive-a 86/0/0, on both the fixed copy (12.26 s) and the start control (12.38 s).
  - Verifier: the same files plus sg-focus-transparency, with no skips (32.20 s and 45.97 s).
  - `helper-threshold-sync.R` and `test-dina-gaussian-floor.R` reach the edited remap.
  - Expectations updated: none.
- **What a user sees:** `effect_measure = "IRD"`, `effect.threshold = 1.5`, count data (`evidence/C2/unaffected/compare.txt`). The unchanged first line of each warning is left out. Both trees record 0.01.
  ```
  # start
    A risk difference cannot exceed 1.0.
    Did you mean effect_measure = 'OR' (ratio scale)?
    Did you mean effect.threshold = 1.50 on the OR scale?
    Remapping to default: 0.05.
    To use a custom IRD threshold, pass a value in [-1, 1] (e.g., effect.threshold = 0.07 for a 7 pct-point threshold).
  # fixed
    A rate-difference threshold above 1.0 is outside that range.
    Did you mean effect_measure = 'IRR' (ratio scale)?
    Did you mean effect.threshold = 1.50 on the IRR scale?
    Remapping to default: 0.01.
    To use a custom IRD threshold, pass a value in [-1, 1] (e.g., effect.threshold = 0.02 for a rate difference of 0.02).
  ```
- **Verification:** CONFIRMED in round 1. There was no repair round and no polish (C2 is not in `polish_result.json`). Seven minor problems:
  - The validator names no default: a stated departure (below).
  - The implementer's proposed NEWS wording gave 0.01 to both warnings: not carried into `fb80dd8f`. The committed bullet said the validator names "the IRR scale" (it suggests `effect_measure = 'IRR'`) and left out the changed example line; both corrected in `2a8be3bc`.
  - The hunk count and two off-by-one pointers in NOTES: corrected in the detail file.
  - Left open: the start comments "bounded in [-1, 1]" (head `:296`, `:2369`), which contradict the new comment at `:302-304`; the IRD test asserts fragments, not full strings; and with `outcome_type = "survival"` and a supplied IRD, the validator suggests `'IRR'`, a pre-existing kind of oddness.
- **Departures:**
  - **The validator message names no default.** It also fires on survival with a supplied RD or IRD, where no remap follows and a default would be false; the verifier ran fits that showed this (screening log(1.5) = 0.4054651). The remap message, which follows on every GLM call, names the default.
  - **The IRD line says "A rate-difference threshold above 1.0 is outside that range"**, not "cannot exceed 1.0". A rate difference is not bounded by 1: IRD = `lam1 - lam0` (start `R/glm_effect_estimators.R:896-905`).
  - **The remap's IRD example changed** from "0.07 for a 7 pct-point threshold" to "0.02 for a rate difference of 0.02". The task does not mention it.
  - **The assignment reads the same local the message prints.** The value is unchanged.
- Full evidence: detail file, section C2.

### C3 — `sg_tables()` GRF CATE note names the count scale given to the forest

- **Status:** PASS — commit `c3323a59` (parent `9ec67b7e`) — `R/summary_utility_functions.R` +25/-2, new `tests/testthat/test-grf-footnote-scale.R` +201 (2 files, +226/-2). No `man/` change. The commit's changed lines equal `patches/C3.patch` (after polish).
- **Site, before → after:** `sg_tables()`, `grf_glm_result` branch, `R/summary_utility_functions.R:912-921` at `b06de4e8` → `:946-978` at `fb80dd8f`. The pre-set condition holds: every return path of `grf.subg.harm.glm()` records `grf_y_transform` and `offset.name` (start `R/grf_subg_harm_glm.R:546-551`, `:569-574`, `:640-647`, `:704-710`). For a count the forest gets the per-subject count, as `log(Y + 0.5)` (default) or unchanged (`"identity"`), never divided by exposure (`:401`, `:410-412`, `:421-423`).

  Before (start `:917`):
  ```r
            "**GRF CATE:** H&#x302; = %.4f, H&#x302;<sup>c</sup> = %.4f (risk-difference scale)",
  ```
  After (head `:956-968`, `:974`; 8 comment lines `:948-955` omitted):
  ```r
      cate_scale <- "risk-difference scale"
      if (identical(outcome_type, "count")) {
        count_scales <- c("log(Y + 0.5)" = "log(Y + 0.5) scale, Y = count",
                          "identity"     = "count scale")
        y_tr <- fs$grf_y_transform
        if (is.character(y_tr) && length(y_tr) == 1L &&
            y_tr %in% names(count_scales)) {
          cate_scale <- count_scales[[y_tr]]
          if (!is.null(offset.name)) {
            cate_scale <- paste0(cate_scale, ", not offset-adjusted")
          }
        } else {
          cate_scale <- "scale not recorded"
        ...
            "**GRF CATE:** H&#x302; = %.4f, H&#x302;<sup>c</sup> = %.4f (%s)",
  ```
- **Red → green:** Main session, repository at `9ec67b7e` + the test file only: red `....12...3...4....56...`, 23 expectations, 17 pass, 6 failures, 0 errors. First: `:127:3`, `"... = -1.0322, ... = 1.0631 (risk-difference scale)"` vs `"... (log(Y + 0.5) scale, Y = count, not offset-adjusted)"` (5.3 s). Green with the code: 23 pass (5.0 s).
  - Reproduced: implementer red 17 pass / 6 failures (5.47 s), green 23/23 (5.62 s); verifier red 6 failures of 23 (5.63 s), green 0 failures (5.73 s); after polish, red 17 / 6 (5.30 s), green 23/23 (5.18 s; polish check 5.23 s). The binary test passes on every red tree.
- **Expected values:** By hand from the code at the start SHA (test header `fb80dd8f:tests/testthat/test-grf-footnote-scale.R:9-44`, comments beside the asserts). Exact strings, no tolerance.
  - `log(Y + 0.5) scale, Y = count, not offset-adjusted` (`:129`, and `:144` with `adverse_outcome = TRUE`: negation runs after `grf_y_transform` is set and does not change the scale); `count scale, not offset-adjusted` (`:157`); `log(Y + 0.5) scale, Y = count` without an offset (`:173`); `scale not recorded` with the transform removed (`:182`); binary `risk-difference scale` (`:200`).
  - Leading text and `%.4f` copied from the start format string (`:917`). The numbers are `fit$cate_sg` / `fit$cate_sgc` from `grf.subg.harm.glm()`, which the fix does not touch.
  - The verifier checked the printed numbers sit on the named scale: log fit -1.0322 / 1.0631 vs sample differences in `log(Y + 0.5)` of -1.0108 / 1.0502 (population, by hand, -0.9316 / 1.0960); identity fit -2.1554 / 5.8495 vs count differences -2.1160 / 5.7920. A risk difference of -2.1554 is impossible.
- **Unaffected case:** Task: "Other outcome types are unchanged."
  - Implementer: 8 binary and continuous `grf.subg.harm.glm()` cases (tree, adverse, frontier, harm notation, no subgroup), `identical()` TRUE for the whole `sg_tables()` list and the HTML; count with no subgroup identical; 4 count cases differ only in `_source_notes`; a survival `forestsearch()` fit identical (Table 1 whole, `sg10_out` without its `_formats` closures, which also differ start vs start2). Walls start 8.01 s, start2 8.07 s, fixed 7.95 s; survival 5.41 / 5.29 / 5.49 s.
  - Verifier, own data: 13 binary and continuous cases identical (whole list and HTML), count no-subgroup identical, 8 count cases differ only in `_source_notes`, survival identical as above. Walls 9.46 s start, 9.94 s start2, 9.79 s fixed (`/usr/bin/time`); compare 1.54 s.
  - At head (`cumulative/C3/RESULT.md`): of 37 cases, 9 identical, 16 differ only through A1 (`_data`, effect column, `adverse_outcome = FALSE` fits), 12 only through C3 (`_source_notes`). With A1 and C3 reverse-applied, all 37 are identical to start.
- **Existing tests:** `test-no-subgroup-dina-grf-summary.R` (the only one that calls `sg_tables()`) 21/21, `test-grf-count-orientation.R` 14/14, `test-glm-pipeline-integration.R` 10/10, `test-grf-inclusion-band.R` 2875/2875, 0 failures (implementer 5.44 / 4.40 / 5.88 / 9.33 s; verifier 5.36 / 4.51 / 5.81 / 9.48 s). No test referenced the note. Expectations updated: none.
- **What a user sees:** count fit, `effect_measure = "log_IRR"`, `offset.name = "person_time"`, source note of `sg_tables(fit)$tab_estimates` (`evidence/C3/unaffected/run_start.txt`, `run_fixed.txt`):
  ```
  start: **GRF CATE:** H&#x302; = -1.0322, H&#x302;<sup>c</sup> = 1.0631 (risk-difference scale)
  fixed: **GRF CATE:** H&#x302; = -1.0322, H&#x302;<sup>c</sup> = 1.0631 (log(Y + 0.5) scale, Y = count, not offset-adjusted)
  ```
  With `grf_count_transform = "identity"`: `(count scale, not offset-adjusted)`. Binary and continuous notes are unchanged. No package article reaches this note.
- **Verification:** CONFIRMED; no repair round. Four minor problems:
  - NOTES gave the R file as +27/-2 (numstat +25/-2): evidence prose, not edited; this report uses +25/-2.
  - The test header said `offset.name` is used only by the GLM estimates; it also sets the complete-case rows given to the forest (start `:372-373`): fixed in the polish (test `:16-21`).
  - The note does not state orientation (with `adverse_outcome = TRUE` the values are on `-log(Y + 0.5)` or `-Y`): left open, as for binary and continuous.
  - The R comment's version boundary "before 0.4.1.9010" was loose (`f2e8de47` already negates under 0.4.1.9009): fixed in the polish (head `:954` names commit `f2e8de47`).
  - Polish: comment-only, 2 files +8/-4; polish check CONFIRMED with no problems and two wording notes, not acted on.
- **Departures:**
  - `, not offset-adjusted` is added when an offset is recorded, so the CATE is not read as a rate difference.
  - `scale not recorded` for a count object without a recognised `grf_y_transform` (hand-built or modified objects only).
  - The sign is not named: negation changes orientation, not scale, and the object does not record whether a count was negated (objects from before `f2e8de47` carry `adverse_outcome = TRUE` un-negated). Implementer and verifier decided this independently.
  - The task says the label is wrong "for counts as well"; it is also wrong for continuous outcomes (a mean difference). Left unchanged by the task's rule; open item.
- Full evidence: detail file, section C3.

## 4. The sweep (A4) and the blocked item (B2)

### A4 — sweep of make_effect_estimator() call sites

- **Status:** PASS (no change). No commit, no patch, no test, no NEWS entry (`patches/` has no `A4*` file, there is no `main/A4/`, and `NEWS.md` at `fb80dd8f` has no A4 bullet). No call site beyond A1, A2, A3 and part 4's F12 builds a fit's effect without the fit's orientation.
- **Verification:** round 1 PROBLEMS. The call-site table was correct, but the other-routes list had 2 blocking problems (one route missing, one wrong function name) and 3 minor ones. One repair round changed only the sweep record (`evidence/A4/sweep.md`; the round-1 copy is kept as `sweep.r1.md`). Round 2 CONFIRMED, with 6 minor problems (wording and line pointers). The detail file applies them.
- **Sweep method:**
  - `git grep -n "make_effect_estimator(" b06de4e8 -- R/` finds 12 calls outside comments. The name has 40 matches in all: the 12 calls, the definition `R/glm_effect_estimators.R:100`, and 27 comment lines.
  - It also finds 3 calls of the bit-identical fast twins `.make_lm_estimator_fast(` and `.make_glm_or_estimator_fast(`. That makes 15 sites.
  - No call goes through `do.call`, `match.fun` or `get`. The enclosing functions were found by parsing the files (`verify/enclosing.txt`, `verify_r2/enclosing_r2.txt`).
- **At head `fb80dd8f`:** the same 15 sites, and no new call. The definition is at `:118`. The 4 extra matches of the name (44 in all) are roxygen lines added by C4 `ba7c4ea2`. Parsing each call's argument list at both SHAs shows that `adverse_outcome` is missing at 4 sites at the start (rows 2, 3, 14, 15) and present at all 15 at head.
- Line numbers below: start `b06de4e8` → head `fb80dd8f` (`R/` is identical at `835fac1b` and `fb80dd8f`).

| # | file:line, start → head | function | passes `adverse_outcome`? (source; start → head lines) | purpose | action |
|---|---|---|---|---|---|
| 1 | `R/bootstrap_dofuture_main.R:310` → `:310` (file unchanged) | `forestsearch_bootstrap_dofuture()` | **yes:** the fit's `args_call_all` (`args_forestsearch_call <- fs.est$args_call_all`, `:268`), with `forestsearch()`'s rule if NULL (`:306-309`), passed at `:317` | resampling: bootstrap estimates for the bias correction. The same closure builds the displayed bootstrap table (`SG_tab_estimates_glm()`, `:714-719`) | none (already passes) |
| 2 | `R/forestsearch_cross_validation.R:1510` → `:1534` | `forestsearch_KfoldOut()`, branch with no original subgroup | start: **no** (the estimator's default `TRUE`). Head: **yes**, `res$cv_args$adverse_outcome` (the fit's `args_call_all`), rule if NULL (`:1517-1520`), passed at `:1540` | display: ITT table (`SG_tab_estimates_glm()`, `:1518-1523` → `:1543-1548`) | fixed by A2 (`a2cd36db`) |
| 3 | `R/forestsearch_cross_validation.R:1566` → `:1591` | `forestsearch_KfoldOut()`, main GLM branch | start: **no**. Head: **yes**, the same resolved value (`:1517-1520`), passed at `:1597` | display: ITT, original-subgroup and K-fold tables (`:1583`, `:1602`, `:1624` → `:1609`, `:1628`, `:1650`) | fixed by A2 (`a2cd36db`) |
| 4 | `R/forestsearch_main.R:2263` → `:2300` | `forestsearch()` | **yes:** the resolved formal (`:1840-1841`, again at `:2241-2243` → `:1877-1878`, `:2278-2280`), passed at `:2270` → `:2307` | fitting: search and consistency closure | none |
| 5 | `R/forestsearch_main.R:2287` → `:2324` (`.make_lm_estimator_fast`) | `forestsearch()` | **yes**, passed at `:2290` → `:2327` | fitting: fast continuous twin | none |
| 6 | `R/forestsearch_main.R:2297` → `:2334` (`.make_glm_or_estimator_fast`) | `forestsearch()` | **yes**, passed at `:2300` → `:2337` | fitting: fast binary-OR twin | none |
| 7 | `R/forestsearch_main.R:2748` → `:2798` | `forestsearch()`, propensity rebuild | **yes**, passed at `:2754` → `:2804` | fitting: PS-adjusted closure | none |
| 8 | `R/frontier_cis.R:143` → `:151` | `compute_frontier_cis()` | start: **yes, with fallback `TRUE`**: `args$adverse_outcome %\|\|% TRUE` (`:150`). Head: `args$adverse_outcome`, with `forestsearch()`'s rule if NULL (`:144-147`), passed at `:158` | display: the frontier CI table, which is the `ci_table` input of `pareto_frontier_table()` and `plot_pareto_frontier()` | fixed by A3 (`fb21396f`) |
| 9 | `R/mr_estimates_table.R:211` → `:229` | `.mr_table_from_fs()`, called by `mr_estimates_table()` when `boot_results = NULL` (`:82` → `:100`) | **yes:** `a$adverse_outcome` with `a <- fs$args_call_all` (`:190` → `:208`), rule if NULL (`:207-210` → `:225-228`), passed at `:215` → `:233` | display | none (F12, part 4) |
| 10 | `R/plot_sg_glm_outcomes.R:262` → `:277` | `plot_sg_glm_outcomes()` | **yes:** `fs.est$args_call_all$adverse_outcome`, rule if NULL (`:169-172` → `:184-187`), passed at `:268` → `:283` | display | none (F12) |
| 11 | `R/plot_subgroup_results_forestplot.R:444` → `:468` | local `create_glm_row()` (`:404-485` → `:428-509`) inside `plot_subgroup_results_forestplot()` | **yes:** `fs.est$args_call_all$adverse_outcome`, rule if NULL, and `TRUE` when `fs.est = NULL` (no fit) (`:316-322` → `:340-346`), passed at `:450` → `:474` | display | none (F12) |
| 12 | `R/run_subgroup_sims.R:321` → `:321` | `subgroup_glm()` | **yes:** its own argument `adverse_outcome` (default `TRUE`, `:256`), passed at `:327`; same lines at head | other: a simulation fitter factory. Its input is a DGM sample, not a `forestsearch()` fit. Its only package caller (`:698` → `:720`) passes DGM fields | none (not a fit) |
| 13 | `R/run_subgroup_sims.R:336` → `:336` (`.make_lm_estimator_fast`) | `subgroup_glm()` | **yes**, passed at `:339` | other: same as row 12 | none (not a fit) |
| 14 | `R/summary_utility_functions.R:830` → `:863` | `sg_tables()`, GRF branch (`fs` of class `grf_glm_result`, `:811` → `:833`) | start: **no**. Head: **yes**, `fs$adverse_outcome` (the GRF result's top-level value), rule if NULL (`:856-859`), passed at `:869` | display: table 1 of a GRF result | fixed by A1 (`f4c7dd2d`) |
| 15 | `R/summary_utility_functions.R:983` → `:1049` | `sg_tables()`, main GLM branch | start: **no**. Head: **yes**, `args_fs$adverse_outcome`, rule if NULL (`:1042-1045`), passed at `:1055` | display: table 1 (ITT and subgroup estimates) | fixed by A1 (`f4c7dd2d`) |

- **Estimator layer.** Inside `make_effect_estimator()`, `.make_cox_estimator()` (`:142` → `:160`) and `.make_poisson_rate_estimator()` (`:163`, `:204` → `:181`, `:222`) receive no `adverse_outcome`. So the orientation changes only continuous MD and binary OR, RD and RR. The trace confirmed this: TRUE and FALSE give the same estimate for HR, IRR and IRD, and different estimates for MD, OR, RD and RR (`evidence/A4/trace_displays.txt`).
- **Dynamic check (start tree):** `trace_displays.R` traced `make_effect_estimator()` while 16 exported displays and entry points were run on the F12 continuous fit (stored `adverse_outcome = FALSE`): 14 completed; `plot.forestsearch` and `fs_declaration_table` stopped with an error (detail file). Process wall 10.00 s.
  - Only `sg_tables()` (A1) and `forestsearch_KfoldOut(outall = TRUE)` (A2) reached the factory without `adverse_outcome`. There the default `TRUE` applied.
  - The F12 displays, `compute_frontier_cis()` (A3, value stored, so the fallback was not reached), the bootstrap, and the refits all passed `FALSE`.
- **Other routes to a displayed effect (listed, not changed).** These build an effect without `make_effect_estimator()` carrying the orientation. Lines are at the start SHA. The full table is in the detail file.
  - **The one route that takes a `forestsearch()` fit and has no orientation** is `.extract_fs_estimates_gen()` (`R/run_simulation_analysis.R:963`; the file is unchanged at head).
    - It is reached from the exported `run_simulation_analysis()` (`:389`) through `.run_fs_analysis_gen()` (`:774`; fit at `:820`, passed as `fs_full` at `:846`, in the call at `:843`).
    - Its GLM `hr.itt`, `hr.H.hat`, `hr.Hc.hat`, `hr.H.true` and `hr.Hc.true` are a raw `stats::glm(y ~ treat)` coefficient (`.glm_effect()`, `:1005-1022`). The fit's `adverse_outcome` is never read.
    - The simulation summaries, `build_estimation_table()` and `plot_effect_distribution()` show these values.
    - Measured through the real caller on the F12 continuous fit: `hr.H.hat` = -1.094112, while the fit's own H estimate is +1.094112 (`repair_route_e2e.txt`, 3.29 s). Round 2 confirmed this on its own data: the route gives the same value under `adverse_outcome = FALSE` and `TRUE` (5.07 s).
    - It is not changed, for two reasons. It does not call `make_effect_estimator()`. And the simulation tables use one convention throughout (the truth columns use the same function; `compute_aor()` and `compute_cde_glm()` compare arm 1 with arm 0).
  - **GRF twin:** `.extract_grf_estimates_gen()` (`:1360`) uses the same raw coefficient. Its input is a GRF result, not a `forestsearch()` fit.
  - **No estimator of their own:** `analyze_subgroup_glm()` / `SG_tab_estimates_glm()`. A NULL closure gives the effect cell `"NA"`, and `Diff` is the raw arm difference.
  - **A NULL closure means the Cox path:** `bootstrap_results()`, `subgroup.search()`, `subgroup.consistency()`, `evaluate_subgroup_consistency()` and `evaluate_consistency_twostage()`.
  - **Fit-time inference with its own flip:** `.consistency_glm_pieces()`. Every spec that `forestsearch()` builds carries the resolved value.
  - **Not a `forestsearch()` fit:**
    - exported `create_glm_row()` (`R/grf_subg_harm_glm.R:760-866`), which has no package caller;
    - `grf.subg.harm.glm()` (`:284`, GRF fitting with its own orientation);
    - `.fs_region_effect()`, `fs_mr_oc_summary()`, `glm_effect_profile()`, `glm_cs_fit()` and `dina_subgroup_refit()`.
  - **Survival only, so no orientation applies:** `plot_sg_results()` → `compute_sg_hr_estimates()`, `fs_to_guohe()`, `fs_fdr_report()`, and the MRCT helpers (`.mrct_analyze_trial()`, `.mrct_record_fit()`).
  - **Raw arm contrasts of a hypothesis test:** `test_hte_from_forestsearch()`, which uses the fit only for covariate names.
  - **Values stored at fit time, no estimator call (trace):** `print`/`summary` of the fit, `mr_values()`, `pareto_frontier_table()`, `plot_pareto_frontier()`, `explain_pareto_selection()` and `frontier_member_flags()`; `fs_declaration_table()` by reading the code (no estimator call in `R/fs_declaration_table.R`; under the trace it stopped before any work).
  - **Refits that pass `args_call_all`:** `forestsearch_Kfold()`, `forestsearch_tenfold()` and the bootstrap.
- **Left open by A4 (listed, not changed):**
  - The simulation route's orientation (above).
  - The simulation scale labels (`R/run_simulation_analysis.R:995-1003`, `:1392-1399`): the value labelled RD is the log OR, IRD is the log IRR, and RR is the OR. Measured: `hr.H.hat` (RD) = -0.200671 = the hand log OR; the hand RD is -0.049536 (`verify/sim_rd_scale.txt`).
  - `fs_to_guohe()` has no `outcome_type` guard.
- **Differences between the task and the code:**
  - The task's "near" pointers hold. The calls are at `:830` and `:983` (A1; the `tryCatch(` lines are `:829` and `:982`), `:1510` and `:1566` (A2), and `:143` with the fallback at `:150` (A3).
  - The task's A4 bullet (task lines 90-94) asks for the call-site list and for fixes only. The evidence calls the other-routes list the item's "point 3". That point and the example routes quoted in `NOTES.md` §8 are not in the task document; they came from the run's item prompt, which is not saved under `$S`.
- Full evidence: detail file, section A4.

### B2 — rate-difference branch (dr_gcomp)

- **Status:** BLOCKED. No commit. Nothing changed: no patch (`patches/B2.patch` does not exist), no test, no NEWS entry. The task's pre-set condition (Part B, B2 step 2) does not hold, so step 3 applies.
- **Verification:** CONFIRMED by the verifier (round 1). No blocking problem, one minor problem (below). No repair round.
- **Site:** `R/glm_effect_estimators.R:849-862` at the start SHA `b06de4e8` = `:938-951` at head `fb80dd8f`. It is the `dr_gcomp` IRD branch of `.make_poisson_rate_estimator()`. The task's pointer "near 856" is the SE line itself (`:856`; head `:945`).
- **Same code at both SHAs:** the whole function (`:724-927` at start, `:813-1016` at head) is byte-identical (`diff` empty). It moved down 89 lines because of edits above it: B1 `8e87cf30` (+75/-4) and C4 `ba7c4ea2` (+23/-5).
- **How the branch is reached:**
  - Through `make_effect_estimator()` with `effect_measure = "IRD"` and `ps_adjust_method = "dr_gcomp"`. Both `outcome_type = "count"` (`:198-213`; head `:216-231`) and `outcome_type = "binary"` (`:162-172`; head `:180-190`) route there.
  - The slice must carry both an `ips_covar` and a `ps_hat` column (`:812-814`; head `:901-903`). Without them, the call falls through to the unadjusted branch (`poisson_offset_delta`, measured).
  - In `forestsearch()`, this needs `ps_method != "none"` and `ps_adjust_method = "dr_gcomp"`. The default is `"none"`. The closure is rebuilt at `R/forestsearch_main.R:2747-2756` (head `:2797-2806`).

**The form found.** Line numbers are at the start SHA. For head, add 89 in `R/glm_effect_estimators.R` and 50 in `R/forestsearch_main.R`.

```
:815  fmla <- stats::as.formula(
:816    paste0(outcome.name, " ~ ", treat.name,
:817           " + ips_covar + offset(.log_offset)")
:819  fit <- stats::glm(
:822    family = stats::poisson(link = "log")
:827  ps <- data_slice$ps_hat
:828  nd1 <- nd0 <- data_slice
:829  nd1[[treat.name]] <- 1L
:830  nd1$ips_covar     <- 1 / ps
:831  nd0[[treat.name]] <- 0L
:832  nd0$ips_covar     <- 1 / (1 - ps)
:834  pred1 <- stats::predict(fit, newdata = nd1, type = "response")
:835  pred0 <- stats::predict(fit, newdata = nd0, type = "response")
:851  est <- mean(pred1) - mean(pred0)
:852  se_treat <- sqrt(diag(stats::vcov(fit)))[[treat.name]]
:853  avg_deriv <- mean(pred1)
:856  se          = abs(avg_deriv) * se_treat,
```

- **Estimate.** Let `b = (b0, bT, bI)` be the coefficients of the single Poisson fit `glm(y ~ treat + ips_covar + offset(log t))` on the slice. Then, averaged over all `n` rows of the slice (both arms):

  ```
  pred1_i = t_i * exp(b0 + bT + bI / ps_i)
  pred0_i = t_i * exp(b0      + bI / (1 - ps_i))
  est     = mean(pred1) - mean(pred0)
  ```

| Quantity | Where it is estimated |
|---|---|
| `b0`, `bT`, `bI` | the slice's Poisson fit (`:815-823`), once per slice |
| `t_i` | the slice's `offset.name` column, which is data. Values `<= 0` are floored at `.Machine$double.eps` (`:775-781`). |
| `ps_i` (`ps_hat`) | a **propensity model** fitted once to the full analysis data before the search: `estimate_propensity_scores()` (`:1419-1423`, clamped at `:1426-1434`; GRF `:1482-1494`, LASSO `:1502-1518` or logistic `:1526-1533`), called at `R/forestsearch_main.R:2731-2738`. Alternatively it is the user's `ps_hat` (`R/forestsearch_main.R:2700-2720`). It is carried into every slice and not refitted there. |
| `ips_covar` (a regressor in the fit) | made from `ps_hat` as `ifelse(W == 1, 1/ps_hat, 1/(1-ps_hat))` (`:1459`, or `R/forestsearch_main.R:2726-2728`). This is the Bang-Robins augmentation covariate (comments `:809-810`, `:825-826`, `:1445-1451`). |

- **Current SE:** `|mean(pred1)| * sqrt(vcov(fit)[treat, treat])` (`:852-856`).
  - `d est / d bT = mean(pred1)`, so this SE keeps only the `bT` term of the delta method.
  - It drops the `b0` term (`d est / d b0 = est`) and the `bI` term (`d est / d bI = mean(pred1/ps) - mean(pred0/(1-ps))`).
  - It drops the covariances between coefficients.
  - It ignores the propensity model's variability.
  - It uses the model-based `vcov()` whatever `robust_se` is set to.
- **Why BLOCKED:** The condition is "a smooth function of the coefficients of one fitted outcome model, with nothing estimated elsewhere entering it (no propensity model, no augmentation term)". Both exclusions apply:
  - **Propensity model.** `ps_hat` enters the prediction directly, as `1/ps_i` and `1/(1-ps_i)` (`:830`, `:832`). With `b` held fixed, `d est / d ps_i` is non-zero for every subject:
    - implementer: range [0.0002572, 0.003724], n = 400;
    - verifier: range [0.0001213, 0.03743], 600 of 600.
  - Changing only the propensity model changes `est`, with the outcome data the same:
    - implementer: 0.326687 / 0.327447 / 0.325843;
    - verifier: 0.68293652 vs 0.69648025.
  - **Augmentation term.** `ips_covar` is a regressor in the fit, so `b` itself depends on `ps_hat`. The verifier measured `bI` = -0.128067, z = -2.965.
  - `sqrt(g' V g)` over the outcome fit alone would treat `ps_hat` as known. The task does not authorize that. Values for information only (not proposed):

| Run | `sqrt(g' V g)`, `ps_hat` held fixed | Current SE | Bootstrap SD, whole two-stage procedure (PS refitted per replicate) |
|---|---|---|---|
| implementer, n = 400 | 0.082362 | 0.109610 | not run |
| verifier, n = 600 | 0.0940735 | 0.118239 | 0.122884 (B = 300) |

- **The form was confirmed by measurement on the start tree:**
  - Implementer (`evidence/B2/form_check.R`): the closure's estimate 0.3266870139 equals the written-out formula. Its SE 0.1096101930 equals `|mean(pred1)| * se(bT)`. The analytic gradient (0.326687, 0.854623, 0.410181) equals central finite differences. Wall 4.42 s in the script, 4.62 s for the process.
  - Verifier (own data and script): estimate 0.682936519003 by hand. SE 0.118238831259 is `identical()` to `|mean(pred1)| * sqrt(vcov(glm)[trt, trt])`. Walls 5.84 s / 6.01 s; addendum 2.23 s / 2.38 s.
- **Other branches with a scaled-coefficient SE (listed, not changed).** Search: `git grep -n -E 'se_treat|avg_deriv' b06de4e8 -- R/` returns only these sites in `R/glm_effect_estimators.R`:

| Branch | Start | Head | Estimate | SE | Status |
|---|---|---|---|---|---|
| count, `dr_gcomp` IRR | `:837-848` | `:926-937` | `log(mean(pred1) / mean(pred0))` (`:838`) | `se(bT)`, model-based (`:840`, `:843`). Factor 1. The comment `:839` says "Delta-method SE approximation". | Not changed. The estimate is not `bT`: 0.481686 vs 0.461615 (implementer), 0.53956396 vs 0.58291894 (verifier). Same `ps_hat` dependence as B2. |
| binary `.estimate_or()`, `dr_gcomp` | `:602-628` | `:691-717` | marginal log-OR of the G-computation means (`:623-624`) | `se_treat` (`:625-626`). Factor 1. | Not changed. Same `ps_hat` dependence. Found by reading the code, not measured. |
| binary `.estimate_rd()`, logistic tier, adjusted | `:546` | `:576` | `mean(pred1) - mean(pred0)` | `abs(avg_deriv) * se_treat` | Changed by B1 (`8e87cf30`). At head it is only the fallback when the full delta-method value is not finite and non-negative. |

  - The other count-estimator branches do not scale a coefficient SE:
    - `iptw` (`:797-798`) and the unadjusted IRR (`:885-886`) report the coefficient with its own SE.
    - The unadjusted IRD (`:897-904`) uses the delta method over `(b0, bT)`.
- **What a user sees:** no change.
- **Verifier's minor problem:** NOTES says only `outcome_type = "count"` reaches the branch. `outcome_type = "binary"` with `effect_measure = "IRD"` reaches it too (`:162-172`). The verifier measured `poisson_dr_gcomp`, with estimate and SE identical to the count route. This section states it.
- **Left open by B2:**
  - The `dr_gcomp` IRD and IRR SEs, and the `.estimate_or()` `dr_gcomp` SE, stay as they are. An SE that matches these estimators has to include the propensity model's variability, either through joint estimating equations or by resampling the whole two-stage procedure. GRF and LASSO scores have no closed-form `V`. This is a decision for Larry.
  - No test covers any `dr_gcomp` branch: `git grep -n dr_gcomp -- tests/` finds nothing at `b06de4e8` or at `fb80dd8f`.
  - Found while listing the branches, not changed:
    - `iptw` ignores `effect_measure`: an IRD request returns the log-IRR.
    - Both `dr_gcomp` branches ignore `robust_se`.
    - `estimate_propensity_scores()` documents `trim` as "Quantile bounds", but the code clamps at absolute values.
    - The covariate-adjusted unadjusted-branch IRD is evaluated at covariates = 0.
  - Full evidence: detail file, section B2.

## 5. Documentation

### C4 — help pages and source comments

- **Status:** PASS — commit `ba7c4ea2` (parent `c3323a59`) — 18 `R/` files (roxygen `#'`, plus four `#` comment blocks: the task's three and the `use_dina` screening comment) and 17 `man/` pages regenerated by `devtools::document()`; 35 files, +699/-96. No file under `tests/`.
- **Made of three verified patches:** C4a (display orientation, `adverse_outcome`, `CV_sgs()`, bootstrap bins, HTE mapping; 19 files, +405/-18), C4b (`?forestsearch` and three `#` comments; 6 files, +111/-30), C4c (calibration, CDE, bootstrap noise, GBSG DGM; 10 files, +183/-48). Applying the three final patches to `c3323a59` gives `ba7c4ea2` exactly (`R/`, `man/`).
- **Verification:** each patch CONFIRMED in round 1 with no blocking problem (4, 7 and 7 minor). Each was then polished in comments only and the polish CONFIRMED with no problem. No repair round.
- **Order:** the `?sg_tables`, `?forestsearch_KfoldOut` and `?cv_summary_tables` sentences describe A1/A2's behaviour. A1 (`f4c7dd2d`) and A2 (`a2cd36db`) are committed before C4, and the rule is in the head code.
- Line numbers: start `b06de4e8` → head `fb80dd8f` (`R/` source). Full old and new text, the code behind each sentence, probes and every verifier problem: detail file, section C4.

| # | Site (start → head) | Old text (short) | New text (short) | Differs from the task's wording? |
|---|---|---|---|---|
| 1 | Display pages: `?sg_tables` (`R/summary_utility_functions.R`, after `:735` → `:736-755`), `?forestsearch_KfoldOut` (`R/forestsearch_cross_validation.R`, after `:1325` → `:1326-1339`), `?cv_summary_tables` (`R/cv_summary_tables.R`, after `:34` → `:35-50`), `?plot_sg_glm_outcomes` (after `:25` → `:26-40`), `?plot_subgroup_results_forestplot` (after `:132` → `:133-156`), `?mr_estimates_table` (after `:66` → `:67-84`) | Nothing on orientation | The effect follows the fit's stored `adverse_outcome` (else `forestsearch()`'s rule); arm columns and `Diff` stay raw. For an unadjusted fit with `adverse_outcome = FALSE` they have opposite signs (MD = `-Diff`, RD = `-Diff / 100`, log OR / log RR opposite). Count, IRR, IRD are not reoriented. Five pages are scoped to `est.scale = "hr"`. | Yes. Adds "unadjusted", the sign per measure, the count/rate exception and the `est.scale = "hr"` scope (also on the three pages with no `est.scale` argument, which read the fit's). The plot has no `Diff` and the forest plot shows arm sizes, so both are worded for that. |
| 2 | `make_effect_estimator()` `@param adverse_outcome` (`R/glm_effect_estimators.R:54-60` → `:54-78`) | "Only applies to binary outcomes; ignored for survival, continuous, and count." | Per type: binary OR/RR/RD fit on `1 - Y`; continuous on `-Y`; count and binary IRR/IRD ignore it; survival ignores it. `forestsearch()` resolves its own default. | No. Adds that binary IRR/IRD ignore it and that the log RR on `1 - Y` is the complementary-event ratio, not `-log(p1/p0)`. |
| 3 | `forestsearch()` "Resolved defaults" (`R/forestsearch_main.R:1313-1315` → `:1335-1352`) | "suppresses the corresponding remap and the supplied value is used" | For RD/IRD a supplied value greater than 1 is replaced by 0.05 / 0.01 with two warnings; `consistency.threshold` by 0 with one. The fit records the remapped value. It is the DINA floor on the two routes where DINA selects a cut, gaussian only there. `selected_only = FALSE` has neither. | Yes. "Greater than 1", not "out-of-range" (a value below -1 is kept, no warning). Two warnings, not one. Adds `consistency.threshold`. The floor applies on two routes only. |
| 4 | `forestsearch()` `@param dina_args` (`R/forestsearch_main.R:440-445` → `:447-459`), plus the screening `#` comment (`:3273-3275` → `:3323-3329`) | "i.e. the same cut `subgroup_method = "dina"` would select" | Ranked by `sg_focus` on DINA's tau-hat. It is that cut only under `dina_select_statistic = "dina"`. Under the default `"effect"` the DINA route re-ranks and can select a different cut; screening does not re-rank. | Yes, inverted. The task says the claim holds for the default statistic; the code's default is `"effect"` (start `:1688`, head `:1725`), so it holds only for the non-default `"dina"`. The `#` comment is an extra site. |
| 5 | `forestsearch()` `@param consistency_method` (`:1076-1078` → `:1090-1100`) and `.glm_resample_supported()` (`R/subgroup_consistency_helpers.R:336-338` → `:336-342`) | "IRD, propensity-adjusted effects, ... fall back to split" | The route is decided by the effect measure alone. IRD, or a candidate without a closed form, goes to `"split"`. Under propensity adjustment the closed form is fitted without weights and covariates, so the rate is that of the unadjusted effect. | No. Adds what follows under propensity adjustment. |
| 6 | `CV_sgs()` `@return` (`R/forestsearch_cross_validation.R:1813-1814` → `:1839-1865`) | "List with indicators for any match, exact match, one match, and covariate-specific matches." | Defines each indicator per fold. One-factor "Exact": the fold's subgroup is that factor alone. NA counts as 0. With no factor or more than two, four indicators are 0. | No. Adds the other indicators and the empty / more-than-two case. |
| 7 | `summarize_bootstrap_subgroups()` bins (`R/summarize_bootstrap_subgroups.R:22-23` → `:22-35`) | "Data.table of Pcons distribution by bins" | Bins closed on the left (`"a-b"` holds `a <= Pcons < b`; `">=0.95"` holds `0.95 <= Pcons <= 1`); size bins the same; `Percent = 100 * Count / n_found` | No. |
| 8 | `test_hte_from_forestsearch()` `@return` (`R/test_hte_crump.R:386` → `:386-417`) | "List with hte_test, covariates_used, source." | Lists every element. A label that is a column name is that column; other labels are parsed. Gives the two warnings. Says what `source` holds (the string is fixed). | No. Adds `n_covariates`, `n_candidates`, the third message "No covariates found." and factor coding. |
| 9 | `calibrate_k_treat()` (`R/calibrate_k_treat.R:9-15` → `:9-18`; `@seealso` `:90` → `:93-94`) | "natural sibling of `calibrate_k_inter` (harm-subgroup HR under `model = "alt"`)" | The sibling is `find_k_inter_for_target_hr()`. `calibrate_k_inter()` calibrates the GBSG mechanism of `setup_gbsg_dgm()`. Link added to `@seealso`. | Yes, in form. The page said "natural sibling", not "the general-generator counterpart". The substance is as the task says. |
| 10 | `find_k_inter_for_target_hr()` (`R/find_k_inter_main.R:10-11`, `:35`, `:42-43`, `:51-52` → `:10-13`, `:37`, `:44-45`, `:53-56`) | Class `"k_inter_result"`; a `method` item; "standardized"; "evaluates the boundaries and provides diagnostic information" | A plain list with no `method`. `convergence` counts `uniroot` iterations, including extensions. Covariates stay on their original scale (no argument to change it). A failed calibration stops with an error. | No on the two named points. Beyond the task: the false class is removed and the `@details` failure-path sentence is replaced. |
| 11 | `compute_dgm_cde()` `@param dgm` (`R/oc_analyses.R:109-114` → `:109-118`) | "Must contain `df_super_rand`" | Reads `df_super_rand` when present, else `df_super`. Needs `theta_0`, `theta_1`; otherwise it warns and returns `dgm`. | No. Adds `@param harm_col` and the warn-and-return case. |
| 12 | `generate_bootstrap_with_noise()` (`R/synthetic_data_perturbation.R:429`, `:449-450`, `:461`, `:463`, `:517-518` → `:429-431`, `:451-456`, `:467`, `:469-471`, `:525-537`) | "Preserves integer type"; "same structure as the input"; "NA ... not imputed" | `continuous_vars` come back double. `n` rows with row names `1:n`; `id_var` becomes `1:n`. NA rules: kept in a column used only as continuous. In a categorical column a selected NA becomes the first level (two levels) or a random level (more than two). It stops on an ordered or `grade` column with more than two levels. It stays NA with fewer than two levels. `continuous_vars` must be numeric (logical becomes 0/1). | No on the four points. Adds the logical case, the `grade` column, the fewer-than-two-levels case and the column in both sets. `×` becomes "times". |
| 13 | `setup_gbsg_dgm()` Details (`R/setup_gbsg_dgm.R:74` → `:74-111`) | Noise paragraph only | Three paragraphs (Subgroup, Model, Interaction) from `d93e6621:R/sim_aft_gbsg.R:134-169`, checked against the builder | Yes. The old `HR(H)`/`HR(Hc)` formulas and "`k_inter = 0` gives HR(H) = HR(Hc)" hold for each patient's HR and `AHR_*_true`, not for the Cox `hr_*_true` (1.387 vs 1.504; 0.661 vs 0.585). The text says which is which. |
| 14 | `#` comments: `R/dina_subgroup.R:404` → `:404-406`; `R/forestsearch_main.R:2495-2500` → `:2540-2550`; `R/get_fsdata.R:327` → `:327-329` | "log-IRR for poisson"; "event count threshold ... would filter out all candidates ... (0--1)"; "GLM for binary/continuous" | Log mean-count ratio (`dina()` has no offset). `dmin.grf` is a DR-scale harm threshold (RMST for survival, outcome scale for GLM; default 0.0, any non-negative effect). LASSO family per outcome type, including count. | Only the line pointers ("near 402" → block :403-404, changed line :404; "near 2457-2458" → block :2494-2500, changed lines :2495-2500). The false old `dmin.grf` claims are dropped. |

- **Gates:**
  - Comment-only gate: `main/C4/comment_gate.txt` shows 18 files comment-only, PASS. The drafting agent's re-run of the task's command, `c3323a59` vs `ba7c4ea2`, is empty for 18 of 18. No added line is blank or non-ASCII (the one non-ASCII character, `×` at start `R/synthetic_data_perturbation.R:461`, is removed), and no changed line is an `@examples`/`\dontrun` line. `parse(keep.source = FALSE)` is identical for 18 of 18 files (0.12 s).
  - `devtools::document()`: the main session measured 11.77 s and 11.64 s for the re-run. Those files record wall time only. The drafting agent's runs on `git archive` exports wrote no file at `ba7c4ea2` (13.18 s) or at head (14.46 s). The verifiers regenerated `man/` from `R/` alone, byte-identical to each patch.
  - `tools::checkRd()`: `main/C4/checkRd.txt` shows 17 of 17 OK. The drafting agent's run with the declared encoding (`encoding = "UTF-8", def_enc = TRUE`) gives 0 problems on all 17 at start, at `ba7c4ea2` and at head. The plain call gives only the pre-existing "León" notes (14, 2, 1, 1 on four pages), the same as at start. `generate_bootstrap_with_noise.Rd` goes from 1 note to 0.
- **NEWS** (`fb80dd8f:NEWS.md:536-562`) had three inexact C4 claims; all three were reworded in `2a8be3bc` (§8). The claims as committed in `fb80dd8f`:
  - "an out-of-range ... threshold is remapped with a warning": only values greater than 1 are remapped, and `effect.threshold` gets two warnings.
  - "the remapped value is also the DINA floor": only `effect.threshold`, and only where DINA selects a cut.
  - "the `use_dina` screening cut matches `subgroup_method = "dina"` only under `dina_select_statistic = "dina"`": under the default the cut can still match (5 of 10 verifier runs). The help says "can select a different cut".
  - The display bullet also drops the pages' "unadjusted" and `est.scale = "hr"` qualifiers. It says "can", so it is not false.
- **Open items:**
  - `test_hte_from_forestsearch()` `source` string "confounders.evaluated (GRF/LASSO selected)" (head `R/test_hte_crump.R:447`) is wrong for a default fit and the DINA/GRF routes. It is a code string, and `test-hte-crump-confounders.R:76` asserts it. The same fault is in `?forestsearch` `@return` `confounders.evaluated` (head `R/forestsearch_main.R:1225`).
  - `est.scale = "1/hr"` GLM fits: the subgroup rows' arm columns and `Diff` are computed on swapped treatment while the effect is not (`prepare_subgroup_data()`, `R/summary_utility_functions.R:51-58`). The help is scoped to `"hr"`; the code is unchanged. Separately, `mr_estimates_table(est.scale = "1/hr")` inverts the H cell but not the complement's naive cell (observation).
  - NEWS: the inexact claims above were corrected in `2a8be3bc`.
  - Pre-existing faults found while documenting, unchanged (code, or outside the site list):
    - `args_call_all` keeps the supplied `hr.threshold` beside the remapped `effect.threshold`.
    - `R/forestsearch_helpers.R:1607-1610` calls `"dina"` the default statistic.
    - The `consistency_method = "resample"` + propensity-adjustment design question (`R/consistency_resample.R:216-221`).
    - `CV_sgs()` gives 0 for a full-data subgroup of three or more factors, and `cv_metrics_tables()` prints "Exact match on all factors" (head `R/cv_summary_tables.R:434`).
    - `?plot_sg_glm_outcomes` `@param outcome_type` omits count.
    - `grf.subg.harm.glm()` defaults `adverse_outcome = FALSE` (`R/grf_subg_harm_glm.R:304`).
    - The `dmin.grf` `missing()` reset has no effect (head `:2551-2553`).
    - The `compute_dgm_cde()` warning names `df_super_rand` when `df_super` was read (head `R/oc_analyses.R:144`).
    - "typically converges within 5-10 iterations" (head `R/find_k_inter_main.R:52`; 11 and 19 observed).
    - `@param k_treat_range` "(uniroot will error otherwise)" (head `R/calibrate_k_treat.R:29-30`).
    - `generate_bootstrap_with_noise()` on a one-column data frame returns a vector with no noise.
  - Count/IRD fits with `fs.splits = 10`, `use_twostage = TRUE` fail inside `subgroup.consistency` ("argument must be coercible to non-negative integer"), as seen on the start tree by the C4b verifier. Not investigated; not caused by C4.
- Full evidence: detail file, section C4.

### C5 — examples that cannot run as written

- **Status:** PASS — commit `5073fd67` (parent `ba7c4ea2`), documentation only.
  - 11 `R/` files, `@examples` roxygen lines only, +99/-26.
  - 16 `man/` pages regenerated, +99/-26.
  - 27 files in all, +198/-52.
  - NEWS line in `fb80dd8f`, reworded in `2a8be3bc`.
- **Result:** 16 help pages.
  - At the start SHA, 15 examples stopped with an error and 1 (`sg_tables`) ended without an error but printed no table.
  - Fixed: 15 RAN and 1 NOT RUN (`fs_fdr_report`, over the 2-minute cap).
- **How each example was run.**
  - The example was taken from its Rd with `tools::Rd2ex(commentDontrun = FALSE)`, with `library(forestsearch)` first.
  - It ran once in a fresh `Rscript --vanilla` under `timeout 120`.
  - The library was a temporary install of the C5 copy. Its R code parses identically to the start's, so the start Rd was run the same way.
  - Implementer, verifier and head runs were made the same way. The head runs were made by a read-only after-the-fixes agent from a library with the head's `R/` and `man/`.
- **Line numbers:** start `b06de4e8` → head `fb80dd8f`. `R/` and `man/` are the same at `5073fd67` and `fb80dd8f`.
- **The shared fix, "L+G+P" below**, on twelve pages:
  - `gbsg <- survival::gbsg`. A bare `gbsg` is `object 'gbsg' not found` after `library(forestsearch)`.
  - `gbsg$grade3 <- as.integer(gbsg$grade == 3)`. `grade3` is not a GBSG column.
  - `id.name = "pid"` in the call. The default `id.name = "id"` names no GBSG column, and without it `forestsearch()` stops with `undefined columns selected`.

| Rd page | file:line, start → head | What was wrong (start run) | What changed | Result | Wall, s (impl / verifier / head) |
|---|---|---|---|---|---|
| `format_bootstrap_table` | `R/bootstrap_summaries_helpers.R:28`, `:32-33` → `:28-30`, `:34-35` | `object 'gbsg' not found`; also `summaries$table_data` does not exist and the required `nb_boots` was not passed | L+G+P; `summarize_bootstrap_results()` line dropped; `format_bootstrap_table(fs_bc$FSsg_tab, nb_boots = 100, sg_definition = paste(fs$sg.harm, collapse = " & "))` | RAN (gt table) | 17.52 / 17.23 / 17.75 |
| `print.fs_kfold` | `R/forestsearch_cross_validation.R:2040` → `:2091-2093` | `object 'gbsg' not found` | L+G+P | RAN (686 folds, found in 99.4 %) | 67.07 / 67.25 / 68.54 |
| `print.fs_tenfold` | `R/forestsearch_cross_validation.R:2069` → `:2122-2124` | `object 'gbsg' not found` | L+G+P | RAN (5/5 simulations) | 14.40 / 14.64 / 14.61 |
| `plot_sg_distribution` | `R/plot_sg_distribution.R:75-79`, `:81` → `:75-83`, `:85-92` | `object 'gbsg' not found`; `fs$grp.consistency$result_new` does not exist; `results_alt` never defined | `forestsearch()` call dropped; two stand-in tables in the two documented schemas, `sims_mrct` and `results_alt`; plot calls unchanged | RAN (3 bars; 15 bars at 2.5 %, `fig_h` 8.25, as derived by hand; one pdf-device glyph warning) | 1.78 / 1.72 / 1.71 |
| `print.fs_sg_plot` | `R/plot_sg_results.R:866`, `:869` → `:866-868`, `:871-872` | `object 'gbsg' not found`; `plot_subgroup(fs)` calls the internal KM helper, not the `fs_sg_plot` constructor | L+G+P; `p <- plot_sg_results(fs, outcome.name = "rfstime", event.name = "status", treat.name = "hormon")` | RAN (empty "Available plots:" line, open item 4) | 6.49 / 6.36 / 6.44 |
| `plot.fs_sg_plot` | `R/plot_sg_results.R:912`, `:915` → `:915-917`, `:920-921` | as `print.fs_sg_plot` | as `print.fs_sg_plot` | RAN ("No plots available", open item 4) | 7.06 / 6.16 / 6.42 |
| `print.fs_forestplot` | `R/render_forestplot.R:372` → `:372-374` | `object 'gbsg' not found` | L+G+P | RAN | 6.39 / 6.19 / 6.44 |
| `plot.fs_forestplot` | `R/render_forestplot.R:410` → `:412-414` | `object 'gbsg' not found` | L+G+P | RAN | 6.44 / 6.32 / 6.76 |
| `run_subgroup_sims` | `R/run_subgroup_sims.R:652` (insert after), `:663` → `:653-674`, `:685` | `object 'dgm_uniform' not found`; `subgroups`, `cal_uniform`, `k_treat_uniform`, `gbsg` also undefined; `saveRDS()` to `results/` | 22 lines build them as in `extreme_subgroups.qmd`, shortened to 3 subgroups; `saveRDS()` to `tempdir()`; the `run_subgroup_sims()` call (`n_sims = 1000`) unchanged | RAN (1000 simulations; prints nothing) | 13.01 / 12.93 / 13.35 |
| `sg_consistency_out` | `R/subgroup_consistency_helpers.R:1031` → `:1035-1037` | `object 'gbsg' not found` | L+G+P | RAN | 6.51 / 6.28 / 6.50 |
| `evaluate_subgroup_consistency` | `R/subgroup_consistency_helpers.R:1361` → `:1367-1369` | `object 'gbsg' not found` | L+G+P | RAN | 6.38 / 6.26 / 6.66 |
| `evaluate_consistency_twostage` | `R/subgroup_consistency_helpers.R:1645` → `:1653-1655` | `object 'gbsg' not found` | L+G+P; `use_twostage = TRUE` kept | RAN | 6.34 / 6.32 / 6.39 |
| `fs_fdr_report` | `R/fs_fdr_report.R:151` (insert after), `:165`, `:170` → `:152-158`, `:172`, `:177` | `object 'gbsg_df' not found`; `n_workers = 100` in both calls is an error under parallelly's limit (3 x cores) below 34 cores | `gbsg_df <- within(survival::gbsg, { id <- pid; time <- rfstime / 30.4375; treat <- hormon })`; `n_workers = 100` dropped from both calls (default is about 90 % of physical cores, 115 here) | **NOT RUN: over the cap** (2 calls x 2 nulls x 1000 MR fits). Killed while "Running 1000 reps/null x 2 null(s) at n = 686 on 115 workers". At `nsims = 4, n_workers = 4` both calls ran to the end, rc 0, verdict CLEAN | 120.76 killed (rc 124) / 120.76 killed / not run |
| `sg_tables` | `R/summary_utility_functions.R:742`, insert after `:744` → `:762-763`, `:766` | rc 0, no table: `library(survival)` finds `gbsg`, but `forestsearch()` warns about `grade3` and returns a NULL result; `sg_tables()` warns; prints `NULL` | `library(survival)` → L+G; `id.name = "pid",` added | RAN (gt table) | 6.89 / 6.76 / 6.89 |
| `pareto_frontier_table` | `R/pareto_frontier_table.R:125` → `:125-128` | placeholder `forestsearch(gbsg, ...)`: `'...' used in an incorrect context` | G+P and the full call (the existing `gbsg <- survival::gbsg`, `:124`, kept) | RAN (one-member frontier, HR 1.72) | 10.06 / 9.88 / 10.66 |
| `summarize_simulation_results` | `R/oc_analyses.R:414` → `:418-424` | `unused argument (nsim = 10)`: `run_simulation_analysis()` has no `nsim` and runs one replicate per call | ten stacked calls: `data.table::rbindlist(lapply(1:10, function(i) run_simulation_analysis(sim_id = i, dgm = dgm, n_sample = 500, confounders_base = c("z1", ..., "size", "grade3"))), fill = TRUE)` | RAN (FS / FSlg / GRF table) | 18.08 / 17.77 / 18.81 |

- **Start-run walls:** 1.48-1.64 s (implementer) and 1.45-1.60 s (verifier).
- **Head re-run** (after-the-fixes agent, `$S/examples_head/`): the 15 RAN pages all exit 0.
  - After gt's random table ids and the K-fold time line are normalised, each output equals the implementer's and the verifier's fixed output, 15/15.
  - `fs_fdr_report` was not re-run at head. Its extracted head code equals the code that was run.
- **Gates:**
  - **Comment-only gate** (the task's command, each of the 11 `R/` files): empty for every file. Main session: `main/C5/comment_gate.txt` says `COMMENT GATE: PASS`. Re-run for this report on `git diff -U0 ba7c4ea2 5073fd67`: also empty. No added line is non-ASCII.
  - **Only `@examples` lines changed.** The verifier found 0 changed lines outside `@examples`.
  - **Parsed code identical** start vs fixed for 11/11 files (`parse(keep.source = FALSE)`, implementer and verifier).
  - **`tools::checkRd()`:** 16/16 `OK (0 problems)` (main session, `main/C5/checkRd.txt`). The implementer and verifier got 16/16 OK with `encoding = "UTF-8"`. A plain `checkRd()` gives a non-ASCII note ("León") on two pages, outside the examples; the start pages give the same note.
  - **`man/` equals the verified patch.**
    - The commit's 250 changed lines in 27 files equal those of `patches/C5.patch`: `diff` is empty, checked for this report. The patch is not repaired or polished, and it predates `VERIFY.md`.
    - The verifier ran `document()` on start + patch. It rewrote exactly the 16 pages, byte-identical to the patch's `man/`, and the Rd text outside `\examples{}` is unchanged.
    - At head, 15 of the 16 pages are byte-identical to the C5 copy. `sg_tables.Rd` differs only by C4's 20-line `\details{}` block.
    - `\dontrun{}` is kept on all 16 pages.
    - The main session's `main/C5/document.txt` is empty, so no `document()` wall was recorded for this commit.
  - **Scope:** `git diff --stat ba7c4ea2 5073fd67` lists only the 11 `R/` files and 16 `man/` pages. No tests changed, so no tests were run.
- **Verification:** CONFIRMED in round 1, with no repair and no polish. Three minor points:
  1. NOTES calls dropping `n_workers = 100` "workers lowered", but the default here is 115. The hard-limit error below 34 cores is the reason that holds.
  2. The proposed NEWS bullet overstated the result (see open item 6).
  3. `plot.fs_sg_plot`'s `@return` does not match what the corrected example shows (open item 4).
- **Departures from the task's wording:**
  - `id.name = "pid"` is added to every GBSG call, because the task's two fixes alone still stop.
  - `plot_sg_distribution` uses stand-in tables, since a fit has no per-replicate table.
  - Lines that call a non-existent element (`summaries$table_data`) or the wrong function (`plot_subgroup(fs)`) are replaced.
  - `run_subgroup_sims` is shortened to 3 subgroups.
  - `fs_fdr_report` drops `n_workers = 100`; the task says nothing about workers.
  - `summarize_simulation_results` passes `confounders_base`. The default `v1`-`v7` are factors in the `setup_gbsg_dgm()` data, and every FS replicate fails inside.
- **Pointers:**
  - The task's "near 2002, 2031" are `2040` and `2069` at start.
  - The task's "near 148-164" is the block `:150-172`.
  - Part 3 said the `sg_tables` example makes `forestsearch()` stop. It warns and returns a NULL result.
- **Open items:**
  1. `fs_fdr_report` was not run in full (over the cap) and not at head. Its comment "published reference (n = 700): A ~ 0.138 (hom), 0.278 (harm)" (head `R/fs_fdr_report.R:176`) is not verified.
  2. Three GBSG examples outside the C5 list stop with `undefined columns selected` because their call has no `id.name` (measured):
     - `summarize_bootstrap_results`: `R/summarize_bootstrap_results.R:52-57`.
     - `forestsearch_KfoldOut`: `R/forestsearch_cross_validation.R:1331-1335` → `:1345-1349`.
     - `summarize_cv_results`: `R/summarize_cv_results.R:148-159`.
  3. `run_simulation_analysis()`'s default `confounders_base = c("v1", ..., "v7")` (`R/run_simulation_analysis.R:398`) fails inside every FS replicate on a `setup_gbsg_dgm()` DGM. Its own example uses that default and was not run.
  4. `plot_sg_results()` always returns `plots = list()`, because its three panel helpers return `invisible(NULL)` (`R/plot_sg_results.R:661`, `:775`, `:849`). So `plot.fs_sg_plot` prints "No plots available", unlike its `@return` (head `:909-911`).
  5. Examples with placeholder data (`my_data`, `mydata`, `trial_data`) were not in this item.
  6. NEWS (`fb80dd8f:NEWS.md:556-559`) said "Sixteen `\dontrun{}` examples now run as written", but `fs_fdr_report` ran only at reduced size, and the list of causes omitted the wrong call, the missing list elements and `n_workers = 100`. Corrected in `2a8be3bc`: sixteen corrected, fifteen run to the end; the causes listed "among other defects".
- Full evidence: detail file, section C5.

### C6 — comments on test expectations

- **Status:** PASS. Commit `ffc55486` (parent `5073fd67`), `test: derivation comments on part 4's test expectations (part 5, C6)`.
  - It changes 9 files under `tests/testthat/` (+111/-8), comment lines only. No file under `R/`, `man/` or `vignettes/` changes.
  - The commit equals `patches/C6.patch` (the polished patch); `main/C6/patch_equal.txt` is empty.
- **Requests:**
  - Part 4 report §8, OPEN ITEMS 17, has 13 requests in the test files of F1, F5, F6, F8, F9, F10, F12, F13 and F14.
  - Part 4's detail file, F6 problem 2, adds 3 optional ones (`:65`, `:71`, the survival case of test 4).
  - All 16 are done.
- **Line numbers** are at start `b06de4e8` → head `fb80dd8f`.
  - Each of the nine files is the same blob at its part-4 fix commit and at `b06de4e8` (9 of 9). So item 17's lines are start lines.
  - No other part-5 commit touches these files, so head lines are `ffc55486` lines.
- **Verification:**
  - CONFIRMED, with no repair round. The verifier found three minor wording points (F10 `:111-112`, F6 test 4 header, F10 frontier upper bound).
  - All three were fixed in a polish round: comment lines only, +15/-10 against the verified patch.
  - The polish check gave CONFIRMED with no problems.

| # | File (`tests/testthat/`) start → head | Request | Comment added (short) |
|---|---|---|---|
| 1 | `test-bootstrap-summary-bins.R:76-83` → comment `:76-81`, checks `:82-89` | Label-order checks pass on both trees | They pass before and after; the labels come from the `cut()` calls and `table()` keeps level order; they anchor the row order the Count/Percent checks rely on |
| 2 | `test-hte-crump-confounders.R:97-99` → `:97-105`, checks `:106-108` | Derive 0 warnings, non-NULL, class `hte_test` | By hand from the labels (`all.vars()` gives `age`/`bmi`, both columns); pre-fix, only the warning count passed |
| 3 | `test-lasso-intercept.R:88` (test 2) → `:94-99`, test 2 at `:100` | Test 2 passes before and after | It mocks the function being fixed, so it cannot go red; it is a lock on `get_FSdata()`'s rule; the red is test 3 |
| 4 | `test-lasso-intercept.R:65` → `:65-66`, check `:67` | Optional: passes on both trees | The fix leaves `coefficients` (every row of `coef(fit)`) as it was |
| 5 | `test-lasso-intercept.R:69-71` → `:70-75`, checks `:76-77` | Optional: `:71` passes on both trees | `selected` is the red; `omitted` passes before and after (the intercept is non-zero) |
| 6 | `test-lasso-intercept.R:181-185` (test 4 header) → `:198-205`, checks `:228-232` | Optional: survival passes on both trees | Survival passes before and after; pre-fix, binary printed no line (both checks fail), continuous/count printed "Cox-LASSO ..." (match fails). Rewritten in polish |
| 7 | `test-duplicate-cut-names.R:162-178` → `:163-172`, checks `:178`, `:187` | The guard does not pin the `%in% confs` clause | `unique()` drops the exact duplicate anyway; a copy without the clause passes every check; the clause only affects `protect` (`collapse_cuts = TRUE`) |
| 8 | `test-mrct-default-confounders.R:40-49` → `:50-56`, `expected` `:57` | Why the other `df_super` columns are absent | Only `^z_\|^ecog\|^strat` names enter; `zcens_*` do not match `^z_`; only the three helpers are removed by name |
| 9 | `test-mrct-default-confounders.R:63` → `:70-72`, check `:73` | One line: one replicate, `verbose_n_sims = 1` | `n_sims = 1`; printed for `sim <= verbose_n_sims` inside `if (!is.na(hr_test))` |
| 10 | `test-mrct-default-confounders.R:82` → `:92-94`, check `:95` | `:82` is a regression guard | Passes before and after; the identifier `!{z_male}` uses no helper column |
| 11 | `test-grf-count-orientation.R:103-106` → `:103-112`, checks `:113-116` | Derive the tree existence checks | One region has positive harm score; ~500 members, >= `n.min` 60, < n = 1000, score > `dmin.grf` 0; pre-fix both fits found the same z == 0 subgroup (polish) |
| 12 | `test-grf-count-orientation.R:129-132` → `:139-145`, checks `:146-149` | Derive the frontier existence checks | Eligible: 60 to 999 members, effect >= 0; the z candidate has 500 members, effect ~6 / ~2; pre-fix both selected z <= 0. Upper bound added in polish |
| 13 | `test-display-orientation.R:27-30` → `:27-38` (rewritten); 1e-4 check `:151` → `:166` | 1e-4 justified by a z difference the test does not have | Both the printed and the hand MR limits use 1.96; only 4-decimal rounding (<= 5e-5) separates them; 1e-4 is twice that |
| 14 | `test-display-orientation.R:127-130` → `:135-141`, checks `:142-145` | Explain the forest-plot 1e-10 | Unrounded estimator values vs the same OLS by hand; 1e-10 is a relative numerical tolerance, far below the 0.19-2.2 a sign flip moves |
| 15 | `test-rd-standard-error.R:78-81` → `:78-83`, checks `:84-87` | No comment for `converged`, `method_used`, `measure`, `n0`/`n1` | Tier 2 values with no covariates; arm sizes 40 by construction; pass before and after (the fix changes only the SE) |
| 16 | `test-dina-gaussian-floor.R:159-161` → `:159-168` (rewritten) | Cite the remap "to default" step | The branch `user_set_threshold && effect_threshold > 1.0` sets `switch(effect_measure, RD = 0.05, IRD = 0.01)`; these equal the no-threshold defaults; `.gf_capture()` muffles both warnings |

- **Existing comments rewritten** (the 8 deleted lines): F6 `:69` (1 line, extended), F12 `:27-30` (4 lines; the old text claimed a z difference that is not there), F14 `:159-161` (3 lines; the old text cited only the no-threshold table).
- **Gates:**
  - **Comment-only:** `main/C6/comment_gate.txt` gives "COMMENT GATE: PASS" for all 9 files. Re-run on `git diff -U0 5073fd67 ffc55486`, the task's command prints 0 lines for each file.
  - **Code unchanged:** with comment and blank lines stripped, each file is identical to start (9 of 9). Parsed code is identical (verifier, 9 of 9).
  - **Added lines:** 0 blank, 0 non-ASCII, 0 over 80 characters.
  - **No changed values:** no expected value, seed or tolerance changed.
- **Test counts, before vs after** (expectations; blocks from the scratch runs, because `tests_after.txt` does not show blocks):

| File | Part 4 (§6) blocks / exp | Main after (`main/C6/tests_after.txt`) exp |
|---|---|---|
| `test-bootstrap-summary-bins.R` | 3 / 8 | 8 |
| `test-hte-crump-confounders.R` | 5 / 35 | 35 |
| `test-lasso-intercept.R` | 4 / 62 | 62 |
| `test-duplicate-cut-names.R` | 7 / 13 | 13 |
| `test-mrct-default-confounders.R` | 1 / 8 | 8 |
| `test-grf-count-orientation.R` | 2 / 14 | 14 |
| `test-display-orientation.R` | 6 / 93 | 93 |
| `test-rd-standard-error.R` | 2 / 14 | 14 |
| `test-dina-gaussian-floor.R` | 6 / 55 | 55 |
| **Total** | **36 / 302** | **302** |

  - **Main run:** all nine files in one process, real 23.079 s, run before the commit. It shows no failures, skips or warnings. **PASS.**
  - **Implementer:** each file in its own process, on the start copy and the patched copy. Blocks and expectations are identical (36 / 302 both), all pass. Walls 3.04-14.85 s before and 3.13-15.50 s after.
  - **Verifier:** the same per-file check gave identical counts. Walls 3.02-15.18 s on start and 3.03-15.20 s patched.
  - **Polish re-run** of the two files it changed: `test-lasso-intercept.R` 4 / 62 (5.93 s) and `test-grf-count-orientation.R` 2 / 14 (4.47 s).
- **How the claims were checked:**
  - Every "passes before and after" comment was checked on a pre-fix copy: start plus the fix commit's `R/` diff reverse-applied (F1, F5, F6, F9, F10, F13), with each expectation's outcome recorded. F8 was checked with a mutant without the `%in% confs` clause, which passes 13 of 13.
  - The verifier repeated these runs. The table of outcomes is in the detail file.
- **Departures:**
  - The three optional F6 comments were added.
  - Three existing comments were rewritten instead of extended.
  - The F14 comment names the remap branch and the `switch()` rather than quoting the warning text, which C2 changes.
  - At head, after C2 (`2984dbee`), the `switch()` is assigned to `.fs_remap_default` (`R/forestsearch_main.R:2376`) and then to `effect_threshold` (`:2394`), inside the same branch (`:2371`). The comment still holds.
- **NEWS:** at `fb80dd8f` nothing named C6, although that commit's message says its Documentation entry is "for C4-C7". `2a8be3bc` added a C6 line to the entry (nine test files carry the comments their verifiers asked for; comments only). C6 changes nothing a user sees.
- Full evidence: detail file, section C6.

### C7 — two article sources

- **Status:** PASS. Commit `835fac1b` (parent `ffc55486`): `vignettes/articles/biomarker_effects.qmd` +3/-1 and `vignettes/articles/paper_simulations.qmd` +46/-16 (2 files, +49/-17). The change is prose only and nothing was re-rendered. The commit equals `patches/C7.patch`; the two differ only in the length of the `index` hashes.
- **Line numbers:** start `b06de4e8` → head `fb80dd8f`. Both files are identical at `835fac1b` and `fb80dd8f`.
- **Verification:** round 1 CONFIRMED, with 0 blocking and 6 minor problems. A polish pass fixed minors 2 and 3. The polish check then found 2 problems in the new treatment-column sentence, and the main session reworded that sentence. No verifier ran on the reworded sentence (see open item 4).

**Old → new text per site** (the full text is in the detail file, C7.2):

| Site | Start → head | Old | New |
|---|---|---|---|
| B | `biomarker_effects.qmd:66` → `:66`, `:69-70` | "built from the GBSG data by `create_gbsg_dgm()` and `generate_aft_dgm_flex()`" | "… by `setup_gbsg_dgm()` and `generate_aft_dgm_flex()`" + "`setup_gbsg_dgm()` uses the same internal builder as the removed `create_gbsg_dgm()`." |
| P1 | `paper_simulations.qmd:366-367` → `:366-369` | "the binary discretisation of continuous covariates" | "the discretisation of continuous covariates (age, progesterone receptor and nodes in quartile groups; tumour size in three groups cut at 20 and 50 mm)" |
| P2 | `:388`, `:391-393` → `:390`, `:393-396` | Heading "Replicating `setup_gbsg_dgm()` with `generate_aft_dgm_flex()`"; "The GBSG DGM in `setup_gbsg_dgm()` is equivalent to the following call — it simply automates this preparation" | Heading "A comparable construction with `generate_aft_dgm_flex()`"; "a comparable construction of the GBSG mechanism with the general generator. It is not numerically equivalent to `setup_gbsg_dgm()`; the differences are listed after the code." |
| P3 | `:442-445` → `:445-475` | "`setup_gbsg_dgm()` is preferred … because it guarantees exact replication of the León et al. (2024) covariate construction and subgroup definition …" | One paragraph and three bullets. Shared: source data, outcome in months, $H$, `k_treat`, `k_inter`, `n_super`, `seed`, and the same resampled rows. Different: the super-population treatment column (0.5 vs about 0.36, which matters only when treatment is not redrawn); **covariates and cut points** (quartile groups vs median splits, plus size and grade); **interaction** (same `k_inter` multiplication of a different fitted coefficient: HR in $H$ about 4.4 vs 4.1); **censoring** (Weibull on the super-population with the scale-ratio factor vs an AIC choice on the source data, here intercept-only Weibull). |
| P4 | `:454` → `:484` | Table row "GBSG-based simulation with custom parameters → `setup_gbsg_dgm()`" | "New GBSG-based simulation with custom parameters → `generate_aft_dgm_flex()`" |
| P5 | `:459-463` → `:489-493` | "`setup_gbsg_dgm()` is not deprecated … it remains the correct entry point for the simulation studies in this package." | "`setup_gbsg_dgm()` is used in this article to reproduce the published simulations, and the package keeps it for that purpose. It is superseded for new work and will be deprecated in a future release: new simulations, on the GBSG data or another dataset, should use `generate_aft_dgm_flex()` with `simulate_from_dgm()`." |

- **What the new claims rest on:**
  - Code at `b06de4e8`, unchanged at head except that `R/setup_gbsg_dgm.R` lines from 75 on moved by +37 (C4).
  - Builders compared on the start tree, running the article's own block and `setup_gbsg_dgm()` call: implementer 2.94 s, verifier 3.33 s, polish-check probe 3.27 s.
  - Measured: HR($H$) 4.3808 vs 4.0743; super-population treatment share 0.5000 vs 0.3598; censoring fit mu 4.0492 / scale 0.4363 vs intercept-only Weibull 4.0312 / 0.4246; resampled `pid` sequence `identical()` TRUE.
  - "Superseded for new work" and "deprecated in a future release" come from `?setup_gbsg_dgm`, section Superseded (`R/setup_gbsg_dgm.R:80-100`; head `:117-137`).
  - The claim-by-claim table is in the detail file, C7.3.
- **Departures from the task's wording:**
  - B adds one sentence. The cited paper names `create_gbsg_dgm()`: fs-post-selection `manuscript/postselection_subgroup_inference_supplement.qmd:1182`.
  - P1 and P4 are outside the named ranges (388–445, 459–463). Each removes text that contradicted the new prose or the code.
  - P3 adds the super-population treatment column, a difference the code shows.
  - `k_inter` and the seed are stated as **not** differences.
- **Main-session correction:**
  - Polisher's sentence: "`simulate_from_dgm()` uses this column only when it keeps the super-population's arms, with `draw_treatment = FALSE` or with `n = NULL` …; otherwise it redraws treatment in every simulated trial."
  - The polish check found two problems with it. Under `baseline = "fixed"`, the arms come from `df_source`, so "otherwise it redraws" was false there. And `n = NULL` stops on a `setup_gbsg_dgm()` object.
  - Committed sentence (head 453–455): "In a simulated trial this column matters only when treatment is not redrawn; `simulate_from_dgm()` redraws it in every trial by default (`draw_treatment = TRUE`)."

**Gate** (`main/C7/chunk_gate.txt`: `C7 GATE: PASS`)

| Check | Result |
|---|---|
| Only the two files change under `vignettes/` | Commit: the 2 `.qmd` files only. Whole part (`git diff --name-only b06de4e8 fb80dd8f`): under `vignettes/` only these two; nothing under `quarto/` |
| No changed line inside a code chunk | `biomarker_effects.qmd` 4 changed lines (old + new), 0 in a chunk or fence; `paper_simulations.qmd` 62, 0. Every chunk is byte-identical before and after (5/5, 64/64; re-checked from git) |
| No `.html` edited | `html edited: []`; no `.html` in the part's diff |
| Nothing rendered | No render step was run. `vignettes/articles` is in `.Rbuildignore` (line 20), so neither file reaches build or check |
| `R/`, `tests/`, `man/` | Unchanged by C7; no tests apply |

- Earlier gate runs, all PASS: implementer 0.05 s, verifier's own fence check 0.04 s, polish 0.05 s. The polish run came before the main-session correction.
- **Articles:** no code chunk changed, so no evaluated call and no article output can change from C7.

**Open items (C7)**

1. The chunk label `generate-aft-equivalent` (head `paper_simulations.qmd:398`) still says "equivalent". The in-chunk comment "Analyst-visible covariates (what FS will search over)" (`:430`) is misleading, because `factor_vars` enter the outcome model. Both are inside a chunk, so the gate kept them.
2. The committed `vignettes/articles/biomarker_effects.html` still shows `create_gbsg_dgm()` (1 occurrence) until the article is re-rendered. By rule it was not edited.
3. "About 4.4 … 4.1" (head `:467`) and "about 0.36" (`:453`) were measured on the current builders. They go stale if either builder changes.
4. No verifier ran on the reworded sentence (head `:453-455`). Read against the code it holds for trial simulation, but there is a gap: with the default `n = NULL`, `simulate_from_dgm()` keeps the super-population's arms whatever `draw_treatment` is (`R/simulate_from_dgm.R:221-234`). The article's own calls pass `n = 700` (head `:509`, `:757`).
5. Code problem found by the polish check and left unchanged because it is outside C7: `simulate_from_dgm(dgm)` with `n = NULL` stops on a `setup_gbsg_dgm()` object with "replacement has 0 rows, data has 5000" (no `lin_pred_obs` in `df_super`). The unchanged line at head `:499-501` says `simulate_from_dgm()` "works identically with either".
6. Out of scope and unchanged: a stray `---` / `title:` / `---` block mid-document (`paper_simulations.qmd:358-360`).

## 6. Existing expectations that were updated

None. No existing test failed at any step, and no existing expectation was changed:
- every code fix's test is a new file, and each verifier ran the existing test files that exercise the changed functions on its fixed copy (detail file, "Existing tests" of each item);
- C6 changed comment lines only; its nine files pass with the same counts as part 4 (8, 55, 93, 13, 14, 35, 62, 8, 14 expectations; `$S/main/C6/tests_after.txt`, 23.1 s);
- in the full suite every existing file has the same blocks, expectations, skips and warnings as in part 4 (§7.5).

## 7. After the parts

### 7.1 Default-path gate (step 1) — **PASS**

- The six equivalence cases, run with part 4's scripts unchanged (`REPORT_cran_prep_code_fixes_2026-11-12_gate_run.R` and `_gate_compare.R`): the package test helpers' synthetic data, `set.seed(123)` before each call, sequential plan, `pkgload::load_all()` of each tree in its own `Rscript --vanilla` process.
- Trees: the start-SHA export `$S/start` and `$S/head_fixes`, a `git archive 835fac1b` export of the package paths (`R/` is identical at `fb80dd8f` and `2a8be3bc`).
- Output: `REPORT_cran_prep_displays_se_docs_2026-11-12_gate_compare.txt`.

| # | Call | `identical()` after dropping timing fields | Non-timing leaves that differ | Timing leaves dropped | Wall, start / head |
|---|---|---|---|---|---|
| 1 | `forestsearch()`, survival, defaults | TRUE | 0 | `$find.grps$time_search`, `$minutes_all` | 0.6 / 0.5 s |
| 2 | same, `mr_inference = TRUE` | TRUE | 0 | the two above + three `timing_seconds` | 7.8 / 7.3 s |
| 3 | same, with `df.test` | TRUE | 0 | 2 (as row 1) | 0.2 / 0.1 s |
| 4 | `forestsearch_bootstrap_dofuture()` on fit 1, `nb_boots = 10` | TRUE | 0 | `$results$tmins_search`, `$results$tmins_iteration` | 5.6 / 5.2 s |
| 5 | `forestsearch_Kfold()` on fit 1, `Kfolds = 5` | TRUE | 0 | `$timing_minutes` | 1.3 / 1.3 s |
| 6 | `forestsearch()`, binary outcome, defaults | TRUE | 0 | 2 (as row 1) | 0.6 / 0.6 s |

- `.Random.seed` after each call equals the state before it in both trees. Whole script: 18.5 s (start), 17.3 s (head).

### 7.2 Every code fix's unaffected case at the head (all fixes in)

Each fix was gated on a tree holding only its own patch. Read-only agents re-ran every code fix's unaffected-case scripts, the implementer's and the verifier's, with the fixed tree replaced by `$S/head_fixes`. The start tree was unchanged, both sides were re-run, and the comparison was the original script's own. Any case that differed at the head was traced to a named fix by reverse-applying that fix's `R/` hunks to a copy of the head tree and re-running it. Output: `$S/cumulative/<ID>/RESULT.md`.

| Fix | Cases | Identical to start at head | Differences, traced |
|---|---|---|---|
| A1 | 31 | 28 | 3 GRF count cases differ only in the CATE source note, "(risk-difference scale)" → the count scale: **C3**. With C3's `R/` hunk reverse-applied to the head tree, 31/31 are identical to start. |
| A2 | 44 | 44 | — |
| A3 | 38 | 35 | 3 cases are A3's own intended change (continuous fit with the stored value removed: every estimate and limit changes sign); each equals A3's own gated output, and with A3 reverse-applied all 38 are identical to start. |
| B1 | 2679 | 2253 | 426 differ, all B1's own intended change in the covariate-adjusted logistic tier: 385 estimator results and 27 probe results differ in `$se` only; 12 displays differ in CI limits (e.g. `0.8044 (-830.7797, 832.3885)` → `0.8044 (0.6838, 0.9251)`); 2 `forestsearch()` RD + `adjust_covariates = "site"` runs differ only in the candidate table's `L(HR)`/`U(HR)` (same rows, HR and selected subgroup). The unaffected sets themselves (unadjusted RD incl. F13, OR, RR, IRR/IRD, identity tier) are identical. With B1 reverse-applied, identical to start in full. |
| C1 | 12 | 12 | — (the only difference start vs head is C1's intended highlight on the `">=0.95"` row, which is C1's gate; head equals C1 alone; with C1 reverse-applied, identical to start) |
| C2 | 366 | 349 | 17 cases differ only in C2's two IRD warning texts (its intended change; 8 of 336 validator-grid cells, exactly {IRD, threshold > 1, user-set}); fits, recorded thresholds and every other warning identical. With C2 reverse-applied, 366/366 identical to start. |
| C3 | 37 | 9 | 16 cases (binary and continuous fits with `adverse_outcome = FALSE`) differ only in the effect column: **A1** (OR → 1/OR, MD → −MD, RD → −RD); with A1 reverse-applied they are identical to start. 12 count cases are C3's own intended note change. With A1 and C3 reverse-applied, all identical to start. |

- No unexplained difference, and nothing was repaired. The only cross-fix effects are between A1 and C3, which both change `sg_tables()`:
  - C3 changes the GRF count footnote that A1's cases print.
  - A1 changes the effect column that C3's binary and continuous `adverse_outcome = FALSE` cases print.

### 7.3 Articles (step 2) — not re-rendered

One read-only agent per code fix listed the evaluated chunks of all 13 documents under `vignettes/` at `835fac1b`, followed child documents and `source()`, and traced the call graph to the fix's changed lines. Where a document reaches them, it ran the chunks with the article's own arguments and seed on both trees, sequentially and without rendering. Output: `$S/articles/<ID>.md`.

| Fix | Documents swept | Reaches changed lines | Calls and branch | Can the output change |
|---|---|---|---|---|
| A1 | 13 | no | — | no |
| A2 | 13 | no | — | no |
| A3 | 13 | yes: `vignettes/articles/mrct-consistency.qmd` | Chunk `frontier`: `compute_frontier_cis(fs, n_splits = 1000L, seed = seed_rep5)` on the survival fit of chunk `case-fs`, which stores `adverse_outcome = FALSE` (`forestsearch()` resolves it before capturing `args_call_all`). The new `is.null()` block (head `R/frontier_cis.R:144-147`) runs, the fallback is not taken, and the Cox estimator does not read the value. | No. Chunks `case-read`, `case-region`, `case-fs` and `frontier` were run with the article's own arguments and seed (613800005) on both trees in separate processes: `ci_front` (15 × 12) `identical()`; plot layer data identical (labels identical with the ggplot `S7_class` closure dropped); fits identical apart from 5 timing leaves. |
| B1 | 13 | no | — | no |
| C1 | 13 | no | — | no |
| C2 | 13 | no | — | no |
| C3 | 13 | no | — | no |

- **Only A3 reaches its changed lines, in `mrct-consistency.qmd`, and the output cannot change.** No article's output can change through any part-5 code fix.
- The documentation commits change no evaluated code. C7 changes prose outside code chunks in two articles (§5).
- **`quarto/`** is read-only and was not rendered; these findings are recorded, not gated.
  - A1: ten files call `sg_tables()` on a `grf.subg.harm.glm()` result built at that function's default `adverse_outcome = FALSE`, so their effect column would change if re-rendered (MD changes sign; binary OR becomes the OR of the non-event). They are `applications/actg175/template_actg175_continuous.qmd:182`, the only one outside `_archive/` and `_legacy/`; six continuous and two binary files under `applications/actg175/_archive/`; and `simulations/_legacy/actg175_binary_harm_simulations_1simexample.qmd:292`.
  - C1: `resampling/test_bootstrap_grf_cuts.qmd:314`, `:401` calls `forestsearch:::format_subgroup_summary_tables()`, so its consistency table now carries the highlight, but the document never displays that table.
  - A2, A3, B1, C2, C3: no `quarto/` file would change.
  - The full notes are in the detail file.

### 7.4 C5 examples at the head

The head was installed into `$S/lib/head` (`R CMD INSTALL --no-test-load`, 17.96 s), and the 16 examples were extracted from the head's `man/` with `tools::Rd2ex(commentDontrun = FALSE)`. The extracted code is identical to what the C5 agent ran. The 15 that RAN in C5 also run at the head (rc 0, walls within 1.5 s of C5's). Their output is identical to C5's once random gt table ids and one k-fold "Time:" line are masked. `fs_fdr_report` is NOT RUN, as in C5 (over the cap). No part-5 code fix changes any of these displays. Output: `$S/examples_head/RESULT.md`.

### 7.5 Full suite and `check()` (steps 4 and 5)

| Run | Part 4 (`REPORT_cran_prep_code_fixes_2026-11-12.md` §6) | This branch | Result |
|---|---|---|---|
| Full suite (`devtools::test()`, `NOT_CRAN=true`, in the repository at `fb80dd8f`) | 703 blocks, 102 files, 8,853 expectations, 0 failures, 0 errors, 3 skips, 38 warnings; 730 s | **736 blocks, 109 files, 9,242 expectations, 0 failures, 0 errors, 3 skips, 38 warnings**; wall **754 s** (17:10:36–17:23:10; cap 30 min) | **PASS** |
| `devtools::check(document = FALSE, error_on = "never")` | 0 errors / 0 warnings / 0 notes; 1,092 s | **0 errors / 0 warnings / 0 notes**, Status OK; wall **1,103 s** (`check()` 17:10:36–17:28:59, after the export that began at 17:10:30; tests 11 m; examples 38 s; `--run-donttest` 67 s; vignettes 76 s; cap 60 min) | **PASS** |

- **Difference, accounted for** (`REPORT_cran_prep_displays_se_docs_2026-11-12_suite_diff.txt`, per file against part 4's saved results): +33 blocks, +7 files, +389 expectations. These are exactly the seven new test files:

  | File | Blocks | Expectations |
  |---|---|---|
  | `test-consistency-highlight-labels.R` (C1) | 2 | 7 |
  | `test-cv-tables-orientation.R` (A2) | 5 | 111 |
  | `test-frontier-cis-orientation.R` (A3) | 4 | 31 |
  | `test-grf-footnote-scale.R` (C3) | 6 | 23 |
  | `test-rd-adjusted-se.R` (B1) | 7 | 35 |
  | `test-sg-tables-orientation.R` (A1) | 6 | 155 |
  | `test-threshold-remap-messages.R` (C2) | 3 | 27 |

  No other file changed counts. The 3 skips are the same (`test-cross-outcome-parity.R` ×2, `test-cv-no-subgroup-edges.R`), and warnings are unchanged at 38.
- **Check route** (`REPORT_cran_prep_2026-11-11.md` §6 with the `-z` correction): an export of the working tree at `fb80dd8f` (`git ls-files -z -co --exclude-standard`; 9,715 files listed, 9,715 copied) into `$S/check/forestsearch`. All 79 files changed since the start SHA were `cmp`-identical to the repository's. `forestsearch_0.4.1.9011.tar.gz` was built there. The local `legacy_vignette/` (git-excluded) was not exported, so its known NOTE did not arise.
- **Contention and order.** The suite was started right after the check export finished, so the export could not pick up files the suite writes. The two then ran together, so both wall clocks include contention. The suite left the working tree clean (`git status --porcelain` empty after it).
- **The NEWS wording commit `2a8be3bc` came after both runs.** It changes `NEWS.md` only, which no test, Rd page or code reads, so neither run was repeated (§8).

## 8. Departures and decisions

- **Polish before commit (comments only).** The verifiers' minor findings on comments and help text were fixed before each commit, not left open as in part 4. Each edit was comment-, roxygen- or prose-only, was checked by an independent agent, and the comment-only gate passed against the verified patch. This is the part-4 practice for F3/F11, applied more widely; part 4's own leftover requests were the reason for item C6. The verified patches are kept as `$S/patches/<ID>.verified.patch`.
- **Commit order.** C2 was committed first: it had no polish pending. The other code fixes followed in task order, then C4 to C7. All fixes are independent, and every patch applied cleanly.
- **A1, the GRF branch.** A `grf.subg.harm.glm()` result has no `args_call_all`; it stores `adverse_outcome` at the top level. The GRF branch therefore reads `fs$adverse_outcome`, which is still "the fit's stored `adverse_outcome`" as the task says. The default of `grf.subg.harm.glm()` is `FALSE` for every outcome type, so binary GRF results at their default also change (OR of the non-event).
- **A2, scope.** Only a `forestsearch_Kfold()` result builds these tables. A `forestsearch_tenfold()` result has no `cv_args`, and `forestsearch_tenfold()` calls `forestsearch_KfoldOut(outall = FALSE)`. Count and binary IRR/IRD rows cannot change, because the rate estimator takes no `adverse_outcome`.
- **B1, the safety rule.** With the literal rule, a slice where the treatment effect is not estimable would get `sqrt(g' V g)` = 0, because the aliased treatment column drops out. That value is finite and non-negative, so the rule would accept it, and an NA interval would become zero width. In the repair round the helper was changed to return NA whenever a column that treatment changes is aliased, so the existing expression, NA there, is kept as at the start. As a result, an aliased user-supplied treatment-interaction column with a nonzero constant also keeps the old value.
- **B1, other notes.** The task asked for B1's hunks inside `.estimate_rd()` plus the helper. The one-line helper-list comment in the file header (`R/glm_effect_estimators.R:19-23`) was added at the polish stage, comment only. The tolerance-reason comment of test (b) was edited after the first green run in round 1; the tolerance did not change, and this was disclosed.
- **C2.** The validator message names the measure and its ratio counterpart but no numeric default. It can fire for `outcome_type = "survival"` with a supplied `effect_measure = "RD"`/`"IRD"`, where no remap follows, so a default there would be wrong. The remap message names the default it assigns. The RD messages are byte-identical.
- **C3.** Besides naming the transform, the note adds ", not offset-adjusted" when an offset is recorded, and "(scale not recorded)" for a count object without a recognised `grf_y_transform`, which only a hand-built object can be. The note names the scale and not the sign. Negation changes orientation, not scale, and the object does not record whether the count was negated.
- **C4, sentences that differ from the task's wording.** The full list is in §5 (C4 table) and the detail file. Three of them:
  - The task says the `use_dina` "same cut" claim holds for the default selection statistic. The code's default `dina_select_statistic` is `"effect"`, and the claim holds only under `"dina"`; the help now says so.
  - Only values greater than 1 are remapped.
  - The `setup_gbsg_dgm()` HR formulas hold for patient-level hazard ratios and `AHR_*_true`, not for the Cox `hr_*_true`.
- **C4, edits beyond the listed sites.** One `find_k_inter_for_target_hr()` details sentence on failure, now stop-with-error. One non-ASCII character in `generate_bootstrap_with_noise()`'s details, so the page passes `checkRd()` with 0 problems. The `use_dina` screening `#` comment, which made the same claim as site 2.
- **C5.**
  - `n_workers = 100` was dropped from the `fs_fdr_report` example. It exceeds `parallelly`'s limit on machines with fewer than 34 cores.
  - The `summarize_simulation_results` example passes `confounders_base` as `paper_simulations.qmd` does: with the default, every replicate fails inside on `setup_gbsg_dgm()` data.
  - `fs_fdr_report` is NOT RUN (over the cap). Its lines ran to the end at `nsims = 4, n_workers = 4` (63.2 s), as a smoke check, not the example.
- **C6.** The verifier's optional F6 comments were added too. Three existing comment blocks were rewritten rather than appended to (F12 `:27-30`, F14 `:159-161`, one F6 line).
- **C7.**
  - One sentence was added to `biomarker_effects.qmd`, mapping the removed function, which the cited paper names, to `setup_gbsg_dgm()`.
  - In `paper_simulations.qmd`, the summary table's row for new GBSG work now recommends `generate_aft_dgm_flex()`, and "binary discretisation" was corrected to quartile groups. Both are in the same passage, and both say what the code does.
- **A second NEWS commit.** The task puts version and NEWS in one commit (`fb80dd8f`). The report drafting checked each NEWS entry against the code and found inexact claims; one further independent agent checked the rewording. Examples: the K-fold bullet named `forestsearch_tenfold()`; the IRD bullet said the validator names "the IRR scale"; the documentation entry said sixteen examples run, when fifteen were run to the end. The entries were corrected in `2a8be3bc`, `NEWS.md` only, wording only, rather than left wrong. History was not rewritten. The suite and `check()` ran at `fb80dd8f`; `2a8be3bc` changes no code, test or Rd file.
- **Unaffected-case normalisation.** A1 and C3 compare the display data (`_data` of every gt table) and whole objects with closures compared by their code; gt `_formats` closure environments differ between any two runs of the same tree. A3, B1, C1 and C2 compare whole objects (C2 drops the two timing leaves of each fit). The ggplot `S7_class` attribute is dropped as in part 4.
- **Report split.** The per-item evidence is long, so it is in `REPORT_cran_prep_displays_se_docs_2026-11-12_detail.md`, as in part 4.

## 9. OPEN ITEMS

Limited to what this task leaves open: B2, the scope of the fixes, and faults found while doing them. None of these changes the gates.

**Left open by this task's items**

1. **B2 (BLOCKED).** The `dr_gcomp` rate-difference SE stays `abs(mean(pred1)) * se(bT)` (head `R/glm_effect_estimators.R:945`). A correct SE has to include the propensity model's variability, through joint estimating equations or by resampling the whole two-stage procedure; for GRF or LASSO scores there is no closed-form `V`. The same branch is reached by `outcome_type = "binary"` with `effect_measure = "IRD"`. Related branches, not changed:
   - the `dr_gcomp` IRR branch reports the model-based `se(bT)` for `log(mean(pred1)/mean(pred0))`;
   - the binary `.estimate_or()` `dr_gcomp` branch does the same;
   - both `dr_gcomp` count branches use the model-based `vcov()` whatever `robust_se` is.
2. **B1, cases that keep the old expression.**
   - When both arms are at 0% or 100% the adjusted SE stays near 0 (about 2.3e-07 to 2.4e-07, the same on both trees); the unadjusted SE is exactly 0 there (F13).
   - An aliased treatment-interaction column with a nonzero constant keeps the old value (§8).
   - The tier-1 `invokeRestart("muffleWarning")` inside `tryCatch()` (`.estimate_rd()`, `.estimate_rr()`) still stops on a non-"converge" warning instead of falling back (code audit A-056).
3. **Orientation, outside this task's displays.**
   - The count and rate-measure estimator takes no `adverse_outcome`, by design, so count and IRR/IRD columns are never oriented.
   - The simulation extractor `.extract_fs_estimates_gen()` (head `R/run_simulation_analysis.R`, from the A4 sweep) builds raw GLM effects without orientation for `run_simulation_analysis()` and its summaries; listed, not changed.
   - The GRF CATE note prints treated-minus-control on GRF's working outcome, so with `adverse_outcome = TRUE` its sign is opposite to a difference on the named scale. Since A1, the note and the effect column have opposite signs in both orientations (before A1, only with `TRUE`): for example a continuous result at the default `FALSE` shows MD 0.558 beside H^ = -0.5838 (`$S/evidence/A1/grf_footnote.txt`). `?sg_tables` does not explain the note. The continuous note still says "(risk-difference scale)", which the task kept unchanged.
   - `?compute_frontier_cis` does not say the estimate takes the fit's orientation.
4. **C2 comments.** The comments at head `R/forestsearch_main.R:296` and `:2369` still say rate differences are bounded in [-1, 1]. Check 3 of the validator suggests `effect_measure = "RD"` for a small IRR threshold on a count outcome. Both are outside C2's two messages.
5. **C5.**
   - `fs_fdr_report`'s example is NOT RUN (over the cap).
   - Three GBSG examples outside the C5 list stop with "undefined columns selected", because their call has no `id.name`: `summarize_bootstrap_results`, `forestsearch_KfoldOut` and `summarize_cv_results`.
   - `plot_sg_results()` always returns `plots = list()`, so `plot.fs_sg_plot` prints "No plots available", while its `@return` describes a different message.
   - `R/subgroup_consistency_helpers.R:1223` at the start (head `:1229`, in a block ending `@noRd` at head `:1247`) has the bare-`gbsg` pattern.
6. **C7.**
   - `vignettes/articles/biomarker_effects.html` still shows `create_gbsg_dgm()` until the site is re-rendered.
   - The chunk label `generate-aft-equivalent` and an in-chunk comment are inside a code chunk, so the gate left them.
   - The quoted hazard ratios (about 4.4 and 4.1) go stale if either builder changes.
   - `paper_simulations.qmd:358-360` (start and head) holds a stray `---`/`title:` block.
7. **Articles not re-rendered.** No article's output can change (§7.3). The ten `quarto/` files of §7.3 (A1) would show a different effect column if re-rendered.

**Found while working; pre-existing; not changed**

8. **C4 documentation of code that is itself wrong.**
   - `test_hte_from_forestsearch()`'s `source` string "confounders.evaluated (GRF/LASSO selected)" (head `R/test_hte_crump.R:447`) is wrong for a default fit, and so is `?forestsearch`'s `confounders.evaluated` entry.
   - `cv_metrics_tables()` prints "Exact match on all factors" (head `R/cv_summary_tables.R:434`), and `CV_sgs()` gives 0 for a full-data subgroup of three or more factors.
   - With `est.scale = "1/hr"` on a GLM fit, the subgroup rows' arm columns and `Diff` are computed on swapped treatment while the effect is not. Separately, `mr_estimates_table(est.scale = "1/hr")` inverts the H cell but not the complement's naive cell.
   - `compute_dgm_cde()`'s warning names `df_super_rand` when `df_super` was read.
   - `args_call_all` keeps the supplied `hr.threshold` beside the remapped `effect.threshold`.
   - The `dmin.grf` `missing()` reset has no effect.
9. **C4, help text outside the listed sites.**
   - `find_k_inter_for_target_hr()` "typically converges within 5-10 iterations" (head `R/find_k_inter_main.R:52`; 11 and 19 observed).
   - `calibrate_k_treat()` `@param k_treat_range` "(uniroot will error otherwise)" (head `R/calibrate_k_treat.R:30`); `extendInt = "yes"` widens the range.
   - `?plot_sg_glm_outcomes` `@param outcome_type` omits count.
   - `R/forestsearch_helpers.R:1607-1610` calls `"dina"` the default statistic.
   - `estimate_propensity_scores()` documents `trim` as quantile bounds; the code clamps at absolute values.
10. **Displays and summaries.**
    - `compute_frontier_cis()` and `forestsearch_KfoldOut()` fall back to `"RD"` for a binary object with no stored `effect_measure` (`forestsearch()` uses `"OR"`).
    - `.label_to_expr()` turns a frontier label on a multi-level factor (`"!{stage}"`) into `stage == 0`, so `compute_frontier_cis()` returns an all-NA row for that member with no warning.
    - Count CV tables print IRR "NA", because `forestsearch_Kfold()` does not carry the offset into `resCV`.
    - `format_subgroup_summary_tables()` Table 5 highlights `rows = 2:4` by position, which are not the match rows.
    - `sg_tables()` stops on a DINA continuous fit whose consistency result has 0 rows.
    - `generate_bootstrap_with_noise()` on a one-column data frame returns a bare vector with no noise.
11. **`run_simulation_analysis()`'s default `confounders_base`** (`v1`–`v7`) are factors in `setup_gbsg_dgm()` data, so every FS replicate fails inside `get_FSdata` (reported as a warning; `any.H` 0). The value labelled RD there is the log OR, IRD the log IRR and RR the OR, so its `ahr.*` and `hr.*` columns disagree for those measures (A4 sweep).
12. **Count/IRD fits with `fs.splits = 10`, `use_twostage = TRUE`** failed inside `subgroup.consistency()` ("argument must be coercible to non-negative integer") on the start tree in one verifier probe. Not investigated; not caused by this task.

## 10. Final gates (step 8)

`REPORT_cran_prep_displays_se_docs_2026-11-12_final_gates.sh` recomputes the gates from git alone:
- G1: the working tree.
- G2: both sibling repositories (branch, SHA, empty status, no file newer than the start commit).
- G3: `git diff b06de4e8 --stat` lists only files this report names, under `vignettes/` only the two C7 files, and nothing under `quarto/`.
- G4: each commit touches only its own files. The lists for the item commits were taken from the gated patches, not from the commits; the start, version/NEWS and NEWS-wording commits' lists are those the task's rules fix.
- G5: nothing pushed or merged.
- G6: the documentation gates, recomputed per commit against its parent. That is comment-only for every `R/` and `tests/` file of C4, C5 and C6, and for C7 no changed line inside a code chunk or on a fence.

The script was run on the staged tree just before the report commit (output `..._final_gates.txt`, committed) and again after the commit.

| Gate | Result on the staged tree |
|---|---|
| G1 | PASS: only the 9 report files are staged; nothing else modified or untracked |
| G2 | PASS: fs-post-selection `main` `5d165acc…`, fs-glms-interpretable `main` `dc94a17e…`, both `git status --porcelain` empty, 0 files newer than the start commit (2026-11-12 15:45:11) |
| G3 | PASS: 88 paths in `git diff b06de4e8 --stat` including the staged report files; every one is named in this report or its detail file; under `vignettes/` only the two C7 files; nothing under `quarto/` |
| G4 | PASS: each of the 14 commits since the start, from `b06de4e8` to `2a8be3bc`, touches exactly its own files (the 11 item commits as their gated patches; `b06de4e8` the task file; `fb80dd8f` `DESCRIPTION` + `NEWS.md`; `2a8be3bc` `NEWS.md`) |
| G5 | PASS: HEAD is on no remote-tracking ref; 0 merge commits since `b06de4e8` |
| G6 | PASS: `ba7c4ea2` comment-only in 18 `R/` and `tests/` files, `5073fd67` in 11, `ffc55486` in 9; `835fac1b` 66 changed lines (old + new sides), none inside a code chunk or on a fence |

- With every report file staged (a placeholder in `_final_gates.txt` the first time), two runs gave the same output, all PASS; that output is the committed `REPORT_cran_prep_displays_se_docs_2026-11-12_final_gates.txt`. The gates were run again after the commit, when G1 requires an empty `git status --porcelain`.

## 11. Files

**Changed since the start SHA by the items and the two NEWS commits** (`git diff b06de4e8 2a8be3bc --stat`):

| Kind | Files (item) |
|---|---|
| R (27) | `R/bootstrap_summaries_helpers.R` (C5), `R/calibrate_k_treat.R` (C4), `R/cv_summary_tables.R` (C4), `R/dina_subgroup.R` (C4), `R/find_k_inter_main.R` (C4), `R/forestsearch_cross_validation.R` (A2, C4, C5), `R/forestsearch_main.R` (C2, C4), `R/format_subgroup_summary_tables.R` (C1), `R/frontier_cis.R` (A3), `R/fs_fdr_report.R` (C5), `R/get_fsdata.R` (C4), `R/glm_effect_estimators.R` (B1, C4), `R/mr_estimates_table.R` (C4), `R/oc_analyses.R` (C4, C5), `R/pareto_frontier_table.R` (C5), `R/plot_sg_distribution.R` (C5), `R/plot_sg_glm_outcomes.R` (C4), `R/plot_sg_results.R` (C5), `R/plot_subgroup_results_forestplot.R` (C4), `R/render_forestplot.R` (C5), `R/run_subgroup_sims.R` (C5), `R/setup_gbsg_dgm.R` (C4), `R/subgroup_consistency_helpers.R` (C4, C5), `R/summarize_bootstrap_subgroups.R` (C4), `R/summary_utility_functions.R` (A1, C3, C4, C5), `R/synthetic_data_perturbation.R` (C4), `R/test_hte_crump.R` (C4) |
| man (32) | `man/CV_sgs.Rd` (C4), `man/calibrate_k_treat.Rd` (C4), `man/compute_dgm_cde.Rd` (C4), `man/cv_summary_tables.Rd` (C4), `man/dot-glm_resample_supported.Rd` (C4), `man/evaluate_consistency_twostage.Rd` (C5), `man/evaluate_subgroup_consistency.Rd` (C5), `man/find_k_inter_for_target_hr.Rd` (C4), `man/forestsearch.Rd` (C4), `man/forestsearch_KfoldOut.Rd` (C4), `man/format_bootstrap_table.Rd` (C5), `man/fs_fdr_report.Rd` (C5), `man/generate_bootstrap_with_noise.Rd` (C4), `man/make_effect_estimator.Rd` (C4), `man/mr_estimates_table.Rd` (C4), `man/pareto_frontier_table.Rd` (C5), `man/plot.fs_forestplot.Rd` (C5), `man/plot.fs_sg_plot.Rd` (C5), `man/plot_sg_distribution.Rd` (C5), `man/plot_sg_glm_outcomes.Rd` (C4), `man/plot_subgroup_results_forestplot.Rd` (C4), `man/print.fs_forestplot.Rd` (C5), `man/print.fs_kfold.Rd` (C5), `man/print.fs_sg_plot.Rd` (C5), `man/print.fs_tenfold.Rd` (C5), `man/run_subgroup_sims.Rd` (C5), `man/setup_gbsg_dgm.Rd` (C4), `man/sg_consistency_out.Rd` (C5), `man/sg_tables.Rd` (C4, C5), `man/summarize_bootstrap_subgroups.Rd` (C4), `man/summarize_simulation_results.Rd` (C5), `man/test_hte_from_forestsearch.Rd` (C4) |
| tests (16) | `tests/testthat/test-bootstrap-summary-bins.R` (C6), `tests/testthat/test-consistency-highlight-labels.R` (C1), `tests/testthat/test-cv-tables-orientation.R` (A2), `tests/testthat/test-dina-gaussian-floor.R` (C6), `tests/testthat/test-display-orientation.R` (C6), `tests/testthat/test-duplicate-cut-names.R` (C6), `tests/testthat/test-frontier-cis-orientation.R` (A3), `tests/testthat/test-grf-count-orientation.R` (C6), `tests/testthat/test-grf-footnote-scale.R` (C3), `tests/testthat/test-hte-crump-confounders.R` (C6), `tests/testthat/test-lasso-intercept.R` (C6), `tests/testthat/test-mrct-default-confounders.R` (C6), `tests/testthat/test-rd-adjusted-se.R` (B1), `tests/testthat/test-rd-standard-error.R` (C6), `tests/testthat/test-sg-tables-orientation.R` (A1), `tests/testthat/test-threshold-remap-messages.R` (C2) |
| vignettes (2) | `vignettes/articles/biomarker_effects.qmd` (C7), `vignettes/articles/paper_simulations.qmd` (C7) |
| Version and NEWS (2) | `DESCRIPTION` (version), `NEWS.md` (`fb80dd8f`, `2a8be3bc`) |

The task file `dev/tasks/TASK_cran_prep_displays_se_docs_2026-11-12.md` is in the start commit itself.

**Beside this report** (committed with it):

| File | What it is |
|---|---|
| `REPORT_cran_prep_displays_se_docs_2026-11-12.md` | this report (also copied to `~/Downloads`) |
| `REPORT_cran_prep_displays_se_docs_2026-11-12_detail.md` | the full per-item evidence: sites, every red/green run, every expected value, unaffected cases, existing tests, repair and polish rounds, verifier problems, departures, the A4 sweep, the B2 form, C4 old/new text, the C5 runs, the cumulative and articles checks (also copied to `~/Downloads`) |
| `REPORT_cran_prep_displays_se_docs_2026-11-12_apply_fix.sh` | main-session per-fix gate: test file only → red; code → green |
| `REPORT_cran_prep_displays_se_docs_2026-11-12_commit_fix.sh` | main-session per-fix commit: staged diff equals the gated patch; `--stat`; commit |
| `REPORT_cran_prep_displays_se_docs_2026-11-12_comment_gate.sh` | the task's comment-only gate over every changed `R/` and `tests/` file |
| `REPORT_cran_prep_displays_se_docs_2026-11-12_gate_compare.txt` | default-path gate output (§7.1); the scripts are part 4's, unchanged |
| `REPORT_cran_prep_displays_se_docs_2026-11-12_suite_diff.txt` | per-file suite counts against part 4 (§7.5) |
| `REPORT_cran_prep_displays_se_docs_2026-11-12_final_gates.sh`, `REPORT_cran_prep_displays_se_docs_2026-11-12_final_gates.txt` | final gates (§10) and their output on the staged tree |

**Scratch, not in the repository** (`/tmp/cran_prep_displays_se_docs_2026-11-12/`):
- `start/`, the start export; `head_fixes/`, the head export;
- `bin/` (`mkcopy.sh`, the run scripts);
- `wt/<ID>/`, the scratch copies; `evidence/<ID>/`, implementer and verifier runs and notes; `verify/`;
- `patches/`, the gated patches plus `*.verified.patch` and `*.r1.patch`;
- `main/<ID>/`, the main-session red/green;
- `gate/`, `cumulative/`, `articles/`, `examples_head/`, `suite/`, `check/`;
- `report_sections/`, the drafts of this report and the detail file;
- the workflow results (`workflow1_result.json`, `polish_result.json`, `after_fixes_result.json`, `drafts_result.json`).
