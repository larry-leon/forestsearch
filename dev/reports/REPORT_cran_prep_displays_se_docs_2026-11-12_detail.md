# REPORT detail — CRAN preparation, part 5 (2026-11-12)

Detail file for `REPORT_cran_prep_displays_se_docs_2026-11-12.md`: the full per-item evidence, drafted from the scratch evidence (`/tmp/cran_prep_displays_se_docs_2026-11-12/`, `$S` below) and checked against git. Line numbers are at the start SHA `b06de4e8` and at head `fb80dd8f` unless stated. Where an earlier evidence file says something that a later repair, polish or main-session correction superseded, the later text is what was committed. NEWS remarks refer to `fb80dd8f`; its wording was corrected in `2a8be3bc` (main report §8).

## A1 — `sg_tables()` uses the fit's own orientation

`S` = `/tmp/cran_prep_displays_se_docs_2026-11-12`. Line numbers are at the start SHA `b06de4e8`, at the commit `f4c7dd2d` and at head `fb80dd8f`. `R/` is identical at `835fac1b` and `fb80dd8f`.

- The test file is the same blob (`d9dc110`) at `f4c7dd2d` and `fb80dd8f`, so its line numbers at the commit and at head are the same.
- `R/summary_utility_functions.R` changed after A1 in three commits: C4 `ba7c4ea2` (+20 roxygen lines at `:733`), C5 `5073fd67` (+2 net at `:759`) and C3 `c3323a59` (+23 net near `:925`). So A1's GRF-branch lines are 22 lower at head than at the commit, and its main-branch lines are 45 lower.

### A1.1 Status, commit, provenance

| | |
|---|---|
| Status | **PASS** |
| Commit | `f4c7dd2d` "fix: sg_tables() uses the fit's own orientation" |
| Parent | `2984dbee` (C2). `R/summary_utility_functions.R` there is the same blob as at the start SHA (`b2de758d`). |
| Files (`git diff --stat 2984dbee f4c7dd2d`, `S/main/A1/diffstat.txt`) | `R/summary_utility_functions.R` +22/-0; `tests/testthat/test-sg-tables-orientation.R` +268 (new). 2 files, +290/-0. No `man/` change. |
| Commit vs gated patch | `git diff 2984dbee f4c7dd2d` equals `S/patches/A1.patch` once the `index` lines are dropped. Those lines name the same blobs (`b2de758`→`10fc287`, new `d9dc110`) at a different abbreviation length. `S/main/A1/patch_equal.txt` is empty. |
| Earlier patch | `S/patches/A1.verified.patch` (the verified patch, before polish): R +20, test +248. |
| Committed tree | The committed R file and test file equal those in `S/wt/A1_polish/` (polish plus the main-session correction). Checked here with `diff`. |
| Workflow | `S/workflow1_result.json`: status READY, final verdict CONFIRMED, `repaired: false`, no blocking problem, 7 minor. `S/polish_result.json`: polish DONE. Its check returned PROBLEMS with one problem, which the main session corrected. |

### A1.2 Site, before → after

Function `sg_tables()`, `R/summary_utility_functions.R`. The task's pointers (829–839 and 982–992) are exact at the start SHA.

**What reaches each branch** (code at `b06de4e8`, checked by implementer and verifier, and re-checked here):

- **GRF branch.** It starts at `if (inherits(fs, "grf_glm_result"))` (start `:811`; head `:833`). `fs` is a `grf.subg.harm.glm()` result.
  - It has no `args_call_all`. The test asserts this (`expect_null(g$args_call_all)`, `:187`).
  - It stores `adverse_outcome` at the top level in all four return paths: `R/grf_subg_harm_glm.R:551`, `:574`, `:647`, `:710`.
  - `grf.subg.harm.glm()`'s own default is `adverse_outcome = FALSE` for every outcome type (`:304`).
  - In this branch `outcome_type` is `fs$outcome_type %||% "binary"` (start `:817`; head `:839`).
- **Main branch.** It reads `args_fs <- fs$args_call_all` (start `:937`; head `:994`).
  - `forestsearch()` resolves `adverse_outcome` with the rule at `R/forestsearch_main.R:1840-1841`. That happens before `mget()` at `:1924`, so the stored value is the resolved one.
  - The changed block sits inside `if (is_glm)` (start `:976`; head `:1033`), so survival fits never reach it.
- **The closure.** `make_effect_estimator()` has the formal `adverse_outcome = TRUE` (start `R/glm_effect_estimators.R:107`).
  - It passes the value to the binary OR/RR/RD estimator (`:181`), which recodes `Y -> 1 - Y` under FALSE (`:390-391`).
  - It passes the value to the linear-model estimator (`:193`), which negates Y under FALSE (`:956-957`).
  - Count outcomes, and binary IRR/IRD, go to `.make_poisson_rate_estimator()` without it (`:162-172`, `:198` on).

**Before.** Neither closure passed `adverse_outcome`, so the formal `TRUE` applied whatever the fit used.

GRF branch, start `:828-839`:

```r
    # Build estimator closure
    estimator_fn_grf <- tryCatch(
      make_effect_estimator(
        outcome_type      = outcome_type,
        treat.name        = treat.name_grf,
        outcome.name      = outcome.name,
        offset.name       = offset.name,
        effect_measure    = effect_measure,
        adjust_covariates = fs$args_call_all$adjust_covariates
      ),
      error = function(e) NULL
    )
```

Main branch, start `:981-992`:

```r
    # Build estimator closure for the table
    estimator_fn_tab <- tryCatch(
      make_effect_estimator(
        outcome_type      = outcome_type,
        treat.name        = args_fs$treat.name,
        outcome.name      = args_fs$outcome.name,
        offset.name       = args_fs$offset.name,
        effect_measure    = effect_measure,
        adjust_covariates = args_fs$adjust_covariates
      ),
      error = function(e) NULL
    )
```

**After, GRF branch.** Head `:850-873`: comment `:850-855`, value `:856-859`, call `:862-873`, new argument `:869`. At the commit this block is `:828-851`, with the argument at `:847`.

```r
    # Orientation: the GRF result's own adverse_outcome (stored at the top
    # level; a grf_glm_result has no args_call_all), so the effect column has
    # the fit's own orientation (MD: the sign; binary OR, RR, RD: which level
    # is coded as the event; count, and binary IRR or IRD: no change, as the
    # rate estimator takes no adverse_outcome).  NULL only on an object not
    # built by grf.subg.harm.glm(); resolved by forestsearch()'s rule.
    adverse_outcome <- fs$adverse_outcome
    if (is.null(adverse_outcome)) {
      adverse_outcome <- outcome_type %in% c("binary", "count")
    }

    # Build estimator closure
    estimator_fn_grf <- tryCatch(
      make_effect_estimator(
        ...
        effect_measure    = effect_measure,
        adverse_outcome   = adverse_outcome,
        adjust_covariates = fs$args_call_all$adjust_covariates
      ),
```

**After, main branch.** Head `:1038-1059`: comment `:1038-1041`, value `:1042-1045`, call `:1048-1059`, new argument `:1055`. At the commit this block is `:993-1014`, with the argument at `:1010`.

```r
    # Orientation: the fit's resolved adverse_outcome, so the effect column
    # is the fit's own estimate, sign included (the arm means and Diff stay
    # raw).  NULL only on an object not built by forestsearch(); resolved by
    # forestsearch()'s rule.
    adverse_outcome <- args_fs$adverse_outcome
    if (is.null(adverse_outcome)) {
      adverse_outcome <- outcome_type %in% c("binary", "count")
    }
    ...
        effect_measure    = effect_measure,
        adverse_outcome   = adverse_outcome,
        adjust_covariates = args_fs$adjust_covariates
```

- The commit has 4 hunks. All are in the two closures or in the lines that resolve the value just before each.
- The pattern is part 4's F12 (start `R/plot_sg_glm_outcomes.R:166-172`, call `:268`). The one difference is that the GRF branch reads the top-level value.
- **Not touched:**
  - Table 2 (`sg10_out`).
  - The arm means and `Diff`.
  - The GRF CATE note (start `:917`, which is C3's site; head `:974`).
  - `adjust_covariates = fs$args_call_all$adjust_covariates` (start `:836`; head `:870`).

### A1.3 Red and green runs

**Command.** Implementer, verifier and polish ran `NOT_CRAN=true timeout 300 Rscript --vanilla scripts/run_test.R <tree> sg-tables-orientation`.

- The script calls `devtools::test(<tree>, filter = ..., reporter = "summary", stop_on_failure = FALSE)` with `testthat.summary.max_reports = 1000`, then prints a count table per block.
- The main session ran the same filter on the repository. Its `red.txt` uses testthat's default cap, so it lists only the first 10 failures ("Maximum number of 10 failures reached") and prints no totals.
- Its counts below come from the progress line: 155 marks, of which 41 are failure labels (`1`–`9`, `a`–`z`, `A`–`F`) and 114 are dots.

| Run | Tree | Expectations | Failed | Errors | Skipped | Per block (nb / failed) | Wall | File (under `S/`) |
|---|---|---|---|---|---|---|---|---|
| Main red | repository at `2984dbee` + the test file only | 155 | 41 | 0 | – | not printed | 5.8 s | `main/A1/red.txt`, `wall.txt` |
| Main green | repository at `f4c7dd2d` | 155 | 0 | 0 | – | – | 5.5 s | `main/A1/green.txt`, `wall.txt` |
| Implementer red | start copy + test file (verified version) | 155 | 41 | 0 | 0 | 13/0, 24/11, 23/10, 46/18, 45/0, 4/2 | 6.84 s | `evidence/A1/red.txt` |
| Implementer green | start copy + fix | 155 | 0 | 0 | 0 | 13/0, 24/0, 23/0, 46/0, 45/0, 4/0 | 5.85 s | `evidence/A1/green.txt` |
| Verifier red | `verify/A1/red` (start + `tests/*` of the patch) | 155 | 41 | 0 | 0 | 13/0, 24/11, 23/10, 46/18, 45/0, 4/2 | 6.83 s | `evidence/A1/verify/red.txt` |
| Verifier green | `verify/A1/fixed` (start + whole patch) | 155 | 0 | 0 | 0 | all 0 | 8.34 s | `evidence/A1/verify/green.txt` |
| Polish red | `wt/A1_polish_red` (start + polished `tests/` only) | 155 | 41 | 0 | 0 | 13/0, 24/11, 23/10, 46/18, 45/0, 4/2 | 7.03 s | `evidence/A1/polish/red_polish.txt` |
| Polish green | `wt/A1_polish` | 155 | 0 | 0 | 0 | all 0 | 5.93 s | `evidence/A1/polish/green_polish.txt` |
| Polish check green | `verify/A1_polish` | 155 | 0 | 0 | 0 | 13/24/23/46/45/4, all 0 | 5.98 s | `evidence/A1/verify_polish/green.txt` |

Blocks: (1) fixtures, 13; (2) `forestsearch()` fit, 24; (3) GRF result, 23; (4) no stored value, 46; (5) `adverse_outcome = TRUE`, 45; (6) mirror twins, 4.

- **First failure, main red:** `test-sg-tables-orientation.R:199:3`. This is block 2, where `.sgo_check()` is called at `:199` and `expect_lt()` sits at `:140`.
  - Message: "Expected MD est ITT < `5e-07 + 1e-12`. Actual comparison: 0.4573760 >= 0.0000005".
  - 0.4573760 = 2 × 0.228688. That is the sign flip.
  - The next listed failures are "MD sign ITT" (`actual` -1.0, `expected` 1.0; `:141`), "MD CI ITT < 1e-04" (0.4574; `:142`) and "MD est H" (2.1882242 = 2 × 1.094112).
- **Implementer and verifier:** same expectation and same message, at `:182:3` in the pre-polish test file.
  - The polish added 7 comment lines to the header and 10 to the fixture block, which moves it to `:199:3`.
  - The polish red run shows the same message at `:199:3`.
- **Red failures by block:**
  - `forestsearch()` fit, 11: 3 sets × (estimate, sign, CI), plus the two `max_sg_est` ties.
  - GRF result, 10: 9, plus the `-cate_sg` sign tie.
  - No stored value, 18: 9 per branch.
  - Mirror twins, 2.
  - The fixture block and the TRUE block pass on both trees, as they should. Every raw check (Mean(C), Mean(T), Diff) passes on both trees.
- **Replaced run:** the implementer's first red run used the default reporter cap (wall 6.36 s; NOTES §3). The run in the table replaced it: same tree, same test file.

### A1.4 Expected values

All are in `tests/testthat/test-sg-tables-orientation.R`; line numbers are at head (equal to the commit). None is copied from the fixed code's output. The verifier re-derived them all (`S/evidence/A1/verify/scripts/hand.R`, `hand.txt`).

| Value | How derived | Where (head) |
|---|---|---|
| Effect-cell estimate per set (ITT; Questionable = `treat.recommend == 0`; Recommend = `== 1`) | **By hand, from the documentation:** the unadjusted OLS treatment coefficient is `mean(y | T = 1) - mean(y | T = 0)`. Orientation follows `@param adverse_outcome` of `?forestsearch` and `?grf.subg.harm.glm` (TRUE = higher is worse, FALSE = higher is better) and the convention effect > 0 = harm. So under FALSE the value is `-(mean1 - mean0)`. | header `:21-30`; `.sgo_md()` `:102-111`; checks `:140` (estimate), `:141` (sign) |
| Effect-cell limits | By hand: estimate ± `qnorm(0.975)` × pooled SE `sqrt(s2 (1/n0 + 1/n1))`, with `s2 = (SS0 + SS1) / (n0 + n1 - 2)`. | header `:23-25`; `.sgo_md()` `:107-109`; `:137-138`, `:142` |
| Mean(C), Mean(T), Diff | By hand: the raw arm means and the raw `mean1 - mean0`, the same under either orientation. | header `:31`; `:143-148` |
| Tie to the fit, `forestsearch()` | `max_sg_est` equals the hand Questionable value (`expect_equal`, tolerance 1e-10) and is > 0. The printed Questionable cell is within 5e-7 + 1e-12 of `max_sg_est`, with the same sign. | fixture `:176-177` (comment `:166-175`); block 2 `:202-203`; TRUE block `:247` (comment `:244-246`) |
| Tie to the fit, GRF result | The sign of the Questionable cell equals the sign of the result's harm score `-cate_sg`. Only the sign is compared, because `cate_sg` is a doubly robust mean, not the arm difference. | header `:32-35`; `:217` |
| Premises (pass on both trees) | The `forestsearch()` fit stores FALSE by the rule (`:163`) and finds a subgroup (`:164`); `fs_t` stores TRUE (`:178`). The GRF result has class `grf_glm_result` (`:185`), stores FALSE as the GRF default (`:186`), has `args_call_all` NULL (`:187`), finds a subgroup (`:188`), and has `-cate_sg > 0` (`:189`). | fixture block `:154-192` |
| Fallback | The same hand values under FALSE once the stored value is removed (rule: continuous gives FALSE), in both branches. | `:220-234` |
| Mirror twins | An invariant, not an expected value. The FALSE fit on `-y` and the TRUE fit on `y` select the same subgroup (`:260-261`) and show identical effect columns (`:266-267`). | `:252-268` |

**Data.** The data is part 4's F12 construction: `.make_continuous_data(N = 300, MD_harm = 1.5, seed = 42)`, with `biomarker_hi` and `sex` made numeric and a mirror twin `ben` with `y` negated (`.sgo_data()` `:58-64`). The sets come from `fs$df.est` for a `forestsearch()` fit, or `fs$data` for a GRF result (`.sgo_sets()` `:115-117`).

**Hand values** (verifier `hand.txt`). The verifier used its own copy of the generator (`identical()` to the helper's) and checked each value against `lm()` on the oriented outcome (`stopifnot`, 1e-12):

| Fit | ITT | Questionable | Recommend | SE (ITT / Q / R) |
|---|---|---|---|---|
| `forestsearch()`, FALSE on `ben` (same for TRUE on `harm`) | 0.2286880 | 1.0941122 | 0.0957932 | 0.120 / 0.291 / 0.124 |
| `grf.subg.harm.glm()`, FALSE on `ben` (same for TRUE on `harm`) | 0.2286880 | 0.5580767 | -0.3413504 | 0.120 / 0.149 / 0.193 |

- **Fixed tree:** max |table − hand| is 4.3e-7 for estimates, 1.1e-5 for limits and 4.3e-7 for the raw columns.
- **Start tree:** the TRUE fits match the fixed tree exactly. The FALSE fits are off by 2 × |estimate| (2.188224 at Questionable).
- **`max_sg_est` vs hand Questionable:** the difference is 2.22e-16 for both `forestsearch()` fits. Three runs measured it: verifier `hand.txt`, polish `polish/measure_maxsg.R`, polish check `verify_polish/probe1.txt`.
- **Set membership:** `(age <= 47) & (biomarker_hi == 1)` evaluated on the raw columns equals `treat.recommend == 0`. TRUE for both `forestsearch()` fits.

**Tolerances**, set before the first run (header `:36-48`):

- **5e-7.** The table is printed at `ndecimals = 6`, so a cell is within half a unit of the last digit.
- **`+ 1e-12`.** Floating-point slack: the printed decimal is parsed back, the reference may differ by rounding, and `expect_lt()` is strict. The header has explained this since the polish (`:38-45`).
- **1e-4 for limits.** The table uses the literal `z = 1.96` (`analyze_subgroup_glm()`), while the hand value uses `qnorm(0.975)`.
  - The gap is 3.6e-5 × SE, which is below 2e-5 at SE < 0.5.
  - The verifier measured `(1.96 - qnorm(0.975)) × max SE` = 1.05e-5.
- **1e-10 on `max_sg_est`.** It is relative, as edition-3 `expect_equal()` applies it. The reason was added in the polish (`:166-175`, `:244-246`); the measured difference is 2.22e-16.
- **No changes after the first run.** No tolerance, seed or value changed (NOTES §4). The polish changed comments only. The polish check confirmed this with `parse(keep.source = FALSE)`: the code is identical, verified vs polished.

### A1.5 Unaffected case

The task's case: "survival fits, and non-survival fits with `adverse_outcome = TRUE`, give identical display data."

**Implementer** (`S/evidence/A1/scripts/unaffected_run.R`, `unaffected_compare.R`; output in `unaffected.txt`, `unaffected_compare.txt`):

- 13 fits per tree, each tree in its own `Rscript --vanilla`.
- `sg_tables()` ran at `ndecimals` 3 and 6.
- The comparison ran in a third process. A start-vs-start control ran the same way.

| Case | Branch | Stored `adverse_outcome` | Measure | Subgroup |
|---|---|---|---|---|
| M_surv | main | FALSE (not used for survival) | HR | `{age <= 55}` |
| M_cont_TRUE | main | TRUE | MD | `{age <= 47} & {biomarker_hi}` |
| M_bin_OR_TRUE | main | TRUE | OR | `{age <= 62} & !{biomarker <= -0.8}` |
| M_bin_OR_unset | main | TRUE (rule) | OR | same |
| M_bin_RD_TRUE | main | TRUE | RD | same |
| M_cnt_IRR_TRUE | main | TRUE | IRR | `!{biomarker <= -1}` |
| M_cnt_IRR_unset | main | TRUE (rule) | IRR | same |
| M_cnt_IRR_FALSE (extra) | main | FALSE | IRR | same |
| G_cont_TRUE | GRF | TRUE | MD | two-leaf union |
| G_bin_TRUE | GRF | TRUE | log_OR | two-leaf union |
| G_bin_RR_TRUE | GRF | TRUE | log_RR | two-leaf union |
| G_cnt_TRUE | GRF | TRUE | log_IRR | `age <= 48` |
| G_cnt_FALSE (extra) | GRF | FALSE (GRF default) | log_IRR | two-leaf union |

- **Display data.** This is `[["_data"]]` of every gt table returned (`tab_estimates`, `sg10_out`), at both `ndecimals`. `identical()` is TRUE in 13/13 start vs fixed, and in 13/13 start vs start.
- **Whole returned objects.**
  - `identical()` is FALSE for the 8 `forestsearch()` cases, in the start-vs-start control as well.
  - The only leaves that differ are 192 gt closures, `sg10_out$_formats[[i]]$func$*` (8 cases × 2 `ndecimals` × 2 formats × 6 renderers). Their environments differ between any two processes.
  - With each closure replaced by its deparsed code (environment dropped), whole objects are identical in 13/13, start vs fixed. Only closure environments were dropped.
- **Count cases with FALSE.** They are identical because the count route of `make_effect_estimator()` (and binary IRR/IRD) does not pass `adverse_outcome` to the Poisson estimator.
- **Walls:**
  - start 9.51 s (fits + load 7.889 s, tables 1.148 s);
  - fixed 9.68 s (8.021 s / 1.177 s);
  - start control 9.47 s (7.883 s / 1.113 s).

**Verifier** (`S/evidence/A1/verify/scripts/unaff_run.R`, `unaff_compare.R`; output in `unaff_run.txt`, `unaff_compare.txt`):

- 15 configurations of its own, with data seeds and settings different from the implementer's.
- `ndecimals` 3 and 6, plus three `subgroup_notation` runs.
- Three processes: start control, fixed, and a second start control.

| Case | Branch | Stored | Measure |
|---|---|---|---|
| V_surv_gbsg (`survival::gbsg`) | main | FALSE (not used) | HR |
| V_cont_TRUE_s7 | main | TRUE | MD |
| V_bin_RR_TRUE | main | TRUE | RR |
| V_bin_IRR_TRUE | main (binary → Poisson) | TRUE | IRR |
| V_bin_OR_unset | main | TRUE (rule) | OR |
| V_bin_OR_TRUE_adj (`adjust_covariates = "age"`) | main | TRUE | OR |
| V_cnt_IRD_TRUE | main | TRUE | IRD |
| V_cnt_IRR_TRUE_g (`use_grf = TRUE`) | main | TRUE | IRR |
| V_cont_TRUE_dina (`subgroup_method = "dina"`) | main | TRUE | MD (`sg_tables()` errors on both trees; A1.10, item 5) |
| V_bin_grf_inside (`fit$grf_res`) | GRF | TRUE | log_OR |
| G_cont_TRUE_front (`grf_selection = "frontier"`) | GRF | TRUE | MD |
| G_cont_TRUE_nosg (no subgroup found) | GRF | TRUE | MD |
| G_bin_OR_TRUE | GRF | TRUE | log_OR |
| G_bin_RD_TRUE | GRF | TRUE | RD |
| G_cnt_TRUE_id (`grf_count_transform = "identity"`) | GRF | TRUE | log_IRR |

- **Results.** Display data is identical in 15/15 at both `ndecimals` and in the 3 notation runs. Whole objects are identical in 15/15 after closures are replaced by their code. The gt `_formats` closure environments differ start vs start as well.
- **Walls:** start control 13.98 s, fixed 13.26 s, second start control 13.68 s.
- **Count check** (`changed.txt`): a count fit with FALSE gives the same effect column as with TRUE, on both trees.
- **Safety** (`safety.txt`): stored values `NA`, `logical(0)`, `"yes"` or `1`, in either branch, never make `sg_tables()` error.
  - Invalid values give `NA` cells, because both the construction and the closure call sit inside the existing `tryCatch`.
  - `1` behaves as TRUE.

### A1.6 Cases that change beyond the test (intended by the rule)

- Binary OR/RR/RD results that store FALSE now show the effect on `1 - Y`. Two kinds of fit do this:
  - `forestsearch(adverse_outcome = FALSE)` on a binary outcome (verifier `changed.txt`). Questionable OR goes from 0.728877 to 1.371974 (= 1/0.728877). That equals the hand OR of the adverse event and the TRUE fit on `y`.
  - `grf.subg.harm.glm()` on a binary outcome at its default FALSE. The verifier's run shows Questionable 0.376989 → 2.652597. The implementer's run (`binary_false_example.txt`) shows OR: ITT 1.057 → 0.946, Questionable 0.540 (0.288, 1.010) → 1.853 (0.990, 3.468), Recommend 2.306 → 0.434.
- Count displays never change (A1.5).

### A1.7 Existing tests

- No R function calls `sg_tables()`; at `b06de4e8`, `git grep` finds it only in roxygen and comments.
- The only existing test file that calls it is `tests/testthat/test-no-subgroup-dina-grf-summary.R`.

| File | Run by | Expectations | Failed | Errors | Skipped | Wall |
|---|---|---|---|---|---|---|
| `test-no-subgroup-dina-grf-summary.R` | implementer, fixed copy | 21 | 0 | 0 | 0 | 5.74 s |
| `test-no-subgroup-dina-grf-summary.R` | verifier, fixed copy | 21 | 0 | 0 | 0 | 5.68 s |
| `test-display-orientation.R` (part 4 F12; context run, not a changed function) | implementer | 93 | 0 | 0 | 0 | 7.07 s |
| `test-display-orientation.R` | verifier | 93 | 0 | 0 | 0 | 6.91 s |

- Expectations updated: **none**.
- Output: `S/evidence/A1/existing_tests.txt`, `S/evidence/A1/verify/existing_tests.txt`.

### A1.8 What a user sees

Source: `S/evidence/A1/example.txt`. The data is the test's `ben` (`y` negated). The output is `sg_tables(fit)$tab_estimates` at the default `ndecimals = 3`.

`forestsearch()` fit: `adverse_outcome` resolved to FALSE, `sg.harm = {age <= 47} & {biomarker_hi}`, `max_sg_est = 1.09411`.

| Row | Mean(C) | Mean(T) | Diff | `MD (95% CI)`, start tree | `MD (95% CI)`, fixed tree |
|---|---|---|---|---|---|
| ITT | 0.003 | -0.226 | -0.229 | `-0.229 (-0.464, 0.007)` | `0.229 (-0.007, 0.464)` |
| Questionable | -0.330 | -1.424 | -1.094 | `-1.094 (-1.665, -0.523)` | `1.094 (0.523, 1.665)` |
| Recommend | 0.045 | -0.051 | -0.096 | `-0.096 (-0.338, 0.146)` | `0.096 (-0.146, 0.338)` |

`grf.subg.harm.glm()` result on the same data: default FALSE, harm score `-cate_sg = 0.583784`.

| Row | Mean(C) | Mean(T) | Diff | start | fixed |
|---|---|---|---|---|---|
| ITT | 0.003 | -0.226 | -0.229 | `-0.229 (-0.464, 0.007)` | `0.229 (-0.007, 0.464)` |
| Questionable | 0.153 | -0.405 | -0.558 | `-0.558 (-0.850, -0.266)` | `0.558 (0.266, 0.850)` |
| Recommend | -0.252 | 0.089 | +0.341 | `0.341 (-0.037, 0.720)` | `-0.341 (-0.720, 0.037)` |

- Mean(C), Mean(T) and Diff are the same on both trees. So for these fits, Diff = −MD after the fix.

### A1.9 Articles

- `vignettes/articles/forestsearch-2024.Rmd:415` calls `sg_tables()` on a survival fit (`outcome_type = "survival"`, `:376`). The survival path does not reach the changed lines, so the article's output cannot change.
- No other file under `vignettes/` calls `sg_tables()` (`git grep` at `b06de4e8`).
- One read-only file outside `vignettes/` is affected: `quarto/applications/actg175/template_actg175_continuous.qmd`, which is not rendered.
  - At `:182` it calls `sg_tables(grf_result, ...)`.
  - `grf_result` is a continuous `grf.subg.harm.glm()` result at the default FALSE (`:147`, `:151`, `:157`).
  - Its effect column would flip if the template were re-rendered (verifier problem 5).

### A1.10 Verification, repair, polish

- **Verifier:** CONFIRMED (`S/evidence/A1/verify/VERIFY.md`). No blocking problem, 7 minor.
- **Repair round:** none (`workflow1_result.json`: `repaired: false`, `repair: null`).
- **Scope check (verifier):** its fixed copy re-diffs to `A1.verified.patch` byte for byte (`rediff.patch`). The patch is ASCII only and has no `man/` change.
- **Polish:** a comment-only edit for verifier items 1, 2 and 4.
  - The polish worked in `S/wt/A1_polish/`, which was `A1.verified.patch` committed and tagged `verified`.
  - The comment-only gate against `verified` printed nothing.
  - The polished test file gave red 155/41 (7.03 s) and green 155/0 (5.93 s), the same per-block counts as the verified runs.
- **Polish check** (`polish_result.json`, `check`): verdict PROBLEMS, one problem.
  - The reworded R comment ("binary: which level is coded as the event") overclaimed for count results.
  - The Poisson route takes no `adverse_outcome`, so a GRF count result shows the same IRR column at TRUE and FALSE (`verify_polish/probe2_count.txt`: identical TRUE).
  - Everything else checked out. The gate passed, the code is identical under `parse(keep.source = FALSE)`, and green was 155/0 (5.98 s).
- **Main-session correction:** the comment now names count, and binary IRR or IRD, as "no change, as the rate estimator takes no adverse_outcome" (head `:852-854`).
  - Comment only. NOTES records the gate re-run as empty.
  - The final delta from the verified patch is `S/evidence/A1/polish.diff`: 2 files, +25/-3.
  - The committed files equal `S/wt/A1_polish/`.

| # | Verifier problem (VERIFY §8) | What happened |
|---|---|---|
| 1 | `expect_equal(fit$max_sg_est, h$est, tolerance = 1e-10)` (verified test `:159`, `:227`) had no stated reason. | **Polish.** Comments added at head `:166-175` and `:244-246`. The polisher checked the reason against the code: `max_sg_est = max(hr.out$HR)` (`R/subgroup_search.R:1264`); the unadjusted fast path, `lm.fit(cbind(1, t), y)` (`R/glm_effect_estimators.R:1039-1042` at start); the relative tolerance of `waldo::num_equal()`. Measured difference: 2.22e-16. |
| 2 | The `+ 1e-12` added to 5e-7 (verified `:133`, `:141`, `:185`) was not explained. | **Polish.** Header `:38-45`. The check sites (head `:140`, `:148`, `:202`) are unchanged code. The polish check confirmed that all 48 six-decimal cells parse to the correctly rounded double (`verify_polish/cells6.txt`). |
| 3 | NOTES §2 has the GRF footnote at fixed `:926`; the string is at `:927` (the `sprintf(` call is at `:926`). | **Not corrected.** It is evidence only, and NOTES §2 still says `:926`. The string is at `:929` at the commit and `:974` at head (C3 changed it). |
| 4 | The R comment "so the effect column is on the fit's own scale" is loose: the fix aligns orientation, not scale. | **Polish, then main-session correction.** The polish wrote "orientation (MD: the sign; binary: which level is coded as the event)". The polish check found this wrong for count results. The main session added "count, and binary IRR or IRD: no change, as the rate estimator takes no adverse_outcome" (head `:850-855`). The main-branch comment (head `:1038-1041`) was not listed and is unchanged. |
| 5 | `quarto/applications/actg175/template_actg175_continuous.qmd:182` (read-only) would flip if re-rendered. | **Not changed** (`quarto/` is read-only). It is an open item for the articles statement (A1.9). |
| 6 | Pre-existing, outside A1, identical on both trees: `sg_tables()` stops at start `:1108` (`sg10$d1 <- NA`, "replacement has 1 row, data has 0") on a DINA continuous fit that found a subgroup but whose `grp.consistency$out_sg$result` has 0 rows (`verify/dina_probe.txt`). | **Not changed.** It is an open item. The line is at head `:1175`. |
| 7 | The NEWS draft "Diff has the opposite sign to the effect" is exact for MD and RD only. For OR/RR it means ratio > 1 where Diff < 0. | **The final NEWS bullet drops the phrase.** `NEWS.md:404-412` at head says only that Mean(C), Mean(T) and Diff "stay raw". C4's help text gives the relation per measure (`R/summary_utility_functions.R:748-752` at head). |

### A1.11 Departures

- **The GRF branch reads `fs$adverse_outcome`, not `fs$args_call_all$adverse_outcome`.**
  - A `grf_glm_result` has no `args_call_all` and stores the value at the top level (A1.2).
  - Following F12's pattern would give NULL and then the rule. That means TRUE for a binary GRF result that was fitted at FALSE.
  - The task's rule says "the fit's stored `adverse_outcome`" and names no slot (`TASK:80-82`). So this departs from F12's code, not from the task's words.
  - When GRF runs inside `forestsearch()`, the value is the parent fit's resolved one (`R/forestsearch_helpers.R:1895`, `R/grf_args.R:65`, `R/forestsearch_main.R:3155`).
- **The GRF-branch fallback is `forestsearch()`'s rule**, applied to `fs$outcome_type %||% "binary"`, as the task states.
  - Only an object not built by `grf.subg.harm.glm()` reaches it, because every return path stores the value.
  - `grf.subg.harm.glm()`'s own default would be FALSE. The rule and that default differ only for a binary or count object without the field; there the rule keeps today's TRUE.
- **The test ties the GRF result to `-cate_sg`, not `max_sg_est`.** A GRF result has no `max_sg_est`. The tie is on the sign only, because `cate_sg` is a doubly robust mean.
- **The test goes beyond the task's one case.** It adds the fallback in both branches, the TRUE case in both branches (which passes on both trees) and the mirror twins.
  - The `forestsearch()` fixture is F12's call without `mr_inference`, which `sg_tables()` does not use.
  - It selects the same subgroup and `max_sg_est` as F12's fit.
- **The unaffected case is extended.**
  - It adds binary and count with the value unset (rule gives TRUE), count with FALSE in both branches, and GRF binary `log_RR`.
  - "Identical" is judged on display data, and on whole objects after closure environments are dropped (A1.5).
- **Cases that change beyond continuous** (intended by the rule, not in the test): binary OR/RR/RD results that store FALSE (A1.6).

### A1.12 Open items

1. **The GRF branch's `adjust_covariates`.** `adjust_covariates = fs$args_call_all$adjust_covariates` (head `:870`) is always NULL for a `grf_glm_result`, so the estimate is unadjusted. This predates A1, is harmless, and is unchanged.
2. **The GRF CATE note's sign.** The note prints `cate_sg`, the treated-minus-control CATE on GRF's working outcome. Its sign is opposite to the effect column's "> 0 = harm".
   - Before A1 this was so only for TRUE results; after A1 it is so for both orientations. `grf_footnote.txt` shows H-hat = -0.5838 beside MD 0.558 in both.
   - C3 changed only the scale label for counts. Its comment says the sign is not named (head `:948-955`).
   - C4's help text for `sg_tables()` does not mention the note. Still open.
3. **The quarto ACTG175 continuous template** (A1.9) would show a flipped GRF effect column if re-rendered.
4. **The pre-existing DINA stop** in `sg_tables()` (head `:1175`; A1.10, item 6).
5. **The help page** (NOTES §9.3) is closed by C4 `ba7c4ea2`: `?sg_tables` `@details`, head `R/summary_utility_functions.R:736-754`.

### A1.13 NEWS entry

The entry is in `NEWS.md` at head `fb80dd8f`, under "Behaviour changes" (`:154`), at lines `:404-412`:

> `sg_tables()` now builds the effect column of its GLM tables with the fit's own orientation: the `adverse_outcome` that `forestsearch()` stored (for a `grf.subg.harm.glm()` result, the value stored on it), or `forestsearch()`'s rule when none is stored. For a continuous fit at the default orientation the subgroup effect reads MD 1.094, the fit's `max_sg_est`, where it read -1.094; a binary OR, RR or RD fit with `adverse_outcome = FALSE` now shows the effect for the non-event, as the fit does. `Mean(C)`, `Mean(T)` and `Diff` stay raw. Survival tables, count and rate-measure (IRR, IRD) rows, and fits with `adverse_outcome = TRUE` are unchanged.

| Claim | Check | Result |
|---|---|---|
| Effect column built with the stored value, or the rule when none is stored | Head `:856-859`, `:1042-1045`; tests `:220-234` | True |
| For a GRF result, the value stored on it | Head `:856`; `R/grf_subg_harm_glm.R:551/574/647/710` | True |
| "For a continuous fit at the default orientation the subgroup effect reads MD 1.094, the fit's `max_sg_est`, where it read -1.094" | `example.txt`, Questionable -1.094 → 1.094, `max_sg_est = 1.09411` | **Inexact.** The numbers come from one dataset (the test fixture), but the sentence reads as general and has no "for example". `max_sg_est` is the search's largest candidate (`max(hr.out$HR)`, `R/subgroup_search.R:1264`). It equals the selected subgroup's estimate in this example (test `:176`), not by definition. A GRF result has no `max_sg_est`. |
| Binary OR, RR or RD with FALSE shows the effect for the non-event, as the fit does | Start `R/glm_effect_estimators.R:390-391` (`Y -> 1 - Y`); `changed.txt`; `binary_false_example.txt` | True. Not said: GRF's default is FALSE for every type (`R/grf_subg_harm_glm.R:304`), so every binary GRF result at its default changes. |
| `Mean(C)`, `Mean(T)` and `Diff` stay raw | Tests `:143-148`; `example.txt` | True for continuous and count. **Inexact for binary**, which the sentence before is about: the binary arm columns are `Rate(C)` and `Rate(T)` (head `R/summary_utility_functions.R:345`; `binary_false_example.txt`). |
| Survival tables unchanged | `is_glm` guard; M_surv and V_surv_gbsg identical | True |
| Count and rate-measure (IRR, IRD) "rows" unchanged | Poisson route (start `R/glm_effect_estimators.R:162-172`, `:198` on); M_cnt_IRR_FALSE, G_cnt_FALSE identical; `changed.txt` | True. "Rows" is loose: what is unchanged is the whole table (its effect column) for count and IRR/IRD fits, not some rows. Count IRR with FALSE was run; binary IRR/IRD and count IRD with FALSE were not, and are unchanged by code reading. |
| Fits with `adverse_outcome = TRUE` unchanged | A1.5, 13/13 and 15/15 identical | True |

### A1.14 Evidence notes

- **The task's wording.** NOTES §8 ("The rule names `<fit>$args_call_all$adverse_outcome`") and VERIFY §1 ("departs from the task's wording (`args_call_all`)") do not match the task.
  - The task text (`b06de4e8:dev/tasks/TASK_cran_prep_displays_se_docs_2026-11-12.md:80-82`) says "the fit's stored `adverse_outcome`" and names no slot.
  - `args_call_all` comes from F12's code (start `R/plot_sg_glm_outcomes.R:169`).
- **The footnote line.** NOTES §2 puts the GRF CATE footnote at fixed `:926`; in the verified copy the string is at `:927` (verifier problem 3). Not corrected.
- **The polish counts.**
  - NOTES "Polish before commit" and `polish_result.json` give the polish delta as 2 files, +25/-4, and the regenerated patch as +289/-0.
  - After the main-session correction, `S/evidence/A1/polish.diff` is 2 files, +25/-3 (`git apply --stat`), and `S/patches/A1.patch` is +290/-0 (R +22), equal to `f4c7dd2d`.
  - The correction section gives no new counts.
  - NOTES polish item 3's "New (polished `:829-832`)" text is not what was committed; the committed text is at `:828-833` at the commit and `:850-855` at head.
- **Stale line numbers in NOTES.**
  - NOTES' header ("2 files, +268/-0") and §2's fixed-copy lines describe `A1.verified.patch`: GRF `:828-835`, call `:838-849`, argument `:845`; main `:991-998`, call `:1001-1012`, argument `:1008`. At the commit these are GRF `:828-837`, `:840-851`, `:847` and main `:993-1000`, `:1003-1014`, `:1010`.
  - NOTES §3–§4 give test-file lines for the pre-polish file: first failure `:182:3`; `.sgo_md()` `:97-104`; `:124-131`, `:136-141`, `:159`, `:185-186`, `:200`. At head these are `:199:3`; `:102-111`; `:126-151`; `:176`; `:202-203`; `:217`. The polish section records the shift only for the first failure.
- **Main-session red.** `S/main/A1/red.txt` lists 10 of the 41 failures and prints no totals. The 155/41 counts come from its progress line, and they agree with the implementer, verifier and polish red runs.

## A2 — CV summary tables use the fit's own orientation

`S` = `/tmp/cran_prep_displays_se_docs_2026-11-12`. Line numbers are at the start SHA `b06de4e8` and at head `fb80dd8f`. `R/` is identical at `835fac1b` and `fb80dd8f`. The test file does not change between `a2cd36db` and `fb80dd8f` (`git diff` empty), so its head line numbers equal those at the commit. `R/forestsearch_cross_validation.R` does change after the commit, in roxygen lines only (C4 `ba7c4ea2`, C5 `5073fd67`); the A2 code lines are unchanged and sit 14 lines lower at head.

### A2.1 Status, commit, provenance

| | |
|---|---|
| Status | **PASS** |
| Commit | `a2cd36db` "fix: CV summary tables use the fit's own orientation" |
| Parent | `f4c7dd2d` (A1) |
| Files (`git show --stat a2cd36db`) | `R/forestsearch_cross_validation.R` +12/-0; `tests/testthat/test-cv-tables-orientation.R` +234 (new). 2 files, +246/-0. No `man/`, `NEWS.md` or `DESCRIPTION` change. |
| Commit vs gated patch | `git diff f4c7dd2d a2cd36db` equals `S/patches/A2.patch` once the `index` lines are dropped (checked for this section). `S/main/A2/patch_equal.txt` is empty. |
| Patch history | `A2.r1.patch` (round 1; test +229; sha256 `46e5e27d...`) → `A2.verified.patch` (repair round 1, CONFIRMED in round 2; test +234; sha256 `495d45ae...`) → `A2.patch` (polish, committed; test +234; sha256 `62818222...`). The `R/` part is the same 12 lines in all three. |
| Later changes to the files | `R/forestsearch_cross_validation.R`: C4 `ba7c4ea2` adds 14 `#'` lines to `forestsearch_KfoldOut()`'s `@return` (head `:1326-1339`) and rewrites `CV_sgs()`'s `@return`; C5 `5073fd67` edits two `\dontrun` examples. All changed lines are roxygen. Test file: unchanged. |
| Origin | Part 4 report, F12 verifier minor problem (1): the CV summary tables kept the default `TRUE` (`b06de4e8:dev/reports/REPORT_cran_prep_code_fixes_2026-11-12.md`, F12 entry). Task: `b06de4e8:dev/tasks/TASK_cran_prep_displays_se_docs_2026-11-12.md:87` (A2) and the Part A rule (`:80-84`). |

### A2.2 Site, before → after

Function `forestsearch_KfoldOut()` (start `:1341`, head `:1355`). Both sites are inside `if (outall)` (start `:1493`, head `:1507`). The task's pointers (1509-1517, 1565-1575) are exact at the start SHA. These are the only two `make_effect_estimator(` calls in the CV code (verifier grep of `forestsearch_cross_validation.R`, `cv_summary_tables.R`, `summarize_cv_results.R`).

| piece | start `b06de4e8` | commit `a2cd36db` | head `fb80dd8f` |
|---|---|---|---|
| `is_glm <- ...` (unchanged) | `:1496` | `:1496` | `:1510` |
| new resolution block | none | `:1498-1506` | `:1512-1520` |
| site 1: estimator, no-original-subgroup ITT branch | `:1509-1517` | `:1519-1528`, new line `:1526` | `:1533-1542`, new line `:1540` |
| site 2: estimator, main GLM branch | `:1565-1575` | `:1576-1587`, new line `:1583` | `:1590-1601`, new line `:1597` |

Before, `b06de4e8:R/forestsearch_cross_validation.R:1565-1575` (site 2; site 1 `:1509-1517` has the same two argument lines):

```r
      estimator_fn_cv <- tryCatch(
        make_effect_estimator(
          outcome_type      = outcome_type,
          treat.name        = treat.name,
          outcome.name      = outcome.name,
          offset.name       = offset.name,
          effect_measure    = effect_measure,
          adjust_covariates = res$cv_args$adjust_covariates
        ),
        error = function(e) NULL
      )
```

After, `fb80dd8f:R/forestsearch_cross_validation.R:1512-1520` and `:1596-1598` (site 1: `:1539-1541`):

```r
    # Orientation for the GLM effect column: the fit's resolved
    # adverse_outcome (cv_args is the fit's args_call_all), so the table's
    # estimate is the fit's own, sign included.  NULL only on an object not
    # built by forestsearch(); resolved by forestsearch()'s rule.  Arm means
    # and Diff stay raw.
    adverse_outcome <- res$cv_args$adverse_outcome
    if (is_glm && is.null(adverse_outcome)) {
      adverse_outcome <- outcome_type %in% c("binary", "count")
    }
...
          effect_measure    = effect_measure,
          adverse_outcome   = adverse_outcome,
          adjust_covariates = res$cv_args$adjust_covariates
```

- Before the fix, neither call passed `adverse_outcome`, so both used `make_effect_estimator()`'s formal `TRUE` (`R/glm_effect_estimators.R:107` start, `:125` head).
- Where the stored value comes from: `forestsearch_Kfold()` sets `fs_args <- fs.est$args_call_all` (`:347`) and `cv_args <- fs_args` (`:389`); `forestsearch()` resolves `adverse_outcome` by its rule (`R/forestsearch_main.R:1840-1841` start, `:1877-1878` head) before `args_call_all <- mget(...)` (`:1924` start, `:1961` head). So every `forestsearch()` fit stores a logical. NULL arises only on an object not built by `forestsearch()`.
- Pattern: part 4's F12 (`R/plot_sg_glm_outcomes.R:166-172`, call `:268`, at `b06de4e8`). The one difference is the `is_glm` guard: the survival path never evaluates `outcome_type %in% ...` (a legacy object with `cv_args$outcome_type` removed has `is_glm` FALSE), and the survival tables (`SG_tab_estimates()`) do not take the value.
- Arm means and `Diff` come from `analyze_subgroup_glm()`, which the fix does not touch.
- Which tables can change: binary with OR, RD or RR (`.make_glm_binary_estimator()` flips `y -> 1 - y` when `adverse_outcome` is FALSE, `R/glm_effect_estimators.R:390-392` start, `:408-410` head) and continuous (`.make_lm_estimator()` negates y under FALSE). Count (IRR, IRD) and binary IRR/IRD go to `.make_poisson_rate_estimator()` (branches `:162-172` and `:198-212` start, `:180-190` and `:216-230` head), which has no `adverse_outcome` formal (`:737-746` start, `:826-835` head). Their tables cannot change.
- Reach: the package's own callers use `outall = FALSE` only (`forestsearch_Kfold()` `:621`, `forestsearch_tenfold()` `:1188`, same lines at head), so they never reach the change. A `forestsearch_tenfold()` result (class `fs_tenfold`, built at `:1284-1298`) has no `cv_args`, `sg0.name`, `sg1.name` or `resCV`, so `forestsearch_KfoldOut()` stops at its input check (head `:1365-1369`); the tables are built only from a `forestsearch_Kfold()` result.

### A2.3 Red and green

Filter `cv-tables-orientation`. It matches only the new file. Five `test_that()` blocks: (1) fixture premises, (2) all tables of the default-orientation fit, (3) no original subgroup, (4) no stored `adverse_outcome` (rule), (5) `adverse_outcome = TRUE`.

| run | tree | result | wall |
|---|---|---|---|
| main-session red | repository at `f4c7dd2d` plus only the patch's `tests/` hunks (`S/bin/apply_fix.sh`); `S/main/A2/red.txt` | 111 expectations: progress line `..........12...34..56....78..9a...b.......cd.......ef...gh..ij....kl..mn....op.................................` = 86 passes and 25 failures (markers 1-9, a-p). The summary reporter lists 10, then "... and 15 more". First: `test-cv-tables-orientation.R:183:3`, via `.cvo_check_row()` `:124`: "Expected MD cell ITT < `.cvo_tol`. Actual comparison: 0.457375922 >= 0.000000005" | 7.7 s |
| main-session green | same tree after the code hunks; `S/main/A2/green.txt` | 111 passes, `DONE`, no failures | 7.3 s |
| implementer red, round 1 | `wt/A2` copy of `start` + round-1 test file; `evidence/A2/red.txt` | 111 expectations, 25 failed, 0 errors, 86 passed; first failure at round-1 lines `:178:3` / `:119` (same message) | 8.34 s (test 8.17 s) |
| implementer green, round 1 | `wt/A2`; `evidence/A2/green.txt` | 111, 0 failed, 111 passed | 7.82 s (test 7.64 s) |
| implementer red, repair round | `wt/A2_r2_red`; `evidence/A2/r2/red.txt` | 111, 25 failed, 86 passed; first failure `:183:3` / `:124` | 8.53 s (test 8.36 s) |
| implementer green, repair round | `wt/A2_r2`; `evidence/A2/r2/green.txt` | 111, 0 failed | 7.90 s (test 7.75 s) |
| verifier red / green, round 1 | `verify/A2/red`, `verify/A2/fixed`; `evidence/A2/verify/{red,green}.txt` | 25 of 111 failed / 0 of 111 failed | 8.43 s / 7.66 s |
| verifier red / green, round 2 | `verify/A2_r2/{red,fixed}`; `evidence/A2/verify_r2/{red,green}.txt` | 25 of 111 failed / 0 of 111 failed | 7.92 s / 7.36 s |
| polish red / green | `wt/A2_polish_red`, `wt/A2_polish`; `evidence/A2/polish/{red,green}.txt` | 25 of 111 failed / 0 of 111 failed | 8.37 s / 7.79 s |
| polish controls (verified test file) | `wt/A2_polish_redctl`, `wt/A2_polish_vctl`; `polish/{red,green}_verified.txt` | 25 / 0 of 111, same per-block counts | 8.39 s / 7.96 s |
| polish check green | `verify/A2_polish`; `evidence/A2/verify_polish/green.txt` | 111, 0 failed (blocks 6/29/11/33/32) | 7.68 s (test 7.51 s) |
| full suite at head (in progress when written) | `S/suite/suite.log:344` | `✔ \|        111 \| cv-tables-orientation [3.1s]` | — |

- Red failures per block (every red run): fixture 0 of 6, main 11 of 29, no subgroup 2 of 11, rule fallback 12 of 33, `TRUE` 0 of 32.
- Every failure is a sign mismatch. Each row check has two failing expectations (cell, sign) and two passing ones (length, raw `Diff`). Main block: 5 rows x 2 + the printed tie `:186` = 11. Rule block: 5 x 2 + the no-subgroup ITT row = 12. No-subgroup block: 2.
- The cell differences are 2|est| or the mirrored limit: ITT 0.457, FA H 2.188, FA Hc 0.192, KfA H 1.603, KfA Hc 0.064, tie 2.188 (verifier r1 and r2).
- Main-session walls are `date` differences around each `Rscript` process (`apply_fix.sh`), rounded to 0.1 s. The others are `/usr/bin/time` or `proc.time()` in the runner.

### A2.4 Expected values

No expected number is a literal in the test. Each is computed at run time from the data by `.cvo_md()` (`:103-110`), a closed-form OLS on the subset, or is the fit's own `max_sg_est`. The values below are what those computations give on the fixture, as re-derived outside the package.

Fixture (`:52-58`, `:80-99`): `.make_continuous_data(N = 300L, MD_harm = 1.5, seed = 42L)`; `ben` is the same data with y negated. `fit_f` = `ben`, `adverse_outcome` unset (resolves to FALSE); `fit_t` = harm twin, `adverse_outcome = TRUE`; `fit_0` = `ben` with `effect.threshold = 10`, `pconsistency.threshold = 0.99` (no subgroup). Each is followed by `forestsearch_Kfold(Kfolds = 3L)`.

| expected value | how it was derived | test line (head) |
|---|---|---|
| MD = OLS treatment coefficient = `mean(y\|T=1) - mean(y\|T=0)`; SE `sqrt(s2 * (1/n0 + 1/n1))`, `s2 = (SS0 + SS1) / (n0 + n1 - 2)` | By hand (closed form of OLS on a 0/1 regressor) | comment `:21-23`; `.cvo_md()` `:103-110` |
| Orientation: under FALSE the fit's MD is `-(mean(y\|T=1) - mean(y\|T=0))`, same SE | Documentation: `?forestsearch` `@param adverse_outcome` (`R/forestsearch_main.R:1104-1107` start, `:1126-1129` head); package convention effect > 0 = harm | comment `:24-28`; `:108` |
| Limits `est -/+ 1.96 * se` | The literal z of `analyze_subgroup_glm()` (`R/summary_utility_functions.R:223`, same at head) | comment `:29-30`; `:122` |
| ITT 0.22868796 (-0.00700868, 0.46438460), se 0.12025339, n 300 | Above, all of `resCV` | `:139` → `:124-125` |
| FA H 1.09411222 (0.52314346, 1.66508098), se 0.29131059, n 36 | Above, `treat.recommend.original == 0`; row order sg0 then sg1 from `SG_tab_estimates_glm()`'s documented return | `:144`, comment `:140-141` |
| FA Hc 0.09579320 (-0.14638670, 0.33797310), se 0.12356117, n 264 | Above, `treat.recommend.original == 1` | `:145` |
| KfA H 0.80136279 (0.36076387, 1.24196171), se 0.22479537, n 71 | Above, `treat.recommend == 0` (held-out assignments) | `:150`, premise both levels `:147` |
| KfA Hc 0.03182661 (-0.23692130, 0.30057453), se 0.13711628, n 229 | Above, `treat.recommend == 1` | `:151` |
| `Diff` raw: -0.22868796, -1.09411222, -0.09579320, -0.80136279, -0.03182661 | `mean(y\|T=1) - mean(y\|T=0)`, not oriented (task: arm means and their raw difference are not touched) | `:126-128` |
| `tab_all` carries the same cells, rows `Overall, FA_0, KfA_0, FA_1, KfA_1` | Structure of the start code (unchanged) | `:154-157` |
| Hand H = `max_sg_est` (1.09411222056755), tolerance 1e-10, FALSE fit and TRUE twin | `max_sg_est = max(hr.out$HR)` (`R/subgroup_search.R:1264`) is the closed form of the same coefficient on the same subjects; `fit_glm_for_subgroup_fast()` negates y first under FALSE (`:1132`) | `:176`, `:230`; reason `:36-42` |
| Printed H cell = `max_sg_est` within `.cvo_tol` | Same tie, at the printed digit | `:186` |
| Planted harm reads as harm: hand H > 0 | Sign convention above | `:177` |
| No subgroup: ITT row of `fit_0` = the ITT values above | Same formula on the whole frame; premise `abs(raw) > 1e-3` (raw -0.22868796) so the two signs are far apart | `:204-205`, reason `:201-202` |
| Rule fallback: `cv_args$adverse_outcome` removed gives the same values (continuous → FALSE) | `forestsearch()`'s rule `outcome_type %in% c("binary", "count")` | `:213-215`, `:217-219` |
| `TRUE` twin: raw MD, same magnitudes and signs | Orientation TRUE = no negation | `:227-231`; mirror premise `:233` |

- Tolerances, each set before the first run with its reason in a comment:
  - `.cvo_tol = 5e-9 + 1e-12` (`:46`): tables built at `digits = 8`, so half a unit of the last digit plus the parse (`:34-36`). Largest fixed-table deviation observed: 4.922e-09 (`verify_r2/hand_r2.txt`).
  - 1e-10 for the `max_sg_est` tie (`:176`, `:230`; reason `:36-42`). Observed: 2.22e-16 (both fits).
  - 1e-3 premise (`:204`).
- Independent re-derivations, all equal to every printed digit:
  - Implementer `evidence/A2/hand_values.R`: base R, data regenerated, H by the printed definition `{age <= 47} & {biomarker_hi}` (ITT, H, Hc).
  - Verifier round 1 `verify/hand_check.R`: `stats::t.test(var.equal = TRUE)`, se = |diff|/|t| (wall 5.14 s for the fixed run).
  - Verifier round 2 `verify_r2/hand_r2.R`: normal equations `b = (X'X)^-1 X'y_o`, package not loaded (wall 0.19 s). It also confirms H (n 36) equals both the `resCV` set and the fit's own H, and that the subset ids are `identical()` start vs fixed.
- The verifiers found no value copied from the fixed code's output.

### A2.5 Unaffected case

- **Task:** survival fits, and non-survival fits with `adverse_outcome = TRUE`, give identical display data.
- **Method:** each tree in its own `Rscript --vanilla` process, results saved with `saveRDS()`, compared with `identical()` in a third process. Each case runs `forestsearch()` → `forestsearch_Kfold(Kfolds = 3)` → `forestsearch_KfoldOut(outall = TRUE)`, `(outall = TRUE, digits = 8)` and `(outall = FALSE)`, and saves the Kfold object. The only leaf dropped is `timing_minutes` (wall clock).

| script (author) | cases | result | walls start / fixed (compare) |
|---|---|---|---|
| `evidence/A2/unaffected.R` (implementer) | 11: `surv_sg`, `surv_null` (survival, stored FALSE by the rule, unused); `cont_true_sg`, `cont_true_null`, `cont_true_adj` (continuous, explicit TRUE; adjusted by `biomarker`); `bin_or_rule_sg` (OR, rule TRUE); `bin_rd_true_sg` (RD, explicit TRUE); `bin_null`; `bin_or_stored_NULL` (`cv_args$adverse_outcome` removed → rule TRUE); `count_rule_sg`, `count_null` | `identical(start, fixed)` TRUE per case and overall; start control also identical (`unaffected_compare.txt`) | round 1: 18.10 s / 17.42 s, control 17.81 s (0.18 s); repair round: 17.03 s / 17.31 s |
| `verify/unaffected_v.R` (verifier r1) | 20: survival rule, explicit TRUE, legacy (no `outcome_type`), null; continuous TRUE x 4 (other seed, no subgroup, no `effect_measure`); binary x 8 (OR, RR, RD; rule and explicit TRUE; adjusted; RD fallback); count x 4 incl. `n2_irr_FALSE`, `n3_irr_FALSE_off` (explicit FALSE; offset merged into `resCV`) | all TRUE (`verify/unaff_compare.txt`) | 31.85 s / 31.96 s (0.16 s); repair re-run 31.55 s / 31.17 s (0.19 s) |
| `verify_r2/unaffected_r2v.R` (verifier r2) | 13 new: survival explicit FALSE, Kfolds 5, tenfold; continuous TRUE adjusted no subgroup, Kfolds 5, stored TRUE by edit, tenfold; binary OR explicit TRUE, RD rule null, legacy (both trees stop with the same Lapack "exactly singular" error), OR tenfold; count IRD TRUE and IRR rule null with offset. Also saves `details = TRUE` printout and `cv_summary_tables(use_gt = FALSE)` | all TRUE (`verify_r2/unaff_r2v_compare.txt`) | 25.11 s / 25.46 s (0.16 s) |
| re-runs by verifier r1 / r2 | implementer's 11; verifier r1's 20 | all TRUE; start RDS `identical()` across rounds (`verify_r2/unaff_compare_all.txt`) | r1: 17.69 s / 18.05 s; r2: 17.88 s / 17.64 s and 31.94 s / 30.81 s |

- 44 distinct configurations, all identical start vs fixed.
- Some go beyond the task's case: `n2_irr_FALSE`, `n3_irr_FALSE_off` (count at FALSE) show that count tables do not move even at FALSE with the IRR populated (`IRR 1.3797 (1.1084, 1.7173)` on both trees, `r2/compare.txt`).
- Positive control (`verify_r2/unaff_compare_all.txt`): the same compare on `extract_*.rds` gives FALSE for F, F_rule, Z, Z_rule and default digits, TRUE for T. The comparison detects the change.
- The two implementer count cases print `IRR` "NA" in every row on both trees (open item 1), so they only show that nothing else in the count tables moved.

### A2.6 Existing tests

One process per file, `NOT_CRAN=true timeout 300`, fixed copy and start control, repair round (`evidence/A2/r2/existing/`). Round 1 and both verifier rounds give the same counts.

| file | expectations | result (both trees) | wall fixed / start |
|---|---|---|---|
| cv-exact-match | 6 | 6 passed | 3.18 / 3.03 s |
| cross-outcome-parity | 49 | 35 passed, 2 skipped (primary fit did not identify), 12 warnings, 0 failed | 46.91 / 46.37 s |
| cv-no-subgroup-dina-grf | 48 | 48 passed | 14.46 / 14.57 s |
| cv-fold-na-assignments | 6 | 6 passed | 3.34 / 3.38 s |
| cv-no-subgroup-edges | 45 | 44 passed, 1 skipped (multisession dev-load) | 34.36 / 34.78 s |
| cv-mixed-and-fatal | 58 | 58 passed | 18.48 / 18.34 s |
| downstream-consumers | 43 | 38 passed, 5 warnings | 22.95 / 22.84 s |
| input-validation | 19 | 18 passed, 1 warning | 4.63 / 4.63 s |
| rng-state-restored | 16 | 16 passed | 7.20 / 7.38 s |
| missing-membership-messages | 60 | 60 passed | 8.11 / 8.03 s |
| no-subgroup-cv (the only existing `outall = TRUE` file) | 72 | 72 passed | 30.23 / 29.84 s |
| display-orientation (F12) | 93 | 93 passed | 7.42 / 7.53 s |
| mr-inference (whole file; `outall = FALSE` only) | 148 | 148 passed, 0 failed, 0 errors (fixed only) | 181.61 s (round 1 188.14 s; verifier r1 189.39 s, r2 187.38 s) |

- No failure or error. Skips and warnings are the same on both trees.
- Expectations updated: none.

### A2.7 What a user sees

Continuous fit at the default orientation (`adverse_outcome` unset → FALSE), the test's `ben` data, `sg.harm` `{age <= 47} & {biomarker_hi}`, `max_sg_est` 1.09411. `forestsearch_KfoldOut(forestsearch_Kfold(fit, Kfolds = 3), outall = TRUE)$tab_all`, default `digits = 4` (`evidence/A2/example.R`, `example.txt`):

```
before (start tree)                          after (fixed tree)
                                MD                                         MD
Overall  -0.2287 (-0.4644, 0.0070)           Overall 0.2287 (-0.0070, 0.4644)
FA_0    -1.0941 (-1.6651, -0.5231)           FA_0     1.0941 (0.5231, 1.6651)
KfA_0   -0.8014 (-1.2420, -0.3608)           KfA_0    0.8014 (0.3608, 1.2420)
FA_1     -0.0958 (-0.3380, 0.1464)           FA_1    0.0958 (-0.1464, 0.3380)
KfA_1    -0.0318 (-0.3006, 0.2369)           KfA_1   0.0318 (-0.2369, 0.3006)
```

- `FA_0`, the selected subgroup, now reads as the fit's `max_sg_est`, 1.0941.
- `Mean(C)`, `Mean(T)` and `Diff` are the same on both trees (`Diff` `FA_0` -1.0941), so `Diff = -MD` for this fit, as in F12's displays.
- Reproduced by verifier r1 (`verify/hand_check_{start,fixed}.txt`) and verifier r2 (`verify_r2/extract.R`, default digits).
- Binary neighbours at `adverse_outcome = FALSE` (verifier probes, not in the test):
  - OR, H cell: 0.60000000 → 1.66666667, the hand OR on `1 - y` (`verify/affected_probe.txt`; walls 6.64 s / 6.43 s).
  - RR, ITT: 1.00751880 → 0.99319728; RD, ITT: 0.00357143 → -0.00357143, each the hand value on `1 - y` within 1.4e-9 (`verify_r2/affected_r2_probe.txt`; walls 8.30 s / 8.26 s).
  - Binary IRR with an offset: `1.04760702 (0.81325660, 1.34948855)` on both trees.

### A2.8 Verification

- **Round 1 (verifier, `evidence/A2/verify/VERIFY.md`): PROBLEMS.** The code and test passed every gate. Three blocking problems, all text about count fits:
  1. NOTES open item 1 said the fix "applies to count once the offset is carried". False: `make_effect_estimator()`'s count branch never passes `adverse_outcome` to `.make_poisson_rate_estimator()`, which has no such formal. Case `n3_irr_FALSE_off` (offset merged, FALSE) is identical on both trees.
  2. NOTES section 8 said an explicit FALSE "for binary or count is changed". False for count (`n2`, `n3` identical).
  3. The proposed NEWS bullet named count fits as changed.
  - Minor: (1) test comment "both sides are the same arithmetic on the same numbers" not literal; (2) NOTES open item 3 missed that the binary `"RD"` fallback is also in the main GLM branch; (3) NEWS "matches the fit's estimate": for an adjusted continuous fit the H cell is the adjusted MD on H (0.91517984), not `max_sg_est` (0.97140465) (`verify/adj_probe.txt`).
- **Repair round 1 (implementer):** all three blocking problems accepted; no code change.
  - NOTES sections 1, 8, 9 (items 1 and 3) and 10 (NEWS) corrected in place.
  - Test header: the tie-tolerance reason, 2 comment lines → 7 (+5; `A2.r1.patch` → `A2.verified.patch`). No expectation, seed or tolerance changed.
  - New count probe (`r2/count_probe.txt`; walls 4.99 s / 4.94 s): `make_effect_estimator()` at TRUE and FALSE gives identical estimate and SE for count IRR (0.75837683, se 0.18291257), count IRD, binary IRR and binary IRD; binary OR flips (0.81093022 vs -0.81093022). A count `forestsearch()` fit at TRUE and at FALSE gives the same `sg.harm` and `max_sg_est` (0.81655414).
  - Every gate re-run (A2.3, A2.5, A2.6).
- **Round 2 (verifier, `evidence/A2/verify_r2/VERIFY.md`): CONFIRMED.** No blocking problem. Two minor:
  1. Test `:38`: "an lm.fit of y on treatment within the candidate" is imprecise for the FALSE fixture; the search negates y first (`R/subgroup_search.R:1132`). **Fixed by the polish.**
  2. Stale round-1 test-file pointers in NOTES (`:34-37` for the tolerance comment, which is `:34-42`; `:171`, `:181`, `:225`). **Not fixed**; this section uses head numbers (see A2.12).
- **Polish (main session, comment only):** test `:36-42` rewrapped, 7 comment lines out and 7 in, file stays 234 lines; the parenthetical now says "with y negated first under FALSE by fit_glm_for_subgroup_fast()".
  - Basis: `R/subgroup_search.R:1132` at the start SHA; probe `polish/fastpath_probe.txt` (106 fast-path calls, 0 slow, closure `adverse_outcome` FALSE, `max_sg_est` 1.09411222056755; wall 3.12 s).
  - Comment gate (`git diff -U0 verified -- R tests` minus comment lines) printed nothing; the `R/` part is byte-identical to the verified patch's.
- **Polish check (`polish_result.json`, A2): CONFIRMED, no problems.** Its probe (`verify_polish/probe_negation.txt`, wall 4.32 s): FALSE fixture, every fast-path estimate equals the coefficient of -y within 1.33e-15 (up to 2.19 from the coefficient of y); TRUE fixture the reverse; `max_sg_est` `identical()` to `max(hr.subgroups$HR)`.

### A2.9 Departures and interpretations

- **One resolution block for both sites**, placed before the `no_original_sg` branch, instead of one copy per site. Same values as F12's per-display block.
- **`is_glm` guard** on the fallback. The survival path never reads the value, and a legacy object without `cv_args$outcome_type` is not touched (verifier case `s3_legacy`, identical).
- **The test goes beyond the task's minimum:** every GLM table `KfoldOut` builds (ITT, full-data, K-fold, `tab_all`), the no-original-subgroup branch (reached deterministically with `effect.threshold = 10`, premise `expect_null(fit$sg.harm)` `:195`), the rule fallback with `cv_args$adverse_outcome` removed, and the `TRUE` twin as a guard. The fixture and `TRUE` blocks pass on both trees by design.
- **Not tested here:** a binary fit with explicit `adverse_outcome = FALSE` and OR, RD or RR changes under this fix; the task's test names only the continuous default. Covered by the verifier probes in A2.7.
- **Unaffected case extended** beyond the task's (count at FALSE, tenfold, Kfolds 5, `cv_summary_tables()`); all identical.
- **Task text vs code:** the pointers are exact at the start SHA; no other difference recorded.
- **Extra check, not asked by the task:** NOTES section 8 says the task asked whether `cv_args` holds `adverse_outcome` for Kfold and tenfold results. The task's A2 item (`:87`) does not ask this (see A2.12). The check itself is right: both build `cv_args` from `fs.est$args_call_all` (`:347`/`:389`, `:874`/`:892`).

### A2.10 Open items

1. **Count CV tables print `IRR` "NA"** in every row, before and after: `forestsearch_Kfold()`'s `keep_cols` (`R/forestsearch_cross_validation.R:526-529`, same at head) does not carry the offset column into `resCV`, so the count estimator fails inside `tryCatch`. Not part of A2. Carrying the offset would fill the cells, but A2 would still not orient them.
2. **Count and binary IRR/IRD are never oriented**, by the CV tables or by `forestsearch()` itself (`R/forestsearch_main.R:2263-2272` at the start SHA builds the fit's own count estimator through the same branch). Whether they should honour `adverse_outcome` is a question for the estimator, not for the CV tables, which follow the fit.
3. **Adjusted continuous fit:** the table's H cell is the adjusted MD on H (0.91517984), while `max_sg_est` is the search maximum (0.97140465) (`verify/adj_probe.txt`). Predates A2; A2 aligns the orientation only.
4. **Binary `effect_measure` fallback "RD"** in both GLM branches of `forestsearch_KfoldOut()` (start `:1506-1507`, `:1555-1560`; head `:1530-1531`, `:1580-1585`), while `forestsearch()` uses "OR". Unreachable for `forestsearch()` fits (the resolved value is stored). Not changed.
5. **quarto/ (read-only, not re-rendered):** by literal grep, the non-archive callers of `forestsearch_KfoldOut()` are survival, binary at TRUE, or count (in `_broken/`); continuous CV analyses are only under `quarto/applications/actg175/_archive/`. Files that set `outcome_type` through variables were not traced (the five actg175 files that do set `adverse_outcome = TRUE`, verifier r2).
6. **Pre-existing error, noted:** a binary result with `cv_args$outcome_type` removed stops in `KfoldOut(outall = TRUE)` with "Lapack routine dgesv: system is exactly singular" on both trees (verifier r2 case `b11_rule_legacy`).

### A2.11 NEWS entry (`fb80dd8f:NEWS.md:414-423`)

Text: "`forestsearch_KfoldOut(outall = TRUE)`, the tables of a `forestsearch_Kfold()` or `forestsearch_tenfold()` result, now builds the effect column of its GLM tables (ITT, full-data and fold subgroups) with the fit's own `adverse_outcome`, or `forestsearch()`'s rule when none is stored. A continuous fit at the default orientation showed the selected subgroup as MD -1.0941 where the fit's estimate is 1.0941; it now shows 1.0941. Binary OR, RR and RD fits with `adverse_outcome = FALSE` now show the effect for the non-event. Arm means and `Diff` stay raw. Survival tables, count and rate-measure (IRR, IRD) rows, and fits with `adverse_outcome = TRUE` are unchanged."

| claim | check | verdict |
|---|---|---|
| "the tables of a ... `forestsearch_tenfold()` result" | A `forestsearch_tenfold()` result (`fs_tenfold`, `:1284-1298`) has no `cv_args`, `sg0.name`, `sg1.name` or `resCV`; `forestsearch_KfoldOut()` stops at its input check (head `:1365-1369`). `forestsearch_tenfold()` calls it with `outall = FALSE` only (`:1188`). NOTES sections 2 and 8: "tenfold never builds them". | **Inexact.** Only a `forestsearch_Kfold()` result gets these tables. |
| effect column of the GLM tables (ITT, full-data and fold subgroups) with the fit's `adverse_outcome`, rule when none stored | Head `:1512-1520`, `:1540`, `:1597`; test blocks 2-4 | Correct |
| continuous default: selected subgroup MD -1.0941 → 1.0941, the fit's estimate | `example.txt` (`FA_0`), `max_sg_est` 1.09411 | Correct for that fit. Minor: the numbers are one synthetic fit's (the test fixture) and the bullet does not say "for example". For an adjusted continuous fit the H cell is not `max_sg_est` (open item 3), so "the fit's estimate" holds in general for the sign only. |
| binary OR, RR, RD at FALSE now show the effect for the non-event | `y -> 1 - y` under FALSE (`R/glm_effect_estimators.R:408-410` head); probes A2.7 | Correct |
| arm means and `Diff` stay raw | Fix does not touch `analyze_subgroup_glm()`; test `:126-128`; example | Correct (for binary the arm columns are `Rate(C)`/`Rate(T)`) |
| survival, count and IRR/IRD rows, `TRUE` fits unchanged | 44 unaffected configurations; `n2`, `n3`; count probe; binary IRR probe | Correct |

### A2.12 Corrections to the notes

- **Stale test-file pointers (NOTES sections 3, 4, 6; verifier r1).** They are round-1 numbers; the committed file has 5 more comment lines after line 37. At head: derivation comment `:19-33`, tolerance comment `:34-42` (NOTES `:34-37`), `.cvo_md()` `:101-110` (NOTES `:96-104`), ties `:176` and `:230` (NOTES `:171`, `:225`), printed tie `:186` (NOTES `:181`), premise `:204` (NOTES `:199`), rule fallback `:215`, `:219` (NOTES `:210`, `:214`), first red failure `:183` / `:124` (NOTES `:178` / `:119`). NOTES gives a "+5 after `:37`" rule, which cannot map `:34-37` (verifier r2 minor 2). Head numbers are used here.
- **Binary IRR/IRD branch range.** NOTES (sections 8, 9) and verifier r1 cite `R/glm_effect_estimators.R:160-169` at `b06de4e8`. Line 160 is the `choices` line of `match.arg()`; the branch is `:162-172` (`if` `:162`, call `:163-172`), as verifier r2 gives it.
- **Round-1 start wall 18.10 s (NOTES section 5).** No log backs it: `evidence/A2/` has `unaffected_start.rds` but no `unaffected_start.txt`. The fixed (17.42 s) and control (17.81 s) walls are in their `.txt` files. Not contradicted; the repair-round start wall (17.03 s) is backed by `r2/unaffected_start.txt`.
- **"red equals `verify_r2/red.txt`" (NOTES polish section).** The two files differ in their two wall lines (8.21 s / 8.37 s vs 7.74 s / 7.92 s); all else is equal.
- **Tenfold cases (verifier r2 Check 4a).** "Each case saves ... `KfoldOut(outall = TRUE)` ... `cv_summary_tables()`" does not hold for `s7_tenfold`, `c8_true_tenfold`, `b12_or_rule_tenfold`: the script's `ten()` saves only the `forestsearch_tenfold()` object.
- **"task asked to check" (NOTES section 8).** The task document's A2 item (`:87`) names only the two sites; no request about `cv_args` in Kfold and tenfold results is in it. The phrase "tenfold results" there means the per-simulation `res` built inside `forestsearch_tenfold()` (`:1169-1181`), which holds `cv_args`; the returned `fs_tenfold` object does not.
- **quarto caller count (verifier r1 Check 5).** "38 callers are outside `_archive/` and `_broken/`" is not reproduced by `git grep -l forestsearch_KfoldOut b06de4e8 -- quarto`: 70 files outside both (37 non-HTML, 35 `.qmd`). Verifier r2's "82 files outside `_archive/`" matches git. The conclusion (literal grep finds no continuous or FALSE call outside `_archive/` and `_broken/`) does not depend on the count.
- **Head re-run.** `S/cumulative/A2/` holds a head-tree re-run of the three unaffected scripts (written 17:12-17:15); it had no compare output when this section was written, so it is not reported here.

## A3 — `compute_frontier_cis()` falls back to `forestsearch()`'s orientation rule

`S` = `/tmp/cran_prep_displays_se_docs_2026-11-12`. Line numbers are at the start SHA `b06de4e8` and at head `fb80dd8f`. `R/` is identical at `835fac1b` and `fb80dd8f`. The two files of this fix do not change between `fb21396f` and `fb80dd8f` (`git diff fb21396f fb80dd8f -- R/frontier_cis.R tests/testthat/test-frontier-cis-orientation.R` is empty), so head line numbers equal those at the commit. Test-file lines before the polish ("verified" numbering) map to head as follows: verified `:n` is head `:n+1` for `12 <= n <= 186`, and verified `:187-189` is head `:191-193`.

### A3.1 Status, commit, provenance

| | |
|---|---|
| Status | **PASS** |
| Commit | `fb21396f` "fix: compute_frontier_cis() falls back to forestsearch()'s orientation rule" |
| Parent | `a2cd36db` (A2) |
| Files (`git show --numstat fb21396f`) | `R/frontier_cis.R` +9/-1; `tests/testthat/test-frontier-cis-orientation.R` +194 (new). 2 files, +203/-1 (`S/main/A3/diffstat.txt`). No roxygen, `man/`, `NEWS.md` or `DESCRIPTION` change. |
| Commit vs gated patch | `git diff a2cd36db fb21396f` equals `S/patches/A3.patch` once the `index` lines are dropped (checked for this section). `S/main/A3/patch_equal.txt` is empty. |
| Earlier patch | `S/patches/A3.verified.patch` (224 lines, before the polish): same `R/` hunk; test file 190 lines (+199/-1 in all). Final `A3.patch`: 228 lines. |
| Added lines | 0 non-ASCII, none over 80 characters (`git show fb21396f`, checked for this section). The em dashes in `R/frontier_cis.R` roxygen (`:31`, `:37`, `:52`) predate the fix. |
| Task text | `b06de4e8:dev/tasks/TASK_cran_prep_displays_se_docs_2026-11-12.md`, Part A, A3: "replace the `%||% TRUE` fallback with the rule above. A stored value is used exactly as now." |
| Workflow | `S/workflow1_result.json`: `status` READY, `final_verdict` CONFIRMED, `repaired` false, `blocking_r1` []. `S/polish_result.json`: polish DONE, check CONFIRMED, `problems` []. |

### A3.2 Site, before → after

Function `compute_frontier_cis()`, step "1. Build the per-data-slice estimator closure". `outcome_type` is resolved at `:136` (`fs$outcome_type %||% "survival"`, both trees). The `%||%` bound in the namespace is `function(a, b) if (is.null(a)) b else a` (verifier, `get('%||%', asNamespace('forestsearch'))` under `load_all()`), so the old expression was `if (is.null(x)) TRUE else x`.

Before, `b06de4e8:R/frontier_cis.R:150` (the `make_effect_estimator()` call is `:143-153`, inside `tryCatch()` `:142-158`):

```r
      adverse_outcome   = args$adverse_outcome %||% TRUE,
```

After, `fb80dd8f:R/frontier_cis.R:140-147` (new block after `effect_measure` is resolved) and `:158` (the call is `:151-161`, inside `tryCatch()` `:150-166`):

```r
  # Orientation: the fit's stored adverse_outcome.  A fit returned by
  # forestsearch() always stores it (resolved before args_call_all is
  # captured); if none is stored, use forestsearch()'s own rule on the
  # resolved outcome_type, not TRUE for every outcome type.
  adverse_outcome <- args$adverse_outcome
  if (is.null(adverse_outcome)) {
    adverse_outcome <- outcome_type %in% c("binary", "count")
  }
...
      adverse_outcome   = adverse_outcome,
```

- The rule is `forestsearch()`'s: `b06de4e8:R/forestsearch_main.R:1840-1842` (Section 1A, header `:1766`), again at `:2241-2243` under the "Resolve adverse_outcome default" comment (`:2236-2240`). F12 used the same form (`b06de4e8:R/plot_sg_glm_outcomes.R:169-172`).
- A stored value is read by the same `args$adverse_outcome` as before. Any non-NULL value (TRUE, FALSE, NA, 0, 1, `"TRUE"`, `logical(0)`) reaches `make_effect_estimator()` unchanged; only a NULL entry changes path.

**When the fallback is reached.** `forestsearch()` resolves `adverse_outcome` (`b06de4e8:R/forestsearch_main.R:1840-1842`) before it captures its formals (`args_call_all <- mget(names(formals()), ...)`, `:1923-1924`). `.sync_args_call_all()` (`:21-30`) writes with `[n] <- list(...)` and never deletes. So every fit `forestsearch()` returns stores a logical. Only an object without the entry reaches the fallback. Checked for this section with `git log -S`: every `forestsearch()` since `adverse_outcome` was added (`99ced431`, 2026-03-24) writes the resolved value back (`args_call_all$adverse_outcome <- adverse_outcome`, e.g. `99ced431:R/forestsearch_main.R:372`, `bc4690b3^:R/forestsearch_main.R:693`). Older objects have no `pareto_frontier` (first in `bdef2595`, 2026-04-27), and `compute_frontier_cis()` returns an empty table before it builds the estimator (`b06de4e8:R/frontier_cis.R:113-118`). So in practice only an object built or edited by hand reaches the changed branch (git and code reading; not run).

### A3.3 What changes, by outcome type (no value stored)

`make_effect_estimator()` reads `adverse_outcome` only on the binary non-rate path (`.make_glm_binary_estimator`, `b06de4e8:R/glm_effect_estimators.R:174-183`) and the continuous path (`.make_lm_estimator`, `:187-195`, which negates `y` when the value is FALSE, `:956-958`). The Cox and Poisson estimators do not receive it.

| outcome type | old fallback | rule | value read? | output |
|---|---|---|---|---|
| binary | TRUE | TRUE | yes (OR, RR, RD) | unchanged |
| count | TRUE | TRUE | no | unchanged |
| survival | TRUE | FALSE | no | unchanged |
| continuous | TRUE | **FALSE** | yes | **changes: every estimate and limit negated, limits swapped** |

### A3.4 Red and green

Command (all runs): `NOT_CRAN=true timeout 300 Rscript --vanilla -e "devtools::test('<tree>', filter = 'frontier-cis-orientation', reporter = 'summary', stop_on_failure = FALSE)"`. Per-test counts from a second run with the silent reporter (`S/evidence/A3/counts.R`, verifier `verify/counts.R`). Red tree = start (or parent) + the test file only.

| run | tree | expectations | failed (premise / selected / every row / identical) | errors | skipped | wall | file |
|---|---|---|---|---|---|---|---|
| implementer red | `S/wt/A3` before the code (start + test) | 31 | 18 (0/5/12/1) | 0 | 0 | 3.96 s | `evidence/A3/red.txt`, `red_counts.txt` |
| implementer green | `S/wt/A3` | 31 | 0 | 0 | 0 | 3.91 s | `green.txt`, `green_counts.txt` |
| verifier red | `S/verify/A3/red` | 31 | 18 (0/5/12/1) | 0 | 0 | 4.17 s | `verify/red_run.txt`, `red_counts.txt`, `red_wall.txt` |
| verifier green | `S/verify/A3/fixed` | 31 | 0 | 0 | 0 | 3.83 s | `verify/fixed_run.txt`, `fixed_counts.txt`, `fixed_wall.txt` |
| polish red | `S/wt/A3_polish_red` | 31 | 18 (0/5/12/1) | 0 | 0 | 4.15 s (summary), 3.95 s (counts) | `evidence/A3/polish/polish_red*.txt` |
| polish green | `S/wt/A3_polish` | 31 | 0 | 0 | 0 | 3.86 s (summary), 3.78 s (counts) | `polish/polish_green*.txt` |
| polish check red | start + polished test | 31 | 18 (0/5/12/1) | 0 | 0 | 4.17 s | `evidence/A3/verify_polish/red.txt` |
| polish check green | `S/verify/A3_polish` | 31 | 0 | 0 | 0 | 3.83 s | `verify_polish/green.txt` |
| **main red** | repository at `a2cd36db` + test file only | 31 | 18 not passed (progress line: 13 dots, 18 marks; the 10 listed are all "Failure", then "... and 8 more") | not counted | not counted | 3.8 s | `S/main/A3/red.txt`, `wall.txt` |
| **main green** | repository with the code (committed as `fb21396f`) | 31 | 0 (31 dots) | 0 | 0 | 3.5 s | `S/main/A3/green.txt`, `wall.txt` |

- First failure, main session: `test-frontier-cis-orientation.R:147:3`, "Expected `row$estimate` to equal `ref$est`." `actual`: -1.1, `expected`: 1.1. The implementer and first-round verifier report the same failure at `:146:3` (verified numbering, one line lower).
- All 18 failures are sign mismatches: `:147-151` (selected subgroup: estimate, sign, naive limits, Split-CI midpoint), `:173`, `:175`, `:177`, `:179` for each of the 3 frontier rows, and the relational `:191`.
- The 13 expectations that pass in red do not depend on the fix: 6 premises (`:116-125`), `:136` (selected `m` is in the table), `:144` (hand value = `max_sg_est`), `:145` (`ref$est > 0`), `:163` (row count) and 3 label-form checks (`:167`).
- The polish check's walls (4.17 s, 3.83 s) equal the first-round verifier's to 0.01 s. They are separate runs: the first failure is at `:147:3` in the polish check and at `:146:3` in the first round.

### A3.5 Expected values

Fixture: `.make_continuous_data(N = 300L, MD_harm = 1.5, seed = 42L)`, `y` negated, so the planted region is where treatment harms under the continuous default (`fb80dd8f:tests/testthat/test-frontier-cis-orientation.R:48-57`). One `forestsearch()` fit, continuous, MD, `adverse_outcome` unset (`:59-78`). `.fco_unstored()` removes the stored entry (`:80-84`).

| value | derivation | comment at head | checked at head |
|---|---|---|---|
| Subgroup MD, oriented | OLS of `y` on 0/1 treatment = `mean(y \| T=1) - mean(y \| T=0)`. Under `adverse_outcome = FALSE` (`?forestsearch` `@param adverse_outcome`: continuous default, higher is better; effect > 0 = harm) the fit's MD is on `-y`, so `-(mean_T - mean_C)`. By hand in `.fco_md()` (`:87-93`). | header `:23-31`; `:86` | `:147`, `:173-174` |
| Sign | `sign()` of the hand value, and of `max_sg_est` for the selected subgroup | `:143` | `:148`, `:175-176` (`expect_identical`) |
| Standard error | pooled variance `s2 = (SS0 + SS1)/(n0 + n1 - 2)`, `se = sqrt(s2 (1/n0 + 1/n1))`; not changed by negation | `:24-26`, `:31` | used in `:149-150`, `:177-180` |
| Naive CI | `?compute_frontier_cis` "Naive CI": full-sample Wald, natural scale for differences, `est -/+ qnorm(0.975) se` | `:32-34` | `:149-150`, `:177-180` |
| Split-CI midpoint | `?compute_frontier_cis` "Split CI": centred at the full-sample estimate, so `(split_lcl + split_ucl)/2 = est` | `:35-36` | `:151-152` |
| Tie to the fit | the selected subgroup's hand MD equals `fit$max_sg_est` (members `treat.recommend == 0`, the first row of `out_sg$result`) | `:37-38`, `:132-134` | `:144` |
| Frontier membership | each label evaluated by hand (`.fco_mask()`: `"{name <= number}"` → `name <= number`; `"{name}"` → `name == 1`; leading `!` negates); every label is asserted to be one of these forms first | `:95-98` | `:167-168` |
| Relational | unstored table `identical()` to the stored-FALSE table (both from `compute_frontier_cis()`, not a hand value) | `:185-186`, `:188-190` | `:191-193` |

Hand numbers (`S/evidence/A3/hand_values.R` → `hand_values.txt`; base R, data from the test helper, no package `R/` code loaded):

| m | labels | n (n1/n0) | mean_T | mean_C | est = -(mean_T - mean_C) | se | naive CI |
|---|---|---|---|---|---|---|---|
| 1 (selected) | `{age <= 47} & {biomarker_hi}` | 36 (19/17) | -1.42388211 | -0.329769888 | 1.09411222 | 0.29131059 | (0.52315395, 1.66507049) |
| 3 | `{biomarker_hi} & {sex}` | 80 (44/36) | -0.58731687 | 0.197397793 | 0.78471467 | 0.24044426 | (0.31345257, 1.25597676) |
| 5 | `{age <= 62} & {biomarker_hi}` | 111 (56/55) | -0.55120304 | 0.034751579 | 0.58595462 | 0.20987338 | (0.17461035, 0.99729889) |

- Tolerance `1e-10` (`:46`), set before the first run, with its reason at `:40-42`: OLS by QR and the hand difference of means agree to rounding (order 1e-15); an orientation error is `2 |MD|`, of order 1.
- Verifier re-derivation (`verify/hand_verify.R`, its own label parser, base-R arithmetic): the same numbers to 10 digits. m = 1 hand − `max_sg_est` = -2.2e-16; `lm(I(-y) ~ treat)` gives 1.09411222057. Fixed tree, entry removed: max |table − hand| 1.3e-15 (estimate), 9.4e-16 / 1.8e-15 (limits), 1.4e-15 (Split midpoint). Start tree: every difference 2.188224 (= 2 × 1.094112). The set from the `sg.harm` labels is `identical()` to `treat.recommend == 0` (n = 36). Walls 3.05 s (fixed), 3.26 s (start). It found no expected value taken from `compute_frontier_cis()`; the cited help texts match the roxygen at `b06de4e8`.

### A3.6 Unaffected case

Task: "For A3, any object that stores the value gives identical output." Each tree in its own `Rscript --vanilla` process, saved with `saveRDS()`, compared with `identical()` in a third process. Per case: the returned data.table, the `verbose = TRUE` printed text, the warnings (and the error text where one is raised). No timing leaf in these objects; nothing dropped.

**Implementer** (`evidence/A3/unaffected.R`, `compute_frontier_cis(n_splits = 50L, seed = 1L, verbose = TRUE)`; start export vs `S/wt/A3`, plus a second start run as control):

| case | data | stored | frontier rows | start vs fixed |
|---|---|---|---|---|
| survival | `.make_survival_data(N = 400, HR_harm = 3)` | FALSE (rule) | 3 (1 finite) | TRUE |
| continuous, default | N 300, seed 42, `y` negated | FALSE | 3 | TRUE |
| continuous, `adverse_outcome = TRUE` | N 300, seed 42 | TRUE | 3 | TRUE |
| binary OR, default | `.make_binary_data(N = 400, OR_harm = 6, seed = 42)` | TRUE | 10 | TRUE |
| binary RD, `adverse_outcome = TRUE` | same | TRUE | 10 | TRUE |
| count IRR, default (offset `ftime`) | `.make_count_data(N = 400, IRR_harm = 3, seed = 42)` | TRUE | 8 | TRUE |
| extension: binary OR, `adverse_outcome = FALSE`, `y` flipped | same | FALSE | 9 | TRUE |
| extension: survival / binary OR / count, entry removed | as above | none | 3 / 10 / 8 | TRUE / TRUE / TRUE |
| intended: continuous, entry removed | N 300, seed 42 | none | 3 | FALSE |

- 10 of 10 non-intended cases identical (`unaffected_cmp.log`); start vs start control identical in all 11 (`unaffected_cmp_control.log`).
- Intended difference (`unaffected_detail.txt`): an exact sign flip of `estimate`, limits swapped and negated. On the fixed tree the entry-removed table is `identical()` to the stored-FALSE table; on the start tree it is not.
- Walls: start 16.26 s, fixed 15.59 s, start control 15.82 s (in R 16.07 / 15.41 / 15.63 s).

**Verifier** (re-ran the implementer's script on its own copies: start control 16.21 s, fixed 16.89 s; 10 of 10 identical, `verify/impl_unaff_cmp.log`; its fixed RDS `identical()` to the implementer's). Its own script `verify/unaff_v.R` (`n_splits = 40L`, `seed = 3L`), start control 45.66 s, fixed 43.29 s, compare 0.13 s (`verify/unaff_cmp.log`):

| case | stored | rows / finite est | start vs fixed |
|---|---|---|---|
| survival on `survival::gbsg` (treat = `hormon`), default / `adverse_outcome = TRUE` | FALSE / TRUE | 10/10, 10/10 | TRUE, TRUE |
| continuous (N 400, `MD_harm` 1.2, seed 7), default / TRUE | FALSE / TRUE | 10/10, 10/10 | TRUE, TRUE |
| continuous, stored `NA` / `"TRUE"` / `logical(0)` | odd | error, same text on both trees | TRUE ×3 |
| continuous, stored `0` / `1` | 0 / 1 | 10/10 | TRUE ×2 |
| continuous stored FALSE, `fs$outcome_type` removed | FALSE | 10/10 | TRUE |
| binary RD (seed 9), `adverse_outcome = FALSE` | FALSE | 10/10 | TRUE |
| binary RR, `adverse_outcome = TRUE` | TRUE | error "no 'restart' 'muffleWarning' found" (pre-existing A-056) on both trees | TRUE |
| binary OR, default | TRUE | 14/14 | TRUE |
| count IRR (seed 13), `adverse_outcome = FALSE` / count IRD default | FALSE / TRUE | 10/10, 15/15 | TRUE, TRUE |
| unstored: survival gbsg, binary OR, binary RD (was FALSE), count IRD | none | full | TRUE ×4 |
| unstored: continuous (intended) | none | 10/10 | FALSE |

- "ALL stored cases identical: TRUE". Intended (`verify/unaff_intended.log`): start estimates are the exact negation of fixed on all 10 rows; start `naive_lcl` = −fixed `naive_ucl`.
- Try-to-break (`verify/break.R`, entry removed, both trees, `break_cmp.log`): `outcome_type` removed (with and without `effect_measure`), the unresolved 4-vector with `"HR"`, `NA`, `"Continuous"`, a continuous fit relabelled binary/OR. All identical, including error texts already raised at the start SHA. The only difference is the intended continuous case. No case becomes an error.

**Main-session cumulative check at head** (`S/cumulative/A3/`; `S/head_fixes` has `R/` and `tests/` equal to `fb80dd8f`, checked for this section):

| script | start wall | head wall | result |
|---|---|---|---|
| `unaffected.R` | 16.92 s | 16.68 s | 10 non-intended identical; intended FALSE (`unaffected_cmp.log`, 0.15 s) |
| `unaff_v.R` | 48.54 s | 43.74 s | all stored identical; 4 non-continuous unstored identical; intended FALSE (`unaff_v_cmp.log`, 0.15 s) |
| `break.R` | 3.68 s | 3.94 s | all non-intended identical (`break_cmp.log`, 0.15 s) |

- `cross_cmp.log`: head vs the A3-only fixed tree, whole objects identical (11/11, 20/20, 7/7 cases); start vs the earlier start runs identical. So the later part-5 commits do not change A3's outputs.
- A head-minus-A3 tree was also built; `unaffected.R` ran on it (17.35 s). No comparison of it is saved, and its `unaff_v` log is empty.

### A3.7 Existing tests

- No existing test file calls `compute_frontier_cis()` or its only package caller `compare_selection_rules()` (`b06de4e8:R/compare_selection_rules.R:597`): `git grep` over `tests/` at `b06de4e8` finds nothing; at `fb80dd8f` only the new file. `make_effect_estimator()` is not changed.
- Neighbour checks on the fixed copies: implementer `test-display-orientation.R` 93 expectations, 0 failed, 6.86 s (`evidence/A3/existing.txt`); verifier `test-display-orientation.R` 93/0 (6 tests), 7.20 s, and `test-sg-focus-transparency.R` 64/0 (13 tests), 8.24 s (`verify/existing.txt`).
- Full suite at `fb80dd8f`: `✔ | 31 | frontier-cis-orientation` (`S/suite/suite.log:561`). The log was still running (in `mr-inference`) when this section was drafted.
- Expectations updated: none.

### A3.8 What a user sees

`evidence/A3/example.R` on the start export and on `S/wt/A3` (`example.txt`). Continuous fit at the default orientation, `max_sg_est = 1.09411`, `sg.harm = {age <= 47} & {biomarker_hi}`, entry removed, `compute_frontier_cis(fit, n_splits = 200L, seed = 1L)`:

```
tree: /tmp/cran_prep_displays_se_docs_2026-11-12/start
       m estimate naive_lcl naive_ucl split_lcl split_ucl fsbc_estimate
1:     1  -1.0941   -1.6651   -0.5232   -1.7340   -0.4542       -1.0925
2:     3  -0.7847   -1.2560   -0.3135   -1.2837   -0.2857       -0.7870
3:     5  -0.5860   -0.9973   -0.1746   -0.9654   -0.2065       -0.5888

tree: /tmp/cran_prep_displays_se_docs_2026-11-12/wt/A3
       m estimate naive_lcl naive_ucl split_lcl split_ucl fsbc_estimate
1:     1   1.0941    0.5232    1.6651    0.4542    1.7340        1.0925
2:     3   0.7847    0.3135    1.2560    0.2857    1.2837        0.7870
3:     5   0.5860    0.1746    0.9973    0.2065    0.9654        0.5888
```

(The `Key: <m>` and column-type lines of the printout are left out.)

Before, the selected subgroup's row (m = 1) was the negative of the fit's own estimate; after, it equals it.

### A3.9 Verification, repair, polish

- **Verdict:** CONFIRMED (`evidence/A3/verify/VERIFY.md`). No blocking problem. No repair round (`workflow1_result.json`: `repaired` false).
- **Polish** (one round, comments only; `evidence/A3/NOTES.md` "Polish before commit"; delta `evidence/A3/polish.diff`, one file, +6/-2): fixed problems 1 and 2. Comment gate printed nothing; `git diff verified -- R` empty. Polish check: CONFIRMED, `problems` [] (`polish_result.json`; `verify_polish/delta_tree.diff` shows only the test file differs).

| # | verifier problem (minor) | what happened |
|---|---|---|
| 1 | Test `:10` (verified) said "the 'Resolve adverse_outcome default' blocks" (plural); only `R/forestsearch_main.R:2236-2243` carries that comment, the Section 1A resolution `:1840-1842` has none | Fixed in the polish: header now reads "in Section 1A, and again in the block under the 'Resolve adverse_outcome default' comment" (`fb80dd8f` test `:10-12`) |
| 2 | Test `:187` (verified) compares two `compute_frontier_cis()` outputs, not a hand value; acceptable as supplementary | Fixed in the polish: three-line comment added (`:188-190`) saying it is relational and that the hand-value checks above carry the gate |
| 3 | The implementer's survival unaffected case has only 1 of 3 finite rows (pre-existing `.label_to_expr()` `!{stage}` → `stage == 0`) | No change; the verifier's `survival::gbsg` case (10/10 finite rows, identical) closes the gap; the label issue is an open item |
| 4 | `compute_frontier_cis()` on a binary RR fit errors "no 'restart' 'muffleWarning' found" (`.estimate_rr()`, `b06de4e8:R/glm_effect_estimators.R:681-687`); pre-existing A-056, identical on both trees | No change; already recorded as A-056 |

### A3.10 Departures from the task's wording

- **The test removes the stored entry.** The task's test is "a continuous-outcome fit at the default orientation". Every `forestsearch()` fit stores the resolved value, so such a fit takes the stored path and would pass on both trees. The test removes the entry (`fit$args_call_all$adverse_outcome <- NULL`) and first asserts the premise that `FALSE` was stored (`:116`) and that the entry is gone (`:124-125`).
- **Placement.** The fallback is a block before the call (`:140-147`) rather than inline, as in F12. The stored value is read by the same `args$adverse_outcome`.
- **Rule input.** The rule is applied to the function's resolved `outcome_type` (`fs$outcome_type %||% "survival"`, `:136`). With no stored `outcome_type` it gives FALSE, which the Cox estimator does not read.
- **More checks than named.** Besides the selected subgroup's estimate, the test checks every frontier row (estimate, sign, naive limits), the Split-CI midpoint, and one relational `identical()`.
- **Unaffected case extended.** Binary with stored FALSE; survival, binary and count with the entry removed (implementer); `survival::gbsg`, odd stored values, rate measures and try-to-break cases (verifier).
- **No roxygen change.** The fix needs none; the help page does not say how the estimate is oriented (open item).
- The task's site and wording match the code: `b06de4e8:R/frontier_cis.R:150` is `adverse_outcome   = args$adverse_outcome %||% TRUE,`.

### A3.11 Articles and examples (step 2)

- `vignettes/articles/mrct-consistency.qmd:632` (both SHAs) evaluates `compute_frontier_cis(fs, n_splits = 1000L, seed = seed_rep5)`. `fs` is a survival `forestsearch()` fit (`:510-523`: `event.name`, no `outcome_type`), which stores `adverse_outcome`. The call runs `:144-145` and `:158` but not the new branch `:146`. Its output cannot change. No other article source calls `compute_frontier_cis()` or `compare_selection_rules()`.
- `\dontrun` examples at `fb80dd8f`: `R/pareto_frontier_table.R:135` passes the `forestsearch()` fit built at `:126`, which stores the value. `R/explain_pareto_selection.R:71` passes an `fs` the example does not build.

### A3.12 NEWS entry (`fb80dd8f:NEWS.md:425-430`, "Behaviour changes"; the inexact sentence is `:427-429`)

> `compute_frontier_cis()` on a fit object that stores no `adverse_outcome` now orients the effect as `forestsearch()` does (`TRUE` for binary and count outcomes, `FALSE` for continuous); it used `TRUE`. For a continuous outcome the table's estimates and limits now have the sign of the fit's own estimate. A stored value is used as before, so objects built by `forestsearch()` are unchanged.

| claim | check | verdict |
|---|---|---|
| "on a fit object that stores no `adverse_outcome`" | branch is `is.null(args$adverse_outcome)` (`:144-145`) | accurate |
| "orients the effect as `forestsearch()` does (TRUE for binary and count, FALSE for continuous)" | `:146` = `R/forestsearch_main.R:1840-1842`. Survival also gets FALSE and count gets TRUE, but neither estimator reads it; only continuous output changes | accurate |
| "it used `TRUE`" | `b06de4e8:R/frontier_cis.R:150` | accurate |
| "For a continuous outcome the table's estimates and limits now have the sign of the fit's own estimate" | The fit's own estimate is `max_sg_est`, for the selected subgroup only; other rows need not share its sign. A limit need not share any estimate's sign (an interval can cross 0). What changes is that each estimate and limit is negated (limits swapped): the table is now oriented like the fit's estimate. It also holds only when the fit's estimate was made at the default orientation (FALSE). The implementer's proposed bullet said "at the default orientation"; the NEWS dropped it | **inexact** |
| "A stored value is used as before" | same `args$adverse_outcome`; verifier's odd stored values identical | accurate |
| "so objects built by `forestsearch()` are unchanged" | every `forestsearch()` since `99ced431` stores a logical; older objects have no frontier and return an empty table before the estimator on both trees (git and code reading, not run); current fits: unaffected cases identical | accurate |

### A3.13 Evidence discrepancies

- **History of the stored value.** `VERIFY.md` §1: "before `bc4690b3` (0.2.0 dev), capture (`:571`) preceded resolution (`:691`). Unset fits from those versions stored `adverse_outcome = NULL`. The NEWS wording 'one saved by an older version' is therefore accurate." Git: `bc4690b3^:R/forestsearch_main.R:693` is `args_call_all$adverse_outcome <- adverse_outcome`, right after the resolution at `:690-692`. The same write-back is in every version since `adverse_outcome` was added (`99ced431:R/forestsearch_main.R:372`). Objects from before `99ced431` have no `pareto_frontier` (first in `bdef2595`), so `compute_frontier_cis()` returns early (`b06de4e8:R/frontier_cis.R:113-118`). The "saved by an older version" wording is in `NOTES.md` §2 and §10, `example.R:2`, the proposed NEWS bullet in `workflow1_result.json`, and the committed test comment `fb80dd8f:tests/testthat/test-frontier-cis-orientation.R:19-20`. The final NEWS does not carry it.
- **Line range of the hand-value checks.** `VERIFY.md` (problem 2) gives ":146-153 and :172-179" (verified numbering). The polish (`NOTES.md`, `polish_result.json`) gives "verified :143-151 and :172-179". In the verified file the selected-subgroup hand checks are `:143-151` and the test closes at `:152`; `:146-153` starts at the first table check and runs one line past the closing brace. Head: `:144-152` and `:173-180`.
- **Number of `%||%` definitions.** `VERIFY.md` §1: "Ten files under `R/` define it." `git grep` at `b06de4e8`: 12 files, 9 at top level and 3 inside functions (`R/compare_selection_rules.R:203`, `R/plot_pareto_combined.R:419`, `R/run_simulation_analysis.R:28`). The verifier's conclusion stands: the bound definition is `if (is.null(a)) b else a`, as its runtime check and the identical error texts for `logical(0)` show.

### A3.14 Open items

- `?compute_frontier_cis` does not say how the estimate is oriented (the stored `adverse_outcome`, or `forestsearch()`'s rule). The roxygen of `R/frontier_cis.R` is unchanged from `b06de4e8` to `fb80dd8f` apart from A3's code lines; C4 did not cover it. Documentation only.
- Comment: the test header (`fb80dd8f:tests/testthat/test-frontier-cis-orientation.R:19-20`) says the fallback is reached by "an object saved by an older version"; per git (A3.13), only an object built or edited by hand reaches it. Comment only.
- Out of scope, pre-existing: with no stored `effect_measure`, the binary fallback is `"RD"` (`R/frontier_cis.R:137-139`, both SHAs), while `forestsearch()`'s binary default is `"OR"` (`b06de4e8:R/forestsearch_main.R:1832-1838`). `forestsearch()` always stores `effect_measure` for non-survival fits.
- Out of scope, pre-existing: `.label_to_expr()` turns a label on a multi-level factor (`"!{stage}"`, levels 1/2/3) into `stage == 0`, which matches no row. `compute_frontier_cis()` then returns an all-NA row with no warning (the `n_sub < 10` branch, `b06de4e8:R/frontier_cis.R:214-217` → `fb80dd8f:R/frontier_cis.R:222-225`). Reproducer `evidence/A3/explore/surv_na.R`; verifier `verify/explore/surv_na_rerun.txt` (2 of 3 survival frontier rows NA, including the selected subgroup). Identical on both trees.
- Out of scope, pre-existing: binary RR fits error "no 'restart' 'muffleWarning' found" (A-056, `b06de4e8:R/glm_effect_estimators.R:681-687`).

## B1 — Risk-difference standard error by the full delta method, covariate-adjusted logistic tier (change of method)

`S` = `/tmp/cran_prep_displays_se_docs_2026-11-12`. Line numbers are at the start SHA `b06de4e8` and at head `fb80dd8f`. `R/` is identical at `835fac1b` and `fb80dd8f`.

- `R/glm_effect_estimators.R` is the same at `b06de4e8` and `fb21396f` (B1's parent): `git diff --stat b06de4e8 fb21396f` lists neither B1 file.
- After B1, one commit touches the R file: C4 `ba7c4ea2` rewrites the `@param adverse_outcome` roxygen of `make_effect_estimator()` (`:53-`, +23/-5). Every B1 line therefore sits **18 lines lower at head** than at the commit `8e87cf30`. Both are given below where they differ.
- The test file is unchanged from `8e87cf30` to `fb80dd8f` (`git diff --stat 8e87cf30 fb80dd8f` lists only the R file), so the test line numbers below hold at the commit and at head.

### B1.1 Status, commit, provenance

| | |
|---|---|
| Status | **PASS** (change of method, formula approved by Larry in the task) |
| Commit | `8e87cf30` "fix: covariate-adjusted risk-difference SE by the full delta method" |
| Parent | `fb21396f` (A3) |
| Files (`git show --numstat 8e87cf30`) | `R/glm_effect_estimators.R` +75/-4; `tests/testthat/test-rd-adjusted-se.R` +303 (new). 2 files, +378/-4 (`S/main/B1/diffstat.txt`). No `man/` or `NAMESPACE` change: the helper is `@noRd`, and `devtools::document()` wrote nothing in six separate runs (B1.8). |
| Commit vs gated patch | `git diff b06de4e8 8e87cf30 -- <both files>` equals `S/patches/B1.patch` once the `index` lines are dropped (checked for this section; the only difference is the hash abbreviation length). `S/main/B1/patch_equal.txt` is empty. |
| Patch history | `S/patches/B1.r1.patch` (round 1, sha256 `fcf50a8f…0a458`, R +56/-3, test 254 lines) → `S/patches/B1.verified.patch` (repair round, sha256 `f1d92a25…cad8fac`, R +68/-3, test 303 lines) → `S/patches/B1.patch` (polish plus the main-session correction, sha256 `9af8cdd3…ecc921c`, R +75/-4, test 303 lines; the test blob `ad7473e` is the same in the verified and the final patch). |
| Evidence | Implementer `S/evidence/B1/NOTES.md` (round 1; repair round `r2/`; polish `polish/`, `polish.diff`; main-session correction). Verifiers `S/evidence/B1/verify/VERIFY.md` (round 1), `S/evidence/B1/verify_r2/VERIFY.md` (round 2). Polish check `S/evidence/B1/verify_polish/` and `S/polish_result.json` (no VERIFY.md was written). Main session `S/main/B1/`. Workflow `S/workflow1_result.json` (status READY, final verdict CONFIRMED, repaired). |

### B1.2 Site, before → after

Function `.estimate_rd()`, tier 2 (logistic g-computation). At `b06de4e8` the function starts at `:457`; tier 2 is `:514-561`; the comment is `:531-532`; the adjusted branch is `:545-547`, with the expression at `:546`. The task's pointers ("near 534-547", "near 531-532") match by content.

Before, `b06de4e8:R/glm_effect_estimators.R:531-532` and `:545-547`:

```r
    # SE via coefficient SE as approximation (exact delta method is complex
    # with covariates; this is conservative for subgroup estimation)
...
    } else {
      abs(avg_deriv) * se_treat
    }
```

After, `fb80dd8f:R/glm_effect_estimators.R:549-553` and `:566-578` (`8e87cf30:531-535`, `:548-560`):

```r
    # SE of rd_est by the delta method: the two-sample binomial SE when
    # unadjusted, sqrt(g' V g) over every estimable coefficient when
    # adjusted (.fs_rd_gcomp_se()).  The scaled-coefficient form
    # abs(avg_deriv) * se_treat, which uses the treatment coefficient alone,
    # is kept only as the adjusted branch's fallback.
...
    } else {
      # Adjusted (B1): the full delta method.  A value that is not a finite,
      # non-negative number keeps the scaled-coefficient form (NA when the
      # treatment coefficient is aliased, i.e. the effect is not
      # estimable); the helper never errors, so no slice is routed to tier
      # 3 by it.
      se_dm <- .fs_rd_gcomp_se(fit_logit, nd0, nd1, pred0, pred1)
      if (length(se_dm) == 1L && is.finite(se_dm) && se_dm >= 0) {
        se_dm
      } else {
        abs(avg_deriv) * se_treat
      }
    }
```

New helper `.fs_rd_gcomp_se()`, `fb80dd8f:R/glm_effect_estimators.R:611-666` (roxygen `:611-639`, `@noRd`; body `:640-666`) = `8e87cf30:593-648`. It sits directly after `.estimate_rd()`, before the OR section. Body:

```r
.fs_rd_gcomp_se <- function(fit, nd0, nd1, p0, p1) {
  tryCatch({
    b   <- stats::coef(fit)
    est <- names(b)[!is.na(b)]
    al  <- names(b)[is.na(b)]
    tt  <- stats::delete.response(stats::terms(fit))
    design <- function(nd) {
      mf <- stats::model.frame(tt, nd, na.action = stats::na.pass,
                               xlev = fit$xlevels)
      stats::model.matrix(tt, mf, contrasts.arg = fit$contrasts)
    }
    X0 <- design(nd0)
    X1 <- design(nd1)
    # Treatment changes an aliased column (stop -> NA; see above).
    stopifnot(!any(X1[, al, drop = FALSE] != X0[, al, drop = FALSE],
                   na.rm = TRUE))
    X0 <- X0[, est, drop = FALSE]
    X1 <- X1[, est, drop = FALSE]
    stopifnot(length(p0) == nrow(X0), length(p1) == nrow(X1))
    V <- stats::vcov(fit)[est, est, drop = FALSE]
    g <- colMeans(p1 * (1 - p1) * X1) - colMeans(p0 * (1 - p0) * X0)
    v <- sum(g * drop(V %*% g))
    if (is.finite(v) && v >= 0) sqrt(v) else NA_real_
  },
  error   = function(e) NA_real_,
  warning = function(w) NA_real_)
}
```

- `g` and `V` are the task's formula. "Estimable columns only" is done by dropping the `NA` (aliased) coefficients by name from `b`, `X0`, `X1` and `V`.
- `X0`, `X1` are built as `predict.lm()` builds them (the fit's terms, `xlevels`, `contrasts`), so their rows match `pred0`/`pred1`, the probabilities the estimate averages.
- The helper returns `NA_real_` (and the caller keeps the old expression) when: the value is not finite and non-negative; an aliased column differs between `X1` and `X0` (`:654-655`, added in the repair round); or any step errors or warns.
- File-header list of internal helpers, `fb80dd8f:21` (same line at `8e87cf30`): `#   - .estimate_rd, .estimate_or, .estimate_rr` → `#   - .estimate_rd, .estimate_or, .estimate_rr, .fs_rd_gcomp_se` (added in polish).
- Unchanged (checked against the code at both SHAs): the unadjusted branch (F13, start `:536-544`, head `:557-565`); the estimate `rd_est`; `converged`; `method_used`; tiers 1 and 3; the IPTW path; `.estimate_or()`, `.estimate_rr()` and the Poisson estimators. `se_treat` and `avg_deriv` are still computed before the branch (head `:554`, `:556`), so an error there still routes to tier 3 as before.

**When the changed lines are reached** (implementer NOTES §8, checked by both verifiers at `b06de4e8`):

- Only through `make_effect_estimator(outcome_type = "binary", effect_measure = "RD", adjust_covariates = <non-empty>, ps_adjust_method = "none")`. With `"dr_gcomp"`, or `"iptw"` without an `sw` column, the covariates are dropped before the call, so `n_adj` is 0 and the slice gets F13's unadjusted branch.
- A slice then reaches tier 2 when the identity-link fit fails: an arm at 0% or 100% events; a factor or character covariate when the slice holds 3 or more of its levels (tier 1's start vector has one value per covariate name, the design one column per non-reference level); or an identity-link failure reported as an error or as a warning containing "converge". (The round-1 wording "on every slice" was corrected in the repair round; see B1.9.)
- No default path reaches it: `forestsearch()` has `adjust_covariates = NULL`, and the binary `effect_measure` resolves to `"OR"`.
- What the SE feeds (callers table, NOTES §8, line numbers at `b06de4e8`):
  - `forestsearch()` (`forestsearch_main.R:2263-2272`): candidate `L(HR)`/`U(HR)` (`subgroup_search.R:1080-1089`). These are also near-duplicate keys (`subgroup_consistency_helpers.R:1120-1131`) and exact-tie keys in `remove_redundant_subgroups()` (`:1159`; verifier round 2).
  - Bootstrap: observed H/Hc SE (`bootstrap_dofuture_main.R:310-320`, `:346-348`) → `ci_est()` (`bootstrap_calculations_helpers.R:286`). The replicates' inner `forestsearch()` runs (`bootstrap_analysis_dofuture.R:615`) use it through the candidate `L(HR)`/`U(HR)` (verifier round 2; NOTES says "Replicates use the estimate only").
  - `SG_tab_estimates_glm()` "RD (95% CI)" cell via `sg_tables()` (`summary_utility_functions.R:830`, `:983`), the CV tables (`forestsearch_cross_validation.R:1510`, `:1566`) and `mr_estimates_table()` (`mr_estimates_table.R:211`).
  - `plot_sg_glm_outcomes()` (`:262`) and `plot_subgroup_results_forestplot()` (`:444`).
  - The exported `subgroup_glm()` (`run_subgroup_sims.R:321`), whose upper bound is `est + z * se` (`:364-373`; head `:371`).
  - Not reached: `frontier_cis.R:143` (`adjust_covariates = NULL`), MR inference, the consistency resample, Guo-He.
- Screening and the consistency splits use the estimate only (`subgroup_consistency_helpers.R:300-306`).

### B1.3 Red and green runs

Test file `tests/testthat/test-rd-adjusted-se.R` (new, 303 lines, 7 `test_that` blocks, 35 expectations). Command for every summary run: `NOT_CRAN=true timeout 300 Rscript --vanilla -e "devtools::test('<copy>', filter = 'rd-adjusted-se', reporter = 'summary', stop_on_failure = FALSE)"`. Per-block counts: `counts.R` (silent reporter).

**Main session, on the repository** (`S/main/B1/`): red = the commit's parent `fb21396f` plus the test file only; green = with the code.

| Run | Result | Wall |
|---|---|---|
| Red (`red.txt`) | `..1.234..56.78.9`: 7 passes, 5 failures, 4 errors (12 expectations) | 3.0 s |
| Green (`green.txt`) | 35 passes, 0 failures, 0 errors | 2.9 s |

- First failure, `test-rd-adjusted-se.R:111:3`: "Expected `res$se` to equal `.rdadj_fd_se(d, ~treat + site + x)`." `actual`: 0.1291, `expected`: 0.1217.
- Other failures: `:124:3` and `:127:3` (aliased covariate, 0.1291 vs 0.1217); `:207:3` and `:210:3` (treated arm at 100%, actual 0.000 vs 0.078).
- Errors: `:161:3`, `:243:3` and `:294:5`, `could not find function ".fs_rd_gcomp_se"`; `:232:5`, `local_mocked_bindings(...)`: "Can't find binding for `.fs_rd_gcomp_se`".

**Implementer and verifiers:**

| Round | Run | Copy | Expectations | Failed | Errors | Warnings | Wall (summary / counts) | Evidence |
|---|---|---|---|---|---|---|---|---|
| Implementer, round 1 | red | `wt/B1_red` (start + test) | 11 | 5 | 3 | 0 | 3.30 s / 2.85 s | `red.txt`, `red_counts.txt` |
| Implementer, round 1 | red, first draft of the test | — | 11 | 5 | 3 | 0 | — | `red_first_draft.txt`, `red_counts_first_draft.txt` |
| Implementer, round 1 | green | `wt/B1` | 27 | 0 | 0 | 0 | 2.81 s / 3.13 s | `green.txt`, `green_counts.txt` |
| Verifier, round 1 | red | `verify/B1/red` | 11 | 5 | 3 | 0 | 3.14 s / 3.11 s | `verify/red.txt`, `red_counts.txt` |
| Verifier, round 1 | green | `verify/B1/fixed` | 27 | 0 | 0 | 0 | 2.88 s / 2.76 s | `verify/green.txt`, `green_counts.txt` |
| Implementer, repair | red | `wt/B1_r2_red` | 12 | 5 | 4 | 0 | 3.00 s / 2.85 s | `r2/red.txt`, `r2/red_counts.txt` |
| Implementer, repair | new test on round-1 code | `wt/B1_r2` before the repair | 35 | 4 | 0 | 0 | 2.97 s / 2.94 s | `r2/r1code_newtest*.txt` |
| Implementer, repair | green | `wt/B1_r2` | 35 | 0 | 0 | 0 | 2.82 s / 2.80 s | `r2/green.txt`, `r2/green_counts.txt` |
| Verifier, round 2 | red | `verify/B1_r2/red` | 12 | 5 | 4 | 0 | 2.92 s / 2.66 s | `verify_r2/red*.txt` |
| Verifier, round 2 | new test on round-1 code | `verify/B1_r2/r1code_newtest` | 35 | 4 | 0 | 0 | 2.89 s / 2.74 s | `verify_r2/r1code_newtest*.txt` |
| Verifier, round 2 | green | `verify/B1_r2/fixed` | 35 | 0 | 0 | 0 | 2.87 s / 2.67 s | `verify_r2/green*.txt` |
| Polish | polished | `wt/B1_polish` | 35 | 0 | 0 | 0 | 3.06 s / 3.16 s | `polish/summary_B1_polish.txt`, `counts_B1_polish.txt` |
| Polish | verified (tag `verified`) | `wt/B1_polish_verified` | 35 | 0 | 0 | 0 | 2.84 s / 2.84 s | `polish/*_verified.txt` |
| Polish | red (start + polished tests) | `wt/B1_polish_red` | 12 | 5 | 4 | 0 | 3.02 s / 2.89 s | `polish/*_red.txt` |
| Polish check | polished | `verify/B1_polish` | 35 | 0 | 0 | 0 | 3.00 s / 2.87 s | `verify_polish/tests_polish.txt`, `counts_polish.txt` |

- Per-block red (repair round on): (a) 3 expectations, 1 failed; (a2) 3, 2 failed; (b) 0, error; (c) 5, 2 failed; safety 0, error; helper-never-errors 0, error; not-estimable 1 (its premise passes), then error at `:294`.
- Per-block green: (a) 3, (a2) 3, (b) 2, (c) 5, safety 8, helper-never-errors 6, not-estimable 8 = 35.
- Round-1 counts (11 red / 27 green) are for the round-1 test file (254 lines, 6 blocks). The repair added the not-estimable block (+8 expectations, +1 error on the start tree).
- The new block fails on the round-1 code, both slices: `:294` helper `0.0` vs `NA`, and `:301` `is.na(res$se)` FALSE. So the test catches the round-1 blocking defect.

### B1.4 Expected values

Every value is worked out by hand or recomputed in the test from a separately fitted `glm()`; none is taken from the fixed code. Each tolerance was written before the first run, with its reason in a comment. Line numbers are in `tests/testthat/test-rd-adjusted-se.R` (same at `8e87cf30` and `fb80dd8f`).

| Block | Assertion (line) | Expected value | Derivation | Tolerance (reason) | Comment |
|---|---|---|---|---|---|
| (a) `:101`, factor `site` (3 levels) + `x`, through `make_effect_estimator()` | premise `method_used == "logistic_gcomp_adj"` (`:108`), `converged` TRUE (`:109`) | `"logistic_gcomp_adj"`, TRUE | Tier 1's start vector has 4 values for 5 design columns, so `glm()` stops and the slice goes to tier 2 | exact | `:105-107`; data `:40-42` |
| (a) | `res$se` (`:111`) | `.rdadj_fd_se(d, ~ treat + site + x)`: `sqrt(g' V g)`, `g` by central finite differences (h = 1e-5) of `mean(plogis(X1 b)) - mean(plogis(X0 b))`, explicit `model.matrix()`, `V = vcov(fit, complete = FALSE)` | Recomputed in the test (`:76-95`) | 1e-6 relative. Truncation `h^2/6 * |f'''|` with `|f'''| <= 2 * 0.125 * 3.375 = 0.84`, so <= 1.4e-11; rounding about 1e-11; gradient entries 1e-2 to 1e-1 | `:58-74` |
| (a2) `:115`, as (a) + aliased `x2 = 2x` | premise (`:123`); `res$se` vs FD on `~ treat + site + x + x2` (`:124`) and on `~ treat + site + x` (`:127`) | the (a) value both times | `x2` is in the span of `x`: same fitted model, same SE; the FD helper drops the NA column | 1e-6 | `:116-118`, `:126` |
| (b) `:135`, helper on an unadjusted fit, control 12/40, treated 15/30 | `.fs_rd_gcomp_se()` (`:161`) | `sqrt(0.3*0.7/40 + 0.5*0.5/30) = sqrt(0.0135833) = 0.1165476` | By hand: `g' V g` written out for the saturated two-arm model (`V_aa = 1/(n0 p0 q0)`, `V_at = -V_aa`, `V_tt = V_aa + 1/(n1 p1 q1)`) gives `p0 q0/n0 + p1 q1/n1` | 1e-6 (glm convergence; `vcov()` uses the IRLS weights from before the final update; measured 1.9e-7, stated in the comment) | `:146-160` |
| (b) | `round(., 5)` (`:163`) | `0.11655` | From the hand value 0.11654756 | exact | `:154-155` |
| (c) `:170`, binary `z`, treated 40/40, through `make_effect_estimator()` | premise (`:183`), `converged` (`:184`) | `"logistic_gcomp_adj"`, TRUE | Tier 1's start puts the treated arm at probability exactly 1, which `binomial()` rejects | exact | `:181-182` |
| (c) | `res$se` (`:207`) | `sqrt(0.375^2*0.2*0.8/20 + 0.625^2*0.5*0.5/20) = sqrt(0.0060078125) = 0.0775101` | By hand. With the treated arm at 100%, `p1 (1 - p1) -> 0`, so the treated arm adds nothing; the control model in `z` is saturated; weights `w_z` = share of all 80 at `z` (0.375, 0.625) | 1e-6 (treated fitted probabilities about 1e-8 from 1 add about 1e-9 to SE^2 = 6e-3, so at most about 1e-7 in the SE) | `:186-206` |
| (c) | `round(res$se, 6)` (`:210`) | `0.07751` | From the hand value | exact | `:197` |
| (c) | `res$estimate` (`:211`) | `1 - (0.375*0.2 + 0.625*0.5) = 0.6125` | By hand | 1e-6 | `:198` |
| Safety `:218`, helper mocked to `NA`, `NaN`, `Inf`, `-0.5` | premise (`:234`) and `res$se` (`:235`), 4 times | `abs(mean(p1 (1 - p1))) * se(treat)` on a separately fitted `glm(y ~ treat + site + x)` | The old expression, by its definition (`:222-225`) | 1e-12 (same model, same data, same arithmetic) | `:220-221`, `:227-230` |
| Helper never errors `:239` | `fit = NULL` (`:243`); probabilities of the wrong length (`:245`); NA covariate: no error (`:250`), premise (`:251`), NA estimate (`:252`), NA SE (`:253`) | `NA_real_`, `NA_real_`, no error, `"logistic_gcomp_adj"`, NA, NA | The helper's contract; the old expression is NA too | exact | `:242`, `:244`, `:247-248` |
| Not estimable `:260`, slices `empty_arm` (no treated patient) and `trt_all_na` (10 treated, all outcomes NA) | premise 1 `coef(fit)[["treat"]]` NA (`:279`); helper `NA_real_` (`:294`); premise 2 (`:300`); `is.na(res$se)` (`:301`) | NA | By definition: the contrast is not estimable, so it has no SE; the start tree's old expression is NA as well (`vcov()` has an NA row) | exact | `:261-268`, `:276-277`, `:281-289`, `:297-298` |

Measured agreement after green (not used to set anything):

| Value | Fixed code | Independent value | Relative difference | Evidence |
|---|---|---|---|---|
| (a), (a2) | 0.121747444697 | FD 0.121747444694 | 1.85e-11 | `agreement.txt` |
| (a), (a2) | 0.1217474447 | verifier's own Newton-Raphson MLE, no `glm()`/`vcov()`: 0.1217474716 | 2.2e-07 (glm convergence) | `verify_r2/values.txt`, `closure_values.txt` |
| (b) | 0.116547536455 | hand 0.116547558247 | 1.87e-07 | `agreement.txt`, `verify_r2/closure_values.txt` |
| (c) | 0.077510080285 | hand 0.077510079990 | 3.82e-09 (verifier: 3.68e-09 against 0.0775100800) | `agreement.txt`, `verify_r2/closure_values.txt` |
| Start-tree (a) | 0.1291168089 | — | the old expression, the red value | `verify_r2/values.txt` |
| Start-tree (c) | 6.79e-06 | — | the old expression | `verify_r2/closure_values.txt` |

Verifier checks of the derivations:

- Both verifiers re-derived (a), (a2), (b) and (c) without the fixed code: own Newton-Raphson MLE, V the inverse Fisher information, analytic gradient (`verify/values.R` 2.43 s; `verify_r2/values.R` 0.10 s).
- Data claims: control 10/30, treated 16/30, every treat × site cell has both outcomes, `x` in [-1.5, 1.5]; (c) cells control `z=0` 4/20, `z=1` 10/20, treated 10/10 and 30/30. The not-estimable design `(1, treat, z)` has rank 2 of 3 on the fitted rows.
- (c) limit: `glm()` stops at `b_t` = 20.059 with `1 - min p1` = 7.77e-09, iter 18 (`closure_values.txt`). Holding `b_t` finite gives relative errors 1.8e-6 (15), 8.7e-8 (18), 1.2e-8 (20), 8e-11 (25).
- Gradient check by finite differences of the estimate exactly as the closure computes it (`mean(predict(fit_b, nd1)) - mean(predict(fit_b, nd0))`, coefficients perturbed, no design matrix built):
  - round 1, 10 cases, agreement 1e-12 to 2.5e-09 (`verify/fd_predict.txt`, 2.84 s);
  - round 2, 14 cases (treatment interaction, covariate equal to treatment, `data.table`/tibble slices, NA in `y` and in `treat`, double treatment, `splines::ns()`, n = 20000, within-cell separation), agreement 1.1e-11 to 2.5e-10 (`verify_r2/probe_fixed.txt`, 5.94 s).
  - `one_row_per_arm` (4 rows, saturated): fixed 5.0e-06, FD 0, start 4.9e-06, all numerically 0.

### B1.5 Unaffected case

Task: "unadjusted risk difference (the F13 values), odds ratio, risk ratio, and an adjusted risk difference that the identity-link tier fits, are identical."

Each tree ran in its own `Rscript --vanilla` process, results saved with `saveRDS()`, compared with `identical()` in a third process. Dropped before comparing: the timing leaves `time_search`, `minutes_all`, `timing_seconds`, `timing_minutes`, `tmins_search`, `tmins_iteration`, and `row.names` attributes (verifier round 2 also `S7_class`).

| Run | Inputs | Unaffected set, `identical()` | Affected set (start `method_used = "logistic_gcomp_adj"`) | `forestsearch()` | Walls start / fixed / compare |
|---|---|---|---|---|---|
| Implementer, round 1 (`unaffected.R`, `unaffected.txt`) | 260 `make_effect_estimator(outcome_type = "binary")` results (24 configurations × 12 slices, `adverse_outcome` TRUE/FALSE; F13 slices included) + 3 `forestsearch()` | **216 TRUE**: unadjusted RD 48 (F13 slices included), OR/RR 132, `identity_binomial_adj` 8, plus `iptw_rd` 20, `dr_gcomp` RD and 20 pre-existing errors | 44, identical with `se` removed | OR default, RD unadjusted, RD + `age`: TRUE (no slice of RD + `age` reached tier 2) | 5.18 s / 5.32 s / 0.27 s |
| Verifier, round 1 (`verify/unaffected.R`) | 960 results (`survival::gbsg` slices, F13 slices, degenerate and interior; RD/OR/RR/IRR/IRD configurations × `adverse_outcome` × `robust_se`) + 3 `forestsearch()` | **852 TRUE**, `identity_binomial_adj` 92 included; 40 pre-existing errors identical | 108, identical with `se` removed | OR, RD unadjusted: TRUE. RD + `site` (3-level factor, reaches tier 2): differs only in `find.grps$out.found$hr.subgroups` columns `L(HR)`, `U(HR)`; same 32 rows, same HR, same `sg.harm {age <= 64}` | 8.42 s / 8.19 s / 0.22 s |
| Implementer, repair round, both scripts re-run (`r2/unaffected_impl*`, `r2/unaffected_verif*`) | as above | 216 TRUE; 852 TRUE | 44 and 108, identical with `se` removed | as above | 5.13 / 5.20 / 0.26 s; 7.76 / 8.14 / 0.23 s |
| Verifier, round 2 (`verify_r2/unaffected.R`) | 1360 results (synthetic n = 80, 200, 500, 14, rare-event 150; F13; degenerate; gbsg; empty treated arm; all-NA treated outcomes) + 50 display results (`subgroup_glm()`, `SG_tab_estimates_glm()` ITT cell) + 4 `forestsearch()` (N = 320, seed 23) | **1044 TRUE**: unadjusted RD 68 (F13 included), OR 340, RR 340, `identity_binomial_adj` 132; 48 pre-existing errors identical. Not-estimable `logistic_gcomp_adj` (start `se` NA): 56, identical in full | 316 (260 estimable), identical with `se` removed; every fixed `se` finite and `>= 0` where the start `se` was finite | OR default, RD unadjusted, RD + `age`: TRUE. RD + `site` (4-level): differs only in `L(HR)`, `U(HR)`; `sg.harm {sex}` on both | 8.21 s / 8.86 s / 0.29 s |

- Round 1 against round 2 (`r2/r1_vs_r2.txt`, 0.38 s): both unaffected sets, `forestsearch()` runs included, are `identical()` between the round-1 and round-2 fixed results. The probe differs only in `$empty_trt_adj_z$se`, `$empty_trt_adj_site$se`, `$empty_ctl_adj_z$se` (0 in round 1, NA in round 2). The start re-runs of all three scripts are `identical()` to the round-1 start results.
- The polish changed comments and roxygen only (gate: no non-comment line), so the unaffected case was not re-run for it.
- Probes for errors (no case became an error; warning lists the same on both trees):
  - verifier round 1: 21 edge cases (2.61 s / 2.47 s), plus `options(warn = 2)`;
  - verifier round 2: 18 further cases (5.74 s / 5.94 s);
  - polish check: P1-P5 (2.54 s).
- Affected examples (start → fixed SE):
  - degenerate single arm: 451.8 → 0.0648 (`ctl0_trt80` + `x`), 5.48e-06 → 0.0748 (`ctl30_trt100` + `x`);
  - F13 slices with a covariate: 448.5 → 0.0632456 = sqrt(0.004) (`f13_ctl0_trt80`), 5.4e-06 → 0.0724569 (`f13_ctl30_trt100`), the F13 values;
  - 3-level factor, interior slices: 0.1120 → 0.1099, 0.0656 → 0.0626, 0.2233 → 0.2253, 0.0469 → 0.0464;
  - both arms degenerate (`ctl0_trt100`, `all_zero`, `both0`): about 2.3e-07 on both trees.

### B1.6 Existing tests

| Run | Files | Blocks | Expectations | Failed | Errors | Skipped | Warnings | Wall | Evidence |
|---|---|---|---|---|---|---|---|---|---|
| Implementer, round 1, fixed | 15 | 142 | 1085 | 0 | 0 | 5 | 12 | 91.8 s | `existing.txt` |
| Implementer, round 1, start | 15 | 142 | 1085 | 0 | 0 | 5 | 12 | 88.4 s | `existing_start.txt` |
| Verifier, round 1, fixed | 17 | 153 | 1217 | 0 | 0 | 5 | 12 | 87.3 s | `verify/existing_fixed.txt` |
| Verifier, round 1, start (`cross-outcome-parity` only) | 1 | — | 49 | 0 | 0 | 2 | 12 | 46.32 s | `verify/cross_parity_start.txt` |
| Implementer, repair, fixed | 17 | 154 | 1225 | 0 | 0 | 5 | 12 | 83.9 s | `r2/existing_fixed.txt` |
| Implementer, repair, start | 16 | 147 | 1190 | 0 | 0 | 5 | 12 | 83.6 s | `r2/existing_start.txt` |
| Verifier, round 2, fixed | 28 | 260 | 2390 | 0 | 0 | 8 | 16 | 103.8 s | `verify_r2/existing_fixed.txt` |
| Verifier, round 2, start | 27 | 253 | 2355 | 0 | 0 | 8 | 16 | 104.0 s | `verify_r2/existing_start.txt` |

- The 15 files: rd-standard-error, estimability-boundary, binary-default-or, binary-default-or-entry-points, lean-fits, directive-c, threshold-pair-directive-a, fs_dgm_scale, fs-dgm-feasibility, dina-gaussian-floor, glm-pipeline-integration, display-orientation, iptw-continuous, cross-outcome-parity, bootstrap-glm-parent-scale.
- The 17: those plus fs-family-report and rd-adjusted-se.
- The 28 add: no-subgroup-dina-grf-plot, mr-admission-set, mr-admission-rounded, mr_values, betaHhat-exact, betaHhat-contract, get-fsdata-diagnostics, missing-membership-messages, fpr-calibration-fallback, knoise-population, save-leff-calibration. The verifier chose these as every test file that calls `make_effect_estimator()`, `.estimate_rd()`, `subgroup_glm()` or `SG_tab_estimates_glm()`, uses `effect_measure = "RD"`, or runs a binary `forestsearch()`.
- Per-file counts of the common files are identical start vs fixed in every pair. The warnings are the same on both trees: 12 from `cross-outcome-parity`, 4 from `get-fsdata-diagnostics`.
- 1225 = 1217 + 8 new expectations; 2390 = 2355 + 35.
- **Expectations updated: none.**
- The main session ran no per-item existing-test run for B1. The full suite at head is reported under step 4.

### B1.7 What a user sees

`S/evidence/B1/example.R`, output `example.txt`. Each tree in its own process. Start = the start export; fixed = `wt/B1` (round 1). The repair changed only the not-estimable case (`r2/r1_vs_r2.txt`) and the polish only comments, so these values hold at the commit. Test (c) pins case (1) at head.

```
(1) control 14/40, treated 40/40, adjust_covariates = "z"
start:  estimate 0.61249999587013  se 6.79127662230129e-06  logistic_gcomp_adj
        95% Wald CI: 0.6125 (0.6125, 0.6125)
fixed:  estimate 0.61249999587013  se 0.0775100802854232    logistic_gcomp_adj
        95% Wald CI: 0.6125 (0.4606, 0.7644)

(2) control 10/30, treated 16/30, adjust_covariates = c("site", "x");
    SG_tab_estimates_glm() ITT row, "RD (95% CI)"
start:  se 0.129116774183919   1 0.200 (-0.053, 0.453)
fixed:  se 0.121747444696635   1 0.200 (-0.039, 0.439)
```

Walls 2.42 s (start), 2.59 s (fixed).

- A control arm at 0% blew the old SE up instead (verifier round 2, display results, `f13_ctl0_trt80` + `age`):
  - `SG_tab_estimates_glm()` cell `0.8044 (-830.7797, 832.3885)` → `0.8044 (0.6838, 0.9251)`;
  - `subgroup_glm()` `0.8044286 832.3732325` → `0.8044286 0.9251011`.
- Not-estimable slice, 30 patients, none treated, `adjust_covariates = "z"` (`r2/empty_arm.txt`; verifier `verify/empty_arm.txt`):
  - start: `se NA`, cell `0.000 (NA, NA)`, `subgroup_glm()` `NA NA`;
  - round 1: `se 0`, cell `0.000 (0.000, 0.000)`, `subgroup_glm()` `0 0`;
  - final: `se NA`, `0.000 (NA, NA)`, `NA NA`, the same as start. Walls 2.52 s / 2.59 s.

### B1.8 Gates

| Gate | Result | Evidence |
|---|---|---|
| Red, then green | Red on the parent + test only (5 failures, 4 errors); green 35/35 | `S/main/B1/red.txt`, `green.txt` |
| Independent expected values | Every value by hand or from a separately fitted `glm()` with finite differences; comments beside each; re-derived by both verifiers | B1.4 |
| Unaffected case | Identical in four runs (216, 852, 216 + 852, 1044 results) | B1.5 |
| Scope | 2 files: `R/glm_effect_estimators.R`, the new test | `S/main/B1/diffstat.txt` |
| `document()` | Wrote nothing: 12.85 s (round 1), 12.51 s (verifier 1), 12.39 s (repair), 12.49 s (verifier 2), 12.56 s (polish), 12.46 s (polish check) | `document_check.txt` files |
| ASCII | Added lines ASCII only. The file's two pre-existing non-ASCII lines (`:65`, `:71`, "León") are in start and untouched. | VERIFY r1 §5, r2 §1 |

### B1.9 Verification, repair and polish

**Round 1 verifier: PROBLEMS, one blocking** (`verify/VERIFY.md`):

- **Blocking.** `.fs_rd_gcomp_se()` dropped an aliased **treatment** coefficient like any aliased covariate (fixed round-1 copy `:612`). On a slice whose effect is not estimable (one arm empty, or all of one arm's outcomes NA), `X1 = X0` on the kept columns, so `g = 0`, the helper returned 0 and the closure accepted it: `se` 0 where start gave NA. Effects: cell `0.000 (0.000, 0.000)`; `subgroup_glm()` `c(0, 0)` instead of `c(NA, NA)` (`verify/empty_arm.txt`, `empty_arm2.txt`, `probe.txt`).
  - **Repair:** the design matrices are built with every column, and the helper stops (`tryCatch` → `NA_real_`) when any aliased column differs between `X1` and `X0`. The caller then keeps the old expression, NA there, as at start. New test block `:256-303`, written first; it fails 4 times on the round-1 code. Verifier round 2 confirmed.
- Minor 1: NOTES said a 3+-level factor sends adjusted RD to tier 2 "on every slice" / "always". This is false when the slice holds fewer levels (probe `site_2of3`). NOTES corrected in the repair section (supersedes §8 and §9 item 2). `workflow1_result.json` `impl.open_items[1]` still carries "always".
- Minor 2: test (b)'s tolerance comment and its redundant rounding assertion were edited after the first green run (disclosed; tolerance unchanged). Not changed.
- Minor 3: the call-site comment said "over every coefficient". Fixed in repair ("every estimable coefficient").
- Minor 4: the file-header helper list did not name `.fs_rd_gcomp_se()`. Fixed in polish (`:21`).
- Minor 5: the proposed NEWS bullet was silent on the empty-arm change. Moot after the repair.

**Round 2 verifier: CONFIRMED, no blocking problem, 6 minor** (`verify_r2/VERIFY.md`):

1. The helper roxygen said an aliased column that treatment changes means "the risk difference is not estimable, g would be 0". That is false for a user-supplied aliased treatment-interaction term with the treatment coefficient estimable (`adjust_covariates = c("site", "treat:age")`, `age = 60`: helper NA, closure keeps the old 0.0936; the delta method on the estimable columns gives 0.0884; `probe_interaction_aliased.txt`). **Polish:** roxygen (head `:618-627`), `@return` (`:636-638`) and the inline comment (`:653`) reworded to say what the code does.
2. File-header list → **polish**, `:21`.
3. NOTES "the new block's closure assertions pass on the start tree" is not shown by the red run, which stops at `:294`. It is true by the verifier's start run (`closure_values.txt`: `se` NA for both slices). NOTES not changed.
4. NOTES callers table, bootstrap row, "Replicates use the estimate only": the replicates' inner `forestsearch()` (`bootstrap_analysis_dofuture.R:615`) use the SE through `L(HR)`/`U(HR)`, also exact-tie keys in `remove_redundant_subgroups()` (`subgroup_consistency_helpers.R:1159`). NOTES not changed; recorded here (B1.2) and as an open item.
5. Test (b) headroom about 5× (1.87e-07 against 1e-6). Carried over; no change.
6. NOTES §2 still quotes the round-1 comment and lines. It is superseded by the repair section. Not changed.

**Polish** (main session's minor-finding pass, `NOTES.md` "Polish before commit", `polish.diff`, `polish_result.json`):

- Comment and roxygen only, in `R/glm_effect_estimators.R`: `:21` header list; helper roxygen description; `@return`; the inline comment above the `stopifnot`. No code line, `tests/` file, expected value, seed or tolerance changed.
- Comment gate `git diff -U0 verified -- R tests | ... | grep -vE "^[+-][[:space:]]*#"` printed nothing; `parse()` OK; `document()` wrote nothing (12.56 s); tests 35/0 (B1.3).

**Polish check: PROBLEMS, one minor** (`polish_result.json` `check`; probes `verify_polish/probe_polish.txt`, 2.54 s):

- The new roxygen claim "with age constant on the slice" fails when the constant is 0: `treat:age` is then an all-zero aliased column on which `X1` and `X0` agree, so the helper returns the delta-method SE (P2: 0.0884271104; old expression 0.0935974056). P1 (age = 60), P3/P3b (empty arm), P4 (all-NA treated outcomes) and P5 (never an error) confirmed the other claims.
- **Main-session correction:** "with age constant on the slice" → "with age constant and nonzero on the slice" (head `:625-627`). Comment only. `patches/B1.patch` was regenerated, and `polish.diff` (17:04) carries the corrected text. The gate was re-run, empty (NOTES last section).

### B1.10 Departures and differences from the task's text

- **Safety rule, wider than "not a finite, non-negative number".** The helper also returns `NA_real_`, and the caller then keeps the old expression, in two more cases:
  - on any error or warning while computing the value, so the fix adds no warning and routes no slice to tier 3;
  - when an aliased column differs between `X1` and `X0` (the repair). Without that, a finite 0 would replace start's NA on a slice whose effect is not estimable.

  In the second case the old expression gives NA when the treatment coefficient is aliased, and the old finite value for an aliased treatment-interaction column.
- **"Estimable columns only"** = drop the NA (aliased) coefficients by name from `b`, `X0`, `X1`, `V`.
- **Helper signature.** It takes the fitted probabilities `pred0`/`pred1` the estimate averages, not `plogis(X b)`: `binomial()`'s `linkinv` clamps |eta| > 30, `plogis` does not.
- **Test route.** Deterministic datasets reach the tier through `make_effect_estimator()` (a 3-level factor; an arm at 100%), so the task's fallback "test the calculation directly" was not needed. Test (b) calls the helper directly because the closure's unadjusted branch does not use it.
- **Extra tests** beyond the task's three: the mocked safety rule, the helper never erroring, the not-estimable block.
- **Comment and header edits.** The comment at start `:531-532` was replaced (head `:549-553`), as the task asks. The header list at `:21`, outside `.estimate_rd()` and the helper, was edited in polish.
- **Test-comment wording changed after the first green run** (round 1), with no tolerance or expected value changed:
  - test (b)'s tolerance reason, "far below 1e-6" → the IRLS-lag mechanism with the measured 1.9e-7;
  - `round(., 6) == 0.116548` → `round(., 5) == 0.11655` (the value sits 3.6e-8 from a 6-digit rounding boundary; both pass);
  - "near 1e-9" / "near 1e-7" → "at most about".

  Red was re-run on the final file.
- **Which `g' V g`.** The approved formula is the delta method conditional on the sample's covariates. It leaves out the variance from averaging over the empirical covariate distribution, about `var(p1i - p0i)/n`. Recorded, not changed.
- **Line pointers.** The task says "near 534-547"; the adjusted expression is at `:546` and the comment at `:531-532`, with the quoted text matching.
- **Exploration before the test** (`scratch/explore1-3.R`): it printed the start-tree SE (6.79e-06), a hand-coded `g'Vg` and the (c) closed form, found the factor route to tier 2 and the tier-1 warning error, and checked the factor data's cells. The expected values in the test are by hand or by finite differences.

### B1.11 Open items

1. **Pre-existing, outside B1's scope.** Tier 1's warning handler calls `invokeRestart("muffleWarning")` inside a `tryCatch()` handler, after the stack has unwound (`.estimate_rd()` start `:509`, head `:527`; `.estimate_rr()` start `:686`, head `:775`).
   - Any identity-link or log-binomial warning without "converge" (e.g. "step size truncated: out of bounds") makes the closure stop with "no 'restart' 'muffleWarning' found" instead of falling to the next tier.
   - Examples: adjusted RD with a numeric covariate at n = 200, 400 (`tier1_warning_probe.txt`); unadjusted RR with an arm at 0% or 100%.
   - 20 / 40 / 48 such errors in the three unaffected runs, identical on both trees.
2. Both arms degenerate: adjusted SE about 2.3e-07, unadjusted (F13) exactly 0. Both reflect zero binomial variance at the boundary.
3. Conditional delta method (B1.10).
4. An aliased user-supplied treatment-interaction column with a nonzero constant keeps the old SE (0.0936 against 0.0884 by the delta method on the estimable columns). The value is the same as at start and there is no error; the helper roxygen says so. It lies outside `.build_adjusted_formula()`'s `y ~ treat + covariates`.
5. The SE enters `L(HR)`/`U(HR)`, which are near-duplicate keys (`subgroup_consistency_helpers.R:1120-1131` start) and exact-tie keys in `remove_redundant_subgroups()` (`:1159` start, `:1165` head). It is also used inside each bootstrap replicate's `forestsearch()`. A changed SE can therefore in principle change which candidates are removed as duplicates. In both verifiers' `forestsearch()` runs with an adjusted RD the candidate rows and the selected subgroup were unchanged. NOTES' callers table still says "Replicates use the estimate only".
6. Test (b) has about 5× headroom against its tolerance (measured 1.87e-07, tolerance 1e-6).
7. Evidence wording not corrected in NOTES: §2's round-1 quote, the "closure assertions pass on the start tree" sentence, and the bootstrap row (verifier round 2 minors 3, 4, 6).

### B1.12 NEWS entry (`fb80dd8f:NEWS.md:432-447`, "Behaviour changes")

Each claim checked against the code at `fb80dd8f` and the evidence:

| Claim | Check |
|---|---|
| Applies to `make_effect_estimator(outcome_type = "binary", effect_measure = "RD")` with `adjust_covariates` when `method_used = "logistic_gcomp_adj"` | True (`:566-578`; `method_used` `:587`). |
| Now `sqrt(g' V g)` of the averaged predictions over every estimable coefficient | True (`:640-666`). |
| It was the treatment coefficient's SE scaled by the average derivative | True (start `:533-535`, `:546`: `avg_deriv = mean(pred1 * (1 - pred1))`). |
| "could collapse to about 0 when an arm had 0% or 100% events" | **Incomplete.** The collapse (about 5e-06) happens when the **treated** arm is at 0% or 100%. With the **control** arm at 0% or 100%, the old SE blew up instead: 424.5-451.8 in the implementer's run (`unaffected.txt`) and 387.8-448.5 in verifier round 2, e.g. the cell `0.8044 (-830.7797, 832.3885)`. The implementer's proposed bullet said "collapse towards 0 or reach hundreds"; the committed bullet keeps only the first half. |
| Control 14/40, treated 40/40, one covariate: 0.0775, was 6.8e-06 | True (`example.txt`: 6.79127662230129e-06 → 0.0775100802854232). |
| "giving a zero-width interval" | Inexact. The interval is 0.6125 ± 1.33e-05 (width 2.7e-05). It prints as `0.6125 (0.6125, 0.6125)` at 4 decimals. |
| Previous value kept if the new one is not finite and non-negative, or a column that treatment changes is aliased | True (`:573`, `:654-655`). The helper also gives NA, and the old value is kept, on any error or warning. Read as covered by "not a finite ... number". |
| Estimate, unadjusted RD, identity-link tier and every other measure unchanged | True (B1.5). |
| The SE reaches the candidate table's confidence limits, the bootstrap's observed-estimate interval, and the GLM tables and plots of a fit with `adjust_covariates` | True. **Omits** the exported `subgroup_glm()`, whose returned upper bound is `est + z * se` (`R/run_subgroup_sims.R:371`). Also omits that the candidate `L(HR)`/`U(HR)` are near-duplicate and exact-tie keys (`subgroup_consistency_helpers.R:1127`, `:1165`), inside `forestsearch()` and inside each bootstrap replicate's `forestsearch()`. |
| No default call reaches this tier | True (`forestsearch()` and `make_effect_estimator()` default `adjust_covariates = NULL`; the binary `effect_measure` resolves to `"OR"`). |

## C1 — `format_subgroup_summary_tables()` highlights the `">=0.95"` consistency row

`S` = `/tmp/cran_prep_displays_se_docs_2026-11-12`. Line numbers are at the start SHA `b06de4e8` and at head `fb80dd8f`. `R/` is identical at `835fac1b` and `fb80dd8f`. The two files of this fix do not change between `9ec67b7e` and `fb80dd8f` (`git diff 9ec67b7e fb80dd8f -- <both files>` is empty), so head line numbers equal those at the commit.

### C1.1 Status, commit, provenance

| | |
|---|---|
| Status | **PASS** |
| Commit | `9ec67b7e` "fix: highlight the >=0.95 row of the bootstrap consistency table" |
| Parent | `8e87cf30` (B1) |
| Files (`git diff --numstat 8e87cf30 9ec67b7e`) | `R/format_subgroup_summary_tables.R` +4/-2; `tests/testthat/test-consistency-highlight-labels.R` +121 (new). 2 files, +125/-2. No `man/` change (no roxygen changed). |
| Commit vs gated patch | `git diff 8e87cf30 9ec67b7e` equals `S/patches/C1.patch` once the `index` lines are dropped (checked for this section). `S/main/C1/patch_equal.txt` is empty. |
| Earlier patch | `S/patches/C1.verified.patch` (before polish): same `R/` hunk; test file 116 lines (+120/-2 in all). |
| Origin | Part 4 report (`b06de4e8:dev/reports/REPORT_cran_prep_code_fixes_2026-11-12.md`): F1 verifier minor problem (`:111`) and OPEN ITEMS 7 (`:800`). |

- History (verifier, checked here with `git log -S`): the formatter has had `">= 0.95"` since `eae57f82` (2025-11-04). The summary label was `"≥0.95"` (non-ASCII) until `b440785e` (2025-11-20) made it `">=0.95"`. The two spellings never matched.

### C1.2 Site, before → after

Function `format_subgroup_summary_tables()`, Table 4 "Consistency (P<sub>cons</sub>) Distribution", the green `gt::tab_style()` call.

Before, `b06de4e8:R/format_subgroup_summary_tables.R:215`, `:222`:

```r
        # Color coding: High consistency (green)
...
            rows = `Consistency Range` %in% c("0.85-0.95", ">= 0.95")
```

After, `fb80dd8f:R/format_subgroup_summary_tables.R:215-217`, `:224`:

```r
        # Color coding: High consistency (green).  The labels in the three
        # highlight sets are the Pcons bin labels of
        # summarize_bootstrap_subgroups() (no space in ">=0.95").
...
            rows = `Consistency Range` %in% c("0.85-0.95", ">=0.95")
```

- One code token changes: `">= 0.95"` → `">=0.95"`. One comment line is reworded and two comment lines are added.
- The label the summary produces: `R/summarize_bootstrap_subgroups.R:420` at `b06de4e8`, `:432` at `fb80dd8f` (moved by C4's roxygen; text unchanged): `labels = c("<0.5", "0.5-0.7", "0.7-0.8", "0.8-0.85", "0.85-0.95", ">=0.95")`, with `include.lowest = TRUE, right = FALSE` (`:421-422` / `:433-434`). `table()` of that factor keeps every level, so the distribution always has these six rows.

Every highlight set in the function, checked against what the summary produces (implementer's table, re-checked at both SHAs):

| Site (start → head) | Highlight | Summary produces | Result |
|---|---|---|---|
| Table 4 green, `:222` → `:224` | `c("0.85-0.95", ">= 0.95")` | `"0.85-0.95"`, `">=0.95"` | **mismatch, fixed** |
| Table 4 yellow, `:232` → `:234` | `c("0.7-0.8", "0.8-0.85")` | both | match |
| Table 4 red, `:242` → `:244` | `c("<0.5", "0.5-0.7")` | both | match |
| Table 4 footnote, `:248` → `:250` | "Color coding: Green (>=0.85), Yellow (0.7-0.85), Red (<0.7)" | consistent with the three sets after the fix | match |
| Table 1 basic stats, `:45-50` | four `Metric` labels | `summarize_bootstrap_subgroups.R:195, 227, 259, 287` (start) | match |
| `create_factor_summary_tables()`, `:455` → `:457` | `Position == "M.1"` | `Position = paste0("M.", i)` (`summarize_bootstrap_subgroups.R:342`, start) | match |
| Table 5 original agreement, `:276` → `:278` | `rows = 2:4` (positional) | rows 2-4 are "Successful iterations", "Failed iterations (no subgroup)" and `""` | not a label set; listed, not changed (C1.11) |
| Tables 2, 2C, 2D (`:86`, `:122`, `:159`); `create_factor_summary_tables()` `:522` → `:524` | `Rank == 1`, `Percent >= 20`, whole column, `Percent_of_successful > 50` | numeric or column-wide | not label sets |

- No other `">= 0` label under `R/`: `git grep -n '">= 0' b06de4e8 -- R` finds only `:222`; the same grep at `fb80dd8f` finds nothing.

### C1.3 Red and green runs

Command: `NOT_CRAN=true timeout 300 Rscript --vanilla -e "devtools::test('<tree>', filter = 'consistency-highlight-labels', reporter = 'summary', stop_on_failure = FALSE)"`. The test file has 2 `test_that` blocks and 7 expectations.

| Run | Tree | Result | Failure lines | Wall | File |
|---|---|---|---|---|---|
| Main session, red | repository at `8e87cf30` + test file only | `.123.45`: 2 pass, 5 fail | `:98, :103, :106, :119, :120` | 3.1 s | `S/main/C1/red.txt`, `wall.txt` |
| Main session, green | with the code | `.......`: 7 pass, 0 fail | — | 3.0 s | `S/main/C1/green.txt` |
| Implementer, red | start copy + verified test file | 2 pass, 5 fail | `:93, :98, :101, :114, :115` | 3.27 s | `S/evidence/C1/red.txt` |
| Implementer, green | `S/wt/C1` (verified patch) | 7 pass, 0 fail | — | 3.00 s | `S/evidence/C1/green.txt` |
| Implementer, polish red | `S/wt/C1_polish_red` (start + polished test file) | 2 pass, 5 fail | `:98, :103, :106, :119, :120` | 3.32 s | `S/evidence/C1/polish_red.txt` |
| Implementer, polish green | `S/wt/C1_polish` | 7 pass, 0 fail | — | 3.26 s | `S/evidence/C1/polish_green.txt` |
| Verifier, red | `S/verify/C1/red` (start + verified test file) | 2 pass, 5 fail | `:93, :98, :101, :114, :115` | 3.41 s | `verify/test_red.txt` |
| Verifier, green | `S/verify/C1/fixed` | 7 pass, 0 fail | — | 3.34 s | `verify/test_fixed.txt` |
| Polish check, start control | `S/verify/C1_polish_startctl` (start + polished test file) | 2 pass, 5 fail | `:98, :103, :106, :119, :120` | not recorded | `verify_polish/test_startctl.txt` |
| Polish check, green | `S/verify/C1_polish` | 7 pass, 0 fail | — | 3.36 s | `verify_polish/test_polished.txt` |

- First failure, the same in every red run (head line `:98`, verified-file line `:93`):

  ```
  Expected `rs$fill` to be identical to `c(.c1_red, .c1_red, .c1_yellow, .c1_yellow, .c1_green, .c1_green)`.
  `actual[3:6]`:   "#FFF3CD" "#FFF3CD" "#D4EDDA" ""
  `expected[3:6]`: "#FFF3CD" "#FFF3CD" "#D4EDDA" "#D4EDDA"
  ```

- All five failures are the `">=0.95"` row: no fill (`""`) and bold `FALSE`. The two passes on the start tree are the premises (row labels `:92`; counts `:114`).
- The +5 shift of the failure lines between the verified and polished files comes from the polish (header 5 lines longer).

Verifier mutants (round 1, of the verified copy; `verify/mutations.txt`):

| Mutant | Change | Result | Wall |
|---|---|---|---|
| mut1 | green set `c("0.85-0.95")` | 5 fail (`.123.45`) | 3.26 s |
| mut2 | green set `c("0.85-0.95", ">=0.95", ">= 0.95")` (extra unmatched label) | **7 pass** | 3.29 s |
| mut3 | yellow `"0.8-0.85"` → `"0.80-0.85"` | 1 fail (`:93`) | 3.32 s |

Polish-check mutants (of the polished copy; `verify_polish/test_m1..m4.txt`; walls not recorded):

| Mutant | Change | Result |
|---|---|---|
| m1 | green set gains `">= 0.95"` (= round-1 mut2) | 7 pass |
| m2 | `"0.7-0.8"` added to the green set | 2 fail (`:98`, `:103`; that row `"#D4EDDA+#FFF3CD"`, bold) |
| m3 | `"0.85-0.95"` added to the red set | 1 fail (`:98`; that row `"#D4EDDA+#F8D7DA"`) |
| m4 | `"0.5-0.7"` removed from the red set | 1 fail (`:98`; that row `""`) |

### C1.4 Expected values

Test file `fb80dd8f:tests/testthat/test-consistency-highlight-labels.R` (121 lines). The table is built by the tree's own `summarize_bootstrap_subgroups()` from synthetic rows (`.c1_summary()`, `:51-60`) and formatted by `format_subgroup_summary_tables()`. The test reads which rows of the gt object carry a style from its `_styles` record (`.c1_row_styles()`, `:65-81`). It does not read the source text.

| Line | Value | Derivation | Comment |
|---|---|---|---|
| `:92-95` | labels `"<0.5","0.5-0.7","0.7-0.8","0.8-0.85","0.85-0.95",">=0.95"` | documented bins, left-closed (as in `test-bootstrap-summary-bins.R`); a premise, passes on both trees | header `:14-18`, `:91` |
| `:87` (input) | Pcons 0.30, 0.60, 0.75, 0.82, 0.90, 0.97 → one value per bin | by hand, a <= v < b | `:84-86` |
| `:98-101` | fills red, red, yellow, yellow, green, green | the table's own footnote "Green (>=0.85), Yellow (0.7-0.85), Red (<0.7)" applied to each bin's range | header `:19-27`, `:97` |
| `:44-46` | colours `#D4EDDA`, `#FFF3CD`, `#F8D7DA` | the three `cell_fill()` literals of the green / yellow / red calls (not touched by the fix), upper-cased because gt stores them so (seen on the start tree in the red run) | header `:28-32` |
| `:103` | bold `FALSE FALSE FALSE FALSE TRUE TRUE` | only the green call adds `cell_text(weight = "bold")` | header `:32`, `:102` |
| `:106` | `">=0.95"` fill = green | the row the task names | `:105` |
| `:114` | Count `c(0L,0L,0L,0L,0L,5L)` for Pcons 0.95, 0.96, 0.98, 0.99, 1.00 | by hand: all >= 0.95; 1.00 kept by the closed top end | `:110-112` |
| `:119`, `:120` | `">=0.95"` fill green, bold `TRUE` | as `:98-103` | `:118` |

- No tolerance: every check is `expect_identical()` or `expect_true()`.
- No value was taken from the fixed code's output (header `:42`). The verifier re-derived every value with base R only, no package code (`verify/scripts/v_hand.R`, `verify/hand_derivation.txt`, wall 0.15 s): bins, counts `0 0 0 0 0 5`, colour per row, bold per row and the upper-cased fills all equal the expectations.
- What the tests check and do not check is stated in the header after the polish (`:34-41`): a row with no fill or a second fill fails; a label in a highlight set that matches no summary row is not seen.

### C1.5 Unaffected case

The task names no unaffected case for C1. The implementer's case: everything in the formatted output other than the highlight on the `">=0.95"` row is identical, start vs fixed.

Method (implementer and verifier alike):

- Each tree in its own `Rscript --vanilla` process; tables built with that tree's `summarize_bootstrap_subgroups()` and `format_subgroup_summary_tables()`; compared in a third process.
- Summaries: `identical()`.
- Every gt table: (a) `identical(x, y, ignore.environment = TRUE)`; (b) `identical()` of the rendered HTML (`gt::as_raw_html(inline_css = FALSE)`, `set.seed(1)` before each render because gt draws a random table id).
- `consistency_dist`: the fixed object's `_styles` rows on the `">=0.95"` row are removed (row names reset on both) before (a) and (b).
- Leaf relaxed: closure environments inside gt's `_formats` (the `fmt_number()` formatters) are new in every process, so plain `identical()` is FALSE even start vs start. (a) ignores them; (b) re-renders, which runs the formatters, so the formatted cell text is still compared exactly.

Implementer (`S/evidence/C1/scripts/unaffected_run.R`, `unaffected_compare.R`; outputs in `S/evidence/C1/unaffected/`):

- Four inputs: `per_bin` (one Pcons per bin), `all_top` (all Pcons >= 0.95), `none_top` (`">=0.95"` Count 0), `full` (60 synthetic rows, seed 20261112, 80% found, with `original_sg` and GRF cuts, so all nine tables are built).
- Result: `ALL IDENTICAL (outside the >=0.95 row's highlight): TRUE`. Summaries identical in all four cases; all 30 tables (7 + 7 + 7 + 9) identical in structure and HTML.
- `consistency_dist`, every case: the only differing component before removal is `_styles`; the removed rows are exactly 3 (cells `Consistency Range` / `Count` / `Percent` of row 6, fill `#D4EDDA`, weight `bold`); the start tree has 0 style rows there.
- HTML diff of the consistency table, all four cases: lines 480-482 only (`html_diff_start_vs_fixed.txt`).
- Control, start vs start (`start_a.rds` vs `start_b.rds`): all TRUE (`control_start_vs_start.txt`).
- Walls (`/usr/bin/time -f %e`): run start 3.17 s, run fixed 3.01 s, compare 3.79 s. Control compare 3.88 s (NOTES §5; no file records it).

Verifier (`verify/scripts/v_run.R`, `v_compare.R`, `v_control.R`; outputs in `S/evidence/C1/verify/unaffected/`):

- Eight other inputs: `on_breaks` (Pcons exactly 0, .5, .7, .8, .85, .95, 1), `below_breaks` (.4999 ... .9999), `na_and_failed` (NA Pcons and `any_found = 0` rows), `out_of_range` (1.02, -0.1, 0.96), `single_top` (nb = 1), `no_pcons` (no consistency table), `none_found` (early return), `big` (120 rows, seed 7, `original_sg`, 30 GRF cuts; all nine tables).
- Result: `ALL IDENTICAL outside the >=0.95 highlight: TRUE`. In the six cases with a consistency table, start fills end `>=0.95=none` and fixed fills end `>=0.95=#D4EDDA`; 3 fixed-only style rows, all green + bold; HTML differs on lines 480-482 only. `no_pcons` and `none_found` identical; no error on either tree.
- Control (two start runs): plain `identical()` FALSE on 42 of 50 tables, only `_formats` differing; TRUE on the 8 without `fmt_number()`; with `ignore.environment = TRUE`, 0 FALSE.
- Walls: run start 3.53 s, run fixed 3.32 s, compare 5.24 s; second start run 3.48 s, control compare 1.83 s.

### C1.6 Existing tests

- No existing test calls `format_subgroup_summary_tables()` or `create_factor_summary_tables()` (`git grep` over `b06de4e8:tests` finds none). No package function calls the formatter.
- The test of the label producer was run anyway:
  - implementer, fixed copy: `bootstrap-summary-bins` 8/8, `consistency-highlight-labels` 7/7, 0 fail, 3.33 s (`S/evidence/C1/existing_tests.txt`);
  - verifier, fixed copy: the same counts, 3.40 s (`verify/existing_tests.txt`); start copy: `bootstrap-summary-bins` 8/8, 3.12 s (`verify/existing_tests_start.txt`).
- Full suite at `fb80dd8f` (log still being written when this section was drafted): `✔ |          8 | bootstrap-summary-bins` (`S/suite/suite.log:90`) and `✔ |          7 | consistency-highlight-labels` (`:105`).
- Expectations updated: none.

### C1.7 What a user sees

`S/evidence/C1/scripts/example.R`: seven bootstraps with Pcons 0.95, 0.96, 0.98, 0.99, 1.00, 0.90, 0.88; `S/evidence/C1/example_before_after.txt`. The fixed tree is `S/wt/C1`, whose `R/` equals the commit's.

```
tree: /tmp/cran_prep_displays_se_docs_2026-11-12/start
 Consistency Range Count Percent    highlight
              <0.5     0    0.00      #F8D7DA
           0.5-0.7     0    0.00      #F8D7DA
           0.7-0.8     0    0.00      #FFF3CD
          0.8-0.85     0    0.00      #FFF3CD
         0.85-0.95     2   28.57 #D4EDDA bold
            >=0.95     5   71.43       (none)
tree: /tmp/cran_prep_displays_se_docs_2026-11-12/wt/C1
 Consistency Range Count Percent    highlight
              <0.5     0    0.00      #F8D7DA
           0.5-0.7     0    0.00      #F8D7DA
           0.7-0.8     0    0.00      #FFF3CD
          0.8-0.85     0    0.00      #FFF3CD
         0.85-0.95     2   28.57 #D4EDDA bold
            >=0.95     5   71.43 #D4EDDA bold
```

- In the rendered table the row holding most bootstraps is now green and bold, as the footnote "Green (>=0.85)" says. Before, it was the only row without colour.
- Rendered HTML of that row (`all_top` case, `html_diff_start_vs_fixed.txt`), start then fixed:

  ```
  <tr><td headers="Consistency-Range" class="gt_row gt_left">&gt;=0.95</td>
  <tr><td headers="Consistency-Range" class="gt_row gt_left" style="background-color: #D4EDDA; font-weight: bold;">&gt;=0.95</td>
  ```

### C1.8 Reach and articles

- `format_subgroup_summary_tables()` is internal (`@keywords internal`, `fb80dd8f:R/format_subgroup_summary_tables.R:10`), not in `NAMESPACE`, and has no example (`man/format_subgroup_summary_tables.Rd`).
- No call under `R/`, `tests/` (other than the new test), `inst/` or `vignettes/` at either SHA. `summarize_bootstrap_results()` returns the unformatted summary.
- The only calls are `quarto/resampling/test_bootstrap_grf_cuts.qmd:314, 401` (`forestsearch:::`), read-only, unchanged and not rendered.
- No article evaluates a call that reaches the changed line, so no article output can change.

### C1.9 Verification, repair, polish

- **Round 1 verdict: CONFIRMED** (`S/evidence/C1/verify/VERIFY.md`). No blocking problem. No repair round (`workflow1_result.json`: `"repaired": false`, `blocking_r1` and `blocking_final` empty).
- The verifier re-ran red and green, built three mutants, re-derived every value by hand, ran an eight-case unaffected comparison with a control, checked scope (2 files; `R/` changes are 3 comment lines and one `%in%` vector), ASCII, callers, and that `gt` and `data.table` are in `Imports`.

Minor problems and what happened to each:

| # | Problem (VERIFY.md) | Outcome |
|---|---|---|
| 1 | Test header (verified file `:33-36`) concluded that every highlight label is a summary label, from "two labels each (six in all)", a fact about the source the tests never check; mut2 (extra `">= 0.95"`) passes 7/7. | **Fixed in polish**: the header now says what the tests check and what they do not (head `:33-41`). |
| 2 | NOTES §2 says the CV binning "has no highlight"; `.format_cv_pconsistency_gt()` applies a `tab_style()` on `rows = is_summary` (`b06de4e8:R/summarize_cv_results.R:1174-1180`), never on a bin label. | Not in the polish list; NOTES not changed. No code or test effect. Exact statement: the CV Pcons table (labels `:761-762`, the same six including `">=0.95"`) has no bin-label highlight. |
| 3 | The test reads gt's internal `_styles` record (gt 1.4.0); a gt that restructures it would break the test with no change in behaviour. | Left as is; recorded in NOTES §8 and in C1.11. |

Polish (main session; `S/evidence/C1/NOTES.md`, "Polish before commit"; `polish_result.json` id C1):

- What changed: `tests/testthat/test-consistency-highlight-labels.R` header only, verified `:33-36` → polished `:33-41` (one hunk, -4/+9, all `#` lines; `S/evidence/C1/polish.diff`). `R/format_subgroup_summary_tables.R` is the same as in the verified patch.
- Old:

  ```
  #   - So every one of the six rows carries exactly one fill.  Since the
  #     three highlight sets name two labels each (six in all) and each of
  #     the six summary labels is highlighted exactly once, every label in
  #     the highlight sets is one of the labels the summary produces.
  ```

- New:

  ```
  #   - So every one of the six rows carries exactly one fill.
  # What the tests check: which rows of the formatted table get each fill
  # and which get bold text.  In the first test every row's fills must equal
  # the one above (the second checks the ">=0.95" row only), so a row with
  # no fill or with a second fill fails; a summary label missing from its
  # highlight set (as ">=0.95" was) shows as a row without its fill.
  # What they do not check: a label in a highlight set that matches no
  # summary row (say an extra ">= 0.95" next to ">=0.95") styles no row, so
  # it changes nothing these tests can see, and they still pass.
  ```

- Why: verifier minor 1; the old conclusion rested on the source, not on what the tests see.
- Gate: comment-only grep against the verified copy prints nothing; ASCII only; added lines <= 80 characters. Red and green re-run (C1.3): same counts and failures, lines +5.
- **Polish check: CONFIRMED**, no problems (`polish_result.json`, `check`). It confirmed each new sentence with a probe (`verify_polish/probe.R`) and the four mutants m1-m4 (C1.3), and the polished test passes 7/7 (3.36 s).

### C1.10 Departures

- **Unaffected case.** The task names none for C1. The implementer defined one (C1.5).
- **Direction of the test.** The task asks: "every label in the highlight set is one of the labels of the consistency distribution". The test checks the converse: every distribution row carries exactly the fill its range calls for, so every distribution label is in the right highlight set. A label in a highlight set that matches no row is not seen (round-1 mut2, polish m1). The header says so (`:39-41`). The task's property holds in the code at head by reading: the three sets name `"0.85-0.95"`, `">=0.95"`, `"0.7-0.8"`, `"0.8-0.85"`, `"<0.5"`, `"0.5-0.7"` (`:224, :234, :244`), all six labels at `R/summarize_bootstrap_subgroups.R:432`. No test checks it.
- **Comment lines.** Three comment lines above the call (one reworded, two new) say where the labels come from. No roxygen change.
- **gt internals.** The test reads `_styles` rather than rendered HTML, to stay independent of gt's markup; HTML is used in the unaffected case and the example.
- **Commit subject.** The implementer proposed "... bootstrap consistency distribution table"; the commit says "... bootstrap consistency table".
- The task's site "near 222" is exactly `b06de4e8:R/format_subgroup_summary_tables.R:222`; the code matches the task's description.

### C1.11 Open items

- **Table 5 (`original_agreement`) highlights by position.** `rows = 2:4` (`b06de4e8:R/format_subgroup_summary_tables.R:276`, `fb80dd8f:...:278`) styles "Successful iterations", "Failed iterations (no subgroup)" and the blank separator row of the table built at `R/summarize_bootstrap_subgroups.R:548-557` (start) / `:560-569` (head). The match rows ("Exact match with original", "Different from original", "Partial match (shared factor)") are rows 6-8. Not a label mismatch of C1's kind; listed for a decision, not changed.
- **The test depends on gt's internal `_styles`** (gt 1.4.0). A gt release that restructures it would break the test with no change in behaviour.
- **The test cannot see an extra unmatched highlight label** (header `:39-41`).

### C1.12 NEWS entry

`fb80dd8f:NEWS.md:449-451`, under "### Behaviour changes" (`:154`):

> * The internal `format_subgroup_summary_tables()` now highlights the `">=0.95"` row of the bootstrap consistency distribution, as its footnote says. The label was written `">= 0.95"`, which matched no row.

| Claim | Checked against | Holds |
|---|---|---|
| "internal" | `@keywords internal` (`:10`); not in `NAMESPACE` at either SHA | yes |
| "now highlights the `">=0.95"` row" | head `:224`; red/green runs; example (green fill `#D4EDDA`, bold) | yes |
| "of the bootstrap consistency distribution" | Table 4 "Consistency (P<sub>cons</sub>) Distribution" from `summarize_bootstrap_subgroups()`'s `consistency_dist` (head `:191-194`) | yes |
| "as its footnote says" | footnote "Green (>=0.85)" (head `:250`); every value in `">=0.95"` is >= 0.85 | yes |
| "The label was written `">= 0.95"`, which matched no row" | start `:222`; summary labels at start `:420` contain no `">= 0.95"` | yes |

- No inexact claim. The entry does not say "and bold"; that is an omission, not an error.

### C1.13 Evidence notes

- NOTES §2 (last paragraph) and `workflow1_result.json` `impl.open_items[1]` say the CV binning "has no highlight". The code applies a `tab_style()` (fill `#f0f0f0`, bold) on `rows = is_summary` (`b06de4e8:R/summarize_cv_results.R:1174-1180`, `fb80dd8f:...:1179`). The exact statement is "no bin-label highlight" (verifier minor 2; NOTES not corrected).
- VERIFY.md §3 places the verified file's expectations at "lines 85-115" and comments at ":79-81, :84, :91, :97, :106-108, :112". In the verified test file (rebuilt from `S/patches/C1.verified.patch`) the expectations are at `:87-115` and the comments at `:79-81, :86, :92, :97, :100, :105-107, :113`. The values are not affected.
- NOTES §3-§4 and VERIFY.md §2, §5 give line numbers (`:93 ...`) and counts (116 lines, +120/-2) for the verified file. The polish section and the commit have the polished file (121 lines, +125/-2, failures at `:98 ...`). This is superseded, not contradicted.
- The polish check has no `VERIFY.md` in `S/evidence/C1/verify_polish/`; its verdict is in `polish_result.json`. Its mutants are numbered m1-m4 independently of round 1's mut1-mut3 (polish m1 = round-1 mut2).
- The implementer's control compare wall (3.88 s) is only in NOTES §5; no output file records it.

## C2 — Threshold warnings name the measure and the default that apply (rate difference)

**Status:** PASS. Commit `2984dbee` ("fix: name the rate difference, IRR scale and 0.01 default in the IRD threshold warnings"), parent `b06de4e8` (start SHA).

| file (`git show --stat 2984dbee`, `--numstat`) | change |
|---|---|
| `R/forestsearch_main.R` | `29 +++--` (numstat: 22 insertions, 7 deletions) |
| `tests/testthat/test-threshold-remap-messages.R` | new, 145 lines |
| total | 2 files, 167 insertions, 7 deletions |

- Scratch root `S=/tmp/cran_prep_displays_se_docs_2026-11-12`. Evidence: `$S/evidence/C2/` (implementer), `$S/evidence/C2/verify/` (verifier), `$S/main/C2/` (main session), `$S/patches/C2.patch`.
- No `man/`, `NAMESPACE`, `DESCRIPTION` or `NEWS.md` change in this commit. No roxygen line is touched.
- The committed diff equals `patches/C2.patch`. `main/C2/patch_equal.txt` is empty (0 bytes). I checked it again: `git diff b06de4e8 2984dbee` and the patch differ only in the length of the abbreviated `index` hashes (patch `9738114..c03a8b9` and new file `3dfaee3`; commit diff `97381140..c03a8b98` and `3dfaee3a`). The blobs are the same.
- There is only one patch for C2. No `C2.verified.patch` or `C2.r1.patch` exists, because there was no repair and no polish.
- **After the commit:**
  - `R/forestsearch_main.R` changed once more, in `ba7c4ea2` (C4). Its hunks start at lines 447, 1080, 1318, 2507 and 3286 of `2984dbee`. None of them is inside a C2 region.
  - C4 adds 30 lines above the remap block, so the remap moves from `:2341-2365` to `:2371-2395`. The validator does not move.
  - The test file is unchanged from `2984dbee` to `fb80dd8f`: only `2984dbee` touches it.
  - `R/` is identical at `835fac1b` and `fb80dd8f` (`git diff --quiet` exits 0).

### Site, before and after

| site | function | start `b06de4e8` | `2984dbee` | head `fb80dd8f` |
|---|---|---|---|---|
| validator, check 4 (changed) | `.validate_outcome_threshold_config()` | `R/forestsearch_main.R:301-310` | `:301-317` | `:301-317` |
| warned remap (changed) | `forestsearch()`, SECTION 2B | `:2334-2350` | `:2341-2365` | `:2371-2395` |
| silent default remap (unchanged) | `forestsearch()` | `:2325` | `:2332` | `:2362` |
| consistency remap message (unchanged) | `forestsearch()` | `:2361` | `:2376` | `:2406` |
| comments "RD/IRD are probability/rate differences bounded in [-1, 1]" (unchanged) | both | `:296`, `:2332` | `:296`, `:2339` | `:296`, `:2369` |

- The diff has three hunk headers, `@@ -299,13 +299,20 @@`, `@@ -2332,6 +2339,14 @@` and `@@ -2341,12 +2356,12 @@`, in two regions.
- The task's pointers ("near 301–309" and "near 2334–2350") match the code at the start SHA.

**Validator.** Before, `b06de4e8`, `R/forestsearch_main.R:301-310`:

```r
  if (is_bounded_identity && user_set_threshold && hr.threshold > 1.0) {
    warning(
      sprintf(paste0(
        "effect.threshold = %.2f but effect_measure = '%s' ",
        "(identity scale, range [-1, 1]).\n",
        "  A risk difference cannot exceed 1.0.\n",
        "  Did you mean effect_measure = 'OR' (ratio scale)?"),
        hr.threshold, effect_measure),
      call. = FALSE)
  }
```

After, `fb80dd8f`, `R/forestsearch_main.R:301-317`:

```r
  if (is_bounded_identity && user_set_threshold && hr.threshold > 1.0) {
    # Name the measure and its ratio-scale counterpart: RD -> OR, IRD -> IRR
    # (the ratio measure of the same rate model).  A rate difference is not
    # bounded by 1 in itself; [-1, 1] is the range the threshold remap uses.
    bound_line <- switch(effect_measure,
      RD  = "  A risk difference cannot exceed 1.0.\n",
      IRD = "  A rate-difference threshold above 1.0 is outside that range.\n")
    ratio_measure <- switch(effect_measure, RD = "OR", IRD = "IRR")
    warning(
      sprintf(paste0(
        "effect.threshold = %.2f but effect_measure = '%s' ",
        "(identity scale, range [-1, 1]).\n",
        "%s",
        "  Did you mean effect_measure = '%s' (ratio scale)?"),
        hr.threshold, effect_measure, bound_line, ratio_measure),
      call. = FALSE)
  }
```

**Warned remap.** Before, `b06de4e8`, `R/forestsearch_main.R:2334-2350`:

```r
      if (user_set_threshold && effect_threshold > 1.0) {
        warning(
          sprintf(
            paste0(
              "effect.threshold = %.2f appears to be on a ratio scale, ",
              "but effect_measure = '%s' uses the identity scale ",
              "(range [-1, 1]).\n",
              "  Did you mean effect.threshold = %.2f on the %s scale?\n",
              "  Remapping to default: %.2f.\n",
              "  To use a custom %s threshold, pass a value in [-1, 1] ",
              "(e.g., effect.threshold = 0.07 for a 7 pct-point threshold)."),
            effect_threshold, effect_measure,
            effect_threshold, "OR", 0.05, effect_measure),
          call. = FALSE)
        effect_threshold <- switch(effect_measure,
          RD = 0.05, IRD = 0.01)
      }
```

After, `fb80dd8f`, `R/forestsearch_main.R:2371-2395` (the four unchanged format lines `:2383-2388` are shortened to `...`):

```r
      if (user_set_threshold && effect_threshold > 1.0) {
        # The measure's ratio-scale counterpart, its default and an example:
        # RD -> OR, 0.05; IRD -> IRR, 0.01 (?forestsearch, "Resolved defaults,
        # per estimand").  The message prints the value the remap assigns.
        .fs_remap_ratio   <- switch(effect_measure, RD = "OR", IRD = "IRR")
        .fs_remap_default <- switch(effect_measure, RD = 0.05, IRD = 0.01)
        .fs_remap_example <- switch(effect_measure,
          RD  = "effect.threshold = 0.07 for a 7 pct-point threshold",
          IRD = "effect.threshold = 0.02 for a rate difference of 0.02")
        warning(
          sprintf(
            paste0(
              ...
              "  To use a custom %s threshold, pass a value in [-1, 1] ",
              "(e.g., %s)."),
            effect_threshold, effect_measure,
            effect_threshold, .fs_remap_ratio, .fs_remap_default,
            effect_measure, .fs_remap_example),
          call. = FALSE)
        effect_threshold <- .fs_remap_default
      }
```

- **The defect at the start.** The remap message passed the literals `"OR"` and `0.05` (start `:2346`), while the assignment used `switch(effect_measure, RD = 0.05, IRD = 0.01)` (`:2348-2349`). For IRD, the message said 0.05 while 0.01 was used. The validator said "A risk difference cannot exceed 1.0." and suggested `'OR'` for both measures.
- **The fix.** The value assigned is the same `switch()` as before. It is now computed once, so the message prints the value that is assigned.
- **Every new `switch()` returns a value.** The validator switches run only when `is_bounded_identity` is true (`effect_measure %in% c("RD", "IRD")`, start `:298-299`). The remap switches sit inside `if (effect_measure %in% c("RD", "IRD"))` (start `:2321`, head `:2358`).
- **The new locals do not reach the fit.** They are dot-prefixed, like `.fs_c1_legacy`. `args_call_all` holds formals only: `args_names <- names(formals())` (start `:1923`), then `args_call_all <- mget(args_names, envir = environment())` (`:1924`).
- **Callers.** The validator has one caller, `forestsearch()` at start `:2224`. It receives the measure through `effect_measure = if (exists("effect_measure")) effect_measure else NULL` (`:2229`).

### Red and green

Filter `threshold-remap-messages`; it matches only the new file. Command, all six runs: `NOT_CRAN=true timeout 300 Rscript --vanilla -e "devtools::test('<tree>', filter = 'threshold-remap-messages', reporter = 'summary', stop_on_failure = FALSE)"`.

| run | where | result | wall |
|---|---|---|---|
| implementer red | `wt/C2`, a copy of the start tree plus the test file only (`git diff start -- R/` empty); `evidence/C2/red.txt` | `...1234.56789..............`: 27 expectations, 18 pass, 9 fail | 8.18 s |
| implementer green | `wt/C2` with the patch; `evidence/C2/green.txt` | 27 pass, 0 fail, 0 skip | 7.90 s |
| verifier red | `verify/C2/red`, a fresh copy plus `git apply --include='tests/*'`; `evidence/C2/verify/red.txt` | 18 pass, 9 fail, the same failures | 8.44 s (`red_wall.txt`) |
| verifier green | `verify/C2/fixed`, a fresh copy with the whole patch; `evidence/C2/verify/fixed.txt` | 27 pass, 0 fail | 8.12 s (`fixed_wall.txt`) |
| main-session red | repository at `b06de4e8` plus the test file only; `main/C2/red.txt` | `...1234.56789..............`: 18 pass, 9 fail | 8.1 s (`main/C2/wall.txt`) |
| main-session green | repository with the code (`2984dbee`); `main/C2/green.txt` | `...........................`: 27 pass, 0 fail | 7.7 s |

- **First failure, every red run:** `test-threshold-remap-messages.R:84:3`, Expected `val` to match regexp "rate[- ]difference". The actual text was:
  ```
  effect.threshold = 1.50 but effect_measure = 'IRD' (identity scale, range [-1, 1]).
    A risk difference cannot exceed 1.0.
    Did you mean effect_measure = 'OR' (ratio scale)?
  ```
- **The nine failures, all in the IRD block:**
  - Validator: `:84` (no "rate difference"), `:85` (no `effect_measure = 'IRR'`), `:86` (says "risk difference"), `:87` (says `'OR'`).
  - Remap: `:92` (no "on the IRR scale"), `:93` (no "Remapping to default: 0.01."), `:94` (says "OR scale"), `:95` (prints "0.05"), `:97` (has "pct-point").
- **The 18 that pass on red:**
  - Six in the IRD block: the two `expect_length` checks (`:78`, `:79`), the two "effect_measure = 'IRD'" matches (`:83`, `:91`), and the two recorded values, 0.01 (`:101`, `:102`). The start tree already assigned 0.01; only the text was wrong.
  - The RD block (6) and the in-range block (6).
- **Identical outputs.** `main/C2/red.txt` and `green.txt` are byte-identical to the verifier's `red.txt` and `fixed.txt` (`cmp`). The summary reporter prints no path and no time, so this is expected. The walls are recorded separately and differ.

### Expected values

All expected values are in `tests/testthat/test-threshold-remap-messages.R`. The file is the same at `2984dbee` and `fb80dd8f`, so the line numbers hold at head.

| expected value (test line) | how it was derived | comment |
|---|---|---|
| IRD validator names the measure: `effect_measure = 'IRD'` (`:83`) | By hand: `'%s'` prints the supplied measure. It passes on red too. | `:81-82` |
| IRD validator speaks of a rate difference: matches `"rate[- ]difference"` (`:84`), has no "risk difference" (`:86`) | The measure's name; the task's defect statement. | `:81-82` |
| IRD validator suggests IRR: `effect_measure = 'IRR'` (`:85`), has no `'OR'` (`:87`) | Documentation: `@param effect_measure` lists IRR with IRD, for count (IRR the default) and for binary (start `R/forestsearch_main.R:1089-1093`, head `:1111-1115`). | `:81-82`; header `:13-16` |
| IRD remap names the measure and the IRR scale: `effect_measure = 'IRD'` (`:91`), "on the IRR scale" (`:92`), no "OR scale" (`:94`) | The same documentation; the format line "on the %s scale" is unchanged. | header `:13-16` |
| IRD remap prints the IRD default: "Remapping to default: 0.01." (`:93`), no "0.05" (`:95`) | Documentation: `?forestsearch`, "Resolved defaults, per estimand", IRD c1 = 0.01 (`man/forestsearch.Rd:1208` at `b06de4e8`, `:1223` at `fb80dd8f`). `%.2f` prints "0.01". | `:89-90`; header `:11-12` |
| IRD remap has no percentage-point example: no "pct-point" (`:97`) | That example is in risk-difference units. | `:96` |
| IRD recorded values: `threshold_config$screening` 0.01 (`:101`), `args_call_all$effect.threshold` 0.01 (`:102`) | Documentation: IRD c1 = 0.01. The scale is identity, so the screening value is the threshold itself. Both pass on red. | `:99-100` |
| one validator message and one remap message per fit (`:78-79`, `:109-110`) | By hand: one validator call (start `:2224`) and one remap block per fit. | header `:4-8` |
| RD validator, full text (`:117-121`) | By hand from the start format strings (start `:302-308`): 1.5 printed by `%.2f` is "1.50", plus "OR". | `:112-116` |
| RD remap, full text (`:122-129`) | By hand from the start format strings (start `:2335-2347`): "1.50", "OR" and the RD default 0.05 (`man/forestsearch.Rd:1207` at start, `:1222` at head). | `:112-116` |
| RD recorded values 0.05 (`:132`, `:133`) | Documentation: RD c1 = 0.05. | `:131` |
| in range, RD 0.07 and IRD 0.02 (loop `:139`): no validator message (`:141`), no remap message (`:142`), screening equals the supplied value (`:143`) | By hand from the code's condition `> 1.0` (start `:301`, `:2334`). | `:137-138` |

- **Count:** 27 expectations. The IRD block has 15, the RD block 6, and the in-range block 6 (3 per measure).
- **No tolerance.** Every check is a string match or `identical()` on doubles.
- **Seeds.** 11 (`:22`) and 12 (`:30`) are in the data builders. NOTES §4 says they were set before the first run and never changed.
- **The verifier's re-derivation** (`evidence/C2/verify/v_derive.R`, `v_derive.txt`), in base R with no package code:
  - It parsed RD 0.05 (Rd `:1207`) and IRD 0.01 (Rd `:1208`) from the start Rd: `identical(c1_ird, 0.01): TRUE`, `identical(c1_rd, 0.05): TRUE`.
  - `sprintf('%.2f', 1.5) = 1.50`.
  - It rebuilt the two RD texts by hand from the start format strings. Both equal the test's expressions: `test RD validator text == hand-built: TRUE`, `test RD remap text == hand-built: TRUE`.
- **Independence.** No expected value depends on the fixed code's output. The IRD expectations are fragments built from the unchanged format lines plus the documented IRR and 0.01. The RD texts are what the start tree prints: the RD block passes on red.

### Unaffected case

- **The task's case:** "The remapped values do not change." The implementer and the verifier also checked that the RD text is unchanged and that no other warning or message changes.
- **Implementer** (`evidence/C2/unaffected/run_case.R`, `compare.R`, `compare.txt`):
  - Trees: the start export `$S/start` and `wt/C2`. Each ran in its own `Rscript --vanilla` process and saved with `saveRDS()`. The results were compared with `identical()` in a third process.
  - Data: synthetic, N = 200. Binary with `set.seed(11)`; count with an `ftime` offset and `set.seed(12)`.
  - Fit settings: `forestsearch()` with `maxk = 2`, `fs.splits = 50`, `n.min = 40`, `use_grf = FALSE`, `use_lasso = FALSE`, sequential plan.
  - Timing leaves dropped (wall-clock values only): `$find.grps$time_search` and `$minutes_all`. A recursive leaf diff showed these were the only differing leaves.

  | case | fit identical (timing leaves dropped) | `threshold_config` identical | screening start / fixed | warnings identical | C2 messages start / fixed | messages identical |
  |---|---|---|---|---|---|---|
  | RD 1.5 | TRUE | TRUE | 0.05 / 0.05 | TRUE | 2 / 2 | TRUE |
  | RD 2 | TRUE | TRUE | 0.05 / 0.05 | TRUE | 2 / 2 | TRUE |
  | IRD 1.5 | TRUE | TRUE | 0.01 / 0.01 | FALSE (the two C2 messages only; other warnings identical TRUE) | 2 / 2 | TRUE |
  | IRD 3 | TRUE | TRUE | 0.01 / 0.01 | FALSE (the two C2 messages only; other warnings identical TRUE) | 2 / 2 | TRUE |
  | RD 0.07 | TRUE | TRUE | 0.07 / 0.07 | TRUE | 0 / 0 | TRUE |
  | IRD 0.02 | TRUE | TRUE | 0.02 / 0.02 | TRUE | 0 / 0 | TRUE |

  - `args_call_all$effect.threshold` equals the screening value on both trees in every case.
  - Direct calls to `.validate_outcome_threshold_config()` are identical for RD 1.5, RD 0.07 and IRD 0.02 (`TRUE`). IRD 1.5 differs, as intended.
  - Subgroups found: RD `{x1 <= 0.7}`; IRD `{x1 <= 0} & {x2 <= 5}`.
  - Walls: start process 23.12 s, fixed process 23.16 s, compare 0.19 s (`run_start.txt`, `run_fixed.txt`, `compare.txt`).
- **Verifier** (`evidence/C2/verify/unaffected/v_run.R`, `v_compare.R`, `v_compare.txt`):
  - Its own data: seeds 21, 22 and 23, N = 240, `fs.splits = 40`, `maxk = 2` (`maxk = 1` for the survival cases).
  - Trees: `verify/C2/startctl` (a fresh start copy) and `verify/C2/fixed`, each in its own process, compared in a third.
  - The leaf diff ran before any drop. In every case the only differing leaves were the two timing leaves. In the two DINA cases the result object had no differing leaf, because both trees stop with the same error ("subgroup_method = "dina" does not support identity-scale estimands (RD, IRD)").

  | case | fit identical | screening start / fixed | warnings | C2 messages |
  |---|---|---|---|---|
  | RD 1.5, RD 1.0001, RD `hr.threshold = 1.5`, RD 1.5 + `consistency.threshold = 1.2` | TRUE | 0.05 / 0.05 | identical | 2 / 2 |
  | RD 1.0, RD 0.07, RD -0.5 | TRUE | 1 / 1, 0.07 / 0.07, -0.5 / -0.5 | identical | 0 / 0 |
  | RD unset (silent remap) | TRUE | 0.05 / 0.05 | identical | 0 / 0 |
  | RD 1.5, `subgroup_method = "dina"` | same error on both trees | — | identical | 2 / 2 |
  | IRD 1.5, IRD 3, IRD `hr.threshold = 1.5` (count) | TRUE | 0.01 / 0.01 | only the 2 C2 messages differ | 2 / 2 |
  | IRD 1.5 on a binary outcome with an offset | TRUE | 0.01 / 0.01 | only the 2 C2 messages differ | 2 / 2 |
  | IRD 1.0, IRD 0.02 | TRUE | 1 / 1, 0.02 / 0.02 | identical | 0 / 0 |
  | IRD unset (silent remap) | TRUE | 0.01 / 0.01 | identical | 0 / 0 |
  | IRD 1.5, `subgroup_method = "dina"` | same error on both trees | — | only the 2 C2 messages differ | 2 / 2 |
  | survival + `effect_measure = "IRD"`, 1.5 | TRUE | 0.4054651 / 0.4054651 | only the validator message differs | 1 / 1 |
  | survival + `effect_measure = "RD"`, 1.5 | TRUE | 0.4054651 / 0.4054651 | identical | 1 / 1 |

  - In every case `args_call_all$effect.threshold` is identical, and so are the non-C2 warnings and the messages.
  - **Validator grid** (direct calls): 7 measures (RD, IRD, OR, RR, IRR, MD, NULL), 4 outcome types, 6 thresholds (1.5, 1.0001, 1.0, 0.3, 0.07, -2) and `user_set_threshold` TRUE/FALSE, 336 cells.
    - 328 cells are identical.
    - The 8 that differ are exactly {IRD, threshold > 1, user-set}, on all four outcome types.
    - "non-C2 warnings identical in every cell: TRUE". The same check also compares the return values (`v_compare.R:46-50`), so they are identical too.
  - **Threshold-resolution probe.** This is the 175-cell extended matrix of `helper-threshold-sync.R`, built from statements lifted out of `body(forestsearch)`, so it includes the edited remap block. The result: "probe matrix identical start vs fixed: TRUE". 8 of its 50 RD/IRD cells carry the remap tag (`evidence/C2/verify/probe_compare.txt`).
  - **Walls:**
    - Fits: start 49.50 s, fixed 49.51 s.
    - Validator grid: 2.41 s / 2.61 s.
    - Compare: 0.22 s.
    - Probe: 2.86 s / 2.83 s. These two are stated only in VERIFY.md; no log in `evidence/C2/verify/` records them.

### Existing tests

All runs used the conventions command with `NOT_CRAN=true`. The results:

| run | files: pass / fail / skip | wall | evidence |
|---|---|---|---|
| implementer batch 1, fixed copy | binary-default-or-entry-points 25/0/0, binary-default-or 59/0/1, dina-gaussian-floor 55/0/0, threshold-pair-directive-a 85/0/1, threshold-sync 16/0/0 | 25.08 s | `evidence/C2/existing_batch1.txt` |
| the same batch, start-tree control | the same counts and skips | 24.74 s | `existing_batch1_start.txt` |
| implementer batch 2, fixed copy | declaration-c0 71/0/0, directive-c 154/0/1, estimability-boundary 60/0/0, fs_dgm_scale 45/0/0, fs-dgm-feasibility 40/0/0, lean-fits 362/0/0, rd-standard-error 14/0/0 | 48.29 s | `existing_batch2.txt` |
| implementer, with the two `dev/reports/` baselines placed, fixed | directive-c 182/0/0, threshold-pair-directive-a 86/0/0 | 12.26 s | `existing_baselines_fixed.txt` |
| the same, start control | directive-c 182/0/0, threshold-pair-directive-a 86/0/0 | 12.38 s | `existing_baselines_start.txt` |
| verifier batch 1, fixed (baselines and the campaign template placed) | binary-default-or-entry-points 25, binary-default-or 66, dina-gaussian-floor 55, sg-focus-transparency 64, threshold-pair-directive-a 86, threshold-remap-messages 27, threshold-sync 16; 0 fail, 0 skip | 32.20 s | `evidence/C2/verify/existing_b1.txt` |
| verifier batch 2, fixed | declaration-c0 71, directive-c 182, estimability-boundary 60, fs_dgm_scale 45, fs-dgm-feasibility 40, lean-fits 362, rd-standard-error 14; 0 fail, 0 skip | 45.97 s | `existing_b2.txt` |

- **Skips.** The implementer's three skips are files missing from the start export:
  - `test-binary-default-or.R:188:3`: the campaign template is not in the tree.
  - `test-threshold-pair-directive-a.R:148:3` and `test-directive-c.R:667:3`: the `dev/reports/` baseline CSVs.
  - The baselines came from `git show b06de4e8:dev/reports/...`. They were placed in the scratch copies only and removed afterwards, so they are not in the patch. The verifier also placed the campaign template. That is why it has binary-default-or at 66 with no skip, where the implementer has 59 with one skip.
- **Existing tests reach the changed lines.**
  - `helper-threshold-sync.R` lifts the RD/IRD branch out of `body(forestsearch)` (start `:63`). Its "default" cells supply c1 = 1.25 for RD and IRD, so they run the edited remap.
  - `probe_resolver_reach.txt` shows the raw messages on both trees. The tag stays `c1-ratio-scale-remap` (`.probe_warn_tag()`, start `:215-216`, keys on the opening word "effect.threshold"). Screening is 0.05 / 0.01 on both trees.
  - `test-dina-gaussian-floor.R` runs RD and IRD with `effect.threshold = 1.5` through `forestsearch()` (start `:142`, `:195`).
- **The deparsed-body tests pass on the fixed tree** (binary-default-or check 5, sg-focus-transparency, threshold-pair-directive-a). The new `switch(effect_measure, RD = "OR", ...)` does not match their patterns (VERIFY Check 5).
- **`test-mr-inference.R` was not run.** Its only RD reference is `.fs_mr_confirm_null("RD")` at start `:135`, and it supplies no RD or IRD threshold. I confirmed this with `git grep` at `b06de4e8`.
- **Expectations updated: none.**

### What a user sees

`forestsearch()` on synthetic count data with `effect_measure = "IRD"` and `effect.threshold = 1.5` (`evidence/C2/unaffected/compare.txt`; the two warnings are separated by `--`).

Start tree:
```
effect.threshold = 1.50 but effect_measure = 'IRD' (identity scale, range [-1, 1]).
  A risk difference cannot exceed 1.0.
  Did you mean effect_measure = 'OR' (ratio scale)?
--
effect.threshold = 1.50 appears to be on a ratio scale, but effect_measure = 'IRD' uses the identity scale (range [-1, 1]).
  Did you mean effect.threshold = 1.50 on the OR scale?
  Remapping to default: 0.05.
  To use a custom IRD threshold, pass a value in [-1, 1] (e.g., effect.threshold = 0.07 for a 7 pct-point threshold).
```

Fixed tree:
```
effect.threshold = 1.50 but effect_measure = 'IRD' (identity scale, range [-1, 1]).
  A rate-difference threshold above 1.0 is outside that range.
  Did you mean effect_measure = 'IRR' (ratio scale)?
--
effect.threshold = 1.50 appears to be on a ratio scale, but effect_measure = 'IRD' uses the identity scale (range [-1, 1]).
  Did you mean effect.threshold = 1.50 on the IRR scale?
  Remapping to default: 0.01.
  To use a custom IRD threshold, pass a value in [-1, 1] (e.g., effect.threshold = 0.02 for a rate difference of 0.02).
```

- Both trees record `threshold_config$screening = 0.01` and `args_call_all$effect.threshold = 0.01`.
- For RD 1.5, the two messages are byte-identical on both trees (`compare.txt`, the "RD_out warnings" blocks).

### Verification

- **Verdict:** CONFIRMED, round 1 (`evidence/C2/verify/VERIFY.md`).
- `workflow1_result.json` (entry C2): `status` READY, `final_verdict` CONFIRMED, `repaired` false, `blocking_r1` [], `blocking_final` [], 7 entries in `minor_final`.
- **No repair round.**
- **No polish.** C2 is not among the 11 entries of `polish_result.json` (A1, A2, A3, B1, C1, C3, C4a, C4b, C4c, C6, C7).
- **The commit equals the gated patch.**
- **The seven minor problems, and what happened to each:**
  1. **The validator names no default** (`R/forestsearch_main.R:301-317`). The verifier rated it a justified departure: it ran the fits SRV_RD_1.5 and SRV_IRD_1.5, where the validator fires on survival, no remap follows, and the fit screens at log(1.5) = 0.4054651. **Outcome:** unchanged at head; reported as a departure (below).
  2. **Contradictory adjacent comments.** The start comment "RD/IRD are probability/rate differences bounded in [-1, 1]" sits directly above the new comment saying a rate difference is not bounded by 1. The same premise appears at the remap NOTE. **Outcome:** unchanged at head (`fb80dd8f:296` against `:302-304`, and `:2369`). C4 did not touch these lines, and they are not in C4's table. Open item.
  3. **NEWS wording.** The proposed bullet said "the two warnings" print 0.01, but only the remap warning does. **Outcome:** resolved. The bullet at `fb80dd8f`, `NEWS.md:453-458`, gives the 0.01 to the remap warning only.
  4. **Hunk count.** NOTES §2 says two hunks; the patch has three hunk headers in two regions. **Outcome:** corrected here (see Site).
  5. **Off-by-one pointers in NOTES.** **Outcome:** corrected here (see Corrections to the notes).
  6. **The IRD block asserts fragments, not the full strings** (test `:74-100`), so a malformed IRD example clause would pass. The RD block asserts full text. **Outcome:** the test is unchanged at head. Open item.
  7. **The survival suggestion.** With `outcome_type = "survival"` and a supplied `effect_measure = "IRD"`, the validator now suggests `'IRR'`, which survival also ignores. The start suggested `'OR'`, which had the same pre-existing oddness, so this is not a regression. **Outcome:** unchanged. Open item.
- **Other verifier checks:**
  - The start export matches `b06de4e8` byte for byte (`cmp`).
  - The names `.fs_remap_*`, `bound_line` and `ratio_measure` occur nowhere else under `R/` or `tests/` at the start SHA.
  - C2 and C4b apply cleanly in either order (`verify/C2/combo_c4b`, `combo_c4b_rev`).
  - The new lines are ASCII only.

### Departures and interpretations

- **The validator names no default.** The task asks that "each message name the measure and the default that apply". The validator message never carried a default, and none was added.
  - The reason: the validator also fires where no remap follows, with `outcome_type = "survival"` and a supplied `effect_measure` of RD or IRD. The remap sits inside `if (outcome_type != "survival")` (start `:2251`). A default printed there would be false.
  - The verifier confirmed this by running the fits (above).
  - On every GLM call, the remap message follows immediately and names the default.
  - The task's test sentence ("the messages for a risk difference and for a rate difference carry the defaults") is met by the remap message for each measure.
- **The IRD validator line.** The RD line "A risk difference cannot exceed 1.0." is kept. For IRD, "A rate difference cannot exceed 1.0" would be false, because a rate difference is not bounded by 1: `poisson_offset_delta` computes IRD as `lam1 - lam0`, with `lam = exp(b0 [+ b1] + mean(log t))` (start `R/glm_effect_estimators.R:896-905`). The IRD line says only what the code does: "A rate-difference threshold above 1.0 is outside that range."
- **The IRD example in the remap's last line.** The task does not mention it. The RD example ("0.07 for a 7 pct-point threshold") is in risk-difference units. For IRD it is now "effect.threshold = 0.02 for a rate difference of 0.02", which makes no unit claim.
- **The ratio counterpart for IRD is IRR.** The task says only "name the measure". IRR comes from `@param effect_measure`, where IRR and IRD are listed together as the ratio and difference measures of the same rate model.
- **The assignment reads the local** (`effect_threshold <- .fs_remap_default`) instead of repeating the `switch()`. The value is unchanged: the fits, `threshold_config`, `args_call_all` and the 175-cell probe are identical.
- **The consistency remap message is unchanged.** "Remapping to default: 0.0." (start `:2361`) is right for both RD and IRD: c2 = 0.0 for both in the help table.
- **The test calls `forestsearch()` itself** (`maxk = 1`, `fs.splits = 20`, N = 200). No smaller entry point was needed: the 27 expectations ran in 7.7-8.1 s in the green runs.

### Articles

- The C2 evidence contains no article statement.
- My own check, by `git grep` at `fb80dd8f`: no file under `vignettes/` passes `effect_measure = "RD"` or `"IRD"`, and none contains the string `IRD`. The `effect_measure` values given there are `"OR"`, `"MD"`, `"HR"` and `NULL`.
- `NULL` resolves to OR, MD or IRR, never to RD or IRD.
- Both changed sites need RD or IRD: `is_bounded_identity` for the validator, `effect_measure %in% c("RD", "IRD")` for the remap.
- So no evaluated article call reaches the changed lines, and no article output can change. This comes from reading the code; nothing was rendered.

### Open items

1. **Check 3 of the same validator** (`R/forestsearch_main.R:257-293`, the same at `b06de4e8` and `fb80dd8f`) has the same kind of wrong-measure text for count outcomes. For a ratio threshold below 0.5 with `effect_measure = "IRR"`, it suggests `outcome_type = 'binary', effect_measure = 'RD'` (c1, `:272`) and `effect_measure = 'RD'` (c2, `:288`), where IRD would apply. This is outside C2, which names two messages.
2. **The `[-1, 1]` premise for IRD.**
   - The comments at `fb80dd8f:296` and `:2369` ("RD/IRD are probability/rate differences bounded in [-1, 1]") are not true for a rate difference. They contradict the new comment at `:302-304`.
   - The remap treats any IRD threshold above 1 as a ratio-scale value. A legitimate IRD threshold above 1 (possible for high rates or long follow-up) is replaced by 0.01, with a warning.
   - Behaviour is unchanged by design, because the remapped values must not change.
3. **The IRD block asserts fragments** (`tests/testthat/test-threshold-remap-messages.R:83-97`), not the full strings. A malformed IRD example clause would still pass.
4. **Survival with a supplied IRD** gets the suggestion `effect_measure = 'IRR'`, which survival ignores (it was `'OR'` at the start). This is pre-existing.

### NEWS entry

The bullet is at `fb80dd8f`, `NEWS.md:453-458`, under "Behaviour changes" (`:154`). Each claim, checked:

| claim | check | result |
|---|---|---|
| "`forestsearch()` with `effect_measure = "IRD"` and an `effect.threshold` above 1" | Both conditions are `> 1.0` (head `:301`, `:2371`). 1.0 itself draws no message (verifier IRD_1.0: 0/0). | correct; see note 2 |
| "the validator warning now names the rate difference" | head `:307`: "A rate-difference threshold above 1.0 is outside that range." | correct |
| "... and the IRR scale" (validator) | head `:308`, `:314`: "Did you mean effect_measure = 'IRR' (ratio scale)?". The words "IRR scale" are only in the remap warning. | loose; see note 1 |
| "(it said "A risk difference cannot exceed 1.0" and suggested OR)" | start `:306-307` | correct |
| "the remap warning names the IRR scale and prints the 0.01 it assigns" | head `:2375-2376`, `:2386-2387`, `:2391`, `:2394`; `compare.txt` | correct |
| "(it suggested the OR scale and printed "Remapping to default: 0.05")" | start `:2341-2342`, `:2346`; `compare.txt` | correct |
| "The risk-difference messages and every remapped value are unchanged." | RD texts byte-identical (the test's RD block passes on red; `compare.txt`; `v_compare.txt`). Fits, `threshold_config`, `args_call_all` and the 175-cell probe are identical. | correct |

- **Note 1.** The validator suggests the measure IRR ("effect_measure = 'IRR' (ratio scale)"); it does not say "IRR scale". A closer wording: "the validator warning now speaks of a rate-difference threshold and suggests `effect_measure = 'IRR'`".
- **Note 2.** The bullet names only `effect.threshold`. The legacy `hr.threshold` above 1 draws the same two changed messages (verifier IRD_hr1.5). With `outcome_type = "survival"` only the validator message fires (verifier SRV_IRD_1.5). The bullet does not claim otherwise, so this is not an error.
- **Omission.** The bullet does not mention that the remap's last line changes for IRD: "(e.g., effect.threshold = 0.07 for a 7 pct-point threshold)" becomes "(e.g., effect.threshold = 0.02 for a rate difference of 0.02)". A user sees this change.

### Corrections to the notes

- **Hunk count.** NOTES §2 says "Hunks: two, at start 299–311 and 2332–2352". `git diff b06de4e8 2984dbee` has three hunk headers (`@@ -299`, `@@ -2332`, `@@ -2341`), in two regions.
- **`args_call_all`.** NOTES §2 gives `mget(names(formals()))` at start `:1923`. That expression is two lines: `args_names <- names(formals())` at `:1923` and `args_call_all <- mget(args_names, envir = environment())` at `:1924`. VERIFY cites `:1924`.
- **Survival guard.** NOTES §8 puts `if (outcome_type != "survival")` at start `:2250`. It is at `:2251`.
- **The red passes.** NOTES §3 lists 14 expectations as "passing on red as well as green" (RD 6, in-range 6, IRD recorded values 2). There are 18. The other four are `:78`, `:79` (`expect_length`) and `:83`, `:91` ("effect_measure = 'IRD'").
- **The NEWS bullet.** NOTES §10's bullet was reworded in `fb80dd8f`, as the verifier asked (problem 3).

## C3 — `sg_tables()` GRF CATE note names the count scale given to the forest

`S` = `/tmp/cran_prep_displays_se_docs_2026-11-12`. Line numbers are at the start SHA `b06de4e8` and at head `fb80dd8f`. `R/` is identical at `835fac1b` and `fb80dd8f` (both trees `22509e66`). The test file is the same blob (`621e495e`) at `c3323a59` and `fb80dd8f`, so head test line numbers equal those at the commit. The R block moved from `:924-956` at `c3323a59` to `:946-978` at head, because C4 (`ba7c4ea2`) and C5 (`5073fd67`) added 22 roxygen lines above it (`git diff c3323a59 fb80dd8f -- R/summary_utility_functions.R`: one hunk at `:733`, +23/-1). The block itself is byte-identical at the commit and at head.

### C3.1 Status, commit, provenance

| | |
|---|---|
| Status | **PASS** |
| Commit | `c3323a59` "fix: sg_tables() GRF CATE note names the count scale given to the forest" |
| Parent | `9ec67b7e` (C1) |
| Files (`git diff --numstat 9ec67b7e c3323a59`) | `R/summary_utility_functions.R` +25/-2; `tests/testthat/test-grf-footnote-scale.R` +201 (new). 2 files, +226/-2. (`--stat` shows 27 changed lines for the R file: 25 + 2.) No `man/` change (no roxygen edited). No non-ASCII character in the added lines (checked for this section). |
| Commit vs gated patch | The changed lines of `git diff 9ec67b7e c3323a59` equal those of `S/patches/C3.patch` (checked for this section; `S/main/C3/patch_equal.txt` is empty). The hunk header differs only by offset: `-923` in the commit, `-911` in the patch, because the committed A1 (`f4c7dd2d`) adds 12 lines above the site. |
| Final patch | `S/patches/C3.patch`, sha256 `e6e0df96...18b517f3` (after polish). |
| Earlier patch | `S/patches/C3.verified.patch`, sha256 `f214eaf8...d78df84453` (verified; test file 197 lines). |
| Pre-set condition | "If the object does not record enough to say, report C3 BLOCKED." It holds: every return path of `grf.subg.harm.glm()` records `outcome_type`, `offset.name`, `grf_count_transform`, `grf_y_transform` and `adverse_outcome` (start `R/grf_subg_harm_glm.R:546-551`, `:569-574`, `:640-647`, `:704-710`). `grf_y_transform` is documented as the "transformation actually applied to Y for GRF screening" (`@return`, `:157-159`). `9409d2cc`, the commit that created the file, already sets it. |

### C3.2 What GRF is given for a count (start `R/grf_subg_harm_glm.R`)

Traced by the implementer and independently by the verifier. Re-read for this section.

| Question | Code |
|---|---|
| Raw count, rate, or offset-adjusted? | The per-subject count: `Y_raw <- data[[outcome.name]]` (`:401`). It is not divided by exposure, and no log-offset enters it. |
| Transform | `grf_count_transform = "log"` (default): `Y_grf <- log(Y_raw + 0.5)`, `grf_y_transform <- "log(Y + 0.5)"` (`:410-412`). `"identity"`: `Y_grf <- Y_raw`, `grf_y_transform <- "identity"` (`:421-423`). Other outcome types: `Y_grf <- Y_raw`; `grf_y_transform` stays `"none"` (`:407`, `:433`). |
| What `offset.name` does | Argument checks (`:338-345` warning, `:354-356` stop, `:362` required columns); the complete-case filter (`:372-373`, applied at `:380`), so rows with a missing exposure are dropped before the forest; a verbose message (`:390-391`); the four return lists (`:547`, `:570`, `:643`, `:706`). It never enters `Y_grf`. |
| Negated? | With `adverse_outcome = TRUE` a count gets `Y_grf <- -Y_grf` (`:452-457`), part 4's F10 (`f2e8de47`). This runs after `grf_y_transform` is set, so that field does not record the flip. Before `f2e8de47` no count was negated. |
| What is printed | Tree path: `cate_sg`, `cate_sgc` are means of DR-score contrasts of the forest fitted to `Y_grf` (`:605-606`). Frontier path: `cate_sg = -sel$effect`, `cate_sgc = NA_real_` (`:576`). |
| `effect_measure` | Selects the GLM only. The verifier's run with `log_RR` gave the same `Y` and the same CATE as `log_IRR` (`S/evidence/C3/verify/capture_y_fixed.txt`). |

Runtime confirmation (verifier, `verify/scripts/capture_y.R`; `grf::causal_forest` traced, exposure `pt ~ U(0.5, 2)`): the forest receives `log(events + 0.5)` for the log fits (not `log(events/pt + 0.5)`, not `events/pt`), `events` for identity, `-log(events + 0.5)` and `-events` under `adverse_outcome = TRUE`. The CATE with and without the offset is the same (-1.0640 / 1.1389). Walls (`proc.time`): fixed 4.93 s, start 5.20 s.

### C3.3 Site, before → after

Function `sg_tables()`, the `if (inherits(fs, "grf_glm_result"))` branch (start `:811`, head `:833`). The locals it reads are set at start `:817-818` / head `:839-840`: `outcome_type <- fs$outcome_type %||% "binary"`, `offset.name <- fs$offset.name`. Nothing reassigns them before the note.

Before, `b06de4e8:R/summary_utility_functions.R:912-921`:

```r
    # CATE footnote (GRF-specific context)
    if (!is.null(fs$cate_sg)) {
      tab_estimates_grf <- tab_estimates_grf |>
        gt::tab_source_note(
          source_note = gt::md(sprintf(
            "**GRF CATE:** H&#x302; = %.4f, H&#x302;<sup>c</sup> = %.4f (risk-difference scale)",
            fs$cate_sg, fs$cate_sgc
          ))
        )
    }
```

After, `fb80dd8f:R/summary_utility_functions.R:946-978` (comment lines `:948-955` shortened here to the first and last):

```r
    # CATE footnote (GRF-specific context)
    if (!is.null(fs$cate_sg)) {
      # The CATE is on the scale of the outcome given to causal_forest().
      ...
      # (counts were not negated before part 4's fix F10, commit f2e8de47).
      # Other outcome types keep the original label.
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
        }
      }
      tab_estimates_grf <- tab_estimates_grf |>
        gt::tab_source_note(
          source_note = gt::md(sprintf(
            "**GRF CATE:** H&#x302; = %.4f, H&#x302;<sup>c</sup> = %.4f (%s)",
            fs$cate_sg, fs$cate_sgc, cate_scale
          ))
        )
    }
```

- For every non-count type `cate_scale` is `"risk-difference scale"`, so `sprintf()` returns the start string.
- For a count, the parenthesis is one of: `log(Y + 0.5) scale, Y = count`; `count scale`; either with `, not offset-adjusted` when `offset.name` is recorded; `scale not recorded` when `grf_y_transform` is missing or not one of the two values.
- Only `grf.subg.harm.glm()` builds class `grf_glm_result`. The only `GRF CATE` string under `R/` is this one (start `:917`, head `:974`). The survival GRF does not return that class, so survival never prints this note.

### C3.4 Red and green runs

Command for every run: `NOT_CRAN=true ... devtools::test('<tree>', filter = 'grf-footnote-scale', reporter = 'summary', stop_on_failure = FALSE)`. The implementer and verifier ran the verified test file (197 lines); the polish runs and the main session ran the committed one (201 lines, 4 more header comment lines, so failure lines are +4).

| Run | Tree | Result | Wall |
|---|---|---|---|
| Main, red | repository at `9ec67b7e` + the committed test file only | `....12...3...4....56...`: 23 expectations, 17 pass, 6 failures, 0 errors. Failures at `:127`, `:131`, `:142`, `:155`, `:171`, `:182` (`S/main/C3/red.txt`) | 5.3 s |
| Main, green | the same + the code (= `c3323a59`) | `.......................`: 23 pass (`S/main/C3/green.txt`) | 5.0 s |
| Implementer, red | `S/wt/C3` = start + test file only | 6 tests, 23 expectations, 17 pass, 6 failures, 0 errors; failures at `:123`, `:127`, `:138`, `:151`, `:167`, `:178` (`evidence/C3/red.txt`) | 5.47 s |
| Implementer, green | `S/wt/C3` with the fix | 6 tests, 23 pass, 0 failures, 0 errors, 0 skips (`green.txt`) | 5.62 s |
| Implementer, green again | same | 6 tests, 23 pass (`existing_counts.txt`) | 5.37 s |
| Verifier, red | `S/verify/C3/red` = start + test file only | 6 tests, 23 expectations, 6 failures, 0 errors, 0 skips; same six lines (`verify/red.txt`, `redgreen_counts.txt`) | 5.63 s |
| Verifier, green | `S/verify/C3/fixed` = start + whole patch | 6 tests, 23 expectations, 0 failures (`verify/fixed.txt`) | 5.73 s |
| Verifier, A1 + C3 | `S/verify/C3/combo_a1c3` | `test-grf-footnote-scale.R` 23/23, `test-sg-tables-orientation.R` 155/155 (`verify/combo_apply.txt`) | 7.79 s |
| Polish, red | `S/wt/C3_polish_red` = start + polished test file only | `....12...3...4....56...`: 6 tests, 23 expectations, 17 pass, 6 failures, 0 errors (`polish/red.txt`) | 5.30 s |
| Polish, green | `S/wt/C3_polish` | 6 tests, 23 expectations, 23 pass (`polish/polish.txt`) | 5.18 s |
| Polish check, green | `S/verify/C3_polish` = fresh start + final patch | 6 tests, 23 expectations, 0 failed, 0 skipped, 0 errors (`verify_polish/green.txt`) | 5.23 s |

First failure, main session (`S/main/C3/red.txt`):

```
── 1. Failure ('test-grf-footnote-scale.R:127:3'): count fit, log transform, off
Expected `notes` to be identical to `.c3_expected(fit, "log(Y + 0.5) scale, Y = count, not offset-adjusted")`.
- "**GRF CATE:** H&#x302; = -1.0322, H&#x302;<sup>c</sup> = 1.0631 (risk-difference scale)"
+ "**GRF CATE:** H&#x302; = -1.0322, H&#x302;<sup>c</sup> = 1.0631 (log(Y + 0.5) scale, Y = count, not offset-adjusted)"
```

- The other five red failures: `:131` (`"risk-difference"` still present), `:142` (count, `adverse_outcome = TRUE`), `:155` (identity), `:171` (no offset), `:182` (unrecorded transform). Each count note still reads `(risk-difference scale)`.
- The binary test (`:185-201`) passes on both trees, as an unchanged case should.
- The red numbers equal the start tree's: the parent `9ec67b7e` includes A1, which does not touch `cate_sg` or `cate_sgc`.

### C3.5 Expected values

All expected strings are exact; no tolerance is used (`expect_identical()`). The working is in the test file header, `fb80dd8f:tests/testthat/test-grf-footnote-scale.R:1-52`: code trace `:9-28`, expected text `:30-40`, why the sign is not named `:41-44`, data `:46-52`. Comments also sit beside the assertions (`:120`, `:138`, `:162`, `:179-180`, `:198-199`) and the helper (`:109-110`).

| Test (head lines) | Case | Premises asserted | Expected text in the parenthesis | Derivation (by hand, from the code at the start SHA) |
|---|---|---|---|---|
| `:118-132` | count, `grf_count_transform = "log"`, offset `person_time` | `:121-124`: `outcome_type`, `grf_y_transform`, `offset.name`, `cate_sg` not NULL | `log(Y + 0.5) scale, Y = count, not offset-adjusted` (`:129`); no `"risk-difference"` (`:131`) | GRF gets `log(Y_raw + 0.5)` of the per-subject count (`:401`, `:410-412`). The offset never enters `Y_grf`, so the clause is added. |
| `:134-146` | count, `adverse_outcome = TRUE` | `:137`, `:139`, `:140` | the same (`:144`) | Negating `Y_grf` (`:452-457`) does not change the scale. It runs after `grf_y_transform` is set, so that field is unchanged (`:138-139`). |
| `:148-159` | count, `"identity"`, offset | `:151-153` | `count scale, not offset-adjusted` (`:157`) | `Y_grf <- Y_raw` (`:421-423`). |
| `:161-175` | count, no offset | `:163-166` (`expect_warning(..., "offset.name = NULL")`, the warning at `:338-345`), `:167-169` | `log(Y + 0.5) scale, Y = count` (`:173`) | No offset recorded, so no clause. |
| `:177-183` | count, `grf_y_transform` set to NULL by hand | none | `scale not recorded` (`:182`) | The note must not name a scale the object does not record (`:179-180`). |
| `:185-201` | binary, default | `:195-196` | `risk-difference scale` (`:200`) | The start-SHA string (`:917`). GRF on a 0/1 outcome estimates a difference in risks. Task: other outcome types unchanged. |

- The leading part and the `%.4f` format are copied from the start SHA's format string (`:917`) into the helper `.c3_expected()` (`:111-116`).
- The two numbers are `fit$cate_sg` and `fit$cate_sgc`, computed by `grf.subg.harm.glm()`, which this fix does not touch. They are inputs of the note, not values of the fix. The verifier accepted this ("None is copied from the fixed display").
- Data: the count design of part 4's `test-grf-count-orientation.R` (n = 1000; region `z`; Poisson means 3 control, 9 treated in `z == 1`, 1 treated in `z == 0`; seed `20261112L`, `:56-71`) and a binary analogue (risks 0.4 control, 0.8 / 0.1 treated; seed `20261113L`, `:73-87`). Seeds were fixed once. No assertion depends on a CATE value.
- Six tests, 23 expectations: 17 premise and unchanged-case checks, 6 that need the fix.

**The verifier's check that the printed numbers are on the named scale** (`verify/exact_scale.txt`, 0.15 s; `verify/sample_scale.txt`, 3.58 s; same generator and seed as the test):

| | Difference on `log(Y + 0.5)` | Difference on the count | Log ratio of means |
|---|---|---|---|
| Population, by hand (`E[log(Y + 0.5)]` = 0.1689, 1.1005, 2.1965 for means 1, 3, 9) | z = 0: -0.9316; z = 1: 1.0960 | -2 and 6 | -1.0986 and 1.0986 |
| Sample, z = 0 | -1.0108 | -2.1160 | -1.2283 |
| Sample, z = 1 | 1.0502 | 5.7920 | 1.0348 |
| Printed by the fits | log fit -1.0322 / 1.0631 | identity fit -2.1554 / 5.8495 | — |

- The log fit's numbers match the `log(Y + 0.5)` column, and the identity fit's match the count column. The policy-tree `H` is close to, not equal to, `z == 0`.
- A "risk difference" of -2.1554 is impossible (it must lie in [-1, 1]), so the old label was false.

### C3.6 Unaffected case

The task names it: "Other outcome types are unchanged." Each tree ran in its own `Rscript --vanilla` process and saved RDS; `identical()` ran in a separate process. A second start run is the run-to-run control.

**Implementer** (`evidence/C3/unaffected/run.R`, `run_surv.R`, `compare.R`; output `compare.txt`, `run_*.txt`):

| Cases | Compared | start vs fixed |
|---|---|---|
| Binary: tree default, `adverse_outcome = TRUE`, frontier, `subgroup_notation = "harm"`, no subgroup (`dmin.grf = 10`). Continuous: tree default, `adverse_outcome = TRUE`, frontier. (8) | whole `sg_tables()` list, and raw HTML of Table 1 | all `identical()` TRUE; HTML identical |
| Count, no subgroup (no note printed) (1) | whole list | TRUE |
| Count: log + offset, log + adverse, identity, frontier (4) | every component of Table 1 except `_source_notes`; `sg10_out` | identical; `_source_notes` is the only difference |
| Survival `forestsearch()` fit, FS path (helper DGM, N = 250, HR_harm = 2.5, seed 42; subgroup `{age <= 55}`) (1) | Table 1 whole; `sg10_out` without `_formats`; raw HTML of both tables | identical |

- Control (start vs start2): all 13 GRF cases identical; survival identical except `sg10_out[["_formats"]]`.
- Dropped leaf: `sg10_out[["_formats"]]` (gt `fmt_number()` closures), survival case only. It differs between two runs of the same start tree.
- Walls: `run.R` start 8.01 s, start2 8.07 s, fixed 7.95 s; `run_surv.R` start 5.41 s, start2 5.29 s, fixed 5.49 s; `compare.R` 0.17 s (stated in NOTES §5; no file in `unaffected/` records it).

**Verifier** (own generators, n = 800, harm region `x1 > 0.3 & x2 == 1`; `verify/scripts/run_unaff.R`, `compare_unaff.R`; output `verify/compare_unaff.txt`):

| Cases | start vs start2 | start vs fixed |
|---|---|---|
| Binary: tree (notation none / harm / benefit), `adverse_outcome = TRUE`, frontier, `RD`, `log_RR` + adverse, no subgroup (8) | all TRUE | whole list `identical()` TRUE and raw HTML identical |
| Continuous: tree (none / harm / benefit), adverse, frontier (5) | all TRUE | whole list and HTML identical |
| Count, no subgroup (1) | TRUE | whole list identical |
| Count with a note: log + offset (none / harm / benefit), identity, log + adverse, log without offset, frontier, `log_RR` (8) | TRUE | everything except `_source_notes` identical; `sg10_out` identical |
| Survival `forestsearch()`, FS path (N = 300, HR_harm = 3, seed 7; subgroup `!{age <= 63} & !{stage}`) (1) | `tab_estimates` TRUE; `sg10_out` FALSE through `_formats` | `tab_estimates` identical; `sg10_out` identical without `_formats`; HTML of both tables identical; no `GRF CATE` text |

- Walls (`/usr/bin/time`, `verify/run_unaff_*.txt`): start 9.46 s, start2 9.94 s, fixed 9.79 s; compare 1.54 s. (VERIFY.md gives start as 9.63 s; see C3.13.)

**Cumulative re-run at head** (`S/cumulative/C3/RESULT.md`, written after the version commit; the same scripts, byte-identical except two path lines in the verifier's compare script):

- 37 cases (implementer 13 + 1 survival, verifier 22 + 1 survival), head = `git archive 835fac1b` (`R/` equal to `fb80dd8f`).
- 9 identical to start (all `adverse_outcome = TRUE` fits, `count_nosg`, `cnt_nosg`, survival under the `_formats` normalisation).
- 16 differ through A1 only: binary and continuous fits with `adverse_outcome = FALSE`; only `_data` differs, in the effect column. The GRF note text is unchanged.
- 12 differ through C3 only: the count cases with a note; only `_source_notes` differs. The head prints exactly C3's fixed-tree note.
- Reverse-applying A1 and C3 from head (`trees/minus_A1_C3`) makes all 37 cases identical to start. Reverting only A1 leaves the 12 C3 cases; reverting only C3 leaves the 16 A1 cases.
- Walls (`/usr/bin/time`): `run.R` start 8.32 s, start2 8.31 s, head 8.27 s; `run_surv.R` 5.72 s, 5.64 s, 7.69 s; `run_unaff.R` 10.92 s, 9.91 s, 9.89 s; compares 0.19 s, 1.61 s, 1.51 s.

### C3.7 Existing tests

No existing test calls `sg_tables()` on a `grf_glm_result`. The only test file that calls `sg_tables()` at the start SHA is `test-no-subgroup-dina-grf-summary.R` (FS path). No test references `_source_notes`, `GRF CATE` or `risk-difference scale` (`git grep` at `b06de4e8`). All runs with `NOT_CRAN=true`.

| File | Why | Implementer (`existing_counts.txt`) | Verifier (`verify/existing_*.txt`) |
|---|---|---|---|
| `test-no-subgroup-dina-grf-summary.R` | calls `sg_tables()` | 5 tests, 21 pass, 0 fail; 5.44 s | 5 / 21, 0 / 0 / 0; 5.36 s |
| `test-grf-count-orientation.R` | count `grf.subg.harm.glm()` fits | 2 tests, 14 pass; 4.40 s | 2 / 14, 0 / 0 / 0; 4.51 s |
| `test-glm-pipeline-integration.R` | GRF GLM pipeline | 4 tests, 10 pass; 5.88 s | 4 / 10, 0 / 0 / 0; 5.81 s |
| `test-grf-inclusion-band.R` | `grf.subg.harm.glm()` selection | 6 tests, 2875 pass; 9.33 s | 6 / 2875, 0 / 0 / 0; 9.48 s |

- The implementer's first pass over the same four files (`existing.txt`, rc 0 each) has walls 5.68 s, 4.54 s, 5.79 s, 8.93 s.
- Expectations updated: none.
- Full suite at head (`S/suite/suite.log`, still running when this section was drafted): `✔ | 23 | grf-footnote-scale [1.8s]`, `✔ | 10 | glm-pipeline-integration`, `✔ | 14 | grf-count-orientation`, `✔ | 2875 | grf-inclusion-band`. The final suite result is in the report's suite section.

### C3.8 What a user sees

Count fit, `grf.subg.harm.glm(outcome_type = "count", effect_measure = "log_IRR", offset.name = "person_time")`, default log transform; the source note of `sg_tables(fit)$tab_estimates` (`evidence/C3/unaffected/run_start.txt`, `run_fixed.txt`):

```
start: **GRF CATE:** H&#x302; = -1.0322, H&#x302;<sup>c</sup> = 1.0631 (risk-difference scale)
fixed: **GRF CATE:** H&#x302; = -1.0322, H&#x302;<sup>c</sup> = 1.0631 (log(Y + 0.5) scale, Y = count, not offset-adjusted)
```

- Rendered in the table footer (verifier, `verify/html_notes.txt`, another data set): `GRF CATE: Ĥ = -0.4194, Ĥc = 0.1766 (log(Y + 0.5) scale, Y = count, not offset-adjusted)`.
- With `grf_count_transform = "identity"`: `... = -2.1554, ... = 5.8495 (count scale, not offset-adjusted)`; before, `(risk-difference scale)`.
- Binary and continuous fits print the same note as before, for example `H = -0.3335, Hc = 0.3437 (risk-difference scale)` on both trees.

### C3.9 Reach and articles

- Only `sg_tables()` on a `grf.subg.harm.glm()` result with a count outcome and a found subgroup reaches the changed text.
- Articles: the articles sweep on the head tree (`S/articles/C3_work/sweep_head_fixes.txt`, 13 sources under `vignettes/`) finds one that reaches `sg_tables()`: `vignettes/articles/forestsearch-2024.Rmd`, chunk `fs-results` (`:412-433`), on a `forestsearch()` fit. That takes the FS path, which prints no GRF note. No article output can change.
- `quarto/` (read-only, not a package article): `quarto/applications/count_data_demo.qmd:747` calls `sg_tables(fs)` only when `forestsearch()` found a subgroup; its count GRF results `grf_pois` and `grf_id` (`:349`, `:380`) are not tabulated. The ACTG175 files that call `sg_tables(grf_result)` are binary or continuous.

### C3.10 Verification, repair, polish

- **Verdict:** CONFIRMED (`verify/VERIFY.md`; `S/workflow1_result.json` entry C3: `final_verdict` CONFIRMED, `repaired` false, no blocking problem in either list).
- **Repair rounds:** none.
- **The verifier's own decision on the sign:** the note need not state it. Negation changes orientation, not scale. The object does not record whether a count was negated: `adverse_outcome` is the request, and objects from before `f2e8de47` carry `adverse_outcome = TRUE` with an un-negated count (checked at `f2e8de47^`). The binary (`1 - Y`) and continuous (`-Y`) notes, kept by the task, are in the same position.

Four minor problems:

| # | Problem (verifier) | What happened |
|---|---|---|
| 1 | NOTES §2 gives the R file as "+27/-2"; `git diff --numstat` gives +25/-2. | Not edited (evidence prose; outside the polish list). This section uses +25/-2, checked against the commit. |
| 2 | The test file (verified copy `:17`) and NOTES §1 say `offset.name` "is used only by the GLM effect estimates"; at the start SHA it also selects the complete-case rows given to the forest (`:372-373`). | Fixed in the polish: test header rewritten, now `:16-21` (+6/-2 comment lines). |
| 3 | The note does not state orientation: with `adverse_outcome = TRUE` the printed count CATEs are differences in `-log(Y + 0.5)` or `-Y`, the opposite sign to a difference on the named scale (`H = -1.1480` printed where treatment raises the count). Exists for binary and continuous at the start SHA. | Left open (NOTES §9 item 2; C3.12 item 2). |
| 4 | The R comment "counts were not negated before 0.4.1.9010" and NOTES' "0.4.1.9009 or earlier" are loose: `f2e8de47` carries version 0.4.1.9009 and already negates (`88302212` bumps to 0.4.1.9010). | Fixed in the polish: the comment now reads `(counts were not negated before part 4's fix F10, commit f2e8de47).` (head `:954`; 2 lines → 2 lines). |

**Polish** (main session, comment-only; `evidence/C3/polish.diff`, 2 files, +8/-4):

- Test file `:16-17` → `:16-21`. Old: `(offset.name is used only by the GLM effect estimates).` New: apart from argument checks (`:338-366`) and a verbose message (`:390-391`), `offset.name` adds the exposure column to the complete-case filter (`:372-373`), so rows with a missing exposure are dropped before the forest, and is recorded for the GLM effect estimates.
- R comment `:920-921` (fixed copy; head `:954-955`): the version boundary replaced by the commit, as in minor 4.
- Gate: the comment-only check over the polish delta prints nothing; added lines are ASCII. No code line, expected value, seed or tolerance changed.
- Tests after polish: green 23/23 (5.18 s); red with the polished test only 17 pass / 6 fail, same messages at lines +4 (5.30 s).

**Polish check** (`S/polish_result.json`, C3 `check`): CONFIRMED, no problems.

- The fresh copy `S/verify/C3_polish` (start + final patch) equals `S/wt/C3_polish` (`diff -r -x .git` empty).
- Both new comments are true at the start SHA. Probe on the polished copy (`verify_polish/probe_offset.R`, `.txt`; 4.64 s): with no missing exposure, `offset.name = "person_time"` and `NULL` give results that differ only in `offset.name`; varying the exposure leaves the forest's `Y.orig` (= `log(events + 0.5)`) and `cate_sg` identical; with 5 exposures set to NA, 995 rows reach the forest (ids 3, 50, 77, 400, 999 dropped) against 1000 with `offset.name = NULL`.
- Green 23/23, 5.23 s.
- Two non-blocking notes, neither acted on (the commit equals the polished patch):
  - The recorded `offset.name` is also read by C3's own clause (fixed copy `:930`, head `:964`), so the test header's "recorded in the result for the GLM effect estimates" (head `:20-21`) is incomplete, not false. The header's expected-text section (`:36-37`) already says it.
  - "part 4's fix F10" is an internal task label in shipped R source (head `:954`); the commit hash beside it makes it checkable.

**Main session:** committed the final patch unchanged (`S/main/C3/patch_equal.txt` empty; changed lines re-checked here). Red and green on the repository as in C3.4.

### C3.11 Departures

From the task's wording ("For a count outcome, name the scale of the outcome given to GRF, read from what the object records. Other outcome types are unchanged."):

- **`, not offset-adjusted`** is added when `offset.name` is recorded. GRF is given the count, not a rate. With an offset recorded a reader could take the CATE as a rate difference. Read from `offset.name`.
- **`scale not recorded`** for a count object without a recognised `grf_y_transform`. Every object `grf.subg.harm.glm()` builds records it, so only a hand-built or modified object reaches this. The alternative would name a scale the object does not support.
- **The sign is not named.** Reasons in C3.10 (implementer and verifier reached the same decision independently).
- **The task's text vs the code:** the task says the label is wrong "for counts as well". It is also printed for continuous outcomes, where the CATE is a mean difference of `Y` (or `-Y`). Left unchanged, as the task says; open item 1.
- **Tests beyond the minimum:** a binary unchanged-case test (`:185-201`) and premise checks before each note.
- **Unaffected case beyond the task:** a survival `forestsearch()` fit (FS path) was added by both the implementer and the verifier.
- Site: "near 917" is the format string at start `:917`, in the block `:912-921`. Confirmed.
- No departure by the main session.

### C3.12 Open items

1. The continuous note still says `(risk-difference scale)`; its CATE is a difference in means of `Y` (or `-Y`). Kept by the task's rule. A1's evidence shows it beside the MD column: `H = -0.5838` next to MD 0.558 (`S/evidence/A1/grf_footnote.txt`).
2. Orientation of the printed CATE values, all outcome types: with `adverse_outcome = TRUE` they are differences on the negated outcome (`1 - Y`, `-Y`, `-log(Y + 0.5)`, `-Y` for counts since `f2e8de47`), and no note says so. A1's open item adds: the note's sign is opposite to the effect column's "> 0 = harm"; before A1 this held only for `adverse_outcome = TRUE`, after A1 for both orientations. The object does not record whether the count flip was applied; a field such as `grf_y_negated` would let a display state it.
3. Frontier selection prints `Hc = NA` (`cate_sgc = NA_real_`, start `R/grf_subg_harm_glm.R:576`). Pre-existing, unchanged.
4. Comment wording, not acted on (polish check notes): test header `:20-21` omits that C3's own clause reads `offset.name`; R comment `:954` uses the internal label "part 4's fix F10".

### C3.13 NEWS entry

`fb80dd8f:NEWS.md:460-464`, under "Behaviour changes" (`:154`):

> `sg_tables()` on a `grf.subg.harm.glm()` result with a count outcome no longer labels the GRF CATE note "(risk-difference scale)". It names the scale the forest was fitted on, "(log(Y + 0.5) scale, Y = count)" or "(count scale)", adding ", not offset-adjusted" when the fit has an offset. Binary and continuous results print the same note as before.

| Claim | Check | Result |
|---|---|---|
| A count result no longer gets "(risk-difference scale)" | head `:956-970`: every count path sets another label | Exact |
| The two labels, verbatim | head `:958-959` | Exact |
| ", not offset-adjusted" added when the fit has an offset | head `:964-966`, on `fs$offset.name` not NULL | Exact |
| "the scale the forest was fitted on" | Exact for the transform. Under `adverse_outcome = TRUE` the forest was fitted on `-log(Y + 0.5)` or `-Y`; the note names the unsigned scale (the recorded decision; open item 2). | Borderline, not wrong |
| Binary and continuous print the same note as before | `cate_scale` stays `"risk-difference scale"` for non-count types; unaffected case identical | Exact |
| Not mentioned | `(scale not recorded)` for a count object without a recognised `grf_y_transform`. Every `grf.subg.harm.glm()` result records it, so the bullet's scope ("a `grf.subg.harm.glm()` result") is right. | Omission, harmless |

### C3.14 Evidence notes

Contradictions between evidence files, or with git:

- **R file line counts.** NOTES §2: "+27/-2". Git: `git diff --numstat 9ec67b7e c3323a59` gives +25/-2 (27 is `--stat`'s changed-line total). Raised as verifier minor 1; NOTES §2 was not corrected.
- **Verifier start wall.** VERIFY.md §3: "run start 9.63 s". `verify/run_unaff_start.txt`: `WALL(time) 9.46 s` (in-script 9.29 s). No file records 9.63 s. The start2 (9.94 s) and fixed (9.79 s) walls match their files.
- **Version boundary of the count flip.** NOTES §1 ("released in 0.4.1.9010", "0.4.1.9009 or earlier") and `workflow1_result.json` C3 `impl.open_items[2]` ("objects from 0.4.1.9009 or earlier ... were not negated"). Git: `f2e8de47` has `Version: 0.4.1.9009` and already negates counts; `88302212` bumps to 0.4.1.9010. Superseded by the NOTES polish section and the R comment fix.
- **What `offset.name` does.** NOTES §1 table: "used only by the GLM effect estimates". Code: it also sets the complete-case rows given to the forest (start `:372-373`). Superseded by the polish section and the test header fix.

Other notes (no contradiction):

- Test line numbers in NOTES §3-§4 and VERIFY.md (`:123`, `:127`, `:138`, `:151`, `:167`, `:178`; literals `:125`, `:140`, `:153`, `:169`, `:178`, `:196`; header `:1-48`) are those of the verified 197-line file. At the commit and head they are +4 after line 17 (`:127`, `:131`, `:142`, `:155`, `:171`, `:182`; literals `:129`, `:144`, `:157`, `:173`, `:182`, `:200`; header `:1-52`), as the NOTES polish section says.
- The A1 + C3 coexistence checks (`evidence/C3/combo_apply.txt`, `verify/combo_apply.txt`) found C3's hunk at offset +10 after A1. That was the verified A1 patch (`patches/A1.verified.patch`, +9 and +1 lines above the site). The committed A1 (`f4c7dd2d`, final `patches/A1.patch`) adds +11 and +1, so C3's hunk sits at +12 in the commit (`-923` vs `-911`). The commit applied cleanly and its changed lines equal the patch's.
- NOTES §5's `compare.R` wall (0.17 s) is not recorded in any file under `evidence/C3/unaffected/`.
- NOTES §6 quotes the walls of `existing_counts.txt` (5.44, 4.40, 5.88, 9.33 s); `existing.txt` holds a second run of the same files (5.68, 4.54, 5.79, 8.93 s).
- The polish check has no `VERIFY.md`; its verdict and summary are in `S/polish_result.json` (C3 `check`), its runs in `evidence/C3/verify_polish/`.

## A4 — sweep of `make_effect_estimator()` call sites (detail)

- **Line numbers.** They are at the start SHA `b06de4e8` and at head `fb80dd8f`. `R/` is identical at `835fac1b` and `fb80dd8f`.
- **Paths.** Scratch paths are relative to `$S = /tmp/cran_prep_displays_se_docs_2026-11-12`.
- **Re-check.** Every file:line below was re-read on 2026-11-12 with `git show b06de4e8:<file>`, and every head line with `git show fb80dd8f:<file>`. Where the evidence's pointer differs from the code, the code's line is given, and A4.11 lists the difference.

### A4.1 Status and provenance

- **Status:** NO_CHANGE. Every display site that builds a `forestsearch()` fit's effect without the orientation is an A1, A2 or A3 site, and part 4's F12 had already fixed three more. No further site qualifies.
- **No artefacts:**
  - no patch (`$S/patches/` has no `A4*` file);
  - no `$S/main/A4/`;
  - no commit (`git log b06de4e8..fb80dd8f` has no A4 commit);
  - no NEWS entry (`NEWS.md` at `fb80dd8f` has no A4 bullet);
  - no `wt/A4` copy in round 1. The repair round made `$S/wt/A4_r2` (`mkcopy.sh`, no patch); its `git status` stayed clean.
- **Rounds:**

  | Round | Record | Outcome |
  |---|---|---|
  | Implement | `evidence/A4/NOTES.md` §1-10, `sweep.r1.md`, `trace_displays.*` | NO_CHANGE, 15 sites |
  | Verify 1 | `evidence/A4/verify/VERIFY.md` (16:11) | **PROBLEMS**: 2 blocking, 3 minor. All are in the other-routes list; the call-site table is correct |
  | Repair | `NOTES.md` "Repair round (round 2)"; `sweep.md` (16:16); `repair_route_e2e.*`, `repair_glm_calls.*` | NO_CHANGE. Only `sweep.md` changed |
  | Verify 2 | `evidence/A4/verify_r2/VERIFY.md` (16:25) | **CONFIRMED**: 0 blocking, 6 minor |

- **Workflow record.** In `workflow1_result.json`, the entry `id = "A4"` has `status` NO_CHANGE, `final_verdict` CONFIRMED, `repaired` true, `patches` [], `blocking_final` [], and 6 `minor_final`. `polish_result.json` has no A4 entry.
- **Round-2 minor problems.** `sweep.md` predates `verify_r2/VERIFY.md`, so it does not carry the round-2 minor problems. This file applies them (A4.6, A4.11).
- **Task wording** (task lines 90-94): "List every `make_effect_estimator(` call under `R/` with: file:line, the function it is in, whether it passes `adverse_outcome`, and what it is for (fitting, resampling, display, other). Fix a further site only if it is in a function that takes a `forestsearch()` fit and builds an effect for a table or plot without passing the orientation; ... Anything else is listed, not changed."

### A4.2 How the sweep was done

**Grep at the start SHA** (implementer; repeated by both verifiers):

    git grep -n "make_effect_estimator" b06de4e8 -- R/ | grep -vE '^[^:]*:[^:]*:[0-9]+:[[:space:]]*#'
    git grep -nE "\.make_lm_estimator_fast\(|\.make_glm_or_estimator_fast\(" b06de4e8 -- R/

- **Matches.** `make_effect_estimator` has 40 matches under `R/`. 13 are code: 12 calls and the definition (`R/glm_effect_estimators.R:100`). The other 27 are `#` / `#'` lines.
- **Fast twins.** There are 3 calls: `.make_lm_estimator_fast(` at `R/forestsearch_main.R:2287` and `R/run_subgroup_sims.R:336`, and `.make_glm_or_estimator_fast(` at `R/forestsearch_main.R:2297`. They are defined at `R/glm_effect_estimators.R:1083` and `:1190`.
- **No indirect route.** No call goes through `do.call`, `match.fun`, `get` or a function-valued use of the factory (both verifiers grepped for these). `.estimate_*()` is called only inside the factory (`:417`, `:424`, `:432`).
- **Enclosing functions.** They were found by parsing the files: `verify/enclosing.txt` (0.42 s) and `verify_r2/enclosing_r2.txt` (1.10 s). All 15 function names match the table, including the local `create_glm_row` (`R/plot_subgroup_results_forestplot.R:404-485`).
- **Call graph** (round 2, `verify_r2/callgraph_routes.R`, 7.43 s).
  - It took the 787 namespace functions and found the 77 that read `args_call_all`, `grp.consistency`, `df.est`, `sg.harm` or `mr_inference`. It then walked the codetools call graph.
  - Through the factory, it reaches only A1-A3 and F12 functions and their callers.

**Re-check for this report** (both SHAs):

- `git grep` at `fb80dd8f` finds the same 12 calls and 3 twin calls. There are 44 matches of the name; the 4 new ones are roxygen lines added by C4 `ba7c4ea2`. The definition is at `:118`.
- A small parser walked each call's argument list, from the opening parenthesis to the matching close. At `b06de4e8`, `adverse_outcome` is missing from 4 calls: rows 2, 3, 14 and 15 (row 8 passes `args$adverse_outcome %||% TRUE`). At `fb80dd8f` all 15 calls pass it.

### A4.3 Call sites (full table, start → head)

Copied from `evidence/A4/sweep.md` §1, re-checked at `b06de4e8`, with head lines added. All 15 rows agree with the code at `b06de4e8`.

| # | file:line, start → head | function | passes adverse_outcome? (source) | purpose | action |
|---|---|---|---|---|---|
| 1 | `R/bootstrap_dofuture_main.R:310` → `:310` (file unchanged in part 5) | `forestsearch_bootstrap_dofuture()` (`:211`) | **yes**: `args_forestsearch_call$adverse_outcome`, where `args_forestsearch_call <- fs.est$args_call_all` (`:268`), with `forestsearch()`'s rule if NULL (`:306-309`), passed at `:317` | resampling: bootstrap estimates for the bias correction. The same closure builds the displayed bootstrap table (`SG_tab_estimates_glm()`, `:714-719`) | none (already passes) |
| 2 | `R/forestsearch_cross_validation.R:1510` → `:1534` | `forestsearch_KfoldOut()` (`:1341` → `:1355`), branch with no original subgroup | start **no** (estimator default `TRUE`); head **yes** (A4.4) | display: ITT table (`SG_tab_estimates_glm()`, `:1518-1523` → `:1543-1548`) | fixed by A2 |
| 3 | `R/forestsearch_cross_validation.R:1566` → `:1591` | `forestsearch_KfoldOut()`, main GLM branch | start **no**; head **yes** (A4.4) | display: ITT, original-subgroup and K-fold tables (`:1583`, `:1602`, `:1624` → `:1609`, `:1628`, `:1650`) | fixed by A2 |
| 4 | `R/forestsearch_main.R:2263` → `:2300` | `forestsearch()` (`:1664`) | **yes**: the resolved formal (`:1840-1841`, again at `:2241-2243` → `:1877-1878`, `:2278-2280`), passed at `:2270` → `:2307` | fitting: search and consistency closure | none |
| 5 | `R/forestsearch_main.R:2287` → `:2324` (`.make_lm_estimator_fast`, alias) | `forestsearch()` | **yes**, `:2290` → `:2327` | fitting: fast continuous twin | none |
| 6 | `R/forestsearch_main.R:2297` → `:2334` (`.make_glm_or_estimator_fast`, alias) | `forestsearch()` | **yes**, `:2300` → `:2337` | fitting: fast binary-OR twin | none |
| 7 | `R/forestsearch_main.R:2748` → `:2798` | `forestsearch()`, propensity rebuild | **yes**, `:2754` → `:2804` | fitting: PS-adjusted closure | none |
| 8 | `R/frontier_cis.R:143` → `:151` | `compute_frontier_cis()` (`:106`) | start **yes, with fallback `TRUE`**: `args$adverse_outcome %\|\|% TRUE` (`:150`); head: the rule (A4.4) | display: the frontier CI table, which is the `ci_table` input of `pareto_frontier_table()` (`:152`) and `plot_pareto_frontier()` (`:78`) | fixed by A3 |
| 9 | `R/mr_estimates_table.R:211` → `:229` | `.mr_table_from_fs()` (`:174` → `:192`), called by `mr_estimates_table(boot_results = NULL)` at `:82` → `:100` | **yes**: `a$adverse_outcome` with `a <- fs$args_call_all` (`:190` → `:208`), with the rule if NULL (`:207-210` → `:225-228`), passed at `:215` → `:233` | display | none (fixed by part 4's F12) |
| 10 | `R/plot_sg_glm_outcomes.R:262` → `:277` | `plot_sg_glm_outcomes()` (`:92`) | **yes**: `fs.est$args_call_all$adverse_outcome`, with the rule if NULL (`:169-172` → `:184-187`), passed at `:268` → `:283` | display | none (F12) |
| 11 | `R/plot_subgroup_results_forestplot.R:444` → `:468` | local `create_glm_row()` (`:404-485` → `:428-509`) inside `plot_subgroup_results_forestplot()` (`:196`) | **yes**: `fs.est$args_call_all$adverse_outcome`, with the rule if NULL; `TRUE` when `fs.est = NULL`, i.e. no fit (`:316-322` → `:340-346`), passed at `:450` → `:474` | display | none (F12) |
| 12 | `R/run_subgroup_sims.R:321` → `:321` | `subgroup_glm()` (`:252`) | **yes**: its own argument `adverse_outcome` (default `TRUE`, `:256`), passed at `:327` (same at head) | other: a simulation fitter factory. Its input is a DGM sample, not a `forestsearch()` fit. Its only package caller (`:698` → `:720`) passes DGM fields | none (not a fit) |
| 13 | `R/run_subgroup_sims.R:336` → `:336` (`.make_lm_estimator_fast`, alias) | `subgroup_glm()` | **yes**, `:339` | other: same as row 12 | none |
| 14 | `R/summary_utility_functions.R:830` → `:863` | `sg_tables()` (`:757` → `:779`), GRF branch (`fs` of class `grf_glm_result`, `:811` → `:833`) | start **no**; head **yes** (A4.4) | display: table 1 of a GRF result | fixed by A1 |
| 15 | `R/summary_utility_functions.R:983` → `:1049` | `sg_tables()`, main GLM branch | start **no**; head **yes** (A4.4) | display: table 1 (ITT and subgroup estimates) | fixed by A1 |

- **Where the head lines moved from.** The other code commits did not touch these calls:
  - C2 `2984dbee` changes `R/forestsearch_main.R` only in the threshold warnings (hunks at start `:299-311` and `:2332-2352`).
  - C3 `c3323a59` changes `sg_tables()` only in the GRF CATE note (hunk at start `:923-933`).
  - B1 `8e87cf30` changes only `R/glm_effect_estimators.R`.
  - The remaining shifts come from A1-A3 themselves and from C4 `ba7c4ea2` and C5 `5073fd67`, which are comment-only by their gate.
  - The `adverse_outcome` lines quoted in rows 1, 4-7 and 9-13 have the same text at both SHAs.
- **Result.** The sites that do not pass the fit's orientation are rows 2, 3, 14 and 15 (no argument) and row 8 (fallback `TRUE`): A2, A1 and A3, all named in the task. Rows 9-11 were fixed by F12. Rows 1 and 4-7 pass the resolved value. Rows 12-13 take a DGM sample, not a fit. So no further site qualifies, and A4 makes **no change**.
- **Note on row 14.** Its input is a `grf.subg.harm.glm()` result, not a `forestsearch()` fit. It is fixed because the task names `sg_tables()`'s GRF branch under A1, not because of A4's criterion.

### A4.4 The A1, A2 and A3 rows: start → head

| Row | Start `b06de4e8` | Head `fb80dd8f` | Commit |
|---|---|---|---|
| 14 (A1, GRF branch) | `tryCatch(` `:829`; call `:830-837` without `adverse_outcome` | rule block `:856-859` (`adverse_outcome <- fs$adverse_outcome`; if NULL, `outcome_type %in% c("binary", "count")`); call `:863-871`; `adverse_outcome   = adverse_outcome,` at `:869` | `f4c7dd2d` (parent `2984dbee`) |
| 15 (A1, main branch) | `tryCatch(` `:982`; call `:983-990` without `adverse_outcome` | rule block `:1042-1045` (`args_fs$adverse_outcome`, rule if NULL); call `:1049-1057`; passed at `:1055` | `f4c7dd2d` |
| 2 (A2, no original subgroup) | call `:1510-1516` without `adverse_outcome` | rule block `:1517-1520` (`res$cv_args$adverse_outcome`; rule if NULL and `is_glm`), before the branch; call `:1534-1541`; passed at `:1540` | `a2cd36db` (parent `f4c7dd2d`) |
| 3 (A2, main branch) | call `:1566-1573` without `adverse_outcome` | the same block `:1517-1520`; call `:1591-1599`; passed at `:1597` | `a2cd36db` |
| 8 (A3) | call `:143-153`; `:150` `adverse_outcome   = args$adverse_outcome %\|\|% TRUE,` | rule block `:144-147` (`args$adverse_outcome`; rule if NULL); call `:151-161`; `:158` `adverse_outcome   = adverse_outcome,` | `fb21396f` (parent `a2cd36db`) |

- **The commits' hunks.** They are the sites the table assigns to each item:
  - `f4c7dd2d`: `@@ -825`, `-833`, `-978`, `-986`;
  - `a2cd36db`: `-1495`, `-1513`, `-1569`;
  - `fb21396f`: `-137`, `-147`.

  Round 2 found the same hunks in `$S/patches/A1.patch`, `A2.patch` and `A3.patch`.

### A4.5 Estimator layer inside `make_effect_estimator()`

The factory's own dispatch. These are not call sites of a fit. Lines are start → head; B1 and C4 moved them.

| Dispatch | Line | `adverse_outcome` |
|---|---|---|
| `.make_cox_estimator()` (survival; defined `:286` → `:304`) | `:142` → `:160` | not an argument |
| `.make_poisson_rate_estimator()` (binary IRR/IRD; defined `:737` → `:826`) | `:163` → `:181` | not an argument |
| `.make_glm_binary_estimator()` | `:174` → `:192` | passed, `:181` → `:199` |
| `.make_lm_estimator()` | `:188` → `:206` | passed, `:193` → `:211` |
| `.make_poisson_rate_estimator()` (count) | `:204` → `:222` | not an argument |

Measured on the start tree (`evidence/A4/trace_displays.txt`): the estimate on one slice with `adverse_outcome = TRUE` and with `FALSE`.

| Estimator | TRUE | FALSE |
|---|---|---|
| continuous MD | 0.484494 | -0.484494 |
| binary OR | 0.735707 | -0.735707 |
| binary RD | 0.180000 | -0.180000 |
| binary RR | 0.405465 | -0.330242 |
| count IRR | 0.257418 | 0.257418 |
| count IRD | 0.344772 | 0.344772 |
| survival HR | 0.375625 | 0.375625 |

So the orientation changes only continuous MD and binary OR, RD and RR. For count and survival fits, A1-A3 change nothing, whatever value is passed.

### A4.6 Other routes to a displayed effect (listed, not changed)

- **Where the table comes from.** It is copied from `evidence/A4/sweep.md` §2 (the round-2 version) and re-checked at `b06de4e8`.
- **What changed from `sweep.md`.** Corrected pointers are marked "(corr.)" and listed in A4.11. Two items are new: the MRCT row (round-2 minor 6) and the `:1252` note.
- **Line numbers** are at the start SHA. At head, these files are unchanged: `R/run_simulation_analysis.R`, `R/simulation_tables.R`, `R/plot_effect_distribution.R`, `R/grf_subg_harm_glm.R`, `R/fs_to_guohe.R` and `R/betaHhat_truth.R`. `R/oc_analyses.R` moved by +10 at the cited lines.

| file:line | function | what the code does | relevant to a fit's display? | action |
|---|---|---|---|---|
| `R/summary_utility_functions.R:218-240` (corr.) | `analyze_subgroup_glm()` (`:187`), via `SG_tab_estimates_glm()` (`:328`) | It has an effect only through the closure passed in. With `estimator_fn = NULL` the effect cell is `"NA"` (`:238-240`). A failed estimate also gives `"NA"` (`:235-237`). **No default estimator is built.** The `Diff` column is the raw `rate1 - rate0`: `:251` (binary, percent points) and `:255` (continuous, count) (corr.). The rule leaves it alone. | no route of its own; the orientation comes from the caller's closure (rows 1-3, 9, 14, 15) | none |
| `R/bootstrap_analysis_dofuture.R:445` | `bootstrap_results()` | `estimator_fn = NULL` means the Cox path. The GLM path always gets `estimator_fn_boot` (row 1, oriented). | no | none |
| `R/subgroup_search.R:207`; `R/subgroup_consistency_main.R:512`; `R/subgroup_consistency_helpers.R:1485`, `:1863` | `subgroup.search()`, `subgroup.consistency()`, `evaluate_subgroup_consistency()`, `evaluate_consistency_twostage()` | `estimator_fn = NULL` means the Cox path. There is no default GLM estimator. | no (fitting) | none |
| `R/consistency_resample.R:254-261` | `.consistency_glm_pieces()` (`:239`) | It flips the outcome itself from its own `adverse_outcome` argument (continuous `-y`, binary `1 - y`). Callers: `consistency_resample()` (exported, `:404`; input a data frame; default `TRUE` at `:417`), which calls it at `:432` (corr.: the call line); and `.fs_mr_pieces()` (`R/fs_mr_inference.R:67`, fallback `TRUE` at `:71`). Every spec `forestsearch()` builds carries the resolved value: `R/forestsearch_main.R:2891-2895`, `:3079-3083`, `:3900-3904` (the `TRUE` at `:3898` is the survival spec), and `R/forestsearch_helpers.R:1215-1219`, `:1737-1741`. | no: inference at fit time. The results are stored in `fs$mr_inference` and displayed as stored | none |
| `R/grf_subg_harm_glm.R:760-866` (GLM fit at `:833` / `:835`) | `create_glm_row()` (exported) | A raw GLM treatment coefficient, with no orientation argument. Its input is a data frame. It has no caller under `R/`: the forest plot's `create_glm_row` (`R/plot_subgroup_results_forestplot.R:404`, called at `:491`) is a local closure that shadows it, and the trace confirms the local one runs. | no (not a fit) | none |
| `R/grf_subg_harm_glm.R:284` | `grf.subg.harm.glm()` | GRF fitting. It orients the GRF response itself from its own `adverse_outcome` argument (default `FALSE`, `:304`). The flip block is `:445-457` (corr.): when `adverse_outcome` is TRUE, binary `1 - Y` (`:447`), continuous `-Y` (`:450`), count `-Y` (`:456`). It stores the value in its result (`:710`). It builds no displayed effect; its table display is `sg_tables()`'s GRF branch (row 14, A1). Its input is a data frame, not a `forestsearch()` fit. | no (fitting) | none |
| `R/run_simulation_analysis.R:963` | `.extract_fs_estimates_gen()`, called by `.run_fs_analysis_gen()` (`:774`; fit built by `do.call(forestsearch, fs_args)` at `:820`, passed as `fs_full` at `:843`), which the exported `run_simulation_analysis()` (`:389`) calls at `:592`, `:617`, `:679` | **Takes a `forestsearch()` fit** (`fs_full`) but reads only its membership (`fs_full$df.est$treat.recommend`, `:1049-1054`, `:1090-1092`) and, for the NA-membership message, `fs_full$sg.harm` (`:1093`) (corr.). `.run_fs_analysis_gen()` also reads `sg.harm` and `grp.consistency` diagnostics (`:864-886`). None of these is an effect. For GLM outcomes the extractor builds `hr.itt` (`:1025`), `hr.H.hat` / `hr.Hc.hat` (`:1132-1133`) and `hr.H.true` / `hr.Hc.true` (`:1188-1189`) with `.glm_effect()` (`:1005-1022`). That is the `stats::glm(y ~ treat)` coefficient (exponentiated for OR/IRR/RR), treatment versus control, with **no orientation**: the fit's stored `adverse_outcome` is never read. Displays: `summarize_simulation_results()` (`:418`) through `summarize_single_analysis()` (`R/oc_analyses.R:450`, called at `:431`; lines `:482-485`) (corr.); `format_oc_results()` (`:618`; `:702-705`); `build_estimation_table()` (`R/simulation_tables.R:584`; `:828-831`, row "theta-hat(sg-hat)"); `plot_effect_distribution()` (`R/plot_effect_distribution.R:199-201`). One convention holds within the table: the truth columns use the same `.glm_effect()`, and `ahr.*` / `cde.*` come from `compute_aor()` (`R/oc_analyses.R:276`) and `compute_cde_glm()` (`:332`), also arm 1 versus arm 0. `plot_effect_distribution()` flips only on the user's `subgroup_notation = "benefit"` (`:143`, `:184-192`). | **yes**: a fit-taking route whose GLM effect is raw. For continuous MD at the default orientation, it is opposite in sign to the fit's own estimate (A4.7) | listed, not changed |
| `R/run_simulation_analysis.R:1360` | `.extract_grf_estimates_gen()`, called by `.run_grf_analysis_gen()` (`:1226`) at `:1252` (only when the GRF function is not found, with `grf_est = NULL`) and at `:1317` | The twin of the row above: `.glm_effect_grf()` (`:1401-1417`), the same raw `glm()` coefficient, used at `:1419`, `:1507-1509` and `:1602-1604`. Its input is a `grf.subg.harm.glm()` / `grf.subg.harm.survival()` result. That result's GRF search is oriented by the `adverse_outcome` that `.run_grf_analysis_gen()` resolves with `forestsearch()`'s rule (`:1300-1305`); the displayed effect is not. | no (not a `forestsearch()` fit) | listed, not changed |
| `R/plot_sg_results.R:386` (inside `safe_cox()`, `:377-400`) | `compute_sg_hr_estimates()` (`:365-417`). `plot_sg_results()` calls it at `:245`, and `plot.forestsearch()` (`R/forestsearch_methods.R:599`) calls `plot_sg_results()` at `:622` | Cox HR and KM curves, survival only. The survival estimator has no orientation. | no orientation applies | none |
| `R/fs_to_guohe.R:338-373` (call at `:357`) | `fs_to_guohe()` | Survival only **by construction, not by a guard**. It always calls `guohe_algorithm3(outcome = "survival", orient = +1)` (`:357-362`; the maximum-HR harm subgroup) and never checks `fit$outcome_type`, although `guohe_algorithm3()` accepts `"survival"`, `"binary"` and `"continuous"` (`R/guohe_algorithm3.R:233`). It first calls `fs_export_maxeff_family()` (`:342`). That function stops on a fit that is not `sg_focus = "maxeff"` or has no family (`:140-150`), or on misaligned labels or missing dummy columns (`:154-163`), before the column check at `:345-350` (corr.). | no orientation applies (survival HR) | none |
| `R/test_hte_crump.R:388` | `test_hte_from_forestsearch()` | Uses the fit only for covariate names (`confounders.evaluated` / `confounders.candidate`) and runs `test_hte()` on the caller's `df`. The zero-CATE and constant-CATE statistics are quadratic forms (`T_n` at `:179`, `:226`), so they do not depend on the sign. The zero-ATE test is signed: `ate = mean(Y[W == 1]) - mean(Y[W == 0])` (`:256`), `z = ate / se` (`:259`), and `print.hte_test()` prints that signed `z` (`:296-298`; -1.90 on the F12 `ben` data, `verify/route_check.txt`). The result also carries `zero_ate$ate` (raw: -0.228688 there, against the fit's ITT orientation of +0.228688) and `constant_cate$ate_diff` (`:233`). All are treatment-minus-control contrasts with no orientation. | no: raw arm contrasts of a hypothesis test, which the rule leaves alone | none |
| `R/betaHhat_truth.R:151` | `.fs_region_effect()`, used by `fs_betaHhat_one()` (`:434`), `fs_betaHhat_table()` (`:495`) and `fs_attach_betaHhat()` (`:602`) | A raw region effect (Cox HR, logit OR at `:176`, or `compute_aor()`) on a DGM evaluation frame. The harm/benefit convention is a required `focus` argument, checked by `.fs_check_focus()` (defined at `:193`) (corr.). The input is a rule and a frame (or simulation results), not a fit. | no (not a fit) | none |
| `R/fs_mr_oc_summary.R:122` | `fs_mr_oc_summary()` | `orient` comes from a simulation payload's `meta$adverse_outcome` | no (not a fit) | none |
| `R/glm_effect_profile.R:79`; `R/glm_spline_fit.R:133` | `glm_effect_profile()`, `glm_cs_fit()` | Take a data frame | no (not a fit) | none |
| `R/dina_subgroup_refit.R:151` | `dina_subgroup_refit()` | Takes a `dina_subgroup` object | no (not a `forestsearch()` fit) | none |
| `R/mrct_trial_helpers.R:30` (fit at `:47`); `R/mrct_simulation.R:836-859` (new, round-2 minor 6) | `.mrct_analyze_trial()`, `.mrct_record_fit()` | They fit `forestsearch()` on a simulated trial and report Cox HRs (`survival::coxph(Surv(y_sim, event_sim) ~ treat_sim)`, `:840-841`). Survival only. | no orientation applies | none |
| `R/forestsearch_methods.R:245`, `:368`; `R/mr_values.R:121`; `R/pareto_frontier_table.R:143`; `R/plot_pareto_frontier.R:77`; `R/explain_pareto_selection.R:84`; `R/fs_declaration_table.R:58`; `R/frontier_member_flags.R:44` | `print`/`summary.forestsearch`, `mr_values()`, `pareto_frontier_table()`, `plot_pareto_frontier()`, `explain_pareto_selection()`, `fs_declaration_table()`, `frontier_member_flags()` | Show values computed at fit time by the fit's own (oriented) closure. None calls `make_effect_estimator()`; the trace confirms this. | already oriented | none |
| `forestsearch_Kfold()`, `forestsearch_tenfold()`, `fs_fdr_report()`, the bootstrap refits | re-run `forestsearch()` | Re-run `forestsearch()` with `args_call_all`, which carries `adverse_outcome`. `fs_fdr_report()` is survival only (guard `R/fs_fdr_report.R:223-228`). | no (fitting) | none |

**The `glm(` / `lm(` sweep** (repair round, `repair_glm_calls.R` / `.txt`). It parsed every code line under `R/` that calls `glm(`, `lm(`, `glm.nb(`, `glm.fit(` or `lm.fit(`, with its enclosing top-level function: 55 call lines in 27 functions (re-counted from the `.txt`). The groups:

- the estimator layer (`.estimate_*`, `.make_*_estimator`, `.fs_*_two_group`, `.ps_logistic`);
- fitting (`.consistency_glm_pieces`, `.dfbeta_glm`, the `dina.R` nuisance fits);
- DGM generation (`generate_glm_dgm`);
- inputs that are a data frame or a non-`forestsearch()` object (`create_glm_row`, the `glm_effect_profile` internals, `glm_cs_fit`, `.dina_fit_subgroup_model`, the Guo-He internals `.g3_coef` / `.ar_coef_se` with an explicit `orient`, `.fit_arms`, `.fs_region_effect`);
- `calibrate_L_eff` (a power-law fit, not an effect);
- the two simulation extractors.

Two of these take a `forestsearch()` fit: `.extract_fs_estimates_gen()` and `test_hte_from_forestsearch()` over `.fit_arms`. Only the first shows an effect without the orientation.

### A4.7 Dynamic checks (all on the start tree)

| Run | By | Copy | What it shows | Wall |
|---|---|---|---|---|
| `trace_displays.R` | implementer | private copy of `$S/start` | `make_effect_estimator()` traced while 16 exported displays and entry points ran on the F12 continuous fit (`tests/testthat/test-display-orientation.R`: data `.dor_data()` `:39-45`, the `ben` twin; fit `.dor_fit()` `:47-60`; stored `adverse_outcome = FALSE`). Without `adverse_outcome`: `sg_tables` (A1) `passed=FALSE value=TRUE`, n = 1; `forestsearch_KfoldOut(outall = TRUE)` (A2) `passed=FALSE value=TRUE`, n = 1. With `FALSE`: `plot_sg_glm_outcomes`, `plot_subgroup_results_forestplot` (n = 3, via the local `create_glm_row`), `mr_estimates_table` (via `.mr_table_from_fs`), `compute_frontier_cis`, `forestsearch_bootstrap_dofuture`, and the refits inside the bootstrap (n = 2) and `forestsearch_Kfold` (n = 2). No call: `print`, `summary`, `mr_values`, `pareto_frontier_table`, `plot_pareto_frontier`, `explain_pareto_selection`, `frontier_member_flags`. Errors: `fs_declaration_table` (the fit has no declaration field); `plot.forestsearch` (a survival-only KM/Cox display, given the continuous column as `event.name`) | 10.00 s process (fit 0.68 s, R 9.82 s) |
| `verify/route_check.R` | verifier 1 | `$S/verify/A4` | 0 factory calls for `test_hte_from_forestsearch`, `fs_family_report`, `plot_pareto_combined` and `.extract_fs_estimates_gen`. The `sg_tables` control: 1 call, `passed=FALSE value=TRUE`. Simulation extractor called directly: H = -1.094112 = hand raw, while the fit's own H is +1.094112. `test_hte` `zero_ate$ate` = -0.228688, z = -1.90 | 3.84 s (R 3.66 s) |
| `devtools::test(filter = "display-orientation")` | verifier 1 | `$S/verify/A4` | F12's tests on the start copy: 93 passes, 0 failures | 6.82 s |
| `verify/sim_rd_scale.R` | verifier 1 | `$S/verify/A4` | `hr.H.hat` labelled RD = -0.200671; hand RD = -0.049536; hand log OR = -0.200671 | 3.02 s (R 2.85 s) |
| `repair_route_e2e.R` | implementer, repair | `$S/wt/A4_r2` (no patch) | Through the real caller `.run_fs_analysis_gen()` with the F12 `ben` data and F12's `forestsearch()` arguments. The fit inside stores `adverse_outcome = FALSE`; subgroup `{age <= 47} & {biomarker_hi}`, n = 36. Route = hand raw on all five (tolerance 1e-10, set before the run): itt -0.228688, H.hat -1.094112, Hc.hat -0.0957932, H.true -1.286941, Hc.true -0.04755184. The hand oriented values are the negatives. The fit's own H is +1.094112 | 3.29 s process (route 0.70 s, R 3.10 s) |
| `repair_glm_calls.R` | implementer, repair | parse of `R/` | 55 call lines, 27 functions | 0.79 s (stated in `NOTES.md` only; not in the `.txt`) |
| `verify_r2/route_orientation_invariance.R` | verifier 2 | `$S/verify/A4_r2` | Own data (seed 2026, N = 240, three covariates), the real caller run twice. `hr.itt` = -0.168826 and `hr.H.true` = -1.216737 under both `adverse_outcome = FALSE` and `TRUE` (`identical()` TRUE), equal to hand raw (tolerance 1e-10, set before the run). The fit's own oriented ITT is +0.168826 under FALSE and -0.168826 under TRUE | 5.07 s (R 4.46 s) |
| `verify_r2/callgraph_routes.R` | verifier 2 | `$S/verify/A4_r2` | 787 namespace functions; 77 fit readers; through the factory, only A1-A3, F12 and their callers are reached | 7.43 s (R 7.26 s) |
| `enclosing.R`, `enclosing_routes.R`, `enclosing_r2.R` | verifiers | parse | enclosing function per site and per route | 0.42 s, 0.41 s, 1.10 s |

### A4.8 Verification and the repair round

**Round 1: PROBLEMS** (`verify/VERIFY.md`). The 15-row call-site table is correct, and NO_CHANGE for the call sites holds.

- **Blocking 1, a missing route.** `.extract_fs_estimates_gen()` takes a `forestsearch()` fit (`fs_full`) and builds the displayed simulation effects with a raw `glm()` coefficient, without the orientation. The verifier measured `hr.H.hat` = -1.094112 against the fit's +1.094112. It asked for the route to be listed, not changed, with the GRF twin listed as "not a fit".
- **Blocking 2, a wrong function name.** The row "`R/grf_subg_harm_glm.R:835` | `grf.subg.harm.glm()`" was wrong: `:835` is inside `create_glm_row()` (`:760-866`), and `grf.subg.harm.glm()` starts at `:284`.
- **Minor 1-3:**
  - The `plot_sg_results` row named only `plot_sg_results()`; `:386` is inside `compute_sg_hr_estimates()` > `safe_cox()`.
  - The `test_hte_from_forestsearch()` row said "no effect is shown". `print.hte_test()` prints the signed zero-ATE `z`, and the result carries `zero_ate$ate`. The "none" action still holds.
  - `fs_to_guohe()` has no `outcome_type` guard; it is survival only by construction.

**Repair** (`NOTES.md` "Repair round (round 2)"). Both blocking problems were accepted as correct. Only `sweep.md` changed; the round-1 copy is kept as `sweep.r1.md`. No code change and no patch.

- The `.extract_fs_estimates_gen()` row and the `.extract_grf_estimates_gen()` row were added.
- The `create_glm_row()` row now reads `:760-866` (GLM fit at `:833` / `:835`). The `grf.subg.harm.glm()` row now points to `:284` and describes its own orientation.
- The three minor rows were rewritten.
- Two runs were added: `repair_route_e2e.R` (through the real caller) and `repair_glm_calls.R` (the `glm(`/`lm(` sweep, which added the `.fs_region_effect()` row).

**Round 2: CONFIRMED** (`verify_r2/VERIFY.md`). No missing row, no wrong function name, no wrong purpose, no wrong "passes" cell. Point 3's list has no missing orientation-relevant route.

**Round-2 minor problems.** None was repaired in `sweep.md`; this file applies all six:

| # | Problem (verifier) | Re-check at `b06de4e8` | Applied here |
|---|---|---|---|
| 1 | `R/oc_analyses.R:482-485` is credited to `summarize_simulation_results()`; the lines are in `summarize_single_analysis()` (`:450`), which `summarize_simulation_results()` calls at `:431` | confirmed: `summarize_simulation_results` defined `:418`, call `:431`, `summarize_single_analysis` `:450` | A4.6 simulation row |
| 2 | "reads only its membership" is too narrow: `:1093` reads `fs_full$sg.harm` for the NA message, and `.run_fs_analysis_gen()` reads `sg.harm` / `grp.consistency` at `:864-886`. None is an effect | confirmed (`:1093` `sg_def_txt <- fs_full$sg.harm`; `:864-865`) | A4.6 simulation row |
| 3 | The `grf.subg.harm.glm()` flip is `:445-457` (count branch `:452-456`), not `:445-451` | confirmed | A4.6 GRF row |
| 4 | `R/betaHhat_truth.R:193` is cited as "the call site" of the `focus` convention; it is the definition of `.fs_check_focus()` | confirmed (`:193` `.fs_check_focus <- function(focus) {`) | A4.6 `.fs_region_effect` row |
| 5 | The open item says `fs_to_guohe()` on a GLM fit "stops only" when the named columns are absent (`:345-350`); it also stops earlier, in `fs_export_maxeff_family()` | confirmed: `:342` calls it, and it stops at `:140-150` and `:154-163` | A4.6 `fs_to_guohe` row; A4.10 |
| 6 | The MRCT survival routes (`.mrct_analyze_trial()`, `.mrct_record_fit()`) are not listed | confirmed (`R/mrct_trial_helpers.R:30`, fit `:47`; `R/mrct_simulation.R:836-859`, Cox) | A4.6 new row |

### A4.9 Differences between the task's text and the code

- **A1 pointers.** The task gives "near 829-839 and 982-992". The `tryCatch(` lines are `:829` and `:982`, and the calls `:830` and `:983`. Head: `:863` and `:1049`.
- **A2 pointers.** The task gives "near 1509-1517 and 1565-1575". The calls are `:1510` and `:1566`. Head: `:1534` and `:1591`.
- **A3 pointer.** The task gives "near 150". The call is `:143` and the fallback `:150`. Head: the call is `:151`, the rule block `:144-147`, and the argument `:158`.
- **"Point 2" / "point 3" and the quotes in `NOTES.md` §8.** The task's A4 bullet (lines 90-94) has no numbered points.
  - The evidence's "point 2" is the fix rule, and its "point 3" is the other-routes list.
  - `NOTES.md` §8 quotes as "the task" two phrases: "a default estimator used when the closure is NULL, or `SG_tab_estimates_glm()`/`analyze_subgroup_glm()` building their own", and F12's pattern "`R/plot_sg_glm_outcomes.R:166-172` and its call `:268`".
  - `git show b06de4e8:dev/tasks/TASK_cran_prep_displays_se_docs_2026-11-12.md` contains none of these.
  - They came from the run's item prompt, which is not saved under `$S`.
- **The code against those quotes:**
  - No default estimator exists. A NULL closure gives `"NA"` in `analyze_subgroup_glm()` or takes the Cox path elsewhere.
  - At `b06de4e8`, F12's call in `plot_sg_glm_outcomes()` opens at `:262`, and `:268` is its `adverse_outcome =` line. The rule block is `:166-172`, comment included (`:169-172` is the code).

### A4.10 Observations and open items (not caused by A4; listed, not changed)

- **Simulation route orientation.** `.extract_fs_estimates_gen()` and `.extract_grf_estimates_gen()` show raw treatment-minus-control GLM effects, whatever the fit's stored `adverse_outcome`.
  - For continuous MD (and binary OR/RD/RR) at a non-default orientation, the simulation tables and plots carry the opposite sign (or the reciprocal) of the fit's own estimate.
  - The tables are internally consistent. Any change would be a convention decision for the simulation layer, not a display fix.
- **Simulation scale labels** (`R/run_simulation_analysis.R:995-1003`, `:1392-1399`):
  - `"RD"` maps to `binomial()` (logit) and is not exponentiated, so the value labelled RD is the log OR.
  - `"IRD"` maps to `poisson()` and is not exponentiated, so it is the log IRR.
  - `"RR"` falls to the `binomial()` default and is exponentiated, so it is the OR.
  - `compute_aor()` returns `m1 - m0` for RD/IRD and `m1 / m0` for RR/IRR, so the `ahr.*` and `hr.*` columns disagree for these measures.
  - Measured: RD `hr.H.hat` = -0.200671 = the hand log OR; the hand RD is -0.049536.
- **`fs_to_guohe()`** has no `outcome_type` guard.
  - On a GLM fit, it stops if `fs_export_maxeff_family()` stops (not `sg_focus = "maxeff"`, no family, misaligned labels) or if the named `time` / `event` / `treatment` columns are absent.
  - Otherwise it runs a survival analysis on those columns.
- **Propensity adjustment.** The displays do not carry `ps_adjust_method`. This is part 4's open item, and it also holds at the A1 and A2 sites, which use the default `"none"`.
- **Run hygiene.** At 16:00:33 an empty directory `$S/start/tests/testthat/_snaps/` appeared. It is not in `b06de4e8`, A4 did not create it, and no file content changed. Recorded in `NOTES.md` §9 and by both verifiers.

### A4.11 Evidence against git: differences found

| Evidence says | Git at `b06de4e8` says |
|---|---|
| `sweep.md` §2, `NOTES.md` §8, `workflow1_result.json` A4 `impl.open_items`: NULL closure gives `"NA"` at `:218-238`; raw `Diff` at `:245-253`. Verifier 2 confirms `:245-253`; verifier 1 gave `:245-255` | The NULL-closure `"NA"` is assigned at `:239` (`else` block `:238-240`). The raw Diff is at `:251` (binary) and `:255` (continuous and count), so `:245-253` misses `:255` |
| `sweep.md` §2 simulation row: `R/oc_analyses.R:482-485` belongs to `summarize_simulation_results()` | The lines are in `summarize_single_analysis()` (`:450`), called at `:431` by `summarize_simulation_results()` (`:418`) (verifier 2, minor 1) |
| `sweep.md` §2, `NOTES.md` repair, and the `workflow1_result.json` repair open item: the extractor "reads only its membership" | `:1093` also reads `fs_full$sg.harm`, for the NA message (not an effect) (verifier 2, minor 2) |
| `sweep.md` §2, `NOTES.md` repair, `workflow1_result.json` `blocking_r1`: the `grf.subg.harm.glm()` flip is at `:445-451` | The flip block is `:445-457`; the count branch is `:452-456` (verifier 2, minor 3) |
| `sweep.md` §2: `R/betaHhat_truth.R:193` is the `focus` "call site" | `:193` is `.fs_check_focus <- function(focus) {`, a definition (verifier 2, minor 4) |
| `NOTES.md` repair open items and `workflow1_result.json` repair open item: `fs_to_guohe()` "on a GLM fit stops only when the named columns are absent (`:345-350`)" | `:342` first calls `fs_export_maxeff_family()`, which stops at `:140-150` / `:154-163` (verifier 2, minor 5) |
| `NOTES.md` §8: the cox and poisson dispatch calls are at `R/glm_effect_estimators.R:142-152`, `:163-172`, `:204-213` | The cox call is `:142-149` (`:150-152` are `},`, a blank line and `binary = {`). The poisson calls `:163-172` and `:204-213` are exact |
| `NOTES.md` §7: the fit is "the F12 continuous fixture (`tests/testthat/test-display-orientation.R:39-58`)" | `.dor_data()` is `:39-45` and `.dor_fit()` is `:47-60`, so `:39-58` leaves out the last two lines of `.dor_fit()` (`sweep.md` §4 cites `:39-45` for the data, which is exact) |
| `NOTES.md` §8 quotes two passages and the F12 pointer `:166-172` / `:268` as "the task" | The task document at `b06de4e8` contains none of them (A4.9) |
| `NOTES.md` repair and `workflow1_result.json` repair: `repair_glm_calls.R` wall 0.79 s | `repair_glm_calls.txt` records no wall line, so the figure cannot be checked against an output file (not a contradiction) |

**Other pointers re-checked, with no difference:**

- every row of A4.3;
- the dispatch lines of A4.5;
- in A4.6: `:445`, `:512`, `:1485`, `:1863`, `:254-261`, `:404`, `:417`, `R/fs_mr_inference.R:71`, the spec lines, `:284`, `:304`, `:710`, `:760-866`, `:833`, `:835`;
- in `R/run_simulation_analysis.R`: `:389`, `:592`, `:617`, `:679`, `:774`, `:820`, `:843`, `:963`, `:995-1003`, `:1005-1022`, `:1025`, `:1049-1054`, `:1132-1133`, `:1188-1189`, `:1226`, `:1252`, `:1300-1305`, `:1317`, `:1360`, `:1392-1399`, `:1401-1417`, `:1419`, `:1507-1509`, `:1602-1604`;
- `R/oc_analyses.R:276`, `:332`, `:702-705`; `R/simulation_tables.R:828-831`; `R/plot_effect_distribution.R:143`, `:184-192`, `:199-201`;
- `R/plot_sg_results.R:245`, `:365`, `:377`, `:386`, `:400`, `:417`; `R/forestsearch_methods.R:245`, `:368`, `:599`, `:622`;
- `R/fs_to_guohe.R:338`, `:345-350`, `:357-362`, `:373`; `R/guohe_algorithm3.R:233`;
- `R/test_hte_crump.R:179`, `:226`, `:233`, `:256`, `:259`, `:296-298`, `:388`;
- `R/betaHhat_truth.R:151`, `:176`, `:434`, `:495`, `:602`; `R/fs_mr_oc_summary.R:122`;
- the stored-value display definitions; `R/fs_fdr_report.R:223-228`.

## B2 — rate-difference branch of `.make_poisson_rate_estimator()` (`dr_gcomp`, IRD)

Line numbers are at the start SHA `b06de4e8` and at head `fb80dd8f`. `R/` is identical at `835fac1b` and `fb80dd8f`. Every quoted line was re-read on 2026-11-12 with `git show b06de4e8:R/glm_effect_estimators.R`, and the head lines with `git show fb80dd8f:...`.

### B2.1 Status and provenance

- **Status:** BLOCKED. Task Part B, B2 step 2 sets a condition, and it is not met. Step 3 applied: nothing was changed.
- **No artefacts:** no patch (`$S/patches/` has no `B2*.patch`), no `$S/main/B2/` directory, no commit, no NEWS entry. `git diff b06de4e8 fb80dd8f -- NEWS.md` adds no line on `dr_gcomp` or on the rate-difference SE.
- **Verifier:** CONFIRMED (round 1). Blocking problems: none. Minor problems: one (B2.9). No repair round. `polish_result.json` has no B2 entry.
- **Workflow record** (`workflow1_result.json`, entry `id = "B2"`): `status` BLOCKED, `final_verdict` CONFIRMED, `repaired` false, `patches` [], `blocking_final` [], one `minor_final`.
- **Copies used:**
  - Implementer: `$S/wt/B2/`, unmodified; `git status --porcelain` empty.
  - Verifier: `$S/verify/B2/start_ctrl/`. Four R files were `cmp`-identical to `b06de4e8`, and the copy stayed clean.
- The repository was only read.

### B2.2 Same code at both SHAs

- `git show b06de4e8:R/glm_effect_estimators.R | sed -n '724,927p'` and `git show fb80dd8f:R/glm_effect_estimators.R | sed -n '813,1016p'` are byte-identical (`diff` empty).
- The 89-line shift comes from two part-5 commits. All of their hunks in this file are above line 611 at head:

| Commit | Item | Lines in `R/glm_effect_estimators.R` |
|---|---|---|
| `8e87cf30` | B1 | +75/-4 |
| `ba7c4ea2` | C4 | +23/-5 |

- `.estimate_or()` (`:585-648` at start, `:674-737` at head) and the propensity helpers (`:1257-1545` at start, `:1346-1634` at head) are byte-identical too.
- `R/forestsearch_main.R:2678-2758` (start) = `:2728-2808` (head), byte-identical, offset +50.
- `R/bootstrap_analysis_dofuture.R` is unchanged (`git diff --stat` empty).
- The `cv_args$ps_hat <- NULL` lines in `R/forestsearch_cross_validation.R` are at `:398` and `:901` at both SHAs.

### B2.3 Line map (start → head)

| What | Start `b06de4e8` | Head `fb80dd8f` |
|---|---|---|
| Binary routing to the rate estimator (`effect_measure %in% c("IRR", "IRD")`) | `glm_effect_estimators.R:162-172` | `:180-190` |
| Count routing | `:198-213` | `:216-231` |
| `.estimate_rd()` logistic tier, adjusted SE (B1 site) | `:546` | `:576` (fallback only) |
| `.estimate_or()` `dr_gcomp` branch | `:602-628` | `:691-717` |
| `.make_poisson_rate_estimator()` roxygen | `:724-735` | `:813-824` |
| `.make_poisson_rate_estimator()` | `:737-927` | `:826-1016` |
| IRR existence check | `:759-772` | `:848-861` |
| Offset floor and `.log_offset` | `:775-781` | `:864-870` |
| `iptw` branch | `:787-806` | `:876-895` |
| `dr_gcomp` condition | `:812-814` | `:901-903` |
| `dr_gcomp` outcome fit | `:815-823` | `:904-912` |
| `ps`, `nd1`, `nd0` | `:827-832` | `:916-921` |
| `pred1`, `pred0` | `:834-835` | `:923-924` |
| `dr_gcomp` IRR | `:837-848` | `:926-937` |
| **`dr_gcomp` IRD (B2)** | **`:849-862`** | **`:938-951`** |
| Unadjusted branch (IRR `:884-894`, IRD `:895-913`) | `:867-913` | `:956-1002` |
| Error handler | `:916-923` | `:1005-1012` |
| `.robust_or_model_se()` | `:1259-1272` | `:1348-1361` |
| `estimate_propensity_scores()` `@param trim` ("Quantile bounds") | `:1357` | `:1446` |
| PS `switch()` (grf / lasso / logistic) | `:1419-1423` | `:1508-1512` |
| Clamp to `trim` | `:1426-1434` | `:1515-1523` |
| Bang-Robins comment | `:1445-1451` | `:1534-1540` |
| `ips_covar <- ifelse(...)` | `:1459` | `:1548` |
| Columns written (`ps_hat`, `sw`, `ips_covar`) | `:1461-1463` | `:1550-1552` |
| `.ps_grf()`, `.ps_lasso()`, `.ps_logistic()` | `:1482-1494`, `:1502-1518`, `:1526-1533` | `:1571-1583`, `:1591-1607`, `:1615-1622` |
| `forestsearch()` `@param effect_measure`, binary list includes `"IRD"` | `forestsearch_main.R:1090-1091` | `:1112-1113` |
| `forestsearch()` `@param ps_adjust_method` | `:1160-1166` | `:1182-1188` |
| `ps_adjust_method = c("none", "iptw", "dr_gcomp")` | `:1758` | `:1795` |
| IRR/IRD require `offset.name` | `:2254-2260` | `:2291-2297` |
| `ps_adjust_resolved <- if (ps_method == "none") "none" else ps_adjust_method` | `:2683` | `:2733` |
| User-supplied `ps_hat` (checks, then `df.analysis$ps_hat <- ps_hat`) | `:2700-2720` | `:2750-2770` |
| User path `ips_covar` | `:2726-2728` | `:2776-2778` |
| `estimate_propensity_scores()` call | `:2731-2738` | `:2781-2788` |
| Closure rebuilt with `ps_adjust_method = ps_adjust_resolved` (GLM only) | `:2747-2756` | `:2797-2806` |
| CV nulls `ps_hat` | `forestsearch_cross_validation.R:398`, `:901` | same |
| Bootstrap nulls `ps_hat` | `bootstrap_analysis_dofuture.R:593` | same |
| `trim` help text | `man/estimate_propensity_scores.Rd:33` | same |

### B2.4 The code, quoted

`R/glm_effect_estimators.R` at `b06de4e8`. Add 89 for head. Each line is verbatim; only the leading indentation is trimmed.

```
:775  time_vec <- data_slice[[offset.name]]
:776  if (any(time_vec <= 0, na.rm = TRUE)) {
:777    time_vec <- pmax(time_vec, .Machine$double.eps)
:781  data_slice$.log_offset <- log(time_vec)
:809  # Option B: Bang & Robins DR -- IPS covariate + G-computation
:810  # (Section 2.2, Bang & Robins 2005)
:812  } else if (ps_adjust_method == "dr_gcomp" &&
:813             "ips_covar" %in% names(data_slice) &&
:814             "ps_hat" %in% names(data_slice)) {
:815    fmla <- stats::as.formula(
:816      paste0(outcome.name, " ~ ", treat.name,
:817             " + ips_covar + offset(.log_offset)")
:819    fit <- stats::glm(
:820      fmla,
:821      data   = data_slice,
:822      family = stats::poisson(link = "log")
:825    # G-computation: predict under both treatment scenarios
:826    # Under W=1: ips = 1/e(x);  Under W=0: ips = 1/(1-e(x))
:827    ps <- data_slice$ps_hat
:828    nd1 <- nd0 <- data_slice
:829    nd1[[treat.name]] <- 1L
:830    nd1$ips_covar     <- 1 / ps
:831    nd0[[treat.name]] <- 0L
:832    nd0$ips_covar     <- 1 / (1 - ps)
:834    pred1 <- stats::predict(fit, newdata = nd1, type = "response")
:835    pred0 <- stats::predict(fit, newdata = nd0, type = "response")
:849    } else {
:850      # IRD
:851      est <- mean(pred1) - mean(pred0)
:852      se_treat <- sqrt(diag(stats::vcov(fit)))[[treat.name]]
:853      avg_deriv <- mean(pred1)
:854      list(
:855        estimate    = est,
:856        se          = abs(avg_deriv) * se_treat,
:857        converged   = fit$converged,
:858        n0 = n0, n1 = n1,
:859        measure     = "IRD",
:860        method_used = "poisson_dr_gcomp"
```

- `.log_offset` is added to `data_slice` at `:781`, before the copies at `:828`. So `nd1` and `nd0` carry the offset, and `predict(type = "response")` returns `t_i * exp(linear predictor)`.

### B2.5 The estimate as a formula, and its inputs

- `b = (b0, bT, bI)` are the coefficients of the single fit `glm(y ~ treat + ips_covar + offset(.log_offset), family = poisson(link = "log"))` on the slice (`:815-823`).
- `t_i` is the slice's follow-up time after the floor.
- `ps_i` is the slice's `ps_hat` column.
- `n` is the number of rows of the slice. All rows of both arms are predicted under both treatments.

```
pred1_i = t_i * exp(b0 + bT + bI / ps_i)          (:829-830, :834)
pred0_i = t_i * exp(b0      + bI / (1 - ps_i))    (:831-832, :835)
est     = (1/n) sum_i pred1_i - (1/n) sum_i pred0_i   (:851)
```

| Input | Computed by | From what | When |
|---|---|---|---|
| `b0`, `bT`, `bI` | `stats::glm()` (`:819-823`) | the slice: outcome, treatment, `ips_covar`, `.log_offset` | each closure call, so each slice |
| `t_i` (`.log_offset = log(t_i)`) | the closure (`:775-781`) | the slice's `offset.name` column. This is data. Values `<= 0` are floored at `.Machine$double.eps`. | each closure call |
| `ips_covar` (a regressor in the fit) | `estimate_propensity_scores()` `:1459`: `ifelse(W == 1, 1 / ps_hat, 1 / (1 - ps_hat))`. For a user-supplied `ps_hat`, `forestsearch()` at `R/forestsearch_main.R:2726-2728`. | `ps_hat` | once, on the full analysis data, before the search. The column is carried into every slice. |
| `ps_i` (`ps_hat`, used in the prediction as `1/ps`, `1/(1-ps)`, `:827-832`) | `estimate_propensity_scores()` `:1419-1423`, then clamped to `trim` at `:1426-1434` and written at `:1461` | a **propensity model** fitted to the full analysis data (see the list below), or the user's `ps_hat` argument (`R/forestsearch_main.R:2700-2720`). The user's values are not clamped; values outside (0, 1) are an error (`:2712-2719`). | once, before the search (`R/forestsearch_main.R:2731-2738`). Not refitted within a slice. Refitted per CV training fold (`cv_args$ps_hat <- NULL`, `R/forestsearch_cross_validation.R:398`, `:901`) and per bootstrap replicate (`args_FS_boot$ps_hat <- NULL`, `R/bootstrap_analysis_dofuture.R:593`). |

- The propensity model is one of:
  - `.ps_grf()` (`:1482-1494`): `grf::regression_forest(X, W)` with OOB predictions, or a supplied forest's `W.hat`. `forestsearch()` passes no forest at `:2731-2737`.
  - `.ps_lasso()` (`:1502-1518`): `glmnet::cv.glmnet`, binomial, `lambda.min`.
  - `.ps_logistic()` (`:1526-1533`): `glm(W ~ confounders, binomial)`.
- **When the branch runs:**
  - Direct call: `make_effect_estimator()` with `effect_measure = "IRD"` and `ps_adjust_method = "dr_gcomp"`, for `outcome_type = "count"` (`:198-213`) or `"binary"` (`:162-172`). The slice must carry both `ips_covar` and `ps_hat`.
  - Through `forestsearch()`: a GLM outcome, `ps_method != "none"` (`ps_adjust_resolved`, `:2683`) and `ps_adjust_method = "dr_gcomp"`. The default is `"none"` (`:1758`). The help calls the option implemented but not tested against the manuscript analyses (`:1160-1163`). The closure is rebuilt at `:2747-2756`.
  - `forestsearch()` documents `"IRD"` for binary outcomes (`:1090-1091`) and requires `offset.name` for it (`:2254-2260`).
  - Without the two columns the call falls through to the unadjusted branch: measured `method_used = "poisson_offset_delta"` (implementer F4, verifier 6e).

### B2.6 The current standard error

```
se = |mean(pred1)| * sqrt(vcov(fit)[treat, treat])        (:852-856)
```

- The gradient of `est` in `b`, derived by hand from B2.5:
  - `d est / d b0 = est`
  - `d est / d bT = mean(pred1)`
  - `d est / d bI = mean(pred1 / ps) - mean(pred0 / (1 - ps))`
- The implementer checked it: analytic = central finite differences, h = 1e-6, (0.326687, 0.854623, 0.410181). The verifier checked it on its own data: (0.682937, 1.63775, 1.92436), and the first component equals `est`.
- The current SE is therefore `|d est / d bT| * se(bT)`, the `bT` term of the delta method alone. It drops:
  - the `b0` and `bI` terms and every covariance between coefficients;
  - the variability of `ps_hat`, which enters both the fit (through `ips_covar`) and the prediction.
- It is model-based. `stats::vcov(fit)` is used whatever `robust_se` is (default `TRUE`); `.robust_or_model_se()` is not called. Verifier 6b: the IRD SE is `identical()` with `robust_se` `TRUE` and `FALSE`.
- Consistency check across the evidence:
  - Implementer: 0.8546228 x `se(bT)` 0.1282556317 (the F1 IRR SE on the same data) = 0.109610, the closure SE.
  - Verifier: 1.63775 x 0.07219603387 (its 6a) = 0.118239, its closure SE.

### B2.7 Why the pre-set condition fails

The condition (task B2 step 2) is: "Only if the estimate is a smooth function of the coefficients of one fitted outcome model, with nothing estimated elsewhere entering it (no propensity model, no augmentation term)". Both exclusions apply, whether `ps_hat` is estimated or supplied by the user.

1. **A propensity model enters.** `ps_hat` comes from a separate fit, or from the user. It enters the prediction directly (`:830`, `:832`). With `b` fixed:
   - `d est / d ps_i = (1/n) * (-pred1_i * bI / ps_i^2 - pred0_i * bI / (1 - ps_i)^2)`. This is non-zero for every subject.
   - Implementer: range [0.0002572, 0.003724], n = 400.
   - Verifier: central differences, h = 1e-4, range [0.0001213, 0.03743], non-zero for 600 of 600.
   - Same outcome data, only the propensity model changed:
     - implementer: est 0.326687 (logistic x1 + x2), 0.327447 (logistic x1), 0.325843 (logistic x1 + x2, `trim = c(0.30, 0.70)`, 65 clamped);
     - verifier: 0.68293652 (z1 + z2 + z3) vs 0.69648025 (z1).
2. **An augmentation term enters.** `ips_covar` is the Bang-Robins inverse-propensity covariate (`:809-810`, `:825-826`, `:1445-1451`) and a regressor of the fit. So `b`, `bI` included, is itself a function of `ps_hat`. Verifier: `bI` = -0.128067, z = -2.965.
3. **Within a slice, `ps_hat` is the full-data score.** It is carried, not refitted. Implementer E, slice `x2 == 1`, n = 195: closure 0.2620550034 = by hand with the slice's carried column.

- `sqrt(g' V g)` with `V` from the outcome fit alone would treat `ps_hat` as known. The task's rule excludes that case, so no such change was made.
- Values for information only, not proposed:

| | Implementer (n = 400, logistic PS) | Verifier (n = 600, logistic PS, three confounders) |
|---|---|---|
| Current closure SE | 0.109610 | 0.118239 |
| `sqrt(g' V g)`, outcome fit only, `ps_hat` fixed | 0.082362 | 0.0940735 |
| Bootstrap SD, rows resampled, `ps_hat`/`ips_covar` carried | not run | 0.124714 (B = 300, seed 7) |
| Bootstrap SD, whole two-stage procedure (PS refitted per replicate) | not run | 0.122884 (B = 300, seed 7) |

### B2.8 Measurements (start tree)

| Check | Implementer (`evidence/B2/form_check.R`) | Verifier (`verify/v_form.R`, `v_se_addendum.R`) |
|---|---|---|
| Data | synthetic, n = 400, seed 20261112, PS `estimate_propensity_scores(method = "logistic")` on x1 + x2 | own design, n = 600, three confounders, seed 20261112, logistic PS; written without the implementer's script |
| Closure | `poisson_dr_gcomp`, `IRD`, est 0.3266870139, se 0.1096101930 | `poisson_dr_gcomp`, `IRD`, est 0.682936519003, se 0.118238831259 |
| Estimate by hand (written-out formula, separate fit) | 0.3266870139, `all.equal` (tol 1e-10) TRUE | 0.682936519003 (`glm.fit` on `cbind(1, trt, ips_covar)`), `all.equal` TRUE |
| SE form `|mean(pred1)| * se(bT)` | 0.1096101930, `all.equal` (tol 1e-10) TRUE | with `vcov()` of a separate `glm()`: 0.118238831259, `identical()` TRUE (addendum). With `V = (X'WX)^{-1}` at the converged `b`: 0.118239060557, rel. 2e-6, `all.equal` FALSE in `v_form.txt`. With `glm(epsilon = 1e-14)`: 0.118239060555. The gap comes from `vcov.glm`'s last-iteration IRLS weights, not from the form. |
| Gradient in `(b0, bT, bI)` | (0.326687, 0.854623, 0.410181), analytic = central FD | (0.682937, 1.63775, 1.92436) |
| `dr_gcomp` IRR | est 0.481686 = `log(mean(pred1)/mean(pred0))`; `bT` 0.461615; se 0.1282556317 = `se(bT)` | est 0.53956396; `bT` 0.58291894; se 0.07219603387 = model `se(bT)` |
| `iptw` with `effect_measure = "IRD"` | `poisson_iptw`, `log-IRR`, est 0.355011 | `poisson_iptw`, `log-IRR`, est 0.334301 |
| Unadjusted IRD | `poisson_offset_delta`, est 0.304270, se 0.077614 | est 0.48325391, se 0.062978855, both = hand delta method over `(b0, bT)` with `vcov()` |
| `dr_gcomp` IRD without the PS columns | `poisson_offset_delta` | `poisson_offset_delta` |
| Binary route, IRD, `dr_gcomp` | not run | `poisson_dr_gcomp`, estimate and SE identical to the count route (6f) |
| `trim` | `c(0.30, 0.70)`: 65 of 400 clamped | `c(0.3, 0.7)`: `trimmed` 230 = raw scores outside [0.3, 0.7]; `ps_hat` range exactly 0.3 .. 0.7; raw 0.3/0.7 quantiles 0.339373 / 0.579148 |
| Adjusted unadjusted-branch IRD | code reading only | closure 0.20119603 = rate difference at covariates 0; marginal value 0.26505089 |
| Walls | 4.42 s in the script, 4.62 s for the process | `v_form`: 5.84 s / 6.01 s. Addendum: 2.23 s / 2.38 s. Implementer script re-run on the control copy: 2.6 s / 2.77 s, output identical apart from the wall lines. |

### B2.9 Branches and scaled-coefficient SE sites

**Every branch of `.make_poisson_rate_estimator()`** (start lines; add 89 for head). None was changed.

| Branch | Lines | Estimate | SE | Scales a coefficient SE? |
|---|---|---|---|---|
| IRR existence check | `:759-772` | none (`.fs_nonestimable`) | NA | no |
| `iptw` with `sw` (IRR and IRD requests) | `:787-806` | `coef(fit)[[treat]]` of the `sw`-weighted Poisson fit (`:797`) | `.robust_or_model_se()` (`:798`; `:1259-1272`) | No: the coefficient's own SE. An IRD request returns the log-IRR (open item). |
| `dr_gcomp` IRR | `:837-848` | `log(mean(pred1) / mean(pred0))` (`:838`) | model `se(bT)` (`:840`, `:843`) | **Yes, factor 1.** The estimate is not `bT`. |
| `dr_gcomp` IRD | `:849-862` | `mean(pred1) - mean(pred0)` (`:851`) | `abs(mean(pred1)) * se(bT)` (`:852-856`) | **Yes. This is B2.** |
| Unadjusted (or fall-through) IRR | `:884-894` | `coef(fit)[[treat]]` (`:885`) | `.robust_or_model_se()` (`:886`) | no |
| Unadjusted (or fall-through) IRD | `:895-913` | `lam1 - lam0`, `lam = exp(b0 [+ bT] + mean(log t))` (`:899-902`) | `sqrt(grad' V[1:2, 1:2] grad)`, `grad = c(lam1 - lam0, lam1)` (`:903-904`) | No: delta method over `(b0, bT)`, model `V` |
| Error handler | `:916-923` | NA | NA | no |

**Package-wide.** `git grep -n -E 'se_treat|avg_deriv' b06de4e8 -- R/` hits only `R/glm_effect_estimators.R` lines 533, 535, 546, 625, 626, 840, 843, 852, 853 and 856. `git grep -n '"dr_gcomp"' b06de4e8 -- R/` (roxygen lines excluded) hits `R/glm_effect_estimators.R:110` (the argument), `:602` and `:812`, plus `R/forestsearch_main.R:1758` (the argument) and `:2635` (a comment). The scaled-coefficient SE sites:

| Site | Start | Head | Estimate | SE | In this task |
|---|---|---|---|---|---|
| Count `dr_gcomp` IRD | `:849-862` | `:938-951` | G-computation rate difference | `abs(mean(pred1)) * se(bT)` | **B2: BLOCKED, not changed** |
| Count `dr_gcomp` IRR | `:837-848` | `:926-937` | marginal log rate ratio | `se(bT)` (factor 1). The comment `:839` "Delta-method SE approximation" does not describe it. | Listed, not changed |
| Binary `.estimate_or()` `dr_gcomp` | `:602-628` | `:691-717` | marginal log-OR from `mean(pred1)`, `mean(pred0)` (`:623-624`) | `se_treat` (`:625-626`, factor 1) | Listed, not changed. Same `ps_hat` dependence. Found by reading the code, not measured. |
| Binary `.estimate_rd()` logistic tier, adjusted | `:531-547` | `:549-578` | `mean(pred1) - mean(pred0)` | `abs(avg_deriv) * se_treat` at start (`:546`). At head it is only the fallback (`:576`) when `.fs_rd_gcomp_se()` (`:640`) does not return a finite, non-negative value. | Changed by B1, `8e87cf30` |

### B2.10 Verification and the verifier's problems

- **Verdict:** CONFIRMED. The verifier's work:
  - It traced the code independently and found every quoted line in NOTES §1 verbatim at `b06de4e8`.
  - It checked the input table and the branch table line by line.
  - It ran its own reproducer on its own data, and re-ran the implementer's script on a control copy.
- **Blocking:** none.
- **Minor (one):** NOTES §1 says `outcome_type = "count"` with `effect_measure = "IRD"` reaches `:849-862`. `outcome_type = "binary"` with `effect_measure = "IRD"` also routes to `.make_poisson_rate_estimator()` (`:162-172`; head `:180-190`) and reaches the same branch. Verifier 6f: `poisson_dr_gcomp`, estimate and SE identical to the count route. NOTES leaves this out.
  - Outcome: this report states both routes (B2.5). No code change was involved.
  - The binary IRD route is also documented for `forestsearch()` (`forestsearch_main.R:1090-1091`; head `:1112-1113`).
- **Open items the verifier checked:**
  - `dr_gcomp` ignores `robust_se`: measured.
  - `trim` is absolute, not a quantile: measured.
  - The adjusted unadjusted-branch IRD is evaluated at covariates 0: measured. NOTES had it from code reading only.
  - `.estimate_or()` lines: correct.
  - No `dr_gcomp` test: `git grep` exit 1.
  - The roxygen at `:724-735` describes only the unadjusted fit: correct.

### B2.11 The task's text against the code

- "near 856 ... `abs(avg_deriv) * se_treat`": exact. `:853` `avg_deriv <- mean(pred1)`, `:856` `se = abs(avg_deriv) * se_treat`.
- "Reported to use `abs(avg_deriv) * se_treat` as well": true. The B1 site used the same expression with `avg_deriv <- mean(pred1 * (1 - pred1))` (`:535`, `:546`).
- Step 2's condition is not met; step 3 applied.
- The roxygen of `.make_poisson_rate_estimator()` (`:724-735`; head `:813-824`, `@noRd`) describes only the unadjusted fit (`glm(event ~ treat, ..., offset = log(time))`). It does not describe the `dr_gcomp` or `iptw` fits. Noted, not changed (outside B2).

### B2.12 Open items B2 leaves

1. **B2 itself.** The `dr_gcomp` IRD SE stays `abs(mean(pred1)) * se(bT)` (`:852-856`; head `:941-945`). An SE that reflects this estimator needs the propensity model's variability. Two routes: joint estimating equations of the propensity and outcome fits, or resampling the whole two-stage procedure. GRF and LASSO scores have no closed-form `V`. A decision for Larry, outside this task's pre-set rule.
2. **`dr_gcomp` IRR** (`:837-848`; head `:926-937`) and **`.estimate_or()` `dr_gcomp`** (`:602-628`; head `:691-717`) are the same class: `se(bT)` is reported for a marginal G-computation contrast. The `:839` comment "Delta-method SE approximation" does not describe the code.
3. **`iptw` ignores `effect_measure`** (`:787-806`). With `effect_measure = "IRD"` and an `sw` column, the closure returns `measure = "log-IRR"`, `method_used = "poisson_iptw"`. Measured est 0.355011 (implementer), 0.334301 (verifier).
4. **Both `dr_gcomp` branches use the model-based `vcov()`** (`:840`, `:852`) whatever `robust_se` is (default `TRUE`). Verifier 6b.
5. **`trim` documentation.** `estimate_propensity_scores()` documents `trim` as "Quantile bounds for PS trimming" (`:1357`; `man/estimate_propensity_scores.Rd:33`, unchanged at head). The code clamps at the absolute values (`pmax(pmin(ps_hat, hi), lo)`, `:1427-1434`). Measured: 65 of 400 (implementer), 230 of 600 (verifier).
6. **No test covers any `dr_gcomp` branch.** `git grep -n dr_gcomp -- tests/` returns nothing at `b06de4e8` or at `fb80dd8f` (exit 1).
7. **Adjusted unadjusted-branch IRD.** With `adjust_covariates` (`:872-877`, `:897-904`), both rates are evaluated at covariates = 0 and the mean of `log(t)`. That is a conditional rate difference at the origin, not a marginal one. Verifier: 0.20119603 vs marginal 0.26505089.

### B2.13 Evidence files

- Implementer:
  - `/tmp/cran_prep_displays_se_docs_2026-11-12/evidence/B2/NOTES.md`. It has no repair, polish or correction section.
  - `/tmp/cran_prep_displays_se_docs_2026-11-12/evidence/B2/form_check.R`, `form_check.txt`.
- Verifier:
  - `/tmp/cran_prep_displays_se_docs_2026-11-12/evidence/B2/verify/VERIFY.md`
  - `v_form.R`, `v_form.txt`
  - `v_se_addendum.R`, `v_se_addendum.txt`
  - `impl_form_check_on_ctrl.R`, `impl_form_check_on_ctrl.txt`
- Workflow: `/tmp/cran_prep_displays_se_docs_2026-11-12/workflow1_result.json`, entry `B2`.

## C4 — help pages and source comments (detail)

Line numbers: "start" = `b06de4e8`, "head" = `fb80dd8f` (`R/` identical at `835fac1b`). Code-support lines are given at start, with the head line in parentheses where it differs. Files not changed between start and head keep their lines (`R/bootstrap_dofuture_main.R`, `R/forestsearch_helpers.R`, `R/consistency_resample.R`, `R/dina.R`, `R/dina_df.R`, `R/grf_helpers.R`, `R/grf_main.R`, `R/grf_subg_harm_glm.R`, `R/get_FSdata_helpers.R`, `R/calibrate_helpers.R`, `R/generate_aft_dgm_main.R`, `R/generate_aft_dgm_helpers.R`, `R/sim_aft_gbsg.R`, `R/simulate_from_dgm.R`). Evidence root `S=/tmp/cran_prep_displays_se_docs_2026-11-12`.

### C4.0 Provenance

- **Commit:** `ba7c4ea2` "docs: help pages and source comments say what the code does (part 5, C4)", parent `c3323a59` (C3). `main/C4/commit.txt`, `main/C4_pre.txt` agree.
- **Diffstat:** 35 files, +699/-96 (`main/C4/diffstat.txt` = `git show --stat ba7c4ea2`). 18 `R/` files, 17 `man/` pages. No file under `tests/`, `vignettes/`, `quarto/`.
- The commit message ends "three source comments". The commit changes four `#` blocks: the task's three (Site 14) and the `use_dina` screening comment (Site 4, C4b polish P4).
- **Three verified patches, one commit:**

| Part | Final patch (sha256) | Files | +/- | Verified patch (sha256, +/-) | Polish delta |
|---|---|---|---|---|---|
| C4a | `patches/C4a.patch` (`2cb414fe...ef2204`) | 9 R + 10 man | +405/-18 | `C4a.verified.patch` (`60ec92c5...65925dc6`), +393/-18 | `evidence/C4a/polish.diff`: 11 files, +40/-28 (5 R, 6 man) |
| C4b | `patches/C4b.patch` (`594c828d...`) | 4 R + 2 man | +111/-30 | `C4b.verified.patch` (`cf783fe3...`), +93/-28 | `evidence/C4b/polish.diff`: `R/forestsearch_main.R` +24/-14, `man/forestsearch.Rd` +18/-12 |
| C4c | `patches/C4c.patch` (`9e331f17...`) | 5 R + 5 man | +183/-48 | `C4c.verified.patch` (`adf70151...`), +171/-40 | `evidence/C4c/polish.diff`: 2 R +15/-13, 2 Rd +15/-13 |

- 405 + 111 + 183 = 699 and 18 + 30 + 48 = 96, the commit's counts.
- **Commit = parent + the three final patches.** I exported `c3323a59` (R, man), applied `C4a.patch`, `C4b.patch`, `C4c.patch` with `git apply` (all clean), and compared with `git archive ba7c4ea2`: `diff -r` empty for `R/` and `man/`.
- **Order condition met.** C4a's sentences on `?sg_tables`, `?forestsearch_KfoldOut`, `?cv_summary_tables` need A1 and A2 (C4a NOTES open item 8; C4a VERIFY "Condition"). A1 `f4c7dd2d` and A2 `a2cd36db` precede `ba7c4ea2`. At head the code has the rule: `R/summary_utility_functions.R:856-858` (GRF branch, `fs$adverse_outcome`, fallback `outcome_type %in% c("binary", "count")`), `:1042-1044` (main branch, `args_fs$adverse_outcome`, same fallback); `R/forestsearch_cross_validation.R:1517-1519` (`res$cv_args$adverse_outcome`, same fallback). C4a VERIFY minor 3 (C4a's `man/sg_tables.Rd` hunk does not apply after C5) is moot: C4 is committed before C5 (`5073fd67`).
- **Workflow results:** `workflow1_result.json`: C4a, C4b, C4c each `status` READY, `final_verdict` CONFIRMED, `repaired` false, no blocking problem in round 1. `polish_result.json`: each polish `status` DONE, polish check `verdict` CONFIRMED, `problems` [].

### C4.1 Gates

| Gate | Evidence | Result |
|---|---|---|
| Comment-only (task Rules), each changed `R/` file | `main/C4/comment_gate.txt`: 18 files "comment-only", "COMMENT GATE: PASS". Re-run by me on git: `git diff -U0 c3323a59 ba7c4ea2 -- <file> \| grep -E '^[+-]' \| grep -vE '^(\+\+\+\|---)' \| grep -vE "^[+-][[:space:]]*#"` | Empty for 18 of 18. 0 blank added lines in `R/`. 0 non-ASCII added lines. 0 changed lines with `@examples`, `\dontrun` or `\donttest`. |
| Code identity (extra) | My run: `parse(keep.source = FALSE)`, deparse, `identical()`, `c3323a59` vs `ba7c4ea2` | TRUE for 18 of 18 files (1 to 28 expressions each), 0.12 s |
| `devtools::document()` | `main/C4/document.txt`: `real 0m11.765s`; `document_noop.txt`: `real 0m11.641s` (time output only). Implementers: C4a wrote exactly 10 pages (13.04 s), C4b 2 (13.00 s), C4c 5 (12.89 s); after polish 6 (12.70 s), 1 (12.91 s), 2 (13.03 s); each re-run wrote nothing. Verifiers: `R/`-only copies regenerate `man/` byte-identical to the patches (C4a 12.95 s, polish 11.71 s; C4b 12.56 s; C4c 12.95 s, polish 12.35 s). | My run on `git archive` exports: `document()` writes no file at `ba7c4ea2` (13.18 s) or at `fb80dd8f` (14.46 s); `git status --porcelain` empty. |
| `tools::checkRd()` on every rewritten page (17) | `main/C4/checkRd.txt`: 17 pages "OK (0 problems)" | My run, start / `ba7c4ea2` / head: with the declared encoding (`encoding = "UTF-8", def_enc = TRUE`; `DESCRIPTION` `Encoding: UTF-8`) 0 problems on all 17 at all three. Plain call: `forestsearch.Rd` 14, `make_effect_estimator.Rd` 2, `plot_sg_glm_outcomes.Rd` 1, `plot_subgroup_results_forestplot.Rd` 1, all "Non-ASCII contents without declared encoding" on the "León" lines, identical to start after line numbers are stripped; `generate_bootstrap_with_noise.Rd` 1 at start (`×`) and 0 at `ba7c4ea2` and head. Wall 0.32 to 0.36 s. |
| Scope | `git show --name-only ba7c4ea2` | 18 `R/` files (roxygen `#'`, plus four `#` comment blocks: Site 14's three and Site 4's screening comment), 17 `man/` pages from `document()`. |

### C4.2 Sites

Old and new text are the roxygen / comment source as committed. The `man/` page carries the same text.

#### Site 1 — display orientation (six pages)

| Page | R source, start → head | Rd, start → head |
|---|---|---|
| `?sg_tables` | `R/summary_utility_functions.R`: new `@details`, inserted after start `:735` → head `:736-755` | `man/sg_tables.Rd` after `:82` → `:83-102` |
| `?forestsearch_KfoldOut` | `R/forestsearch_cross_validation.R`: `@return` start `:1324-1325`, extended after `:1325` → head `:1326-1339` | `man/forestsearch_KfoldOut.Rd` after `:29` → `:30-43` |
| `?cv_summary_tables` | `R/cv_summary_tables.R`: new `@details` after start `:34` → head `:35-50` | `man/cv_summary_tables.Rd` after `:52` → `:53-68` |
| `?plot_sg_glm_outcomes` | `R/plot_sg_glm_outcomes.R`: new description paragraph after start `:25` → head `:26-40` | `man/plot_sg_glm_outcomes.Rd` after `:126` → `:127-141` |
| `?plot_subgroup_results_forestplot` | `R/plot_subgroup_results_forestplot.R`: new subsection after start `:132` → head `:133-156` | `man/plot_subgroup_results_forestplot.Rd` after `:197` → `:198-222` |
| `?mr_estimates_table` | `R/mr_estimates_table.R`: new `@details` paragraph after start `:66` → head `:67-84` | `man/mr_estimates_table.Rd` after `:77` → `:78-95` |

**Old text.** None of the six pages said anything on orientation. The nearest start text: `?sg_tables` `@return` "List with gt tables for estimates, subgroups, and optionally search info." (`:736-737`); `?forestsearch_KfoldOut` `@return` "If `outall=TRUE`, a list with summary tables and metrics." (`:1324-1325`).

**New text, `?sg_tables`** (head `R/summary_utility_functions.R:736-755`):
```
#' @details
#' \strong{Orientation of the GLM estimates table.}  For a non-survival fit
#' with the default \code{est.scale = "hr"} (and for a
#' \code{grf_glm_result}), the estimates table shows the per-arm
#' \code{Rate(C)}, \code{Rate(T)} (binary) or \code{Mean(C)}, \code{Mean(T)}
#' (continuous, count), their raw difference \code{Diff} (treated minus
#' control; in percentage points for binary outcomes) and the effect column.
#' The effect column follows the fit's orientation: its estimator is built
#' with the fit's stored \code{adverse_outcome} or, when none is stored,
#' with \code{\link{forestsearch}}'s rule (\code{TRUE} for binary and count
#' outcomes, \code{FALSE} otherwise), so it is on the scale of the fit's own
#' estimates.  The arm columns and \code{Diff} do not follow it: they are
#' the raw outcome.  So for an unadjusted fit where a lower outcome is
#' adverse (\code{adverse_outcome = FALSE}, the default for continuous
#' outcomes) the effect and \code{Diff} have opposite signs: an MD is
#' \code{-Diff}, an RD is \code{-Diff / 100}, and the log of an OR or RR has
#' the sign opposite to \code{Diff}.  Count outcomes and the rate measures
#' (IRR, IRD) are not reoriented (see \code{\link{make_effect_estimator}}),
#' so for them \code{adverse_outcome} does not change the effect column.
```

**New text, `?forestsearch_KfoldOut`** (appended to `@return`, head `R/forestsearch_cross_validation.R:1326-1339`):
```
#'   For GLM outcomes the tables (\code{itt_tab}, \code{SG_tab_original},
#'   \code{SG_tab_Kfold}, \code{tab_all}) are built by
#'   \code{\link{SG_tab_estimates_glm}}.  Their effect column follows the
#'   fit's orientation (the fit's stored \code{adverse_outcome} or, when
#'   none is stored, \code{\link{forestsearch}}'s rule: \code{TRUE} for
#'   binary and count outcomes, \code{FALSE} otherwise), while, for a fit
#'   with the default \code{est.scale = "hr"}, the arm columns
#'   (\code{Rate(C)} or \code{Mean(C)}, \code{Rate(T)} or \code{Mean(T)})
#'   and their raw difference \code{Diff} (treated minus control) do not.
#'   So for an unadjusted fit where a lower outcome is
#'   adverse (\code{adverse_outcome = FALSE}, the default for continuous
#'   outcomes) the effect and \code{Diff} have opposite signs.  Count
#'   outcomes and the rate measures (IRR, IRD) are not reoriented.  See
#'   \code{\link{sg_tables}} for the details.
```

**New text, `?cv_summary_tables`** (head `R/cv_summary_tables.R:35-50`):
```
#' @details
#' The estimate tables are those of \code{kfold_out}, formatted without
#' recomputation.  For GLM outcomes \code{\link{forestsearch_KfoldOut}}
#' builds their effect column with the fit's orientation (the fit's stored
#' \code{adverse_outcome} or, when none is stored,
#' \code{\link{forestsearch}}'s rule: \code{TRUE} for binary and count
#' outcomes, \code{FALSE} otherwise), while, for a fit with the default
#' \code{est.scale = "hr"}, the arm columns (\code{Rate(C)} or
#' \code{Mean(C)}, \code{Rate(T)} or \code{Mean(T)}) and their raw
#' difference \code{Diff} (treated minus control) are the raw outcome.  So
#' for an unadjusted fit where a lower outcome is adverse
#' (\code{adverse_outcome = FALSE}, the default for continuous outcomes)
#' the effect and \code{Diff} have opposite signs.  Count outcomes and the
#' rate measures (IRR, IRD) are not reoriented.  See
#' \code{\link{sg_tables}} for the details.
```

**New text, `?plot_sg_glm_outcomes`** (head `R/plot_sg_glm_outcomes.R:26-40`):
```
#' The effect annotation follows the fit's orientation: its estimator is
#' built with the fit's stored \code{adverse_outcome}
#' (\code{fs.est$args_call_all$adverse_outcome}) or, when none is stored,
#' with \code{\link{forestsearch}}'s rule (\code{TRUE} for binary and count
#' outcomes, \code{FALSE} otherwise), so it is on the scale of the fit's
#' own estimates, sign included.  The plotted per-arm rates or means (and
#' the returned \code{data}) are the raw outcome and do not follow it.  So
#' for an unadjusted fit where a lower outcome is adverse
#' (\code{adverse_outcome = FALSE}, the default for continuous outcomes)
#' the annotated effect and the treated-minus-control difference of the
#' plotted values have opposite signs (for OR and RR, the log of the
#' effect): an MD annotation is the control mean minus the treated mean.
#' Count outcomes and the rate measures (IRR, IRD) are not reoriented (see
#' \code{\link{make_effect_estimator}}).
```

**New text, `?plot_subgroup_results_forestplot`** (head `R/plot_subgroup_results_forestplot.R:133-156`):
```
#' ## Orientation of the GLM rows
#'
#' For a non-survival outcome every row estimated from the data (ITT,
#' reference and post-hoc subgroups) is a GLM effect.  What follows is for
#' the default \code{est.scale = "hr"}.  With a fit
#' (\code{fs_results$fs.est}) its estimator is built with the fit's stored
#' \code{adverse_outcome} or, when none is stored, with
#' \code{\link{forestsearch}}'s rule (\code{TRUE} for binary and count
#' outcomes, \code{FALSE} otherwise), so each row is on the scale of the
#' fit's own estimates, sign included; without a fit the estimator's
#' default (\code{adverse_outcome = TRUE}) is used.  The rows show the arm
#' sizes, not the arm means; the arm means and their raw difference
#' (\code{Diff} in \code{\link{sg_tables}} and
#' \code{\link{mr_estimates_table}}) do not follow the orientation.  So for
#' an unadjusted fit where a lower outcome is adverse
#' (\code{adverse_outcome = FALSE}, the default for continuous outcomes) a
#' row's effect and the raw treated-minus-control difference have opposite
#' signs (for OR and RR, the log of the effect): an MD row is the control
#' mean minus the treated mean.  Count outcomes and the rate measures (IRR,
#' IRD) are not reoriented (see \code{\link{make_effect_estimator}}).
#' Bias-corrected rows are taken from \code{fs_results$fs_bc}, whose
#' estimator \code{\link{forestsearch_bootstrap_dofuture}} builds with the
#' same rule.
```

**New text, `?mr_estimates_table`** (head `R/mr_estimates_table.R:67-84`):
```
#' For GLM outcomes the effect columns follow the fit's orientation: the
#' naive and MR values are those of `fs$mr_inference`, on the scale of the
#' fit's own estimates; the descriptive table's estimator (which gives the
#' complement's naive cell when MR did not cover the complement) is built
#' with the fit's stored `adverse_outcome` or, when none is stored, with
#' [forestsearch()]'s rule (`TRUE` for binary and count outcomes, `FALSE`
#' otherwise); and a template's effect column comes from
#' [forestsearch_bootstrap_dofuture()], which builds its estimator with the
#' same rule.  With the default `est.scale = "hr"`, the descriptive columns
#' `Rate(C)` or `Mean(C)`, `Rate(T)` or `Mean(T)`, and their raw difference
#' `Diff` (treated minus control) are the raw outcome and do not follow it.
#' So for an unadjusted fit where a
#' lower outcome is adverse (`adverse_outcome = FALSE`, the default for
#' continuous outcomes) the effect and `Diff` have opposite signs: an MD is
#' `-Diff`, an RD is `-Diff / 100`, and the log of an OR or RR has the sign
#' opposite to `Diff`.  Count outcomes and the rate measures (IRR, IRD) are
#' not reoriented (see [make_effect_estimator()]).
```

**Code each sentence rests on** (start; head in parentheses):

| Sentence | Code |
|---|---|
| Arm columns are the raw outcome | `R/summary_utility_functions.R:203-204` `rate0 <- mean(Y[Treat == 0])`, `rate1 <- mean(Y[Treat == 1])`; labels `:344-346` (head same) |
| `Diff` is treated minus control, x 100 for binary | `:251` `sprintf(fmt_signed_pct, 100 * (rate1 - rate0))`; `:255` raw otherwise (head same) |
| Effect column is the estimator's | `:219` `estimator_fn(df_sub)` (head same) |
| Stored `adverse_outcome`, else `forestsearch()`'s rule | `R/forestsearch_main.R:1840-1842` (head `:1877-1879`) `adverse_outcome <- (outcome_type %in% c("binary", "count"))`, stored by `mget()` at `:1924` (head `:1961`); A1/A2 at head (C4.0); F12 displays `R/plot_sg_glm_outcomes.R:169-172` (head `:184-187`), `R/plot_subgroup_results_forestplot.R:316-322` (head `:340-346`), `R/mr_estimates_table.R:204-216` (head `:222-234`); bootstrap `R/bootstrap_dofuture_main.R:306-320` |
| "on the scale of the fit's own estimates" | `forestsearch()` builds its estimator with the same value, `R/forestsearch_main.R:2270` (head `:2307`) |
| MD is `-Diff` | `R/glm_effect_estimators.R:956-958` (head `:1045-1047`) fits on `-Y`; fast closure `:1093-1094` (head `:1182-1183`) |
| RD is `-Diff / 100`; log OR, log RR opposite sign | `:390-392` (head `:408-410`) fits on `1L - Y` |
| "unadjusted" | an adjusted effect is a model coefficient, not a function of the raw arm means (`.build_adjusted_formula()`, `:1247`, head `:1336`) |
| Count, IRR, IRD not reoriented | `:162-171` (head `:180-189`), `:198-213` (head `:216-231`) call `.make_poisson_rate_estimator()` without `adverse_outcome`; its formals `:737-746` (head `:826-835`) |
| `est.scale = "hr"` scope | `prepare_subgroup_data()` `R/summary_utility_functions.R:51-58` swaps `treat` on `"1/hr"` (called `:366`) while the closure keeps the original (`:219`); `sg_tables()` passes `args_fs$est.scale` (start `:1002`, `:1018`; head `:1069`, `:1085`), GRF branch `"hr"` (start `:849`, `:865`; head `:883`, `:899`); `forestsearch_KfoldOut()` `res$cv_args$est.scale` (start `:1376`, head `:1390`); forest plot swaps `df_fs$treat` (start `:727-731`, head `:751-755`) |
| CV tables built by `SG_tab_estimates_glm()`, formatted without recomputation | `R/forestsearch_cross_validation.R:1519`, `:1583`, `:1602`, `:1624` (head `:1544`, `:1609`, `:1628`, `:1650`); `R/cv_summary_tables.R:84-87`, `:97-113` (head `:100-103`, `:113-129`) |
| Forest plot rows show arm sizes; no fit → `TRUE`; bias-corrected rows from `fs_bc` | `row[[E.name]] <- ntreat` `:477` (head `:501`); `:316` (head `:340`); `create_bc_row()` `:599-617` (head `:623-641`) |
| MR table cells | naive/MR from `fs$mr_inference` `:290-292` (head `:308-310`); complement naive cell `:286-293` (head `:304-311`); template `FSsg_tab` `:89` (head `:107`), built `R/bootstrap_dofuture_main.R:714-727` |

**Measured** (`evidence/C4a/confirm_*.txt`, `verify/probe_*.txt`, `polish/probe_estscale_*.txt`):
- Continuous fit, default orientation (stored `FALSE`), `max_sg_est` 1.094112: plot annotation `MD = 1.094`, plotted means Control -0.330 / Treatment -1.424; forest plot H row `est` 1.0941; `mr_estimates_table()` `Diff` -1.0941, `MD` 1.0941. `sg_tables()` at start: `MD` -1.094 (same sign as `Diff`); start + A1 + A2: `Diff` -1.094, `MD` 1.094; `forestsearch_KfoldOut()` FA_0 `Diff` -1.0941, `MD` 1.0941.
- Binary, estimator oriented `FALSE`: Questionable `Diff` +6.1982 (points), `RD` -0.0620, `OR` 0.7688, `RR` 0.8500.
- Verifier, own generator: `sg_tables()` `Diff` -2.118692, `MD` 2.118692; binary RD 0.245274 vs `Diff` -24.527363%; `cv_summary_tables(use_gt = FALSE)` `identical()` to `tab_all` and `SG_tab_original`; forest plot without a fit uses `TRUE` (ITT -0.2816907 = hand -0.281691).
- Fit with `est.scale = "1/hr"` (stack start + A1 + A2): Questionable `Mean(C)` 0.246039 (the treated mean), `Diff` +2.118692, `MD` +2.118692 (same sign). Polish verifier: with `"hr"` `MD = -Diff` (-2.060310); with `"1/hr"` same sign.

**Differs from the task's wording** ("State that the effect column follows the fit's orientation, while the arm means and their raw difference do not, so for a fit where a lower outcome is adverse the two have opposite signs."):
- "for an unadjusted fit": with `adjust_covariates` the effect is not a function of the raw arm means, so the sign relation holds only unadjusted.
- Signs per measure added: MD is `-Diff`; RD is `-Diff / 100` (`Diff` is in points for binary); for OR and RR the sign is that of the log.
- "Count outcomes and the rate measures (IRR, IRD) are not reoriented": the Poisson estimator takes no `adverse_outcome`.
- "with the default `est.scale = "hr"`" (polish, five pages): on a `"1/hr"` fit the arm columns and `Diff` are swapped and `Diff` has the effect's sign. The polish list said to qualify only pages whose function has an `est.scale` argument (2 of the 5: `mr_estimates_table()`, `plot_subgroup_results_forestplot()`). The polisher also qualified `sg_tables()`, `forestsearch_KfoldOut()` and `cv_summary_tables()`, worded for the fit, because they read the fit's own `est.scale` and the probes show the sentence fails there. C4a NOTES gives the revert path. Not on `?plot_sg_glm_outcomes` (the file reads no `est.scale`).
- `?plot_sg_glm_outcomes` has no `Diff`: the sentence compares the annotation with the treated-minus-control difference of the plotted values. The forest plot shows arm sizes, not means; the page adds the no-fit default `TRUE` and the source of the bias-corrected rows.
- "the CV summary tables": written on `?forestsearch_KfoldOut` (builder) and `?cv_summary_tables` (formatter). `cv_metrics_tables()`, `cv_summary_text()`, `cv_compare_results()` show no effect column and are unchanged.
- Three of the six pages describe behaviour after A1/A2 (met: C4.0).

**Polish** (C4a, `evidence/C4a/polish.diff`): the five `est.scale = "hr"` qualifiers above (old → new in C4a NOTES "Polish" table). Polish check: CONFIRMED, no problem.

**Verifier problems:** C4a round 1 minor 2 (est.scale) → fixed by polish; the `"1/hr"` code behaviour → open item. Minor 3 (C5 integration) → moot by commit order. Minor 4 (needs A1/A2) → met by commit order.

#### Site 2 — `make_effect_estimator()` `@param adverse_outcome`

- `R/glm_effect_estimators.R` start `:54-60` (changed `:56-60`) → head `:54-78` (changed `:56-78`); `man/make_effect_estimator.Rd` start `:45-49` → head `:45-67`.

Old (start `:54-60`):
```
#' @param adverse_outcome Logical. If \code{TRUE} (default), higher
#'   outcome values indicate harm and the effect is computed on the
#'   raw outcome.  If \code{FALSE}, the outcome is beneficial (e.g.,
#'   improvement) and the effect is computed on \code{1 - Y} internally
#'   so that effect > 0 (or OR > 1) consistently indicates treatment
#'   harm.  Only applies to binary outcomes; ignored for survival,
#'   continuous, and count.
```
New (head `:54-78`):
```
#' @param adverse_outcome Logical. If \code{TRUE} (default), higher
#'   outcome values indicate harm and the effect is computed on the
#'   raw outcome.  If \code{FALSE}, higher values indicate benefit
#'   (e.g., improvement) and, where it applies, the effect is computed
#'   on a reversed outcome so that effect > 0 (or a ratio > 1) indicates
#'   treatment harm.  By outcome type and effect measure:
#'   \itemize{
#'     \item \code{"binary"} with \code{"OR"}, \code{"RR"} or \code{"RD"}:
#'       with \code{FALSE} the outcome is replaced by \code{1 - Y} before
#'       fitting.  The log-OR and the RD change sign; the log-RR becomes
#'       the log risk ratio of the complementary event,
#'       \code{log((1 - p1) / (1 - p0))}, which is not \code{-log(p1 / p0)}
#'       (\code{p1}, \code{p0}: the treated and control risks).
#'     \item \code{"continuous"} (\code{"MD"}): with \code{FALSE} the
#'       outcome is replaced by \code{-Y} before fitting, so the MD
#'       changes sign.
#'     \item \code{"count"} (\code{"IRR"}, \code{"IRD"}), and
#'       \code{"binary"} with \code{"IRR"} or \code{"IRD"}: ignored; the
#'       Poisson rate estimator always uses the raw outcome.
#'     \item \code{"survival"}: ignored; the Cox estimator always uses the
#'       raw time and event.
#'   }
#'   \code{\link{forestsearch}} resolves its own default (\code{TRUE} for
#'   binary and count outcomes, \code{FALSE} otherwise) and passes it here;
#'   this function's default is \code{TRUE} for every outcome type.
```

| Sentence | Code |
|---|---|
| binary OR/RR/RD: `1 - Y` | `:173-183` passes it (`:181`; head `:191-201`, `:199`) → `:390-392` (head `:408-410`) before the RD/OR/RR fits `:414-434` (head `:432-452`) |
| continuous: `-Y` | `:187-196` passes it (`:193`; head `:205-214`, `:211`) → `:956-958` (head `:1045-1047`); fast closure `:1093-1094` (head `:1182-1183`) |
| count, binary IRR/IRD ignored | `:162-171`, `:198-213` (head `:180-189`, `:216-231`); Poisson formals `:737-746` (head `:826-835`) |
| survival ignored | `:141-150` (head `:159-168`) calls `.make_cox_estimator()` without it |
| function default `TRUE`; `forestsearch()` resolves its own | `:107` (head `:125`); `R/forestsearch_main.R:1840-1842` (head `:1877-1879`), passed `:2270` (head `:2307`) |

- Measured, start tree (`confirm_claims.txt`): `TRUE` → `FALSE`: OR -0.0763 → +0.0763; RD -0.0188 → +0.0188 (hand p1 - p0 = -0.018767); RR -0.0333 → +0.0429 (hand log(p1/p0) = -0.033350, log((1-p1)/(1-p0)) = 0.042930); MD 0.1508 → -0.1508; binary IRR, IRD, count IRR, IRD, survival `identical()` TRUE. Verifier (p1 0.55, p0 0.555): RR -0.009050 → +0.011173 = hand complementary ratio.
- **Differs from the task's wording** ("the linear-model estimators negate the outcome. Say what the code does for each outcome type"): no, that part is as the task says. Added from the code: binary IRR/IRD also ignore it (the old "Only applies to binary outcomes" was wrong there), and the RR on `1 - Y` is the complementary-event ratio, not the reciprocal.
- Polish: none. Verifier problems: none.

#### Site 3 — `forestsearch()` "Resolved defaults": a supplied threshold

- `R/forestsearch_main.R` start `:1313-1315` (changed `:1314-1315`) → head `:1335-1352` (changed `:1336-1352`); `man/forestsearch.Rd` start `:1213-1214` → head `:1228-1244`.

Old (start `:1313-1315`):
```
#' Supplying \code{effect.threshold} or \code{consistency.threshold} (or the
#' legacy names) suppresses the corresponding remap and the supplied value is
#' used on the scale shown above.
```
New (head `:1335-1352`):
```
#' Supplying \code{effect.threshold} or \code{consistency.threshold} (or the
#' legacy names) suppresses the corresponding default remap, and the supplied
#' value is used on the scale shown above, with one exception: for \code{RD}
#' and \code{IRD} a supplied value greater than 1 is taken to be on a ratio
#' scale and is replaced by the default in the table: \code{effect.threshold}
#' by 0.05 for \code{RD} or 0.01 for \code{IRD}, with two warnings (one from
#' the input check, one from the remap), and \code{consistency.threshold} by
#' 0, with one warning.  The remapped \code{effect.threshold}
#' is the value the fit records (in \code{threshold_config} and
#' \code{args_call_all}, which a bootstrap or cross-validation replay
#' re-uses), and it is also the DINA harm floor wherever DINA selects a cut
#' (\code{subgroup_method = "dina"}, or \code{use_dina = TRUE} screening with
#' \code{selected_only = TRUE}): on these two routes \code{RD} and \code{IRD}
#' run with DINA only on the \code{gaussian} family, whose floor is the
#' resolved \code{effect.threshold} (see \code{dina_args}).  Screening with
#' \code{selected_only = FALSE} has neither the restriction nor this floor:
#' it contributes the frontier candidate set on any family, for \code{RD} and
#' \code{IRD} too.
```

| Clause | Code |
|---|---|
| a supplied value suppresses the default remap | `user_set_*` `:1848-1849` (head `:1885-1886`); default remap only when not user-set `:2323-2329` (head `:2360-2366`) |
| greater than 1, replaced by 0.05 / 0.01 | `:2334` `user_set_threshold && effect_threshold > 1.0` (head `:2371`); `:2348-2349` `switch(effect_measure, RD = 0.05, IRD = 0.01)` (head `:2376`, after C2) |
| two warnings | input check `:301-309` (head `:301-316`, text changed by C2) and the remap warning `:2335-2347` (head within `:2371-2395`) |
| consistency to 0, one warning | `:2355` (head `:2400`), `:2366` (head `:2411`); the input check has no identity-scale check on `hr.consistency` |
| recorded in `threshold_config`, `args_call_all`; replay | `:2479-2483` (head `:2524-2528`); `:2573-2580` (head `:2623-2630`); `R/bootstrap_dofuture_main.R:268`, `R/forestsearch_cross_validation.R:99` |
| DINA floor where DINA selects a cut | `subgroup_method = "dina"`: `:2823` (head `:2873`) → `R/forestsearch_helpers.R:1528-1529` `m_diff <- if (identical(da$fit$family, "gaussian")) effect.threshold`; `use_dina` + `selected_only`: `:3272` (head `:3322`), `:3287-3289` (head `:3341-3343`) |
| gaussian only on those two routes | `.dina_assert_ratio_estimand()` `R/forestsearch_helpers.R:1446-1463`, called at `:1527` and `R/forestsearch_main.R:3286` (head `:3340`) |
| `selected_only = FALSE`: no restriction, no floor | returns `fr$cut_expr` `:3314` (head `:3368`); `dina_frontier()` with `da$frontier` `:3263-3267` (head `:3313-3317`), `m_diff` = `dina_args$m_diff` or `NULL` (`R/forestsearch_helpers.R:1165`) |

- Probes: implementer `repro_site1.txt` (8.72 s fixed / 9.20 s start, identical output): RD `effect.threshold = 1.5` → `threshold_config$screening` 0.05, `args_call_all$effect.threshold` 0.05, `m_diff` 0.05 on both DINA routes; IRD 0.01; `consistency.threshold = 1.5` → 0. Verifier `probe_site1` (10.33 s start / 10.51 s fixed): RD 1.5, legacy `hr.threshold = 2`, IRD 3 remapped with 2 warnings each; 1 kept; -1.5 kept with 0 warnings; replay screening 0.05 with 0 warnings. Polish `probe_p2_start.txt` (9.29 s): `selected_only = FALSE` gives 6 frontier cuts, `m_diff = NULL`, 0 `dina_subgroup()` calls, 0 assert calls. Polish verifier `probe_s1` (5.7 s): exactly 2 warnings (effect) and 1 (consistency); 0.08 / 0.02 kept with 0 warnings.
- **Differs from the task's wording** ("An out-of-range risk- or rate-difference value is remapped with a warning, and the remapped value is also the DINA floor"):
  - "greater than 1", not "out-of-range": the code tests `> 1.0` only; a value below -1 is kept with no warning.
  - "two warnings" for `effect.threshold`, one for `consistency.threshold`.
  - `consistency.threshold` remap to 0 added (the task does not name it).
  - The DINA floor is limited to the two routes where DINA selects a cut; `selected_only = FALSE` has no floor.
- **Polish** (C4b P1, P2): "read as a ratio-scale value and is remapped, with a warning" → "taken to be on a ratio scale and is replaced ... with two warnings ..."; "RD and IRD run with DINA only on the gaussian family" → "on these two routes ..." plus the `selected_only = FALSE` sentence.
- **Verifier problems:** C4b minor 1 (two warnings; "read as") → fixed (P1). Minor 3 (gaussian-only holds only on those routes) → fixed (P2). Minor 2 (`args_call_all` also keeps the supplied `hr.threshold` beside the remapped `effect.threshold`; the replay still uses 0.05) → not in the polish list, unchanged; open item.

#### Site 4 — `forestsearch()` `@param dina_args`: "the same cut `subgroup_method = "dina"` would select"

- `R/forestsearch_main.R` start `:440-445` (changed `:443-445`) → head `:447-459` (changed `:450-459`); `man/forestsearch.Rd` start `:219-221` → head `:219-228`.
- Extra site (not in the task's table): the `#` comment in the `use_dina` screening block, start `:3273-3275` (changed `:3274-3275`) → head `:3323-3329` (changed `:3324-3329`).

Old (start `:440-445`):
```
#'   When \code{selected_only = TRUE} (the default), \code{use_dina} screening
#'   contributes the cut(s) chosen by \code{\link{dina_subgroup}} -- using
#'   this call's \code{sg_focus} / \code{selection_rule} /
#'   \code{effect_neighborhood} / \code{n.min} / \code{hr.threshold}, i.e. the
#'   same cut \code{subgroup_method = "dina"} would select. Set \code{FALSE} to
#'   contribute the full frontier candidate set instead.
```
New (head `:447-459`):
```
#'   When \code{selected_only = TRUE} (the default), \code{use_dina} screening
#'   contributes the cut(s) chosen by \code{\link{dina_subgroup}} -- using
#'   this call's \code{sg_focus} / \code{selection_rule} /
#'   \code{effect_neighborhood} / \code{n.min} / \code{hr.threshold}, ranked by
#'   \code{sg_focus} with DINA's native tau-hat as the effect statistic (under
#'   \code{"maxSG"} / \code{"minSG"} size ranks first and tau-hat breaks
#'   ties).  That is the cut \code{subgroup_method = "dina"}
#'   would select only under \code{dina_select_statistic = "dina"}; under the
#'   default \code{dina_select_statistic = "effect"},
#'   \code{subgroup_method = "dina"} re-ranks the qualifying family on the
#'   inferential effect and can select a different cut, while this screening
#'   step does not re-rank.  Set \code{FALSE} to contribute the full frontier
#'   candidate set instead.
```
Old comment (start `:3273-3275`):
```
        # Selected-cut screening: contribute the SINGLE cut dina_subgroup()
        # selects, using forestsearch's own selection criteria (so this
        # matches what subgroup_method = "dina" would choose).  m_diff is the
```
New comment (head `:3323-3329`):
```
        # Selected-cut screening: contribute the SINGLE cut dina_subgroup()
        # selects, using forestsearch's own selection criteria with DINA's
        # native tau-hat as the effect statistic.  That is the cut
        # subgroup_method = "dina" would choose only under
        # dina_select_statistic = "dina"; under the default "effect" that
        # route re-ranks the qualifying family on the inferential effect
        # (.dina_reselect_on_effect()), and this step does not.  m_diff is the
```

| Clause | Code |
|---|---|
| default is `"effect"` | `:1688` `dina_select_statistic = c("effect", "dina")` (head `:1725`); `match.arg` `:1819` (head `:1856`); injected into `dina_args` `:1830` (head `:1867`) |
| `subgroup_method = "dina"` re-ranks under `"effect"` | `R/forestsearch_helpers.R:1611-1624`; `.dina_reselect_on_effect(` called only at `:1614` |
| screening does not re-rank | screening block `R/forestsearch_main.R:3272-3312` (head `:3322-3366`) emits the `dina_subgroup()` winner, rounded by `signif(, digits)` at `:3309` (head `:3363`) |
| ranked by `sg_focus` on tau-hat; size first under maxSG/minSG | `R/dina_subgroup.R:535-560` (head `:537-562`), `eff` = tau-hat `:519` (head `:521`) |

- Probes: implementer `repro_site2.txt` (22.87 s): continuous MD, 8 seeds, `"effect"` and `"dina"` select different cuts in 7 of 8; the screening cut equals the `"dina"` cut. Verifier `probe_site2_start.txt` (57.31 s): binary OR, 10 runs; screening = `"dina"` cut in 10 of 10; `"effect"` differs in 5 of 10, all with `sg_focus = "hr"`. Polish verifier `probe_s45` (18.0 s): `.dina_reselect_on_effect` called once on the DINA route and 0 times in screening; under maxSG / minSG ties the winner has the largest tau-hat.
- **Differs from the task's wording** ("holds for the default selection statistic only"): yes, inverted. The code's default is `"effect"`, under which `subgroup_method = "dina"` re-ranks; the claim holds only for the non-default `"dina"`. Also added: the ranking statistic (polish P5). Not added: the `digits` rounding of the screening threshold, which the same paragraph already states (start `:471-472`, head `:485-486`).
- **Polish:** P4 (the `#` comment, verifier minor 5) and P5 ("ranked on DINA's native tau-hat" → "ranked by `sg_focus` with DINA's native tau-hat as the effect statistic (...)", verifier minor 7).
- **Verifier problems:** C4b minor 5 → comment fixed (P4); its second part, `R/forestsearch_helpers.R:1607-1610` calling `"dina"` "the default", is outside C4b's files and unchanged (open item). Minor 6 (proposed NEWS bullet stronger than "can select a different cut") → not addressed; the committed NEWS keeps the stronger form (C4.4). Minor 7 → fixed (P5).

#### Site 5 — `forestsearch()` `@param consistency_method` and `.glm_resample_supported()`

- (a) `R/forestsearch_main.R` start `:1074-1078` (changed `:1076-1078`) → head `:1088-1100` (changed `:1090-1100`); `man/forestsearch.Rd` start `:813-815` → head `:820-830`.
- (b) `R/subgroup_consistency_helpers.R` start `:336-338` (changed `:337-338`) → head `:336-342` (changed `:337-342`); `man/dot-glm_resample_supported.Rd` (internal page) start `:20-21` → head `:20-25`.

Old (a) (start `:1076-1078`):
```
#'   coefficient (OR, RR, RD, MD, IRR); configurations it cannot represent that
#'   way (IRD, propensity-adjusted effects, or a non-convergent fit on a given
#'   candidate) fall back to \code{"split"} automatically. When
```
New (a) (head `:1090-1100`):
```
#'   coefficient (OR, RR, RD, MD, IRR); an IRD effect, or a candidate whose
#'   closed form is not available (for example a non-convergent fit), falls
#'   back to \code{"split"} automatically.  The route is decided by the
#'   effect measure alone: propensity adjustment
#'   (\code{ps_adjust_method = "iptw"} or \code{"dr_gcomp"}, with
#'   \code{ps_method} other than \code{"none"}) does not send a supported
#'   measure to \code{"split"}.  Under propensity adjustment the
#'   closed form fits the treatment effect without the propensity weights,
#'   and without covariates (\code{adjust_covariates} cannot be combined with
#'   \code{ps_method} other than \code{"none"}), so the consistency rate is
#'   that of the unadjusted effect. When
```
Old (b) (start `:336-338`):
```
#' IRR (Poisson with offset). It does not hold for IRD (a delta-method rate
#' difference) or propensity-adjusted (IPTW / G-computation) effects, for which
#' the caller falls back to literal splitting.
```
New (b) (head `:336-342`):
```
#' IRR (Poisson with offset). It does not hold for IRD (a delta-method rate
#' difference), for which the caller falls back to literal splitting.
#' The decision reads `spec$effect_measure` only.  The spec carries no
#' propensity setting, so a propensity-adjusted (IPTW / G-computation)
#' analysis with a supported measure also takes the resampling path, whose
#' model is fitted without the propensity weights (see
#' `.consistency_glm_pieces()`).
```

| Clause | Code |
|---|---|
| decided by effect measure alone | `.glm_resample_supported()` `R/subgroup_consistency_helpers.R:343-346` (head `:347-350`); spec fields `R/forestsearch_main.R:3675-3684` (head `:3729-3738`), no propensity field; routing `:1502-1504`, `:1863-1866` (head `:1510-1512`, `:1873-1876`) |
| IRD, or no closed form → split | `R/consistency_resample.R:249` (IRD → `NULL`); NA closed form `R/subgroup_consistency_helpers.R:1519-1527` (head `:1527-1535`), two-stage `:1886-1891` (head `:1896-1901`); failed fit `R/consistency_resample.R:302-305` |
| propensity adjustment defined | `R/forestsearch_main.R:2683` `ps_adjust_resolved <- if (ps_method == "none") "none" else ps_adjust_method` (head `:2733`) |
| closed form without weights | `R/consistency_resample.R:267-298` (`rhs <- c(treat.name, adj)`, no `weights`) |
| without covariates under propensity adjustment | `R/forestsearch_main.R:2693-2697` (head `:2743-2747`) stops `adjust_covariates` with `ps_method != "none"` |

- Probes: implementer `repro_site3.txt` (3.23 s / 3.55 s start, identical): iptw OR, `consistency_resample()` 10 calls, splits 0. Verifier `probe_site3` (9.35 s / 9.43 s): OR+iptw 8/0, OR+dr_gcomp 12/0, RD+iptw 10/0, MD+iptw 1/0; closed-form rates `identical()` to the `ps_method = "none"` run (8 of 8, 11 of 11); `probe_site3b` IRD 0 closed form, 480 split calls. Polish verifier `probe_s3`: `beta_hat` equals an unweighted `glm(y ~ treat)` (max abs difference 0).
- **Differs from the task's wording** ("the code decides by effect measure only"): no. Added: what follows under propensity adjustment (unadjusted closed form).
- **Polish** P3 (verifier minor 4): "The closed form fits ... without covariates" → "Under propensity adjustment the closed form fits ... without covariates", since without propensity adjustment the closed form includes `adjust_covariates` (`R/consistency_resample.R:267-268`).

#### Site 6 — `CV_sgs()` `@return`

- `R/forestsearch_cross_validation.R` start `:1813-1814` (changed `:1814`) → head `:1839-1865` (changed `:1840-1865`); `man/CV_sgs.Rd` start `:20` → head `:20-45`.

Old (start `:1813-1814`):
```
#' @return List with indicators for any match, exact match, one match, and
#'   covariate-specific matches.
```
New (head `:1839-1865`):
```
#' @return List with indicators for any match, exact match, one match, and
#'   covariate-specific matches, one value per fold (\code{sg1[i]},
#'   \code{sg2[i]} are the factor labels of fold \code{i}'s subgroup,
#'   \code{NA} for an empty slot):
#'   \describe{
#'     \item{any_found}{1 if the fold found a subgroup (either slot not
#'       \code{NA}), else 0.}
#'     \item{exact_match}{1 if the fold's subgroup is the full-data
#'       subgroup \code{sg_analysis}, else 0.  For a two-factor
#'       \code{sg_analysis}, both factors appear in the fold's two slots,
#'       in either order.  For a one-factor \code{sg_analysis}, the fold's
#'       subgroup is that factor alone: one slot holds it and the other is
#'       \code{NA}.  A fold subgroup that combines it with a second factor
#'       is not exact; it counts in \code{one_match}.}
#'     \item{one_match}{1 if at least one full-data factor appears in
#'       either slot, else 0.}
#'     \item{cov1_any, cov2_any}{1 if the fold's subgroup uses the
#'       covariate of the first (second) full-data factor, with any cut
#'       (see \code{\link{find_covariate_any_match}}), else 0.}
#'     \item{cov1_exact, cov2_exact}{1 if the first (second) full-data
#'       factor appears in either slot, else 0.}
#'   }
#'   \code{cov2_any} and \code{cov2_exact} are \code{NA} unless
#'   \code{sg_analysis} has two factors.  A comparison with \code{NA}
#'   counts as 0.  When \code{sg_analysis} has no element, or more than
#'   two, \code{exact_match}, \code{one_match}, \code{cov1_any} and
#'   \code{cov1_exact} are 0 for every fold.
```

| Item | Code (start; head) |
|---|---|
| `any_found` | `:1826` (`:1877`) |
| two-factor exact, either order | `:1834` (`:1885`) |
| one-factor exact: that factor alone | `:1860` `(sg1 == sg1a & is.na(sg2)) \| (sg2 == sg1a & is.na(sg1))` (`:1911`) (part 4's F4) |
| `one_match`, `cov1_exact` | `:1838`, `:1842`, `:1845`; one-factor `:1856`, `:1863` (`:1889`, `:1893`, `:1896`, `:1907`, `:1914`) |
| `cov*_any` | `:1849-1850` (`:1900-1901`) |
| `cov2_*` NA unless two factors | `:1865`, `:1869` (`:1916`, `:1920`) |
| no element or more than two | `else` branch `:1871-1879` (`:1922-1930`) |

- Probes: implementer one-factor `exact_match` (A,NA) 1, (A,B) 0, (B,A) 0, (NA,A) 1, (A,A) 0, (NA,NA) 0; `one_match` 1 1 1 1 1 0. Verifier 9 folds, two-factor either order. Polish `probe_cvsgs_start.txt` (2.44 s) and polish verifier (2.33 s): `character(0)`, `NULL`, three and four factors give the four indicators 0 and `cov2_*` `NA`.
- **Differs from the task's wording:** no for "Exact". Added: every indicator, the NA rule, and the empty / more-than-two branch.
- **Polish** (verifier minor 1): "every indicator except `any_found` is 0" → names the four indicators that are 0 (`cov2_*` are `NA` there).

#### Site 7 — `summarize_bootstrap_subgroups()` bins

- `R/summarize_bootstrap_subgroups.R` start `:22-23` → head `:22-35`; `man/summarize_bootstrap_subgroups.Rd` start `:24-25` → head `:24-37`.

Old (start `:22-23`):
```
#'     \item{consistency_dist}{Data.table of Pcons distribution by bins}
#'     \item{size_dist}{Data.table of subgroup size distribution}
```
New (head `:22-35`):
```
#'     \item{consistency_dist}{Data.table of Pcons distribution by bins
#'       (columns \code{Consistency Range}, \code{Count}, \code{Percent} =
#'       \code{100 * Count / n_found}).  The bins are closed on the left:
#'       a label \code{"a-b"} holds \code{a <= Pcons < b} (for example
#'       \code{"0.85-0.95"} holds \code{0.85 <= Pcons < 0.95}),
#'       \code{"<0.5"} holds \code{Pcons < 0.5}, and the top bin
#'       \code{">=0.95"} holds \code{0.95 <= Pcons <= 1}.  A value on a
#'       break counts in the bin that starts there (0.95 is in
#'       \code{">=0.95"}).}
#'     \item{size_dist}{Data.table of subgroup size distribution
#'       (columns \code{Size Range}, \code{Count}, \code{Percent}), with
#'       bins closed on the left in the same way: \code{"<50"} holds
#'       \code{N_sg < 50}, \code{"50-99"} holds \code{50 <= N_sg < 100},
#'       and so on to \code{">=300"}, which holds \code{N_sg >= 300}.}
```
- Code: `:417-423` (head `:429-435`) `breaks = c(0, 0.5, 0.7, 0.8, 0.85, 0.95, 1.0)`, `include.lowest = TRUE`, `right = FALSE` (part 4's F1); size `:442-451` (head `:454-463`); `Percent` `:427-431`, `:455-459` (head `:439-443`, `:467-471`).
- Probes: implementer Pcons 0, 0.5, 0.85, 0.95, 1.0, 0.949 → `<0.5` 1, `0.5-0.7` 1, `0.85-0.95` 2, `>=0.95` 2. Verifier `probe_bins_fn.txt` (2.44 s), through the function: 0.95 and 1.0 → `>=0.95`; 0.4999 → `<0.5`; N_sg 50, 99 → `50-99`; 300 → `>=300`; `Percent` = 100 x Count / 11 (`n_found`, including the NA-Pcons row).
- Differs from the task's wording: no (labels, columns and size bins added). Polish: none. Verifier problems: none.

#### Site 8 — `test_hte_from_forestsearch()` `@return`

- `R/test_hte_crump.R` start `:386` → head `:386-417`; `man/test_hte_from_forestsearch.Rd` start `:33` → head `:33-64`.

Old (start `:386`): `#' @return List with hte_test, covariates_used, source.`

New (head `:386-417`):
```
#' @return List with \code{hte_test} (the \code{\link{test_hte}} result),
#'   \code{covariates_used} (the \code{df} columns entered in the test),
#'   \code{n_covariates} (their number), \code{n_candidates}
#'   (\code{length(fs_result$find.grps$max_sg_est)}, \code{NA} when that
#'   is absent) and \code{source}; or \code{NULL} with a warning when no
#'   covariate is found.
#'
#'   The covariates come from \code{fs_result$confounders.evaluated} when
#'   it is non-empty, else from \code{fs_result$confounders.candidate},
#'   taken as column names.  Each \code{confounders.evaluated} label is
#'   mapped to columns of \code{df}: a label that is exactly a column name
#'   (a categorical confounder such as \code{"meno"}, including a
#'   non-syntactic name such as \code{"ER-pos"}) is that column; any other
#'   label (a cut such as \code{"age <= 53"}) is parsed, and its variables
#'   that are columns of \code{df} are used.  Each column is used once, in
#'   the order first met.  Two warnings: a label that resolves to no column
#'   of \code{df} warns that it "cannot be resolved to a column of df" and
#'   is skipped; and when no covariate column is left, the function warns
#'   "No covariate found in df for the ForestSearch covariates" and returns
#'   \code{NULL}.  (When neither field is present it warns "No covariates
#'   found." and returns \code{NULL}.)  Factor and character columns enter
#'   the test as their level codes minus 1, other columns as numbers.
#'
#'   \code{source} names the field used:
#'   \code{"confounders.evaluated (GRF/LASSO selected)"} or
#'   \code{"confounders.candidate (all candidates)"}.  The first string is
#'   fixed.  On the default (consistency) route
#'   \code{confounders.evaluated} holds the fit's candidate factor labels,
#'   which are LASSO-screened or include GRF cuts only when the fit used
#'   \code{use_lasso} or \code{use_grf} (both \code{FALSE} by default); on
#'   the \code{subgroup_method = "dina"} and \code{"grf"} routes it is the
#'   fit's \code{confounders.name}.
```

| Clause | Code (start; head) |
|---|---|
| elements; `n_candidates` | list `:440-449` (`:471-480`); `:443-445` (`:474-476`) |
| `confounders.evaluated` first, else candidates as column names | `:397-398` (`:428-429`); `:417-420` (`:448-451`) |
| exact column name; else parsed | `:406` (`:437`); `:407-409` (`:438-440`) |
| once each, first-met order | `:413`, `:415` (`:444`, `:446`) |
| warning 1 | `:410-412` (`:441-443`) |
| warning 2 | `:432-436` (`:463-467`); candidate names not in `df` skipped silently `:427` (`:458`) |
| "No covariates found." | `:422` (`:453`) |
| factor coding | `:428-430` (`:459-461`) |
| `source` string fixed | `:416` (`:447`) and `:420` (`:451`) |
| what `confounders.evaluated` holds | `R/forestsearch_main.R:4009` (head `:4063`), `:3420` (head `:3474`); `R/get_fsdata.R:752` (head `:754`); `use_lasso`/`use_grf` defaults `:1680-1681` (head `:1717-1718`); DINA/GRF routes `:2840`, `:3029` (head `:2890`, `:3079`) |

- Probes: implementer GBSG default fit: `use_lasso`, `use_grf` FALSE; `source` "confounders.evaluated (GRF/LASSO selected)"; `covariates_used` age size nodes pgr er meno; `n_candidates` 1; each warning text reproduced. Verifier `probe_hte.txt` (3.41 s): `ER-pos` used as a column; factor enters as 0, 1, 2; existing test `hte-crump-confounders` 35 passed (5.77 s).
- **Differs from the task's wording** ("Describe the mapping of cut labels to data columns and the two warnings"): no. Added: `n_covariates`, `n_candidates`, the third message "No covariates found.", factor coding, and what `source` means. The `source` string is wrong for a default fit and the DINA/GRF routes; changing it is code, outside the comment-only gate; the page now says what the field holds (open item).
- Polish: none. Verifier problems: none.

#### Site 9 — `calibrate_k_treat()`

- `R/calibrate_k_treat.R` start `:9-15` → head `:9-18`; `@seealso` start `:90` → head `:93-94`. `man/calibrate_k_treat.Rd` start `:59-65` → head `:59-68`; `\seealso` line added at head `:118`.

Old (start `:9-15`, `:90`):
```
#' This is the natural sibling of \code{\link{calibrate_k_inter}} (which
#' targets the harm-subgroup HR under \code{model = "alt"}) and
#' \code{\link{calibrate_cens_adjust}} (which targets the censoring rate).
#' It is typically used with \code{model = "null"} to fix the overall HR
#' for a uniform-effect DGM, but also works with \code{model = "alt"}
#' when the marginal HR is the calibration target rather than the
#' subgroup-specific HR.
...
#' @seealso \code{\link{calibrate_k_inter}},
```
New (head `:9-18`, `:93-94`):
```
#' This is the natural sibling of \code{\link{find_k_inter_for_target_hr}}
#' (which targets the harm-subgroup HR of a \code{generate_aft_dgm_flex()}
#' DGM under \code{model = "alt"}) and \code{\link{calibrate_cens_adjust}}
#' (which targets the censoring rate). \code{\link{calibrate_k_inter}}
#' targets the harm-subgroup HR of the GBSG mechanism of
#' \code{\link{setup_gbsg_dgm}} instead.
#' \code{calibrate_k_treat()} is typically used with \code{model = "null"}
#' to fix the overall HR for a uniform-effect DGM, but also works with
#' \code{model = "alt"} when the marginal HR is the calibration target
#' rather than the subgroup-specific HR.
...
#' @seealso \code{\link{find_k_inter_for_target_hr}},
#'   \code{\link{calibrate_k_inter}},
```
- Code: `calibrate_k_treat()` builds with `generate_aft_dgm_flex()` `:133-136` (head `:137-140`); `find_k_inter_for_target_hr()` with `generate_aft_dgm_flex(..., model = "alt")`, returns `$hazard_ratios$harm_subgroup` (`R/find_k_inter_main.R:105-121`, head `:109-125`); `calibrate_k_inter()` with `.create_gbsg_dgm_()` (`R/sim_aft_gbsg.R:673-679`), requires `model = "alt"` (`:668-669`); `setup_gbsg_dgm()` calls the same builder (`R/setup_gbsg_dgm.R:154-165`, head `:191-202`); `calibrate_cens_adjust()` default target `"rate"` (`R/simulate_from_dgm.R:790`).
- **Differs from the task's wording** ("It calls `calibrate_k_inter()` the general-generator counterpart"): the page said "the natural sibling", not "the general-generator counterpart". The substance is as the task says; the text also says what `calibrate_k_inter()` calibrates, and `@seealso` gains the link.
- **Polish** (verifier minor 1): "It is typically used ..." → "`calibrate_k_treat()` is typically used ...", since the preceding sentence's subject is `calibrate_k_inter()`, which requires `model = "alt"`.
- **Verifier problems:** minor 1 → fixed. Minor 7 (`@param k_treat_range` "(uniroot will error otherwise)", start `:26-27`, head `:29-30`; `extendInt = "yes"` widens automatically) → outside the site list, unchanged; open item.

#### Site 10 — `find_k_inter_for_target_hr()`

- `R/find_k_inter_main.R`: `@param continuous_vars` start `:10-11` (changed `:11`) → head `:10-13`; `@return` start `:35`, `:42-43` → head `:37`, `:44-45`; `@details` start `:51-52` → head `:53-56`. `man/find_k_inter_for_target_hr.Rd` start `:32`, `:69`, `:76-77`, `:90-91` → head `:32-34`, `:71`, `:78-79`, `:92-95`.

Old (start):
```
#' @param continuous_vars Character vector of continuous variable names to be
#'   standardized and included as covariates.
...
#' @return A list of class "k_inter_result" containing:
...
#'   \item{convergence}{Integer number of iterations to convergence}
#'   \item{method}{Character string "root-finding" indicating method used}
...
#' precision (within the specified tolerance). If the root-finding fails, the
#' function evaluates the boundaries and provides diagnostic information.
```
New (head):
```
#' @param continuous_vars Character vector of continuous variable names to be
#'   included as covariates on their original scale: the function calls
#'   \code{\link{generate_aft_dgm_flex}} at its default
#'   \code{standardize = FALSE} and has no argument to change it.
...
#' @return A plain list (no class attribute) containing:
...
#'   \item{convergence}{Integer number of \code{uniroot} iterations, including
#'     any automatic extensions of the search range}
...
#' precision (within the specified tolerance). The search range is widened
#' automatically until it brackets the target; if no root can be bracketed, or
#' the achieved HR is not within \code{tol_rel} percent of the target, the call
#' stops with an error.
```
- Code: no `standardize` passed (`:106-120`, `:139-153`; head `:110-124`, `:143-157`), no `standardize` or `...` formal (`:88-102`, head `:92-106`); generator default `standardize = FALSE` (`R/generate_aft_dgm_main.R:240`); raw column copied (`R/generate_aft_dgm_helpers.R:495-499`). Result `list(...)` with no class and no `method` (`:166-173`, head `:170-177`); `convergence = sol$iter` from `stats::uniroot(..., extendInt = "yes")` (`R/calibrate_helpers.R:61`); bracket failure → `stop()` `:62-71`; tolerance gate → `stop()` `:77-84`.
- Probes: implementer `run_s2_s3_find_k_inter.txt`: `class` "list"; names `k_inter achieved_hr_harm target_hr_harm error dgm convergence`; `convergence` 11; `z_age`, `z_er`, `z_pgr` ranges [21, 80], [0, 1144], [0, 2380]; `run_s2b_failure_path.txt`: `tol_rel = 1e-9` stops. Verifier `run_p1` (8.43 s): 19 iterations from a start range `c(-0.2, -0.1)`.
- **Differs from the task's wording** (`method` element; "standardized" where the default does not standardize): as the task says, plus: the false class `"k_inter_result"` removed; "no argument to change it" (the function has none, so not only the default); `@details` failure-path sentence replaced (beyond the task's list; separable 3-line edit, C4c NOTES departure 3).
- Polish: none. Verifier problems: minor 6 ("typically converges within 5-10 iterations", unchanged at head `:52`; 11 and 19 observed) → unchanged; open item.

#### Site 11 — `compute_dgm_cde()` `@param dgm`

- `R/oc_analyses.R` start `:109-114` (changed `:110-111`, `:113`) → head `:109-118` (changed `:110-115`, `:117`); `man/compute_dgm_cde.Rd` start `:11-12`, `:15` → head `:11-16`, `:19`.

Old (start `:109-113`):
```
#' @param dgm A DGM object (e.g., from \code{generate_aft_dgm_flex} or
#'   \code{setup_gbsg_dgm}). Must contain \code{df_super_rand}
#'   with columns \code{theta_0} and \code{theta_1}.
#' @param harm_col Character. Name of the subgroup indicator column in
#'   \code{dgm$df_super_rand}. If \code{NULL} (default), auto-detected
```
New (head `:109-117`):
```
#' @param dgm A DGM object (e.g., from \code{generate_aft_dgm_flex} or
#'   \code{setup_gbsg_dgm}). The super-population is read from
#'   \code{dgm$df_super_rand} when present (as from \code{setup_gbsg_dgm}),
#'   otherwise from \code{dgm$df_super} (as from
#'   \code{generate_aft_dgm_flex}); it must have columns \code{theta_0} and
#'   \code{theta_1}. If neither data frame is present, or those columns are
#'   missing, the function warns and returns \code{dgm} unchanged.
#' @param harm_col Character. Name of the subgroup indicator column in
#'   that super-population data frame. If \code{NULL} (default), auto-detected
```
- Code: `df <- dgm$df_super_rand %||% dgm$df_super` `:134` (head `:138`), `%||%` `:27`; warn and return `:135-138`, `:139-142` (head `:139-142`, `:143-146`); indicator `:146-147`, `:160` (head `:150-151`, `:164`); `setup_gbsg_dgm()` objects carry both frames (`R/sim_aft_gbsg.R:529`, `R/setup_gbsg_dgm.R:209`, head `:246`); `generate_aft_dgm_flex()` objects `df_super` only (`R/generate_aft_dgm_helpers.R:1320-1330`).
- Probes: implementer: `aft_dgm_flex` object `CDE_harm` 2.122295 = hand from `df_super`; `compute_dgm_cde(list(a = 1))` warns and returns its input `identical()`. Verifier `run_p2` (2.53 s): with both frames and `df_super$theta_1` perturbed, `CDE_harm` 1.504105 = hand from `df_super_rand` (perturbed `df_super` would give 4.088582).
- Differs from the task's wording: no. Added: `@param harm_col` (same wording fault) and the warn-and-return behaviour. Polish: none. Verifier problems: none.

#### Site 12 — `generate_bootstrap_with_noise()`

- `R/synthetic_data_perturbation.R`: `@param continuous_vars` start `:427-429` (changed `:429`) → head `:427-431`; `@return` start `:449-450` → head `:451-456`; Details bullets start `:461`, `:463` → head `:467`, `:469-471`; `@note` start `:517-518` → head `:525-537`. `man/generate_bootstrap_with_noise.Rd` start `:25`, `:50-51`, `:74`, `:76`, `:100-101` → head `:25-27`, `:52-57`, `:80`, `:82-84`, `:108-120`.

Old (start):
```
#'   used, including numeric columns with few distinct values.
...
#' @return A data frame with the same structure as the input data, containing
#'   bootstrap sampled observations with added noise.
...
#'     \item Adds Gaussian noise with standard deviation = original SD × noise_level
#'     \item Constrains values to remain within original variable bounds
#'     \item Preserves integer type for variables that appear to be integers
...
#'   \item Missing values (NA) are handled appropriately in calculations but
#'     are not imputed
```
New (head):
```
#'   used, including numeric columns with few distinct values. A column named
#'   here must be numeric: a factor or character column stops the call with an
#'   error, and a logical column is returned as a numeric 0/1 column.
...
#' @return A data frame of \code{n} rows, with the columns of \code{data} in
#'   the same order and row names \code{1:n}, containing bootstrap sampled
#'   observations with added noise. Column types are kept, except that every
#'   column in \code{continuous_vars} is returned as double (an integer column
#'   comes back double, holding whole numbers) and \code{id_var}, if present,
#'   is replaced by the integer sequence \code{1:n}.
...
#'     \item Adds Gaussian noise with standard deviation = original SD times noise_level
#'     \item Constrains values to remain within original variable bounds
#'     \item Rounds to whole numbers the variables that appear to be integers
#'       (see below); the column is returned as double, not integer
#'     \item Leaves missing values missing
...
#'   \item Missing values (NA) in a column treated only as continuous (not
#'     also in \code{cat_vars}) stay missing. In a
#'     categorical column a missing value is selected for change with the same
#'     probability, \code{noise_level / 2}, as any other value; if it is not
#'     selected it stays missing. A selected missing value in an ordered
#'     factor, or a column named \code{grade}, with more than two levels stops
#'     the call with the error "replacement has length zero". Otherwise, in a
#'     column with two levels (for a non-factor column, two distinct
#'     non-missing values in \code{data}) it is replaced by the first level
#'     (for a non-factor column, the smallest value); in a column with more
#'     than two levels, by a level drawn at random from all its levels; and in
#'     a column with fewer than two levels (for a non-factor column, fewer than
#'     two distinct non-missing values) it stays missing.
```

| Clause | Code (start; head) |
|---|---|
| factor `continuous_vars` stops; character stops; logical → 0/1 | `sd()` `:605` (`:624`) "Calling var(x) on a factor x is defunct."; `+ noise` `:609` (`:628`) "non-numeric argument to binary operator"; logical coerced at `:609`, bounded `:610-612`, rounded `:615-617` (`:629-631`, `:634-636`) |
| `n` rows, columns, row names `1:n`, `id_var` → `1:n` | `:598-599` (`:617-618`); `:684` (`:703`); `:687-689` (`:706-708`) |
| continuous → double, whole numbers for integers | `+ rnorm()` `:606`, `:609` (`:625`, `:628`); `round()` `:615-617` |
| NA in continuous-only column stays NA | `NA + noise`, `pmax`/`pmin` without `na.rm`, `round(NA)` `:609-617` (`:628-636`) |
| categorical NA selected with prob `noise_level / 2` | `which(runif(n) < flip_prob)` `:634-635` (`:653-654`) |
| ordered / `grade`, more than two levels → error | `which(levels_vec == current)` is `integer(0)` for NA, then `weights[j] <- ...` `:653-660` (`:672-679`) |
| two levels → first level (non-factor: smallest value) | levels `:625-631` (`:644-650`); `setdiff(levels_vec, NA)[1]` `:640-643` (`:659-662`) |
| more than two → random level | `sample(other_levels, 1)` `:665` (`:684`) |
| fewer than two → stays NA | only `n_levels == 2` (`:638`) and `n_levels > 2` (`:645`) change values (`:657`, `:664`) |
| column in both sets loses selected NAs | auto `cat_vars` drop `continuous_vars` `:584` (`:603`); categorical step `:622-681` (`:641-700`) |
| SD times `noise_level` | `:605` (`:624`) |

- Probes: implementer `run_s4_bootstrap_noise.txt` (2.42 s): integer column → double, whole; `survival::gbsg` defaults: 10 of 11 integer columns come back double, `pid` integer `1:686`; two-level factor 3 selected NAs → `"lo"`; 0/1 numeric 4 → 0; three-level factor 10 → x 5, y 3, z 2; ordered three-level 8 selected → "replacement has length zero"; numeric `grade` 1/2/3/NA → same error. Polish `run_polish_probe.txt` (2.65 s): 1-level and 0-level columns keep selected NAs (7 of 7; 24 of 24); column in both sets 12 NAs in → 7 out. Verifier `run_p3` (2.54 s) and polish verifier (2.56 s) reproduce these.
- **Differs from the task's wording** (integer type; "same structure"; NA not imputed: a selected NA in a categorical column is replaced, an ordered factor with more than two levels and an NA stops; `continuous_vars` must be numeric):
  - A logical `continuous_vars` column is accepted (returned 0/1), so "must be numeric" is stated with that exception.
  - A column named `grade` with more than two levels also stops; the stop happens only when an NA is selected (`noise_level = 0` runs).
  - Two-level → first level; more than two (unordered) → random level; fewer than two → stays NA (polish).
  - NA kept only in a column treated only as continuous (polish).
  - `×` replaced by "times" (non-ASCII, start `:461`). With the declared encoding the start page already had 0 `checkRd` problems; the plain call reported 1 (my run).
- **Polish** (verifier minors 2, 3, 5): fewer-than-two-levels clause added; "in continuous columns" → "in a column treated only as continuous (not also in `cat_vars`)"; "holding whole numbers when every non-missing input value is whole" → "holding whole numbers".
- **Verifier problems:** minor 4 (a one-column data frame drops to a vector, no noise) → behaviour, unchanged; open item. Polish verifier's non-problem: an `id_var` named in `continuous_vars` loses its NAs to `1:n` (`@return` states the `1:n` reset).

#### Site 13 — `setup_gbsg_dgm()` Details

- `R/setup_gbsg_dgm.R` start `:74` (`#' @details Noise columns ...`) → head `:74-111` (new `:74-110`; the noise paragraph now starts at head `:111`). `man/setup_gbsg_dgm.Rd`: inserted after start `:105` → head `:106-141`.
- Source recovered: `git show d93e6621:R/sim_aft_gbsg.R`, `create_gbsg_dgm()` `@details`, lines 134-169 (headings "Subgroup Definition", "Model Specification", "Interaction Effect (k_inter)", "Alignment with generate_aft_dgm_flex"; checked by me).

Old (start `:74`): `#' @details Noise columns requested via \code{k_random_noise} are inert` (the noise paragraph `:74-78` is otherwise unchanged).

New (head `:74-111`):
```
#' @details
#' **Subgroup.** The harm subgroup H is \code{z1 = 1} and \code{z3 = 1}:
#' \code{z1} is 1 when the estrogen receptor value \code{er} is at or below
#' its \code{z1_quantile} sample quantile in \code{survival::gbsg} (the 25th
#' percentile by default), and \code{z3} is 1 for premenopausal patients
#' (\code{meno == 0}). Its indicator is \code{flag.harm} in
#' \code{df_super_rand} and \code{flag_harm} in \code{df_super}. Under
#' \code{model = "null"} the indicator is 0 for every patient.
#'
#' **Model.** A Weibull accelerated failure time model is fitted with
#' \code{survival::survreg()} to \code{survival::gbsg}, with time
#' \code{rfstime / 30.4375} (months), event \code{status} and treatment
#' \code{treat = hormon}. Its covariates are \code{treat}, \code{z1},
#' \code{z2}, \code{z3}, \code{z4} and \code{z5}, where \code{z2}, \code{z4}
#' and \code{z5} are the quartile groups of \code{age}, \code{pgr} and
#' \code{nodes}, entered as scores 1 to 4; under \code{model = "alt"} the
#' interaction \code{zh = treat * z1 * z3} is added. \code{k_treat} and
#' \code{k_z3} multiply the fitted \code{treat} and \code{z3} coefficients.
#' With AFT coefficients \code{gamma} and scale \code{sigma}, the
#' hazard-scale coefficients are \code{-gamma / sigma}
#' (\code{model_params$b_hr}).
#'
#' **Interaction (`k_inter`).** Under \code{model = "alt"}, \code{k_inter}
#' multiplies the fitted \code{zh} coefficient,
#' `gamma[zh] <- k_inter * gamma[zh]` (the default \code{k_inter = 1} keeps
#' the fitted value). Each patient's log hazard ratio
#' (\code{loghr_po = theta_1 - theta_0} in \code{df_super_rand}) is then
#' `-(gamma[treat] + gamma[zh]) / sigma` in H and `-gamma[treat] / sigma` in
#' its complement, so \code{k_inter = 0} gives every patient the same hazard
#' ratio. The exponentials of these two values are the stored
#' \code{AHR_H_true} and \code{AHR_Hc_true}. The Cox-based
#' \code{hr_H_true} and \code{hr_Hc_true} are fitted to stacked potential
#' outcomes of the \code{n_super} patients drawn with replacement from the
#' GBSG rows, and differ from them; at \code{k_inter = 0} they still differ
#' slightly from each other. Under \code{model = "null"}, \code{k_inter} is
#' not used.
#'
#' Noise columns requested via \code{k_random_noise} are inert
```

| Sentence | Code (`R/sim_aft_gbsg.R`, same at head, unless noted) |
|---|---|
| H is `z1 = 1` and `z3 = 1`; `z1` from the `er` quantile; `z3` premenopausal | `:139` `flag.harm <- ifelse(z1 == 1 & z3 == 1, 1L, 0L)`; `:128-129`, `:131`; default 0.25 `:90`, `R/setup_gbsg_dgm.R:141` (head `:178`) |
| `flag.harm` / `flag_harm` | `:529`; rename `R/setup_gbsg_dgm.R:207-209` (head `:244-246`) |
| null: indicator 0 | `:160` |
| Weibull AFT on gbsg; time in months; `treat = hormon` | `survreg` `:224-228` on `dfa` `:119`; `:122`, `DAYS_PER_MONTH <- 30.4375` `:32`; `:124-125` |
| covariates; quartile scores | `covs_main` `:181`; `:130`, `:132`, `:133` with `cut_numeric()` `:51-58` |
| `zh` under alt | `:136`, `:183-188` |
| `k_treat`, `k_z3` | `:265-266` |
| `-gamma / sigma` = `b_hr` | `b0 <- -gamma / sigma` `:290`; `b_hr = b0` `:551` |
| `k_inter` multiplies `gamma[zh]` | `:268-270` `if (model == "alt") gamma["zh"] <- k_inter * gamma["zh"]` |
| patient log HR in H and Hc | `theta_0`, `theta_1`, `loghr_po` `:316-318`, carried `:337-339` |
| exp = `AHR_H_true`, `AHR_Hc_true` | `:408-409` `exp(mean(df_samp$loghr_po[...]))` |
| Cox `hr_*_true` from stacked potential outcomes of `n_super` resampled patients | `:327-331`, `:347-394` |
| null: `k_inter` unused | `:268` guard; `covs_true <- covs_main` `:190` |

- Probes: implementer `run_s5_gbsg_model.txt` (3.25 s): `er` threshold 8 = `quantile(gbsg$er, 0.25)`; `loghr_po` H 0.4081982, Hc -0.5365232 (one value each); `AHR_H_true` 1.504105 vs Cox `hr_H_true` 1.387292; `AHR_Hc_true` 0.584778 vs `hr_Hc_true` 0.661188; at `k_inter = 0` AHR 0.584778 both, Cox 0.655173 vs 0.661188. Verifier `run_p4` (3.02 s): an independent `survreg` refit (covariates by `findInterval`) gives the same `gamma` (tolerance 1e-8) and `sigma`.
- **Differs from the task's wording** ("Add a Details paragraph carrying the model description ... subgroup definition, model, `k_inter`"):
  - Three paragraphs with bold lead-ins, not one paragraph and not `##` headings (so the existing noise paragraph does not fall under a subsection).
  - The old page's `HR(H) = exp(-gamma[treat]/sigma - gamma[zh]/sigma)`, `HR(Hc) = exp(-gamma[treat]/sigma)` and "When k_inter = 0, HR(H) = HR(Hc)" hold for each patient's HR and `AHR_*_true`, not for the Cox `hr_*_true`; the new text says which is which.
  - Added from the builder: null-model indicator and `k_inter`; quartile scores; months; `k_treat`/`k_z3`; `b_hr`; `n_super` resampling. The old "Alignment" section is carried only as far as `loghr_po` and AHR.
- Polish: none. Verifier problems: none.

#### Site 14 — three `#` comments

(a) `R/dina_subgroup.R` start `:403-404` (changed `:404`) → head `:403-406` (changed `:404-406`). Task pointer "near 402".

Old: `  # gaussian, log-OR for binomial, log-IRR for poisson, log-HR for cox.`
New:
```
  # gaussian, log-OR for binomial, the log mean-count ratio for poisson
  # (dina() takes no exposure offset, so this is log-IRR only when mean
  # follow-up does not depend on arm given the covariates), log-HR for cox.
```
Code: `dina()` has no offset argument (`R/dina_df.R:84-104`); `dina_fit()` formals none (`R/dina.R:157-167`); log link (`R/dina.R:394-396`); same wording as `@param m_diff` (`R/dina_subgroup.R:138-143`).

(b) `R/forestsearch_main.R` start `:2494-2500` (changed `:2495-2500`) → head `:2539-2550` (changed `:2540-2550`). Task pointer "near 2457-2458".

Old:
```
    # Resolve dmin.grf for GLM outcomes.
    # The survival default for dmin.grf is on the RMST scale (event count
    # threshold over a window) and would filter out all candidates for
    # binary/continuous outcomes where DR score differences are on the
    # probability/mean scale (0--1).  When the user has not supplied a
    # value, default to 0.0 (any positive harm qualifies); when the user
    # supplied a value explicitly, respect it.
```
New:
```
    # Resolve dmin.grf for GLM outcomes.
    # dmin.grf is a harm threshold on GRF's doubly robust (DR) scale, not a
    # number of events: a subgroup's DR harm effect (mean control DR score
    # minus mean treated DR score) must reach it (>= dmin.grf).  For survival
    # that effect is a difference in restricted mean survival time up to the
    # GRF horizon; for GLM outcomes it is a difference on the scale of the
    # outcome given to GRF (a probability for binary, outcome units for
    # continuous, log(Y + 0.5) or raw counts for count, per
    # grf_count_transform), so a value chosen on the RMST scale does not
    # carry over.  When the user has not supplied a value, use 0.0 (the
    # formal default: any non-negative DR harm effect qualifies); a supplied
    # value is used as given.
```
Code: DR harm effect `R/grf_helpers.R:132` (`control mean - treated mean`), `R/grf_subg_harm_glm.R:600-610`; `>=` `R/grf_helpers.R:207`, `:221`, `R/grf_subg_harm_glm.R:917`; RMST horizon `R/grf_main.R:241-245`, `R/grf_helpers.R:74-86`; count transform `R/grf_subg_harm_glm.R:408-433`; formal default `R/forestsearch_main.R:1694` (head `:1731`); `missing()` reset `:2501-2503` (head `:2551-2553`). The old claims "event count threshold", "would filter out all candidates" (default is 0.0 on every path), "any positive harm" (the test is `>=`) and "(0--1)" are false and are dropped.

(c) `R/get_fsdata.R` start `:327` → head `:327-329`.

Old: `    # Reduce dimension via LASSO (Cox for survival, GLM for binary/continuous)`
New:
```
    # Reduce dimension via LASSO: Cox for survival; a GLM otherwise --
    # logistic for binary (Poisson when offset.name is given), linear for
    # continuous, Poisson for count (log offset when offset.name is given).
```
Code: `lasso_selection()` family by `outcome_type` (`R/get_FSdata_helpers.R:286-311`); `R/get_fsdata.R:328-335` (head `:330-337`) passes `outcome_type` and `offset.name`.

- Verifier (C4b VERIFY §3 site 4): all three correct. Differs from the task's wording: only the line pointers (402 → 403-404; 2457-2458 → 2494-2500). Polish: none.

### C4.3 Verifier problems and what happened

| Item | Problem (round 1, minor) | Outcome |
|---|---|---|
| C4a 1 | `?CV_sgs` "every indicator except `any_found` is 0" (`cov2_*` are `NA`) | Fixed by polish (Site 6) |
| C4a 2 | Orientation sentences hold for `est.scale = "hr"` only | Fixed by polish on five pages (Site 1); `"1/hr"` code behaviour → open item |
| C4a 3 | `man/sg_tables.Rd` hunk does not apply after C5 | Moot: C4 committed before C5 |
| C4a 4 | Three pages true only after A1/A2 | Met: A1, A2 committed first; rule present at head |
| C4b 1 | "remapped, with a warning": two warnings; "read as" | Fixed (P1) |
| C4b 2 | `args_call_all` keeps the supplied `hr.threshold` beside the remapped value | Not in the polish list; unchanged; open item |
| C4b 3 | Gaussian-only holds only on the two DINA-selects routes | Fixed (P2) |
| C4b 4 | "without covariates" read as general | Fixed (P3) |
| C4b 5 | Screening `#` comment still says "matches"; `R/forestsearch_helpers.R:1607-1610` calls `"dina"` the default | Comment fixed (P4); helpers comment outside C4b's files, unchanged; open item |
| C4b 6 | Proposed NEWS bullet stronger than "can select a different cut" | Not addressed; committed NEWS keeps it (C4.4) |
| C4b 7 | "ranked on DINA's native tau-hat" under maxSG/minSG | Fixed (P5) |
| C4b (outside) | Count/IRD `use_twostage = TRUE`, `fs.splits = 10`: "argument must be coercible to non-negative integer" | Not investigated; open item |
| C4c 1 | "It is typically used" subject | Fixed (Site 9) |
| C4c 2 | One-level column: selected NA stays NA | Fixed ("fewer than two levels") |
| C4c 3 | Column in both sets loses NAs | Fixed ("only continuous") |
| C4c 4 | One-column data frame drops to a vector | Behaviour; unchanged; open item |
| C4c 5 | Redundant "when" clause | Fixed |
| C4c 6 | "5-10 iterations" | Unchanged; open item |
| C4c 7 | `@param k_treat_range` "uniroot will error otherwise" | Outside the site list; unchanged; open item |

- No blocking problem in any round. No repair round. Each polish was checked by a second verifier: CONFIRMED, `problems` [] (C4a, C4b, C4c).
- The polish verifiers recomputed each polish delta from the trees: 68, 68 and 56 changed lines, equal to `evidence/C4*/polish.diff`.

### C4.4 NEWS check (`fb80dd8f:NEWS.md:536-562`, "Documentation pass (part 5 ...)", under `### Documentation`)

| NEWS line | Claim | Check |
|---|---|---|
| 536 | "no behaviour change" | True for C4: parse-identical code, 18 of 18 files |
| 537-541 | six display pages "say that the effect column follows the fit's orientation while the arm means or rates and `Diff` stay raw, so for a fit where a lower outcome is adverse the two can have opposite signs" | Pages named correctly. Inexact in detail: the pages state it for an **unadjusted** fit and (five pages) for the default `est.scale = "hr"`; NEWS drops both. "can have" keeps it from being false. `?plot_sg_glm_outcomes` shows no `Diff` and the forest plot shows arm sizes, not means. |
| 541-542 | `?make_effect_estimator` says what `adverse_outcome` does for each outcome type | Correct |
| 544-545 | "an out-of-range risk- or rate-difference threshold is remapped with a warning" | **Inexact.** Only values greater than 1 are remapped (start `R/forestsearch_main.R:2334`, `:2355`; head `:2371`, `:2400`); a value below -1 is kept with no warning (C4b verifier probe, -1.5). A supplied `effect.threshold` gets **two** warnings; `consistency.threshold` one. The help says "greater than 1" and "two warnings". |
| 545 | "the remapped value is also the DINA floor" | **Inexact.** Only the remapped `effect.threshold`, and only where DINA selects a cut (`subgroup_method = "dina"`, or `use_dina` with `selected_only = TRUE`); `selected_only = FALSE` has no floor (help, polish P2). |
| 546-547 | "the `use_dina` screening cut matches `subgroup_method = "dina"` only under `dina_select_statistic = "dina"`" | **Stronger than the help and the probes.** Under the default `"effect"` the route "can select a different cut"; the cuts were the same in 5 of 10 verifier runs and 1 of 8 implementer seeds. C4b round-1 minor 6 flagged this for the proposed bullet; C4b's polish NOTES recommended "can select a different cut under the default `"effect"`". |
| 547-548 | "what `consistency_method` does under propensity adjustment" | Correct |
| 549-555 | `CV_sgs()`, bins, `test_hte_from_forestsearch()`, the four corrected pages, `setup_gbsg_dgm()` | Correct |
| 556-562 | C5 and C7 bullets | Not C4. Light check only: C5 `5073fd67` regenerates 16 pages and its subject says "sixteen"; C7 `835fac1b` adds "It is not numerically equivalent to `setup_gbsg_dgm()`" (head `paper_simulations.qmd:395`). |

- Not mentioned in NEWS: the three `#` comments and the internal `.glm_resample_supported` page (not user-facing).

### C4.5 Open items (left by this item)

1. `test_hte_from_forestsearch()` `source` string "confounders.evaluated (GRF/LASSO selected)" (start `R/test_hte_crump.R:416`, head `:447`) is wrong for a default fit and on the DINA/GRF routes. Code string, outside the comment-only gate; `tests/testthat/test-hte-crump-confounders.R:76` asserts it. `n_candidates` is `length(max_sg_est)`, so 1 or `NA`.
2. `est.scale = "1/hr"` GLM fits: `prepare_subgroup_data()` (`R/summary_utility_functions.R:51-58`) swaps the treatment for the arm columns and `Diff` while the estimator keeps the original. `sg_tables()` and the CV tables then label the treated mean `Mean(C)` and show `Diff` with the effect's sign. Help scoped to `"hr"`; code unchanged.
3. `mr_estimates_table(est.scale = "1/hr")` without `include_complement`: H cell inverted (2.7238 → 0.3671), the complement's naive cell not (1.1172). Observation, not analysed.
4. NEWS documentation entry: the inexact claims in C4.4.
5. `args_call_all` keeps the supplied `hr.threshold` (e.g. 1.5) beside the remapped `effect.threshold` (0.05). The replay uses 0.05.
6. `R/forestsearch_helpers.R:1607-1610` (start = head) says "The default "dina" path leaves `sg` untouched"; the default statistic is `"effect"`.
7. `consistency_method = "resample"` under propensity adjustment gives the unadjusted effect's consistency rate (now documented); whether to refuse, route to `"split"` or extend the path is unresolved (`R/consistency_resample.R:216-221`).
8. `CV_sgs()` with a full-data subgroup of three or more factors gives 0 for every indicator but `any_found` (head `R/forestsearch_cross_validation.R:1922-1930`); `cv_metrics_tables()` still prints "Exact match on all factors" (start `R/cv_summary_tables.R:418`, head `:434`).
9. `?forestsearch` `@return` "confounders.evaluated: Confounders after variable selection" (start `R/forestsearch_main.R:1203`, head `:1225`): same fault as item 1.
10. `?plot_sg_glm_outcomes` `@param outcome_type` lists binary and continuous only (start `:37`, head `:52`); count outcomes are drawn too.
11. `grf.subg.harm.glm()` defaults `adverse_outcome = FALSE` for every outcome type (`R/grf_subg_harm_glm.R:304`), unlike `forestsearch()`'s rule.
12. The `dmin.grf` `missing()` reset re-assigns the formal default 0.0 (start `:2501-2503`, head `:2551-2553`): no effect.
13. `compute_dgm_cde()`'s second warning says "df_super_rand lacks theta_0/theta_1" also when `df_super` was read (start `R/oc_analyses.R:140`, head `:144`).
14. `find_k_inter_for_target_hr()` `@details` "typically converges within 5-10 iterations" (head `R/find_k_inter_main.R:52`): 11 (help example) and 19 (narrow start range) observed. The help example prints two `coxph.fit` "Loglik converged before variable 1" warnings.
15. `calibrate_k_treat()` `@param k_treat_range` "(uniroot will error otherwise)" (start `:26-27`, head `:29-30`): `extendInt = "yes"` widens the range.
16. `generate_bootstrap_with_noise()`: a one-column data frame drops to a vector with no noise; an `id_var` also in `continuous_vars` loses its NAs to `1:n`.
17. Count/IRD fits: "argument must be coercible to non-negative integer" in `subgroup.consistency` (caught as "Consistency analysis failed") with `fs.splits = 10`, `use_twostage = TRUE` (C4b verifier, start tree; also every IRD count run of the C4b polish verifier's synthetic data). Not investigated; not caused by C4.
18. Plain `tools::checkRd()` still reports the pre-existing "León" notes (14 + 2 + 1 + 1 on four pages); 0 with the declared encoding.

### C4.6 Evidence problems

1. **C4b NOTES §2, Site 1 table** cites `R/forestsearch_helpers.R:1554-1556` as the gaussian floor and `:1554` as a call of `.dina_assert_ratio_estimand()`. At `b06de4e8` the call is at `:1527` and the floor `m_diff <- if (identical(da$fit$family, "gaussian")) effect.threshold` at `:1528-1529`; `:1554-1556` is the printed "Harm floor: m_diff = ..." line. C4b's polish section (P2, ":1527") and C4b VERIFY §3 (":1527", ":1528") give the right lines.
2. **`main/C4/document.txt` and `document_noop.txt`** hold only `time` output (`real 0m11.765s`, `real 0m11.641s`); they do not say which files `document()` wrote. My re-run: no file written at `ba7c4ea2` or head.
3. **`main/C4/checkRd.txt`** reports "OK (0 problems)" on all 17 pages without naming the call. That matches only the declared-encoding call; the plain call gives 14 notes on `forestsearch.Rd` and 2, 1, 1 on three other pages (all pre-existing).
4. **C4a stack checks used different A1/A2 patch versions.** Round-1 verifier: A1 `6539fa4e`, A2 `46e5e27d`. Polish implementer: A2 `495d45ae`. Polish verifier: A1 `efe668ab`, A2 `62818222`, called "current". `patches/A1.patch` is now `b12a282a` (modified 17:04:13, after the polish check). The dependency is met at head (rule checked directly in C4.0); only the "current" hash statement is stale.
5. **C4b NOTES polish P2** cites `R/forestsearch_main.R:3310-3312` for "the screening returns `fr$cut_expr`" under `selected_only = FALSE`. At `b06de4e8` `:3310-3312` is `} else { character(0) }` (no cut selected); `fr$cut_expr` is at `:3314` (head `:3368`), and the `dina_frontier()` call starts at `:3263`, not `:3262`. The claim itself holds.
6. **C4c NOTES §8 item 7** says the `×` was replaced "so the rewritten page passes `checkRd` with no problems, as the gate requires". That holds for the plain call and for `encoding = "UTF-8"` alone; with `def_enc = TRUE` (the call C4a/C4b NOTES name as R CMD check's) the start page already had 0 problems (my run). Not a contradiction; the change is harmless.

## C5 — examples that cannot run as written (detail)

Commit `5073fd67` (parent `ba7c4ea2`). Line numbers: start `b06de4e8` and head `fb80dd8f`. `R/` and `man/` are identical at `5073fd67`, `835fac1b` and `fb80dd8f` (`git diff --stat 5073fd67 fb80dd8f` lists no `R/` or `man/` file), so head lines are `5073fd67` lines. `$S = /tmp/cran_prep_displays_se_docs_2026-11-12`.

- **Sources:**
  - Implementer: `$S/evidence/C5/NOTES.md` (no repair, polish or correction section), the RAN / NOT RUN table `$S/evidence/C5/examples.md`, and per-page runs `$S/evidence/C5/examples/<page>/{before,after}.{R,out,rc,wall}`.
  - Verifier: `$S/evidence/C5/verify/VERIFY.md`, with its own runs in `$S/evidence/C5/verify/ex/<page>/` and extra checks in `.../verify/checks/`.
  - Main session: `$S/main/C5/{comment_gate,checkRd,commit,diffstat,document}.txt`, and the head re-runs in `$S/examples_head/`.
  - Workflow record: `$S/workflow1_result.json` entry `C5`: status READY, final verdict CONFIRMED, `repaired: false`. C5 is not in `$S/polish_result.json`.
- **Patch:** `$S/patches/C5.patch`, sha256 `91bf9790...b75eac1`, mtime 16:12:09. That is before `VERIFY.md` (16:27:32). It is the only C5 patch: there is no `.verified.patch` or `.r1.patch`.
- **The fixed copy's line numbers are not the head's.** `NOTES.md` and `examples.md` give line numbers in the fixed copy, which is the start plus C5 only. A1, A2, C3 and C4 were committed before C5 and move four files at head:
  - `R/forestsearch_cross_validation.R`: +51 lines.
  - `R/oc_analyses.R`: +4.
  - `R/subgroup_consistency_helpers.R`: +4.
  - `R/summary_utility_functions.R`: +20.
  - The other seven files have the same lines in the copy and at head.

### 1. How the examples were run

- **Extraction.** Each example was taken from its Rd with `tools::Rd2ex(..., commentDontrun = FALSE)`, so the `\dontrun{}` body is kept, and `library(forestsearch)` was put first.
  - "Start" runs use the start Rd. "Fixed" runs use the Rd that `document()` generated in the C5 copy.
  - The verifier's extracted scripts are byte-identical to the implementer's, 16/16 each way (VERIFY §2).
- **Run.** Each example ran once in a fresh `Rscript --vanilla`, capped by `timeout 120`, with the wall from `/usr/bin/time -f %e`.
  - The start and fixed runs both use one temporary library built from the C5 copy. Only the Rd text differs between them, because the parsed R code is identical (section 4). A start run is therefore the start example as a user would run it.
  - Implementer: `devtools::document()` 12.68 s; `R CMD INSTALL --no-test-load --library=$S/lib/C5` 17.15 s.
  - Verifier: `document()` 11.93 s; install into `$S/lib/C5_verify` 16.68 s; forestsearch 0.4.1.9010.
  - The machine has 128 physical cores (VERIFY §2).
- **Head re-run (main session, `$S/examples_head/`).**
  - **Library.** It was built from `$S/head_fixes`. Its `R/` and `man/` are identical to `835fac1b` (checked for this report with `git archive` and `diff -r`), and `835fac1b` has the same `R/` and `man/` as `fb80dd8f`. Its `DESCRIPTION` still says 0.4.1.9010. Install took 17.96 s.
  - **Extracted code.** The extracted head example equals the C5 run's `after.R` on 16/16 pages (`extract.txt`).
  - **Runs.** 15 pages were run. `fs_fdr_report` was not run at head: its folder holds `head.R` only.
- **Head output equals the earlier runs.** For each of the 15 pages, the output of the head run was compared with the implementer's and the verifier's fixed runs.
  - Two things were normalised first: gt's random 10-letter table ids, and the K-fold `Time: <x> minutes` line (the main session's `normalize.sh`).
  - Result: 15/15 equal (`normalized_compare.txt` says IDENTICAL). Recomputed for this report with the same script (md5): 15/15 equal to both runs.

### 2. The shared correction ("L+G+P")

On twelve pages the start example begins:

```r
#' fs <- forestsearch(gbsg,
#'   confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
#'   outcome.name = "rfstime", treat.name = "hormon", event.name = "status")
```

It is now:

```r
#' gbsg <- survival::gbsg
#' gbsg$grade3 <- as.integer(gbsg$grade == 3)
#' fs <- forestsearch(gbsg, id.name = "pid",
#'   confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
#'   outcome.name = "rfstime", treat.name = "hormon", event.name = "status")
```

Three things were wrong, so the fix has three parts:

- **L (load).** `survival` is imported, not attached, so a bare `gbsg` is `object 'gbsg' not found` after `library(forestsearch)`.
- **G (`grade3`).** `grade3` is not a `survival::gbsg` column.
- **P (`id.name = "pid"`).** The default `id.name = "id"` (`R/forestsearch_main.R:1668` at start, `:1705` at head) names no GBSG column.
  - With L and G alone, `forestsearch()` stops with ``Error in `[.data.frame`(df, , c(id.name)) : undefined columns selected``. The verifier measured this in `verify/checks/misc_check.out`, 6.21 s.
  - `pid` is GBSG's id column.

`sg_tables` and `pareto_frontier_table` take the same lines in their own layout (sections 3.14 and 3.15).

### 3. Per example

Walls are in seconds, given as implementer / verifier / head. A start run used the start Rd; a fixed run used the C5 Rd.

#### 3.1 `format_bootstrap_table`

- **Site:** `R/bootstrap_summaries_helpers.R:28` and `:32-33` → `:28-30` and `:34-35`. The example block is `:26-34` → `:26-36`.
- **Wrong:**
  - The `forestsearch(gbsg, ...)` call has all three defects of section 2.
  - The start ends `summaries <- summarize_bootstrap_results(fs$sg.harm, fs_bc)` / `format_bootstrap_table(summaries$table_data)`. But `summarize_bootstrap_results()` returns `table`, `diagnostics`, `diagnostics_table_gt`, `plots`, `timing` and `subgroup_summary`, and no `table_data` (`@return`, `R/summarize_bootstrap_results.R:28-36`).
  - `nb_boots` has no default (`R/bootstrap_summaries_helpers.R:36-37` at start).
- **Changed:**
  - L+G+P.
  - The `summarize_bootstrap_results()` line is dropped.
  - The last call is now `format_bootstrap_table(fs_bc$FSsg_tab, nb_boots = 100,` / `sg_definition = paste(fs$sg.harm, collapse = " & "))`. These are the arguments `summarize_bootstrap_results()` itself passes (`R/summarize_bootstrap_results.R:169-175`).
- **Start:** rc 1, `Error: object 'gbsg' not found`. 1.56 / 1.47 s.
- **Fixed:** rc 0, prints the gt table. Its footnote is "Identified subgroup: !{size <= 29} & {er <= 8}" (VERIFY). 17.52 / 17.23 / 17.75 s.

#### 3.2 `print.fs_kfold`

- **Site:** `R/forestsearch_cross_validation.R:2040` → `:2091-2093`. The block is `:2038-2045` → `:2089-2098`. The task's pointer was "near 2002".
- **Wrong:** L, G and P.
- **Changed:** L+G+P. `forestsearch_Kfold(fs.est = fs)` and `print(oob)` are unchanged.
- **Start:** rc 1, `object 'gbsg' not found`. 1.55 / 1.49 s.
- **Fixed:** rc 0. 67.07 / 67.25 / 68.54 s. Output:
  - `Folds: 686`, `Observations: 686`. The default `Kfolds = .fs_cv_nobs(fs.est)` is the number of observations, so this is leave-one-out, on the default 6 workers.
  - `Subgroup found in: 99.4 % of folds`, `Original subgroup: !{size <= 29} & {er <= 8}`.

#### 3.3 `print.fs_tenfold`

- **Site:** `R/forestsearch_cross_validation.R:2069` → `:2122-2124`. The block is `:2067-2074` → `:2120-2129`. The task's pointer was "near 2031".
- **Wrong / changed:** As in 3.2. `forestsearch_tenfold(fs.est = fs, sims = 5)` is unchanged.
- **Start:** rc 1, `object 'gbsg' not found`. 1.51 / 1.52 s.
- **Fixed:** rc 0. 14.40 / 14.64 / 14.61 s. Output: `Successful simulations: 5 / 5`, `Folds with subgroup identified: 74.0% (37 of 50 sim x fold pairs)`.

#### 3.4 `plot_sg_distribution`

- **Site:** `R/plot_sg_distribution.R:75-79` and `:81` → `:75-83` and `:85-92`. The block is `:73-96` → `:73-107`.
- **Wrong:**
  - Bare `gbsg`.
  - `plot_sg_distribution(fs$grp.consistency$result_new)`: a `forestsearch()` fit has no `result_new` in `grp.consistency`. Verifier `misc_check.out`: `result_new present: FALSE`. A single fit has no per-replicate table at all.
  - `results_alt` is used but never defined.
- **Changed:** The `forestsearch()` call is dropped. Two small stand-in tables are built in the two column schemas the page documents, and each is labelled as a stand-in in the example:
  - `sims_mrct` (`any_found`, `sg_found`), for `mrct_region_sims()` output: 40 of 50 found, with labels in counts 24, 10 and 6, and `"none"` for misses (as in `R/mrct_simulation.R:508`).
  - `results_alt` (`analysis`, `any.H`, `sg.def`), for `run_simulation_analysis()` output: 30 of 50 found, 30 distinct labels.
  - The `fs_rows <- ...` line and the top-K / auto-halving call are unchanged.
- **Start:** rc 1, `object 'gbsg' not found`. 1.50 / 1.49 s.
- **Fixed:** rc 0. 1.78 / 1.72 / 1.71 s.
  - One graphics-device warning: the package's own `≥` in the subtitle (`R/plot_sg_distribution.R:245`) is replaced by `>=` on the pdf device.
- **Values worked by hand** (NOTES §4, VERIFY §2a):
  - `sims_mrct`: 24/10/6 of 40 found is 60/25/15 %, all at least `min_pct = 5`, so 3 bars at threshold 5.
  - `results_alt`: each label is 1/30 = 3.33 % < 5, so one halving to 2.5, after which all 30 pass. `top_k = 15` caps the plot at 15 bars. `sgdist_fig_height(15)` = 3 + 0.35 x 15 = 8.25.
  - Measured: `n_bars` 3 / `effective_min_pct` 5; `n_bars` 15 / 2.5 / `fig_h` 8.25.
  - Neither `drafts/sgdist_attr.R` nor `verify/checks/sgdist_attr.R` saved its output. The verifier's script, re-run for this report from `$S/lib/head`, printed the same values (1.88 s).

#### 3.5 `print.fs_sg_plot`

- **Site:** `R/plot_sg_results.R:866` and `:869` → `:866-868` and `:871-872`. The block is `:864-871` → `:864-874`.
- **Wrong:**
  - L, G and P.
  - `p <- plot_subgroup(fs)`: `plot_subgroup()` is the KM helper `function(df.sub, df.subC, by.risk, confs_labels, this.1_label, top_result)` (`R/subgroup_consistency_helpers.R:967` at start). An `fs_sg_plot` is made by `plot_sg_results()` (class set at `R/plot_sg_results.R:339`).
- **Changed:** L+G+P, and `p <- plot_sg_results(fs, outcome.name = "rfstime", event.name = "status",` / `treat.name = "hormon")`.
- **Start:** rc 1, `object 'gbsg' not found`. 1.48 / 1.51 s.
- **Fixed:** rc 0. 6.49 / 6.36 / 6.44 s.
  - Prints the subgroup `!{size <= 29} & {er <= 8}`, the summary table and the HR table (H: n 76, 47 events, HR 1.722).
  - The `Available plots:` line is empty (open item 4).

#### 3.6 `plot.fs_sg_plot`

- **Site:** `R/plot_sg_results.R:912` and `:915` → `:915-917` and `:920-921`. The block is `:910-917` → `:913-923`.
- **Wrong / changed:** As in 3.5.
- **Start:** rc 1, `object 'gbsg' not found`. 1.57 / 1.51 s.
- **Fixed:** rc 0. Prints `No plots available` (open item 4). 7.06 / 6.16 / 6.42 s.

#### 3.7 `print.fs_forestplot`

- **Site:** `R/render_forestplot.R:372` → `:372-374`. The block is `:370-378` → `:370-380`.
- **Wrong:** L, G and P. The bare `gbsg` is also the data argument of `plot_subgroup_results_forestplot(list(fs.est = fs), gbsg, ...)` (`:375` → `:377`).
- **Changed:** L+G+P. The forest-plot lines are unchanged and now resolve.
- **Start:** rc 1, `object 'gbsg' not found`. 1.56 / 1.53 s.
- **Fixed:** rc 0. Prints the forest-plot object's help lines. 6.39 / 6.19 / 6.44 s.

#### 3.8 `plot.fs_forestplot`

- **Site:** `R/render_forestplot.R:410` → `:412-414`. The block is `:408-416` → `:410-420`.
- **Wrong / changed:** As in 3.7.
- **Start:** rc 1, `object 'gbsg' not found`. 1.51 / 1.49 s.
- **Fixed:** rc 0, and `Rplots.pdf` is written. 6.44 / 6.32 / 6.76 s.

#### 3.9 `run_subgroup_sims`

- **Site:** 22 lines are inserted after `R/run_subgroup_sims.R:652` (`\dontrun{`), becoming `:653-674` at head, and `:663` → `:685`. The block is `:651-664` → `:651-686`. The task's pointer was "near 660".
- **Wrong:**
  - `dgm_uniform`, `subgroups`, `cal_uniform`, `k_treat_uniform` and `gbsg` are never defined.
  - `saveRDS(sims, "results/extreme_sims_fixed_1000_payload.rds")` writes under the working directory. This fails where no `results/` exists.
- **Changed:** The objects are built as in `vignettes/articles/extreme_subgroups.qmd` (chunks build-dgm, calibrate and define-uniform-subgroups), shortened:
  - `gbsg` plus `time_months`.
  - `base_args`, the uniform-effect GBSG DGM. Arguments the article sets to their defaults are omitted.
  - `k_treat_uniform <- calibrate_k_treat(target_hr_overall = 0.70, base_args = base_args, use_ahr = TRUE)`.
  - `dgm_uniform <- do.call(generate_aft_dgm_flex, ...)`.
  - `cal_uniform <- calibrate_cens_adjust(dgm_uniform, target = "rate", ..., analysis_time = 84, max_entry = 24, ..., baseline = "fixed")`.
  - A three-subgroup list (`meno == 0`, `age <= age_med`, `random60 == 1`), verbatim from the article's `:978`, `:983` and `:1050`.
  - `saveRDS()` now writes to `file.path(tempdir(), ...)`.
  - The `run_subgroup_sims()` call itself is unchanged, `n_sims = 1000` included.
  - The omitted calibration arguments are defaults: `k_treat_range = c(-5, 5)` and `tol = 1e-6` (`R/calibrate_k_treat.R:98-99`); `seed = 42`, `interval = c(-3, 3)` and `tol = 1e-4` (`R/simulate_from_dgm.R:795-797`).
- **Start:** rc 1, `object 'dgm_uniform' not found`. 1.53 / 1.47 s.
- **Fixed:** rc 0. The example prints nothing (empty `after.out` in all three runs). 13.01 / 12.93 / 13.35 s.
- **Verifier's check** (`verify/checks/rss_check.out`: the example plus print lines, rc 0, 12.90 s):
  - `k_treat_uniform` 0.956, AHR 0.70.
  - `cens_adjust` 0.0412, with reference rate 0.5756 and simulated rate 0.5758.
  - 1000 sims x 3 subgroups, no NA.
  - mHR 0.72 in each subgroup.
  - The file was written to `tempdir()`.

#### 3.10 `sg_consistency_out`

- **Site:** `R/subgroup_consistency_helpers.R:1031` → `:1035-1037`. The block is `:1027-1035` → `:1031-1041`.
- **Wrong / changed:** L, G and P; fixed by L+G+P. `fs$grp.consistency$out_sg` is unchanged.
- **Start:** rc 1, `object 'gbsg' not found`. 1.53 / 1.48 s.
- **Fixed:** rc 0. 6.51 / 6.28 / 6.50 s.

#### 3.11 `evaluate_subgroup_consistency`

- **Site:** `R/subgroup_consistency_helpers.R:1361` → `:1367-1369`. The block is `:1357-1364` → `:1363-1372`.
- **Wrong / changed:** L, G and P; fixed by L+G+P.
- **Start:** rc 1, `object 'gbsg' not found`. 1.56 / 1.56 s.
- **Fixed:** rc 0. 6.38 / 6.26 / 6.66 s.

#### 3.12 `evaluate_consistency_twostage`

- **Site:** `R/subgroup_consistency_helpers.R:1645` → `:1653-1655`. The block is `:1641-1649` → `:1649-1659`.
- **Wrong / changed:** L, G and P; fixed by L+G+P. `use_twostage = TRUE` is kept.
- **Start:** rc 1, `object 'gbsg' not found`. 1.52 / 1.57 s.
- **Fixed:** rc 0. 6.34 / 6.32 / 6.39 s.

#### 3.13 `fs_fdr_report`

- **Site:** 7 lines are inserted after `R/fs_fdr_report.R:151`, becoming `:152-158` at head; `:165` → `:172` and `:170` → `:177`. The block is `:150-172` → `:150-179`. The task's pointer was "near 148-164"; the `\dontrun{}` body is `:151-171` at start.
- **Wrong:**
  - `gbsg_df` is never built: `Error: object 'gbsg_df' not found`.
  - Both calls pass `n_workers = 100`. parallelly's hard limit is 3 x the available cores, so 100 workers is an error on a machine with fewer than 34 cores. The verifier checked this with `options(mc.cores = 8)` (parallelly 1.48.0): "100 localhost parallel workers with only 8 CPU cores ... hard limit is set to 300%".
- **Changed:**
  - New lines build `gbsg_df <- within(survival::gbsg, {` / `id <- pid` / `time <- rfstime / 30.4375` / `treat <- hormon` / `})`, plus a comment and a `#'` spacer. Time is in months.
  - `, n_workers = 100` is dropped from both calls. The default is then `max(1L, ceiling(0.90 * (parallel::detectCores(logical = FALSE) - 1L)))` (`R/fs_fdr_report.R:232-234` at start, `:239-241` at head), which is 115 on this machine.
- **Start:** rc 1, `object 'gbsg_df' not found`. 1.52 / 1.45 s.
- **Fixed:** **NOT RUN, over the cap.** Killed by `timeout` at 120.76 s (rc 124) in both runs.
  - By then the example had built `gbsg_df`, fitted the model (MR on, "MR harm confirmation: HR = 1.002") and built both nulls (`H0-harm k_inter = 0.5679`).
  - It had printed `Running 1000 reps/null x 2 null(s) at n = 686 on 115 workers ...`.
  - The example runs 2 calls x 2 nulls x 1000 null `forestsearch()` fits with MR.
- **Smoke checks** (these are not the example):
  - Implementer `examples/fs_fdr_report/smoke.*`: the example's lines with the first call at `nsims = 4, n_workers = 4`. rc 0, 63.20 s, verdict CLEAN, six `coxph.fit` "Loglik converged before variable 1" warnings.
  - Verifier `fdr_smoke1`: the same reduction. rc 0, 63.43 s, CLEAN, n_valid 4/4 per null, 6 `fdr` rows, realised H0-harm HR_H 1.002.
  - Verifier `fdr_smoke2`: the last call, `fs_fdr_report(fit, gbsg_df, n = 700, ...)`, at `nsims = 4, n_workers = 4`. rc 0, 21.11 s, "n per replicate: 700", CLEAN.

#### 3.14 `sg_tables`

- **Site:** `R/summary_utility_functions.R:742`, plus a line inserted after `:744` → `:762-763` and `:766`. The block is `:740-755` → `:760-777`.
- **Wrong:** `library(survival)` finds `gbsg`, but `grade3` is not a column. The start run (rc 0) shows:
  - `forestsearch()` warns "Error in get_FSdata: confounders.name references columns not in df.analysis: 'grade3'." and "FSdata failure - returning NULL result".
  - `sg_tables()` warns "no harm subgroup was identified (fs$sg.harm is NULL); returning empty tables."
  - `tabs$tab_estimates` prints `NULL`.
- **Changed:** `library(survival)` is replaced by the L and G lines, and `id.name = "pid",` is added as its own argument line (head `:766`).
- **Start:** rc 0, no table. 1.56 / 1.51 s.
- **Fixed:** rc 0, prints the gt estimates table. 6.89 / 6.76 / 6.89 s.

#### 3.15 `pareto_frontier_table`

- **Site:** `R/pareto_frontier_table.R:125` → `:125-128`. The block is `:121-136` → `:121-139`.
- **Wrong:** The placeholder `fs <- forestsearch(gbsg, ...)` gives `Error: '...' used in an incorrect context`. The existing `gbsg <- survival::gbsg` at `:124` was already right.
- **Changed:** `:124` is kept. G is added, and the full call uses P with the confounders and the outcome, treatment and event names of the other GBSG examples.
- **Start:** rc 1, `'...' used in an incorrect context`. 1.50 / 1.45 s.
- **Fixed:** rc 0. 10.06 / 9.88 / 10.66 s.
  - A one-member frontier: `!{size <= 29}`, `{er <= 8}`, N 76, E 47, HR 1.72.
  - ggplot prints "Each group consists of only one observation" from `plot_pareto_frontier()`.

#### 3.16 `summarize_simulation_results`

- **Site:** `R/oc_analyses.R:414` → `:418-424`. The block is `:411-416` → `:415-426`.
- **Wrong:** `run_simulation_analysis(dgm, nsim = 10)` fails because there is no `nsim` argument: `Error in run_simulation_analysis(dgm, nsim = 10) : unused argument (nsim = 10)`. The function runs one replicate per call and needs `sim_id` and `n_sample` (`R/run_simulation_analysis.R:389-419`).
- **Changed:** Three comment lines, then:
  ```r
  sim_res <- data.table::rbindlist(lapply(1:10, function(i)
    run_simulation_analysis(sim_id = i, dgm = dgm, n_sample = 500,
      confounders_base = c("z1", "z2", "z3", "z4", "z5", "size", "grade3"))),
    fill = TRUE)
  ```
  - `confounders_base` is the one used in `vignettes/articles/paper_simulations.qmd:883`.
  - The default `c("v1", ..., "v7")` (`R/run_simulation_analysis.R:398`) holds factors in the `setup_gbsg_dgm()` data. With it, every FS replicate fails inside: the verifier's `ssr_default.out` (2 replicates, 4.13 s) shows 8 warnings, among them "Error in get_FSdata: Calling var(x) on a factor x is defunct" and "FSdata failure - returning NULL result", and `any.H` is 0 on every row.
- **Start:** rc 1, `unused argument (nsim = 10)`. 1.64 / 1.60 s.
- **Fixed:** rc 0. Prints the FS / FSlg / GRF table, for example `any.H 0.100 0.300 0.200`. 18.08 / 17.77 / 18.81 s.

#### 3.17 All runs

| Rd page | start rc, impl / ver wall | fixed rc | fixed wall impl / ver / head | result |
|---|---|---|---|---|
| `format_bootstrap_table` | 1, 1.56 / 1.47 | 0 | 17.52 / 17.23 / 17.75 | RAN |
| `print.fs_kfold` | 1, 1.55 / 1.49 | 0 | 67.07 / 67.25 / 68.54 | RAN |
| `print.fs_tenfold` | 1, 1.51 / 1.52 | 0 | 14.40 / 14.64 / 14.61 | RAN |
| `plot_sg_distribution` | 1, 1.50 / 1.49 | 0 | 1.78 / 1.72 / 1.71 | RAN |
| `print.fs_sg_plot` | 1, 1.48 / 1.51 | 0 | 6.49 / 6.36 / 6.44 | RAN |
| `plot.fs_sg_plot` | 1, 1.57 / 1.51 | 0 | 7.06 / 6.16 / 6.42 | RAN |
| `print.fs_forestplot` | 1, 1.56 / 1.53 | 0 | 6.39 / 6.19 / 6.44 | RAN |
| `plot.fs_forestplot` | 1, 1.51 / 1.49 | 0 | 6.44 / 6.32 / 6.76 | RAN |
| `run_subgroup_sims` | 1, 1.53 / 1.47 | 0 | 13.01 / 12.93 / 13.35 | RAN |
| `sg_consistency_out` | 1, 1.53 / 1.48 | 0 | 6.51 / 6.28 / 6.50 | RAN |
| `evaluate_subgroup_consistency` | 1, 1.56 / 1.56 | 0 | 6.38 / 6.26 / 6.66 | RAN |
| `evaluate_consistency_twostage` | 1, 1.52 / 1.57 | 0 | 6.34 / 6.32 / 6.39 | RAN |
| `fs_fdr_report` | 1, 1.52 / 1.45 | 124 | 120.76 / 120.76 / not run | NOT RUN (over the cap) |
| `sg_tables` | **0** (no table), 1.56 / 1.51 | 0 | 6.89 / 6.76 / 6.89 | RAN |
| `pareto_frontier_table` | 1, 1.50 / 1.45 | 0 | 10.06 / 9.88 / 10.66 | RAN |
| `summarize_simulation_results` | 1, 1.64 / 1.60 | 0 | 18.08 / 17.77 / 18.81 | RAN |

- **Start:** 15 runs stopped with an error. 11 gave `object 'gbsg' not found`; the other four gave `dgm_uniform`, `gbsg_df`, `'...'` and `nsim`. `sg_tables` ended with rc 0 and no table. Walls 1.48-1.64 s (implementer) and 1.45-1.60 s (verifier).
- **Fixed:** 15 RAN and 1 NOT RUN. No fixed output contains "Error" (VERIFY §2). The only warning is the `plot_sg_distribution` device warning.
- **Agreement:** The verifier's walls are within 0.90 s of the implementer's; the largest gap is `plot.fs_sg_plot`, 7.06 vs 6.16. The head walls are within 1.47 s of the implementer's; the largest gap is `print.fs_kfold`, 68.54 vs 67.07.

### 4. Gates

- **Comment-only (the task's command), every changed `R/` file:**
  - Main session `main/C5/comment_gate.txt`: "comment-only" for all 11 files, `COMMENT GATE: PASS`.
  - Re-run for this report on `git diff -U0 ba7c4ea2 5073fd67`: empty output for all 11.
  - 0 added lines are non-ASCII: implementer `gate.txt`, verifier `gate.txt`, and re-checked on git. The verifier also found 0 blank lines added in `R/`.
- **Only `@examples` lines changed:** The verifier's `examples_only.py` found `TOTAL outside @examples: 0`. `R/` is +99/-26, the same as `git diff --numstat ba7c4ea2 5073fd67 -- R`.
- **Parsed code identical:** `parse(keep.source = FALSE)` start vs fixed, then `identical(lapply(exprs, deparse))`, is TRUE for 11/11 files with equal expression counts (implementer and verifier `parse_identity.txt`).
- **`checkRd`:**
  - Main session `main/C5/checkRd.txt`: 16/16 `OK (0 problems)`.
  - Implementer and verifier, with `encoding = "UTF-8", def_enc = TRUE` (the package declares UTF-8): 16/16 OK.
  - A plain `checkRd(f)` gives one note, "Non-ASCII contents without declared encoding", at `evaluate_consistency_twostage.Rd:96` and `evaluate_subgroup_consistency.Rd:82`. The cause is the accent in "León" in a paragraph outside the examples. The start pages give the same note.
- **`man/` equals the verified patch:**
  - **Commit vs patch.** For each of the 27 files, the commit's changed lines (`git diff ba7c4ea2 5073fd67`) were listed with their file and compared with those of `patches/C5.patch` (built against the start). They are the same 250 lines (+198/-52); the `diff` is empty (checked for this report).
  - **Verifier, `document()`.** `document()` on start + patch rewrote exactly the 16 pages, byte-identical to the patch's `man/`. On a start control copy it changed 0 files. On all 16 pages the Rd text outside `\examples{}` is identical start vs fixed; the negative control (`sg_tables` `\examples`) compares as different.
  - **Head.** 15 of the 16 pages are byte-identical (md5) to the C5 copy after `document()`.
    - `sg_tables.Rd` differs only by C4's `\details{}` block: 20 lines at head `:83-102` (`examples_head/sg_tables_rd.diff`). This block is also the whole 20-line `b06de4e8 → ba7c4ea2` change to that page.
    - The extracted example code is identical on all 16 pages (`examples_head/extract.txt`).
  - **No `document()` record for C5.** The main session's `main/C5/document.txt` is empty, so there is no wall for a `document()` run in this commit (C4's `main/C4/document.txt` records `real 0m11.765s`).
- **`\dontrun{}` kept** on all 16 pages (16 occurrences in the 16 head pages).
- **Scope:** `main/C5/diffstat.txt` and `git diff --stat ba7c4ea2 5073fd67` both give 27 files changed, 198 insertions(+), 52 deletions(-). That is 11 `R/` files and 16 `man/` pages, and nothing under `tests/` or `vignettes/`.
- **Existing tests:** None were run, because no executable line changed.

### 5. Departures from the task's wording; task text vs code

- **`id.name = "pid"`.** The task names loading and `grade3` only. Without `id.name`, `forestsearch()` stops (section 2), so it is added to every GBSG `forestsearch()` call.
- **`plot_sg_distribution`.** No simulation is run. A fit has no per-replicate table, so stand-in tables are used.
- **`format_bootstrap_table`.** A line that uses a non-existent element is dropped, and the required `nb_boots` is passed.
- **`print.fs_sg_plot` / `plot.fs_sg_plot`.** A wrong function call is replaced by the constructor.
- **`run_subgroup_sims`.**
  - Three subgroups instead of the article's list.
  - Arguments left at their defaults are omitted.
  - The output goes to `tempdir()`.
- **`fs_fdr_report`.** `n_workers = 100` is dropped.
  - The task's C5 text says nothing about workers (no "worker" in the task). NOTES §8 says "workers lowered as the task allows".
  - On this machine the default (115) is more than 100 (verifier minor 1). The reason that holds is the parallelly hard-limit error below 34 cores.
- **`summarize_simulation_results`.** Besides removing `nsim`, the example passes `confounders_base` (section 3.16).
- **Pointers.** The task's "near 2002, 2031" are `2040` and `2069` at start. "near 148-164" is the block `:150-172` (`\dontrun{}` body `:151-171`). "near 660" is inside the block `:651-664`. The other pointers match the start lines.
- **Part 3 report vs code.** Part 3 says (`REPORT_cran_prep_docs_legacy_2026-11-12.md:364`, "by code reading ... not run") that `forestsearch()` stops in the `sg_tables` example. The code warns and returns a NULL result, and `sg_tables()` returns empty tables (section 3.14).
- **Same pattern, not edited.** `R/subgroup_consistency_helpers.R:1223` (head `:1229`) has the bare-`gbsg` pattern in an `@noRd` block (`:1241`). It is on no help page.
- **Commit subject.** The commit subject is "docs: make sixteen \dontrun examples runnable as written (part 5, C5)". The implementer proposed "docs: make sixteen GBSG and simulation help-page examples runnable as written".

### 6. Verification record

- Verifier `VERIFY.md`: **CONFIRMED** in round 1, with no blocking problem. There was no repair round and no polish.
- Minor points:
  1. `examples.md` and NOTES §8 call dropping `n_workers = 100` "workers lowered". On this 128-core machine the default is 115. The fix is still correct: 100 workers is a hard error below 34 cores.
  2. The proposed NEWS bullet says all sixteen examples "failed" and "now run", which does not fit `sg_tables` (it ended with rc 0) or `fs_fdr_report` (run only at reduced size). Its list of causes omits `plot_subgroup(fs)`, `summaries$table_data` and `n_workers = 100`.
     - The committed bullet (`fb80dd8f:NEWS.md:556-559`) is a different text. It still says "Sixteen `\dontrun{}` examples now run as written", and its list still omits those three.
  3. `plot.fs_sg_plot`'s `@return` (start `R/plot_sg_results.R:906-908`, head `:909-911`) promises "a message names the `plot_sg_results()` call", but the corrected example prints "No plots available". C5 may not edit `@return`, so this is open item 4.
  4. The implementer's open items were re-checked as stated.

### 7. Open items

1. **`fs_fdr_report` not run in full.** It goes over the cap (2 calls x 2 nulls x 1000 MR replicates) and was not re-run at head. Its comment "Reproduce the published reference (n = 700): A ~ 0.138 (hom), 0.278 (harm)" (start `R/fs_fdr_report.R:169`, head `:176`) is not verified.
2. **Three GBSG examples outside the C5 list stop.** Their `forestsearch()` call has no `id.name`, and each stops with ``Error in `[.data.frame`(df, , c(id.name)) : undefined columns selected``. The implementer measured this with the examples' prefixes, `evidence/C5/outside_list/*.prefix.out`. The examples:
   - `summarize_bootstrap_results`: `R/summarize_bootstrap_results.R:52-57`, unchanged at head.
   - `forestsearch_KfoldOut`: `R/forestsearch_cross_validation.R:1331-1335` → `:1345-1349`.
   - `summarize_cv_results`: `R/summarize_cv_results.R:148-159`, unchanged at head.
3. **`run_simulation_analysis()`'s default confounders fail on a `setup_gbsg_dgm()` DGM.** The default `confounders_base = c("v1", ..., "v7")` (`R/run_simulation_analysis.R:398`, unchanged) holds factors there, so every FS replicate fails inside with a warning (section 3.16). The function's own example uses the default and was not run.
4. **`plot_sg_results()` always returns `plots = list()`.** `plot_sg_km()`, `plot_sg_forest()` and `plot_sg_summary_panel()` return `invisible(NULL)` (`R/plot_sg_results.R:661`, `:775`, `:849`, the same at start and head), and assigning NULL drops the element. So:
   - `print.fs_sg_plot` shows an empty "Available plots:" line.
   - `plot.fs_sg_plot` prints "No plots available", while its `@return` (head `:909-911`) describes a different message.
   - The plots are drawn while `plot_sg_results()` runs.
5. **Other examples use placeholder data.** Examples with `my_data`, `mydata` or `trial_data` (for example in `R/plot_sg_results.R:95`, `R/forestsearch_methods.R:236`, `R/bootstrap_dofuture_main.R:164` at head) were not in this item and were not checked here.
6. **NEWS overstates the result.** `fb80dd8f:NEWS.md:556-559` says "Sixteen `\dontrun{}` examples now run as written". The evidence has 15 RAN and `fs_fdr_report` NOT RUN in full: its lines ran only at `nsims = 4, n_workers = 4`.

### 8. Evidence notes (conflicts and gaps)

- **"Workers lowered."** NOTES §8 calls this "as the task allows", but the task text has no worker wording. "Lowered" is also wrong for the measured runs, which printed "on 115 workers".
- **`run_subgroup_sims` summary.** `examples.md` gives the result as "RAN (1000 simulations; summary: 3 subgroups, mHR 0.72)". The example's own output is empty: implementer `after.out` is 0 bytes, and the verifier's and head outputs are empty too. The summary values are in the verifier's `checks/rss_check.out`, which adds print lines. The implementer's `drafts/rss.R` has no saved output.
- **`plot_sg_distribution` attributes.** The measured values (3 / 5; 15 / 2.5 / 8.25) are stated in NOTES §4, `examples.md` and VERIFY §2a, but no saved output exists for either `sgdist_attr.R`. Re-run for this report: same values.
- **No C5 `document()` record.** `main/C5/document.txt` is empty (0 bytes), so there is no main-session record of a `document()` run or its wall for C5. Nothing contradicts the commit's `man/`, which equals the verified patch (section 4).
- **Line-number frames.** NOTES §2 and `examples.md` give fixed-copy lines. They differ from head in four files (front matter). This is not a contradiction.

## C6 — comments on test expectations (detail)

Commit `ffc55486` (parent `5073fd67`). Line numbers: start `b06de4e8` and head `fb80dd8f`.

- **Where the requests come from:**
  - Part 4 report §8, OPEN ITEMS 17 (`git show b06de4e8:dev/reports/REPORT_cran_prep_code_fixes_2026-11-12.md`) has 13 requests in nine test files.
  - The F6 "Verification" table in part 4's detail file, problem 2 (`..._detail.md:1335`), adds three optional ones: `:65`, `:71` and the survival case of test 4. All 16 are done.
- **Line numbers:**
  - Item 17 gives its line numbers at each fix's commit. Each of the nine files is the same blob at its part-4 fix commit and at `b06de4e8` (`git rev-parse`, 9 of 9), so item 17's lines are start lines.
  - No part-5 commit before or after C6 touches the nine files. Both `git diff --stat b06de4e8 5073fd67` and `git diff --stat ffc55486 fb80dd8f` are empty for them. So head lines are `ffc55486` lines.
- **Commit equals patch:** `git diff 5073fd67 ffc55486` equals `patches/C6.patch` (the polished patch). The only difference is how long the abbreviated blob hashes on the `index` lines are; the hashes themselves match. `main/C6/patch_equal.txt` is empty.

### Requests, comments and polish

All paths are `tests/testthat/`. "Checks" are the expectations the comment is about. The text in the "Comment added" column is the committed comment, with its lines joined.

| # | Fix | Request (source) | Start `b06de4e8` | Head `fb80dd8f`: comment; checks | Comment added (quoted) | Polish |
|---|---|---|---|---|---|---|
| 1 | F1 | The label-order checks pass on both trees; say so, and say they anchor the row order (item 17; detail F1 problem 1) | `test-bootstrap-summary-bins.R:76-83` (comment `:75`) | `:76-81`; checks `:82-89` | "These two label checks pass before and after the A-147 fix (on the pre-fix code they pass while the per-value, Count and Percent checks fail): the labels are the \`labels\` of the two cut() calls, and table() of a factor keeps every level, zero counts included, in level order; the fix changed only right = FALSE.  They anchor the row order that the positional Count and Percent checks below rely on." | none |
| 2 | F5 | Derivation comment for 0 warnings, non-NULL value, class `hte_test` (item 17; detail F5 round 2 problem 1) | `test-hte-crump-confounders.R:97-99` | `:97-105`; checks `:106-108` | "By hand: every label in ev is a cut label "age <= v" or "bmi <= v" (precondition above); all.vars() of the parsed label is "age" or "bmi", both columns of df, so no label is skipped (no "cannot be resolved" warning) and X gets two columns (no "No covariate found" warning, no NULL).  hte_test is test_hte()'s return value, whose class test_hte() sets to "hte_test". Before the A-151 fix the labels were looked up as column names, none matched, and the call returned NULL without a warning: the warning count passed then too; the NULL and class checks did not." | none |
| 3 | F6 | A note that test 2 passes before and after (item 17; detail F6 problem 1) | `test-lasso-intercept.R:88` (test 2, `:88-136`) | `:94-99`, above `test_that` (now `:100`; test 2 `:100-148`) | "This test passes before and after the fix (on the pre-fix code every check in it passes).  It replaces lasso_selection(), the function the fix changes, so it cannot show the defect, and get_FSdata()'s rule "if LASSO selects nothing, keep the full pool" was already right.  It is kept as a lock on that rule.  The red for this case is the next test, which runs the real lasso_selection()." | none |
| 4 | F6 | Optional: `:65` passes on both trees (detail F6 problem 2; not in item 17) | `test-lasso-intercept.R:65` | `:65-66`; check `:67-68` | "This check passes before and after the fix, which leaves \`coefficients\` as it was (every row of coef(fit))." | none |
| 5 | F6 | Optional: `:71` passes on both trees (detail F6 problem 2; not in item 17) | `test-lasso-intercept.R:69-71` (comment `:68-69`) | `:70-75` (start `:69` extended); checks `:76-77` | Appended to "...neither \`selected\` nor \`omitted\`.": "The \`selected\` check is the red: before the fix the non-zero fitted intercept was listed in \`selected\`. The \`omitted\` check passes before and after: \`omitted\` held only the zero coefficients, and the fitted intercept is non-zero in these three fits." | none |
| 6 | F6 | Optional: the survival case of test 4 passes on both trees (detail F6 problem 2; not in item 17) | `test-lasso-intercept.R:180-214` (header `:181-185`; checks `:208`, `:209-212`) | `:198-205`; checks `:228` (`expect_length`), `:229-232` (`expect_match`) | "The survival case passes before and after the fix (the line always said "Cox-LASSO" before it).  Before the fix the binary fit printed no LASSO summary line at all: LASSO selected only "(Intercept)", get_FSdata() built no candidate cuts and forestsearch() returned before the summary, so both checks fail on the empty grep() result.  The continuous and count fits printed "Cox-LASSO selected: 4 of 5" and "2 of 5 candidate factors" (old label, intercept counted): one line, so only the match check fails for them." | **Rewritten.** Verified text (`:198-200`): "...(the line always said "Cox-LASSO" before it); the binary, continuous and count labels are the red." Verifier minor 2: before the fix the binary fit printed no line (length 0). |
| 7 | F8 | The exact-duplicate guard also holds through `unique()` and does not by itself pin the `%in% confs` clause (item 17; detail F8 problem 1) | `test-duplicate-cut-names.R:162-178` (checks `:168`, `:177`) | `:163-172`; checks `:178`, `:187` | "What this block pins: the outcome "the duplicate appears once", not the \`grf_cuts %in% confs\` / \`dina_cuts %in% confs\` clause of the merges by itself.  Both merges end in confs <- unique(c(confs, grf_cuts_keep)) or unique(c(confs, dina_cuts_keep)) (get_fsdata.R), which drops an exact duplicate anyway: a copy of get_fsdata.R with that clause removed from both merges still passes every check in this file.  The clause's only other effect is that the duplicate stays out of grf_cuts_keep / dina_cuts_keep, and so out of the \`protect\` list passed to collapse_redundant_cuts(), which runs only with collapse_cuts = TRUE (FALSE here)." | none |
| 8 | F9 | Say why df_super's raw copies and `zcens_age` / `zcens_prior_treat` are absent (item 17; detail F9 problem 2) | `test-mrct-default-confounders.R:40-49` (`expected` `:50`) | `:50-56`; `expected` `:57` | "The other columns of df_super are outside the default for another reason: both functions keep only names matching "^z_\|^ecog\|^strat", so the raw copies age, bm, male, histology, prior_treat and regA, and zcens_age and zcens_prior_treat ("zcens_" does not match "^z_"), never enter it, nor do id, y, treat, event, flag_harm, lin_pred_\*, theta_\* and loghr_po. Among the matching names only the three spline helpers are removed, by exact name." | none |
| 9 | F9 | `expect_length(line, 1L)` needs "(one replicate, verbose_n_sims = 1)" (item 17; detail F9 problem 2) | `test-mrct-default-confounders.R:63` | `:70-72`; check `:73` | "Exactly one line: one replicate (n_sims = 1), and with details = TRUE the line is printed for each replicate sim <= verbose_n_sims (= 1) that runs forestsearch() (inside if (!is.na(hr_test)) in mrct_region_sims())." | none |
| 10 | F9 | Mark `:82` as a regression guard that passes on both trees (item 17; detail F9 problem 1) | `test-mrct-default-confounders.R:82` | `:92-94`; check `:95` | "A regression guard, not a check of the fix: it passes before and after the fix here, since this trial's identifier (!{z_male} on both trees) uses no spline helper column." | none |
| 11 | F10 | Derivation for the tree-path existence checks (item 17; detail F10 problem 1) | `test-grf-count-orientation.R:103-106` | `:103-112`; checks `:113-116` | "Existence checks, so that a missing subgroup fails here and not as an empty (NaN) share below.  By hand: in each orientation one region has a positive harm score -CATE on the forest's scale (header: 1.10 in z == 1 with adverse_outcome = TRUE, 0.93 in z == 0 with FALSE), so the policy tree that gives control to that region yields a subgroup of about 500 subjects, at least n.min (60, the default) and fewer than n = 1000, with a score above dmin.grf (0, the default); sg.harm.id is NULL, with no members, only when no tree depth passes these (.select_best_sg_glm()). They pass before and after the fix: before it both fits found the same z == 0 subgroup, the wrong region for adverse_outcome = TRUE." | **Last sentence rewritten** (`:111-112`). Verified text: "...before it a subgroup was found too, in the wrong region." Verifier minor 1: true only for the TRUE fit. |
| 12 | F10 | Derivation for the frontier-path existence checks (item 17; detail F10 problem 1) | `test-grf-count-orientation.R:129-132` | `:139-145`; checks `:146-149` | "Existence checks, as in the tree test.  By hand: a frontier candidate is eligible when it has at least n.min (60) and at most n - 1 (999) members (.grf_dr_candidates()) and a harm-effect of at least dmin.grf (0) (.grf_frontier_select()); the candidate on z has 500 members and a harm-effect of about 6 (z > 0, TRUE) or 2 (z <= 0, FALSE) (header), so a subgroup is selected in both orientations.  They pass before and after the fix: before it both fits selected z <= 0." | **`:140-145` rewritten and reflowed.** Verified text: "...eligible when it has at least n.min (60) members and a harm-effect of at least dmin.grf (0) (.grf_dr_candidates(), .grf_frontier_select())...". Verifier minor 3: the upper bound `nS <= n - 1` was missing. |
| 13 | F12 | The header justifies the 1e-4 MR-limit tolerance by a z difference (1.96 vs `qnorm(0.975)`) that the test does not have (item 17; detail F12 problem 2a) | `test-display-orientation.R:27-30` (1e-4 check `:151`) | `:27-38`, rewritten; 1e-4 check `:166` | "Tolerances.  For printed numbers the tolerance is the display rounding, half a unit of the last printed digit: labels at 3 decimals, 5e-4 (the labels' limits and the hand limits both use z = qnorm(0.975)); the MR table at digits = 4, 5e-5 for the estimate.  The MR limits checked (the complement row's) are printed from est +/- 1.96 \* se with the literal z = 1.96 (analyze_subgroup_glm()), and the hand limits below use the same 1.96, so no z difference enters: only the 4-decimal rounding (at most 5e-5) separates them, and the 1e-4 used is twice that bound.  (With qnorm(0.975) in the hand limits the difference would grow by (1.96 - qnorm(0.975)) \* se, below 1e-5 at the complement's se of about 0.12.)  The forest plot's est and se are unrounded numbers, not printed text: see the comment at that check." | none |
| 14 | F12 | The forest-plot 1e-10 tolerance (unrounded columns) is not explained in the test (item 17; detail F12 problem 2b) | `test-display-orientation.R:127-130` | `:135-141`; checks `:142-145` | "dt$est and dt$se are the estimator's unrounded estimate and standard error (MD is on the identity scale, so est is not transformed), not printed text; they and the hand formulas compute the same OLS quantities, so they differ only by floating-point error.  1e-10 is therefore a numerical tolerance (testthat 3e: a relative difference, as in all.equal()), far above that error and far below the 0.19 to 2.2 by which a sign flip moves these estimates (twice \|est\|)." | none |
| 15 | F13 | No comment for `converged`, `method_used`, `measure`, `n0`/`n1` in the second block (item 17; detail F13 flag) | `test-rd-standard-error.R:78-81` | `:78-83`; checks `:84-87` | "As in the first test: tier 2 returns converged = TRUE and, with no adjustment covariates, method_used "logistic_margins" (tier 3 would give FALSE and "raw_means"); measure is "RD", the effect_measure requested; n0 and n1 are the arm sizes, 40 each by construction of .rdse_slice(), counted as integers (sum() of a logical).  These four checks pass before and after the fix, which changes only the standard error." | none |
| 16 | F14 | The comment cites only the "Resolved defaults" table, which is stated for neither threshold supplied; cite the remap "to default" step (item 17; detail F14 problem 1) | `test-dina-gaussian-floor.R:159-161` | `:159-168`, rewritten (the sentence ending "The" runs on into `:169`) | "Expected recorded values: a supplied RD/IRD effect.threshold above 1 is remapped "to default" (forestsearch_main.R, the branch \`user_set_threshold && effect_threshold > 1.0\`, which warns and sets switch(effect_measure, RD = 0.05, IRD = 0.01)).  These are the defaults ?forestsearch "Resolved defaults, per estimand" lists for a run with neither threshold supplied, RD -> c1 = 0.05, IRD -> c1 = 0.01 (identity scale, so threshold_config$screening is that value itself, not its log); here they come from the remap step, since a threshold is supplied. .gf_capture() muffles the remap warning and the earlier validator warning on the same value." | none |

### Existing comment text that was rewritten (8 deleted lines)

| File | Start lines | Old text (quoted) | Head lines |
|---|---|---|---|
| `test-lasso-intercept.R` | `:69` (1 line) | "# neither \`selected\` nor \`omitted\`." | `:71-75` (the sentence is kept, and four sentences follow) |
| `test-display-orientation.R` | `:27-30` (4 lines) | "Tolerances are the display rounding: half a unit of the last printed digit (labels at 3 decimals: 5e-4; the MR table at digits = 4: 5e-5 for the estimate, 1e-4 for its limits, whose z is the literal 1.96 rather than qnorm(0.975), a difference below 1e-5 at these standard errors)." | `:27-38` |
| `test-dina-gaussian-floor.R` | `:159-161` (3 lines) | "remaps.  Expected recorded values, from ?forestsearch "Resolved defaults, per estimand": RD -> c1 = 0.05, IRD -> c1 = 0.01 (identity scale, so threshold_config$screening is that value itself, not its log).  The" | `:159-168` |

1 + 4 + 3 = 8 lines, which matches the commit's 8 deletions.

### How each comment's facts were checked

- **Method:**
  - Each "passes before and after" claim was checked on a pre-fix copy. That copy is the start tree with the fix commit's `R/` diff reverse-applied (`git show <fix> -- R/ | git apply -R`) for F1 `0dc99947`, F5 `4aaab2f0`, F6 `f5b65516`, F9 `a30d6e7c`, F10 `f2e8de47` and F13 `89ad9c19`.
  - Each expectation's outcome was recorded with a `ListReporter`. The implementer's results are in `evidence/C6/prefix/<F>.per_expect.txt` and the verifier's in `evidence/C6/verify/prefix/`.
  - F8 used mutant M2 instead: both `| (*_cuts %in% confs)` clauses removed.

| Fix | Pre-fix outcome (verifier, patched-file lines) | Code and probe basis |
|---|---|---|
| F1 | `:82`, `:86` (labels) pass; `:47`, `:62`, `:95`, `:98`, `:104`, `:105` fail (3.37 s) | `cut(..., labels = ..., include.lowest = TRUE, right = FALSE)`, `R/summarize_bootstrap_subgroups.R:417-423`, `:445-451` (start) |
| F5 | `:106` passes; `:107`, `:108` fail (5.91 s) | `probes/f5_probe.txt`: labels `age <= 55` ... `bmi <= 30`, 0 warnings; `class(result) <- "hte_test"` at `R/test_hte_crump.R:268` (start) |
| F6 | `:67` passes x3; `:76` fails x3; `:77` passes x3; test 2 passes in full x3; test 4: survival passes both checks, binary fails both (length 0), continuous and count fail the match only (5.75 s) | `polish/f6_test4_msgs_pre_F6.txt`; `polish/f6_test4_probe_pre_F6.txt` (4.56 s): the binary fit stops in `get_FSdata` ("LASSO selected (1)"), continuous prints "Cox-LASSO selected: 4 of 5 candidate factors", count prints "... 2 of 5 ..." |
| F8 | M2 passes all 13 expectations (3.32 s) | Both merges end in `unique(c(confs, *_cuts_keep))`; `protect` only under `collapse_cuts = TRUE` (`R/get_fsdata.R`, start) |
| F9 | `:73`, `:95` pass; `:75`, `:86` fail (7.46 s) | `probes/f9_probe_C6.txt`: `df_super` has 30 columns, 9 match `^z_\|^ecog\|^strat`; one "Confounders:" line; `Hhat` = `sg_found` = `!{z_male}` on both trees |
| F10 | `:113-116`, `:146-149` pass; `:120`, `:124`, `:151`, `:153` fail (4.68 s) | Pre-fix: both tree fits give the same 497-member subgroup, all in z == 0; both frontier fits select `z <= 0` (500). Fixed: tree TRUE `z > 0` (500), frontier TRUE `z > 0` (500), frontier FALSE `z <= 0` (500). Bounds `nS < n_min \|\| nS > n - 1L` (`R/grf_subgroup_labels.R:269`), `effect >= dmin` (`:361`), same lines at both SHAs |
| F12 | not a pre-fix claim | `probes/f12_probe.txt`: MR complement limits vs hand 1.96 limits at most 2.69e-05; (1.96 - qnorm(0.975)) x 0.1236 = 4.45e-06; forest est/se relative error 1.1e-16 to 1.7e-15; 2\|est\| = 0.1916, 0.4574, 2.1882 |
| F13 | `:84-87` pass; `:72-73` (SE) fail (2.77 s) | Tier 2 `converged = TRUE`, `measure = "RD"`, `"logistic_margins"` when `n_adj == 0`; tier 3 `FALSE`, `"raw_means"` |
| F14 | not a pre-fix claim | `probes/f14_probe.txt`: recorded 0.05 (RD) and 0.01 (IRD) on both routes; validator warning first, then the remap warning |

- **The verifier's three minor points (`verify/VERIFY.md` §5) were all fixed in the polish round** (rows 6, 11, 12):
  - The polish changed 2 files: +15/-10 against the verified patch, comment lines only.
  - The polish checker (`polish_result.json`, check) gave CONFIRMED with no problems. It rebuilt the delta (`verify_polish/polish_delta_recomputed.diff`; the hunk bodies equal `evidence/C6/polish.diff`) and re-ran the probes on pre-F6 and pre-F10 copies.
- **F14 at head:**
  - Commit `2984dbee` (C2, before C6) moved the remap value into a helper. Start: `effect_threshold <- switch(effect_measure, RD = 0.05, IRD = 0.01)` at `R/forestsearch_main.R:2348-2349`, inside the branch at `:2334`. Head: `.fs_remap_default <- switch(effect_measure, RD = 0.05, IRD = 0.01)` at `:2376`, then `effect_threshold <- .fs_remap_default` at `:2394`, inside the same branch at `:2371`.
  - The branch condition and the values are unchanged, and the branch still warns. Only the warning's text changed (C2). So the comment still describes the code.
  - `?forestsearch` at head still states the "Resolved defaults, per estimand" table "when neither threshold is supplied" (`man/forestsearch.Rd:1213-1215`).
- **F13 at head:** B1 (`8e87cf30`) changed the adjusted SE only. Tier 2 still returns `converged = TRUE`, `measure = "RD"` and `"logistic_margins"` when `n_adj == 0` (`R/glm_effect_estimators.R:583-587` at head). Tier 3 still returns `FALSE` and `"raw_means"` (`:602-606`).
- **F12 at head:** `analyze_subgroup_glm()` still has `z <- 1.96` (`R/summary_utility_functions.R:223`, both SHAs).

### Gates and runs

| Check | Evidence | Result |
|---|---|---|
| Comment-only gate, per file, on the commit | `main/C6/comment_gate.txt` ("COMMENT GATE: PASS", 9 files). Re-run here on `git diff -U0 5073fd67 ffc55486`: 0 non-comment lines in each of the 9 files | PASS |
| Code unchanged | Implementer `evidence/C6/gate.txt`: with comment and blank lines stripped, each file is identical to start (9 of 9). Verifier: `deparse(parse(keep.source = FALSE))` identical (9 of 9). Polish: `polish/parse_check.txt` identical to start and to the verified copy (2 of 2) | PASS |
| Hygiene of added lines | On the commit: 0 blank, 0 non-ASCII, 0 over 80 characters | PASS |
| Scope | `main/C6/diffstat.txt` = `git show --stat ffc55486`: 9 files under `tests/testthat/`, +111/-8 | PASS |
| Values | No expected value, seed or tolerance changed (follows from the gate) | PASS |

Expectation counts. The part 4 column is from its report §6. Main "after" is `main/C6/tests_after.txt`, counting the dots per file: one process for all nine files, real 23.079 s, run before the commit (file written 17:07:55; commit 17:08:07). Blocks come from the scratch runs, because `tests_after.txt` does not show them.

| File | Part 4 blocks / exp | Main after: exp | Implementer before → after: blocks / exp (wall s) | Verifier start → patched (wall s) |
|---|---|---|---|---|
| `test-bootstrap-summary-bins.R` | 3 / 8 | 8 | 3 / 8 → 3 / 8 (3.21 → 3.36) | 3 / 8 → 3 / 8 (3.15 → 3.21) |
| `test-hte-crump-confounders.R` | 5 / 35 | 35 | 5 / 35 → 5 / 35 (5.67 → 5.62) | same (5.50 → 5.48) |
| `test-lasso-intercept.R` | 4 / 62 | 62 | 4 / 62 → 4 / 62 (6.15 → 6.27) | same (6.00 → 5.69) |
| `test-duplicate-cut-names.R` | 7 / 13 | 13 | 7 / 13 → 7 / 13 (3.37 → 3.36) | same (3.21 → 3.24) |
| `test-mrct-default-confounders.R` | 1 / 8 | 8 | 1 / 8 → 1 / 8 (7.49 → 7.54) | same (7.73 → 7.64) |
| `test-grf-count-orientation.R` | 2 / 14 | 14 | 2 / 14 → 2 / 14 (4.44 → 4.60) | same (4.58 → 4.43) |
| `test-display-orientation.R` | 6 / 93 | 93 | 6 / 93 → 6 / 93 (7.32 → 7.23) | same (7.16 → 7.12) |
| `test-rd-standard-error.R` | 2 / 14 | 14 | 2 / 14 → 2 / 14 (3.04 → 3.13) | same (3.02 → 3.03) |
| `test-dina-gaussian-floor.R` | 6 / 55 | 55 | 6 / 55 → 6 / 55 (14.85 → 15.50) | same (15.18 → 15.20) |
| **Total** | **36 / 302** | **302** | 36 / 302 → 36 / 302 | 36 / 302 → 36 / 302 |

- Every run had 0 failures, 0 errors, 0 skips and 0 warnings. `tests_after.txt` has no Failed, Skipped or Warnings section.
- The implementer's "after" ran on the pre-polish copy. The polish re-ran the two files it changed:
  - `test-lasso-intercept.R`: 4 / 62, 5.93 s (polish checker 5.71 s);
  - `test-grf-count-orientation.R`: 2 / 14, 4.47 s (polish checker 4.35 s).

### Departures and notes

- **Three optional F6 comments added** (rows 4-6). They come from part 4's F6 Verification paragraph, not from item 17.
- **Three existing comments rewritten instead of only extended:** F6 `:69`, F12 `:27-30` and F14 `:159-161` (table above). In F12 and F14 the old text was wrong or incomplete.
- **The F14 comment does not quote the warning text.** C2 changes that text, and C4 changes the `?forestsearch` sentence on supplied thresholds. The comment names the branch and the `switch()` instead.
- **No line numbers in the comments.** Functions and code expressions are named instead, because line numbers move.
- **Process:** the implementer's first application put the F10 frontier comment into the tree block, because the anchor occurred twice. The implementer caught it before writing the patch, restored the files and re-applied the comments with anchors unique to each block (`evidence/C6/NOTES.md` §8).
- **Commit message:** `test: derivation comments on part 4's test expectations (part 5, C6)`. The implementer had proposed `docs: derivation comments on the part-4 test expectations (C6)`.
- **NEWS:** no part-5 bullet names the test comments.
  - The `fb80dd8f` commit message says it adds "one Documentation entry for C4-C7".
  - That entry ("Documentation pass (part 5 of the CRAN preparation; no behaviour change)") has sub-bullets for help pages, `\dontrun{}` examples and the two articles, and none for C6.
  - The implementer's NOTES §10 proposed that the documentation-pass entry cover C6, with an optional line if one is wanted.

### Open items C6 leaves

- **F12:** the 1e-4 tolerance on the MR limits (`test-display-orientation.R:166` at head) is twice the 5e-5 rounding bound. It was kept because C6 changes no tolerance, and the new header comment says so.
- **F13:** the first block's `measure` and `n0`/`n1` checks (`test-rd-standard-error.R:61-62`, same lines at start and head) have no comment of their own. Its comment at `:56-58` covers `converged` and `method_used`. This was not requested.

## C7 — two article sources (detail)

Line numbers: start `b06de4e8` → head `fb80dd8f`. Both `.qmd` files are identical at `835fac1b` and `fb80dd8f` (`git diff --stat 835fac1b fb80dd8f` lists only `DESCRIPTION` and `NEWS.md`). Scratch root `S=/tmp/cran_prep_displays_se_docs_2026-11-12`.

### C7.1 Provenance

| Step | Record | Result |
|---|---|---|
| Implement | `$S/evidence/C7/NOTES.md` §1–10; copy `$S/wt/C7/` | READY; patch +42/-17 (kept as `$S/patches/C7.verified.patch`) |
| Verify, round 1 | `$S/evidence/C7/verify/VERIFY.md` | CONFIRMED; 0 blocking, 6 minor |
| Polish (minor 2 and 3) | NOTES "Polish before commit"; copy `$S/wt/C7_polish/` | +10/-2 on the verified tree; patch +50/-17 |
| Polish check | `$S/polish_result.json` (entry `C7`, `check`); probes in `$S/evidence/C7/verify_polish/` | PROBLEMS: 2 |
| Main-session correction | NOTES "Main-session correction after the polish check"; `$S/evidence/C7/polish.diff` | sentence reworded; delta from the verified tree now +9/-2; final patch +49/-17 |
| Commit | `835fac1b` (parent `ffc55486`); `$S/main/C7/{commit,diffstat,chunk_gate,patch_equal}.txt` | 2 files, +49/-17; `C7 GATE: PASS` |
| NEWS | `fb80dd8f:NEWS.md:560-562` (Documentation pass entry) | one bullet for C7 |

- Commit subject: `docs: two article sources name setup_gbsg_dgm() and say what the flex block is (part 5, C7)`. The implementer had proposed a different subject (NOTES §10); the main session used its own (`$S/msg/C7.txt`).
- The commit equals `$S/patches/C7.patch`. Checked for this section: `git diff --binary ffc55486 835fac1b` and the patch differ only in the `index` lines (7- vs 8-character hashes). `$S/main/C7/patch_equal.txt` is empty (0 bytes).
- Per file: `biomarker_effects.qmd` +3/-1, `paper_simulations.qmd` +46/-16 (`git diff --numstat ffc55486 835fac1b`).

### C7.2 Sites: old and new text in full

Site labels used here, with the implementer's numbers (NOTES §2):

| Label | NOTES site | File | Start lines | Head lines | Inside the task's named ranges |
|---|---|---|---|---|---|
| B | 1 | `biomarker_effects.qmd` | 66 (paragraph 63–68) | 66, 69–70 (paragraph 63–70) | yes ("near 66") |
| P1 | 6 | `paper_simulations.qmd` | 366–367 | 366–369 | no (added) |
| P2 | 2 | `paper_simulations.qmd` | 388, 391–393 | 390, 393–396 | yes (388–445) |
| P3 | 3 | `paper_simulations.qmd` | 442–445 | 445–475 | yes (388–445) |
| P4 | 5 | `paper_simulations.qmd` | 454 | 484 | no (added; between the two ranges) |
| P5 | 4 | `paper_simulations.qmd` | 459–463 | 489–493 | yes (459–463) |

#### B — `biomarker_effects.qmd`

Old, `b06de4e8:vignettes/articles/biomarker_effects.qmd:63-68` (line 66 changes):

~~~~text
This DGM family is the one behind the published simulation studies.
León and Anderson (2026a, Section 7) evaluate the post-selection correction
on an accelerated-failure-time super-population built from the GBSG data by
`create_gbsg_dgm()` and `generate_aft_dgm_flex()` (their Supplementary
Section S9), and 2026b (Section 5.1) builds its survival design from the
same trial's covariate structure, censoring and event rates.
~~~~

New, `fb80dd8f:vignettes/articles/biomarker_effects.qmd:63-70` (lines 66, 69, 70 changed or added):

~~~~text
This DGM family is the one behind the published simulation studies.
León and Anderson (2026a, Section 7) evaluate the post-selection correction
on an accelerated-failure-time super-population built from the GBSG data by
`setup_gbsg_dgm()` and `generate_aft_dgm_flex()` (their Supplementary
Section S9), and 2026b (Section 5.1) builds its survival design from the
same trial's covariate structure, censoring and event rates.
`setup_gbsg_dgm()` uses the same internal builder as the removed
`create_gbsg_dgm()`.
~~~~

#### P1 — "What `setup_gbsg_dgm()` does"

Old, start `paper_simulations.qmd:364-369` (366–367 change):

~~~~text
`setup_gbsg_dgm()` is a **GBSG-specific convenience wrapper**. It encodes
the particular data preparation choices made in León et al. (2024) — the
subgroup definition (ER ≤ q25, premenopausal), the binary discretisation of
continuous covariates, the Weibull outcome and censoring models — so that a
calibrated DGM for the GBSG simulation setting can be created with a single
call:
~~~~

New, head `paper_simulations.qmd:364-371` (366–369 changed):

~~~~text
`setup_gbsg_dgm()` is a **GBSG-specific convenience wrapper**. It encodes
the particular data preparation choices made in León et al. (2024) — the
subgroup definition (ER ≤ q25, premenopausal), the discretisation of
continuous covariates (age, progesterone receptor and nodes in quartile
groups; tumour size in three groups cut at 20 and 50 mm), the Weibull
outcome and censoring models — so that a
calibrated DGM for the GBSG simulation setting can be created with a single
call:
~~~~

#### P2 — heading and lead-in to the `generate-aft-equivalent` chunk

Old, start `paper_simulations.qmd:388-393`:

~~~~text
## Replicating `setup_gbsg_dgm()` with `generate_aft_dgm_flex()`

`generate_aft_dgm_flex()` is the general DGM builder. It requires the analyst
to supply the dataset and make every modelling choice explicitly. The GBSG
DGM in `setup_gbsg_dgm()` is equivalent to the following call — it simply
automates this preparation:
~~~~

New, head `paper_simulations.qmd:390-396`:

~~~~text
## A comparable construction with `generate_aft_dgm_flex()`

`generate_aft_dgm_flex()` is the general DGM builder. It requires the analyst
to supply the dataset and make every modelling choice explicitly. The call
below is a comparable construction of the GBSG mechanism with the general
generator. It is not numerically equivalent to `setup_gbsg_dgm()`; the
differences are listed after the code.
~~~~

#### P3 — paragraph after the chunk

Old, start `paper_simulations.qmd:442-445`:

~~~~text
`setup_gbsg_dgm()` is preferred for the GBSG simulation setting because it
guarantees exact replication of the León et al. (2024) covariate
construction and subgroup definition without requiring the analyst to
reproduce all preparation steps by hand.
~~~~

New, head `paper_simulations.qmd:445-475`:

~~~~text
The block and `setup_gbsg_dgm()` share the source data, the outcome time in
months, the subgroup $H$ (ER at or below its 25th percentile and
premenopausal), `k_treat`, `k_inter`, `n_super` and `seed`. With the same
seed they resample the same source rows into the super-population, so the
random draws are not what differs. The super-population's treatment column
does differ: `setup_gbsg_dgm()` assigns treatment to half of the resampled
rows (share 0.5), while the block, with the default
`draw_treatment = FALSE` of `generate_aft_dgm_flex()`, keeps the observed
`hormon` (share about 0.36). In a simulated trial this column matters only
when treatment is not redrawn; `simulate_from_dgm()` redraws it in every
trial by default (`draw_treatment = TRUE`). The models differ:

- **Covariates and cut points.** `setup_gbsg_dgm()` fits the Weibull outcome
  model on treatment, the ER and premenopausal indicators, and the quartile
  group (1 to 4, entered as a number) of age, progesterone receptor and
  nodes, plus the treatment-by-$H$ interaction; tumour size and grade are
  analysis factors only. The block uses median splits (0/1) of age,
  progesterone receptor and nodes, and adds tumour size (median split) and
  grade 3 to the outcome model.
- **Interaction.** Both builders multiply the fitted interaction coefficient
  by `k_inter`, but the coefficient comes from a different outcome model, so
  `k_inter = 2.5` does not give the same interaction: the true hazard ratio
  in $H$ is about 4.4 from `setup_gbsg_dgm()` and 4.1 from the block.
- **Censoring.** `setup_gbsg_dgm()` fits a Weibull censoring model to the
  resampled super-population, in which treatment is the arm assigned for
  simulation, and multiplies its covariate term by the ratio of the
  event-model scale to the censoring-model scale (`?setup_gbsg_dgm`, section
  Superseded). `generate_aft_dgm_flex()`, with its default
  `select_censoring = TRUE`, fits the censoring model to the source data and
  chooses by AIC among Weibull and log-normal models with and without
  covariates; for this call it chooses an intercept-only Weibull model.
~~~~

#### P4 — "When to use each function" table, one row

Old, start `paper_simulations.qmd:454`:

~~~~text
| GBSG-based simulation with custom parameters | `setup_gbsg_dgm()` |
~~~~

New, head `paper_simulations.qmd:484`:

~~~~text
| New GBSG-based simulation with custom parameters | `generate_aft_dgm_flex()` |
~~~~

The row above it (start 453, head 483), "Replicating León et al. (2024) simulations → `setup_gbsg_dgm()`", is unchanged.

#### P5 — paragraph after the table

Old, start `paper_simulations.qmd:459-463`:

~~~~text
`setup_gbsg_dgm()` is not deprecated in the sense of being removed or
discouraged for GBSG work — it remains the correct entry point for the
simulation studies in this package. The general function `generate_aft_dgm_flex()`
is the right choice when the simulation is based on a **different dataset or
subgroup structure**.
~~~~

New, head `paper_simulations.qmd:489-493`:

~~~~text
`setup_gbsg_dgm()` is used in this article to reproduce the published
simulations, and the package keeps it for that purpose. It is superseded
for new work and will be deprecated in a future release: new simulations,
on the GBSG data or another dataset, should use `generate_aft_dgm_flex()`
with `simulate_from_dgm()`.
~~~~

- Wording removed from the prose (`git show <sha>:...paper_simulations.qmd | grep -n -i -E "equivalent|correct entry point|exact replication|binary discretisation"`): at start, prose hits at 366, 392, 443, 460; at head, the only prose hit is the new "not numerically equivalent" (395). Hits left at head: the chunk label (398), an unrelated table header (2397) and an unrelated "exact replication" at 2845 (start 2815), both outside C7's sites.

### C7.3 Code facts behind each new claim

Code files read are unchanged from start to head except `R/setup_gbsg_dgm.R` (C4, `ba7c4ea2`, replaced line 74 with 38 lines, so its lines from 75 on move by +37). `git diff --stat b06de4e8 fb80dd8f` is empty for `R/sim_aft_gbsg.R`, `R/simulate_from_dgm.R`, `R/generate_aft_dgm_main.R`, `R/generate_aft_dgm_helpers.R` and `R/generate_aft_dgm_helpers_outcome.R`; their lines are the same at both SHAs.

Measurements were made on the start tree (no code changed): implementer `$S/evidence/C7/compare_builders.R` → `compare_builders.txt` (wall 2.94 s, `/usr/bin/time`); verifier's own `$S/evidence/C7/verify/v_compare.R` → `v_compare.txt` (wall 3.33 s, `/usr/bin/time`; 3.15 s `proc.time()`); polish check `$S/evidence/C7/verify_polish/probe.R` → `probe.txt` (wall 3.27 s, `/usr/bin/time`; 3.1 s `proc.time()`) and `probe_nnull.R` → `probe_nnull.txt` (no wall recorded). Each runs the article's `generate-aft-equivalent` chunk (start 396–439, `#|` and `library(forestsearch)` lines dropped) and the article's `setup_gbsg_dgm(model = "alt", k_inter = 2.5, k_treat = 1.0, seed = 8316951)`.

| Claim (head line) | Code, start `b06de4e8` | Code, head `fb80dd8f` | Measured |
|---|---|---|---|
| B: AFT super-population "built from the GBSG data by `setup_gbsg_dgm()`" (66) | `R/setup_gbsg_dgm.R:154` calls `.create_gbsg_dgm_()`; it reads `survival::gbsg` (`R/sim_aft_gbsg.R:119`), fits a Weibull `survreg` (`:224-228`), resamples `n_super` rows (`:327-331`) | `R/setup_gbsg_dgm.R:191`; `sim_aft_gbsg.R` same lines | — |
| B: "uses the same internal builder as the removed `create_gbsg_dgm()`" (69–70) | removed body `add6cabc^:R/sim_aft_gbsg.R:199-231` passes every argument to `.create_gbsg_dgm_()` (`:226-230`); `R/sim_aft_gbsg.R` changed after `add6cabc` only in `dee93483`, comment lines only (comment gate on `git diff -U0 add6cabc b06de4e8` prints nothing); `NEWS.md:410-412` ("same internal builder") | `NEWS.md:472-474` | — |
| B: the cited paper names `create_gbsg_dgm()` (reason for the added sentence) | fs-post-selection `HEAD:manuscript/postselection_subgroup_inference_supplement.qmd:1182` (read only) | same | — |
| P1: age, progesterone receptor, nodes in quartile groups (367–368) | `R/sim_aft_gbsg.R:130,132,133` `cut_numeric()`; `cut_numeric()` `:51-58` (probs 0.25/0.5/0.75 → codes 1–4) | same | z2, z4, z5 take 1 2 3 4 (`v_compare.txt`) |
| P1: tumour size in three groups cut at 20 and 50 mm (368) | `:147` `v6 <- as.factor(cut_size(size))`; `cut_size()` `:69-75`, `breaks = c(20, 50)`: 1, 2 for `> 20`, 3 for `> 50` | same | 180 / 453 / 53 of 686; `cut_size(20)` = 1, `(50)` = 2, `(51)` = 3 (`probe.txt`); `?survival::gbsg` "tumor size, mm" (NOTES polish) |
| P2: "comparable construction … not numerically equivalent" (393–396) | rows below | — | rows below |
| P3: shared source data, months, $H$, `k_treat`, `k_inter`, `n_super`, `seed` (445–447) | `DAYS_PER_MONTH <- 30.4375` (`sim_aft_gbsg.R:32`, used `:122`); block `dfa$y <- dfa$rfstime / 30.4375` (article start 406, head 409); $H$ `:139` (`z1 == 1 & z3 == 1`), `z1` = ER ≤ 25th percentile (`:128-129`), `z3` = `meno == 0` (`:131`); default `n_super` 5000 (`:35`) | same; block at head 409 | $H$ 84 of 686 source rows; 634 of 5000 super-population rows in both; flag sequences `identical()` TRUE (`v_compare.txt`) |
| P3: same seed resamples the same source rows (447–449) | setup `set.seed(seed)` then `sample()` (`sim_aft_gbsg.R:327-330`); flex `sample()` (`R/generate_aft_dgm_helpers.R:858`) | same | `pid` sequence `identical()` TRUE (`compare_builders.txt`, `check_resample.txt` wall 3.03 s, `v_compare.txt`); first `rexp()` at the same RNG state (`v_compare.txt`) |
| P3: setup assigns treatment to half the rows, share 0.5 (450–451) | `n_treat <- round(n_super / 2)` and first `n_treat` rows `treat <- 1L`, rest `0L` (`sim_aft_gbsg.R:435-447`); stored as `df_super_rand` (`:529`); `R/setup_gbsg_dgm.R:195,209` copies it to `df_super` (treat kept) | `R/setup_gbsg_dgm.R:232,246` | 2500 of 5000; first 2500 are 1 (`probe.txt`) |
| P3: block keeps observed `hormon` with default `draw_treatment = FALSE`, share about 0.36 (451–453) | block `dfa$treat <- dfa$hormon` (article start 407, head 410); `generate_aft_dgm_flex(draw_treatment = FALSE)` (`R/generate_aft_dgm_main.R:229`, roxygen `:25-26`), passed on at `:338`; halves only `if (draw_treatment)` (`R/generate_aft_dgm_helpers.R:867-873`) | same | 0.3598 in the super-population; 246/686 = 0.3586 in the source (`probe.txt`) |
| P3: column matters only when treatment is not redrawn; redrawn by default, `draw_treatment = TRUE` (453–455) | `simulate_from_dgm(draw_treatment = TRUE)` (`R/simulate_from_dgm.R:179`); trial path sets `treat_sim <- treat` (`:266`) and overwrites it only `if (draw_treatment)` (`:271-291`); `baseline = "fixed"` takes rows from `dgm$df_source` (`:245`); `n = NULL` with `baseline = "resample"` keeps `df_super$treat` (`:221-234`) | same | n = 700: redrawn, share 0.5 for both; `draw_treatment = FALSE`: equals `df_super$treat[idx]` for both (0.5257 setup, 0.3529 block) (`probe.txt`) |
| P3: setup outcome model = treat, ER and premenopausal indicators, quartile groups of age/pgr/nodes, treat-by-$H$; size and grade analysis factors only (457–461) | `covs_main <- c("treat","z1","z2","z3","z4","z5")`, `+ "zh"` for `alt` (`sim_aft_gbsg.R:181-184`); `zh` `:136`; `v6`, `v7` built `:147-149` | same | gamma names `treat z1 z2 z3 z4 z5 zh`; z2 integer 1–4 |
| P3: block uses median splits and adds size (median split) and grade 3 (461–463) | `fit_aft_model()` takes every `^z_` column (`R/generate_aft_dgm_helpers_outcome.R:8`, formula `:46-50`, fit `:52-54`); block article start 413, 417–418, `factor_vars` 429 (head 416, 420–421, 432) | same | gamma names `treat z_z1 … z_z5 z_v6 z_v7 treat_harm`; `z_z2 == (age <= median)`, `z_v6 == (size <= median)`, `z_v7 == (grade == 3)` TRUE |
| P3: both multiply the fitted interaction by `k_inter` (464–465) | setup `gamma["zh"] <- k_inter * gamma["zh"]` (`sim_aft_gbsg.R:269`); flex `b0 <- -gamma / tau`, `b0["treat_harm"] <- k_inter * b0["treat_harm"]`, `gamma <- -b0 * tau` (`..._outcome.R:62,66-68,76`) | same | setup −0.673744 → −1.684360; flex −0.640353 → −1.600883; ratio 2.5 in both (`v_compare.txt`) |
| P3: true HR in $H$ about 4.4 (setup) and 4.1 (block) (466–467) | — | — | 4.3808 vs 4.0743; $H^c$ 0.6612 vs 0.6483; overall 0.7411 vs 0.7242; AHR($H$) 6.2045 vs 5.4087 (`compare_builders.txt`, reproduced in `v_compare.txt`) |
| P3: setup censoring = Weibull on the super-population with assigned treatment, covariate term × event scale / censoring scale (468–472) | `survreg(Surv(y, 1 - event) ~ z_cens, data = df_big)` with `z_cens <- df_big[, covs_main]` (`sim_aft_gbsg.R:474-481`); `b_cens <- -b0_cens * sigma` with `b0_cens <- -gamma_cens / sigma_cens` (`:486,490`); `?setup_gbsg_dgm` Superseded `R/setup_gbsg_dgm.R:88-93`, `man/setup_gbsg_dgm.Rd:113` | `R/setup_gbsg_dgm.R:125-130`; `man/setup_gbsg_dgm.Rd:149` | mu 4.0492, scale 0.4363; `linC1.conf == Z %*% (gamma_cens * sigma / sigma_cens)` TRUE; refit on `df_super_rand` reproduces the fit (`v_compare.txt`) |
| P3: flex censoring, default `select_censoring = TRUE`, on the source data, AIC among Weibull / log-normal with and without covariates; intercept-only Weibull here (472–475) | default `R/generate_aft_dgm_main.R:234`; four `survreg` fits on `df_work` (`R/generate_aft_dgm_helpers.R:312-317`), compared and best taken (`:320-326`) | same | chose `Weibull0`: mu 4.0312, scale 0.4246, 0 covariate coefficients; equals `survreg(Surv(y, 1 - event) ~ 1, dist = "weibull")` on the 686 source rows (`compare_builders.txt`, `v_compare.txt`) |
| P4, P5: kept to reproduce published simulations; superseded for new work; deprecated in a future release; new work → `generate_aft_dgm_flex()` with `simulate_from_dgm()` (484, 489–493) | `?setup_gbsg_dgm` Superseded: "kept so that published simulations can be reproduced", "New work should use …" (`R/setup_gbsg_dgm.R:80-83`), "will be deprecated in a future release" (`:99-100`) | `R/setup_gbsg_dgm.R:117-120,136-137` | — |
| P5: "used in this article to reproduce the published simulations" (489) | article DGMs built with `setup_gbsg_dgm()` at start 274, 507, 614, 645, 699, 1756 | head 274, 537, 644, 675, 729, 1786 | — |

- The prose's "about 4.4 … 4.1" are 4.3808 and 4.0743 rounded to one decimal. "About 0.36" covers 0.3598 (super-population) and 0.3586 (source).

### C7.4 Departures from the task's wording

- **B, added sentence.** The task asked only to name `setup_gbsg_dgm()`. The added sentence (head 69–70) maps the paper's `create_gbsg_dgm()` to `setup_gbsg_dgm()`. Basis in C7.3.
- **P1, outside the named ranges.** Start 366–367 said "binary discretisation"; the code uses quartile groups (age, pgr, nodes) and three size groups. The P3 bullet states the quartile groups, so the two passages would have contradicted each other.
- **P4, outside the named ranges.** Start 454 recommended `setup_gbsg_dgm()` for "GBSG-based simulation with custom parameters", which is new work and contradicts "superseded for new work".
- **P3, treatment column.** Added in the polish (verifier minor 2). The task does not name it; it is a difference the code shows.
- **`k_inter` and seeds are not differences.** NOTES §8 records that the same multiplication applies in both builders (only the fitted coefficient differs) and that the same seed draws the same rows; the prose says both.
- Verifier: the heading, the "comparable construction … not numerically equivalent" sentence, the listed differences and the "used in this article … superseded for new work" paragraph are all present; it judged the B, P1, P4 additions correct (VERIFY.md, Check 2).

### C7.5 Polish and main-session corrections

**Round 1 → polish (verifier minor 2 and 3).** Delta on the verified tree, +10/-2 (NOTES "Polish before commit"):

- Minor 3 → P1: "groups)," became "groups; tumour size in three groups cut at 20 and 50 mm),". Kept in the commit (head 368).
- Minor 2 → P3: the polisher's text for the treatment column (verified tree 448 → polished 448–456):

~~~~text
random draws are not what differs. The super-population's treatment column
does differ: `setup_gbsg_dgm()` assigns treatment to half of the resampled
rows (share 0.5), while the block, with the default
`draw_treatment = FALSE` of `generate_aft_dgm_flex()`, keeps the observed
`hormon` (share about 0.36). `simulate_from_dgm()` uses this column only
when it keeps the super-population's arms, with `draw_treatment = FALSE` or
with `n = NULL` (the whole super-population as it is); otherwise it redraws
treatment in every simulated trial. The models differ:
~~~~

- Not addressed by the polish: verifier minor 1 (chunk label and in-chunk comment; inside a chunk), 4 (4.4 / 4.1 tied to current builders), 5 (`biomarker_effects.html` not re-rendered), 6 (`docs:` prefix).

**Polish check (`polish_result.json`, `C7.check`): PROBLEMS, two.**

1. `baseline = "fixed"`: `simulate_from_dgm()` takes `dgm$df_source` (`R/simulate_from_dgm.R:245,266`) and, with `draw_treatment = FALSE`, keeps the source arms; so "otherwise it redraws" was false there (probe C4: block DGM, 686 rows, `treat_sim == df_source$treat`, share 0.3586).
2. `n = NULL` works only for the block's DGM: on the `setup_gbsg_dgm()` object, `simulate_from_dgm(dgm)` stops with "replacement has 0 rows, data has 5000" (`probe_nnull.txt`), because its `df_super` has no `lin_pred_obs` (`:310`). The polisher's text named `n = NULL` as a working path for both builders.

**Main-session correction (17:04; `polish.diff`, now +9/-2 on the verified tree).** The last three lines of the polisher's text were replaced. Committed text, head 449–455:

~~~~text
random draws are not what differs. The super-population's treatment column
does differ: `setup_gbsg_dgm()` assigns treatment to half of the resampled
rows (share 0.5), while the block, with the default
`draw_treatment = FALSE` of `generate_aft_dgm_flex()`, keeps the observed
`hormon` (share about 0.36). In a simulated trial this column matters only
when treatment is not redrawn; `simulate_from_dgm()` redraws it in every
trial by default (`draw_treatment = TRUE`). The models differ:
~~~~

- Problem 2: `n = NULL` is no longer named.
- Problem 1: "matters only when treatment is not redrawn" is a necessary condition. It still holds under `baseline = "fixed"`, where the super-population column is not used at all. With `baseline = "fixed"` and the default `draw_treatment = TRUE`, treatment is redrawn (probe C5, block: share 0.5).
- No verifier ran on the corrected sentence. The `verify_polish/` files date 16:53–16:55 and test the polisher's wording. The checks above are read from the code at `b06de4e8` (unchanged at head) and from the existing probes.
- The `n = NULL` path is still there. With `n = NULL` (the default of `n`) and `baseline = "resample"`, `simulate_from_dgm()` uses `df_super` as it is and sets `treat_sim <- treat` (`R/simulate_from_dgm.R:221-234`), whatever `draw_treatment` is. So "redraws it in every trial by default" holds for calls that give `n`, as the article's calls do (head 509 and 757, `n = 700`). See open items.
- The final patch, start → committed, is +49/-17 (`patches/C7.patch` = commit `835fac1b`).

### C7.6 Gate

Task gate: only the two files change under `vignettes/`; no changed line is inside a code chunk; no `.html` is edited; nothing is rendered.

- **Main session, on the commit** (`$S/main/C7/chunk_gate.txt`, no wall recorded):

~~~~text
changed files: ['vignettes/articles/biomarker_effects.qmd', 'vignettes/articles/paper_simulations.qmd']
vignettes/articles/biomarker_effects.qmd: 4 changed lines (old+new), inside a chunk or a fence: 0
vignettes/articles/paper_simulations.qmd: 62 changed lines (old+new), inside a chunk or a fence: 0
html edited: []
C7 GATE: PASS
~~~~

- **Re-checked from git for this section** (`ffc55486` → `835fac1b`; fence lines and every line between an opening and closing fence count as chunk lines):
  - `biomarker_effects.qmd`: −1 (start 66), +3 (head 66, 69, 70); 0 in a chunk; 5 chunks before and after, identical.
  - `paper_simulations.qmd`: −16 (start 366, 367, 388, 391–393, 442–445, 454, 459–463), +46 (head 366–369, 390, 393–396, 445–475, 484, 489–493); 0 in a chunk; 64 chunks before and after, identical; no unclosed fence.
  - Nearest fences, start → head: B between 34 and 171 → 34 and 173; P1 between 340 and 371 → 340 and 373 (`setup-gbsg-dgm-usage`); P2 between 381 and 395 → 383 and 398 (`generate-aft-equivalent`); P3–P5 between 440 and 473 → 443 and 503 (`simulate-from-either`).
  - `git diff --check ffc55486 835fac1b` clean. Non-ASCII in added lines: `≤` (head 366) and `—` (head 369), both carried over from the start sentence.
- **Scope for the whole part:** `git diff --name-only b06de4e8 fb80dd8f` lists under `vignettes/` only the two C7 files, nothing under `quarto/`, and no `.html`.
- **Earlier gate runs:** implementer `gate_c7.py` → `gate.txt`, GATE: PASS, 0.05 s; verifier's own `verify/fence_check.py` → `fence_check.txt`, FENCE CHECK: PASS, 0.04 s; polish `polish_gate.txt`, GATE: PASS, 0.05 s (before the main-session correction). Verifier copies by `mkcopy.sh`: 0.38 s and 0.39 s.
- **R/tests comment gate:** no file under `R/`, `tests/` or `man/` changes. `devtools::document()` on the polish-check copy changed nothing (`verify_polish/document.txt`, wall 12.09 s).
- **Build:** `vignettes/articles` is in `.Rbuildignore` (`^vignettes/articles$`, line 20 at start and head), so neither file reaches the package build or check.
- **Articles statement:** every code chunk in both files is byte-identical, so no evaluated call changes and no article output can change from C7. Nothing was re-rendered.
- **Existing tests:** none apply; none run.

### C7.7 Open items

1. The chunk label `generate-aft-equivalent` (head 398) still says "equivalent". The in-chunk comment "Analyst-visible covariates (what FS will search over)" (head 430) is misleading: `fit_aft_model()` puts every `^z_` column, including those from `factor_vars`, in the outcome model (`R/generate_aft_dgm_helpers_outcome.R:8`). Both are inside a chunk, so the gate kept them.
2. The committed render `vignettes/articles/biomarker_effects.html` still has 1 `create_gbsg_dgm` and 0 `setup_gbsg_dgm` at head until the site is re-rendered (not edited, by rule). There is no committed `paper_simulations.html`.
3. "About 4.4 … 4.1" (head 467) and "share about 0.36" (head 453) were measured on the current builders. They go stale if either builder changes.
4. The corrected sentence (head 453–455) has no independent verifier record. Read against the code, it holds on the trial path. On the default `n = NULL` path, `simulate_from_dgm()` keeps the super-population's arms whatever `draw_treatment` is (`R/simulate_from_dgm.R:221-234`).
5. Found by the polish check, not changed (code, outside C7): `simulate_from_dgm(dgm)` with `n = NULL` on a `setup_gbsg_dgm()` object stops with "replacement has 0 rows, data has 5000". The cause is that `df_super` has no `lin_pred_obs` (`R/simulate_from_dgm.R:310`; `R/setup_gbsg_dgm.R:195-209` at start, `:232-246` at head). The unchanged prose at head 499–501 says `simulate_from_dgm()` "works identically with either".
6. Out of scope, unchanged: a stray `---` / `title:` / `---` block mid-document (`paper_simulations.qmd:358-360`, start and head).

## After the parts: every code fix's unaffected case at the head (full)

### A1 — 31 cases, 28 identical at head; all explained: True

- Output: `/tmp/cran_prep_displays_se_docs_2026-11-12/cumulative/A1/RESULT.md`
- Walls: Implementer script: start 12.00 s, head 11.58 s, head minus C3 10.06 s, compare 2.01 s. Verifier script: start 13.87 s, head 13.34 s, head minus C3 14.41 s, compare 2.21 s.
- G_cnt_TRUE (implementer; GRF count result, log transform, offset ftime): only tab_estimates$_source_notes[[1]] differs, at ndecimals 3 and 6. Start says 'H^ = -0.7798, H^c = 0.0508 (risk-difference scale)'; head says '(log(Y + 0.5) scale, Y = count, not offset-adjusted)'. Display data are identical. Traced to C3 (c3323a59), which changes the GRF CATE note's scale label for count outcomes.
- G_cnt_FALSE (implementer; GRF count result at its default adverse_outcome = FALSE, log transform, offset): same leaf and same change, 'risk-difference scale' -> 'log(Y + 0.5) scale, Y = count, not offset-adjusted' (H^ = -0.3129, H^c = 0.6072 unchanged). Display data are identical. Traced to C3.
- G_cnt_TRUE_id (verifier; GRF count result, grf_count_transform = 'identity', offset): same leaf, 'risk-difference scale' -> 'count scale, not offset-adjusted' (H^ = -1.3386, H^c = 0.1281 unchanged). Display data are identical. Traced to C3.
- Confirmation: head_fixes with C3's R/ hunk reverse-applied (git apply -R --include='R/*' patches/C3.patch) applied cleanly, hunk at line 945, no fuzz; C4 and C5 later touched the file but not the hunk's context. Start vs head minus C3 is identical in 31/31 cases: display data and whole objects with closure environments dropped; leafdiff.R shows only the 192 gt _formats closure leaves, as at gating.

### A2 — 44 cases, 44 identical at head; all explained: True

- Output: `/tmp/cran_prep_displays_se_docs_2026-11-12/cumulative/A2/RESULT.md`
- Walls: Process walls, start / head: unaffected.R (11 cases) 21.03 s / 17.81 s; unaffected_v.R (20 cases) 31.80 s / 32.43 s; unaffected_r2v.R (13 cases) 26.81 s / 27.22 s. Compare walls: 0.15 s, 0.31 s, 0.15 s. Extra positive control extract.R: 7.37 s / 7.44 s.

### A3 — 38 cases, 35 identical at head; all explained: True

- Output: `/tmp/cran_prep_displays_se_docs_2026-11-12/cumulative/A3/RESULT.md`
- Walls: Separate timeout 120 Rscript --vanilla per tree and script (/usr/bin/time %e), all exit 0. unaffected.R: start 16.92 s, head 16.68 s, head_minus_A3 17.35 s. unaff_v.R: start 48.54 s, head 43.74 s, head_minus_A3 45.35 s. break.R: start 3.68 s, head 3.94 s, head_minus_A3 3.89 s. Each compare in a third process: 0.15-0.16 s. Normalisation: none; the outputs hold compute_frontier_cis() tables, printed text, warnings and error strings, with no timing leaf and no ggplot or gt object, as in the original scripts. Scripts were copied byte-identical, with trees and outputs passed as arguments. break_cmp.R was reconstructed because the verifier saved only break_cmp.log. Cross-checks: head is identical() to the A3-only gate outputs (11/11, 20/20 whole object, 7/7), and start is identical to the earlier start runs (11/11, 20/20, 7/7). In summary, 38 cases (implementer 11, verifier unaff_v 20, verifier break 7): 35 identical at head. The 3 that differ are A3's own intended unstored-continuous sign flip, and all 38 are identical to start once A3's R/ hunk is reverse-applied to the head. No other part-5 fix changes any of these cases.
- continuous_unstored_INTENDED (implementer unaffected.R): A3 itself, its intended change. For a continuous fit with the adverse_outcome entry removed, the fallback goes from TRUE to forestsearch()'s rule (FALSE), so every estimate flips sign and the limits swap and negate. The head output is identical() to A3's own gate output (evidence/A3/unaffected_fixed.rds). On head_minus_A3 (head with A3's R/ hunk reverse-applied, which applied cleanly) it is identical to start.
- UNSTORED_cont_INTENDED (verifier unaff_v.R): A3 itself, its intended change. unaff_intended.log: on head the unstored table is identical() to the stored-FALSE table; start estimate equals -(head estimate) exactly (m=1: -1.1298341 vs 1.1298341). The head output is identical() to verify/unaff_fixed.rds. On head_minus_A3 it is identical to start.
- cont_no_measure_INTENDED (verifier break.R): A3 itself, its intended change. Estimates are -1.094, -0.7847, -0.586 on start and 1.094, 0.7847, 0.586 on head. The head output is identical() to verify/break_fixed.rds. On head_minus_A3 it is identical to start.

### B1 — 2679 cases, 2253 identical at head; all explained: True

- Output: `/tmp/cran_prep_displays_se_docs_2026-11-12/cumulative/B1/RESULT.md`
- Walls: Run walls in seconds, start/head/head-minus-B1/start-plus-B1: unaffected_impl 5.11/5.06/5.69/5.57; unaffected_verif 8.21/8.24/8.47/8.72; unaffected_verif_r2 9.02/9.72/9.34/11.20; probe_verif 2.62/2.72/3.01/2.97; probe_verif_r2 5.89/6.04/6.37/6.45. Compares start vs head: impl 0.34, verif 0.22, verif_r2 0.29, probes 0.22. Full compares: start vs head 0.45, confirmation 0.68, evidence 0.66. Timeout 120 s each; none hit it, all exited 0.
- Unaffected sets as defined by every B1 script are identical at head: impl 216/216, verifier r1 852/852, verifier r2 1044/1044 (unadjusted RD incl. F13 slices, OR, RR, IRR/IRD, identity_binomial_adj, iptw/dr_gcomp, pre-existing muffleWarning errors); OR-default, RD-unadjusted and RD+age forestsearch() runs identical; 38/50 displays identical.
- 385 estimator results whose start method_used is logistic_gcomp_adj (impl 41, verif 92, verif_r2 252) differ in $se only (identical() with se removed TRUE; 56 not-estimable NA-se cases identical in full). Traced to B1 itself (8e87cf30, .estimate_rd() adjusted logistic tier -> .fs_rd_gcomp_se()).
- 27 probe cases (verify/probe.R 14 of 21, verify_r2/probe.R 13 of 18), all logistic_gcomp_adj, differ in $se only. Traced to B1 itself.
- 12 verif_r2 displays (subgroup_glm() and SG_tab_estimates_glm() ITT cell for RD|age on the two F13 slices and RD|site+sex on syn200, the two F13 slices, gb_all) differ in CI bounds, e.g. f13_ctl0_trt80+age 0.8044 (-830.7797, 832.3885) -> 0.8044 (0.6838, 0.9251). Traced to B1 itself.
- 2 forestsearch() RD + adjust_covariates='site' runs (verif r1, verif r2) differ only in find.grps$out.found$hr.subgroups columns L(HR)/U(HR); same rows, HR and sg.harm; fs_diff output byte-identical to the verifier's evidence. Traced to B1 itself.
- Confirmation: head_fixes with B1's R/ hunks reverse-applied (clean; 18-line offset from C4a's +23/-5 roxygen lines at :53, no fuzz) is identical() in full to start for all 5 scripts; start + B1 R/ hunks is identical() in full to head for all 5; head results identical in full to B1's own evidence fixed runs (r2 probe after removing the verifier's fd_se). No other part-5 fix (A1, A2, A3, C1, C2, C3, C4, C5, C6, C7) contributes: C2's changed branches need a user-set threshold, not set in these runs.

### C1 — 12 cases, 12 identical at head; all explained: True

- Output: `/tmp/cran_prep_displays_se_docs_2026-11-12/cumulative/C1/RESULT.md`
- Walls: Implementer scripts: start run 3.30 s, head run 3.35 s, compare 4.20 s; start-vs-start control: second run 3.70 s, compare 4.03 s. Verifier scripts: start run 3.60 s, head run 3.59 s, compare 5.62 s; control: second run 3.57 s, compare 1.91 s. Head-minus-C1 tree: implementer run 3.10 s and plain compare 3.65 s; verifier run 3.58 s and plain compare 5.26 s. All measured with /usr/bin/time -f %e. Plain compares of the C1-only outputs vs head, measured with proc.time() inside the process: 3.84 s implementer, 5.20 s verifier.
- None. All 12 cases are identical at head by C1's own gate: 4 implementer cases with 30 gt tables, and 8 verifier cases with 50 gt tables. The only difference start vs head is C1's intended highlight on the >=0.95 row, the same one seen at gating: 3 green #D4EDDA + bold _styles rows and HTML lines 480-482 in each of the 10 consistency tables. The verifier's no_pcons and none_found cases have no consistency table and match exactly, with nothing removed.
- Head equals C1 alone. The gated C1-only outputs (evidence/C1/unaffected/fixed.rds and verify/unaffected/v_fixed.rds) match the head outputs in all 80 tables, with nothing removed. All 10 rendered consistency-table HTML files are byte-identical to the C1-only files from gating.
- Reverse-apply check. C1.patch's R/ hunk reverses cleanly on a copy of head_fixes, since only C1 changes R/format_subgroup_summary_tables.R. That tree matches start with nothing removed (30/30 and 50/50 tables). The verifier's v_compare.R reports FALSE on this tree only because it requires C1's 3 highlight rows to be present.
- The only other part-5 edit to a file these cases execute is C4a's roxygen (#') lines 22-35 in R/summarize_bootstrap_subgroups.R. It changes no output: summaries are identical in all 12 cases.

### C2 — 366 cases, 349 identical at head; all explained: True

- Output: `/tmp/cran_prep_displays_se_docs_2026-11-12/cumulative/C2/RESULT.md`
- Walls: Whole-process walls, /usr/bin/time, start / head / head-minus-C2: run_case.R 23.84 / 23.36 / 23.74 s; v_run.R fits 52.99 / 53.99 / 51.93 s; v_run.R val 2.85 / 2.88 / 2.61 s; v_probe.R 3.27 / 3.16 / 3.28 s. Compares: compare.R 0.19 s, v_compare.R 0.23 s, probe_compare.R 0.13 s, aux 0.21 s.
- Implementer IRD_out (effect.threshold 1.5) and IRD_out2 (3) fits: only the two C2 IRD warning texts differ. Fit (timing leaves dropped), threshold_config, recorded effect.threshold (0.01), messages and other warnings are identical. Cause: C2 (2984dbee), its intended rewording.
- Implementer direct validator call IRD_out: the C2 validator message differs. Cause: C2.
- Verifier fits IRD_1.5, IRD_3, IRD_hr1.5, IRD_1.5_dina, IRDbin_1.5, SRV_IRD_1.5: only the C2 messages differ (2 per case; 1 for SRV_IRD_1.5, where only the validator fires). Fits, screening values, args_call_all, messages, non-C2 warnings and the DINA errors are identical. Cause: C2.
- Verifier validator grid: 8 of 336 cells differ, exactly {IRD, threshold > 1, user-set} across the 4 outcome types. Return values and non-C2 warnings are identical in every cell. Cause: C2.
- No other part-5 fix (A1, A2, A3, B1, C1, C3, C4, C5, C6, C7) changes any of the 366 cases. The head compare outputs are byte-identical to C2's original start-vs-C2 outputs apart from the wall line, and the head results are identical to the original C2-alone results (aux_head_vs_c2only.txt).
- Confirmation: C2.patch R/ hunks reverse-applied cleanly to head_fixes ('Applied patch R/forestsearch_main.R cleanly'; hunks 2 and 3 at offset -30 because C4 edited comments in the same file). Re-run on that tree: all 366 cases are identical to the start tree (validator grid 336/336, probe identical).
- Normalisation, as in the originals: only $find.grps$time_search and $minutes_all are dropped from each fit, and the pre-drop leaf diff shows only those two leaves. No ggplot or gt objects are present, so no S7_class or closure normalisation was done. The probe matrix was compared whole. The probe-compare script was reconstructed, because the verifier saved only its output; it reproduces the original two output lines.

### C3 — 37 cases, 9 identical at head; all explained: True

- Output: `/tmp/cran_prep_displays_se_docs_2026-11-12/cumulative/C3/RESULT.md`
- Walls: run.R (13 GRF GLM cases): start 8.32 s, start2 8.31 s, head 8.27 s, minus_A1 8.56 s, minus_C3 8.38 s, minus_A1_C3 8.45 s. run_surv.R: 5.72 / 5.64 / 7.69 / 5.47 / 5.58 / 5.52 s. Verifier run_unaff.R (22 GRF cases + survival): 10.92 / 9.91 / 9.89 / 10.56 / 10.57 / 10.75 s. Comparisons: compare.R GRF 0.19 s, survival 1.61 s; compare_unaff.R 1.51 s. All timed with /usr/bin/time; each tree ran in its own Rscript --vanilla process with timeout 120. With both A1 and C3 reversed (minus_A1_C3), all 37 cases are identical to start. Survival is compared without sg10_out$_formats (gt fmt_number closures), the same leaf the originals dropped, which also differs between start and start2.
- A1 (f4c7dd2d, sg_tables() uses the fit's own orientation), 16 cases. Implementer: binary_tree, binary_frontier, binary_harm_notation, binary_nosg, continuous_tree, continuous_frontier. Verifier: bin_tree.{none,harm,benefit}, bin_frontier, bin_RD, bin_nosg, con_tree.{none,harm,benefit}, con_frontier. All are adverse_outcome=FALSE fits. Only tab_estimates _data differs, in the effect column alone: OR becomes 1/OR, MD becomes -MD, RD becomes -RD. At start, sg_tables() used make_effect_estimator()'s default adverse_outcome=TRUE; A1 passes the fit's FALSE. Confirmed by reverse-applying A1's R hunks to head (clean, offsets -22/-22/-45/-45): with A1 reversed, all 16 cases are whole- and HTML-identical to start. bin_RD is not B1: a grf_glm_result has no args_call_all, so the RD is unadjusted, and B1 changes only the adjusted branch.
- C3 itself (c3323a59), 12 count cases with a note. Implementer: count_log_offset, count_log_adverse, count_identity, count_frontier. Verifier: cnt_log_off.{none,harm,benefit}, cnt_id_off, cnt_log_adv, cnt_log_nooff, cnt_frontier, cnt_logRR. Only _source_notes differs; this is the intended scale label and the same result as in C3's own gate. The head's notes are identical to C3's fixed-tree RDS. With C3 reversed (clean, offset -34), these 12 cases are identical to start.

## After the parts: articles, per code fix (full)

### A1

A1's changed lines are in sg_tables() only, at head_fixes :850-869 (GRF branch) and :1038-1055 (inside if (is_glm)). No function in R/ calls sg_tables. Of the 13 .qmd/.Rmd documents under vignettes/ at 835fac1b, only forestsearch-2024.Rmd calls it (chunk fs-results), on a survival forestsearch() fit, so neither changed block runs. The article's chunks were run on start and head_fixes with its own arguments and seed (sequential plan): the counted changed statements executed 0 times, positive controls confirmed the counters fire, and the table data is identical; the only differences are gt closure environments and timing leaves, which a start-vs-start control shows too. No article reaches the changed lines and no article output can change. In quarto/ (not a gate), 10 files would change if re-rendered: each calls sg_tables() on a grf.subg.harm.glm() result at its default adverse_outcome = FALSE. Only one is outside _archive/_legacy: actg175/template_actg175_continuous.qmd:182. Report: /tmp/cran_prep_displays_se_docs_2026-11-12/articles/A1.md; scripts and outputs in /tmp/cran_prep_displays_se_docs_2026-11-12/articles/A1_work/.

| Document | Reaches | Calls and branch | Can the output change |
|---|---|---|---|
| `vignettes/articles/actg175-glm-outcomes.qmd` | no | 20 R chunks, all evaluated. Package calls: dina, forestsearch, forestsearch_bootstrap_dofuture, fs_declaration_table, fs_in_subgroup, mr_cost_table, mr_method_comparison, mr_values. No sg_tables() call, and no package function calls sg_tables (its reverse call closure at head_fixes is {sg_tables}). | no |
| `vignettes/articles/biomarker_effects.qmd` | no | 4 R chunks, all evaluated (plus one css chunk). No forestsearch function is called. | no |
| `vignettes/articles/causal_effects_brief_review.qmd` | no | 4 R chunks, all evaluated. No forestsearch function is called. | no |
| `vignettes/articles/extreme_subgroups.qmd` | no | 36 R chunks, all evaluated. Package calls: calibrate_cens_adjust, calibrate_k_treat, check_censoring_dgm, create_summary_table, forestsearch, generate_aft_dgm_flex, gg_forest, simulate_from_dgm. No sg_tables() call; none of these reaches it. | no |
| `vignettes/articles/forestsearch-2024.Rmd` | no | 26 R chunks, 25 evaluated (oob-cv has eval=FALSE). This is the only document that calls sg_tables(): chunk fs-results (lines 412-433), with default arguments, on fs from forestsearch-main (outcome_type = "survival", gbsg, seedit 8316951). The A1 lines are not reached on either branch: fs has class c("forestsearch","list"), so the GRF block (:833-984) is skipped; args_call_all$outcome_type is "survival", so is_glm is FALSE and the if (is_glm) block holding :1038-1055 is skipped. Measured on both trees (sequential plan instead of the article's multisession): the counted A1 statements ran 0 of 6 times on head and 0 of 2 on start, with no make_effect_estimator() calls; positive controls fired on both branches. | no. tab_estimates and sg10_out _data are identical() start vs head_fixes (HR 0.69/1.95/0.61). Whole res_tabs differs only in 12 gt formatter closure environments, the same 12 as in the start-vs-start control, and is identical once closures are deparsed. fs differs only in its timing leaves. |
| `vignettes/articles/guohe-nested-continuum.qmd` | no | 22 R chunks, all evaluated. Package calls: forestsearch, fs_mr_inference, guohe_algorithm3. No sg_tables() call. | no |
| `vignettes/articles/methodology.qmd` | no | 7 R chunks, all evaluated (plus dot chunks). No forestsearch function is called. | no |
| `vignettes/articles/mrct-consistency-simulations.qmd` | no | 38 R chunks, 34 evaluated (own-cells, own-campaign, own-summary and campaign-code have eval=FALSE). Package calls: ahr_* functions, analysis_covariates, forestsearch, get_dfpred, inject_mrct_structure, mr_estimates_table, mr_values, mrct_* functions, scale_dgm_hazard. No sg_tables() call; none of these reaches it. | no |
| `vignettes/articles/mrct-consistency.qmd` | no | 19 R chunks, all evaluated. Package calls: ahr_design, ahr_limit, ahr_population, compute_frontier_cis, forestsearch, fs_in_subgroup, get_dfpred, inject_mrct_structure, mr_estimates_table, mr_values, mrct_four_sets, mrct_record, mrct_trial_analysis, plot_pareto_frontier, scale_dgm_hazard. No sg_tables() call. | no |
| `vignettes/articles/paper_simulations.qmd` | no | 47 R chunks, 35 evaluated (12 have eval=FALSE). Package calls: build_classification_table, build_estimation_table, calibrate_k_inter, compute_detection_probability, compute_dgm_cde, create_summary_table, forestsearch, format_oc_results, generate_detection_curve, interpret_estimation_table, plot_detection_curve, render_reference_table, run_simulation_analysis, setup_gbsg_dgm, simulate_from_dgm, validate_k_inter_effect. No sg_tables() call. | no |
| `vignettes/articles/survival-post-selection.qmd` | no | 17 R chunks, all evaluated. Package calls: dina, forestsearch, forestsearch_bootstrap_dofuture, mr_cost_table, mr_estimates_table, mr_method_comparison, summarize_bootstrap_results. No sg_tables() call. | no |
| `vignettes/articles/treatment_effect_definitions.qmd` | no | 5 R chunks, 4 evaluated (worked-example has eval=FALSE). No forestsearch function is called. | no |
| `vignettes/forestsearch.qmd` | no | 17 R chunks, all evaluated. Package calls: dina, forestsearch, fs_declaration_calibration. No sg_tables() call. | no |

**`quarto/` (read-only, not a gate):** 159 non-HTML files under quarto/ call sg_tables() (grep at HEAD; quarto/ is unchanged since b06de4e8). 10 of them would show a different effect column if re-rendered. Each calls sg_tables() on a grf.subg.harm.glm() result built without adverse_outcome, so the result stores that function's default, FALSE. The change is to the effect column in the ITT and subgroup rows; Mean(C), Mean(T) and Diff stay raw.

Continuous, MD sign flips:
- applications/actg175/template_actg175_continuous.qmd:182 (the only one not in _archive or _legacy)
- applications/actg175/_archive/20260520_actg175_continuous_analysis.qmd:149
- applications/actg175/_archive/20260520_analysis_actg175_continuous.qmd:1088
- applications/actg175/_archive/20260520_analysis_actg175_continuous_effMaxSG-both.qmd:1112
- applications/actg175/_archive/20260520_analysis_actg175_continuous_hrMaxSG-both.qmd:1088
- applications/actg175/_archive/20260520_analysis_actg175_continuous_hrMaxSG-pareto.qmd:1088
- applications/actg175/_archive/20260520_analysis_actg175_continuous_hr.qmd:1088

Binary, the OR becomes the OR of the non-event:
- applications/actg175/_archive/20260520_analysis_actg175_binary_effMaxSG-both.qmd:1152
- applications/actg175/_archive/20260520_actg175_binary_analysis_working.qmd:142
- simulations/_legacy/actg175_binary_harm_simulations_1simexample.qmd:292

The other 149 are unchanged:
- survival fits: the gbsg files, qc, resampling, mrct_slides;
- count fits: count_data_demo and the cox_poisson files;
- binary or continuous forestsearch() fits with adverse_outcome = TRUE;
- binary fits with adverse_outcome unset, which resolves to TRUE: the 106 simulation files with sg_tables(fs_single), and _utils/restructure_bootstrap_chunks.R, which only writes those calls as strings.

### A2

No article reaches A2's changed lines, so no article's output can change. The changed lines are in R/forestsearch_cross_validation.R: :1517-1520, :1540 and :1597 at head (:1503-1506, :1526 and :1583 in a2cd36db). All of them sit inside forestsearch_KfoldOut()'s `if (outall)` block. In R/, the only callers are forestsearch_Kfold() (:621) and forestsearch_tenfold() (:1188), and both pass a literal outall = FALSE; outall is never reassigned inside the function. A static call graph over 786 functions agrees: only those two functions reach forestsearch_KfoldOut. A static sweep covered the evaluated chunks and inline code of all 13 documents; every chunk parses, and there are no child documents or source() calls. The only evaluated route in is forestsearch-2024.Rmd's kfold-cv chunk (forestsearch_tenfold -> outall = FALSE). Its one outall = TRUE call is in oob-cv, which has eval = FALSE, and that fit is survival. I ran the article's own call on both trees, each in a separate process with a sequential plan. The trace shows 1 entry with outall FALSE, and identical(start, head) is TRUE (timing_minutes dropped). Walls: start 8.43 s, head 9.13 s. Statement: /tmp/cran_prep_displays_se_docs_2026-11-12/articles/A2.md; scripts and logs are in /tmp/cran_prep_displays_se_docs_2026-11-12/articles/A2_work/.

| Document | Reaches | Calls and branch | Can the output change |
|---|---|---|---|
| `vignettes/forestsearch.qmd` | no | 17/17 chunks evaluated; package calls: forestsearch, fs_declaration_calibration. No evaluated call reaches forestsearch_KfoldOut(). | no |
| `vignettes/articles/actg175-glm-outcomes.qmd` | no | 20/20 evaluated; forestsearch, forestsearch_bootstrap_dofuture, fs_declaration_table, fs_in_subgroup, mr_cost_table, mr_method_comparison, mr_values. A GLM article, but it makes no cross-validation call. | no |
| `vignettes/articles/biomarker_effects.qmd` | no | 4/4 evaluated; no package function calls. | no |
| `vignettes/articles/causal_effects_brief_review.qmd` | no | 4/4 evaluated; no package function calls. | no |
| `vignettes/articles/extreme_subgroups.qmd` | no | 36/36 evaluated; calibrate_cens_adjust, calibrate_k_treat, check_censoring_dgm, create_summary_table, forestsearch, generate_aft_dgm_flex, gg_forest, simulate_from_dgm. None reaches forestsearch_KfoldOut(). | no |
| `vignettes/articles/forestsearch-2024.Rmd` | no | 25/26 evaluated. Chunk kfold-cv (line 598): forestsearch_tenfold() -> forestsearch_KfoldOut(outall = FALSE) at R/forestsearch_cross_validation.R:1188 -> if (outall) at :1507 not taken, so changed lines :1517-1520, :1540 and :1597 (head numbering; :1503-1506, :1526, :1583 in a2cd36db) do not run. The trace on the article's own call shows 1 entry with outall FALSE on both trees. The only outall = TRUE call (line 660) is in chunk oob-cv, which has eval = FALSE, and its fit is survival. | no: identical(start, head) is TRUE for the article's own forestsearch_tenfold() call (seeds 8316951/8316951L, sims 1, Kfolds 2, sequential; only timing_minutes dropped) |
| `vignettes/articles/guohe-nested-continuum.qmd` | no | 22/22 evaluated; forestsearch, fs_mr_inference, guohe_algorithm3. | no |
| `vignettes/articles/methodology.qmd` | no | 7/7 evaluated; no package function calls. | no |
| `vignettes/articles/mrct-consistency-simulations.qmd` | no | 34/38 evaluated (own-cells, own-campaign, own-summary, campaign-code have eval=FALSE); ahr_*, analysis_covariates, forestsearch, get_dfpred, inject_mrct_structure, mr_estimates_table, mrct_* (build_dgm, case_study, core_summary, record, region_sims, select_case, trial_analysis), scale_dgm_hazard. None reaches forestsearch_KfoldOut(). | no |
| `vignettes/articles/mrct-consistency.qmd` | no | 19/19 evaluated; ahr_design, ahr_limit, ahr_population, compute_frontier_cis, forestsearch, fs_in_subgroup, inject_mrct_structure, mr_estimates_table, mr_values, mrct_four_sets, mrct_record, mrct_trial_analysis, plot_pareto_frontier, scale_dgm_hazard. | no |
| `vignettes/articles/paper_simulations.qmd` | no | 35/47 evaluated; build_*_table, calibrate_k_inter, compute_detection_probability, compute_dgm_cde, create_summary_table, forestsearch, format_oc_results, generate_detection_curve, interpret_estimation_table, plot_detection_curve, render_reference_table, run_simulation_analysis, setup_gbsg_dgm, simulate_from_dgm, validate_k_inter_effect. None reaches forestsearch_KfoldOut(). | no |
| `vignettes/articles/survival-post-selection.qmd` | no | 17/17 evaluated; forestsearch, forestsearch_bootstrap_dofuture, mr_cost_table, mr_estimates_table, mr_method_comparison, summarize_bootstrap_results. forestsearch_Kfold() appears only in prose. | no |
| `vignettes/articles/treatment_effect_definitions.qmd` | no | 4/5 evaluated (worked-example has eval=FALSE); no package function calls. | no |

**`quarto/` (read-only, not a gate):** None of the quarto/ sources would change if re-rendered. 68 source files (.qmd/.Rmd/.rmarkdown; 115 files mention KfoldOut, including 45 .html and 2 .md) call forestsearch_KfoldOut(outall = TRUE). They fall into five groups:
- 27 survival fits (GBSG applications, GuoHe, qc): the new value is assigned but never used.
- 5 GBSG Cox+Poisson fits (count IRR): the rate estimator takes no orientation.
- 10 ACTG175 binary-OR fits with forestsearch(adverse_outcome = TRUE).
- 5 ACTG175 continuous-MD _archive/20260520 fits with forestsearch(adverse_outcome = TRUE) set explicitly.
- 21 ACTG175 binary-OR simulations: the sourced _fsparams bundles do not set adverse_outcome, so forestsearch()'s rule stores TRUE.
In the last three groups the stored TRUE is the same value as the estimator default TRUE that the start code used. No KfoldOut caller passes or reassigns adverse_outcome = FALSE, and no GLM caller loads a cached CV object.

### A3

Only one article reaches A3's changed lines, and its output cannot change. That article is mrct-consistency.qmd, chunk frontier, which calls compute_frontier_cis(fs, n_splits = 1000L, seed = seed_rep5). Its survival fit from forestsearch() stores adverse_outcome = FALSE, so the new is.null() fallback is not taken. make_effect_estimator() receives FALSE on both trees, and its Cox branch does not read the value anyway. Running the article's chunks with its own arguments and seed on start and head_fixes in separate processes gave identical() ci_front and plot layer data; labels were identical with S7_class dropped. The other 12 documents do not name compute_frontier_cis() or compare_selection_rules(), its only namespace caller, which nothing in the package calls. No document uses source() or child documents. Report: /tmp/cran_prep_displays_se_docs_2026-11-12/articles/A3.md. Runs, scripts and logs: /tmp/cran_prep_displays_se_docs_2026-11-12/articles/A3_runs/.

| Document | Reaches | Calls and branch | Can the output change |
|---|---|---|---|
| `vignettes/articles/actg175-glm-outcomes.qmd` | no | 20 evaluated chunks + 21 inline; forestsearch functions: forestsearch, forestsearch_bootstrap_dofuture, fs_declaration_table, fs_in_subgroup, mr_cost_table, mr_method_comparison, mr_values. None is or reaches compute_frontier_cis() in the static call graph (its only namespace caller is compare_selection_rules(), which nothing in the package calls). | no |
| `vignettes/articles/biomarker_effects.qmd` | no | 4 evaluated chunks, display tables only; no forestsearch function named. | no |
| `vignettes/articles/causal_effects_brief_review.qmd` | no | 4 evaluated chunks (gt tables); no forestsearch function named. | no |
| `vignettes/articles/extreme_subgroups.qmd` | no | 36 evaluated chunks + 7 inline; calibrate_cens_adjust, calibrate_k_treat, check_censoring_dgm, create_summary_table, forestsearch, generate_aft_dgm_flex, gg_forest, simulate_from_dgm. None reaches compute_frontier_cis(). | no |
| `vignettes/articles/forestsearch-2024.Rmd` | no | 25 evaluated chunks (oob-cv eval=FALSE excluded) + 2 inline; forestsearch, forestsearch_bootstrap_dofuture, forestsearch_tenfold, grf.subg.harm.survival, sg_tables, summarize_bootstrap_*, cv_metrics_tables, plot/render forest helpers. None reaches compute_frontier_cis(). | no |
| `vignettes/articles/guohe-nested-continuum.qmd` | no | 22 evaluated chunks + 18 inline; forestsearch, fs_mr_inference, guohe_algorithm3. None reaches compute_frontier_cis(). | no |
| `vignettes/articles/methodology.qmd` | no | 7 evaluated chunks, base R / ggplot2 only; no forestsearch function named. | no |
| `vignettes/articles/mrct-consistency-simulations.qmd` | no | 34 evaluated chunks (own-cells, own-campaign, own-summary, campaign-code eval=FALSE excluded) + 78 inline; ahr_*, inject_mrct_structure, mrct_build_dgm, mrct_case_study, mrct_region_sims, mrct_select_case, mrct_trial_analysis, mr_estimates_table, forestsearch, etc. None reaches compute_frontier_cis(). | no |
| `vignettes/articles/mrct-consistency.qmd` | yes | Chunk case-fs: forestsearch(... survival default, no adverse_outcome passed) resolves adverse_outcome = FALSE (R/forestsearch_main.R:1877-1879 head) before capturing args_call_all (:1961), so the fit stores FALSE. Chunk frontier: compute_frontier_cis(fs, n_splits = 1000L, seed = seed_rep5) with 15 frontier rows, so it reaches the new block R/frontier_cis.R:144-147 (head). The if (is.null(adverse_outcome)) is FALSE and the fallback is not taken. :158 passes FALSE to make_effect_estimator(outcome_type = 'survival', effect_measure = 'HR'), whose Cox branch does not read it. At start, :150 args$adverse_outcome %\|\|% TRUE also gives FALSE. Traced: A3 if reached 1x per call, fallback 0x; make_effect_estimator args identical on both trees. | No. Ran case-read/case-region/case-fs/frontier with the article's own arguments and seed (613800005) on start and head_fixes in separate processes. ci_front (15x12) was identical() start vs head, and head code vs start code on the same fit. ggplot_build layer data were identical(); labels were identical with the S7_class closure dropped. The fits were identical except 5 timing leaves (dropped, named in A3.md). |
| `vignettes/articles/paper_simulations.qmd` | no | 35 evaluated chunks (12 eval-false excluded) + 4 inline; run_simulation_analysis, setup_gbsg_dgm, simulate_from_dgm, forestsearch, calibrate_k_inter, table/plot builders. None reaches compute_frontier_cis(). | no |
| `vignettes/articles/survival-post-selection.qmd` | no | 17 evaluated chunks + 7 inline; forestsearch, forestsearch_bootstrap_dofuture, mr_cost_table, mr_estimates_table, mr_method_comparison, summarize_bootstrap_results. None reaches compute_frontier_cis(). | no |
| `vignettes/articles/treatment_effect_definitions.qmd` | no | 4 evaluated chunks (gt tables; worked-example eval=FALSE excluded); no forestsearch function named. | no |
| `vignettes/forestsearch.qmd` | no | 17 evaluated chunks + 34 inline; forestsearch, fs_declaration_calibration. None reaches compute_frontier_cis(). | no |

**`quarto/` (read-only, not a gate):** No file under quarto/ would change output if re-rendered. At 835fac1b, 125 tracked non-report files (124 .qmd plus _utils/restructure_bootstrap_chunks.R, which writes the call as a string) name compute_frontier_cis( or compare_selection_rules(. A3 changes output only for a continuous fit that stores no adverse_outcome. The 12 continuous files (applications/actg175/analysis_actg175_continuous_compare_all.qmd and 11 _archive/20260520_analysis_actg175_continuous*.qmd) all pass adverse_outcome = TRUE explicitly, so it is stored, and they compute their fits in the document with no readRDS. In the 104 binary files the fallback (TRUE) equals the old TRUE, including cached bundles read with readRDS in simulations/actg175/binary_legacy. In the 9 survival/default files (gbsg, dina, _utils) the Cox estimator ignores the value, and the dina hits are prose, not calls.

### B1

No article reaches B1's changed lines, so no rendered output can change. Nothing was rendered and no start-vs-head comparison was needed.

B1's executable changes are the adjusted branch at R/glm_effect_estimators.R:567-577 and the new helper .fs_rd_gcomp_se() at :611-666; the other changed lines are comments. The changed code runs only when all of these hold: make_effect_estimator(outcome_type="binary", effect_measure="RD") is used, with ps_adjust_method "none" and a non-empty adjust_covariates (n_adj > 0); tier 1 fails; and the tier-2 logistic fit converges.

Every make_effect_estimator() call site passes the user's adjust_covariates (default NULL), a stored copy of it, or NULL. Nothing in R/ creates a non-NULL value. None of the 13 documents (all evaluated chunks and inline code; no child documents or source() calls) supplies adjust_covariates or RD. They are survival throughout, or binary OR / continuous MD in actg175-glm-outcomes.qmd.

A run-time trace of actg175's four binary fits and the DINA bootstrap and displays confirmed this: 31 estimator constructions, all binary/OR/NULL, with .estimate_rd() called 0 times and .fs_rd_gcomp_se() 0 times. The write-up is /tmp/cran_prep_displays_se_docs_2026-11-12/articles/B1.md; scripts and evidence are in /tmp/cran_prep_displays_se_docs_2026-11-12/articles/B1_work/.

| Document | Reaches | Calls and branch | Can the output change |
|---|---|---|---|
| `vignettes/forestsearch.qmd` | no | 17 R chunks, all evaluated. forestsearch() x3 (fit-fs, fit-dina, fit-grf) with no outcome_type, so survival; fs_declaration_calibration(). The test at forestsearch_main.R:2288 (outcome_type != "survival") is FALSE, so no GLM estimator closure is built (Cox path). No adjust_covariates anywhere. .estimate_rd() is not reached. | no |
| `vignettes/articles/actg175-glm-outcomes.qmd` | no | 20 R chunks, all evaluated. forestsearch(outcome_type="binary", effect_measure="OR") x4 (fit-bin, fit-s3-fs/-dina/-grf), forestsearch(continuous, MD) x1, plus fs_in_subgroup, mr_values, fs_declaration_table, forestsearch_bootstrap_dofuture, mr_method_comparison and mr_cost_table. adjust_covariates is absent (NULL). OR goes to the fast unadjusted OR closure (forestsearch_main.R:2322-2330) or .estimate_or() (glm_effect_estimators.R:440-443); MD goes to the lm estimator. .estimate_rd() is never called, and even if it were, n_adj=0 would take the unadjusted branch :557-565. Run-time trace on the head tree (fit-bin, fit-s3-fs, fit-s3-dina, fit-s3-grf, and the DINA bootstrap with nb_boots=25 plus mr_values and mr_method_comparison): every make_effect_estimator() construction was binary/OR/NULL (31 in all); .estimate_rd() was called 0 times and .fs_rd_gcomp_se() 0 times. | no |
| `vignettes/articles/biomarker_effects.qmd` | no | 4 R chunks, all evaluated: static tables only. No package function is called. | no |
| `vignettes/articles/causal_effects_brief_review.qmd` | no | 4 R chunks, all evaluated: setup and static tables. No package function is called. | no |
| `vignettes/articles/extreme_subgroups.qmd` | no | 36 R chunks, all evaluated. Calls calibrate_k_treat, generate_aft_dgm_flex, simulate_from_dgm, check_censoring_dgm, calibrate_cens_adjust, create_summary_table, gg_forest, and S3 summary/print/plot. It is an AFT/Cox survival simulation with no forestsearch() call and no binary make_effect_estimator(). | no |
| `vignettes/articles/forestsearch-2024.Rmd` | no | 25 R chunks evaluated, 1 eval=FALSE (oob-cv). Calls forestsearch(outcome_type="survival", effect_measure=NULL, ps_adjust_method="none") with no adjust_covariates, plus grf.subg.harm.survival, sg_tables, forestsearch_bootstrap_dofuture, forestsearch_tenfold, cv_metrics_tables, plot_subgroup_results_forestplot(outcome_type=NULL), plot_sg_weighted_km and plot_km_band_forestsearch; inline code calls figure_note(). Survival throughout: is_glm is FALSE in sg_tables (summary_utility_functions.R:1025-1027), the bootstrap (:296), CV (:1510) and the forest plot (:348). | no |
| `vignettes/articles/guohe-nested-continuum.qmd` | no | 22 R chunks, all evaluated. Calls guohe_algorithm3 and fs_mr_inference(spec: survival/HR/adjust_covariates=NULL), plus plot. There is no forestsearch() call. MR and Guo-He fit their own Cox models and do not use .estimate_rd(). | no |
| `vignettes/articles/methodology.qmd` | no | 7 R chunks, all evaluated: setup, a base-R power curve, static tables and session info. No package function is called. The forestsearch() listings at :1138 and :1247 are plain code blocks, not chunks. | no |
| `vignettes/articles/mrct-consistency-simulations.qmd` | no | 34 R chunks evaluated, 4 eval=FALSE. Calls mrct_trial_analysis, mrct_region_sims and mrct_case_study (their fs_args carry no outcome_type, effect_measure or adjust_covariates), mr_estimates_table, inject_mrct_structure, mrct_build_dgm and ahr_*. All survival; mr_estimates_table takes its survival branch (:214). No adjust_covariates anywhere. | no |
| `vignettes/articles/mrct-consistency.qmd` | no | 19 R chunks, all evaluated. Calls forestsearch() (case-fs, survival default), fs_in_subgroup, mrct_four_sets, compute_frontier_cis (adjust_covariates=NULL by design, frontier_cis.R:159), mr_estimates_table, mrct_trial_analysis, mr_values and plot_pareto_frontier. All survival. | no |
| `vignettes/articles/paper_simulations.qmd` | no | 35 R chunks evaluated, 12 eval=FALSE. Calls run_simulation_analysis(fs_params=fs_params) x2; fs_params has no outcome_type, effect_measure or adjust_covariates, so is_glm is FALSE (run_simulation_analysis.R:972/1241/1370). Also setup_gbsg_dgm, simulate_from_dgm, calibrate_k_inter, and tables. GBSG survival only. | no |
| `vignettes/articles/survival-post-selection.qmd` | no | 17 R chunks, all evaluated. Calls forestsearch() x3 (survival default), forestsearch_bootstrap_dofuture, summarize_bootstrap_results, mr_estimates_table, mr_method_comparison and mr_cost_table. Cox paths only; no adjust_covariates. | no |
| `vignettes/articles/treatment_effect_definitions.qmd` | no | 4 R chunks evaluated, 1 eval=FALSE (worked-example): setup and static tables. No package function is called. | no |

**`quarto/` (read-only, not a gate):** No file under quarto/ (git grep at 835fac1b) calls the changed function in a way that would change its output if re-rendered: none combines an "RD" estimator with a non-empty adjust_covariates.
- Files that set adjust_covariates are OR or survival. OR: the actg175 binary_methods/*_adjusted*.qmd files (fs_params from _fsparams/fs*_params.R, all effect_measure="OR"), fs_mr_inference_minimal_test.qmd, verify_glm_adjust_layers.R and fdr_*.R; they reach .estimate_or(), and "RD" appears in them only in comments. Survival: smoke_test_adjust_covariates.R, verify_adjust_selection_change.R, gbsg_consistency_resample.qmd, GuoHe/* and gbsg_020 declcal*.
- Files that build an RD estimator pass no adjust_covariates, so they take the unadjusted branch, which B1 did not change: qc/qc_exhaustive_glm_pipeline.R:98 and qc/test_adverse_outcome.R:129-135, 223-226, 270, 280.
- No file references .estimate_rd, .fs_rd_gcomp_se or logistic_gcomp_adj outside .md and .html files.

### C1

C1 (commit 9ec67b7e) changes only R/format_subgroup_summary_tables.R: the green highlight set in Table 4 of format_subgroup_summary_tables() (start :222, head :224), plus comment lines. No document under vignettes/ at 835fac1b reaches the changed lines (13 documents: 12 under articles/ plus forestsearch.qmd), so no article output can change. The reasons:
- The function is internal and not exported, and it is not an S3 method.
- No package function references it by call, symbol or string (801 functions parsed at head, 800 at start), and no get, match.fun or eval(parse) site in R/ can produce its name.
- No vignette names it, and there are no child documents or source() calls.
- summarize_bootstrap_results() (used by forestsearch-2024.Rmd and survival-post-selection.qmd) returns the unformatted subgroup_summary and never calls the formatter.

The check was a static sweep of every evaluated chunk and inline expression, with a conservative closure over the package call graph (calls, symbols, strings and S3 methods). It ran on both trees in separate Rscript --vanilla processes (walls 2.65 s head, 2.60 s start). Nothing was rendered, and no runtime comparison was needed. This agrees with the implementer's NOTES section 8.

Report: /tmp/cran_prep_displays_se_docs_2026-11-12/articles/C1.md. Script and outputs: /tmp/cran_prep_displays_se_docs_2026-11-12/articles/C1_work/ (static_sweep.R, head_static.txt/.rds, start_static.txt/.rds).

| Document | Reaches | Calls and branch | Can the output change |
|---|---|---|---|
| `vignettes/forestsearch.qmd` | no | 17 chunks, 17 evaluated, 34 inline. forestsearch calls: forestsearch, fs_declaration_calibration. The document never names format_subgroup_summary_tables, and the conservative closure (calls, symbols, strings and S3 methods) does not contain it. | no |
| `vignettes/articles/actg175-glm-outcomes.qmd` | no | 20 chunks, 20 evaluated, 21 inline. forestsearch calls: forestsearch, forestsearch_bootstrap_dofuture, fs_declaration_table, fs_in_subgroup, mr_cost_table, mr_method_comparison, mr_values. No package function references the formatter, so none of these calls reaches it. | no |
| `vignettes/articles/biomarker_effects.qmd` | no | 4 chunks, 4 evaluated, 0 inline. No forestsearch function is called. | no |
| `vignettes/articles/causal_effects_brief_review.qmd` | no | 4 chunks, 4 evaluated, 0 inline. No forestsearch function is called. | no |
| `vignettes/articles/extreme_subgroups.qmd` | no | 36 chunks, 36 evaluated, 7 inline. forestsearch calls: calibrate_cens_adjust, calibrate_k_treat, check_censoring_dgm, create_summary_table, generate_aft_dgm_flex, gg_forest, simulate_from_dgm. Its eval(parse(text = s$id / sg_ids[[gg]])) lines only evaluate subgroup expressions on data. The closure does not contain the formatter. | no |
| `vignettes/articles/forestsearch-2024.Rmd` | no | 26 chunks, 25 evaluated (oob-cv has eval=FALSE), 2 inline. forestsearch calls include forestsearch, forestsearch_bootstrap_dofuture, forestsearch_tenfold, summarize_bootstrap_results, plot_*, render_forestplot, sg_tables. summarize_bootstrap_results() (:501) calls summarize_bootstrap_subgroups() and returns the unformatted subgroup_summary; it never calls format_subgroup_summary_tables(). Lines :573-579 print the raw data.tables original_agreement and factor_presence, which C1 does not change. | no |
| `vignettes/articles/guohe-nested-continuum.qmd` | no | 22 chunks, 22 evaluated, 18 inline. forestsearch calls: fs_mr_inference, guohe_algorithm3. The closure does not contain the formatter. | no |
| `vignettes/articles/methodology.qmd` | no | 7 chunks, 7 evaluated, 0 inline. No forestsearch function is called; knit_code$get appears only in the code-folding hook. | no |
| `vignettes/articles/mrct-consistency-simulations.qmd` | no | 38 chunks, 34 evaluated (own-cells, own-campaign, own-summary and campaign-code are off), 78 inline. forestsearch calls are the ahr_* and mrct_* functions plus inject_mrct_structure, scale_dgm_hazard, mr_estimates_table, analysis_covariates and get_dfpred. getNamespaceExports() only checks a fixed list of names. The closure does not contain the formatter. | no |
| `vignettes/articles/mrct-consistency.qmd` | no | 19 chunks, 19 evaluated, 50 inline. forestsearch calls: ahr_design, ahr_limit, ahr_population, compute_frontier_cis, forestsearch, fs_in_subgroup, inject_mrct_structure, mr_estimates_table, mr_values, mrct_four_sets, mrct_record, mrct_trial_analysis, plot_pareto_frontier, scale_dgm_hazard. The closure does not contain the formatter. | no |
| `vignettes/articles/paper_simulations.qmd` | no | 47 chunks, 35 evaluated (12 have eval false), 4 inline. forestsearch calls: run_simulation_analysis, setup_gbsg_dgm, simulate_from_dgm, calibrate_k_inter, build_*_table, format_oc_results, render_reference_table, the detection-curve functions and others. No bootstrap-subgroup formatting is involved, and the closure does not contain the formatter. | no |
| `vignettes/articles/survival-post-selection.qmd` | no | 17 chunks, 17 evaluated, 7 inline. forestsearch calls: forestsearch, forestsearch_bootstrap_dofuture, mr_cost_table, mr_estimates_table, mr_method_comparison, summarize_bootstrap_results. summarize_bootstrap_results() (:285) calls summarize_bootstrap_subgroups() but not the formatter, and the article prints only summ_fs$table (:287). | no |
| `vignettes/articles/treatment_effect_definitions.qmd` | no | 5 chunks, 4 evaluated (worked-example has eval=FALSE), 0 inline. No forestsearch function is called. | no |

**`quarto/` (read-only, not a gate):** quarto/resampling/test_bootstrap_grf_cuts.qmd calls forestsearch:::format_subgroup_summary_tables() at :314 (chunk test3) and :401 (chunk test5), so it reaches the changed line. The built consistency_dist gt object now styles the >=0.95 row green and bold, a _styles change only. The document never displays that table: it shows only fmt_tables$grf_cut_freq (:338), and its .tc() checks test only the presence and class of grf_cut_freq, basic_stats and agreement. Its rendered output would therefore not change if re-rendered. quarto/guides/forestsearch_diagnostic_scoping.qmd:309 mentions the function in prose only.

### C2

No article reaches C2's changed lines, so C2 cannot change any rendered article.

C2 changes two hunks in R/forestsearch_main.R:
- validator check 4: start :301-310, head :301-317;
- the warned remap: start :2334-2350, head :2371-2395.

Their guards are the same on both trees. To reach either hunk, forestsearch() must get effect_measure = "RD" or "IRD" explicitly, with a user-set threshold above 1; the remap also needs a non-survival outcome. No default resolves to RD or IRD (head :1869-1875).

What was checked:
- All 13 .qmd/.Rmd files at 835fac1b, with eval=FALSE and #| eval: false chunks excluded. There are no child documents, source() calls or data loads.
- A findGlobals call graph from forestsearch() found 26 functions that can reach it.
- Every evaluated call into those functions passes survival (effect_measure NULL), OR, or MD. Bootstrap and CV replays inherit the parent fit's measure.
- "RD" and "IRD" appear nowhere in the 13 documents.

Calls that come closest:
- actg175 fit-cont: MD with effect.threshold = 10. It takes the MD branch, not the RD/IRD branch.
- forestsearch-2024 and paper_simulations: hr.threshold = 1.25 on survival. is_bounded_identity is FALSE and the remap block is skipped.

Nothing was re-run on the two trees, because no document reaches the lines. Nothing was rendered or installed, and nothing was written to the repository.

Report: /tmp/cran_prep_displays_se_docs_2026-11-12/articles/C2.md. Evidence: /tmp/cran_prep_displays_se_docs_2026-11-12/articles/C2_work/ (list_calls.txt, excluded.txt, reach_graph_nostr.txt, reach_graph.txt).

| Document | Reaches | Calls and branch | Can the output change |
|---|---|---|---|
| `vignettes/forestsearch.qmd` | no | forestsearch() in fit-fs, fit-dina and fit-grf: survival, effect_measure NULL, hr.threshold = 1.00. The validator runs with is_bounded_identity FALSE (effect_measure NULL), so check 4 (head :301-317) is false. The remap (head :2371-2395) sits inside outcome_type != "survival" (:2288), so it is skipped. fs_declaration_calibration() (declaration) only reads the fit. | no |
| `vignettes/articles/actg175-glm-outcomes.qmd` | no | forestsearch() in fit-bin, fit-s3-fs and fit-s3-dina (binary, OR, effect.threshold = 1.0) and in fit-s3-grf (binary, OR, threshold not set). fit-cont: continuous, MD, effect.threshold = 10. MD is not in c("RD","IRD"), so validator check 4 is false and the code takes the MD else-if branch (head :2414), not the remap. forestsearch_bootstrap_dofuture() (bootstrap) replays the three OR fits with effect_measure OR. The RD/IRD guard is never true. | no |
| `vignettes/articles/biomarker_effects.qmd` | no | 4/4 chunks evaluated. No package function is called. | no |
| `vignettes/articles/causal_effects_brief_review.qmd` | no | 4/4 chunks evaluated. No package function is called. | no |
| `vignettes/articles/extreme_subgroups.qmd` | no | 36/36 chunks evaluated. Only library(forestsearch). generate_aft_dgm_flex, simulate_from_dgm, calibrate_k_treat, calibrate_cens_adjust, check_censoring_dgm, create_summary_table and gg_forest do not reach forestsearch(). | no |
| `vignettes/articles/forestsearch-2024.Rmd` | no | forestsearch() in forestsearch-main: survival, effect_measure NULL, hr.threshold = 1.25. That is user-set and above 1, but is_bounded_identity is FALSE, and the remap is skipped (survival). forestsearch_bootstrap_dofuture() (bootstrap, NB = 2) and forestsearch_tenfold() (kfold-cv) replay that survival fit. oob-cv (eval=FALSE) is excluded. | no |
| `vignettes/articles/guohe-nested-continuum.qmd` | no | Only library(forestsearch). fs_mr_inference() and guohe_algorithm3() (mr_spec: survival, HR) do not call forestsearch(). | no |
| `vignettes/articles/methodology.qmd` | no | 7/7 chunks evaluated. No package function is called. The full forestsearch() argument listings are plain fenced r display blocks and are not evaluated. | no |
| `vignettes/articles/mrct-consistency-simulations.qmd` | no | mrct_trial_analysis() (bridge-case) and mrct_case_study() (find-case) go through .mrct_analyze_trial() to forestsearch(). mrct_region_sims() (find-sims) calls forestsearch() directly. All are survival, with no effect_measure and hr.threshold = 0.90. mrct_record() only reads CSV records. Four eval=FALSE campaign chunks are excluded. | no |
| `vignettes/articles/mrct-consistency.qmd` | no | forestsearch() (case-fs) and mrct_trial_analysis() (one-call): survival, no effect_measure, hr.threshold = 0.90. The RD/IRD guards are false. | no |
| `vignettes/articles/paper_simulations.qmd` | no | run_simulation_analysis() (run-alt-sims, run-null-sims) goes through .run_fs_analysis_gen() to forestsearch(). fs_params: survival, no outcome_type and no effect_measure, hr.threshold = 1.25. That is above 1, but is_bounded_identity is FALSE and the remap is skipped (survival). run_fs_grf = FALSE, and the GRF arm does not call forestsearch(). 35 of 47 chunks are evaluated. | no |
| `vignettes/articles/survival-post-selection.qmd` | no | forestsearch() in fit-fs (hr.threshold = 1.00), fit-dina (0.90) and fit-grf (threshold not set): all survival, effect_measure NULL. forestsearch_bootstrap_dofuture() (bootstrap) replays them. The RD/IRD guards are false. | no |
| `vignettes/articles/treatment_effect_definitions.qmd` | no | 4 of 5 chunks evaluated (worked-example, eval=FALSE, is excluded). No package function is called. | no |

**`quarto/` (read-only, not a gate):** No file under quarto/ at 835fac1b calls forestsearch(), or a wrapper that reaches it, with effect_measure "RD" or "IRD", so re-rendering any of them would not change its output through C2.

- RD/IRD string literals in code that runs occur only in make_effect_estimator() calls: quarto/qc/test_adverse_outcome.R and quarto/qc/qc_exhaustive_glm_pipeline.R:98. A GLM-family switch() in the old package-code copy 'quarto/legacy calibration files/run_simulation_analysis.R' also has them. None of these calls forestsearch().
- Elsewhere RD/IRD is prose or display text.
- Non-literal effect_measure arguments (fs_effect_measure, pois_effect_measure, grf_effect_measure, em, emeasure, spec$effect_measure, fs_params_list$effect_measure) are assigned OR, MD, IRR or HR.
- No yml, yaml, csv, json or txt file contains the word RD or IRD.

### C3

C3 changes only the GRF CATE source note inside sg_tables()'s grf_glm_result branch (R/summary_utility_functions.R :948-970 and :974-975 at 835fac1b), and the text changes only for count fits. No package function calls sg_tables(): the codetools reverse call graph over all 786 functions finds no caller and there are no string references. So an article reaches the lines only through a direct, evaluated sg_tables() call on a grf_glm_result. Of the 13 documents under vignettes/ at 835fac1b, only forestsearch-2024.Rmd calls sg_tables(). It passes a survival forestsearch() fit, so the branch is not taken. I ran its chunks with counters on both trees: branch 0, CATE block 0, and the positive control fired. The article's tables have no CATE note and are identical between trees once gt's random container id is masked. Verdict: no article reaches C3's lines and no article output can change. Report: /tmp/cran_prep_displays_se_docs_2026-11-12/articles/C3.md. Evidence: /tmp/cran_prep_displays_se_docs_2026-11-12/articles/C3_work/ (sweep.R, sweep_*.txt, run_fs2024.R, fs2024_*.txt/.rds, compare_fs2024.txt).

| Document | Reaches | Calls and branch | Can the output change |
|---|---|---|---|
| `vignettes/forestsearch.qmd` | no | 17 R chunks, all evaluated (34 inline). forestsearch calls: forestsearch, fs_declaration_calibration. No evaluated call to sg_tables(), and no package function calls sg_tables() (codetools reverse call graph over 786 functions: reach set = {sg_tables}; 0 string references). No route to C3's lines. | no |
| `vignettes/articles/actg175-glm-outcomes.qmd` | no | 20 R chunks, all evaluated (21 inline). forestsearch calls: forestsearch, forestsearch_bootstrap_dofuture, fs_declaration_table, fs_in_subgroup, mr_cost_table, mr_method_comparison, mr_values. No sg_tables() call; no grf.subg.harm.glm / grf_glm_result. | no |
| `vignettes/articles/biomarker_effects.qmd` | no | 4 R chunks, all evaluated; no forestsearch function called. | no |
| `vignettes/articles/causal_effects_brief_review.qmd` | no | 4 R chunks, all evaluated; no forestsearch function called. | no |
| `vignettes/articles/extreme_subgroups.qmd` | no | 36 R chunks, all evaluated (7 inline). forestsearch calls: calibrate_cens_adjust, calibrate_k_treat, check_censoring_dgm, create_summary_table, generate_aft_dgm_flex, gg_forest, simulate_from_dgm. None reaches sg_tables(). | no |
| `vignettes/articles/forestsearch-2024.Rmd` | no | 26 R chunks, 1 eval=FALSE (oob-cv), 25 evaluated. The only document that calls the changed function: chunk fs-results (Rmd :412-433) calls sg_tables(fs = fs, ...). fs comes from forestsearch(..., outcome_type = "survival") (:301), which has class c("forestsearch","list") (R/forestsearch_main.R :2974/:3159/:4098). So inherits(fs, "grf_glm_result") (R/summary_utility_functions.R :833 at 835fac1b) is FALSE and the call takes the FORESTSEARCH PATH. C3's lines (:948-970, :974-975, inside the cate_sg block at :947 in the grf branch) are not executed. I ran the article's chunks data-setup, forestsearch-main (sequential plan) and fs-results with counters on both trees: sg_tables calls 1, grf branch 0, CATE block 0. In the positive control (count grf_glm_result), the counters fired 1/1. | no: the article's tables have 0 source notes in both trees (no GRF CATE note); tab_estimates and sg10_out raw HTML are identical once gt's random container id is masked (1 line each) |
| `vignettes/articles/guohe-nested-continuum.qmd` | no | 22 R chunks, all evaluated (18 inline). forestsearch calls: fs_mr_inference, guohe_algorithm3. None reaches sg_tables(). | no |
| `vignettes/articles/methodology.qmd` | no | 7 R chunks, all evaluated; no forestsearch function called. | no |
| `vignettes/articles/mrct-consistency-simulations.qmd` | no | 38 R chunks, 4 eval false, 34 evaluated (78 inline). forestsearch calls: ahr_* (design, design_table, limit, population, solve, truths_table), analysis_covariates, get_dfpred, inject_mrct_structure, mr_estimates_table, mrct_build_dgm, mrct_case_study, mrct_core_summary, mrct_record, mrct_region_sims, mrct_select_case, mrct_trial_analysis, scale_dgm_hazard. None reaches sg_tables(). | no |
| `vignettes/articles/mrct-consistency.qmd` | no | 19 R chunks, all evaluated (50 inline). forestsearch calls: ahr_design, ahr_limit, ahr_population, compute_frontier_cis, forestsearch, fs_in_subgroup, inject_mrct_structure, mr_estimates_table, mr_values, mrct_four_sets, mrct_record, mrct_trial_analysis, plot_pareto_frontier, scale_dgm_hazard. None reaches sg_tables(). | no |
| `vignettes/articles/paper_simulations.qmd` | no | 47 R chunks, 12 eval false, 35 evaluated (4 inline). forestsearch calls: build_classification_table, build_estimation_table, calibrate_k_inter, compute_detection_probability, compute_dgm_cde, create_summary_table, format_oc_results, generate_detection_curve, interpret_estimation_table, plot_detection_curve, render_reference_table, run_simulation_analysis, setup_gbsg_dgm, simulate_from_dgm, validate_k_inter_effect. None reaches sg_tables(). | no |
| `vignettes/articles/survival-post-selection.qmd` | no | 17 R chunks, all evaluated (7 inline). forestsearch calls: forestsearch, forestsearch_bootstrap_dofuture, mr_cost_table, mr_estimates_table, mr_method_comparison, summarize_bootstrap_results. None reaches sg_tables(). | no |
| `vignettes/articles/treatment_effect_definitions.qmd` | no | 5 R chunks, 1 eval false, 4 evaluated; no forestsearch function called. | no |

**`quarto/` (read-only, not a gate):** No quarto/ file would change if re-rendered. Ten files call sg_tables() on a grf.subg.harm.glm() result, so they reach C3's block, but all ten fits are continuous (template_actg175_continuous.qmd and six _archive/20260520_*continuous*.qmd, as grf_result) or binary (two _archive/20260520 binary files as grf_result, and _legacy/actg175_binary_harm_simulations_1simexample.qmd as grf_sim). For those types the note keeps \"(risk-difference scale)\", the same string as before. The count GRF objects in applications/count_data_demo.qmd (grf_pois :349, grf_id :380) are never passed to sg_tables(): its :747 call tabulates the forestsearch(outcome_type = \"count\") fit (:653), which takes the FS path. The GBSG Cox/Poisson files (_archive/2026-05-27_gbsg_survival_cox_poisson.qmd :986 and seven _broken/ files) also tabulate forestsearch() fits. Every other sg_tables() call in quarto/ is on a forestsearch() fit.

## After the parts: C5 examples at the head (full)

All 15 examples that ran in C5 also ran on the head with rc 0 and no errors. The 16th, fs_fdr_report, was NOT RUN because it is over the 120 s cap, as in C5. Every page has the same result as the C5 run.

The head example code is the same as the C5 copy's. tools::Rd2ex(commentDontrun = FALSE) on $S/head_fixes/man and on $S/wt/C5/man (after document()) gives identical code for 16/16 pages, and it also equals the after.R scripts C5 actually ran. The Rd files are byte-identical for 15 pages. sg_tables.Rd differs only by C4's \details section; its \examples block is unchanged.

$S/head_fixes is the same as git archive 835fac1b for R/, man/, DESCRIPTION and NAMESPACE, and R/ and man/ are unchanged at fb80dd8f. It was installed from the copy $S/lib_src_head into $S/lib/head (R CMD INSTALL --no-test-load, rc 0, 17.96 s). Under R_LIBS, $S/lib/head comes first in .libPaths().

Outputs were compared after masking only the random gt table ids and the k-fold "Time:" line. All 15 outputs are identical to C5's after.out. The 6 Rplots.pdf files are byte-identical with dates masked. Walls are within 1.5 s of C5's.

No part-5 code fix changes any of these displays. Each fix's code path is either not reached by these survival-outcome examples or gives the same result (details in RESULT.md section 5): C2 IRD warnings, A2 CV orientation, A3 frontier fallback (the fit stores adverse_outcome), B1 binary RD SE, C1 >=0.95 highlight (format_subgroup_summary_tables() has no callers), and A1/C3 sg_tables for non-survival fits and the GRF count note.

run_subgroup_sims prints nothing. As an extra check, its result was saved from lib/C5 and lib/head; identical() is TRUE once the timing leaves are dropped (t_plan_secs, t_sims_secs, created).

The repository was not written to (git status clean, HEAD fb80dd8f). Report: /tmp/cran_prep_displays_se_docs_2026-11-12/examples_head/RESULT.md. Runs: /tmp/cran_prep_displays_se_docs_2026-11-12/examples_head/<page>/head.{R,out,rc,wall}.

| Page | Result | Wall | Same as the C5 run |
|---|---|---|---|
| `format_bootstrap_table` | RAN, rc 0 (prints the gt table) | 17.75 s (C5 17.52 s) | yes: output identical once the random gt table id is masked |
| `print.fs_kfold` | RAN, rc 0 | 68.54 s (C5 67.07 s) | yes: output identical except the timing line (Time: 1 -> 1.02 minutes) |
| `print.fs_tenfold` | RAN, rc 0 | 14.61 s (C5 14.40 s) | yes: output byte-identical |
| `plot_sg_distribution` | RAN, rc 0 (the same single warning: the pdf device substitutes >= for U+2265) | 1.71 s (C5 1.78 s) | yes: output byte-identical; Rplots.pdf byte-identical with dates masked |
| `print.fs_sg_plot` | RAN, rc 0 | 6.44 s (C5 6.49 s) | yes: output byte-identical; Rplots.pdf identical |
| `plot.fs_sg_plot` | RAN, rc 0 (prints 'No plots available', as in C5) | 6.42 s (C5 7.06 s) | yes: output byte-identical; Rplots.pdf identical |
| `print.fs_forestplot` | RAN, rc 0 | 6.44 s (C5 6.39 s) | yes: output byte-identical; Rplots.pdf identical |
| `plot.fs_forestplot` | RAN, rc 0 | 6.76 s (C5 6.44 s) | yes: output byte-identical; Rplots.pdf identical |
| `run_subgroup_sims` | RAN, rc 0 (prints nothing: the example only assigns and saves to tempdir()) | 13.35 s (C5 13.01 s) | yes: output empty in both runs; extra check saving the result from lib/C5 and lib/head: identical() TRUE once the three timing leaves are dropped (sim_config$t_plan_secs, sim_config$t_sims_secs, created) |
| `sg_consistency_out` | RAN, rc 0 | 6.50 s (C5 6.51 s) | yes: output byte-identical |
| `evaluate_subgroup_consistency` | RAN, rc 0 | 6.66 s (C5 6.38 s) | yes: output byte-identical |
| `evaluate_consistency_twostage` | RAN, rc 0 | 6.39 s (C5 6.34 s) | yes: output byte-identical |
| `fs_fdr_report` | NOT RUN (over the 120 s cap): the example code is identical to C5's, which the cap killed inside its first call (2 calls x 2 nulls x 1000 MR fits); no part-5 code change makes it faster | not run (C5 120.76 s, killed, rc 124) | yes: NOT RUN in both |
| `sg_tables` | RAN, rc 0 (prints the gt table) | 6.89 s (C5 6.89 s) | yes: output identical once the gt table id is masked (A1/C3 change only non-survival fits and the GRF count note; this example is a survival fit with use_grf = FALSE) |
| `pareto_frontier_table` | RAN, rc 0 (the same ggplot 'only one observation' message) | 10.66 s (C5 10.06 s) | yes: output identical once the 2 gt table ids are masked; Rplots.pdf identical (the A3 fallback is not reached: the forestsearch() fit stores adverse_outcome) |
| `summarize_simulation_results` | RAN, rc 0 | 18.81 s (C5 18.08 s) | yes: output byte-identical |
