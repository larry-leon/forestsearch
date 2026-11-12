# REPORT — CRAN preparation, part 4: approved code fixes, 2026-11-12

Task: `dev/tasks/TASK_cran_prep_code_fixes_2026-11-12.md` (committed as received at the start SHA). Run unattended from `~/Documents/GitHub/forestsearch` on branch `cran-prep`, after part 3 (`REPORT_cran_prep_docs_legacy_2026-11-12.md`). Fourteen fixes, each authorized by Larry on 2026-11-12; one fix, one commit.

**Result: 14 PASS, 0 FAILED, 0 BLOCKED.** The default-path gate passes (6/6 calls identical to the start SHA). Full suite: 0 failures, 0 errors. `check()`: 0 errors, 0 warnings, 0 notes. Full per-fix evidence: `REPORT_cran_prep_code_fixes_2026-11-12_detail.md` (beside this report).

## 1. Provenance

| Item | Value |
|---|---|
| Branch | `cran-prep`; `git status --porcelain` empty at the start (gate) |
| Part 3 gate | `dev/reports/REPORT_cran_prep_docs_legacy_2026-11-12.md` exists on `cran-prep` (commit `9e877866`). Its status lines: Step 0 PASS; Parts 1-4 PASS; Part 5 PASS (read-only); Part 6 PASS (0.4.1.9009, suite 0 failures with 652 blocks and 8,478 expectations, `check()` 0/0/0). |
| Start SHA | `76a62a71` — `dev: add CRAN preparation task, part 4` (parent `9e877866`). Exported with `git archive` to `/tmp/cran_prep_code_fixes_2026-11-12/start/`. |
| `DESCRIPTION` version | 0.4.1.9009 at the start → 0.4.1.9010 |
| fs-post-selection | `main` at `5d165acce104a336e11093f1881cfdb572348462`, `git status --porcelain` empty, at the start and at the end |
| fs-glms-interpretable | `main` at `dc94a17e0e93653d5d4a9fb83bacca93bb8f87ea`, `git status --porcelain` empty, at the start and at the end |
| R / tools | R 4.6.2; devtools 2.5.2, roxygen2 8.1.1, testthat 3.3.2, pkgload 1.5.3; survival 3.8-12, grf 2.6.1, glmnet 5.1 |
| Scratch | `/tmp/cran_prep_code_fixes_2026-11-12/` (start export, per-fix scratch copies `wt/F*`, evidence, verifier copies, patches, gate runs, check export). Nothing installed into the user or site library; every scratch reproducer used `pkgload::load_all()`, a sequential plan and `timeout 120`; the main-session red/green runs (`..._apply_fix.sh`) used `timeout 300`. |

**Commits** (every `git add` named explicit paths; nothing fetched, pulled, pushed or merged):

| SHA | Time | Commit | Fix |
|---|---|---|---|
| `76a62a71` | 11:37 | `dev: add CRAN preparation task, part 4` | Step 0 (start SHA) |
| `0dc99947` | 12:13 | `fix: bootstrap summary bins match their labels (A-147)` | F1 |
| `1c2565b3` | 12:13 | `fix: CV fold assignments {1, NA} or {0, NA} no longer stop` | F2 |
| `3612d0ba` | 12:13 | `fix: CV "Exact" means the same subgroup (A-103)` | F4 |
| `4aaab2f0` | 12:14 | `fix: test_hte_from_forestsearch() includes continuous covariates (A-151)` | F5 |
| `f5b65516` | 12:14 | `fix: LASSO intercept is not a selected variable (A-042)` | F6 |
| `65831b54` | 12:14 | `fix: fpr_calibration() outcome-type fallback is survival (A-107)` | F7 |
| `ea5f82d9` | 12:14 | `fix: duplicate-cut test compares names exactly (A-047)` | F8 |
| `a30d6e7c` | 12:15 | `fix: MRCT default confounders exclude internal spline columns (A-125)` | F9 |
| `f2e8de47` | 12:15 | `fix: count outcomes with GRF are oriented (A-064)` | F10 |
| `25ecaa0f` | 12:15 | `fix: three displays use the fit's own orientation` | F12 |
| `89ad9c19` | 12:16 | `fix: risk-difference standard error in the unadjusted case (A-055)` | F13 |
| `685e0c5d` | 12:16 | `fix: DINA floor: fit, record and replay agree in the forced-gaussian corner` | F14 |
| `a4cce90a` | 12:31 | `fix: messages when subgroup membership is missing (three sites)` | F3 |
| `61c15088` | 12:31 | `fix: DINA orientation sign on the screening route` | F11 |
| `88302212` | 12:34 | `chore: version and NEWS for part 4` | After the fixes, step 3 |
| `a577b59d` | 13:43 | `chore: NEWS wording for part 4` | Report review (§7) |
| (this commit) | | `dev: CRAN preparation report, part 4` | After the fixes, step 7 |

F3 and F11 were committed last because their verifiers finished last. The order changes nothing, because the fixes are independent; every patch applied cleanly on top of the earlier commits.

**How the work was done.**
- **Implement (scratch, one agent per fix).** Each fix was built in a private copy of the start SHA (`git archive 76a62a71` of the package paths) under `/tmp/cran_prep_code_fixes_2026-11-12/wt/F<n>/`. The agent wrote the test first, recorded red on the unfixed copy, applied the fix, recorded green, ran the existing test files that exercise the changed functions, ran the unaffected case (start export and fixed copy, each in its own `Rscript --vanilla` process, `identical()`), and wrote a patch and notes.
- **Verify (one independent adversarial verifier per fix).** It re-applied the patch to fresh copies, re-ran red and green, re-derived every expected value by hand, re-ran the unaffected case and checked scope and callers. 13 fixes were CONFIRMED in the first round. F5 had one blocking problem (non-syntactic categorical names); it was repaired and CONFIRMED by a second verifier.
- **Commit (main session, one fix at a time).** `..._apply_fix.sh` applied only the patch's test file to the repository and ran it (red), then applied the code and ran it again (green). `..._commit_fix.sh` checked that the staged diff equals the gated patch (ignoring `index` and hunk-header lines), printed `--stat` and committed. Both scripts are beside this report; their outputs are summarized per fix in §3.
- **Comment-only additions before commit (F3, F11).** The verifiers of F3 and F11 asked for derivation comments on three expected values. The main session added them (no code or value change) and regenerated the patch before commit: F3 one line (`any.H`); F11 one citation and one `expect_length` comment. The verified patches are kept as `patches/F3.verified.patch` and `patches/F11.verified.patch`. For the other fixes, verifier comment requests are listed as open items: the commits are the gated patches unchanged.
- **After the fixes.** Read-only agents re-ran every fix's unaffected case at the head with all 14 fixes in (§5.2), derived the `paper_simulations` reach per fix (§5.3), and swept the other articles (§5.4). Fourteen agents drafted the per-fix detail sections from the evidence, checked against git; the main session condensed them for §3.
- **Measured wall clocks.** Implement and verify workflow (30 agents, up to 16 at a time): 2,934 s. Post-fix checks (6 agents): 3,036 s. Main-session red/green runs per fix: 2.7–14.0 s (listed in §3). Default-path gate: 18.1 s (start) and 17.9 s (head). `devtools::document()`: 12.3 s. Suite 730 s, `check()` 1,092 s (§6). Scratch reproducers ran under `timeout 120` and the main-session red/green runs under `timeout 300`; no gate run hit its cap. Several whole-file runs of the existing `test-mr-inference.R` did hit the 120 s cap in the scratch copies; it was then run block by block, and it passes in the full suite (182.4 s there).

## 2. Status per fix

| Fix | Item | Status | Commit | One line |
|---|---|---|---|---|
| F1 | A-147 | **PASS** | `0dc99947` | Both `cut()` calls take `right = FALSE`; edge values land in the bin their label names. |
| F2 | CV {1, NA} / {0, NA} | **PASS** | `1c2565b3` | Table built before the branch; two-level branch only for a 2×2 table; no more "subscript out of bounds". |
| F3 | NA membership messages | **PASS** | `a4cce90a` | F3a stops in exactly the same cases with a message; F3b and F3c warn once; returned values identical. |
| F4 | A-103 | **PASS** | `3612d0ba` | One-factor "Exact" requires the fold's subgroup to be that factor alone. |
| F5 | A-151 | **PASS** | `4aaab2f0` | Cut labels map to their columns; continuous covariates now used; unresolved labels and an empty set warn. Repaired once. |
| F6 | A-042 | **PASS** | `f5b65516` | `(Intercept)` out of `selected`/`omitted`; full pool kept when LASSO selects nothing; label by outcome type. |
| F7 | A-107 | **PASS** | `65831b54` | Fallback type is `"survival"`; an unrecognised type stops. |
| F8 | A-047 | **PASS** | `ea5f82d9` | GRF/DINA duplicate test and median removal compare names exactly. |
| F9 | A-125 | **PASS** | `a30d6e7c` | `<var>_treat`, `<var>_k`, `<var>_k_treat` excluded by exact name in both MRCT functions. |
| F10 | A-064 | **PASS** | `f2e8de47` | Count outcome negated before the forest when `adverse_outcome = TRUE`. |
| F11 | DINA screening sign | **PASS** | `61c15088` | `tau_sign = .dina_tau_sign(outcome_type, adverse_outcome)` on the `use_dina = TRUE` route. |
| F12 | Display orientation | **PASS** | `25ecaa0f` | Three displays pass the fit's `adverse_outcome` (forestsearch()'s rule if absent). |
| F13 | A-055 | **PASS** | `89ad9c19` | Change of method: unadjusted RD logistic-tier SE is the two-sample binomial SE. |
| F14 | DINA gaussian floor | **PASS** | `685e0c5d` | Floor is the resolved threshold at both sites; corner fit = record = replay; every other configuration unchanged. |
| Version/NEWS | | **PASS** | `88302212`, `a577b59d` | 0.4.1.9010; one NEWS entry per fix; `document()` changed nothing else; four bullets reworded after the report review. |

## 3. The fixes

Line numbers are at the start SHA `76a62a71` and at the head of the fixes `61c15088` (`R/` is identical at `88302212`). The pointers in the task were taken at `28e71ffa`. Each entry gives the site before and after, the red and green runs, how the expected values were derived, the unaffected case and one before/after example. The detail file has the full tables (every expected value, every run, every verifier problem).

### F1 — A-147: bootstrap summary bins match their labels

- **Status:** PASS — commit `0dc99947` (parent `76a62a71`) — `R/summarize_bootstrap_subgroups.R` +4/-2, new `tests/testthat/test-bootstrap-summary-bins.R` +100 (2 files, +104/-2).
- **Site, before → after:** `R/summarize_bootstrap_subgroups.R:417-422` (Pcons `cut()`) and `:444-449` (N_sg `cut()`) at the start SHA `76a62a71` → `:417-423` and `:445-451` at head `61c15088`. Both calls get `right = FALSE` (new lines `:422`, `:450`). `include.lowest = TRUE` stays. The bins go from right-closed `(a, b]` to left-closed `[a, b)`, as in `R/summarize_cv_results.R:764-765`. Pcons site shown; the N_sg site is the same change.

Before:
```r
      include.lowest = TRUE
```
After:
```r
      include.lowest = TRUE,
      right = FALSE
```

- **Red → green:** Main session, repository at `76a62a71` plus the test only: red 6 failures of 8 expectations, 2 passes. Example: `:47:3` `actual` `"0.85-0.95" "0.85-0.95" "0.8-0.85" ...` vs `expected` `">=0.95" "0.85-0.95" "0.85-0.95" ...` (3.0 s). Green after the code hunks: 8 pass, no failures (2.9 s).
- The implementer (red 6 failed/2 passed, 3.04 s; green 8 passed, 2.84 s) and verifier (red `failed= 6 passed= 2`, 3.60 s; green `failed= 0 passed= 8`, 3.08 s) reproduced it.
- **Expected values:**
  - The ten per-value bins (`:43`, `:58`) are the spec's values, derived by hand from the labels ("a-b" holds a <= v < b; header comment `:8-11`). Counts (`:89`, `:92`) are hand tallies. Percents (`:98`, `:99`) are 100 x Count / 5, the function's own definition.
  - Label-order checks (`:76-83`) use the labels in `R/summarize_bootstrap_subgroups.R:420`, `:442` at `76a62a71`.
  - The verifier re-derived all values by hand: every value has a derivation comment and none needs the fixed code. It flagged that the label-order checks pass before and after, and their comment `:75` does not say so.
- **Unaffected case:** Spec: "values not on a bin edge give identical tables." Three synthetic cases (`consistency_wide`, `consistency_default_grid`, `dina_mode`; 400 rows, 332 found; `stopifnot()` that nothing is on a break) ran on the start export and the fixed copy, all 11 return slots compared: `OVERALL identical(start, fixed): TRUE`. Walls: start 2.88 s, fixed 2.87 s, driver 6.02 s.
- The verifier re-ran it (identical, also without normalisation; 3.10 s / 3.09 s, driver 6.46 s). Its grid gave identical bins for all 994 off-edge Pcons and 445 off-edge N_sg values; only the interior breaks moved up one bin (8.46 s / 8.65 s).
- **Existing tests:** none reach the changed lines (`grep`/`git grep` exit 1; only caller `R/summarize_bootstrap_results.R:394`), so none ran. Expectations updated: none. The full suite at `88302212` passed it: `✔ |          8 | bootstrap-summary-bins`.
- **What a user sees:** `consistency_dist` for Pcons 0.95, 0.90, 0.85, 1.00, 0.949 (`evidence/F1/example.txt`):
```
before (start tree)                          after (fixed tree)
4:          0.8-0.85     1      20           4:          0.8-0.85     0       0
5:         0.85-0.95     3      60           5:         0.85-0.95     3      60
6:            >=0.95     1      20           6:            >=0.95     2      40
```
  A consistency of exactly 0.95 moves from "0.85-0.95" to ">=0.95". In `size_dist` a subgroup of 50 moves from "<50" to "50-99".
- **Verification:** CONFIRMED. No repair or re-verify round. Two minor problems, both left as open items:
  - The label-order checks pass on both trees. The suggested comment was not added (the commit equals `patches/F1.patch`).
  - `R/format_subgroup_summary_tables.R:222` highlights `">= 0.95"` (with a space), which never matches `">=0.95"`. It is out of F1's scope and unchanged at `61c15088:222`.
- **Departures:** none, per the implementer. The interpretations (extra test checks, how the unaffected case was read, no roxygen change) are in the detail file.
- Full evidence: detail file, section F1.

### F2 — CV fold assignments {1, NA} or {0, NA} no longer stop

- **Status:** PASS — commit `1c2565b3` (parent `0dc99947`) — `R/forestsearch_cross_validation.R` (+7/-3), `tests/testthat/test-cv-fold-na-assignments.R` (new, 78 lines); 2 files, 85 insertions, 3 deletions.
- **Site, before → after:** `forestsearch_KfoldOut()`, `R/forestsearch_cross_validation.R:1417-1420` (stop at `:1420`) and `:1423-1425` at the start SHA `76a62a71` → `:1447-1454` and `:1458-1459` at head `61c15088`. The table is now built once before the branch, the two-level branch also requires a 2x2 table, and the duplicate `tabit <-` line in `else` (start `:1424`) is removed. Before (start `:1417-1420`):
  ```r
      if (n_sgfound == 2) {
        tabit <- with(df_CV, table(treat.recommend, treat.recommend.original))
        sensH <- tabit[1, 1] / sum(tabit[, 1])
        sensHc <- tabit[2, 2] / sum(tabit[, 2])
  ```
  After (head `:1447-1454`, without the five-line comment and the unchanged `sensH` line):
  ```r
      tabit <- with(df_CV, table(treat.recommend, treat.recommend.original))
      if (n_sgfound == 2 && nrow(tabit) == 2L && ncol(tabit) == 2L) {
  ```
- **Red → green:** Main session, repository at `0dc99947` + test file only: red `.1.2`, 2 passes (premise checks) and 2 errors (test `:57:3`, `:70:3`), ``Error in `[.default`(tabit, 2, 2): subscript out of bounds`` at `R/forestsearch_cross_validation.R:1420:7`, 3.0 s; with the code, green 6 passes, 0 failures, 2.9 s.
  Implementer and verifier reproduced it: red 2 passes / 2 errors (2.95 s; verifier `[ FAIL 2 | WARN 0 | SKIP 0 | PASS 2 ]`, 3.19 s), green 6 passes (2.97 s; verifier `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 6 ]`, 3.08 s).
- **Expected values:** By hand from the test data (100 subjects; `table()` drops the 20 NA folds, leaving 80) and the single-level rule (test header lines 15-23; code comment at start `:1427-1431`): `{1, NA}` → `c(0, 1, 0, 53 / 80)` (line 64), `{0, NA}` → `c(1, 0, 27 / 80, 0)` (line 77). The verifier re-derived both with base-R arithmetic (`TRUE`, `TRUE`).
  - Also two premise checks (`2L`, lines 55, 68, by hand) and two names checks (lines 58, 71, from `names(sens_metrics_original) <-` at start `:1448`). Comments: lines 59-63 and 72-76 for the values, line 54 for the premise check at 55; none of its own at 68, none for the names checks. The verifier flagged no value as lacking a comment.
- **Unaffected case:** Spec: {0, 1} and {0, 1, NA} give identical values. Ran `forestsearch_KfoldOut()` on hand-built `res` (8 cases) and `forestsearch()` + `forestsearch_Kfold(Kfolds = 5L)` on N = 200 data, complete and with 20 NA ages. `identical()` TRUE for both (implementer start vs scratch; verifier start vs `post`). Walls start/fixed: implementer 3.12/3.03 s and 8.47/8.72 s; verifier 3.08/3.08 s and 8.33/8.68 s. Verifier sweep, 84 cases (4.24/5.01 s): 56 identical, the 28 that differ all stopped on the start tree.
- **Existing tests:** All CV-calling test files pass with 0 failures for implementer and verifier (`test-mr-inference.R` in four chunks; `test-input-validation.R` verifier only); warnings and skips predate the fix. Existing expectations updated: none.
- **What a user sees** (`evidence/F2/example.txt`, folds {1, NA}); `forestsearch_tenfold()` now completes too. `forestsearch_Kfold()` did not stop before: it caught the error and returned all-NA `sens_summary`/`find_summary`, with a warning only under `details = TRUE`. It now returns these values (inferred from the code, not measured):
  ```
  # start SHA
  > forestsearch_KfoldOut(res_1na)$sens_metrics_original   # folds {1, NA}
  Error in `[.default`(tabit, 2, 2): subscript out of bounds
  # fixed
  > forestsearch_KfoldOut(res_1na)$sens_metrics_original   # folds {1, NA}
   sens_H sens_Hc   ppv_H  ppv_Hc
   0.0000  1.0000  0.0000  0.6625
  ```
- **Verification:** CONFIRMED; no repair round; the commit equals the gated patch. Three minor problems, no code change: 1x1 and 1x2-with-zero-cell tables also reach the position-labelled fallback (noted as an open item; they stopped before and the spec routes them there); survival `outall = TRUE` stops when `treat.recommend.original` has NA, for every fold shape on both trees (noted, pre-existing); the unmeasured caller-level claim was reworded.
- **Departures:**
  - The guard keeps `n_sgfound == 2` and adds the 2x2 check: `table()` drops NA, so a shape-only guard would send {0, 1, NA} (n_sgfound 3, 2x2 table) to the two-level branch and turn its NA `sens_Hc`/`ppv_Hc` into values (code reading; not run), breaking the named unaffected case. Both earlier reports recommend this variant.
  - The duplicate `tabit <-` line was removed because the table is now built before the branch, as the spec asks.
  - The test has two premise checks and two names checks besides the four values per shape.
  - `test-mr-inference.R` was run in four chunks because the whole file exceeds the 120 s per-run cap.
- Full evidence: detail file, section F2.

### F3 — Missing-covariate messages at three sites

- **Status:** PASS — commit `a4cce90a` (parent `685e0c5d`) — `R/forestsearch_cross_validation.R` +31/-1, `R/forestsearch_helpers.R` +15, `R/run_simulation_analysis.R` +50, new `tests/testthat/test-missing-membership-messages.R` (363 lines); 4 files, +459/-1. No `man/` change.
- **Site, before → after** (start `76a62a71` → head `61c15088`):
  - F3a: FS extractor `R/run_simulation_analysis.R:1129-1131` → `:1172-1177`, new stop `:1175` (quoted); GRF extractor `:1441-1444` → `:1487-1494`, stop `:1491-1492`; new helpers `.stop_na_membership()` `:931-941`, `.sg_def_columns()` `:947-955`.
  - F3b: GRF route `R/forestsearch_helpers.R:1990-1991` → `:1995-2006`, a warning when `tr_est` has NA, before `df.est` is built.
  - F3c: `R/forestsearch_cross_validation.R` Kfold `:552-555` → `:558-564`; tenfold per-simulation counts `:1184-1185` → `:1196-1201`, warning in the caller `:1206-1209` → `:1228-1237`.

  Before (start `:1129-1131`):
  ```r
      hat_H  <- df$sg_hat == 1L

      tp <- sum(true_H & hat_H);   fp <- sum(!true_H & hat_H)
  ```
  After (head `:1172-1175`):
  ```r
      hat_H  <- df$sg_hat == 1L
      # NA membership: tp and fn, or tn and fp, below would be NA and the first
      # `if` would stop.  Stop here with the reason (same cases stop).
      if (anyNA(hat_H)) .stop_na_membership(df$sg_hat, sg_def_txt, names(df))
  ```
- **Red → green:**
  - Main, repository: red at `685e0c5d` with the committed test file, 9 not passed (3 errors, 6 failures), 42 passed, ``Error in `if ((tp + fn) > 0) tp/(tp + fn) else NA_real_`: missing value where TRUE/FALSE needed``, F3b/F3c tests get no warning, 7.3 s; green at `a4cce90a`, 60 passed, 0 failed, 7.3 s.
  - Reproduced: implementer red 9 not passed / 42 passed (7.64 s), green 60 expectations, 0 failed (7.50 s); verifier red `FAIL 9 | PASS 42` (7.40 s), green `FAIL 0 | PASS 60` (7.28 s), and 60, 0 failed on `685e0c5d` + patch (7.14 s).
- **Expected values:**
  - By hand from the hand-built inputs: NA rows per shape (A rows 9-10, B rows 20-22), row counts (tp 8, fn 2, fp 4, tn 26; `size.H = 12`, `size.Hc = 28`), three-valued logic for F3b (rows 3, 7, 11 NA), `x1` NA at rows 5, 12, 33 for F3c (2 x 3 = 6 of 80 for tenfold). Tenfold result names from the `@return` list of `forestsearch_tenfold()` (`R/forestsearch_cross_validation.R:726-777`, start); `call. = FALSE` from the spec.
  - Comments: the verifier flagged `any.H = 1L` (test `:157`) as lacking one; the main session added it at `:156` before the commit. `:139` and `:342` repeat commented checks, and the complete-`x1` block (357, 358, 362) has no comment of its own.
  - The wording is the implementer's design, accepted by the verifier, who re-derived every value by hand and found none that could only come from running the fixed code.
- **Unaffected case:** the Tests bullet, "No returned value changes at any of the three sites". `evidence/F3/unaffected.R`, 15 cases on the start and fixed trees: all 11 complete cases `identical()` TRUE after dropping timing leaves only; all 4 F3a stop shapes stop on both trees. Walls 9.28 s start / 9.22 s fixed (implementer), 8.69 s / 8.84 s (verifier).
- **Existing tests:** 27 whole files, 244 tests, 4351 expectations, 0 failed, 9 skipped (implementer and verifier); `test-mr-inference.R` in parts, 21 tests, 148 expectations, 0 failed; verifier probe over all 89 existing files plus the new one: F3 code fired only in the new file. Expectations updated: none.
- **What a user sees:** `forestsearch(subgroup_method = "grf", grf_selection = "frontier", outcome_type = "binary", ...)`, 500 rows, 30 NA in each of `x1`, `x2`, `x3`. Before: no message, and 34 subjects silently have NA membership (`<NA> 34`). After: the same returned values plus
  ```
  Warning message:
  34 of 500 subjects have no subgroup assignment (treat.recommend is NA): a covariate that defines the GRF-selected subgroup (x1, x2) has missing values.
  ```
- **Verification:** CONFIRMED; the verifier reproduced red, green and the unaffected case; no repair round is recorded, and the only change after verification is one test comment line. Three minor problems: F3c gives the count only → recorded as departure 7, no code change; `.mrct_analyze_trial()` `$warnings` can now hold the F3b text → open item, not changed; `any.H` lacked a comment → comment-only line at test `:156` added by the main session before the commit (`patches/F3.verified.patch` and `patches/F3.patch` differ only there).
- **Departures:**
  - F3a names the covariates: the general bullet asks for them where the definition is at hand; looked up inside `tryCatch`, so which cases stop never depends on them.
  - GRF extractor has one check, not two: NA cannot reach the classification block (start `:1536-1546`), since `m_class$sg_hat` is copied from `harm_ind`.
  - F3b applies to any outcome type: "(binary, continuous or count)" was read as where the NA arises, not as a filter.
  - F3b warns for `df` only: the spec names the analysis frame; `df.predict` and `df.test` do not warn.
  - Tenfold counts inside each simulation and the caller warns: the spec requires the warning in the calling session; the counts do not reach the result.
  - `forestsearch_KfoldOut()` does not warn: Kfold and tenfold call it after collection, so a warning there would break "one warning per exported call".
  - F3c gives the count alone: the F3c bullet asks for "the count"; each fold has its own definition and `resCV` holds only the first two fold labels.
- Full evidence: detail file, section F3.

### F4 — A-103: CV "Exact" means the same subgroup

- **Status:** PASS — commit `3612d0ba` — `R/forestsearch_cross_validation.R` +8/-4 (`CV_sgs()`, one-factor branch only); `tests/testthat/test-cv-exact-match.R` new, +71.
- **Site, before → after:** `R/forestsearch_cross_validation.R:1818-1827` at the start SHA `76a62a71` → `:1852-1865` at head `61c15088`. The containment indicator is kept for `one_match`/`cov1_match`; "Exact" now requires the other slot NA. Depth-2 and no-subgroup branches unchanged.

  Before:
  ```r
  exact_match <- ifelse((sg1 == sg1a | sg2 == sg1a), 1, 0)
  one_match <- exact_match
  cov1_match <- exact_match
  ```
  After:
  ```r
  contains_match <- ifelse((sg1 == sg1a | sg2 == sg1a), 1, 0)
  exact_match <- ifelse((sg1 == sg1a & is.na(sg2)) | (sg2 == sg1a & is.na(sg1)), 1, 0)
  one_match <- contains_match
  cov1_match <- contains_match
  ```
- **Red → green:** main session, repository at `1c2565b3` with only the test file: red 2 failures of 6 — `test-cv-exact-match.R:32:3` `actual`: 1.0 1.0 1.0 0.0 0.0 vs `expected`: 1.0 0.0 0.0 0.0 0.0, and `:68:3` 0.60 vs 0.20 — 3.0 s; green with `3612d0ba` code: 6 passes, 2.8 s.
- Implementer (2 fail/4 pass, 2.89 s → 6 pass, 2.78 s) and verifier (FAIL 2/PASS 4, 2.93 s → PASS 6, 2.90 s) reproduced it.
- **Expected values:** by hand from the spec's F4 definition, folds (A, NA), (A, B), (B, A), (B, NA), (NA, NA), full-data subgroup {A}.
  - Per fold: `exact_match` `c(1, 0, 0, 0, 0)`; `one_match` and `cov1_exact` `c(1, 1, 1, 0, 0)` — the spec's Test-line vectors. `find_metrics` Exact 1/5, At least 1 3/5, Cov1 exact 3/5 are their fold means.
  - Every value carries a comment (:28-31, :34-35, :65-67); the verifier flagged none missing and re-derived all six by hand.
- **Unaffected case:** two-factor and no full-data subgroup (`evidence/F4/unaffected.R`): `CV_sgs()` on 16 (sg1, sg2) pairs with `sg_analysis` `c(A, B)`, `c(B, A)`, `NULL`, `character(0)`; two-factor (N = 400) and null (N = 160) survival fits through `forestsearch_Kfold`, `forestsearch_KfoldOut` (outall F/T), `forestsearch_tenfold`. All 14 objects `identical: TRUE` (six timing leaves dropped). Walls: implementer 13.00 s / 12.73 s; verifier 12.89 s / 12.90 s.
- **Existing tests:** 8 CV-related files plus `test-mr-inference.R` (one block per process) ran with 0 failures and no errors for implementer and verifier; expectations updated: none.
- **What a user sees:** one-factor fit (`sg.harm: {age <= 55}`, N = 300, seed 31), `forestsearch_KfoldOut(details = TRUE)`: `Exact match: 0.4` → `Exact match: 0.2`. Fold 1 (`{age <= 55} & !{sex}`) no longer counts as exact; "At least 1 match" and "Cov 1 exact" stay 0.4.
- **Verification:** CONFIRMED; no repair round; committed diff equals `patches/F4.patch` apart from index and hunk-header lines. Three minor verifier problems: the `forestsearch-2024.Rmd` `kfold-cv` chunk may show a lower "Exact Match" at next render, and the roxygen `@return`/"Exact match on all factors" wording is unchanged — both left as open items (`vignettes/` read-only; NEWS states the new meaning, `88302212:NEWS.md:317-322`); the (NA, A) clause has no test — optional, not added.
- **Departures:**
  - "The other slot is NA" read as either slot, so (NA, A) is exact — the audit's proposed form; package folds never produce (NA, A).
  - A fold (A, A) now scores 0 (was 1), by the same wording; agrees with canonical equality; package folds do not produce it.
  - Spec names `cov1_match`; the returned field is `cov1_exact`, which the test checks, as the spec's Test line does.
  - Extra test block via exported `forestsearch_KfoldOut()`; its values follow arithmetically from the spec's vectors.
  - `test-mr-inference.R` run one block per process, because the whole file exceeds the 120 s per-run cap.
  - Spec's site pointer "near 1818-1826" (at `28e71ffa`); at the start SHA the block is :1818-1827.
- Full evidence: detail file, section F4.

### F5 — A-151: test_hte_from_forestsearch() includes continuous covariates

- **Status:** PASS (after one repair and a second verification) — commit `4aaab2f0` — `R/test_hte_crump.R` +22 / -2; `tests/testthat/test-hte-crump-confounders.R` new, +191 (2 files, 213 insertions, 2 deletions). No roxygen, `man/`, `NAMESPACE`, `DESCRIPTION` or `NEWS.md` change.
- **Site, before → after:** `R/test_hte_crump.R:397-400` and `:416` at the start SHA `76a62a71` → `:397-416` and `:432-436` at head `61c15088`. Site 1: each `confounders.evaluated` label is mapped to `df` columns (a label that is exactly a column name is that column; any other label is parsed), an unresolved label warns, unique columns are used. Site 2: an empty `X_cols` now warns before returning NULL (it was silent).
```r
# Before, R/test_hte_crump.R:397-400 (cut labels such as "age <= 55" are not columns, skipped silently at :411)
    covs <- fs_result$confounders.evaluated
```
```r
# After, R/test_hte_crump.R:397-416 (excerpt)
      lab_vars <- if (lab %in% names(df)) lab else
        intersect(tryCatch(all.vars(parse(text = lab)),
                           error = function(e) character(0)),
                  names(df))
      if (length(lab_vars) == 0)
        warning("Label '", lab, "' in confounders.evaluated cannot be ",
```
- **Red → green:** main session, repository at `3612d0ba` plus only the test file: 19 of 35 failed (all printed reports `Failure`, at :71, :73, :81, :98, ...; e.g. :71 `` `actual`: "bm" `` vs `` `expected`: "age" "bmi" "bm" `` as in the implementer's red), 5.9 s; after the code was applied, 35 of 35 pass, 5.4 s.
- Implementer reproduced it (red 19 of 35 failed, 0 errors; green 35 of 35, also on a fresh export), and verifier round 2 too (pre 19 failures, 0 errors; post 35 of 35).
- **Expected values:**
  - Derived by hand from get_fsdata() docs and code (a confounder with 4 or more distinct values gets `var <= value` cut labels, a binary one keeps its bare name; cuts first, then categorical), from the spec text, and from the code path for warning counts. None copied from the fixed code's output; both verifiers re-derived every value and agreed.
  - Every value carries a comment except block 2 :97-99 (0 warnings, non-NULL, class `hte_test`), flagged by verifier round 2 (minor) and left open.
- **Unaffected case:** spec: "a fit whose evaluated labels are all column names gives identical output." `evidence/F5/unaffected.R` ran `forestsearch()` (no LASSO or GRF) then `test_hte_from_forestsearch()` (OLS and logistic) per tree, for s1 (syntactic names), s2 (`"bm status"`, `"ER-pos"`) and s2x (s2 on a df that also has `ER`, `pos`). `start` vs `wt/F5`: all `identical()` TRUE (4.61 s / 4.48 s / compare 0.18 s). Verifier round 2: all TRUE (4.32 s / 4.46 s), own check 24 of 24 TRUE (2.51 s / 2.43 s). Cumulative at head `61c15088`: all TRUE (4.37 s / 4.44 s).
- **Existing tests:** none exist for this function (grep of `tests/testthat` finds nothing), so none were run; expectations updated: none.
- **What a user sees** (`evidence/F5/example.txt`; labels `age <= 55 | ... | bmi <= 30 | bm`; `is.null(result)` is the fit with only continuous confounders):
```
# Before (start)
covariates_used: bm
is.null(result): TRUE
# After (F5)
covariates_used: age bmi bm
is.null(result): FALSE
```
  New warnings where there was silence: `Label 'zzz <= 1' in confounders.evaluated cannot be resolved to a column of df; skipped.` and `No covariate found in df for the ForestSearch covariates; returning NULL.`
- **Verification:** round 1 PROBLEMS, one repair, round 2 CONFIRMED.
  - Round 1 blocking: the first version parsed every label, so `"bm status"` was dropped and `"ER-pos"` became `ER`, `pos`, breaking the unaffected case. Fixed: exact column-name labels kept (head :406), test block 5 added, unaffected members s2, s2x added.
  - Round 1 minor: block 4 warning checks too loose. Fixed (:140, :150, :151). Round 1 minor (harness recipe only, `git apply --index` on an unstaged export): no patch change.
  - Round 2 minor: no comment at block 2 :97-99, left open; notes cite the deleted round-1 `'my age'` probe, noted, an equivalent probe survives at `verify/F5/probe_nonsyn_cont.R`.
- **Departures:**
  - Column-name labels are resolved before parsing: the spec's literal `all.vars()` method breaks its own unaffected case for non-syntactic names. Verifier round 2 judged this acceptable.
  - Mapping only on the `confounders.evaluated` branch, as the spec names only that field; the `confounders.candidate` fallback shares the new no-covariate warning, which the spec states generally.
  - Warning wording is the implementer's (the spec gives none); one warning per unresolved label, same message for a parse failure and a non-column.
  - Three extra test blocks (3-5) for the spec's two warning clauses and the column-name unaffected member; block 1 tests order exactly, with continuous confounders passed first on purpose.
  - Roxygen not changed: the spec does not ask for it.
- Full evidence: detail file, section F5.

### F6 — A-042: LASSO, the intercept is not a selected variable

- **Status:** PASS — commit `f5b65516` (parent `4aaab2f0`) — `R/get_FSdata_helpers.R` (+15/-5), `R/forestsearch_main.R` (+8/-1), `man/lasso_selection.Rd` (+8/-2), new `tests/testthat/test-lasso-intercept.R` (+214); 245 insertions, 8 deletions. `get_FSdata()` is not edited.
- **Site, before → after:** `R/get_FSdata_helpers.R:333-338` at the start SHA `76a62a71` → `339-348` at head `61c15088`: `lasso_selection()` counts only confounder rows as selected/omitted. The roxygen title (:214) and description (:216 → :216-219) drop the Cox-only claim, `@return` gains :231-233, and `man/lasso_selection.Rd` matches. The summary label at `R/forestsearch_main.R:3394-3397` → `3393-3403` now uses the outcome-type switch of `R/get_fsdata.R:343-348`.
  ```r
  # Before
  selected <- names(coefs)[coefs != 0]
  omitted <- names(coefs)[coefs == 0]
  ```
  ```r
  # After
  conf_coefs <- coefs[names(coefs) %in% confounders.name]
  selected <- names(conf_coefs)[conf_coefs != 0]
  omitted <- names(conf_coefs)[conf_coefs == 0]
  ```
- **Red → green:** main session, `devtools::test(filter = "lasso-intercept")`: red at `4aaab2f0` + only the `tests/` hunks, 37 pass / 14 not (the same progress line as the implementer's; the body lists only the first 10), key message ``Expected `"(Intercept)" %in% res$selected` to be FALSE.``, 5.6 s; green with the whole patch, 62 pass / 0 failures, 5.8 s.
- Reproduced by the implementer (red 37 pass / 13 failures / 1 error, 5.78 s; green 62 pass, 5.97 s) and the verifier (the same red progress line, `failed=13 errors=1`, 5.83 s; green `expectations=62 failed=0 errors=0`, 5.93 s).
- **Expected values:** derived by hand from `?glmnet` (the intercept default and the `coef()` row order), the F6 spec (the intercept is not a confounder), the roxygen, the definition of selection at `lambda.min`, `get_FSdata()`/`cut_var()` default cut construction, and the label switch in `R/get_fsdata.R`.
  - The verifier said every expected value carries a comment and none was copied from the fixed output; it flagged none. The one-summary-line check (:208) has no comment of its own; the header at :181-185 covers it.
  - :65, :71, all of test 2 and the survival case of test 4 pass on both trees (invariant locks, not red evidence).
- **Unaffected case:** spec: the survival path, and a non-survival fit where LASSO selects at least one covariate, give identical candidate cuts. `evidence/F6/unaffected.R` ran `lasso_selection()`, `get_FSdata(use_lasso = TRUE)` and `forestsearch(use_lasso = TRUE, use_grf = FALSE)` for all four outcome types (survival also `use_grf = TRUE`): 952 leaves per side, 6 differ, all the intended `"(Intercept)"` removal from `selected`/`lassokeep`. `identical()` is TRUE for every survival leaf, and for binary/continuous/count for the candidate cuts (`df`, `confs_names`, `confs`), the whole `forestsearch` object without timing leaves, and `coefficients`; `ALL GATES: TRUE`. Wall 6.31 s / 6.35 s, compare 0.95 s. The verifier's re-run gave the same (6.22 s / 6.30 s / 1.01 s). Its own check (n = 400, four `get_FSdata()` configurations plus `survival::gbsg`) gave `ALL: TRUE` (5.31 s / 5.16 s).
- **Existing tests:** none updated (the commit touches only the new test file). The 7 related files, plus the verifier's `test-plan-restored.R`, have 0 failures and 0 errors; the pre-fix counts are the same.
- **What a user sees:** a binary outcome with no covariate signal, `forestsearch(outcome_type = "binary", use_lasso = TRUE, use_grf = FALSE, details = TRUE)` (`evidence/F6/example.txt`); the first line of each is the value of `lasso_selection(...)$selected`. When LASSO keeps a covariate, only the label and count change: `Cox-LASSO selected: 4 of 5 candidate factors` → `Linear-LASSO selected: 3 of 4 candidate factors`.
  ```
  # Before
  [1] "(Intercept)"
  Error in get_FSdata: get_FSdata: no candidate cut expressions could be constructed from the supplied inputs. ...
  # After
  character(0)
  Logistic-LASSO selected: 0 of 4 candidate factors
  > number of candidate factors evaluated: 8
  ```
- **Verification:** CONFIRMED; no repair round. `tools::checkRd` parses, roxygen regenerates the `.Rd` byte-identically, and the committed diff equals the gated patch. The verifier's 3 minor problems: (1) the spec's mocked test 2 cannot go red: kept as a lock, and the red comes from test 3; (2) :65, :71 and the survival case of test 4 pass on both trees: noted, no change; (3) the docs go beyond the literal spec: accepted as Departure 3. The suggested comments for (1) and (2) were not added.
- **Departures:**
  - The spec's mocked test passes on the start tree: the defect is in `lasso_selection()`, not `get_FSdata()`'s guard. Test 2 is kept as a lock, and the added test 3 (real `lasso_selection()`, only `coef()` mocked) gives the red.
  - Added test 4 for the printed label, because the label is in the spec's change list.
  - The roxygen description and `@return` were corrected as well as the title: the description repeated the Cox-only claim, and `@return` now states the changed contract.
  - "The existing rule in `get_FSdata()` then applies" was read as no edit to `get_FSdata()`.
  - "Identical candidate cuts" was read as identical `df`, `confs_names`, `confs` and the whole `forestsearch` object. `selected`/`lassokeep` may lose `"(Intercept)"`, which is the fix.
  - The label mirrors the existing switch exactly, so a binary outcome with an offset (a Poisson fit) still prints `Logistic-LASSO` (open item 14). The spec said to mirror that switch.
- Full evidence: detail file, section F6.

### F7 — A-107: fpr_calibration() outcome-type fallback

- **Status:** PASS — commit `65831b54` — `R/fpr_calibration.R` (+14/-4), `man/fpr_calibration.Rd` (+3/-1), `tests/testthat/test-fpr-calibration-fallback.R` (new, +128).
- **Site, before → after:** `R/fpr_calibration.R:236-237` and `:245-247` at the start SHA `76a62a71` → `238-248` and `256-257` at head `61c15088`. The no-type fallback changes from `"binary"` to `"survival"` (forestsearch()'s default), an unrecognised type now stops, and the switch's binary fallback arm is deleted. Roxygen `:91-92` → `:91-94` and `man/fpr_calibration.Rd:80` → `:80-82` gain one sentence stating the fallback.

Before:
```r
  otype <- outcome_type_approx %||% fs_params$outcome_type %||% "binary"
      d_eff_binary(n_min, p_event)  # fallback
```
After:
```r
  otype <- outcome_type_approx %||% fs_params$outcome_type %||% "survival"
  otype_ok <- c("survival", "binary", "continuous", "count")
  if (!is.character(otype) || length(otype) != 1L || is.na(otype) ||
      !(otype %in% otype_ok))
    stop(sprintf(paste0(
      "Unrecognised outcome type %s (from outcome_type_approx or ",
```

- **Red → green:** main session, repository at `f5b65516` (parent) + test hunks only: red 11 failures, 5 passes of 16 expectations, first `` Expected `res$call_args$otype` to be identical to "survival". ``, 3.8 s; green with the rest of the patch: 16 passes, 0 failures, 3.2 s.
- Implementer (red `[ FAIL 11 | WARN 0 | SKIP 0 | PASS 5 ]`, 3.87 s; green `PASS 16`, 3.42 s) and verifier (red 11 failed, 0 errors, 5 passed, 3.97 s; green 16 passed, 3.40 s) reproduced it.
- **Expected values:**
  - By hand: `"survival"` from forestsearch()'s `outcome_type` formals and `match.arg()` (`R/forestsearch_main.R:1746-1747`, `:2215`); d_eff 45, 60, 15 and 80/3 from `?d_eff_survival`, `?d_eff_binary` and `?d_eff_continuous` (e.g. 60 x (1 - 0.25) = 45); the message checks from the documented type values. None was taken from the fixed code's output.
  - Every value has a derivation comment; the verifier re-derived every value and flagged none.
- **Unaffected case:** spec: valid type → identical output; no type → identical simulated FPR. Implementer (`evidence/F7/unaffected.R`): start and fixed trees in separate `Rscript --vanilla` processes, real forestsearch() replicates (N = 240, `nsims = 20`, `set.seed(777)`). Four valid-type cases (object, print, warnings, `.Random.seed`) and the no-type case (FPR fields, warnings, RNG) all identical() TRUE: `UNAFFECTED CASE RESULT: IDENTICAL`. Walls 22.35 s start, 22.07 s fixed, 0.18 s compare. The verifier's re-run (21.86 / 21.65 / 0.19 s) and its own script (continuous, count, DGM mode; 11.31 / 10.64 / 0.18 s) also gave identical() TRUE.
- **Existing tests:** `tests/testthat/test-plan-restored.R` (block at line 153 now takes the survival path): implementer `[ FAIL 0 | WARN 0 | SKIP 0 | PASS 20 ]`, 20.15 s; verifier 20 passed, 20.09 s. Expectations updated: none.
- **What a user sees:** `print()` for a survival analysis with no type given (`evidence/F7/unaffected.txt`, case `none`); the Sim FPR lines are unchanged.
```
before:   P1 (per-subgroup): 23.1%  (d_eff = 10.0, uncorrected)
          Implied L_eff:     3.035
after:    P1 (per-subgroup): 20.1%  (d_eff = 28.0, uncorrected)
          Implied L_eff:     3.567
```
- **Verification:** CONFIRMED; red, green and the unaffected case reproduced. No repair round (the commit equals `patches/F7.patch`). Three minor problems:
  - Exact matching rejects `"surv"`, which forestsearch() accepts through `match.arg()`: left as specified, open item.
  - The notes' `Error:` prefix on the misspelt-type example came from the example script's `cat()`: noted, no code change.
  - `fpr_calibration_transfer()` drops the nested error, so a misspelt type now gives NA silently: open item, already disclosed by the implementer.
- **Departures:**
  - Validation is unconditional: the spec's stop has no condition, so the type is checked before the replicates even when `d_eff_approx` is supplied.
  - Exact matching (`%in%`), not `match.arg()`: the documented values are the four exact strings.
  - The switch's binary fallback arm was deleted: it cannot be reached after the up-front check.
  - The test stubs the `%dofuture%` loop with `local_mocked_bindings()` (as `test-plan-restored.R` does): type and `d_eff` are resolved before the loop; the unaffected case runs real replicates.
  - The misspelt-type checks were split into two `expect_error()` calls before green after a cascade error in the first red run; no expected value changed.
- Full evidence: detail file, section F7.

### F8 — A-047: duplicate-cut test compares names exactly

- **Status:** PASS — commit `ea5f82d9` — `R/get_fsdata.R` (+17/-11); new `tests/testthat/test-duplicate-cut-names.R` (+194).
- **Site, before → after:** `R/get_fsdata.R`, three sites in `get_FSdata()`, start `76a62a71` → head `61c15088`. Each now takes a cut's name as the text before `<=`/`>=`, trimmed, and compares it exactly.
  - Median removal (`replace_med_grf`): :448-453 → :448-452; `any(grepl(x, grf_cuts))` becomes `conf.cont_medians %in% grf_cut_names`.
  - GRF merge: :498-499 → :497-502 (quoted below).
  - DINA merge: comment :515-516, code :520 → comment :518-521, code :525-526; same two-clause test as the GRF merge.
```r
# Before, 76a62a71 :497-500 (GRF merge)
flag_omit <- sapply(grf_cuts, function(cut) any(sapply(confs, function(x) grepl(x, cut))))
```
```r
# After, 61c15088 :496-503 (GRF merge)
grf_cut_names <- trimws(sub("\\s*[<>]=.*$", "", grf_cuts))
flag_omit <- (grf_cut_names %in% conf.categorical) | (grf_cuts %in% confs)
```
- **Red → green:** main session, repository at `65831b54` + test file: red 8 failures, 0 errors of 13 (e.g. `:105:3` `` `actual`: 0 `` / `` `expected`: 1 ``; `:142:3` `Absent: "x <= 25.3"`), 3.3 s; green 13 pass, 3.1 s.
  - Implementer (8 failures, 3.25 s; 13 pass, 3.16 s) and verifier (8 failures, 3.34 s; 13 pass, 3.36 s) reproduced it, same lines and messages.
- **Expected values:**
  - All 13 derived by hand on three datasets built without random numbers, from `?is.continuous` (at least 4 distinct values is continuous) and the `cut_type` rules ("median": `round(median, 2)`; "default": mean, median, Q1, Q3, type-7).
  - The verifier re-derived all 13 before looking at the fixed output. Every value has a derivation comment; none is flagged.
- **Unaffected case:** non-overlapping confounder names (n = 300, `set.seed(20261112)`), 7 `get_FSdata()` configs + 2 `forestsearch()` runs. Start vs fixed: `OVERALL: IDENTICAL`, 9 of 9 (timing leaves dropped), 5.98 s / 5.90 s. Plain `identical()`: TRUE for c1-c7; FALSE for f1, f2 only in `find.grps$time_search` and `minutes_all`. Cumulative (all 14 fixes): also IDENTICAL, 6.11 s / 5.98 s.
- **Existing tests:** implementer 26 of 27 files and verifier 28 of 29 files completed with 0 failures, 0 errors; `mr-inference` hit the 120 s cap (it cannot reach the changed lines: `use_grf = FALSE`). Expectations updated: none.
- **What a user sees:** `confounders.name = c("bm", "bmi", "age")`, `grf_cuts = "bmi <= 25.3"`, `dina_cuts = "bmi >= 25.3"`, `bm` binary. Before, both `bmi` cuts were dropped only because "bm" is a prefix of "bmi":
```
== start tree
[1] "bmi <= 26" "bmi <= 23" "bmi <= 29" "age <= 60" "age <= 50" "age <= 70"
[7] "bm"
== fixed tree
[1] "bmi <= 25.3" "bmi >= 25.3" "bmi <= 26"   "bmi <= 23"   "bmi <= 29"
[6] "age <= 60"   "age <= 50"   "age <= 70"   "bm"
```
- **Verification:** CONFIRMED; no repair round recorded. Mutants M1 and M3 fail the new tests; M2 passes all 13.
  - Problem 1 (minor): the exact-duplicate guard (test :162-178) does not pin the `%in% confs` clause, because `unique()` already removes exact duplicates. Noted, not fixed: the gated patch was committed unchanged, and the caveat is stated as the verifier asked.
  - Problem 2 (minor): the unaffected evidence prints `confs=0` for the `forestsearch()` rows (the pool is `fit$confounders.evaluated`, 11 entries). Noted, not fixed: an evidence label only.
- **Departures:**
  - Names are parsed from `<=`/`>=` only, as the spec words it (the audit's sketch used `[<>]=?`); GRF and DINA write only these operators.
  - "The current pool" was read as `confs` at each merge, and "a categorical confounder" as `conf.categorical`, the vector that built the pool.
  - The DINA comment now says `"var <= value" or "var >= value"`, because the new parse depends on it.
  - Tests beyond the spec's list: a median-removal block (:181-194), the `x <= 25.3` vs `x <= 2` case at both merges, and the categorical guard at both.
- Full evidence: detail file, section F8.

### F9 — A-125: MRCT default confounders exclude internal spline columns

- **Status:** PASS — commit `a30d6e7c` (parent `ea5f82d9`) — `R/mrct_case_study.R` +18/-2, `R/mrct_simulation.R` +20/-2, `man/mrct_case_study.Rd` +6/-1, `man/mrct_region_sims.Rd` +9/-1, `tests/testthat/test-mrct-default-confounders.R` (new) +83/-0. Commit diff equals the gated `patches/F9.patch`.
- **Site, before → after:** start `76a62a71` (F9 R files identical at `ea5f82d9`) → head `61c15088` (files unchanged since `a30d6e7c`).
  - `mrct_region_sims()` NULL branch: `R/mrct_simulation.R:321-324` → `329-342` (new 329-338, changed 341). Quoted below.
  - `mrct_case_study()` defaults list: `R/mrct_case_study.R:291-293` → `291-300` (same `spline_helpers` block) and `306-309` (`spline_helpers` added to the excluded names).
  - Roxygen: `R/mrct_simulation.R:52` → `52-60`; `R/mrct_case_study.R:172` → `172-177` (`@param fs_args`); Rd regenerated, `tools::checkRd()` 0 problems.

Before (`R/mrct_simulation.R:321-324`):
```r
    confounders.name <- setdiff(
      all_vars[grepl("^z_|^ecog|^strat", all_vars)],
      exclude_vars
    )
```
After (`R/mrct_simulation.R:333-341`, excerpt):
```r
    spline_info <- dgm$model_params$spline_info
    spline_helpers <- if (!is.null(spline_info)) {
      paste0(spline_info$var, c("_treat", "_k", "_k_treat"))
    ...
      c(exclude_vars, spline_helpers)
```
- **Red → green:** main session, repo at `ea5f82d9` + test hunk only: 2 failures of 8, "Expected `rs_set` to be identical to `expected`" at `:65:3` (extra `z_bm_treat`, `z_bm_k`, `z_bm_k_treat`), same for `cs_set` at `:76:3`, 6.9 s; green after the rest of the patch: 8 of 8 pass, 6.7 s.
  - Implementer (2 fail / 6 pass → 8 pass; 7.58 / 7.18 s) and verifier (`FAIL 2 | PASS 6` → `PASS 8`; 7.55 / 7.22 s) reproduced it.
- **Expected values:**
  - By hand from source: the six-name list (:50) from `create_dgm_for_mrct()` (`continuous_vars`, `factor_vars`) and the column order of `process_continuous_vars()` / `process_factor_vars()`; helper names from `create_spline_variables()`; `z_prior_treat` and the region indicator from the spec. Verifier re-derived all from source at `76a62a71`; none copied from fixed output.
  - Comments: all but `expect_length(line, 1L)` (:63, one replicate) carry one; verifier also flagged (minor) that :40-49 does not say why other `df_super` columns are absent.
- **Unaffected case:** spec "with `confounders.name` given, output is identical": both functions with the six and with the nine names, the printed line, plus 3 no-spline extras; `identical()` TRUE for every case (`ALL IDENTICAL: TRUE`) for implementer and verifier. Walls start / fixed: 15.89 / 16.20 s and 15.75 / 15.96 s. Not run by the main session.
- **Existing tests:** 7 files using the two functions, all FAIL 0 for implementer and verifier (verifier's exit 1 came from its own summary `cat()`); expectations updated: none.
- **What a user sees:** `mrct_region_sims(..., details = TRUE)`, `confounders.name` omitted; `mrct_case_study()`'s `fs_args$confounders.name` also goes from 9 to 6 names; `Hhat` `!{z_male}` on both.
```
Before:    Confounders:  z_age, z_bm, z_male, z_histology, z_prior_treat, z_regA, z_bm_treat, z_bm_k, z_bm_k_treat
After:     Confounders:  z_age, z_bm, z_male, z_histology, z_prior_treat, z_regA
```
- **Verification:** CONFIRMED; no repair round. Two minor problems not applied, kept as open items: :82 (`cs$Hhat` vs `rs$sg_found`) passes on both trees; missing comments at :40-49 and :63. No main-session comment added.
- **Departures:** none in the code. Interpretations: spline info read as `!is.null(dgm$model_params$spline_info)`, exact names (no suffix match); exclusion stated in `@param fs_args` because `mrct_case_study()` has no `@param confounders.name`. Additions (identifier check :82, no-spline unaffected case) do not change the scope.
- Full evidence: detail file, section F9.

### F10 — A-064: count outcomes with GRF are oriented

- **Status:** PASS — commit `f2e8de47` — `R/grf_subg_harm_glm.R` (+15 / -8); new `tests/testthat/test-grf-count-orientation.R` (138 lines); 2 files, +153 / -8.
- **Site, before → after:** `R/grf_subg_harm_glm.R` :448-452 at start `76a62a71` → :450-458 at head `61c15088` (new count branch :452-457, after the count transform, before `causal_forest()`). Header comment :437-441 → :437-443 adds count to the flip list; frontier comment :515-517 → :521-524 says the flip covers count too.

```r
    # Before (start :448-452)
    }
    # count: log transform already applied above; flip not needed
```
```r
    # After (head :452-457)
    } else if (outcome_type == "count") {
      # count: the log(Y + 0.5) or identity transform above keeps "higher
      # count = higher Y_grf", so negate it: a higher count is then worse,
      # as for binary and continuous.
      Y_grf <- -Y_grf
    }
```
- **Red → green:** Main, repository at `a30d6e7c` with only the test file: red, 4 failures, 0 errors of 14 expectations at test :110, :114, :134, :136 ("Actual comparison: 0.00 < 0.80"; "Expected `identical(res_adv$data$treat.recommend, res_ben$data$treat.recommend)` to be FALSE"), wall 4.3 s; green 14 pass, 0 failures, 4.2 s.
- The implementer (red 10 pass / 4 failures, 4.60 s; green 14 pass, 4.29 s) and the verifier (same 4 failures, 4.60 s; 14 pass, 4.54 s) reproduced it.
- **Expected values:**
  - By hand from the roxygen: members are `treat.recommend == 0L` (`@return`, :142-144); `adverse_outcome` (:100-103) puts harm in `z == 1` for TRUE and `z == 0` for FALSE, so the two memberships differ.
  - The 0.8 share margin is worked out by formula in the test header (:28-53; about 7 SE); the verifier re-derived every number independently (`verify/F10/handcalc.txt`).
  - Every value carries a comment except the existence assertions (:103-106, :129-132), which the verifier flagged as lacking an inline comment (minor).
- **Unaffected case:** `adverse_outcome = FALSE` for counts, plus binary and continuous outcomes. 12 cases (direct calls and 3 `forestsearch()` cases), start vs fixed tree in separate `Rscript --vanilla` processes: all `identical = TRUE` after dropping timing leaves; walls start 14.94 s, fixed 14.81 s, compare 5.01 s. Verifier: same 12 identical (15.20 / 15.26 / 5.05 s); its own 7 new cases identical, positive control `identical=FALSE`.
- **Existing tests:** 5 files reach count + `adverse_outcome = TRUE`; 16 files pass on the fixed tree (3 baseline-skipped files re-run in full, pass); verifier 30 of 31 pass, `mr-inference` hit the 120 s cap (cannot reach the function). Expectations updated: none.
- **What a user sees:** direct verbose call on the test data, where treatment raises the count in `z == 1`. Before (`76a62a71`, first line) the "harm" subgroup was the benefit region; after (fixed tree, last two lines) it is where treatment raises the adverse count.
```
n(H) = 497 | share of H in z == 1 = 0.00 | IRR(H) = 0.290 (GLM)
[grf.subg.harm.glm] Subgroup found: n(H) = 500 (50.0%) | definition: z > 0
n(H) = 500 | share of H in z == 1 = 1.00 | IRR(H) = 2.815 (GLM)
```
- **Verification:** CONFIRMED; no repair round. Three minor problems: existence assertions lack inline comments (not changed, open item); NEWS hand-off omitted the GRF arm of `run_simulation_analysis()` (fixed in NEWS commit `88302212`, `NEWS.md` :361-365); CATE footnote "(risk-difference scale)" in `R/summary_utility_functions.R` :913-918 is wrong for counts (pre-existing, open item).
- **Departures:**
  - Header comment reworded to "a higher Y_grf is better and positive -CATE = ...": the old wording contradicted the next sentence (comment only; verifier judged it correct).
  - "Negate the transformed outcome" read as negating whichever count transform applies (log or identity), placed after the transform.
  - "Mostly" read as a share of at least 0.8; "memberships differ" read as `!identical(treat.recommend)`.
  - Test runs two paths (tree + log, frontier + identity) and also checks FALSE, so "differ" compares two correctly oriented answers.
  - `skip_on_cran()` at the top: the file fits four forests, as `test-grf-inclusion-band.R` does.
  - No verbose message for the count flip and no roxygen change: the spec did not ask for them.
- Full evidence: detail file, section F10.

### F11 — DINA orientation sign on the screening route

- **Status:** PASS — commit `61c15088` (start `76a62a71`) — `R/forestsearch_main.R` (3 +-), `tests/testthat/test-dina-screening-tau-sign.R` (new, +132); 2 files, +134/−1. R/ and the test file are unchanged to final head `88302212`.
- **Site, before → after:** `R/forestsearch_main.R:3292-3303` at `76a62a71` → `:3290-3302` at head `61c15088` (added line :3301). The screening call inside `if (use_dina && is.null(dina_cuts))` / `if (isTRUE(da$selected_only))` now passes the sign the `dina` route already passes (`R/forestsearch_helpers.R:1542`, `:1604`).

Before (from :3299-3303 at `76a62a71`):
```r
          effect_neighborhood = effect_neighborhood
        )
```
After (from :3298-3303 at `61c15088`):
```r
          effect_neighborhood = effect_neighborhood,
          tau_sign            = .dina_tau_sign(outcome_type, adverse_outcome)
        )
```

- **Red → green:** Main session, repository at `a4cce90a` plus test hunks only: red 13 pass, 3 failures of 16, 0 errors (:96:5 twice, "Expected \`run$signs\` to be identical to \`cs$expected\`", actual `1.0`, expected `-1.0`; :125:3, `x1 >= -0.481` does not match `"^x1 <= "`), 6.8 s; green after the code hunk 16 pass, 6.7 s.
- Implementer (red 13/3, 7.17 s; green 16/0, 6.92 s) and verifier (red 13/3, 7.22 s; green 16/0, 6.81 s) reproduced it, at :95/:122 in the pre-comment test file.
- **Expected values:**
  - Signs derived by hand from `?forestsearch` (`adverse_outcome` NULL resolution), `?dina_subgroup` (`tau_sign`: -1 for continuous and binary with `adverse_outcome = FALSE`) and the `@noRd` roxygen of `.dina_tau_sign()` (`R/forestsearch_helpers.R:1406-1407`: survival, count and `TRUE` keep 1).
  - Test 2 cut direction by construction (:107: y shifts −2 where `x1 <= 0`, +2 where `x1 > 0`); the `adverse_outcome = TRUE` mirror is the control. Inputs fixed before the first run; the verifier found no value came from running the fixed code.
  - Every value carries a comment at commit. The verifier flagged both `expect_length(..., 1L)` in test 2 (:124, :130) as lacking one; the comment at :122-123 was added and :130 shares it.
- **Unaffected case:** Spec names survival, count, binary `adverse_outcome = TRUE` and every `subgroup_method = "dina"` fit. `evidence/F11/unaffected.R` ran U1-U4, D1-D9 and extras X1-X3: `identical()` TRUE in all 16 cases for implementer and verifier. Walls: implementer 10.24 s start / 10.21 s fixed (compare 1.65 s); verifier in-R 9.73 s / 9.30 s.
- **Existing tests:** 0 failures in every run (implementer: 8 DINA files + 39 broad `forestsearch()` files, `test-mr-inference.R` block by block after the 120 s cap; verifier: 8 + 12 + 9 files); expectations updated: none.
- **What a user sees:** continuous outcome, default orientation, `use_dina = TRUE`, `max_depth = 1L` (`evidence/F11/example.txt`; verifier reproduced exactly):
```
# start (76a62a71)
fs$dina_cuts: x1 >= -0.481
fs$sg.harm:   {x1 <= 0}
# fixed
fs$dina_cuts: x1 <= 0.352
fs$sg.harm:   {x1 <= 0.352}
```
- **Verification:** CONFIRMED; no repair round recorded. Two minor problems, both resolved by the main session before commit as comment-only test edits: the "keeps 1" source now cites both `?dina_subgroup` and the `.dina_tau_sign()` roxygen (:60-63); the `expect_length` comment was added (:122-123), as the verifier suggested. Code hunk identical to the verified patch; red/green re-run on the repository afterward.
- **Departures:**
  - Code: none; the one argument the spec names. No roxygen change, so `document()` was not run.
  - Site found by content: the spec's "near 3254–3265" was taken at `28e71ffa`.
  - Test 1 reads the effective sign from a `local_mocked_bindings` wrapper, so the four sign-1 cases are controls that pass before and after.
  - Test 2 "points at that region" read as a single x1 cut checked by operator, not cut value; "every dina fit" read as D1-D9 (D7-D9 added after a first 12-case run).
- Full evidence: detail file, section F11.

### F12 — Three displays use the fit's own orientation

- **Status:** PASS — commit `25ecaa0f` (parent `f2e8de47`) — `R/mr_estimates_table.R` +8/-0, `R/plot_sg_glm_outcomes.R` +9/-0, `R/plot_subgroup_results_forestplot.R` +13/-0, `tests/testthat/test-display-orientation.R` (new) +227/-0. No roxygen, `man/`, `NEWS.md` or `DESCRIPTION` change.
- **Site, before → after:** each display now resolves `adverse_outcome` from the fit (fallback: `forestsearch()`'s rule, `R/forestsearch_main.R:1840-1841`) and passes it to `make_effect_estimator()`, whose default `TRUE` applied before.
  - `R/plot_sg_glm_outcomes.R:163-166`, call `:254-262` at `76a62a71` → block `:166-172`, call `:262-271` gains `:268` at head `61c15088` (quoted below).
  - `R/plot_subgroup_results_forestplot.R:309-312`, call `:432-440` → block `:312-322`, call `:444-453` gains `:450` (inside `create_glm_row()`).
  - `R/mr_estimates_table.R:203-208` (`boot_results = NULL` route) → block `:204-210`, call `:211-216` gains `:215`.

Before
```r
      effect_measure = effect_measure,
      adjust_covariates = adjust_covariates,
```
After
```r
  adverse_outcome <- fs.est$args_call_all$adverse_outcome
  if (is.null(adverse_outcome)) {
    adverse_outcome <- outcome_type %in% c("binary", "count")
  }
...
      adverse_outcome = adverse_outcome,
```

- **Red → green:** main session, repository at `f2e8de47` with the test file only: 24 of 93 not passed (69 passes; 10 listed, all "Failure", then "... and 14 more"), e.g. "Expected forest est H to equal `ref[[k]]$est`." actual -1.1, expected 1.1; 6.8 s. Green (tree committed as `25ecaa0f`): 93 passes, 6.4 s.
- Implementer reproduced it (24 failed, 0 errors of 93 → 0 failed; 7.28 s → 7.08 s), and so did the verifier (24 failed, 0 errors of 93 → 0 failed; 7.6 s → 6.68 s in R). Every failure is a sign mismatch.
- **Expected values:**
  - By hand: oriented MD `mean(y | T = 1) - mean(y | T = 0)`, negated when `adverse = FALSE`, with pooled SE. Orientation from `?forestsearch` `@param adverse_outcome` (`R/forestsearch_main.R:1104-1107`); z from the documented default `conf.level = 0.95`; MR limits use the 1.96 of `analyze_subgroup_glm()` (`R/summary_utility_functions.R:223`). The hand H value is tied to `max_sg_est` within `1e-10`.
  - The verifier rederived them in base R (H 1.0941, se 0.2913) and found no value copied from the fixed code's output.
  - Comments: the verifier flagged the forest-plot `1e-10` tolerance (`:127-130`) as unexplained in the test, and the `:29-30` comment on the `1e-4` MR-limit tolerance as muddled. The MR row locator `c("Qstnbl", "Recmnd")` (`:138`) has no comment at the line.
- **Unaffected case:** spec: survival fits and non-survival fits with `adverse_outcome = TRUE`. `evidence/F12/unaffected.R` ran 6 cases (survival, continuous TRUE, binary OR and RD, count, `forest_no_fit`) on start vs fixed plus a start control. `survival` and `forest_no_fit` `identical()` = TRUE; the four others differ only at ggplot2's labels `S7_class` attribute, as in the control; with it dropped all six are TRUE. Walls 10.26 / 10.55 / 10.21 s (implementer), 10.42 / 10.18 / 10.27 s (verifier, same output).
- **Existing tests:** the five files that call the displays all pass on the fixed tree (0 failed / 0 errors / 0 skipped); expectations updated: none.
- **What a user sees:** continuous fit at the default (`adverse_outcome = FALSE`, `max_sg_est = 1.09411`), `mr_estimates_table(fit)`, `MD (95% CI)`. Before, the table mixed orientations: `Qstnbl` `1.09 (0.50, 1.68)` but `Recmnd` `-0.10 (-0.34, 0.15)`. After, `Recmnd` reads `0.10 (-0.15, 0.34)`. The plot label for Questionable goes from `MD = -1.094 (-1.665, -0.523)` to `MD = 1.094 (0.523, 1.665)`. `Mean(C)`, `Mean(T)` and `Diff` stay raw, so `Diff = -MD` for this fit.
- **Verification:** CONFIRMED; no repair round; the commit equals the gated patch (`main/F12/patch_equal.txt` empty). Two minor problems, both left as open items with no code change: (1) `sg_tables()`, the CV summary tables and `compute_frontier_cis()`'s `%||% TRUE` fallback keep the default `TRUE` (outside F12's named sites); (2) the two tolerance comments above (not changed, to keep the commit equal to the gated patch; neither affects whether a tolerance separates the trees).
- **Departures:**
  - Forest plot without a fit (`fs.est = NULL`): the default `TRUE` is kept and passed explicitly (`:316` at head), since applying the rule would flip continuous rows drawn without a fit.
  - With a fit, the forest plot orients all GLM rows, not only H and Hc, because `create_glm_row()` builds every GLM row; the verifier judged this to match the spec's intent.
  - The fallback uses each display's resolved `outcome_type`, following `R/bootstrap_dofuture_main.R:306-309`.
  - The fit stores its own estimate only for H, so the test checks every set against a hand value and ties H to `max_sg_est`.
  - MR table H-row limits are not checked (documented as robust sandwich limits); narrowed before the first red run.
  - Unaffected case extended (binary and count with `adverse_outcome` unset, `forest_no_fit`); "identical" judged after dropping the `S7_class` attribute, which differs between two runs of the same tree.
  - No roxygen change: the spec names a code change only.
- Full evidence: detail file, section F12.

### F13 — A-055: risk-difference standard error, unadjusted case (change of method)

- **Status:** PASS — commit `89ad9c19` (parent `25ecaa0f`) — `R/glm_effect_estimators.R` (+12/-1), `tests/testthat/test-rd-standard-error.R` (new, +82); 2 files, +94/-1. No `man/` change (`.estimate_rd` is `@noRd`).
- **Site, before → after:** `R/glm_effect_estimators.R:534-536` (`.estimate_rd()`, tier 2, logistic g-computation) at the start SHA `76a62a71` → `:534-547` at head `61c15088` (new SE `:544`; old expression kept at `:546` for the adjusted case; the new explanatory comment is omitted below).
  Before:
  ```r
  se_rd <- abs(avg_deriv) * se_treat
  ```
  After:
  ```r
  se_rd <- if (n_adj == 0L) {
    sqrt(p1 * (1 - p1) / n1 + p0 * (1 - p0) / n0)
  } else {
    abs(avg_deriv) * se_treat
  }
  ```
- **Red → green:** main session on the repository (red at `25ecaa0f`, test file only): 4 failures of 14 expectations, e.g. `:50` "Expected `res$se` to equal `sqrt(0.3 * 0.7/40 + 1 * 0/40)`." (`actual: 0.000`, `expected: 0.072`), 2.8 s; green: 14 of 14 passed, 2.7 s.
  Implementer (14 expectations: 4 failed, 0 errors, 3.10 s → 0 failed, 2.89 s) and verifier (14: 4 failed, 0 errors, 2.77 s → 0 failed, 2.76 s) reproduced it.
- **Expected values:**
  - SEs by hand from the spec's two-sample binomial formula on 40-per-arm slices: p0 = 12/40, p1 = 40/40 → 0.0724569 (`:50`; rounded 0.07246, `:51`); p0 = 0/40, p1 = 32/40 → 0.0632456 (`:72`; rounded 0.06325, `:73`). Estimates, `converged`, `method_used` (`"logistic_margins"`), `measure`, `n0`/`n1` from the spec and the start-SHA source; estimate tolerance 1e-6 set before any run from glm's 1e-8 stopping rule.
  - Comments: the verifier flagged that the second block's `converged`, `method_used`, `measure`, `n0`/`n1` (`:78-81`) have no comment of their own (minor; not added). `measure` `"RD"` has no source comment in either block.
- **Unaffected case:** spec: OR, RR, and RD with adjustment covariates give identical output. `evidence/F13/unaffected.R`: 112 `make_effect_estimator()` results (14 configurations x 8 slices, 5 with a degenerate arm) plus two `forestsearch()` runs (OR default; RD with `"age"`): `OVERALL identical: TRUE`; walls 4.75 s start, 4.85 s fixed. Verifier re-run TRUE (4.82 s, 4.87 s); its own 324-result probe: the 288 spec-named results identical, the 12 that differ are all unadjusted RD (`se` only). Scope: RD with covariates is identical under `ps_adjust_method = "none"`; with `"dr_gcomp"`, or `"iptw"` without `sw`, the covariates are dropped and the new SE applies.
- **Existing tests:** 14 files, no failures; 13 identical start vs fixed, mr-inference capped at 120 s on both trees after the same 133 passing expectations. Expectations updated: none (estimability-boundary reaches the new branch, SE 1.5299e-05 → 0.21909, but asserts only a non-NA estimate and NULL `reason`).
- **What a user sees:** `make_effect_estimator("binary", effect_measure = "RD")`, 40 per arm, control 30%, treated 100%; the ITT row of `SG_tab_estimates_glm()` goes from `0.7000 (0.7000, 0.7000)` to `0.7000 (0.5580, 0.8420)`.
  ```
  Before: control  30%, treated 100%: estimate 0.7000  se 5.4089e-06  95% CI (0.7000, 0.7000)  converged TRUE  method logistic_margins
  After:  control  30%, treated 100%: estimate 0.7000  se 0.072457  95% CI (0.5580, 0.8420)  converged TRUE  method logistic_margins
  ```
- **Verification:** CONFIRMED (red, green and unaffected reproduced; expected values re-derived by hand). No repair round. Three minor problems: (1) the weights trace missed `adjust_covariates` with `dr_gcomp` or `iptw` without `sw` → row added to the report's weights table, no code change; (2) the kept comment at head `:531-532` ("conservative for subgroup estimation") contradicts the new one at `:540-541` → left open; (3) missing test comments at `:78-81` → not added (commit equals `F13.patch`; no `F13.verified.patch`).
- **Departures:**
  - "No adjustment covariates and no weights" read as `n_adj == 0L` in tier 2, because no weighted fit reaches tier 2 (IPTW with `sw` returns earlier; neither `glm()` there takes weights).
  - Covariates dropped before the call (`dr_gcomp`, or `iptw` without `sw`) get the new SE, read as within the spec; so the unaffected clause holds for RD with covariates under `"none"` only (plus `iptw` with `sw`). `forestsearch()` never builds the other combination.
  - No routing to the raw-proportion tier: `se_treat` and `avg_deriv` are still computed in both branches, so failures route as before; tests assert `"logistic_margins"` and `converged = TRUE`.
  - Extra test assertions beyond the spec's two values (`round(se, 5)`, `method_used`, `measure`, `n0`/`n1`), all from the spec or the start-SHA source.
  - mr-inference did not finish within the 2-minute cap on either tree; judged out of reach of `.estimate_rd` (binary OR only, 0 hits in the instrumented run).
- Full evidence: detail file, section F13.

### F14 — DINA floor: fit, record and replay agree in the forced-gaussian corner

- **Status:** PASS — commit `685e0c5d` — `R/forestsearch_main.R` (+7 −9), `tests/testthat/test-dina-gaussian-floor.R` (+77 −7); no roxygen, `man/`, `NEWS.md` or `DESCRIPTION` change.
- **Site, before → after:** `R/forestsearch_main.R:2820-2824` and `:3283-3291` at start `76a62a71` → `:2820-2823` and `:3282-3289` at head `61c15088`. Site 1 (`subgroup_method = "dina"`, last argument of `.forestsearch_dina_select(`) and site 2 (`use_dina = TRUE`, `m_diff_sel`) used the raw supplied `hr.threshold` as the gaussian floor; both now use the resolved `effect.threshold`. Site 2:

Before
```r
m_diff_sel <- if (identical(da$fit$family, "gaussian")) {
                if (user_set_threshold) hr.threshold else effect.threshold
```
After
```r
m_diff_sel <- if (identical(da$fit$family, "gaussian")) {
                effect.threshold
```
- **Red → green:** main session, `89ad9c19` + new test file: 10 failures, 0 errors of 55 ("Expected `unique(fit$floors)` to be identical to 0.05." `actual`: 1.5), 10.8 s; at `685e0c5d`: 55 pass, 0 failures, 14.0 s.
- Implementer and verifier reproduced it: red 45 pass / 10 failures of 55 each (byte-identical outputs); green 55 / 0 failures each.
- **Expected values:**
  - RD 0.05 and IRD 0.01 by hand from `?forestsearch` "Resolved defaults, per estimand" (`man/forestsearch.Rd:1207-1208`) and the remap warning text; floor rule from `?forestsearch` `dina_args` (`man/forestsearch.Rd:233-234`). Fit floor == recorded and replay == fit are relational (the spec's own test).
  - Every value carries a comment; the verifier re-derived all, none copied from output. One minor flag: the IRD comment (`:159-161`) cites the no-threshold table, which does not literally cover a supplied out-of-range value.
- **Unaffected case:** spec: every configuration the existing test file covers, plus continuous outcome with a supplied threshold. `evidence/F14/unaffected.R`, 26 rows (fit + `args_call_all` replay) on start vs start + `F14.patch`: `identical()` TRUE for all 20 non-corner rows (floors and whole objects), `UNAFFECTED GATE: PASS`; walls 60.58 s / 64.04 s. Verifier re-run same (60.84 s / 64.29 s); its 18 extra rows change only at corner rows X8, X9.

| group (rows) | start floor | fixed floor | identical |
|---|---|---|---|
| continuous MD (R01-R06, R11-R12) | 0 / 0.4 / 1.5 | 0 / 0.4 / 1.5 | yes |
| corner RD/IRD, gaussian, 1.5 (R07-R08, R13-R16) | 1.5 | 0.05 (RD) / 0.01 (IRD) | no (intended) |
| other gaussian RD/IRD (R17-R22) | 0.07 / 0.05 / 0.02 | 0.07 / 0.05 / 0.02 | yes |
| binary OR (R23-R24) | 0.405465108108164 | 0.405465108108164 | yes |
| survival HR (R25-R26) | 0.22314355131421 | 0.22314355131421 | yes |
| helper `.forestsearch_dina_select()` (R09-R10) | 0.3 / 0 | 0.3 / 0 | yes |

- Cumulative run at head `61c15088`: every non-corner floor unchanged; R02/R05/R06 whole objects differ from F11's `tau_sign` (gate PASS with F11 removed).
- **Existing tests:** implementer 52 files pass (`test-mr-inference.R` 21/21 blocks, one per process); verifier 51 files, nb 5456, failed 0. Updated: `:144` `1.5` → `:150` `0.05` (title and comment reworded).
- **What a user sees** (RD, gaussian forced, `effect.threshold = 1.5`, `details = TRUE`; both trees print the same two remap warnings and record 0.05):
```
# start SHA
DINA selection: no subgroup met the criteria
  Harm floor:          m_diff = 1.5000
# fixed
DINA-selected subgroup: {x1 <= 2.75}
  Harm floor:          m_diff = 0.0500
```
- **Verification:** CONFIRMED; no repair round. Three minor problems: (1) IRD comment derivation — left as noted (committed test equals gated patch); (2a) `?forestsearch` sentence `man/forestsearch.Rd:1212-1214` inexact in the corner — open; (2b) A-001 NEWS bullet — fixed in NEWS commit `88302212` (`NEWS.md:390-396`); (3) "same warning" was two warnings — fixed in the section.
- **Departures:**
  - Test-only: `.gf_replay()` gained `df = .gf_df()` so the corner test can replay binary and count fits; bootstrap/CV replay checked via `do.call(forestsearch, args_call_all)`, `forestsearch_bootstrap_dofuture()` and `forestsearch_Kfold()` not run.
  - Spec scope read as R01-R10 plus MD with 1.5 (R11-R12); legacy `hr.threshold = 1.5` (R15-R16) and X8/X9 counted as corner, since the same remap fires. No roxygen change: `?forestsearch` `dina_args` already said the floor comes from the resolved `effect.threshold`.
  - `test-mr-inference.R` run block by block (120 s cap); skip-gated blocks run on a scratch copy holding three start-export files.
- Full evidence: detail file, section F14.


## 4. Existing expectations that were updated

| Fix | File:line (start → head) | Old | New | Why |
|---|---|---|---|---|
| F14 | `tests/testthat/test-dina-gaussian-floor.R:144` → `:150` | `expect_identical(unique(fit$floors), 1.5, info = route)` | `expect_identical(unique(fit$floors), 0.05, info = route)` | The pinning block named in the task: it encoded the old corner floor (the raw 1.5). The block's title (`:123` → `:124`, "… is still the floor as given" → "… is the floor as remapped") and its comment were reworded to match. |

- No other existing expectation changed. No existing test failed for any other reason; the existing test files that exercise each fix's functions were run on each fix's scratch copy (detail file, "Existing tests" of each fix), and the full suite at the head has 0 failures (§6).
- F3: no existing test sees one of the new warnings or the new stop message. The F3 verifier ran all 89 existing test files on the patched tree with a probe on the three new conditions; they fired only in the new test file. So no existing test needed `expect_warning()`.
- The full suite's warning count is unchanged at 38 (§6).

## 5. After the fixes

### 5.1 Default-path gate (step 1) — **PASS**

- The six equivalence cases of `REPORT_cran_prep_seed_legacy_2026-11-11.md` §3.3, run by `REPORT_cran_prep_code_fixes_2026-11-12_gate_run.R` (that report's `equiv_run.R` without its step-5 block): the package test helpers' synthetic data (survival N = 600, seed 42; `df.test` N = 300, seed 43; binary N = 600, seed 44), `set.seed(123)` before each call, sequential plan, `pkgload::load_all()` of each tree in its own `Rscript --vanilla` process.
- Trees: the start-SHA export `/tmp/cran_prep_code_fixes_2026-11-12/start` and a `git archive 61c15088` export of the package paths (`head_fixes`; `R/` is identical at `88302212`).
- Comparison: `REPORT_cran_prep_code_fixes_2026-11-12_gate_compare.R` (that report's `equiv_compare.R`, same timing-field rule), output `..._gate_compare.txt`.

| # | Call | `identical()` after dropping timing fields | Non-timing leaves that differ | Timing leaves dropped | Wall, start / head |
|---|---|---|---|---|---|
| 1 | `forestsearch()`, survival, defaults | TRUE | 0 | `$find.grps$time_search`, `$minutes_all` | 0.5 / 0.6 s |
| 2 | same, `mr_inference = TRUE` | TRUE | 0 | the two above + three `timing_seconds` | 7.7 / 7.6 s |
| 3 | same, with `df.test` | TRUE | 0 | 2 (as row 1) | 0.2 / 0.2 s |
| 4 | `forestsearch_bootstrap_dofuture()` on fit 1, `nb_boots = 10` | TRUE | 0 | `$results$tmins_search`, `$results$tmins_iteration` | 5.3 / 5.3 s |
| 5 | `forestsearch_Kfold()` on fit 1, `Kfolds = 5` | TRUE | 0 | `$timing_minutes` | 1.4 / 1.3 s |
| 6 | `forestsearch()`, binary outcome, defaults | TRUE | 0 | 2 (as row 1) | 0.6 / 0.6 s |

- A default survival `forestsearch()` run is identical. No difference anywhere in the six, so nothing needed tracing. `.Random.seed` after each call equals the state before it in both trees.
- Walls: 18.1 s (start) and 17.9 s (head), whole script.

### 5.2 Every fix's unaffected case at the head (all 14 fixes in)

Each fix was gated on a tree holding only its own patch. Read-only agents re-ran every fix's unaffected-case script with the fixed tree replaced by `head_fixes` (start tree unchanged; both sides re-run; scripts copied with only the tree and output paths changed, or byte-identical where the paths are arguments). Output: `/tmp/cran_prep_code_fixes_2026-11-12/cumulative/<F>/`.

| Fix | Cases | Identical to start at head | Differences, traced |
|---|---|---|---|
| F1 | 3 | 3 | — |
| F2 | 10 | 6 | 4 synthetic cases whose full-data subgroup has one factor. **F4**: fold 1 is `{age <= 50}` & `{sex >= 1}` against full-data `{age <= 50}`, so `exact_match` is 0 at head `R/forestsearch_cross_validation.R:1860` (1 at start `:1822`); Exact 0.25 → 0. Confirmed: on `head_fixes` with `F4.patch` reverse-applied, all 8 synthetic keys are identical to start. |
| F3 | 15 | 15 | — (the four stop cases stop on both trees, with the gated messages) |
| F4 | 14 | 14 | — |
| F5 | 3 | 3 | — |
| F6 | 4 | 4 | — (by F6's own spec gate: the only differing leaves are F6's intended `selected`/`lassokeep`, as at gating) |
| F7 | 5 | 5 | — (by F7's spec gate: the only differing leaves are F7's intended `P1`, `d_eff`, `L_eff_implied`, `otype`) |
| F8 | 9 | 9 | — |
| F9 | 8 | 8 | — |
| F10 | 12 | 12 | — |
| F11 | 16 | 16 | — |
| F12 | 6 | 6 | — (one ggplot `S7_class` closure attribute differs between any two runs, start vs start too; excluded as at gating) |
| F13 | 3 | 3 | — |
| F14 | 26 | 17 | 6 corner rows (R07, R08, R13–R16): F14 itself, intended. 3 non-corner rows R02, R05, R06 (continuous MD, `use_dina = TRUE`): whole objects differ through **F11**'s `tau_sign` on the screening route; their DINA floors are unchanged. Confirmed: identical to start on `head_fixes` minus `F11.patch`, and start + F11 alone reproduces head exactly. |

- 134 cases: 121 identical; 13 differ, each traced to a named fix (F4, F11, or F14's own corner). The F4 and F11 differences were confirmed by reverse-applying that fix. The 6 F14 corner rows were confirmed by matching the gated F14-only output (`evidence/F14/unaffected_fixed.rds`); F14 was not reverse-applied. No unexplained difference; nothing was repaired.

### 5.3 `paper_simulations` article (step 2) — not re-rendered

An independent agent listed every evaluated chunk of `vignettes/articles/paper_simulations.qmd` (35 chunks; no child documents or `source()`), followed the call graph from the evaluated calls to each fix's changed lines, and probed the article's own arguments on all 88 of its seeds (alt 1–44, null 1–44) with `get_FSdata` wrapped, sequential plan, nothing rendered. Its verdicts agree with each implementer's.

| Fix | Reaches changed lines | Calls and branch | Can the article's output change |
|---|---|---|---|
| F1 | no | — | no |
| F2 | no | no cross-validation call | no |
| **F3** | **yes** | `run_simulation_analysis()` (`run-alt-sims`, `run-null-sims`): the F3a checks `R/run_simulation_analysis.R:1093`/`:1101`, `:1175` (FS arm, when a subgroup is found) and `:1491-1492` (GRF arm) | No: the stop fires only with NA membership; the simulated GBSG covariates have no NA (`anyNA` FALSE on all 88 seeds); and it would stop in exactly the cases that stopped before. |
| F4 | no | — | no |
| F5 | no | — | no |
| **F6** | **yes** | both FS arms: `forestsearch(use_lasso = TRUE, outcome_type = "survival")` → `lasso_selection()` (`R/get_FSdata_helpers.R:346-348`); the label (`R/forestsearch_main.R:3393-3403`) under `details` | No: the Cox LASSO has no `(Intercept)` row and its row names equal the confounders, so `selected`/`omitted` are the same sets (all 88 seeds); the survival label stays "Cox-LASSO". |
| F7 | no | — | no |
| **F8** | **yes** | FS arm: `fs_params$use_grf = TRUE` overrides the arm's `use_grf = FALSE` (`modifyList`, `R/run_simulation_analysis.R:587-591`), so `get_FSdata()`'s GRF merge (`R/get_fsdata.R:501-502`) runs (77 of 88 seeds; 11 null seeds had no GRF cuts). Median-removal and DINA hunks not reached. | No: at the merge the pool is the LASSO-filtered categorical names, and no categorical name is a substring of another confounder's cut; the kept cuts are the same under both rules. |
| F9 | no | no MRCT call | no |
| F10 | no | survival only (`grf.subg.harm.survival()`) | no |
| F11 | no | `use_dina` is FALSE | no |
| F12 | no | no display of a fit | no |
| F13 | no | survival only | no |
| F14 | no | `subgroup_method = "consistency"`, `use_dina = FALSE` | no |

F3, F6 and F8 reach changed lines, so they are listed under OPEN ITEMS (1).

### 5.4 Other articles (read-only sweep, not re-rendered)

The other 12 documents under `vignettes/` were swept the same way. Where a fix's lines are reached, the evaluated chunks were run in stages on both trees (sequential plan, counters on the changed statements) and the outputs compared. **In every case the rendered output is unchanged.**

| Document | Fixes reached | Why the output does not change |
|---|---|---|
| `articles/forestsearch-2024.Rmd` | F1, F2, F3, F4, F6, F8, F12 | F1: the bin tables are not printed and the values are off the edges. F2/F3: GBSG has no NA. F4: its full-data subgroup `{er <= 0}` is one factor, so the new Exact rule runs; with seed 8316951 neither fold contains `{er <= 0}`, so Exact is 0 in both trees (a different fold result could lower it). F6: survival LASSO, same sets and label. F8: the GRF cut `er <= 0` is kept under both rules. F12: survival rows do not use the GLM estimator. |
| `articles/survival-post-selection.qmd` | F1, F3, F14 | F1: one replicate has `N_sg = 100`, which moves from "50-99" to "100-149" in the hidden `size_dist` (21/2 instead of 22/1); the article prints only `summ_fs$table`, which is identical. F3: no NA. F14: survival, the old and new floor arguments are the same value, and the cox family does not use it. |
| `articles/actg175-glm-outcomes.qmd` | F3, F14 | No NA; OR is a ratio measure, so the same value passes and the binomial family does not use the gaussian floor. |
| `forestsearch.qmd` | F3, F14 | No NA; survival, same value. |
| `articles/mrct-consistency-simulations.qmd` | F9 | `mrct_case_study()` builds the new default (it drops `z_cd40_treat`, `z_cd40_k`, `z_cd40_k_treat`), but `fs_args` supplies `confounders.name`, which replaces it: 0 differing leaves. |
| The other 7 documents | none | No evaluated call reaches a changed line. |

## 6. Full suite and `check()` (steps 4 and 5)

| Run | Part 3 (`REPORT_cran_prep_docs_legacy_2026-11-12.md` §8.4–8.5) | This branch | Result |
|---|---|---|---|
| Full suite (`devtools::test()`, `NOT_CRAN=true`, in the repository at `88302212`) | 652 blocks, 89 files, 8,478 expectations, 0 failures, 3 skips, 38 warnings; 694 s | **703 blocks, 102 files, 8,853 expectations, 0 failures, 0 errors, 3 skips, 38 warnings**; wall **730 s** (12:34:28–12:46:39; cap 30 min) | **PASS** |
| `devtools::check(document = FALSE, error_on = "never")` | 0 errors / 0 warnings / 0 notes; 1,035 s | **0 errors / 0 warnings / 0 notes**, Status OK; wall **1,092 s** (12:34:28–12:52:41; tests 11 m; examples 36 s; `--run-donttest` 64 s; vignettes 79 s; cap 60 min) | **PASS** |

- **Difference, accounted for** (`..._suite_diff.txt`, per file against part 3's saved results): +51 blocks, +13 files, +375 expectations. These are exactly the 13 new test files (blocks/expectations: bootstrap-summary-bins 3/8, cv-exact-match 2/6, cv-fold-na-assignments 2/6, dina-screening-tau-sign 2/16, display-orientation 6/93, duplicate-cut-names 7/13, fpr-calibration-fallback 4/16, grf-count-orientation 2/14, hte-crump-confounders 5/35, lasso-intercept 4/62, missing-membership-messages 10/60, mrct-default-confounders 1/8, rd-standard-error 2/14) plus F14's corner block in `test-dina-gaussian-floor.R` (5/31 → 6/55). No other file changed counts. The 3 skips are the same (`test-cross-outcome-parity.R` ×2, `test-cv-no-subgroup-edges.R`); warnings are unchanged at 38.
- **Check route** (`REPORT_cran_prep_2026-11-11.md` §6 with the `-z` correction): an export of the working tree at `88302212`, `git ls-files -z -co --exclude-standard`, 9,697 files listed and 9,697 copied, into `/tmp/cran_prep_code_fixes_2026-11-12/check/forestsearch`. All 36 files changed since the start SHA were `cmp`-identical to the repository's. `forestsearch_0.4.1.9010.tar.gz` was built there. The local `legacy_vignette/` (git-excluded) was not exported, so its known NOTE did not arise.
- **Contention.** The suite and `check()` started together, and six read-only agents and fourteen drafting agents were running, so both wall clocks include contention.
- The suite left the working tree clean (`git status --porcelain` empty after it).

## 7. Departures and decisions

- **F2, the branch condition.** The task says to take the two-level branch "only when the table has two rows and two columns". The commit keeps `n_sgfound == 2` and adds the 2×2 condition. `n_sgfound` counts NA as a value, so {0, 1, NA} gives 3 and took the single-level branch at the start SHA; a condition on the table shape alone would send it to the two-level branch and break the unaffected case the task names. The task also forbids changing how `n_sgfound` is counted. Both earlier reports recommend this variant.
- **F3 and F11, comments added after verification.** The main session added the three comments the verifiers asked for (no code, no value) before committing: one line in F3, and in F11 a reworded citation (2 lines become 3) and a new two-line comment. The F11 test grew from 129 to 132 lines and the F3 test from 362 to 363. The verified patches are kept in scratch (§1).
- **F3c gives the count only.** The task's F3c bullet asks for "one warning per exported call with the count". Each fold has its own subgroup definition, so the covariates are not named there; F3a and F3b name them.
- **F5, repaired once.** The first version parsed every label, which dropped or mis-mapped categorical confounders with non-syntactic names (`"bm status"`, `"ER-pos"`). The repair maps a label that is exactly a column name to that column and parses the rest; a second verifier confirmed it, and the unaffected case covers both kinds of name.
- **F6, the mocked-LASSO test.** The task's second test ("with `lasso_selection()` mocked to select nothing, `get_FSdata()` keeps the full pool") cannot go red: it replaces the function being fixed, and `get_FSdata()`'s rule was already right. It is kept as a lock on that rule. The red for that bullet comes from test 3, which runs the real `lasso_selection()` with `coef()` mocked to keep only the intercept; on the start tree it stops with "no candidate cut expressions could be constructed … LASSO selected (1)". F6 also corrected the roxygen description and added three `@return` lines describing the changed contract, beyond the title the task names.
- **F14, the corner.** The legacy spelling `hr.threshold = 1.5` (rows R15, R16) is counted as part of the corner, because the same remap fires.
- **A second NEWS commit.** The task puts version and NEWS in one commit (`88302212`). The report review then found three NEWS bullets that said something the code does not do (F2, F10, F13; OPEN ITEMS 16), and the task asks each entry to say what a user sees change. They were corrected in `a577b59d`, `NEWS.md` only, wording only, rather than left wrong; history was not rewritten. The suite and `check()` ran at `88302212`; `a577b59d` changes no code, test or Rd file.
- **Commit order.** F3 and F11 were committed after F14 (their verifiers finished last). The fixes are independent and every patch applied cleanly.
- **Report split.** The per-fix evidence is long, so it is in `REPORT_cran_prep_code_fixes_2026-11-12_detail.md`, as part 3 of the 2026-11-11 work did with its detail file. Each entry in §3 above has every item the task asks for; the detail file has every table.
- **Verifier comment requests on committed tests.** For F1, F5, F6, F8, F9, F10, F12, F13 and F14 the verifiers asked for extra comments on a few expected values or assertions. The values are correct (each verifier re-derived them), and the commits were kept equal to the gated patches; the requests are listed in OPEN ITEMS (17).

## 8. OPEN ITEMS

1. **`paper_simulations` reaches three fixes' changed lines (step 2).** None can change the article's output (§5.3); it was not re-rendered.
   - F3: the F3a checks in `run_simulation_analysis()` (head `R/run_simulation_analysis.R:1093`/`:1101`, `:1175`, `:1491-1492`) run in every replicate that finds a subgroup. They stop only with NA membership; the simulated covariates have none.
   - F6: `lasso_selection()` runs in every FS fit (`use_lasso = TRUE`, survival). Cox LASSO has no intercept, so `selected`/`omitted` are unchanged.
   - F8: `fs_params$use_grf = TRUE` overrides the FS arm's `use_grf = FALSE`, so the GRF merge (head `R/get_fsdata.R:501-502`) runs in 77 of 88 seeds. The kept cuts are the same under both rules.
2. **F13 (recorded as OPEN by the task).** The covariate-adjusted logistic tier still uses `se_rd <- abs(avg_deriv) * se_treat` (head `R/glm_effect_estimators.R:546`). It leaves out the control arm's variance and degenerates when an arm is at 0% or 100%. With both arms degenerate (control 0%, treated 100%) the new unadjusted SE is also exactly 0. Related:
   - the comment "conservative for subgroup estimation" (head `:531-532`) no longer holds for the form that remains and contradicts the new comment;
   - `.make_poisson_rate_estimator()`'s `dr_gcomp` IRD branch has the same shortcut (head `:856`);
   - with `adjust_covariates` and `ps_adjust_method = "dr_gcomp"`, or `"iptw"` without an `sw` column, the closure drops the covariates before `.estimate_rd` (`:415-419`), so those direct calls now get the new SE. `forestsearch()` never builds this. Weights never reach `.estimate_rd`'s logistic tier: the IPTW path returns earlier.
3. **F12: other displays still use the estimator's default orientation.** `sg_tables()` (`R/summary_utility_functions.R:829-839` and `:982-992`) and the CV summary tables (head `R/forestsearch_cross_validation.R:1509-1517`, `:1565-1575`) pass no `adverse_outcome`. For a fit whose resolved value is FALSE (a continuous fit at the default), their effect column now has the opposite sign to the fit and to the three fixed displays (example: `sg_tables()` Questionable MD −1.094 against `max_sg_est` 1.094). `compute_frontier_cis()` falls back to `%||% TRUE` (`R/frontier_cis.R:150`). The same one-line change would apply.
4. **F14, documentation and messages.**
   - `?forestsearch` (`man/forestsearch.Rd:1212-1214`, roxygen `R/forestsearch_main.R:1313-1315`) says a supplied threshold "suppresses the corresponding remap and the supplied value is used"; an out-of-range RD/IRD value (> 1) is remapped with a warning, and since F14 the remapped value is also the DINA floor.
   - The remap warning hard-codes "Remapping to default: 0.05" and "OR scale" for IRD too, whose default is 0.01 (head `R/forestsearch_main.R:2334-2350`); the earlier validator warning (`:301-309`) says "A risk difference cannot exceed 1.0" for IRD. Message text only.
   - With `dina_args = list(family = "gaussian")` forced on a ratio estimand (survival, OR, IRR), the raw ratio threshold is used as an identity-scale floor. Same on both trees; outside F14.
5. **F2: shapes that now reach the single-level fallback.** Besides {1, NA} and {0, NA} against a two-level original, these shapes stopped before and now return values the fallback labels by position (head `R/forestsearch_cross_validation.R:1473-1477`): 2×1 (full data all Hc) `0.4 NA 1.0 NA`; 1×1 `1 NA 1 NA`, labelled H; 1×2 with a zero cell `NaN NA 0 NA`. The `outall = TRUE` assembly still branches on `n_sgfound == 2` (head `:1622`, `:1669`, `:1729`, `:1762`); with NA in `treat.recommend.original` it stops with "incorrect number of dimensions" for every fold shape on both trees (pre-existing).
6. **F7.**
   - Exact matching rejects abbreviations that `forestsearch()` accepts through `match.arg()` (for example `fs_params$outcome_type = "surv"`); at the start SHA such a call ran with the binary approximation. A later fix could use `match.arg(otype, otype_ok)`.
   - `fpr_calibration_transfer()` discards the nested call's error (head `R/fpr_calibration.R:569`), so a misspelt type passed through `...` now gives NA results and `n_sim` 0 with no message (it gave `fpr_sim` 0, `n_sim` 2).
   - `"count"` uses `d_eff_binary()` although `d_eff_count()` exists (pre-existing). `@param prop_cens` and `@param p_event` describe the explicit type only.
7. **F1: the top consistency bin is never highlighted.** `R/format_subgroup_summary_tables.R:222` highlights `c("0.85-0.95", ">= 0.95")`, with a space; the label is `">=0.95"`. Pre-existing; after F1 more iterations land in that row. The roxygen does not say the bins are left-closed.
8. **F11.**
   - With `dina_args = list(selected_only = FALSE)`, the screening pool is `dina_frontier()`'s `cut_expr`, which reads the raw tau-hat (head `R/forestsearch_main.R:3314`; A-020 of the 2026-11-11 report).
   - `dina_select_statistic = "effect"` re-selects only on the `dina` route (`R/forestsearch_helpers.R:1611-1624`), so "the same cut `subgroup_method = "dina"` would select" (`R/forestsearch_main.R:440-445`) holds for the default statistic only.
   - `quarto/dina/smoke-tests/dina_screening_smoketest.qmd` runs a continuous battery with `use_dina = TRUE` (:112, :126, :134, :169); its cuts and pass/fail checks may change. `quarto/` is read-only here.
9. **F10.** The GRF CATE footnote is labelled "(risk-difference scale)" (`R/summary_utility_functions.R:917`); for counts the label is wrong, and with `adverse_outcome = TRUE` the printed sign is now reversed relative to before. The count flip prints no verbose message, and `grf_y_transform` does not record the flip (as for binary and continuous).
10. **F9.** `mrct_region_sims()` still builds `fs_call_args <- modifyList(fs_args, fs_base_args)` (head `R/mrct_simulation.R:561`), which replaces a `confounders.name` passed in `fs_args` with the default. `analysis_covariates()` still drops columns by the suffix pattern `_treat$|_k$` (`R/mrct_inject_dgm.R:635`, A-124), which would also drop `z_prior_treat`.
11. **F3.** A loop that calls `forestsearch(subgroup_method = "grf")` on data with missing covariates repeats the F3b warning once per replicate (`run_simulation_analysis()` with a GLM DGM, `fpr_calibration()`, `mrct_simulation()`; `mrct_trial_analysis()`/`mrct_case_study()` now carry the text in `$warnings` when a `grf_res` is supplied). No message for NA in `df.predict`/`df.test` on the GRF route, or from a direct `forestsearch_KfoldOut()` call. With a count of 1 the messages read "1 of N subjects have".
12. **F5.** The roxygen `@return` (`R/test_hte_crump.R:386`) does not describe the mapping or the warnings; the kept `source` text "confounders.evaluated (GRF/LASSO selected)" is wrong on fits that use neither. `quarto/applications/validation_hte_tests_crump.qmd:1366` calls the function; its printed covariates would now include continuous ones (not re-rendered).
13. **F4.** The wording "Exact match on all factors" (`R/cv_summary_tables.R:413`, `:418`) and `CV_sgs()`'s `@return` are unchanged. With `maxk >= 3` the fold refits keep only two factors (`R/forestsearch_cross_validation.R:505-506` and `:1034-1035` at start; `:505-506` and `:1046-1047` at head), so a depth-2 "Exact" can over-count. `vignettes/articles/forestsearch-2024.Rmd` (chunk `kfold-cv`) runs the new rule; its displayed Exact is unchanged with the current folds and could be lower at a future render.
14. **F6.** A binary outcome with `offset.name` fits a Poisson LASSO, but both printed labels say "Logistic-LASSO" (`R/get_fsdata.R:343-348`, head `R/forestsearch_main.R:3396-3401`); the task said to mirror the existing switch. The comment at `R/get_fsdata.R:327` omits count.
15. **F8.** `conf.cont_medians_force` is appended to `conf.cont_medians` twice (`R/get_fsdata.R:434-436` and head `:469`), so with GRF/DINA off and no `conf_force` a forced median cut appears twice in the pool (same on both trees). `get_cut_name()` still uses an unanchored `grepl()` (head `R/get_FSdata_helpers.R:463-466`); its only caller `is_flag_drop()` is never called.
16. **NEWS wording — corrected in `a577b59d`.** The report review found three part-4 NEWS bullets inexact at `88302212`: F2 said `forestsearch_Kfold()` used to stop (it caught the error and returned all-NA summaries; only `forestsearch_tenfold()` stopped); F10 called `adverse_outcome = TRUE` "the default for counts" (it is the value `forestsearch()` and `run_simulation_analysis()` resolve; a direct `grf.subg.harm.glm()` call defaults to FALSE, `R/grf_subg_harm_glm.R:304`); F13 said the covariate-adjusted case is unchanged (a direct call with `"dr_gcomp"`, or `"iptw"` without `sw`, drops the covariates and now gets the new form; item 2). F9's "`confounders.name = NULL`" was also loose. All four were reworded in `a577b59d` (`NEWS.md` only). Nothing left open.
17. **Verifier comment requests on committed tests (no value changes).** F1 `test-bootstrap-summary-bins.R:76-83` (label-order checks pass on both trees; anchor row order); F5 `test-hte-crump-confounders.R:97-99`; F6 a note that test 2 passes before and after; F8 `test-duplicate-cut-names.R:162-178` (the exact-duplicate guard also holds through `unique()` and does not by itself pin the `%in% confs` clause); F9 `test-mrct-default-confounders.R:40-49`, `:63`, `:82`; F10 `test-grf-count-orientation.R:103-106`, `:129-132`; F12 `test-display-orientation.R:27-30`, `:127-130` (tolerance comments); F13 `test-rd-standard-error.R:78-81`; F14 `test-dina-gaussian-floor.R:159-161` (cite the remap "to default" step). Line numbers at each fix's commit.
18. **`tests/testthat/test-mr-inference.R`** takes more than 120 s as one process (scratch runs under contention hit the reproducer cap and were run block by block). It passes in the full suite (182.4 s).
19. **Other articles, no current output change (§5.4).** `survival-post-selection.qmd`'s hidden `size_dist` moves one replicate (`N_sg = 100`) from "50-99" to "100-149" (F1); the article prints neither distribution.

## 9. Final gates (step 8)

`REPORT_cran_prep_code_fixes_2026-11-12_final_gates.sh` recomputes the gates from git alone: G1 working tree, G2 both sibling repositories (branch, SHA, empty status, no file newer than the start commit), G3 `git diff 76a62a71 --stat` lists only files this report names and nothing under `vignettes/` or `quarto/`, G4 each fix commit touches only its own files, G5 nothing pushed or merged. It was run on the staged tree just before the report commit (output `..._final_gates.txt`, committed) and again after the commit.

| Gate | Result on the staged tree |
|---|---|
| G1 | PASS: only the 10 report files are staged; nothing else modified or untracked |
| G2 | PASS: fs-post-selection `main` `5d165acc…`, fs-glms-interpretable `main` `dc94a17e…`, both `git status --porcelain` empty, 0 files newer than the start commit (2026-11-12 11:37:04) |
| G3 | PASS: 46 paths in `git diff 76a62a71 --stat` including the staged report files; every one is named in this report (§10); nothing under `vignettes/` or `quarto/` |
| G4 | PASS: each of the 17 commits since the start touches exactly its own files (the 14 fixes as in §10, `88302212` `DESCRIPTION` + `NEWS.md`, `a577b59d` `NEWS.md`, `76a62a71` the task file) |
| G5 | PASS: HEAD is on no remote-tracking ref; 0 merge commits since `76a62a71` |

The first run on the staged tree failed G3 only: two file names in §10 were abbreviated (`..._gate_compare.txt`, `..._final_gates.txt`). They were written out in full. The next run failed G1 only, because the edited report had not been staged again. With every report file staged, two consecutive runs gave byte-identical output, all PASS; that output is the committed `REPORT_cran_prep_code_fixes_2026-11-12_final_gates.txt`. The gates were run again after the commit (G1 then requires an empty `git status --porcelain`).

## 10. Files

**Changed since the start SHA by the fixes and the two NEWS commits** (36 files; `git diff 76a62a71 a577b59d --stat`, the same list as at `88302212`):

| Kind | Files |
|---|---|
| Code (16) | `R/forestsearch_cross_validation.R` (F2, F3, F4), `R/forestsearch_helpers.R` (F3), `R/forestsearch_main.R` (F6, F11, F14), `R/fpr_calibration.R` (F7), `R/get_FSdata_helpers.R` (F6), `R/get_fsdata.R` (F8), `R/glm_effect_estimators.R` (F13), `R/grf_subg_harm_glm.R` (F10), `R/mr_estimates_table.R` (F12), `R/mrct_case_study.R` (F9), `R/mrct_simulation.R` (F9), `R/plot_sg_glm_outcomes.R` (F12), `R/plot_subgroup_results_forestplot.R` (F12), `R/run_simulation_analysis.R` (F3), `R/summarize_bootstrap_subgroups.R` (F1), `R/test_hte_crump.R` (F5) |
| Help pages (4) | `man/fpr_calibration.Rd` (F7), `man/lasso_selection.Rd` (F6), `man/mrct_case_study.Rd` (F9), `man/mrct_region_sims.Rd` (F9) |
| Tests (14) | new: `tests/testthat/test-bootstrap-summary-bins.R` (F1), `tests/testthat/test-cv-fold-na-assignments.R` (F2), `tests/testthat/test-missing-membership-messages.R` (F3), `tests/testthat/test-cv-exact-match.R` (F4), `tests/testthat/test-hte-crump-confounders.R` (F5), `tests/testthat/test-lasso-intercept.R` (F6), `tests/testthat/test-fpr-calibration-fallback.R` (F7), `tests/testthat/test-duplicate-cut-names.R` (F8), `tests/testthat/test-mrct-default-confounders.R` (F9), `tests/testthat/test-grf-count-orientation.R` (F10), `tests/testthat/test-dina-screening-tau-sign.R` (F11), `tests/testthat/test-display-orientation.R` (F12), `tests/testthat/test-rd-standard-error.R` (F13); changed: `tests/testthat/test-dina-gaussian-floor.R` (F14) |
| Version and NEWS (2) | `DESCRIPTION`, `NEWS.md` |

The task file `dev/tasks/TASK_cran_prep_code_fixes_2026-11-12.md` is in the start commit itself.

**Beside this report** (committed with it):

| File | What it is |
|---|---|
| `REPORT_cran_prep_code_fixes_2026-11-12.md` | this report (also copied to `~/Downloads`) |
| `REPORT_cran_prep_code_fixes_2026-11-12_detail.md` | the full per-fix evidence: sites, every red/green run, every expected value, unaffected cases, existing tests, verifier problems, departures, open items (also copied to `~/Downloads`) |
| `REPORT_cran_prep_code_fixes_2026-11-12_apply_fix.sh` | main-session per-fix gate: test file only → red; code → green |
| `REPORT_cran_prep_code_fixes_2026-11-12_commit_fix.sh` | main-session per-fix commit: staged diff equals the gated patch; `--stat`; commit |
| `REPORT_cran_prep_code_fixes_2026-11-12_gate_run.R`, `REPORT_cran_prep_code_fixes_2026-11-12_gate_compare.R`, `REPORT_cran_prep_code_fixes_2026-11-12_gate_compare.txt` | default-path gate (§5.1) and its output |
| `REPORT_cran_prep_code_fixes_2026-11-12_suite_diff.txt` | per-file suite counts against part 3 (§6) |
| `REPORT_cran_prep_code_fixes_2026-11-12_final_gates.sh`, `REPORT_cran_prep_code_fixes_2026-11-12_final_gates.txt` | final gates (§9) and their output on the staged tree |

**Scratch, not in the repository** (`/tmp/cran_prep_code_fixes_2026-11-12/`): `start/` (start export), `wt/F*/` (per-fix scratch copies), `evidence/F*/` (implementer runs and notes), `verify/F*/` (verifier copies and runs), `patches/` (gated patches; `F3.verified.patch`, `F11.verified.patch`), `main/F*/` (main-session red/green), `cumulative/` (§5.2), `gate/`, `suite/`, `check/`, `report_sections/` (drafts of the detail file).
