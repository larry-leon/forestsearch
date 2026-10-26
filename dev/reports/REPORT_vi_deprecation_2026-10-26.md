# REPORT — deprecate `vi.grf.min` and `max_n_confounders`; remove GRF variable-importance screening

**Task:** `dev/tasks/cc_task_vi_deprecation_2026-10-26.md` · **Executed:** 2026-10-26 by CC, unattended
**Category:** `R/` change, CHANGES BEHAVIOUR (documented API deprecation), decided by Larry
**Record directory:** `dev/reports/` exists; battery under `dev/reports/vi_deprecation_2026-10-26/`

## 1. Provenance

| item | value |
|---|---|
| start branch | `master` |
| base SHA (master HEAD) | `a0025013f103d0921455102129caeb33127bb139` |
| `git status --porcelain` at start | `?? quarto/extreme_subgroups/fixed_random/` (pre-existing, never staged) and `?? dev/tasks/cc_task_vi_deprecation_2026-10-26.md` (this document, copied in just before the status was taken) |
| branch created | `feature/vi-deprecation-lean-fits` (from base SHA) |
| task-document filenames | `~/Downloads/cc_task_vi_deprecation_2026-10-26.md` (hyphens intact) → `dev/tasks/cc_task_vi_deprecation_2026-10-26.md` |
| MD5 at base | `R/forestsearch_main.R` 11d80e74c4fcdba43fe8bbf9a1287d3f, `R/get_fsdata.R` a7f71d659558706ae121d72e82099181 (match Stage-0b record) |
| installed build before (Stage A) | Version 0.4.0, Packaged 2026-10-10 04:12:05 UTC, Repository CRAN, Built R 4.6.1; 2026-10-26 16:12:28 UTC — the CRAN build |
| installed build after (§3, Stage C) | Version 0.4.0.9000, Built R 4.6.1; 2026-10-26 17:05:25 UTC — **the installed package is now the branch build** |

Commits, in order:

1. `4654b4cc` tasks: vi.grf.min / max_n_confounders deprecation (remove GRF VI screening) — task document alone
2. `94738533` dev: Stage A baselines, vi deprecation
3. `41347eb4` Deprecate vi.grf.min and max_n_confounders; remove VI screening — all §3 edits, one commit
4. Stage C battery + this report (see §7 closing block for the SHA)
5. Report closing block (post-conditions)

## 2. The edits

Every named 0.4.0 line was verified by content before it was edited. Line drift against the document:

- `R/fs_family_report.R`: the block "L268–L283 (`vi_on <- …` through the `add("confounder cap", …)` call)" actually ends at **L282**. The `add("confounder cap", …)` call closes on L282, and L283 is `add("dummy expansion", …)`. The replacement covers L268–L282 by content, and the `dummy expansion` row is untouched.
- Every other named line was at its stated 0.4.0 number. After the edits, the formals sit at L1519 / L1569 (the roxygen above them is shorter). Their text and position in the formals list are unchanged.

**`R/forestsearch_main.R`**
- Formals `max_n_confounders = 1000` and `vi.grf.min = NULL` are unchanged.
- The deprecation guard is inserted verbatim (document §3 item 2) right after `consistency_method <- match.arg(consistency_method)`.
- Deleted: the inert-cap warning block (0.4.0 L3213–L3225) and the whole SECTION 5 VI block (L3227–L3313), 101 lines in total. They are replaced by the 5-line `conf.screen <- FSconfounders.name` block (verbatim). Post-edit grep: none of `X`, `Y_vi`, `cs.forest`, `vi.cs`, `vi.cs2`, `vi.order`, `vi_max`, `vi_ratio`, `selected.vars`, `tau.rmst`, `tune_arg`, `y_evt_c`, `y_evt_t` appears on any code line of the file. `conf.screen` is assigned once (L3223) and used at L3235 (`df.confounders <- df[, conf.screen, drop = FALSE]`).
- Roxygen (Rd markup): `@param vi.grf.min` and `@param max_n_confounders` are replaced with the deprecation text; `max_n_confounders` also gets "It capped the importance-ordered column list, which no longer exists." The `@param use_grf` sentence and the `@param frac.tau` sentence are replaced as specified.
- Sweep `grep -n -i "variable.importance\|section 5\|vi_ratio\|vi_max"`: what remains is the new deprecation text (L388, L507, L917, the guard L1608/L1616/L1625, the replacement comment L3220) plus the `@importFrom` line (L1484, kept). The L2051 comment (`# Survival: not used (causal_survival_forest handles sign correctly)`) is in the `adverse_outcome` default block and does not concern the VI step, so it was left as is.
- `#' @importFrom grf causal_survival_forest variable_importance` is kept; `NAMESPACE` is byte-identical.
- `fit_causal_forest_glm()` (`R/grf_helpers.R`) is left in place; see OPEN ITEMS.

**`R/fs_family_report.R`**
- The two rows ("GRF variable-importance ordering", "confounder cap") are replaced with the single `"GRF variable-importance screening (deprecated)"` row, status `INE`, verbatim. `vi_on` and `cap_status` are referenced nowhere now (grep empty). The `.arg()` reads (L142–L143) and `n_cut_columns` are kept.
- Roxygen L14, L38–L40 and both examples are edited as specified.

**Other files**
- `R/run_simulation_analysis.R`: deleted `vi.grf.min = -0.2,` from the default list.
- `R/get_fsdata.R`: the diagnostic string now reads `"(4) GRF-related thresholds (dmin.grf, frac.tau) "`.
- `tests/testthat/helper-synthetic-dgm.R`: deleted `vi.grf.min = -0.2,` from `.fs_args_for()`.
- `tests/testthat/test-vi-deprecation.R`: new. It holds the four test blocks of §3 item 14: deprecation warnings; inert (`sg.harm`, sorted membership, `out_sg$result`, `hr.subgroups` all `identical()`); no forest (0 fits through `count_grf_fits()`); and the echo in `fs_family_report()`.
- `vignettes/forestsearch.Rmd`: deleted the `max_n_confounders` and `vi.grf.min` lines. Not built.
- `inst/doc/`: **does not exist** in the working tree, and nothing under it is tracked. Nothing to leave alone.
- `NEWS.md`: master already had an empty `# forestsearch (development version)` heading above `# forestsearch 0.4.0`, so no `0.4.0.9000` heading was added. The `### Deprecated` entry sits under the existing heading and covers every point in §3 item 16.

**Regenerated by `devtools::document()`:** `man/forestsearch.Rd`, `man/fs_family_report.Rd` (only these two). `NAMESPACE` is unchanged.

**Verbatim test-file edits (§3 item 13),** from `git diff a0025013 HEAD -- tests/testthat/test-fs-family-report.R`:

```diff
-test_that("resolution agreement (drift guard): stop_threshold and the max_n_confounders cap match the engine", {
+test_that("resolution agreement (drift guard): stop_threshold matches the engine; the deprecated VI row is inert", {
@@
-  # max_n_confounders inert under vi.grf.min = NULL, on the fitted object and on its args
   expect_null(fit$args_call_all$vi.grf.min)
-  expect_identical(r[r$stage == "confounder cap", "status"], "inert")
-  expect_identical(r[r$stage == "GRF variable-importance ordering", "status"], "disabled")
+  expect_identical(r[r$stage == "GRF variable-importance screening (deprecated)", "status"], "inert")
@@
-  r_null <- fs_family_report(.fr_driver_args(vi.grf.min = NULL))
-  expect_identical(row_status(r_null, "GRF variable-importance ordering"), "disabled")
-  expect_identical(row_status(r_null, "confounder cap"), "inert")
-  r_num <- fs_family_report(.fr_driver_args(vi.grf.min = -0.2))
-  expect_identical(row_status(r_num, "GRF variable-importance ordering"), "data-dependent")
-  expect_identical(row_status(r_num, "confounder cap"), "data-dependent")
+  # deprecated VI row: inert whatever is passed
+  expect_identical(row_status(fs_family_report(.fr_driver_args(vi.grf.min = NULL)), "GRF variable-importance screening (deprecated)"), "inert")
+  expect_identical(row_status(fs_family_report(.fr_driver_args(vi.grf.min = -0.2)), "GRF variable-importance screening (deprecated)"), "inert")
+  expect_identical(row_status(fs_family_report(.fr_driver_args(max_n_confounders = 5)), "GRF variable-importance screening (deprecated)"), "inert")
```

L17 (`.fr_driver_args()` keeps `vi.grf.min = -0.2`), L147–L148, L159 and L177 are untouched.

**`git diff --stat a0025013 41347eb4 -- . ':!dev'`** (package files only; the dev/ battery is listed in §7):

```
 NEWS.md                                |  21 ++
 R/forestsearch_main.R                  | 178 +++--------
 R/fs_family_report.R                   |  31 +-
 R/get_fsdata.R                         |   2 +-
 R/run_simulation_analysis.R            |   1 -
 man/forestsearch.Rd                    |  49 +--
 man/fs_family_report.Rd                |  12 +-
 tests/testthat/helper-synthetic-dgm.R  |   1 -
 tests/testthat/test-fs-family-report.R |  16 +-
 tests/testthat/test-vi-deprecation.R   |  76 +++++
 vignettes/forestsearch.Rmd             |   2 -
```

## 3. Stage A vs Stage C

Fixtures are transplanted verbatim from `dev/glm-continuous-sims/stage0b_glm_2026-10-26/` into `vi_common.R`. Run (i) passes the args as committed. Run (ii) passes `vi.grf.min = NULL` explicitly and omits `max_n_confounders`. Every comparison is made on saved `.rds` against saved `.rds` (see §6, deviation 2).

**Stage A (CRAN 0.4.0)** — GATE PASS (binary (i) = 1 fit, all others 0):

| fixture | run | n_fits | warnings | sg.harm |
|---|---|---|---|---|
| survival | i | 0 | 0 | {er <= 0} & {pgr <= 26} |
| survival | ii | 0 | 0 | {er <= 0} & {pgr <= 26} |
| continuous | i | 0 | 0 | !{cd40 <= 507} & {gender} |
| continuous | ii | 0 | 0 | !{cd40 <= 507} & {gender} |
| binary | i | **1** | 0 | !{wtkg <= 86} & !{cd40 <= 380} |
| binary | ii | 0 | 0 | !{cd40 <= 380} & !{wtkg <= 86} |

**Stage C (branch 0.4.0.9000):** every run has n_fits = 0. Warnings are 0 everywhere except binary (i), which has exactly 1 (the `vi.grf.min` deprecation warning, verbatim in `stageC_summary.txt`).

Each component compared against its Stage-A counterpart (`identical()`). The first five rows compare to the matching Stage-A run; binary (i) is shown against both Stage-C (ii) and Stage-A (ii):

| component | surv i | surv ii | cont i | cont ii | bin ii | bin i vs C(ii) | bin i vs A(ii) |
|---|---|---|---|---|---|---|---|
| sg.harm | identical | identical | identical | identical | identical | identical | identical |
| membership (sorted) | identical | identical | identical | identical | identical | identical | identical |
| out_sg$result | identical | identical | identical | identical | identical | identical | identical |
| hr.subgroups | identical | identical | identical | identical | identical | identical | identical |
| filter_counts | identical | identical | identical | identical | identical | identical | identical |
| n_candidates_* | identical | identical | identical | identical | identical | identical | identical |
| deprecation warnings | 0 | 0 | 0 | 0 | 0 | 1 | — |
| forest fits (A → C) | 0 → 0 | 0 → 0 | 0 → 0 | 0 → 0 | 0 → 0 | 1 → 0 | — |

**Survival `nb = 20` bootstrap payload:** found `results`, `SG_CIs`, `H_estimates`, `Hc_estimates`, `Ystar_mat`, `FSsg_tab`, and all six are `identical()` to Stage A (`results` excludes `tmins_search`/`tmins_iteration` and any `timing` attribute). Full object names: results, SG_CIs, FSsg_tab, Ystar_mat, H_estimates, Hc_estimates, nb_boots, original_sg, outcome_type, effect_measure, est.scale, mr_replicates.

**62/62 Stage-C gates PASS.**

**Binary (i), Stage C against the old screened Stage A (i)** (report only, not a gate):

| check | result |
|---|---|
| `setequal()` of the `sg.harm` clauses | TRUE |
| sorted membership `identical()` | **TRUE** (as expected) |
| `sg.harm` clause order differs | TRUE: C `!{cd40 <= 380} & !{wtkg <= 86}` vs old `!{wtkg <= 86} & !{cd40 <= 380}` |
| candidate table equal after sorting by (HR, K, label) | **FALSE: 160 rows now vs 145 under the screen** |
| `filter_counts` | `n_passed_events` 2659 → 2666, `n_passed_sample_size` and `n_passed_cox` 2577 → 2584, `n_passed_hr` 145 → 160; `n_evaluated` 3240, `n_passed_redundancy` 2975 and earlier counts unchanged |
| `n_candidates_evaluated` / `_total` | 136 → 137 |
| `out_sg$result` | top rows have the same Pcons and HR; the `g`/`m` indices and the first-clause labels differ (column order) |

The "label" used for the candidate-table sort is the sorted set of active `q*` indicator columns. The table has no label column, and `grp` is dropped because it is a column-order index. The `q*` names come from `get_FSdata()` before the VI step, so they are stable across the two runs.

Reading: membership of the selected subgroup is unchanged, as expected. But ordering-only VI did more than reorder clauses. Through the redundancy column-walk it changed *which* combinations survived: the redundancy count is the same (2975), but different combinations make up that set, so 15 more candidates now pass the effect screen and the consistency pool grows by one. The 0.4.0 `@param vi.grf.min` text ("Measured effect on `n_candidates_evaluated` … zero") was therefore not exact on this fixture. The 0.4.0 docs only claimed that for the record's own measurement; this branch removes that text. See OPEN ITEMS.

## 4. Test results

Run with `devtools::load_all()` + `testthat::test_file()`, each under `timeout 10m`. Logs: `tests_gating_log.txt`, `tests_dependents_log.txt`, `tests_dependents_notcran_log.txt`, and per-file `dep_*.txt` / `depNC_*.txt`.

**Gating** — both pass, no skips:

```
TOTALS test-vi-deprecation: tests=4 failed=0 skipped=0 warnings=0 error=0
TOTALS test-fs-family-report: tests=5 failed=0 skipped=0 warnings=0 error=0
```

**Non-gating dependents**, pass 1 (as specified: `NOT_CRAN` unset, so `skip_on_cran()` applies; about 30 s in total):

```
TOTALS test-cross-outcome-parity: tests=9 failed=0 skipped=7 warnings=5 error=0
TOTALS test-return-shape-contract: tests=0 failed=0 skipped=0 warnings=0 error=0   (file-level skip_on_cran(), L22)
TOTALS test-sg-focus-transparency: tests=13 failed=0 skipped=3 warnings=0 error=0
TOTALS test-get-fsdata-diagnostics: tests=19 failed=0 skipped=0 warnings=4 error=0
TOTALS test-dfbeta-lm-influence-path: tests=5 failed=0 skipped=0 warnings=0 error=0
TOTALS test-downstream-consumers: tests=11 failed=0 skipped=4 warnings=1 error=0
```

Pass 2, the same six files with `NOT_CRAN=true` (same caps; 98 s in total; this deviation is noted in §6):

```
TOTALS test-cross-outcome-parity: tests=9 failed=0 skipped=2 warnings=12 error=0
TOTALS test-return-shape-contract: tests=10 failed=0 skipped=0 warnings=9 error=0
TOTALS test-sg-focus-transparency: tests=13 failed=0 skipped=0 warnings=0 error=0
TOTALS test-get-fsdata-diagnostics: tests=19 failed=0 skipped=0 warnings=4 error=0
TOTALS test-dfbeta-lm-influence-path: tests=5 failed=0 skipped=0 warnings=0 error=0
TOTALS test-downstream-consumers: tests=11 failed=0 skipped=0 warnings=5 error=0
```

There are no failures in either pass, and none of the test warnings is a deprecation warning: `grep deprecated` over every `dep*_*.txt` is empty. The two remaining parity skips are conditional on the data: `primary fit did not identify; test N/A` (test-cross-outcome-parity.R L262, L358). All six files ran inside the 45-minute cap.

## 5. §6 audit: explicit settings outside this task's edits (read-only, nothing edited)

Full table: `dev/reports/vi_deprecation_2026-10-26/audit_vi_settings.tsv`. It has 491 rows, one per line matching `(vi.grf.min|max_n_confounders) =` in `quarto/ dev/ inst/ R/ tests/ vignettes/ NEWS.md`. Rendered `.html`, this task's own battery and this task's own document are excluded. Columns are file, line, value, kind, and whether a rendered sibling `.html` exists and is committed.

By kind: 200 rows in application documents (`.qmd`/`.Rmd`; 118 files, 34 of them with a committed rendered `.html`), 164 in drivers/scripts, 108 in task/report records, 16 in tests (mostly this task's own), and 3 in package source (the two formals and the `fs_family_report` sprintf). `inst/doc/` does not exist. Largest groups: `quarto/applications` (85), `dev/applications_legacy` (71), `dev/glm-continuous-sims` (114 across drivers and records), `quarto/simulations` (27).

Entries Larry should look at first:

| file | line | value | kind |
|---|---|---|---|
| `quarto/applications/actg175/analysis_actg175_binary_multimethod_fixed_family.qmd` | 577, 1040 | `vi.grf.min = -0.2` | committed application document (rendered HTML committed) |
| `tests/testthat/test-fs-family-report.R` | 17 | `vi.grf.min = -0.2` (`.fr_driver_args()` driver mirror) | test (left as instructed) |
| `dev/applications_legacy/systolic-heart-failure/shf_noisey_EXTREME-K=100_VI=0.10_mins=15_benefit_nolasso_v0B.Rmd` | 218 | `vi.grf.min = 0.10` — **binding** | legacy application document (HTML committed) |
| `dev/applications_legacy/systolic-heart-failure/shf_noisey_EXTREME-K=100_VI=0.05_mins=15_benefit_nolasso_v0B.Rmd` | 218 | `vi.grf.min = 0.05` — **binding** | legacy application document (HTML committed) |
| `dev/applications_legacy/systolic-heart-failure/shf_noisey_EXTREME-K=100_VI-min=10_benefit_nolasso_v0B.Rmd` | 218 | `vi.grf.min = 0.05` — **binding** | legacy application document (HTML committed) |
| `dev/applications_legacy/systolic-heart-failure/.Rproj.user/…/A7D3F111-contents` | 218 | `vi.grf.min = 0.05` | RStudio session cache (not source) |

Values seen: `-0.2` (about 160 literal), `NULL` (about 85), indirections (`fs_vi_grf_min`, `vi_grf_min`, `grf_vi_min`, `cfg$fs_vi_grf_min`, `args$vi.grf.min`, …), and three binding values `> 0` (above, all legacy). After this change, every non-NULL call emits one warning and behaves as `NULL`.

**External repo.** `~/GitHub/fs-glms-interpretable/` does not exist. The sibling `~/Documents/GitHub/fs-glms-interpretable/` exists and was grepped read-only: 15 lines, all `max_n_confounders = 1000` (the default value, so the guard does not warn), and no `vi.grf.min`. List: `audit_vi_settings_fs-glms-interpretable.txt`.

Nothing in this audit was edited. Whether to clean the source lines (and re-render the records) is Larry's call.

## 6. Script deviations and OPEN ITEMS

**Deviations (verbatim):**

1. **`devtools::document()` also rewrote `DESCRIPTION`**:
   ```
   -Config/roxygen2/version: 8.1.0
   +Config/roxygen2/version: 8.1.1
   ```
   This is roxygen's own version stamp, not one of this task's edits. The task forbids `DESCRIPTION` changes, so it was restored with `git checkout -- DESCRIPTION` **before** `devtools::install()`. `DESCRIPTION` is byte-identical to the base. Expect roxygen 8.1.1 to re-stamp it on every future `document()` until Larry accepts the bump.
2. **Stage C run 1 compared an in-memory object against a reloaded one.** `stageC.R` run 1 compared each fresh in-memory Stage-C result against the Stage-A `.rds`. All 12 comparisons of the two `data.table` components (`out_sg_result`, `hr.subgroups`) returned FALSE, and every other component passed. Diagnosis, read-only: comparing the saved Stage-C `.rds` against the saved Stage-A `.rds` gives TRUE for all 12. An in-memory `data.table` compared with a `readRDS()` copy of itself is not `identical()`, because the `.internal.selfref` pointer is nil after `readRDS()`, while `all.equal(check.attributes = TRUE)` is TRUE. The gate is defined on *saved* components, so `stageC.R` was corrected to reload each Stage-C result from its own `.rds` before comparing (comment in the script) and re-run: 62/62 PASS. Run 1's record is kept as `stageC_summary_run1_selfref_artifact.txt` / `stageC_log_run1_selfref_artifact.txt`. The bootstrap gates already compared saved against saved and passed in both runs.
3. `devtools::install(upgrade = "never")` failed on its own argument check (`upgrade` must be TRUE, FALSE or NA) before doing anything. It was re-run with `upgrade = FALSE` and succeeded.
4. The dependents were run a second time with `NOT_CRAN=true`, because `skip_on_cran()` hid most of their tests in the specified run (§4). Both passes are reported.
5. In `fs_family_report.R` the block ends at L282, not L283 (§2).

**OPEN ITEMS for Larry:**

- `fit_causal_forest_glm()` (`R/grf_helpers.R` L117) has no callers now. It was left in place, as instructed.
- §5 audit: the committed ACTG175 binary fixed-family document (`vi.grf.min = -0.2`, L577/L1040) and the `.fr_driver_args()` mirror (L17) now warn and behave as `NULL`. Three legacy SHF documents pass **binding** values (0.05, 0.10) and would lose their screen if re-rendered.
- The binary (i) finding (§3): ordering-only VI changed the candidate pool (145 → 160 effect-passing candidates, 136 → 137 in the consistency pool) while leaving the selected membership unchanged. Any written record that says VI at `<= 0` is "ordering only, zero effect on `n_candidates_evaluated`" is not exact for the ACTG175 binary application.
- `Config/roxygen2/version` 8.1.0 vs the installed roxygen 8.1.1 (deviation 1): a one-line `DESCRIPTION` decision.
- The two data-conditional skips in `test-cross-outcome-parity.R` (L262, L358, "primary fit did not identify"). `.fs_args_for()` no longer passes `vi.grf.min = -0.2`, and it keeps `use_grf = TRUE`. Whether those two tests also skipped under 0.4.0 was not checked; they are in the full suite, Larry's pre-release gate.
- The installed package is the branch build `0.4.0.9000`, not CRAN 0.4.0.

**Verdict:**

1. The deprecation is implemented as specified. Both formals are kept, a non-default value warns once, and either value has no effect.
2. SECTION 5 (GRF VI screening) and the inert-cap warning are deleted, and no removed symbol is used afterwards.
3. Default callers are bit-identical to CRAN 0.4.0 on all three fixtures (survival, continuous, binary) and on the survival nb = 20 bootstrap payload.
4. A caller passing `vi.grf.min = -0.2` now gets exactly the default fit, with one warning and zero forest fits (CRAN fitted one).
5. Against the old screened fit, membership is identical and the clause order differs, as expected. The candidate pool also differs (see OPEN ITEMS).
6. `NAMESPACE` and `DESCRIPTION` are unchanged. `man/` regenerated `forestsearch.Rd` and `fs_family_report.Rd` only.
7. Both gating test files pass with no skips. All six dependents pass in both passes.
8. Two process deviations occurred (the roxygen stamp, the `data.table` selfref comparison); both were resolved without touching package code.
9. No R CMD check, no full suite, no vignette build, no fetch, pull, push or merge.
10. Ready for Larry's review and the lean-fit (A+F) task on this branch.

## 7. Closing block

Run at HEAD `65ed179e` (the Stage C + report commit). The closing-block commit that follows touches only this report and `postconditions.txt`, both inside the allowed list. Full output, including the complete `git diff --name-only` list: `dev/reports/vi_deprecation_2026-10-26/postconditions.txt`.

Commits: `4654b4cc` (task document) · `94738533` (Stage A) · `41347eb4` (edits) · `65ed179e` (Stage C, audit, report) · closing block (this commit).

```
HEAD: 65ed179e260480393ae11ecc8d27dc912c754084  branch: feature/vi-deprecation-lean-fits  base: a0025013f103d0921455102129caeb33127bb139
$ git diff --name-only <base> HEAD | grep -v -E <allowed>   # must print nothing
[end]
$ grep -n "variable_importance|causal_survival_forest|causal_forest|fit_causal_forest_glm|vi_ratio|vi_max|vi.cs\b|vi.order" R/forestsearch_main.R | grep -v comments   # must print nothing
[end]
$ grep -n "vi.grf.min|max_n_confounders" R/forestsearch_main.R | grep -v comments   # formals + guard only
1519:                         max_n_confounders = 1000,
1569:                         vi.grf.min = NULL,
1614:  if (!is.null(vi.grf.min)) {
1615:    warning("forestsearch(): 'vi.grf.min' is deprecated and ignored. GRF ",
1622:  if (!missing(max_n_confounders) &&
1623:      !identical(as.numeric(max_n_confounders), 1000)) {
1624:    warning("forestsearch(): 'max_n_confounders' is deprecated and ignored ",
[end]
$ grep -rn "vi.grf.min|max_n_confounders" R/ | grep -v forestsearch_main|fs_family_report   # must print nothing
[end]
$ git status --porcelain | grep -v '^?? '   # must print nothing
[end]
```

Every command printed what it should. The allowed-list filter, the removed-step grep, the other-`R/`-files grep and the dirty-tree check are all empty, and the deprecated-argument grep shows exactly the two formals and the guard's own lines. Neither `NAMESPACE` nor `DESCRIPTION` appears in the diff.
