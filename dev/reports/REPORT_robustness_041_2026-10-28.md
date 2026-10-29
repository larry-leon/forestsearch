# REPORT — robustness for 0.4.1: survival floor, continuous IPTW weights fix, lean Cox core with probe

**Task:** `dev/tasks/cc_task_robustness_041_2026-10-28.md` · **Executed:** 2026-10-29 by CC, unattended
**Branch:** `feature/robustness-0.4.1` from `master` at `606a78ea` · **Record directory:** `dev/reports/robustness_041_2026-10-28/`

---

## 1. Provenance

- Starting HEAD: `606a78ea55b0ec772ae0c99e1462b866d95731f1` on `master` (= the audit basis); `git status --porcelain` empty (no untracked files either).
- §1 MD5 gate: all nine files matched the pinned hashes.
- `packageVersion("survival")`: **3.8.12** (≥ 3.5-8, gate passed). `R version 4.6.1 (2026-06-24)`.
- Installed `forestsearch` before the task: 0.4.0.9000, Built `R 4.6.1; ; 2026-10-27 21:01:28 UTC; unix`.
- Installed for Stage A (§2, = `master`): 0.4.0.9000, Built `2026-10-29 14:19:31 UTC`.
- Installed after the edits (§6, the branch build — **this is the build installed at the end of the task**): 0.4.0.9000, Built `2026-10-29 14:24:17 UTC`; installed `Imports` carries `survival (>= 3.5-8)` (checked: `TRUE`).
- Task document: `~/Downloads/cc_task_robustness_041_2026-10-28.md` → `dev/tasks/cc_task_robustness_041_2026-10-28.md` (same filename).

Commits, in order:

| SHA | Message |
|---|---|
| `a3c36b6f` | Task document: robustness for 0.4.1 |
| `ff6cd9e7` | Stage A baseline on master, robustness task |
| `255361ad` | Reproduce the continuous IPTW weights defect on master |
| `f6735639` | Fix continuous IPTW weights resolution in .make_lm_estimator() |
| `11937932` | Lean Cox core .fs_cox_two_group() with once-per-session probe and formula fallback |
| `fe2b3342` | DESCRIPTION: survival (>= 3.5-8) |
| `ddaacbd2` | Stage C: identity against Stage A, robustness task |
| `342b1742` | NEWS: survival floor, lean-core probe, continuous IPTW fix |
| (this commit) | Report: robustness for 0.4.1 |

## 2. The edits

```
$ git diff --stat 606a78ea HEAD -- R/ DESCRIPTION tests/ NEWS.md
 DESCRIPTION                           |   2 +-
 NEWS.md                               |  18 ++++
 R/forestsearch_main.R                 |   2 +-
 R/glm_effect_estimators.R             |   8 +-
 R/globals.R                           |   8 +-
 R/subgroup_search.R                   | 160 ++++++++++++++++++++++++++--------
 tests/testthat/test-iptw-continuous.R |  47 ++++++++++
 tests/testthat/test-lean-fits.R       |  36 ++++++++
 8 files changed, 243 insertions(+), 38 deletions(-)
```

Hunk headers of the `R/` files:

- `R/glm_effect_estimators.R`: `@@ -954,7 +954,13 @@ make_effect_estimator <- function(` — the hunk (new lines 954–966) sits inside `.make_lm_estimator()` (L929–L989); git labels it with `make_effect_estimator` only because its default function-name heuristic skips definitions starting with `.`. It is the closure's IPTW branch: six comment lines and `weights = data_slice$sw` → `weights = sw`.
- `R/globals.R`: `@@ -263,7 +263,13 @@ utils::globalVariables(c(` — inside the `globalVariables` vector (adds `"sw"` with its section header).
- `R/forestsearch_main.R`: `@@ -2573,7 +2573,7 @@ forestsearch <- function(df.analysis,` — the PS-block comment line only.
- `R/subgroup_search.R`: `@@ -798,6 +798,121 @@ meets_event_criteria <- ...` — the new helpers (`.fs_cox_lean_state`, `.fs_cox_lean_conventions()`, `.fs_cox_two_group()`, `.fs_cox_lean_ok()`) inserted between `meets_event_criteria()` and the `fit_cox_for_subgroup()` roxygen header; `@@ -843,40 +958,17 @@ fit_cox_for_subgroup <- ...` — the unadjusted branch's `try({ ... })` block only.

Every hunk sits in one of the permitted places. All §4/§5 anchors were at the line numbers the task gives (L956–957, L266–267, L2576, L799/L801, L845–L890); none moved. `devtools::document()` produced no `man/`, `NAMESPACE` or `DESCRIPTION` change (not even the roxygen re-stamp). `DESCRIPTION` numstat `1 1`.

§6.4 (informational): survival's NEWS database is readable (printing it errors inside `format.news_db`, so it was read via `as.data.frame()`). The earliest entry mentioning `nocenter` is **3.2-9** ("The coxph routine no longer recenters 0/1 predictors ... per the default value of the new nocenter argument"); the earliest mentioning the `coxph.fit` residual argument is **2.41-4** ("Add the residual and concordance argument to coxph.fit and agreg.fit"). The floor stays **3.5-8**, the lower end of the verified range.

## 3. Reproduce (§3 and §8 item 4)

Before the fix (`reproduce_iptw_log.txt`, Stage-A build):

```
closure: estimate = NA, se = NA, converged = FALSE, method_used = ols_failed
PASS: closure-level weighted fit returns NA (the defect is present)
end to end: sw in df_clean = TRUE; fitted candidates = 0; sg.harm =
```

After the fix (`stageC_summary.txt`, branch build):

```
INFO: [4] closure after the fix: estimate = 0.8358334, se = 0.1334824, converged = TRUE, method_used = ols
PASS: [4] closure-level weighted fit now returns an estimate
```

End to end on the branch build, `test-iptw-continuous.R` test 2 asserts a non-empty `hr.subgroups` for the same call — passed.

## 4. Identity

**Stage C vs Stage A (GATE), saved against saved, `identical()`:**

| Fixture | End to end (10 components: `sg.harm`, `membership`, `out_sg_result`, `hr.subgroups`, `filter_counts`, `n_candidates`, `n_fits`, `n_warnings`, `warnings`, `args_call_all`) | Bootstrap (`names_found` + 6 components: `results`, `SG_CIs`, `H_estimates`, `Hc_estimates`, `Ystar_mat`, `FSsg_tab`) |
|---|---|---|
| survival (gbsg) | 10/10 PASS | 7/7 PASS |
| continuous (ACTG175 MD) | 10/10 PASS | 7/7 PASS |
| binary (ACTG175 OR) | 10/10 PASS | 7/7 PASS |

**Stage A vs committed lean-fits records (§2.3, INFO):** all 48 lines `identical` — every end-to-end component against `../lean_fits_2026-10-26/identF_<fixture>_ii.rds` and every bootstrap component against `../lean_fits_2026-10-26/stageC_<fixture>_boot.rds`, for all three fixtures. No component differs, including `args_call_all` (the `max.minutes` deprecation did not change what the echo carries for these fixtures).

## 5. Tests

Gating (`FAIL 0`, `PASS > 0` required):

```
test-lean-fits.R        [ FAIL 0 | WARN 0 | SKIP 0 | PASS 362 ]
test-iptw-continuous.R  [ FAIL 0 | WARN 0 | SKIP 0 | PASS 7 ]
```

Non-gating:

```
test-fs-family-report.R [ FAIL 0 | WARN 0 | SKIP 0 | PASS 105 ]
test-vi-deprecation.R   [ FAIL 0 | WARN 0 | SKIP 0 | PASS 15 ]
```

`test_file(..., package = "forestsearch")` loaded the helpers on this testthat (the IPTW end-to-end test uses `.make_continuous_data()` from `helper-synthetic-dgm.R` and ran, SKIP 0); the fallback invocation was not needed.

## 6. Timing (informational)

Sequential, one warm-up + three timed calls + one `nb_boots = 20` bootstrap per fixture.

| Fixture | Call median s A → C | Replicate mean s A → C | Boot outer s A → C |
|---|---|---|---|
| survival | 1.53 → 1.56 (+2.0%) | 2.07 → 2.07 (−0.2%) | 44.0 → 43.8 (−0.4%) |
| continuous | 1.64 → 1.48 (−9.9%) | 1.83 → 1.83 (−0.1%) | 38.9 → 38.9 (−0.1%) |
| binary | 2.36 → 2.26 (−4.4%) | 2.82 → 2.82 (−0.1%) | 58.8 → 58.7 (−0.1%) |

The replicate means agree to 0.2%; the call-median movements are three-call noise (the continuous fixture's Stage A median includes one 1.88 s outlier; the continuous path is not touched by any edit). The probe adds one 30-row fit per session, not visible at this resolution.

Probe after the Stage C run: `.fs_cox_lean_state$ok = TRUE`; conventions read from `coxph()`'s formals on survival 3.8.12: `ties = "efron"`, `nocenter = c(-1, 0, 1)` — the literal values the 0.4.0.9000 lean block hard-coded.

## 7. Behaviour statement for Larry's record

1. **Floor.** `survival (>= 3.5-8)` changes no code; it only refuses installation against an older survival. This machine has 3.8.12.
2. **IPTW fix.** The changed statement is reachable only with `is.RCT = FALSE` (or `ps_method` supplied), `ps_adjust_method = "iptw"` and a continuous outcome. The repository grep for `ps_adjust_method = "iptw"` over `quarto/applications`, `vignettes`, `tests`, `inst` finds only the new regression test `tests/testthat/test-iptw-continuous.R` (closing block): no committed analysis uses the path; all are `is.RCT = TRUE`. The `else` branch (unweighted `lm()`) is byte-identical, the fast continuous closure is untouched, and the survival/OR/count IPTW sites are untouched (three non-comment `weights = data_slice$sw` lines remain). Before the fix the path returned NA for every fit and an empty search; after it the closure equals `lm(weights = ...)` exactly (`expect_identical` on estimate and SE) and the search yields candidates.
3. **Lean core.** On this survival the probe passes and the lean core executes the same `coxph.fit()` call with the same argument values as the 0.4.0.9000 block (`ties`/`nocenter` now read from the formals and equal to the former literals). Stage C equals Stage A on every end-to-end component of every fixture and on every bootstrap component; the forced-fallback test shows the formula branch reproduces the 0.4.0 reference as well. The only behavioural difference is in the failure case (a survival whose `coxph.fit()` no longer reproduces `coxph()`), which previously produced a silent empty search and now falls back to the formula interface with one message.

## 8. Script deviations and open items

**Script deviations**

1. Stage A and Stage C were run as `timeout 15m Rscript stage?.R > stage?_log.txt 2>&1` rather than `2>&1 | tee stage?_log.txt`; the log content is the same (console not echoed live).
2. §6.4: `print()` of the `news()` result errors in `format.news_db` (`values must be length 1`); the database was read via `as.data.frame()` instead. Findings in §2 above.
3. **§11 expected values corrected (accepted by Larry, 2026-10-29).** Three of the task's post-condition checks printed values other than the document's expectations, because text the task itself prescribed matches the patterns. The committed task document is not edited. Original lines and outputs:
   ```
   grep -c 'weights = data_slice\$sw' R/glm_effect_estimators.R
       # must print 3          -> printed 4 (the fourth is the prescribed §4 comment line 959:
       #                           "# whose frame does not hold `data_slice`, so `weights = data_slice$sw`")
   grep -n "coxph.fit(" R/subgroup_search.R | grep -v -E '^[0-9]+:[[:space:]]*#' | wc -l
       # must print 1          -> printed 2 (the second is the prescribed §5 probe message string, line 906:
       #                           message("forestsearch: survival::coxph.fit() did not reproduce coxph() ", ...)
   grep -rn 'ps_adjust_method *= *"iptw"' --include=*.qmd --include=*.Rmd --include=*.R quarto/applications vignettes tests inst 2>/dev/null
       # expected empty        -> printed two lines, both in the new tests/testthat/test-iptw-continuous.R (L16, L43)
   ```
   CC stopped at the gate and reported; Larry accepted these as defects in the expected values and supplied the corrected checks, which replace the three originals in the closing block:
   ```
   grep -n 'weights = data_slice\$sw' R/glm_effect_estimators.R | grep -v -E '^[0-9]+:[[:space:]]*#' | wc -l
       # must print 3 (non-comment lines only: survival, OR, count sites)
   grep -n "coxph.fit(" R/subgroup_search.R | grep -v -E '^[0-9]+:[[:space:]]*#' | grep -v 'message(' | wc -l
       # must print 1 (the call inside .fs_cox_two_group(); the probe's message string is the other match)
   grep -rln 'ps_adjust_method *= *"iptw"' --include=*.qmd --include=*.Rmd --include=*.R quarto/applications vignettes tests inst 2>/dev/null
       # must print exactly one file: tests/testthat/test-iptw-continuous.R
   ```

**OPEN ITEMS (Larry's decisions; nothing changed here)**

- (a) The consistency-stage pieces (`.consistency_glm_pieces()`, `R/consistency_resample.R:239–`) take no weights for any `ps_adjust_method`, so an IPTW analysis is weighted in the search stage and unweighted in the consistency screen. Now that the continuous IPTW search stage actually fits, this asymmetry is live for continuous outcomes too — a method question.
- (b) At `is.RCT = FALSE` with the defaults (`ps_method = "grf"`, `ps_adjust_method = "none"`) a propensity forest is fitted once per call and once per bootstrap replicate and never consumed, and the rebuilt estimator closure is untagged, so the lean continuous/OR path is lost (`R/forestsearch_main.R:2626–2697`).
- (c) `subgroup_glm()` in `R/run_subgroup_sims.R:331–340` still routes only the continuous case to a fast closure.
- (d) The §11 grep patterns in future task documents should exclude comment and string matches of the prescribed text (deviation 3).

**Verdict**

1. All §1 provenance gates passed at the audit basis `606a78ea`.
2. Stage A on `master` reproduces the committed lean-fits records on every component (48/48).
3. The IPTW defect was reproduced on the unmodified build before any edit (closure NA; 0 candidates; `sw` present in the frame).
4. The one-statement fix makes the weighted closure equal `lm(weights = ...)` exactly and the end-to-end path yields candidates.
5. The lean Cox core moved into `.fs_cox_two_group()` behind a once-per-session probe; the probe passes on survival 3.8.12.
6. Stage C equals Stage A on all 30 end-to-end and 21 bootstrap checks across the three fixtures.
7. Gating tests: 362 + 7 passes, 0 failures; non-gating: 105 + 15 passes, 0 failures.
8. Timing is unchanged within noise (replicate means within 0.2%).
9. No `man/`, `NAMESPACE`, or version change; `DESCRIPTION` changes only the survival floor.
10. The §11 gate stopped on three mis-specified expected values; with Larry's corrected checks, every post-condition holds. The branch is ready for Larry's `R CMD check --as-cran`.

## 9. Closing block

Run after the report commit; see the record of the run below.

CLOSING_BLOCK_PLACEHOLDER
