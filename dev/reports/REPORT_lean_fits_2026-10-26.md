# REPORT — lean search-stage fits A (survival, `coxph.fit`) and F (binary OR, `glm.fit`)

**Task:** `dev/tasks/cc_task_lean_fits_2026-10-26.md` · **Executed:** 2026-10-26 by CC, unattended · **Record:** `dev/reports/lean_fits_2026-10-26/`

## 1. Provenance

- Branch `feature/vi-deprecation-lean-fits`; starting HEAD `60a6496bbe24f531519546004ede03d6bdc079d2` (as required). Starting `git status --porcelain`: `?? dev/tasks/cc_task_lean_fits_2026-10-26.md` (the copy just made) and `?? quarto/extreme_subgroups/fixed_random/` (pre-existing, never staged).
- MD5 gate: `R/subgroup_search.R` `f06d6bba8135f41a57962a659ebbe601`, `R/glm_effect_estimators.R` `8e8239f94e5f832af33e4d9cc0ff64d3` — PASS.
- Transport: `~/Downloads/cc_task_lean_fits_2026-10-26.md` (hyphens intact) → `dev/tasks/cc_task_lean_fits_2026-10-26.md`.
- Commits, in order:

| SHA | message |
|---|---|
| `679cdda1` | tasks: lean search-stage fits A (coxph.fit) and F (glm.fit), bit-identical |
| `d0b2fab4` | dev: Stage A timing baseline, lean fits |
| `7ff98f3d` | Lean survival search fit: coxph.fit() in fit_cox_for_subgroup() |
| `d9761a71` | Lean binary OR search fit: glm.fit closure and vector path |
| `a3e5d1d4` | Lean-fit identity tests; NEWS |
| `6cb98bc0` | dev: Stage C timing and bootstrap payloads, lean fits |
| (this report, then the closing block) | |

- Installed builds (`packageDescription("forestsearch")$Built`, all 0.4.0.9000, R 4.6.1, survival 3.8.12):
  - before (§1): `2026-10-26 17:05:25 UTC` — the deprecation-task build, as expected; no reinstall needed.
  - after A (§3): `2026-10-26 17:44:17 UTC`.
  - after F (§4): `2026-10-26 17:45:59 UTC`. **The installed package is now the A+F branch build.**
- `devtools::document()` was run after A and after F; each time it re-stamped `Config/roxygen2/version` in `DESCRIPTION` (restored with `git checkout -- DESCRIPTION` before installing) and changed nothing under `man/` or `NAMESPACE`.

## 2. The edits

```
 NEWS.md                         |  17 ++++
 R/forestsearch_main.R           |  29 +++++--
 R/glm_effect_estimators.R       | 111 +++++++++++++++++++++++++
 R/subgroup_search.R             |  93 ++++++++++++++++++++-
 tests/testthat/test-lean-fits.R | 177 ++++++++++++++++++++++++++++++++++++++++
 5 files changed, 415 insertions(+), 12 deletions(-)
```

Hunk headers (`git diff -U0 60a6496b HEAD`):

```
== R/subgroup_search.R
@@ -810,2 +810,3 @@ fit_cox_for_subgroup <- ...       comment phrase, re-wrapped (inside fit_cox_for_subgroup)
@@ -835 +836,34 @@ fit_cox_for_subgroup <- ...       coxph() -> coxph.fit() block (inside fit_cox_for_subgroup)
@@ -836,0 +871 @@ fit_cox_for_subgroup <- ...        names(beta) <- "Treat" (inside fit_cox_for_subgroup)
@@ -968,0 +1004,5 @@ fit_glm_for_subgroup <- ...       roxygen sentence of fit_glm_for_subgroup_fast (git labels the hunk with the previous top-level function)
@@ -981,0 +1022,6 @@ fit_glm_for_subgroup_fast <- ...  the fs_fast_core dispatch (inside fit_glm_for_subgroup_fast)
@@ -1004,0 +1051,41 @@ fit_glm_for_subgroup_fast <- ... new .fit_or_for_subgroup_vector() after it
== R/glm_effect_estimators.R
@@ -1100,0 +1101,111 @@ make_effect_estimator <- ...      new .fs_logit_two_group() + .make_glm_or_estimator_fast(), pure addition before the Utilities banner (git's context label is the last unindented top-level line it found)
== R/forestsearch_main.R
@@ -2091 +2091 @@ forestsearch <- ...               comment: "continuous/MD and binary/OR cases"
@@ -2101,2 +2101 @@ forestsearch <- ...
@@ -2104,6 +2103,18 @@ forestsearch <- ...             the routing predicate
```

Every hunk sits inside `fit_cox_for_subgroup()`, inside `fit_glm_for_subgroup_fast()` (body or roxygen), in a new function, or in the routing predicate/comment.

All line anchors were found where the task said: `R/subgroup_search.R` L810 / L834–L845 / L981–L982 / L1003–L1005 and `R/glm_effect_estimators.R` L1098–L1103 at their 0.4.0 numbers; the `R/forestsearch_main.R` block (anchored by content) now sits at L2091–L2104. None of them had moved.

Code was transplanted verbatim from the task, with these small layout choices:

- The edited comment line in `forestsearch_main.R` L2091 is now 86 columns. The task named a phrase swap with no re-wrap; I left it unwrapped.
- The roxygen sentence in `fit_glm_for_subgroup_fast()` was added as its own paragraph after the description paragraph, wrapped to 80 columns.

## 3. Identity

### Part A (§3) — `identA_summary.txt`

| check | result |
|---|---|
| `fit_cox_for_subgroup()` calls captured | **1410** (= Stage-0 count), 0 adjusted |
| per candidate, value `identical()` to the 0.4.0 formula block | **1410 / 1410** |
| per candidate, warning sets identical | **1410 / 1410** (no candidate warns; 0 NULL returns) |
| end to end, survival `"ii"` vs committed `stageC_survival_ii.rds` | sg.harm, membership, out_sg_result, hr.subgroups, filter_counts, n_candidates, n_fits, n_warnings, warnings: **all identical** (args_call_all also identical, report only) |

### Part F (§4) — `identF_summary.txt`

| check | result |
|---|---|
| `fit_glm_for_subgroup_fast()` calls captured | **2584** (see deviation §6.1 on the task's 2577) |
| slow `fit_glm_for_subgroup()` calls during the run | **0** |
| `attr(estimator_fn, "fs_fast_unadjusted")` / `fs_fast_core` | `TRUE` / `"logit_or"` |
| entry level, value `identical()` fast vs slow chain | **2584 / 2584** |
| entry level, warning sets identical | **2584 / 2584** (no candidate warns) |
| kinds (fast × slow) | fit/fit = 2584 (no NULL or nonestimable candidates in this application; those kinds are covered by the test file) |
| closure level (every 50th + every non-fit) | **51 / 51** identical (full list + warnings) |
| end to end, binary / continuous / survival `"ii"` vs committed `stageC_<fixture>_ii.rds` | all nine components **identical** for each fixture (args_call_all also identical, report only) |

### Bootstrap payloads (§6), `stageC_<fixture>_boot.rds` vs `stageA_<fixture>_boot.rds` (saved vs saved)

| fixture | components | results | SG_CIs | H_estimates | Hc_estimates | Ystar_mat | FSsg_tab |
|---|---|---|---|---|---|---|---|
| survival | same 6 | identical | identical | identical | identical | identical | identical |
| continuous | same 6 | identical | identical | identical | identical | identical | identical |
| binary | same 6 | identical | identical | identical | identical | identical | identical |

Stage A sanity gate: `stageA_survival_boot.rds` was identical, component by component, to the committed `vi_deprecation_2026-10-26/stageC_survival_boot.rds`. The machine reproduces its own record.

## 4. Timing

Sequential, the same three fixtures, `"ii"` args. Each cell has one warm-up, three timed calls (median taken) and one `nb_boots = 20` bootstrap. The replicate time is `mean(tmins_iteration) * 60`. Load average was about 0.4 at both starts and nothing else was running.

| fixture | call s (3 raw): before → after | call median s: before → after (%) | replicate mean s: before → after (%) | boot outer s: before → after (%) | Stage-0b predicted saving / replicate |
|---|---|---|---|---|---|
| survival | 3.06, 3.19, 3.04 → 1.56, 1.51, 1.64 | 3.06 → 1.56 (−49.1 %) | 3.60 → 2.02 (−43.9 %; −1.58 s) | 74.6 → 42.9 (−42.5 %) | 1.51 s |
| continuous | 1.48, 1.45, 1.61 → 1.57, 1.56, 1.45 | 1.48 → 1.56 (+5.3 %) | 1.87 → 1.84 (−1.5 %; −0.03 s) | 39.7 → 39.1 (−1.5 %) | 0 |
| binary | 5.52, 5.51, 5.34 → 2.20, 2.31, 2.15 | 5.51 → 2.20 (−60.1 %) | 5.93 → 2.75 (−53.6 %; −3.18 s) | 120.9 → 57.3 (−52.6 %) | 3.13 s |

**How the results compare with the predictions.** The absolute savings land on the Stage-0b predictions:

- **Survival:** saved 1.58 s per replicate against a predicted 1.51 s.
- **Binary:** saved 3.18 s per replicate against a predicted 3.13 s.
- **Continuous:** unchanged within run-to-run noise. Its call median rose 5 % while its replicate and bootstrap times fell 1.5 %; the three raw calls overlap.

The baseline replicate times here (3.60 s survival, 5.93 s binary) differ from Stage 0b's (3.14 s, 6.37 s). Stage 0b's binary figure was measured on the VI-screened candidate family, and a different measurement session explains the rest. As a result, the percentage savings (−44 % survival, −54 % binary) are near Stage 0b's projected halving, while the absolute seconds match almost exactly. **Both survival and binary clear Larry's 30 % bar.** No `below the 30% bar` FINDING was raised.

## 5. Tests (`NOT_CRAN=true`, `timeout 10m`, `devtools::load_all()` + `testthat::test_file()`)

Gating:
```
test-lean-fits.R            [ FAIL 0 | WARN 0 | SKIP 0 | PASS 344 ]   (3.2 s)
```
Non-gating:
```
test-fs-family-report.R     [ FAIL 0 | WARN 0 | SKIP 0 | PASS 105 ]
test-vi-deprecation.R       [ FAIL 0 | WARN 0 | SKIP 0 | PASS 15 ]
```
Logs: `tests_gating_log.txt`, `tests_nongating_log.txt`.

The test file follows the §5 battery: A on gbsg (25 random subsets, all rows, near-tie set, monotone likelihood — which asserts that the "Loglik converged before variable" warning occurs on both paths — n = 8, and single arm); F at both `adverse_outcome` settings at closure and entry level (25 random subsets, four zero-cell slices, `y == treat`, one event, one event per arm, 5 rows, single arm); and routing (unadjusted OR: 0 slow / > 0 fast; RD: > 0 slow / 0 fast; OR with `adjust_covariates = "age"`: 0 fast).

Additions beyond the literal list:
- The "one event among many non-events" item is covered twice: as written, which is non-estimable by the existence check, and as one event per arm, which reaches `glm.fit()` with extreme fitted probabilities.
- The routing test's call shape is transplanted from `test-binary-default-or.R`'s `.bd_fit()`.

## 6. Script deviations and OPEN ITEMS

**6.1 Part-F candidate count (verbatim).** Run 1 of `identF.R` printed:
```
INFO: captured 2584 fit_glm_for_subgroup_fast() calls (Stage-0b count 2577); slow fit_glm_for_subgroup() calls = 0
FAIL: captured candidate count equals the Stage-0b count 2577
```
Every other check in that run passed: all 2584 candidates were identical, warnings were identical, routing was correct, and all end-to-end comparisons matched. The count check was my own addition to the script. The task says "expected 2577 calls — the Stage-0b count; print it"; its §4 gates are identity, warning sets and routing.

The committed record resolves the difference:
- 2577 is `filter_counts$n_passed_cox` of `vi_deprecation_2026-10-26/stageA_binary_i.rds`. That is the VI-screened run (`vi.grf.min = -0.2`, CRAN 0.4.0) — the config Stage 0b measured.
- The `"ii"` family has `n_passed_cox = 2584` on both CRAN 0.4.0 (`stageA_binary_ii.rds`) and the deprecation build (`stageC_binary_ii.rds`). It is the deprecation-task change in the candidate family (145 → 160 passing the effect screen), not this task.

What I did:
- Preserved run 1 as `identF_log_run1_count_premise.txt` / `identF_summary_run1_count_premise.txt`.
- Re-anchored the count check to the committed `"ii"` record's `n_passed_cox`, with the Stage-0b number printed as INFO.
- Reran. All gates passed.

The Part-A count (1410) matched its stated expectation exactly.

**6.2 Order of §5 and §6.** The §5 NEWS entry must carry the §6 measured numbers, so Stage C (§6) ran before the §5 commit.
- No code or installed build changed between the test runs and Stage C.
- The test logs and the §5 outcomes copied into `stageC_summary.txt` are from runs made before Stage C.
- Commits: `a3e5d1d4` (tests + NEWS), then `6cb98bc0` (Stage C).

**6.3 Record layout.**
- `timing_common.R` holds the §2/§6 protocol shared by `stageA.R` and `stageC.R`, so both stages run identical timing code. It sources `vi_common.R` unchanged.
- The end-to-end comparisons in `identA.R` / `identF.R` use a separate, untraced `run_collect()` call; the traced run is used only for the per-candidate capture.

**OPEN ITEMS**
1. `run_subgroup_sims.R` still routes only the continuous case to a fast closure. Giving binary-OR simulations the same gain is a one-predicate follow-up, if Larry wants it.
2. The `DESCRIPTION` `Config/roxygen2/version` stamp (roxygen 8.1.1 re-stamps it on every `document()`; restored each time) is still Larry's call.

**Verdict.**
1. Both lean fits are in and bit-identical on every measured surface.
2. A: 1410/1410 survival candidates identical in value and warnings.
3. F: 2584/2584 binary candidates identical at entry level, 51/51 at closure level, with zero slow-path calls.
4. End-to-end runs of all three fixtures match the committed deprecation record, warnings included.
5. All three `nb = 20` bootstrap payloads are identical before and after.
6. Gating tests: 344/344. Non-gating: 105/105 and 15/15.
7. Per-replicate time: survival −44 % (−1.58 s), binary −54 % (−3.18 s), continuous unchanged.
8. The savings match Stage 0b's predictions to within 0.07 s; both clear the 30 % bar.
9. No `man/`, `NAMESPACE` or `DESCRIPTION` change.
10. One script deviation (§6.1): a count expectation carried over from the pre-deprecation candidate family. It is documented, not a defect.

## 7. Closing block

(appended after the report commit)
