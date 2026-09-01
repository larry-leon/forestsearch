# REPORT — survival medians on the passing set: the `survfit` computation moved behind the effect screen (0.3.4)

**Task:** `dev/tasks/cc_task_medians_on_survivors_2026-09-02.md` · **Executed:** 2026-09-02 (session date 2026-09-01) by CC
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`

**Verdict up front: all gates pass.** The medians computation moved bit-identically: every retained component of every fixture and the 20-replicate gbsg bootstrap payload is `identical()` across 0.3.3 → 0.3.4, `survfit` dispatches on gbsg fell 1,410 → 121 (exactly the passer count, zero residual calls), the F2 medians bucket fell 2.83 s / 37.8 % → 0.22 s / 4.6 %, and the gbsg bootstrap replicate fell 7.47 → 5.49 s (B = 1000: 124.5 → 91.5 min).

---

## 1. Provenance and §2 re-verification

**Machine/repo:** pop-os, `/home/larryleon/Documents/GitHub/forestsearch`, branch `feature/glm-extension`, HEAD at task start `5d8f7f68`, `git status --porcelain` clean. `packageVersion("forestsearch")` **0.3.3** ✓. R 4.6.1 (2026-06-24), 128 cores. *GATE:* branch ✓; the four fast-path commits `c79f797d`, `b9d7b622`, `4dbf9f26`, `5d8f7f68` all in `git log` ✓. `vi.grf.min` default in force: **`NULL`** (`R/forestsearch_main.R:1249`) ✓.

**Task doc committed alone:** `~/Downloads/cc_task_medians_on_survivors_2026-09-02.md` → `dev/tasks/cc_task_medians_on_survivors_2026-09-02.md`, commit `ff2f9ca2`.

### §2.1 `fit_cox_for_subgroup()` at HEAD (`5d8f7f68`, `R/subgroup_search.R:730–815`)

Signature `fit_cox_for_subgroup(yy, dd, tt, id.x, df_clean = NULL, adjust_covariates = NULL)`. The `df.x` assembly (L736–749):

```r
data.x <- data.table::data.table(Y = yy, E = dd, Treat = tt, id.x = id.x)   # L736
...                                       # adjusted arm cbinds raw adj columns, L740-747
df.x <- data.x[id.x == 1]                                                   # L749
```

The 0.3.3 unadjusted fast branch (L753–782, direct `conf.int` arithmetic) and the adjusted branch (L783–788, `summary(coxph(...))$conf.int`) both sit inside `try(..., silent = TRUE)` with a NULL-on-error contract (L790). **The medians computation** — after the arm branch, so common to both arms (L800–806):

```r
  # Get median survival times (descriptive; treatment-only by design)
  meds <- try(
    summary(survival::survfit(survival::Surv(Y, E) ~ Treat, data = df.x))$table[, "median"],
    silent = TRUE
  )

  if (inherits(meds, "try-error")) return(NULL)
```

and the extraction into the return (L812–813): `med0 = meds[1], med1 = meds[2]`. Mechanics that decided the §4 frame rebuild: the `survfit` formula is `survival::Surv(Y, E) ~ Treat` with `data = df.x`, no other arguments; `$table` is indexed by **column name `"median"`** and the two rows are taken **positionally** (`meds[1]` = `"Treat=0"`, `meds[2]` = `"Treat=1"` — row order and the names carried on `meds` come from the `Treat` strata labels). The computation consumes **only the `Y`, `E`, `Treat` columns of `df.x`** (the adjusted arm's extra columns are not in the formula) — so the rebuild must reproduce exactly those names, and does.

### §2.2 The per-candidate sequence (`evaluate_combination_with_status`, survival path L650–682)

Cox fit at L665; `NULL` → `return(list(status = 5L, result = NULL))` (L668–670). **The effect screen is per candidate, in-loop** (L672–677):

```r
  # Status 6: Check HR threshold. ...
  if (!disable_effect_floor && cox_result$hr <= hr.threshold) {
    return(list(status = 6L, result = NULL))
  }

  # Status 7: Passed all criteria - return result
  result_row <- create_result_row(kk, covs.in, nx, event_counts, cox_result)  # L680
  return(list(status = 7L, result = result_row))
```

Status codes (roxygen L540–548): 5 = passed sample size, failed model fit; 6 = passed model fit, failed effect threshold; 7 = passed all. *GATE:* the verdict is known per candidate inside the iteration, before row creation, with `yy`, `dd`, `tt`, `id.x` all in scope (they are formals of the function) — the §4 design applies as written. ✓

### §2.3 Row provenance

`result` is non-NULL only on the status-7 return; the collection loop (`R/subgroup_search.R:290–304`) appends to `results_list` only `if (!is.null(res$result))`, and `format_search_results()` builds `hr.subgroups` from exactly that list (L964–981). **Rows reaching `hr.subgroups` exist only for effect-screen passers.** ✓ The record's counts reproduced at Stage A: gbsg 1,410 fits → 121 rows; continuous F1 1,975 → 852; anchor 4,702 → 146.

### §2.4 Callers of `fit_cox_for_subgroup()`

Package-wide grep (R/ and tests/): exactly one call site — `R/subgroup_search.R:665`, inside the per-candidate iteration; it consumes `hr` (screen) and passes the whole list to `create_result_row()` (`m1`/`m0` enter the row at L924–925). *GATE:* no caller outside the iteration consumes `med0`/`med1`. ✓

### §2.5 The two consumer sites at HEAD

Dedup key — `R/subgroup_consistency_helpers.R:1108–1122` (`remove_near_duplicate_subgroups`):

```r
  # Columns to check: K, n, E, d1, m1, m0, HR, L(HR), U(HR)
  cols_to_check <- 2:min(10, ncol(df))                       # L1109
  ...
  key_cols <- df_rounded[, cols_to_check, drop = FALSE]      # L1118
  dup_key <- apply(key_cols, 1, function(x) paste(x, collapse = "_"))
```

`m1.threshold` filter — `R/subgroup_consistency_main.R:531–546`:

```r
  if (is.finite(m1.threshold)) {                             # L531
    hr.subgroups <- hr.subgroups[!is.na(hr.subgroups$m1), ]  # L532
    ...
    found.hrs <- hr.subgroups[hr.subgroups$HR >= hr.threshold &
                                hr.subgroups$m1 <= m1.threshold, ]   # L545-546
```

Formal defaults `Inf` at `subgroup.consistency()` L357 and `forestsearch()` `R/forestsearch_main.R:1228` — the two sites the 09-01 report recorded; `forestsearch()`'s default flows to `subgroup.consistency` through `filter_call_args(args_call_all, ...)`, and a package-wide grep finds **no call site setting it finite** (the 09-01 report's "two live call sites passing `Inf`" are these two formal-default paths). Both blocks unchanged in substance from the 09-01 quotes. ✓ These files are untouched by this task (`git diff` at close: only `R/subgroup_search.R`, `DESCRIPTION`, `NEWS.md`).

### §2.6 No-other-reader survey at HEAD

Grep of every `m1`/`m0`/`med0`/`med1` read in R/, classified:

- **Consumer-of-`hr.subgroups`:** the dedup key (helpers L1108–1118); the `m1.threshold` filter (main L531–546); the column-name bookkeeping `subgroup_consistency_main.R:482`.
- **Display-for-selected-subgroup:** `summary_utility_functions.R` (L128–139, 445–449, 556–557, 621–635 — computes its **own** medians via `survival::survfit` at L74 for the selected subgroup's tables); `bootstrap_summaries_helpers.R:99–167` (gt column labels); `forestsearch_cross_validation.R:1702–1718` (table column names).
- **Other/unrelated:** `oc_analyses.R:300–306` (`m1`/`m0` are local arm means of a GLM outcome, nothing to do with medians); `guohe_adaptive_r.R:41` (comment); `fs_family_report.R` (reports the `m1.threshold` *setting*, reads no candidate medians).

*GATE:* nothing reads medians for candidates outside `hr.subgroups`. ✓

### §2.7 Fixture arguments

Reused **verbatim** from `dev/glm-continuous-sims/fastpath_baseline_2026-09-02.R`: the F1/F2/F3/F5 argument blocks, the E4 (heavy-ties) and E5 (monotone-likelihood) constructions, the pruner and exclusion list, the instrumented-runner architecture.

---

## 2. The move as implemented (commit `3921ffdd`)

Files changed: `R/subgroup_search.R`, `DESCRIPTION` (0.3.3 → **0.3.4**), `NEWS.md`. Nothing else; `devtools::document()` produced no NAMESPACE/man changes (the touched functions are internal with unchanged signatures).

**Removed** from `fit_cox_for_subgroup()` (pre-move L800–806 + the `meds[1]`/`meds[2]` extraction at L812–813): the seven medians lines quoted in §2.1 above, and the return now carries `med0 = NA_real_, med1 = NA_real_` (post-move `R/subgroup_search.R:823–833`). The rest of the body — `df.x` assembly, the 0.3.3 fast branch, the adjusted branch, both `try()` contracts — is untouched apart from one stale comment clause inside the fast-branch comment block: "the medians computation below consumes it" → "the coxph() fit consumes it" (the `df.x` assembly is still consumed, by the fit itself; recorded here since it is a line inside the otherwise-untouched body). `fit_cox_for_subgroup()`'s roxygen never mentioned medians; no roxygen change was needed.

**Added** at the post-screen point of the survival path in `evaluate_combination_with_status()` (post-move `R/subgroup_search.R:679–704`):

```r
  idx <- id.x == 1
  med_df <- data.frame(Y = yy[idx], E = dd[idx], Treat = tt[idx])
  meds <- try(
    summary(survival::survfit(survival::Surv(Y, E) ~ Treat, data = med_df))$table[, "median"],
    silent = TRUE
  )
  # Pre-move contract, replicated: a survfit failure inside
  # fit_cox_for_subgroup() made it return NULL, i.e. status 5 (failed
  # model fit).
  if (inherits(meds, "try-error")) {
    return(list(status = 5L, result = NULL))
  }
  cox_result$med0 <- meds[1]
  cox_result$med1 <- meds[2]
```

then the unchanged `create_result_row()` call. The `survfit` call and extraction are the removed lines operation for operation; the rebuilt frame carries the same `Y`/`E`/`Treat` names so the `"Treat=0"`/`"Treat=1"` strata labels, row order, and the names on `meds` are identical, and the patched `cox_result$med0/med1` are the same named scalars that entered the row before. Non-passers return at status ≤ 6 before this block — no medians call, and (per §2.3) no row. One recorded theoretical asymmetry: a candidate whose `coxph` succeeds, whose `survfit` errors, **and** which would fail the screen was status 5 at 0.3.3 and is status 6 at 0.3.4 (the screen now runs first) — no such candidate exists in any fixture (`survfit` does not error where `coxph` fit the same frame), and the equality gates cover the counters that would expose one.

The GLM path (continuous/binary) returns at L647, before the block — bitwise untouched, as gated.

---

## 3. Stage A / Stage C mechanics and the equality matrix

**Scripts:** `medians_baseline_2026-09-02.R` (one script, both stages via `FS_MEDIANS_OUT`; Stage A at `4c923dfd` on 0.3.3, Stage C on installed 0.3.4), `medians_compare_2026-09-02.R` (the gates), logs and `.rds` artefacts alongside. **The volatile-field exclusion list, reused verbatim** from `fastpath_baseline_2026-09-02.R`: list fields `time_search` (find.grps) and `minutes_all` (top level); bootstrap frame columns `tmins_search`, `tmins_iteration`, `tmins_all`; data.tables normalised to data.frame (`.internal.selfref`); environments/functions/calls/formulas replaced by tags, every dropped path recorded in the rds (26 unique paths, all of the listed kinds). Self-consistency (F2, E-ties, bootstrap run twice per stage): 3/3 identical at Stage A **and** at Stage C.

**Fixtures and findings.**

- **E-finite value chosen: `m1.threshold = 60` months** on the F2 gbsg configuration. Verified at Stage A: of the 121 screen-passing rows, it **excludes 5** (finite `m1` > 60 or `m1` NA among passers) and **keeps 116**; the filter path (`is.finite` branch, NA-drop, `m1 <=` comparison) runs live. Selection under the filter: `{er <= 0} & {size <= 35}` (same subgroup as F2 here), search-stage frames at 116 post-filter rows.
- **E-named:** engineered control arm never reaching median survival inside a passing candidate (z1 == 1: controls 6 events / 40, KM floor ≈ 0.85; treated 36/40 events overlapping in time — no monotone likelihood). Stage A check "NA median in a passer" OK; the NA flows into the dedup key string and reproduced identically.
- **E-zero:** E-ties data with `hr.threshold = 100` — the screen passes nobody, `hr.subgroups` absent, no-find return, no error, warning count 0 both stages, and **0 `survfit` calls at 0.3.4** (the no-op elision).
- **E-degen:** the fast-path E5 monotone-likelihood construction **passes the screen as-built** (diverging HR > 100 verified in Stage A `hr.subgroups`; no threshold adjustment was needed). 2 warnings (coxph "Loglik converged before variable") at parity both stages; medians computed for it identically.
- **E-ties:** heavy-ties survfit reproduced identically (2 passers, moved calls match).

**Equality matrix (`medians_compare_2026-09-02.log`):** every gate × case × component `identical() == TRUE` — for each of F2, F5 (E-1), F1, F3 (E-2), Eties, Enamed, Efinite, Ezero, Edegen (E-3): `pruned` (the full pruned `forestsearch` object — `hr.subgroups`/`out.found` at full precision including `m1`/`m0`, selected subgroup, consistency outputs), `warnings`, `counters` (fit/consistency/dedup), `sg.harm`; plus the E-4 bootstrap payload. **37/37 TRUE. ALL EQUALITY GATES PASS.**

---

## 4. Elision proof — the `survfit` count (§6)

`trace()` counter on the `survfit` generic in the survival namespace, innermost forestsearch caller recorded per dispatch, per single fixture call (battery instrumentation, both stages):

| Fixture | 0.3.3 count (attribution) | 0.3.4 count (attribution) | passers |
|---|---|---|---|
| F2 gbsg | 1,410 — all `fit_cox_for_subgroup` | **121** — all `evaluate_combination_with_status` | 121 |
| F5 gbsg adjusted | 1,410 — all `fit_cox_for_subgroup` | **111** — all `evaluate_combination_with_status` | 111 |
| E-finite | 1,410 — all `fit_cox_for_subgroup` | 121 — all `evaluate_combination_with_status` | 121 |
| E-ties / E-named / E-degen | 4 / 5 / 4 | 2 / 3 / 4 | 2 / 3 / 4 |
| E-zero | 4 | **0** | 0 |
| F1 / F3 (continuous) | **0 / 0** | **0 / 0** | — |

After-count = passer count **exactly** on every fixture; there are **no** counts above the passer set to attribute — the display sites (`summary_utility_functions.R:74` etc.) are not exercised under this configuration (`plot.sg = FALSE`, quiet), consistent with the before-counts equalling the fit counts exactly. Continuous guards: zero `survfit` calls at both versions. Before-ratio on F2: 1,410/121 = **11.7×**, the recorded measurement.

---

## 5. Realized recovery (§7)

**F2 single call, `Rprof(interval = 0.010, line.profiling = TRUE)`** (`medians_profile_2026-09-02.R`, established buckets verbatim; fresh before-profile at Stage A):

- Medians bucket: **2.83 s / 37.8 %** (0.3.3; record: 2.69 s / 35.5 % at 0.3.2) → **0.22 s / 4.6 %** (0.3.4). Expected ≈ 0.2–0.3 s: on target.
- Profiled gbsg wall: **7.50 s → 4.82 s** (record's expectation ≈ 7.4–7.6 → ≈ 5 s ✓); battery walls F2 8.10 → 5.38 s, E-finite 7.90 → 4.91 s; F5 adjusted 40.42 → 38.14 s (its wall is dominated by the adjusted `summary.coxph` fits, as recorded).
- **Bootstrap:** per-replicate **7.47 s → 5.49 s** (walls 149.4 → 109.9 s / 20 reps); **B = 1000 projection for the survival configuration: 124.5 → 91.5 min**.
- Realized vs predicted: medians-bucket recovery **2.61 s** on the profiled call (total-wall Δ 2.68 s) against the predicted ≈ 2.5 s — **realized in full, no shortfall**.
- Residual: the gbsg samples now sit in the **fit bucket, 3.69 s / 76.7 %** — `coxph()` wrapper/internals (`survival::coxph` total 2.52 s, `coxph.fit` self 0.13 s) plus model-frame/assembly machinery — exactly the material for the §9 unblocked assembly-skip follow-on.

---

## 6. Package health gates (§8)

- `devtools::test()` on 0.3.4 source: **0 failures | WARN 31 | 3 skipped | 4,864 passed** — warning-count parity with the recorded 0.3.3 baseline (31) and identical pass/skip counts; `tests/` untouched by this task, so `test-search-reproducibility.R` ran **unmodified** and passed. *GATE ✓*
- `devtools::check(document = FALSE, vignettes = FALSE, manual = FALSE)` (RSTUDIO_PANDOC set per the standing note): **0 errors | 0 warnings | 0 notes**, matching the recorded 0.3.3 result. *GATE ✓*
- Standing guards after an `R/` edit, both run on 0.3.4 source: `fidelity_fs_oc_predict_2026-08-28.R` — FIDELITY GATE: PASS (bit-identical, exit 0); `prerefactor_reference_2026-08-29.R check` — REFACTOR GUARD: PASS (identical to the 0.2.4 reference, exit 0).

---

## 7. Ten-line verdict

1. The `survfit` medians computation moved from inside `fit_cox_for_subgroup()` to the post-screen point of the same per-candidate iteration — a move, not a removal; call and extraction operation for operation.
2. Every equality gate passes: 37/37 `identical()` across nine fixtures × four components plus the 20-replicate gbsg bootstrap payload, under the reused volatile-field exclusion, self-consistency 3/3 at both stages.
3. `survfit` dispatches on gbsg: 1,410 → 121, exactly the passer set; adjusted arm 1,410 → 111; zero residual calls to attribute; continuous paths at zero both versions.
4. E-named reproduced the NA-median key string bit-identically; E-finite (`m1.threshold = 60`, excludes 5 / keeps 116) exercised the finite-filter consumer identically; E-zero no-opped with zero `survfit` calls; E-degen passed the screen as-built with warning parity.
5. Realized recovery: medians bucket 2.83 → 0.22 s (37.8 % → 4.6 %), gbsg call ≈ 7.5 → ≈ 4.8–5.4 s — the predicted ≈ 2.5 s realized in full.
6. Bootstrap replicate 7.47 → 5.49 s; B = 1000 survival projection 124.5 → 91.5 min.
7. Residual gbsg cost now sits in `coxph` internals + per-candidate assembly (fit bucket 76.7 %) — the §9 assembly-skip follow-on is the next lever, recorded, not implemented.
8. Tests 0 fail / WARN 31 / 4,864 pass (parity); `R CMD check` 0/0/0 (parity); `test-search-reproducibility.R` unmodified.
9. Touched surfaces: `R/subgroup_search.R`, `DESCRIPTION` (0.3.4), `NEWS.md` only; consumer sites read-only as gated; no push, no renders.
10. Deviations: none beyond the two recorded notes — the one-clause comment truthfulness fix inside `fit_cox_for_subgroup()`'s fast-branch comment, and the theoretical status-5/6 asymmetry for a screen-failing candidate whose `survfit` (but not `coxph`) errors, unreachable in every fixture.

## Commits

1. `ff2f9ca2` — task doc alone (`dev/tasks/cc_task_medians_on_survivors_2026-09-02.md`).
2. `4c923dfd` — Stage A: `medians_baseline_2026-09-02.{R,rds,log}`, `medians_compare_2026-09-02.R`, `medians_profile_2026-09-02.R`, `medians_profile_F2_before_2026-09-02.{Rprof,rds}`, `medians_profile_before_2026-09-02.log`.
3. `3921ffdd` — the implementation: `R/subgroup_search.R`, `DESCRIPTION` → 0.3.4, `NEWS.md`.
4. (this commit) — Stage C/verification: `medians_postchange_2026-09-02.{rds,log}`, `medians_compare_2026-09-02.log`, `medians_profile_F2_after_2026-09-02.{Rprof,rds}`, `medians_profile_after_2026-09-02.log`, `medians_check_034_2026-09-02.log`, `medians_test_034_2026-09-02.log`, this report.

No push. Tree clean at close.
