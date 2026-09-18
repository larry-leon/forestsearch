# REPORT — threshold sync: replicates resolve what the parent fit resolved

**Task:** `dev/tasks/TASK_threshold_sync_2026-09-18.md`, committed alone at `d0708011`.
**Date:** 2026-09-18. **Branch:** `feature/glm-extension`. **Machine:** pop-os. **Evidence base
re-verified from source:** the audit's §3.9 — pointers hold, line numbers had drifted. **Compute:** none;
the probe evaluates argument lists only. **Tests:** this task's two acceptance files, nothing else. No
`R CMD check`, no full suite. Not pushed. This starts the `dev/reports/` convention.

| SHA | What |
|---|---|
| `d0708011` | the task document, as received |
| `8bc8052b` | baseline replicate-equality table at `d0708011` |
| `111bd410` | Gate 1 pin: `quarto/qc/smoke_forestsearch_robustness.qmd` |
| `571b202c` | the sync (`forestsearch_main.R`, SECTION 2B-ii) |
| `a611928e` | acceptance tests + the probe moved beside them |
| `85847aa0` | the rider: binary default → `"OR"` at the estimation entry points |
| this one | NEWS + this report |

## 1. Source (Step 1)

`forestsearch_main.R`: `effect_measure` default `:1446-1452`; `user_set_*` detection + alias merge
`:1462-1465`; `args_call_all <- mget(names(formals()))` `:1512-1513`, *after* the merge; GLM resolution
`:1896-1990`; `threshold_config` `:2053-2067` (GLM) / `:2103-2117` (survival); `.sync_args_call_all()`
`:21-30`, called at `:1781, 2079, 2212, 2296, 2870`. Replays: `bootstrap_analysis_dofuture.R:406, 558, 614`
and `forestsearch_cross_validation.R:345, 388-421, 479` / `:859, 885-916, 1001`. The resolved thresholds
live in the locals `effect_threshold` / `consistency_threshold`, which are not formals and so were never in
`args_call_all`; `dmin.grf` escapes the same trap only because it *is* a formal, synced at `:2079`.

## 2. Gate 1 — exposure inventory: NOT empty

**Genuinely exposed, one file.** `quarto/qc/smoke_forestsearch_robustness.qmd`, `continuous` family
(MD). `fit_fs()` (`:164-170`) builds the call from `cfg$args` + `thr_*` (`:129-133`), which supply
`effect.threshold` but **no** consistency spelling; it then runs `forestsearch_tenfold()` (`:202`, `:239`)
and `forestsearch_Kfold()` (`:217`) on that fit. The fit resolved `c2 = 0.0`; every fold resolved `1.0`.

**Matches criterion (b) as written but provably unaffected, ten files** — explicit legacy supply sets
`user_set_*` TRUE in the parent too, so both sides resolve the same value:
`quarto/resampling/test_continuous_glm_compatibility.qmd` (MD, `hr.threshold = 0.5` / `hr.consistency = 0.0`,
+ tenfold `:325` and bootstrap `:539`); the four
`quarto/applications/actg175/_archive/20260520_analysis_actg175_continuous{,_hr,_hrMaxSG-both,_hrMaxSG-pareto}.qmd`
(MD, `10` / `5`, + bootstrap, K-fold, LOO); and `tests/testthat/helper-synthetic-dgm.R` `.fs_args_for()`
(`1.25` / `1.00` for **all** outcome types), consumed by `test-cross-outcome-parity.R`,
`test-no-subgroup-cv.R`, `test-cv-no-subgroup-edges.R`, `test-mr-inference.R`.

Verified clean: the MD campaigns and the `dev/glm-continuous-sims/*` drivers pass **both** new spellings —
their legacy-spelled argument lists are the survival fixtures; `dev/glm/*_test_suite.qmd` and the
graduation-step copy never bootstrap or CV their RD/MD fits; `_archive/…continuous_effMaxSG-both.qmd` uses
the new spellings; the OR campaign is ratio-scale throughout.

**Decision at the gate (Larry): pin first, then proceed.** `111bd410` adds `consistency.threshold = 0.0`
to that file's continuous `thr_found` / `thr_nosg` — the value its fit already resolved — so the fit is
unmoved, its folds now agree with it, and the sync changes nothing anywhere in committed work.

## 3. Baseline (Step 2) — the bug, measured

63 cells: estimand {survival, binary-unset(→OR), binary RD, binary OR, MD, IRR, IRD} × spelling
{none, legacy, new} × subset {neither, c1 only, both} × value {default, custom}, degenerate combinations
collapsed. `dev/reports/baseline_threshold_sync_2026-09-18.csv`, committed with the HEAD SHA.

**15 violations, all identity-scale**, so the "STOP if ratio or survival violate" branch did not fire:
binary/RD 5 of 9, continuous/MD 5 of 9, count/IRD 5 of 9; survival, binary-unset, binary/OR and count/IRR
zero. In every violating cell `c2` goes `0` → `1`; on the MD cell with nothing supplied `c1` also goes
`0` → `1.25`. The violating cells are exactly those leaving `c2` unsupplied — supplying **both**
thresholds, in either spelling, is already safe. Two cells (RD and IRD, nothing supplied) additionally
raised a spurious ratio-scale warning **once per replicate**. This reproduces the audit's §3.9 table
independently, cell for cell.

## 4. The sync (Step 3)

`forestsearch_main.R`, new **SECTION 2B-ii**, immediately after `threshold_config` is built: the resolved
values, on the natural scale, are written into `effect.threshold` / `consistency.threshold` and synced with
`.sync_args_call_all()`. Those spellings are `is.null()`-detected, so they survive any wrapper, and the
alias merge at `:1464-1465` makes the replay adopt them. The block runs **after** resolution and writes
only to `args_call_all`. The naturals are taken pre-`log()` (`hr.threshold` / `hr.consistency`, which the
ratio branch leaves alone) rather than as `exp(threshold_config$screening)`: a replay applies `log()`
again, and an `exp(log(x))` round-trip is not guaranteed to return `x` bit-for-bit. The resolution logic
is untouched — 36 inserted lines, **zero** deleted.

## 5. Gates (Step 5)

The probe was rewritten to be **source-driven**: it lifts the `effect_measure` default, the `user_set_*`
detection, the alias merge, the resolution block and the sync out of `body(forestsearch)` and rebuilds them
as a function whose formals are `forestsearch()`'s own — so `missing()` behaves exactly as in a real call
and the probe cannot drift from the source. Run with the sync excluded it reproduces the committed
transcribed baseline cell for cell.

| Gate | Result |
|---|---|
| **5a** replicate resolution == parent, every cell | **PASS** — 15 violations → 0 |
| **5b** parent resolution byte-identical to baseline | **PASS** — all 63 cells |
| **5c** ratio + survival replicate cells byte-identical | **PASS** — all 36 cells |

`dev/reports/postsync_threshold_sync_2026-09-18.csv`. The only change outside the 15 violating cells is
that six per-replicate ratio-scale warnings on RD / IRD stop firing; those cells already resolved equally.

## 6. The rider (Step 4)

`make_effect_estimator()` (`glm_effect_estimators.R`) and `.consistency_glm_pieces()`
(`consistency_resample.R`, the resolution `consistency_resample()` hands it) now default an unset binary
`effect_measure` to `"OR"`. Continuous, count and survival defaults are unmoved; an explicit measure is
untouched. `@param` blocks updated, `man/` regenerated (three `.Rd`; `NAMESPACE` unchanged).

**Exposure gate re-run**, not merely cited: every git-tracked `.R` / `.qmd` / `.Rmd` under `R/`, `tests/`,
`vignettes/`, `quarto/` parsed, each call matched the way R matches arguments — positional and partial
included. **61 call sites; 10 omit `effect_measure`; every one is survival or continuous.** No committed
caller reaches a binary default. This reproduces Directive B's §6
(`dev/directive_b/caller_inventory.txt`) and is stricter than a named-argument scan, which misreports the
positional call at `consistency_resample.R:426`.

## 7. Acceptance tests (Step 6)

`test-threshold-sync.R` (16 checks) and `test-binary-default-or-entry-points.R` (25). **41 pass, 0 fail,
0 skip. Wall clock 1.2 s** for the two files, 3.3 s including `load_all()`, against a 3-minute abort. No
fit, no resample, no fold. The probe now lives at `tests/testthat/helper-threshold-sync.R` so the tests and
the dev-side runs share one copy; `dev/tasks/probe_threshold_sync_2026-09-18.R` is a shim that sources it.
One test deliberately asserts the **pre-fix** configuration still produces exactly the 15 identity-scale
violations, so the probe cannot quietly stop measuring anything.

## 8. Findings (no tasks attached)

**S1 — criterion (b) is broader than the mechanism.** An explicitly supplied legacy threshold is safe: it
sets `user_set_*` TRUE on both sides. Only an *unsupplied* threshold on an identity-scale fit is exposed.
Ten of the eleven inventory hits were of the safe kind.

**S2 — `consistency_resample_compare()` has no binary default to flip.** Directive B's F2 names it with the
other two, but it has neither an `outcome_type` nor an `effect_measure` formal: it is survival-only and
calls `consistency_resample()` without `outcome_type`. B's inventory column recorded `outcome_type`, not
`effect_measure`, for its two call sites. A test records this.

**S3 — a second, currently unreachable `"RD"` default.** `make_effect_estimator()`'s binary branch calls
`match.arg(effect_measure, choices = c("RD", "OR", "RR", "IRR", "IRD"))`, and `match.arg()` returns the
first choice for a `NULL` argument — so the choice order is itself an `"RD"` default. Unreachable today
because the resolution above it always yields a length-1 value. Flagged, not changed.

**S4 — five further `binary = "RD"` resolution sites remain**: `frontier_cis.R:138`,
`forestsearch_cross_validation.R:1449` (and its `"RD"` fallback) and `:1498`, `plot_sg_glm_outcomes.R:147`.
Directive B's F1; outside this task's scope.

**S5 — the sync removes a per-replicate warning.** On RD / IRD with nothing supplied each replicate used to
emit "`effect.threshold` = 1.25 appears to be on a ratio scale". Resolution was already equal there; only
the warning stops. Worth knowing when comparing replicate logs across the change.

## Out of scope, untouched

Directive A (validation, derivation, `c2 <= c1` documentation). Directive C (DINA guard, frontier warnings,
display defaults). Floors. Anything in `fs-glms-interpretable`.
