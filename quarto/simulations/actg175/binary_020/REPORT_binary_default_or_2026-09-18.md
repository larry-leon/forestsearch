# REPORT — Directive B: the binary default estimand becomes OR

**Date:** 2026-09-18 · **Task:** `dev/tasks/TASK_binary_default_or_2026-09-18.md` (`6d0ce4df`)
**Start SHA:** `02491f4f` (docs-task tree) · **Nothing pushed.**

## 1. Gates

| Gate | Result |
|---|---|
| **1** — site count / reachability | count claim **refuted**, reachability **confirmed**; dispositioned by Larry (option 1), scope unchanged |
| **5** — probe diff | **PASS** — 3 of 19 cells moved, all `binary`/`<unset>` |
| **6** — fit digests | **PASS** — 4 parity cells byte-identical, `binary_unset` moved to the OR path |
| **7** — task acceptance tests | **PASS** — 66 passing, 0 failures, 9.5 s |

## 2. Step 1 — the proof, and the refuted count

The task asserted two resolution sites in `forestsearch_main.R` and "no third anywhere in `R/`".
Searching all of `R/` **refutes the count**: eight more resolve a binary default — `"RD"` in
`make_effect_estimator()` (`glm_effect_estimators.R`), `.consistency_glm_pieces()`
(`consistency_resample.R`), two sites in `forestsearch_cross_validation.R`, `frontier_cis()`,
`plot_sg_glm_outcomes()` and two in `summary_utility_functions.R`; `"OR"` already in
`subgroup_glm()` (`run_subgroup_sims.R`), `generate_glm_dgm()` and `calibrate_glm_interaction()`.
The package held two opposed binary defaults across layers.

**Reachability confirms the audit's substance.** No extra site is reachable with `effect_measure`
still `NULL` from `forestsearch()`: both `make_effect_estimator()` calls and the MR `gspec` pass the
resolved value; the CV, summary and plot sites read it back off `args_call_all` / `fs$effect_measure`.

**Site 2 unreachable — proven.** `match.arg(outcome_type)` precedes both sites, so the `switch`
always returns non-`NULL` for the three non-survival types; the live site runs unconditionally for
`outcome_type != "survival"`; between the two sites `effect_measure` appears once, read-only, into
`.validate_outcome_threshold_config()`; the only `return()` in that span is the zero-row guard, which
exits; and no nested function definition lies between `forestsearch()`'s opening brace and site 2, so
`forestsearch()` is the only way in. Deleted on this proof in its own commit.

**`args_call_all`** is captured (`mget(names(formals()))`) **after** the live site, so bootstrap and
CV replays carry the resolved estimand and never re-resolve it.

## 3. Gate 5 — resolution probe, 19 cells

Baseline `6d0ce4df`, after `7e4e4c87`; artefacts in `dev/directive_b/`. The probe runs **no fits**:
it extracts the real top-level expressions from `body(forestsearch)` and evaluates them in a function
whose formals *are* `forestsearch()`'s, so `missing(hr.threshold)` behaves as in a real call and the
probe cannot drift from the source it measures.

Only the three `binary` / `<unset>` cells (one per `subgroup_method`) differ: resolved measure
`RD`→`OR`; `c1` `0.05`→log(1.25); screening natural `0.05`→`1.25`; consistency natural `0`→`1`;
scale `identity`→`log`; admission floor `0.05`→log(1.25); `args_call_all$effect_measure` `RD`→`OR`.
Every explicit cell (`OR`, `RD`, `RR`), survival, `MD` and `IRR` is **identical**. `c2` is numerically
`0` on both sides — log(1.0) = 0.0 — so the move shows in `consistency_natural` and `scale`.

## 4. Gate 6 — five seeded fits, bootstrap and MR off (4.5 s / 4.8 s)

| cell | baseline | after | verdict |
|---|---|---|---|
| `binary_OR` / `binary_RD` / `survival_default` / `continuous_MD` | — | — | all four byte-identical |
| `binary_unset` | `8a96d9cd…` (= `binary_RD`) | `d8099ef8…` (= `binary_OR`) | moved, as required |

The baseline `binary_unset` digest was byte-identical to `binary_RD`; the after digest is
byte-identical to `binary_OR`. `continuous_MD` selects nothing on both trees by construction
(`adverse_outcome` is `FALSE` for continuous and the factory's subgroup effect is beneficial) — still
a valid parity digest.

## 5. Acceptance checks

| # | Check | Result |
|---|---|---|
| 1 | binary unset → OR, log(1.25)/log(1.0), all three methods | **PASS** |
| 2 | binary unset ≡ binary `"OR"` explicit | **PASS** (probe, digest and test) |
| 3 | binary `"RD"` explicit unchanged (0.05/0.0, identity) | **PASS** |
| 4 | survival, MD, IRR byte-identical | **PASS** |
| 5 | exactly one resolution site remains | **PASS** — one `switch` in the deparsed body, yielding `"OR"` |
| 6 | campaign template re-resolves identically | **PASS** — passes `"OR"` explicitly; log(0.90)/log(0.80) |
| 7 | replayed `args_call_all` resolves the parent's estimand | **PASS** |
| 8 | identifier agreement across the three methods | **PASS** — below |
| 9 | `hr.threshold` roxygen fixed; rider roxygen-only | **PASS** — asserted mechanically |
| 10 | `NEWS.md` entry present | **PASS** — check findings: §7 |

**Check 8, verified from source, not repeated from the audit.** `.map_dina_family()` maps
`outcome_type = "binary"` to family `binomial`, and `forestsearch_helpers.R` builds the DINA harm
floor as `m_diff <- if (family == "gaussian") hr.threshold else log(hr.threshold)` — from the raw
`hr.threshold` formal, never from the resolved `effect_measure`. **DINA already screened on log-OR.**
Before the change, a default binary call had consistency and GRF screening `RD >= 0.05` while DINA
screened `log-OR >= log(1.25)`; they now agree. **DINA's own behaviour is unchanged.**

**Step 4 rider, verified before writing.** `fit_cox_for_subgroup()` returns the `"exp(coef)"` column
— the natural HR — and the survival path leaves `effect_threshold` `NULL`, so the natural
`hr.threshold` reaches the comparison: natural against natural. The old roxygen claimed log scale
"for ratio measures (OR, HR)", wrong for the measure the parameter is named after. The same wrong
claim appeared a second time in the same file on `evaluate_combination_with_status()`; corrected
identically.

## 6. Caller inventory (Larry's addition 2)

Read-only over git-tracked `R/`, `tests/`, `vignettes/`, `quarto/`; prose and `\link{}` excluded,
`dev/` excluded as archival. Full listing: `dev/directive_b/caller_inventory.txt`.

| entry point | call sites | leave `effect_measure` NULL |
|---|---|---|
| `make_effect_estimator()` | 48 (11 in `R/`) | 2 |
| `consistency_resample()` | 15 (6 in `R/`) | 6 |
| `consistency_resample_compare()` | 3 (0 in `R/`) | 3 |

**Not one committed caller can reach a binary default at any of the three.** Every site omitting
`effect_measure` is on a survival or continuous path — the `R/` ones are the `is.null(estimator_fn)` /
`cox_resample` branches, where `outcome_type` defaults to `"survival"`, while the GLM branches beside
them pass `glm_resample_spec$effect_measure`. The RD default there is reachable only by a direct user
call.

## 7. Check policy, amended mid-task

Larry amended the verification policy during this task, converging on: **per task, run the task's own
acceptance tests and nothing else; the full suite or any CRAN check runs only when Larry explicitly
asks.** Recorded in `CLAUDE.md` (`e500e5ff`, superseding `3427ca07` and `e68c5505`, kept in history).

**No `R CMD check` result exists for this task** — two `--as-cran` runs and one reduced run were
started and killed before completion, producing no artifact. The docs task's reference set (0 errors
| 1 warning | 2 notes) stands untouched; this task edited no file named in it.

**One-off data point:** a full `devtools::test()` did complete here, in **7.12 min** — 361 files,
**FAIL 0 | ERROR 0 | PASS 5160 | SKIP 3 | WARN 32**. It finished seconds before the instruction to
kill it, so nothing was killed. Under the final policy it is not run again.

**Gate as finally defined:** `tests/testthat/test-binary-default-or.R` — **66 passing, 0 failures,
9.5 s.** Confirmed failing pre-change: against the step-2 baseline tree (`49893dec`) in a detached
worktree it gives **FAIL 30 | PASS 36**. Checks 1, 2, 5, 7, 9 move; 3, 4, 6 pass on both trees and
exist to catch collateral movement.

## 8. Commits

| SHA | Class | Subject |
|---|---|---|
| `6d0ce4df` | record | task document into `dev/tasks/`, alone |
| `49893dec` | record | step 2 baseline, before any edit |
| `5b1dac43` | **changes behaviour** | binary default `"RD"` → `"OR"` at the live site |
| `1ba22ded` | **removes code** | delete the unreachable second site, on the proof |
| `fda9dcb0` | doc only | `effect_measure` param states OR as the binary default |
| `7e4e4c87` | doc only | `hr.threshold` scale rider |
| `f391b816` | tests + docs | acceptance tests, `NEWS.md`, regenerated `man/` |
| `3427ca07`→`e68c5505`→`e500e5ff` | doc only | `CLAUDE.md` check policy, amended to final text |

Both `R/` doc-only commits were asserted mechanically: every changed line is a roxygen line.
`document()` regenerated `forestsearch.Rd`, `subgroup.search.Rd`,
`evaluate_combination_with_status.Rd`. The threshold-vocabulary table needed no change — its rows are
keyed by *resolved estimand*, not outcome type, so its RD row already describes explicit RD.

## 9. Findings, no task attached

**F1 — the eight extra resolution sites** (§2). Two opposed binary defaults across layers; not
reachable with `NULL` from `forestsearch()`, so not a live defect today.

**F2 — exported entry points now disagree with `forestsearch()`.** `make_effect_estimator()`,
`consistency_resample()` and `consistency_resample_compare()` still default binary to `"RD"` while
`forestsearch()` defaults to `"OR"`. No committed caller reaches it (§6); a direct user call does.
**A follow-up directive will disposition this separately.**

**F3 — an explicit RD threshold becomes a log-ratio floor under `subgroup_method = "dina"`.**
`m_diff` is `log(hr.threshold)` whenever the DINA family is not `gaussian`, and family is mapped from
`outcome_type` alone, so binary + `effect_measure = "RD"` + `effect.threshold = 0.07` gives a floor of
`log(0.07) = -2.659` against DINA's log-OR — admitting nearly everything. Pre-existing.

**F4 — vignette build fails on an untouched tree.** Verbatim, no diagnosis:

```
* creating vignettes ... ERROR
--- re-building ‘forestsearch.Rmd’ using rmarkdown
Error: processing vignette 'forestsearch.Rmd' failed with diagnostics:
Pandoc is required to build R Markdown vignettes but not available. Please make sure it is installed.
--- failed re-building ‘forestsearch.Rmd’
SUMMARY: processing the following file failed:
  ‘forestsearch.Rmd’
Error: Vignette re-building failed.
```

Reproduced identically on the untouched pre-change tree at `02491f4f`. It restates the environment
note already in `REPORT_threshold_docs_2026-09-18.md`, which this session rediscovered from scratch —
the reason the `RSTUDIO_PANDOC` export now lives in `CLAUDE.md`.

**F5 — pre-existing roxygen defect in `R/fs_bias_coverage.R`.** `document()` reports `@description
failed to evaluate inline markdown code` / `Failed to parse the inline R code: \`r = se_mean /
sd_emp\``. Same file as the docs task's non-ASCII warning. Not touched here.

**F6 — the truth layer is already unconditionally OR for binary.** `.fs_region_effect()`
(`R/betaHhat_truth.R`) ignores `effect_measure` on the binary branch and always returns
`exp(logit coefficient)`.

## 10. Compute

Two fit batches of five (4.5 s, 4.8 s), two probe runs (no fits), one replay inside the tests, tests
9.5 s. Well inside the two-minute target; the ten-minute abort was never approached. Bootstrap and CV
never ran.
