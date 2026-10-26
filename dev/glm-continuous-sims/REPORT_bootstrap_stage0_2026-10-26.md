# REPORT — bootstrap efficiency Stage 0: measurement only, against CRAN 0.4.0

**Task:** `dev/tasks/cc_task_bootstrap_stage0_2026-10-26.md` (committed alone as `b428b94e`; received in `~/Downloads` as `cc_task_bootstrap_stage0_2026-10-26.md`, same name). **Bundle:** `dev/glm-continuous-sims/stage0_bootstrap_2026-10-26/` (from `~/Downloads/stage0_bootstrap_2026-10-26.tgz`, tarball route; all seven sha256 sums matched). **Run:** 2026-10-26, plain terminal on `pop-os` (not RStudio), `bash run_stage0.sh`, all five steps exit 0 (s0 05:16:12 → s4 05:18, ≈2.5 min wall).

## 1. Audit basis

- **Audited release:** CRAN `forestsearch` 0.4.0 (Date/Publication 2026-10-10).
- **Installed build:** `Version 0.4.0`, `Packaged: 2026-10-10 04:12:05 UTC; larryleon`, `Repository: CRAN`, `RemoteType: standard`, `RemoteSha: 0.4.0`, `Built: R 4.6.1; ; 2026-10-25 23:17:25 UTC; unix`, library `/home/larryleon/R/x86_64-pc-linux-gnu-library/4.6`. **The installed build is the CRAN build** (`Repository: CRAN`).
- **Environment:** survival 3.8.12 · future 1.80.0 · doFuture 1.4.0 · R 4.6.1 (2026-06-24) · Linux x86_64 · 128 available cores · `supportsMulticore() = TRUE` · `RSTUDIO` unset.
- **Repo:** branch `master`; HEAD at §1 = `4e77c37e830082c6258f0431695d6e179b2b6b2f`; s0 recorded `b428b94e…` because it ran after the task-document commit (that commit touches only `dev/tasks/`). `git status --porcelain` at §1, verbatim: `?? quarto/extreme_subgroups/fixed_random/` (pre-existing dirt; never staged).
- **Source of the 0.3.5 figures:** `dev/glm-continuous-sims/REPORT_bootstrap_closeout_2026-09-02.md` is the in-repo source. `HANDOFF_bootstrap_efficiency_2026-09-02.md` is a chat-side document, not a repo file. Closeout cites the single-call, replicate and candidate-count figures from "the assembly-skip record". That record is `REPORT_assembly_skip_2026-09-02.md`, where those numbers first appear, and §3 names it per row.

**Repo files under study vs CRAN 0.4.0 MD5** (`results/provenance_md5.csv`):

| file | CRAN 0.4.0 | repo | match |
|---|---|---|---|
| DESCRIPTION | 19b9459a03ed5c4a0c12943eb57d82b6 | 8c1cec4441bf8c2a8dff1ef720ec34b5 | DIFFERS |
| NAMESPACE | 00257d44d6f21d10f9af68b9f0dc227c | same | match |
| R/Cox_estimation_helpers.R | ffa9744bf4b3264041aa4c787cd5958f | same | match |
| R/bootstrap_analysis_dofuture.R | a891d66e7337b9a1ebdf524463102978 | same | match |
| R/bootstrap_calculations_helpers.R | 2c25b9b63263653c599f7143df6b854c | same | match |
| R/bootstrap_dofuture_main.R | 65cb60f5b23303578362172bc3cb9291 | same | match |
| R/bootstrap_dofuture_setup_helpers.R | 97cac83ad17bff5bf16efe7f437e164a | same | match |
| R/consistency_resample.R | 861c6b4e7d00c8176bbab6b7f863f724 | same | match |
| R/forestsearch_main.R | 11d80e74c4fcdba43fe8bbf9a1287d3f | same | match |
| R/get_fsdata.R | a7f71d659558706ae121d72e82099181 | same | match |
| R/guohe_algorithm3.R | 5a6e842f871c61ddef41dbfcc93bf6c5 | same | match |
| R/subgroup_consistency_helpers.R | f32b55b62d683337f6e1b36b3aa3bc17 | same | match |
| R/subgroup_consistency_main.R | 1285150ed4891883468d01d5ee981b78 | same | match |
| R/subgroup_search.R | f06d6bba8135f41a57962a659ebbe601 | same | match |

DESCRIPTION, git log since 2026-10-09: `99d0f987 chore: begin 0.4.0.9000 development`, `2405f11f docs: cite arXiv:2609.38361 via its arXiv DOI (CRAN incoming NOTE)`, `c373fd83 docs: cite Leon and Anderson (2026) in DESCRIPTION; NEWS entry`, `48674761 chore(release): DESCRIPTION 0.4.0, CRAN weightedsurv messages, ignore inst/patches`. The CRAN tarball's DESCRIPTION also carries CRAN-appended fields, so its MD5 can never match a repo copy. The survey makes no line reference to DESCRIPTION. **All twelve R files and NAMESPACE match CRAN 0.4.0 byte for byte, so every survey line reference still holds on `master`.**

## 2. Measurement record — `results/stage0_summary.md`, verbatim

# Stage 0 -- bootstrap efficiency measurements (CRAN forestsearch 0.4.0)

Installed: forestsearch 0.4.0 (2026-10-10 04:12:05 UTC; larryleon; repository CRAN) · survival 3.8.12 · future 1.80.0 · R version 4.6.1 (2026-06-24) · Linux/x86_64 · cores 128 · multicore TRUE
Repo: /home/larryleon/Documents/GitHub/forestsearch @ b428b94ec23fed6acc9275e97b1d86701ed539e8 (master); files matching CRAN 0.4.0: 13, differing: 1
Fixture: V: survival-post-selection vignette gbsg config, mr_inference = FALSE (n = 686)

## 1. The single sequential call (fixture V)

| quantity | value |
|---|---|
| call wall-clock, median of 3 (s) | 3.158 (reps 3.045, 3.197, 3.158) |
| candidates evaluated in the search | 2211 |
| candidates reaching the Cox fit | 1410 |
| candidates passing the effect screen | 121 |
| consistency stage: candidates total / evaluated / passed | 120 / 120 / 16 (early stop: FALSE) |
| identified subgroup | {er <= 0} & {pgr <= 26} |

Rprof, median of 3 profiles (seconds of total time; see rprof_watched_median.csv for the full list):

| function | total s |
|---|---|
| coxph (formula interface, all callers) | 2.04 |
| coxph.fit (+ agreg.fit) -- the fit proper | 0.33 |
| formula-interface machinery = coxph - coxph.fit | 1.71 (83.8% of coxph) |
| of which concordancefit | 0.45 |
| search_combinations_parallel (whole search stage) | 2.56 |
| fit_cox_for_subgroup (search-stage fits) | 1.98 |
| survfit (medians on effect-screen passers) | 0.11 |
| evaluate_consistency_twostage (whole consistency evaluation) | 0.51 |
| .consistency_cox_pieces / residuals.coxph | 0.26 / 0.10 |
| coxph time outside the search stage (Rprof coxph total - direct-timed search fits; candidate B's territory) | 0.15 |

## 2. Candidate A -- lean search-stage Cox fit (coxph.fit instead of formula coxph)

| quantity | value |
|---|---|
| captured search-stage fits (= candidates reaching the fit) | 1410 (matches filter counts: TRUE; adjusted: 0) |
| CI rows identical, lean vs formula (identical()) | 1410 of 1410 |
| CI values identical ignoring dimnames | 1410 |
| failure sets agree (try-errors) | TRUE (formula 0, lean 0) |
| warning sets agree | TRUE (candidates warned: 0; identical among them: NA) |
| formula path, all captured fits, median of 5 (s) | 1.887 (1.338 ms/fit) |
| lean path, all captured fits, median of 5 (s) | 0.378 (0.268 ms/fit) |
| **machinery removed per call (s)** | **1.509** (direct timers) vs 1.710 (Rprof coxph - coxph.fit, all callers) |
| machinery share of the formula path | 80.0% |
| sequential replicate now, mean of 5 (s) | 3.211 |
| **A: share of a sequential replicate** | **47.0%** (of the single call: 47.8%) |
| **A clears the 30.0% bar on the per-replicate (campaign) shape?** | **YES** |

Captured coxph.fit() conventions (the lean recipe, from the installed survival):

- coxph.fit formals: `x,y,strata,offset,init,control,weights,method,rownames,resid,nocenter`
- x: storage `double`, colnames `Treat`, values = treatment column: TRUE
- y: type `right`; equals Surv(Y,E): TRUE; equals aeqSurv(Surv(Y,E)): TRUE; control$timefix = TRUE
- strata NULL: TRUE · offset: numeric length 624, all zero: TRUE · init: NULL · control = coxph.control(): TRUE · weights NULL: TRUE · method: efron · rownames: character length 624 · resid: TRUE · nocenter: -1,0,1

## 3. Candidates C and D -- application shape (bootstrap at W workers)

W = 48 workers, nb = 192 for the M/C runs; identity pair at nb = 20.

| run | plan | nb | outer s | inner s | spawn s | overhead s | sum rep s | mean rep s | contiguous makespan s | balanced s | imbalance s |
|---|---|---|---|---|---|---|---|---|---|---|---|
| S | sequential | 5 | 18.6 | 18.4 | 0.2 | 2.4 | 16.1 | 3.21 | 16.1 | 16.1 | 0.0 |
| M | multisession | 192 | 30.7 | 28.9 | 1.8 | 4.3 | 964.4 | 5.02 | 24.6 | 20.1 | 4.5 |
| C | multicore | 192 | 28.2 | 28.1 | 0.2 | 2.4 | 987.9 | 5.15 | 25.6 | 20.6 | 5.1 |
| Mi | multisession | 20 | 11.3 | 10.3 | 1.0 | 3.7 | 85.3 | 4.27 | 6.6 | 4.3 | 2.3 |
| Ci | multicore | 20 | 9.6 | 9.5 | 0.1 | 2.5 | 88.8 | 4.44 | 6.9 | 4.4 | 2.5 |

| quantity | value |
|---|---|
| identity across plans, nb = 20 (all non-volatile components) | TRUE (results table only: TRUE) |
| per-replicate slowdown under load (M / sequential, C / sequential) | 1.56 / 1.60 |
| startup (spawn + load): multisession vs multicore (s) | 6.1 vs 2.6; removed: 3.6 |
| **C: measured wall reduction at nb = 192** | **8.1%** |
| **C: projected wall reduction at B = 1000** (startup + B·rep/W + imbalance) | **2.6%** (M 117 s -> C 114 s) |
| **D: static-chunk imbalance, share of the M run's wall** | **14.6%** (4.5 s of 30.7 s) |
| A inside the application shape at B = 1000: multisession 117 s -> 67 s; multicore 114 s -> 64 s | A+C combined: 45.3% |

## 4. Reading

- Per-replicate (campaign) shape: A removes 1.51 s of a 3.21 s replicate = 47.0%; bar 30.0%: cleared.
- Bit-identity of the lean fit: 1410 of 1410 CI rows identical; failure sets agree: TRUE.
- Application shape at B = 1000 (projection from measured components): C alone 2.6%; D (imbalance) 14.6% of the measured M wall; A+C 45.3%.
- Projections carry the per-replicate slowdown under load as measured and the imbalance as measured at nb = 4W; they are not a substitute for a timed B = 1000 run.


---

## §3 Against the 0.3.5 record

Record source: `REPORT_bootstrap_closeout_2026-09-02.md` (in-repo), with the figures it cites from `REPORT_assembly_skip_2026-09-02.md`. The "now" column quotes `results/keyvalues.txt` and `results/stage0_summary.md`.

| quantity | measured now (0.4.0, Stage 0) | record (0.3.5) | record figure: measured or extrapolated, outcome type | like-for-like? |
|---|---|---|---|---|
| single gbsg `forestsearch()` call | `call_s_median = 3.158` s (unprofiled, median of 3 after warm-up) | 3.58 s | **measured**, survival (F2 gbsg); an **Rprof-profiled** wall (`interval = 0.010`), `REPORT_assembly_skip_2026-09-02.md` §4 | partly. Same candidate flow, but the record is a profiled wall and this one is not, so profiler overhead is in the record figure only |
| sequential bootstrap replicate | `boot_S_mean_rep_s = 3.2112` s (mean of 5) | 4.18 s | **measured**, survival; 83.6 s / 20 sequential replicates (`REPORT_assembly_skip_2026-09-02.md` §4) | partly. Same configuration flow; nb 5 vs 20 |
| candidates reaching the fit / passing the screen | `search_n_reached_cox = 1410` / `search_n_passed_effect = 121` | 1,410 / 121 | **measured** (counted), survival F2 (`REPORT_assembly_skip_2026-09-02.md` §4 counts table) | **yes**. Fixture V reproduces the record's candidate flow exactly (s1 gate) |
| "coxph internals" per call | Rprof `coxph` total `rprof_coxph_total_s = 2.04` s, of which machinery `rprof_formula_machinery_s = 1.71` s; direct timers `machinery_s_per_call = 1.509` s | ≈2.4 s (survey's figure). Closeout states a 2.29 s fit bucket: wrapper/internals 1.26, model-frame 0.60, solve 0.26 | **measured** (Rprof bucket), survival F2 (`assembly_rebucket_035_2026-09-02.log`, cited in closeout item 5) | no. The record buckets by a re-bucketing scheme; Stage 0 splits `coxph` − `coxph.fit` by function. Different partition, different survival version |
| fixed startup under multisession at 48 workers | spawn + overhead, M run: 6.1 s (`boot_M_spawn_s = 1.832` + `boot_M_overhead_s = 4.317`) | ≈69 s | **inferred**, **continuous**: 75.2 s measured wall at nb = 100 minus ≈6 s ideal compute (closeout §4, line 56). It lumps every non-ideal second, including load slowdown and imbalance, into "startup" | **no**. Different outcome type and a different definition. The survival record has no startup figure |
| B = 1000 at 48 workers | projection 117 s (M), 114 s (C), from measured components (`stage0_summary.md` §3) | ≈3.9 min | **extrapolated**, survival: ×10 of a single measured nb = 100 wall of 23.3 s (closeout §4, line 59: "implied") | both are projections, not timed B = 1000 runs. The record scales a wall linearly, including its fixed part; Stage 0 adds startup + B·rep/W + imbalance |

**Not like-for-like across the board:**
- Survival version: 3.8.12 now; the records do not state a version.
- Machine: `pop-os` now; the records do not state a host.
- W = 48 in both, so workers match.
- The 0.3.5 → 0.4.0 code change between the two measurements.

## §4 Verdicts — numbers only

**A — lean search-stage Cox fit.**
- Machinery per call: `machinery_s_per_call = 1.509` s from direct timers over the 1,410 captured fits. Rprof corroborates it at 1.71 s, measured as `coxph` − `coxph.fit` across all callers.
- Share of the measured sequential replicate: **47.0%** of `boot_S_mean_rep_s = 3.2112` s (47.8% of the single call).
- 30% bar on the per-replicate shape: **cleared (YES)**.
- Identity: `ident_ci_row_identical = 1410` of `ident_n_candidates = 1410` (`identical()` on the consumed CI row, dimnames included). Failure sets agree: TRUE (0 / 0). Warning sets agree: TRUE (0 / 0).
- Captured `coxph.fit()` conventions (the recipe):
  - `x`: double, one column named `Treat`.
  - `y`: `Surv` type `right`. `control$timefix = TRUE`, so `aeqSurv()` is applied; on the probe candidate it equals the raw `Surv` (see OPEN ITEMS).
  - `strata`: NULL. `weights`: NULL. `init`: NULL. `control`: `coxph.control()`. `method`: `efron`.
  - `offset`: numeric zeros of length n. `rownames`: character of length n.
  - `resid`: TRUE. `nocenter`: `-1,0,1`.
  - No other formals.

**C — `multicore` vs `multisession` outer plan.**
- Startup removed: `startup_removed_by_multicore_s = 3.563` s (6.1 → 2.6 s).
- Measured wall reduction at nb = 192: **8.1%** (30.7 → 28.2 s).
- Projected reduction at B = 1000: **2.6%** (M 117 s → C 114 s; a projection).
- Identity across plans at nb = 20: `identity_multicore_vs_multisession_nb20 = TRUE` (results table only: TRUE).
- Per-replicate slowdown under 48-worker load vs sequential: `rep_slowdown_M = 1.56426`, `rep_slowdown_C = 1.60225`.

**D — static-chunk tail imbalance, judged from the Stage-0 M run only.**
- Imbalance: `imbalance_M_s = 4.47735` s (contiguous makespan 24.6 s vs balanced 20.1 s at nb = 192, W = 48, 4 replicates per chunk).
- Share of the M run's wall: **14.6%** (4.5 s of 30.7 s).

**B — not separately measured.**
- The Rprof line "coxph time outside the search stage" is **0.15** s: Rprof `coxph` total (2.04 s) minus the direct-timed formula path over the captured search fits (`time_formula_s_median = 1.887`).
- It subtracts a direct timer from an Rprof total, two different instruments, so it is not a precise figure. Rprof's own `.consistency_cox_pieces` total is 0.26 s, with `residuals.coxph` at 0.10 s.
- On fixture V the consistency stage ran without early stop (120 / 120 evaluated, `consistency_early_stop = FALSE`), so this is the full-evaluation shape.

No recommendation is made; whether to commission A, C or anything else is Larry's decision.

## §5 Script deviations

None. All five steps ran unmodified (exit 0); the seven sources are the bundle bytes as verified by sha256.

## §6 OPEN ITEMS

1. **Error and warning paths not exercised by A's identity test.** None of the 1,410 captured fits errored or warned: `ident_errors_formula = 0`, `ident_n_warned = 0`, `ident_identical_among_warned = NA`. The 1410/1410 identity therefore covers only the clean-fit path. The edge cases the survey lists are untested: zero-event arm under `d1.min = 0`, the monotone-likelihood warning, and n just above `n.min`.
2. **`aeqSurv` not discriminated.** On the probe candidate, `conv_y_values_equal_raw = TRUE` and `conv_y_values_equal_aeqSurv = TRUE`. That candidate has no near-ties that `aeqSurv()` would change, so the probe cannot show which form `coxph()` passes. The lean block applies `aeqSurv()` because `timefix = TRUE`, and 1,410/1,410 identity holds. Whether any fixture V candidate actually exercises the tie-correction path is not recorded.
3. **The record's ≈69 s startup is not reproduced for survival.** Stage 0 measures 6.1 s of spawn + overhead under multisession at 48 workers. The record's figure is an inferred residual from a continuous run, so the two are not comparable (§3). Whether continuous startup differs is outside this task (no continuous fixture).
4. **The B = 1000 figures on both sides are projections.** The record's ≈3.9 min is a ×10 extrapolation of an nb = 100 wall; Stage 0's 117 / 114 s come from summed components. No timed B = 1000 run exists for 0.4.0 (out of scope).
5. **No `current_status.md`.** `dev/glm-continuous-sims/current_status.md` does not exist; it was not created, and no pin commit was made.
6. **The 0.3.5 records state no survival version or host.** Machine and package-version comparability to the record cannot be established from the repo.

## §7 Ten-line verdict

1. Audit basis is CRAN `forestsearch` 0.4.0. The installed build is the CRAN build (`Repository: CRAN`), and all 12 studied R files plus NAMESPACE match the CRAN MD5s.
2. Fixture V reproduces the 0.3.5 candidate flow exactly: 1410 reach the fit, 121 pass the screen, 120 enter consistency, and the subgroup is `{er <= 0} & {pgr <= 26}`.
3. A single sequential call takes 3.158 s now (record: 3.58 s profiled). A sequential replicate takes 3.2112 s (record: 4.18 s).
4. A: the lean `coxph.fit()` path removes 1.509 s per call (80.0% of the formula path), which is 47.0% of a sequential replicate. The 30% bar is cleared.
5. A identity: 1410 of 1410 CI rows are `identical()`, with failure and warning sets agreeing. The error and warning paths were not exercised (OPEN ITEM 1).
6. C: multicore removes 3.563 s of startup. Wall drops 8.1% at nb = 192 and is projected to drop 2.6% at B = 1000. Payloads are identical across plans.
7. D: static-chunk imbalance is 4.47735 s, 14.6% of the M run's 30.7 s wall.
8. Under 48-worker load each replicate runs 1.56× (multisession) to 1.60× (multicore) slower than sequential.
9. The record's ≈69 s startup (continuous, inferred) and ≈3.9 min B = 1000 (survival, extrapolated) are not like-for-like with these measurements. Stage 0 projects 117 s at B = 1000.
10. No package code was modified (closing block). This report decides nothing; whether to commission A is Larry's call.

## Closing block — post-conditions

The checks were evaluated on the staged tree immediately before this report's commit. Check 1 was run over the committed range plus the staged index, which equals the post-commit range. All three were re-run after the commit with the same output. HEAD recorded in §1 = `4e77c37e`.

```
$ git diff --name-only 4e77c37e HEAD | grep -v '^dev/'
(nothing)
$ git status --porcelain | grep -v '^?? ' | grep -v '^.. dev/'
(nothing)
$ git status --porcelain | grep '^?? ' | grep -v ' dev/'
?? quarto/extreme_subgroups/fixed_random/      # pre-existing dirt, listed in §1
```

Commits (no push, no fetch): `b428b94e` (task document alone), then this one (bundle, `results/`, `logs/`, report). No `current_status.md` commit (OPEN ITEM 5).
