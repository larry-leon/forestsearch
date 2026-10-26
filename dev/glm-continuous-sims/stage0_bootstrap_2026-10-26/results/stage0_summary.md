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

