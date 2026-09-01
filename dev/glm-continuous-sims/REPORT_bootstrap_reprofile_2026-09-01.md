# REPORT — re-profiling the search at 0.3.2: where the bootstrap's cost sits now that the dispatch overhead is gone

**Task:** `dev/tasks/cc_task_bootstrap_reprofile_2026-09-01.md` (commit `d7d5aae2`)
**Repo / branch:** `~/Documents/GitHub/forestsearch` (pop-os), `feature/glm-extension`
**Category:** measurement only. No `R/` change, no optimisation, no edit to any estimator, driver or document. Artefacts, all under `dev/glm-continuous-sims/` (chosen to sit beside the 08-30 baseline this task supersedes): `bootstrap_reprofile_2026-09-01.R`, `bootstrap_reprofile_2026-09-01.log`, `bootstrap_reprofile_2026-09-01.rds`, three raw `bootstrap_reprofile_<config>_2026-09-01.Rprof` sample files, this report. No renders, no push, no install.

---

## 1. Provenance (§1 — GATE passed)

```
pop-os · /home/larryleon/Documents/GitHub/forestsearch · feature/glm-extension
HEAD at task start 2a9ceeb6 · git status --porcelain at task start: empty
0.3.1 dispatch commits in the log: fbd564de (the fix: sequential plan -> plain loop,
  31.1 -> 5.4 s per call, bootstrap 32.2 -> 4.95 s per replicate), with 911d73ba
  (the narrowed task doc) and 774d768c (the reverted first attempt) beside it
0.3.2 orientation commit in the log: eb136a35 (signed orientation in
  fs_oc_family_enumerate(); report 15ff5640)
packageVersion("forestsearch") 0.3.2 = DESCRIPTION Version: 0.3.2 · R 4.6.1 · 128 cores
```

`~/Downloads/cc_task_bootstrap_reprofile_2026-09-01.md` was absent under that name; the stem arrives hyphen-stripped as **`cc_task_bootstrap_reprofile_20260901.md`** (as with the 08-30 task) → `dev/tasks/cc_task_bootstrap_reprofile_2026-09-01.md`, committed alone as **`d7d5aae2`**. The measurement script ran at HEAD `d7d5aae2`.

**`vi.grf.min` default in force: `NULL`** (`R/forestsearch_main.R:1249`; confirmed `formals(forestsearch)$vi.grf.min` at run time). Every call below uses `vi.grf.min = NULL` — no fixture passes it, so the installed default is what runs.

**Configuration note.** The task says "the current default `parallel_args` (i.e. the plain loop)". The *formal* default of `forestsearch()` is `list(plan = "multisession", workers = .default_parallel_workers(), show_message = TRUE)` — not a plain loop — so the parenthetical was read as governing: the continuous and survival fixtures are profiled at the **replicate configuration** `list(plan = "sequential", workers = 1L)`, which at 0.3.2 takes the plain loop (the 0.3.1 fix short-circuits on the plan alone, `R/subgroup_consistency_main.R` Section 6, L693–700) and is both the configuration every bootstrap replicate runs (`bootstrap_analysis_dofuture.R`) and the configuration the 08-30 baseline profiled, so then/now is like-for-like. The anchor fixture runs its document's own `parallel_args = list(plan = "sequential")` verbatim — at 0.3.2 the two spellings take the identical path.

## 2. Fixtures (§2)

The first two are **the 08-30 profile's, reproduced from its script verbatim** (fixture blocks copied unchanged from `bootstrap_profile_2026-08-30.R`). Quoting that report's §3.1:

- **Continuous / MD.** "The MD40 DGM the OC-wrapper scripts rebuild (`generate_glm_dgm` on ACTG175 arms 1/3, `k_inter = truth$beta_inter` from the tracked n = 500 payload, gated on `effect_Q` and `prevalence_Q`), one trial `simulate_from_glm_dgm(dgm, n = 500, seed = 8316951 + 1)` (sim_id 1), and the driver's `forestsearch()` arguments verbatim (`sim_fs_maxeffCons_mr_md40_knoise0_n500_batch_1_1000.qmd` L561–590): 13 confounders incl. `str2`, `effect.threshold = 30`, `consistency.threshold = 10`, `pconsistency = 0.90`, `fs.splits = 400`, `n.min = 60`, `d0.min = d1.min = 12`, `maxk = 2`, `sg_focus = "maxeffCons"`, `consistency_method = "resample"`, `conf.cont_jcuts = list(age = 10, preanti = 10)`, `use_twostage = TRUE`, `is.RCT = TRUE`, `seedit = 8316952`." (Also `stop_threshold = NULL` explicitly — this matters for the anchor contrast below.)
- **Survival / gbsg.** "`survival::gbsg` (686 rows, 299 events), the effMaxSG application's arguments (`quarto/applications/gbsg/analysis_gbsg_survival_effMaxSG.qmd` L539–585 with its focus profile resolved): confounders `age, meno, size, grade3, nodes, pgr, er`, `hormon` treatment, `time_months = rfstime / 30.4375`, `est.scale = "hr"`, `conf.cont_jcuts = list(er = 10, pgr = 10)`, `conf_force = c("er <= 0", "pgr <= 0")`, `sg_focus = "effMaxSG"`, `n.min = 60`, `d0.min = d1.min = 10`, `maxk = 2`, `hr.threshold = hr.consistency = 1.0`, `pconsistency = 0.90`, `minp = 0.025`, `use_twostage = TRUE`, `consistency_method = "resample"`, `fs.splits = 1000`, `max.minutes = 3`, `vi.grf.min = NULL`."
- **Anchor (scale only, new).** The ACTG175 applied anchor: the fixed-family `maxeffCons` call of `quarto/applications/actg175/analysis_actg175_continuous_oc.qmd` §2 (`{r anchor}` chunk, L89–120) with its data-prep (L54–64: arms 1/3, `y_decline = cd40 − cd420`, N = 1083, six continuous + six binary confounders) and `params$seed = 8316951`, verbatim: `adverse_outcome = TRUE`, `effect.threshold = 10`, `consistency.threshold = 5`, `conf.cont_jcuts` = 10 on all six continuous confounders, `fs.splits = 500`, `parallel_args = list(plan = "sequential")`. One profiled call; no bootstrap on it. Its selected subgroup reproduced the document's gate: `{age <= 37} & !{cd40 <= 507}`.

**One inherited difference surfaced by the anchor.** The anchor omits `stop_threshold`, whose formal default is `pconsistency.threshold` (`R/forestsearch_main.R:1235`) — and first-passer stopping is sound for `maxeffCons`, so the anchor's consistency screen evaluated **1** candidate and stopped. The continuous driver fixture passes `stop_threshold = NULL` explicitly, so all 749 are evaluated. The anchor therefore shows the fit bucket almost undiluted.

## 3. Profile — the same buckets, re-shared (§3)

Procedure per fixture as in the 08-30 script: an instrumented call (`trace()` counters, per-fit timing), one un-timed warm-up after `untrace()`, then the profiled call under `Rprof(interval = 0.010, line.profiling = TRUE)` (10 ms is this platform's minimum; Rprof warns below it). **Bucket definitions are the 08-30 script's `BUCKETS` list verbatim**, priority order medians › GRF › fit › cuts › consistency evaluator › future machinery › dedup/selection › enumeration › rest of `subgroup.consistency()` › everything else.

### 3.1 Wall-clock

| config | `parallel_args` | instrumented wall | **profiled wall** | samples | sampled / wall |
|---|---|---:|---:|---:|---:|
| cont (replicate configuration) | `plan = "sequential", workers = 1L` | 5.21 s | **4.74 s** | 473 | 99.8 % |
| surv gbsg | `plan = "sequential", workers = 1L` | 7.72 s | **7.59 s** | 757 | 99.8 % |
| anchor (ACTG175 applied) | `plan = "sequential"` | 8.21 s | **8.19 s** | 817 | 99.8 % |

Against the narrow-dispatch report's post-fix measurements: continuous ≈ 5.4 s per call (its C1 variant), gbsg 7.41 s — this run's 4.74 s and 7.59 s sit beside them.

### 3.2 Bucket tables — share of the current total beside the 08-30 share, absolute seconds then and now

**Continuous, replicate configuration** (0.3.0 wall 29.42 s → 0.3.2 wall 4.74 s):

| bucket | secs 08-30 | secs now | share 08-30 | **share now** |
|---|---:|---:|---:|---:|
| per-candidate effect fit | 2.60 | 2.93 | 8.8 % | **61.9 %** |
| candidate medians | 0 | 0 | 0 | 0 |
| GRF variable importance | 0 | 0 | 0 | 0 |
| cut construction | 0 | 0 | 0 | 0 |
| enumeration / floors | 0.01 | 0.01 | 0.03 % | 0.2 % |
| consistency screen (evaluator) | 1.56 | 1.40 | 5.3 % | **29.6 %** |
| dedup / selection | 0.01 | 0.01 | 0.03 % | 0.2 % |
| future machinery | 25.17 | 0.31 | 85.6 % | 6.6 % |
| everything else | 0.05 | 0.07 | 0.2 % | 1.5 % |

**Survival gbsg** (8.19 s → 7.59 s):

| bucket | secs 08-30 | secs now | share 08-30 | **share now** |
|---|---:|---:|---:|---:|
| per-candidate effect fit (excl. survfit) | 3.53 | 3.97 | 43.2 % | **52.4 %** |
| candidate medians (`survfit` + `summary.survfit`) | 2.54 | 2.69 | 31.1 % | **35.5 %** |
| consistency screen (evaluator) | 0.60 | 0.61 | 7.3 % | 8.1 % |
| future machinery | 1.45 | 0.26 | 17.7 % | 3.4 % |
| enumeration / floors | 0.01 | 0.00 | 0.1 % | 0 |
| dedup / selection | 0.00 | 0.01 | 0 | 0.1 % |
| GRF / cuts | 0 | 0 | 0 | 0 |
| everything else | 0.04 | 0.03 | 0.5 % | 0.4 % |

**Anchor, ACTG175 applied** (8.19 s; no 08-30 baseline — scale reference only):

| bucket | secs | share |
|---|---:|---:|
| per-candidate effect fit | 7.23 | **88.5 %** |
| future machinery (the search's own in-process `future_lapply`) | 0.85 | 10.4 % |
| consistency screen (evaluator) | 0.00 | 0 (1 candidate evaluated — default `stop_threshold`, §2) |
| medians / GRF / cuts / enumeration / dedup | 0.00 each | 0 |
| everything else | 0.09 | 1.1 % |

**Absolute seconds where the fix did not touch the work.** The fit buckets read +0.33 s (cont, 2.60 → 2.93) and +0.44 s (gbsg, 3.53 → 3.97), medians +0.15 s, consistency evaluator −0.16 s / +0.01 s. The `trace()`-timed fit totals are, by contrast, essentially identical run-to-run: cont **2.94 s (08-30) → 2.96 s (now)**, gbsg 6.48 → 6.76 s (that timer wraps the whole of `fit_cox_for_subgroup`, i.e. fit + medians buckets, 3.97 + 2.69 = 6.66 s of samples). So the real work is unchanged, as expected — the ±0.3–0.4 s bucket deltas are 10 ms-sampling attribution and run-to-run spread, not new cost. The one bucket the fix touched, future machinery, fell 25.17 → 0.31 s (cont) and 1.45 → 0.26 s (gbsg); the ~0.3–0.9 s that remains is the *search phase's* single in-process `future_lapply` (`search_combinations_parallel`), not the consistency batching, which is gone.

*GATE (§3):* the bucket parser assigns every sample, so the bucket sum equals the sampled total by construction; the sampled total is **99.8 % of the measured wall on all three fixtures**, and the unattributed "everything else" bucket is ≤ 1.5 %. **Passed.**

### 3.3 Counts and per-fit times (`trace()`; per-call `proc.time()` deltas)

| fixture | candidates reaching a fit | fits NULL | mean / median per fit | passed effect screen | dedup before → after | reaching the consistency screen | passing it |
|---|---:|---:|---:|---:|---:|---:|---:|
| cont | 1 975 | 0 | 1.5 ms / 1.0 ms | 852 | 852 → 749 | 749 | 319 |
| surv gbsg | 1 410 | 0 | 4.8 ms / 5.0 ms | 121 | 121 → 120 | 120 | 9 |
| anchor | **4 702** | 0 | 1.5 ms / 1.0 ms | 146 | 146 → 131 | **1** (early stop) | 1 |

Counts are identical to 08-30 where comparable (cont: 1 975 / 852 / 749 / 319; gbsg: 1 410 / 121 / 120 / 9), and every fixture selected the same subgroup on the instrumented and profiled calls. The closed-form consistency never returned `NA` and the literal-split fallback fired 0 times (as at 08-30). The anchor's candidate space is 2.4× the MD40 fixture's (4 702 vs 1 975 fits at the same 1.5 ms mean), which is what pushes its fit bucket to 88.5 %: **as M grows, the call converges to pure fit cost** — the buckets move exactly the way the task supposed.

## 4. The two source questions (§4)

### 4(a) The per-candidate fit: machinery versus arithmetic

**Estimator path, unadjusted two-group case, from current source.**

*Continuous:* `subgroup_search.R:618` — `glm_result <- fit_glm_for_subgroup(df_clean, id.x, estimator_fn)`. `fit_glm_for_subgroup()` (`R/subgroup_search.R:798`) does `df_sg <- df_clean[id.x == 1, , drop = FALSE]` and calls the closure. The closure is built once per `forestsearch()` call (`R/forestsearch_main.R:1721–1729`, `make_effect_estimator(... adjust_covariates = adjust_covariates)`) and for `outcome_type = "continuous"` is `.make_lm_estimator()`'s function (`R/glm_effect_estimators.R:803`), which **per candidate** rebuilds the formula, refits, and inverts:

```r
fmla <- .build_adjusted_formula(outcome.name, treat.name,                  # L826
          if (ps_adjust_method == "none") adjust_covariates else NULL)
...
      fit <- stats::lm(fmla, data = data_slice)                            # L833
      coef_val <- stats::coef(fit)[[treat.name]]                          # L835
      se_val   <- sqrt(diag(stats::vcov(fit)))[[treat.name]]              # L836
```

(`vcov.lm` runs `summary.lm`.) *Survival:* `subgroup_search.R:657–659` — `cox_result <- fit_cox_for_subgroup(yy, dd, tt, id.x, df_clean = df_clean, adjust_covariates = adjust_covariates)`; the body (`R/subgroup_search.R:722–748`) assembles a `data.table`, keeps `rhs <- "Treat"` unless `.fs_adjust_terms(adjust_covariates)` is non-empty (L725, L731–739), then per candidate:

```r
cox_fmla <- stats::as.formula(paste0("survival::Surv(Y, E) ~ ", rhs))      # L744
hr.cox <- try(
  summary(survival::coxph(cox_fmla, data = df.x, robust = FALSE))$conf.int, # L745-746
  silent = TRUE)
```

**The measured split** (Rprof sub-classification of the fit bucket's samples, priority solve › model.frame/model.matrix › vcov/summary › formula construction › subset/assembly › wrapper/other; solve = `lm.fit` / `glm.fit` / `coxph.fit` / `agreg.fit`):

| part of the fit bucket | cont (2.93 s) | surv (3.97 s) | anchor (7.23 s) |
|---|---:|---:|---:|
| **the actual solve** | **4.8 %** (0.14 s) | **6.8 %** (0.27 s) | **2.1 %** (0.15 s) |
| subset / data assembly (`[.data.frame`, `data.table`, …) | 58.7 % (1.72 s) | 25.9 % (1.03 s) | 63.9 % (4.62 s) |
| `model.frame` / `model.matrix` | 19.5 % (0.57 s) | 16.9 % (0.67 s) | 17.8 % (1.29 s) |
| `vcov` / `summary.*` | 8.2 % (0.24 s) | 40.3 % (1.60 s) | 8.6 % (0.62 s) |
| formula construction | 3.1 % (0.09 s) | 0.5 % (0.02 s) | 1.0 % (0.07 s) |
| wrapper / other | 5.8 % (0.17 s) | 9.6 % (0.38 s) | 6.6 % (0.48 s) |

**The arithmetic is 2–7 % of the fit bucket; 93–98 % is machinery** — dominated on the continuous path by data-frame subsetting/assembly plus `model.frame`/`model.matrix` (the top self-time functions of the whole call are `NextMethod`, `[[.data.frame`, `sys.call`, `[.factor`, `[.data.frame`; `lm.fit` is 0.11 s self of 4.73 s), and on the survival path by `summary.coxph` (1.60 s — computing the full summary table to read one confidence interval) plus assembly.

**Is an unadjusted fast path routable by condition?** Yes, on both paths, by a predicate available at the call site with no runtime cost:

- *Continuous:* every branch inside the closure is decided by values **fixed at closure construction** (`force()`d): the weighted branch requires `ps_adjust_method == "iptw"` (L830), the adjusted formula requires `adjust_covariates` non-empty under `ps_adjust_method == "none"` (L826–827, `.build_adjusted_formula` L870–878), and the lm path takes no offset (offsets exist only on the Poisson rate path, enforced at `make_effect_estimator()` L121–128). So the single predicate `ps_adjust_method == "none" && length(.fs_adjust_terms(adjust_covariates)) == 0L`, evaluable once at `R/forestsearch_main.R:1721` where both values are in scope, characterises exactly the calls that fit `y ~ treat`, unweighted, no offset — every call in all three fixtures. The adjusted/weighted paths would be untouched.
- *Survival:* `fit_cox_for_subgroup()` accepts no weights, no offset, no IPTW at all; its only branch is `length(.fs_adjust_terms(adjust_covariates)) > 0L` (L732). The predicate `is.null(adjust_covariates)` (equivalently `length(.fs_adjust_terms(adjust_covariates)) == 0L`) at the call site `subgroup_search.R:657–659` selects the treatment-only model `Surv(Y, E) ~ Treat` exactly. (`forestsearch()` additionally guards that GRF-adjustment and `adjust_covariates` interplay at L2057, all resolved before the search starts.)

No change is designed here; the predicate finding is the deliverable.

### 4(b) The survival medians: how many candidates actually need them?

**The two consumer sites, quoted from current source.** (i) The survival dedup key — `R/subgroup_consistency_helpers.R:1108–1121`, called at `R/subgroup_consistency_main.R:594`:

```r
# Columns to check: K, n, E, d1, m1, m0, HR, L(HR), U(HR)
cols_to_check <- 2:min(10, ncol(df))                                       # L1109
...
key_cols <- df_rounded[, cols_to_check, drop = FALSE]                      # L1118
dup_key <- apply(key_cols, 1, function(x) paste(x, collapse = "_"))        # L1119
keep_rows <- !duplicated(dup_key)                                          # L1121
```

`m1`/`m0` are columns 6–7 of the nine key columns. (ii) The `m1.threshold` filter — `R/subgroup_consistency_main.R:531–546`, behind a formal whose default is `Inf` (`R/forestsearch_main.R:1228`; `subgroup.consistency()` L357):

```r
if (is.finite(m1.threshold)) {                                             # L531
  hr.subgroups <- hr.subgroups[!is.na(hr.subgroups$m1), ]                  # L532
  ...
  found.hrs <- hr.subgroups[hr.subgroups$HR >= hr.threshold &              # L545
                              hr.subgroups$m1 <= m1.threshold, ]           # L546
```

No non-archived driver or application sets it finite: the two live call sites that name it pass `Inf` explicitly (`quarto/resampling/gbsg_analysis_cox_effMinSG_split_vs_resample.qmd:128`; `quarto/applications/gbsg/analysis_gbsg_survival_multimethod.qmd:759, 1325`), unchanged from the 08-30 survey.

**The measured ratio.** Medians are computed inside `fit_cox_for_subgroup()` (L761–763) for **every candidate whose Cox fit succeeds** — on gbsg, all **1 410** fitted candidates — while both consumers act on `hr.subgroups`, the effect-screen survivors: the `m1.threshold` filter reads that frame (L531–546), and the dedup receives it at **121 rows** (→ 120). **The ratio is 1 410 / 121 = 11.7×**: medians are computed for 11.7 candidates for every one that reaches a consumer. The medians bucket is 2.69 s (35.5 % of the gbsg call); computed on the 121-row set instead it would be ≈ 0.23 s — **the prize is ≈ 2.5 s per gbsg call, ≈ 32 % of it**. Both consumers *could* be served from the later, smaller set: they read `m1`/`m0` only from `hr.subgroups` rows, all present at the point the effect screen has finished, and the remaining reads in the package are display tables for the selected subgroup only (08-30 §2.1 survey; no site reads medians for candidates that failed the effect screen). The anchor shows the same shape at larger M for the general row-builder pattern (4 702 fits → 146 screen survivors, ratio 32.2 — though continuous fits set `m1 = m0 = NA` and pay nothing). No change designed; the ratio and the served-from-later finding are the deliverable.

## 5. The extrapolation, confirmed (§5)

`forestsearch_bootstrap_dofuture()` at **`nb_boots = 20`** (the low end of the task's 20–50, matching the 08-30 and dispatch-task checks), continuous fixture, `plan = "sequential", workers = 1L`, seed 8316952:

| quantity | value |
|---|---:|
| bootstrap wall (20 replicates, sequential) | 100.4 s |
| **per replicate** | **5.02 s** |
| the bootstrap's own `tmins_search × 60`: mean / median | 4.86 s / 4.95 s |
| the bootstrap's own `tmins_iteration × 60`: mean | 4.87 s |
| replicates with an identified subgroup | 20 of 20 |
| profiled single call | 4.74 s |
| ratio per-replicate / single call | **1.059** |
| ratio vs the narrow-dispatch task's 4.95 s per replicate | **1.015** |
| **projection, B = 1000, sequential** | **5 022 s = 83.7 min** (narrow-dispatch projected ≈ 82 min) |

Well inside the ~20 % band: nothing outside the search dominates — the replicate's iteration time (4.87 s) is its search time (4.86 s), and the per-replicate mean agrees with the dispatch task's 4.95 s to 1.5 %. The profile's ranking transfers to the bootstrap unchanged.

## 6. The batch-size question, settled by measurement (§6)

Two direct `forestsearch()` calls on the continuous fixture under `parallel_args = list(plan = "multisession", workers = 6L)` — the configuration class the applications run, no longer a bootstrap item (replicates take the plain loop). 749 candidates reach the screen; both calls selected the same subgroup (`!{wtkg <= 84} & !{cd40 <= 320}`).

| batch size | `future_lapply` rounds | wall |
|---|---:|---:|
| default `min(n_workers * 2L, n_candidates)` = 12 (`R/subgroup_consistency_main.R:893`) | 63 | **20.21 s** |
| `parallel_args$batch_size = 125` (one batch per worker, ⌈749/6⌉) | 6 | **12.80 s** |

The 57 extra rounds cost 7.4 s — ≈ 0.13 s of globals resolution and future dispatch per round — making the default-batch call 1.6× the one-batch-per-worker call on identical work. (Both walls include the multisession pool start-up and the search phase, common to the two.) **Measurement only; no change proposed.**

## 7. Ranked by measured current share (§7)

What is worth changing, ranked by the share of a **bootstrap replicate** (continuous, 5.02 s; gbsg call in brackets) it would recover. **Nothing is changed here.**

| rank | what | measured current share | recoverable | risk note |
|---|---|---:|---|---|
| **1** | **The per-candidate fit's machinery** — 95.2 % of the fit bucket is subsetting/assembly, `model.frame`/`model.matrix`, `vcov`→`summary`, formula rebuild, against 4.8 % `lm.fit` (§4a); an unadjusted fast path is routable by a single construction-time predicate on both paths | **61.9 %** (2.93 s) [gbsg fit 52.4 %, 3.97 s, solve 6.8 %] | up to ≈ 2.6–2.8 s per continuous replicate if the fast path approaches the solve's cost — B = 1000 from ≈ 84 min toward ≈ 40 min; on gbsg up to ≈ 3.5 s per call (with rank 2, ≈ 6 s of its 7.6) | Medium: the closure's `estimate`/`se` must reproduce to the digits the payloads store; the closure is shared with MR, the bootstrap's H/Hc fits and CV, so the predicate must gate at construction, not by changing the shared path |
| **2** | **Survival medians computed for all fitted candidates** — needed by consumers only on the effect-screen survivors: 1 410 computed vs 121 consumed, ratio **11.7×** (§4b) | [gbsg **35.5 %**, 2.69 s] — 0 on continuous | ≈ 2.5 s per gbsg call (compute on the 121-row set: ≈ 0.23 s) | Low-medium: the dedup key keeps its `m1`/`m0` columns if medians are attached to `hr.subgroups` before dedup — a move, not a removal; `m1.threshold` (default `Inf`, nothing sets it finite) reads the same frame; display sites need the selected subgroup only |
| 3 | The consistency evaluator — 749 closed-form calls at ≈ 1.9 ms | 29.6 % (1.40 s) [gbsg 8.1 %, 0.61 s] | little per call (already the cheap closed form); the lever is the **count** — the anchor's default `stop_threshold = pconsistency.threshold` evaluated 1 candidate where the drivers' explicit `stop_threshold = NULL` evaluates all 749 | The drivers disable early stop deliberately (the OC family needs every candidate's Pcons); any revisit is a semantics decision, not an optimisation |
| 4 | The search phase's own single in-process `future_lapply` | 6.6 % (0.31 s) [gbsg 3.4 %; anchor 10.4 %, 0.85 s] | ≤ 0.3–0.9 s per call | Small; grows with M (anchor), so worth a look only after ranks 1–2 |
| 5 | Batch size for **direct** parallel calls (§6; not a bootstrap item) | 63 rounds vs 6: 20.2 → 12.8 s at workers = 6 | ≈ 7.4 s per application-style multisession call | Low mechanically, but the default batching exists for early stopping granularity; measurement only here |
| — | Cuts, enumeration, dedup/selection, GRF (at the `NULL` default) | ≤ 0.2 % each | nothing | — |

**Where the measurement stands against chat's re-scaling arithmetic** (the task's §7 asks for this explicitly):

1. **"The fit is now roughly half a continuous replicate" — close, slightly understated, and it misses the operative fact.** Measured: 61.9 % of the profiled call, 58 % of the 5.02 s replicate — more than half, not roughly half. But the share was never the actionable number: **only ~5 % of that bucket is arithmetic**. The re-scaling could not have said this — it is a composition fact the old profile did not resolve — and it reverses the natural reading. "The fit dominates" suggests the estimator is expensive; the measurement says the estimator is nearly free and its scaffolding (subsetting, `model.frame`, `summary`) is the cost.
2. **"Fits plus `survfit` medians are ~80 % of a survival call" — understated: measured 87.9 %** (52.4 + 35.5). Direction and rank order correct.
3. So the re-scaling's *ranking* survives measurement — no contradiction in what comes first — but both shares were low by 6–10 points (the removed overhead was slightly larger than the arithmetic assumed relative to the remaining work), and the split inside rank 1 is the finding the arithmetic could not produce.

## 8. Close-out

`git status --porcelain` before staging: this task's own artefacts only. Staged by explicit path: the script, log, `.rds`, three `.Rprof` files, this report. **No push. No install. No `R/` change.** Script wall-clock ≈ 4.5 min (three fixtures × 3 calls, the 100 s bootstrap, two multisession calls).

Commits: `d7d5aae2` task doc; the artefacts-and-report commit follows (hash in `git log`).

## 9. Verdict (ten lines)

1. Installed 0.3.2 (`fbd564de` dispatch fix and `eb136a35` orientation both in the log), `vi.grf.min` default `NULL`; all timings at that default, sequential plain-loop configuration.
2. The continuous replicate-configuration call is **29.42 → 4.74 s**; the dispatch bucket fell 25.17 → 0.31 s and every other bucket's absolute seconds are unchanged within sampling noise (traced fit totals 2.94 → 2.96 s).
3. Re-shared, the continuous call is **fit 61.9 % · consistency evaluator 29.6 % · future 6.6 %**; gbsg is **fit 52.4 % · medians 35.5 % · evaluator 8.1 %**; the anchor (4 702 fits, 2.4× M) is **fit 88.5 %** — buckets converge to fit as M grows. Bucket-sum gate passed (99.8 % of wall attributed, ≤ 1.5 % unassigned).
4. Inside the fit bucket the **solve is 2–7 %**; 93–98 % is machinery — subsetting/assembly and `model.frame`/`model.matrix` on the continuous path, `summary.coxph` (40 %) on the survival path.
5. An unadjusted fast path is routable by a construction-time predicate on both paths: `ps_adjust_method == "none" && length(.fs_adjust_terms(adjust_covariates)) == 0L` (continuous, closure build at `forestsearch_main.R:1721`); `is.null(adjust_covariates)` (survival, call site `subgroup_search.R:657`). Quoted, not designed.
6. Medians: computed for all 1 410 gbsg fits, consumed (dedup key cols 6–7; `m1.threshold` filter, `Inf` default nothing sets) only on the 121 effect-screen survivors — **ratio 11.7×, prize ≈ 2.5 s ≈ 32 % of the call**; both consumers could be served from the later set.
7. Bootstrap check: 20 replicates, **5.02 s each** — 1.059× the profiled call, 1.015× the dispatch task's 4.95 s; **B = 1000 = 83.7 min sequential**, confirming ≈ 82.
8. Batch size (direct multisession calls, workers = 6): default 12 → 63 rounds, 20.21 s; one batch per worker (125) → 6 rounds, 12.80 s. Measurement only.
9. Chat's re-scaling: ranking confirmed, shares understated (fit 62 % not ~50; gbsg fit+medians 88 % not ~80); its real gap is composition — the fit bucket is machinery, not arithmetic.
10. No `R/` change, no optimisation written; one script, one `.rds`, three `.Rprof`, one report; anchor reproduced its document's subgroup exactly.
