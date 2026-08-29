# REPORT — profiling `forestsearch()` to locate the bootstrap's cost

**Task:** `dev/tasks/cc_task_bootstrap_profile_2026-08-30.md` (commit `9d53f465`)
**Repo / branch:** `~/Documents/GitHub/forestsearch` (pop-os), `feature/glm-extension`
**Category:** measurement only. No `R/` change, no optimisation, no edit to any estimator, driver or document. Artefacts, all under `dev/glm-continuous-sims/` (chosen over `dev/profiling/` to sit beside the OC-wrapper reports the fixture comes from): `bootstrap_profile_2026-08-30.R`, `bootstrap_profile_2026-08-30.log`, `bootstrap_profile_2026-08-30.rds`, five raw `bootstrap_profile_<config>_2026-08-30.Rprof` sample files, this report. No renders, no push, no install.

**Chat's memo.** `memo_bootstrap_performance_2026-08-30.md` is **not on this machine** — not in the repository, `~/Downloads`, or anywhere under `$HOME` (`find` by name). The only record of its ranking available here is the task's own §3 bucket list, which names the suspected hot spots in this order: per-candidate effect fit, candidate medians, GRF variable importance, cut construction, enumeration/floors, consistency screen, dedup/selection. §6 ranks against *that* list and says where the measurement contradicts it; if the memo ranked differently, this report cannot know.

---

## 1. Provenance (§1 — GATE passed)

```
pop-os · /home/larryleon/Documents/GitHub/forestsearch · feature/glm-extension
HEAD at task start fd115e3f (the residual task's last commit; this task ran second in the session)
git status --porcelain at task start: the profile script and its in-flight log only (this task's own scratch, untracked); nothing else
git log -4: fd115e3f 243801e6 80544be1 5ff1487a
packageVersion forestsearch 0.3.0 · R 4.x · 128 cores, 251 GB
```

`~/Downloads/cc_task_bootstrap_profile_2026-08-30.md` was absent under that name; the stem arrives hyphen-stripped as **`cc_task_bootstrap_profile_20260830.md`** (plus an identical `(1)` copy) → `dev/tasks/cc_task_bootstrap_profile_2026-08-30.md`, committed alone as **`9d53f465`**.

**`vi.grf.min` default in force: `NULL`** (`formals(forestsearch)$vi.grf.min`; the queued default task landed as `202ef15b`, 0.3.0). Every timing below states which value the call used.

**Operator note.** Three launches. The first died on an `Rprof` header regex (this R writes `line profiling: sample.interval=10000`) after one instrumented call. The second ran to the end, but I patched the script while it was running — `Rscript` reads a file expression by expression, so the tail was mis-parsed after the bootstrap check had already printed; **every number in that run is superseded**. The third launch is the clean one and is the only one reported; its bootstrap ratio (1.096) agrees with the second's (1.038) to within the replicate-to-replicate spread. Sampling interval: 10 ms is this platform's minimum (`Rprof` warns below it); 451–3 977 samples per call.

---

## 2. Two source questions, answered before profiling

### 2.1 Are `med0` / `med1` consumed anywhere?

They are computed in `fit_cox_for_subgroup()` (`R/subgroup_search.R` L758–773: `summary(survfit(Surv(Y, E) ~ Treat, data = df.x))$table[, "median"]`, "descriptive; treatment-only by design") and set to `NA_real_` on the GLM path (`fit_glm_for_subgroup()` L819–820, "no median survival for GLM outcomes"). `create_result_row()` (L827–841) places them as columns 6 and 7 — `m1`, `m0` — of every candidate row, named in `format_search_results()` L878–879.

Every read of those columns in `R/`:

| site | what it does | selection or display |
|---|---|---|
| `remove_near_duplicate_subgroups()`, `R/subgroup_consistency_helpers.R` L1110 `cols_to_check <- 2:min(10, ncol(df))` → columns `K, n, E, d1, m1, m0, HR, L(HR), U(HR)`, rounded to 0.001 and pasted into the dedup key (L1114–1122); called at `subgroup_consistency_main.R` L579 for every focus except `maxeff` | **`m1`/`m0` are two of the nine key columns** — two survival candidates that agree on `(K, n, E, d1, HR, L, U)` to 0.001 but differ in a median are kept as distinct; on the GLM path both are `NA` and inert | **selection** (survival only) |
| `subgroup_consistency_main.R` L516–531: `if (is.finite(m1.threshold)) { hr.subgroups <- hr.subgroups[!is.na(hr.subgroups$m1), ]; … found.hrs <- hr.subgroups[hr.subgroups$HR >= hr.threshold & hr.subgroups$m1 <= m1.threshold, ] }` | drops candidates whose treatment-arm median exceeds `m1.threshold` | **selection, gated by a formal whose default is `Inf`** (`forestsearch()` L1220; `subgroup.consistency()` L342). No non-archived driver or application sets it finite: the three that name it pass `m1.threshold = Inf` (`quarto/resampling/gbsg_analysis_cox_effMinSG_split_vs_resample.qmd` L128; `quarto/applications/gbsg/analysis_gbsg_survival_multimethod.qmd` L759, L1325) |
| `sort_subgroups()` L540–590 and `sort_subgroups_preview()` L641–690 | keys: `(-HR, K)`, `(-Pcons, -hr, K)`, `(-N, -Pcons, K)`, `(N, …)`, and for `hrMaxSG`/`hrMinSG` `(-in_band, [-on_frontier,] ∓N, -Pcons, -hr, K)` | **no read** |
| the `hrMaxSG`/`hrMinSG` band, `.compute_inclusion_band()` via `forestsearch_helpers.R` L1240–1275 and the sort keys above | effect and size only | **no read** |
| `evaluate_subgroup_consistency()` / `evaluate_consistency_twostage()` (`subgroup_consistency_helpers.R` L1858–1861): `resultk <- c(p.consistency, found.hrs$HR[m], found.hrs$n[m], found.hrs$E[m], found.hrs$grp[m], m, k, mfound)` | carries `HR, n, E, grp` forward; the medians are dropped here | **no read** |
| `bootstrap_summaries_helpers.R` L102–111, L167–170 | `gt` column labels `Med<sub>T</sub>` / `Med<sub>C</sub>` and the "Survival" spanner, when the columns exist | display |
| `forestsearch_cross_validation.R` L1702–1718 | `estimates_toget <- c("Subgroup", "n", "n1", "m1", "m0", "RMST", "HR (95% CI)")` for the CV comparison table | display |
| `fs_mr_oc_summary()`, `fs_mr_inference*.R`, `print_candidate_summary.R`, `forestsearch_methods.R`, `forestsearch_main.R`, the bootstrap replicate loop | — | **no read** (`grep` for `m1`/`m0`/`med0`/`med1` finds only the `m1.threshold` formal and the `oc_analyses.R` arm-mean variables `m1`/`m0`, which are unrelated) |

**Verdict: consumed by selection on two sites, neither of which acts on the paths profiled here.** (i) On the survival path the medians are *part of the near-duplicate key* — removing them would change which candidates the dedup collapses only for pairs identical on the other seven columns to 0.001, which the residual task's analysis of the same key shows means identical sample membership, where the medians coincide anyway; so the key's behaviour would be unchanged in practice, but the dependence is real and must be stated. (ii) The `m1.threshold` filter is a live selection path behind an `Inf` default that nothing in the repository sets. Everywhere else — every sort key, the band, the consistency carry-forward, MR, the bootstrap — they are unread; the remaining reads are two table-formatting sites. On the continuous path they are `NA` and unread by anything.

### 2.2 Does the closed-form consistency ever fall back to `"split"`?

Two kinds of fallback exist in source.

**Global (whole call), `subgroup.consistency()` L386–398:** `consistency_method = "resample"` is replaced by `"split"` with a warning when `consistency_resample()` is not a function in scope, or on the GLM path when `glm_resample_spec` is `NULL`. `forestsearch()` builds the spec for supported measures (MD included), so neither fires on these fixtures.

**Per candidate, when the closed-form rate is `NA`:** `evaluate_subgroup_consistency()` L1450–1466 (Cox) and L1470–1492 (GLM) and `evaluate_consistency_twostage()` L1842–1848 (both families) call `.consistency_via_splits()` — the literal repeated 50/50 split-and-refit, `fs.splits` refits — when `consistency_resample(..., method = "closed")$rate_closed` is `NA`. That is `NA` exactly when the "pieces" builder returns `NULL` (`R/consistency_resample.R` L412–430, L446): Cox — fewer than 2 rows or 2 events (L62), `coxph()` error (L73–74), a missing/`NA` treatment coefficient (L76), or a non-finite or zero `sigma_D = sqrt(sum(dfbeta²))` (L88); GLM — a fit error (L299), a missing/`NA` coefficient or `converged == FALSE` (L304–305), a `.dfbeta_glm()` failure (L309), or non-finite `beta_hat`/`sigma_D` (L320).

**Counted in one search on each fixture** (`trace()` on `consistency_resample()` and `.consistency_via_splits()` in the loaded namespace; sequential plan so every call is in-process):

| fixture | candidates reaching the screen | `consistency_resample()` calls | closed form `NA` | literal-split fallbacks |
|---|---:|---:|---:|---:|
| continuous MD40, n = 500 (`vi.grf.min = NULL`) | 749 | 749 | **0** | **0** |
| continuous, `vi.grf.min = −0.2` | 747 | 747 | 0 | 0 |
| survival gbsg | 120 | 120 | **0** | **0** |

**Verdict: the fallback exists on three code paths and fired on none of 1 616 candidate evaluations.** A candidate that hit it would cost `fs.splits` (400 / 1 000) refits, i.e. ≈ 0.4–4 s on these fixtures per candidate; at the observed rate of zero it contributes nothing to the bootstrap.

---

## 3. Profile

### 3.1 Fixtures and configuration

- **Continuous / MD.** The MD40 DGM the OC-wrapper scripts rebuild (`generate_glm_dgm` on ACTG175 arms 1/3, `k_inter = truth$beta_inter` from the tracked n = 500 payload, gated on `effect_Q` and `prevalence_Q`), one trial `simulate_from_glm_dgm(dgm, n = 500, seed = 8316951 + 1)` (sim_id 1), and the driver's `forestsearch()` arguments verbatim (`sim_fs_maxeffCons_mr_md40_knoise0_n500_batch_1_1000.qmd` L561–590): 13 confounders incl. `str2`, `effect.threshold = 30`, `consistency.threshold = 10`, `pconsistency = 0.90`, `fs.splits = 400`, `n.min = 60`, `d0.min = d1.min = 12`, `maxk = 2`, `sg_focus = "maxeffCons"`, `consistency_method = "resample"`, `conf.cont_jcuts = list(age = 10, preanti = 10)`, `use_twostage = TRUE`, `is.RCT = TRUE`, `seedit = 8316952`.
- **Survival.** `survival::gbsg` (686 rows, 299 events), the effMaxSG application's arguments (`quarto/applications/gbsg/analysis_gbsg_survival_effMaxSG.qmd` L539–585 with its focus profile resolved): confounders `age, meno, size, grade3, nodes, pgr, er`, `hormon` treatment, `time_months = rfstime / 30.4375`, `est.scale = "hr"`, `conf.cont_jcuts = list(er = 10, pgr = 10)`, `conf_force = c("er <= 0", "pgr <= 0")`, `sg_focus = "effMaxSG"`, `n.min = 60`, `d0.min = d1.min = 10`, `maxk = 2`, `hr.threshold = hr.consistency = 1.0`, `pconsistency = 0.90`, `minp = 0.025`, `use_twostage = TRUE`, `consistency_method = "resample"`, `fs.splits = 1000`, `max.minutes = 3`, `vi.grf.min = NULL` (the application does not pass it → installed default).
- **Replicate configuration.** Both are run as `forestsearch_bootstrap_dofuture()` runs each replicate (`bootstrap_analysis_dofuture.R` L573–585): `parallel_args$plan <- "sequential"`, `$workers <- 1L`, `mr_inference = FALSE`, no plots, quiet. The driver's own single call differs in two ways that are each measured separately below: it passes `parallel_args = list(plan = "sequential")` *without* `workers`, and `mr_inference = TRUE` (5 000 draws).
- **Procedure per configuration:** an instrumented call (`trace()` counters and per-fit timing), one un-timed warm-up after `untrace()` (so the JIT recompilation of the restored closures — a once-per-worker cost — is not in the profile), then the profiled call with `Rprof(interval = 0.010, line.profiling = TRUE, memory.profiling = FALSE)`. Buckets are assigned per sample from the whole stack, first match in priority order (medians › GRF › fit › cuts › consistency evaluator › future machinery › dedup/selection › enumeration › rest of `subgroup.consistency()` › everything else), so a `survfit()` under `fit_cox_for_subgroup()` counts as "medians", not "fit", and the consistency evaluator's own time is separated from the `future_lapply` machinery around it.

### 3.2 Wall-clock and the bucket table

| config | `vi.grf.min` | `parallel_args` | MR | instrumented wall | **profiled wall** | samples |
|---|---|---|---|---:|---:|---:|
| cont_viNULL (**replicate configuration**) | `NULL` | `plan = "sequential", workers = 1` | off | 29.55 s | **29.42 s** | 2 940 |
| cont_viNULL_seqloop (**the driver's `parallel_args`**) | `NULL` | `plan = "sequential"` | off | 4.49 s | **4.52 s** | 451 |
| cont_vim02 | `−0.2` | `plan = "sequential", workers = 1` | off | 29.23 s | **29.35 s** | 2 971 |
| surv_viNULL (gbsg) | `NULL` | `plan = "sequential", workers = 1` | off | 8.56 s | **8.19 s** | 817 |
| cont_viNULL_mr (the driver's single call, MR on) | `NULL` | `plan = "sequential", workers = 1` | on, 5 000 draws | 39.92 s | **39.79 s** | 3 977 |

Share of sampled time by bucket (seconds in parentheses):

| bucket | cont, replicate config (29.4 s) | cont, driver config (4.5 s) | cont, `vi.grf.min = −0.2` (29.7 s) | **surv gbsg** (8.2 s) | cont, MR on (39.8 s) |
|---|---:|---:|---:|---:|---:|
| per-candidate effect fit (`fit_glm/cox_for_subgroup` and below, excl. survfit) | 8.8 % (2.60) | **60.8 %** (2.74) | 9.2 % (2.72) | **43.2 %** (3.53) | 6.7 % (2.67) |
| candidate medians (`survfit`, `summary.survfit`) | 0 | 0 | 0 | **31.1 %** (2.54) | 0 |
| GRF variable importance (`causal_forest`, `fit_causal_forest_glm`, `variable_importance`) | 0 | 0 | **1.5 %** (0.46) | 0 | 0 |
| cut construction (`get_FSdata`) | 0 (< 1 sample) | 0 | 0 | 0 | 0 |
| enumeration / floors (excl. the fit) | 0.03 % (0.01) | 0 | 0.03 % (0.01) | 0.1 % (0.01) | 0.2 % (0.06) |
| consistency screen — the evaluator itself (`evaluate_consistency_twostage`, `consistency_resample` and below) | 5.3 % (1.56) | **30.6 %** (1.38) | 4.4 % (1.32) | 7.3 % (0.60) | 3.5 % (1.41) |
| dedup / selection (`remove_near_duplicate_subgroups`, `sort_subgroups*`) | 0.03 % (0.01) | 0.2 % (0.01) | 0 | 0 | 0.03 % (0.01) |
| **future machinery around the consistency screen** (`future_lapply` → `getGlobalsAndPackages` → `findGlobals`, `serializedSize`, future creation/resolution; not the evaluator) | **85.6 %** (25.17) | 7.1 % (0.32) | **84.5 %** (25.10) | **17.7 %** (1.45) | 62.4 % (24.83) |
| everything else | 0.2 % (0.05) | 1.3 % (0.06) | 0.3 % (0.10) | 0.5 % (0.04) | **27.1 %** (10.79) — MR: `crossprod` 7.10 s self |

Top functions (self / total, seconds) — the full top-20 tables for every configuration are in the log and `.rds`:

- **cont, replicate config:** total — `future.apply::future_lapply` 29.33, `future_xapply` 28.95, `getGlobalsAndPackagesXApply` 22.79, `getGlobalsAndPackages` 22.29, `globalsOf` 21.29, `findGlobals` 20.38, `findGlobals_dfs_call` 18.80, `dframe` 16.76, `data.frame` 14.62, `as.data.frame` 10.30; self — `deparse` 1.71, `%in%` 1.70, `data.frame` 1.56, `.deparseOpts` 1.41, `structure` 1.10, `paste` 0.74, `pmatch` 0.73, `serializedSize` 0.70. Not one package function appears in the top 20 by self time.
- **cont, driver config:** total — `search_combinations_parallel` 3.07 (68 %), `future_lapply` 3.06 (the *search's* own future path, evaluating in-process), `run.Future` 2.97; self — `[.data.frame` 0.29, `NextMethod` 0.23, `[[.data.frame` 0.18, `%in%` 0.18, `lm.fit` 0.08, `summary.lm` 0.07: the fit is data-frame subsetting plus `lm`.
- **surv gbsg:** total — `future_lapply` 8.12, `FUN` 7.13, `survival::coxph` 2.27 (27.8 %); self — `order` 0.23, `ifelse` 0.17, `coxph.fit` 0.14, `model.frame.default` 0.10, `[.data.table` 0.10.
- **cont, MR on:** the MR block adds 10.4 s, of which `crossprod` is 7.10 s self (the multiplier-draw products); the rest of the profile is the replicate profile unchanged.

### 3.3 What the 85 % is

`subgroup.consistency()` validates `parallel_args` against `c("plan", "workers")` (`subgroup_consistency_main.R` L663–670). With both present — the replicate configuration — it takes the "parallel" path even for `plan = "sequential"`: `setup_parallel_SGcons()` sets the future plan, `batch_size_parallel <- min(n_workers * 2L, n_candidates)` = **2** (L841; no `stop_threshold`, so neither of the other branches applies), and the loop at L906–924 issues **one `future.apply::future_lapply()` per batch — 375 calls for 749 candidates** — each of which re-derives the closure's globals (`getGlobalsAndPackages` → `findGlobals` walking `eval_fun`'s code with `deparse`, then `serializedSize` on the captured `df`, `found.hrs`, `index.Z`) before evaluating two candidates in-process. That costs ≈ 67 ms per batch, 25 s per call, against 1.6 s of actual consistency evaluation. With `workers` absent — the driver's configuration — the same function warns `parallel_args missing required elements. Using sequential.` and runs the plain loop (L668–669): 4.5 s for the identical search, identical selected subgroup (`!{wtkg <= 84} & !{cd40 <= 320}`), identical counters. On gbsg (120 candidates → 60 batches) the same overhead is 1.45 s of 8.2 s.

The bootstrap sets exactly the expensive configuration on every replicate (`args_FS_boot$parallel_args$plan <- "sequential"; $workers <- 1L`, `bootstrap_analysis_dofuture.R` L573–575). **The replicate therefore costs 6.5× the driver's search on the same data.** The 7.1 % "future machinery" left in the driver configuration is the *search's* own `search_combinations_parallel()` → `future_lapply` path (one call, in-process), not the consistency batching.

### 3.4 Per-candidate fit counts and times (`trace()` on the fit functions; per-call `proc.time()` deltas)

| fixture | candidates reaching a fit | fits returning `NULL` | mean per fit | median per fit | fits total | share of wall | passed the effect screen (`hr.subgroups` rows) | dedup before → after | reaching the consistency screen | passing it |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| cont, `NULL` | **1 975** | 0 | **1.5 ms** | **1.0 ms** | 2.94 s | 9.9 % of 29.6 s / **62 %** of 4.5 s | 852 | 852 → **749** (−103) | 749 | 319 |
| cont, `−0.2` | 1 962 | 0 | 1.4 ms | 1.0 ms | 2.65 s | 9.1 % | 829 | 829 → 747 | 747 | 319 |
| surv gbsg | **1 410** | 0 | **4.6 ms** | **4.0 ms** | 6.48 s | **76 %** of 8.6 s (of which `survfit` ≈ 2.5 s, 31 %) | 121 | 121 → 120 | 120 | 9 |
| cont, MR on | 1 975 | 0 | 1.4 ms | 1.0 ms | 2.78 s | 7.0 % | 852 | 852 → 749 | 749 | 319 |

Two facts worth recording beside the timings. (i) **The near-duplicate dedup removes 12 % of the continuous candidates** (103 of 852) in a sample — identical sample memberships within K under sample-quantile cuts — where the population family of the residual task had 189 identical memberships of 1 885; the mechanism is the one that report characterises. (ii) **749 of 852 effect-screen survivors reach the consistency screen** on the continuous fixture (`maxeffCons`, two-stage, no `stop_threshold`), and 319 pass it; on gbsg 120 reach it and 9 pass. The screen is evaluated on nearly every survivor, so its per-candidate cost — 2 ms of evaluator, 67 ms of batching — is multiplied by ~750.

### 3.5 The `vi.grf.min` toggle (continuous fixture, replicate configuration)

| `vi.grf.min` | forest fits | GRF bucket | profiled wall | candidates fitted | reaching the screen | selected |
|---|---:|---:|---:|---:|---:|---|
| `NULL` (installed default) | 0 | 0 | 29.42 s | 1 975 | 749 | `!{wtkg <= 84} & !{cd40 <= 320}` |
| `−0.2` (the drivers' value) | 1 | **0.46 s, 1.5 %** | 29.35 s | 1 962 | 747 | same |

**The forest is worth ≈ 0.46 s per call, 1.5 % of a replicate at the replicate configuration — and would be ≈ 10 % of the 4.5 s the search costs in the driver configuration.** (The screen at −0.2 is a no-op on membership, as the smoke report found; here it trimmed 13 fits and 2 screen candidates through `max_n_confounders`-independent ordering effects, and selected the same subgroup.) The queued default change buys the bootstrap 1.5 % per replicate as things stand; it is not where the time is.

---

## 4. Extrapolation to `nb_boots = 1000`, and the check

**Projection from the profiled single call** (replicate configuration, `vi.grf.min = NULL`, MR off — the replicate does not carry MR unless `mr_in_replicates = TRUE`): 1 000 × 29.4 s = **29 400 s = 490 min sequential**; on *W* workers ≈ 490 / *W* min (the bootstrap's outer `future` plan; `resolve_bootstrap_parallel_args()` takes the observed fit's `parallel_args` when none are supplied). At the driver's search cost the same bootstrap would be 1 000 × 4.5 s = **75 min sequential**.

**The check:** `forestsearch_bootstrap_dofuture(fs.est, nb_boots = 20, seed = 8316952, parallel_args = list(plan = "sequential", workers = 1))` on the continuous fixture, `fs.est` = the profiled `vi.grf.min = NULL` call (MR off in the fit as well, so `args_call_all` carries the replicate configuration):

| quantity | value |
|---|---:|
| bootstrap wall (20 replicates, sequential) | 644.7 s |
| per replicate (wall / 20) | **32.24 s** |
| the bootstrap's own `tmins_search × 60`: mean / median | 32.09 s / 31.09 s |
| the bootstrap's own `tmins_iteration × 60`: mean | 32.10 s |
| replicates with an identified subgroup | 20 of 20 |
| profiled single call | 29.42 s |
| **ratio per-replicate / single call** | **1.096** |

Within 10 %, i.e. inside the ~20 % band the task set: **the profile's ranking transfers to the bootstrap.** The bootstrap's per-replicate search (32.09 s) accounts for the whole iteration (32.10 s); the four `fit_subgroup_effect()` calls at H/Hc and the subsetting outside the search are ≈ 0.01 s. The 2.8 s excess over the profiled call is the replicate-to-replicate spread of the search (median 31.1 s) plus the bootstrap's `dfnew_predict` scoring; nothing outside the search dominates. Measured projection: 1 000 × 32.24 s = **537 min sequential**.

---

## 5. Ranked by measured share

What is worth changing, ranked by the share of a **bootstrap replicate** it would recover (continuous fixture, replicate configuration, 29.4 s; gbsg in brackets), with the risk in one line. **Nothing is changed here.**

| rank | what | measured share | recoverable | risk note |
|---|---|---:|---|---|
| **1** | The consistency screen's `future_lapply` batching when `plan = "sequential"` — `subgroup_consistency_main.R` L663–670 / L841 / L906–924: with `workers = 1` the screen runs 375 two-candidate `future_lapply` calls, each re-resolving globals; the bootstrap forces this on every replicate (`bootstrap_analysis_dofuture.R` L573–575) | **85.6 %** (25.2 s of 29.4 s) [gbsg **17.7 %**, 1.45 s of 8.2 s] | ≈ 25 s per replicate: the same search runs in 4.5 s under the plain loop the driver already gets — a **6.5× replicate speed-up**, ≈ 490 → 75 min per 1 000 sequential replicates | Low in principle — the plain loop and the batched path evaluate the same closure with the same per-candidate seeds (`.make_candidate_rng_seeds`), but the plain loop's seeding must be checked for replicate-level reproducibility before any switch; the batching also exists for early stopping under `stop_threshold`, which the sequential plan does not need |
| 2 | Candidate medians on the survival path — `survfit()` + `summary.survfit()` inside `fit_cox_for_subgroup()` for every candidate (1 410 calls), consumed only by the dedup key and two display tables (§2.1) | [gbsg **31.1 %**, 2.54 s] — 0 on continuous | ≈ 2.5 s per gbsg replicate, ≈ 55 % of its fit cost | Low: the medians never enter a sort key or the consistency screen; the dedup key would lose two of nine columns, which in practice changes nothing (identical-sample-membership pairs share medians), but that is an argued equivalence, not a measured one, and the display tables would need the medians computed for the *selected* subgroup only |
| 3 | The per-candidate effect fit itself — `lm`/`summary.lm` on a data-frame subset, 1 975 × 1.5 ms (continuous); `coxph()` 1 410 × ~2.5 ms (gbsg, medians excluded) | 8.8 % (2.6 s) [gbsg 43 %, 3.5 s] | Only by a cheaper estimator (e.g. closed-form MD from arm sums instead of `lm` + `summary.lm`; `coxph.fit` on prebuilt matrices instead of `coxph` + `model.frame`): perhaps half of it | Medium: the estimator closure is shared with MR and the bootstrap's H/Hc fits; any replacement must reproduce `estimate`, `se` to the digits the payloads store |
| 4 | The consistency evaluator (`consistency_resample` closed form: one `lm` + `.dfbeta_glm` per candidate, 749 calls) | 5.3 % (1.6 s) [gbsg 7.3 %, 0.6 s] | Marginal; it is already the cheap closed form, 2 ms per candidate | Low; nothing to do until rank 1 is gone |
| 5 | GRF variable importance (`vi.grf.min = −0.2`) | 1.5 % (0.46 s) — already 0 at the installed default `NULL` | The queued default change has recovered this: 0.46 s per replicate, **1.5 %** at the replicate configuration (10 % of the 4.5 s driver-configuration search) | None (landed) |
| — | Cut construction, enumeration/floors, dedup/selection | each ≤ 0.2 % (≤ 0.06 s) | nothing | — |
| — | MR inference (`crossprod` of the multiplier draws) | 10.4 s on the driver's *single* call (26 % of its 39.8 s); **not in the replicates** (`mr_inference <- FALSE` unless `mr_in_replicates`) | irrelevant to the bootstrap as configured | — |
| — | Literal-split fallback | fired 0 times in 1 616 evaluations | nothing | — |

**Where the measurement contradicts the task's (and, as far as it is recorded, chat's) ranking.**

1. **The top item is not on the list at all.** The bucket table the task specified — fit, medians, GRF, cuts, enumeration, consistency, dedup — accounts for **14 %** of a continuous replicate; the other **86 %** is `future`'s globals resolution around the consistency screen, an artefact of `parallel_args = list(plan = "sequential", workers = 1L)` that the bootstrap itself injects. Ranking the listed buckets against each other, as the task asks, would rank the small parts of a small remainder.
2. **The per-candidate fit is first among the listed buckets, but it is not "the" cost.** 2.6–2.9 s per replicate on the continuous fixture, 60 % of the *driver-configuration* search and 9 % of the replicate. Its per-fit cost (1.0 ms median) is close to the floor for `lm` on a data frame; the count (1 975) is the lever, and the count is the enumeration's business, not the fit's.
3. **GRF is 1.5 %, not a leading item**: the queued default change was worth ≈ 0.5 s per replicate. Real, small, done.
4. **Cuts, enumeration and dedup are nil** — below one 10 ms sample each on the continuous fixture. Chat's suspicion of `get_FSdata` and the combination loop does not survive measurement.
5. **The medians are the survival path's genuine second item (31 %)** and are consumed only by the dedup key and display — that part of the memo's list is confirmed, on survival only.
6. **The consistency screen's *evaluator* is cheap (5 %); its *scaffolding* is not.** A memo that reads "consistency screen" as one bucket would attribute 91 % of the replicate to it and draw the wrong conclusion about what to change; the split in the bucket table above is the operative finding.
7. **The literal-split fallback never fires** on either fixture, so the `fs.splits`-refits scenario the task raised is not a contributor.
8. **The bootstrap check passes (1.096)** — nothing outside the search dominates; the replicate *is* the search, and the search inside a replicate is 6.5× the search the driver runs.

---

## 6. Close-out

`git status --porcelain` before staging: the script, log, `.rds`, the five `.Rprof` files and this report; nothing else. Staged by explicit path. **No push. No install** (nothing in `R/` changed). Wall-clock of the clean run: 22 min (four continuous configurations 3 × 5–40 s each, gbsg 3 × 8.5 s, the 20-replicate bootstrap 645 s).

Commits: `9d53f465` task doc; `b8daf89c` script, outputs and this report; the hash line in the follow-up commit.

---

## 7. Verdict (ten lines)

1. Installed 0.3.0, `vi.grf.min` default `NULL`; every timing states its value.
2. `m1`/`m0`: **consumed by selection on two sites** — the survival dedup key (columns 6–7 of nine) and the `m1.threshold` filter behind an `Inf` default nothing sets — and by two display tables; unread by every sort key, the band, the consistency carry-forward, MR and the bootstrap; `NA` and inert on the GLM path.
3. The `"split"` fallback exists on three code paths (whole-call and two per-candidate) and **fired 0 times in 1 616 candidate evaluations**.
4. One continuous replicate (n = 500, MD40, replicate configuration) is **29.4 s**, of which **85.6 % is `future_lapply` globals resolution** — 375 two-candidate batches around a 1.6 s consistency evaluation — because the bootstrap passes `plan = "sequential", workers = 1L`; the driver's `list(plan = "sequential")` fails the `c("plan", "workers")` check and runs the plain loop in **4.5 s** for the identical search.
5. Fits: 1 975 × 1.5 ms (median 1.0 ms) = 2.9 s, 9 % of the replicate, 61 % of the plain-loop search; gbsg 1 410 × 4.6 ms = 6.5 s, of which `survfit` medians are 2.5 s (31 % of the 8.2 s call).
6. Consistency screen reached by 749 / 120 candidates (of 852 / 121 effect-screen survivors; the dedup removed 103 / 1); evaluator 5 % / 7 %.
7. `vi.grf.min = −0.2` costs **0.46 s, 1.5 %** of a replicate; the queued default change has already recovered it. Cuts, enumeration, dedup: ≤ 0.2 % each.
8. Bootstrap check: 20 sequential replicates 644.7 s, **32.2 s per replicate vs 29.4 s single call, ratio 1.096** — the profile transfers; the replicate is the search.
9. Projection: 1 000 replicates = 490–537 min sequential at the replicate configuration; **≈ 75 min if the replicate ran the search the way the driver does.**
10. **Where the measurement contradicts the memo:** the dominant cost is not any of the seven suspected buckets — those sum to 14 % — but the future scaffolding the bootstrap itself configures; the fit is first *among* them at 9 %; GRF is 1.5 %; cuts/enumeration/dedup are nil; medians are real on survival only. No `R/` change, no optimisation written; one script, one `.rds`, five `.Rprof`, one report.
