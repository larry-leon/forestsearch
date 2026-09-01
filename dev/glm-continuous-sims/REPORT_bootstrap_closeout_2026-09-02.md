# REPORT — bootstrap-efficiency close-out: two verifications, production-shape confirmation, and the register stamped

**Task:** `dev/tasks/cc_task_bootstrap_closeout_2026-09-02.md` · **Executed:** 2026-09-02 (session date 2026-09-01) by CC
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension` · **Read-only: zero edits outside `dev/`.**

**Verdict up front.** Verification A: the replicate payload consumes only winner-level content plus named search-stage bookkeeping, and a caller's `stop_threshold` rides into every replicate unmodified — early-stop enablement is a driver-arg decision, no wrapper change needed. Verification B: direct-call results are `identical()` between sequential and multisession plans for both outcome types — the seeding-unification item closes as moot. Production shape: nb = 100 payloads identical sequential vs 48 workers; continuous B = 1000 ≈ 12.5 min (linear upper bound; ≈ 2–3 min amortized), survival ≈ 3.9 min at 48 workers; +≈ 17.7 GB for 48 workers on this machine. Register items E, F, and memo §4 stamped closed by citation; deep Cox surgery recommend-closed; the binary limb recorded as handed off.

---

## 1. Provenance

pop-os, `~/Documents/GitHub/forestsearch`, branch `feature/glm-extension`, HEAD at task start `e31e612d`, porcelain clean. `packageVersion("forestsearch")` **0.3.5** ✓; R 4.6.1, 128 cores. *GATE:* assembly-skip commits `4b826c98`, `934f31ed`, `f3975b99`, `e31e612d` in the log ✓. `vi.grf.min` default **`NULL`** ✓. Task doc committed alone: `dev/tasks/cc_task_bootstrap_closeout_2026-09-02.md`, commit `91268179`.

## 2. Verification A — replicate consumption and `stop_threshold` forwardability

### 2.1 What the replicate consumes from the inner `forestsearch()` result

Every read of `run_bootstrap` inside the `%dofuture%` body (`R/bootstrap_analysis_dofuture.R`):

- **Identification gate and flag:** `!is.null(run_bootstrap$sg.harm)` (L626) gates the bias-correction block; `any_found` is derived from the same signal (L939–943); DINA mode reads `sg.harm` again for fallback labels (L843–867, consistency mode leaves it untouched).
- **Winner membership frames:** `df_PredBoot <- run_bootstrap$df.predict; dfboot_PredBoot <- run_bootstrap$df.est` (L629–630) — the H*/Hc* splits (`treat.recommend`) that feed all four bias-correction fits (`get_Cox_sg` / `fit_subgroup_effect`, L687–773). Winner-level only.
- **Search-stage bookkeeping:** `fg <- run_bootstrap$find.grps`; `max_sg_est`, `max_count`, `L` (L638–642) → payload columns. These come from the **search** (`format_search_results()`), which runs before and independently of the consistency stage — untouched by consistency early stop.
- **The consistency table — first row only:** L797–828 reads `run_bootstrap$grp.consistency$out_sg$result`, takes `first_row <- sg_result[1, ]`, and extracts `Pcons, hr, N, E, K, g, m, M.1–M.7`. **No other row, no per-candidate vector, no count from the table is read**; `grp.consistency`'s bookkeeping (`n_candidates_evaluated`, `n_passed`, `early_stop_triggered`, `early_stop_candidate`, `stop_threshold` echo) is never consumed.
- **GRF cuts:** `run_bootstrap$grf_cuts` (L898–902); the MR object as an attribute only when `mr_in_replicates` (L962+).

The stored payload row (L905–960) is built from exactly these plus the replicate's own H/Hc fits, event counts, and volatile `tmins_*`.

**Verdict (i): winner-level content plus named bookkeeping.** The bookkeeping fields that would legitimately differ under early stop are precisely the **un-consumed** `grp.consistency` fields (`n_candidates_evaluated`, `n_passed`, `early_stop_triggered`, `early_stop_candidate`, the `stop_threshold` echo) and the row count of `out_sg$result` beyond row 1 — nothing that reaches the payload. The consumed `find.grps` fields (`max_sg_est`, `L`, `max_count`) are search-stage and stop-invariant. Given Settlement B's row-level result (the winner row `identical()` between early-stop and full evaluation, `g`/`m` included — `REPORT_assembly_skip_2026-09-02.md` §6), **early stop under `maxeffCons` is payload-invariant on this branch's wrapper**: a bookkeeping change, not a content change.

### 2.2 Forwardability

`args_FS_boot <- args_FS_template` (L558), where `args_FS_template <- fs.est$args_call_all` (L406), and `args_call_all <- mget(names(formals()), ...)` captures every `forestsearch()` formal — `stop_threshold` included (`R/forestsearch_main.R:1372–1373`). The wrapper's per-replicate mutations (L559–607) touch only: `df.analysis`, `df.predict`, `details`, `show_candidate_summary`, `plot.sg`, `plot.grf`, `quiet`, `grf_res`, `grf_cuts`, `dina_res`, `dina_cuts`, `ps_hat`, `parallel_args${plan,workers,show_message}`, and conditionally `mr_inference`. A grep of both bootstrap files finds **zero** `stop_threshold` occurrences (the doc line aside): **a caller-supplied `stop_threshold` rides through to the inner `forestsearch()` unmodified.** The only rewrite sites live inside `forestsearch()` itself: the `sg_focus = "maxeff"` override (L1513–1519) and the non-`maxeffCons` soundness guard (L1608–1643), both of which NULL it and sync `args_call_all` — for `maxeffCons` neither fires. **Enablement is driver-arg-only; no wrapper edit required.** Two corollaries for the decision: (a) since the formal's default is `stop_threshold = pconsistency.threshold`, a `maxeffCons` analysis that never passes the argument already early-stops in its replicates today; the drivers that pass `stop_threshold = NULL` explicitly (the committed fixture blocks do) are the ones opting out; (b) what's actually being decided is whether those drivers drop their explicit `NULL`.

### 2.3 Branch caveat

This verdict holds for `feature/glm-extension`'s wrapper. The `feature/mr-in-replicates` branch adds in-replicate machinery with its own consumption; the early-stop decision must be re-checked at merge time. (Recorded; no work on that branch here.)

## 3. Verification B — direct-call plan reproducibility

`closeout_planrepro_2026-09-02.R` (fixture blocks and pruner verbatim from the committed battery; F1 continuous `seedit = 8316952`, F2 gbsg `seedit = 8316951`): each run once under `parallel_args = list(plan = "sequential", workers = 1L)` and once under `list(plan = "multisession", workers = 6L)`, pruned and compared with `identical()`. The two calls differ in exactly one argument by construction, and `forestsearch()` echoes its arguments in `args_call_all`, so the `args_call_all$parallel_args` echo is excluded before comparison (it is the experimental manipulation, not a result — the exclusion is in the committed script with that comment; nothing else is excluded beyond the battery's volatile list).

| fixture | sequential wall | multisession/6 wall | selected | pruned results `identical()` |
|---|---|---|---|---|
| F1 continuous | 2.31 s | 20.86 s | `!{wtkg <= 84} & !{cd40 <= 320}` (both) | **TRUE** |
| F2 gbsg survival | 3.35 s | 6.47 s | `{er <= 0} & {size <= 35}` (both) | **TRUE** |

**Verdict (i): identical, per outcome type — the seeding-unification item is moot; closed.** (Consistent with the mechanism: the search loop's `future_lapply` evaluation is RNG-free per candidate, and the consistency stage re-seeds deterministically.) Side observation for the record: a parallel plan on a **direct** call of this size only costs (F1 2.3 → 20.9 s — globals export dwarfs the 2-second search); outer-plan parallelism belongs at the bootstrap wrapper level, not inside single calls.

## 4. Production-shape confirmation and campaign guidance

**Transplant:** `closeout_parallel_2026-09-02.R` = `assembly_settleA_2026-09-02.R` with only the named lines changed — `nb_boots = 40L → 100L` (and its two `wall / 40` per-replicate divisors + the label print), `workers = 20L → 48L` (and its label/print occurrences), output filename → `closeout_parallel_2026-09-02.rds`. Diff verified to contain nothing else (shown in the commit-2 log and reproducible via `git diff` of the two committed files). *GATE ✓*

**Continuous F1, nb = 100, seed 8316952 (0.3.5):**

- Sequential: **204.9 s** (2.05 s/replicate). `multisession, workers = 48`: **75.2 s**. **Payloads `identical()` — TRUE.** *GATE ✓*
- Throughput and amortization: ideal 48-worker compute for 100 replicates ≈ ⌈100/48⌉ × 2.05 ≈ 6 s; observed 75.2 s ⇒ **fixed startup + globals-export overhead ≈ 69 s** (48 worker sessions spun up, package load, per-chunk serialization). Implied **B = 1000 at 48 workers: 12.5 min by linear extrapolation — an upper bound**; with the ≈ 69 s treated as fixed, the amortized estimate is ≈ 69 s + 1000 × 2.05/48 ≈ **2 min** (sequential comparator: 34.1 min).
- **Memory** (`closeout_parallel_mem_2026-09-02.log`, `free -m` sampled every 10 s): baseline before launch 32.6 GB used / 224.9 GB available (of 257.5 GB); peak during the 48-worker phase **50.4 GB used** (207.2 GB still available) — **≈ +17.7 GB for 48 workers, ≈ 370 MB/worker** on this fixture. Worker counts on this machine are memory-uncritical at this footprint; the survival OC-grid memory note (~5–6 GB/process) is the binding case, not these bootstraps.

**Survival F2, nb = 100, 48 workers (wall only):** **23.3 s** (0.23 s/replicate effective) ⇒ **implied survival B = 1000 at 48 workers ≈ 3.9 min** (sequential comparator from the assembly-skip record: 69.7 min). Reproducibility follows from the same mechanism — main-process index matrix + per-replicate deterministic seeding — plus Settlement A and this section's continuous `identical()`; the ~7-minute sequential comparator was not paid, per the task.

**Campaign guidance (record only — no driver or document edited):**

- Invocation, as the wrapper's source defines it: `forestsearch_bootstrap_dofuture(fs.est, nb_boots = B, seed = <seed>, parallel_args = list(plan = "multisession", workers = W, show_message = FALSE))` — the plan enters at `setup_parallel_SGcons(parallel_args)` (`bootstrap_dofuture_main.R:366`; allowed plans `multisession` / `multicore` / `callr` / `sequential`, workers capped at `detectCores()`).
- Measured throughputs at W = 48, nb = 100: continuous 75.2 s (≈ 69 s of it fixed startup), survival 23.3 s. B = 1000: continuous ≈ 2 min amortized (≤ 12.5 min linear bound), survival ≈ 3.9 min.
- Startup amortization means small nb undersells parallelism (Settlement A's 2.6× at nb = 40) — size B or batch campaigns so the ~1 min spin-up is paid once per launch, not per small run.
- Inner searches stay sequential **by construction** (`args_FS_boot$parallel_args$plan <- "sequential"`, `bootstrap_analysis_dofuture.R:595`) — outer workers do not nest.
- Memory: ≈ 370 MB/worker observed at 48 workers on the continuous fixture (+17.7 GB total, 207 GB headroom on this machine); worker count can go well above 48 here before memory matters for bootstraps.
- Payload identity to sequential is proven at nb = 40/20w and nb = 100/48w (continuous) and follows mechanically for survival; a campaign can adopt the parallel plan without changing any recorded result.

## 5. The register, stamped

1. **Memo item E (combination-index hoisting): closed — dead by measurement.** Enumeration/floors bucket: 08-30 report bucket table — "enumeration / floors (excl. the fit): 0.03 % (0.01) | 0 | 0.03 % (0.01) | 0.1 % (0.01) | 0.2 % (0.06)" across the five profiled configurations (`REPORT_bootstrap_profile_2026-08-30.md`, bucket table); 09-01 tables — 0.01 s / ≤ 0.2 % per fixture and "Cuts, enumeration, dedup/selection, GRF … ≤ 0.2 % each — nothing" (`REPORT_bootstrap_reprofile_2026-09-01.md` L62/L76/L88/L218); 0.3.4/0.3.5 re-bucketing — enumeration 0 samples (`assembly_rebucket_034/035_2026-09-02.log`).
2. **Memo item F (membership/floor checks): closed — dead by measurement.** Its targets live inside the same enumeration and dedup/selection buckets, ≤ 0.2 % per call at every profile through 0.3.5 (same citations as item E).
3. **Memo §4 (the `"split"` consistency fallback inside replicates): closed.** "Literal-split fallback: fired 0 times in 1 616 evaluations" (`REPORT_bootstrap_profile_2026-08-30.md`, §counters table and finding 7); "the closed-form consistency never returned `NA` and the literal-split fallback fired 0 times (as at 08-30)" (`REPORT_bootstrap_reprofile_2026-09-01.md` L103).
4. **Memo item A, binary limb: handed off, not open here.** Deferred by design to the OR workstream: "binary/count outcomes are **out of scope** (§9): the closure fast variant is built only under `outcome_type == "continuous"`; no glm-family fast path of any kind" (`dev/tasks/cc_task_fastpath_unadjusted_2026-09-02.md` L94; NEWS 0.3.3: "binary/count outcomes are out of scope"). IRLS bit-identity is uncertifiable at this layer; glm.fit-level routing is that track's call.
5. **Deep Cox surgery (`coxph.fit`/`agreg.fit`): recommend-closed as mooted by reproducible outer parallelism.** The residual per-call cost at 0.3.5 is the fit bucket's 2.29 s on gbsg (wrapper/internals 1.26 s, model-frame 0.60 s, solve 0.26 s — `assembly_rebucket_035_2026-09-02.log`); §4's campaign configuration covers it in wall-clock (survival B = 1000 ≈ 3.9 min at 48 workers). Larry may reopen; the recommendation and the residual seconds are the record.
6. **Stale doc line `bootstrap_analysis_dofuture.R:200`** ("Sets reproducible seeds: 8316951 + boot * 100" — no such code exists): assigned to the queued pre-CRAN housekeeping pass; untouched here per the task's read-only category.

## 6. The decision hand-back

Two decisions now sit with Larry; findings only, the recommendation layer is chat's:

- **Early-stop enablement in bootstrap replicates.** §2's verdict: payload-invariant under `maxeffCons` on this branch (winner-level consumption only; the fields that differ are un-consumed bookkeeping), and `stop_threshold` forwards through the wrapper unmodified — **driver-arg-only if commissioned** (concretely: the drivers' explicit `stop_threshold = NULL` is the opt-out; the formal default already early-stops). Caveat: re-check at the `feature/mr-in-replicates` merge (§2.3).
- **Medians removal (NA-out + display recompute).** Hygiene, not speed — the medians bucket is ~0.22 s per gbsg call at 0.3.4/0.3.5 (`REPORT_medians_on_survivors_2026-09-02.md` §7; `assembly_rebucket_035_2026-09-02.log` 0.23 s). If commissioned: a separate task — NA-out in the search payload plus display-site recompute, per the medians record.

## 7. Ten-line verdict

1. Read-only close-out executed: zero edits to `R/`, `DESCRIPTION`, `NEWS.md`, tests, drivers, or documents; writes confined to `dev/tasks/` and `dev/glm-continuous-sims/`.
2. Verification A(i): the replicate consumes winner-level content only — `sg.harm`, the two membership frames, the first row of `out_sg$result`, `grf_cuts` — plus stop-invariant search bookkeeping (`max_sg_est`, `L`, `max_count`); the fields early stop would change are never consumed.
3. `stop_threshold` rides from `args_call_all` into every replicate unmodified; enablement needs no wrapper edit — the drivers' explicit `NULL` is what's actually being decided; mr-in-replicates re-check recorded.
4. Verification B(i): direct-call sequential vs multisession/6 results `identical()` for F1 continuous and F2 survival (argument echo excluded) — the seeding-unification item closes as moot.
5. Production shape confirmed: nb = 100 continuous payloads `identical()` sequential vs 48 workers (204.9 → 75.2 s).
6. B = 1000 at 48 workers: continuous ≈ 2 min amortized (12.5 min linear bound), survival ≈ 3.9 min measured at nb = 100 (23.3 s) — the whole survival campaign cost collapses from the 69.7-min sequential projection.
7. Startup ≈ 69 s fixed at 48 workers; ≈ 370 MB/worker, +17.7 GB total, 207 GB headroom — memory-uncritical for bootstraps on this machine.
8. Register stamped: items E and F and memo §4 closed by measurement with committed citations; the binary limb handed to the OR workstream; deep Cox surgery recommend-closed with residual seconds stated; the stale doc line assigned to pre-CRAN housekeeping.
9. Direct-call parallel plans cost rather than help at fixture scale (F1 2.3 → 20.9 s) — outer-plan parallelism belongs at the wrapper, a guidance line for drivers.
10. Deviation of process: the §3 script's comparison initially flagged the `args_call_all$parallel_args` argument echo; the committed script excludes that one echo with an explanatory comment (commit 2 amended) — no result field is excluded.

## Commits

1. `91268179` — task doc alone (`dev/tasks/cc_task_bootstrap_closeout_2026-09-02.md`).
2. `e964b129` — scripts: `closeout_planrepro_2026-09-02.R`, `closeout_parallel_2026-09-02.R` (transplant), `closeout_parallel_surv_2026-09-02.R` (amended once: the argument-echo exclusion in the §3 script's verdict).
3. (this commit) — outputs + report: `closeout_planrepro_2026-09-02.{log,rds}`, `closeout_parallel_2026-09-02.{log,rds}`, `closeout_parallel_mem_2026-09-02.log`, `closeout_parallel_surv_2026-09-02.{log,rds}`, this report.

No push. Tree clean at close.
