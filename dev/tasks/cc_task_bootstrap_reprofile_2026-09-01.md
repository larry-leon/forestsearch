# CC task — re-profile the search at 0.3.2: where the bootstrap's cost sits now that the dispatch overhead is gone

**File:** `dev/tasks/cc_task_bootstrap_reprofile_2026-09-01.md` · **Issued:** 2026-09-01 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first (§1).
**Reads (in-repo):** `dev/glm-continuous-sims/REPORT_bootstrap_profile_2026-08-30.md` (the 0.3.0 profile — the baseline this task supersedes), `REPORT_parallel_dispatch_2026-08-30.md` and `REPORT_parallel_dispatch_narrow_2026-08-30.md` (the fix and its measured effect).

**Why now.** The 08-30 profile was taken **before** the dispatch fix, when 85.6% of a continuous call was `future_lapply` globals resolution. That overhead is gone at 0.3.1/0.3.2 (31.1 → 5.4 s per continuous call; 8.9 → 7.4 s on gbsg; 32.2 → 4.95 s per bootstrap replicate; B = 1000 from 490–537 min to ≈ 82 min). Removing overhead does not change the absolute cost of the real work, so **every remaining bucket's share has risen sharply** and the old percentages no longer rank anything. Chat's re-scaling of the old numbers suggests the per-candidate fit is now roughly half a continuous replicate and that fits plus `survfit` medians are ~80% of a survival one — **but that is arithmetic on superseded measurements, not a reading.** This task measures the current state before any optimisation is designed.

The standing lesson from this exact topic applies: **profile before ranking, not after.** No optimisation is written here.

---

## ⚠ CATEGORY

**No `R/` change of any kind.** No edits to any package file, driver, application document or payload. Writes: one scratch script, its profile outputs and logs, and the report, under `dev/glm-continuous-sims/` (or `dev/profiling/` — say which). Plus this task document.

**Compute:** a handful of `forestsearch()` calls under `Rprof`, plus one short bootstrap. **Estimate 30–60 minutes.** No simulation study, no full bootstrap campaign, no renders.

**Unattended.** Gates stop the task; never ask, never work around. **Provenance gate is dirt-tolerant**; stop only if dirt touches this task's own paths.

---

## 1. Provenance and first commit — GATE

Standard block plus the installed version (expect **0.3.2**). *GATE:* branch `feature/glm-extension`, and the 0.3.1 dispatch commits plus the 0.3.2 orientation commit in the log. Copy this document into `dev/tasks/` and commit it alone; report both filenames. **Record the `vi.grf.min` default in force** (expect `NULL` since 0.3.0) — every timing must state it.

## 2. Fixtures

Use the **same two** the 08-30 profile used, so the comparison is like-for-like — quote from that report exactly what they were and reproduce them:

- **Continuous / MD40**, one trial at n = 500, the drivers' `forestsearch_args`.
- **Survival / gbsg-based**, the gbsg application's arguments.

Add one **third, for scale only**: the ACTG175 applied anchor configuration (the fixed-family `maxeffCons` call in `quarto/applications/actg175/analysis_actg175_continuous_oc.qmd` §2, measured at ≈ 8 s sequential). Its candidate space is far larger than MD40's, so it shows how the buckets move with M. One profiled call; no bootstrap on it.

## 3. Profile — the same buckets, re-shared

Profile **one `forestsearch()` call** per fixture with `Rprof(line.profiling = TRUE)` at a fine interval, under the current default `parallel_args` (i.e. the plain loop). Report per fixture:

- total wall-clock, and the top 20 functions by **self** and by **total** time;
- the bucket table, using the 08-30 report's bucket definitions verbatim (per-candidate effect fit · candidate medians (survival) · GRF variable importance · cut construction · enumeration/floors · consistency screen · dedup/selection · everything else);
- **each bucket as a share of the current total**, printed beside the 08-30 share, with a column for the *absolute* seconds then and now. The absolute seconds should be roughly unchanged where the fix did not touch the work — say plainly where they are not;
- candidates reaching a fit, mean and median time per fit, candidates reaching the consistency screen, candidates reaching dedup.

*GATE:* the sum of buckets is within 10% of the measured total; if not, report what is unattributed before drawing conclusions.

## 4. Two source questions for the two live items

**(a) The per-candidate fit — what fraction is machinery rather than arithmetic?** For each outcome type, quote the estimator path taken for the **unadjusted** two-group case (`fit_glm_for_subgroup` / `fit_cox_for_subgroup` and what they call), and from the profile attribute the fit bucket between `model.frame` / `model.matrix` / `vcov` / `summary.*` / formula construction on one side and the actual solve (`lm.fit`, `glm.fit`, `coxph.fit` iterations) on the other. Report the split. Then establish from source **whether an unadjusted fast path can be routed by condition** — is there a single predicate (no adjustment covariates, no weights, no IPTW, no offset) available at the call site that would select it, leaving the adjusted/weighted paths untouched? Quote it or say why not.

**(b) The survival medians — how many candidates actually need them?** The 08-30 report established `m1`/`m0` are consumed, not display-only: the survival dedup key (columns 6–7 of nine) and the `m1.threshold` filter behind an `Inf` default no driver sets. Quote both sites from current source. Then measure the ratio that sizes the opportunity: **how many candidates have medians computed, versus how many reach the dedup step** (the gbsg profile recorded dedup acting on 121 → 120). If medians are computed in the per-candidate row builder for every candidate reaching a fit while dedup sees only the passing set, state that ratio explicitly — it is the size of the prize for computing them on the passing set instead. **Do not design the change; measure the ratio and report whether the two consumers could be served from a later, smaller set.**

## 5. Confirm the extrapolation

Run `forestsearch_bootstrap_dofuture()` at a small `nb_boots` (20–50, state which) on the continuous fixture at 0.3.2. Report the per-replicate mean against the profiled single call and against the 4.95 s the narrow-dispatch task measured, and project B = 1000 against its ≈ 82 min. If the per-replicate mean and the single call disagree by more than ~20%, something outside the search dominates — say so prominently and locate it.

## 6. One recorded item, now cheaply settled

The batch-size question (`min(n_workers * 2L, n_candidates)`) is **no longer a bootstrap item** — after the dispatch fix the bootstrap's inner search takes the plain loop — but it still governs direct `forestsearch()` calls under a parallel plan, which is how the applications run. Settle it with two timed calls on the continuous fixture under `list(plan = "multisession", workers = 6)`: the default batch size, and an explicit `parallel_args$batch_size` large enough to make one batch per worker. Report both wall-clocks and the number of `future_lapply` rounds each implies. **Measurement only** — no change proposed here.

## 7. Report

`REPORT_bootstrap_reprofile_2026-09-01.md`:

1. Provenance, including the `vi.grf.min` default and the installed version.
2. §3's bucket tables, three fixtures, with then/now shares and absolute seconds side by side.
3. §4(a)'s machinery-versus-solve split and the fast-path predicate finding.
4. §4(b)'s quoted consumer sites and the medians-computed / medians-needed ratio.
5. §5's per-replicate check and the B = 1000 projection.
6. §6's batch-size measurement.
7. **A ranked list of what is worth changing**, by *measured current* share, each with the seconds it would recover per replicate, the projected effect on B = 1000, and a one-line risk note. **State explicitly where the measurement contradicts chat's re-scaling arithmetic** (the claim that the fit is now ~half a continuous replicate and that fits plus medians are ~80% of a survival one) — that contradiction, if any, is the most useful thing this report can contain.
8. Ten-line verdict.

Commit the script, its outputs and the report by explicit path. **No push. No install. No `R/` change.**

## 8. Out of scope

- No optimisation, no `R/` edit, no change to any estimator, driver, application or document. The profile decides what gets written; nothing gets written here.
- No full bootstrap campaign — §5 is 20–50 replicates.
- No change to `parallel_args` defaults or batch sizing; §6 is measurement.
- No renders, no push.
