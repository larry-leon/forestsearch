# CC task — survival medians on the passing set: move the `survfit` computation behind the effect screen, under a bit-identity equality guard

**File:** `dev/tasks/cc_task_medians_on_survivors_2026-09-02.md` · **Issued:** 2026-09-02 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first (§1). If absent under the exact name, the hyphen-stripped stem `cc_task_medians_on_survivors_20260902.md` is the same file.
**Reads (in-repo):** `dev/glm-continuous-sims/REPORT_bootstrap_reprofile_2026-09-01.md` §4(b) (the measurement and the two consumer sites), `dev/glm-continuous-sims/REPORT_fastpath_unadjusted_2026-09-02.md` (the 0.3.3 state of `fit_cox_for_subgroup`, the volatile-field exclusion list, the recorded check/test baselines), `dev/glm-continuous-sims/fastpath_baseline_2026-09-02.R` (the battery pattern and fixture blocks to reuse).

**Why.** The 09-01 re-profile measured that survival candidate medians (`survfit` + `summary.survfit`) are computed inside the per-candidate fit for **every candidate whose Cox fit succeeds** — 1,410 on gbsg — while both consumers act only on `hr.subgroups`, the effect-screen survivors: 121 rows. Ratio **11.7×**; the medians bucket is **2.69 s, 35.5% of the gbsg call**; computed on the 121-row set it would be ≈0.23 s — **prize ≈2.5 s per gbsg call, ≈32%**. The 09-02 fast-path report sharpened the stakes: the survival fit-machinery prize was a stack-attribution artifact and realized only ~0.1–0.2 s, so **this move is essentially the whole survival lever**. The anchor showed the row-builder ratio at 32.2× for general M; the continuous path sets `m1 = m0 = NA` and pays nothing — it must remain bitwise untouched.

**The shape of the change, from the rank-2 risk note:** *a move, not a removal.* Every row that exists today exists because its candidate passed the effect screen; computing medians only for passers therefore changes **no value on any row** — it elides computation whose results never reach any output. That is the equality argument, and the bit-identity gates check it rather than assume it.

**The governing rule:** a performance change that alters a result is a defect regardless of its speed. Guard: `identical()` — no tolerance.

---

## ⚠ CATEGORY — this task edits `R/`

- **Moves existing code** (the substance): the `survfit` medians computation relocates from inside `fit_cox_for_subgroup()` (currently ≈L761–763 pre-0.3.3 numbering; re-anchor at HEAD) to the post-screen point of the same per-candidate iteration. The computation itself — call and extraction — is preserved operation for operation.
- **Adds:** only the minimal plumbing the move requires (a small local frame rebuild at the new site; `NA` medians in the fit's return, patched for passers).
- **Changes behaviour:** none on any output surface. Non-passing candidates lose a computation whose results were never consumed; passing candidates get identical values. Adjusted and unadjusted survival arms are **both** in scope — the medians lines sit after the arm branch and move for both.
- **Changes the method:** no.

Files expected to change: `R/subgroup_search.R` only, plus `DESCRIPTION` (version → 0.3.4) and `NEWS.md`. The two consumer sites — the dedup key in `R/subgroup_consistency_helpers.R` and the `m1.threshold` filter in `R/subgroup_consistency_main.R` — are **read-only for this task**: they must keep receiving identical frames from untouched code. If implementation genuinely requires touching them or any other `R/` file, **STOP** and report why.

**Compute:** Stage-A battery ≈5–6 min (dominated by a 20-replicate gbsg bootstrap, ≈150 s); two install/test cycles; `R CMD check` ≈9 min; Stage-C battery ≈4–5 min; verification + one profiled gbsg call ≈3 min. **Estimate 45–75 minutes wall total.** No simulation campaign, no renders.

**Unattended.** Gates stop the task; never ask, never work around. On any equality failure: report the case, the component, and the first differing values verbatim, leave the tree committed and documented, stop. Provenance gate is dirt-tolerant; stop only if dirt touches this task's own paths.

---

## 1. Provenance and first commit — GATE

Standard block: machine, repo path, branch, HEAD, `git status --porcelain`, `packageVersion("forestsearch")` (expect **0.3.3**), R version, cores. *GATE:* branch `feature/glm-extension`; the four fast-path commits `c79f797d`, `b9d7b622`, `4dbf9f26`, `5d8f7f68` in the log. Copy this document into `dev/tasks/`, commit it **alone**, report both filenames and the hash. Record the `vi.grf.min` default in force (expect `NULL`).

## 2. Re-verify the design basis at HEAD — GATE

Line numbers have moved (0.3.3 added ≈96 lines to `subgroup_search.R`); re-anchor everything. Quote line-anchored at current HEAD:

1. **`fit_cox_for_subgroup()` as it now stands** — the full body: signature, the `df.x` assembly, the 0.3.3 unadjusted fast branch and the adjusted `summary.coxph` branch, the `try()` contracts, and **the medians computation**: the exact `survfit` call (formula, `data =`, any arguments), the exact extraction of `m1`/`m0` (which object, which indices or names — if extraction indexes `$table` by **row name**, say so: it decides the §4 frame-rebuild requirement), and precisely which inputs it consumes (`df.x` columns and nothing else?).
2. **The per-candidate sequence in the search loop** — from the Cox fit's return to row creation: where the **effect screen** is applied (the `hr.threshold` / `minp` conditions), as which status/branch, and whether the screen verdict is known **per candidate, inside the iteration, before the row is created**. Quote the status codes and the row-builder call. *GATE:* the design below (§4) requires the verdict per-candidate in-loop, with `yy`, `dd`, `tt`, `id.x` still in scope at that point. If the source shows the screen is instead applied later, vectorized on an assembled frame, **STOP and report the actual structure** — the move then needs a different design and goes back to chat.
3. **Row provenance:** confirm from source that rows reaching `hr.subgroups` exist **only** for screen-passers (the record's counts — gbsg 1,410 fits → 121 rows; continuous 1,975 → 852 — say so, but verify the code path, not the counts). This is the premise of the equality argument.
4. **Every caller of `fit_cox_for_subgroup()`** in the package, by file:line, and what each consumes from its return. *GATE:* no caller outside the search's per-candidate iteration consumes `med0`/`med1`. If one does, STOP and report it.
5. **The two consumer sites, re-quoted at HEAD:** the dedup key (`subgroup_consistency_helpers.R`, the `cols_to_check <- 2:min(10, ncol(df))` block) and the `m1.threshold` filter (`subgroup_consistency_main.R`, the `is.finite(m1.threshold)` block), plus the formal's `Inf` default and the two live call sites passing `Inf` explicitly. Confirm unchanged in substance from the 09-01 report's quotes.
6. **The no-other-reader survey, re-run at HEAD:** grep the package for every read of `m1`/`m0`/`med0`/`med1`; classify each as (consumer-of-`hr.subgroups`) or (display-for-selected-subgroup) or (other). *GATE:* nothing reads medians for candidates outside `hr.subgroups`. Quote anything in the "other" class and stop if it exists.
7. **Where the fixture arguments live**, unchanged: reuse the F1/F2/F3/F5 blocks and the volatile-field exclusion from `fastpath_baseline_2026-09-02.R` verbatim.

## 3. Stage A — pre-change baselines, at HEAD, before any edit

One script, `dev/glm-continuous-sims/medians_baseline_2026-09-02.R`, same architecture as the fast-path battery (same stage-switch mechanism, same pruning, **the same volatile-field exclusion list — reused verbatim, not re-derived**). Cases:

**Core fixtures:**
- **F2 gbsg** — the primary; this is where the prize lives.
- **F5 gbsg adjusted** (`adjust_covariates = "age"`) — the adjusted arm is touched by the move; it needs its own equality row.
- **F1 continuous MD40** and **F3 anchor** — untouched-path guards: the continuous fast path must remain bitwise inert (`m1 = m0 = NA` throughout).

**Edge fixtures (synthetic survival, fixed seeds, small):**
- **E-ties** — reuse the fast-path task's heavy-ties construction (integer times 1:6, n = 120): `survfit` on tied data, moved call must match.
- **E-named** — a set where at least one **passing** subgroup's arm never reaches median survival (heavy censoring in one arm, engineered so the candidate still clears `hr.threshold` and the event floors): `m1` or `m0` is `NA` today; the moved computation must produce the identical `NA`, and the dedup key must paste the identical key string.
- **E-finite** — a gbsg-configuration call with **`m1.threshold` set finite** (e.g., 60 months): no live driver exercises this consumer, so this task must — the filter path runs, and selection/frames must be `identical()` pre/post. Choose a value that actually excludes at least one row at 0.3.3 (verify in Stage A; adjust the value if not, and record it).
- **E-zero** — a configuration whose effect screen passes **nobody** (e.g., `hr.threshold` set prohibitively high on the E-ties data): `hr.subgroups` is empty or absent; the moved medians step must no-op identically — same no-find return, no error, no warning delta.
- **E-degen** — the monotone-likelihood construction from the fast-path battery, provided its degenerate candidate **passes** the screen: warnings and counters at parity, medians computed for it identically. If it cannot pass, adjust thresholds so it does, and record the adjustment.

**End-to-end:** a **20-replicate gbsg bootstrap** (`forestsearch_bootstrap_dofuture`, `nb_boots = 20`, fixed seed, sequential) — the survival analogue of the fast-path task's F1 bootstrap; this is the pre-change wall-clock and payload baseline both.

Save the pruned comparison sets to `medians_baseline_2026-09-02.rds`, with the log; commit script + artefacts as the second commit.

## 4. Implement — the design

**One move, arm-agnostic, inside the per-candidate iteration:**

1. `fit_cox_for_subgroup()` no longer computes medians: the `survfit` lines are removed from it, and its return carries `med0 = NA_real_, med1 = NA_real_` (both arms; the rest of the body, including the 0.3.3 fast branch and both `try()` contracts, byte-for-byte untouched apart from the removed lines).
2. At the post-screen point of the iteration (the §2.2-verified location where the candidate is known to pass and `yy`, `dd`, `tt`, `id.x` are in scope), compute medians **for passers only**, by rebuilding the minimal frame and replicating the removed lines verbatim:

```
idx    <- id.x == 1
med_df <- data.frame(Y = yy[idx], E = dd[idx], Treat = tt[idx])
<the exact survfit call the removed lines made, with data = med_df>
<the exact m1/m0 extraction, operation for operation>
```

   The frame rebuild exists precisely so the `survfit` formula and any **name-indexed** extraction (`"Treat=0"` / `"Treat=1"` row labels — per the §2.1 finding) see identical variable names and therefore identical labels. If §2.1 shows the current call's `data` argument or column names differ from `Y`/`E`/`Treat`, match **whatever the current lines use** — the removed code is the specification, not this sketch.
3. The computed `m1`/`m0` enter the row exactly where the fit's values entered it before — patch the result list before the row builder, or pass alongside, whichever requires fewer touched lines; the row's fields, names, types, and order are unchanged.
4. Non-passers: no medians call, and — per §2.3 — no row, so no output surface exists on which they could differ.

**Housekeeping in the same commit:** `DESCRIPTION` → **0.3.4**; `NEWS.md` entry under 0.3.4 (performance: survival candidate medians computed only for effect-screen survivors — a relocation of the existing computation, verified bit-identical; ~11.7× fewer `survfit` calls on the gbsg application); if `fit_cox_for_subgroup()`'s roxygen mentions medians, update it truthfully (`@noRd` internals, markdown-on conventions); `devtools::document()` clean; tidyverse style. Commit as the third commit; reinstall.

## 5. Stage C — the equality gates, post-change

Rerun the same script against 0.3.4 → `medians_postchange_2026-09-02.rds`; compare with `identical()` under the reused exclusion list:

- *GATE E-1 (F2, F5):* every retained component identical — `hr.subgroups`/`out.found` at full precision (all columns, `m1`/`m0` included), selected subgroup, consistency outputs, warnings/counters.
- *GATE E-2 (F1, F3):* the continuous guards identical — the move must be invisible there.
- *GATE E-3 (E-ties, E-named, E-finite, E-zero, E-degen):* identical per case, including the `NA` key strings (E-named), the finite-filter selection (E-finite), and the empty no-op (E-zero).
- *GATE E-4 (gbsg bootstrap):* the 20-replicate payload identical against Stage A, volatile fields excluded.

On any failure: case, component, first differing values verbatim, both sides, full precision — then STOP.

## 6. Elision proof — the `survfit` count

Equality proves harmlessness; this proves the move happened.

- With 0.3.4 installed, count `survival::survfit` dispatches (`trace()` counter) during one F2 call and one F5 call. Expected: ≈ the passer count (121-ish for F2's configuration) plus attributable display-site calls for the selected subgroup — attribute every count above the passer set by source line. Repeat at 0.3.3 (or cite Stage A instrumentation) for the before-count ≈1,410+. Report the before/after table.
- One line per fixture confirming zero `survfit` calls on F1/F3 (continuous) at both versions.

## 7. Realized recovery — measure, don't extrapolate

- Re-profile **one F2 single call** under `Rprof(interval = 0.010, line.profiling = TRUE)`, the established procedure and buckets: medians-bucket seconds and share before → after (2.69 s / 35.5% at 0.3.2–0.3.3 → expect ≈0.2–0.3 s), total gbsg wall before → after (≈7.4–7.6 s → expect ≈5 s).
- The Stage-C gbsg bootstrap: per-replicate mean against Stage A's, and the B = 1000 projection for the survival configuration, before → after.
- One line comparing realized to the predicted ≈2.5 s, and one line stating where the residual gbsg samples now sit (expected: coxph internals + assembly — the material for the unblocked follow-on in §9). A shortfall is a finding for the record, not a trigger for edits.

## 8. Package health gates

- `devtools::test()` on 0.3.4 source: 0 failures; warning-count parity with the recorded 0.3.3 baseline (31); `test-search-reproducibility.R` passes **unmodified**. *GATE.*
- `R CMD check` (`devtools::check(document = FALSE, vignettes = FALSE, manual = FALSE)`): **0 errors / 0 warnings / 0 notes**, matching the recorded 0.3.3 result. *GATE.*

## 9. Out of scope — explicit

- **The now-unblocked assembly skip.** With medians no longer consuming `df.x`, the unadjusted fast branch could fit on vectors directly (the fast-path task's unexercised tier). It is a **separate change**: recorded here, decided in the workstream after §7's residual numbers land. Do not implement it.
- **The IPTW continuous defect** (recorded at 0.3.3): untouched, its own future task.
- **Deeper Cox surgery, `stop_threshold` semantics, batch-size and `parallel_args` defaults:** recorded items, untouched.
- **No renders, no push, no payload, application, or driver changes.**

## 10. Report

`dev/glm-continuous-sims/REPORT_medians_on_survivors_2026-09-02.md`:

1. Provenance; the §2 re-verification with HEAD-anchored quotes — the medians lines and their extraction mechanics, the per-candidate screen location and status codes, the row-provenance confirmation, the caller list, the two consumer sites, and the no-other-reader survey's classification.
2. The move as implemented: exact lines removed, exact lines added, file:line, and the before/after of the roxygen if touched.
3. Stage A/C mechanics: the reused exclusion list (stated verbatim), the equality matrix — every fixture × every retained component — and the E-case findings, including the E-finite value chosen and what it excluded, and any threshold adjustments made for E-degen.
4. The `survfit` count table with residual-call attribution.
5. §7's bucket line, walls, bootstrap per-replicate and B = 1000 before/after, realized-vs-predicted, and the residual-cost statement.
6. §8's results against the recorded baselines.
7. Ten-line verdict.

Commits, in order: (1) this task doc alone; (2) Stage-A script + artefacts; (3) the implementation + version/NEWS/docs; (4) Stage-C/verification artefacts + profile outputs + the report. Explicit paths; tree clean at close. **No push. No render. Nothing outside `R/subgroup_search.R`, `DESCRIPTION`, `NEWS.md`, and `dev/` artefacts.**
