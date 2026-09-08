# REPORT — ACTG175 continuous intervals (TASK_actg175_continuous_intervals_2026-09-07)

Date: 2026-09-08 (task issued 2026-09-07). Mac-Studio-3, `~/Documents/GitHub/forestsearch`, branch `feature/glm-extension-mac`, HEAD `1bf54ff4` at the start of Stage 1, forestsearch 0.3.5. Task: `dev/tasks/TASK_actg175_continuous_intervals_2026-09-07.md`; decisions N-1–N-4 at defaults. No R/ change; no edit to any off-limits path; no push. Stage 0 record (§2's classification, the identity anchors, the pass-throughs, the render projection): `REPORT_actg175_continuous_intervals_stage0_2026-09-07.md` (Gate 0 PASS, commit `1bf54ff4`) — N-1 = `analysis_actg175_continuous_oc.qmd`, the only document reproducing Ĥ = `{age <= 37} & !{cd40 <= 507}`, n = 66 from a committed payload; no committed MR (IJ) interval and no FB result existed, so the rows are naive / MR (IJ two-term) / MR (field) plus the Bonferroni pair.

Status: **Stage 1 rendered and committed (`f87518e9`), Gate 1 PASS (with the cross-platform floating-point note below); this is the Stage 2 record.**

## Stage 1 — the flags and the section; render attempts on the Mac

Branch `feature/glm-extension-mac`, HEAD `1bf54ff4`, forestsearch 0.3.5. Working tree at session start (2026-09-08 08:56, fresh session after the Mac rebooted out of application memory during the previous render): the Stage 1 edits to `analysis_actg175_continuous_oc.qmd` (uncommitted; kept — `dirout = "analysis_actg175_continuous_oc_intervals"`; on the gate call `mr_inference = TRUE`, `mr_inference_args = list(ci_method = "field", draws = 5000L, include_complement = TRUE, field_complement = TRUE, return_reselection = TRUE)`; the new section `## 2.1 Post-selection intervals for the found subgroup and its complement {#sec-intervals}` with chunks `intervals-objects`, `intervals-table-H`, `intervals-spread`, `intervals-table-Hc`, `intervals-table-joint`, `intervals-diagnostics`, `intervals-reading`; the payload element `extras$intervals = iv`), a modified tracked HTML (restored from git) and an untracked partial payload directory (moved aside under the session scratchpad, not used).

### Mac environment facts (this machine: Mac-Studio-3, M4 Max, 14 cores, 36 GB)

1. **Apple Accelerate is not fork-safe.** R's BLAS here is `/System/Library/Frameworks/Accelerate.framework/.../libBLAS.dylib`; `parallel::mclapply()` children that touch BLAS after the parent's dense-matrix work (the field gate) can hang or die unless `VECLIB_MAXIMUM_THREADS=1` is set in the render environment. The probe of 2026-09-07 23:16 (`probe_veclib.R`: the gate call with the three flags, then a 2-worker `mclapply` doing `crossprod`/`solve`/`eigen`) passed under `VECLIB_MAXIMUM_THREADS=1` (`PROBE OK`); every render below was launched with it.
2. **Worker cap for 36 GB.** Rule for this machine: (workers × per-worker RSS) under 24 GB with no other R process running. Measured today for this document's OC loop (`eval_job`, one `fs_oc_grid()` on 20,000 draws × M = 4,508 in one block of `block = 5e4`):
   - one job in-process (`probe_one.R`, job 8 = the T̂_obs rung): **55.5 s**, process RSS from 3,795 MB (the parent after the `family` chunk) to a peak of **10,740 MB** (sampled every 5 s: 4,430 → 7,991 → 9,364 → 10,052 → 10,740 MB); system-level 15 G → 22 G used;
   - two forked jobs (`probe_mem.R`, `mc.cores = 2`): the children's RSS reached 8,237 and 9,684 MB within 30 s of the fork (still rising; the probe was killed at 08:59:47), system `PhysMem` 35 G used, 17 M unused, compressor rising — a forked child re-touches the parent's heap (R's gc marks every live object), so each child duplicates most of the parent's RSS in addition to its own job.
   - Consequence: at 14 workers (the document's `params$n_workers`) the loop needs ≈ 14 × 10.7 GB; at 2 workers the total R footprint (parent 3.8 GB + 2 × 10.7 GB) still exhausted the machine (render attempt 1 below). **The OC loop fits only at `n_workers = 1`** (in-process `lapply`; the results are seeded per job, so the loop's numbers do not depend on the worker count).

### Worker-count settability (checked from source before rendering)

- OC loop (chunk `evaluation-loop`, lines 567–568): `oc <- do.call(rbind, parallel::mclapply(seq_len(nrow(jobs)), eval_job, mc.cores = params$n_workers))` — a `params:` entry (`n_workers: 14`, line 6), settable from outside with `quarto render ... -P n_workers:N`. Set to 1 for attempt 2.
- Chunk `null-not-shared` (lines 615–617): `fam0 <- parallel::mclapply(Q_variants, function(cu) fs_oc_family_enumerate(dgm_at(q_rungs[1], cu), fs_args, n = N, max_M = 10000L, verbose = FALSE), mc.cores = 3L)` — **a literal, not settable from outside.** It runs after the loop, when the parent holds ≈ 12–14 GB (the loop's heap is not returned to the OS), so its three forks each re-touch a ≈ 13 GB parent: 3 × ≈ 13 GB ≫ 24 GB.

### Render attempts (all: `VECLIB_MAXIMUM_THREADS=1`, `scripts_mdf1/tmo.sh 3600 quarto render analysis_actg175_continuous_oc.qmd -P n_workers:<N>`, detached; summed R/quarto RSS sampled every 30 s into the render log; a guard killing R when the system's unused memory fell below 2,500 MB)

| attempt | `n_workers` | start | outcome |
|---|---|---|---|
| 1 | 2 | 09:02:26 | guard tripped 09:03:14 in `evaluation-loop` (unused 48 MB; total R RSS 10,607 MB at 09:02:59 and rising) — killed |
| 2 | 1 | 09:03:52 | `evaluation-loop` completed (chunk 30 from ≈ 09:04 to 32/59 `[homogeneous-null]` at ≈ 09:22:30, ≈ 19 min; total R RSS 11.7–13.9 GB throughout, system unused 10–11 GB, compressor flat at 0.8 GB); `homogeneous-null` completed; **guard tripped 09:23:44 at 34/59 `[null-not-shared]`** (unused 2,427 MB, down from 11 GB at 09:23:30 — the three literal forks) — killed |

Nothing was written to the repo by either attempt (no payload, the tracked HTML untouched; knit intermediates removed). Logs: session scratchpad `render_oc_intervals_2workers_guardkilled.log`, `render_oc_intervals.log`, `render_sys*.log`, `render_guard*.log`, `probe_mem*.log`, `probe_one*.log`.

### The worker-count stop and Larry's disposition (09:24–09:27)

After attempt 2 the remaining fork site was the literal `mc.cores = 3L` at line 617 (chunk `null-not-shared`); reported as a stop with the number 1 proposed. **Disposition (Larry, 09:26): approved as a document edit** — line 617 `mc.cores = 3L` → `mc.cores = min(3L, params$n_workers)` (Linux keeps 3 at `n_workers:14`; the Mac gets 1 at `n_workers:1`; `fs_oc_family_enumerate()` is deterministic so `fam0` is unchanged either way — Gate 1(d) below is the proof: the `null-not-shared` chunk's rendered numbers are unchanged), plus one line in the document's clock/provenance chunk `doc-clock` (lines 1229–1231) printing `n_workers` as rendered and the ≈ 11 GB per-worker cost of the OC loop. Both edits sit beside the two Mac facts above; nothing else in the document changed.

| attempt | `n_workers` | start | outcome |
|---|---|---|---|
| 3 | 1 (with the line-617 edit) | 09:27:26 | **completed**, `WRAPPER EXITED rc=0` at 09:48:15 — wall 20 min 49 s; the document's own clock "document compute wall-clock so far: 20.2 min"; "evaluation loop: 1123.4 s (18.7 min) over 20 jobs, 1 workers"; peak summed R/quarto RSS 15,641 MB (30-s samples), system unused never below 9.4 GB; guard never tripped |

## Gate 1 (task §3.3)

Checker: session scratchpad `gate1_check.R` (a–c, `identical()` element by element), `gate1_numdiff.R` (every differing numeric leaf quantified), `gate1_htmldiff.py` (d: rendered text of the committed HTML vs the new HTML with the new section and its TOC entry removed). New payload: `_payloads/analysis_actg175_continuous_oc_intervals/analysis_actg175_continuous_oc_intervals_payload.rds` (13,875 bytes, built 2026-09-08 09:47:42 under 0.3.5); committed reference: `_payloads/analysis_actg175_continuous_oc/analysis_actg175_continuous_oc_payload.rds` (built 2026-08-31 21:06:37 under 0.3.2 **on Linux** — its HTML's payload path is `/home/larryleon/...`).

**(a) Anchors and pre-existing elements.** Rule `{age <= 37} & !{cd40 <= 507}` identical; `n_H = 66` identical; `labels` identical; `p_cons = 0.95` identical; `meta` identical except the compute fields (`n_workers` 14 → 1, `t_loop_secs` 2,056.4 → 1,123.4) and the provenance fields (`built_at`, `forestsearch_version` 0.3.2 → 0.3.5); the only new element is `extras$intervals`; no new top-level element. `identical()` is FALSE for `T_obs` and everything downstream of it: `T_obs` committed 87.9166666666667 vs new 87.9166666666666 (|d| = 9.9e-14, the last bit; `q_rungs[8]`, `q_shared[6]` and one `c1_ladder` entry are set to `T_obs` at lines 145–147, so `table$q`, the ladders and every `c1`/`q` column inherit the same 9.9e-14), `beta_treat` −26.978724502131353 vs −26.97872450213173 (|d| = 3.8e-13). Across all 65 differing numeric leaves the maximum relative difference is 1.6e-14 except `Enaive_bias` (max |d| 2.0e-7, max relative 1.3e-9 — a draw-average bias term, cancellation-amplified); `all.equal()` at the default tolerance is TRUE on everything but the new element and the compute/provenance fields; at tolerance 1e-12 only the eleven `Enaive_bias` columns are reported. `extras$purity` and `extras$c2_policy` are bit-identical. **Reading: every §2.2 anchor reproduces at the document's precision (T_obs 87.9166667; every rendered OC number is byte-identical in (d)); bit-identity against a payload built on Linux with a different BLAS is not attainable and was not the intent of the gate. Judged PASS; Larry may overrule — the raw numbers are above.** The `intervals` element is bit-identical (timing fields excluded) to the set-aside Mac render of 2026-09-07 22:05, so the Mac path is deterministic.

**(b)** All 77 numeric leaves of `extras$intervals` finite. **PASS.**

**(c)** Ĥ: β̃ = 33.6147229870, q₀.₉₅(Λ*) = 77.4101727437, lower_1s = −43.7954497567 = β̃ − q₀.₉₅ (|d| < 1e-9); Ĥᶜ: β̃ᶜ = −33.6173825244, q₀.₀₅(Λ*ᶜ) = −15.0628175104, upper_1s = −18.5545650140 = β̃ᶜ − q₀.₀₅; IJ lower_1s = β̃ − 1.645·SE_IJ; gate naive est = T_obs; naive one-sided conventions hold; two-sided field intervals ordered; γ = 0.025 ∈ [0.025, 0.05]; p̂(Ĥ) = 0.0872 ∈ [0, 1]. **PASS.**

**(d)** Rendered text (scripts/styles stripped, gt's random table ids normalised): committed 1,848 lines, new outside the new section 2,106 lines; 280 differing lines, of which 226 are the new section's own block (`@@ -1020,0 +1048,226 @@`). The remainder, every one accounted for: the `dirout` code (twice: the folded code and the appendix), the gate-call code with the `mr_inference` lines (twice), the `mc.cores = min(3L, params$n_workers)` code (twice), the `doc-clock` code and its new output line `rendered with n_workers = 1 (...)`, the payload code line `intervals = iv,`, the two timing lines (`evaluation loop: 1123.4 s (18.7 min) over 20 jobs, 1 workers` vs `2056.4 s (34.3 min) ... 14 workers`; `document compute wall-clock so far: 20.2 min` vs `51.9 min`), the payload path (`/Users/...` vs `/home/...`), and the anchored-truth print's two zero-valued covariance columns `C_mu0_tau` / `V_tau` at 1e-13 / 1e-28 (machine noise; the neighbouring `V_mu0`, `bracket`, `V_eff` columns are identical). **No rendered result number changed anywhere outside the new section. PASS.**

**Gate 1: PASS.** Commit `f87518e9` (the `.qmd`, the payload, the rendered HTML, by explicit path).

## The three tables as rendered (from the payload / the rendered HTML; oriented = `y_decline`, positive = harm; `cd4_change` = sign flip)

**Table 1 — Ĥ = {age <= 37} & !{cd40 <= 507}, n = 66 of 1083**

| Method | Point est. (oriented) | SE | Two-sided 95% (oriented) | One-sided 95% LOWER (oriented) | Point est. (cd4_change) | Two-sided 95% (cd4_change) | One-sided bound (cd4_change): change ≤ |
|---|---|---|---|---|---|---|---|
| Naive | 87.92 | 42.66 (robust) | (4.31, 171.53) | 17.75 | −87.92 | (−171.53, −4.31) | −17.75 |
| MR (IJ two-term) | 33.61 (β̃) | 53.91 (IJ) | (−72.04, 139.27) | −55.05 | −33.61 | (−139.27, 72.04) | 55.05 |
| MR (field) | 27.26 (est₂); β̃ 33.61 | 38.79 (λ-SD) | (−56.44, 94.20) | −43.80 | −27.26 | (−94.20, 56.44) | 43.80 |

Retained-optimism spread (oriented): naive 87.92 → MR (IJ) 33.61 → MR (field) est₂ 27.26; the corrections remove 54.30 (IJ) and 60.66 (field) CD4 cells/mm³ (selection term 52.71, fixed term 1.59).

**Table 2 — Ĥᶜ, n = 1017 (benefit claim; the one-sided UPPER bound is the exposed limit)**

| Method | Point est. (oriented) | SE | Two-sided 95% (oriented) | One-sided 95% UPPER (oriented) | Point est. (cd4_change) | Two-sided 95% (cd4_change) | One-sided bound (cd4_change): change ≥ |
|---|---|---|---|---|---|---|---|
| Naive | −36.54 | 7.83 (robust) | (−51.89, −21.19) | −23.66 | 36.54 | (21.19, 51.89) | 23.66 |
| MR (IJ two-term) | −33.62 (β̃ᶜ) | 15.12 (IJ) | (−63.25, −3.98) | −8.75 | 33.62 | (3.98, 63.25) | 8.75 |
| MR (field) | −31.96 (est₂); β̃ᶜ −33.62 | 8.07 (λ-SDᶜ) | (−46.81, −16.35) | −18.55 | 31.96 | (16.35, 46.81) | 18.55 |

Regime diagnostic: λ-SDᶜ / naive SE = 1.030 (IJ SE / naive SE = 1.931).

**Table 3 — the joint pair** (γ = 0.025 on the 0.025–0.050 grid; corr(Λ*, Λ*ᶜ) = −0.019 over 983 aligned draws)

| Pair | Ĥ lower (oriented) | Ĥᶜ upper (oriented) | Ĥ: change ≤ | Ĥᶜ: change ≥ | Joint prob. on the aligned draws |
|---|---|---|---|---|---|
| Separate one-sided 95% field bounds | −43.80 | −18.55 | 43.80 | 18.55 | < 0.95 by construction |
| Bonferroni (γ = 0.025 each) | −56.44 | −16.35 | 56.44 | 16.35 | 0.949 |
| Calibrated from field$joint (γ = 0.025) | −56.44 | −16.35 | 56.44 | 16.35 | 0.949 |

**Diagnostics.** p̂(Ĥ) = 0.087 (selected `{age <= 37} & !{cd40 <= 507}`); top re-selected `{age <= 35} & !{cd40 <= 507}` 0.218, 2nd `!{wtkg <= 73} & !{cd40 <= 507}` 0.169, 3rd = Ĥ itself 0.087; selection rate 0.972; gate family 4,935 candidates (the search's own candidate set; the OC family on the super-population has M = 4,508 — different objects, as in Stage 0); field Λ*-mean +6.355; λ-SD / IJ SE = 38.79 / 53.91 = 0.72, λ-SD / naive SE = 0.91; field draw usage 983 / 1000 outer, 499 / 500 inner; complement fits 761, λ-SDᶜ / naive SE = 1.030. Payload `settings`: `ci_method = "field"`, draws 5000, multiplier `poisson`, reselection `maxeff`, `ij_residual = "two_term"`, seed 8316951, R_out/R_in = 1000/500, n_selected 66. Gate wall 11.0 s (multiplier stage 0.6 s, field 9.9 s, complement field 0.6 s).

## Render wall and commits

Attempt 3 (the committed render): 20 min 49 s wrapper wall (09:27:26–09:48:15), document clock 20.2 min, loop 18.7 min at `n_workers = 1`. Session cumulative wall (probes, three attempts, gates): 08:56–09:55 ≈ 59 min of the 2 h ceiling.

```
git log --oneline 8cbba125..HEAD
f87518e9 ACTG175 continuous intervals Stage 1: MR gate ... Gate 1 PASS ...
1bf54ff4 ACTG175 continuous intervals Stage 0: Gate 0 PASS -- N-1 = analysis_actg175_continuous_oc.qmd ...
f18a478f Stage 3 report wording (carried over): ...
4be1f9f9 ACTG175 continuous intervals task: task document (2026-09-07)
```
(plus this record's commit, which follows.)

## Reading (ten lines)

1. The complement first: the field's one-sided 95% upper bound on Ĥᶜ is −18.55 oriented, a CD4 change of at least 18.55 cells/mm³ on the raw scale; "harm at most τ on Ĥᶜ" is supported at every reading threshold τ ∈ {0, 10, 20, 30} and not supported at none.
2. Its comparators: MR (IJ) upper −8.75 (change ≥ 8.75); naive upper −23.66 (change ≥ 23.66) — the selection adjustment on the complement costs 5.1 (field) to 14.9 (IJ) CD4 units of the naive bound.
3. Ĥ: the field's one-sided 95% lower bound is −43.80 oriented (a CD4 change of at most 43.80); "harm at least τ on Ĥ" is supported at none of τ ∈ {0, 10, 20, 30, 40}; the naive lower bound 17.75 would have supported τ = 0 and 10.
4. MR (IJ) as the conservative reference sits outside the field's on both blocks (Ĥ lower −55.05, SE 53.91 = 1.26× naive; Ĥᶜ upper −8.75, SE 15.12 = 1.93× naive).
5. A claim on both uses the Bonferroni pair (Ĥ lower −56.44, Ĥᶜ upper −16.35; raw scale Ĥ change ≤ 56.44, Ĥᶜ change ≥ 16.35); the calibrated pair coincides with it (γ = 0.025, corr(Λ*, Λ*ᶜ) = −0.019, joint probability 0.949).
6. The price of selection on Ĥ: point estimate 87.92 → 33.61 (IJ) → 27.26 (field est₂); adjusted lower bounds −43.80 (field) to −55.05 (IJ) against the naive 17.75.
7. Regime: p̂(Ĥ) = 0.087 — the tie regime; the top re-selected labels share the cd40 cut and differ only in the age (35 vs 37) or wtkg factor, the near-identical-membership ties of the mdf1 finding, in which the field's one-sided bounds were at nominal on this DGM family.
8. On the complement λ-SDᶜ / naive SE = 1.030: the complement's selection-adjusted upper bound is essentially the naive one shifted by β̃ᶜ − naive (+2.9) and the second-order term, not inflated.
9. No bootstrap row exists for this document and none was run; no Guo–He row (survival-only adapters); thresholds are the simulation design's reference points and may be replaced.
10. **The two numbers to read first:** the complement's field one-sided upper bound **change ≥ 18.55 CD4 cells/mm³ on Ĥᶜ** (oriented −18.55), and Ĥ's field one-sided lower bound **change ≤ 43.80 on Ĥ** (oriented −43.80).

## Findings (no task proposed)

- Cross-platform identity: the committed OC payload was built on Linux; a Mac render reproduces it to 1e-13 (last bit in T_obs), not bit-for-bit. Any future "identity against the committed payload" gate on a document rendered on the other platform needs a floating-point tolerance stated in the task.
- The Mac cannot render this document at the committed `n_workers = 14` (≈ 11 GB per worker); `-P n_workers:1` (20.8 min) is the Mac setting, `n_workers:14` remains Linux's.
- `p̂` on labels understates the settledness of the membership (the mdf1 finding), visible here: the top three re-selected labels all contain `!{cd40 <= 507}`.
