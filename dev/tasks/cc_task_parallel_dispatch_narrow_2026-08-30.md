# CC task — the consistency screen's parallel dispatch, narrowed: a sequential plan takes the plain loop; worker count does not

**File:** `dev/tasks/cc_task_parallel_dispatch_narrow_2026-08-30.md` · **Issued:** 2026-08-30 by chat, commissioned by Larry
**Repo:** `~/Documents/GitHub/forestsearch` (pop-os) · **Branch:** `feature/glm-extension`
**Transport:** `~/Downloads` → `dev/tasks/`, committed first (§1).
**Follows:** `REPORT_parallel_dispatch_2026-08-30.md` (stopped at its §6 gate) and the reverted edit it preserved, `dev/glm-continuous-sims/parallel_dispatch_edit_2026-08-30.patch`.

**Where it stands.** The first dispatch task wrote the short-circuit as *sequential plan OR `workers <= 1`*. Every identity gate passed — C1, C2, C4, C5 `identical()`; C3 the same selection with Monte-Carlo-scale `p.consistency` differences — and the measured effect was 31.1 → 5.5 s per continuous call, 8.9 → 7.5 s on gbsg, 32.2 → 4.97 s per bootstrap replicate. One test failed: `tests/testthat/test-search-reproducibility.R:86`, the documented invariance of `p.consistency` to worker count under `split`. The `workers <= 1` clause sent `multisession, workers = 1` to the plain loop, which seeds a single stream, while `workers = 2` kept the batched path's per-candidate streams. The edit was reverted and preserved as a patch.

**Larry's decision, 2026-08-30: the narrow exit.** The short-circuit fires on the plan alone — `plan = "sequential"` takes the plain loop whatever `workers` says; a parallel plan with one worker stays on the batched path, invariant and slow. The bootstrap's forced `list(plan = "sequential", workers = 1L)` and the driver's `list(plan = "sequential")` both satisfy the plan test, so the whole measured saving is retained. **The test is not edited. Neither path's seeding is touched.** The contract the test asserts — invariance to worker count within a plan — holds for every plan; the sequential plan's single stream under `split` is documented as the one place plans differ.

---

## ⚠ `R/` CALLOUT — **CHANGES BEHAVIOUR**

Same estimand, same results on every `resample` path; a different execution path for a sequential plan; `split` under a sequential plan uses the plain loop's single stream (Monte-Carlo-scale differences from the batched path, as the first report characterised at C3).

Editable: `R/subgroup_consistency_main.R` (Section 6's validation block and the `@param parallel_args` roxygen), `man/`, `NAMESPACE` as regenerated, `DESCRIPTION`, `NEWS.md`. **Nothing else.**

**Explicitly NOT changed:** `tests/` (no test is edited, re-pinned, skipped or given a tolerance); the plain loop's seeding; the batched path's seeding; the batch-size logic at L832–858; `valid_plans`; `R/bootstrap_dofuture_main.R`; any driver, application document, payload or OC-wrapper file.

**Compute:** one extra pre-change capture (§3, under a minute), six post-change captures (§5, a few minutes), the 20-replicate bootstrap re-run (§6, ≈ 2 min), `devtools::test()` and `R CMD check` (≈ 15 min). **Estimate 30–45 minutes.** No simulation study, no renders.

**Unattended.** Gates stop the task; never ask, never work around.

---

## 1. Provenance and first commit — GATE

Standard block (`hostname; pwd; git rev-parse --abbrev-ref HEAD; git rev-parse --short HEAD; git status --porcelain; git log --oneline -8`) plus the installed version (expect **0.3.0**). *GATE:* branch `feature/glm-extension`, clean tree, and the log shows both `5196ea35` (Task 1's close) and `acb49b82` (Task 2's close). Then `git apply --check dev/glm-continuous-sims/parallel_dispatch_edit_2026-08-30.patch` — *GATE:* it must apply cleanly.

Copy this document into `dev/tasks/` and commit it alone; the `~/Downloads` stem arrives with hyphens stripped, so report both names.

## 2. Source verification — GATE

Quote, from the current working tree, not from the first report:

- the test at `test-search-reproducibility.R:86` — both `parallel_args` lists, the `consistency_method` it runs, `fs.splits`, and every expectation it makes;
- the roxygen or `NEWS.md` sentence that documents the invariance the test asserts (the first report called it "documented" — say where);
- the plain loop's seeding (the `use_parallel <- FALSE` branch: where `set.seed()` is called and how many times) and the batched path's per-candidate seeding (how the per-candidate seeds are generated and whether they depend on batch layout);
- how `R/bootstrap_dofuture_main.R` sets the inner `parallel_args`, and what the drivers pass.

*GATE:* the test's two arms must name the **same non-sequential plan** and differ only in `workers`. If either arm is a sequential plan, or the arms differ in plan, the narrow condition cannot pass the test — **stop and report.**
*GATE:* both callers' plans must be the string `"sequential"`, so that `identical(parallel_args$plan, "sequential")` is `TRUE` for each. If a caller passes the plan any other way, stop and report.

## 3. Pre-change capture of a sixth configuration — before any edit

The five configurations' pre-change reference is committed (`parallel_dispatch_pre_2026-08-30.rds`), and `R/` has not changed since it was captured — **do not re-capture it.** Add one configuration, captured with the same script and the same fields:

| # | fixture | `parallel_args` | consistency | early stopping | what it tests |
|---|---|---|---|---|---|
| C6 | continuous | `list(plan = "multisession", workers = 1)` | `resample` | off | **the configuration the narrow condition must leave on the batched path** |

Save it alongside the others and record its wall-clock (expect the batched path's order, ≈ 30 s, not the plain loop's ≈ 5 s).

## 4. The edit

`git apply` the patch, then change the short-circuit so its condition is `identical(parallel_args$plan, "sequential")` and nothing else — delete the `workers` clause. Everything else the patch did stands: placement before the required-elements check, no `warning()`, the `DESCRIPTION` bump, the `NEWS.md` entry (reworded per §7).

Update the `@param parallel_args` roxygen to say: a sequential plan runs the plain loop regardless of `workers`; a parallel plan with `workers = 1` is honoured as parallel (pass `plan = "sequential"` to run the plain loop); `resample` results are identical on every path; under `split` the plain loop uses a single RNG stream, so its `p.consistency` differs from the batched path within Monte-Carlo noise at the stated `fs.splits`, while invariance to worker count holds within any plan; `batch_size` is the knob for genuinely parallel runs.

Quote the full diff of `R/subgroup_consistency_main.R` in the report.

## 5. The identity gates

Re-run C1–C6 post-change with the same script.

**Hard gates — `identical()` on every captured component:**

- **C1, C2, C4, C5** against the committed pre-change reference — as before, `resample` consumes no RNG, so any difference means the edit did more than change which branch runs.
- **C6** against §3 — bit-identical, **and** its post-change wall-clock must be of the batched path's order (≈ 30 s). A plain-loop-fast C6 (≈ 5 s) means the short-circuit still catches a parallel plan — **stop.**
- **C2 additionally:** `early_stop_candidate` and `n_candidates_evaluated` exact.
- **If the first task's post-change capture is in the tree** (`parallel_dispatch_post_2026-08-30.rds`, or whatever name the first report gives it): C1–C5 post-narrow must be `identical()` to it on every component, **C3 included** — the two edits route these five configurations identically, so any difference means this edit did something other than drop one clause. If that file is not in the tree, say so and skip this line.

**C3 — characterised, not gated on identity**, exactly as in the first task: `set.seed()` immediately before each call; identical, or a Monte-Carlo-scale difference with the same selected subgroup, is benign; anything larger, or a changed selection, is a stop.

*GATE:* all of the above.

## 6. Measure what was retained

Per configuration, wall-clock pre and post — the five from the first report, C6 from §3. Then re-run the bootstrap check at the same `nb_boots` the first report used (state it): per-replicate mean against 32.2 s and 4.97 s, and the B = 1000 projection. The narrow condition should retain the whole saving; if the bootstrap does not land near 5 s per replicate, its inner plan is not reaching the short-circuit — report why, and stop.

## 7. Close-out

- `devtools::document()`; `devtools::test()`. *GATE:* **FAIL 0** — in particular `test-search-reproducibility.R:86` passes **unedited**. The first edit's only failure was that test, and this edit's routing changes are a subset of the first's, so any other failure is new and is a stop.
- Warning parity: baseline **31**, and the first task measured no movement. A decrease must be accounted for test by test and attributable only to `"parallel_args missing required elements. Using sequential."`; an increase, or a decrease from any other warning, is a stop.
- `R CMD check` (`RSTUDIO_PANDOC` on Rscript's PATH). Status verbatim; **0/0/0** is the gate.
- `DESCRIPTION` → **0.3.1**.
- `NEWS.md` (rewrite the patch's entry): a sequential plan now runs the plain loop regardless of `workers`; the measured effect (the first report's numbers and this run's); no change to results under `resample`; under `split`, a sequential plan uses a single stream — Monte-Carlo-scale differences from 0.3.0 for callers who passed a sequential plan *with* `workers`; invariance to worker count unchanged within any plan; a parallel plan with `workers = 1` is still parallel; `parallel_args$batch_size` for genuinely parallel runs.
- Commit, **explicit paths only**, never `git add -A`. **No push.** `devtools::install()` (`upgrade = FALSE`); confirm **0.3.1**.
- Report `REPORT_parallel_dispatch_narrow_2026-08-30.md`: §2's quotations · the diff · §5's six results with C3 in full · §6's timings and projection · test and check output raw · the parity accounting · commits · ten-line verdict.

## 8. Out of scope

- No edit to any test — no re-pin, no `skip()`, no tolerance. No change to either path's seeding. No change to the batch-size logic, `valid_plans`, or the bootstrap caller.
- No re-run of the five-configuration pre-change capture, of Task 2, or of anything already committed.
- No push, no simulation study, no renders.
