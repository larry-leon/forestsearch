# REPORT — the consistency screen's parallel dispatch, narrowed: a sequential plan takes the plain loop; worker count does not

**Task:** `dev/tasks/cc_task_parallel_dispatch_narrow_2026-08-30.md` (commit `911d73ba`)
**Repo / branch:** `~/Documents/GitHub/forestsearch` (pop-os), `feature/glm-extension`
**Follows:** `REPORT_parallel_dispatch_2026-08-30.md` (stopped at its §6 gate) and the patch it preserved, `dev/glm-continuous-sims/parallel_dispatch_edit_2026-08-30.patch`.
**Category:** `R/` change, behaviour change — the short-circuit of the first task with its condition narrowed to `identical(parallel_args$plan, "sequential")`. Edited: `R/subgroup_consistency_main.R` (Section 6 and the `@param parallel_args` roxygen), `man/subgroup.consistency.Rd` (regenerated), `DESCRIPTION` (0.3.1), `NEWS.md`. Not touched: any test, either path's seeding, the batch-size logic, `valid_plans`, the bootstrap caller, any driver or document. Scratch: the capture script (extended with C6 and four `*_narrow` tags), its outputs and logs, this report.

---

## 1. Provenance (§1 — GATES passed)

```
pop-os · /home/larryleon/Documents/GitHub/forestsearch · feature/glm-extension
HEAD at start acb49b82 · git status --porcelain: empty
git log -8: acb49b82 2a4c0126 d9ccc24c 5196ea35 774d768c 1e45219f 80ffa2da d77702ce
  -> 5196ea35 (Task 1's close) and acb49b82 (Task 2's close) both present
packageVersion forestsearch 0.3.0
git apply --check dev/glm-continuous-sims/parallel_dispatch_edit_2026-08-30.patch: applies cleanly
```

`~/Downloads/cc_task_parallel_dispatch_narrow_2026-08-30.md` was absent under that name; the stem arrives hyphen-stripped as **`cc_task_parallel_dispatch_narrow_20260830.md`**, with a second copy **`cc_task_parallel_dispatch_narrow_20260830 (1).md`** that differs in two sentences of §1/§4 (it adds "locate the patch the first report preserved … if the report names another path, use that one"; the patch is at the expected path, so the two agree in effect). The plain-named copy → `dev/tasks/cc_task_parallel_dispatch_narrow_2026-08-30.md`, committed alone as **`911d73ba`**.

---

## 2. Source verification (§2 — GATES passed)

**The test**, `tests/testthat/test-search-reproducibility.R` — `.repro_run()` L34–55 and the test L65–95, quoted from the working tree:

```r
.repro_run <- function(workers) {
  fs <- suppressWarnings(forestsearch(
    df.analysis = .repro_df(),
    confounders.name = c("x1", "x2", "x3", "b1", "b2"),
    outcome.name = "y", treat.name = "treat", id.name = "id",
    outcome_type = "continuous", effect_measure = "MD",
    effect.threshold = 2, consistency.threshold = 2,
    pconsistency.threshold = 0.05, fs.splits = 40L,
    n.min = 30, d0.min = 5, d1.min = 5, maxk = 2L,
    sg_focus = "maxeffCons", consistency_method = "split",
    stop_threshold = NULL,   # evaluate every candidate; no early stop
    is.RCT = TRUE, adverse_outcome = FALSE,
    details = FALSE, quiet = TRUE, seedit = 20260827,
    parallel_args = list(plan = "multisession", workers = workers,
                         show_message = FALSE),
    mr_inference = FALSE))
  future::plan(future::sequential)
  fs
}
test_that("search + consistency results are invariant to worker count", {
  skip_on_cran()
  fs1 <- .repro_run(1L)
  fs2 <- .repro_run(2L)
  n_cand <- fs1$grp.consistency$n_candidates_total
  expect_gt(n_cand, 4L)
  pc1 <- fs1$grp.consistency$out_sg$result
  expect_gt(nrow(pc1), 3L)
  expect_true(any(as.numeric(pc1$Pcons) > 0.1 & as.numeric(pc1$Pcons) < 0.9))
  expect_identical(.repro_flat(fs1$find.grps$out.found$hr.subgroups),
                   .repro_flat(fs2$find.grps$out.found$hr.subgroups))
  expect_identical(.repro_flat(pc1),
                   .repro_flat(fs2$grp.consistency$out_sg$result))
  expect_identical(fs1$sg.harm, fs2$sg.harm)
})
```

Both arms: `parallel_args = list(plan = "multisession", workers = workers, show_message = FALSE)` with `workers` = 1 and 2 — **the same non-sequential plan, differing only in `workers`**; `consistency_method = "split"`, `fs.splits = 40`, `seedit = 20260827`; six expectations (candidate count > 4, ≥ 4 result rows, some mid-range `Pcons`, `hr.subgroups` identical, the `Pcons` table identical, `sg.harm` identical). **GATE passed.**

**Where the invariance is documented.** Not in `NEWS.md` and not in a user-facing roxygen sentence: the first report's "documented" pointed at (i) the test file's own header (L1–24: "every candidate's resampling RNG a function of (master seed, global index) only — invariant to worker count and batch size"), (ii) the internal roxygen of `.make_candidate_rng_seeds()` (`R/subgroup_consistency_main.R` L17–38, L31–32: "makes each candidate's RNG depend only on its index and the master seed, invariant to batch size and worker count"), and (iii) the `@param seed` roxygen L181–184: "The seed is used both for sequential execution (via `set.seed()`) and parallel execution (via `future.seed`)" — which already states, in one line, that the two paths seed differently. The 0.3.1 roxygen and NEWS now state the plan-level distinction explicitly.

**The plain loop's seeding** (`use_parallel <- FALSE` branch): `set.seed(seed)` is called **once**, at L471 inside `if (!is.null(seed))`, before Section 8; the sequential loop then draws every candidate's splits from that single stream in candidate order. **The batched path's seeding:** L944 `candidate_seeds <- .make_candidate_rng_seeds(seed, n_candidates)` — L40–65 switch to L'Ecuyer-CMRG, `set.seed(seed)`, and step `nextRNGStream()` once per candidate to build one stream per **global** candidate index, restoring the caller's RNG; each `future_lapply()` batch receives `candidate_seeds[batch_indices]` (L968–969), so a candidate's stream depends on its index and the master seed only, **not on the batch layout or worker count**. That is why arms at `workers` 1 and 2 agree within a plan, and why the plain loop (one stream) differs from either.

**The callers.** `R/bootstrap_analysis_dofuture.R` L573–575: `args_FS_boot$parallel_args$plan <- "sequential"; args_FS_boot$parallel_args$workers <- 1L; …$show_message <- FALSE` — the string `"sequential"`. The drivers (`sim_fs_maxeffCons_mr_md40_knoise0_n{500,700}_batch_1_1000.qmd` L528–532): `inner_parallel <- if (identical(parallel_mode, "sims")) list(plan = "sequential") else list(plan = "multisession", workers = n_workers)` — the string `"sequential"` under `"sims"`. (The gbsg application passes `list(plan = "multisession", workers = n_cores)`, L583 — a parallel plan, unaffected.) `identical(parallel_args$plan, "sequential")` is `TRUE` for both callers' sequential inputs. **GATE passed.**

---

## 3. The sixth configuration, captured pre-change (§3)

`parallel_dispatch_capture_2026-08-30.R pre6` at HEAD `911d73ba`, `R/` unedited (the committed five-configuration reference `parallel_dispatch_pre_2026-08-30.rds` was **not** re-captured):

| # | fixture | `parallel_args` | consistency | early stop | selected (members) | `hr.subgroups` | evaluated / total, passed | wall |
|---|---|---|---|---|---|---:|---|---:|
| C6 | cont | `list(plan = "multisession", workers = 1L)` | resample | off | `!{wtkg <= 84} & !{cd40 <= 320}` (64) | 852 | 749 / 749, 319 | **30.80 s** — the batched path's order |

---

## 4. The edit (§4)

`git apply` of the patch, then the condition narrowed and the roxygen rewritten. Full diff of `R/subgroup_consistency_main.R` against `911d73ba`:

```diff
@@ -230,9 +230,24 @@
 #' @param parallel_args List. Parallel processing configuration:
 #'   \describe{
-#'     \item{plan}{Future plan: "multisession", "multicore", or "sequential"}
-#'     \item{workers}{Number of parallel workers}
-#'     \item{batch_size}{Number of candidates to evaluate per batch. Smaller
+#'     \item{plan}{Future plan: "multisession", "multicore", "callr", or
+#'       "sequential".  A sequential plan runs the consistency screen as a
+#'       plain in-process loop regardless of \code{workers}: no future backend
+#'       is set up and no per-batch globals resolution is paid.  A parallel
+#'       plan with \code{workers = 1} is honoured as parallel (one-worker
+#'       batched backend); pass \code{plan = "sequential"} to run the plain
+#'       loop.  Results under \code{consistency_method = "resample"} are
+#'       identical on every path (the closed form consumes no RNG).  Under
+#'       \code{"split"} the plain loop draws from a single RNG stream seeded
+#'       once, whereas the batched backend gives each candidate its own
+#'       L'Ecuyer-CMRG stream indexed by global position, so a sequential
+#'       plan's \code{p.consistency} values differ from a parallel plan's
+#'       within Monte-Carlo noise at the stated \code{fs.splits}; invariance
+#'       to worker count holds within any plan.}
+#'     \item{workers}{Number of parallel workers on a parallel plan.}
+#'     \item{batch_size}{Number of candidates to evaluate per batch on a
+#'       genuinely parallel plan -- the knob for the per-batch overhead
+#'       described above.  Smaller
 #'       values provide finer-grained early stopping but may increase overhead.
@@ -662,6 +677,27 @@ subgroup.consistency <- function(df,
   use_parallel <- length(parallel_args) > 0 && !is.null(parallel_args[[1]])
 
+  # A sequential plan is the caller asking NOT to parallelise: run the plain
+  # loop (Section 8, sequential branch) instead of standing up a one-worker
+  # future backend.  The backend path batches the candidates into
+  # future_lapply() calls of min(2 * workers, n) and pays a full round of
+  # globals resolution per call -- measured at 29.4 s against 4.5 s for the
+  # identical 749-candidate screen, 86% of a bootstrap replicate
+  # (REPORT_bootstrap_profile_2026-08-30.md).  The plan alone decides: a
+  # parallel plan with workers = 1 stays on the batched path, so results are
+  # invariant to worker count WITHIN any plan (the per-candidate RNG streams of
+  # .make_candidate_rng_seeds() are indexed by global position; the plain loop
+  # seeds a single stream once).  Under consistency_method = "resample" no RNG
+  # is consumed and every path gives identical results.  No warning: this is
+  # the code doing what it was asked, once per replicate in a bootstrap.
+  if (use_parallel && identical(parallel_args$plan, "sequential")) {
+    if (details) {
+      cat("Parallel config: sequential plan -- running the consistency ",
+          "screen as a plain loop\n", sep = "")
+    }
+    use_parallel <- FALSE
+  }
+
   if (use_parallel) {
     required_parallel <- c("plan", "workers")
```

`valid_plans`, the required-elements / plan / workers validations, the batch-size block and the bootstrap caller are untouched. `devtools::document()` regenerated `man/subgroup.consistency.Rd` only. `DESCRIPTION` → 0.3.1; `NEWS.md` entry rewritten per §7.

---

## 5. The identity gates (§5 — GATE passed; `parallel_dispatch_compare_narrow_2026-08-30.log`)

Post-change capture `parallel_dispatch_post_narrow_2026-08-30.rds` (C1–C6) against the committed pre-change reference (C1–C5) and §3 (C6); `identical()` on `sg.harm`, membership, `treat.recommend`, `hr.subgroups`, the counters, the full result table, the top row:

| # | identical on all seven | wall pre → post | note |
|---|---|---|---|
| **C1** | **yes** | 31.08 → **5.38 s** (5.8×) | sequential plan, batch size 2 → plain loop |
| **C2** | **yes** — `early_stop_candidate` 1 = 1, `n_candidates_evaluated` 1 = 1 | 4.44 → 4.29 s | batch size 1 halts where the plain loop halts |
| C3 | selection identical (same 64 members, same `sg.harm`, same `hr.subgroups`); `Pcons` table not identical | 87.32 → 57.81 s | below |
| **C4** | **yes** | 8.87 → **7.41 s** (1.20×) | survival |
| **C5** | **yes** (bit-identical) | 26.04 → 30.27 s | `multisession, 2` — untouched |
| **C6** | **yes** (bit-identical) | 30.80 → **30.65 s** | `multisession, 1` — **stays on the batched path** (batched order, not plain-loop order) |

**C3 characterisation** (`set.seed(20260830)` before each call, `split`, 50 splits): 348 vs 346 passing rows, 296 in both; per-candidate `Pcons` difference mean +0.0005, **sd 0.028**, max |Δ| 0.10, against a binomial sd of 0.042 at 50 splits — Monte-Carlo scale, selected subgroup unchanged: benign, and **numerically identical to the first task's C3** (same figures to the digit).

**Against the first task's post-change capture** (`parallel_dispatch_post_2026-08-30.rds`, in the tree): C1, C2, C3, C4, C5 post-narrow are **identical on all seven components** to the first edit's post-change values, **C3 included** — the two edits route these five configurations identically; this edit did nothing beyond dropping the `workers` clause.

---

## 6. What was retained (§6)

| configuration | pre | first edit | **narrow edit** |
|---|---:|---:|---:|
| C1 continuous, 749 candidates screened | 31.08 s | 5.51 s | **5.38 s** |
| C2 early stop at candidate 1 | 4.44 s | 4.28 s | 4.29 s |
| C3 split, 50 splits | 87.32 s | 58.34 s | 57.81 s |
| C4 survival gbsg | 8.87 s | 7.49 s | **7.41 s** |
| C5 multisession 2 | 26.04 s | 28.18 s | 30.27 s (noise) |
| C6 multisession 1 | 30.80 s | — | **30.65 s** (unchanged, as required) |

**Bootstrap check** (`parallel_dispatch_boot_narrow_2026-08-30.log`): `forestsearch_bootstrap_dofuture()` at **`nb_boots = 20`** (as in the profile task and the first task), `plan = "sequential"`, on the C1 fit: wall **98.9 s, 4.95 s per replicate** (the bootstrap's own `tmins_search`: mean 4.79 s, median 4.88 s; 20 / 20 identified) — against the profile's 32.2 s and the first edit's 4.97 s. **B = 1000 projects to 4 946 s = 82 min sequential** (490–537 min before). The inner plan reaches the short-circuit; the whole saving is retained.

---

## 7. Close-out (§7)

- `devtools::document()`: `man/subgroup.consistency.Rd` only.
- `devtools::test()`: **`[ FAIL 0 | WARN 31 | SKIP 3 | PASS 4844 ]`** against the baseline `[ FAIL 0 | WARN 31 | SKIP 3 | PASS 4844 ]` (`parallel_dispatch_tests_post_narrow_2026-08-30.log`; baseline `parallel_dispatch_tests_pre_2026-08-30.log`). **`test-search-reproducibility.R:86` passes unedited, 6 / 6 expectations.** GATE passed.
- Warning parity: **31 = 31, and no test moved in any count** — the per-test comparison of all 313 tests (failed / warning / skipped / passed, pre vs post) shows zero differences. As the first task found, the `"parallel_args missing required elements. Using sequential."` warning is not among the 31 (every test that passes `parallel_args` passes `workers`), so the decrease the task allowed for had nothing to act on. GATE passed.
- `R CMD check` (`devtools::check(document = FALSE, args = "--no-manual")`, pandoc on PATH): **`Status: OK` — 0 errors | 0 warnings | 0 notes** (`parallel_dispatch_check_narrow_2026-08-30.log`). GATE passed.
- `DESCRIPTION` → **0.3.1**; `NEWS.md` entry rewritten (sequential plan → plain loop regardless of `workers`; the measured effect from both runs; `resample` results unchanged; `split` under a sequential plan on a single stream — Monte-Carlo-scale differences from 0.3.0 only for callers who passed a sequential plan *with* `workers`; invariance to worker count unchanged within any plan; a parallel plan with `workers = 1` still parallel; `batch_size` for genuinely parallel runs).
- Commit by explicit path; **no push**; `devtools::install(".", upgrade = FALSE)` — **installed 0.3.1**, confirmed by `packageVersion()`.

Commits: `911d73ba` task doc; the `R/`, `man/`, `DESCRIPTION`, `NEWS.md` edit with the capture script, its outputs and logs and this report in the next commit; its hash in the follow-up line.

---

## 8. Verdict (ten lines)

1. Provenance gates passed (both prior closes in the log; the patch applies cleanly); the task doc arrived as `cc_task_parallel_dispatch_narrow_20260830.md` (+ a `(1)` copy differing only in patch-location wording).
2. §2 gates passed: the test's arms are both `multisession` differing only in `workers` (1 / 2) under `split` at 40 splits; both callers pass the string `"sequential"`.
3. Seeding, from source: the plain loop calls `set.seed(seed)` once (L471); the batched path gives each candidate a global-index L'Ecuyer stream (`.make_candidate_rng_seeds`, L944–969) independent of batch layout — worker-count invariance holds within a plan by construction, and no seeding was changed.
4. The edit: the first task's patch with the condition reduced to `identical(parallel_args$plan, "sequential")`; roxygen and NEWS rewritten; nothing else.
5. **§5 gates all passed:** C1, C2, C4, C5, C6 `identical()`; C6 stays on the batched path (30.8 → 30.7 s); C1–C5 identical, C3 included, to the first task's post-change capture.
6. C3: same selection, `Pcons` sd 0.028 vs binomial 0.042 — benign, the figures of the first task reproduced exactly.
7. Retained: 31.1 → 5.4 s per continuous call, 8.9 → 7.4 s on gbsg, bootstrap 32.2 → **4.95 s** per replicate, B = 1000 → **82 min**.
8. `devtools::test()` `[ FAIL 0 | WARN 31 | SKIP 3 | PASS 4844 ]` = baseline; the worker-count invariance test passes unedited; per-test parity exact across all 313 tests.
9. `R CMD check` `Status: OK`, 0 / 0 / 0.
10. Installed 0.3.1; no test edited, no seeding changed, no push.
