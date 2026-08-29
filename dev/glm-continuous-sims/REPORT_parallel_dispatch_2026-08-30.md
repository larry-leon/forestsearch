# REPORT — the consistency screen's parallel dispatch: `plan = "sequential"` / `workers <= 1` as sequential

**Task:** `dev/tasks/cc_task_parallel_dispatch_2026-08-30.md` (commit `1e45219f`)
**Repo / branch:** `~/Documents/GitHub/forestsearch` (pop-os), `feature/glm-extension`
**Follows:** `REPORT_bootstrap_profile_2026-08-30.md`.

---

## ⛔ GATE FAILED at §6 — the edit is NOT committed; the tree is clean; `R/` is unchanged; the installed package is 0.3.0

**What failed.** With the §3 edit in place exactly as specified — sequential plan **or** `workers <= 1` under any plan → plain loop — `devtools::test()` returns **`[ FAIL 1 | WARN 31 | SKIP 3 | PASS 4843 ]`** against the baseline **`[ FAIL 0 | WARN 31 | SKIP 3 | PASS 4844 ]`**. The one failure is `tests/testthat/test-search-reproducibility.R:86`, "search + consistency results are invariant to worker count". That test pins a documented package invariant (its header, L1–24): under `consistency_method = "split"` the per-candidate `Pcons` table is **identical between `multisession` runs at `workers = 1` and `workers = 2`**, because the batched path hands `future_lapply` one L'Ecuyer-CMRG stream per candidate indexed by *global* position (`.make_candidate_rng_seeds`, `subgroup_consistency_main.R` L944–970). The edit routes `workers = 1` — plan `multisession` included, as §3 asked — to the plain loop, which seeds once with `set.seed(seed)` (L471) and draws the splits from a single stream, so workers 1 and 2 now draw different splits. The failing rows (`Pcons` 0.90 vs 0.97, 0.80 vs 0.73, 0.50 vs 0.57 … at `fs.splits = 40`) are Monte-Carlo-scale differences — the same phenomenon §4 characterised as benign for C3 — but the test is the executable form of a contract the package states, and the task's §6 forbids adjusting a test to pass without a stated reason, while §3 forbids narrowing the condition and the editable scope excludes the plain loop's seeding. Every route out of the failure is a decision the task did not authorise:

1. **Narrow the short-circuit to `identical(plan, "sequential")`** and drop the `workers <= 1` clause. The bootstrap (`plan = "sequential", workers = 1L`) and the driver still reach the plain loop; `multisession, workers = 1` keeps the batched path and the invariant holds. This is the smallest change and the one I would recommend if asked — but it deviates from §3's specification.
2. **Re-pin the test** to compare `workers = 2` against `workers = 3` (or state the invariant as "among genuinely parallel worker counts"). Adjusts a test to pass; the invariant narrows.
3. **Seed the plain loop with the same per-candidate streams** the batched path uses, so both branches are RNG-identical and C3 becomes `identical()` too. The principled fix, and it edits Section 8 — outside the task's editable scope.

**What is committed.** The task document (`1e45219f`), the capture script, the pre/post/compare/bootstrap outputs and logs, the two `devtools::test()` logs, this report, and the edit itself as a patch file (`parallel_dispatch_edit_2026-08-30.patch`, 154 lines — `R/subgroup_consistency_main.R`, `man/subgroup.consistency.Rd`, `DESCRIPTION` 0.3.1, `NEWS.md`), so it can be applied under whichever of the three decisions is taken. Nothing under `R/`, `man/`, `DESCRIPTION` or `NEWS.md` differs from `1e45219f`. No install, no push.

**What was established before the gate closed** (all sections below): the source verification passed; the §4 identity gates **all passed** (C1, C2, C4, C5 `identical()`; C3 same selection, MC-scale `Pcons` differences); the timings and the bootstrap check were measured (29.4 → 5.5 s per call; 32.2 → 4.97 s per replicate; B = 1000: 490–537 → 83 min); warning parity is exact (31 = 31, no per-test movement). The edit does what the task said it would; the one thing it also does is break the worker-count invariance test on the `split` path at `workers = 1`.

---

## 1. Provenance (§1 — GATE passed)

```
pop-os · /home/larryleon/Documents/GitHub/forestsearch · feature/glm-extension
HEAD at start 80ffa2da · git status --porcelain: empty
git log -4: 80ffa2da d77702ce b8daf89c 9d53f465
packageVersion forestsearch 0.3.0
```

`~/Downloads/cc_task_parallel_dispatch_2026-08-30.md` was absent under that name; the stem arrives hyphen-stripped as **`cc_task_parallel_dispatch_20260830.md`**, alongside chat's companion **`parallel_dispatch_sites_20260830.md`** (its three source quotations). → `dev/tasks/cc_task_parallel_dispatch_2026-08-30.md`, committed alone as **`1e45219f`**.

### 1.1 Source verification against the working tree (GATE passed — the premise has not moved)

**Section 6 validation block**, `R/subgroup_consistency_main.R` **L663–683** — identical to chat's quotation, same line numbers:

```r
  use_parallel <- length(parallel_args) > 0 && !is.null(parallel_args[[1]])

  if (use_parallel) {
    required_parallel <- c("plan", "workers")
    if (!all(required_parallel %in% names(parallel_args))) {
      warning("parallel_args missing required elements. Using sequential.")
      use_parallel <- FALSE
    }

    valid_plans <- c("multisession", "multicore", "callr", "sequential")
    if (use_parallel && !parallel_args$plan %in% valid_plans) {
      warning("Invalid parallel plan. Using sequential.")
      use_parallel <- FALSE
    }

    if (use_parallel && (!is.numeric(parallel_args$workers) ||
                         parallel_args$workers < 1)) {
      warning("Invalid workers value. Using sequential.")
      use_parallel <- FALSE
    }
  }
```

**Batch-size block**, **L832–858** — identical to chat's quotation (chat: L831–858):

```r
    n_workers <- parallel_args$workers

    if (!is.null(parallel_args$batch_size)) {
      batch_size_parallel <- min(as.integer(parallel_args$batch_size), n_candidates)
    } else if (!is.null(stop_threshold)) {
      if (sg_focus %in% c("maxeffCons")) {
        batch_size_parallel <- 1L
      } else {
        batch_size_parallel <- min(n_workers, max(1L, n_candidates %/% 4), n_candidates)
      }
    } else {
      batch_size_parallel <- min(n_workers * 2L, n_candidates)      # L857
    }
```

**The caller**, `R/bootstrap_dofuture_main.R`: L279 `parallel_args <- resolve_bootstrap_parallel_args(parallel_args, args_forestsearch_call)`; L366 `setup_parallel_SGcons(parallel_args)`; L443 `cat("  - parallel_args: sequential (prevents nested parallelism)\n")` — all three at chat's line numbers. The inner replicate's arguments are set in `R/bootstrap_analysis_dofuture.R` L573–575: `args_FS_boot$parallel_args$plan <- "sequential"; $workers <- 1L; $show_message <- FALSE`.

---

## 2. The pre-change reference (§2)

Fixtures as in the profile task: **continuous** — the MD40 DGM rebuilt and gated on the payload's `truth`, one trial `simulate_from_glm_dgm(dgm, n = 500, seed = 8316952)`, the drivers' `forestsearch()` arguments (13 confounders, `maxeffCons`, `resample`, `fs.splits = 400`, `vi.grf.min = NULL`, `mr_inference = FALSE`); **survival** — `survival::gbsg` (686 rows, 299 events) with the effMaxSG application's arguments. Script `parallel_dispatch_capture_2026-08-30.R`; pre-change values in `parallel_dispatch_pre_2026-08-30.rds`, captured at HEAD `1e45219f` with `R/` unedited.

| # | fixture | `parallel_args` | consistency | early stop | selected (membership size) | `hr.subgroups` | evaluated / total, passed | early stop at | top row: effect, Pcons, N | wall |
|---|---|---|---|---|---|---:|---|---|---|---:|
| C1 | cont | `sequential, workers 1` | resample | off | `!{wtkg <= 84} & !{cd40 <= 320}` (64) | 852 | 749 / 749, 319 | — | 101.05, 1.00, 64 | **31.08 s** |
| C2 | cont | `sequential, workers 1` | resample | `stop_threshold = 0.95` | same (64) | 852 | **1 / 749, 1** | **candidate 1** | 101.05, 1.00, 64 | 4.44 s |
| C3 | cont | `sequential, workers 1` | **split, `fs.splits = 50`**, `set.seed(20260830)` before the call | off | same (64) | 852 | 749 / 749, **348** | — | 101.05, 1.00, 64 | 87.32 s |
| C4 | surv | `sequential, workers 1` | resample | off | `{er <= 0} & {size <= 35}` (61) | 121 | 120 / 120, 9 | — | HR 2.537, 0.99, 61 | **8.87 s** |
| C5 | cont | `multisession, workers 2` | resample | off | same as C1 (64) | 852 | 749 / 749, 319 | — | 101.05, 1.00, 64 | 26.04 s |

Every run was warning-free. Captured per configuration: `sg.harm`, `sort(which(sg.harm.id == 1))`, `df.est$treat.recommend`, the full `hr.subgroups` table, the `grp.consistency` counters, the full `out_sg$result` table and its top row, wall-clock.

---

## 3. The edit (§3) — written, verified, then reverted at the §6 gate; preserved as `parallel_dispatch_edit_2026-08-30.patch`

`R/subgroup_consistency_main.R`, inserted between L663 and the `required_parallel` check:

```r
  if (use_parallel) {
    w <- parallel_args$workers
    single_worker <- !is.null(w) && is.numeric(w) && length(w) == 1L && !is.na(w) && w <= 1
    if (identical(parallel_args$plan, "sequential") || single_worker) {
      if (details) {
        cat("Parallel config: sequential plan or a single worker -- ",
            "running the consistency screen as a plain loop\n", sep = "")
      }
      use_parallel <- FALSE
    }
  }
```

with a comment block citing the measurement, no `warning()`, and the `@param parallel_args` roxygen rewritten (plan / workers / batch_size items: a sequential plan or a single worker runs the plain loop; performance property, not results property; identical under `resample`, RNG re-laid-out under `split`; `batch_size` is the knob for genuinely parallel runs). `valid_plans`, the batch-size logic and the bootstrap caller untouched. `devtools::document()` regenerated `man/subgroup.consistency.Rd` only. `DESCRIPTION` → 0.3.1 and the `NEWS.md` entry were written (both in the patch) and reverted with the rest.

---

## 4. The identity gates (§4) — **all PASS** (`parallel_dispatch_compare_2026-08-30.log`)

Post-change capture at the same HEAD with the edit in the working tree (`parallel_dispatch_post_2026-08-30.rds`). `identical()` on `sg.harm`, membership, `treat.recommend`, `hr.subgroups`, counters, the full result table, the top row:

| # | identical on all seven | wall pre → post | note |
|---|---|---|---|
| **C1** | **yes** | 31.08 → **5.51 s** (5.6×) | the main path, batch size 2 → plain loop |
| **C2** | **yes** — `early_stop_candidate` 1 = 1, `n_candidates_evaluated` 1 = 1 | 4.44 → 4.28 s | batch size 1 halted where the plain loop halts |
| C3 | selection **identical** (same 64 members, same `sg.harm`, same `hr.subgroups`); `Pcons` table **not** identical | 87.32 → 58.34 s | characterised below |
| **C4** | **yes** | 8.87 → **7.49 s** (1.18×) | survival |
| **C5** | **yes** (bit-identical) | 26.04 → 28.18 s | the genuinely parallel path untouched |

**C3 characterisation.** 348 vs 346 candidates passing the screen; 296 in both tables (52 pre-only, 50 post-only — candidates near the 0.90 threshold that crossed it one way or the other). Per-candidate `Pcons` difference over the 296: mean +0.0005, **sd 0.028**, max |Δ| 0.10, quartiles −0.02 / 0.00 / +0.02, against a binomial sd of **0.042** for a rate near 0.9 at 50 splits. Within Monte-Carlo noise for the stated `fs.splits`, selected subgroup unchanged: **benign and expected** under §4's rule. The mechanism is the one the failing test formalises: the plain loop seeds once (`set.seed(seed)`, L471) where the batched path assigns per-candidate streams, so the same seed lays the splits out differently.

---

## 5. What was bought (§5)

| configuration | pre | post | factor |
|---|---:|---:|---:|
| C1 continuous, 749 candidates screened (the profile's 29.4 s case) | 31.08 s | **5.51 s** | 5.6× |
| C4 survival gbsg, 120 candidates screened — **pinned: `parallel_args = list(plan = "sequential", workers = 1L)` on both sides** | 8.87 s | **7.49 s** | 1.18× (1.4 s, 16 % — the profile's 18 % "future overhead" bucket) |
| C2 (early stop at candidate 1) | 4.44 s | 4.28 s | — |
| C3 (split, 50 splits) | 87.32 s | 58.34 s | 1.5× |
| C5 (multisession 2) | 26.04 s | 28.18 s | unchanged (noise) |

**Bootstrap check** (`parallel_dispatch_boot_2026-08-30.log`): `forestsearch_bootstrap_dofuture()` at `nb_boots = 20`, `plan = "sequential"`, on the C1 fit, with the edit in place: wall **99.3 s, 4.97 s per replicate** (the bootstrap's own `tmins_search`: mean 4.82 s, median 4.88 s; 20 / 20 identified) against the profile's **32.2 s**. **B = 1000 projects to 4 967 s = 83 min sequential**, against 490–537 min.

---

## 6. Close-out (§6) — where the gate closed

- `devtools::document()`: `man/subgroup.consistency.Rd` only.
- `devtools::test()` **baseline** (unchanged tree, `1e45219f`): `[ FAIL 0 | WARN 31 | SKIP 3 | PASS 4844 ]` (`parallel_dispatch_tests_pre_2026-08-30.log`). **Post-change:** `[ FAIL 1 | WARN 31 | SKIP 3 | PASS 4843 ]` (`parallel_dispatch_tests_post_2026-08-30.log`).
- **Warning-parity accounting, test by test:** per-test warning counts were compared across all 313 tests — **no test's warning count moved** (31 = 31). The `"parallel_args missing required elements. Using sequential."` warning was never among the 31: every test that passes `parallel_args` passes `workers` too (`test-glm-pipeline-integration.R`, `test-cross-outcome-parity.R`, `helper-synthetic-dgm.R`, the CV/MR helpers), so no test reached the plain loop by that route before, and none loses the warning now. Parity holds; the decrease the task anticipated did not occur because it had nothing to act on in the suite.
- **The failure:** one test, one expectation — the `Pcons`-table `expect_identical()` at `test-search-reproducibility.R:86` — described at the top. The other five expectations of that test pass (candidate count, `hr.subgroups` identical, `sg.harm` identical). No other test changed in any count.
- `R CMD check` (`devtools::check(document = FALSE, args = "--no-manual")`, pandoc on PATH), run on the edited tree before the revert: **`Status: 1 ERROR`** — the single ERROR is `checking tests`, the same `test-search-reproducibility.R:86` failure; no WARNINGs, no NOTEs (`parallel_dispatch_check_2026-08-30.log`). Not 0/0/0.
- `DESCRIPTION` 0.3.1 and `NEWS.md`: written, in the patch, reverted. **No `devtools::install()`; installed remains 0.3.0.**
- Tree after the revert: `git status --porcelain` shows only the new `dev/glm-continuous-sims/` files; `R/`, `man/`, `DESCRIPTION`, `NEWS.md` at `1e45219f`.

Commits: `1e45219f` task doc; the capture script, its outputs and logs, the patch and this report in the next commit; its hash in the follow-up line.

---

## 7. Verdict (ten lines)

1. Source verification: the dispatch block, the batch-size block and the caller are exactly where and what chat quoted; premise intact.
2. The §3 edit was written as specified (sequential plan **or** `workers <= 1` → plain loop; no warning; roxygen updated) and is preserved as a patch.
3. §4 gates: **C1, C2, C4, C5 `identical()` on every captured component**; C2 halts at the same candidate; C5 bit-identical.
4. C3 (`split`, 50 splits): same selection; per-candidate `Pcons` sd 0.028 vs binomial 0.042 — benign, as §4 allows.
5. Bought: 31.1 → 5.5 s per continuous call (5.6×); 8.9 → 7.5 s on gbsg; bootstrap 32.2 → 4.97 s per replicate; B = 1000 from 490–537 min to 83 min.
6. Warning parity exact: 31 = 31, no per-test movement; the anticipated decrease had no test to act on.
7. **§6 gate failed:** `test-search-reproducibility.R:86` — worker-count invariance of `Pcons` under `split` — fails because `multisession, workers = 1` now takes the plain loop (single `set.seed` stream) while `workers = 2` keeps per-candidate streams.
8. Every exit is an unauthorised decision: narrow to `plan == "sequential"` only (deviates from §3), re-pin the test (adjusts a test), or seed the plain loop per candidate (Section 8, out of scope). Option 1 is the smallest and keeps the bootstrap's and the driver's fast path.
9. Tree left clean: edit reverted, `R/` unchanged, 0.3.0 installed, no push; evidence and patch committed.
10. The measurements stand and transfer: whichever option is chosen, the identity results of §4 for the sequential-plan configurations do not depend on the `workers <= 1` clause.
