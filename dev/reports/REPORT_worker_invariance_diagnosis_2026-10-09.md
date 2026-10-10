# REPORT — Worker-count invariance failure: diagnosis (2026-10-09)

- Task: `dev/tasks/TASK_worker_invariance_diagnosis_2026-10-09.md` (committed alone at `f1d86f50`).
- Type: diagnosis only. No edits to `R/`, `tests/`, `man/`, `DESCRIPTION`, `NAMESPACE`.
- Outcome: **Gate 2 passed (3/3) → classification H0. Steps 3–5 not run, as the task specifies.**

## 1. Classification — **H0: harness artifact**

- The failure is caused by **`multisession` workers loading a stale installed build**. Under `devtools::test()` the main session runs HEAD source through `load_all()`, but each `multisession` worker runs `loadNamespace("forestsearch")` from the user library. That copy is **0.2.2** (Packaged 2026-08-19), not HEAD's 0.4.0.
- Evidence:
  - **Installed-HEAD route (Step 2): PASS 3/3.** Main session and a worker both resolved `forestsearch` 0.4.0 from the scratch library.
  - **Positive control: FAIL, reproducing the Stage 3 failure exactly.** This was `load_all()` on the same HEAD tree, with workers on the default library path; the worker reported `/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/forestsearch`, version 0.2.2. The failures were:
    - `test-search-reproducibility.R:86`: the `out_sg$result` tables differ. One run has 5 rows (top row Pcons 0.97, `!{x1 <= 0.6}` & `{b1}`); the other is empty.
    - `test-search-reproducibility.R:90`: `sg.harm` is `c('!{x1 <= 0.6}', '{b1}')` vs `NULL`.
  - Only the code the workers load differs between the two routes. Fixture, seed and HEAD source are identical. The worker-invariance property itself holds at HEAD.
- Why only one of the two runs breaks (an inference from `future`'s documented behaviour; I didn't instrument it):
  - `future::multisession` with `workers = 1` falls back to sequential, in-process evaluation, so the `workers = 1L` run executes the `load_all()` HEAD code.
  - The `workers = 2L` run sends `eval_fun` to background R processes, which run 0.2.2.
  - The 0.2.2 result yields an empty consistency table. The `suppressWarnings({ future.apply::future_lapply(...) })` wrapper at `R/subgroup_consistency_main.R:986` hides `future`'s package-version-mismatch warning, which would otherwise have pointed straight at the cause.

## 2. Repo and installed-package state

- P1: tree clean at start. Pre-task HEAD `1ddd8a37`; task doc commit `f1d86f50`; branch `feature/glm-extension`; `DESCRIPTION` Version 0.4.0.
- P2: the user library's forestsearch is version **0.2.2**:
  - `RemoteSha`: none (not a remote install).
  - Packaged: `2026-08-19 06:39:06 UTC; larryleon`. Built: `R 4.5.2; ; 2026-08-19 06:39:07 UTC; unix`.
  - Path: `/Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/forestsearch`.
- PC3: re-checked after cleanup; unchanged (0.2.2, same Packaged stamp, same path).

## 3. Findings

### Step 1 — What the test does (`tests/testthat/test-search-reproducibility.R`)

- The block is `test_that("search + consistency results are invariant to worker count", ...)`, lines 65–91, with `skip_on_cran()` at line 66.
- Helpers:
  - `.repro_df()`, lines 26–35: n = 360, `set.seed(42)`, continuous outcome with true harm region `x1 > 0 & b1 == 1`.
  - `.repro_run(workers)`, lines 37–55.
  - `.repro_flat()`, lines 59–63: drops data.table attributes before `identical()`.
- Settings (lines 38–52):
  - Outcome: `outcome_type = "continuous"`, `effect_measure = "MD"`.
  - Thresholds: `effect.threshold = 2`, `consistency.threshold = 2`, `pconsistency.threshold = 0.05`.
  - Consistency: `consistency_method = "split"` (literal random splits, chosen on purpose; see the header comment at lines 19–24), `fs.splits = 40L`, `stop_threshold = NULL` (every candidate evaluated), `sg_focus = "maxeffCons"`.
  - Size limits: `n.min = 30`, `d0.min = 5`, `d1.min = 5`, `maxk = 2L`.
  - Seed: `seedit = 20260827`. `pconsistency.digits` is not set, so the default applies.
- Plans compared: both are `plan = "multisession"`, with `workers = 1L` vs `workers = 2L` (lines 50, 68–69). After each run, line 53 resets the plan with `future::plan(future::sequential)`. Inside the call, `R/subgroup_consistency_main.R:879-880` saves the plan and restores it with `on.exit`.
- Assertions:
  - Lines 75–78 check that the fixture is in the multi-batch regime.
  - Line 81: the search table is identical across runs.
  - Line 86: the per-candidate consistency table `out_sg$result` is identical across runs.
  - Line 90: `sg.harm` is identical across runs.
- `M.2` is the second-factor column of the `out_sg$result` table that line 86 compares.

### Step 2 — Harness hypothesis

- Installed build vs HEAD: **0.2.2 (2026-08-19) vs 0.4.0 at `f1d86f50`.** The installed build predates the 0.3.x per-candidate RNG-stream work and the GLM consistency changes.
- Procedure:
  - Installed HEAD from a detached scratch worktree into `<scratchpad>/rlib` with `R CMD INSTALL --no-docs`.
  - Ran 3 times in fresh `Rscript` processes with `R_LIBS=<scratch>/rlib`, `NOT_CRAN=true`, `library(forestsearch)`, and `testthat::test_file(..., package = "forestsearch", load_package = "installed")`.
- Load check, every run:
  - Main session: `<scratch>/rlib/forestsearch`, 0.4.0.
  - `future::value(future::future(...))` worker: `<scratch>/rlib/forestsearch`, 0.4.0.
  - The `Packaged` field is absent for a source-directory install, so I checked the library path and version instead.
- Results: run 1 PASS (6 expectations, 0 failed); run 2 PASS; run 3 PASS.
- Positive control (the `load_all()` route described in §1): FAIL, matching the Stage 3 failure.
- **Gate 2: 3/3 pass → H0.**

### Steps 3–5

- Not run. Gate 2 classified the failure as H0, and the task says to skip these steps in that case.
- For reference, the invariance mechanism at HEAD:
  - `R/subgroup_consistency_main.R:40-66`: `.make_candidate_rng_seeds()` creates one L'Ecuyer-CMRG stream per global candidate index.
  - `R/subgroup_consistency_main.R:967, 987-994`: each `future_lapply` call gets `future.seed = candidate_seeds[batch_indices]`.
  - `R/subgroup_consistency_main.R:713`: a sequential plan short-circuits to the plain loop.
- The 0.3.1 bisection anchor wasn't needed and was not identified.

## 4. Proposals (none applied)

H0 needs no change to package behaviour or method. All four proposals below are harness or test changes; none moves any committed result.

1. **(Recommended, no code change)** Before any `devtools::test()` run that includes multisession tests, reinstall HEAD into the user library (`devtools::install(dependencies = FALSE)`), or run the test file against an installed build as in Step 2. This brings the user library from 0.2.2 to 0.4.0, so it needs Larry's go-ahead.
2. **(Test hygiene; changes no package behaviour)** Add a guard to `test-search-reproducibility.R`, and to any other test that uses `multisession`: skip with an explicit message when `pkgload::is_dev_package("forestsearch")` is TRUE and the installed version differs from the source `DESCRIPTION` version. Without it, a stale library produces a misleading invariance failure.
3. **(Diagnostics; changes warning output only)** The `suppressWarnings()` around `future_lapply` (`R/subgroup_consistency_main.R:986`) hides `future`'s package-version-mismatch warning, which would have named the cause at once. Narrowing the suppression to the specific warnings it targets would let that one through. This changes user-visible warnings, not results.
4. Keep the CRAN path in mind:
   - The test is `skip_on_cran()`.
   - Under `R CMD check` the package is installed before tests run, so workers load the same build.
   - This failure therefore cannot occur on win-builder or CRAN, and doesn't block submission.

## Post-conditions

- PC1: this task changes only the task doc (`f1d86f50`) and this report.
- PC2: the scratch worktree has been removed (`git worktree list` shows only the main tree). The scratch library was in the session scratchpad and has been deleted. Nothing was added inside the repo apart from the two documents.
- PC3: the user library is unchanged (0.2.2, Packaged `2026-08-19 06:39:06 UTC`).
- Wall-clock: about 3 min (17:07–17:11 PDT).
