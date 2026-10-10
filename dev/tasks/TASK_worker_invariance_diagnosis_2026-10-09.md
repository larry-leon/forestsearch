# TASK — Diagnose the worker-count invariance failure (read-only)

- Date: 2026-10-09
- Repo: `larry-leon/forestsearch`, branch `feature/glm-extension`, at the HEAD left by the Stage 3 task (`TASK_cran_release_wb_followup_2026-10-09.md`). Run in a fresh CC session after Stage 3 has finished.
- Why: `tests/testthat/test-search-reproducibility.R`, test "search + consistency results are invariant to worker count", fails under `devtools::test()` (`NOT_CRAN=true`): one run returned `sg.harm = NULL`, the other `c('!{x1 <= 0.6}', '{b1}')`; `M.2` differs (Stage 3 halt note, identical failure at `50b1b9a5`). The package documents that consistency results are invariant to worker count; that guarantee passed at 0.3.1 (2026-08-30). Larry requires this understood before the CRAN submission.
- Type: **diagnosis only.** No edits to `R/`, `tests/`, `man/`, `DESCRIPTION`, `NAMESPACE`. No fixes. Anything that would change behaviour is reported as a proposal, never applied.
- R/ callout: none.
- Authorizations (kickoff): installing forestsearch builds of specific commits into a **scratch library only** (never the user library); temporary git worktrees under the scratchpad (removed at the end); running the single test file and small diagnostic scripts. No `R CMD check`, no full suite, no uploads, no `git fetch`/`pull`/`push`.
- Hard cap: 90 min wall-clock for the whole task; each step has its own cap. A timeout is a STOP with the report written to the point reached.
- Gates are stop-on-failure.

## First action

1. Copy this file to `dev/tasks/TASK_worker_invariance_diagnosis_2026-10-09.md` and commit it alone.

## Preconditions

- P1. `git status --porcelain` shows no tracked-file modifications. Record HEAD SHA.
- P2. Record the user library's installed forestsearch: `packageVersion("forestsearch")`, `packageDescription("forestsearch")$RemoteSha` (or `Packaged`/`Built` if not a remote install), and its library path. Do not modify it.

## Steps

### Step 1 — What the test does (≤ 10 min, read-only)
Record verbatim with file:line:
- The full failing `test_that` block and any helpers it calls.
- Consistency mode used (split / resample / closed-form), `n.splits` or equivalent, `pconsistency.threshold`, `pconsistency.digits`, seed arguments.
- The two plans compared (e.g. `sequential` vs `multisession`), `workers`, and how the plan is set and restored.
- The assertion(s) that fail and what `M.2` is.

### Step 2 — Harness hypothesis: stale installed package in workers (≤ 15 min)
- `devtools::test()` loads source via `load_all()` in the main session, but `multisession` workers load the forestsearch **installed** in the library path. If that build differs from HEAD, the two routes run different code.
- Compare P2's installed build with HEAD (version and SHA). Record.
- Install HEAD into `<scratch>/rlib` (`R CMD INSTALL --library=<scratch>/rlib .`). Run only this test file in a fresh R process with `R_LIBS=<scratch>/rlib` set in the environment (so workers inherit it) and `NOT_CRAN=true`, using `library(forestsearch)` from that library (not `load_all`): `testthat::test_file("tests/testthat/test-search-reproducibility.R", package = "forestsearch", load_package = "installed")`. Confirm inside the test session and inside one worker (`future::value(future::future(packageDescription("forestsearch")$Packaged))`) that both load the scratch build.
- Run it 3 times. Record pass/fail each time.
- **Gate 2:** if all 3 pass → classification **H0 (harness artifact)**; skip to the Deliverable (Steps 3–5 not needed). If any fail → continue.

### Step 3 — Random-number routing (≤ 15 min, read-only)
From source at HEAD, trace how random numbers are generated for the consistency splits on each route, with file:line:
- The sequential route (including the 0.3.1 sequential-plan short-circuit / plain loop) and the parallel route (`future`/`foreach`/`doFuture` with `seed = TRUE`, `%dorng%`, `set.seed` calls, `future.seed`, L'Ecuyer streams).
- Whether streams are assigned per task (candidate/split) or per worker / by execution order; whether both routes use the same generator; whether any filtering, ordering or tie-breaking happens before draws on one route but not the other.

### Step 4 — Bisect (≤ 45 min)
- Identify the 0.3.1 commit (the 2026-08-30 "narrow exit" parallel-dispatch commit; search `git log` for it and confirm `DESCRIPTION` reads 0.3.1 there). Record SHA.
- Confirm the test passes there using the Step 2 procedure (install that commit from a worktree into a separate scratch library). If it does not pass at 0.3.1 either, record that and skip bisection.
- Otherwise `git bisect` between that commit (good) and HEAD (bad) with a script that, per commit: installs the worktree into a per-commit scratch library and runs the test file as in Step 2 (once; NOT_CRAN=true; result = pass/fail of this one `test_that` block only). Skip commits that do not install (`git bisect skip`).
- Record the first bad commit: SHA, subject, files changed in `R/`.
- Remove worktrees and scratch libraries afterwards.

### Step 5 — The flipping candidate (≤ 15 min)
At HEAD (scratch library), reproduce the test's data and settings in a script and, on each route, record for the candidate `!{x1 <= 0.6} & {b1}`: its effect estimate, the per-split consistency indicators or proportion, the rounded proportion, and the admission decision; also the top 5 candidates of each route's `M.2`. Identify the first quantity that differs between routes.

## Deliverable

`dev/reports/REPORT_worker_invariance_diagnosis_2026-10-09.md` (bullets), sections:
1. Classification — one of: **H0 harness artifact** (workers load a stale installed build); **H1 random-number routing** (routes draw differently by design or by defect); **H2 code-path defect** (non-random difference between routes); **H3 knife-edge** (same random draws, tiny numeric difference flips a threshold) — with the evidence for it.
2. Repo and installed-package state (P1–P2).
3. Steps 1–5 findings with file:line.
4. If H1/H2/H3: the minimal change that would restore invariance, as a **proposal** (classified: changes behaviour / changes method), with which committed results it could move. Not applied.

Commit the report; copy it to `~/Downloads`.

## Post-conditions

- PC1. `git diff --name-only <task-doc-commit>^ HEAD` lists only this task doc and the report.
- PC2. `git status --porcelain` clean; no worktrees left (`git worktree list` shows only the main tree); no scratch libraries inside the repo.
- PC3. The user library's forestsearch matches P2 (same version and SHA/Packaged).
