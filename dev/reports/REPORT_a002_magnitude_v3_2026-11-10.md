# REPORT — A-002 magnitude check, v3, 2026-11-10

**Outcome: STOPPED at Step 2.** Step 2 says to choose `parallel_args` for
`forestsearch_bootstrap_dofuture()` that leave the arm's cluster plan in place, and to stop and
report the lines if every value resets the plan. Every value does: the function either replaces
the plan or stops with an error before the replicate loop. Neither arm was run (Steps 3–5 not
reached), so there is no results table and no verdict.

## 1. Header

| Item | Value |
|---|---|
| `master` SHA | `217656e2b7142ad62932d12da91698752f174a96` |
| `fs-a002-check` before this report | `9c2404710aba2e6264c0a4272660bcd679be9846` (v3 task commit); `R/` identical to `master` (0 diff lines) |
| Version (DESCRIPTION) | 0.4.1.9002 |
| R / future / doFuture / pkgload | 4.6.2 / 1.80.0 / 1.4.0 / 1.5.3 |
| B (`nb_boots`) | 500 (planned; not run) |
| Bootstrap seed | 8316951 (planned; not run) |
| Wall clock, arm A | not run |
| Wall clock, arm B | not run |
| Scratch | `/tmp/a002_2026-11-10/v3/` |

Changes from the document (`quarto/applications/gbsg/analysis_gbsg_survival_effMaxSG.qmd`), as
carried in `v3/app_v3.R` (a copy of v1's `app.R`; v1's marks `# A002:`, v3's marks `# A002v3:`):

- v1 (kept): `library(forestsearch)` → `pkgload::load_all(<tree>)`; the `forestsearch()` call's
  `parallel_args` plan `"multisession"` → `"sequential"`; the bootstrap chunk stops after its
  timing line.
- v3: the arm's tree path and arm label are command-line arguments (`Rscript --vanilla app_v3.R
  <tree> <A|B>`), loaded with `pkgload::load_all()`.
- v3: `mr_inference = FALSE` in the `forestsearch()` call (`mr_inference_args` left as is).
- v3: `NB` 1000L → 500L; v1's fit-only mode removed.
- v3: before the bootstrap, `cl <- parallel::makeCluster(100)`, `clusterEvalQ(cl,
  pkgload::load_all(<tree>, quiet = TRUE))`, `future::plan(future::cluster, workers = cl)`, cluster
  stopped on exit; then the Step 3 worker gate (`deparse(forestsearch::forestsearch)` on every
  worker, checked for `get_dfpred(df.predict, sg.harm`).
- v3: the bootstrap call gets `seed = 8316951L`; its `parallel_args` could not be chosen (§3), so
  the script stops with an error just before the call. The output path would have been `v3/<arm>.rds`.

## 2. Gate checklist

| Gate | Result | Evidence |
|---|---|---|
| Step 1: `git status --porcelain` empty on `fs-a002-check` | PASS | 0 lines |
| Step 1: task copied to `dev/tasks/` and committed | PASS | `9c240471 dev: add A-002 magnitude task v3`; `cmp` with `~/Downloads` copy identical |
| Step 2: a `parallel_args` value that leaves the cluster plan in place | **STOP** | none exists (§3) |
| Step 3: worker gate, arm A / arm B | not run | stopped at Step 2 |
| Step 4: arm A, arm B ≤ 30 min | not run | stopped at Step 2 |
| Step 5: both fits select `{er <= 0} & {size <= 35}`; per-replicate found/not-found, `M.1`, `M.2` equal | not run | stopped at Step 2 (the arm A fit in the §3 probe, `mr_inference = FALSE`, selected `{er <= 0} & {size <= 35}` in 2.4 s) |
| Step 7: `git diff --name-status master..fs-a002-check` only added files under `dev/tasks/`, `dev/reports/`; `git status --porcelain` empty | run after this report's commit | result in the reply and in `v3/VERIFY_a002_magnitude_v3_2026-11-10.txt` |

## 3. Why Step 2 stops: the lines

`forestsearch_bootstrap_dofuture()` (R/bootstrap_dofuture_main.R) handles the plan in three places.

1. The supplied value is resolved; an empty value (`list()`, `NULL`) inherits the fit's
   `parallel_args` (here `plan = "sequential"`):

   ```r
   # R/bootstrap_dofuture_main.R:279
   parallel_args <- resolve_bootstrap_parallel_args(parallel_args, args_forestsearch_call)
   # R/bootstrap_dofuture_setup_helpers.R:21-22
   if (length(parallel_args) == 0) {
     resolved_args <- as.list(forestsearch_call_args$parallel_args)
   ```

2. The plan is then set unconditionally, before the index matrix and the replicate loop:

   ```r
   # R/bootstrap_dofuture_main.R:366
   setup_parallel_SGcons(parallel_args)
   ```

   `setup_parallel_SGcons()` (R/subgroup_consistency_helpers.R:32-74) admits four plans and calls
   `future::plan()` in each; anything else stops:

   ```r
   # :41   if (is.null(plan_type)) stop("parallel_args$plan must be specified.")
   # :43   allowed_plans <- c("multisession", "multicore", "callr", "sequential")
   # :45-47 if (!plan_type %in% allowed_plans) stop("plan_type must be one of: ", ...)
   # :57   future::plan(future::multisession, workers = n_workers)
   # :60   future::plan(future::multicore, workers = n_workers)
   # :66   future::plan(future.callr::callr, workers = n_workers)
   # :69   future::plan(future::sequential)
   ```

   `workers` must be numeric (:50-54), so a cluster object passed as `workers` becomes 1 worker.

3. On exit the plan is set again, to `sequential`, because `.Last.future.plan` does not exist
   (it is declared only in R/globals.R:30 as a global for R CMD check):

   ```r
   # R/bootstrap_dofuture_main.R:349-356
   on.exit({
     options(warn = old_warn)  # Restore first
     if (exists(".Last.future.plan")) {
       future::plan(.Last.future.plan)
     } else {
       future::plan("sequential")
     }
   ```

**Live probe** (`v3/plan_probe.R`, log `v3/probe/plan_probe.log`, table `v3/plan_probe.rds`; 8 s, no
replicate run). Arm A tree loaded with `load_all()`; the document's fit with `mr_inference = FALSE`;
then for each candidate a 2-worker `future::cluster` plan was set, the real
`forestsearch_bootstrap_dofuture(fs, nb_boots = 500L, seed = 8316951L, ...)` was called, the plan
was recorded as the next call after setup (`.make_boot_index_matrix()`, :383) started, and the call
was aborted there by a trace.

| `parallel_args` | plan before | plan after setup (:366) | plan after exit | outcome |
|---|---|---|---|---|
| `list()` (inherits the fit's) | cluster (2 w) | sequential (1 w) | sequential (1 w) | reached :383 |
| `NULL` | cluster (2 w) | sequential (1 w) | sequential (1 w) | reached :383 |
| `plan = "sequential"` | cluster (2 w) | sequential (1 w) | sequential (1 w) | reached :383 |
| `plan = "multisession"`, 100 | cluster (2 w) | multisession (100 w) | sequential (1 w) | reached :383 |
| `plan = "multicore"`, 100 | cluster (2 w) | multicore (100 w) | sequential (1 w) | reached :383 |
| `plan = "callr"`, 100 | cluster (2 w) | callr (100 w) | sequential (1 w) | reached :383 |
| `plan = "cluster"` | cluster (2 w) | — | sequential (1 w) | error: `plan_type must be one of: ...` |
| `plan = "multisession"`, `workers = cl` | cluster (2 w) | multisession (1 w) | sequential (1 w) | reached :383 |
| no `plan` element | cluster (2 w) | — | sequential (1 w) | error: `parallel_args$plan must be specified.` |

In no row does the replicate loop see the caller's cluster: the admissible values replace it at
:366, and the rest stop the call.

**Independent check.** A separate read-only review of the same path, not given the probe, was asked
to find any value that keeps the cluster plan through the loop. It found none. It also tested
unusual types against `resolve_bootstrap_parallel_args()` and `setup_parallel_SGcons()` under a
2-node cluster plan:

- These ran and left `sequential`: `factor`, `list` and `FutureStrategy`-classed plan values, an
  environment, an environment with an active binding, a `data.frame`, and
  `workers = cl`.
- These errored inside setup, after `on.exit` is registered: `plan` of length 2 or 0, `NA`, the
  cluster object itself as `X`, and a list without `plan`.
- These errored in the resolver at :279: closures (`future::cluster`, `future::plan()`) and an
  atomic `X`. That is before `on.exit` is registered, so the cluster plan survives, but the call fails.

It confirmed that the `%dofuture%` loop (R/bootstrap_analysis_dofuture.R:408-415) sets no plan of
its own and uses the one in force. Inside each replicate the inner `forestsearch()` is pinned to
`sequential` (:595-597). It also observed that `plan("sequential")` at :355 does not stop a
user-supplied cluster: the processes stay alive, unused.

## 4. Results

Not available: no arm was run.

## 5. Verdict

Not reached: the task stopped at Step 2, so neither MATERIAL nor IMMATERIAL can be stated.

## 6. Notes

- Nothing under `R/`, `tests/`, `man/`, `vignettes/`, `inst/`, `quarto/`, `DESCRIPTION` or
  `NAMESPACE` changed; no install. The repo receives this report and the v3 task file only.
- The v1 patched copy `/tmp/a002_2026-11-10/fs_patched` and v1's `app.R` were not modified.
- A second-route check of every gate above (static parse of the three functions plus the probe
  table and the git state; no replicates) is in `/tmp/a002_2026-11-10/v3/verify_v3.R` →
  `VERIFY_a002_magnitude_v3_2026-11-10.txt`, kept in scratch because the task allows only the
  report in the repo.
