# CC task — bootstrap efficiency Stage 0: measurement only, against CRAN 0.4.0

**File:** `dev/tasks/cc_task_bootstrap_stage0_2026-10-26.md` · **Issued:** 2026-10-26 by chat, commissioned by Larry (decision (b) on `survey_bootstrap_efficiency_revisit_2026-10-25.md`)
**Repo:** `forestsearch` · **Branch:** `master` · **Transport:** `~/Downloads` → `dev/tasks/` (this document, committed first and alone) and `~/Downloads` → `dev/glm-continuous-sims/stage0_bootstrap_2026-10-26/` (the script bundle, §2). Filenames may arrive with hyphens stripped; locate by pattern.

**Audit basis — stated once, carried into the report.** The survey and this task are based on **CRAN `forestsearch` 0.4.0** (Date/Publication 2026-10-10; `Packaged: 2026-10-10 04:12:05 UTC`). Every line reference in the survey is to the 0.4.0 sources as published. The measurements below run against the **installed** package, which `s0` gates at version 0.4.0 and records (`Packaged`, `Repository`), and `s0` additionally compares the working tree's R files under study to the CRAN 0.4.0 MD5 sums so the report can say whether the survey's line references still hold on `master`.

**What this task answers.** Three numbers the survey could not supply from the record:
- **A** — of the per-candidate Cox fit time in the search stage, how much is formula-interface machinery that a direct `survival::coxph.fit()` call would remove, and is that call bit-identical to the current one on every candidate the search fits;
- **C** — how much of the bootstrap's fixed startup a `multicore` plan removes relative to `multisession`, with an `identical()` check of the payload across plans;
- **D** — whether the outer static chunking leaves measurable tail imbalance.

**What this task decides: nothing.** It produces a report with numbers. Whether to commission A (an `R/` change) is Larry's decision afterwards.

---

## ⚠ CATEGORY

**ABSOLUTELY NO PACKAGE CODE IS MODIFIED.** No edit to any file under `R/`, `tests/`, `man/`, `vignettes/`, `inst/`, `src/`, `DESCRIPTION`, `NAMESPACE`, `NEWS.md`, or any driver, application document or payload. No `devtools::install()`, no `load_all()`, no `R CMD check`, no test suite, no vignette build. The scripts call the installed package through its public entry points exactly as a user would.

The one thing to know about and say in the report: `s2` uses **`base::trace()` on two functions in the running R session only** — `forestsearch:::fit_cox_for_subgroup` (to record each candidate's subset vectors as the search fits them) and `survival::coxph.fit` (to record the arguments `coxph()` passes to it) — and `untrace()`s both before any timing. That instruments the session's in-memory copies; it writes nothing to the installed library or the repo. Machine-checkable post-condition in §6 (`git diff` empty outside `dev/`).

**Writes:** this document under `dev/tasks/`; the bundle's seven files, their `results/` and `logs/` under `dev/glm-continuous-sims/stage0_bootstrap_2026-10-26/`; the report `dev/glm-continuous-sims/REPORT_bootstrap_stage0_2026-10-26.md`; and, if `dev/glm-continuous-sims/current_status.md` exists, its regenerated entry.

**Compute (approved with decision (b); stated here as its own line):** on fixture V (gbsg, n = 686): one warm-up + 3 timed + 3 profiled + 1 traced single sequential `forestsearch()` calls; one formula-vs-lean timing loop over the ≈1,410 captured candidate fits (5 reps each path, seconds); and five `forestsearch_bootstrap_dofuture()` runs — sequential nb = 5, multisession nb = 4·W, multicore nb = 4·W, multisession nb = 20, multicore nb = 20, with W = min(48, available cores). Bound from the record: a B = 1000 survival bootstrap at 48 workers measured 3.9 min (`HANDOFF_bootstrap_efficiency_2026-09-02.md`); the five runs total 8·W + 45 ≤ 429 replicates, so `s3` is bounded by roughly that figure plus startups. Hard caps per step are in the driver (`s3`: 40 min). No B = 1000 run. No continuous or binary fixture.

**Unattended.** Gates stop the task; never ask, never work around. The provenance gate is dirt-tolerant: pre-existing untracked files are never staged and never a stop unless they sit under this task's own paths.

**Run from a plain terminal, not RStudio** — `multicore` futures are disabled under RStudio and `s0` stops (exit 3) if the session cannot fork.

---

## 1. Provenance and first commit — GATE

Standard block: `git rev-parse HEAD`, branch (expect **`master`**), `git status --porcelain` (recorded verbatim; dirt tolerated as above), installed `packageVersion("forestsearch")` (expect **0.4.0**; `s0` re-checks and stops otherwise), `packageVersion("survival")`, `R.version.string`. Record the `Packaged:` and `Repository:` fields of the installed package's DESCRIPTION — `Repository: CRAN` means the installed build is the CRAN build; anything else means a local build, which the report must say.

Copy this document to `dev/tasks/cc_task_bootstrap_stage0_2026-10-26.md`, `git add` it by path, commit it alone. Report both filenames (as received in `~/Downloads`, as committed).

## 2. Place and verify the bundle — GATE

The bundle `stage0_bootstrap_2026-10-26.tgz` (also delivered to `~/Downloads`) extracts to one directory. Steps, all machine-checked:

1. `sha256sum ~/Downloads/stage0_bootstrap_2026-10-26.tgz` must print
   `679e984130d3ccbb426792c5be4bb8377edb1b6b4e4efe546d35e70705d17d4c`.
   If it does not (or the file is absent), **do not stop yet**: save the seven sources from Appendix A verbatim into `dev/glm-continuous-sims/stage0_bootstrap_2026-10-26/` and continue with step 3 on them.
2. `tar -xzf ~/Downloads/stage0_bootstrap_2026-10-26.tgz -C dev/glm-continuous-sims/` → creates `dev/glm-continuous-sims/stage0_bootstrap_2026-10-26/` with exactly these seven files.
3. In that directory, `sha256sum *` must print exactly:

```
e1a7b9f97ca2eca363bda5cb2d9ea0eaa5ae0878815b3acf275c8ac96bad16eb  run_stage0.sh
4b45953952511d3a069101b20c5bf3b6d53f921800cee947513c89491fcfa33c  s0_provenance.R
ae8ff345edabe4d39888f445ceb90b076aa1377e99dbcb7477f1f7f9c5488975  s1_profile.R
59f5a4dd4d52007c3df7bef914debda93ed1772aa76f807039e442316be25923  s2_lean_identity.R
4370bf3c3bae7baceda1d50c8ca1993e4ee2efa12c3f75d2d6e6d4a76ea69ee9  s3_bootstrap_runs.R
df583baa090b3d03bd8a36d1c24e5bdecab1c6615974b17a28c0a8ea7aa4557f  s4_summary.R
2af6c94aa0c605e05647f7e1231be05884a5e979b8acf1a469657a295a6d6e30  stage0_common.R
```

   *GATE:* any mismatch after both routes (tarball, then Appendix A) is a stop — report which file and which route.

`dev/glm-continuous-sims/` is where every prior bootstrap report lives (`REPORT_bootstrap_profile_2026-08-30.md` … `REPORT_bootstrap_closeout_2026-09-02.md`). If that directory does not exist on `master`, stop and report — do not invent a location.

## 3. Run

```
cd dev/glm-continuous-sims/stage0_bootstrap_2026-10-26 && bash run_stage0.sh
```

The driver runs `s0` → `s4` in order, each under `timeout` (5 / 15 / 15 / 40 / 5 min), logging to `logs/<step>.log`, writing results to `results/`. A non-zero exit stops the driver; the tail of the failing log says why. The gates inside the scripts:

| step | gate | exit | meaning |
|---|---|---|---|
| s0 | installed `forestsearch` is 0.4.0 | 2 | wrong version installed — this audit is of CRAN 0.4.0 |
| s0 | `parallelly::supportsMulticore()` | 3 | not a plain terminal (RStudio/Windows); relaunch from a terminal |
| s1 | candidate flow on fixture V equals CRAN 0.4.0's: 1,410 candidates reach the Cox fit, 121 pass the effect screen, 120 enter the consistency stage, subgroup `{er <= 0} & {pgr <= 26}` | 5 | the installed package is not computing what 0.4.0 computes on this data |
| s2 | at least one search-stage fit captured | 4 | the search did not run in-process; nothing to measure |

**Script deviations.** The scripts ran end to end in chat's dry run (CRAN 0.4.0 from CRAN, R 4.3.3, survival 3.5.8, future 1.80.0, 2 cores, Linux). If a step fails for a reason that is a script defect rather than a gate (a changed `survival` signature, a missing field), fix the script with the smallest edit that restores the stated measurement, re-run from that step, and list every edit verbatim (file, old line, new line) in the report's "Script deviations" section. Never change what is measured or how identity is judged; if that would be needed, stop and report.

## 4. What each step measures (the script headers say the same)

- **s0** provenance; the two gates; repo R files under study vs CRAN 0.4.0 MD5 (recorded, never a stop; mismatches get their `git log --since=2026-10-09`).
- **s1** fixture V = the package's own `survival-post-selection` vignette gbsg call with `mr_inference = FALSE` and `parallel_args = list(plan = "sequential")` (the per-replicate shape). Warm-up, candidate-flow counts, 3 timed calls, 3 `Rprof` profiles; reports total time by function with `pkg::` prefixes stripped. `coxph` − `coxph.fit` is the formula-interface machinery; `concordancefit` is reported separately because `coxph()` computes it on every fit and nothing in the search consumes it.
- **s2** candidate A: captures the (Y, E, Treat) subset of every candidate reaching `fit_cox_for_subgroup()` during one call; captures the arguments `coxph()` passes to `coxph.fit()`; for every captured candidate runs the package's own block (`R/subgroup_search.R:834–845` in 0.4.0, transplanted verbatim) and a lean `coxph.fit()` block built from the captured conventions; `identical()` on the consumed CI row, agreement of error and warning sets, and 5 alternating timing reps of each path. Output: machinery seconds per call and the identity count.
- **s3** candidates C and D: `fs.est` fitted once; five bootstrap runs as listed under Compute; per run the decomposition outer = spawn + inner, inner = overhead + static-chunk makespan (from the bootstrap's own `tmins_iteration` column and `parallel::splitIndices`), balanced makespan, imbalance; per-replicate slowdown under load vs sequential; `identical()` of the nb = 20 payloads across plans with the volatile timing columns/attribute excluded (and only those).
- **s4** assembles `results/stage0_summary.md` from `results/keyvalues.txt` — every figure in the summary is a line in that file. The 30% bar is applied to the per-replicate (campaign) shape; the application-shape (B = 1000) figures are projections from measured components and are labelled as such.

## 5. Report — `dev/glm-continuous-sims/REPORT_bootstrap_stage0_2026-10-26.md`

1. **Audit basis** (first lines): CRAN `forestsearch` 0.4.0; installed build fields; whether the installed build is the CRAN build; repo HEAD and branch; the MD5 table from `results/provenance_md5.csv` and, for any file that differs from CRAN 0.4.0, its git log since the release and one line on whether the survey's line references to it still hold (`git diff` the file against HEAD~N if needed — read-only).
2. `results/stage0_summary.md` **transplanted verbatim** (it is the measurement record).
3. **Against the 0.3.5 record** (`HANDOFF_bootstrap_efficiency_2026-09-02.md` §1 and §4), one table, measured-now beside recorded: single gbsg call (record 3.58 s), sequential replicate (4.18 s), candidates reaching the fit / passing the screen (1,410 / 121), "coxph internals" (≈2.4 s per call), startup under multisession (≈69 s fixed at 48 workers), B = 1000 at 48 workers (3.9 min measured). Note which of these are like-for-like (fixture V reproduces the record's 1,410 / 121 candidate flow) and which are not (machine, `survival` version, W if below 48).
4. **Verdicts, numbers only** — A: machinery seconds per call, share of the measured sequential replicate, bar cleared yes/no, identity count, the captured `coxph.fit()` conventions (the recipe an implementation task would have to reproduce: `nocenter`, `aeqSurv` under `timefix`, `init`, `offset`, `rownames`, `resid`). C: startup removed, measured reduction at nb = 4·W, projected reduction at B = 1000, identity across plans. D: imbalance seconds and share. B (not separately measured): the Rprof line "coxph time outside the search stage" sizes the consistency-stage fits under full evaluation — report it as what it is. **No recommendation to implement anything**; that is the next decision, Larry's.
5. **Script deviations** (verbatim edits, or "none").
6. **OPEN ITEMS** (documentation gaps; never a stop).
7. **Ten-line verdict.**

Stop-point figures must be quoted verbatim from `results/keyvalues.txt` / `results/stage0_summary.md`; an empty or paraphrased figure is ABSENT.

## 6. Commit and post-conditions — GATE

Commit by explicit path, in this order, no push, no fetch:
1. (already done in §1) the task document alone;
2. `dev/glm-continuous-sims/stage0_bootstrap_2026-10-26/` — the seven sources, `results/` (all files, including the three `Rprof_call_*.out` and the `.rds`; total is a few MB), `logs/`; and the report `dev/glm-continuous-sims/REPORT_bootstrap_stage0_2026-10-26.md`;
3. if `dev/glm-continuous-sims/current_status.md` exists: its entry for this task and the regenerated pin, as the final commit, with the standing post-condition that the stated pin equals `HEAD` at commit time. If it does not exist, note that in OPEN ITEMS and do not create one.

Before each commit: `git diff --cached --name-only` must list only paths under `dev/`; abort otherwise. After the final commit, these must all hold and be printed in the report's closing block:

```
git diff --name-only <HEAD recorded in §1> HEAD | grep -v '^dev/'      # must print nothing
git status --porcelain | grep -v '^?? ' | grep -v '^.. dev/'            # must print nothing (no tracked file modified)
git status --porcelain | grep '^?? ' | grep -v ' dev/'                   # only pre-existing dirt, listed in §1
```

The first line is the machine check of "no package code modified". If it prints anything, do not fix it — stop and report what it prints.

## 7. Out of scope

- No `R/` edit, no estimator change, no driver, application, document or payload change. No install, no check, no tests, no vignette build. No push, no fetch.
- No B = 1000 run; no continuous (MD) or binary (OR) fixture; no re-profile of either path.
- No change to `parallel_args` defaults, chunking, or batch sizing — D is measured from recorded timings, not by altering the loop.
- No re-running of anything in the 09-02 record.

---

## Appendix A — script sources, verbatim (fallback for §2 only)

Save each block to the named file under `dev/glm-continuous-sims/stage0_bootstrap_2026-10-26/`, byte for byte (no trailing-whitespace cleanup, Unix newlines, a final newline), then re-check §2 step 3. The tarball and these blocks are the same bytes.

### `run_stage0.sh`

```bash
#!/usr/bin/env bash
# run_stage0.sh -- driver for the bootstrap-efficiency Stage-0 measurement.
# Measurement only: no package file is read for writing, let alone changed.
# Run from a plain terminal (multicore is disabled under RStudio):
#     cd dev/glm-continuous-sims/stage0_bootstrap_2026-10-26 && bash run_stage0.sh
# Each step has a hard time cap; a non-zero exit stops the driver.
set -uo pipefail
cd "$(dirname "$0")"
export STAGE0_OUT=results
mkdir -p results logs
run() {
  local name=$1 cap=$2; shift 2
  echo "== $name (cap $cap)"
  timeout -k 15s "$cap" "$@" > "logs/$name.log" 2>&1
  local rc=$?
  tail -n 3 "logs/$name.log" | sed 's/^/   /'
  if [ $rc -ne 0 ]; then echo "STOP: $name exited $rc -- see logs/$name.log"; exit $rc; fi
}
run s0_provenance 5m  Rscript s0_provenance.R
run s1_profile    15m Rscript s1_profile.R
run s2_lean       15m Rscript s2_lean_identity.R
run s3_bootstrap  40m Rscript s3_bootstrap_runs.R
run s4_summary    5m  Rscript s4_summary.R
echo "== all steps complete; summary at results/stage0_summary.md"
```

### `stage0_common.R`

```r
# stage0_common.R -- shared setup for the bootstrap-efficiency Stage-0
# measurement (task cc_task_bootstrap_stage0_2026-10-26.md).
#
# MEASUREMENT ONLY.  Nothing here edits any file of the forestsearch package,
# on disk or in the installed library.  s2 uses base::trace() on two functions
# IN THE RUNNING SESSION (forestsearch:::fit_cox_for_subgroup and
# survival::coxph.fit) to record their arguments, and untrace()s them before
# any timing; the installed package files are untouched.
#
# Audit basis: CRAN forestsearch 0.4.0 (Date/Publication 2026-10-10).

suppressPackageStartupMessages({
  library(forestsearch)
  library(survival)
  library(future)
  library(doFuture)
  library(data.table)
})

OUT <- Sys.getenv("STAGE0_OUT", "results")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

out_path <- function(...) file.path(OUT, ...)

save_result <- function(x, name) {
  saveRDS(x, out_path(paste0(name, ".rds")))
  invisible(x)
}

# key = value lines, appended to results/keyvalues.txt (the report quotes these)
kv <- function(...) {
  vals <- list(...)
  stopifnot(!is.null(names(vals)), all(nzchar(names(vals))))
  lines <- vapply(names(vals), function(k) {
    v <- vals[[k]]
    v <- if (is.numeric(v)) paste(format(v, digits = 6, trim = TRUE), collapse = ", ")
         else paste(as.character(v), collapse = ", ")
    sprintf("%s = %s", k, v)
  }, character(1))
  cat(lines, sep = "\n")
  cat(lines, file = out_path("keyvalues.txt"), sep = "\n", append = TRUE)
  invisible(lines)
}

# --- Fixture V: the gbsg configuration of the package's own vignette ---------
# inst/doc/survival-post-selection.qmd (CRAN 0.4.0), chunks `data` and `fit`,
# with two changes that give the per-replicate shape of the bootstrap:
#   * mr_inference = FALSE (bootstrap_results() strips mr_inference from every
#     replicate under the default mr_in_replicates = FALSE), and
#     mr_inference_args dropped with it;
#   * parallel_args = list(plan = "sequential") -- as the vignette already has,
#     and as every bootstrap replicate forces (CATEGORY 3 in bootstrap_results()).
# Everything else is the vignette's call verbatim.
stage0_fixture <- function() {
  df.analysis <- survival::gbsg
  df.analysis <- within(df.analysis, {
    id          <- seq_len(nrow(df.analysis))
    time_months <- rfstime / 30.4375
    grade3      <- ifelse(grade == "3", 1, 0)
  })
  confounders.name <- c("age", "meno", "size", "grade3", "nodes", "pgr", "er")
  args <- list(
    df.analysis = df.analysis,
    outcome.name = "time_months", event.name = "status",
    treat.name = "hormon", id.name = "id",
    confounders.name = confounders.name,
    is.RCT = TRUE, seedit = 8316951, est.scale = "hr",
    use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE,
    conf.cont_jcuts = list(er = 10, pgr = 10), collapse_cuts = TRUE,
    conf_force = c("er <= 0", "pgr <= 0"),
    subgroup_method = "consistency", sg_focus = "effMaxSG",
    selection_rule = "neighborhood", effect_neighborhood = 0.20,
    mr_inference = FALSE,
    n.min = 60, d0.min = 10, d1.min = 10, maxk = 2,
    hr.threshold = 1.0, hr.consistency = 1.0, pconsistency.threshold = 0.80,
    fs.splits = 200, max.minutes = 3,
    details = FALSE, quiet = TRUE, plot.sg = FALSE,
    parallel_args = list(plan = "sequential")
  )
  list(label = "V: survival-post-selection vignette gbsg config, mr_inference = FALSE",
       df = df.analysis, args = args)
}

run_fixture_call <- function(fx) {
  suppressWarnings(suppressMessages(do.call(forestsearch::forestsearch, fx$args)))
}

# Candidate-flow counters of one forestsearch() result (survival consistency path)
candidate_counts <- function(fs) {
  fc <- fs$find.grps$filter_counts
  gc <- fs$grp.consistency
  list(
    n_evaluated_search   = fc$n_evaluated,
    n_reached_cox_fit    = fc$n_passed_sample_size,   # status >= 5: fit attempted
    n_cox_fit_ok         = fc$n_passed_cox,           # status >= 6
    n_passed_effect      = fc$n_passed_hr,            # status 7 (medians computed here)
    n_consistency_total  = gc$n_candidates_total,
    n_consistency_evaluated = gc$n_candidates_evaluated,
    n_consistency_passed = gc$n_passed,
    early_stop_triggered = isTRUE(gc$early_stop_triggered),
    sg_harm = paste(fs$sg.harm, collapse = " & ")
  )
}
```

### `s0_provenance.R`

```r
# s0_provenance.R -- provenance record and the two gates that must hold before
# anything is measured.  Exits non-zero on a failed gate (the driver stops).
#
# GATE 1  installed forestsearch is version 0.4.0 (the CRAN release this audit
#         is based on).
# GATE 2  multicore is usable in this session (candidate C cannot be measured
#         otherwise; relaunch from a plain terminal).
# RECORD  repo R/ files under study vs their CRAN 0.4.0 MD5 sums (from the
#         CRAN tarball's MD5 file).  A mismatch is NOT a stop: the measurement
#         runs against the installed package, and GATE 1 covers that.  Every
#         mismatch is listed with its git log since the release so the report
#         can say whether the survey's line references still hold.

source("stage0_common.R")
unlink(out_path("keyvalues.txt"))

# --- installed package -------------------------------------------------------
pd <- utils::packageDescription("forestsearch")
kv(installed_version   = as.character(pd$Version),
   installed_packaged  = if (is.null(pd$Packaged)) "<absent>" else pd$Packaged,
   installed_built     = if (is.null(pd$Built)) "<absent>" else pd$Built,
   installed_repository = if (is.null(pd$Repository)) "<absent>" else pd$Repository,
   installed_remotetype = if (is.null(pd$RemoteType)) "<absent>" else pd$RemoteType,
   installed_remotesha  = if (is.null(pd$RemoteSha)) "<absent>" else pd$RemoteSha,
   installed_libpath   = dirname(dirname(dirname(attr(pd, "file")))),
   survival_version    = as.character(utils::packageVersion("survival")),
   future_version      = as.character(utils::packageVersion("future")),
   doFuture_version    = as.character(utils::packageVersion("doFuture")),
   r_version           = R.version.string,
   sysname             = unname(Sys.info()[["sysname"]]),
   machine             = unname(Sys.info()[["machine"]]),
   available_cores     = parallelly::availableCores(),
   supports_multicore  = parallelly::supportsMulticore(),
   in_rstudio          = nzchar(Sys.getenv("RSTUDIO")))

if (as.character(pd$Version) != "0.4.0") {
  cat("GATE 1 FAILED: installed forestsearch is", as.character(pd$Version),
      "- this audit is based on CRAN 0.4.0.\n")
  quit(status = 2)
}
if (!isTRUE(parallelly::supportsMulticore())) {
  cat("GATE 2 FAILED: multicore futures are not supported in this session",
      "(RStudio or Windows).  Relaunch run_stage0.sh from a plain terminal.\n")
  quit(status = 3)
}

# --- repo files vs CRAN 0.4.0 ------------------------------------------------
cran_md5 <- c(
  "DESCRIPTION"                           = "19b9459a03ed5c4a0c12943eb57d82b6",
  "NAMESPACE"                             = "00257d44d6f21d10f9af68b9f0dc227c",
  "R/Cox_estimation_helpers.R"            = "ffa9744bf4b3264041aa4c787cd5958f",
  "R/bootstrap_analysis_dofuture.R"       = "a891d66e7337b9a1ebdf524463102978",
  "R/bootstrap_calculations_helpers.R"    = "2c25b9b63263653c599f7143df6b854c",
  "R/bootstrap_dofuture_main.R"           = "65cb60f5b23303578362172bc3cb9291",
  "R/bootstrap_dofuture_setup_helpers.R"  = "97cac83ad17bff5bf16efe7f437e164a",
  "R/consistency_resample.R"              = "861c6b4e7d00c8176bbab6b7f863f724",
  "R/forestsearch_main.R"                 = "11d80e74c4fcdba43fe8bbf9a1287d3f",
  "R/get_fsdata.R"                        = "a7f71d659558706ae121d72e82099181",
  "R/guohe_algorithm3.R"                  = "5a6e842f871c61ddef41dbfcc93bf6c5",
  "R/subgroup_consistency_helpers.R"      = "f32b55b62d683337f6e1b36b3aa3bc17",
  "R/subgroup_consistency_main.R"         = "1285150ed4891883468d01d5ee981b78",
  "R/subgroup_search.R"                   = "f06d6bba8135f41a57962a659ebbe601"
)
root <- tryCatch(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE, stderr = FALSE),
                 error = function(e) character(0), warning = function(w) character(0))
if (length(root) != 1L || !nzchar(root)) {
  cat("not inside a git repository; repo-vs-CRAN comparison skipped\n")
  kv(repo_root = "<not a git repo>")
} else {
  head_sha <- system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE)
  branch   <- system2("git", c("-C", shQuote(root), "rev-parse", "--abbrev-ref", "HEAD"), stdout = TRUE)
  kv(repo_root = root, repo_head = head_sha, repo_branch = branch)
  files <- file.path(root, names(cran_md5))
  have  <- file.exists(files)
  got   <- ifelse(have, unname(tools::md5sum(files)), NA_character_)
  tab <- data.frame(file = names(cran_md5), cran_0_4_0 = unname(cran_md5),
                    repo_head = got,
                    match = ifelse(is.na(got), "MISSING", ifelse(got == unname(cran_md5), "match", "DIFFERS")),
                    stringsAsFactors = FALSE)
  utils::write.csv(tab, out_path("provenance_md5.csv"), row.names = FALSE)
  print(tab[, c("file", "match")], row.names = FALSE)
  kv(repo_files_matching_cran = sum(tab$match == "match"),
     repo_files_differing     = sum(tab$match != "match"))
  diff_files <- tab$file[tab$match != "match"]
  if (length(diff_files)) {
    cat("\nFiles differing from CRAN 0.4.0 -- git log since the release (2026-10-09):\n")
    for (f in diff_files) {
      cat("--", f, "\n")
      lg <- system2("git", c("-C", shQuote(root), "log", "--oneline", "--since=2026-10-09", "--", shQuote(f)),
                    stdout = TRUE)
      cat(if (length(lg)) paste0("   ", lg) else "   (no commits since 2026-10-09; file predates or is untracked)", sep = "\n")
    }
    cat("NOTE: measurements run against the INSTALLED package (GATE 1 = 0.4.0).",
        "If the installed build is not the CRAN build (installed_repository != CRAN),",
        "say so in the report's provenance section.\n")
  }
}
cat("\ns0 complete\n")
```

### `s1_profile.R`

```r
# s1_profile.R -- the single sequential forestsearch() call on fixture V:
# wall-clock (3 timed calls after a warm-up) and an Rprof split (3 profiles)
# attributing the call to the Cox fit machinery vs the fit proper.
#
# Nothing is modified; the installed package is called as a user would.

source("stage0_common.R")
fx <- stage0_fixture()
kv(fixture = fx$label, fixture_n = nrow(fx$df))

# warm-up (package code loading, lazy-load DB, JIT) -- recorded, not used
t_warm <- system.time(fs <- run_fixture_call(fx))[["elapsed"]]
kv(call_warmup_s = t_warm)
cc <- candidate_counts(fs)
kv(search_n_evaluated     = cc$n_evaluated_search,
   search_n_reached_cox   = cc$n_reached_cox_fit,
   search_n_cox_ok        = cc$n_cox_fit_ok,
   search_n_passed_effect = cc$n_passed_effect,
   consistency_n_total    = cc$n_consistency_total,
   consistency_n_evaluated = cc$n_consistency_evaluated,
   consistency_n_passed   = cc$n_consistency_passed,
   consistency_early_stop = cc$early_stop_triggered,
   identified_subgroup    = cc$sg_harm)

# GATE: the candidate flow is fixed by the data and the 0.4.0 code, not by the
# machine.  These are the counts CRAN 0.4.0 produced on fixture V in the chat's
# dry run (and the 1,410 / 121 the 0.3.5 record carries).  A mismatch means the
# installed package is not computing what 0.4.0 computes here -- stop.
expected <- list(n_reached_cox_fit = 1410L, n_passed_effect = 121L,
                 n_consistency_total = 120L, sg_harm = "{er <= 0} & {pgr <= 26}")
bad <- names(expected)[!mapply(function(k, v) identical(as.character(cc[[k]]), as.character(v)),
                               names(expected), expected)]
if (length(bad)) {
  cat("GATE FAILED: candidate flow differs from CRAN 0.4.0 on fixture V for:",
      paste(bad, collapse = ", "), "\n")
  quit(status = 5)
}

# --- plain timings -----------------------------------------------------------
t_call <- vapply(1:3, function(i) system.time(run_fixture_call(fx))[["elapsed"]], numeric(1))
kv(call_s_rep = t_call, call_s_median = median(t_call))

# --- Rprof split -------------------------------------------------------------
# Functions whose total time is reported.  "coxph" is the formula interface;
# "coxph.fit" is the fit proper it delegates to (the lean path of candidate A
# calls it directly).  Their difference is the formula-interface machinery:
# model.frame / model.matrix / Surv / aeqSurv / concordancefit / object assembly.
watch <- c("forestsearch", "subgroup.search", "search_combinations_parallel",
           "evaluate_combination_with_status", "fit_cox_for_subgroup",
           "coxph", "coxph.fit", "agreg.fit", "model.frame.default",
           "model.matrix.default", "Surv", "aeqSurv", "concordancefit",
           "concordance", "survfit", "survfit.formula", "summary.survfit",
           "subgroup.consistency", "evaluate_consistency_twostage",
           "evaluate_subgroup_consistency", "consistency_resample",
           ".consistency_cox_pieces", "residuals.coxph",
           "remove_near_duplicate_subgroups", "sort_subgroups_preview",
           "get_FSdata", "dummy", "get_dfpred", "merge", "future_lapply")

prof_one <- function(i) {
  f <- out_path(sprintf("Rprof_call_%d.out", i))
  # the platform may raise a too-short interval (e.g. to 0.01 s); harmless
  suppressWarnings(Rprof(f, interval = 0.002))
  run_fixture_call(fx)
  Rprof(NULL)
  sp <- summaryRprof(f)
  bt <- sp$by.total
  # Rprof labels namespaced calls "pkg::fun"; strip the prefix so "survival::coxph"
  # and "coxph" are the same row.  Functions reached through do.call() show as
  # "<Anonymous>" and are not attributable by name (forestsearch(), get_FSdata(),
  # subgroup.consistency()) -- their children are.
  bt$fun <- sub("^[A-Za-z0-9._]+::", "", gsub('"', "", rownames(bt)))
  bt <- aggregate(cbind(total.time, total.pct, self.time, self.pct) ~ fun, data = bt, FUN = sum)
  bt$rep <- i
  bt$sampling_time_s <- sp$sampling.time
  bt
}
prof <- do.call(rbind, lapply(1:3, prof_one))
utils::write.csv(prof, out_path("rprof_by_total_all.csv"), row.names = FALSE)

# median over reps for the watched functions
w <- prof[prof$fun %in% watch, c("fun", "total.time", "total.pct", "rep")]
med <- aggregate(cbind(total.time, total.pct) ~ fun, data = w, FUN = median)
med <- med[order(-med$total.time), ]
med$fun <- factor(med$fun, levels = watch)
med <- med[order(med$fun), ]; med$fun <- as.character(med$fun)
utils::write.csv(med, out_path("rprof_watched_median.csv"), row.names = FALSE)
cat("\nRprof by.total, median of 3 profiles (seconds, % of samples):\n")
print(med, row.names = FALSE)

get_t <- function(f) { v <- med$total.time[med$fun == f]; if (length(v)) v else 0 }
cox_total <- get_t("coxph")
fit_total <- get_t("coxph.fit") + get_t("agreg.fit")
kv(rprof_sampling_s_median   = median(prof$sampling_time_s[!duplicated(prof$rep)]),
   rprof_coxph_total_s       = cox_total,
   rprof_coxph_fit_total_s   = fit_total,
   rprof_formula_machinery_s = cox_total - fit_total,
   rprof_machinery_share_of_coxph = if (cox_total > 0) (cox_total - fit_total) / cox_total else NA,
   rprof_fit_cox_for_subgroup_s = get_t("fit_cox_for_subgroup"),
   rprof_search_stage_s      = get_t("search_combinations_parallel"),
   rprof_survfit_s           = get_t("survfit"),
   rprof_concordancefit_s    = get_t("concordancefit"),
   rprof_consistency_stage_s = get_t("evaluate_consistency_twostage") + get_t("evaluate_subgroup_consistency"),
   rprof_consistency_cox_pieces_s = get_t(".consistency_cox_pieces"),
   rprof_residuals_coxph_s   = get_t("residuals.coxph"))
save_result(list(t_warm = t_warm, t_call = t_call, counts = cc, prof = prof, med = med), "s1_profile")
cat("\ns1 complete\n")
```

### `s2_lean_identity.R`

```r
# s2_lean_identity.R -- candidate A sized and its bit-identity tested, with
# NO change to the package:
#
#  1. base::trace() on forestsearch:::fit_cox_for_subgroup (in this session
#     only) records the subset vectors (Y, E, Treat) of every candidate that
#     reaches the search-stage Cox fit during ONE forestsearch() call on
#     fixture V.  untrace() follows immediately.
#  2. base::trace() on survival::coxph.fit records the arguments coxph()
#     passes to it during ONE formula fit, so the lean call below uses the
#     installed survival's conventions rather than assumptions.  untrace().
#  3. For every captured candidate: the package's formula block
#     (R/subgroup_search.R lines 834-845, transplanted verbatim) and a lean
#     block that calls survival::coxph.fit() directly.  identical() on the
#     consumed CI row; agreement of the failure sets; timings (5 reps,
#     alternating order).
#
# Output: the machinery seconds the lean path would remove per call, the
# identity count, and the captured coxph.fit() conventions (the recipe).

source("stage0_common.R")
fx <- stage0_fixture()

# ---- 1. capture the search-stage Cox-fit inputs ------------------------------
.cap <- new.env()
.cap$n <- 0L
.cap$items <- vector("list", 20000L)
.cap$push <- function(item) {
  .cap$n <- .cap$n + 1L
  if (.cap$n > length(.cap$items)) length(.cap$items) <- 2L * length(.cap$items)
  .cap$items[[.cap$n]] <- item
}
assign(".cap", .cap, envir = globalenv())

stopifnot(exists("fit_cox_for_subgroup", envir = asNamespace("forestsearch"), inherits = FALSE))
invisible(capture.output(
  trace("fit_cox_for_subgroup", where = asNamespace("forestsearch"), print = FALSE,
        tracer = quote({
          .GlobalEnv$.cap$push(list(
            Y = yy[id.x == 1], E = dd[id.x == 1], Treat = tt[id.x == 1],
            adj = forestsearch:::.fs_adjust_terms(adjust_covariates)))
        }))
))
on.exit(try(suppressMessages(untrace("fit_cox_for_subgroup", where = asNamespace("forestsearch"))), silent = TRUE), add = TRUE)
t_traced <- system.time(fs <- run_fixture_call(fx))[["elapsed"]]
invisible(capture.output(untrace("fit_cox_for_subgroup", where = asNamespace("forestsearch"))))

K <- .cap$n
cands <- .cap$items[seq_len(K)]
cc <- candidate_counts(fs)
kv(capture_n_fits = K,
   capture_matches_filter_counts = identical(K, as.integer(cc$n_reached_cox_fit)),
   capture_n_adjusted = sum(vapply(cands, function(z) length(z$adj) > 0L, logical(1))),
   capture_median_n = median(vapply(cands, function(z) length(z$Y), numeric(1))),
   capture_call_s_under_trace = t_traced)
if (K == 0L) {
  cat("GATE FAILED: no fit_cox_for_subgroup() call was captured -- the search did not run",
      "in this process (future fallback changed?).  Nothing to measure.\n")
  quit(status = 4)
}
if (any(vapply(cands, function(z) length(z$adj) > 0L, logical(1)))) {
  cat("NOTE: adjusted candidates present; the lean path applies to the unadjusted branch only.",
      "They are excluded from the comparison below and counted above.\n")
  cands <- Filter(function(z) length(z$adj) == 0L, cands)
}
cands <- Filter(function(z) length(z$Y) > 0L, cands)
K <- length(cands)

# ---- 2. capture what coxph() passes to coxph.fit() -----------------------------
.capfit <- new.env(); assign(".capfit", .capfit, envir = globalenv())
# Every formal of the installed coxph.fit() is captured by name, so a formal
# this script does not know about (a newer survival) is still recorded and
# reused below rather than silently left at its default.
invisible(capture.output(
  trace("coxph.fit", where = asNamespace("survival"), print = FALSE,
        tracer = quote({
          .fr <- environment()
          .nm <- setdiff(names(formals(survival::coxph.fit)), "...")
          .GlobalEnv$.capfit$args <- lapply(setNames(.nm, .nm), function(.a)
            tryCatch(get(.a, envir = .fr, inherits = FALSE),
                     error = function(e) structure(list(), class = "stage0_capture_failed")))
        }))
))
probe <- cands[[which.max(vapply(cands, function(z) length(z$Y), numeric(1)))]]
Y <- probe$Y; E <- probe$E; Treat <- probe$Treat
fit_probe <- survival::coxph(survival::Surv(Y, E) ~ Treat, robust = FALSE)
invisible(capture.output(untrace("coxph.fit", where = asNamespace("survival"))))
a <- .capfit$args
stopifnot(!is.null(a), all(c("x", "y") %in% names(a)))
failed_caps <- names(a)[vapply(a, inherits, logical(1), "stage0_capture_failed")]
a <- a[setdiff(names(a), failed_caps)]
kv(conv_formals_not_capturable = if (length(failed_caps)) paste(failed_caps, collapse = ",") else "<none>")
ctl <- survival::coxph.control()
fit_formals <- names(formals(survival::coxph.fit))
y_raw <- survival::Surv(Y, E)
y_fix <- if (isTRUE(ctl$timefix)) survival::aeqSurv(y_raw) else y_raw
same_vals <- function(u, v) isTRUE(all.equal(unclass(unname(u)), unclass(unname(v)), check.attributes = FALSE))
conv <- list(
  coxph_fit_formals   = paste(fit_formals, collapse = ","),
  x_storage           = storage.mode(a$x),
  x_colnames          = paste(colnames(a$x), collapse = ","),
  x_values_equal_treat = same_vals(a$x[, 1L], as.double(Treat)),
  y_type              = attr(a$y, "type"),
  y_values_equal_raw  = same_vals(a$y, y_raw),
  y_values_equal_aeqSurv = same_vals(a$y, y_fix),
  timefix             = isTRUE(ctl$timefix),
  strata_is_null      = is.null(a$strata),
  offset_desc         = if (is.null(a$offset)) "NULL" else sprintf("numeric length %d, all zero: %s", length(a$offset), all(a$offset == 0)),
  init_desc           = if (is.null(a$init)) "NULL" else paste(format(a$init), collapse = ","),
  control_is_default  = identical(a$control, ctl),
  weights_is_null     = is.null(a$weights),
  method              = a$method,
  rownames_desc       = if (is.null(a$rownames)) "NULL" else sprintf("%s length %d", class(a$rownames)[1], length(a$rownames)),
  resid               = if (!"resid" %in% names(a)) "<no formal>" else as.character(a$resid),
  nocenter            = if (!"nocenter" %in% names(a)) "<no formal>" else if (is.null(a$nocenter)) "NULL" else paste(a$nocenter, collapse = ","),
  other_formals       = {
    oth <- setdiff(names(a), c("x", "y", "strata", "offset", "init", "control", "weights",
                               "method", "rownames", "resid", "nocenter"))
    if (length(oth)) paste(vapply(oth, function(k) paste0(k, "=", paste(deparse(a[[k]]), collapse = "")), character(1)), collapse = "; ") else "<none>"
  }
)
do.call(kv, setNames(conv, paste0("conv_", names(conv))))

# ---- 3. the two blocks ---------------------------------------------------------
# (a) the package's own block, R/subgroup_search.R:834-845 verbatim (CRAN 0.4.0)
formula_block <- function(Y, E, Treat) {
  try({
    fit.cox <- survival::coxph(survival::Surv(Y, E) ~ Treat, robust = FALSE)
    beta <- fit.cox$coefficients
    se_beta <- sqrt(diag(fit.cox$var))
    z_ci <- stats::qnorm((1 + 0.95) / 2, 0, 1)
    ci <- cbind(exp(beta), exp(-beta),
                exp(beta - z_ci * se_beta), exp(beta + z_ci * se_beta))
    dimnames(ci) <- list(names(beta),
                         c("exp(coef)", "exp(-coef)",
                           "lower .95", "upper .95"))
    ci
  }, silent = TRUE)
}
# (b) the lean block: coxph.fit() with the conventions captured in step 2.
#     Data-dependent arguments (x, y, offset, rownames) are rebuilt per
#     candidate; every other captured formal (strata, init, control, weights,
#     method, resid, nocenter, and anything a newer survival adds) is reused
#     exactly as captured.  names(beta) <- "Treat" restores the coefficient
#     name coxph() takes from the model-matrix column, so the CI row's
#     dimnames agree.
lean_args_fixed <- a[setdiff(names(a), c("x", "y", "offset", "rownames"))]
use_offset_zero <- !is.null(a$offset)
use_rownames    <- !is.null(a$rownames)
use_aeqSurv     <- isTRUE(ctl$timefix)
lean_block <- function(Y, E, Treat) {
  try({
    n <- length(Y)
    ys <- survival::Surv(Y, E)
    if (use_aeqSurv) ys <- survival::aeqSurv(ys)
    args <- c(list(x = matrix(as.double(Treat), ncol = 1L),
                   y = ys,
                   offset = if (use_offset_zero) rep(0, n) else NULL,
                   rownames = if (use_rownames) as.character(seq_len(n)) else NULL),
              lean_args_fixed)
    fit.cox <- do.call(survival::coxph.fit, args)
    beta <- fit.cox$coefficients
    names(beta) <- "Treat"
    se_beta <- sqrt(diag(fit.cox$var))
    z_ci <- stats::qnorm((1 + 0.95) / 2, 0, 1)
    ci <- cbind(exp(beta), exp(-beta),
                exp(beta - z_ci * se_beta), exp(beta + z_ci * se_beta))
    dimnames(ci) <- list(names(beta),
                         c("exp(coef)", "exp(-coef)",
                           "lower .95", "upper .95"))
    ci
  }, silent = TRUE)
}

run_all <- function(block) {
  res <- vector("list", K); nwarn <- integer(K)
  for (k in seq_len(K)) {
    z <- cands[[k]]
    res[[k]] <- withCallingHandlers(block(z$Y, z$E, z$Treat),
                 warning = function(w) { nwarn[k] <<- nwarn[k] + 1L; invokeRestart("muffleWarning") })
  }
  list(res = res, nwarn = nwarn)
}

# identity pass (unordered timing -- separate from the timed reps below)
ra <- run_all(formula_block); rb <- run_all(lean_block)
err_a <- vapply(ra$res, inherits, logical(1), "try-error")
err_b <- vapply(rb$res, inherits, logical(1), "try-error")
ident <- mapply(function(u, v) identical(u, v), ra$res, rb$res)
ident_vals <- mapply(function(u, v) {
  if (inherits(u, "try-error") || inherits(v, "try-error")) return(NA)
  identical(unname(u), unname(v)) }, ra$res, rb$res)
kv(ident_n_candidates = K,
   ident_ci_row_identical = sum(ident),
   ident_ci_values_identical = sum(ident_vals, na.rm = TRUE),
   ident_errors_formula = sum(err_a), ident_errors_lean = sum(err_b),
   ident_error_sets_agree = identical(err_a, err_b),
   ident_warn_formula_total = sum(ra$nwarn), ident_warn_lean_total = sum(rb$nwarn),
   ident_warn_sets_agree = identical(ra$nwarn > 0L, rb$nwarn > 0L),
   ident_identical_among_warned = if (any(ra$nwarn > 0L)) sum(ident[ra$nwarn > 0L]) else NA,
   ident_n_warned = sum(ra$nwarn > 0L))
if (any(!ident)) {
  k <- which(!ident)[1L]
  cat("\nFirst non-identical candidate (k =", k, "): n =", length(cands[[k]]$Y), "\n")
  cat("formula:\n"); print(ra$res[[k]]); cat("lean:\n"); print(rb$res[[k]])
  if (!inherits(ra$res[[k]], "try-error") && !inherits(rb$res[[k]], "try-error"))
    cat("max abs diff:", max(abs(ra$res[[k]] - rb$res[[k]])), "\n")
}

# ---- timings: 5 reps, alternating order --------------------------------------
time_block <- function(block) {
  system.time(for (k in seq_len(K)) { z <- cands[[k]]; suppressWarnings(block(z$Y, z$E, z$Treat)) })[["elapsed"]]
}
ta <- tb <- numeric(5)
for (r in 1:5) {
  if (r %% 2L == 1L) { ta[r] <- time_block(formula_block); tb[r] <- time_block(lean_block) }
  else               { tb[r] <- time_block(lean_block);    ta[r] <- time_block(formula_block) }
}
kv(time_formula_s_reps = ta, time_lean_s_reps = tb,
   time_formula_s_median = median(ta), time_lean_s_median = median(tb),
   time_ms_per_fit_formula = 1000 * median(ta) / K,
   time_ms_per_fit_lean = 1000 * median(tb) / K,
   machinery_s_per_call = median(ta) - median(tb),
   machinery_share_of_formula = (median(ta) - median(tb)) / median(ta))

save_result(list(K = K, conv = conv, capfit = a, ident = ident, ident_vals = ident_vals,
                 err_a = err_a, err_b = err_b, nwarn_a = ra$nwarn, nwarn_b = rb$nwarn,
                 ta = ta, tb = tb, cand_sizes = vapply(cands, function(z) length(z$Y), numeric(1))),
            "s2_lean_identity")
cat("\ns2 complete\n")
```

### `s3_bootstrap_runs.R`

```r
# s3_bootstrap_runs.R -- the application-shape measurements for candidates C
# (multicore vs multisession startup) and D (tail imbalance of the outer
# chunking), from the bootstrap's own recorded per-replicate timings.
#
# Runs (fixture V fit once as fs.est):
#   S   sequential,    nb = 5            -> sequential per-replicate seconds now
#   M   multisession,  nb = 4*W          -> startup + makespan decomposition
#   C   multicore,     nb = 4*W          -> same, forked workers
#   Mi  multisession,  nb = 20, 20 wrk   -> identity pair across plans
#   Ci  multicore,     nb = 20, 20 wrk
# W = min(48, parallelly::availableCores()).  Every run uses the package's
# public entry point with user-level arguments; nothing is modified.
#
# Decomposition per run, all from recorded fields:
#   outer_s      = system.time() around the call
#   inner_s      = attr(results, "timing")$total_minutes * 60  (bootstrap_results()
#                  only: worker export/package load + the replicate loop)
#   spawn_s      = outer_s - inner_s                    (plan setup, main-process fits, Ystar)
#   sum_rep_s    = sum(tmins_iteration) * 60            (total replicate CPU-seconds)
#   contig_s     = max over contiguous chunks (parallel::splitIndices(nb, W)) of the
#                  chunk's replicate seconds  -- the static-chunking makespan
#   balanced_s   = sum_rep_s / W
#   overhead_s   = inner_s - contig_s  (worker export + package load + globals discovery + collection)

source("stage0_common.R")
fx <- stage0_fixture()
fs <- run_fixture_call(fx)
stopifnot(!is.null(fs$sg.harm))
kv(boot_fixture_sg = paste(fs$sg.harm, collapse = " & "))

W  <- as.integer(min(48L, parallelly::availableCores()))
NB <- 4L * W
kv(boot_workers_W = W, boot_nb_main = NB)

volatile_cols <- c("tmins_search", "tmins_iteration")
strip_volatile <- function(b) {
  r <- as.data.frame(b$results)
  r <- r[, setdiff(names(r), volatile_cols), drop = FALSE]
  list(results = r, SG_CIs = b$SG_CIs, FSsg_tab = b$FSsg_tab, Ystar_mat = b$Ystar_mat,
       H_estimates = b$H_estimates, Hc_estimates = b$Hc_estimates,
       nb_boots = b$nb_boots, original_sg = b$original_sg)
}

run_boot <- function(plan, workers, nb, tag) {
  pa <- list(plan = plan, workers = workers, show_message = FALSE)
  t0 <- proc.time()[["elapsed"]]
  b <- suppressWarnings(suppressMessages(
    forestsearch::forestsearch_bootstrap_dofuture(fs.est = fs, nb_boots = nb,
                                                  seed = 8316951L, parallel_args = pa)))
  outer_s <- proc.time()[["elapsed"]] - t0
  future::plan("sequential")
  tm <- attr(b$results, "timing")
  inner_s <- tm$total_minutes * 60
  rep_s <- b$results$tmins_iteration * 60
  rep_s <- rep_s[order(b$results$boot_id)]
  idx <- parallel::splitIndices(nb, min(workers, nb))
  contig_s <- max(vapply(idx, function(ii) sum(rep_s[ii], na.rm = TRUE), numeric(1)))
  sum_rep_s <- sum(rep_s, na.rm = TRUE)
  stats <- list(tag = tag, plan = plan, workers = workers, nb = nb,
                outer_s = outer_s, inner_s = inner_s, spawn_s = outer_s - inner_s,
                sum_rep_s = sum_rep_s, mean_rep_s = mean(rep_s, na.rm = TRUE),
                max_rep_s = max(rep_s, na.rm = TRUE),
                contig_s = contig_s, balanced_s = sum_rep_s / min(workers, nb),
                overhead_s = inner_s - contig_s,
                n_found = sum(b$results$any_found == 1L, na.rm = TRUE),
                n_rows = nrow(b$results))
  do.call(kv, setNames(stats[-1L], paste0("boot_", tag, "_", names(stats)[-1L])))
  list(stats = stats, obj = b, rep_s = rep_s)
}

runs <- list()
runs$S  <- run_boot("sequential",   1L,  5L, "S")
runs$M  <- run_boot("multisession", W,   NB, "M")
runs$C  <- run_boot("multicore",    W,   NB, "C")
W2 <- as.integer(min(20L, W))
runs$Mi <- run_boot("multisession", W2, 20L, "Mi")
runs$Ci <- run_boot("multicore",    W2, 20L, "Ci")

# identity across plans (volatile timing columns/attributes excluded)
id_pair <- identical(strip_volatile(runs$Mi$obj), strip_volatile(runs$Ci$obj))
id_results_only <- identical(strip_volatile(runs$Mi$obj)$results, strip_volatile(runs$Ci$obj)$results)
kv(identity_multicore_vs_multisession_nb20 = id_pair,
   identity_results_table_only = id_results_only)
if (!id_pair) {
  s1 <- strip_volatile(runs$Mi$obj); s2 <- strip_volatile(runs$Ci$obj)
  cat("components differing:", paste(names(s1)[!mapply(identical, s1, s2)], collapse = ", "), "\n")
}

# per-replicate slowdown under load vs sequential
kv(rep_slowdown_M = runs$M$stats$mean_rep_s / runs$S$stats$mean_rep_s,
   rep_slowdown_C = runs$C$stats$mean_rep_s / runs$S$stats$mean_rep_s,
   startup_removed_by_multicore_s = (runs$M$stats$spawn_s + runs$M$stats$overhead_s) -
                                    (runs$C$stats$spawn_s + runs$C$stats$overhead_s),
   imbalance_M_s = runs$M$stats$contig_s - runs$M$stats$balanced_s,
   imbalance_C_s = runs$C$stats$contig_s - runs$C$stats$balanced_s)

stats_tab <- do.call(rbind, lapply(runs, function(r) as.data.frame(r$stats, stringsAsFactors = FALSE)))
utils::write.csv(stats_tab, out_path("bootstrap_runs.csv"), row.names = FALSE)
print(stats_tab, row.names = FALSE)
save_result(list(stats = stats_tab, rep_s = lapply(runs, `[[`, "rep_s"), W = W, NB = NB,
                 id_pair = id_pair), "s3_bootstrap_runs")
cat("\ns3 complete\n")
```

### `s4_summary.R`

```r
# s4_summary.R -- assembles results/stage0_summary.md from the key = value
# record written by s0-s3.  Pure arithmetic on measured numbers; every figure
# in the summary is traceable to a line in results/keyvalues.txt.

source("stage0_common.R")
lines <- readLines(out_path("keyvalues.txt"))
kvs <- strsplit(lines, " = ", fixed = TRUE)
vals <- setNames(lapply(kvs, function(z) paste(z[-1L], collapse = " = ")), vapply(kvs, `[[`, "", 1L))
# later writes of the same key win
vals <- vals[!duplicated(names(vals), fromLast = TRUE)]
num <- function(k) { v <- vals[[k]]; if (is.null(v)) NA_real_ else suppressWarnings(as.numeric(strsplit(v, ", ")[[1L]][1L])) }
chr <- function(k) { v <- vals[[k]]; if (is.null(v)) "<missing>" else v }
pct <- function(x) if (is.na(x)) "n/a" else sprintf("%.1f%%", 100 * x)
s   <- function(x, d = 2) if (is.na(x)) "n/a" else formatC(x, format = "f", digits = d)

BAR <- 0.30

# ---- A ------------------------------------------------------------------------
rep_s     <- num("boot_S_mean_rep_s")          # sequential replicate seconds, measured now
call_s    <- num("call_s_median")              # single sequential call, measured now
mach_s    <- num("machinery_s_per_call")       # direct timers (s2)
mach_rp   <- num("rprof_formula_machinery_s")  # Rprof corroboration (s1)
A_share_rep  <- mach_s / rep_s
A_share_call <- mach_s / call_s
A_clears  <- isTRUE(A_share_rep >= BAR)

# ---- C / D --------------------------------------------------------------------
W <- num("boot_workers_W"); NB <- num("boot_nb_main")
g <- function(tag, f) num(sprintf("boot_%s_%s", tag, f))
startup_M <- g("M", "spawn_s") + g("M", "overhead_s")
startup_C <- g("C", "spawn_s") + g("C", "overhead_s")
# Projected application-shape wall at B replicates:
#   startup (plan-specific) + B * rep / W + static-chunk imbalance (plan-specific,
#   as measured at nb = 4W).  The replicate term is COMMON to both plans (mean of
#   the two loaded per-replicate means) so that the C comparison isolates startup
#   and is not swamped by run-to-run noise in the replicate mean.
rep_common <- mean(c(g("M", "mean_rep_s"), g("C", "mean_rep_s")))
proj <- function(tag, B, rep_mult = 1) {
  (g(tag, "spawn_s") + g(tag, "overhead_s")) + B * rep_common * rep_mult / W +
    (g(tag, "contig_s") - g(tag, "balanced_s"))
}
B <- 1000
wall_M_1000 <- proj("M", B); wall_C_1000 <- proj("C", B)
C_share_1000 <- 1 - wall_C_1000 / wall_M_1000
C_share_meas <- 1 - g("C", "outer_s") / g("M", "outer_s")
D_share_M <- (g("M", "contig_s") - g("M", "balanced_s")) / g("M", "outer_s")
# A inside the application shape: replicate seconds shrink by the machinery share
A_rep_mult <- 1 - A_share_rep
wall_M_1000_A  <- proj("M", B, A_rep_mult)
wall_C_1000_A  <- proj("C", B, A_rep_mult)

md <- c(
"# Stage 0 -- bootstrap efficiency measurements (CRAN forestsearch 0.4.0)",
"",
sprintf("Installed: forestsearch %s (%s; repository %s) · survival %s · future %s · %s · %s/%s · cores %s · multicore %s",
        chr("installed_version"), chr("installed_packaged"), chr("installed_repository"),
        chr("survival_version"), chr("future_version"), chr("r_version"), chr("sysname"), chr("machine"),
        chr("available_cores"), chr("supports_multicore")),
sprintf("Repo: %s @ %s (%s); files matching CRAN 0.4.0: %s, differing: %s",
        chr("repo_root"), chr("repo_head"), chr("repo_branch"),
        chr("repo_files_matching_cran"), chr("repo_files_differing")),
sprintf("Fixture: %s (n = %s)", chr("fixture"), chr("fixture_n")),
"",
"## 1. The single sequential call (fixture V)",
"",
"| quantity | value |",
"|---|---|",
sprintf("| call wall-clock, median of 3 (s) | %s (reps %s) |", s(call_s, 3), chr("call_s_rep")),
sprintf("| candidates evaluated in the search | %s |", chr("search_n_evaluated")),
sprintf("| candidates reaching the Cox fit | %s |", chr("search_n_reached_cox")),
sprintf("| candidates passing the effect screen | %s |", chr("search_n_passed_effect")),
sprintf("| consistency stage: candidates total / evaluated / passed | %s / %s / %s (early stop: %s) |",
        chr("consistency_n_total"), chr("consistency_n_evaluated"), chr("consistency_n_passed"), chr("consistency_early_stop")),
sprintf("| identified subgroup | %s |", chr("identified_subgroup")),
"",
"Rprof, median of 3 profiles (seconds of total time; see rprof_watched_median.csv for the full list):",
"",
"| function | total s |",
"|---|---|",
sprintf("| coxph (formula interface, all callers) | %s |", s(num("rprof_coxph_total_s"))),
sprintf("| coxph.fit (+ agreg.fit) -- the fit proper | %s |", s(num("rprof_coxph_fit_total_s"))),
sprintf("| formula-interface machinery = coxph - coxph.fit | %s (%s of coxph) |", s(mach_rp), pct(num("rprof_machinery_share_of_coxph"))),
sprintf("| of which concordancefit | %s |", s(num("rprof_concordancefit_s"))),
sprintf("| search_combinations_parallel (whole search stage) | %s |", s(num("rprof_search_stage_s"))),
sprintf("| fit_cox_for_subgroup (search-stage fits) | %s |", s(num("rprof_fit_cox_for_subgroup_s"))),
sprintf("| survfit (medians on effect-screen passers) | %s |", s(num("rprof_survfit_s"))),
sprintf("| evaluate_consistency_twostage (whole consistency evaluation) | %s |", s(num("rprof_consistency_stage_s"))),
sprintf("| .consistency_cox_pieces / residuals.coxph | %s / %s |", s(num("rprof_consistency_cox_pieces_s")), s(num("rprof_residuals_coxph_s"))),
sprintf("| coxph time outside the search stage (Rprof coxph total - direct-timed search fits; candidate B's territory) | %s |",
        s(num("rprof_coxph_total_s") - num("time_formula_s_median"))),
"",
"## 2. Candidate A -- lean search-stage Cox fit (coxph.fit instead of formula coxph)",
"",
"| quantity | value |",
"|---|---|",
sprintf("| captured search-stage fits (= candidates reaching the fit) | %s (matches filter counts: %s; adjusted: %s) |",
        chr("capture_n_fits"), chr("capture_matches_filter_counts"), chr("capture_n_adjusted")),
sprintf("| CI rows identical, lean vs formula (identical()) | %s of %s |", chr("ident_ci_row_identical"), chr("ident_n_candidates")),
sprintf("| CI values identical ignoring dimnames | %s |", chr("ident_ci_values_identical")),
sprintf("| failure sets agree (try-errors) | %s (formula %s, lean %s) |", chr("ident_error_sets_agree"), chr("ident_errors_formula"), chr("ident_errors_lean")),
sprintf("| warning sets agree | %s (candidates warned: %s; identical among them: %s) |", chr("ident_warn_sets_agree"), chr("ident_n_warned"), chr("ident_identical_among_warned")),
sprintf("| formula path, all captured fits, median of 5 (s) | %s (%s ms/fit) |", s(num("time_formula_s_median"), 3), s(num("time_ms_per_fit_formula"), 3)),
sprintf("| lean path, all captured fits, median of 5 (s) | %s (%s ms/fit) |", s(num("time_lean_s_median"), 3), s(num("time_ms_per_fit_lean"), 3)),
sprintf("| **machinery removed per call (s)** | **%s** (direct timers) vs %s (Rprof coxph - coxph.fit, all callers) |", s(mach_s, 3), s(mach_rp, 3)),
sprintf("| machinery share of the formula path | %s |", pct(num("machinery_share_of_formula"))),
sprintf("| sequential replicate now, mean of 5 (s) | %s |", s(rep_s, 3)),
sprintf("| **A: share of a sequential replicate** | **%s** (of the single call: %s) |", pct(A_share_rep), pct(A_share_call)),
sprintf("| **A clears the %s bar on the per-replicate (campaign) shape?** | **%s** |", pct(BAR), if (A_clears) "YES" else "NO"),
"",
"Captured coxph.fit() conventions (the lean recipe, from the installed survival):",
"",
sprintf("- coxph.fit formals: `%s`", chr("conv_coxph_fit_formals")),
sprintf("- x: storage `%s`, colnames `%s`, values = treatment column: %s", chr("conv_x_storage"), chr("conv_x_colnames"), chr("conv_x_values_equal_treat")),
sprintf("- y: type `%s`; equals Surv(Y,E): %s; equals aeqSurv(Surv(Y,E)): %s; control$timefix = %s",
        chr("conv_y_type"), chr("conv_y_values_equal_raw"), chr("conv_y_values_equal_aeqSurv"), chr("conv_timefix")),
sprintf("- strata NULL: %s · offset: %s · init: %s · control = coxph.control(): %s · weights NULL: %s · method: %s · rownames: %s · resid: %s · nocenter: %s",
        chr("conv_strata_is_null"), chr("conv_offset_desc"), chr("conv_init_desc"), chr("conv_control_is_default"),
        chr("conv_weights_is_null"), chr("conv_method"), chr("conv_rownames_desc"), chr("conv_resid"), chr("conv_nocenter")),
"",
"## 3. Candidates C and D -- application shape (bootstrap at W workers)",
"",
sprintf("W = %s workers, nb = %s for the M/C runs; identity pair at nb = 20.", chr("boot_workers_W"), chr("boot_nb_main")),
"",
"| run | plan | nb | outer s | inner s | spawn s | overhead s | sum rep s | mean rep s | contiguous makespan s | balanced s | imbalance s |",
"|---|---|---|---|---|---|---|---|---|---|---|---|",
vapply(c("S", "M", "C", "Mi", "Ci"), function(t) sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
        t, chr(sprintf("boot_%s_plan", t)), chr(sprintf("boot_%s_nb", t)), s(g(t, "outer_s"), 1), s(g(t, "inner_s"), 1),
        s(g(t, "spawn_s"), 1), s(g(t, "overhead_s"), 1), s(g(t, "sum_rep_s"), 1), s(g(t, "mean_rep_s"), 2),
        s(g(t, "contig_s"), 1), s(g(t, "balanced_s"), 1), s(g(t, "contig_s") - g(t, "balanced_s"), 1)), character(1)),
"",
"| quantity | value |",
"|---|---|",
sprintf("| identity across plans, nb = 20 (all non-volatile components) | %s (results table only: %s) |", chr("identity_multicore_vs_multisession_nb20"), chr("identity_results_table_only")),
sprintf("| per-replicate slowdown under load (M / sequential, C / sequential) | %s / %s |", s(num("rep_slowdown_M"), 2), s(num("rep_slowdown_C"), 2)),
sprintf("| startup (spawn + load): multisession vs multicore (s) | %s vs %s; removed: %s |", s(startup_M, 1), s(startup_C, 1), s(num("startup_removed_by_multicore_s"), 1)),
sprintf("| **C: measured wall reduction at nb = %s** | **%s** |", chr("boot_nb_main"), pct(C_share_meas)),
sprintf("| **C: projected wall reduction at B = 1000** (startup + B·rep/W + imbalance) | **%s** (M %s s -> C %s s) |", pct(C_share_1000), s(wall_M_1000, 0), s(wall_C_1000, 0)),
sprintf("| **D: static-chunk imbalance, share of the M run's wall** | **%s** (%s s of %s s) |", pct(D_share_M), s(g("M", "contig_s") - g("M", "balanced_s"), 1), s(g("M", "outer_s"), 1)),
sprintf("| A inside the application shape at B = 1000: multisession %s s -> %s s; multicore %s s -> %s s | A+C combined: %s |",
        s(wall_M_1000, 0), s(wall_M_1000_A, 0), s(wall_C_1000, 0), s(wall_C_1000_A, 0), pct(1 - wall_C_1000_A / wall_M_1000)),
"",
"## 4. Reading",
"",
sprintf("- Per-replicate (campaign) shape: A removes %s s of a %s s replicate = %s; bar %s: %s.",
        s(mach_s, 2), s(rep_s, 2), pct(A_share_rep), pct(BAR), if (A_clears) "cleared" else "not cleared"),
sprintf("- Bit-identity of the lean fit: %s of %s CI rows identical; failure sets agree: %s.",
        chr("ident_ci_row_identical"), chr("ident_n_candidates"), chr("ident_error_sets_agree")),
sprintf("- Application shape at B = 1000 (projection from measured components): C alone %s; D (imbalance) %s of the measured M wall; A+C %s.",
        pct(C_share_1000), pct(D_share_M), pct(1 - wall_C_1000_A / wall_M_1000)),
"- Projections carry the per-replicate slowdown under load as measured and the imbalance as measured at nb = 4W; they are not a substitute for a timed B = 1000 run.",
""
)
writeLines(md, out_path("stage0_summary.md"))
cat(md, sep = "\n")
cat("\ns4 complete\n")
```
