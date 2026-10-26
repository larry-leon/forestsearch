# CC task — bootstrap efficiency Stage 0b: the GLM paths and the consistency stage, measured, with the modification map — no package code modified

**File:** `dev/tasks/cc_task_bootstrap_stage0b_2026-10-26.md` · **Issued:** 2026-10-26 by chat, commissioned by Larry ("a comprehensive measurement to gain a clear quantification of the efficiency gains and where the package functions will be modified — do NOT modify any R package functions")
**Repo:** `forestsearch` · **Branch:** `master` · **Transport:** `~/Downloads` → `dev/tasks/` (this document, committed first and alone) and `~/Downloads` → `dev/glm-continuous-sims/stage0b_glm_2026-10-26/` (the script bundle, §2). Filenames may arrive with hyphens stripped; locate by pattern.
**Predecessor:** Stage 0 (`dev/tasks/cc_task_bootstrap_stage0_2026-10-26.md`, `dev/glm-continuous-sims/REPORT_bootstrap_stage0_2026-10-26.md`) measured the survival search-stage fit (candidate A) and the outer plan (C, D). This task extends the same method to everything else and produces the map of where the package would change.

**Audit basis.** CRAN `forestsearch` 0.4.0 (Date/Publication 2026-10-10; `Packaged: 2026-10-10 04:12:05 UTC`). Every line reference below is to the 0.4.0 sources as published; `b0` gates the installed version at 0.4.0 and compares the working tree's files under study to the CRAN MD5 sums.

**What this task answers.**
- For each outcome type (survival, continuous/MD, binary/OR) and each stage of a `forestsearch()` call (search-stage per-candidate fit; consistency-stage per-candidate pieces; consistency-stage member extraction): the seconds a lean equivalent would remove per call, as a share of a measured sequential bootstrap replicate, against the 30% bar; and whether the lean equivalent is bit-identical on every candidate the package evaluates.
- The modification map: which functions, at which 0.4.0 lines, each candidate would change, and with what; each row carries its measured status.

**What this task decides: nothing.** Whether to commission any change is Larry's decision afterwards.

---

## ⚠ CATEGORY

**ABSOLUTELY NO PACKAGE CODE IS MODIFIED.** No edit to any file under `R/`, `tests/`, `man/`, `vignettes/`, `inst/`, `src/`, `DESCRIPTION`, `NAMESPACE`, `NEWS.md`, or any driver, application document or payload. No `devtools::install()`, no `load_all()`, no `R CMD check`, no test suite, no vignette build, no render. The scripts call the installed package through its public entry points; every "package path" they time is the installed package's own function, unmodified, run on recorded inputs.

As in Stage 0, `b2` uses **`base::trace()` in the running session only** — on `forestsearch::evaluate_consistency_twostage`, `forestsearch::evaluate_subgroup_consistency`, `forestsearch:::.consistency_glm_pieces` / `.consistency_cox_pieces`, `forestsearch:::fit_glm_for_subgroup` / `fit_glm_for_subgroup_fast`, and `survival::coxph.fit` — to record their arguments during one call per fixture, and `untrace()`s before any timing. Nothing on disk is touched. Machine-checkable post-condition in §6.

**Writes:** this document under `dev/tasks/`; the bundle's eight files (two of them the fixture files, as transplanted), their `results/` and `logs/`, under `dev/glm-continuous-sims/stage0b_glm_2026-10-26/`; the report `dev/glm-continuous-sims/REPORT_bootstrap_stage0b_2026-10-26.md`.

**Compute (approved with the commission; stated as its own line).** Per fixture: 1 warm-up + 3 timed + 3 profiled + 1 traced sequential `forestsearch()` calls; one sequential `nb_boots = 5` bootstrap; and the timing loops of `b2` (each candidate's package function and lean equivalent, 3 alternating reps over all recorded candidates). Three fixtures. Prior measured figures: a sequential gbsg call is 3.2 s (Stage 0); the ACTG175 continuous application call was measured at ≈8 s sequential in the 0.3.2 record (reprofile task §2); no prior figure for the binary application. Hard caps: `b1` 60 min, `b2` 60 min (driver). No parallel bootstrap runs, no B = 1000, no campaigns.

**Unattended.** Gates stop the task; never ask, never work around. Provenance is dirt-tolerant (pre-existing untracked files are never staged and never a stop unless under this task's own paths).

**Run from a plain terminal.**

---

## 1. Provenance and first commit — GATE

Standard block: `git rev-parse HEAD`, branch (expect **`master`**), `git status --porcelain` (verbatim), installed `packageVersion("forestsearch")` (expect **0.4.0**), `packageVersion("survival")`, `packageVersion("data.table")`, `R.version.string`; the installed DESCRIPTION's `Packaged:` and `Repository:` fields. Copy this document to `dev/tasks/cc_task_bootstrap_stage0b_2026-10-26.md`, `git add` it by path, commit it alone. Report both filenames.

## 2. Place the bundle, then transplant the two GLM fixtures — GATE

1. `sha256sum ~/Downloads/stage0b_glm_2026-10-26.tgz` must print `5f61d5182116b0e1668a9687025dcccb5bea5c6c2d0cea09009816478aab5a3c`. If not (or absent), save the eight sources from Appendix A verbatim into `dev/glm-continuous-sims/stage0b_glm_2026-10-26/` and continue with step 3.
2. `tar -xzf ~/Downloads/stage0b_glm_2026-10-26.tgz -C dev/glm-continuous-sims/` → creates `dev/glm-continuous-sims/stage0b_glm_2026-10-26/` with exactly eight files.
3. In that directory, `sha256sum *` must print exactly:

```
98f5bc1af8bb91e556cc60477ed760a2649d062bd95aaa71c0d2f199b032903b  b0_provenance.R
08001177c58c43bb9c2223ccda64bfbcac2bf0f0acd29e791c9742c06212c8e6  b1_profiles.R
453573052f8e73d19902749da7999398d1539e8776d54f0317db662cb06ecfea  b2_lean.R
fcfa83c8775a4c74c9abc5010f7676fc36e9956015aa0c766f6404549696419d  b3_summary.R
08daeb07132eb05abf1b65d53f4a2584a9cbe3d00214f822c3d042336a2ef423  fixture_binary.R
16044dfb49a9e10d5532f66552a21ffd4e9314d4dacb58b99795b25b749166a3  fixture_continuous.R
024ae519c9d80081c7a21e0672b2d6521f594337e67b65246cd7e7996d4438b4  run_stage0b.sh
9429af0edb65b30c28d50991f43b0a497454dcd9e228ae42c14763fb91b9b101  stage0b_common.R
```

   *GATE:* any mismatch after both routes is a stop.

4. **Transplant the fixtures (the only authoring in this task, and it is a copy).** `fixture_continuous.R` and `fixture_binary.R` ship with chat's APPROXIMATE ACTG175 configurations as defaults. Replace each `args` block with the committed application's `forestsearch()` parameter block, verbatim, found in this order of preference, recording which:
   - the ACTG175 continuous and binary analysis documents under `quarto/applications/actg175/` in this repo (the continuous one is the `y_decline` / MD analysis with `adverse_outcome = TRUE`, `effect.threshold = 10`, `consistency.threshold = 5`, twelve confounders; the binary one is the OR analysis);
   - failing that, the corresponding documents under `~/GitHub/fs-glms-interpretable/quarto/actg175/` — **read-only**, nothing written there;
   - failing both for one outcome type, keep that default and say so.

   Apply exactly these overrides to the transplanted block, nothing else: `mr_inference = FALSE` (drop `mr_inference_args`), `details = FALSE`, `quiet = TRUE`, `plot.sg = FALSE`, `plot.grf = FALSE`, `parallel_args = list(plan = "sequential")`. Keep `use_grf`, `use_lasso`, `vi.grf.min`, `use_twostage`, `consistency_method` and every threshold exactly as the document has them — the bootstrap replicate inherits them. The data-preparation lines (arm subset, treatment coding, outcome construction, derived covariates) are transplanted too; `speff2trial::ACTG175` is the source (`b0` gates that the package is installed). Set each fixture's `source` string to the document path and chunk name. The report quotes both transplanted blocks and the diff against the shipped defaults.

   The fixture files are the only files in the bundle whose checksum is expected to change; the other seven must still match the manifest after step 4 (`sha256sum` them again and print the result in the report).

## 3. Run

```
cd dev/glm-continuous-sims/stage0b_glm_2026-10-26 && bash run_stage0b.sh
```

`b0` → `b3` in order, each under `timeout` (5 / 60 / 60 / 5 min), logs in `logs/`, results in `results/`. Gates inside the scripts:

| step | gate | exit | meaning |
|---|---|---|---|
| b0 | installed `forestsearch` is 0.4.0 | 2 | wrong version — this audit is of CRAN 0.4.0 |
| b0 | `speff2trial` installed | 3 | ACTG175 unavailable |
| b1 | survival candidate flow equals CRAN 0.4.0's (1,410 reach the fit, 121 pass, 120 enter consistency, `{er <= 0} & {pgr <= 26}`) | 5 | the installed package is not computing what 0.4.0 computes |
| b2 | a consistency-stage call captured for every fixture | 4 | the evaluator did not run in-process; nothing to measure |

A fixture that identifies no subgroup is not a stop: `b1` records it and the bootstrap denominator for that fixture is `NA`; the candidate shares for it are then reported against the single-call time (say so in the report).

**Script deviations.** The scripts ran end to end in chat's dry run (CRAN 0.4.0 from CRAN, R 4.3.3, survival 3.5.8, 2 cores, Linux, the shipped approximate fixtures). If a step fails for a reason that is a script defect rather than a gate, fix the script with the smallest edit that restores the stated measurement, re-run from that step, and list every edit verbatim in the report's "Script deviations". Never change what is measured or how identity is judged; if that would be needed, stop and report.

## 4. What each step measures (the script headers say the same)

- **b0** provenance, gates, repo files under study vs CRAN 0.4.0 MD5, fixture sources.
- **b1** per fixture: warm-up, candidate-flow counts, 3 timed calls, 3 `Rprof` profiles bucketed by stage (search stage / per-candidate fit / consistency stage / pieces / extraction remainder / dedup / GRF / LASSO / medians), and a sequential `nb = 5` bootstrap for the replicate denominator.
- **b2** per fixture, from one traced call: the consistency evaluator's arguments (shared objects once, `m` per candidate), the pieces inputs per evaluated candidate (needed columns only), the search-stage fit inputs per candidate. Then, per candidate, package function vs lean equivalent with `identical()` and 3 alternating timing reps:
  - **B / E / G** — consistency pieces (`.consistency_cox_pieces` / `.consistency_glm_pieces`): lean = `coxph.fit` ×2 + `residuals.coxph` on a constructed object (survival), `lm.fit` + hat values from the fit's QR (MD), `glm.fit` + hat values from the fit's QR (OR); `.dfbeta_glm()`'s arithmetic line for line. Identity on `beta_hat`, `sigma_D` and the closed-form rate.
  - **F** — binary search-stage fit (`fit_glm_for_subgroup` → `.estimate_or`): lean = vectors + `glm.fit` + SE from the fit's QR, with the package's own `.fs_existence_reason()` / `.fs_binary_cells()` for the existence check. Identity on the returned `hr/lower/upper` list and on the NULL / non-estimable pattern.
  - **H** — the whole consistency evaluator (`evaluate_consistency_twostage`, or `evaluate_subgroup_consistency` if the fixture has `use_twostage = FALSE`): package function per candidate vs a lean evaluator (precomputed logical column cache + index matrix instead of `as.data.table(df)` + sequential row filters per candidate; the needed vectors only; then the lean pieces), and an "H alone" variant (lean extraction, package pieces). Identity on `Pcons` (rounded as the package rounds) and on the pass/fail pattern.
  - Survival's search-stage fit (A) is Stage 0's and is quoted, not re-measured; the continuous search-stage fit is already the 0.3.3 fast path and is timed as a bucket only.
- **b3** assembles `results/stage0b_summary.md`: per-fixture call/replicate/bucket tables; the candidate table (package s, lean s, machinery s, share of replicate, identity, bar); the combined per-fixture share; the modification map with measured status.

## 5. Report — `dev/glm-continuous-sims/REPORT_bootstrap_stage0b_2026-10-26.md`

1. **Audit basis** (first lines) and provenance; the MD5 table; fixture sources with both transplanted parameter blocks quoted verbatim and the diff against the shipped defaults; whether each GLM fixture's identified subgroup matches the document's committed payload (read-only; a difference is a finding, not a stop).
2. `results/stage0b_summary.md` **transplanted verbatim**.
3. **Reading, numbers only**: per fixture, which candidates clear 30% alone, which only in combination, and the identity status of each; where the consistency stage's share comes from (number of candidates evaluated under full evaluation vs early stop). State plainly which figures are for the application configurations and which for chat's defaults, if any default was kept.
4. **Modification map** — the table from the summary, plus for each row: category (re-implements the same computation with a bit-identity target / moves code / changes behaviour — expected "same computation" for all; say if any measured identity failed), exported vs internal, whether any signature changes (expected none).
5. **Script deviations** (verbatim edits, or "none").
6. **OPEN ITEMS**.
7. **Ten-line verdict.** No recommendation to implement anything.

Stop-point figures must be quoted verbatim from `results/keyvalues.txt` / `results/stage0b_summary.md`.

## 6. Commit and post-conditions — GATE

Commit by explicit path, in order, no push, no fetch: (1) the task document alone (§1); (2) `dev/glm-continuous-sims/stage0b_glm_2026-10-26/` — the eight sources as they stand after §2.4, `results/` (all files), `logs/`; and the report. If `dev/glm-continuous-sims/current_status.md` exists, update it per the standing rule as the final commit; it did not exist at Stage 0 — if still absent, note it and do not create one.

Before each commit `git diff --cached --name-only` must list only paths under `dev/`. After the final commit, print in the report's closing block:

```
git diff --name-only <HEAD recorded in §1> HEAD | grep -v '^dev/'      # must print nothing
git status --porcelain | grep -v '^?? ' | grep -v '^.. dev/'            # must print nothing
git status --porcelain | grep '^?? ' | grep -v ' dev/'                   # only pre-existing dirt, listed in §1
```

The first line is the machine check of "no package code modified". If it prints anything, do not fix it — stop and report what it prints.

## 7. Out of scope

- No `R/` edit, no estimator change, no driver, application, document or payload change. No install, check, tests, vignette build, render, push or fetch.
- No parallel bootstrap runs, no B = 1000, no campaign cells. No count/IRR, RR, RD, adjusted, propensity-weighted, DINA or GRF measurements (the summary says why).
- No re-running of Stage 0's measurements.

---

## Appendix A — script sources, verbatim (fallback for §2 only)

Save each block to the named file under `dev/glm-continuous-sims/stage0b_glm_2026-10-26/`, byte for byte, then re-check §2 step 3. The tarball and these blocks are the same bytes.

### `run_stage0b.sh`

```bash
#!/usr/bin/env bash
# run_stage0b.sh -- driver for the bootstrap-efficiency Stage-0b measurement
# (GLM paths and the consistency stage).  Measurement only: no package file is
# changed.  Run from a plain terminal:
#     cd dev/glm-continuous-sims/stage0b_glm_2026-10-26 && bash run_stage0b.sh
# Each step has a hard time cap; a non-zero exit stops the driver.
set -uo pipefail
cd "$(dirname "$0")"
export STAGE0B_OUT=results
mkdir -p results logs
run() {
  local name=$1 cap=$2; shift 2
  echo "== $name (cap $cap)"
  timeout -k 15s "$cap" "$@" > "logs/$name.log" 2>&1
  local rc=$?
  tail -n 3 "logs/$name.log" | sed 's/^/   /'
  if [ $rc -ne 0 ]; then echo "STOP: $name exited $rc -- see logs/$name.log"; exit $rc; fi
}
run b0_provenance 5m  Rscript b0_provenance.R
run b1_profiles   60m Rscript b1_profiles.R
run b2_lean       60m Rscript b2_lean.R
run b3_summary    5m  Rscript b3_summary.R
echo "== all steps complete; summary at results/stage0b_summary.md"
```

### `stage0b_common.R`

```r
# stage0b_common.R -- shared setup for the bootstrap-efficiency Stage-0b
# measurement (task cc_task_bootstrap_stage0b_2026-10-26.md): the GLM paths
# and the consistency stage, all three outcome types.
#
# MEASUREMENT ONLY.  Nothing here edits any file of the forestsearch package,
# on disk or in the installed library.  b2 uses base::trace() on package
# functions IN THE RUNNING SESSION to record their arguments, and untrace()s
# them before any timing.  Every "formula path" timed below is the installed
# package's own function, called unmodified on the recorded arguments.
#
# Audit basis: CRAN forestsearch 0.4.0 (Date/Publication 2026-10-10).

suppressPackageStartupMessages({
  library(forestsearch)
  library(survival)
  library(future)
  library(doFuture)
  library(data.table)
})

OUT <- Sys.getenv("STAGE0B_OUT", "results")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
out_path <- function(...) file.path(OUT, ...)
save_result <- function(x, name) { saveRDS(x, out_path(paste0(name, ".rds"))); invisible(x) }

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

# ---- fixtures ----------------------------------------------------------------
# Survival: fixture V of Stage 0 (the package's own vignette gbsg configuration,
# mr_inference = FALSE, sequential).  Continuous and binary: the committed
# ACTG175 application configurations, transplanted by CC into
# fixture_continuous.R / fixture_binary.R (see the task document, section 2).
# Each fixture function returns list(label, df, args, source).
fixture_survival <- function() {
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
  list(label = "survival: vignette gbsg config (Stage-0 fixture V)", df = df.analysis, args = args,
       source = "inst/doc/survival-post-selection.qmd (CRAN 0.4.0), mr_inference = FALSE")
}

source("fixture_continuous.R")   # defines fixture_continuous()
source("fixture_binary.R")       # defines fixture_binary()

FIXTURES <- list(survival = fixture_survival, continuous = fixture_continuous, binary = fixture_binary)

run_fixture_call <- function(fx) {
  suppressWarnings(suppressMessages(do.call(forestsearch::forestsearch, fx$args)))
}

candidate_counts <- function(fs) {
  fc <- fs$find.grps$filter_counts
  gc <- fs$grp.consistency
  list(
    n_evaluated_search      = fc$n_evaluated,
    n_reached_fit           = fc$n_passed_sample_size,
    n_fit_ok                = fc$n_passed_cox,
    n_passed_effect         = fc$n_passed_hr,
    n_consistency_total     = gc$n_candidates_total,
    n_consistency_evaluated = gc$n_candidates_evaluated,
    n_consistency_passed    = gc$n_passed,
    early_stop_triggered    = isTRUE(gc$early_stop_triggered),
    sg_harm                 = paste(fs$sg.harm, collapse = " & ")
  )
}

# capture environment used by the tracers (lives in .GlobalEnv so a tracer
# evaluated inside a package function's frame can reach it)
new_cap <- function() {
  e <- new.env()
  e$n <- 0L
  e$items <- vector("list", 65536L)
  e$push <- function(item) {
    e$n <- e$n + 1L
    if (e$n > length(e$items)) length(e$items) <- 2L * length(e$items)
    e$items[[e$n]] <- item
  }
  e
}
cap_items <- function(e) e$items[seq_len(e$n)]

# median of timing reps, alternating order to avoid order bias
time_pair <- function(f_a, f_b, reps = 3L) {
  ta <- tb <- numeric(reps)
  for (r in seq_len(reps)) {
    if (r %% 2L == 1L) { ta[r] <- system.time(f_a())[["elapsed"]]; tb[r] <- system.time(f_b())[["elapsed"]] }
    else               { tb[r] <- system.time(f_b())[["elapsed"]]; ta[r] <- system.time(f_a())[["elapsed"]] }
  }
  list(a = ta, b = tb, a_med = median(ta), b_med = median(tb))
}
```

### `fixture_continuous.R`

```r
# fixture_continuous.R -- the continuous (MD) fixture for Stage 0b.
#
# CC: replace the `args` block below with the committed ACTG175 continuous
# application's forestsearch() parameter block, transplanted VERBATIM from the
# document named in the task (quarto/applications/actg175/, the continuous
# analysis), with exactly these overrides for the per-replicate shape:
#   mr_inference = FALSE, details = FALSE, quiet = TRUE, plot.sg = FALSE,
#   plot.grf = FALSE, parallel_args = list(plan = "sequential"),
# and drop mr_inference_args.  Keep use_grf / use_lasso / vi.grf.min exactly as
# the document has them (the bootstrap re-runs them per replicate).  Set
# `source` to the document path and the chunk name.  Record the diff against
# this default in the report.
#
# DEFAULT (chat's approximation, used only if no committed configuration is
# located): ACTG175 arms 1 vs 3, outcome y_decline = cd40 - cd420 (positive =
# CD4 decline, harm-oriented), the twelve confounders of the application
# (six continuous, six binary, no str2), ten quantile cuts on each continuous
# confounder, c1 = 10, c2 = 5, p* = 0.80, effMaxSG, resample consistency.

fixture_continuous <- function() {
  if (!requireNamespace("speff2trial", quietly = TRUE))
    stop("speff2trial is required for the ACTG175 fixtures")
  d <- speff2trial::ACTG175
  d <- d[d$arms %in% c(1, 3), ]
  d$treat <- 1L - as.integer(d$arms == 1)
  d$id <- seq_len(nrow(d))
  d$y_decline <- d$cd40 - d$cd420
  cont <- c("age", "wtkg", "karnof", "preanti", "cd40", "cd80")
  bin  <- c("hemo", "homo", "drugs", "race", "gender", "symptom")
  args <- list(
    df.analysis = d,
    outcome.name = "y_decline", treat.name = "treat", id.name = "id",
    confounders.name = c(cont, bin),
    outcome_type = "continuous", effect_measure = "MD", adverse_outcome = TRUE,
    is.RCT = TRUE, seedit = 8316951,
    use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE,
    conf.cont_jcuts = as.list(setNames(rep(10L, length(cont)), cont)), collapse_cuts = TRUE,
    subgroup_method = "consistency", sg_focus = "effMaxSG",
    selection_rule = "neighborhood", effect_neighborhood = 0.20,
    consistency_method = "resample",
    n.min = 60, maxk = 2, max.minutes = 3,
    hr.threshold = 10, hr.consistency = 5, pconsistency.threshold = 0.80,
    mr_inference = FALSE,
    details = FALSE, quiet = TRUE, plot.sg = FALSE, plot.grf = FALSE,
    parallel_args = list(plan = "sequential")
  )
  list(label = "continuous: ACTG175 MD (APPROXIMATE DEFAULT -- replace with the committed application block)",
       df = d, args = args,
       source = "APPROXIMATE DEFAULT (chat), not a committed configuration")
}
```

### `fixture_binary.R`

```r
# fixture_binary.R -- the binary (OR) fixture for Stage 0b.
#
# CC: replace the `args` block below with the committed ACTG175 binary
# application's forestsearch() parameter block, transplanted VERBATIM from the
# document named in the task (quarto/applications/actg175/, the binary/OR
# analysis), with exactly these overrides for the per-replicate shape:
#   mr_inference = FALSE, details = FALSE, quiet = TRUE, plot.sg = FALSE,
#   plot.grf = FALSE, parallel_args = list(plan = "sequential"),
# and drop mr_inference_args.  Keep use_grf / use_lasso / vi.grf.min exactly as
# the document has them.  Set `source` to the document path and chunk name.
# Record the diff against this default in the report.
#
# DEFAULT (chat's approximation, used only if no committed configuration is
# located): ACTG175 arms 1 vs 3, outcome y_bin = 1 if CD4 did not improve
# (cd420 <= cd40; adverse), the same twelve confounders and cuts as the
# continuous default, c1 = 1.0 (OR), c2 = 0.8, p* = 0.80, effMaxSG,
# resample consistency, d0.min = d1.min = 10.

fixture_binary <- function() {
  if (!requireNamespace("speff2trial", quietly = TRUE))
    stop("speff2trial is required for the ACTG175 fixtures")
  d <- speff2trial::ACTG175
  d <- d[d$arms %in% c(1, 3), ]
  d$treat <- 1L - as.integer(d$arms == 1)
  d$id <- seq_len(nrow(d))
  d$y_bin <- as.integer(d$cd420 <= d$cd40)
  cont <- c("age", "wtkg", "karnof", "preanti", "cd40", "cd80")
  bin  <- c("hemo", "homo", "drugs", "race", "gender", "symptom")
  args <- list(
    df.analysis = d,
    outcome.name = "y_bin", treat.name = "treat", id.name = "id",
    confounders.name = c(cont, bin),
    outcome_type = "binary", effect_measure = "OR", adverse_outcome = TRUE,
    is.RCT = TRUE, seedit = 8316951,
    use_lasso = FALSE, use_grf = FALSE, use_dina = FALSE,
    conf.cont_jcuts = as.list(setNames(rep(10L, length(cont)), cont)), collapse_cuts = TRUE,
    subgroup_method = "consistency", sg_focus = "effMaxSG",
    selection_rule = "neighborhood", effect_neighborhood = 0.20,
    consistency_method = "resample",
    n.min = 60, d0.min = 10, d1.min = 10, maxk = 2, max.minutes = 3,
    hr.threshold = 1.0, hr.consistency = 0.8, pconsistency.threshold = 0.80,
    mr_inference = FALSE,
    details = FALSE, quiet = TRUE, plot.sg = FALSE, plot.grf = FALSE,
    parallel_args = list(plan = "sequential")
  )
  list(label = "binary: ACTG175 OR (APPROXIMATE DEFAULT -- replace with the committed application block)",
       df = d, args = args,
       source = "APPROXIMATE DEFAULT (chat), not a committed configuration")
}
```

### `b0_provenance.R`

```r
# b0_provenance.R -- provenance record and gates for Stage 0b.
#   GATE 1  installed forestsearch is 0.4.0 (exit 2)
#   GATE 2  speff2trial (ACTG175) is installed (exit 3)
#   RECORD  repo files under study vs CRAN 0.4.0 MD5; fixture sources

source("stage0b_common.R")
unlink(out_path("keyvalues.txt"))

pd <- utils::packageDescription("forestsearch")
kv(installed_version    = as.character(pd$Version),
   installed_packaged   = if (is.null(pd$Packaged)) "<absent>" else pd$Packaged,
   installed_repository = if (is.null(pd$Repository)) "<absent>" else pd$Repository,
   installed_libpath    = dirname(dirname(dirname(attr(pd, "file")))),
   survival_version     = as.character(utils::packageVersion("survival")),
   future_version       = as.character(utils::packageVersion("future")),
   data_table_version   = as.character(utils::packageVersion("data.table")),
   r_version            = R.version.string,
   sysname              = unname(Sys.info()[["sysname"]]),
   available_cores      = parallelly::availableCores(),
   speff2trial_installed = requireNamespace("speff2trial", quietly = TRUE))
if (as.character(pd$Version) != "0.4.0") {
  cat("GATE 1 FAILED: installed forestsearch is", as.character(pd$Version), "- this audit is based on CRAN 0.4.0.\n")
  quit(status = 2)
}
if (!requireNamespace("speff2trial", quietly = TRUE)) {
  cat("GATE 2 FAILED: speff2trial (ACTG175) is not installed.\n")
  quit(status = 3)
}

cran_md5 <- c(
  "DESCRIPTION"                          = "19b9459a03ed5c4a0c12943eb57d82b6",
  "NAMESPACE"                            = "00257d44d6f21d10f9af68b9f0dc227c",
  "R/Cox_estimation_helpers.R"           = "ffa9744bf4b3264041aa4c787cd5958f",
  "R/consistency_resample.R"             = "861c6b4e7d00c8176bbab6b7f863f724",
  "R/forestsearch_main.R"                = "11d80e74c4fcdba43fe8bbf9a1287d3f",
  "R/glm_effect_estimators.R"            = "8e8239f94e5f832af33e4d9cc0ff64d3",
  "R/subgroup_consistency_helpers.R"     = "f32b55b62d683337f6e1b36b3aa3bc17",
  "R/subgroup_consistency_main.R"        = "1285150ed4891883468d01d5ee981b78",
  "R/subgroup_search.R"                  = "f06d6bba8135f41a57962a659ebbe601",
  "R/bootstrap_analysis_dofuture.R"      = "a891d66e7337b9a1ebdf524463102978",
  "R/forestsearch_cross_validation.R"    = "359955f9ba9bd850c978230aee3b4429"
)
root <- tryCatch(system2("git", c("rev-parse", "--show-toplevel"), stdout = TRUE, stderr = FALSE),
                 error = function(e) character(0), warning = function(w) character(0))
if (length(root) != 1L || !nzchar(root)) {
  kv(repo_root = "<not a git repo>")
} else {
  kv(repo_root = root,
     repo_head = system2("git", c("-C", shQuote(root), "rev-parse", "HEAD"), stdout = TRUE),
     repo_branch = system2("git", c("-C", shQuote(root), "rev-parse", "--abbrev-ref", "HEAD"), stdout = TRUE))
  files <- file.path(root, names(cran_md5)); have <- file.exists(files)
  got <- ifelse(have, unname(tools::md5sum(files)), NA_character_)
  tab <- data.frame(file = names(cran_md5), cran_0_4_0 = unname(cran_md5), repo_head = got,
                    match = ifelse(is.na(got), "MISSING", ifelse(got == unname(cran_md5), "match", "DIFFERS")),
                    stringsAsFactors = FALSE)
  utils::write.csv(tab, out_path("provenance_md5.csv"), row.names = FALSE)
  print(tab[, c("file", "match")], row.names = FALSE)
  kv(repo_files_matching_cran = sum(tab$match == "match"), repo_files_differing = sum(tab$match != "match"))
  for (f in tab$file[tab$match != "match"]) {
    cat("--", f, "\n")
    lg <- system2("git", c("-C", shQuote(root), "log", "--oneline", "--since=2026-10-09", "--", shQuote(f)), stdout = TRUE)
    cat(if (length(lg)) paste0("   ", lg) else "   (no commits since 2026-10-09)", sep = "\n")
  }
}

# fixture sources (recorded; the task document says what they must be)
for (nm in names(FIXTURES)) {
  fx <- FIXTURES[[nm]]()
  do.call(kv, setNames(list(fx$label, fx$source, nrow(fx$df), paste(names(fx$args), collapse = ",")),
                       paste0("fixture_", nm, c("_label", "_source", "_n", "_argnames"))))
}
cat("\nb0 complete\n")
```

### `b1_profiles.R`

```r
# b1_profiles.R -- per fixture (survival V, continuous, binary): one warm-up
# call, candidate-flow counts, 3 timed sequential calls, 3 Rprof profiles
# bucketed by stage, and a sequential nb = 5 bootstrap for the per-replicate
# denominator.  Nothing is modified; the installed package is called as a
# user would.

source("stage0b_common.R")

watch <- c(
  # search stage
  "search_combinations_parallel", "evaluate_combination_with_status",
  "fit_cox_for_subgroup", "fit_glm_for_subgroup", "fit_glm_for_subgroup_fast",
  ".fs_lm_two_group", ".estimate_or", ".estimate_rd", ".estimate_rr",
  "coxph", "coxph.fit", "glm", "glm.fit", "lm", "lm.fit",
  "model.frame.default", "model.matrix.default", "summary.glm", "vcov", "concordancefit",
  "survfit", "survfit.formula",
  # consistency stage
  "evaluate_consistency_twostage", "evaluate_subgroup_consistency",
  "consistency_resample", ".consistency_glm_pieces", ".consistency_cox_pieces",
  ".dfbeta_glm", "lm.influence", "residuals.coxph",
  "as.data.table", "[.data.table", "FS_labels",
  # selection / dedup
  "remove_near_duplicate_subgroups", ".maxeff_membership_dedup",
  "sort_subgroups_preview", "sort_subgroups",
  # screens / cuts
  "get_FSdata", "lasso_selection", "glmnet", "cv.glmnet",
  "causal_forest", "causal_survival_forest", "fit_causal_forest_glm",
  "variable_importance", "policy_tree", "dummy", "get_dfpred", "merge")

prof_one <- function(fx, tag, i) {
  f <- out_path(sprintf("Rprof_%s_%d.out", tag, i))
  suppressWarnings(Rprof(f, interval = 0.002))
  run_fixture_call(fx)
  Rprof(NULL)
  sp <- summaryRprof(f)
  bt <- sp$by.total
  bt$fun <- sub("^[A-Za-z0-9._]+::", "", gsub('"', "", rownames(bt)))
  bt <- aggregate(cbind(total.time, total.pct, self.time, self.pct) ~ fun, data = bt, FUN = sum)
  bt$rep <- i; bt$fixture <- tag; bt$sampling_time_s <- sp$sampling.time
  bt
}

all_med <- list(); all_prof <- list()
for (tag in names(FIXTURES)) {
  cat("\n==========", tag, "==========\n")
  fx <- FIXTURES[[tag]]()
  t_warm <- system.time(fs <- run_fixture_call(fx))[["elapsed"]]
  cc <- candidate_counts(fs)
  do.call(kv, setNames(c(list(t_warm), cc), paste0(tag, "_", c("call_warmup_s", names(cc)))))
  if (tag == "survival") {
    # GATE (as in Stage 0): the survival candidate flow is fixed by the data and
    # the 0.4.0 code, not by the machine.
    expected <- list(n_reached_fit = 1410L, n_passed_effect = 121L,
                     n_consistency_total = 120L, sg_harm = "{er <= 0} & {pgr <= 26}")
    bad <- names(expected)[!mapply(function(k, v) identical(as.character(cc[[k]]), as.character(v)),
                                   names(expected), expected)]
    if (length(bad)) {
      cat("GATE FAILED: survival candidate flow differs from CRAN 0.4.0 on fixture V for:",
          paste(bad, collapse = ", "), "\n")
      quit(status = 5)
    }
  }

  t_call <- vapply(1:3, function(i) system.time(run_fixture_call(fx))[["elapsed"]], numeric(1))
  do.call(kv, setNames(list(t_call, median(t_call)), paste0(tag, c("_call_s_rep", "_call_s_median"))))

  prof <- do.call(rbind, lapply(1:3, function(i) prof_one(fx, tag, i)))
  all_prof[[tag]] <- prof
  w <- prof[prof$fun %in% watch, c("fun", "total.time", "total.pct", "rep")]
  med <- aggregate(cbind(total.time, total.pct) ~ fun, data = w, FUN = median)
  med <- med[order(-med$total.time), ]
  all_med[[tag]] <- cbind(fixture = tag, med)
  cat("\nRprof by.total, median of 3 (s):\n"); print(med, row.names = FALSE)
  get_t <- function(f) { v <- med$total.time[med$fun == f]; if (length(v)) v else 0 }
  do.call(kv, setNames(list(
    median(prof$sampling_time_s[!duplicated(prof$rep)]),
    get_t("search_combinations_parallel"),
    get_t("fit_cox_for_subgroup") + get_t("fit_glm_for_subgroup") + get_t("fit_glm_for_subgroup_fast"),
    get_t("evaluate_consistency_twostage") + get_t("evaluate_subgroup_consistency"),
    get_t("consistency_resample"),
    get_t(".consistency_glm_pieces") + get_t(".consistency_cox_pieces"),
    get_t("remove_near_duplicate_subgroups") + get_t(".maxeff_membership_dedup") + get_t("sort_subgroups_preview"),
    get_t("causal_forest") + get_t("causal_survival_forest") + get_t("fit_causal_forest_glm") + get_t("policy_tree"),
    get_t("lasso_selection"), get_t("get_FSdata"), get_t("survfit"),
    get_t("coxph"), get_t("coxph.fit"), get_t("glm"), get_t("glm.fit"), get_t("lm"), get_t("lm.fit"),
    get_t("lm.influence"), get_t(".dfbeta_glm"), get_t("model.frame.default"), get_t("model.matrix.default"),
    get_t("as.data.table")),
    paste0(tag, "_rprof_", c("sampling_s", "search_stage_s", "search_fit_s", "consistency_stage_s",
                             "consistency_resample_s", "consistency_pieces_s", "dedup_sort_s", "grf_s",
                             "lasso_s", "get_FSdata_s", "survfit_s", "coxph_s", "coxph_fit_s", "glm_s", "glm_fit_s",
                             "lm_s", "lm_fit_s", "lm_influence_s", "dfbeta_glm_s", "model_frame_s", "model_matrix_s",
                             "as_data_table_s"))))

  # sequential replicate denominator: nb = 5 bootstrap, plan sequential
  if (!is.null(fs$sg.harm)) {
    t0 <- proc.time()[["elapsed"]]
    b <- suppressWarnings(suppressMessages(forestsearch::forestsearch_bootstrap_dofuture(
      fs.est = fs, nb_boots = 5L, seed = 8316951L,
      parallel_args = list(plan = "sequential", workers = 1L, show_message = FALSE))))
    t_outer <- proc.time()[["elapsed"]] - t0
    future::plan("sequential")
    rep_s <- b$results$tmins_iteration * 60
    do.call(kv, setNames(list(t_outer, mean(rep_s), sum(b$results$any_found == 1L)),
                         paste0(tag, c("_boot5_outer_s", "_rep_s_mean", "_boot5_n_found"))))
  } else {
    do.call(kv, setNames(list(NA, NA, NA), paste0(tag, c("_boot5_outer_s", "_rep_s_mean", "_boot5_n_found"))))
    cat("NOTE: no subgroup identified on this fixture; bootstrap denominator unavailable\n")
  }
}
utils::write.csv(do.call(rbind, all_prof), out_path("rprof_by_total_all.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, all_med), out_path("rprof_watched_median.csv"), row.names = FALSE)
cat("\nb1 complete\n")
```

### `b2_lean.R`

```r
# b2_lean.R -- the candidates sized and their bit-identity tested, per fixture,
# with NO change to the package.  One traced forestsearch() call per fixture
# records (in this session only, untraced before any timing):
#   * the consistency-stage evaluator's arguments
#     (evaluate_consistency_twostage / evaluate_subgroup_consistency): the
#     shared objects once, the candidate index m per call;
#   * the consistency "pieces" inputs per evaluated candidate
#     (.consistency_glm_pieces / .consistency_cox_pieces): the needed columns
#     of the subgroup frame, and cox_init for survival;
#   * the search-stage fit inputs per candidate (binary:
#     fit_glm_for_subgroup; continuous: fit_glm_for_subgroup_fast).
# Then, per candidate, the package's own function (unmodified) is run on the
# recorded arguments and timed against a lean equivalent, and the outputs are
# compared with identical().
#
# Candidates measured here (letters as in the survey / Stage 0):
#   B  survival consistency pieces: coxph.fit + residuals.coxph on a constructed
#      object instead of two formula coxph() fits + a model-frame re-evaluation
#   E  continuous consistency pieces: lm.fit + QR hat values instead of
#      lm() + lm.influence()
#   F  binary search-stage fit: glm.fit on vectors instead of glm() on a
#      data-frame slice, SE from the fit's QR instead of vcov()/summary.glm()
#   G  binary consistency pieces: glm.fit + QR hat values instead of
#      glm() + lm.influence()
#   H  consistency-stage member extraction: a logical index over the cut
#      columns instead of as.data.table(df) + sequential row filters, per
#      candidate; measured end to end as the whole evaluator (package) vs a
#      lean evaluator (H + the lean pieces of B/E/G).
# A (survival search fit) is Stage 0's and is quoted, not re-measured.

source("stage0b_common.R")

# ------------------------------------------------------------------------------
# lean pieces (vectors in; beta_hat, sigma_D out; NULL where the package's
# function returns NULL).  The dfbeta arithmetic is .dfbeta_glm()'s, line for
# line, with lm.influence()'s hat values recomputed from the fit's own QR.
# ------------------------------------------------------------------------------
hat_from_qr <- function(qr, n, k) {
  h <- rowSums(base::qr.qy(qr, diag(1, n, k))^2)
  h[h > 1 - 10 * .Machine$double.eps] <- 1
  h
}
lean_glm_pieces <- function(y, z, effect_measure, adverse_outcome, outcome_type, treat.name) {
  if (!adverse_outcome) { if (outcome_type == "continuous") y <- -y else if (outcome_type == "binary") y <- 1L - y }
  n0 <- sum(z == 0L, na.rm = TRUE); n1 <- sum(z == 1L, na.rm = TRUE)
  if (n0 < 3L || n1 < 3L) return(NULL)
  X <- cbind("(Intercept)" = 1, as.double(z)); colnames(X)[2L] <- treat.name
  if (effect_measure == "MD") {
    fit <- stats::lm.fit(X, as.double(y)); w <- rep(1, nrow(X))
  } else if (effect_measure == "OR") {
    fit <- tryCatch(suppressWarnings(stats::glm.fit(X, as.double(y), family = stats::binomial("logit"))),
                    error = function(e) NULL)
    if (is.null(fit) || isFALSE(fit$converged)) return(NULL)
    w <- fit$weights
  } else return(NULL)
  cf <- fit$coefficients
  if (any(is.na(cf))) return(NULL)
  beta_hat <- unname(cf[[treat.name]])
  h <- hat_from_qr(fit$qr, nrow(X), fit$rank)
  A <- solve(crossprod(X, w * X))
  r <- sqrt(w) * fit$residuals
  dfb <- t(A %*% t(X * (sqrt(w) * r / (1 - h))))
  dfb <- as.numeric(dfb[, 2L])
  sigma_D <- sqrt(sum(dfb^2))
  if (!is.finite(beta_hat) || !is.finite(sigma_D) || sigma_D <= 0) return(NULL)
  list(beta_hat = beta_hat, sigma_D = sigma_D)
}

# survival: the conventions coxph() passes to coxph.fit() are read at run time
# (as in Stage 0); `init` is taken from the captured second fit (the warm start).
.coxfit_conv <- NULL
.cox_terms <- stats::terms(survival::Surv(Y, E) ~ Treat, specials = c("strata", "cluster", "tt"))
lean_cox_fit <- function(Y, E, Treat, init) {
  a <- .coxfit_conv
  n <- length(Y)
  ys <- survival::Surv(Y, E)
  if (isTRUE(a$control$timefix)) ys <- survival::aeqSurv(ys)
  args <- c(list(x = matrix(as.double(Treat), ncol = 1L, dimnames = list(NULL, "Treat")), y = ys,
                 offset = if (is.null(a$offset)) NULL else rep(0, n),
                 rownames = if (is.null(a$rownames)) NULL else as.character(seq_len(n)),
                 init = init),
            a[setdiff(names(a), c("x", "y", "offset", "rownames", "init"))])
  do.call(survival::coxph.fit, args)
}
lean_cox_pieces <- function(Y, E, Treat) {
  n <- length(Y)
  if (n < 2L || sum(E) < 2L) return(NULL)
  fit0 <- tryCatch(suppressWarnings(lean_cox_fit(Y, E, Treat, NULL)), error = function(e) NULL)
  cox_init <- if (is.null(fit0)) 0 else fit0$coefficients[1]
  fit <- tryCatch(suppressWarnings(lean_cox_fit(Y, E, Treat, cox_init)), error = function(e) NULL)
  if (is.null(fit)) return(NULL)
  cf <- fit$coefficients
  if (length(cf) == 0L || any(is.na(cf))) return(NULL)
  beta_hat <- unname(cf[[1L]])
  obj <- fit
  obj$x <- matrix(as.double(Treat), ncol = 1L, dimnames = list(NULL, "Treat"))
  ys <- survival::Surv(Y, E); if (isTRUE(.coxfit_conv$control$timefix)) ys <- survival::aeqSurv(ys)
  obj$y <- ys; obj$strata <- NULL; obj$method <- .coxfit_conv$method; obj$weights <- NULL
  obj$terms <- .cox_terms      # residuals.coxph() requires a terms component
  class(obj) <- "coxph"
  dfbeta <- tryCatch(stats::residuals(obj, type = "dfbeta"), error = function(e) NULL)
  if (is.null(dfbeta)) return(NULL)
  dfbeta <- as.numeric(dfbeta)
  sigma_D <- sqrt(sum(dfbeta^2))
  if (!is.finite(beta_hat) || !is.finite(sigma_D) || sigma_D <= 0) return(NULL)
  list(beta_hat = beta_hat, sigma_D = sigma_D, cox_init = cox_init)
}

# ------------------------------------------------------------------------------
# per fixture
# ------------------------------------------------------------------------------
summ <- list()
fixture_sel <- strsplit(Sys.getenv("STAGE0B_FIXTURES", "survival,continuous,binary"), ",")[[1L]]
for (tag in intersect(names(FIXTURES), fixture_sel)) {
  cat("\n==========", tag, "==========\n")
  fx <- FIXTURES[[tag]]()
  is_surv <- tag == "survival"
  out_nm <- fx$args$outcome.name; trt_nm <- fx$args$treat.name

  # ---- 1. traced call: capture -------------------------------------------------
  .stage <- new.env(); .stage$args <- NULL; .stage$m <- integer(0); .stage$fun <- NA_character_
  assign(".stage", .stage, envir = globalenv())
  stage_tracer <- function(fname) bquote({
    .fr <- environment()
    if (is.null(.GlobalEnv$.stage$args)) {
      .nm <- setdiff(names(formals(get(.(fname), envir = asNamespace("forestsearch")))), c("m", "..."))
      .GlobalEnv$.stage$args <- lapply(setNames(.nm, .nm), function(.a) get(.a, envir = .fr, inherits = FALSE))
      .GlobalEnv$.stage$fun <- .(fname)
    }
    .GlobalEnv$.stage$m <- c(.GlobalEnv$.stage$m, m)
  })
  invisible(capture.output({
    trace("evaluate_consistency_twostage", where = asNamespace("forestsearch"), print = FALSE,
          tracer = stage_tracer("evaluate_consistency_twostage"))
    trace("evaluate_subgroup_consistency", where = asNamespace("forestsearch"), print = FALSE,
          tracer = stage_tracer("evaluate_subgroup_consistency"))
  }))
  .pieces <- new_cap(); assign(".pieces", .pieces, envir = globalenv())
  if (is_surv) {
    invisible(capture.output(trace(".consistency_cox_pieces", where = asNamespace("forestsearch"), print = FALSE,
      tracer = quote({ .d <- as.data.frame(df)
        .GlobalEnv$.pieces$push(list(Y = .d[[tte.name]], E = .d[[event.name]], Treat = .d[[treat.name]],
                                     cox_init = cox_init, adj = adjust_covariates, ids = .d[["id"]])) }))))
  } else {
    invisible(capture.output(trace(".consistency_glm_pieces", where = asNamespace("forestsearch"), print = FALSE,
      tracer = quote({ .d <- as.data.frame(df)
        .GlobalEnv$.pieces$push(list(y = .d[[outcome.name]], z = .d[[treat.name]], ids = .d[["id"]],
                                     offset = if (!is.null(offset.name)) .d[[offset.name]] else NULL))
        if (is.null(.GlobalEnv$.pieces$conv))
          .GlobalEnv$.pieces$conv <- list(outcome_type = outcome_type, effect_measure = effect_measure,
            treat.name = treat.name, outcome.name = outcome.name, offset.name = offset.name,
            adjust_covariates = adjust_covariates, adverse_outcome = adverse_outcome) }))))
  }
  .sfit <- new_cap(); assign(".sfit", .sfit, envir = globalenv())
  if (tag == "binary") {
    invisible(capture.output(trace("fit_glm_for_subgroup", where = asNamespace("forestsearch"), print = FALSE,
      tracer = quote({ .GlobalEnv$.sfit$push(list(id.x = id.x))
        if (is.null(.GlobalEnv$.sfit$df_clean)) { .GlobalEnv$.sfit$df_clean <- df_clean; .GlobalEnv$.sfit$estimator_fn <- estimator_fn } }))))
  } else if (tag == "continuous") {
    invisible(capture.output(trace("fit_glm_for_subgroup_fast", where = asNamespace("forestsearch"), print = FALSE,
      tracer = quote({ .GlobalEnv$.sfit$push(list(id.x = id.x))
        if (is.null(.GlobalEnv$.sfit$yy)) { .GlobalEnv$.sfit$yy <- yy; .GlobalEnv$.sfit$tt <- tt; .GlobalEnv$.sfit$estimator_fn <- estimator_fn } }))))
  }
  t_traced <- system.time(fs <- run_fixture_call(fx))[["elapsed"]]
  invisible(capture.output({
    untrace("evaluate_consistency_twostage", where = asNamespace("forestsearch"))
    untrace("evaluate_subgroup_consistency", where = asNamespace("forestsearch"))
    if (is_surv) untrace(".consistency_cox_pieces", where = asNamespace("forestsearch"))
    else untrace(".consistency_glm_pieces", where = asNamespace("forestsearch"))
    if (tag == "binary") untrace("fit_glm_for_subgroup", where = asNamespace("forestsearch"))
    if (tag == "continuous") untrace("fit_glm_for_subgroup_fast", where = asNamespace("forestsearch"))
  }))
  cc <- candidate_counts(fs)
  P <- cap_items(.pieces); S <- cap_items(.sfit); M <- .stage$m; SA <- .stage$args
  kv_tag <- function(...) { v <- list(...); do.call(kv, setNames(v, paste0(tag, "_", names(v)))) }
  kv_tag(capture_call_s = t_traced, capture_stage_fun = .stage$fun,
         capture_n_stage_calls = length(M), capture_n_pieces_calls = length(P),
         capture_n_searchfit_calls = length(S),
         capture_stage_calls_match_counts = identical(length(M), as.integer(cc$n_consistency_evaluated)))
  if (length(M) == 0L || is.null(SA)) {
    cat("GATE FAILED: no consistency-stage call captured for", tag, "\n"); quit(status = 4)
  }
  if (is_surv) {
    # the conventions coxph() passes to coxph.fit(): captured during ONE package pieces call
    .capfit <- new.env(); .capfit$calls <- list(); assign(".capfit", .capfit, envir = globalenv())
    invisible(capture.output(trace("coxph.fit", where = asNamespace("survival"), print = FALSE, tracer = quote({
      .fr <- environment(); .nm <- setdiff(names(formals(survival::coxph.fit)), "...")
      .GlobalEnv$.capfit$calls[[length(.GlobalEnv$.capfit$calls) + 1L]] <-
        lapply(setNames(.nm, .nm), function(.a) tryCatch(get(.a, envir = .fr, inherits = FALSE), error = function(e) NULL)) }))))
    p1 <- P[[1L]]
    d1 <- data.frame(Y = p1$Y, Event = p1$E, Treat = p1$Treat)
    invisible(forestsearch:::.consistency_cox_pieces(d1, "Y", "Event", "Treat", cox_init = p1$cox_init, adjust_covariates = NULL))
    invisible(capture.output(untrace("coxph.fit", where = asNamespace("survival"))))
    cf <- .capfit$calls
    kv_tag(coxfit_calls_per_pieces = length(cf),
           coxfit_init_as_passed = if (length(cf)) paste(deparse(cf[[length(cf)]]$init), collapse = "") else "<none>",
           coxfit_init_equals_cox_init = if (length(cf)) isTRUE(all.equal(unname(cf[[length(cf)]]$init), unname(p1$cox_init))) else NA,
           coxfit_nocenter = if (length(cf)) paste(deparse(cf[[1]]$nocenter), collapse = "") else "<none>",
           coxfit_timefix = isTRUE(cf[[1]]$control$timefix))
    .coxfit_conv <<- cf[[length(cf)]]
  }

  # ---- 2. search-stage fit (binary: F; continuous: bucket only) --------------
  if (tag == "binary") {
    dfc <- .sfit$df_clean; efn <- .sfit$estimator_fn
    yy <- dfc[[out_nm]]; tt <- dfc[[trt_nm]]
    fn_env <- environment(efn)
    adv <- isTRUE(fn_env$adverse_outcome)
    pkg_fit <- forestsearch:::fit_glm_for_subgroup
    lean_fit <- function(id.x) {
      idx <- id.x == 1
      if (sum(idx) < 6L) return(NULL)
      y <- yy[idx]; z <- tt[idx]
      if (!adv) y <- 1L - y
      n0 <- sum(z == 0L, na.rm = TRUE); n1 <- sum(z == 1L, na.rm = TRUE)
      why <- forestsearch:::.fs_existence_reason("OR", forestsearch:::.fs_binary_cells(y, z))
      if (!is.null(why)) return(list(nonestimable = TRUE, reason = why))
      X <- cbind("(Intercept)" = 1, as.double(z)); colnames(X)[2L] <- trt_nm
      fit <- tryCatch(suppressWarnings(stats::glm.fit(X, as.double(y), family = stats::binomial("logit"))),
                      error = function(e) NULL)
      if (is.null(fit)) return(NULL)
      est <- fit$coefficients[[trt_nm]]
      if (isFALSE(fit$converged) && !is.na(est)) return(list(nonestimable = TRUE, reason = "non-convergent fit"))
      if (is.na(est)) return(NULL)
      p1 <- seq_len(fit$rank)
      se <- unname(sqrt(diag(1 * chol2inv(fit$qr$qr[p1, p1, drop = FALSE])))[2L])
      z_crit <- 1.96
      list(hr = est, lower = est - z_crit * se, upper = est + z_crit * se, med0 = NA_real_, med1 = NA_real_)
    }
    res_f <- lapply(S, function(s) pkg_fit(dfc, s$id.x, efn))
    res_l <- lapply(S, function(s) lean_fit(s$id.x))
    same <- mapply(function(u, v) identical(u, v), res_f, res_l)
    kind <- function(u) if (is.null(u)) "NULL" else if (isTRUE(u$nonestimable)) "nonestimable" else "fit"
    kf <- vapply(res_f, kind, ""); kl <- vapply(res_l, kind, "")
    tp <- time_pair(function() for (s in S) pkg_fit(dfc, s$id.x, efn), function() for (s in S) lean_fit(s$id.x))
    kv_tag(F_n = length(S), F_identical = sum(same), F_kinds_agree = identical(kf, kl),
           F_n_fit = sum(kf == "fit"), F_n_nonestimable = sum(kf == "nonestimable"), F_n_null = sum(kf == "NULL"),
           F_formula_s_reps = tp$a, F_lean_s_reps = tp$b, F_formula_s = tp$a_med, F_lean_s = tp$b_med,
           F_machinery_s = tp$a_med - tp$b_med, F_machinery_share = (tp$a_med - tp$b_med) / tp$a_med)
    if (any(!same)) { k <- which(!same)[1]; cat("first F mismatch at", k, "\n"); str(res_f[[k]]); str(res_l[[k]]) }
  } else if (tag == "continuous") {
    yy <- .sfit$yy; tt <- .sfit$tt; efn <- .sfit$estimator_fn
    t_fast <- vapply(1:3, function(i) system.time(for (s in S) forestsearch:::fit_glm_for_subgroup_fast(yy, tt, s$id.x, efn))[["elapsed"]], numeric(1))
    kv_tag(searchfit_fast_n = length(S), searchfit_fast_s = median(t_fast))
  }

  # ---- 3. pieces alone (B / E / G) -------------------------------------------
  if (is_surv) {
    pkg_pieces <- function(p) forestsearch:::.consistency_cox_pieces(
      data.frame(Y = p$Y, Event = p$E, Treat = p$Treat), "Y", "Event", "Treat",
      cox_init = p$cox_init, adjust_covariates = p$adj)
    lean_pieces <- function(p) lean_cox_pieces(p$Y, p$E, p$Treat)
    c_cmp <- log(SA$hr.consistency)
  } else {
    cv <- .pieces$conv
    pkg_pieces <- function(p) {
      d <- data.frame(y = p$y, z = p$z); names(d) <- c(cv$outcome.name, cv$treat.name)
      if (!is.null(p$offset)) d[[cv$offset.name]] <- p$offset
      forestsearch:::.consistency_glm_pieces(d, cv$outcome_type, cv$effect_measure, cv$treat.name,
                                             cv$outcome.name, cv$offset.name, cv$adjust_covariates, cv$adverse_outcome)
    }
    lean_pieces <- function(p) lean_glm_pieces(p$y, p$z, cv$effect_measure, cv$adverse_outcome, cv$outcome_type, cv$treat.name)
    c_cmp <- SA$glm_resample_spec$comparison_threshold
  }
  res_f <- lapply(P, pkg_pieces); res_l <- lapply(P, lean_pieces)
  nf <- vapply(res_f, is.null, logical(1)); nl <- vapply(res_l, is.null, logical(1))
  id_beta  <- mapply(function(u, v) !is.null(u) && !is.null(v) && identical(u$beta_hat, v$beta_hat), res_f, res_l)
  id_sigma <- mapply(function(u, v) !is.null(u) && !is.null(v) && identical(u$sigma_D, v$sigma_D), res_f, res_l)
  rate_of <- function(r) if (is.null(r)) NA_real_ else max(0, 2 * stats::pnorm((r$beta_hat - c_cmp) / r$sigma_D) - 1)
  id_rate <- mapply(function(u, v) identical(rate_of(u), rate_of(v)), res_f, res_l)
  tp <- time_pair(function() for (p in P) pkg_pieces(p), function() for (p in P) lean_pieces(p))
  L <- if (is_surv) "B" else if (tag == "continuous") "E" else "G"
  kv_tag(pieces_letter = L, pieces_n = length(P), pieces_null_formula = sum(nf), pieces_null_lean = sum(nl),
         pieces_null_sets_agree = identical(nf, nl), pieces_beta_identical = sum(id_beta),
         pieces_sigma_identical = sum(id_sigma), pieces_rate_identical = sum(id_rate),
         pieces_formula_s_reps = tp$a, pieces_lean_s_reps = tp$b,
         pieces_formula_s = tp$a_med, pieces_lean_s = tp$b_med,
         pieces_machinery_s = tp$a_med - tp$b_med, pieces_machinery_share = (tp$a_med - tp$b_med) / tp$a_med,
         pieces_lean_identical_all = all(id_rate))
  if (any(!id_rate)) { k <- which(!id_rate)[1]; cat("first pieces mismatch at", k, ":\n"); str(res_f[[k]]); str(res_l[[k]]) }
  pieces_identical_all <- all(id_rate)

  # ---- 4. the whole consistency evaluator: package vs lean (H + pieces) -------
  stage_fun <- get(.stage$fun, envir = asNamespace("forestsearch"))
  pkg_stage <- function(m) do.call(stage_fun, c(list(m = m), SA))
  df_all <- as.data.frame(SA$df)
  index_mat <- as.matrix(SA$index.Z); storage.mode(index_mat) <- "double"
  names.Z <- SA$names.Z
  col_cache <- lapply(setNames(names.Z, names.Z), function(v) df_all[[v]] == 1)
  pdig <- SA$pconsistency.digits; pthr <- SA$pconsistency.threshold
  if (is_surv) { Y_all <- df_all$Y; E_all <- df_all$Event; T_all <- df_all$Treat } else {
    y_all <- df_all[[cv$outcome.name]]; z_all <- df_all[[cv$treat.name]]
    off_all <- if (!is.null(cv$offset.name)) df_all[[cv$offset.name]] else NULL }
  lean_stage <- function(m, use_lean_pieces = TRUE) {
    this.m <- names.Z[index_mat[m, ] == 1]
    keep <- col_cache[[this.m[1L]]]
    for (v in this.m[-1L]) keep <- keep & col_cache[[v]]
    N.x <- sum(keep)
    if (N.x < 10) return(NULL)
    rr <- if (is_surv) {
      if (use_lean_pieces) lean_cox_pieces(Y_all[keep], E_all[keep], T_all[keep])
      else pkg_pieces(list(Y = Y_all[keep], E = E_all[keep], Treat = T_all[keep],
                           cox_init = lean_cox_fit(Y_all[keep], E_all[keep], T_all[keep], NULL)$coefficients[1], adj = NULL))
    } else {
      if (use_lean_pieces) lean_glm_pieces(y_all[keep], z_all[keep], cv$effect_measure, cv$adverse_outcome, cv$outcome_type, cv$treat.name)
      else pkg_pieces(list(y = y_all[keep], z = z_all[keep], offset = if (is.null(off_all)) NULL else off_all[keep]))
    }
    if (is.null(rr)) return(structure(list(), class = "stage0b_fallback"))  # package would fall through to splitting
    p <- round(rate_of(rr), pdig)
    if (p < pthr) return(NULL)
    p
  }
  res_pkg <- lapply(M, pkg_stage)
  res_lean <- lapply(M, lean_stage)
  pcons_pkg <- vapply(res_pkg, function(r) if (is.null(r)) NA_real_ else as.numeric(r[["Pcons"]]), numeric(1))
  pcons_lean <- vapply(res_lean, function(r) if (is.null(r) || inherits(r, "stage0b_fallback")) NA_real_ else r, numeric(1))
  n_fallback <- sum(vapply(res_lean, inherits, logical(1), "stage0b_fallback"))
  stage_ident <- identical(pcons_pkg, pcons_lean)
  tp <- time_pair(function() for (m in M) pkg_stage(m), function() for (m in M) lean_stage(m))
  # variant: lean extraction with the package's own pieces (what H alone gives)
  tH <- vapply(1:3, function(i) system.time(for (m in M) lean_stage(m, use_lean_pieces = FALSE))[["elapsed"]], numeric(1))
  res_H <- lapply(M, function(m) lean_stage(m, use_lean_pieces = FALSE))
  pcons_H <- vapply(res_H, function(r) if (is.null(r) || inherits(r, "stage0b_fallback")) NA_real_ else r, numeric(1))
  kv_tag(stage_n = length(M), stage_n_passed_pkg = sum(!is.na(pcons_pkg)), stage_n_passed_lean = sum(!is.na(pcons_lean)),
         stage_pcons_identical_all = stage_ident, stage_n_lean_fallback = n_fallback,
         stage_pcons_identical_H_only = identical(pcons_pkg, pcons_H),
         stage_formula_s_reps = tp$a, stage_lean_s_reps = tp$b, stage_formula_s = tp$a_med, stage_lean_s = tp$b_med,
         stage_lean_Honly_s = median(tH),
         stage_machinery_s = tp$a_med - tp$b_med, stage_machinery_share = (tp$a_med - tp$b_med) / tp$a_med,
         stage_H_only_machinery_s = tp$a_med - median(tH))
  # extraction overhead = whole package evaluator minus the package pieces alone (both on the same candidates)
  kv_tag(stage_extraction_overhead_s = tp$a_med - median(vapply(1:3, function(i) system.time(for (p in P) pkg_pieces(p))[["elapsed"]], numeric(1))))
  if (!stage_ident) { k <- which(!(is.na(pcons_pkg) == is.na(pcons_lean)) | (pcons_pkg != pcons_lean))[1]; cat("first stage mismatch at m =", M[k], ": pkg", pcons_pkg[k], "lean", pcons_lean[k], "\n") }

  summ[[tag]] <- list(cc = cc, P_n = length(P), S_n = length(S), M_n = length(M),
                      pieces_identical_all = pieces_identical_all, stage_identical = stage_ident)
  save_result(list(cc = cc, M = M, pcons_pkg = pcons_pkg, pcons_lean = pcons_lean,
                   stage_args_names = names(SA), coxfit_conv = if (is_surv) .coxfit_conv else NULL), paste0("b2_", tag))
}
save_result(summ, "b2_summary")
cat("\nb2 complete\n")
```

### `b3_summary.R`

```r
# b3_summary.R -- assembles results/stage0b_summary.md from results/keyvalues.txt:
# per fixture the call and replicate times, the Rprof buckets, each candidate's
# measured machinery seconds and share of a replicate against the 30% bar, the
# identity counts, and the modification map (where package functions would
# change) with the measured status of each entry.  Pure arithmetic on
# recorded numbers.

source("stage0b_common.R")
lines <- readLines(out_path("keyvalues.txt"))
kvs <- strsplit(lines, " = ", fixed = TRUE)
vals <- setNames(lapply(kvs, function(z) paste(z[-1L], collapse = " = ")), vapply(kvs, `[[`, "", 1L))
vals <- vals[!duplicated(names(vals), fromLast = TRUE)]
num <- function(k) { v <- vals[[k]]; if (is.null(v)) NA_real_ else suppressWarnings(as.numeric(strsplit(v, ", ")[[1L]][1L])) }
chr <- function(k) { v <- vals[[k]]; if (is.null(v)) "<missing>" else v }
pct <- function(x) if (is.na(x)) "n/a" else sprintf("%.1f%%", 100 * x)
s   <- function(x, d = 2) if (is.na(x)) "n/a" else formatC(x, format = "f", digits = d)
BAR <- 0.30

# Stage-0 (committed) survival A figures, quoted for the survival row
A_mach_s <- 1.509; A_rep_s <- 3.2112

md <- c("# Stage 0b -- GLM paths and the consistency stage (CRAN forestsearch 0.4.0)", "",
  sprintf("Installed: forestsearch %s (%s; repository %s) · survival %s · data.table %s · %s · cores %s",
          chr("installed_version"), chr("installed_packaged"), chr("installed_repository"), chr("survival_version"),
          chr("data_table_version"), chr("r_version"), chr("available_cores")),
  sprintf("Repo: %s @ %s (%s); files matching CRAN 0.4.0: %s, differing: %s",
          chr("repo_root"), chr("repo_head"), chr("repo_branch"), chr("repo_files_matching_cran"), chr("repo_files_differing")),
  "")

# ---- per-fixture profile table ----------------------------------------------------
md <- c(md, "## 1. The single sequential call and the replicate, per fixture", "",
  "| fixture | source | n | call s (median of 3) | sequential replicate s (nb = 5) | search: evaluated / reached fit / passed effect | consistency: total / evaluated / passed | early stop | identified subgroup |",
  "|---|---|---|---|---|---|---|---|---|")
for (t in names(FIXTURES)) md <- c(md, sprintf("| %s | %s | %s | %s | %s | %s / %s / %s | %s / %s / %s | %s | %s |",
  t, chr(paste0("fixture_", t, "_source")), chr(paste0("fixture_", t, "_n")), s(num(paste0(t, "_call_s_median")), 3),
  s(num(paste0(t, "_rep_s_mean")), 3), chr(paste0(t, "_n_evaluated_search")), chr(paste0(t, "_n_reached_fit")),
  chr(paste0(t, "_n_passed_effect")), chr(paste0(t, "_n_consistency_total")), chr(paste0(t, "_n_consistency_evaluated")),
  chr(paste0(t, "_n_consistency_passed")), chr(paste0(t, "_early_stop_triggered")), chr(paste0(t, "_sg_harm"))))

md <- c(md, "", "Rprof buckets, median of 3 profiles, seconds of total time (a function reached through do.call() is not attributable by name; its children are):", "",
  "| fixture | call s | search stage | of which per-candidate fit | consistency stage | of which pieces (fit + dfbeta) | consistency stage minus pieces (extraction, bookkeeping) | dedup / sort | GRF | LASSO | survfit medians |",
  "|---|---|---|---|---|---|---|---|---|---|---|")
for (t in names(FIXTURES)) md <- c(md, sprintf("| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |",
  t, s(num(paste0(t, "_rprof_sampling_s"))), s(num(paste0(t, "_rprof_search_stage_s"))), s(num(paste0(t, "_rprof_search_fit_s"))),
  s(num(paste0(t, "_rprof_consistency_stage_s"))), s(num(paste0(t, "_rprof_consistency_pieces_s"))),
  s(num(paste0(t, "_rprof_consistency_stage_s")) - num(paste0(t, "_rprof_consistency_pieces_s"))),
  s(num(paste0(t, "_rprof_dedup_sort_s"))), s(num(paste0(t, "_rprof_grf_s"))), s(num(paste0(t, "_rprof_lasso_s"))), s(num(paste0(t, "_rprof_survfit_s")))))

# ---- candidates --------------------------------------------------------------------
md <- c(md, "", "## 2. Candidates: machinery removed per call, identity, and the 30% bar (per-replicate shape)", "",
  "Machinery = (package function on the recorded inputs) - (lean equivalent on the same inputs), direct timers, median of 3 alternating reps. Share = machinery / measured sequential replicate. Identity = identical() on what the package consumes.", "",
  "| candidate | fixture | what | n calls | package s | lean s | machinery s | share of replicate | identity | clears 30%? |",
  "|---|---|---|---|---|---|---|---|---|---|")
row <- function(letter, t, what, n, a, b, ident) {
  rep_s <- num(paste0(t, "_rep_s_mean")); mach <- a - b; sh <- mach / rep_s
  sprintf("| **%s** | %s | %s | %s | %s | %s | **%s** | **%s** | %s | %s |", letter, t, what, n, s(a, 3), s(b, 3), s(mach, 3), pct(sh), ident,
          if (is.na(sh)) "n/a" else if (sh >= BAR) "YES" else "no")
}
# survival A (Stage 0, quoted)
md <- c(md, sprintf("| **A** | survival | search-stage Cox fit (Stage 0, committed) | 1410 | 1.887 | 0.378 | **1.509** | **%s** | 1410/1410 | YES |", pct(A_mach_s / A_rep_s)))
for (t in names(FIXTURES)) {
  L <- chr(paste0(t, "_pieces_letter"))
  md <- c(md, row(L, t, "consistency pieces (fit + dfbeta), narrow frame",
                  chr(paste0(t, "_pieces_n")), num(paste0(t, "_pieces_formula_s")), num(paste0(t, "_pieces_lean_s")),
                  sprintf("rate identical %s/%s; NULL sets agree %s", chr(paste0(t, "_pieces_rate_identical")), chr(paste0(t, "_pieces_n")), chr(paste0(t, "_pieces_null_sets_agree")))))
  if (t == "binary") md <- c(md, row("F", t, "search-stage OR fit (glm.fit on vectors, SE from QR)",
                  chr("binary_F_n"), num("binary_F_formula_s"), num("binary_F_lean_s"),
                  sprintf("identical %s/%s; kinds agree %s", chr("binary_F_identical"), chr("binary_F_n"), chr("binary_F_kinds_agree"))))
  md <- c(md, row(if (t == "survival") "H (+B)" else if (t == "continuous") "H+E" else "H+G", t,
                  "whole consistency evaluator: lean extraction + lean pieces",
                  chr(paste0(t, "_stage_n")), num(paste0(t, "_stage_formula_s")), num(paste0(t, "_stage_lean_s")),
                  sprintf("Pcons identical: %s (passed %s vs %s; lean fallbacks %s)", chr(paste0(t, "_stage_pcons_identical_all")),
                          chr(paste0(t, "_stage_n_passed_pkg")), chr(paste0(t, "_stage_n_passed_lean")), chr(paste0(t, "_stage_n_lean_fallback")))))
  md <- c(md, row("H alone", t, if (t == "survival") "lean extraction + lean warm start, package pieces" else "lean extraction, package pieces",
                  chr(paste0(t, "_stage_n")), num(paste0(t, "_stage_formula_s")), num(paste0(t, "_stage_lean_Honly_s")),
                  sprintf("Pcons identical: %s", chr(paste0(t, "_stage_pcons_identical_H_only")))))
}
md <- c(md, "", sprintf("Continuous search-stage fit (already lean since 0.3.3): %s s per call over %s candidates -- no candidate.",
                        s(num("continuous_searchfit_fast_s"), 3), chr("continuous_searchfit_fast_n")))

# ---- combined per fixture ----------------------------------------------------------
md <- c(md, "", "## 3. Combined per fixture (per-replicate shape; machinery seconds are additive across stages)", "",
  "| fixture | replicate s | search-stage machinery s | consistency-stage machinery s (whole evaluator) | total s | share | clears 30%? |", "|---|---|---|---|---|---|---|")
for (t in names(FIXTURES)) {
  rep_s <- num(paste0(t, "_rep_s_mean"))
  srch <- if (t == "survival") A_mach_s else if (t == "binary") num("binary_F_machinery_s") else 0
  cons <- num(paste0(t, "_stage_machinery_s"))
  tot <- srch + cons; sh <- tot / rep_s
  md <- c(md, sprintf("| %s | %s | %s | %s | **%s** | **%s** | %s |", t, s(rep_s, 3), s(srch, 3), s(cons, 3), s(tot, 3), pct(sh),
                      if (is.na(sh)) "n/a" else if (sh >= BAR) "YES" else "no"))
}

md <- c(md, "", "The survival search-stage figure (A) is Stage 0's committed measurement on this machine (1.509 s of a 3.2112 s replicate); the survival replicate measured here is the denominator of the combined row.", "")

# ---- modification map --------------------------------------------------------------
st <- function(k) chr(k)
md <- c(md, "", "## 4. Modification map (CRAN 0.4.0 line references; nothing here was changed)", "",
  "| candidate | function(s) that would change | file:lines (0.4.0) | what replaces what | new helper | measured status |", "|---|---|---|---|---|---|",
  "| A | `fit_cox_for_subgroup()` unadjusted branch | R/subgroup_search.R:809-845 | formula `coxph()` -> `survival::coxph.fit()` with captured conventions; fallback to the formula block on error | `.fs_cox_lean_ci()` | Stage 0: 1410/1410 identical; 1.509 s/call |",
  sprintf("| B | `.consistency_cox_pieces()` (and the warm-start fits in `evaluate_subgroup_consistency()` / `evaluate_consistency_twostage()`) | R/consistency_resample.R:55-91; R/subgroup_consistency_helpers.R:1457-1467, 1827-1837 | two formula `coxph()` + `residuals(type = 'dfbeta')` with model-frame re-evaluation -> `coxph.fit()` x2 + `residuals.coxph()` on a constructed object | shares `.fs_cox_lean_ci()`'s fit | pieces rate identical %s/%s; machinery %s s/call |",
          st("survival_pieces_rate_identical"), st("survival_pieces_n"), s(num("survival_pieces_machinery_s"), 3)),
  sprintf("| E | `.consistency_glm_pieces()` MD branch + `.dfbeta_glm()` | R/consistency_resample.R:292, 168-176 | `lm()` + `lm.influence()` -> `lm.fit()` + hat values from the fit's QR | `.fs_glm_pieces_lean()` | pieces rate identical %s/%s; machinery %s s/call |",
          st("continuous_pieces_rate_identical"), st("continuous_pieces_n"), s(num("continuous_pieces_machinery_s"), 3)),
  sprintf("| F | routing predicate in `forestsearch()` + a fast OR closure + `fit_glm_for_subgroup_fast()` dispatch | R/forestsearch_main.R:2085-2104; R/glm_effect_estimators.R (new closure beside `.make_lm_estimator_fast`, :1058); R/subgroup_search.R:975-1003 | data-frame slice + `glm()` + `vcov()` -> vectors + `glm.fit()` + SE from the QR; existence check via the package's own `.fs_existence_reason()` | `.make_or_estimator_fast()`, `.fs_glm_two_group()` | search fits identical %s/%s; machinery %s s/call |",
          st("binary_F_identical"), st("binary_F_n"), s(num("binary_F_machinery_s"), 3)),
  sprintf("| G | `.consistency_glm_pieces()` OR branch + `.dfbeta_glm()` | R/consistency_resample.R:271-272, 168-176 | `glm()` + `lm.influence()` -> `glm.fit()` + hat values from the fit's QR | shares `.fs_glm_pieces_lean()` | pieces rate identical %s/%s; machinery %s s/call |",
          st("binary_pieces_rate_identical"), st("binary_pieces_n"), s(num("binary_pieces_machinery_s"), 3)),
  sprintf("| H | member extraction in `evaluate_consistency_twostage()` and `evaluate_subgroup_consistency()` | R/subgroup_consistency_helpers.R:1803-1811 and 1431-1439 (+ index row extraction :1778, :1391) | `as.data.table(df)` + sequential row filters per candidate -> one precomputed logical column cache and an index matrix, subset only the columns the pieces need | none (hoisting into `subgroup.consistency()`) | Pcons identical, H alone: survival %s, continuous %s, binary %s; machinery (H alone) %s / %s / %s s/call |",
          st("survival_stage_pcons_identical_H_only"), st("continuous_stage_pcons_identical_H_only"), st("binary_stage_pcons_identical_H_only"),
          s(num("survival_stage_H_only_machinery_s"), 3), s(num("continuous_stage_H_only_machinery_s"), 3), s(num("binary_stage_H_only_machinery_s"), 3)),
  "| C, D | none | -- | closed by Stage 0 measurement | -- | -- |",
  "", "Count/IRR is not measured: its search fit and pieces are the binary code paths with `poisson(log)` and an offset; the binary figures are the proxy. RR/RD (multi-tier fallbacks), adjusted and propensity-weighted fits, DINA and GRF identifiers are outside every candidate above.",
  "", "## 5. Reading", "")
for (t in names(FIXTURES)) {
  rep_s <- num(paste0(t, "_rep_s_mean"))
  srch <- if (t == "survival") A_mach_s else if (t == "binary") num("binary_F_machinery_s") else 0
  cons <- num(paste0(t, "_stage_machinery_s"))
  md <- c(md, sprintf("- %s: replicate %s s; search-stage machinery %s s, consistency-stage machinery %s s; combined %s (bar %s: %s). Consistency-stage identity: %s.",
    t, s(rep_s, 2), s(srch, 2), s(cons, 2), pct((srch + cons) / rep_s), pct(BAR), if ((srch + cons) / rep_s >= BAR) "cleared" else "not cleared",
    chr(paste0(t, "_stage_pcons_identical_all"))))
}
md <- c(md, "- Shares depend on each fixture's candidate flow (the consistency stage's cost is proportional to the number of candidates evaluated under full evaluation); the fixture sources above say which configuration each row describes.", "")
writeLines(md, out_path("stage0b_summary.md"))
cat(md, sep = "\n")
cat("\nb3 complete\n")
```
