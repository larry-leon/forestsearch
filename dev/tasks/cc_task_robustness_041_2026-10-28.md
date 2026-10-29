# CC task — robustness for 0.4.1: survival floor, continuous IPTW weights fix, lean Cox core with probe

**File:** `dev/tasks/cc_task_robustness_041_2026-10-28.md` · **Issued:** 2026-10-28 by chat, commissioned by Larry (decision: items (1)–(3) of the robustness review in the efficiency-routes chat; recorded in `decision_efficiency_routes_0.4.1_2026-10-28.md`, chat-side)
**Repo:** `forestsearch` · **Start on `master`**, HEAD `606a78ea` ("commit html vignettes", 2026-10-28) or a descendant with no change under `R/`, `DESCRIPTION`, `NAMESPACE`, `tests/` since then (§1). CC creates **`feature/robustness-0.4.1`** from it; everything in this task goes on that branch. CC never fetches, pulls, pushes or merges; Larry runs `R CMD check` on the branch and merges in GitHub Desktop before the release task.
**Transport:** `~/Downloads` → `dev/tasks/` (this document, committed first and alone). The filename has no dots other than `.md`; locate by pattern if punctuation was altered.
**Audit basis:** `master` at `606a78ea`. Line numbers below are from that commit, and §1 pins the files by MD5; verify every anchor by content before editing and report any line that moved.

---

## ⚠ CATEGORY — three `R/`-adjacent changes, each called out

1. **`DESCRIPTION` floor `survival (>= 3.5-8)`.** No code change. The lean unadjusted search-stage Cox fit calls `survival::coxph.fit()` with the calling conventions verified on survival 3.5-8 → 3.8.x; on an older survival that call errors inside `try()` and every candidate silently becomes non-fittable. The floor turns that into an install-time refusal. **Changes nothing for any analysis.**
2. **Bug fix — continuous IPTW weights (`R/glm_effect_estimators.R:957`). CHANGES BEHAVIOUR in one unreachable-so-far branch.** Inside `.make_lm_estimator()`'s closure, `stats::lm(fmla, data = data_slice, weights = data_slice$sw)` never resolved `data_slice` (the formula comes from `.build_adjusted_formula()`, so `model.frame()` evaluates `weights` in `data` and then in that helper's frame), the `tryCatch` returned `estimate = NA`, and `fit_glm_for_subgroup()` turned every candidate into status 5: an IPTW continuous analysis returned an empty search with no error. This is the defect recorded, deliberately unfixed, in `dev/glm-continuous-sims/REPORT_fastpath_unadjusted_2026-09-02.md` (§ "A pre-existing defect surfaced by F6": all 1,975 candidate fits fail at 0.3.2 and 0.3.3, `method_used = "ols_failed"`, reproduced minimally outside the search) and handed on as a correctness item in the 2026-09-02 close-out. The statement sits behind `if (ps_adjust_method == "iptw" && "sw" %in% names(data_slice))`, which is entered only when `is.RCT = FALSE` (or `ps_method` supplied), `ps_adjust_method = "iptw"` and the outcome is continuous. **Every committed analysis is `is.RCT = TRUE`** (`ps_method` resolves to `"none"`, no `sw` column is attached, and the fast continuous closure — a different function — is used); the `else` branch (unweighted `lm()`) is not touched; the three other IPTW sites (survival `:322`, OR `:579`, count `:776`) build their formulas in the closure frame, work today, and are not touched. §3 reproduces the defect on the unmodified build before anything is edited, and §11 greps the repository for any committed use of the path.
3. **Lean Cox core extracted with a once-per-session probe (`R/subgroup_search.R`). MOVES CODE; output bit-identical.** The `coxph.fit()` block of `fit_cox_for_subgroup()` (L845–L890) becomes the internal helper `.fs_cox_two_group()` (the twin of `.fs_logit_two_group()`), `coxph()`'s `ties`/`nocenter` defaults are read from its formals once per session instead of being hard-coded (with the literal values as fallback), and a once-per-session probe `.fs_cox_lean_ok()` fits a 30-row near-tie fixture through both `coxph.fit()` and `coxph()`: if they are not `identical()` (or the lean call errors — a changed signature in a future survival), the search-stage fit uses the 0.4.0 formula block for the session and says so once. On a healthy survival the same statements run with the same values; the only behavioural difference is in the failure case, which today yields a silent empty search. Identity is re-proved three ways (§5 tests, §7 gating tests, §8 Stage C against Stage A).

- **Writes:** `DESCRIPTION`; `R/glm_effect_estimators.R`; `R/globals.R`; `R/forestsearch_main.R` (one comment line); `R/subgroup_search.R`; `tests/testthat/test-lean-fits.R` (additions); new `tests/testthat/test-iptw-continuous.R`; `NEWS.md`; this document; the record under `dev/reports/`. **Nothing else.** No `man/` change (new functions are `@noRd`; no existing roxygen header changes), no `NAMESPACE` change.
- **Installs:** `devtools::install(upgrade = FALSE)` from the working tree twice: at §2 (so Stage A measures `master` itself) and at §7 (the edited tree). The installed package at the end is the branch build; the report says so.
- **Check policy:** no `R CMD check`, no full test suite, no vignette build. CC runs `test-lean-fits.R` and `test-iptw-continuous.R` (gating) and `test-fs-family-report.R`, `test-vi-deprecation.R` (non-gating), each under `timeout 10m` with `NOT_CRAN = "true"`. Larry's `R CMD check --as-cran` on the branch is the release gate.
- **Compute (the go):** Stage A and Stage C each run, per fixture (survival gbsg, continuous ACTG175, binary ACTG175), one warm-up call, three timed sequential calls and one `nb_boots = 20` sequential bootstrap — about 3–4 minutes per stage on the post-lean build (`stageC_timing.csv` of 2026-10-26: bootstraps 43 / 39 / 57 s); the reproduce step and the four test files a few minutes more. Hard caps per section below. Nothing else may run on the machine during the timed runs.
- **Unattended.** Gates stop the task; never ask, never work around. Provenance is dirt-tolerant (pre-existing untracked files are never staged).

---

## 1. Provenance, branch and first commit — GATE

1. `git rev-parse --abbrev-ref HEAD` must be `master`. Record `git rev-parse HEAD` (full SHA) and `git status --porcelain` verbatim. Tracked modifications (lines not starting `??`) → stop and report. Untracked files are tolerated and never staged.
2. If HEAD ≠ `606a78ea…`: `git diff --quiet 606a78ea HEAD -- R/ DESCRIPTION NAMESPACE tests/` must exit 0 (nothing under study changed since the audit basis); otherwise stop and report the `--stat`.
3. `md5sum` must print (GATE):
   - `51b866688442e727c9aaf1f7b777dad9  R/subgroup_search.R`
   - `e2fca0031b2e50577e23bea48775fa79  R/glm_effect_estimators.R`
   - `16a7771d31a15ae6f93ff6f906db4636  R/forestsearch_main.R`
   - `e8afa0f667951a17349a52d6774e3560  R/globals.R`
   - `dfb518541601e277871ebe85757efaa2  DESCRIPTION`
   - `6ea0c98f967922b2cd7556d54073e139  NEWS.md`
   - `0c76102bf26b18933d33b48b40915468  tests/testthat/test-lean-fits.R`
   - `7fddbf6e7d8740762baa14a7b459a208  dev/reports/vi_deprecation_2026-10-26/vi_common.R`
   - `6897eb239ceb3bc5d0d7a62d518bca2b  dev/reports/lean_fits_2026-10-26/timing_common.R`
4. Record `packageVersion("survival")` (must be ≥ 3.5-8 — GATE, since the floor is being pinned on this machine), `packageDescription("forestsearch")[c("Version", "Built")]`, `R.version.string`.
5. `git checkout -b feature/robustness-0.4.1`.
6. Copy this document to `dev/tasks/cc_task_robustness_041_2026-10-28.md`; `git add` by path; commit alone ("Task document: robustness for 0.4.1"). Report the `~/Downloads` filename and the repo filename.
7. Record directory: `dev/reports/robustness_041_2026-10-28/`. Report: `dev/reports/REPORT_robustness_041_2026-10-28.md`.

## 2. Stage A — baseline on the unmodified tree (install first)

1. `devtools::install(upgrade = FALSE)` from the working tree (= `master`). Record `Built`.
2. Into `dev/reports/robustness_041_2026-10-28/` copy `dev/reports/lean_fits_2026-10-26/timing_common.R` (MD5 above) and make **one** addition inside `run_timing()`, immediately after the line `cat(sprintf("[%s] warm-up done (warnings = %d)\n", fx_name, length(warm$warnings)))`:
   ```r
    saveRDS(extract_run(warm, fx$args$id.name),
            file.path(OUT, sprintf("stage%s_%s_ii.rds", stage, fx_name)))
   ```
   Nothing else in it changes (it still sources `../vi_deprecation_2026-10-26/vi_common.R` **unchanged**, so fixtures, `variant_args()`, `run_collect()`, `extract_run()`, `run_boot()`, `extract_boot()` are the committed ones and objects compare like with like).
3. `stageA.R` (transplant of `dev/reports/lean_fits_2026-10-26/stageA.R`, minus its deprecation-record gate): print and save `packageDescription()` + `sessionInfo()` as `stageA_sessionInfo.rds`; `tab <- run_timing("A")`; print it. Then, **informational, not a gate**: compare, saved against saved and component by component with `identical()`, `stageA_<fixture>_ii.rds` against the committed `../lean_fits_2026-10-26/identF_<fixture>_ii.rds` and `stageA_<fixture>_boot.rds` against `../lean_fits_2026-10-26/stageC_<fixture>_boot.rds`, printing one `INFO: ... identical / DIFFERS` line per component. (`args_call_all` may legitimately differ — the `max.minutes` deprecation of 2026-10-27 changed what the argument echo carries; a result component that differs is a finding for the report, not a stop.)
4. Run from the record directory: `timeout 15m Rscript stageA.R 2>&1 | tee stageA_log.txt`. Outputs: `stageA_<fixture>_ii.rds`, `stageA_<fixture>_boot.rds` (three each), `stageA_timing.csv`, `stageA_sessionInfo.rds`, `stageA_log.txt`.
5. Commit `timing_common.R`, `stageA.R` and the outputs by explicit path ("Stage A baseline on master, robustness task"). Hard cap 15 min.

## 3. Reproduce the IPTW defect on the Stage-A build — GATE (premise)

`dev/reports/robustness_041_2026-10-28/reproduce_iptw.R`, run **before any edit**, from the repo root with `timeout 10m Rscript ... 2>&1 | tee reproduce_iptw_log.txt`:

```r
library(forestsearch)
# (a) closure level: an explicit sw column, the slow continuous closure built as
#     forestsearch() builds it for ps_adjust_method = "iptw"
set.seed(11)
n <- 200L
d <- data.frame(treat = rep(0:1, each = n / 2))
d$y  <- 1 + 0.5 * d$treat + stats::rnorm(n)
ps   <- stats::runif(n, 0.2, 0.8)
d$sw <- ifelse(d$treat == 1L, mean(d$treat) / ps, (1 - mean(d$treat)) / (1 - ps))
est <- make_effect_estimator(outcome_type = "continuous", treat.name = "treat",
                             outcome.name = "y", ps_adjust_method = "iptw")
r <- est(d)
cat(sprintf("closure: estimate = %s, se = %s, converged = %s, method_used = %s\n",
            format(r$estimate), format(r$se), r$converged, r$method_used))
cat(sprintf("%s: closure-level weighted fit returns NA (the defect is present)\n",
            if (is.na(r$estimate)) "PASS" else "FAIL"))
# (b) end to end: the same call shape as tests/testthat/test-vi-deprecation.R's
#     .run(), with is.RCT = FALSE, ps_method = "logistic", ps_adjust_method = "iptw";
#     numeric confounders for the propensity model.
source("tests/testthat/helper-synthetic-dgm.R")
df <- .make_continuous_data(N = 300L, MD_harm = 1.5, seed = 42L)
df$biomarker_hi <- as.numeric(as.character(df$biomarker_hi))
df$sex <- as.numeric(as.character(df$sex))
rec <- new.env(); rec$names <- NULL
suppressMessages(trace("fit_glm_for_subgroup",
  tracer = bquote(if (is.null(.(rec)$names)) assign("names", names(df_clean), envir = .(rec))),
  where = asNamespace("forestsearch"), print = FALSE))
fit <- suppressWarnings(suppressMessages(forestsearch(
  df.analysis = df, confounders.name = c("age", "biomarker", "biomarker_hi", "sex"),
  outcome.name = "y", event.name = "y", treat.name = "treat", id.name = "id",
  outcome_type = "continuous", effect_measure = "MD", adverse_outcome = TRUE,
  effect.threshold = 0.5, consistency.threshold = 0.0, pconsistency.threshold = 0.5,
  fs.splits = 20, n.min = 30, maxk = 2, use_lasso = FALSE, use_grf = FALSE,
  sg_focus = "eff", is.RCT = FALSE, ps_method = "logistic", ps_adjust_method = "iptw",
  details = FALSE, quiet = TRUE, seedit = 1L, parallel_args = list(plan = "sequential"))))
suppressMessages(untrace("fit_glm_for_subgroup", where = asNamespace("forestsearch")))
tab <- fit$find.grps$out.found$hr.subgroups
cat(sprintf("end to end: sw in df_clean = %s; fitted candidates = %d; sg.harm = %s\n",
            "sw" %in% rec$names, if (is.null(tab)) 0L else nrow(tab),
            paste(fit$sg.harm, collapse = " & ")))
```

**GATE:** the closure-level line must print `PASS` (the defect is present on the unmodified build). If it prints `FAIL`, the premise of item 2 is wrong: stop and report. The end-to-end line is recorded, not gated; expected `sw in df_clean = TRUE` and `fitted candidates = 0`. If `sw in df_clean = FALSE`, also stop and report: the fix would then need to carry `sw` into the search frame, which is a different change and Larry's decision.

Commit `reproduce_iptw.R` and `reproduce_iptw_log.txt` by path ("Reproduce the continuous IPTW weights defect on master"). Hard cap 10 min.

## 4. Edit A — the IPTW fix

**`R/glm_effect_estimators.R`**, inside `.make_lm_estimator()`'s closure, L956–L957:

```r
      if (ps_adjust_method == "iptw" && "sw" %in% names(data_slice)) {
        fit <- stats::lm(fmla, data = data_slice, weights = data_slice$sw)
```
→
```r
      if (ps_adjust_method == "iptw" && "sw" %in% names(data_slice)) {
        # `weights` is evaluated by model.frame() in `data` and then in the
        # formula's environment.  `fmla` comes from .build_adjusted_formula(),
        # whose frame does not hold `data_slice`, so `weights = data_slice$sw`
        # never resolved and every weighted fit failed (0.4.0 and earlier).
        # `sw` is a column of data_slice (the guard above), so it resolves in
        # `data`.
        fit <- stats::lm(fmla, data = data_slice, weights = sw)
```
Nothing else in the file changes (`git diff --numstat` for this file must show 7 insertions, 1 deletion).

**`R/globals.R`**: the `utils::globalVariables(c(` vector ends at L266–L267:
```r
  "b", "r", "cell", "estimator", "cov", "cov1", "cov2", "obs", "ref"
))
```
→
```r
  "b", "r", "cell", "estimator", "cov", "cov1", "cov2", "obs", "ref",

# ============================================================================
# Propensity-weighted estimator closures: `weights = sw` is evaluated by
# model.frame() in the data slice (.make_lm_estimator())
# ============================================================================
  "sw"
))
```

**`R/forestsearch_main.R`** L2576 (a comment in the PS block header): `  #   ps_adjust_method = "iptw" (default)` → `  #   ps_adjust_method = "none" (default; "iptw" or "dr_gcomp" on request)`. The only change in this file.

**New `tests/testthat/test-iptw-continuous.R`** (gating):

```r
# =============================================================================
# Continuous IPTW weighted fit (regression, 0.4.0.9000): the weights did not
# resolve inside .make_lm_estimator()'s closure, so every candidate failed.
# Closure-level identity against lm(weights = ...), and the path end to end.
# =============================================================================

test_that("continuous IPTW closure reproduces the weighted OLS fit", {
  set.seed(11)
  n <- 200L
  d <- data.frame(treat = rep(0:1, each = n / 2))
  d$y  <- 1 + 0.5 * d$treat + stats::rnorm(n)
  ps   <- stats::runif(n, 0.2, 0.8)
  d$sw <- ifelse(d$treat == 1L, mean(d$treat) / ps, (1 - mean(d$treat)) / (1 - ps))

  est <- make_effect_estimator(outcome_type = "continuous", treat.name = "treat",
                               outcome.name = "y", ps_adjust_method = "iptw")
  r   <- est(d)
  ref <- stats::lm(y ~ treat, data = d, weights = d$sw)   # d is visible here
  expect_true(isTRUE(r$converged))
  expect_identical(r$method_used, "ols")
  expect_identical(r$estimate, stats::coef(ref)[["treat"]])
  expect_identical(r$se, sqrt(diag(stats::vcov(ref)))[["treat"]])

  # Without an sw column the same closure is the unweighted fit (unchanged).
  d0 <- d; d0$sw <- NULL
  r0 <- est(d0)
  ref0 <- stats::lm(y ~ treat, data = d0)
  expect_identical(r0$estimate, stats::coef(ref0)[["treat"]])
  expect_identical(r0$se, sqrt(diag(stats::vcov(ref0)))[["treat"]])
})

test_that("continuous IPTW path is reachable end to end and yields fitted candidates", {
  skip_on_cran()   # the closure-level test above is the regression guard everywhere
  df <- .make_continuous_data(N = 300L, MD_harm = 1.5, seed = 42L)
  df$biomarker_hi <- as.numeric(as.character(df$biomarker_hi))  # numeric confounders
  df$sex <- as.numeric(as.character(df$sex))                    # for the PS model
  fit <- suppressWarnings(suppressMessages(forestsearch(
    df.analysis = df, confounders.name = c("age", "biomarker", "biomarker_hi", "sex"),
    outcome.name = "y", event.name = "y", treat.name = "treat", id.name = "id",
    outcome_type = "continuous", effect_measure = "MD", adverse_outcome = TRUE,
    effect.threshold = 0.5, consistency.threshold = 0.0, pconsistency.threshold = 0.5,
    fs.splits = 20, n.min = 30, maxk = 2, use_lasso = FALSE, use_grf = FALSE,
    sg_focus = "eff", is.RCT = FALSE, ps_method = "logistic", ps_adjust_method = "iptw",
    details = FALSE, quiet = TRUE, seedit = 1L, parallel_args = list(plan = "sequential"))))
  tab <- fit$find.grps$out.found$hr.subgroups
  expect_true(!is.null(tab) && nrow(tab) > 0L)
})
```

Commit `R/glm_effect_estimators.R`, `R/globals.R`, `R/forestsearch_main.R`, `tests/testthat/test-iptw-continuous.R` by path ("Fix continuous IPTW weights resolution in .make_lm_estimator()"). (The test is run at §7, after the install.)

## 5. Edit B — the lean Cox core and the probe (`R/subgroup_search.R`)

**Insert**, immediately before L801 `#' Fit Cox Model for Subgroup` (and after the closing `}` of `meets_event_criteria()` at L799, keeping one blank line each side):

```r
# ---------------------------------------------------------------------------
# Lean unadjusted Cox core (0.4.0.9000)
# ---------------------------------------------------------------------------

# Per-session state for the lean core: the probe verdict and the calling
# conventions read from coxph().  One environment per R process, so every
# worker probes once.
.fs_cox_lean_state <- new.env(parent = emptyenv())

#' Calling conventions of `coxph()` that the lean core reproduces
#'
#' `coxph()` passes `method = match.arg(ties)` and its own `nocenter` default
#' to `coxph.fit()`.  Both are read from `coxph()`'s formals once per session,
#' so the lean core follows survival if either default moves; the literal
#' values of survival 3.5-8 .. 3.8.x (`"efron"`, `c(-1, 0, 1)`) are the
#' fallback when the formals cannot be read.  Whatever is read here is
#' checked against `coxph()` by `.fs_cox_lean_ok()` before the lean core is
#' used.
#'
#' @return List with `ties` (character scalar) and `nocenter` (numeric).
#' @keywords internal
#' @noRd
.fs_cox_lean_conventions <- function() {
  cv <- .fs_cox_lean_state$conventions
  if (!is.null(cv)) return(cv)
  f <- formals(survival::coxph)
  ties <- tryCatch(eval(f$ties)[1L], error = function(e) NULL)
  if (!is.character(ties) || length(ties) != 1L || is.na(ties)) ties <- "efron"
  nocenter <- tryCatch(eval(f$nocenter), error = function(e) NULL)
  if (!is.numeric(nocenter) || length(nocenter) == 0L) nocenter <- c(-1, 0, 1)
  cv <- list(ties = ties, nocenter = nocenter)
  .fs_cox_lean_state$conventions <- cv
  cv
}

#' Two-group Cox fit through `survival::coxph.fit()`, bit-identical to `coxph()`
#'
#' The fitting-engine call `coxph()` makes for `Surv(Y, E) ~ Treat`
#' (survival 3.5-8 .. 3.8.x): the one-column model matrix without intercept,
#' the response through `aeqSurv()` when `control$timefix`, no strata, a zero
#' offset, `init = NULL`, the default control, no weights, `coxph()`'s ties
#' method and `nocenter` (so the 0/1 treatment column is not centred, exactly
#' as in the formula fit), the model-frame row names, residuals on.
#' `coxph()` additionally builds the model frame, the full `coxph` object and
#' the concordance, none of which the search reads; `coefficients` and `var`
#' are returned by `coxph.fit()` unchanged, and the same warnings are raised
#' (monotone likelihood, non-convergence) -- nothing is suppressed.  Verified
#' bit-identical on every candidate of the gbsg application
#' (`dev/reports/REPORT_bootstrap_stage0_2026-10-26.md`, 1410/1410) and in
#' `tests/testthat/test-lean-fits.R`.  Guarded by `.fs_cox_lean_ok()`.
#'
#' @param Y,E,Treat Numeric vectors: time, event indicator and 0/1 treatment,
#'   already restricted to the candidate subgroup.
#' @return The `coxph.fit()` result (a list with `coefficients`, `var`, ...).
#' @keywords internal
#' @noRd
.fs_cox_two_group <- function(Y, E, Treat) {
  n_sg <- length(Y)
  cv   <- .fs_cox_lean_conventions()
  ctl  <- survival::coxph.control()
  ys   <- survival::Surv(Y, E)
  if (isTRUE(ctl$timefix)) ys <- survival::aeqSurv(ys)
  survival::coxph.fit(
    x        = matrix(as.double(Treat), ncol = 1L,
                      dimnames = list(NULL, "Treat")),
    y        = ys,
    strata   = NULL,
    offset   = rep(0, n_sg),
    init     = NULL,
    control  = ctl,
    weights  = NULL,
    method   = cv$ties,
    rownames = as.character(seq_len(n_sg)),
    resid    = TRUE,
    nocenter = cv$nocenter
  )
}

#' Once-per-session probe: does the lean core reproduce `coxph()`?
#'
#' Fits a fixed 30-row fixture both ways -- ten clusters of three near-tied
#' times that `aeqSurv()` collapses to exact ties (the tripwire for the
#' `timefix` convention), both arms present, events in both arms -- and
#' compares `coefficients` and `var` with `identical()`.  An error (for
#' example a changed `coxph.fit()` signature in a future survival) or any
#' difference turns the lean core off for the session: the search-stage fit
#' then goes through the formula interface, which is the same estimator, and
#' one message says so.  The verdict is cached in `.fs_cox_lean_state`.
#'
#' @return `TRUE` if the lean core may be used in this session.
#' @keywords internal
#' @noRd
.fs_cox_lean_ok <- function() {
  ok <- .fs_cox_lean_state$ok
  if (!is.null(ok)) return(ok)
  ok <- isTRUE(tryCatch({
    Y <- rep(seq(5, 50, by = 5), each = 3) + rep(c(0, 1e-9, 2e-9), 10)
    E <- rep(c(1L, 1L, 0L), 10)
    Treat <- rep(c(0, 1), 15)
    lean <- .fs_cox_two_group(Y, E, Treat)
    ref  <- survival::coxph(survival::Surv(Y, E) ~ Treat, robust = FALSE)
    identical(unname(lean$coefficients), unname(ref$coefficients)) &&
      identical(unname(as.vector(lean$var)), unname(as.vector(ref$var)))
  }, error = function(e) FALSE))
  if (!ok) {
    message("forestsearch: survival::coxph.fit() did not reproduce coxph() ",
            "on the probe fixture (survival ",
            as.character(utils::packageVersion("survival")),
            "); the unadjusted search-stage Cox fit uses the formula ",
            "interface in this session.")
  }
  .fs_cox_lean_state$ok <- ok
  ok
}

```

**Replace** L845–L890 of `fit_cox_for_subgroup()` (from `    hr.cox <- try({` through `    }, silent = TRUE)`, the whole lean block including its comment) with:

```r
    hr.cox <- try({
      # Unadjusted fit (0.4.0.9000): the lean core .fs_cox_two_group() --
      # survival::coxph.fit() with exactly the arguments coxph() passes for
      # Surv(Y, E) ~ Treat -- once the per-session probe .fs_cox_lean_ok()
      # has confirmed it reproduces coxph(); otherwise the 0.4.0 formula
      # fit.  The same estimator either way; coefficients and var are the
      # only fields consumed below.  Identity: tests/testthat/test-lean-fits.R.
      fit.cox <- if (.fs_cox_lean_ok()) {
        .fs_cox_two_group(Y, E, Treat)
      } else {
        survival::coxph(survival::Surv(Y, E) ~ Treat, robust = FALSE)
      }
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
```

(The eight lines from `beta <- fit.cox$coefficients` to `ci` are the current lines verbatim. In the formula branch `fit.cox$coefficients` is already named `"Treat"`, so `names(beta) <- "Treat"` is a no-op there and the `ci` matrix equals the 0.4.0 one.) The comment block above the branch (L820–L841, "Unadjusted fast branch (0.3.3) …") stays as it is; the formals, roxygen header, the three subset locals, the adjusted branch and everything after stay byte-identical. `@importFrom survival coxph Surv survfit` stays (the fallback and the adjusted branch use `coxph()`).

**Append to `tests/testthat/test-lean-fits.R`** (at the end of the file):

```r

test_that("A: the lean-core probe passes on the installed survival, and the formula fallback is identical", {
  state  <- forestsearch:::.fs_cox_lean_state
  old_ok <- state$ok
  old_cv <- state$conventions
  on.exit({ state$ok <- old_ok; state$conventions <- old_cv }, add = TRUE)

  # A fresh probe on this survival: TRUE and silent; conventions well-formed.
  state$ok <- NULL; state$conventions <- NULL
  expect_silent(ok <- forestsearch:::.fs_cox_lean_ok())
  expect_true(ok)
  cv <- forestsearch:::.fs_cox_lean_conventions()
  expect_true(is.character(cv$ties) && length(cv$ties) == 1L)
  expect_true(is.numeric(cv$nocenter) && length(cv$nocenter) >= 1L)

  # The core against coxph() on the probe fixture.
  Y_p <- rep(seq(5, 50, by = 5), each = 3) + rep(c(0, 1e-9, 2e-9), 10)
  E_p <- rep(c(1L, 1L, 0L), 10)
  T_p <- rep(c(0, 1), 15)
  lean <- forestsearch:::.fs_cox_two_group(Y_p, E_p, T_p)
  ref  <- survival::coxph(survival::Surv(Y_p, E_p) ~ T_p, robust = FALSE)
  expect_identical(unname(lean$coefficients), unname(ref$coefficients))
  expect_identical(unname(as.vector(lean$var)), unname(as.vector(ref$var)))

  # Forced fallback: fit_cox_for_subgroup() through the formula branch is
  # still identical to the 0.4.0 reference (gbsg subsets and the near-tie set).
  state$ok <- FALSE
  g <- survival::gbsg
  Yg <- g$rfstime / 30.4375; Eg <- g$status; Tg <- g$hormon
  set.seed(3)
  for (i in 1:5) {
    idx <- sort(sample(nrow(g), sample(40:300, 1L)))
    .lf_expect_cox_same(Yg[idx], Eg[idx], Tg[idx], sprintf("fallback subset %d", i))
  }
  .lf_expect_cox_same(Y_p, Eg[1:30], Tg[1:30], "fallback near-tie set")
})
```

Commit `R/subgroup_search.R` and `tests/testthat/test-lean-fits.R` by path ("Lean Cox core .fs_cox_two_group() with once-per-session probe and formula fallback").

## 6. `devtools::document()`, the `DESCRIPTION` floor, install

1. `devtools::document()`. Then `git status --porcelain` must show **no** `man/` and **no** `NAMESPACE` change. If `DESCRIPTION` shows only the `Config/roxygen2/version` re-stamp, `git checkout -- DESCRIPTION` (expected, not a deviation); if it shows anything else, stop and report.
2. `DESCRIPTION` L53 `    survival,` → `    survival (>= 3.5-8),` (the one edit to this file; `git diff --numstat` must show 1 1).
3. `devtools::install(upgrade = FALSE)`; record `packageDescription("forestsearch")[c("Version", "Built")]` and confirm the installed `DESCRIPTION` carries `survival (>= 3.5-8)` (`packageDescription("forestsearch")$Imports`).
4. Informational: `news(grepl("nocenter|resid", Text, ignore.case = TRUE), package = "survival")` — record the earliest version mentioning `nocenter` and the `coxph.fit` `resid` argument, if the NEWS database is readable; otherwise say so. (The floor stays 3.5-8, the verified range's lower end, whatever this says.)

Commit `DESCRIPTION` by path ("DESCRIPTION: survival (>= 3.5-8)").

## 7. Gating tests

From the repo root, each `timeout 10m`, with the console `tee`'d into the record directory:

```sh
timeout 10m Rscript -e 'Sys.setenv(NOT_CRAN = "true"); testthat::test_file("tests/testthat/test-lean-fits.R", package = "forestsearch")'          2>&1 | tee dev/reports/robustness_041_2026-10-28/tests_gating_lean_fits.txt
timeout 10m Rscript -e 'Sys.setenv(NOT_CRAN = "true"); testthat::test_file("tests/testthat/test-iptw-continuous.R", package = "forestsearch")'    2>&1 | tee dev/reports/robustness_041_2026-10-28/tests_gating_iptw.txt
timeout 10m Rscript -e 'Sys.setenv(NOT_CRAN = "true"); testthat::test_file("tests/testthat/test-fs-family-report.R", package = "forestsearch")'    2>&1 | tee dev/reports/robustness_041_2026-10-28/tests_nongating_family_report.txt
timeout 10m Rscript -e 'Sys.setenv(NOT_CRAN = "true"); testthat::test_file("tests/testthat/test-vi-deprecation.R", package = "forestsearch")'      2>&1 | tee dev/reports/robustness_041_2026-10-28/tests_nongating_vi_deprecation.txt
```

**GATE:** the two gating files report `FAIL 0` with `PASS > 0`. The two non-gating files are recorded (their `TOTALS` lines go in the report). If `test_file(..., package = )` does not load the helpers on this testthat, use the invocation recorded for the lean-fits task (`dev/reports/lean_fits_2026-10-26/tests_gating_log.txt`) and say so.

## 8. Stage C — identity against Stage A (GATE), records, timing

`stageC.R` in the record directory (transplant of `dev/reports/lean_fits_2026-10-26/stageC.R`: same `run_timing("C")` protocol via the §2 `timing_common.R`; drop its sections-3/4 copy-in and its prediction table; keep the gate/summary mechanics). It writes `stageC_summary.txt` with:

1. **GATE — end to end:** for each fixture, `stageC_<fixture>_ii.rds` vs `stageA_<fixture>_ii.rds`, saved against saved, `identical()` per component (`sg.harm`, `membership`, `out_sg_result`, `hr.subgroups`, `filter_counts`, `n_candidates`, `n_fits`, `n_warnings`, `warnings`, `args_call_all`) — one `PASS/FAIL` line each.
2. **GATE — bootstrap:** `stageC_<fixture>_boot.rds` vs `stageA_<fixture>_boot.rds`, per component as the lean-fits `stageC.R` does.
3. **GATE — tests:** the `TOTALS` of the two gating logs (`FAIL 0`, `PASS > 0`); the two non-gating `TOTALS` as `INFO` lines.
4. **GATE — reproduce:** the `PASS` line of `reproduce_iptw_log.txt` copied in; and, on the new build, re-run the §3 closure-level block (source the same script's part (a) or re-execute it inline) and assert `!is.na(r$estimate)` — `PASS: closure-level weighted fit now returns an estimate`.
5. **INFO — timing:** per fixture, call median and replicate mean, Stage A → Stage C, with the percentage; expected within noise (the probe adds one 30-row fit per session). A difference is a finding, not a gate.
6. **INFO — probe:** `forestsearch:::.fs_cox_lean_state$ok` after the run (`TRUE`) and `.fs_cox_lean_conventions()` (`ties`, `nocenter`) on this survival.

Run `timeout 15m Rscript stageC.R 2>&1 | tee stageC_log.txt` from the record directory. **GATE:** `grep -c FAIL stageC_summary.txt` is 0. Commit `stageC.R`, the outputs, the four test logs and `stageC_summary.txt` by path ("Stage C: identity against Stage A, robustness task"). Hard cap 15 min.

## 9. `NEWS.md`

Under `# forestsearch (development version)`:

- `### Internal` (after the `fit_causal_forest_glm()` bullet):
  ```
  * `survival (>= 3.5-8)` is now required: the lean unadjusted search-stage
    Cox fit calls `survival::coxph.fit()` with the calling conventions
    verified on survival 3.5-8 and later.
  * The lean unadjusted Cox fit lives in the internal helper
    `.fs_cox_two_group()` and is guarded by a once-per-session probe: if
    `coxph.fit()` does not reproduce `coxph()` on a near-tie fixture (for
    example after a change in a future survival), the search-stage fit falls
    back to the formula interface for the session, with one message. Results
    are identical either way (`tests/testthat/test-lean-fits.R`; record
    `dev/reports/REPORT_robustness_041_2026-10-28.md`).
  ```
- `### Bug fixes` (after the `xlog` bullet):
  ```
  * `ps_adjust_method = "iptw"` with a continuous outcome: the weighted
    least-squares fit failed for every candidate, because the `weights`
    argument did not resolve inside the estimator closure, so the search
    returned no candidates without an error. The weights now resolve from the
    data slice. The unadjusted and covariate-adjusted continuous fits, and the
    binary, count and survival paths, are unchanged; no committed analysis
    used this path (all are randomized-trial analyses with `is.RCT = TRUE`).
    Regression test `tests/testthat/test-iptw-continuous.R`.
  ```

Commit `NEWS.md` by path ("NEWS: survival floor, lean-core probe, continuous IPTW fix").

## 10. Report — `dev/reports/REPORT_robustness_041_2026-10-28.md`

1. Provenance: starting HEAD, branch, every commit SHA in order, installed build before (§2) and after (§6), `packageVersion("survival")`, R version.
2. The edits, file by file: `git diff --stat 606a78ea HEAD -- R/ DESCRIPTION tests/ NEWS.md` and the `@@` hunk headers of each `R/` file (each hunk must sit inside `.make_lm_estimator()`'s closure, the `globalVariables` vector, the PS-block comment, the new helpers, or `fit_cox_for_subgroup()`'s unadjusted branch — say so).
3. Reproduce (§3): the two recorded lines before the fix; the closure-level line after it (§8 item 4).
4. Identity: the Stage C vs Stage A component tables (end to end and bootstrap), and the Stage A vs committed-record `INFO` lines (§2.3) with any differing component named.
5. Tests: the four `TOTALS` lines, gating and non-gating separately.
6. Timing: the Stage A → Stage C table (informational).
7. **Behaviour statement for Larry's record**, written from the evidence above: (1) no code runs differently under the floor; (2) the IPTW statement is unreachable in every committed analysis (`is.RCT = TRUE` throughout; the §11 grep result quoted), the `else` branch is byte-identical, and the fast continuous closure is untouched; (3) the lean core executes the same statements with the same values on this survival (probe `TRUE`), and Stage C equals Stage A on every component of every fixture and bootstrap.
8. Script deviations (verbatim if any). OPEN ITEMS, at least: (a) the consistency-stage pieces (`.consistency_glm_pieces()`, `R/consistency_resample.R:239–`) take no weights for any `ps_adjust_method`, so an IPTW analysis is weighted in the search stage and unweighted in the consistency screen — a method question for Larry, not changed here; (b) at `is.RCT = FALSE` with the defaults (`ps_method = "grf"`, `ps_adjust_method = "none"`) a propensity forest is fitted once per call and once per bootstrap replicate and never consumed, and the rebuilt estimator closure is untagged, so the lean continuous/OR path is lost (`R/forestsearch_main.R:2626–2697`); (c) `subgroup_glm()` in `R/run_subgroup_sims.R:331–340` still routes only the continuous case to a fast closure. Ten-line verdict.
9. Closing block: the §11 commands and their output.

## 11. Commit and post-conditions — GATE

Allowed paths: `dev/tasks/cc_task_robustness_041_2026-10-28.md`; `DESCRIPTION`; `R/glm_effect_estimators.R`; `R/globals.R`; `R/forestsearch_main.R`; `R/subgroup_search.R`; `tests/testthat/test-lean-fits.R`; `tests/testthat/test-iptw-continuous.R`; `NEWS.md`; `dev/reports/REPORT_robustness_041_2026-10-28.md`; `dev/reports/robustness_041_2026-10-28/`. Before every commit `git diff --cached --name-only` must be a subset of that list. After the final commit, print in the closing block:

```
git diff --name-only 606a78ea HEAD | grep -v -E '^(dev/|NEWS\.md$|DESCRIPTION$|R/(subgroup_search|glm_effect_estimators|forestsearch_main|globals)\.R$|tests/testthat/(test-lean-fits|test-iptw-continuous)\.R$)'
    # must print nothing (in particular: no man/, no NAMESPACE)
git diff 606a78ea HEAD -- DESCRIPTION | grep -E '^[-+] '
    # must print exactly:  -    survival,   and   +    survival (>= 3.5-8),
git diff 606a78ea HEAD -- R/forestsearch_main.R | grep -E '^[-+] '
    # must print exactly the two comment lines of section 4
git diff --numstat 606a78ea HEAD -- R/glm_effect_estimators.R
    # must print 7 1 (the comment lines and the one changed statement)
grep -c 'weights = data_slice\$sw' R/glm_effect_estimators.R
    # must print 3 (survival, OR, count sites untouched; the continuous site is gone)
grep -c 'weights = sw)' R/glm_effect_estimators.R
    # must print 1
grep -n "coxph.fit(" R/subgroup_search.R | grep -v -E '^[0-9]+:[[:space:]]*#' | wc -l
    # must print 1 (inside .fs_cox_two_group())
grep -n "survival::coxph(" R/subgroup_search.R | grep -v -E '^[0-9]+:[[:space:]]*#' | wc -l
    # must print 3 (probe reference, formula fallback, adjusted branch)
grep -n "aeqSurv\|coxph.control" R/subgroup_search.R | grep -v -E '^[0-9]+:[[:space:]]*#' | wc -l
    # must print 2 (both inside .fs_cox_two_group())
grep -ln "fs_cox_lean_state\|fs_cox_lean_conventions\|fs_cox_two_group\|fs_cox_lean_ok" R/*.R
    # must print only R/subgroup_search.R
grep -rn 'ps_adjust_method *= *"iptw"' --include=*.qmd --include=*.Rmd --include=*.R quarto/applications vignettes tests inst 2>/dev/null
    # record the output: expected empty for the application documents, vignettes, tests and inst (quarto/grf/ uses the binary OR path, which is untouched)
grep -c FAIL dev/reports/robustness_041_2026-10-28/stageC_summary.txt
    # must print 0
git status --porcelain | grep -v '^?? '
    # must print nothing
```

If any command prints something unexpected, stop and report what it prints.

Finally, copy `dev/reports/REPORT_robustness_041_2026-10-28.md`, `dev/reports/robustness_041_2026-10-28/stageC_summary.txt`, `dev/reports/robustness_041_2026-10-28/stageC_timing.csv` and `dev/reports/robustness_041_2026-10-28/reproduce_iptw_log.txt` to `~/Downloads` and print `ls -la` of the four copies.

## 12. Out of scope

- The consistency-stage items B/E/G/H and every other efficiency route (closed for 0.4.1; decision memo).
- Weights in the consistency-stage pieces; the observational-path defaults; `subgroup_glm()`'s binary routing — OPEN ITEMS in the report, Larry's decisions.
- The three other IPTW sites (survival, OR, count): they work today and are not touched.
- Adjusted, weighted, stratified, RD and RR search fits: unchanged.
- No `R CMD check`, no full test suite, no vignette build, no version bump, no push, no merge.
