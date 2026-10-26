# CC task — lean search-stage fits: A (survival, `coxph.fit`) and F (binary OR, `glm.fit`), bit-identical

**File:** `dev/tasks/cc_task_lean_fits_2026-10-26.md` · **Issued:** 2026-10-26 by chat, commissioned by Larry (decision: proceed with A + F after Stage 0b, `dev/reports/REPORT_bootstrap_stage0b_2026-10-26.md`, and after the deprecation task on this branch, `dev/reports/REPORT_vi_deprecation_2026-10-26.md`)
**Repo:** `forestsearch` · **Branch:** `feature/vi-deprecation-lean-fits`, already checked out, HEAD `60a6496b` (the deprecation closing-block commit). Everything in this task goes on it. CC never fetches, pulls, pushes or merges; Larry merges in GitHub Desktop after his own `R CMD check`.
**Transport:** `~/Downloads` → `dev/tasks/` (this document, committed first and alone). Filenames may arrive with hyphens stripped; locate by pattern.
**Audit basis:** CRAN 0.4.0 line numbers for `R/subgroup_search.R` (MD5 f06d6bba…) and `R/glm_effect_estimators.R` (MD5 8e8239f9…) — neither was touched by the deprecation task, so those are the current numbers; `R/forestsearch_main.R` was edited by that task, so its one anchor below is given by content. Verify every named line by content before editing; if a line has moved, find it by content and say so in the report.

---

## ⚠ CATEGORY — `R/` change: MOVES EXISTING CODE TO A LEANER EQUIVALENT, OUTPUT BIT-IDENTICAL; no behaviour change, no method change

- **What changes.** Two per-candidate search-stage fits stop going through the formula interface and call the fitting engine directly with exactly the arguments the formula interface passes:
  - **A (survival):** `fit_cox_for_subgroup()`'s unadjusted branch calls `survival::coxph.fit()` instead of `survival::coxph(Surv(Y, E) ~ Treat, robust = FALSE)`. Same engine, same arguments, same `coefficients` and `var`; what is skipped is the model frame, the full `coxph` object and the concordance computation the search never reads.
  - **F (binary, `effect_measure = "OR"`, no covariate or propensity adjustment):** a tagged fast estimator closure (the twin of the existing continuous `.make_lm_estimator_fast()`) and a vector path in `fit_glm_for_subgroup_fast()` call `stats::glm.fit()` on the two-column model matrix instead of `stats::glm()` on a data-frame slice, with the standard error taken from the fit's QR exactly as `summary.glm()`/`vcov.glm()` compute it.
- **Evidence it is identical.** Stage 0 (A: 1410/1410 candidates of the gbsg survival application, `identical()` on the consumed CI row, warning sets equal) and Stage 0b (F: 2577/2577 candidates of the ACTG175 binary application, `identical()` on the returned list, failure kinds equal). Both measured against unmodified package code with the recipes below. This task re-proves it in three layers: per candidate on captured inputs, end to end on the three committed fixtures against the committed deprecation-task record, and in a permanent testthat file with edge fixtures.
- **Expected gain (Stage 0b, sequential replicate):** survival 3.14 s → about 1.6 s (A removes 1.51 s); binary 6.37 s → about 3.2 s (F removes 3.13 s); continuous unchanged (its search fit is already lean). Realized numbers are measured here (§2, §6) and reported against those predictions. **A shortfall is a finding, not a gate.**
- **Dependence introduced.** A hard-codes `coxph()`'s calling conventions for `coxph.fit()` (survival 3.5.8 → 3.8.x: `x` model matrix without intercept, `y` through `aeqSurv()` when `control$timefix`, `strata = NULL`, `offset = rep(0, n)`, `init = NULL`, `control = coxph.control()`, `weights = NULL`, `method = "efron"`, `rownames`, `resid = TRUE`, `nocenter = c(-1, 0, 1)`). `coxph.fit`, `aeqSurv` and `coxph.control` are exported; the package already calls `coxph.fit` (`.g3_cox_lean()`). The new testthat identity test is the tripwire if a future survival changes those conventions.
- **Writes:** `R/subgroup_search.R`, `R/glm_effect_estimators.R`, `R/forestsearch_main.R`; new `tests/testthat/test-lean-fits.R`; `NEWS.md`; this document; the battery and report under `dev/reports/`. **Nothing else.** No `man/` change (every new function is `@noRd`; `fit_cox_for_subgroup()`'s roxygen is not touched), no `NAMESPACE` change, no `DESCRIPTION` change (`devtools::document()` re-stamps `Config/roxygen2/version` — restore it with `git checkout -- DESCRIPTION` before installing, as the deprecation task did; that is expected, not a deviation).
- **Installs:** `devtools::install(upgrade = FALSE)` from the working tree twice (after A, after F). The installed package at the end is the A+F branch build; the report says so.
- **Check policy:** no `R CMD check`, no full test suite, no vignette build. CC runs the new test file and, non-gating, `test-fs-family-report.R` and `test-vi-deprecation.R` (they exercise `forestsearch()` end to end on the continuous path) — each with `Sys.setenv(NOT_CRAN = "true")`, under `timeout 10m`. Larry's `R CMD check` on the branch is the release gate.
- **Compute:** per fixture and stage, one warm-up call, three timed sequential `forestsearch()` calls and one `nb = 20` sequential bootstrap (≈ 4–5 min per stage over the three fixtures); two per-part identity passes (≈ 2 min each); three test files. Nothing else may run on the machine during the timed runs. Hard caps in §2, §4, §5.
- **Unattended.** Gates stop the task; never ask, never work around. Provenance is dirt-tolerant.

---

## 1. Provenance and first commit — GATE

1. `git rev-parse --abbrev-ref HEAD` must be `feature/vi-deprecation-lean-fits` and `git rev-parse HEAD` must be `60a6496b…` (record the full SHA; if HEAD differs, stop and report — Larry may have committed something). `git status --porcelain` verbatim (pre-existing untracked dirt tolerated, never staged).
2. `md5sum R/subgroup_search.R R/glm_effect_estimators.R` must print `f06d6bba…` and `8e8239f9…` (the CRAN 0.4.0 files). **GATE.**
3. `packageDescription("forestsearch")[c("Version", "Built")]` — expected the deprecation-task build (0.4.0.9000, built 2026-10-26 17:05:25 UTC). If it is anything else, `devtools::install(upgrade = FALSE)` from the current working tree first and say so: Stage A must be measured on the deprecation build.
4. Copy this document to `dev/tasks/cc_task_lean_fits_2026-10-26.md`, `git add` by path, commit alone. Report both filenames.
5. Record directory: `dev/reports/lean_fits_2026-10-26/`.

## 2. Stage A — timing baseline and bootstrap payloads on the deprecation build (no edits yet)

Identity baselines for `forestsearch()` itself are **already committed**: `dev/reports/vi_deprecation_2026-10-26/stageC_<fixture>_ii.rds` (the default-argument runs on this very build) and `stageC_survival_boot.rds`. Do not regenerate them. What Stage A adds is timing, and the binary and continuous bootstrap payloads that record lacks.

Script `dev/reports/lean_fits_2026-10-26/stageA.R`, which `source()`s `../vi_deprecation_2026-10-26/vi_common.R` **unchanged** (fixtures, `variant_args()`, `run_collect()`, `extract_run()`, `run_boot()`, `extract_boot()`, `count_grf_fits()` — the same code that built the committed record, so objects compare like with like). Use the `"ii"` argument variant for every fixture (no deprecated argument, no warning). For each fixture in order survival, continuous, binary:

1. one untimed warm-up `forestsearch()` call;
2. three timed sequential calls (`system.time()` elapsed each; keep all three, report the median);
3. one `nb_boots = 20` sequential bootstrap via `run_boot()` on the warm-up fit, timed as a whole (`elapsed`), with `mean(results$tmins_iteration) * 60` as the per-replicate seconds; save `extract_boot()` of it as `stageA_<fixture>_boot.rds` (the survival one must be `identical()` component by component to the committed `stageC_survival_boot.rds` — a sanity gate that the machine reproduces its own record; continuous and binary have no prior record and become the baseline).

Save `stageA_timing.csv` (fixture, call_s_1..3, call_median_s, boot_outer_s, rep_mean_s) and `sessionInfo()`. Commit script and outputs by path ("Stage A timing baseline, lean fits"). **GATE:** the survival bootstrap payload identity. Hard cap 20 min.

## 3. Part A — the survival search fit (`R/subgroup_search.R`)

Edit **only** the `hr.cox <- try({ ... }, silent = TRUE)` block of the unadjusted branch of `fit_cox_for_subgroup()`, 0.4.0 **L834–L845** (from `hr.cox <- try({` through `}, silent = TRUE)`), and one phrase of the comment block above it: L810 `# Unadjusted fast branch (0.3.3): same coxph() fit, but the consumed` → `# Unadjusted fast branch (0.3.3): the same Cox fit (through coxph.fit() since 0.4.0.9000, see below), but the consumed` (re-wrap to 80 columns). The rest of that comment (L811–L830, the summary.coxph arithmetic) still holds and stays. The formals, the roxygen header (L791–L801), the three subset locals (L831–L833), the adjusted branch and everything after L845 stay byte-identical.

Replace L834–L845 with:

```r
    hr.cox <- try({
      # Lean unadjusted fit (0.4.0.9000): survival::coxph.fit() called
      # directly with exactly the arguments coxph() itself passes for
      # Surv(Y, E) ~ Treat (survival 3.5.8 .. 3.8.x): the one-column model
      # matrix without intercept, the response through aeqSurv() when
      # control$timefix, no strata, a zero offset, init = NULL, the default
      # control, no weights, Efron ties, the model-frame row names,
      # residuals on, and coxph()'s nocenter = c(-1, 0, 1) (so the 0/1
      # treatment column is not centred, exactly as in the formula fit).
      # coxph() additionally builds the model frame, the full coxph object
      # and the concordance, none of which the search reads; coefficients
      # and var -- the only fields consumed below -- are returned by
      # coxph.fit() unchanged.  Verified bit-identical on every candidate
      # of the gbsg application (REPORT_bootstrap_stage0_2026-10-26.md,
      # 1410/1410) and re-verified in tests/testthat/test-lean-fits.R.
      # coxph.fit() raises the same warnings coxph() would (monotone
      # likelihood, non-convergence): nothing is suppressed.
      n_sg <- length(Y)
      ctl  <- survival::coxph.control()
      ys   <- survival::Surv(Y, E)
      if (isTRUE(ctl$timefix)) ys <- survival::aeqSurv(ys)
      fit.cox <- survival::coxph.fit(
        x        = matrix(as.double(Treat), ncol = 1L,
                          dimnames = list(NULL, "Treat")),
        y        = ys,
        strata   = NULL,
        offset   = rep(0, n_sg),
        init     = NULL,
        control  = ctl,
        weights  = NULL,
        method   = "efron",
        rownames = as.character(seq_len(n_sg)),
        resid    = TRUE,
        nocenter = c(-1, 0, 1)
      )
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

(The eight lines from `beta <- fit.cox$coefficients` to `ci` are the 0.4.0 lines verbatim, plus `names(beta) <- "Treat"`.) Do not change `@importFrom survival coxph Surv survfit` (the adjusted branch still uses `coxph()`).

Then `devtools::document()` → `git status --porcelain` must show **no** `man/` or `NAMESPACE` change (if `DESCRIPTION` shows the roxygen re-stamp, `git checkout -- DESCRIPTION`); `devtools::install(upgrade = FALSE)`; record `packageDescription("forestsearch")$Built`.

**Part-A identity pass** — `dev/reports/lean_fits_2026-10-26/identA.R` (sources `vi_common.R`), survival fixture, `"ii"` args:

- **per candidate:** `base::trace()` on `forestsearch:::fit_cox_for_subgroup` (session only, `untrace()` right after the call) records `list(Y = yy[id.x == 1], E = dd[id.x == 1], Treat = tt[id.x == 1], adj = forestsearch:::.fs_adjust_terms(adjust_covariates))` for every call during one `forestsearch()` run (expected 1410 calls, all unadjusted — the Stage-0 count; print it). Then for every captured candidate compare, with `identical()`, (a) the installed function `forestsearch:::fit_cox_for_subgroup(Y, E, Treat, rep(1L, length(Y)))` against (b) the reference below, which is the 0.4.0 formula block plus the 0.4.0 tail, transplanted verbatim; and compare the warning texts collected around each (`withCallingHandlers`, muffled). **GATE:** all candidates identical in value and in warning set.

  ```r
  ref_fit_cox <- function(Y, E, Treat) {           # CRAN 0.4.0 R/subgroup_search.R L834-L845, L873, L875-L881, L888-L894, verbatim
    hr.cox <- try({
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
    if (inherits(hr.cox, "try-error")) return(NULL)
    if (is.null(dim(hr.cox))) {
      hr.cox <- matrix(hr.cox, nrow = 1L, dimnames = list("Treat", names(hr.cox)))
    }
    trow <- if ("Treat" %in% rownames(hr.cox)) hr.cox["Treat", ] else hr.cox[1, ]
    list(hr = trow[1], lower = trow[3], upper = trow[4], med0 = NA_real_, med1 = NA_real_)
  }
  ```

- **end to end:** `run_collect()` + `extract_run()` on the survival `"ii"` args, saved as `identA_survival_ii.rds`, reloaded, and compared component by component (`identical()`, saved against saved) to the committed `../vi_deprecation_2026-10-26/stageC_survival_ii.rds`, **including the `warnings` component**. **GATE.**

Commit `R/subgroup_search.R`, `identA.R` and its outputs by explicit path ("Lean survival search fit: coxph.fit() in fit_cox_for_subgroup()"). Hard cap 15 min.

## 4. Part F — the binary OR search fit

**`R/glm_effect_estimators.R`** — append, immediately before the `# Utilities` banner block (0.4.0 L1101–L1103; `.make_lm_estimator_fast()` closes at L1098), as pure additions:

```r
# ---------------------------------------------------------------------------
# Unadjusted fast path, binary OR (0.4.0.9000)
# ---------------------------------------------------------------------------

#' Two-group logistic core, bit-identical to `glm()` + `coef()` + `vcov()`
#'
#' Fits `y ~ t` for an unadjusted, unweighted, offset-free two-group
#' comparison via `stats::glm.fit()` on the explicit two-column model matrix
#' -- the call `stats::glm()` makes for that formula (weights, start,
#' etastart, mustart and offset at their defaults, `intercept = TRUE`,
#' `singular.ok = TRUE`, `glm.control()`) -- then takes the treatment
#' coefficient and its standard error by replicating `stats::summary.glm()`
#' / `stats::vcov.glm()` operation for operation (R 4.6.1 sources):
#' `p1 <- 1L:p`; `covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])`;
#' `covmat <- dispersion * covmat.unscaled` with `dispersion = 1` for the
#' binomial family (kept as `1 * ...` so the arithmetic is the same);
#' `sqrt(diag(covmat))` at the treatment column.  A rank-deficient fit
#' (single treatment arm) returns `estimate = NA` and `se = NA` without
#' error, as the slow chain does.  `glm.fit()`'s own warnings ("fitted
#' probabilities numerically 0 or 1 occurred", "algorithm did not converge")
#' are raised exactly as through `glm()`: nothing is suppressed here.
#'
#' @param y Numeric vector. 0/1 outcome, already sliced, already on the
#'   adverse scale.
#' @param t Numeric or integer vector. Treatment indicator (0/1), same length.
#' @param treat.name Character. The treatment column name `model.matrix()`
#'   would use.
#' @return List with `estimate`, `se` (doubles, possibly `NA`) and
#'   `converged` (`glm.fit()`'s flag).
#' @keywords internal
#' @noRd
.fs_logit_two_group <- function(y, t, treat.name) {
  X <- cbind(1, t)
  colnames(X) <- c("(Intercept)", treat.name)
  fit <- stats::glm.fit(X, as.double(y),
                        family = stats::binomial(link = "logit"))
  coef_val <- fit$coefficients[[treat.name]]
  p1 <- seq_len(fit$rank)
  covmat <- 1 * chol2inv(fit$qr$qr[p1, p1, drop = FALSE])
  # Column 2 is the treatment term; with the treatment column aliased
  # (rank 1) the index is out of range and gives NA, as vcov(complete = TRUE)
  # does for an aliased coefficient.
  se_val <- unname(sqrt(diag(covmat))[2L])
  list(estimate = coef_val, se = se_val, converged = fit$converged)
}


#' Fast unadjusted binary odds-ratio estimator closure
#'
#' Drop-in replacement for the `.make_glm_binary_estimator()` closure when
#' `effect_measure = "OR"`, there is no covariate adjustment and no
#' propensity adjustment: the routing predicate at the construction site
#' (`forestsearch_main.R`) selects it, and every consumer of the closure
#' works unchanged because the signature and the returned structure --
#' fields, names, types, failure behaviour, the `reason` of a non-estimable
#' slice -- are identical to the slow closure's, and the numbers are
#' bit-identical (see `.fs_logit_two_group()`).  The search loop routes its
#' per-candidate fit to the vector path on the `fs_fast_unadjusted`
#' attribute; `fs_fast_core` in this closure's environment tells
#' `fit_glm_for_subgroup_fast()` which core to run (the continuous closure
#' carries no marker and keeps its OLS path).
#'
#' @keywords internal
#' @noRd
.make_glm_or_estimator_fast <- function(treat.name, outcome.name,
                                        adverse_outcome = TRUE, ...) {
  force(treat.name)
  force(outcome.name)
  force(adverse_outcome)
  fs_fast_core <- "logit_or"

  function(data_slice) {
    # Sign convention as in .make_glm_binary_estimator(): Y -> 1 - Y when the
    # outcome is beneficial, so that OR > 1 indicates treatment harm.
    if (!adverse_outcome) {
      data_slice[[outcome.name]] <- 1L - data_slice[[outcome.name]]
    }

    n0 <- sum(data_slice[[treat.name]] == 0L, na.rm = TRUE)
    n1 <- sum(data_slice[[treat.name]] == 1L, na.rm = TRUE)

    why <- .fs_existence_reason(
      "OR",
      .fs_binary_cells(data_slice[[outcome.name]], data_slice[[treat.name]]))
    if (!is.null(why)) {
      return(.fs_nonestimable(why, n0, n1, "log-OR", "nonestimable"))
    }

    result <- tryCatch({
      core <- .fs_logit_two_group(data_slice[[outcome.name]],
                                  data_slice[[treat.name]], treat.name)
      list(estimate = core$estimate, se = core$se, converged = core$converged,
           n0 = n0, n1 = n1, measure = "log-OR", method_used = "logistic")
    },
    error = function(e) {
      list(estimate = NA_real_, se = NA_real_, converged = FALSE,
           n0 = n0, n1 = n1, measure = "log-OR",
           method_used = "logistic_failed_none")
    })

    # The slow closure's non-convergence rule, verbatim.
    if (isFALSE(result$converged) && !is.na(result$estimate)) {
      return(.fs_nonestimable("non-convergent fit", n0, n1,
                              result$measure, result$method_used))
    }

    result
  }
}
```

**`R/subgroup_search.R`** — `fit_glm_for_subgroup_fast()` (0.4.0 L955–L1003). Keep its roxygen (it is `@noRd`), adding one sentence after the first paragraph: `The tagged unadjusted binary OR closure (\code{.make_glm_or_estimator_fast()}) is routed here too; its environment carries \code{fs_fast_core = "logit_or"}, on which the OR vector path \code{.fit_or_for_subgroup_vector()} is taken instead of the OLS core.` In the body, insert after `t <- tt[idx]` (L981) and before `if (!isTRUE(fn_env$adverse_outcome)) y <- -y` (L982):

```r
  # Binary OR core (0.4.0.9000): the tagged fast OR closure marks its
  # environment; the OLS closure carries no marker and takes the path below.
  if (identical(fn_env$fs_fast_core, "logit_or")) {
    return(.fit_or_for_subgroup_vector(y, t, fn_env))
  }

```

and add, immediately after `fit_glm_for_subgroup_fast()` closes (L1003) and before `#' Create Result Row` (L1005):

```r
#' Vector-path binary OR fit for a candidate subgroup
#'
#' The vector twin of `fit_glm_for_subgroup()` with the unadjusted OR
#' closure: the same sign flip, the same existence check
#' (`.fs_existence_reason()` on the 2x2 cells), the same non-convergence
#' rule and the same `NULL` / `nonestimable` contract, with the fit done by
#' `.fs_logit_two_group()` on the already-sliced vectors instead of
#' `stats::glm()` on a data-frame slice.  Called by
#' `fit_glm_for_subgroup_fast()` after its `< 6 rows -> NULL` check.
#'
#' @param y,t Numeric vectors: outcome and treatment, already sliced.
#' @param fn_env The tagged OR closure's environment (`adverse_outcome`,
#'   `treat.name`).
#' @keywords internal
#' @noRd
.fit_or_for_subgroup_vector <- function(y, t, fn_env) {
  if (!isTRUE(fn_env$adverse_outcome)) y <- 1L - y

  why <- .fs_existence_reason("OR", .fs_binary_cells(y, t))
  if (!is.null(why)) return(list(nonestimable = TRUE, reason = why))

  core <- tryCatch(.fs_logit_two_group(y, t, fn_env$treat.name),
                   error = function(e) NULL)
  # The slow chain's error branch: NA estimate with no reason -> NULL.
  if (is.null(core)) return(NULL)
  if (isFALSE(core$converged) && !is.na(core$estimate)) {
    return(list(nonestimable = TRUE, reason = "non-convergent fit"))
  }
  if (is.na(core$estimate)) return(NULL)

  # Wald CI: estimate +/- 1.96 * SE (as in fit_glm_for_subgroup)
  z_crit <- 1.96
  list(
    hr    = core$estimate,
    lower = core$estimate - z_crit * core$se,
    upper = core$estimate + z_crit * core$se,
    med0  = NA_real_,
    med1  = NA_real_
  )
}

```

**`R/forestsearch_main.R`** — the routing predicate. Anchor by content: the comment beginning `# Unadjusted fast path (0.3.3): for the continuous/MD case with no` and the block

```r
    if (outcome_type == "continuous" &&
        identical(ps_adjust_method[[1L]], "none") &&
        length(.fs_adjust_terms(adjust_covariates)) == 0L) {
      estimator_fn <- .make_lm_estimator_fast(
        treat.name      = treat.name,
        outcome.name    = outcome.name,
        adverse_outcome = adverse_outcome
      )
      attr(estimator_fn, "fs_fast_unadjusted") <- TRUE
    }
```

(0.4.0 L2095–L2104; the deprecation task shortened the roxygen above and inserted its guard below the formals, so the block now sits a few lines from there — anchor by content). Replace that block with:

```r
    if (identical(ps_adjust_method[[1L]], "none") &&
        length(.fs_adjust_terms(adjust_covariates)) == 0L) {
      if (outcome_type == "continuous") {
        estimator_fn <- .make_lm_estimator_fast(
          treat.name      = treat.name,
          outcome.name    = outcome.name,
          adverse_outcome = adverse_outcome
        )
        attr(estimator_fn, "fs_fast_unadjusted") <- TRUE
      } else if (outcome_type == "binary" && identical(effect_measure, "OR")) {
        # Binary OR twin (0.4.0.9000): same predicate, same tag, bit-identical
        # glm.fit closure (.make_glm_or_estimator_fast).  RD and RR keep the
        # slow closure (their fallback tiers are untouched).
        estimator_fn <- .make_glm_or_estimator_fast(
          treat.name      = treat.name,
          outcome.name    = outcome.name,
          adverse_outcome = adverse_outcome
        )
        attr(estimator_fn, "fs_fast_unadjusted") <- TRUE
      }
    }
```

and in the comment above it change `for the continuous/MD case` to `for the continuous/MD and binary/OR cases`. `effect_measure` is already resolved to `"OR"` for binary by this point (the `switch(outcome_type, binary = "OR", ...)` default block in Section 1A). Nothing else in the file changes; the `ps_adjust_resolved` rebuild below it (which replaces the closure with an untagged slow one when PS adjustment resolves active) is untouched and now covers both tagged closures.

Then `devtools::document()` → no `man/`/`NAMESPACE` change, restore `DESCRIPTION` if re-stamped; `devtools::install(upgrade = FALSE)`; record `Built`.

**Part-F identity pass** — `dev/reports/lean_fits_2026-10-26/identF.R` (sources `vi_common.R`), binary fixture, `"ii"` args:

- **per candidate:** `base::trace()` on `forestsearch:::fit_glm_for_subgroup_fast` records `id.x` per call and, once, `yy`, `tt` and `estimator_fn` (expected 2577 calls — the Stage-0b count; print it). Reference: `df_ref <- data.frame(yy, tt)` named `c(outcome.name, treat.name)` and the **slow** closure `slow_fn <- forestsearch:::.make_glm_binary_estimator(treat.name, outcome.name, effect_measure = "OR", adverse_outcome = <the fixture's adverse_outcome>)`; for every captured `id.x`, `identical(forestsearch:::fit_glm_for_subgroup_fast(yy, tt, id.x, estimator_fn), forestsearch:::fit_glm_for_subgroup(df_ref, id.x, slow_fn))`, warnings collected and compared, and the kinds (`NULL` / `nonestimable` / fit) tabulated. Also `identical(estimator_fn(df_ref[id.x == 1, , drop = FALSE]), slow_fn(df_ref[id.x == 1, , drop = FALSE]))` for every 50th candidate and every non-fit candidate (closure-level identity of the full returned list). **GATE:** all identical, warning sets equal. Confirm the routing: `attr(estimator_fn, "fs_fast_unadjusted")` is `TRUE` and `environment(estimator_fn)$fs_fast_core` is `"logit_or"`, and `trace()` on `forestsearch:::fit_glm_for_subgroup` (the slow entry) counted **zero** calls during the run.
- **end to end:** binary `"ii"` saved as `identF_binary_ii.rds`, reloaded, compared component by component (including `warnings`) to the committed `stageC_binary_ii.rds`. **GATE.**
- **regression on the other two:** continuous and survival `"ii"` end to end, saved and compared the same way to their committed records. **GATE.**

Commit `R/glm_effect_estimators.R`, `R/subgroup_search.R`, `R/forestsearch_main.R`, `identF.R` and its outputs by explicit path ("Lean binary OR search fit: glm.fit closure and vector path"). Hard cap 15 min.

## 5. Tests — new `tests/testthat/test-lean-fits.R` (gating), then NEWS

Fixtures from the existing helpers (`.make_survival_data()`, `.make_binary_data()`); every `forestsearch()` call well under 30 s; `Sys.setenv(NOT_CRAN = "true")` when running under `Rscript`. Tests:

1. **A, function level.** On `survival::gbsg` (Y = `rfstime / 30.4375`, E = `status`, Treat = `hormon`) with `set.seed(1)`: 25 random subsets of size 40–300, plus edge subsets: all rows; a near-tie set — `Y <- rep(seq(5, 50, by = 5), each = 3) + rep(c(0, 1e-9, 2e-9), 10)` with E and Treat from rows 1:30 (ten clusters of three times that `aeqSurv()` collapses to exact ties — chat confirmed the lean fit without `aeqSurv()` gives a different coefficient here, so this fixture is the tripwire for the `timefix` convention; an ordinary subset is likewise the tripwire for `nocenter`, since the centred fit is not bit-identical); a monotone-likelihood set (all treated events earlier than every control event — both paths warn "Loglik converged before variable 1"); n = 8 with at least one event per arm; a single-arm set (Treat constant — both paths return an all-`NA` row with no warning). For each, `expect_identical(forestsearch:::fit_cox_for_subgroup(Y, E, Treat, rep(1L, length(Y))), ref_fit_cox(Y, E, Treat))` with `ref_fit_cox` from §3 copied verbatim into the test, and `expect_identical()` of the warning texts collected around each side.
2. **F, closure and vector level.** `df <- .make_binary_data(N = 300L, OR_harm = 2.5, seed = 7L)`; `slow <- forestsearch:::.make_glm_binary_estimator("treat", "y", effect_measure = "OR", adverse_outcome = TRUE)`, `fast <- forestsearch:::.make_glm_or_estimator_fast("treat", "y", adverse_outcome = TRUE)`, and the same pair with `adverse_outcome = FALSE`. Slices: 25 random subsets of size 20–300, plus: zero-cell slices (each of the four 2x2 cells empty in turn — `y == treat` on 40 rows is one of them; these return the `nonestimable` reason from the existence check, which is the only way a two-group logistic fit fails here, since the MLE of a 2x2 table with four non-empty cells always exists); a slice with one event among many non-events (extreme fitted probabilities, `glm.fit()` may warn — warnings must match); a 5-row slice (`NULL`); a single-arm slice (`nonestimable`). For each, `expect_identical(fast(slice), slow(slice))` (whole list) and, through the search entry points, `expect_identical(forestsearch:::fit_glm_for_subgroup_fast(df$y, df$treat, id.x, fast), forestsearch:::fit_glm_for_subgroup(df, id.x, slow))`, plus warning-text identity on both comparisons. (Chat ran exactly this battery against the installed 0.4.0 on 2026-10-26: 30/30 identical at both `adverse_outcome` settings, closure and vector level, warnings included; the A battery of item 1 likewise 30/30 against the 0.4.0 formula block, including the monotone-likelihood warning.)
3. **Routing.** With the `trace()`-counter pattern of `test-vi-deprecation.R` (counters on `forestsearch:::fit_glm_for_subgroup` and `forestsearch:::fit_glm_for_subgroup_fast`, `untrace()` on exit): a small `forestsearch()` binary run with `effect_measure = "OR"`, `use_lasso = use_grf = FALSE` (transplant the call shape of the fastest existing binary test; `fs.splits = 20, n.min = 30, maxk = 2, quiet = TRUE, parallel_args = list(plan = "sequential")`) makes zero slow calls and > 0 fast calls; the same call with `effect_measure = "RD"` makes > 0 slow calls and zero fast calls; the same OR call with `adjust_covariates = "age"` makes zero fast calls.

Run `test-lean-fits.R` (gating, `timeout 10m`), then `test-fs-family-report.R` and `test-vi-deprecation.R` (non-gating, same cap each). All with `NOT_CRAN=true`; record the `TOTALS` lines.

**`NEWS.md`** — under `# forestsearch (development version)`, a `### Performance` entry: the two search-stage fits now call the fitting engines directly (`survival::coxph.fit()` for the unadjusted survival fit; `stats::glm.fit()` through a tagged closure for the unadjusted binary odds-ratio fit, the twin of the continuous fast path); results are bit-identical (identity tests in `test-lean-fits.R`; verification record `dev/reports/REPORT_lean_fits_2026-10-26.md`); the measured per-replicate bootstrap time on the committed applications (fill in from §6); nothing changes for adjusted, weighted, RD or RR fits.

Commit `tests/testthat/test-lean-fits.R` and `NEWS.md` by path ("Lean-fit identity tests; NEWS").

## 6. Stage C — timing on the A+F build, and the bootstrap payloads

`dev/reports/lean_fits_2026-10-26/stageC.R`, same protocol as §2 on the installed A+F build (the machine otherwise idle; same order of fixtures): warm-up, three timed calls, one `nb = 20` sequential bootstrap per fixture. Save `stageC_<fixture>_boot.rds` and compare, saved against saved, to `stageA_<fixture>_boot.rds` component by component. **GATE:** all three payloads identical. Write `stageC_timing.csv` and `stageC_summary.txt` (one line per gate — §3, §4, §5 outcomes copied in, §6 gates — and the timing table as `FINDING` lines):

| fixture | call median s: before → after (%) | replicate mean s: before → after (%) | boot outer s: before → after (%) | Stage-0b predicted saving per replicate |
|---|---|---|---|---|
| survival | | | | 1.51 s |
| continuous | | | | 0 |
| binary | | | | 3.13 s |

If a realized per-replicate saving is below 30 % on survival or binary, print `FINDING: below the 30% bar` on its line — report it, do not stop. Hard cap 20 min.

## 7. Report — `dev/reports/REPORT_lean_fits_2026-10-26.md`

1. Provenance: branch, starting HEAD, every commit SHA in order, installed build before (§1), after A (§3), after F (§4); the installed package is now the A+F branch build.
2. The edits, file by file, with `git diff --stat 60a6496b HEAD -- R/ tests/ NEWS.md` and the `@@` hunk headers of each `R/` file (every hunk must sit inside `fit_cox_for_subgroup()`, `fit_glm_for_subgroup_fast()`, the new functions, or the routing predicate — say so).
3. Identity: the per-candidate counts and results of §3 and §4 (candidates captured, identical, warning sets, kinds table), the end-to-end component tables against the committed deprecation record, and the §6 bootstrap-payload table.
4. Timing: the §6 table, with all three raw call times per cell, and one paragraph reading it against Stage 0b's predictions and Larry's 30 % bar.
5. Test results, verbatim `TOTALS` lines, gating and non-gating separately.
6. Script deviations (verbatim if any). OPEN ITEMS (at least: `run_subgroup_sims.R` still routes only the continuous case to a fast closure — a one-predicate follow-up if Larry wants binary simulations to share the gain; the `DESCRIPTION` roxygen stamp, still Larry's call). Ten-line verdict.
7. Closing block: the §8 post-condition commands and their output.

## 8. Commit and post-conditions — GATE

Allowed paths: `dev/tasks/cc_task_lean_fits_2026-10-26.md`; `R/subgroup_search.R`; `R/glm_effect_estimators.R`; `R/forestsearch_main.R`; `tests/testthat/test-lean-fits.R`; `NEWS.md`; `dev/reports/REPORT_lean_fits_2026-10-26.md`; `dev/reports/lean_fits_2026-10-26/`. Before every commit `git diff --cached --name-only` must be a subset of that list. After the final commit print in the closing block:

```
git diff --name-only 60a6496b HEAD | grep -v -E '^(dev/|NEWS\.md$|R/(subgroup_search|glm_effect_estimators|forestsearch_main)\.R$|tests/testthat/test-lean-fits\.R$)'
    # must print nothing (in particular: no man/, no NAMESPACE, no DESCRIPTION)
git diff --numstat 60a6496b HEAD -- R/glm_effect_estimators.R
    # deletions must be 0: the slow OR path is untouched, the new functions are pure additions
grep -n "survival::coxph(" R/subgroup_search.R | grep -v -E "^[0-9]+:[[:space:]]*#"
    # must print exactly one line: the adjusted branch's summary(survival::coxph(cox_fmla, ...)) call
grep -n "coxph.fit\|aeqSurv\|coxph.control" R/subgroup_search.R | grep -v -E "^[0-9]+:[[:space:]]*#" | wc -l
    # must print 3 (one call each, all inside fit_cox_for_subgroup())
grep -n "fs_fast_core\|.make_glm_or_estimator_fast\|.fit_or_for_subgroup_vector\|.fs_logit_two_group" R/*.R | grep -v -E "^[^:]+:[0-9]+:[[:space:]]*#"
    # must print only lines inside R/glm_effect_estimators.R (definitions), R/subgroup_search.R (the dispatch and the helper) and R/forestsearch_main.R (the routing call)
grep -c "PASS" dev/reports/lean_fits_2026-10-26/stageC_summary.txt; grep -c "FAIL" dev/reports/lean_fits_2026-10-26/stageC_summary.txt
    # the second count must be 0
git status --porcelain | grep -v '^?? '
    # must print nothing
```

If any command prints something unexpected, stop and report what it prints.

Finally, copy `dev/reports/REPORT_lean_fits_2026-10-26.md`, `dev/reports/lean_fits_2026-10-26/stageC_summary.txt` and `dev/reports/lean_fits_2026-10-26/stageC_timing.csv` to `~/Downloads` and print `ls -la` of the three copies.

## 9. Out of scope

- The consistency-stage items B/E/G/H of Stage 0b (recorded as deferred; each under 15 % of a replicate in the application configurations).
- Adjusted, weighted (propensity), stratified, RD and RR fits: they keep their formula paths unchanged. A lean adjusted/weighted `coxph.fit()`/`glm.fit()` extension is a later proposal, not this task.
- `run_subgroup_sims.R`'s own routing predicate (OPEN ITEM, §7).
- `.g3_cox_lean()` (Guo–He path): untouched; it is a different, non-bit-identical estimator by design.
- No `DESCRIPTION` version bump; no `R CMD check`; no full test suite; no vignette build; no push; no merge.
