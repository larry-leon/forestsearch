# File: R/dina_subgroup_bootstrap.R
# Part of the forestsearch package.
#
# Selection-adjusted bootstrap inference for dina_subgroup().
#
# For each of n_boot bootstrap resamples of the original data frame:
#   1. Refit DINA via dina()
#   2. Re-run dina_subgroup() on that fit
#   3. Record (covariate, direction, threshold, n_subgroup, mean_tau_hat)
#
# Aggregate per-iteration results into a non-parametric bootstrap
# distribution.  The percentile CI on mean_tau_hat is selection-adjusted
# because each bootstrap iteration re-runs the full search procedure,
# allowing different (covariate, direction, threshold) selections to
# contribute to the distribution.  If the procedure is unstable (the
# chosen subgroup varies substantially across bootstrap samples), the
# resulting CI is appropriately wide.


# ---------------------------------------------------------------------------
# Main exported function
# ---------------------------------------------------------------------------

#' Selection-adjusted bootstrap CI for dina_subgroup()
#'
#' Performs selection-adjusted inference for the harm subgroup discovered
#' by `dina_subgroup()`.  For each of `n_boot` bootstrap resamples of the
#' input data frame, the function refits DINA via `dina()` and re-runs
#' the subgroup search via `dina_subgroup()`, then reports three classes
#' of output:
#'
#'   * An **effect CI** on `a*^T beta_b`, where `a* = (1, x_bar_S)` is
#'     held FIXED from the original-data subgroup `S*` and `beta_b` is
#'     each iteration's DINA coefficient vector.  This is the
#'     bootstrap analog of the Wald CI returned by
#'     `dina_subgroup()` -- conditional on the discovered subgroup
#'     definition, with width reflecting sampling variance of
#'     `beta-hat`.
#'   * **Stability CIs** on the bootstrap-selected `n_subgroup` and
#'     `threshold`, restricted to iterations that selected the same
#'     `(covariate, direction)` as the original-data subgroup.
#'   * A **selection-frequency table** counting how often each
#'     `(covariate, direction)` was chosen.  The covariate-selection
#'     stability diagnostic.
#'
#' Each bootstrap iteration calls `dina()` (single-pass, sandwich
#' variance) under the hood, not `dina_bagged()`.
#'
#' @section Why the bootstrap target is the *fixed-subgroup* effect:
#' A naive bootstrap of `dina_subgroup()`'s `mean_tau_hat` is
#' uninformative.  The search rule selects the *largest* subgroup
#' whose mean tau-hat exceeds `m_diff`, so the chosen subgroup lies
#' at the boundary by construction: adding one more patient would
#' push the mean below `m_diff`.  Each bootstrap iteration's
#' selected subgroup is similarly pinned, so the bootstrap
#' distribution of `mean_tau_hat` clusters tightly above `m_diff`
#' and the resulting CI quantifies threshold-grid granularity
#' rather than effect uncertainty.  The fixed-subgroup linear
#' functional `a*^T beta_b` instead targets the clinically
#' meaningful estimand "treatment effect for patients whose
#' covariates lie in the discovered subgroup".  Because DINA's
#' `beta-hat` is cross-fit (and hence asymptotically unbiased even
#' for data-driven contrasts), the percentile distribution of
#' `a*^T beta_b` gives a non-parametric CI that complements the
#' sandwich Wald CI of `dina_subgroup()`.
#'
#' @section Performance:
#' With `parallel = "boots"`, each bootstrap iteration runs in a worker
#' via `foreach::foreach() %dofuture%`.  Set `future::plan()` to a
#' parallel backend (e.g., `future::multisession`) before calling, and
#' reset to `"sequential"` afterwards.  Workers run from the installed
#' package; `devtools::install()` is required after code edits before
#' the parallel branch will see them.
#'
#' @param df data frame containing the outcome, treatment, and
#'   covariate columns.
#' @param outcome character(1); name of the response column in `df`.
#'   For `family = "cox"` this is the event/censoring time column.
#' @param treatment character(1); name of the binary treatment
#'   indicator column.
#' @param covariates character vector; names of the covariate columns
#'   to search over.  All must be numeric.
#' @param family one of `"gaussian"`, `"binomial"`, `"poisson"`, `"cox"`.
#' @param status character(1) or `NULL`; for `family = "cox"` only,
#'   the event indicator column name.
#' @param m_diff scalar harm threshold on the natural-parameter scale.
#' @param n_min positive integer; minimum subgroup size.  Default `60L`.
#' @param direction one of `"both"` (default), `"left"`, `"right"`.
#' @param alpha confidence level for the percentile CI.  Default `0.05`
#'   (percentile endpoints at the `alpha/2` and `1 - alpha/2` quantiles).
#' @param n_boot positive integer; number of bootstrap iterations.
#'   Default `200L`.  Warns if `n_boot < 100`.
#' @param parallel one of `"none"` (default) or `"boots"`.  When
#'   `"boots"`, bootstrap iterations dispatch via
#'   `foreach::foreach() %dofuture%` against the active
#'   `future::plan()`.
#' @param seed optional integer seed for reproducibility of the
#'   bootstrap indices.  Default `NULL`.
#' @param verbose logical; if `TRUE`, print a progress message every
#'   10 iterations (serial mode) or a start summary (parallel mode).
#'   Default `FALSE`.
#' @param ... further arguments forwarded to `dina()` for each
#'   per-iteration DINA fit (e.g., `n_folds`, `cens_type`,
#'   `cens_params`, `propensity_method`, `baseline_method`).
#'
#' @return An object of class `"dina_subgroup_bootstrap"`, a list with
#'   components:
#'   \describe{
#'     \item{point}{the `"dina_subgroup"` object computed on the
#'       original (unbootstrapped) data.}
#'     \item{effect_ci}{named length-2 numeric vector (`lower`,
#'       `upper`) giving the bootstrap percentile CI on the
#'       fixed-subgroup linear functional `a*^T beta_b`, where
#'       `a* = (1, x_bar_S)` is the covariate-mean vector of the
#'       original-data subgroup and `beta_b` is each bootstrap
#'       iteration's DINA coefficient vector.  This is the
#'       meaningful effect-uncertainty CI; width reflects sampling
#'       variance of `beta-hat`.}
#'     \item{effect_dist}{numeric vector of length `n_boot`
#'       containing `a*^T beta_b` for each iteration (`NA` for
#'       iterations whose DINA fit failed).}
#'     \item{n_subgroup_ci, threshold_ci}{named length-2 numeric
#'       vectors giving percentile CIs for the structural quantities
#'       of bootstrap-selected subgroups, restricted to iterations
#'       that selected the same `(covariate, direction)` as the
#'       original-data subgroup.  Subgroup-size and boundary-location
#'       stability diagnostics.}
#'     \item{n_modal_match}{integer; number of bootstrap iterations
#'       whose selection matched the original-data
#'       `(covariate, direction)`.  Denominator for the stability CIs.}
#'     \item{selection_frequency}{a `table` of
#'       `(covariate, direction)` selection counts across found
#'       iterations.  Covariate-selection stability diagnostic.}
#'     \item{boot_results}{data frame with one row per bootstrap
#'       iteration; columns `found`, `covariate`, `direction`,
#'       `threshold`, `n_subgroup`, `mean_tau_hat`, `failed`.  The
#'       `mean_tau_hat` column records the bootstrap-selected
#'       subgroup's boundary value (pinned near `m_diff` by the
#'       search rule) and is retained for diagnostic transparency
#'       rather than as a meaningful effect estimate.}
#'     \item{n_boot, n_boot_found, n_boot_failed}{convergence
#'       diagnostics.}
#'     \item{alpha, m_diff, n_min, family}{echoed inputs.}
#'     \item{call}{the matched call.}
#'   }
#'
#' @seealso `dina_subgroup()` for the underlying threshold search;
#'   `dina()` for the per-iteration DINA fit.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' n <- 400
#' df_demo <- data.frame(
#'   w  = stats::rbinom(n, 1, 0.5),
#'   x1 = stats::runif(n, -1, 1),
#'   x2 = stats::runif(n, -1, 1)
#' )
#' tau_x      <- 0.4 + 1.2 * df_demo$x1
#' df_demo$y  <- 0.5 * df_demo$x1 + df_demo$w * tau_x + stats::rnorm(n)
#'
#' bs <- dina_subgroup_bootstrap(
#'   df         = df_demo,
#'   outcome    = "y",
#'   treatment  = "w",
#'   covariates = c("x1", "x2"),
#'   family     = "gaussian",
#'   m_diff     = 0.5,
#'   n_boot     = 50L,
#'   seed       = 1L
#' )
#' print(bs)
#' }
#'
#' @importFrom foreach foreach
#' @importFrom future nbrOfWorkers
#' @export
dina_subgroup_bootstrap <- function(df,
                                    outcome,
                                    treatment,
                                    covariates,
                                    family = c("gaussian", "binomial",
                                               "poisson", "cox"),
                                    status = NULL,
                                    m_diff,
                                    n_min = 60L,
                                    direction = c("both", "left", "right"),
                                    alpha = 0.05,
                                    n_boot = 200L,
                                    parallel = c("none", "boots"),
                                    seed = NULL,
                                    verbose = FALSE,
                                    ...) {

  family    <- match.arg(family)
  parallel  <- match.arg(parallel)
  direction <- match.arg(direction)
  call      <- match.call()

  # ---- Argument validation ----------------------------------------------
  if (missing(m_diff) || length(m_diff) != 1L || !is.numeric(m_diff) ||
      !is.finite(m_diff)) {
    stop("`m_diff` must be a single finite numeric value.")
  }
  n_boot <- as.integer(n_boot)
  if (length(n_boot) != 1L || is.na(n_boot) || n_boot < 2L) {
    stop("`n_boot` must be a single integer >= 2.")
  }
  if (n_boot < 100L) {
    warning("`n_boot` < 100; bootstrap percentile CI will be unstable. ",
            "Consider increasing `n_boot`.", call. = FALSE)
  }
  if (length(alpha) != 1L || !is.numeric(alpha) ||
      alpha <= 0 || alpha >= 1) {
    stop("`alpha` must be a single numeric in (0, 1).")
  }

  if (parallel == "boots") {
    n_workers <- future::nbrOfWorkers()
    if (n_workers <= 1L) {
      warning(
        "`parallel = \"boots\"` requested but `future::plan()` reports ",
        "only 1 worker; bootstraps will run serially. ",
        "Call e.g. `future::plan(future::multisession, workers = N)` ",
        "before this function to enable actual parallelism.",
        call. = FALSE
      )
    } else if (verbose) {
      message(sprintf(
        "dina_subgroup_bootstrap: running %d bootstraps across %d workers",
        n_boot, n_workers
      ))
    }
  }

  dina_args <- list(...)
  # Prevent the user's `seed` (if passed via ...) from being applied to
  # every per-iteration DINA fit identically.  The bootstrap-index seed
  # is governed by the explicit `seed` argument above.
  dina_args[["seed"]] <- NULL

  # ---- Point estimate on the original data ------------------------------
  fit_point <- do.call(dina, c(
    list(df = df, outcome = outcome, treatment = treatment,
         covariates = covariates, family = family, status = status),
    dina_args
  ))
  sg_point <- dina_subgroup(
    fit = fit_point, df = df, covariates = covariates,
    m_diff = m_diff, n_min = n_min,
    direction = direction, alpha = alpha
  )

  # ---- Bootstrap indices ------------------------------------------------
  n <- nrow(df)
  boot <- .dina_make_bootstrap_indices(n, n_boot, seed = seed)

  # ---- Bootstrap loop ---------------------------------------------------
  if (parallel == "boots") {
    b <- NULL  # silence R CMD check "no visible binding"
    boot_list <- foreach::foreach(
      b = seq_len(n_boot),
      .options.future = list(
        packages = "forestsearch",
        seed     = TRUE
      )
    ) %dofuture% {
      .dina_safe_one_bootstrap(
        idx = boot$indices[, b],
        df = df, outcome = outcome, treatment = treatment,
        covariates = covariates, family = family, status = status,
        m_diff = m_diff, n_min = n_min, direction = direction,
        alpha = alpha, dina_args = dina_args
      )
    }
  } else {
    boot_list <- lapply(seq_len(n_boot), function(b) {
      if (verbose && b %% 10L == 0L) {
        message(sprintf("dina_subgroup_bootstrap: iteration %d of %d",
                        b, n_boot))
      }
      .dina_safe_one_bootstrap(
        idx = boot$indices[, b],
        df = df, outcome = outcome, treatment = treatment,
        covariates = covariates, family = family, status = status,
        m_diff = m_diff, n_min = n_min, direction = direction,
        alpha = alpha, dina_args = dina_args
      )
    })
  }

  # ---- Aggregate --------------------------------------------------------
  # Per-iteration bootstrap-selected subgroup statistics (for stability
  # diagnostics).  The effect-uncertainty target -- the fixed-subgroup
  # linear functional a*^T beta_b -- is computed separately below using
  # each iteration's full coefficient vector and the FIXED a* from the
  # original-data subgroup.
  boot_df <- data.frame(
    found        = vapply(boot_list, function(r) {
                    if (is.na(r$found)) NA else as.logical(r$found)
                   }, logical(1L)),
    covariate    = vapply(boot_list, `[[`, character(1L), "covariate"),
    direction    = vapply(boot_list, `[[`, character(1L), "direction"),
    threshold    = vapply(boot_list, `[[`, numeric(1L),   "threshold"),
    n_subgroup   = vapply(boot_list, `[[`, integer(1L),   "n_subgroup"),
    mean_tau_hat = vapply(boot_list, `[[`, numeric(1L),   "mean_tau_hat"),
    failed       = vapply(boot_list, `[[`, logical(1L),   "failed"),
    stringsAsFactors = FALSE
  )

  n_boot_failed <- sum(boot_df$failed)
  fit_mask      <- !boot_df$failed
  found_mask    <- !is.na(boot_df$found) & boot_df$found
  n_boot_found  <- sum(found_mask)

  # ---- Fixed-subgroup effect distribution -------------------------------
  # The point-estimate subgroup S* and its covariate-mean vector a* are
  # held FIXED from the original-data fit.  For each bootstrap iteration
  # whose DINA fit succeeded, compute T_b = a*^T beta_b.  Per the design
  # rationale (see dina_inference_design_notes.md, Section 4), this is the
  # meaningful effect-uncertainty target -- width reflects sampling
  # variance of beta-hat, not boundary slack.
  effect_dist <- rep(NA_real_, n_boot)
  effect_ci   <- c(lower = NA_real_, upper = NA_real_)

  if (isTRUE(sg_point$found)) {
    X_all <- .dina_extract_X_from_df(df, covariates)
    x_bar_S <- colMeans(X_all[sg_point$mask, , drop = FALSE])
    a_star  <- c(1, x_bar_S)

    # Stack per-iteration beta vectors into a B x (d+1) matrix; rows for
    # failed iterations remain NA.  Matrix-vector product yields a_star^T
    # beta_b for every iteration in one shot.
    beta_mat <- do.call(rbind, lapply(boot_list, `[[`, "beta"))
    effect_dist <- as.numeric(beta_mat %*% a_star)

    effect_finite <- effect_dist[is.finite(effect_dist)]
    if (length(effect_finite) >= 2L) {
      qs <- stats::quantile(
        effect_finite, probs = c(alpha / 2, 1 - alpha / 2),
        names = FALSE, na.rm = TRUE
      )
      effect_ci <- c(lower = qs[1L], upper = qs[2L])
    }
  }

  # ---- Stability CIs on subgroup structure ------------------------------
  # n_subgroup and threshold are restricted to iterations that selected
  # the SAME (covariate, direction) as the original-data subgroup, since
  # thresholds on different covariates are not comparable.
  n_subgroup_ci <- c(lower = NA_real_, upper = NA_real_)
  threshold_ci  <- c(lower = NA_real_, upper = NA_real_)
  n_modal_match <- NA_integer_

  if (isTRUE(sg_point$found) && n_boot_found > 0L) {
    modal_match <- found_mask &
      boot_df$covariate == sg_point$covariate &
      boot_df$direction == sg_point$direction
    n_modal_match <- as.integer(sum(modal_match))

    if (n_modal_match >= 2L) {
      n_sg_qs <- stats::quantile(
        boot_df$n_subgroup[modal_match],
        probs = c(alpha / 2, 1 - alpha / 2),
        names = FALSE, na.rm = TRUE
      )
      n_subgroup_ci <- c(lower = n_sg_qs[1L], upper = n_sg_qs[2L])

      thr_qs <- stats::quantile(
        boot_df$threshold[modal_match],
        probs = c(alpha / 2, 1 - alpha / 2),
        names = FALSE, na.rm = TRUE
      )
      threshold_ci <- c(lower = thr_qs[1L], upper = thr_qs[2L])
    }
  }

  # ---- Selection-frequency diagnostic -----------------------------------
  selection_frequency <- if (n_boot_found > 0L) {
    table(
      covariate = boot_df$covariate[found_mask],
      direction = boot_df$direction[found_mask]
    )
  } else {
    NULL
  }

  if (n_boot_failed > 0L) {
    warning(sprintf("%d of %d bootstrap iterations errored and were excluded.",
                    n_boot_failed, n_boot), call. = FALSE)
  }
  if (sum(fit_mask) < n_boot %/% 2L) {
    warning(sprintf(
      "Only %d of %d bootstrap iterations produced a usable DINA fit; ",
      sum(fit_mask), n_boot
    ), "effect CI is based on a small sample.", call. = FALSE)
  }

  out <- list(
    point                = sg_point,
    effect_ci            = effect_ci,
    effect_dist          = effect_dist,
    n_subgroup_ci        = n_subgroup_ci,
    threshold_ci         = threshold_ci,
    n_modal_match        = n_modal_match,
    selection_frequency  = selection_frequency,
    boot_results         = boot_df,
    n_boot               = n_boot,
    n_boot_found         = n_boot_found,
    n_boot_failed        = n_boot_failed,
    alpha                = alpha,
    m_diff               = m_diff,
    n_min                = n_min,
    family               = family,
    call                 = call
  )
  class(out) <- "dina_subgroup_bootstrap"
  out
}


# ---------------------------------------------------------------------------
# Internal helper: one bootstrap iteration with error handling
# ---------------------------------------------------------------------------

#' Run one bootstrap iteration with error handling, returning a uniform
#' list result.
#'
#' Wraps the dina() + dina_subgroup() pipeline in tryCatch so that
#' per-iteration failures (e.g., from a pathological bootstrap sample
#' that triggers a fit error) surface as
#' `list(found = NA, ..., failed = TRUE, message = ...)` instead of
#' propagating as conditions.  Both parallel and serial code paths
#' consume the same list contract downstream.
#'
#' @noRd
.dina_safe_one_bootstrap <- function(idx, df, outcome, treatment, covariates,
                                     family, status,
                                     m_diff, n_min, direction, alpha,
                                     dina_args) {
  d <- length(covariates)
  na_beta <- rep(NA_real_, d + 1L)

  tryCatch({
    df_b <- df[idx, , drop = FALSE]

    fit_b <- do.call(dina, c(
      list(df = df_b, outcome = outcome, treatment = treatment,
           covariates = covariates, family = family, status = status),
      dina_args
    ))

    # Full coefficient vector for downstream fixed-subgroup linear-functional
    # computation (a*^T beta_b with a* fixed from the original data).
    beta_b <- as.numeric(stats::coef(fit_b))

    sg_b <- dina_subgroup(
      fit = fit_b, df = df_b, covariates = covariates,
      m_diff = m_diff, n_min = n_min,
      direction = direction, alpha = alpha
    )

    if (isTRUE(sg_b$found)) {
      list(
        beta         = beta_b,
        found        = TRUE,
        covariate    = sg_b$covariate,
        direction    = sg_b$direction,
        threshold    = as.numeric(sg_b$threshold),
        n_subgroup   = as.integer(sg_b$n_subgroup),
        mean_tau_hat = as.numeric(sg_b$mean_tau_hat),
        failed       = FALSE,
        message      = ""
      )
    } else {
      list(
        beta         = beta_b,
        found        = FALSE,
        covariate    = NA_character_,
        direction    = NA_character_,
        threshold    = NA_real_,
        n_subgroup   = NA_integer_,
        mean_tau_hat = NA_real_,
        failed       = FALSE,
        message      = ""
      )
    }
  }, error = function(e) {
    list(
      beta         = na_beta,
      found        = NA,
      covariate    = NA_character_,
      direction    = NA_character_,
      threshold    = NA_real_,
      n_subgroup   = NA_integer_,
      mean_tau_hat = NA_real_,
      failed       = TRUE,
      message      = conditionMessage(e)
    )
  })
}


# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' Print a `dina_subgroup_bootstrap` result.
#'
#' @param x a `"dina_subgroup_bootstrap"` object.
#' @param digits number of digits for numeric summary.
#' @param max_select maximum number of rows of the selection-frequency
#'   table to display.  Default `5L`.
#' @param ... unused.
#' @return invisibly returns `x`.
#' @export
print.dina_subgroup_bootstrap <- function(x,
                                          digits = max(3L,
                                                       getOption("digits") - 3L),
                                          max_select = 5L,
                                          ...) {
  cat("Bootstrap inference for dina_subgroup()\n")
  cat("  Family:           ", x$family, "\n", sep = "")
  cat("  m_diff:           ", format(x$m_diff, digits = digits), "\n", sep = "")
  cat("  n_min:            ", x$n_min, "\n", sep = "")
  cat("  Bootstrap iters:  ", x$n_boot_found, " found / ",
      sum(!is.na(x$boot_results$found) & !x$boot_results$found),
      " not-found / ", x$n_boot_failed, " failed (of ",
      x$n_boot, " total)\n", sep = "")

  cat("\nPoint estimate (original data):\n")
  if (isTRUE(x$point$found)) {
    cmp <- if (x$point$direction == "left") "<=" else ">="
    cat("  Signature:    ", x$point$covariate, " ", cmp, " ",
        format(x$point$threshold, digits = digits), "\n", sep = "")
    cat("  n_subgroup:   ", x$point$n_subgroup, " / ",
        x$point$n_total, "\n", sep = "")
    cat("  mean tau-hat: ", format(x$point$mean_tau_hat, digits = digits),
        "\n", sep = "")
  } else {
    cat("  No qualifying subgroup on the original data.\n")
  }

  ci_pct <- format(100 * (1 - x$alpha), digits = 3L)

  cat("\nFixed-subgroup effect CI (", ci_pct,
      " pct percentile):  [",
      format(x$effect_ci[["lower"]], digits = digits), ", ",
      format(x$effect_ci[["upper"]], digits = digits), "]\n", sep = "")
  cat("  Target:       a*^T beta_b with a* = (1, x_bar_S) fixed from\n")
  cat("                the original-data subgroup S*.\n")
  cat("  Interpretation: bootstrap percentile CI for the effect in\n")
  cat("                the discovered subgroup; complements the Wald\n")
  cat("                CI returned by dina_subgroup().\n")

  if (!is.na(x$n_modal_match) && x$n_modal_match >= 2L) {
    cat("\nStructural stability CIs (", ci_pct,
        " pct percentile, ", x$n_modal_match,
        " iterations matching the modal (covariate, direction)):\n", sep = "")
    cat("  n_subgroup:   [",
        format(x$n_subgroup_ci[["lower"]], digits = digits), ", ",
        format(x$n_subgroup_ci[["upper"]], digits = digits), "]\n", sep = "")
    cat("  threshold:    [",
        format(x$threshold_ci[["lower"]], digits = digits), ", ",
        format(x$threshold_ci[["upper"]], digits = digits), "]\n", sep = "")
  }

  if (!is.null(x$selection_frequency)) {
    sel_long <- as.data.frame(x$selection_frequency,
                              stringsAsFactors = FALSE)
    sel_long <- sel_long[sel_long$Freq > 0L, , drop = FALSE]
    sel_long <- sel_long[order(-sel_long$Freq), , drop = FALSE]
    n_show <- min(max_select, nrow(sel_long))
    cat("\nTop ", n_show, " (covariate, direction) selections across found iterations:\n",
        sep = "")
    for (i in seq_len(n_show)) {
      cat(sprintf("  %s (%s):  %d of %d  (%.0f pct)\n",
                  sel_long$covariate[i], sel_long$direction[i],
                  sel_long$Freq[i], x$n_boot_found,
                  100 * sel_long$Freq[i] / x$n_boot_found))
    }
    if (nrow(sel_long) > n_show) {
      cat(sprintf("  ... %d more (covariate, direction) combinations.\n",
                  nrow(sel_long) - n_show))
    }
  }

  invisible(x)
}
