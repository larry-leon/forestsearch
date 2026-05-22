# Silence R CMD check NOTE for non-standard evaluation in data.table.
# The bare `m` symbol below is a data.table column name, not a global.
utils::globalVariables(c("m"))

#' Compute Naive and Split-Derived 95% CIs for Pareto Frontier Members
#'
#' For each candidate subgroup on \code{fs$grp.consistency$out_sg$pareto_frontier},
#' computes three 95\% confidence intervals on the effect estimate:
#'
#' \itemize{
#'   \item \strong{Naive CI}: full-sample Wald CI from a Cox or GLM model
#'     refit on the subgroup's data.  This is the textbook CI that ignores
#'     the multiple-testing / subgroup-search context.  Anti-conservative
#'     by construction.
#'   \item \strong{Split CI}: a 95\% interval centered at the full-sample
#'     estimate, with SE estimated by the empirical standard deviation of
#'     the \eqn{2S} individual half-sample effect estimates produced by
#'     \code{n_splits} random 50/50 splits.  This is a half-jackknife SE
#'     (\cite{Shao1996}).  Naming: the \eqn{\sim} in column labels marks
#'     this as a resampling-derived approximation, not a model-based CI.
#'   \item \strong{FSBC-mimic CI}: a bias-corrected interval inspired by
#'     the bootstrap bias-correction algorithm of \cite{Leon2024fs}
#'     (eq 7), with the selection-on-bootstrap-data term
#'     \eqn{\eta_b^*(\hat H_b^*)} zeroed out (the selected subgroup is
#'     treated as fixed across half-jackknife replicates).  The
#'     bias-corrected estimate is
#'     \eqn{\hat\beta^{\mathrm{FSBC}} = 2\hat\beta - \overline{\hat\beta^{(h)}}}
#'     where \eqn{\overline{\hat\beta^{(h)}}} is the mean of the
#'     \eqn{2S} half-sample estimates; the SE is the same
#'     half-jackknife SE as for the Split CI.  See
#'     \dQuote{Details — FSBC-mimic interpretation} below.
#' }
#'
#' All three CIs are computed post-hoc on the returned \code{forestsearch}
#' object; this function does not modify any internal state.
#'
#' @section Details — Half-jackknife SE:
#' For each of \code{n_splits} random 50/50 splits we fit the effect model
#' separately on each half, obtaining \eqn{2S} per-half estimates
#' \eqn{\{\hat\beta^{(s,h)} : s=1,\ldots,S; h=1,2\}}.  The half-jackknife
#' SE estimator is
#' \deqn{\widehat{\mathrm{SE}} = \mathrm{sd}\bigl\{\hat\beta^{(s,h)}\bigr\},}
#' i.e., the plain empirical SD of the \eqn{2S} estimates.  No
#' \eqn{\sqrt{2}} correction is applied: simulation studies confirm that
#' this estimator's CI achieves nominal coverage in a no-selection
#' setting (see the package's pressure-test script).  Earlier versions
#' of this function used the average \eqn{(\hat\beta^{(s,1)} +
#' \hat\beta^{(s,2)})/2} per split; that approach is invalid because
#' the two halves are complementary and their average is nearly
#' constant.
#'
#' @section Details — FSBC-mimic interpretation:
#' The FSBC-mimic CI is a \emph{quick approximation} to the bootstrap
#' bias-corrected estimator of \cite{Leon2024fs}.  Substituting the
#' \eqn{2S} per-half estimates for bootstrap replicates, and treating
#' the selected subgroup \eqn{\hat H} as fixed (i.e., omitting the
#' \eqn{\eta_b^*(\hat H_b^*) = \hat\beta_b^*(\hat H_b^*) -
#' \hat\beta(\hat H_b^*)} term that captures variability in FS selection
#' on the bootstrap data), the bias-corrected estimator from paper
#' eq (7) collapses to
#' \deqn{\hat\beta^{\mathrm{FSBC}}(\hat H) = \hat\beta(\hat H) -
#' \frac{1}{2S}\sum_{s,h}\bigl[\hat\beta^{(s,h)} - \hat\beta(\hat H)\bigr]
#' = 2\hat\beta(\hat H) - \overline{\hat\beta^{(h)}}.}
#' The SE is the same half-jackknife SE as for the Split CI.
#'
#' The FSBC-mimic CI is \emph{not} a full FSBC CI: it omits the
#' selection-on-bootstrap-data variance.  Use this column as a quick
#' bias-corrected diagnostic, not for inference.  For full FSBC CIs use
#' \code{\link{forestsearch_bootstrap_dofuture}}.
#'
#' @param fs A \code{forestsearch} result object.
#' @param n_splits Integer.  Number of 50/50 random splits per frontier
#'   member for the Split CI.  Default \code{1000L}.
#' @param ci_level Numeric in \code{(0, 1)}.  Default \code{0.95}.
#' @param seed Integer or \code{NULL}.  Seed for reproducibility of the
#'   per-split halving.  Default \code{NULL} (no explicit seed).
#' @param verbose Logical.  Default \code{FALSE}.  When \code{TRUE},
#'   prints per-frontier-member progress.
#'
#' @return A \code{data.table} keyed by frontier row index (\code{m}) with
#'   columns:
#'   \describe{
#'     \item{\code{m}}{Candidate id (matches \code{pareto_frontier$m}).}
#'     \item{\code{estimate}}{Refit effect estimate on the natural scale.}
#'     \item{\code{naive_lcl}, \code{naive_ucl}}{Naive Wald CI bounds.}
#'     \item{\code{split_sd}}{Half-jackknife SE estimate: empirical SD of
#'       the \eqn{2S} per-half effect estimates on the SE scale (log for
#'       ratio measures, natural for differences).}
#'     \item{\code{split_lcl}, \code{split_ucl}}{Split-derived CI bounds.}
#'     \item{\code{fsbc_estimate}}{Bias-corrected effect estimate on the
#'       natural scale (\eqn{2\hat\beta - \overline{\hat\beta^{(h)}}} on
#'       the SE scale, then back-transformed).}
#'     \item{\code{fsbc_lcl}, \code{fsbc_ucl}}{FSBC-mimic CI bounds.}
#'     \item{\code{n_valid_splits}}{Number of splits that produced valid
#'       per-half estimates.}
#'   }
#'   Rows where the refit fails (e.g., zero events in one arm) have
#'   \code{NA} in all CI columns.
#'
#' @seealso \code{\link{pareto_frontier_table}},
#'   \code{\link{frontier_member_flags}}, \code{\link{plot_pareto_frontier}}.
#'
#' @importFrom data.table data.table setkey is.data.table
#' @importFrom stats sd quantile coef vcov
#' @export
compute_frontier_cis <- function(fs,
                                 n_splits = 1000L,
                                 ci_level = 0.95,
                                 seed     = NULL,
                                 verbose  = FALSE) {

  # --- 0. Resolve frontier and data ---------------------------------------
  out_sg   <- tryCatch(fs$grp.consistency$out_sg, error = function(e) NULL)
  frontier <- tryCatch(out_sg$pareto_frontier,    error = function(e) NULL)
  if (is.null(frontier) || !data.table::is.data.table(frontier) ||
      nrow(frontier) == 0L) {
    return(data.table::data.table())
  }

  df.est <- tryCatch(fs$df.est, error = function(e) NULL)
  if (is.null(df.est) || !nrow(df.est)) {
    stop("compute_frontier_cis(): fs$df.est is missing or empty.",
         call. = FALSE)
  }

  args <- fs$args_call_all
  if (is.null(args)) {
    stop("compute_frontier_cis(): fs$args_call_all is missing.",
         call. = FALSE)
  }

  treat.name   <- args$treat.name
  outcome.name <- args$outcome.name
  event.name   <- args$event.name
  offset.name  <- args$offset.name
  outcome_type <- fs$outcome_type %||% "survival"
  effect_measure <- fs$effect_measure %||%
    switch(outcome_type, survival = "HR", binary = "RD",
           continuous = "MD", count = "IRR")

  # --- 1. Build the per-data-slice estimator closure ----------------------
  est_fn <- tryCatch(
    make_effect_estimator(
      outcome_type   = outcome_type,
      treat.name     = treat.name,
      outcome.name   = outcome.name,
      event.name     = event.name,
      offset.name    = offset.name,
      effect_measure = effect_measure,
      adverse_outcome   = args$adverse_outcome %||% TRUE,
      adjust_covariates = NULL,        # frontier CIs are unadjusted by design
      ps_adjust_method  = "none"
    ),
    error = function(e) {
      stop("compute_frontier_cis(): could not construct estimator closure: ",
           conditionMessage(e), call. = FALSE)
    }
  )

  is_log_scale <- effect_measure %in% c("HR", "OR", "RR", "IRR")
  z_crit       <- stats::qnorm(0.5 + ci_level / 2)

  # --- 2. Iterate over frontier members -----------------------------------
  if (!is.null(seed)) set.seed(seed)

  m_cols <- grep("^M\\.", names(frontier), value = TRUE)
  out    <- vector("list", nrow(frontier))

  for (k in seq_len(nrow(frontier))) {
    row_k <- frontier[k, ]

    # Cuts for this frontier member: M.* values, excluding empty/NA.
    # Each cut is a human-readable label like "{er <= 0}" -- NOT a column
    # name on df.est.  Strip the outer braces and evaluate the resulting
    # expression ("er <= 0") against df.est to recover per-subject
    # membership.  This handles simple inequalities ("{er <= 0}"),
    # compound cuts ("{age <= 47} & !{age <= 30}"), and any other syntax
    # the package emits in human-label form.
    cuts <- unlist(row_k[, m_cols, with = FALSE], use.names = FALSE)
    cuts <- cuts[!is.na(cuts) & nzchar(cuts)]

    keep <- rep(TRUE, nrow(df.est))
    bad_cuts <- character(0)
    for (c_label in cuts) {
      # Translate the human-readable label into an evaluable expression.
      # See .label_to_expr() in forestsearch_helpers.R for the rules; the
      # helper errors loudly on malformed labels rather than silently
      # mangling them, and produces operationally correct complements
      # for both numeric comparisons and binary factor identifiers
      # (the old gsub-based path mishandled the latter, producing an
      # Ops.factor warning plus an all-NA mask that silently zeroed
      # out the subgroup).
      expr_text <- .label_to_expr(c_label)
      mask <- tryCatch(
        eval(parse(text = expr_text), envir = df.est),
        error = function(e) NULL)
      if (is.null(mask) || length(mask) != nrow(df.est)) {
        bad_cuts <- c(bad_cuts, c_label)
        next
      }
      mask <- as.logical(mask)
      keep <- keep & !is.na(mask) & mask
    }
    if (length(bad_cuts) > 0L) {
      warning(sprintf(
        "compute_frontier_cis(): frontier row %d could not evaluate cut(s) on df.est: %s. Returning NA for this row.",
        k, paste(bad_cuts, collapse = "; ")), call. = FALSE)
      out[[k]] <- .empty_ci_row(row_k)
      next
    }

    sub <- df.est[keep, , drop = FALSE]
    n_sub <- nrow(sub)
    if (n_sub < 10L) {
      out[[k]] <- .empty_ci_row(row_k)
      next
    }

    # --- 2a. Naive (full-sample) CI ---------------------------------------
    # For survival outcomes we use the robust (sandwich) SE to match the
    # CI reported by sg_tables() / cox_summary(), which the rest of the
    # package treats as the canonical Naive CI.  For GLM outcomes we fall
    # back to make_effect_estimator()'s default SE.
    naive_est_se <- NA_real_
    naive_se     <- NA_real_
    if (outcome_type == "survival") {
      cox_res <- tryCatch({
        fit <- survival::coxph(
          survival::Surv(sub[[outcome.name]], sub[[event.name]]) ~
            sub[[treat.name]],
          robust = TRUE, model = FALSE, x = FALSE, y = FALSE
        )
        list(
          estimate  = as.numeric(stats::coef(fit))[1],
          se        = sqrt(stats::vcov(fit)[1, 1]),
          converged = TRUE
        )
      },
      error = function(e) list(estimate = NA_real_, se = NA_real_, converged = FALSE)
      )
      if (!isTRUE(cox_res$converged) ||
          is.na(cox_res$estimate) || is.na(cox_res$se)) {
        out[[k]] <- .empty_ci_row(row_k); next
      }
      naive_est_se <- cox_res$estimate
      naive_se     <- cox_res$se
    } else {
      full <- est_fn(sub)
      if (!isTRUE(full$converged) ||
          is.na(full$estimate) || is.na(full$se)) {
        out[[k]] <- .empty_ci_row(row_k); next
      }
      naive_est_se <- full$estimate
      naive_se     <- full$se
    }
    naive_lcl_se <- naive_est_se - z_crit * naive_se
    naive_ucl_se <- naive_est_se + z_crit * naive_se

    # --- 2b. Per-split half estimates -------------------------------------
    # The pressure-test revealed that AVERAGED halves
    # (theta_h1 + theta_h2) / 2 have near-zero variance: the two halves are
    # complementary, so for a near-linear statistic their sum is ~ 2*theta_full
    # by construction and the average is ~ theta_full regardless of split.
    #
    # The correct approach is the half-jackknife of Shao (1996): collect the
    # 2*S individual half-estimates and use their empirical SD directly as
    # the SE of the full-sample estimator.  No averaging, no sqrt(2)
    # correction, no IJ formula.
    half_estimates <- numeric(0L)
    n_valid_splits <- 0L
    for (s in seq_len(n_splits)) {
      in1 <- sample(c(TRUE, FALSE), n_sub, replace = TRUE,
                    prob = c(0.5, 0.5))
      h1 <- sub[in1, ,  drop = FALSE]
      h2 <- sub[!in1, , drop = FALSE]
      if (nrow(h1) < 5L || nrow(h2) < 5L) next
      r1 <- est_fn(h1)
      r2 <- est_fn(h2)
      if (!isTRUE(r1$converged) || !isTRUE(r2$converged)) next
      if (is.na(r1$estimate) || is.na(r2$estimate)) next
      n_valid_splits <- n_valid_splits + 1L
      half_estimates <- c(half_estimates, r1$estimate, r2$estimate)
    }

    # --- 2c. Compute Split CI and FSBC-mimic CI ---------------------------
    # Both rely on sd(half_estimates) as a half-jackknife SE estimate.
    # Bias-corrected estimator follows eq 7 of Leon et al. 2024 with the
    # eta(H_b*) selection term zeroed out:
    #   beta_FSBC = beta_hat - (1/(2S)) sum_h (beta_h - beta_hat)
    #             = 2 * beta_hat - mean(half_estimates).
    if (length(half_estimates) < 10L) {
      split_sd     <- NA_real_
      split_lcl_se <- NA_real_
      split_ucl_se <- NA_real_
      fsbc_est_se  <- NA_real_
      fsbc_lcl_se  <- NA_real_
      fsbc_ucl_se  <- NA_real_
    } else {
      split_sd     <- stats::sd(half_estimates)
      split_lcl_se <- naive_est_se - z_crit * split_sd
      split_ucl_se <- naive_est_se + z_crit * split_sd
      fsbc_est_se  <- 2 * naive_est_se - mean(half_estimates)
      fsbc_lcl_se  <- fsbc_est_se - z_crit * split_sd
      fsbc_ucl_se  <- fsbc_est_se + z_crit * split_sd
    }

    # --- 2d. Back-transform to natural scale for display ------------------
    if (is_log_scale) {
      est_nat       <- exp(naive_est_se)
      naive_lcl_nat <- exp(naive_lcl_se)
      naive_ucl_nat <- exp(naive_ucl_se)
      split_lcl_nat <- if (is.na(split_lcl_se)) NA_real_ else exp(split_lcl_se)
      split_ucl_nat <- if (is.na(split_ucl_se)) NA_real_ else exp(split_ucl_se)
      fsbc_est_nat  <- if (is.na(fsbc_est_se)) NA_real_ else exp(fsbc_est_se)
      fsbc_lcl_nat  <- if (is.na(fsbc_lcl_se)) NA_real_ else exp(fsbc_lcl_se)
      fsbc_ucl_nat  <- if (is.na(fsbc_ucl_se)) NA_real_ else exp(fsbc_ucl_se)
    } else {
      est_nat       <- naive_est_se
      naive_lcl_nat <- naive_lcl_se
      naive_ucl_nat <- naive_ucl_se
      split_lcl_nat <- split_lcl_se
      split_ucl_nat <- split_ucl_se
      fsbc_est_nat  <- fsbc_est_se
      fsbc_lcl_nat  <- fsbc_lcl_se
      fsbc_ucl_nat  <- fsbc_ucl_se
    }

    out[[k]] <- data.table::data.table(
      m              = as.integer(row_k$m),
      estimate       = est_nat,
      naive_lcl      = naive_lcl_nat,
      naive_ucl      = naive_ucl_nat,
      split_sd       = split_sd,
      split_lcl      = split_lcl_nat,
      split_ucl      = split_ucl_nat,
      fsbc_estimate  = fsbc_est_nat,
      fsbc_lcl       = fsbc_lcl_nat,
      fsbc_ucl       = fsbc_ucl_nat,
      n_valid_splits = n_valid_splits
    )

    if (verbose) {
      cat(sprintf(
        paste0("  Frontier row %d/%d (m=%d): naive (%.3f, %.3f), ",
               "split (%.3f, %.3f), FSBC est %.3f (%.3f, %.3f), n_valid=%d\n"),
        k, nrow(frontier), as.integer(row_k$m),
        naive_lcl_nat, naive_ucl_nat,
        split_lcl_nat, split_ucl_nat,
        fsbc_est_nat, fsbc_lcl_nat, fsbc_ucl_nat,
        n_valid_splits))
    }
  }

  result <- data.table::rbindlist(out, fill = TRUE)
  if (nrow(result) > 0L) data.table::setkey(result, m)
  result[]
}


# Internal helpers ------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

.empty_ci_row <- function(row_k) {
  data.table::data.table(
    m              = as.integer(row_k$m),
    estimate       = NA_real_,
    naive_lcl      = NA_real_,
    naive_ucl      = NA_real_,
    split_sd       = NA_real_,
    split_lcl      = NA_real_,
    split_ucl      = NA_real_,
    fsbc_estimate  = NA_real_,
    fsbc_lcl       = NA_real_,
    fsbc_ucl       = NA_real_,
    n_valid_splits = 0L
  )
}
