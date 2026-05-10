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
#'   \item \strong{Split CI}: a 95\% interval derived from the variance of
#'     averaged half-sample effect estimates over \code{n_splits} random
#'     50/50 splits.  For each split, the per-half effects are averaged,
#'     and the empirical SD of those averaged values (computed on the log
#'     scale for ratio measures) approximates the SE of the full-sample
#'     estimate.  See \dQuote{Details — Why averaged halves work} below.
#'   \item \strong{FSBC-mimic CI}: a bias-corrected interval based on the
#'     bootstrap bias-correction algorithm of \cite{Leon2024fs} (eq 7
#'     and 9) with the \eqn{\eta_b^*(\hat H_b^*)} term zeroed out
#'     (treating the selected subgroup as fixed across splits).  Each
#'     averaged-halves estimator \eqn{\tilde\beta_s} plays the role of a
#'     bootstrap replicate; the bias correction is
#'     \eqn{\hat\beta^{\mathrm{FSBC}} = \hat\beta - (1/S)\sum_s \eta_s
#'     = 2\hat\beta - \bar{\tilde\beta}}, and the variance is the
#'     infinitesimal-jackknife formula based on the realized
#'     half-assignment matrix.  See \dQuote{Details — FSBC-mimic
#'     interpretation} below.
#' }
#'
#' Both CIs are computed post-hoc on the returned \code{forestsearch}
#' object; this function does not modify any internal state.
#'
#' @section Details — Why averaged halves \enquote{work}:
#' Let \eqn{\hat\theta^{(s,1)}} and \eqn{\hat\theta^{(s,2)}} be the
#' per-half effect estimates on disjoint random halves of the data, each
#' having variance \eqn{\approx 2\sigma^2} (relative to a full-sample
#' estimate's \eqn{\sigma^2}).  The averaged value
#' \eqn{\tilde\theta^{(s)} = \tfrac{1}{2}(\hat\theta^{(s,1)} + \hat\theta^{(s,2)})}
#' has variance \eqn{\tfrac{1}{4}(2\sigma^2 + 2\sigma^2) = \sigma^2},
#' since the two halves are disjoint and approximately independent.  Thus
#' \eqn{\mathrm{sd}\{\tilde\theta^{(s)}\}} approximates the SE of the
#' full-sample estimator without requiring an explicit \eqn{\sqrt{2}}
#' correction.
#'
#' @section Details — FSBC-mimic interpretation:
#' The FSBC-mimic CI is a \emph{quick approximation} to the bootstrap
#' bias-corrected estimator of \cite{Leon2024fs}.  Substituting splits
#' \eqn{s=1,\ldots,S} for bootstrap replicates \eqn{b=1,\ldots,B}, and
#' treating the selected subgroup \eqn{\hat H} as fixed (i.e., omitting
#' the \eqn{\eta_b^*(\hat H_b^*) = \hat\beta_b^*(\hat H_b^*) - \hat\beta(\hat H_b^*)}
#' term that captures variability in FS selection on the bootstrap
#' data):
#' \deqn{\hat\beta^{\mathrm{FSBC}}(\hat H) = \hat\beta(\hat H) - \frac{1}{S}\sum_s \eta_s}
#' where \eqn{\eta_s = \tilde\beta_s - \hat\beta(\hat H)}.  This collapses
#' to \eqn{2\hat\beta(\hat H) - \bar{\tilde\beta}}.  The IJ variance
#' (\cite{Leon2024fs} eq 9) with the same omission reduces to
#' \deqn{\widetilde V = \sum_{i=1}^{n_\mathrm{sub}} \tilde c_i^2,
#' \quad \tilde c_i = \frac{1}{S}\sum_s (K_{si} - \bar K_i)(\bar{\tilde\beta} - \tilde\beta_s),}
#' \deqn{\hat V = \widetilde V - (n_\mathrm{sub}/S)\,\widetilde\sigma^2,
#' \quad \widetilde\sigma^2 = \frac{1}{S}\sum_s (\bar{\tilde\beta} - \tilde\beta_s)^2.}
#' Here \eqn{K_{si} \in \{0,1\}} indicates subject \eqn{i}'s half
#' assignment in split \eqn{s}.  When \eqn{\hat V \le 0} (which can
#' happen because of Monte Carlo bias correction in a low-replicate
#' regime), the CI is reported with \eqn{V_\mathrm{used} = 0} — the
#' result is degenerate but informative as a flag.  The
#' \code{fsbc_var_pos} column indicates whether the bias-corrected
#' variance was positive.
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
#'     \item{\code{split_sd}}{Empirical SD of averaged-half effects on
#'       the SE scale (log for ratio measures, natural for differences).}
#'     \item{\code{split_lcl}, \code{split_ucl}}{Split-derived CI bounds.}
#'     \item{\code{fsbc_estimate}}{Bias-corrected effect estimate on the
#'       natural scale (\eqn{2\hat\beta - \bar{\tilde\beta}} on the SE
#'       scale, then back-transformed).}
#'     \item{\code{fsbc_lcl}, \code{fsbc_ucl}}{FSBC-mimic CI bounds.}
#'     \item{\code{fsbc_var_pos}}{Logical; \code{TRUE} when the
#'       bias-corrected variance \eqn{\hat V} was positive.  When
#'       \code{FALSE}, the CI is degenerate (\eqn{V} clipped to 0) and
#'       should be interpreted cautiously.}
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

    # Cuts for this frontier member: M.* values, excluding empty/NA
    cuts <- unlist(row_k[, m_cols, with = FALSE], use.names = FALSE)
    cuts <- cuts[!is.na(cuts) & nzchar(cuts)]

    # Identify the subgroup on df.est:
    # each cut name is itself a binary column on df.est
    keep <- rep(TRUE, nrow(df.est))
    missing_cuts <- character(0)
    for (c_name in cuts) {
      if (!c_name %in% names(df.est)) {
        missing_cuts <- c(missing_cuts, c_name)
        next
      }
      col <- df.est[[c_name]]
      keep <- keep & !is.na(col) & col == 1L
    }
    if (length(missing_cuts) > 0L) {
      warning(sprintf(
        "compute_frontier_cis(): frontier row %d references cut(s) not found on df.est: %s. Returning NA for this row.",
        k, paste(missing_cuts, collapse = ", ")), call. = FALSE)
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
    full <- est_fn(sub)
    if (!isTRUE(full$converged) || is.na(full$estimate) || is.na(full$se)) {
      out[[k]] <- .empty_ci_row(row_k)
      next
    }
    naive_est_se <- full$estimate       # SE scale (log for ratio measures)
    naive_se     <- full$se
    naive_lcl_se <- naive_est_se - z_crit * naive_se
    naive_ucl_se <- naive_est_se + z_crit * naive_se

    # --- 2b. Split-derived CI: per-split averaged halves ------------------
    # We collect TWO things per valid split:
    #   averaged[s] = (theta_hat_h1 + theta_hat_h2) / 2     (averaged-half estimator)
    #   K_list[[s]] = logical vector of length n_sub indicating half-1 membership
    # The K_list is needed for the FSBC-mimic IJ variance (Section 2c).
    averaged <- numeric(0L)
    K_list   <- vector("list", n_splits)
    n_valid  <- 0L
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
      n_valid <- n_valid + 1L
      averaged[n_valid] <- (r1$estimate + r2$estimate) / 2
      K_list[[n_valid]] <- as.integer(in1)
    }

    if (length(averaged) < 10L) {
      split_sd <- NA_real_
      split_lcl_se <- NA_real_
      split_ucl_se <- NA_real_
    } else {
      split_sd     <- stats::sd(averaged)
      split_lcl_se <- naive_est_se - z_crit * split_sd
      split_ucl_se <- naive_est_se + z_crit * split_sd
    }

    # --- 2c. FSBC-mimic CI: bias-corrected + IJ-style variance ------------
    # Treats each averaged-half estimator tilde_beta_s as a "split-bootstrap"
    # replicate; eta_s = tilde_beta_s - beta_hat (the per-split bias term);
    # the bias-correction follows paper eq (7) with the eta(H_b*) term
    # zeroed out (H is treated as fixed across splits):
    #   beta_FSBC = beta_hat - mean(eta_s) = 2*beta_hat - mean(tilde_beta_s)
    # Variance is the paper's IJ formula (eq 9) with the same omission:
    #   c_tilde_i = mean_s (K_si - mean_s K_si) * (mean(tilde) - tilde_s)
    #   V_tilde   = sum_i c_tilde_i^2
    #   V_hat     = V_tilde - (n_sub / S) * mean_s ((mean(tilde) - tilde_s)^2)
    # All quantities on the SE scale (log for ratio measures, natural for
    # difference measures); back-transformed in Section 2d below.
    if (length(averaged) < 10L) {
      fsbc_est_se  <- NA_real_
      fsbc_lcl_se  <- NA_real_
      fsbc_ucl_se  <- NA_real_
      fsbc_var_pos <- NA
    } else {
      K_mat <- do.call(rbind, K_list[seq_len(n_valid)])  # n_valid x n_sub
      mean_tilde <- mean(averaged)
      r_vec      <- mean_tilde - averaged                # length n_valid
      Kbar       <- colMeans(K_mat)                      # length n_sub
      # Per-subject covariance: c_i = (1/S) sum_s (K_si - Kbar_i) * r_s
      Kc         <- sweep(K_mat, 2L, Kbar, FUN = "-")    # n_valid x n_sub centered
      c_vec      <- as.numeric(crossprod(Kc, r_vec)) / n_valid  # length n_sub
      V_tilde    <- sum(c_vec * c_vec)
      sigma2_S   <- mean(r_vec * r_vec)
      V_hat      <- V_tilde - (n_sub / n_valid) * sigma2_S
      fsbc_var_pos <- isTRUE(V_hat > 0)
      V_use      <- max(V_hat, 0)
      fsbc_se    <- sqrt(V_use)
      fsbc_est_se <- 2 * naive_est_se - mean_tilde
      fsbc_lcl_se <- fsbc_est_se - z_crit * fsbc_se
      fsbc_ucl_se <- fsbc_est_se + z_crit * fsbc_se
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
      fsbc_var_pos   = fsbc_var_pos,
      n_valid_splits = n_valid
    )

    if (verbose) {
      cat(sprintf(
        paste0("  Frontier row %d/%d (m=%d): naive (%.3f, %.3f), ",
               "split (%.3f, %.3f), FSBC est %.3f (%.3f, %.3f) [var_pos=%s], ",
               "n_valid=%d\n"),
        k, nrow(frontier), as.integer(row_k$m),
        naive_lcl_nat, naive_ucl_nat,
        split_lcl_nat, split_ucl_nat,
        fsbc_est_nat, fsbc_lcl_nat, fsbc_ucl_nat,
        as.character(fsbc_var_pos), n_valid))
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
    fsbc_var_pos   = NA,
    n_valid_splits = 0L
  )
}
