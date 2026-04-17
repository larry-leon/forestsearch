# =============================================================================
# Nonparametric Tests for Treatment Effect Heterogeneity
# =============================================================================
#
# Crump, R.K., Hotz, V.J., Imbens, G.W., Mitnik, O.A. (2008).
# "Nonparametric Tests for Treatment Effect Heterogeneity."
# Review of Economics and Statistics, 90(3), 389-405.
#
# Extensions: GLM-based tests (logistic, Poisson+offset),
# variance diagnostic output, selected-covariate reporting.
# =============================================================================

#' Build Polynomial Basis Matrix
#' @param X Numeric matrix (N x d). Must be numeric (no factors).
#' @param poly_order Integer. 1 = linear, 2 = squares + interactions.
#' @param include_intercept Logical. Default: TRUE.
#' @return Numeric matrix (N x K).
#' @importFrom stats glm.fit var binomial
#' @keywords internal
build_basis <- function(X, poly_order = 1L,
                        include_intercept = TRUE) {
  X <- as.matrix(X)
  n <- nrow(X); d <- ncol(X)
  X_std <- scale(X)
  sc <- attr(X_std, "scaled:scale")
  sc[is.na(sc) | sc == 0] <- 1
  X_std <- scale(X, center = attr(X_std, "scaled:center"),
                 scale = sc)
  basis <- list()
  if (include_intercept) basis[["intercept"]] <- rep(1, n)
  for (j in seq_len(d)) basis[[paste0("x", j)]] <- X_std[, j]
  if (poly_order >= 2) {
    for (j in seq_len(d))
      basis[[paste0("x", j, "_sq")]] <- X_std[, j]^2
    if (d >= 2)
      for (j in seq_len(d - 1))
        for (k in (j + 1):d)
          basis[[paste0("x", j, "_x", k)]] <-
            X_std[, j] * X_std[, k]
  }
  do.call(cbind, basis)
}


#' OLS Sandwich Variance (HC0)
#' @keywords internal
sandwich_var_ols <- function(P, Y, beta_hat) {
  N <- nrow(P)
  resid <- as.numeric(Y - P %*% beta_hat)
  PtP_inv <- tryCatch(solve(crossprod(P)),
    error = function(e) MASS::ginv(crossprod(P)))
  N * PtP_inv %*% crossprod(P * resid) %*% PtP_inv
}


#' GLM Sandwich Variance
#' @keywords internal
sandwich_var_glm <- function(P, Y, mu, family) {
  N <- nrow(P)
  dmu_deta <- family$mu.eta(family$linkfun(mu))
  var_mu <- family$variance(mu)
  var_mu[var_mu < 1e-10] <- 1e-10
  w <- dmu_deta^2 / var_mu
  score <- P * as.numeric(dmu_deta * (Y - mu) / var_mu)
  PwP_inv <- tryCatch(solve(crossprod(P * sqrt(w))),
    error = function(e) MASS::ginv(crossprod(P * sqrt(w))))
  N * PwP_inv %*% crossprod(score) %*% PwP_inv
}


#' Fit Separate Regressions by Treatment Arm
#' @keywords internal
.fit_arms <- function(Y, W, X, poly_order, covariate_select,
                      t_threshold, regression = "ols",
                      offset = NULL) {
  X <- as.matrix(X)
  W <- as.integer(W)
  stopifnot(length(Y) == length(W), length(Y) == nrow(X))

  sel_result <- .select_covariates(Y, W, X, covariate_select,
                                    t_threshold)
  X_sel <- sel_result$X
  selected_names <- sel_result$selected_names

  P_full <- build_basis(X_sel, poly_order = poly_order)

  # Rank check: drop linearly dependent columns
  qr_full <- qr(P_full)
  if (qr_full$rank < ncol(P_full)) {
    keep <- qr_full$pivot[seq_len(qr_full$rank)]
    P_full <- P_full[, keep, drop = FALSE]
  }
  K <- ncol(P_full)

  idx0 <- which(W == 0); idx1 <- which(W == 1)
  N0 <- length(idx0); N1 <- length(idx1)
  P0 <- P_full[idx0, , drop = FALSE]
  P1 <- P_full[idx1, , drop = FALSE]
  Y0 <- Y[idx0]; Y1 <- Y[idx1]
  off0 <- if (!is.null(offset)) offset[idx0] else NULL
  off1 <- if (!is.null(offset)) offset[idx1] else NULL

  if (regression == "ols") {
    beta0 <- tryCatch(
      as.numeric(solve(crossprod(P0), crossprod(P0, Y0))),
      error = function(e) as.numeric(
        MASS::ginv(crossprod(P0)) %*% crossprod(P0, Y0)))
    beta1 <- tryCatch(
      as.numeric(solve(crossprod(P1), crossprod(P1, Y1))),
      error = function(e) as.numeric(
        MASS::ginv(crossprod(P1)) %*% crossprod(P1, Y1)))
    Omega0 <- sandwich_var_ols(P0, Y0, beta0)
    Omega1 <- sandwich_var_ols(P1, Y1, beta1)
  } else {
    fam <- if (regression == "logistic") binomial("logit")
           else poisson("log")
    fit0 <- tryCatch(
      glm.fit(P0, Y0, family = fam, offset = off0,
              intercept = FALSE),
      error = function(e) NULL)
    fit1 <- tryCatch(
      glm.fit(P1, Y1, family = fam, offset = off1,
              intercept = FALSE),
      error = function(e) NULL)
    if (is.null(fit0) || is.null(fit1)) return(NULL)
    if (!fit0$converged || !fit1$converged) return(NULL)
    beta0 <- as.numeric(fit0$coefficients)
    beta1 <- as.numeric(fit1$coefficients)
    if (any(is.na(beta0)) || any(is.na(beta1))) return(NULL)
    Omega0 <- sandwich_var_glm(P0, Y0, fit0$fitted.values, fam)
    Omega1 <- sandwich_var_glm(P1, Y1, fit1$fitted.values, fam)
  }

  V_hat <- Omega0 / N0 + Omega1 / N1
  list(beta0 = beta0, beta1 = beta1, V_hat = V_hat,
       K = K, N0 = N0, N1 = N1,
       selected_names = selected_names)
}


#' Test for Zero Conditional Average Treatment Effect
#'
#' H0: tau(x) = 0 for all x.
#'
#' @param Y Numeric outcome vector.
#' @param W Treatment indicator (0/1).
#' @param X Numeric covariate matrix (no factors; use
#'   \code{as.numeric(factor) - 1} to convert).
#' @param poly_order Integer. Default: 1.
#' @param covariate_select "all", "top_down", or "bottom_up".
#' @param t_threshold Numeric. Default: 2.0.
#' @param regression "ols", "logistic", or "poisson".
#' @param offset Numeric offset vector (Poisson only).
#' @return List with test results including \code{covariates_selected}.
#' @export
test_zero_cate <- function(Y, W, X, poly_order = 1L,
                           covariate_select = c("all", "top_down",
                                                 "bottom_up"),
                           t_threshold = 2.0,
                           regression = c("ols", "logistic",
                                           "poisson"),
                           offset = NULL) {
  covariate_select <- match.arg(covariate_select)
  regression <- match.arg(regression)
  arms <- .fit_arms(Y, W, X, poly_order, covariate_select,
                     t_threshold, regression, offset)
  if (is.null(arms))
    return(list(test = "Zero CATE", chi_sq = NA, df = NA,
                p_value_chi = NA, normal = NA,
                p_value_normal = NA, regression = regression,
                K = NA, N0 = NA, N1 = NA,
                covariates_selected = character(0)))
  K <- arms$K
  diff <- arms$beta1 - arms$beta0
  V_inv <- tryCatch(solve(arms$V_hat),
    error = function(e) MASS::ginv(arms$V_hat))
  Q <- as.numeric(t(diff) %*% V_inv %*% diff)
  p_chi <- pchisq(Q, df = K, lower.tail = FALSE)
  T_n <- (Q - K) / sqrt(2 * K)
  list(test = "Zero CATE", chi_sq = Q, df = K,
       p_value_chi = p_chi, normal = T_n,
       p_value_normal = pnorm(T_n, lower.tail = FALSE),
       K = K, N0 = arms$N0, N1 = arms$N1,
       regression = regression,
       covariates_selected = arms$selected_names,
       coefs_0 = arms$beta0, coefs_1 = arms$beta1,
       diff = diff, V_hat = arms$V_hat)
}


#' Test for Constant Conditional Average Treatment Effect
#'
#' H0': tau(x) = tau for some tau and all x.
#'
#' @inheritParams test_zero_cate
#' @return A list with components \code{test}, \code{chi_sq}, \code{df},
#'   \code{p_value_chi}, \code{normal} (normal-approximation statistic),
#'   \code{p_value_normal}, \code{K} (number of covariates), \code{N0},
#'   \code{N1}, \code{regression}, \code{covariates_selected},
#'   \code{ate_diff} (intercept-difference point estimate), and
#'   \code{diff_slope}, \code{V_slope}.  Returns \code{NULL} if per-arm
#'   fits fail or fewer than two covariates are retained.
#' @export
test_constant_cate <- function(Y, W, X, poly_order = 1L,
                               covariate_select = c("all",
                                 "top_down", "bottom_up"),
                               t_threshold = 2.0,
                               regression = c("ols", "logistic",
                                               "poisson"),
                               offset = NULL) {
  covariate_select <- match.arg(covariate_select)
  regression <- match.arg(regression)
  arms <- .fit_arms(Y, W, X, poly_order, covariate_select,
                     t_threshold, regression, offset)
  if (is.null(arms)) return(NULL)
  K <- arms$K
  if (K < 2) return(NULL)
  idx <- 2:K
  d_slope <- arms$beta1[idx] - arms$beta0[idx]
  V_s <- arms$V_hat[idx, idx, drop = FALSE]
  V_s_inv <- tryCatch(solve(V_s),
    error = function(e) MASS::ginv(V_s))
  Q <- as.numeric(t(d_slope) %*% V_s_inv %*% d_slope)
  df <- K - 1
  p_chi <- pchisq(Q, df = df, lower.tail = FALSE)
  T_n <- (Q - df) / sqrt(2 * df)
  list(test = "Constant CATE", chi_sq = Q, df = df,
       p_value_chi = p_chi, normal = T_n,
       p_value_normal = pnorm(T_n, lower.tail = FALSE),
       K = K, N0 = arms$N0, N1 = arms$N1,
       regression = regression,
       covariates_selected = arms$selected_names,
       ate_diff = arms$beta1[1] - arms$beta0[1],
       diff_slope = d_slope, V_slope = V_s)
}


#' Run Both Crump et al. Tests
#' @inheritParams test_zero_cate
#' @return Object of class \code{hte_test}.
#' @export
test_hte <- function(Y, W, X, poly_order = 1L,
                     covariate_select = c("all", "top_down",
                                           "bottom_up"),
                     t_threshold = 2.0,
                     regression = c("ols", "logistic", "poisson"),
                     offset = NULL) {
  covariate_select <- match.arg(covariate_select)
  regression <- match.arg(regression)
  t0 <- test_zero_cate(Y, W, X, poly_order, covariate_select,
                        t_threshold, regression, offset)
  t1 <- test_constant_cate(Y, W, X, poly_order,
                            covariate_select, t_threshold,
                            regression, offset)
  W <- as.integer(W)
  ate <- mean(Y[W == 1]) - mean(Y[W == 0])
  se <- sqrt(var(Y[W == 1]) / sum(W == 1) +
              var(Y[W == 0]) / sum(W == 0))
  z <- ate / se
  t2 <- list(test = "Zero ATE", ate = ate, se = se,
             chi_sq = z^2, df = 1,
             p_value = 2 * pnorm(abs(z), lower.tail = FALSE),
             normal = z,
             p_value_normal = 2 * pnorm(abs(z),
                                         lower.tail = FALSE))
  result <- list(zero_cate = t0, constant_cate = t1,
                 zero_ate = t2)
  class(result) <- "hte_test"
  result
}


#' Print method for \code{hte_test} objects
#' @param x An \code{hte_test} object.
#' @param ... Unused; present for S3 compatibility.
#' @return The input \code{x}, invisibly.
#' @export
print.hte_test <- function(x, ...) {
  cat("=== Crump et al. (2008) HTE Tests ===\n")
  reg <- x$zero_cate$regression
  if (!is.null(reg) && reg != "ols")
    cat(sprintf("    Regression: %s\n", reg))
  covs <- x$zero_cate$covariates_selected
  if (!is.null(covs) && length(covs) > 0)
    cat(sprintf("    Covariates: %s\n", paste(covs, collapse = ", ")))
  cat("\n")
  t0 <- x$zero_cate
  cat(sprintf("1. Zero Cond. Ave. TE:     chi-sq=%6.2f (df=%d) p=%.4f | z=%6.2f p=%.4f\n",
              t0$chi_sq, t0$df, t0$p_value_chi,
              t0$normal, t0$p_value_normal))
  t1 <- x$constant_cate
  if (!is.null(t1))
    cat(sprintf("2. Constant Cond. Ave. TE: chi-sq=%6.2f (df=%d) p=%.4f | z=%6.2f p=%.4f\n",
                t1$chi_sq, t1$df, t1$p_value_chi,
                t1$normal, t1$p_value_normal))
  t2 <- x$zero_ate
  cat(sprintf("3. Zero Ave. TE:           chi-sq=%6.2f (df=%d) p=%.4f | z=%6.2f p=%.4f\n",
              t2$chi_sq, t2$df, t2$p_value,
              t2$normal, t2$p_value_normal))
  cat(sprintf("\nN0=%d, N1=%d, K=%d\n", t0$N0, t0$N1, t0$K))
  invisible(x)
}


# ---- Covariate selection (internal) ----

#' Select Covariates for the Crump Test
#'
#' Returns a list with X (selected columns) and selected_names.
#' @keywords internal
.select_covariates <- function(Y, W, X, method, threshold) {
  X <- as.matrix(X)
  all_names <- colnames(X)
  if (is.null(all_names))
    all_names <- paste0("V", seq_len(ncol(X)))

  if (method == "all")
    return(list(X = X, selected_names = all_names))

  idx0 <- which(W == 0)
  Y0 <- Y[idx0]; X0 <- X[idx0, , drop = FALSE]; d <- ncol(X)

  if (method == "top_down") {
    sel <- seq_len(d)
    repeat {
      if (length(sel) == 0) break
      P0 <- cbind(1, X0[, sel, drop = FALSE])
      fit <- tryCatch({
        b <- solve(crossprod(P0), crossprod(P0, Y0))
        r <- Y0 - P0 %*% b
        se <- sqrt(diag(solve(crossprod(P0)) *
                         sum(r^2) / (length(Y0) - ncol(P0))))
        abs(b / se)[-1]
      }, error = function(e) rep(Inf, length(sel)))
      mi <- which.min(fit)
      if (fit[mi] >= threshold) break
      sel <- sel[-mi]
    }
    return(list(X = X[, sel, drop = FALSE],
                selected_names = all_names[sel]))
  }

  if (method == "bottom_up") {
    sel <- integer(0); rem <- seq_len(d)
    repeat {
      if (length(rem) == 0) break
      best_t <- -Inf; best_j <- NULL
      for (j in rem) {
        cands <- c(sel, j)
        P0 <- cbind(1, X0[, cands, drop = FALSE])
        fit <- tryCatch({
          b <- solve(crossprod(P0), crossprod(P0, Y0))
          r <- Y0 - P0 %*% b
          se <- sqrt(diag(solve(crossprod(P0)) *
                           sum(r^2) / (length(Y0) - ncol(P0))))
          abs(b[length(b)] / se[length(se)])
        }, error = function(e) 0)
        if (fit > best_t) { best_t <- fit; best_j <- j }
      }
      if (best_t < threshold) break
      sel <- c(sel, best_j); rem <- setdiff(rem, best_j)
    }
    if (length(sel) == 0)
      return(list(X = X, selected_names = all_names))
    return(list(X = X[, sel, drop = FALSE],
                selected_names = all_names[sel]))
  }

  list(X = X, selected_names = all_names)
}


# =============================================================================
# ForestSearch Integration
# =============================================================================

#' HTE Test Using ForestSearch-Selected Covariates
#'
#' @param fs_result A \code{forestsearch} result object.
#' @param df Data frame used in the \code{forestsearch()} call.
#' @param outcome.name Character. Outcome variable name.
#' @param treat.name Character. Treatment variable name.
#' @param regression Character: "ols", "logistic", or "poisson".
#' @param offset Numeric offset vector (Poisson only).
#' @param poly_order Integer. Default: 1.
#' @return List with hte_test, covariates_used, source.
#' @export
test_hte_from_forestsearch <- function(fs_result, df,
                                        outcome.name,
                                        treat.name,
                                        regression = c("ols",
                                          "logistic", "poisson"),
                                        offset = NULL,
                                        poly_order = 1L) {
  regression <- match.arg(regression)
  covs <- NULL; source_desc <- NULL
  if (!is.null(fs_result$confounders.evaluated) &&
      length(fs_result$confounders.evaluated) > 0) {
    covs <- fs_result$confounders.evaluated
    source_desc <- "confounders.evaluated (GRF/LASSO selected)"
  } else if (!is.null(fs_result$confounders.candidate) &&
             length(fs_result$confounders.candidate) > 0) {
    covs <- fs_result$confounders.candidate
    source_desc <- "confounders.candidate (all candidates)"
  }
  if (is.null(covs)) { warning("No covariates found."); return(NULL) }

  X_cols <- list()
  for (v in covs) {
    col <- df[[v]]
    if (is.null(col)) next
    if (is.factor(col)) X_cols[[v]] <- as.numeric(col) - 1
    else if (is.character(col)) X_cols[[v]] <- as.numeric(as.factor(col)) - 1
    else X_cols[[v]] <- as.numeric(col)
  }
  if (length(X_cols) == 0) return(NULL)

  X <- do.call(cbind, X_cols)
  colnames(X) <- names(X_cols)
  hte <- test_hte(df[[outcome.name]], df[[treat.name]], X,
                   poly_order = poly_order,
                   regression = regression, offset = offset)
  n_cands <- NA
  if (!is.null(fs_result$find.grps$max_sg_est))
    n_cands <- length(fs_result$find.grps$max_sg_est)

  list(hte_test = hte, covariates_used = names(X_cols),
       n_covariates = ncol(X), n_candidates = n_cands,
       source = source_desc)
}
