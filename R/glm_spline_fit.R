# =============================================================================
# glm_spline_fit.R
#
# Fit GLM with natural-spline treatment-biomarker interaction.
# Analog of cox_cs_fit() for binary, continuous, and count outcomes.
#
# Model:
#   g(E[Y|Z,A]) = beta_0 * A + f(Z) + g(Z) * A  [+ offset]
#
# where g() is the link function (logit for binary, log for Poisson,
# identity for Gaussian), f(Z) is a prognostic spline, and g(Z) is a
# treatment-interaction spline.
#
# The treatment-effect profile is:
#   beta(Z) = beta_0 + g(Z)
#
# which represents log-OR (binary/logistic), log-IRR (Poisson), or
# mean difference (Gaussian) as a function of the biomarker Z.
# =============================================================================


#' Fit GLM with Spline Treatment-Biomarker Interaction
#'
#' Estimates treatment effects as a function of a continuous covariate using
#' a generalized linear model with natural cubic splines.  This is the GLM
#' analog of \code{\link{cox_cs_fit}} for binary, continuous, and count
#' outcomes.
#'
#' @param df Data frame containing outcome data.
#' @param outcome_name Character. Name of the outcome variable.
#'   Default: \code{"event"}.
#' @param treat_name Character. Name of the treatment variable (0/1).
#'   Default: \code{"treat"}.
#' @param z_name Character. Name of the continuous biomarker covariate.
#'   Default: \code{"bm"}.
#' @param family Character or family object.  One of \code{"binomial"}
#'   (default), \code{"poisson"}, or \code{"gaussian"}.  Determines the
#'   link function and the scale of the treatment-effect profile:
#'   \code{"binomial"} uses logit link (profile is log-OR);
#'   \code{"poisson"} uses log link (profile is log-IRR, requires
#'   \code{offset_name}); \code{"gaussian"} uses identity link (profile
#'   is mean difference).
#' @param offset_name Character or \code{NULL}.  Name of the offset
#'   variable (e.g., log follow-up time for Poisson models).  The offset
#'   is applied as \code{log(offset_var)} internally.  Required when
#'   \code{family = "poisson"}.
#' @param alpha Numeric. Significance level for confidence intervals
#'   (two-sided).  Default: 0.05 (95 percent CI).
#' @param spline_df Integer. Degrees of freedom for the natural spline.
#'   Default: 3.
#' @param z_by Numeric. Increment for the biomarker prediction grid.
#'   Default: 0.05.
#' @param z_quantile Numeric (0--1). Upper quantile for the prediction
#'   grid.  Default: 0.95 (5th to 95th percentile).
#' @param z_max Numeric. Maximum z value for predictions.  Default:
#'   \code{Inf} (no truncation).
#' @param z_window Numeric. Half-width for counting observations near
#'   each grid point.  Default: 0.0.
#' @param show_plot Logical. Display a base-R diagnostic plot.
#'   Default: \code{FALSE}.
#' @param plot_params List. Optional plot parameter overrides (see
#'   \code{\link{cox_cs_fit}} for format).
#' @param truebeta_name Character or \code{NULL}. Name of a column
#'   containing true log-effect values (for simulation validation).
#' @param verbose Logical. Print diagnostic information.
#'   Default: \code{TRUE}.
#'
#' @return A list of class \code{"glm_cs_fit"} containing:
#' \describe{
#'   \item{z_profile}{Numeric vector. Biomarker grid values.}
#'   \item{loghr_est}{Numeric vector. Point estimates of the log
#'     treatment effect (log-OR, log-IRR, or MD) at each grid point.}
#'   \item{loghr_lower}{Numeric vector. Lower confidence bounds.}
#'   \item{loghr_upper}{Numeric vector. Upper confidence bounds.}
#'   \item{se_loghr}{Numeric vector. Standard errors (delta method).}
#'   \item{counts_profile}{Integer vector. Observation counts near each
#'     grid point.}
#'   \item{glm_primary}{Numeric. Overall treatment effect from the
#'     no-interaction GLM.}
#'   \item{model_fit}{The fitted \code{glm} object.}
#'   \item{spline_basis}{The natural spline basis object.}
#'   \item{family}{Character. The GLM family used.}
#'   \item{effect_label}{Character. Human-readable label for the effect
#'     scale (e.g., "log(OR)", "log(IRR)", "Mean Difference").}
#'   \item{lrt_pvalue}{Numeric. P-value from the likelihood ratio test
#'     comparing the interaction model to the main-effects-only model.}
#'   \item{alpha}{Numeric. Significance level used.}
#'   \item{ci_level}{Numeric. Confidence level (1 - alpha).}
#' }
#'
#' @details
#' ## Model Structure
#'
#' For binary outcomes (\code{family = "binomial"}):
#' \deqn{\text{logit}(P(Y=1|Z,A)) = \beta_0 A + f(Z) + g(Z) \cdot A}
#'
#' For count outcomes with offset (\code{family = "poisson"}):
#' \deqn{\log(E[Y|Z,A]) = \beta_0 A + f(Z) + g(Z) \cdot A + \log(t)}
#'
#' For continuous outcomes (\code{family = "gaussian"}):
#' \deqn{E[Y|Z,A] = \beta_0 A + f(Z) + g(Z) \cdot A}
#'
#' The treatment-effect profile \eqn{\beta(Z) = \beta_0 + g(Z)} is
#' estimated on the link scale (log-OR, log-IRR, or identity) with
#' pointwise delta-method confidence intervals.
#'
#' @seealso \code{\link{cox_cs_fit}} for the survival analog,
#'   \code{\link{glm_effect_profile}} for the extended interface with
#'   sandwich SEs and overdispersion correction.
#'
#' @examples
#' \dontrun{
#' # Binary outcome: log-OR profile over a biomarker
#' set.seed(42)
#' df <- data.frame(
#'   event = rbinom(500, 1, 0.3),
#'   treat = rbinom(500, 1, 0.5),
#'   bm    = rnorm(500, 0, 1)
#' )
#' result <- glm_cs_fit(df, outcome_name = "event", z_name = "bm",
#'                       family = "binomial")
#'
#' # Poisson outcome with offset: log-IRR profile
#' df$ftime <- rexp(500, 0.01) + 1
#' result_irr <- glm_cs_fit(df, outcome_name = "event", z_name = "bm",
#'                            family = "poisson", offset_name = "ftime")
#' }
#'
#' @export
#' @importFrom splines ns
#' @importFrom stats glm coef vcov qnorm quantile binomial poisson gaussian
#'   as.formula predict
#' @importFrom graphics plot lines abline axis text legend par box
glm_cs_fit <- function(df,
                       outcome_name = "event",
                       treat_name   = "treat",
                       z_name       = "bm",
                       family       = "binomial",
                       offset_name  = NULL,
                       alpha        = 0.05,
                       spline_df    = 3,
                       z_by         = 0.05,
                       z_quantile   = 0.95,
                       z_max        = Inf,
                       z_window     = 0.0,
                       show_plot    = FALSE,
                       plot_params  = NULL,
                       truebeta_name = NULL,
                       verbose      = TRUE) {

  # ==========================================================================
  # Input Validation
  # ==========================================================================

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar), add = TRUE)

  required_vars <- c(outcome_name, treat_name, z_name)
  if (!is.null(offset_name)) required_vars <- c(required_vars, offset_name)
  missing_vars <- setdiff(required_vars, names(df))
  if (length(missing_vars) > 0) {
    stop("Missing required variables: ",
         paste(missing_vars, collapse = ", "), call. = FALSE)
  }

  if (!is.null(truebeta_name) && !truebeta_name %in% names(df)) {
    stop("True beta variable '", truebeta_name, "' not found in data",
         call. = FALSE)
  }

  if (alpha <= 0 || alpha >= 1) {
    stop("alpha must be between 0 and 1", call. = FALSE)
  }

  if (spline_df < 1) {
    stop("spline_df must be at least 1", call. = FALSE)
  }

  # Resolve family
  family_name <- if (is.character(family)) family else family$family
  if (!family_name %in% c("binomial", "poisson", "gaussian")) {
    stop("family must be 'binomial', 'poisson', or 'gaussian'",
         call. = FALSE)
  }
  family_obj <- switch(family_name,
    binomial = stats::binomial(),
    poisson  = stats::poisson(),
    gaussian = stats::gaussian()
  )

  # Poisson requires offset

  if (family_name == "poisson" && is.null(offset_name)) {
    stop("family = 'poisson' requires offset_name (follow-up time column)",
         call. = FALSE)
  }

  # Effect label
  effect_label <- switch(family_name,
    binomial = "log(OR)",
    poisson  = "log(IRR)",
    gaussian = "Mean Difference"
  )

  # ==========================================================================
  # Extract Variables
  # ==========================================================================

  df <- df[order(df[[z_name]]), ]

  Y     <- df[[outcome_name]]
  Treat <- df[[treat_name]]
  z     <- df[[z_name]]

  if (!all(Treat %in% c(0, 1))) {
    warning("Treatment indicator should be 0 or 1")
  }

  if (family_name == "binomial" && !all(Y %in% c(0, 1))) {
    warning("Binary outcome should be 0 or 1")
  }

  # Offset (Poisson)
  offset_vec <- NULL
  if (!is.null(offset_name)) {
    offset_raw <- df[[offset_name]]
    if (any(offset_raw <= 0, na.rm = TRUE)) {
      warning("Non-positive offset values detected; adding small constant")
      offset_raw <- pmax(offset_raw, 1e-6)
    }
    offset_vec <- log(offset_raw)
  }

  # True beta (simulation)
  beta_true <- NULL
  if (!is.null(truebeta_name)) {
    beta_true <- df[[truebeta_name]]
  }

  # ==========================================================================
  # Define Prediction Grid
  # ==========================================================================

  z_lower <- quantile(z, probs = 1 - z_quantile, na.rm = TRUE)
  z_upper <- quantile(z, probs = z_quantile, na.rm = TRUE)
  z_profile <- seq(z_lower, z_upper, by = z_by)
  z_profile <- unique(pmin(z_profile, z_max))

  if (verbose) {
    cat("\n=== GLM Spline Fit ===\n")
    cat("Family:", family_name, " | Effect:", effect_label, "\n")
    cat("Z range:", round(range(z, na.rm = TRUE), 3), "\n")
    cat("Prediction grid:", length(z_profile), "points from",
        round(min(z_profile), 3), "to", round(max(z_profile), 3), "\n")
    cat("N =", nrow(df), " | Events =", sum(Y), "\n")
  }

  # ==========================================================================
  # Fit Primary (No Interaction) GLM
  # ==========================================================================

  if (!is.null(offset_vec)) {
    fit_primary <- stats::glm(Y ~ Treat, family = family_obj,
                              offset = offset_vec)
  } else {
    fit_primary <- stats::glm(Y ~ Treat, family = family_obj)
  }
  effect_primary <- coef(fit_primary)["Treat"]

  if (verbose) {
    cat("\n--- Primary GLM (no interaction) ---\n")
    cat(effect_label, ":", round(effect_primary, 4), "\n")
    if (family_name != "gaussian") {
      cat("Exponentiated:", round(exp(effect_primary), 4), "\n")
    }
  }

  # ==========================================================================
  # Create Spline Basis
  # ==========================================================================

  z_basis  <- splines::ns(z, df = spline_df)
  z_matrix <- as.matrix(z_basis)

  # Treatment-spline interactions
  z_treat_matrix <- Treat * z_matrix

  # Combined model matrix
  model_matrix <- cbind(Treat, z_matrix, z_treat_matrix)
  colnames(model_matrix) <- c(
    "treat",
    paste0("z_spline_", seq_len(spline_df)),
    paste0("treat_z_spline_", seq_len(spline_df))
  )

  # ==========================================================================
  # Fit GLM with Spline Interaction
  # ==========================================================================

  if (!is.null(offset_vec)) {
    fit_spline <- stats::glm(Y ~ model_matrix, family = family_obj,
                              offset = offset_vec)
  } else {
    fit_spline <- stats::glm(Y ~ model_matrix, family = family_obj)
  }

  beta_hat  <- coef(fit_spline)
  vcov_beta <- stats::vcov(fit_spline)

  # ==========================================================================
  # Fit Main-Effects-Only Model for LRT
  # ==========================================================================

  if (!is.null(offset_vec)) {
    fit_main <- stats::glm(Y ~ Treat + z_matrix, family = family_obj,
                            offset = offset_vec)
  } else {
    fit_main <- stats::glm(Y ~ Treat + z_matrix, family = family_obj)
  }

  lrt_stat  <- fit_main$deviance - fit_spline$deviance
  lrt_df    <- spline_df  # interaction terms added
  lrt_pval  <- stats::pchisq(lrt_stat, df = lrt_df, lower.tail = FALSE)

  if (verbose) {
    cat("\n--- LRT for treat x ns(", z_name, ") interaction ---\n")
    cat("Chi-sq =", round(lrt_stat, 2), ", df =", lrt_df,
        ", p =", formatC(lrt_pval, format = "g", digits = 4), "\n")
  }

  # ==========================================================================
  # Calculate Treatment-Effect Profile
  # ==========================================================================

  # Spline basis at profile z values
  z_basis_profile  <- predict(z_basis, newx = z_profile)
  z_matrix_profile <- as.matrix(z_basis_profile)

  # The treatment effect (on link scale) is:
  #   beta(Z) = beta_treat + sum(beta_treat_spline_k * spline_k(Z))
  # which is the difference in linear predictors: lp(A=1) - lp(A=0)

  # Indices: [1] = intercept, [2] = treat, [3..spline_df+2] = z_spline,
  #          [spline_df+3..2*spline_df+2] = treat_z_spline
  idx_treat    <- 2  # treat coefficient (offset by intercept)
  idx_interact <- (spline_df + 3):(2 * spline_df + 2)
  idx_diff     <- c(idx_treat, idx_interact)

  # Design matrix for the treatment effect difference
  diff_matrix <- cbind(
    rep(1, length(z_profile)),
    z_matrix_profile
  )

  # Point estimates
  loghr_z <- as.numeric(diff_matrix %*% beta_hat[idx_diff])

  # Delta-method standard errors
  vcov_diff <- vcov_beta[idx_diff, idx_diff]
  var_loghr <- diag(diff_matrix %*% vcov_diff %*% t(diff_matrix))
  se_loghr  <- sqrt(var_loghr)

  # Confidence intervals
  z_crit      <- qnorm(1 - alpha / 2)
  loghr_lower <- loghr_z - z_crit * se_loghr
  loghr_upper <- loghr_z + z_crit * se_loghr

  if (verbose) {
    cat("\n--- Treatment-Effect Profile ---\n")
    cat(effect_label, "range:", round(min(loghr_z), 4), "to",
        round(max(loghr_z), 4), "\n")
  }

  # ==========================================================================
  # Count Observations Near Each Profile Point
  # ==========================================================================

  if (z_window > 0) {
    counts_profile <- vapply(z_profile, function(z_val) {
      sum(z >= (z_val - z_window) & z <= (z_val + z_window), na.rm = TRUE)
    }, integer(1))
  } else {
    z_breaks <- c(z_profile - z_by / 2, max(z_profile) + z_by / 2)
    counts_profile <- as.vector(table(cut(z, breaks = z_breaks)))
    if (length(counts_profile) != length(z_profile)) {
      counts_profile <- rep(0L, length(z_profile))
    }
  }

  # ==========================================================================
  # Optional Base-R Plot
  # ==========================================================================

  if (show_plot) {
    default_pp <- list(
      xlab = z_name,
      ylab = effect_label,
      main_title = paste0("Treatment-Effect Profile (", effect_label, ")"),
      ylimit = NULL,
      y_pad_zero = 0.01,
      y_delta = 0.25,
      cex_legend = 0.9,
      show_primary = TRUE,
      show_null = TRUE
    )
    if (!is.null(plot_params)) {
      for (pn in names(plot_params)) default_pp[[pn]] <- plot_params[[pn]]
    }
    pp <- default_pp

    ymin <- min(loghr_lower, na.rm = TRUE)
    ymax <- max(loghr_upper, na.rm = TRUE)
    if (!is.null(beta_true)) {
      ymin <- min(ymin, beta_true, na.rm = TRUE)
      ymax <- max(ymax, beta_true, na.rm = TRUE)
    }
    if (is.null(pp$ylimit)) {
      pp$ylimit <- c(ymin - pp$y_delta, ymax + pp$y_delta)
    }

    null_ref <- if (family_name == "gaussian") 0 else 0  # log-scale: 0

    plot(z_profile, loghr_z,
         type = "l", lty = 1, col = "black", lwd = 3,
         xlab = pp$xlab, ylab = pp$ylab,
         ylim = pp$ylimit, main = pp$main_title)
    lines(z_profile, loghr_lower, lty = 2, col = "black", lwd = 1)
    lines(z_profile, loghr_upper, lty = 2, col = "black", lwd = 1)

    if (!is.null(beta_true)) {
      lines(z, beta_true, col = "grey", lwd = 3)
    }
    if (pp$show_null) abline(h = null_ref, col = "red", lty = 2, lwd = 0.5)
    if (pp$show_primary) {
      abline(h = effect_primary, col = "blue", lty = 2, lwd = 2)
    }

    ci_pct <- paste0(round((1 - alpha) * 100), "%")
    legend("topleft", horiz = TRUE,
           legend = c(effect_label, paste0(ci_pct, " CI")),
           lty = c(1, 2), lwd = c(3, 1), col = "black",
           bty = "n", cex = pp$cex_legend)
  }

  # ==========================================================================
  # Return
  # ==========================================================================

  results <- list(
    z_profile      = z_profile,
    loghr_est      = loghr_z,
    loghr_lower    = loghr_lower,
    loghr_upper    = loghr_upper,
    se_loghr       = se_loghr,
    counts_profile = counts_profile,
    glm_primary    = effect_primary,
    model_fit      = fit_spline,
    model_main     = fit_main,
    spline_basis   = z_basis,
    family         = family_name,
    effect_label   = effect_label,
    lrt_pvalue     = lrt_pval,
    alpha          = alpha,
    ci_level       = 1 - alpha,
    z_name         = z_name,
    treat_name     = treat_name,
    outcome_name   = outcome_name,
    offset_name    = offset_name
  )

  class(results) <- c("glm_cs_fit", "list")
  return(results)
}


#' Print method for glm_cs_fit objects
#'
#' @param x A \code{glm_cs_fit} object.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns \code{x}.
#' @rdname glm_cs_fit
#' @export
print.glm_cs_fit <- function(x, ...) {
  cat("GLM Spline Fit (", x$family, " / ", x$effect_label, ")\n", sep = "")
  cat("  Biomarker:", x$z_name, "\n")
  cat("  Profile points:", length(x$z_profile), "\n")
  cat("  Effect range:", round(min(x$loghr_est), 3), "to",
      round(max(x$loghr_est), 3), "\n")
  cat("  Overall effect:", round(x$glm_primary, 3), "\n")
  cat("  LRT interaction p:", formatC(x$lrt_pvalue, format = "g",
                                       digits = 4), "\n")
  cat("  CI level:", x$ci_level * 100, "%\n")
  invisible(x)
}
