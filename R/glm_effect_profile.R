# =============================================================================
# glm_effect_profile.R
# GLM delta-method treatment-effect profile for continuous biomarker
# interaction. Supports binary (log-OR, log-RR, RD), continuous (MD),
# and count/rate (log-IRR with optional person-time offset, quasi-Poisson,
# neg-binomial).
#
# Distinct from glm_cs_fit() in glm_spline_fit.R:
#   - glm_cs_fit()          : family-based interface, LRT, base-R plot
#   - glm_effect_profile()  : effect_measure-based interface, sandwich SEs,
#                             overdispersion (quasi / negbin), strata support
# =============================================================================

#' GLM Treatment-Effect Profile for Continuous Biomarker Interaction
#'
#' Estimates a treatment effect profile across a continuous biomarker using
#' a generalized linear model with natural cubic spline interactions.
#' Supports binary, continuous, and count or rate outcomes. Standard errors
#' are computed via the delta method applied to the full model covariance,
#' with optional sandwich variance estimation for modified-Poisson models.
#'
#' @param df Data frame containing the analysis data.
#' @param outcome.name Character. Name of the outcome variable. Binary
#'   outcomes must be coded 0/1; count outcomes must be non-negative integers;
#'   continuous outcomes may be any numeric.
#' @param treat.name Character. Name of the binary treatment indicator
#'   (0 = control, 1 = treated). Default: \code{"treat"}.
#' @param z.name Character. Name of the continuous biomarker (effect modifier).
#'   Default: \code{"bm"}.
#' @param effect_measure Character. Effect measure to estimate on the link
#'   scale. Allowed values: \code{"log_OR"} (log odds ratio),
#'   \code{"log_RR"} (log relative risk, modified Poisson),
#'   \code{"log_IRR"} (log incidence rate ratio, requires offset),
#'   \code{"RD"} (risk difference), \code{"MD"} (mean difference).
#'   Default: \code{"log_OR"}.
#' @param offset.name Character or \code{NULL}. Column name for person-time
#'   or exposure. When provided, \code{log(offset)} enters the linear
#'   predictor; required for \code{effect_measure = "log_IRR"}.
#'   Default: \code{NULL}.
#' @param overdispersion Character. Overdispersion correction for Poisson
#'   family models. \code{"none"} uses standard Poisson (model-based SEs or
#'   sandwich if \code{effect_measure = "log_RR"}); \code{"quasi"} uses
#'   quasi-Poisson; \code{"negbin"} uses negative binomial via
#'   \code{MASS::glm.nb} (requires \pkg{MASS}).
#'   Default: \code{"none"}.
#' @param spline_df Integer. Degrees of freedom for natural spline basis.
#'   Default: \code{3}.
#' @param z_by Numeric. Increment for the biomarker prediction grid.
#'   Default: \code{0.05}.
#' @param z_quantile Numeric. Upper quantile of the biomarker used as the
#'   grid endpoint (avoids extrapolation into sparse tails).
#'   Default: \code{0.95}.
#' @param alpha Numeric. Two-sided significance level for confidence
#'   intervals. Default: \code{0.05}.
#' @param conf.level Numeric or \code{NULL}. Confidence level; overrides
#'   \code{alpha} when supplied. Default: \code{NULL}.
#' @param strata.name Character or \code{NULL}. Optional stratification
#'   variable entered as a factor covariate. Default: \code{NULL}.
#' @param show_plot Logical. Display a base-R diagnostic plot of the
#'   estimated treatment effect profile. Default: \code{FALSE}.
#' @param verbose Logical. Print diagnostic messages. Default: \code{FALSE}.
#'
#' @return An object of class \code{"glm_effect_profile"}: a named list
#'   containing \code{z_profile} (biomarker grid), \code{est} (point
#'   estimates on the link scale), \code{lower} and \code{upper}
#'   (confidence bounds), \code{se} (standard errors),
#'   \code{effect_measure}, \code{family_used}, \code{overdispersion},
#'   \code{dispersion} (quasi or negbin dispersion estimate, else NA),
#'   \code{model_fit}, \code{spline_basis}, and \code{alpha}.
#'
#' @seealso \code{\link{glm_cs_fit}} for a simpler family-based interface,
#'   \code{\link{cox_cs_fit}} for survival outcomes,
#'   \code{\link{grf.subg.harm.glm}} for GRF-based subgroup identification.
#'
#' @importFrom splines ns
#' @importFrom stats coef vcov qnorm quantile glm binomial poisson gaussian quasipoisson model.matrix as.formula update complete.cases
#' @importFrom graphics plot lines abline
#' @export
glm_effect_profile <- function(
    df,
    outcome.name,
    treat.name      = "treat",
    z.name          = "bm",
    effect_measure  = c("log_OR", "log_RR", "log_IRR", "RD", "MD"),
    offset.name     = NULL,
    overdispersion  = c("none", "quasi", "negbin"),
    spline_df       = 3L,
    z_by            = 0.05,
    z_quantile      = 0.95,
    alpha           = 0.05,
    conf.level      = NULL,
    strata.name     = NULL,
    show_plot       = FALSE,
    verbose         = FALSE
) {

  # ---------------------------------------------------------------------------
  # Argument matching
  # ---------------------------------------------------------------------------
  effect_measure <- match.arg(effect_measure)
  overdispersion <- match.arg(overdispersion)

  # conf.level overrides alpha
  if (!is.null(conf.level)) {
    alpha <- 1 - conf.level
  }

  z_mult <- stats::qnorm(1 - alpha / 2)

  # ---------------------------------------------------------------------------
  # Input validation
  # ---------------------------------------------------------------------------
  required_cols <- c(outcome.name, treat.name, z.name)
  if (!is.null(offset.name))  required_cols <- c(required_cols, offset.name)
  if (!is.null(strata.name))  required_cols <- c(required_cols, strata.name)

  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0L) {
    stop("Columns not found in df: ", paste(missing_cols, collapse = ", "))
  }

  if (effect_measure == "log_IRR" && is.null(offset.name)) {
    stop(
      "'offset.name' must be supplied when effect_measure = 'log_IRR'. ",
      "Provide the name of the person-time (exposure) column."
    )
  }

  if (!is.null(offset.name)) {
    offset_vals <- df[[offset.name]]
    if (any(is.na(offset_vals) | offset_vals <= 0, na.rm = TRUE)) {
      stop("'", offset.name, "' contains non-positive or missing values. ",
           "Offset must be strictly positive (person-time / exposure).")
    }
  }

  if (overdispersion != "none" && !effect_measure %in% c("log_IRR", "log_RR")) {
    warning(
      "overdispersion = '", overdispersion, "' is only applicable to ",
      "Poisson-family models (log_IRR, log_RR). Ignoring for effect_measure = '",
      effect_measure, "'."
    )
    overdispersion <- "none"
  }

  # ---------------------------------------------------------------------------
  # Resolve GLM family
  # ---------------------------------------------------------------------------
  family_obj <- .gep_resolve_family(effect_measure, overdispersion)
  family_used <- family_obj$label  # human-readable label

  if (verbose) {
    message(sprintf(
      "[glm_effect_profile] outcome = '%s' | effect_measure = '%s' | family = '%s' | ",
      outcome.name, effect_measure, family_used
    ), appendLF = FALSE)
    message(sprintf(
      "overdispersion = '%s' | offset = '%s' | n = %d",
      overdispersion, if (is.null(offset.name)) "none" else offset.name, nrow(df)
    ))
  }

  # ---------------------------------------------------------------------------
  # Prepare vectors
  # ---------------------------------------------------------------------------
  Y     <- df[[outcome.name]]
  Treat <- df[[treat.name]]
  Z     <- df[[z.name]]

  # Remove incomplete cases on key columns
  keep <- complete.cases(Y, Treat, Z)
  if (!is.null(offset.name)) keep <- keep & !is.na(df[[offset.name]])
  if (!is.null(strata.name)) keep <- keep & !is.na(df[[strata.name]])
  if (sum(keep) < 10L) stop("Fewer than 10 complete observations.")
  Y     <- Y[keep]
  Treat <- Treat[keep]
  Z     <- Z[keep]

  log_offset <- if (!is.null(offset.name)) log(df[[offset.name]][keep]) else NULL

  # ---------------------------------------------------------------------------
  # Biomarker profile grid
  # ---------------------------------------------------------------------------
  z_lower <- min(Z, na.rm = TRUE)
  z_upper <- stats::quantile(Z, probs = z_quantile, na.rm = TRUE)
  z_profile <- seq(z_lower, z_upper, by = z_by)

  if (verbose) {
    message(sprintf(
      "[glm_effect_profile] Biomarker range: [%.3f, %.3f] | profile points: %d",
      z_lower, z_upper, length(z_profile)
    ))
  }

  # ---------------------------------------------------------------------------
  # Natural spline basis
  # ---------------------------------------------------------------------------
  z_basis   <- splines::ns(Z, df = spline_df)
  z_matrix  <- as.matrix(z_basis)
  z_treat   <- Treat * z_matrix

  col_treat  <- treat.name
  col_z      <- paste0("z_spline_", seq_len(spline_df))
  col_zt     <- paste0("treat_z_spline_", seq_len(spline_df))

  # Build data frame for glm()  -- keeps formula interface simple
  df_fit <- data.frame(
    .Y     = Y,
    .Treat = Treat,
    z_matrix,
    z_treat
  )
  names(df_fit)[3L:(2L + spline_df)]                             <- col_z
  names(df_fit)[(3L + spline_df):(2L + 2L * spline_df)]         <- col_zt
  if (!is.null(log_offset)) df_fit[[".log_offset"]] <- log_offset
  if (!is.null(strata.name)) df_fit[[".strata"]] <- df[[strata.name]][keep]

  rhs <- paste(c(".Treat", col_z, col_zt), collapse = " + ")
  if (!is.null(strata.name)) rhs <- paste0(rhs, " + factor(.strata)")
  fml <- stats::as.formula(paste(".Y ~", rhs))

  # ---------------------------------------------------------------------------
  # Fit GLM
  # ---------------------------------------------------------------------------
  fit <- .gep_fit_model(
    fml           = fml,
    df_fit        = df_fit,
    log_offset    = log_offset,
    family_obj    = family_obj,
    overdispersion = overdispersion,
    verbose       = verbose
  )

  beta_hat <- stats::coef(fit$model)
  vcov_mat <- fit$vcov  # may be sandwich-corrected

  # ---------------------------------------------------------------------------
  # Delta-method profile
  # ---------------------------------------------------------------------------
  z_basis_profile <- predict(z_basis, newx = z_profile)
  z_mat_prof      <- as.matrix(z_basis_profile)

  n_params <- length(beta_hat)
  n_z      <- length(z_profile)

  idx_treat <- which(names(beta_hat) == ".Treat")
  idx_zt    <- which(names(beta_hat) %in% col_zt)

  est   <- numeric(n_z)
  se    <- numeric(n_z)
  lower <- numeric(n_z)
  upper <- numeric(n_z)

  for (j in seq_len(n_z)) {
    zb <- z_mat_prof[j, ]  # length spline_df

    # Contrast: treat effect at biomarker value z_profile[j]
    cvec <- numeric(n_params)
    cvec[idx_treat] <- 1
    cvec[idx_zt]    <- zb

    est_j  <- sum(cvec * beta_hat)
    var_j  <- as.numeric(t(cvec) %*% vcov_mat %*% cvec)
    se_j   <- sqrt(max(var_j, 0))

    est[j]   <- est_j
    se[j]    <- se_j
    lower[j] <- est_j - z_mult * se_j
    upper[j] <- est_j + z_mult * se_j
  }

  # ---------------------------------------------------------------------------
  # Optional plot
  # ---------------------------------------------------------------------------
  if (show_plot) {
    ylab <- switch(effect_measure,
                   log_OR  = "log(OR)",
                   log_RR  = "log(RR)",
                   log_IRR = "log(IRR)",
                   RD      = "Risk Difference",
                   MD      = "Mean Difference",
                   "Effect estimate")
    graphics::plot(
      z_profile, est, type = "l", lwd = 2,
      xlab = z.name, ylab = ylab,
      ylim = range(c(lower, upper), na.rm = TRUE),
      main = paste("Treatment Effect Profile:", effect_measure)
    )
    graphics::lines(z_profile, lower, lty = 2)
    graphics::lines(z_profile, upper, lty = 2)
    graphics::abline(h = 0, lty = 3, col = "grey50")
  }

  # ---------------------------------------------------------------------------
  # Return (with class for S3 dispatch)
  # ---------------------------------------------------------------------------
  result <- list(
    z_profile      = z_profile,
    est            = est,
    lower          = lower,
    upper          = upper,
    se             = se,
    effect_measure = effect_measure,
    family_used    = family_used,
    overdispersion = overdispersion,
    dispersion     = fit$dispersion,
    model_fit      = fit$model,
    spline_basis   = z_basis,
    alpha          = alpha
  )
  class(result) <- "glm_effect_profile"
  result
}


# =============================================================================
# S3 print method
# =============================================================================

#' Print method for glm_effect_profile objects
#'
#' @param x A \code{glm_effect_profile} object.
#' @param ... Additional arguments (unused).
#' @return Invisibly returns \code{x}.
#' @export
print.glm_effect_profile <- function(x, ...) {
  cat("GLM Effect Profile (", x$family_used, " / ",
      x$effect_measure, ")\n", sep = "")
  cat("  Profile points:", length(x$z_profile), "\n")
  if (length(x$est) > 0L) {
    cat("  Effect range:", round(min(x$est, na.rm = TRUE), 3), "to",
        round(max(x$est, na.rm = TRUE), 3), "\n")
  }
  if (!is.na(x$dispersion)) {
    cat("  Dispersion:", round(x$dispersion, 4), "\n")
  }
  cat("  CI level:", (1 - x$alpha) * 100, "%\n")
  invisible(x)
}


# =============================================================================
# Internal helpers (prefixed .gep_ to avoid collision with glm_cs_fit helpers)
# =============================================================================

# Resolve GLM family object from effect_measure and overdispersion
#' @noRd
.gep_resolve_family <- function(effect_measure, overdispersion) {

  switch(effect_measure,

    log_OR = list(
      family = stats::binomial(link = "logit"),
      label  = "binomial(logit)",
      is_poisson_type = FALSE
    ),

    RD = list(
      family = stats::binomial(link = "identity"),
      label  = "binomial(identity)",
      is_poisson_type = FALSE
    ),

    MD = list(
      family = stats::gaussian(link = "identity"),
      label  = "gaussian(identity)",
      is_poisson_type = FALSE
    ),

    log_RR = {
      # Modified Poisson (Zou 2004) for binary outcome
      # Sandwich SEs handle variance misspecification
      list(
        family          = stats::poisson(link = "log"),
        label           = "poisson(log) [modified Poisson for RR]",
        is_poisson_type = TRUE,
        use_sandwich    = (overdispersion == "none")
      )
    },

    log_IRR = {
      # True count / rate outcome
      if (overdispersion == "quasi") {
        list(
          family          = stats::quasipoisson(link = "log"),
          label           = "quasipoisson(log)",
          is_poisson_type = TRUE,
          use_sandwich    = FALSE
        )
      } else if (overdispersion == "negbin") {
        list(
          family          = NULL,  # MASS::glm.nb() used directly
          label           = "negative_binomial",
          is_poisson_type = TRUE,
          use_sandwich    = FALSE
        )
      } else {
        list(
          family          = stats::poisson(link = "log"),
          label           = "poisson(log)",
          is_poisson_type = TRUE,
          use_sandwich    = FALSE  # model-based for proper Poisson
        )
      }
    },

    stop("Unsupported effect_measure: ", effect_measure)
  )
}


# Fit the GLM model with appropriate family and variance estimation
#' @noRd
.gep_fit_model <- function(fml, df_fit, log_offset, family_obj,
                            overdispersion, verbose) {

  is_negbin <- identical(family_obj$label, "negative_binomial")

  # ---- Negative binomial (MASS::glm.nb) ------------------------------------
  if (is_negbin) {
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop(
        "Package 'MASS' is required for overdispersion = 'negbin'. ",
        "Install it with: install.packages('MASS')"
      )
    }
    if (!is.null(log_offset)) {
      fml_nb <- stats::update(fml, paste(". ~ . + offset(.log_offset)"))
    } else {
      fml_nb <- fml
    }
    model <- tryCatch(
      MASS::glm.nb(fml_nb, data = df_fit),
      error = function(e) {
        warning("MASS::glm.nb() failed: ", e$message,
                ". Falling back to poisson().")
        stats::glm(fml, data = df_fit, family = stats::poisson(link = "log"),
                   offset = log_offset)
      }
    )
    dispersion <- if (inherits(model, "negbin")) 1 / model$theta else NA_real_

    return(list(
      model      = model,
      vcov       = stats::vcov(model),
      dispersion = dispersion
    ))
  }

  # ---- Standard GLM --------------------------------------------------------
  if (!is.null(log_offset)) {
    fml_off <- stats::update(fml, ". ~ . + offset(.log_offset)")
  } else {
    fml_off <- fml
  }

  model <- stats::glm(fml_off, data = df_fit, family = family_obj$family)

  # Variance: sandwich or model-based
  use_sandwich <- isTRUE(family_obj$use_sandwich)

  if (use_sandwich) {
    if (!requireNamespace("sandwich", quietly = TRUE)) {
      warning(
        "Package 'sandwich' is recommended for effect_measure = 'log_RR' ",
        "(modified Poisson). Falling back to model-based SEs."
      )
      vcov_mat <- stats::vcov(model)
    } else {
      vcov_mat <- sandwich::vcovHC(model, type = "HC0")
    }
  } else {
    vcov_mat <- stats::vcov(model)
  }

  # Dispersion (quasi only; else NA)
  dispersion <- if (grepl("quasi", family_obj$label, fixed = TRUE)) {
    summary(model)$dispersion
  } else {
    NA_real_
  }

  if (verbose && !is.na(dispersion)) {
    message(sprintf("[glm_effect_profile] Estimated dispersion: %.4f",
                    dispersion))
  }

  list(
    model      = model,
    vcov       = vcov_mat,
    dispersion = dispersion
  )
}
