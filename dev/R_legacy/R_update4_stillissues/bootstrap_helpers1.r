#' Target Estimate and Standard Error for Bootstrap (Fixed)
#'
#' Calculates target estimate and standard error for bootstrap samples with proper NA handling.
#'
#' @param x Numeric vector of estimates (may contain NA).
#' @param ystar Matrix of bootstrap samples.
#' @param cov_method Character. Covariance method ("standard" or "nocorrect").
#' @param cov_trim Numeric. Trimming proportion for covariance (default: 0.0).
#' @return List with target estimate, standard errors, and correction term.
#' @export

get_targetEst <- function(x, ystar, cov_method = "standard", cov_trim = 0.0) {
  # Input validation
  if (!is.numeric(x)) stop("'x' must be numeric")
  if (!is.matrix(ystar) && !is.array(ystar)) stop("'ystar' must be a matrix")
  if (!is.character(cov_method) || length(cov_method) != 1) {
    stop("'cov_method' must be a single character string")
  }
  if (!is.numeric(cov_trim) || length(cov_trim) != 1 || cov_trim < 0 || cov_trim > 1) {
    stop("'cov_trim' must be a single numeric value between 0 and 1")
  }

  # CRITICAL FIX: Return NA results for all-NA or insufficient data
  n_valid <- sum(!is.na(x))

  # Return NA structure if all values are NA
  if (n_valid == 0) {
    return(list(
      target_est = NA_real_,
      sehat = NA_real_,
      sehat_new = NA_real_,
      term_correct = NA_real_,
      varhat = NA_real_,
      n_valid = 0,
      n_total = length(x)
    ))
  }

  # Return NA structure if insufficient non-NA values (need at least 2 for variance)
  if (n_valid < 2) {
    # Can calculate mean with 1 value, but not variance
    mx <- mean(x, na.rm = TRUE)
    return(list(
      target_est = mx,
      sehat = NA_real_,
      sehat_new = NA_real_,
      term_correct = NA_real_,
      varhat = NA_real_,
      n_valid = n_valid,
      n_total = length(x)
    ))
  }

  # Calculate mean excluding NAs
  mx <- mean(x, na.rm = TRUE)

  # Center the estimates (excluding NAs)
  xc <- x - mx

  # Get dimensions
  N <- ncol(ystar)
  B.eval <- sum(!is.na(xc))

  # Calculate covariance for each observation
  if (cov_method == "standard" || cov_method == "nocorrect") {
    # Modified calc_cov to handle NAs
    calc_cov_safe <- function(ystar_col, Est) {
      # Remove NAs from both ystar_col and Est
      valid_idx <- !is.na(Est)
      if (sum(valid_idx) < 2) return(0)

      ystar_valid <- ystar_col[valid_idx]
      Est_valid <- Est[valid_idx]

      # Calculate covariance
      mean((ystar_valid - mean(ystar_valid)) * Est_valid, na.rm = TRUE)
    }

    # Apply to each column of ystar
    cov_i <- apply(ystar, 2, calc_cov_safe, Est = xc)
  }

  # Handle different covariance methods
  if (cov_method != "nocorrect") {
    # Standard method with correction
    varhat <- N * mean(cov_i^2, trim = cov_trim, na.rm = TRUE)
    seH <- sqrt(varhat)

    # Calculate correction term
    nb_ratio <- N / (B.eval^2)
    termc <- nb_ratio * B.eval * mean(xc^2, na.rm = TRUE, trim = cov_trim)

    # Apply correction
    if (varhat <= termc) {
      varhat_new <- varhat
      seH_new <- sqrt(varhat)
    } else {
      varhat_new <- varhat - termc
      seH_new <- sqrt(varhat_new)
    }
  } else {
    # No correction method
    varhat <- N * mean(cov_i^2, trim = cov_trim, na.rm = TRUE)
    seH <- sqrt(varhat)
    termc <- 0.0
    varhat_new <- varhat
    seH_new <- seH
  }

  # Return results
  list(
    target_est = mx,
    sehat = seH,
    sehat_new = seH_new,
    term_correct = termc,
    varhat = varhat_new,
    n_valid = n_valid,
    n_total = length(x)
  )
}


#' Calculate Covariance for Bootstrap Estimates (Fixed)
#'
#' Calculates the covariance between a vector and bootstrap estimates with NA handling.
#'
#' @param x Numeric vector.
#' @param Est Numeric vector of bootstrap estimates (may contain NA).
#' @return Numeric value of covariance.
#' @export
calc_cov <- function(x, Est) {
  # CRITICAL FIX: Add input validation
  if (!is.numeric(x)) stop("'x' must be numeric")
  if (!is.numeric(Est)) stop("'Est' must be numeric")
  if (length(x) != length(Est)) stop("'x' and 'Est' must have the same length")

  # Remove NA values from Est and corresponding x values
  valid_idx <- !is.na(Est) & !is.na(x)

  # Need at least 2 valid values
  if (sum(valid_idx) < 2) {
    return(0)
  }

  x_valid <- x[valid_idx]
  Est_valid <- Est[valid_idx]

  # CRITICAL FIX: Check for zero variance
  if (var(x_valid) == 0 || var(Est_valid) == 0) {
    return(0)
  }

  # Calculate covariance with additional safety
  cov_result <- tryCatch(
    mean((x_valid - mean(x_valid, na.rm = TRUE)) * Est_valid, na.rm = TRUE),
    error = function(e) 0,
    warning = function(w) 0
  )

  # Ensure finite result
  if (!is.finite(cov_result)) {
    return(0)
  }

  return(cov_result)
}


#' Bootstrap Confidence Interval and Bias Correction Results (Fixed)
#'
#' Calculates confidence intervals and bias-corrected estimates for bootstrap results with NA handling.
#'
#' @param Hobs Numeric. Observed estimate.
#' @param seHobs Numeric. Standard error of observed estimate.
#' @param H1_adj Numeric vector. Bias-corrected estimate 1 (may contain NA).
#' @param H2_adj Numeric vector. Bias-corrected estimate 2 (may contain NA).
#' @param ystar Matrix of bootstrap samples.
#' @param cov_method Character. Covariance method ("standard" or "nocorrect").
#' @param cov_trim Numeric. Trimming proportion for covariance (default: 0.0).
#' @param est.scale Character. "hr" or "1/hr".
#' @param est.loghr Logical. Is estimate on log(HR) scale?
#' @return Data.table with confidence intervals and estimates.
#' @importFrom data.table data.table
#' @export
get_dfRes <- function(Hobs, seHobs, H1_adj, H2_adj = NULL, ystar,
                      cov_method = "standard", cov_trim = 0.0,
                      est.scale = "hr", est.loghr = TRUE) {

  # Load required package
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required.")
  }

  # Input validation
  if (!is.numeric(Hobs) || length(Hobs) != 1 || is.na(Hobs) || is.infinite(Hobs)) {
    stop("'Hobs' must be a single, finite numeric value")
  }
  if (!is.numeric(seHobs) || length(seHobs) != 1 || is.na(seHobs) || is.infinite(seHobs) || seHobs < 0) {
    stop("'seHobs' must be a single, non-negative, finite numeric value")
  }
  if (missing(ystar)) {
    stop("'ystar' must be provided")
  }
  if (!is.character(cov_method) || length(cov_method) != 1) {
    stop("'cov_method' must be a single character string")
  }
  if (!is.numeric(cov_trim) || length(cov_trim) != 1 || is.na(cov_trim) || cov_trim < 0 || cov_trim > 1) {
    stop("'cov_trim' must be a single numeric value between 0 and 1")
  }
  if (!est.scale %in% c("hr", "1/hr")) {
    stop("'est.scale' must be either 'hr' or '1/hr'")
  }
  if (!is.logical(est.loghr) || length(est.loghr) != 1 || is.na(est.loghr)) {
    stop("'est.loghr' must be a single logical value")
  }

  # CRITICAL FIX: Handle all-NA cases gracefully
  # Un-adjusted estimates (these don't depend on bootstrap)
  cest <- ci_est(x = Hobs, sd = seHobs, scale = est.scale, est.loghr = est.loghr)
  H0_lower <- cest$lower
  H0_upper <- cest$upper
  sdH0 <- cest$sd
  H0 <- cest$est

  # Initialize H1 results with NA values
  H1 <- NA_real_
  sdH1 <- NA_real_
  H1_lower <- NA_real_
  H1_upper <- NA_real_

  # Calculate H1 adjusted estimates only if we have valid data
  if (!all(is.na(H1_adj))) {
    est <- tryCatch(
      get_targetEst(x = H1_adj, ystar = ystar, cov_method = cov_method, cov_trim = cov_trim),
      error = function(e) {
        warning("Failed to calculate H1 estimates: ", e$message)
        list(target_est = NA, sehat_new = NA)
      }
    )

    # Only calculate CI if we have valid estimates
    if (!is.na(est$target_est) && !is.na(est$sehat_new)) {
      q <- est$target_est
      se_new <- est$sehat_new

      cest <- ci_est(x = q, sd = se_new, scale = est.scale, est.loghr = est.loghr)
      H1_lower <- cest$lower
      H1_upper <- cest$upper
      sdH1 <- cest$sd
      H1 <- cest$est
    }
  }

  # Create basic output
  outres <- data.table::data.table(
    H0, sdH0, H0_lower, H0_upper,
    H1, sdH1, H1_lower, H1_upper
  )

  # Handle H2_adj if provided
  if (!is.null(H2_adj)) {
    # Initialize H2 results with NA values
    H2 <- NA_real_
    sdH2 <- NA_real_
    H2_lower <- NA_real_
    H2_upper <- NA_real_

    if (!all(is.na(H2_adj))) {
      # Calculate H2 adjusted estimates
      est <- tryCatch(
        get_targetEst(x = H2_adj, ystar = ystar, cov_method = cov_method, cov_trim = cov_trim),
        error = function(e) {
          warning("Failed to calculate H2 estimates: ", e$message)
          list(target_est = NA, sehat_new = NA)
        }
      )

      if (!is.na(est$target_est) && !is.na(est$sehat_new)) {
        q <- est$target_est
        se_new <- est$sehat_new
        cest <- ci_est(x = q, sd = se_new, scale = est.scale, est.loghr = est.loghr)
        H2_lower <- cest$lower
        H2_upper <- cest$upper
        sdH2 <- cest$sd
        H2 <- cest$est
      }
    }

    outres <- data.table::data.table(
      H0, sdH0, H0_lower, H0_upper,
      H1, sdH1, H1_lower, H1_upper,
      H2, sdH2, H2_lower, H2_upper
    )
  }

  return(outres)
}

#' Find integer pairs (x, y) such that x * y = z and y >= x
#'
#' Given an integer z, this function finds all integer pairs (x, y) such that x * y = z and y >= x.
#' Optionally, you can return only the pair with the largest value of x or y.
#'
#' @param z Integer. The target product.
#' @param return_largest Character. If \"x\", returns the pair with the largest x. If \"y\", returns the pair with the largest y. If NULL, returns all pairs.
#'
#' @return A matrix of integer pairs (x, y) satisfying the conditions, or a single pair if return_largest is specified.
#' @examples
#' find_xy_given_z(12)
#' find_xy_given_z(12, return_largest = \"x\")
#' find_xy_given_z(12, return_largest = \"y\")
find_xy_given_z <- function(z, return_largest = NULL) {
  pairs <- list()
  for (x in 1:z) {
    if (z %% x == 0) {
      y <- z / x
      if (y >= x && y %% 1 == 0) {
        pairs[[length(pairs) + 1]] <- c(x, y)
      }
    }
  }
  result <- do.call(rbind, pairs)
  if (!is.null(return_largest)) {
    if (return_largest == "x") {
      idx <- which.max(result[,1])
      return(result[idx, , drop = FALSE])
    } else if (return_largest == "y") {
      idx <- which.max(result[,2])
      return(result[idx, , drop = FALSE])
    }
  }
  return(result)
}


#' Ensure Required Packages Are Installed and Loaded
#'
#' Installs and loads required packages if not already available.
#'
#' @param pkgs Character vector of package names.
#' @return None. Packages are loaded into the session.
#' @export
ensure_packages <- function(pkgs) {
  for (pkg in pkgs) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

#' Build Cox Model Formula
#'
#' Constructs a Cox model formula from variable names.
#'
#' @param outcome.name Character. Name of outcome variable.
#' @param event.name Character. Name of event indicator variable.
#' @param treat.name Character. Name of treatment variable.
#' @return An R formula object for Cox regression.
#' @export
build_cox_formula <- function(outcome.name, event.name, treat.name) {
  sf <- paste0("Surv(", outcome.name, ",", event.name, ") ~ ", treat.name)
  as.formula(sf)
}

#' Fit Cox Models for Subgroups
#'
#' Fits Cox models for two subgroups defined by treatment recommendation.
#'
#' @param df Data frame.
#' @param formula Cox model formula.
#' @return List with HR and SE for each subgroup.
#' @export
fit_cox_models <- function(df, formula) {
  fitH <- get_Cox_sg(df_sg = subset(df, treat.recommend == 0), cox.formula = formula)
  fitHc <- get_Cox_sg(df_sg = subset(df, treat.recommend == 1), cox.formula = formula)
  list(H_obs = fitH$est_obs, seH_obs = fitH$se_obs, Hc_obs = fitHc$est_obs, seHc_obs = fitHc$se_obs)
}


#' Bootstrap Ystar Matrix
#'
#' Generates a bootstrap matrix for Ystar using parallel processing.
#'
#' @param df Data frame.
#' @param nb_boots Integer. Number of bootstrap samples.
#' @return Matrix of bootstrap samples.
#' @importFrom foreach foreach
#' @export
bootstrap_ystar <- function(df, nb_boots) {
  NN <- nrow(df)
  foreach::foreach(
    boot = seq_len(nb_boots),
    .options.future = list(seed = TRUE),
    .combine = "rbind",
    .errorhandling = "pass"
  ) %dofuture% {
    set.seed(8316951 + boot * 100)
    in_boot <- sample.int(NN, size = NN, replace = TRUE)
    df_boot <- df[in_boot, ]
    df_boot$id_boot <- seq_len(nrow(df_boot))
    ystar <- unlist(lapply(df$id, count.id, dfb = df_boot))
    return(ystar)
  }
}

#' Format Confidence Interval for Estimates
#'
#' Formats confidence interval for estimates.
#'
#' @param estimates Data frame or data.table of estimates.
#' @param col_names Character vector of column names for estimate, lower, upper.
#' @return Character string formatted as \"estimate (lower, upper)\".
#' @export
format_CI <- function(estimates, col_names) {
  resH <- estimates[, ..col_names]
  Hstat <- round(unlist(resH[1, ]), 2)
  paste0(Hstat[1], " (", Hstat[2], ",", Hstat[3], ")")
}

