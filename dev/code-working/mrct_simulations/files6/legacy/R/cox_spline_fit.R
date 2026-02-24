#' Fit Cox Model with Cubic Spline for Treatment Effect Heterogeneity
#'
#' Estimates treatment effects as a function of a continuous covariate using
#' a Cox proportional hazards model with natural cubic splines. The function
#' models treatment-by-covariate interactions to detect effect modification.
#'
#' @param df Data frame containing survival data
#' @param tte_name Character string specifying time-to-event variable name.
#'   Default: "os_time"
#' @param event_name Character string specifying event indicator variable name
#'   (1=event, 0=censored). Default: "os_event"
#' @param treat_name Character string specifying treatment variable name
#'   (1=treated, 0=control). Default: "treat"
#' @param strata_name Character string specifying stratification variable name.
#'   If NULL, no stratification is used. Default: NULL
#' @param z_name Character string specifying continuous covariate name for
#'   effect modification. Default: "bm"
#' @param alpha Numeric value for confidence level (two-sided). Default: 0.20
#'   (80% confidence intervals)
#' @param spline_df Integer specifying degrees of freedom for natural spline.
#'   Default: 3
#' @param z_max Numeric maximum value for z in predictions. Values beyond this
#'   are truncated. Default: Inf (no truncation)
#' @param z_by Numeric increment for z values in prediction grid. Default: 1
#' @param z_window Numeric half-width for counting observations near each z value.
#'   Default: 0.0 (exact matches only)
#' @param z_quantile Numeric quantile (0-1) for upper limit of z profile.
#'   Default: 0.90 (90th percentile)
#' @param show_plot Logical indicating whether to display plot. Default: TRUE
#' @param plot_params List of plotting parameters (see Details). Default: NULL
#' @param truebeta_name Character string specifying variable containing true
#'   log(HR) values for validation/simulation. Default: NULL
#' @param verbose Logical indicating whether to print diagnostic information.
#'   Default: TRUE
#'
#' @details
#' ## Model Structure
#'
#' The function fits:
#' \deqn{h(t|Z,A) = h_0(t) \exp(\beta_0 A + f(Z) + g(Z) \cdot A)}
#'
#' Where:
#' - A is treatment (0/1)
#' - Z is the continuous effect modifier
#' - f(Z) is modeled with natural splines (main effect)
#' - g(Z) is modeled with natural splines (interaction)
#' - The log hazard ratio is: \eqn{\beta(Z) = \beta_0 + g(Z)}
#'
#' ## Plot Parameters
#'
#' The `plot_params` argument accepts a list with:
#' - `xlab`: x-axis label
#' - `main_title`: plot title
#' - `ylimit`: y-axis limits c(min, max)
#' - `y_pad_zero`: padding below zero line
#' - `y_delta`: extra space for count labels
#' - `cex_legend`: legend text size
#' - `cex_count`: count text size
#' - `show_cox_primary`: show standard Cox estimate line
#' - `show_null`: show null effect line (log(HR)=0)
#' - `show_target`: show target effect line (e.g., log(0.80))
#'
#' @return List containing:
#' \describe{
#'   \item{z_profile}{Vector of z values where treatment effect is estimated}
#'   \item{loghr_est}{Point estimates of log(HR) at each z value}
#'   \item{loghr_lower}{Lower confidence bound}
#'   \item{loghr_upper}{Upper confidence bound}
#'   \item{se_loghr}{Standard errors of log(HR) estimates}
#'   \item{counts_profile}{Number of observations near each z value}
#'   \item{cox_primary}{Log(HR) from standard Cox model (no interaction)}
#'   \item{model_fit}{The fitted coxph model object}
#'   \item{spline_basis}{The natural spline basis object}
#' }
#'
#' @examples
#' \dontrun{
#' # Simulate data
#' set.seed(123)
#' df <- data.frame(
#'   os_time = rexp(500, 0.01),
#'   os_event = rbinom(500, 1, 0.7),
#'   treat = rbinom(500, 1, 0.5),
#'   bm = rnorm(500, 50, 10)
#' )
#'
#' # Fit model
#' result <- cox_cs_fit(df, z_name = "bm", alpha = 0.20)
#'
#' # Custom plotting
#' result <- cox_cs_fit(
#'   df,
#'   z_name = "bm",
#'   plot_params = list(
#'     xlab = "Biomarker Level",
#'     main_title = "Treatment Effect by Biomarker",
#'     cex_legend = 1.2
#'   )
#' )
#' }
#'
#' @export
#' @importFrom survival coxph Surv
#' @importFrom splines ns
#' @importFrom stats coef quantile qnorm vcov
#' @importFrom graphics plot lines abline axis text legend par box
cox_cs_fit <- function(df,
                       tte_name = "os_time",
                       event_name = "os_event",
                       treat_name = "treat",
                       strata_name = NULL,
                       z_name = "bm",
                       alpha = 0.20,
                       spline_df = 3,
                       z_max = Inf,
                       z_by = 1,
                       z_window = 0.0,
                       z_quantile = 0.90,
                       show_plot = TRUE,
                       plot_params = NULL,
                       truebeta_name = NULL,
                       verbose = TRUE) {

  # ==========================================================================
  # Input Validation
  # ==========================================================================

  required_vars <- c(tte_name, event_name, treat_name, z_name)
  missing_vars <- setdiff(required_vars, names(df))
  if (length(missing_vars) > 0) {
    stop("Missing required variables: ", paste(missing_vars, collapse = ", "))
  }

  if (!is.null(strata_name) && !strata_name %in% names(df)) {
    stop("Stratification variable '", strata_name, "' not found in data")
  }

  if (!is.null(truebeta_name) && !truebeta_name %in% names(df)) {
    stop("True beta variable '", truebeta_name, "' not found in data")
  }

  if (alpha <= 0 || alpha >= 1) {
    stop("alpha must be between 0 and 1")
  }

  if (spline_df < 1) {
    stop("spline_df must be at least 1")
  }

  # ==========================================================================
  # Extract Variables
  # ==========================================================================

  # Sort by biomarker
  df <- df[order(df[[z_name]]),]

  Y <- df[[tte_name]]
  E <- df[[event_name]]
  Treat <- df[[treat_name]]
  z <- df[[z_name]]

  # Check for valid values
  if (any(Y <= 0, na.rm = TRUE)) {
    warning("Non-positive survival times detected")
  }

  if (!all(E %in% c(0, 1))) {
    warning("Event indicator should be 0 or 1")
  }

  if (!all(Treat %in% c(0, 1))) {
    warning("Treatment indicator should be 0 or 1")
  }

  # Handle stratification
  if (!is.null(strata_name)) {
    stratum <- df[[strata_name]]
    # Find most common stratum for prediction reference
    stratum_table <- table(stratum)
    stratum_reference <- names(stratum_table)[which.max(stratum_table)]
  } else {
    stratum <- rep("1", length(Y))
    stratum_reference <- "1"
  }

  # Extract true beta if provided (for validation/simulation studies)
  beta_true <- NULL
  if (!is.null(truebeta_name)) {
    beta_true <- df[[truebeta_name]]
  }

  # ==========================================================================
  # Define Prediction Grid
  # ==========================================================================

  z_quantiles <- quantile(z, probs = c(0.75, 0.80, 0.90), na.rm = TRUE)

  if (verbose) {
    cat("\n=== Z Variable Summary ===\n")
    cat("Range:", round(range(z, na.rm = TRUE), 2), "\n")
    cat("Quantiles (75%, 80%, 90%):", round(z_quantiles, 2), "\n")
  }

  # Create prediction grid
  z_min <- min(z, na.rm = TRUE)
  z_upper <- quantile(z, probs = z_quantile, na.rm = TRUE)
  z_profile <- seq(z_min, z_upper, by = z_by)

  # Apply maximum cutoff
  z_profile <- unique(pmin(z_profile, z_max))

  if (verbose) {
    cat("Prediction grid:", length(z_profile), "points from",
        round(min(z_profile), 2), "to", round(max(z_profile), 2), "\n")
  }

  # ==========================================================================
  # Fit Standard Cox Model (No Interaction)
  # ==========================================================================

  if (!is.null(strata_name)) {
    cox_formula_primary <- as.formula(paste0(
      "Surv(Y, E) ~ Treat + strata(stratum)"
    ))
  } else {
    cox_formula_primary <- as.formula("Surv(Y, E) ~ Treat")
  }

  fit_primary <- coxph(cox_formula_primary)
  loghr_primary <- coef(fit_primary)["Treat"]

  if (verbose) {
    cat("\n=== Primary Cox Model ===\n")
    cat("Treatment log(HR):", round(loghr_primary, 3), "\n")
    cat("Treatment HR:", round(exp(loghr_primary), 3), "\n")
  }

  # ==========================================================================
  # Create Spline Basis
  # ==========================================================================

  # Create natural spline basis for z
  z_basis <- ns(z, df = spline_df)
  z_matrix <- as.matrix(z_basis)

  # Create treatment-spline interactions
  z_treat_matrix <- Treat * z_matrix

  # Combine all covariates
  model_matrix <- cbind(Treat, z_matrix, z_treat_matrix)

  # Give informative column names
  colnames(model_matrix) <- c(
    "treat",
    paste0("z_spline_", 1:spline_df),
    paste0("treat_z_spline_", 1:spline_df)
  )

  # ==========================================================================
  # Fit Cox Model with Spline Interaction
  # ==========================================================================

  if (!is.null(strata_name)) {
    fit_spline <- coxph(Surv(Y, E) ~ model_matrix + strata(stratum))
  } else {
    fit_spline <- coxph(Surv(Y, E) ~ model_matrix)
  }

  beta_hat <- coef(fit_spline)
  vcov_beta <- stats::vcov(fit_spline)

  if (verbose) {
    cat("\n=== Spline Model ===\n")
    cat("Number of parameters:", length(beta_hat), "\n")
    cat("Treatment main effect:", round(beta_hat[1], 3), "\n")
  }

  # ==========================================================================
  # Calculate Predictions
  # ==========================================================================

  # Predict spline basis at profile z values
  z_basis_profile <- predict(z_basis, newx = z_profile)
  z_matrix_profile <- as.matrix(z_basis_profile)

  # Design matrix for treated (A=1)
  treat_1 <- rep(1, length(z_profile))
  z_treat_matrix_1 <- treat_1 * z_matrix_profile
  model_matrix_1 <- cbind(treat_1, z_matrix_profile, z_treat_matrix_1)
  colnames(model_matrix_1) <- colnames(model_matrix)

  # Design matrix for control (A=0)
  treat_0 <- rep(0, length(z_profile))
  z_treat_matrix_0 <- treat_0 * z_matrix_profile
  model_matrix_0 <- cbind(treat_0, z_matrix_profile, z_treat_matrix_0)
  colnames(model_matrix_0) <- colnames(model_matrix)

  # Linear predictors
  lp_1 <- as.vector(model_matrix_1 %*% beta_hat)
  lp_0 <- as.vector(model_matrix_0 %*% beta_hat)

  # Treatment effect (log hazard ratio) as function of z
  loghr_z <- lp_1 - lp_0

  # ==========================================================================
  # Calculate Standard Errors
  # ==========================================================================

  # The log(HR) depends only on treatment and treatment-spline interactions
  # log(HR)(z) = beta_treat + sum(beta_treat_spline_k * spline_k(z))

  # Index for treatment and interaction terms
  idx_treat <- 1
  idx_interact <- (spline_df + 2):(2 * spline_df + 1)
  idx_diff <- c(idx_treat, idx_interact)

  # Design matrix for difference (log HR)
  # Each row is: [1, spline_1(z_i), spline_2(z_i), ..., spline_df(z_i)]
  diff_matrix <- cbind(
    rep(1, length(z_profile)),
    z_matrix_profile
  )

  # Extract relevant covariance submatrix
  vcov_diff <- vcov_beta[idx_diff, idx_diff]

  # Calculate variance of log(HR) estimates
  var_loghr <- diag(diff_matrix %*% vcov_diff %*% t(diff_matrix))
  se_loghr <- sqrt(var_loghr)

  # Confidence intervals
  z_crit <- qnorm(1 - alpha / 2)
  loghr_lower <- loghr_z - z_crit * se_loghr
  loghr_upper <- loghr_z + z_crit * se_loghr

  # ==========================================================================
  # Count Observations Near Each Profile Point
  # ==========================================================================

  if (z_window > 0) {
    counts_profile <- sapply(z_profile, function(z_val) {
      sum(z >= (z_val - z_window) & z <= (z_val + z_window), na.rm = TRUE)
    })
  } else {
    # Count exact matches or use binning
    counts_profile <- sapply(z_profile, function(z_val) {
      sum(abs(z - z_val) < 1e-10, na.rm = TRUE)
    })
    # If no exact matches, use histogram approach
    if (all(counts_profile == 0)) {
      z_breaks <- c(z_profile - z_by/2, max(z_profile) + z_by/2)
      counts_profile <- as.vector(table(cut(z, breaks = z_breaks)))
    }
  }

  # ==========================================================================
  # Create Plot
  # ==========================================================================

  if (show_plot) {
    # Set default plot parameters
    default_plot_params <- list(
      xlab = z_name,
      main_title = NULL,
      ylimit = NULL,
      y_pad_zero = 0.01,
      y_delta = 0.25,
      cex_legend = 0.9,
      cex_count = 0.7,
      show_cox_primary = TRUE,
      show_null = TRUE,
      show_target = TRUE,
      target_loghr = log(0.80)
    )

    # Override with user-provided parameters
    if (!is.null(plot_params)) {
      for (param_name in names(plot_params)) {
        default_plot_params[[param_name]] <- plot_params[[param_name]]
      }
    }
    plot_params <- default_plot_params

    # Calculate y-axis limits
    ymin <- min(loghr_lower)
    ymax <- max(loghr_upper)

    if (!is.null(beta_true)) {
      ymin <- min(ymin, beta_true, na.rm = TRUE)
      ymax <- max(ymax, beta_true, na.rm = TRUE)
    }

    if (is.null(plot_params$ylimit)) {
      y_range_min <- ymin - plot_params$y_delta
      y_range_max <- ymax + plot_params$y_delta
      plot_params$ylimit <- c(y_range_min, y_range_max)
    }

    # Create y-axis tick marks
    y_start <- floor(log(0.5) * 10) / 10
    y_end <- ceiling(ymax * 10) / 10
    y_ticks <- seq(y_start, y_end, length.out = 10)
    y_ticks <- round(y_ticks, 1)
    y_ticks <- sort(unique(c(y_ticks, 0, round(ymin, 1))))

    # Position for count labels
    y_count_pos <- plot_params$ylimit[1] +
      (round(ymin, 1) - plot_params$y_pad_zero -
         plot_params$ylimit[1]) / 2

    # Create main plot
    plot(z_profile, loghr_z,
         type = "l", lty = 1, col = "black", lwd = 3,
         xlab = plot_params$xlab,
         ylab = "log(Hazard Ratio)",
         ylim = plot_params$ylimit,
         axes = FALSE,
         main = plot_params$main_title)

    # Add confidence bands
    lines(z_profile, loghr_lower, lty = 2, col = "black", lwd = 1)
    lines(z_profile, loghr_upper, lty = 2, col = "black", lwd = 1)

    # Add true values if provided
    if (!is.null(beta_true)) {
      lines(z, beta_true, type = "l", lty = 1, col = "grey", lwd = 3)
    }

    # Add reference lines
    if (plot_params$show_null) {
      abline(h = 0, col = "red", lty = 2, lwd = 0.5)
    }

    if (plot_params$show_cox_primary) {
      abline(h = loghr_primary, col = "blue", lty = 2, lwd = 2)
    }

    # if (plot_params$show_target) {
    #   abline(h = plot_params$target_loghr, col = "red", lty = 1, lwd = 1)
    # }

    # Add axes
    axis(2, at = y_ticks, las = 1)
    axis(1, at = z_profile, labels = z_profile)

    # Add horizontal line for count baseline
    y_count_baseline <- round(ymin, 1) - plot_params$y_pad_zero
    abline(h = y_count_baseline, lty = 1, col = "black")

    # Add box
    graphics::box()

    # Add count labels
    text(z_profile, y_count_pos, counts_profile,
         col = "black", cex = plot_params$cex_count)

    # Add legends
    ci_level <- paste0(round((1 - alpha) * 100), "%")
    legend("topleft", horiz = TRUE,
           legend = c("log(HR)", paste0(ci_level, " CI"),"True (causal)"),
           lty = c(1, 2, 1),
           lwd = c(3, 1, 3),
           col = c("black", "black","grey"),
           bty = "n",
          cex = plot_params$cex_legend)

    legend_items <- character()
    legend_lty <- numeric()
    legend_lwd <- numeric()
    legend_col <- character()

    # if (plot_params$show_target) {
    #   legend_items <- c(legend_items, "log(0.80)")
    #   legend_lty <- c(legend_lty, 1)
    #   legend_lwd <- c(legend_lwd, 1)
    #   legend_col <- c(legend_col, "red")
    # }

    if (plot_params$show_cox_primary) {
      legend_items <- c(legend_items, "Cox primary")
      legend_lty <- c(legend_lty, 2)
      legend_lwd <- c(legend_lwd, 2)
      legend_col <- c(legend_col, "blue")
    }

    if (length(legend_items) > 0) {
      legend("topright",
             legend = legend_items,
             lty = legend_lty,
             lwd = legend_lwd,
             col = legend_col,
             bty = "n",
             cex = plot_params$cex_legend)
    }
  }

  # ==========================================================================
  # Return Results
  # ==========================================================================

  results <- list(
    z_profile = z_profile,
    loghr_est = loghr_z,
    loghr_lower = loghr_lower,
    loghr_upper = loghr_upper,
    se_loghr = se_loghr,
    counts_profile = counts_profile,
    cox_primary = loghr_primary,
    model_fit = fit_spline,
    spline_basis = z_basis,
    alpha = alpha,
    ci_level = 1 - alpha
  )

  class(results) <- c("cox_cs_fit", "list")

  return(results)
}





#' Plot Subgroup Analysis Results
#'
#' Creates diagnostic plots for subgroup treatment effects from df_super object
#'
#' @param df_super A data frame containing subgroup analysis results with columns:
#'   loghr_po (log hazard ratios), and optionally theta_1 and theta_0
#'   (treatment-specific parameters)
#' @param z Character string specifying the column name to use as the subgroup
#'   score (e.g., "z_age", "z_size", "subgroup"). Required.
#' @param hrz_crit Critical hazard ratio threshold for defining optimal subgroup.
#'   Default is 1 (HR=1 on log scale is 0).
#' @param log.hrs Optional vector of reference log hazard ratios to display as
#'   horizontal lines. Default is NULL.
#' @param ahr_empirical Optional empirical average hazard ratio to display.
#'   If NULL, calculated from data. Default is NULL.
#' @param plot_type Character string specifying plot type: "both" (default),
#'   "profile", or "ahr".
#' @param add_rug Logical indicating whether to add rug plot of z values.
#'   Default is TRUE.
#' @param zpoints_by Step size for z-axis grid when calculating AHR curves.
#'   Default is 1.
#' @param ... Additional graphical parameters passed to plot()
#'
#' @return A list containing:
#'   \item{cut.zero}{The minimum z value where loghr_po < hrz_crit}
#'   \item{AHR_opt}{Average hazard ratio for optimal subgroup (z >= cut.zero)}
#'   \item{zpoints}{Grid of z values used for AHR calculations}
#'   \item{HR.zpoints}{AHR for population with z >= zpoints}
#'   \item{HRminus.zpoints}{AHR for population with z <= zpoints}
#'   \item{HR2.zpoints}{Alternative AHR calculation for z >= zpoints}
#'   \item{HRminus2.zpoints}{Alternative AHR calculation for z <= zpoints}
#'
#' @details
#' The function creates up to two plots:
#' 1. Treatment effect profile: Shows log hazard ratio as function of z
#' 2. Average hazard ratio curve: Shows AHR for subgroups z >= threshold
#'
#' The "optimal" subgroup is defined as patients with z >= cut.zero, where
#' cut.zero is the minimum z value with favorable treatment effect (loghr < hrz_crit).
#'
#' @examples
#' \dontrun{
#' # Using z_age as the subgroup score
#' results <- plot_subgroup_effects(dgm_spline$df_super, z = "z_age", hrz_crit = 0)
#'
#' # Using subgroup identifier
#' results <- plot_subgroup_effects(dgm_spline$df_super, z = "subgroup", hrz_crit = 0)
#'
#' # With reference lines
#' results <- plot_subgroup_effects(dgm_spline$df_super, z = "z_size",
#'                                   hrz_crit = 0,
#'                                   log.hrs = c(-0.5, 0, 0.5))
#'
#' # Only AHR plot
#' results <- plot_subgroup_effects(dgm_spline$df_super, z = "z_pgr",
#'                                   plot_type = "ahr")
#' }
#'
#' @importFrom stats quantile
#' @importFrom graphics par plot rug abline legend
#' @export
plot_subgroup_effects <- function(df_super,
                                  z,
                                  hrz_crit = 0,
                                  log.hrs = NULL,
                                  ahr_empirical = NULL,
                                  plot_type = c("both", "profile", "ahr"),
                                  add_rug = TRUE,
                                  zpoints_by = 1,
                                  ...) {

  # Input validation
  if (!is.data.frame(df_super)) {
    stop("df_super must be a data frame")
  }

  if (missing(z)) {
    stop("Argument 'z' is required. Please specify the column name to use as subgroup score.")
  }

  if (!z %in% names(df_super)) {
    stop("Column '", z, "' not found in df_super. Available columns: ",
         paste(names(df_super), collapse = ", "))
  }

  if (!"loghr_po" %in% names(df_super)) {
    stop("df_super missing required column: loghr_po")
  }

  plot_type <- match.arg(plot_type)

  # Sort by biomarker
  df_super <- df_super[order(df_super[[z]]),]


  # Extract the z variable
  z_values <- df_super[[z]]
  loghr_values <- df_super$loghr_po

  # Remove NA values
  valid_idx <- !is.na(z_values) & !is.na(loghr_values)
  if (sum(valid_idx) == 0) {
    stop("No valid observations after removing NAs")
  }

  df_work <- df_super[valid_idx, ]
  z_values <- df_work[[z]]
  loghr_values <- df_work$loghr_po

  # Calculate optimal cutpoint
  valid_indices <- which(loghr_values < hrz_crit)
  if (length(valid_indices) == 0) {
    warning("No observations with loghr_po < hrz_crit. Using minimum z value.")
    cut.zero <- min(z_values)
  } else {
    cut.zero <- min(z_values[valid_indices])
  }

  # Calculate AHR for optimal subgroup
  dfs_opt <- df_work[z_values >= cut.zero, ]
  if (nrow(dfs_opt) == 0) {
    stop("No observations in optimal subgroup (z >= cut.zero)")
  }
  AHR_opt <- exp(mean(dfs_opt$loghr_po))

  # Calculate empirical AHR if not provided
  if (is.null(ahr_empirical)) {
    ahr_empirical <- exp(mean(df_work$loghr_po))
  }

  # Generate z-axis grid
  zpoints <- seq(min(z_values), max(z_values), by = zpoints_by)

  # Initialize result vectors
  HR_zpoints <- rep(NA, length(zpoints))
  HRminus_zpoints <- rep(NA, length(zpoints))
  HR2_zpoints <- rep(NA, length(zpoints))
  HRminus2_zpoints <- rep(NA, length(zpoints))

  # Check if theta columns exist for alternative calculations
  has_theta <- all(c("theta_1", "theta_0") %in% names(df_work))

  # Calculate AHR for each threshold
  for (zindex in seq_along(zpoints)) {
    zz <- zpoints[zindex]

    # Population with z >= zz
    dfz <- df_work[z_values >= zz, ]
    if (nrow(dfz) > 0) {
      HR_zpoints[zindex] <- exp(mean(dfz$loghr_po))

      if (has_theta) {
        aa <- mean(exp(dfz$theta_1))
        bb <- mean(exp(dfz$theta_0))
        HR2_zpoints[zindex] <- aa / bb
      }
    }

    # Population with z <= zz
    dfz_minus <- df_work[z_values <= zz, ]
    if (nrow(dfz_minus) > 0) {
      HRminus_zpoints[zindex] <- exp(mean(dfz_minus$loghr_po))

      if (has_theta) {
        aa <- mean(exp(dfz_minus$theta_1))
        bb <- mean(exp(dfz_minus$theta_0))
        HRminus2_zpoints[zindex] <- aa / bb
      }
    }
  }

  # Set up plotting layout
  if (plot_type == "both") {
    par(mfrow = c(1, 2))
  } else {
    par(mfrow = c(1, 1))
  }

  # Plot 1: Treatment effect profile
  if (plot_type %in% c("both", "profile")) {
    plot(z_values, loghr_values,
         type = "s",
         lty = 1,
         xlab = paste("Biomarker (", z, ")", sep = ""),
         ylab = expression(paste("Log Hazard Ratio ", psi[0], "(z)")),
         main = "Treatment Effect Profile",
         ...)

    if (add_rug) {
      rug(z_values)
    }

    # Add reference lines
    if (!is.null(log.hrs)) {
      abline(h = log.hrs, col = "gray", lty = 2)
    }

    abline(h = log(ahr_empirical), lwd = 2, col = "orange", lty = 1)
    abline(h = hrz_crit, lwd = 2, col = "red", lty = 2)
    abline(v = cut.zero, lwd = 2, col = "blue", lty = 2)

    legend("topright",
           legend = c("Empirical AHR", "HR threshold", "Optimal cutpoint"),
           col = c("orange", "red", "blue"),
           lty = c(1, 2, 2),
           lwd = c(2, 2, 2),
           cex = 0.8, bty = "n")
  }

  # Plot 2: Average hazard ratio curve
  if (plot_type %in% c("both", "ahr")) {
    plot(zpoints, HR_zpoints,
         type = "s",
         xlab = paste("BM (", z, ")", sep = ""),
         ylab = "Average Hazard Ratio",
         main = paste("AHR for BM:", z, ">= z"),
         lty = 1,
         col = "darkblue",
         lwd = 2,
         ...)

    abline(h = 1, col = "red", lty = 2)
    abline(h = exp(hrz_crit), col = "red", lty = 2, lwd = 2)
    abline(v = cut.zero, col = "blue", lty = 2, lwd = 2)

    legend("topright",
           legend = c("AHR(z+)", "HR = 1", "Optimal cutpoint"),
           col = c("darkblue", "red", "blue"),
           lty = c(1, 2, 2),
           lwd = c(2, 1, 2),
           cex = 0.8, bty = "n")
  }

  # Reset plotting parameters
  par(mfrow = c(1, 1))

  # Return results
  results <- list(
    z_variable = z,
    cut.zero = cut.zero,
    AHR_opt = AHR_opt,
    zpoints = zpoints,
    AHR_zpos = HR_zpoints,
    AHR_zneg = HRminus_zpoints
  )

  if (has_theta) {
    results$CDE_zpos <- HR2_zpoints
    results$CDE_zneg <- HRminus2_zpoints
  }

  invisible(results)
}
