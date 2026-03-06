
#' Fit AFT Model with Optional Spline Treatment Effect
#' @keywords internal
fit_aft_model <- function(df_work, interaction_term, k_treat, k_inter,
                          verbose, spline_spec = NULL, set_var = NULL, beta_var = NULL) {

  # Prepare model matrix
  covariate_cols <- grep("^z_", names(df_work), value = TRUE)
  X <- as.matrix(df_work[, c("treat", covariate_cols)])

  # Add interaction term if needed
  if (!is.null(interaction_term)) {
    X <- cbind(X, treat_harm = interaction_term)
  }

  # Check if spline specification is provided
  if (!is.null(spline_spec)) {
    # Use spline-based treatment effect
    result <- fit_aft_model_spline(df_work, covariate_cols, interaction_term,
                                   k_treat, k_inter, spline_spec, verbose, set_var, beta_var)
  } else {
    # Use standard (non-spline) model
    result <- fit_aft_model_standard(df_work, covariate_cols, interaction_term,
                                     k_treat, k_inter, verbose, set_var, beta_var)
  }

  return(result)
}


#' Fit Standard AFT Model (Non-Spline)
#' @keywords internal
fit_aft_model_standard <- function(df_work, covariate_cols, interaction_term,
                                   k_treat, k_inter, verbose, set_var = NULL, beta_var = NULL) {

  # Prepare model matrix
  X <- as.matrix(df_work[, c("treat", covariate_cols)])

  # Add interaction term if needed
  if (!is.null(interaction_term)) {
    df_work$treat_harm <- interaction_term
    X <- cbind(X, treat_harm = interaction_term)
  }

  # Fit Weibull AFT model
  formula_str <- paste("Surv(y, event) ~ ",
                       paste(c("treat", covariate_cols), collapse = " + "))
  if (!is.null(interaction_term)) {
    formula_str <- paste(formula_str, "+ treat_harm")
  }

  fit_aft <- survreg(as.formula(formula_str),
                     data = df_work,
                     dist = "weibull")

  # Extract parameters
  mu <- coef(fit_aft)[1]  # Intercept
  tau <- fit_aft$scale
  gamma <- coef(fit_aft)[-1]  # Coefficients (excluding intercept)

  # Weibull parameterization
  b0 <- -gamma / tau

  # Apply effect modifiers to Weibull log(hazard-ratio) parameterization
  b0["treat"] <- k_treat * b0["treat"]
  if ("treat_harm" %in% names(b0)) {
    b0["treat_harm"] <- k_inter * b0["treat_harm"]
  }

  if(!is.null(set_var)){
  b0[set_var] <- beta_var
  }


  # Transform back to corresponding revised gamma
  gamma <- -b0 * tau

  if (verbose) {
    cat("\n=== Model Parameters (AFT, log(T)) ===\n")
    cat("Intercept (mu):", round(mu, 3), "\n")
    cat("Scale (tau):", round(tau, 3), "\n")
    cat("Treatment effect:", round(gamma["treat"], 3), "\n")
    if ("treat_harm" %in% names(gamma)) {
      cat("Interaction effect:", round(gamma["treat_harm"], 3), "\n")
    }
  }

  return(list(
    mu = mu,
    tau = tau,
    gamma = gamma,
    b0 = b0,
    spline_info = NULL
  ))
}


#' Fit AFT Model with Spline Treatment Effect
#' @keywords internal
fit_aft_model_spline <- function(df_work, covariate_cols, interaction_term,
                                 k_treat, k_inter, spline_spec, verbose, set_var = NULL, beta_var = NULL) {

  # Validate spline specification
  spline_spec <- validate_spline_spec(spline_spec, df_work)

  # Extract spline parameters
  spline_var <- spline_spec$var
  knot <- spline_spec$knot
  log_hrs <- spline_spec$log_hrs
  zeta <- spline_spec$zeta

  # Create spline variables
  df_work <- create_spline_variables(df_work, spline_var, knot)

  # Build formula with spline terms
  spline_terms <- c(spline_var,
                    paste0(spline_var, "_treat"),
                    paste0(spline_var, "_k"),
                    paste0(spline_var, "_k_treat"))

  formula_terms <- c("treat", covariate_cols, spline_terms)

  # Add subgroup interaction if present
  if (!is.null(interaction_term)) {
    df_work$treat_harm <- interaction_term
    formula_terms <- c(formula_terms, "treat_harm")
  }

  formula_str <- paste("Surv(y, event) ~",
                       paste(formula_terms, collapse = " + "))

  # Fit Weibull AFT model
  fit_aft <- survreg(as.formula(formula_str),
                     data = df_work,
                     dist = "weibull")

  # Extract parameters
  mu <- coef(fit_aft)[1]
  tau <- fit_aft$scale
  gamma <- coef(fit_aft)[-1]

  # Convert to Weibull hazard parameterization
  b0 <- -gamma / tau

  # Apply spline constraint to treatment effect
  b0 <- apply_spline_constraint(b0, spline_var, knot, zeta, log_hrs,
                                k_treat, verbose)

  # Apply k_inter to subgroup interaction if present
  if ("treat_harm" %in% names(b0)) {
    b0["treat_harm"] <- k_inter * b0["treat_harm"]
  }

  if(!is.null(set_var)){
    b0[set_var] <- beta_var
  }


  # Convert back to AFT parameterization
  gamma <- -b0 * tau

  if (verbose) {
    cat("\n=== Model Parameters (AFT with Spline, log(T)) ===\n")
    cat("Intercept (mu):", round(mu, 3), "\n")
    cat("Scale (tau):", round(tau, 3), "\n")
    cat("Spline variable:", spline_var, "\n")
    cat("  Knot:", round(knot, 3), "\n")
    cat("  Zeta:", round(zeta, 3), "\n")
    cat("Treatment effect at", spline_var, "= 0:",
        round(gamma["treat"], 3), "\n")
    cat("Treatment slope before knot:",
        round(gamma[paste0(spline_var, "_treat")], 3), "\n")
    cat("Treatment slope change after knot:",
        round(gamma[paste0(spline_var, "_k_treat")], 3), "\n")
    if ("treat_harm" %in% names(gamma)) {
      cat("Subgroup interaction effect:", round(gamma["treat_harm"], 3), "\n")
    }
  }

  # Store spline information
  spline_info <- list(
    var = spline_var,
    knot = knot,
    zeta = zeta,
    log_hrs = log_hrs,
    df_work = df_work  # Contains spline variables
  )

  return(list(
    mu = mu,
    tau = tau,
    gamma = gamma,
    b0 = b0,
    spline_info = spline_info
  ))
}


#' Validate Spline Specification
#' @keywords internal
validate_spline_spec <- function(spline_spec, df_work) {

  # Check required fields
  if (is.null(spline_spec$var)) {
    stop("spline_spec must include 'var' (variable name)")
  }

  if (!spline_spec$var %in% names(df_work)) {
    stop("Spline variable '", spline_spec$var, "' not found in data")
  }

  # Set defaults if not provided
  if (is.null(spline_spec$knot)) {
    spline_spec$knot <- median(df_work[[spline_spec$var]], na.rm = TRUE)
    message("Using median as default knot: ", round(spline_spec$knot, 2))
  }

  if (is.null(spline_spec$zeta)) {
    spline_spec$zeta <- quantile(df_work[[spline_spec$var]], 0.75,
                                 na.rm = TRUE)
    message("Using 75th percentile as default zeta: ",
            round(spline_spec$zeta, 2))
  }

  if (is.null(spline_spec$log_hrs)) {
    spline_spec$log_hrs <- log(c(0.75, 0.75, 0.75))
    message("Using default log_hrs: ",
            paste(round(spline_spec$log_hrs, 3), collapse = ", "))
  }

  if (length(spline_spec$log_hrs) != 3) {
    stop("log_hrs must be a vector of length 3")
  }

  return(spline_spec)
}


#' Create Spline Variables
#' @keywords internal
create_spline_variables <- function(df_work, spline_var, knot) {

  # Main variable should already exist (either original or z_ transformed)
  # Create spline basis variables
  df_work[[paste0(spline_var, "_treat")]] <-
    df_work[[spline_var]] * df_work$treat

  df_work[[paste0(spline_var, "_k")]] <-
    (df_work[[spline_var]] - knot) * ifelse(df_work[[spline_var]] > knot, 1, 0)

  df_work[[paste0(spline_var, "_k_treat")]] <-
    df_work[[paste0(spline_var, "_k")]] * df_work$treat

  return(df_work)
}


#' Apply Spline Constraint to Treatment Effect Coefficients
#' @keywords internal
apply_spline_constraint <- function(b0, spline_var, knot, zeta, log_hrs,
                                    k_treat, verbose) {

  # Extract target log-HRs
  loghr_0 <- log_hrs[1]
  loghr_knot <- log_hrs[2]
  loghr_zeta <- log_hrs[3]

  # Get coefficient names
  treat_name <- "treat"
  interact_name <- paste0(spline_var, "_treat")
  spline_interact_name <- paste0(spline_var, "_k_treat")

  # Calculate constrained coefficients
  # b0[treat] = log(HR) at spline_var = 0
  # b0[spline_var:treat] = slope of log(HR) before knot
  # b0[spline_var_k:treat] = change in slope after knot

  b0_new <- b0
  b0_new[treat_name] <- k_treat * loghr_0
  b0_new[interact_name] <- k_treat * (loghr_knot - loghr_0) / knot
  b0_new[spline_interact_name] <- k_treat *
    (loghr_zeta - loghr_0 - zeta * (loghr_knot - loghr_0) / knot) /
    (zeta - knot)

  if (verbose) {
    cat("\n=== Spline Constraint Applied ===\n")
    cat("Target log(HR) at", spline_var, "= 0:", round(loghr_0, 3), "\n")
    cat("Target log(HR) at", spline_var, "= knot (", knot, "):",
        round(loghr_knot, 3), "\n")
    cat("Target log(HR) at", spline_var, "= zeta (", zeta, "):",
        round(loghr_zeta, 3), "\n")
    cat("k_treat multiplier:", k_treat, "\n")
  }

  return(b0_new)
}


#' Plot Spline Treatment Effect Function
#'
#' @param dgm_result Result object from generate_aft_dgm_flex with spline
#' @param add_points Logical; add observed data points. Default TRUE
#' @export
plot_spline_treatment_effect <- function(dgm_result, add_points = TRUE) {

  if (is.null(dgm_result$model_params$spline_info)) {
    stop("Result object does not contain spline information")
  }

  spline_info <- dgm_result$model_params$spline_info
  b0 <- dgm_result$model_params$b0

  spline_var <- spline_info$var
  knot <- spline_info$knot
  zeta <- spline_info$zeta
  log_hrs <- spline_info$log_hrs

  # Get data range
  var_data <- dgm_result$df_super[[spline_var]]
  var_range <- range(var_data, na.rm = TRUE)

  # Create sequence for plotting
  z_seq <- seq(var_range[1], var_range[2], length.out = 200)

  # Calculate log(HR) across range
  treat_name <- "treat"
  interact_name <- paste0(spline_var, "_treat")
  spline_interact_name <- paste0(spline_var, "_k_treat")

  loghr_seq <- b0[treat_name] +
    b0[interact_name] * z_seq +
    b0[spline_interact_name] * (z_seq - knot) *
    ifelse(z_seq > knot, 1, 0)

  # Create plot
  plot(z_seq, loghr_seq, type = "l", lwd = 2,
       xlab = spline_var, ylab = "log(HR)",
       main = paste("Treatment Effect Function:", spline_var))

  # Add reference lines
  abline(h = log_hrs, lty = 2, col = "gray60", lwd = 1)
  abline(v = c(knot, zeta), lty = 2, col = "blue", lwd = 1)
  abline(h = 0, lty = 1, col = "red", lwd = 0.5)

  # Add points if requested
  if (add_points && !is.null(dgm_result$df_super$loghr_po)) {
    # Calculate actual log(HR) for each observation
    df_plot <- dgm_result$df_super
    points(df_plot[[spline_var]], df_plot$loghr_po,
           pch = 16, cex = 0.3, col = rgb(0, 0, 0, 0.3))
  }

  # Add legend
  legend("topright",
         legend = c("Treatment effect", "Target log(HRs)",
                    "Knot & Zeta", "Null effect"),
         lty = c(1, 2, 2, 1),
         lwd = c(2, 1, 1, 0.5),
         col = c("black", "gray60", "blue", "red"), bty = "n")

  # Add labels for key points
  text(knot, par("usr")[3] + 0.05 * diff(par("usr")[3:4]),
       labels = "knot", pos = 3, cex = 0.8, col = "blue")
  text(zeta, par("usr")[3] + 0.05 * diff(par("usr")[3:4]),
       labels = "zeta", pos = 3, cex = 0.8, col = "blue")
}



#' Fit AFT Model and Apply Effect Modifiers
#' @keywords internal
fit_aft_model_legacy <- function(df_work, interaction_term, k_treat, k_inter,
                          verbose) {

  # Prepare model matrix
  covariate_cols <- grep("^z_", names(df_work), value = TRUE)
  X <- as.matrix(df_work[, c("treat", covariate_cols)])

  # Add interaction term if needed
  if (!is.null(interaction_term)) {
    X <- cbind(X, treat_harm = interaction_term)
  }

  # Fit Weibull AFT model
  formula_str <- paste("Surv(y, event) ~ ",
                       paste(c("treat", covariate_cols), collapse = " + "))
  if (!is.null(interaction_term)) {
    df_work$treat_harm <- interaction_term
    formula_str <- paste(formula_str, "+ treat_harm")
  }

  fit_aft <- survreg(as.formula(formula_str),
                     data = df_work,
                     dist = "weibull")

  # Extract parameters
  mu <- coef(fit_aft)[1]  # Intercept
  tau <- fit_aft$scale
  gamma <- coef(fit_aft)[-1]  # Coefficients (excluding intercept)

  # Weibull parameterization
  b0 <- -gamma / tau

  # Apply effect modifiers Weibull log(hazard-ratio) parameterization
  b0["treat"] <- k_treat * b0["treat"]
  if ("treat_harm" %in% names(b0)) {
    b0["treat_harm"] <- k_inter * b0["treat_harm"]
  }

  # Transform to corresponding revised gamma
  gamma <- -b0 * tau

  if (verbose) {
    cat("\n=== Model Parameters (AFT, log(T)) ===\n")
    cat("Intercept (mu):", round(mu, 3), "\n")
    cat("Scale (tau):", round(tau, 3), "\n")
    cat("Treatment effect:", round(gamma["treat"], 3), "\n")
    if ("treat_harm" %in% names(gamma)) {
      cat("Interaction effect:", round(gamma["treat_harm"], 3), "\n")
    }
  }

  return(list(mu = mu, tau = tau, gamma = gamma, b0 = b0))
}
