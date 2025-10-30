#' Additional S3 Methods for AFT DGM Objects
#'
#' This file provides S3 methods for objects created by the AFT DGM functions.
#' These methods enhance the usability of the package by providing standard
#' print and summary functions for the data generating mechanism objects.
#'
#' @name aft_dgm_methods
#' @keywords internal

# ================================================================================
# S3 Methods for aft_dgm_flex Objects
# ================================================================================

#' Print Method for AFT DGM Flex Objects
#' 
#' Provides a formatted summary of an AFT data generating mechanism object.
#'
#' @param x An object of class "aft_dgm_flex"
#' @param digits Number of digits for numeric output. Default is 3.
#' @param ... Additional arguments (currently unused)
#'
#' @return Invisibly returns the input object
#' 
#' @examples
#' \dontrun{
#' dgm <- generate_aft_dgm_flex(data = gbsg, ...)
#' print(dgm)
#' }
#' 
#' @export
#' @method print aft_dgm_flex
print.aft_dgm_flex <- function(x, digits = 3, ...) {
  cat("AFT Data Generating Mechanism\n")
  cat("=============================\n")
  cat("Model type:", x$model_type, "\n")
  cat("Super population size:", x$n_super, "\n")
  cat("Number of covariates:", length(x$analysis_vars$covariates), "\n")
  
  if (x$model_type == "alt" && !is.null(x$subgroup_info$vars)) {
    cat("\nSubgroup Information:\n")
    cat("  Variables:", paste(x$subgroup_info$vars, collapse = ", "), "\n")
    cat("  Size:", x$subgroup_info$size, "\n")
    cat("  Proportion:", round(x$subgroup_info$proportion, digits), "\n")
    
    if (!is.null(x$subgroup_info$definitions)) {
      cat("  Definitions:\n")
      for (def in x$subgroup_info$definitions) {
        cat("    -", def, "\n")
      }
    }
  }
  
  cat("\nHazard Ratios:\n")
  for (name in names(x$hazard_ratios)) {
    cat("  ", gsub("_", " ", name), ":", round(x$hazard_ratios[[name]], digits), "\n")
  }
  
  cat("\nModel Parameters:\n")
  cat("  Intercept (mu):", round(x$model_params$mu, digits), "\n")
  cat("  Scale (sigma):", round(x$model_params$sigma, digits), "\n")
  if ("treat" %in% names(x$model_params$gamma)) {
    cat("  Treatment effect:", round(x$model_params$gamma["treat"], digits), "\n")
  }
  if ("treat_harm" %in% names(x$model_params$gamma)) {
    cat("  Interaction effect:", round(x$model_params$gamma["treat_harm"], digits), "\n")
  }
  
  invisible(x)
}

#' Summary Method for AFT DGM Flex Objects
#' 
#' Provides a comprehensive summary of an AFT data generating mechanism object.
#'
#' @param object An object of class "aft_dgm_flex"
#' @param ... Additional arguments (currently unused)
#'
#' @return A list of class "summary.aft_dgm_flex" containing summary information
#' 
#' @examples
#' \dontrun{
#' dgm <- generate_aft_dgm_flex(data = gbsg, ...)
#' summary(dgm)
#' }
#' 
#' @export
#' @method summary aft_dgm_flex
summary.aft_dgm_flex <- function(object, ...) {
  summary_list <- list(
    model_type = object$model_type,
    n_super = object$n_super,
    n_covariates = length(object$analysis_vars$covariates),
    continuous_vars = object$analysis_vars$continuous,
    factor_vars = object$analysis_vars$factor,
    subgroup_info = object$subgroup_info,
    hazard_ratios = object$hazard_ratios,
    model_params = object$model_params,
    seed = object$seed
  )
  
  class(summary_list) <- "summary.aft_dgm_flex"
  return(summary_list)
}

#' Print Method for Summary of AFT DGM Flex
#' 
#' @param x An object of class "summary.aft_dgm_flex"
#' @param digits Number of digits for numeric output
#' @param ... Additional arguments (ignored)
#' @return Invisible NULL
#' @export
#' @method print summary.aft_dgm_flex
print.summary.aft_dgm_flex <- function(x, digits = 3, ...) {
  cat("\nSummary of AFT Data Generating Mechanism\n")
  cat("=========================================\n\n")
  
  cat("Model Configuration:\n")
  cat("  Type:", x$model_type, "\n")
  cat("  Super population:", x$n_super, "\n")
  cat("  Number of covariates:", x$n_covariates, "\n")
  cat("  Random seed:", x$seed, "\n")
  
  cat("\nVariables:\n")
  cat("  Continuous:", paste(x$continuous_vars, collapse = ", "), "\n")
  cat("  Factor:", paste(x$factor_vars, collapse = ", "), "\n")
  
  if (!is.null(x$subgroup_info$vars)) {
    cat("\nSubgroup Definition:\n")
    cat("  Variables:", paste(x$subgroup_info$vars, collapse = ", "), "\n")
    cat("  Size:", x$subgroup_info$size, "(", 
        round(100 * x$subgroup_info$proportion, 1), "%)\n", sep = "")
  }
  
  cat("\nHazard Ratios:\n")
  hr_df <- data.frame(
    Population = gsub("_", " ", names(x$hazard_ratios)),
    HR = round(unlist(x$hazard_ratios), digits)
  )
  print(hr_df, row.names = FALSE)
  
  cat("\nKey Model Parameters:\n")
  cat("  Intercept (μ):", round(x$model_params$mu, digits), "\n")
  cat("  Scale (σ):", round(x$model_params$sigma, digits), "\n")
  
  invisible(NULL)
}

# ================================================================================
# Additional Utility S3 Methods
# ================================================================================

#' Coef Method for k_inter Results
#' 
#' Extract the optimal k_inter coefficient from calibration results.
#'
#' @param object An object of class "k_inter_result", "k_inter_grid_result", 
#'   or "k_inter_optim_result"
#' @param ... Additional arguments (ignored)
#' @return Numeric value of optimal k_inter
#' @export
#' @method coef k_inter_result
coef.k_inter_result <- function(object, ...) {
  return(object$k_inter)
}

#' @export
#' @method coef k_inter_grid_result
coef.k_inter_grid_result <- function(object, ...) {
  return(object$k_inter)
}

#' @export
#' @method coef k_inter_optim_result
coef.k_inter_optim_result <- function(object, ...) {
  return(object$k_inter)
}

#' Predict Method for k_inter Sensitivity Results
#' 
#' Predict hazard ratios for new k_inter values based on sensitivity analysis.
#'
#' @param object An object of class "k_inter_sensitivity"
#' @param newdata Numeric vector of new k_inter values
#' @param type Character string specifying which HR to predict: 
#'   "harm", "no_harm", or "overall" (default)
#' @param ... Additional arguments (ignored)
#' @return Numeric vector of predicted hazard ratios
#' @export
#' @method predict k_inter_sensitivity
predict.k_inter_sensitivity <- function(object, newdata, type = "overall", ...) {
  if (!type %in% c("harm", "no_harm", "overall")) {
    stop("type must be 'harm', 'no_harm', or 'overall'")
  }
  
  # Use linear interpolation
  hr_column <- switch(type,
                     harm = "hr_harm",
                     no_harm = "hr_no_harm",
                     overall = "hr_overall")
  
  predicted <- approx(x = object$k_inter, 
                     y = object[[hr_column]], 
                     xout = newdata)$y
  
  return(predicted)
}

#' Plot Method for k_inter Sensitivity Results
#' 
#' Create visualization of sensitivity analysis results.
#'
#' @param x An object of class "k_inter_sensitivity"
#' @param type Character string: "all" (default), "harm", "ratio", or "comparison"
#' @param ... Additional graphical parameters
#' @return Invisible NULL
#' @export
#' @method plot k_inter_sensitivity
plot.k_inter_sensitivity <- function(x, type = "all", ...) {
  
  if (type == "all") {
    # Create 4-panel plot as in original
    oldpar <- par(mfrow = c(2, 2))
    on.exit(par(oldpar))
    
    # Panel 1: Harm subgroup
    plot(x$k_inter, x$hr_harm, type = "l", lwd = 2, col = "red",
         xlab = "k_inter", ylab = "Hazard Ratio",
         main = "Harm Subgroup HR vs k_inter", ...)
    abline(h = 1, lty = 2, col = "gray")
    grid()
    
    # Panel 2: All groups
    plot(x$k_inter, x$hr_overall, type = "l", lwd = 2,
         xlab = "k_inter", ylab = "Hazard Ratio",
         main = "All Hazard Ratios vs k_inter",
         ylim = range(c(x$hr_harm, x$hr_no_harm, x$hr_overall)), ...)
    lines(x$k_inter, x$hr_harm, col = "red", lwd = 2)
    lines(x$k_inter, x$hr_no_harm, col = "blue", lwd = 2)
    abline(h = 1, lty = 2, col = "gray")
    legend("topright", legend = c("Overall", "Harm", "No-harm"),
           col = c("black", "red", "blue"), lty = 1, lwd = 2)
    grid()
    
    # Panel 3: Ratio
    plot(x$k_inter, x$hr_harm / x$hr_no_harm, type = "l", lwd = 2,
         xlab = "k_inter", ylab = "HR Ratio (Harm/No-harm)",
         main = "Effect Modification", ...)
    abline(h = 1, lty = 2, col = "gray")
    grid()
    
    # Panel 4: Summary table
    plot.new()
    text(0.5, 0.9, "Summary", cex = 1.2, font = 2)
    text(0.5, 0.75, paste("k_inter range:", 
                         round(range(x$k_inter), 2)[1], "to",
                         round(range(x$k_inter), 2)[2]), cex = 0.9)
    text(0.5, 0.65, paste("HR harm range:",
                         round(range(x$hr_harm), 2)[1], "to",
                         round(range(x$hr_harm), 2)[2]), cex = 0.9)
    text(0.5, 0.55, paste("HR overall range:",
                         round(range(x$hr_overall), 2)[1], "to",
                         round(range(x$hr_overall), 2)[2]), cex = 0.9)
    
  } else if (type == "harm") {
    plot(x$k_inter, x$hr_harm, type = "l", lwd = 2, col = "red",
         xlab = "k_inter", ylab = "Hazard Ratio",
         main = "Harm Subgroup HR vs k_inter", ...)
    abline(h = 1, lty = 2, col = "gray")
    grid()
    
  } else if (type == "ratio") {
    plot(x$k_inter, x$hr_harm / x$hr_no_harm, type = "l", lwd = 2,
         xlab = "k_inter", ylab = "HR Ratio (Harm/No-harm)",
         main = "Effect Modification", ...)
    abline(h = 1, lty = 2, col = "gray")
    grid()
    
  } else if (type == "comparison") {
    plot(x$k_inter, x$hr_overall, type = "l", lwd = 2,
         xlab = "k_inter", ylab = "Hazard Ratio",
         main = "Hazard Ratios Comparison",
         ylim = range(c(x$hr_harm, x$hr_no_harm, x$hr_overall)), ...)
    lines(x$k_inter, x$hr_harm, col = "red", lwd = 2)
    lines(x$k_inter, x$hr_no_harm, col = "blue", lwd = 2)
    abline(h = 1, lty = 2, col = "gray")
    legend("topright", legend = c("Overall", "Harm", "No-harm"),
           col = c("black", "red", "blue"), lty = 1, lwd = 2)
    grid()
  }
  
  invisible(NULL)
}

#' As.data.frame Method for k_inter Batch Results
#' 
#' Convert k_inter batch results to a standard data frame.
#'
#' @param x An object of class "k_inter_batch_result"
#' @param row.names NULL or character vector giving row names
#' @param optional Logical. If TRUE, setting row names is optional
#' @param ... Additional arguments (ignored)
#' @return A data.frame
#' @export
#' @method as.data.frame k_inter_batch_result
as.data.frame.k_inter_batch_result <- function(x, row.names = NULL, 
                                               optional = FALSE, ...) {
  class(x) <- "data.frame"
  if (!is.null(row.names)) {
    rownames(x) <- row.names
  }
  return(x)
}
