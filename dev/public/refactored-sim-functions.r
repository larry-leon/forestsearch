# ==============================================================================
# REFACTORED SIMULATION FUNCTIONS
# Complete implementation with modular design
# ==============================================================================

# ------------------------------------------------------------------------------
# CONSTANTS AND CONFIGURATION
# ------------------------------------------------------------------------------

SIMULATION_CONSTANTS <- list(
  BASE_SEED = 8316951,
  SEED_INCREMENT = 1000,
  CENSORING_SEED_OFFSET = 100,
  MAX_KM_PLOTS = 10,
  DEFAULT_Z_INCREMENT = 1,
  DEFAULT_KNOT = 5,
  DEFAULT_ZETA = 10
)

# ------------------------------------------------------------------------------
# MAIN SIMULATION FUNCTION (ORCHESTRATOR)
# ------------------------------------------------------------------------------

#' Draw Stratified Simulation (Refactored)
#'
#' @param dgm Data generating mechanism from get_dgm_stratified()
#' @param ss Simulation seed/index
#' @param Ndraw Number of observations to draw
#' @param strata_rand Randomization stratification variable
#' @param wname Prognostic covariate name
#' @param bw Log hazard ratio for prognostic covariate
#' @param checking Logical. Validate simulation parameters?
#' @param details Logical. Print detailed output?
#' @param return_df Logical. Return dataframe (TRUE) or summaries (FALSE)?
#' @param time_eos Time of end-of-study censoring
#' @param keep_rand Logical. Keep original randomization?
#' @param hrz_crit Log HR threshold for optimal biomarker cutpoint
#'
#' @return Either simulated dataframe or list of population summaries
#' @export
draw_sim_stratified <- function(dgm, 
                                ss = 1, 
                                Ndraw = nrow(dgm$df_super),
                                strata_rand = "stratum", 
                                wname = "ecog", 
                                bw = 0,
                                checking = FALSE, 
                                details = FALSE,
                                return_df = TRUE, 
                                time_eos = Inf, 
                                keep_rand = FALSE,
                                hrz_crit = log(1.2)) {
  
  # Validate inputs
  validate_simulation_inputs(dgm, strata_rand, wname)
  
  # 1. Prepare simulation data
  sim_data <- prepare_simulation_data(dgm, Ndraw, ss)
  
  # 2. Calculate adjusted gamma parameter for prognostic factor
  gamma_w <- calculate_gamma_w(bw, dgm$tau.approx)
  
  # 3. Generate potential outcomes
  po_data <- generate_potential_outcomes(
    df = sim_data,
    gamma_true = dgm$gamma.true,
    mu = dgm$mu,
    tau_strata = sim_data$tau.strataO,
    gamma_w = gamma_w,
    wname = wname,
    seed = ss
  )
  
  # 4. Apply randomization and treatment assignment
  obs_data <- apply_randomization_treatment(
    po_data = po_data,
    strata_rand = strata_rand,
    keep_rand = keep_rand,
    seed = ss
  )
  
  # 5. Apply censoring
  obs_data <- apply_censoring(
    obs_data = obs_data,
    muC = dgm$muC,
    tauC = dgm$tauC,
    time_eos = time_eos,
    seed = ss
  )
  
  # 6. Add stratification variables
  obs_data <- add_stratification_variables(
    obs_data = obs_data,
    strata_rand = strata_rand,
    strata_tte = dgm$strata_tte
  )
  
  # 7. Optional: Validation
  if (checking) {
    validate_simulation_parameters(obs_data, dgm, details)
  }
  
  # 8. Optional: Report metrics
  if (details && ss <= SIMULATION_CONSTANTS$MAX_KM_PLOTS) {
    report_simulation_metrics(obs_data, ss)
  }
  
  # 9. Format and return output
  return(format_simulation_output(obs_data, return_df, hrz_crit))
}

# ------------------------------------------------------------------------------
# VALIDATION FUNCTIONS
# ------------------------------------------------------------------------------

#' Validate simulation inputs
#' @keywords internal
validate_simulation_inputs <- function(dgm, strata_rand, wname) {
  required_vars <- c(strata_rand, wname)
  
  if (!all(required_vars %in% names(dgm$df_super))) {
    missing_vars <- required_vars[!required_vars %in% names(dgm$df_super)]
    stop(sprintf(
      "Required variables not in dgm$df_super: %s",
      paste(missing_vars, collapse = ", ")
    ), call. = FALSE)
  }
  
  required_dgm_elements <- c("gamma.true", "mu", "tau", "muC", "tauC", "tau.approx")
  if (!all(required_dgm_elements %in% names(dgm))) {
    missing_elements <- required_dgm_elements[!required_dgm_elements %in% names(dgm)]
    stop(sprintf(
      "Required DGM elements missing: %s",
      paste(missing_elements, collapse = ", ")
    ), call. = FALSE)
  }
  
  invisible(TRUE)
}

# ------------------------------------------------------------------------------
# DATA PREPARATION
# ------------------------------------------------------------------------------

#' Prepare simulation dataset with sampling
#' @keywords internal
prepare_simulation_data <- function(dgm, Ndraw, seed) {
  set.seed(SIMULATION_CONSTANTS$BASE_SEED + seed * SIMULATION_CONSTANTS$SEED_INCREMENT)
  
  df_super <- dgm$df_super
  
  # Sample if Ndraw differs from population size
  if (Ndraw != nrow(df_super)) {
    id_sample <- sample(seq_len(nrow(df_super)), size = Ndraw, replace = TRUE)
    df_super <- df_super[id_sample, ]
  }
  
  return(df_super)
}

#' Calculate gamma parameter for prognostic covariate
#' @keywords internal
calculate_gamma_w <- function(bw, tau_approx) {
  -bw * tau_approx
}

# ------------------------------------------------------------------------------
# POTENTIAL OUTCOMES GENERATION
# ------------------------------------------------------------------------------

#' Generate potential outcomes under both treatment conditions
#' @keywords internal
generate_potential_outcomes <- function(df, gamma_true, mu, tau_strata, 
                                       gamma_w, wname, seed) {
  N <- nrow(df)
  
  # Extract prognostic covariate
  w <- df[[wname]]
  
  # Build covariate matrices for both treatment conditions
  zmat_1 <- build_covariate_matrix(df, treat = 1)
  zmat_0 <- build_covariate_matrix(df, treat = 0)
  
  # Generate error term (on log scale)
  set.seed(SIMULATION_CONSTANTS$BASE_SEED + seed * SIMULATION_CONSTANTS$SEED_INCREMENT + 1)
  epsilon <- log(rexp(N))
  
  # Calculate log survival times under treatment
  df$log.Y1 <- mu + c(zmat_1 %*% gamma_true) + w * gamma_w + tau_strata * epsilon
  
  # Calculate log survival times under control
  df$log.Y0 <- mu + c(zmat_0 %*% gamma_true) + w * gamma_w + tau_strata * epsilon
  
  # Calculate treatment effects (log hazard ratios)
  df$loghr.po <- calculate_log_hazard_ratio(zmat_1, zmat_0, gamma_true, tau_strata)
  
  # Calculate hazard components (for CDE calculations)
  df <- calculate_hazard_components(df, zmat_1, zmat_0, gamma_true, gamma_w, w, tau_strata)
  
  return(df)
}

#' Build covariate matrix for given treatment condition
#' @keywords internal
build_covariate_matrix <- function(df, treat) {
  # Extract base covariate matrix
  zmat <- as.matrix(df[, c("treat", "z", "z.treat", "z.k", "z.k.treat")])
  
  # Set treatment and interactions
  zmat[, "treat"] <- treat
  zmat[, "z.treat"] <- treat * df$z
  zmat[, "z.k.treat"] <- treat * zmat[, "z.k"]
  
  return(zmat)
}

#' Calculate log hazard ratio from covariate matrices
#' @keywords internal
calculate_log_hazard_ratio <- function(zmat_1, zmat_0, gamma_true, tau_strata) {
  (-1) * c((zmat_1 - zmat_0) %*% gamma_true) / tau_strata
}

#' Calculate hazard components for CDE (Controlled Direct Effect)
#' @keywords internal
calculate_hazard_components <- function(df, zmat_1, zmat_0, gamma_true, 
                                       gamma_w, w, tau_strata) {
  # Theta = exp(linear predictor / tau)
  # These are log(theta) values
  df$theta1.po <- -c(zmat_1 %*% gamma_true + w * gamma_w) / tau_strata
  df$theta0.po <- -c(zmat_0 %*% gamma_true + w * gamma_w) / tau_strata
  
  # Also calculate with W set to specific values (for subgroup analyses)
  df$theta1_w1.po <- -c(zmat_1 %*% gamma_true + 1 * gamma_w) / tau_strata
  df$theta0_w1.po <- -c(zmat_0 %*% gamma_true + 1 * gamma_w) / tau_strata
  
  df$theta1_w0.po <- -c(zmat_1 %*% gamma_true + 0 * gamma_w) / tau_strata
  df$theta0_w0.po <- -c(zmat_0 %*% gamma_true + 0 * gamma_w) / tau_strata
  
  return(df)
}

# ------------------------------------------------------------------------------
# RANDOMIZATION AND TREATMENT ASSIGNMENT
# ------------------------------------------------------------------------------

#' Apply randomization and assign treatment
#' @keywords internal
apply_randomization_treatment <- function(po_data, strata_rand, keep_rand, seed) {
  if (!keep_rand) {
    # Perform block randomization by strata
    blocks <- po_data[[strata_rand]]
    
    set.seed(SIMULATION_CONSTANTS$BASE_SEED + seed * SIMULATION_CONSTANTS$SEED_INCREMENT + 2)
    po_data$treat.sim <- randomizr::block_ra(blocks = blocks)
  } else {
    # Keep original treatment assignment
    po_data$treat.sim <- ifelse(po_data$treat == 1, 1, 0)
  }
  
  # Generate observed outcomes based on treatment assignment
  po_data$y.sim <- exp(
    po_data$treat.sim * po_data$log.Y1 + 
    (1 - po_data$treat.sim) * po_data$log.Y0
  )
  
  return(po_data)
}

# ------------------------------------------------------------------------------
# CENSORING
# ------------------------------------------------------------------------------

#' Apply censoring mechanism
#' @keywords internal
apply_censoring <- function(obs_data, muC, tauC, time_eos, seed) {
  N <- nrow(obs_data)
  
  # Generate censoring times
  set.seed(SIMULATION_CONSTANTS$BASE_SEED + 
           seed * SIMULATION_CONSTANTS$SEED_INCREMENT + 
           SIMULATION_CONSTANTS$CENSORING_SEED_OFFSET)
  
  epsilonC <- log(rexp(N))
  censor_time <- exp(muC + tauC * epsilonC)
  
  # Apply administrative censoring at end of study
  censor_time <- pmin(censor_time, time_eos)
  
  # Determine event indicator and observed time
  obs_data$event.sim <- ifelse(obs_data$y.sim <= censor_time, 1, 0)
  obs_data$y.sim <- pmin(obs_data$y.sim, censor_time)
  
  return(obs_data)
}

# ------------------------------------------------------------------------------
# STRATIFICATION VARIABLES
# ------------------------------------------------------------------------------

#' Add stratification variables to output
#' @keywords internal
add_stratification_variables <- function(obs_data, strata_rand, strata_tte) {
  # Randomization stratification
  obs_data$strata.simR <- obs_data[[strata_rand]]
  
  # Outcome stratification
  if (!is.null(strata_tte)) {
    obs_data$strata.simO <- obs_data[[strata_tte]]
  } else {
    obs_data$strata.simO <- "All"
  }
  
  return(obs_data)
}

# ------------------------------------------------------------------------------
# VALIDATION AND CHECKING
# ------------------------------------------------------------------------------

#' Validate simulation parameters against super-population
#' @keywords internal
validate_simulation_parameters <- function(obs_data, dgm, details) {
  if (details) {
    cat("\n=== Validation Checks ===\n")
    cat("Stratification parameters (tau) from dgm:\n")
    print(dgm$tau)
  }
  
  # Fit Weibull model to simulated data
  if (!is.null(dgm$strata_tte)) {
    strata_formula <- paste("strata(", "strata.simO", ")")
    formula_str <- paste(
      "Surv(y.sim, event.sim) ~ treat.sim + z + z.treat + z.k + z.k.treat + w +",
      strata_formula
    )
  } else {
    formula_str <- "Surv(y.sim, event.sim) ~ treat.sim + z + z.treat + z.k + z.k.treat + w"
  }
  
  fit_weib <- survival::survreg(
    as.formula(formula_str),
    dist = 'weibull',
    data = obs_data
  )
  
  tau_sim <- fit_weib$scale
  
  if (details) {
    cat("Stratification parameters (tau) from simulated data:\n")
    print(tau_sim)
    
    # Check log(HR) calculation
    max_diff <- max(abs(with(obs_data, loghr.po - (log.Y0 - log.Y1) / tau.strataO)))
    cat(sprintf("\nMax |loghr.po - (log.Y0-log.Y1)/tau| = %.12f\n", max_diff))
    
    # Compare Weibull and Cox estimates
    bhat_weib <- -(1) * coef(fit_weib)[-1] / tau_sim
    
    fit_cox <- survival::coxph(as.formula(formula_str), data = obs_data)
    
    comparison <- data.frame(
      Parameter = names(bhat_weib),
      Weibull = as.numeric(bhat_weib),
      Cox = as.numeric(coef(fit_cox))
    )
    
    cat("\nWeibull vs Cox parameter estimates:\n")
    print(comparison, row.names = FALSE, digits = 6)
  }
  
  invisible(TRUE)
}

# ------------------------------------------------------------------------------
# REPORTING
# ------------------------------------------------------------------------------

#' Report simulation metrics
#' @keywords internal
report_simulation_metrics <- function(obs_data, ss) {
  censoring_rate <- mean(1 - obs_data$event.sim)
  
  cat(sprintf("\n=== Simulation %d Metrics ===\n", ss))
  cat(sprintf("Sample size: %d\n", nrow(obs_data)))
  cat(sprintf("Censoring rate: %.3f\n", censoring_rate))
  cat(sprintf("Events: %d (%.1f%%)\n", 
              sum(obs_data$event.sim),
              100 * mean(obs_data$event.sim)))
  
  # Treatment assignment balance
  cat(sprintf("Treatment assignment: %d (%.1f%%) vs %d (%.1f%%)\n",
              sum(obs_data$treat.sim),
              100 * mean(obs_data$treat.sim),
              sum(1 - obs_data$treat.sim),
              100 * mean(1 - obs_data$treat.sim)))
  
  invisible(NULL)
}

# ------------------------------------------------------------------------------
# OUTPUT FORMATTING
# ------------------------------------------------------------------------------

#' Format simulation output
#' @keywords internal
format_simulation_output <- function(obs_data, return_df, hrz_crit) {
  # Sort by biomarker
  obs_data <- data.table::setorder(obs_data, z)
  
  if (!return_df) {
    # Return population-level summaries
    return(calculate_population_summaries(obs_data, hrz_crit))
  }
  
  # Return full dataset
  return(obs_data)
}

#' Calculate population-level summaries
#' @keywords internal
calculate_population_summaries <- function(df, hrz_crit) {
  # Overall AHR (Average Hazard Ratio)
  ahr_overall <- exp(mean(df$loghr.po))
  
  # AHR by prognostic factor W
  ahr_w1 <- with(subset(df, w == 1), exp(mean(loghr.po)))
  ahr_w0 <- with(subset(df, w == 0), exp(mean(loghr.po)))
  
  # CDE (Controlled Direct Effect) - averaging hazards not log(hazards)
  cde_overall <- mean(exp(df$theta1.po)) / mean(exp(df$theta0.po))
  cde_w1 <- with(subset(df, w == 1), mean(exp(theta1_w1.po)) / mean(exp(theta0_w1.po)))
  cde_w0 <- with(subset(df, w == 0), mean(exp(theta1_w0.po)) / mean(exp(theta0_w0.po)))
  
  # Fit Cox models for ITT analyses
  cox_results <- fit_itt_cox_models(df)
  
  # Calculate optimal biomarker cutpoint
  cut_optimal <- with(df, min(z[which(loghr.po < hrz_crit)]))
  ahr_optimal <- with(subset(df, z >= cut_optimal), exp(mean(loghr.po)))
  
  # Calculate AHR across biomarker thresholds
  z_thresholds <- seq(min(df$z), max(df$z), by = SIMULATION_CONSTANTS$DEFAULT_Z_INCREMENT)
  ahr_by_threshold <- calculate_ahr_by_threshold(df, z_thresholds)
  
  # Compile results
  list(
    # Overall measures
    AHR = ahr_overall,
    CDE = cde_overall,
    
    # By prognostic factor
    AHR_W1 = ahr_w1,
    AHR_W0 = ahr_w0,
    CDE_W1 = cde_w1,
    CDE_W0 = cde_w0,
    
    # Cox ITT analyses
    ITT_unadjusted = cox_results$unadjusted,
    ITT_strata_R = cox_results$strata_R,
    ITT_strata_R_W = cox_results$strata_R_W,
    W1_subpop = cox_results$W1_subpop,
    W0_subpop = cox_results$W0_subpop,
    
    # Optimal biomarker
    cut_optimal = cut_optimal,
    AHR_optimal = ahr_optimal,
    
    # Threshold-based analyses
    zpoints = z_thresholds,
    HR.zpoints = ahr_by_threshold$ahr_above,
    HRminus.zpoints = ahr_by_threshold$ahr_below,
    HR2.zpoints = ahr_by_threshold$cde_above,
    HRminus2.zpoints = ahr_by_threshold$cde_below
  )
}

#' Fit ITT Cox models with different adjustment strategies
#' @keywords internal
fit_itt_cox_models <- function(df) {
  # Unadjusted
  fit_unadj <- survival::coxph(
    Surv(y.sim, event.sim) ~ treat.sim,
    data = df
  )
  
  # Stratified by randomization strata
  fit_sR <- survival::coxph(
    Surv(y.sim, event.sim) ~ treat.sim + strata(strata.simR),
    data = df
  )
  
  # Adjusted for W and stratified by R
  fit_sRW <- survival::coxph(
    Surv(y.sim, event.sim) ~ treat.sim + w + strata(strata.simR),
    data = df
  )
  
  # Subpopulation analyses
  fit_w1 <- survival::coxph(
    Surv(y.sim, event.sim) ~ treat.sim,
    data = subset(df, w == 1)
  )
  
  fit_w0 <- survival::coxph(
    Surv(y.sim, event.sim) ~ treat.sim,
    data = subset(df, w == 0)
  )
  
  list(
    unadjusted = summary(fit_unadj)$conf.int[1, 1],
    strata_R = summary(fit_sR)$conf.int[1, 1],
    strata_R_W = summary(fit_sRW)$conf.int[1, 1],
    W1_subpop = summary(fit_w1)$conf.int[1, 1],
    W0_subpop = summary(fit_w0)$conf.int[1, 1]
  )
}

#' Calculate AHR for populations above/below biomarker thresholds
#' @keywords internal
calculate_ahr_by_threshold <- function(df, z_thresholds) {
  ahr_above <- numeric(length(z_thresholds))
  ahr_below <- numeric(length(z_thresholds))
  cde_above <- numeric(length(z_thresholds))
  cde_below <- numeric(length(z_thresholds))
  
  for (i in seq_along(z_thresholds)) {
    z_cut <- z_thresholds[i]
    
    # Populations
    df_above <- subset(df, z >= z_cut)
    df_below <- subset(df, z <= z_cut)
    
    # AHR (averaging log hazard ratios)
    ahr_above[i] <- with(df_above, exp(mean(loghr.po)))
    ahr_below[i] <- with(df_below, exp(mean(loghr.po)))
    
    # CDE (averaging hazards)
    cde_above[i] <- with(df_above, mean(exp(theta1.po)) / mean(exp(theta0.po)))
    cde_below[i] <- with(df_below, mean(exp(theta1.po)) / mean(exp(theta0.po)))
  }
  
  list(
    ahr_above = ahr_above,
    ahr_below = ahr_below,
    cde_above = cde_above,
    cde_below = cde_below
  )
}

# ------------------------------------------------------------------------------
# PLOTTING HELPER
# ------------------------------------------------------------------------------

#' Plot AHR profiles
#' @param popsummary Output from calculate_population_summaries
#' @param dfcase Original case study dataset
#' @export
plot_AHRs <- function(popsummary, dfcase) {
  par(mfrow = c(1, 2))
  
  # AHR(z+) - populations with biomarker >= z
  ymin <- min(c(popsummary$HR.zpoints, popsummary$HR2.zpoints), na.rm = TRUE)
  ymax <- max(c(popsummary$HR.zpoints, popsummary$HR2.zpoints), na.rm = TRUE)
  
  plot(popsummary$zpoints, popsummary$HR.zpoints,
       xlab = "Biomarker (z)", 
       ylab = "Average Hazard Ratio",
       type = "s", lty = 1, col = "black", lwd = 2,
       ylim = c(ymin, ymax),
       main = "AHR(z+): Populations with z >= threshold")
  
  lines(popsummary$zpoints, popsummary$HR2.zpoints,
        type = "s", lty = 2, col = "blue", lwd = 2)
  
  rug(jitter(dfcase$z), col = "gray")
  
  legend("topright", 
         c("AHR (avg log HR)", "CDE (avg HR)"),
         lty = c(1, 2), 
         col = c("black", "blue"), 
         lwd = 2, 
         bty = "n")
  
  # AHR(z-) - populations with biomarker <= z
  ymin <- min(c(popsummary$HRminus.zpoints, popsummary$HRminus2.zpoints), na.rm = TRUE)
  ymax <- max(c(popsummary$HRminus.zpoints, popsummary$HRminus2.zpoints), na.rm = TRUE)
  
  plot(popsummary$zpoints, popsummary$HRminus.zpoints,
       xlab = "Biomarker (z)", 
       ylab = "Average Hazard Ratio",
       type = "s", lty = 1, col = "black", lwd = 2,
       ylim = c(ymin, ymax),
       main = "AHR(z-): Populations with z <= threshold")
  
  lines(popsummary$zpoints, popsummary$HRminus2.zpoints,
        type = "s", lty = 2, col = "blue", lwd = 2)
  
  rug(jitter(dfcase$z), col = "gray")
  
  legend("topright", 
         c("AHR (avg log HR)", "CDE (avg HR)"),
         lty = c(1, 2), 
         col = c("black", "blue"), 
         lwd = 2, 
         bty = "n")
}

# ------------------------------------------------------------------------------
# EXAMPLE USAGE
# ------------------------------------------------------------------------------

#' Example: Generate simulation with refactored functions
#' 
#' @examples
#' \dontrun{
#' # Assuming you have dgm from get_dgm_stratified()
#' 
#' # Generate single simulation
#' df_sim <- draw_sim_stratified(
#'   dgm = dgm,
#'   ss = 1,
#'   wname = "meno",
#'   bw = -log(5),
#'   strata_rand = "stratum",
#'   checking = TRUE,
#'   details = TRUE,
#'   return_df = TRUE
#' )
#' 
#' # Generate population summaries
#' pop_summary <- draw_sim_stratified(
#'   dgm = dgm,
#'   ss = 1,
#'   Ndraw = 10000,  # Large sample for approximation
#'   wname = "meno",
#'   bw = -log(5),
#'   strata_rand = "stratum",
#'   return_df = FALSE
#' )
#' 
#' # Plot AHR profiles
#' plot_AHRs(pop_summary, df.case)
#' }
#' @export
example_simulation_workflow <- function() {
  message("See examples in function documentation")
}
