# =============================================================================
# GBSG-Based AFT Data Generating Mechanism (Refactored & Aligned)
# =============================================================================
#
# Refactored simulation functions for survival data based on the German Breast
# Cancer Study Group (GBSG) dataset. Aligned with forestsearch package conventions
# and generate_aft_dgm_flex() / calculate_hazard_ratios() methodology.
#
# Key functions:
#   - create_gbsg_dgm(): Create DGM with treatment effect heterogeneity
#   - simulate_from_gbsg_dgm(): Generate simulated trial data
#   - calibrate_k_inter(): Find k_inter for target subgroup HR
#   - get_dgm_with_output(): Wrapper for DGM creation with file output
#
# Alignment with generate_aft_dgm_flex():
#   - Computes theta_0, theta_1, loghr_po for each subject
#   - Uses stacked potential outcomes with same epsilon for causal HR
#   - Computes AHR, AHR_harm (AHR_H), AHR_no_harm (AHR_Hc)
#   - Returns hazard_ratios list matching generate_aft_dgm_flex output
#
# =============================================================================


# =============================================================================
# Constants
# =============================================================================

# Default seed for reproducibility
SEED_BASE <- 8316951L

# Days per month conversion
DAYS_PER_MONTH <- 30.4375

# Default super-population size
DEFAULT_N_SUPER <- 5000L


# =============================================================================
# Helper Functions
# =============================================================================

#' Discretize Continuous Variable into Quantile-Based Categories
#'
#' @param x Numeric vector to discretize
#' @param probs Numeric vector of probabilities for quantile breaks.
#'   Default: c(0.25, 0.5, 0.75) creates quartiles coded as 1, 2, 3, 4
#'
#' @return Integer vector with category codes (1 = lowest, max = highest)
#'
#' @keywords internal
cut_numeric <- function(x, probs = c(0.25, 0.5, 0.75)) {
  breaks <- stats::quantile(x, probs = probs, na.rm = TRUE)
  result <- rep(1L, length(x))
  for (i in seq_along(breaks)) {
    result[x > breaks[i]] <- i + 1L
  }
  result
}


#' Discretize Continuous Variable by Size Categories
#'
#' @param x Numeric vector (typically tumor size)
#' @param breaks Numeric vector of breakpoints. Default: c(20, 50)
#'
#' @return Integer vector with category codes
#'
#' @keywords internal
cut_size <- function(x, breaks = c(20, 50)) {
  result <- rep(1L, length(x))
  for (i in seq_along(breaks)) {
    result[x > breaks[i]] <- i + 1L
  }
  result
}


# =============================================================================
# Main DGM Creation Function
# =============================================================================

#' Create GBSG-Based AFT Data Generating Mechanism
#'
#' Creates a data generating mechanism (DGM) for survival simulations based on
#' the German Breast Cancer Study Group (GBSG) dataset. Supports heterogeneous
#' treatment effects via treatment-subgroup interactions.
#'
#' This version is aligned with \code{generate_aft_dgm_flex()} and
#' \code{calculate_hazard_ratios()} methodology, computing individual-level
#' potential outcomes and average hazard ratios (AHR).
#'
#' @param model Character. Either "alt" for alternative hypothesis with
#'   heterogeneous treatment effects, or "null" for uniform treatment effect.
#'   Default: "alt"
#' @param k_treat Numeric. Treatment effect multiplier applied to the treatment
#'   coefficient from the fitted AFT model. Values > 1 strengthen the treatment
#'   effect. Default: 1
#' @param k_inter Numeric. Interaction effect multiplier for the
#'   treatment-subgroup interaction (z1 * z3). Only used when model = "alt".
#'   Higher values create more heterogeneity between HR(H) and HR(Hc).
#'   Default: 1
#' @param k_z3 Numeric. Effect multiplier for the z3 (menopausal status)
#'   coefficient. Default: 1
#' @param z1_quantile Numeric. Quantile threshold for z1 (estrogen receptor).
#'   Observations with ER <= quantile are coded as z1 = 1. Default: 0.25
#' @param n_super Integer. Size of super-population for empirical HR estimation.
#'   Default: 5000
#' @param cens_type Character. Censoring distribution type: "weibull" or
#'   "uniform". Default: "weibull"
#' @param use_rand_params Logical. If TRUE, modifies confounder coefficients
#'   using estimates from randomized subset (meno == 0). Default: FALSE
#' @param seed Integer. Random seed for super-population generation.
#'   Default: 8316951
#' @param verbose Logical. Print diagnostic information. Default: FALSE
#'
#' @return A list of class "gbsg_dgm" containing:
#' \describe{
#'   \item{df_super_rand}{Data frame with randomized super-population including
#'     potential outcomes (theta_0, theta_1, loghr_po)}
#'   \item{hr_H_true}{Empirical hazard ratio in harm subgroup (Cox-based)}
#'   \item{hr_Hc_true}{Empirical hazard ratio in complement subgroup (Cox-based)}
#'   \item{hr_causal}{Overall causal (ITT) hazard ratio (Cox-based)}
#'   \item{AHR}{Overall average hazard ratio (from loghr_po)}
#'   \item{AHR_H_true}{Average hazard ratio in harm subgroup}
#'   \item{AHR_Hc_true}{Average hazard ratio in complement subgroup}
#'   \item{hazard_ratios}{List matching generate_aft_dgm_flex output format}
#'   \item{model_params}{List with AFT model parameters (mu, sigma, gamma, etc.)}
#'   \item{cens_params}{List with censoring model parameters}
#'   \item{subgroup_info}{List with subgroup definitions and true factor names}
#'   \item{analysis_vars}{Character vector of analysis variable names}
#'   \item{model_type}{Character indicating "alt" or "null"}
#' }
#'
#' @details
#' ## Subgroup Definition
#'
#' The harm subgroup H is defined as: z1 = 1 AND z3 = 1, where:
#' \itemize{
#'   \item z1: Low estrogen receptor (ER <= 25th percentile by default)
#'   \item z3: Premenopausal status (meno == 0)
#' }
#'
#' ## Model Specification
#'
#' The AFT model uses covariates: treat, z1, z2, z3, z4, z5, and (for "alt")
#' the interaction zh = treat * z1 * z3.
#'
#' ## Interaction Effect (k_inter)
#'
#' The k_inter parameter modifies the zh coefficient in the AFT model:
#' \preformatted{gamma[zh] <- k_inter * gamma[zh]}
#'
#' This affects the hazard ratio for the harm subgroup:
#' \itemize{
#'   \item HR(H) = exp(-gamma\[treat\]/sigma - gamma\[zh\]/sigma)
#'   \item HR(Hc) = exp(-gamma\[treat\]/sigma)
#' }
#'
#' When k_inter = 0, HR(H) = HR(Hc) (no heterogeneity).
#'
#' ## Alignment with generate_aft_dgm_flex
#'
#' This function now computes:
#' \itemize{
#'   \item theta_0: Log-hazard contribution under control
#'   \item theta_1: Log-hazard contribution under treatment
#'   \item loghr_po: Individual causal log hazard ratio (theta_1 - theta_0)
#'   \item AHR metrics: exp(mean(loghr_po)) for overall and subgroups
#' }
#'
#' @examples
#' \dontrun{
#' # Alternative hypothesis with default parameters
#' dgm_alt <- create_gbsg_dgm(model = "alt", verbose = TRUE)
#'
#' # Null hypothesis
#' dgm_null <- create_gbsg_dgm(model = "null", verbose = TRUE)
#'
#' # Custom subgroup HR via k_inter
#' dgm_custom <- create_gbsg_dgm(
#'   model = "alt",
#'   k_treat = 1.2,
#'   k_inter = 2.0,
#'   verbose = TRUE
#' )
#'
#' # Access AHR metrics (aligned with generate_aft_dgm_flex)
#' dgm_alt$hazard_ratios$AHR_harm
#' dgm_alt$hazard_ratios$AHR_no_harm
#' }
#'
#' @seealso
#' \code{\link{simulate_from_gbsg_dgm}} for generating data from the DGM
#' \code{\link{calibrate_k_inter}} for finding k_inter to achieve target HR
#'
#' @importFrom survival survreg coxph Surv
#' @importFrom stats coef quantile rexp
#' @export
create_gbsg_dgm <- function(
    model = c("alt", "null"),
    k_treat = 1,
    k_inter = 1,
    k_z3 = 1,
    z1_quantile = 0.25,
    n_super = DEFAULT_N_SUPER,
    cens_type = c("weibull", "uniform"),
    use_rand_params = FALSE,
    seed = SEED_BASE,
    verbose = FALSE
) {

  # -------------------------------------------------------------------------
  # Input Validation
  # -------------------------------------------------------------------------
  model <- match.arg(model)
  cens_type <- match.arg(cens_type)

  stopifnot(
    "k_treat must be positive" = k_treat > 0,
    "k_z3 must be numeric" = is.numeric(k_z3),
    "z1_quantile must be between 0 and 1" = z1_quantile > 0 && z1_quantile < 1,
    "n_super must be positive integer" = n_super > 0
  )

  # -------------------------------------------------------------------------
  # Load and Prepare GBSG Data
  # -------------------------------------------------------------------------
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' required for GBSG dataset.")
  }

  # Get GBSG data
  dfa <- survival::gbsg

  # Create derived variables
  dfa$y <- dfa$rfstime / DAYS_PER_MONTH
  dfa$id <- seq_len(nrow(dfa))
  dfa$event <- ifelse(dfa$status == 1, 1L, 0L)
  dfa$treat <- dfa$hormon

  # Create subgroup-defining variables
  er_threshold <- stats::quantile(dfa$er, probs = z1_quantile)
  dfa$z1 <- ifelse(dfa$er <= er_threshold, 1L, 0L)
  dfa$z2 <- cut_numeric(dfa$age)
  dfa$z3 <- ifelse(dfa$meno == 0, 1L, 0L)  # Premenopausal
  dfa$z4 <- cut_numeric(dfa$pgr)
  dfa$z5 <- cut_numeric(dfa$nodes)

  # Interaction term (treatment * subgroup indicator)
  dfa$zh <- dfa$treat * dfa$z1 * dfa$z3

  # Harm subgroup indicator (based on covariate pattern, NOT treatment)
  dfa$flag.harm <- ifelse(dfa$z1 == 1 & dfa$z3 == 1, 1L, 0L)

  # Analysis factors (observed by analyst)
  dfa$v1 <- as.factor(dfa$z1)
  dfa$v2 <- as.factor(dfa$z2)
  dfa$v3 <- as.factor(dfa$z3)
  dfa$v4 <- as.factor(dfa$z4)
  dfa$v5 <- as.factor(dfa$z5)
  dfa$v6 <- as.factor(cut_size(dfa$size))
  dfa$grade3 <- ifelse(dfa$grade == 3, 1L, 0L)
  dfa$v7 <- as.factor(dfa$grade3)

  # -------------------------------------------------------------------------
  # Define Subgroup Identities
  # -------------------------------------------------------------------------
  if (model == "alt") {
    fs_harm_true <- c("v1.1", "v3.1")
    grf_harm_true <- c("v1=1", "v3=1")
  } else {
    fs_harm_true <- NULL
    grf_harm_true <- NULL
    dfa$flag.harm <- 0L
  }

  # -------------------------------------------------------------------------
  # Verbose Output: Initial Data Summary
  # -------------------------------------------------------------------------
  if (verbose) {
    message("=== GBSG DGM Creation (Aligned) ===")
    message(sprintf("Model type: %s", model))
    message(sprintf("Effect modifiers: k_treat = %.3f, k_inter = %.3f, k_z3 = %.3f",
                    k_treat, k_inter, k_z3))
    message(sprintf("Sample size: %d", nrow(dfa)))
    message(sprintf("Events: %d (%.1f%%)", sum(dfa$event), 100 * mean(dfa$event)))
    message(sprintf("Harm subgroup size: %d (%.1f%%)",
                    sum(dfa$flag.harm), 100 * mean(dfa$flag.harm)))
    message(sprintf("ER threshold (z1_quantile = %.2f): %.1f", z1_quantile, er_threshold))
  }

  # -------------------------------------------------------------------------
  # Fit AFT Model
  # -------------------------------------------------------------------------
  covs_main <- c("treat", "z1", "z2", "z3", "z4", "z5")

  if (model == "alt") {
    covs_true <- c(covs_main, "zh")
    col_treat <- 1
    col_z1 <- 2
    col_z3 <- 4
    col_zh <- 7
  } else {
    covs_true <- covs_main
    col_treat <- 1
  }

  # Create design matrix from original data
  z_true <- as.matrix(dfa[, covs_true])

  # -------------------------------------------------------------------------
  # Create Potential Outcome Design Matrices
  # -------------------------------------------------------------------------
  # Get z1 and z3 as numeric vectors
  z1_vec <- as.numeric(dfa$z1)
  z3_vec <- as.numeric(dfa$z3)

  if (model == "alt") {
    # Potential outcome under treatment (treat = 1)
    z_true_1 <- z_true
    z_true_1[, col_treat] <- 1
    z_true_1[, col_zh] <- z1_vec * z3_vec  # zh = 1 * z1 * z3

    # Potential outcome under control (treat = 0)
    z_true_0 <- z_true
    z_true_0[, col_treat] <- 0
    z_true_0[, col_zh] <- 0  # zh = 0 * z1 * z3 = 0

  } else {
    z_true_1 <- z_true
    z_true_1[, col_treat] <- 1

    z_true_0 <- z_true
    z_true_0[, col_treat] <- 0
  }

  # Fit Weibull AFT model
  fit_aft <- survival::survreg(
    survival::Surv(y, event) ~ z_true,
    data = dfa,
    dist = "weibull"
  )
  names(fit_aft$coefficients) <- c("(Intercept)", colnames(z_true))

  # Extract parameters
  sigma <- fit_aft$scale
  mu <- stats::coef(fit_aft)[1]
  gamma <- stats::coef(fit_aft)[-1]

  if (verbose) {
    message("\n--- Original AFT Coefficients ---")
    message(sprintf("sigma (scale): %.4f", sigma))
    message(sprintf("mu (intercept): %.4f", mu))
    message("gamma coefficients:")
    for (nm in names(gamma)) {
      message(sprintf("  %s: %.4f", nm, gamma[nm]))
    }
  }

  # -------------------------------------------------------------------------
  # Optionally Modify Parameters Using Randomized Subset
  # -------------------------------------------------------------------------
  if (use_rand_params) {
    df_rand <- subset(dfa, meno == 0)
    fit_rand <- survival::survreg(
      survival::Surv(y, event) ~ treat + z1 + z2 + z3 + z4 + z5,
      data = df_rand,
      dist = "weibull"
    )
    gamma_rand <- stats::coef(fit_rand)[-1]
    gamma[c("z1", "z2", "z4", "z5")] <- gamma_rand[c("z1", "z2", "z4", "z5")]
  }

  # -------------------------------------------------------------------------
  # Apply Effect Modifiers
  # -------------------------------------------------------------------------
  gamma_orig <- gamma

  gamma["z3"] <- k_z3 * gamma["z3"]
  gamma["treat"] <- k_treat * gamma["treat"]

  if (model == "alt") {
    gamma["zh"] <- k_inter * gamma["zh"]
  }

  if (verbose) {
    message("\n--- Modified AFT Coefficients ---")
    message(sprintf("gamma['treat']: %.4f -> %.4f (k_treat = %.3f)",
                    gamma_orig["treat"], gamma["treat"], k_treat))
    message(sprintf("gamma['z3']: %.4f -> %.4f (k_z3 = %.3f)",
                    gamma_orig["z3"], gamma["z3"], k_z3))
    if (model == "alt") {
      message(sprintf("gamma['zh']: %.4f -> %.4f (k_inter = %.3f)",
                      gamma_orig["zh"], gamma["zh"], k_inter))
    }
  }

  # -------------------------------------------------------------------------
  # Convert to Hazard Scale Coefficients
  # -------------------------------------------------------------------------
  # In AFT: log(T) = mu + X*gamma + sigma*epsilon
  # Hazard: h(t|X) = h_0(t) * exp(X * b0) where b0 = -gamma/sigma
  b_true <- gamma
  b0 <- -gamma / sigma

  if (verbose && model == "alt") {
    hr_H_theoretical <- exp(b0["treat"] + b0["zh"])
    hr_Hc_theoretical <- exp(b0["treat"])
    message("\n--- Theoretical HRs (from coefficients) ---")
    message(sprintf("HR(H) = exp(b0['treat'] + b0['zh']) = exp(%.4f + %.4f) = %.4f",
                    b0["treat"], b0["zh"], hr_H_theoretical))
    message(sprintf("HR(Hc) = exp(b0['treat']) = exp(%.4f) = %.4f",
                    b0["treat"], hr_Hc_theoretical))
  }

  # -------------------------------------------------------------------------
  # Compute Linear Predictors for Potential Outcomes (AFT scale)
  # -------------------------------------------------------------------------
  lin_conf <- as.vector(z_true %*% b_true)
  lin1_conf <- as.vector(z_true_1 %*% b_true)
  lin0_conf <- as.vector(z_true_0 %*% b_true)

  # -------------------------------------------------------------------------
  # Compute Individual-Level Potential Outcome Log-Hazards (ALIGNED)
  # -------------------------------------------------------------------------
  # theta_0: log-hazard contribution under control (treat=0)
  # theta_1: log-hazard contribution under treatment (treat=1)
  # These use HAZARD-SCALE coefficients (b0), matching generate_aft_dgm_flex

  dfa$theta_0 <- as.vector(z_true_0 %*% b0)
  dfa$theta_1 <- as.vector(z_true_1 %*% b0)
  dfa$loghr_po <- dfa$theta_1 - dfa$theta_0

  # Also store linear predictors (AFT scale) for simulation
  dfa$lin_pred_0 <- lin0_conf
  dfa$lin_pred_1 <- lin1_conf

  # -------------------------------------------------------------------------
  # Generate Super-Population Sample
  # -------------------------------------------------------------------------
  set.seed(seed)

  # Sample with replacement from original data
  id_sample <- sample(seq_len(nrow(dfa)), size = n_super, replace = TRUE)
  df_samp <- dfa[id_sample, ]
  df_samp$id <- seq_len(n_super)

  # Carry forward linear predictors and potential outcomes
  df_samp$lin1.conf <- lin1_conf[id_sample]
  df_samp$lin0.conf <- lin0_conf[id_sample]
  df_samp$theta_0 <- dfa$theta_0[id_sample]
  df_samp$theta_1 <- dfa$theta_1[id_sample]
  df_samp$loghr_po <- dfa$loghr_po[id_sample]

  # -------------------------------------------------------------------------
  # Compute Hazard Ratios Using STACKED Potential Outcomes (ALIGNED)
  # -------------------------------------------------------------------------
  # This matches the methodology in calculate_hazard_ratios()
  # Use SAME epsilon for both potential outcomes (causal framework)

  epsilon <- log(stats::rexp(n_super))

  # Potential survival time under treatment
  logT_1 <- mu + sigma * epsilon + df_samp$lin1.conf
  T_1 <- exp(logT_1)

  # Potential survival time under control
  logT_0 <- mu + sigma * epsilon + df_samp$lin0.conf
  T_0 <- exp(logT_0)

  # Stack for Cox model (2 * n_super rows)
  df_po <- data.frame(
    time = c(T_1, T_0),
    event = 1L,
    treat = c(rep(1L, n_super), rep(0L, n_super)),
    flag.harm = rep(df_samp$flag.harm, 2)
  )

  # Cox-based HRs from stacked potential outcomes
  hr_causal <- exp(survival::coxph(
    survival::Surv(time, event) ~ treat,
    data = df_po
  )$coefficients)

  if (model == "alt") {
    df_po_H <- subset(df_po, flag.harm == 1)
    df_po_Hc <- subset(df_po, flag.harm == 0)

    if (nrow(df_po_H) > 20 && length(unique(df_po_H$treat)) == 2) {
      hr_H_true <- exp(survival::coxph(
        survival::Surv(time, event) ~ treat,
        data = df_po_H
      )$coefficients)
    } else {
      warning("Insufficient data in harm subgroup for HR calculation")
      hr_H_true <- NA_real_
    }

    if (nrow(df_po_Hc) > 20 && length(unique(df_po_Hc$treat)) == 2) {
      hr_Hc_true <- exp(survival::coxph(
        survival::Surv(time, event) ~ treat,
        data = df_po_Hc
      )$coefficients)
    } else {
      warning("Insufficient data in complement subgroup for HR calculation")
      hr_Hc_true <- NA_real_
    }

  } else {
    hr_H_true <- NA_real_
    hr_Hc_true <- hr_causal
  }

  # -------------------------------------------------------------------------
  # Compute Average Hazard Ratios (AHR) from loghr_po (ALIGNED)
  # -------------------------------------------------------------------------
  # This matches the AHR calculation in calculate_hazard_ratios()

  AHR <- exp(mean(df_samp$loghr_po))

  if (model == "alt" && sum(df_samp$flag.harm) > 0) {
    AHR_H_true <- exp(mean(df_samp$loghr_po[df_samp$flag.harm == 1]))
    AHR_Hc_true <- exp(mean(df_samp$loghr_po[df_samp$flag.harm == 0]))
  } else {
    AHR_H_true <- NA_real_
    AHR_Hc_true <- AHR
  }

  if (verbose) {
    message("\n--- Hazard Ratios (Stacked Potential Outcomes) ---")
    message(sprintf("HR (overall/causal): %.4f", hr_causal))
    message(sprintf("HR (harm subgroup, H): %.4f", hr_H_true))
    message(sprintf("HR (complement, Hc): %.4f", hr_Hc_true))

    message("\n--- Average Hazard Ratios (from loghr_po) ---")
    message(sprintf("AHR (overall): %.4f", AHR))
    message(sprintf("AHR (harm subgroup, H): %.4f", AHR_H_true))
    message(sprintf("AHR (complement, Hc): %.4f", AHR_Hc_true))

    if (model == "alt" && !is.na(hr_H_true) && !is.na(hr_Hc_true)) {
      message(sprintf("\nHR(H) / HR(Hc) ratio: %.4f", hr_H_true / hr_Hc_true))
      message(sprintf("AHR(H) / AHR(Hc) ratio: %.4f", AHR_H_true / AHR_Hc_true))
    }
  }

  # -------------------------------------------------------------------------
  # Create Super-Population with Assigned Treatment (for simulation)
  # -------------------------------------------------------------------------
  n_treat <- round(n_super / 2)
  n_control <- n_super - n_treat

  # Treatment group
  df_treat <- df_samp[seq_len(n_treat), ]
  df_treat$treat <- 1L
  df_treat$lin.conf.true <- df_treat$lin1.conf
  df_treat$zh <- df_treat$z1 * df_treat$z3

  # Control group
  df_control <- df_samp[(n_treat + 1):n_super, ]
  df_control$treat <- 0L
  df_control$lin.conf.true <- df_control$lin0.conf
  df_control$zh <- 0L

  df_big <- rbind(df_treat, df_control)

  # Generate observed survival times for the assigned-treatment super-population
  set.seed(seed + 1)
  epsilon_obs <- log(stats::rexp(n_super))
  log_Ts_obs <- mu + sigma * epsilon_obs + df_big$lin.conf.true
  df_big$Ts <- exp(log_Ts_obs)
  df_big$es <- 1L

  # Add potential outcome hazard ratios (for compatibility)
  df_big$hlin.conf.1 <- exp(df_big$theta_1)
  df_big$hlin.conf.0 <- exp(df_big$theta_0)
  df_big$hlin.ratio <- df_big$hlin.conf.1 / df_big$hlin.conf.0
  df_big$h1.potential <- df_big$hlin.conf.1
  df_big$h0.potential <- df_big$hlin.conf.0

  # -------------------------------------------------------------------------
  # Fit Censoring Model
  # -------------------------------------------------------------------------
  cens_model <- NULL
  mu_cens <- NULL
  sigma_cens <- NULL
  gamma_cens <- NULL

  if (cens_type == "weibull") {
    z_cens <- as.matrix(df_big[, covs_main])

    fit_cens <- survival::survreg(
      survival::Surv(y, 1 - event) ~ z_cens,
      data = df_big,
      dist = "weibull"
    )

    sigma_cens <- fit_cens$scale
    mu_cens <- stats::coef(fit_cens)[1]
    gamma_cens <- stats::coef(fit_cens)[-1]
    b0_cens <- -gamma_cens / sigma_cens
    b_cens <- -b0_cens * sigma

    # Censoring linear predictors for potential outcomes
    z_cens_1 <- z_cens
    z_cens_1[, 1] <- 1
    z_cens_0 <- z_cens
    z_cens_0[, 1] <- 0

    df_big$linC1.conf <- as.vector(z_cens_1 %*% b_cens)
    df_big$linC0.conf <- as.vector(z_cens_0 %*% b_cens)
    df_big$linC.conf <- ifelse(df_big$treat == 1,
                                df_big$linC1.conf,
                                df_big$linC0.conf)

    cens_model <- list(
      type = "weibull",
      mu = mu_cens,
      sigma = sigma_cens,
      gamma = gamma_cens
    )
  }

  # -------------------------------------------------------------------------
  # Assemble Output (ALIGNED with generate_aft_dgm_flex)
  # -------------------------------------------------------------------------
  analysis_vars <- c("v1", "v2", "v3", "v4", "v5", "v6", "v7")

  # Hazard ratios list matching generate_aft_dgm_flex output
  hazard_ratios <- list(
    overall = hr_causal,
    AHR = AHR,
    AHR_harm = AHR_H_true,
    AHR_no_harm = AHR_Hc_true,
    harm_subgroup = hr_H_true,
    no_harm_subgroup = hr_Hc_true
  )

  result <- list(
    # Super-population data
    df_super_rand = df_big,

    # Cox-based HRs (backward compatible)
    hr_H_true = hr_H_true,
    hr_Hc_true = hr_Hc_true,
    hr_causal = hr_causal,

    # AHR metrics (aligned with generate_aft_dgm_flex)
    AHR = AHR,
    AHR_H_true = AHR_H_true,
    AHR_Hc_true = AHR_Hc_true,

    # Combined hazard_ratios list (matches generate_aft_dgm_flex)
    hazard_ratios = hazard_ratios,

    # Model parameters
    model_params = list(
      mu = mu,
      sigma = sigma,
      gamma = gamma,
      gamma_orig = gamma_orig,
      b_true = b_true,
      b_hr = b0
    ),

    # Censoring parameters
    cens_params = list(
      type = cens_type,
      mu = mu_cens,
      sigma = sigma_cens,
      model = cens_model
    ),

    # Subgroup information
    subgroup_info = list(
      fs_harm_true = fs_harm_true,
      grf_harm_true = grf_harm_true,
      z1_quantile = z1_quantile,
      er_threshold = er_threshold,
      definition = "z1 == 1 & z3 == 1 (low ER & premenopausal)",
      size = sum(df_big$flag.harm),
      proportion = mean(df_big$flag.harm)
    ),

    # Effect modifiers
    effect_modifiers = list(
      k_treat = k_treat,
      k_inter = k_inter,
      k_z3 = k_z3
    ),

    # Analysis variables
    analysis_vars = analysis_vars,
    model_type = model,
    n_super = n_super,
    seed = seed
  )

  class(result) <- c("gbsg_dgm", "list")

  return(result)
}


# =============================================================================
# Simulation Function
# =============================================================================

#' Simulate Trial Data from GBSG DGM
#'
#' Generates simulated clinical trial data from a GBSG-based data generating
#' mechanism.
#'
#' @param dgm A "gbsg_dgm" object from \code{\link{create_gbsg_dgm}}
#' @param n Integer. Sample size. If NULL, uses full super-population.
#'   Default: NULL
#' @param rand_ratio Numeric. Randomization ratio (treatment:control).
#'   Default: 1 (1:1 randomization)
#' @param sim_id Integer. Simulation ID used for seed offset. Default: 1
#' @param max_follow Numeric. Administrative censoring time (months).
#'   Default: Inf (no administrative censoring)
#' @param muC_adj Numeric. Adjustment to censoring distribution location
#'   parameter. Positive values increase censoring. Default: 0
#' @param min_cens Numeric. Minimum censoring time for uniform censoring.
#'   Required if cens_type = "uniform"
#' @param max_cens Numeric. Maximum censoring time for uniform censoring.
#'   Required if cens_type = "uniform"
#' @param draw_treatment Logical. If TRUE, randomly assigns treatment.
#'   If FALSE, samples from existing treatment arms. Default: TRUE
#'
#' @return Data frame with simulated trial data including:
#' \describe{
#'   \item{id}{Subject identifier}
#'   \item{y.sim}{Observed follow-up time}
#'   \item{event.sim}{Event indicator (1 = event, 0 = censored)}
#'   \item{t.sim}{True event time (before censoring)}
#'   \item{treat}{Treatment indicator}
#'   \item{flag.harm}{Harm subgroup indicator}
#'   \item{loghr_po}{Individual log hazard ratio (potential outcome)}
#'   \item{v1-v7}{Analysis factors}
#' }
#'
#' @examples
#' \dontrun{
#' dgm <- create_gbsg_dgm(model = "alt", k_inter = 2.0)
#' sim_data <- simulate_from_gbsg_dgm(dgm, n = 500, sim_id = 1)
#'
#' # Check AHR in simulated data
#' exp(mean(sim_data$loghr_po))
#' }
#'
#' @export
simulate_from_gbsg_dgm <- function(
    dgm,
    n = NULL,
    rand_ratio = 1,
    sim_id = 1,
    max_follow = Inf,
    muC_adj = 0,
    min_cens = NULL,
    max_cens = NULL,
    draw_treatment = TRUE
) {

  stopifnot(
    "dgm must be a gbsg_dgm object" = inherits(dgm, "gbsg_dgm"),
    "rand_ratio must be positive" = rand_ratio > 0
  )

  cens_type <- dgm$cens_params$type

  # -------------------------------------------------------------------------
  # Sample from Super-Population
  # -------------------------------------------------------------------------
  df_super <- dgm$df_super_rand

  if (is.null(n)) {
    df_sim <- df_super
  } else {
    # Calculate group sizes
    n0 <- round(n / (1 + rand_ratio))
    n1 <- n - n0

    if (draw_treatment) {
      # Random sampling with treatment assignment
      set.seed(dgm$seed + sim_id)
      id_sample <- sample(seq_len(nrow(df_super)), size = n, replace = TRUE)
      df_sample <- df_super[id_sample, ]

      # Assign treatment
      df_treat <- df_sample[seq_len(n1), ]
      df_treat$treat <- 1L
      df_treat$lin.conf.true <- df_treat$lin1.conf
      df_treat$zh <- df_treat$z1 * df_treat$z3
      if ("linC1.conf" %in% names(df_treat)) {
        df_treat$linC.conf <- df_treat$linC1.conf
      }

      df_control <- df_sample[(n1 + 1):n, ]
      df_control$treat <- 0L
      df_control$lin.conf.true <- df_control$lin0.conf
      df_control$zh <- 0L
      if ("linC0.conf" %in% names(df_control)) {
        df_control$linC.conf <- df_control$linC0.conf
      }

      df_sim <- rbind(df_control, df_treat)

    } else {
      # Sample from existing treatment arms
      df1_super <- subset(df_super, treat == 1)
      df0_super <- subset(df_super, treat == 0)

      set.seed(dgm$seed + sim_id)
      id1 <- sample(seq_len(nrow(df1_super)), size = n1, replace = TRUE)
      id0 <- sample(seq_len(nrow(df0_super)), size = n0, replace = TRUE)

      df_sim <- rbind(df0_super[id0, ], df1_super[id1, ])
    }
  }

  n_obs <- nrow(df_sim)

  # -------------------------------------------------------------------------
  # Generate Event Times
  # -------------------------------------------------------------------------
  mu <- dgm$model_params$mu
  sigma <- dgm$model_params$sigma

  set.seed(dgm$seed + sim_id + 1000)
  epsilon <- log(stats::rexp(n_obs))
  log_T <- mu + sigma * epsilon + df_sim$lin.conf.true
  T_sim <- exp(log_T)

  # -------------------------------------------------------------------------
  # Generate Censoring Times
  # -------------------------------------------------------------------------
  set.seed(dgm$seed + sim_id + 2000)

  if (cens_type == "weibull" && !is.null(dgm$cens_params$mu)) {
    mu_cens <- dgm$cens_params$mu + muC_adj
    sigma_cens <- dgm$cens_params$sigma

    epsilon_cens <- log(stats::rexp(n_obs))
    log_C <- mu_cens + sigma_cens * epsilon_cens + df_sim$linC.conf
    C_sim <- exp(log_C)

  } else {
    if (is.null(min_cens) || is.null(max_cens)) {
      min_cens <- 0
      max_cens <- max(T_sim) * 1.5
    }
    C_sim <- stats::runif(n_obs, min = min_cens, max = max_cens)
  }

  # Apply administrative censoring
  C_sim <- pmin(C_sim, max_follow)

  # -------------------------------------------------------------------------
  # Create Observed Outcomes
  # -------------------------------------------------------------------------
  event_sim <- ifelse(T_sim <= C_sim, 1L, 0L)
  y_sim <- pmin(T_sim, C_sim)

  # -------------------------------------------------------------------------
  # Finalize Output
  # -------------------------------------------------------------------------
  df_sim$y.sim <- y_sim
  df_sim$event.sim <- event_sim
  df_sim$t.sim <- T_sim
  df_sim$id <- seq_len(n_obs)

  # Select output columns (including new aligned variables)
  keep_cols <- c(
    "id", "y.sim", "event.sim", "t.sim", "treat", "flag.harm",
    "v1", "v2", "v3", "v4", "v5", "v6", "v7",
    "z1", "z2", "z3", "z4", "z5", "zh", "grade3", "size",
    "theta_0", "theta_1", "loghr_po",
    "hlin.ratio", "h1.potential", "h0.potential"
  )
  keep_cols <- intersect(keep_cols, names(df_sim))

  as.data.frame(df_sim[, keep_cols])
}


# =============================================================================
# Calibration Functions
# =============================================================================

#' Calibrate k_inter for Target Subgroup Hazard Ratio
#'
#' Finds the interaction effect multiplier (k_inter) that achieves a target
#' hazard ratio in the harm subgroup.
#'
#' @param target_hr_harm Numeric. Target hazard ratio for the harm subgroup
#' @param model Character. Model type ("alt" only). Default: "alt"
#' @param k_treat Numeric. Treatment effect multiplier. Default: 1
#' @param cens_type Character. Censoring type. Default: "weibull"
#' @param k_inter_range Numeric vector of length 2. Search range for k_inter.
#'   Default: c(-100, 100)
#' @param tol Numeric. Tolerance for root finding. Default: 1e-6
#' @param use_ahr Logical. If TRUE, calibrate to AHR instead of Cox-based HR.
#'   Default: FALSE
#' @param verbose Logical. Print diagnostic information. Default: FALSE
#' @param ... Additional arguments passed to \code{\link{create_gbsg_dgm}}
#'
#' @return Numeric value of k_inter that achieves the target HR
#'
#' @details
#' This function uses \code{uniroot} to find the k_inter value such that
#' the empirical HR (or AHR) in the harm subgroup equals target_hr_harm.
#'
#' @examples
#' \dontrun{
#' # Find k_inter for HR = 1.5 in harm subgroup
#' k <- calibrate_k_inter(target_hr_harm = 1.5, verbose = TRUE)
#'
#' # Verify
#' dgm <- create_gbsg_dgm(model = "alt", k_inter = k, verbose = TRUE)
#' print(dgm$hr_H_true)  # Should be close to 1.5
#'
#' # Calibrate to AHR instead
#' k_ahr <- calibrate_k_inter(target_hr_harm = 1.5, use_ahr = TRUE, verbose = TRUE)
#' dgm_ahr <- create_gbsg_dgm(model = "alt", k_inter = k_ahr, verbose = TRUE)
#' print(dgm_ahr$AHR_H_true)  # Should be close to 1.5
#' }
#'
#' @importFrom stats uniroot
#' @export
calibrate_k_inter <- function(
    target_hr_harm,
    model = "alt",
    k_treat = 1,
    cens_type = "weibull",
    k_inter_range = c(-100, 100),
    tol = 1e-6,
    use_ahr = FALSE,
    verbose = FALSE,
    ...
) {

  stopifnot(
    "target_hr_harm must be positive" = target_hr_harm > 0,
    "model must be 'alt' for calibration" = model == "alt"
  )

  # Objective function
  objective <- function(k_val) {
    dgm <- create_gbsg_dgm(
      model = model,
      k_treat = k_treat,
      k_inter = k_val,
      cens_type = cens_type,
      verbose = FALSE,
      ...
    )
    if (use_ahr) {
      dgm$AHR_H_true - target_hr_harm
    } else {
      dgm$hr_H_true - target_hr_harm
    }
  }

  hr_type <- if (use_ahr) "AHR(H)" else "HR(H)"

  if (verbose) {
    message(sprintf("Calibrating k_inter to achieve %s = %.4f", hr_type, target_hr_harm))
    message(sprintf("Search range: [%.1f, %.1f]", k_inter_range[1], k_inter_range[2]))

    # Show HR at boundaries
    dgm_lower <- create_gbsg_dgm(model = model, k_treat = k_treat,
                                  k_inter = k_inter_range[1], verbose = FALSE, ...)
    dgm_upper <- create_gbsg_dgm(model = model, k_treat = k_treat,
                                  k_inter = k_inter_range[2], verbose = FALSE, ...)

    if (use_ahr) {
      message(sprintf("%s at k_inter = %.1f: %.4f", hr_type, k_inter_range[1], dgm_lower$AHR_H_true))
      message(sprintf("%s at k_inter = %.1f: %.4f", hr_type, k_inter_range[2], dgm_upper$AHR_H_true))
    } else {
      message(sprintf("%s at k_inter = %.1f: %.4f", hr_type, k_inter_range[1], dgm_lower$hr_H_true))
      message(sprintf("%s at k_inter = %.1f: %.4f", hr_type, k_inter_range[2], dgm_upper$hr_H_true))
    }
  }

  result <- tryCatch({
    stats::uniroot(objective, interval = k_inter_range, tol = tol)
  }, error = function(e) {
    warning("Root finding failed: ", e$message)
    return(NULL)
  })

  if (is.null(result)) {
    return(NA_real_)
  }

  k_inter <- result$root

  if (verbose) {
    dgm_verify <- create_gbsg_dgm(
      model = model,
      k_treat = k_treat,
      k_inter = k_inter,
      cens_type = cens_type,
      verbose = FALSE,
      ...
    )
    achieved <- if (use_ahr) dgm_verify$AHR_H_true else dgm_verify$hr_H_true
    message(sprintf("\n=== Calibration Result ==="))
    message(sprintf("Found k_inter = %.6f", k_inter))
    message(sprintf("Achieved %s = %.4f (target: %.4f)", hr_type, achieved, target_hr_harm))
    message(sprintf("Error: %.6f", abs(achieved - target_hr_harm)))
    message(sprintf("Iterations: %d", result$iter))
  }

  k_inter
}


#' Create DGM with Output File Path
#'
#' Wrapper function that creates a GBSG DGM and generates a standardized
#' output file path for saving results.
#'
#' @param model_harm Character. Model type ("alt" or "null")
#' @param n Integer. Planned sample size (for filename)
#' @param k_treat Numeric. Treatment effect multiplier
#' @param target_hr_harm Numeric. Target HR for harm subgroup (used for
#'   calibration when model = "alt")
#' @param cens_type Character. Censoring type
#' @param out_dir Character. Output directory path. If NULL, no file path
#'   is generated
#' @param file_prefix Character. Prefix for output filename
#' @param file_suffix Character. Suffix for output filename
#' @param include_hr_in_name Logical. Include achieved HR in filename.
#'   Default: FALSE
#' @param verbose Logical. Print diagnostic information. Default: FALSE
#' @param ... Additional arguments passed to \code{\link{create_gbsg_dgm}}
#'
#' @return List with components:
#' \describe{
#'   \item{dgm}{The gbsg_dgm object}
#'   \item{out_file}{Character path to output file (NULL if out_dir is NULL)}
#'   \item{k_inter}{The k_inter value used (either calibrated or default)}
#' }
#'
#' @keywords internal
get_dgm_with_output <- function(
    model_harm,
    n,
    k_treat = 1,
    target_hr_harm = NULL,
    cens_type = "weibull",
    out_dir = NULL,
    file_prefix = "sim",
    file_suffix = "",
    include_hr_in_name = FALSE,
    verbose = FALSE,
    ...
) {

  # Calibrate k_inter if needed
  k_inter <- 1
  if (model_harm != "null" && !is.null(target_hr_harm)) {
    k_inter <- calibrate_k_inter(
      target_hr_harm = target_hr_harm,
      model = model_harm,
      k_treat = k_treat,
      cens_type = cens_type,
      verbose = verbose,
      ...
    )
    if (is.na(k_inter)) {
      warning("Calibration failed, using k_inter = 1")
      k_inter <- 1
    }
  }

  # Create DGM
  dgm <- create_gbsg_dgm(
    model = model_harm,
    k_treat = k_treat,
    k_inter = k_inter,
    cens_type = cens_type,
    verbose = verbose,
    ...
  )

  # Generate output file path
  out_file <- NULL
  if (!is.null(out_dir)) {
    out_file <- file.path(
      out_dir,
      paste0(
        file_prefix, "_",
        "N=", n, "_",
        model_harm, "_",
        "ktreat=", k_treat,
        if (include_hr_in_name && model_harm == "alt") {
          paste0("_hrH=", round(dgm$hr_H_true, 2))
        } else "",
        file_suffix,
        ".Rdata"
      )
    )
  }

  list(dgm = dgm, out_file = out_file, k_inter = k_inter)
}


# =============================================================================
# Print Method
# =============================================================================

#' Print Method for gbsg_dgm Objects
#'
#' @param x A gbsg_dgm object
#' @param ... Additional arguments (unused)
#'
#' @export
print.gbsg_dgm <- function(x, ...) {
  cat("GBSG-Based AFT Data Generating Mechanism (Aligned)\n")
  cat("===================================================\n\n")

  cat("Model type:", x$model_type, "\n")
  cat("Super-population size:", x$n_super, "\n\n")

  cat("Effect Modifiers:\n")
  cat("  k_treat:", x$effect_modifiers$k_treat, "\n")
  cat("  k_inter:", x$effect_modifiers$k_inter, "\n")
  cat("  k_z3:", x$effect_modifiers$k_z3, "\n\n")

  cat("Hazard Ratios (Cox-based, stacked PO):\n")
  cat("  Overall (causal):", round(x$hr_causal, 4), "\n")
  if (!is.na(x$hr_H_true)) {
    cat("  Harm subgroup (H):", round(x$hr_H_true, 4), "\n")
  }
  cat("  Complement (Hc):", round(x$hr_Hc_true, 4), "\n")
  if (!is.na(x$hr_H_true) && !is.na(x$hr_Hc_true)) {
    cat("  Ratio HR(H)/HR(Hc):", round(x$hr_H_true / x$hr_Hc_true, 4), "\n")
  }
  cat("\n")

  cat("Average Hazard Ratios (from loghr_po):\n")
  cat("  AHR (overall):", round(x$AHR, 4), "\n")
  if (!is.na(x$AHR_H_true)) {
    cat("  AHR_harm (H):", round(x$AHR_H_true, 4), "\n")
  }
  if (!is.na(x$AHR_Hc_true)) {
    cat("  AHR_no_harm (Hc):", round(x$AHR_Hc_true, 4), "\n")
  }
  if (!is.na(x$AHR_H_true) && !is.na(x$AHR_Hc_true)) {
    cat("  Ratio AHR(H)/AHR(Hc):", round(x$AHR_H_true / x$AHR_Hc_true, 4), "\n")
  }
  cat("\n")

  cat("Subgroup definition:", x$subgroup_info$definition, "\n")
  cat("ER threshold:", round(x$subgroup_info$er_threshold, 2),
      sprintf("(quantile = %.2f)\n", x$subgroup_info$z1_quantile))
  cat("Subgroup size:", x$subgroup_info$size,
      sprintf("(%.1f%%)\n", 100 * x$subgroup_info$proportion))
  cat("Analysis variables:", paste(x$analysis_vars, collapse = ", "), "\n")

  invisible(x)
}


# =============================================================================
# Validation Function
# =============================================================================

#' Validate k_inter Effect on HR Heterogeneity
#'
#' Test function to verify that k_inter properly modulates the difference
#' between HR(H) and HR(Hc), and that AHR metrics align with Cox-based HRs.
#'
#' @param k_inter_values Numeric vector of k_inter values to test.
#'   Default: c(-2, -1, 0, 1, 2, 3)
#' @param verbose Logical. Print results. Default: TRUE
#' @param ... Additional arguments passed to create_gbsg_dgm
#'
#' @return Data frame with k_inter, hr_H, hr_Hc, AHR_H, AHR_Hc, and ratio columns
#'
#' @examples
#' \dontrun{
#' # Test k_inter effect
#' results <- validate_k_inter_effect()
#'
#' # k_inter = 0 should give hr_H approximately equals hr_Hc (ratio approximately 1)
#' }
#'
#' @export
validate_k_inter_effect <- function(
    k_inter_values = c(-2, -1, 0, 1, 2, 3),
    verbose = TRUE,
    ...
) {

  results <- data.frame(
    k_inter = numeric(),
    hr_H = numeric(),
    hr_Hc = numeric(),
    AHR_H = numeric(),
    AHR_Hc = numeric(),
    ratio_cox = numeric(),
    ratio_ahr = numeric()
  )

  if (verbose) {
    cat("Testing k_inter effect on HR heterogeneity...\n\n")
    cat(sprintf("%-8s %-8s %-8s %-8s %-8s %-10s %-10s\n",
                "k_inter", "HR(H)", "HR(Hc)", "AHR(H)", "AHR(Hc)",
                "Ratio(Cox)", "Ratio(AHR)"))
    cat(paste(rep("-", 70), collapse = ""), "\n")
  }

  for (k in k_inter_values) {
    dgm <- create_gbsg_dgm(model = "alt", k_inter = k, verbose = FALSE, ...)

    ratio_cox <- if (!is.na(dgm$hr_H_true) && !is.na(dgm$hr_Hc_true)) {
      dgm$hr_H_true / dgm$hr_Hc_true
    } else NA

    ratio_ahr <- if (!is.na(dgm$AHR_H_true) && !is.na(dgm$AHR_Hc_true)) {
      dgm$AHR_H_true / dgm$AHR_Hc_true
    } else NA

    results <- rbind(results, data.frame(
      k_inter = k,
      hr_H = dgm$hr_H_true,
      hr_Hc = dgm$hr_Hc_true,
      AHR_H = dgm$AHR_H_true,
      AHR_Hc = dgm$AHR_Hc_true,
      ratio_cox = ratio_cox,
      ratio_ahr = ratio_ahr
    ))

    if (verbose) {
      cat(sprintf("%-8.1f %-8.4f %-8.4f %-8.4f %-8.4f %-10.4f %-10.4f\n",
                  k, dgm$hr_H_true, dgm$hr_Hc_true,
                  dgm$AHR_H_true, dgm$AHR_Hc_true,
                  ratio_cox, ratio_ahr))
    }
  }

  if (verbose) {
    cat("\n")
    # Check if k_inter = 0 gives ratio approximately 1
    k0_row <- results[results$k_inter == 0, ]
    if (nrow(k0_row) > 0) {
      if (abs(k0_row$ratio_cox - 1) < 0.05) {
        cat("PASS: k_inter = 0 gives Cox ratio ~= 1 (no heterogeneity)\n")
      } else {
        cat("CHECK: k_inter = 0 gives Cox ratio =", round(k0_row$ratio_cox, 4),
            "- expected ~= 1\n")
      }
      if (abs(k0_row$ratio_ahr - 1) < 0.05) {
        cat("PASS: k_inter = 0 gives AHR ratio ~= 1 (no heterogeneity)\n")
      } else {
        cat("CHECK: k_inter = 0 gives AHR ratio =", round(k0_row$ratio_ahr, 4),
            "- expected ~= 1\n")
      }
    }

    # Check AHR vs Cox alignment
    cat("\nAHR vs Cox HR alignment:\n")
    for (i in seq_len(nrow(results))) {
      hr_diff <- abs(results$hr_H[i] - results$AHR_H[i])
      if (hr_diff < 0.02) {
        cat(sprintf("  k_inter = %.1f: HR(H) vs AHR(H) diff = %.4f (ALIGNED)\n",
                    results$k_inter[i], hr_diff))
      } else {
        cat(sprintf("  k_inter = %.1f: HR(H) vs AHR(H) diff = %.4f\n",
                    results$k_inter[i], hr_diff))
      }
    }
  }

  invisible(results)
}
