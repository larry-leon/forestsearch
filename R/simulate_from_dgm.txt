#' Simulate Survival Data from AFT Data Generating Mechanism
#'
#' Generates simulated survival data from a previously created AFT data generating
#' mechanism (DGM). Samples from the super population and generates survival times
#' with specified censoring.
#'
#' @param dgm An object of class "aft_dgm_flex" created by
#'   \code{\link{generate_aft_dgm_flex}}
#' @param entry_var Character string specifying the name of an entry time variable
#'   in the data. If NULL, entry times are simulated uniformly. Default NULL
#' @param max_entry Numeric specifying maximum entry time for staggered entry.
#'   Entry times are simulated as Uniform(0, max_entry). Default 24
#' @param analysis_time Numeric specifying the calendar time at which analysis
#'   occurs. Follow-up time is calculated as analysis_time - entry_time. Default 48
#' @param n Integer specifying the sample size. If NULL (default), uses the
#'   entire super population
#' @param rand_ratio Numeric specifying the randomization ratio (treatment:control).
#'   Default is 1 (1:1 allocation)
#' @param cens_adjust Numeric adjustment to censoring distribution on log scale.
#'   Positive values increase censoring, negative values decrease it. Default is 0
#' @param draw_treatment Logical indicating whether to redraw treatment assignment.
#'   If TRUE (default), reassigns treatment according to rand_ratio. If FALSE,
#'   keeps original treatment assignments from super population
#' @param seed Integer random seed for reproducibility. Default is NULL (no seed set)
#'
#' @return A data.frame containing simulated survival data with columns:
#' \describe{
#'   \item{id}{Subject identifier}
#'   \item{treat}{Treatment assignment (0 or 1)}
#'   \item{treat_sim}{Simulated treatment assignment (may differ from treat if
#'     draw_treatment = TRUE)}
#'   \item{flag_harm}{Subgroup indicator (1 if all subgroup conditions met, 0 otherwise)}
#'   \item{z_*}{Standardized covariate values}
#'   \item{y_sim}{Observed survival time (minimum of true time and censoring time)}
#'   \item{event_sim}{Event indicator (1 = event observed, 0 = censored)}
#'   \item{t_true}{True underlying survival time (before censoring)}
#'   \item{c_time}{Censoring time}
#' }
#'
#' @details
#' ## Simulation Process
#'
#' 1. **Sampling**: Draws n observations with replacement from the super population
#' 2. **Treatment Assignment**:
#'    - If `draw_treatment = TRUE`: Reassigns treatment based on `rand_ratio`
#'    - If `draw_treatment = FALSE`: Keeps original treatment assignments
#' 3. **Survival Times**: Generates from Weibull AFT model:
#'    \deqn{\log(T) = \mu + \sigma \epsilon + X'\gamma}
#'    where \eqn{\epsilon} ~ extreme value distribution
#' 4. **Censoring**: Applies specified censoring distribution (Weibull or uniform)
#' 5. **Administrative Censoring**: Applies max_follow cutoff if specified
#'
#' ## Censoring Adjustment
#'
#' The `cens_adjust` parameter modifies the censoring distribution:
#' - `cens_adjust = log(2)` doubles expected censoring times
#' - `cens_adjust = log(0.5)` halves expected censoring times
#'
#' @examples
#' \dontrun{
#' # Create DGM first
#' dgm <- generate_aft_dgm_flex(
#'   data = gbsg,
#'   continuous_vars = c("age", "size", "nodes", "pgr", "er"),
#'   factor_vars = c("meno", "grade"),
#'   outcome_var = "rfstime",
#'   event_var = "status",
#'   treatment_var = "hormon",
#'   subgroup_vars = c("er", "meno"),
#'   subgroup_cuts = list(er = 20, meno = 0),
#'   model = "alt"
#' )
#'
#' # Simulate data with 1:1 randomization
#' sim_data <- simulate_from_dgm(
#'   dgm = dgm,
#'   n = 1000,
#'   rand_ratio = 1,
#'   max_follow = 84,
#'   cens_adjust = log(1.5),
#'   seed = 123
#' )
#'
#' # Check results
#' table(sim_data$treat_sim)
#' mean(sim_data$event_sim)
#'
#' # Simulate with 2:1 randomization
#' sim_data_2to1 <- simulate_from_dgm(
#'   dgm = dgm,
#'   n = 900,
#'   rand_ratio = 2,  # 2:1 treatment:control
#'   seed = 456
#' )
#' }
#'
#' @seealso
#' \code{\link{generate_aft_dgm_flex}} for creating the DGM
#'
#' @export
#' @importFrom stats rexp runif

simulate_from_dgm <- function(dgm,
                              n = NULL,
                              rand_ratio = 1,
                              entry_var = NULL,
                              max_entry = 24,
                              analysis_time = 48,
                              cens_adjust = 0,
                              draw_treatment = TRUE,
                              seed = NULL) {

  if (!inherits(dgm, c("aft_dgm_flex", "aft_dgm"))) {
    stop("dgm must be an object created by generate_aft_dgm_flex()")
  }

  if (!is.null(seed)) set.seed(seed)

  df_super <- dgm$df_super
  params <- dgm$model_params

  # Determine sample size
  if (is.null(n)) {
    df_sim <- df_super
    n <- nrow(df_sim)
  } else {
    # Sample from super population

    idx_sample <- sample(1:nrow(df_super), size = n, replace = TRUE)
    df_sim <- df_super[idx_sample, ]

    if(!is.null(entry_var)){
    entry_time <- df_sim[,c(entry_var)]
    } else {
    entry_time <- runif(n, min = 0, max = max_entry)
    }

    follow_up <- analysis_time - entry_time
    # Data at analysis_time: follow_up > 0

    # Reassign treatment if draw_treatment, otherwise retain per super_population
    # Set to original assignment in df_super
    df_sim$treat_sim <- df_sim$treat

    if(draw_treatment){
      n_treat <- round(n * rand_ratio / (1 + rand_ratio))
      n_control <- n - n_treat
      df_sim$treat_sim[1:n_treat] <- 1
      df_sim$treat_sim[(n_treat + 1):n] <- 0
    }

    # Update linear predictors based on assigned treatment
    df_sim$lin_pred_obs <- ifelse(df_sim$treat_sim == 1,
                                  df_sim$lin_pred_1,
                                  df_sim$lin_pred_0)

    # Reset IDs
    df_sim$id <- 1:n
  }

  # Generate survival times
  epsilon <- log(rexp(n))  # Extreme value distribution
  logT_sim <- params$mu + params$tau * epsilon + df_sim$lin_pred_obs
  T_sim <- exp(logT_sim)

  # Generate censoring times
  if (params$censoring$type == "weibull") {
    # Weibull censoring
    lin_pred_cens <- ifelse(df_sim$treat_sim == 1,
                            df_sim$lin_pred_cens_1,
                            df_sim$lin_pred_cens_0)
    epsilon_cens <- log(rexp(n))  # Extreme value distribution
    logC_sim <- params$censoring$mu + cens_adjust +
      params$censoring$tau * epsilon_cens + lin_pred_cens
    C_sim <- exp(logC_sim)

  } else if (params$censoring$type == "lognormal") {
    # Lognormal censoring
    lin_pred_cens <- ifelse(df_sim$treat_sim == 1,
                            df_sim$lin_pred_cens_1,
                            df_sim$lin_pred_cens_0)
    epsilon_cens <- rnorm(n)  # Normal errors for lognormal
    logC_sim <- params$censoring$mu + cens_adjust +
      params$censoring$tau * epsilon_cens + lin_pred_cens
    C_sim <- exp(logC_sim)

  } else if (params$censoring$type == "uniform") {
    # Uniform censoring
    C_sim <- runif(n,
                   min = params$censoring$min,
                   max = params$censoring$max)
  }

  # Apply administrative censoring
  C_sim <- pmin(C_sim, follow_up)

  # Observed times and events
  df_sim$y_sim <- pmin(T_sim, C_sim)
  df_sim$event_sim <- ifelse(T_sim <= C_sim, 1, 0)
  df_sim$t_true <- T_sim
  df_sim$c_time <- C_sim

  # Data at analysis time "x = analysis_time"
  df_sim <- subset(df_sim, follow_up > 0)

  return(df_sim)
}
