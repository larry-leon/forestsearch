# =============================================================================
# simulate_from_glm_dgm() -- Generate trial data from a GLM DGM
# =============================================================================
#
# Mirrors simulate_from_dgm() for non-survival outcomes.
# Simpler because there is no censoring, entry times, or AFT structure.
# =============================================================================


#' Simulate Trial Data from a GLM Data Generating Mechanism
#'
#' Generates a simulated clinical trial dataset by sampling from the
#' super-population in a \code{"glm_dgm"} object.  Assigns treatment
#' and generates outcomes from the individual-level potential outcomes
#' stored in the DGM.
#'
#' @param dgm An object of class \code{"glm_dgm"} created by
#'   \code{\link{generate_glm_dgm}}.
#' @param n Integer. Sample size.  If \code{NULL} (default), uses the
#'   entire super-population without resampling.
#' @param rand_ratio Numeric. Treatment:control randomisation ratio.
#'   Default \code{1} (1:1 allocation).
#' @param seed Integer. Random seed.  Default \code{NULL} (no explicit
#'   seed; use the calling environment's RNG state).
#' @param draw_treatment Logical.  If \code{TRUE} (default), randomly
#'   assigns treatment.  If \code{FALSE}, retains the treatment column
#'   from the super-population.
#'
#' @return A data frame with columns:
#'   \describe{
#'     \item{\code{id}}{Integer subject identifier.}
#'     \item{\code{y_sim}}{Simulated outcome.}
#'     \item{\code{treat_sim}}{Treatment indicator (0/1).}
#'     \item{\code{flag_harm}}{True subgroup membership indicator.}
#'   }
#'   Plus all covariate columns from the super-population.
#'
#'   For binary outcomes, \code{y_sim} is a 0/1 Bernoulli draw from
#'   the potential outcome probabilities (\code{p0} or \code{p1}).
#'   There is no \code{event_sim} column -- for ForestSearch, set
#'   \code{event.name = "y_sim"}.
#'
#' @seealso \code{\link{generate_glm_dgm}},
#'   \code{\link{run_simulation_analysis}}
#'
#' @examples
#' \dontrun{
#' dgm <- generate_glm_dgm(...)
#' df  <- simulate_from_glm_dgm(dgm, n = 500, seed = 42)
#' table(df$y_sim, df$treat_sim)
#' }
#'
#' @export
#' @importFrom stats rbinom rnorm rpois runif
simulate_from_glm_dgm <- function(
    dgm,
    n               = NULL,
    rand_ratio      = 1,
    seed            = NULL,
    draw_treatment  = TRUE
) {

  # -- Validate --------------------------------------------------------------
  if (!inherits(dgm, "glm_dgm")) {
    stop("'dgm' must be an object of class 'glm_dgm'.")
  }

  df_super <- dgm$df_super

  if (!is.null(seed)) set.seed(seed)

  # -- Sample ----------------------------------------------------------------
  if (is.null(n)) {
    df_sim <- df_super
  } else {
    idx <- sample(nrow(df_super), n, replace = TRUE)
    df_sim <- df_super[idx, , drop = FALSE]
  }

  n_obs <- nrow(df_sim)
  df_sim$id <- seq_len(n_obs)
  rownames(df_sim) <- NULL

  # -- Assign treatment ------------------------------------------------------
  if (draw_treatment) {
    p_treat <- rand_ratio / (1 + rand_ratio)
    df_sim$treat_sim <- stats::rbinom(n_obs, 1L, p_treat)
  } else {
    # Retain treatment from super-population
    treat_col <- names(dgm$model_params$coefficients)[2]
    if (treat_col %in% names(df_sim)) {
      df_sim$treat_sim <- as.integer(df_sim[[treat_col]])
    } else {
      stop("Cannot find treatment column in super-population. ",
           "Set draw_treatment = TRUE.")
    }
  }

  # -- Generate outcome ------------------------------------------------------
  outcome_type <- dgm$outcome_type

  if (outcome_type == "binary") {
    p_obs <- ifelse(df_sim$treat_sim == 1L, df_sim$p1, df_sim$p0)
    df_sim$y_sim <- stats::rbinom(n_obs, 1L, p_obs)
  } else if (outcome_type == "continuous") {
    mu_obs <- ifelse(df_sim$treat_sim == 1L, df_sim$mu1, df_sim$mu0)
    sigma <- dgm$model_params$sigma
    if (is.null(sigma)) sigma <- 1
    df_sim$y_sim <- mu_obs + stats::rnorm(n_obs, sd = sigma)
  } else if (outcome_type == "count") {
    mu_obs <- ifelse(df_sim$treat_sim == 1L, df_sim$mu1, df_sim$mu0)
    df_sim$y_sim <- stats::rpois(n_obs, lambda = mu_obs)
  } else {
    stop("Unknown outcome_type: '", outcome_type, "'")
  }

  df_sim
}
