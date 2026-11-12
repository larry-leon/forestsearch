# =============================================================================
# setup_gbsg_dgm() — GBSG DGM wrapper returning an aft_dgm_flex object
# =============================================================================
#
# This function is the new public API for creating a GBSG-based DGM.
# Internally it calls the builder .create_gbsg_dgm_() (unchanged), guaranteeing
# identical numeric output.  The returned object is reshaped to class
# "aft_dgm_flex" so it can be passed directly to simulate_from_dgm().
#
# Backward compatibility:
#   - df_super_rand is retained, so compute_dgm_cde() and print.gbsg_dgm work.
#   - All legacy top-level fields (hr_H_true, AHR, etc.) are preserved.
#   - class is c("aft_dgm_flex", "gbsg_dgm", "list"), so both dispatch methods
#     and inherits() checks work.
#
# Migration path:
#   Step 1 (this file): setup_gbsg_dgm() wraps .create_gbsg_dgm_().
#   Later steps:        switch internals to generate_aft_dgm_flex() once
#                       simulate_from_dgm() compatibility is confirmed.
# =============================================================================


#' Set Up a GBSG-Based AFT Data Generating Mechanism
#'
#' Creates a GBSG-based data generating mechanism that is fully compatible with
#' \code{\link{simulate_from_dgm}}.  It builds the mechanism with the
#' package's internal GBSG builder and returns it as an object of class
#' \code{c("aft_dgm_flex", "gbsg_dgm", "list")}.
#'
#' Internally the function calls that builder and then:
#' \enumerate{
#'   \item Adds a \code{df_super} field with column names aligned to
#'     \code{simulate_from_dgm()} conventions
#'     (\code{lin_pred_1}, \code{lin_pred_0}, \code{lin_pred_cens_1},
#'     \code{lin_pred_cens_0}, \code{flag_harm}).
#'   \item Adds a \code{model_params$tau} field (= \code{model_params$sigma})
#'     and a \code{model_params$censoring} sub-list.
#'   \item Sets class to \code{c("aft_dgm_flex", "gbsg_dgm", "list")}.
#' }
#' The original \code{df_super_rand} field is kept so that
#' \code{compute_dgm_cde()} and \code{print.gbsg_dgm} continue to work.
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
#' @param k_random_noise Integer. Number of standard-normal noise columns
#'   (\code{noise1}, \code{noise2}, ...) drawn once onto the super-population
#'   at construction. Default \code{0L} (no noise columns; output identical
#'   to previous versions).
#' @param noise_seed Integer. Seed for the population-level noise draw, used
#'   only when \code{k_random_noise > 0}. Default \code{20260807L}, the
#'   committed population-noise seed (commit \code{2a4787bc}).
#'
#' @details Noise columns requested via \code{k_random_noise} are inert
#'   N(0,1) population attributes intended as candidate confounders; the
#'   builder never puts them in the outcome model. They are drawn once onto
#'   the super-population (both \code{df_super_rand} and \code{df_super}),
#'   so simulated trials and evaluation frames inherit them by row resampling.
#'
#' @section Superseded:
#' This function is kept so that published simulations can be reproduced.
#' New work should use \code{\link{generate_aft_dgm_flex}()} with
#' \code{\link{simulate_from_dgm}()}, which take any survival dataset with a
#' 0/1 treatment column. On the GBSG data they give the same outcome model,
#' true hazard ratios and \code{k_inter} calibration once the caller has built
#' the GBSG covariates, at \code{k_z3 = 1} and \code{use_rand_params = FALSE}
#' (\code{generate_aft_dgm_flex()} has neither argument); there is no drop-in
#' replacement call. The censoring step differs: here (with the default
#' \code{cens_type = "weibull"}) a Weibull censoring model is fitted to the
#' resampled super-population, in which treatment is the arm assigned for
#' simulation, not the observed treatment, and its covariate term is the
#' fitted term multiplied by the ratio of the event-model scale to the
#' censoring-model scale. The true hazard ratios stored in the object and the
#' \code{k_inter} calibration are not affected; simulated censoring times are,
#' and with them the observed times and event indicators and anything computed
#' from them, such as the beta(Hhat) targets scored on a frame from
#' \code{\link{fs_build_eval_frame}()}. The "Extreme Subgroup
#' Effects" article on the package website shows the general machinery on the
#' GBSG data, with a uniform treatment effect. This function will be deprecated
#' in a future release.
#'
#' @return An object of class \code{c("aft_dgm_flex", "gbsg_dgm", "list")}
#'   with the builder's fields (\code{df_super_rand}; the true hazard ratios
#'   \code{hr_causal}, \code{hr_H_true}, \code{hr_Hc_true}, \code{AHR},
#'   \code{AHR_H_true}, \code{AHR_Hc_true} and the list
#'   \code{hazard_ratios}; \code{model_params}, \code{cens_params},
#'   \code{subgroup_info}, \code{effect_modifiers}, \code{analysis_vars},
#'   \code{model_type}, \code{n_super} and \code{seed}) plus:
#'   \describe{
#'     \item{\code{df_super}}{Super-population data frame with
#'       \code{simulate_from_dgm()}-compatible column names.}
#'     \item{\code{model_params$tau}}{Copy of \code{model_params$sigma}.}
#'     \item{\code{model_params$censoring}}{Sub-list with
#'       \code{type}, \code{mu}, \code{tau} for the censoring model.}
#'     \item{\code{noise_names}}{Character vector of noise column names
#'       (\code{character(0)} when \code{k_random_noise = 0}).}
#'     \item{\code{noise_scheme}}{\code{"population"} when noise was drawn,
#'       otherwise \code{"none"}.}
#'     \item{\code{noise_seed}}{The seed used for the noise draw, or
#'       \code{NA_integer_} when none was drawn.}
#'   }
#'
#' @seealso \code{\link{generate_aft_dgm_flex}}, \code{\link{simulate_from_dgm}},
#'   \code{\link{compute_dgm_cde}}
#'
#' @examples
#' \dontrun{
#' dgm <- setup_gbsg_dgm(model = "alt", k_inter = 2, verbose = FALSE)
#' dgm <- compute_dgm_cde(dgm)
#' print(dgm)
#' sim <- simulate_from_dgm(dgm, n = 400, seed = 1)
#' }
#'
#' @export
#' @importFrom stats rnorm
setup_gbsg_dgm <- function(
    model           = c("alt", "null"),
    k_treat         = 1,
    k_inter         = 1,
    k_z3            = 1,
    z1_quantile     = 0.25,
    n_super         = 5000L,
    cens_type       = c("weibull", "uniform"),
    use_rand_params = FALSE,
    seed            = 8316951L,
    verbose         = FALSE,
    k_random_noise  = 0L,
    noise_seed      = 20260807L
) {
  model     <- match.arg(model)
  cens_type <- match.arg(cens_type)

  # ── Step 1: Run the legacy DGM (unchanged logic, guaranteed exact numbers) ──
  dgm <- .create_gbsg_dgm_(
    model           = model,
    k_treat         = k_treat,
    k_inter         = k_inter,
    k_z3            = k_z3,
    z1_quantile     = z1_quantile,
    n_super         = n_super,
    cens_type       = cens_type,
    use_rand_params = use_rand_params,
    seed            = seed,
    verbose         = verbose
  )

  # ── Step 1b: Population-level noise columns (scheme of commit 2a4787bc) ────
  #
  # Drawn once onto df_super_rand BEFORE the Step-2 copy, so df_super inherits
  # the identical values and both frames stay consistent.  The global RNG
  # state is saved/restored so a k_random_noise > 0 build leaves the caller's
  # RNG stream exactly as a k_random_noise = 0 build does.
  noise_names <- character(0)
  if (k_random_noise > 0) {
    noise_names <- paste0("noise", seq_len(k_random_noise))
    has_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    if (has_seed) old_seed <- get(".Random.seed", envir = globalenv())
    set.seed(noise_seed)
    for (nm in noise_names)
      dgm$df_super_rand[[nm]] <- stats::rnorm(nrow(dgm$df_super_rand))
    if (has_seed) assign(".Random.seed", old_seed, envir = globalenv())
    else rm(".Random.seed", envir = globalenv())
  }

  # ── Step 2: Build df_super with simulate_from_dgm()-compatible column names ──
  #
  # simulate_from_dgm() reads:
  #   df_super$lin_pred_1  / lin_pred_0        (AFT linear predictors)
  #   df_super$lin_pred_cens_1 / lin_pred_cens_0  (censoring linear predictors)
  #   df_super$flag_harm                       (subgroup indicator)
  #
  # df_super_rand uses dot-notation equivalents:
  #   lin1.conf, lin0.conf, linC1.conf, linC0.conf, flag.harm
  #
  df_s <- dgm$df_super_rand

  .rename <- function(df, old, new) {
    if (old %in% names(df) && !new %in% names(df))
      names(df)[names(df) == old] <- new
    df
  }

  df_s <- .rename(df_s, "lin1.conf",  "lin_pred_1")
  df_s <- .rename(df_s, "lin0.conf",  "lin_pred_0")
  df_s <- .rename(df_s, "linC1.conf", "lin_pred_cens_1")
  df_s <- .rename(df_s, "linC0.conf", "lin_pred_cens_0")
  df_s <- .rename(df_s, "flag.harm",  "flag_harm")

  dgm$df_super <- df_s          # add aligned field
  # df_super_rand is kept for backward compat (compute_dgm_cde, print.gbsg_dgm)

  # ── Step 3: Extend model_params for simulate_from_dgm() ────────────────────
  #
  # simulate_from_dgm() reads:
  #   params$mu               (AFT intercept)          -- already present
  #   params$tau              (AFT scale)               -- legacy uses $sigma
  #   params$censoring$type   ("weibull" / "uniform")
  #   params$censoring$mu     (censoring intercept)
  #   params$censoring$tau    (censoring scale)
  #
  dgm$model_params$tau <- dgm$model_params$sigma   # alias

  dgm$model_params$censoring <- list(
    type = dgm$cens_params$type,
    mu   = dgm$cens_params$mu,
    tau  = dgm$cens_params$sigma   # cens_params uses $sigma; sim_from_dgm needs $tau
  )

  # ── Step 4: Promote class to aft_dgm_flex (keeps gbsg_dgm for dispatch) ────
  class(dgm) <- c("aft_dgm_flex", "gbsg_dgm", "list")

  # ── Step 5: Record population-noise provenance ─────────────────────────────
  dgm$noise_names  <- noise_names
  dgm$noise_scheme <- if (k_random_noise > 0) "population" else "none"
  dgm$noise_seed   <- if (k_random_noise > 0) noise_seed else NA_integer_

  dgm
}
