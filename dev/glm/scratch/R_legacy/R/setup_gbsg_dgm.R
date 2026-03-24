# =============================================================================
# setup_gbsg_dgm() — GBSG DGM wrapper returning an aft_dgm_flex object
# =============================================================================
#
# This function is the new public API for creating a GBSG-based DGM.
# Internally it calls create_gbsg_dgm() (unchanged), guaranteeing identical
# numeric output.  The returned object is reshaped to class "aft_dgm_flex" so
# it can be passed directly to simulate_from_dgm().
#
# Backward compatibility:
#   - df_super_rand is retained, so compute_dgm_cde() and print.gbsg_dgm work.
#   - All legacy top-level fields (hr_H_true, AHR, etc.) are preserved.
#   - class is c("aft_dgm_flex", "gbsg_dgm", "list"), so both dispatch methods
#     and inherits() checks work.
#
# Migration path:
#   Step 1 (this file): setup_gbsg_dgm() wraps create_gbsg_dgm().
#   Later steps:        switch internals to generate_aft_dgm_flex() once
#                       simulate_from_dgm() compatibility is confirmed.
# =============================================================================


#' Set Up a GBSG-Based AFT Data Generating Mechanism
#'
#' Creates a GBSG-based data generating mechanism that is fully compatible with
#' \code{\link{simulate_from_dgm}}.  This is the replacement for
#' \code{create_gbsg_dgm()}: it accepts exactly the same arguments and produces
#' the same numeric output, but returns an object of class
#' \code{"aft_dgm_flex"} instead of \code{"gbsg_dgm"}.
#'
#' Internally the function calls \code{create_gbsg_dgm()} and then:
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
#' @inheritParams create_gbsg_dgm
#'
#' @return An object of class \code{c("aft_dgm_flex", "gbsg_dgm", "list")}
#'   with all fields from \code{create_gbsg_dgm()} plus:
#'   \describe{
#'     \item{\code{df_super}}{Super-population data frame with
#'       \code{simulate_from_dgm()}-compatible column names.}
#'     \item{\code{model_params$tau}}{Copy of \code{model_params$sigma}.}
#'     \item{\code{model_params$censoring}}{Sub-list with
#'       \code{type}, \code{mu}, \code{tau} for the censoring model.}
#'   }
#'
#' @seealso \code{\link{create_gbsg_dgm}}, \code{\link{simulate_from_dgm}},
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
    verbose         = FALSE
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

  dgm
}
