# =============================================================================
# run_simulation_analysis() — General simulation wrapper
# =============================================================================
#
# Replaces oc_analyses_gbsg.R::run_simulation_analysis(), which was coupled
# to simulate_from_gbsg_dgm() and hardcoded column names (y.sim, event.sim,
# flag.harm, treat).
#
# Design:
#   - Calls simulate_from_dgm() (the general simulator).
#   - Column names are explicit parameters; no GBSG-specific names anywhere.
#   - Helpers (extract_fs_estimates_gen, extract_grf_estimates_gen,
#     run_forestsearch_analysis_gen, run_grf_analysis_gen) are contained here
#     and receive column names from the top-level call.
#   - default_sim_params() / default_grf_params_gen() use general names.
#   - GBSG usage is just: run_simulation_analysis(dgm, ...) with defaults.
# =============================================================================


# =============================================================================
# Default parameter lists
# =============================================================================

# Null-coalescing helper (local to this file).
# Defined here rather than importing rlang::%||% because this file
# runs inside parallel workers where rlang may not be loaded.
`%||%` <- function(a, b) if (!is.null(a)) a else b


#' Default ForestSearch parameters (general)
#' @keywords internal
default_sim_params <- function() {
  list(
    outcome.name             = "y_sim",
    event.name               = "event_sim",
    treat.name               = "treat_sim",
    id.name                  = "id",
    use_lasso                = TRUE,
    use_grf                  = FALSE,
    hr.threshold             = 1.25,
    hr.consistency           = 1.0,
    pconsistency.threshold   = 0.90,
    fs.splits                = 400,
    n.min                    = 60,
    d0.min                   = 12,
    d1.min                   = 12,
    maxk                     = 2,
    max.minutes              = 5,
    by.risk                  = 12,
    vi.grf.min               = -0.2,
    use_twostage             = FALSE,
    twostage_args            = list()
  )
}

#' Default GRF parameters (general)
#' @keywords internal
default_grf_params_gen <- function() {
  list(
    outcome.name = "y_sim",
    event.name   = "event_sim",
    treat.name   = "treat_sim",
    id.name      = "id",
    n.min        = 60,
    dmin.grf     = 12,
    frac.tau     = 0.60,
    maxdepth     = 2,
    RCT          = TRUE,
    sg.criterion = "mDiff",
    seedit       = 8316951L
  )
}


# =============================================================================
# Main function
# =============================================================================

#' Run One Simulation Replicate
#'
#' General replacement for the legacy \code{run_simulation_analysis()} that
#' was coupled to \code{simulate_from_gbsg_dgm()} and GBSG-specific column
#' names.  This version calls \code{\link{simulate_from_dgm}} and accepts
#' explicit column-name parameters, making it applicable to any DGM built
#' with \code{\link{generate_aft_dgm_flex}}.
#'
#' @param sim_id Integer. Simulation replicate index (used as seed offset).
#' @param dgm An \code{"aft_dgm_flex"} object from
#'   \code{\link{generate_aft_dgm_flex}} or \code{\link{setup_gbsg_dgm}}.
#' @param n_sample Integer. Per-replicate sample size.
#' @param analysis_time Numeric. Calendar time of analysis on the DGM time
#'   scale.  Use \code{Inf} (default) for no administrative censoring —
#'   equivalent to the legacy \code{max_follow = Inf}.
#' @param cens_adjust Numeric. Log-scale shift to censoring times passed to
#'   \code{simulate_from_dgm(cens_adjust = ...)}. Replaces legacy
#'   \code{muC_adj}. Default \code{0}.
#' @param max_follow \strong{Deprecated.} Use \code{analysis_time} instead.
#'   If supplied, its value is forwarded to \code{analysis_time} with a
#'   warning. Retained for backward compatibility with legacy scripts.
#' @param muC_adj \strong{Deprecated.} Use \code{cens_adjust} instead.
#'   If supplied, its value is forwarded to \code{cens_adjust} with a
#'   warning. Retained for backward compatibility with legacy scripts.
#' @param confounders_base Character vector of base confounder names.
#' @param n_add_noise Integer. Number of independent N(0,1) noise variables
#'   to append. Default \code{0L}.
#' @param outcome_name Name of the observed time column in simulated data.
#'   Default \code{"y_sim"}.
#' @param event_name Name of the event indicator column. Default
#'   \code{"event_sim"}.
#' @param treat_name Name of the treatment column. Default \code{"treat_sim"}.
#' @param harm_col Name of the true-subgroup indicator column. Default
#'   \code{"flag_harm"}.
#' @param run_fs Logical. Run ForestSearch (LASSO). Default \code{TRUE}.
#' @param run_fs_grf Logical. Run ForestSearch (LASSO + GRF). Default
#'   \code{TRUE}.
#' @param run_grf Logical. Run standalone GRF. Default \code{TRUE}.
#' @param fs_params Named list of ForestSearch parameter overrides.
#' @param grf_params Named list of GRF parameter overrides.
#' @param cox_formula Optional Cox formula for unadjusted ITT.
#' @param cox_formula_adj Optional adjusted Cox formula.
#' @param n_sims_total Integer. Total simulations (for progress messages).
#' @param seed_base Integer. Base seed; replicate seed = \code{seed_base +
#'   sim_id}. Default \code{8316951L}.
#' @param verbose Logical. Print progress messages. Default \code{FALSE}.
#' @param verbose_n Integer. If set, only print for \code{sim_id <=
#'   verbose_n}. Default \code{NULL}.
#' @param debug Logical. Print detailed debug output. Default \code{FALSE}.
#'
#' @return A \code{data.table} with one row per analysis method, containing
#'   subgroup size, HR, AHR, CDE, and classification metrics.
#'
#' @seealso \code{\link{simulate_from_dgm}},
#'   \code{\link{generate_aft_dgm_flex}}, \code{\link{setup_gbsg_dgm}}
#'
#' @section GLM Parameters:
#' GLM-specific parameters (\code{outcome_type}, \code{effect_measure},
#' \code{offset.name}) must be passed inside \code{fs_params}, not as
#' top-level arguments.  They route only to the estimation step
#' (\code{.extract_fs_estimates_gen}, \code{.extract_grf_estimates_gen}),
#' not to \code{forestsearch()} itself, which uses Cox PH for subgroup
#' identification in v0.1.x.  Passing these as top-level arguments will
#' result in them being silently ignored.
#'
#' @importFrom data.table data.table rbindlist
#' @importFrom survival coxph Surv
#' @importFrom stats rnorm
#' @examples
#' \dontrun{
#' dgm <- setup_gbsg_dgm(model = "null", verbose = FALSE)
#' results <- run_simulation_analysis(dgm, nsim = 5,
#'   parallel_args = list(plan = "sequential"))
#' summarize_simulation_results(results)
#' }
#' @export
run_simulation_analysis <- function(
    sim_id,
    dgm,
    n_sample,
    analysis_time    = Inf,
    cens_adjust      = 0,
    # Deprecated legacy parameter names — emit a warning if used
    max_follow       = NULL,
    muC_adj          = NULL,
    confounders_base = c("v1", "v2", "v3", "v4", "v5", "v6", "v7"),
    n_add_noise      = 0L,
    outcome_name     = "y_sim",
    event_name       = "event_sim",
    treat_name       = "treat_sim",
    harm_col         = "flag_harm",
    run_fs           = TRUE,
    run_fs_grf       = TRUE,
    run_grf          = TRUE,
    fs_params        = list(),
    grf_params       = list(),
    cox_formula      = NULL,
    cox_formula_adj  = NULL,
    n_sims_total     = NULL,
    seed_base        = 8316951L,
    verbose          = FALSE,
    verbose_n        = NULL,
    debug            = FALSE
) {

  show_verbose <- verbose && (is.null(verbose_n) || sim_id <= verbose_n)

  # ── Deprecated parameter shims ─────────────────────────────────────────────
  if (!is.null(max_follow)) {
    warning(
      "'max_follow' is deprecated. Use 'analysis_time' instead.\n",
      "  Mapping: analysis_time = max_follow (", max_follow, ")",
      call. = FALSE
    )
    analysis_time <- max_follow
  }
  if (!is.null(muC_adj)) {
    warning(
      "'muC_adj' is deprecated. Use 'cens_adjust' instead.\n",
      "  Mapping: cens_adjust = muC_adj (", muC_adj, ")",
      call. = FALSE
    )
    cens_adjust <- muC_adj
  }

  if (show_verbose) {
    message("\n", paste(rep("=", 60), collapse = ""))
    message(sprintf("Simulation %d", sim_id))
    if (!is.null(n_sims_total))
      message(sprintf("  Progress: %d / %d (%.1f%%)",
                      sim_id, n_sims_total, 100 * sim_id / n_sims_total))
    message(paste(rep("=", 60), collapse = ""))
  }

  # ── Simulate data ──────────────────────────────────────────────────────────
  if (show_verbose)
    message(sprintf("\n[1] Simulating data (n=%d, analysis_time=%s, cens_adjust=%g)...",
                    n_sample, analysis_time, cens_adjust))

  sim_data <- tryCatch(
    simulate_from_dgm(
      dgm           = dgm,
      n             = n_sample,
      seed          = seed_base + sim_id,
      analysis_time = analysis_time,
      cens_adjust   = cens_adjust
    ),
    error = function(e) stop("simulate_from_dgm failed: ", e$message)
  )

  if (show_verbose)
    message(sprintf("    n=%d  events=%d (%.1f%%)",
                    nrow(sim_data), sum(sim_data[[event_name]]),
                    100 * mean(sim_data[[event_name]])))

  # ── Optional noise variables ───────────────────────────────────────────────
  confounders_name <- confounders_base
  if (n_add_noise > 0L) {
    set.seed(seed_base + 1000L * sim_id)
    noise_nms <- paste0("noise", seq_len(n_add_noise))
    for (nm in noise_nms) sim_data[[nm]] <- stats::rnorm(nrow(sim_data))
    confounders_name <- c(confounders_base, noise_nms)
    if (show_verbose)
      message(sprintf("    Added %d noise variable(s)", n_add_noise))
  }

  # ── True subgroup properties ───────────────────────────────────────────────
  has_harm <- harm_col %in% names(sim_data)
  size_H_true  <- if (has_harm) sum(sim_data[[harm_col]])        else NA_integer_
  prop_H_true  <- if (has_harm) mean(sim_data[[harm_col]])       else NA_real_
  size_Hc_true <- if (has_harm) nrow(sim_data) - size_H_true    else NA_integer_
  prop_Hc_true <- if (has_harm) 1 - prop_H_true                 else NA_real_

  hr_H_dgm   <- get_dgm_hr(dgm, "hr_H")
  hr_Hc_dgm  <- get_dgm_hr(dgm, "hr_Hc")
  ahr_H_dgm  <- get_dgm_hr(dgm, "ahr_H")
  ahr_Hc_dgm <- get_dgm_hr(dgm, "ahr_Hc")

  if (show_verbose) {
    message(sprintf("\n[2] True subgroup: H n=%s (%.1f%%), Hc n=%s (%.1f%%)",
                    ifelse(is.na(size_H_true), "NA", size_H_true),
                    100 * prop_H_true,
                    ifelse(is.na(size_Hc_true), "NA", size_Hc_true),
                    100 * prop_Hc_true))
    message(sprintf("    DGM HR_H=%.3f  HR_Hc=%.3f  AHR_H=%.3f  AHR_Hc=%.3f",
                    hr_H_dgm, hr_Hc_dgm,
                    ifelse(is.na(ahr_H_dgm), NA, ahr_H_dgm),
                    ifelse(is.na(ahr_Hc_dgm), NA, ahr_Hc_dgm)))
  }

  df_pop <- data.table::data.table(
    sim          = sim_id,
    sizeH_true   = size_H_true,
    propH_true   = prop_H_true,
    sizeHc_true  = size_Hc_true,
    propHc_true  = prop_Hc_true
  )

  # ── Merge parameter defaults ───────────────────────────────────────────────
  fs_defaults  <- default_sim_params()
  grf_defaults <- default_grf_params_gen()
  grf_merged   <- modifyList(grf_defaults, grf_params)

  # Ensure outcome/event/treat names in FS defaults match what simulate produced
  fs_defaults$outcome.name <- outcome_name
  fs_defaults$event.name   <- event_name
  fs_defaults$treat.name   <- treat_name
  grf_merged$outcome.name  <- outcome_name
  grf_merged$event.name    <- event_name
  grf_merged$treat.name    <- treat_name

  # ── Run analyses ───────────────────────────────────────────────────────────
  results_list <- list()

  if (show_verbose) {
    to_run <- c(if (run_fs) "FS", if (run_fs_grf) "FSlg", if (run_grf) "GRF")
    message(sprintf("\n[3] Running: %s", paste(to_run, collapse = ", ")))
  }

  if (run_fs) {
    fs_p <- modifyList(
      modifyList(fs_defaults, list(use_lasso = TRUE, use_grf = FALSE)),
      fs_params
    )
    results_list[["FS"]] <- cbind(df_pop,
      .run_fs_analysis_gen(
        data = sim_data, confounders_name = confounders_name,
        params = fs_p, dgm = dgm,
        cox_formula = cox_formula, cox_formula_adj = cox_formula_adj,
        outcome_name = outcome_name, event_name = event_name,
        treat_name = treat_name, harm_col = harm_col,
        label = "FS", verbose = show_verbose
      ))
  }

  if (run_fs_grf) {
    fs_p <- modifyList(
      modifyList(fs_defaults, list(use_lasso = TRUE, use_grf = TRUE)),
      fs_params
    )
    results_list[["FSlg"]] <- cbind(df_pop,
      .run_fs_analysis_gen(
        data = sim_data, confounders_name = confounders_name,
        params = fs_p, dgm = dgm,
        cox_formula = cox_formula, cox_formula_adj = cox_formula_adj,
        outcome_name = outcome_name, event_name = event_name,
        treat_name = treat_name, harm_col = harm_col,
        label = "FSlg", verbose = show_verbose
      ))
  }

  if (run_grf) {
    # Resolve id_name from merged GRF params (falls back to "id")
    id_name_resolved <- grf_merged$id.name %||% "id"
    results_list[["GRF"]] <- cbind(df_pop,
      .run_grf_analysis_gen(
        data = sim_data, confounders_name = confounders_name,
        params = grf_merged, dgm = dgm,
        cox_formula = cox_formula, cox_formula_adj = cox_formula_adj,
        outcome_name = outcome_name, event_name = event_name,
        treat_name = treat_name, harm_col = harm_col,
        label = "GRF", verbose = show_verbose, debug = debug,
        outcome_type   = fs_params$outcome_type,
        effect_measure = fs_params$effect_measure,
        offset_name    = fs_params$offset.name,
        id_name        = id_name_resolved
      ))
  }

  if (length(results_list) == 0) {
    warning("No analyses were run. Check run_fs / run_fs_grf / run_grf.")
    return(NULL)
  }

  result <- data.table::rbindlist(results_list, fill = TRUE)

  if (show_verbose)
    message(sprintf("\n[4] Done: %d rows x %d cols\n%s",
                    nrow(result), ncol(result),
                    paste(rep("=", 60), collapse = "")))

  result
}


# =============================================================================
# Internal: ForestSearch analysis
# =============================================================================

#' @keywords internal
.run_fs_analysis_gen <- function(
    data, confounders_name, params, dgm,
    cox_formula, cox_formula_adj,
    outcome_name, event_name, treat_name, harm_col,
    label, verbose
) {
  if (verbose) {
    message(sprintf("\n  [%s] ForestSearch: n=%d events=%d (%.1f%%)",
                    label, nrow(data),
                    sum(data[[event_name]]), 100 * mean(data[[event_name]])))
  }

  fs_args <- list(
    df.analysis      = data,
    confounders.name = confounders_name,
    details          = verbose,
    plot.sg          = FALSE
  )

  valid_pnames <- c(
    "outcome.name", "event.name", "treat.name", "id.name",
    "use_lasso", "use_grf", "conf_force",
    "hr.threshold", "hr.consistency", "pconsistency.threshold",
    "fs.splits", "n.min", "d0.min", "d1.min",
    "maxk", "max.minutes", "by.risk", "vi.grf.min",
    "frac.tau", "dmin.grf", "grf_depth",
    "use_twostage", "twostage_args",
    "adverse_outcome", "ps_method", "ps_adjust_method", "ps_hat",
    "parallel_args"
    # NOTE: outcome_type, effect_measure, offset.name are intentionally
    # EXCLUDED from the forestsearch() whitelist.  They affect only the
    # *estimation* step (.extract_fs_estimates_gen), not subgroup
    # *identification* — which always uses Cox PH in forestsearch v0.1.x.
    # Add them back here once forestsearch() natively supports GLM.
  )
  for (pn in valid_pnames)
    if (!is.null(params[[pn]])) fs_args[[pn]] <- params[[pn]]

  fs_args$quiet <- TRUE
  fs_result <- tryCatch(
    do.call(forestsearch, fs_args),
    error = function(e) { warning(label, " failed: ", e$message); NULL }
  )

  has_result <- !is.null(fs_result) &&
    !is.null(fs_result$grp.consistency$out_sg$result) &&
    nrow(fs_result$grp.consistency$out_sg$result) > 0

  .extract_fs_estimates_gen(
    df           = data,
    fs_res       = if (has_result) fs_result$grp.consistency$out_sg$result else NULL,
    fs_full      = if (has_result) fs_result else NULL,
    dgm          = dgm,
    cox_formula  = cox_formula,
    cox_formula_adj = cox_formula_adj,
    outcome_name = outcome_name,
    event_name   = event_name,
    treat_name   = treat_name,
    harm_col     = harm_col,
    analysis     = label,
    verbose      = verbose,
    outcome_type   = params$outcome_type,
    effect_measure = params$effect_measure,
    offset_name    = params$offset.name
  )
}


# =============================================================================
# Internal: extract FS estimates  (column-name-general)
# =============================================================================

#' @keywords internal
.extract_fs_estimates_gen <- function(
    df, fs_res, fs_full, dgm,
    cox_formula, cox_formula_adj,
    outcome_name, event_name, treat_name, harm_col,
    analysis, verbose,
    outcome_type = NULL, effect_measure = NULL, offset_name = NULL
) {

  # ── Detect GLM mode ──────────────────────────────────────────────────────
  is_glm <- !is.null(outcome_type) && outcome_type != "survival"

  out <- data.table::data.table(
    analysis    = analysis,
    any.H       = 0L,
    size.H      = NA_integer_,
    size.Hc     = nrow(df),
    hr.H.true   = NA_real_, hr.H.hat   = NA_real_,
    hr.Hc.true  = NA_real_, hr.Hc.hat  = NA_real_,
    ahr.H.true  = NA_real_, ahr.Hc.true = NA_real_,
    ahr.H.hat   = NA_real_, ahr.Hc.hat  = NA_real_,
    cde.H.true  = NA_real_, cde.Hc.true = NA_real_,
    cde.H.hat   = NA_real_, cde.Hc.hat  = NA_real_,
    hr.itt      = NA_real_,
    hr.adj.itt  = NA_real_,
    ppv = NA_real_, npv = NA_real_, sens = NA_real_, spec = NA_real_,
    p.cens = 1 - mean(df[[event_name]]),
    taumax = max(df[[outcome_name]])
  )

  # ── Build estimation closures ─────────────────────────────────────────────
  if (is_glm) {
    # Resolve GLM family from effect_measure
    glm_family <- switch(
      effect_measure %||% "OR",
      "OR"  = , "RD"  = stats::binomial(),
      "IRR" = , "IRD" = stats::poisson(),
      "MD"  =          stats::gaussian(),
      stats::binomial()
    )
    # Whether the effect is on log scale (exponentiate) or identity
    log_scale <- effect_measure %in% c("OR", "IRR", "RR")

    .glm_effect <- function(sub_df) {
      tryCatch({
        if (!is.null(offset_name) && offset_name %in% names(sub_df)) {
          off_vec <- log(pmax(sub_df[[offset_name]], 1e-6))
          fit <- stats::glm(sub_df[[event_name]] ~ sub_df[[treat_name]],
                            family = glm_family, offset = off_vec)
        } else {
          fit <- stats::glm(sub_df[[event_name]] ~ sub_df[[treat_name]],
                            family = glm_family)
        }
        b <- coef(fit)[[2]]  # treatment coefficient
        if (log_scale) exp(b) else b
      }, error = function(e) NA_real_)
    }

    # ITT
    out$hr.itt <- .glm_effect(df)
  } else {
    # ── Survival (Cox) ITT ─────────────────────────────────────────────────
    out$hr.itt <- tryCatch(
      exp(survival::coxph(
        survival::Surv(df[[outcome_name]], df[[event_name]]) ~ df[[treat_name]]
      )$coefficients),
      error = function(e) NA_real_
    )
  }

  # Adjusted ITT (Cox only; skip for GLM)
  if (!is_glm && !is.null(cox_formula_adj)) {
    out$hr.adj.itt <- tryCatch(
      exp(survival::coxph(cox_formula_adj, data = df)$coefficients[treat_name]),
      error = function(e) NA_real_
    )
  }

  if (is.null(fs_res) || nrow(fs_res) == 0) {
    out$hr.Hc.hat <- out$hr.itt
    return(out)
  }

  out$any.H <- 1L

  # Build subgroup indicator
  df_sg <- NULL
  if (!is.null(fs_full) && !is.null(fs_full$df.est) &&
      "treat.recommend" %in% names(fs_full$df.est)) {
    df$sg_hat <- as.integer(fs_full$df.est$treat.recommend == 0)
  } else {
    sg_factors <- character(0)
    for (i in seq_len(min(7, ncol(fs_res)))) {
      cn <- paste0("M.", i)
      if (cn %in% names(fs_res) && !is.na(fs_res[[cn]][1]))
        sg_factors <- c(sg_factors, fs_res[[cn]][1])
    }
    if (length(sg_factors) > 0)
      df$sg_hat <- create_subgroup_indicator(df, sg_factors)
    else {
      out$hr.Hc.hat <- out$hr.itt
      return(out)
    }
  }

  out$size.H  <- sum(df$sg_hat, na.rm = TRUE)
  out$size.Hc <- sum(df$sg_hat == 0L, na.rm = TRUE)

  # ── Subgroup effect estimates ────────────────────────────────────────────
  if (is_glm) {
    if (out$size.H  > 10) out$hr.H.hat  <- .glm_effect(df[df$sg_hat == 1, ])
    if (out$size.Hc > 10) out$hr.Hc.hat <- .glm_effect(df[df$sg_hat == 0, ])
  } else {
    .cox_hr <- function(sub_df) tryCatch(
      exp(survival::coxph(
        survival::Surv(sub_df[[outcome_name]], sub_df[[event_name]]) ~
          sub_df[[treat_name]]
      )$coefficients),
      error = function(e) NA_real_
    )

    if (out$size.H  > 10) out$hr.H.hat  <- .cox_hr(df[df$sg_hat == 1, ])
    if (out$size.Hc > 10) out$hr.Hc.hat <- .cox_hr(df[df$sg_hat == 0, ])
  }

  # ── AHR and CDE (survival-only estimands) ────────────────────────────────
  if (!is_glm) {
    if ("loghr_po" %in% names(df)) {
      out$ahr.H.hat  <- compute_ahr(df, df$sg_hat)
      out$ahr.Hc.hat <- compute_ahr(df, 1L - df$sg_hat)
    }
    if (all(c("theta_0", "theta_1") %in% names(df))) {
      out$cde.H.hat  <- compute_cde(df, df$sg_hat)
      out$cde.Hc.hat <- compute_cde(df, 1L - df$sg_hat)
    }
  }

  # ── True-subgroup classification metrics ─────────────────────────────────
  if (harm_col %in% names(df)) {
    true_H <- df[[harm_col]] == 1L
    hat_H  <- df$sg_hat == 1L

    tp <- sum(true_H & hat_H);   fp <- sum(!true_H & hat_H)
    tn <- sum(!true_H & !hat_H); fn <- sum(true_H & !hat_H)

    out$sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    out$spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
    out$ppv  <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    out$npv  <- if ((tn + fn) > 0) tn / (tn + fn) else NA_real_

    # True-subgroup effect estimates
    if (is_glm) {
      if (sum(true_H) > 10)  out$hr.H.true  <- .glm_effect(df[true_H, ])
      if (sum(!true_H) > 10) out$hr.Hc.true <- .glm_effect(df[!true_H, ])
    } else {
      if (sum(true_H) > 10)  out$hr.H.true  <- .cox_hr(df[true_H, ])
      if (sum(!true_H) > 10) out$hr.Hc.true <- .cox_hr(df[!true_H, ])

      if ("loghr_po" %in% names(df)) {
        out$ahr.H.true  <- compute_ahr(df, df[[harm_col]])
        out$ahr.Hc.true <- compute_ahr(df, 1L - df[[harm_col]])
      }
      if (all(c("theta_0", "theta_1") %in% names(df))) {
        out$cde.H.true  <- compute_cde(df, df[[harm_col]])
        out$cde.Hc.true <- compute_cde(df, 1L - df[[harm_col]])
      }
    }
  }

  out
}


# =============================================================================
# Internal: GRF analysis
# =============================================================================

#' @keywords internal
.run_grf_analysis_gen <- function(
    data, confounders_name, params, dgm,
    cox_formula, cox_formula_adj,
    outcome_name, event_name, treat_name, harm_col,
    label, verbose, debug,
    outcome_type = NULL, effect_measure = NULL, offset_name = NULL,
    id_name = "id"
) {
  grf_fun <- tryCatch(
    get("grf.subg.harm.survival", mode = "function",
        envir = asNamespace("forestsearch")),
    error = function(e) NULL
  )

  if (is.null(grf_fun)) {
    warning("grf.subg.harm.survival not found. Skipping GRF analysis.")
    return(.extract_grf_estimates_gen(
      df = data, grf_est = NULL, dgm = dgm,
      cox_formula = cox_formula, cox_formula_adj = cox_formula_adj,
      outcome_name = outcome_name, event_name = event_name,
      treat_name = treat_name, harm_col = harm_col,
      analysis = label, verbose = verbose, debug = debug,
      outcome_type = outcome_type, effect_measure = effect_measure,
      offset_name = offset_name, id_name = id_name
    ))
  }

  grf_args <- list(data = data, confounders.name = confounders_name,
                   details = verbose)
  grf_pnames <- c("outcome.name", "event.name", "id.name", "treat.name",
                  "frac.tau", "n.min", "dmin.grf", "RCT", "sg.criterion",
                  "maxdepth", "seedit")
  for (pn in grf_pnames)
    if (!is.null(params[[pn]])) grf_args[[pn]] <- params[[pn]]

  grf_result <- tryCatch(
    do.call(grf_fun, grf_args),
    error = function(e) { warning(label, " failed: ", e$message); NULL }
  )

  .extract_grf_estimates_gen(
    df = data, grf_est = grf_result, dgm = dgm,
    cox_formula = cox_formula, cox_formula_adj = cox_formula_adj,
    outcome_name = outcome_name, event_name = event_name,
    treat_name = treat_name, harm_col = harm_col,
    analysis = label, verbose = verbose, debug = debug,
    outcome_type = outcome_type, effect_measure = effect_measure,
    offset_name = offset_name, id_name = id_name
  )
}


# =============================================================================
# Internal: extract GRF estimates  (column-name-general)
# =============================================================================

#' @keywords internal
.extract_grf_estimates_gen <- function(
    df, grf_est, dgm,
    cox_formula, cox_formula_adj,
    outcome_name, event_name, treat_name, harm_col,
    analysis, verbose, debug,
    outcome_type = NULL, effect_measure = NULL, offset_name = NULL,
    id_name = "id"
) {

  # ── Detect GLM mode ──────────────────────────────────────────────────────
  is_glm <- !is.null(outcome_type) && outcome_type != "survival"

  out <- data.table::data.table(
    analysis    = analysis,
    any.H       = 0L,
    size.H      = NA_integer_,
    size.Hc     = nrow(df),
    hr.H.true   = NA_real_, hr.H.hat   = NA_real_,
    hr.Hc.true  = NA_real_, hr.Hc.hat  = NA_real_,
    ahr.H.true  = NA_real_, ahr.Hc.true = NA_real_,
    ahr.H.hat   = NA_real_, ahr.Hc.hat  = NA_real_,
    cde.H.true  = NA_real_, cde.Hc.true = NA_real_,
    cde.H.hat   = NA_real_, cde.Hc.hat  = NA_real_,
    hr.itt      = NA_real_,
    hr.adj.itt  = NA_real_,
    ppv = NA_real_, npv = NA_real_, sens = NA_real_, spec = NA_real_,
    p.cens = 1 - mean(df[[event_name]]),
    taumax = max(df[[outcome_name]])
  )

  # ── Build estimation closures ─────────────────────────────────────────────
  if (is_glm) {
    glm_family <- switch(
      effect_measure %||% "OR",
      "OR"  = , "RD"  = stats::binomial(),
      "IRR" = , "IRD" = stats::poisson(),
      "MD"  =          stats::gaussian(),
      stats::binomial()
    )
    log_scale <- effect_measure %in% c("OR", "IRR", "RR")

    .glm_effect_grf <- function(sub_df, oc, ec, tc) {
      tryCatch({
        if (!is.null(offset_name) && offset_name %in% names(sub_df)) {
          off_vec <- log(pmax(sub_df[[offset_name]], 1e-6))
          fit <- stats::glm(sub_df[[ec]] ~ sub_df[[tc]],
                            family = glm_family, offset = off_vec)
        } else {
          fit <- stats::glm(sub_df[[ec]] ~ sub_df[[tc]],
                            family = glm_family)
        }
        b <- coef(fit)[[2]]
        if (log_scale) exp(b) else b
      }, error = function(e) NA_real_)
    }

    out$hr.itt <- .glm_effect_grf(df, outcome_name, event_name, treat_name)
  } else {
    out$hr.itt <- tryCatch(
      exp(survival::coxph(
        survival::Surv(df[[outcome_name]], df[[event_name]]) ~ df[[treat_name]]
      )$coefficients),
      error = function(e) NA_real_
    )
  }

  if (!is_glm && !is.null(cox_formula_adj)) {
    out$hr.adj.itt <- tryCatch(
      exp(survival::coxph(cox_formula_adj, data = df)$coefficients[treat_name]),
      error = function(e) NA_real_
    )
  }

  if (is.null(grf_est)) {
    out$hr.Hc.hat <- out$hr.itt
    return(out)
  }

  # GRF subgroup detection
  sg_harm_id <- grf_est$sg.harm.id %||% grf_est$sg_harm_id
  has_sg     <- !is.null(sg_harm_id) && length(sg_harm_id) > 0 &&
                !all(is.na(sg_harm_id))
  has_tr     <- !is.null(grf_est$data) &&
                "treat.recommend" %in% names(grf_est$data)

  if (!has_sg || !has_tr) {
    out$hr.Hc.hat <- out$hr.itt
    return(out)
  }

  out$any.H <- 1L
  grf_data  <- grf_est$data
  harm_ind  <- as.integer(grf_data$treat.recommend == 0)
  grf_data$sg_hat <- harm_ind

  out$size.H  <- sum(harm_ind)
  out$size.Hc <- sum(harm_ind == 0L)

  # Detect column names in grf_data (may differ from sim_data names)
  o_col <- if (outcome_name %in% names(grf_data)) outcome_name else
           if ("y.sim" %in% names(grf_data)) "y.sim" else outcome_name
  e_col <- if (event_name %in% names(grf_data)) event_name else
           if ("event.sim" %in% names(grf_data)) "event.sim" else event_name
  t_col <- if (treat_name %in% names(grf_data)) treat_name else
           if ("treat" %in% names(grf_data)) "treat" else treat_name

  if (is_glm) {
    if (out$size.H  > 10) out$hr.H.hat  <-
      .glm_effect_grf(grf_data[grf_data$sg_hat == 1, ], o_col, e_col, t_col)
    if (out$size.Hc > 10) out$hr.Hc.hat <-
      .glm_effect_grf(grf_data[grf_data$sg_hat == 0, ], o_col, e_col, t_col)
  } else {
    .cox_hr <- function(sub_df, oc, ec, tc) tryCatch(
      exp(survival::coxph(
        survival::Surv(sub_df[[oc]], sub_df[[ec]]) ~ sub_df[[tc]]
      )$coefficients),
      error = function(e) NA_real_
    )

    if (out$size.H  > 10) out$hr.H.hat  <-
      .cox_hr(grf_data[grf_data$sg_hat == 1, ], o_col, e_col, t_col)
    if (out$size.Hc > 10) out$hr.Hc.hat <-
      .cox_hr(grf_data[grf_data$sg_hat == 0, ], o_col, e_col, t_col)
  }

  # AHR / CDE (survival-only)
  if (!is_glm) {
    .merge_sg <- function(base_df, sg_df, cols) {
      if (id_name %in% names(base_df) && id_name %in% names(sg_df))
        merge(base_df[, c(id_name, cols), drop = FALSE],
              sg_df[, c(id_name, "sg_hat")], by = id_name, all.y = TRUE)
      else if (nrow(base_df) == nrow(sg_df)) {
        base_df$sg_hat <- sg_df$sg_hat; base_df
      } else NULL
    }

    if ("loghr_po" %in% names(df)) {
      m <- .merge_sg(df, grf_data, "loghr_po")
      if (!is.null(m)) {
        out$ahr.H.hat  <- compute_ahr(m, m$sg_hat)
        out$ahr.Hc.hat <- compute_ahr(m, 1L - m$sg_hat)
      }
    }
    if (all(c("theta_0", "theta_1") %in% names(df))) {
      m <- .merge_sg(df, grf_data, c("theta_0", "theta_1"))
      if (!is.null(m)) {
        out$cde.H.hat  <- compute_cde(m, m$sg_hat)
        out$cde.Hc.hat <- compute_cde(m, 1L - m$sg_hat)
      }
    }
  }

  # True-subgroup metrics
  if (harm_col %in% names(df)) {
    .merge_sg_class <- function(base_df, sg_df, hcol) {
      if (id_name %in% names(base_df) && id_name %in% names(sg_df))
        merge(base_df[, c(id_name, hcol), drop = FALSE],
              sg_df[, c(id_name, "sg_hat")], by = id_name, all.y = TRUE)
      else if (nrow(base_df) == nrow(sg_df)) {
        base_df$sg_hat <- sg_df$sg_hat; base_df
      } else NULL
    }

    m_class <- .merge_sg_class(df, grf_data, harm_col)
    if (!is.null(m_class)) {
      true_H <- m_class[[harm_col]] == 1L
      hat_H  <- m_class$sg_hat == 1L

      tp <- sum(true_H & hat_H);   fp <- sum(!true_H & hat_H)
      tn <- sum(!true_H & !hat_H); fn <- sum(true_H & !hat_H)

      out$sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
      out$spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
      out$ppv  <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
      out$npv  <- if ((tn + fn) > 0) tn / (tn + fn) else NA_real_
    }

    true_H_df <- df[[harm_col]] == 1L
    if (is_glm) {
      if (sum(true_H_df) > 10)  out$hr.H.true  <-
        .glm_effect_grf(df[true_H_df, ], outcome_name, event_name, treat_name)
      if (sum(!true_H_df) > 10) out$hr.Hc.true <-
        .glm_effect_grf(df[!true_H_df, ], outcome_name, event_name, treat_name)
    } else {
      if (sum(true_H_df) > 10)  out$hr.H.true  <-
        .cox_hr(df[true_H_df, ], outcome_name, event_name, treat_name)
      if (sum(!true_H_df) > 10) out$hr.Hc.true <-
        .cox_hr(df[!true_H_df, ], outcome_name, event_name, treat_name)

      if ("loghr_po" %in% names(df)) {
        out$ahr.H.true  <- compute_ahr(df, df[[harm_col]])
        out$ahr.Hc.true <- compute_ahr(df, 1L - df[[harm_col]])
      }
      if (all(c("theta_0", "theta_1") %in% names(df))) {
        out$cde.H.true  <- compute_cde(df, df[[harm_col]])
        out$cde.Hc.true <- compute_cde(df, 1L - df[[harm_col]])
      }
    }
  }

  out
}
