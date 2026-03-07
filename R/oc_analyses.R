# =============================================================================
# Operating Characteristics Analysis Functions
# =============================================================================
#
# Works with DGMs produced by generate_aft_dgm_flex().
#
# Canonical column convention (underscore): y_sim, event_sim, flag_harm,
# loghr_po, theta_0, theta_1.
#
# Key public functions:
#   - run_simulation_analysis()      Main simulation wrapper
#   - run_forestsearch_analysis()    ForestSearch analysis helper
#   - run_grf_analysis()             GRF analysis helper
#   - extract_fs_estimates()         Extract estimates from FS results
#   - extract_grf_estimates()        Extract estimates from GRF results
#   - summarize_simulation_results() Summary tables
#   - format_oc_results()            GT-formatted OC table
#   - compute_dgm_cde()              Attach CDE values to a DGM
#
# =============================================================================

#' @import survival
#' @import data.table
NULL


# =============================================================================
# Configuration Defaults
# =============================================================================

#' Default ForestSearch Parameters
#'
#' Returns a list of default parameters for ForestSearch analysis.
#' Column names use the canonical underscore convention (y_sim, event_sim).
#'
#' @return List of default ForestSearch parameters
#' @keywords internal
default_fs_params <- function() {
  list(
    outcome.name           = "y_sim",
    event.name             = "event_sim",
    treat.name             = "treat",
    id.name                = "id",
    use_lasso              = TRUE,
    use_grf                = FALSE,
    hr.threshold           = 1.25,
    hr.consistency         = 1.0,
    pconsistency.threshold = 0.90,
    fs.splits              = 400,
    n.min                  = 60,
    d0.min                 = 12,
    d1.min                 = 12,
    maxk                   = 2,
    max.minutes            = 5,
    by.risk                = 12,
    vi.grf.min             = -0.2,
    use_twostage           = FALSE,
    twostage_args          = list()
  )
}


#' Default GRF Parameters
#'
#' Returns a list of default parameters for GRF analysis.
#' Column names use the canonical underscore convention.
#'
#' @return List of default GRF parameters
#' @keywords internal
default_grf_params <- function() {
  list(
    outcome.name  = "y_sim",
    event.name    = "event_sim",
    treat.name    = "treat",
    id.name       = "id",
    n.min         = 60,
    dmin.grf      = 12,
    frac.tau      = 0.60,
    maxdepth      = 2,
    RCT           = TRUE,
    sg.criterion  = "mDiff",
    seedit        = 8316951
  )
}


# =============================================================================
# DGM Simulation Dispatcher
# =============================================================================

#' Dispatch Simulation
#'
#' Calls \code{simulate_from_dgm()} for \code{"aft_dgm_flex"} objects and
#' returns the resulting data frame.
#'
#' @param dgm A DGM object of class \code{"aft_dgm_flex"}
#' @param n Integer sample size
#' @param sim_id Integer simulation index (used for seed offset)
#' @param seed_base Integer base random seed. Default 8316951L
#' @param sim_args Named list of additional arguments forwarded to
#'   \code{simulate_from_dgm()} — commonly \code{analysis_time},
#'   \code{max_entry}, \code{cens_adjust}, \code{draw_treatment},
#'   \code{time_eos}, \code{rand_ratio}.
#' @param verbose Logical. Print dispatch details. Default FALSE
#'
#' @return Data frame from \code{simulate_from_dgm()}
#' @keywords internal
dispatch_simulate <- function(dgm,
                              n,
                              sim_id    = 1L,
                              seed_base = 8316951L,
                              sim_args  = list(),
                              verbose   = FALSE) {

  if (!inherits(dgm, "aft_dgm_flex"))
    stop("Unsupported DGM class: ", paste(class(dgm), collapse = ", "),
         "\nExpected 'aft_dgm_flex' from generate_aft_dgm_flex().")

  args <- list(dgm = dgm, n = n, seed = seed_base + sim_id)
  for (nm in names(sim_args)) {
    if (nm != "seed") args[[nm]] <- sim_args[[nm]]
  }

  if (verbose) {
    message(sprintf("  [dispatch] simulate_from_dgm(n=%d, seed=%d)",
                    n, args$seed))
    if ("analysis_time" %in% names(args))
      message(sprintf("  [dispatch]   analysis_time = %s", args$analysis_time))
    if ("cens_adjust" %in% names(args))
      message(sprintf("  [dispatch]   cens_adjust   = %.4f", args$cens_adjust))
  }

  tryCatch(
    do.call(simulate_from_dgm, args),
    error = function(e) stop("simulate_from_dgm() failed: ", e$message)
  )
}


# =============================================================================
# Helper Functions
# =============================================================================

#' Extract HR from DGM
#'
#' Extracts hazard ratios from an \code{aft_dgm_flex} object's
#' \code{hazard_ratios} list.
#'
#' @param dgm DGM object
#' @param which Character. Which HR to extract: \code{"hr_H"},
#'   \code{"hr_Hc"}, \code{"ahr_H"}, \code{"ahr_Hc"},
#'   \code{"hr_overall"}, \code{"ahr"}, \code{"cde_H"},
#'   \code{"cde_Hc"}, \code{"cde"}.
#' @return Numeric hazard ratio value
#' @keywords internal
get_dgm_hr <- function(dgm, which = "hr_H") {
  hr_list <- dgm$hazard_ratios
  if (is.null(hr_list)) return(NA_real_)
  result <- switch(which,
    "hr_H"       = hr_list$harm_subgroup,
    "hr_Hc"      = hr_list$no_harm_subgroup,
    "ahr_H"      = hr_list$AHR_harm,
    "ahr_Hc"     = hr_list$AHR_no_harm,
    "hr_overall" = hr_list$overall,
    "ahr"        = hr_list$AHR,
    "cde_H"      = hr_list$CDE_harm,
    "cde_Hc"     = hr_list$CDE_no_harm,
    "cde"        = hr_list$CDE,
    NA_real_
  )
  if (is.null(result)) NA_real_ else result
}


#' Resolve Super-Population Data Frame from a DGM
#'
#' Returns \code{dgm$df_super} or \code{NULL} if absent.
#'
#' @param dgm A DGM object
#' @return Data frame or NULL
#' @keywords internal
resolve_df_super <- function(dgm) {
  dgm$df_super
}


#' Compute and Attach CDE Values to a DGM Object
#'
#' Calculates Controlled Direct Effect (CDE) hazard ratios from the
#' super-population potential outcomes (\code{theta_0}, \code{theta_1})
#' and attaches them to the DGM's \code{hazard_ratios} list.
#'
#' The subgroup indicator column is auto-detected from \code{flag_harm}
#' or \code{H} (in that order of preference).
#'
#' @param dgm An \code{aft_dgm_flex} object from \code{\link{generate_aft_dgm_flex}}.
#' @param harm_col Character. Name of the subgroup indicator column.
#'   If \code{NULL} (default), auto-detected.
#' @return The DGM object with CDE values added to \code{dgm$hazard_ratios}
#'   and top-level fields \code{dgm$CDE}, \code{dgm$cde_H}, \code{dgm$cde_Hc}.
#' @seealso \code{\link{get_dgm_hr}}
#' @export
compute_dgm_cde <- function(dgm, harm_col = NULL) {

  df <- resolve_df_super(dgm)

  if (is.null(df)) {
    warning("DGM has no super-population data frame; cannot compute CDE.")
    return(dgm)
  }
  if (!all(c("theta_0", "theta_1") %in% names(df))) {
    warning("Super-population lacks theta_0/theta_1; cannot compute CDE.")
    return(dgm)
  }

  # Auto-detect subgroup indicator (canonical name first)
  if (is.null(harm_col)) {
    candidates <- c("flag_harm", "H")
    found <- intersect(candidates, names(df))
    if (length(found) == 0) {
      warning("No subgroup indicator found (tried: flag_harm, H).")
      cde_overall <- mean(exp(df$theta_1)) / mean(exp(df$theta_0))
      dgm$CDE <- cde_overall
      if (is.null(dgm$hazard_ratios)) dgm$hazard_ratios <- list()
      dgm$hazard_ratios$CDE <- cde_overall
      return(dgm)
    }
    harm_col <- found[1]
  }

  in_H <- df[[harm_col]] == 1

  cde_H       <- mean(exp(df$theta_1[in_H]))   / mean(exp(df$theta_0[in_H]))
  cde_Hc      <- mean(exp(df$theta_1[!in_H]))  / mean(exp(df$theta_0[!in_H]))
  cde_overall <- mean(exp(df$theta_1))          / mean(exp(df$theta_0))

  if (is.null(dgm$hazard_ratios)) dgm$hazard_ratios <- list()
  dgm$hazard_ratios$CDE         <- cde_overall
  dgm$hazard_ratios$CDE_harm    <- cde_H
  dgm$hazard_ratios$CDE_no_harm <- cde_Hc

  dgm$CDE    <- cde_overall
  dgm$cde_H  <- cde_H
  dgm$cde_Hc <- cde_Hc

  dgm
}


#' Create Subgroup Indicator from Factor Definitions
#'
#' @param df Data frame containing the variables
#' @param sg_factors Character vector of factor definitions (e.g. "v1.1")
#' @return Integer vector (1 = in subgroup, 0 = not)
#' @keywords internal
create_subgroup_indicator <- function(df, sg_factors) {
  indicator <- rep(TRUE, nrow(df))
  for (factor_def in sg_factors) {
    if (is.na(factor_def) || factor_def == "") next
    parts <- strsplit(factor_def, "\\.")[[1]]
    if (length(parts) >= 2) {
      var_name <- parts[1]
      level    <- parts[2]
      if (var_name %in% names(df)) {
        indicator <- indicator & (as.character(df[[var_name]]) == level)
      }
    }
  }
  as.integer(indicator)
}


#' Compute AHR from loghr_po
#' @keywords internal
compute_ahr <- function(df, subset_indicator = NULL) {
  if (!"loghr_po" %in% names(df)) return(NA_real_)
  loghr <- df$loghr_po
  if (!is.null(subset_indicator)) loghr <- loghr[subset_indicator == 1]
  if (length(loghr) == 0 || all(is.na(loghr))) return(NA_real_)
  exp(mean(loghr, na.rm = TRUE))
}


#' Compute CDE from theta_0 and theta_1
#' @keywords internal
compute_cde <- function(df, subset_indicator = NULL) {
  if (!all(c("theta_0", "theta_1") %in% names(df))) return(NA_real_)
  t0 <- df$theta_0
  t1 <- df$theta_1
  if (!is.null(subset_indicator)) {
    idx <- subset_indicator == 1
    t0  <- t0[idx]
    t1  <- t1[idx]
  }
  if (length(t0) == 0 || all(is.na(t0))) return(NA_real_)
  mean(exp(t1), na.rm = TRUE) / mean(exp(t0), na.rm = TRUE)
}


# =============================================================================
# Extract ForestSearch Estimates
# =============================================================================

#' Extract Estimates from ForestSearch Results
#'
#' Extracts operating characteristics from ForestSearch analysis results.
#' Expects the canonical underscore column convention in \code{df}:
#' \code{y_sim}, \code{event_sim}, \code{flag_harm}, \code{loghr_po},
#' \code{theta_0}, \code{theta_1}.
#'
#' @param df Simulated data frame (canonical column names)
#' @param fs_res ForestSearch result table or NULL
#' @param dgm DGM object
#' @param cox_formula Optional Cox formula for estimation
#' @param cox_formula_adj Optional adjusted Cox formula
#' @param analysis Analysis label (e.g. "FS", "FSlg")
#' @param fs_full Full forestsearch result object (for df.est access)
#' @param verbose Logical. Print extraction details. Default FALSE
#' @return data.table with extracted estimates
#' @importFrom data.table data.table
#' @importFrom survival coxph Surv
#' @keywords internal
extract_fs_estimates <- function(df,
                                 fs_res,
                                 dgm,
                                 cox_formula     = NULL,
                                 cox_formula_adj = NULL,
                                 analysis        = "FS",
                                 fs_full         = NULL,
                                 verbose         = FALSE) {

  out <- data.table::data.table(
    analysis    = analysis,
    any.H       = 0L,
    size.H      = NA_integer_,
    size.Hc     = nrow(df),
    hr.H.true   = NA_real_,
    hr.H.hat    = NA_real_,
    hr.Hc.true  = NA_real_,
    hr.Hc.hat   = NA_real_,
    ahr.H.true  = NA_real_,
    ahr.Hc.true = NA_real_,
    ahr.H.hat   = NA_real_,
    ahr.Hc.hat  = NA_real_,
    cde.H.true  = NA_real_,
    cde.Hc.true = NA_real_,
    cde.H.hat   = NA_real_,
    cde.Hc.hat  = NA_real_,
    hr.itt      = NA_real_,
    hr.adj.itt  = NA_real_,
    ppv         = NA_real_,
    npv         = NA_real_,
    sens        = NA_real_,
    spec        = NA_real_,
    p.cens      = 1 - mean(df$event_sim),
    taumax      = max(df$y_sim)
  )

  # ITT estimate
  out$hr.itt <- tryCatch(
    exp(survival::coxph(
      survival::Surv(y_sim, event_sim) ~ treat,
      data = df)$coefficients),
    error = function(e) NA_real_)

  # Adjusted ITT
  if (!is.null(cox_formula_adj)) {
    out$hr.adj.itt <- tryCatch(
      exp(survival::coxph(cox_formula_adj, data = df)$coefficients["treat"]),
      error = function(e) NA_real_)
  }

  # No subgroup found
  if (is.null(fs_res) || nrow(fs_res) == 0) {
    if (verbose) message(sprintf("  [%s] No subgroup identified", analysis))
    out$hr.Hc.hat <- out$hr.itt
    return(out)
  }

  out$any.H <- 1L

  # Resolve subgroup indicator
  df_with_sg <- if (!is.null(fs_full) && !is.null(fs_full$df.est))
                  fs_full$df.est else NULL

  if (!is.null(df_with_sg) && "treat.recommend" %in% names(df_with_sg)) {
    df$sg_hat <- as.integer(df_with_sg$treat.recommend == 0)
    if (verbose)
      message(sprintf("  [%s] Using treat.recommend from df.est", analysis))
  } else {
    sg_factors <- character(0)
    for (i in seq_len(min(7, ncol(fs_res)))) {
      col_name <- paste0("M.", i)
      if (col_name %in% names(fs_res) && !is.na(fs_res[[col_name]][1]))
        sg_factors <- c(sg_factors, fs_res[[col_name]][1])
    }
    if (length(sg_factors) > 0) {
      df$sg_hat <- create_subgroup_indicator(df, sg_factors)
      if (verbose)
        message(sprintf("  [%s] Parsed subgroup: %s",
                        analysis, paste(sg_factors, collapse = " & ")))
    } else {
      out$hr.Hc.hat <- out$hr.itt
      return(out)
    }
  }

  out$size.H  <- sum(df$sg_hat, na.rm = TRUE)
  out$size.Hc <- sum(!df$sg_hat, na.rm = TRUE)

  if (verbose)
    message(sprintf("  [%s] Subgroup: n_H = %d (%.1f%%), n_Hc = %d",
                    analysis, out$size.H,
                    100 * out$size.H / nrow(df), out$size.Hc))

  # Cox HRs in identified subgroups
  if (out$size.H > 10) {
    out$hr.H.hat <- tryCatch(
      exp(survival::coxph(
        survival::Surv(y_sim, event_sim) ~ treat,
        data = subset(df, sg_hat == 1))$coefficients),
      error = function(e) NA_real_)
  }
  if (out$size.Hc > 10) {
    out$hr.Hc.hat <- tryCatch(
      exp(survival::coxph(
        survival::Surv(y_sim, event_sim) ~ treat,
        data = subset(df, sg_hat == 0))$coefficients),
      error = function(e) NA_real_)
  }

  # AHR in identified subgroups
  if ("loghr_po" %in% names(df)) {
    out$ahr.H.hat   <- compute_ahr(df, df$sg_hat)
    out$ahr.Hc.hat  <- compute_ahr(df, 1L - df$sg_hat)
  }

  # CDE in identified subgroups
  if (all(c("theta_0", "theta_1") %in% names(df))) {
    out$cde.H.hat   <- compute_cde(df, df$sg_hat)
    out$cde.Hc.hat  <- compute_cde(df, 1L - df$sg_hat)
  }

  if (verbose)
    message(sprintf("  [%s] HR estimates: H = %.3f, Hc = %.3f",
                    analysis,
                    ifelse(is.na(out$hr.H.hat), NA, out$hr.H.hat),
                    ifelse(is.na(out$hr.Hc.hat), NA, out$hr.Hc.hat)))

  # Classification metrics (against canonical flag_harm)
  if ("flag_harm" %in% names(df)) {
    true_H <- df$flag_harm == 1
    hat_H  <- df$sg_hat == 1

    tp <- sum(true_H & hat_H)
    fp <- sum(!true_H & hat_H)
    tn <- sum(!true_H & !hat_H)
    fn <- sum(true_H & !hat_H)

    out$sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    out$spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
    out$ppv  <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    out$npv  <- if ((tn + fn) > 0) tn / (tn + fn) else NA_real_

    if (verbose)
      message(sprintf("  [%s] Classification: Sens=%.3f, Spec=%.3f, PPV=%.3f, NPV=%.3f",
                      analysis, out$sens, out$spec, out$ppv, out$npv))

    # Cox HRs in true subgroups
    if (out$size.H > 10)
      out$hr.H.true <- tryCatch(
        exp(survival::coxph(
          survival::Surv(y_sim, event_sim) ~ treat,
          data = subset(df, flag_harm == 1))$coefficients),
        error = function(e) NA_real_)

    if (out$size.Hc > 10)
      out$hr.Hc.true <- tryCatch(
        exp(survival::coxph(
          survival::Surv(y_sim, event_sim) ~ treat,
          data = subset(df, flag_harm == 0))$coefficients),
        error = function(e) NA_real_)

    # AHR in true subgroups
    if ("loghr_po" %in% names(df)) {
      out$ahr.H.true  <- compute_ahr(df, df$flag_harm)
      out$ahr.Hc.true <- compute_ahr(df, 1L - df$flag_harm)
    }
    # CDE in true subgroups
    if (all(c("theta_0", "theta_1") %in% names(df))) {
      out$cde.H.true  <- compute_cde(df, df$flag_harm)
      out$cde.Hc.true <- compute_cde(df, 1L - df$flag_harm)
    }
  }

  out
}


# =============================================================================
# Run ForestSearch Analysis
# =============================================================================

#' Run ForestSearch Analysis
#'
#' Helper function to run ForestSearch and extract estimates.
#'
#' @param data Data frame with simulated trial data (canonical column names)
#' @param confounders_name Character vector of confounder names
#' @param params List of ForestSearch parameters
#' @param dgm DGM object
#' @param cox_formula Cox formula for estimation
#' @param cox_formula_adj Adjusted Cox formula
#' @param analysis_label Character label for this analysis
#' @param verbose Print details
#' @return data.table with analysis estimates
#' @keywords internal
run_forestsearch_analysis <- function(data,
                                      confounders_name,
                                      params,
                                      dgm,
                                      cox_formula     = NULL,
                                      cox_formula_adj = NULL,
                                      analysis_label  = "FS",
                                      verbose         = FALSE) {

  if (verbose) {
    message(sprintf("\n  [%s] Starting ForestSearch...", analysis_label))
    message(sprintf("  [%s] n = %d, events = %d (%.1f%%)",
                    analysis_label, nrow(data), sum(data$event_sim),
                    100 * mean(data$event_sim)))
    message(sprintf("  [%s] Confounders: %s",
                    analysis_label, paste(confounders_name, collapse = ", ")))
    message(sprintf(
      "  [%s] use_lasso=%s, use_grf=%s, hr.threshold=%.2f, use_twostage=%s",
      analysis_label, params$use_lasso, params$use_grf,
      params$hr.threshold, params$use_twostage))
  }

  fs_args <- list(
    df.analysis      = data,
    confounders.name = confounders_name,
    details          = verbose,
    plot.sg          = FALSE
  )

  param_names <- c(
    "outcome.name", "event.name", "treat.name", "id.name",
    "use_lasso", "use_grf", "conf_force",
    "hr.threshold", "hr.consistency", "pconsistency.threshold",
    "fs.splits", "n.min", "d0.min", "d1.min",
    "maxk", "max.minutes", "by.risk", "vi.grf.min",
    "frac.tau", "dmin.grf", "grf_depth",
    "use_twostage", "twostage_args"
  )
  for (pn in param_names) {
    if (!is.null(params[[pn]])) fs_args[[pn]] <- params[[pn]]
  }

  if (verbose) message(sprintf("  [%s] Calling forestsearch()...", analysis_label))

  fs_result <- tryCatch(
    do.call(forestsearch, fs_args),
    error = function(e) {
      warning(sprintf("%s failed: %s", analysis_label, e$message))
      NULL
    })

  if (verbose) {
    if (is.null(fs_result)) {
      message(sprintf("  [%s] forestsearch() returned NULL", analysis_label))
    } else {
      sg_def <- if (!is.null(fs_result$sg.harm))
                  paste(fs_result$sg.harm, collapse = " & ") else "none"
      message(sprintf("  [%s] Subgroup: %s", analysis_label, sg_def))
    }
  }

  has_result <- !is.null(fs_result) &&
                !is.null(fs_result$grp.consistency) &&
                !is.null(fs_result$grp.consistency$out_sg) &&
                !is.null(fs_result$grp.consistency$out_sg$result) &&
                nrow(fs_result$grp.consistency$out_sg$result) > 0

  extract_fs_estimates(
    df              = data,
    fs_res          = if (has_result) fs_result$grp.consistency$out_sg$result else NULL,
    fs_full         = if (has_result) fs_result else NULL,
    dgm             = dgm,
    cox_formula     = cox_formula,
    cox_formula_adj = cox_formula_adj,
    analysis        = analysis_label,
    verbose         = verbose
  )
}


# =============================================================================
# Extract GRF Estimates
# =============================================================================

#' Extract Estimates from GRF Results
#'
#' Extracts operating characteristics from GRF analysis results.
#' Expects canonical underscore column names in \code{df}.
#'
#' @param df Simulated data frame (canonical column names)
#' @param grf_est GRF result from grf.subg.harm.survival()
#' @param dgm DGM object
#' @param cox_formula Cox formula
#' @param cox_formula_adj Adjusted Cox formula
#' @param analysis Analysis label
#' @param frac_tau Fraction of tau used
#' @param verbose Print extraction details
#' @param debug Print detailed GRF result structure
#' @return data.table with extracted estimates
#' @keywords internal
extract_grf_estimates <- function(df,
                                  grf_est,
                                  dgm,
                                  cox_formula     = NULL,
                                  cox_formula_adj = NULL,
                                  analysis        = "GRF",
                                  frac_tau        = 1.0,
                                  verbose         = FALSE,
                                  debug           = FALSE) {

  if (debug && !is.null(grf_est)) {
    message(sprintf("  [%s] DEBUG grf_est names: %s",
                    analysis, paste(names(grf_est), collapse = ", ")))
    message(sprintf("  [%s] DEBUG sg.harm.id = '%s'",
                    analysis, paste(grf_est$sg.harm.id, collapse = ", ")))
    if ("data" %in% names(grf_est) && "treat.recommend" %in% names(grf_est$data)) {
      tr <- grf_est$data$treat.recommend
      message(sprintf("  [%s] DEBUG treat.recommend: 0=%d, 1=%d",
                      analysis, sum(tr == 0, na.rm = TRUE),
                      sum(tr == 1, na.rm = TRUE)))
    }
  }

  out <- data.table::data.table(
    analysis    = analysis,
    any.H       = 0L,
    size.H      = NA_integer_,
    size.Hc     = nrow(df),
    hr.H.true   = NA_real_,
    hr.H.hat    = NA_real_,
    hr.Hc.true  = NA_real_,
    hr.Hc.hat   = NA_real_,
    ahr.H.true  = NA_real_,
    ahr.Hc.true = NA_real_,
    ahr.H.hat   = NA_real_,
    ahr.Hc.hat  = NA_real_,
    cde.H.true  = NA_real_,
    cde.Hc.true = NA_real_,
    cde.H.hat   = NA_real_,
    cde.Hc.hat  = NA_real_,
    hr.itt      = NA_real_,
    hr.adj.itt  = NA_real_,
    ppv         = NA_real_,
    npv         = NA_real_,
    sens        = NA_real_,
    spec        = NA_real_,
    p.cens      = 1 - mean(df$event_sim),
    taumax      = max(df$y_sim)
  )

  # ITT
  out$hr.itt <- tryCatch(
    exp(survival::coxph(
      survival::Surv(y_sim, event_sim) ~ treat,
      data = df)$coefficients),
    error = function(e) NA_real_)

  if (!is.null(cox_formula_adj)) {
    out$hr.adj.itt <- tryCatch(
      exp(survival::coxph(cox_formula_adj, data = df)$coefficients["treat"]),
      error = function(e) NA_real_)
  }

  if (is.null(grf_est)) {
    if (verbose)
      message(sprintf("  [%s] No GRF result - returning ITT", analysis))
    out$hr.Hc.hat <- out$hr.itt
    return(out)
  }

  # Extract subgroup definition
  sg_harm_id <- grf_est$sg.harm.id %||% grf_est$sg_harm_id

  has_sg_harm_id <- !is.null(sg_harm_id) &&
                    length(sg_harm_id) > 0 &&
                    !all(is.na(sg_harm_id)) &&
                    any(nchar(as.character(sg_harm_id)) > 0, na.rm = TRUE)

  has_treat_recommend <- !is.null(grf_est$data) &&
                         "treat.recommend" %in% names(grf_est$data)

  if (debug)
    message(sprintf("  [%s] DEBUG has_sg_harm_id=%s, has_treat_recommend=%s",
                    analysis, has_sg_harm_id, has_treat_recommend))

  if (!has_sg_harm_id || !has_treat_recommend) {
    if (verbose)
      message(sprintf("  [%s] No subgroup identified from GRF", analysis))
    out$hr.Hc.hat <- out$hr.itt
    return(out)
  }

  out$any.H <- 1L
  grf_data  <- grf_est$data

  # treat.recommend == 0 => harm subgroup
  harm_indicator <- as.integer(grf_data$treat.recommend == 0)
  out$size.H  <- sum(harm_indicator, na.rm = TRUE)
  out$size.Hc <- sum(harm_indicator == 0, na.rm = TRUE)

  if (verbose) {
    message(sprintf("  [%s] Subgroup: %s",
                    analysis, paste(sg_harm_id, collapse = " & ")))
    message(sprintf("  [%s] n_H=%d (%.1f%%), n_Hc=%d",
                    analysis, out$size.H,
                    100 * out$size.H / nrow(grf_data), out$size.Hc))
  }

  # Determine column names in grf_data (may use either convention)
  outcome_col <- if ("y_sim"    %in% names(grf_data)) "y_sim"    else
                 if ("y.sim"    %in% names(grf_data)) "y.sim"    else NULL
  event_col   <- if ("event_sim" %in% names(grf_data)) "event_sim" else
                 if ("event.sim" %in% names(grf_data)) "event.sim" else NULL
  treat_col   <- if ("treat" %in% names(grf_data)) "treat" else NULL

  if (is.null(outcome_col) || is.null(event_col) || is.null(treat_col)) {
    if (debug)
      message(sprintf("  [%s] DEBUG Missing cols in grf_data, falling back to df",
                      analysis))
    outcome_col <- "y_sim"
    event_col   <- "event_sim"
    treat_col   <- "treat"
    grf_data    <- df
  }
  grf_data$sg_hat <- harm_indicator

  # Cox HRs in identified subgroups
  cox_form_str <- sprintf("survival::Surv(%s, %s) ~ %s",
                          outcome_col, event_col, treat_col)

  if (out$size.H > 10 &&
      sum(grf_data$sg_hat == 1 & grf_data[[event_col]] == 1) >= 5) {
    out$hr.H.hat <- tryCatch({
      df_H <- grf_data[grf_data$sg_hat == 1, ]
      if (length(unique(df_H[[treat_col]])) == 2)
        exp(survival::coxph(as.formula(cox_form_str), data = df_H)$coefficients)
      else NA_real_
    }, error = function(e) NA_real_)
  }

  if (out$size.Hc > 10 &&
      sum(grf_data$sg_hat == 0 & grf_data[[event_col]] == 1) >= 5) {
    out$hr.Hc.hat <- tryCatch({
      df_Hc <- grf_data[grf_data$sg_hat == 0, ]
      if (length(unique(df_Hc[[treat_col]])) == 2)
        exp(survival::coxph(as.formula(cox_form_str), data = df_Hc)$coefficients)
      else NA_real_
    }, error = function(e) NA_real_)
  }

  if (verbose)
    message(sprintf("  [%s] HR: H=%.3f, Hc=%.3f",
                    analysis,
                    ifelse(is.na(out$hr.H.hat), NA, out$hr.H.hat),
                    ifelse(is.na(out$hr.Hc.hat), NA, out$hr.Hc.hat)))

  # AHR / CDE in identified subgroups (need to match grf_data rows to df rows)
  if ("loghr_po" %in% names(df)) {
    if ("id" %in% names(grf_data) && "id" %in% names(df)) {
      dm <- merge(df[, c("id", "loghr_po")],
                  grf_data[, c("id", "sg_hat")], by = "id", all.y = TRUE)
      out$ahr.H.hat   <- compute_ahr(dm, dm$sg_hat)
      out$ahr.Hc.hat  <- compute_ahr(dm, 1L - dm$sg_hat)
    } else if (nrow(df) == nrow(grf_data)) {
      df$sg_hat <- harm_indicator
      out$ahr.H.hat  <- compute_ahr(df, df$sg_hat)
      out$ahr.Hc.hat <- compute_ahr(df, 1L - df$sg_hat)
    }
  }

  if (all(c("theta_0", "theta_1") %in% names(df))) {
    if ("id" %in% names(grf_data) && "id" %in% names(df)) {
      dm2 <- merge(df[, c("id", "theta_0", "theta_1")],
                   grf_data[, c("id", "sg_hat")], by = "id", all.y = TRUE)
      out$cde.H.hat   <- compute_cde(dm2, dm2$sg_hat)
      out$cde.Hc.hat  <- compute_cde(dm2, 1L - dm2$sg_hat)
    } else if (nrow(df) == nrow(grf_data)) {
      if (!"sg_hat" %in% names(df)) df$sg_hat <- harm_indicator
      out$cde.H.hat  <- compute_cde(df, df$sg_hat)
      out$cde.Hc.hat <- compute_cde(df, 1L - df$sg_hat)
    }
  }

  # Classification metrics
  if ("flag_harm" %in% names(df)) {
    if ("id" %in% names(grf_data) && "id" %in% names(df)) {
      dm3     <- merge(df[, c("id", "flag_harm")],
                       grf_data[, c("id", "sg_hat")], by = "id")
      true_H  <- dm3$flag_harm == 1
      hat_H   <- dm3$sg_hat == 1
    } else if (nrow(df) == nrow(grf_data)) {
      true_H  <- df$flag_harm == 1
      hat_H   <- harm_indicator == 1
    } else {
      true_H  <- NULL
      hat_H   <- NULL
    }

    if (!is.null(true_H)) {
      tp <- sum(true_H & hat_H)
      fp <- sum(!true_H & hat_H)
      tn <- sum(!true_H & !hat_H)
      fn <- sum(true_H & !hat_H)

      out$sens <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
      out$spec <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
      out$ppv  <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
      out$npv  <- if ((tn + fn) > 0) tn / (tn + fn) else NA_real_

      if (verbose)
        message(sprintf("  [%s] Classification: Sens=%.3f Spec=%.3f PPV=%.3f NPV=%.3f",
                        analysis, out$sens, out$spec, out$ppv, out$npv))
    }

    if ("loghr_po" %in% names(df)) {
      out$ahr.H.true  <- compute_ahr(df, df$flag_harm)
      out$ahr.Hc.true <- compute_ahr(df, 1L - df$flag_harm)
    }
    if (all(c("theta_0", "theta_1") %in% names(df))) {
      out$cde.H.true  <- compute_cde(df, df$flag_harm)
      out$cde.Hc.true <- compute_cde(df, 1L - df$flag_harm)
    }
  }

  out
}


# =============================================================================
# Run GRF Analysis
# =============================================================================

#' Run Standalone GRF Analysis
#'
#' @param data Data frame with simulated trial data (canonical column names)
#' @param confounders_name Character vector of confounder names
#' @param params List of GRF parameters
#' @param dgm DGM object
#' @param cox_formula Cox formula
#' @param cox_formula_adj Adjusted Cox formula
#' @param analysis_label Character label
#' @param verbose Print details
#' @param debug Print detailed debugging
#' @return data.table with analysis estimates
#' @keywords internal
run_grf_analysis <- function(data,
                             confounders_name,
                             params,
                             dgm,
                             cox_formula     = NULL,
                             cox_formula_adj = NULL,
                             analysis_label  = "GRF",
                             verbose         = FALSE,
                             debug           = FALSE) {

  if (verbose)
    message(sprintf(
      "\n  [%s] GRF: n.min=%d, dmin.grf=%.1f, frac.tau=%.2f, maxdepth=%d",
      analysis_label, params$n.min, params$dmin.grf,
      params$frac.tau, params$maxdepth))

  grf_fun <- tryCatch(
    get("grf.subg.harm.survival", mode = "function", envir = parent.frame()),
    error = function(e) tryCatch(
      get("grf.subg.harm.survival", mode = "function", envir = globalenv()),
      error = function(e2) NULL))

  if (is.null(grf_fun)) {
    if (verbose)
      message(sprintf("  [%s] grf.subg.harm.survival not found. Skipping.",
                      analysis_label))
    warning("grf.subg.harm.survival not found. Skipping GRF analysis.")
    return(extract_grf_estimates(
      df = data, grf_est = NULL, dgm = dgm,
      cox_formula = cox_formula, cox_formula_adj = cox_formula_adj,
      analysis = analysis_label, verbose = verbose, debug = debug))
  }

  grf_args <- list(data = data, confounders.name = confounders_name,
                   details = verbose)

  grf_param_names <- c(
    "outcome.name", "event.name", "id.name", "treat.name",
    "frac.tau", "n.min", "dmin.grf", "RCT", "sg.criterion",
    "maxdepth", "seedit"
  )
  for (pn in grf_param_names) {
    if (!is.null(params[[pn]])) grf_args[[pn]] <- params[[pn]]
  }

  grf_result <- tryCatch(
    do.call(grf_fun, grf_args),
    error = function(e) {
      if (verbose)
        message(sprintf("  [%s] GRF failed: %s", analysis_label, e$message))
      warning(sprintf("%s failed: %s", analysis_label, e$message))
      NULL
    })

  extract_grf_estimates(
    df              = data,
    grf_est         = grf_result,
    dgm             = dgm,
    cox_formula     = cox_formula,
    cox_formula_adj = cox_formula_adj,
    analysis        = analysis_label,
    frac_tau        = params$frac.tau,
    verbose         = verbose,
    debug           = debug
  )
}


# =============================================================================
# Main Simulation Analysis Function
# =============================================================================

#' Run Single Simulation Analysis
#'
#' Executes ForestSearch and/or GRF analysis on a single simulated dataset
#' from a DGM produced by \code{\link{generate_aft_dgm_flex}}.
#'
#' @param sim_id Integer. Simulation index (used for seed offset and progress)
#' @param dgm A DGM object of class \code{"aft_dgm_flex"}
#' @param n_sample Integer. Sample size for simulation
#' @param confounders_base Character vector. Confounder names presented to the
#'   analysis algorithms.
#' @param n_add_noise Integer. Number of standard-normal noise variables to
#'   add to \code{confounders_base}. Default 0L
#' @param run_fs Logical. Run ForestSearch with LASSO selection. Default TRUE
#' @param run_fs_grf Logical. Run ForestSearch with LASSO+GRF selection.
#'   Default TRUE
#' @param run_grf Logical. Run standalone GRF analysis. Default TRUE
#' @param fs_params List. ForestSearch parameter overrides (merged on top of
#'   \code{default_fs_params()}; user values always win)
#' @param grf_params List. GRF parameter overrides (merged on top of
#'   \code{default_grf_params()})
#' @param sim_args Named list of additional arguments forwarded verbatim to
#'   \code{\link{simulate_from_dgm}}. Commonly: \code{analysis_time},
#'   \code{max_entry}, \code{cens_adjust}, \code{rand_ratio},
#'   \code{draw_treatment}, \code{time_eos}. Default: \code{list()}
#' @param cox_formula Optional Cox formula for ITT estimation
#' @param cox_formula_adj Optional adjusted Cox formula
#' @param n_sims_total Integer. Total simulations (for progress messages)
#' @param seed_base Integer. Base random seed. Default 8316951L
#' @param verbose Logical. Print progress. Default FALSE
#' @param verbose_n Integer. Limit verbose output to the first N sims.
#'   Default NULL (all sims when verbose = TRUE)
#' @param debug Logical. Print detailed GRF debugging. Default FALSE
#'
#' @return A data.table with one row per analysis method run, containing
#'   simulation ID, true subgroup properties, and all OC metrics.
#'
#' @details
#' ## Analysis methods
#'
#' Three methods can be run:
#' \itemize{
#'   \item \strong{FS}: ForestSearch with LASSO only
#'     (use_lasso=TRUE, use_grf=FALSE)
#'   \item \strong{FSlg}: ForestSearch with LASSO+GRF
#'     (use_lasso=TRUE, use_grf=TRUE)
#'   \item \strong{GRF}: Standalone GRF via grf.subg.harm.survival()
#' }
#'
#' ## Parameter merging
#'
#' For FS/FSlg: defaults -> analysis-type-specific defaults ->
#' user's \code{fs_params} (user always wins).
#'
#' @importFrom data.table data.table rbindlist
#' @importFrom stats rnorm
#' @export
run_simulation_analysis <- function(
    sim_id,
    dgm,
    n_sample,
    confounders_base,
    n_add_noise     = 0L,
    run_fs          = TRUE,
    run_fs_grf      = TRUE,
    run_grf         = TRUE,
    fs_params       = list(),
    grf_params      = list(),
    sim_args        = list(),
    cox_formula     = NULL,
    cox_formula_adj = NULL,
    n_sims_total    = NULL,
    seed_base       = 8316951L,
    verbose         = FALSE,
    verbose_n       = NULL,
    debug           = FALSE
) {

  # Effective verbosity
  show_verbose <- verbose && (is.null(verbose_n) || sim_id <= verbose_n)

  if (show_verbose) {
    message("\n", paste(rep("=", 60), collapse = ""))
    message(sprintf("Simulation %d", sim_id))
    if (!is.null(n_sims_total))
      message(sprintf("  Progress: %d / %d (%.1f%%)",
                      sim_id, n_sims_total, 100 * sim_id / n_sims_total))
    message(sprintf("  DGM class: %s", paste(class(dgm), collapse = ", ")))
    message(paste(rep("=", 60), collapse = ""))
  }

  # ---------------------------------------------------------------------------
  # Simulate Data
  # ---------------------------------------------------------------------------
  if (show_verbose)
    message(sprintf("\n[1] Simulating data (n=%d)...", n_sample))

  sim_data <- dispatch_simulate(
    dgm       = dgm,
    n         = n_sample,
    sim_id    = sim_id,
    seed_base = seed_base,
    sim_args  = sim_args,
    verbose   = show_verbose
  )

  if (show_verbose)
    message(sprintf("    Simulated: n=%d, events=%d (%.1f%%)",
                    nrow(sim_data), sum(sim_data$event_sim),
                    100 * mean(sim_data$event_sim)))

  # ---------------------------------------------------------------------------
  # Add Noise Variables
  # ---------------------------------------------------------------------------
  confounders_name <- confounders_base

  if (n_add_noise > 0) {
    set.seed(seed_base + 1000L * sim_id)
    noise_names <- paste0("noise", seq_len(n_add_noise))
    for (nm in noise_names) sim_data[[nm]] <- stats::rnorm(nrow(sim_data))
    confounders_name <- c(confounders_base, noise_names)
    if (show_verbose)
      message(sprintf("    Added %d noise variables", n_add_noise))
  }

  # ---------------------------------------------------------------------------
  # True Subgroup Properties
  # ---------------------------------------------------------------------------
  size_H_true  <- if ("flag_harm" %in% names(sim_data)) sum(sim_data$flag_harm)   else NA_integer_
  prop_H_true  <- if ("flag_harm" %in% names(sim_data)) mean(sim_data$flag_harm)  else NA_real_
  size_Hc_true <- if ("flag_harm" %in% names(sim_data)) sum(!sim_data$flag_harm)  else NA_integer_
  prop_Hc_true <- if ("flag_harm" %in% names(sim_data)) mean(!sim_data$flag_harm) else NA_real_

  hr_H_dgm   <- get_dgm_hr(dgm, "hr_H")
  hr_Hc_dgm  <- get_dgm_hr(dgm, "hr_Hc")
  ahr_H_dgm  <- get_dgm_hr(dgm, "ahr_H")
  ahr_Hc_dgm <- get_dgm_hr(dgm, "ahr_Hc")

  if (show_verbose) {
    message("\n[2] True subgroup properties:")
    message(sprintf("    H:  n=%d (%.1f%%)", size_H_true,  100 * prop_H_true))
    message(sprintf("    Hc: n=%d (%.1f%%)", size_Hc_true, 100 * prop_Hc_true))
    message(sprintf("    DGM HR_H=%.3f, HR_Hc=%.3f", hr_H_dgm, hr_Hc_dgm))
    if (!is.na(ahr_H_dgm))
      message(sprintf("    DGM AHR_H=%.3f, AHR_Hc=%.3f", ahr_H_dgm, ahr_Hc_dgm))
  }

  df_pop <- data.table::data.table(
    sim          = sim_id,
    sizeH_true   = size_H_true,
    propH_true   = prop_H_true,
    sizeHc_true  = size_Hc_true,
    propHc_true  = prop_Hc_true
  )

  # ---------------------------------------------------------------------------
  # Merge Parameters with Defaults
  # ---------------------------------------------------------------------------
  fs_defaults  <- default_fs_params()
  grf_defaults <- default_grf_params()
  grf_merged   <- modifyList(grf_defaults, grf_params)
  # Ensure GRF uses a unique seed per simulation so that forest splits are
  # independently randomised.  A fixed seedit (the default) would make every
  # GRF call use identical random partitions, inflating or deflating detection
  # rates in a non-stochastic way across simulations.
  grf_merged$seedit <- seed_base + sim_id

  # ---------------------------------------------------------------------------
  # Run Analyses
  # ---------------------------------------------------------------------------
  results_list <- list()

  if (show_verbose) {
    to_run <- c(
      if (run_fs)      "FS (LASSO)"       else NULL,
      if (run_fs_grf)  "FSlg (LASSO+GRF)" else NULL,
      if (run_grf)     "GRF (standalone)"  else NULL
    )
    message(sprintf("\n[3] Running: %s", paste(to_run, collapse = ", ")))
  }

  # FS: LASSO only
  if (run_fs) {
    # Apply user overrides first, then enforce method-level flags last so that
    # fs_params$use_grf can never accidentally make FS run as FSlg.
    params_fs             <- modifyList(fs_defaults, fs_params)
    params_fs$use_lasso   <- TRUE
    params_fs$use_grf     <- FALSE   # always forced: FS = LASSO only

    if (show_verbose)
      message(sprintf("    [FS] use_lasso=%s use_grf=%s use_twostage=%s",
                      params_fs$use_lasso, params_fs$use_grf,
                      params_fs$use_twostage))

    results_list[["FS"]] <- cbind(df_pop, run_forestsearch_analysis(
      data            = sim_data,
      confounders_name = confounders_name,
      params          = params_fs,
      dgm             = dgm,
      cox_formula     = cox_formula,
      cox_formula_adj = cox_formula_adj,
      analysis_label  = "FS",
      verbose         = show_verbose
    ))
  }

  # FSlg: LASSO + GRF
  if (run_fs_grf) {
    params_fslg             <- modifyList(fs_defaults, fs_params)
    params_fslg$use_lasso   <- TRUE
    params_fslg$use_grf     <- TRUE   # always forced: FSlg = LASSO + GRF

    if (show_verbose)
      message(sprintf("    [FSlg] use_lasso=%s use_grf=%s use_twostage=%s",
                      params_fslg$use_lasso, params_fslg$use_grf,
                      params_fslg$use_twostage))

    results_list[["FSlg"]] <- cbind(df_pop, run_forestsearch_analysis(
      data             = sim_data,
      confounders_name = confounders_name,
      params           = params_fslg,
      dgm              = dgm,
      cox_formula      = cox_formula,
      cox_formula_adj  = cox_formula_adj,
      analysis_label   = "FSlg",
      verbose          = show_verbose
    ))
  }

  # GRF: standalone
  if (run_grf) {
    if (show_verbose)
      message(sprintf("    [GRF] n.min=%d dmin.grf=%.1f frac.tau=%.2f maxdepth=%d",
                      grf_merged$n.min, grf_merged$dmin.grf,
                      grf_merged$frac.tau, grf_merged$maxdepth))

    results_list[["GRF"]] <- cbind(df_pop, run_grf_analysis(
      data             = sim_data,
      confounders_name = confounders_name,
      params           = grf_merged,
      dgm              = dgm,
      cox_formula      = cox_formula,
      cox_formula_adj  = cox_formula_adj,
      analysis_label   = "GRF",
      verbose          = show_verbose,
      debug            = debug
    ))
  }

  if (length(results_list) == 0) {
    warning("No analyses were run. Check run_fs, run_fs_grf, run_grf settings.")
    return(NULL)
  }

  result <- data.table::rbindlist(results_list, fill = TRUE)

  if (show_verbose) {
    message("\n[4] Done.")
    message(sprintf("    Results: %d rows x %d columns",
                    nrow(result), ncol(result)))
    message(paste(rep("=", 60), collapse = ""), "\n")
  }

  result
}


# =============================================================================
# Summary Functions
# =============================================================================

#' Summarize Simulation Results
#'
#' Creates a summary table of operating characteristics across all simulations.
#' Includes Cox HR, AHR, and CDE metrics.
#'
#' @param results data.table from \code{\link{run_simulation_analysis}}
#' @param analyses Character vector. Analysis methods to include. Default: all
#' @param digits Integer. Decimal places for proportions. Default 2
#' @param digits_hr Integer. Decimal places for hazard ratios. Default 3
#' @return Data frame with summary statistics
#' @export
summarize_simulation_results <- function(results,
                                         analyses  = NULL,
                                         digits    = 2,
                                         digits_hr = 3) {
  if (is.null(analyses)) analyses <- unique(results$analysis)

  summaries <- lapply(analyses, function(a) {
    res <- subset(results, analysis == a)
    summarize_single_analysis(res, digits = digits, digits_hr = digits_hr)
  })

  summary_df <- do.call(cbind, summaries)
  colnames(summary_df) <- analyses
  summary_df
}


#' Summarize Single Analysis Results
#' @keywords internal
summarize_single_analysis <- function(result, digits = 2, digits_hr = 3) {

  class_cols  <- intersect(c("any.H", "sens", "spec", "ppv", "npv"), names(result))
  class_means <- round(sapply(result[, class_cols, with = FALSE],
                              mean, na.rm = TRUE), digits)

  res_found <- subset(result, any.H == 1)

  if (nrow(res_found) > 0) {
    avg_H <- round(mean(res_found$size.H, na.rm = TRUE), 0)
    min_H <- round(min(res_found$size.H,  na.rm = TRUE), 0)
    max_H <- round(max(res_found$size.H,  na.rm = TRUE), 0)
  } else {
    avg_H <- min_H <- max_H <- NA
  }

  avg_Hc <- round(mean(result$size.Hc, na.rm = TRUE), 0)
  min_Hc <- round(min(result$size.Hc,  na.rm = TRUE), 0)
  max_Hc <- round(max(result$size.Hc,  na.rm = TRUE), 0)

  safe_mean <- function(x, ...) if (length(x) == 0) NA else mean(x, ...)

  if (nrow(res_found) > 0) {
    hr_H_true   <- round(safe_mean(res_found$hr.H.true,   na.rm = TRUE), digits_hr)
    hr_H_hat    <- round(safe_mean(res_found$hr.H.hat,    na.rm = TRUE), digits_hr)
    hr_Hc_true  <- round(safe_mean(res_found$hr.Hc.true,  na.rm = TRUE), digits_hr)
    hr_Hc_hat   <- round(safe_mean(res_found$hr.Hc.hat,   na.rm = TRUE), digits_hr)
    ahr_H_true  <- round(safe_mean(res_found$ahr.H.true,  na.rm = TRUE), digits_hr)
    ahr_Hc_true <- round(safe_mean(res_found$ahr.Hc.true, na.rm = TRUE), digits_hr)
    ahr_H_hat   <- round(safe_mean(res_found$ahr.H.hat,   na.rm = TRUE), digits_hr)
    ahr_Hc_hat  <- round(safe_mean(res_found$ahr.Hc.hat,  na.rm = TRUE), digits_hr)
  } else {
    hr_H_true <- hr_H_hat <- hr_Hc_true <- hr_Hc_hat <- NA
    ahr_H_true <- ahr_Hc_true <- ahr_H_hat <- ahr_Hc_hat <- NA
  }

  hr_H_true_all  <- round(safe_mean(result$hr.H.true,   na.rm = TRUE), digits_hr)
  hr_Hc_true_all <- round(safe_mean(result$hr.Hc.true,  na.rm = TRUE), digits_hr)
  hr_itt         <- round(safe_mean(result$hr.itt,       na.rm = TRUE), digits_hr)
  hr_adj_itt     <- if ("hr.adj.itt" %in% names(result))
                      round(safe_mean(result$hr.adj.itt, na.rm = TRUE), digits_hr) else NA

  ahr_H_true_all  <- if ("ahr.H.true"  %in% names(result))
                       round(safe_mean(result$ahr.H.true,  na.rm = TRUE), digits_hr) else NA
  ahr_Hc_true_all <- if ("ahr.Hc.true" %in% names(result))
                       round(safe_mean(result$ahr.Hc.true, na.rm = TRUE), digits_hr) else NA

  values <- c(
    class_means,
    avg_H, min_H, max_H, avg_Hc, min_Hc, max_Hc,
    hr_H_true,  hr_H_hat,  hr_Hc_true,  hr_Hc_hat,
    hr_H_true_all, hr_Hc_true_all, hr_itt, hr_adj_itt,
    ahr_H_true, ahr_H_hat, ahr_Hc_true, ahr_Hc_hat,
    ahr_H_true_all, ahr_Hc_true_all
  )

  row_names <- c(
    names(class_means),
    "Avg(#H)", "minH", "maxH", "Avg(#Hc)", "minHc", "maxHc",
    "theta-dag(H)", "theta-hat(H-hat)", "theta-dag(Hc)", "theta-hat(Hc-hat)",
    "theta-dag(H)all", "theta-dag(Hc)all", "theta-hat(ITT)", "theta-hat(ITTadj)",
    "ahr-dag(H)", "ahr-hat(H-hat)", "ahr-dag(Hc)", "ahr-hat(Hc-hat)",
    "ahr-dag(H)all", "ahr-dag(Hc)all"
  )

  out <- data.frame(value = values, stringsAsFactors = FALSE)
  rownames(out) <- row_names
  out
}


# =============================================================================
# Format Operating Characteristics Table
# =============================================================================

#' Format Operating Characteristics Results as GT Table
#'
#' Creates a formatted gt table from simulation operating characteristics.
#'
#' @param results data.table or data.frame from
#'   \code{\link{run_simulation_analysis}}
#' @param analyses Character vector. Analysis methods to include.
#'   Default NULL (all)
#' @param metrics Character vector. Metrics to display. Options:
#'   \code{"detection"}, \code{"classification"}, \code{"hr_estimates"},
#'   \code{"ahr_estimates"}, \code{"cde_estimates"}, \code{"subgroup_size"},
#'   \code{"all"} (default)
#' @param digits Integer. Decimal places for proportions. Default 3
#' @param digits_hr Integer. Decimal places for hazard ratios. Default 3
#' @param title Character. Table title
#' @param subtitle Character. Table subtitle
#' @param use_gt Logical. Return gt table if TRUE. Default TRUE
#' @return A gt table or data.frame
#' @importFrom data.table is.data.table as.data.table
#' @export
format_oc_results <- function(results,
                              analyses  = NULL,
                              metrics   = "all",
                              digits    = 3,
                              digits_hr = 3,
                              title     = "Operating Characteristics Summary",
                              subtitle  = NULL,
                              use_gt    = TRUE) {

  if (!data.table::is.data.table(results))
    results <- data.table::as.data.table(results)

  if (is.null(analyses)) analyses <- unique(results$analysis)
  results <- results[results$analysis %in% analyses, ]

  summary_list <- lapply(analyses, function(a) {
    res       <- results[results$analysis == a, ]
    n_sims    <- nrow(res)
    res_found <- res[res$any.H == 1, ]

    # Detection and classification
    detection_rate <- mean(res$any.H,  na.rm = TRUE)
    sens <- mean(res$sens, na.rm = TRUE)
    spec <- mean(res$spec, na.rm = TRUE)
    ppv  <- mean(res$ppv,  na.rm = TRUE)
    npv  <- mean(res$npv,  na.rm = TRUE)

    # Conditional estimates (when subgroup found)
    cond <- function(col) {
      if (nrow(res_found) > 0 && col %in% names(res_found))
        mean(res_found[[col]], na.rm = TRUE) else NA
    }

    hr_itt <- mean(res$hr.itt, na.rm = TRUE)

    data.frame(
      Analysis     = a,
      N_sims       = n_sims,
      Detection    = detection_rate,
      Sen          = sens, Spec = spec, PPV = ppv, NPV = npv,
      HR_H_hat     = cond("hr.H.hat"),
      HR_Hc_hat    = cond("hr.Hc.hat"),
      HR_H_true    = cond("hr.H.true"),
      HR_Hc_true   = cond("hr.Hc.true"),
      HR_ITT       = hr_itt,
      AHR_H_true   = cond("ahr.H.true"),
      AHR_Hc_true  = cond("ahr.Hc.true"),
      AHR_H_hat    = cond("ahr.H.hat"),
      AHR_Hc_hat   = cond("ahr.Hc.hat"),
      CDE_H_true   = cond("cde.H.true"),
      CDE_Hc_true  = cond("cde.Hc.true"),
      CDE_H_hat    = cond("cde.H.hat"),
      CDE_Hc_hat   = cond("cde.Hc.hat"),
      Size_H_mean  = if (nrow(res_found) > 0) mean(res_found$size.H,  na.rm = TRUE) else NA,
      Size_H_min   = if (nrow(res_found) > 0) min(res_found$size.H,   na.rm = TRUE) else NA,
      Size_H_max   = if (nrow(res_found) > 0) max(res_found$size.H,   na.rm = TRUE) else NA,
      stringsAsFactors = FALSE
    )
  })

  summary_df <- do.call(rbind, summary_list)

  # Column filtering by metrics
  if (!"all" %in% metrics) {
    cols_to_keep <- c("Analysis", "N_sims")
    if ("detection"       %in% metrics) cols_to_keep <- c(cols_to_keep, "Detection")
    if ("classification"  %in% metrics) cols_to_keep <- c(cols_to_keep, "Sen","Spec","PPV","NPV")
    if ("hr_estimates"    %in% metrics) cols_to_keep <- c(cols_to_keep,
      "HR_H_hat","HR_Hc_hat","HR_H_true","HR_Hc_true","HR_ITT")
    if ("ahr_estimates"   %in% metrics) cols_to_keep <- c(cols_to_keep,
      "AHR_H_true","AHR_Hc_true","AHR_H_hat","AHR_Hc_hat")
    if ("cde_estimates"   %in% metrics) cols_to_keep <- c(cols_to_keep,
      "CDE_H_true","CDE_Hc_true","CDE_H_hat","CDE_Hc_hat")
    if ("subgroup_size"   %in% metrics) cols_to_keep <- c(cols_to_keep,
      "Size_H_mean","Size_H_min","Size_H_max")
    summary_df <- summary_df[, intersect(cols_to_keep, names(summary_df)), drop = FALSE]
  }

  if (!use_gt || !requireNamespace("gt", quietly = TRUE))
    return(summary_df)

  # Build GT table
  gt_table <- gt::gt(summary_df)
  gt_table <- gt::tab_header(gt_table, title = title, subtitle = subtitle)

  # Format proportions
  prop_cols <- intersect(c("Detection","Sen","Spec","PPV","NPV"), names(summary_df))
  if (length(prop_cols) > 0)
    gt_table <- gt::fmt_number(gt_table,
                               columns  = gt::all_of(prop_cols),
                               decimals = digits)

  # Format HR columns
  hr_cols <- intersect(c("HR_H_hat","HR_Hc_hat","HR_H_true","HR_Hc_true","HR_ITT"),
                       names(summary_df))
  if (length(hr_cols) > 0)
    gt_table <- gt::fmt_number(gt_table,
                               columns  = gt::all_of(hr_cols),
                               decimals = digits_hr)

  # Format AHR columns
  ahr_cols <- intersect(c("AHR_H_true","AHR_Hc_true","AHR_H_hat","AHR_Hc_hat"),
                        names(summary_df))
  if (length(ahr_cols) > 0)
    gt_table <- gt::fmt_number(gt_table,
                               columns  = gt::all_of(ahr_cols),
                               decimals = digits_hr)

  # Format CDE columns
  cde_cols <- intersect(c("CDE_H_true","CDE_Hc_true","CDE_H_hat","CDE_Hc_hat"),
                        names(summary_df))
  if (length(cde_cols) > 0)
    gt_table <- gt::fmt_number(gt_table,
                               columns  = gt::all_of(cde_cols),
                               decimals = digits_hr)

  # Format size columns
  size_cols <- intersect(c("Size_H_mean","Size_H_min","Size_H_max"), names(summary_df))
  if (length(size_cols) > 0)
    gt_table <- gt::fmt_number(gt_table,
                               columns  = gt::all_of(size_cols),
                               decimals = 0)

  # Column labels (Unicode math notation from Leon et al. 2024)
  label_list <- list(Analysis = "Method", N_sims = "Sims", Detection = "Found")

  map_labels <- function(nm_map) {
    for (col in names(nm_map)) {
      if (col %in% names(summary_df)) label_list[[col]] <<- nm_map[[col]]
    }
  }
  map_labels(list(
    HR_H_hat   = "\u03b8\u0302(\u0124)",
    HR_Hc_hat  = "\u03b8\u0302(\u0124\u1d9c)",
    HR_H_true  = "\u03b8\u0302(H)",
    HR_Hc_true = "\u03b8\u0302(H\u1d9c)",
    HR_ITT     = "\u03b8\u0302(ITT)"
  ))
  map_labels(list(
    AHR_H_true  = "\u00e2hr(H)",
    AHR_Hc_true = "\u00e2hr(H\u1d9c)",
    AHR_H_hat   = "\u00e2hr(\u0124)",
    AHR_Hc_hat  = "\u00e2hr(\u0124\u1d9c)"
  ))
  map_labels(list(
    CDE_H_true  = "\u03b8\u2021(H)",
    CDE_Hc_true = "\u03b8\u2021(H\u1d9c)",
    CDE_H_hat   = "\u03b8\u2021(\u0124)",
    CDE_Hc_hat  = "\u03b8\u2021(\u0124\u1d9c)"
  ))
  gt_table <- gt::cols_label(gt_table, .list = label_list)

  # Column spanners
  add_spanner <- function(tbl, label, cols) {
    present <- intersect(cols, names(summary_df))
    if (length(present) > 0)
      gt::tab_spanner(tbl, label = label, columns = gt::all_of(present))
    else
      tbl
  }

  if ("all" %in% metrics || "classification"  %in% metrics)
    gt_table <- add_spanner(gt_table, "Classification",
                            c("Sen","Spec","PPV","NPV"))
  if ("all" %in% metrics || "hr_estimates"    %in% metrics)
    gt_table <- add_spanner(gt_table, "Cox Hazard Ratios",
                            c("HR_H_hat","HR_Hc_hat","HR_H_true","HR_Hc_true","HR_ITT"))
  if ("all" %in% metrics || "ahr_estimates"   %in% metrics)
    gt_table <- add_spanner(gt_table, "Average Hazard Ratios",
                            c("AHR_H_true","AHR_Hc_true","AHR_H_hat","AHR_Hc_hat"))
  if ("all" %in% metrics || "cde_estimates"   %in% metrics)
    gt_table <- add_spanner(gt_table, "Controlled Direct Effects",
                            c("CDE_H_true","CDE_Hc_true","CDE_H_hat","CDE_Hc_hat"))

  gt_table <- gt::tab_style(
    gt_table,
    style     = gt::cell_text(weight = "bold"),
    locations = gt::cells_column_labels())

  if (utils::packageVersion("gt") >= "0.6.0")
    gt_table <- gt::sub_missing(gt_table,
                                columns      = gt::everything(),
                                missing_text = "-")
  gt_table
}


# =============================================================================
# Utility: null-coalescing operator (internal)
# =============================================================================
`%||%` <- function(a, b) if (!is.null(a)) a else b
