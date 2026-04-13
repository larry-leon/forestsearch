# =============================================================================
# DGM and Simulation Reporting Utilities
# =============================================================================
#
# Utility functions for extracting true effects from DGM objects,
# computing subgroup-level estimands, and formatting operating
# characteristics tables.
#
# These are internal helpers used by:
#   - run_simulation_analysis.R (effect extraction, classification)
#   - simulation_tables.R (build_estimation_table, format_oc_results)
#   - setup_gbsg_dgm.R (compute_dgm_cde)
#   - generate_glm_dgm.R (get_dgm_hr compatibility)
#
# Key functions:
#   - get_dgm_hr()                  DGM effect accessor
#   - compute_dgm_cde()             Causal direct effect from DGM
#   - create_subgroup_indicator()   Build subgroup flag from cut expressions
#   - compute_ahr()                 Average hazard ratio (survival)
#   - compute_cde()                 Causal direct effect (survival)
#   - summarize_simulation_results() Summary table from stacked results
#   - summarize_single_analysis()   Single-analysis summary
#   - format_oc_results()           Operating characteristics gt table
# =============================================================================

# Null-coalescing operator (used by compute_dgm_cde)
`%||%` <- function(a, b) if (!is.null(a)) a else b



# =============================================================================
# Helper Functions
# =============================================================================

#' Extract HR from DGM (Backward Compatible)
#'
#' Extracts hazard ratios from DGM object, supporting both old and new formats.
#' Also supports CDE (controlled direct effect) extraction for Table 5 of
#' Leon et al. (2024) alignment (theta-ddagger).
#'
#' @param dgm DGM object (gbsg_dgm or aft_dgm_flex)
#' @param which Character. Which HR to extract: \code{"hr_H"}, \code{"hr_Hc"},
#'   \code{"ahr_H"}, \code{"ahr_Hc"}, \code{"hr_overall"}, \code{"ahr"},
#'   \code{"cde_H"}, \code{"cde_Hc"}, \code{"cde"}.
#'
#' @return Numeric hazard ratio value
#'
#' @keywords internal
get_dgm_hr <- function(dgm, which = "hr_H") {

  # Try new hazard_ratios list first (aligned format)
  if (!is.null(dgm$hazard_ratios)) {
    hr_list <- dgm$hazard_ratios
    result <- switch(which,
      "hr_H" = hr_list$harm_subgroup,
      "hr_Hc" = hr_list$no_harm_subgroup,
      "ahr_H" = hr_list$AHR_harm,
      "ahr_Hc" = hr_list$AHR_no_harm,
      "hr_overall" = hr_list$overall,
      "ahr" = hr_list$AHR,
      "cde_H" = hr_list$CDE_harm,
      "cde_Hc" = hr_list$CDE_no_harm,
      "cde" = hr_list$CDE,
      NA_real_
    )
    if (!is.null(result) && !is.na(result)) return(result)
  }


  # Fall back to old format (direct properties)
  result <- switch(which,
    "hr_H" = dgm$hr_H_true,
    "hr_Hc" = dgm$hr_Hc_true,
    "ahr_H" = dgm$AHR_H_true,
    "ahr_Hc" = dgm$AHR_Hc_true,
    "hr_overall" = dgm$hr_causal,
    "ahr" = dgm$AHR,
    "cde_H" = dgm$cde_H,
    "cde_Hc" = dgm$cde_Hc,
    "cde" = dgm$CDE,
    NA_real_
  )

  if (is.null(result)) NA_real_ else result
}


#' Compute and Attach CDE Values to a DGM Object
#'
#' Calculates Controlled Direct Effect (CDE) hazard ratios from the
#' super-population potential outcomes (\code{theta_0}, \code{theta_1})
#' and attaches them to the DGM's \code{hazard_ratios} list. This enables
#' automatic CDE detection by \code{\link{build_estimation_table}}.
#'
#' @details
#' The CDE for subgroup S is defined as:
#'
#' \code{CDE(S) = mean(exp(theta_1[S])) / mean(exp(theta_0[S]))}
#'
#' which is the ratio of average hazard contributions on the natural scale.
#' This differs from the AHR (\code{exp(mean(loghr_po))}) due to Jensen's
#' inequality. In the notation of Leon et al. (2024), CDE corresponds to
#' theta-ddagger.
#'
#' The function detects the subgroup indicator column automatically,
#' checking for \code{flag.harm}, \code{flag_harm}, and \code{H} in
#' the super-population data frame.
#'
#' @param dgm A DGM object (e.g., from \code{create_gbsg_dgm} or
#'   \code{generate_aft_dgm_flex}). Must contain \code{df_super_rand}
#'   with columns \code{theta_0} and \code{theta_1}.
#' @param harm_col Character. Name of the subgroup indicator column in
#'   \code{dgm$df_super_rand}. If \code{NULL} (default), auto-detected
#'   from \code{flag.harm}, \code{flag_harm}, or \code{H}.
#'
#' @return The DGM object with CDE values added to
#'   \code{dgm$hazard_ratios} (\code{CDE}, \code{CDE_harm},
#'   \code{CDE_no_harm}) and to top-level fields (\code{dgm$CDE},
#'   \code{dgm$cde_H}, \code{dgm$cde_Hc}).
#'
#' @examples
#' \dontrun{
#' dgm <- create_gbsg_dgm(model = "alt", k_inter = 2.0)
#' dgm <- compute_dgm_cde(dgm)
#' dgm$hazard_ratios$CDE_harm   # theta-ddagger(H)
#' dgm$hazard_ratios$CDE         # theta-ddagger overall
#' }
#'
#' @seealso \code{\link{build_estimation_table}}, \code{\link{get_dgm_hr}}
#' @export
compute_dgm_cde <- function(dgm, harm_col = NULL) {

  # Support both gbsg_dgm (df_super_rand) and aft_dgm_flex (df_super)
  df <- dgm$df_super_rand %||% dgm$df_super
  if (is.null(df)) {
    warning("DGM has no df_super_rand or df_super; cannot compute CDE.")
    return(dgm)
  }
  if (!all(c("theta_0", "theta_1") %in% names(df))) {
    warning("df_super_rand lacks theta_0/theta_1; cannot compute CDE.")
    return(dgm)
  }

  # Auto-detect subgroup indicator
  if (is.null(harm_col)) {
    candidates <- c("flag.harm", "flag_harm", "H")
    found <- intersect(candidates, names(df))
    if (length(found) == 0) {
      warning("No subgroup indicator found (tried: flag.harm, flag_harm, H).")
      # Compute overall CDE only
      cde_overall <- mean(exp(df$theta_1)) / mean(exp(df$theta_0))
      dgm$CDE <- cde_overall
      if (is.null(dgm$hazard_ratios)) dgm$hazard_ratios <- list()
      dgm$hazard_ratios$CDE <- cde_overall
      return(dgm)
    }
    harm_col <- found[1]
  }

  in_H <- df[[harm_col]] == 1

  # CDE(S) = mean(exp(theta_1[S])) / mean(exp(theta_0[S]))
  cde_H  <- mean(exp(df$theta_1[in_H]))  / mean(exp(df$theta_0[in_H]))
  cde_Hc <- mean(exp(df$theta_1[!in_H])) / mean(exp(df$theta_0[!in_H]))
  cde_overall <- mean(exp(df$theta_1))    / mean(exp(df$theta_0))

  # Attach to hazard_ratios list (for get_dgm_hr auto-detection)
  if (is.null(dgm$hazard_ratios)) dgm$hazard_ratios <- list()
  dgm$hazard_ratios$CDE         <- cde_overall
  dgm$hazard_ratios$CDE_harm    <- cde_H
  dgm$hazard_ratios$CDE_no_harm <- cde_Hc

  # Also attach as top-level fields (for old-format fallback)
  dgm$CDE    <- cde_overall
  dgm$cde_H  <- cde_H
  dgm$cde_Hc <- cde_Hc

  dgm
}


#' Create Subgroup Indicator from Factor Definitions
#'
#' Parses factor definitions (e.g., "v1.1", "grade3.1") and creates
#' a binary indicator for subgroup membership.
#'
#' @param df Data frame containing the variables
#' @param sg_factors Character vector of factor definitions
#'
#' @return Integer vector (1 = in subgroup, 0 = not in subgroup)
#'
#' @keywords internal
create_subgroup_indicator <- function(df, sg_factors) {
  indicator <- rep(TRUE, nrow(df))

  for (factor_def in sg_factors) {
    if (is.na(factor_def) || factor_def == "") next

    parts <- strsplit(factor_def, "\\.")[[1]]
    if (length(parts) >= 2) {
      var_name <- parts[1]
      level <- parts[2]

      if (var_name %in% names(df)) {
        indicator <- indicator & (as.character(df[[var_name]]) == level)
      }
    }
  }

  as.integer(indicator)
}


#' Compute AHR from loghr_po
#'
#' Computes Average Hazard Ratio from individual log hazard ratios.
#'
#' @param df Data frame with loghr_po column
#' @param subset_indicator Optional logical/integer vector for subsetting
#'
#' @return Numeric AHR value
#'
#' @keywords internal
compute_ahr <- function(df, subset_indicator = NULL) {
  if (!"loghr_po" %in% names(df)) {
    return(NA_real_)
  }

  loghr <- df$loghr_po

  if (!is.null(subset_indicator)) {
    loghr <- loghr[subset_indicator == 1]
  }

  if (length(loghr) == 0 || all(is.na(loghr))) {
    return(NA_real_)
  }

  exp(mean(loghr, na.rm = TRUE))
}


#' Compute CDE from theta_0 and theta_1
#'
#' Computes Controlled Direct Effect as the ratio of average hazard
#' contributions on the natural scale:
#' \code{CDE(S) = mean(exp(theta_1[S])) / mean(exp(theta_0[S]))}.
#'
#' @param df Data frame with \code{theta_0} and \code{theta_1} columns.
#' @param subset_indicator Optional logical/integer vector for subsetting.
#'   If provided, only rows where \code{subset_indicator == 1} are used.
#'
#' @return Numeric CDE value, or \code{NA_real_} if columns are missing.
#'
#' @keywords internal
compute_cde <- function(df, subset_indicator = NULL) {
  if (!all(c("theta_0", "theta_1") %in% names(df))) {
    return(NA_real_)
  }

  t0 <- df$theta_0
  t1 <- df$theta_1

  if (!is.null(subset_indicator)) {
    idx <- subset_indicator == 1
    t0 <- t0[idx]
    t1 <- t1[idx]
  }

  if (length(t0) == 0 || all(is.na(t0))) {
    return(NA_real_)
  }

  mean(exp(t1), na.rm = TRUE) / mean(exp(t0), na.rm = TRUE)
}




# =============================================================================
# Summary Functions
# =============================================================================

#' Summarize Simulation Results
#'
#' Creates a summary table of operating characteristics across all simulations.
#' Includes both HR and AHR metrics.
#'
#' @param results data.table with simulation results from run_simulation_analysis
#' @param analyses Character vector. Analysis methods to include. Default: all
#' @param digits Integer. Decimal places for proportions. Default: 2
#' @param digits_hr Integer. Decimal places for hazard ratios. Default: 3
#'
#' @return Data frame with summary statistics
#'
#' @examples
#' \dontrun{
#' dgm <- create_gbsg_dgm()
#' sim_res <- run_simulation_analysis(dgm, nsim = 10)
#' summarize_simulation_results(sim_res)
#' }
#' @export
summarize_simulation_results <- function(
    results,
    analyses = NULL,
    digits = 2,
    digits_hr = 3
) {

  if (is.null(analyses)) {
    analyses <- unique(results$analysis)
  }

  summaries <- lapply(analyses, function(a) {
    res <- subset(results, analysis == a)
    summarize_single_analysis(res, digits = digits, digits_hr = digits_hr)
  })

  summary_df <- do.call(cbind, summaries)
  colnames(summary_df) <- analyses

  summary_df
}


#' Summarize Single Analysis Results
#'
#' @param result data.table with results for a single analysis method
#' @param digits Integer. Decimal places for proportions
#' @param digits_hr Integer. Decimal places for hazard ratios
#'
#' @return Data frame with summary statistics
#'
#' @keywords internal
summarize_single_analysis <- function(result, digits = 2, digits_hr = 3) {

  # ---------------------------------------------------------------------------
  # Classification metrics (all replicates)
  # ---------------------------------------------------------------------------
  class_cols <- c("any.H", "sens", "spec", "ppv", "npv")
  class_cols <- intersect(class_cols, names(result))
  class_means <- sapply(result[, class_cols, with = FALSE], mean, na.rm = TRUE)
  class_means <- round(class_means, digits)

  # ---------------------------------------------------------------------------
  # Subgroup sizes (H: conditional on detection; Hc: all replicates)
  # ---------------------------------------------------------------------------
  res_found <- subset(result, any.H == 1)

  if (nrow(res_found) > 0) {
    avg_H <- round(mean(res_found$size.H, na.rm = TRUE), 0)
    min_H <- round(min(res_found$size.H, na.rm = TRUE), 0)
    max_H <- round(max(res_found$size.H, na.rm = TRUE), 0)
  } else {
    avg_H <- min_H <- max_H <- NA
  }

  avg_Hc <- round(mean(result$size.Hc, na.rm = TRUE), 0)
  min_Hc <- round(min(result$size.Hc, na.rm = TRUE), 0)
  max_Hc <- round(max(result$size.Hc, na.rm = TRUE), 0)

  # ---------------------------------------------------------------------------
  # HR and AHR estimates (conditional on detection)
  # ---------------------------------------------------------------------------
  if (nrow(res_found) > 0) {
    # Cox-based HRs
    hr_H_true <- round(mean(res_found$hr.H.true, na.rm = TRUE), digits_hr)
    hr_H_hat <- round(mean(res_found$hr.H.hat, na.rm = TRUE), digits_hr)
    hr_Hc_true <- round(mean(res_found$hr.Hc.true, na.rm = TRUE), digits_hr)
    hr_Hc_hat <- round(mean(res_found$hr.Hc.hat, na.rm = TRUE), digits_hr)

    # AHR estimates (true subgroup and identified subgroup)
    ahr_H_true <- if ("ahr.H.true" %in% names(res_found)) {
      round(mean(res_found$ahr.H.true, na.rm = TRUE), digits_hr)
    } else NA
    ahr_Hc_true <- if ("ahr.Hc.true" %in% names(res_found)) {
      round(mean(res_found$ahr.Hc.true, na.rm = TRUE), digits_hr)
    } else NA
    ahr_H_hat <- if ("ahr.H.hat" %in% names(res_found)) {
      round(mean(res_found$ahr.H.hat, na.rm = TRUE), digits_hr)
    } else NA
    ahr_Hc_hat <- if ("ahr.Hc.hat" %in% names(res_found)) {
      round(mean(res_found$ahr.Hc.hat, na.rm = TRUE), digits_hr)
    } else NA
  } else {
    hr_H_true <- hr_H_hat <- hr_Hc_true <- hr_Hc_hat <- NA
    ahr_H_true <- ahr_Hc_true <- ahr_H_hat <- ahr_Hc_hat <- NA
  }

  # ---------------------------------------------------------------------------
  # Unconditional estimates (all replicates)
  # ---------------------------------------------------------------------------
  hr_H_true_all <- round(mean(result$hr.H.true, na.rm = TRUE), digits_hr)
  hr_Hc_true_all <- round(mean(result$hr.Hc.true, na.rm = TRUE), digits_hr)
  hr_itt <- round(mean(result$hr.itt, na.rm = TRUE), digits_hr)
  hr_adj_itt <- if ("hr.adj.itt" %in% names(result)) {
    round(mean(result$hr.adj.itt, na.rm = TRUE), digits_hr)
  } else NA

  # AHR unconditional (all replicates)
  ahr_H_true_all <- if ("ahr.H.true" %in% names(result)) {
    round(mean(result$ahr.H.true, na.rm = TRUE), digits_hr)
  } else NA
  ahr_Hc_true_all <- if ("ahr.Hc.true" %in% names(result)) {
    round(mean(result$ahr.Hc.true, na.rm = TRUE), digits_hr)
  } else NA

  # ---------------------------------------------------------------------------
  # Combine into output data frame
  # ---------------------------------------------------------------------------
  values <- c(
    class_means,
    avg_H, min_H, max_H, avg_Hc, min_Hc, max_Hc,
    hr_H_true, hr_H_hat, hr_Hc_true, hr_Hc_hat,
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
#' Creates a formatted gt table from simulation operating characteristics results.
#'
#' @param results data.table or data.frame. Simulation results from
#'   \code{\link{run_simulation_analysis}} or combined results from multiple simulations.
#' @param analyses Character vector. Analysis methods to include.
#'   Default: NULL (all analyses in results)
#' @param metrics Character vector. Metrics to display. Options include:
#'   "detection", "classification", "hr_estimates", "ahr_estimates",
#'   "cde_estimates", "subgroup_size", "all". Default: "all"
#' @param digits Integer. Decimal places for proportions. Default: 3
#' @param digits_hr Integer. Decimal places for hazard ratios. Default: 3
#' @param title Character. Table title. Default: "Operating Characteristics Summary"
#' @param subtitle Character. Table subtitle. Default: NULL
#' @param use_gt Logical. Return gt table if TRUE, data.frame if FALSE. Default: TRUE
#' @param subgroup_notation Character. \code{"harm"} (default) labels
#'   subgroups as H/Hc; \code{"benefit"} labels as Q/Qc for
#'   benefit-search analyses (treatment switching).
#'
#' @return A gt table object (if use_gt = TRUE and gt package available) or data.frame
#'
#' @details
#' The function summarizes simulation results across multiple metrics:
#' \itemize{
#'   \item \strong{Found}: Proportion of simulations finding a subgroup (any.H)
#'   \item \strong{Classification}: Sen, spec, PPV, NPV
#'   \item \strong{HR Estimates}: Mean Cox hazard ratios in true (H) and
#'     identified (H-hat) subgroups and their complements
#'   \item \strong{AHR Estimates}: Mean average hazard ratios (from loghr_po)
#'     in true and identified subgroups
#'   \item \strong{CDE Estimates}: Controlled direct effects (from
#'     theta_0/theta_1) in true and identified subgroups
#'   \item \strong{Subgroup Size}: Average, min, max sizes
#' }
#'
#' Column notation aligns with \code{\link{build_estimation_table}} and
#' Leon et al. (2024): \code{H} = true (oracle) subgroup, \code{H-hat} =
#' identified subgroup. The asterisk (\code{*}) is reserved for bootstrap
#' bias-corrected estimates and is not used in this summary table.
#'
#' @importFrom data.table is.data.table as.data.table
#' @examples
#' \dontrun{
#' # format_oc_results() is called by summarize_simulation_results().
#' # See run_simulation_analysis() for the standard entry point.
#' }
#' @export
format_oc_results <- function(
    results,
    analyses = NULL,
    metrics = "all",
    digits = 3,
    digits_hr = 3,
    title = "Operating Characteristics Summary",
    subtitle = NULL,
    use_gt = TRUE,
    subgroup_notation = c("harm", "benefit")
) {

  subgroup_notation <- match.arg(subgroup_notation)
  L <- get_sg_labels(subgroup_notation)

  # Benefit search: invert HRs from switched -> original scale
  if (subgroup_notation == "benefit") {
    results <- invert_hr_columns(results)
  }

  # Convert to data.table if needed
  if (!data.table::is.data.table(results)) {
    results <- data.table::as.data.table(results)
  }

  # Get analyses if not specified
  if (is.null(analyses)) {
    analyses <- unique(results$analysis)
  }

  # Filter to requested analyses
  results <- results[results$analysis %in% analyses, ]

  # Compute summary statistics for each analysis
  summary_list <- lapply(analyses, function(a) {
    res <- results[results$analysis == a, ]
    n_sims <- nrow(res)

    # Detection rate
    detection_rate <- mean(res$any.H, na.rm = TRUE)

    # Classification metrics (averaged across sims with detection).
    # Under null models (León et al. 2024, Sec 3): Q = emptyset,
    # so sens = NA and ppv = 0.  NaN from mean(all-NA) is replaced
    # with NA for clean table display.
    sens <- mean(res$sens, na.rm = TRUE)
    spec <- mean(res$spec, na.rm = TRUE)
    ppv <- mean(res$ppv, na.rm = TRUE)
    npv <- mean(res$npv, na.rm = TRUE)

    if (is.nan(sens)) sens <- NA
    if (is.nan(spec)) spec <- NA
    if (is.nan(ppv))  ppv  <- NA
    if (is.nan(npv))  npv  <- NA

    # HR estimates (only when subgroup found)
    res_found <- res[res$any.H == 1, ]
    if (nrow(res_found) > 0) {
      hr_H_hat <- mean(res_found$hr.H.hat, na.rm = TRUE)
      hr_Hc_hat <- mean(res_found$hr.Hc.hat, na.rm = TRUE)
      hr_H_true <- mean(res_found$hr.H.true, na.rm = TRUE)
      hr_Hc_true <- mean(res_found$hr.Hc.true, na.rm = TRUE)
      size_H_mean <- mean(res_found$size.H, na.rm = TRUE)
      size_H_min <- min(res_found$size.H, na.rm = TRUE)
      size_H_max <- max(res_found$size.H, na.rm = TRUE)
    } else {
      hr_H_hat <- hr_Hc_hat <- hr_H_true <- hr_Hc_true <- NA
      size_H_mean <- size_H_min <- size_H_max <- NA
    }

    # AHR estimates (conditional on detection)
    if (nrow(res_found) > 0) {
      ahr_H_true <- if ("ahr.H.true" %in% names(res_found)) {
        mean(res_found$ahr.H.true, na.rm = TRUE)
      } else NA
      ahr_Hc_true <- if ("ahr.Hc.true" %in% names(res_found)) {
        mean(res_found$ahr.Hc.true, na.rm = TRUE)
      } else NA
      ahr_H_hat <- if ("ahr.H.hat" %in% names(res_found)) {
        mean(res_found$ahr.H.hat, na.rm = TRUE)
      } else NA
      ahr_Hc_hat <- if ("ahr.Hc.hat" %in% names(res_found)) {
        mean(res_found$ahr.Hc.hat, na.rm = TRUE)
      } else NA
    } else {
      ahr_H_true <- ahr_Hc_true <- ahr_H_hat <- ahr_Hc_hat <- NA
    }

    # CDE estimates (conditional on detection)
    if (nrow(res_found) > 0) {
      cde_H_true <- if ("cde.H.true" %in% names(res_found)) {
        mean(res_found$cde.H.true, na.rm = TRUE)
      } else NA
      cde_Hc_true <- if ("cde.Hc.true" %in% names(res_found)) {
        mean(res_found$cde.Hc.true, na.rm = TRUE)
      } else NA
      cde_H_hat <- if ("cde.H.hat" %in% names(res_found)) {
        mean(res_found$cde.H.hat, na.rm = TRUE)
      } else NA
      cde_Hc_hat <- if ("cde.Hc.hat" %in% names(res_found)) {
        mean(res_found$cde.Hc.hat, na.rm = TRUE)
      } else NA
    } else {
      cde_H_true <- cde_Hc_true <- cde_H_hat <- cde_Hc_hat <- NA
    }

    # ITT estimate (all simulations)
    hr_itt <- mean(res$hr.itt, na.rm = TRUE)

    data.frame(
      Analysis = a,
      N_sims = n_sims,
      Detection = detection_rate,
      Sen = sens,
      Spec = spec,
      PPV = ppv,
      NPV = npv,
      HR_H_hat = hr_H_hat,
      HR_Hc_hat = hr_Hc_hat,
      HR_H_true = hr_H_true,
      HR_Hc_true = hr_Hc_true,
      HR_ITT = hr_itt,
      AHR_H_true = ahr_H_true,
      AHR_Hc_true = ahr_Hc_true,
      AHR_H_hat = ahr_H_hat,
      AHR_Hc_hat = ahr_Hc_hat,
      CDE_H_true = cde_H_true,
      CDE_Hc_true = cde_Hc_true,
      CDE_H_hat = cde_H_hat,
      CDE_Hc_hat = cde_Hc_hat,
      Size_H_mean = size_H_mean,
      Size_H_min = size_H_min,
      Size_H_max = size_H_max,
      stringsAsFactors = FALSE
    )
  })

  summary_df <- do.call(rbind, summary_list)

  # Filter columns based on metrics
  if (!"all" %in% metrics) {
    cols_to_keep <- c("Analysis", "N_sims")
    if ("detection" %in% metrics) {
      cols_to_keep <- c(cols_to_keep, "Detection")
    }
    if ("classification" %in% metrics) {
      cols_to_keep <- c(cols_to_keep, "Sen", "Spec", "PPV", "NPV")
    }
    if ("hr_estimates" %in% metrics) {
      cols_to_keep <- c(cols_to_keep, "HR_H_hat", "HR_Hc_hat",
                        "HR_H_true", "HR_Hc_true", "HR_ITT")
    }
    if ("ahr_estimates" %in% metrics) {
      cols_to_keep <- c(cols_to_keep, "AHR_H_true", "AHR_Hc_true",
                        "AHR_H_hat", "AHR_Hc_hat")
    }
    if ("cde_estimates" %in% metrics) {
      cols_to_keep <- c(cols_to_keep, "CDE_H_true", "CDE_Hc_true",
                        "CDE_H_hat", "CDE_Hc_hat")
    }
    if ("subgroup_size" %in% metrics) {
      cols_to_keep <- c(cols_to_keep, "Size_H_mean", "Size_H_min", "Size_H_max")
    }
    # Only keep columns that exist (AHR may be all-NA and still present)
    cols_to_keep <- intersect(cols_to_keep, names(summary_df))
    summary_df <- summary_df[, cols_to_keep, drop = FALSE]
  }

  # Format as gt table if requested and available
  if (use_gt && requireNamespace("gt", quietly = TRUE)) {
    gt_table <- gt::gt(summary_df)

    # Add title
    gt_table <- gt::tab_header(
      gt_table,
      title = title,
      subtitle = subtitle
    )

    # Proportion columns (0-1 scale)
    prop_cols <- intersect(c("Detection", "Sen", "Spec", "PPV", "NPV"),
                           names(summary_df))
    if (length(prop_cols) > 0) {
      gt_table <- gt::fmt_number(
        gt_table,
        columns = gt::all_of(prop_cols),
        decimals = digits
      )
    }

    # Cox HR columns
    hr_cols <- intersect(c("HR_H_hat", "HR_Hc_hat", "HR_H_true",
                           "HR_Hc_true", "HR_ITT"),
                         names(summary_df))
    if (length(hr_cols) > 0) {
      gt_table <- gt::fmt_number(
        gt_table,
        columns = gt::all_of(hr_cols),
        decimals = digits_hr
      )
    }

    # AHR columns
    ahr_cols <- intersect(c("AHR_H_true", "AHR_Hc_true",
                            "AHR_H_hat", "AHR_Hc_hat"),
                          names(summary_df))
    if (length(ahr_cols) > 0) {
      gt_table <- gt::fmt_number(
        gt_table,
        columns = gt::all_of(ahr_cols),
        decimals = digits_hr
      )
    }

    # CDE columns
    cde_cols <- intersect(c("CDE_H_true", "CDE_Hc_true",
                            "CDE_H_hat", "CDE_Hc_hat"),
                          names(summary_df))
    if (length(cde_cols) > 0) {
      gt_table <- gt::fmt_number(
        gt_table,
        columns = gt::all_of(cde_cols),
        decimals = digits_hr
      )
    }

    # Size columns
    size_cols <- intersect(c("Size_H_mean", "Size_H_min", "Size_H_max"),
                           names(summary_df))
    if (length(size_cols) > 0) {
      gt_table <- gt::fmt_number(
        gt_table,
        columns = gt::all_of(size_cols),
        decimals = 0
      )
    }

    # Rename columns for display using math notation (pure Unicode)
    # Notation-aware: H/Hc for harm search, Q/Qc for benefit search
    # Unicode: \u0302 = combining circumflex, \u1D9C = superscript-c,
    #          \u2020 = dagger, \u2021 = double-dagger,
    #          \u03b8 = theta, \u00e2 = a-hat

    label_list <- list(
      Analysis  = "Method",
      N_sims    = "Sims",
      Detection = "Found"
    )

    # HR column labels (conditional on existence)
    hr_label_map <- list(
      HR_H_hat   = sprintf("\u03b8\u0302(%s)", L$sg_hat),
      HR_Hc_hat  = sprintf("\u03b8\u0302(%s)", L$sg_hat_c),
      HR_H_true  = sprintf("\u03b8\u0302(%s)", L$sg),
      HR_Hc_true = sprintf("\u03b8\u0302(%s)", L$sg_c),
      HR_ITT     = "\u03b8\u0302(ITT)"
    )
    for (col_nm in names(hr_label_map)) {
      if (col_nm %in% names(summary_df)) {
        label_list[[col_nm]] <- hr_label_map[[col_nm]]
      }
    }

    # AHR column labels (conditional on existence)
    ahr_label_map <- list(
      AHR_H_true  = sprintf("\u00e2hr(%s)", L$sg),
      AHR_Hc_true = sprintf("\u00e2hr(%s)", L$sg_c),
      AHR_H_hat   = sprintf("\u00e2hr(%s)", L$sg_hat),
      AHR_Hc_hat  = sprintf("\u00e2hr(%s)", L$sg_hat_c)
    )
    for (col_nm in names(ahr_label_map)) {
      if (col_nm %in% names(summary_df)) {
        label_list[[col_nm]] <- ahr_label_map[[col_nm]]
      }
    }

    # CDE column labels (conditional on existence)
    cde_label_map <- list(
      CDE_H_true  = sprintf("\u03b8\u2021(%s)", L$sg),
      CDE_Hc_true = sprintf("\u03b8\u2021(%s)", L$sg_c),
      CDE_H_hat   = sprintf("\u03b8\u2021(%s)", L$sg_hat),
      CDE_Hc_hat  = sprintf("\u03b8\u2021(%s)", L$sg_hat_c)
    )
    for (col_nm in names(cde_label_map)) {
      if (col_nm %in% names(summary_df)) {
        label_list[[col_nm]] <- cde_label_map[[col_nm]]
      }
    }

    gt_table <- gt::cols_label(gt_table, .list = label_list)

    # Add column spanners
    if ("all" %in% metrics || "classification" %in% metrics) {
      class_cols <- intersect(c("Sen", "Spec", "PPV", "NPV"), names(summary_df))
      if (length(class_cols) > 0) {
        gt_table <- gt::tab_spanner(
          gt_table,
          label = "Classification",
          columns = gt::all_of(class_cols)
        )
      }
    }

    if ("all" %in% metrics || "hr_estimates" %in% metrics) {
      hr_span_cols <- intersect(c("HR_H_hat", "HR_Hc_hat", "HR_H_true",
                                  "HR_Hc_true", "HR_ITT"),
                                names(summary_df))
      if (length(hr_span_cols) > 0) {
        gt_table <- gt::tab_spanner(
          gt_table,
          label = "Cox Hazard Ratios",
          columns = gt::all_of(hr_span_cols)
        )
      }
    }

    if ("all" %in% metrics || "ahr_estimates" %in% metrics) {
      ahr_span_cols <- intersect(c("AHR_H_true", "AHR_Hc_true",
                                   "AHR_H_hat", "AHR_Hc_hat"),
                                 names(summary_df))
      if (length(ahr_span_cols) > 0) {
        gt_table <- gt::tab_spanner(
          gt_table,
          label = "Average Hazard Ratios",
          columns = gt::all_of(ahr_span_cols)
        )
      }
    }

    if ("all" %in% metrics || "cde_estimates" %in% metrics) {
      cde_span_cols <- intersect(c("CDE_H_true", "CDE_Hc_true",
                                   "CDE_H_hat", "CDE_Hc_hat"),
                                 names(summary_df))
      if (length(cde_span_cols) > 0) {
        gt_table <- gt::tab_spanner(
          gt_table,
          label = "Controlled Direct Effects",
          columns = gt::all_of(cde_span_cols)
        )
      }
    }

    # Style
    gt_table <- gt::tab_style(
      gt_table,
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_column_labels()
    )

    # Handle missing values
    if (utils::packageVersion("gt") >= "0.6.0") {
      gt_table <- gt::sub_missing(
        gt_table,
        columns = gt::everything(),
        missing_text = "-"
      )
    }

    return(gt_table)
  }

  # Return data.frame if gt not available or not requested
  summary_df
}
