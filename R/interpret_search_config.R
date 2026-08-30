# =============================================================================
# interpret_search_config() -- Search alignment diagnostic
#
# Prints a human-readable interpretation of the user's ForestSearch
# parameter combination, explaining:
#   - What the search will find (effect direction)
#   - How GRF pre-screening is aligned
#   - How LASSO screening is aligned
#   - Clinical interpretation of the threshold
#   - Potential misalignment warnings
#
# Called from forestsearch() after parameter resolution.
# =============================================================================


#' Interpret and Diagnose Search Configuration
#'
#' Generates a human-readable diagnostic of how ForestSearch will
#' interpret the user's parameter combination.  Prints via
#' \code{message()} and returns an invisible list of diagnostics.
#'
#' @param outcome_type Character. One of \code{"survival"}, \code{"binary"},
#'   \code{"continuous"}, \code{"count"}.
#' @param effect_measure Character. \code{"HR"}, \code{"OR"}, \code{"RR"},
#'   \code{"RD"}, \code{"MD"}, \code{"IRR"}, \code{"IRD"}.
#' @param adverse_outcome Logical. As resolved by forestsearch().
#' @param effect_threshold Numeric. Resolved screening threshold (on
#'   log scale for ratio measures, identity for others).
#' @param consistency_threshold Numeric. Resolved consistency threshold.
#' @param use_lasso Logical.
#' @param use_grf Logical.  Whether GRF candidate-cut generation is on;
#'   see \code{\link{forestsearch}}.
#' @param outcome.name Character. Name of outcome column.
#' @param event.name Character. Name of event column.
#' @param treat.name Character. Name of treatment column.
#' @param offset.name Character or NULL.
#' @param quiet Logical. If TRUE, suppress output.
#'
#' @return Invisible list with diagnostic fields.
#' @export
interpret_search_config <- function(
    outcome_type,
    effect_measure,
    adverse_outcome,
    effect_threshold,
    consistency_threshold,
    use_lasso,
    use_grf,
    outcome.name,
    event.name,
    treat.name,
    offset.name   = NULL,
    quiet         = FALSE
) {

  is_survival  <- outcome_type == "survival"
  is_binary    <- outcome_type == "binary"
  is_cont      <- outcome_type == "continuous"
  is_count     <- outcome_type == "count"
  is_ratio     <- effect_measure %in% c("HR", "OR", "RR", "IRR")
  is_identity  <- effect_measure %in% c("RD", "IRD", "MD")

  # -- Core search direction ----------------------------------------------
  # ForestSearch always uses: estimate >= threshold
  # Determine what this means clinically.

  if (is_survival) {
    measure_desc  <- "Hazard Ratio (HR)"
    direction_pos <- "treatment INCREASES hazard (more/faster events)"
    direction_neg <- "treatment DECREASES hazard (fewer/slower events)"
    search_finds  <- paste0(
      "Subgroups where ", treat.name, " increases the event hazard\n",
      "    (HR >= threshold in the candidate subgroup)")
    grf_note <- paste0(
      "causal_survival_forest() estimates treatment effect on\n",
      "    survival probability. Positive CATE = treatment increases\n",
      "    survival = treatment HELPS (events assumed adverse).")

  } else if (is_binary) {
    if (effect_measure == "OR") {
      ao_label <- if (adverse_outcome) "as-is (no flip)" else "flipped (1-Y)"
      measure_desc <- paste0("Odds Ratio (OR) -- estimator sees Y ", ao_label)

      if (adverse_outcome) {
        # No flip: OR computed on Y directly
        search_finds <- paste0(
          "Subgroups where ", treat.name, " increases P(",
          outcome.name, "=1)\n",
          "    (OR >= threshold on the ORIGINAL outcome scale)")
        direction_pos <- paste0(
          treat.name, " increases P(", outcome.name, "=1)")
      } else {
        # Flip: OR computed on 1-Y
        search_finds <- paste0(
          "Subgroups where ", treat.name, " increases P(",
          outcome.name, "=0)\n",
          "    (OR >= threshold on the FAILURE/complement scale)")
        direction_pos <- paste0(
          treat.name, " increases P(1-", outcome.name, "=1) = P(",
          outcome.name, "=0)")
      }

    } else if (effect_measure == "RD") {
      ao_label <- if (adverse_outcome) "as-is" else "flipped (1-Y)"
      measure_desc <- paste0("Risk Difference (RD) -- estimator sees Y ", ao_label)
      search_finds <- paste0(
        "Subgroups where ", treat.name, " increases the event/outcome rate\n",
        "    (RD >= threshold)")
      direction_pos <- paste0(treat.name, " increases event rate")

    } else {
      measure_desc <- paste0(effect_measure, " (binary)")
      search_finds <- paste0("Subgroups where ", effect_measure, " >= threshold")
      direction_pos <- paste0(treat.name, " increases outcome")
    }

    # GRF note for binary
    if (adverse_outcome) {
      grf_note <- paste0(
        "GRF sees 1-", outcome.name, " (flipped for causal_forest).\n",
        "    Positive CATE = treatment REDUCES ", outcome.name,
        " = treatment HELPS.")
    } else {
      grf_note <- paste0(
        "GRF sees ", outcome.name, " directly (no flip).\n",
        "    Positive CATE = treatment INCREASES ", outcome.name,
        " = treatment HELPS.")
    }

  } else if (is_cont) {
    measure_desc <- "Mean Difference (MD)"
    search_finds <- paste0(
      "Subgroups where ", treat.name, " increases mean(",
      outcome.name, ")\n",
      "    (MD >= threshold)")
    direction_pos <- paste0(treat.name, " increases mean(", outcome.name, ")")
    direction_neg <- paste0(treat.name, " decreases mean(", outcome.name, ")")

    if (adverse_outcome) {
      grf_note <- paste0(
        "GRF sees 1-", outcome.name,
        " (flipped: higher ", outcome.name, " = worse).\n",
        "    Positive CATE = treatment REDUCES ", outcome.name, ".")
    } else {
      grf_note <- paste0(
        "GRF sees ", outcome.name,
        " directly (higher = better).\n",
        "    Positive CATE = treatment INCREASES ", outcome.name, ".")
    }

  } else if (is_count) {
    measure_desc <- if (effect_measure == "IRR") {
      "Incidence Rate Ratio (IRR)"
    } else {
      "Incidence Rate Difference (IRD)"
    }
    off_desc <- if (!is.null(offset.name)) {
      paste0(" with offset = log(", offset.name, ")")
    } else ""
    search_finds <- paste0(
      "Subgroups where ", treat.name, " increases the event rate",
      off_desc, "\n",
      "    (", effect_measure, " >= threshold)")
    direction_pos <- paste0(treat.name, " increases event rate")

    grf_note <- paste0(
      "GRF sees log(", outcome.name,
      " + 0.5) (variance-stabilising transform).\n",
      "    Variable importance identifies covariates with\n",
      "    heterogeneous treatment effects on the count scale.")
  }

  # -- Threshold interpretation --------------------------------------------
  if (is_ratio) {
    thresh_nat <- exp(effect_threshold)
    thresh_desc <- sprintf(
      "Screening: %s >= %.3f in the candidate subgroup\n    Consistency: %s >= %.3f in each bootstrap split",
      effect_measure, thresh_nat, effect_measure, exp(consistency_threshold))
  } else {
    thresh_desc <- sprintf(
      "Screening: %s >= %.4f in the candidate subgroup\n    Consistency: %s >= %.4f in each bootstrap split",
      effect_measure, effect_threshold, effect_measure, consistency_threshold)
  }

  # -- LASSO note ----------------------------------------------------------
  if (use_lasso) {
    lasso_family <- switch(outcome_type,
      survival   = "Cox (coxnet)",
      binary     = "binomial (logistic)",
      continuous = "gaussian",
      count      = "poisson"
    )
    lasso_note <- paste0(
      "LASSO screening uses glmnet with family = \"", lasso_family,
      "\".\n    Variables with non-zero coefficients are retained.")
  } else {
    lasso_note <- "LASSO screening: DISABLED"
  }

  # -- Alignment warnings -------------------------------------------------
  warnings_list <- character(0)

  # Warn: binary + adverse_outcome=FALSE + benefit-like threshold
  if (is_binary && !adverse_outcome) {
    warnings_list <- c(warnings_list, paste0(
      "adverse_outcome = FALSE: the binary estimator flips Y -> 1-Y\n",
      "      internally. This means OR/RD is computed on the COMPLEMENT\n",
      "      of ", outcome.name, ". If you are doing a benefit search with\n",
      "      switched treatment, you likely want adverse_outcome = TRUE\n",
      "      (the default for binary)."))
  }

  # Warn: continuous + adverse_outcome=TRUE but GRF default is FALSE
  if (is_cont && adverse_outcome) {
    warnings_list <- c(warnings_list, paste0(
      "adverse_outcome = TRUE for continuous outcome: GRF will see\n",
      "      1-", outcome.name, ". Confirm that higher ", outcome.name,
      " = worse\n",
      "      (e.g., pain score, CD4 loss). If higher = better, use\n",
      "      adverse_outcome = FALSE (the default for continuous)."))
  }

  # Warn: survival with GRF but no LASSO
  if (is_survival && use_grf && !use_lasso) {
    warnings_list <- c(warnings_list,
      "Survival with GRF only (no LASSO): causal_survival_forest()\n      assumes events are adverse. If the event is positive (recovery),\n      variable importance signs may be reversed (magnitude still valid).")
  }

  # -- Print --------------------------------------------------------------
  if (!quiet) {
    message("\n", strrep("=", 62))
    message("  Search Alignment Diagnostic")
    message(strrep("=", 62))

    message(sprintf("\n  Outcome:    %s (%s)", outcome.name, outcome_type))
    if (!is.null(offset.name))
      message(sprintf("  Offset:     %s", offset.name))
    message(sprintf("  Treatment:  %s", treat.name))
    message(sprintf("  Measure:    %s", measure_desc))

    message(sprintf("\n  ForestSearch searches for:\n    %s", search_finds))
    message(sprintf("\n  Thresholds:\n    %s", thresh_desc))

    message(sprintf("\n  Variable screening:"))
    message(sprintf("    LASSO: %s", lasso_note))
    if (use_grf) {
      message(sprintf("    GRF:   %s", grf_note))
    } else {
      message("    GRF:   DISABLED")
    }

    if (length(warnings_list) > 0) {
      message(sprintf("\n  %s ALIGNMENT NOTES:", length(warnings_list)))
      for (i in seq_along(warnings_list)) {
        message(sprintf("    [%d] %s", i, warnings_list[[i]]))
      }
    } else {
      message("\n  No alignment concerns detected.")
    }

    message(strrep("=", 62), "\n")
  }

  # -- Return -------------------------------------------------------------
  invisible(list(
    outcome_type         = outcome_type,
    effect_measure       = effect_measure,
    adverse_outcome      = adverse_outcome,
    search_finds         = search_finds,
    threshold_desc       = thresh_desc,
    lasso_note           = lasso_note,
    grf_note             = if (use_grf) grf_note else "GRF disabled",
    warnings             = warnings_list,
    n_warnings           = length(warnings_list)
  ))
}
