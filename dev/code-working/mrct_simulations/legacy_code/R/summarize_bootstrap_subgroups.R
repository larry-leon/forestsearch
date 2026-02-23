#' Summarize Bootstrap Subgroup Analysis Results
#'
#' Comprehensive summary of bootstrap subgroup identification results including
#' basic statistics, factor frequencies, consistency distributions, and agreement
#' with the original analysis subgroup.
#'
#' @param results Data.table or data.frame. Bootstrap results with subgroup
#'   characteristics including columns like Pcons, hr_sg, N_sg, K_sg, and M.1-M.k
#' @param nb_boots Integer. Total number of bootstrap iterations
#' @param original_sg Character vector. Original subgroup definition from main

#'   analysis (e.g., c("\{age>=50\}", "\{nodes>=3\}") for a 2-factor subgroup)
#' @param maxk Integer. Maximum number of factors allowed in subgroup definition
#'
#' @return List with summary components:
#'   \describe{
#'     \item{basic_stats}{Data.table of summary statistics}
#'     \item{consistency_dist}{Data.table of Pcons distribution by bins}
#'     \item{size_dist}{Data.table of subgroup size distribution}
#'     \item{factor_freq}{Data.table of factor frequencies by position}
#'     \item{agreement}{Data.table of subgroup definition agreement counts}
#'     \item{factor_presence}{Data.table of base factor presence counts}
#'     \item{factor_presence_specific}{Data.table of specific factor definitions}
#'     \item{original_agreement}{Data.table comparing to original analysis subgroup}
#'     \item{n_found}{Integer. Number of successful iterations}
#'     \item{pct_found}{Numeric. Percentage of successful iterations}
#'   }
#'
#' @importFrom data.table data.table .N setnames setcolorder copy as.data.table
#'   rbindlist
#' @importFrom stats median sd quantile
#' @keywords internal
summarize_bootstrap_subgroups <- function(results,
                                          nb_boots,
                                          original_sg = NULL,
                                          maxk = 2) {


  # ---------------------------------------------------------------------------

# Declare data.table variables to avoid R CMD check NOTEs
  # ---------------------------------------------------------------------------
  Pcons <- K_sg <- Subgroup <- N <- Percent_of_successful <- NULL
  M.1 <- M.2 <- hr_sg <- N_sg <- Pcons_bin <- N_sg_bin <- NULL

  # ===========================================================================
  # SECTION 1: INPUT VALIDATION AND DATA PREPARATION
  # ===========================================================================

  # Convert to data.table if needed
  if (!inherits(results, "data.table")) {
    if (is.data.frame(results)) {
      results <- data.table::as.data.table(results)
    } else if (is.matrix(results)) {
      results <- data.table::as.data.table(as.data.frame(results))
    } else {
      stop("results must be a data.frame, data.table, or matrix")
    }
  }

  # Make a copy to avoid modifying the original
  results <- data.table::copy(results)

  # Filter to successful iterations (where subgroup was found)
  if ("Pcons" %in% names(results)) {
    sg_found <- results[!is.na(results$Pcons), ]
  } else {
    warning("Column 'Pcons' not found in results")
    sg_found <- results
  }

  n_found <- nrow(sg_found)
  pct_found <- 100 * n_found / nb_boots

  # Early return if no subgroups found
  if (n_found == 0) {
    warning("No subgroups identified in any bootstrap iteration")
    return(list(
      basic_stats = NULL,
      consistency_dist = NULL,
      size_dist = NULL,
      factor_freq = NULL,
      agreement = NULL,
      factor_presence = NULL,
      factor_presence_specific = NULL,
      original_agreement = NULL,
      n_found = 0,
      pct_found = 0
    ))
  }

  # ===========================================================================
  # SECTION 2: BASIC STATISTICS TABLE
  # ===========================================================================

  metric_vec <- c(
    "Total bootstrap iterations",
    "Subgroups identified",
    "Success rate (%)"
  )

  value_vec <- c(
    as.character(nb_boots),
    as.character(n_found),
    sprintf("%.1f%%", pct_found)
  )

  # Add Pcons statistics if available
  if ("Pcons" %in% names(sg_found) && sum(!is.na(sg_found$Pcons)) > 0) {
    pcons_vals <- sg_found$Pcons[!is.na(sg_found$Pcons)]

    pcons_metrics <- c(
      "",
      "Consistency (Pcons)",
      "  Mean",
      "  Median",
      "  SD",
      "  Min",
      "  Max",
      "  Q25",
      "  Q75"
    )

    pcons_values <- c(
      "",
      "",
      sprintf("%.3f", mean(pcons_vals)),
      sprintf("%.3f", median(pcons_vals)),
      sprintf("%.3f", sd(pcons_vals)),
      sprintf("%.3f", min(pcons_vals)),
      sprintf("%.3f", max(pcons_vals)),
      sprintf("%.3f", quantile(pcons_vals, 0.25)),
      sprintf("%.3f", quantile(pcons_vals, 0.75))
    )

    metric_vec <- c(metric_vec, pcons_metrics)
    value_vec <- c(value_vec, pcons_values)
  }

  # Add hr_sg statistics if available
  if ("hr_sg" %in% names(sg_found) && sum(!is.na(sg_found$hr_sg)) > 0) {
    hr_vals <- sg_found$hr_sg[!is.na(sg_found$hr_sg)]

    hr_metrics <- c(
      "",
      "Hazard Ratio (hr_sg)",
      "  Mean",
      "  Median",
      "  SD",
      "  Min",
      "  Max",
      "  Q25",
      "  Q75"
    )

    hr_values <- c(
      "",
      "",
      sprintf("%.2f", mean(hr_vals)),
      sprintf("%.2f", median(hr_vals)),
      sprintf("%.2f", sd(hr_vals)),
      sprintf("%.2f", min(hr_vals)),
      sprintf("%.2f", max(hr_vals)),
      sprintf("%.2f", quantile(hr_vals, 0.25)),
      sprintf("%.2f", quantile(hr_vals, 0.75))
    )

    metric_vec <- c(metric_vec, hr_metrics)
    value_vec <- c(value_vec, hr_values)
  }

  # Add N_sg statistics if available
  if ("N_sg" %in% names(sg_found) && sum(!is.na(sg_found$N_sg)) > 0) {
    n_vals <- sg_found$N_sg[!is.na(sg_found$N_sg)]

    n_metrics <- c(
      "",
      "Subgroup Size (N_sg)",
      "  Mean",
      "  Median",
      "  SD",
      "  Min",
      "  Max"
    )

    n_values <- c(
      "",
      "",
      sprintf("%.1f", mean(n_vals)),
      sprintf("%.0f", median(n_vals)),
      sprintf("%.1f", sd(n_vals)),
      sprintf("%d", min(n_vals)),
      sprintf("%d", max(n_vals))
    )

    metric_vec <- c(metric_vec, n_metrics)
    value_vec <- c(value_vec, n_values)
  }

  # Add K_sg statistics if available
  if ("K_sg" %in% names(sg_found) && sum(!is.na(sg_found$K_sg)) > 0) {
    k_vals <- sg_found$K_sg[!is.na(sg_found$K_sg)]

    k_metrics <- c(
      "",
      "Number of Factors (K_sg)",
      "  Mean",
      "  Mode",
      "  Distribution"
    )

    # Calculate mode
    k_table <- table(k_vals)
    k_mode <- as.integer(names(k_table)[which.max(k_table)])

    # Distribution string
    k_dist <- paste(
      vapply(sort(unique(k_vals)), function(k) {
        sprintf("K=%d: %d (%.1f%%)", k, sum(k_vals == k), 100 * mean(k_vals == k))
      }, character(1)),
      collapse = "; "
    )

    k_values <- c(
      "",
      "",
      sprintf("%.2f", mean(k_vals)),
      as.character(k_mode),
      k_dist
    )

    metric_vec <- c(metric_vec, k_metrics)
    value_vec <- c(value_vec, k_values)
  }

  # Create basic_stats data.table
  basic_stats <- data.table::data.table(
    Metric = metric_vec,
    Value = value_vec
  )

  # ===========================================================================
  # SECTION 3: FACTOR FREQUENCY TABLE
  # ===========================================================================

  factor_cols <- paste0("M.", seq_len(maxk))
  factor_freq_list <- list()

  for (i in seq_len(maxk)) {
    col <- paste0("M.", i)
    if (col %in% names(sg_found)) {
      valid_rows <- !is.na(sg_found[[col]]) & sg_found[[col]] != ""
      if (sum(valid_rows) > 0) {
        factor_vals <- sg_found[[col]][valid_rows]
        freq_table <- table(factor_vals)
        n_position <- sum(valid_rows)

        freq <- data.table::data.table(
          Factor = names(freq_table),
          N = as.integer(freq_table),
          Position = paste0("M.", i),
          Percent = 100 * as.integer(freq_table) / n_position
        )

        freq <- freq[order(-N)]
        factor_freq_list[[i]] <- freq
      }
    }
  }

  if (length(factor_freq_list) > 0) {
    factor_freq <- data.table::rbindlist(factor_freq_list, fill = TRUE)
    data.table::setcolorder(factor_freq, c("Position", "Factor", "N", "Percent"))
  } else {
    factor_freq <- data.table::data.table(
      Position = character(),
      Factor = character(),
      N = integer(),
      Percent = numeric()
    )
  }

  # ===========================================================================
  # SECTION 4: SUBGROUP DEFINITION AGREEMENT
  # ===========================================================================

  agreement <- NULL

  if (maxk == 1) {
    if ("M.1" %in% names(sg_found)) {
      valid_rows <- !is.na(sg_found$M.1) & sg_found$M.1 != ""
      if (sum(valid_rows) > 0) {
        freq_table <- table(sg_found$M.1[valid_rows])
        agreement <- data.table::data.table(
          Subgroup = names(freq_table),
          K_sg = 1L,
          N = as.integer(freq_table)
        )
      }
    }
  } else if (maxk >= 2) {
    if (all(c("M.1", "M.2") %in% names(sg_found))) {
      valid_rows <- !is.na(sg_found$M.1) & sg_found$M.1 != ""
      if (sum(valid_rows) > 0) {
        sg_temp <- sg_found[valid_rows, ]
        sg_temp$Subgroup <- ifelse(
          is.na(sg_temp$M.2) | sg_temp$M.2 == "",
          sg_temp$M.1,
          paste(sg_temp$M.1, "&", sg_temp$M.2)
        )

        if ("K_sg" %in% names(sg_temp)) {
          agreement <- sg_temp[, .N, by = .(Subgroup, K_sg)]
        } else {
          agreement <- sg_temp[, .N, by = Subgroup]
          agreement$K_sg <- NA_integer_
        }
      }
    }
  }

  # Add percentage and sort if agreement exists
  if (!is.null(agreement) && nrow(agreement) > 0) {
    agreement$Percent_of_successful <- 100 * agreement$N / n_found
    agreement <- agreement[order(-agreement$N)]
    agreement$Rank <- seq_len(nrow(agreement))
  }

  # ===========================================================================
  # SECTION 5: CONSISTENCY DISTRIBUTION
  # ===========================================================================

  consistency_dist <- NULL

  if ("Pcons" %in% names(sg_found) && sum(!is.na(sg_found$Pcons)) > 0) {
    sg_found$Pcons_bin <- cut(
      sg_found$Pcons,
      breaks = c(0, 0.5, 0.7, 0.8, 0.85, 0.95, 1.0),
      labels = c("<0.5", "0.5-0.7", "0.7-0.8", "0.8-0.85", "0.85-0.95", ">=0.95"),
      include.lowest = TRUE
    )

    if (sum(!is.na(sg_found$Pcons_bin)) > 0) {
      freq_table <- table(sg_found$Pcons_bin[!is.na(sg_found$Pcons_bin)])
      consistency_dist <- data.table::data.table(
        `Consistency Range` = names(freq_table),
        Count = as.integer(freq_table),
        Percent = 100 * as.integer(freq_table) / n_found
      )
    }
  }

  # ===========================================================================
  # SECTION 6: SIZE DISTRIBUTION
  # ===========================================================================

  size_dist <- NULL

  if ("N_sg" %in% names(sg_found) && sum(!is.na(sg_found$N_sg)) > 0) {
    size_breaks <- c(0, 50, 100, 150, 200, 300, Inf)
    size_labels <- c("<50", "50-99", "100-149", "150-199", "200-299", ">=300")

    sg_found$N_sg_bin <- cut(
      sg_found$N_sg,
      breaks = size_breaks,
      labels = size_labels,
      include.lowest = TRUE
    )

    if (sum(!is.na(sg_found$N_sg_bin)) > 0) {
      freq_table <- table(sg_found$N_sg_bin[!is.na(sg_found$N_sg_bin)])
      size_dist <- data.table::data.table(
        `Size Range` = names(freq_table),
        Count = as.integer(freq_table),
        Percent = 100 * as.integer(freq_table) / n_found
      )
    }
  }

  # ===========================================================================
  # SECTION 7: FACTOR PRESENCE
  # ===========================================================================

  factor_presence_results <- NULL

  if (n_found > 0) {
    tryCatch({
      factor_presence_results <- summarize_factor_presence_robust(
        sg_found,
        maxk = maxk,
        threshold = 10,
        as_gt = FALSE
      )
    }, error = function(e) {
      warning("Could not summarize factor presence: ", e$message)
      NULL
    })
  }

  # ===========================================================================
  # SECTION 8: ORIGINAL AGREEMENT WITH MAIN ANALYSIS SUBGROUP (FIXED)
  # ===========================================================================

  original_agreement <- NULL

  if (!is.null(original_sg) && n_found > 0) {

    # -------------------------------------------------------------------------
    # FIX: Construct single comparable string from original_sg
    # -------------------------------------------------------------------------
    # original_sg may be a multi-element vector like c("{age>=50}", "{nodes>=3}")
    # We need to create a single comparable string

    orig_sg_char <- as.character(original_sg)

    # Remove NA and empty strings, then collapse to single string
    orig_sg_char <- orig_sg_char[!is.na(orig_sg_char) & orig_sg_char != ""]

    if (length(orig_sg_char) == 0) {
      # No valid original subgroup definition
      original_agreement <- NULL

    } else {
      # Sort and collapse to create canonical form
      orig_sg_combined <- paste(sort(orig_sg_char), collapse = " & ")

      # -----------------------------------------------------------------------
      # Construct comparable strings for each bootstrap iteration
      # -----------------------------------------------------------------------
      # Combine M.1, M.2, etc. columns into single strings per row

      m_cols <- intersect(paste0("M.", seq_len(maxk)), names(sg_found))

      if (length(m_cols) > 0) {
        # Build combined subgroup string for each bootstrap iteration
        boot_subgroups <- vapply(seq_len(nrow(sg_found)), function(i) {
          factors <- character(0)
          for (col in m_cols) {
            val <- sg_found[[col]][i]
            if (!is.na(val) && val != "") {
              factors <- c(factors, val)
            }
          }
          if (length(factors) == 0) {
            return(NA_character_)
          }
          # Sort and collapse to create canonical form (same as original)
          paste(sort(factors), collapse = " & ")
        }, character(1))

        # Count exact matches (comparing single strings)
        matches <- sum(boot_subgroups == orig_sg_combined, na.rm = TRUE)

        # Also count partial matches (at least one factor in common)
        partial_matches <- 0
        if (length(orig_sg_char) > 0) {
          partial_matches <- sum(vapply(boot_subgroups, function(bs) {
            if (is.na(bs)) return(FALSE)
            boot_factors <- strsplit(bs, " & ", fixed = TRUE)[[1]]
            any(boot_factors %in% orig_sg_char)
          }, logical(1)), na.rm = TRUE)
        }

        # Create agreement summary
        original_agreement <- data.table::data.table(
          Metric = c(
            "Total bootstrap iterations",
            "Successful iterations",
            "Failed iterations (no subgroup)",
            "",
            "Original subgroup definition",
            "Exact match with original",
            "Different from original",
            "Partial match (shared factor)"
          ),
          Value = c(
            as.character(nb_boots),
            as.character(n_found),
            as.character(nb_boots - n_found),
            "",
            orig_sg_combined,
            sprintf("%d (%.1f%%)", matches, 100 * matches / n_found),
            sprintf("%d (%.1f%%)", n_found - matches,
                    100 * (n_found - matches) / n_found),
            sprintf("%d (%.1f%%)", partial_matches,
                    100 * partial_matches / n_found)
          )
        )
      }
    }
  }

  # ===========================================================================
  # RETURN COMPILED RESULTS
  # ===========================================================================

  list(
    basic_stats = basic_stats,
    consistency_dist = consistency_dist,
    size_dist = size_dist,
    factor_freq = factor_freq,
    agreement = agreement,
    factor_presence = if (!is.null(factor_presence_results)) {
      factor_presence_results$base_factors
    } else {
      NULL
    },
    factor_presence_specific = if (!is.null(factor_presence_results)) {
      factor_presence_results$specific_factors
    } else {
      NULL
    },
    original_agreement = original_agreement,
    n_found = n_found,
    pct_found = pct_found
  )
}
