#' Fixed Version of summarize_bootstrap_subgroups
#'
#' This version fixes the data.table ordering error by ensuring proper
#' data structure handling throughout.
#'
#' @param results Data.table or data.frame. Bootstrap results with subgroup characteristics
#' @param nb_boots Integer. Total number of bootstrap iterations
#' @param original_sg Character vector. Original subgroup definition (M.1, M.2 from fs.est)
#' @param maxk Integer. Maximum number of factors allowed
#'
#' @return List with multiple summary components
#' @importFrom data.table data.table .N setnames setcolorder
#' @export

summarize_bootstrap_subgroups <- function(results, nb_boots,
                                          original_sg = NULL,
                                          maxk = 2) {

  # =========================================================================
  # FIX 1: Ensure results is a proper data.table
  # =========================================================================

  # Convert to data.table if it's not already
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

  # =========================================================================
  # SECTION 1: BASIC STATISTICS
  # =========================================================================

  # Filter to successful iterations (where subgroup was found)
  # FIX: Use proper data.table syntax
  if ("Pcons" %in% names(results)) {
    sg_found <- results[!is.na(results$Pcons), ]
  } else {
    warning("Column 'Pcons' not found in results")
    sg_found <- results
  }

  n_found <- nrow(sg_found)
  pct_found <- 100 * n_found / nb_boots

  if (n_found == 0) {
    warning("No subgroups identified in any bootstrap iteration")
    return(list(
      basic_stats = NULL,
      consistency_dist = NULL,
      size_dist = NULL,
      factor_freq = NULL,
      agreement = NULL
    ))
  }

  # =========================================================================
  # SECTION 2: BASIC STATISTICS TABLE
  # =========================================================================

  # Create basic stats with safe extraction
  basic_stats_list <- list(
    Metric = character(),
    Value = character()
  )

  # Add basic counts
  basic_stats_list$Metric <- c(
    "Total bootstrap iterations",
    "Subgroups identified",
    "Success rate (%)"
  )
  basic_stats_list$Value <- c(
    as.character(nb_boots),
    as.character(n_found),
    sprintf("%.1f%%", pct_found)
  )

  # Add Pcons statistics if available
  if ("Pcons" %in% names(sg_found) && sum(!is.na(sg_found$Pcons)) > 0) {
    basic_stats_list$Metric <- c(basic_stats_list$Metric,
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
    basic_stats_list$Value <- c(basic_stats_list$Value,
                                "",
                                "",
                                sprintf("%.3f", mean(sg_found$Pcons, na.rm = TRUE)),
                                sprintf("%.3f", median(sg_found$Pcons, na.rm = TRUE)),
                                sprintf("%.3f", sd(sg_found$Pcons, na.rm = TRUE)),
                                sprintf("%.3f", min(sg_found$Pcons, na.rm = TRUE)),
                                sprintf("%.3f", max(sg_found$Pcons, na.rm = TRUE)),
                                sprintf("%.3f", quantile(sg_found$Pcons, 0.25, na.rm = TRUE)),
                                sprintf("%.3f", quantile(sg_found$Pcons, 0.75, na.rm = TRUE))
    )
  }

  # Add hr_sg statistics if available
  if ("hr_sg" %in% names(sg_found) && sum(!is.na(sg_found$hr_sg)) > 0) {
    basic_stats_list$Metric <- c(basic_stats_list$Metric,
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
    basic_stats_list$Value <- c(basic_stats_list$Value,
                                "",
                                "",
                                sprintf("%.2f", mean(sg_found$hr_sg, na.rm = TRUE)),
                                sprintf("%.2f", median(sg_found$hr_sg, na.rm = TRUE)),
                                sprintf("%.2f", sd(sg_found$hr_sg, na.rm = TRUE)),
                                sprintf("%.2f", min(sg_found$hr_sg, na.rm = TRUE)),
                                sprintf("%.2f", max(sg_found$hr_sg, na.rm = TRUE)),
                                sprintf("%.2f", quantile(sg_found$hr_sg, 0.25, na.rm = TRUE)),
                                sprintf("%.2f", quantile(sg_found$hr_sg, 0.75, na.rm = TRUE))
    )
  }

  # Convert to data.table
  basic_stats <- data.table::data.table(basic_stats_list)

  # =========================================================================
  # SECTION 3: FACTOR FREQUENCY TABLE - FIXED
  # =========================================================================

  # Count frequency of each factor appearing
  factor_cols <- paste0("M.", 1:maxk)
  factor_freq_list <- list()

  for (i in 1:maxk) {
    col <- paste0("M.", i)
    if (col %in% names(sg_found)) {
      # FIX: Use proper data.table subsetting
      valid_rows <- !is.na(sg_found[[col]]) & sg_found[[col]] != ""
      if (sum(valid_rows) > 0) {
        # Create frequency table
        factor_vals <- sg_found[[col]][valid_rows]
        freq_table <- table(factor_vals)

        # Convert to data.table
        freq <- data.table::data.table(
          Factor = names(freq_table),
          N = as.integer(freq_table),
          Position = paste0("M.", i),
          Percent = 100 * as.integer(freq_table) / n_found
        )

        # Sort by N
        freq <- freq[order(-N)]
        factor_freq_list[[i]] <- freq
      }
    }
  }

  # FIX: Properly combine data.tables
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

  # =========================================================================
  # SECTION 4: SUBGROUP DEFINITION AGREEMENT - FIXED
  # =========================================================================

  agreement <- NULL

  if (maxk == 1) {
    if ("M.1" %in% names(sg_found)) {
      valid_rows <- !is.na(sg_found$M.1) & sg_found$M.1 != ""
      if (sum(valid_rows) > 0) {
        freq_table <- table(sg_found$M.1[valid_rows])
        agreement <- data.table::data.table(
          Subgroup = names(freq_table),
          K_sg = 1,
          N = as.integer(freq_table)
        )
      }
    }
  } else if (maxk == 2) {
    if (all(c("M.1", "M.2") %in% names(sg_found))) {
      # Create subgroup combinations
      valid_rows <- !is.na(sg_found$M.1) & sg_found$M.1 != ""
      if (sum(valid_rows) > 0) {
        sg_temp <- sg_found[valid_rows, ]
        sg_temp$Subgroup <- ifelse(
          is.na(sg_temp$M.2) | sg_temp$M.2 == "",
          sg_temp$M.1,
          paste(sg_temp$M.1, "&", sg_temp$M.2)
        )

        # Count combinations
        if ("K_sg" %in% names(sg_temp)) {
          agreement <- sg_temp[, .N, by = .(Subgroup, K_sg)]
        } else {
          agreement <- sg_temp[, .N, by = Subgroup]
          agreement$K_sg <- NA
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

  # =========================================================================
  # SECTION 5: CONSISTENCY DISTRIBUTION - FIXED
  # =========================================================================

  consistency_dist <- NULL

  if ("Pcons" %in% names(sg_found) && sum(!is.na(sg_found$Pcons)) > 0) {
    # Create bins
    sg_found$Pcons_bin <- cut(
      sg_found$Pcons,
      breaks = c(0, 0.5, 0.7, 0.8, 0.9, 0.95, 1.0),
      labels = c("<0.5", "0.5-0.7", "0.7-0.8", "0.8-0.9", "0.9-0.95", "≥0.95"),
      include.lowest = TRUE
    )

    # Count by bin
    if (sum(!is.na(sg_found$Pcons_bin)) > 0) {
      freq_table <- table(sg_found$Pcons_bin[!is.na(sg_found$Pcons_bin)])
      consistency_dist <- data.table::data.table(
        `Consistency Range` = names(freq_table),
        Count = as.integer(freq_table),
        Percent = 100 * as.integer(freq_table) / n_found
      )
    }
  }

  # =========================================================================
  # SECTION 6: SIZE DISTRIBUTION - FIXED
  # =========================================================================

  size_dist <- NULL

  if ("N_sg" %in% names(sg_found) && sum(!is.na(sg_found$N_sg)) > 0) {
    # Create bins
    size_breaks <- c(0, 50, 100, 150, 200, 300, Inf)
    size_labels <- c("<50", "50-99", "100-149", "150-199", "200-299", "≥300")

    sg_found$N_sg_bin <- cut(
      sg_found$N_sg,
      breaks = size_breaks,
      labels = size_labels,
      include.lowest = TRUE
    )

    # Count by bin
    if (sum(!is.na(sg_found$N_sg_bin)) > 0) {
      freq_table <- table(sg_found$N_sg_bin[!is.na(sg_found$N_sg_bin)])
      size_dist <- data.table::data.table(
        `Size Range` = names(freq_table),
        Count = as.integer(freq_table),
        Percent = 100 * as.integer(freq_table) / n_found
      )
    }
  }

  # =========================================================================
  # SECTION 7: FACTOR PRESENCE - FIXED
  # =========================================================================

  factor_presence_results <- NULL

  if (n_found > 0) {
    # Try to summarize factor presence
    tryCatch({
      factor_presence_results <- summarize_factor_presence_robust(sg_found, maxk = maxk, threshold = 20)
    }, error = function(e) {
      warning("Could not summarize factor presence: ", e$message)
      NULL
    })
  }

  # list(
  #   basic_stats = basic_stats,
  #   consistency_dist = consistency_dist,
  #   size_dist = size_dist,
  #   factor_freq = factor_freq,
  #   agreement = agreement,
  #   factor_presence = if (!is.null(factor_presence_results)) factor_presence_results$base_factors else NULL,
  #   factor_presence_specific = if (!is.null(factor_presence_results)) factor_presence_results$specific_factors else NULL,
  #   original_agreement = NULL,  # Simplified for now
  #   n_found = n_found,
  #   pct_found = pct_found
  # )
  #

  # =========================================================================
  # SECTION 8: ORIGINAL AGREEMENT WITH MAIN ANALYSIS SUBGROUP
  # =========================================================================

  original_agreement <- NULL

  if (!is.null(original_sg) && n_found > 0) {
    # Try to find the bootstrap subgroup column
    subgroup_col <- NULL

    # Check for common subgroup column names
    for (possible_col in c("M.1", "Subgroup", "H_bootstrap", "H", "H_star")) {
      if (possible_col %in% names(sg_found)) {
        subgroup_col <- possible_col
        break
      }
    }

    # Calculate agreement if subgroup column found
    if (!is.null(subgroup_col)) {
      # Extract original subgroup as character
      orig_sg_char <- as.character(original_sg)

      # Count exact matches
      matches <- sum(as.character(sg_found[[subgroup_col]]) == orig_sg_char, na.rm = TRUE)

      # Create agreement summary
      original_agreement <- data.table::data.table(
        Metric = c(
          "Total bootstrap iterations",
          "Successful iterations",
          "Failed iterations (no subgroup)",
          "Exact match with original",
          "Different from original"
        ),
        Value = c(
          as.character(nb_boots),
          as.character(n_found),
          as.character(nb_boots - n_found),
          sprintf("%d (%.1f%%)", matches, 100 * matches / n_found),
          sprintf("%d (%.1f%%)", n_found - matches, 100 * (n_found - matches) / n_found)
        ),
        stringsAsFactors = FALSE
      )
    }
  }

  # =========================================================================
  # RETURN COMPILED RESULTS
  # =========================================================================

  list(
    basic_stats = basic_stats,
    consistency_dist = consistency_dist,
    size_dist = size_dist,
    factor_freq = factor_freq,
    agreement = agreement,
    factor_presence = if (!is.null(factor_presence_results)) factor_presence_results$base_factors else NULL,
    factor_presence_specific = if (!is.null(factor_presence_results)) factor_presence_results$specific_factors else NULL,
    original_agreement = original_agreement,  # NOW CALCULATED!
    n_found = n_found,
    pct_found = pct_found
  )
}



