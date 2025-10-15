# Fixed implementation for the merge section in forestsearch function

if (!is.null(grp.consistency$sg.harm)) {
  sg.harm <- grp.consistency$sg.harm
  
  if (details) {
    cat("\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("✓ SUBGROUP IDENTIFIED AND VALIDATED!\n")
    cat(rep("=", 70), "\n", sep = "")
    cat("Definition:", paste(sg.harm, collapse = " & "), "\n")
    cat("Total time:", round(t.min_all, 2), "minutes\n\n")
  }
  
  # Merge treatment recommendations to datasets
  # data containing id and treatment flag
  temp <- grp.consistency$df_flag
  
  # ============================================================================
  # FIX: Use safe merge function that handles bootstrap duplicates
  # ============================================================================
  
  # Helper function for safe merging (should be defined at package level)
  merge_safe_local <- function(df_main, df_flags) {
    # Convert to data.frames to ensure consistent behavior
    df_main <- as.data.frame(df_main)
    df_flags <- as.data.frame(df_flags)
    
    # Validate inputs
    if (!"id" %in% names(df_main)) stop("df_main missing 'id' column")
    if (!"id" %in% names(df_flags)) stop("df_flags missing 'id' column")
    
    # Ensure df_flags has unique IDs
    if (any(duplicated(df_flags$id))) {
      df_flags <- df_flags[!duplicated(df_flags$id), ]
    }
    
    # Add row index to preserve order (especially important for bootstrap)
    df_main$._row_idx_. <- seq_len(nrow(df_main))
    
    # Perform merge
    result <- merge(
      x = df_main,
      y = df_flags,
      by = "id",
      all.x = TRUE,
      all.y = FALSE,
      sort = FALSE
    )
    
    # Restore original order
    result <- result[order(result$._row_idx_.), ]
    result$._row_idx_. <- NULL
    rownames(result) <- NULL
    
    # Validate
    if (nrow(result) != nrow(df_main)) {
      warning(sprintf("Merge changed row count: %d -> %d", 
                     nrow(df_main), nrow(result)))
    }
    
    return(as.data.frame(result))
  }
  
  # ============================================================================
  # Apply safe merge to each dataset
  # ============================================================================
  
  # Merge to analysis data
  df.est_out <- merge_safe_local(df, temp)
  
  # Debug info for bootstrap context
  if (details) {
    n_unique_ids <- length(unique(df$id))
    n_total_rows <- nrow(df)
    if (n_unique_ids < n_total_rows) {
      cat("Note: Bootstrap sample detected\n")
      cat("  Unique IDs:", n_unique_ids, "\n")
      cat("  Total rows:", n_total_rows, "\n")
      cat("  Duplicates:", n_total_rows - n_unique_ids, "\n\n")
    }
  }
  
  # Handle df.predict
  if (!is.null(df.predict)) {
    # Check if df.predict has same IDs as analysis data
    predict_ids_in_analysis <- all(df.predict$id %in% df$id)
    
    if (predict_ids_in_analysis) {
      # IDs match - can merge directly
      df.predict_out <- merge_safe_local(df.predict, temp)
    } else {
      # Different IDs (e.g., test set) - need to apply subgroup definition
      if (details) {
        cat("Note: df.predict has different IDs, applying subgroup definition directly\n")
      }
      df.predict_out <- get_dfpred(
        df.predict = df.predict,
        sg.harm = grp.consistency$sg.harm,
        version = 2
      )
    }
  }
  
  # Handle df.test (always different IDs)
  if (!is.null(df.test)) {
    df.test_out <- get_dfpred(
      df.predict = df.test,
      sg.harm = grp.consistency$sg.harm,
      version = 2
    )
  }
  
} else {
  # No subgroup found
  if (details) {
    cat("\n✗ No subgroup met consistency criteria\n")
    cat("  Possible reasons:\n")
    cat("    - Consistency rate below threshold (", pconsistency.threshold, ")\n")
    cat("    - HR in splits below threshold (", hr.consistency, ")\n\n")
  }
}


# ============================================================================
# Alternative: Complete replacement with better error handling
# ============================================================================

#' Process subgroup consistency results and merge with data
#'
#' This function handles the merging of subgroup treatment recommendations
#' back to the original and prediction datasets, with proper handling of
#' bootstrap samples that contain duplicate IDs.
#'
#' @param grp.consistency Subgroup consistency results
#' @param df Analysis data frame (may have duplicate IDs in bootstrap)
#' @param df.predict Prediction data frame (optional)
#' @param df.test Test data frame (optional)
#' @param details Logical for detailed output
#' @param t.min_all Total computation time
#' @param pconsistency.threshold Consistency threshold
#' @param hr.consistency HR consistency threshold
#' @return List with processed datasets

process_subgroup_results <- function(grp.consistency, df, df.predict = NULL, 
                                     df.test = NULL, details = FALSE,
                                     t.min_all = NULL, 
                                     pconsistency.threshold = 0.9,
                                     hr.consistency = 1.0) {
  
  # Initialize outputs
  sg.harm <- NULL
  df.est_out <- NULL
  df.predict_out <- NULL
  df.test_out <- NULL
  
  if (!is.null(grp.consistency$sg.harm)) {
    sg.harm <- grp.consistency$sg.harm
    
    if (details) {
      cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
      cat("✓ SUBGROUP IDENTIFIED AND VALIDATED!\n")
      cat(paste(rep("=", 70), collapse = ""), "\n")
      cat("Definition:", paste(sg.harm, collapse = " & "), "\n")
      if (!is.null(t.min_all)) {
        cat("Total time:", round(t.min_all, 2), "minutes\n")
      }
      cat("\n")
    }
    
    # Get treatment recommendations
    temp <- grp.consistency$df_flag
    
    # Check for bootstrap context
    is_bootstrap <- length(unique(df$id)) < nrow(df)
    
    if (is_bootstrap && details) {
      cat("Bootstrap context detected:\n")
      cat("  Unique IDs:", length(unique(df$id)), "\n")
      cat("  Total rows:", nrow(df), "\n")
      cat("  Bootstrap duplicates:", sum(duplicated(df$id)), "\n\n")
    }
    
    # Safe merge for analysis data
    df.est_out <- tryCatch({
      merge_safe(df, temp)
    }, error = function(e) {
      warning("Merge failed for df.est, using fallback: ", e$message)
      # Fallback: apply subgroup definition directly
      get_dfpred(df.predict = df, sg.harm = grp.consistency$sg.harm, version = 2)
    })
    
    # Handle prediction dataset
    if (!is.null(df.predict)) {
      # Determine if IDs overlap
      ids_overlap <- any(df.predict$id %in% unique(df$id))
      
      if (ids_overlap) {
        # Try merge first
        df.predict_out <- tryCatch({
          merge_safe(df.predict, temp)
        }, error = function(e) {
          if (details) {
            cat("Note: Merge failed for df.predict, applying subgroup definition\n")
          }
          get_dfpred(df.predict = df.predict, 
                    sg.harm = grp.consistency$sg.harm, 
                    version = 2)
        })
      } else {
        # No ID overlap - apply definition directly
        if (details) {
          cat("df.predict has independent IDs, applying subgroup definition\n")
        }
        df.predict_out <- get_dfpred(df.predict = df.predict,
                                     sg.harm = grp.consistency$sg.harm,
                                     version = 2)
      }
    }
    
    # Handle test dataset (always independent IDs)
    if (!is.null(df.test)) {
      df.test_out <- get_dfpred(df.predict = df.test,
                               sg.harm = grp.consistency$sg.harm,
                               version = 2)
    }
    
    # Validate results
    if (!is.null(df.est_out)) {
      n_assigned <- sum(!is.na(df.est_out$treat.recommend))
      if (details) {
        cat("Treatment recommendations assigned:\n")
        cat("  df.est:", n_assigned, "of", nrow(df.est_out), "rows\n")
        if (!is.null(df.predict_out)) {
          n_pred <- sum(!is.na(df.predict_out$treat.recommend))
          cat("  df.predict:", n_pred, "of", nrow(df.predict_out), "rows\n")
        }
        if (!is.null(df.test_out)) {
          n_test <- sum(!is.na(df.test_out$treat.recommend))
          cat("  df.test:", n_test, "of", nrow(df.test_out), "rows\n")
        }
      }
    }
    
  } else {
    # No subgroup found
    if (details) {
      cat("\n✗ No subgroup met consistency criteria\n")
      cat("  Possible reasons:\n")
      cat("    - Consistency rate below threshold (", pconsistency.threshold, ")\n")
      cat("    - HR in splits below threshold (", hr.consistency, ")\n\n")
    }
  }
  
  return(list(
    sg.harm = sg.harm,
    df.est = df.est_out,
    df.predict = df.predict_out,
    df.test = df.test_out
  ))
}


# ============================================================================
# Define merge_safe at package level (in a utilities file)
# ============================================================================

#' Safe merge for ForestSearch data
#'
#' Handles merging with proper support for bootstrap samples containing
#' duplicate IDs from sampling with replacement.
#'
#' @param df_main Main data frame (may have duplicate IDs)
#' @param df_flags Flags data frame with treatment recommendations
#' @return Merged data frame preserving all rows from df_main

merge_safe <- function(df_main, df_flags) {
  # Convert to data.frame for consistency
  df_main <- as.data.frame(df_main)
  df_flags <- as.data.frame(df_flags)
  
  # Ensure unique IDs in flags
  if (any(duplicated(df_flags$id))) {
    df_flags <- df_flags[!duplicated(df_flags$id), ]
  }
  
  # Preserve row order
  df_main$._merge_order_. <- seq_len(nrow(df_main))
  
  # Merge
  result <- merge(
    x = df_main,
    y = df_flags,
    by = "id",
    all.x = TRUE,
    all.y = FALSE,
    sort = FALSE
  )
  
  # Restore order
  result <- result[order(result$._merge_order_.), ]
  result$._merge_order_. <- NULL
  rownames(result) <- NULL
  
  # Validate
  if (nrow(result) != nrow(df_main)) {
    stop(sprintf("Critical: Merge changed rows from %d to %d", 
                nrow(df_main), nrow(result)))
  }
  
  return(result)
}