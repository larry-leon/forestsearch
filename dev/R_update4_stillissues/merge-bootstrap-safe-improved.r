#' Safe merge handling for bootstrap and regular samples
#'
#' This function properly handles merging data frames in both regular and bootstrap
#' scenarios, preserving duplicate IDs when necessary. It converts inputs to data.frame
#' to ensure consistent behavior across different object types.
#'
#' @param df_main Main data frame (may have duplicate IDs in bootstrap)
#' @param df_flags Flags data frame with treat.recommend (should have unique IDs)
#' @return Merged data frame with preserved row order and duplicates
#' @keywords internal

merge_bootstrap_safe <- function(df_main, df_flags) {
  # ============================================================================
  # STEP 1: INPUT VALIDATION AND CONVERSION
  # ============================================================================

  # Convert to data.frame to ensure consistent behavior
  # This handles data.table, tibble, and other data.frame-like objects
  df_main <- as.data.frame(df_main)
  df_flags <- as.data.frame(df_flags)

  # Validate required columns
  if (!"id" %in% names(df_main)) {
    stop("df_main missing 'id' column")
  }
  if (!"id" %in% names(df_flags)) {
    stop("df_flags missing 'id' column")
  }
  if (!"treat.recommend" %in% names(df_flags)) {
    stop("df_flags missing 'treat.recommend' column")
  }

  # ============================================================================
  # STEP 2: ENSURE UNIQUE IDS IN FLAGS
  # ============================================================================

  # df_flags should have unique IDs (one recommendation per unique ID)
  # If duplicates exist, keep only the first occurrence
  if (any(duplicated(df_flags$id))) {
    warning("Duplicate IDs found in df_flags, keeping first occurrence only")
    df_flags <- df_flags[!duplicated(df_flags$id), ]
  }

  # ============================================================================
  # STEP 3: PRESERVE ORIGINAL ORDER
  # ============================================================================

  # Add row index to preserve original order of df_main
  # This is critical for bootstrap samples where order matters
  df_main$._row_idx_. <- seq_len(nrow(df_main))

  # ============================================================================
  # STEP 4: PERFORM MERGE
  # ============================================================================

  # Use base R merge for reliability
  # all.x = TRUE: keeps all rows from df_main (including duplicates)
  # all.y = FALSE: don't add rows from df_flags that aren't in df_main
  # sort = FALSE: don't sort by ID (preserves relative order within ID groups)

  result <- merge(
    x = df_main,
    y = df_flags[, c("id", "treat.recommend"), drop = FALSE],
    by = "id",
    all.x = TRUE,
    all.y = FALSE,
    sort = FALSE
  )

  # ============================================================================
  # STEP 5: RESTORE ORIGINAL ORDER
  # ============================================================================

  # Sort back to original order using the row index
  result <- result[order(result$._row_idx_.), ]

  # Remove the temporary row index column
  result$._row_idx_. <- NULL

  # Reset row names to sequential
  rownames(result) <- NULL

  # ============================================================================
  # STEP 6: VALIDATION
  # ============================================================================

  # Check that merge preserved the correct number of rows
  if (nrow(result) != nrow(df_main)) {
    warning(sprintf(
      "Merge resulted in different number of rows: expected %d, got %d",
      nrow(df_main),
      nrow(result)
    ))
  }

  # Check for any IDs that didn't get a treatment recommendation
  n_missing <- sum(is.na(result$treat.recommend))
  if (n_missing > 0) {
    warning(sprintf(
      "%d rows in df_main did not match any ID in df_flags (treat.recommend set to NA)",
      n_missing
    ))
  }

  # ============================================================================
  # STEP 7: RETURN AS DATA.FRAME
  # ============================================================================

  # Ensure output is a data.frame (not data.table or tibble)
  return(as.data.frame(result))
}


#' Alternative version with more defensive programming
#'
#' This version includes additional safety checks and handles edge cases
#'
#' @param df_main Main data frame
#' @param df_flags Flags data frame
#' @param id_col Name of ID column (default: "id")
#' @param flag_col Name of flag column (default: "treat.recommend")
#' @return Merged data frame

merge_bootstrap_safe_v2 <- function(df_main, df_flags,
                                    id_col = "id",
                                    flag_col = "treat.recommend") {

  # Convert to data.frame
  df_main <- as.data.frame(df_main)
  df_flags <- as.data.frame(df_flags)

  # Store original attributes
  orig_names <- names(df_main)
  orig_nrow <- nrow(df_main)

  # Validate columns exist
  if (!id_col %in% names(df_main)) {
    stop(sprintf("Column '%s' not found in df_main", id_col))
  }
  if (!id_col %in% names(df_flags)) {
    stop(sprintf("Column '%s' not found in df_flags", id_col))
  }
  if (!flag_col %in% names(df_flags)) {
    stop(sprintf("Column '%s' not found in df_flags", flag_col))
  }

  # Check for column name conflicts
  temp_idx_name <- "._merge_idx_."
  while (temp_idx_name %in% c(names(df_main), names(df_flags))) {
    temp_idx_name <- paste0(temp_idx_name, ".")
  }

  # Handle empty data frames
  if (nrow(df_main) == 0) {
    df_main[[flag_col]] <- numeric(0)
    return(df_main)
  }
  if (nrow(df_flags) == 0) {
    df_main[[flag_col]] <- NA
    return(df_main)
  }

  # Remove duplicates from flags
  df_flags_unique <- df_flags[!duplicated(df_flags[[id_col]]), ]

  # Add temporary index
  df_main[[temp_idx_name]] <- seq_len(nrow(df_main))

  # Prepare flags for merge (only ID and flag columns)
  df_flags_merge <- df_flags_unique[, c(id_col, flag_col), drop = FALSE]

  # Perform merge
  result <- tryCatch({
    merge(
      x = df_main,
      y = df_flags_merge,
      by = id_col,
      all.x = TRUE,
      all.y = FALSE,
      sort = FALSE
    )
  }, error = function(e) {
    stop(sprintf("Merge failed: %s", e$message))
  })

  # Restore order
  result <- result[order(result[[temp_idx_name]]), ]
  result[[temp_idx_name]] <- NULL
  rownames(result) <- NULL

  # Validate result
  if (nrow(result) != orig_nrow) {
    stop(sprintf(
      "Critical error: Merge changed number of rows from %d to %d. ",
      "This should not happen with all.x=TRUE. Check for data issues.",
      orig_nrow, nrow(result)
    ))
  }

  # Ensure column order (flag column at the end)
  col_order <- c(orig_names, flag_col)
  col_order <- unique(col_order)  # Remove duplicates if flag_col already existed
  result <- result[, col_order, drop = FALSE]

  return(as.data.frame(result))
}


#' Test function for merge_bootstrap_safe
#'
#' Tests the merge function with various scenarios
#'
#' @examples
#' test_merge_bootstrap_safe()

test_merge_bootstrap_safe <- function() {
  cat("Testing merge_bootstrap_safe function...\n\n")

  # Test 1: Regular merge (no duplicates)
  cat("Test 1: Regular merge\n")
  df_main <- data.frame(id = 1:5, value = letters[1:5])
  df_flags <- data.frame(id = c(1, 3, 5), treat.recommend = c(0, 1, 0))
  result1 <- merge_bootstrap_safe(df_main, df_flags)
  print(result1)

  # Test 2: Bootstrap with duplicates
  cat("\nTest 2: Bootstrap with duplicate IDs\n")
  df_main <- data.frame(id = c(1, 2, 1, 3, 2), value = letters[1:5])
  df_flags <- data.frame(id = 1:3, treat.recommend = c(0, 1, 0))
  result2 <- merge_bootstrap_safe(df_main, df_flags)
  print(result2)

  # Test 3: Data.table input
  cat("\nTest 3: Data.table input\n")
  if (requireNamespace("data.table", quietly = TRUE)) {
    dt_main <- data.table::data.table(id = c(1, 1, 2), value = 1:3)
    dt_flags <- data.table::data.table(id = 1:2, treat.recommend = c(1, 0))
    result3 <- merge_bootstrap_safe(dt_main, dt_flags)
    print(result3)
    cat("Output class:", class(result3), "\n")
  }

  # Test 4: Missing IDs in flags
  cat("\nTest 4: Missing IDs in flags\n")
  df_main <- data.frame(id = 1:5, value = letters[1:5])
  df_flags <- data.frame(id = c(2, 4), treat.recommend = c(0, 1))
  result4 <- merge_bootstrap_safe(df_main, df_flags)
  print(result4)

  # Test 5: Duplicate IDs in flags (should warn)
  cat("\nTest 5: Duplicate IDs in flags\n")
  df_main <- data.frame(id = 1:3, value = letters[1:3])
  df_flags <- data.frame(id = c(1, 2, 2, 3), treat.recommend = c(0, 1, 0, 1))
  result5 <- merge_bootstrap_safe(df_main, df_flags)
  print(result5)

  cat("\nAll tests completed!\n")
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
