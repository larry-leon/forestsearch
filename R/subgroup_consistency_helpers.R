#' Quietly Evaluate an Expression
#'
#' Suppresses console output for the evaluation of an expression.
#'
#' @param x An R expression to evaluate.
#' @return Invisibly returns the result of the expression.
#' @export

quiet <- function(x) {
  sink(tempfile())
  on.exit(sink())
  invisible(force(x))
}

#' Sort Subgroups by Focus
#'
#' Sorts a data.table of subgroup results according to the specified focus.
#'
#' @param result_new A data.table of subgroup results.
#' @param sg_focus Sorting focus: "hr", "hrMaxSG", "maxSG", "hrMinSG", "minSG".
#' @return A sorted data.table.
#' @importFrom data.table setorder
#' @export

sort_subgroups <- function(result_new, sg_focus) {
  if (sg_focus == "hr") data.table::setorder(result_new, -Pcons, -hr, K)
  if (sg_focus %in% c("hrMaxSG", "maxSG")) data.table::setorder(result_new, -N, -Pcons, K)
  if (sg_focus %in% c("hrMinSG", "minSG")) data.table::setorder(result_new, N, -Pcons, K)
  result_new
}

#' Extract Subgroup Definition and Treatment Recommendation
#'
#' Identifies the subgroup definition and recommends treatment for each subject.
#'
#' @param df The original data.frame.
#' @param top_result The top subgroup result row.
#' @param index.Z Matrix of subgroup indices.
#' @param names.Z Names of covariates.
#' @param confs_labels Label mapping for covariates.
#' @return A list with subgroup definition, labels, recommended treatment, and group id.
#' @export

extract_subgroup <- function(df, top_result, index.Z, names.Z, confs_labels) {
  grp1 <- as.numeric(top_result$g)
  m1 <- as.numeric(top_result$m)
  index1 <- as.numeric(unlist(index.Z[m1, ]))
  this.1 <- names.Z[index1 == 1]
  this.1_label <- unlist(lapply(this.1, FS_labels, confs_labels = confs_labels))
  id.harm <- paste(paste(this.1, collapse = "==1 & "), "==1")
  df.sub <- subset(df, eval(parse(text = id.harm)))
  df.sub$treat.recommend <- 0
  id.noharm <- paste(paste(this.1, collapse = "!=1 | "), "!=1")
  df.subC <- subset(df, eval(parse(text = id.noharm)))
  df.subC$treat.recommend <- 1
  data.found <- rbind(df.sub, df.subC)
  list(
    sg.harm = this.1,
    sg.harm_label = this.1_label,
    df_flag = data.found[, c("id", "treat.recommend")],
    sg.harm.id = grp1
  )
}

#' Plot Subgroup Survival Curves
#'
#' Plots weighted Kaplan-Meier survival curves for a specified subgroup and its complement using the \pkg{weightedSurv} package.
#'
#' @param df.sub A data frame containing data for the subgroup of interest.
#' @param df.subC A data frame containing data for the complement subgroup.
#' @param by.risk Numeric. The risk interval for plotting (passed to \code{weightedSurv::df_counting}).
#' @param confs_labels Named character vector. Covariate label mapping (not used directly in this function, but may be used for labeling).
#' @param this.1_label Character. Label for the subgroup being plotted.
#' @param top_result Data frame row. The top subgroup result row, expected to contain a \code{Pcons} column for consistency criteria.
#'
#' @importFrom weightedsurv df_counting plot_weighted_km
#' @export

plot_subgroup <- function(df.sub, df.subC, by.risk, confs_labels, this.1_label, top_result) {
  if (requireNamespace("weightedsurv", quietly = TRUE)) {
    tte.name <- "Y"
    event.name <- "Event"
    treat.name <- "Treat"
    con.lab <- "control"
    exp.lab <- "treat"
    dfcount <- weightedsurv::df_counting(df.sub, tte.name = tte.name, event.name = event.name, treat.name = treat.name, arms = c(exp.lab, con.lab), by.risk = by.risk)
    dfcountC <- weightedsurv::df_counting(df.subC, tte.name = tte.name, event.name = event.name, treat.name = treat.name, arms = c(exp.lab, con.lab), by.risk = by.risk)
    par(mfrow = c(1, 2))
    weightedsurv::plot_weighted_km(dfcount, conf.int = TRUE, show.logrank = TRUE, put.legend.lr = "topleft", ymax = 1.05, xmed.fraction = 0.65)
    weightedsurv::plot_weighted_km(dfcountC, conf.int = TRUE, show.logrank = TRUE, put.legend.lr = "topleft", ymax = 1.05, xmed.fraction = 0.65)
    cat("*** Subgroup found:", c(this.1_label), "\n")
    cat("% consistency criteria met=", c(top_result$Pcons), "\n")
  }  else {
    message("Package 'weightedsurv' not available: skipping weighted KM plots.")
  }
}




#' Output Subgroup Consistency Results
#'
#' Returns the top subgroup(s) and recommended treatment flags.
#'
#' @param df The original data.frame.
#' @param result_new Sorted subgroup results.
#' @param sg_focus Sorting focus.
#' @param index.Z Matrix of subgroup indices.
#' @param names.Z Names of covariates.
#' @param details Logical; print details.
#' @param plot.sg Logical; plot subgroup curves.
#' @param by.risk Risk interval for plotting.
#' @param confs_labels Covariate label mapping.
#' @return A list with results, subgroup definition, labels, flags, and group id.
#' @importFrom data.table copy
#' @export

sg_consistency_out <- function(df, result_new, sg_focus, index.Z, names.Z, details = FALSE, plot.sg = FALSE, by.risk = 12, confs_labels) {
  result_new <- sort_subgroups(result_new, sg_focus)
  top_result <- result_new[1, ]
  subgroup_info <- extract_subgroup(df, top_result, index.Z, names.Z, confs_labels)
  if (details && plot.sg) {
    plot_subgroup(
      subset(df, eval(parse(text = paste(paste(subgroup_info$sg.harm, collapse = "==1 & "), "==1")))),
      subset(df, eval(parse(text = paste(paste(subgroup_info$sg.harm, collapse = "!=1 | "), "!=1")))),
      by.risk, confs_labels, subgroup_info$sg.harm_label, top_result
    )
  }
  result_out <- data.table::copy(result_new)
  list(
    result = result_out,
    sg.harm = subgroup_info$sg.harm,
    sg.harm_label = subgroup_info$sg.harm_label,
    df_flag = subgroup_info$df_flag,
    sg.harm.id = subgroup_info$sg.harm.id
  )
}

#' Remove Redundant Subgroups
#'
#' Removes redundant subgroups from the results by checking for ties.
#'
#' @param found.hrs Data.table of found subgroups.
#' @return Data.table of non-redundant subgroups.
#' @export

remove_redundant_subgroups <- function(found.hrs) {
  found.new <- found.hrs[order(found.hrs$HR, decreasing = TRUE), ]
  f1.hrs <- found.new[1, ]
  temp <- found.new[-c(1), ]
  temp2 <- as.matrix(found.new[, c("HR", "n", "E", "K", "L(HR)", "U(HR)")])
  id_keep <- which(round(diff(temp2, lag = 1), 6) != 0)
  fkeep.hrs <- temp[id_keep, ]
  na.omit(data.table::data.table(rbind(f1.hrs, fkeep.hrs)))
}


#' Identify Duplicate Rows Based on Near-Identical Numeric Columns
#'
#' Finds rows that are duplicates based on columns 2-10 being nearly identical.
#' Useful for removing redundant subgroups with same HR, sample size, events, etc.
#'
#' @param dt Data.table or data.frame to check for duplicates
#' @param cols_to_check Integer vector. Column indices to check (default: 2:10)
#' @param tolerance Numeric. Tolerance for numeric comparison (default: 1e-6)
#' @param keep Character. Which duplicates to keep: "first", "last", "none", or "all" (default: "first")
#' @param return_indices Logical. Return row indices instead of filtered data (default: FALSE)
#'
#' @return Data.table with duplicates removed (or indices if return_indices=TRUE)
#' @importFrom data.table data.table copy
#' @export
identify_near_duplicates <- function(dt,
                                     cols_to_check = 2:10,
                                     tolerance = 1e-6,
                                     keep = "first",
                                     return_indices = FALSE) {

  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' required")
  }

  # Ensure it's a data.table
  if (!is.data.table(dt)) {
    dt <- data.table::as.data.table(dt)
  } else {
    dt <- data.table::copy(dt)
  }

  # Validate inputs
  if (max(cols_to_check) > ncol(dt)) {
    stop("Column indices exceed number of columns in data")
  }

  # Get the columns to compare
  compare_cols <- names(dt)[cols_to_check]

  # Round numeric columns to handle near-identical values
  dt_rounded <- data.table::copy(dt)
  for (col in compare_cols) {
    if (is.numeric(dt_rounded[[col]])) {
      # Round to eliminate floating point differences
      dt_rounded[[col]] <- round(dt_rounded[[col]] / tolerance) * tolerance
    }
  }

  # Create a composite key for comparison
  dt_rounded[, dup_key := do.call(paste, c(.SD, sep = "_")), .SDcols = compare_cols]

  # Identify duplicates
  dt_rounded[, dup_group := .GRP, by = dup_key]
  dt_rounded[, dup_count := .N, by = dup_key]
  dt_rounded[, dup_rank := seq_len(.N), by = dup_key]

  # Determine which rows to keep
  if (keep == "first") {
    keep_rows <- dt_rounded$dup_rank == 1
  } else if (keep == "last") {
    keep_rows <- dt_rounded$dup_rank == dt_rounded$dup_count
  } else if (keep == "none") {
    keep_rows <- dt_rounded$dup_count == 1
  } else if (keep == "all") {
    keep_rows <- rep(TRUE, nrow(dt_rounded))
  } else {
    stop("keep must be one of: 'first', 'last', 'none', 'all'")
  }

  if (return_indices) {
    # Return information about duplicates
    return(list(
      keep_indices = which(keep_rows),
      duplicate_groups = dt_rounded[, .(grp, dup_group, dup_count, dup_rank)],
      n_duplicates = sum(!keep_rows)
    ))
  } else {
    # Return filtered data
    return(dt[keep_rows])
  }
}

#' Remove Near-Duplicate Subgroups (Warning-Free Version)
#'
#' Simplified version that avoids data.table modification warnings
#'
#' @param hr_subgroups Data.table of subgroup results
#' @param tolerance Numeric. Tolerance for numeric comparison (default: 0.001)
#' @param details Logical. Print details about removed duplicates
#'
#' @return Data.table with near-duplicate rows removed
#' @export
remove_near_duplicate_subgroups <- function(hr_subgroups,
                                            tolerance = 0.001,
                                            details = FALSE) {

  # Convert to regular data.frame to avoid data.table copy warnings
  df <- as.data.frame(hr_subgroups)

  # Columns to check: K, n, E, d1, m1, m0, HR, L(HR), U(HR)
  cols_to_check <- 2:10

  # Round numeric columns
  df_rounded <- df
  for (i in cols_to_check) {
    if (is.numeric(df_rounded[, i])) {
      df_rounded[, i] <- round(df_rounded[, i] / tolerance) * tolerance
    }
  }

  # Create composite key
  key_cols <- df_rounded[, cols_to_check, drop = FALSE]
  dup_key <- apply(key_cols, 1, function(x) paste(x, collapse = "_"))

  # Find duplicates (keep first occurrence)
  keep_rows <- !duplicated(dup_key)

  n_removed <- sum(!keep_rows)

  if (details && n_removed > 0) {
    cat("Removed", n_removed, "near-duplicate subgroups\n")
    cat("Original rows:", nrow(hr_subgroups), "\n")
    cat("After removal:", sum(keep_rows), "\n")
  }

  # Return as data.table
  return(data.table::as.data.table(hr_subgroups[keep_rows, ]))
}



#' Show Duplicate Subgroups (Diagnostic Function)
#'
#' Display which subgroups are near-duplicates of each other
#'
#' @param hr_subgroups Data.table of subgroup results
#' @param tolerance Numeric. Tolerance for numeric comparison (default: 0.001)
#'
#' @return Data.table showing duplicate groups
#' @export
show_duplicate_subgroups <- function(hr_subgroups, tolerance = 0.001) {

  result <- identify_near_duplicates(
    dt = hr_subgroups,
    cols_to_check = 2:10,
    tolerance = tolerance,
    keep = "all",
    return_indices = TRUE
  )

  dup_info <- result$duplicate_groups
  dup_info <- dup_info[dup_count > 1][order(dup_group, dup_rank)]

  if (nrow(dup_info) == 0) {
    cat("No near-duplicate subgroups found\n")
    return(NULL)
  }

  cat("Found", length(unique(dup_info$dup_group)), "groups of duplicates\n")
  cat("Total duplicate rows:", nrow(dup_info), "\n\n")

  # Show the actual duplicates with their key columns
  hr_with_dup <- cbind(hr_subgroups[dup_info$grp],
                       dup_info[, .(dup_group, dup_rank)])

  return(hr_with_dup[order(dup_group, dup_rank)])
}


#' Set up parallel processing for subgroup consistency
#'
#' Sets up parallel processing using the specified approach and number of workers.
#'
#' @param parallel_args List with \code{plan} (character), \code{workers} (integer), and \code{show_message} (logical).
#' @return None. Sets up parallel backend.
#' @importFrom future plan
#' @export

setup_parallel_SGcons <- function(parallel_args = list(plan = "multisession", workers = 4, show_message = TRUE)) {
  plan_type <- parallel_args$plan
  n_workers <- parallel_args$workers
  show_message <- parallel_args$show_message
  if (is.null(plan_type)) stop("parallel_args$plan must be specified.")

  allowed_plans <- c("multisession", "multicore", "callr", "sequential")

  if (!plan_type %in% allowed_plans) {
    stop("plan_type must be one of: ", paste(allowed_plans, collapse = ", "))
  }

  max_cores <- parallel::detectCores()
  if (is.null(n_workers) || !is.numeric(n_workers) || n_workers < 1) {
    n_workers <- 1
  } else {
    n_workers <- min(n_workers, max_cores)
  }
  if (plan_type == "multisession") {
    plan(multisession, workers = n_workers)
 if(show_message) message("Parallel plan: multisession with ", n_workers, " workers.")
  } else if (plan_type == "multicore") {
    plan(multicore, workers = n_workers)
  if(show_message)  message("Parallel plan: multicore with ", n_workers, " workers.")
  } else if (plan_type == "callr") {
    if (!requireNamespace("future.callr", quietly = TRUE)) {
      stop("The 'future.callr' package is required for the callr approach.")
    }
    plan(future.callr::callr, workers = n_workers)
if(show_message) message("Parallel plan: callr with ", n_workers, " workers.")
  }
  else if (plan_type == "sequential") {
  plan(sequential)
  if(show_message)  message("Sequential plan")
  }
  else {
    stop("Unknown parallel plan: ", plan_type)
  }
}

