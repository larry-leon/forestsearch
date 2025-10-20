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
#' @importFrom weightedSurv df_counting plot_weighted_km
#' @export

plot_subgroup <- function(df.sub, df.subC, by.risk, confs_labels, this.1_label, top_result) {
  if (requireNamespace("weightedSurv", quietly = TRUE)) {
    tte.name <- "Y"
    event.name <- "Event"
    treat.name <- "Treat"
    con.lab <- "control"
    exp.lab <- "treat"
    dfcount <- weightedSurv::df_counting(df.sub, tte.name = tte.name, event.name = event.name, treat.name = treat.name, arms = c(exp.lab, con.lab), by.risk = by.risk)
    dfcountC <- weightedSurv::df_counting(df.subC, tte.name = tte.name, event.name = event.name, treat.name = treat.name, arms = c(exp.lab, con.lab), by.risk = by.risk)
    par(mfrow = c(1, 2))
    weightedSurv::plot_weighted_km(dfcount, conf.int = TRUE, show.logrank = TRUE, put.legend.lr = "topleft", ymax = 1.05, xmed.fraction = 0.65)
    weightedSurv::plot_weighted_km(dfcountC, conf.int = TRUE, show.logrank = TRUE, put.legend.lr = "topleft", ymax = 1.05, xmed.fraction = 0.65)
    cat("*** Subgroup found:", c(this.1_label), "\n")
    cat("% consistency criteria met=", c(top_result$Pcons), "\n")
  }  else {
    message("Package 'weightedSurv' not available: skipping weighted KM plots.")
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

#' Get Hazard Ratio from Split Data
#'
#' Calculates the hazard ratio from a split data set using Cox regression.
#'
#' @param df Data frame with survival data.
#' @return Numeric hazard ratio or NA if error.
#' @importFrom survival coxph Surv
get_split_hr_legacy <- function(df, cox_initial = log(1)) {
  hr <- try((summary(coxph(Surv(Y, Event) ~ Treat, data = df, robust = FALSE, init = cox_initial))$conf.int), TRUE)
  if (inherits(hr, "try-error")) return(NA_real_)
  hr[1,1]
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

#' Subgroup Consistency Evaluation (WITH COMPREHENSIVE ERROR CHECKS)
#'
#' Evaluates consistency of subgroups found in a survival analysis, using random splits and hazard ratio criteria.
#'
#' @param df The original data.frame.
#' @param hr.subgroups Data.table of subgroup hazard ratio results.
#' @param hr.threshold Minimum hazard ratio for subgroup inclusion.
#' @param hr.consistency Minimum hazard ratio for consistency in splits.
#' @param pconsistency.threshold Minimum proportion of splits meeting consistency.
#' @param pconsistency.digits Significant digits for pconsistency.threshold
#' @param m1.threshold Maximum median survival for treatment arm.
#' @param n.splits Number of random splits for consistency evaluation.
#' @param details Logical; print details.
#' @param stop.threshold Consistency threshold for early stopping.
#' @param by.risk Risk interval for plotting.
#' @param plot.sg Logical; plot subgroup curves.
#' @param maxk Maximum number of covariates in a subgroup.
#' @param Lsg Number of covariates.
#' @param confs_labels Covariate label mapping.
#' @param sg_focus Sorting focus.
#' @param stop_Kgroups Maximum number of subgroups to identify.
#' @param checking Logical; for debugging.
#' @param parallel_args List for parallel processing.
#' @return A list with results for top subgroups and recommended treatment flags.
#' @importFrom data.table copy
#' @importFrom survival coxph Surv
#' @importFrom future.apply future_lapply
#' @export

subgroup.consistency <- function(df, hr.subgroups,
                                 hr.threshold = 1.0,
                                 hr.consistency = 1.0,
                                 pconsistency.threshold = 0.9,
                                 m1.threshold = Inf,
                                 n.splits = 100,
                                 details = FALSE,
                                 stop.threshold = 1.1,
                                 by.risk = 12,
                                 plot.sg = FALSE,
                                 maxk = 7,
                                 Lsg,
                                 confs_labels,
                                 sg_focus = "hr",
                                 stop_Kgroups = 10,
                                 pconsistency.digits = 2,
                                 checking = FALSE,
                                 parallel_args = list(NULL)) {

  # =========================================================================
  # SECTION 1: INPUT VALIDATION
  # =========================================================================

  # Check required inputs exist
  if (missing(df) || !is.data.frame(df)) {
    stop("'df' must be provided and must be a data.frame")
  }

  if (missing(hr.subgroups) || is.null(hr.subgroups)) {
    stop("'hr.subgroups' must be provided")
  }

  if (missing(Lsg) || !is.numeric(Lsg) || Lsg < 1) {
    stop("'Lsg' must be a positive integer representing number of covariates")
  }

  if (missing(confs_labels) || !is.character(confs_labels)) {
    stop("'confs_labels' must be a character vector of covariate labels")
  }

  # Check data.table format
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required")
  }

  if (!data.table::is.data.table(hr.subgroups)) {
    if (is.data.frame(hr.subgroups)) {
      hr.subgroups <- data.table::as.data.table(hr.subgroups)
      if (details) message("Converting hr.subgroups to data.table")
    } else {
      stop("'hr.subgroups' must be a data.frame or data.table")
    }
  }

  # Validate required columns in df
  required_cols <- c("Y", "Event", "Treat", "id")
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("df is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Validate required columns in hr.subgroups
  required_hr_cols <- c("grp", "K", "n", "E", "d1", "m1", "m0", "HR", "L(HR)", "U(HR)")
  missing_hr_cols <- setdiff(required_hr_cols, names(hr.subgroups))
  if (length(missing_hr_cols) > 0) {
    stop("hr.subgroups is missing required columns: ", paste(missing_hr_cols, collapse = ", "))
  }

  # Check numeric parameters
  if (!is.numeric(hr.threshold) || hr.threshold < 0) {
    stop("'hr.threshold' must be a non-negative numeric value")
  }

  if (!is.numeric(hr.consistency) || hr.consistency < 0) {
    stop("'hr.consistency' must be a non-negative numeric value")
  }

  if (!is.numeric(pconsistency.threshold) ||
      pconsistency.threshold < 0 || pconsistency.threshold > 1) {
    stop("'pconsistency.threshold' must be between 0 and 1")
  }

  if (!is.numeric(n.splits) || n.splits < 1) {
    stop("'n.splits' must be a positive integer")
  }
  n.splits <- as.integer(n.splits)

  if (!is.numeric(stop.threshold) || stop.threshold < 0 || stop.threshold > 1.1) {
    warning("'stop.threshold' is typically between 0 and 1.1")
  }

  if (!is.numeric(maxk) || maxk < 1 || maxk > 7) {
    stop("'maxk' must be between 1 and 7")
  }
  maxk <- as.integer(maxk)

  if (!is.numeric(stop_Kgroups) || stop_Kgroups < 1) {
    stop("'stop_Kgroups' must be a positive integer")
  }
  stop_Kgroups <- as.integer(stop_Kgroups)

  if (!is.numeric(pconsistency.digits) || pconsistency.digits < 0) {
    stop("'pconsistency.digits' must be a non-negative integer")
  }
  pconsistency.digits <- as.integer(pconsistency.digits)

  # Check logical parameters
  if (!is.logical(details)) {
    stop("'details' must be logical (TRUE/FALSE)")
  }

  if (!is.logical(plot.sg)) {
    stop("'plot.sg' must be logical (TRUE/FALSE)")
  }

  if (!is.logical(checking)) {
    stop("'checking' must be logical (TRUE/FALSE)")
  }

  # Validate sg_focus
  valid_sg_focus <- c("hr", "hrMaxSG", "hrMinSG", "maxSG", "minSG")
  if (!sg_focus %in% valid_sg_focus) {
    stop("'sg_focus' must be one of: ", paste(valid_sg_focus, collapse = ", "))
  }

  # Check data quality in df
  if (nrow(df) == 0) {
    stop("'df' has no rows")
  }

  if (any(is.na(df$Y)) || any(is.na(df$Event)) || any(is.na(df$Treat))) {
    stop("df contains missing values in Y, Event, or Treat columns")
  }

  if (!all(df$Event %in% c(0, 1))) {
    stop("Event column must contain only 0 and 1")
  }

  if (!all(df$Treat %in% c(0, 1))) {
    stop("Treat column must contain only 0 and 1")
  }

  if (any(df$Y < 0)) {
    stop("Y (outcome) must be non-negative")
  }

  # Check sufficient data
  min_events <- 10  # Reasonable minimum for Cox models
  if (sum(df$Event) < min_events) {
    stop("Insufficient events (< ", min_events, ") for analysis")
  }

  events_by_treat <- table(df$Event, df$Treat)
  if (any(events_by_treat[2, ] < 5)) {
    warning("Very few events in one or both treatment arms. Results may be unstable.")
  }

  # =========================================================================
  # SECTION 2: DEFINE get_split_hr HELPER FUNCTION
  # =========================================================================

  get_split_hr <- function(df, cox_initial = NULL) {
    if (nrow(df) < 2 || sum(df$Event) < 2) {
      return(NA_real_)
    }

    hr <- try({
      fit <- survival::coxph(survival::Surv(Y, Event) ~ Treat,
                             data = df,
                             init = cox_initial,
                             robust = FALSE)
      summary(fit)$conf.int[1, 1]
    }, silent = TRUE)

    if (inherits(hr, "try-error")) return(NA_real_)
    return(hr)
  }

  # =========================================================================
  # SECTION 3: EXTRACT AND VALIDATE COVARIATE NAMES
  # =========================================================================

  # Extract Z names (factor columns)
  exclude_cols <- c("grp", "K", "n", "E", "d1", "m1", "m0", "HR", "L(HR)", "U(HR)")
  names.Z <- setdiff(names(hr.subgroups), exclude_cols)

  if (length(names.Z) == 0) {
    stop("No covariate columns found in hr.subgroups after excluding metric columns")
  }

  # Validate Lsg matches actual number of covariates
  if (length(names.Z) != Lsg) {
    stop("Lsg (", Lsg, ") does not match actual number of covariates (",
         length(names.Z), ") in hr.subgroups")
  }

  # Validate that covariate columns contain only 0/1
  z_cols <- hr.subgroups[, names.Z, with = FALSE]
  for (col_name in names.Z) {
    unique_vals <- unique(z_cols[[col_name]])
    if (!all(unique_vals %in% c(0, 1, NA))) {
      stop("Covariate column '", col_name, "' contains values other than 0, 1, or NA")
    }
  }

  # =========================================================================
  # SECTION 4: FILTER SUBGROUPS BY CRITERIA
  # =========================================================================

  # Check if any subgroups to evaluate
  if (nrow(hr.subgroups) == 0) {
    warning("No subgroups provided in hr.subgroups")
    return(list(
      out_hr = NULL, out_maxSG = NULL, out_minSG = NULL,
      df_flag = NULL, sg.harm = NULL, sg.harm.id = NULL
    ))
  }

  # Filter by m1 if finite threshold
  if (is.finite(m1.threshold)) {
    # Check if m1 column has values
    if (all(is.na(hr.subgroups$m1))) {
      warning("m1.threshold specified but all m1 values are NA")
    } else {
      # Remove rows where m1 is NA before filtering
      hr.subgroups <- hr.subgroups[!is.na(hr.subgroups$m1), ]
      if (nrow(hr.subgroups) == 0) {
        warning("All subgroups removed after filtering NA m1 values")
        return(list(
          out_hr = NULL, out_maxSG = NULL, out_minSG = NULL,
          df_flag = NULL, sg.harm = NULL, sg.harm.id = NULL
        ))
      }
    }
  }

  # Apply HR and m1 filters
  if (is.finite(m1.threshold)) {
    found.hrs <- hr.subgroups[hr.subgroups$HR >= hr.threshold &
                                hr.subgroups$m1 <= m1.threshold, ]
  } else {
    found.hrs <- hr.subgroups[hr.subgroups$HR >= hr.threshold, ]
  }

  # Check if any subgroups pass filters
  if (nrow(found.hrs) == 0) {
    if (details) {
      cat("No subgroups meet criteria (HR >=", hr.threshold)
      if (is.finite(m1.threshold)) cat(" and m1 <=", m1.threshold)
      cat(")\n")
    }
    return(list(
      out_hr = NULL, out_maxSG = NULL, out_minSG = NULL,
      df_flag = NULL, sg.harm = NULL, sg.harm.id = NULL
    ))
  }

  if (details) {
    cat("# of subgroups meeting initial criteria:", nrow(found.hrs), "\n")
  }

  # =========================================================================
  # SECTION 5: REMOVE NEAR-DUPLICATE SUBGROUPS
  # =========================================================================

  if (nrow(found.hrs) > 1) {
    n_before <- nrow(found.hrs)

    tryCatch({
      found.hrs <- remove_near_duplicate_subgroups(found.hrs, details = details)
    }, error = function(e) {
      warning("Error removing duplicates: ", e$message, ". Proceeding with original subgroups.")
    })

    if (nrow(found.hrs) == 0) {
      stop("All subgroups removed during duplicate removal. This should not happen.")
    }

    if (details && nrow(found.hrs) < n_before) {
      cat("Removed", n_before - nrow(found.hrs), "near-duplicate subgroups\n")
    }
  }

  # =========================================================================
  # SECTION 6: SORT AND LIMIT SUBGROUPS
  # =========================================================================

  # Sort based on sg_focus
  if (sg_focus == "maxSG") {
    found.hrs <- found.hrs[order(found.hrs$n, decreasing = TRUE), ]
  } else if (sg_focus == "minSG") {
    found.hrs <- found.hrs[order(found.hrs$n, decreasing = FALSE), ]
  }
  # For "hr", "hrMaxSG", "hrMinSG": sorting happens later in sg_consistency_out

  # Extract index matrix
  index.Z <- found.hrs[, names.Z, with = FALSE]

  # Validate index.Z dimensions
  if (ncol(index.Z) != Lsg) {
    stop("index.Z has ", ncol(index.Z), " columns but Lsg=", Lsg)
  }

  if (details) {
    cat("# of unique initial candidates:", nrow(found.hrs), "\n")
  }

  # Limit to top stop_Kgroups
  maxsgs <- min(nrow(found.hrs), stop_Kgroups)
  found.hrs <- found.hrs[seq_len(maxsgs), ]

  if (details) {
    cat("# Restricting to top stop_Kgroups =", stop_Kgroups, "\n")
    cat("# of candidates restricted to 'top K':", nrow(found.hrs), "\n")
  }

  # Final check
  if (nrow(found.hrs) == 0) {
    stop("No subgroups remaining after restriction. This should not happen.")
  }

  # =========================================================================
  # SECTION 7: VALIDATE PARALLEL CONFIGURATION
  # =========================================================================

  use_parallel <- length(parallel_args) > 0 && !is.null(parallel_args[[1]])

  if (use_parallel) {
    # Validate parallel_args structure
    required_parallel <- c("plan", "workers")
    if (!all(required_parallel %in% names(parallel_args))) {
      warning("parallel_args missing required elements. Using sequential processing.")
      use_parallel <- FALSE
    }

    # Check plan validity
    valid_plans <- c("multisession", "multicore", "callr", "sequential")
    if (!parallel_args$plan %in% valid_plans) {
      warning("Invalid parallel plan '", parallel_args$plan,
              "'. Using sequential processing.")
      use_parallel <- FALSE
    }

    # Check workers
    if (!is.numeric(parallel_args$workers) || parallel_args$workers < 1) {
      warning("Invalid workers value. Using sequential processing.")
      use_parallel <- FALSE
    }
  }

  if (details) {
    if (use_parallel) {
      cat("Using parallel processing:", parallel_args$plan,
          "with", parallel_args$workers, "workers\n")
    } else {
      cat("Using sequential processing\n")
    }
  }

  # =========================================================================
  # SECTION 8: INITIALIZE RESULTS STORAGE
  # =========================================================================

  res <- NULL
  any.found <- 0
  resultk <- NULL
  if (details) t.start <- proc.time()[3]

  # =========================================================================
  # SECTION 9: EVALUATE EACH SUBGROUP (SEQUENTIAL OR PARALLEL)
  # =========================================================================

  # Standard lapply (sequential)
  if (!use_parallel) {
    results_list <- lapply(seq_len(nrow(found.hrs)), function(m) {
      # --- Begin: code for each subgroup ---

      # Validate subgroup extraction
      if (m < 1 || m > nrow(index.Z)) {
        warning("Invalid subgroup index m=", m, ". Skipping.")
        return(NULL)
      }

      indexm <- as.numeric(unlist(index.Z[m, ]))

      if (length(indexm) != length(names.Z)) {
        warning("Subgroup ", m, ": index length mismatch. Skipping.")
        return(NULL)
      }

      if (!all(indexm %in% c(0, 1, NA))) {
        warning("Subgroup ", m, ": invalid index values. Skipping.")
        return(NULL)
      }

      this.m <- names.Z[indexm == 1]

      if (length(this.m) == 0) {
        warning("Subgroup ", m, ": no factors selected. Skipping.")
        return(NULL)
      }

      # if (length(this.m) > maxk) {
      #   warning("Subgroup ", m, ": ", length(this.m), " factors exceeds maxk=", maxk, ". Skipping.")
      #   return(NULL)
      # }

      # Validate label conversion
      this.m_label <- tryCatch({
        unlist(lapply(this.m, FS_labels, confs_labels = confs_labels))
      }, error = function(e) {
        warning("Subgroup ", m, ": error in FS_labels: ", e$message, ". Skipping.")
        return(NULL)
      })

      if (is.null(this.m_label)) return(NULL)

      if (length(this.m_label) != length(this.m)) {
        warning("Subgroup ", m, ": label length mismatch. Skipping.")
        return(NULL)
      }

      # Create subgroup definition
      id.m <- paste(paste(this.m, collapse = "==1 & "), "==1")

      # Extract subgroup data
      df.sub <- tryCatch({
        subset(df, eval(parse(text = id.m)))
      }, error = function(e) {
        warning("Subgroup ", m, ": error extracting data: ", e$message, ". Skipping.")
        return(NULL)
      })

      if (is.null(df.sub)) return(NULL)

      if (nrow(df.sub) == 0) {
        warning("Subgroup ", m, ": no observations match criteria. Skipping.")
        return(NULL)
      }

      df.x <- data.table::data.table(df.sub)
      N.x <- nrow(df.x)

      # Get cox_init safely
      cox_init <- log(found.hrs$HR[m])
      if (is.na(cox_init) || is.infinite(cox_init)) {
        cox_init <- 0  # Fallback to null effect
        if (details) {
          cat("Subgroup ", m, ": using default cox_init=0\n")
        }
      }

      # Perform splits
      flag.consistency <- sapply(seq_len(n.splits), function(bb) {
        # Create splits with error handling
        in.split1 <- tryCatch({
          sample(c(TRUE, FALSE), N.x, replace = TRUE, prob = c(0.5, 0.5))
        }, error = function(e) {
          warning("Split ", bb, " of subgroup ", m, ": sampling error. Returning NA.")
          return(NULL)
        })

        if (is.null(in.split1)) return(NA_real_)

        df.x$insplit1 <- in.split1

        df.x.split1 <- subset(df.x, insplit1 == 1)
        df.x.split2 <- subset(df.x, insplit1 == 0)

        # Check split sizes
        if (nrow(df.x.split1) < 5 || nrow(df.x.split2) < 5) {
          return(NA_real_)  # Skip this split
        }

        # Check events in splits
        if (sum(df.x.split1$Event) < 2 || sum(df.x.split2$Event) < 2) {
          return(NA_real_)  # Need at least 2 events per split
        }

        # Fit models with error suppression
        hr.split1 <- suppressWarnings(get_split_hr(df = df.x.split1, cox_initial = cox_init))
        hr.split2 <- suppressWarnings(get_split_hr(df = df.x.split2, cox_initial = cox_init))

        # Return consistency flag
        if (!is.na(hr.split1) && !is.na(hr.split2)) {
          as.numeric(hr.split1 > hr.consistency && hr.split2 > hr.consistency)
        } else {
          NA_real_
        }
      })

      # Check if we got any valid splits
      n_valid_splits <- sum(!is.na(flag.consistency))

      if (n_valid_splits == 0) {
        # NO valid splits - skip this subgroup entirely
        if (details) {
          cat("Bootstrap", boot, ": No valid consistency splits for subgroup\n")
        }
        # Don't calculate p.consistency, just skip
        p.consistency <- NA

      } else {
      # At least some valid splits
      # Calculate consistency with error handling
      p.consistency <- tryCatch({
        round(mean(flag.consistency, na.rm = TRUE), pconsistency.digits)
      }, error = function(e) {
        warning("Subgroup ", m, ": error calculating consistency: ", e$message)
        return(NA_real_)
      })
      if (is.na(p.consistency)) {
        warning("Subgroup ", m, ": consistency calculation failed. Skipping.")
        return(NULL)
      }

      if (isTRUE(p.consistency < pconsistency.threshold)) {
        if(details) cat("*** Not met: Subgroup, % Consistency =", c(this.m_label, p.consistency), "\n")
      } else {
        k <- length(this.m)
        covsm <- rep(m, maxk)
        mindex <- c(1:maxk)
        Mnames <- paste(covsm, mindex, sep = ".")
        mfound <- matrix(rep("", maxk))
        mfound[c(1:k)] <- this.m_label
        resultk <- c(p.consistency, found.hrs$HR[m], found.hrs$n[m], found.hrs$E[m], found.hrs$grp[m], m, k, mfound)
        names(resultk) <- c("Pcons", "hr", "N", "E", "g", "m", "K", Mnames)
        if (details) {
          cat("Consistency met!\n")
          cat("# of splits =", c(n.splits), "\n")
          cat("**** Subgroup, % Consistency Met=", c(this.m_label, p.consistency), "\n")
        }
        return(resultk)
      }
      }
      # --- End: code for each subgroup ---
      return(NULL)
    })
  } else {
    # Parallel processing
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)  # Restore plan on exit
    setup_parallel_SGcons(parallel_args)

    results_list <- future.apply::future_lapply(seq_len(nrow(found.hrs)), function(m) {
      # --- Begin: code for each subgroup (same as above) ---

      # Validate subgroup extraction
      if (m < 1 || m > nrow(index.Z)) {
        warning("Invalid subgroup index m=", m, ". Skipping.")
        return(NULL)
      }

      indexm <- as.numeric(unlist(index.Z[m, ]))

      if (length(indexm) != length(names.Z)) {
        warning("Subgroup ", m, ": index length mismatch. Skipping.")
        return(NULL)
      }

      if (!all(indexm %in% c(0, 1, NA))) {
        warning("Subgroup ", m, ": invalid index values. Skipping.")
        return(NULL)
      }

      this.m <- names.Z[indexm == 1]

      if (length(this.m) == 0) {
        warning("Subgroup ", m, ": no factors selected. Skipping.")
        return(NULL)
      }

      # Validate label conversion
      this.m_label <- tryCatch({
        unlist(lapply(this.m, FS_labels, confs_labels = confs_labels))
      }, error = function(e) {
        warning("Subgroup ", m, ": error in FS_labels: ", e$message, ". Skipping.")
        return(NULL)
      })

      if (is.null(this.m_label)) return(NULL)

      if (length(this.m_label) != length(this.m)) {
        warning("Subgroup ", m, ": label length mismatch. Skipping.")
        return(NULL)
      }

      # Create subgroup definition
      id.m <- paste(paste(this.m, collapse = "==1 & "), "==1")

      # Extract subgroup data
      df.sub <- tryCatch({
        subset(df, eval(parse(text = id.m)))
      }, error = function(e) {
        warning("Subgroup ", m, ": error extracting data: ", e$message, ". Skipping.")
        return(NULL)
      })

      if (is.null(df.sub)) return(NULL)

      if (nrow(df.sub) == 0) {
        warning("Subgroup ", m, ": no observations match criteria. Skipping.")
        return(NULL)
      }

      df.x <- data.table::data.table(df.sub)
      N.x <- nrow(df.x)

      # Get cox_init safely
      cox_init <- log(found.hrs$HR[m])
      if (is.na(cox_init) || is.infinite(cox_init)) {
        cox_init <- 0
        if (details) {
          cat("Subgroup ", m, ": using default cox_init=0\n")
        }
      }

      # Perform splits
      flag.consistency <- sapply(seq_len(n.splits), function(bb) {
        in.split1 <- tryCatch({
          sample(c(TRUE, FALSE), N.x, replace = TRUE, prob = c(0.5, 0.5))
        }, error = function(e) {
          warning("Split ", bb, " of subgroup ", m, ": sampling error. Returning NA.")
          return(NULL)
        })

        if (is.null(in.split1)) return(NA_real_)

        df.x$insplit1 <- in.split1

        df.x.split1 <- subset(df.x, insplit1 == 1)
        df.x.split2 <- subset(df.x, insplit1 == 0)

        if (nrow(df.x.split1) < 5 || nrow(df.x.split2) < 5) {
          return(NA_real_)
        }

        if (sum(df.x.split1$Event) < 2 || sum(df.x.split2$Event) < 2) {
          return(NA_real_)
        }

        hr.split1 <- suppressWarnings(get_split_hr(df = df.x.split1, cox_initial = cox_init))
        hr.split2 <- suppressWarnings(get_split_hr(df = df.x.split2, cox_initial = cox_init))

        if (!is.na(hr.split1) && !is.na(hr.split2)) {
          as.numeric(hr.split1 > hr.consistency && hr.split2 > hr.consistency)
        } else {
          NA_real_
        }
      })

      n_valid_splits <- sum(!is.na(flag.consistency))

      if (n_valid_splits < 10) {
        warning("Subgroup ", m, ": only ", n_valid_splits, " valid splits out of ",
                n.splits, ". Results may be unreliable.")
      }

      if (n_valid_splits == 0) {
        # NO valid splits - skip this subgroup entirely
        if (details) {
          cat("Bootstrap", boot, ": No valid consistency splits for subgroup\n")
        }
        # Don't calculate p.consistency, just skip
        p.consistency <- NA

      } else {
        # At least some valid splits
        # Calculate consistency with error handling
        p.consistency <- tryCatch({
          round(mean(flag.consistency, na.rm = TRUE), pconsistency.digits)
        }, error = function(e) {
          warning("Subgroup ", m, ": error calculating consistency: ", e$message)
          return(NA_real_)
        })
        if (is.na(p.consistency)) {
          warning("Subgroup ", m, ": consistency calculation failed. Skipping.")
          return(NULL)
        }

        if (isTRUE(p.consistency < pconsistency.threshold)) {
          if(details) cat("*** Not met: Subgroup, % Consistency =", c(this.m_label, p.consistency), "\n")
        } else {
        k <- length(this.m)
        covsm <- rep("M", maxk)
        mindex <- c(1:maxk)
        Mnames <- paste(covsm, mindex, sep = ".")
        mfound <- matrix(rep("", maxk))
        mfound[c(1:k)] <- this.m_label
        resultk <- c(p.consistency, found.hrs$HR[m], found.hrs$n[m], found.hrs$E[m], found.hrs$grp[m], m, k, mfound)
        names(resultk) <- c("Pcons", "hr", "N", "E", "g", "m", "K", Mnames)
        if (details) {
          cat("Consistency met!\n")
          cat("# of splits =", c(n.splits), "\n")
          cat("**** Subgroup, % Consistency Met=", c(this.m_label, p.consistency), "\n")
        }
        return(resultk)
        }
      }
      # --- End: code for each subgroup ---
      return(NULL)
    },
    future.seed = TRUE,
    future.packages = c("survival", "data.table"),
    future.globals = structure(
      TRUE,
      add = c("get_split_hr", "FS_labels", "hr.consistency", "pconsistency.threshold")
    )
    )
  }

  # =========================================================================
  # SECTION 10: COMPILE RESULTS
  # =========================================================================

  # Check if any results returned
  if (length(results_list) == 0) {
    if (details) {
      cat("No subgroups met consistency criteria\n")
    }
    return(list(
      out_hr = NULL, out_maxSG = NULL, out_minSG = NULL,
      df_flag = NULL, sg.harm = NULL, sg.harm.id = NULL
    ))
  }

  # Filter NULL results
  results_list <- Filter(Negate(is.null), results_list)

  if (length(results_list) == 0) {
    if (details) {
      cat("All subgroup evaluations returned NULL\n")
    }
    return(list(
      out_hr = NULL, out_maxSG = NULL, out_minSG = NULL,
      df_flag = NULL, sg.harm = NULL, sg.harm.id = NULL
    ))
  }

  # Combine results
  res <- tryCatch({
    data.table::as.data.table(do.call(rbind, results_list))
  }, error = function(e) {
    stop("Error combining results: ", e$message,
         "\nThis may indicate inconsistent result structure across subgroups.")
  })

  any.found <- nrow(res)

  if (any.found == 0) {
    if (details) {
      cat("No subgroups found meeting consistency threshold\n")
    }
    return(list(
      out_hr = NULL, out_maxSG = NULL, out_minSG = NULL,
      df_flag = NULL, sg.harm = NULL, sg.harm.id = NULL
    ))
  }

  # Convert columns to numeric
  cols_to_numeric <- c("Pcons", "hr", "N", "E", "K")

  # Check columns exist
  missing_cols <- setdiff(cols_to_numeric, names(res))
  if (length(missing_cols) > 0) {
    stop("Result data.table missing expected columns: ",
         paste(missing_cols, collapse = ", "))
  }

  tryCatch({
    res[, (cols_to_numeric) := lapply(.SD, as.numeric), .SDcols = cols_to_numeric]
  }, error = function(e) {
    stop("Error converting result columns to numeric: ", e$message)
  })

  # =========================================================================
  # SECTION 11: GENERATE OUTPUTS FOR DIFFERENT sg_focus VALUES
  # =========================================================================

  out_hr <- out_maxSG <- out_minSG <- NULL
  df_flag <- sg.harm <- sg.harm.id <- NULL

  if (any.found > 0) {
    result_new <- data.table::copy(res)

    # Generate output for 'hr' focus
    sgdetails <- ifelse(plot.sg && sg_focus == "hr", TRUE, FALSE)
    out_hr <- tryCatch({
      sg_consistency_out(df = df, result_new = result_new, sg_focus = "hr", details = sgdetails,
                         plot.sg = sgdetails, index.Z = index.Z, names.Z = names.Z, by.risk = by.risk, confs_labels = confs_labels)
    }, error = function(e) {
      warning("Error in sg_consistency_out for 'hr': ", e$message)
      NULL
    })

    # Generate output for 'maxSG' focus
    sgdetails <- ifelse(plot.sg && sg_focus %in% c("hrMaxSG", "maxSG"), TRUE, FALSE)
    out_maxSG <- tryCatch({
      sg_consistency_out(df = df, result_new = result_new, sg_focus = "maxSG", details = sgdetails,
                         plot.sg = sgdetails, index.Z = index.Z, names.Z = names.Z, by.risk = by.risk, confs_labels = confs_labels)
    }, error = function(e) {
      warning("Error in sg_consistency_out for 'maxSG': ", e$message)
      NULL
    })

    # Generate output for 'minSG' focus
    sgdetails <- ifelse(plot.sg && sg_focus %in% c("hrMinSG", "minSG"), TRUE, FALSE)
    out_minSG <- tryCatch({
      sg_consistency_out(df = df, result_new = result_new, sg_focus = "minSG", details = sgdetails,
                         plot.sg = sgdetails, index.Z = index.Z, names.Z = names.Z, by.risk = by.risk, confs_labels = confs_labels)
    }, error = function(e) {
      warning("Error in sg_consistency_out for 'minSG': ", e$message)
      NULL
    })

    # Create a mapping from sg_focus to the corresponding object
    sg_map <- list(
      hr      = out_hr,
      hrMaxSG = out_maxSG,
      maxSG   = out_maxSG,
      hrMinSG = out_minSG,
      minSG   = out_minSG
    )

    # Check if sg_focus is valid
    if (!sg_focus %in% names(sg_map)) {
      stop(sprintf("Unknown sg_focus value: %s", sg_focus))
    }

    # Extract the relevant object
    sg_obj <- sg_map[[sg_focus]]

    # Check if selected output exists
    if (is.null(sg_obj)) {
      warning("No valid output for sg_focus='", sg_focus, "'")
      df_flag <- NULL
      sg.harm <- NULL
      sg.harm.id <- NULL
    } else {
      # Validate sg_obj structure
      required_fields <- c("df_flag", "sg.harm_label", "sg.harm.id")
      missing_fields <- setdiff(required_fields, names(sg_obj))
      if (length(missing_fields) > 0) {
        stop("sg_consistency_out result missing fields: ",
             paste(missing_fields, collapse = ", "))
      }

      # Assign variables
      df_flag    <- sg_obj$df_flag
      sg.harm    <- sg_obj$sg.harm_label
      sg.harm.id <- sg_obj$sg.harm.id
    }

    if (details) cat("SG focus=", c(sg_focus), "\n")
  }

  # =========================================================================
  # SECTION 12: VALIDATE AND RETURN OUTPUT
  # =========================================================================

  if (details) {
    t.end <- proc.time()[3]
    t.min <- (t.end - t.start) / 60
    cat("Subgroup Consistency Minutes=", c(t.min), "\n")
    if (any.found > 0) {
      cat("Subgroup found (FS)\n")
      if (!is.null(sg.harm)) {
        cat("Selected subgroup:", paste(sg.harm, collapse = " & "), "\n")
      }
    } else {
      cat("NO subgroup found (FS)\n")
    }
  }

  # Build output list
  output <- list(
    out_hr = out_hr,
    out_maxSG = out_maxSG,
    out_minSG = out_minSG,
    df_flag = df_flag,
    sg.harm = sg.harm,
    sg.harm.id = sg.harm.id
  )

  # Validate output structure
  if (!is.null(df_flag)) {
    if (!is.data.frame(df_flag) && !data.table::is.data.table(df_flag)) {
      warning("df_flag is not a data.frame or data.table")
    }

    required_flag_cols <- c("id", "treat.recommend")
    missing_flag_cols <- setdiff(required_flag_cols, names(df_flag))
    if (length(missing_flag_cols) > 0) {
      warning("df_flag missing expected columns: ",
              paste(missing_flag_cols, collapse = ", "))
    }
  }

  if (!is.null(sg.harm)) {
    if (!is.character(sg.harm)) {
      warning("sg.harm should be a character vector")
    }
  }

  return(output)
}
