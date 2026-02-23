# =============================================================================
# Subgroup Consistency Helper Functions
# =============================================================================
#
# Helper functions for subgroup consistency evaluation that need to be
# exported for parallel worker access. When using callr or multisession
# parallelization, workers only see exported functions from installed packages.
#
# =============================================================================

# =============================================================================
# PARALLEL SETUP
# =============================================================================

#' Set up parallel processing for subgroup consistency
#'
#' Sets up parallel processing using the specified approach and number of workers.
#'
#' @param parallel_args List with \code{plan} (character), \code{workers} (integer),
#'   and \code{show_message} (logical).
#'
#' @return None. Sets up parallel backend as side effect.
#'
#' @importFrom future plan multisession multicore sequential
#' @export
setup_parallel_SGcons <- function(parallel_args = list(
    plan = "multisession",
    workers = 4,
    show_message = TRUE
)) {
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
    future::plan(future::multisession, workers = n_workers)
    if (show_message) message("Parallel plan: multisession with ", n_workers, " workers.")
  } else if (plan_type == "multicore") {
    future::plan(future::multicore, workers = n_workers)
    if (show_message) message("Parallel plan: multicore with ", n_workers, " workers.")
  } else if (plan_type == "callr") {
    if (!requireNamespace("future.callr", quietly = TRUE)) {
      stop("The 'future.callr' package is required for the callr approach.")
    }
    future::plan(future.callr::callr, workers = n_workers)
    if (show_message) message("Parallel plan: callr with ", n_workers, " workers.")
  } else if (plan_type == "sequential") {
    future::plan(future::sequential)
    if (show_message) message("Sequential plan")
  } else {
    stop("Unknown parallel plan: ", plan_type)
  }
}


# =============================================================================
# STATISTICAL HELPER FUNCTIONS
# =============================================================================

#' Wilson Score Confidence Interval
#'
#' Computes Wilson score confidence interval for a proportion, which has
#' better coverage properties than the normal approximation for small samples
#' and proportions near 0 or 1.
#'
#' @param x Integer. Number of successes.
#' @param n Integer. Number of trials.
#' @param conf.level Numeric. Confidence level (default 0.95).
#'
#' @return Named numeric vector with elements: estimate, lower, upper.
#'
#' @export
wilson_ci <- function(x, n, conf.level = 0.95) {
  if (n == 0) {
    return(c(estimate = NA_real_, lower = 0, upper = 1))
  }

  z <- qnorm(1 - (1 - conf.level) / 2)
  p_hat <- x / n

  denominator <- 1 + z^2 / n
  center <- (p_hat + z^2 / (2 * n)) / denominator
  margin <- (z / denominator) * sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2))

  lower <- max(0, center - margin)
  upper <- min(1, center + margin)

  c(estimate = p_hat, lower = lower, upper = upper)
}


#' Early Stopping Decision
#'
#' Evaluates whether enough evidence exists to stop early based on
#' confidence interval for consistency proportion.
#'
#' @param n_success Integer. Number of splits meeting consistency.
#' @param n_total Integer. Total number of valid splits.
#' @param threshold Numeric. Target consistency threshold.
#' @param conf.level Numeric. Confidence level for decision (default 0.95).
#' @param min_samples Integer. Minimum samples before allowing early stop.
#'
#' @return Character. One of "continue", "pass", or "fail".
#'
#' @export
early_stop_decision <- function(n_success, n_total, threshold,
                                conf.level = 0.95, min_samples = 20) {
  if (n_total < min_samples) {
    return("continue")
  }

  ci <- wilson_ci(n_success, n_total, conf.level)

  if (ci["lower"] >= threshold) {
    return("pass")
  }

  if (ci["upper"] < threshold) {
    return("fail")
  }

  return("continue")
}


# =============================================================================
# COX MODEL HELPERS
# =============================================================================

#' Fast Cox Model HR Estimation
#'
#' Fits a minimal Cox model to estimate hazard ratio with reduced overhead.
#' Disables robust variance, model matrix storage, and other extras for speed.
#'
#' @param df data.frame or data.table with Y, Event, Treat columns.
#' @param cox_init Numeric. Initial value for coefficient (default 0).
#'
#' @return Numeric. Estimated hazard ratio, or NA if model fails.
#'
#' @importFrom survival coxph Surv
#' @export
get_split_hr_fast <- function(df, cox_init = 0) {
  if (nrow(df) < 2 || sum(df$Event) < 2) {
    return(NA_real_)
  }

  fit <- tryCatch(
    suppressWarnings(
      survival::coxph(
        survival::Surv(Y, Event) ~ Treat,
        data = df,
        init = cox_init,
        robust = FALSE,
        model = FALSE,
        x = FALSE,
        y = FALSE
      )
    ),
    error = function(e) NULL
  )

  if (is.null(fit)) return(NA_real_)
  return(exp(fit$coefficients[1]))
}


#' Run Single Consistency Split
#'
#' Performs one random 50/50 split and evaluates whether both halves
#' meet the HR consistency threshold.
#'
#' @param df.x data.table. Subgroup data with columns Y, Event, Treat.
#' @param N.x Integer. Number of observations in subgroup.
#' @param hr.consistency Numeric. Minimum HR threshold for consistency.
#' @param cox_init Numeric. Initial value for Cox model (log HR).
#'
#' @return Numeric. 1 if both splits meet threshold, 0 if not, NA if error.
#'
#' @export
run_single_consistency_split <- function(df.x, N.x, hr.consistency, cox_init = 0) {

  in.split1 <- tryCatch({
    sample(c(TRUE, FALSE), N.x, replace = TRUE, prob = c(0.5, 0.5))
  }, error = function(e) {
    return(NULL)
  })

  if (is.null(in.split1)) return(NA_real_)

  df.x$insplit1 <- in.split1
  df.x.split1 <- df.x[df.x$insplit1 == TRUE, ]
  df.x.split2 <- df.x[df.x$insplit1 == FALSE, ]

  if (nrow(df.x.split1) < 5 || nrow(df.x.split2) < 5) {
    return(NA_real_)
  }

  if (sum(df.x.split1$Event) < 2 || sum(df.x.split2$Event) < 2) {
    return(NA_real_)
  }

  hr.split1 <- get_split_hr_fast(df.x.split1, cox_init)
  hr.split2 <- get_split_hr_fast(df.x.split2, cox_init)

  if (!is.na(hr.split1) && !is.na(hr.split2)) {
    as.numeric(hr.split1 > hr.consistency && hr.split2 > hr.consistency)
  } else {
    NA_real_
  }
}


# =============================================================================
# LABEL HELPERS
# =============================================================================

#' Convert Factor Code to Label
#'
#' Converts q-indexed codes to human-readable labels using the confs_labels
#' mapping. Supports both full format ("q1.1", "q3.0") and short format
#' ("q1", "q3"). Handles vector input via recursion.
#'
#' @param Qsg Character. Factor code in format `"q<index>.<action>"` or
#'   `"q<index>"`. For the full format, action 0 = NOT (complement),
#'   action 1 = IN (member). Short format defaults to action 1.
#'   Can also be a character vector for vectorized input.
#' @param confs_labels Character vector. Labels for each factor, indexed by
#'   factor number.
#'
#' @return Character. Human-readable label wrapped in braces, e.g.,
#'   `"{age <= 50}"` or `"!{age <= 50}"` for complement. Returns the
#'   original code if no match is found.
#'
#' @examples
#' \dontrun{
#' labels <- c("age <= 50", "tumor size <= 20", "nodes <= 3")
#' FS_labels("q1.1", labels)
#' FS_labels("q1.0", labels)
#' FS_labels("q2", labels)
#' FS_labels(c("q1", "q3"), labels)
#' }
#' @export

FS_labels <- function(Qsg, confs_labels) {
  # Handle vector input
  if (length(Qsg) > 1) {
    return(vapply(Qsg, FS_labels, character(1), confs_labels = confs_labels,
                  USE.NAMES = FALSE))
  }

  # Try full format: q<index>.<action> (e.g., "q1.1", "q3.0")
  pattern_full <- "^q(\\d+)\\.(\\d)$"
  matches <- regmatches(Qsg, regexec(pattern_full, Qsg))[[1]]

  if (length(matches) >= 3) {
    idx <- as.integer(matches[2])
    action <- matches[3]
  } else {
    # Try short format: q<index> (e.g., "q1") — default action = 1
    pattern_short <- "^q(\\d+)$"
    matches <- regmatches(Qsg, regexec(pattern_short, Qsg))[[1]]

    if (length(matches) >= 2) {
      idx <- as.integer(matches[2])
      action <- "1"
    } else {
      return(Qsg)
    }
  }

  if (idx < 1 || idx > length(confs_labels)) return(Qsg)

  base_label <- confs_labels[idx]

  if (action == "0") {
    paste0("!{", base_label, "}")
  } else {
    paste0("{", base_label, "}")
  }
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
  if (sg_focus == "maxSG") data.table::setorder(result_new, -N, -Pcons, K)
  if (sg_focus == "hrMaxSG") data.table::setorder(result_new, -N, -hr, K)
  if (sg_focus == "minSG") data.table::setorder(result_new, N, -Pcons, K)
  if (sg_focus == "hrMinSG") data.table::setorder(result_new, N, -hr, K)
  result_new
}

#' Sort Subgroups by Focus at consistency stage (consistency not available at this point)
#'
#' Sorts a data.table of subgroup results according to the specified focus.
#'
#' @param result_new A data.table of subgroup results.
#' @param sg_focus Sorting focus: "hr", "hrMaxSG", "maxSG", "hrMinSG", "minSG".
#' @return A sorted data.table.
#' @importFrom data.table setorder
#' @keywords internal
sort_subgroups_preview <- function(result_new, sg_focus) {
  if (sg_focus == "hr") data.table::setorder(result_new, -HR, K)
  if (sg_focus == "maxSG") data.table::setorder(result_new, -n, K)
  if (sg_focus == "hrMaxSG") data.table::setorder(result_new, -HR, -n, K)
  if (sg_focus == "minSG") data.table::setorder(result_new, n, K)
  if (sg_focus == "hrMinSG") data.table::setorder(result_new, -HR, n, K)
  result_new
}

#' Extract Subgroup Information
#'
#' Extracts subgroup definition and membership from results.
#'
#' @param df Data.frame. Original analysis data.
#' @param top_result Data.table row. Top subgroup result.
#' @param index.Z Matrix. Factor indicators for all subgroups.
#' @param names.Z Character vector. Factor column names.
#' @param confs_labels Character vector. Human-readable labels.
#'
#' @return List with sg.harm, sg.harm_label, df_flag, sg.harm.id.
#'
#' @keywords internal
extract_subgroup <- function(df, top_result, index.Z, names.Z, confs_labels) {
  m <- as.integer(top_result$m)
  indexm <- as.numeric(unlist(index.Z[m, ]))
  sg.harm <- names.Z[indexm == 1]
  sg.harm_label <- vapply(sg.harm, function(q) FS_labels(q, confs_labels), character(1))

  # Create membership flag
  df_temp <- data.table::as.data.table(df)
  for (varname in sg.harm) {
    df_temp <- df_temp[df_temp[[varname]] == 1, ]
  }

  sg.harm.id <- rep(0, nrow(df))
  if (nrow(df_temp) > 0 && "id" %in% names(df_temp)) {
    sg.harm.id[df$id %in% df_temp$id] <- 1
  }

  df_flag <- data.frame(
    id = df$id,
    treat.recommend = ifelse(sg.harm.id == 1, 0, 1)
  )

  list(
    sg.harm = sg.harm,
    sg.harm_label = sg.harm_label,
    df_flag = df_flag,
    sg.harm.id = sg.harm.id
  )
}


#' Plot Subgroup Survival Curves
#'
#' Plots weighted Kaplan-Meier survival curves for a specified subgroup and its complement using the \pkg{weightedsurv} package.
#'
#' @param df.sub A data frame containing data for the subgroup of interest.
#' @param df.subC A data frame containing data for the complement subgroup.
#' @param by.risk Numeric. The risk interval for plotting (passed to \code{weightedsurv::df_counting}).
#' @param confs_labels Named character vector. Covariate label mapping (not used directly in this function, but may be used for labeling).
#' @param this.1_label Character. Label for the subgroup being plotted.
#' @param top_result Data frame row. The top subgroup result row, expected to contain a \code{Pcons} column for consistency criteria.
#'
#' @importFrom weightedsurv df_counting plot_weighted_km
#' @keywords internal

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
#' @param df Data.frame. Original analysis data.
#' @param result_new Data.table. Sorted subgroup results.
#' @param sg_focus Character. Sorting focus criterion.
#' @param index.Z Matrix. Subgroup factor indicators.
#' @param names.Z Character vector. Factor column names.
#' @param details Logical. Print details.
#' @param plot.sg Logical. Plot subgroup curves.
#' @param by.risk Numeric. Risk interval for plotting.
#' @param confs_labels Character vector. Human-readable labels.
#'
#' @return List with results, subgroup definition, labels, flags, and group id.
#'
#' @importFrom data.table copy
#' @export
sg_consistency_out <- function(df, result_new, sg_focus, index.Z, names.Z,
                               details = FALSE, plot.sg = FALSE,
                               by.risk = 12, confs_labels) {

  result_new <- sort_subgroups(result_new, sg_focus)
  top_result <- result_new[1, ]
  subgroup_info <- extract_subgroup(df, top_result, index.Z, names.Z, confs_labels)

  # === ADD THIS PLOTTING SECTION ===


  if (details && plot.sg) {
    sg.harm <- subgroup_info$sg.harm

    in_harm <- Reduce(`&`, lapply(sg.harm, function(v) df[[v]] == 1L))
    df.sub  <- df[in_harm, , drop = FALSE]
    df.subC <- df[!in_harm, , drop = FALSE]

    plot_subgroup(
      df.sub = df.sub,
      df.subC = df.subC,
      by.risk = by.risk,
      confs_labels = confs_labels,
      this.1_label = subgroup_info$sg.harm_label,
      top_result = top_result
    )
  }


  # === END PLOTTING SECTION ===

  result_out <- data.table::copy(result_new)

  list(
    result = result_out,
    sg.harm = subgroup_info$sg.harm,
    sg.harm_label = subgroup_info$sg.harm_label,
    df_flag = subgroup_info$df_flag,
    sg.harm.id = subgroup_info$sg.harm.id
  )
}

# =============================================================================
# DUPLICATE REMOVAL
# =============================================================================

#' Remove Near-Duplicate Subgroups
#'
#' Removes subgroups with nearly identical statistics (HR, n, E, etc.)
#' to reduce redundancy in candidate list.
#'
#' @param hr_subgroups Data.table of subgroup results.
#' @param tolerance Numeric. Tolerance for numeric comparison (default 0.001).
#' @param details Logical. Print details about removed duplicates.
#'
#' @return Data.table with near-duplicate rows removed.
#'
#' @importFrom data.table as.data.table
#' @keywords internal
remove_near_duplicate_subgroups <- function(hr_subgroups,
                                            tolerance = 0.001,
                                            details = FALSE) {

  df <- as.data.frame(hr_subgroups)

  # Columns to check: K, n, E, d1, m1, m0, HR, L(HR), U(HR)
  cols_to_check <- 2:min(10, ncol(df))

  df_rounded <- df
  for (i in cols_to_check) {
    if (is.numeric(df_rounded[, i])) {
      df_rounded[, i] <- round(df_rounded[, i] / tolerance) * tolerance
    }
  }

  key_cols <- df_rounded[, cols_to_check, drop = FALSE]
  dup_key <- apply(key_cols, 1, function(x) paste(x, collapse = "_"))

  keep_rows <- !duplicated(dup_key)

  n_removed <- sum(!keep_rows)

  if (details && n_removed > 0) {
    cat("Removed", n_removed, "near-duplicate subgroups\n")
  }

  return(data.table::as.data.table(hr_subgroups[keep_rows, ]))
}


#' Remove Redundant Subgroups
#'
#' Removes redundant subgroups by checking for exact ties in key columns.
#'
#' @param found.hrs Data.table of found subgroups.
#'
#' @return Data.table of non-redundant subgroups.
#'
#' @importFrom data.table data.table
#' @keywords internal
remove_redundant_subgroups <- function(found.hrs) {
  found.new <- found.hrs[order(found.hrs$HR, decreasing = TRUE), ]
  f1.hrs <- found.new[1, ]
  temp <- found.new[-c(1), ]
  temp2 <- as.matrix(found.new[, c("HR", "n", "E", "K", "L(HR)", "U(HR)")])
  id_keep <- which(round(diff(temp2, lag = 1), 6) != 0)
  fkeep.hrs <- temp[id_keep, ]
  na.omit(data.table::data.table(rbind(f1.hrs, fkeep.hrs)))
}


# =============================================================================
# FIXED-SAMPLE CONSISTENCY EVALUATION
# =============================================================================

#' Evaluate Single Subgroup for Consistency (Fixed-Sample)
#'
#' Evaluates a single subgroup for consistency across random splits using
#' a fixed number of splits.
#'
#' @param m Integer. Index of the subgroup to evaluate.
#' @param index.Z Data.table or matrix. Factor indicators for all subgroups.
#' @param names.Z Character vector. Names of factor columns.
#' @param df Data.frame. Original data with Y, Event, Treat, id columns.
#' @param found.hrs Data.table. Subgroup hazard ratio results.
#' @param n.splits Integer. Number of random splits for consistency evaluation.
#' @param hr.consistency Numeric. Minimum HR threshold for consistency.
#' @param pconsistency.threshold Numeric. Minimum proportion of splits meeting
#'   consistency.
#' @param pconsistency.digits Integer. Rounding digits for consistency proportion.
#' @param maxk Integer. Maximum number of factors in a subgroup.
#' @param confs_labels Character vector. Labels for confounders.
#' @param details Logical. Print details during execution.
#'
#' @return Named numeric vector with consistency results, or NULL if criteria
#'   not met.
#'
#' @importFrom data.table data.table
#' @importFrom survival coxph Surv
#' @export
evaluate_subgroup_consistency <- function(
    m,
    index.Z,
    names.Z,
    df,
    found.hrs,
    n.splits,
    hr.consistency,
    pconsistency.threshold,
    pconsistency.digits = 2,
    maxk,
    confs_labels,
    details = FALSE
) {

  # -------------------------------------------------------------------------
  # SECTION 1: VALIDATE INPUT
  # -------------------------------------------------------------------------

  if (m < 1 || m > nrow(index.Z)) {
    warning("Invalid subgroup index m=", m, ". Skipping.")
    return(NULL)
  }

  # -------------------------------------------------------------------------
  # SECTION 2: EXTRACT SUBGROUP DEFINITION
  # -------------------------------------------------------------------------

  indexm <- as.numeric(unlist(index.Z[m, ]))

  if (length(indexm) != length(names.Z)) {
    warning("Subgroup ", m, ": index length mismatch. Skipping.")
    return(NULL)
  }

  this.m <- names.Z[indexm == 1]

  if (length(this.m) == 0) {
    warning("Subgroup ", m, ": no factors selected. Skipping.")
    return(NULL)
  }

  # Convert to labels
  this.m_label <- vapply(this.m, function(q) FS_labels(q, confs_labels), character(1))

  # -------------------------------------------------------------------------
  # SECTION 3: IDENTIFY SUBGROUP MEMBERS
  # -------------------------------------------------------------------------

  df.x <- data.table::as.data.table(df)

  for (varname in this.m) {
    if (!varname %in% names(df.x)) {
      warning("Subgroup ", m, ": variable '", varname, "' not in df. Skipping.")
      return(NULL)
    }
    df.x <- df.x[df.x[[varname]] == 1, ]
  }

  N.x <- nrow(df.x)

  if (N.x < 10) {
    if (details) cat("Subgroup ", m, ": too few observations (", N.x, ")\n", sep = "")
    return(NULL)
  }

  # -------------------------------------------------------------------------
  # SECTION 4: INITIALIZE COX MODEL
  # -------------------------------------------------------------------------

  cox_init <- tryCatch({
    fit0 <- survival::coxph(
      survival::Surv(Y, Event) ~ Treat,
      data = df.x,
      model = FALSE,
      x = FALSE,
      y = FALSE
    )
    fit0$coefficients[1]
  }, error = function(e) 0)

  # -------------------------------------------------------------------------
  # SECTION 5: RUN CONSISTENCY SPLITS
  # -------------------------------------------------------------------------

  flag.consistency <- numeric(n.splits)

  for (i in seq_len(n.splits)) {
    flag.consistency[i] <- run_single_consistency_split(df.x, N.x, hr.consistency, cox_init)
  }

  # -------------------------------------------------------------------------
  # SECTION 6: CALCULATE CONSISTENCY PROPORTION
  # -------------------------------------------------------------------------

  n_valid <- sum(!is.na(flag.consistency))

  if (n_valid < 10) {
    if (details) cat("Subgroup ", m, ": too few valid splits (", n_valid, ")\n", sep = "")
    return(NULL)
  }

  p.consistency <- tryCatch({
    round(mean(flag.consistency, na.rm = TRUE), pconsistency.digits)
  }, error = function(e) {
    return(NA_real_)
  })

  if (is.na(p.consistency)) {
    return(NULL)
  }

  # -------------------------------------------------------------------------
  # SECTION 7: CHECK CONSISTENCY THRESHOLD
  # -------------------------------------------------------------------------

  if (isTRUE(p.consistency < pconsistency.threshold)) {
    if (details) {
      cat("*** Not met: Subgroup, % Consistency =",
          c(this.m_label, p.consistency), "\n")
    }
    return(NULL)
  }

  # -------------------------------------------------------------------------
  # SECTION 8: FORMAT AND RETURN RESULT
  # -------------------------------------------------------------------------

  k <- length(this.m)
  covsm <- rep("M", maxk)
  mindex <- seq_len(maxk)
  Mnames <- paste(covsm, mindex, sep = ".")

  mfound <- matrix(rep("", maxk))
  mfound[seq_len(k)] <- this.m_label

  resultk <- c(
    p.consistency,
    found.hrs$HR[m],
    found.hrs$n[m],
    found.hrs$E[m],
    found.hrs$grp[m],
    m,
    k,
    mfound
  )

  names(resultk) <- c("Pcons", "hr", "N", "E", "g", "m", "K", Mnames)

  if (details) {
    cat("*** Met: Subgroup, % Consistency =",
        c(this.m_label, p.consistency), "\n")
  }

  return(resultk)
}


# =============================================================================
# TWO-STAGE CONSISTENCY EVALUATION
# =============================================================================

#' Evaluate Consistency (Two-Stage Algorithm)
#'
#' Evaluates a single subgroup for consistency using a two-stage approach:
#' Stage 1 screens with fewer splits, Stage 2 uses sequential batched
#' evaluation with early stopping for efficient evaluation.
#'
#' @param m Integer. Index of subgroup to evaluate.
#' @param index.Z data.table or matrix. Factor indicators for all subgroups.
#' @param names.Z Character vector. Names of factor columns.
#' @param df data.frame. Original data with Y, Event, Treat, id columns.
#' @param found.hrs data.table. Subgroup hazard ratio results.
#' @param hr.consistency Numeric. Minimum HR threshold for consistency.
#' @param pconsistency.threshold Numeric. Final consistency threshold.
#' @param pconsistency.digits Integer. Rounding digits for output.
#' @param maxk Integer. Maximum number of factors in a subgroup.
#' @param confs_labels Character vector. Labels for confounders.
#' @param details Logical. Print progress details.
#' @param n.splits.screen Integer. Number of splits for Stage 1 (default 30).
#' @param screen.threshold Numeric. Screening threshold for Stage 1
#'   (default auto-calculated).
#' @param n.splits.max Integer. Maximum total splits (default 400).
#' @param batch.size Integer. Splits per batch in Stage 2 (default 20).
#' @param conf.level Numeric. Confidence level for early stopping (default 0.95).
#' @param min.valid.screen Integer. Minimum valid splits in Stage 1 (default 10).
#'
#' @return Named numeric vector with consistency results, or NULL if not met.
#'
#' @importFrom data.table data.table as.data.table
#' @importFrom survival coxph Surv
#' @export
evaluate_consistency_twostage <- function(
    m,
    index.Z,
    names.Z,
    df,
    found.hrs,
    hr.consistency,
    pconsistency.threshold,
    pconsistency.digits = 2,
    maxk,
    confs_labels,
    details = FALSE,
    n.splits.screen = 30,
    screen.threshold = NULL,
    n.splits.max = 400,
    batch.size = 20,
    conf.level = 0.95,
    min.valid.screen = 10
) {

  # ===========================================================================
  # EMBEDDED HELPER FUNCTIONS (for parallel execution compatibility)
  # ===========================================================================

  .wilson_ci <- function(x, n, conf.level = 0.95) {
    if (n == 0) {
      return(c(estimate = NA_real_, lower = NA_real_, upper = NA_real_))
    }
    z <- qnorm(1 - (1 - conf.level) / 2)
    p_hat <- x / n
    denom <- 1 + z^2 / n
    center <- (p_hat + z^2 / (2 * n)) / denom
    margin <- (z / denom) * sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2))
    c(estimate = p_hat, lower = max(0, center - margin), upper = min(1, center + margin))
  }

  .early_stop_decision <- function(n_success, n_total, threshold, conf.level = 0.95, min_samples = 20) {
    if (n_total < min_samples) return("continue")
    ci <- .wilson_ci(n_success, n_total, conf.level)
    if (is.na(ci["lower"]) || is.na(ci["upper"])) return("continue")
    if (ci["lower"] >= threshold) return("pass")
    if (ci["upper"] < threshold) return("fail")
    return("continue")
  }

  .get_split_hr_fast <- function(df_split, cox_initial = NULL) {
    if (nrow(df_split) < 2 || sum(df_split$Event) < 2) return(NA_real_)
    fit <- tryCatch(
      suppressWarnings(
        survival::coxph(
          survival::Surv(Y, Event) ~ Treat,
          data = df_split,
          init = cox_initial,
          robust = FALSE,
          model = FALSE,
          x = FALSE,
          y = FALSE
        )
      ),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NA_real_)
    return(exp(fit$coefficients[1]))
  }

  .run_single_consistency_split <- function(df.x, N.x, hr.cons, cox_init) {
    in.split1 <- tryCatch({
      sample(c(TRUE, FALSE), N.x, replace = TRUE, prob = c(0.5, 0.5))
    }, error = function(e) NULL)
    if (is.null(in.split1)) return(NA_real_)

    df.x$insplit1 <- in.split1
    df.x.split1 <- df.x[df.x$insplit1 == TRUE, ]
    df.x.split2 <- df.x[df.x$insplit1 == FALSE, ]

    if (nrow(df.x.split1) < 5 || nrow(df.x.split2) < 5) return(NA_real_)
    if (sum(df.x.split1$Event) < 2 || sum(df.x.split2$Event) < 2) return(NA_real_)

    hr.split1 <- .get_split_hr_fast(df_split = df.x.split1, cox_initial = cox_init)
    hr.split2 <- .get_split_hr_fast(df_split = df.x.split2, cox_initial = cox_init)

    if (!is.na(hr.split1) && !is.na(hr.split2)) {
      as.numeric(hr.split1 > hr.cons && hr.split2 > hr.cons)
    } else {
      NA_real_
    }
  }


  # ---------------------------------------------------------------------------
  # Parameter initialization
  # ---------------------------------------------------------------------------

  if (is.null(screen.threshold)) {
    se_estimate <- sqrt(pconsistency.threshold * (1 - pconsistency.threshold) / n.splits.screen)
    screen.threshold <- max(0.5, pconsistency.threshold - 2.5 * se_estimate)
  }

  # ---------------------------------------------------------------------------
  # Validate and extract subgroup
  # ---------------------------------------------------------------------------

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

  this.m_label <- vapply(this.m, function(q) FS_labels(q, confs_labels), character(1))

  # ---------------------------------------------------------------------------
  # Identify subgroup members
  # ---------------------------------------------------------------------------

  df.x <- data.table::as.data.table(df)

  for (varname in this.m) {
    if (!varname %in% names(df.x)) {
      warning("Subgroup ", m, ": variable '", varname, "' not in df. Skipping.")
      return(NULL)
    }
    df.x <- df.x[df.x[[varname]] == 1, ]
  }

  N.x <- nrow(df.x)

  if (N.x < 10) {
    if (details) cat("Subgroup ", m, ": too few observations (", N.x, ")\n", sep = "")
    return(NULL)
  }

  # ---------------------------------------------------------------------------
  # Initialize Cox model
  # ---------------------------------------------------------------------------

  cox_init <- tryCatch({
    fit0 <- survival::coxph(
      survival::Surv(Y, Event) ~ Treat,
      data = df.x,
      model = FALSE,
      x = FALSE,
      y = FALSE
    )
    fit0$coefficients[1]
  }, error = function(e) 0)

  # ---------------------------------------------------------------------------
  # Stage 1: Quick screening
  # ---------------------------------------------------------------------------

  stage1_flags <- numeric(n.splits.screen)
  for (i in seq_len(n.splits.screen)) {
    stage1_flags[i] <- .run_single_consistency_split(df.x, N.x, hr.consistency, cox_init)
  }

  n_valid_stage1 <- sum(!is.na(stage1_flags))
  n_success_stage1 <- sum(stage1_flags == 1, na.rm = TRUE)

  if (n_valid_stage1 < min.valid.screen) {
    if (details) {
      cat("Subgroup ", m, ": insufficient valid Stage 1 splits (", n_valid_stage1, ")\n", sep = "")
    }
    return(NULL)
  }

  p_stage1 <- n_success_stage1 / n_valid_stage1

  if (p_stage1 < screen.threshold) {
    if (details) {
      cat("Subgroup ", m, ": SCREENED OUT (p=", round(p_stage1, 3), ")\n", sep = "")
    }
    return(NULL)
  }

  # ---------------------------------------------------------------------------
  # Stage 2: Sequential evaluation with early stopping
  # ---------------------------------------------------------------------------

  if (details) {
    cat("Subgroup ", m, ": Stage 2 sequential\n", sep = "")
  }

  all_flags <- stage1_flags
  n_total_valid <- n_valid_stage1
  n_total_success <- n_success_stage1

  n_remaining <- n.splits.max - n.splits.screen
  n_batches <- ceiling(n_remaining / batch.size)

  final_decision <- "continue"

  for (batch_num in seq_len(n_batches)) {

    decision <- .early_stop_decision(
      n_success = n_total_success,
      n_total = n_total_valid,
      threshold = pconsistency.threshold,
      conf.level = conf.level,
      min_samples = max(20, min.valid.screen)
    )

    if (decision != "continue") {
      final_decision <- decision
      if (details) {
        cat("Subgroup ", m, ": EARLY ", toupper(decision),
            " at n=", n_total_valid, "\n", sep = "")
      }
      break
    }

    batch_flags <- numeric(batch.size)
    for (i in seq_len(batch.size)) {
      batch_flags[i] <- .run_single_consistency_split(df.x, N.x, hr.consistency, cox_init)
    }

    n_batch_valid <- sum(!is.na(batch_flags))
    n_batch_success <- sum(batch_flags == 1, na.rm = TRUE)

    n_total_valid <- n_total_valid + n_batch_valid
    n_total_success <- n_total_success + n_batch_success
  }

  # ---------------------------------------------------------------------------
  # Final evaluation
  # ---------------------------------------------------------------------------

  if (final_decision == "continue") {
    final_decision <- .early_stop_decision(
      n_success = n_total_success,
      n_total = n_total_valid,
      threshold = pconsistency.threshold,
      conf.level = conf.level,
      min_samples = 20
    )
  }

  p.consistency <- tryCatch({
    round(n_total_success / n_total_valid, pconsistency.digits)
  }, error = function(e) {
    return(NA_real_)
  })

  if (is.na(p.consistency)) {
    return(NULL)
  }

  if (final_decision == "fail" || p.consistency < pconsistency.threshold) {
    if (details) {
      cat("Subgroup ", m, ": FAILED (Pcons=", p.consistency, ")\n", sep = "")
    }
    return(NULL)
  }

  # ---------------------------------------------------------------------------
  # Format and return result
  # ---------------------------------------------------------------------------

  k <- length(this.m)
  covsm <- rep("M", maxk)
  mindex <- seq_len(maxk)
  Mnames <- paste(covsm, mindex, sep = ".")

  mfound <- matrix(rep("", maxk))
  mfound[seq_len(k)] <- this.m_label

  resultk <- c(
    p.consistency,
    found.hrs$HR[m],
    found.hrs$n[m],
    found.hrs$E[m],
    found.hrs$grp[m],
    m,
    k,
    mfound
  )

  names(resultk) <- c("Pcons", "hr", "N", "E", "g", "m", "K", Mnames)

  if (details) {
    cat("Subgroup ", m, ": PASSED (Pcons=", p.consistency,
        ", n_splits=", n_total_valid, ")\n", sep = "")
  }

  return(resultk)
}
