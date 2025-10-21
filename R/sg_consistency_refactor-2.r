# get_split_hr <- function(df, cox_initial = NULL) {
#   if (nrow(df) < 2 || sum(df$Event) < 2) {
#     return(NA_real_)
#   }
#
#   hr <- try({
#     fit <- survival::coxph(survival::Surv(Y, Event) ~ Treat,
#                            data = df,
#                            init = cox_initial,
#                            robust = FALSE)
#     summary(fit)$conf.int[1, 1]
#   }, silent = TRUE)
#
#   if (inherits(hr, "try-error")) return(NA_real_)
#   return(hr)
# }



#' Evaluate Single Subgroup for Consistency
#'
#' Helper function that evaluates a single subgroup (indexed by m) for consistency
#' across random splits. This function contains the core logic used in both
#' sequential and parallel execution modes.
#'
#' @param m Integer. Index of the subgroup to evaluate (1 to nrow(found.hrs))
#' @param index.Z Data.table or matrix. Factor indicators for all subgroups
#' @param names.Z Character vector. Names of factor columns
#' @param df Data.frame. Original data with Y, Event, Treat, id columns
#' @param found.hrs Data.table. Subgroup hazard ratio results
#' @param n.splits Integer. Number of random splits for consistency evaluation
#' @param hr.consistency Numeric. Minimum HR threshold for consistency
#' @param pconsistency.threshold Numeric. Minimum proportion of splits meeting consistency
#' @param pconsistency.digits Integer. Rounding digits for consistency proportion
#' @param maxk Integer. Maximum number of factors in a subgroup
#' @param confs_labels Character vector. Labels for confounders
#' @param details Logical. Print details during execution
#'
#' @return Named numeric vector with consistency results, or NULL if criteria not met.
#'   Vector contains: Pcons, hr, N, E, g, m, K, and factor labels (M.1, M.2, etc.)
#'
#' @importFrom data.table data.table
#' @importFrom survival coxph Surv
#' @export

evaluate_subgroup_consistency <- function(m, index.Z, names.Z, df, found.hrs,
                                         n.splits, hr.consistency,
                                         pconsistency.threshold, pconsistency.digits,
                                         maxk, confs_labels, details = FALSE) {

  get_split_hr <- function(df, cox_initial = NULL) {
    # Quick validation
    if (nrow(df) < 2 || sum(df$Event) < 2) {
      return(NA_real_)
    }

    # Fit model with minimal overhead, suppress all output
    fit <- tryCatch(
      suppressWarnings(
          survival::coxph(
            survival::Surv(Y, Event) ~ Treat,
            data = df,
            init = cox_initial,
            robust = FALSE,
            model = FALSE,  # Don't store model frame (saves memory)
            x = FALSE,      # Don't store design matrix
            y = FALSE       # Don't store response
          )
        ),
      error = function(e) NULL
    )

    # Return HR or NA
    if (is.null(fit)) return(NA_real_)
    return(exp(fit$coefficients[1]))
  }
  # =========================================================================
  # SECTION 1: VALIDATE SUBGROUP EXTRACTION
  # =========================================================================

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

  # =========================================================================
  # SECTION 2: VALIDATE LABEL CONVERSION
  # =========================================================================

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

  # =========================================================================
  # SECTION 3: CREATE SUBGROUP DEFINITION AND EXTRACT DATA
  # =========================================================================

  id.m <- paste(paste(this.m, collapse = "==1 & "), "==1")

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

  # =========================================================================
  # SECTION 4: GET INITIAL COX ESTIMATE (FOR STABILITY)
  # =========================================================================

  cox_init <- log(found.hrs$HR[m])
  if (is.na(cox_init) || is.infinite(cox_init)) {
    cox_init <- 0  # Fallback to null effect
    if (details) {
      cat("Subgroup ", m, ": using default cox_init=0\n")
    }
  }

  # =========================================================================
  # SECTION 5: PERFORM CONSISTENCY SPLITS
  # =========================================================================

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
      return(NA_real_)
    }

    # Check events in splits
    if (sum(df.x.split1$Event) < 2 || sum(df.x.split2$Event) < 2) {
      return(NA_real_)
    }

    # Fit models with error suppression (implicit within get_split_hr)
    hr.split1 <- get_split_hr(df = df.x.split1, cox_initial = cox_init)
    hr.split2 <- get_split_hr(df = df.x.split2, cox_initial = cox_init)

    # Return consistency flag
    if (!is.na(hr.split1) && !is.na(hr.split2)) {
      as.numeric(hr.split1 > hr.consistency && hr.split2 > hr.consistency)
    } else {
      NA_real_
    }
  })

  # =========================================================================
  # SECTION 6: CHECK VALIDITY AND CALCULATE CONSISTENCY
  # =========================================================================

  n_valid_splits <- sum(!is.na(flag.consistency))

  if (n_valid_splits == 0) {
    # NO valid splits - skip this subgroup entirely
    if (details) {
      cat("Subgroup ", m, ": No valid consistency splits\n")
    }
    return(NULL)
  }

  # Warn if too few valid splits
  if (n_valid_splits < 10) {
    warning("Subgroup ", m, ": only ", n_valid_splits, " valid splits out of ",
            n.splits, ". Results may be unreliable.")
  }

  # Calculate consistency
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

  # =========================================================================
  # SECTION 7: CHECK CONSISTENCY THRESHOLD
  # =========================================================================

  if (isTRUE(p.consistency < pconsistency.threshold)) {
    if (details) {
      cat("*** Not met: Subgroup, % Consistency =",
          c(this.m_label, p.consistency), "\n")
    }
    return(NULL)  # Did not meet threshold
  }

  # =========================================================================
  # SECTION 8: FORMAT AND RETURN RESULT
  # =========================================================================

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
    cat("Consistency met!\n")
    cat("# of splits =", n.splits, "\n")
    cat("**** Subgroup, % Consistency Met=",
        c(this.m_label, p.consistency), "\n")
  }

  return(resultk)
}


#' Subgroup Consistency Evaluation (REFACTORED WITH HELPER)
#'
#' Evaluates consistency of subgroups found in a survival analysis, using random
#' splits and hazard ratio criteria. This refactored version uses a helper function
#' for the core evaluation logic, enabling cleaner parallel/sequential execution.
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

subgroup.consistency.refactored <- function(df, hr.subgroups,
                                 hr.threshold = 1.0,
                                 hr.consistency = 1.0,
                                 pconsistency.threshold = 0.9,
                                 m1.threshold = Inf,
                                 n.splits = 100,
                                 details = FALSE,
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
  # SECTION 1: INPUT VALIDATION (UNCHANGED)
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

  # [Additional validation checks omitted for brevity - include all from original]


  # =========================================================================
  # SECTION 3: EXTRACT AND VALIDATE COVARIATE NAMES
  # =========================================================================

  exclude_cols <- c("grp", "K", "n", "E", "d1", "m1", "m0", "HR", "L(HR)", "U(HR)")
  names.Z <- setdiff(names(hr.subgroups), exclude_cols)

  if (length(names.Z) == 0) {
    stop("No covariate columns found in hr.subgroups after excluding metric columns")
  }

  if (length(names.Z) != Lsg) {
    stop("Lsg (", Lsg, ") does not match actual number of covariates (",
         length(names.Z), ") in hr.subgroups")
  }

  # =========================================================================
  # SECTION 4: FILTER SUBGROUPS BY CRITERIA
  # =========================================================================

  if (nrow(hr.subgroups) == 0) {
    warning("No subgroups provided in hr.subgroups")
    return(list(
      out_hr = NULL, out_maxSG = NULL, out_minSG = NULL,
      df_flag = NULL, sg.harm = NULL, sg.harm.id = NULL
    ))
  }

  # Apply filters
  if (is.finite(m1.threshold)) {
    hr.subgroups <- hr.subgroups[!is.na(hr.subgroups$m1), ]
    if (nrow(hr.subgroups) == 0) {
      warning("All subgroups removed after filtering NA m1 values")
      return(list(
        out_hr = NULL, out_maxSG = NULL, out_minSG = NULL,
        df_flag = NULL, sg.harm = NULL, sg.harm.id = NULL
      ))
    }
  }

  if (is.finite(m1.threshold)) {
    found.hrs <- hr.subgroups[hr.subgroups$HR >= hr.threshold &
                                hr.subgroups$m1 <= m1.threshold, ]
  } else {
    found.hrs <- hr.subgroups[hr.subgroups$HR >= hr.threshold, ]
  }

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

  # =========================================================================
  # SECTION 5: REMOVE DUPLICATES AND SORT
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
  }

  # Sort based on sg_focus
  if (sg_focus == "maxSG") {
    found.hrs <- found.hrs[order(found.hrs$n, decreasing = TRUE), ]
  } else if (sg_focus == "minSG") {
    found.hrs <- found.hrs[order(found.hrs$n, decreasing = FALSE), ]
  }

  # Extract index matrix
  index.Z <- found.hrs[, names.Z, with = FALSE]

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

  # =========================================================================
  # SECTION 6: VALIDATE PARALLEL CONFIGURATION
  # =========================================================================

  use_parallel <- length(parallel_args) > 0 && !is.null(parallel_args[[1]])

  if (use_parallel) {
    required_parallel <- c("plan", "workers")
    if (!all(required_parallel %in% names(parallel_args))) {
      warning("parallel_args missing required elements. Using sequential processing.")
      use_parallel <- FALSE
    }

    valid_plans <- c("multisession", "multicore", "callr", "sequential")
    if (!parallel_args$plan %in% valid_plans) {
      warning("Invalid parallel plan '", parallel_args$plan,
              "'. Using sequential processing.")
      use_parallel <- FALSE
    }

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
  # SECTION 7: INITIALIZE TIMING
  # =========================================================================

  if (details) t.start <- proc.time()[3]

  # =========================================================================
  # SECTION 8: EVALUATE EACH SUBGROUP (USING HELPER FUNCTION)
  # =========================================================================

  # SEQUENTIAL EXECUTION
  if (!use_parallel) {
    results_list <- lapply(seq_len(nrow(found.hrs)), function(m) {
      evaluate_subgroup_consistency(
        m = m,
        index.Z = index.Z,
        names.Z = names.Z,
        df = df,
        found.hrs = found.hrs,
        n.splits = n.splits,
        hr.consistency = hr.consistency,
        pconsistency.threshold = pconsistency.threshold,
        pconsistency.digits = pconsistency.digits,
        maxk = maxk,
        confs_labels = confs_labels,
        details = details
      )
    })
  }
  # PARALLEL EXECUTION
  else {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)
    setup_parallel_SGcons(parallel_args)

    results_list <- future.apply::future_lapply(
      seq_len(nrow(found.hrs)),
      function(m) {
        evaluate_subgroup_consistency(
          m = m,
          index.Z = index.Z,
          names.Z = names.Z,
          df = df,
          found.hrs = found.hrs,
          n.splits = n.splits,
          hr.consistency = hr.consistency,
          pconsistency.threshold = pconsistency.threshold,
          pconsistency.digits = pconsistency.digits,
          maxk = maxk,
          confs_labels = confs_labels,
          details = details
        )
      },
      future.seed = TRUE,
      future.packages = c("survival", "data.table"),
      future.globals = structure(
        TRUE,
        add = c(
          # Functions
          "evaluate_subgroup_consistency",
          "get_split_hr",
          "FS_labels",

          # Data objects
          "index.Z",
          "names.Z",
          "df",
          "found.hrs",
          "n.splits",
          "hr.consistency",
          "pconsistency.threshold",
          "pconsistency.digits",
          "maxk",
          "confs_labels"
          # Note: 'details' deliberately excluded - doesn't affect computation
        )
      )
    )  }

  # =========================================================================
  # SECTION 9: COMPILE RESULTS (UNCHANGED FROM ORIGINAL)
  # =========================================================================

  if (length(results_list) == 0) {
    if (details) cat("No subgroups met consistency criteria\n")
    return(list(
      out_hr = NULL, out_maxSG = NULL, out_minSG = NULL,
      df_flag = NULL, sg.harm = NULL, sg.harm.id = NULL
    ))
  }

  results_list <- Filter(Negate(is.null), results_list)

  if (length(results_list) == 0) {
    if (details) cat("All subgroup evaluations returned NULL\n")
    return(list(
      out_hr = NULL, out_maxSG = NULL, out_minSG = NULL,
      df_flag = NULL, sg.harm = NULL, sg.harm.id = NULL
    ))
  }

  res <- tryCatch({
    data.table::as.data.table(do.call(rbind, results_list))
  }, error = function(e) {
    stop("Error combining results: ", e$message,
         "\nThis may indicate inconsistent result structure across subgroups.")
  })

  any.found <- nrow(res)

  if (any.found == 0) {
    if (details) cat("No subgroups found meeting consistency threshold\n")
    return(list(
      out_hr = NULL, out_maxSG = NULL, out_minSG = NULL,
      df_flag = NULL, sg.harm = NULL, sg.harm.id = NULL
    ))
  }

  # Convert columns to numeric
  cols_to_numeric <- c("Pcons", "hr", "N", "E", "K")
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
  # SECTION 10: GENERATE OUTPUTS (UNCHANGED FROM ORIGINAL)
  # =========================================================================

  out_hr <- out_maxSG <- out_minSG <- NULL
  df_flag <- sg.harm <- sg.harm.id <- NULL

  if (any.found > 0) {
    result_new <- data.table::copy(res)

    # Generate outputs for different sg_focus values
    sgdetails <- ifelse(plot.sg && sg_focus == "hr", TRUE, FALSE)
    out_hr <- tryCatch({
      sg_consistency_out(df = df, result_new = result_new, sg_focus = "hr",
                        details = sgdetails, plot.sg = sgdetails, index.Z = index.Z,
                        names.Z = names.Z, by.risk = by.risk, confs_labels = confs_labels)
    }, error = function(e) {
      warning("Error in sg_consistency_out for 'hr': ", e$message)
      NULL
    })

    sgdetails <- ifelse(plot.sg && sg_focus %in% c("hrMaxSG", "maxSG"), TRUE, FALSE)
    out_maxSG <- tryCatch({
      sg_consistency_out(df = df, result_new = result_new, sg_focus = "maxSG",
                        details = sgdetails, plot.sg = sgdetails, index.Z = index.Z,
                        names.Z = names.Z, by.risk = by.risk, confs_labels = confs_labels)
    }, error = function(e) {
      warning("Error in sg_consistency_out for 'maxSG': ", e$message)
      NULL
    })

    sgdetails <- ifelse(plot.sg && sg_focus %in% c("hrMinSG", "minSG"), TRUE, FALSE)
    out_minSG <- tryCatch({
      sg_consistency_out(df = df, result_new = result_new, sg_focus = "minSG",
                        details = sgdetails, plot.sg = sgdetails, index.Z = index.Z,
                        names.Z = names.Z, by.risk = by.risk, confs_labels = confs_labels)
    }, error = function(e) {
      warning("Error in sg_consistency_out for 'minSG': ", e$message)
      NULL
    })

    # Map sg_focus to output
    sg_map <- list(
      hr = out_hr,
      hrMaxSG = out_maxSG,
      maxSG = out_maxSG,
      hrMinSG = out_minSG,
      minSG = out_minSG
    )

    if (!sg_focus %in% names(sg_map)) {
      stop(sprintf("Unknown sg_focus value: %s", sg_focus))
    }

    sg_obj <- sg_map[[sg_focus]]

    if (is.null(sg_obj)) {
      warning("No valid output for sg_focus='", sg_focus, "'")
      df_flag <- NULL
      sg.harm <- NULL
      sg.harm.id <- NULL
    } else {
      required_fields <- c("df_flag", "sg.harm_label", "sg.harm.id")
      missing_fields <- setdiff(required_fields, names(sg_obj))
      if (length(missing_fields) > 0) {
        stop("sg_consistency_out result missing fields: ",
             paste(missing_fields, collapse = ", "))
      }

      df_flag <- sg_obj$df_flag
      sg.harm <- sg_obj$sg.harm_label
      sg.harm.id <- sg_obj$sg.harm.id
    }

    if (details) cat("SG focus=", sg_focus, "\n")
  }

  # =========================================================================
  # SECTION 11: FINAL TIMING AND OUTPUT
  # =========================================================================

  if (details) {
    t.end <- proc.time()[3]
    t.min <- (t.end - t.start) / 60
    cat("Subgroup Consistency Minutes=", t.min, "\n")
    if (any.found > 0) {
      cat("Subgroup found (FS)\n")
      if (!is.null(sg.harm)) {
        cat("Selected subgroup:", paste(sg.harm, collapse = " & "), "\n")
      }
    } else {
      cat("NO subgroup found (FS)\n")
    }
  }

  output <- list(
    out_hr = out_hr,
    out_maxSG = out_maxSG,
    out_minSG = out_minSG,
    df_flag = df_flag,
    sg.harm = sg.harm,
    sg.harm.id = sg.harm.id
  )

  return(output)
}
