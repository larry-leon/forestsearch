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

#' Get Hazard Ratio from Split Data
#'
#' Calculates the hazard ratio from a split data set using Cox regression.
#'
#' @param df Data frame with survival data.
#' @return Numeric hazard ratio or NA if error.
#' @importFrom survival coxph Surv
#' @export

get_split_hr <- function(df) {
  hr <- try((summary(coxph(Surv(Y, Event) ~ Treat, data = df, robust = FALSE))$conf.int), TRUE)
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


#' Subgroup Consistency Evaluation
#'
#' Evaluates consistency of subgroups found in a survival analysis, using random splits and hazard ratio criteria.
#'
#' @param df The original data.frame.
#' @param hr.subgroups Data.table of subgroup hazard ratio results.
#' @param hr.threshold Minimum hazard ratio for subgroup inclusion.
#' @param hr.consistency Minimum hazard ratio for consistency in splits.
#' @param pconsistency.threshold Minimum proportion of splits meeting consistency.
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

subgroup.consistency <- function(df, hr.subgroups, hr.threshold = 1.0, hr.consistency = 1.0, pconsistency.threshold = 0.9, m1.threshold = Inf, n.splits = 100,
details = FALSE, stop.threshold = 1.1, by.risk = 12, plot.sg = FALSE, maxk = 7, Lsg, confs_labels, sg_focus = "hr", stop_Kgroups = 10,
checking = FALSE, parallel_args = list(NULL)) {
  names.Z <- c(names(hr.subgroups[, -c("grp", "K", "n", "E", "d1", "m1", "m0", "HR", "L(HR)", "U(HR)")]))
  if (length(names.Z) != Lsg) stop("HR subgroup results not matching L, check subgroup search function")
  if (m1.threshold == Inf) {
    found.hrs <- hr.subgroups[hr.subgroups$HR >= hr.threshold, ]
  } else {
    hr.subgroups <- hr.subgroups[which(!is.na(hr.subgroups$m1)), ]
    found.hrs <- hr.subgroups[hr.subgroups$HR >= hr.threshold & hr.subgroups$m1 <= m1.threshold, ]
  }
  if (nrow(found.hrs) > 1) {
    found.hrs <- remove_redundant_subgroups(found.hrs)
  }
  if (sg_focus == "maxSG") found.hrs <- found.hrs[order(found.hrs$n, decreasing = TRUE), ]
  if (sg_focus == "minSG") found.hrs <- found.hrs[order(found.hrs$n, decreasing = FALSE), ]
  index.Z <- found.hrs[, -c("grp", "K", "n", "E", "d1", "m1", "m0", "HR", "L(HR)", "U(HR)")]
  if (dim(index.Z)[2] != Lsg) stop("HR subgroup results not matching L, check subgroup search function")

if(details) cat("# of initial candidates",c(nrow(found.hrs)),"\n")

  maxsgs <- min(c(nrow(found.hrs),stop_Kgroups))
  # Restrict to top "stop_Kgroups"
  found.hrs <- found.hrs[c(1:maxsgs)]

if(details) cat("# of candidates restricted to 'top 10'", c(nrow(found.hrs)),"\n")

  res <- NULL
  any.found <- 0
  resultk <- NULL
  if (details) t.start <- proc.time()[3]
  # Standard lapply
  if (length(parallel_args) == 0) {
    results_list <- lapply(seq_len(nrow(found.hrs)), function(m) {
      # --- Begin: code for each subgroup ---
      indexm <- as.numeric(unlist(index.Z[m, ]))
      this.m <- names.Z[indexm == 1]
      this.m_label <- unlist(lapply(this.m, FS_labels, confs_labels = confs_labels))
      id.m <- paste(paste(this.m, collapse = "==1 & "), "==1")
      df.sub <- subset(df, eval(parse(text = id.m)))
      df.x <- data.table::data.table(df.sub)
      set.seed(8316951)
      N.x <- nrow(df.x)
      flag.consistency <- sapply(seq_len(n.splits), function(bb) {
        in.split1 <- sample(c(TRUE, FALSE), N.x, replace = TRUE, prob = c(0.5, 0.5))
        df.x$insplit1 <- in.split1
        df.x.split1 <- subset(df.x, insplit1 == 1)
        df.x.split2 <- subset(df.x, insplit1 == 0)
        hr.split1 <- suppressWarnings(get_split_hr(df.x.split1))
        hr.split2 <- suppressWarnings(get_split_hr(df.x.split2))
        if (!is.na(hr.split1) && !is.na(hr.split2)) {
          as.numeric(hr.split1 > hr.consistency & hr.split2 > hr.consistency)
        } else {
          NA_real_
        }
      })
      p.consistency <- mean(flag.consistency, na.rm = TRUE)
      if (p.consistency < pconsistency.threshold && details){
        cat("*** Not met: Subgroup, % Consistency =", c(this.m_label, p.consistency), "\n")
      }
      if (p.consistency >= pconsistency.threshold) {
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
      # --- End: code for each subgroup ---
      return(NULL)
    })
  } else {
    old_plan <- future::plan()
    on.exit(future::plan(old_plan), add = TRUE)  # Restore plan on exit
    setup_parallel_SGcons(parallel_args)
    results_list <- future.apply::future_lapply(seq_len(nrow(found.hrs)), function(m) {
      # --- Begin: code for each subgroup ---
      indexm <- as.numeric(unlist(index.Z[m, ]))
      this.m <- names.Z[indexm == 1]
      this.m_label <- unlist(lapply(this.m, FS_labels, confs_labels = confs_labels))
      id.m <- paste(paste(this.m, collapse = "==1 & "), "==1")
      df.sub <- subset(df, eval(parse(text = id.m)))
      df.x <- data.table::data.table(df.sub)
      set.seed(8316951)
      N.x <- nrow(df.x)
      flag.consistency <- sapply(seq_len(n.splits), function(bb) {
        in.split1 <- sample(c(TRUE, FALSE), N.x, replace = TRUE, prob = c(0.5, 0.5))
        df.x$insplit1 <- in.split1
        df.x.split1 <- subset(df.x, insplit1 == 1)
        df.x.split2 <- subset(df.x, insplit1 == 0)
        hr.split1 <- suppressWarnings(get_split_hr(df.x.split1))
        hr.split2 <- suppressWarnings(get_split_hr(df.x.split2))
        if (!is.na(hr.split1) && !is.na(hr.split2)) {
          as.numeric(hr.split1 > hr.consistency & hr.split2 > hr.consistency)
        } else {
          NA_real_
        }
      })
      p.consistency <- mean(flag.consistency, na.rm = TRUE)
      if (p.consistency < pconsistency.threshold & details){
        cat("*** Not met: Subgroup, % Consistency =", c(this.m_label, p.consistency), "\n")
      }
      if (p.consistency >= pconsistency.threshold) {
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
      # --- End: code for each subgroup ---
      return(NULL)
    }, future.seed = TRUE)
  }
  res <- data.table::as.data.table(do.call(rbind, results_list))
  any.found <- nrow(res)
  out_hr <- out_maxSG <- out_minSG <- NULL
  df_flag <- sg.harm <- sg.harm.id <- NULL
  if (any.found > 0) {
    cols_to_numeric <- c("Pcons", "hr", "N", "E", "K")
    res[, (cols_to_numeric) := lapply(.SD, as.numeric), .SDcols = cols_to_numeric]
    result_new <- data.table::copy(res)

    sgdetails <- ifelse(plot.sg & sg_focus == "hr", TRUE, FALSE)
    out_hr <- sg_consistency_out(df = df, result_new = result_new, sg_focus = "hr", details = sgdetails,
                                 plot.sg = sgdetails, index.Z = index.Z, names.Z = names.Z, by.risk = by.risk, confs_labels = confs_labels)
    sgdetails <- ifelse(plot.sg & sg_focus %in% c("hrMaxSG", "maxSG"), TRUE, FALSE)
    out_maxSG <- sg_consistency_out(df = df, result_new = result_new, sg_focus = "maxSG", details = sgdetails,
                                    plot.sg = sgdetails, index.Z = index.Z, names.Z = names.Z, by.risk = by.risk, confs_labels = confs_labels)
    sgdetails <- ifelse(plot.sg & sg_focus %in% c("hrMinSG", "minSG"), TRUE, FALSE)
    out_minSG <- sg_consistency_out(df = df, result_new = result_new, sg_focus = "minSG", details = sgdetails,
                                    plot.sg = sgdetails, index.Z = index.Z, names.Z = names.Z, by.risk = by.risk, confs_labels = confs_labels)

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
    # Assign variables
    df_flag    <- sg_obj$df_flag
    sg.harm    <- sg_obj$sg.harm_label
    sg.harm.id <- sg_obj$sg.harm.id

    if (details) cat("SG focus=", c(sg_focus), "\n")
  }
  if (details) {
    t.end <- proc.time()[3]
    t.min <- (t.end - t.start) / 60
    cat("Subgroup Consistency Minutes=", c(t.min), "\n")
    if (any.found > 0) cat("Subgroup found (FS)\n")
    if (any.found == 0) cat("NO subgroup found (FS)\n")
  }
  list(out_hr = out_hr, out_maxSG = out_maxSG, out_minSG = out_minSG, df_flag = df_flag, sg.harm = sg.harm, sg.harm.id = sg.harm.id)
}
