# minp = Minimum prevalence rate
# Only combinations where all factors have prevalance at least minp are evaluated
# Minimum difference in subgroup sample size
# Max subgroup size is n: After first factor x1 is added, the sample size for
# the subgroup is then n(x1), say.   If combining with next factor x2
# does not reduce the sample size by rmin, then consider this combination "redundant"
# Adding d.min to require min(Events) per treatment arm (d0.min for control, d1.min for treatment)
# covs.in[isin,]<-unlist(lapply(covs.in[isin,],function(x){ifelse(x==1,0,1)})) # slower



#' Get all combinations of subgroup factors up to maxk
#'
#' Generates all possible combinations of subgroup factors up to a specified maximum size.
#'
#' @param L Integer. Number of subgroup factors.
#' @param maxk Integer. Maximum number of factors in a combination.
#' @return List with \code{max_count} (total combinations) and \code{indices_list} (indices for each k).
#' @importFrom utils combn
#' @export

get_combinations_info <- function(L, maxk) {
  # Calculate max_count
  max_count <- 0
  indices_list <- list()
  for (k in 1:maxk) {
    indices_k <- base::t(combn(L, k))
    indices_list[[paste0("k", k)]] <- indices_k
    max_count <- max_count + nrow(indices_k)
  }
  # Consistency check
  tot_counts <- sum(sapply(indices_list, nrow))
  if (tot_counts != max_count) stop("Error with maximum combinations for kmax = ", maxk)
  list(
    max_count = max_count,
    indices_list = indices_list
  )
}


#' Get subgroup membership vector (old, legacy version)
#'
#' Returns a vector indicating subgroup membership (1 if all selected factors are present, 0 otherwise).
#'
#' @param zz Matrix or data frame of subgroup factor indicators.
#' @param covs.in Numeric vector indicating which factors are selected (1 = included).
#' @return Numeric vector of subgroup membership (1/0).
#' @export

get_subgroup_membership_old <- function(zz, covs.in) {
  x <- zz[, which(covs.in == 1), drop = FALSE]
  if (ncol(x) == 0) return(rep(1, nrow(zz)))
  apply(x, 1, function(row) all(row == 1))
}

#' Get subgroup membership vector
#'
#' Returns a vector indicating subgroup membership (1 if all selected factors are present, 0 otherwise).
#'
#' @param zz Matrix or data frame of subgroup factor indicators.
#' @param covs.in Numeric vector indicating which factors are selected (1 = included).
#' @return Numeric vector of subgroup membership (1/0).
#' @export

get_subgroup_membership <- function(zz, covs.in) {
  selected_cols <- which(covs.in == 1)
  if (length(selected_cols) == 0) return(rep(1, nrow(zz)))
  x <- zz[, selected_cols, drop = FALSE]
  rowSums(x) == length(selected_cols)
}


#' Get indicator vector for selected subgroup factors
#'
#' Returns a vector indicating which factors are included in a subgroup combination.
#'
#' @param kk Integer. Index of the combination.
#' @param maxk Integer. Maximum number of factors in a combination.
#' @param L Integer. Number of subgroup factors.
#' @param counts_1factor Integer. Number of single-factor combinations.
#' @param index_1factor Matrix of indices for single-factor combinations.
#' @param counts_2factor Integer. Number of two-factor combinations.
#' @param index_2factor Matrix of indices for two-factor combinations.
#' @param counts_3factor Integer. Number of three-factor combinations.
#' @param index_3factor Matrix of indices for three-factor combinations.
#' @return Numeric vector indicating selected factors (1 = included, 0 = not included).
#' @export

get_covs_in <- function(kk, maxk, L,
                        counts_1factor, index_1factor,
                        counts_2factor = NULL, index_2factor = NULL,
                        counts_3factor = NULL, index_3factor = NULL) {
  covs.in <- rep(0, L)
  if (kk <= counts_1factor) {
    which1 <- index_1factor[kk]
    covs.in[which1] <- 1.0
  } else if (maxk == 2 && kk > counts_1factor) {
    kk_new <- kk - counts_1factor
    which1 <- index_2factor[kk_new, 1]
    which2 <- index_2factor[kk_new, 2]
    covs.in[c(which1, which2)] <- 1.0
  } else if (maxk == 3) {
    if (kk > counts_1factor && kk <= (counts_1factor + counts_2factor)) {
      kk_new <- kk - counts_1factor
      which1 <- index_2factor[kk_new, 1]
      which2 <- index_2factor[kk_new, 2]
      covs.in[c(which1, which2)] <- 1.0
    }
    if (kk > (counts_1factor + counts_2factor)) {
      kk_new <- kk - (counts_1factor + counts_2factor)
      which1 <- index_3factor[kk_new, 1]
      which2 <- index_3factor[kk_new, 2]
      which3 <- index_3factor[kk_new, 3]
      covs.in[c(which1, which2, which3)] <- 1.0
    }
  }
  covs.in
}

#' Extract redundancy flag for subgroup combinations
#'
#' Checks if adding each factor to a subgroup reduces the sample size by at least \code{rmin}.
#'
#' @param x Matrix of subgroup factor indicators.
#' @param rmin Integer. Minimum required reduction in sample size.
#' @return List with \code{id.x} (membership vector) and \code{flag.redundant} (logical).
#' @export

extract_idx_flagredundancy <- function(x, rmin) {
  n <- nrow(x)
  id.x <- rep(1, n)
  flag.redundant <- FALSE
  nx.last <- n
  for (m in 1:ncol(x)) {
    if (!flag.redundant) {
      id.x <- id.x * x[, m]
      nx.this <- sum(id.x)
      if (nx.last - nx.this <= rmin) flag.redundant <- TRUE
      nx.last <- nx.this
    }
  }
  list(id.x = id.x, flag.redundant = flag.redundant)
}

#' Subgroup search for treatment effect heterogeneity
#'
#' Searches for subgroups with treatment effect heterogeneity using combinations of candidate factors.
#' Evaluates subgroups for minimum prevalence, event counts, and hazard ratio threshold.
#'
#' @param Y Numeric vector of outcome (e.g., time-to-event).
#' @param Event Numeric vector of event indicators (0/1).
#' @param Treat Numeric vector of treatment group indicators (0/1).
#' @param ID Optional vector of subject IDs.
#' @param Z Matrix or data frame of candidate subgroup factors (binary indicators).
#' @param n.min Integer. Minimum subgroup size.
#' @param d0.min Integer. Minimum number of events in control.
#' @param d1.min Integer. Minimum number of events in treatment.
#' @param hr.threshold Numeric. Hazard ratio threshold for subgroup selection.
#' @param max.minutes Numeric. Maximum minutes for search.
#' @param minp Numeric. Minimum prevalence rate for each factor.
#' @param rmin Integer. Minimum required reduction in sample size when adding a factor.
#' @param details Logical. Print details during execution.
#' @param maxk Integer. Maximum number of factors in a subgroup.
#'
#' @return List with found subgroups, maximum HR, search time, and configuration info.
#'
#' @importFrom data.table data.table setorder
#' @importFrom survival coxph Surv survfit
#' @importFrom utils combn
#' @export

subgroup.search_legacy <- function(Y,Event,Treat,ID = NULL,Z,n.min = 30,d0.min = 15,d1.min = 15,
hr.threshold = 1.0,max.minutes = 30,minp = 0.05,rmin = 5,details = FALSE,maxk = 2){
  temp <- cbind(Y,Event,Treat,Z)
  temp <- na.exclude(temp)
  yy <- temp[,1]
  dd <- temp[,2]
  tt <- temp[,3]
  zz <- temp[,-c(1,2,3)]
  n <- length(yy)
  L <- ncol(zz)
  covs.in <- as.matrix(rep(0,L)) # Initially no members
  # Store results appending to HR.model.k
  HR.model.k <- NULL
  limit_all<-(2**L)-1 # All-possible configurations including redundancies
  # Maximum number of combinations among factors <= maxk
  if(maxk == 1) max_count <- L
  if(maxk == 2){
  # choose(L,2)+choose(L,1)
  max_count <- (L*(L-1)/2) + L
  }
  if(maxk==3){
  choose(L,3)+choose(L,2)+choose(L,1)
  max_count <- ((L*(L-2)*(L-1))/6)+(L*(L-1)/2) + L
  }
  if(maxk == 1){
    index_1factor <- t(combn(L, 1))
    counts_1factor <- nrow(index_1factor)
    tot_counts <- counts_1factor
    if(tot_counts != max_count) stop("Error with maximum combinations kmax=2")
  }
  if(maxk == 2){
  index_2factor <- t(combn(L, 2))
  counts_2factor <- nrow(index_2factor)
  index_1factor <- t(combn(L, 1))
  counts_1factor <- nrow(index_1factor)
  tot_counts <- counts_2factor+counts_1factor
  if(tot_counts != max_count) stop("Error with maximum combinations kmax=2")
  }
  if(maxk == 3){
    index_3factor <- t(combn(L, 3))
    counts_3factor <- nrow(index_3factor)
    index_2factor <- t(combn(L, 2))
    counts_2factor <- nrow(index_2factor)
    index_1factor <- t(combn(L, 1))
    counts_1factor <- nrow(index_1factor)
    tot_counts <- counts_3factor+counts_2factor+counts_1factor
    if(tot_counts != max_count) stop("Error with maximum combinations kmax=2")
  }
  if(details){
    cat("Number of possible configurations (<= maxk): maxk, # <= maxk",c(maxk,max_count),"\n")
  }
  t.start<-proc.time()[3]
  t.sofar<-0
  # Loop through single-factor subgroups first
  results_list <- lapply(seq_len(tot_counts), function(kk) {
    covs.in <- get_covs_in(
      kk, maxk, L,
      counts_1factor, index_1factor,
      counts_2factor, index_2factor,
      counts_3factor, index_3factor
    )
    k.in <- sum(covs.in)
    t.now <- proc.time()[3]
    t.sofar <- (t.now - t.start) / 60
    if (t.sofar > max.minutes) return(NULL)
    if (k.in > maxk) return(NULL)

    selected_cols <- which(covs.in == 1)
    x <- zz[, selected_cols, drop = FALSE]

    xpx <- t(x) %*% x
    if (!all(xpx > 0)) return(NULL)

    get_idx_rflag <- extract_idx_flagredundancy(x, rmin)

    id.x <- get_idx_rflag$id.x
    flag.redundant <- get_idx_rflag$flag.redundant

    rm("get_idx_rflag")

    d1 <- sum((dd * tt)[which(id.x == 1)])
    d0 <- sum((dd * (1 - tt))[which(id.x == 1)])
    nx <- sum(id.x)

    crit3 <- all(apply(x, 2, mean) >= minp)

    if (!flag.redundant && d1 >= d1.min && d0 >= d0.min && nx > n.min && crit3) {
      data.x <- data.table::data.table(Y = yy, E = dd, Treat = tt, id.x = id.x)
      df.x <- data.x[id.x == 1]
      hr.cox <- try(summary(coxph(Surv(Y, E) ~ Treat, data = df.x, robust = FALSE))$conf.int, TRUE)
      if (!inherits(hr.cox, "try-error")) {
        if (hr.cox[1] > hr.threshold) {
          meds <- summary(survfit(Surv(Y, E) ~ Treat, data = df.x))$table[, "median"]
          m0 <- meds[1]
          m1 <- meds[2]
          id.group <- kk
          hr.this <- c(id.group, sum(covs.in), nx, d0 + d1, d1, m1, m0, hr.cox[c(1, 3, 4)], covs.in)
          return(hr.this)
        }
      }
    }
    return(NULL)
  })

  # Filter out NULLs and combine results
  HR.model.k <- base::do.call(rbind, Filter(Negate(is.null), results_list))

   t.end<-proc.time()[3]
   t.min<-(t.end-t.start)/60

  if(details){
    cat("Events criteria for control,exp=",c(d0.min,d1.min),"\n")
    cat("*Subgroup Searching Minutes=*",c(t.min),"\n")
    }

  out.found<-NULL
  if(!is.null(HR.model.k)){
    hr.out <- data.table::data.table(HR.model.k)
    names(hr.out)<-c("grp","K","n","E","d1","m1","m0","HR","L(HR)","U(HR)",colnames(Z))
    rownames(hr.out)<-NULL
    hr.out <- data.table::setorder(hr.out,-HR,K)
    out.found <- list(hr.subgroups = hr.out)
  }
  if(is.null(out.found)){
  if(details)    cat("NO subgroup candidate found (FS)","\n")
  return(NULL)
  } else{
  if(details) cat("Subgroup candidate(s) found (FS)","\n")
  max_sg_est <- max(hr.out[["HR"]])
  return(list(out.found = out.found, max_sg_est = max_sg_est, time_search = t.sofar,L = L,max_count = max_count, prop_max_count = 1))
  }
}



#' Subgroup Search for Treatment Effect Heterogeneity (Improved)
#'
#' Searches for subgroups with treatment effect heterogeneity using combinations
#' of candidate factors. Evaluates subgroups for minimum prevalence, event counts,
#' and hazard ratio threshold.
#'
#' @param Y Numeric vector of outcome (e.g., time-to-event).
#' @param Event Numeric vector of event indicators (0/1).
#' @param Treat Numeric vector of treatment group indicators (0/1).
#' @param ID Optional vector of subject IDs.
#' @param Z Matrix or data frame of candidate subgroup factors (binary indicators).
#' @param n.min Integer. Minimum subgroup size.
#' @param d0.min Integer. Minimum number of events in control.
#' @param d1.min Integer. Minimum number of events in treatment.
#' @param hr.threshold Numeric. Hazard ratio threshold for subgroup selection.
#' @param max.minutes Numeric. Maximum minutes for search.
#' @param minp Numeric. Minimum prevalence rate for each factor.
#' @param rmin Integer. Minimum required reduction in sample size when adding a factor.
#' @param details Logical. Print details during execution.
#' @param maxk Integer. Maximum number of factors in a subgroup.
#'
#' @return List with found subgroups, maximum HR, search time, and configuration info.
#'
#' @importFrom data.table data.table setorder
#' @importFrom survival coxph Surv survfit
#' @importFrom utils combn
#' @export

subgroup.search <- function(Y, Event, Treat, ID = NULL, Z,
                            n.min = 30, d0.min = 15, d1.min = 15,
                            hr.threshold = 1.0, max.minutes = 30,
                            minp = 0.05, rmin = 5,
                            details = FALSE, maxk = 2) {

  # =========================================================================
  # SECTION 1: DATA PREPARATION AND VALIDATION
  # =========================================================================

  # Clean and prepare data
  prepared_data <- prepare_search_data(Y, Event, Treat, Z)
  yy <- prepared_data$Y
  dd <- prepared_data$Event
  tt <- prepared_data$Treat
  zz <- prepared_data$Z

  n <- length(yy)
  L <- ncol(zz)

  # =========================================================================
  # SECTION 2: GENERATE COMBINATION INDICES
  # =========================================================================

  combo_info <- generate_combination_indices(L, maxk)

  if (details) {
    cat("Number of possible configurations (<= maxk): maxk =", maxk,
        ", # combinations =", combo_info$max_count, "\n")
  }

  # =========================================================================
  # SECTION 3: SEARCH THROUGH COMBINATIONS
  # =========================================================================

  t.start <- proc.time()[3]

  results_list <- search_combinations(
    yy = yy, dd = dd, tt = tt, zz = zz,
    combo_info = combo_info,
    n.min = n.min, d0.min = d0.min, d1.min = d1.min,
    hr.threshold = hr.threshold, minp = minp, rmin = rmin,
    max.minutes = max.minutes, t.start = t.start,
    maxk = maxk, L = L
  )

  # =========================================================================
  # SECTION 4: COMPILE AND FORMAT RESULTS
  # =========================================================================

  t.end <- proc.time()[3]
  t.min <- (t.end - t.start) / 60
  t.sofar <- t.min

  if (details) {
    cat("Events criteria: control >=", d0.min, ", treatment >=", d1.min, "\n")
    cat("Subgroup search completed in", round(t.min, 2), "minutes\n")
  }

  # Format results
  output <- format_search_results(
    results_list = results_list,
    Z = Z,
    details = details,
    t.sofar = t.sofar,
    L = L,
    max_count = combo_info$max_count
  )

  return(output)
}


# =========================================================================
# HELPER FUNCTIONS FOR IMPROVED READABILITY
# =========================================================================

#' Prepare Data for Subgroup Search
#'
#' Cleans data by removing missing values and extracting components
#'
#' @keywords internal
prepare_search_data <- function(Y, Event, Treat, Z) {
  temp <- cbind(Y, Event, Treat, Z)
  temp <- na.exclude(temp)

  list(
    Y = temp[, 1],
    Event = temp[, 2],
    Treat = temp[, 3],
    Z = temp[, -c(1, 2, 3)]
  )
}


#' Generate Combination Indices
#'
#' Creates indices for all factor combinations up to maxk
#'
#' @keywords internal
generate_combination_indices <- function(L, maxk) {

  # Calculate maximum possible combinations
  max_count <- calculate_max_combinations(L, maxk)

  # Generate indices for each level
  indices <- list()
  counts <- numeric(maxk)

  for (k in 1:maxk) {
    indices[[k]] <- t(combn(L, k))
    counts[k] <- nrow(indices[[k]])
  }

  # Validate total count
  if (sum(counts) != max_count) {
    stop("Error: combination count mismatch for maxk = ", maxk)
  }

  list(
    max_count = max_count,
    indices_1 = indices[[1]],
    counts_1 = counts[1],
    indices_2 = if (maxk >= 2) indices[[2]] else NULL,
    counts_2 = if (maxk >= 2) counts[2] else 0,
    indices_3 = if (maxk >= 3) indices[[3]] else NULL,
    counts_3 = if (maxk >= 3) counts[3] else 0
  )
}


#' Calculate Maximum Combinations
#'
#' @keywords internal
calculate_max_combinations <- function(L, maxk) {
  if (maxk == 1) return(L)
  if (maxk == 2) return(L + (L * (L - 1) / 2))
  if (maxk == 3) return(L + (L * (L - 1) / 2) + (L * (L - 2) * (L - 1) / 6))
  stop("maxk must be 1, 2, or 3")
}


#' Search Through All Combinations
#'
#' Main search loop evaluating each factor combination
#'
#' @keywords internal
search_combinations <- function(yy, dd, tt, zz, combo_info,
                                n.min, d0.min, d1.min, hr.threshold,
                                minp, rmin, max.minutes, t.start,
                                maxk, L) {

  tot_counts <- combo_info$max_count
  results_list <- vector("list", tot_counts)
  n_results <- 0

  for (kk in seq_len(tot_counts)) {

    # Check time limit
    if (is_time_exceeded(t.start, max.minutes)) break

    # Get factor selection for this combination
    covs.in <- get_covs_in(
      kk, maxk, L,
      combo_info$counts_1, combo_info$indices_1,
      combo_info$counts_2, combo_info$indices_2,
      combo_info$counts_3, combo_info$indices_3
    )

    # Skip if too many factors
    if (sum(covs.in) > maxk) next

    # Evaluate this combination
    result <- evaluate_combination(
      covs.in = covs.in,
      yy = yy, dd = dd, tt = tt, zz = zz,
      n.min = n.min, d0.min = d0.min, d1.min = d1.min,
      hr.threshold = hr.threshold, minp = minp, rmin = rmin,
      kk = kk
    )

    # Store if valid result found
    if (!is.null(result)) {
      n_results <- n_results + 1
      results_list[[n_results]] <- result
    }
  }

  # Return only non-NULL results
  Filter(Negate(is.null), results_list[1:n_results])
}


#' Check if Time Limit Exceeded
#'
#' @keywords internal
is_time_exceeded <- function(t.start, max.minutes) {
  t.now <- proc.time()[3]
  t.elapsed <- (t.now - t.start) / 60
  return(t.elapsed > max.minutes)
}


#' Evaluate a Single Factor Combination
#'
#' Tests whether a specific combination meets all criteria
#'
#' @keywords internal
evaluate_combination <- function(covs.in, yy, dd, tt, zz,
                                 n.min, d0.min, d1.min, hr.threshold,
                                 minp, rmin, kk) {

  # Extract selected factors
  selected_cols <- which(covs.in == 1)
  x <- zz[, selected_cols, drop = FALSE]

  # Quick checks that can fail fast
  if (!has_positive_variance(x)) return(NULL)
  if (!meets_prevalence_threshold(x, minp)) return(NULL)

  # Check for redundancy
  redundancy_check <- extract_idx_flagredundancy(x, rmin)
  if (redundancy_check$flag.redundant) return(NULL)

  id.x <- redundancy_check$id.x

  # Check event counts
  event_counts <- calculate_event_counts(dd, tt, id.x)
  if (!meets_event_criteria(event_counts, d0.min, d1.min)) return(NULL)

  # Check sample size
  nx <- sum(id.x)
  if (nx <= n.min) return(NULL)

  # Fit Cox model and check HR threshold
  cox_result <- fit_cox_for_subgroup(yy, dd, tt, id.x)
  if (is.null(cox_result)) return(NULL)
  if (cox_result$hr <= hr.threshold) return(NULL)

  # Passed all criteria - return result
  create_result_row(kk, covs.in, nx, event_counts, cox_result)
}


#' Check if Matrix Has Positive Variance
#'
#' @keywords internal
has_positive_variance <- function(x) {
  xpx <- t(x) %*% x
  return(all(xpx > 0))
}


#' Check Prevalence Threshold
#'
#' @keywords internal
meets_prevalence_threshold <- function(x, minp) {
  return(all(colMeans(x) >= minp))
}


#' Calculate Event Counts by Treatment Arm
#'
#' @keywords internal
calculate_event_counts <- function(dd, tt, id.x) {
  list(
    d0 = sum(dd[id.x == 1 & tt == 0]),
    d1 = sum(dd[id.x == 1 & tt == 1]),
    total = sum(dd[id.x == 1])
  )
}


#' Check Event Count Criteria
#'
#' @keywords internal
meets_event_criteria <- function(event_counts, d0.min, d1.min) {
  return(event_counts$d0 >= d0.min && event_counts$d1 >= d1.min)
}


#' Fit Cox Model for Subgroup
#'
#' @keywords internal
#' @importFrom survival coxph Surv survfit
fit_cox_for_subgroup <- function(yy, dd, tt, id.x) {

  # Create subgroup data
  data.x <- data.table::data.table(Y = yy, E = dd, Treat = tt, id.x = id.x)
  df.x <- data.x[id.x == 1]

  # Fit Cox model
  hr.cox <- try(
    summary(coxph(Surv(Y, E) ~ Treat, data = df.x, robust = FALSE))$conf.int,
    silent = TRUE
  )

  if (inherits(hr.cox, "try-error")) return(NULL)

  # Get median survival times
  meds <- try(
    summary(survfit(Surv(Y, E) ~ Treat, data = df.x))$table[, "median"],
    silent = TRUE
  )

  if (inherits(meds, "try-error")) return(NULL)

  list(
    hr = hr.cox[1],
    lower = hr.cox[3],
    upper = hr.cox[4],
    med0 = meds[1],
    med1 = meds[2]
  )
}


#' Create Result Row
#'
#' @keywords internal
create_result_row <- function(kk, covs.in, nx, event_counts, cox_result) {
  c(
    kk,                           # group id
    sum(covs.in),                 # K (number of factors)
    nx,                           # n (subgroup size)
    event_counts$total,           # E (total events)
    event_counts$d1,              # d1 (events in treatment)
    cox_result$med1,              # m1 (median survival treatment)
    cox_result$med0,              # m0 (median survival control)
    cox_result$hr,                # HR
    cox_result$lower,             # Lower CI
    cox_result$upper,             # Upper CI
    covs.in                       # Factor indicators
  )
}


#' Format Search Results
#'
#' @keywords internal
#' @importFrom data.table data.table setorder
format_search_results <- function(results_list, Z, details, t.sofar, L, max_count) {

  # No results found
  if (length(results_list) == 0) {
    if (details) cat("NO subgroup candidates found\n")
    return(NULL)
  }

  # Combine results into matrix
  HR.model.k <- do.call(rbind, results_list)

  # Convert to data.table with proper names
  hr.out <- data.table::data.table(HR.model.k)
  names(hr.out) <- c("grp", "K", "n", "E", "d1", "m1", "m0",
                     "HR", "L(HR)", "U(HR)", colnames(Z))
  rownames(hr.out) <- NULL

  # Sort by HR (descending) and K (ascending)
  hr.out <- data.table::setorder(hr.out, -HR, K)

  if (details) {
    cat("Found", nrow(hr.out), "subgroup candidate(s)\n")
  }

  # Return formatted output
  list(
    out.found = list(hr.subgroups = hr.out),
    max_sg_est = max(hr.out$HR),
    time_search = t.sofar,
    L = L,
    max_count = max_count,
    prop_max_count = 1
  )
}


