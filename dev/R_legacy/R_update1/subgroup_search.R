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


#' Get subgroup membership vector
#'
#' Returns a vector indicating subgroup membership (1 if all selected factors are present, 0 otherwise).
#'
#' @param zz Matrix or data frame of subgroup factor indicators.
#' @param covs.in Numeric vector indicating which factors are selected (1 = included).
#' @return Numeric vector of subgroup membership (1/0).
#' @export

get_subgroup_membership <- function(zz, covs.in) {
  x <- zz[, which(covs.in == 1), drop = FALSE]
  if (ncol(x) == 0) return(rep(1, nrow(zz)))
  apply(x, 1, function(row) all(row == 1))
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

subgroup.search <- function(Y,Event,Treat,ID = NULL,Z,n.min = 30,d0.min = 15,d1.min = 15,
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
  results_list <- lapply(1:tot_counts, function(kk) {
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
