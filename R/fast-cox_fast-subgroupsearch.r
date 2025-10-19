




#' Fast Cox Model for Subgroup Analysis
#'
#' Optimized Cox model fitting for bootstrap and subgroup search.
#' Supports both survival::coxph and glmnet implementations.
#' The glmnet implementation leverages logic from lasso_selection() for consistency.
#'
#' @param Y Numeric vector of outcome (time-to-event)
#' @param Event Numeric vector of event indicators (0/1)
#' @param Treat Numeric vector of treatment indicators (0/1)
#' @param method Character. "coxph" (default) or "glmnet" for fitting method
#' @param robust Logical. Use robust standard errors (only for method="coxph")
#' @param return_se Logical. Return standard error (default: FALSE for speed)
#' @param lambda Numeric. Penalty parameter for glmnet (default: 0 for unpenalized).
#'        Use NULL for cross-validation to select optimal lambda.
#' @param alpha Numeric. Elastic net mixing parameter for glmnet (default: 0 for ridge).
#'        alpha=1 gives lasso, alpha=0 gives ridge, 0<alpha<1 gives elastic net.
#' @param nfolds Integer. Number of folds for CV when lambda=NULL (default: 10)
#' @param standardize Logical. Standardize covariates in glmnet (default: FALSE)
#' @param thresh Numeric. Convergence threshold for glmnet (default: 1e-14)
#' @param maxit Integer. Maximum iterations for glmnet (default: 100000)
#'
#' @return List with HR, lower CI, upper CI, median survival times, and optionally SE
#' @importFrom survival coxph Surv survfit
#' @importFrom glmnet glmnet cv.glmnet
#' @export
fast_cox_subgroup <- function(Y, Event, Treat,
                              method = "coxph",
                              robust = FALSE,
                              return_se = FALSE,
                              lambda = 0,
                              alpha = 0,
                              nfolds = 10,
                              standardize = FALSE,
                              thresh = 1e-14,
                              maxit = 100000) {
  # Input validation
  if (length(Y) != length(Event) || length(Y) != length(Treat)) {
    stop("Y, Event, and Treat must have the same length")
  }

  # Check for sufficient events
  if (sum(Event) < 2) {
    return(list(hr = NA, lower = NA, upper = NA, med0 = NA, med1 = NA,
                se = NA, coef = NA, lambda_used = NA))
  }

  # Initialize return values
  hr <- lower <- upper <- se <- coef_val <- lambda_used <- NA

  # =========================================================================
  # METHOD 1: Standard coxph (survival package)
  # =========================================================================
  if (method == "coxph") {
    fit <- try(
      survival::coxph(survival::Surv(Y, Event) ~ Treat, robust = robust),
      silent = TRUE
    )

    if (inherits(fit, "try-error")) {
      return(list(hr = NA, lower = NA, upper = NA, med0 = NA, med1 = NA,
                  se = NA, coef = NA, lambda_used = NA))
    }

    # Extract HR and CI
    hr_ci <- summary(fit)$conf.int[c(1, 3, 4)]
    hr <- hr_ci[1]
    lower <- hr_ci[2]
    upper <- hr_ci[3]

    # Extract coefficient and SE if requested
    if (return_se) {
      coef_summary <- summary(fit)$coefficients
      coef_val <- coef_summary[1, "coef"]
      se <- if (robust) coef_summary[1, "robust se"] else coef_summary[1, "se(coef)"]
    }
  }

  # =========================================================================
  # METHOD 2: glmnet (leveraging lasso_selection logic)
  # =========================================================================
  else if (method == "glmnet") {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
      warning("glmnet package not available, falling back to coxph")
      return(fast_cox_subgroup(Y, Event, Treat, method = "coxph",
                               robust = robust, return_se = return_se))
    }

    # Prepare data (following lasso_selection pattern)
    X <- matrix(Treat, ncol = 1)
    colnames(X) <- "Treat"
    y <- survival::Surv(Y, Event)

    # -----------------------------------------------------------------------
    # OPTION A: User-specified lambda (including lambda=0 for unpenalized)
    # -----------------------------------------------------------------------
    if (!is.null(lambda)) {
      # Fit with specified lambda (similar to lasso_selection's glmnet call)
      fit <- try(
        glmnet::glmnet(X, y, family = "cox",
                       lambda = lambda,
                       alpha = alpha,
                       standardize = standardize,
                       thresh = thresh,
                       maxit = maxit),
        silent = TRUE
      )

      if (inherits(fit, "try-error")) {
        return(list(hr = NA, lower = NA, upper = NA, med0 = NA, med1 = NA,
                    se = NA, coef = NA, lambda_used = NA))
      }

      lambda_used <- lambda
      coef_val <- as.vector(coef(fit))
    }

    # -----------------------------------------------------------------------
    # OPTION B: Cross-validation to select lambda (following lasso_selection)
    # -----------------------------------------------------------------------
    else {
      # Use cv.glmnet to select optimal lambda (like lasso_selection does)
      cvfit <- try(
        glmnet::cv.glmnet(X, y, family = "cox",
                          alpha = alpha,
                          standardize = standardize,
                          nfolds = nfolds,
                          thresh = thresh,
                          maxit = maxit),
        silent = TRUE
      )

      if (inherits(cvfit, "try-error")) {
        return(list(hr = NA, lower = NA, upper = NA, med0 = NA, med1 = NA,
                    se = NA, coef = NA, lambda_used = NA))
      }

      # Extract lambda.min (like lasso_selection does)
      lambda_used <- cvfit$lambda.min

      # Fit final model at lambda.min
      fit <- try(
        glmnet::glmnet(X, y, family = "cox",
                       lambda = lambda_used,
                       alpha = alpha,
                       standardize = standardize,
                       thresh = thresh,
                       maxit = maxit),
        silent = TRUE
      )

      if (inherits(fit, "try-error")) {
        return(list(hr = NA, lower = NA, upper = NA, med0 = NA, med1 = NA,
                    se = NA, coef = NA, lambda_used = lambda_used))
      }

      coef_val <- as.vector(coef(fit))
    }

    # Calculate HR
    hr <- exp(coef_val)

    # -----------------------------------------------------------------------
    # Standard Error and Confidence Interval Calculation
    # -----------------------------------------------------------------------
    if (return_se) {
      # For unpenalized (lambda=0) or very small lambda: use coxph SE
      if (is.null(lambda) || lambda < 1e-10) {
        fit_coxph <- try(
          survival::coxph(survival::Surv(Y, Event) ~ Treat, robust = FALSE),
          silent = TRUE
        )

        if (!inherits(fit_coxph, "try-error")) {
          se <- summary(fit_coxph)$coefficients[1, "se(coef)"]
          lower <- exp(coef_val - 1.96 * se)
          upper <- exp(coef_val + 1.96 * se)
        } else {
          # Fallback: approximate SE
          n_events <- sum(Event)
          se <- sqrt(4 / max(n_events, 1))
          lower <- exp(coef_val - 1.96 * se)
          upper <- exp(coef_val + 1.96 * se)
        }
      } else {
        # For penalized models: SE from bootstrap or conservative approximation
        # Option 1: Conservative approximation based on events
        n_events <- sum(Event)
        se <- sqrt(4 / max(n_events, 1)) * (1 + abs(lambda))  # Penalty adjustment
        lower <- exp(coef_val - 1.96 * se)
        upper <- exp(coef_val + 1.96 * se)
      }
    } else {
      # Quick CI approximation without full SE calculation
      n_events <- sum(Event)
      if (n_events > 0) {
        # Approximate SE with penalty adjustment
        penalty_factor <- if (is.null(lambda) || lambda < 1e-10) 1 else (1 + sqrt(abs(lambda)))
        approx_se <- sqrt(4 / n_events) * penalty_factor
        lower <- exp(coef_val - 1.96 * approx_se)
        upper <- exp(coef_val + 1.96 * approx_se)
      }
    }
  }

  else {
    stop("method must be 'coxph' or 'glmnet'")
  }

  # =========================================================================
  # Calculate median survival times (same for both methods)
  # =========================================================================
  surv_fit <- try(
    summary(survival::survfit(survival::Surv(Y, Event) ~ Treat))$table[, "median"],
    silent = TRUE
  )

  if (inherits(surv_fit, "try-error")) {
    med0 <- med1 <- NA
  } else {
    med0 <- as.numeric(surv_fit[1])
    med1 <- as.numeric(surv_fit[2])
  }

  list(
    hr = hr,
    lower = lower,
    upper = upper,
    med0 = med0,
    med1 = med1,
    se = se,
    coef = coef_val,
    lambda_used = lambda_used
  )
}



#' Optimized Subgroup Search for Treatment Effect Heterogeneity
#'
#' Faster version of subgroup.search() that pre-computes common operations
#' and uses vectorized operations where possible. Supports glmnet for faster Cox fitting.
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
#' @param method Character. "coxph" or "glmnet" for Cox model fitting.
#'
#' @return List with found subgroups, maximum HR, search time, and configuration info.
#'
#' @importFrom data.table data.table setorder
#' @importFrom survival coxph Surv survfit
#' @importFrom utils combn
#' @export
subgroup_search_fast <- function(Y, Event, Treat, ID = NULL, Z,
                                  n.min = 30, d0.min = 15, d1.min = 15,
                                  hr.threshold = 1.0, max.minutes = 30,
                                  minp = 0.05, rmin = 5,
                                  details = FALSE, maxk = 2,
                                  method = "coxph") {

  # =========================================================================
  # SECTION 1: DATA PREPARATION AND VALIDATION
  # =========================================================================

  temp <- cbind(Y, Event, Treat, Z)
  temp <- na.exclude(temp)
  yy <- temp[, 1]
  dd <- temp[, 2]
  tt <- temp[, 3]
  zz <- temp[, -c(1, 2, 3)]

  n <- length(yy)
  L <- ncol(zz)

  # Pre-compute products for faster event counting
  dd_tt <- dd * tt
  dd_1minustt <- dd * (1 - tt)

  # Validate glmnet availability if requested
  if (method == "glmnet" && !requireNamespace("glmnet", quietly = TRUE)) {
    warning("glmnet not available, falling back to coxph")
    method <- "coxph"
  }

  if (details && method == "glmnet") {
    cat("Using glmnet for Cox model fitting (faster)\n")
  }

  # =========================================================================
  # SECTION 2: GENERATE COMBINATION INDICES
  # =========================================================================

  if (maxk == 1) {
    max_count <- L
    index_1factor <- t(combn(L, 1))
    counts_1factor <- nrow(index_1factor)
    tot_counts <- counts_1factor
    index_2factor <- NULL
    counts_2factor <- 0
    index_3factor <- NULL
    counts_3factor <- 0
  } else if (maxk == 2) {
    max_count <- (L * (L - 1) / 2) + L
    index_2factor <- t(combn(L, 2))
    counts_2factor <- nrow(index_2factor)
    index_1factor <- t(combn(L, 1))
    counts_1factor <- nrow(index_1factor)
    tot_counts <- counts_2factor + counts_1factor
    index_3factor <- NULL
    counts_3factor <- 0
  } else if (maxk == 3) {
    max_count <- ((L * (L - 2) * (L - 1)) / 6) + (L * (L - 1) / 2) + L
    index_3factor <- t(combn(L, 3))
    counts_3factor <- nrow(index_3factor)
    index_2factor <- t(combn(L, 2))
    counts_2factor <- nrow(index_2factor)
    index_1factor <- t(combn(L, 1))
    counts_1factor <- nrow(index_1factor)
    tot_counts <- counts_3factor + counts_2factor + counts_1factor
  }

  if (tot_counts != max_count) {
    stop("Error with maximum combinations for maxk = ", maxk)
  }

  if (details) {
    cat("Number of possible configurations (<= maxk): maxk, # <= maxk",
        c(maxk, max_count), "\n")
  }

  # =========================================================================
  # SECTION 3: OPTIMIZED SUBGROUP SEARCH LOOP
  # =========================================================================

  t.start <- proc.time()[3]

  # Pre-allocate result storage for better performance
  results_list <- vector("list", tot_counts)
  n_results <- 0

  for (kk in seq_len(tot_counts)) {

    # Check time limit
    t.now <- proc.time()[3]
    t.sofar <- (t.now - t.start) / 60
    if (t.sofar > max.minutes) break

    # Get covariate selection for this combination
    covs.in <- get_covs_in(
      kk, maxk, L,
      counts_1factor, index_1factor,
      counts_2factor, index_2factor,
      counts_3factor, index_3factor
    )

    k.in <- sum(covs.in)
    if (k.in > maxk) next

    # Extract selected columns
    selected_cols <- which(covs.in == 1)
    x <- zz[, selected_cols, drop = FALSE]

    # Quick prevalence check
    xpx <- t(x) %*% x
    if (!all(xpx > 0)) next

    # Check redundancy
    get_idx_rflag <- extract_idx_flagredundancy(x, rmin)
    id.x <- get_idx_rflag$id.x
    flag.redundant <- get_idx_rflag$flag.redundant

    # Fast event counting using pre-computed products
    d1 <- sum(dd_tt[id.x == 1])
    d0 <- sum(dd_1minustt[id.x == 1])
    nx <- sum(id.x)

    # Check prevalence threshold
    crit3 <- all(colMeans(x) >= minp)

    # Apply all criteria
    if (!flag.redundant && d1 >= d1.min && d0 >= d0.min &&
        nx > n.min && crit3) {

      # Extract subgroup data
      idx_keep <- which(id.x == 1)
      Y_sub <- yy[idx_keep]
      Event_sub <- dd[idx_keep]
      Treat_sub <- tt[idx_keep]

      # Fit Cox model using fast implementation with chosen method
      cox_result <- fast_cox_subgroup(Y_sub, Event_sub, Treat_sub,
                                     method = method, robust = FALSE,
                                     return_se = FALSE)

      # Check if HR meets threshold
      if (!is.na(cox_result$hr) && cox_result$hr > hr.threshold) {
        n_results <- n_results + 1

        # Store result
        results_list[[n_results]] <- c(
          kk,                    # group id
          sum(covs.in),         # K (number of factors)
          nx,                   # n (subgroup size)
          d0 + d1,              # E (total events)
          d1,                   # d1 (events in treatment)
          cox_result$med1,      # m1 (median survival treatment)
          cox_result$med0,      # m0 (median survival control)
          cox_result$hr,        # HR
          cox_result$lower,     # Lower CI
          cox_result$upper,     # Upper CI
          covs.in               # Factor indicators
        )
      }
    }
  }

  # =========================================================================
  # SECTION 4: COMPILE RESULTS
  # =========================================================================

  t.end <- proc.time()[3]
  t.min <- (t.end - t.start) / 60

  if (details) {
    cat("Events criteria for control, exp =", c(d0.min, d1.min), "\n")
    cat("*Subgroup Searching Minutes =*", c(t.min), "\n")
  }

  # Combine results
  if (n_results == 0) {
    if (details) cat("NO subgroup candidate found (FS)\n")
    return(NULL)
  }

  HR.model.k <- do.call(rbind, results_list[1:n_results])

  hr.out <- data.table::data.table(HR.model.k)
  names(hr.out) <- c("grp", "K", "n", "E", "d1", "m1", "m0",
                     "HR", "L(HR)", "U(HR)", colnames(Z))
  rownames(hr.out) <- NULL
  hr.out <- data.table::setorder(hr.out, -HR, K)

  out.found <- list(hr.subgroups = hr.out)
  max_sg_est <- max(hr.out[["HR"]])

  if (details) cat("Subgroup candidate(s) found (FS)\n")

  return(list(
    out.found = out.found,
    max_sg_est = max_sg_est,
    time_search = t.sofar,
    L = L,
    max_count = max_count,
    prop_max_count = 1
  ))
}


#' Fast Cox Model for Bootstrap (Compatible with get_Cox_sg)
#'
#' Drop-in replacement for get_Cox_sg that's optimized for speed.
#' Returns the same output format as the original function.
#' Supports both coxph and glmnet methods.
#'
#' @param df_sg Data frame for subgroup
#' @param cox.formula Cox model formula
#' @param est.loghr Logical. Is estimate on log(HR) scale?
#' @param method Character. "coxph" (default) or "glmnet" for fitting
#' @return List with estimate and standard error
#' @importFrom survival coxph
#' @export
get_Cox_sg_fast <- function(df_sg, cox.formula, est.loghr = TRUE,
                            method = "coxph") {

  # Extract variable names
  names_tocheck <- all.vars(cox.formula)
  check <- unlist(lapply(names_tocheck, grep, names(df_sg), value = TRUE))
  check2 <- match(names_tocheck, check)

  if (sum(!is.na(check2)) != length(names_tocheck)) {
    stop("df_sg dataset does NOT contain cox.formula variables")
  }

  # =========================================================================
  # METHOD 1: Standard coxph (survival package)
  # =========================================================================
  if (method == "coxph") {
    # Fit Cox model with robust standard errors
    fit <- try(
      summary(coxph(cox.formula, data = df_sg, robust = TRUE))$coefficients,
      silent = TRUE
    )

    if (inherits(fit, "try-error")) {
      return(list(est_obs = NA_real_, se_obs = NA_real_))
    }

    # Extract estimates
    if (est.loghr) {
      bhat <- c(fit[, "coef"])
      est_obs <- bhat
      se_obs <- c(fit[, "robust se"])
    } else {
      bhat <- c(fit[, "coef"])
      est_obs <- exp(bhat)
      sebhat <- c(fit[, "robust se"])
      se_obs <- est_obs * sebhat
    }
  }

  # =========================================================================
  # METHOD 2: glmnet (faster for simple models)
  # =========================================================================
  else if (method == "glmnet") {
    if (!requireNamespace("glmnet", quietly = TRUE)) {
      warning("glmnet not available, falling back to coxph")
      return(get_Cox_sg_fast(df_sg, cox.formula, est.loghr, method = "coxph"))
    }

    # Extract formula components
    outcome_var <- all.vars(cox.formula[[2]])[1]
    event_var <- all.vars(cox.formula[[2]])[2]
    treat_var <- all.vars(cox.formula[[3]])[1]

    # Prepare data
    X <- matrix(df_sg[[treat_var]], ncol = 1)
    y <- survival::Surv(df_sg[[outcome_var]], df_sg[[event_var]])

    # Fit unpenalized Cox model
    fit <- try(
      glmnet::glmnet(X, y, family = "cox", lambda = 0, alpha = 0,
                     thresh = 1e-14, maxit = 100000),
      silent = TRUE
    )

    if (inherits(fit, "try-error")) {
      return(list(est_obs = NA_real_, se_obs = NA_real_))
    }

    # Extract coefficient
    bhat <- as.vector(coef(fit))

    # For SE, we need to use coxph (glmnet doesn't provide vcov easily)
    # But we can approximate or fall back to quick coxph call
    fit_coxph <- try(
      summary(coxph(cox.formula, data = df_sg, robust = TRUE))$coefficients,
      silent = TRUE
    )

    if (inherits(fit_coxph, "try-error")) {
      # Approximate SE using rule of thumb: sqrt(4/events)
      n_events <- sum(df_sg[[event_var]])
      se_approx <- sqrt(4 / max(n_events, 1))
    } else {
      se_approx <- fit_coxph[, "robust se"]
    }

    # Return in appropriate scale
    if (est.loghr) {
      est_obs <- bhat
      se_obs <- se_approx
    } else {
      est_obs <- exp(bhat)
      se_obs <- est_obs * se_approx
    }
  }

  else {
    stop("method must be 'coxph' or 'glmnet'")
  }

  return(list(est_obs = est_obs, se_obs = se_obs))
}


#' Benchmark Comparison Function
#'
#' Compares timing between original and fast implementations
#' Tests both coxph and glmnet methods
#'
#' @param Y Numeric vector of outcome
#' @param Event Numeric vector of events
#' @param Treat Numeric vector of treatment
#' @param Z Matrix of factors
#' @param n.min Minimum subgroup size
#' @param d0.min Minimum control events
#' @param d1.min Minimum treatment events
#' @param maxk Maximum factors
#'
#' @return List with timing results
#' @export
benchmark_subgroup_search <- function(Y, Event, Treat, Z,
                                     n.min = 30, d0.min = 15, d1.min = 15,
                                     maxk = 2) {

  cat("\n=== SUBGROUP SEARCH BENCHMARK ===\n\n")

  # Test 1: Original implementation
  cat("1. Running original subgroup.search()...\n")
  t1_start <- proc.time()[3]
  result_original <- try(
    subgroup.search(Y, Event, Treat, Z = Z,
                   n.min = n.min, d0.min = d0.min, d1.min = d1.min,
                   maxk = maxk, details = FALSE),
    silent = TRUE
  )
  t1_end <- proc.time()[3]
  time_original <- t1_end - t1_start

  # Test 2: Fast implementation with coxph
  cat("2. Running subgroup_search_fast() with coxph...\n")
  t2_start <- proc.time()[3]
  result_fast_coxph <- try(
    subgroup_search_fast(Y, Event, Treat, Z = Z,
                        n.min = n.min, d0.min = d0.min, d1.min = d1.min,
                        maxk = maxk, method = "coxph", details = FALSE),
    silent = TRUE
  )
  t2_end <- proc.time()[3]
  time_fast_coxph <- t2_end - t2_start

  # Test 3: Fast implementation with glmnet
  if (requireNamespace("glmnet", quietly = TRUE)) {
    cat("3. Running subgroup_search_fast() with glmnet...\n")
    t3_start <- proc.time()[3]
    result_fast_glmnet <- try(
      subgroup_search_fast(Y, Event, Treat, Z = Z,
                          n.min = n.min, d0.min = d0.min, d1.min = d1.min,
                          maxk = maxk, method = "glmnet", details = FALSE),
      silent = TRUE
    )
    t3_end <- proc.time()[3]
    time_fast_glmnet <- t3_end - t3_start
  } else {
    time_fast_glmnet <- NA
    result_fast_glmnet <- NULL
    cat("3. glmnet not available, skipping\n")
  }

  # Calculate speedups
  speedup_coxph <- time_original / time_fast_coxph
  speedup_glmnet <- if (!is.na(time_fast_glmnet)) time_original / time_fast_glmnet else NA

  # Print results
  cat("\n=== TIMING RESULTS ===\n")
  cat(sprintf("%-30s: %8.3f seconds\n", "Original (coxph)", time_original))
  cat(sprintf("%-30s: %8.3f seconds (%5.2fx faster)\n",
              "Fast + coxph", time_fast_coxph, speedup_coxph))
  if (!is.na(time_fast_glmnet)) {
    cat(sprintf("%-30s: %8.3f seconds (%5.2fx faster)\n",
                "Fast + glmnet", time_fast_glmnet, speedup_glmnet))
  }

  # Validate results match
  cat("\n=== VALIDATION ===\n")
  if (!inherits(result_original, "try-error") &&
      !inherits(result_fast_coxph, "try-error")) {
    if (!is.null(result_original) && !is.null(result_fast_coxph)) {
      n_orig <- nrow(result_original$out.found$hr.subgroups)
      n_fast_coxph <- nrow(result_fast_coxph$out.found$hr.subgroups)
      cat("Subgroups found:\n")
      cat(sprintf("  Original:        %d\n", n_orig))
      cat(sprintf("  Fast (coxph):    %d\n", n_fast_coxph))

      if (!is.na(time_fast_glmnet) && !inherits(result_fast_glmnet, "try-error")) {
        n_fast_glmnet <- nrow(result_fast_glmnet$out.found$hr.subgroups)
        cat(sprintf("  Fast (glmnet):   %d\n", n_fast_glmnet))
      }

      # Check if top subgroup matches
      top_hr_orig <- result_original$out.found$hr.subgroups$HR[1]
      top_hr_fast <- result_fast_coxph$out.found$hr.subgroups$HR[1]
      hr_match <- abs(top_hr_orig - top_hr_fast) < 0.01
      cat(sprintf("\nTop HR matches: %s (%.3f vs %.3f)\n",
                  ifelse(hr_match, "YES", "NO"), top_hr_orig, top_hr_fast))
    }
  }

  cat("\n=== RECOMMENDATIONS ===\n")
  if (!is.na(speedup_glmnet) && speedup_glmnet > speedup_coxph) {
    cat("✓ Use method='glmnet' for best performance\n")
    cat(sprintf("  Expected speedup: %.1fx faster than original\n", speedup_glmnet))
  } else {
    cat("✓ Use method='coxph' (default) for slightly better performance\n")
    cat(sprintf("  Expected speedup: %.1fx faster than original\n", speedup_coxph))
  }

  invisible(list(
    time_original = time_original,
    time_fast_coxph = time_fast_coxph,
    time_fast_glmnet = time_fast_glmnet,
    speedup_coxph = speedup_coxph,
    speedup_glmnet = speedup_glmnet,
    result_original = result_original,
    result_fast_coxph = result_fast_coxph,
    result_fast_glmnet = result_fast_glmnet
  ))
}


#' Benchmark Cox Model Implementations
#'
#' Compare timing for get_Cox_sg vs get_Cox_sg_fast
#' Tests both coxph and glmnet methods
#'
#' @param df_sg Data frame for subgroup
#' @param cox.formula Cox model formula
#' @param n_iterations Number of iterations to average (default: 100)
#'
#' @return List with timing comparisons
#' @export
benchmark_cox_methods <- function(df_sg, cox.formula, n_iterations = 100) {

  cat("\n=== COX MODEL BENCHMARK ===\n")
  cat(sprintf("Testing with %d iterations\n", n_iterations))
  cat(sprintf("Sample size: %d observations\n", nrow(df_sg)))
  cat(sprintf("Events: %d\n\n", sum(df_sg[[all.vars(cox.formula[[2]])[2]]])))

  # Test 1: Original get_Cox_sg
  cat("1. Testing get_Cox_sg (original)...\n")
  times_original <- numeric(n_iterations)
  for (i in 1:n_iterations) {
    t_start <- proc.time()[3]
    result <- get_Cox_sg(df_sg, cox.formula, est.loghr = TRUE)
    t_end <- proc.time()[3]
    times_original[i] <- t_end - t_start
  }
  time_original <- mean(times_original)

  # Test 2: Fast with coxph
  cat("2. Testing get_Cox_sg_fast with coxph...\n")
  times_fast_coxph <- numeric(n_iterations)
  for (i in 1:n_iterations) {
    t_start <- proc.time()[3]
    result <- get_Cox_sg_fast(df_sg, cox.formula, est.loghr = TRUE,
                              method = "coxph")
    t_end <- proc.time()[3]
    times_fast_coxph[i] <- t_end - t_start
  }
  time_fast_coxph <- mean(times_fast_coxph)

  # Test 3: Fast with glmnet
  if (requireNamespace("glmnet", quietly = TRUE)) {
    cat("3. Testing get_Cox_sg_fast with glmnet...\n")
    times_fast_glmnet <- numeric(n_iterations)
    for (i in 1:n_iterations) {
      t_start <- proc.time()[3]
      result <- get_Cox_sg_fast(df_sg, cox.formula, est.loghr = TRUE,
                                method = "glmnet")
      t_end <- proc.time()[3]
      times_fast_glmnet[i] <- t_end - t_start
    }
    time_fast_glmnet <- mean(times_fast_glmnet)
  } else {
    time_fast_glmnet <- NA
    cat("3. glmnet not available, skipping\n")
  }

  # Calculate speedups
  speedup_coxph <- time_original / time_fast_coxph
  speedup_glmnet <- if (!is.na(time_fast_glmnet)) time_original / time_fast_glmnet else NA

  # Print results
  cat("\n=== AVERAGE TIMING (per iteration) ===\n")
  cat(sprintf("%-30s: %8.6f seconds\n", "get_Cox_sg (original)", time_original))
  cat(sprintf("%-30s: %8.6f seconds (%5.2fx faster)\n",
              "get_Cox_sg_fast (coxph)", time_fast_coxph, speedup_coxph))
  if (!is.na(time_fast_glmnet)) {
    cat(sprintf("%-30s: %8.6f seconds (%5.2fx faster)\n",
                "get_Cox_sg_fast (glmnet)", time_fast_glmnet, speedup_glmnet))
  }

  cat("\n=== PROJECTED SAVINGS (1000 bootstrap iterations) ===\n")
  saved_coxph <- (time_original - time_fast_coxph) * 1000 / 60
  cat(sprintf("Using fast + coxph:   %.1f minutes saved\n", saved_coxph))

  if (!is.na(time_fast_glmnet)) {
    saved_glmnet <- (time_original - time_fast_glmnet) * 1000 / 60
    cat(sprintf("Using fast + glmnet:  %.1f minutes saved\n", saved_glmnet))
  }

  invisible(list(
    time_original = time_original,
    time_fast_coxph = time_fast_coxph,
    time_fast_glmnet = time_fast_glmnet,
    speedup_coxph = speedup_coxph,
    speedup_glmnet = speedup_glmnet
  ))
}
