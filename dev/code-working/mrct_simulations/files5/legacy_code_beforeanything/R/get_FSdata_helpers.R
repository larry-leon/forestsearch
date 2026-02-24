
#' Check if a variable is continuous
#'
#' Determines if a variable is continuous based on the number of unique values.
#'
#' @param x A vector.
#' @param cutoff Integer. Minimum number of unique values to be considered continuous.
#' @return 1 if continuous, 2 if not.
#' @keywords internal

is.continuous <- function(x,cutoff = 4){ifelse(length(unique(x))>=cutoff,1,2)}

#' 25th Percentile (Quantile Low)
#'
#' Returns the 25th percentile of a numeric vector.
#'
#' @param x A numeric vector.
#' @return Numeric value of the 25th percentile.
#' @importFrom stats quantile
#' @keywords internal

qlow <- function(x) c(quantile(x,0.25))


#' 75th Percentile (Quantile High)
#'
#' Returns the 75th percentile of a numeric vector.
#'
#' @param x A numeric vector.
#' @return Numeric value of the 75th percentile.
#' @importFrom stats quantile
#' @keywords internal

qhigh <- function(x) c(quantile(x,0.75))

# For continuous variables in conf_force_names
# setup for mean, median, qlow, and qhigh

#' Generate cut expressions for a variable
#'
#' For a continuous variable, returns expressions for mean, median, qlow, and qhigh cuts.
#'
#' @param x Character. Variable name.
#' @return Character vector of cut expressions.
#' @keywords internal

cut_var <- function(x){
  mx <- paste0("mean(",x,")")
  a <- paste0(x," <= ",mx)
  mdx <- paste0("median(",x,")")
  b <- paste0(x," <= ",mdx)
  qlx <- paste0("qlow(",x,")")
  c <- paste0(x," <= ",qlx)
  qhx <- paste0("qhigh(",x,")")
  d <- paste0(x," <= ",qhx)
  return(c(a,b,c,d))
}

#' Get forced cut expressions for variables
#'
#' For each variable in \code{conf.force.names}, returns cut expressions if continuous.
#'
#' @param df Data frame.
#' @param conf.force.names Character vector of variable names.
#' @param cont.cutoff Integer. Cutoff for continuous.
#' @return Character vector of cut expressions.
#' @export

get_conf_force <- function(df, conf.force.names, cont.cutoff = 4) {
  # Validate input
  if (!is.data.frame(df)) stop("df must be a data.frame.")
  if (!is.character(conf.force.names)) stop("conf.force.names must be a character vector.")
  res <- list()
  for (name in conf.force.names) {
    if (!name %in% names(df)) {
      warning(paste("Variable", name, "not found in data frame. Skipping."))
      next
    }
    var_data <- df[[name]]
    # Check if variable is continuous
    flag_cont <- is.continuous(var_data, cutoff = cont.cutoff)
    if (flag_cont == 1) {
      # Create mean, median, qlow, and qhigh cuts
      cuts <- cut_var(x = name)
      res[[name]] <- cuts
    } else {
      res[[name]] <- NULL
    }
  }
  # Flatten to character vector if needed
  unlist(res)
}

#' LASSO selection for Cox model
#'
#' Performs LASSO variable selection using Cox regression.
#'
#' @param df Data frame.
#' @param confounders.name Character vector of confounder names.
#' @param outcome.name Character. Name of outcome variable.
#' @param event.name Character. Name of event indicator variable.
#' @param seedit Integer. Random seed.
#' @return List with selected, omitted variables, coefficients, lambda, and fits.
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom survival Surv
#' @export

lasso_selection <- function(df, confounders.name, outcome.name, event.name, seedit = 8316951) {
set.seed(seedit)
# Package checks
if (!requireNamespace("glmnet", quietly = TRUE)) stop("Package 'glmnet' is required.")
if (!requireNamespace("survival", quietly = TRUE)) stop("Package 'survival' is required.")

  # Prepare the design matrix and response
  x <- as.matrix(df[, confounders.name, drop = FALSE])

  y <- survival::Surv(df[[outcome.name]], df[[event.name]])

  # Fit Cox LASSO with cross-validation
  cvfit <- glmnet::cv.glmnet(x, y, family = "cox", alpha = 1)
  lambda_min <- cvfit$lambda.min
  fit <- glmnet::glmnet(x, y, family = "cox", alpha = 1, lambda = lambda_min)

  # Extract coefficients at lambda.min
  coefs <- as.vector(coef(fit))
  names(coefs) <- rownames(coef(fit))

  # Determine selected and omitted variables
  selected <- names(coefs)[coefs != 0]
  omitted <- names(coefs)[coefs == 0]

  # Return results as a list
  list(
    selected = selected,
    omitted = omitted,
    coefficients = coefs,
    lambda_min = lambda_min,
    cvfit = cvfit,
    fit = fit
  )
}


#' Filter a vector by LASSO-selected variables
#'
#' Returns elements of \code{x} that are in \code{lassokeep}.
#'
#' @param x Character vector.
#' @param lassokeep Character vector of selected variables.
#' @return Filtered character vector or NULL.
#' @export

filter_by_lassokeep <- function(x, lassokeep) {
  if (!is.null(x)) {
    filtered <- x[x %in% lassokeep]
    if (length(filtered) > 0) {
      return(filtered)
    } else {
      return(NULL)
    }
  }
  return(NULL)
}


#' Process forced cut expression for a variable
#'
#' Evaluates a cut expression (e.g., "age <= mean(age)") and returns the expression with the value.
#'
#' @param expr Character string of the cut expression.
#' @param df Data frame.
#' @return Character string with evaluated value.
#' @export

process_conf_force_expr <- function(expr, df) {
  # Match pattern: variable <= function(variable)
  # Examples: "age <= mean(age)", "size <= qlow(size)"
  pattern <- "^\\s*([a-zA-Z0-9_.]+)\\s*<=\\s*([a-zA-Z]+)\\(([^)]+)\\)\\s*$"
  m <- regexec(pattern, expr)
  matches <- regmatches(expr, m)[[1]]
  if (length(matches) == 0) {
    # If not matching, return as is
    return(expr)
  }
  var <- matches[2]
  fun <- matches[3]
  arg <- matches[4]
  # Only proceed if var and arg match
  if (!(var %in% colnames(df)) || !(arg %in% colnames(df))) return(expr)
  # Evaluate the function
  if (fun == "mean") {
    val <- round(mean(df[[arg]], na.rm = TRUE), 1)
  } else if (fun == "median") {
    val <- round(median(df[[arg]], na.rm = TRUE), 1)
  } else if (fun == "qlow") {
    val <- round(quantile(df[[arg]], 0.25, na.rm = TRUE), 1)
  } else if (fun == "qhigh") {
    val <- round(quantile(df[[arg]], 0.75, na.rm = TRUE), 1)
  } else {
    # Unknown function, return as is
    return(expr)
  }
  # Return the evaluated expression
  paste0(var, " <= ", val)
}


#' Get variable name from cut expression
#'
#' Extracts the variable name from a cut expression.
#'
#' @param thiscut Character string of the cut expression.
#' @param confounders.name Character vector of confounder names.
#' @return Character vector of variable names.
#' @keywords internal

get_cut_name <- function(thiscut, confounders.name) {
  cov_index <- vapply(confounders.name, function(x) grepl(x, thiscut), logical(1))
  confounders.name[cov_index]
}

#' Check if cut expression is for a continuous variable (OPTIMIZED)
#'
#' Determines if a cut expression refers to a continuous variable.
#' This optimized version avoids redundant lookups by using word boundary
#' matching instead of partial string matching.
#'
#' @param thiscut Character string of the cut expression.
#' @param confounders.name Character vector of confounder names.
#' @param df Data frame.
#' @param cont.cutoff Integer. Cutoff for continuous.
#'
#' @return Logical; TRUE if continuous, FALSE otherwise.
#' @keywords internal

is_flag_continuous <- function(thiscut, confounders.name, df, cont.cutoff) {
  # Use word boundaries to avoid partial matches
  # e.g., "z1" won't accidentally match "z11"
  for (conf in confounders.name) {
    # Pattern: word boundary + variable name + word boundary
    pattern <- paste0("\\b", conf, "\\b")
    if (grepl(pattern, thiscut)) {
      # Found the confounder, now check if it's continuous
      return(is.continuous(df[[conf]], cutoff = cont.cutoff) == 1)
    }
  }
  FALSE
}


#' Check if cut expression should be dropped
#'
#' Determines if a cut expression should be dropped (e.g., variable has <=1 unique value).
#'
#' @param thiscut Character string of the cut expression.
#' @param confounders.name Character vector of confounder names.
#' @param df Data frame.
#' @return Logical; TRUE if should be dropped, FALSE otherwise.
#' @keywords internal

is_flag_drop <- function(thiscut, confounders.name, df) {
  cut_name <- get_cut_name(thiscut, confounders.name)
  aa <- df[, cut_name, drop = FALSE]
  # If multiple columns, check all; if any has <=1 unique value, return TRUE
  any(vapply(aa, function(col) length(unique(col)) <= 1, logical(1)))
}


#' Disjunctive (dummy) coding for factor columns
#'
#' @param df Data frame with factor variables.
#' @return Data frame with dummy-coded columns.
#' @keywords internal
acm.disjctif <- function(df) {
  encode_col <- function(i) {
    cl <- as.factor(df[, i])
    cha <- colnames(df)[i]
    n <- length(cl)
    x <- matrix(0L, n, nlevels(cl))
    x[cbind(seq_len(n), as.integer(cl))] <- 1L
    dimnames(x) <- list(rownames(df), paste(cha, levels(cl), sep = "."))
    x
  }
  parts <- lapply(seq_len(ncol(df)), encode_col)
  data.frame(do.call(cbind, parts), check.names = FALSE)
  }


#' Dummy-code a data frame (numeric pass-through, factors expanded)
#'
#' @param df Data frame with numeric and/or factor columns.
#' @return Data frame with numeric columns unchanged and factor columns
#'   expanded via \code{\link{acm.disjctif}}.
#' @export
dummy_encode <- function(df) {
  stopifnot(is.data.frame(df))
  is_num <- vapply(df, is.numeric, logical(1))
  is_fac <- vapply(df, is.factor, logical(1))
  parts <- list()
  if (any(is_num)) parts[[length(parts) + 1L]] <- df[, is_num, drop = FALSE]
  if (any(is_fac)) parts[[length(parts) + 1L]] <- acm.disjctif(df[, is_fac, drop = FALSE])
  if (length(parts) == 0L) stop("df contains no numeric or factor columns")
  do.call(data.frame, c(parts, list(check.names = FALSE)))
}

#' @rdname dummy_encode
#' @noRd
dummy <- dummy_encode

#' @rdname dummy_encode
#' @noRd
dummy2 <- dummy_encode



#' Trailing zeros in binary representation
#'
#' Returns the number of trailing zeros in the binary representation of an integer.
#'
#' @param kk Integer.
#' @return Integer count of trailing zeros.
#' @noRd

ztrail <- function(kk){
  ii <- 1
  zz <- kk
  while(zz%%2 == 0){
    ii <- ii+1
    zz <- zz/2
  }
  return(ii)
}


#' Flip binary value(s)
#'
#' @param x Integer vector of 0s and 1s.
#' @return Integer vector with values flipped.
#' @noRd
one.zero <- function(x) 1L - x


#' Cache and validate cut expressions efficiently
#'
#' Evaluates all cut expressions once and caches results to avoid
#' redundant evaluation. Much faster than evaluating repeatedly.
#'
#' @param confs Character vector of cut expressions.
#' @param df Data frame to evaluate expressions against.
#' @param details Logical. Print details during execution.
#'
#' @return List with:
#'   - evaluations: List of evaluated vectors (logical TRUE/FALSE) for each cut
#'   - is_valid: Logical vector indicating which cuts produced >1 unique value
#'   - has_error: Logical vector indicating which cuts failed to evaluate
#'
#' @details
#' This replaces multiple eval(parse()) calls scattered throughout get_FSdata.
#' By caching results, we avoid:
#' 1. Repeated parsing of expressions
#' 2. Repeated evaluation on dataframe
#' 3. Redundant uniqueness checks
#'
#' @export

evaluate_cuts_once <- function(confs, df, details = FALSE) {
  n_confs <- length(confs)
  evaluations <- vector("list", n_confs)
  is_valid <- logical(n_confs)
  has_error <- logical(n_confs)

  for (i in seq_along(confs)) {
    thiscut <- confs[i]

    tryCatch({
      # Use evaluate_comparison() — no eval(parse()) needed
      # Cut expressions are always single comparisons like "er <= 0"
      result <- evaluate_comparison(thiscut, df)

      evaluations[[i]] <- as.logical(result)
      is_valid[i] <- length(unique(result)) > 1

    }, error = function(e) {
      has_error[i] <<- TRUE
      is_valid[i] <<- FALSE
      if (details) {
        cat("Error evaluating cut '", thiscut, "': ", e$message, "\n", sep = "")
      }
    })
  }

  if (details) {
    cat("Cut evaluation summary:\n")
    cat("  Total cuts: ", n_confs, "\n")
    cat("  Valid cuts: ", sum(is_valid), "\n")
    cat("  Errors: ", sum(has_error), "\n")
  }

  list(
    evaluations = evaluations,
    is_valid = is_valid,
    has_error = has_error
  )
}



