
#' Check if a variable is continuous
#'
#' Determines if a variable is continuous based on the number of unique values.
#'
#' @param x A vector.
#' @param cutoff Integer. Minimum number of unique values to be considered continuous.
#' @return 1 if continuous, 2 if not.
#' @export

is.continuous <- function(x,cutoff = 4){ifelse(length(unique(x))>=cutoff,1,2)}

#' 25th Percentile (Quantile Low)
#'
#' Returns the 25th percentile of a numeric vector.
#'
#' @param x A numeric vector.
#' @return Numeric value of the 25th percentile.
#' @importFrom stats quantile
#' @export

qlow <- function(x) c(quantile(x,0.25))


#' 75th Percentile (Quantile High)
#'
#' Returns the 75th percentile of a numeric vector.
#'
#' @param x A numeric vector.
#' @return Numeric value of the 75th percentile.
#' @importFrom stats quantile
#' @export

qhigh <- function(x) c(quantile(x,0.75))

# For continuous variables in conf_force_names
# setup for mean, median, qlow, and qhigh

#' Generate cut expressions for a variable
#'
#' For a continuous variable, returns expressions for mean, median, qlow, and qhigh cuts.
#'
#' @param x Character. Variable name.
#' @return Character vector of cut expressions.
#' @export

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

#' Format subgroup labels for ForestSearch
#'
#' Converts q-indexed subgroup codes to human-readable labels.
#'
#' @param Qsg Character vector of codes (e.g., 'q1.0', 'q2.1').
#' @param confs_labels Character vector of confounder labels.
#' @return Character vector of formatted subgroup labels.
#' @importFrom stringr str_length str_sub
#' @export

FS_labels <- function(Qsg, confs_labels) {
  # Validate inputs
  if (!is.character(Qsg)) stop("Qsg must be a character vector.")
  if (!is.character(confs_labels)) stop("confs_labels must be a character vector.")

  # Use regex to extract index and action
  # Pattern: q<index>.<action>, e.g., q1.0, q12.1, q123.0, q1234.1
  pattern <- "^q(\\d+)\\.(\\d+)$"
  matches <- regexec(pattern, Qsg)
  parts <- regmatches(Qsg, matches)

  sg_labels <- character(length(Qsg))

  for (i in seq_along(parts)) {
    if (length(parts[[i]]) == 3) {
      idx <- as.numeric(parts[[i]][2])
      action <- as.numeric(parts[[i]][3])
      if (!is.na(idx) && idx >= 1 && idx <= length(confs_labels)) {
        label <- confs_labels[idx]
        if (action == 0) {
          sg_labels[i] <- paste0("!{", label, "}")
        } else if (action == 1) {
          sg_labels[i] <- paste0("{", label, "}")
        } else {
          sg_labels[i] <- NA
        }
      } else {
        sg_labels[i] <- NA
      }
    } else {
      sg_labels[i] <- NA
    }
  }
  return(sg_labels)
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
#' @export

get_cut_name <- function(thiscut, confounders.name) {
  cov_index <- vapply(confounders.name, function(x) grepl(x, thiscut), logical(1))
  confounders.name[cov_index]
}


#' Check if cut expression is for a continuous variable
#'
#' Determines if a cut expression refers to a continuous variable.
#'
#' @param thiscut Character string of the cut expression.
#' @param confounders.name Character vector of confounder names.
#' @param df Data frame.
#' @param cont.cutoff Integer. Cutoff for continuous.
#' @return Logical; TRUE if continuous, FALSE otherwise.
#' @export

is_flag_continuous <- function(thiscut, confounders.name, df, cont.cutoff) {
  cut_name <- get_cut_name(thiscut, confounders.name)
  aa <- df[, cut_name, drop = FALSE]
  # If multiple columns, check all; if any is continuous, return TRUE
  any(vapply(aa, function(col) is.continuous(col, cutoff = cont.cutoff) == 1, logical(1)))
}


#' Check if cut expression should be dropped
#'
#' Determines if a cut expression should be dropped (e.g., variable has <=1 unique value).
#'
#' @param thiscut Character string of the cut expression.
#' @param confounders.name Character vector of confounder names.
#' @param df Data frame.
#' @return Logical; TRUE if should be dropped, FALSE otherwise.
#' @export

is_flag_drop <- function(thiscut, confounders.name, df) {
  cut_name <- get_cut_name(thiscut, confounders.name)
  aa <- df[, cut_name, drop = FALSE]
  # If multiple columns, check all; if any has <=1 unique value, return TRUE
  any(vapply(aa, function(col) length(unique(col)) <= 1, logical(1)))
}


#' Disjunctive coding for factors (modified from ade4)
#'
#' Converts factor variables in a data frame to disjunctive (dummy) coding.
#'
#' @param df Data frame with factor variables.
#' @return Data frame with dummy-coded variables.
#' @export

acm.disjctif<-function (df)
{
  acm.util.df <- function(i) {
    cl <- df[, i]
    cha <- names(df)[i]
    n <- length(cl)
    cl <- as.factor(cl)
    x <- matrix(0, n, length(levels(cl)))
    x[(1:n) + n * (unclass(cl) - 1)] <- 1
    dimnames(x) <- list(row.names(df), paste(cha, levels(cl),
                                             sep = "."))
    return(x)
  }
  # For FAC(df) with only single variable
  acm.util.df2 <- function(i) {
    cl <- df[, i]
    cha <- colnames(df)[i]
    n <- length(cl)
    cl <- as.factor(cl)
    x <- matrix(0, n, length(levels(cl)))
    x[(1:n) + n * (unclass(cl) - 1)] <- 1
    dimnames(x) <- list(row.names(df), paste(cha, levels(cl),
                                             sep = "."))
    return(x)
  }
  if(!is.null(ncol(df)) & ncol(df)>1) G <- lapply(1:ncol(df), acm.util.df)
  if(!is.null(ncol(df)) & ncol(df)==1) G <- lapply(1, acm.util.df2)
  G <- data.frame(G, check.names = FALSE)
  return(G)
}


#' Dummy coding for data frame (numeric and factor)
#'
#' Converts numeric and factor variables in a data frame to dummy-coded format.
#'
#' @param df Data frame.
#' @return Data frame with numeric and dummy-coded factor variables.
#' @export

dummy2 <- function(df) {
  NUM <- function(dataframe)dataframe[,sapply(dataframe,is.numeric)]
  FAC <- function(dataframe)dataframe[,sapply(dataframe,is.factor)]
  if (is.null(ncol(NUM(df)))){
    DF <- data.frame(NUM(df), acm.disjctif(FAC(df)))
    names(DF)[1] <- colnames(df)[which(sapply(df, is.numeric))]
  } else {
    if (!is.null(ncol(FAC(df))) && ncol(FAC(df))>0) DF <- data.frame(NUM(df), acm.disjctif(FAC(df)))
    if (!is.null(ncol(FAC(df))) | ncol(FAC(df))==0) DF <- data.frame(NUM(df))
    if (is.null(ncol(FAC(df)))){
      temp <- as.matrix(FAC(df))
      colnames(temp)[1] <- colnames(df)[which(sapply(df, is.factor))]
      df.fac <- acm.disjctif(temp)
      DF <- data.frame(NUM(df), df.fac)
    }
  }
  return(DF)
}

#' Dummy coding for data frame (numeric and factor)
#'
#' Converts numeric and factor variables in a data frame to dummy-coded format.
#'
#' @param df Data frame.
#' @return Data frame with numeric and dummy-coded factor variables.
#' @export

dummy <- function(df) {
  NUM <- function(dataframe)dataframe[,sapply(dataframe,is.numeric)]
  FAC <- function(dataframe)dataframe[,sapply(dataframe,is.factor)]
  if (is.null(ncol(NUM(df)))){
    DF <- data.frame(NUM(df), acm.disjctif(FAC(df)))
    names(DF)[1] <- colnames(df)[which(sapply(df, is.numeric))]
  } else {
    if (!is.null(ncol(FAC(df)))) DF <- data.frame(NUM(df), acm.disjctif(FAC(df)))
    if (is.null(ncol(FAC(df)))){
      temp <- as.matrix(FAC(df))
      colnames(temp)[1] <- colnames(df)[which(sapply(df, is.factor))]
      df.fac <- acm.disjctif(temp)
      DF <- data.frame(NUM(df), df.fac)
    }
  }
  return(DF)
}

#' Trailing zeros in binary representation
#'
#' Returns the number of trailing zeros in the binary representation of an integer.
#'
#' @param kk Integer.
#' @return Integer count of trailing zeros.
#' @export

ztrail <- function(kk){
  ii <- 1
  zz <- kk
  while(zz%%2 == 0){
    ii <- ii+1
    zz <- zz/2
  }
  return(ii)
}

#' Flip 1/0 value
#'
#' Flips a binary value: returns 0 if input is 1, returns 1 if input is 0.
#'
#' @param x Integer (0 or 1).
#' @return Integer (0 or 1).
#' @export

one.zero <- function(x){
  if (x == 1){
    x <- 0}
  else {
    x <- 1}
  return(x)
}
