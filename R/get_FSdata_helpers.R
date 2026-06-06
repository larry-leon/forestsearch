
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

#' k-th J-Quantile
#'
#' Returns the (k/J)-th quantile of a numeric vector.  Used by
#' \code{cut_var_jq()} to emit deferred cut expressions of the form
#' \code{"X <= qj(X, k, J)"}, which are then resolved to literal
#' numerics by \code{process_conf_force_expr()}.
#'
#' @param x A numeric vector.
#' @param k Integer in \code{1, ..., J - 1L}.  Index of the cut point
#'   (where \code{J} here is the third argument to \code{qj()}, not the
#'   user-facing \code{conf.cont_jcuts} value -- see \code{cut_var_jq()}
#'   for how these relate).
#' @param J Integer >= 2.  Total number of intervals implied by the
#'   probability \code{k/J}.  In the cut expressions emitted by
#'   \code{cut_var_jq(x, J_user)}, this argument is \code{J_user + 1}.
#' @return Numeric value of the (k/J)-th quantile of \code{x}.
#' @importFrom stats quantile
#' @keywords internal

qj <- function(x, k, J) c(quantile(x, probs = k / J, na.rm = TRUE))

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

#' Generate J-quantile cut expressions for a continuous variable
#'
#' For a continuous variable, returns \code{J} cut expressions of the
#' form \code{"X <= qj(X, k, J + 1)"} for \code{k = 1, ..., J}.  These
#' \code{J} cut points are placed at the (k/(J+1))-th empirical
#' quantiles of X and partition its range into \code{J + 1}
#' non-overlapping intervals
#' \deqn{[\min(X), c_1),\ [c_1, c_2),\ \ldots,\ [c_J, \max(X)]}
#' where \eqn{c_k} is the (k/(J+1))-th quantile of X.  Note the cut
#' expressions themselves are nested half-spaces (each is a subset of
#' the next); ForestSearch consumes them as binary candidate factors
#' and combines them via intersection during the search.
#'
#' Expressions are emitted in deferred form (with literal \code{qj(...)}
#' calls inside the string) so that they are correctly recomputed when
#' the same expression is processed against a different data subset
#' (e.g., a bootstrap replicate).  They are subsequently resolved to
#' literal numerics by \code{process_conf_force_expr()}.
#'
#' @param x Character.  Variable name.
#' @param J Integer >= 1.  Number of binary cut expressions to emit.
#'   The resulting partition has \code{J + 1} non-overlapping intervals.
#'   Cuts are placed at the (k/(J+1))-th empirical quantiles for
#'   \code{k = 1, ..., J} (no cut at the boundaries).
#' @return Character vector of \code{J} cut expressions.
#' @keywords internal

cut_var_jq <- function(x, J) {
  J <- as.integer(J)
  if (length(J) != 1L || is.na(J) || J < 1L) {
    stop("cut_var_jq: J must be a single positive integer.", call. = FALSE)
  }
  ks <- seq_len(J)
  paste0(x, " <= qj(", x, ", ", ks, ", ", J + 1L, ")")
}

#' Get forced cut expressions for variables
#'
#' For each variable in \code{conf.force.names}, returns cut expressions if continuous.
#'
#' @param df Data frame.
#' @param conf.force.names Character vector of variable names.
#' @param cont.cutoff Integer. Cutoff for continuous.
#' @return Character vector of cut expressions.
#' @examples
#' df <- data.frame(age = c(45, 60, 35), size = c(20, 35, 15), grade = c(1, 2, 3))
#' get_conf_force(df, conf.force.names = c("age", "size"))
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

#' Get J-quantile cut expressions for variables
#'
#' For each named variable in \code{conf_jcuts}, returns \code{J}
#' deferred cut expressions of the form \code{"X <= qj(X, k, J + 1)"}
#' via \code{cut_var_jq()}.  Variables that are not continuous (per
#' \code{cont.cutoff}) are skipped with a warning.  This is the
#' \code{cut_var_jq}/J-quantile analog of \code{get_conf_force()}.
#'
#' Conflicts with other forced-cut mechanisms (\code{defaultcut_names},
#' \code{conf.cont_medians_force}) are validated upstream by
#' \code{get_FSdata()}; this helper performs only structural checks.
#'
#' @param df Data frame.
#' @param conf_jcuts Named list of positive integers.  Each name is a
#'   continuous-variable column in \code{df}; each value \code{J} is
#'   the number of binary cut expressions to emit for that variable.
#'   The cuts are placed at the (k/(J+1))-th empirical quantiles for
#'   \code{k = 1, ..., J}, partitioning the variable's range into
#'   \code{J + 1} non-overlapping intervals.
#' @param cont.cutoff Integer.  Cutoff for continuous determination
#'   (passed through to \code{is.continuous()}).
#' @return Character vector of cut expressions (length \code{sum(J_v)}
#'   over continuous \code{v}), or \code{character(0)} if none.
#' @examples
#' df <- data.frame(age = c(45, 60, 35, 50, 70, 25, 80, 55))
#' # 5 binary cuts at the 1/6, 2/6, ..., 5/6 quantiles of age:
#' get_conf_force_jq(df, conf_jcuts = list(age = 5))
#' @export

get_conf_force_jq <- function(df, conf_jcuts, cont.cutoff = 4) {
  if (!is.data.frame(df)) stop("df must be a data.frame.", call. = FALSE)
  if (is.null(conf_jcuts) || length(conf_jcuts) == 0L) return(character(0))
  if (!is.list(conf_jcuts) || is.null(names(conf_jcuts)) ||
      any(!nzchar(names(conf_jcuts)))) {
    stop("conf_jcuts must be a NAMED list, e.g. list(X8 = 10).",
         call. = FALSE)
  }
  res <- list()
  for (name in names(conf_jcuts)) {
    if (!name %in% names(df)) {
      warning(sprintf(
        "Variable '%s' not found in data frame. Skipping.", name),
        call. = FALSE)
      next
    }
    var_data <- df[[name]]
    flag_cont <- is.continuous(var_data, cutoff = cont.cutoff)
    if (flag_cont != 1) {
      warning(sprintf(
        "Variable '%s' is not continuous (cont.cutoff = %d). Skipping J-quantile cuts.",
        name, cont.cutoff), call. = FALSE)
      next
    }
    J <- conf_jcuts[[name]]
    res[[name]] <- cut_var_jq(x = name, J = J)
  }
  out <- unlist(res, use.names = FALSE)
  if (is.null(out)) character(0) else out
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
#' @param outcome_type Character. One of \code{"survival"} (default),
#'   \code{"binary"}, \code{"continuous"}, or \code{"count"}.
#' @param offset.name Character or \code{NULL}. Name of the follow-up time
#'   column for rate-based measures (IRR, IRD).
#' @return List with selected, omitted variables, coefficients, lambda, and fits.
#' @importFrom glmnet cv.glmnet glmnet
#' @importFrom survival Surv
#' @examples
#' \donttest{
#' library(survival)
#' df <- survival::gbsg
#' df$grade3 <- as.integer(df$grade == "3")
#' lasso_selection(df,
#'   confounders.name = c("age", "meno", "size", "grade3", "nodes", "pgr", "er"),
#'   outcome.name = "rfstime", event.name = "status")
#' }
#' @export

lasso_selection <- function(df, confounders.name, outcome.name, event.name,
                           seedit = 8316951,
                           outcome_type = "survival",
                           offset.name = NULL) {
# Scope the seed locally so this function does not pollute the caller's
# RNG stream.  Save the current state, set the local seed, and restore
# on exit.  If seedit is NULL, leave the RNG untouched entirely.
if (!is.null(seedit)) {
  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv)) {
    get(".Random.seed", envir = .GlobalEnv)
  } else NULL
  on.exit({
    if (!is.null(old_seed)) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  set.seed(seedit)
}
# Package checks
if (!requireNamespace("glmnet", quietly = TRUE)) stop("Package 'glmnet' is required.")
if (!requireNamespace("survival", quietly = TRUE)) stop("Package 'survival' is required.")

  # Prepare the design matrix and response
  # Coerce factor columns with all-numeric levels to integer so that
  # as.matrix() produces a numeric matrix (not a character one).
  # as.integer(as.character(f)) preserves the numeric value of each level
  # (e.g. factor "1" -> 1L), unlike plain as.integer(f) which gives 1-based codes.
  x_df <- df[, confounders.name, drop = FALSE]
  for (.nm in names(x_df)) {
    if (is.factor(x_df[[.nm]])) {
      .lvls <- levels(x_df[[.nm]])
      if (!anyNA(suppressWarnings(as.numeric(.lvls))))
        x_df[[.nm]] <- as.integer(as.character(x_df[[.nm]]))
    }
  }
  x <- as.matrix(x_df)

  # Branch response and family by outcome_type
  if (outcome_type == "survival") {
    y <- survival::Surv(df[[outcome.name]], df[[event.name]])
    glm_family <- "cox"
    offset_vec <- NULL
  } else if (outcome_type == "binary") {
    y <- df[[outcome.name]]
    # Detect rate-based (Poisson) if offset is provided
    if (!is.null(offset.name)) {
      glm_family <- "poisson"
      offset_vec <- log(df[[offset.name]])
    } else {
      glm_family <- "binomial"
      offset_vec <- NULL
    }
  } else if (outcome_type == "continuous") {
    y <- df[[outcome.name]]
    glm_family <- "gaussian"
    offset_vec <- NULL
  } else if (outcome_type == "count") {
    y <- df[[outcome.name]]
    glm_family <- "poisson"
    offset_vec <- if (!is.null(offset.name)) log(df[[offset.name]]) else NULL
  } else {
    stop("Unknown outcome_type: ", outcome_type)
  }

  # Fit LASSO with cross-validation
  # For family = "cox", pass cox.ties = "efron" to (i) align tie-handling with
  # survival::coxph (which is what downstream HR estimation uses), and
  # (ii) silence the glmnet v5.0 -> v5.1 default-change transition warning.
  # Argument is no-op for other families; only built for cox to keep intent
  # explicit and to avoid passing it through if glmnet ever tightens its API.
  glmnet_extra <- if (identical(glm_family, "cox")) {
    list(cox.ties = "efron")
  } else {
    list()
  }
  cvfit <- do.call(
    glmnet::cv.glmnet,
    c(list(x = x, y = y, family = glm_family, alpha = 1,
           offset = offset_vec),
      glmnet_extra)
  )
  lambda_min <- cvfit$lambda.min
  fit <- do.call(
    glmnet::glmnet,
    c(list(x = x, y = y, family = glm_family, alpha = 1,
           lambda = lambda_min, offset = offset_vec),
      glmnet_extra)
  )

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
#' @examples
#' filter_by_lassokeep(c("age", "size", "nodes"), lassokeep = c("age", "nodes"))
#' filter_by_lassokeep(c("pgr", "er"), lassokeep = c("age", "nodes"))  # returns NULL
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
#' @examples
#' df <- data.frame(age = c(40, 55, 70, 35), size = c(20, 30, 25, 15))
#' process_conf_force_expr("age <= mean(age)", df)
#' process_conf_force_expr("size <= median(size)", df)
#' @export

process_conf_force_expr <- function(expr, df) {
  # Try the 3-arg qj() pattern first.
  # Examples: "X8 <= qj(X8, 1, 10)"
  pattern_qj <- "^\\s*([a-zA-Z0-9_.]+)\\s*<=\\s*qj\\(\\s*([a-zA-Z0-9_.]+)\\s*,\\s*([0-9]+)\\s*,\\s*([0-9]+)\\s*\\)\\s*$"
  m_qj <- regexec(pattern_qj, expr)
  matches_qj <- regmatches(expr, m_qj)[[1]]
  if (length(matches_qj) > 0) {
    var <- matches_qj[2]
    arg <- matches_qj[3]
    k <- as.integer(matches_qj[4])
    J <- as.integer(matches_qj[5])
    if (!(var %in% colnames(df)) || !(arg %in% colnames(df))) return(expr)
    if (is.na(k) || is.na(J) || J < 2L || k < 1L || k >= J) return(expr)
    col_vals <- df[[arg]]
    if (is.factor(col_vals)) col_vals <- as.numeric(as.character(col_vals))
    val <- round(quantile(col_vals, probs = k / J, na.rm = TRUE), 1)
    return(paste0(var, " <= ", val))
  }

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
  # Coerce factor columns to numeric before computing summary statistics
  col_vals <- df[[arg]]
  if (is.factor(col_vals)) col_vals <- as.numeric(as.character(col_vals))
  # Evaluate the function
  if (fun == "mean") {
    val <- round(mean(col_vals, na.rm = TRUE), 1)
  } else if (fun == "median") {
    val <- round(median(col_vals, na.rm = TRUE), 1)
  } else if (fun == "qlow") {
    val <- round(quantile(col_vals, 0.25, na.rm = TRUE), 1)
  } else if (fun == "qhigh") {
    val <- round(quantile(col_vals, 0.75, na.rm = TRUE), 1)
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


#' Round half away from zero
#'
#' Base R \code{round()} uses round-half-to-even ("banker's rounding"), so
#' \code{round(0.5) == 0} and \code{round(74.5) == 74}.  Cut coarsening rounds
#' representative thresholds half *up* (away from zero) so that, e.g., a cluster
#' centroid of 74.5 maps to 75, matching the intuitive "nearest integer" rule.
#'
#' @param x Numeric vector.
#' @param digits Integer.  Number of decimal places.
#' @return Numeric vector rounded half away from zero.
#' @keywords internal

collapse_cuts_round <- function(x, digits = 0L) {
  f <- 10^as.integer(digits)
  sign(x) * floor(abs(x) * f + 0.5) / f
}

#' Format a numeric threshold for a cut expression
#'
#' Renders a (already-rounded) numeric threshold as a clean string with no
#' trailing zeros or floating-point noise, for use in "var <= value" cut
#' expressions.  Integer-valued thresholds print without a decimal point.
#'
#' @param v Numeric scalar.
#' @param digits Integer.  Decimal places used during rounding.
#' @return Character scalar.
#' @keywords internal

collapse_cuts_fmt <- function(v, digits = 0L) {
  s <- formatC(v, format = "f", digits = max(0L, as.integer(digits)))
  if (as.integer(digits) > 0L) {
    s <- sub("0+$", "", s)
    s <- sub("\\.$", "", s)
  }
  trimws(s)
}

#' Tidy a cut / subgroup expression for display
#'
#' Display-only cleanup for diagnostic printing of cut and subgroup-definition
#' strings.  Collapses runs of whitespace to single spaces (so padded GRF
#' definitions like \code{"cd40 <=                320"} render cleanly) and
#' rounds each post-operator numeric threshold half-up to \code{digits} places
#' (matching the \code{collapse_cuts} representative rounding, e.g.
#' \code{"wtkg > 79.833600000000004"} -> \code{"wtkg > 80"} at \code{digits = 0}).
#'
#' Only numbers that follow a comparison operator (\code{<=}, \code{>=},
#' \code{==}, \code{<}, \code{>}) are rounded, so integers embedded in variable
#' names (\code{cd40}, \code{prior_6mo}) are never altered.  This function does
#' NOT touch the canonical, membership-exact labels used in the search; it is
#' intended only for the diagnostic strings passed to \code{cat()}.
#'
#' @param s Character vector of expressions, e.g.
#'   \code{"cd40 <=    320 & wtkg > 79.8336"}.
#' @param digits Integer >= 0.  Rounding for embedded thresholds.  Default 0
#'   (nearest integer), matching the \code{collapse_cuts} default.
#' @return Character vector, tidied for display (same length as \code{s}).
#' @keywords internal

tidy_cut_display <- function(s, digits = 0L) {
  if (length(s) == 0L) return(s)
  pat <- "(<=|>=|==|<|>)\\s*(-?[0-9.]+([eE][-+]?[0-9]+)?)"
  vapply(s, function(x) {
    if (is.na(x) || !nzchar(x)) return(x)
    x <- gsub("[[:space:]]+", " ", x)            # collapse padding
    m <- gregexpr(pat, x, perl = TRUE)[[1L]]
    if (m[1L] == -1L) return(trimws(x))
    L <- attr(m, "match.length")
    for (i in rev(seq_along(m))) {               # right-to-left: positions stay valid
      tok <- substr(x, m[i], m[i] + L[i] - 1L)
      gp  <- regmatches(tok, regexec(
        "^(<=|>=|==|<|>)\\s*(-?[0-9.]+([eE][-+]?[0-9]+)?)$", tok))[[1L]]
      if (length(gp) < 3L) next
      num <- suppressWarnings(as.numeric(gp[3L]))
      if (is.na(num)) next
      rep <- paste0(gp[2L], " ",
                    collapse_cuts_fmt(collapse_cuts_round(num, digits), digits))
      x <- paste0(substr(x, 1L, m[i] - 1L), rep, substr(x, m[i] + L[i], nchar(x)))
    }
    trimws(x)
  }, character(1), USE.NAMES = FALSE)
}

#' Collapse near-redundant continuous candidate cuts
#'
#' Continuous covariates can generate many candidate threshold cuts that are
#' practically redundant -- e.g. "age <= 35.7" and "age <= 35", or
#' "wtkg <= 75.2" and "wtkg <= 74.5" -- because the candidate pool unions
#' quantile cuts, GRF / DINA splits, and calibrated-DGM cuts at different
#' precisions.  This helper merges such near-duplicates to a single
#' representative threshold, leaving categorical and indicator cuts untouched.
#'
#' For each (variable, operator) group of continuous-variable cuts, distinct
#' thresholds are single-linkage clustered using a per-variable band
#' \code{band = c_band * sd(x) / sqrt(n)} (the standard error of the variable:
#' thresholds finer than ~1 SE are not statistically resolvable).  Each cluster
#' collapses to one representative, the cluster mean rounded half-up to
#' \code{digits} places (see \code{collapse_cuts_round()}).  Rounding alone can
#' also merge singletons that round to the same value (e.g. 75.2 and 74.5 both
#' round to 75 at \code{digits = 0}); de-duplication then drops the duplicate.
#'
#' Safety check: a cluster is collapsed only if replacing each member threshold
#' by the representative changes subgroup membership for no more than
#' \code{safety_tol} subjects (a fraction of n when \code{safety_tol < 1}, an
#' absolute count when \code{safety_tol >= 1}).  Clusters that would move more
#' than that are kept unchanged, so coarsening never silently redefines a
#' candidate subgroup by more than the tolerance.
#'
#' Cut expressions that are categorical, indicator-valued (fewer than
#' \code{cont.cutoff} unique values), bare variable names, or equality tests
#' (\code{==}) pass through unchanged.
#'
#' @param cuts Character vector of candidate cut expressions, e.g.
#'   \code{"age <= 35.7"}, \code{"wtkg >= 87.2"}, \code{"gender"}.
#' @param df Data frame supplying the covariate columns, used for sd(x), the
#'   sample size n, and the membership safety check.
#' @param confounders.name Character vector of confounder names (currently
#'   informational; variable identity is parsed from each cut expression).
#' @param c_band Numeric >= 0.  Band multiplier; \code{band = c_band *
#'   sd(x) / sqrt(n)}.  Larger values merge more aggressively.  Default 1.
#' @param safety_tol Numeric > 0.  Membership safety tolerance: a fraction of
#'   n when < 1, an absolute subject count when >= 1.  Default 0.05.
#' @param digits Integer >= 0.  Decimal places for the representative
#'   threshold.  Default 0 (nearest integer).  Raise for variables measured on
#'   a sub-unit scale.
#' @param cont.cutoff Integer.  A variable with at least this many unique values
#'   is treated as continuous (see \code{is.continuous()}).  Default 4.
#' @param details Logical.  If TRUE, print a per-cluster collapse report.
#' @return Character vector of de-duplicated cut expressions (length <= input).
#' @keywords internal

collapse_redundant_cuts <- function(cuts, df, confounders.name = NULL,
                                    c_band = 1.0, safety_tol = 0.05,
                                    digits = 0L, cont.cutoff = 4,
                                    details = FALSE) {
  if (length(cuts) == 0L) return(cuts)
  n <- nrow(df)
  if (is.null(n) || n == 0L) return(unique(cuts))

  pat <- "^\\s*([A-Za-z0-9_.]+)\\s*(<=|>=|<|>)\\s*(-?[0-9.]+([eE][-+]?[0-9]+)?)\\s*$"
  ops <- list("<=" = function(a, b) a <= b,
              ">=" = function(a, b) a >= b,
              "<"  = function(a, b) a <  b,
              ">"  = function(a, b) a >  b)

  # Parse each cut: classify as coarsenable (continuous numeric threshold) or
  # passthrough (categorical / indicator / bare name / equality test).
  mm     <- regmatches(cuts, regexec(pat, cuts))
  is_num <- vapply(mm, function(z) length(z) > 0L, logical(1))
  var_of <- vapply(seq_along(cuts),
                   function(i) if (is_num[i]) mm[[i]][2L] else NA_character_,
                   character(1))
  op_of  <- vapply(seq_along(cuts),
                   function(i) if (is_num[i]) mm[[i]][3L] else NA_character_,
                   character(1))
  val_of <- vapply(seq_along(cuts),
                   function(i) if (is_num[i]) as.numeric(mm[[i]][4L]) else NA_real_,
                   numeric(1))

  # Restrict coarsening to continuous variables present in df.
  coarsen <- is_num & !is.na(var_of) & var_of %in% names(df)
  for (i in which(coarsen)) {
    if (is.continuous(df[[var_of[i]]], cutoff = cont.cutoff) != 1L) {
      coarsen[i] <- FALSE
    }
  }

  out          <- cuts   # default: every cut unchanged
  safety_bound <- if (safety_tol < 1) ceiling(safety_tol * n) else as.integer(safety_tol)
  report       <- character(0)

  grp_key <- paste(var_of, op_of, sep = "\r")   # group by (variable, operator)
  for (key in unique(grp_key[coarsen])) {
    idx  <- which(coarsen & grp_key == key)
    v    <- var_of[idx[1L]]
    op   <- op_of[idx[1L]]
    x    <- df[[v]]
    sdx  <- stats::sd(x, na.rm = TRUE)
    band <- if (is.finite(sdx) && sdx > 0) c_band * sdx / sqrt(n) else 0
    opf  <- ops[[op]]

    vals  <- sort(unique(val_of[idx]))
    cl_id <- if (length(vals) == 1L) 1L
             else cumsum(c(1L, as.integer(diff(vals) > band)))

    for (cid in unique(cl_id)) {
      members <- vals[cl_id == cid]
      rep_val <- collapse_cuts_round(mean(members), digits)
      mr      <- opf(x, rep_val)
      flips   <- max(vapply(members,
                            function(vv) sum(opf(x, vv) != mr, na.rm = TRUE),
                            numeric(1)))
      changed <- (length(members) > 1L) ||
                 !isTRUE(all.equal(rep_val, members[1L]))
      if (flips <= safety_bound) {
        rep_str <- paste0(v, " ", op, " ", collapse_cuts_fmt(rep_val, digits))
        for (vv in members) {
          hit <- idx[abs(val_of[idx] - vv) < .Machine$double.eps^0.5]
          out[hit] <- rep_str
        }
        if (changed) {
          report <- c(report, sprintf(
            "  %s %s {%s} -> %s  (%d subj moved)",
            v, op, paste(format(members), collapse = ", "), rep_str, flips))
        }
      } else {
        report <- c(report, sprintf(
          "  SAFETY-BLOCKED %s %s {%s} -> %s  (%d > %d subj); kept unchanged",
          v, op, paste(format(members), collapse = ", "),
          collapse_cuts_fmt(rep_val, digits), flips, safety_bound))
      }
    }
  }

  collapsed <- out[!duplicated(out)]
  if (isTRUE(details)) {
    cat(sprintf(
      "\nCut coarsening: %d -> %d candidates  (c=%.3g, safety_tol=%.3g, digits=%d)\n",
      length(cuts), length(collapsed), c_band, safety_tol, as.integer(digits)))
    if (length(report) > 0L) cat(paste(report, collapse = "\n"), "\n")
  }
  collapsed
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
#' @examples
#' df <- data.frame(age = c(40, 55, 70), grade = factor(c("1", "2", "3")))
#' dummy_encode(df)
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
#' @examples
#' df <- data.frame(age = c(40, 55, 70), size = c(20, 30, 25))
#' confs <- c("age <= 50", "size > 25")
#' evaluate_cuts_once(confs, df)
#' @export

evaluate_cuts_once <- function(confs, df, details = FALSE) {
  n_confs <- length(confs)
  evaluations <- vector("list", n_confs)
  is_valid <- logical(n_confs)
  has_error <- logical(n_confs)

  for (i in seq_along(confs)) {
    thiscut <- confs[i]

    result_i <- tryCatch({
      # Use evaluate_comparison() -- no eval(parse()) needed
      # Cut expressions are always single comparisons like "er <= 0"
      result <- evaluate_comparison(thiscut, df)
      list(
        evaluation = as.logical(result),
        is_valid   = length(unique(result)) > 1,
        has_error  = FALSE
      )
    }, error = function(e) {
      if (details) {
        cat("Error evaluating cut '", thiscut, "': ", e$message, "\n", sep = "")
      }
      list(evaluation = NULL, is_valid = FALSE, has_error = TRUE)
    })

    evaluations[[i]] <- result_i$evaluation
    is_valid[i]      <- result_i$is_valid
    has_error[i]     <- result_i$has_error
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



